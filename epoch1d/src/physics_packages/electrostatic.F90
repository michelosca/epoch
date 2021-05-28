! Copyright (C) 2009-2019 University of Warwick
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE electrostatic
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
  USE petscvec
  USE petscmat
  USE petscksp
  USE shared_data
  USE custom_electrostatic
  USE evaluator
  USE boundary
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: es_update_e_field, setup_electrostatic
  PUBLIC :: init_potential, es_wall_current_diagnostic
  PUBLIC :: attach_potential, setup_petsc_variables, finalize_petsc
  PUBLIC :: potential_update_spatial_profile, min_init_electrostatic
  PUBLIC :: initialize_petsc, finalize_electrostatic_solver

  Vec, SAVE :: es_potential_vec, charge_dens_vec
  Mat, SAVE :: transform_mtrx
  MatNullSpace, SAVE :: matnull
  PetscErrorCode, SAVE :: perr
  KSP, SAVE :: ksp


CONTAINS

  SUBROUTINE es_update_e_field
    
    REAL(num) :: rhodx_min, rhodx_max
    REAL(num), DIMENSION(:), ALLOCATABLE :: es_charge_density

    ! Charge weighting from particles to the grid, i.e. charge density
    ALLOCATE(es_charge_density(1-ng:nx+ng))
    CALL es_calc_charge_density(es_charge_density)

    ! Charge density to electrostatic potential
    !  - This subroutine calculates the electric potential on es_potential
    CALL es_calc_potential(es_charge_density)

    ! Calculate electric field in x-direction
    !  - In case of non-periodic boundaries
    IF (x_min_boundary_open) THEN
      ! chargdens_1/2 = (rho_0+rho_1)/2*dx: charge density near wall
      rhodx_min = (es_charge_density(0)+es_charge_density(1)) * dx * 0.5_num
    END IF
    IF (x_max_boundary_open) THEN
      ! chargdens_nx-1/2 = rho_nx*dx/2: charge density near wall
      rhodx_max = (es_charge_density(nx)+es_charge_density(nx-1)) * dx * 0.5_num
    END IF
    DEALLOCATE(es_charge_density)
    
    ! - Calculate electric field in x-direction
    CALL es_calc_ex(rhodx_min, rhodx_max)
    
    ! - Update electric field in y- and z-direction
    CALL set_ez
    CALL set_ey


  END SUBROUTINE es_update_e_field



  SUBROUTINE es_calc_charge_density(charge_density)

    REAL(num), DIMENSION(1-ng:nx+ng), INTENT(INOUT) :: charge_density
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_q
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: part_weight, idx
    INTEGER :: ispecies, ix, i
    TYPE(particle), POINTER :: current
#include "particle_head.inc"

    idx = 1.0_num / dx
    charge_density = 0._num

    DO ispecies = 1, n_species
      IF (species_list(ispecies)%species_type == c_species_id_photon) CYCLE

      ! Copy the particle properties out for speed
      part_q  = species_list(ispecies)%charge
      IF (ABS(part_q) < TINY(0._num)) CYCLE
#ifdef PER_SPECIES_WEIGHT      
      part_weight = species_list(ispecies)%weight
      wdata = part_q * part_weight
#endif

      current => species_list(ispecies)%attached_list%head
      DO WHILE (ASSOCIATED(current))
#ifndef PER_SPECIES_WEIGHT
        part_weight = current%weight
        wdata = part_q * part_weight
#endif

#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = (current%part_pos - x_min_local) * idx - 0.5_num
#else
        cell_x_r = (current%part_pos - x_min_local) * idx
#endif
        cell_x = FLOOR(cell_x_r + 0.5_num)
        cell_frac_x = REAL(cell_x, num) - cell_x_r

#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gxfac.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gxfac.inc"
#else
#include "triangle/gxfac.inc"
#endif

        DO ix = sf_min, sf_max
          charge_density(cell_x+ix) = charge_density(cell_x+ix) + gx(ix) * wdata
        END DO

        current => current%next
      END DO
    END DO

    ! Adds charge density values from adjacent processors
    CALL processor_summation_bcs(charge_density, ng)
    ! Places charge density values from adjacent processor in ghost cells
    CALL field_bc(charge_density, ng)

    charge_density = charge_density * idx
    IF (x_max_boundary_open) THEN
      DO i = 1,ng-1
        charge_density(nx-i) = charge_density(nx-i) + charge_density(nx+i)
      END DO
      charge_density(nx+1:) = 0._num
      charge_density(nx) = charge_density(nx) * 2._num
    END IF
    IF (x_min_boundary_open) THEN
      DO i = 1,ng-1
        charge_density(i) = charge_density(i) + charge_density(-i)
      END DO
      charge_density(:-1) = 0._num
      charge_density(0) = charge_density(0) * 2._num
    END IF

  END SUBROUTINE es_calc_charge_density



  SUBROUTINE es_calc_potential(charge_density)
  
    REAL(num), DIMENSION(1-ng:nx+ng), INTENT(IN) :: charge_density
    INTEGER :: nx_end
    
    ! The Poisson solver resolves the potential inside the grid with
    !given boundary conditions.
    nx_end = nx
    IF (x_max_boundary) nx_end = nx_end - 1
    
    CALL poisson_solver(charge_density(1:nx_end), nx_end)
    
    ! Potential boundaries (between MPI processors, valid for periodic bc)
    CALL field_bc(es_potential, ng)

    ! Potential (open) boundary conditions
    IF (x_min_boundary) es_potential(0) = set_potential_x_min()
    IF (x_max_boundary) es_potential(nx) = set_potential_x_max()

  END SUBROUTINE es_calc_potential



  SUBROUTINE poisson_solver(charge_density, nx_end)

    ! This subroutine solves Poisson's equation
    ! Input: charge_density array and the length of the data array
    ! Output: electric potential array (data_array)
    INTEGER, INTENT(IN) :: nx_end
    REAL(num), DIMENSION(1:nx_end), INTENT(IN) :: charge_density
    REAL(num), DIMENSION(1:nx_end) :: b_array
    REAL(num) :: constant
    INTEGER :: i
    REAL(num), DIMENSION(:), POINTER :: vec_pointer

    ! Get array ready for solver
    constant = -dx*dx/epsilon0
    b_array = charge_density * constant
    b_array(1) = b_array(1) - set_potential_x_min()
    b_array(nx_end) = b_array(nx_end) - set_potential_x_max()

    ! Get charge density data (data_array) into the Petsc vector
    CALL VecGetArrayF90(charge_dens_vec, vec_pointer, perr)
    DO i = 1, nx_end
      vec_pointer(i) = b_array(i)
    END DO
    CALL VecRestoreArrayF90(charge_dens_vec, vec_pointer, perr)

    ! Solve linear system: transform_mtrx*es_potential_vec = charge_dens_vec
    !                      A * x = b
    CALL KSPSolve(ksp, charge_dens_vec, es_potential_vec, perr)

    ! Pass electric potential data from PETSc to Fortran
    CALL VecGetArrayReadF90(es_potential_vec, vec_pointer, perr)
    DO i = 1, nx_end
      es_potential(i) = vec_pointer(i)
    END DO
    CALL VecRestoreArrayReadF90(es_potential_vec, vec_pointer, perr)

  END SUBROUTINE poisson_solver



  SUBROUTINE es_calc_ex(rhodx_min, rhodx_max)

    INTEGER :: ix
    REAL(num) :: idx
    REAL(num), INTENT(IN) :: rhodx_min, rhodx_max

    idx = 1._num/dx
    ex = 0._num
    
    ! Calculate the electric field: Div(E) = -Phi
    DO ix = 0, nx
      ex(ix) = es_potential(ix-1) - es_potential(ix+1)
    END DO
    ex = ex * 0.5_num * idx

    ! E-field boundaries
    ! For between-processors boundaries and external (periodic) bc
    CALL field_bc(ex, ng)

    ! For non-periodic boundaries
    IF (x_min_boundary_open) THEN
      ! (E_{1} - E_0) / dx = rho_1/2/epsilon0
      ! Wall charge and difference with respect to last time step
      ! Previous wall charge surface density value
      dwcharge_min = wcharge_min  
      ! Charge surface density at wall
      wcharge_min = ex(1) * epsilon0 - rhodx_min 
      ! This is used for diagnostics
      dwcharge_min = wcharge_min - dwcharge_min 

      !Update E-field at wall
      ex(0) = wcharge_min/epsilon0
      ! Option two would be ex(0) = ex_half*2._num - ex(1)
    END IF

    IF (x_max_boundary_open) THEN
      ! (E_{nx} - E_{nx-1}) / dx = rho_{nx-1/2}/epsilon0
      ! Wall charge and difference with respect to last time step
      ! Previous wall charge surface density value
      dwcharge_max = wcharge_max
      ! Charge surface density at wall
      wcharge_max = ex(nx-1) * epsilon0 + rhodx_max 
      ! This is used for diagnostics
      dwcharge_max = wcharge_max - dwcharge_max

      !Update E-field at wall
      ex(nx) = wcharge_max/epsilon0
      ! Option two would be ex(nx) = ex_half*2._num - ex(nx-1)
    END IF

    DO ix = 1, 2*c_ndims
      CALL field_clamp_wall(ex, ng, c_stagger_ex, ix)
    END DO

  END SUBROUTINE es_calc_ex



  SUBROUTINE setup_petsc_vector(ncells_local, ncells_global)
    
    INTEGER, INTENT(IN) :: ncells_local, ncells_global
    PetscInt :: nx_local, nx_glob
      
    ! Birdsall and Langdon, Appendix D:
    ! vector size is (nx_glog-1), hence last processor
    ! must do without the cell at the boundary. The potential at grid=0
    ! and grid=nx are given by boundary conditions phiL and phiR
    nx_local = ncells_local
    IF (x_max_boundary) nx_local = nx_local - 1
    nx_glob = ncells_global - 1 ! don't count for last global cell

    ! Charge density vector
    CALL VecCreateMPI(comm, nx_local, nx_glob, charge_dens_vec, perr)
    CALL VecSetFromOptions(charge_dens_vec, perr)
    ! Electric potential vector
    CALL VecDuplicate(charge_dens_vec, es_potential_vec, perr)
    
  END SUBROUTINE setup_petsc_vector



  SUBROUTINE setup_petsc_matrix(ncells_local, ncells_global)
    
    ! ncells_local: number of cells in the processor
    ! ncells_global: number of all cells (all processors)
    INTEGER, INTENT(IN) :: ncells_local, ncells_global
    INTEGER :: n_first, n_last
    PetscInt :: nx_local, nx_glob, zero
    PetscScalar :: values(3)
    PetscInt :: col(3), row

    ! Birdsall and Langdon, Appendix D:
    ! matrix size is (nx_glob-1)x(nx_glob-1), hence last processor
    ! must do without the cell at the boundary.
    nx_local = ncells_local
    IF (x_max_boundary) nx_local = nx_local - 1
    nx_glob = ncells_global - 1

    ! Petsc subroutines setting up the Matrix
    ! MatCreateAIJ(MPI communicator,
    !             number of local rows,
    !             number of local columns,
    !             number of global rows,
    !             number of global columns,
    !             d_nz: number of non-zeros per row in DIAGONAL portion of
    !                   local submatrix
    !             d_nnz-array -> same as d_nz but specifying each row
    !             o_nz: number of non-zeros per row in the OFF-DIAGONAL portion
    !                   of local submatrix
    !             o_nnz-array -> same as o_nz but specifying each row
    CALL MatCreate(comm, transform_mtrx, perr)
    CALL MatSetSizes(transform_mtrx, nx_local, nx_local, nx_glob, nx_glob, perr)
    CALL MatSetFromOptions(transform_mtrx, perr)
    CALL MatSetUp(transform_mtrx, perr)

    ! Petsc matrix indexes (global). First index in Petsc is zero, hence
    ! Shift cell indexing so that it starts at 0
    n_first = nx_global_min-1
    n_last = nx_global_max-1

    ! Matrix diagonal values
    values(1) = 1._num
    values(2) = -2._num
    values(3) = 1._num

    ! Set matrix values
    col = 0
    IF (x_min_boundary) THEN
      ! Top matrix row
      row = n_first
      col(1) = row
      col(2) = row+1
      !MatSetValues(matrix, #rows,rows-global-indexes, 
      !  #columns, col-global-indexes, insertion_values, insert/add,err )
      CALL MatSetValues(transform_mtrx, 1, row, 2, col(1:2), values(2:3), &
              INSERT_VALUES, perr)
      ! Move to next row to be set
      n_first = n_first + 1
    END IF

    IF (x_max_boundary) THEN
      n_last = n_last - 1
      ! Bottom matrix row
      row = n_last
      col(1) = row-1
      col(2) = row
      !MatSetValues(matrix, #rows,rows-global-indexes, 
      !  #columns, col-global-indexes, insertion_values, insert/add,err )
      CALL MatSetValues(transform_mtrx, 1, row, 2, col(1:2), values(1:2), &
              INSERT_VALUES, perr)
      ! Move to last row to be set
      n_last = n_last - 1
    END IF


    DO row = n_first, n_last
      col(1) = row-1
      col(2) = row
      col(3) = row+1
      CALL MatSetValues(transform_mtrx, 1, row, 3, col, values, INSERT_VALUES,&
        perr)
    END DO

    CALL MatAssemblyBegin(transform_mtrx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(transform_mtrx, MAT_FINAL_ASSEMBLY, perr)

   ! Generate nullspace required for GAMG preconditioner
    zero = 0
    CALL MatNullSpaceCreate( comm, PETSC_TRUE, zero, &
      PETSC_NULL_VEC, matnull, perr)
    CALL MatSetNearNullSpace(transform_mtrx, matnull, perr)

    CALL MatSetOption(transform_mtrx, MAT_SYMMETRIC, PETSC_TRUE, perr)

    CALL MatSetOption(transform_mtrx, MAT_SYMMETRIC, PETSC_TRUE, perr)

  END SUBROUTINE setup_petsc_matrix



  SUBROUTINE setup_petsc_ksp

    ! The actual KSP setup happens through the file
    ! 'housekeeping/petsc_runtime_options.opt'

    ! Create linear solver context
    CALL KSPCreate(comm, ksp, perr)

    ! Set operators
    CALL KSPSetOperators(ksp, transform_mtrx, transform_mtrx, perr)

    ! Setup solver
    CALL KSPSetFromOptions(ksp,perr)
    CALL KSPSetUp(ksp, perr)


  END SUBROUTINE setup_petsc_ksp



  SUBROUTINE destroy_petsc

    CALL VecDestroy(es_potential_vec, perr)
    CALL VecDestroy(charge_dens_vec, perr)
    CALL MatDestroy(transform_mtrx, perr)
    CALL MatNullSpaceDestroy(matnull, perr)
    CALL KSPDestroy(ksp, perr)

  END SUBROUTINE destroy_petsc



  SUBROUTINE init_potential(boundary, potential)

    INTEGER, INTENT(IN) :: boundary
    TYPE(potential_block), INTENT(INOUT) :: potential

    potential%boundary = boundary
    potential%id = -1
    potential%use_time_function = .FALSE.
    potential%use_profile_function = .FALSE.
    potential%amp = 0.0_num
    potential%t_start = 0.0_num
    potential%t_end = HUGE(0._num)
    NULLIFY(potential%next)

    potential%profile = 1.0_num

  END SUBROUTINE init_potential
  
  
  
  SUBROUTINE deallocate_potential(potential)

    TYPE(potential_block), POINTER :: potential

    IF (potential%use_profile_function) &
        CALL deallocate_stack(potential%profile_function)
    IF (potential%use_time_function) &
        CALL deallocate_stack(potential%time_function)
    DEALLOCATE(potential)

  END SUBROUTINE deallocate_potential
  
  
  
  SUBROUTINE deallocate_potentials

    TYPE(potential_block), POINTER :: current, next

    current => potential_list_x_min
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL deallocate_potential(current)
      current => next
    END DO

    current => potential_list_x_max
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL deallocate_potential(current)
      current => next
    END DO

  END SUBROUTINE deallocate_potentials
  
  

  SUBROUTINE deallocate_efield_profiles

    IF (ey_profile%use_profile_function) THEN
      CALL deallocate_stack(ey_profile%profile_function)
    END IF
    DEALLOCATE(ey_profile)

    IF (ez_profile%use_profile_function) THEN
      CALL deallocate_stack(ez_profile%profile_function)
    END IF
    DEALLOCATE(ez_profile)

  END SUBROUTINE deallocate_efield_profiles



  SUBROUTINE attach_potential(potential)

    INTEGER :: boundary
    TYPE(potential_block), POINTER :: potential

    boundary = potential%boundary

    IF (boundary == c_bd_x_min) THEN
      n_potential_source_x_min = n_potential_source_x_min + 1
      CALL attach_potential_to_list(potential_list_x_min, potential)
    ELSE IF (boundary == c_bd_x_max) THEN
      n_potential_source_x_max = n_potential_source_x_max + 1
      CALL attach_potential_to_list(potential_list_x_max, potential)
    END IF

  END SUBROUTINE attach_potential
  
  
  
  SUBROUTINE attach_potential_to_list(list, potential)

    TYPE(potential_block), POINTER :: list
    TYPE(potential_block), POINTER :: potential
    TYPE(potential_block), POINTER :: current

    IF (ASSOCIATED(list)) THEN
      current => list
      DO WHILE(ASSOCIATED(current%next))
        current => current%next
      END DO
      current%next => potential
    ELSE
      list => potential
    END IF

  END SUBROUTINE attach_potential_to_list
  
  
  
  SUBROUTINE populate_pack_from_potential(potential, parameters)

    TYPE(potential_block), POINTER :: potential
    TYPE(parameter_pack), INTENT(INOUT) :: parameters

    parameters%pack_ix = 0

    SELECT CASE(potential%boundary)
      CASE(c_bd_x_min)
        parameters%pack_ix = 0
      CASE(c_bd_x_max)
        parameters%pack_ix = nx
    END SELECT

  END SUBROUTINE populate_pack_from_potential
  
  
  
  SUBROUTINE potential_update_spatial_profile(potential)

    TYPE(potential_block), POINTER :: potential
    INTEGER :: err
    TYPE(parameter_pack) :: parameters

    err = 0
    CALL populate_pack_from_potential(potential, parameters)
    IF (potential%use_profile_function) THEN
      potential%profile = &
        evaluate_with_parameters(potential%profile_function, parameters, err)
      RETURN
    END IF
    
  END SUBROUTINE potential_update_spatial_profile
  
  
  
  FUNCTION potential_time_profile(potential)

    TYPE(potential_block), POINTER :: potential
    REAL(num) :: potential_time_profile
    INTEGER :: err
    TYPE(parameter_pack) :: parameters

    err = 0
    CALL populate_pack_from_potential(potential, parameters)
    IF (potential%use_time_function) THEN
      potential_time_profile = evaluate_with_parameters( &
        potential%time_function, parameters, err)
      RETURN
    END IF

    potential_time_profile = custom_electric_potential_time_profile(potential)

  END FUNCTION potential_time_profile
  
  
  
  FUNCTION set_potential_x_min()
  
    REAL(num) :: set_potential_x_min
    REAL(num) :: pot_x_min, source
    TYPE(potential_block), POINTER :: current
  
    pot_x_min = 0._num
    
    IF (x_min_boundary_open) THEN
      IF (add_potential_source(c_bd_x_min)) THEN
        current => potential_list_x_min
          DO WHILE(ASSOCIATED(current))
            ! evaluate the temporal evolution of the potential source
            IF (time >= current%t_start .AND. time <= current%t_end) THEN
              IF (current%use_profile_function) &
                CALL potential_update_spatial_profile(current)
              source = potential_time_profile(current) * current%amp
              source = source * current%profile
              pot_x_min = pot_x_min + source
            END IF
            current => current%next
          END DO
      END IF
    END IF
    
    set_potential_x_min = pot_x_min
   
  END FUNCTION set_potential_x_min
  
  
  
  FUNCTION set_potential_x_max()

    REAL(num) :: set_potential_x_max
    REAL(num) :: pot_x_max, source
    TYPE(potential_block), POINTER :: current
  
    pot_x_max = 0._num
    
    IF (x_max_boundary_open) THEN
      IF (add_potential_source(c_bd_x_max)) THEN
        current => potential_list_x_max
          DO WHILE(ASSOCIATED(current))
            ! evaluate the temporal evolution of the potential source
            IF (time >= current%t_start .AND. time <= current%t_end) THEN
              IF (current%use_profile_function) &
                CALL potential_update_spatial_profile(current)
              source = potential_time_profile(current) * current%amp
              source = source * current%profile
              pot_x_max = pot_x_max + source
            END IF
            current => current%next
          END DO
      END IF
    END IF
    
    set_potential_x_max = pot_x_max
   
  END FUNCTION set_potential_x_max



  SUBROUTINE min_init_electrostatic

    es_dt_fact = 0.1_num
    NULLIFY(potential_list_x_min)
    NULLIFY(potential_list_x_max)

  END SUBROUTINE min_init_electrostatic



  SUBROUTINE setup_electrostatic

    IF (x_min_boundary .AND. bc_field(c_bd_x_min) /= c_bc_periodic ) THEN
      x_min_boundary_open = .TRUE.
    ELSE
      x_min_boundary_open = .FALSE.
    END IF
    IF (x_max_boundary .AND. bc_field(c_bd_x_max) /= c_bc_periodic ) THEN
      x_max_boundary_open = .TRUE.
    ELSE
      x_max_boundary_open = .FALSE.
    END IF

    wcharge_min = 0._num
    wcharge_max = 0._num
    dwcharge_min = 0._num
    dwcharge_max = 0._num
    convect_curr_min = 0._num
    convect_curr_max = 0._num

    es_current = 0._num

  END SUBROUTINE setup_electrostatic



  SUBROUTINE es_wall_current_diagnostic

    IF (x_min_boundary_open) THEN
      es_current(1) = (dwcharge_min - convect_curr_min) / dt
    END IF

    IF (x_max_boundary_open) THEN
      es_current(nx) = (dwcharge_max - convect_curr_max) / dt
    END IF

  END SUBROUTINE es_wall_current_diagnostic



  SUBROUTINE set_ey
  
    IF (ey_profile%use_profile_function) THEN
      CALL efield_update_profile(ey_profile, ey)
    END IF

  END SUBROUTINE set_ey



  SUBROUTINE set_ez

    IF (ez_profile%use_profile_function) THEN
      CALL efield_update_profile(ez_profile, ez)
    END IF

  END SUBROUTINE set_ez



  SUBROUTINE efield_update_profile(e_block, e_array)

    TYPE(efield_block), INTENT(INOUT) :: e_block
    REAL(num), DIMENSION(1-ng:nx+ng), INTENT(OUT) :: e_array
    INTEGER :: err, i
    TYPE(parameter_pack) :: parameters

    err = 0
    DO i = 1-ng, nx+ng
      parameters%pack_ix = i
      e_array(i) = evaluate_with_parameters(e_block%profile_function, &
        parameters, err)
    END DO

  END SUBROUTINE efield_update_profile


  SUBROUTINE setup_petsc_variables(ncells_local, ncells_global)

    INTEGER, INTENT(IN) :: ncells_local, ncells_global

    CALL setup_petsc_vector(ncells_local, ncells_global)
    CALL setup_petsc_matrix(ncells_local, ncells_global)
    CALL setup_petsc_ksp

  END SUBROUTINE setup_petsc_variables



  SUBROUTINE initialize_petsc(communicator)

    INTEGER, INTENT(IN) :: communicator

    PETSC_COMM_WORLD = communicator
    CALL PetscInitialize('./src/housekeeping/petsc_runtime_options.opt', perr)

  END SUBROUTINE initialize_petsc


  SUBROUTINE finalize_petsc

    CALL destroy_petsc
    CALL PetscFinalize(perr)

  END SUBROUTINE finalize_petsc


  SUBROUTINE finalize_electrostatic_solver

    CALL deallocate_potentials
    CALL deallocate_efield_profiles
    CALL finalize_petsc

  END SUBROUTINE finalize_electrostatic_solver

END MODULE electrostatic
