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
#ifdef ELECTROSTATIC
#ifndef TRIDIAG
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
  USE petscvec
  USE petscmat
  USE petscksp
#endif
  USE shared_data
  USE custom_electrostatic
  USE evaluator
  USE boundary
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: es_update_e_field, setup_electrostatic
  PUBLIC :: init_potential, es_wall_current_diagnostic
  PUBLIC :: attach_potential, update_cell_count
  PUBLIC :: potential_update_spatial_profile, min_init_electrostatic
  PUBLIC :: finalize_electrostatic_solver
#ifndef TRIDIAG
  PUBLIC :: initialize_petsc, setup_petsc_variables, finalize_petsc

  Vec, SAVE :: es_potential_vec, charge_dens_vec
  Mat, SAVE :: transform_mtrx
  MatNullSpace, SAVE :: matnull
  PetscErrorCode, SAVE :: perr
  KSP, SAVE :: ksp
#else
  INTEGER, ALLOCATABLE, DIMENSION(:) :: cell_start_each_rank
#endif
  INTEGER :: nx_start, nx_end, nx_all
  ! Wall charge surface density
  REAL(num) :: wcharge_min_prev, wcharge_min_now
  REAL(num) :: wcharge_max_prev, wcharge_max_now
  REAL(num) :: wcharge_min_diff, wcharge_max_diff
  REAL(num) :: Q_min_now, Q_min_prev
  REAL(num) :: Q_max_now, Q_max_prev
  REAL(num) :: Q_conv_max, Q_conv_min
  REAL(num) :: pot_ext_max, pot_ext_min


CONTAINS

  SUBROUTINE es_update_e_field
    
    REAL(num), DIMENSION(:), ALLOCATABLE :: es_charge_density
    REAL(num) :: rho_max, rho_min

    ! Charge weighting from particles to the grid, i.e. charge density
    ALLOCATE(es_charge_density(1-ng:nx+ng))
    CALL es_calc_charge_density(es_charge_density)

    ! This subroutine updates the charge surface density values at wall
    CALL es_calc_charge_density_at_wall

    ! Charge density to electrostatic potential
    !  - This subroutine calculates the electric potential on es_potential
    CALL es_calc_potential(es_charge_density(nx_start:nx_end))

    ! Save charge density at dx/2 from boundary and deallocate
    IF (x_min_boundary_open) rho_min = SUM(es_charge_density(0:1)) * 0.5_num
    IF (x_max_boundary_open) rho_max = SUM(es_charge_density(nx-1:nx)) * 0.5_num
    DEALLOCATE(es_charge_density)

    ! Calculate electric field in x-direction
    CALL es_calc_ex(rho_min, rho_max)

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

    ! Same as in particle_head.inc but without integers defined for electrostatic
    REAL(num), DIMENSION(sf_min:sf_max) :: gx
    REAL(num) :: cell_x_r, cell_frac_x
    INTEGER :: cell_x
#ifndef PARTICLE_SHAPE_TOPHAT
    REAL(num) :: cx2
#endif
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num), PARAMETER :: third = 1.0_num / 3.0_num
    REAL(num), PARAMETER :: fac1 = 0.125_num * third
    REAL(num), PARAMETER :: fac2 = 0.5_num * third
    REAL(num), PARAMETER :: fac3 = 7.1875_num * third
#endif

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
  
    REAL(num), DIMENSION(nx_start:nx_end), INTENT(IN) :: charge_density
    
    CALL poisson_solver(charge_density(nx_start:nx_end))
    
    ! Potential boundaries (between MPI processors, valid for periodic bc)
    CALL field_bc(es_potential, ng)

    ! Potential (open) boundary conditions
    IF (x_min_boundary) THEN
      IF (capacitor_min) THEN
        es_potential(0) = 0._num
      ELSE IF (.NOT.capacitor_flag) THEN
        es_potential(0) = set_potential_x_min()
      END IF
    END IF
    IF (x_max_boundary) THEN
      IF (capacitor_max) THEN
        es_potential(nx) = 0._num
      ELSE IF (.NOT.capacitor_flag) THEN
        es_potential(nx) = set_potential_x_max()
      END IF
    END IF

  END SUBROUTINE es_calc_potential


#ifdef TRIDIAG
  SUBROUTINE poisson_solver(rho)

    ! This subroutine solves Poisson's equation
    ! Input: charge_density array and the length of the data array
    ! Output: electric potential array (data_array)
    REAL(num), DIMENSION(nx_start:nx_end), INTENT(IN) :: rho
    REAL(num), DIMENSION(nx_start:nx_end) :: solver_rho
    INTEGER, PARAMETER :: dp = selected_real_kind(30,300)
    INTEGER :: i
    REAL(num), ALLOCATABLE, DIMENSION(:) :: rho_all, pot_all
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: b
    REAL(dp) :: w
    REAL(num) :: fac

    ! Get charge density array ready for solver
    fac = -dx*dx/epsilon0
    solver_rho = rho * fac
    IF (x_min_boundary) THEN
        CALL set_poissonsolver_min_bc(rho(nx_start), solver_rho(nx_start))
    END IF
    IF (x_max_boundary) THEN
        CALL set_poissonsolver_max_bc(rho(nx_end), solver_rho(nx_end))
    END IF

    IF (rank==0) THEN
      ALLOCATE(rho_all(nx_start:nx_all), pot_all(nx_start:nx_all))
      IF (nproc > 1) THEN
        CALL MPI_GATHERV(solver_rho(nx_start), nx_each_rank(rank+1), &
          MPI_DOUBLE_PRECISION, rho_all(nx_start), nx_each_rank, &
          cell_start_each_rank, MPI_DOUBLE_PRECISION, 0, comm, errcode)
      ELSE
        rho_all = solver_rho
      END IF

      ! Tridiagonal solver method
      ALLOCATE(b(nx_start:nx_all))
      b = -2._dp
      IF (capacitor_max) b(nx_start) = -1._dp - REAL(dx*capacitor/epsilon0,dp)
      IF (capacitor_min) b(nx_all) = -1._dp - REAL(dx*capacitor/epsilon0,dp)

      DO i = nx_start+1, nx_all
        w = 1._dp / b(i-1)
        b(i) = b(i) - w
        rho_all(i) = rho_all(i) - REAL(w,num) * rho_all(i-1)
      END DO

      pot_all(nx_all) = rho_all(nx_all) / REAL(b(nx_all),num)
      DO i = nx_all-1, nx_start, -1
        pot_all(i) = (rho_all(i) - pot_all(i+1)) / REAL(b(i),num)
      END DO

      CALL MPI_SCATTERV(pot_all(nx_start), nx_each_rank, cell_start_each_rank, &
        MPI_DOUBLE_PRECISION, es_potential(nx_start), nx_each_rank(rank+1), &
        MPI_DOUBLE_PRECISION, 0, comm, errcode)
      DEALLOCATE(rho_all, pot_all, b)

    ELSE
      ALLOCATE(rho_all(0), pot_all(0))
      CALL MPI_GATHERV(solver_rho(nx_start), nx_each_rank(rank+1), &
        MPI_DOUBLE_PRECISION, &
        rho_all, nx_each_rank, cell_start_each_rank, &
        MPI_DOUBLE_PRECISION, 0, comm, errcode)

      CALL MPI_SCATTERV(pot_all, nx_each_rank, cell_start_each_rank, &
        MPI_DOUBLE_PRECISION, es_potential(nx_start), nx_each_rank(rank+1), &
        MPI_DOUBLE_PRECISION, 0, comm, errcode)
      DEALLOCATE(rho_all, pot_all)
    END IF

  END SUBROUTINE poisson_solver

#else

  SUBROUTINE poisson_solver(rho)

    ! This subroutine solves Poisson's equation
    ! Input: charge_density array
    ! Output: electric potential array (data_array)
    REAL(num), DIMENSION(nx_start:nx_end), INTENT(IN) :: rho
    REAL(num), ALLOCATABLE, DIMENSION(:) :: solver_rho
    REAL(num) :: fac
    INTEGER :: i, i_start, i_end, n0
    REAL(num), DIMENSION(:), POINTER :: vec_pointer

    IF (capacitor_max .AND. x_min_boundary) THEN
      n0 = 1
    ELSE
      n0 = 0
    END IF
    i_start = nx_start + n0
    i_end = nx_end + n0

    ALLOCATE(solver_rho(i_start:i_end))

    ! Get array ready for solver
    fac = -dx*dx/epsilon0
    solver_rho(i_start:i_end) = rho(nx_start:nx_end) * fac

    ! Set up boundaries
    IF (x_min_boundary) THEN
      CALL set_poissonsolver_min_bc(rho(nx_start), solver_rho(i_start))
    END IF
    IF (x_max_boundary) THEN
      CALL set_poissonsolver_max_bc(rho(nx_end), solver_rho(i_end))
    END IF

    ! Get charge density data (data_array) into the Petsc vector
    CALL VecGetArrayF90(charge_dens_vec, vec_pointer, perr)

    DO i = i_start, i_end
      vec_pointer(i) = solver_rho(i)
    END DO
    CALL VecRestoreArrayF90(charge_dens_vec, vec_pointer, perr)
    DEALLOCATE(solver_rho)

    ! Solve linear system: transform_mtrx*es_potential_vec = charge_dens_vec
    CALL KSPSolve(ksp, charge_dens_vec, es_potential_vec, perr)

    ! Pass electric potential data from PETSc to Fortran
    CALL VecGetArrayReadF90(es_potential_vec, vec_pointer, perr)
    DO i = i_start, i_end
      es_potential(i-n0) = vec_pointer(i)
    END DO
    CALL VecRestoreArrayReadF90(es_potential_vec, vec_pointer, perr)

  END SUBROUTINE poisson_solver



  SUBROUTINE setup_petsc_vector
    
    INTEGER :: n_local, n_global

    n_local = nx_end - nx_start + 1
    n_global = nx_all

    ! Charge density vector
    CALL VecCreateMPI(comm, n_local, n_global, charge_dens_vec, perr)
    CALL VecSetFromOptions(charge_dens_vec, perr)
    ! Electric potential vector
    CALL VecDuplicate(charge_dens_vec, es_potential_vec, perr)

  END SUBROUTINE setup_petsc_vector



  SUBROUTINE setup_petsc_matrix
    
    INTEGER :: n_first, n_last, n_local, n_global
    PetscInt :: zero
    PetscScalar :: values(3)
    PetscInt :: col(3), row

    n_local = nx_end - nx_start + 1
    n_global = nx_all

    ! Petsc subroutines setting up the Matrix
    CALL MatCreate(comm, transform_mtrx, perr)
    CALL MatSetSizes(transform_mtrx, n_local, n_local, n_global, n_global, perr)
    CALL MatSetFromOptions(transform_mtrx, perr)
    CALL MatSetUp(transform_mtrx, perr)

    ! Petsc matrix indexes (global). First index in Petsc is zero, hence
    ! Shift cell indexing so that it starts at 0
    n_first = nx_global_min - 1
    n_last = nx_global_max - 1
    IF (capacitor_max) n_last = nx_global_max

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
      IF (capacitor_max) values(2) = -1._num - dx*capacitor/epsilon0
      !MatSetValues(matrix, #rows,rows-global-indexes, 
      !  #columns, col-global-indexes, insertion_values, insert/add,err )
      CALL MatSetValues(transform_mtrx, 1, row, 2, col(1:2), values(2:3), &
              INSERT_VALUES, perr)
      IF (capacitor_max) values(2) = -2._num
      ! Move to next row to be set
      n_first = n_first + 1
    END IF

    IF (x_max_boundary) THEN
      IF (.NOT.capacitor_min) THEN
        n_last = n_last - 1
      END IF
      ! Bottom matrix row
      row = n_last
      col(1) = row-1
      col(2) = row
      IF (capacitor_min) values(2) = -1._num - dx*capacitor/epsilon0
      !MatSetValues(matrix, #rows,rows-global-indexes, 
      !  #columns, col-global-indexes, insertion_values, insert/add,err )
      CALL MatSetValues(transform_mtrx, 1, row, 2, col(1:2), values(1:2), &
              INSERT_VALUES, perr)
      IF (capacitor_min) values(2) = -2._num
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
#endif


  SUBROUTINE set_poissonsolver_min_bc(rho_min, solver_rho_min)

    REAL(num), INTENT(IN) :: rho_min
    REAL(num), INTENT(INOUT) :: solver_rho_min
    REAL(num) :: term0, term1, term2
    REAL(num) :: fac

    IF (capacitor_max) THEN
      fac = -dx*dx/epsilon0
      term0 = rho_min * 0.5_num
      term1 = wcharge_min_prev / dx
      term2 = Q_conv_min - Q_min_prev + pot_ext_min * capacitor
      solver_rho_min = term0 + term1 + term2 / dx
      solver_rho_min = solver_rho_min * fac
    ELSE IF (.NOT.capacitor_flag) THEN
      solver_rho_min = solver_rho_min - set_potential_x_min()
    END IF

  END SUBROUTINE set_poissonsolver_min_bc


  SUBROUTINE set_poissonsolver_max_bc(rho_max, solver_rho_max)

    REAL(num), INTENT(IN) :: rho_max
    REAL(num), INTENT(INOUT) :: solver_rho_max
    REAL(num) :: term0, term1, term2
    REAL(num) :: fac

    IF (capacitor_min) THEN
      fac = -dx*dx/epsilon0
      term0 = rho_max * 0.5_num
      term1 = wcharge_max_prev / dx
      term2 = Q_conv_max - Q_max_prev + pot_ext_max * capacitor
      solver_rho_max = term0 + term1 + term2 / dx
      solver_rho_max = solver_rho_max * fac
    ELSE IF (.NOT.capacitor_flag) THEN
      solver_rho_max = solver_rho_max - set_potential_x_max()
    END IF

  END SUBROUTINE set_poissonsolver_max_bc


  SUBROUTINE es_calc_ex(rho_min, rho_max)

    REAL(num), INTENT(IN) :: rho_min, rho_max
    INTEGER :: ix
    REAL(num) :: idx

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

    IF (x_min_boundary_open) THEN
      !Update E-field at wall
      ex(0) = ex(1) - rho_min * dx / epsilon0
    END IF

    IF (x_max_boundary_open) THEN
      !Update E-field at wall
      ex(nx) = ex(nx-1) + rho_max * dx / epsilon0
    END IF

    DO ix = 1, 2*c_ndims
      CALL field_clamp_wall(ex, ng, c_stagger_ex, ix)
    END DO

  END SUBROUTINE es_calc_ex



  SUBROUTINE es_calc_charge_density_at_wall

    ! This subroutine only makes sense for open boundary conditions
    IF (x_min_boundary_open) THEN
      ! Buffer previous values
      Q_min_prev = Q_min_now
      wcharge_min_prev = wcharge_min_now

      ! Calculate new surface charge density
      pot_ext_min = set_potential_x_min()
      Q_min_now = capacitor * (pot_ext_min - es_potential(0))
      Q_conv_min = convect_curr_min
      wcharge_min_now = wcharge_min_prev + Q_conv_min + Q_min_now - Q_min_prev

      ! Charge wall difference
      wcharge_min_diff = wcharge_min_now - wcharge_min_prev
    END IF

    IF (x_max_boundary_open) THEN
      ! Buffer previous values
      Q_max_prev = Q_max_now
      wcharge_max_prev = wcharge_max_now

      ! Calculate new surface charge density
      pot_ext_max = set_potential_x_max()
      Q_max_now = capacitor * (pot_ext_max + es_potential(nx_end))
      Q_conv_max = convect_curr_max
      wcharge_max_now = wcharge_max_prev + Q_conv_max + Q_max_now - Q_max_prev

      ! Charge wall difference
      wcharge_max_diff = wcharge_max_now - wcharge_max_prev
    END IF

  END SUBROUTINE es_calc_charge_density_at_wall



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
    capacitor_flag = .FALSE.
    capacitor_min = .FALSE.
    capacitor_max = .FALSE.
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

    wcharge_min_now = 0._num
    wcharge_max_prev = 0._num
    wcharge_min_now = 0._num
    wcharge_max_prev = 0._num
    wcharge_min_diff = 0._num
    wcharge_max_diff = 0._num
    convect_curr_min = 0._num
    convect_curr_max = 0._num
    Q_min_now = 0._num
    Q_min_prev = 0._num
    Q_max_now = 0._num
    Q_max_prev = 0._num
    Q_conv_min = 0._num
    Q_conv_max = 0._num
    pot_ext_max = 0._num
    pot_ext_min = 0._num

    es_current = 0._num

#ifdef TRIDIAG
    ALLOCATE(nx_each_rank(nproc), cell_start_each_rank(nproc))
    nx_each_rank = cell_x_max - cell_x_min + 1
    IF (capacitor_flag) THEN
      IF (capacitor_max) THEN
        cell_start_each_rank = cell_x_min
        cell_start_each_rank(1) = cell_x_min(1) - 1
        nx_each_rank(1) = nx_each_rank(1) + 1
        nx_each_rank(nproc) = nx_each_rank(nproc) - 1
      END IF
      IF (capacitor_min) THEN
        cell_start_each_rank = cell_x_min - 1
      END IF
    ELSE
      nx_each_rank(nproc) = nx_each_rank(nproc) - 1
      cell_start_each_rank = cell_x_min - 1
    END IF
#endif

    CALL update_cell_count(nx)

  END SUBROUTINE setup_electrostatic



  SUBROUTINE update_cell_count(nx)

    INTEGER, INTENT(IN) :: nx

    ! Number of cells to be solved on each processors
    nx_start = 1
    nx_end = nx
    IF (x_min_boundary) THEN
      IF (capacitor_max) nx_start = 0
    END IF
    IF (x_max_boundary) THEN
      IF (.NOT.capacitor_min) nx_end = nx_end - 1
    END IF

#ifdef TRIDIAG
    IF (capacitor_min) THEN
#else
    IF (capacitor_flag) THEN
#endif
      nx_all = nx_global
    ELSE
      nx_all = nx_global - 1
    END IF

  END SUBROUTINE update_cell_count



  SUBROUTINE es_wall_current_diagnostic

    IF (x_min_boundary_open) THEN
      es_current(1) = (wcharge_min_diff - convect_curr_min) / dt
      convect_curr_min = 0._num
    END IF

    IF (x_max_boundary_open) THEN
      es_current(nx) = (wcharge_max_diff - convect_curr_max) / dt
      convect_curr_max = 0._num
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



  SUBROUTINE finalize_electrostatic_solver

#ifdef TRIDIAG
    DEALLOCATE(nx_each_rank, cell_start_each_rank)
#endif
    CALL deallocate_potentials
    CALL deallocate_efield_profiles
#ifndef TRIDIAG
    CALL finalize_petsc
#endif
  END SUBROUTINE finalize_electrostatic_solver


#ifndef TRIDIAG
  SUBROUTINE setup_petsc_variables(nx)

    INTEGER, INTENT(IN) :: nx

    CALL update_cell_count(nx)
    CALL setup_petsc_vector
    CALL setup_petsc_matrix
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
#endif
#endif
END MODULE electrostatic
