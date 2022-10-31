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
!------------------------------------------------------------------------------
! Inductive heating
!   - Based on: A. Meige et al., Phys. Plasmas 12, 052317 (2005)
!                https://doi.org/10.1063/1.1897390
!   - Written by M. Osca Engelbrecht
!------------------------------------------------------------------------------

MODULE inductive_heating 
#ifdef ELECTROSTATIC
  USE shared_data
  USE evaluator
  USE strings_advanced
  
  IMPLICIT NONE

  PRIVATE

  TYPE inductive_heating_block 

    ! Input parameters
    INTEGER :: id
    TYPE(primitive_stack) :: time_function
    REAL(num) :: j0_amp, t_start, t_end
    REAL(num) :: x_min, x_max

    ! Simulation paramters
    REAL(num) :: j_total, j_cond_fac
    REAL(num) :: dEy_dt, sum_vy_e
    INTEGER :: ix_min, ix_max

    LOGICAL :: source_on

    CHARACTER(LEN=string_length) :: e_name
    INTEGER :: id_electrons

  END TYPE inductive_heating_block 

  TYPE(inductive_heating_block), DIMENSION(:), POINTER :: inductive_sources 
  INTEGER, SAVE :: n_heating_sources = 0
  LOGICAL :: inductive_heating_flag = .FALSE.

  PUBLIC :: inductive_heating_y, inductive_heating_prepare_sources
  PUBLIC :: inductive_heating_init_inductive_block
  PUBLIC :: inductive_heating_allocate_sources_list
  PUBLIC :: inductive_heating_setup, inductive_heating_add_source_to_list
  PUBLIC :: inductive_heating_add_particle_velocity
  PUBLIC :: inductive_heating_flag, inductive_heating_block

CONTAINS

  SUBROUTINE inductive_heating_y

    INTEGER :: i
    TYPE(inductive_heating_block), POINTER :: inductive_source

    IF (inductive_heating_flag) THEN
      DO i = 1, n_heating_sources
        inductive_source => inductive_sources(i)
        CALL update_total_current(inductive_source)
        CALL update_sum_vy_e(inductive_source)
        CALL update_dEy_dt(inductive_source)
      END DO
      CALL update_ey
    END IF
  
  END SUBROUTINE inductive_heating_y



  SUBROUTINE update_total_current(inductive_source_block)

    TYPE(inductive_heating_block),POINTER,INTENT(INOUT)::inductive_source_block
    REAL(num) :: j0_amp
    INTEGER :: err
    TYPE(parameter_pack) :: parameters

    err = 0
    IF (inductive_source_block%source_on) THEN
      j0_amp = inductive_source_block%j0_amp
      inductive_source_block%j_total = j0_amp * evaluate_with_parameters( &
        inductive_source_block%time_function, parameters, err)
    ELSE
      inductive_source_block%j_total = 0._num
    END IF

  END SUBROUTINE update_total_current 



  SUBROUTINE update_sum_vy_e(inductive_source_block)

    TYPE(inductive_heating_block),POINTER,INTENT(INOUT)::inductive_source_block
    REAL(num) :: sendbuf, recvbuf
    INTEGER :: ierror

    sendbuf = inductive_source_block%sum_vy_e
    CALL MPI_ALLREDUCE(sendbuf, recvbuf, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
      comm, ierror) 
    inductive_source_block%sum_vy_e = recvbuf
    
  END SUBROUTINE update_sum_vy_e

  

  SUBROUTINE update_dEy_dt(inductive_source_block)

    TYPE(inductive_heating_block),POINTER,INTENT(INOUT)::inductive_source_block
    REAL(num) :: j_total, j_cond_fac, sum_vy_e

    j_total = inductive_source_block%j_total
    j_cond_fac = inductive_source_block%j_cond_fac  ! = electron_weight * e / L_source
    sum_vy_e = inductive_source_block%sum_vy_e
    inductive_source_block%dEy_dt = (j_total - j_cond_fac*sum_vy_e) / epsilon0

  END SUBROUTINE update_dEy_dt



  SUBROUTINE update_ey

    TYPE(inductive_heating_block), POINTER :: inductive_source
    REAL(num) :: dEy_dt
    INTEGER :: i, j, ix_min, ix_max


    DO j = 1, n_heating_sources
      inductive_source => inductive_sources(j)
      dEy_dt = inductive_source%dEy_dt
      ix_min = inductive_source%ix_min
      ix_max = inductive_source%ix_max
      DO i = ix_min, ix_max 
        ey(i) = ey(i) + dEy_dt * dt
      END DO
    END DO

  END SUBROUTINE update_ey

  

  SUBROUTINE inductive_heating_allocate_sources_list(count)
    INTEGER, INTENT(IN) :: count

    ALLOCATE(inductive_sources(count))

  END SUBROUTINE inductive_heating_allocate_sources_list



  SUBROUTINE inductive_heating_init_inductive_block(inductive_source_block)

    TYPE(inductive_heating_block), POINTER :: inductive_source_block
    
    n_heating_sources = n_heating_sources + 1

    inductive_source_block%id = n_heating_sources
    inductive_source_block%x_min = 0._num
    inductive_source_block%x_max = 0._num 
    inductive_source_block%ix_min = -100000000
    inductive_source_block%ix_max = -100000000 
    inductive_source_block%t_start = 0._num
    inductive_source_block%t_end = HUGE(0._num) 
    
    inductive_source_block%j0_amp = 0._num
    inductive_source_block%j_total = 0._num
    
    inductive_source_block%j_cond_fac = 0._num
    inductive_source_block%sum_vy_e = 0._num
    
    inductive_source_block%dEy_dt = 0._num
    
    inductive_source_block%source_on = .FALSE.
    inductive_source_block%e_name = "None"
    inductive_source_block%id_electrons = -1
      
  END SUBROUTINE inductive_heating_init_inductive_block 



  SUBROUTINE inductive_heating_add_source_to_list(inductive_source_block)

    TYPE(inductive_heating_block),POINTER,INTENT(IN) :: inductive_source_block
    INTEGER :: id

    id = inductive_source_block%id
    inductive_sources(id) = inductive_source_block

  END SUBROUTINE inductive_heating_add_source_to_list



  SUBROUTINE inductive_heating_prepare_sources

    INTEGER :: i
    TYPE(inductive_heating_block), POINTER :: inductive_source

    DO i = 1, n_heating_sources
      inductive_source => inductive_sources(i)
      IF ((time >= inductive_source%t_start) .AND. &
          (time <= inductive_source%t_end)) THEN
        inductive_source%source_on = .TRUE.
      ELSE
        inductive_source%source_on = .FALSE.
      END IF
      inductive_source%sum_vy_e = 0._num
    END DO

  END SUBROUTINE inductive_heating_prepare_sources



  SUBROUTINE inductive_heating_add_particle_velocity(ispecies, part_pos, vy)

    INTEGER, INTENT(IN) :: ispecies
    REAL(num), INTENT(IN) :: vy, part_pos
    INTEGER :: i
    TYPE(inductive_heating_block), POINTER :: inductive_source
    
    DO i = 1, n_heating_sources
      inductive_source => inductive_sources(i)
      IF (inductive_source%source_on .AND. &
        inductive_source%id_electrons == ispecies) THEN
        IF ( (part_pos >= inductive_source%x_min) .AND. &
          (part_pos <= inductive_source%x_max) ) THEN
          inductive_source%sum_vy_e = inductive_source%sum_vy_e + vy
        END IF
      END IF
    END DO

  END SUBROUTINE inductive_heating_add_particle_velocity



  SUBROUTINE inductive_heating_setup 

    INTEGER :: ispecies, i, j
    REAL(num) :: x_min, x_max, L_source, electron_weight
    REAL(num) :: x_cell_left
    CHARACTER(LEN=string_length) :: e_name, species_name
    TYPE(inductive_heating_block), POINTER :: inductive_source

    DO j = 1, n_heating_sources
      inductive_source => inductive_sources(j)
      IF (inductive_heating_flag) THEN 
          ! Set electron species parameters
          e_name = TRIM(ADJUSTL(inductive_source%e_name))
          DO ispecies = 1, n_species
            species_name = TRIM(ADJUSTL(species_list(ispecies)%name)) 
            IF (str_cmp(species_name, e_name)) THEN
              inductive_source%id_electrons = ispecies
              electron_weight = species_list(ispecies)%weight
              EXIT 
            END IF
          END DO
          
          ! Set conduction current factor
          x_max = inductive_source%x_max
          x_min = inductive_source%x_min
          L_source = x_max - x_min
          inductive_source%j_cond_fac = -q0 / L_source * electron_weight 

          ! Set source start cell
          IF (x_min > x_max_local) THEN
            inductive_source%ix_min = nx+ng+1
          ELSE
            DO i = 1-ng, nx+ng
              x_cell_left = i * dx + x_min_local
              IF (x_cell_left >= x_min) THEN
                inductive_source%ix_min = i
                EXIT
              END IF
            END DO
          END IF

          ! Set source end cell
          IF (x_max < x_min_local) THEN
            inductive_source%ix_max = 1-ng-1
          ELSE
            DO i = nx+ng,1-ng,-1
              x_cell_left = i *dx + x_min_local
              IF (x_cell_left <= x_max) THEN
                inductive_source%ix_max = i
                EXIT
              END IF
            END DO
          END IF
      END IF
    END DO
  
  END SUBROUTINE inductive_heating_setup 

#endif
END MODULE inductive_heating