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

  TYPE inductive_block 

    ! Input parameters
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

  END TYPE inductive_block 

  TYPE(inductive_block), POINTER :: inductive_source 
  LOGICAL :: inductive_heating_flag = .FALSE.

  PUBLIC :: inductive_heating_y, inductive_heating_prepare_source
  PUBLIC :: inductive_heating_flag, inductive_source
  PUBLIC :: init_inductive_source_block
  PUBLIC :: setup_inductive_heating

CONTAINS

  SUBROUTINE inductive_heating_y
    CALL update_total_current(inductive_source)
    CALL update_sum_vy_e(inductive_source)
    CALL update_dEy_dt(inductive_source)
    CALL update_ey(inductive_source)
  END SUBROUTINE inductive_heating_y



  SUBROUTINE init_inductive_source_block(inductive_source_block)

    TYPE(inductive_block), POINTER :: inductive_source_block
    
    ALLOCATE(inductive_source_block)

    inductive_source_block%x_min = 0._num
    inductive_source_block%x_max = 0._num 
    inductive_source_block%ix_min = -1
    inductive_source_block%ix_max = -1 
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
      
  END SUBROUTINE init_inductive_source_block



  SUBROUTINE update_total_current(inductive_source_block)

    TYPE(inductive_block), POINTER :: inductive_source_block
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
    print*, rank,time,step, &
      'current_amplitude', inductive_source_block%j_total

  END SUBROUTINE update_total_current 



  SUBROUTINE update_sum_vy_e(inductive_source_block)

    TYPE(inductive_block), POINTER :: inductive_source_block
    REAL(num) :: sendbuf, recvbuf
    INTEGER :: ierror

    sendbuf = inductive_source_block%sum_vy_e
    CALL MPI_ALLREDUCE(sendbuf, recvbuf, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
      comm, ierror) 
    inductive_source_block%sum_vy_e = recvbuf
    print*, rank, time, step, &
      'vy_sum', inductive_source_block%sum_vy_e
    
  END SUBROUTINE update_sum_vy_e

  

  SUBROUTINE update_dEy_dt(inductive_source_block)

    TYPE(inductive_block), POINTER :: inductive_source_block
    REAL(num) :: j_total, j_cond_fac, sum_vy_e

    j_total = inductive_source_block%j_total
    j_cond_fac = inductive_source_block%j_cond_fac  ! = electron_weight * e / L_source
    sum_vy_e = inductive_source_block%sum_vy_e
    inductive_source_block%dEy_dt = (j_total - j_cond_fac*sum_vy_e) / epsilon0
    print*, rank, time, step, &
      ' dEydt ', inductive_source_block%dEy_dt


  END SUBROUTINE update_dEy_dt



  SUBROUTINE update_ey(inductive_source_block)

    TYPE(inductive_block), POINTER :: inductive_source_block
    REAL(num) :: dEy_dt
    INTEGER :: i, ix_min, ix_max

    dEy_dt = inductive_source_block%dEy_dt
    ix_min = inductive_source%ix_min
    ix_max = inductive_source%ix_max
    DO i = ix_min, ix_max 
      ey(i) = ey(i) + dEy_dt * dt
      print*, rank, time, step, 'ey', i*dx + x_min_local, ey(i)
    END DO


  END SUBROUTINE update_ey

  

  SUBROUTINE inductive_heating_prepare_source
    IF ((time >= inductive_source%t_start) .AND. &
        (time <= inductive_source%t_end)) THEN
      inductive_source%source_on = .TRUE.
    ELSE
      inductive_source%source_on = .TRUE.
    END IF
    inductive_source%sum_vy_e = 0._num
  END SUBROUTINE inductive_heating_prepare_source

  SUBROUTINE setup_inductive_heating

    INTEGER :: ispecies, i
    REAL(num) :: x_min, x_max, L_source, electron_weight
    REAL(num) :: x_cell_left
    CHARACTER(LEN=string_length) :: e_name, species_name

    IF (deck_state /= c_ds_first .AND. inductive_heating_flag) THEN 
        ! Set electron species parameters
        e_name = TRIM(ADJUSTL(inductive_source%e_name))
        DO ispecies = 1, n_species
          species_name = TRIM(ADJUSTL(species_list(ispecies)%name)) 
          IF (str_cmp(species_name, e_name)) THEN
            inductive_source%id_electrons = ispecies
            electron_weight = species_list(inductive_source%id_electrons)%weight
            EXIT 
          END IF
        END DO
        
        ! Set conduction current factor
        x_max = inductive_source%x_max
        x_min = inductive_source%x_min
        L_source = x_max - x_min
        inductive_source%j_cond_fac = q0 / L_source * electron_weight 

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

  print*,'Input data'
  print*,rank, ' - current amplitude', inductive_source%j0_amp
  print*,rank, ' - t start and end', inductive_source%t_start,&
   inductive_source%t_end
  print*,rank, ' - x min and max', inductive_source%x_min, &
  inductive_source%x_max
  print*,rank, ' - cell min and max', inductive_source%ix_min, &
  inductive_source%ix_max
  print*,rank, ' - j_total ', inductive_source%j_total
  print*,rank, ' - j_cond_fac ', inductive_source%j_cond_fac
  print*,rank, ' - dEy/dt', inductive_source%dEy_dt
  print*,rank, ' - sum(v_y_e) ', inductive_source%sum_vy_e
  print*,rank, ' - source on/off', inductive_source%source_on
  print*,rank, ' - electron name', inductive_source%e_name
  print*,rank, ' - electron species', inductive_source%id_electrons

    END IF
  
  END SUBROUTINE setup_inductive_heating
#endif
END MODULE inductive_heating