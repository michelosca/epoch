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

    LOGICAL :: source_on

    CHARACTER(LEN=string_length) :: e_name
    INTEGER :: id_electrons

  END TYPE inductive_block 

  TYPE(inductive_block), POINTER :: inductive_source 
  LOGICAL :: inductive_heating_flag = .FALSE.

  PUBLIC :: inductive_heating_y, inductive_heating_prepare_source
  PUBLIC :: inductive_heating_flag, inductive_source
  PUBLIC :: init_inductive_source_block

CONTAINS

  SUBROUTINE inductive_heating_y
    CALL update_total_current(inductive_source)
    CALL update_sum_vy_e(inductive_source)
    CALL update_dEy_dt(inductive_source)
  END SUBROUTINE inductive_heating_y



  SUBROUTINE init_inductive_source_block(inductive_source_block)

    TYPE(inductive_block), POINTER :: inductive_source_block
    
    inductive_source_block%x_min = 0._num
    inductive_source_block%x_max = 0._num 
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

  END SUBROUTINE update_total_current 



  SUBROUTINE update_sum_vy_e(inductive_source_block)

    TYPE(inductive_block), POINTER :: inductive_source_block
    REAL(num) :: sendbuf, recvbuf
    INTEGER :: ierror

    sendbuf = inductive_source_block%sum_vy_e
    CALL MPI_ALLREDUCE(sendbuf, recvbuf, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
      comm, ierror) 
    inductive_source_block%sum_vy_e = recvbuf
    
  END SUBROUTINE update_sum_vy_e

  

  SUBROUTINE update_dEy_dt(inductive_source_block)

    TYPE(inductive_block), POINTER :: inductive_source_block
    REAL(num) :: j_total, j_cond_fac, sum_vy_e

    j_total = inductive_source_block%j_total
    j_cond_fac = inductive_source_block%j_cond_fac  ! = electron_weight * e / L_source
    sum_vy_e = inductive_source_block%sum_vy_e
    inductive_source_block%dEy_dt = (j_total - j_cond_fac*sum_vy_e) / epsilon0


  END SUBROUTINE update_dEy_dt

  

  SUBROUTINE inductive_heating_prepare_source
    IF ((time >= inductive_source%t_start) .AND. &
        (time <= inductive_source%t_end)) THEN
      inductive_source%source_on = .TRUE.
    ELSE
      inductive_source%source_on = .TRUE.
    END IF
    inductive_source%sum_vy_e = 0._num
  END SUBROUTINE inductive_heating_prepare_source

#endif
END MODULE inductive_heating