
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

MODULE deck_psd_diagnostics_block

  USE strings_advanced
  USE shared_data
  
  IMPLICIT NONE
  SAVE
  
  PRIVATE
  PUBLIC :: psd_diagnostics_deck_initialise, psd_diagnostics_deck_finalise
  PUBLIC :: psd_diagnostics_block_start, psd_diagnostics_block_end
  PUBLIC :: psd_diagnostics_block_handle_element, psd_diagnostics_block_check
  
  INTEGER :: array_index
CONTAINS
  
  SUBROUTINE psd_diagnostics_deck_initialise
    
    LOGICAL :: max_boundary, min_boundary

    IF (deck_state /= c_ds_first) THEN
      ! Determine processor and cell
      IF (rank == 0) THEN
        min_boundary = psd_diag_xpos >= x_min_local
      ELSE
        min_boundary = psd_diag_xpos > x_min_local
      END IF
      max_boundary = psd_diag_xpos <= x_max_local

      IF (min_boundary .AND. max_boundary) THEN
        ! Processor involved
        psd_diag_rank = rank
        ! Cell where diagnostics take place
        psd_diag_cellx = FLOOR((psd_diag_xpos - x_min_local) / dx ) + 1

        ! Load PSD dump ID array
        IF (psd_diag_n_dumps > 0) THEN
          ALLOCATE(psd_diag_dump_id(psd_diag_n_dumps))
        END IF
      END IF
    END IF

  END SUBROUTINE psd_diagnostics_deck_initialise
  
  
  
  SUBROUTINE psd_diagnostics_deck_finalise

  END SUBROUTINE psd_diagnostics_deck_finalise
  
  
  
  SUBROUTINE psd_diagnostics_block_start
    
    IF (deck_state == c_ds_first) THEN
      psd_diag_n_dumps = 0
      av_over_sampling_period = .TRUE.
      psd_diag_period_is_dt = .FALSE.
      psd_diag_concatenated = .FALSE.
      psd_diag_print_setup = .FALSE.
      psd_diag_sampling_period = -1._num
      psd_diag_averaging_time = -1._num
      psd_diag_rank = -1
      psd_diag_cellx = -1
      psd_diag_xpos = HUGE(1._num)
    ELSE
      array_index = 1
    END IF
  
  END SUBROUTINE psd_diagnostics_block_start
  
  
  
  SUBROUTINE psd_diagnostics_block_end
  
  END SUBROUTINE psd_diagnostics_block_end
  
  
  
  FUNCTION psd_diagnostics_block_handle_element(element, value) RESULT(errcode)
  
    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    INTEGER :: ispecies 
    REAL(num) :: sampling_frequency

    errcode = c_err_none
    IF (element == blank) RETURN
    
    IF (deck_state == c_ds_first) THEN
      IF (str_cmp(element, 'sampling_period')) THEN
        psd_diag_sampling_period = as_real_print(value, element, errcode)
        ! By default averaging_time = sampling_period unless defined otherwise
        psd_diag_averaging_time = psd_diag_sampling_period

        ! If negative psd_diag_sampling_period, then = dt
        IF (psd_diag_sampling_period < 0._num) THEN
          psd_diag_period_is_dt = .TRUE.
        END IF

        RETURN
      END IF
      
      IF (str_cmp(element, 'averaging_time')) THEN
        psd_diag_averaging_time = as_real_print(value, element, errcode)
        RETURN
      END IF

      IF (str_cmp(element, 'sampling_frequency')) THEN
        sampling_frequency = as_real_print(value, element, errcode)
        psd_diag_sampling_period = 1._num / sampling_frequency
        RETURN
      END IF

      IF (str_cmp(element, 'position_x')) THEN
        psd_diag_xpos = as_real_print(value, element, errcode)
        RETURN
      END IF

      IF (str_cmp(element, 'average_over_sampling_period')) THEN
        av_over_sampling_period = as_logical_print(value, element, errcode)
        RETURN
      END IF

      IF (str_cmp(element, 'show_setup')) THEN
        psd_diag_print_setup = as_logical_print(value, element, errcode)
        RETURN
      END IF

      IF (str_cmp(element, 'number_density')) THEN
        psd_diag_n_dumps = psd_diag_n_dumps + n_species
      END IF
      IF (str_cmp(element, 'electric_potential')) THEN
        psd_diag_n_dumps = psd_diag_n_dumps + 1
      END IF
      IF (str_cmp(element, 'electric_field')) THEN
        psd_diag_n_dumps = psd_diag_n_dumps + 1
      END IF

    ELSE
      IF (.NOT.psd_diag_rank == rank) RETURN

      IF (str_cmp(element, 'number_density')) THEN
        DO ispecies = 1, n_species
          psd_diag_dump_id(array_index) = c_dump_number_density
          array_index = array_index + 1
        END DO

      ELSE IF (str_cmp(element, 'electric_potential')) THEN
        psd_diag_dump_id(array_index) = c_dump_es_potential 
        array_index = array_index + 1

      ELSE IF (str_cmp(element, 'electric_field')) THEN
        psd_diag_dump_id(array_index) = c_dump_ex
        array_index = array_index + 1
      END IF
    END IF

    
  END FUNCTION psd_diagnostics_block_handle_element
  
  
  
  FUNCTION psd_diagnostics_block_check() RESULT(errcode)
  
    INTEGER :: errcode

    errcode = c_err_none

    IF (psd_diag_n_dumps > 0 ) THEN
      IF (rank == psd_diag_rank) THEN
        IF (psd_diag_cellx < 1) THEN
          errcode = c_err_bad_setup
          PRINT*, 'ERROR: PSD diagnostic cell location is incorrect'      
        ELSE IF (psd_diag_cellx > nx) THEN
          errcode = c_err_bad_setup
          PRINT*, 'ERROR: PSD diagnostic cell location is incorrect'      
        ELSE IF (psd_diag_sampling_period < dt) THEN
          psd_diag_sampling_period = dt
          PRINT*,'WARNING: PSD diagnostic sampling period must be at least = dt' 
          PRINT*,' - sampling period is set T = dt'
        END IF
      END IF
      IF (psd_diag_xpos > x_max .OR. psd_diag_xpos < x_min) THEN
        errcode = c_err_bad_setup
        IF (rank == 0) THEN
          PRINT*, 'ERROR: PSD diagnotic x-position in out of bounds'
        END IF
      END IF
    END IF

  END FUNCTION psd_diagnostics_block_check
  
  
END MODULE deck_psd_diagnostics_block
