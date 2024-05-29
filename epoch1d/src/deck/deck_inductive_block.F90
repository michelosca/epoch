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

MODULE deck_inductive_block
#ifdef ELECTROSTATIC
  USE strings_advanced
  USE inductive_heating
  
  IMPLICIT NONE
  SAVE
  
  PRIVATE

  TYPE(inductive_heating_block), POINTER :: inductive_source
  INTEGER :: inductive_block_counter

  PUBLIC :: inductive_deck_initialise, inductive_deck_finalise
  PUBLIC :: inductive_block_start, inductive_block_end
  PUBLIC :: inductive_block_handle_element, inductive_block_check
  
CONTAINS
  
  SUBROUTINE inductive_deck_initialise
    
    IF (deck_state /= c_ds_first) THEN 
      CALL inductive_heating_allocate_sources_list(inductive_block_counter)
    END IF
  END SUBROUTINE inductive_deck_initialise
  
  
  
  SUBROUTINE inductive_deck_finalise
    
  END SUBROUTINE inductive_deck_finalise
  
  
  
  SUBROUTINE inductive_block_start
  
    IF (deck_state == c_ds_first) THEN 
      inductive_block_counter = inductive_block_counter + 1
    ELSE
      ALLOCATE(inductive_source)
      CALL inductive_heating_init_inductive_block(inductive_source)
      inductive_heating_flag = .TRUE.
    END IF
  
  END SUBROUTINE inductive_block_start
  
  
  
  SUBROUTINE inductive_block_end

    IF (deck_state == c_ds_first) RETURN
    CALL inductive_heating_add_source_to_list(inductive_source)
    NULLIFY(inductive_source)
  
  END SUBROUTINE inductive_block_end
  
  
  
  FUNCTION inductive_block_handle_element(element, value) RESULT(errcode)
  
    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    REAL(num) :: dummy

    errcode = c_err_none
    
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN
    
    IF (str_cmp(element, 'amp')) THEN
      inductive_source%j0_amp = as_real_print(value, element, errcode)
      RETURN
    END IF
    
    IF (str_cmp(element, 'x_min')) THEN
      inductive_source%x_min = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'x_max')) THEN
      inductive_source%x_max = as_real_print(value, element, errcode)
      RETURN
    END IF
    
    IF (str_cmp(element, 't_start')) THEN
      inductive_source%t_start = as_time_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 't_end')) THEN
      inductive_source%t_end = as_time_print(value, element, errcode)
      RETURN
    END IF
    
    IF (str_cmp(element, 'profile')) THEN
      CALL initialise_stack(inductive_source%function)
      CALL tokenize(value, inductive_source%function, errcode)
      ! evaluate it once to check that it's a valid block
      dummy = evaluate(inductive_source%function, errcode)
      RETURN
    END IF

  END FUNCTION inductive_block_handle_element
  
  
  
  FUNCTION inductive_block_check() RESULT(errcode)
  
    INTEGER :: errcode

    errcode = c_err_none

  END FUNCTION inductive_block_check
#endif

END MODULE deck_inductive_block