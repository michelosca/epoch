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
  PUBLIC :: inductive_deck_initialise, inductive_deck_finalise
  PUBLIC :: inductive_block_start, inductive_block_end
  PUBLIC :: inductive_block_handle_element, inductive_block_check
  
CONTAINS
  
  SUBROUTINE inductive_deck_initialise
    
  END SUBROUTINE inductive_deck_initialise
  
  
  
  SUBROUTINE inductive_deck_finalise

    INTEGER :: ispecies
    REAL(num) :: x_min, x_max, L_source, electron_weight
    CHARACTER(LEN=string_length) :: e_name, species_name

    IF (deck_state /= c_ds_first .AND. inductive_heating_flag) THEN 
        ! Set electron species parameters
        e_name = inductive_source%e_name
        DO ispecies = 1, n_species
          species_name = TRIM(ADJUSTL(species_list(ispecies)%name)) 
          IF (str_cmp(species_name, e_name)) THEN
            electron_weight = species_list(ispecies)%weight
            inductive_source%id_electrons = ispecies
          END IF
        END DO
        
        ! Set source length
        x_max = inductive_source%x_max
        x_min = inductive_source%x_min

        ! Set conduction current factor
        L_source = x_max - x_min

        inductive_source%j_cond_fac = q0 * electron_weight / L_source 
    END IF
  
  END SUBROUTINE inductive_deck_finalise
  
  
  
  SUBROUTINE inductive_block_start
  
    IF (deck_state == c_ds_first) RETURN
    ALLOCATE(inductive_source)
    CALL init_inductive_source_block(inductive_source)
    inductive_heating_flag = .TRUE.
  
  END SUBROUTINE inductive_block_start
  
  
  
  SUBROUTINE inductive_block_end
  
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
    
    IF (str_cmp(element, 't_profile')) THEN
      CALL initialise_stack(inductive_source%time_function)
      CALL tokenize(value, inductive_source%time_function, errcode)
      ! evaluate it once to check that it's a valid block
      dummy = evaluate(inductive_source%time_function, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'electron_species_name')) THEN
      inductive_source%e_name = TRIM(ADJUSTL(value))
      RETURN
    END IF
    
  END FUNCTION inductive_block_handle_element
  
  
  
  FUNCTION inductive_block_check() RESULT(errcode)
  
    INTEGER :: errcode

    errcode = c_err_none
    IF (inductive_source%x_min >= inductive_source%x_max) THEN
      errcode = c_err_missing_elements
      IF (rank == 0) THEN
        WRITE(*,*) '*** ERROR ***'
        WRITE(*,*) 'Source boundaries must be x_max > x_min'
      END IF
    END IF
    IF (inductive_source%id_electrons <= 0) THEN
      errcode = c_err_missing_elements
      IF (rank == 0) THEN
        WRITE(*,*) '*** ERROR ***'
        WRITE(*,*) 'Electron species was not identified'
      END IF
    END IF


  END FUNCTION inductive_block_check
#endif

END MODULE deck_inductive_block