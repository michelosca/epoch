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
!
! Background neutral gas block
! written by M. Osca Engelbrecht

MODULE deck_see_block
#ifdef SEE
  USE strings_advanced
  USE see

  IMPLICIT NONE
  SAVE
  
  PRIVATE
  PUBLIC :: see_deck_initialise, see_deck_finalise
  PUBLIC :: see_block_start, see_block_end
  PUBLIC :: see_block_handle_element, see_block_check
  
  TYPE(see_type), POINTER :: see_block
  LOGICAL :: model_set, name_set
  INTEGER :: species_id
  REAL(num) :: see_rate_buffer
  
CONTAINS

  SUBROUTINE see_deck_initialise

  
  END SUBROUTINE see_deck_initialise

  
  
  SUBROUTINE see_deck_finalise
  
  END SUBROUTINE see_deck_finalise

  
  
  SUBROUTINE see_block_start

    IF (deck_state /= c_ds_first) THEN
      name_set = .FALSE.
      see_block => NULL()
      model_set = .FALSE.
      see_rate_buffer = -1._num
    END IF
  
  END SUBROUTINE see_block_start
  
  

  SUBROUTINE see_block_end

    INTEGER, DIMENSION(2) :: bc_species

    IF (deck_state /= c_ds_first) THEN
      bc_species = species_list(species_id)%bc_particle
      IF (bc_species(1) /= c_bc_open) THEN
        IF (rank == 0) THEN
          WRITE(*,*) '*** WARNING ***'
          WRITE(*,*) 'SEE module wont work at x-min unless open bc is set'
        END IF
      END IF

      IF (bc_species(2) /= c_bc_open) THEN
        IF (rank == 0) THEN
          WRITE(*,*) '*** WARNING ***'
          WRITE(*,*) 'SEE module wont work at x-max unless open bc is set'
        END IF
      END IF

      IF (.NOT.model_set) THEN
        IF (rank == 0) THEN
          WRITE(*,*) '*** ERROR ***'
          WRITE(*,*) 'SEE model has not been defined'
          CALL abort_code(c_err_required_element_not_set)
        END IF
      END IF

    END IF

  END SUBROUTINE see_block_end
  
  
  
  FUNCTION see_block_handle_element(element, value) RESULT(errcode)
  
    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode, i
    CHARACTER(string_length) :: species_name

    errcode = c_err_none

    IF (deck_state == c_ds_first) RETURN

    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'impact_species')) THEN
      DO i = 1, n_species
        species_name = TRIM(ADJUSTL(species_list(i)%name))
        IF (str_cmp(value, species_name)) THEN
          species_id = i
          CALL init_see_block(i)
          see_block => species_list(i)%see
          name_set = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT.name_set) THEN
        IF (rank == 0) THEN
          WRITE(*,*) '*** ERROR ***'
          WRITE(*,*) 'SEE species "impact_species" is missing/not recognized'
          CALL abort_code(c_err_required_element_not_set)
        END IF
        extended_error_string = 'SEE impact species'
        errcode = c_err_required_element_not_set
      END IF
    END IF

    IF (.NOT. name_set) THEN
      IF (rank == 0) THEN
        WRITE(*,*) '*** ERROR ***'
        WRITE(*,*) 'SEE "impact_species" must be defined first'
        CALL abort_code(c_err_required_element_not_set)
      END IF
      extended_error_string = 'SEE impact species'
      errcode = c_err_required_element_not_set
      RETURN
    END IF
    
    IF (str_cmp(element, 'see_species')) THEN
      DO i = 1, n_species
        species_name = TRIM(ADJUSTL(species_list(i)%name))
        IF (str_cmp(value, species_name)) THEN
          see_block%species_electron = i 
        END IF
        EXIT
      END DO
      RETURN
    END IF
    
    IF (str_cmp(element, 'see_model')) THEN
      ! Input in electron mass units
      CALL setup_see_model(see_block, value, errcode)

      IF (see_rate_buffer > 0._num .AND. see_block%const_yield_model) THEN
        see_block%cross_section(1) = see_rate_buffer
        see_rate_buffer = -1._num
      END IF

      model_set = .TRUE.
      ! Check a see_model has been detected
      IF (errcode /= c_err_none) THEN
        IF (rank == 0) THEN
          WRITE(*,*) '*** ERROR ***'
          WRITE(*,*) 'SEE model definition is not recognized'
          CALL abort_code(c_err_required_element_not_set)
        END IF
        extended_error_string = 'SEE model'
        errcode = c_err_required_element_not_set
        RETURN
      END IF
    END IF
    
    IF (str_cmp(element, 'see_rate')) THEN
      IF (see_block%const_yield_model) THEN
        see_block%cross_section(1) = as_real_print(value, element, errcode)
      ELSE
        see_rate_buffer = as_real_print(value, element, errcode)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'see_temp')) THEN
      see_block%see_temp = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'see_temp_ev')) THEN
      see_block%see_temp = as_real_print(value, element, errcode) * 1.16e4_num
      RETURN
    END IF

  END FUNCTION see_block_handle_element
  
  
  
  FUNCTION see_block_check() RESULT(errcode)
  
    INTEGER :: errcode

    errcode = c_err_none

    IF (.NOT.ASSOCIATED(see_block)) RETURN

    IF (see_block%const_yield_model) THEN
      IF (see_block%cross_section(1) < 0._num) THEN
        IF (rank == 0) THEN
          PRINT*, '*** ERROR ***'
          PRINT*, 'SEE rate has not been defined for the constant yield model'
        END IF
        errcode = c_err_missing_elements
      END IF

      IF (see_block%cross_section(1) > 1._num) THEN
        IF (rank == 0) THEN
          PRINT*, '*** ERROR ***'
          PRINT*, 'SEE rate for the constant yield model must be <= 1'
        END IF
        errcode = c_err_missing_elements
      END IF
    END IF

    IF (see_block%species_electron <= 0) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'Electron species has not been defined for SEE'
      END IF
      errcode = c_err_missing_elements
    END IF
    
    IF (see_block%see_temp <= 0._num) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'Secondary electron tempererature, "see_temp(_ev)", ', &
          'has not been defined'
      END IF
      errcode = c_err_missing_elements
    END IF
  END FUNCTION see_block_check
#endif
END MODULE deck_see_block
