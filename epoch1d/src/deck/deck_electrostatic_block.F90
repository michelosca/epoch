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

MODULE deck_electrostatic_block
#ifdef ELECTROSTATIC
  USE strings_advanced
  USE electrostatic
  
  IMPLICIT NONE
  SAVE
  
  PRIVATE
  PUBLIC :: electrostatic_deck_initialise, electrostatic_deck_finalise
  PUBLIC :: electrostatic_block_start, electrostatic_block_end
  PUBLIC :: electrostatic_block_handle_element, electrostatic_block_check
  
  TYPE(potential_block), POINTER :: potential_source
  LOGICAL :: boundary_set
  INTEGER :: boundary_es
  
CONTAINS
  
  SUBROUTINE electrostatic_deck_initialise
    
    n_potential_source_x_min = 0
    n_potential_source_x_max = 0
    
  END SUBROUTINE electrostatic_deck_initialise
  
  
  
  SUBROUTINE electrostatic_deck_finalise
  
  END SUBROUTINE electrostatic_deck_finalise
  
  
  
  SUBROUTINE electrostatic_block_start
  
    IF (deck_state == c_ds_first) RETURN
  
    boundary_set = .FALSE.
  
    ! Every new laser uses the internal time function
    ALLOCATE(potential_source)
  
  END SUBROUTINE electrostatic_block_start
  
  
  
  SUBROUTINE electrostatic_block_end
  
    IF (deck_state == c_ds_first) RETURN

    CALL attach_potential(potential_source)
  
  END SUBROUTINE electrostatic_block_end
  
  
  
  FUNCTION electrostatic_block_handle_element(element, value) RESULT(errcode)
  
    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    INTEGER :: io, iu
    REAL(num) :: dummy

    errcode = c_err_none
    
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN
    
    IF (str_cmp(element, 'boundary')) THEN
      ! If the boundary has already been set, simply ignore further calls to it
      IF (boundary_set) RETURN
      boundary_es = as_boundary_print(value, element, errcode)
      boundary_set = .TRUE.
      add_potential_source(boundary_es) = .TRUE.
      CALL init_potential(boundary_es, potential_source)
      RETURN
    END IF

    IF (.NOT. boundary_set) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Cannot set potential source properties', & 
            ' before boundary is set'
        END DO
        CALL abort_code(c_err_required_element_not_set)
      END IF
      extended_error_string = 'boundary'
      errcode = c_err_required_element_not_set
      RETURN
    END IF
    
    IF (str_cmp(element, 'amp')) THEN
      potential_source%amp = as_real_print(value, element, errcode)
      RETURN
    END IF
    
    IF (str_cmp(element, 'profile')) THEN
      CALL initialise_stack(potential_source%profile_function)
      CALL tokenize(value, potential_source%profile_function, errcode)
      CALL potential_update_spatial_profile(potential_source)
      potential_source%use_profile_function = .TRUE.
      ! evaluate it once to check that it's a valid block
      dummy = evaluate(potential_source%profile_function, errcode)
      RETURN
    END IF
    
    IF (str_cmp(element, 't_start')) THEN
      potential_source%t_start = as_time_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 't_end')) THEN
      potential_source%t_end = as_time_print(value, element, errcode)
      RETURN
    END IF
    
    IF (str_cmp(element, 't_profile')) THEN
      CALL initialise_stack(potential_source%time_function)
      CALL tokenize(value, potential_source%time_function, errcode)
      potential_source%use_time_function = .TRUE.
      ! evaluate it once to check that it's a valid block
      dummy = evaluate(potential_source%time_function, errcode)
      RETURN
    END IF
    
    IF (str_cmp(element, 'id')) THEN
      potential_source%id = as_integer_print(value, element, errcode)
      RETURN
    END IF
    
  END FUNCTION electrostatic_block_handle_element
  
  
  
  FUNCTION electrostatic_block_check() RESULT(errcode)
  
    INTEGER :: errcode
    TYPE(potential_block), POINTER :: current
    INTEGER :: error, io, iu

    errcode = c_err_none

    error = 0
    current => potential_list_x_min
    DO WHILE(ASSOCIATED(current))
      IF (ABS(current%amp) < TINY(0._num)) error = IOR(error, 1)
      current => current%next
    END DO

    current => potential_list_x_max
    DO WHILE(ASSOCIATED(current))
      IF (ABS(current%amp) < TINY(0._num)) error = IOR(error, 1)
      current => current%next
    END DO

    IF (IAND(error, 1) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define an "amp" for every potential source.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF
  
  END FUNCTION electrostatic_block_check
#endif  
END MODULE deck_electrostatic_block
