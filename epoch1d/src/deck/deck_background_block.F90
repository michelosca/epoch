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

MODULE deck_background_block

  USE strings_advanced
  USE background
  USE constants

  IMPLICIT NONE
  SAVE
  
  PRIVATE
  PUBLIC :: background_deck_initialise, background_deck_finalise
  PUBLIC :: background_block_start, background_block_end
  PUBLIC :: background_block_handle_element, background_block_check
  
  LOGICAL :: name_set
  INTEGER :: bg_counter
  TYPE(background_block), POINTER :: working_background
  
CONTAINS

  SUBROUTINE background_deck_initialise
  
    IF (deck_state == c_ds_first) THEN
      n_backgrounds = 0
    ELSE
      ALLOCATE(background_list(n_backgrounds))
      bg_counter = 0
    END IF
  
  END SUBROUTINE background_deck_initialise

  
  
  SUBROUTINE background_deck_finalise
  
    INTEGER :: jspecies, ibg, bg_species
    
    IF (deck_state == c_ds_first) THEN
      n_species_bg = n_species + n_backgrounds
      RETURN
    END IF
    
    ! Add background gas' collision factors to collision matrices
    DO ibg = 1, n_backgrounds
      bg_species = n_species + ibg
      DO jspecies = 1, n_species
        IF (background_list(ibg)%colliding_species(jspecies)) THEN
          coll_pairs(bg_species, jspecies) = 1._num
          coll_pairs(jspecies, bg_species) = 1._num
          neutral_coll(bg_species, jspecies) = .TRUE.
          neutral_coll(jspecies, bg_species) = .TRUE.
        END IF
      END DO
    END DO
  
  END SUBROUTINE background_deck_finalise

  
  
  SUBROUTINE background_block_start
  
    name_set = .FALSE.
    IF (deck_state == c_ds_first) RETURN

    IF (n_backgrounds > 0) THEN
      bg_counter = bg_counter + 1
      working_background => background_list(bg_counter)
    END IF

  END SUBROUTINE background_block_start
  
  

  SUBROUTINE background_block_end
  
    IF (deck_state == c_ds_first) RETURN
    NULLIFY(working_background)
    
  END SUBROUTINE background_block_end
  
  
  
  FUNCTION background_block_handle_element(element, value) RESULT(errcode)
  
    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    REAL(num) :: dummy, electron_masses, temperature_eV
    INTEGER :: io, iu, ispecies
    LOGICAL :: check_species
    
    errcode = c_err_none
    
    IF (deck_state == c_ds_first) THEN
      IF (str_cmp(element, 'name')) THEN
        IF (name_set) RETURN
        name_set = .TRUE.
        n_backgrounds = n_backgrounds + 1
      END IF
      RETURN
    END IF
    
    IF (element == blank .OR. value == blank) RETURN
    
    IF (str_cmp(element, 'name')) THEN
      IF (name_set) RETURN
      name_set = .TRUE.
      CALL init_background(value, working_background)
      RETURN
    END IF
    
    IF (.NOT. name_set) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Cannot set background properties before name is set'
        END DO
        CALL abort_code(c_err_required_element_not_set)
      END IF
      extended_error_string = 'background name'
      errcode = c_err_required_element_not_set
      RETURN
    END IF
    
    IF (str_cmp(element, 'id')) THEN
      working_background%id = as_integer_print(value, element, errcode)
      RETURN
    END IF
    
    IF (str_cmp(element, 'mass')) THEN
      ! Input in electron mass units
      electron_masses = as_real_print(value, element, errcode)
      working_background%mass = electron_masses * m0
      RETURN
    END IF
    
    IF (str_cmp(element, 'density')) THEN
      ! Input in number density [m^-3]
      working_background%dens = as_real_print(value, element, errcode)
      RETURN
    END IF
    
    IF (str_cmp(element, 'mean_velocity_x')) THEN
      ! Input in number density [m^-3]
      working_background%mean_velocity(1) = as_real_print(value, element, &
        errcode)
      RETURN
    END IF
    
    IF (str_cmp(element, 'mean_velocity_y')) THEN
      ! Input in number density [m^-3]
      working_background%mean_velocity(2) = as_real_print(value, element, &
        errcode)
      RETURN
    END IF
    
    IF (str_cmp(element, 'mean_velocity_z')) THEN
      ! Input in number density [m^-3]
      working_background%mean_velocity(3) = as_real_print(value, element, &
        errcode)
      RETURN
    END IF
    
    IF (str_cmp(element, 'temperature_eV')) THEN
      ! Input in eV
      temperature_eV = as_real_print(value, element, errcode)
      working_background%temp = temperature_eV * 1.16e4_num
      RETURN
    END IF
    
    IF (str_cmp(element, 'temperature')) THEN
      ! Input in K
      working_background%temp = as_real_print(value, element, errcode)
      RETURN
    END IF
    
    IF (str_cmp(element, 'colliding_species')) THEN
      IF (str_cmp(value, 'all')) THEN
        working_background%colliding_species = .TRUE.
        RETURN
      END IF
      
      ! Colliding partner's name
      check_species = .FALSE.
      DO ispecies = 1, n_species
        IF (str_cmp(value, species_list(ispecies)%name)) THEN
          working_background%colliding_species(ispecies) = .TRUE.
          check_species = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT.(check_species)) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) 'Background collision partner ', TRIM(ADJUSTL(value)), &
           ' does not match any defined species'
          WRITE(io,*)
        END DO
      END IF
      RETURN
    END IF
    
    IF (str_cmp(element, 'density_profile')) THEN
      working_background%use_dens_profile_fn = .TRUE.
      CALL initialise_stack(working_background%dens_profile_fn)
      CALL tokenize(value, working_background%dens_profile_fn, errcode)
      ! evaluate it once to check that it's a valid block
      dummy = evaluate(working_background%dens_profile_fn, errcode)
      RETURN
    END IF
    
    IF (str_cmp(element, 'temperature_profile')) THEN
      working_background%use_temp_profile_fn = .TRUE.
      CALL initialise_stack(working_background%temp_profile_fn)
      CALL tokenize(value, working_background%temp_profile_fn, errcode)
      ! evaluate it once to check that it's a valid block
      dummy = evaluate(working_background%temp_profile_fn, errcode)
      RETURN
    END IF
    
    IF (str_cmp(element, 't_start')) THEN
      working_background%t_start = as_time_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 't_end')) THEN
      working_background%t_end = as_time_print(value, element, errcode)
      RETURN
    END IF
  
    errcode = c_err_unknown_element
  
  END FUNCTION background_block_handle_element
  
  
  
  FUNCTION background_block_check() RESULT(errcode)
  
    INTEGER :: errcode, ibg
    INTEGER :: error, io, iu

    errcode = c_err_none

    error = 0
    DO ibg = 1, n_backgrounds
      IF (background_list(ibg)%mass <= 0.0_num) error = 1
      IF (background_list(ibg)%dens <= 0.0_num) error = 2
      IF (background_list(ibg)%temp <= 0.0_num) error = 3
    END DO

    IF (error==1) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define a "mass" for every background.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (error==2) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define a "density" for every background.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF
    
    IF (error==3) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define a "temperature">0 for every background.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF
    
  END FUNCTION background_block_check

END MODULE deck_background_block
