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
! written by M. Osca Engelbrecht

MODULE see

  USE strings_advanced
  USE shared_data

CONTAINS

  SUBROUTINE init_see_block(species_id)

    INTEGER, INTENT(IN) :: species_id
    TYPE(see_type), POINTER :: see_block

    IF (ASSOCIATED(species_list(species_id)%see)) THEN
      IF (rank == 0) THEN
        WRITE(*,*) '*** WARNING ***'
        WRITE(*,*) 'Species ',TRIM(ADJUSTL(species_list(species_id)%name)), &
          ' has already a SEE model associated. First defined models prevails.'
      END IF
      RETURN
    END IF

    ALLOCATE(species_list(species_id)%see)
    see_block => species_list(species_id)%see
    see_block%species_impact = species_id
    see_block%species_electron = -1
    see_block%const_yield_model = .FALSE.
    see_block%see_subroutine => NULL()

  END SUBROUTINE init_see_block

  
  SUBROUTINE setup_see_model(see_block, model_str, errcode)

    TYPE(see_type), POINTER, INTENT(INOUT) :: see_block
    CHARACTER(*), INTENT(IN) :: model_str
    INTEGER, INTENT(INOUT) :: errcode

    IF (str_cmp(model_str, 'constant_yield')) THEN
      see_block%const_yield_model = .TRUE.
      ALLOCATE(see_block%cross_section(1))
      see_block%cross_section = -1._num
      see_block%see_subroutine => see_constant_yield 
    ELSE
      errcode = c_err_required_element_not_set 
    END IF

  END SUBROUTINE setup_see_model


  SUBROUTINE see_constant_yield(current_part, species)

    TYPE(particle), POINTER, INTENT(INOUT) :: current_part
    TYPE(particle_species), POINTER, INTENT(IN) :: species

  END SUBROUTINE

END MODULE see 
