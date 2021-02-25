! Copyright (C) 2010      Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009      Chris Brady <C.S.Brady@warwick.ac.uk>
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

MODULE custom_electrostatic

  USE shared_data
  USE constants

  IMPLICIT NONE

CONTAINS

  FUNCTION custom_electric_potential_time_profile(potential_source)

    TYPE(potential_block), INTENT(IN) :: potential_source
    REAL(num) :: custom_electric_potential_time_profile

    custom_electric_potential_time_profile = 1.0_num

  END FUNCTION custom_electric_potential_time_profile

END MODULE custom_electrostatic
