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

MODULE background

  USE evaluator
  USE nc_auxiliary
  USE random_generator

  INTEGER, PRIVATE :: id_counter = 1
  
CONTAINS

  SUBROUTINE init_background(name, background)

    CHARACTER(LEN=string_length), INTENT(IN) :: name
    TYPE(background_block), INTENT(INOUT) :: background

    background%name = name
    background%id = id_counter
    background%use_dens_profile_fn = .FALSE.
    background%use_temp_profile_fn = .FALSE.
    background%mass = -1.0_num
    background%dens = -1.0_num
    background%temp = -1.0_num
    background%mean_velocity = 0.0_num
    
    background%t_start = 0.0_num
    background%t_end = t_end
    
    ALLOCATE(background%colliding_species(n_species))
    background%colliding_species = .FALSE.
    
    ALLOCATE(background%dens_profile(nx))
    background%dens_profile = 1._num
    
    ALLOCATE(background%temp_profile(nx))
    background%temp_profile = 1._num

    id_counter = id_counter + 1

  END SUBROUTINE init_background



  SUBROUTINE deallocate_backgrounds

    INTEGER :: ibg
    TYPE(background_block), POINTER :: background
    
    DO ibg = 1, n_backgrounds
      background => background_list(ibg)
      
      IF (background%use_dens_profile_fn) &
        CALL deallocate_stack(background%dens_profile_fn)

      IF (background%use_temp_profile_fn) &
        CALL deallocate_stack(background%temp_profile_fn)

      DEALLOCATE(background%colliding_species)
      DEALLOCATE(background%dens_profile)
      DEALLOCATE(background%temp_profile)
    END DO
    
    IF (n_backgrounds>0) DEALLOCATE(background_list)

  END SUBROUTINE deallocate_backgrounds
  
  
  
  SUBROUTINE bg_update_dens_profile(background)

    TYPE(background_block), INTENT(INOUT):: background
    INTEGER :: ix, err
    TYPE(parameter_pack) :: parameters

    err = 0
    IF (background%use_dens_profile_fn) THEN
      DO ix = 1,nx
        parameters%pack_ix = ix
        background%dens_profile(ix) = &
          evaluate_with_parameters(background%dens_profile_fn, parameters, err)
      END DO
      RETURN
    END IF
    
!     DO ix = 1,nx
!       background%dens_profile(ix) = custom_dens_profile(background, ix)
!     END DO
   
   END SUBROUTINE bg_update_dens_profile
  
  
  
  SUBROUTINE bg_update_temp_profile(background)

    TYPE(background_block), INTENT(INOUT) :: background
    INTEGER :: ix, err
    TYPE(parameter_pack) :: parameters

    err = 0
    IF (background%use_temp_profile_fn) THEN
      DO ix = 1,nx
        parameters%pack_ix = ix
        background%temp_profile(ix) = &
          evaluate_with_parameters(background%temp_profile_fn, parameters, err)
      END DO
      RETURN
    END IF
    
!     DO ix = 1,nx
!       background%temp_profile(ix) = custom_temp_profile(background, ix)
!     END DO
   
   END SUBROUTINE bg_update_temp_profile

   
   
  FUNCTION random_background_velocity(collision)
 
    REAL(num), DIMENSION(3) :: random_background_velocity
    TYPE(current_collision_block) :: collision
    REAL(num) :: temp, mass
    REAL(num) :: ran1, ran2, variance, speed
    TYPE(background_block), POINTER :: background

    background => collision%collision_block%background
    ! this is to ensure that 0 < ran1 < 1
    ! ran1=0 gives NaN in logarithm
    ! ran1=1 could give positive logarithm due to rounding errors
    ran1 = (1.0_num - 1.0e-10_num) * random() + 0.5e-10_num
    ran2 = 2.0_num * pi * random()
    ! Box Muller method for random from Gaussian distribution,
    ! mean is mean_velocity, variance from background gas temperature
    ! Possible place for speed up by caching the second Box Muller number
    ! and using it later
    ! variance * SQRT(-2.0_num * LOG(ran1)) * COS(ran2) + mean_velocity

    mass = background%mass
    temp = collision%ix_temp
    variance = SQRT(2._num*kb*temp/mass)
    speed = variance * SQRT(-2.0_num * LOG(ran1)) * SIN(ran2)

    random_background_velocity = random_unit_vector(1._num)
    random_background_velocity = random_background_velocity * speed
    random_background_velocity = random_background_velocity + &
      background%mean_velocity

  END FUNCTION random_background_velocity  
  
END MODULE background
