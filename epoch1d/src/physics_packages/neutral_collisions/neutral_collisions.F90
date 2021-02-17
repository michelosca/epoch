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

! Non-relativistic neutral-neutral and charged-neutral collision module
! written by M. Osca Engelbrecht
! particle-particle collisions based on G. A. Bird (Oxford Science 
! Publications, 1994) and K. Nanbu (IEEE transactions on plasma science, 
! Vol. 28, No. 3, 2000). Background-particle interactions based on V. Vahedi
! and M. Surendra (Computer Physics Communications, Vol. 87, 179-198, 1995)

MODULE neutral_collisions

  USE shared_data
  USE background
  USE nc_subroutines
  USE nc_auxiliary
  USE random_generator
  USE collisions
!   USE helper

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: particle_neutral_collisions

CONTAINS

  SUBROUTINE particle_neutral_collisions

    INTEGER :: species1, species2, ibg, ix
    REAL(num) :: idx, gsigma_max, w_max
    TYPE(current_collision_block), POINTER :: collision => NULL()
    TYPE(particle_list), POINTER :: p_list1, p_list2
    TYPE(neutrals_block), POINTER :: collision_block

    idx = 1._num/dx

    ! Update background gas spatial & temporal profile
    DO ibg = 1, n_backgrounds
      ! Update density profile
      CALL bg_update_dens_profile(background_list(ibg))
      ! Update temperature profile
      CALL bg_update_temp_profile(background_list(ibg))
    END DO

    DO species1 = 1, n_species
      ! Thic checks whether the species has any collision partner at all
      IF ( .NOT.ANY(neutral_coll(species1,:)) ) CYCLE

      CALL init_current_collision_block(collision)
      ! species1 mass & weight
      collision%m1 = species_list(species1)%mass
      collision%im1 = 1._num/collision%m1
#ifdef PER_SPECIES_WEIGHT
      collision%w1 = species_list(species1)%weight
#endif
      !Randomize species1's particle list
      DO ix = 1, nx
        p_list1 => species_list(species1)%secondary_list(ix)
        CALL shuffle_particle_list_random(p_list1)
      END DO ! ix
      
      
      DO species2 = species1, n_species_bg
        ! Neutral collision?
        IF (.NOT.neutral_coll(species1, species2)) CYCLE
        
        ! Collision factor defined in the input.deck if <=0 then no-collision
        IF (coll_pairs(species1, species2) <= 0._num) CYCLE

        ! Get collision block
        collision_block => species_list(species1)%neutrals(species2)

        ! Set current collision block
        collision%species1 = species1
        collision%species2 = species2
        IF (species1 == species2) THEN
          collision%same_species = .TRUE.
        ELSE
          collision%same_species = .FALSE.
        END IF
        collision%collision_block => collision_block

        ! Background collisions
        IF (collision_block%is_background) THEN
          collision%bg_species = species2 - n_species
          CALL setup_background_collisions(collision)
          CYCLE
        END IF

        !Maximum(g*cross-section)
        gsigma_max = collision_block%gsigma_max_total
        IF (gsigma_max <= 0._num) CYCLE

        ! species2 mass & weight
        collision%m2 = species_list(species2)%mass
        collision%im2 = 1._num/collision%m2
        collision%m12 = collision%m1 + collision%m2
        collision%im12 = 1._num/collision%m12
        collision%reducedm = collision%m1*collision%m2*collision%im12
        collision%ireducedm = 1._num/collision%reducedm
        w_max = collision_block%max_weight
#ifdef PER_SPECIES_WEIGHT
        collision%w2 = species_list(species2)%weight
        collision%w1_ratio = collision%w1 / w_max
        collision%w2_ratio = collision%w2 / w_max
#endif

        !Maximum collision probability
        ! after the following lines pmax is not fully calculated, the number
        !of target particles 'jcount' needs to be multiplied
        collision%pmax = dt*gsigma_max*idx*w_max
        ! Collisions per grid
        DO ix = 1, nx
          collision%ix = ix
          collision%p_list1 => species_list(species1)%secondary_list(ix)

          IF (collision%same_species) THEN
            CALL intra_species_collisions_neutrals(collision)

          ELSE
            p_list2 => species_list(species2)%secondary_list(ix)
            ! Shuffle target particle list, but only if not been done before
            IF (p_list2%coll_counter==0) &
              CALL shuffle_particle_list_random(p_list2)

            collision%p_list2 => p_list2
            CALL inter_species_collisions_neutrals(collision)

            collision%p_list2 => NULL()
          END IF
          collision%p_list1 => NULL()
        END DO ! ix

        collision%collision_block => NULL()
      END DO ! species2
      CALL end_current_collision_block(collision)
    END DO ! species1


  END SUBROUTINE particle_neutral_collisions



  SUBROUTINE inter_species_collisions_neutrals(collision)

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    INTEGER(8) :: n_part1, n_part2, counter1, counter2, ipair, int_coll_pairs
    REAL(num) :: n_coll_pairs
    TYPE(particle), POINTER :: current, impact

    ! Number of particles in the cell
    n_part1 = collision%p_list1%count
    IF (n_part1 == 0) RETURN

    n_part2 = collision%p_list2%count
    IF (n_part2 == 0) RETURN

    ! Collision pairs
    n_coll_pairs = n_part1*n_part2 * collision%pmax ! Real
    int_coll_pairs = CEILING(n_coll_pairs)  ! Integer

    ! Number of collisions undergone by this list in this dt cycle
    counter1 = collision%p_list1%coll_counter
    counter2 = collision%p_list2%coll_counter

    n_part1 = n_part1 - counter1
    n_part2 = n_part2 - counter2
    IF (int_coll_pairs > MIN(n_part1, n_part2)) THEN
      ! Checking whether pairing is performing correctly
      CALL warning_coll_pairs(collision, int_coll_pairs)
      int_coll_pairs = MIN(n_part1, n_part2)
      n_coll_pairs = int_coll_pairs
    END IF
    
    ! The conversion from real to integer of n_coll_pairs induces an error
    !that is prevented with boyd_factor (Boyd 2017, Chapter 6)
    collision%prob_factor = n_coll_pairs/REAL(int_coll_pairs, num)

    ! Particle pairing
    current => collision%p_list1%head
    impact => collision%p_list2%head
    ! Point to the first particle in the list that is going to collide.
    !counter1/2 discards the ones that may have collided already with other
    !species
    DO ipair = 1, counter1
      current => current%next
    END DO
    DO ipair = 1, counter2
      impact => impact%next
    END DO

    ! Execute collisions
    collision%part1 => current
    collision%part2 => impact
    DO ipair = 1, int_coll_pairs

      CALL particle_collision_dynamics(collision)

      collision%part1 => collision%part1%next
      collision%part2 => collision%part2%next

    END DO

    collision%p_list1%coll_counter = counter1 + int_coll_pairs
    collision%p_list2%coll_counter = counter2 + int_coll_pairs

  END SUBROUTINE inter_species_collisions_neutrals



  SUBROUTINE intra_species_collisions_neutrals( collision )

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    INTEGER(8) :: n_part1, counter1, ipair, int_coll_pairs
    REAL(num) :: n_coll_pairs
    TYPE(particle), POINTER :: current, impact

    ! Number of particles in the cell
    n_part1 = collision%p_list1%count
    IF (n_part1 == 0) RETURN

    ! Collision pairs
    n_coll_pairs = 0.5_num * n_part1 * (n_part1-1) * collision%pmax ! Real
    int_coll_pairs = CEILING(n_coll_pairs)  ! Integer

    ! Number of collisions undergone by this list in this dt cycle
    counter1 = collision%p_list1%coll_counter

    n_part1 = n_part1 - counter1
    IF (int_coll_pairs > n_part1) THEN
      ! Checkin whether pairing is performing correctly
      CALL warning_coll_pairs(collision, int_coll_pairs)
      int_coll_pairs = n_part1
      n_coll_pairs = int_coll_pairs
    END IF

    ! The conversion from real to integer of n_coll_pairs induces an error
    !that is prevented with boyd_factor (Boyd 2017, Chapter 6)
    collision%prob_factor = n_coll_pairs/REAL(int_coll_pairs, num)

    ! Particle pairing
    current => collision%p_list1%head
    ! Point to the first particle in the list that is going to collide.
    !counter1 discards the ones that may have collided already with other
    !species
    DO ipair = 1, counter1
      current => current%next
    END DO
    impact => current%next
    
    ! Execute collisions
    collision%part1 => current
    collision%part2 => impact
    DO ipair = 1, int_coll_pairs
      CALL particle_collision_dynamics(collision)
      
      collision%part1 => collision%part2%next
      collision%part2 => collision%part1%next
      
    END DO
    
    collision%p_list1%coll_counter = counter1 + int_coll_pairs*2

  END SUBROUTINE intra_species_collisions_neutrals



  SUBROUTINE particle_collision_dynamics(collision)

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    INTEGER :: ncolltypes, coll_type, the_coll_type, ix
    REAL(num) :: p_total, p_low, p_top, igsigma_max
    REAL(num) :: ran1
    REAL(num) :: g_mag, g2, ke
    REAL(num), DIMENSION(3) :: u_1, u_2, g
    REAL(num), ALLOCATABLE, DIMENSION(:) :: sigma_arr, prob_arr
    TYPE(neutrals_block), POINTER :: collision_block
    TYPE(collision_type_block), POINTER :: coll_type_block

    collision_block => collision%collision_block

    igsigma_max = collision_block%igsigma_max_total
    ncolltypes = collision_block%ncolltypes
    ALLOCATE(sigma_arr(ncolltypes), prob_arr(ncolltypes))

    !Particles momenta
    u_1 = collision%part1%part_p * collision%im1
    IF (collision_block%is_background) THEN
      u_2 = random_background_velocity(collision)
      collision%u_2 = u_2
#ifndef PER_SPECIES_WEIGHT
      collision%w1 = collision%part1%weight
#endif
    ELSE
      u_2 = collision%part2%part_p * collision%im2
#ifndef PER_SPECIES_WEIGHT
      collision%w1 = collision%part1%weight
      collision%w2 = collision%part2%weight
      collision%w1_ratio = collision%w1 / collision_block%max_weight
      collision%w2_ratio = collision%w2 / collision_block%max_weight
#endif
    END IF
    !Relative velocity
    g = u_1 - u_2                 ! vector
    collision%g = g
    g2 = DOT_PRODUCT(g, g)        ! scalar^2
    g_mag = SQRT(g2)              ! scalar

    ! Relative kinetic energy
    ke = 0.5_num * collision%reducedm * g2

    ! Specific cross-section values for the calculated energies
    CALL get_coll_cross_section(collision, sigma_arr, ke)
    
    ! Collision probabilities
    prob_arr = sigma_arr * g_mag * igsigma_max  ! Probability array
    prob_arr = prob_arr * collision%prob_factor! User and Boyd(2017) factors
    p_total = SUM(prob_arr)

    ! Checkin whether MAX(g*sigma) has been chosen correctly
    IF (p_total>1._num) THEN
      CALL warning_pmax(collision, SUM(sigma_arr), ke)
      ! What should be done in this case? Maybe have a second collision?
    END IF

    ! Generate a random number
    ran1 = random()

    ! First check if there is a collision
    IF (ran1 <= p_total) THEN !if TRUE then -> collision

      ! Determines which type of collision (including null collision)
      the_coll_type = ncolltypes
      p_low = 0._num
      p_top = prob_arr(1)
      DO coll_type = 1, ncolltypes - 1
        IF ( ran1 <= p_top .AND. ran1 >= p_low ) THEN
          the_coll_type = coll_type
          EXIT
        END IF
        p_low = p_top
        p_top = p_top + prob_arr(coll_type+1)
      END DO
      coll_type_block => collision_block%collision_set(the_coll_type)
      collision%type_block => coll_type_block
      collision%type = the_coll_type

      ! IO collision data
      IF (neutral_collision_counter) THEN
        ix = collision%ix
        coll_type_block%coll_counter(ix) = coll_type_block%coll_counter(ix)+1
      END IF

      ! Store velocity values (used later in coll_subroutine)
      collision%g_mag = g_mag
      collision%u_cm = (collision%part1%part_p+collision%m2*u_2)*collision%im12

      ! Calculate post-collision momenta
      CALL coll_type_block%coll_subroutine(collision)

      collision%type_block => NULL()
    END IF

    DEALLOCATE(sigma_arr, prob_arr)

  END SUBROUTINE particle_collision_dynamics



  SUBROUTINE get_coll_cross_section(collision, sigma_arr, ke_cm)
    ! This routine outputs the cross-sections for the different collision types

    TYPE(current_collision_block), POINTER, INTENT(IN) :: collision

    REAL(num), DIMENSION(:), INTENT(OUT) :: sigma_arr
    REAL(num), INTENT(IN) :: ke_cm
    
    REAL(num) :: sigma, ke_min, ke_max, sigma_min, sigma_max, ethreshold
    REAL(num) :: user_factor
    INTEGER :: table_len, coll_type, ix
    TYPE(collision_type_block), POINTER :: coll_type_block
    TYPE(neutrals_block), POINTER :: collision_block

    sigma_arr = 0._num
    collision_block => collision%collision_block
    
    DO coll_type = 1, collision_block%ncolltypes
    
      coll_type_block => collision_block%collision_set(coll_type)
      user_factor = coll_type_block%user_factor

      ! Collision threshold energy
      ethreshold = coll_type_block%ethreshold
      
      !Check whether collision types with e_threshold!=0 are possible
      IF (ke_cm < ethreshold) CYCLE
      
      ! Cross-section
      !Negative cross-section is initially assigned
      sigma = -1._num
      !Length of the cross-section table
      table_len = coll_type_block%table_len
      IF (table_len == 1) THEN
        ! Only one line: cross-section is constant
        sigma = coll_type_block%cross_section(1)
      ELSE
        ! cross-section is energy dependent
        DO ix = 2, table_len
          ke_min = coll_type_block%energy(ix-1)
          ke_max = coll_type_block%energy(ix)
          IF (ke_cm < ke_max .AND. ke_cm >= ke_min ) THEN
            sigma_min = coll_type_block%cross_section(ix-1)
            sigma_max = coll_type_block%cross_section(ix)
            ! Interpolation
            sigma = interpolation(ke_min, ke_max, sigma_min, sigma_max, ke_cm)
            sigma = sigma * user_factor
            EXIT
          END IF
        END DO
        CALL crosssection_out_of_range(sigma, ke_cm, coll_type_block)
      END IF

      !Assig calculated value to the output array
      sigma_arr(coll_type) = sigma
    END DO

  END SUBROUTINE get_coll_cross_section



  SUBROUTINE setup_background_collisions(collision)

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num) :: gsigma_max, density, temp
    REAL(num) :: ix_density
    INTEGER :: ix, species1
    TYPE(background_block), POINTER :: background

    species1 = collision%species1
    background => collision%collision_block%background
    IF (time < background%t_start .OR. time > background%t_end) RETURN 

    !Maximum(g*cross-section)
    gsigma_max = collision%collision_block%gsigma_max_total
    IF (gsigma_max <= 0._num) RETURN
        
    ! Background gas mass
    collision%m2 = background%mass
    collision%im2 = 1._num/collision%m2
    collision%m12 = collision%m1 + collision%m2
    collision%im12 = 1._num/collision%m12
    collision%reducedm = collision%m1*collision%m2*collision%im12
    collision%ireducedm = 1._num/collision%reducedm

    ! Background gass density & temperature 'amplitude'
    density = background%dens
    temp = background%temp
          
    DO ix = 1, nx
      collision%ix = ix
      collision%p_list1 => species_list(species1)%secondary_list(ix)

      ! Background density
      ix_density = background%dens_profile(ix) * density
      ! Background temperature
      collision%ix_temp = background%temp_profile(ix) * temp
      collision%pmax = ix_density*gsigma_max*dt

      CALL background_collisions_neutrals(collision)

      collision%p_list1 => NULL()
    END DO ! ix

    collision%collision_block => NULL()

  END SUBROUTINE setup_background_collisions



  SUBROUTINE background_collisions_neutrals(collision)

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    INTEGER(8) :: n_part1, counter1, int_coll_pairs, ipair
    REAL(num) :: n_coll_pairs

    TYPE(particle), POINTER :: current
    
    ! Number of particles in the cell
    n_part1 = collision%p_list1%count
    IF (n_part1 == 0) RETURN
   
    ! Collision pairs
    n_coll_pairs = n_part1 * collision%pmax ! Real
    int_coll_pairs = CEILING(n_coll_pairs)  ! Integer

    ! Number of collisions undergone by this list in this dt cycle
    counter1 = collision%p_list1%coll_counter

    n_part1 = n_part1 - counter1
    IF (int_coll_pairs > n_part1) THEN
      ! Checking whether pairing is performing correctly
      CALL warning_coll_pairs(collision, int_coll_pairs)
      int_coll_pairs = n_part1
      n_coll_pairs = int_coll_pairs
    END IF
    
    ! The conversion from real to integer of n_coll_pairs induces an error
    !that is prevented with boyd_factor (Boyd 2017, Chapter 6)
    collision%prob_factor = n_coll_pairs/REAL(int_coll_pairs, num)
    
    ! Particle pairing
    current => collision%p_list1%head
    ! Point to the first particle in the list that is going to collide
    !counter1 discards the ones that may have collided already with other
    !species
    DO ipair = 1, counter1
      current => current%next
    END DO

    ! Execute collisions
    collision%part1 => current
    DO ipair = 1, int_coll_pairs
      CALL particle_collision_dynamics(collision)
      collision%part1 => collision%part1%next
    END DO
    
    collision%p_list1%coll_counter = counter1 + int_coll_pairs
    collision%part1 => NULL()

  END SUBROUTINE background_collisions_neutrals



  SUBROUTINE init_current_collision_block(collision)

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    ALLOCATE(collision)

    ! Cell position
    collision%ix = -1

    ! Species
    collision%species1 = -1
    collision%species2 = -1
    collision%same_species = .FALSE.
    collision%bg_species = -1
    collision%p_list1 => NULL()
    collision%p_list2 => NULL()
    collision%part1 => NULL()
    collision%part2 => NULL()

    ! Collision
    collision%collision_block => NULL()
    collision%type_block => NULL()
    collision%pmax = 0._num
    collision%prob_factor = 1._num
    collision%type = -1

    ! Collision dynamics
    collision%g_mag = 0._num
    collision%u_cm = 0._num
    collision%g = 0._num
    collision%u_2 = 0._num

    ! Particle's mass
    collision%m1 = 0._num
    collision%m2 = 0._num
    collision%im1 = 0._num
    collision%im2 = 0._num
    collision%m12 = 0._num
    collision%im12 = 0._num
    collision%reducedm = 0._num
    collision%ireducedm = 0._num

    ! Particle's weight
    collision%w1 = 0._num
    collision%w2 = 0._num
    collision%w1_ratio = 0._num
    collision%w2_ratio = 0._num
    
    ! Background temperature
    collision%ix_temp = 0._num

  END SUBROUTINE init_current_collision_block



  SUBROUTINE end_current_collision_block(collision)

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    collision%p_list1 => NULL()
    collision%p_list2 => NULL()
    collision%part1 => NULL()
    collision%part2 => NULL()

    ! Collision
    collision%collision_block => NULL()
    collision%type_block => NULL()
    
    DEALLOCATE(collision)

  END SUBROUTINE end_current_collision_block

END MODULE neutral_collisions