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

MODULE nc_subroutines

  USE shared_data
  USE partlist
  USE nc_auxiliary
  USE random_generator

  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: bird_elastic_scattering, bird_ionisation, bird_excitation
  PUBLIC :: bird_charge_exchange
  PUBLIC :: boyd_elastic_scattering, boyd_excitation
  PUBLIC :: bird_elastic_scattering_bg, bird_ionisation_bg, bird_excitation_bg
  PUBLIC :: bird_charge_exchange_bg
  PUBLIC :: vahedi_electron_elastic_scattering_bg, vahedi_excitation_bg
  PUBLIC :: vahedi_ion_elastic_scattering_bg, vahedi_ionisation_bg
#ifndef PER_SPECIES_WEIGHT
  PUBLIC :: split_ionisation, split_charge_exchange
#endif
CONTAINS

  SUBROUTINE set_particle_properties(collision, species1, species2, &
    m1, im1, m2, im2, w1, w2, g_mag, u_cm, part1, part2, p_list)

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision
    REAL(num), INTENT(OUT) :: m1, im1, m2, im2, w1, w2, g_mag
    REAL(num), DIMENSION(3), INTENT(OUT) :: u_cm
    INTEGER, INTENT(OUT) :: species1, species2
    TYPE(particle), POINTER, INTENT(OUT) :: part1, part2
    TYPE(particle_list), POINTER, INTENT(OUT), OPTIONAL :: p_list

    REAL(num) :: w_ratio_temp

    IF (collision%species1 == collision%type_block%source_species_id) THEN
      species1 = collision%species2
      species2 = collision%species1
      m1 = collision%m2
      im1 = collision%im2
      m2 = collision%m1
      im2 = collision%im1
      w1 = collision%w2
      w2 = collision%w1
      w_ratio_temp = collision%w1_ratio
      collision%w1_ratio = collision%w2_ratio
      collision%w2_ratio = w_ratio_temp
      part1 => collision%part2
      part2 => collision%part1
      IF (PRESENT(p_list)) p_list => collision%p_list1
    ELSE
      species1 = collision%species1
      species2 = collision%species2
      m1 = collision%m1
      im1 = collision%im1
      m2 = collision%m2
      im2 = collision%im2
      w1 = collision%w1
      w2 = collision%w2
      part1 => collision%part1
      part2 => collision%part2
      IF (PRESENT(p_list)) p_list => collision%p_list2
    END IF
    g_mag = collision%g_mag
    u_cm = collision%u_cm

  END SUBROUTINE set_particle_properties


#ifndef PER_SPECIES_WEIGHT
  SUBROUTINE check_weight_difference(w2, w1, part2, p_list)

    REAL(num), INTENT(INOUT) :: w2
    REAL(num), INTENT(IN) :: w1
    TYPE(particle), POINTER, INTENT(INOUT) :: part2
    TYPE(particle_list), POINTER , INTENT(INOUT) :: p_list

    REAL(num) :: weight_difference, new_w2
    REAL(num), DIMENSION(3) :: new_p
    TYPE(particle), POINTER :: current

    weight_difference = w2 - w1 ! Neutrals weight should always be larger
    ! Check that w2 - w1 > w1
    IF (weight_difference < 0._num) THEN
      WRITE(*,*) 'weight_difference is below zero: ', weight_difference
    END IF
    IF (weight_difference <= w1) THEN
      current => part2 ! Low weight particle
      part2 => part2%next ! Next particle
      ! Merge current to next particle
      new_w2 = part2%weight
      ! Average momentum
      new_p = part2%part_p*new_w2 + current%part_p*weight_difference
      w2 = new_w2 + weight_difference
      new_p = new_p/w2
      part2%part_p = new_p
      part2%weight = w2
      CALL remove_particle_from_partlist(p_list, current)
      CALL destroy_particle(current)
      NULLIFY(current)
    END IF

  END SUBROUTINE check_weight_difference
#endif


  SUBROUTINE link_particle_pointers(collision, part1, part2)

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision
    TYPE(particle), POINTER, INTENT(IN) :: part1, part2

    REAL(num) :: w_ratio_temp

    IF (collision%species1 == collision%type_block%source_species_id) THEN
      collision%part1 => part2
      collision%part2 => part1
      w_ratio_temp = collision%w1_ratio
      collision%w1_ratio = collision%w2_ratio
      collision%w2_ratio = w_ratio_temp
    ELSE
      collision%part1 => part1
      collision%part2 => part2
    END IF

  END SUBROUTINE link_particle_pointers



  SUBROUTINE release_energy_loss(collision, u_cm)

    ! This subroutine only applies to Boyd's method

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision
    REAL(num), DIMENSION(3), INTENT(INOUT) :: u_cm

    REAL(num), DIMENSION(3) :: u_cm_unit
    REAL(num) :: delta_E
    REAL(num) :: u_cm_mod, u_cm_mod2,  u_cm_mod_corrected
    INTEGER :: ix

    ix = collision%ix
    delta_E = energy_loss(ix)
    IF (delta_E > 0._num .AND. collision%type_block%energy_correction) THEN
      ! u_cm vector is corrected -> energy conservation
      u_cm_mod2 = DOT_PRODUCT(u_cm, u_cm)
      u_cm_mod = SQRT(u_cm_mod2) ! u_cm magnitude
      u_cm_unit = u_cm/u_cm_mod ! u_cm unitary vector
      u_cm_mod_corrected = SQRT(u_cm_mod2 + 2._num*delta_E*collision%im12)
      u_cm = u_cm_unit*u_cm_mod_corrected ! new vector
      ! Update energy_loss array
      energy_loss(ix) = 0._num
    END IF

  END SUBROUTINE release_energy_loss



  SUBROUTINE save_energy_loss(ix, p, p_precoll, w_ratio, im)

    ! This subroutine only applies to Boyd's method

    INTEGER, INTENT(IN) :: ix
    REAL(num), DIMENSION(3), INTENT(IN) :: p_precoll, p
    REAL(num), INTENT(IN) :: w_ratio, im

    REAL(num), DIMENSION(3) :: p_diff
    REAL(num) :: delta_E, p_diff2
    ! Energy loss
    p_diff = p_precoll - p
    p_diff2 = DOT_PRODUCT(p_diff, p_diff)
    delta_E = 0.5_num*(1._num - w_ratio)*w_ratio*p_diff2*im
    energy_loss(ix) = energy_loss(ix) + delta_E

  END SUBROUTINE save_energy_loss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                 BIRD METHOD                                 !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE bird_elastic_scattering(collision)

    ! Collision process: electron scattering
    ! e + N -> e + N

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num) :: reducedm, g_mag, m, ran1
    REAL(num), DIMENSION(3) :: random_direction, p_scat, u_cm

    reducedm = collision%reducedm
    g_mag = collision%g_mag
    u_cm = collision%u_cm

    ! Post-collision momentum
    random_direction = random_unit_vector()
    p_scat = reducedm * g_mag * random_direction

    ran1 = random()
    IF (ran1 <= collision%w2_ratio) THEN
      m = collision%m1
      collision%part1%part_p = u_cm*m + p_scat ! Post-collision momentum
    END IF

    IF (ran1 <= collision%w1_ratio) THEN
      m = collision%m2
      collision%part2%part_p = u_cm*m - p_scat ! Post-collision momentum
    END IF

  END SUBROUTINE bird_elastic_scattering



  SUBROUTINE bird_ionisation(collision)

    ! Collision process: electron impact ionisation
    ! e + N -> 2e + I
    ! species1 = electron (e)
    ! species2 = neutral (N)

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num), DIMENSION(3) :: u_cm, u_e1, u_e2, random_direction
    REAL(num) :: m1, m2, im1, im2, w1, w2
    REAL(num) :: ran_e, ran_w
    REAL(num) :: g_mag, e_threshold, e_excess, u_excess, u_excess_total
    INTEGER :: ion_id, species1, species2
    TYPE(particle), POINTER :: new_part, part1, part2
    TYPE(particle_list) , POINTER :: p_list

    CALL set_particle_properties(collision, species1, species2, &
      m1, im1, m2, im2, w1, w2, g_mag, u_cm, part1, part2, p_list)
    e_threshold = collision%type_block%ethreshold


    ! Excess energy distribution between the two electrons
    e_excess = 0.5_num*collision%reducedm*g_mag*g_mag - e_threshold

    ! Speed associated to the energy excess
    u_excess_total = SQRT(2._num*im1*e_excess)

    ran_e = 0.5_num !random()

    ! Random weight ratio
    ran_w = random()

    ! Current particle: the existing electron
    u_excess = u_excess_total*SQRT(ran_e)
    random_direction = random_unit_vector() ! Random vector
    u_e1 = u_excess*random_direction ! cm frame of reference
    IF (ran_w <= collision%w2_ratio) THEN
      part1%part_p = (u_cm + u_e1)*m1
    END IF

    ! Impact particle: neutral -> ion + electron
    IF (ran_w <= collision%w1_ratio) THEN

      ! The new electron
      CALL create_particle(new_part)
      u_excess = u_excess_total*SQRT(1._num - ran_e)
      random_direction = random_unit_vector() ! Random unitary vector
      u_e2 = u_excess*random_direction  ! cm frame of reference
      new_part%part_p = (u_cm + u_e2)*m1
      new_part%part_pos = part1%part_pos
#ifndef PER_SPECIES_WEIGHT
      new_part%weight = w2 ! Electron's weight
#endif
      CALL add_particle_to_partlist( &
        species_list(species1)%attached_list, new_part)
      NULLIFY(new_part)

      ! New ion: move neutral particle to ion species list
      new_part => part2
      new_part%part_p = u_cm * (m2 - m1) - m1*(u_e1 + u_e2)
      IF (.NOT.ASSOCIATED(p_list%head, TARGET=part2)) THEN
        ! If particle is not the head of the list
        part2 => part2%prev
      ELSE
        ! If particle is the head of the list then one particle is skipped
        p_list%coll_counter = p_list%coll_counter + 1
        part2 => part2%next
      END IF
      CALL remove_particle_from_partlist(p_list, new_part)
      ion_id = collision%type_block%new_species_id ! ions
      CALL add_particle_to_partlist( species_list(ion_id)%attached_list, &
        new_part)
      NULLIFY(new_part)
      
    END IF

    CALL link_particle_pointers(collision, part1, part2)
    NULLIFY(part1, part2)

  END SUBROUTINE bird_ionisation



  SUBROUTINE bird_excitation(collision)

    ! Collision process: excitation
    ! e + N -> e + N*

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    INTEGER :: excited_id, species1, species2
    REAL(num) :: g_mag, e_threshold
    REAL(num) :: m1, im1, m2, im2, w1, w2, reducedm, ireducedm
    REAL(num) :: ran_w
    REAL(num), DIMENSION(3) :: random_direction, p_2, u_cm, p_scat
    TYPE(particle), POINTER :: part1, part2, new_part
    TYPE(particle_list), POINTER :: p_list

    CALL set_particle_properties(collision, species1, species2, &
      m1, im1, m2, im2, w1, w2, g_mag, u_cm, part1, part2, p_list)
    excited_id = collision%type_block%new_species_id
    reducedm = collision%reducedm
    ireducedm = collision%ireducedm
    e_threshold = collision%type_block%ethreshold

    !Random vector
    random_direction = random_unit_vector()
    g_mag = SQRT(g_mag*g_mag - 2._num*e_threshold*ireducedm)
    p_scat = reducedm * g_mag * random_direction

    ran_w = random()
    IF (ran_w <= collision%w2_ratio) THEN
      ! Post-collision momentum
      part1%part_p = u_cm*m1 + p_scat
    END IF

    IF (ran_w <= collision%w1_ratio) THEN

      ! Post-collision momentum
      p_2 = u_cm*m2 - p_scat

      IF (excited_id>0) THEN !Move particle to excited particle list
        new_part => part2
        ! Link part1 pointer to next colliding particle
        IF (.NOT.ASSOCIATED(p_list%head, TARGET=part2)) THEN
          part2 => part2%prev
        ELSE
          part2 => part2%next
        END IF
        new_part%part_p = p_2
        CALL remove_particle_from_partlist(p_list, new_part)
        CALL add_particle_to_partlist( &
          species_list(excited_id)%attached_list, new_part)
        NULLIFY(new_part)

      ELSE ! Just carry on with the same particle
        part2%part_p = p_2
      END IF

    END IF

    CALL link_particle_pointers(collision, part1, part2)
    NULLIFY(part1, part2)

  END SUBROUTINE bird_excitation



  SUBROUTINE bird_charge_exchange(collision)

    ! Collision process: backscattering
    ! Elastic 180 degrees scattering

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num) :: ran_w
    REAL(num) :: reducedm, m
    REAL(num), DIMENSION(3) :: g, p_scat, u_cm

    ran_w = random()
    reducedm = collision%reducedm
    g = collision%g
    u_cm = collision%u_cm

    ! Post-collision momentum
    p_scat = reducedm * g

    IF (ran_w <= collision%w2_ratio) THEN
      m = collision%m1
      collision%part1%part_p = u_cm*m - p_scat
    END IF

    IF (ran_w <= collision%w1_ratio) THEN
      m = collision%m2
      collision%part2%part_p = u_cm*m + p_scat
    END IF
  END SUBROUTINE bird_charge_exchange

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                 BOYD METHOD                                 !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE boyd_elastic_scattering(collision)

    ! Collision process: electron scattering
    ! e + N -> e + N

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num), DIMENSION(3) :: p_1, p_2, p_precoll
    REAL(num), DIMENSION(3) :: u_cm, random_direction
    REAL(num) :: w_ratio
    REAL(num) :: m1, m2, im1, im2, reducedm
    REAL(num) :: w1, w2, w_diff
    REAL(num) :: g_mag
    INTEGER :: species1, species2
    TYPE(particle), POINTER :: impact, current

    CALL set_particle_properties(collision, species1, species2, &
      m1, im1, m2, im2, w1, w2, g_mag, u_cm, current, impact)
    reducedm = collision%reducedm

    random_direction = random_unit_vector()
    random_direction = reducedm * g_mag * random_direction

    w_diff = w1 - w2
    IF ( ABS(w_diff) < TINY(0._num)) THEN
      ! Add energy add energy loss
      CALL release_energy_loss(collision, u_cm)
      
      ! Post collision momenta
      p_1 = u_cm*m1 + random_direction
      p_2 = u_cm*m2 - random_direction
      current%part_p = p_1
      impact%part_p = p_2

    ELSE IF ( w_diff < 0._num ) THEN ! w2 > w1
      w_ratio = w1/w2
      p_precoll = impact%part_p
      ! Post collision momenta
      p_1 = u_cm*m1 + random_direction
      p_2 = u_cm*m2 - random_direction
      current%part_p = p_1
      impact%part_p = p_precoll*(1._num - w_ratio) + p_2*w_ratio
      ! Energy loss
      CALL save_energy_loss(collision%ix, p_2, p_precoll, w_ratio, im2)

    ELSE IF ( w_diff > 0._num ) THEN ! w1 > w2
      w_ratio = w2/w1
      p_precoll = current%part_p
      ! Post collision momenta
      p_1 = u_cm*m1 + random_direction
      p_2 = u_cm*m2 - random_direction
      current%part_p = p_precoll*(1._num - w_ratio) + p_1*w_ratio
      impact%part_p = p_2
      CALL save_energy_loss(collision%ix, p_1, p_precoll, w_ratio, im1)

    END IF

    CALL link_particle_pointers(collision, current, impact)

  END SUBROUTINE boyd_elastic_scattering



  SUBROUTINE boyd_excitation(collision)

    ! Collision process: excitation
    ! e + N -> e + N*

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num), DIMENSION(3) :: p_1, p_2, p_precoll
    REAL(num), DIMENSION(3) :: u_cm, random_direction
    REAL(num) :: w_ratio
    REAL(num) :: m1, m2, im1, im2, reducedm, ireducedm
    REAL(num) :: w1, w2, w_diff
    REAL(num) :: g_mag, e_threshold
    INTEGER :: species1, species2
    TYPE(particle), POINTER :: impact, current

    CALL set_particle_properties(collision, species1, species2, &
      m1, im1, m2, im2, w1, w2, g_mag, u_cm, current, impact)
    reducedm = collision%reducedm
    ireducedm = collision%ireducedm
    e_threshold = collision%type_block%ethreshold

    !Random vector
    random_direction = random_unit_vector()
    g_mag = SQRT(g_mag*g_mag - 2._num*e_threshold*ireducedm)
    random_direction = reducedm * g_mag * random_direction

    w_diff = w1 - w2
    IF ( ABS(w_diff) < TINY(0._num)) THEN
      ! Add energy add energy loss
      CALL release_energy_loss(collision, u_cm)

      ! Post collision momenta
      p_1 = u_cm*m1 + random_direction
      p_2 = u_cm*m2 - random_direction
      current%part_p = p_1
      impact%part_p = p_2

    ELSE IF ( w_diff < 0._num ) THEN ! w2 > w1
      w_ratio = w1/w2
      p_precoll = impact%part_p
      ! Post collision momenta
      p_1 = u_cm*m1 + random_direction
      p_2 = u_cm*m2 - random_direction
      current%part_p = p_1
      impact%part_p = p_precoll*(1._num - w_ratio) + p_2*w_ratio
      ! Energy loss
      CALL save_energy_loss(collision%ix, p_2, p_precoll, w_ratio, im2)

    ELSE IF ( w_diff > 0._num ) THEN ! w1 > w2
      w_ratio = w2/w1
      p_precoll = current%part_p
      ! Post collision momenta
      p_1 = u_cm*m1 + random_direction
      p_2 = u_cm*m2 - random_direction
      current%part_p = p_precoll*(1._num - w_ratio) + p_1*w_ratio
      impact%part_p = p_2
      CALL save_energy_loss(collision%ix, p_1, p_precoll, w_ratio, im1)

    END IF

    CALL link_particle_pointers(collision, current, impact)

  END SUBROUTINE boyd_excitation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                SPLIT METHOD                                 !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef PER_SPECIES_WEIGHT
  SUBROUTINE split_ionisation(collision)

    ! Collision process: electron impact ionisation
    ! e + N -> 2e + I
    ! species1 = electron (e)
    ! species2 = neutral (N)

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num), DIMENSION(3) :: u_cm, u_e1, u_e2, random_direction
    REAL(num) :: m1, m2, im1, im2, w1, w2, reducedm
    REAL(num) :: ran_e
    REAL(num) :: g_mag, e_threshold, e_excess, u_excess, u_excess_total
    INTEGER :: species1, species2, ion_id
    TYPE(particle), POINTER :: part1, part2, new_part
    TYPE(particle_list), POINTER :: p_list

    CALL set_particle_properties(collision, species1, species2, &
      m1, im1, m2, im2, w1, w2, g_mag, u_cm, part1, part2, p_list)
    reducedm = collision%reducedm
    e_threshold = collision%type_block%ethreshold


    ! Excess energy distribution between the two electrons
    e_excess = 0.5_num*reducedm*g_mag*g_mag - e_threshold

    ! Speed associated to the energy excess (electrons)
    u_excess_total = SQRT(2._num*im1*e_excess)
    ran_e = 0.5_num ! 0.5 for equal energy distribution between electrons,
                    ! Otherwise: random()
    u_excess = u_excess_total*SQRT(ran_e)

    ! Impact electron
    random_direction = random_unit_vector() ! Random vector
    u_e1 = u_excess*random_direction ! cm frame of reference
    part1%part_p = (u_cm + u_e1)*m1


    ! Neutral (part2 pointer) is split into two
    ! - Ion and electron of weight w1
    ! - Untouched neutral of weight w2 - w1

    ! New electron
    CALL create_particle(new_part)
    u_excess = u_excess_total*SQRT(1._num - ran_e)
    random_direction = random_unit_vector() ! New random vector
    u_e2 = u_excess*random_direction  ! cm frame of reference
    new_part%part_p = (u_cm + u_e2)*m1
    new_part%part_pos = part2%part_pos
    new_part%weight = w1 ! Electron's weight
    CALL add_particle_to_partlist( &
      species_list(species1)%attached_list, new_part)
    NULLIFY(new_part)
    
    ! New ion
    ion_id = collision%type_block%new_species_id ! ions
    CALL create_particle(new_part)
    new_part%part_p = u_cm * (m2 - m1) - m1*(u_e1 + u_e2)
    new_part%part_pos = part2%part_pos
    new_part%weight = w1 ! Electron's weight
    CALL add_particle_to_partlist(species_list(ion_id)%attached_list, &
      new_part)
    NULLIFY(new_part)

    ! Neutral: substract ion's weight from neutral particle
    part2%weight = w2 - w1
    CALL check_weight_difference(w2, w1, part2, p_list)

    CALL link_particle_pointers(collision, part1, part2)
    NULLIFY(part1, part2)

  END SUBROUTINE split_ionisation



  SUBROUTINE split_charge_exchange(collision)

    ! Collision process: charge exchange
    ! Assumes that: N and I have the same mass (same species)
    !               N's weight > I's weight
    ! I + N -> N + I
    ! species1 -> ion
    ! species2 -> neutral

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num) :: m1, im1, m2, im2, w1, w2, g_mag
    REAL(num), DIMENSION(3) :: u_cm
    INTEGER :: species1, species2
    TYPE(particle), POINTER :: part1, part2, new_part
    TYPE(particle_list), POINTER :: p_list

    CALL set_particle_properties(collision, species1, species2, &
      m1, im1, m2, im2, w1, w2, g_mag, u_cm, part1, part2, p_list)

    ! Create new neutral (former ion)
    CALL create_particle(new_part)
    new_part%part_p = 2._num * m2 * u_cm - part1%part_p
    new_part%part_pos = part1%part_pos
    new_part%weight = w1 ! Ion's weight
    CALL add_particle_to_partlist( &
      species_list(species2)%attached_list, new_part)
    NULLIFY(new_part)

    ! Neutral becomes ion
    part1%part_p = 2._num * m1 * u_cm - part2%part_p
    part1%part_pos = part2%part_pos

    ! High weight neutral: substract ions weight
    part2%weight = w2 - w1
    CALL check_weight_difference(w2, w1, part2, p_list)

    CALL link_particle_pointers(collision, part1, part2)
    NULLIFY(part1, part2)

  END SUBROUTINE split_charge_exchange
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                 BACKGROUND                                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE bird_elastic_scattering_bg(collision)

    ! Collision process: electron scattering
    ! e + N -> e + N

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num) :: reducedm, g_mag, m1
    REAL(num), DIMENSION(3) :: random_direction, p_scat, u_cm

    reducedm = collision%reducedm
    g_mag = collision%g_mag
    u_cm = collision%u_cm
    m1 = collision%m1

    ! Post-collision momentum
    random_direction = random_unit_vector()
    p_scat = reducedm * g_mag * random_direction
    collision%part1%part_p = u_cm*m1 + p_scat 

  END SUBROUTINE bird_elastic_scattering_bg



  SUBROUTINE bird_ionisation_bg(collision)

    ! Collision process: electron impact ionisation
    ! e + N -> 2e + I
    ! species1 = electron (e)
    ! species2 = neutral background (N)

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num), DIMENSION(3) :: u_cm, u_e1, u_e2, random_direction
    REAL(num) :: m1, m2, im1,im2, reducedm
#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: w1
#endif
    REAL(num) :: ran_e
    REAL(num) :: g_mag, e_threshold, e_excess, u_excess, u_excess_total
    INTEGER :: species1, ion_id
    TYPE(particle), POINTER :: part1, new_part

    m1 = collision%m1
    im1 = collision%im1
    m2 = collision%m2
    im2 = collision%im2
    u_cm = collision%u_cm
#ifndef PER_SPECIES_WEIGHT
    w1 = collision%w1
#endif
    g_mag = collision%g_mag
    reducedm = collision%reducedm
    e_threshold = collision%type_block%ethreshold
    species1 = collision%species1
    part1 => collision%part1


    ! Excess energy distribution between the two electrons
    e_excess = 0.5_num*reducedm*g_mag*g_mag - e_threshold

    ! Speed associated to the energy excess (electrons)
    u_excess_total = SQRT(2._num*im1*e_excess)
    ran_e = 0.5_num ! 0.5 for equal energy distribution between electrons,
                    ! Otherwise: random()
    u_excess = u_excess_total*SQRT(ran_e)

    ! Impact electron
    random_direction = random_unit_vector() ! Random unitary vector
    u_e1 = u_excess*random_direction ! cm frame of reference
    part1%part_p = (u_cm + u_e1)*m1

    ! Background generates: ion and electron of weight = w1
    ! New electron
    CALL create_particle(new_part)
    u_excess = u_excess_total*SQRT(1._num - ran_e)
    random_direction = random_unit_vector() ! Random unitary vector
    u_e2 = u_excess*random_direction  ! cm frame of reference
    new_part%part_p = (u_cm + u_e2)*m1
    new_part%part_pos = part1%part_pos
#ifndef PER_SPECIES_WEIGHT
    new_part%weight = w1 ! Electron's weight
#endif
    CALL add_particle_to_partlist( &
      species_list(species1)%attached_list, new_part)
    NULLIFY(new_part)
    
    ! New ion
    ion_id = collision%type_block%new_species_id ! ions
    CALL create_particle(new_part)
    new_part%part_p = u_cm * (m2 - m1) - m1*(u_e1 + u_e2)
    new_part%part_pos = part1%part_pos
#ifndef PER_SPECIES_WEIGHT
    new_part%weight = w1 ! Electron's weight
#endif
    CALL add_particle_to_partlist(species_list(ion_id)%attached_list,new_part)
    NULLIFY(new_part, part1)

  END SUBROUTINE bird_ionisation_bg



  SUBROUTINE bird_excitation_bg(collision)

    ! Collision process: excitation
    ! e + N -> e + N*

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num) :: g_mag, e_threshold
    REAL(num) :: m1, m2, reducedm, ireducedm
    REAL(num), DIMENSION(3) :: random_direction, u_cm, p_scat
    INTEGER :: excited_id
    TYPE(particle), POINTER :: part1, new_part

    m1 = collision%m1
    reducedm = collision%reducedm
    ireducedm = collision%ireducedm
    u_cm = collision%u_cm
    g_mag = collision%g_mag
    e_threshold = collision%type_block%ethreshold
    excited_id = collision%type_block%new_species_id
    part1 => collision%part1

    !Random vector
    random_direction = random_unit_vector()
    g_mag = SQRT(g_mag*g_mag - 2._num*e_threshold*ireducedm)
    p_scat = reducedm * g_mag * random_direction

    ! Post-collision momentum
    part1%part_p = u_cm*m1 + p_scat

    IF (excited_id>0) THEN !Create particle to excited particle list

      m2 = collision%m2
      ! Post-collision momentum
      CALL create_particle(new_part)
      new_part%part_pos = part1%part_pos
      new_part%part_p = u_cm*m2 - p_scat
#ifndef PER_SPECIES_WEIGHT
      new_part%weight = collision%w1
#endif
      CALL add_particle_to_partlist( &
        species_list(excited_id)%attached_list, new_part)
      NULLIFY(new_part)

    END IF

    NULLIFY(part1)

  END SUBROUTINE bird_excitation_bg



  SUBROUTINE bird_charge_exchange_bg(collision)

    ! Collision process: charge_exchange
    ! I + N -> N + I
    ! Where the ions (I) are species1 and neutrals (N) are background

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num) :: reducedm, m1
    REAL(num), DIMENSION(3) :: g, p_scat, u_cm

    reducedm = collision%reducedm
    g = collision%g
    u_cm = collision%u_cm
    m1 = collision%m1

    ! Post-collision momentum
    p_scat = reducedm * g
    collision%part1%part_p = u_cm*m1 - p_scat 

  END SUBROUTINE bird_charge_exchange_bg



  SUBROUTINE vahedi_electron_elastic_scattering_bg(collision)

    ! Collision process: electron elastic scattering
    ! e + N -> e + N
    REAL(num), DIMENSION(3) :: v_inc, v_scat, v_inc_i
    REAL(num) :: costheta, sintheta, coschi, sinchi, cosphi, sinphi, phi
    REAL(num) :: e_inc, delta_e, e_scat, g_scat, g_mag
    REAL(num) :: sinratio, m1, m2

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    ! Incoming normalised velocity vector
    v_inc = vector_normalisation(collision%g)

    !Theta angle
    costheta = v_inc(1)
    sintheta = SQRT(1._num - costheta*costheta)
    ! Chi angle
    coschi = 1._num - 2._num * random()
    sinchi = SQRT(1._num - coschi*coschi)
    ! Phi angle
    phi = 2._num * pi * random()
    cosphi = COS(phi)
    sinphi = SIN(phi)

    ! Scattered normalised vector
    sinratio = sinchi / sintheta
    v_inc_i = crossproduct(v_inc,(/1._num, 0._num, 0._num/))
    v_scat = v_inc * coschi + v_inc_i * sinratio * sinphi + &
      crossproduct(v_inc_i,v_inc) * sinratio * cosphi

    m1 = collision%m1
    m2 = collision%m2
    g_mag = collision%g_mag
    e_inc = 0.5_num * m1 * g_mag * g_mag
    delta_e = 2._num * m1 / m2 * (1._num - coschi)
    e_scat = e_inc * (1._num - delta_e)
    g_scat = SQRT(2._num * e_scat / m1)

    ! Post-collision momentum
    collision%part1%part_p = v_scat * g_scat * m1

  END SUBROUTINE vahedi_electron_elastic_scattering_bg



  SUBROUTINE vahedi_ion_elastic_scattering_bg(collision)

    ! Collision process: ion elastic scattering
    ! e + N -> e + N
    REAL(num), DIMENSION(3) :: v_inc, v_scat, v_inc_i
    REAL(num) :: costheta, sintheta, coschi, sinchi, cosphi, sinphi, phi
    REAL(num) :: e_inc, e_scat, g_scat, g_mag
    REAL(num) :: sinratio, m1

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    ! Incoming normalised velocity vector
    v_inc = vector_normalisation(collision%g)

    !Theta angle
    costheta = v_inc(1)
    sintheta = SQRT(1._num - costheta*costheta)
    ! Chi angle
    coschi = SQRT(1._num - random())
    sinchi = SQRT(1._num - coschi*coschi)
    ! Phi angle
    phi = 2._num * pi * random()
    cosphi = COS(phi)
    sinphi = SIN(phi)

    ! Scattered normalised vector
    sinratio = sinchi / sintheta
    v_inc_i = crossproduct(v_inc,(/1._num, 0._num, 0._num/))
    v_scat = v_inc * coschi + v_inc_i * sinratio * sinphi + &
      crossproduct(v_inc_i,v_inc) * sinratio * cosphi

    mu = collision%reducedm
    g_mag = collision%g_mag
    e_inc = 0.5_num * mu * g_mag * g_mag
    e_scat = e_inc * coschi * coschi
    g_scat = SQRT(2._num * e_scat / mu)

    ! Post-collision momentum
    collision%part1%part_p = collision%u_cm*collision%m1 + v_scat*g_scat*mu

  END SUBROUTINE vahedi_ion_elastic_scattering_bg



  SUBROUTINE vahedi_excitation_bg(collision)

    ! Collision process: electron impact excitation
    ! e + N -> e + N
    REAL(num), DIMENSION(3) :: v_inc, v_scat, v_inc_i
    REAL(num) :: costheta, sintheta, coschi, sinchi, cosphi, sinphi, phi
    REAL(num) :: e_inc, delta_e, e_scat, g_scat, g_mag
    REAL(num) :: sinratio, m1

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    ! Incoming normalised velocity vector
    v_inc = vector_normalisation(collision%g)

    !Theta angle
    costheta = v_inc(1)
    sintheta = SQRT(1._num - costheta*costheta)
    ! Chi angle
    coschi = 1._num - 2._num * random()
    sinchi = SQRT(1._num - coschi*coschi)
    ! Phi angle
    phi = 2._num * pi * random()
    cosphi = COS(phi)
    sinphi = SIN(phi)

    ! Scattered normalised vector
    sinratio = sinchi / sintheta
    v_inc_i = crossproduct(v_inc,(/1._num, 0._num, 0._num/))
    v_scat = v_inc * coschi + v_inc_i * sinratio * sinphi + &
      crossproduct(v_inc_i,v_inc) * sinratio * cosphi

    m1 = collision%m1
    g_mag = collision%g_mag
    e_inc = 0.5_num * m1 * g_mag * g_mag
    delta_e = collision%type_block%ethreshold
    e_scat = e_inc - delta_e
    g_scat = SQRT(2._num * e_scat / m1)

    ! Post-collision momentum
    collision%part1%part_p = v_scat * g_scat * m1

  END SUBROUTINE vahedi_excitation_bg



  SUBROUTINE vahedi_ionisation_bg(collision)

    ! Collision process: electron impact ionisation
    ! e + N -> e + N
    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision
    REAL(num), DIMENSION(3) :: v_inc, v_scat, v_inc_i
    REAL(num) :: costheta, sintheta, coschi, sinchi, cosphi, sinphi, phi
    REAL(num) :: e_inc, delta_e, e_scat, g_scat, g_mag
    REAL(num) :: sinratio, m1, part_pos
    INTEGER :: species_id

    TYPE(particle), POINTER :: new_part

    part_pos = collision%part1%part_pos

    m1 = collision%m1
    g_mag = collision%g_mag
    e_inc = 0.5_num * m1 * g_mag * g_mag
    delta_e = collision%type_block%ethreshold
    e_scat = (e_inc - delta_e) * 0.5_num
    g_scat = SQRT(2._num * e_scat / m1)

    ! Incoming normalised velocity vector
    v_inc = vector_normalisation(collision%g)
    v_inc_i = crossproduct(v_inc,(/1._num, 0._num, 0._num/))
    !Theta angle
    costheta = v_inc(1)
    sintheta = SQRT(1._num - costheta*costheta)

    ! Electron #1:
    ! Chi angle
    coschi = 1._num - 2._num * random()
    sinchi = SQRT(1._num - coschi*coschi)
    ! Phi angle
    phi = 2._num * pi * random()
    cosphi = COS(phi)
    sinphi = SIN(phi)
    ! Scattered normalised vector
    sinratio = sinchi / sintheta
    v_scat = v_inc * coschi + v_inc_i * sinratio * sinphi + &
      crossproduct(v_inc_i,v_inc) * sinratio * cosphi
    ! Post-collision momentum
    collision%part1%part_p = v_scat * g_scat * m1


    ! Electron #2:
    ! Chi angle
    coschi = 1._num - 2._num * random()
    sinchi = SQRT(1._num - coschi*coschi)
    ! Phi angle
    phi = 2._num * pi * random()
    cosphi = COS(phi)
    sinphi = SIN(phi)
    ! Scattered normalised vector
    sinratio = sinchi / sintheta
    v_scat = v_inc * coschi + v_inc_i * sinratio * sinphi + &
      crossproduct(v_inc_i,v_inc) * sinratio * cosphi
    ! Create a new electron particle
    species_id = collision%species1
    CALL create_particle(new_part)
    new_part%part_p = v_scat * g_scat * m1
    new_part%part_pos = part_pos
#ifndef PER_SPECIES_WEIGHT
    new_part%weight = collision%w1 ! Electron's weight
#endif
    CALL add_particle_to_partlist(species_list(species_id)%attached_list, &
      new_part)
    NULLIFY(new_part)


    ! Ion
    ! Create a new ion particle and assign neutral particle's momentum
    species_id = collision%type_block%new_species_id ! ions
    CALL create_particle(new_part)
    new_part%part_p = collision%u_2 * collision%m2
    new_part%part_pos = part_pos
#ifndef PER_SPECIES_WEIGHT
    new_part%weight = w1 ! Electron's weight
#endif
    CALL add_particle_to_partlist(species_list(species_id)%attached_list, &
      new_part)
    NULLIFY(new_part)

  END SUBROUTINE vahedi_ionisation_bg



  FUNCTION crossproduct(a, b)
    ! Computes: a x b = crossproduct
    REAL(num), DIMENSION(3) :: crossproduct
    REAL(num),  DIMENSION(3), INTENT(IN) :: a, b

    crossproduct(1) = a(2) * b(3) - a(3) * b(2)
    crossproduct(2) = a(3) * b(1) - a(1) * b(3)
    crossproduct(3) = a(1) * b(2) - a(2) * b(1)

  END FUNCTION crossproduct


  FUNCTION vector_normalisation(vec) RESULT(vec_norm)
    ! Normalised the input vector
    REAL(num), DIMENSION(3) :: vec, vec_norm
    REAL(num) :: norm

    norm = SQRT(DOT_PRODUCT(vec,vec))
    vec_norm = vec / norm

  END FUNCTION vector_normalisation


END MODULE nc_subroutines
