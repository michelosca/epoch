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
#ifdef NEUTRAL_COLLISIONS

  USE shared_data
  USE partlist
  USE nc_auxiliary
  USE random_generator

  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: nanbu_elastic_scattering, nanbu_ionisation, nanbu_excitation
  PUBLIC :: nanbu_charge_exchange, nanbu_excitation_bg
  PUBLIC :: nanbu_elastic_scattering_bg, nanbu_ionisation_bg
  PUBLIC :: nanbu_charge_exchange_bg
  PUBLIC :: vahedi_electron_elastic_scattering_bg, vahedi_excitation_bg
  PUBLIC :: vahedi_ion_elastic_scattering_bg, vahedi_ionisation_bg
  PUBLIC :: vahedi_electron_elastic_scattering, vahedi_excitation
  PUBLIC :: vahedi_ionisation, vahedi_ion_elastic_scattering
#ifndef PER_SPECIES_WEIGHT
  PUBLIC :: vahedi_split_ionisation!, vahedi_split_charge_exchange
  PUBLIC :: nanbu_split_ionisation, nanbu_split_charge_exchange
#endif
CONTAINS


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



  SUBROUTINE set_particle_properties(collision, species1, species2, &
    m1, im1, m2, im2, w1, w2, w1_ratio, w2_ratio, g, part1, part2, &
#ifndef PER_SPECIES_WEIGHT
    p_list1, p_list2, w1max_rat, w2max_rat)
#else
    p_list1, p_list2 )
#endif

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision
    REAL(num), INTENT(OUT), OPTIONAL :: m1, im1, m2, im2
    REAL(num), INTENT(OUT), OPTIONAL :: w1, w2, w1_ratio, w2_ratio
#ifndef PER_SPECIES_WEIGHT
    REAL(num), INTENT(OUT), OPTIONAL :: w1max_rat, w2max_rat
#endif
    REAL(num), DIMENSION(3), INTENT(OUT), OPTIONAL :: g
    INTEGER, INTENT(OUT), OPTIONAL :: species1, species2
    TYPE(particle), POINTER, INTENT(OUT), OPTIONAL :: part1, part2
    TYPE(particle_list), POINTER, INTENT(OUT), OPTIONAL :: p_list1, p_list2

    IF (collision%species1 == collision%type_block%source_species_id) THEN
      IF (PRESENT(species1)) species1 = collision%species2
      IF (PRESENT(species2)) species2 = collision%species1
      IF (PRESENT(m1)) m1 = collision%m2
      IF (PRESENT(im1))im1 = collision%im2
      IF (PRESENT(m2)) m2 = collision%m1
      IF (PRESENT(im2)) im2 = collision%im1
      IF (PRESENT(w1)) w1 = collision%w2
      IF (PRESENT(w2)) w2 = collision%w1
      IF (PRESENT(w1_ratio)) w1_ratio = collision%w2_ratio
      IF (PRESENT(w2_ratio)) w2_ratio = collision%w1_ratio
#ifndef PER_SPECIES_WEIGHT
      IF (PRESENT(w1max_rat)) &
        w1max_rat = collision%w2/collision%collision_block%max_w2
      IF (PRESENT(w2max_rat)) &
        w2max_rat = collision%w1/collision%collision_block%max_w1
#endif
      IF (PRESENT(part1)) part1 => collision%part2
      IF (PRESENT(part2)) part2 => collision%part1
      IF (PRESENT(p_list2)) p_list2 => collision%p_list1
      IF (PRESENT(p_list1)) p_list1 => collision%p_list2
      IF (PRESENT(g)) g = -collision%g
    ELSE
      IF (PRESENT(species1)) species1 = collision%species1
      IF (PRESENT(species2)) species2 = collision%species2
      IF (PRESENT(m1)) m1 = collision%m1
      IF (PRESENT(im1)) im1 = collision%im1
      IF (PRESENT(m2)) m2 = collision%m2
      IF (PRESENT(im2)) im2 = collision%im2
      IF (PRESENT(w1)) w1 = collision%w1
      IF (PRESENT(w2)) w2 = collision%w2
      IF (PRESENT(w1_ratio)) w1_ratio = collision%w1_ratio
      IF (PRESENT(w2_ratio)) w2_ratio = collision%w2_ratio
#ifndef PER_SPECIES_WEIGHT
      IF (PRESENT(w2max_rat)) &
        w2max_rat = collision%w2/collision%collision_block%max_w2
      IF (PRESENT(w1max_rat)) &
        w1max_rat = collision%w1/collision%collision_block%max_w1
#endif
      IF (PRESENT(part1)) part1 => collision%part1
      IF (PRESENT(part2)) part2 => collision%part2
      IF (PRESENT(p_list1)) p_list1 => collision%p_list1
      IF (PRESENT(p_list2)) p_list2 => collision%p_list2
      IF (PRESENT(g)) g = collision%g
    END IF

  END SUBROUTINE set_particle_properties


  SUBROUTINE link_particle_pointers(collision, part1, part2)

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision
    TYPE(particle), POINTER, INTENT(INOUT) :: part1, part2

    IF (collision%species1 == collision%type_block%source_species_id) THEN
      collision%part1 => part2
      collision%part2 => part1
    ELSE
      collision%part1 => part1
      collision%part2 => part2
    ENDIF

    NULLIFY(part1, part2)
  END SUBROUTINE link_particle_pointers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                 NANBU METHOD                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE nanbu_elastic_scattering(collision)

    ! Collision process: electron scattering
    ! e + N -> e + N

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num) :: reducedm, g_mag, m, ran1
    REAL(num), DIMENSION(3) :: random_direction, p_scat, u_cm

    reducedm = collision%reducedm
    g_mag = collision%g_mag
    u_cm = cm_velocity(collision)

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

  END SUBROUTINE nanbu_elastic_scattering



  SUBROUTINE nanbu_excitation(collision)

    ! Collision process: excitation
    ! e + N -> e + N*

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    INTEGER :: excited_id
    REAL(num) :: g_mag, e_threshold
    REAL(num) :: m1, m2, reducedm, ireducedm
    REAL(num) :: ran_w, w1rat, w2rat
    REAL(num), DIMENSION(3) :: random_direction, p_2, u_cm, p_scat
    TYPE(particle), POINTER :: part1, part2, new_part
    TYPE(particle_list), POINTER :: p_list2

    CALL set_particle_properties(collision, &
      m1 = m1, m2 = m2, w1_ratio = w1rat, w2_ratio = w2rat, &
      part1 = part1, part2 = part2, p_list2 = p_list2)
    excited_id = collision%type_block%new_species_id
    reducedm = collision%reducedm
    ireducedm = collision%ireducedm
    e_threshold = collision%type_block%ethreshold
    g_mag = collision%g_mag
    u_cm = cm_velocity(collision) 

    !Random vector
    random_direction = random_unit_vector()
    g_mag = SQRT(g_mag*g_mag - 2._num*e_threshold*ireducedm)
    p_scat = reducedm * g_mag * random_direction

    ran_w = random()
    IF (ran_w <= w2rat) THEN
      ! Post-collision momentum
      part1%part_p = u_cm*m1 + p_scat
    END IF

    IF (ran_w <= w1rat) THEN

      ! Post-collision momentum
      p_2 = u_cm*m2 - p_scat

      IF (excited_id>0) THEN !Move particle to excited particle list
        new_part => part2
        ! Link part2 pointer to next colliding particle
        IF (.NOT.ASSOCIATED(p_list2%head, TARGET=part2)) THEN
          part2 => part2%prev
        ELSE
          part2 => part2%next
        END IF
        new_part%part_p = p_2
        CALL remove_particle_from_partlist(p_list2, new_part)
        CALL add_particle_to_partlist( &
          species_list(excited_id)%attached_list, new_part)
        NULLIFY(new_part)

      ELSE ! Just carry on with the same particle
        part2%part_p = p_2
      END IF

    END IF

    CALL link_particle_pointers(collision, part1, part2)
    NULLIFY(p_list2)

  END SUBROUTINE nanbu_excitation


#ifndef PER_SPECIES_WEIGHT
  SUBROUTINE merge_particles(part1, part2, p_list1)

    ! Merge particles 1 and 2 and remove particle 1
    TYPE(particle), POINTER, INTENT(INOUT) :: part1, part2
    TYPE(particle_list), POINTER, INTENT(INOUT) :: p_list1
    REAL(num), DIMENSION(3) :: p
    REAL(num) :: w1, w2, wtot

    ! Merge current to particle2
    w1 = part1%weight
    w2 = part2%weight
    wtot = w1 + w2
    ! Average momentum
    p = (part1%part_p*w1 + part2%part_p*w2 ) / wtot
    part2%part_p = p
    part2%weight = wtot
    CALL remove_particle_from_partlist(p_list1, part1)
    CALL destroy_particle(part1)

  END SUBROUTINE merge_particles



  SUBROUTINE nanbu_split_ionisation(collision)

    ! Collision process: electron impact ionisation
    ! e + N -> 2e + I
    ! species1 = electron (e)
    ! species2 = neutral (N)

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num), DIMENSION(3) :: u_cm, u_e1, u_e2, random_direction
    REAL(num) :: m1, m2, im1, w1, w2, w_diff, w_min, reducedm
    REAL(num) :: ran_e, part2_pos
    REAL(num) :: g_mag, e_threshold, e_excess, u_excess, u_excess_total
    INTEGER :: species1, ion_id
    LOGICAL :: merge_electrons
    TYPE(particle), POINTER :: part1, part2, new_part, part_next
    TYPE(particle_list), POINTER :: p_list1, p_list2

    CALL set_particle_properties(collision, species1 = species1, &
      m1 = m1, im1 = im1, m2 = m2, w1 = w1, w2 = w2, part1 = part1, &
      part2 = part2, p_list1 = p_list1, p_list2 = p_list2)
    g_mag = collision%g_mag
    u_cm = cm_velocity(collision) 
    reducedm = collision%reducedm
    e_threshold = collision%type_block%ethreshold
    part2_pos = part2%part_pos

    ! Excess energy distribution between the two electrons
    e_excess = 0.5_num*reducedm*g_mag*g_mag - e_threshold
    ! Speed associated to the energy excess (electrons)
    u_excess_total = SQRT(2._num*im1*e_excess)
    ran_e = 0.5_num ! 0.5 for equal energy distribution between electrons,
                    ! Otherwise: random()

    ! Check particle's weight
    merge_electrons = .FALSE.
    IF (w1 >= w2) THEN
      ! Electron's weight is larger than neutral
      w_min = w2
      w_diff = w1 - w2
      !Update impact electron weight
      part1%weight = w_diff
      !Remove neutral particle from list
      new_part => part2%next ! This is just a buffer
      CALL remove_particle_from_partlist(p_list2, part2)
      CALL destroy_particle(part2)
      part2 => new_part
      NULLIFY(new_part)
      IF (w_diff < w1) merge_electrons = .TRUE.
    ELSE
      ! Neutral's weight is larger
      w_min = w1
      w_diff = w2 - w1
      part2%weight = w_diff
      IF (w_diff < w1) THEN
        ! If w_diff is small merge neutral with next particle
        new_part => part2 ! Part to be merged and removed
        part2 => part2%next
        CALL merge_particles(new_part, part2, p_list2)
        NULLIFY(new_part)
      END IF
    ENDIF

    ! Impact electron
    random_direction = random_unit_vector() ! Random vector
    u_excess = u_excess_total*SQRT(ran_e)
    u_e1 = u_excess*random_direction ! cm frame of reference
    part1%part_p = (u_cm + u_e1)*m1

    ! Neutral (part2 pointer) is split into
    ! - New ion
    ! - New electron
    ! - Untouched neutral

    ! New electron
    CALL create_particle(new_part)
    random_direction = random_unit_vector() ! New random vector
    u_excess = u_excess_total*SQRT(1._num - ran_e)
    u_e2 = u_excess*random_direction  ! cm frame of reference
    new_part%part_p = (u_cm + u_e2)*m1
    new_part%part_pos = part2_pos
    new_part%weight = w_min
    CALL add_particle_to_partlist( &
      species_list(species1)%attached_list, new_part)
    IF (merge_electrons) THEN
      part_next => part1%next
      CALL merge_particles(part1, new_part, p_list1)
      part1 => part_next
      NULLIFY(part_next)
    END IF
    NULLIFY(new_part)

    ! New ion
    ion_id = collision%type_block%new_species_id ! ions
    CALL create_particle(new_part)
    new_part%part_p = u_cm * (m2 - m1) - m1*(u_e1 + u_e2)
    new_part%part_pos = part2_pos
    new_part%weight = w_min
    CALL add_particle_to_partlist(species_list(ion_id)%attached_list, &
      new_part)
    NULLIFY(new_part)

    CALL link_particle_pointers(collision, part1, part2)
    NULLIFY(p_list1, p_list2)

  END SUBROUTINE nanbu_split_ionisation



  SUBROUTINE nanbu_split_charge_exchange(collision)

    ! Collision process: charge exchange
    ! I + N -> N + I

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num) :: m1, m2, mu, w1, w2, w_diff
    REAL(num), DIMENSION(3) :: u_cm, p_scat
    TYPE(particle), POINTER :: part1, part2, new_part, part_next
    TYPE(particle_list), POINTER :: p_list

    w1 = collision%w1
    w2 = collision%w2
    m1 = collision%m1
    m2 = collision%m2
    part1 => collision%part1
    part2 => collision%part2
    u_cm = cm_velocity(collision) 
    mu = collision%reducedm
    p_scat = mu * collision%g

    ! Check particle's weight
    IF (w1 >= w2) THEN
      ! Particle 1 is split into two
      ! Particle 1: unperturbed (only weight loss)
      w_diff = w1 - w2
      part1%weight = w_diff
      ! Species 1: create new particle
      CALL create_particle(new_part)
      new_part%part_p = u_cm * m1 - p_scat
      new_part%part_pos = part2%part_pos
      new_part%weight = w2
      CALL add_particle_to_partlist( &
        species_list(collision%species1)%attached_list, new_part)
      ! Particle 2
      part2%part_p = u_cm * m2 + p_scat
      part2%part_pos = part1%part_pos
      ! Check weight difference
      IF (w_diff < w2) THEN
        ! If w_diff small merge particles of species1
        p_list => collision%p_list1
        part_next => part1%next
        CALL merge_particles(part1, new_part, p_list)
        part1 => part_next
        NULLIFY(p_list, part_next)
      END IF
      NULLIFY(new_part)
    ELSE
      ! Particle 2 is split into two
      ! Particle 1
      part1%part_p = u_cm * m1 - p_scat
      part1%part_pos = part2%part_pos
      ! Particle 2: unperturbed (only weight loss)
      w_diff = w2 - w1
      part2%weight = w_diff
      ! Species 2: create new particle
      CALL create_particle(new_part)
      new_part%part_p = u_cm * m2 + p_scat
      new_part%part_pos = part1%part_pos
      new_part%weight = w1
      CALL add_particle_to_partlist( &
        species_list(collision%species2)%attached_list, new_part)
      ! Check weight difference
      IF (w_diff < w1) THEN
        ! If w_diff small merge particles of species2
        p_list => collision%p_list2
        part_next => part2%next
        CALL merge_particles(part2, new_part, p_list)
        part2 => part_next
        NULLIFY(p_list, part_next)
      END IF
      NULLIFY(new_part)
    ENDIF

    collision%part1 => part1
    collision%part2 => part2
    NULLIFY(part1, part2)

  END SUBROUTINE nanbu_split_charge_exchange
#endif


  SUBROUTINE nanbu_ionisation(collision)

    ! Collision process: electron impact ionisation
    ! e + N -> 2e + I
    ! species1 = electron (e)
    ! species2 = neutral (N)

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num), DIMENSION(3) :: u_cm, u_e1, u_e2, random_direction
    REAL(num) :: m1, m2, im1, mu
#ifndef PER_SPECIES_WEIGHT
    REAL(num) ::  w2
#endif
    REAL(num) :: ran_e, ran_w, w1rat, w2rat
    REAL(num) :: g_mag, e_threshold, e_excess, u_excess, u_excess_total
    INTEGER :: ion_id, species1
    TYPE(particle), POINTER :: new_part, part1, part2
    TYPE(particle_list) , POINTER :: p_list2

    CALL set_particle_properties(collision, species1 = species1, &
      m1 = m1, im1 = im1, m2 = m2, w1_ratio = w1rat, w2_ratio = w2rat, &
#ifndef PER_SPECIES_WEIGHT
      w2 = w2, &
#endif
      part1 = part1, part2 = part2, p_list2 = p_list2)
    u_cm = cm_velocity(collision) 
    g_mag = collision%g_mag
    mu = collision%reducedm

    ! Excess energy distribution between the two electrons
    e_threshold = collision%type_block%ethreshold
    e_excess = 0.5_num*mu*g_mag*g_mag - e_threshold

    ! Speed associated to the energy excess
    u_excess_total = SQRT(2._num*im1*e_excess)
    ran_e = 0.5_num !random()

    ! Random weight ratio
    ran_w = random()

    ! Current particle: the existing electron
    random_direction = random_unit_vector()
    u_excess = u_excess_total*SQRT(ran_e)
    u_e1 = u_excess*random_direction ! cm frame of reference
    IF (ran_w <= w2rat) THEN
      part1%part_p = (u_cm + u_e1)*m1
    END IF

    ! Impact particle: neutral -> ion + electron
    IF (ran_w <= w1rat) THEN
      ! The new electron
      CALL create_particle(new_part)
      random_direction = random_unit_vector()
      u_excess = u_excess_total*SQRT(1._num - ran_e)
      u_e2 = u_excess*random_direction  ! cm frame of reference
      new_part%part_p = (u_cm + u_e2)*m1
      new_part%part_pos = part2%part_pos
#ifndef PER_SPECIES_WEIGHT
      new_part%weight = w2
#endif
      CALL add_particle_to_partlist( &
        species_list(species1)%attached_list, new_part)
      NULLIFY(new_part)

      ! New ion: move neutral particle to ion species list
      new_part => part2
      new_part%part_p = u_cm * (m2 - m1) - m1*(u_e1 + u_e2)
      IF (.NOT.ASSOCIATED(p_list2%head, TARGET=part2)) THEN
        ! If particle is not the head of the list
        part2 => part2%prev
      ELSE
        ! If particle is the head of the list then one particle is skipped
        p_list2%coll_counter = p_list2%coll_counter + 1
        part2 => part2%next
      END IF
      CALL remove_particle_from_partlist(p_list2, new_part)
      ion_id = collision%type_block%new_species_id ! ions
      CALL add_particle_to_partlist( species_list(ion_id)%attached_list, &
        new_part)
      NULLIFY(new_part)
    END IF

    CALL link_particle_pointers(collision, part1, part2)
    NULLIFY(p_list2)

  END SUBROUTINE nanbu_ionisation



  SUBROUTINE nanbu_charge_exchange(collision)

    ! Collision process: charge_exchange
    ! I + N -> N + I
    ! Where the ions (I) are species1 and neutrals (N) are background

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num) :: mu, m, ranw
    REAL(num), DIMENSION(3) :: g, p_scat, u_cm

    mu = collision%reducedm
    g = collision%g
    u_cm = cm_velocity(collision) 

    ! Post-collision momentum
    p_scat = mu * g

    ranw = random()
    IF (ranw <= collision%w2_ratio) THEN
      m = collision%m1
      collision%part1%part_p = u_cm*m - p_scat
    END IF

    IF (ranw <= collision%w1_ratio) THEN
      m = collision%m2
      collision%part2%part_p = u_cm*m + p_scat
    END IF

  END SUBROUTINE nanbu_charge_exchange

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                VAHEDI METHOD                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE vahedi_electron_elastic_scattering(collision)

    ! Collision process: electron elastic scattering
    ! e + N -> e + N
    REAL(num), DIMENSION(3) :: g, v_inc, v_scat, v_inc_i
    REAL(num) :: costheta, sintheta, coschi, sinchi, cosphi, sinphi, phi
    REAL(num) :: e_inc, delta_e, e_scat, g_scat_m1, g_mag
    REAL(num) :: sinratio, m1, im2
    REAL(num) :: ranw, w2rat
    TYPE(particle), POINTER :: part1

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    ! Electrons must be species 1, and neutrals => species2
    CALL set_particle_properties(collision, &
      m1 = m1, im2 = im2, w2_ratio = w2rat, &
      part1 = part1, g = g)

    ranw = random()
    IF (ranw <= w2rat) THEN
      ! Incoming normalised velocity vector
      v_inc = vector_normalisation(g)

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

      ! Post-collision speed (g_scat)
      g_mag = collision%g_mag
      e_inc = 0.5_num * m1 * g_mag * g_mag
      delta_e = 2._num * m1 * im2 * (1._num - coschi)
      e_scat = e_inc * (1._num - delta_e)
      g_scat_m1 = SQRT(2._num * e_scat * m1)

      ! Post-collision momentum
      part1%part_p = v_scat * g_scat_m1
    END IF
    !IF (ran1 <= collision%w1_ratio) THEN
    !  collision%part2%part_p =
    !END IF

  END SUBROUTINE vahedi_electron_elastic_scattering



  SUBROUTINE vahedi_ion_elastic_scattering(collision)

    ! Collision process: electron elastic scattering
    ! e + N -> e + N
    REAL(num), DIMENSION(3) :: g, v_inc, v_scat, v_inc_i, u_cm
    REAL(num) :: costheta, sintheta, coschi, sinchi, cosphi, sinphi, phi
    REAL(num) :: e_inc, e_scat, g_scat_m1, g_mag
    REAL(num) :: sinratio, m1, m12, ucm_mag2
    REAL(num) :: ranw, w2rat
    TYPE(particle), POINTER :: part1

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    ! Ions must be species 1, and neutrals => species2
    CALL set_particle_properties(collision, &
      m1 = m1, w2_ratio = w2rat, &
      part1 = part1, g = g)

    ranw = random()
    IF (ranw <= w2rat) THEN
      m12 = collision%m12

      ! Incoming normalised velocity vector
      v_inc = vector_normalisation(g)

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

      ! Post-collision speed (g_scat)
      u_cm = cm_velocity(collision)
      g_mag = collision%g_mag
      ucm_mag2 = DOT_PRODUCT(u_cm, u_cm)
      e_inc = 0.5_num * m1 * g_mag * g_mag
      e_scat = e_inc * coschi * coschi + 0.5_num * m12 * ucm_mag2
      g_scat_m1 = SQRT(2._num * e_scat * m1)

      ! Post-collision momentum
      part1%part_p = v_scat * g_scat_m1
    END IF
    !IF (ranw <= collision%w1_ratio) THEN
    !  collision%part2%part_p =
    !END IF

  END SUBROUTINE vahedi_ion_elastic_scattering



  SUBROUTINE vahedi_excitation(collision)

    ! Collision process: electron impact excitation
    ! e + N -> e + N

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num) :: g_mag, e_threshold
    REAL(num) :: e_inc, e_scat, g_scat_m1
    REAL(num) :: costheta, sintheta, coschi, sinchi, cosphi, sinphi, phi
    REAL(num) :: m1, w2rat
    REAL(num) :: ran_w, sinratio
    REAL(num), DIMENSION(3) :: g
    REAL(num), DIMENSION(3) :: v_inc, v_scat, v_inc_i
    TYPE(particle), POINTER :: part1

    CALL set_particle_properties(collision, &
      m1 = m1, w2_ratio = w2rat, &
      part1 = part1, g = g)

    ran_w = random()
    IF (ran_w <= w2rat) THEN
      e_threshold = collision%type_block%ethreshold
      g_mag = collision%g_mag

      ! Incoming normalised velocity vector
      v_inc = vector_normalisation(g)

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

      e_inc = 0.5_num * m1 * g_mag * g_mag
      e_scat = e_inc - e_threshold
      g_scat_m1 = SQRT(2._num * e_scat * m1)

      ! Post-collision momentum
      part1%part_p = v_scat * g_scat_m1
    END IF

    !IF (ran_w <= w1rat) THEN
    !  ! Post-collision momentum
    !  p_2 = u_cm*m2 - p_scat
    !  excited_id = collision%type_block%new_species_id
    !  IF (excited_id>0) THEN !Move particle to excited particle list
    !    new_part => part2
    !    ! Link part1 pointer to next colliding particle
    !    IF (.NOT.ASSOCIATED(p_list%head, TARGET=part2)) THEN
    !      part2 => part2%prev
    !    ELSE
    !      part2 => part2%next
    !    END IF
    !    new_part%part_p = p_2
    !    CALL remove_particle_from_partlist(p_list, new_part)
    !    CALL add_particle_to_partlist( &
    !      species_list(excited_id)%attached_list, new_part)
    !    NULLIFY(new_part)
    !  ELSE ! Just carry on with the same particle
    !    part2%part_p = p_2
    !  END IF
    !END IF

  END SUBROUTINE vahedi_excitation



  SUBROUTINE vahedi_ionisation(collision)

    ! Collision process: electron impact ionisation
    ! e + N -> 2e + I
    ! species1 = electron (e)
    ! species2 = neutral (N)

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision
    REAL(num), DIMENSION(3) :: g, v_inc, v_scat, v_inc_i
    REAL(num) :: costheta, sintheta, coschi, sinchi, cosphi, sinphi, phi
    REAL(num) :: e_inc, delta_e, e_scat, g_scat_m1, g_mag
    REAL(num) :: sinratio, m1, m2, part_pos, ran_e, ranw, w2rat
    INTEGER :: species_id, species1
    TYPE(particle), POINTER :: part1, part2, new_part

    ! Electron is species 1 and neutral => species 2
    CALL set_particle_properties(collision, m1 = m1, m2 = m2, &
      g = g, species1 = species1, w2_ratio = w2rat, &
#ifndef PER_SPECIES_WEIGHT
      w1 = w1, &
#endif
      part1 = part1, part2 = part2)

    ranw = random()
    IF (ranw <= w2rat) THEN

      part_pos = part1%part_pos
      g_mag = collision%g_mag

      e_inc = 0.5_num * m1 * g_mag * g_mag
      delta_e = collision%type_block%ethreshold
      e_scat = e_inc - delta_e
      g_scat_m1 = SQRT(2._num * e_scat * m1)
      ran_e = 0.5_num ! random()

      ! Incoming normalised velocity vector
      v_inc = vector_normalisation(g)
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
      part1%part_p = v_scat*g_scat_m1*SQRT(ran_e)

      ! Collision product from ionisation
      species_id = collision%type_block%new_species_id ! ions
      IF (species_id > 0) THEN
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
        CALL create_particle(new_part)
        new_part%part_p = v_scat*g_scat_m1*SQRT(1._num - ran_e)
        new_part%part_pos = part_pos
#ifndef PER_SPECIES_WEIGHT
        new_part%weight = w1 ! Electron's weight
#endif
        CALL add_particle_to_partlist(&
          species_list(species1)%attached_list, &
          new_part)
        NULLIFY(new_part)
      END IF


      ! Ion
      ! Create a new ion particle and assign neutral particle's momentum
      CALL create_particle(new_part)
      new_part%part_p = part2%part_p
      new_part%part_pos = part_pos
#ifndef PER_SPECIES_WEIGHT
      new_part%weight = w1 ! Electron's weight
#endif
      CALL add_particle_to_partlist(species_list(species_id)%attached_list, &
        new_part)
      NULLIFY(new_part)
    END IF

    ! Impact particle: neutral -> ion + electron
    !IF (ran_w <= w1rat) THEN
    !  ! The new electron
    !  ! Chi angle
    !  coschi = 1._num - 2._num * random()
    !  sinchi = SQRT(1._num - coschi*coschi)
    !  ! Phi angle
    !  phi = 2._num * pi * random()
    !  cosphi = COS(phi)
    !  sinphi = SIN(phi)
    !  ! Scattered normalised vector
    !  sinratio = sinchi / sintheta
    !  v_scat = v_inc * coschi + v_inc_i * sinratio * sinphi + &
    !    crossproduct(v_inc_i,v_inc) * sinratio * cosphi
    !  ! The new electron
    !  CALL create_particle(new_part)
    !  u_e2 = v_scat * g_scat * SQRT(1._num - ran_e)
    !  new_part%part_p = (u_cm + u_e2)*m1
    !  new_part%part_pos = part2%part_pos
    !#ifndef PER_SPECIES_WEIGHT
    !  new_part%weight = w2 ! Electron's weight
    !#endif
    !  CALL add_particle_to_partlist( &
    !    species_list(species1)%attached_list, new_part)
    !  NULLIFY(new_part)
    !  ! New ion: move neutral particle to ion species list
    !  new_part => part2
    !  new_part%part_p = collision%u_2 * m2
    !  IF (.NOT.ASSOCIATED(p_list%head, TARGET=part2)) THEN
    !    ! If particle is not the head of the list
    !    part2 => part2%prev
    !  ELSE
    !    ! If particle is the head of the list then one particle is skipped
    !    p_list%coll_counter = p_list%coll_counter + 1
    !    part2 => part2%next
    !  END IF
    !  CALL remove_particle_from_partlist(p_list, new_part)
    !  ion_id = collision%type_block%new_species_id ! ions
    !  CALL add_particle_to_partlist(species_list(ion_id)%attached_list, &
    !    new_part)
    !  NULLIFY(new_part)
    !END IF

  END SUBROUTINE vahedi_ionisation


#ifndef PER_SPECIES_WEIGHT
  SUBROUTINE vahedi_split_ionisation(collision)

    ! Collision process: electron impact ionisation
    ! e + N -> 2e + I
    ! species1 = electron (e)
    ! species2 = neutral (N)

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num), DIMENSION(3) :: u_cm, u_e1, u_e2, v_inc, v_inc_i, v_scat, g
    REAL(num) :: m1, m2, im1, mu, part2_pos
    REAL(num) :: w1, w2, w_diff, w_min, w1max_rat, w2max_rat
    REAL(num) :: ran_e, phi
    REAL(num) :: g_mag, e_threshold
    REAL(num) :: costheta, sintheta, cosphi, sinphi, coschi, sinchi
    REAL(num) :: sinratio, e_inc, e_scat, g_scat
    INTEGER :: ion_id, species1
    LOGICAL :: merge_electrons, electron_collides, neutral_collides
    TYPE(particle), POINTER :: new_part, part_next, part1, part2
    TYPE(particle_list) , POINTER :: p_list1, p_list2

    ! Link indexes so that 1 is electron and 2 is neutral
    CALL set_particle_properties(collision, species1 = species1, &
      m1 = m1, im1 = im1, m2 = m2, w1 = w1, w2 = w2, part1 = part1, &
      w1max_rat = w1max_rat, w2max_rat = w2max_rat, &
      g = g, part2 = part2, p_list1 = p_list1, p_list2 = p_list2)
    g_mag = collision%g_mag
    u_cm = cm_velocity(collision) 
    mu = collision%reducedm
    e_threshold = collision%type_block%ethreshold
    part2_pos = part2%part_pos

    ! Energy balance and post collision speed
    e_inc = 0.5_num * mu * g_mag * g_mag
    e_scat = e_inc - e_threshold
    g_scat = SQRT(2._num * e_scat * im1)
    ran_e = 0.5_num !random() ! Energy split ratio

    ! Incoming normalised velocity vector
    v_inc = vector_normalisation(g)
    v_inc_i = crossproduct(v_inc,(/1._num, 0._num, 0._num/))
    !Theta angle
    costheta = v_inc(1)
    sintheta = SQRT(1._num - costheta*costheta)

    ! Check particle's weight
    electron_collides = .FALSE.
    neutral_collides = .FALSE.
    IF (random() <= w1max_rat) electron_collides = .TRUE.
    IF (random() <= w2max_rat) neutral_collides = .TRUE.
    merge_electrons = .FALSE.
    IF (w1 >= w2) THEN
      ! Electron's weight is larger than neutral
      w_min = w2
      w_diff = w1 - w2
      !Update impact electron weight
      IF (electron_collides) part1%weight = w_diff
      !Remove neutral particle from list
      IF (neutral_collides) THEN
        new_part => part2%next ! This is just a buffer
        CALL remove_particle_from_partlist(p_list2, part2)
        CALL destroy_particle(part2)
        part2 => new_part
        NULLIFY(new_part)
      END IF
      IF (w_diff < w1) merge_electrons = .TRUE.
    ELSE
      ! Neutral's weight is larger
      w_min = w1
      w_diff = w2 - w1
      IF (neutral_collides) THEN
        part2%weight = w_diff
        IF (w_diff < w1) THEN
          ! If w_diff is small merge neutral with next particle
          new_part => part2 ! Part to be merged and removed
          part2 => part2%next
          CALL merge_particles(new_part, part2, p_list2)
          NULLIFY(new_part)
        END IF
      END IF
    END IF

    ! Impact electron
    IF (electron_collides) THEN
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
      u_e1 = v_scat * g_scat * SQRT(ran_e)
      part1%part_p = u_e1 * m1
    END IF

    ! Neutral (part2 pointer) is split into
    ! - New ion
    ! - New electron
    ! - Untouched neutral
    IF (neutral_collides) THEN
      ! New electron
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
      CALL create_particle(new_part)
      u_e2 = v_scat * g_scat * SQRT(1._num - ran_e)
      new_part%part_p = u_e2 * m1
      new_part%part_pos = part2_pos
      new_part%weight = w_min
      CALL add_particle_to_partlist( &
        species_list(species1)%attached_list, new_part)
      IF (merge_electrons) THEN
        part_next => part1%next
        CALL merge_particles(part1, new_part, p_list1)
        part1 => part_next
        NULLIFY(part_next)
      END IF
      NULLIFY(new_part)

      ! New ion
      ion_id = collision%type_block%new_species_id ! ions
      CALL create_particle(new_part)
      new_part%part_p = collision%u_2 * m2
      new_part%part_pos = part2_pos
      new_part%weight = w_min
      CALL add_particle_to_partlist(species_list(ion_id)%attached_list, &
        new_part)
      NULLIFY(new_part)
    END IF

    CALL link_particle_pointers(collision, part1, part2)
    NULLIFY(p_list1, p_list2)

  END SUBROUTINE vahedi_split_ionisation
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                             NANBU BACKGROUND                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE nanbu_elastic_scattering_bg(collision)

    ! Collision process: electron scattering
    ! e + N -> e + N

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num) :: reducedm, g_mag, m1
    REAL(num), DIMENSION(3) :: random_direction, p_scat, u_cm
#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: ran_w

    ran_w = random()
    IF (ran_w > collision%w1_ratio) RETURN
#endif

    reducedm = collision%reducedm
    g_mag = collision%g_mag
    u_cm = cm_velocity(collision) 
    m1 = collision%m1

    ! Post-collision momentum
    random_direction = random_unit_vector()
    p_scat = reducedm * g_mag * random_direction
    collision%part1%part_p = u_cm*m1 + p_scat 

  END SUBROUTINE nanbu_elastic_scattering_bg



  SUBROUTINE nanbu_ionisation_bg(collision)

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
#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: ran_w

    ran_w = random()
    IF (ran_w > collision%w1_ratio) RETURN
#endif

    m1 = collision%m1
    im1 = collision%im1
    m2 = collision%m2
    im2 = collision%im2
    u_cm = cm_velocity(collision) 
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

  END SUBROUTINE nanbu_ionisation_bg



  SUBROUTINE nanbu_excitation_bg(collision)

    ! Collision process: excitation
    ! e + N -> e + N*

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num) :: g_mag, e_threshold
    REAL(num) :: m1, m2, reducedm, ireducedm
    REAL(num), DIMENSION(3) :: random_direction, u_cm, p_scat
    INTEGER :: excited_id
    TYPE(particle), POINTER :: part1, new_part
#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: ran_w

    ran_w = random()
    IF (ran_w > collision%w1_ratio) RETURN
#endif

    m1 = collision%m1
    reducedm = collision%reducedm
    ireducedm = collision%ireducedm
    u_cm = cm_velocity(collision) 
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

  END SUBROUTINE nanbu_excitation_bg



  SUBROUTINE nanbu_charge_exchange_bg(collision)

    ! Collision process: charge_exchange
    ! I + N -> N + I
    ! Where the ions (I) are species1 and neutrals (N) are background

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num) :: mu, m1
    REAL(num), DIMENSION(3) :: g, p_scat, u_cm
#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: ranw

    ranw = random()
    IF (ranw > collision%w1_ratio) RETURN
#endif

    mu = collision%reducedm
    g = collision%g
    u_cm = cm_velocity(collision) 
    m1 = collision%m1

    ! Post-collision momentum
    p_scat = mu * g
    collision%part1%part_p = u_cm*m1 - p_scat 

  END SUBROUTINE nanbu_charge_exchange_bg


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                            VAHEDI BACKGROUND                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE vahedi_electron_elastic_scattering_bg(collision)

    ! Collision process: electron elastic scattering
    ! e + N -> e + N
    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num), DIMENSION(3) :: v_inc, v_scat, v_inc_i
    REAL(num) :: costheta, sintheta, coschi, sinchi, cosphi, sinphi, phi
    REAL(num) :: e_inc, delta_e, e_scat, g_scat_m1, g_mag
    REAL(num) :: sinratio, m1, im2
#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: ran_w

    ran_w = random()
    IF (ran_w > collision%w1_ratio) RETURN
#endif

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

    im2 = collision%im2
    m1 = collision%m1
    g_mag = collision%g_mag
    e_inc = 0.5_num * m1 * g_mag * g_mag
    delta_e = 2._num * m1 * im2 * (1._num - coschi)
    e_scat = e_inc * (1._num - delta_e)
    g_scat_m1 = SQRT(2._num * e_scat * m1)

    ! Post-collision momentum
    collision%part1%part_p = v_scat * g_scat_m1

  END SUBROUTINE vahedi_electron_elastic_scattering_bg



  SUBROUTINE vahedi_ion_elastic_scattering_bg(collision)

    ! Collision process: ion elastic scattering
    ! e + N -> e + N
    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num), DIMENSION(3) :: v_inc, v_scat, v_inc_i, u_cm
    REAL(num) :: costheta, sintheta, coschi, sinchi, cosphi, sinphi, phi
    REAL(num) :: e_inc, e_scat, g_scat_m1, g_mag, ucm_mag2
    REAL(num) :: sinratio, m1, m12
#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: ran_w

    ran_w = random()
    IF (ran_w > collision%w1_ratio) RETURN
#endif

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

    m1 = collision%m1
    m12 = collision%m12
    u_cm = cm_velocity(collision) 
    g_mag = collision%g_mag
    ucm_mag2 = DOT_PRODUCT(u_cm, u_cm)
    e_inc = 0.5_num * m1 * g_mag * g_mag
    e_scat = e_inc * coschi*coschi + 0.5_num * m12 * ucm_mag2
    g_scat_m1 = SQRT(2._num * e_scat * m1)

    ! Post-collision momentum
    collision%part1%part_p = v_scat * g_scat_m1

  END SUBROUTINE vahedi_ion_elastic_scattering_bg



  SUBROUTINE vahedi_excitation_bg(collision)

    ! Collision process: electron impact excitation
    ! e + N -> e + N
    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision

    REAL(num), DIMENSION(3) :: v_inc, v_scat, v_inc_i
    REAL(num) :: costheta, sintheta, coschi, sinchi, cosphi, sinphi, phi
    REAL(num) :: e_inc, delta_e, e_scat, g_scat_m1, g_mag
    REAL(num) :: sinratio, m1
#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: ran_w

    ran_w = random()
    IF (ran_w > collision%w1_ratio) RETURN
#endif

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
    g_scat_m1 = SQRT(2._num * e_scat * m1)

    ! Post-collision momentum
    collision%part1%part_p = v_scat * g_scat_m1

  END SUBROUTINE vahedi_excitation_bg



  SUBROUTINE vahedi_ionisation_bg(collision)

    ! Collision process: electron impact ionisation
    ! e + N -> e + N
    TYPE(particle), POINTER :: new_part

    TYPE(current_collision_block), POINTER, INTENT(INOUT) :: collision
    REAL(num), DIMENSION(3) :: v_inc, v_scat, v_inc_i
    REAL(num) :: costheta, sintheta, coschi, sinchi, cosphi, sinphi, phi
    REAL(num) :: e_inc, delta_e, e_scat, g_scat_m1, g_mag
    REAL(num) :: sinratio, m1, part_pos, ran_e
    INTEGER :: species_id
#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: ran_w

    ran_w = random()
    IF (ran_w > collision%w1_ratio) RETURN
#endif

    part_pos = collision%part1%part_pos
    m1 = collision%m1
    g_mag = collision%g_mag

    e_inc = 0.5_num * m1 * g_mag * g_mag
    delta_e = collision%type_block%ethreshold
    e_scat = e_inc - delta_e
    g_scat_m1 = SQRT(2._num * e_scat * m1)
    ran_e = 0.5_num ! random()

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
    collision%part1%part_p = v_scat*g_scat_m1*SQRT(ran_e)

    ! Collision product from ionisation
    species_id = collision%type_block%new_species_id ! ions
    IF (species_id > 0) THEN
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
      CALL create_particle(new_part)
      new_part%part_p = v_scat*g_scat_m1*SQRT(1._num - ran_e)
      new_part%part_pos = part_pos
#ifndef PER_SPECIES_WEIGHT
      new_part%weight = collision%w1 ! Electron's weight
#endif
      CALL add_particle_to_partlist(&
        species_list(collision%species1)%attached_list, &
        new_part)
      NULLIFY(new_part)


      ! Ion
      ! Create a new ion particle and assign neutral particle's momentum
      CALL create_particle(new_part)
      new_part%part_p = collision%u_2 * collision%m2
      new_part%part_pos = part_pos
#ifndef PER_SPECIES_WEIGHT
      new_part%weight = collision%w1 ! Electron's weight
#endif
      CALL add_particle_to_partlist(species_list(species_id)%attached_list, &
        new_part)
      NULLIFY(new_part)
    END IF

  END SUBROUTINE vahedi_ionisation_bg

  

  FUNCTION cm_velocity(collision)

    REAL(num), DIMENSION(3) :: cm_velocity
    TYPE(current_collision_block) :: collision

    cm_velocity = (collision%part1%part_p+collision%m2*collision%u_2) &
      * collision%im12
  END FUNCTION cm_velocity

#endif
END MODULE nc_subroutines
