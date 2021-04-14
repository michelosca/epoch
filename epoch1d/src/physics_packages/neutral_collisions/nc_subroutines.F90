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
  
  PUBLIC :: nanbu_elastic_scattering, nanbu_ionisation, nanbu_excitation
  PUBLIC :: nanbu_charge_exchange, nanbu_excitation_bg
  PUBLIC :: nanbu_elastic_scattering_bg, nanbu_ionisation_bg
  PUBLIC :: nanbu_charge_exchange_bg
  PUBLIC :: vahedi_electron_elastic_scattering_bg, vahedi_excitation_bg
  PUBLIC :: vahedi_ion_elastic_scattering_bg, vahedi_ionisation_bg
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

  END SUBROUTINE nanbu_elastic_scattering



  SUBROUTINE nanbu_ionisation(collision)

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

  END SUBROUTINE nanbu_ionisation



  SUBROUTINE nanbu_excitation(collision)

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

  END SUBROUTINE nanbu_excitation



  SUBROUTINE nanbu_charge_exchange(collision)

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
  END SUBROUTINE nanbu_charge_exchange


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                 BACKGROUND                                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE nanbu_elastic_scattering_bg(collision)

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

  END SUBROUTINE nanbu_excitation_bg



  SUBROUTINE nanbu_charge_exchange_bg(collision)

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

  END SUBROUTINE nanbu_charge_exchange_bg



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
    REAL(num) :: sinratio, mu

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
