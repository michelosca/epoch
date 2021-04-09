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

MODULE nc_auxiliary

  USE shared_data
  USE utilities

  IMPLICIT NONE

CONTAINS

  FUNCTION interpolation(x1, x2, y1, y2, x_in)

    REAL(num) :: interpolation
    REAL(num) :: x1, x2, y1, y2, x_in

    interpolation = y1 + (x_in - x1) * (y2 -y1) / (x2 - x1)

  END FUNCTION interpolation



  SUBROUTINE crosssection_out_of_range(gsigma, g, coll_type_block)

    REAL(num), INTENT(IN) :: g
    REAL(num), INTENT(INOUT) :: gsigma
    TYPE(collision_type_block), POINTER, INTENT(IN) :: coll_type_block
    REAL(num) :: g_max, g_min

    IF ( gsigma < 0._num ) THEN

      g_min = coll_type_block%energy(1)
      g_max = coll_type_block%energy(coll_type_block%table_len)

      IF (g <= g_min) THEN
        gsigma = coll_type_block%cross_section(1)
      ELSE IF (g >= g_max) THEN
        gsigma = coll_type_block%cross_section(coll_type_block%table_len)
      ELSE
        WRITE(*,*) 'ERROR: Cross-section out of range but energy is not!'
      END IF

    END IF

  END SUBROUTINE crosssection_out_of_range



  SUBROUTINE warning_pmax(collision, gsigma_total, g)

    TYPE(current_collision_block), POINTER ,INTENT(IN) :: collision
    REAL(num), INTENT(IN) :: gsigma_total, g

    INTEGER :: species1, species2, iu, io
    CHARACTER(14) :: time_str, pos_str, gsigma_str, gsigma_max_str
    CHARACTER(14) :: g_str, boyd_factor_str
    CHARACTER(len=string_length) :: ispecies_str, jspecies_str
    TYPE(neutrals_block), POINTER:: collision_block

    species1 = collision%species1
    species2 = collision%species2

    collision_block => species_list(species1)%neutrals(species2)

    WRITE(time_str, 987) time
    WRITE(pos_str, 987) collision%part1%part_pos
    WRITE(gsigma_str, 987) gsigma_total
    WRITE(gsigma_max_str, 987) collision_block%gsigma_max_total
    WRITE(g_str, 987) g
    WRITE(boyd_factor_str, 987) collision%prob_factor
    ispecies_str = species_list(species1)%name
    IF (species2 > n_species) THEN
      jspecies_str = background_list(species2 - n_species)%name
    ELSE
      jspecies_str = species_list(species2)%name
    END IF
    
    DO iu = 1, nio_units ! Print to stdout and to file
      io = io_units(iu)
      WRITE(io,*)
      WRITE(io,*) '*** WARNING ***'
      WRITE(io,*) 'MAX(g*sigma) of species ', TRIM(ADJUSTL(ispecies_str)), &
        ' and ', TRIM(ADJUSTL(jspecies_str)), ' is lower than expected.'
      WRITE(io,*) 'Position [m]:         ', TRIM(ADJUSTL(pos_str))
      WRITE(io,*) 'Time [s]:             ', TRIM(ADJUSTL(time_str))
      WRITE(io,*) 'g*sigma [m^3/s]:      ', TRIM(ADJUSTL(gsigma_str))
      WRITE(io,*) 'MAX(g*sigma) [m^3/s]: ', TRIM(ADJUSTL(gsigma_max_str))
      WRITE(io,*) 'g [m/s]:              ', TRIM(ADJUSTL(g_str))
      WRITE(io,*) 'Boyd factor []:       ', TRIM(ADJUSTL(boyd_factor_str))
      WRITE(io,*)
    END DO

987 FORMAT (ES14.6)

  END SUBROUTINE warning_pmax



  SUBROUTINE warning_coll_pairs(collision, n_pairs)

    TYPE(current_collision_block), POINTER ,INTENT(IN) :: collision
    INTEGER(8), INTENT(IN) :: n_pairs

    INTEGER :: species1, species2, iu, io
    INTEGER(8) :: n_part1, n_part2, coll_counter1, coll_counter2
    REAL(num) :: time, position
    CHARACTER(14) :: time_str, pos_str, n_ipart_str, n_jpart_str
    CHARACTER(14) :: n_pairs_str, coll_icounter_str, coll_jcounter_str
    CHARACTER(LEN=string_length) :: ispecies_str, jspecies_str

    species1 = collision%species1
    species2 = collision%species2
    position = collision%p_list1%head%part_pos
    n_part1 = collision%p_list1%count
    IF (collision%collision_block%is_background) THEN
      n_part2 = -1
    ELSE
      n_part2 = collision%p_list2%count
    END IF

    coll_counter1 = collision%p_list1%coll_counter
    IF (collision%collision_block%is_background) THEN
      coll_counter2 = -1
    ELSE
      coll_counter2 = collision%p_list2%coll_counter
    END IF

    WRITE(time_str, 987) time
    WRITE(pos_str, 987) position
    WRITE(n_ipart_str, 765) n_part1
    WRITE(n_jpart_str, 765) n_part2
    WRITE(n_pairs_str, 765) n_pairs
    WRITE(coll_icounter_str, 765) coll_counter1
    IF (coll_counter2 == -1) THEN
      WRITE(coll_jcounter_str,*) 'N/A'
    ELSE
      WRITE(coll_jcounter_str, 765) coll_counter2
    END IF
    ispecies_str = species_list(species1)%name
    IF (species2 > n_species) THEN
      jspecies_str = background_list(species2 - n_species)%name
    ELSE
      jspecies_str = species_list(species2)%name
    END IF

    DO iu = 1, nio_units ! Print to stdout and to file
      io = io_units(iu)
      WRITE(io,*)
      WRITE(io,*) '*** WARNING ***'
      WRITE(io,*) 'Too many collision pairs between species ', &
        TRIM(ADJUSTL(ispecies_str)), ' and ', TRIM(ADJUSTL(jspecies_str)),'.'
      WRITE(io,*) 'Position [m]:    ', TRIM(ADJUSTL(pos_str))
      WRITE(io,*) 'Time [s]:        ', TRIM(ADJUSTL(time_str))
      WRITE(io,*) 'Number Sup-pcls: ', TRIM(ADJUSTL(n_ipart_str)), &
        ' and ', TRIM(ADJUSTL(n_jpart_str))
      WRITE(io,*) 'Coll. counter:   ', TRIM(ADJUSTL(coll_icounter_str)), &
        ' and ', TRIM(ADJUSTL(coll_jcounter_str))
      WRITE(io,*) 'Collision pairs: ', TRIM(ADJUSTL(n_pairs_str))
      WRITE(io,*)
    END DO

987 FORMAT (ES14.6)
765 FORMAT (I14)

  END SUBROUTINE warning_coll_pairs



  SUBROUTINE test_neutral_collision_setup


    LOGICAL :: boyd_collision, boyd_same_species, is_background
    INTEGER :: ncerr
    INTEGER :: ispecies, jspecies, nc_type, ibg
    REAL(num) :: m1, m2, m_new, m_source
    CHARACTER(len=string_length) :: iname, jname, type_name
    TYPE(collision_type_block), POINTER :: coll_type
    TYPE(neutrals_block), POINTER :: coll_block

    ncerr = 0

    boyd_collision = .FALSE.
    boyd_same_species = .FALSE.

    DO ispecies = 1, n_species
      iname = species_list(ispecies)%name
      
      DO jspecies = 1, n_species_bg

        IF (.NOT.neutral_coll(ispecies,jspecies)) CYCLE
        coll_block => species_list(ispecies)%neutrals(jspecies)

        is_background = coll_block%is_background
        m1 = species_list(ispecies)%mass
        IF (is_background) THEN
          ibg = jspecies - n_species
          jname = coll_block%background%name
          m2 = coll_block%background%mass
        ELSE
          jname = species_list(jspecies)%name
          m2 = species_list(jspecies)%mass
        END IF

        DO nc_type = 1, coll_block%ncolltypes

          coll_type => coll_block%collision_set(nc_type)
          type_name = coll_type%name

          ! In case no collision types are found
          IF (coll_type%id < 0) ncerr = 8

          IF (coll_type%wbird) THEN
            IF (coll_type%id == c_nc_ionisation) THEN

              ! Bird-ionisation requires species_source_id
              IF ((coll_type%source_species_id <= 0 .OR. &
                  coll_type%source_species_id > n_species_bg) .AND. &
                  .NOT.is_background) ncerr = 2

              ! Bird-ionisation requires species_target_id
              IF (coll_type%new_species_id <= 0 .OR. &
                 coll_type%new_species_id > n_species) ncerr = 3

              ! Bird-ionisation requires same mass for source and new species 
              IF ( ncerr == 0 ) THEN
                IF (is_background) THEN
                  m_source = background_list(ibg)%mass
                ELSE
                  m_source = species_list(coll_type%source_species_id)%mass
                END IF
                m_new = species_list(coll_type%new_species_id)%mass
                IF (ABS(m_source - m_new) > TINY(0._num)) ncerr = 11
              END IF

            ELSE IF (coll_type%id == c_nc_excitation) THEN

              IF (coll_type%source_species_id > 0 .AND. &
                  coll_type%source_species_id <= n_species_bg &
                  .AND. .NOT.is_background) THEN

                ! Bird-excitation has new_species_id when source_species_id
                IF (coll_type%new_species_id < 0 .OR. &
                    coll_type%new_species_id > n_species) THEN
                  ncerr = 4
                  CYCLE
                END IF

                ! Bird-excitation requires same mass for source and new species
                m_source = species_list(coll_type%source_species_id)%mass
                m_new = species_list(coll_type%new_species_id)%mass
                IF (ABS(m_source - m_new) > TINY(0._num)) ncerr = 10
              END IF

            ELSE IF (coll_type%id == c_nc_charge_exchange) THEN
              IF (ABS(m1 - m2) > TINY(0._num)) ncerr = 9
            END IF

          ELSE IF (coll_type%wboyd) THEN
          
            !No background collision available for split method
            IF (is_background) ncerr = 12

            ! Boyd's method requires a target species for energy balance
            boyd_collision = .TRUE.
            IF (energy_correction_species <= 0) ncerr = 1

            ! An intra-collision type must be defined for energy balance
            IF (ispecies==jspecies) boyd_same_species = .TRUE.

          ELSE IF (coll_type%wsplit) THEN
            !No background collision available for split method
            IF (is_background) ncerr = 13

            IF (coll_type%id == c_nc_ionisation) THEN

              ! Valid source_species_id
              IF (coll_type%source_species_id <= 0 .OR. &
                  coll_type%source_species_id > n_species .AND. &
                  .NOT.is_background) ncerr = 5

              ! Valir new_species_id
              IF (coll_type%new_species_id <= 0 .OR. &
                  coll_type%new_species_id > n_species) ncerr = 6

              ! Ionisation requires same mass for source and new species 
              m_source = species_list(coll_type%source_species_id)%mass
              m_new = species_list(coll_type%new_species_id)%mass
              IF (ABS(m_source - m_new) > TINY(0._num)) ncerr = 11

            ELSE IF (coll_type%id == c_nc_charge_exchange) THEN
              ! Must point which particle has the largest weight
              IF (coll_type%source_species_id <= 0 .OR. &
                  coll_type%source_species_id > n_species) ncerr = 5

              ! Collision species must have the same mass
              m1 = species_list(ispecies)%mass
              m2 = species_list(jspecies)%mass
              IF (ABS(m1 - m2) > TINY(0._num)) ncerr = 9

            END IF
          END IF
          IF (ncerr /= 0) EXIT
        END DO
        IF (ncerr /= 0) EXIT
      END DO
      IF (ncerr /= 0) EXIT
    END DO

    IF (boyd_collision .AND. .NOT.boyd_same_species) ncerr = 7

    IF (ncerr /= 0) THEN
      IF (rank == 0) THEN
        WRITE(*,*) ''
        WRITE(*,*) '*** ERROR ***'
        WRITE(*,*) 'Collision pair: ', TRIM(ADJUSTL(iname)),' - ', &
          TRIM(ADJUSTL(jname))
        IF (ncerr == 1) THEN
          WRITE(*,*) 'Boyd requires sepecific "energy_correction_species"'
        ELSE IF (ncerr == 2) THEN
          WRITE(*,*) 'Bird-ionisation requires requires a valid ', &
            '"source_species"'
        ELSE IF (ncerr ==3) THEN
          WRITE(*,*) 'Bird-ionisation requires "new_species"'
        ELSE IF (ncerr ==4) THEN
          WRITE(*,*) 'Bird-excitation requires a valid ', &
            '"new_species" when "source_species" is specified'
        ELSE IF (ncerr ==5) THEN
          WRITE(*,*) 'Split methods requires a valid ', &
            '"source_species"'
        ELSE IF (ncerr ==6) THEN
          WRITE(*,*) 'Split methods requires a valid ', &
            '"new_species"'
        ELSE IF (ncerr ==7) THEN
          WRITE(*,*) 'Boyd intra-species collision type must be defined ', &
            'for energy balance'
        ELSE IF (ncerr ==8) THEN
          WRITE(*,*) 'Collision type is not recognised'
        ELSE IF (ncerr ==9) THEN
          WRITE(*,*) 'CX collisions require same mass species'
        ELSE IF (ncerr ==10) THEN
          WRITE(*,*) 'Excitation requires that source and new species', &
            ' have the same mass'
        ELSE IF (ncerr ==11) THEN
          WRITE(*,*) 'Ionisation requires that source and new species', &
            ' have the same mass'
        ELSE IF (ncerr ==12) THEN
          WRITE(*,*) 'Boyd method does not support background collisions'
        ELSE IF (ncerr ==13) THEN
          WRITE(*,*) 'Split method does not support background collisions'
        END IF
      CALL print_collision_type(ispecies, jspecies, nc_type)
      END IF
      CALL MPI_BARRIER(comm, ncerr)
      CALL abort_code(c_err_bad_setup)
    END IF

  END SUBROUTINE test_neutral_collision_setup



  SUBROUTINE print_collisions

    INTEGER :: i, j, ctype
    CHARACTER(LEN=string_length) :: iname, jname
    TYPE(neutrals_block), POINTER :: coll_block
    TYPE(background_block), POINTER :: background

    IF (rank == 0 .AND. use_collisions) THEN
      WRITE(*,*) ' '
      WRITE(*,*) '**********************NEUTRAL COLLISIONS*********************'
      DO i = 1, n_species
        iname = species_list(i)%name
        DO j = i, n_species_bg
          IF (.NOT.neutral_coll(i, j)) CYCLE

          coll_block => species_list(i)%neutrals(j)
          IF (coll_block%is_background) THEN
            jname = coll_block%background%name
          ELSE
            jname = species_list(j)%name
          END IF

          WRITE(*,*)
          WRITE(*,333) 'Species:', TRIM(iname), 'and', TRIM(jname)
          WRITE(*,888) 'Max. g-sigma [m^3/s]:', coll_block%gsigma_max_total
          WRITE(*,666) 'User max. g-sigma:', coll_block%user_gsigma_max
          WRITE(*,777) 'Collision types:', coll_block%ncolltypes
          WRITE(*,888) 'Max. Weight:', coll_block%max_weight

          IF (coll_block%is_background) THEN
            background => coll_block%background
            WRITE(*,900) 'Background:', TRIM(ADJUSTL(background%name))
            WRITE(*,777) 'ID:', background%id
            WRITE(*,888) 'Mass:', background%mass
            WRITE(*,888) 'Density:', background%dens
            WRITE(*,888) 'Temperature:', background%temp
            WRITE(*,500) 'Mean velocity:', background%mean_velocity
            WRITE(*,888) 'Start time:', background%t_start
            WRITE(*,888) 'End time:', background%t_end
          END IF

          DO ctype = 1, coll_block%ncolltypes
            CALL print_collision_type(i, j, ctype)
          END DO
          WRITE(*,*)
        END DO
      END DO
    END IF

900 FORMAT(A, 1X, A)
500 FORMAT(A, 3(1X, G10.4))
333 FORMAT(4(A, 1X))
888 FORMAT(A, 1X, G10.4)
777 FORMAT(A, 1X, I2)
666 FORMAT(A, 1X, L1)

  END SUBROUTINE print_collisions



  SUBROUTINE print_collision_type(ispecies, jspecies, ctype)

    INTEGER, INTENT(IN) :: ispecies, jspecies, ctype
    INTEGER :: species_id
    TYPE(neutrals_block), POINTER :: coll_block
    TYPE(collision_type_block), POINTER :: coll_type_block

    coll_block => species_list(ispecies)%neutrals(jspecies)
    coll_type_block => coll_block%collision_set(ctype)

    WRITE(*,999) 'Collision type:', TRIM(coll_type_block%name)

    WRITE(*,888) 'Threshold energy:', &
      coll_type_block%ethreshold / coll_type_block%energy_units

    WRITE(*,888) 'Energy units:', coll_type_block%energy_units

    WRITE(*,888) 'Sigma units:', coll_type_block%cross_section_units

    WRITE(*,777) 'Id:', coll_type_block%id

    WRITE(*,888) 'Collision user factor:', coll_type_block%user_factor

    species_id = coll_type_block%source_species_id
    IF ( species_id > 0 .AND..NOT.coll_block%is_background) THEN
      WRITE(*,555) 'Source species:', &
       TRIM(ADJUSTL(species_list(species_id)%name))
    END IF

    species_id = coll_type_block%new_species_id
    IF ( species_id > 0 ) THEN
      WRITE(*,555) 'Target species:', &
        TRIM(ADJUSTL(species_list(species_id)%name))
    END IF

    WRITE(*,333) 'Table_len:', coll_type_block%table_len

!     WRITE(*,555) 'Energy vs. Cross-section'
!     DO species_id = 1, coll_type_block%table_len
!       WRITE(*,*) &
!       coll_type_block%energy(species_id) / coll_type_block%energy_units, &
!       coll_type_block%cross_section(species_id) / &
!       coll_type_block%cross_section_units
!     END DO

    IF (coll_type_block%wboyd) THEN
      WRITE(*,666) 'Boyd method'
      WRITE(*,444) 'Energy correction: ', coll_type_block%energy_correction
    ELSE IF (coll_type_block%wbird) THEN
      WRITE(*,666) 'Bird method'
    ELSE IF (coll_type_block%wsplit) THEN
      WRITE(*,666) 'Split method'
    ELSE IF (coll_type_block%wvahedi) THEN
      WRITE(*,666) 'Vahedi method'
    END IF
    WRITE(*,555) 'Output name:', TRIM(ADJUSTL(coll_type_block%io_name))

999 FORMAT(2X, A, 1X, A)
888 FORMAT(6X, A, 1X, G10.4)
777 FORMAT(6X, A, 1X, I2)
666 FORMAT(6X, A)
555 FORMAT(6X, A, 1X, A)
444 FORMAT(6X, A, L2)
333 FORMAT(6X, A, 1X, I4)

  END SUBROUTINE print_collision_type

END MODULE nc_auxiliary
