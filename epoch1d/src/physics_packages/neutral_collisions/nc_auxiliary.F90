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
#ifdef NEUTRAL_COLLISIONS
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

    INTEGER :: species1, species2
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
    
      WRITE(*,*)
      WRITE(*,*) '*** WARNING ***'
      WRITE(*,*) 'MAX(g*sigma) of species ', TRIM(ADJUSTL(ispecies_str)), &
        ' and ', TRIM(ADJUSTL(jspecies_str)), ' is lower than expected.'
      WRITE(*,*) 'Position [m]:         ', TRIM(ADJUSTL(pos_str))
      WRITE(*,*) 'Time [s]:             ', TRIM(ADJUSTL(time_str))
      WRITE(*,*) 'g*sigma [m^3/s]:      ', TRIM(ADJUSTL(gsigma_str))
      WRITE(*,*) 'MAX(g*sigma) [m^3/s]: ', TRIM(ADJUSTL(gsigma_max_str))
      WRITE(*,*) 'g [m/s]:              ', TRIM(ADJUSTL(g_str))
      WRITE(*,*) 'Boyd factor []:       ', TRIM(ADJUSTL(boyd_factor_str))
      WRITE(*,*)

987 FORMAT (ES14.6)

  END SUBROUTINE warning_pmax



  SUBROUTINE test_neutral_collision_setup


    LOGICAL :: is_background
    INTEGER :: ncerr
    INTEGER :: ispecies, jspecies, nc_type, ibg
    REAL(num) :: m1, m2, m_new, m_source
    CHARACTER(len=string_length) :: iname, jname, type_name
    TYPE(collision_type_block), POINTER :: coll_type
    TYPE(neutrals_block), POINTER :: coll_block

    ncerr = 0

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

          IF (coll_type%wnanbu .OR. coll_type%wvahedi) THEN
            IF (coll_type%id == c_nc_ionisation) THEN

              ! Nanbu-ionisation requires species_source_id
              IF ((coll_type%source_species_id <= 0 .OR. &
                  coll_type%source_species_id > n_species_bg) .AND. &
                  .NOT.is_background) ncerr = 2

              ! Nanbu-ionisation requires species_target_id
              IF ((coll_type%new_species_id <= 0 .OR. &
                 coll_type%new_species_id > n_species) .AND. &
                 .NOT.is_background) ncerr = 3

              ! Nanbu-ionisation requires same mass for source and new species 
              IF ( ncerr == 0 ) THEN
                IF (is_background) THEN
                  m_source = background_list(ibg)%mass
                ELSE
                  m_source = species_list(coll_type%source_species_id)%mass
                END IF
                IF (coll_type%new_species_id > 0) THEN
                  m_new = species_list(coll_type%new_species_id)%mass
                  IF (ABS(m_source - m_new) > TINY(0._num)) ncerr = 11
                END IF
              END IF

            ELSE IF (coll_type%id == c_nc_excitation) THEN

              IF (coll_type%source_species_id > 0 .AND. &
                  coll_type%source_species_id <= n_species_bg &
                  .AND. .NOT.is_background) THEN

                ! Nanbu-excitation has new_species_id when source_species_id
                IF (coll_type%new_species_id < 0 .OR. &
                    coll_type%new_species_id > n_species) THEN
                  ncerr = 4
                  CYCLE
                END IF

                ! Nanbu-excitation requires same mass for source and new species
                m_source = species_list(coll_type%source_species_id)%mass
                m_new = species_list(coll_type%new_species_id)%mass
                IF (ABS(m_source - m_new) > TINY(0._num)) ncerr = 10
              END IF

            END IF
          ELSE IF (coll_type%wnanbusplit .OR. coll_type%wvahedisplit) THEN
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

            END IF

          END IF
          IF (ncerr /= 0) EXIT
        END DO
        IF (ncerr /= 0) EXIT
      END DO
      IF (ncerr /= 0) EXIT
    END DO

    IF (ncerr /= 0) THEN
      IF (rank == 0) THEN
        WRITE(*,*) ''
        WRITE(*,*) '*** ERROR ***'
        WRITE(*,*) 'Collision pair: ', TRIM(ADJUSTL(iname)),' - ', &
          TRIM(ADJUSTL(jname))
        IF (ncerr == 2) THEN
          WRITE(*,*) 'Ionisation requires requires a valid ', &
            '"source_species"'
        ELSE IF (ncerr ==3) THEN
          WRITE(*,*) 'Ionisation requires "new_species"'
        ELSE IF (ncerr ==4) THEN
          WRITE(*,*) 'Excitation requires a valid ', &
            '"new_species" when "source_species" is specified'
        ELSE IF (ncerr ==5) THEN
          WRITE(*,*) 'Split methods requires a valid "source_species"'
        ELSE IF (ncerr ==6) THEN
          WRITE(*,*) 'Split methods requires a valid "new_species"'
        ELSE IF (ncerr ==8) THEN
          WRITE(*,*) 'Collision type is not recognised'
        ELSE IF (ncerr ==10) THEN
          WRITE(*,*) 'Excitation requires that source and new species', &
            ' have the same mass'
        ELSE IF (ncerr ==11) THEN
          WRITE(*,*) 'Ionisation requires that source and new species', &
            ' have the same mass'
        ELSE IF (ncerr ==13) THEN
          WRITE(*,*) 'Split methods does not support background collisions'
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

    IF (coll_type_block%wnanbu) THEN
      WRITE(*,666) 'Nanbu method'
    ELSE IF (coll_type_block%wvahedi) THEN
      WRITE(*,666) 'Vahedi method'
    ELSE IF (coll_type_block%wnanbusplit) THEN
      WRITE(*,666) 'Nanbu(split) method'
    ELSE IF (coll_type_block%wvahedisplit) THEN
      WRITE(*,666) 'Vahedi(split) method'
    END IF
    WRITE(*,555) 'Output name:', TRIM(ADJUSTL(coll_type_block%io_name))

999 FORMAT(2X, A, 1X, A)
888 FORMAT(6X, A, 1X, G10.4)
777 FORMAT(6X, A, 1X, I2)
666 FORMAT(6X, A)
555 FORMAT(6X, A, 1X, A)
333 FORMAT(6X, A, 1X, I4)

  END SUBROUTINE print_collision_type

#endif
END MODULE nc_auxiliary
