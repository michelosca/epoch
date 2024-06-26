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

MODULE nc_setup
#ifdef NEUTRAL_COLLISIONS

  USE strings_advanced
  USE shared_data
  USE background
  USE nc_subroutines
  USE nc_auxiliary

  IMPLICIT NONE

  LOGICAL :: allocated_energy_loss

CONTAINS
  
  SUBROUTINE neutral_collisions_start
  
    ALLOCATE(neutral_coll(n_species_bg, n_species_bg))
    ALLOCATE(coulomb_coll(n_species_bg, n_species_bg))
    neutral_coll = .FALSE.
    coulomb_coll = .FALSE.
    
  END SUBROUTINE neutral_collisions_start
  
  
  
  SUBROUTINE load_neutral_collisions

    INTEGER :: ispecies, jspecies, n_coll_type, lines, l, species
    INTEGER :: errcode, iread, file_id, header_lines
    LOGICAL :: ij_exists, ji_exists, header, is_background
    
    CHARACTER(string_length) :: element, value
    CHARACTER(string_length) :: ispecies_name, jspecies_name
    CHARACTER(string_length) :: ij_filename, ij_type_k_filename
    CHARACTER(string_length) :: ji_filename, ji_type_k_filename
    CHARACTER(string_length) :: ij_full_path_file, ji_full_path_file
    CHARACTER(string_length) :: full_path_file
    CHARACTER(4) :: type_str

    TYPE(neutrals_block), POINTER :: coll_block
    TYPE(collision_type_block), POINTER :: coll_type
    TYPE(background_block), POINTER :: background

    file_id = 140892
    errcode = c_err_none

    ! This is used for output purposes
    ALLOCATE(total_collision_types(n_species))
    total_collision_types = 0

    ! Look for cross section files
    DO ispecies = 1, n_species
      ispecies_name = species_list(ispecies)%name

      DO jspecies = ispecies, n_species_bg

        IF ( .NOT.neutral_coll(ispecies,jspecies) ) CYCLE

        IF (jspecies > n_species) THEN
          background => background_list(jspecies - n_species)
          jspecies_name = background%name
          is_background = .TRUE.
        ELSE
          jspecies_name = species_list(jspecies)%name
          is_background = .FALSE.
        END IF
        ij_filename = TRIM(ispecies_name) // '_' // TRIM(jspecies_name)
        ji_filename = TRIM(jspecies_name) // '_' // TRIM(ispecies_name)

        ! Read each table-file twice, the first time it identifies how many
        ! collision types there are per species pair
        DO iread = 1, 2
          !After reading the number of collision types, allocate memory
          IF (iread == 2) THEN
            IF (n_coll_type > 0) THEN
              IF (.NOT.ASSOCIATED(species_list(ispecies)%neutrals)) THEN
                ALLOCATE(species_list(ispecies)%neutrals(n_species_bg))
                  DO species = 1, n_species_bg
                    coll_block => species_list(ispecies)%neutrals(species)
                    CALL init_dummy_collision_block(coll_block)
                  END DO
              END IF
              coll_block => species_list(ispecies)%neutrals(jspecies)
              CALL init_collision_block(coll_block, n_coll_type, &
                ispecies, jspecies)
              total_collision_types(ispecies) = &
                total_collision_types(ispecies) + n_coll_type
              IF (is_background) coll_block%background => background
            ELSE ! No table is found
              IF (rank == 0) THEN
                WRITE(*,*)
                WRITE(*,*) '*** ERROR ***'
                WRITE(*,*) 'No energy-cross-section table is found for', &
                  ' species pair ', TRIM(ADJUSTL(ij_filename))
                WRITE(*,*) 'Please add table and rerun code'
                WRITE(*,*) ''
              END IF
              CALL MPI_BARRIER(comm, errcode)
              CALL abort_code(c_err_missing_elements)
            END IF
          END IF

          n_coll_type = 0
          
          ! Loop over as many energy-cross-section tables there are
          DO
            ! Collision type file
            WRITE(type_str, '(i2.2)') n_coll_type + 1

            ! Table file name
            ij_type_k_filename = TRIM(ij_filename) // '_' &
                // TRIM(type_str) //'.table'
            ji_type_k_filename = TRIM(ji_filename) // '_' &
                // TRIM(type_str) //'.table'

            ! Table file name + full path
            ij_full_path_file = TRIM(cross_section_table_location) // '/' &
                // TRIM(ij_type_k_filename)
            ji_full_path_file = TRIM(cross_section_table_location) // '/' &
                // TRIM(ji_type_k_filename)

            ! Checks whether collision type k in i-j species is found
            INQUIRE(FILE= TRIM(ij_full_path_file), EXIST=ij_exists)
            INQUIRE(FILE= TRIM(ji_full_path_file), EXIST=ji_exists)

            ! Read file
            IF (ij_exists .OR. ji_exists) THEN
              ! Update de collision-type counter
              n_coll_type = n_coll_type + 1

              ! This first read only counts the number of collision types
              IF ( iread == 1 ) CYCLE
              ! Only one of the two options should exits
              IF (ji_exists) full_path_file = ji_full_path_file
              IF (ij_exists) full_path_file = ij_full_path_file

              ! Point to the collision type
              coll_type => coll_block%collision_set(n_coll_type)
              CALL init_collision_type_block(coll_type, ispecies, jspecies)

              ! Open table file and check in case error
              errcode = 0
              OPEN(file_id,FILE=TRIM(full_path_file), STATUS='OLD', &
                IOSTAT=errcode)
              IF (errcode /= 0 .AND. rank==0) THEN
                WRITE(*,*)
                WRITE(*,*) '*** ERROR ***'
                WRITE(*,*) 'Error reading cross section files'
                WRITE(*,*)
                CALL abort_code(c_err_bad_value)
              END IF

              ! Read header lines (threshold energy, type of collision, etc.)
              header_lines = 0
              header = .TRUE.
              DO WHILE (header)
                CALL read_table_header(file_id, header, header_lines, element, &
                  value)
                CALL allocate_table_header(element,value,coll_type,coll_block)
              END DO
              BACKSPACE(file_id) ! The file-reader moves a line back
              header_lines = header_lines - 1

              ! Count length of the energy-cross_section table
              lines = 0
              DO
                READ(file_id,*,IOSTAT=errcode)
                IF (errcode/=0) EXIT
                lines = lines + 1
              END DO
              
              ! Check the number of lines read. Minimum is one line
              IF (rank==0 .AND. lines == 0) THEN
                WRITE(*,*) '*** ERROR ***'
                WRITE(*,*) 'Table in file ', TRIM(ij_type_k_filename), &
                    ' is empty'
                WRITE(*,*)
                CALL abort_code(c_err_io_error)
              END IF

              ! Allocate energy and cross_section arrays and store table data
              coll_type%table_len = lines
              ALLOCATE(coll_type%energy(lines), coll_type%cross_section(lines))
              
              ! Start reading energy-cross-section data and save to arrays
              ! Jump back to the start of the file
              REWIND(file_id)
              
              ! Skip the header lines
              DO l = 1, header_lines
                READ(file_id,*,IOSTAT=errcode)
              END DO
              
              ! Next read table data
              DO l = 1, lines
                READ(file_id,*,IOSTAT=errcode) coll_type%energy(l), &
                  coll_type%cross_section(l)
              END DO

              ! Close file
              CLOSE(file_id)

              ! Fill collision type block
              CALL customise_collision_type_block(coll_block, coll_type, &
                ispecies, jspecies)

              NULLIFY(coll_type)

            ELSE ! If no file is found it is assumed there aren't any more files
              EXIT
            END IF
          END DO ! as many tables as there are
        END DO ! 1st/2nd read

        ! Inter species collision -> link collision block to the pair species
        IF (ispecies /= jspecies .AND. .NOT.is_background) THEN
          IF (.NOT.ASSOCIATED(species_list(jspecies)%neutrals)) THEN
            ALLOCATE(species_list(jspecies)%neutrals(n_species_bg))
            DO species = 1, n_species_bg
              coll_block => species_list(jspecies)%neutrals(species)
              CALL init_dummy_collision_block( coll_block )
            END DO
            coll_block => species_list(ispecies)%neutrals(jspecies)
          END IF
          species_list(jspecies)%neutrals(ispecies) = coll_block
        END IF

        IF (is_background) NULLIFY(background)
        NULLIFY(coll_block)

      END DO ! jspecies
    END DO ! ispecies

    errcode = c_err_none

  END SUBROUTINE load_neutral_collisions



  SUBROUTINE customise_collision_type_block(coll_block, coll_type, &
    ispecies, jspecies)

    TYPE(neutrals_block), POINTER, INTENT(INOUT) :: coll_block
    TYPE(collision_type_block), POINTER, INTENT(INOUT) :: coll_type
    INTEGER, INTENT(IN) :: ispecies, jspecies
    REAL(num) :: mi, mj, mu

    ! Convert energy data into Joules
    coll_type%energy = coll_type%energy * coll_type%energy_units
    coll_type%ethreshold = coll_type%ethreshold * coll_type%energy_units
    ! Convert energy data into speed difference, i.e. g = |u1-u2|
    mi = species_list(ispecies)%mass
    IF (coll_block%is_background) THEN
      mj = coll_block%background%mass
    ELSE
      mj = species_list(jspecies)%mass
    END IF
    mu = mj * mi / (mi + mj)
    coll_type%energy = SQRT(coll_type%energy * 2._num / mu)
    coll_type%gthreshold = SQRT(coll_type%ethreshold * 2._num / mu)


    ! Convert cross-section data into m^2
    coll_type%cross_section = coll_type%cross_section * &
      coll_type%cross_section_units
    ! Convert cross-section data into m^3/s, i.e. g*sigma
    coll_type%cross_section = coll_type%cross_section * coll_type%energy

    ! If collision with a background gas, either Bird or Vahedi
    IF (coll_block%is_background) THEN
      IF (.NOT.coll_type%wnanbu .AND. .NOT.coll_type%wvahedi) THEN
        IF (rank == 0) THEN
          WRITE(*,*)
          WRITE(*,*) '*** WARNING ***'
          WRITE(*,*) "Collisions with background gas are only possible", &
            " with Nanbu's or Vahedi's method"
          WRITE(*,*) 'Collision between ', &
            TRIM(ADJUSTL(species_list(ispecies)%name)), ' and ', &
            TRIM(ADJUSTL(coll_block%background%name))
          WRITE(*,*) ' - Type ', TRIM(ADJUSTL(coll_type%name)), &
            " switched to Nanbu's method"
          WRITE(*,*) ''
        END IF
        coll_type%wnanbu = .TRUE.
        coll_type%wvahedi = .FALSE.
        coll_type%wnanbusplit = .FALSE.
        coll_type%wvahedisplit = .FALSE.
      END IF
    END IF

    ! Link collision subroutine
    CALL link_collision_subroutine(coll_type, coll_block%is_background)

  END SUBROUTINE customise_collision_type_block



  SUBROUTINE link_collision_subroutine(coll_type, is_background)

    TYPE(collision_type_block), POINTER, INTENT(INOUT) :: coll_type
    LOGICAL, INTENT(IN) :: is_background
    LOGICAL :: linked, background_mismatch
    CHARACTER(len=20) :: method_str

    linked = .FALSE.
    background_mismatch = .FALSE.

    ! NANBU COLLISIONS
    IF (coll_type%wnanbu) THEN
      IF (is_background) THEN
        IF ((coll_type%id == c_nc_elastic) .OR. &
          (coll_type%id == c_nc_elastic_electron) .OR. &
          (coll_type%id == c_nc_elastic_ion)) THEN
          coll_type%coll_subroutine => nanbu_elastic_scattering_bg
          linked = .TRUE.
        ELSE IF (coll_type%id == c_nc_excitation) THEN
          coll_type%coll_subroutine => nanbu_excitation_bg
          linked = .TRUE.
        ELSE IF (coll_type%id == c_nc_ionisation) THEN
          coll_type%coll_subroutine => nanbu_ionisation_bg
          linked = .TRUE.
        ELSE IF (coll_type%id == c_nc_charge_exchange) THEN
          coll_type%coll_subroutine => nanbu_charge_exchange_bg
          linked = .TRUE.
        END IF
      ELSE
        IF ((coll_type%id == c_nc_elastic) .OR. &
          (coll_type%id == c_nc_elastic_electron) .OR. &
          (coll_type%id == c_nc_elastic_ion)) THEN
          coll_type%coll_subroutine => nanbu_elastic_scattering
          linked = .TRUE.
        ELSEIF (coll_type%id == c_nc_excitation) THEN
          coll_type%coll_subroutine => nanbu_excitation
          linked = .TRUE.
        ELSEIF (coll_type%id == c_nc_ionisation) THEN
          coll_type%coll_subroutine => nanbu_ionisation
          linked = .TRUE.
        ELSEIF (coll_type%id == c_nc_charge_exchange) THEN
          coll_type%coll_subroutine => nanbu_charge_exchange
          linked = .TRUE.
        END IF
      END IF

    ! NANBU-SPLIT COLLISIONS
    ELSEIF (coll_type%wnanbusplit) THEN
#ifndef PER_SPECIES_WEIGHT
      IF (is_background) THEN
        background_mismatch = .TRUE.
      ELSEIF (coll_type%id == c_nc_ionisation) THEN
        coll_type%coll_subroutine => nanbu_split_ionisation
        linked = .TRUE.
      ELSEIF (coll_type%id == c_nc_charge_exchange) THEN
        coll_type%coll_subroutine => nanbu_split_charge_exchange
        linked = .TRUE.
      END IF
#else
      IF (rank == 0) THEN
        WRITE(*,*)
        WRITE(*,*) '*** ERROR ***'
        WRITE(*,*) 'Collision split method is only possible with', &
          ' per-particle weight set-up.'
        WRITE(*,*) 'Please recompile, and rerun code'
        WRITE(*,*) ''
      END IF
      CALL MPI_BARRIER(comm, errcode)
      CALL abort_code(c_err_bad_setup)
#endif

    ! VAHEDI COLLISIONS
    ELSEIF (coll_type%wvahedi) THEN
      IF (is_background) THEN
        IF (coll_type%id == c_nc_elastic_electron) THEN
          coll_type%coll_subroutine => vahedi_electron_elastic_scattering_bg
          linked = .TRUE.
        ELSEIF (coll_type%id == c_nc_elastic_ion) THEN
          coll_type%coll_subroutine => vahedi_ion_elastic_scattering_bg
          linked = .TRUE.
        ELSEIF (coll_type%id == c_nc_excitation) THEN
          coll_type%coll_subroutine => vahedi_excitation_bg
          linked = .TRUE.
        ELSEIF (coll_type%id == c_nc_ionisation) THEN
          coll_type%coll_subroutine => vahedi_ionisation_bg
          linked = .TRUE.
        ELSEIF (coll_type%id == c_nc_charge_exchange) THEN
          coll_type%coll_subroutine => nanbu_charge_exchange_bg
          linked = .TRUE.
        END IF
      ELSE
        IF (coll_type%id == c_nc_elastic_electron) THEN
          coll_type%coll_subroutine => vahedi_electron_elastic_scattering
          linked = .TRUE.
        ELSEIF (coll_type%id == c_nc_elastic_ion) THEN
          coll_type%coll_subroutine => vahedi_ion_elastic_scattering
          linked = .TRUE.
        ELSEIF (coll_type%id == c_nc_excitation) THEN
          coll_type%coll_subroutine => vahedi_excitation
          linked = .TRUE.
        ELSEIF (coll_type%id == c_nc_ionisation) THEN
          coll_type%coll_subroutine => vahedi_ionisation
          linked = .TRUE.
        ELSEIF (coll_type%id == c_nc_charge_exchange) THEN
          coll_type%coll_subroutine => nanbu_charge_exchange
          linked = .TRUE.
        END IF
      END IF

    ! VAHEDI SPLIT COLLISIONS
    ELSEIF (coll_type%wvahedisplit) THEN
#ifndef PER_SPECIES_WEIGHT
      IF (is_background) THEN
        background_mismatch = .TRUE.
      ELSEIF (coll_type%id == c_nc_ionisation) THEN
        coll_type%coll_subroutine => vahedi_split_ionisation
        linked = .TRUE.
!      ELSEIF (coll_type%id == c_nc_charge_exchange) THEN
!        coll_type%coll_subroutine => nanbu_split_charge_exchange
!        linked = .TRUE.
      END IF
#else
      IF (rank == 0) THEN
        WRITE(*,*)
        WRITE(*,*) '*** ERROR ***'
        WRITE(*,*) 'Collision split method is only possible with', &
          ' per-particle weight set-up.'
        WRITE(*,*) 'Please recompile, and rerun code'
        WRITE(*,*) ''
      END IF
      CALL MPI_BARRIER(comm, errcode)
      CALL abort_code(c_err_bad_setup)
#endif
    END IF

    IF (.NOT.linked) THEN
      IF (rank == 0) THEN
        IF (coll_type%wnanbu) method_str = 'Nanbu'
        IF (coll_type%wvahedi) method_str = 'Vahedi'
        IF (coll_type%wvahedisplit) method_str = 'Vahedi split'
        IF (coll_type%wnanbusplit) method_str = 'Nanbu split'
        WRITE(*,*)
        WRITE(*,*) '*** ERROR ***'
        WRITE(*,*) TRIM(coll_type%name), ' is not available for ', &
          TRIM(ADJUSTL(method_str)), ' method'
        WRITE(*,*) ''
      END IF
      CALL MPI_BARRIER(comm, errcode)
      CALL abort_code(c_err_bad_setup)
    END IF

    IF (background_mismatch) THEN
      IF (rank == 0) THEN
        WRITE(*,*)
        WRITE(*,*) '*** ERROR ***'
        WRITE(*,*) TRIM(coll_type%name), ' method is not available for ', &
          'collisions with background gas'
      END IF
      CALL MPI_BARRIER(comm, errcode)
      CALL abort_code(c_err_bad_setup)
    END IF

  END SUBROUTINE link_collision_subroutine



  SUBROUTINE get_max_weight(coll_block, ispecies, jspecies)

    INTEGER, INTENT(IN) :: ispecies, jspecies
    TYPE(neutrals_block), POINTER, INTENT(INOUT) :: coll_block
#ifndef PER_SPECIES_WEIGHT
    INTEGER(8) :: npart, ipart
    REAL(num) :: all_max_weight, all_min_weight
    REAL(num) :: max_w1, max_w2
    REAL(num) :: min_w1, min_w2
    TYPE(particle), POINTER :: current
#endif
    REAL(num) :: max_weight, min_weight, wi, wj
    
    max_weight = -1._num
    min_weight = 1.e100_num
#ifndef PER_SPECIES_WEIGHT
    current => species_list(ispecies)%attached_list%head
    npart = species_list(ispecies)%attached_list%count
    max_w1 = max_weight
    min_w1 = min_weight
    DO ipart = 1, npart
      wi = current%weight
      max_w1 = MAX(wi, max_w1)
      min_w1 = MIN(wi, min_w1)
      current => current%next
    END DO
    CALL MPI_ALLREDUCE(max_w1,all_max_weight,1,MPIREAL,MPI_MAX,comm,errcode)
    coll_block%max_w1 = all_max_weight
    CALL MPI_ALLREDUCE(min_w1,all_min_weight,1,MPIREAL,MPI_MIN,comm,errcode)
    coll_block%min_w1 = all_min_weight
    
    max_w2 = max_weight
    min_w2 = min_weight
    IF (ispecies /= jspecies) THEN
      IF (jspecies <= n_species) THEN
        current => species_list(jspecies)%attached_list%head
        npart = species_list(jspecies)%attached_list%count
        DO ipart = 1, npart
          wj = current%weight
          max_w2 = MAX(wj, max_w2)
          min_w2 = MIN(wj, min_w2)
          current => current%next
        END DO
      END IF
    END IF
    CALL MPI_ALLREDUCE(max_w2,all_max_weight,1,MPIREAL,MPI_MAX,comm,errcode)
    coll_block%max_w2 = all_max_weight
    CALL MPI_ALLREDUCE(min_w2,all_min_weight,1,MPIREAL,MPI_MIN,comm,errcode)
    coll_block%min_w2 = all_min_weight

    max_weight = MAX(coll_block%max_w1, coll_block%max_w2)
    min_weight = MIN(coll_block%min_w1, coll_block%min_w2)
#else
    wi = species_list(ispecies)%weight
    IF (jspecies > n_species) THEN
      wj = 0._num
    ELSE
      wj = species_list(jspecies)%weight
    END IF
    max_weight = MAX(wi, wj)

    IF (jspecies > n_species) THEN
      wj = 1.e100_num
    ELSE
      wj = species_list(jspecies)%weight
    END IF
    min_weight = MIN(wi, wj)
#endif
    coll_block%max_weight = max_weight
    coll_block%min_weight = min_weight
    
  END SUBROUTINE get_max_weight



  SUBROUTINE set_max_weight_in_coll_block

    INTEGER :: ispecies, jspecies
    TYPE(neutrals_block), POINTER :: coll_block

    DO ispecies = 1, n_species
      DO jspecies = 1, n_species_bg
        IF (.NOT.neutral_coll(ispecies, jspecies)) CYCLE
        coll_block => species_list(ispecies)%neutrals(jspecies)
        CALL get_max_weight(coll_block, ispecies, jspecies)
        NULLIFY(coll_block)
      END DO
    END DO

  END SUBROUTINE set_max_weight_in_coll_block



  SUBROUTINE final_nc_setup

    ! Get max. g*sigma values
    CALL get_gsigma_max

    ! Get max. super-particle weight values
    CALL set_max_weight_in_coll_block

    ! Test neutral collision setup: tables, collision types, etc.
    CALL test_neutral_collision_setup

    ! Print the different features from each collision block and type
    CALL print_collisions

  END SUBROUTINE final_nc_setup



  SUBROUTINE get_gsigma_max
    

    REAL(num) :: g, g_min, g_max
    REAL(num) :: gsigma,gs_min,gs_max,gsigma_total,global_gsigma_max
    REAL(num) :: user_factor
    INTEGER :: ispecies, jspecies, ncolltypes, tlen
    INTEGER :: ctype1, ctype2, line, g_entries, i
    TYPE(neutrals_block), POINTER :: coll_block
    TYPE(collision_type_block), POINTER :: coll_type1, coll_type2

    DO ispecies = 1, n_species
      DO jspecies = 1, n_species_bg
        IF (.NOT.neutral_coll(ispecies, jspecies)) CYCLE
        coll_block => species_list(ispecies)%neutrals(jspecies)

        IF (coll_block%user_gsigma_max) THEN
          coll_block%igsigma_max_total = 1._num/coll_block%gsigma_max_total
          CYCLE
        END IF

        ! Get MAX(g*sigma_total)
        global_gsigma_max = 0._num

        ! Go through all energy values present in this coll_block
        ncolltypes = coll_block%ncolltypes

        ! Loop over collision types
        DO ctype1 = 1, ncolltypes
          coll_type1 => coll_block%collision_set(ctype1)
          g_entries = coll_type1%table_len

          ! Loop over g-values in cross-section tables
          DO i = 1, g_entries
            g = coll_type1%energy(i)

            ! Get gsigma total for the current value of g
            gsigma_total = 0._num

            ! Loop over collision types
            DO ctype2 = 1, ncolltypes

              coll_type2 => coll_block%collision_set(ctype2)
              user_factor = coll_type2%user_factor
              ! Get g-cross-section value
              IF (ctype1 == ctype2) THEN
                gsigma = coll_type1%cross_section(i) * user_factor
              ELSE
                gsigma = -1._num
                tlen = coll_type2%table_len
                IF (tlen == 1) THEN
                  ! Only one line: cross-section is constant
                  gsigma = coll_type2%cross_section(1)
                ELSE
                  ! cross-section is energy dependent
                  DO line = 2, tlen
                    g_min = coll_type2%energy(line-1)
                    g_max = coll_type2%energy(line)
                    IF (g < g_max .AND. g >= g_min ) THEN
                      gs_min = coll_type2%cross_section(line-1)
                      gs_max = coll_type2%cross_section(line)
                      ! Interpolation
                      gsigma = interpolation(g_min, g_max, gs_min, gs_max, g)
                      gsigma = gsigma * user_factor
                      EXIT
                    END IF
                  END DO
                  CALL crosssection_out_of_range(gsigma, g, coll_type2)
                END IF
              END IF

              gsigma_total = gsigma_total + gsigma
            END DO ! collision types 2
            NULLIFY(coll_type2)

            ! Checker whether current g-sigma is max.
            global_gsigma_max = MAX(global_gsigma_max, gsigma_total)

          END DO ! g-entries loop
        END DO ! collision types 1

        ! Update coll_block
        coll_block%gsigma_max_total = global_gsigma_max
        coll_block%igsigma_max_total = 1._num/coll_block%gsigma_max_total

      END DO ! jspecies
    END DO ! ispecies

  END SUBROUTINE get_gsigma_max



  SUBROUTINE read_table_header(file_id, header, header_lines, element, value)
  
    ! This read one line of the file (file_id)
    INTEGER, INTENT(IN) :: file_id
    INTEGER, INTENT(INOUT) :: header_lines
    LOGICAL, INTENT(OUT) :: header
    CHARACTER(len=string_length), INTENT(OUT) :: element, value
    INTEGER :: pos, ic, s, f
    CHARACTER :: c
    CHARACTER(len=string_length) :: buffer
    
    header = .FALSE.
    pos = 1
    ic = 1
    DO WHILE (ic <= string_length)
      READ(file_id, '(A1)', advance='no', size=s, iostat=f, eor=10) c

      ic = ic + 1
      buffer(pos:pos) = c
      
      IF (c == '=') THEN
        element = buffer(1:pos-1)
        header = .TRUE.
        pos = 0
      ELSEIF (c == '#') THEN
        IF (header) THEN
          READ(file_id,*)
          EXIT
        ELSE
          READ(file_id,*)
          pos = 1
          ic = 1
          header_lines = header_lines + 1
          CYCLE
        END IF
      END IF
      
      pos = pos + 1
    END DO
    
10  IF (header) THEN
      value = buffer(1:pos-1)
      header_lines = header_lines + 1
!       WRITE(*,*) TRIM(ADJUSTL(element)), ' = ', TRIM(ADJUSTL(value))
    ELSE
      header_lines = header_lines + 1
      RETURN
    END IF

  END SUBROUTINE read_table_header



  SUBROUTINE allocate_table_header(element, value, coll_type, coll_block)
  
    INTEGER :: errcode, ibg
    LOGICAL :: flag
    REAL(num) :: conversion_factor
    CHARACTER(*), INTENT(IN) :: element, value
    TYPE(background_block), POINTER :: background
    TYPE(collision_type_block), INTENT(INOUT) :: coll_type
    TYPE(neutrals_block), INTENT(INOUT) :: coll_block

    errcode = c_err_none
    flag = .FALSE.

    IF (str_cmp(element, 'energy_threshold')) THEN
      coll_type%ethreshold = as_real_print(value, element, errcode)

    ELSEIF (str_cmp(element, 'output_name')) THEN
      coll_type%io_name = TRIM(ADJUSTL(value))

    ELSEIF (str_cmp(element, 'gsigma_max')) THEN
      coll_block%gsigma_max_total = as_real(value, errcode)
      coll_block%user_gsigma_max = .TRUE.

    ELSEIF (str_cmp(element, 'user_factor')) THEN
      coll_type%user_factor = as_real(value, errcode)

    ELSEIF (str_cmp(element, 'cross_section_units')) THEN
      IF (str_cmp(value, 'cm2')) THEN
        coll_type%cross_section_units = 1.e4_num
      ELSEIF (str_cmp(value, 'cm^2')) THEN
        coll_type%cross_section_units = 1.e4_num
      ELSEIF (str_cmp(value, 'm2')) THEN
        coll_type%cross_section_units = 1._num
      ELSEIF (str_cmp(value, 'm^2')) THEN
        coll_type%cross_section_units = 1._num
      ELSE
        conversion_factor = as_real(value, errcode)
        IF (errcode == c_err_none) &
          coll_type%cross_section_units = conversion_factor
        errcode = c_err_none
      END IF

    ELSEIF (str_cmp(element, 'energy_units')) THEN
      IF (str_cmp(value, 'eV')) THEN
        coll_type%energy_units = q0
      ELSEIF (str_cmp(value, 'ev')) THEN
        coll_type%energy_units = q0
      ELSE
        conversion_factor = as_real(value, errcode)
        IF (errcode == c_err_none) &
          coll_type%energy_units = conversion_factor
        errcode = c_err_none
      END IF

    ELSEIF (str_cmp(element, 'collision_method')) THEN
      IF (str_cmp(value, 'nanbu')) THEN ! This one is default
        coll_type%wnanbu = .TRUE.
      ELSEIF (str_cmp(value, 'vahedi')) THEN
        coll_type%wnanbu = .FALSE.
        coll_type%wvahedi = .TRUE.
      ELSEIF (str_cmp(value, 'vahedisplit')) THEN
        coll_type%wnanbu = .FALSE.
        coll_type%wvahedisplit = .TRUE.
      ELSEIF (str_cmp(value, 'nanbusplit')) THEN
        coll_type%wnanbu = .FALSE.
        coll_type%wnanbusplit = .TRUE.
      END IF

    ELSEIF (str_cmp(element, 'collision_type')) THEN        
      IF (str_cmp(value, 'elastic')) THEN
        coll_type%id = c_nc_elastic
        coll_type%name = TRIM(ADJUSTL(value))
      ELSEIF (str_cmp(value, 'excitation')) THEN
        coll_type%id = c_nc_excitation
        coll_type%name = TRIM(ADJUSTL(value))
      ELSEIF (str_cmp(value, 'ionisation') .OR. &
          str_cmp(value, 'ionization')) THEN
        coll_type%id = c_nc_ionisation
        coll_type%name = TRIM(ADJUSTL(value))
      ELSEIF (str_cmp(value, 'charge exchange') .OR. &
          str_cmp(value, 'charge_exchange')) THEN
        coll_type%id = c_nc_charge_exchange
        coll_type%name = TRIM(ADJUSTL(value))
      ELSEIF (str_cmp(value, 'elastic_electron') .OR. &
          str_cmp(value, 'elastic electron')) THEN
        coll_type%id = c_nc_elastic_electron
        coll_type%name = TRIM(ADJUSTL(value))
      ELSEIF (str_cmp(value, 'elastic_ion') .OR. &
          str_cmp(value, 'elastic ion')) THEN
        coll_type%id = c_nc_elastic_ion
        coll_type%name = TRIM(ADJUSTL(value))
      END IF

    ELSEIF (str_cmp(element, 'source_species')) THEN
      IF (.NOT.str_cmp(value, 'none')) THEN
        DO ibg = 1, n_backgrounds
          background => background_list(ibg)
          IF ( str_cmp(value, TRIM(ADJUSTL(background%name))) ) THEN
            coll_type%source_species_id = n_species + ibg
            flag = .TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT.flag) coll_type%source_species_id = as_integer(value, errcode)
      END IF

    ELSEIF (str_cmp(element, 'new_species')) THEN
      IF (.NOT.str_cmp(value, 'none')) THEN
        DO ibg = 1, n_backgrounds
          background => background_list(ibg)
          IF ( str_cmp(value, TRIM(ADJUSTL(background%name))) ) THEN
            coll_type%new_species_id = n_species + ibg
            flag = .TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT.flag) coll_type%new_species_id = as_integer(value, errcode)
      END IF
    END IF

    IF (errcode /= c_err_none) THEN
      IF (rank == 0) THEN
        WRITE(*,*)
        WRITE(*,*) '*** ERROR ***'
        WRITE(*,*) 'Neutral collision table value "species_source_name"', &
          ' "species_target_name" has not been defined correctly'
        WRITE(*,*) 'Please fix value(s) and rerun code'
        WRITE(*,*) ''
      END IF
      CALL MPI_BARRIER(comm, errcode)
      CALL abort_code(c_err_bad_setup)
    END IF

  END SUBROUTINE allocate_table_header



  SUBROUTINE init_collision_type_block(coll_type, ispecies, jspecies)

    TYPE(collision_type_block), POINTER, INTENT(INOUT) :: coll_type
    INTEGER, INTENT(IN) :: ispecies, jspecies
    
    coll_type%ethreshold = 0._num
    coll_type%gthreshold = 0._num
    coll_type%energy_units = 1._num
    coll_type%cross_section_units = 1._num
    coll_type%user_factor = coll_pairs(ispecies, jspecies)

    coll_type%id = -1
    coll_type%source_species_id = -1
    coll_type%new_species_id = -1
    coll_type%table_len = 0
    
    coll_type%io_name = 'none'

    coll_type%wnanbu = .TRUE.
    coll_type%wvahedi = .FALSE.
    coll_type%wnanbusplit = .FALSE.
    coll_type%wvahedisplit = .FALSE.
    coll_type%coll_subroutine => NULL()

    IF (neutral_collision_counter) THEN
      ALLOCATE(coll_type%coll_counter(nx))
#ifdef PER_SPECIES_WEIGHT
      coll_type%coll_counter = 0
#else
      coll_type%coll_counter = 0._num
#endif
    END IF

  END SUBROUTINE init_collision_type_block



  SUBROUTINE init_collision_block(coll_block, n_coll_types, ispecies, &
    jspecies)
  
    TYPE(neutrals_block), POINTER, INTENT(INOUT) :: coll_block
    INTEGER, INTENT(IN) :: n_coll_types, ispecies, jspecies
    
    ! Initialised the variables from the block
    coll_block%ncolltypes = n_coll_types
    ALLOCATE(coll_block%collision_set(n_coll_types))

    IF (jspecies > n_species) THEN
      ALLOCATE(coll_block%background)
      coll_block%is_background = .TRUE.
    ELSE
      coll_block%is_background = .FALSE.
    END IF

  END SUBROUTINE init_collision_block



  SUBROUTINE init_dummy_collision_block(coll_block)

    TYPE(neutrals_block), POINTER, INTENT(INOUT) :: coll_block
    
    ! Initialised the variables from the block
    coll_block%gsigma_max_total = 0._num
    coll_block%igsigma_max_total = 0._num
    coll_block%user_gsigma_max = .FALSE.
    coll_block%ncolltypes = 0
    coll_block%max_weight = 0._num
#ifndef PER_SPECIES_WEIGHT
    coll_block%max_w1 = 0._num
    coll_block%max_w2 = 0._num
#endif
    coll_block%is_background = .FALSE.


  END SUBROUTINE init_dummy_collision_block



  SUBROUTINE deallocate_neutral_collisions

    INTEGER :: stat

    CALL deallocate_collision_block
    CALL deallocate_backgrounds
    DEALLOCATE(neutral_coll, coulomb_coll, STAT=stat)
    
  END SUBROUTINE deallocate_neutral_collisions



  SUBROUTINE deallocate_collision_block

    TYPE(collision_type_block), POINTER :: coll_type
    TYPE(neutrals_block), POINTER :: coll_block
    INTEGER :: ispecies,jspecies,ctype
    
    DO ispecies = 1, n_species
      IF (ASSOCIATED(species_list(ispecies)%neutrals)) THEN
        DO jspecies = ispecies, n_species_bg
          IF ( .NOT.neutral_coll(ispecies,jspecies) ) CYCLE
          coll_block => species_list(ispecies)%neutrals(jspecies)
          DO ctype = 1, coll_block%ncolltypes
            coll_type => coll_block%collision_set(ctype)
            CALL deallocate_collision_type_block(coll_type)
          END DO ! ctype
          IF (coll_block%is_background) THEN
            NULLIFY(coll_block%background)
          END IF
        END DO ! jspecies
        DEALLOCATE(species_list(ispecies)%neutrals)
      END IF
    END DO ! ispecies

  END SUBROUTINE deallocate_collision_block



  SUBROUTINE deallocate_collision_type_block(coll_type)

    TYPE(collision_type_block), POINTER, INTENT(INOUT) :: coll_type

    DEALLOCATE(coll_type%energy)
    DEALLOCATE(coll_type%cross_section)
    IF(neutral_collision_counter) DEALLOCATE(coll_type%coll_counter)

  END SUBROUTINE deallocate_collision_type_block

#endif
END MODULE nc_setup
