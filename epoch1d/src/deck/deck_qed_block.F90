MODULE deck_qed_block

  USE mpi
  USE strings_advanced

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: qed_deck_initialise, qed_deck_finalise
  PUBLIC :: qed_block_start, qed_block_end
  PUBLIC :: qed_block_handle_element, qed_block_check

CONTAINS

  SUBROUTINE qed_deck_initialise

#ifdef PHOTONS
    IF (deck_state .EQ. c_ds_first) THEN
      qed_table_location = 'src/physics_packages/TABLES'
    ENDIF
#endif

  END SUBROUTINE qed_deck_initialise



  SUBROUTINE qed_deck_finalise

#ifdef PHOTONS
    LOGICAL :: exists
    INTEGER :: io, ierr

    IF (deck_state .EQ. c_ds_first) RETURN

    IF (rank .EQ. 0) THEN
      INQUIRE(file=TRIM(qed_table_location)//'/hsokolov.table', exist=exists)
      IF (.NOT.exists) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to find QED tables in the ', &
              'directory "' // TRIM(qed_table_location) // '"'
        ENDDO
        CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      ENDIF
    ENDIF
#endif

  END SUBROUTINE qed_deck_finalise



  SUBROUTINE qed_block_start

  END SUBROUTINE qed_block_start



  SUBROUTINE qed_block_end

  END SUBROUTINE qed_block_end



  FUNCTION qed_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode

    errcode = c_err_none
    IF (deck_state .EQ. c_ds_first) RETURN
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

#ifdef PHOTONS
    IF (str_cmp(element, 'use_qed') .OR. str_cmp(element, 'qed')) THEN
      use_qed = as_logical(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'qed_start_time')) THEN
      qed_start_time = as_real(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'produce_photons')) THEN
      produce_photons = as_logical(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'photon_energy_min') &
        .OR. str_cmp(element, 'min_photon_energy')) THEN
      photon_energy_min = as_real(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'produce_pairs')) THEN
      produce_pairs = as_logical(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'qed_table_location')) THEN
      qed_table_location = TRIM(ADJUSTL(value))
      RETURN
    ENDIF

    IF (str_cmp(element, 'photon_dynamics')) THEN
      photon_dynamics = as_logical(value, errcode)
      RETURN
    ENDIF

    errcode = c_err_unknown_element
#else
    extended_error_string = '-DPHOTONS'
    errcode = c_err_pp_options_wrong
#endif

  END FUNCTION qed_block_handle_element



  FUNCTION qed_block_check() RESULT(errcode)

    INTEGER :: errcode, io

    errcode = c_err_none

#ifdef PHOTONS
    IF (produce_pairs .AND. .NOT. photon_dynamics) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'You cannot set photon_dynamics=F when ', &
            'produce_pairs=T. Without ', 'photon motion, pair ', &
            'creation will be incorrect.'
          WRITE(io,*) 'Code will terminate.'
        ENDDO
      ENDIF
      errcode = c_err_bad_value + c_err_terminate
    ENDIF
#endif

  END FUNCTION qed_block_check

END MODULE deck_qed_block