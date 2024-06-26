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

PROGRAM pic

  ! EPOCH1D is a Birdsall and Langdon type PIC code derived from the PSC
  ! written by Hartmut Ruhl.

  ! The particle pusher (particles.F90) and the field solver (fields.f90) are
  ! almost exact copies of the equivalent routines from PSC, modified slightly
  ! to allow interaction with the changed portions of the code and for
  ! readability. The MPI routines are exactly equivalent to those in PSC, but
  ! are completely rewritten in a form which is easier to extend with arbitrary
  ! fields and particle properties. The support code is entirely new and is not
  ! equivalent to PSC.

  ! EPOCH1D written by C.S.Brady, Centre for Fusion, Space and Astrophysics,
  ! University of Warwick, UK
  ! PSC written by Hartmut Ruhl

  USE balance
  USE deck
  USE diagnostics
  USE fields
  USE helper
  USE ic_module
  USE mpi_routines
  USE particles
  USE setup
  USE finish
  USE welcome
  USE window
  USE split_particle
  USE collisions
  USE nc_setup
  USE particle_migration
  USE ionise
  USE calc_df
  USE injectors
  USE current_smooth
#ifdef PHOTONS
  USE photons
#endif
#ifdef BREMSSTRAHLUNG
  USE bremsstrahlung
#endif
#ifdef ELECTROSTATIC
  USE electrostatic
  USE inductive_heating
#endif
#ifdef NEUTRAL_COLLISIONS
  USE neutral_collisions
#endif

  IMPLICIT NONE

  INTEGER :: ispecies, ierr
  LOGICAL :: halt = .FALSE., push = .TRUE.
  LOGICAL :: force_dump = .FALSE.
  CHARACTER(LEN=64) :: deck_file = 'input.deck'
  CHARACTER(LEN=*), PARAMETER :: data_dir_file = 'USE_DATA_DIRECTORY'
  CHARACTER(LEN=64) :: timestring
  REAL(num) :: runtime
#ifndef ELECTROSTATIC
  REAL(num) :: dt_store
#endif

  step = 0
  time = 0.0_num

  CALL mpi_minimal_init ! mpi_routines.f90
  real_walltime_start = MPI_WTIME()
  CALL minimal_init     ! setup.f90
  CALL get_job_id(jobid)
  CALL welcome_message  ! welcome.f90

  IF (rank == 0) THEN
    OPEN(unit=lu, status='OLD', file=TRIM(data_dir_file), iostat=ierr)
    IF (ierr == 0) THEN
      READ(lu,'(A)') data_dir
      CLOSE(lu)
      PRINT*, 'Using data directory "' // TRIM(data_dir) // '"'
    ELSE
      PRINT*, 'Specify output directory'
      READ(*,'(A)') data_dir
    END IF
    CALL cleanup_stop_files
  END IF

  CALL MPI_BCAST(data_dir, c_max_path_length, MPI_CHARACTER, 0, comm, errcode)

  ! version check only, exit silently
  IF (TRIM(data_dir) == 'VERSION_INFO') CALL finalise

  CALL register_objects ! custom.f90
  CALL read_deck(deck_file, .TRUE., c_ds_first)

  CALL setup_partlists  ! partlist.f90
  CALL timer_init

  IF (use_exact_restart) CALL read_cpu_split
  CALL setup_boundaries ! boundary.f90
  CALL mpi_initialise  ! mpi_routines.f90
  CALL after_control   ! setup.f90
  CALL open_files      ! setup.f90

  ! Re-scan the input deck for items which require allocated memory
  CALL read_deck(deck_file, .TRUE., c_ds_last)
  CALL after_deck_last

  ! restart flag is set
  IF (ic_from_restart) THEN
    CALL restart_data(step)
  ELSE
    ! auto_load particles
    CALL pre_load_balance
    CALL auto_load
    time = 0.0_num
  END IF

  CALL custom_particle_load
  CALL manual_load
  CALL finish_injector_setup
#ifdef ELECTROSTATIC
#ifdef PETSC
  CALL setup_petsc_variables(nx)
#endif
#endif


#ifdef NEUTRAL_COLLISIONS
  IF (ANY(neutral_coll)) CALL final_nc_setup
#endif
  IF (inductive_heating_flag) CALL inductive_heating_setup
  
  CALL initialise_window ! window.f90
  CALL set_dt
  CALL set_maxwell_solver
  CALL deallocate_ic
  CALL update_particle_count

  npart_global = 0
  DO ispecies = 1, n_species
    npart_global = npart_global + species_list(ispecies)%count
  END DO

  ! .TRUE. to over_ride balance fraction check
  IF (npart_global > 0) CALL balance_workload(.TRUE.)

  IF (use_current_correction) CALL calc_initial_current
  CALL particle_bcs
#ifndef ELECTROSTATIC
  CALL efield_bcs
#endif

  IF (ic_from_restart) THEN
    IF (dt_from_restart > 0) dt = dt_from_restart
    IF (step == 0) THEN
      CALL bfield_final_bcs
    ELSE
#ifndef ELECTROSTATIC
      time = time + dt / 2.0_num
      CALL update_eb_fields_final
      CALL moving_window
#else
      CALL es_initialize_e_field
#endif
    END IF
  ELSE
#ifndef ELECTROSTATIC
    dt_store = dt
    dt = dt / 2.0_num
    time = time + dt
    CALL bfield_final_bcs
    dt = dt_store
#else
    CALL bfield_final_bcs
    CALL es_initialize_e_field
    CALL es_particles_push_back
#endif
  END IF
  CALL count_n_zeros

  ! Setup particle migration between species
  IF (use_particle_migration) CALL initialise_migration
  CALL build_persistent_subsets
#ifdef PHOTONS
  IF (use_qed) CALL setup_qed_module()
#endif
#ifdef BREMSSTRAHLUNG
  IF (use_bremsstrahlung) CALL setup_bremsstrahlung_module()
#endif

  IF (rank == 0) THEN
    PRINT*
    PRINT*, 'Equilibrium set up OK, running code'
    PRINT*
  END IF

  walltime_started = MPI_WTIME()
  IF (.NOT.ic_from_restart) CALL output_routines(step) ! diagnostics.f90
  IF (use_field_ionisation) CALL initialise_ionisation

  IF (timer_collect) CALL timer_start(c_timer_step)

  DO
    IF ((step >= nsteps .AND. nsteps >= 0) &
        .OR. (time >= t_end) .OR. halt) EXIT

    IF (timer_collect) THEN
      CALL timer_stop(c_timer_step)
      CALL timer_reset
      timer_first(c_timer_step) = timer_walltime
    END IF

    push = (time >= particle_push_start_time)
#ifdef PHOTONS
    IF (push .AND. use_qed .AND. time > qed_start_time) THEN
      CALL qed_update_optical_depth()
    END IF
#endif

#ifdef BREMSSTRAHLUNG
    IF (push .AND. use_bremsstrahlung &
        .AND. time > bremsstrahlung_start_time) THEN
      CALL bremsstrahlung_update_optical_depth()
    END IF
#endif

#ifndef ELECTROSTATIC
    CALL update_eb_fields_half
#endif
    IF (push) THEN
      CALL run_injectors
      ! .FALSE. this time to use load balancing threshold
      IF (use_balance) CALL balance_workload(.FALSE.)
#ifdef ELECTROSTATIC
      CALL es_wall_current_diagnostic
      CALL es_push_particles
      CALL es_update_e_field
#else
      CALL push_particles
#endif
      IF (use_particle_lists) THEN
        ! After this line, the particles can be accessed on a cell by cell basis
        ! Using the particle_species%secondary_list property
        CALL reorder_particles_to_grid

        ! call collision operator
        IF (use_collisions) THEN
          IF (use_collisional_ionisation) THEN
            CALL collisional_ionisation
          ELSE
#ifdef NEUTRAL_COLLISIONS
            IF (ANY(coulomb_coll)) CALL particle_collisions
            IF (ANY(neutral_coll)) CALL run_neutral_collisions
#else
            CALL particle_collisions
#endif
          END IF
        END IF

        ! Early beta version of particle splitting operator
        IF (use_split) CALL split_particles

        CALL reattach_particles_to_mainlist
      END IF
      IF (use_particle_migration) CALL migrate_particles(step)
      IF (use_field_ionisation) CALL ionise_particles
#ifndef ELECTROSTATIC
      CALL current_finish
#endif
      CALL update_particle_count
    END IF

    CALL check_for_stop_condition(halt, force_dump)
    IF (halt) EXIT
    step = step + 1
#ifndef ELECTROSTATIC
    time = time + dt / 2.0_num
    CALL output_routines(step)
    time = time + dt / 2.0_num
    CALL update_eb_fields_final
    CALL moving_window
#else
    time = time + dt
    CALL output_routines(step)
#endif
  END DO

  IF (rank == 0) runtime = MPI_WTIME() - walltime_started

#ifdef PHOTONS
  IF (use_qed) CALL shutdown_qed_module()
#endif

#ifdef BREMSSTRAHLUNG
  IF (use_bremsstrahlung) CALL shutdown_bremsstrahlung_module()
#endif

  CALL output_routines(step, force_dump)

  IF (rank == 0) THEN
    CALL create_full_timestring(runtime, timestring)
    WRITE(*,*) 'Final runtime of core = ' // TRIM(timestring)
  END IF

  CALL finalise

END PROGRAM pic
