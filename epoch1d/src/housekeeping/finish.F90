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

MODULE finish

  USE constants
  USE diagnostics
  USE setup
  USE deck
  USE window
  USE laser
  USE collisions
  USE nc_setup
  USE dist_fn
  USE ionise
  USE injectors
  USE probes
#ifdef ELECTROSTATIC
  USE electrostatic
#endif
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: finalise

CONTAINS

  SUBROUTINE finalise

    CALL close_files
    CALL close_psd_diagnostics_files
    IF (done_mpi_initialise) CALL deallocate_memory
    CALL MPI_FINALIZE(errcode)
    STOP

  END SUBROUTINE finalise



  SUBROUTINE deallocate_memory

    INTEGER :: i, n, stat
    CLASS(particle_id_hash), POINTER :: current_hash

    DEALLOCATE(x, xb, x_global, xb_global, xb_offset_global)
    DEALLOCATE(ex, ey, ez, bx, by, bz, jx, jy, jz)

    DEALLOCATE(npart_each_rank)
    DEALLOCATE(x_grid_mins, x_grid_maxs, cell_x_min, cell_x_max)

    DEALLOCATE(total_particle_energy_species)
    CALL deallocate_probes

    DO i = 1, n_species
      CALL deallocate_stack(species_list(i)%density_function)
      DO n = 1, 3
        CALL deallocate_stack(species_list(i)%temperature_function(n))
        CALL deallocate_stack(species_list(i)%drift_function(n))
        CALL deallocate_stack(species_list(i)%dist_fn_range(n))
      END DO
      CALL destroy_partlist(species_list(i)%attached_list)
      DEALLOCATE(species_list(i)%ext_temp_x_min, STAT=stat)
      DEALLOCATE(species_list(i)%ext_temp_x_max, STAT=stat)
      IF (ASSOCIATED(species_list(i)%background_density)) &
          DEALLOCATE(species_list(i)%background_density, STAT=stat)
    END DO

    CALL deallocate_neutral_collisions

    DEALLOCATE(species_list, STAT=stat)

    DO i = 1, n_io_blocks
      IF (ASSOCIATED(io_block_list(i)%dump_at_times)) &
          DEALLOCATE(io_block_list(i)%dump_at_times, STAT=stat)
      IF (ASSOCIATED(io_block_list(i)%dump_at_nsteps)) &
          DEALLOCATE(io_block_list(i)%dump_at_nsteps, STAT=stat)
    END DO
    DEALLOCATE(io_block_list, STAT=stat)
    DEALLOCATE(io_list_data, STAT=stat)
    DEALLOCATE(file_prefixes, STAT=stat)
    DEALLOCATE(file_numbers, STAT=stat)

    DO i = 1, n_subsets
      DEALLOCATE(subset_list(i)%dumpmask, STAT=stat)
      DEALLOCATE(subset_list(i)%use_species, STAT=stat)
      IF (subset_list(i)%persistent) THEN
        current_hash => id_registry%get_hash(subset_list(i)%name)
        DEALLOCATE(current_hash)
      END IF
      IF (.NOT. subset_list(i)%time_varying) CYCLE
      DO n = 1, c_subset_max
        IF (subset_list(i)%use_restriction_function(n)) THEN
          CALL deallocate_stack(subset_list(i)%restriction_function(n))
        END IF
      END DO
    END DO
    DEALLOCATE(subset_list, STAT=stat)
    CALL id_registry%reset

    DO i = 1, n_deck_constants
      CALL deallocate_stack(deck_constant_list(i)%execution_stream)
    END DO
    DEALLOCATE(deck_constant_list, STAT=stat)

    CALL deallocate_input_deck_buffer
    CALL deallocate_window
    CALL deallocate_lasers
    CALL deallocate_collisions
    CALL deallocate_file_list
    CALL deallocate_dist_fns
    CALL deallocate_ionisation
    CALL deallocate_partlists
    CALL deallocate_eval_stack
    CALL deallocate_injectors
#ifdef ELECTROSTATIC
    CALL deallocate_potentials
    CALL deallocate_efield_profiles
    CALL destroy_petsc
    CALL PetscFinalize(perr)
#endif
    CALL MPI_COMM_FREE(comm, errcode)

  END SUBROUTINE deallocate_memory



  SUBROUTINE close_psd_diagnostics_files

    INTEGER :: i, dump_id, species_counter, file_unit

    IF (psd_diag_rank == rank) THEN
      species_counter = 1
      DO i = 1, psd_diag_n_dumps
        dump_id = psd_diag_dump_id(i)
        IF (dump_id == c_dump_number_density) THEN
          file_unit = dump_id + 100 + species_counter
          species_counter = species_counter + 1
        ELSE
          file_unit = dump_id + 100
        END IF

        CLOSE(unit=file_unit)

      END DO
    END IF

 END SUBROUTINE close_psd_diagnostics_files

END MODULE finish
