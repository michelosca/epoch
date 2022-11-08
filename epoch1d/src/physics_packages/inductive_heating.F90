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
!------------------------------------------------------------------------------
! Inductive heating
!   - Based on: A. Meige et al., Phys. Plasmas 12, 052317 (2005)
!                https://doi.org/10.1063/1.1897390
!   - Written by M. Osca Engelbrecht
!------------------------------------------------------------------------------

MODULE inductive_heating 
#ifdef ELECTROSTATIC
  USE shared_data
  USE evaluator
  USE strings_advanced
  
  IMPLICIT NONE

  PRIVATE

  TYPE inductive_heating_block 

    ! Input parameters
    INTEGER :: id
    TYPE(primitive_stack) :: function
    REAL(num) :: j0_amp, t_start, t_end
    REAL(num) :: x_min, x_max

    ! Simulation paramters
    REAL(num), DIMENSION(:), ALLOCATABLE :: j_total, j_cond, j_disp
    INTEGER :: ix_min, ix_max

    LOGICAL :: source_on
    LOGICAL :: involved
    LOGICAL :: x_min_sendrecv, x_max_sendrecv

  END TYPE inductive_heating_block 

  TYPE(inductive_heating_block), DIMENSION(:), POINTER :: inductive_sources 
  INTEGER, SAVE :: n_heating_sources = 0
  LOGICAL :: inductive_heating_flag = .FALSE.

  PUBLIC :: inductive_heating_y, inductive_heating_prepare_sources
  PUBLIC :: inductive_heating_init_inductive_block
  PUBLIC :: inductive_heating_allocate_sources_list
  PUBLIC :: inductive_heating_setup, inductive_heating_add_source_to_list
  PUBLIC :: inductive_heating_add_particle_velocity
  PUBLIC :: inductive_heating_flag, inductive_heating_block
  PUBLIC :: inductive_heating_finalize

CONTAINS

  SUBROUTINE inductive_heating_y

    INTEGER :: i
    TYPE(inductive_heating_block), POINTER :: inductive_source

    IF (inductive_heating_flag) THEN
      DO i = 1, n_heating_sources
        inductive_source => inductive_sources(i)
        
        IF (inductive_source%source_on .AND. inductive_source%involved) THEN
          ! Update total current (time and position)
          CALL update_total_current(inductive_source)

          ! Update conduction current
          CALL update_conduction_current(inductive_source)

          ! Update displacement current from total and conduction currents
          CALL update_displacement_current(inductive_source)
        END IF
      END DO

      CALL update_ey
    END IF
  
  END SUBROUTINE inductive_heating_y



  SUBROUTINE update_total_current(inductive_source)

    TYPE(inductive_heating_block),POINTER,INTENT(INOUT)::inductive_source
    REAL(num) :: j0_amp
    INTEGER :: err, i, ix_min, ix_max
    LOGICAL :: return_flag
    TYPE(parameter_pack) :: parameters

    err = 0
    ! The current amplitude defined in the input block is per unit length [A/m]
    j0_amp = inductive_source%j0_amp * dx

    ! Cell boundaries. Ghost cells are discarded
    ix_min = inductive_source%ix_min
    ix_max = inductive_source%ix_max
    
    CALL inductive_heating_check_cell_boundaries(ix_min, ix_max, return_flag)
    IF (return_flag) RETURN

    ! Update total current
    DO i = ix_min, ix_max
      parameters%pack_ix = i
      inductive_source%j_total(i) = j0_amp * evaluate_with_parameters( &
        inductive_source%function, parameters, err)
      print*, step, rank, 'j_total',inductive_source%id, &
        x_min_local + dx*i, inductive_source%j_total(i)
    END DO

  END SUBROUTINE update_total_current 



  SUBROUTINE update_conduction_current(inductive_source)

    TYPE(inductive_heating_block),POINTER,INTENT(INOUT)::inductive_source
    REAL(num) :: idx
    REAL(num), DIMENSION(:), ALLOCATABLE :: temp
    REAL(num), DIMENSION(:), ALLOCATABLE :: j_cond
    INTEGER :: i, j, ix, ix_min, ix_max, n_send, n_recv, start_ix, end_ix

    ix_min = inductive_source%ix_min
    ix_max = inductive_source%ix_max
    ! Spatial factor is still required 
    idx = 1._num/dx
    ALLOCATE(j_cond(ix_min:ix_max))
    j_cond = inductive_source%j_cond * idx

    ! Ghost cells must be shared between processors: send/recv cases
    !_________|_|_|_|  : this line shows the domain boundary and ghost cells 
    !x x x x  |        : the 'x' shows the cells where the heating source is applied
    ! CASE A
    ! PROCESSOR LEFT from BOUNDARY (e.g. ng = 3, nx = 100, ix_max = 99 )
    !_________|_|_|_|   send = MAX(0, ix_max - nx)       ->   MAX(0, 99 - 100 = -1) = 0
    !x x x x  |         recv = MIN(ix_max-(nx-ng),ng)    ->   MIN(99-(100-3) = 2, 3) = 2

    ! PROCESSOR RIGHT from BOUNDARY (e.g. ix_min = -2, ix_max = -1 )
    !   |_|_|_|_______  send = MIN(ix_max-ix_min+1,ng)   ->   MIN(-1-(-2)+1=2, 3) = 2
    !    x x  |         recv = MAX(0, MIN(ix_max,ng))    ->   MAX(0,MIN(-1,3)) = MAX(0,-1) = 0
    
    ! CASE B
    ! PROCESSOR LEFT from BOUNDARY (e.g. ng = 3, nx = 100, ix_max = 102)
    !_________|_|_|_|   send = MAX(0, ix_max - nx)       ->   MAX(0,102 - 100 = 2) = 2
    !x x x x x|x x      recv = MIN(ix_max-(nx-ng),ng)    ->   MIN(102-(100-3)= 5, 3) = 3

    ! PROCESSOR RIGHT from BOUNDARY (e.g. ix_min=-2, ix_max = 2)
    !   |_|_|_|_______  send = MIN(ix_max-ix_min+1,ng)   ->   MIN(2-(-2)+1=5, 3) = 3      
    !    x x x|x x      recv = MAX(0, MIN(ix_max,ng))    ->   MAX(0, MIN(2, 3)) = MAX(0,2) = 2
    
    ! CASE C
    ! PROCESSOR LEFT from BOUNDARY (e.g., ng = 3, nx = 100, ix_max = 103)
    !_________|_|_|_|   send = MAX(0,ix_max-nx)          ->   MAX(0, 103-100=3) = 3
    !x x x x x|x x x    recv = MIN(ix_max-(nx-ng),ng)    ->   MIN(103-(100-3)=6,3) = 3 

    ! PROCESSOR RIGHT from BOUNDARY (e.g. ix_min = -2. ix_max = 50)
    !   |_|_|_|_______  send = MIN(ix_max-ix_min+1,ng)   ->   MIN(50-(-2)+1 = 53, 3) = 3
    !    x x x|x x x x  recv = MAX(0, MIN(ix_max,ng))    ->   MAX(0, MIN(50,3)) = MAX(0,3) = 3
    
    IF (inductive_source%x_min_sendrecv) THEN
      n_send = MIN(ix_max-ix_min+1, ng)
      n_recv = MAX(0, MIN(ix_max, ng))

      print*, step, rank, 'x_min_sendrecv', n_send, n_recv
      ALLOCATE(temp(n_recv))

      ! Send ghost cells at x_min to left processor
      CALL MPI_SENDRECV(j_cond(ix_min), n_send, mpireal, proc_x_min, &
          tag, temp, n_recv, mpireal, proc_x_min, tag, comm, status, errcode)

      ! Load data received from left processor on x_min side (not in ghost cells!)
      DO i = 1, n_recv
        j_cond(i) = j_cond(i) + temp(i)
      END DO
      ! Remove data from ghot cells
      IF (ix_min < 1) j_cond(ix_min:MIN(0,ix_max)) = 0._num
      DEALLOCATE(temp)

    END IF

    IF (inductive_source%x_max_sendrecv) THEN
      n_send = MAX(0, ix_max-nx)
      n_recv = MIN(ix_max-(nx-ng), ng)
      
      print*, step, rank, 'x_max_sendrecv', n_send, n_recv
      ALLOCATE(temp(n_recv))

      ! Send ghost cells at x_max to right processor
      IF (n_send == 0) THEN
        ix = ix_max
      ELSE
        ix = ix_max - n_send + 1
      END IF

      CALL MPI_SENDRECV(j_cond(ix), n_send, mpireal, proc_x_max, &
          tag, temp, n_recv, mpireal, proc_x_max, tag, comm, status, errcode)

      ! Load data received from right processor on x_max side (not in ghost cells!)
      j = 1
      start_ix = MAX(nx-ng+1, ix_min)
      end_ix = nx-ng+n_recv
      DO i = start_ix, end_ix 
        j_cond(i) = j_cond(i) + temp(j)
        j = j + 1
      END DO
      ! Remove data from ghot cells
      IF (ix_max > nx) j_cond(MAX(ix_min,nx+1):ix_max) = 0._num
      DEALLOCATE(temp)
    END IF

    inductive_source%j_cond = j_cond
    DEALLOCATE(j_cond)

    do i = ix_min, ix_max
      print*, step, rank, 'j_cond',inductive_source%id, &
        x_min_local + dx*i, inductive_source%j_cond(i)
    end do

  END SUBROUTINE update_conduction_current



  SUBROUTINE inductive_heating_add_particle_velocity(x_pos,y_vel,charge_weight)
    ! Each particle with position x_pos and y-velocity y_vel is
    ! added to conduction current

    REAL(num), INTENT(IN) :: x_pos, y_vel, charge_weight 
    REAL(num) :: x_min, x_max
    INTEGER :: i, ix
    TYPE(inductive_heating_block), POINTER :: inductive_source
    
    DO i = 1, n_heating_sources
      inductive_source => inductive_sources(i)
      
      ! First check: heating source is on/involved
      IF (.NOT.inductive_source%involved) CYCLE
      IF (.NOT.inductive_source%source_on) CYCLE 

      ! Second check: particle is within the heating source domain
      IF ( (x_pos >= inductive_source%x_min) .AND. &
        (x_pos <= inductive_source%x_max) ) THEN
        
        ! Find particle cell
        DO ix = inductive_source%ix_min, inductive_source%ix_max
          x_min = (ix-1)*dx + x_min_local
          x_max = ix*dx + x_min_local
          IF ( (x_pos >= x_min) .AND. (x_pos <= x_max) ) THEN
            ! Add particle to conduction current
            inductive_source%j_cond(ix) = inductive_source%j_cond(ix) + &
              y_vel * charge_weight
            EXIT
          END IF
        END DO

      END IF
    END DO

  END SUBROUTINE inductive_heating_add_particle_velocity

  

  SUBROUTINE update_displacement_current(ind_source)

    TYPE(inductive_heating_block),POINTER,INTENT(INOUT)::ind_source
    INTEGER :: i, ix_min, ix_max
    LOGICAL :: return_flag

    ix_min = ind_source%ix_min
    ix_max = ind_source%ix_max
    CALL inductive_heating_check_cell_boundaries(ix_min, ix_max, return_flag)
    IF (return_flag) RETURN

    DO i = ix_min, ix_max
      ind_source%j_disp(i) = ind_source%j_total(i) - ind_source%j_cond(i)
      print*, step, rank, 'j_disp',ind_source%id, &
        x_min_local + dx*i, ind_source%j_disp(i)
    END DO

  END SUBROUTINE update_displacement_current



  SUBROUTINE update_ey

    TYPE(inductive_heating_block), POINTER :: inductive_source
    REAL(num), DIMENSION(:), ALLOCATABLE :: dEy_dt
    REAL(num) :: ieps0
    INTEGER :: i, j, ix_min, ix_max
    LOGICAL :: return_flag

    ieps0 = 1._num / epsilon0

    DO j = 1, n_heating_sources
      inductive_source => inductive_sources(j)

      IF (.NOT.inductive_source%involved) CYCLE
      IF (.NOT.inductive_source%source_on) CYCLE

      ! Set cell boundaries.
      ! - Ghost cells are discarded
      ix_min = inductive_source%ix_min
      ix_max = inductive_source%ix_max

      CALL inductive_heating_check_cell_boundaries(ix_min, ix_max, return_flag)
      IF (return_flag) RETURN

      ! Electric field time derivative
      ALLOCATE(dEy_dt(ix_min:ix_max))
      dEy_dt = inductive_source%j_disp(ix_min:ix_max) * ieps0

      ! Update electric field
      DO i = ix_min, ix_max 
        ey(i) = ey(i) + dEy_dt(i) * dt
      END DO
      DEALLOCATE(dEy_dt)
    END DO

  END SUBROUTINE update_ey



  SUBROUTINE inductive_heating_check_cell_boundaries(ix_min,ix_max,return_flag)

    INTEGER, INTENT(INOUT) :: ix_min, ix_max
    LOGICAL, INTENT(INOUT) :: return_flag

      return_flag = .FALSE.
      ! If heating source is ONLY at the ghost region -> RETURN
      IF (ix_min > nx) return_flag = .TRUE. 
      IF (ix_max < 1) return_flag = .TRUE. 

      ! Get rid of cells at the ghost region 
      IF (ix_min<1) ix_min = 1
      IF (ix_max>nx) ix_max = nx

  END SUBROUTINE inductive_heating_check_cell_boundaries

  

  SUBROUTINE inductive_heating_allocate_sources_list(count)
    INTEGER, INTENT(IN) :: count

    ALLOCATE(inductive_sources(count))

  END SUBROUTINE inductive_heating_allocate_sources_list



  SUBROUTINE inductive_heating_init_inductive_block(inductive_source_block)

    TYPE(inductive_heating_block), POINTER :: inductive_source_block
    
    n_heating_sources = n_heating_sources + 1

    inductive_source_block%id = n_heating_sources
    inductive_source_block%x_min = 0._num
    inductive_source_block%x_max = 0._num 
    inductive_source_block%ix_min = -100000000
    inductive_source_block%ix_max = -100000000 
    inductive_source_block%t_start = 0._num
    inductive_source_block%t_end = HUGE(0._num) 
    
    inductive_source_block%j0_amp = 0._num
    
    inductive_source_block%source_on = .FALSE.
    inductive_source_block%involved = .FALSE.
    inductive_source_block%x_min_sendrecv = .FALSE.
    inductive_source_block%x_max_sendrecv = .FALSE.
      
  END SUBROUTINE inductive_heating_init_inductive_block 



  SUBROUTINE inductive_heating_add_source_to_list(inductive_source_block)

    TYPE(inductive_heating_block),POINTER,INTENT(IN) :: inductive_source_block
    INTEGER :: id

    id = inductive_source_block%id
    inductive_sources(id) = inductive_source_block

  END SUBROUTINE inductive_heating_add_source_to_list



  SUBROUTINE inductive_heating_prepare_sources

    INTEGER :: i
    TYPE(inductive_heating_block), POINTER :: inductive_source

    DO i = 1, n_heating_sources
      inductive_source => inductive_sources(i)
      IF (inductive_source%involved) THEN
        IF ((time >= inductive_source%t_start) .AND. &
            (time <= inductive_source%t_end)) THEN
          inductive_source%source_on = .TRUE.
        ELSE
          inductive_source%source_on = .FALSE.
          inductive_source%j_total = 0._num
          inductive_source%j_disp = 0._num
        END IF
        inductive_source%j_cond = 0._num
      END IF
    END DO

  END SUBROUTINE inductive_heating_prepare_sources



  SUBROUTINE inductive_heating_setup 

    INTEGER :: i, j, ix_min, ix_max
    REAL(num) :: source_x_min, source_x_max
    REAL(num) :: x_cell_left, x_cell_right
    TYPE(inductive_heating_block), POINTER :: inductive_source

    IF (inductive_heating_flag) THEN 
    DO j = 1, n_heating_sources
      inductive_source => inductive_sources(j)

        ! Set source START(MIN) cell
        source_x_min = inductive_source%x_min
        source_x_max = inductive_source%x_max
        IF (source_x_min > x_max_local + ng * dx) THEN
          ! In this case the heating source is right of this processor
          ! i.e. this processor does not have a heating source
          inductive_source%involved = .FALSE.
        ELSEIF (source_x_min < x_min_local - ng *dx) THEN
          inductive_source%involved = .TRUE.
          inductive_source%ix_min = 1 - ng
        ELSE
          inductive_source%involved = .TRUE.
          DO i = 1-ng, nx+ng
            x_cell_left = (i-1) * dx + x_min_local
            x_cell_right = i * dx + x_min_local
            IF ((x_cell_left <= source_x_min) .AND. &
              (x_cell_right > source_x_min)) THEN
              inductive_source%ix_min = i
              EXIT
            END IF
          END DO
        END IF

        ! Set source end cell
        IF (source_x_max < x_min_local - ng * dx) THEN
          ! In this case the heating source is left of this processor
          ! i.e. this processor does not have a heating source
          inductive_source%involved = .FALSE.
        ELSEIF (source_x_max > x_max_local + ng *dx) THEN
          inductive_source%involved = .TRUE.
          inductive_source%ix_max = nx + ng
        ELSE
          inductive_source%involved = .TRUE.
          DO i = 1-ng, nx+ng
            x_cell_left = (i-1) * dx + x_min_local
            x_cell_right = i * dx + x_min_local
            IF ((x_cell_left < source_x_max) .AND. &
              (x_cell_right >= source_x_max)) THEN
              inductive_source%ix_max = i
              EXIT
            END IF
          END DO
        END IF

        ! Setup IH block in involved processes 
        IF (inductive_source%involved) THEN
          ix_min = inductive_source%ix_min
          ix_max = inductive_source%ix_max
            
          ! Set flag for MPI send/recv subroutines
          IF (ix_min <= ng .AND. .NOT.x_min_boundary_open) THEN
            inductive_source%x_min_sendrecv = .TRUE.
          END IF
          IF (ix_max >= nx-ng+1 .AND. .NOT.x_max_boundary_open) THEN
            inductive_source%x_max_sendrecv = .TRUE.
          END IF
          
          print*,step, rank, 'setup', ix_min, ix_max, &
            inductive_source%x_min_sendrecv, inductive_source%x_max_sendrecv
          ! Allocate cunrrent arrays
          ALLOCATE(inductive_source%j_cond(ix_min:ix_max))

          IF (ix_min < 1) ix_min = 1
          IF (ix_max > nx) ix_max = nx
          ALLOCATE(inductive_source%j_total(ix_min:ix_max))
          ALLOCATE(inductive_source%j_disp(ix_min:ix_max))

        END IF

      END DO
    END IF

    CALL inductive_heating_prepare_sources
  
  END SUBROUTINE inductive_heating_setup 


  SUBROUTINE inductive_heating_finalize

    TYPE(inductive_heating_block), POINTER :: inductive_source
    INTEGER :: j

    DO j = 1, n_heating_sources
      inductive_source => inductive_sources(j)

      IF (inductive_source%involved) THEN
        DEALLOCATE(inductive_source%j_total)
        DEALLOCATE(inductive_source%j_cond)
        DEALLOCATE(inductive_source%j_disp)
      END IF 
      CALL deallocate_stack(inductive_source%function)

    END DO

    DEALLOCATE(inductive_sources)

  END SUBROUTINE inductive_heating_finalize

#endif
END MODULE inductive_heating