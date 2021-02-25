#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = (current%part_pos - x_grid_min_local) / dx - 0.5_num
#else
        cell_x_r = (current%part_pos - x_grid_min_local) / dx
#endif
        cell_x = FLOOR(cell_x_r + 0.5_num)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_x = cell_x + 1

#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gxfac.inc"
#ifdef ELECTROSTATIC
        IF (x_min_boundary_open .AND. cell_x == 1) THEN
          gx(0) = gx(0) + gx(-1) + gx(-2)
          gx(-1) = 0._num
          gx(-2) = 0._num
        END IF
        IF (x_max_boundary_open .AND. cell_x == nx) THEN
          gx(0) = gx(0) + gx(1) + gx(2)
          gx(1) = 0._num
          gx(2) = 0._num
        END IF
#endif

#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gxfac.inc"
#ifdef ELECTROSTATIC
        IF (x_max_boundary_open .AND. cell_x == nx) THEN
          gx(0) = gx(0) + gx(1)
          gx(1) = 0._num
        END IF
#endif

#else
#include "triangle/gxfac.inc"
#ifdef ELECTROSTATIC
        IF (x_min_boundary_open .AND. cell_x == 1) THEN
          gx(0) = gx(0) + gx(-1)
          gx(-1) = 0._num
        END IF
        IF (x_max_boundary_open .AND. cell_x == nx) THEN
          gx(0) = gx(0) + gx(1)
          gx(1) = 0._num
        END IF
#endif
#endif
