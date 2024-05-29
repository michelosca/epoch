#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = (parameters%pack_pos - x_grid_min_local) / dx - 0.5_num
#else
        cell_x_r = (parameters%pack_pos - x_grid_min_local) / dx
#endif
        cell_x = FLOOR(cell_x_r + 0.5_num)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_x = cell_x + 1

#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gxfac.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gxfac.inc"
#else
#include "triangle/gxfac.inc"
#endif

#ifdef ELECTROSTATIC
    IF (x_max_boundary_open .AND. cell_x == nx) THEN
      DO j = 1, sf_max
        gx(0) = gx(0) + gx(j)
        gx(j) = 0._num
      END DO
    END IF
#ifdef PARTICLE_SHAPE_TOPHAT
    IF (x_min_boundary_open .AND. cell_x == 0) THEN
      gx(1) = gx(1) + gx(0)
      gx(0) = 0._num
#else
    IF (x_min_boundary_open .AND. cell_x == 1) THEN
      DO j = 1, minus_sf_min
        gx(0) = gx(0) + gx(-j)
        gx(-j) = 0._num
      END DO
#endif
    END IF
#endif
