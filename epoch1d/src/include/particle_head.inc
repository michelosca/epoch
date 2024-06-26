    REAL(num), DIMENSION(sf_min:sf_max) :: gx
    REAL(num) :: cell_x_r, cell_frac_x
    INTEGER :: cell_x
#ifndef PARTICLE_SHAPE_TOPHAT
    REAL(num) :: cx2
#endif
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num), PARAMETER :: third = 1.0_num / 3.0_num
    REAL(num), PARAMETER :: fac1 = 0.125_num * third
    REAL(num), PARAMETER :: fac2 = 0.5_num * third
    REAL(num), PARAMETER :: fac3 = 7.1875_num * third
#endif
#ifdef ELECTROSTATIC
    INTEGER :: j
#ifndef PARTICLE_SHAPE_TOPHAT
    INTEGER :: minus_sf_min = -sf_min
#endif
#endif
