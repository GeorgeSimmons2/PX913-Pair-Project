MODULE VERLET_MOD
    
    USE ISO_FORTRAN_ENV
    
    IMPLICIT NONE

    CONTAINS
    
    !Verlet subprogram is for calculating positions and velocities of the particle
    !using the velocity verlet algorithm
    SUBROUTINE VERLET(gauss_seidel, r_init, v_init, dt, steps, n_x, n_y, dx, dy)
        INTEGER(INT64) :: n_x, n_y, x_cell, y_cell, steps
        REAL(REAL64), INTENT(IN) :: dx, dy, dt
        REAL(REAL64) :: E_x, E_y
        REAL(REAL64), DIMENSION(n_x + 2, n_y + 2), INTENT(IN) :: gauss_seidel
        REAL(REAL64), DIMENSION(2), INTENT(IN) :: r_init, v_init
        INTENT(IN) :: steps, n_x, n_y

        x_cell = FLOOR((r_init(1) - 1.0) / dx) + 1
        y_cell = FLOOR((r_init(2) - 1.0) / dy) + 1

    END SUBROUTINE VERLET



END MODULE VERLET_MOD