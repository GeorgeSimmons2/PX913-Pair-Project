MODULE VERLET_MOD
    
    USE ISO_FORTRAN_ENV
    
    IMPLICIT NONE

    
    
    !Type used to pass all trajectories of position, velocity, accelerations, and E-field as the return 
    !argument of verlet function
    TYPE :: trajectory
        REAL(REAL64), DIMENSION(:), ALLOCATABLE :: x_traj, y_traj, vx_traj, vy_traj, ax_traj, ay_traj
        REAL(REAL64), DIMENSION(:, :), ALLOCATABLE :: E_x, E_y
    END TYPE

    CONTAINS

    !Verlet subprogram is for calculating positions and velocities of the particle
    !using the velocity verlet algorithm; also calculates E-field
    FUNCTION VERLET(gauss_seidel, r_init, v_init, dt, steps, n_x, n_y, dx, dy) RESULT(particle_traj)
        INTEGER(INT64) :: n_x, n_y, x_cell, y_cell, steps, i, j
        INTENT(IN) :: steps, n_x, n_y
        REAL(REAL64), INTENT(IN) :: dx, dy, dt
        REAL(REAL64), DIMENSION(n_x, n_y) :: E_x, E_y
        REAL(REAL64), DIMENSION(n_x + 2, n_y + 2), INTENT(IN) :: gauss_seidel !(n_x+2)x(n_y+2) array of potentials (ghost nodes on boundaries)
        REAL(REAL64), DIMENSION(2), INTENT(IN) :: r_init, v_init
        TYPE(trajectory) :: particle_traj
        
        !Allocate and initialise the trajectories with the initial conditions
        !(we will initialise acceleration after calculating the E-field)
        ALLOCATE(particle_traj%x_traj(0:steps))
        particle_traj%x_traj(0) = r_init(1)
        ALLOCATE(particle_traj%y_traj(0:steps))
        particle_traj%y_traj(0) = r_init(2)
        ALLOCATE(particle_traj%vx_traj(0:steps))
        particle_traj%vx_traj(0) = v_init(1)
        ALLOCATE(particle_traj%vy_traj(0:steps))
        particle_traj%vy_traj(0) = v_init(2)
        ALLOCATE(particle_traj%ax_traj(0:steps))
        ALLOCATE(particle_traj%ay_traj(0:steps))
        ALLOCATE(particle_traj%E_x(n_x, n_y))
        ALLOCATE(particle_traj%E_y(n_x, n_y))

        !Obtaining the electric field so we can calculate accelerations
        DO i = 1, n_x
            DO j = 1, n_y
                particle_traj%E_x(i, j) = (gauss_seidel(i + 1, j) - gauss_seidel(i - 1, j) / (2 * dx))
                particle_traj%E_y(i, j) = (gauss_seidel(i, j + 1) - gauss_seidel(i, j - 1) / (2 * dy))
            END DO
        END DO

        !Initialise acceleration now using cell corresponding to E-field cell we are starting in
        x_cell = FLOOR((r_init(1) - 1.0) / dx) + 1
        y_cell = FLOOR((r_init(2) - 1.0) / dy) + 1
        particle_traj%ax_traj(0) = - 1. * particle_traj%E_x(x_cell, y_cell)
        particle_traj%ay_traj(0) = - 1. * particle_traj%E_y(x_cell, y_cell)
        
        !Verlet algorithm
        


    END FUNCTION VERLET



END MODULE VERLET_MOD