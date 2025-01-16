MODULE VERLET_MOD
    
    USE ISO_FORTRAN_ENV
    
    IMPLICIT NONE

    
    
    !Type used to pass all trajectories of position, velocity, accelerations, and E-field as the return 
    !argument of verlet function. We also want names for when we define these variables in the netCDF writer
    TYPE :: TRAJECTORY
        REAL(REAL64), DIMENSION(:), ALLOCATABLE     :: x_traj, y_traj, vx_traj, vy_traj, ax_traj, ay_traj
        REAL(REAL64), DIMENSION(:, :), ALLOCATABLE  :: E_x, E_y, rho, phi
        CHARACTER(LEN=20)                            :: x_name, y_name, vx_name, vy_name, ax_name, ay_name, Ex_name, Ey_name, &
                                                       rho_name, phi_name
    END TYPE

    CONTAINS

    !Verlet subprogram is for calculating positions and velocities of the particle
    !using the velocity verlet algorithm; also calculates E-field
    FUNCTION VERLET(gauss_seidel, r_init, v_init, dt, steps, n_x, n_y, dx, dy) RESULT(particle_traj)
        INTEGER(INT64) :: n_x, n_y, x_cell, y_cell, steps, i, j
        INTENT(IN) :: steps, n_x, n_y
        REAL(REAL64), INTENT(IN) :: dx, dy, dt
        REAL(REAL64), DIMENSION(0:n_y+1, 0:n_x+1), INTENT(IN) :: gauss_seidel !(n_x+2)x(n_y+2) array of potentials (ghost nodes on boundaries)
        REAL(REAL64), DIMENSION(2), INTENT(IN) :: r_init, v_init
        TYPE(TRAJECTORY) :: particle_traj
        
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
        ALLOCATE(particle_traj%E_x(1:n_y, 1:n_x)) !n_y:1 because we want y idecreasing as we go down array
        ALLOCATE(particle_traj%E_y(1:n_y, 1:n_x))

        !Obtaining the electric field so we can calculate accelerations
        DO i = 1, n_y
            DO j = 1, n_x
                particle_traj%E_x(i, j) = (gauss_seidel(i, j + 1) - gauss_seidel(i, j - 1)) / (2 * dy)
                particle_traj%E_y(i, j) = (gauss_seidel(i + 1, j) - gauss_seidel(i - 1, j)) / (2 * dx)
            END DO
        END DO


        !Initialise acceleration now using cell corresponding to E-field cell we are starting in
        x_cell = FLOOR((r_init(1) - 1.0) / dx) + n_x !This equation and future uses are from the briefing sheet (credit to C Brady & H Ratcliffe)
        y_cell = FLOOR((r_init(2) - 1.0) / dy) + n_y
        particle_traj%ax_traj(0) = - 1. * particle_traj%E_x(y_cell, x_cell)
        particle_traj%ay_traj(0) = - 1. * particle_traj%E_y(y_cell, x_cell)
        
        !Verlet algorithm
        DO i = 1, steps
            particle_traj%x_traj(i)  = particle_traj%x_traj(i - 1) + particle_traj%vx_traj(i - 1) * dt &
                                     + 0.5 * particle_traj%ax_traj(i - 1) * dt ** 2
            particle_traj%y_traj(i)  = particle_traj%y_traj(i - 1) + particle_traj%vy_traj(i - 1) * dt &
                                     + 0.5 * particle_traj%ay_traj(i - 1) * dt ** 2
            x_cell = FLOOR((particle_traj%x_traj(i) - 1.0) / dx) + n_x
            y_cell = FLOOR((particle_traj%y_traj(i) - 1.0) / dy) + n_y
            !We do not want the particle going outside the bounds of the box, so if this occurs we will keep it at its final
            !postition and velocity with this if statement
            IF (particle_traj%x_traj(i) > 1 .OR. particle_traj%x_traj(i) < - 1. .OR. particle_traj%y_traj(i) > 1 .OR. &
                particle_traj%y_traj(i) < - 1.) THEN
                    particle_traj%x_traj(i:) = particle_traj%x_traj(i - 1)
                    particle_traj%y_traj(i:) = particle_traj%y_traj(i - 1)
                    particle_traj%vx_traj(i:) = particle_traj%vx_traj(i - 1)
                    particle_traj%vy_traj(i:) = particle_traj%vy_traj(i - 1)
                    particle_traj%ax_traj(i:) = particle_traj%ax_traj(i - 1)
                    particle_traj%ay_traj(i:) = particle_traj%ay_traj(i - 1)
                    EXIT
            END IF
            particle_traj%ax_traj(i) = -1. * particle_traj%E_x(y_cell, x_cell)
            particle_traj%ay_traj(i) = -1. * particle_traj%E_y(y_cell, x_cell)
            particle_traj%vx_traj(i) = particle_traj%vx_traj(i - 1) + dt * (particle_traj%ax_traj(i) &
                                     + particle_traj%ax_traj(i - 1)) / 2
            particle_traj%vy_traj(i) = particle_traj%vy_traj(i - 1) + dt * (particle_traj%ay_traj(i) &
                                     + particle_traj%ay_traj(i - 1)) / 2
        END DO


    END FUNCTION VERLET

END MODULE VERLET_MOD
