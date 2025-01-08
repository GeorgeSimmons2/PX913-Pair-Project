PROGRAM MAIN
    USE ISO_FORTRAN_ENV
    USE VERLET_MOD
    USE grid_initialisation
    USE NETCDF_WRITER
    USE command_line
    IMPLICIT NONE

    TYPE(TRAJECTORY) :: trajectory_data
    TYPE(DIMENSIONS) :: dimension_data
    INTEGER          :: nx, ny
    INTEGER(INT64)   :: steps = 1000,ierr
    REAL(REAL64), DIMENSION(:, :), ALLOCATABLE :: rho, phi
    REAL(REAL64)     :: dx, dy
    CHARACTER(LEN=10):: problem, filename = 'traj.nc'
    REAL(REAL64), DIMENSION(2) :: r_init,v_init
    REAL(REAL64)     :: dt = 0.01
    LOGICAL          :: arg

    CALL parse_args

    arg = get_arg('nx', nx)
    arg = get_arg('ny', ny)
    arg = get_arg('problem', problem)

    dimension_data%steps = steps
    dimension_data%n_x = nx
    dimension_data%n_y = ny
    dimension_data%n_x_name = 'x_axis'
    dimension_data%n_y_name = 'y_axis'
    dimension_data%steps_name = 'Time'

    
   
    ! TODO: Include run data in the netCDF file, i.e number of steps, run file name, problem type, etc.

    ! TODO: Final finishing touches on the plotter file. 
    CALL Initial_Conditions(problem,r_init,v_init)
    CALL define_rho(rho, nx, ny, problem, dx, dy)
    CALL calc_potential(rho, phi, dx, dy)

 

    trajectory_data = VERLET(phi, r_init, v_init, dt, steps, INT(nx,INT64), INT(ny,INT64), dx, dy)
    trajectory_data%x_name = 'x_trajectory'
    trajectory_data%y_name = 'y_trajectory'
    trajectory_data%vx_name = 'x_velocity'
    trajectory_data%vy_name = 'y_velocity'
    trajectory_data%ax_name = 'x_acceleration'
    trajectory_data%ay_name = 'y_acceleration'
    trajectory_data%Ex_name = 'x_Electric_field'
    trajectory_data%Ey_name = 'y_Electric_field'
    trajectory_data%rho_name = 'Charge_density'
    trajectory_data%phi_name = 'Electric_potential'

    CALL WRITER(trajectory_data, rho, phi, filename, dimension_data, ierr, INT(nx,INT64),INT(ny,INT64))

END PROGRAM MAIN