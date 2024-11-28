MODULE NETCDF_WRITER

    USE ISO_FORTRAN_ENV
    USE VERLET_MOD
    USE NETCDF

    IMPLICIT NONE
    
    TYPE :: DIMENSIONS
        INTEGER(INT64)   :: n_x, n_y, steps
        CHARACTER(LEN=*) :: n_x_name, n_y_name, steps_name !These are just to give appropriate names to each dimension
    END TYPE

    CONTAINS
    
    SUBROUTINE WRITER(particle_traj, filename, dims, ierr)
        TYPE(DIMENSIONS), INTENT(IN)                 :: dims
        TYPE(TRAJECTORY), INTENT(IN)                 :: particle_traj !Easiest to use derived type for all trajectories
        INTEGER(INT64)                               :: ierr, file_id
        INTEGER(INT64), DIMENSION(:), ALLOCATABLE    :: dimension_ids, variable_ids
        INTENT(OUT)                                  :: ierr
        CHARACTER(LEN=*), INTENT(IN)                 :: filename

        !We are allocating for storing each id number of the variables and dimensions we have
        ALLOCATE(dimension_ids(3))
        ALLOCATE(variable_ids(10))
        ALLOCATE(particle_traj%x_traj(0:dims%steps))
        ALLOCATE(particle_traj%y_traj(0:dims%steps))
        ALLOCATE(particle_traj%vx_traj(0:dims%steps))
        ALLOCATE(particle_traj%vy_traj(0:dims%steps))
        ALLOCATE(particle_traj%ax_traj(0:dims%steps))
        ALLOCATE(particle_traj%ay_traj(0:dims%steps))
        ALLOCATE(particle_traj%E_x(1:dims%n_x, dims%n_y:1)) !n_y:1 because we want y idecreasing as we go down array
        ALLOCATE(particle_traj%E_y(1:dims%n_x, dims%n_y:1))


        ierr = nf90_create(filename, NF90_CLOBBER, file_id)
        PRINT *, ierr

        !Dimensions n_x_name and n_y_name are for the E-field variables
        ierr = nf90_def_dim(file_id, dims%n_x_name, dims%n_x, dimension_ids(1))
        PRINT *, ierr

        ierr = nf90_def_dim(file_id, dims%n_y_name, dims%n_y, dimension_ids(2))
        PRINT *, ierr

        !Dimension steps_name are for the trajectory variables
        ierr = nf90_def_dim(file_id, dims%steps_name, dims%steps + 1, dimension_ids(3))
        PRINT *, ierr

        !Defining all our variables to store all trajectories, E-fields, charge densities and scalar fields
        ierr = nf90_def_var(file_id, particle_traj%Ex_name, NF90_REAL, dimension_ids(1:2), variable_ids(1))
        PRINT *, ierr

        ierr = nf90_def_var(file_id, particle_traj%Ey_name, NF90_REAL, dimension_ids(1:2), variable_ids(2))
        PRINT *, ierr

        ierr = nf90_def_var(file_id, particle_traj%x_traj, NF90_REAL, dimension_ids(3), variable_ids(3))
        PRINT *, ierr

        ierr = nf90_def_var(file_id, particle_traj%y_traj, NF90_REAL, dimension_ids(3), variable_ids(4))
        PRINT *, ierr

        ierr = nf90_def_var(file_id, particle_traj%vx_traj, NF90_REAL, dimension_ids(3), variable_ids(5))
        PRINT *, ierr

        ierr = nf90_def_var(file_id, particle_traj%vy_traj, NF90_REAL, dimension_ids(3), variable_ids(6))
        PRINT *, ierr

        ierr = nf90_def_var(file_id, particle_traj%ax_traj, NF90_REAL, dimension_ids(3), variable_ids(7))
        PRINT *, ierr

        ierr = nf90_def_var(file_id, particle_traj%ay_traj, NF90_REAL, dimension_ids(3), variable_ids(8))
        PRINT *, ierr

        !These two remaining variables are for charge density and scalar potential
        !ierr = nf90_def_var(file_id, NAME, NF90_REAL, dimension_ids(1:2), variable_ids(9))
        !PRINT *, ierr

        !ierr = nf90_def_var(file_id, NAME, NF90_REAL, dimension_ids(1:2), variable_ids(10))
        !PRINT *, ierr



        ierr = nf90_close(file_id)
        PRINT *, ierr


    END SUBROUTINE WRITER

END MODULE NETCDF_WRITER