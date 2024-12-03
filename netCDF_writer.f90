MODULE NETCDF_WRITER

    USE ISO_FORTRAN_ENV
    USE VERLET_MOD
    USE NETCDF

    IMPLICIT NONE
    
    TYPE :: DIMENSIONS
        INTEGER(INT32)   :: n_x, n_y, steps
        CHARACTER(LEN=30) :: n_x_name, n_y_name, steps_name !These are just to give appropriate names to each dimension
    END TYPE

    CONTAINS
    
    SUBROUTINE WRITER(particle_traj, rho, phi, filename, dims, ierr, N)
        TYPE(DIMENSIONS), INTENT(INOUT)                 :: dims
        TYPE(TRAJECTORY), INTENT(INOUT)                 :: particle_traj !Easiest to use derived type for all trajectories
        INTEGER(INT32)                                  :: file_id
        INTEGER(INT32), DIMENSION(:), ALLOCATABLE       :: dimension_ids, variable_ids
        REAL(REAL64), DIMENSION(:, :), ALLOCATABLE      :: rho, phi, rho_corrected, phi_corrected
        INTENT(INOUT)                                   :: rho, phi
        INTEGER(INT64), INTENT(OUT)                     :: ierr
        CHARACTER(LEN=*), INTENT(IN)                    :: filename
        INTEGER(INT64)                                  :: i, j, N

        ALLOCATE(rho_corrected(N, N))
        ALLOCATE(phi_corrected(N, N))

        DO i = 2, N-1
            DO j = 2, N-1
                rho_corrected(i, j) = rho(i, j)
                phi_corrected(i, j) = phi(i, j)
            END DO
        END DO

        !We are allocating for storing each id number of the variables and dimensions we have
        ALLOCATE(dimension_ids(3))
        ALLOCATE(variable_ids(10))

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

        ierr = nf90_def_var(file_id, particle_traj%x_name, NF90_REAL, dimension_ids(3), variable_ids(3))
        PRINT *, ierr

        ierr = nf90_def_var(file_id, particle_traj%y_name, NF90_REAL, dimension_ids(3), variable_ids(4))
        PRINT *, ierr

        ierr = nf90_def_var(file_id, particle_traj%vx_name, NF90_REAL, dimension_ids(3), variable_ids(5))
        PRINT *, ierr

        ierr = nf90_def_var(file_id, particle_traj%vy_name, NF90_REAL, dimension_ids(3), variable_ids(6))
        PRINT *, ierr

        ierr = nf90_def_var(file_id, particle_traj%ax_name, NF90_REAL, dimension_ids(3), variable_ids(7))
        PRINT *, ierr

        ierr = nf90_def_var(file_id, particle_traj%ay_name, NF90_REAL, dimension_ids(3), variable_ids(8))
        PRINT *, ierr

        ierr = nf90_def_var(file_id, particle_traj%rho_name, NF90_REAL, dimension_ids(1:2), variable_ids(9))
        PRINT *, ierr

        ierr = nf90_def_var(file_id, particle_traj%phi_name, NF90_REAL, dimension_ids(1:2), variable_ids(10))
        PRINT *, ierr

        !Also add attributes to variables/global
        ierr = nf90_enddef(file_id)
        PRINT *, ierr

        !We need to put all our data into the variables we defined now
        !All output data is containted within the derived type "trajectory"
        ierr = nf90_put_var(file_id, variable_ids(1), particle_traj%E_x)
        PRINT *, ierr

        ierr = nf90_put_var(file_id, variable_ids(2), particle_traj%E_y)
        PRINT *, ierr

        ierr = nf90_put_var(file_id, variable_ids(3), particle_traj%x_traj)
        PRINT *, ierr
        
        ierr = nf90_put_var(file_id, variable_ids(4), particle_traj%y_traj)
        PRINT *, ierr

        ierr = nf90_put_var(file_id, variable_ids(5), particle_traj%vx_traj)
        PRINT *, ierr
        
        ierr = nf90_put_var(file_id, variable_ids(6), particle_traj%vy_traj)
        PRINT *, ierr

        ierr = nf90_put_var(file_id, variable_ids(7), particle_traj%ax_traj)
        PRINT *, ierr
        
        ierr = nf90_put_var(file_id, variable_ids(8), particle_traj%ay_traj)
        PRINT *, ierr

        ierr = nf90_put_var(file_id, variable_ids(9), rho_corrected)
        PRINT *, ierr

        ierr = nf90_put_var(file_id, variable_ids(10), phi_corrected)
        PRINT *, ierr

        ierr = nf90_close(file_id)
        PRINT *, ierr


    END SUBROUTINE WRITER

END MODULE NETCDF_WRITER