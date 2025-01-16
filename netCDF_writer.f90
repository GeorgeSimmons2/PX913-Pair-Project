MODULE NETCDF_WRITER

    USE ISO_FORTRAN_ENV
    USE VERLET_MOD
    USE NETCDF

    IMPLICIT NONE
    
    TYPE :: DIMENSIONS
        INTEGER(INT64)   :: n_x, n_y, steps
        CHARACTER(LEN=30) :: n_x_name, n_y_name, steps_name !These are just to give appropriate names to each dimension
    END TYPE

    CONTAINS
    
    SUBROUTINE WRITER(particle_traj, rho, phi, filename, dims, ierr, nx, ny, problem, r_init, v_init)
        TYPE(DIMENSIONS), INTENT(INOUT)                 :: dims
        TYPE(TRAJECTORY), INTENT(INOUT)                 :: particle_traj !Easiest to use derived type for all trajectories
        INTEGER(INT32)                                  :: file_id
        INTEGER(INT32), DIMENSION(:), ALLOCATABLE       :: dimension_ids, variable_ids
        REAL(REAL64), DIMENSION(2), INTENT(IN)          :: r_init, v_init
        REAL(REAL64), DIMENSION(:, :), ALLOCATABLE      :: rho, phi, rho_corrected, phi_corrected
        INTENT(INOUT)                                   :: rho, phi
        INTEGER, INTENT(OUT)                     :: ierr
        CHARACTER(LEN=*), INTENT(IN)                    :: filename, problem
        INTEGER(INT64)                                  :: i, j, nx,ny

        ALLOCATE(rho_corrected(ny, nx))
        ALLOCATE(phi_corrected(ny, nx))
  
        DO i = 2, ny-1
            DO j = 2, nx-1
                rho_corrected(i, j) = rho(i, j)
                phi_corrected(i, j) = phi(i, j)
            END DO
        END DO
    
        !We are allocating for storing each id number of the variables and dimensions we have
        ALLOCATE(dimension_ids(4))
        ALLOCATE(variable_ids(12))

        ierr = nf90_create(filename, NF90_CLOBBER, file_id)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        !Dimensions n_x_name and n_y_name are for the E-field variables
        ierr = nf90_def_dim(file_id, dims%n_x_name, INT(dims%n_x, kind=4), dimension_ids(2))
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_def_dim(file_id, dims%n_y_name, INT(dims%n_y, kind=4), dimension_ids(1))
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        !Dimension steps_name are for the trajectory variables
        ierr = nf90_def_dim(file_id, dims%steps_name, INT(dims%steps + 1, kind=4), dimension_ids(3))
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        !Dimension steps_name are for the trajectory variables
        ierr = nf90_def_dim(file_id, "Initial Conditions", INT(2, kind=4), dimension_ids(4))
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        !Defining all our variables to store all trajectories, E-fields, charge densities and scalar fields
        ierr = nf90_def_var(file_id, particle_traj%Ex_name, NF90_REAL, dimension_ids(1:2), variable_ids(1))
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_def_var(file_id, particle_traj%Ey_name, NF90_REAL, dimension_ids(1:2), variable_ids(2))
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_def_var(file_id, particle_traj%x_name, NF90_REAL, dimension_ids(3), variable_ids(3))
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_def_var(file_id, particle_traj%y_name, NF90_REAL, dimension_ids(3), variable_ids(4))
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_def_var(file_id, particle_traj%vx_name, NF90_REAL, dimension_ids(3), variable_ids(5))
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_def_var(file_id, particle_traj%vy_name, NF90_REAL, dimension_ids(3), variable_ids(6))
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_def_var(file_id, particle_traj%ax_name, NF90_REAL, dimension_ids(3), variable_ids(7))
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_def_var(file_id, particle_traj%ay_name, NF90_REAL, dimension_ids(3), variable_ids(8))
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_def_var(file_id, particle_traj%rho_name, NF90_REAL, dimension_ids(1:2), variable_ids(9))
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_def_var(file_id, particle_traj%phi_name, NF90_REAL, dimension_ids(1:2), variable_ids(10))
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_def_var(file_id, "Initial Position", NF90_REAL, dimension_ids(4), variable_ids(11))
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_def_var(file_id, "Initial Velocity", NF90_REAL, dimension_ids(4), variable_ids(12))
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        !Meta-data about the problem described by global attribute and two above variables
        ierr = nf90_put_att(file_id, NF90_GLOBAL, "Problem", problem)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_enddef(file_id)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        !We need to put all our data into the variables we defined now
        !All output data is containted within the derived type "trajectory"
        ierr = nf90_put_var(file_id, variable_ids(1), particle_traj%E_x)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_put_var(file_id, variable_ids(2), particle_traj%E_y)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_put_var(file_id, variable_ids(3), particle_traj%x_traj)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
        
        ierr = nf90_put_var(file_id, variable_ids(4), particle_traj%y_traj)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_put_var(file_id, variable_ids(5), particle_traj%vx_traj)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
        
        ierr = nf90_put_var(file_id, variable_ids(6), particle_traj%vy_traj)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_put_var(file_id, variable_ids(7), particle_traj%ax_traj)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
        
        ierr = nf90_put_var(file_id, variable_ids(8), particle_traj%ay_traj)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_put_var(file_id, variable_ids(9), rho_corrected)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_put_var(file_id, variable_ids(10), phi_corrected)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_put_var(file_id, variable_ids(11), r_init)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_put_var(file_id, variable_ids(12), v_init)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_close(file_id)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF


    END SUBROUTINE WRITER

END MODULE NETCDF_WRITER