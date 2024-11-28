MODULE NETCDF_WRITER

    USE ISO_FORTRAN_ENV
    USE VERLET_MOD
    USE NETCDF

    IMPLICIT NONE
    
    CONTAINS
    
    SUBROUTINE WRITER(particle_traj, filename, ierr)
        TYPE(trajectory), INTENT(IN) :: particle_traj !Easiest to use derived type for all trajectories
        INTEGER(INT64) :: ierr, file_id
        INTENT(OUT) :: ierr
        CHARACTER(LEN=*), INTENT(IN) :: filename


        ierr = nf90_create(filename, NF90_CLOBBER, file_id)
        PRINT *, ierr

        !Insert all information about file here
        ierr = nf90_def_dim(file_id, )

        ierr = nf90_close(file_id)
        PRINT *, ierr


    END SUBROUTINE WRITER

END MODULE NETCDF_WRITER