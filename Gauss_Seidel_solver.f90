MODULE grid_initialisation
    USE ISO_FORTRAN_ENV
    USE domain_tools
    IMPLICIT NONE

    CONTAINS
    
    ! Based on the sizes nx and ny and the type of problem,
    ! the subroutine creates the correct charge distribution 
    ! rho. The charge distribution values are defined as being in the 
    ! centre of the grid squares, as in the notes provided.

    SUBROUTINE define_rho(rho, nx,ny,problem)
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: rho
        REAL(REAL64),DIMENSION(:),ALLOCATABLE :: x_axis,y_axis

        REAL(REAL64),DIMENSION(2) :: axis_range
        
        REAL(REAL64) :: dx,dy
        INTEGER(INT64) :: x_dim,y_dim, i,j
        
        
        INTEGER, INTENT(IN) :: nx,ny
        CHARACTER(LEN=10), INTENT(IN) :: problem
        
        REAL(REAL64), DIMENSION(2) :: pos
        REAL(REAL64) :: x,y,mean_x,mean_y,std



        axis_range(1) = -1.0
        axis_range(2) = 1.0
        
        
        
        CALL create_axis(x_axis, nx, axis_range,1,dx)
        CALL create_axis(y_axis,ny,axis_range,1,dy)

        ! Gets the dimensions of axes. Unnecessary since 
        ! it's directly defined by nx, ny and the fact the 
        ! number of ghost cells on each side of the grid is 1.
        x_dim = INT(SIZE(x_axis))
        y_dim = INT(SIZE(y_axis))
        
       
        ! From above, could equivelantly have the limits as 
        ! 0:nx + 1 and similarly for y
        ALLOCATE(rho(0:x_dim -1, 0 : y_dim-1))

        ! TODO: Trim the input to ensure that it matches.
        ! Check correctness of user input.
        

        ! Sets up the charge distribution based on the pre-defined problems.
        
        IF (problem == 'null') THEN 
            rho =0.0
        ELSE IF (problem == 'single') THEN
            ! Sets the paramters of the problem.
            mean_x = 0.0
            mean_y = 0.0
            std = 0.1

            ! Loops over the grid and sets the correct 
            ! charge distribution. Omitts the ghost cells. i.e
            ! 0 and nx + 1 (similarly for y).
            DO i=1,nx
                DO j=1,ny
                    pos = position_converter(x_axis,y_axis,i,j)
                    x = pos(1)
                    y = pos(2)
    
                    rho(i,j) = Gaussian(x,y,mean_x,mean_y,std)
                END DO
            END DO
        END IF 

                    
        
    END SUBROUTINE
    
    ! Converts the grid (i,j) coordinates to the physical 
    ! (x,y) coordinates, given the x and y axes. i runs from 
    ! (1,nx) and j runs from (1,ny). Returns (x,y)
    
    FUNCTION position_converter(x_axis,y_axis,i,j)
        INTEGER(INT64), INTENT(IN) :: i,j
        REAL(REAL64),DIMENSION(:),ALLOCATABLE,INTENT(IN) :: x_axis,y_axis

        REAL(REAL64) :: x,y

        REAL(REAL64),DIMENSION(2) :: position_converter


        x = x_axis(i)
        y = y_axis(j)

        position_converter(1) = x
        position_converter(2) = y

    END FUNCTION

    ! Returns the value of a 2D Gaussian at a given point (x,y),
    ! with respective x and y means and a common standard deviation,
    ! to ensure that each indidividual charge distribution is 
    ! spherically symmetric.
    
    FUNCTION Gaussian(x,y,mean_x,mean_y,std)
        REAL(REAL64) :: Gaussian
        REAL(REAL64), INTENT(IN) :: x,y,mean_x,mean_y,std

        Gaussian = EXP(-( (x + mean_x)/std ) ** 2 - ( (y + mean_y)/std ) ** 2 )
    END FUNCTION



END MODULE



PROGRAM Solver
    USE ISO_FORTRAN_ENV
    USE domain_tools
    USE grid_initialisation
    USE command_line

    IMPLICIT NONE

    REAL(REAL64),DIMENSION(:,:), ALLOCATABLE :: rho

    INTEGER :: nx,ny
    CHARACTER(LEN=10) :: problem
    LOGICAL :: arg

    ! TODO: Check user input. Allocate exact string length for problem.
    
    CALL  parse_args
    

    ! Gets the requierd user inputs.
    arg = get_arg('nx',nx)
    arg = get_arg('ny',ny)
    arg = get_arg('problem',problem)


    CALL define_rho(rho,nx,ny,problem)
    
    ! Prints the charge distribution, depending on the user 
    ! input 'problem'
    
    PRINT*,rho


END PROGRAM


        
          
   
