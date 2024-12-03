MODULE grid_initialisation
    USE ISO_FORTRAN_ENV
    USE domain_tools
    IMPLICIT NONE

    CONTAINS
    
    ! Based on the sizes nx and ny and the type of problem,
    ! the subroutine creates the correct charge distribution 
    ! rho. The charge distribution values are defined as being in the 
    ! centre of the grid squares, as in the notes provided.

    SUBROUTINE define_rho(rho,nx,ny,problem,dx,dy)
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: rho
        REAL(REAL64),DIMENSION(:),ALLOCATABLE :: x_axis,y_axis

        REAL(REAL64),DIMENSION(2) :: axis_range
        
        REAL(REAL64),INTENT(OUT) :: dx,dy
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
            ! 0 and nx + 1 (similarly for y). Loops from the top to bottom to 
            ! correctly, such that it uses the same convention as in the notes. 
            ! i.e the origin of the grid is in the bottom left, as per te he top
            ! left


            ! TODO: Set the input of the Gaussian using literals maybe.
            ! TODO: Make sure that the looping is correct.
            DO j=ny,1,-1
                DO i=1,nx
                    pos = position_converter(x_axis,y_axis,i,j)
                    x = pos(1)
                    y = pos(2)
    
                    rho(i,j) = Gaussian(x,y,mean_x,mean_y,std)
                END DO
            END DO
        
        ELSE IF (problem == 'double') THEN
            ! Uses the same looping convention as above. 
            ! This time the paramaters are set directly in the function call 
            ! to stop from creating more variables.
 
              DO j=ny,1,-1
                DO i=1,nx
                    pos = position_converter(x_axis,y_axis,i,j)
                    x = pos(1)
                    y = pos(2)
                    rho(i,j) = Gaussian(x,y,mean_x = 0.25_REAL64,mean_y = 0.25_REAL64,std = 0.1_REAL64) & 
                     + Gaussian(x,y,mean_x = -0.75_REAL64, mean_y = -0.75_REAL64, std = 0.2_REAL64)
                END DO
            END DO
        END IF 

    END SUBROUTINE
    
    SUBROUTINE calc_potential(rho,phi,dx,dy)
        REAL(REAL64),DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: rho
        REAL(REAL64),DIMENSION(:,:), ALLOCATABLE,INTENT(INOUT) :: phi
        
        REAL(REAL64), INTENT(IN) :: dx,dy
        INTEGER(INT64) :: i,j ,nx,ny,counter
        
        REAL(REAL64) :: error,factor,dphi_x,dphi_y

        LOGICAL :: loop_flag

        loop_flag = .TRUE.
        
        error = 0.0
        
        ! Need the -2 to account for the 2 ghost cells on each side of the grid.
        
       
        nx = SIZE(rho(:,1)) -2
        ny = SIZE(rho(1,:)) -2 

        ALLOCATE(phi(0:nx + 1,0:ny + 1))
        
        
        phi(0:nx +1,0:ny + 1) = 0.0_REAL64



        !phi(1:nx,1:ny)  = 1.0_REAL64
        
        counter = 0
        DO WHILE (loop_flag)
           
              DO j=ny,1,-1
                DO i=1,nx
                    dphi_x = (phi(i+1,j) + phi(i-1,j))/(dx**2)
                    dphi_y = (phi(i,j+1) + phi(i,j-1))/(dy**2)
                    
                    phi(i,j) = (-rho(i,j) + dphi_x + dphi_y)/( 2.0_REAL64/(dx**2) + 2.0_REAL64/(dy**2))

                END DO
            END DO
           
            error = calc_error(rho,phi,dx,dy)
            
            counter = counter + 1
            !PRINT*,phi
            !PRINT*,error
            IF (error < 10E-5) THEN
               
                loop_flag = .FALSE.
                PRINT '(A," ",G0," ",A)',' The algorithm took', counter, 'iterations to converge'
            END IF
           
            
        END DO
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

    FUNCTION calc_error(rho,phi,dx,dy)
        REAL(REAL64) :: calc_error
        REAL(REAL64),DIMENSION(:,:),ALLOCATABLE, INTENT(IN) :: rho,phi
        INTEGER(INT64) :: i,j,nx,ny
        REAL(REAL64), INTENT(IN) :: dx,dy
        REAL(REAL64) :: error_tot, error_rms




        nx = SIZE(rho(:,1)) -2
        ny = SIZE(rho(1,:)) -2
        
        error_tot = 0.0
        error_rms = 0.0
      

        DO j=ny,1,-1
            DO i=1,nx
                error_tot = error_tot + ABS( -rho(i,j) & 
                + (phi(i+1,j) -2.0_REAL64 * phi(i,j) + phi(i-1,j))/(dx**2)  & 
                + (phi(i,j+1) - 2.0_REAL64 * phi(i,j)  + phi(i,j-1))/(dy**2))
                
            END DO
        END DO
   
        DO j = ny,1,-1
            DO i=1,nx
                error_rms = error_rms + ( (phi(i+1,j) -2.0_REAL64 * phi(i,j) + phi(i-1,j))/(dx**2) &
                + (phi(i,j+1) - 2.0_REAL64 * phi(i,j)  + phi(i,j-1))/(dy**2))**2
             
            END DO
        END DO
        

        ! Either error, may be 0. This is the case when there is no charge denisty,
        ! and the second derivatives of the potential dissappear. However, there still is 
        ! a potential. 
      

        IF (error_rms <10E-14) THEN 
            error_rms = 1.0
        END IF


        error_rms = SQRT(error_rms/(nx * ny))
     
        calc_error = error_tot/error_rms
    END FUNCTION
END MODULE



PROGRAM Solver
    USE ISO_FORTRAN_ENV
    USE domain_tools
    USE grid_initialisation
    USE command_line

    IMPLICIT NONE

    REAL(REAL64),DIMENSION(:,:), ALLOCATABLE :: rho
    REAL(REAL64),DIMENSION(:,:), ALLOCATABLE :: phi
    INTEGER :: nx,ny,i,fu,j
    REAL(REAL64) :: dx,dy
    CHARACTER(LEN=10) :: problem
    LOGICAL :: arg

    CHARACTER(LEN=*),PARAMETER :: OUT_FILE = 'potential.csv'
    
        
    ! TODO: Check user input. Allocate exact string length for problem.
    
    CALL  parse_args
    

    ! Gets the requierd user inputs.
    arg = get_arg('nx',nx)
    arg = get_arg('ny',ny)
    arg = get_arg('problem',problem)

    fu = 50

    CALL define_rho(rho,nx,ny,problem,dx,dy)
    CALL calc_potential(rho,phi,dx,dy)
    
    ! Prints the charge distribution, depending on the user 
    ! input 'problem'
    
    !PRINT*,rho(4,:)
    
    OPEN(fu,action = 'write',file = OUT_FILE,status = 'replace')

    DO j=ny,1,-1
        DO i=1,nx
            WRITE(fu,'(F10.5)',ADVANCE = 'NO') phi(i,j)
            IF (i < nx) THEN
                WRITE(fu, '(A)', ADVANCE='NO') ','
            END IF
        END DO
        WRITE(fu,*)
    END DO
        
      

    CLOSE(fu)

END PROGRAM


        
          
   
