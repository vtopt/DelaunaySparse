PROGRAM SAMPLE_MAIN_S
! Driver code that reads a set P of data points from a file and computes
! the containing simplices and interpolation weights for a set Q of
! user-specified interpolation points using DELAUNAYSPARSES. If response
! values are provided, the interpolant f_{DT}(q) is also computed for all
! q \in Q.
!
! Usage: ./samples $(filepath)
!
! where $(filepath) is the relative or absolute path to the input file,
! formatted as follows.
!
! D,N,M,IR
! [Data/training points]
! [Response/function values]
! [Interpolation points]
!
! where
! D is the dimension of problem, 
! N is the number of data/training points (contained in lines 2 -- N+1),
! M is the number of interpolation points (contained in lines 2N+2 -- 2N+1+M),
! IR is the dimension of the output f(x) (the corresponding f(p) for p \in P
!   are stored in lines N+2 -- 2N+1).
!
! If IR = 0, then no interpolation will be done (and the M interpolation points
! are stored in lines N+2 -- N+1+M).
!
! A sample input file with D=2, N=43, M=101, IR=1 is provided by
! sample_input2d.dat.
! A sample input file with D=4, N=432, M=432, IR=1 is provided by
! sample_input4d.dat.
!
! Last Update: March, 2020
! Primary Author: Tyler Chang
USE DELSPARSE_MOD
USE OMP_LIB
IMPLICIT NONE

! Declare arguments and local data.
! Problem dimensions.
INTEGER :: D ! Problem dimension.
INTEGER :: N ! Number of data points.
INTEGER :: M ! Number of interpolation points.
INTEGER :: IR ! Response values (i.e., the dimension of the output).
! DELAUNAYSPARSE argument arrays.
REAL(KIND=R8), ALLOCATABLE :: PTS(:,:) ! The input data points.
REAL(KIND=R8), ALLOCATABLE :: Q(:,:) ! The interpolation points.
REAL(KIND=R8), ALLOCATABLE :: WEIGHTS(:,:) ! The interpolation weights.
INTEGER, ALLOCATABLE :: SIMPS(:,:) ! The indices of the simplex vertices.
INTEGER, ALLOCATABLE :: IERR(:) ! Array of integer error flags.
! Optional argument arrays.
REAL(KIND=R8), ALLOCATABLE :: INTERP_IN(:,:) ! Response value array.
REAL(KIND=R8), ALLOCATABLE :: INTERP_OUT(:,:) ! Output array for f_DT(q).
REAL(KIND=R8), ALLOCATABLE :: RNORM(:) ! Array of extrapolation residuals.
! Local variables.
INTEGER :: I ! Loop index/temp value.
REAL(KIND=R8) :: TICK ! The current clock time/total walltime.
CHARACTER(LEN=80) :: FILEPATH ! Input filepath.

! Open the file path $(filepath), and get the metadata from the
! first line (D, N, M, and IR).
CALL GET_COMMAND_ARGUMENT(1, FILEPATH)
OPEN(1, FILE=TRIM(FILEPATH))
READ(1, *) D, N, M, IR
IF(D .LE. 0 .OR. N .LE. 0 .OR. M .LE. 0) THEN
   WRITE(*,*) "Illegal input dimensions in input file, line 1."; STOP
END IF

! Allocate all necessarry arrays.
ALLOCATE(PTS(D,N), WEIGHTS(D+1,M), Q(D,M), SIMPS(D+1,M), IERR(M), &
  & RNORM(M), STAT=I)
IF(I .NE. 0) THEN
   WRITE(*,*) "Memory allocation error."; STOP
END IF

! Read the input data/training points into PTS.
DO I = 1, N
   READ(1, *) PTS(:, I)
END DO
! Check if there are any response values.
IF (IR > 0) THEN
   ! If so, allocate INTERP_IN and INTERP_OUT.
   ALLOCATE(INTERP_IN(IR,N), INTERP_OUT(IR,M), STAT=I)
   IF(I .NE. 0) THEN
      WRITE(*,*) "Memory allocation error."; STOP
   END IF
   ! Then, read the response values into INTERP_IN.
   DO I = 1, N
      READ(1, *) INTERP_IN(:,I)
   END DO
END IF
! Read the interpolation points into Q.
DO I = 1, M
   READ(1, *) Q(:, I)
END DO
CLOSE(1)

! Compute the interpolation results and time.
! If response values are provided, compute the outputs f_{DT}(q).
IF (IR > 0) THEN
   TICK = OMP_GET_WTIME()
   ! Call DELAUNAYSPARSES with INTERP_IN and INTERP_OUT.
   CALL DELAUNAYSPARSES(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, &
     ! Optional argument list.
     & INTERP_IN=INTERP_IN, INTERP_OUT=INTERP_OUT,             &
     & EPS=SQRT(EPSILON(0.0_R8)), EXTRAP=0.1_R8, RNORM=RNORM,  &
     & IBUDGET = 50000, CHAIN=.FALSE., EXACT=.TRUE.)
   TICK = OMP_GET_WTIME() - TICK
! Otherwise, just compute the simplices and weights.
ELSE
   TICK = OMP_GET_WTIME()
   ! Call DELAUNAYSPARSES without INTERP_IN and INTERP_OUT.
   CALL DELAUNAYSPARSES(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, &
     ! Optional argument list. Note that INTERP_IN and INTERP_OUT
     ! have been excluded.
     & EPS=SQRT(EPSILON(0.0_R8)), EXTRAP=0.1_R8, RNORM=RNORM,  &
     & IBUDGET = 50000, CHAIN=.FALSE., EXACT=.TRUE.)
   TICK = OMP_GET_WTIME() - TICK
   
END IF

! Display the results of the interpolation.
DO I = 1, M
   IF(IERR(I) .EQ. 0) THEN
      WRITE(*,10) 'Interpolation point: ', Q(:,I)
      WRITE(*,11) 'Simplex: ', SIMPS(:,I)
      WRITE(*,10) 'Weights: ', WEIGHTS(:,I)
      IF (IR > 0) THEN
         WRITE(*,12) 'f(x) = ', INTERP_OUT(:,I)
      END IF
   ELSE IF(IERR(I) .EQ. 1) THEN
      WRITE(*,10) 'Extrapolation point: ', Q(:,I)
      WRITE(*,11) 'Simplex: ', SIMPS(:,I)
      WRITE(*,10) 'Weights: ', WEIGHTS(:,I)
      IF (IR > 0) THEN
         WRITE(*,12) 'f(x) = ', INTERP_OUT(:,I)
      END IF
      WRITE(*,13) 'Residual: ', RNORM(I)
   ELSE IF(IERR(I) .EQ. 2) THEN
      WRITE(*,10) 'Extrapolation point: ', Q(:,I)
      WRITE(*,13) 'Residual: ', RNORM(I)
   ELSE
      WRITE(*,14) 'Error at point ', I, '. IERR(I) = ', IERR(I)
   END IF
END DO
! Print the timing data.
WRITE(*,15) M, ' points interpolated in ', TICK, ' seconds.'
10 FORMAT(1X,A,/,(1X,5ES15.7))
11 FORMAT(1X,A,/,(10I7))
12 FORMAT(1X,A,4ES15.7,/,(1X,5ES15.7))
13 FORMAT(1X,A,ES16.8)
14 FORMAT(1X,A,I7,A,I2)
15 FORMAT(/,I7,A,ES16.8,A,/)

END PROGRAM SAMPLE_MAIN_S
