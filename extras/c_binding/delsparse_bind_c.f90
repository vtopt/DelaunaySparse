

SUBROUTINE C_DELAUNAYSPARSES_NOOPTS(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR) &
           BIND(C, NAME="c_delaunaysparses")
   ! This is a wrapper for DELAUNAYSPARSES with no optional arguments.
   !
   !
   ! On input:
   !
   ! D is the dimension of the space for PTS and Q.
   !
   ! N is the number of data points in PTS.
   !
   ! PTS(1:D,1:N) is a real valued matrix with N columns, each containing the
   !    coordinates of a single data point in R^D.
   !
   ! M is the number of interpolation points in Q.
   !
   ! Q(1:D,1:M) is a real valued matrix with M columns, each containing the
   !    coordinates of a single interpolation point in R^D.
   !
   !
   ! On output:
   !
   ! PTS and Q have been rescaled and shifted. All the data points in PTS
   !    are now contained in the unit hyperball in R^D, and the points in Q
   !    have been shifted and scaled accordingly in relation to PTS.
   !
   ! SIMPS(1:D+1,1:M) contains the D+1 integer indices (corresponding to columns
   !    in PTS) for the D+1 vertices of the Delaunay simplex containing each
   !    interpolation point in Q.
   !
   ! WEIGHTS(1:D+1,1:M) contains the D+1 real valued weights for expressing each
   !    point in Q as a convex combination of the D+1 corresponding vertices
   !    in SIMPS.
   !
   ! IERR(1:M) contains integer valued error flags associated with the
   !    computation of each of the M interpolation points in Q. The error
   !    codes are given in the definition of DELAUNAYSPARSES in delsparse.f90.
   !
   !
   ! LAST UPDATE:
   !   11/2020 by THC
   !
   USE REAL_PRECISION , ONLY : R8
   USE ISO_C_BINDING

   IMPLICIT NONE
  
   INTEGER(C_INT), INTENT(IN) :: D
   INTEGER(C_INT), INTENT(IN) :: N
   REAL(C_DOUBLE), INTENT(INOUT) :: PTS(D,N)
   INTEGER(C_INT), INTENT(IN) :: M
   REAL(C_DOUBLE), INTENT(INOUT) :: Q(D,M)
   INTEGER(C_INT), INTENT(OUT) :: SIMPS(D+1,M)
   REAL(C_DOUBLE), INTENT(OUT) :: WEIGHTS(D+1,M)
   INTEGER(C_INT), INTENT(OUT) :: IERR(M)
  
   INTERFACE
      SUBROUTINE DELAUNAYSPARSES(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, &
                                 INTERP_IN, INTERP_OUT, EPS, EXTRAP,    &
                                 RNORM, IBUDGET, CHAIN, EXACT)
         USE REAL_PRECISION , ONLY : R8
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: D
         INTEGER, INTENT(IN) :: N
         REAL(KIND=R8), INTENT(INOUT) :: PTS(:,:)
         INTEGER, INTENT(IN) :: M
         REAL(KIND=R8), INTENT(INOUT) :: Q(:,:)
         INTEGER, INTENT(OUT) :: SIMPS(:,:)
         REAL(KIND=R8), INTENT(OUT) :: WEIGHTS(:,:)
         INTEGER, INTENT(OUT) :: IERR(:)
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: INTERP_IN(:,:)
         REAL(KIND=R8), INTENT(OUT), OPTIONAL :: INTERP_OUT(:,:)
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: EPS
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: EXTRAP
         REAL(KIND=R8), INTENT(OUT), OPTIONAL :: RNORM(:)
         INTEGER, INTENT(IN), OPTIONAL :: IBUDGET
         LOGICAL, INTENT(IN), OPTIONAL :: CHAIN
         LOGICAL, INTENT(IN), OPTIONAL :: EXACT
      END SUBROUTINE DELAUNAYSPARSES
   END INTERFACE

   INTEGER :: D_LOC
   INTEGER :: N_LOC
   REAL(KIND=R8) :: PTS_LOC(D, N)
   INTEGER :: M_LOC
   REAL(KIND=R8) :: Q_LOC(D, M)
   INTEGER :: SIMPS_LOC(D+1, M)
   REAL(KIND=R8) :: WEIGHTS_LOC(D+1, M)
   INTEGER :: IERR_LOC(M)

   D_LOC = INT(D)
   N_LOC = INT(N)
   PTS_LOC = REAL(PTS, KIND=R8)
   M_LOC = INT(M)
   Q_LOC = REAL(Q, KIND=R8)

   CALL DELAUNAYSPARSES(D_LOC, N_LOC, PTS_LOC, M_LOC, Q_LOC, SIMPS_LOC, &
                        WEIGHTS_LOC, IERR_LOC)

   PTS = REAL(PTS_LOC, KIND=C_DOUBLE)
   Q = REAL(Q_LOC, KIND=C_DOUBLE)
   SIMPS = INT(SIMPS_LOC, KIND=C_INT)
   WEIGHTS = REAL(WEIGHTS_LOC, KIND=C_DOUBLE)
   IERR = INT(IERR_LOC, KIND=C_INT)

   RETURN
END SUBROUTINE C_DELAUNAYSPARSES_NOOPTS


SUBROUTINE C_DELAUNAYSPARSES_INTERP(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, &
                                    IR, INTERP_IN, INTERP_OUT)             &
           BIND(C, NAME="c_delaunaysparses_interp")
   ! This is a wrapper for DELAUNAYSPARSES with INTERP_IN and INTERP_OUT
   ! specified, but no other optional arguments. Unlike the Fortran interface,
   ! in this interface the dimension of the response variables (IR) must
   ! be explicitly specified by an additional input, IR.
   !
   !
   ! On input:
   !
   ! D is the dimension of the space for PTS and Q.
   !
   ! N is the number of data points in PTS.
   !
   ! PTS(1:D,1:N) is a real valued matrix with N columns, each containing the
   !    coordinates of a single data point in R^D.
   !
   ! M is the number of interpolation points in Q.
   !
   ! Q(1:D,1:M) is a real valued matrix with M columns, each containing the
   !    coordinates of a single interpolation point in R^D.
   !
   ! IR is the dimension of the response variables.
   !
   ! INTERP_IN(1:IR,1:N) contains real valued response vectors for each of
   !    the data points in PTS on input. The first dimension of INTERP_IN is
   !    inferred to be the dimension of these response vectors, and the
   !    second dimension must match N.
   !
   !
   ! On output:
   !
   ! PTS and Q have been rescaled and shifted. All the data points in PTS
   !    are now contained in the unit hyperball in R^D, and the points in Q
   !    have been shifted and scaled accordingly in relation to PTS.
   !
   ! SIMPS(1:D+1,1:M) contains the D+1 integer indices (corresponding to columns
   !    in PTS) for the D+1 vertices of the Delaunay simplex containing each
   !    interpolation point in Q.
   !
   ! WEIGHTS(1:D+1,1:M) contains the D+1 real valued weights for expressing each
   !    point in Q as a convex combination of the D+1 corresponding vertices
   !    in SIMPS.
   !
   ! IERR(1:M) contains integer valued error flags associated with the
   !    computation of each of the M interpolation points in Q. The error
   !    codes are given in the definition of DELAUNAYSPARSES in delsparse.f90.
   ! 
   ! INTERP_OUT(1:IR,1:M) contains real valued response vectors for each
   !    interpolation point in Q on output. The first dimension of INTERP_OU
   !    must match the first dimension of INTERP_IN, and the second dimension
   !    must match M.
   !
   !
   ! LAST UPDATE:
   !   11/2020 by THC
   !
   USE REAL_PRECISION , ONLY : R8
   USE ISO_C_BINDING

   IMPLICIT NONE
  
   INTEGER(C_INT), INTENT(IN) :: D
   INTEGER(C_INT), INTENT(IN) :: N
   REAL(C_DOUBLE), INTENT(INOUT) :: PTS(D,N)
   INTEGER(C_INT), INTENT(IN) :: M
   REAL(C_DOUBLE), INTENT(INOUT) :: Q(D,M)
   INTEGER(C_INT), INTENT(OUT) :: SIMPS(D+1,M)
   REAL(C_DOUBLE), INTENT(OUT) :: WEIGHTS(D+1,M)
   INTEGER(C_INT), INTENT(OUT) :: IERR(M)
   INTEGER(C_INT), INTENT(IN) :: IR
   REAL(C_DOUBLE), INTENT(IN) :: INTERP_IN(IR, N)
   REAL(C_DOUBLE), INTENT(OUT) :: INTERP_OUT(IR, M)
  
   INTERFACE
      SUBROUTINE DELAUNAYSPARSES(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, &
                                 INTERP_IN, INTERP_OUT, EPS, EXTRAP,    &
                                 RNORM, IBUDGET, CHAIN, EXACT)
         USE REAL_PRECISION , ONLY : R8
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: D
         INTEGER, INTENT(IN) :: N
         REAL(KIND=R8), INTENT(INOUT) :: PTS(:,:)
         INTEGER, INTENT(IN) :: M
         REAL(KIND=R8), INTENT(INOUT) :: Q(:,:)
         INTEGER, INTENT(OUT) :: SIMPS(:,:)
         REAL(KIND=R8), INTENT(OUT) :: WEIGHTS(:,:)
         INTEGER, INTENT(OUT) :: IERR(:)
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: INTERP_IN(:,:)
         REAL(KIND=R8), INTENT(OUT), OPTIONAL :: INTERP_OUT(:,:)
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: EPS
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: EXTRAP
         REAL(KIND=R8), INTENT(OUT), OPTIONAL :: RNORM(:)
         INTEGER, INTENT(IN), OPTIONAL :: IBUDGET
         LOGICAL, INTENT(IN), OPTIONAL :: CHAIN
         LOGICAL, INTENT(IN), OPTIONAL :: EXACT
      END SUBROUTINE DELAUNAYSPARSES
   END INTERFACE

   INTEGER :: D_LOC
   INTEGER :: N_LOC
   REAL(KIND=R8) :: PTS_LOC(D, N)
   INTEGER :: M_LOC
   REAL(KIND=R8) :: Q_LOC(D, M)
   INTEGER :: SIMPS_LOC(D+1, M)
   REAL(KIND=R8) :: WEIGHTS_LOC(D+1, M)
   INTEGER :: IERR_LOC(M)
   REAL(KIND=R8) :: INTERP_IN_LOC(IR, N)
   REAL(KIND=R8) :: INTERP_OUT_LOC(IR, M)

   D_LOC = INT(D)
   N_LOC = INT(N)
   PTS_LOC = REAL(PTS, KIND=R8)
   M_LOC = INT(M)
   Q_LOC = REAL(Q, KIND=R8)
   INTERP_IN_LOC = REAL(INTERP_IN, KIND=R8)

   CALL DELAUNAYSPARSES(D_LOC, N_LOC, PTS_LOC, M_LOC, Q_LOC, SIMPS_LOC, &
                        WEIGHTS_LOC, IERR_LOC, INTERP_IN=INTERP_IN_LOC, &
                        INTERP_OUT=INTERP_OUT_LOC)

   PTS = REAL(PTS_LOC, KIND=C_DOUBLE)
   Q = REAL(Q_LOC, KIND=C_DOUBLE)
   SIMPS = INT(SIMPS_LOC, KIND=C_INT)
   WEIGHTS = REAL(WEIGHTS_LOC, KIND=C_DOUBLE)
   IERR = INT(IERR_LOC, KIND=C_INT)
   INTERP_OUT = REAL(INTERP_OUT_LOC, KIND=C_DOUBLE)

   RETURN
END SUBROUTINE C_DELAUNAYSPARSES_INTERP


SUBROUTINE C_DELAUNAYSPARSES_OPTS(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, EPS, &
                                  EXTRAP, RNORM, IBUDGET, CHAIN, EXACT)       &
           BIND(C, NAME="c_delaunaysparses_opts")
   ! This is a wrapper for DELAUNAYSPARSES without INTERP_IN and INTERP_OUT,
   ! but all other optional arguments present.
   !
   !
   ! On input:
   !
   ! D is the dimension of the space for PTS and Q.
   !
   ! N is the number of data points in PTS.
   !
   ! PTS(1:D,1:N) is a real valued matrix with N columns, each containing the
   !    coordinates of a single data point in R^D.
   !
   ! M is the number of interpolation points in Q.
   !
   ! Q(1:D,1:M) is a real valued matrix with M columns, each containing the
   !    coordinates of a single interpolation point in R^D.
   !
   ! EXTRAP contains the real maximum extrapolation distance (relative to the
   !    diameter of PTS) on input. Interpolation at a point outside the convex
   !    hull of PTS is done by projecting that point onto the convex hull, and
   !    then doing normal Delaunay interpolation at that projection.
   !    Interpolation at any point in Q that is more than EXTRAP * DIAMETER(PTS)
   !    units outside the convex hull of PTS will not be done and an error code
   !    of 2 will be returned. Note that computing the projection can be
   !    expensive. Setting EXTRAP=0 will cause all extrapolation points to be
   !    ignored without ever computing a projection.
   ! 
   ! IBUDGET on input contains the integer budget for performing flips while
   !    iterating toward the simplex containing each interpolation point in Q.
   !    This prevents DELAUNAYSPARSES from falling into an infinite loop when
   !    an inappropriate value of EPS is given with respect to the problem
   !    conditioning. For most cases, the default value of 50000 should be
   !    more than sufficient.
   !
   ! CHAIN is a logical input argument that determines whether a new first
   !    simplex should be constructed for each interpolation point
   !    (CHAIN=.FALSE.), or whether the simplex walks should be "daisy-chained."
   !    Setting CHAIN=.TRUE. is generally not recommended, unless the size of
   !    the triangulation is relatively small or the interpolation points are
   !    known to be tightly clustered.
   !
   ! EXACT is a logical input argument that determines whether the exact
   !    diameter should be computed and whether a check for duplicate data
   !    points should be performed in advance. When EXACT=.FALSE., the
   !    diameter of PTS is approximated by twice the distance from the
   !    barycenter of PTS to the farthest point in PTS, and no check is
   !    done to find the closest pair of points, which could result in hard
   !    to find bugs later on. When EXACT=.TRUE., the exact diameter is
   !    computed and an error is returned whenever PTS contains duplicate
   !    values up to the precision EPS. Setting EXACT=.FALSE. could result
   !    in significant speedup when N is large, but it is strongly
   !    recommended that most users leave EXACT=.TRUE., as setting
   !    EXACT=.FALSE. could result in input errors that are difficult
   !    to identify. Also, the diameter approximation could be wrong by up
   !    to a factor of two.
   ! 
   !
   ! On output:
   !
   ! PTS and Q have been rescaled and shifted. All the data points in PTS
   !    are now contained in the unit hyperball in R^D, and the points in Q
   !    have been shifted and scaled accordingly in relation to PTS.
   !
   ! SIMPS(1:D+1,1:M) contains the D+1 integer indices (corresponding to columns
   !    in PTS) for the D+1 vertices of the Delaunay simplex containing each
   !    interpolation point in Q.
   !
   ! WEIGHTS(1:D+1,1:M) contains the D+1 real valued weights for expressing each
   !    point in Q as a convex combination of the D+1 corresponding vertices
   !    in SIMPS.
   !
   ! IERR(1:M) contains integer valued error flags associated with the
   !    computation of each of the M interpolation points in Q. The error
   !    codes are given in the definition of DELAUNAYSPARSES in delsparse.f90.
   !
   ! RNORM(1:M) contains the real unscaled projection (2-norm) distances from
   !    any projection computations on output.
   !
   !
   ! LAST UPDATE:
   !   11/2020 by THC
   !
   USE REAL_PRECISION , ONLY : R8
   USE ISO_C_BINDING

   IMPLICIT NONE
  
   INTEGER(C_INT), INTENT(IN) :: D
   INTEGER(C_INT), INTENT(IN) :: N
   REAL(C_DOUBLE), INTENT(INOUT) :: PTS(D,N)
   INTEGER(C_INT), INTENT(IN) :: M
   REAL(C_DOUBLE), INTENT(INOUT) :: Q(D,M)
   INTEGER(C_INT), INTENT(OUT) :: SIMPS(D+1,M)
   REAL(C_DOUBLE), INTENT(OUT) :: WEIGHTS(D+1,M)
   INTEGER(C_INT), INTENT(OUT) :: IERR(M)
   REAL(C_DOUBLE), INTENT(IN) :: EPS
   REAL(C_DOUBLE), INTENT(IN) :: EXTRAP
   REAL(C_DOUBLE), INTENT(OUT) :: RNORM(M)
   INTEGER(C_INT), INTENT(IN) :: IBUDGET
   LOGICAL(C_BOOL), INTENT(IN) :: CHAIN
   LOGICAL(C_BOOL), INTENT(IN) :: EXACT
  
   INTERFACE
      SUBROUTINE DELAUNAYSPARSES(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, &
                                 INTERP_IN, INTERP_OUT, EPS, EXTRAP,    &
                                 RNORM, IBUDGET, CHAIN, EXACT)
         USE REAL_PRECISION , ONLY : R8
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: D
         INTEGER, INTENT(IN) :: N
         REAL(KIND=R8), INTENT(INOUT) :: PTS(:,:)
         INTEGER, INTENT(IN) :: M
         REAL(KIND=R8), INTENT(INOUT) :: Q(:,:)
         INTEGER, INTENT(OUT) :: SIMPS(:,:)
         REAL(KIND=R8), INTENT(OUT) :: WEIGHTS(:,:)
         INTEGER, INTENT(OUT) :: IERR(:)
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: INTERP_IN(:,:)
         REAL(KIND=R8), INTENT(OUT), OPTIONAL :: INTERP_OUT(:,:)
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: EPS
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: EXTRAP
         REAL(KIND=R8), INTENT(OUT), OPTIONAL :: RNORM(:)
         INTEGER, INTENT(IN), OPTIONAL :: IBUDGET
         LOGICAL, INTENT(IN), OPTIONAL :: CHAIN
         LOGICAL, INTENT(IN), OPTIONAL :: EXACT
      END SUBROUTINE DELAUNAYSPARSES
   END INTERFACE

   INTEGER :: D_LOC
   INTEGER :: N_LOC
   REAL(KIND=R8) :: PTS_LOC(D, N)
   INTEGER :: M_LOC
   REAL(KIND=R8) :: Q_LOC(D, M)
   INTEGER :: SIMPS_LOC(D+1, M)
   REAL(KIND=R8) :: WEIGHTS_LOC(D+1, M)
   INTEGER :: IERR_LOC(M)
   REAL(KIND=R8) :: EPS_LOC
   REAL(KIND=R8) :: EXTRAP_LOC
   REAL(KIND=R8) :: RNORM_LOC(M)
   INTEGER :: IBUDGET_LOC
   LOGICAL :: CHAIN_LOC
   LOGICAL :: EXACT_LOC

   D_LOC = INT(D)
   N_LOC = INT(N)
   PTS_LOC = REAL(PTS, KIND=R8)
   M_LOC = INT(M)
   Q_LOC = REAL(Q, KIND=R8)
   EPS_LOC = REAL(EPS, KIND=R8)
   EXTRAP_LOC = REAL(EXTRAP, KIND=R8)
   IBUDGET_LOC = INT(IBUDGET)
   CHAIN_LOC = LOGICAL(CHAIN)
   EXACT_LOC = LOGICAL(EXACT)

   CALL DELAUNAYSPARSES(D_LOC, N_LOC, PTS_LOC, M_LOC, Q_LOC, SIMPS_LOC, &
                        WEIGHTS_LOC, IERR_LOC, EPS=EPS_LOC,             &
                        EXTRAP=EXTRAP_LOC, RNORM=RNORM_LOC,             &
                        IBUDGET=IBUDGET_LOC, CHAIN=CHAIN_LOC,           &
                        EXACT=EXACT_LOC)

   PTS = REAL(PTS_LOC, KIND=C_DOUBLE)
   Q = REAL(Q_LOC, KIND=C_DOUBLE)
   SIMPS = INT(SIMPS_LOC, KIND=C_INT)
   WEIGHTS = REAL(WEIGHTS_LOC, KIND=C_DOUBLE)
   IERR = INT(IERR_LOC, KIND=C_INT)
   RNORM = REAL(RNORM_LOC, KIND=C_DOUBLE)

   RETURN
END SUBROUTINE C_DELAUNAYSPARSES_OPTS


SUBROUTINE C_DELAUNAYSPARSES_INTERP_OPTS(D, N, PTS, M, Q, SIMPS, WEIGHTS,    &
                                         IERR, IR, INTERP_IN, INTERP_OUT,    &
                                         EPS, EXTRAP, RNORM, IBUDGET, CHAIN, &
                                         EXACT, PMODE)                       &
           BIND(C, NAME="c_delaunaysparses_interp_opts")
   ! This is a wrapper for DELAUNAYSPARSES with all optional arguments present.
   !
   !
   ! On input:
   !
   ! D is the dimension of the space for PTS and Q.
   !
   ! N is the number of data points in PTS.
   !
   ! PTS(1:D,1:N) is a real valued matrix with N columns, each containing the
   !    coordinates of a single data point in R^D.
   !
   ! M is the number of interpolation points in Q.
   !
   ! Q(1:D,1:M) is a real valued matrix with M columns, each containing the
   !    coordinates of a single interpolation point in R^D.
   !
   ! IR is the dimension of the response variables.
   !
   ! INTERP_IN(1:IR,1:N) contains real valued response vectors for each of
   !    the data points in PTS on input. The first dimension of INTERP_IN is
   !    inferred to be the dimension of these response vectors, and the
   !    second dimension must match N.
   !
   ! EXTRAP contains the real maximum extrapolation distance (relative to the
   !    diameter of PTS) on input. Interpolation at a point outside the convex
   !    hull of PTS is done by projecting that point onto the convex hull, and
   !    then doing normal Delaunay interpolation at that projection.
   !    Interpolation at any point in Q that is more than EXTRAP * DIAMETER(PTS)
   !    units outside the convex hull of PTS will not be done and an error code
   !    of 2 will be returned. Note that computing the projection can be
   !    expensive. Setting EXTRAP=0 will cause all extrapolation points to be
   !    ignored without ever computing a projection.
   ! 
   ! IBUDGET on input contains the integer budget for performing flips while
   !    iterating toward the simplex containing each interpolation point in Q.
   !    This prevents DELAUNAYSPARSES from falling into an infinite loop when
   !    an inappropriate value of EPS is given with respect to the problem
   !    conditioning. For most cases, the default value of 50000 should be
   !    more than sufficient.
   !
   ! CHAIN is a logical input argument that determines whether a new first
   !    simplex should be constructed for each interpolation point
   !    (CHAIN=.FALSE.), or whether the simplex walks should be "daisy-chained."
   !    Setting CHAIN=.TRUE. is generally not recommended, unless the size of
   !    the triangulation is relatively small or the interpolation points are
   !    known to be tightly clustered.
   !
   ! EXACT is a logical input argument that determines whether the exact
   !    diameter should be computed and whether a check for duplicate data
   !    points should be performed in advance. When EXACT=.FALSE., the
   !    diameter of PTS is approximated by twice the distance from the
   !    barycenter of PTS to the farthest point in PTS, and no check is
   !    done to find the closest pair of points, which could result in hard
   !    to find bugs later on. When EXACT=.TRUE., the exact diameter is
   !    computed and an error is returned whenever PTS contains duplicate
   !    values up to the precision EPS. Setting EXACT=.FALSE. could result
   !    in significant speedup when N is large, but it is strongly
   !    recommended that most users leave EXACT=.TRUE., as setting
   !    EXACT=.FALSE. could result in input errors that are difficult
   !    to identify. Also, the diameter approximation could be wrong by up
   !    to a factor of two.
   !
   ! PMODE is an integer specifying the level of parallelism to be exploited.
   !    If PMODE = 1, then parallelism is exploited at the level of the loop
   !    over all interpolation points (Level 1 parallelism).
   !    If PMODE = 2, then parallelism is exploited at the level of the loops
   !    over data points when constructing/flipping simplices (Level 2
   !    parallelism).
   !    If PMODE = 3, then parallelism is exploited at both levels. Note: this
   !    implies that the total number of threads active at any time could be up
   !    to OMP_NUM_THREADS^2.
   !
   !
   ! On output:
   !
   ! PTS and Q have been rescaled and shifted. All the data points in PTS
   !    are now contained in the unit hyperball in R^D, and the points in Q
   !    have been shifted and scaled accordingly in relation to PTS.
   !
   ! SIMPS(1:D+1,1:M) contains the D+1 integer indices (corresponding to columns
   !    in PTS) for the D+1 vertices of the Delaunay simplex containing each
   !    interpolation point in Q.
   !
   ! WEIGHTS(1:D+1,1:M) contains the D+1 real valued weights for expressing each
   !    point in Q as a convex combination of the D+1 corresponding vertices
   !    in SIMPS.
   !
   ! IERR(1:M) contains integer valued error flags associated with the
   !    computation of each of the M interpolation points in Q. The error
   !    codes are given in the definition of DELAUNAYSPARSES in delsparse.f90.
   ! 
   ! INTERP_OUT(1:IR,1:M) contains real valued response vectors for each
   !    interpolation point in Q on output. The first dimension of INTERP_OUT
   !    must match the first dimension of INTERP_IN, and the second dimension
   !    must match M.
   ! 
   ! RNORM(1:M) contains the real unscaled projection (2-norm) distances from
   !    any projection computations on output.
   !
   !
   ! LAST UPDATE:
   !   11/2020 by THC
   !
   USE REAL_PRECISION , ONLY : R8
   USE ISO_C_BINDING

   IMPLICIT NONE
  
   INTEGER(C_INT), INTENT(IN) :: D
   INTEGER(C_INT), INTENT(IN) :: N
   REAL(C_DOUBLE), INTENT(INOUT) :: PTS(D,N)
   INTEGER(C_INT), INTENT(IN) :: M
   REAL(C_DOUBLE), INTENT(INOUT) :: Q(D,M)
   INTEGER(C_INT), INTENT(OUT) :: SIMPS(D+1,M)
   REAL(C_DOUBLE), INTENT(OUT) :: WEIGHTS(D+1,M)
   INTEGER(C_INT), INTENT(OUT) :: IERR(M)
   INTEGER(C_INT), INTENT(IN) :: IR
   REAL(C_DOUBLE), INTENT(IN) :: INTERP_IN(IR, N)
   REAL(C_DOUBLE), INTENT(OUT) :: INTERP_OUT(IR, M)
   REAL(C_DOUBLE), INTENT(IN) :: EPS
   REAL(C_DOUBLE), INTENT(IN) :: EXTRAP
   REAL(C_DOUBLE), INTENT(OUT) :: RNORM(M)
   INTEGER(C_INT), INTENT(IN) :: IBUDGET
   LOGICAL(C_BOOL), INTENT(IN) :: CHAIN
   LOGICAL(C_BOOL), INTENT(IN) :: EXACT
   INTEGER(C_INT), INTENT(IN) :: PMODE
  
   INTERFACE
      SUBROUTINE DELAUNAYSPARSES(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, &
                                 INTERP_IN, INTERP_OUT, EPS, EXTRAP,    &
                                 RNORM, IBUDGET, CHAIN, EXACT)
         USE REAL_PRECISION , ONLY : R8
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: D
         INTEGER, INTENT(IN) :: N
         REAL(KIND=R8), INTENT(INOUT) :: PTS(:,:)
         INTEGER, INTENT(IN) :: M
         REAL(KIND=R8), INTENT(INOUT) :: Q(:,:)
         INTEGER, INTENT(OUT) :: SIMPS(:,:)
         REAL(KIND=R8), INTENT(OUT) :: WEIGHTS(:,:)
         INTEGER, INTENT(OUT) :: IERR(:)
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: INTERP_IN(:,:)
         REAL(KIND=R8), INTENT(OUT), OPTIONAL :: INTERP_OUT(:,:)
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: EPS
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: EXTRAP
         REAL(KIND=R8), INTENT(OUT), OPTIONAL :: RNORM(:)
         INTEGER, INTENT(IN), OPTIONAL :: IBUDGET
         LOGICAL, INTENT(IN), OPTIONAL :: CHAIN
         LOGICAL, INTENT(IN), OPTIONAL :: EXACT
      END SUBROUTINE DELAUNAYSPARSES
   END INTERFACE

   INTEGER :: D_LOC
   INTEGER :: N_LOC
   REAL(KIND=R8) :: PTS_LOC(D, N)
   INTEGER :: M_LOC
   REAL(KIND=R8) :: Q_LOC(D, M)
   INTEGER :: SIMPS_LOC(D+1, M)
   REAL(KIND=R8) :: WEIGHTS_LOC(D+1, M)
   INTEGER :: IERR_LOC(M)
   REAL(KIND=R8) :: INTERP_IN_LOC(IR, N)
   REAL(KIND=R8) :: INTERP_OUT_LOC(IR, M)
   REAL(KIND=R8) :: EPS_LOC
   REAL(KIND=R8) :: EXTRAP_LOC
   REAL(KIND=R8) :: RNORM_LOC(M)
   INTEGER :: IBUDGET_LOC
   LOGICAL :: CHAIN_LOC
   LOGICAL :: EXACT_LOC
   INTEGER :: PMODE_LOC

   D_LOC = INT(D)
   N_LOC = INT(N)
   PTS_LOC = REAL(PTS, KIND=R8)
   M_LOC = INT(M)
   Q_LOC = REAL(Q, KIND=R8)
   INTERP_IN_LOC = REAL(INTERP_IN, KIND=R8)
   EPS_LOC = REAL(EPS, KIND=R8)
   EXTRAP_LOC = REAL(EXTRAP, KIND=R8)
   IBUDGET_LOC = INT(IBUDGET)
   CHAIN_LOC = LOGICAL(CHAIN)
   EXACT_LOC = LOGICAL(EXACT)
   PMODE_LOC = INT(PMODE)

   CALL DELAUNAYSPARSES(D_LOC, N_LOC, PTS_LOC, M_LOC, Q_LOC, SIMPS_LOC, &
                        WEIGHTS_LOC, IERR_LOC, INTERP_IN=INTERP_IN_LOC, &
                        INTERP_OUT=INTERP_OUT_LOC, EPS=EPS_LOC,         &
                        EXTRAP=EXTRAP_LOC, RNORM=RNORM_LOC,             &
                        IBUDGET=IBUDGET_LOC, CHAIN=CHAIN_LOC,           &
                        EXACT=EXACT_LOC)

   PTS = REAL(PTS_LOC, KIND=C_DOUBLE)
   Q = REAL(Q_LOC, KIND=C_DOUBLE)
   SIMPS = INT(SIMPS_LOC, KIND=C_INT)
   WEIGHTS = REAL(WEIGHTS_LOC, KIND=C_DOUBLE)
   IERR = INT(IERR_LOC, KIND=C_INT)
   INTERP_OUT = REAL(INTERP_OUT_LOC, C_DOUBLE)
   RNORM = REAL(RNORM_LOC, KIND=C_DOUBLE)

   RETURN
END SUBROUTINE C_DELAUNAYSPARSES_INTERP_OPTS


SUBROUTINE C_DELAUNAYSPARSEP_NOOPTS(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR) &
           BIND(C, NAME="c_delaunaysparsep")
   ! This is a wrapper for DELAUNAYSPARSEP with no optional arguments.
   !
   !
   ! On input:
   !
   ! D is the dimension of the space for PTS and Q.
   !
   ! N is the number of data points in PTS.
   !
   ! PTS(1:D,1:N) is a real valued matrix with N columns, each containing the
   !    coordinates of a single data point in R^D.
   !
   ! M is the number of interpolation points in Q.
   !
   ! Q(1:D,1:M) is a real valued matrix with M columns, each containing the
   !    coordinates of a single interpolation point in R^D.
   !
   !
   ! On output:
   !
   ! PTS and Q have been rescaled and shifted. All the data points in PTS
   !    are now contained in the unit hyperball in R^D, and the points in Q
   !    have been shifted and scaled accordingly in relation to PTS.
   !
   ! SIMPS(1:D+1,1:M) contains the D+1 integer indices (corresponding to columns
   !    in PTS) for the D+1 vertices of the Delaunay simplex containing each
   !    interpolation point in Q.
   !
   ! WEIGHTS(1:D+1,1:M) contains the D+1 real valued weights for expressing each
   !    point in Q as a convex combination of the D+1 corresponding vertices
   !    in SIMPS.
   !
   ! IERR(1:M) contains integer valued error flags associated with the
   !    computation of each of the M interpolation points in Q. The error
   !    codes are given in the definition of DELAUNAYSPARSEP in delsparse.f90.
   !
   !
   ! LAST UPDATE:
   !   11/2020 by THC
   !
   USE REAL_PRECISION , ONLY : R8
   USE ISO_C_BINDING
   IMPLICIT NONE
  
   INTEGER(C_INT), INTENT(IN) :: D
   INTEGER(C_INT), INTENT(IN) :: N
   REAL(C_DOUBLE), INTENT(INOUT) :: PTS(D,N)
   INTEGER(C_INT), INTENT(IN) :: M
   REAL(C_DOUBLE), INTENT(INOUT) :: Q(D,M)
   INTEGER(C_INT), INTENT(OUT) :: SIMPS(D+1,M)
   REAL(C_DOUBLE), INTENT(OUT) :: WEIGHTS(D+1,M)
   INTEGER(C_INT), INTENT(OUT) :: IERR(M)
  
   INTERFACE
      SUBROUTINE DELAUNAYSPARSEP(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, &
                                 INTERP_IN, INTERP_OUT, EPS, EXTRAP,    &
                                 RNORM, IBUDGET, CHAIN, EXACT, PMODE)
         USE REAL_PRECISION , ONLY : R8
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: D
         INTEGER, INTENT(IN) :: N
         REAL(KIND=R8), INTENT(INOUT) :: PTS(:,:)
         INTEGER, INTENT(IN) :: M
         REAL(KIND=R8), INTENT(INOUT) :: Q(:,:)
         INTEGER, INTENT(OUT) :: SIMPS(:,:)
         REAL(KIND=R8), INTENT(OUT) :: WEIGHTS(:,:)
         INTEGER, INTENT(OUT) :: IERR(:)
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: INTERP_IN(:,:)
         REAL(KIND=R8), INTENT(OUT), OPTIONAL :: INTERP_OUT(:,:)
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: EPS
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: EXTRAP
         REAL(KIND=R8), INTENT(OUT), OPTIONAL :: RNORM(:)
         INTEGER, INTENT(IN), OPTIONAL :: IBUDGET
         LOGICAL, INTENT(IN), OPTIONAL :: CHAIN
         LOGICAL, INTENT(IN), OPTIONAL :: EXACT
         INTEGER, INTENT(IN), OPTIONAL :: PMODE
      END SUBROUTINE DELAUNAYSPARSEP
   END INTERFACE
  
   INTEGER :: D_LOC
   INTEGER :: N_LOC
   REAL(KIND=R8) :: PTS_LOC(D, N)
   INTEGER :: M_LOC
   REAL(KIND=R8) :: Q_LOC(D, M)
   INTEGER :: SIMPS_LOC(D+1, M)
   REAL(KIND=R8) :: WEIGHTS_LOC(D+1, M)
   INTEGER :: IERR_LOC(M)

   D_LOC = INT(D)
   N_LOC = INT(N)
   PTS_LOC = REAL(PTS, KIND=R8)
   M_LOC = INT(M)
   Q_LOC = REAL(Q, KIND=R8)

   CALL DELAUNAYSPARSEP(D_LOC, N_LOC, PTS_LOC, M_LOC, Q_LOC, SIMPS_LOC, &
                        WEIGHTS_LOC, IERR_LOC)

   PTS = REAL(PTS_LOC, KIND=C_DOUBLE)
   Q = REAL(Q_LOC, KIND=C_DOUBLE)
   SIMPS = INT(SIMPS_LOC, KIND=C_INT)
   WEIGHTS = REAL(WEIGHTS_LOC, KIND=C_DOUBLE)
   IERR = INT(IERR_LOC, KIND=C_INT)

   RETURN
END SUBROUTINE C_DELAUNAYSPARSEP_NOOPTS


SUBROUTINE C_DELAUNAYSPARSEP_INTERP(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, &
                                    IR, INTERP_IN, INTERP_OUT)             &
           BIND(C, NAME="c_delaunaysparsep_interp")
   ! This is a wrapper for DELAUNAYSPARSEP with INTERP_IN and INTERP_OUT
   ! specified, but no other optional arguments. Unlike the Fortran interface,
   ! in this interface the dimension of the response variables (IR) must
   ! be explicitly specified by an additional input, IR.
   !
   !
   ! On input:
   !
   ! D is the dimension of the space for PTS and Q.
   !
   ! N is the number of data points in PTS.
   !
   ! PTS(1:D,1:N) is a real valued matrix with N columns, each containing the
   !    coordinates of a single data point in R^D.
   !
   ! M is the number of interpolation points in Q.
   !
   ! Q(1:D,1:M) is a real valued matrix with M columns, each containing the
   !    coordinates of a single interpolation point in R^D.
   !
   ! IR is the dimension of the response variables.
   !
   ! INTERP_IN(1:IR,1:N) contains real valued response vectors for each of
   !    the data points in PTS on input. The first dimension of INTERP_IN is
   !    inferred to be the dimension of these response vectors, and the
   !    second dimension must match N.
   !
   !
   ! On output:
   !
   ! PTS and Q have been rescaled and shifted. All the data points in PTS
   !    are now contained in the unit hyperball in R^D, and the points in Q
   !    have been shifted and scaled accordingly in relation to PTS.
   !
   ! SIMPS(1:D+1,1:M) contains the D+1 integer indices (corresponding to columns
   !    in PTS) for the D+1 vertices of the Delaunay simplex containing each
   !    interpolation point in Q.
   !
   ! WEIGHTS(1:D+1,1:M) contains the D+1 real valued weights for expressing each
   !    point in Q as a convex combination of the D+1 corresponding vertices
   !    in SIMPS.
   !
   ! IERR(1:M) contains integer valued error flags associated with the
   !    computation of each of the M interpolation points in Q. The error
   !    codes are given in the definition of DELAUNAYSPARSEP in delsparse.f90.
   ! 
   ! INTERP_OUT(1:IR,1:M) contains real valued response vectors for each
   !    interpolation point in Q on output. The first dimension of INTERP_OU
   !    must match the first dimension of INTERP_IN, and the second dimension
   !    must match M.
   !
   !
   ! LAST UPDATE:
   !   11/2020 by THC
   !
   USE REAL_PRECISION , ONLY : R8
   USE ISO_C_BINDING

   IMPLICIT NONE
  
   INTEGER(C_INT), INTENT(IN) :: D
   INTEGER(C_INT), INTENT(IN) :: N
   REAL(C_DOUBLE), INTENT(INOUT) :: PTS(D,N)
   INTEGER(C_INT), INTENT(IN) :: M
   REAL(C_DOUBLE), INTENT(INOUT) :: Q(D,M)
   INTEGER(C_INT), INTENT(OUT) :: SIMPS(D+1,M)
   REAL(C_DOUBLE), INTENT(OUT) :: WEIGHTS(D+1,M)
   INTEGER(C_INT), INTENT(OUT) :: IERR(M)
   INTEGER(C_INT), INTENT(IN) :: IR
   REAL(C_DOUBLE), INTENT(IN) :: INTERP_IN(IR, N)
   REAL(C_DOUBLE), INTENT(OUT) :: INTERP_OUT(IR, M)
  
   INTERFACE
      SUBROUTINE DELAUNAYSPARSEP(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, &
                                 INTERP_IN, INTERP_OUT, EPS, EXTRAP,    &
                                 RNORM, IBUDGET, CHAIN, EXACT, PMODE)
         USE REAL_PRECISION , ONLY : R8
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: D
         INTEGER, INTENT(IN) :: N
         REAL(KIND=R8), INTENT(INOUT) :: PTS(:,:)
         INTEGER, INTENT(IN) :: M
         REAL(KIND=R8), INTENT(INOUT) :: Q(:,:)
         INTEGER, INTENT(OUT) :: SIMPS(:,:)
         REAL(KIND=R8), INTENT(OUT) :: WEIGHTS(:,:)
         INTEGER, INTENT(OUT) :: IERR(:)
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: INTERP_IN(:,:)
         REAL(KIND=R8), INTENT(OUT), OPTIONAL :: INTERP_OUT(:,:)
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: EPS
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: EXTRAP
         REAL(KIND=R8), INTENT(OUT), OPTIONAL :: RNORM(:)
         INTEGER, INTENT(IN), OPTIONAL :: IBUDGET
         LOGICAL, INTENT(IN), OPTIONAL :: CHAIN
         LOGICAL, INTENT(IN), OPTIONAL :: EXACT
         INTEGER, INTENT(IN), OPTIONAL :: PMODE
      END SUBROUTINE DELAUNAYSPARSEP
   END INTERFACE

   INTEGER :: D_LOC
   INTEGER :: N_LOC
   REAL(KIND=R8) :: PTS_LOC(D, N)
   INTEGER :: M_LOC
   REAL(KIND=R8) :: Q_LOC(D, M)
   INTEGER :: SIMPS_LOC(D+1, M)
   REAL(KIND=R8) :: WEIGHTS_LOC(D+1, M)
   INTEGER :: IERR_LOC(M)
   REAL(KIND=R8) :: INTERP_IN_LOC(IR, N)
   REAL(KIND=R8) :: INTERP_OUT_LOC(IR, M)

   D_LOC = INT(D)
   N_LOC = INT(N)
   PTS_LOC = REAL(PTS, KIND=R8)
   M_LOC = INT(M)
   Q_LOC = REAL(Q, KIND=R8)
   INTERP_IN_LOC = REAL(INTERP_IN, KIND=R8)

   CALL DELAUNAYSPARSEP(D_LOC, N_LOC, PTS_LOC, M_LOC, Q_LOC, SIMPS_LOC, &
                        WEIGHTS_LOC, IERR_LOC, INTERP_IN=INTERP_IN_LOC, &
                        INTERP_OUT=INTERP_OUT_LOC)

   PTS = REAL(PTS_LOC, KIND=C_DOUBLE)
   Q = REAL(Q_LOC, KIND=C_DOUBLE)
   SIMPS = INT(SIMPS_LOC, KIND=C_INT)
   WEIGHTS = REAL(WEIGHTS_LOC, KIND=C_DOUBLE)
   IERR = INT(IERR_LOC, KIND=C_INT)
   INTERP_OUT = REAL(INTERP_OUT_LOC, KIND=C_DOUBLE)

   RETURN
END SUBROUTINE C_DELAUNAYSPARSEP_INTERP


SUBROUTINE C_DELAUNAYSPARSEP_OPTS(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, EPS,  &
                                  EXTRAP, RNORM, IBUDGET, CHAIN, EXACT, PMODE) &
           BIND(C, NAME="c_delaunaysparsep_opts")
   ! This is a wrapper for DELAUNAYSPARSEP without INTERP_IN and INTERP_OUT,
   ! but all other optional arguments present.
   !
   !
   ! On input:
   !
   ! D is the dimension of the space for PTS and Q.
   !
   ! N is the number of data points in PTS.
   !
   ! PTS(1:D,1:N) is a real valued matrix with N columns, each containing the
   !    coordinates of a single data point in R^D.
   !
   ! M is the number of interpolation points in Q.
   !
   ! Q(1:D,1:M) is a real valued matrix with M columns, each containing the
   !    coordinates of a single interpolation point in R^D.
   !
   ! EXTRAP contains the real maximum extrapolation distance (relative to the
   !    diameter of PTS) on input. Interpolation at a point outside the convex
   !    hull of PTS is done by projecting that point onto the convex hull, and
   !    then doing normal Delaunay interpolation at that projection.
   !    Interpolation at any point in Q that is more than EXTRAP * DIAMETER(PTS)
   !    units outside the convex hull of PTS will not be done and an error code
   !    of 2 will be returned. Note that computing the projection can be
   !    expensive. Setting EXTRAP=0 will cause all extrapolation points to be
   !    ignored without ever computing a projection.
   ! 
   ! IBUDGET on input contains the integer budget for performing flips while
   !    iterating toward the simplex containing each interpolation point in Q.
   !    This prevents DELAUNAYSPARSEP from falling into an infinite loop when
   !    an inappropriate value of EPS is given with respect to the problem
   !    conditioning. For most cases, the default value of 50000 should be
   !    more than sufficient.
   !
   ! CHAIN is a logical input argument that determines whether a new first
   !    simplex should be constructed for each interpolation point
   !    (CHAIN=.FALSE.), or whether the simplex walks should be "daisy-chained."
   !    Setting CHAIN=.TRUE. is generally not recommended, unless the size of
   !    the triangulation is relatively small or the interpolation points are
   !    known to be tightly clustered.
   !
   ! EXACT is a logical input argument that determines whether the exact
   !    diameter should be computed and whether a check for duplicate data
   !    points should be performed in advance. When EXACT=.FALSE., the
   !    diameter of PTS is approximated by twice the distance from the
   !    barycenter of PTS to the farthest point in PTS, and no check is
   !    done to find the closest pair of points, which could result in hard
   !    to find bugs later on. When EXACT=.TRUE., the exact diameter is
   !    computed and an error is returned whenever PTS contains duplicate
   !    values up to the precision EPS. Setting EXACT=.FALSE. could result
   !    in significant speedup when N is large, but it is strongly
   !    recommended that most users leave EXACT=.TRUE., as setting
   !    EXACT=.FALSE. could result in input errors that are difficult
   !    to identify. Also, the diameter approximation could be wrong by up
   !    to a factor of two.
   !
   ! PMODE is an integer specifying the level of parallelism to be exploited.
   !    If PMODE = 1, then parallelism is exploited at the level of the loop
   !    over all interpolation points (Level 1 parallelism).
   !    If PMODE = 2, then parallelism is exploited at the level of the loops
   !    over data points when constructing/flipping simplices (Level 2
   !    parallelism).
   !    If PMODE = 3, then parallelism is exploited at both levels. Note: this
   !    implies that the total number of threads active at any time could be up
   !    to OMP_NUM_THREADS^2.
   !
   !
   ! On output:
   !
   ! PTS and Q have been rescaled and shifted. All the data points in PTS
   !    are now contained in the unit hyperball in R^D, and the points in Q
   !    have been shifted and scaled accordingly in relation to PTS.
   !
   ! SIMPS(1:D+1,1:M) contains the D+1 integer indices (corresponding to columns
   !    in PTS) for the D+1 vertices of the Delaunay simplex containing each
   !    interpolation point in Q.
   !
   ! WEIGHTS(1:D+1,1:M) contains the D+1 real valued weights for expressing each
   !    point in Q as a convex combination of the D+1 corresponding vertices
   !    in SIMPS.
   !
   ! IERR(1:M) contains integer valued error flags associated with the
   !    computation of each of the M interpolation points in Q. The error
   !    codes are given in the definition of DELAUNAYSPARSEP in delsparse.f90.
   ! 
   ! RNORM(1:M) contains the real unscaled projection (2-norm) distances from
   !    any projection computations on output.
   !
   !
   ! LAST UPDATE:
   !   11/2020 by THC
   !
   USE REAL_PRECISION , ONLY : R8
   USE ISO_C_BINDING

   IMPLICIT NONE
  
   INTEGER(C_INT), INTENT(IN) :: D
   INTEGER(C_INT), INTENT(IN) :: N
   REAL(C_DOUBLE), INTENT(INOUT) :: PTS(D,N)
   INTEGER(C_INT), INTENT(IN) :: M
   REAL(C_DOUBLE), INTENT(INOUT) :: Q(D,M)
   INTEGER(C_INT), INTENT(OUT) :: SIMPS(D+1,M)
   REAL(C_DOUBLE), INTENT(OUT) :: WEIGHTS(D+1,M)
   INTEGER(C_INT), INTENT(OUT) :: IERR(M)
   REAL(C_DOUBLE), INTENT(IN) :: EPS
   REAL(C_DOUBLE), INTENT(IN) :: EXTRAP
   REAL(C_DOUBLE), INTENT(OUT) :: RNORM(M)
   INTEGER(C_INT), INTENT(IN) :: IBUDGET
   LOGICAL(C_BOOL), INTENT(IN) :: CHAIN
   LOGICAL(C_BOOL), INTENT(IN) :: EXACT
   INTEGER(C_INT), INTENT(IN) :: PMODE
  
   INTERFACE
      SUBROUTINE DELAUNAYSPARSEP(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, &
                                 INTERP_IN, INTERP_OUT, EPS, EXTRAP,    &
                                 RNORM, IBUDGET, CHAIN, EXACT, PMODE)
         USE REAL_PRECISION , ONLY : R8
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: D
         INTEGER, INTENT(IN) :: N
         REAL(KIND=R8), INTENT(INOUT) :: PTS(:,:)
         INTEGER, INTENT(IN) :: M
         REAL(KIND=R8), INTENT(INOUT) :: Q(:,:)
         INTEGER, INTENT(OUT) :: SIMPS(:,:)
         REAL(KIND=R8), INTENT(OUT) :: WEIGHTS(:,:)
         INTEGER, INTENT(OUT) :: IERR(:)
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: INTERP_IN(:,:)
         REAL(KIND=R8), INTENT(OUT), OPTIONAL :: INTERP_OUT(:,:)
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: EPS
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: EXTRAP
         REAL(KIND=R8), INTENT(OUT), OPTIONAL :: RNORM(:)
         INTEGER, INTENT(IN), OPTIONAL :: IBUDGET
         LOGICAL, INTENT(IN), OPTIONAL :: CHAIN
         LOGICAL, INTENT(IN), OPTIONAL :: EXACT
         INTEGER, INTENT(IN), OPTIONAL :: PMODE
      END SUBROUTINE DELAUNAYSPARSEP
   END INTERFACE

   INTEGER :: D_LOC
   INTEGER :: N_LOC
   REAL(KIND=R8) :: PTS_LOC(D, N)
   INTEGER :: M_LOC
   REAL(KIND=R8) :: Q_LOC(D, M)
   INTEGER :: SIMPS_LOC(D+1, M)
   REAL(KIND=R8) :: WEIGHTS_LOC(D+1, M)
   INTEGER :: IERR_LOC(M)
   REAL(KIND=R8) :: EPS_LOC
   REAL(KIND=R8) :: EXTRAP_LOC
   REAL(KIND=R8) :: RNORM_LOC(M)
   INTEGER :: IBUDGET_LOC
   LOGICAL :: CHAIN_LOC
   LOGICAL :: EXACT_LOC
   INTEGER :: PMODE_LOC

   D_LOC = INT(D)
   N_LOC = INT(N)
   PTS_LOC = REAL(PTS, KIND=R8)
   M_LOC = INT(M)
   Q_LOC = REAL(Q, KIND=R8)
   EPS_LOC = REAL(EPS, KIND=R8)
   EXTRAP_LOC = REAL(EXTRAP, KIND=R8)
   IBUDGET_LOC = INT(IBUDGET)
   CHAIN_LOC = LOGICAL(CHAIN)
   EXACT_LOC = LOGICAL(EXACT)
   PMODE_LOC = INT(PMODE)

   CALL DELAUNAYSPARSEP(D_LOC, N_LOC, PTS_LOC, M_LOC, Q_LOC, SIMPS_LOC, &
                        WEIGHTS_LOC, IERR_LOC, EPS=EPS_LOC,             &
                        EXTRAP=EXTRAP_LOC, RNORM=RNORM_LOC,             &
                        IBUDGET=IBUDGET_LOC, CHAIN=CHAIN_LOC,           &
                        EXACT=EXACT_LOC, PMODE=PMODE_LOC)

   PTS = REAL(PTS_LOC, KIND=C_DOUBLE)
   Q = REAL(Q_LOC, KIND=C_DOUBLE)
   SIMPS = INT(SIMPS_LOC, KIND=C_INT)
   WEIGHTS = REAL(WEIGHTS_LOC, KIND=C_DOUBLE)
   IERR = INT(IERR_LOC, KIND=C_INT)
   RNORM = REAL(RNORM_LOC, KIND=C_DOUBLE)

   RETURN
END SUBROUTINE C_DELAUNAYSPARSEP_OPTS


SUBROUTINE C_DELAUNAYSPARSEP_INTERP_OPTS(D, N, PTS, M, Q, SIMPS, WEIGHTS,    &
                                         IERR, IR, INTERP_IN, INTERP_OUT,    &
                                         EPS, EXTRAP, RNORM, IBUDGET, CHAIN, &
                                         EXACT, PMODE)                       &
           BIND(C, NAME="c_delaunaysparsep_interp_opts")
   ! This is a wrapper for DELAUNAYSPARSEP with all optional arguments present.
   !
   !
   ! On input:
   !
   ! D is the dimension of the space for PTS and Q.
   !
   ! N is the number of data points in PTS.
   !
   ! PTS(1:D,1:N) is a real valued matrix with N columns, each containing the
   !    coordinates of a single data point in R^D.
   !
   ! M is the number of interpolation points in Q.
   !
   ! Q(1:D,1:M) is a real valued matrix with M columns, each containing the
   !    coordinates of a single interpolation point in R^D.
   !
   ! IR is the dimension of the response variables.
   !
   ! INTERP_IN(1:IR,1:N) contains real valued response vectors for each of
   !    the data points in PTS on input. The first dimension of INTERP_IN is
   !    inferred to be the dimension of these response vectors, and the
   !    second dimension must match N.
   !
   ! EXTRAP contains the real maximum extrapolation distance (relative to the
   !    diameter of PTS) on input. Interpolation at a point outside the convex
   !    hull of PTS is done by projecting that point onto the convex hull, and
   !    then doing normal Delaunay interpolation at that projection.
   !    Interpolation at any point in Q that is more than EXTRAP * DIAMETER(PTS)
   !    units outside the convex hull of PTS will not be done and an error code
   !    of 2 will be returned. Note that computing the projection can be
   !    expensive. Setting EXTRAP=0 will cause all extrapolation points to be
   !    ignored without ever computing a projection.
   ! 
   ! IBUDGET on input contains the integer budget for performing flips while
   !    iterating toward the simplex containing each interpolation point in Q.
   !    This prevents DELAUNAYSPARSEP from falling into an infinite loop when
   !    an inappropriate value of EPS is given with respect to the problem
   !    conditioning. For most cases, the default value of 50000 should be
   !    more than sufficient.
   !
   ! CHAIN is a logical input argument that determines whether a new first
   !    simplex should be constructed for each interpolation point
   !    (CHAIN=.FALSE.), or whether the simplex walks should be "daisy-chained."
   !    Setting CHAIN=.TRUE. is generally not recommended, unless the size of
   !    the triangulation is relatively small or the interpolation points are
   !    known to be tightly clustered.
   !
   ! EXACT is a logical input argument that determines whether the exact
   !    diameter should be computed and whether a check for duplicate data
   !    points should be performed in advance. When EXACT=.FALSE., the
   !    diameter of PTS is approximated by twice the distance from the
   !    barycenter of PTS to the farthest point in PTS, and no check is
   !    done to find the closest pair of points, which could result in hard
   !    to find bugs later on. When EXACT=.TRUE., the exact diameter is
   !    computed and an error is returned whenever PTS contains duplicate
   !    values up to the precision EPS. Setting EXACT=.FALSE. could result
   !    in significant speedup when N is large, but it is strongly
   !    recommended that most users leave EXACT=.TRUE., as setting
   !    EXACT=.FALSE. could result in input errors that are difficult
   !    to identify. Also, the diameter approximation could be wrong by up
   !    to a factor of two.
   !
   ! PMODE is an integer specifying the level of parallelism to be exploited.
   !    If PMODE = 1, then parallelism is exploited at the level of the loop
   !    over all interpolation points (Level 1 parallelism).
   !    If PMODE = 2, then parallelism is exploited at the level of the loops
   !    over data points when constructing/flipping simplices (Level 2
   !    parallelism).
   !    If PMODE = 3, then parallelism is exploited at both levels. Note: this
   !    implies that the total number of threads active at any time could be up
   !    to OMP_NUM_THREADS^2.
   !
   !
   ! On output:
   !
   ! PTS and Q have been rescaled and shifted. All the data points in PTS
   !    are now contained in the unit hyperball in R^D, and the points in Q
   !    have been shifted and scaled accordingly in relation to PTS.
   !
   ! SIMPS(1:D+1,1:M) contains the D+1 integer indices (corresponding to columns
   !    in PTS) for the D+1 vertices of the Delaunay simplex containing each
   !    interpolation point in Q.
   !
   ! WEIGHTS(1:D+1,1:M) contains the D+1 real valued weights for expressing each
   !    point in Q as a convex combination of the D+1 corresponding vertices
   !    in SIMPS.
   !
   ! IERR(1:M) contains integer valued error flags associated with the
   !    computation of each of the M interpolation points in Q. The error
   !    codes are given in the definition of DELAUNAYSPARSEP in delsparse.f90.
   ! 
   ! INTERP_OUT(1:IR,1:M) contains real valued response vectors for each
   !    interpolation point in Q on output. The first dimension of INTERP_OUT
   !    must match the first dimension of INTERP_IN, and the second dimension
   !    must match M.
   ! 
   ! RNORM(1:M) contains the real unscaled projection (2-norm) distances from
   !    any projection computations on output.
   !
   !
   ! LAST UPDATE:
   !   11/2020 by THC
   !
   USE REAL_PRECISION , ONLY : R8
   USE ISO_C_BINDING

   IMPLICIT NONE
  
   INTEGER(C_INT), INTENT(IN) :: D
   INTEGER(C_INT), INTENT(IN) :: N
   REAL(C_DOUBLE), INTENT(INOUT) :: PTS(D,N)
   INTEGER(C_INT), INTENT(IN) :: M
   REAL(C_DOUBLE), INTENT(INOUT) :: Q(D,M)
   INTEGER(C_INT), INTENT(OUT) :: SIMPS(D+1,M)
   REAL(C_DOUBLE), INTENT(OUT) :: WEIGHTS(D+1,M)
   INTEGER(C_INT), INTENT(OUT) :: IERR(M)
   INTEGER(C_INT), INTENT(IN) :: IR
   REAL(C_DOUBLE), INTENT(IN) :: INTERP_IN(IR, N)
   REAL(C_DOUBLE), INTENT(OUT) :: INTERP_OUT(IR, M)
   REAL(C_DOUBLE), INTENT(IN) :: EPS
   REAL(C_DOUBLE), INTENT(IN) :: EXTRAP
   REAL(C_DOUBLE), INTENT(OUT) :: RNORM(M)
   INTEGER(C_INT), INTENT(IN) :: IBUDGET
   LOGICAL(C_BOOL), INTENT(IN) :: CHAIN
   LOGICAL(C_BOOL), INTENT(IN) :: EXACT
   INTEGER(C_INT), INTENT(IN) :: PMODE
  
   INTERFACE
      SUBROUTINE DELAUNAYSPARSEP(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, &
                                 INTERP_IN, INTERP_OUT, EPS, EXTRAP,    &
                                 RNORM, IBUDGET, CHAIN, EXACT, PMODE)
         USE REAL_PRECISION , ONLY : R8
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: D
         INTEGER, INTENT(IN) :: N
         REAL(KIND=R8), INTENT(INOUT) :: PTS(:,:)
         INTEGER, INTENT(IN) :: M
         REAL(KIND=R8), INTENT(INOUT) :: Q(:,:)
         INTEGER, INTENT(OUT) :: SIMPS(:,:)
         REAL(KIND=R8), INTENT(OUT) :: WEIGHTS(:,:)
         INTEGER, INTENT(OUT) :: IERR(:)
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: INTERP_IN(:,:)
         REAL(KIND=R8), INTENT(OUT), OPTIONAL :: INTERP_OUT(:,:)
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: EPS
         REAL(KIND=R8), INTENT(IN), OPTIONAL :: EXTRAP
         REAL(KIND=R8), INTENT(OUT), OPTIONAL :: RNORM(:)
         INTEGER, INTENT(IN), OPTIONAL :: IBUDGET
         LOGICAL, INTENT(IN), OPTIONAL :: CHAIN
         LOGICAL, INTENT(IN), OPTIONAL :: EXACT
         INTEGER, INTENT(IN), OPTIONAL :: PMODE
      END SUBROUTINE DELAUNAYSPARSEP
   END INTERFACE

   INTEGER :: D_LOC
   INTEGER :: N_LOC
   REAL(KIND=R8) :: PTS_LOC(D, N)
   INTEGER :: M_LOC
   REAL(KIND=R8) :: Q_LOC(D, M)
   INTEGER :: SIMPS_LOC(D+1, M)
   REAL(KIND=R8) :: WEIGHTS_LOC(D+1, M)
   INTEGER :: IERR_LOC(M)
   REAL(KIND=R8) :: INTERP_IN_LOC(IR, N)
   REAL(KIND=R8) :: INTERP_OUT_LOC(IR, M)
   REAL(KIND=R8) :: EPS_LOC
   REAL(KIND=R8) :: EXTRAP_LOC
   REAL(KIND=R8) :: RNORM_LOC(M)
   INTEGER :: IBUDGET_LOC
   LOGICAL :: CHAIN_LOC
   LOGICAL :: EXACT_LOC
   INTEGER :: PMODE_LOC

   D_LOC = INT(D)
   N_LOC = INT(N)
   PTS_LOC = REAL(PTS, KIND=R8)
   M_LOC = INT(M)
   Q_LOC = REAL(Q, KIND=R8)
   INTERP_IN_LOC = REAL(INTERP_IN, KIND=R8)
   EPS_LOC = REAL(EPS, KIND=R8)
   EXTRAP_LOC = REAL(EXTRAP, KIND=R8)
   IBUDGET_LOC = INT(IBUDGET)
   CHAIN_LOC = LOGICAL(CHAIN)
   EXACT_LOC = LOGICAL(EXACT)
   PMODE_LOC = INT(PMODE)

   CALL DELAUNAYSPARSEP(D_LOC, N_LOC, PTS_LOC, M_LOC, Q_LOC, SIMPS_LOC, &
                        WEIGHTS_LOC, IERR_LOC, INTERP_IN=INTERP_IN_LOC, &
                        INTERP_OUT=INTERP_OUT_LOC, EPS=EPS_LOC,         &
                        EXTRAP=EXTRAP_LOC, RNORM=RNORM_LOC,             &
                        IBUDGET=IBUDGET_LOC, CHAIN=CHAIN_LOC,           &
                        EXACT=EXACT_LOC, PMODE=PMODE_LOC)

   PTS = REAL(PTS_LOC, KIND=C_DOUBLE)
   Q = REAL(Q_LOC, KIND=C_DOUBLE)
   SIMPS = INT(SIMPS_LOC, KIND=C_INT)
   WEIGHTS = REAL(WEIGHTS_LOC, KIND=C_DOUBLE)
   IERR = INT(IERR_LOC, KIND=C_INT)
   INTERP_OUT = REAL(INTERP_OUT_LOC, C_DOUBLE)
   RNORM = REAL(RNORM_LOC, KIND=C_DOUBLE)

   RETURN
END SUBROUTINE C_DELAUNAYSPARSEP_INTERP_OPTS

