# Usage Information for DELAUNAYSPARSE.

DELAUNAYSPARSE solves the multivariate interpolation problem:

Given a set of `N` points `PTS` and a set of `M` interpolation points
`Q` in `R^D`, for each interpolation point `Q_i` in `Q`, identify the set
of `D+1` data points in `PTS` that are the vertices of a Delaunay simplex
containing `Q_i`.

These vertices can be used to calculate the Delaunay interpolant using
a piecewise linear model.

For more information on the underlying algorithm, see

Chang et al. 2018. A polynomial time algorithm for multivariate interpolation
in arbitrary dimension via the Delaunay triangulation. In Proc. ACMSE 2018
Conference.

For more information on this software, see

Chang et al. 2020. Algorithm 1012: DELAUNAYSPARSE: Interpolation via a sparse
subset of the Delaunay triangulation in medium to high dimensions.
ACM Trans. Math. Softw. 46(4). Article No. 38.

DELAUNAYSPARSE contains a Fortran module
 * `delsparse`;

as well as C bindings
 * `delsparsec`;

two command-line drivers
 * `delsparses` and
 * `delsparsep`;

and a Python 3 wrapper
 * `delsparsepy`.

These interfaces are described in the following sections.

## Fortran interface

DELAUNAYSPARSE is written in Fortran 2003, and this is its native interface.
The Fortran interface contains two drivers:
 * `DELAUNAYSPARSES` (serial driver) and
 * `DELAUNAYSPARSEP` (OpenMP parallel driver).

### DELAUNAYSPARSES

The interface for DELAUNAYSPARSES is

```
SUBROUTINE DELAUNAYSPARSES( D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR,     &
                            INTERP_IN, INTERP_OUT, EPS, EXTRAP, RNORM, &
                            IBUDGET, CHAIN, EXACT                      )
```

Each of the above parameters is described below.


On input:

 * `D` is the dimension of the space for `PTS` and `Q`.

 * `N` is the number of data points in `PTS`.

 * `PTS(1:D,1:N)` is a real valued matrix with `N` columns, each containing the
   coordinates of a single data point in R^D.

 * `M` is the number of interpolation points in `Q`.

 * `Q(1:D,1:M)` is a real valued matrix with `M` columns, each containing the
   coordinates of a single interpolation point in R^D.


On output:

 * `PTS` and `Q` have been rescaled and shifted. All the data points in `PTS`
   are now contained in the unit hyperball in R^D, and the points in `Q`
   have been shifted and scaled accordingly in relation to `PTS`.

 * `SIMPS(1:D+1,1:M)` contains the `D+1` integer indices (corresponding to
   columns in `PTS`) for the `D+1` vertices of the Delaunay simplex
   containing each interpolation point in `Q`.

 * `WEIGHTS(1:D+1,1:M)` contains the `D+1` real-valued weights for expressing
   each point in `Q` as a convex combination of the `D+1` corresponding vertices
   in `SIMPS`.

 * `IERR(1:M)` contains integer valued error flags associated with the
   computation of each of the `M` interpolation points in `Q`. The error
   codes are:

   Codes 0, 1, 2 are expected to occur during normal execution.

    - 00 : Succesful interpolation.
    - 01 : Succesful extrapolation (up to the allowed extrapolation distance).
    - 02 : This point was outside the allowed extrapolation distance; the
           corresponding entries in SIMPS and WEIGHTS contain zero values.

    Error codes 10--28 indicate that one or more inputs contain illegal
    values or are incompatible with each other.

    - 10 : The dimension `D` must be positive.
    - 11 : Too few data points to construct a triangulation (i.e., `N < D+1`).
    - 12 : No interpolation points given (i.e., `M < 1`).
    - 13 : The first dimension of `PTS` does not agree with the dimension `D`.
    - 14 : The second dimension of `PTS` does not agree with the number of
           points `N`.
    - 15 : The first dimension of `Q` does not agree with the dimension `D`.
    - 16 : The second dimension of `Q` does not agree with the number of
           interpolation points `M`.
    - 17 : The first dimension of the output array `SIMPS` does not match the
           number of vertices needed for a `D`-simplex (`D+1`).
    - 18 : The second dimension of the output array `SIMPS` does not match the
           number of interpolation points `M`.
    - 19 : The first dimension of the output array `WEIGHTS` does not match the
           number of vertices for a a `D`-simplex (`D+1`).
    - 20 : The second dimension of the output array `WEIGHTS` does not match
           the number of interpolation points `M`.
    - 21 : The size of the error array `IERR` does not match the number of
           interpolation points `M`.
    - 22 : `INTERP_IN` cannot be present without `INTERP_OUT` or vice versa.
    - 23 : The first dimension of `INTERP_IN` does not match the first
           dimension of `INTERP_OUT`.
    - 24 : The second dimension of `INTERP_IN` does not match the number of
           data points `PTS`.
    - 25 : The second dimension of `INTERP_OUT` does not match the number of
           interpolation points `M`.
    - 26 : The budget supplied in `IBUDGET` does not contain a positive
           integer.
    - 27 : The extrapolation distance supplied in `EXTRAP` cannot be negative.
    - 28 : The size of the `RNORM` output array does not match the number of
           interpolation points `M`.

    The errors 30, 31 typically indicate that DELAUNAYSPARSE has been given
    an unclean dataset. These errors can be fixed by preprocessing your
    data (remove duplicate points and apply PCA or other dimension reduction
    technique).

    - 30 : Two or more points in the data set `PTS` are too close together with
           respect to the working precision (`EPS`), which would result in a
           numerically degenerate simplex.
    - 31 : All the data points in `PTS` lie in some lower dimensional linear
           manifold (up to the working precision), and no valid triangulation
           exists.

    The error code 40 occurs when another earlier error prevented this point
    from ever being evaluated.

    - 40 : An error caused `DELAUNAYSPARSES` to terminate before this value
           could be computed. Note: The corresponding entries in `SIMPS` and
           `WEIGHTS` may contain garbage values.

    The error code 50 corresponds to allocation of the internal WORK array.
    Check your systems internal memory settings and limits, in relation
    to the problem size and DELAUNAYSPARSE's space requirements (see TOMS
    Alg. paper for more details on DELAUNAYSPARSE's space requirements).

    - 50 : A memory allocation error occurred while allocating the work array
           `WORK`.

    The errors 60, 61 should not occur with the default settings. If one of
    these errors is observed, then it is likely that either the value of
    the optional inputs `IBUDGET` or `EPS` has been adjusted in a way that is
    unwise, or there may be another issue with the problem settings, which
    is manifesting in an unusual way.

    - 60 : The budget was exceeded before the algorithm converged on this
           value. If the dimension is high, try increasing `IBUDGET`. This
           error can also be caused by a working precision `EPS` that is too
           small for the conditioning of the problem.
    - 61 : A value that was judged appropriate later caused LAPACK to
           encounter a singularity. Try increasing the value of `EPS`.

    The errors 70+ can occur during extrapolation outside the convex hull.
    Errors 71+ are received from the solver BQPD and should never occur.
    Since extrapolation has a much higher memory overhead than interpolation,
    if one of these errors occurs, it is likely that your system does not have
    sufficient memory to support extrapolation for the given problem size.

    - 70 : Allocation error for the extrapolation work arrays. Note that
           extrapolation has a higher memory overhead than interpolation for
           the current version.
    - 71 : Unbounded problem detected. This error should never occur
           unless there is a user error.
    - 72 : bl(i) > bu(i) for some i This error should never occur
           unless there is a user error.
    - 73 : Infeasible problem detected in Phase 1 (should never occur).
    - 74 : Incorrect setting of m, n, kmax, mlp, mode or tol. This is
           extremely unlikely, but could indicate that the problem is
           too large for usage of BQPD.
    - 75 : Not enough space in lp. This error should not occur.
           Double-check your usage, then contact authors if the issue
           persists.
    - 76 : Not enough space for reduced Hessian matrix (increase kmax).
           This error should not occur. Double-check your usage, then
           contact authors if the issue persists.
    - 77 : Not enough space for sparse factors. This error should not
           occur, but may indicate an un-anticipated problem size.
           Double-check your usage then contact the authors if the issue
           persists.
    - 78 : Maximum number of unsuccessful restarts taken. This issue
           should never occur since the projection problem is convex.

   The errors 80--83 should never occur, and likely indicate a
   compiler bug or hardware failure.

    - 80 : The LAPACK subroutine `DGEQP3` has reported an illegal value.
    - 81 : The LAPACK subroutine `DGETRF` has reported an illegal value.
    - 82 : The LAPACK subroutine `DGETRS` has reported an illegal value.
    - 83 : The LAPACK subroutine `DORMQR` has reported an illegal value.


Optional arguments:

 * `INTERP_IN(1:IR,1:N)` contains real valued response vectors for each of
   the data points in `PTS` on input. The first dimension of `INTERP_IN` is
   inferred to be the dimension of these response vectors, and the
   second dimension must match `N`. If present, the response values will
   be computed for each interpolation point in `Q`, and stored in `INTERP_OUT`,
   which therefore must also be present. If both `INTERP_IN` and `INTERP_OUT`
   are omitted, only the containing simplices and convex combination
   weights are returned.

 * `INTERP_OUT(1:IR,1:M)` contains real valued response vectors for each
   interpolation point in `Q` on output. The first dimension of `INTERP_OUT`
   must match the first dimension of `INTERP_IN`, and the second dimension
   must match `M`. If present, the response values at each interpolation
   point are computed as a convex combination of the response values
   (supplied in `INTERP_IN`) at the vertices of a Delaunay simplex containing
   that interpolation point. Therefore, if `INTERP_OUT` is present, then
   `INTERP_IN` must also be present.  If both are omitted, only the
   simplices and convex combination weights are returned.

 * `EPS` contains the real working precision for the problem on input. By
   default, `EPS` is assigned \sqrt{\mu} where \mu denotes the unit roundoff
   for the machine. In general, any values that differ by less than `EPS`
   are judged as equal, and any weights that are greater than `-EPS` are
   judged as nonnegative. `EPS` cannot take a value less than the default
   value of \sqrt{\mu}. If any value less than \sqrt{\mu} is supplied,
   the default value will be used instead automatically. 

 * `EXTRAP` contains the real maximum extrapolation distance (relative to the
   diameter of `PTS`) on input. Interpolation at a point outside the convex
   hull of `PTS` is done by projecting that point onto the convex hull, and
   then doing normal Delaunay interpolation at that projection.
   Interpolation at any point in `Q` that is more than `EXTRAP * DIAMETER(PTS)`
   units outside the convex hull of `PTS` will not be done and an error code
   of `2` will be returned. Note that computing the projection can be
   expensive. Setting `EXTRAP=0` will cause all extrapolation points to be
   ignored without ever computing a projection. By default, `EXTRAP=0.1`
   (extrapolate by up to 10% of the diameter of `PTS`). 

 * `RNORM(1:M)` contains the real unscaled projection (2-norm) distances from
   any projection computations on output. If not present, these distances
   are still computed for each extrapolation point, but are never returned.

 * `IBUDGET` on input contains the integer budget for performing flips while
   iterating toward the simplex containing each interpolation point in
   `Q`. This prevents `DELAUNAYSPARSES` from falling into an infinite loop when
   an inappropriate value of `EPS` is given with respect to the problem
   conditioning.  By default, `IBUDGET=50000`. However, for extremely
   high-dimensional problems and pathological inputs, the default value
   may be insufficient.

 * `CHAIN` is a logical input argument that determines whether a new first
   simplex should be constructed for each interpolation point
   (`CHAIN=.FALSE.`), or whether the simplex walks should be "daisy-chained."
   By default, `CHAIN=.FALSE.` Setting `CHAIN=.TRUE.` is generally not
   recommended, unless the size of the triangulation is relatively small
   or the interpolation points are known to be tightly clustered.

 * `EXACT` is a logical input argument that determines whether the exact
   diameter should be computed and whether a check for duplicate data
   points should be performed in advance. These checks are $O(N^2 D)$ time
   complexity, while `DELAUNAYSPARSE` tends toward $O(N D^4)$ on average.
   By default, `EXACT=.TRUE.` and the exact diameter is computed and an error
   is returned whenever `PTS` contains duplicate values up to the precision
   `EPS`. When `EXACT=.FALSE.`, the diameter of `PTS` is approximated by twice
   the distance from the barycenter of `PTS` to the farthest point in `PTS`,
   and no check is done to find the closest pair of points.
   When `EXACT=.TRUE.`, `DELAUNAYSPARSE` could spend over 90% of runtime
   calculating these constants, which are not critical to the `DELAUNAYSPARSE`
   algorithm. In particular, this happens for large values of `N`. However,
   setting `EXACT=.FALSE.` could result in input errors that are difficult
   to identify. It is recommended that users verify the input set `PTS`
   and possibly rescale `PTS` manually while `EXACT=.TRUE.` Then, when
   100% sure that `PTS` is valid, users may choose to set `EXACT=.FALSE.`
   in production runs for large values of N to achieve massive speedups.


Subroutines and functions directly referenced from BLAS are
 * `DDOT`,
 * `DGEMV`,
 * `DNRM2`,
 * `DTRSM`,
and from LAPACK are
 * `DGEQP3`,
 * `DGETRF`,
 * `DGETRS`,
 * `DORMQR`.

The BQPD driver
 * `BQPD` is also directly referenced.

`BQPD` and all its dependencies have been flattened into a single file `bqpd.f`.
For a reference to `BQPD`, see
   Annals of Operations Research, 46 : 307--334 (1993).
The module `REAL_PRECISION` from HOMPACK90 (ACM TOMS Algorithm 777) is
used for the real data type. The `REAL_PRECISION` module, `DELAUNAYSPARSEP`,
and `BQPD` and its dependencies comply with the Fortran 2008 standard.  

## DELAUNAYSPARSEP

```
SUBROUTINE DELAUNAYSPARSEP( D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR,  & 
  INTERP_IN, INTERP_OUT, EPS, EXTRAP, RNORM, IBUDGET, CHAIN, EXACT, &
  PMODE )
```

Each of the above parameters is described below.


On input:

 * `D` is the dimension of the space for `PTS` and `Q`.

 * `N` is the number of data points in `PTS`.

 * `PTS(1:D,1:N)` is a real valued matrix with `N` columns, each containing the
   coordinates of a single data point in R^D.

 * `M` is the number of interpolation points in `Q`.

 * `Q(1:D,1:M)` is a real valued matrix with `M` columns, each containing the
   coordinates of a single interpolation point in R^D.


On output:

 * `PTS` and `Q` have been rescaled and shifted. All the data points in `PTS`
   are now contained in the unit hyperball in R^D, and the points in `Q`
   have been shifted and scaled accordingly in relation to `PTS`.

 * `SIMPS(1:D+1,1:M)` contains the `D+1` integer indices (corresponding to
   columns in `PTS`) for the `D+1` vertices of the Delaunay simplex
   containing each interpolation point in `Q`.

 * `WEIGHTS(1:D+1,1:M)` contains the `D+1` real-valued weights for expressing
   each point in `Q` as a convex combination of the `D+1` corresponding vertices
   in `SIMPS`.

 * `IERR(1:M)` contains integer valued error flags associated with the
   computation of each of the `M` interpolation points in `Q`. The error
   codes are:

   Codes 0, 1, 2 are expected to occur during normal execution.

    - 00 : Succesful interpolation.
    - 01 : Succesful extrapolation (up to the allowed extrapolation distance).
    - 02 : This point was outside the allowed extrapolation distance; the
           corresponding entries in SIMPS and WEIGHTS contain zero values.

    Error codes 10--28 indicate that one or more inputs contain illegal
    values or are incompatible with each other.

    - 10 : The dimension `D` must be positive.
    - 11 : Too few data points to construct a triangulation (i.e., `N < D+1`).
    - 12 : No interpolation points given (i.e., `M < 1`).
    - 13 : The first dimension of `PTS` does not agree with the dimension `D`.
    - 14 : The second dimension of `PTS` does not agree with the number of
           points `N`.
    - 15 : The first dimension of `Q` does not agree with the dimension `D`.
    - 16 : The second dimension of `Q` does not agree with the number of
           interpolation points `M`.
    - 17 : The first dimension of the output array `SIMPS` does not match the
           number of vertices needed for a `D`-simplex (`D+1`).
    - 18 : The second dimension of the output array `SIMPS` does not match the
           number of interpolation points `M`.
    - 19 : The first dimension of the output array `WEIGHTS` does not match the
           number of vertices for a a `D`-simplex (`D+1`).
    - 20 : The second dimension of the output array `WEIGHTS` does not match
           the number of interpolation points `M`.
    - 21 : The size of the error array `IERR` does not match the number of
           interpolation points `M`.
    - 22 : `INTERP_IN` cannot be present without `INTERP_OUT` or vice versa.
    - 23 : The first dimension of `INTERP_IN` does not match the first
           dimension of `INTERP_OUT`.
    - 24 : The second dimension of `INTERP_IN` does not match the number of
           data points `PTS`.
    - 25 : The second dimension of `INTERP_OUT` does not match the number of
           interpolation points `M`.
    - 26 : The budget supplied in `IBUDGET` does not contain a positive
           integer.
    - 27 : The extrapolation distance supplied in `EXTRAP` cannot be negative.
    - 28 : The size of the `RNORM` output array does not match the number of
           interpolation points `M`.

    The errors 30, 31 typically indicate that DELAUNAYSPARSE has been given
    an unclean dataset. These errors can be fixed by preprocessing your
    data (remove duplicate points and apply PCA or other dimension reduction
    technique).

    - 30 : Two or more points in the data set `PTS` are too close together with
           respect to the working precision (`EPS`), which would result in a
           numerically degenerate simplex.
    - 31 : All the data points in `PTS` lie in some lower dimensional linear
           manifold (up to the working precision), and no valid triangulation
           exists.

    The error code 40 occurs when another earlier error prevented this point
    from ever being evaluated.

    - 40 : An error caused `DELAUNAYSPARSEP` to terminate before this value
           could be computed. Note: The corresponding entries in `SIMPS` and
           `WEIGHTS` may contain garbage values.

    The error code 50 corresponds to allocation of the internal WORK array.
    Check your systems internal memory settings and limits, in relation
    to the problem size and DELAUNAYSPARSE's space requirements (see TOMS
    Alg. paper for more details on DELAUNAYSPARSE's space requirements).

    - 50 : A memory allocation error occurred while allocating the work array
           `WORK`.

    The errors 60, 61 should not occur with the default settings. If one of
    these errors is observed, then it is likely that either the value of
    the optional inputs `IBUDGET` or `EPS` has been adjusted in a way that is
    unwise, or there may be another issue with the problem settings, which
    is manifesting in an unusual way.

    - 60 : The budget was exceeded before the algorithm converged on this
           value. If the dimension is high, try increasing `IBUDGET`. This
           error can also be caused by a working precision `EPS` that is too
           small for the conditioning of the problem.
    - 61 : A value that was judged appropriate later caused LAPACK to
           encounter a singularity. Try increasing the value of `EPS`.

    The errors 70+ can occur during extrapolation outside the convex hull.
    Errors 71+ are received from the solver BQPD and should never occur.
    Since extrapolation has a much higher memory overhead than interpolation,
    if one of these errors occurs, it is likely that your system does not have
    sufficient memory to support extrapolation for the given problem size.

    - 70 : Allocation error for the extrapolation work arrays. Note that
           extrapolation has a higher memory overhead than interpolation for
           the current version.
    - 71 : Unbounded problem detected. This error should never occur
           unless there is a user error.
    - 72 : bl(i) > bu(i) for some i This error should never occur
           unless there is a user error.
    - 73 : Infeasible problem detected in Phase 1 (should never occur).
    - 74 : Incorrect setting of m, n, kmax, mlp, mode or tol. This is
           extremely unlikely, but could indicate that the problem is
           too large for usage of BQPD.
    - 75 : Not enough space in lp. This error should not occur.
           Double-check your usage, then contact authors if the issue
           persists.
    - 76 : Not enough space for reduced Hessian matrix (increase kmax).
           This error should not occur. Double-check your usage, then
           contact authors if the issue persists.
    - 77 : Not enough space for sparse factors. This error should not
           occur, but may indicate an un-anticipated problem size.
           Double-check your usage then contact the authors if the issue
           persists.
    - 78 : Maximum number of unsuccessful restarts taken. This issue
           should never occur since the projection problem is convex.


   The errors 80--83 should never occur, and likely indicate a
   compiler bug or hardware failure.

    - 80 : The LAPACK subroutine `DGEQP3` has reported an illegal value.
    - 81 : The LAPACK subroutine `DGETRF` has reported an illegal value.
    - 82 : The LAPACK subroutine `DGETRS` has reported an illegal value.
    - 83 : The LAPACK subroutine `DORMQR` has reported an illegal value.

   The error code 90 is unique to DELAUNAYSPARSEP.

    - 90 : The value of `PMODE` is not valid.


Optional arguments:

 * `INTERP_IN(1:IR,1:N)` contains real valued response vectors for each of
   the data points in `PTS` on input. The first dimension of `INTERP_IN` is
   inferred to be the dimension of these response vectors, and the
   second dimension must match `N`. If present, the response values will
   be computed for each interpolation point in `Q`, and stored in `INTERP_OUT`,
   which therefore must also be present. If both `INTERP_IN` and `INTERP_OUT`
   are omitted, only the containing simplices and convex combination
   weights are returned.

 * `INTERP_OUT(1:IR,1:M)` contains real valued response vectors for each
   interpolation point in `Q` on output. The first dimension of `INTERP_OUT`
   must match the first dimension of `INTERP_IN`, and the second dimension
   must match `M`. If present, the response values at each interpolation
   point are computed as a convex combination of the response values
   (supplied in `INTERP_IN`) at the vertices of a Delaunay simplex containing
   that interpolation point. Therefore, if `INTERP_OUT` is present, then
   `INTERP_IN` must also be present.  If both are omitted, only the
   simplices and convex combination weights are returned.

 * `EPS` contains the real working precision for the problem on input. By
   default, `EPS` is assigned \sqrt{\mu} where \mu denotes the unit roundoff
   for the machine. In general, any values that differ by less than `EPS`
   are judged as equal, and any weights that are greater than `-EPS` are
   judged as nonnegative. `EPS` cannot take a value less than the default
   value of \sqrt{\mu}. If any value less than \sqrt{\mu} is supplied,
   the default value will be used instead automatically. 

 * `EXTRAP` contains the real maximum extrapolation distance (relative to the
   diameter of `PTS`) on input. Interpolation at a point outside the convex
   hull of `PTS` is done by projecting that point onto the convex hull, and
   then doing normal Delaunay interpolation at that projection.
   Interpolation at any point in `Q` that is more than `EXTRAP * DIAMETER(PTS)`
   units outside the convex hull of `PTS` will not be done and an error code
   of `2` will be returned. Note that computing the projection can be
   expensive. Setting `EXTRAP=0` will cause all extrapolation points to be
   ignored without ever computing a projection. By default, `EXTRAP=0.1`
   (extrapolate by up to 10% of the diameter of `PTS`). 

 * `RNORM(1:M)` contains the real unscaled projection (2-norm) distances from
   any projection computations on output. If not present, these distances
   are still computed for each extrapolation point, but are never returned.

 * `IBUDGET` on input contains the integer budget for performing flips while
   iterating toward the simplex containing each interpolation point in
   `Q`. This prevents `DELAUNAYSPARSEP` from falling into an infinite loop when
   an inappropriate value of `EPS` is given with respect to the problem
   conditioning.  By default, `IBUDGET=50000`. However, for extremely
   high-dimensional problems and pathological inputs, the default value
   may be insufficient.

 * `CHAIN` is a logical input argument that determines whether a new first
   simplex should be constructed for each interpolation point
   (`CHAIN=.FALSE.`), or whether the simplex walks should be "daisy-chained."
   By default, `CHAIN=.FALSE.` Setting `CHAIN=.TRUE.` is generally not
   recommended, unless the size of the triangulation is relatively small
   or the interpolation points are known to be tightly clustered.

 * `EXACT` is a logical input argument that determines whether the exact
   diameter should be computed and whether a check for duplicate data
   points should be performed in advance. These checks are $O(N^2 D)$ time
   complexity, while `DELAUNAYSPARSE` tends toward $O(N D^4)$ on average.
   By default, `EXACT=.TRUE.` and the exact diameter is computed and an error
   is returned whenever `PTS` contains duplicate values up to the precision
   `EPS`. When `EXACT=.FALSE.`, the diameter of `PTS` is approximated by twice
   the distance from the barycenter of `PTS` to the farthest point in `PTS`,
   and no check is done to find the closest pair of points.
   When `EXACT=.TRUE.`, `DELAUNAYSPARSE` could spend over 90% of runtime
   calculating these constants, which are not critical to the `DELAUNAYSPARSE`
   algorithm. In particular, this happens for large values of `N`. However,
   setting `EXACT=.FALSE.` could result in input errors that are difficult
   to identify. It is recommended that users verify the input set `PTS`
   and possibly rescale `PTS` manually while `EXACT=.TRUE.` Then, when
   100% sure that `PTS` is valid, users may choose to set `EXACT=.FALSE.`
   in production runs for large values of N to achieve massive speedups.

 * `PMODE` is an integer specifying the level of parallelism to be exploited.
    - If `PMODE = 1`, then parallelism is exploited at the level of the loop
      over all interpolation points (Level 1 parallelism).
    - If `PMODE = 2`, then parallelism is exploited at the level of the loops
      over data points when constructing/flipping simplices (Level 2
      parallelism).
    - If `PMODE = 3`, then parallelism is exploited at both levels. Note:
      this implies that the total number of threads active at any time could
      be up to `OMP_NUM_THREADS^2`.
   By default, `PMODE` is set to `1` if there is more than 1 interpolation
   point and `2` otherwise.


Subroutines and functions directly referenced from BLAS are
 * `DDOT`,
 * `DGEMV`,
 * `DNRM2`,
 * `DTRSM`,
and from LAPACK are
 * `DGEQP3`,
 * `DGETRF`,
 * `DGETRS`,
 * `DORMQR`.

The BQPD driver
 * `BQPD` is also directly referenced.

`BQPD` and all its dependencies have been flattened into a single file `bqpd.f`.
For a reference to `BQPD`, see
   Annals of Operations Research, 46 : 307--334 (1993).
The module `REAL_PRECISION` from HOMPACK90 (ACM TOMS Algorithm 777) is
used for the real data type. The `REAL_PRECISION` module, `DELAUNAYSPARSEP`,
and `BQPD` and its dependencies comply with the Fortran 2008 standard.  
