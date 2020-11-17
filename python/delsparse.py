# Python wrapper for DELAUNAYSPARSE using C interface.
import os
import ctypes 
import numpy as np

# --------------------------------------------------------------------
#                        CONFIGURATION
# 
fort_compiler = "gfortran"
shared_object_name = "delsparse_clib.so"
path_to_lib = os.path.join(os.curdir, shared_object_name)
compile_options = "-fPIC -shared -O3 -fopenmp -std=legacy"
# ^^ 'fPIC' and 'shared' are required. 'O3' is for speed and 'fopenmp'
#    is necessary for supporting CPU parallelism during execution.
blas_lapack = "-lblas -llapack"
blas_lapack = "blas.f lapack.f"
# ^^ Use a local BLAS and LAPACK if available by commenting the second line
#    above. The included "blas.f" and "lapack.f" are known to cause error 71
#    during extrapolation, but there is no known resolution.
ordered_dependencies = "real_precision.f90 slatec.f delsparse.f90 delsparse_bind_c.f90"
# 
# --------------------------------------------------------------------


# Try to import the existing object. If that fails, recompile and then try.
try:
    delsparse_clib = ctypes.CDLL(path_to_lib)
except:
    # Remove the shared object if it exists, because it is faulty.
    if os.path.exists(shared_object_name):
        os.remove(shared_object_name)
    # Compile a new shared object.
    command = " ".join([fort_compiler, compile_options, blas_lapack,
                        ordered_dependencies, "-o", path_to_lib])
    print("Running command")
    print("  ", command)
    os.system(command)
    # Remove all ".mod" files that were created to reduce clutter.
    all_mods = [f for f in os.listdir(os.curdir) if f[-4:] == ".mod"]
    for m in all_mods: os.remove(m)

# Import the shared object file as a C library with ctypes.
delsparse_clib = ctypes.CDLL(path_to_lib)

def delaunaysparses(d, n, pts, m, q, simps, weights, ierr, interp_in=None, interp_out=None, eps=None, extrap=None, rnorm=None, ibudget=None, chain=None, exact=None):
    '''! This is a serial implementation of an algorithm for efficiently performing
! interpolation in R^D via the Delaunay triangulation. The algorithm is fully
! described and analyzed in
!
! T. H. Chang, L. T. Watson,  T. C.H. Lux, B. Li, L. Xu, A. R. Butt, K. W.
! Cameron, and Y. Hong. 2018. A polynomial time algorithm for multivariate
! interpolation in arbitrary dimension via the Delaunay triangulation. In
! Proceedings of the ACMSE 2018 Conference (ACMSE '18). ACM, New York, NY,
! USA. Article 12, 8 pages.
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
!    codes are:
!
! 00 : Succesful interpolation.
! 01 : Succesful extrapolation (up to the allowed extrapolation distance).
! 02 : This point was outside the allowed extrapolation distance; the
!      corresponding entries in SIMPS and WEIGHTS contain zero values.
!
! 10 : The dimension D must be positive.
! 11 : Too few data points to construct a triangulation (i.e., N < D+1).
! 12 : No interpolation points given (i.e., M < 1).
! 13 : The first dimension of PTS does not agree with the dimension D.
! 14 : The second dimension of PTS does not agree with the number of points N.
! 15 : The first dimension of Q does not agree with the dimension D.
! 16 : The second dimension of Q does not agree with the number of
!      interpolation points M.
! 17 : The first dimension of the output array SIMPS does not match the number
!      of vertices needed for a D-simplex (D+1).
! 18 : The second dimension of the output array SIMPS does not match the
!      number of interpolation points M.
! 19 : The first dimension of the output array WEIGHTS does not match the
!      number of vertices for a a D-simplex (D+1).
! 20 : The second dimension of the output array WEIGHTS does not match the
!      number of interpolation points M.
! 21 : The size of the error array IERR does not match the number of
!      interpolation points M.
! 22 : INTERP_IN cannot be present without INTERP_OUT or vice versa.
! 23 : The first dimension of INTERP_IN does not match the first
!      dimension of INTERP_OUT.
! 24 : The second dimension of INTERP_IN does not match the number of
!      data points PTS.
! 25 : The second dimension of INTERP_OUT does not match the number of
!      interpolation points M.
! 26 : The budget supplied in IBUDGET does not contain a positive
!      integer.
! 27 : The extrapolation distance supplied in EXTRAP cannot be negative.
! 28 : The size of the RNORM output array does not match the number of
!      interpolation points M.
!
! 30 : Two or more points in the data set PTS are too close together with
!      respect to the working precision (EPS), which would result in a
!      numerically degenerate simplex.
! 31 : All the data points in PTS lie in some lower dimensional linear
!      manifold (up to the working precision), and no valid triangulation
!      exists.
! 40 : An error caused DELAUNAYSPARSES to terminate before this value could
!      be computed. Note: The corresponding entries in SIMPS and WEIGHTS may
!      contain garbage values.
!
! 50 : A memory allocation error occurred while allocating the work array
!      WORK.
!
! 60 : The budget was exceeded before the algorithm converged on this
!      value. If the dimension is high, try increasing IBUDGET. This
!      error can also be caused by a working precision EPS that is too
!      small for the conditioning of the problem.
!
! 61 : A value that was judged appropriate later caused LAPACK to encounter a
!      singularity. Try increasing the value of EPS.
!
! 70 : Allocation error for the extrapolation work arrays.
! 71 : The SLATEC subroutine DWNNLS failed to converge during the projection
!      of an extrapolation point onto the convex hull.
! 72 : The SLATEC subroutine DWNNLS has reported a usage error.
!
!      The errors 72, 80--83 should never occur, and likely indicate a
!      compiler bug or hardware failure.
! 80 : The LAPACK subroutine DGEQP3 has reported an illegal value.
! 81 : The LAPACK subroutine DGETRF has reported an illegal value.
! 82 : The LAPACK subroutine DGETRS has reported an illegal value.
! 83 : The LAPACK subroutine DORMQR has reported an illegal value.
!
!
! Optional arguments:
!
! INTERP_IN(1:IR,1:N) contains real valued response vectors for each of
!    the data points in PTS on input. The first dimension of INTERP_IN is
!    inferred to be the dimension of these response vectors, and the
!    second dimension must match N. If present, the response values will
!    be computed for each interpolation point in Q, and stored in INTERP_OUT,
!    which therefore must also be present. If both INTERP_IN and INTERP_OUT
!    are omitted, only the containing simplices and convex combination
!    weights are returned.
!
! INTERP_OUT(1:IR,1:M) contains real valued response vectors for each
!    interpolation point in Q on output. The first dimension of INTERP_OUT
!    must match the first dimension of INTERP_IN, and the second dimension
!    must match M. If present, the response values at each interpolation
!    point are computed as a convex combination of the response values
!    (supplied in INTERP_IN) at the vertices of a Delaunay simplex containing
!    that interpolation point.  Therefore, if INTERP_OUT is present, then
!    INTERP_IN must also be present.  If both are omitted, only the
!    simplices and convex combination weights are returned.
!
! EPS contains the real working precision for the problem on input. By default,
!    EPS is assigned \sqrt{\mu} where \mu denotes the unit roundoff for the
!    machine. In general, any values that differ by less than EPS are judged
!    as equal, and any weights that are greater than -EPS are judged as
!    nonnegative.  EPS cannot take a value less than the default value of
!    \sqrt{\mu}. If any value less than \sqrt{\mu} is supplied, the default
!    value will be used instead automatically.
!
! EXTRAP contains the real maximum extrapolation distance (relative to the
!    diameter of PTS) on input. Interpolation at a point outside the convex
!    hull of PTS is done by projecting that point onto the convex hull, and
!    then doing normal Delaunay interpolation at that projection.
!    Interpolation at any point in Q that is more than EXTRAP * DIAMETER(PTS)
!    units outside the convex hull of PTS will not be done and an error code
!    of 2 will be returned. Note that computing the projection can be
!    expensive. Setting EXTRAP=0 will cause all extrapolation points to be
!    ignored without ever computing a projection. By default, EXTRAP=0.1
!    (extrapolate by up to 10% of the diameter of PTS).
!
! RNORM(1:M) contains the real unscaled projection (2-norm) distances from
!    any projection computations on output. If not present, these distances
!    are still computed for each extrapolation point, but are never returned.
!
! IBUDGET on input contains the integer budget for performing flips while
!    iterating toward the simplex containing each interpolation point in
!    Q. This prevents DELAUNAYSPARSES from falling into an infinite loop when
!    an inappropriate value of EPS is given with respect to the problem
!    conditioning.  By default, IBUDGET=50000. However, for extremely
!    high-dimensional problems and pathological inputs, the default value
!    may be insufficient.
!
! CHAIN is a logical input argument that determines whether a new first
!    simplex should be constructed for each interpolation point
!    (CHAIN=.FALSE.), or whether the simplex walks should be "daisy-chained."
!    By default, CHAIN=.FALSE. Setting CHAIN=.TRUE. is generally not
!    recommended, unless the size of the triangulation is relatively small
!    or the interpolation points are known to be tightly clustered.
!
! EXACT is a logical input argument that determines whether the exact
!    diameter should be computed and whether a check for duplicate data
!    points should be performed in advance. When EXACT=.FALSE., the
!    diameter of PTS is approximated by twice the distance from the
!    barycenter of PTS to the farthest point in PTS, and no check is
!    done to find the closest pair of points, which could result in hard
!    to find bugs later on. When EXACT=.TRUE., the exact diameter is
!    computed and an error is returned whenever PTS contains duplicate
!    values up to the precision EPS. By default EXACT=.TRUE., but setting
!    EXACT=.FALSE. could result in significant speedup when N is large.
!    It is strongly recommended that most users leave EXACT=.TRUE., as
!    setting EXACT=.FALSE. could result in input errors that are difficult
!    to identify. Also, the diameter approximation could be wrong by up to
!    a factor of two.
!
!
! Subroutines and functions directly referenced from BLAS are
!      DDOT, DGEMV, DNRM2, DTRSM,
! and from LAPACK are
!      DGEQP3, DGETRF, DGETRS, DORMQR.
! The SLATEC subroutine DWNNLS is directly referenced. DWNNLS and all its
! SLATEC dependencies have been slightly edited to comply with the Fortran
! 2008 standard, with all print statements and references to stderr being
! commented out. For a reference to DWNNLS, see ACM TOMS Algorithm 587
! (Hanson and Haskell).  The module REAL_PRECISION from HOMPACK90 (ACM TOMS
! Algorithm 777) is used for the real data type. The REAL_PRECISION module,
! DELAUNAYSPARSES, and DWNNLS and its dependencies comply with the Fortran
! 2008 standard.
!
! Primary Author: Tyler H. Chang
! Last Update: March, 2020
!'''

    # Setting up "d"
    d = ctypes.c_int(d)

    # Setting up "n"
    n = ctypes.c_int(n)

    # Setting up "m"
    m = ctypes.c_int(m)

    # Setting up "pts"
    pts_local = np.asarray(pts, dtype=ctypes.c_double)
    pts_dim_1 = ctypes.c_int(pts_local.shape[0])
    pts_dim_2 = ctypes.c_int(pts_local.shape[1])
    
    # Setting up "q"
    q_local = np.asarray(q, dtype=ctypes.c_double)
    q_dim_1 = ctypes.c_int(q_local.shape[0])
    q_dim_2 = ctypes.c_int(q_local.shape[1])
    
    # Setting up "simps"
    simps_local = np.asarray(simps, dtype=ctypes.c_int)
    simps_dim_1 = ctypes.c_int(simps_local.shape[0])
    simps_dim_2 = ctypes.c_int(simps_local.shape[1])
    
    # Setting up "weights"
    weights_local = np.asarray(weights, dtype=ctypes.c_double)
    weights_dim_1 = ctypes.c_int(weights_local.shape[0])
    weights_dim_2 = ctypes.c_int(weights_local.shape[1])
    
    # Setting up "ierr"
    ierr_local = np.asarray(ierr, dtype=ctypes.c_int)
    ierr_dim_1 = ctypes.c_int(ierr_local.shape[0])
    
    # Setting up "interp_in"
    interp_in_present = ctypes.c_bool(True)
    interp_in_dim_1 = ctypes.c_int(0)
    interp_in_dim_2 = ctypes.c_int(0)
    if (interp_in is None):
        interp_in_present = ctypes.c_bool(False)
        interp_in = np.zeros(shape=(1,1), dtype=ctypes.c_double, order='F')
    elif (type(interp_in) == bool) and (interp_in):
        interp_in = np.zeros(shape=(1,1), dtype=ctypes.c_double, order='F')
        interp_in_dim_1 = ctypes.c_int(interp_in.shape[0])
        interp_in_dim_2 = ctypes.c_int(interp_in.shape[1])
    elif (not np.asarray(interp_in).flags.f_contiguous):
        raise(Exception("The numpy array given as argument 'interp_in' was not f_contiguous."))
    else:
        interp_in_dim_1 = ctypes.c_int(interp_in.shape[0])
        interp_in_dim_2 = ctypes.c_int(interp_in.shape[1])
    interp_in_local = np.asarray(interp_in, dtype=ctypes.c_double)
    
    # Setting up "interp_out"
    interp_out_present = ctypes.c_bool(True)
    interp_out_dim_1 = ctypes.c_int(0)
    interp_out_dim_2 = ctypes.c_int(0)
    if (interp_out is None):
        interp_out_present = ctypes.c_bool(False)
        interp_out = np.zeros(shape=(1,1), dtype=ctypes.c_double, order='F')
    elif (type(interp_out) == bool) and (interp_out):
        interp_out = np.zeros(shape=(1,1), dtype=ctypes.c_double, order='F')
        interp_out_dim_1 = ctypes.c_int(interp_out.shape[0])
        interp_out_dim_2 = ctypes.c_int(interp_out.shape[1])
    elif (not np.asarray(interp_out).flags.f_contiguous):
        raise(Exception("The numpy array given as argument 'interp_out' was not f_contiguous."))
    else:
        interp_out_dim_1 = ctypes.c_int(interp_out.shape[0])
        interp_out_dim_2 = ctypes.c_int(interp_out.shape[1])
    interp_out_local = np.asarray(interp_out, dtype=ctypes.c_double)
    
    # Setting up "eps"
    eps_present = ctypes.c_bool(True)
    if (eps is None):
        eps_present = ctypes.c_bool(False)
        eps = 1
    eps_local = ctypes.c_double(eps)
    
    # Setting up "extrap"
    extrap_present = ctypes.c_bool(True)
    if (extrap is None):
        extrap_present = ctypes.c_bool(False)
        extrap = 1
    extrap_local = ctypes.c_double(extrap)
    
    # Setting up "rnorm"
    rnorm_present = ctypes.c_bool(True)
    rnorm_dim_1 = ctypes.c_int(0)
    if (rnorm is None):
        rnorm_present = ctypes.c_bool(False)
        rnorm = np.zeros(shape=(1), dtype=ctypes.c_double, order='F')
    elif (type(rnorm) == bool) and (rnorm):
        rnorm = np.zeros(shape=(1), dtype=ctypes.c_double, order='F')
        rnorm_dim_1 = ctypes.c_int(rnorm.shape[0])
    elif (not np.asarray(rnorm).flags.f_contiguous):
        raise(Exception("The numpy array given as argument 'rnorm' was not f_contiguous."))
    else:
        rnorm_dim_1 = ctypes.c_int(rnorm.shape[0])
    rnorm_local = np.asarray(rnorm, dtype=ctypes.c_double)
    
    # Setting up "ibudget"
    ibudget_present = ctypes.c_bool(True)
    if (ibudget is None):
        ibudget_present = ctypes.c_bool(False)
        ibudget = 1
    ibudget_local = ctypes.c_int(ibudget)
    
    # Setting up "chain"
    chain_present = ctypes.c_bool(True)
    if (chain is None):
        chain_present = ctypes.c_bool(False)
        chain = 1
    chain_local = ctypes.c_bool(chain)
    
    # Setting up "exact"
    exact_present = ctypes.c_bool(True)
    if (exact is None):
        exact_present = ctypes.c_bool(False)
        exact = 1
    exact_local = ctypes.c_bool(exact)

    # Call C-accessible Fortran wrapper.
    delsparse_clib.c_delaunaysparses(ctypes.byref(d), ctypes.byref(n), ctypes.byref(pts_dim_1), ctypes.byref(pts_dim_2), ctypes.c_void_p(pts_local.ctypes.data), ctypes.byref(m), ctypes.byref(q_dim_1), ctypes.byref(q_dim_2), ctypes.c_void_p(q_local.ctypes.data), ctypes.byref(simps_dim_1), ctypes.byref(simps_dim_2), ctypes.c_void_p(simps_local.ctypes.data), ctypes.byref(weights_dim_1), ctypes.byref(weights_dim_2), ctypes.c_void_p(weights_local.ctypes.data), ctypes.byref(ierr_dim_1), ctypes.c_void_p(ierr_local.ctypes.data), ctypes.byref(interp_in_present), ctypes.byref(interp_in_dim_1), ctypes.byref(interp_in_dim_2), ctypes.c_void_p(interp_in_local.ctypes.data), ctypes.byref(interp_out_present), ctypes.byref(interp_out_dim_1), ctypes.byref(interp_out_dim_2), ctypes.c_void_p(interp_out_local.ctypes.data), ctypes.byref(eps_present), ctypes.byref(eps_local), ctypes.byref(extrap_present), ctypes.byref(extrap_local), ctypes.byref(rnorm_present), ctypes.byref(rnorm_dim_1), ctypes.c_void_p(rnorm_local.ctypes.data), ctypes.byref(ibudget_present), ctypes.byref(ibudget_local), ctypes.byref(chain_present), ctypes.byref(chain_local), ctypes.byref(exact_present), ctypes.byref(exact_local))

    # Return final results, 'INTENT(OUT)' arguments only.
    return np.asarray(pts_local), np.asarray(q_local), np.asarray(simps_local), np.asarray(weights_local), np.asarray(ierr_local), (np.asarray(interp_out_local) if interp_out_present else None), (np.asarray(rnorm_local) if rnorm_present else None)


# ----------------------------------------------
# Wrapper for the Fortran subroutine DELAUNAYSPARSEP

def delaunaysparsep(d, n, pts, m, q, simps, weights, ierr, interp_in=None, interp_out=None, eps=None, extrap=None, rnorm=None, ibudget=None, chain=None, exact=None, pmode=None):
    '''! This is a parallel implementation of an algorithm for efficiently performing
! interpolation in R^D via the Delaunay triangulation. The algorithm is fully
! described and analyzed in
!
! T. H. Chang, L. T. Watson,  T. C.H. Lux, B. Li, L. Xu, A. R. Butt, K. W.
! Cameron, and Y. Hong. 2018. A polynomial time algorithm for multivariate
! interpolation in arbitrary dimension via the Delaunay triangulation. In
! Proceedings of the ACMSE 2018 Conference (ACMSE '18). ACM, New York, NY,
! USA. Article 12, 8 pages.
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
!    codes are:
!
! 00 : Succesful interpolation.
! 01 : Succesful extrapolation (up to the allowed extrapolation distance).
! 02 : This point was outside the allowed extrapolation distance; the
!      corresponding entries in SIMPS and WEIGHTS contain zero values.
!
! 10 : The dimension D must be positive.
! 11 : Too few data points to construct a triangulation (i.e., N < D+1).
! 12 : No interpolation points given (i.e., M < 1).
! 13 : The first dimension of PTS does not agree with the dimension D.
! 14 : The second dimension of PTS does not agree with the number of points N.
! 15 : The first dimension of Q does not agree with the dimension D.
! 16 : The second dimension of Q does not agree with the number of
!      interpolation points M.
! 17 : The first dimension of the output array SIMPS does not match the number
!      of vertices needed for a D-simplex (D+1).
! 18 : The second dimension of the output array SIMPS does not match the
!      number of interpolation points M.
! 19 : The first dimension of the output array WEIGHTS does not match the
!      number of vertices for a a D-simplex (D+1).
! 20 : The second dimension of the output array WEIGHTS does not match the
!      number of interpolation points M.
! 21 : The size of the error array IERR does not match the number of
!      interpolation points M.
! 22 : INTERP_IN cannot be present without INTERP_OUT or vice versa.
! 23 : The first dimension of INTERP_IN does not match the first
!      dimension of INTERP_OUT.
! 24 : The second dimension of INTERP_IN does not match the number of
!      data points PTS.
! 25 : The second dimension of INTERP_OUT does not match the number of
!      interpolation points M.
! 26 : The budget supplied in IBUDGET does not contain a positive
!      integer.
! 27 : The extrapolation distance supplied in EXTRAP cannot be negative.
! 28 : The size of the RNORM output array does not match the number of
!      interpolation points M.
!
! 30 : Two or more points in the data set PTS are too close together with
!      respect to the working precision (EPS), which would result in a
!      numerically degenerate simplex.
! 31 : All the data points in PTS lie in some lower dimensional linear
!      manifold (up to the working precision), and no valid triangulation
!      exists.
! 40 : An error caused DELAUNAYSPARSEP to terminate before this value could
!      be computed. Note: The corresponding entries in SIMPS and WEIGHTS may
!      contain garbage values.
!
! 50 : A memory allocation error occurred while allocating the work array
!      WORK.
!
! 60 : The budget was exceeded before the algorithm converged on this
!      value. If the dimension is high, try increasing IBUDGET. This
!      error can also be caused by a working precision EPS that is too
!      small for the conditioning of the problem.
!
! 61 : A value that was judged appropriate later caused LAPACK to encounter a
!      singularity. Try increasing the value of EPS.
!
! 70 : Allocation error for the extrapolation work arrays.
! 71 : The SLATEC subroutine DWNNLS failed to converge during the projection
!      of an extrapolation point onto the convex hull.
! 72 : The SLATEC subroutine DWNNLS has reported a usage error.
!
!      The errors 72, 80--83 should never occur, and likely indicate a
!      compiler bug or hardware failure.
! 80 : The LAPACK subroutine DGEQP3 has reported an illegal value.
! 81 : The LAPACK subroutine DGETRF has reported an illegal value.
! 82 : The LAPACK subroutine DGETRS has reported an illegal value.
! 83 : The LAPACK subroutine DORMQR has reported an illegal value.
!
! 90 : The value of PMODE is not valid.
!
!
! Optional arguments:
!
! INTERP_IN(1:IR,1:N) contains real valued response vectors for each of
!    the data points in PTS on input. The first dimension of INTERP_IN is
!    inferred to be the dimension of these response vectors, and the
!    second dimension must match N. If present, the response values will
!    be computed for each interpolation point in Q, and stored in INTERP_OUT,
!    which therefore must also be present. If both INTERP_IN and INTERP_OUT
!    are omitted, only the containing simplices and convex combination
!    weights are returned.
!
! INTERP_OUT(1:IR,1:M) contains real valued response vectors for each
!    interpolation point in Q on output. The first dimension of INTERP_OU
!    must match the first dimension of INTERP_IN, and the second dimension
!    must match M. If present, the response values at each interpolation
!    point are computed as a convex combination of the response values
!    (supplied in INTERP_IN) at the vertices of a Delaunay simplex containing
!    that interpolation point.  Therefore, if INTERP_OUT is present, then
!    INTERP_IN must also be present.  If both are omitted, only the
!    simplices and convex combination weights are returned.
!
! EPS contains the real working precision for the problem on input. By
!    default, EPS is assigned \sqrt{\mu} where \mu denotes the unit roundoff
!    for the machine. In general, any values that differ by less than EPS
!    are judged as equal, and any weights that are greater than -EPS are
!    judged as nonnegative.  EPS cannot take a value less than the default
!    value of \sqrt{\mu}. If any value less than \sqrt{\mu} is supplied,
!    the default value will be used instead automatically.
!
! EXTRAP contains the real maximum extrapolation distance (relative to the
!    diameter of PTS) on input. Interpolation at a point outside the convex
!    hull of PTS is done by projecting that point onto the convex hull, and
!    then doing normal Delaunay interpolation at that projection.
!    Interpolation at any point in Q that is more than EXTRAP * DIAMETER(PTS)
!    units outside the convex hull of PTS will not be done and an error code
!    of 2 will be returned. Note that computing the projection can be
!    expensive. Setting EXTRAP=0 will cause all extrapolation points to be
!    ignored without ever computing a projection. By default, EXTRAP=0.1
!    (extrapolate by up to 10% of the diameter of PTS).
!
! RNORM(1:M) contains the real unscaled projection (2-norm) distances from
!    any projection computations on output. If not present, these distances
!    are still computed for each extrapolation point, but are never returned.
!
! IBUDGET on input contains the integer budget for performing flips while
!    iterating toward the simplex containing each interpolation point in Q.
!    This prevents DELAUNAYSPARSEP from falling into an infinite loop when
!    an inappropriate value of EPS is given with respect to the problem
!    conditioning.  By default, IBUDGET=50000. However, for extremely
!    high-dimensional problems and pathological inputs, the default value
!    may be insufficient.
!
! CHAIN is a logical input argument that determines whether a new first
!    simplex should be constructed for each interpolation point
!    (CHAIN=.FALSE.), or whether the simplex walks should be "daisy-chained."
!    By default, CHAIN=.FALSE. Setting CHAIN=.TRUE. is generally not
!    recommended, unless the size of the triangulation is relatively small
!    or the interpolation points are known to be tightly clustered.
!
! EXACT is a logical input argument that determines whether the exact
!    diameter should be computed and whether a check for duplicate data
!    points should be performed in advance. When EXACT=.FALSE., the
!    diameter of PTS is approximated by twice the distance from the
!    barycenter of PTS to the farthest point in PTS, and no check is
!    done to find the closest pair of points, which could result in hard
!    to find bugs later on. When EXACT=.TRUE., the exact diameter is
!    computed and an error is returned whenever PTS contains duplicate
!    values up to the precision EPS. By default EXACT=.TRUE., but setting
!    EXACT=.FALSE. could result in significant speedup when N is large.
!    It is strongly recommended that most users leave EXACT=.TRUE., as
!    setting EXACT=.FALSE. could result in input errors that are difficult
!    to identify. Also, the diameter approximation could be wrong by up to
!    a factor of two.
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
!    By default, PMODE is set to 1 if there is more than 1 interpolation
!    point and 2 otherwise.
!
!
! Subroutines and functions directly referenced from BLAS are
!      DDOT, DGEMV, DNRM2, DTRSM,
! and from LAPACK are
!      DGEQP3, DGETRF, DGETRS, DORMQR.
! The SLATEC subroutine DWNNLS is directly referenced. DWNNLS and all its
! SLATEC dependencies have been slightly edited to comply with the Fortran
! 2008 standard, with all print statements and references to stderr being
! commented out. For a reference to DWNNLS, see ACM TOMS Algorithm 587
! (Hanson and Haskell).  The module REAL_PRECISION from HOMPACK90 (ACM TOMS
! Algorithm 777) is used for the real data type. The REAL_PRECISION module,
! DELAUNAYSPARSEP, and DWNNLS and its dependencies comply with the Fortran
! 2008 standard.
!
! Primary Author: Tyler H. Chang
! Last Update: March, 2020
!'''

    # Setting up "d"
    d = ctypes.c_int(d)

    # Setting up "n"
    n = ctypes.c_int(n)

    # Setting up "m"
    m = ctypes.c_int(m)

    # Setting up "pts"
    pts_local = np.asarray(pts, dtype=ctypes.c_double)
    pts_dim_1 = ctypes.c_int(pts_local.shape[0])
    pts_dim_2 = ctypes.c_int(pts_local.shape[1])
    
    # Setting up "q"
    q_local = np.asarray(q, dtype=ctypes.c_double)
    q_dim_1 = ctypes.c_int(q_local.shape[0])
    q_dim_2 = ctypes.c_int(q_local.shape[1])
    
    # Setting up "simps"
    simps_local = np.asarray(simps, dtype=ctypes.c_int)
    simps_dim_1 = ctypes.c_int(simps_local.shape[0])
    simps_dim_2 = ctypes.c_int(simps_local.shape[1])
    
    # Setting up "weights"
    weights_local = np.asarray(weights, dtype=ctypes.c_double)
    weights_dim_1 = ctypes.c_int(weights_local.shape[0])
    weights_dim_2 = ctypes.c_int(weights_local.shape[1])
    
    # Setting up "ierr"
    ierr_local = np.asarray(ierr, dtype=ctypes.c_int)
    ierr_dim_1 = ctypes.c_int(ierr_local.shape[0])
    
    # Setting up "interp_in"
    interp_in_present = ctypes.c_bool(True)
    interp_in_dim_1 = ctypes.c_int(0)
    interp_in_dim_2 = ctypes.c_int(0)
    if (interp_in is None):
        interp_in_present = ctypes.c_bool(False)
        interp_in = np.zeros(shape=(1,1), dtype=ctypes.c_double, order='F')
    elif (type(interp_in) == bool) and (interp_in):
        interp_in = np.zeros(shape=(1,1), dtype=ctypes.c_double, order='F')
        interp_in_dim_1 = ctypes.c_int(interp_in.shape[0])
        interp_in_dim_2 = ctypes.c_int(interp_in.shape[1])
    elif (not np.asarray(interp_in).flags.f_contiguous):
        raise(Exception("The numpy array given as argument 'interp_in' was not f_contiguous."))
    else:
        interp_in_dim_1 = ctypes.c_int(interp_in.shape[0])
        interp_in_dim_2 = ctypes.c_int(interp_in.shape[1])
    interp_in_local = np.asarray(interp_in, dtype=ctypes.c_double)
    
    # Setting up "interp_out"
    interp_out_present = ctypes.c_bool(True)
    interp_out_dim_1 = ctypes.c_int(0)
    interp_out_dim_2 = ctypes.c_int(0)
    if (interp_out is None):
        interp_out_present = ctypes.c_bool(False)
        interp_out = np.zeros(shape=(1,1), dtype=ctypes.c_double, order='F')
    elif (type(interp_out) == bool) and (interp_out):
        interp_out = np.zeros(shape=(1,1), dtype=ctypes.c_double, order='F')
        interp_out_dim_1 = ctypes.c_int(interp_out.shape[0])
        interp_out_dim_2 = ctypes.c_int(interp_out.shape[1])
    elif (not np.asarray(interp_out).flags.f_contiguous):
        raise(Exception("The numpy array given as argument 'interp_out' was not f_contiguous."))
    else:
        interp_out_dim_1 = ctypes.c_int(interp_out.shape[0])
        interp_out_dim_2 = ctypes.c_int(interp_out.shape[1])
    interp_out_local = np.asarray(interp_out, dtype=ctypes.c_double)
    
    # Setting up "eps"
    eps_present = ctypes.c_bool(True)
    if (eps is None):
        eps_present = ctypes.c_bool(False)
        eps = 1
    eps_local = ctypes.c_double(eps)
    
    # Setting up "extrap"
    extrap_present = ctypes.c_bool(True)
    if (extrap is None):
        extrap_present = ctypes.c_bool(False)
        extrap = 1
    extrap_local = ctypes.c_double(extrap)
    
    # Setting up "rnorm"
    rnorm_present = ctypes.c_bool(True)
    rnorm_dim_1 = ctypes.c_int(0)
    if (rnorm is None):
        rnorm_present = ctypes.c_bool(False)
        rnorm = np.zeros(shape=(1), dtype=ctypes.c_double, order='F')
    elif (type(rnorm) == bool) and (rnorm):
        rnorm = np.zeros(shape=(1), dtype=ctypes.c_double, order='F')
        rnorm_dim_1 = rnorm.shape[0]
    elif (not np.asarray(rnorm).flags.f_contiguous):
        raise(Exception("The numpy array given as argument 'rnorm' was not f_contiguous."))
    else:
        rnorm_dim_1 = ctypes.c_int(rnorm.shape[0])
    rnorm_local = np.asarray(rnorm, dtype=ctypes.c_double)
    
    # Setting up "ibudget"
    ibudget_present = ctypes.c_bool(True)
    if (ibudget is None):
        ibudget_present = ctypes.c_bool(False)
        ibudget = 1
    ibudget_local = ctypes.c_int(ibudget)
    
    # Setting up "chain"
    chain_present = ctypes.c_bool(True)
    if (chain is None):
        chain_present = ctypes.c_bool(False)
        chain = 1
    chain_local = ctypes.c_bool(chain)
    
    # Setting up "exact"
    exact_present = ctypes.c_bool(True)
    if (exact is None):
        exact_present = ctypes.c_bool(False)
        exact = 1
    exact_local = ctypes.c_bool(exact)

    # Setting up "pmode"
    pmode_present = ctypes.c_bool(True)
    if (pmode is None):
        pmode_present = ctypes.c_bool(False)
        pmode = 1
    pmode_local = ctypes.c_int(pmode)

    # Call C-accessible Fortran wrapper.
    delsparse_clib.c_delaunaysparsep(ctypes.byref(d), ctypes.byref(n), ctypes.byref(pts_dim_1), ctypes.byref(pts_dim_2), ctypes.c_void_p(pts_local.ctypes.data), ctypes.byref(m), ctypes.byref(q_dim_1), ctypes.byref(q_dim_2), ctypes.c_void_p(q_local.ctypes.data), ctypes.byref(simps_dim_1), ctypes.byref(simps_dim_2), ctypes.c_void_p(simps_local.ctypes.data), ctypes.byref(weights_dim_1), ctypes.byref(weights_dim_2), ctypes.c_void_p(weights_local.ctypes.data), ctypes.byref(ierr_dim_1), ctypes.c_void_p(ierr_local.ctypes.data), ctypes.byref(interp_in_present), ctypes.byref(interp_in_dim_1), ctypes.byref(interp_in_dim_2), ctypes.c_void_p(interp_in_local.ctypes.data), ctypes.byref(interp_out_present), ctypes.byref(interp_out_dim_1), ctypes.byref(interp_out_dim_2), ctypes.c_void_p(interp_out_local.ctypes.data), ctypes.byref(eps_present), ctypes.byref(eps_local), ctypes.byref(extrap_present), ctypes.byref(extrap_local), ctypes.byref(rnorm_present), ctypes.byref(rnorm_dim_1), ctypes.c_void_p(rnorm_local.ctypes.data), ctypes.byref(ibudget_present), ctypes.byref(ibudget_local), ctypes.byref(chain_present), ctypes.byref(chain_local), ctypes.byref(exact_present), ctypes.byref(exact_local), ctypes.byref(pmode_present), ctypes.byref(pmode_local))

    # Return final results, 'INTENT(OUT)' arguments only.
    return np.asarray(pts_local), np.asarray(q_local), np.asarray(simps_local), np.asarray(weights_local), np.asarray(ierr_local), (np.asarray(interp_out_local) if interp_out_present else None), (np.asarray(rnorm_local) if rnorm_present else None)

