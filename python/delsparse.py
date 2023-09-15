'''This Python code is an automatically generated wrapper
for Fortran code made by 'fmodpy'. The original documentation
for the Fortran source code follows.


'''

import os
import ctypes
import platform
import numpy

# --------------------------------------------------------------------
#               CONFIGURATION
# 
_verbose = True
_fort_compiler = "gfortran"
_shared_object_name = "delsparse." + platform.machine() + ".so"
_this_directory = os.path.dirname(os.path.abspath(__file__))
_src_directory = os.path.join(_this_directory, "delsparse_src")
_path_to_lib = os.path.join(_this_directory, _shared_object_name)
_compile_options = ['-fPIC', '-shared', '-O3']
_ordered_dependencies = ['bqpd_min/bqpd.f', 'lapack.f', 'blas.f', 'delsparse.f90', 'delsparse_bind_c.f90']
_ordered_dependencies = [f"{_src_directory}/{di}" for di in _ordered_dependencies]
_symbol_files = []# 
# --------------------------------------------------------------------
#               AUTO-COMPILING
#
# Try to import the prerequisite symbols for the compiled code.
for _ in _symbol_files:
    _ = ctypes.CDLL(os.path.join(_this_directory, _), mode=ctypes.RTLD_GLOBAL)
# Try to import the existing object. If that fails, recompile and then try.
try:
    # Check to see if the source files have been modified and a recompilation is needed.
    if (max(max([0]+[os.path.getmtime(os.path.realpath(os.path.join(_this_directory,_))) for _ in _symbol_files]),
            max([0]+[os.path.getmtime(os.path.realpath(os.path.join(_this_directory,_))) for _ in _ordered_dependencies]))
        > os.path.getmtime(_path_to_lib)):
        print()
        print("WARNING: Recompiling because the modification time of a source file is newer than the library.", flush=True)
        print()
        if os.path.exists(_path_to_lib):
            os.remove(_path_to_lib)
        raise NotImplementedError(f"The newest library code has not been compiled.")
    # Import the library.
    clib = ctypes.CDLL(_path_to_lib)
except:
    # Remove the shared object if it exists, because it is faulty.
    if os.path.exists(_shared_object_name):
        os.remove(_shared_object_name)
    # Compile a new shared object.
    _command = " ".join([_fort_compiler] + _compile_options + ["-o", _shared_object_name] + _ordered_dependencies)
    if _verbose:
        print("Running system command with arguments")
        print("  ", _command)
    # Run the compilation command.
    import subprocess
    subprocess.run(_command, shell=True, cwd=_this_directory)
    # Import the shared object file as a C library with ctypes.
    clib = ctypes.CDLL(_path_to_lib)
# --------------------------------------------------------------------


# ----------------------------------------------
# Wrapper for the Fortran subroutine DELAUNAYSPARSES

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
! 70 : Allocation error for the extrapolation work arrays. Note that
!      extrapolation has a higher memory overhead than interpolation for
!      the current version.
! 7x : BQPD has reported an error while computing the projection. See the
!      comment block for the PROJECT subroutine for more details.
!
!      The errors 80--83 should never occur, and likely indicate a
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
! The BQPD solver is also used. For more information, see
!      Annals of Operations Research, 46 : 307--334 (1993).
! The module REAL_PRECISION from HOMPACK90 (ACM TOMS Algorithm 777) is used
! for the real data type. The REAL_PRECISION module, DELAUNAYSPARSES, and
! BQPD and its dependencies comply with the Fortran 2008 standard.
!
! Primary Author: Tyler H. Chang (THC)
!
! Version history:
!
! Version 2: Sep 2023 (THC et al.) -- replaced DWNNLS with BQPD (ACM TOMS Rmk)
! Version 1: Mar 2020 (THC et al.) -- original version (ACM TOMS Alg 1012)'''
    
    # Setting up "d"
    if (type(d) is not ctypes.c_int): d = ctypes.c_int(d)
    
    # Setting up "n"
    if (type(n) is not ctypes.c_int): n = ctypes.c_int(n)
    
    # Setting up "pts"
    if ((not issubclass(type(pts), numpy.ndarray)) or
        (not numpy.asarray(pts).flags.f_contiguous) or
        (not (pts.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'pts' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        pts = numpy.asarray(pts, dtype=ctypes.c_double, order='F')
    pts_dim_1 = ctypes.c_long(pts.shape[0])
    pts_dim_2 = ctypes.c_long(pts.shape[1])
    
    # Setting up "m"
    if (type(m) is not ctypes.c_int): m = ctypes.c_int(m)
    
    # Setting up "q"
    if ((not issubclass(type(q), numpy.ndarray)) or
        (not numpy.asarray(q).flags.f_contiguous) or
        (not (q.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'q' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        q = numpy.asarray(q, dtype=ctypes.c_double, order='F')
    q_dim_1 = ctypes.c_long(q.shape[0])
    q_dim_2 = ctypes.c_long(q.shape[1])
    
    # Setting up "simps"
    if ((not issubclass(type(simps), numpy.ndarray)) or
        (not numpy.asarray(simps).flags.f_contiguous) or
        (not (simps.dtype == numpy.dtype(ctypes.c_int)))):
        import warnings
        warnings.warn("The provided argument 'simps' was not an f_contiguous NumPy array of type 'ctypes.c_int' (or equivalent). Automatically converting (probably creating a full copy).")
        simps = numpy.asarray(simps, dtype=ctypes.c_int, order='F')
    simps_dim_1 = ctypes.c_long(simps.shape[0])
    simps_dim_2 = ctypes.c_long(simps.shape[1])
    
    # Setting up "weights"
    if ((not issubclass(type(weights), numpy.ndarray)) or
        (not numpy.asarray(weights).flags.f_contiguous) or
        (not (weights.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'weights' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        weights = numpy.asarray(weights, dtype=ctypes.c_double, order='F')
    weights_dim_1 = ctypes.c_long(weights.shape[0])
    weights_dim_2 = ctypes.c_long(weights.shape[1])
    
    # Setting up "ierr"
    if ((not issubclass(type(ierr), numpy.ndarray)) or
        (not numpy.asarray(ierr).flags.f_contiguous) or
        (not (ierr.dtype == numpy.dtype(ctypes.c_int)))):
        import warnings
        warnings.warn("The provided argument 'ierr' was not an f_contiguous NumPy array of type 'ctypes.c_int' (or equivalent). Automatically converting (probably creating a full copy).")
        ierr = numpy.asarray(ierr, dtype=ctypes.c_int, order='F')
    ierr_dim_1 = ctypes.c_long(ierr.shape[0])
    
    # Setting up "interp_in"
    interp_in_present = ctypes.c_bool(True)
    if (interp_in is None):
        interp_in_present = ctypes.c_bool(False)
        interp_in = numpy.zeros(shape=(1,1), dtype=ctypes.c_double, order='F')
    elif (type(interp_in) == bool) and (interp_in):
        interp_in = numpy.zeros(shape=(1,1), dtype=ctypes.c_double, order='F')
    elif ((not issubclass(type(interp_in), numpy.ndarray)) or
          (not numpy.asarray(interp_in).flags.f_contiguous) or
          (not (interp_in.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'interp_in' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        interp_in = numpy.asarray(interp_in, dtype=ctypes.c_double, order='F')
    if (interp_in_present):
        interp_in_dim_1 = ctypes.c_long(interp_in.shape[0])
        interp_in_dim_2 = ctypes.c_long(interp_in.shape[1])
    else:
        interp_in_dim_1 = ctypes.c_long()
        interp_in_dim_2 = ctypes.c_long()
    
    # Setting up "interp_out"
    interp_out_present = ctypes.c_bool(True)
    if (interp_out is None):
        interp_out_present = ctypes.c_bool(False)
        interp_out = numpy.zeros(shape=(1,1), dtype=ctypes.c_double, order='F')
    elif (type(interp_out) == bool) and (interp_out):
        interp_out = numpy.zeros(shape=(1,1), dtype=ctypes.c_double, order='F')
    elif ((not issubclass(type(interp_out), numpy.ndarray)) or
          (not numpy.asarray(interp_out).flags.f_contiguous) or
          (not (interp_out.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'interp_out' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        interp_out = numpy.asarray(interp_out, dtype=ctypes.c_double, order='F')
    if (interp_out_present):
        interp_out_dim_1 = ctypes.c_long(interp_out.shape[0])
        interp_out_dim_2 = ctypes.c_long(interp_out.shape[1])
    else:
        interp_out_dim_1 = ctypes.c_long()
        interp_out_dim_2 = ctypes.c_long()
    
    # Setting up "eps"
    eps_present = ctypes.c_bool(True)
    if (eps is None):
        eps_present = ctypes.c_bool(False)
        eps = ctypes.c_double()
    else:
        eps = ctypes.c_double(eps)
    if (type(eps) is not ctypes.c_double): eps = ctypes.c_double(eps)
    
    # Setting up "extrap"
    extrap_present = ctypes.c_bool(True)
    if (extrap is None):
        extrap_present = ctypes.c_bool(False)
        extrap = ctypes.c_double()
    else:
        extrap = ctypes.c_double(extrap)
    if (type(extrap) is not ctypes.c_double): extrap = ctypes.c_double(extrap)
    
    # Setting up "rnorm"
    rnorm_present = ctypes.c_bool(True)
    if (rnorm is None):
        rnorm_present = ctypes.c_bool(False)
        rnorm = numpy.zeros(shape=(1), dtype=ctypes.c_double, order='F')
    elif (type(rnorm) == bool) and (rnorm):
        rnorm = numpy.zeros(shape=(1), dtype=ctypes.c_double, order='F')
    elif ((not issubclass(type(rnorm), numpy.ndarray)) or
          (not numpy.asarray(rnorm).flags.f_contiguous) or
          (not (rnorm.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'rnorm' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        rnorm = numpy.asarray(rnorm, dtype=ctypes.c_double, order='F')
    if (rnorm_present):
        rnorm_dim_1 = ctypes.c_long(rnorm.shape[0])
    else:
        rnorm_dim_1 = ctypes.c_long()
    
    # Setting up "ibudget"
    ibudget_present = ctypes.c_bool(True)
    if (ibudget is None):
        ibudget_present = ctypes.c_bool(False)
        ibudget = ctypes.c_int()
    else:
        ibudget = ctypes.c_int(ibudget)
    if (type(ibudget) is not ctypes.c_int): ibudget = ctypes.c_int(ibudget)
    
    # Setting up "chain"
    chain_present = ctypes.c_bool(True)
    if (chain is None):
        chain_present = ctypes.c_bool(False)
        chain = ctypes.c_int()
    else:
        chain = ctypes.c_int(chain)
    if (type(chain) is not ctypes.c_int): chain = ctypes.c_int(chain)
    
    # Setting up "exact"
    exact_present = ctypes.c_bool(True)
    if (exact is None):
        exact_present = ctypes.c_bool(False)
        exact = ctypes.c_int()
    else:
        exact = ctypes.c_int(exact)
    if (type(exact) is not ctypes.c_int): exact = ctypes.c_int(exact)

    # Call C-accessible Fortran wrapper.
    clib.c_delaunaysparses(ctypes.byref(d), ctypes.byref(n), ctypes.byref(pts_dim_1), ctypes.byref(pts_dim_2), ctypes.c_void_p(pts.ctypes.data), ctypes.byref(m), ctypes.byref(q_dim_1), ctypes.byref(q_dim_2), ctypes.c_void_p(q.ctypes.data), ctypes.byref(simps_dim_1), ctypes.byref(simps_dim_2), ctypes.c_void_p(simps.ctypes.data), ctypes.byref(weights_dim_1), ctypes.byref(weights_dim_2), ctypes.c_void_p(weights.ctypes.data), ctypes.byref(ierr_dim_1), ctypes.c_void_p(ierr.ctypes.data), ctypes.byref(interp_in_present), ctypes.byref(interp_in_dim_1), ctypes.byref(interp_in_dim_2), ctypes.c_void_p(interp_in.ctypes.data), ctypes.byref(interp_out_present), ctypes.byref(interp_out_dim_1), ctypes.byref(interp_out_dim_2), ctypes.c_void_p(interp_out.ctypes.data), ctypes.byref(eps_present), ctypes.byref(eps), ctypes.byref(extrap_present), ctypes.byref(extrap), ctypes.byref(rnorm_present), ctypes.byref(rnorm_dim_1), ctypes.c_void_p(rnorm.ctypes.data), ctypes.byref(ibudget_present), ctypes.byref(ibudget), ctypes.byref(chain_present), ctypes.byref(chain), ctypes.byref(exact_present), ctypes.byref(exact))

    # Return final results, 'INTENT(OUT)' arguments only.
    return pts, q, simps, weights, ierr, (interp_out if interp_out_present else None), (rnorm if rnorm_present else None)


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
! 70 : Allocation error for the extrapolation work arrays. Note that
!      extrapolation has a higher memory overhead than interpolation for
!      the current version.
! 7x : BQPD has reported an error while computing the projection. See the
!      comment block for the PROJECT subroutine for more details.
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
! The BQPD solver is also used. For more information, see
!      Annals of Operations Research, 46 : 307--334 (1993).
! The module REAL_PRECISION from HOMPACK90 (ACM TOMS Algorithm 777) is used
! for the real data type. The REAL_PRECISION module, DELAUNAYSPARSES, and
! BQPD and its dependencies comply with the Fortran 2008 standard.
!
! Primary Author: Tyler H. Chang (THC)
!
! Version history:
!
! Version 2: Sep 2023 (THC et al.) -- replaced DWNNLS with BQPD (ACM TOMS Rmk)
! Version 1: Mar 2020 (THC et al.) -- original version (ACM TOMS Alg 1012)'''
    
    # Setting up "d"
    if (type(d) is not ctypes.c_int): d = ctypes.c_int(d)
    
    # Setting up "n"
    if (type(n) is not ctypes.c_int): n = ctypes.c_int(n)
    
    # Setting up "pts"
    if ((not issubclass(type(pts), numpy.ndarray)) or
        (not numpy.asarray(pts).flags.f_contiguous) or
        (not (pts.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'pts' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        pts = numpy.asarray(pts, dtype=ctypes.c_double, order='F')
    pts_dim_1 = ctypes.c_long(pts.shape[0])
    pts_dim_2 = ctypes.c_long(pts.shape[1])
    
    # Setting up "m"
    if (type(m) is not ctypes.c_int): m = ctypes.c_int(m)
    
    # Setting up "q"
    if ((not issubclass(type(q), numpy.ndarray)) or
        (not numpy.asarray(q).flags.f_contiguous) or
        (not (q.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'q' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        q = numpy.asarray(q, dtype=ctypes.c_double, order='F')
    q_dim_1 = ctypes.c_long(q.shape[0])
    q_dim_2 = ctypes.c_long(q.shape[1])
    
    # Setting up "simps"
    if ((not issubclass(type(simps), numpy.ndarray)) or
        (not numpy.asarray(simps).flags.f_contiguous) or
        (not (simps.dtype == numpy.dtype(ctypes.c_int)))):
        import warnings
        warnings.warn("The provided argument 'simps' was not an f_contiguous NumPy array of type 'ctypes.c_int' (or equivalent). Automatically converting (probably creating a full copy).")
        simps = numpy.asarray(simps, dtype=ctypes.c_int, order='F')
    simps_dim_1 = ctypes.c_long(simps.shape[0])
    simps_dim_2 = ctypes.c_long(simps.shape[1])
    
    # Setting up "weights"
    if ((not issubclass(type(weights), numpy.ndarray)) or
        (not numpy.asarray(weights).flags.f_contiguous) or
        (not (weights.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'weights' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        weights = numpy.asarray(weights, dtype=ctypes.c_double, order='F')
    weights_dim_1 = ctypes.c_long(weights.shape[0])
    weights_dim_2 = ctypes.c_long(weights.shape[1])
    
    # Setting up "ierr"
    if ((not issubclass(type(ierr), numpy.ndarray)) or
        (not numpy.asarray(ierr).flags.f_contiguous) or
        (not (ierr.dtype == numpy.dtype(ctypes.c_int)))):
        import warnings
        warnings.warn("The provided argument 'ierr' was not an f_contiguous NumPy array of type 'ctypes.c_int' (or equivalent). Automatically converting (probably creating a full copy).")
        ierr = numpy.asarray(ierr, dtype=ctypes.c_int, order='F')
    ierr_dim_1 = ctypes.c_long(ierr.shape[0])
    
    # Setting up "interp_in"
    interp_in_present = ctypes.c_bool(True)
    if (interp_in is None):
        interp_in_present = ctypes.c_bool(False)
        interp_in = numpy.zeros(shape=(1,1), dtype=ctypes.c_double, order='F')
    elif (type(interp_in) == bool) and (interp_in):
        interp_in = numpy.zeros(shape=(1,1), dtype=ctypes.c_double, order='F')
    elif ((not issubclass(type(interp_in), numpy.ndarray)) or
          (not numpy.asarray(interp_in).flags.f_contiguous) or
          (not (interp_in.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'interp_in' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        interp_in = numpy.asarray(interp_in, dtype=ctypes.c_double, order='F')
    if (interp_in_present):
        interp_in_dim_1 = ctypes.c_long(interp_in.shape[0])
        interp_in_dim_2 = ctypes.c_long(interp_in.shape[1])
    else:
        interp_in_dim_1 = ctypes.c_long()
        interp_in_dim_2 = ctypes.c_long()
    
    # Setting up "interp_out"
    interp_out_present = ctypes.c_bool(True)
    if (interp_out is None):
        interp_out_present = ctypes.c_bool(False)
        interp_out = numpy.zeros(shape=(1,1), dtype=ctypes.c_double, order='F')
    elif (type(interp_out) == bool) and (interp_out):
        interp_out = numpy.zeros(shape=(1,1), dtype=ctypes.c_double, order='F')
    elif ((not issubclass(type(interp_out), numpy.ndarray)) or
          (not numpy.asarray(interp_out).flags.f_contiguous) or
          (not (interp_out.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'interp_out' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        interp_out = numpy.asarray(interp_out, dtype=ctypes.c_double, order='F')
    if (interp_out_present):
        interp_out_dim_1 = ctypes.c_long(interp_out.shape[0])
        interp_out_dim_2 = ctypes.c_long(interp_out.shape[1])
    else:
        interp_out_dim_1 = ctypes.c_long()
        interp_out_dim_2 = ctypes.c_long()
    
    # Setting up "eps"
    eps_present = ctypes.c_bool(True)
    if (eps is None):
        eps_present = ctypes.c_bool(False)
        eps = ctypes.c_double()
    else:
        eps = ctypes.c_double(eps)
    if (type(eps) is not ctypes.c_double): eps = ctypes.c_double(eps)
    
    # Setting up "extrap"
    extrap_present = ctypes.c_bool(True)
    if (extrap is None):
        extrap_present = ctypes.c_bool(False)
        extrap = ctypes.c_double()
    else:
        extrap = ctypes.c_double(extrap)
    if (type(extrap) is not ctypes.c_double): extrap = ctypes.c_double(extrap)
    
    # Setting up "rnorm"
    rnorm_present = ctypes.c_bool(True)
    if (rnorm is None):
        rnorm_present = ctypes.c_bool(False)
        rnorm = numpy.zeros(shape=(1), dtype=ctypes.c_double, order='F')
    elif (type(rnorm) == bool) and (rnorm):
        rnorm = numpy.zeros(shape=(1), dtype=ctypes.c_double, order='F')
    elif ((not issubclass(type(rnorm), numpy.ndarray)) or
          (not numpy.asarray(rnorm).flags.f_contiguous) or
          (not (rnorm.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'rnorm' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        rnorm = numpy.asarray(rnorm, dtype=ctypes.c_double, order='F')
    if (rnorm_present):
        rnorm_dim_1 = ctypes.c_long(rnorm.shape[0])
    else:
        rnorm_dim_1 = ctypes.c_long()
    
    # Setting up "ibudget"
    ibudget_present = ctypes.c_bool(True)
    if (ibudget is None):
        ibudget_present = ctypes.c_bool(False)
        ibudget = ctypes.c_int()
    else:
        ibudget = ctypes.c_int(ibudget)
    if (type(ibudget) is not ctypes.c_int): ibudget = ctypes.c_int(ibudget)
    
    # Setting up "chain"
    chain_present = ctypes.c_bool(True)
    if (chain is None):
        chain_present = ctypes.c_bool(False)
        chain = ctypes.c_int()
    else:
        chain = ctypes.c_int(chain)
    if (type(chain) is not ctypes.c_int): chain = ctypes.c_int(chain)
    
    # Setting up "exact"
    exact_present = ctypes.c_bool(True)
    if (exact is None):
        exact_present = ctypes.c_bool(False)
        exact = ctypes.c_int()
    else:
        exact = ctypes.c_int(exact)
    if (type(exact) is not ctypes.c_int): exact = ctypes.c_int(exact)
    
    # Setting up "pmode"
    pmode_present = ctypes.c_bool(True)
    if (pmode is None):
        pmode_present = ctypes.c_bool(False)
        pmode = ctypes.c_int()
    else:
        pmode = ctypes.c_int(pmode)
    if (type(pmode) is not ctypes.c_int): pmode = ctypes.c_int(pmode)

    # Call C-accessible Fortran wrapper.
    clib.c_delaunaysparsep(ctypes.byref(d), ctypes.byref(n), ctypes.byref(pts_dim_1), ctypes.byref(pts_dim_2), ctypes.c_void_p(pts.ctypes.data), ctypes.byref(m), ctypes.byref(q_dim_1), ctypes.byref(q_dim_2), ctypes.c_void_p(q.ctypes.data), ctypes.byref(simps_dim_1), ctypes.byref(simps_dim_2), ctypes.c_void_p(simps.ctypes.data), ctypes.byref(weights_dim_1), ctypes.byref(weights_dim_2), ctypes.c_void_p(weights.ctypes.data), ctypes.byref(ierr_dim_1), ctypes.c_void_p(ierr.ctypes.data), ctypes.byref(interp_in_present), ctypes.byref(interp_in_dim_1), ctypes.byref(interp_in_dim_2), ctypes.c_void_p(interp_in.ctypes.data), ctypes.byref(interp_out_present), ctypes.byref(interp_out_dim_1), ctypes.byref(interp_out_dim_2), ctypes.c_void_p(interp_out.ctypes.data), ctypes.byref(eps_present), ctypes.byref(eps), ctypes.byref(extrap_present), ctypes.byref(extrap), ctypes.byref(rnorm_present), ctypes.byref(rnorm_dim_1), ctypes.c_void_p(rnorm.ctypes.data), ctypes.byref(ibudget_present), ctypes.byref(ibudget), ctypes.byref(chain_present), ctypes.byref(chain), ctypes.byref(exact_present), ctypes.byref(exact), ctypes.byref(pmode_present), ctypes.byref(pmode))

    # Return final results, 'INTENT(OUT)' arguments only.
    return pts, q, simps, weights, ierr, (interp_out if interp_out_present else None), (rnorm if rnorm_present else None)


# ----------------------------------------------
# Wrapper for the Fortran subroutine PROJECT

def project(d, n, pts, z, eps=None, weights=None):
    '''! Project a point outside the convex hull of the point set onto the convex
! hull by solving a non negatively constrained least squares problem with 1
! equality constraint (an instance of the WNNLS problem):
!
!    min_X   || MATMUL(PTS, X) - Z ||   s.t.   X >= 0, SUM(X) == 1
!
! The solution to the WNNLS problem stated above gives the projection Z_hat as
! a convex combination of the data points:
!
!   Z_hat = MATMUL(PTS, X).
!
! The above WNNLS problem is solved via R. Fletcher's QP solver: BQPD.
! Compared to other existing (D)WNNLS solvers, BQPD's flexible nature allows
! us to exploit the sparsity in the solution X, which should contain at most
! D positive entries (inactive constraints).
!
!
! On input:
!
! D is the dimension of the space for PTS and Z.
!
! N is the number of data points in PTS.
!
! PTS(1:D, 1:N) is a real valued matrix with N columns, each containing the
!    coordinates of a single data point in R^D.
!
! Z(1:D) is a real vector specifying the coordinates of a single
!    extrapolation point in R^D.
!
!
! On output:
!
! Z is overwritten with the result of the projection (labeled Z_hat above).
!
! RNORM contains the norm of the residual vector || Z - Z_hat ||.
!
! IERR contains an integer valued error flag (0=success) forwaded from BQPD.
!    Possible exit codes are listed below:
!
!       0 = solution obtained
!       1 = unbounded problem detected (f(x)<=fmin would occur)
!       2 = bl(i) > bu(i) for some i
!       3 = infeasible problem detected in Phase 1
!       4 = incorrect setting of m, n, kmax, mlp, mode or tol
!       5 = not enough space in lp
!       6 = not enough space for reduced Hessian matrix (increase kmax)
!       7 = not enough space for sparse factors (sparse code only)
!       8 = maximum number of unsuccessful restarts taken
!      -1 = a memory allocation error occurred (indicates the problem size is
!           too large for the additional memory overhead from extrapolation)
!
!
! Optional arguments:
!
! EPS contains the real working precision for the problem on input. By default,
!    EPS is assigned \sqrt{\mu} where \mu denotes the unit roundoff for the
!    machine. In general, any values that differ by less than EPS are judged
!    as equal, and any weights that are greater than -EPS are judged as
!    nonnegative.  EPS cannot take a value less than the default value of
!    \sqrt{\mu}. If any value less than \sqrt{\mu} is supplied, the default
!    value will be used instead automatically. Note that in order to ensure
!    that DELAUNAYSPARSE will be within tolerances, BQPD does not use the
!    value of EPS given here. Instead, BQPD is given a tolerance of
!    EPS ** 1.5.
!
! WEIGHTS(N) is assigned the projection weights on output, when present.
!
!
! Subroutines and functions directly referenced from BLAS are
!      DNRM2, DGEMV.
! BQPD, its utility functions, and its sparse linear algebra library are
! also referenced.
!
!
! This work is from a modification to the driver for BQPD by R. Fletcher,
! with modifications made by Tyler H. Chang (THC) and Sven Leyffer (SL).
!
!
! Version history:
!
! Version 1: Forked from BQPD by SL (Aug 2023)
!            Converted to f90 by THC (Sep 2023)'''
    
    # Setting up "d"
    if (type(d) is not ctypes.c_int): d = ctypes.c_int(d)
    
    # Setting up "n"
    if (type(n) is not ctypes.c_int): n = ctypes.c_int(n)
    
    # Setting up "pts"
    if ((not issubclass(type(pts), numpy.ndarray)) or
        (not numpy.asarray(pts).flags.f_contiguous) or
        (not (pts.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'pts' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        pts = numpy.asarray(pts, dtype=ctypes.c_double, order='F')
    pts_dim_1 = ctypes.c_long(pts.shape[0])
    pts_dim_2 = ctypes.c_long(pts.shape[1])
    
    # Setting up "z"
    if ((not issubclass(type(z), numpy.ndarray)) or
        (not numpy.asarray(z).flags.f_contiguous) or
        (not (z.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'z' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        z = numpy.asarray(z, dtype=ctypes.c_double, order='F')
    z_dim_1 = ctypes.c_long(z.shape[0])
    
    # Setting up "rnorm"
    rnorm = ctypes.c_double()
    
    # Setting up "ierr"
    ierr = ctypes.c_int()
    
    # Setting up "eps"
    eps_present = ctypes.c_bool(True)
    if (eps is None):
        eps_present = ctypes.c_bool(False)
        eps = ctypes.c_double()
    else:
        eps = ctypes.c_double(eps)
    if (type(eps) is not ctypes.c_double): eps = ctypes.c_double(eps)
    
    # Setting up "weights"
    weights_present = ctypes.c_bool(True)
    if (weights is None):
        weights_present = ctypes.c_bool(False)
        weights = numpy.zeros(shape=(1), dtype=ctypes.c_double, order='F')
    elif (type(weights) == bool) and (weights):
        weights = numpy.zeros(shape=(1), dtype=ctypes.c_double, order='F')
    elif ((not issubclass(type(weights), numpy.ndarray)) or
          (not numpy.asarray(weights).flags.f_contiguous) or
          (not (weights.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'weights' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        weights = numpy.asarray(weights, dtype=ctypes.c_double, order='F')
    if (weights_present):
        weights_dim_1 = ctypes.c_long(weights.shape[0])
    else:
        weights_dim_1 = ctypes.c_long()

    # Call C-accessible Fortran wrapper.
    clib.c_project(ctypes.byref(d), ctypes.byref(n), ctypes.byref(pts_dim_1), ctypes.byref(pts_dim_2), ctypes.c_void_p(pts.ctypes.data), ctypes.byref(z_dim_1), ctypes.c_void_p(z.ctypes.data), ctypes.byref(rnorm), ctypes.byref(ierr), ctypes.byref(eps_present), ctypes.byref(eps), ctypes.byref(weights_present), ctypes.byref(weights_dim_1), ctypes.c_void_p(weights.ctypes.data))

    # Return final results, 'INTENT(OUT)' arguments only.
    return z, rnorm.value, ierr.value, (weights if weights_present else None)


# ----------------------------------------------
# Wrapper for the Fortran subroutine GDOTX

def gdotx(n, x, ws, lws, v):
    '''! Auxiliary subroutine needed by BQPD to compute the gradient v = <G, x>
!    where G = <PTS, PTS>. Since PTS is a DxN matrix and N could be large,
!    we will not explicitly form the NxN matrix G.
!
!
! On input:
!
! N is the integer length of the gradient vector.
!
! X is an array of length N containing the current iterate.
!
! WS is a real valued workspace array, whose first D*N entries contain PTS,
!    stored in (flattened) column-major order.
!
! LWS is an integer valued workspace array, with LWS(0) contains the value
!    of D.
!
!
! On output:
!
! V is a real valued vector of length D containing < <PTS, PTS>, x >.
!
!
! Uses the external BLAS subroutine DGEMV.'''
    
    # Setting up "n"
    if (type(n) is not ctypes.c_int): n = ctypes.c_int(n)
    
    # Setting up "x"
    if ((not issubclass(type(x), numpy.ndarray)) or
        (not numpy.asarray(x).flags.f_contiguous) or
        (not (x.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'x' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        x = numpy.asarray(x, dtype=ctypes.c_double, order='F')
    x_dim_1 = ctypes.c_long(x.shape[0])
    
    # Setting up "ws"
    if ((not issubclass(type(ws), numpy.ndarray)) or
        (not numpy.asarray(ws).flags.f_contiguous) or
        (not (ws.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'ws' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        ws = numpy.asarray(ws, dtype=ctypes.c_double, order='F')
    ws_dim_1 = ctypes.c_long(ws.shape[0])
    
    # Setting up "lws"
    if ((not issubclass(type(lws), numpy.ndarray)) or
        (not numpy.asarray(lws).flags.f_contiguous) or
        (not (lws.dtype == numpy.dtype(ctypes.c_int)))):
        import warnings
        warnings.warn("The provided argument 'lws' was not an f_contiguous NumPy array of type 'ctypes.c_int' (or equivalent). Automatically converting (probably creating a full copy).")
        lws = numpy.asarray(lws, dtype=ctypes.c_int, order='F')
    lws_dim_1 = ctypes.c_long(lws.shape[0])
    
    # Setting up "v"
    if ((not issubclass(type(v), numpy.ndarray)) or
        (not numpy.asarray(v).flags.f_contiguous) or
        (not (v.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'v' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        v = numpy.asarray(v, dtype=ctypes.c_double, order='F')
    v_dim_1 = ctypes.c_long(v.shape[0])

    # Call C-accessible Fortran wrapper.
    clib.c_gdotx(ctypes.byref(n), ctypes.byref(x_dim_1), ctypes.c_void_p(x.ctypes.data), ctypes.byref(ws_dim_1), ctypes.c_void_p(ws.ctypes.data), ctypes.byref(lws_dim_1), ctypes.c_void_p(lws.ctypes.data), ctypes.byref(v_dim_1), ctypes.c_void_p(v.ctypes.data))

    # Return final results, 'INTENT(OUT)' arguments only.
    return ws, v


class real_precision:
    ''''''

    # Declare 'r8'
    def get_r8(self):
        r8 = ctypes.c_int()
        clib.real_precision_get_r8(ctypes.byref(r8))
        return r8.value
    def set_r8(self, r8):
        raise(NotImplementedError('Module attributes with PARAMETER status cannot be set.'))
    r8 = property(get_r8, set_r8)

real_precision = real_precision()


class delsparse_mod:
    '''! This module contains the REAL_PRECISION R8 data type for 64-bit arithmetic
! and interface blocks for the DELAUNAYSPARSES and DELAUNAYSPARSEP
! subroutines for computing the Delaunay simplices containing interpolation
! points Q in R^D given data points PTS.'''

delsparse_mod = delsparse_mod()

