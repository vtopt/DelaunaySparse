
# Import the Delaunay Fortran code.
import delsparse

# List out the "help" documentation.
# help(delsparse)

# Return the source point indices and weights associated with a set of
# interpolation points. Takes points in row-major (C style) format.
# 
# INPUTS:
#   pts -- 2D Numpy array of float64 points, where each row is one point.
#   q -- 2D numpy array of float64 points where Delaunay predictions
#        are to be made, where each row is one point.
# 
# OUTPUT:
#   (indices, weights) -- Where "indices" is a 2D NumPy array of integers
#      and each row, i, enumerates the indices of rows in "pts" that are
#      the vertices of the simplex containing q[i], and each corresponding
#      row of weights (a 2D NumPy array of float64) provides the convex
#      weights such that q[i] = np.dot(pts[indices[i]], weights[i]).
# 
def delaunay_simplex(pts, q, allow_extrapolation=True, print_errors=True,
                     parallel=True, pmode=None, chain=None,
                     ibudget=10000, epsilon=2**(-23), check_spacing=False):
    # Enable parallelism.
    if parallel:
        import os
        os.environ["OMP_NESTED"] = "TRUE"
    # Import NumPy.
    import numpy as np
    # Get the predictions from VTdelaunay
    pts_in = np.asarray(pts.T, dtype=np.float64, order="F")
    p_in = np.asarray(q.T, dtype=np.float64, order="F")
    simp_out = np.ones(shape=(p_in.shape[0]+1, p_in.shape[1]), 
                       dtype=np.int32, order="F")
    weights_out = np.ones(shape=(p_in.shape[0]+1, p_in.shape[1]), 
                          dtype=np.float64, order="F")
    error_out = np.ones(shape=(p_in.shape[1],), 
                        dtype=np.int32, order="F")
    if parallel:
        delsparse.delaunaysparsep(pts_in.shape[0], pts_in.shape[1],
                                  pts_in, p_in.shape[1], p_in, simp_out,
                                  weights_out, error_out, extrap=100.0,
                                  pmode=pmode, ibudget=ibudget,
                                  eps=epsilon, chain=chain,
                                  exact=check_spacing)
    else:
        delsparse.delaunaysparses(pts_in.shape[0], pts_in.shape[1],
                                  pts_in, p_in.shape[1], p_in, simp_out,
                                  weights_out, error_out, extrap=100.0, 
                                  ibudget=ibudget, eps=epsilon,
                                  chain=chain, exact=check_spacing)
    # Remove "extrapolation" errors if the user doesn't care.
    if allow_extrapolation: error_out = np.where(error_out == 1, 0, error_out)
    else:
        if 1 in error_out:
            class Extrapolation(Exception): pass
            raise(Extrapolation("Encountered extrapolation point when making Delaunay prediction."))
    # Handle any errors that may have occurred.
    if (sum(error_out) != 0):
        if print_errors:
            unique_errors = sorted(np.unique(error_out))
            print(" [Delaunay errors:",end="")
            for e in unique_errors:
                if (e == 0): continue
                indices = tuple(str(i) for i in range(len(error_out))
                                if (error_out[i] == e))
                if (len(indices) > 5): indices = indices[:2] + ('...',) + indices[-2:]
                print(" %3i"%e,"at","{"+",".join(indices)+"}", end=";")
            print("] ")
        # Reset the errors to simplex of 1s (to be 0) and weights of 0s.
        bad_indices = (error_out > (1 if allow_extrapolation else 0))
        simp_out[:,bad_indices] = 1
        weights_out[:,bad_indices] = 0
    # Adjust the output simplices and weights to be expected shape.
    indices  = simp_out.T - 1
    weights = weights_out.T
    # Return the appropriate shaped pair of points and weights
    return (indices, weights)


# Declare some test function.
import numpy as np
f = lambda x: 3*x[0]+.5*np.cos(8*x[0])+np.sin(5*x[-1])
np.random.seed(0)

# Construct a function that converts indices and weights into a real number prediction.
def delaunay_approx(q, points, values):
    q = np.array(q, dtype=float)
    if len(q.shape) == 1:
        pts, wts = delaunay_simplex(points.copy(), np.reshape(q,(1,len(q))))
        return sum(values[pts[0]]*wts[0])
    else:
        pts, wts = delaunay_simplex(points.copy(), q)
        return np.asarray([sum(values[i]*w) for (i,w) in zip(pts,wts)], dtype=float)

# Generate test data.
test_size = 1000
q = np.random.random(size=(test_size,2))
f_q = np.asarray(list(map(f,q)), dtype=float)

# Show the convergence on the true function.
for train_size in (10, 50, 100, 200, 500, 1000, 5000, 10000):
    # Construct "training" data.
    x = np.random.random(size=(train_size,2))
    y = np.asarray([f(_) for _ in x], dtype=float)
    # Approximate at points.
    f_hat = delaunay_approx(q, x, y)
    # Compute errors.
    abs_error = abs(f_hat - f_q)
    avg_abs_error = sum(abs_error) / test_size
    max_abs_error = max(abs_error)
    # Show errors.
    print()
    print("Train size:", train_size)
    print("  maximum absolute error:  %.2f"%(max_abs_error))
    print("  average absolute error:  %.2f"%(avg_abs_error))

