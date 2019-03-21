cimport cython
import numpy as np
from libc.math cimport log
from cython.parallel cimport prange, parallel

@cython.boundscheck(False)
cdef double entropy_ufun(double p, double q, double tol = 1e-50) nogil:
    if p < tol:
        return 0
    elif q < tol:
        return 0
    else:
        return p * log(p / q)
    
@cython.boundscheck(False)
cpdef double entropy_cy(double[:] p, double[:] q) nogil:
    cdef double epy = 0
    cdef int idx
    
    for idx in range(len(p)):
        epy += entropy_ufun(p[idx], q[idx])
    return epy

@cython.boundscheck(False)
cpdef double jsdiv_cy(double[:] p, double[:] q) nogil:
    return 0.5 * (entropy_cy(p, q) + entropy_cy(q, p))

@cython.boundscheck(False)
def pairwise_jsdiv_cy(double[:,:] mat_out, double[:,:] mat_p1, double[:,:] mat_p2):
    """
    Arg:
        mat_p1, mat_p2: similarity matrix
    """
    #mat_out = np.zeros((mat_p1.shape[1], mat_p2.shape[1]))
    cdef int i, j, m, n
    cdef double[:] p1, p2
    m = mat_p1.shape[1]
    n = mat_p2.shape[1]
    
    for i in range(m):
        for j in range(n):
            p1 = mat_p1[:, i]
            p2 = mat_p2[:, j]
            mat_out[i, j] = jsdiv_cy(p1, p2)