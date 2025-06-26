cimport cython
cimport numpy as np 
import numpy

ctypedef np.float64_t DTYPE_F64_t

from libc.math cimport fabs

cpdef DTYPE_F64_t rot_fit_BB(DTYPE_F64_t depro_radius, list params):
    '''
    BB function to fit rotation curve data to

    PARAMETERS
    ==========

    depro_radius : float
        Deprojected radius as taken from the [PLATE]-[FIBERID] rotation curve 
        data file (in units of kpc); the "x" data of the rotation curve equation
    params : list
        model parameter values

    RETURNS
    =======
    v : float
        The rotation curve equation with the given parameters at the given depro_radius


    '''

    cdef DTYPE_F64_t v_max
    cdef DTYPE_F64_t r_turn
    cdef DTYPE_F64_t alpha
    cdef DTYPE_F64_t v

    v_max, r_turn, alpha = params

    v = v_max * fabs(depro_radius) / (r_turn**alpha + fabs(depro_radius)**alpha)**(1/alpha)

    if depro_radius < 0:
        v = v * -1


    return v


cpdef rot_fit_BB_loop(DTYPE_F64_t [:] depro_radius, list params):
    '''
    BB function to fit rotation curve data to

    PARAMETERS
    ==========

    depro_radius : float
        Deprojected radius as taken from the [PLATE]-[FIBERID] rotation curve 
        data file (in units of kpc); the "x" data of the rotation curve equation
    params : list
        model parameter values

    RETURNS
    =======
    v : float
        The rotation curve equation with the given parameters at the given depro_radius


    '''

    cdef DTYPE_F64_t v_max
    cdef DTYPE_F64_t r_turn
    cdef DTYPE_F64_t alpha
    cdef np.ndarray[DTYPE_F64_t, ndim=1] v = numpy.zeros(depro_radius.shape[0])
    cdef int N = depro_radius.shape[0]
    cdef int i
    

    v_max, r_turn, alpha = params

    for i in range(N):

        v[i] = v_max * fabs(depro_radius[i]) / (r_turn**alpha + fabs(depro_radius[i])**alpha)**(1/alpha)

        if depro_radius[i] < 0:
            v[i] = v[i] * -1


    return v