#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import numpy.ma as ma
import numpy as np

def calc_v(r_array,params):
    '''
    Function
    ^^^^^^^^
    Calculate the velocity at each radius
    =====================================

    Parameters
    ^^^^^^^^^^

    params : length-7 list or array
    model parameters in the following order:
    
        vmax : Quantity object
            maximum velocity (km/s)

        alpha : float
            sharpness of the galaxy rotation curve

        Rturn : float
            radius where the rotation curve changes from rising to flat

        PA : float
            position angle of the galaxy (radians)

        i_angle : float
            inclination angle of the galaxy (radians)

        clean_coords_i : int
            i index of the galaxy center in pixels

        clean_coords_j : int
            j index of the galaxy center in pixels

    r_array : array

    radius array with deprojected length

    ====================================

    Return
    ^^^^^^

    v : array
    the rotational velocity at each radius
    '''
    vmax = params[0]
    Rturn = params[1]
    alpha = params[2]
    
    v = ((vmax*r_array)/((((Rturn)**alpha)+((r_array)**alpha))**(1/alpha)))
    
    return v
