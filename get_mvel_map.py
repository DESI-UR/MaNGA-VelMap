#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import numpy.ma as ma

from calc_v import calc_v
from deproject import deproject_spaxel


def get_mvel_map(map_shape, vmax, alpha, Rturn, checkedPA, i_angle , clean_coords_i, clean_coords_j ,mask_f, r_convert):
    '''
    Function
    ^^^^^^^^
    Create model velocy map for chosen galaxy
    =========================================
    Parameters
    ^^^^^^^^^^
    map_shape : length-2 tuple
    shape of the map (#,#)

    vmax : Quantity object
    maximum velocity (km/s)

    alpha : float
    sharpness of the galaxy rotation curve

    Rturn : float
    the radius where the rotation curve changes from increasing to flat

    checkedPA : float
    the corrected position angle (radians)

    i_angle : float
    the inclination angle of the galaxy (radians)

    clean_coords_i : int
    the i index of flux based center coordinates

    clean_coords_j : int
    the j index of flux based center coordinates

    mask_f : array
    mask to be used on model velocity map

    r_convert : float
    covnersion factor for radius to radius in kpc

    ====================================================================

    Return
    ^^^^^^
    mvel_map : masked array
    velocity map model 
    
    '''
    #make arrays for r and theta in the same shape as map
    r = np.zeros(map_shape)
    theta = np.zeros(map_shape)

    #combine i and j coords for use in deproject_spaxel
    clean_coords = [clean_coords_i,clean_coords_j]

    for i in range(map_shape[0]):
        for j in range(map_shape[1]):

            # De-projected radius for the current point
            r[i,j], theta[i,j] = deproject_spaxel((i,j), clean_coords, checkedPA, i_angle)

    #convert the r out array 
    r_array = r*r_convert

    #list of free parameters for calc_v
    params = [vmax, Rturn, alpha]

    #call calc_v to get velocities at each point
    v = calc_v(r_array, params)

    #adjust the map for our line of sight velocity
    vel_map = v*np.sin(i_angle)*np.cos(theta)
    #apply mask to model velocity map
    mvel_map = ma.array(vel_map, mask = mask_f)

    return mvel_map


def calc_chi2(params, vel_map, vel_ivar, map_shape, r_convert, mask_f):
    '''
    Function
    ^^^^^^^^
    Calculate chi-squared value comparing model and observed velocity maps
    ======================================================================

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

    vel_map : masked array
        observed velocity map

    vel_ivar : array
        inverse variance of the observed velocity values 

    map_shape : length-2 tuple
        shape of the velocity map (#, #)

    r_convert : float
        conversion factor from pixels to physical radius (e.g., kpc)

    mask_f : array
        mask to be applied to the model velocity map

    ======================================================================

    Return
    ^^^^^^
    chi2 : float
        chi-squared value quantifying the fit between model and observed velocities

    '''
    model = get_mvel_map(map_shape, params[0],params[1],params[2],params[3], params[4], params[5],params[6], mask_f, r_convert )      
    
    chi2 = np.sum(((model-vel_map)**2)*vel_ivar)

    return chi2
