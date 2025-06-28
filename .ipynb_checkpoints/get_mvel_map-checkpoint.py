#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import numpy.ma as ma

from calc_v import calc_v
from deproject import deproject_spaxel


def get_mvel_map(map_shape, vmax, alpha, Rturn, checkedPA, i_angle , clean_coords_i, clean_coords_j ,mask_f, r_convert):

    r = np.zeros(map_shape)
    theta = np.zeros(map_shape)
    clean_coords = [clean_coords_i,clean_coords_j]

    for i in range(map_shape[0]):
        for j in range(map_shape[1]):

            # De-projected radius for the current point
            r[i,j], theta[i,j] = deproject_spaxel((i,j), clean_coords, checkedPA, i_angle)
                
    r_array = r*r_convert
    
    params = [vmax, Rturn, alpha]
    
    v = calc_v(r_array, params)


    vel_map = v*np.sin(i_angle)*np.cos(theta)
    mvel_map = ma.array(vel_map, mask = mask_f)

    return mvel_map


def calc_chi2(params, vel_map, vel_ivar, map_shape, r_convert, mask_f):
    '''
    vmax

    alpha

    rturn

    PA

    i_angle

    c_coord_i

    c_coord_j

    '''
    model = get_mvel_map(map_shape, params[0],params[1],params[2],params[3], params[4], params[5],params[6], mask_f, r_convert )      
    
    chi2 = np.sum(((model-vel_map)**2)*vel_ivar)

    return chi2