#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import numpy.ma as ma

import numpy as np

def calc_v(r_array,params  ):
    
    vmax = params[0]
    Rturn = params[1]
    alpha = params[2]
    
    v = ((vmax*r_array)/((((Rturn)**alpha)+((r_array)**alpha))**(1/alpha)))
    
    return v
