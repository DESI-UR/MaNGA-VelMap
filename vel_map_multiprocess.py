################################################################################
# IMPORT MODULES
#-------------------------------------------------------------------------------
import datetime
START = datetime.datetime.now()

import os.path, warnings

import numpy as np

from astropy.table import Table
import astropy.units as u

from multiprocessing import Process, Queue, Value
from queue import Empty

from ctypes import c_long

#importing all packages


import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from astropy.table import Table
from astropy.io import fits
import numpy.ma as ma
# from PIL import Image
import astropy.units as u
# from astropy.coordinates import SkyCoord
# from astropy.wcs import WCS
import astropy.constants as const
# from GenerateCutout import get_cutout,get_cutout_fits
from get_mvel_map import get_mvel_map, calc_chi2
from tqdm import tqdm
from mapSmoothness_functions import how_smooth 
import traceback
import os
import matplotlib.image as mpimg
from deproject import deproject_spaxel
H_0 = 100 *(u.km/u.s)/u.Mpc
q0 = 0.2
c = const.c.to('km/s')


import sys
sys.path.insert(1, '/home/ebenitez/Documents/' ) # CHANGE PATH AS NECESSARY
from mapSmoothness_functions import how_smooth

#warnings.simplefilter('ignore', np.RankWarning)
#warnings.simplefilter('ignore', RuntimeWarning)
################################################################################


def process_1_galaxy(job_queue, i, 
                     return_queue,  
                     drpall,
                     drpall_dict 
                     ):
    '''
    Main body of for-loop for processing one galaxy.
    '''

    ############################################################################
    # Open file to which we can write all the print statements.
    #---------------------------------------------------------------------------
    outfile = open('Process_' + str(i) + '_output.txt', 'wt')
    sys.stdout = outfile
    sys.stderr = outfile
    ############################################################################

    while True:
        try: 
            i_DRP = job_queue.get(timeout=1.0)
        except Empty:
        
            print('Queue is empty!', flush=True)
        
            outfile.close()
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
            
            print('Worker', i, 'redirected stdout and stderr.', flush=True)
            
            return_queue.close()
            
            print('Worker', i, 'closed the return queue.', flush=True)
            
            return_queue.join_thread()
            
            print('Worker', i, 'joined the return queue.', flush=True)
            
            job_queue.close()
            
            print('Worker', i, 'closed the job queue.', flush=True)
            
            job_queue.join_thread()
            
            print('Worker', i, 'returned successfully', datetime.datetime.now(), flush=True)
            return
    
        start = datetime.datetime.now()

        ########################################################################
        # Fitting process here
        #-----------------------------------------------------------------------

        plateifu = drpall['plateifu'][i_DRP] 
        loc = drpall_dict[plateifu]
        
        #import data from plateifu

        plate=plateifu.split('-',1)
        plate1=plate[0]
        plate2=plate[1]

        ########################################################################
        # Check if gal already ran
        #-----------------------------------------------------------------------

        refit = ['10516-1902', '9044-9101', '9187-12702', '8239-1901',
                 '8595-3704', '8444-9101', '8727-12702', '8453-6101',
                 '11834-12701','8259-12704',  '8999-1901', '8549-1901'
                 ]

        if drpall['smoothness_score'][loc] >=0 and plateifu not in refit:
            output_tuple = (None, None, None, None, 
                            None, None, None, None, None, None,
                            None, None, None, 
                            None, None, None, None, None, None, loc)
            return_queue.put(output_tuple)
            print(plate, ' already processed', flush=True)
            continue


        ########################################################################
        # Check galaxy target
        #-----------------------------------------------------------------------
 
        ttype = drpall['TType'][loc]
        mng = drpall['mngtarg1'][loc]

        if ttype <= 0 or mng <=0:
            output_tuple = (None, None, None, None, 
                            None, None, None, None, None, None,
                            None, None, None, 
                            None, None, None, None, None, None, loc)
            return_queue.put(output_tuple)
            print(plate, ' wrong t-type or target type', flush=True)
            continue
            
        print('PLATEIFU : ',plateifu)
        cube_fn = '/scratch/kdougla7/data/SDSS/dr17/manga/spectro/analysis/v3_1_1/3.1.0/HYB10-MILESHC-MASTARSSP/'+plate1+'/'+plate2+'/manga-'+plateifu+'-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz'
        
        ########################################################################
        # Extract maps
        #-----------------------------------------------------------------------

        try:
            cube = fits.open(cube_fn)                                           #Opening cube file to import all neccesary data
            stellar_vel = cube['STELLAR_VEL'].data                              
            stellar_mask = cube['STELLAR_VEL_MASK'].data
            stellar_vel_ivar = cube['STELLAR_VEL_IVAR'].data
            halpha_vel = cube['EMLINE_GVEL'].data[23]        
            halpha_gvel_mask = cube['EMLINE_GVEL_MASK'].data[23]
            halpha_gvel_ivar = cube['EMLINE_GVEL_IVAR'].data[23]
            # ellip_radius = cube['SPX_ELLCOO'].data[3]
            spx_x = cube['SPX_SKYCOO'].data[0]
            spx_y = cube['SPX_SKYCOO'].data[1]
            flux = cube['SPX_MFLUX'].data
            halpha_flux = cube['EMLINE_GFLUX'].data[23]
            halpha_flux_ivar = cube['EMLINE_GFLUX_IVAR'].data[23]
            binarea = cube['BIN_AREA'].data
            cube.close() 

        except:
            print('Could not extract data for ', plate, flush=True)
            output_tuple = (None, None, None, None, 
                            None, None, None, None, None, None,
                            None, None, None, 
                            None, None, None, None, None, None, loc)
            return_queue.put(output_tuple)
            continue     
        
        
        print(plate, ' cube data extracted', flush=True)

        ########################################################################
        # Smoothness score
        #-----------------------------------------------------------------------


        map_smoothness = how_smooth(halpha_vel, halpha_gvel_mask)                                                                            #cut based on smoothness to ensure usability of galaxy by
        map_smoothness_5sigma = how_smooth(halpha_vel, np.logical_or(halpha_gvel_mask>0, np.abs(halpha_flux*np.sqrt(halpha_flux_ivar)<5)))   #avoiding random behavior
        
        if map_smoothness >= 2 and map_smoothness_5sigma >= 2:

            print(plate, ' not smooth enough to fit', flush=True)
            output_tuple = (map_smoothness, None, None, None, 
                            None, None, None, None, None, None,
                            None, None, None, 
                            None, None, None, None, None, None, loc)
            return_queue.put(output_tuple)
            continue
            
        print("smoothness score: ", map_smoothness, flush=True)
        
        ########################################################################
        # Maps and photometry 
        #-----------------------------------------------------------------------

        halpha_mask = np.logical_or(halpha_gvel_mask, np.abs(halpha_flux*np.sqrt(halpha_flux_ivar)) < 5)         #masking halpha velocity map with default mask and 
        mhalpha_vel = ma.array(halpha_vel, mask=halpha_mask)                                                     #the addition of masks on low strength points
        
        #percent spaxel cut
        if (mhalpha_vel.count()/mhalpha_vel.size)*100<5:
            print(plate, ' not enough data to fit', flush=True)
            output_tuple = (map_smoothness, None, None, None, 
                            None, None, None, None, None, None,
                            None, None, None, 
                            None, None, None, None, None, None, loc)
            return_queue.put(output_tuple)
            continue

        mhalpha_ivar = ma.array (halpha_gvel_ivar, mask = halpha_mask )  #masking the halpha inverse variance data
        PA = drpall[loc]['nsa_elpetro_phi']                                   #PA is the position angle of chosen galaxy
        ba = drpall[loc]['nsa_elpetro_ba']                                    #ba is the axes ratio of chosen galaxy
        z = drpall[loc]['z']                                                  #z is red-shift value given by sdss 
        mask_f = halpha_mask
        # x = -spx_x                                                            #x and y are spaxel distance from the center of the galaxy
        # y = spx_y
            
        
        gal_distance = (c*z)/H_0
        gal_dist_kpc = gal_distance.to(u.kpc).value                           #gal_dist_kpc is the distance to the observed galaxy in kpc units
        i_angle = np.arccos(np.sqrt(((ba)**2)-(q0**2)/(1-q0**2)))
        spax_size = 0.5*(1/60)*(1/60)*(np.pi/180)                             #spax_size is the conversion from spaxels to radians
        r_convert = gal_dist_kpc*spax_size                                    #r_convert is the distance of a point from the center of a galaxy in kpc
        map_shape = mhalpha_vel.shape

        mx = np.max(flux)                                         #locating the peak value on the flux map 
        indices = np.where(flux == mx)                            #and assigning that point to the center of the map  
        clean_coords = [int(indices[0][0]), int(indices[1][0])]
        
        # sys_vel = mhalpha_vel[tuple(clean_coords)]*u.km/u.s 

        # theta = np.radians(PA-90)                      #theta is PA reoriented to the positive x axis for calculations
        # for x in range(15,clean_coords[0]):       #goes 15 left to 15 right
        #     if x < 0 or x >= mhalpha_vel.shape[0]:
        #         continue
            
        #     y = clean_coords[1]- round((clean_coords[0]-x) * np.tan(theta))               #y's are found based on the slope given by PA
        #     if y < 0 or y >= mhalpha_vel.shape[1]:
        #         continue
                
        #     if ma.is_masked(mhalpha_vel[x, y]):
        #         continue
        #     break
            
            
        # if (mhalpha_vel[x,y]<0):                        # if velocity comes back negative the position angle will be flipped 180 deg, otherwise left alone
        #     checkedPA = (PA + 180) *(np.pi/180)
        # else:
        #     checkedPA = PA*(np.pi/180)


        theta = np.radians(PA-90)   #theta is PA reoriented to the positive x axis for calculations
        checkedPA = np.radians(PA)
        f = 20

        try: 
            for x in range(int(clean_coords[0]/2),clean_coords[0]):   
                if x < 0 or x >= mhalpha_vel.shape[0]:
                    continue
                
                y = clean_coords[1] - round((clean_coords[0] - x) * np.tan(theta)) #y's are found based on the slope given by PA

                if PA < 15 and PA > 0: # if position angle is near 0
                    x = clean_coords[0]
                    y = clean_coords[1] - f
                    while y < 0 or ma.is_masked(mhalpha_vel[x,y]):
                        f -= 1
                        y = clean_coords[1] - f


                if y < 0 or y >= mhalpha_vel.shape[1]:
                    continue
                    
                if ma.is_masked(mhalpha_vel[x, y]):
                    continue
                
                else:
                    if (mhalpha_vel[x,y]<0):
                        checkedPA = (PA + 180) *(np.pi/180)
                    else:
                        checkedPA = PA*(np.pi/180)
                    break

        except:
            print('couldnt check PA, using NSA value', flush=True)


        ########################################################################
        # fitting guesses and bounds
        #-----------------------------------------------------------------------

        vmax_guess = abs(mhalpha_vel).max()/np.sin(i_angle)                             #set up initial guesses to start chi2 minimization
        alpha_guess = 1
        
        halphamx = np.max(mhalpha_vel)                                                  #defining Rturn_guess with halfway from the center spaxel
        mvel_indices = np.where(mhalpha_vel == halphamx)                                #to the fastest spaxel
        vmax_coords = [int(mvel_indices[0][0]), int(mvel_indices[1][0])]
        
        dep_mvel_dist , _= deproject_spaxel(vmax_coords, clean_coords, checkedPA, i_angle)

        
        dmd_kpc = dep_mvel_dist*r_convert
        Rturn_guess = dmd_kpc/2      

        PA_guess = checkedPA
        i_angle_guess = i_angle
        center_guess = clean_coords
        #sys_vel_guess = mhalpha_vel[tuple(clean_coords)]
        sys_vel_guess =  0
        x0 = [vmax_guess, alpha_guess, Rturn_guess, PA_guess, i_angle_guess, center_guess[0],center_guess[1],sys_vel_guess] #make list for guesses to be input as a list

        args = (mhalpha_vel, mhalpha_ivar, map_shape, r_convert, mask_f)    #non-variable values to go into minimzation function

        vmax_bound = (0, 1000)                                     #setting reasonable bounds for variable used in minimization
        alpha_bound = (0.001,100)
        Rturn_bound = (0.001,100)
        PA_bound = (0,2*np.pi)
        sys_vel_bound = (-100,100)
        
        i_anglelow = np.max([0,i_angle-np.radians(15)])
        i_anglehigh = np.min([np.radians(90), i_angle+np.radians(15)])
        i_angle_bound = (i_anglelow, i_anglehigh)

        centeri_bound = (center_guess[0]-5,center_guess[0]+5)
        centerj_bound = (center_guess[1]-5,center_guess[1]+5)

        bounds = np.array([vmax_bound,alpha_bound,Rturn_bound,PA_bound,i_angle_bound,centeri_bound,centerj_bound, sys_vel_bound ]) 

        print(plate, ' fitting Halpha map', flush=True)
        start = datetime.datetime.now()

        ########################################################################
        # Fit H-alpha map
        #-----------------------------------------------------------------------

        try:

            result = minimize(     #using minimization function on chi2 function
                calc_chi2,
                x0,
                args = args,
                method = 'Powell',
                bounds = bounds
                ) 

        except:
            print(plate, 'CRASHED! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<', 
                            flush=True)
            output_tuple = (map_smoothness, None, None, None, 
                            None, None, None, None, None, None,
                            None, None, None, 
                            None, None, None, None, None, None, loc)
            return_queue.put(output_tuple)
            continue
                        
        
        fit_time = datetime.datetime.now() - start

        print(plate, " Halpha velocity map fit", fit_time, flush=True)

        if result.success == False:
            print(plate, ' fit failed', flush=True)
            output_tuple = (map_smoothness, None, None, None, 
                            None, None, None, None, None, None,
                            None, None, None, 
                            None, None, None, None, None, None, loc)
            return_queue.put(output_tuple)
            continue

        ########################################################################
        # Extract best fit parameters
        #-----------------------------------------------------------------------

        #model based on chi2
        nmvel_map, r, theta = get_mvel_map(map_shape,                       #get_mvel_map is used to get the model vel map
                        result.x[0], 
                        result.x[1],
                        result.x[2],
                        result.x[3],
                        result.x[4],
                        result.x[5],
                        result.x[6],
                        result.x[7],
                        mask_f,r_convert
                       )

        #defining result outputs
        vmax_ideal = result.x[0] 
        alpha_ideal = result.x[1]
        Rturn_ideal = result.x[2]
        PA_ideal = result.x[3]
        i_angle_ideal = result.x[4]
        x_cent_ideal = result.x[5]
        y_cent_ideal = result.x[6]
        sys_vel_ideal = result.x[7]
        chi2_nu = result.fun / (mhalpha_vel.count() - 8)
        print('H alpha map chi2_nu', chi2_nu, flush=True)

        ########################################################################
        # Make H-alpha plots
        #-----------------------------------------------------------------------
        # halpha plot save
        mhalpha_vel = ma.array(halpha_vel, mask=mask_f) #re-applied mask

        val_max_h = mhalpha_vel.max()        #simple plotting section for visualization
        val_min_h = mhalpha_vel.min()
        

        if (val_max_h >= abs(val_min_h)):
            lim = val_max_h
        else:
            lim = abs(val_min_h)

        plt.imshow(mhalpha_vel,cmap = 'bwr',vmin = -lim, vmax = lim)     #creates a color map of data, (data set, cmap = 'desired color code')
        plt.colorbar(label = 'Velocity [km/s]')  #following are just labels
        plt.gca().invert_yaxis()
        plt.xlabel('spaxel')
        plt.ylabel('spaxel')
        plt.title(plateifu + r' H-$\alpha$ Velocity (data)')
        plt.savefig(data_folder + 'Plots/Halpha/'+plateifu+'_Halpha_vel_plot.png')
        plt.close()
        
      
        #model plot save
        val_max_m = nmvel_map.max()        #simple plotting section for visualization
        val_min_m = nmvel_map.min()
        

        if (val_max_m >= abs(val_min_m)):
            lim = val_max_m
        else:
            lim = abs(val_min_m)

        plt.imshow(nmvel_map,cmap = 'bwr',vmin = -lim, vmax = lim)     #creates a color map of data, (data set, cmap = 'desired color code')
        plt.colorbar(label = 'Velocity [km/s]')  #following are just labels
        plt.gca().invert_yaxis()
        plt.xlabel('spaxel')
        plt.ylabel('spaxel')
        plt.title(plateifu+ r' H-$\alpha$ Velocity (model)')
        plt.savefig(data_folder + 'Plots/model/'+plateifu+'_Halpha_vel_model.png')
        plt.close()
        
        
        #flux plot save
        plt.imshow(flux, cmap = 'hot')     #creates a color map of data, (data set, cmap = 'desired color code')
        plt.colorbar(label = '10-17 erg/s/cm2/Å/spaxel')
        plt.gca().invert_yaxis()
        plt.xlabel('spaxel')
        plt.ylabel('spaxel')
        plt.title(plateifu+' Flux')
        plt.savefig(data_folder + 'Plots/flux/'+plateifu+'_Flux.png')
        plt.close()
        
       
        #residual plot save
        residual = nmvel_map-mhalpha_vel
        '''print(nmvel_map[31,31])
        print(mhalpha_vel[31,31])
        print(residual[31,31])'''
        
        mresidual = ma.array(residual, mask=mask_f) #re-applied mask

        val_max_r = mresidual.max()       
        val_min_r = mresidual.min()

        if (val_max_r >= abs(val_min_r)):
            lim = val_max_r
        else:
            lim = abs(val_min_r)

        plt.imshow(mresidual,cmap = 'PiYG',vmin = -lim, vmax = lim)     #creates a color map of data, (data set, cmap = 'desired color code')
        plt.colorbar(label = 'residual(model-data)')  #following are just labels
        plt.gca().invert_yaxis()
        plt.xlabel('spaxel')
        plt.ylabel('spaxel')
        plt.title(plateifu + 'residual')
        plt.savefig(data_folder + 'Plots/residual/'+plateifu+'_residual.png')
        plt.close()
        
        #Rotation curve plot save
        
        radius_vel = []
        for i in range(map_shape[0]):
                for j in range(map_shape[1]):
                    if np.ma.is_masked(mhalpha_vel[i, j]):
                        continue
                   
                    radius = r[i,j]
                    vel = mhalpha_vel[i,j]
                    rot_vel = vel/(np.sin(i_angle_ideal)*np.cos(theta[i,j]))

                    if vel<0:
                        rad = -radius
                        fin_vel = -rot_vel
                    else:
                        rad = radius
                        fin_vel = rot_vel

                    radius_vel.append((rad, fin_vel))    


        r_list = [r for r, v in radius_vel]
        r_array = np.array(r_list)
        r_sorted = np.sort(r_array)


        plt.figure()
        plt.plot( r_sorted, (vmax_ideal*r_sorted)/((abs((Rturn_ideal**alpha_ideal)+(abs(r_sorted)**alpha_ideal)))**(1/alpha_ideal))+result.x[7], linestyle = '-', linewidth = '2', c = 'lime')
        plt.plot( *zip(*radius_vel), linestyle= '', marker = '.', c = 'black', markersize = '1')
        plt.ylim(-vmax_ideal-100, vmax_ideal+100)
        plt.xlabel("Deprojected Radius [kpc/h]")
        plt.ylabel("Rotational Velocity [km/s]")
        plt.title(plateifu+" Rotation curve ")
        plt.savefig(data_folder + 'Plots/rotation_curve/'+plateifu+'_rotation_curve.png')
        plt.close()

        #make subplot (data,model,flux,rotation curve)
        image_paths = [data_folder + 'Plots/Halpha/'+plateifu+'_Halpha_vel_plot.png',data_folder + 'Plots/model/'+plateifu+'_Halpha_vel_model.png',data_folder + 'Plots/flux/'+plateifu+'_Flux.png',data_folder + 'Plots/rotation_curve/'+plateifu+'_rotation_curve.png']
        
        
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize = (10,6)) 
        axs = axs.flatten()

        for i, file_path in enumerate(image_paths):
            img = mpimg.imread(file_path)
            axs[i].imshow(img)
            axs[i].axis('off')
        plt.tight_layout()
        plt.savefig(data_folder+ 'Plots/diagnostic_plots/'+plateifu+'_diagnostic.png')
        plt.close()



        ########################################################################
        # Stellar velocity fitting 
        #-----------------------------------------------------------------------
        stellar_vel_mask = np.logical_or(stellar_mask, binarea>5)                  #updated mask for stellar velocity,
        mstellar_vel = ma.array(stellar_vel, mask = stellar_vel_mask)              #using default mask and binned data points
        # mstellar_ivar = ma.array (stellar_vel_ivar, mask = stellar_vel_mask ) 
        mask_s = stellar_vel_mask
        
        
        #percent spaxel cut
        if (mstellar_vel.count()/mstellar_vel.size)*100<5:                         #same cut based on usable data points
            print(plate, ' not enough data for stellar fit', flush=True)
            output_tuple = (map_smoothness, vmax_ideal, alpha_ideal, Rturn_ideal, 
                        PA_ideal, i_angle_ideal, x_cent_ideal, y_cent_ideal, sys_vel_ideal, chi2_nu,
                        None, None, None, 
                        None, None, None, None, None, None, loc)
            return_queue.put(output_tuple)
            continue
            
        map_shape_s = mstellar_vel.shape

        # theta = np.radians(PA-90)                      #theta is PA reoriented to the positive x axis for calculations
        # for x in range(15,clean_coords[0]):       
        #     if x < 0 or x >= mstellar_vel.shape[0]:
        #         continue
            
        #     y = clean_coords[1]- round((clean_coords[0]-x) * np.tan(theta))               #y's are found based on the slope given by PA
        #     if y < 0 or y >= mstellar_vel.shape[1]:
        #         continue
                
        #     if ma.is_masked(mstellar_vel[x, y]):
        #         continue
        #     break
            
        # if (mstellar_vel[x,y]<0):         # if velocity comes back negative the position angle will be flipped 180 deg, otherwise left alone
        #     checkedPA_s = (PA + 180) *(np.pi/180)
        # else:
        #     checkedPA_s = PA*(np.pi/180)


        theta = np.radians(PA-90)   #theta is PA reoriented to the positive x axis for calculations
        checkedPA = np.radians(PA)
        f = 20
        
        try: 
            for x in range(int(clean_coords[0]/2),clean_coords[0]):   
                if x < 0 or x >= mstellar_vel.shape[0]:
                    continue
                
                y = clean_coords[1] - round((clean_coords[0] - x) * np.tan(theta)) #y's are found based on the slope given by PA

                if PA < 15 and PA > 0: # if position angle is near 0
                    x = clean_coords[0]
                    y = clean_coords[1] - f
                    while y < 0 or ma.is_masked(mstellar_vel[x,y]):
                        f -= 1
                        y = clean_coords[1] - f


                if y < 0 or y >= mstellar_vel.shape[1]:
                    continue
                    
                if ma.is_masked(mstellar_vel[x, y]):
                    continue
                
                else:
                    if (mstellar_vel[x,y]<0):
                        checkedPA = (PA + 180) *(np.pi/180)
                    else:
                        checkedPA = PA*(np.pi/180)
                    break

        except:
            print('couldnt check PA, using NSA value', flush=True)

        ########################################################################
        # Parameter guesses and bounds
        # ----------------------------------------------------------------------    

        vmax_guess_s = abs(mstellar_vel).max()/np.sin(i_angle)                             #set up initial guesses to start chi2 minimization using stellar map information
        alpha_guess_s = 1
        
        stellar_mx = np.max(mstellar_vel)                             
        s_mvel_indices = np.where(mstellar_vel == stellar_mx)                 
        s_vmax_coords = [int(s_mvel_indices[0][0]), int(s_mvel_indices[1][0])]
        
        dep_mvel_dist_stellar , _= deproject_spaxel(s_vmax_coords, clean_coords, checkedPA_s, i_angle)

        
        dmd_kpc_s = dep_mvel_dist_stellar*r_convert
        Rturn_guess_s = dmd_kpc_s/2      

        PA_guess_s = checkedPA_s
        i_angle_guess = i_angle
        center_guess = clean_coords
        #sys_vel_guess = mhalpha_vel[tuple(clean_coords)]
        sys_vel_guess =  0
        x0_s = [vmax_guess_s, alpha_guess_s, Rturn_guess_s, PA_guess_s, i_angle_guess, center_guess[0],center_guess[1],sys_vel_guess] #make list for guesses to be input as a list

        args_s = (mstellar_vel, stellar_vel_ivar, map_shape_s, r_convert, mask_s)    #non-variable values to go into minimzation function

        vmax_bound_s = (0, 1000)                                     #setting reasonable bounds for variable used in minimization
        alpha_bound_s = (0.001,100)
        Rturn_bound_s = (0.001,100)
        PA_bound_s = (0,2*np.pi)
        sys_vel_bound_s = (-100,100)
               
        i_angle_bound_s = (i_anglelow, i_anglehigh) 

        bounds_s = np.array([vmax_bound_s,alpha_bound_s,Rturn_bound_s,PA_bound_s,i_angle_bound,centeri_bound,centerj_bound, sys_vel_bound_s ]) #make list containing bounds to be sent into minimization function


        ########################################################################
        # fit stellar vel
        # ----------------------------------------------------------------------

        print(plate, ' fitting stellar vel', flush=True)
        start = datetime.datetime.now()

        try:
            result = minimize(     #using minimization function on chi2 function
                calc_chi2,
                x0_s,
                args = args_s,
                method = 'Powell',
                bounds = bounds_s
                ) 

        except:
            print(plate, 'CRASHED! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<', 
                            flush=True)
            output_tuple = (map_smoothness, vmax_ideal, alpha_ideal, Rturn_ideal, 
                PA_ideal, i_angle_ideal, x_cent_ideal, y_cent_ideal, sys_vel_ideal, chi2_nu,
                None, None, None, 
                None, None, None, None, None, None, loc)
            return_queue.put(output_tuple)
            continue

        fit_time = datetime.datetime.now() - start
        print(plate, " stellar velocity map fit", fit_time, flush=True)

        if result.success == False:
            print(plate, ' stellar fit failed', flush=True)
            output_tuple = (map_smoothness, vmax_ideal, alpha_ideal, Rturn_ideal, 
                        PA_ideal, i_angle_ideal, x_cent_ideal, y_cent_ideal, sys_vel_ideal, chi2_nu,
                        None, None, None, 
                        None, None, None, None, None, None, loc)
            return_queue.put(output_tuple)
            continue

        ########################################################################
        # Extract fit parameters
        # ----------------------------------------------------------------------

        #model based on chi2
        nmvel_map_stellar, r_stellar, theta_stellar = get_mvel_map(map_shape_s,                       #get_mvel_map is used to get the model vel map
                        result.x[0], 
                        result.x[1],
                        result.x[2],
                        result.x[3],
                        result.x[4],
                        result.x[5],
                        result.x[6],
                        result.x[7],
                        mask_s,r_convert
                       )
        
        #redefining result outputs
        vmax_ideal_s = result.x[0] 
        alpha_ideal_s = result.x[1]
        Rturn_ideal_s = result.x[2]
        PA_ideal_s = result.x[3]
        i_angle_ideal_s = result.x[4]
        x_cent_ideal_s = result.x[5]
        y_cent_ideal_s = result.x[6]
        sys_vel_ideal_s = result.x[7]

        chi2_nu_s = result.fun / (mstellar_vel.count() - 8)
        print('stellar map chi2_nu', chi2_nu_s, flush=True)
        

        ########################################################################
        # Stellar velocity plots
        #-----------------------------------------------------------------------
         # stellar plot save
        mstellar_vel = ma.array(stellar_vel, mask=mask_s) #re-applied mask

        val_max_s = mstellar_vel.max()        #simple plotting section for visualization
        val_min_s = mstellar_vel.min()
        

        if (val_max_s >= abs(val_min_s)):
            lim = val_max_s
        else:
            lim = abs(val_min_s) #gal_distance = (c*z)/H_0
    
    
        map_shape_s = mstellar_vel.shape
                
        
        #Plotting stellar velocity map
        plt.imshow(mstellar_vel,cmap = 'bwr',vmin = -lim, vmax = lim)     #creates a color map of data, (data set, cmap = 'desired color code')
        plt.colorbar(label = 'Velocity [km/s]')  #following are just labels
        plt.gca().invert_yaxis()
        plt.xlabel('spaxel')
        plt.ylabel('spaxel')
        plt.title(plateifu + ' Stellar Velocity (data)')
        plt.savefig(data_folder + 'Plots/stellar_plots/stellar_vel/'+plateifu+'_stellar_vel_plot.png')
        plt.close()
        
        #model plot save
        val_max_sm = nmvel_map_stellar.max()        #simple plotting section for visualization
        val_min_sm = nmvel_map_stellar.min()
        

        if (val_max_sm >= abs(val_min_sm)):
            lim = val_max_sm
        else:
            lim = abs(val_min_sm)

        plt.imshow(nmvel_map_stellar,cmap = 'bwr',vmin = -lim, vmax = lim)     #creates a color map of data, (data set, cmap = 'desired color code')
        plt.colorbar(label = 'Velocity [km/s]')  #following are just labels
        plt.gca().invert_yaxis()
        plt.xlabel('spaxel')
        plt.ylabel('spaxel')
        plt.title(plateifu+ ' Stellar Velocity (model)')
        plt.savefig(data_folder + 'Plots/stellar_plots/stellar_model/'+plateifu+'_stellar_vel_model.png')
        plt.close()
        
        #residual plot save
        residual_s = nmvel_map_stellar-mstellar_vel
        '''print(nmvel_map[31,31])
        print(mhalpha_vel[31,31])
        print(residual[31,31])'''
        
        mresidual_s = ma.array(residual_s, mask=mask_s) #re-applied mask

        val_max_sr = mresidual_s.max()        #simple plotting section for visualization
        val_min_sr = mresidual_s.min()
        

        if (val_max_sr >= abs(val_min_sr)):
            lim = val_max_sr
        else:
            lim = abs(val_min_sr)

        plt.imshow(mresidual_s,cmap = 'PiYG',vmin = -lim, vmax = lim)     #creates a color map of data, (data set, cmap = 'desired color code')
        plt.colorbar(label = 'residual(model-data)')  #following are just labels
        plt.gca().invert_yaxis()
        plt.xlabel('spaxel')
        plt.ylabel('spaxel')
        plt.title(plateifu + 'stellar residual')
        plt.savefig(data_folder + 'Plots/stellar_plots/stellar_residual/'+plateifu+'_stellar_residual.png')
        plt.close()
    
        
        #Rotation curve plot save
      
        radius_vel_s = []

        for i in range(map_shape_s[0]):
                for j in range(map_shape_s[1]):
                    if np.ma.is_masked(mstellar_vel[i, j]):
                        continue
                   
                    radius = r[i,j]
                    vel = mstellar_vel[i,j]
                    rot_vel = vel/(np.sin(i_angle_ideal_s)*np.cos(theta_stellar[i,j]))

                    if vel<0:
                        rad_s = -radius
                        fin_vel_s = -rot_vel
                    else:
                        rad_s = radius
                        fin_vel_s = rot_vel

                    radius_vel_s.append((rad_s, fin_vel_s))    


        r_list_s = [r for r, v in radius_vel_s]
        r_array_s = np.array(r_list_s)
        r_sorted_s = np.sort(r_array_s)


        plt.figure()
        plt.plot( r_sorted_s, (vmax_ideal_s*r_sorted_s)/((abs((Rturn_ideal_s**alpha_ideal_s)+(abs(r_sorted_s)**alpha_ideal_s)))**(1/alpha_ideal_s))+sys_vel_ideal_s, linestyle = '-', linewidth = '2', c = 'lime')
        plt.plot( *zip(*radius_vel_s), linestyle= '', marker = '.', c = 'black', markersize = '1')
        plt.ylim(-vmax_ideal_s-100, vmax_ideal_s+100)
        plt.xlabel("Deprojected Radius [kpc/h]")
        plt.ylabel("Rotational Velocity [km/s]")
        plt.title(plateifu+" Stellar rotation curve ")
        plt.savefig(data_folder + 'Plots/stellar_plots/stellar_RC/'+plateifu+'_stellar_rotation_curve.png')
        plt.close()


        stellar_image_paths = [data_folder + 'Plots/stellar_plots/stellar_vel/'+plateifu+'_stellar_vel_plot.png',data_folder + 'Plots/stellar_plots/stellar_model/'+plateifu+'_stellar_vel_model.png',data_folder + 'Plots/flux/'+plateifu+'_Flux.png',data_folder + 'Plots/stellar_plots/stellar_RC/'+plateifu+'_stellar_rotation_curve.png']
        
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize = (10,6)) 
        axs = axs.flatten()

        for i, file_path in enumerate(stellar_image_paths):
            img = mpimg.imread(file_path)
            axs[i].imshow(img)
            axs[i].axis('off')
        plt.tight_layout()
        plt.savefig(data_folder+ 'Plots/stellar_plots/stellar_diagnosis/'+plateifu+'_stellar_diagnostic.png')
        plt.close()


        ########################################################################
        print('\n', flush=True)

        output_tuple = (map_smoothness, vmax_ideal, alpha_ideal, Rturn_ideal, 
                        PA_ideal, i_angle_ideal, x_cent_ideal, y_cent_ideal, sys_vel_ideal, chi2_nu,
                        vmax_ideal_s, alpha_ideal_s, Rturn_ideal_s, 
                        PA_ideal_s, i_angle_ideal_s, x_cent_ideal_s, y_cent_ideal_s, sys_vel_ideal_s, chi2_nu_s,
                        loc)
    
        return_queue.put(output_tuple)



################################################################################
################################################################################
################################################################################

job_queue = Queue()
return_queue = Queue()


data_folder = "/scratch/ebenitez/"                  #create a variable for directory to common folder


drpall_fn = data_folder + "drpall_ttype_R90_chi2_halpha_stellar.fits"
drpall = Table.read(drpall_fn, format="fits",hdu=1)


drpall_dict = {}                                    #create the dictionary

for i in range(len(drpall)):                        #looping through drpall and redefining the index in terms
    plateifu = drpall['plateifu'][i]                #of the plateifu for simplicity			
    drpall_dict[plateifu]=i
    
print('dictionary created')


#get good galaxies
good_galaxy_folder = '/scratch/ebenitez/Plots/good' #Only want to use galaxies with good data so 
good_galaxies = []                                  #we fetch the good plateifu's and add them to 
                                                    #to our usable list

for filename in os.listdir(good_galaxy_folder):
        plateifu = filename.split('_')[0]
        if plateifu in drpall_dict:
            index = drpall_dict[plateifu]
            good_galaxies.append(index)
            
            
#create new column on table

custom_columns = ['smoothness_score',
    'vmax', 'alpha', 'Rturn', 'PA', 'i_angle', 'center_x', 'center_y', 'sys_vel', 'chi_nu',
    'stellar_vmax', 'stellar_alpha', 'stellar_Rturn', 'stellar_PA', 'stellar_i_angle',
    'stellar_center_x', 'stellar_center_y', 'stellar_chi_nu', 'stellar_sys_vel']


for col in custom_columns:
    if col not in drpall.colnames:
        drpall[col] = np.ones(len(drpall)) * -999

print('New columns created')


num_tasks = len(good_galaxies)

# Load jobs into queue
for i,index in enumerate(good_galaxies):
        
    job_queue.put(index)

    #if i > 10:
    #    num_tasks = 12
    #    break


print('Starting processes', datetime.datetime.now(), flush=True)

processes = []

for i in range(12):

    p = Process(target=process_1_galaxy, args=(job_queue, i, 
                     return_queue,  
                     drpall,
                     drpall_dict 
                     ))
    
    p.start()

    processes.append(p)

print('Populating output table', datetime.datetime.now(), flush=True)

################################################################################
# Iterate through the populated return queue to fill in the table
#-------------------------------------------------------------------------------
num_processed = 0

print(num_tasks)

while num_processed < num_tasks:

    try:
        return_tuple = return_queue.get(timeout=1.0)
    except:
        continue

    ############################################################################
    # Write the best-fit values and calculated parameters to a text file in 
    # ascii format.
    #---------------------------------------------------------------------------


    map_smoothness = return_tuple[0]

    vmax_ideal = return_tuple[1]
    alpha_ideal = return_tuple[2]
    Rturn_ideal = return_tuple[3]

    PA_ideal = return_tuple[4] 
    i_angle_ideal = return_tuple[5] 
    x_cent_ideal = return_tuple[6] 
    y_cent_ideal = return_tuple[7] 
    sys_vel_ideal = return_tuple[8] 

    chi_nu = return_tuple[9]

    vmax_ideal_s = return_tuple[10]
    alpha_ideal_s = return_tuple[11] 
    Rturn_ideal_s = return_tuple[12]

    PA_ideal_s = return_tuple[13] 
    i_angle_ideal_s = return_tuple[14] 
    x_cent_ideal_s = return_tuple[15] 
    y_cent_ideal_s = return_tuple[16] 
    sys_vel_ideal_s = return_tuple[17] 

    chi_nu_s = return_tuple[18]

    loc = return_tuple[19]


    if map_smoothness is not None:

        drpall['smoothness_score'][loc] = map_smoothness


    if vmax_ideal is not None:

        drpall['vmax'][loc] = vmax_ideal
        drpall['alpha'][loc] = alpha_ideal
        drpall['Rturn'][loc] =  Rturn_ideal


        drpall['PA'][loc] = PA_ideal
        drpall['i_angle'][loc] = i_angle_ideal
        drpall['center_x'][loc] = x_cent_ideal
        drpall['center_y'][loc] = y_cent_ideal

        drpall['chi_nu'][loc] = chi_nu

        drpall['sys_vel'][loc] = sys_vel_ideal


    if vmax_ideal_s is not None:

        drpall['stellar_vmax'][loc] = vmax_ideal_s
        drpall['stellar_alpha'][loc] = alpha_ideal_s
        drpall['stellar_Rturn'][loc] =  Rturn_ideal_s


        drpall['stellar_PA'][loc] = PA_ideal_s
        drpall['stellar_i_angle'][loc] = i_angle_ideal_s
        drpall['stellar_center_x'][loc] = x_cent_ideal_s
        drpall['stellar_center_y'][loc] = y_cent_ideal_s

        drpall['stellar_chi_nu'][loc] = chi_nu_s

        drpall['stellar_sys_vel'][loc] = sys_vel_ideal_s



        # put PARAMETER in table

    num_processed += 1

    if num_processed % 5 == 0:
        drpall.write(data_folder +'drpall_ttype_R90_chi2_halpha_stellar_run2.fits',
                format='fits', overwrite=True)
        print('Table written ', num_processed, flush=True)
    
    print(num_processed, ': ', num_processed)


print('Finished populating output table', datetime.datetime.now(), flush=True)
# Go through all the processes and join them back to the parent.
for p in processes:
    p.join(None)

drpall.write(data_folder +'drpall_ttype_R90_chi2_halpha_stellar_run2.fits',
                format='fits', overwrite=True)
print('Table written', flush=True)

FINISH = datetime.datetime.now()
print("Runtime:", FINISH - START, flush=True)
