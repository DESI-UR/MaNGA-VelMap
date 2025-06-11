import matplotlib
matplotlib.use("Agg")

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from PIL import Image
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
import traceback
import numpy.ma as ma
from astropy.wcs import WCS

from GenerateCutout import get_cutout
from MakePoint import make_point
from Gal_plotter_fn import gal_plot



data_folder = "/scratch/ebenitez/" #create a variable for directory to common folder


drpall_fn = data_folder + "drpall_ttype_R90.fits"
drpall = Table.read(drpall_fn, format="fits",hdu=1)


drpall_dict = {}                    #create the dictionary

for i in range(len(drpall)):           #loop that repeats for the length of drpalltt file
    plateifu = drpall['plateifu'][i]    #variable that holds the plateifu withing drpalltt at i(current row) 			
    drpall_dict[plateifu]=i


max_lim = 400              #creates default plot limit

for i in range(9080,9180):

    plateifu = drpall['plateifu'][i]      #calls the specific row by index

    plot_directory = data_folder + 'Plots/2x2_plots/'+ plateifu +'_2x2.png'
    plate=plateifu.split('-',1)
    plate1=plate[0]
    plate2=plate[1]
    try: 
        cube_df = '/scratch/kdougla7/data/SDSS/dr17/manga/spectro/analysis/v3_1_1/3.1.0/HYB10-MILESHC-MASTARSSP/'+plate1+'/'+plate2+'/manga-'+plateifu+'-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz'
        cube = fits.open(cube_df)                        #file requires opening
        stellar_vel = cube['STELLAR_VEL'].data           #for the following we create variable representing data from file HDU's
        stellar_vel_ivar = cube['STELLAR_VEL_IVAR'].data 
        stellar_mask = cube['STELLAR_VEL_MASK'].data
        halpha_vel = cube['EMLINE_GVEL'].data[23]        #the following three have several channels, but were only concerned with the Halpha spectra (channel 23)
        halpha_gvel_ivar = cube['EMLINE_GVEL_IVAR'].data[23]
        halpha_gvel_mask = cube['EMLINE_GVEL_MASK'].data[23]
        flux = cube['SPX_MFLUX'].data

        ttype = drpall['TType'][i]
        ttype>0

        cube.close()
    #creates figure for subplots,  
        fig = plt.figure( figsize=(15,15))
        fig.suptitle(plateifu)
        ra = drpall['objra'][i]
        dec = drpall['objdec'][i]
        r90 = drpall['R90'][i]
        directory = data_folder + 'Plots/cutouts/'
        gal_image, wcs = get_cutout(plateifu,ra,dec,r90,directory) 
       
        img = np.asarray(Image.open(data_folder + 'Plots/cutouts/'+plateifu+'.jpg'))
        index = drpall_dict[plateifu]
        pos_angle = drpall['nsa_elpetro_phi'][index]          #position angle is in degrees, this converts to radians bc math expects radians
        axis_ratio = drpall['nsa_sersic_ba'][index]
        
        alpha1, alpha2, alpha3, alpha4, delta1, delta2, delta3, delta4, center_coords = make_point(ra,dec,pos_angle,r90,axis_ratio) 
        
        gal_plot(halpha_vel, stellar_vel, stellar_mask, flux, max_lim, img, plot_directory, alpha1, alpha2, alpha3, alpha4, delta1, delta2, delta3, delta4, center_coords, wcs, fig)
        
    except Exception as e:
        print(f'Error with {plateifu}: {e}')
        traceback.print_exc()
        continue
        
