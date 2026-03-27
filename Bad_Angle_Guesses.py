import numpy as np
import matplotlib.pyplot as plt
import astropy
from astropy.table import Table
from astropy.io import fits
import numpy.ma as ma
from PIL import Image
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
#custom
import sys 
sys.path.insert(1, '/gpfs/fs1/home/ssampat3/Desktop/MaNGA-VelMap')
from GenerateCutout import get_cutout,get_cutout_fits

GALAXY_DATA = "/scratch/kdougla7/data/SDSS/dr17/manga/spectro/analysis/v3_1_1/3.1.0/HYB10-MILESHC-MASTARSSP/"
EMILIO_DATA = "/scratch/ebenitez/"
DATA_FOLDER = "/scratch/ssampat3/"

with open("bad_angle_practice.txt", "r") as f:
    Plate_IFU = [line.strip() for line in f]

for test_galaxy2 in Plate_IFU:

    test_galaxy1 = test_galaxy2.replace("-", '/') + '/'
    cube_fn = GALAXY_DATA + test_galaxy1 + 'manga-' + test_galaxy2 + '-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz' 

    cube = fits.open(cube_fn)                       
    stellar_vel = cube['STELLAR_VEL'].data           
    stellar_vel_ivar = cube['STELLAR_VEL_IVAR'].data 
    stellar_mask = cube['STELLAR_VEL_MASK'].data
    halpha_vel = cube['EMLINE_GVEL'].data[23]
    halpha_gvel_ivar = cube['EMLINE_GVEL_IVAR'].data[23]
    halpha_gvel_mask = cube['EMLINE_GVEL_MASK'].data[23]
    flux = cube['SPX_MFLUX'].data
    cube.close()

    drpall = EMILIO_DATA + 'drpall_ttype_R90.fits'
    drpalltt = Table.read(drpall, format="fits",hdu=1)
    drpalltt [:5]
    drpalltt_dict = {}
    
    for i in range(len(drpalltt)):           
        plateifu = drpalltt['plateifu'][i] 
        drpalltt_dict[plateifu] = i 

    mstellar_vel = ma.array(stellar_vel, mask = stellar_mask)
    mhalpha_vel = ma.array(halpha_vel, mask = stellar_mask) 
    
    sm_val_max = mstellar_vel.max()
    sm_val_min = mstellar_vel.min()
    sm_vlim = max(np.abs([sm_val_max, sm_val_min]))
    
    hm_val_max = mhalpha_vel.max() 
    hm_val_min = mhalpha_vel.min() 
    hm_vlim = min(400, max(np.abs([hm_val_max, hm_val_min])))
    
    drpalltt_dict[test_galaxy2]
    index = drpalltt_dict[test_galaxy2]
    drpalltt[index]
    ra = drpalltt['objra'][index]
    dec = drpalltt['objdec'][index]
    r90 = drpalltt['R90'][index]
    directory = EMILIO_DATA + 'Plots/cutouts/'
    
    gal_image_jpg, wcs = get_cutout(test_galaxy2,ra,dec,r90,directory)
    img = np.asarray(Image.open(gal_image_jpg))
    center_coords = SkyCoord(ra*u.deg,dec*u.deg)
    pos_angle = drpalltt['nsa_elpetro_phi'][index]
    axis_ratio = drpalltt['nsa_sersic_ba'][index]
    
    point1 = center_coords.directional_offset_by(pos_angle*u.deg,r90*u.arcsec)
    alpha1 = point1.ra.deg
    delta1 = point1.dec.deg

    point2 = center_coords.directional_offset_by((pos_angle+180)*u.deg,r90*u.arcsec)
    alpha2 = point2.ra.deg
    delta2 = point2.dec.deg


    point3 = center_coords.directional_offset_by((pos_angle+90)*u.deg,(r90*axis_ratio)*u.arcsec)
    alpha3 = point3.ra.deg
    delta3 = point3.dec.deg


    point4 = center_coords.directional_offset_by((pos_angle-90)*u.deg,(r90*axis_ratio)*u.arcsec)
    alpha4 = point4.ra.deg
    delta4 = point4.dec.deg

    fig = plt.figure(figsize=(15,15))
    fig.suptitle(test_galaxy2)

    ax1 = fig.add_subplot(221)
    im = ax1.imshow(mstellar_vel,cmap = 'bwr', origin = 'lower', vmin = -sm_vlim, vmax = sm_vlim)
    ax1.set_xlabel('spaxel')
    ax1.set_ylabel('spaxel')
    ax1.set_title('Masked Stellar Velocity')
    fig.colorbar(im, ax = ax1,label = 'Velocity [km/s]' )

    #plotting lines on masked stellar velocity map
    ny, nx = mstellar_vel.shape
    xmin, xmax = ax1.get_xlim()
    ymin, ymax = ax1.get_ylim()
    x = np.linspace(xmin, xmax, 500)
    y1 = (x - nx/2) + ny/2
    y2 = -(x - nx/2) + ny/2
    ax1.axvline(x= nx/2, color = 'black', linewidth=2) #vertical line
    ax1.axhline(y= nx/2, color = 'black', linewidth=2) #horizontal line
    ax1.plot(x, y1, color = 'black', linewidth=2) #45 deg line
    ax1.plot(x, y2, color = 'black', linewidth=2) #135 deg line


    ax2 = fig.add_subplot(222)
    im = ax2.imshow(mhalpha_vel, cmap ='bwr', origin = 'lower', vmin = -hm_vlim, vmax = hm_vlim)    
    ax2.set_xlabel('spaxel') 
    ax2.set_ylabel('spaxel')
    ax2.set_title('Masked H-alpha Velocity')
    fig.colorbar(im, ax = ax2,label = 'Velocity [km/s]' )

    #plotting lines on masked H-alpha velocity map
    ax2.axvline(x= nx/2, color = 'black', linewidth=2) #vertical line
    ax2.axhline(y= nx/2, color = 'black', linewidth=2) #horizontal line
    ax2.plot(x, y1, color = 'black', linewidth=2) #45 deg line
    ax2.plot(x, y2, color = 'black', linewidth=2) #135 deg line


    ax3 = fig.add_subplot(223, projection=wcs)
    im = ax3.imshow(img[::-1])
    ax3.grid(color='white', ls='solid')
    ax3.set(xlabel='Galactic Longitude', ylabel='Galactic Latitude')
    ax3.plot([alpha1,alpha2],[delta1,delta2],color = 'limegreen',transform = ax3.get_transform('world'))
    ax3.plot([alpha3,alpha4],[delta3,delta4],color = 'orange',transform = ax3.get_transform('world'))
    ax3.plot(center_coords.ra.deg,center_coords.dec.deg,'rx',transform = ax3.get_transform('world') )


    ax4 = fig.add_subplot(224)
    im = ax4.imshow(flux, cmap ='hot', origin = 'lower', vmin = 0)
    ax4.set_title('Flux')
    ax4.set_xlabel('spaxel')
    ax4.set_ylabel('spaxel')
    fig.colorbar(im, ax = ax4,label = '10-17 erg/s/cm2/Å/spaxel')

    plt.savefig(DATA_FOLDER + 'bad_angle_practice/' + test_galaxy2, facecolor = "white")
