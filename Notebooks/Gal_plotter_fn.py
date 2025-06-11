import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
import numpy.ma as ma

def gal_plot(halpha_vel, stellar_vel, stellar_mask, flux, max_lim, img, plot_directory,alpha1, alpha2, alpha3, alpha4, delta1, delta2, delta3, delta4, center_coords, wcs, fig):

    mhalpha_vel = ma.array(halpha_vel, mask = stellar_mask)

  #  set equal limits for color bar
    val_max = mhalpha_vel.max()
    val_min = mhalpha_vel.min()
    if (val_max >= abs(val_min)):
        lim = val_max
    else:
        lim = abs(val_min)

    if (lim > max_lim):
        lim = max_lim

    #create subplot for mhalpha_vel
    ax2 = fig.add_subplot(222)
    im = ax2.imshow(mhalpha_vel, cmap ='bwr', origin = 'lower',vmax = lim, vmin = -lim)
    ax2.set_xlabel('spaxel')
    ax2.set_ylabel('spaxel')
    ax2.set_title('Masked H-alpha Velocity')


    fig.colorbar(im, ax = ax2,label = 'Velocity [km/s]' )
    #plt.close()


    mstellar_vel = ma.array(stellar_vel, mask = stellar_mask)

    #create subplot for mhalpha_vel
    val_max = mstellar_vel.max()
    val_min = mstellar_vel.min()
    if (val_max >= abs(val_min)):
        lim = val_max
    else:
        lim = abs(val_min)

    if (lim > max_lim):
        lim = max_lim
    #create subplot for mstellar_vel
    ax1 = fig.add_subplot(221)
    im = ax1.imshow(mstellar_vel, cmap ='bwr', origin = 'lower', vmax = lim, vmin = -lim )
    ax1.set_xlabel('spaxel')
    ax1.set_ylabel('spaxel')
    ax1.set_title('Masked Stellar Velocity')

    fig.colorbar(im, ax = ax1,label = 'Velocity [km/s]')

    #create subplot for flux

    ax4 = fig.add_subplot(224)
    im = ax4.imshow(flux, cmap ='hot', origin = 'lower', vmin = 0)
    ax4.set_title('Flux')
    ax4.set_xlabel('spaxel')
    ax4.set_ylabel('spaxel')



    fig.colorbar(im, ax = ax4,label = '10-17 erg/s/cm2/Å/spaxel' )

    ax3 = fig.add_subplot(223, projection=wcs)


    im = ax3.imshow(img[::-1])
    ax3.grid(color='white', ls='solid')
    ax3.set(xlabel='Right Ascension', ylabel='Declination')
    ax3.plot([alpha1,alpha2],[delta1,delta2],c = 'lime',transform = ax3.get_transform('world'))
    ax3.plot([alpha3,alpha4],[delta3,delta4],c = 'peru',transform = ax3.get_transform('world'))
    ax3.plot(center_coords.ra.deg,center_coords.dec.deg,'rx',transform = ax3.get_transform('world') )


    plt.savefig(plot_directory,bbox_inches='tight')
    plt.close(fig)
