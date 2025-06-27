import os
from astropy.wcs import WCS
import requests
from astropy.io import fits


def get_cutout(plateifu, ra, dec, r90, cache_dir, verbose=False):
    '''Take MaNGA galaxy ra, dec, r90 and generates a cutout from the legacy survey viewer
    with size 4*r90 by 4*r90. cutout is saved at cache_dir with file name plateifu.jpg.
    
    Parameters
    ----------
    plateifu : str
        plateifu for galaxy.
    ra : float
        Right ascension (degrees).
    dec : float
        Declination (degrees).
    r90 : float
        r90 in arcsec for galaxy, to be used for cutout size
    cache_dir : string
        cache location
    verbose : bool
        Add some status messages if true.
        
    Returns
    -------
    img_name : str
        Name of JPG cutout file written after query.
    w : astropy.wcs.WCS
        World coordinate system for the image.
    '''
    
    # Either load an existing image or download a cutout.
    img_name = cache_dir + plateifu + '.jpg'

    
    size = int(4 * r90 / 0.262)
    
    if os.path.exists(img_name):
        if verbose:
            print('{} exists.'.format(img_name))

    else:
        img_url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra={}&dec={}&zoom=14&size={}'.format(ra, dec, size)
        if verbose:
            print('Get {}'.format(img_url))
            
        with open(img_name, 'wb') as handle: 
            response = requests.get(img_url, stream=True) 
            if not response.ok: 
                print(response) 
            for block in response.iter_content(1024): 
                if not block: 
                    break 
                handle.write(block)
                
    # Set up the WCS.

    
    wcs_input_dict = {
        'CTYPE1': 'RA---TAN',
        'CUNIT1': 'deg',
        'CDELT1': -0.262/3600,
        'CRPIX1': size/2 + 0.5,
        'CRVAL1': ra,
        'NAXIS1': size,
        'CTYPE2': 'DEC--TAN',
        'CUNIT2': 'deg',
        'CDELT2': 0.262/3600,
        'CRPIX2': size/2 + 0.5,
        'CRVAL2': dec,
        'NAXIS2': size
    }
    w = WCS(wcs_input_dict)
    
    return img_name, w


def get_cutout_fits(plateifu, ra, dec, r90, cache_dir, verbose=False):
    '''Take MaNGA galaxy ra, dec, r90 and generates a cutout from the legacy survey viewer
    with size 4*r90 by 4*r90. cutout is saved at cache_dir with file name plateifu.jpg.
    
    Parameters
    ----------
    plateifu : str
        plateifu for galaxy.
    ra : float
        Right ascension (degrees).
    dec : float
        Declination (degrees).
    r90 : float
        r90 in arcsec for galaxy, to be used for cutout size
    cache_dir : string
        cache location
    verbose : bool
        Add some status messages if true.
        
    Returns
    -------
    img_name : str
        Name of fits cutout file written after query.
    w : astropy.wcs.WCS
        World coordinate system for the image.
    '''
    
    # Either load an existing image or download a cutout.
    img_name = cache_dir + plateifu + '.fits'

    
    size = int(4 * r90 / 0.262)
    
    if os.path.exists(img_name):
        hdu = fits.open(img_name, format = 'fits')
        if verbose:
            print('{} exists.'.format(img_name))

    else:
        img_url = 'https://www.legacysurvey.org/viewer/cutout.fits?ra={}&dec={}&zoom=14&size={}&layer=ls-dr10'.format(ra, dec, size)
        if verbose:
            print('Get {}'.format(img_url))
            
        hdu = fits.open(img_url, format = 'fits')
        hdu.writeto(img_name) 
            
                
    
    w = WCS(hdu[0].header)
    
    return img_name, w