
from PIL import Image
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS

        
        
def make_point(ra,dec,pos_angle,r90,axis_ratio):
    center_coords = SkyCoord(ra*u.deg,dec*u.deg)

    point1 = center_coords.directional_offset_by(pos_angle*u.deg,(r90*1.5)*u.arcsec)

    alpha1 = point1.ra.deg
    delta1 = point1.dec.deg


    point2 = center_coords.directional_offset_by((pos_angle+180)*u.deg,(r90*1.5)*u.arcsec)

    alpha2 = point2.ra.deg
    delta2 = point2.dec.deg

    point3 = center_coords.directional_offset_by((pos_angle+90)*u.deg,((r90*1.5)*axis_ratio)*u.arcsec)

    alpha3 = point3.ra.deg
    delta3 = point3.dec.deg


    point4 = center_coords.directional_offset_by((pos_angle-90)*u.deg,((r90*1.5)*axis_ratio)*u.arcsec)

    alpha4 = point4.ra.deg
    delta4 = point4.dec.deg

    return alpha1, alpha2, alpha3, alpha4, delta1, delta2, delta3, delta4, center_coords
