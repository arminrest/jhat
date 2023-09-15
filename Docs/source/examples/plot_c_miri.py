"""
=========
JWST MIRI
=========

Aligning JWST/MIRI images with JHAT.
"""
	
###############################################################
# An example MIRI Dataset is downloaded, and then a series of
# alignment methods are used. For more information on the
# key parameters used for alignment see 
# :ref:`params:Useful Parameters`.
   

import sys,os,glob
from astropy.io import fits
from astropy.table import Table
from astropy.nddata import extract_array
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy.wcs.utils import skycoord_to_pixel
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
from astroquery.mast import Observations
from astropy.visualization import (simple_norm,LinearStretch)

import jhat
from jhat import jwst_photclass,st_wcs_align


####################################################################
# 
# ------------------
# Relative Alignment
# ------------------
#
# **Download some Data**
#
# For this example we download 2 MIRI cal images from MAST. They're
# the same field and different filters. Note that 
# the code will also work for level 3 images.

obs_table1 = Observations.query_criteria(obs_id='jw02107-o038_t019_miri_f770w')
data_products_by_obs = Observations.get_product_list(obs_table1)
data_products_by_obs = data_products_by_obs[data_products_by_obs['calib_level']==2]
data_products_by_obs = data_products_by_obs[data_products_by_obs['productSubGroupDescription']=='CAL'][0]
Observations.download_products(data_products_by_obs,extension='fits')

obs_table2 = Observations.query_criteria(obs_id='jw02107-c1018_t019_miri_f1130w')
data_products_by_obs = Observations.get_product_list(obs_table2)
data_products_by_obs = data_products_by_obs[data_products_by_obs['calib_level']==2]
data_products_by_obs = data_products_by_obs[data_products_by_obs['productSubGroupDescription']=='CAL'][0]
Observations.download_products(data_products_by_obs,extension='fits')

####################################################################
# **Examine the Reference Image**
# 

files = glob.glob('mastDownload/JWST/*miri*/*cal.fits')
ref_image = files[0]
print(ref_image)
ref_fits = fits.open(ref_image)
ref_data = fits.open(ref_image)['SCI',1].data
norm1 = simple_norm(ref_data,stretch='log',min_cut=5,max_cut=25)

plt.imshow(ref_data, origin='lower',
                      norm=norm1,cmap='gray')
plt.gca().tick_params(labelcolor='none',axis='both',color='none')
plt.show()

####################################################################
# **Zoom in to see the offset**
# 
# Here add an artificial offset to the wcs, and then we see the 
# same star in both images at the same ra/dec
# location, demonstrating a large offset between
# the images.  
star_location = SkyCoord('23:09:44.0809','-43:26:05.613',unit=(u.hourangle,u.deg))
align_image = files[1]
align_fits = fits.open(align_image)
align_fits['SCI',1].header['CRPIX1']+=2
align_fits['SCI',1].header['CRPIX2']+=2
align_fits.writeto(align_image,overwrite=True)

align_data = fits.open(align_image)['SCI',1].data
ref_y,ref_x = skycoord_to_pixel(star_location,wcs.WCS(ref_fits['SCI',1],ref_fits))
align_y,align_x = skycoord_to_pixel(star_location,wcs.WCS(align_fits['SCI',1],align_fits))

ref_cutout = extract_array(ref_data,(11,11),(ref_x,ref_y))
align_cutout = extract_array(align_data,(11,11),(align_x,align_y))
norm1 = simple_norm(ref_cutout,stretch='log',min_cut=-1,max_cut=200)
norm2 = simple_norm(align_cutout,stretch='log',min_cut=-1,max_cut=200)
fig,axes = plt.subplots(1,2)
axes[0].imshow(ref_cutout, origin='lower',
                      norm=norm1,cmap='gray')
axes[1].imshow(align_cutout, origin='lower',
                      norm=norm2,cmap='gray')
axes[0].set_title('Reference')
axes[1].set_title('To Align')
axes[0].tick_params(labelcolor='none',axis='both',color='none')
axes[1].tick_params(labelcolor='none',axis='both',color='none')

plt.show()

####################################################################
# **Create a Photometric Catalog for Relative Alignment**
# 
# We choose one of the images to be the reference image, and then 
# create a catalog that we will use to align the other image.

jwst_phot = jwst_photclass()
jwst_phot.run_phot(imagename=ref_image,photfilename='auto',overwrite=True)
ref_catname = ref_image.replace('.fits','.phot.txt') # the default
refcat = Table.read(ref_catname,format='ascii')
print(refcat)

####################################################################
# **Align the second image**
# 
# The plots outputted here show the various steps used by jhat to
# determine the true matching sources in the image, and the
# subsequent correction needed for optimal alignment.

wcs_align = st_wcs_align()


wcs_align.run_all(align_image,
		  telescope='jwst',
		  outsubdir='mastDownload',
          refcat_racol='ra',
          refcat_deccol='dec',
          refcat_magcol='mag',
          refcat_magerrcol='dmag',
          overwrite=True,
          d2d_max=1,
          showplots=2,
          refcatname=ref_catname,
          histocut_order='dxdy',
              sharpness_lim=(0.3,0.9),
              roundness1_lim=(-0.7, 0.7),
              SNR_min= 3,
              dmag_max=1.0,
              objmag_lim =(14,24))


####################################################################
# **Check the Output**
# 
# The reference image has not changed, but let's read in the newly
# aligned image and compare with the original. 
# subsequent correction needed for optimal alignment.

aligned_image = os.path.join('mastDownload',os.path.basename(align_image).replace('cal.fits','jhat.fits'))
aligned_fits = fits.open(aligned_image)
aligned_data = fits.open(aligned_image)['SCI',1].data
aligned_y,aligned_x = skycoord_to_pixel(star_location,wcs.WCS(aligned_fits['SCI',1],aligned_fits))
aligned_cutout = extract_array(aligned_data,(11,11),(aligned_x,aligned_y))

norm3 = simple_norm(aligned_cutout,stretch='log',min_cut=-1,max_cut=200)
fig,axes = plt.subplots(1,3)
axes[0].imshow(ref_cutout, origin='lower',
                      norm=norm1,cmap='gray')
axes[1].imshow(align_cutout, origin='lower',
                      norm=norm2,cmap='gray')
axes[2].imshow(aligned_cutout, origin='lower',
                      norm=norm3,cmap='gray')
axes[0].set_title('Reference')
axes[1].set_title('To Align')
axes[2].set_title('Aligned')
for i in range(3):
	axes[i].tick_params(labelcolor='none',axis='both',color='none')


plt.show()