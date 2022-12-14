"""
======
Hubble
======

Aligning HST images with JHAT.
"""
	
###############################################################
# An example HST Dataset is downloaded, and then a series of
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
from jhat import hst_photclass,st_wcs_align


####################################################################
# 
# ------------------
# Relative Alignment
# ------------------
#
# **Download some Data**
#
# For this example we download 2 HST DRZ images from MAST. They're
# the same filter and same field, just separated in time. 

obs_table = Observations.query_criteria(obs_id='hst_16264_12_wfc3_ir_f110w_iebc12')
obs_table1 = obs_table[obs_table['filters']=='F110W']

obs_table = Observations.query_criteria(obs_id='hst_16264_15_wfc3_ir_f110w_iebc15')
obs_table2 = obs_table[obs_table['filters']=='F110W']

data_products_by_obs = Observations.get_product_list(obs_table1)
data_products_by_obs = data_products_by_obs[data_products_by_obs['calib_level']==3]
data_products_by_obs = data_products_by_obs[data_products_by_obs['productSubGroupDescription']=='DRZ'][0]
Observations.download_products(data_products_by_obs,extension='fits')

data_products_by_obs = Observations.get_product_list(obs_table2)
data_products_by_obs = data_products_by_obs[data_products_by_obs['calib_level']==3]
data_products_by_obs = data_products_by_obs[data_products_by_obs['productSubGroupDescription']=='DRZ'][0]
Observations.download_products(data_products_by_obs,extension='fits')

####################################################################
# **Examine the Reference Image**
# 

files = glob.glob('mastDownload/HST/*/*drz.fits')
ref_image = files[0]
ref_fits = fits.open(ref_image)
ref_data = fits.open(ref_image)['SCI',1].data
norm1 = simple_norm(ref_data,stretch='log',min_cut=-1,max_cut=15)

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

star_location = SkyCoord('21:29:40.5351','+0:04:42.697',unit=(u.hourangle,u.deg))
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

hst_phot = hst_photclass(psf_fwhm=1.8,aperture_radius=5)
hst_phot.run_phot(imagename=ref_image,photfilename='auto',overwrite=True)
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
		  telescope='hst',
		  outsubdir='mastDownload',
          refcat_racol='ra',
          refcat_deccol='dec',
          refcat_magcol='mag',
          refcat_magerrcol='dmag',
          overwrite=True,
          d2d_max=.5,
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

aligned_image = os.path.join('mastDownload',os.path.basename(align_image).replace('.fits','.jhat.fits'))
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


####################################################################
# 
# -------------
# Align to Gaia
# -------------
#
# You can also align each image to the Gaia DR3 catalog, or you
# could replace the catalog created in step one with your own
# catalog of the field. 


wcs_align.run_all(align_image,
		  telescope='hst',
		  outsubdir='mastDownload',
          overwrite=True,
          d2d_max=.5,
          showplots=0,
          refcatname='Gaia',
          histocut_order='dxdy',
              sharpness_lim=(0.3,0.9),
              roundness1_lim=(-0.7, 0.7),
              SNR_min= 3,
              dmag_max=1.0,
              objmag_lim =(14,24))

aligned_image = os.path.join('mastDownload',os.path.basename(align_image).replace('.fits','.jhat.fits'))
aligned_fits = fits.open(aligned_image)
aligned_data = fits.open(aligned_image)['SCI',1].data
aligned_y,aligned_x = skycoord_to_pixel(star_location,wcs.WCS(aligned_fits['SCI',1],aligned_fits))
aligned_cutout = extract_array(aligned_data,(11,11),(aligned_x,aligned_y))

norm3 = simple_norm(aligned_cutout,stretch='log',min_cut=-1,max_cut=200)
fig,axes = plt.subplots(1,2)
axes[0].imshow(align_cutout, origin='lower',
                      norm=norm2,cmap='gray')
axes[1].imshow(aligned_cutout, origin='lower',
                      norm=norm3,cmap='gray')
axes[0].set_title('To Align')
axes[1].set_title('Aligned')
for i in range(2):
	axes[i].tick_params(labelcolor='none',axis='both',color='none')


plt.show()

####################################################################
# 
# -------------
# Large Offsets
# -------------
#
# Sometimes the initial images are so poorly aligned, that the code
# fails. Here we read in the same image as in the first example,
# and add an additional 3 pixel offset in the wcs. 

files = glob.glob('mastDownload/HST/*/*drz.fits')
align_image = files[1]
align_fits = fits.open(align_image)
align_fits['SCI',1].header['CRPIX1']+=3
align_fits['SCI',1].header['CRPIX2']+=3
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

wcs_align = st_wcs_align()

try:
	wcs_align.run_all(align_image,
		  telescope='hst',
		  outsubdir='mastDownload',
          refcat_racol='ra',
          refcat_deccol='dec',
          refcat_magcol='mag',
          refcat_magerrcol='dmag',
          overwrite=True,
          d2d_max=.5,
          showplots=2,
          refcatname=ref_catname,
          histocut_order='dxdy',
              sharpness_lim=(0.3,0.9),
              roundness1_lim=(-0.7, 0.7),
              SNR_min= 3,
              dmag_max=1.0,
              objmag_lim =(14,24))

except:
	print('Failed for not enough matches!')

####################################################################
# 
# This is what a failure looks like (compare to the plots above).
# There are now a couple of options here. You can increase the 
# d2d_max parameter, which increases the allowed distance between 
# sources being matched in the reference and target images:

wcs_align = st_wcs_align()


wcs_align.run_all(align_image,
		  telescope='hst',
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

aligned_image = os.path.join('mastDownload',os.path.basename(align_image).replace('.fits','.jhat.fits'))
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

####################################################################
# 
# Or you can apply a rough guess for the offset, and then use a
# smaller d2d_max for matching:

wcs_align = st_wcs_align()


wcs_align.run_all(align_image,
		  telescope='hst',
		  outsubdir='mastDownload',
          refcat_racol='ra',
          refcat_deccol='dec',
          refcat_magcol='mag',
          refcat_magerrcol='dmag',
          overwrite=True,
          d2d_max=.25,
          xshift=5,
          yshift=5,
          showplots=2,
          refcatname=ref_catname,
          histocut_order='dxdy',
              sharpness_lim=(0.3,0.9),
              roundness1_lim=(-0.7, 0.7),
              SNR_min= 3,
              dmag_max=1.0,
              objmag_lim =(14,24))

aligned_image = os.path.join('mastDownload',os.path.basename(align_image).replace('.fits','.jhat.fits'))
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