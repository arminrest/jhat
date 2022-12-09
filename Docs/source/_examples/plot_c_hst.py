"""
======
Hubble
======

Aligning HST images with JHAT.
"""
		
###############################################################
# An example HST Dataset is downloaded, and then a series of
# alignment methods are used.
   

import sys,os,glob
from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from astroquery.mast import Observations
from astropy.visualization import (simple_norm,LinearStretch)

import jhat
from jhat import hst_photclass,st_wcs_align


####################################################################
# **Download some Data**
#
# For this example we download 2 HST FLT images from MAST. They're
# the same filter and same field, just separated in time. Note that 
# the code will also work for drizzled images.

obs_table = Observations.query_criteria(obs_id='hst_16264_13_wfc3_ir_f160w_iebc13')
obs_table1 = obs_table[obs_table['filters']=='F160W']

obs_table = Observations.query_criteria(obs_id='hst_16264_14_wfc3_ir_f160w_iebc14')
obs_table2 = obs_table[obs_table['filters']=='F160W']

data_products_by_obs = Observations.get_product_list(obs_table1)
data_products_by_obs = data_products_by_obs[data_products_by_obs['calib_level']==2]
data_products_by_obs = data_products_by_obs[data_products_by_obs['productSubGroupDescription']=='FLT'][0]
Observations.download_products(data_products_by_obs,extension='fits')

data_products_by_obs = Observations.get_product_list(obs_table2)
data_products_by_obs = data_products_by_obs[data_products_by_obs['calib_level']==2]
data_products_by_obs = data_products_by_obs[data_products_by_obs['productSubGroupDescription']=='FLT'][0]
Observations.download_products(data_products_by_obs,extension='fits')

####################################################################
# **Examine the Reference Image**
# 

files = glob.glob('mastDownload/HST/*/*flt.fits')
ref_image = files[0]
align_image = files[1]
ref_data = fits.open(ref_image)['SCI',1].data
norm1 = simple_norm(ref_data,stretch='log',min_cut=-1,max_cut=15)

plt.imshow(ref_data, origin='lower',
                       #interval=MinMaxInterval(),
                       norm=norm1,cmap='gray')
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

wcs_align = st_wcs_align()
wcs_align.outdir = 'mastDownload'


wcs_align.run_all(align_image,
          refcat_racol='ra',
          refcat_deccol='dec',
          refcat_magcol='mag',
          refcat_magerrcol='dmag',
          overwrite=True,
          #Nbright=1000,
          #Nbright4match=1000,
          d2d_max=.12,
          showplots=2,
          #xshift=xshift,
          #yshift=yshift,
          #refcat_maincolor='gaia_g_rp',
                  refcatname=ref_catname,
                  histocut_order='dxdy',
                      sharpness_lim=(0.3,0.9),
                      roundness1_lim=(-0.7, 0.7),
                      #rough_cut_px_min=0.2,
                      #rough_cut_px_max=1.0,
                      SNR_min= 3,
                      dmag_max=1.0,
                      objmag_lim =(14,24))
