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
#ref_data = fits.open(ref_image)['SCI',1].data
#norm1 = simple_norm(ref_data,stretch='log',min_cut=-1,max_cut=15)

# #plt.imshow(ref_data, origin='lower',
#                       #interval=MinMaxInterval(),
#                       norm=norm1,cmap='gray')
x = np.linspace(0, 2 * np.pi, 100)
y = np.sin(x)
# fig = plt.figure()
plt.plot(x, y)
plt.xlabel(r'$x$')
plt.ylabel(r'$\sin(x)$')
plt.show()