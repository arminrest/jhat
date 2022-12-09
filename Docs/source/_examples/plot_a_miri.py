"""
=========
JWST MIRI
=========

Aligning JWST MIRI images with JHAT.
"""
		
###############################################################
# Explain here what this example does
   
import numpy as np
import matplotlib.pyplot as plt
import sys,os,glob
from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from astroquery.mast import Observations
from astropy.visualization import (simple_norm,LinearStretch)

import jhat
from jhat import hst_photclass,st_wcs_align

x = np.linspace(0, 2 * np.pi, 100)
y = np.sin(x)

plt.plot(x, y)
plt.xlabel(r'$x$')
plt.ylabel(r'$\sin(x)$')
# To avoid matplotlib text output
plt.show()