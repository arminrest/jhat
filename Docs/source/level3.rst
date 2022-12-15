*************************
Level 3 From Aligned Cals
*************************

We run F560W with the best options determined from the :ref:`improve:Improving Alignment` example:

.. code-block:: python

	run_st_wcs_align_batch.py --input_dir '.' --input_files 'miri_example/*_cal.fits' --outrootdir aligned 
		--outsubdir F560W_level2_gaia  --overwr -v --refcat gaia --saveplots -tt -pp --histocut_order dxdy 
		--roundness1_lim -0.5 0.5 --objmag_lim 14 21.5  --refmag_lim 16 25 --delta_mag_lim -2 2 
		--d2d_max 1.5 --xshift 3 --yshift -4 --filter F560W

Now we run the notebook that runs level3 and creates the mosaic and the corresponding catalog 
(`On GitHub here <https://github.com/arminrest/jhat/blob/master/notebooks/st_align_Level3_example.ipynb>`_). 
That notebook produces a catalog (``F560W_snr3_npix10_cat.ecsv``), which we choose
as our secondary astrometric catalog.

Now we run the rest of the filters. However, we remove --delta_mag_lim, since this cut depends on the filter!!!
Note that we define the necessary column names from the catalog, which are different from the defaults. 

.. code-block:: python

	run_st_wcs_align_batch.py --input_dir '.' --input_files 'miri_example/*_cal.fits' --outrootdir miri_example 
		--outsubdir ALL_level2_catF560W  --overwr -v --refcat F560W_snr3_npix10_cat.ecsv --saveplots -tt -pp 
		--histocut_order dxdy --roundness1_lim -0.5 0.5 --objmag_lim 14 21.5  --refmag_lim 16 25 --d2d_max 1.5 
		--xshift 3 --yshift -4 --iterate_with_xyshifts --refcat_racol sky_centroid.ra 
		--refcat_deccol sky_centroid.dec --refcat_magcol aper50_abmag  --refcat_magerrcol aper50_abmag_err  
		--filters F560W F1000W F1280W F1130W F1500W F1800W


We find that F1500W F1800W mostly work, but don't have many stars. Therefore we run it with --d_rotated_Nsigma 0.0: 
Too few stars to do a 3-sigma cut. 

.. code-block:: python

	run_st_wcs_align_batch.py --input_dir '.' --input_files 'miri_example/*_cal.fits' --outrootdir miri_example 
		--outsubdir ALLRED_level2_catF560W  --overwr -vvv --refcat F560W_snr3_npix10_cat.ecsv --saveplots -tt -pp 
		--histocut_order dxdy --roundness1_lim -0.5 0.5 --objmag_lim 14 21.5  --refmag_lim 16 25 --d2d_max 1.5 --xshift 3 
		--yshift -4 --iterate_with_xyshifts --refcat_racol sky_centroid.ra --refcat_deccol sky_centroid.dec --refcat_magcol 
		aper50_abmag  --refcat_magerrcol aper50_abmag_err  --filters F1500W F1800W F2100W --d_rotated_Nsigma 0.0