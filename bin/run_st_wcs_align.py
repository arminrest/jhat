#!/usr/bin/env python3

from jhat import st_wcs_align


wcs_align = st_wcs_align()
parser = wcs_align.define_options()
args = parser.parse_args()

wcs_align.verbose=args.verbose
wcs_align.replace_sip = args.replace_sip
wcs_align.sip_err = args.sip_err
wcs_align.sip_degree = args.sip_degree
wcs_align.sip_points = args.sip_points
wcs_align.rough_cut_px_min = args.rough_cut_px_min
wcs_align.rough_cut_px_max = args.rough_cut_px_max
wcs_align.d_rotated_Nsigma = args.d_rotated_Nsigma


wcs_align.run_all(args.cal_image,
                 telescope = args.telescope,
                 #distortion_file = args.distortion_file,
                 outrootdir = args.outrootdir,
                 outsubdir = args.outsubdir,
                 distortion_file = args.distortion_file,
                 overwrite = args.overwrite,
                 photometry_method = args.photometry_method,
                 find_stars_threshold = args.find_stars_threshold,
                 sci_xy_catalog=args.sci_xy_catalog,
                 #skip_applydistortions_if_exists=args.skip_applydistortions_if_exists,
                 use_dq = not args.skip_use_dq,
                 refcatname = args.refcat,
                 refcat_racol = args.refcat_racol,
                 refcat_deccol = args.refcat_deccol,
                 refcat_magcol = args.refcat_magcol,
                 refcat_magerrcol = args.refcat_magerrcol,
                 refcat_colorcol = args.refcat_colorcol,
                 pmflag = args.refcat_pmflag,
                 pm_median = args.refcat_pmmedian,
                 photfilename = args.photfilename,
                 load_photcat_if_exists=args.load_photcat_if_exists,
                 rematch_refcat=args.rematch_refcat,
                 SNR_min = args.SNR_min,
                 d2d_max = args.d2d_max, # maximum distance refcat to source in image
                 dmag_max = args.dmag_max, # maximum uncertainty of source 
                 sharpness_lim = args.sharpness_lim, # sharpness limits
                 roundness1_lim = args.roundness1_lim, # roundness1 limits 
                 delta_mag_lim =  args.delta_mag_lim, # limits on mag-refcat_mainfilter
                 objmag_lim = args.objmag_lim, # limits on mag, the magnitude of the objects in the image
                 refmag_lim = args.refmag_lim, # limits on refcat_mainfilter, the magnitude of the reference catalog
                 slope_min = args.slope_min,
                 Nbright4match=args.Nbright4match, # Use only the the brightest  Nbright sources from image for the matching with the ref catalog
                 Nbright=args.Nbright,    # U/se only the brightest Nbright sources from image
                 histocut_order=args.histocut_order, # histocut_order defines whether the histogram cut is first done for dx or first for dy
                 xshift=args.xshift,# added to the x coordinate before calculating ra,dec (only impacts ra,dec, not x). This can be used to correct for large shifts before matching!
                 yshift=args.yshift, # added to the y coordinate before calculating ra,dec (only impacts ra,dec, not y). This can be used to correct for large shifts before matching!
                 iterate_with_xyshifts=args.iterate_with_xyshifts, # After the first histogram fit, redo the match with refcat with x/yshift=median(dx/dy) and redo histofit. Use this if the offsets are big, since the second iteration will give you better matching with the refcat 
                 showplots=args.showplots,
                 saveplots=args.saveplots,# 
                 savephottable=args.savephottable,
                 ee_radius=args.ee_radius
                 )