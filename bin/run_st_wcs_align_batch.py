#!/usr/bin/env python

from jhat import align_wcs_batch
import sys,os

usagestring = """JHAT Batch file Alignment

usage: see docs

"""


align_batch = align_wcs_batch()
parser = align_batch.define_options()
args = parser.parse_args()

align_batch.verbose=args.verbose
align_batch.debug=args.debug
align_batch.replace_sip = args.replace_sip
align_batch.sip_err = args.sip_err
align_batch.sip_degree = args.sip_degree
align_batch.sip_points = args.sip_points
 
align_batch.rough_cut_px_min = args.rough_cut_px_min
align_batch.rough_cut_px_max = args.rough_cut_px_max
align_batch.d_rotated_Nsigma = args.d_rotated_Nsigma

# get the input files
align_batch.get_input_files(args.input_files,directory=args.input_dir,
                            detectors=args.detectors,
                            filters=args.filters,pupils=args.pupils)

ixs_all = align_batch.getindices()

if len(ixs_all)==0:
    print('NO IMAGES FOUND!! exiting...')
    sys.exit(0)

if args.distortioncoeffs_dir is not None:
    align_batch.get_distortioncoeff_files(args.distortioncoeffs_dir)
    align_batch.match_distortioncoeffs(ixs_all)


    
# get the output filenames
ixs_exists,ixs_notexists = align_batch.get_output_filenames(ixs=ixs_all,
                                                            outrootdir=args.outrootdir,
                                                            outsubdir=args.outsubdir,
                                                            addfilter2outsubdir=args.addfilter2outsubdir)    


ixs_todo = ixs_notexists[:]
if len(ixs_exists)>0:
    if args.skip_if_exists:
        print(f'{len(ixs_exists)} output images already exist, skipping since --skip_if_exists')
    else:
        if args.overwrite:
           ixs_todo.extend(ixs_exists) 
           print(f'{len(ixs_exists)} output images already exist,overwriting them if continuing!')
        else:
           raise RuntimeError(f'{len(ixs_exists)} output images already exist, exiting! if you want to overwrite them, use the --overwrite option, or if you want to skip them, use the --skip_if_exists option!')

if len(ixs_todo)==0:
    print(f'There are {len(ixs_all)} images, but none of them need to be done!')
    sys.exit(0)

print(f'Output directory:{os.path.dirname(align_batch.t.loc[ixs_todo[0],"outfilename"])}')
do_it = input(f'Do you want to continue and align the wcs for {len(ixs_todo)} images [y/n]?  ')
if do_it.lower() in ['y','yes']:
    pass
elif do_it.lower() in ['n','no']:
    print('OK, stopping....')
    sys.exit(0)
else:
    print(f'Hmm, \'{do_it}\' is neither yes or no. Don\'t know what to do, so stopping ....')
    sys.exit(0)

align_batch.align_wcs(ixs_todo,
                    overwrite = args.overwrite,
                    outrootdir= args.outrootdir,
                    outsubdir = args.outsubdir,
                    addfilter2outsubdir = args.addfilter2outsubdir,
                    telescope = args.telescope,
                    #skip_applydistortions_if_exists=args.skip_applydistortions_if_exists,
                    photometry_method = args.photometry_method,
                    find_stars_threshold = args.find_stars_threshold,
                    sci_xy_catalog=args.sci_xy_catalog,
                    use_dq = not args.skip_use_dq,
                    refcatname = args.refcat,
                    refcat_racol = args.refcat_racol,
                    refcat_deccol = args.refcat_deccol,
                    refcat_magcol = args.refcat_magcol,
                    refcat_magerrcol = args.refcat_magerrcol,
                    refcat_colorcol = args.refcat_colorcol,
                    pmflag = args.refcat_pmflag,
                    pm_median = args.refcat_pmmedian,
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
                    xshift=args.xshift, # added to the x coordinate before calculating ra,dec. This can be used to correct for large shifts before matching!
                    yshift=args.yshift, # added to the y coordinate before calculating ra,dec. This can be used to correct for large shifts before matching!
                    iterate_with_xyshifts=args.iterate_with_xyshifts, # After the first histogram fit, redo the match with refcat with x/yshift=median(dx/dy) and redo histofit. Use this if the offsets are big, since the second iteration will give you better matching with the refcat 
                    showplots=0,
                    saveplots=args.saveplots,# 
                    savephottable=args.savephottable
                    )

align_batch.write()