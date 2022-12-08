#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 09:39:07 2022

@author: arest
"""

import argparse,sys,os,re
from astropy.io import fits
import glob

from .pdastro import pdastroclass,unique
from .st_wcs_align import st_wcs_align
from .simple_jwst_phot import jwst_photclass

class align_wcs_batch(pdastroclass):
    def __init__(self):
        pdastroclass.__init__(self)

        self.verbose=0
        self.debug = False
        
        self.outdir = None
        
        
        self.aperture_col = 'AperName'
        self.detector_col = 'detector'
        self.filter_col = 'filter'
        self.pupil_col = 'pupil'

        
        self.wcs_align = st_wcs_align()

    def define_options(self,parser=None,usage=None,conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)
        
        # default directory for input images
        if 'JWST_INPUT_IMAGEDIR' in os.environ:
            inputdir = os.environ['JWST_INPUT_IMAGEDIR']
        else:
            inputdir = None

        # default directory for output
        if 'JWST_OUTROOTDIR' in os.environ:
            outrootdir = os.environ['JWST_OUTROOTDIR']
        else:
            outrootdir = None


        parser.add_argument('--input_dir', default=inputdir, help='Directory in which the cal or rate images are located, to which the distortions are applied. (default=%(default)s)')
        parser.add_argument('--input_files', nargs='+', default=['*_cal.fits'], help='list of cal or rate file(pattern)s to which the distortion files are applied to. "input_dir" is used if not None (default=%(default)s)')

        parser.add_argument('--outrootdir', default=outrootdir, help='output root directory. The output directoy is the output root directory + the outsubdir if not None (default=%(default)s)')
        parser.add_argument('--outsubdir', default=None, help='outsubdir added to output root directory (default=%(default)s)')
        parser.add_argument('--overwrite', default=False, action='store_true', help='overwrite files if they exist.')

        parser.add_argument('--skip_if_exists', default=False, action='store_true', help='Skip doing the analysis of a given input image if the cal file already exists, assuming the full analysis has been already done')

        parser.add_argument('-v','--verbose', default=0, action='count')

        parser.add_argument('--apertures', nargs='+', default=None, help='constrain the input file list to these apertures (default=%(default)s)')
        parser.add_argument('--filters', nargs='+', default=None, help='constrain the input file list to these filters (default=%(default)s)')
        parser.add_argument('--pupils', nargs='+', default=None, help='constrain the input file list to these pupils (default=%(default)s)')

        
        parser.add_argument('-d','--debug', default=False, action='store_true',help='debug mode: alignment is done outside "try" block!')

        parser = self.wcs_align.default_options(parser)

        return(parser)
    
    #def get_output_filenames(self, suffixmapping = {'cal':'tweakregstep','rate':'tweakregstep'},ixs=None):
    def get_output_filenames(self, outputsuffix = 'tweakregstep',ixs=None):
        ixs = self.getindices(ixs)
        ixs_exists=[]
        ixs_notexists=[]
        for ix in ixs:
            pattern = '(.*)_([a-zA-Z0-9]+)\.fits$'
            m = re.search(pattern,self.t.loc[ix,'filename'])
            if m is None:
                raise RuntimeError(f'cannot determine suffix for file {self.t.loc[ix,"filename"]} with pattern {pattern}!')
            (inputbasename,inputsuffix)=m.groups()
            inputbasename = os.path.basename(inputbasename)
            
            #outfilename = f'{self.outdir}/{inputbasename}_{suffixmapping[inputsuffix]}.fits'
            outfilename = f'{self.outdir}/{inputbasename}_{outputsuffix}.fits'
            self.t.loc[ix,'outfilename']=outfilename
            if os.path.isfile(outfilename):
                self.t.loc[ix,'outfile_exists']='yes'
                ixs_exists.append(ix)
            else:
                self.t.loc[ix,'outfile_exists']='no'
                ixs_notexists.append(ix)
        
        if self.verbose>1:  self.write(indices=ixs)

        return(ixs_exists,ixs_notexists)


    def set_outdir(self,outrootdir=None,outsubdir=None):
        self.outdir = self.wcs_align.set_outdir(outrootdir,outsubdir)
        return(self.outdir)

    def get_files(self,filepatterns,directory=None):
        filenames=[]
        for filepattern in filepatterns:
            #(tmpdir,basename) = os.path.split(filepattern)
            #print(f'{tmpdir},{basename}')
            #if tmpdir =='' and (directory is not None):
            if directory is not None:
                filepattern=os.path.join(directory,filepattern)
            if self.verbose>2: print(f'Looking for filepattern {filepattern}')
            filenames.extend(glob.glob(filepattern))
        
        for i in range(len(filenames)):
            filenames[i] = os.path.abspath(filenames[i])
        if self.verbose: print(f'Found {len(filenames)} files matching filepatterns {filepatterns}')
        return(filenames)
    
    
    def get_input_files(self,filepatterns,directory=None,filters=None,pupils=None):
        self.t['filename'] = self.get_files(filepatterns,directory=directory)
        ixs = self.getindices()
        for ix in ixs:
            hdr = fits.getheader(self.t.loc[ix,'filename'])
            detector = re.sub('long$','5',hdr['DETECTOR'].lower())
            self.t.loc[ix,self.aperture_col]=f'{detector}_{hdr["SUBARRAY"].lower()}'
            self.t.loc[ix,'detector']=f'{detector}'
            self.t.loc[ix,'instrument']=f'{hdr["INSTRUME"].lower()}'
            self.t.loc[ix,'subarray']=f'{hdr["SUBARRAY"].lower()}'
            if self.t.loc[ix,'instrument']=='fgs':
                self.t.loc[ix,self.filter_col]='clear'
                self.t.loc[ix,self.pupil_col]='clear'
            else:   
                self.t.loc[ix,self.filter_col]=f'{hdr["FILTER"].lower()}'
                if "PUPIL" in hdr:
                    self.t.loc[ix,self.pupil_col]=f'{hdr["PUPIL"].lower()}'
                else:
                    self.t.loc[ix,self.pupil_col]=None
        
        if self.verbose:
            print(f'##################\n### Found {len(ixs)} input files with the correct filepatterns {filepatterns}')
            
        # if specified, select on filters
        if filters is not None:
            ixs_filters = []
            for filt in filters:
                ixs_filters.extend(self.ix_equal(self.filter_col,filt.lower()))
            print(f'### after filters cut ({filters}): {len(ixs_filters)} input files left')
            self.t = self.t.loc[ixs_filters]
            
        # if specified, select on pupils
        if pupils is not None:
            ixs_pupils = []
            for pupil in pupils:
                ixs_pupils.extend(self.ix_equal(self.filter_col,pupil.lower()))
            print(f'### after pupils cut ({pupils}): {len(ixs_pupils)} input files left')
            self.t = self.t.loc[ixs_pupils]
            
        if self.verbose>2:
            print('### Input files:')
            self.write()

    
    def align_wcs(self, ixs, 
                  #outrootdir=None, 
                  #outsubdir=None,
                  overwrite = False,
                  skip_if_exists = False,
                  telescope = None,
                  #skip_applydistortions_if_exists = False,
                  # refcat parameters
                  refcatname = 'Gaia',
                  refcat_racol='auto',
                  refcat_deccol='auto',
                  refcat_magcol = None,
                  refcat_magerrcol = None,
                  refcat_colorcol = None,
                  pmflag = False,
                  pm_median=False,
                  load_photcat_if_exists=False,
                  rematch_refcat=False,
                  SNR_min = 10.0, # minimum S/N for photometry
                  # find best matches to refcut
                  d2d_max = None, # maximum distance refcat to source in image
                  dmag_max = 0.1, # maximum uncertainty of source 
                  sharpness_lim = (None, None), # sharpness limits
                  roundness1_lim = (None, None), # roundness1 limits 
                  delta_mag_lim = (None, None), # limits on mag-refcat_mainfilter
                  objmag_lim = (None,None), # limits on mag, the magnitude of the objects in the image
                  refmag_lim = (None,None), # limits on refcat_mainfilter, the magnitude of the reference catalog                    
                  slope_min=-10/2048.0, 
                  Nbright4match=None, # Use only the the brightest  Nbright sources from image for the matching with the ref catalog
                  Nbright=None,    # Use only the brightest Nbright sources from image
                  histocut_order='dxdy', # histocut_order defines whether the histogram cut is first done for dx or first for dy
                  xshift=0.0,# added to the x coordinate before calculating ra,dec. This can be used to correct for large shifts before matching!
                  yshift=0.0, # added to the y coordinate before calculating ra,dec. This can be used to correct for large shifts before matching!
                  showplots=0,
                  saveplots=0,
                  savephottable=0):
        
        self.t.loc[ixs,'errorflag'] = None
        self.t.loc[ixs,'skipflag'] = None
                
        #self.wcs_align.set_outdir(outrootdir,outsubdir)
        
        for ix in ixs:
            #distfile = self.t.loc[ix,'distortion_match']
            inputfile = self.t.loc[ix,'filename']

            self.wcs_align = st_wcs_align()
            self.wcs_align.calphot=jwst_photclass()
            self.wcs_align.outdir = self.outdir

            self.wcs_align.verbose = self.verbose
            self.wcs_align.replace_sip = self.replace_sip
            self.wcs_align.sip_err = self.sip_err
            self.wcs_align.sip_degree = self.sip_degree
            self.wcs_align.sip_points = self.sip_points
            
            self.wcs_align.rough_cut_px_min = self.rough_cut_px_min
            self.wcs_align.rough_cut_px_max = self.rough_cut_px_max
            self.wcs_align.d_rotated_Nsigma = self.d_rotated_Nsigma

            # If debugging: just run one, outside the try block so that we can get real error messages
            if self.debug:
                self.wcs_align.run_all(inputfile,
                                       #distortion_file=distfile,  
                                       telescope = telescope,
                                       overwrite = overwrite,
                                       skip_if_exists = skip_if_exists,
                                       #skip_applydistortions_if_exists=skip_applydistortions_if_exists,
                                       refcatname = refcatname,
                                       refcat_racol = refcat_racol,
                                       refcat_deccol = refcat_deccol,
                                       refcat_magcol=refcat_magcol,
                                       refcat_magerrcol=refcat_magerrcol,
                                       refcat_colorcol=refcat_colorcol,
                                       pmflag = pmflag,
                                       pm_median = pm_median,
                                       load_photcat_if_exists=load_photcat_if_exists,
                                       rematch_refcat=rematch_refcat,
                                       SNR_min = SNR_min,
                                       d2d_max = d2d_max, # maximum distance refcat to source in image
                                       dmag_max = dmag_max, # maximum uncertainty of source 
                                       sharpness_lim = sharpness_lim, # sharpness limits
                                       roundness1_lim = roundness1_lim, # roundness1 limits 
                                       delta_mag_lim =  delta_mag_lim, # limits on mag-refcat_mainfilter
                                       objmag_lim = objmag_lim, # limits on mag, the magnitude of the objects in the image
                                       refmag_lim = refmag_lim, # limits on refcat_mainfilter, the magnitude of the reference catalog
                                       slope_min = slope_min,
                                       Nbright4match=Nbright4match, # Use only the the brightest  Nbright sources from image for the matching with the ref catalog
                                       Nbright=Nbright,    # U/se only the brightest Nbright sources from image
                                       histocut_order=histocut_order, # histocut_order defines whether the histogram cut is first done for dx or first for dy
                                       xshift=xshift,# added to the x coordinate before calculating ra,dec. This can be used to correct for large shifts before matching!
                                       yshift=yshift, # added to the y coordinate before calculating ra,dec. This can be used to correct for large shifts before matching!
                                       showplots=showplots,
                                       saveplots=saveplots,# 
                                       savephottable=savephottable)
                self.t.loc[ix,'errorflag']=False
            else:
                try:
                    self.wcs_align.run_all(inputfile,
                                           #distortion_file=distfile,                     
                                           telescope = telescope,
                                           overwrite = overwrite,
                                           skip_if_exists = skip_if_exists,
                                           #skip_applydistortions_if_exists=skip_applydistortions_if_exists,
                                           refcatname = refcatname,
                                           refcat_racol = refcat_racol,
                                           refcat_deccol = refcat_deccol,
                                           refcat_magcol=refcat_magcol,
                                           refcat_magerrcol=refcat_magerrcol,
                                           refcat_colorcol=refcat_colorcol,
                                           pmflag = pmflag,
                                           pm_median = pm_median,
                                           load_photcat_if_exists=load_photcat_if_exists,
                                           rematch_refcat=rematch_refcat,
                                           SNR_min = SNR_min,
                                           d2d_max = d2d_max, # maximum distance refcat to source in image
                                           dmag_max = dmag_max, # maximum uncertainty of source 
                                           sharpness_lim = sharpness_lim, # sharpness limits
                                           roundness1_lim = roundness1_lim, # roundness1 limits 
                                           delta_mag_lim =  delta_mag_lim, # limits on mag-refcat_mainfilter
                                           objmag_lim = objmag_lim, # limits on mag, the magnitude of the objects in the image
                                           refmag_lim = refmag_lim, # limits on refcat_mainfilter, the magnitude of the reference catalog
                                           slope_min = slope_min,
                                           Nbright4match=Nbright4match, # Use only the the brightest  Nbright sources from image for the matching with the ref catalog
                                           Nbright=Nbright,    # U/se only the brightest Nbright sources from image
                                           histocut_order=histocut_order, # histocut_order defines whether the histogram cut is first done for dx or first for dy
                                           xshift=xshift,# added to the x coordinate before calculating ra,dec. This can be used to correct for large shifts before matching!
                                           yshift=yshift, # added to the y coordinate before calculating ra,dec. This can be used to correct for large shifts before matching!
                                           showplots=showplots,
                                           saveplots=saveplots,# 
                                           savephottable=savephottable)
                    self.t.loc[ix,'errorflag']=False
                except Exception as e:
                    print(f'ERROR while running {inputfile}: {e}')
                    self.t.loc[ix,'errorflag']=True

if __name__ == '__main__':

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


    # set the output directory
    align_batch.set_outdir(outrootdir=args.outrootdir,
                           outsubdir=args.outsubdir)

    # get the input files
    align_batch.get_input_files(args.input_files,directory=args.input_dir,
                                filters=args.filters,pupils=args.pupils)
    
    # apply new distortions?
    #if args.distortion_files is not None:
    #    align_batch.get_distortion_files(args.distortion_files,directory=None)
    #    ixs_matches,ixs_not_matches = align_batch.match_distortion4inputfile(apertures=args.apertures, 
    #                                                                 filts=args.filters, 
    #                                                                 pupils=args.pupils)
    #else:
        
    #align_batch.distortionfiles.t['filename']=None
    #align_batch.t['distortion_match']=None
    ixs_all = align_batch.getindices()
    
    #align_batch.get_output_filenames()

    if len(ixs_all)==0:
        print('NO IMAGES FOUND!! exiting...')
        sys.exit(0)
        
    # get the output filenames
    ixs_exists,ixs_notexists = align_batch.get_output_filenames(ixs=ixs_all)    
    
    
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


    print(f'Output directory:{align_batch.outdir}')
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
                        telescope = args.telescope,
                        #skip_applydistortions_if_exists=args.skip_applydistortions_if_exists,
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
                        showplots=args.showplots,
                        saveplots=args.saveplots,# 
                        savephottable=args.savephottable
                        )

    align_batch.write()