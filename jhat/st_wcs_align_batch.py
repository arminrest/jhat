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

        self.distortionfiles = pdastroclass()        

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
        if 'JWST_OUTPUT_ROOTDIR' in os.environ:
            outrootdir = os.environ['JWST_OUTPUT_ROOTDIR']
        else:
            outrootdir = None


        parser.add_argument('--input_dir', default=inputdir, help='Directory in which the input images are located. If $JWST_INPUT_IMAGEDIR is defined, then this dir is taken as default (default=%(default)s)')
        parser.add_argument('--input_files', nargs='+', default=['*_cal.fits'], help='list of cal or rate file(pattern)s to which the distortion files are applied to. "input_dir" is used if not None (default=%(default)s)')

        #parser.add_argument('--outrootdir', default=outrootdir, help='output root directory. The output directoy is the output root directory + the outsubdir if not None.  If $JWST_OUTPUT_ROOTDIR is defined, then this dir is taken as default (default=%(default)s)')
        #parser.add_argument('--outsubdir', default=None, help='outsubdir added to output root directory (default=%(default)s)')
        parser.add_argument('--addfilter2outsubdir', default=False, action='store_true', help='add the filter to the outsubdir')
        #parser.add_argument('--overwrite', default=False, action='store_true', help='overwrite files if they exist.')

        parser.add_argument('--skip_if_exists', default=False, action='store_true', help='Skip doing the analysis of a given input image if the cal file already exists, assuming the full analysis has been already done')

        parser.add_argument('-v','--verbose', default=0, action='count')

        parser.add_argument('--detectors', nargs='+', default=None, help='constrain the input file list to these detectors (default=%(default)s)')
        parser.add_argument('--filters', nargs='+', default=None, help='constrain the input file list to these filters (default=%(default)s)')
        parser.add_argument('--pupils', nargs='+', default=None, help='constrain the input file list to these pupils (default=%(default)s)')

        parser.add_argument('--distortioncoeffs_dir', default=None, help='Directory in which the distortion coefficients are. If directory is specified, then in this directory matching distortion files of the form <aperture>_<filter>_<pupil>.polycoeff.asdf (e.g., nrcb3_full_f187n_clear.polycoeff.asdf) are looked for and then applied before the WCS alignment.')
        
        parser.add_argument('-d','--debug', default=False, action='store_true',help='debug mode: alignment is done outside "try" block!')

        parser = self.wcs_align.default_options(parser)

        return(parser)

    def addfilter2outsubdir(self,outsubdir,addfilter2outsubdir,ix,filt=None,pupil=None):
        outsubdir_filter = outsubdir
        if addfilter2outsubdir:
            if filt is None: filt = self.t.loc[ix,"filter"].lower()
            outsubdir_filter+=f'/{filt}'
            if (pupil is not None) or  (isinstance(pupil,str) and self.t.loc[ix,"pupil"].lower()!='clear'):
                if pupil is None: pupil = self.t.loc[ix,"pupil"].lower()
                outsubdir_filter+=f'_{pupil}'
        return(outsubdir_filter)

    
    #def get_output_filenames(self, suffixmapping = {'cal':'tweakregstep','rate':'tweakregstep'},ixs=None):
    def get_output_filenames(self, outputsuffix = 'jhat',ixs=None,outrootdir=None,outsubdir=None,addfilter2outsubdir=False):
        ixs = self.getindices(ixs)
        ixs_exists=[]
        ixs_notexists=[]
        for ix in ixs:
            #pattern = '(.*)_([a-zA-Z0-9]+)\.fits$'
            #m = re.search(pattern,self.t.loc[ix,'filename'])
            #if m is None:
            #    raise RuntimeError(f'cannot determine suffix for file {self.t.loc[ix,"filename"]} with pattern {pattern}!')
            #(inputname_nosuffix,inputsuffix)=m.groups()
            #inputbasename = os.path.basename(inputname_nosuffix)
            
            outsubdir_filter = self.addfilter2outsubdir(outsubdir,addfilter2outsubdir,ix)
            
            outfilebasename = self.wcs_align.set_outbasename(outrootdir=outrootdir,
                                                         outsubdir=outsubdir_filter,
                                                         inputname=self.t.loc[ix,'filename'])
            
            outfilename = f'{outfilebasename}_{outputsuffix}.fits'
            
            #outfilename = f'{self.outdir}/{inputbasename}_{suffixmapping[inputsuffix]}.fits'
            #outfilename = f'{self.outdir}/{inputbasename}_{outputsuffix}.fits'
            self.t.loc[ix,'outfilename']=outfilename
            if os.path.isfile(outfilename):
                self.t.loc[ix,'outfile_exists']='yes'
                ixs_exists.append(ix)
            else:
                self.t.loc[ix,'outfile_exists']='no'
                ixs_notexists.append(ix)
        
        if self.verbose>1:  self.write(indices=ixs)

        return(ixs_exists,ixs_notexists)


    #def set_outdir(self,outrootdir=None,outsubdir=None):
    #    self.outdir = self.wcs_align.set_outdir(outrootdir,outsubdir)
    #    return(self.outdir)

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
    
    
    def get_input_files(self,filepatterns,directory=None,detectors=None,filters=None,pupils=None):
        self.t['filename'] = self.get_files(filepatterns,directory=directory)
        ixs = self.getindices()
        for ix in ixs:
            print(self.t.loc[ix,'filename'])
            #print(fits.getheader(self.t.loc[ix,'filename']))
            hdr = fits.getheader(self.t.loc[ix,'filename'])
            #print('CCCC',hdr['TELESCOP'],hdr['DETECTOR'],hdr["INSTRUME"],hdr["FILTER"],hdr["PUPIL"])
            detector = re.sub('long$','5',hdr['DETECTOR'].lower())
            if hdr['TELESCOP'].lower()=='jwst':
                self.t.loc[ix,self.aperture_col]=f'{detector}_{hdr["SUBARRAY"].lower()}'
                self.t.loc[ix,'subarray']=f'{hdr["SUBARRAY"].lower()}'
            elif hdr['TELESCOP'].lower()=='hst':
                self.t.loc[ix,self.aperture_col]=f'{detector}'
            else:
                raise RuntimeError('Cannot identify the telescope {hdr["TELESCOP"]}! Must be hst or jwst!')
            self.t.loc[ix,'telescope']=f'{hdr["TELESCOP"].lower()}'
            self.t.loc[ix,'detector']=f'{detector}'
            self.t.loc[ix,'instrument']=f'{hdr["INSTRUME"].lower()}'
            if self.t.loc[ix,'instrument']=='fgs':
                self.t.loc[ix,self.filter_col]='clear'
                self.t.loc[ix,self.pupil_col]='clear'
            elif self.t.loc[ix,'instrument']=='acs':
                f1=hdr["FILTER1"].lower()
                f2=hdr["FILTER2"].lower()
                if re.search('clear',f1) is not None:
                    self.t.loc[ix,self.filter_col]=f'{f2}'
                else:
                    self.t.loc[ix,self.filter_col]=f'{f1}'
            else:   
                self.t.loc[ix,self.filter_col]=f'{hdr["FILTER"].lower()}'
                if "PUPIL" in hdr:
                    self.t.loc[ix,self.pupil_col]=f'{hdr["PUPIL"].lower()}'
                else:
                    self.t.loc[ix,self.pupil_col]=None
        
        if self.verbose:
            print(f'##################\n### Found {len(ixs)} input files with the correct filepatterns {filepatterns}')
        if len(ixs)<1:
            return(1)
        
        # if specified, select on detectors
        if detectors is not None:
            ixs_detectors = []
            for detector in detectors:
                if re.search('^a[1-5]$|^b[1-5]$', detector.lower()):
                    detector = f'nrc{detector.lower()}'
                ixs_detectors.extend(self.ix_equal(self.detector_col,detector.lower()))
            print(f'### after detectors cut ({detectors}): {len(ixs_detectors)} input files left')
            self.t = self.t.loc[ixs_detectors]
            
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
                ixs_pupils.extend(self.ix_equal(self.pupil_col,pupil.lower()))
            print(f'### after pupils cut ({pupils}): {len(ixs_pupils)} input files left')
            self.t = self.t.loc[ixs_pupils]
            
        if self.verbose>2:
            print('### Input files:')
            self.write()

    def get_distortioncoeff_files(self,coeffdir,filepatterns=['*.polycoeff.asdf']):
        self.distortionfiles.t['filename'] = self.get_files(filepatterns,directory=coeffdir)
        for ix in self.distortionfiles.getindices():
            
            # get the filter and save it in the 'filter' column
            #m = re.search('^([a-zA-Z0-9]+_[a-zA-Z0-9]+)_(f[a-zA-Z0-9]+)_([a-zA-Z0-9]+)',os.path.basename(self.distortionfiles.t.loc[ix,'filename']))
            #if m is None:
            #    raise RuntimeError(f'could not parse filename {os.path.basename(self.distortionfiles.t.loc[ix,"filename"])} for aperture, filter and/or pupil!')
            #aperture,filt,pupil = m.groups()

            # nrcb3_full_f187n_clear.polycoeff.asdf

            m0 = re.search('^([a-zA-Z0-9]+_[a-zA-Z0-9]+)\.polycoeff\.asdf',os.path.basename(self.distortionfiles.t.loc[ix,'filename']))
            if m0 is not None:
                aperture, = m0.groups()    
                filt = pupil = 'clear'
            else:
                m1 = re.search('^([a-zA-Z0-9]+_[a-zA-Z0-9]+).*_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)\.polycoeff\.asdf',os.path.basename(self.distortionfiles.t.loc[ix,'filename']))
                if m1 is not None:
                    aperture,filt,pupil = m1.groups()        
                else:
                    raise RuntimeError(f'could not parse filename {os.path.basename(self.distortionfiles.t.loc[ix,"filename"])} for aperture, filter and/or pupil!')

            m3 = re.search('(^[a-zA-Z0-9]+)',aperture)    
            if m3 is not None:
                detector =  m3.groups()[0]
            else:
                raise RuntimeError(f'could not parse aperture {aperture} for detector!')
                

            self.distortionfiles.t.loc[ix,self.aperture_col]=f'{aperture}'
            self.distortionfiles.t.loc[ix,self.detector_col]=f'{detector}'
            self.distortionfiles.t.loc[ix,self.filter_col]=f'{filt}'
            self.distortionfiles.t.loc[ix,self.pupil_col]=f'{pupil}'
        if self.verbose:
            print('##################\n### Distortion files:')
            ixs = self.distortionfiles.ix_sort_by_cols([self.filter_col,self.pupil_col,self.aperture_col])
            self.distortionfiles.write(indices=ixs)

    def match_distortioncoeffs(self, ixs):
        for ix in ixs:
            ixs_distcoeff = self.distortionfiles.getindices()
            for col in [self.filter_col, self.pupil_col, self.aperture_col]:
                ixs_distcoeff = self.distortionfiles.ix_equal(col, self.t.loc[ix,col], indices=ixs_distcoeff)
            if len(ixs_distcoeff)==0:
                raise RuntimeError(f'no distortion coefficient file found for {self.t.loc[ix,self.aperture_col]}, {self.t.loc[ix,self.filter_col]}, {self.t.loc[ix,self.pupil_col]}')
            elif len(ixs_distcoeff)>1:
                self.distortionfiles.write(indices=ixs_distcoeff)
                raise RuntimeError(f'more than one distortion coefficient files found for {self.t.loc[ix,self.aperture_col]}, {self.t.loc[ix,self.filter_col]}, {self.t.loc[ix,self.pupil_col]}')
            else:
                self.t.loc[ix,'distcoefffile']=self.distortionfiles.t.loc[ixs_distcoeff[0],'filename']
        self.write()
        sys.exit(0)
        return(0)
        
    
    def align_wcs(self, ixs, 
                  outrootdir=None, 
                  outsubdir=None,
                  addfilter2outsubdir=False,
                  overwrite = False,
                  skip_if_exists = False,
                  telescope = None,
                  #skip_applydistortions_if_exists = False,
                  photometry_method = 'aperture',
                  find_stars_threshold = 3.0,
                  sci_xy_catalog=None,
                  # refcat parameters
                  use_dq=False,
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
                  iterate_with_xyshifts = False, # After the first histogram fit, redo the match with refcat with x/yshift=median(dx/dy) and redo histofit. Use this if the offsets are big, since the second iteration will give you better matching with the refcat
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
            #self.wcs_align.outdir = self.outdir

            self.wcs_align.verbose = self.verbose
            self.wcs_align.replace_sip = self.replace_sip
            self.wcs_align.sip_err = self.sip_err
            self.wcs_align.sip_degree = self.sip_degree
            self.wcs_align.sip_points = self.sip_points
            
            self.wcs_align.rough_cut_px_min = self.rough_cut_px_min
            self.wcs_align.rough_cut_px_max = self.rough_cut_px_max
            self.wcs_align.d_rotated_Nsigma = self.d_rotated_Nsigma

            outsubdir_filter = self.addfilter2outsubdir(outsubdir,addfilter2outsubdir,ix)            
            
            if telescope is None:
                if "telescope" in self.t.columns:
                    telescope4image = self.t.loc[ix,'telescope'].upper()
                else:
                    telescope4image = None
            else:
                telescope4image = telescope    
            print(f'Telescope {telescope4image} for image {self.t.loc[ix,"filename"]}')

            # If debugging: just run one, outside the try block so that we can get real error messages
            if self.debug:
                self.wcs_align.run_all(inputfile,
                                       #distortion_file=distfile,  
                                       telescope = telescope4image,
                                       outrootdir = outrootdir,
                                       outsubdir = outsubdir_filter,
                                       overwrite = overwrite,
                                       skip_if_exists = skip_if_exists,
                                       #skip_applydistortions_if_exists=skip_applydistortions_if_exists,
                                       photometry_method = photometry_method,
                                       find_stars_threshold = find_stars_threshold,
                                       sci_xy_catalog=sci_xy_catalog,
                                       use_dq = use_dq,
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
                                       iterate_with_xyshifts=iterate_with_xyshifts, # After the first histogram fit, redo the match with refcat with x/yshift=median(dx/dy) and redo histofit. Use this if the offsets are big, since the second iteration will give you better matching with the refcat
                                       showplots=showplots,
                                       saveplots=saveplots,# 
                                       savephottable=savephottable)
                self.t.loc[ix,'errorflag']=False
            else:
                try:
                    self.wcs_align.run_all(inputfile,
                                           #distortion_file=distfile,                     
                                           telescope = telescope4image,
                                           outrootdir = outrootdir,
                                           outsubdir = outsubdir_filter,
                                           overwrite = overwrite,
                                           skip_if_exists = skip_if_exists,
                                           #skip_applydistortions_if_exists=skip_applydistortions_if_exists,
                                           photometry_method = photometry_method,
                                           find_stars_threshold = find_stars_threshold,
                                           sci_xy_catalog=sci_xy_catalog,
                                           use_dq = use_dq,
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
                                           iterate_with_xyshifts=iterate_with_xyshifts, # After the first histogram fit, redo the match with refcat with x/yshift=median(dx/dy) and redo histofit. Use this if the offsets are big, since the second iteration will give you better matching with the refcat
                                           showplots=showplots,
                                           saveplots=saveplots,# 
                                           savephottable=savephottable)
                    self.t.loc[ix,'errorflag']=False
                except Exception as e:
                    print(f'ERROR while running {inputfile}: {e}')
                    self.t.loc[ix,'errorflag']=True

