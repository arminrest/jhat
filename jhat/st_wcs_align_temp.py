#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 14:32:42 2022

@author: arest, bhilbert, mcorrenti, acanipe, jpierel
"""

import os,re,sys,copy,warnings
#from jwst.tweakreg import TweakRegStep
import tweakreg_hack
import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

from jwst import datamodels

from .simple_jwst_phot import jwst_photclass,hst_photclass
from .pdastro import *

warnings.simplefilter('ignore')
__all__ = ['st_wcs_align']

plot_style={}
plot_style['good_data']={'style':'o','color':'blue', 'ms':5 ,'alpha':0.5}
plot_style['cut_data']={'style':'o','color':'red', 'ms':5 ,'alpha':0.3}
plot_style['do_not_use_data']={'style':'o','color':'gray', 'ms':3 ,'alpha':0.3}

def initplot(nrows=1, ncols=1, figsize4subplot=5, **kwargs):
    sp=[]
    xfigsize=figsize4subplot*ncols
    yfigsize=figsize4subplot*nrows
    plt.figure(figsize=(xfigsize,yfigsize))
    counter=1
    for row in range(nrows):
        for col in range(ncols):
            sp.append(plt.subplot(nrows, ncols, counter,**kwargs))
            counter+=1

    for i in range(len(sp)):
        plt.setp(sp[i].get_xticklabels(),'fontsize',12)
        plt.setp(sp[i].get_yticklabels(),'fontsize',12)
        sp[i].set_xlabel(sp[i].get_xlabel(),fontsize=14)
        sp[i].set_ylabel(sp[i].get_ylabel(),fontsize=14)
        sp[i].set_title(sp[i].get_title(),fontsize=14)

    return(sp)


# These are the info plots to check out dx, dy for good and bad detections
def infoplots(phot,ixs_good,dy_plotlim=(-4,4),dx_plotlim=(-4,4)):
    sp=initplot(2,3)

    if phot.ixs_notuse is not None:
        phot.t.loc[phot.ixs_notuse].plot('y','dx',ax=sp[0],ylim=dx_plotlim, **plot_style['do_not_use_data'])
        phot.t.loc[phot.ixs_notuse].plot('y','dx',ax=sp[0],ylim=dx_plotlim, **plot_style['do_not_use_data'])
        phot.t.loc[phot.ixs_notuse].plot('x','dy',ax=sp[1],ylim=dy_plotlim, **plot_style['do_not_use_data'])
        phot.t.loc[phot.ixs_notuse].plot('x','y',ax=sp[2],**plot_style['do_not_use_data'])
        phot.t.loc[phot.ixs_notuse].plot('sharpness','mag',ax=sp[3],**plot_style['do_not_use_data'])
        phot.t.loc[phot.ixs_notuse].plot('sharpness','dmag',ax=sp[4],**plot_style['do_not_use_data'])
        phot.t.loc[phot.ixs_notuse].plot('sharpness','roundness1',ax=sp[5],**plot_style['do_not_use_data'])

    if phot.ixs_use is not None:
        ixs_cut = AnotB(phot.ixs_use,ixs_good)
        phot.t.loc[ixs_cut].plot('y','dx',ax=sp[0],ylim=dx_plotlim, **plot_style['cut_data'])
        phot.t.loc[ixs_cut].plot('y','dx',ax=sp[0],ylim=dx_plotlim, **plot_style['cut_data'])
        phot.t.loc[ixs_cut].plot('x','dy',ax=sp[1],ylim=dy_plotlim, **plot_style['cut_data'])
        phot.t.loc[ixs_cut].plot('x','y',ax=sp[2],**plot_style['cut_data'])
        phot.t.loc[ixs_cut].plot('sharpness','mag',ax=sp[3],**plot_style['cut_data'])
        phot.t.loc[ixs_cut].plot('sharpness','dmag',ax=sp[4],**plot_style['cut_data'])
        phot.t.loc[ixs_cut].plot('sharpness','roundness1',ax=sp[5],**plot_style['cut_data'])


    
    phot.t.loc[ixs_good].plot('y','dx',ax=sp[0],ylim=dx_plotlim, ylabel='dx in pixels', **plot_style['good_data'])
    phot.t.loc[ixs_good].plot('x','dy',ax=sp[1],ylim=dy_plotlim, ylabel='dy in pixels', **plot_style['good_data'])
    phot.t.loc[ixs_good].plot('x','y',ax=sp[2],**plot_style['good_data'],ylabel='y')
    phot.t.loc[ixs_good].plot('sharpness','mag',ax=sp[3],**plot_style['good_data'],ylabel='mag')
    phot.t.loc[ixs_good].plot('sharpness','dmag',ax=sp[4],**plot_style['good_data'],ylabel='dmag')
    phot.t.loc[ixs_good].plot('sharpness','roundness1',ax=sp[5],**plot_style['good_data'],ylabel='roundness1')
    
    sp[4].set_ylim(0.0,0.15)

    for i in range(6): sp[i].get_legend().remove()
    
    plt.tight_layout()
    return(sp)

# plot the rotated dx or dy versus the original one
def plot_rotated(phot,ixs,d_col,col,
                 histotable,
                 apply_rolling_gaussian=False,
                 d_col_rot='d_rot_tmp',
                 sp=None,
                 spi=[0,1],
                 histolim=[-20,20],
                 bins=None,
                 bin_weights_flag=False,
                 title=None):
    if sp == None:
        sp=initplot(1,2)
        spi=[0,1]
        
    # Get the limits for the plots
    # the maximum limits are passed with histolim
    # however, if the all data is within the limits, then
    # shrink the limits accordingly
    minval = np.min(phot.t.loc[ixs,d_col_rot])
    maxval = np.max(phot.t.loc[ixs,d_col_rot])
    if minval>histolim[0]:histolim[0]=minval
    if maxval<histolim[1]:histolim[1]=maxval
    # add a buffer to the min and max values
    histolim[0]-=0.1*(maxval-minval)
    histolim[1]+=0.1*(maxval-minval)
        
    if phot.ixs_notuse is not None:
        phot.t.loc[phot.ixs_notuse].plot(col,d_col_rot,ax=sp[spi[0]],**plot_style['do_not_use_data'])
    phot.t.loc[ixs].plot(col,d_col,ax=sp[spi[0]],ylim=histolim,title=title,**plot_style['cut_data'])
    phot.t.loc[ixs].plot(col,d_col_rot,ax=sp[spi[0]],ylim=histolim,ylabel=f'{d_col} [pixels]',**plot_style['good_data'])
    sp[spi[0]].get_legend().remove()

    
    if bins is not None:
        if bin_weights_flag:
            phot.t.loc[ixs,d_col_rot].plot.hist(ax=sp[spi[1]],bins=bins,
                                                weights=phot.t.loc[ixs,'__weights'],
                                                xlim=histolim,color='blue',histtype='step')
        else:
            phot.t.loc[ixs,d_col_rot].plot.hist(ax=sp[spi[1]],bins=bins,xlim=histolim,color='blue',histtype='step')
        
        if apply_rolling_gaussian:
            histotable.t.plot('bincenter','histosum',ax=sp[spi[1]],color='red')

        sp[spi[1]].set_xlabel(f'rotated {d_col}')
        #sp[spi[1]].get_legend().remove()
    return(sp)


def initial_dxdy_plot(phot, ixs_use, ixs_notuse, 
                      plots_dxdy_delta_pix_ylim=20,
                      refcat_mainfilter=None,refcat_mainfilter_err=None,refcat_maincolor=None,
                      d2d_max=None,dmag_max=None,Nbright=None,delta_mag_lim=None,verbose=0):
    
    sp = initplot(2,3)
    
    #print('FFFFFFF',len(ixs_notuse),len(ixs_use),len(phot.t))
    #sys.exit(0)
    
    dx_median = phot.t.loc[ixs_use,'dx'].median()
    dy_median = phot.t.loc[ixs_use,'dy'].median()
    if verbose:
        print(f'dx median: {dx_median}\ndy median: {dy_median}')
 
    # these are the general limits for the y-axis for the dx/dy plots
    dy_plotlim = (dy_median-plots_dxdy_delta_pix_ylim,dy_median+plots_dxdy_delta_pix_ylim)
    dx_plotlim = (dx_median-plots_dxdy_delta_pix_ylim,dx_median+plots_dxdy_delta_pix_ylim)


    # plot the residuals
    title = f'Initial cut: d2d_max={d2d_max},\ndmag_max={dmag_max}'
    title_Nbright = f'Nbright={Nbright}'
    title_deltamag = f'delta_mag_lim={delta_mag_lim}'
    phot.t.loc[ixs_notuse].plot('y','dx',ax=sp[0],ylim=dx_plotlim,title=title,**plot_style['do_not_use_data'])
    phot.t.loc[ixs_use].plot('y','dx',ax=sp[0],ylim=dx_plotlim, ylabel='dx [pixel]',**plot_style['good_data'])
    phot.t.loc[ixs_notuse].plot('x','dy',ax=sp[1],ylim=dx_plotlim,**plot_style['do_not_use_data'])
    phot.t.loc[ixs_use].plot('x','dy',ax=sp[1],ylim=dy_plotlim,ylabel='dy [pixel]',**plot_style['good_data'])
    phot.t.loc[ixs_notuse].plot('x','y',ax=sp[2],**plot_style['do_not_use_data'])
    phot.t.loc[ixs_use].plot('x','y',ax=sp[2],ylabel='y [pixel]',**plot_style['good_data'])
    phot.t.loc[ixs_notuse].plot('sharpness','mag',ax=sp[3],**plot_style['do_not_use_data'])
    phot.t.loc[ixs_use].plot('sharpness','mag',ax=sp[3],title=title_Nbright,ylabel='mag',**plot_style['good_data'])
    if phot.refcat_mainfilter is not None:
        if phot.refcat_maincolor is not None:
            phot.t.loc[ixs_notuse].plot(phot.refcat_maincolor,'delta_mag',ax=sp[4],**plot_style['do_not_use_data'])
            phot.t.loc[ixs_use].plot(phot.refcat_maincolor,'delta_mag',title=title_deltamag,ax=sp[4],ylabel=f'mag - {phot.refcat_mainfilter}',**plot_style['good_data'])
            phot.t.loc[ixs_notuse].plot(phot.refcat_maincolor,phot.refcat_mainfilter,ax=sp[5],**plot_style['do_not_use_data'])
            phot.t.loc[ixs_use].plot(phot.refcat_maincolor,phot.refcat_mainfilter,ax=sp[5],ylabel=f'{phot.refcat_mainfilter}',**plot_style['good_data'])
            for i in range(6): sp[i].get_legend().remove()
        else:
            phot.t.loc[ixs_notuse].plot(phot.refcat_mainfilter,'delta_mag',title=title_deltamag,ax=sp[4],**plot_style['do_not_use_data'])
            phot.t.loc[ixs_use].plot(phot.refcat_mainfilter,'delta_mag',ax=sp[4],ylabel=f'mag - {phot.refcat_mainfilter}',**plot_style['good_data'])
            for i in range(5): sp[i].get_legend().remove()
    else:
        for i in range(4): sp[i].get_legend().remove()

    plt.tight_layout()
    print('*** Note: close plot to continue!')
    plt.show()
    
    return(sp)

def dxdy_plot(phot,ixs_selected, sp=None, spi = [0,1,4,5,8,9,2,6,10,3,7,11], title=None,
              refcat_mainfilter=None,refcat_mainfilter_err=None,refcat_maincolor=None):
    if sp is None: sp=initplot(3,4)

    dx_median = phot.t.loc[ixs_selected,'dx'].median()
    dy_median = phot.t.loc[ixs_selected,'dy'].median()

    dlim_big = 15
    #dlim_small= 1
    dx_ylim_big = (dx_median-dlim_big,dx_median+dlim_big)
    dy_ylim_big = (dy_median-dlim_big,dy_median+dlim_big)
    #dx_ylim_small = (dx_median-dlim_small,dx_median+dlim_small)
    #dy_ylim_small = (dy_median-dlim_small,dy_median+dlim_small)

    if (refcat_mainfilter is not None):
        if (refcat_mainfilter in phot.t.columns):
            phot.t['delta_mag'] = phot.t['mag'] - phot.t[refcat_mainfilter]
        else:
            print(f'Warning! cannot make plot with reference filter {refcat_mainfilter}, since it is not a column in the table')
    
    
    if phot.ixs_notuse is not None:
        phot.t.loc[phot.ixs_notuse].plot('y','dx',ax=sp[spi[0]], **plot_style['do_not_use_data'])
        phot.t.loc[phot.ixs_notuse].plot('x','dy',ax=sp[spi[1]], **plot_style['do_not_use_data'])
        phot.t.loc[phot.ixs_notuse].plot('y','dx',ax=sp[spi[2]], **plot_style['do_not_use_data'])
        phot.t.loc[phot.ixs_notuse].plot('x','dx',ax=sp[spi[3]], **plot_style['do_not_use_data'])
        phot.t.loc[phot.ixs_notuse].plot('x','dy',ax=sp[spi[4]], **plot_style['do_not_use_data'])
        phot.t.loc[phot.ixs_notuse].plot('y','dy',ax=sp[spi[5]], **plot_style['do_not_use_data'])
        phot.t.loc[phot.ixs_notuse].plot('x','y' ,ax=sp[spi[6]], **plot_style['do_not_use_data'])
        phot.t.loc[phot.ixs_notuse].plot('sharpness','mag',ax=sp[spi[7]], **plot_style['do_not_use_data'])
        phot.t.loc[phot.ixs_notuse].plot('sharpness','roundness1',ax=sp[spi[8]], **plot_style['do_not_use_data'])
        phot.t.loc[phot.ixs_notuse].plot(refcat_maincolor,'delta_mag',ax=sp[spi[9]], **plot_style['do_not_use_data'])
        phot.t.loc[phot.ixs_notuse].plot(refcat_maincolor,refcat_mainfilter,ax=sp[spi[10]], **plot_style['do_not_use_data'])
        phot.t.loc[phot.ixs_notuse].plot(refcat_mainfilter,'mag',ax=sp[spi[11]], **plot_style['do_not_use_data'])

    if phot.ixs_use is not None:
        ixs_cut = AnotB(phot.ixs_use,ixs_selected)
        phot.t.loc[ixs_cut].plot('y','dx',ax=sp[spi[0]], **plot_style['cut_data'])
        phot.t.loc[ixs_cut].plot('x','dy',ax=sp[spi[1]], **plot_style['cut_data'])
        phot.t.loc[ixs_cut].plot('y','dx',ax=sp[spi[2]], **plot_style['cut_data'])
        phot.t.loc[ixs_cut].plot('x','dx',ax=sp[spi[3]], **plot_style['cut_data'])
        phot.t.loc[ixs_cut].plot('x','dy',ax=sp[spi[4]], **plot_style['cut_data'])
        phot.t.loc[ixs_cut].plot('y','dy',ax=sp[spi[5]], **plot_style['cut_data'])
        phot.t.loc[ixs_cut].plot('x','y' ,ax=sp[spi[6]], **plot_style['cut_data'])
        phot.t.loc[ixs_cut].plot('sharpness','mag',ax=sp[spi[7]], **plot_style['cut_data'])
        phot.t.loc[ixs_cut].plot('sharpness','roundness1',ax=sp[spi[8]], **plot_style['cut_data'])
        phot.t.loc[ixs_cut].plot(refcat_maincolor,'delta_mag',ax=sp[spi[9]], **plot_style['cut_data'])
        phot.t.loc[ixs_cut].plot(refcat_maincolor,refcat_mainfilter,ax=sp[spi[10]], **plot_style['cut_data'])
        phot.t.loc[ixs_cut].plot(refcat_mainfilter,'mag',ax=sp[spi[11]], **plot_style['cut_data'])
        
    phot.t.loc[ixs_selected].plot('y','dx',ylim=dx_ylim_big,ax=sp[spi[0]],ylabel='dx in pixels',title=title, **plot_style['good_data'])
    phot.t.loc[ixs_selected].plot('x','dy',ylim=dy_ylim_big,ax=sp[spi[1]],ylabel='dy in pixels',title=title, **plot_style['good_data'])

    (dx_min,dx_max) = (phot.t.loc[ixs_selected,'dx'].min(),phot.t.loc[ixs_selected,'dx'].max())
    dx_ylim_small = (dx_min - 1.0*(dx_max - dx_min), dx_max + 1.0*(dx_max - dx_min))
    phot.t.loc[ixs_selected].plot('y','dx',ylim=dx_ylim_small,ax=sp[spi[2]],ylabel='dx in pixels',title=title, **plot_style['good_data'])
    phot.t.loc[ixs_selected].plot('x','dx',ylim=dx_ylim_small,ax=sp[spi[3]],ylabel='dx in pixels',title=title, **plot_style['good_data'])

    (dy_min,dy_max) = (phot.t.loc[ixs_selected,'dy'].min(),phot.t.loc[ixs_selected,'dy'].max())
    dy_ylim_small = (dy_min - 1.0*(dy_max - dy_min), dy_max + 1.0*(dy_max - dy_min))
    phot.t.loc[ixs_selected].plot('x','dy',ylim=dy_ylim_small,ax=sp[spi[4]],ylabel='dy in pixels',title=title, **plot_style['good_data'])
    phot.t.loc[ixs_selected].plot('y','dy',ylim=dy_ylim_small,ax=sp[spi[5]],ylabel='dy in pixels',title=title, **plot_style['good_data'])

    phot.t.loc[ixs_selected].plot('x','y',ax=sp[spi[6]],ylabel='y',title=title, **plot_style['good_data'])
    phot.t.loc[ixs_selected].plot('sharpness','mag',ax=sp[spi[7]],ylabel='mag',title=title, **plot_style['good_data'])
    phot.t.loc[ixs_selected].plot('sharpness','roundness1',ax=sp[spi[8]],ylabel='roundness1',title=title, **plot_style['good_data'])
    phot.t.loc[ixs_selected].plot(refcat_maincolor,'delta_mag',ax=sp[spi[9]],ylabel=f'mag - {refcat_mainfilter}',title=title, **plot_style['good_data'])
    phot.t.loc[ixs_selected].plot(refcat_maincolor,refcat_mainfilter,ax=sp[spi[10]],ylabel=f'{refcat_mainfilter}',title=title, **plot_style['good_data'])
    phot.t.loc[ixs_selected].plot(refcat_mainfilter,'mag',ax=sp[spi[11]],ylabel='mag',title=title, **plot_style['good_data'])

    for i in range(12): 
        if sp[spi[i]].get_legend() is not None:
            sp[spi[i]].get_legend().remove()

    plt.tight_layout()
    return(sp)

def find_info_for_maxval(histotable,xcol='bincenter',ycol='histo',use_firstindex_if_multiple=True,):
    # get the max value of the histogram, and its associated bin
    #print(histo)
    xvals = np.array(histotable.t[xcol])
    yvals = np.array(histotable.t[ycol])
    N=len(histotable.t)
    
    maxval = np.max(yvals)
    ixs_max = np.where(yvals==maxval)
    multiple_max = None
    if (len(ixs_max[0])==0):
        raise RuntimeError('BUUUUGGGG!!!!')
    elif (len(ixs_max[0])>1):
        if not use_firstindex_if_multiple:
            raise RuntimeError(f'More than one indices {ixs_max} for maxval={maxval}')
        #print(f'\nWARNING!! More than one indices {ixs_max} for maxval={maxval}!')
        multiple_max=True
        ix_best=ixs_max[0][0]
    else:
        multiple_max=False
        ix_best=ixs_max[0][0]
    # return ()
    
    # get the (rough) fwhm of the peak
    ix_minus = ix_best-1 
    if ix_minus<0: ix_minus=0
    ix_plus = ix_best+1
    if ix_plus>N-1: ix_plus=N-1
    while (ix_minus>0):
        if yvals[ix_minus]<=0.5*yvals[ix_best]:
            break
        ix_minus-=1
    while (ix_plus<N-1):
        if yvals[ix_plus]<=0.5*yvals[ix_best]:
            break
        ix_plus+=1
    if ix_plus==ix_minus:
        print('Warning: problems getting FWHM!!! Setting it to binsize')
        ix_minus=0
        ix_plus=1
    fwhm = xvals[ix_plus]-xvals[ix_minus]
    
    return(xvals[ix_best],yvals[ix_best],ix_best,fwhm,multiple_max)

# straight line.
def f(val,slope,intercept):
    return(val*slope+intercept)

def find_binmax_for_slope(slope,phot,ixs,d_col,col,
                          Naxis_px,
                          gaussian_sigma,
                          windowsize,
                          halfwindowsize,
                          rot_results=None,
                          apply_rolling_gaussian=True,
                          d_col_rot='d_rot_tmp',
                          binsize=0.02,
                          bin_weights_flag=True,
                          showplots=0,
                          sp=None,
                          spi=[0,1]):

    intercept = -0.5*Naxis_px * slope
        
    #phot.t.loc[ixs,d_col_rot] = phot.t.loc[ixs,d_col] - f(phot.t.loc[ixs,col],slope,intercept)
    phot.t[d_col_rot] = phot.t[d_col] - f(phot.t[col],slope,intercept)

    # get the histogram
    d_rotated = phot.t.loc[ixs,d_col_rot]
    bins = np.arange(np.min(d_rotated),np.max(d_rotated)+binsize,binsize)
    if bin_weights_flag:
        histo = np.histogram(d_rotated,bins=bins,weights=phot.t.loc[ixs,'__weights'])
    else:
        histo = np.histogram(d_rotated,bins=bins)

    histotable=pdastroclass()            
    if apply_rolling_gaussian:
        # extend teh bins by +- halfwindowsize, since the kernel has the size windowsize. If this 
        # extension is not done, then the values at the edges within halfwindowsize are ignored!
        binmin = histo[1][0]+0.5*binsize-halfwindowsize*binsize
        bins_extended = np.arange(binmin,binmin+(len(histo[0])+2*halfwindowsize)*binsize,binsize)
        histotable.t['bincenter']=bins_extended

        # fill in the histogram values: make sure the get added +halfwindowsize into the table, 
        # since the bins got extended! The values outside are set to 0.0
        dataindices = range(halfwindowsize,len(histo[0])+halfwindowsize)
        histotable.t['histo']=0.0
        histotable.t.loc[dataindices,'histo']=histo[0]

        # TEST1 and TEST2 should be the same array (within floating point accuracy)
        #print('TEST1',histo[1][:-1]+0.5*binsize)
        #print('TEST2',histotable.t.loc[dataindices,'bincenter'])
        #print('TEST length',len(histo[1][:-1]+0.5*binsize),len(histo[0]),len(dataindices))

        histotable.t['histosum'] = histotable.t['histo'].rolling(windowsize,center=True,win_type='gaussian').sum(std=gaussian_sigma)
        ix_nan = histotable.ix_is_null('histosum')
        histotable.t.loc[ix_nan,'histosum']=0.0
        #sp=initplot(1,1)
        #histotable.t.plot('bincenter','histo',ax=sp[0],color='blue')
        #histotable.t.plot('bincenter','histosum',ax=sp[0],color='red')
        #histotable.write()

        (bincenter4maxval,maxval,index_maxval,fwhm,multiple_max) = find_info_for_maxval(histotable,ycol='histosum')

    else:
        # Note that the bincenter is the value in the bins (left edge of the bin) + the half of the binsize
        histotable.t['bincenter']=histo[1][:-1]+0.5*binsize
        histotable.t['histo']=histo[0]

        (bincenter4maxval,maxval,index_maxval,fwhm,multiple_max) = find_info_for_maxval(histotable)

    # Save the results
    if rot_results is not None:
        rot_results.newrow({'slope':slope,
                            'intercept':intercept,
                            'maxval':maxval,
                            'index':index_maxval,
                            'd_bestguess':bincenter4maxval,
                            'fwhm':fwhm,
                            'multimax':multiple_max
                            })
    # plot it if wanted
    if showplots>2:
        plot_rotated(phot,ixs,
                     d_col,col,
                     histotable,
                     apply_rolling_gaussian=apply_rolling_gaussian,
                     d_col_rot=d_col_rot,
                     bins=bins,
                     bin_weights_flag=bin_weights_flag,
                     sp=sp,
                     spi=spi,
                     histolim = [bincenter4maxval-8,bincenter4maxval+8],
                     title=f'slope:{slope}')

def rotate_d_and_find_binmax(phot,ixs,d_col,col,
                             Naxis_px, # Nx or Ny, depending on col
                             apply_rolling_gaussian=True,
                             gaussian_sigma_px=0.22,
                             d_col_rot='d_rot_tmp',
                             binsize=0.02,
                             bin_weights_flag=False,
                             slope_min=-10.0/2048.0, # 
                             slope_max=10.0/2048.0, # 
                             slope_stepsize=1.0/2048,
                             showplots=0,
                             sp=None,
                             verbose=0,
                             spi=[0,1,2,3,4]):
    if verbose:
        print(f'########################\n### rotate {d_col} versus {col}')
    rot_results = pdastroclass()

    if bin_weights_flag:
        print('weighting with flux!')
        phot.t.loc[ixs,'__weights']=10**(-0.4*phot.t.loc[ixs,'mag'])
    else:
        phot.t.loc[ixs,'__weights']=None
    
    # if rolling gaussian: get window sizes!
    if apply_rolling_gaussian:
        # get the half window size for gaussian kernel in integer bins 
        gaussian_sigma = gaussian_sigma_px/binsize
        # kernel window size: 3*sigma half width
        windowsize = int(3 * gaussian_sigma * 2)
        # make sure the window size is uneven, so that the peak is centered!
        if windowsize % 2 == 0:
            windowsize+=1
        halfwindowsize = int(windowsize*0.5)+1
        if verbose:
            print(f'Applying rolling gaussian:\ngaussian_sigma_px={gaussian_sigma_px}, binsize={binsize}, gaussian_sigma(bins)={gaussian_sigma}, windowsize(bins)={windowsize} halfwindowsize(bins)={halfwindowsize}')
    
    # Loop through the slopes
    if verbose>1:
        print(f'slope min: {slope_min}, slope max: {slope_max}, slope stepsize: slope_stepsize')
    slopes = np.arange(slope_min,slope_max,slope_stepsize)
    for counter in range(len(slopes)):
        #print(f'iteration {counter} out of {len(slopes)}: slope = {slopes[counter]:.6f}')

        find_binmax_for_slope(slopes[counter],phot,ixs,d_col,col,
                              Naxis_px,
                              gaussian_sigma,
                              windowsize,
                              halfwindowsize,
                              rot_results=rot_results,
                              apply_rolling_gaussian=apply_rolling_gaussian,
                              d_col_rot=d_col_rot,
                              binsize=binsize,
                              bin_weights_flag=bin_weights_flag,
                              showplots=showplots,
                              sp=None,
                              spi=[0,1])
    # print the results        
    #rot_results.write()
    
    # find the best rotation
    maxmaxval = np.max(rot_results.t['maxval'])
    ixs_maxmax = np.where(rot_results.t['maxval']==maxmaxval)
    if (len(ixs_maxmax[0])==0):
        raise RuntimeError('BUUUUGGGG!!!!')
    elif (len(ixs_maxmax[0])>1):
        #print(f'\nWARNING!! more than one bin with maxvalue={maxmaxval}!')
        best_index=ixs_maxmax[0][0]
    else:
        best_index=ixs_maxmax[0][0]
    if verbose:
        print('####BEST:')
    rot_results.write(indices=[best_index])
    
    # make the summary plots
    if showplots>1:
        if sp is None:
            sp = initplot(2,3)
            spi=[0,1,2,4,5]
        rot_results.t.plot('slope','maxval',ax=sp[spi[0]],color='blue',title=f'{d_col}',ylabel='histogran peak value')
        rot_results.t.plot.scatter('slope','d_bestguess',ax=sp[spi[1]],color='blue',title=f'{d_col}')
        rot_results.t.plot.scatter('slope','fwhm',ax=sp[spi[2]],color='blue',title=f'{d_col}')
        sp[spi[0]].axvline(rot_results.t.loc[best_index,'slope'],  color='red',linestyle='-', linewidth=2.0)
        sp[spi[1]].axvline(rot_results.t.loc[best_index,'slope'],  color='red',linestyle='-', linewidth=2.0)
        sp[spi[2]].axvline(rot_results.t.loc[best_index,'slope'],  color='red',linestyle='-', linewidth=2.0)
    
        find_binmax_for_slope(rot_results.t.loc[best_index,'slope'],phot,ixs,d_col,col,
                              Naxis_px,
                              gaussian_sigma,
                              windowsize,
                              halfwindowsize,
                              rot_results=rot_results,
                              apply_rolling_gaussian=apply_rolling_gaussian,
                              d_col_rot=d_col_rot,
                              binsize=binsize,
                              bin_weights_flag=bin_weights_flag,
                              showplots=3,
                              sp=sp,
                              spi=spi[3:5])
        
    return(rot_results,best_index)


def sigmacut_d_rot(phot,ixs,
                   d_col,col,
                   slope,intercept,d_rot_bestguess,
                   rough_cut_px = 2.5, #This is the first rough cut:  get rid of everything d_rot_bestguess+-rough_cut_px
                   d_col_rot='d_rot_tmp',
                   binsize=0.02,
                   Nsigma=3.0,
                   percentile_cut_firstiteration=75,
                   bin_weights_flag=True,
                   showplots=0,
                   sp=None,
                   verbose=0,
                   spi=[0]):

    ### recover the slope and intercept of the best binning
    phot.t[d_col_rot] = phot.t[d_col] - f(phot.t[col],slope,intercept)
    
    # Now make the rough cut! only keep data for with dx_rotated within  d_rot_bestguess+-rough_cut_px
    ixs_roughcut = phot.ix_inrange(d_col_rot,d_rot_bestguess-rough_cut_px,d_rot_bestguess+rough_cut_px,indices=ixs)
    #d_rotated = phot.t.loc[ixs,d_col_rot]
    if verbose:
        print(f'\n####################\n### d_rotated cut (Nsigma={Nsigma})')
    if Nsigma is None or Nsigma==0.0:
        # don't do percentile cut if there are no iterations!
        percentile_cut_firstiteration = None
    #ixs_clean4average = phot_clear.ix_inrange(d_col,0,3,indices=ixs_clear_cut)
    phot.calcaverage_sigmacutloop(d_col_rot,verbose=3,indices=ixs_roughcut,Nsigma=Nsigma,percentile_cut_firstiteration=percentile_cut_firstiteration)
    ixs_cut = phot.statparams['ix_good']

    if showplots>1:
        title = f'3-sigma cut: {len(ixs_cut)} out of {len(ixs_roughcut)} left\n'
        title += f'mean = {phot.statparams["mean"]:.3f} px, stdev = {phot.statparams["stdev"]:.3f} px'
        phot.t.loc[AnotB(ixs_roughcut,ixs_cut)].plot(col,d_col_rot,style='o',ax=sp[spi[0]],color='red', ms=5 ,alpha=0.3,title=title)
        phot.t.loc[ixs_cut].plot(col,d_col_rot,style='o',ax=sp[spi[0]],color='blue', 
                                 ms=5 ,alpha=0.3,ylabel=f'{d_col} [pixels]',
                                 title=title)
        if phot.ixs_notuse is not None:
            phot.t.loc[phot.ixs_notuse].plot(col,d_col_rot,style='o',ax=sp[spi[0]],color='gray', ms=1,alpha=0.5)
        sp[spi[0]].get_legend().remove()
    
        # set the appropriate y-axis limits
        (ylim_min,ylim_max) = (phot.t.loc[ixs_roughcut,d_col_rot].min(),phot.t.loc[ixs_roughcut,d_col_rot].max())
        ylim_min -= 0.1*(ylim_max-ylim_min)
        ylim_max += 0.1*(ylim_max-ylim_min)
        sp[spi[0]].set_ylim(ylim_min,ylim_max)
        plt.tight_layout() 

    return(ixs_cut,ixs_roughcut)


def histogram_cut(phot,ixs,d_col,col,
                  Naxis_px, # Nx or Ny, depending on col
                  d_col_rot='d_rot_tmp',
                  binsize=0.02,
                  bin_weights_flag=False,
                  slope_min=-10.0/2048.0, # 
                  slope_max=10.0/2048.0, # 
                  slope_stepsize=1.0/2048,
                  #This is the first rough cut:  get rid of everything d_rot_bestguess+-Nfwhm*fwhm,
                  Nfwhm=2.0,
                  rough_cut_px_min=None,
                  rough_cut_px_max=None,
                  Nsigma=3.0,
                  showplots=0,
                  verbose=0,
                  sp=None
                  ):
    if verbose:
        print(f'### Doing histogram cut for {d_col}, slope_min:{slope_min:.6f} slope_max:{slope_max:.6f} slope_stepsize:{slope_stepsize:.6f}')
        print(f'Nfwhm={Nfwhm}, rough_cut_px_min={rough_cut_px_min}, rough_cut_px_max={rough_cut_px_max}, Nsigma={Nsigma}')
    # initialize plot
    if showplots>1:
        if sp is None:
            sp=initplot(2,3)
    else:
        sp=None

    (rot_results,best_index) = rotate_d_and_find_binmax(phot,ixs,
                                                        d_col,col,
                                                        Naxis_px,
                                                        d_col_rot=d_col_rot,
                                                        binsize=binsize,
                                                        bin_weights_flag=bin_weights_flag,
                                                        slope_min=slope_min,
                                                        slope_max=slope_max,
                                                        slope_stepsize=slope_stepsize,
                                                        showplots=showplots,
                                                        sp=sp,
                                                        verbose=verbose,
                                                        spi=[0,1,2,3,4])
    
    # Using the best dy_rotated, we first remove all entries with dy_rotated outside of dy_bestguess+-Nfwhm*fwhm
    # Note that FWHM ~ 2.355 stdev, so Nfwhm*fwhm should be at least 3*stdev. This is the first ROUGH cut, with 
    # which we just want to remove excessive amounts of outliers. Then a 3-sigma cut is done on the *rotated* dy
    rough_cut_px= Nfwhm*rot_results.t.loc[best_index,'fwhm']
    # when using a rolling gaussian to smooth the histogram (in particular for small numbers of good matches), 
    # the Nfwhm*fwhm method does nto work well. In that case it is better to use upper/lower limits
    if verbose:
        print(f'Setting rough_cut_px={rough_cut_px}. limits: ({rough_cut_px_min}-{rough_cut_px_max})')
    if (rough_cut_px_max is not None) and rough_cut_px>rough_cut_px_max:
        rough_cut_px=rough_cut_px_max
    if (rough_cut_px_min is not None) and rough_cut_px<rough_cut_px_min:
        rough_cut_px=rough_cut_px_min
    if verbose:
        print(f'Setting rough_cut_px={rough_cut_px}')
    
    (ixs_cut,ixs_roughcut) = sigmacut_d_rot(phot,ixs,d_col,col,
                                            rot_results.t.loc[best_index,'slope'],
                                            rot_results.t.loc[best_index,'intercept'],
                                            rot_results.t.loc[best_index,'d_bestguess'],
                                            rough_cut_px =rough_cut_px ,
                                            binsize=binsize,
                                            Nsigma=Nsigma,
                                            bin_weights_flag=bin_weights_flag,
                                            showplots=showplots,
                                            sp=sp,
                                            spi=[5],
                                            verbose=verbose
                                            )
    plt.tight_layout()   
    if showplots>0:
        plt.show()
    else:
        plt.close()
    return(ixs_cut,rot_results)



class st_wcs_align:
    """
    Main class for alignment.

    Attributes
    ----------
    outrootdir : str
        output root directory. The output directoy is the output root directory + the outsubdir if not None
    outsubdir : str
            outsubdir added to output root directory
    overwrite : bool
            overwrite files if they exist.
    telescope : str
            If None, then telescope is determined automatically from the filename ("jw*" and "hst*" for JWST and HST, respectively)
    skip_if_exists : bool
            Skip doing the analysis of a given input image if the cal file already exists, assuming the full analysis has been already done
    verbose : int 
            verbosity count (lower is less verbose)
    SNR_min : float 
            mininum SNR for objects in image to be used for analysis
    use_dq : bool
            use the DQ extensions for masking
    refcat : str
            reference catalog. Can be a filename or Gaia 
    refcat_racol : str
            RA column of reference catalog. If None, then automatically determined 
    refcat_deccol : str
            Dec column of reference catalog. If None, then automatically determined
    refcat_magcol : str
            mag column of reference catalog. If None and not one of the default refcats like gaia, then 3rd column is used
    refcat_magerrcol : str
            magerr column of reference catalog. If None, then not used 
    refcat_colorcol : str:
            color column of reference catalog. If None, then not used
    refcat_pmflag : bool
            Apply the proper motion correction (only for catalogs it is applicable, e.g., gaia
    refcat_pmmedian : bool
            Apply the MEDIAN proper motion correction (only for catalogs it is applicable, e.g., gaia
    photfilename : str
            photometry output filename. if "auto", the fits in the image filename is substituted with phot.txt 
    load_photcat_if_exists : bool
            If the photometric catalog file already exists, skip recreating it.
    rematch_refcat : bool
            if --load_photcat_if_exists and the photcat already exists, load the photcat, but rematch with refcat
    d2d_max : float
            maximum distance between source in image and refcat object, in arcsec
    dmag_max : float
            maximum uncertainty of sources in image 
    sharpness_lim : list of float
            sharpness limits of sources in image (iterable of length 2)
    roundness1_lim : list of float
            roundness1 limits of sources in image (iterable of length 2)
    delta_mag_lim : list of float
            limits on mag - refcat_mainfilter (iterable of length 2)
    objmag_lim : list of float
            limits on mag, the magnitude of the objects in the image (iterable of length 2)
    refmag_lim : list of float
            limits on refcat_mainfilter, the magnitude of the reference catalog (iterable of length 2)
    slope_min : float
            minimum slope for linear correction applied to dx/dy. This effectively accounts for rotation. slopes go from slopemin to -slopemin
    Nbright4match : int 
            Use only Nbright brightest objects for matching to the ref cat 
    Nbright : int
            Use only Nbright brightest objects in image that are matched to refcat for alignment 
    histocut_order : str
            histocut_order defines whether the histogram cut is first done for dx or first for dy (choices are 'dxdy' or 'dydx')
    xshift : float
            added to the x coordinate before calculating ra,dec (only impacts ra,dec, not x). This can be used to correct for large shifts before matching!
    yshift : float
            added to the y coordinate before calculating ra,dec (only impacts ra,dec, not y). This can be used to correct for large shifts before matching! 
    iterate_with_xyshifts: bool
            After the first histogram fit, redo the match with refcat with x/yshift=median(dx/dy) and redo histofit. Use this if the offsets are big, since the second iteration will give you better matching with the refcat
    showplots : int 
            showplots=1: most important plots. showplots=2: all plots (debug/test/finetune)
    saveplots : int 
            saveplots=1: most important plots. saveplots=2: all plots (debug/test/finetune)
    savephottable : bool
            Save the final photometry table
    replace_sip : bool
            Replace the tweaked fits image wcs with the SIP representation.
    sip_err : float
            max_pix_error for SIP transformation.
    sip_degree : int 
            degree for SIP transformation.
    sip_points : int
            npoints for SIP transformation.
    ee_radius : int 
            encircled energy percentage (multiples of 10) for photometry
    rough_cut_px_min : float
            first rough cut: best d_rotated+-rough_cut_pix. This is the lower limit for rough_cut
    rough_cut_px_max : float
            first rough cut: best d_rotated+-rough_cut_pix. This is the upper limit for rough_cut 
    d_rotated_Nsigma : float
            Nsigma for sigma cut of d_rotated. If 0.0, then 3-sigma cut is skipped 
    """
    def __init__(self):
        self.verbose=0
        self.outdir = None
        
        self.telescope = None
        
        self.phot=jwst_photclass()
        self.phot.ixs4use=None
        
        self.replace_sip = True
        self.sip_err = 0.1
        self.sip_degree = 3
        self.sip_points = 128

        self.rough_cut_px_min=0.3
        self.rough_cut_px_max=0.8

        self.d_rotated_Nsigma=3.0        

    def define_options(self,parser=None,usage=None,conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)

        #parser.add_argument('--rate_dir', default=ratedir, help='Directory in which the rate images are located, which will be used to test the distortions. (default=%(default)s)')
        parser.add_argument('cal_image',  help='cal.fits filename or any other image that is at a similar reduction stage.')

        parser = self.default_options(parser)
        return(parser)

    def default_options(self,parser):

        # default directory for input images
        #if 'JWST_DISTORTION_IMAGEDIR' in os.environ:
        #    ratedir = os.environ['JWST_DISTORTION_IMAGEDIR']
        #else:
        #    ratedir = None

        # default directory for output
        if 'JWST_OUTPUT_ROOTDIR' in os.environ:
            outrootdir = os.environ['JWST_OUTPUT_ROOTDIR']
        else:
            outrootdir = None

        parser.add_argument('--outrootdir', default=outrootdir, help='output root directory. The output directoy is the output root directory + the outsubdir if not None.  If $JWST_OUTPUT_ROOTDIR is defined, then this dir is taken as default (default=%(default)s)')
        parser.add_argument('--outsubdir', default=None, help='outsubdir added to output root directory (default=%(default)s)')
        parser.add_argument('--overwrite', default=False, action='store_true', help='overwrite files if they exist.')

        parser.add_argument('--telescope', default=None, help='If None, then telescope is determined automatically from the filename ("jw*" and "hst*" for JWST and HST, respectively) (default=%(default)s)')

        parser.add_argument('--skip_if_exists', default=False, action='store_true', help='Skip doing the analysis of a given input image if the cal file already exists, assuming the full analysis has been already done')

        parser.add_argument('-v','--verbose', default=0, action='count')

        parser.add_argument('--SNR_min', default=None,type=float, help='mininum SNR for object in image to be used for analysis (default=%(default)s)')

        parser.add_argument('--use_dq', default=False, action='store_true', help='use the DQ extensions for masking')


        parser.add_argument('--refcat', default='Gaia', help='reference catalog. Can be a filename or Gaia (default=%(default)s)')
        parser.add_argument('--refcat_racol', default=None, help='RA column of reference catalog. If None, then automatically determined (default=%(default)s)')
        parser.add_argument('--refcat_deccol', default=None, help='Dec column of reference catalog. If None, then automatically determined (default=%(default)s)')
        parser.add_argument('--refcat_magcol', default=None, help='mag column of reference catalog. If None and not one of the default refcats like gaia, then 3rd column is used (default=%(default)s)')
        parser.add_argument('--refcat_magerrcol', default=None, help='magerr column of reference catalog. If None, then not used  (default=%(default)s)')
        parser.add_argument('--refcat_colorcol', default=None, help='color column of reference catalog. If None, then not used (default=%(default)s)')
        parser.add_argument('--refcat_pmflag', default=False, action='store_true', help='Apply the proper motion correction (only for catalogs it is applicable, e.g., gaia')
        parser.add_argument('--refcat_pmmedian', default=False, action='store_true', help='Apply the MEDIAN proper motion correction (only for catalogs it is applicable, e.g., gaia')
        parser.add_argument('--photfilename', default='auto', help='photometry output filename. if "auto", the fits in the image filename is substituted with phot.txt (default=%(default)s)')
#        parser.add_argument('--photfilename', default='auto', help='photometry output filename. if "auto", the fits in the image filename is substituted with phot.txt (default=%(default)s)')

        parser.add_argument('--load_photcat_if_exists', default=False, action='store_true', help='If the photometric catalog file already exists, skip recreating it.')
        parser.add_argument('--rematch_refcat', default=False, action='store_true', help='if --load_photcat_if_exists and the photcat already exists, load the photcat, but rematch with refcat')

        parser.add_argument('--d2d_max', default=None, type=float, help='maximum distance between source in image and refcat object, in arcsec (default=%(default)s)')
        parser.add_argument('--dmag_max', default=None, type=float, help='maximum uncertainty of sources in image (default=%(default)s)')
        parser.add_argument('--sharpness_lim', default=(None,None), nargs=2, type=float, help='sharpness limits of sources in image (default=%(default)s)')
        parser.add_argument('--roundness1_lim', default=(-0.75,0.75), nargs=2, type=float, help='roundness1 limits of sources in image (default=%(default)s)')
        parser.add_argument('--delta_mag_lim', default=(None,None), nargs=2, type=float, help='limits on mag - refcat_mainfilter (default=%(default)s)')
        parser.add_argument('--objmag_lim', default=(None,None), nargs=2, type=float, help='limits on mag, the magnitude of the objects in the image (default=%(default)s)')
        parser.add_argument('--refmag_lim', default=(None,None), nargs=2, type=float, help='limits on refcat_mainfilter, the magnitude of the reference catalog (default=%(default)s)')
        parser.add_argument('--slope_min', default=-0.008, type=float, help='minimum slope for linear correction applied to dx/dy. This effectively accounts for rotation. slopes go from slopemin to -slopemin (default=%(default)s)')
        parser.add_argument('--Nbright4match', default=None, type=int, help='Use only Nbright brightest objects for matching to the ref cat (default=%(default)s)')
        parser.add_argument('--Nbright', default=None, type=int, help='Use only Nbright brightest objects in image that are matched to refcat for alignment (default=%(default)s)')
        parser.add_argument('--histocut_order', default='dxdy', choices=['dxdy','dydx'], help='histocut_order defines whether the histogram cut is first done for dx or first for dy (default=%(default)s)')
        parser.add_argument('--xshift', default=0.0, type=float, help='added to the x coordinate before calculating ra,dec (only impacts ra,dec, not x). This can be used to correct for large shifts before matching! (default=%(default)s)')
        parser.add_argument('--yshift', default=0.0, type=float, help='added to the y coordinate before calculating ra,dec (only impacts ra,dec, not y). This can be used to correct for large shifts before matching! (default=%(default)s)')
        parser.add_argument('--iterate_with_xyshifts', default=False, action='store_true',help='After the first histogram fit, redo the match with refcat with x/yshift=median(dx/dy) and redo histofit. Use this if the offsets are big, since the second iteration will give you better matching with the refcat')
        parser.add_argument('-p','--showplots', default=0, action='count',help='showplots=1: most important plots. showplots=2: all plots (debug/test/finetune)')
        parser.add_argument('--saveplots', default=0, action='count',help='saveplots=1: most important plots. saveplots=2: all plots (debug/test/finetune)')
        parser.add_argument('-t','--savephottable', default=0, action='count',help='Save the final photometry table')
        parser.add_argument('--replace_sip', default=True, action='store_true',help='Replace the tweaked fits image wcs with the SIP representation.')
        parser.add_argument('--sip_err', default=0.3, type=float,help='max_pix_error for SIP transformation.')
        parser.add_argument('--sip_degree', default=3, type=int,help='degree for SIP transformation.')
        parser.add_argument('--sip_points', default=128, type=int,help='npoints for SIP transformation.')
        parser.add_argument('--ee_radius', default=70, type=int, help='encircled energy percentage (multiples of 10) for photometry')
        parser.add_argument('--is_hst', default=False, action='store_true', help='set if your image is from hst not jwst')
        parser.add_argument('--rough_cut_px_min', default=0.3, type=float,help='first rough cut: best d_rotated+-rough_cut_pix. This is the lower limit for rough_cut (default=%(default)s)')
        parser.add_argument('--rough_cut_px_max', default=1.0, type=float,help='first rough cut: best d_rotated+-rough_cut_pix. This is the upper limit for rough_cut (default=%(default)s)')
        parser.add_argument('--d_rotated_Nsigma', default=3.0, type=float,help='Nsigma for sigma cut of d_rotated. If 0.0, then 3-sigma cut is skipped (default=%(default)s)')

        return(parser)
    
    def set_telescope(self,telescope=None,imname=None):
        """
        Figuring out which telescope your data come from (HST or JWST).
        Note that if your filename is non-standard, then you MUST set
        """

        if telescope is None:
            if imname is None:
                raise RuntimeError('neither telescope not image name is specified, cannot determine telescope')
            if re.search('^jw',os.path.basename(imname)):
                telescope = 'JWST'
            elif re.search('^hst',os.path.basename(imname)):
                telescope = 'HST'
            else:
                raise RuntimeError(f'Cannot parse image name {imname} to determine telescope! Please set the "telescope" argument')                
        self.telescope = telescope
        if telescope.lower() == 'jwst':
            self.phot=jwst_photclass()
        elif telescope.lower() == 'hst':
            self.phot=hst_photclass()
            self.replace_sip = True
        else:
            raise RuntimeError(f'unknown telescope {telescope}')
        self.phot.verbose = self.verbose
        self.phot.ixs4use=None            
        if self.verbose: print(f'telescope set to {self.telescope}')

    def set_outdir(self,outrootdir=None,outsubdir=None):
        self.outdir = outrootdir
        if self.outdir is None: self.outdir = '.'
        
        if outsubdir is not None:
            self.outdir+=f'/{outsubdir}'
        return(self.outdir) 
         
    def set_outbasename(self,outrootdir=None,outsubdir=None,outbasename=None,inputname=None):

        if outbasename is not None:
            if outrootdir is not None or outsubdir is not None:
                raise RuntimeError(f'Cannot specify both, outbasename={outbasename} and outrootdir/outsubdir={outrootdir}/{outsubdir}')
            self.set_outdir(outrootdir=os.path.dirname(outbasename))
            self.outbasename = f'{self.outdir}/{os.path.basename(outbasename)}'
        else:
            if inputname is not None:
                inputname = os.path.basename(inputname)
                inputbasename = re.sub('_([a-zA-Z0-9]+)\.fits$','',inputname)
                #if inputbasename == inputname:
                #    raise RuntimeError(f'Cound not strip ".fits" from {inputname}')
                self.set_outdir(outrootdir=outrootdir,outsubdir=outsubdir)
                self.outbasename = f'{self.outdir}/{inputbasename}'
            else:
                RuntimeError('Could not determine outputbasename, neither image nor basename passed to routine')
        
        return(self.outbasename)
        
    def run_align2refcat(self,imfilename,
                         outputfits=None,
                         phot=None,
                         ixs=None,
#                         refcat_short=None,
                         refcat_racol=None,
                         refcat_deccol=None,
                         xcol='x',
                         ycol='y',
                         outdir=None,overwrite=False, 
                         skip_if_exists=False,
                         savephot=True,
                         ):
            
        if phot is None:
            phot=self.phot
        if outdir is None:
            if self.outdir is None:
                self.set_outdir()
            outdir=self.outdir
        
        
#        if (racol is None) or (deccol is None): 
#            if refcat_short is None:
#                refcat_short = phot.refcat.short
#                if refcat_short is None:
#                    raise RuntimeError('need the refcat shortname to identify ra,dec columns!')
#            racol = f'{refcat_short}_ra'
#            deccol = f'{refcat_short}_dec'
            
        if refcat_racol is None: refcat_racol = phot.refcat_racol
        if refcat_deccol is None: refcat_deccol = phot.refcat_deccol

        tweakreg = tweakreg_hack.TweakRegStep()
        if self.telescope=='jwst' or not self.phot.do_driz:
            cal_image = datamodels.open(imfilename)
            tweakreg.do_driz = False
            tweakreg.input_file = imfilename
        else:
            cal_image = datamodels.open(self.phot.driz_name)
            tweakreg.do_driz = True
            tweakreg.true_file = imfilename
            tweakreg.input_file = self.phot.driz_name

        if outputfits is None:
            shortoutputfits = re.sub('_[a-zA-Z0-9]+\.fits$','_jhat.fits',os.path.basename(imfilename))
            if shortoutputfits == os.path.basename(imfilename):
                raise RuntimeError('could not determine output filename for {imfilename}')
            # if outdir is None, try to determine it!
            if outdir is None:
                outdir = self.outdir
            if outdir is None:
                outdir = os.path.dirname(imfilename)
            if outdir is None:
                raise RuntimeError('Could not determine outdir!')
            outputfits = f'{outdir}/{shortoutputfits}'
        else:
            (outdir,shortoutputfits) = os.path.split(outputfits)
            
            
        if self.verbose:    
            print(f'Setting output directory for {shortoutputfits} file to {outdir}')
        tweakreg.output_dir = outdir
        if not os.path.isdir(outdir):
            makepath(outdir)
        
        if os.path.isfile(outputfits):
            if not overwrite:
                if skip_if_exists:
                    # return False means that rate2cal did not run
                    print(f'Image {outputfits} already exists, skipping recreating it...')
                    return(False,outputfits)
                else:
                    raise RuntimeError(f'Image {outputfits} already exists! exiting. If you want to overwrite or skip rate2cal reduction, you can use "overwrite" or "skip_if_exists"')
            else:
                # make sure output filename is deleted
                rmfile(outputfits)


        # It is important to set fitgeometry to rshift for level 2
        tweakreg.pipeline_level = self.phot.pipeline_level
        if self.phot.pipeline_level==2 and not (self.telescope=='hst' and self.phot.do_driz):
            tweakreg.fitgeometry = 'rshift'
        else:    
            tweakreg.fitgeometry = 'general'
        tweakreg.align_to_gaia = False
        tweakreg.save_gaia_catalog = False
        tweakreg.save_results = True
        # minimum number of objects required for fit
        tweakreg.save_catalogs = False
        tweakreg.minobj = 4
        
        # the following parameters should have not impact, since these steps in tweakreg are skipped
        if 1==1:
            tweakreg.snr_threshold = 50
            tweakreg.separation = 9
            tweakreg.searchrad = 0.5
            tweakreg.min_gaia = 30
            tweakreg.xoffset = 0
            tweakreg.yoffset = 0
            tweakreg.brightest = 4000        
        
        tweakreg.already_matched = True
        # phot_cal.t.loc[ixs_cal_good] is the table with the good matches!
        if self.verbose: print(f'{len(ixs)} matches are passed to tweakreg {tweakreg.fitgeometry} fitting')
        t =  Table.from_pandas(phot.t.loc[ixs,[xcol,ycol,refcat_racol,refcat_deccol]])
        #t.write(os.path.join(outdir,shortoutputfits.replace('.fits','_matched.txt')),format='ascii')
        tweakreg.refcat = t
        tweakreg.ref_racol = refcat_racol
        tweakreg.ref_deccol = refcat_deccol
        
        ### Provide your own source catalog, to be used in place of the default daofinder stuff. If you actually have a list
        ### of images, it's okay to provide a source catalog for each. 
        cal_image.source_catalog = t
        cal_data = [cal_image]
        tweakreg.source_xcol = xcol
        tweakreg.source_ycol = ycol
        
        if self.verbose: print(f'Fitting tweakreg fitgeometry={tweakreg.fitgeometry} to xy={xcol},{ycol} to ra,dec={refcat_racol},{refcat_deccol}')
        
        cal_data = [datamodels.open(cal_image)]
        
        tweakreg.output_file = os.path.join(outdir,shortoutputfits)
        tweakreg.telescope = self.telescope
        tweakreg.run(cal_data)

        if not os.path.isfile(outputfits) and not \
            os.path.isfile(outputfits.replace('jhat.fits','tweakregstep.fits')):
            raise RuntimeError(f'Image {outputfits} did not get created!!')
        elif not os.path.isfile(outputfits):
            os.rename(outputfits.replace('jhat.fits','tweakregstep.fits'),outputfits)
        if self.replace_sip and self.phot.pipeline_level==2:
            if self.telescope.lower()=='jwst':
                print('replacing SIP',outputfits)
                dm = datamodels.open(outputfits)
                gwcs_header = dm.meta.wcs.to_fits_sip(max_pix_error=self.sip_err,
                                                       max_inv_pix_error=self.sip_err,
                                                       degree=self.sip_degree,
                                                       npoints=self.sip_points)
                from astropy.io import fits
                dm_fits = fits.open(outputfits)
                dm_fits['SCI',1].header.update(gwcs_header)
                # for key,value in dict(gwcs_header).items():
                #     for k in dm_fits['SCI',1].header.keys():
                #         if k==key:
                #             dm_fits['SCI',1].header[key] = value
                #             break
                #     #astropy.wcs.WCS(header=gwcs_header)
                dm_fits.writeto(outputfits,overwrite=True)
        #print(imcat.meta['image_model'].wcs)
        #print(gwcs_header)
        #imcat.meta['image_model'].wcs = astropy.wcs.WCS(header=gwcs_header)
        #make sure the image got created
        

        # return True means that tweakrun did run
        return(True,outputfits)
    

    def find_good_refcat_matches(self,
                                 phot=None,
                                 ixs=None,
                                 refcat_xcol=None,
                                 refcat_ycol=None,
                                 xcol='x',
                                 ycol='y',
                                 d2d_max = None,
                                 dmag_max =None,
                                 sharpness_lim = (None, None), # sharpness limits
                                 roundness1_lim = (None, None), # roundness1 limits 
                                 delta_mag_lim = (None, None), # limits on mag-refcat_mainfilter
                                 objmag_lim = (None,None), # limits on mag, the magnitude of the objects in the image
                                 refmag_lim = (None,None), # limits on refcat_mainfilter, the magnitude of the reference catalog                    
                                 Nbright=None,
                                 show_initial_plot=1,
                                 show_histofit_plots=1,
                                 savephottable=1,
                                 outbasename=None,
                                 plots_dxdy_delta_pix_ylim=20,
                                 # histogram cut parameters
                                 histocut_order = 'dxdy', # this can only be 'dxdy' or 'dydx'
                                 binsize_px = 0.2, # this is the binsize of the dx/dy histograms
                                 bin_weights_flag=False,# If bin_weights_flag is set to True, 
                                                       #then the dx/dy bins are weighted by 
                                                       # the flux of the detection.
                                 slope_min=-10/2048.0, 
                                 slope_Nsteps = 200, # slope_max=-slope_min, slope_stepsize=(slope_max-slope_min)/slope_Nsteps
                                 Nfwhm = 2.5 
                                 ):
        if phot is None:
            phot=self.phot
            
        if savephottable and (outbasename is None):
            raise RuntimeError('Trying to save plots and/or phot tables, but outbasename is None!')

        if refcat_xcol is None: refcat_xcol = phot.refcat_xcol
        if refcat_ycol is None: refcat_ycol = phot.refcat_ycol

        Nx = phot.scihdr['NAXIS1']
        Ny = phot.scihdr['NAXIS2']
        
        # Calculate dx and dy
        phot.t['dx'] = phot.t[refcat_xcol] - phot.t[xcol]
        phot.t['dy'] = phot.t[refcat_ycol] - phot.t[ycol]

        #phot.write('delme1.txt',columns=['x','y','dx','dy','gaia_x','gaia_y'])
        #phot.write(columns=['x','y','dx','dy','gaia_x','gaia_y'])
        #print(len(phot.ixs_use),len(phot.ixs_notuse),len(phot.t))
        
        # Calculate the difference between JWST mag and main filter of reference catalog
        if phot.refcat_mainfilter is not None:
            phot.t['delta_mag'] = phot.t['mag'] - phot.t[phot.refcat_mainfilter]
        else:
            phot.t['delta_mag'] = np.nan
        
        # do some first very rough cuts on matches
        ixs = self.phot.initial_cut_matches(d2d_max=d2d_max,
                                            delta_mag_lim = delta_mag_lim, # limits on mag-refcat_mainfilter
                                            Nbright=Nbright,
                                            ixs=ixs)
        
        if len(ixs)<4:
            raise RuntimeError(f'Only {len(ixs)} objects pass the initial cut, at least 3 required!')
        
        # do the initial dx,dy plot and other important plots
        # it shows the initial cut.
        if show_initial_plot:
            initial_dxdy_plot(phot, phot.ixs_use, phot.ixs_notuse,
                              plots_dxdy_delta_pix_ylim=plots_dxdy_delta_pix_ylim,
                              refcat_mainfilter=phot.refcat_mainfilter,refcat_mainfilter_err=phot.refcat_mainfilter_err,refcat_maincolor=phot.refcat_maincolor,
                              d2d_max=d2d_max,dmag_max=dmag_max,Nbright=Nbright,delta_mag_lim=delta_mag_lim,verbose=self.verbose)

       

        # Here we correct for rotation so that we can robustly remove outlier matches
        
        # A straight line with a given slope is subtracted from dx versus y, which is 
        # effectively a correction for a rotation.
        
        # Then a histogram of the resulting dx_rotated is calculated with bin size = binsize_px in pixels.
        # The maximum value of the histogram, it's associated dx_rot position, a rough upper limit estimate
        # of the fwhm of the histogram peak are then saved for each slope and intercept in dx_rot_results (pdastro table).
        
        # best_index is the index for this table for which the largest maxval of the histograms, presumably for 
        # the best rotation (or slope)
        
        # slope loops from slope_min to slope_max=-slope_min with slope_stepsize
        # slope_min=-10.0/2048.0 means that slopes with a range 
        # of +-10 pixels of 2048 pixels are probed, i.e. dx does not change
        # by more than 10 pixels over the full detector width
        #  slope_stepsize=(slope_max-slope_min)/slope_Nsteps
 
        slope_max=-slope_min
        slope_stepsize=(slope_max-slope_min)/slope_Nsteps
       
        # histocut_order defines whether the histogram cut is first done for dx or first for dy
        if histocut_order == 'dxdy':
            d_col1,col1,Naxis1_px = 'dx','y',Ny
            d_col2,col2,Naxis2_px = 'dy','x',Nx
        elif histocut_order == 'dydx':
            d_col1,col1,Naxis1_px = 'dy','x',Nx
            d_col2,col2,Naxis2_px = 'dx','y',Ny


        # Do the histogram cut on the first dcol (dx or dy, as selected)
        (ixs_cut1,rot_results1) = histogram_cut(phot,ixs,d_col1,col1,Naxis1_px,
                                                binsize=binsize_px,
                                                bin_weights_flag=bin_weights_flag,
                                                slope_min=slope_min,
                                                slope_max=slope_max,
                                                slope_stepsize=slope_stepsize,
                                                Nfwhm=Nfwhm,
                                                rough_cut_px_min=self.rough_cut_px_min,
                                                rough_cut_px_max=self.rough_cut_px_max,
                                                Nsigma=self.d_rotated_Nsigma,
                                                showplots=show_histofit_plots,
                                                verbose=self.verbose)

        # Do the histogram cut on the second dcol (dx or dy, as selected)
        (ixs_cut2,rot_results2) = histogram_cut(phot,ixs_cut1,d_col2,col2,Naxis2_px,
                                                binsize=binsize_px,
                                                bin_weights_flag=bin_weights_flag,
                                                slope_min=slope_min,
                                                slope_max=slope_max,
                                                slope_stepsize=slope_stepsize,
                                                Nfwhm=Nfwhm,
                                                rough_cut_px_min=self.rough_cut_px_min,
                                                rough_cut_px_max=self.rough_cut_px_max,
                                                Nsigma=self.d_rotated_Nsigma,
                                                showplots=show_histofit_plots,
                                                verbose=self.verbose)


        if savephottable:
            #print(f'Saving {outbasename}.good.phot.txt')
            phot.write(f'{outbasename}.goodmatches.phot.txt',indices=ixs_cut2,verbose=1)
            if savephottable>1:
                #print(f'Saving {outbasename}.all.phot.txt')
                phot.write(f'{outbasename}.allmatches.phot.txt',verbose=1)

        if show_histofit_plots>1:
            plt.tight_layout()
            print('*** Note: close plots to continue!')
            plt.show()

            # get the bad data points
        #    infoplots(phot,ixs_dy_cut,dy_plotlim=dy_plotlim,dx_plotlim=dx_plotlim)
            

        return(ixs_cut2)
                
    def update_phottable_final_wcs(self,outputfits,
                                   ixs_bestmatch,
                                   phot=None,
                                   refcat_racol=None,
                                   refcat_deccol=None,
                                   refcat_xcol=None,
                                   refcat_ycol=None,                                   
                                   showplots=1,
                                   saveplots=1,
                                   savephottable=1,
                                   overwrite=False):
        outbasename = re.sub('\.fits$','',outputfits)
        if (outbasename == outputfits): raise RuntimeError(f'Could not remove .fits from {outputfits}')        
        if phot is None:
            phot=self.phot

        if refcat_racol is None: refcat_racol = phot.refcat_racol
        if refcat_deccol is None: refcat_deccol = phot.refcat_deccol
        if refcat_xcol is None: refcat_xcol = phot.refcat_xcol
        if refcat_ycol is None: refcat_ycol = phot.refcat_ycol


        # show or save dxdy pre WCS correction
        if showplots>=0 or saveplots:
            dxdy_plot(phot, ixs_bestmatch,title='pre WCS correction',
                      refcat_mainfilter=phot.refcat_mainfilter,
                      refcat_mainfilter_err=phot.refcat_mainfilter_err,
                      refcat_maincolor=phot.refcat_maincolor
                     )
            if saveplots:
                outfilename = f'{outbasename}.phot.prewcs.png'
                if os.path.isfile(outfilename):
                    rmfile(outfilename)
                if self.verbose:
                    print(f'Saving {outfilename}')
                plt.savefig(outfilename)
            if showplots>0:
                plt.show()
            plt.close()

        #racol=f'{phot.refcat.short}_ra'
        #deccol=f'{phot.refcat.short}_dec'
        #xcol=f'{phot.refcat.short}_x'
        #ycol=f'{phot.refcat.short}_y'

        phot.t.drop(columns=['dx','dy',refcat_xcol,refcat_ycol],inplace=True)    
        if f'{phot.refcatshort}_d2d' in phot.t.columns:
            phot.t.drop(columns=[f'{phot.refcatshort}_d2d'],inplace=True)    
            

        #image_model = datamodels.ImageModel(outputfits)
        image_model = fits.open(outputfits)
        print(outputfits)
        # recalculate the x,y of the ref cat objects
        #world_to_detector = image_model.meta.wcs.get_transform('world', 'detector')
        #phot.t[refcat_xcol], phot.t[refcat_ycol] = world_to_detector(phot.t[refcat_racol],phot.t[refcat_deccol])
        from astropy import wcs
        from astropy.wcs.utils import skycoord_to_pixel,pixel_to_skycoord
        from astropy import units as u
        imwcs = wcs.WCS(image_model['SCI',1],image_model)
        if self.telescope=='hst' and self.phot.do_driz:
            old_model = fits.open(self.phot.imagename)
            oldwcs1 = wcs.WCS(old_model['SCI',1],old_model)
            oldwcs2 = wcs.WCS(old_model['SCI',2],old_model)
            x1,y1 = skycoord_to_pixel(SkyCoord(phot.t[refcat_racol],
            phot.t[refcat_deccol],unit=u.deg),oldwcs1)
            x2,y2 = skycoord_to_pixel(SkyCoord(phot.t[refcat_racol],
            phot.t[refcat_deccol],unit=u.deg),oldwcs2)
            inds2 = np.where(np.logical_or(np.logical_or(x1<0,x1>old_model['SCI',1].data.shape[1]),
                np.logical_or(y1<0,y1>old_model['SCI',1].data.shape[0])))[0]
            for i in inds2:
                x1[i] = x2[i]
                y1[i] = y2[i]
            #for i in range(len(x1)):
            #    if not 0<x1[i]<old_model['SCI',1].data.shape[1] or not\
            #        0<y1[i]<old_model['SCI',1].data.shape[0]:
            #        x1[i],y1[i] = skycoord_to_pixel(SkyCoord(phot.t[refcat_racol][i],
            #                        phot.t[refcat_deccol][i],unit=u.deg),oldwcs2)
            phot.t['x'] = x1
            phot.t['y'] = y1


            newwcs1 = wcs.WCS(image_model['SCI',1],image_model)
            newwcs2 = wcs.WCS(image_model['SCI',2],image_model)
            x1,y1 = skycoord_to_pixel(SkyCoord(phot.t[refcat_racol],
            phot.t[refcat_deccol],unit=u.deg),newwcs1)
            x2,y2 = skycoord_to_pixel(SkyCoord(phot.t[refcat_racol],
            phot.t[refcat_deccol],unit=u.deg),newwcs2)
            inds2 = np.where(np.logical_or(np.logical_or(x1<0,x1>image_model['SCI',1].data.shape[1]),
                np.logical_or(y1<0,y1>image_model['SCI',1].data.shape[0])))[0]
            sc1 = pixel_to_skycoord(x1,y1,newwcs1)
            sc2 = pixel_to_skycoord(x2,y2,newwcs2)
            phot.t['ra'] = sc1.ra.value
            phot.t['dec'] = sc1.dec.value
            for i in inds2:
                x1[i] = x2[i]
                y1[i] = y2[i]
                phot.t['ra'][i] = sc2[i].ra.value
                phot.t['dec'][i] = sc2[i].dec.value
            #for i in range(len(x1)):
            #    if not 0<x1[i]<old_model['SCI',1].data.shape[1] or not\
            #        0<y1[i]<old_model['SCI',1].data.shape[0]:
            #        x1[i],y1[i] = skycoord_to_pixel(SkyCoord(phot.t[refcat_racol][i],
            #                        phot.t[refcat_deccol][i],unit=u.deg),oldwcs2)
            phot.t[refcat_xcol] = x1
            phot.t[refcat_ycol] = y1
        else:
            phot.t[refcat_xcol], phot.t[refcat_ycol] = skycoord_to_pixel(SkyCoord(phot.t[refcat_racol],
                phot.t[refcat_deccol],unit=u.deg),imwcs)
            sc = pixel_to_skycoord(phot.t['x'],phot.t['y'],imwcs)
            phot.t['ra'] = sc.ra.value
            phot.t['dec'] = sc.dec.value

        # recalculate dx, dy
        phot.t['dx'] = phot.t[refcat_xcol] - phot.t['x']
        phot.t['dy'] = phot.t[refcat_ycol] - phot.t['y']

        # recalculate the RA, Dec of the image objects
        #detector_to_world = image_model.meta.wcs.get_transform('detector', 'world')

        #phot.t['ra'],phot.t['dec'] = detector_to_world(phot.t['x'],phot.t['y'])
        

        if savephottable:
            outfilename = f'{outbasename}.good.phot.txt'
            if os.path.isfile(outfilename):
                if not overwrite:
                    raise RuntimeError(f'Image {outfilename} already exists! exiting. If you want to overwrite or skip rate2cal reduction, you can use "overwrite" or "skip_if_exists"')
                else:
                    # make sure output filename is deleted
                    rmfile(outfilename)
            if self.verbose:
                print(f'Saving {outfilename}')
            phot.write(outfilename,indices=ixs_bestmatch)
            if savephottable>2:
                outfilename = f'{outbasename}.all.phot.txt'
                if os.path.isfile(outfilename):
                    if not overwrite:
                        raise RuntimeError(f'Image {outfilename} already exists! exiting. If you want to overwrite or skip rate2cal reduction, you can use "overwrite" or "skip_if_exists"')
                    else:
                        # make sure output filename is deleted
                        rmfile(outfilename)
                if self.verbose:        
                    print(f'Saving {outfilename}')
                phot.write(outfilename)

        # show or save dxdy post WCS correction
        if showplots>=0 or saveplots:
            print('should be plotting')
            dxdy_plot(phot, ixs_bestmatch,title='after WCS correction',
                      refcat_mainfilter=phot.refcat_mainfilter,
                      refcat_mainfilter_err=phot.refcat_mainfilter_err,
                      refcat_maincolor=phot.refcat_maincolor)
            if saveplots:
                outfilename = f'{outbasename}.phot.finalwcs.png'
                if os.path.isfile(outfilename):
                    rmfile(outfilename)
                if self.verbose:
                    print(f'Saving {outfilename}')
                plt.savefig(outfilename)
            if showplots>0:
                plt.show()
            plt.close()



    def run_all(self,input_image,
                telescope=None,
                #distortion_file=None,
                outrootdir = None,
                outsubdir = None,
                overwrite = False,
                skip_if_exists = False,
                #skip_applydistortions_if_exists = False,
                use_dq=False,
                # refcat parameters
                refcatname = 'Gaia',
                refcat_racol='auto',
                refcat_deccol='auto',
                refcat_magcol = None,
                refcat_magerrcol = None,
                refcat_colorcol = None,
                pmflag = False,
                pm_median=False,
                photfilename=None,
                load_photcat_if_exists=False,
                rematch_refcat=False,
                SNR_min = None, # minimum S/N for photometry
                # find best matches to refcut
                d2d_max = None, # maximum distance refcat to source in image
                dmag_max =None, # maximum uncertainty of source 
                sharpness_lim = (None, None), # sharpness limits
                roundness1_lim = (None, None), # roundness1 limits 
                delta_mag_lim = (None, None), # limits on mag-refcat_mainfilter
                objmag_lim = (None,None), # limits on mag, the magnitude of the objects in the image
                refmag_lim = (None,None), # limits on refcat_mainfilter, the magnitude of the reference catalog                    
                Nbright4match=None, # Use only the the brightest  Nbright sources from image for the matching with the ref catalog
                Nbright=None,    # Use only the brightest Nbright sources from image after all other cuts
                histocut_order='dxdy', # histocut_order defines whether the histogram cut is first done for dx or first for dy
                slope_min=-10/2048.0, 
                slope_Nsteps = 200, # slope_max=-slope_min, slope_stepsize=(slope_max-slope_min)/slope_Nsteps
                Nfwhm = 2.5,
                xshift=0.0,# added to the x coordinate before calculating ra,dec. This can be used to correct for large shifts before matching!
                yshift=0.0, # added to the y coordinate before calculating ra,dec. This can be used to correct for large shifts before matching!
                iterate_with_xyshifts = False,
                showplots=0,
                saveplots=0,
                savephottable=0,
                ee_radius=70,
                check_only=False,
                photometry_method='aperture'                ):
            
        # set self.outbasename based on option
        self.set_outbasename(outrootdir=outrootdir,outsubdir=outsubdir,inputname=input_image)
        
        # set the telescope
        self.set_telescope(telescope=telescope,imname=input_image)
            
        # do the photometry
        self.phot.verbose = self.verbose
        photfilename = f'{self.outbasename}.phot.txt'
        (photfilename_4check,photcat_loaded) = self.phot.run_phot(input_image,
                                                                  use_dq=use_dq,
                                                                  photfilename=photfilename,
                                                                  load_photcat_if_exists=load_photcat_if_exists,
                                                                  overwrite=overwrite,
                                                                  Nbright4match=Nbright4match,
                                                                  SNR_min=SNR_min,
                                                                  xshift=xshift,
                                                                  yshift=yshift,
                                                                  ee_radius=ee_radius,
                                                                  photometry_method=photometry_method)
        if (photfilename!=photfilename_4check):
            raise RuntimeError(f'BUG!!! {photfilename}!={photfilename_4check}')
        
        # make the initial cut on the image photometry catalog on magnitudes, sharpness, roundness etc
        ixs_use = self.phot.initial_cut_photcat(dmag_max = dmag_max,
                                                sharpness_lim = sharpness_lim, # sharpness limits
                                                roundness1_lim = roundness1_lim, # roundness1 limits 
                                                objmag_lim = objmag_lim, # limits on mag, the magnitude of the objects in the image
                                                Nbright = Nbright)
        
        # initialize an (existing!) matched refcat only, instead of loading
        # and matching the refcat for the following reasons:
        # if the photcat is loaded (assuming this means that the loaded
        # photcat has already a matched refcat), and if rematch_refcat=False
        initialize_refcat_only = (not rematch_refcat) and photcat_loaded
        # load the refcat and match it
        self.phot.load_and_match_refcat(ixs_obj=ixs_use,
                                        refcatname=refcatname,
                                        refcat_racol=refcat_racol,
                                        refcat_deccol=refcat_deccol,
                                        refcat_magcol=refcat_magcol,
                                        refcat_magerrcol=refcat_magerrcol,
                                        refcat_colorcol=refcat_colorcol,
                                        refmag_lim=refmag_lim, # limits for initial cut
                                        refmagerr_lim=(None,None), # limits for initial cut, needs to be added to options
                                        refcolor_lim=(None,None), # limits for initial cut, needs to be added to options
                                        pmflag = pmflag, 
                                        pm_median=pm_median,
                                        initialize_only=initialize_refcat_only
                                        )

        ixs_bestmatch= self.find_good_refcat_matches(ixs=ixs_use,
                                                     d2d_max = d2d_max,
                                                     delta_mag_lim = delta_mag_lim, # limits on mag-refcat_mainfilter
                                                     refmag_lim = refmag_lim, # limits on refcat_mainfilter, the magnitude of the reference catalog
                                                     histocut_order=histocut_order,
                                                     slope_min=slope_min, 
                                                     slope_Nsteps = slope_Nsteps, # slope_max=-slope_min, slope_stepsize=(slope_max-slope_min)/slope_Nsteps
                                                     Nfwhm = Nfwhm,
                                                     show_initial_plot=showplots,
                                                     show_histofit_plots=showplots,
                                                     savephottable=savephottable,
                                                     outbasename=self.outbasename
                                                     )     
        
        # If iterate with x/yshift
        if iterate_with_xyshifts:
            # get the median of dx and dy of the best matched objects from the first iteration!
            dx_median = self.phot.t.loc[ixs_bestmatch,'dx'].median()
            dy_median = self.phot.t.loc[ixs_bestmatch,'dy'].median()
            print(f'dx median of best matched objects of 1st iteration: {dx_median}',
                  f'dy median of best matched objects of 1st iteration: {dy_median}')
            
            # remove all columns from previous iteration in the main table
            refcatcols = ['dx','dy','delta_mag']
            for col in self.phot.t.columns:
                if re.search(f'^{self.phot.refcat.short}_',col):
                    refcatcols.append(col)
            self.phot.t.drop(refcatcols,axis=1)
    
            # re-calculate the ra,dec, now with the xy shifts from the first iteration
            self.phot.xy_to_radec(indices=ixs_use,xshift=dx_median,yshift=dy_median)
            
            # rematch the refcat 
            self.phot.match_refcat(ixs_obj=ixs_use,
                                   ixs_refcat=self.phot.ixs_use_refcat)

            # find again the best matches
            ixs_bestmatch= self.find_good_refcat_matches(ixs=ixs_use,
                                                         d2d_max = d2d_max,
                                                         delta_mag_lim = delta_mag_lim, # limits on mag-refcat_mainfilter
                                                         refmag_lim = refmag_lim, # limits on refcat_mainfilter, the magnitude of the reference catalog
                                                         histocut_order=histocut_order,
                                                         slope_min=slope_min, 
                                                         slope_Nsteps = slope_Nsteps, # slope_max=-slope_min, slope_stepsize=(slope_max-slope_min)/slope_Nsteps
                                                         Nfwhm = Nfwhm,
                                                         show_initial_plot=0,
                                                         show_histofit_plots=showplots,
                                                         savephottable=savephottable,
                                                         outbasename=self.outbasename
                                                         )     

        if check_only:
            dxdy_plot(self.phot, ixs_bestmatch,title='pre WCS correction',
                      refcat_mainfilter=self.phot.refcat_mainfilter,
                      refcat_mainfilter_err=self.phot.refcat_mainfilter_err,
                      refcat_maincolor=self.phot.refcat_maincolor
                     )
            if saveplots:
                outfilename = f'{outbasename}.phot.prewcs.png'
                if os.path.isfile(outfilename):
                    rmfile(outfilename)
                if self.verbose:
                    print(f'Saving {outfilename}')
                plt.savefig(outfilename)
            if showplots>0:
                plt.show()
            plt.close()
            return(0)
        jhatfits = f'{self.outbasename}_jhat.fits'
        (runflag,jhatfits) = self.run_align2refcat(input_image,
                                                   outputfits=jhatfits,
                                                   ixs=ixs_bestmatch,
                                                   overwrite=overwrite,skip_if_exists=skip_if_exists)
        #if self.telescope.lower()=='jwst':
        self.update_phottable_final_wcs(jhatfits,
                                            ixs_bestmatch = ixs_bestmatch,
                                            showplots=showplots,
                                            saveplots=saveplots,
                                            savephottable=savephottable,
                                            overwrite=overwrite
                                            )

        if showplots: 
            print('*** Note: close plots to continue!')
            plt.show()

        return(0)

