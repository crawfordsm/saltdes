
"""
SNREDUCE

Reduction script for SALT long slit data for 
follow up of DES Supernova candidates

This includes step that are not yet included in the pipeline 
and can be used for extended reductions of SALT data. 

It does require the pysalt package to be installed 
and up to date.

"""

import os, sys, glob, shutil

import numpy as np
import pyfits
from scipy.ndimage.filters import median_filter


from pyraf import iraf
from iraf import pysalt

from saltobslog import obslog
from saltprepare import  *
from saltbias import saltbias
from saltgain import saltgain
from saltxtalk import saltxtalk
from saltcrclean import saltcrclean
from saltcombine import saltcombine
from saltflat import saltflat
from saltmosaic import saltmosaic
from saltillum import saltillum

from specidentify import specidentify
from specrectify import specrectify
from specsky import skysubtract
from specextract import extract, write_extract
from specsens import specsens
from speccal import speccal

from PySpectrograph.Spectra import findobj

def science_red(rawdir, prodir, imreduce=True, specreduce=True, bpmfile=None, calfile=None, lampfile='Ar', automethod='Matchlines', skysection=[800,1000], cleanup=True):
    print rawdir
    print prodir

    #get the name of the files
    infile_list=glob.glob(rawdir+'P*.fits')
    infiles=','.join(['%s' % x for x in infile_list])
    

    #get the current date for the files
    obsdate=os.path.basename(infile_list[0])[1:9]
    print obsdate

    #set up some files that will be needed
    logfile='spec'+obsdate+'.log'
    flatimage='FLAT%s.fits' % (obsdate)
    dbfile='spec%s.db' % obsdate

    #create the observation log
    obs_dict=obslog(infile_list)

 
    if imreduce:   
      #prepare the data
      saltprepare(infiles, '', 'p', createvar=False, badpixelimage='', clobber=True, logfile=logfile, verbose=True)

      #bias subtract the data
      saltbias('pP*fits', '', 'b', subover=True, trim=True, subbias=False, masterbias='',  
              median=False, function='polynomial', order=5, rej_lo=3.0, rej_hi=5.0, 
              niter=10, plotover=False, turbo=False, 
              clobber=True, logfile=logfile, verbose=True)

      add_variance('bpP*fits', bpmfile)

      #gain correct the data
      saltgain('bpP*fits', '', 'g', usedb=False, mult=True, clobber=True, logfile=logfile, verbose=True)

      #cross talk correct the data
      saltxtalk('gbpP*fits', '', 'x', xtalkfile = "", usedb=False, clobber=True, logfile=logfile, verbose=True)

 
      #flat field correct the data
      flat_imgs=''
      for i in range(len(infile_list)):
        if obs_dict['CCDTYPE'][i].count('FLAT'):
           if flat_imgs: flat_imgs += ','
           flat_imgs += 'xgbp'+os.path.basename(infile_list[i])

      if 0: #len(flat_imgs)!=0:
         saltcombine(flat_imgs,flatimage, method='median', reject=None, mask=False,    \
                weight=False, blank=0, scale=None, statsec='[200:300, 600:800]', lthresh=3,    \
                hthresh=3, clobber=True, logfile=logfile, verbose=True)
         saltillum(flatimage, flatimage, '', mbox=11, clobber=True, logfile=logfile, verbose=True)

         saltflat('xgbpP*fits', '', 'f', flatimage, minflat=0.8, allext=False, clobber=True, logfile=logfile, verbose=True)
      else:
         flats=None
         imfiles=glob.glob('xgbpP*fits')
         for f in imfiles:
             shutil.copy(f, 'f'+f)

      #cosmic ray clean the data
      #only clean the object data
      for i in range(len(infile_list)):
        if obs_dict['CCDTYPE'][i].count('OBJECT') and obs_dict['INSTRUME'][i].count('RSS'):
          img='fxgbp'+os.path.basename(infile_list[i])
          saltcrclean(img, img, '', crtype='edge', thresh=5, mbox=11, bthresh=5.0,
                flux_ratio=0.2, bbox=25, gain=1.0, rdnoise=5.0, fthresh=5.0, bfactor=2,
                gbox=3, maxiter=5, multithread=True,  clobber=True, logfile=logfile, verbose=True)

      #mosaic the data
      geomfile=iraf.osfn("pysalt$data/rss/RSSgeom.dat")
      saltmosaic('fxgbpP*fits', '', 'm', geomfile, interp='linear', cleanup=True, geotran=True, clobber=True, logfile=logfile, verbose=True)

      #clean up the images
      if cleanup:
           for f in glob.glob('p*fits'): os.remove(f)
           for f in glob.glob('bp*fits'): os.remove(f)
           for f in glob.glob('gbp*fits'): os.remove(f)
           for f in glob.glob('xgbp*fits'): os.remove(f)
           for f in glob.glob('fxgbp*fits'): os.remove(f)


    #set up the name of the images
    if specreduce:
       for i in range(len(infile_list)):
           if obs_dict['OBJECT'][i].upper().strip()=='ARC':
               lamp=obs_dict['LAMPID'][i].strip().replace(' ', '')
               arcimage='mfxgbp'+os.path.basename(infile_list[i])
               #lampfile=iraf.osfn("pysalt$data/linelists/%s.salt" % lamp)

               specidentify(arcimage, lampfile, dbfile, guesstype='rss', 
                  guessfile='', automethod=automethod,  function='legendre',  order=3, 
                  rstep=100, rstart='middlerow', mdiff=50, thresh=2, niter=5, smooth=3,
                  inter=True, clobber=True, logfile=logfile, verbose=True)

               specrectify(arcimage, outimages='', outpref='x', solfile=dbfile, caltype='line', 
                   function='legendre',  order=3, inttype='interp', w1=None, w2=None, dw=None, nw=None,
                   blank=0.0, clobber=True, logfile=logfile, verbose=True)
     

    objimages=''
    for i in range(len(infile_list)):
       if obs_dict['CCDTYPE'][i].count('OBJECT') and obs_dict['INSTRUME'][i].count('RSS'):
          if objimages: objimages += ','
          objimages+='mfxgbp'+os.path.basename(infile_list[i])

    if specreduce:
      #run specidentify on the arc files

      specrectify(objimages, outimages='', outpref='x', solfile=dbfile, caltype='line', 
           function='legendre',  order=3, inttype='interp', w1=None, w2=None, dw=None, nw=None,
           blank=0.0, clobber=True, logfile=logfile, verbose=True)


    return

def add_variance(filenames, bpmfile):
    file_list=glob.glob(filenames)
    badpixelstruct = saltio.openfits(bpmfile)
    for f in file_list:
        print f
        struct = pyfits.open(f)
        nsciext=len(struct)-1
        nextend=nsciext
        for i in range(1, nsciext+1):
            hdu=CreateVariance(struct[i], i, nextend+i)
            struct[i].header.set('VAREXT',nextend+i, comment='Extension for Variance Frame')
            struct.append(hdu)
        nextend+=nsciext
        for i in range(1, nsciext+1):
            hdu=createbadpixel(struct, badpixelstruct, i, nextend+i)
            struct[i].header.set('BPMEXT',nextend+i, comment='Extension for Bad Pixel Mask')
            struct.append(hdu)
        nextend+=nsciext
        struct[0].header.set('NEXTEND', nextend)
        if os.path.isfile(f): os.remove(f)
        struct.writeto(f)

if __name__=='__main__':
   rawdir=sys.argv[1]
   prodir=os.path.curdir+'/'
   lampfile = os.path.dirname(sys.argv[0]) + '/Ar.sn'
   bpmfile = os.path.dirname(sys.argv[0]) + '/bpm_sn.fits'
   science_red(rawdir, prodir, imreduce=True, specreduce=True, cleanup=True, bpmfile=bpmfile, lampfile=lampfile)
