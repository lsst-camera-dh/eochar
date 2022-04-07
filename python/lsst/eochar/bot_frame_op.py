#!/usr/bin/env python
# coding: utf-8

###################################################################################################
#
# File : frame_study.py
#
# Auhtor : P.Antilogus
#
# Version : 30 March 2022
#
# Goal : this python file read raw data image , and can be used for specific sensor diagnostic like :
#          - cte
#          - overscan
#          - noise
#
# Example : see notebooks using this package
#
# Remark : it is under the process to be cleaned - simplified ...but this is not my top priority today ;-) 
#


import astropy.io.fits as pyfits
import numpy as np
import glob
import os
import sys
import matplotlib.pyplot as plt
import matplotlib
import time
import pickle
matplotlib.rcParams['axes.formatter.useoffset'] = False

#
def image_area(image) :
    # input : image ==> fits image 
    # output :  image section coordinate to be used in python table:  ymin , ymax ,xmin, xmax
    # -use pyfits to open the file 
    # -extract the image area to be used in python table  from the DATASEC keyword
    #
    r=image[1].header['DATASEC'][1:-1].split(',')
    x=r[0].split(':')
    y=r[1].split(':')
    #
    return int(y[0])-1,int(y[1]),int(x[0])-1,int(x[1])

def smooth253(table) :
    # input 1D table to be smoothed  according to algorithm 353QH, TWICE (see HBOOK)
    # the 4 first and last pixels are left unchanged : this is not a bad idea at least for
    # the first 4 as we can have a high frequecy variation on the bias of the first lines 
    sum353=np.copy(table)
    for k in [3,5,3] :
        sumk=np.copy(sum353[k-1:])
        for i in range(1,k) :
            sumk+=sum353[k-1-i:-i]
        sum353=sumk/k
    result=np.copy(table)
    result[4:-4]=sum353
    #print('smooth353 called with a table of ',len(result),' raws , std table =',np.std(table),' std result =',np.std(result))
    return result

#
class InFile :
   # Handel ( select et all ) file list from LSST ccd test bench 
    def __init__(self,dirall=['/Users/antilog/scratch/20160901'],Pickle=False,root_for_pickle='/sps/lsst/DataBE/lpnws5203',fkey={},verbose=False,Slow=True,Bias='2D',single_t=False,nskip=0,nkeep=-1):
        # dirall : a list of directory/file   to read in  : the file header will be used to select or not the file if fkey is set  or the file will be read from the content of the pickle file (if Pickle is True ) , fkey will also be used for selection. 
        # fkey : {'selection_name' : {'fits_header_name':{'key':value} , ... }  : a triple  dictionary of { 'selection_name' : {header : {key:value}} } , to select a file 
        self.nkept=0
        self.all_file=[]
        self.clap=[]
        self.selection=[]
        self.fkey=fkey
        self.directory=sorted(dirall)
        self.stat={}
        self.nseen=0
        # loop on header 
        if Pickle : 
            self.all_file_from_pickle(dirall=dirall,root_for_pickle=root_for_pickle,fkey=fkey,verbose=verbose,Slow=Slow,Bias=Bias,single_t=single_t,nskip=nskip,nkeep=nkeep)
        else : 
            self.all_file_from_dir(dirall=dirall,fkey=fkey,verbose=verbose,Slow=Slow,Bias=Bias,single_t=single_t,nskip=nskip,nkeep=nkeep)
        return
    def all_file_from_dir(self,dirall,fkey,verbose,Slow,Bias,single_t,nskip,nkeep):
       # dirname : can be a directory name or a file with * at the moment it ends by .fz                                                                        
        #  ex : /Users/antilog/scratch/REB_DATA/20160513/linearity                                                                                               
        #  or : /Users/antilog/scratch/REB_DATA/20160513/linearity/reb3*.fz                                                                                      
        # fkey : dictionary key word for image selection
        # verbose :  if True , print messages about each selected file ( default= False )
        # single_t :  keep only one file per exposure time (default False ) 
        self.all_file=[]
        old_time=[0.]
        fits_is_open=False
        #
        # fill all_file from the header of the files in dirall . 
        for dirname in self.directory :
            # have we allready all the needed file ? 
            if nkeep> 0 and self.nkept==nkeep : 
                break  
            # build the list of file for this directory 
            if (len(os.path.splitext(dirname)[1])>0) :
                file_list=glob.glob(dirname)
            else :
                file_list=glob.glob("%s/*.fz" % (dirname))
            file_list.sort()
            # loop on files to select them if needed  
            for filenamed  in file_list :
                #
                keep=True
                if len(fkey)==0 :  
                    selection='Main'
                else : 
                    fits_is_open=False
                    dname, fname=os.path.split((filenamed))
                    # is there any extra selection based on Header , key , value  ? 
                    for selection, sel_id in fkey.items() :
                        local_keep=False
                        # select of file name if any
                        if 'first' in sel_id.keys() :
                            if  fname < sel_id['first'] : continue
                        if 'last' in sel_id.keys() :
                            if  fname > sel_id['last'] : continue
                        local_keep=True
                        if 'key' in sel_id.keys() : 
                            if not(fits_is_open) :
                                fitsfile=pyfits.open(filenamed)
                                fits_is_open=True
                            for header, key_cur in sel_id['key'].items() :
                                if not ( local_keep ) : break  
                                for key, value in key_cur.items() :
                                    if ( key in fitsfile[header].header ) :
                                        if ( fitsfile[header].header[key]!=value ):
                                            local_keep=False
                                            break 
                                    else : 
                                        local_keep=False
                                        break
                        # this file is ok for the current selection ,  remark : a file is selected in the first selection that is compatible with 
                        if local_keep : break 
                    keep=local_keep
                #
                if (keep and single_t ) :
                    if not(fits_is_open) :
                        fitsfile=pyfits.open(filenamed)
                        fits_is_open=True
                    new_time=fitsfile[0].header['EXPTIME']
                    if new_time in old_time : 
                        keep=False
                    else : 
                        old_time.append(new_time)
                if (keep) :
                    self.nseen+=1
                    if self.nseen>nskip :
                        if not(fits_is_open) :
                            fitsfile=pyfits.open(filenamed)
                            fits_is_open=True
                        self.all_file.append(actfile(fitsfile,Slow,Bias))
                        self.selection.append((selection,filenamed))
                        # to be updated with a call to clap
                        #self.clap.append(new_time)
                        if verbose : print ('%d : Selected %s File %s ' % (self.nkept,selection,filenamed) )
                        self.nkept+=1
                        if self.nkept==nkeep and nkeep > 0  :
                             # we selected the number of files requested 
                            fitsfile.close()
                            break
                if fits_is_open :
                    fitsfile.close()
                    del fitsfile
                    fits_is_open=False
        return
                
class actfile :
    def __init__(self, fitsfile,Slow=True,Bias='2D'):
        #
        # Input Param :
        #    fitsfile : pyfits object with the image to process
        #    Slow      :  True  : return : <> , var , nb_pix used , for each amplifier after bias subtraction  , but don't return the image
        #                 False : return the image ( in Self.Image[amp][:,:] ) for each amplifier after Bias correction , but not the <> , var,...
        #    Bias      : How the Bias correct is done :
        #                 '2D' (default)  : 2D bias correction using the  overscan line and column
        #                 '1D'            : 1D bias correction using the overscan per line 
        #                 'Ct'            : 1 global bias number computed from all serial overscan
        #                 '2D y353' or '1D y353' or '2D x353' or '2D x353 y353'  :  x : smooth <// overscan> along colomn , y : smooth <serial overscan> along line using the 353 smoothing algo before computing the Bias correction       
        #                                 
        # next line need care , as the memory will grow quickly
        # self.fits=np.copy(fitsfile)
        self.Image=[]
        self.OverCol=[]
        self.OverLine=[]
        self.Mean=[]
        self.Var=[]
        self.Nb=[]
        self.Channel=[]
        self.Extname=[]
        self.Datasec=[]
        self.Detsec=[]
        self.Detsize=[]
        # image area
        first_line,first_p_over,first_col,first_s_over=image_area(fitsfile)
        self.first_col=first_col
        self.first_s_over=first_s_over
        self.first_line=first_line
        self.first_p_over=first_p_over
        # 
        try :
            self.exptime=float(fitsfile[0].header['EXPTIME'])
        except :
            try:
                self.exptime=float(fitsfile[0].header['EXPOSURE'])
            except:
                self.exptime=None
        try :
           self.obsid=(fitsfile[0].header['OBSID']).strip()
           self.ccdslot=(fitsfile[0].header['CCDSLOT']).strip()
           self.raftbay=(fitsfile[0].header['RAFTBAY']).strip()
           self.lsst_num=(fitsfile[0].header['LSST_NUM']).strip()
           self.mjd=float(fitsfile[0].header['MJD-TRG'])   
        except :
           self.obsid=''
           self.ccdslot=''
           self.raftbay=''
           self.lsst_num=''
           self.mjd=0.
        # which Bias correction ?
        Bias_val=Bias.split(' ')
        Bias_cor='2D'
        xsmooth=''
        ysmooth=''
        #
        if len(Bias_val) > 0 :
            Bias_cor=Bias_val[0]
            for i in range(1,len(Bias_val)) :
                if 'x' == Bias_val[i][0] :
                    xsmooth=Bias_val[i][1:]
                elif  'y' == Bias_val[i][0] :
                    ysmooth=Bias_val[i][1:]
        #       
#                                                                                                              
# Number of boxes to compute var and mean for each amps 
        nstepy=4
        nstepx=4
# self.Date=JdUtc(fitsfile[0].header['DATE']).Jd for i in range(len(fitsfile)):                                
        for i in range(1,min(17,len(fitsfile))):
            if (fitsfile[i].header['XTENSION'].strip()=='IMAGE' ) :
                self.Channel.append(i-1)
                self.Extname.append((fitsfile[i].header['EXTNAME']).strip())
                self.Datasec.append((fitsfile[i].header['DATASEC']).strip())
                self.Detsec.append((fitsfile[i].header['DETSEC']).strip())
                try :
                    self.Detsize.append((fitsfile[i].header['DETSIZE']).strip())
                except :
                    self.Detsize.append((fitsfile[0].header['DETSIZE']).strip())
                # size
                last_l=len(fitsfile[i].data[:,0])
                last_s=len(fitsfile[i].data[0,:])
                # serial overscan
                mean_over_per_line=np.mean(fitsfile[i].data[:,first_s_over+2:],axis=1)
                rawl=np.zeros((last_l-first_p_over-2,last_s))
                # // ovesrcan per column , corrected by the serial value per line
                for l in range(first_p_over+2,last_l):
                    rawl[l-first_p_over-2,:]=fitsfile[i].data[l,:]-mean_over_per_line[l]
                # from // overscan  correction per colomn 
                mean_over_per_column=np.mean(rawl[:,:],axis=0)
                if xsmooth == '353' :
                     mean_over_per_column=smooth253(mean_over_per_column)
                # from serial overscan correction per line
                if ysmooth == '353' :
                     mean_over_per_line=smooth253(mean_over_per_line)
                #
                linef=mean_over_per_line[:,np.newaxis]
                # generate the 2D correction (thank's to numpy) 
                #  over_cor_mean=mean_over_per_column+linef
                if not(Slow) :
                    if Bias_cor=='Ct' :
                        self.Image.append(np.single(fitsfile[i].data-mean_over_per_line.mean()))
                    elif Bias_cor=='1D' :
                        self.Image.append(np.single(fitsfile[i].data-(mean_over_per_column*0.+linef)))
                    else :
                        # if 2D or what ever else 
                        self.Image.append(np.single(fitsfile[i].data-(mean_over_per_column+linef)))
                    self.OverCol.append(mean_over_per_column)
                    self.OverLine.append(mean_over_per_line)
                else :
                    if Bias_cor=='Ct' :
                        Image_single=fitsfile[i].data-mean_over_per_line.mean()
                    elif Bias_cor=='1D' :
                        Image_single=fitsfile[i].data-(mean_over_per_column*0.+linef)
                    else : 
                        Image_single=fitsfile[i].data-(mean_over_per_column+linef)
                    IMean=np.zeros((nstepy,nstepx))
                    IVar=np.zeros((nstepy,nstepx))
                    INb=np.zeros((nstepy,nstepx))
                    ystart=first_line+10
                    ystep=int((first_p_over-first_line-20)/nstepy)
                    xstart=first_col+10
                    xstep=int((first_s_over-first_col-20)/nstepx)
                    #
                    for iy in range(nstepy) :
                        for ix in range(nstepx) :
                            box=Image_single[ystart+iy*ystep:ystart+(iy+1)*ystep,xstart+ix*xstep:xstart+(ix+1)*xstep]
                            mbox=box.mean()
                            sbox=box.std()
                            #box_trunc=np.array([ box[i,j] for i in range(ystep) for j in range(xstep) if abs(box[i,j]-mbox)<5*sbox])
                            #mbox=box_trunc.mean()
                            #sbox=box_trunc.std()
                            box_trunc=np.array([ box[i,j] for i in range(ystep) for j in range(xstep) if abs(box[i,j]-mbox)<4*sbox])
                            IMean[iy,ix]=box_trunc.mean()
                            IVar[iy,ix]=(box_trunc.std())**2
                            INb[iy,ix]=len(box_trunc)
                    self.Mean.append(IMean)
                    self.Var.append(IVar)
                    self.Nb.append(INb)      
        return
    def WriteImage(self,out_root) :
        # Write the multi-frame image un-biased 
        hdu=pyfits.PrimaryHDU()
        new_hdul = pyfits.HDUList([hdu])
        k=0
        new_hdul[k].header['EXPTIME']=self.exptime
        new_hdul[k].header['OBSID']=   self.obsid
        new_hdul[k].header['CCDSLOT']=   self.ccdslot
        new_hdul[k].header['RAFTBAY']=  self.raftbay
        new_hdul[k].header['LSST_NUM']=   self.lsst_num
        new_hdul[k].header['MJD-TRG']=   self.mjd

        for k in range(len(self.Image)) :
            new_hdul.append(pyfits.ImageHDU(self.Image[k]))
            new_hdul[k+1].header['EXTNAME']=(self.Extname[k])
            new_hdul[k+1].header['DATASEC']=    self.Datasec[k]
            new_hdul[k+1].header['DETSEC']=   self.Detsec[k]
            new_hdul[k+1].header['DETSIZE']=   self.Detsize[k]
            k+=1
        out_file=out_root+self.obsid+"_"+self.raftbay+"_"+self.ccdslot+"_unbiased.fits"
        new_hdul.writeto(out_file)
        new_hdul.close()
        del new_hdul

