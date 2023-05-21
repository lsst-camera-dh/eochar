#!/usr/bin/env python
# coding: utf-8

###################################################################################################
#
# File : frame_study.py
#
# Auhtor : P.Antilogus
#
# Version : 22 Feb 2019
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

try:
    import pyfits
except : 
    import astropy.io.fits as pyfits
import numpy as np
import glob
import os
import sys
import matplotlib.pyplot as plt
import matplotlib
import copy
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

#
class Ifile :
   # Handel ( select et all ) file list from LSST ccd test bench 
    def __init__(self,dirall=['/Users/antilog/scratch/20160901'],Pickle=False,root_for_pickle='/sps/lsst/DataBE/lpnws5203',fkey={},verbose=False,Slow=True,single_t=False,nskip=0,nkeep=-1):
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
            self.all_file_from_pickle(dirall=dirall,root_for_pickle=root_for_pickle,fkey=fkey,verbose=verbose,Slow=Slow,single_t=single_t,nskip=nskip,nkeep=nkeep)
        else : 
            self.all_file_from_dir(dirall=dirall,fkey=fkey,verbose=verbose,Slow=Slow,single_t=single_t,nskip=nskip,nkeep=nkeep)
        return
    def all_file_from_dir(self,dirall,fkey,verbose,Slow,single_t,nskip,nkeep):
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
                        self.all_file.append(datafile(fitsfile,Slow))
                        self.selection.append((selection,filenamed))
                        # to be updated with a call to clap
                        #self.clap.append(new_time)
                        if verbose : print ('%d : Selected %s File %s                                  ' % (self.nkept,selection,filenamed),end='\r' )
                        self.nkept+=1
                        if self.nkept==nkeep and nkeep > 0  :
                             # we selected the number of files requested 
                            fitsfile.close()
                            break
                if fits_is_open :
                    fitsfile.close()
                    del fitsfile
                    fits_is_open=False
        if verbose : print ('\n')
        return
    def all_file_from_pickle(self,dirall,root_for_pickle,fkey,verbose,Slow,single_t,nskip,nkeep):
        old_time=[0.]
         # fill all_file from the header of the files in dirall . 
        for dirname in self.directory :
            # have we allready all the needed file ? 
            if nkeep> 0 and self.nkept==nkeep : 
                break  
            # build the list of file for this directory 
            if (len(os.path.splitext(dirname)[1])>0) :
                file_list=glob.glob(dirname)
            else :
                file_list=glob.glob("%s/*.pkl" % (dirname))
            file_list.sort()
            # loop on files to select them if needed  
            for pickle_file  in file_list :
                # open the pickle file
                input=open(pickle_file,'rb')
                file=pickle.load(input)
                #
                for i_cur in range(len(file)) : 
                    #filename=file[i_cur].filename  
                    dname=file[i_cur].dir  
                    fname=file[i_cur].filename
                    clap=file[i_cur].clap
                    keep=True
                    if len(fkey)==0 :  
                        selection='Main'
                    else :
                        #
                        # is there any extra selection based on Header , key , value  ? 
                        for selection, sel_id in fkey.items() :
                            local_keep=False
                            # select of file name if any
                            if 'first' in sel_id.keys() :
                                if  fname < sel_id['first'] : continue
                            if 'last' in sel_id.keys() :
                                 if  fname > sel_id['last'] : continue
                            local_keep=True
                            # test key (=)  key+ (=>) key- (<=)
                            if 'key' in sel_id.keys() : 
                                for header, key_cur in sel_id['key'].items() :
                                    if not ( local_keep ) : break  
                                    for key, value in key_cur.items() :
                                        if ( key in file[i_cur].header[header] ) :
                                            if (file[i_cur].header[header][key] is None) or ( file[i_cur].header[header][key]!=value ):
                                                local_keep=False
                                                break 
                                        else : 
                                            local_keep=False
                                            break
                            # this file is ok for the current selection ,  remark : a file is selected in the first selection that is compatible with 
                            if not(local_keep) : continue 
                            #
                            if 'key+' in sel_id.keys() : 
                                for header, key_cur in sel_id['key+'].items() :
                                    if not ( local_keep ) : break  
                                    for key, value in key_cur.items() :
                                        if ( key in file[i_cur].header[header] ) :
                                            if (file[i_cur].header[header][key] is None) or ( file[i_cur].header[header][key]<value ):
                                                local_keep=False
                                                break 
                                        else : 
                                            local_keep=False
                                            break
                            # this file is ok for the current selection ,  remark : a file is selected in the first selection that is compatible with 
                            if not(local_keep) : continue
                            #
                            if 'key-' in sel_id.keys() : 
                                for header, key_cur in sel_id['key-'].items() :
                                    if not ( local_keep ) : break  
                                    for key, value in key_cur.items() :
                                        if ( key in file[i_cur].header[header] ) :
                                            if (file[i_cur].header[header][key] is None) or ( file[i_cur].header[header][key]>value ):
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
                        new_time=file[i_cur].header['Primary']['EXPTIME']
                        if new_time in old_time : 
                            keep=False
                        else : 
                            old_time.append(new_time)
                    if (keep) :
                        self.nseen+=1
                        if self.nseen>nskip :
                            fitsfile=pyfits.open(root_for_pickle+'/'+dname+'/'+fname)
                            self.all_file.append(datafile(fitsfile,Slow))
                            fitsfile.close()
                            #
                            self.clap.append(clap)
                            self.selection.append(selection)
                            if verbose : print ('%d : Selected %s File %s ' % (self.nkept,selection,fname) )
                            self.nkept+=1
                            if self.nkept==nkeep and nkeep > 0  :
                                # we selected the number of files requested 
                                break
        return
    def plot(self,plt_name='',on_screen=False) :
        # define last to plot :
        #
        fig=plt.figure(figsize=(15,15))
        title="Noise estimation from Overscan : %s " % (plt_name)
        fig.suptitle(title)
        iplt=1
        ax=fig.add_subplot(3,3,iplt)
        iplt+=1
        
        return
                
class datafile :
    def __init__(self, fitsfile,Slow=True):
        '''
           Construct all the necessary attributes for the datafile object.
  
           Parameters : 
                 fitsfile (list of str ) :  list of file to process , they should be all from the same raft-sensor
                 Slow     (boll)         :  computed extended image properties or not ( Default :  True ) 
                                            remark for CTE analysis alone , Slow can be set to False , will be faster 
        '''
        #
        self.Image=[]
        self.Hdu=[]
        self.HduMax=0
        self.fft=[]
        self.w=[]
        self.Mean=[]
        self.Std=[]
        self.Median=[]
        self.MedPScan=[]
        self.StdPScan=[]
        self.MeanSScan=[]
        self.MedSScan=[]
        self.StdSScan=[]
        self.StdSScanOS=[]
        self.MeanPScan=[]
        self.StdPScanOS=[]
        self.Mean_col=[]
        self.Std_col=[]
        self.Mean_line=[]
        self.Median_line=[]
        self.Std_line=[]
        self.Std_l60=[]
        self.Attenuator=0.
        self.CCD_COND={}
        self.Range=0.
        self.PreExp=0.
        self.PostExp=0.
        # image area
        first_line,first_p_over,first_col,first_s_over=image_area(fitsfile)
        self.first_col=first_col
        self.first_s_over=first_s_over
        self.first_line=first_line
        self.first_p_over=first_p_over
        #
        try :
            # the smallest modification for this code to take into account the CCOB is to replace EXPTIME
            # by CCOBFLUX ... hum sound dangerous ... but we don't have much time to think today ...
            self.exptime=float(fitsfile[0].header['CCOBFLUX'])
        except :
            try:
                self.exptime=float(fitsfile[0].header['EXPTIME'])
            except :
                # Paris test bench key value for exposure time is different
                self.exptime=float(fitsfile[0].header['EXPOSURE'])
        try :
           self.ccdslot=(fitsfile[0].header['CCDSLOT']).strip()
           self.raftbay=(fitsfile[0].header['RAFTBAY']).strip()
        except :
           self.ccdslot=''
           self.raftbay=''
        self.fluxs_last=[]
        self.fluxp_last=[]
        self.fluxs_last_std=[]
        #self.fluxs_last_var=[]
        self.fluxp_last_std=[]
        self.fluxp_used=[]
#
        self.over4_col_std=[]
        self.over4_line_std=[]


#                                                                                                              
#                                                                                                              
# self.Date=JdUtc(fitsfile[0].header['DATE']).Jd for i in range(len(fitsfile)):                                
        for i in range(1,min(17,len(fitsfile))):
            if ( fitsfile[i].header['XTENSION'].strip()=='IMAGE' ) :
                self.Hdu.append(i)
                self.HduMax=i
                # Remark : for the moment we don't which REB slice we are looki 
                #  for e2v and BNL  data it's [8:]
                # self.Channel.append(int(fitsfile[i].header['EXTNAME'][8:]))
                #  for Paris   data it's [5:]
                #self.Channel.append(int(fitsfile[i].header['EXTNAME'][5:]))
                # self.Image.append(np.copy(fitsfile[i].data))                                                 
                self.Image.append(fitsfile[i].header['EXTNAME'].strip())
                # image                                                                                        
                # Mean and noise 
                self.Median.append(np.median(fitsfile[i].data[first_line:first_p_over,first_col:first_s_over])) 
                if Slow :
                    self.Mean.append(fitsfile[i].data[first_line:first_p_over,first_col:first_s_over].mean())   
                    self.Std.append(fitsfile[i].data[first_line:first_p_over,first_col:first_s_over].std())
                    # line OverScan                                                                                      
                    self.MedPScan.append(np.median(fitsfile[i].data[first_p_over+5:,first_col:first_s_over]))
                    self.MeanPScan.append(fitsfile[i].data[first_p_over+5:,first_col:first_s_over].mean())
                    self.StdPScan.append(fitsfile[i].data[first_p_over+5:,first_col:first_s_over].std())

                    # Serial over-scan   + 1 to remove CTE / xxshoot                                                                                  
                    self.MedSScan.append(np.median(fitsfile[i].data[:,first_s_over+5:]))
                    self.MeanSScan.append(fitsfile[i].data[:,first_s_over+5:].mean())
                    self.StdSScan.append(fitsfile[i].data[:,first_s_over+5:].std())

                # information for 2D diagmostic of the overscan : does the overscan is flat in function of the column ? line ?
                # --- data in the overscan corner ( pixels are overscan in line and column ) 
                self.over4_col_std.append(fitsfile[i].data[first_p_over:,first_s_over:].std(axis=1).mean())
                self.over4_line_std.append(fitsfile[i].data[first_p_over:,first_s_over:].std(axis=0).mean())

                if Slow :
                    # Same but bias subtrracted
                    mean_line=np.median(fitsfile[i].data[first_p_over:,:],axis=0)
                    mean_column=np.median(fitsfile[i].data[:,first_s_over:],axis=1)
                    last_l=len(fitsfile[i].data[:,0])
                    last_s=len(fitsfile[i].data[0,:])
                    rawl=np.zeros((last_l-first_p_over,last_s))
                    raws=np.zeros((last_l,last_s-first_s_over))
                    for l in range(first_p_over,last_l) :
                        rawl[l-first_p_over,:]=fitsfile[i].data[l,:]-mean_line
                    self.StdPScanOS.append((rawl[:,first_col:].std(axis=1)).mean())
                    #
                    for c in range(first_s_over,last_s) :
                        raws[:,c-first_s_over]=fitsfile[i].data[:,c]-mean_column
                    self.StdSScanOS.append((raws[first_line:,:].std(axis=0)).mean())
            

                
                # average allong the column and line                                                           
                #self.Mean_col.append(fitsfile[i].data[first_line:first_p_over,:].mean(axis=0))
                #self.Mean_line.append(fitsfile[i].data[:,first_col:first_s_over].mean(axis=1))
                #self.Median_line.append(np.median(fitsfile[i].data[:,first_col:first_s_over],axis=1))
                #self.Std_col.append(fitsfile[i].data[first_line:first_p_over,:].std(axis=0))
                #self.Std_line.append(fitsfile[i].data[:,first_col:first_s_over].std(axis=1))
                #self.Std_l60.append(fitsfile[i].data[:60,first_col:first_s_over].std())
                # For CTE Serie
                #
                # REMARK : The size cut ( 28 , line +/- 10 ... ) are hard wired and used in CTE part of the code to compute statistic !!!!
                #
                self.fluxs_last.append(fitsfile[i].data[first_line+10:first_p_over-10,first_s_over-1:first_s_over+28].mean(axis=0)-fitsfile[i].data[first_line+10:first_p_over-10,first_s_over+15:].mean())
                #
                self.fluxs_last_std.append(fitsfile[i].data[first_line+10:first_p_over-10,first_s_over-1:first_s_over+28].std(axis=0)/np.sqrt(float(first_p_over-10-first_line-10)))
                #
                #self.fluxs_last_var.append(fitsfile[i].data[first_line+100:first_p_over-100,first_s_over-1:first_s_over+28].std(axis=0)**2-fitsfile[i].data[first_line+100:first_p_over-100,first_s_over+15:].std()**2)
                # For CTE //
                #
                # fluxp=np.array([ fitsfile[i].data[first_p_over-1:first_p_over+28,icol] - np.median(fitsfile[i].data[first_p_over+5:,icol ])  for icol in range(first_col+10,first_s_over-10)  ])

                overscan_offset=  np.median(fitsfile[i].data[first_p_over+5:,first_col+10:first_s_over-10],axis=0)
                
                #
                # we do a median per slice , to kill outlier , but keep a statistical precision < 1 adu ...still not that precise compared to a mean 
                #self.fluxp_last.append((np.median(fluxp[0:100],axis=0)+np.median(fluxp[100:200],axis=0)+np.median(fluxp[200:300],axis=0)+np.median(fluxp[300:400],axis=0)+np.median(fluxp[400:500],axis=0))/5.)
                # the correct version : kill outlier  ( there is outlier in case of blooming column ) : to speed up we just kill based on the last physical column
                self.fluxp_last.append(np.zeros((29)))
                self.fluxp_last_std.append(np.zeros((29)))
                self.fluxp_used.append(np.zeros((29)))
                #
                last_line=fitsfile[i].data[first_p_over-1,first_col+10:first_s_over-10]-overscan_offset
                last_line_median=np.median(last_line)
                last_line_std=5*last_line.std()
                
                column_ok=[icol for  icol in range(len(last_line))  if np.abs(last_line[icol]-last_line_median) < last_line_std ]

                #                
                
                for j in range(29):
                    fluxp_truncated=(fitsfile[i].data[first_p_over-1+j,first_col+10:first_s_over-10]-overscan_offset)[column_ok]
                    self.fluxp_last[-1][j]=np.mean(fluxp_truncated)
                    self.fluxp_last_std[-1][j]=np.std(fluxp_truncated)/np.sqrt(len(fluxp_truncated))
                    self.fluxp_used[-1][j]=len(fluxp_truncated)
                #fluxp=np.array([ fitsfile[i].data[first_p_over-1:first_p_over+28,icol] - np.median(fitsfile[i].data[first_p_over+5:,icol ])  for icol in range(first_col+10,first_s_over-10)  ])
                #self.fluxp_last.append(np.median(fluxp,axis=0))
                # self.fluxp_last.append(np.median(fitsfile[i].data[first_p_over-1:first_p_over+28,first_col+10:first_s_over-10],axis=1)-np.median(fitsfile[i].data[first_p_over-1:first_p_over+28,first_s_over+5:first_s_over+15]))
            else:
                # last image section read
                break
        return
class cte :
    def __init__(self, all_file, gain=[0.704650434205,0.68883578783,0.688459358774,0.696697494642,0.689209827484,0.696579402812,0.698973006751,0.689613072912,0.682880384357,0.696206655845,0.690349506621,0.691506176017,0.690763478766,0.689762341309,0.694801544092,0.850025229184 ],serie=True):
        #
        nb_f_max=len(all_file)
        #
        self.cte_flux=np.zeros((16,nb_f_max))
        self.cte_time=np.zeros((nb_f_max))
        self.cte_ftime=np.zeros((nb_f_max))
        self.cte_flux_s=np.zeros((16,nb_f_max))
        self.cte_y=np.zeros((16,nb_f_max,28))
        self.cte_y_s=np.zeros((16,nb_f_max,28))
        self.cte_y_std=np.zeros((16,nb_f_max,28))
        self.cte_y_s_std=np.zeros((16,nb_f_max,28))
        self.ylev=np.zeros((16,nb_f_max))
        self.ylev_std=np.zeros((16,nb_f_max))
        self.nb_file=np.zeros((nb_f_max))
        self.serie=serie
        self.cte_noise_s=np.zeros((16,nb_f_max))
        self.cte_noise_s_std=np.zeros((16,nb_f_max))
        self.overscan_std=np.zeros((16,nb_f_max))
        self.over8_18=np.zeros((16,nb_f_max))
        self.over8_18_std=np.zeros((16,nb_f_max))
        # pixel number in unit of CCD 
        self.i_f=0
        #
        if nb_f_max==0 : return
        #
        self.first_file=copy.deepcopy(all_file[0])
        #
        cte_noise_std=np.zeros((16,nb_f_max,28))
        #
        for f in all_file :
            im_flux=np.median(np.array(f.Median))
            if self.i_f>0 and f.exptime in self.cte_time[0:self.i_f] :
                all_cur=np.flatnonzero(self.cte_time[0:self.i_f] == f.exptime)
                # Attention could be that we have the same exposure time but not the same flux (extra filter)
                found_cur=-1
                for cur_cur in all_cur  :
                    ratio=max(self.cte_ftime[cur_cur]/im_flux,im_flux/self.cte_ftime[cur_cur])
                    if ratio<1.1 :
                        found_cur=cur_cur
                if found_cur > -1 :
                    i_cur=found_cur
                else:
                    i_cur=self.i_f
                    self.cte_ftime[i_cur]=im_flux
                    self.cte_time[i_cur]=f.exptime
                    self.i_f+=1
            else :
                i_cur=self.i_f
                self.cte_ftime[i_cur]=im_flux
                self.cte_time[i_cur]=f.exptime
                self.i_f+=1
            for ch in range(f.HduMax) :
                if serie :
                    # CTE serie
                    if ch==0 :
                         # print ('%s ,time %f, Channel %d,  flux %f (flux last col %f) , image %f ,  signal dispersion %f , scan serie %f , scan serie dispersion %f ' % (f.filename,f.exptime,ch,f.Mean[ch]-f.MeanSScan[ch],f.fluxs_last[ch][0],f.Mean[ch], f.Std[ch],f.MeanSScan[ch],f.StdSScan[ch]))
                        self.first=f.first_s_over
                        # what matter in the CTE def is how many transfer you did for the last column read , which is the size of the pre-scan + size of the image (in the past we subtracted the prescan which is an error) 
                        self.nb_pixel=f.first_s_over
                        self.nb_file[i_cur]+=1
                    flux_last=f.fluxs_last[ch][0]
                    #
                    self.cte_y[ch,i_cur,:]+=f.fluxs_last[ch][1:]
                    self.cte_y_std[ch,i_cur,:]+=(f.fluxs_last_std[ch][1:])**2
                    cte_noise_std[ch,i_cur,:]+=(f.fluxs_last_std[ch][1:])**2*(float(f.first_p_over-10-f.first_line-10))
                    self.overscan_std[ch,i_cur]+=(f.over4_col_std[ch])**2

                else :
                    # CTE //
                    if ch==0 :
                        # print ('%s ,time % f, Channel %d,  flux %f (flux last line %f) , image %f ,  signal dispersion %f , scan //  %f , scan // dispersion %f ' % (f.filename,f.exptime,ch,f.Mean[ch]-f.MedPScan[ch],f.fluxp_last[ch][0],f.Mean[ch], f.Std[ch],f.MedPScan[ch],f.StdPScan[ch]))
                        self.first=f.first_p_over
                        self.nb_pixel=f.first_p_over
                        self.nb_file[i_cur]+=1
                    flux_last=f.fluxp_last[ch][0]
                    #
                    self.cte_y[ch,i_cur,:]+=f.fluxp_last[ch][1:]
                    self.cte_y_std[ch,i_cur,:]+=f.fluxp_last_std[ch][1:]**2
                    cte_noise_std[ch,i_cur,:]+=(f.fluxp_last_std[ch][1:])**2*f.fluxp_used[ch][1:]
                    self.overscan_std[ch,i_cur]+=(f.over4_line_std[ch])**2
                #if flux_last==0. : flux_last=1e-6
                self.cte_flux[ch,i_cur]+=flux_last
        #    self.i_f+=1
        #fl=np.argsort(self.cte_flux,axis=1)
        ft=np.argsort(self.cte_ftime[0:self.i_f])
        l_ft=len(ft)
        #print('order in time ',ft)
        self.lmax=np.zeros((16),dtype=np.int16)
        # we take the number of amplifiers from the last file read 
        for ch in range(f.HduMax) :
            l_k=0
            cte_sig=np.zeros((l_ft))
            #for l in fl[ch,:] :
            for l in ft[:] :
                # protection against divide by 0 improbable ?
                if not(self.cte_flux[ch,l]>0.) : self.cte_flux[ch,l]=1.0e-6
                self.cte_y_s[ch,l_k,:]=self.cte_y[ch,l,:]*gain[ch]/self.nb_file[l]
                # remark that the 1/n below ...it's because sqrt(1/n) **2 is needed to get the error on the mean , and not the dispersion  . ...
                self.cte_y_s_std[ch,l_k,:]=np.sqrt(self.cte_y_std[ch,l,:])*gain[ch]/self.nb_file[l]
                self.cte_noise_s[ch,l_k]=np.sqrt(cte_noise_std[ch,l,2:].mean(axis=0)/(self.nb_file[l]))*gain[ch]
                self.cte_flux_s[ch,l_k]=self.cte_flux[ch,l]*gain[ch]/self.nb_file[l]
                self.cte_noise_s_std[ch,l_k]=np.sqrt(cte_noise_std[ch,l,2:].std(axis=0)/(self.nb_file[l])/(26))*gain[ch]
                l_k+=1
            for l in range(1,l_ft) :
                 if  self.cte_flux_s[ch,l]<self.cte_flux_s[ch,self.lmax[ch]] and self.cte_flux_s[ch,self.lmax[ch]] > 100000    :
                     self.lmax[ch]=l
                     break
                 self.lmax[ch]=l
            if len(self.cte_flux_s[ch,:])==1 : self.lmax[ch]=1
            self.ylev[ch,0:self.lmax[ch]]=(self.cte_y_s[ch,0:self.lmax[ch],0]+self.cte_y_s[ch,0:self.lmax[ch],1])/self.cte_flux_s[ch,0:self.lmax[ch]]/float(self.nb_pixel)
            #self.ylev_std[ch,0:self.lmax[ch]]=self.ylev[ch,0:self.lmax[ch]]*np.sqrt((self.cte_y_s_std[ch,0:self.lmax[ch],0]/self.cte_y_s[ch,0:self.lmax[ch],0])**2+(self.cte_y_s_std[ch,0:self.lmax[ch],1]/self.cte_y_s[ch,0:self.lmax[ch],1])**2)
            self.ylev_std[ch,0:self.lmax[ch]]=np.sqrt(self.cte_y_s_std[ch,0:self.lmax[ch],0]**2+self.cte_y_s_std[ch,0:self.lmax[ch],1]**2)/self.cte_flux_s[ch,0:self.lmax[ch]]/float(self.nb_pixel)
            # re-order and normalize Overscan data 
            self.overscan_std[ch,0:l_ft]=np.sqrt(self.overscan_std[ch,ft]/self.nb_file[ft])*gain[ch]
            # overscan stability
            self.over8_18[ch,0:self.lmax[ch]]=(self.cte_y_s[ch,0:self.lmax[ch],8:18]).mean(axis=1)
            self.over8_18_std[ch,0:self.lmax[ch]]=np.sqrt(np.sum(self.cte_y_s_std[ch,0:self.lmax[ch],8:18]**2,axis=1))/10.

        return
    def print_cte(self,ccd_name,nf=0):
        if self.serie :
            print('Serial CTE for %s ----------------------------------------------------------------' % (ccd_name) )
        else :
            print('  //   CTE for %s ----------------------------------------------------------------' % (ccd_name) ) 
            
        #
        print('Ch |    flux  |          1-CTE        |      Signal in Overscan       | Overscan Noise  Noise in        |')
        print('   |          |                       |   ov 1    ov 2    ov 8 to 18  |                 overscan corner |')
        for n in range(nf,self.i_f) :
            print('---------------------------------------------------------------------------------------------------------')
            for ch in range(16) :
                if n>=self.lmax[ch] :
                    print('%02d |  % 6.0f  | saturation (no eval) | % 6.02f  % 6.02f % 6.02f+/-%5.02f | % 5.02f+/-%5.02f        % 5.02f      |' % (
                        ch,
                        self.cte_flux_s[ch,n],
                        self.cte_y_s[ch,n,0],
                        self.cte_y_s[ch,n,1],
                        self.over8_18[ch,n],self.over8_18_std[ch,n],
                        self.cte_noise_s[ch,n],self.cte_noise_s_std[ch,n],
                        self.overscan_std[ch,n]))
                else : 
                    print('%02d |  % 6.0f  | %9.3g+/-%9.3g | % 6.02f  % 6.02f % 6.02f+/-%5.02f | % 5.02f+/-%5.02f        % 5.02f      |' % (
                      ch,
                      self.cte_flux_s[ch,n],
                      self.ylev[ch,n],self.ylev_std[ch,n],
                      self.cte_y_s[ch,n,0],
                      self.cte_y_s[ch,n,1],
                      self.over8_18[ch,n],self.over8_18_std[ch,n],
                      self.cte_noise_s[ch,n],self.cte_noise_s_std[ch,n],
                      self.overscan_std[ch,n])
                      )
        print('---------------------------------------------------------------------------------------------------------')
        return
                
    def plot_cte(self,ch,ccd_name,nf=0,on_screen=False,root_dir='.',unit='e-') :
        '''
        plot_cte(self,ch,ccd_name,nf=0,on_screen=False,root_dir='.')

        Plot the CTE results from cte class  per channel
        
             Parameters:
             ch        (int)  : channel index ( = hdu number -1 in data file )  to plot
             ccd_name  (str)  : ccd name : extra sting used in caption to identify this plot serie 
                                ( remark the ccd name itself , from file header is automaticaly added , this is more to identify run
                                ( or test level in plots label  
             nf        (int)  : index of first flux entry to plot (default=0)
             on_screen (bool) :  do we plot on display (or just save png on disk ) (default=False)
             root_dir  (str)  : top directory to save directory tree with plots
                                (default = '.' , directory used to save the plots will be ./raft_name/ccd_name/ch/ )
             unit      (str)  : unit of flux used ( e- or ADU ) 
        '''
        #
        root_plt=os.path.join(root_dir,self.first_file.raftbay,self.first_file.ccdslot,str(self.first_file.Hdu[ch])) 
        label_header=ccd_name+' '+self.first_file.raftbay+' '+self.first_file.ccdslot+' '+self.first_file.Image[ch]+' (hdu='+str(self.first_file.Hdu[ch])+')'
        # create the directorty
        os.makedirs(root_plt,exist_ok=True)
        # 
        xx=[max(np.min(self.cte_flux_s[:,nf:self.lmax[ch]])*.9,10.),min(2.0e5,np.max(self.cte_flux_s[:,nf:self.lmax[ch]])*1.1)]
        #

        pix_col=['b','c']
        pix_sym=['<','>']
        fig=plt.figure(figsize=(10,12))
        x=range(self.first,self.first+28)
        if self.serie :
            title="CTI Serial : "+label_header
            yv=5.0e-6
        else : 
            title="CTI //     : "+label_header
            yv=3.0e-6
        yy=[yv,yv]

        fig.suptitle(title,y=0.94)
        #fig.tight_layout()
        iplt=1
        ax=fig.add_subplot(3,3,iplt)
        ax.xaxis.set_label_position('top')
        ax.xaxis.tick_top()
        iplt+=1
        #ylev=(self.cte_y_s[ch,nf:self.lmax[ch],0]+self.cte_y_s[ch,nf:self.lmax[ch],1])/self.cte_flux_s[ch,nf:self.lmax[ch]]/float(self.nb_pixel)
        label='%02d' % self.first_file.Hdu[ch]
        #plt.plot(self.cte_flux_s[ch,nf:self.lmax[ch]],self.ylev[ch,nf:self.lmax[ch]],'o',color='r',label=label)
        plt.errorbar(self.cte_flux_s[ch,nf:self.lmax[ch]],self.ylev[ch,nf:self.lmax[ch]], yerr=self.ylev_std[ch,nf:self.lmax[ch]],fmt='o', ecolor='r',label=label)
        #print(self.ylev[ch,nf:self.lmax[ch]])
        #print(self.ylev_std[ch,nf:self.lmax[ch]])
        plt.plot(xx,yy,'g')
        if self.serie :
            plt.xlabel('<flux> of last column in '+unit)
        else :
            plt.xlabel('<flux> of last line in '+unit)
        plt.xlim(xx[0],xx[1])
        y_min=min(max(int(np.min(self.ylev[ch,nf:self.lmax[ch]])*500)/100.,1.e-8),1e-7)
        y_max=5e-5
        plt.ylim(y_min,y_max)
        plt.ylabel('1-CTE')
        plt.xscale('log')
        if (abs(y_max/y_min>80.)) : plt.yscale('log')
        # plt.locator_params(axis="both", tight=True, nbins=10)
        plt.legend()
        ax=fig.add_subplot(3,3,iplt)
        ax.xaxis.set_label_position('top')
        ax.xaxis.tick_top()
        iplt+=1
        y_min=1.
        y_max=0
        for pix in range(2) :
            ylev=self.cte_y_s[ch,nf:self.lmax[ch],pix]/self.cte_flux_s[ch,nf:self.lmax[ch]]/float(self.nb_pixel)
            #
            label="pix + %d " % (pix+1)
            y_min=min(y_min,np.min(ylev))
            y_max=max(0.,np.max(ylev))
            plt.plot(self.cte_flux_s[ch,nf:self.lmax[ch]],ylev,pix_sym[pix],color=pix_col[pix],label=label)
        plt.plot(xx,yy,'g')
        if self.serie :
            plt.xlabel('<flux> of last column in '+unit)
        else :
            plt.xlabel('<flux> of last line in '+unit)
        #plt.ylabel('1-CTE')
        y_min=min(max(y_min*.5,1.e-8),5e-7)
        y_max=max(y_max*1.5,1e-5)
        plt.ylim(y_min,y_max)
        plt.xscale('log')
        if (y_max/y_min>80.) : plt.yscale('log')
        plt.xlim(xx[0],xx[1])
        plt.legend()
        ax=fig.add_subplot(3,3,iplt)
        ax.xaxis.set_label_position('top')
        ax.xaxis.tick_top()
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        iplt+=1
        y_min=0.
        y_max=0.
        for pix in range(2) : 
            label="pix + %d " % (pix+1)
            plt.plot(self.cte_flux_s[ch,nf:self.lmax[ch]],self.cte_y_s[ch,nf:self.lmax[ch],pix],pix_sym[pix],color=pix_col[pix],label=label)
        y_min=min(-1.,np.min(self.cte_y_s[ch,nf:self.lmax[ch],0:1]))
        y_max=max(1.,np.max(self.cte_y_s[ch,nf:self.lmax[ch],0:1])*1.1)
        if self.serie :
            plt.xlabel('<flux> of last column in '+unit)
        else :
            plt.xlabel('<flux> of last line in '+unit)
        plt.ylabel('signal in overscan pixel(s) in '+unit)
        plt.xscale('log')
        if y_max > 10. :
            plt.yscale('symlog')
        plt.plot([xx[0],xx[1]],[0.,0.],'--',color='black')
        #
        plt.ylim(y_min,y_max)    
        plt.xlim(xx[0],xx[1])
        plt.legend(loc=2)
        #plt.xticks(ticks_flux)
        #
        ax=fig.add_subplot(3,3,iplt)
        iplt+=1
        #
        xx=[self.cte_flux_s[ch,nf]*0.9,self.cte_flux_s[ch,self.lmax[ch]-1]*1.1]
        yy=[0.,0.]
        plt.plot(xx,yy,'b--')
        plt.errorbar(self.cte_flux_s[ch,nf:self.lmax[ch]],
                      self.over8_18[ch,nf:self.lmax[ch]],yerr=self.over8_18_std[ch,nf:self.lmax[ch]],fmt='o',color='r', ecolor='r',label='Signal Overscan[8:18]')
        if self.serie :
            plt.xlabel('<flux> of last column in '+unit)
            plt.ylabel('signal in '+unit+' in serial Overscan')
        else :
            plt.xlabel('<flux> of last line in '+unit)
            plt.ylabel('signal in '+unit+' in // Overscan')
        
        plt.xscale('log')
        plt.ylim(min(np.min(self.over8_18[ch,nf:self.lmax[ch]])*1.2,-0.5),min(10.,max(0.5,np.max(self.over8_18[ch,nf:max(nf+1,self.lmax[ch]-1)])*1.5)))
        plt.legend(loc=2)       
        ax=fig.add_subplot(3,3,iplt)
        ax.yaxis.set_label_position("right")
        iplt+=1
        plt.plot(self.cte_flux_s[ch,nf:self.lmax[ch]],self.overscan_std[ch,nf:self.lmax[ch]],'<',label='Corner Noise')
        plt.errorbar(self.cte_flux_s[ch,nf:self.lmax[ch]],
                     self.cte_noise_s[ch,nf:self.lmax[ch]],yerr=self.cte_noise_s_std[ch,nf:self.lmax[ch]],fmt='o',color='r', ecolor='r',label='Frame Noise')
        try :
            mean_noiseV=np.array([self.cte_noise_s[ch,ii]   for ii in range(nf,self.lmax[ch])   if self.cte_flux_s[ch,ii] > 1000 and self.cte_flux_s[ch,ii] < 50000])
            if len(mean_noiseV)>0 : 
                mean_noise=mean_noiseV.mean()            
                xx=[self.cte_flux_s[ch,nf],self.cte_flux_s[ch,self.lmax[ch]-1]]
                yy=[mean_noise,mean_noise]
                plt.plot(xx,yy,'b--')
        except :
            pass
        if self.serie :
            plt.xlabel('<flux> of last column in '+unit)
            plt.ylabel('Noise from Serial Overscan in '+unit)
        else :
            plt.xlabel('<flux> of last line in '+unit)
            plt.ylabel('Noise from // Overscan in '+unit)
        plt.xscale('log')
        ymin_cc=3.
        ymax_cc=max(min(30.,np.max(self.cte_noise_s[ch,nf:max(nf+1,self.lmax[ch]-2)])*1.2),10.)
        #print(ymax_cc,np.max(self.cte_noise_s[ch,nf:max(nf+1,self.lmax[ch]-5)])*1.5)
        plt.ylim(ymin_cc,ymax_cc)
        #if ymax_cc > 20  :
        plt.legend(loc=2)
        #else :
        #    plt.legend(loc=3)       
        #
        ax=fig.add_subplot(3,3,iplt)
        iplt+=1
        flux=0.
        l_last=nf
        #lmax=len(self.cte_flux_s[ch,:])
        count=0
        im=0
        #
        y_min=0.
        y_max=0.
        #
        step_size=12
        nb_step=3
        last_step=nb_step-1
        nb_plot_max=step_size*nb_step
        # crude estimate of the log progretion , to estimate what step to get not too many plots (we average  fluxes to limit the number of plots ) 
        if (self.lmax[ch]-nf < nb_plot_max ) :
            ratio_step=0.9
        else :
            try :
                ratio_step=(self.cte_flux_s[ch,self.lmax[ch]]/self.cte_flux_s[ch,nf])**(-1./nb_plot_max)
            except :
                ratio_step=0.9
        #max_plt=max(int((self.lmax[ch]-nf)/4)+1,9)
        for l in range(nf,self.lmax[ch]) :
            if ((self.cte_flux_s[ch,l_last]/self.cte_flux_s[ch,l] < ratio_step ) and ( l_last < l )) :
                # first test to only plot result for point different enough , second test to be sure that we have already selected something , third test (l<lamx[ch] )   to avoid to plot too saturated guy
                if im>last_step-1 :
                    if self.serie :
                        label="%5.1f %s in last Col. " % (self.cte_flux_s[ch,l_last:l].mean(axis=0),unit)
                    else : 
                        label="%5.1f %s in last Line" % (self.cte_flux_s[ch,l_last:l].mean(axis=0),unit)
                else :
                    #    
                    label="%5.1f" % (self.cte_flux_s[ch,l_last:l].mean(axis=0))
                yplt=self.cte_y_s[ch,l_last:l,:].mean(axis=0)
                y_min=min(max(min(np.min(yplt)*1.2,0.),-10.),y_min)
                y_max=max(min(np.max(yplt)*1.2,100.),y_max)
                plt.plot(x,yplt,label=label)
                l_last=l
                count+=1
                if count == step_size and im<last_step and l < self.lmax[ch]-1 :
                    count = 0 
                    #plt.yscale('log')
                    if self.serie :
                        plt.xlabel('column number (serial overscan)')
                    else :
                        plt.xlabel('line number (// overscan)')
                    if im==0 or im==1 :
                        plt.ylabel('Overscan Signal in '+unit)
                    ymax=max(y_max,y_min+1.)
                    plt.ylim(y_min,y_max)
                    if im==0 :   plt.plot([x[0],x[-1]],[0.,0.],'--',color='black')
                    if y_max>80. :
                        plt.yscale('symlog')
                    plt.xlim(self.first,self.first+27)
                    if im == 0 :
                        ax.yaxis.set_label_position("right")
                        ax.yaxis.tick_right()
                    #    plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
                    #else :
                    plt.legend(loc=1)
                    ax=fig.add_subplot(3,3,iplt)
                    iplt+=1
                    y_min=0
                    y_max=0
                    im+=1
                    # update the ratio step with what is left to plot
                    try:
                        ratio_step=min(0.9,(self.cte_flux_s[ch,self.lmax[ch]]/self.cte_flux_s[ch,l])**(-1./step_size/(nb_step-im)))
                    except :
                        ratio_step=0.9
        if count !=0 or l==nf : 
            if self.serie :
                plt.xlabel('column number (serial overscan)')
            else :
                plt.xlabel('line number (// overscan)')
            if im<last_step : plt.ylabel('Overscan Signal in '+unit)
            if self.serie :
                label="%5.1f %s in last Col. " % (self.cte_flux_s[ch,l_last:self.lmax[ch]].mean(axis=0),unit)
            else : 
                label="%5.1f %s in last Line" % (self.cte_flux_s[ch,l_last:self.lmax[ch]].mean(axis=0),unit)
            yplt=self.cte_y_s[ch,l_last:self.lmax[ch],:].mean(axis=0)
            y_min=min(max(min(np.min(yplt)*1.2,0.),-10.),y_min)
            y_max=max(min(np.max(yplt)*1.2,100.),y_max)
            plt.plot(x,yplt,label=label)
            plt.xlim(self.first,self.first+27)
            plt.ylim(y_min,y_max)
            #plt.ylim(-10.,min(np.max(yplt)*1.2,100.))
            if y_max>80. :
                plt.yscale('symlog')
            plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
            #plt.legend()
            # Overscan noise 
            #ax=fig.add_subplot(3,3,iplt)
            #iplt+=1
        #
        if on_screen : plt.show()
        if self.serie :
            plotfile=root_plt+'/cte_serial.png' 
        else :
            plotfile=root_plt+'/cte_parallel.png'            
        fig.savefig(plotfile)
        if not(on_screen) : plt.close(fig)
        return
def cte_example():
    get_ipython().magic('matplotlib inline')
    print (' file =',sys.argv[1:-1],' TESTTYPE=',sys.argv[-1])
    selection=sys.argv[-1]
    file=Ifile(dirall=sys.argv[1:-1],fkey={'IMGTYPE':'Acquisition','TESTTYPE':sys.argv[-1]})
    #
    plt.interactive(1)
#file=Ifile(dirall=['/Users/antilog/scratch/e2v_190/20170314102625/*.fits'],fkey={})
#
    cte_data=cte(allfile=file.allfile,gain=[0.704650434205,0.68883578783,0.688459357874,0.696697494642,0.689209827484,0.696579402812,0.698973006751,0.689613072912,0.682880384357,0.696206655845,0.690349506621,0.691506176017,0.690763478766,0.689762341309,0.694801544092,0.850025229184 ])
    for ch in range(16) :
        cte_data.plot_extra(ch=ch,ccd_name=selection,nf=0,on_screen=True)


def fft_noise(h_all,channel=range(1,17),fplot=True,mean=False,start=1,int_pixel=1.8e-6,int_line=30.e-6,verbose=0,legend=True,xboundary=(20,500),yboundary=(30,2000),label='',color_v=None,two=True,axes=None,index=None) :
    cmap=plt.get_cmap('nipy_spectral')
    colors=[cmap(j)[:3] for j in np.linspace(0,1,17)]
    if color_v != None :
        color_val=color_v
    else : 
        color_val=channel
    #
    nb_l=yboundary[1]-yboundary[0]
    nb_c=xboundary[1]-xboundary[0]
    nb_file=len(h_all)
    nb_channel=len(channel)
    freq_x  = np.fft.rfftfreq(nb_c, d=int_pixel)[start:]
    freqf_x = np.flipud(1./np.fft.rfftfreq(nb_c,1)[start:])
    noise=np.zeros((nb_channel))
    #   
    # image area
    first_line,first_p_over,first_col,first_s_over=image_area(h_all[0])
    first_good_overs=first_s_over+2
    first_good_overp=first_p_over+2
    #
    for ich in range(nb_channel) :
        ch=int(channel[ich])
        for i_h  in range(nb_file) :
            h=h_all[i_h]
            #
            if i_h==0 :
                (n_y,n_x)=np.shape(h[1].data)
                 #delta time between 2 pixels from 2 # lines
                delta_line=int_pixel*n_x+int_line
                freq_y  = np.fft.rfftfreq(nb_l, d=delta_line)[start:]
                freqf_y = np.flipud(1./np.fft.rfftfreq(nb_l,d=delta_line/int_pixel)[start:])
                freq=np.append(freq_y,freq_x)
                freqf=np.append(freqf_x,freqf_y)
        #
            mean_line=np.median(h[ch].data[yboundary[0]:yboundary[1],:],axis=0)
            mean_column=np.median(h[ch].data[:,xboundary[0]:xboundary[1]],axis=1)
            if (ich==0 and i_h==0) or ( not(mean) and i_h==0 )  :
                ff_x=np.zeros((int(nb_c/2)))
                ff_y=np.zeros((int(nb_l/2)))
            for l in range(yboundary[0],yboundary[1]) :
                raw=h[ch].data[l,:]-mean_line
                to_fft=raw[xboundary[0]:xboundary[1]]-raw[first_good_overs:].mean()
                ff_x+=np.absolute(np.fft.rfft(to_fft))[start:]
            for c in range(xboundary[0],xboundary[1]) :
                raw=h[ch].data[:,c]-mean_column
                to_fft=raw[yboundary[0]:yboundary[1]]-raw[first_good_overp:].mean()
                ff_y+=np.absolute(np.fft.rfft(to_fft))[start:]
            noise[ich]+=(h[ch].data[yboundary[0]:yboundary[1],first_good_overs:].std())**2
            if verbose>1 : print ('channel %d noise %3.3f  Overscan dispersion = %3.3f '%(ch,h[ch].data[yboundary[0]:yboundary[1],first_good_overs:].std(),(h[ch].data[yboundary[0]:yboundary[1],first_good_overs:].mean(axis=1)).std()))
            if (i_h==nb_file-1)  and ( ich==len(channel)-1 or not(mean) )     :
                if mean :
                    # en fait on doit / par le nombr d bin de la fft , pas du signal ...facteur 2 ? 
                    ff_xn=ff_x/nb_l/nb_c/nb_file/nb_channel/2.
                    ff_yn=ff_y/nb_l/nb_c/nb_file/nb_channel
                    xnorm=np.append(ff_yn,ff_xn)
                    label_ch=label+'<'+','.join(map(str,channel))+'>'
                else :
                    ff_xn=ff_x/nb_l/nb_c/nb_file/2.
                    ff_yn=ff_y/nb_l/nb_c/nb_file
                    xnorm=np.append(ff_yn,ff_xn)
                    label_ch=label+'%d' % (ch)
                if fplot :
                    if two :
                        if index!=None :
                            axes[index[0]].plot(freq_x,ff_xn,label=label_ch,color=colors[color_val[ich]])
                            axes[index[1]].plot(freq_y,ff_yn,label=label_ch,color=colors[color_val[ich]])
                        else : 
                            plt.plot(freq_x,ff_xn,label=label_ch,color=colors[color_val[ich]])
                            plt.plot(freq_y,ff_yn,label=label_ch,color=colors[color_val[ich]])
                    else : 
                        if index!=None :
                            axes[index[0]].plot(freq,xnorm,label=label_ch,color=colors[color_val[ich]])
                        else :
                            plt.plot(freq,xnorm,label=label_ch,color=colors[color_val[ich]])
                else :
                    if two :
                        ff_xnf=np.flipud(ff_xn)
                        ff_ynf=np.flipud(ff_yn)
                        if index!=None :
                            axes[index[0]].plot(freqf_x,ff_xnf,label=label_ch,color=colors[color_val[ich]])
                            axes[index[1]].plot(freqf_y,ff_ynf,label=label_ch,color=colors[color_val[ich]])
                        else :
                            plt.plot(freqf_x,ff_xnf,label=label_ch,color=colors[color_val[ich]])
                            plt.plot(freqf_y,ff_ynf,label=label_ch,color=colors[color_val[ich]])
                    else : 
                        xnormf=np.flipud(xnorm)
                        if index!=None :
                            axes[index[0]].plot(freqf,xnormf,label=label_ch,color=colors[color_val[ich]])
                        else :
                            plt.plot(freqf,xnormf,label=label_ch,color=colors[color_val[ich]])
                if verbose :
                    argsort=np.argsort(xnorm)
                    sort=np.sort(xnorm)
                    ff_mean=np.mean(xnorm)
                    ff_sum=np.sum(xnorm)
                    # do a quick and durty bias on the sigma ...
                    ff_sig=xnorm.std()
                    #for i in range(1,len(argsort)) :
                    #    if xnorm[argsort[-i]] < ff_mean+3*ff_sig : break
                    #ff_sig=sort[0:-i].std()
            # 
                    print(' sum(fft) %g  <fft level>  %g   fft dispersion %g  ' % (ff_sum,ff_mean,ff_sig))
                    if not(fplot) :
                        xnormf=np.flipud(xnorm)
                    for i in range(1,len(argsort)) :
                        if xnorm[argsort[-i]] < ff_mean+2*ff_sig : break
                        if i>1 and ( ( argsort[-i]+1 in  argsort[-i+1:] ) or   ( argsort[-i]-1 in  argsort[-i+1:] ))  : continue
                        print (' fft bin %d ,  delta(pixels) %f (%6.1f Hz) : %g ' % (argsort[-i]+1,freqf[-argsort[-i]-1],freq[argsort[-i]],xnorm[argsort[-i]]))
            #
        if legend :
            if fplot :
                plt.xlabel('Noise in Hz')
            else :
                plt.xlabel('Noise Period in Pixel(s)')
            plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
        #plt.xscale('log')
        #plt.yscale('log')
    noise=np.sqrt(noise/nb_file)
    return freq,xnorm,noise



def for_ever(top_dir='/data/frames',do_fft=False,do_cte=False,xboundary=(20,500),yboundary=(30,2000),gain=[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]) : 
    filename_old=''
    old_dir=''
    while True : 
        dir=top_dir+'/'+time.strftime('%Y%m%d')
        #dir=top_dir
        files = glob.glob(dir+'/*.fz')
        if dir != old_dir : 
            print (' Scan Dir %s ========================================' % (dir) )
            old_dir=dir
        if  files:                
            _,new_file =os.path.split( max(files, key=os.path.getctime) )
            if new_file != filename_old :
                filename=dir+'/'+new_file
                # so it's a new file , but it could be not fully written 
                time.sleep(5)
                all_file=Ifile(dirall=[filename],single_t=False)
                filename_old=new_file
                print ('%s -------------------------------------------------------------------------' % (new_file) )
                print ('Ch |   mean  median  med. -     std    | <ov //>      std       std  | <ov S.>     std       std   |')
                print ('   |  image   image  S_over     image  |            ov //    Fix//o  |            ov S.    Fix S.  |')
                for ch in range(16) : 
                    print ('%02d | % 6.0f  % 6.0f  % 6.0f  % 8.02f  | % 6.0f  % 8.02f  % 8.02f  | % 6.0f  % 8.02f  % 8.02f  |' % (ch,
                                                                            all_file.all_file[0].Mean[ch],
                                                                            all_file.all_file[0].Median[ch],
                                                                            all_file.all_file[0].Median[ch]-all_file.all_file[0].MeanSScan[ch],
                                                                            all_file.all_file[0].Std[ch],
                                                                            all_file.all_file[0].MeanPScan[ch],
                                                                            all_file.all_file[0].StdPScan[ch],
                                                                            all_file.all_file[0].StdPScanOS[ch],
                                                                            all_file.all_file[0].MeanSScan[ch],
                                                                            all_file.all_file[0].StdSScan[ch],
                                                                            all_file.all_file[0].StdSScanOS[ch]) )
                print ('----------------------------------------------------------------------------------------------------')
                if do_cte :
                    cte_s=cte(all_file=all_file.all_file,gain=gain,serie=True)
                    ccd=new_file
                    cte_s.print_cte(ccd_name=ccd,nf=0)
                    #for ch in range(16) :
                    #    cte_s.plot_cte(ch=ch,ccd_name=ccd,nf=0,on_screen=True)
                    #    plt.show()
                    cte_p=cte(all_file=all_file.all_file,gain=gain,serie=False)
                    cte_p.print_cte(ccd_name=ccd,nf=0)
                    #for ch in range(16) :
                    #    cte_p.plot_cte(ch=ch,ccd_name=ccd,nf=0,on_screen=True)
                    #    plt.show()

                if do_fft :
                    fitsfile=pyfits.open(filename)
                    fft_it([fitsfile],channel=range(1,9),xboundary=xboundary,yboundary=yboundary)
                    fft_it([fitsfile],channel=range(9,17),xboundary=xboundary,yboundary=yboundary)
                    plt.show()
                    fitsfile.close()
        time.sleep(2)
    return
 
