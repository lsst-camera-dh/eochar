from astropy.io import fits
import pandas as pd
import numpy as np
from lsst.resources import ResourcePath
from astropy.table import Table
import os 


#
ch16=['C10','C11','C12','C13','C14','C15','C16','C17','C07','C06','C05','C04','C03','C02','C01','C00']
photo_path='s3://lfa@rubinobs-lfa-cp/'
# subroutine to get the information for a given run    (fast ) 
def get_index_info(butler,query,nb_ccd=True):    
    df_exposure = pd.DataFrame(columns= ['science_program','id','obs_id','group','physical_filter','exposure_time','dark_time','observation_type','observation_reason','day_obs','seq_num','timespan','nb_ccd','header','uri','photodiode'])
    df_exposure['uri']=df_exposure['uri'].astype('object')
    df_exposure['header']=df_exposure['header'].astype('object')
    df_exposure['photodiode']=df_exposure['photodiode'].astype('object')
    # where_query="'%s'" % (query)
    for i, ref in enumerate(butler.registry.queryDimensionRecords('exposure',where=query).order_by("timespan.end")):
        #if nb_ccd :
        #    # then each exposure how many detector is there 
        #    where_query="exposure.id = %s  " % (ref.id)
        #    results = butler.registry.queryDimensionRecords( 'detector',
        #                                         datasets='raw' ,
        #                                         where=where_query  )
        #    results = len(list( set(results) ))
        #else :
        #    results = None 
        df_exposure.loc[i] = [ref.science_program,ref.id,ref.obs_id,ref.group,ref.physical_filter,ref.exposure_time,ref.dark_time,ref.observation_type,ref.observation_reason,ref.day_obs,ref.seq_num,ref.timespan,None,None,None,None]      
    return df_exposure

# subroutine to get the information for a run + files uri + photodiode flux information     
def get_index(butler,query_cur,uri_fast=True,photo_use=True,channel='103',verbose=True,header_use=True,header_dm=True,repo_root=None,fsspec_kwargs=None) : 
    header_val={'TEMPAVG':17,'BSSVBS':17}
    #df_ccob= pd.DataFrame(columns= ccob_val]   
    # get all the info on the exposures 
    #dsrefs = get_dsrefs(run_cur,butler)
    if verbose :  
        t0=time.time()
        print('Start queries to identify all exposures of  %s' % (query_cur))
    df = get_index_info(butler,query_cur,nb_ccd=False)
    #
    if verbose :  
        dt=time.time()-t0
        print('Delta t = %s , Start queries to identify all files associated to the %d exposures of  %s '  % (dt,len(df),query_cur))
    #
    nb_ccd=0
    nb_file=0
    # get the name of all fits file 
    for iid in range(len(df)) : 
        detector=None
        uri=np.zeros((205),dtype=np.object_)
        if uri_fast :
            where_query="exposure.id = %s  " % (int(df.loc[iid,'id']))
            results = butler.registry.queryDimensionRecords( 'detector',
                                                 datasets='raw' ,
                                                 where=where_query  )
            detec = list( set(results) )
            nb_ccd=len(detec)
            for detec_cur in detec : 
                detec_dic=detec_cur.toDict()
                Raft_CCD=detec_dic['full_name']
                instrument=detec_dic['instrument']
                detector=int(detec_dic['id'])
                uri[detector]='%s%s/%s/%s/%s_%s.fits' % (repo_root,instrument,df.loc[iid,'day_obs'],df.loc[iid,'obs_id'],df.loc[iid,'obs_id'],Raft_CCD)
        else : 
            nb_ccd=0
            for i in range(205) :
                try :
                    dataId = {'exposure': int(df.loc[iid,'id']), 'detector': i}
                    uri[i]=str(butler.getURI('raw',dataId)) 
                    # TBD : ON DOIT ENLEVER le embargo@ du bucket / file name !!!!!!!!
                    #detector=i
                    nb_ccd+=1
                except: 
                    continue
        df.loc[iid,'nb_ccd']=nb_ccd
        nb_file+=nb_ccd
        #
        df.iat[iid, df.columns.get_loc('uri')]=uri
    # do we read header  info ?
    if header_use  :
        if verbose :  
            dt=time.time()-t0
            print('Delta t = %s , Start read header data for each exposures '  % (dt))
        nb_header=0
        for i in range(len(df)) :
                header_data=np.zeros((205),dtype=np.object_)
                for detector in range(len(df.loc[i,'uri'])) :
                    if df.loc[i,'uri'][detector]!=0 :
                        # This is the minimal open to not load the full file in memory ,
                        # in the past it was taking .025 s per file , DM take ~ 2s per file to access the header
                        # now we are .11 s per file for the direct access , 0.05 s per file for DM ... ????
                        if header_dm :
                            hdu=butler.get('raw.metadata',instrument='LSSTCam',detector=detector,exposure=df.loc[i,'id'])
                            header_values={}
                            #print('3 ',df.loc[i,'uri'][detector],' ',time.time()-t0)
                            for header_cur,header_id in  header_val.items() :
                                try :
                                    header_values[header_cur]=hdu[header_cur]
                                except : 
                                    header_values[header_cur]=None
                        else :
                            hdu=fits.open(df.loc[i,'uri'][detector],cache=False, fsspec_kwargs=fsspec_kwargs ) 
                            header_values={}
                            #print('3 ',df.loc[i,'uri'][detector],' ',time.time()-t0)
                            for header_cur,header_id in  header_val.items() :
                                try :
                                    header_values[header_cur]=hdu[header_id].header[header_cur]
                                except : 
                                    header_values[header_cur]=None
                                    hdu.close()
                        nb_header+=1
                    else :
                        for header_cur,header_id in  header_val.items() :
                            header_values[header_cur]=None
                    #print('4 ',i,' ',time.time()-t0)
                    header_data[detector]=header_values
                df.iat[i, df.columns.get_loc('header')]=header_data
                #print('5 ',i,' ',time.time()-t0)
                if verbose and i%50 == 0 :
                    dt=time.time()-t0
                    print('Delta t = %s , Read fits header for %d events and %d headers '  % (dt,i+1 , nb_header))
    # do we need photodiode info ?
    if photo_use  :
        if verbose :  
            dt=time.time()-t0
            print('Delta t = %s , Start queries to identify photodiode  data for each exposures '  % (dt))
         # select all flat
        iflat=df[df['observation_type']=='flat'].index.to_list()
        icol=df.columns.get_loc('photodiode')
        nb_flat=len(iflat)
        day=[]
        photodiode_start=[]
        photo_data=[]
        if nb_flat!= 0 :
            exposure_start=np.zeros((nb_flat))
            icur=0
            for i in iflat : 
                exposure_start[icur]=df.loc[i,'timespan'].begin.to_value('unix_tai')
                icur+=1
                # for each day select all photodiode data of the day
                date='%d' % (df.loc[i,'day_obs'])
                if not(date in day ) :
                    day.append(date)
                    #
                    basepath='%sElectrometer:%s/fits/%s/%s/%s/' % (photo_path,channel,date[:4],date[4:6],date[6:])
                    rp = ResourcePath(basepath)
                    day_photo=list(rp.walk())
                    if len(day_photo) > 0 and len(day_photo[0]) > 2    :
                        for photo_file in day_photo[0][2] :
                            photo_value={}
                            photo_value['uri']=basepath+photo_file
                            photo_file_name=ResourcePath(basepath+photo_file)
                            with photo_file_name.open("rb") as fd:
                                with fits.open(fd) as hdul :
                                    # read the photodiode start of integration in the header (unit : unix_tai) 
                                    photodiode_start.append(hdul[0].header['DATE-BEG'])
                                    table=Table.read(hdul)
                                    # read / compute a set of data from photodiode 
                                    photo_value['current_mean']=table['Signal'].mean()
                                    photo_value['current_std']=table['Signal'].std()
                                    photo_value['current_nb_sample']=len(table)
                                    photo_value['current_error']=photo_value['current_std']/np.sqrt(photo_value['current_nb_sample'])
                                    photo_value['scan_time']=hdul[0].header['SCANTIME']
                                    photo_value['scan_beg']=hdul[0].header['DATE-BEG']
                                    photo_value['scan_end']=hdul[0].header['DATE-END']
                                    photo_value['uri']=fd
                                    photo_data.append(photo_value)
            if len(photo_data) !=0 :
                photo_start=np.array(photodiode_start)
                icur=0
                # do the match between the exposure and the photo diode file
                for i in iflat : 
                    # for each flat  find get the closest photodiode to the start of the exposures
                    i_start=np.argmin(np.abs(exposure_start[icur]-photo_start[:]))
                    # we got it if it's within 10s 
                    if np.abs(photo_start[i_start]-exposure_start[icur])<10. :
                        #
                        df.iat[i, icol]=photo_data[i_start]
                    icur+=1
    #               
    if verbose :  
        dt=time.time()-t0
        print('Delta t = %s , Done , all data for  %s (nb exposure %d , nb files %d) are collected '  % (dt,query_cur,len(df),nb_file))
    return df
def GetAllDays(butler,verbose=True,instrument='LSSTCam'):
    list_days = []
    nb_event  = {}
    query="instrument='%s'" % (instrument)
    for i, ref in enumerate(butler.registry.queryDimensionRecords('exposure',where=query).order_by("exposure.timespan.end")):
        day = ref.day_obs
        if day not in list_days:
            list_days.append(day) 
            nb_event[day]=1
        else:
            nb_event[day]+=1
            continue
    if verbose : 
        print('number of run',len(list_days))
        print(nb_event)
    return list_runs,nb_event 
def GetDay(butler,day_cur,repo_root=repo_root,instrument='LSSTCam',header_use=True,header_dm=True,fsspec_kwargs=fsspec_kwargs,write_panda=False,panda_path='/home/a/antilog/public_html/LsstCam/Index/'):
    PandaDir=os.path.join(panda_path,day_cur)
    if header_use : 
        PandaFile='%s/PandaDayHeaderIndex.pkl' % (PandaDir)
    else :
        PandaFile='%s/PandaDayIndex.pkl' % (PandaDir)
    #        
    try :
         df_temp=pd.read_pickle(PandaFile)
         df=df_temp
         print('Read Index data for day %s from %s ' % (day_cur,PandaFile))
    except:
         print('file ',PandaFile,' porbably not there , we get it from data')
         query_cur="instrument = '%s' and day_obs = %s" % (instrument,day_cur)
         df=get_index(butler,query_cur,header_use=header_use,header_dm=header_dm,repo_root=repo_root,fsspec_kwargs=fsspec_kwargs)
         if write_panda :
            os.makedirs(PandaDir,exist_ok=True)
            df.to_pickle(PandaFile) 
            print ('Panda PklFile=',PandaFile)
    return df
#
class ImageAna :
    def __init__(self,raw,verbose=True) :
        self.raw=raw
        self.det = raw.getDetector()
        self.nb_amp=len(self.det.getAmplifiers())
        self.vendor = self.det.getPhysicalType()
        det_ref=self.det['C01']
        self.det_name=self.det.getName()
        camera= LsstCam.getCamera()
        self.raw_det=camera[self.det_name]
        self.first_c=det_ref.getRawSerialPrescanBBox().getDimensions()[0]
        self.first_l=0
        self.amp_y_size=self.last_l=det_ref.getRawBBox().getDimensions()[1]
        self.amp_x_size=self.last_c=det_ref.getRawBBox().getDimensions()[0]
        self.im_y_size=det_ref.getRawDataBBox().getDimensions()[1]
        self.im_x_size=det_ref.getRawDataBBox().getDimensions()[0]
        self.first_c_over=self.first_c+self.im_x_size
        self.first_l_over=self.im_y_size
        self.image=np.zeros((self.nb_amp,self.last_l,self.last_c))
        if verbose : 
            print('CCD ',self.det_name)
            print('Number of lines read=',self.amp_y_size,' Number of colomns read=',self.amp_x_size)
            print('Number of lines in Image area=',self.im_y_size,' Number of colomns in Image area=',self.im_x_size)
            print('First line in Image area=',self.first_l,' First column in Image area=',self.first_c)
        # 
    def bias_cor(self,amp_cur,over_c='1D',over_l='1D',skip_c_over=2,skip_l_over=2) :
        # type of correction 
        #     None ==> Nothing 
        #     Cte ==> median over the full overscan , mean per line or column , than median over the mean ,  
        #     1D  ==> Mean per (line or collumn) , 1 correction per line or column
        # get the index of the current amp
        index = ch16.index(amp_cur)
        image_raw=self.raw.getImage()
        image_2d=image_raw[self.det[amp_cur].getRawBBox()].array
        raw_amp=self.raw_det[amp_cur]
        if raw_amp.getRawFlipX() :
            image_2d=image_2d[:,::-1]
        if raw_amp.getRawFlipY() :
            image_2d=image_2d[::-1,:]
        # Overscan correction
        # Over_c 
        self.over_c= over_c
        self.mean_over_per_line=np.mean(image_2d[:,self.first_c_over+skip_c_over:],axis=1)
        if over_c=='1D' :
            cor_l=self.mean_over_per_line[:,np.newaxis]
            cor_over_l=np.median(self.mean_over_per_line[self.first_l_over+skip_l_over:])
        elif over_c=='Cte' :
            cor_over_l=np.median(self.mean_over_per_line[self.first_l_over+skip_l_over:])
            cor_l=np.median(self.mean_over_per_line[self.first_l:self.first_l_over])
        else : 
            cor_over_l=0.
            cor_l=0.
        # over_l
        self.over_l= over_l
        # from // overscan  correction per colomn 
        # Over_l
        self.mean_over_per_column=np.mean(image_2d[self.first_l_over+skip_l_over:self.last_l,:],axis=0)-cor_over_l 
        if over_l=='1D' : 
            cor_c=self.mean_over_per_column
        elif over_l=='Cte' :
            cor_c=np.median(self.mean_over_per_column)
        else :
            cor_c=0.
        #
        self.image[index,:,:]=image_2d-(cor_l+cor_c)
        # create a view of the science part of the CCD image
        # self.science_image=self.image_raw[self.first_l:self.first_l_over,slef.first_c:self.first_c_over]
        #
    def SingleImageNorm(self):
        # display à la ds9 shape 
        # prepare a 4kx4k  CCD image , with amplifiers set at the right place ...there is a DM version which does this better...
        # but here you are in stand alone 
        # the default associated to the image area (pre-overscan excluded) 
        #
        #
        if self.nb_amp==16 :
            spf=np.zeros((self.im_y_size*2,self.im_x_size*8))
        else : 
            spf=np.zeros((self.im_y_size,self.im_x_size*8))
        # compute a relative gain by amplifier 
        relatG=np.zeros((self.nb_amp))
        for i in range(self.nb_amp) :
            relatG[i]=np.median(self.image[i][1700:1950,30:500])
        #
        tot=np.mean(relatG)
        relatG=relatG/tot
        #    
        if(self.vendor == 'E2V'):
             for i in range(self.nb_amp) :
                  if i<8 :
                       xx=i*self.im_x_size-1
                       if xx< 0 :
                            spf[2*self.im_y_size-1:self.im_y_size-1:-1,xx+self.im_x_size::-1]=self.image[i][self.first_l:self.first_l_over,self.first_c:self.first_c_over]/relatG[i]
                       else : 
                            spf[2*self.im_y_size-1:self.im_y_size-1:-1,xx+self.im_x_size:xx:-1]=self.image[i][self.first_l:self.first_l_over,self.first_c:self.first_c_over]/relatG[i]
                  else :
                       xx=(15-i)*self.im_x_size
                       spf[:self.first_l_over-self.first_l,xx:xx+self.im_x_size]=self.image[i][self.first_l:self.first_l_over,self.first_c:self.first_c_over]/relatG[i]
        else:
             for i in range(self.nb_amp) :
                  if i<8 :
                        xx=i*self.im_x_size-1
                        if xx < 0 : 
                            spf[2*self.im_y_size-1:self.im_y_size-1:-1,xx+self.im_x_size::-1]=self.image[i][self.first_l:self.first_l_over,self.first_c:self.first_c_over]/relatG[i]
                        else : 
                            spf[2*self.im_y_size-1:self.im_y_size-1:-1,xx+self.im_x_size:xx:-1]=self.image[i][self.first_l:self.first_l_over,self.first_c:self.first_c_over]/relatG[i]
                  else :
                        xx=(15-i)*self.im_x_size-1
                        if xx<0 : 
                            spf[:self.first_l_over-self.first_l,xx+self.im_x_size::-1]=self.image[i][self.first_l:self.first_l_over,self.first_c:self.first_c_over]/relatG[i]
                        else :
                            spf[:self.first_l_over-self.first_l,xx+self.im_x_size:xx:-1]=self.image[i][self.first_l:self.first_l_over,self.first_c:self.first_c_over]/relatG[i]
        return spf        
    def SingleImage(self):
        # display à la ds9 shape 
        # prepare a 4kx4k  CCD image , with amplifiers set at the right place ...there is a DM version which does this better...
        # but here you are in stand alone 
        # the default associated to the image area (pre-overscan excluded) 
        #
        #
        if self.nb_amp==16 :
            spf=np.zeros((self.im_y_size*2,self.im_x_size*8))
        else : 
            spf=np.zeros((self.im_y_size,self.im_x_size*8))
        #    
        if(self.vendor == 'E2V'):
             for i in range(self.nb_amp) :
                  if i<8 :
                       xx=i*self.im_x_size-1
                       if xx< 0 :
                            spf[2*self.im_y_size-1:self.im_y_size-1:-1,xx+self.im_x_size::-1]=self.image[i][self.first_l:self.first_l_over,self.first_c:self.first_c_over]
                       else : 
                            spf[2*self.im_y_size-1:self.im_y_size-1:-1,xx+self.im_x_size:xx:-1]=self.image[i][self.first_l:self.first_l_over,self.first_c:self.first_c_over]
                  else :
                       xx=(15-i)*self.im_x_size
                       spf[:self.first_l_over-self.first_l,xx:xx+self.im_x_size]=self.image[i][self.first_l:self.first_l_over,self.first_c:self.first_c_over]
        else:
             for i in range(self.nb_amp) :
                  if i<8 :
                        xx=i*self.im_x_size-1
                        if xx < 0 : 
                            spf[2*self.im_y_size-1:self.im_y_size-1:-1,xx+self.im_x_size::-1]=self.image[i][self.first_l:self.first_l_over,self.first_c:self.first_c_over]
                        else : 
                            spf[2*self.im_y_size-1:self.im_y_size-1:-1,xx+self.im_x_size:xx:-1]=self.image[i][self.first_l:self.first_l_over,self.first_c:self.first_c_over]
                  else :
                        xx=(15-i)*self.im_x_size-1
                        if xx<0 : 
                            spf[:self.first_l_over-self.first_l,xx+self.im_x_size::-1]=self.image[i][self.first_l:self.first_l_over,self.first_c:self.first_c_over]
                        else :
                            spf[:self.first_l_over-self.first_l,xx+self.im_x_size:xx:-1]=self.image[i][self.first_l:self.first_l_over,self.first_c:self.first_c_over]
        return spf
