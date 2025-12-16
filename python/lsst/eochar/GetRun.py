from astropy.io import fits
import pandas as pd
import numpy as np
#
ch16=['C10','C11','C12','C13','C14','C15','C16','C17','C07','C06','C05','C04','C03','C02','C01','C00']
# subroutine to get the information for a given run    (fast ) 
def get_run_info(butler,run,nb_ccd=True,instrument='LSSTCam'):    
    df_exposure = pd.DataFrame(columns= ['science_program','id','obs_id','physical_filter','exposure_time','dark_time','observation_type','observation_reason','day_obs','seq_num','timespan','nb_ccd','uri','ccob'])
    df_exposure['uri']=df_exposure['uri'].astype('object')
    df_exposure['ccob']=df_exposure['ccob'].astype('object')
    if len(run)<15 : 
        where_query="exposure.science_program = '%s'  " % (run)
    else :
        where_query=run
    for i, ref in enumerate(butler.registry.queryDimensionRecords('exposure',instrument=instrument,where=where_query).order_by("timespan.end")):
        if nb_ccd :
            # then each exposure how many detector is there 
            where_query="exposure.id = %s  " % (ref.id)
            results = butler.registry.queryDimensionRecords( 'detector',
                                                 datasets='raw' ,
                                                 where=where_query  )
            results = len(list( set(results) ))
        else :
            results = None 
        df_exposure.loc[i] = [ref.science_program,ref.id,ref.obs_id,ref.physical_filter,ref.exposure_time,ref.dark_time,ref.observation_type,ref.observation_reason,ref.day_obs,ref.seq_num,ref.timespan,None,None,None]      
    return df_exposure

# subroutine to get the information for a run + files uri + CCOB flux information     
def get_run(butler,run_cur,uri_fast=True,ccob_use=True,instrument='LSSTCam',verbose=True,repo_root=None,fsspec_kwargs=None) : 
    ccob_val=['TEMPLED1','TEMPLED2','TEMPBRD','CCOBLED','CCOBCURR','CCOBADC','CCOBFLST','PROJTIME','CCOBFLUX','DATEPBEG','MJDPBEG','DATEPEND','MJDPEND']
    #df_ccob= pd.DataFrame(columns= ccob_val]   
    # get all the info on the exposures 
    #dsrefs = get_dsrefs(run_cur,butler)
    if verbose :  
        t0=time.time()
        print('Start queries to identify all exposures of run %s' % (run_cur))
    df = get_run_info(butler,run_cur,nb_ccd=False,instrument='LSSTCam')
    #
    if verbose :  
        dt=time.time()-t0
        print('Delta t = %s , Start queries to identify all files associated to the %d exposures of run %s '  % (dt,len(df),run_cur))
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
    # do we need CCOB info ?
    if ccob_use  :
        if verbose :  
            dt=time.time()-t0
            print('Delta t = %s , Start queries to identify CCOB data for each exposures '  % (dt))
        nb_header=0
        for i in range(len(df)) : 
            if df.loc[i,'observation_type']=='flat' :
                detec_found=False
                for detector in range(len(df.loc[i,'uri'])) :
                    if df.loc[i,'uri'][detector]!=0 :
                        detec_found=True
                        break 
                if detec_found :
                    # This is the minimal open to not load the full file in memory ,
                    # it takes .025 s per file , DM take ~ 2s per file to access the header 
                    hdu=fits.open(df.loc[i,'uri'][detector],cache=False, fsspec_kwargs=fsspec_kwargs ) 
                    ccob_values={}
                    #print('3 ',df.loc[i,'uri'][detector],' ',time.time()-t0)
                    for ccob_cur in ccob_val : 
                        try :
                            ccob_values[ccob_cur]=hdu[0].header[ccob_cur]
                        except : 
                            ccob_values[ccob_cur]=None
                    #print('4 ',i,' ',time.time()-t0)
                    
                    df.iat[i, df.columns.get_loc('ccob')]=ccob_values
                    hdu.close()
                    #print('5 ',i,' ',time.time()-t0)
                    nb_header+=1
                if verbose and nb_header%50 == 0 :
                    dt=time.time()-t0
                    print('Delta t = %s , Read fits header for %d events '  % (dt,nb_header))
    #               
    if verbose :  
        dt=time.time()-t0
        print('Delta t = %s , Done , all data for run %s (nb exposure %d , nb files %d) are collected '  % (dt,run_cur,len(df),nb_file))
    return df
def GetAllRun(butler,verbose=True,instrument='LSSTCam'):
    list_runs = []
    nb_event  = {}
    query="instrument='%s'" % (instrument)
    for i, ref in enumerate(butler.registry.queryDimensionRecords('exposure',where=query).order_by("exposure.timespan.end")):
        run = ref.science_program
        if run not in list_runs:
            list_runs.append(run) 
            nb_event[run]=1
        else:
            nb_event[run]+=1
            continue
    if verbose : 
        print('number of run',len(list_runs))
        print(nb_event)
    return list_runs,nb_event 
def GetRunCur(butler,run_cur,repo_root=repo_root,instrument='LSSTCam',fsspec_kwargs=fsspec_kwargs):
    try :
         Pandafile='/home/a/antilog/public_html/LsstCam/IndexRun7/%s/PandaRun_pkl.pkl' % (run_cur)
         df=pd.read_pickle(Pandafile)
    except:
         print('file ',Pandafile,' porbably not there , we get it from data') 
         df=get_run(butler,run_cur,instrument=instrument,repo_root=repo_root,fsspec_kwargs=fsspec_kwargs)
    return df
