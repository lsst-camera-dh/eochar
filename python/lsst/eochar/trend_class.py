# trend Object definition
# root class , that inherit from pytno dictionary , and will allow to addres the objects like a dictionnary :
#    obj['key']  or obj.key 
class DirDict(dict):
    def __init__(self):
        super().__init__()
        if trend.InTrend :
            self.__dict__ = self
        else :
            raise Exception('Addressing a trend object directly is not allowed', 1 )
            
    def __setitem__(self, key, item):
        if trend.InTrend :
            super().__setitem__( key, item)
        else :
            raise Exception('Addressing a trend object directly is not allowed', 1 )
            
    
# This is the class that contain all the infornmation from a trending entry 
class trend_obj(DirDict):
    def __init__(self,data,time,atype,values=None):
        super().__init__()
        self.data=data
        self.time=time
        self.type=atype
        self.values=values
        self.sync={}
    def fromtimestamp(self) :
    # entry : trending data name , ex : 'focal-plane/R11/Reb1/Temp10'
        # convert the time which is unix time in ms , to a standard datetime format 
        return np.array([datetime.datetime.fromtimestamp(self.time[i]/1000.) for i in range(len(self.time))])
from bs4 import BeautifulSoup as bs 
# the top object , that associate a given trending name , to its trending data :
# obj['trending_data_path'].data ... ( see trend_obj for all available sub object ) 
# 
class trend(DirDict):
    InTrend=False
    def __init__(self,start=None , dur='1h' , get=[] , verbose=False):
        trend.InTrend=True
        super().__init__()
        #self.__sart_time=start
        #self.__duration=dur
        #self.__get_argument=get
        trend.InTrend=False
        for to_get in get :
            # for each entry  'get' query the ternding db 
            if verbose : print('Reading trending for ',to_get,' starting at ',start,' for ',dur )
            # identify first all the db entries name associated to the current 'to_get'
            all_trend = !python  $trender_path --site summit --start $start --dur $dur --match $to_get
            data=[]
            for l in all_trend : 
                if l[0]=="#" : continue 
                all=l.split(' ')
                try :
                    data.append(all[7])
                except :
                    print('Error , unexpected data retruned by trender.py : \n',all_trend)
                    sys.exit(1)
            # read the entries  
            cmd='python %s --site summit  --start %s  --dur %s --xml %s ' % (trender_path,start,dur,to_get)
            result = subprocess.run(cmd ,shell=True,capture_output=True)
            content=result.stdout.decode("utf-8")
            #b=%get_magic_outb run -i $trender  --start $start --dur $dur --xml $to_get
            #content=b.decode("utf-8") 
            bbs_content = bs(content, features="xml") 
            # convert them in a dictionary , time , data , type of data 
            for data_cur in data : 
                con_meta=bbs_content.find_all(path=data_cur)[0].find_all("channelmetadatavalue")
                con=bbs_content.find_all(path=data_cur)[0].find_all("datavalue")
                con_t=bbs_content.find_all(path=data_cur)[0].find_all("axisvalue")
                if verbose : print(data_cur)
                for i in range(len(con_meta)) :
                    if con_meta[i]['name']=='type' :
                        data_type_m=con_meta[i]['value']
                        if 'Enumeration' in data_type_m :
                            data_type_m='Enumeration'
                            data_type='str'
                        elif data_type_m=='boolean' :
                            data_type='str'
                        else :
                            data_type=data_type_m
                    if verbose : print (con_meta[i]['name'],'==>',con_meta[i]['value']) 
                isize=len(con)
                if isize< 1 : continue
                time=np.zeros((isize),dtype=np.int64)
                if data_type=='str' :
                    if data_type_m=='boolean' : 
                        value=np.zeros((isize),dtype=int)
                    else :
                        value=np.zeros((isize),dtype=object)
                else :
                    value=np.zeros((isize))
                for icur in range(len(con)):
                    time[icur]=int(float(con_t[icur]['value']))
                    if data_type=='str' : 
                        val=con[icur]['strValue']
                        if data_type_m=='boolean' :
                            if val=='true' :
                                value[icur]=1
                            elif val=='false' :
                                value[icur]=0
                            else :
                                print('Erreur content in boolean (entry=',icur,') =',val)
                                value[icur]=0
                        else : 
                            value[icur]= val           
                    else : 
                        value[icur]=float(con[icur]['value'])
                #
                if data_type=='str' and data_type_m !='boolean' : 
                    list_entry=[]
                    for data_enum in value :
                        if not(data_enum in list_entry) :
                            list_entry.append(data_enum) 
                    list_entry.sort()
                    data_int=np.zeros((len(value)),dtype=np.int8)
                    i=0
                    for data_enum in value :
                        data_int[i]=list_entry.index(data_enum)
                        i+=1
                    trend.InTrend=True
                    self[data_cur]=trend_obj(data=data_int,time=time,atype=data_type_m,values=list_entry)
                    trend.InTrend=False
                else : 
                    trend.InTrend=True
                    self[data_cur]=trend_obj(data=value,time=time,atype=data_type_m)
                    trend.InTrend=False
    def save(self,out_file):
    # save the object in a pcikle file = out_file 
        out_put=open(out_file,'wb')
        pickle.dump(self,out_put)
        out_put.close()
        return
    def sync_causal_mean(self,s1,s2,tlength=600000):
    # compute a causal mean  of the data     between  ]t-tlenght,t[ 
    # entry : trending data name , ex : 'focal-plane/R11/Reb1/Temp10'
    # tlenght : time in ms ,  10 mins = 600000 ms 
    # so for 10 minutes in ms , to do   0<t-t0<10 , we will do -10/2<(t-to-10/2)<10/2 ==> abs(t-to-10/2) <10/2 
        delta=tlength/2.
        nb_entry=len(self[s2].time)
        cur_mean=np.zeros((nb_entry))
        for i in range(nb_entry) :
            arg=np.argwhere( abs(self[s2].time[i]-self[s1].time-delta)<=delta )
            if len(arg)!=0 : 
                cur_mean[i]=self[s1].data[arg].mean()
            else : 
                # to avoid empty slice , put 0 in it when slice is empty
                cur_mean[i]=0
        trend.InTrend=True
        a=DirDict()
        a['time']=self[s2].time
        a['data']=cur_mean
        a['sync_type']=tlength
        self[s1].sync[s2]=a 
        trend.InTrend=False
        return
    #
    def average(self,keys,tlength=120000):
        # keys a list with all the probe_path to average , should be an entry in the trend dictionary .
        # the time sync will be done on the timing of the first key : keys[0] . 
        for key in keys : 
            if not ( key in self.keys()) : 
                print('ERROR , no key ',key,' in the trend object provided')
                return None
        #
        av_probe=np.zeros(len(self[keys[0]].data))
        for keycur in keys :
            self.sync_causal_mean(keycur,keys[0],tlength=tlength) 
            av_probe+=self[keycur].sync[keys[0]].data
        av_probe/=len(keys)
        return av_probe

# pandas data from DM and images
class pdDM :
        def __init__(self,pickle_in=None , verbose=False):
            # input 
            # pickle_in , input pickle file , the one produced by Yassine , with DM , exposure and photodiode data 
            # verbose , if True print summary of what is read , default False  
            # open the pickle file with the DM data per event 
            self.pd=pd.read_pickle(pickle_in)
            # init per run and per event time for an easy usage
            # runU[run][0] start runU[run][1] end time of the run expressed in Unix in ms (int) 
            self.runU={}
            # runU[run][0] start runU[run][1] end time of the run expressed in datetime (localtime)  
            self.runDT={}
            # len[run] : number of event in the run 
            self.len={}
            # beg[run][0:len[run]]  : begin of the exposure (unclear what is taken as the start in DM , before clear , after ? ) 
            # end[run][0:len[run]]  : end of exposure , start of readout  
            self.beg={}
            self.end={}
            for run_cur in self.pd['run'].unique() :
                rcur=self.pd[self.pd['run']==run_cur]
                self.len[run_cur]=len(rcur)
                sel=rcur['timespan']
                self.runU[run_cur]=(int(((sel.iloc[0]).begin).to_value('unix')*1000.),int(((sel.iloc[-1]).end).to_value('unix')*1000.))
                self.runDT[run_cur]=(datetime.datetime.fromtimestamp(((sel.iloc[0]).begin).to_value('unix')),datetime.datetime.fromtimestamp(((sel.iloc[-1]).end).to_value('unix')))
                self.beg[run_cur]=np.array([int(sel.iloc[i].begin.to_value('unix')*1000.) for i in range(self.len[run_cur]) ])
                self.end[run_cur]=np.array([int(sel.iloc[i].end.to_value('unix')*1000.) for i in range(self.len[run_cur]) ])
            if verbose : 
            # print summary
                print('Data on events read from ',pickle_in)
                print('Total number of events ',len(self.pd),', Number of runs ',len(self.pd['run'].unique()),'(first=',self.pd[self.pd['run']!='unknown']['run'].iloc[0],' last=',self.pd[self.pd['run']!='unknown']['run'].iloc[-1],')')  
                print('First data read produced at ',datetime.datetime.fromtimestamp(self.pd['timespan'].iloc[0].begin.to_value('unix')).astimezone().isoformat())
                print('Flast data read produced at ',datetime.datetime.fromtimestamp(self.pd['timespan'].iloc[0].begin.to_value('unix')).astimezone().isoformat())
            return 
        #
        def probe(self,run,quantity,debug=False):
            offset=15
            # set start from all_run[run] - 15 min 
            tstart=(self.runDT[run][0]-datetime.timedelta(minutes = offset)).isoformat()
            # set end   from all_run[run][1]-[0]  + 15 min 
            tstop=(self.runDT[run][1]+datetime.timedelta(minutes = offset)).isoformat()
            # copute duration in s 
            tdur=int( ((self.runDT[run][1]+datetime.timedelta(minutes = offset))-(self.runDT[run][0]-datetime.timedelta(minutes = offset))).total_seconds())
            #
            #
            dbtr=trend(tstart,tdur,[quantity])        
            if debug : print('start run=',tstart,' end run=',tstop,' dur=',tdur,' probe list=',dbtr.keys())   
            return dbtr
        def average_probe(self,run,q_root,q_var='.',tlength=200000,debug=False):
            offset=15
            # set start from all_run[run] - 15 min 
            tstart=(self.runDT[run][0]-datetime.timedelta(minutes = offset)).isoformat()
            # set end   from all_run[run][1]-[0]  + 15 min 
            tstop=(self.runDT[run][1]+datetime.timedelta(minutes = offset)).isoformat()
            # copute duration in s 
            tdur=int( ((self.runDT[run][1]+datetime.timedelta(minutes = offset))-(self.runDT[run][0]-datetime.timedelta(minutes = offset))).total_seconds())
            #
            #
            quantity = q_root+q_var  
            dbtr=trend(tstart,tdur,[quantity]) 
            tkey=list(dbtr.keys())[0]
            for keycur in dbtr.keys() :
                dbtr.sync_causal_mean(keycur,tkey,tlength=tlength)
            av_probe=np.zeros(len(dbtr[tkey].sync[tkey].data))
            nb=0
            for keycur in dbtr.keys() : 
                av_probe+=dbtr[keycur].sync[tkey].data
            av_probe/=len(dbtr.keys())
            if debug : print('start run=',tstart,' end run=',tstop,' dur=',tdur,' nb measure=',len(av_probe),' first ev=',dbtr[keycur].sync[tkey].time[0],' last ev=',dbtr[keycur].sync[tkey].time[-1],' time between first and last in h=',(dbtr[keycur].sync[tkey].time[-1]-dbtr[keycur].sync[tkey].time[0])/3600000)   
            return dbtr,tkey,av_probe
