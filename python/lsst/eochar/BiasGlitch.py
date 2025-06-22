def ProcessBias(run_cur,raft_cur,ccd_cur,file90,plot=True,show=False,dist=1.5):
#def FileForBiasEstimator(run_cur,raft_cur,ccd_cur,file90):
    #
    #
    if ccd_cur in sensors_8ch :
        nb_amp=8
        ch=ch8
    else : 
        nb_amp=16
        ch=ch16
    #
    nb_file=len(file90)
    # read overscan for ~ nb_file files 
    prescan=np.zeros((nb_file,nb_amp,2048,first_col))
    overser=np.zeros((nb_file,nb_amp,2048))
    overpar=np.zeros((nb_file,nb_amp,576))
    ampnoise=np.zeros((nb_file,nb_amp,3))
    noise_met=['Image','Serial','//']
    # loop on files
    for ifits in range(nb_file):
        fits=pyfits.open(file90[ifits],cache=False, fsspec_kwargs=fsspec_kwargs)
        for iamp in range(nb_amp) :
            # pour itl prescan= col 0,1,2 ...la col 0 ne correspond à rien en pratique , 1 et 2 sont des pixels pre-scan
            prescan[ifits,iamp,:,:]=fits[iamp+1].data[:,0:first_col]
            # overscan dans la direction serie , on prend la moyen de l'overscan de chaque ligne :de la colomne 1er overscan+2 à la fin
            overser[ifits,iamp,:]=fits[iamp+1].data[:,first_cover+2:].mean(axis=1)
            # overscan dans la direction // , on prend la moyen de l'overscan de chaque colome :de la ligne 1er overscan+2 à la fin
            overpar[ifits,iamp,:]=fits[iamp+1].data[first_lover+2:,:].mean(axis=0)
            # noise in most of the image area
            ampnoise[ifits,iamp,0]=np.median(fits[iamp+1].data[first_line+50:first_lover-50,first_col+50:first_cover-50].std(axis=1))
            # noise in  the serial overscan
            ampnoise[ifits,iamp,1]=np.median(fits[iamp+1].data[first_line+3:,first_cover+3:].std(axis=1))
            # noise in  the // overscan 
            ampnoise[ifits,iamp,2]=np.median(fits[iamp+1].data[first_lover+3:,first_col+3:].std(axis=0))
           
        fits.close()
#    return nb_file,prescan,overser,overpar
#        
#def BiasEstimator(run_cur,raft_cur,ccd_cur,nb_file,prescan,overser,overpar,plot=True,dist=1.5):
    # dist=1.5 , distance min in ADU between clusters of bias 
    #
    if ccd_cur in sensors_8ch :
        nb_amp=8
        ch=ch8
    else : 
        nb_amp=16
        ch=ch16
    #

    all_met=['pre1','pre1-pre0','pre1&pre2','pre1&pre2-pre0','over_serie','over_serie-pre0','over_//','over_//-pre0']
    nb_met=len(all_met)
    #
    nb_cluster_max=6
    color=['r','g','b','yellow','purple','black']
    #
    cluster=np.zeros((nb_met,nb_amp,nb_cluster_max,nb_file),dtype=int)
    cluster_size=np.zeros((nb_met,nb_amp,nb_cluster_max),dtype=int)
    #icluster=-1*np.ones((nb_met,nb_amp),dtype=int)
    icluster=np.zeros((nb_met,nb_amp),dtype=int)
    mean=np.zeros((nb_met,nb_amp,nb_cluster_max))
    std=np.zeros((nb_met,nb_amp,nb_cluster_max))
    std_no_cluster=np.zeros((nb_met,nb_amp))
    ref=np.zeros((nb_met,nb_amp,nb_file))
    ref_mean=np.zeros((nb_met,nb_amp))
    ylabel=np.zeros((nb_met),dtype=object)
    #compute for all amplifier the clusters for the different method
    #for imet in [6]  :
    for imet in range(nb_met)  :
        met=all_met[imet]
        for iamp in range(nb_amp) :
            match met :
                case "pre1":
                    ref_mean[imet,iamp]=np.mean(prescan[:,iamp,:,1])
                    ref[imet,iamp,:]=np.mean(prescan[:,iamp,:,1], axis=1)-ref_mean[imet,iamp]
                    ylabel[imet]='Bias from <prescan[1]> \n centred to mean in run'
                case "pre1-pre0":
                    ref_mean[imet,iamp]=np.mean(prescan[:,iamp,:,1])-np.mean(prescan[:,iamp,:,0])
                    ref[imet,iamp,:]=np.mean(prescan[:,iamp,:,1], axis=1)-np.mean(prescan[:,iamp,:,0],axis=1)-ref_mean[imet,iamp]
                    ylabel[imet]='Bias from <prescan[1]>-<prescan[0]> \n centred to mean in run'
                case "pre1&pre2":
                    ref_mean[imet,iamp]=np.mean(prescan[:,iamp,:,1:3])
                    ref[imet,iamp,:]=np.mean(np.mean(prescan[:,iamp,:,1:3], axis=2),axis=1)-ref_mean[imet,iamp]
                    ylabel[imet]='Bias from <prescan[1]&[2]> '
                case "pre1&pre2-pre0":
                    ref_mean[imet,iamp]=np.mean(prescan[:,iamp,:,1:3])-np.mean(prescan[:,iamp,:,0])
                    ref[imet,iamp,:]=np.mean(np.mean(prescan[:,iamp,:,1:3], axis=2),axis=1)-np.mean(prescan[:,iamp,:,0], axis=1)-ref_mean[imet,iamp]
                    ylabel[imet]='Bias from <prescan[1]&[2]>-<prescan[0]> \n centred to mean in run'
                case "over_serie":
                    ref_mean[imet,iamp]=np.mean(overser[:,iamp,:])
                    ref[imet,iamp,:]=np.mean(overser[:,iamp,:],axis=1)-ref_mean[imet,iamp]
                    ylabel[imet]='Bias from <SerialOverscan>'
                case "over_serie-pre0":
                    ref_mean[imet,iamp]=np.mean(overser[:,iamp,:])-np.mean(prescan[:,iamp,:,0])
                    ref[imet,iamp,:]=np.mean(overser[:,iamp,:],axis=1)-np.mean(prescan[:,iamp,:,0], axis=1)-ref_mean[imet,iamp]
                    ylabel[imet]='Bias from <SerialOverscan>-<prescan[0]> \n centred to mean in run'
                case "over_//":
                    ref_mean[imet,iamp]=np.median(overpar[:,iamp,3:])
                    ref[imet,iamp,:]=np.median(overpar[:,iamp,3:],axis=1)-ref_mean[imet,iamp]
                    ylabel[imet]='Bias from med(//Overscan)  '
                case "over_//-pre0":
                    ref_mean[imet,iamp]=np.median(overpar[:,iamp,3:])-np.mean(prescan[:,iamp,:,0])
                    ref[imet,iamp,:]=np.median(overpar[:,iamp,3:],axis=1)-np.mean(prescan[:,iamp,:,0], axis=1)-ref_mean[imet,iamp]
                    ylabel[imet]='Bias from med(//Overscan)-<prescan[0]> \n centred to mean in run' 
                case _:
                    print("unknow config")
                    break  
    #for imet in range(nb_met)  :
    #    met=all_met[imet]
    #    for iamp in range(nb_amp) :                    
            for iref in range(nb_file) : 
                ref_cur=ref[imet,iamp,iref]
                found=False
                for icl in range(icluster[imet,iamp]):
                    if found : break
                    for icl_cur in range(cluster_size[imet,iamp,icl]) : 
                        dist_cur=abs(ref[imet,iamp,cluster[imet,iamp,icl,icl_cur]]-ref_cur)
                        if dist_cur<dist : 
                            cluster[imet,iamp,icl,cluster_size[imet,iamp,icl]]=iref
                            cluster_size[imet,iamp,icl]+=1
                            found=True
                            break
                if not(found) : 
                    # create a new cluster
                    if icluster[imet,iamp] == nb_cluster_max-1 :
                        # too many cluster , there is probably a defect in the image that disturb the process 
                        # so we should identify this clustering as failled for this amp , and fill its result in the last identified cluster
                        cluster[imet,iamp,nb_cluster_max-2,cluster_size[imet,iamp,nb_cluster_max-2]]=iref
                        cluster_size[imet,iamp,nb_cluster_max-2]+=1
                    else :
                        cluster[imet,iamp,icluster[imet,iamp],cluster_size[imet,iamp,icluster[imet,iamp]]]=iref
                        cluster_size[imet,iamp,icluster[imet,iamp]]+=1
                        icluster[imet,iamp]+=1
            #
            for icl in range(icluster[imet,iamp]) : 
                mean[imet,iamp,icl]=ref[imet,iamp,cluster[imet,iamp,icl,0:cluster_size[imet,iamp,icl]]].mean()
                std[imet,iamp,icl]=ref[imet,iamp,cluster[imet,iamp,icl,0:cluster_size[imet,iamp,icl]]].std()
            std_no_cluster[imet,iamp]=ref[imet,iamp,:].std()

    if plot : 
        # plot the results :
        # for each CCD plot the raw dispersion of the // over 
        fig=plt.figure(figsize=[16,16])
        # get the method and amp with the largest sigma 
        arg_max=np.unravel_index(np.argmax(std_no_cluster, axis=None), std_no_cluster.shape)
        imet_max=arg_max[0]
        iamp_max=arg_max[1]
        #                        
        txt='run %s , RAFT %s CCD %s \n median( // Overscan ) per event for each amplifier   ' % (run_cur,raft_cur,ccd_cur) 
        plt.suptitle(txt)
        # select the median verscan method
        imet=6
        # even if some entries are at 0 , I gues this is still ok     
        ymin=np.min(mean[imet,:,:]-5*std[imet,:,:])
        ymax=np.max(mean[imet,:,:]+5*std[imet,:,:])
        for iamp in range(nb_amp) :
            met=all_met[imet]
            plt.subplot(4,int(nb_amp/4),iamp+1)
            label='%s (hdu=%d), std=%6.3f ADU' % (ch[iamp],iamp+1,ref[imet,iamp,:].std())
            plt.gca().set_title(label)
            plt.plot(range(nb_file),ref[imet,iamp,:],color='black')
                #plt.plot(ref[cluster[imet,iamp,icl,0:cluster_size[imet,iamp,icl]]],color=color[icl],label=label)
            if iamp>=nb_amp-4 :
                label='in run event number'
                plt.xlabel(label)
            if iamp%4 == 0 : 
                label='Per event // overscan median (in ADU)  \n centred to mean in run'
                plt.ylabel(label)
            plt.ylim(ymin,ymax)                        
        if show : plt.show() 
        rawPlotFile='BiasOverPar'
        SaveFig(fig,rawPlotFile,run_cur=run_cur,raft_cur=raft_cur,ccd_cur=ccd_cur)
        #
        # for each CCD plot the noise estimated by diferent part of the bias / image 
        fig=plt.figure(figsize=[16,16])
         #                        
        txt='run %s , RAFT %s CCD %s \n noise in ADU  per event for each amplifier   ' % (run_cur,raft_cur,ccd_cur) 
        plt.suptitle(txt)
        # for each type onf noise
        for imet in range(3) : 
            for iamp in range(nb_amp) :
                met=all_met[imet]
                plt.subplot(4,int(nb_amp/4),iamp+1)
                if imet==0 :
                    label='%s,%s=%4.2f,%s=%4.2f' % (ch[iamp],noise_met[1],np.median(ampnoise[:,iamp,1]),noise_met[2],np.median(ampnoise[:,iamp,2]))
                    plt.gca().set_title(label)
                plt.plot(range(nb_file),ampnoise[:,iamp,imet],color=color[imet])
                if iamp>=nb_amp-4 :
                    label='in run event number'
                    plt.xlabel(label)
                    if iamp%4 == 0 : 
                        label='Measured Noise per image '
                        plt.ylabel(label)
                        #plt.ylim(ymin,ymax)                        
        if show :plt.show() 
        rawPlotFile='Noise'
        SaveFig(fig,rawPlotFile,run_cur=run_cur,raft_cur=raft_cur,ccd_cur=ccd_cur)
        #
        # for each CCD plot the amplifier with the largest dispersion plot all cluster methods results
        fig=plt.figure(figsize=[16,16])
        # get the method and amp with the largest sigma 
        arg_max=np.unravel_index(np.argmax(std_no_cluster, axis=None), std_no_cluster.shape)
        imet_max=arg_max[0]
        iamp_max=arg_max[1]
        #                        
        txt='run %s , RAFT %s CCD %s \n amp %s (hdu=%d) with bias estimator with largest sigma: %s  ' % (run_cur,raft_cur,ccd_cur,ch16[iamp_max],iamp_max+1,all_met[imet_max]) 
        plt.suptitle(txt)
        iamp=iamp_max
        imet=imet_max
        # even if some entries are at 0 , I gues this is still ok     
        ymin=np.min(mean[imet,:,:]-5*std[imet,:,:])
        ymax=np.max(mean[imet,:,:]+5*std[imet,:,:])

        for imet in range(nb_met) :
            met=all_met[imet]
            plt.subplot(4,2,imet+1)
            label='%s , amp %s (hdu=%d), mean=%6.1f  std=%6.3f ' % (all_met[imet],ch[iamp],iamp+1,ref_mean[imet,iamp],ref[imet,iamp,:].std())
            plt.gca().set_title(label)
            for icl in range(icluster[imet,iamp]) :
                if icl==icluster[imet,iamp]-1 :
                    if icl!=0 :
                        label='lot=%d, mean = %6.3f , std=%6.3f , dist %d-%d = %6.3f ' % (icl,mean[imet,iamp,icl],std[imet,iamp,icl],0,icl,(mean[imet,iamp,0]-mean[imet,iamp,icl])/np.sqrt(std[imet,iamp,icl]**2+std[imet,iamp,0]**2))
                    else :
                        label='lot=%d, mean = %6.3f , std=%6.3f ' % (icl,mean[imet,iamp,icl],std[imet,iamp,icl])
                else : 
                    label='lot=%d, mean = %6.3f , std=%6.3f , dist %d-%d / sig = %6.3f ' % (icl,mean[imet,iamp,icl],std[imet,iamp,icl],icl+1,icl,(mean[imet,iamp,icl+1]-mean[imet,iamp,icl])/np.sqrt(std[imet,iamp,icl]**2+std[imet,iamp,icl+1]**2))
                plt.plot(cluster[imet,iamp,icl,0:cluster_size[imet,iamp,icl]],ref[imet,iamp,cluster[imet,iamp,icl,0:cluster_size[imet,iamp,icl]]],'.',color=color[icl],label=label)
                #plt.plot(ref[cluster[imet,iamp,icl,0:cluster_size[imet,iamp,icl]]],color=color[icl],label=label)
            if imet>=nb_met-2 :
                label='in run event number'
                plt.xlabel(label)
            plt.ylabel(ylabel[imet])
            plt.ylim(ymin,ymax)                        
            plt.legend()
        if show :plt.show() 
        rawPlotFile='BiasEstimator'
        SaveFig(fig,rawPlotFile,run_cur=run_cur,raft_cur=raft_cur,ccd_cur=ccd_cur,hdu=iamp_max+1)
        # Does the best split in pre1+2 does a correct classification of over_// ? 
        fig=plt.figure(figsize=[16,16])
        # get the method and amp with the largest sigma 
        imet_max=3
        iamp_max=np.argmax(std_no_cluster[imet_max,:])
        #                        
        txt='run %s , RAFT %s CCD %s \n clustering with amp %s (hdu=%d) with bias estimator: %s  ' % (run_cur,raft_cur,ccd_cur,ch16[iamp_max],iamp_max+1,all_met[imet_max]) 
        plt.suptitle(txt)
        iamp=iamp_max
        imet=imet_max
        if nb_amp==8  :
            n_sub=2
        else :
            n_sub=4
        for iamp in range(nb_amp) :
            label='amp %s (hdu=%d)' % (ch16[iamp],iamp+1)
            plt.subplot(n_sub,4,iamp+1,title=label)
            for icl in range(icluster[imet_max,iamp_max]+1):
                plt.plot(ref[imet_max,iamp_max,cluster[imet_max,iamp_max,icl,:cluster_size[imet_max,iamp_max,icl]]],ref[7,iamp,cluster[imet_max,iamp_max,icl,:cluster_size[imet_max,iamp_max,icl]]],'.',color=color[icl])
                if icl==0 and iamp+1==((n_sub-1)*4+2):
                    # last line , we can put the x label 
                    label='Bias ( %s ) of amp %s (hdu=%d) clusterized (color) ' % (all_met[imet_max],ch16[iamp_max],iamp_max+1) 
                    plt.xlabel(label)
                if iamp%4==0 :
                    # first column , we can put the y label 
                    label='Bias ( %s ) for amp %s ' % (all_met[7],ch16[iamp]) 
                    plt.ylabel(label)
        if show :plt.show()    
        rawPlotFile='BiasClusterAppliedToOverPar'
        SaveFig(fig,rawPlotFile,run_cur=run_cur,raft_cur=raft_cur,ccd_cur=ccd_cur)
        #
        # impact of a clusterization on overscan based on a given amp / method 
        if nb_amp==8  :
            n_sub=1
        else :
            n_sub=2
        nb_col=len(overpar[0,0,:])
        fig=plt.figure(figsize=[35,25])
        label='total number of images %d (black) .' % (nb_file)
        for icl in range(icluster[imet_max,iamp_max]):
            label='%s lot %d (%s) : %d images.' %(label,icl+1,color[icl],cluster_size[imet_max,iamp_max,icl])
        icl_increase=np.argsort(cluster_size[imet_max,iamp_max,:icluster[imet_max,iamp_max]])
        txt='run %s , RAFT %s CCD %s \n Overscan // by events clustered using amp %s (hdu=%d) with bias estimator: %s \n %s ' % (run_cur,raft_cur,ccd_cur,ch16[iamp_max],iamp_max+1,all_met[imet_max],label) 
        plt.suptitle(txt,fontsize=26)
        for iamp in range(nb_amp):
            label='amp %s <bias>' % (ch16[iamp])
            plt.subplot(2*n_sub,8,int(iamp/8)*8+iamp+1)
            plt.gca().set_title(label,fontsize=20)

            for icl in icl_increase:
                overparav=overpar[cluster[imet_max,iamp_max,icl,:cluster_size[imet_max,iamp_max,icl]],iamp,3:].mean(axis=0)
                label='data lot %d' % (icl+1)
                plt.plot(range(3,nb_col),overparav,color=color[icl],label=label)
                maxy=np.median(overparav)
                miny=maxy-30
                maxy=maxy+30
                mxy=np.max(overparav)+3
                miy=np.min(overparav)-3
                ylim=False
                if maxy< mxy : 
                    mx=maxy
                    ylim=True
                else :
                    mx=mxy
                if miny>  miy : 
                    mi=miny
                    ylim=True
                else :
                    mi=miy 
                if ylim : 
                    plt.ylim(mi,mx)
                if icl==0 and iamp%8==0 :
                    # first column , we can put the y label 
                    label='<// overscan> for amp %s ' % (ch16[iamp]) 
                    plt.ylabel(label,fontsize=18)
                if icl==0 and iamp+1==((n_sub-1)*8+4):
                    # last line , we can put the x label 
                    label='Column (pre-scan excluded) '  
                    plt.xlabel(label,fontsize=18)

                if iamp==iamp_max:
                    plt.legend()
            label='all data'
            plt.plot(range(3,nb_col),overpar[:,iamp,3:].mean(axis=0),color='black',label=label)
            #label='amp %s (hdu=%d), sig_bias' % (ch16[iamp],iamp+1)
            label='amp %s sig(bias)' % (ch16[iamp])
            plt.subplot(2*n_sub,8,int(iamp/8)*8+iamp+1+8)
            plt.gca().set_title(label,fontsize=20)
    #        for icl in range(icluster[imet_max,iamp_max]):
            for icl in icl_increase :
                overparst=overpar[cluster[imet_max,iamp_max,icl,:cluster_size[imet_max,iamp_max,icl]],iamp,3:].std(axis=0)
                label='data lot %d' % (icl+1)
                plt.plot(range(3,nb_col),overparst,color=color[icl],label=label)
                maxy=np.median(overparst)
                miny=0.
                maxy=maxy+10
                mxy=np.max(overparst)+2
                miy=np.min(overparst)-0.1
                ylim=False
                if maxy< mxy : 
                    mx=maxy
                    ylim=True
                else :
                    mx=mxy
                if miny>  miy : 
                    mi=miny
                    ylim=True
                else :
                    mi=miy 
                if ylim : 
                    plt.ylim(mi,mx)
                if iamp==iamp_max:
                    plt.legend()
                if icl==0 and iamp+1==((n_sub-1)*8+4):
                    # last line , we can put the x label 
                    label='Column (pre-scan excluded) '  
                    plt.xlabel(label,fontsize=18)
                if icl==0 and iamp%8==0 :
                    # first column , we can put the y label 
                    label='// overscan dispersion '
                    plt.ylabel(label,fontsize=18)
            label='all data'
            plt.plot(range(3,nb_col),overpar[:,iamp,3:].std(axis=0),color='black',label=label)
        #
        if show :plt.show()
        rawPlotFile='BiasClusterOverParDisper'
        SaveFig(fig,rawPlotFile,run_cur=run_cur,raft_cur=raft_cur,ccd_cur=ccd_cur)
    # prepare the return dictionary 
    to_return={}
    to_return['icluster']=icluster
    to_return['cluster']=cluster
    to_return['cluster_size']=cluster_size
    to_return['mean']=mean
    to_return['std']=std
    to_return['std_no_cluster']=std_no_cluster
    to_return['ref']=ref
    to_return['ref_mean']=ref_mean
    to_return['noise']=ampnoise

    return to_return
