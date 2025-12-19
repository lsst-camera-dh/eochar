import numpy as np
#
ch16=['C10','C11','C12','C13','C14','C15','C16','C17','C07','C06','C05','C04','C03','C02','C01','C00']
#
class ImageAna :
    def __init__(self,raw,verbose=True) :
        self.raw=raw
        self.det = raw.getDetector()
        self.nb_amp=len(self.det.getAmplifiers())
        self.vendor = self.det.getPhysicalType()
        det_ref=self.det['C10']
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
        self.amp_noise=np.zeros((self.nb_amp,3,2))
        if verbose : 
            print('CCD ',self.det_name)
            print('Number of lines read=',self.amp_y_size,' Number of colomns read=',self.amp_x_size)
            print('Number of lines in Image area=',self.im_y_size,' Number of colomns in Image area=',self.im_x_size)
            print('First line in Image area=',self.first_l,' First column in Image area=',self.first_c)
        # 
    def bias_cor(self,amp_cur,over_c='1D',over_l='1D',skip_c_over=2,skip_l_over=2,noise_analysis=False) :
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
        # Noise analysis
        if noise_analysis :
                # noise in most of the imageself. area per 5x5 bin after bias subtraction
                ny,nx=np.shape(self.image[index,self.first_l:self.first_l_over,self.first_c:self.first_c_over])
                nx5=int(nx/5)
                ny5=int(ny/5)
                nx0=int((nx-5*nx5)/2)+self.first_c
                ny0=int((ny-5*ny5)/2)+self.first_l
                b=self.image[index,ny0:5*ny5+ny0,nx0:5*nx5+nx0].reshape(ny5,5,nx5,5)
                b2=b**2
                E=((b.sum(axis=3).sum(axis=1))/25)
                EX2=((b2.sum(axis=3)).sum(axis=1))/25
                self.amp_noise[index,0,0]=np.median(E)
                self.amp_noise[index,0,1]=np.median(np.sqrt(EX2-E**2))
                # noise in  the serial overscan
                self.amp_noise[index,1,0]=np.median(image_2d[self.first_l+skip_l_over:,self.first_c_over+skip_c_over:].mean(axis=1))
                self.amp_noise[index,1,1]=np.median(image_2d[self.first_l+skip_l_over:,self.first_c_over+skip_c_over:].std(axis=1))
                # noise in  the // overscan 
                self.amp_noise[index,2,0]=np.median(image_2d[self.first_l_over+skip_l_over:,self.first_c+skip_c_over:].mean(axis=0))       
                self.amp_noise[index,2,1]=np.median(image_2d[self.first_l_over+skip_l_over:,self.first_c+skip_c_over:].std(axis=0))       
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
