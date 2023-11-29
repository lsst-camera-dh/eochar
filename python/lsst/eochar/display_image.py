import numpy as np 

def read_dim(a):
    yx=(a[1:-1]).split(',')
    x=np.array(yx[0].split(':')).astype(int)
    y=np.array(yx[1].split(':')).astype(int)
    return y,x

def SingleImageFITS(fits,image_to_display):
    # display à la ds9 shape 
    # fits : a pointer on a fits image/object of the same type than the one you want to display
    # image_to_display  :  a table of amp image  with [nb_amplifier] [:,:]  with the amp in the same order than in the fits file
    try : 
        y_size,x_size=read_dim(fits[1].header['DETSIZE'])
    except :
        y_size,x_size=read_dim(fits[0].header['DETSIZE'])
    nb_amp=len(image_to_display)
    if y_size[1]>4000: 
        print('Image for e2v sensor')
    else : 
        if nb_amp>8 :
            print('Image for itl sensor')
        else : 
            print('Image for itl corner raft sensor')
            # there is a bug in the  header ...fix it
            y_size[1]=int(y_size[1]/2)
    image=np.zeros((y_size[1],x_size[1]))
    #   
    for i in range(1,1+nb_amp) :
        y_datasec,x_datasec = read_dim(fits[i].header['DATASEC'])
        y_detsec,x_detsec   = read_dim(fits[i].header['DETSEC'])   
        if x_detsec[0]>x_detsec[1] :
            xstart=x_detsec[1]-1
            xend=x_detsec[0]
            xdstart=x_datasec[1]-1
            xdend=x_datasec[0]-2
            xsign=-1
        else : 
            xstart=x_detsec[0]-1
            xend=x_detsec[1]
            xdstart=x_datasec[0]-1
            xdend=x_datasec[1]
            xsign=1
        if y_detsec[0]>y_detsec[1] :
            ysign=-1
            ystart=y_detsec[1]-1
            yend=y_detsec[0]  
            ydstart=y_datasec[1]-1
            ydend=y_datasec[0]-2
        else : 
            ystart=y_detsec[0]-1
            yend=y_detsec[1]
            ydstart=y_datasec[0]-1
            ydend=y_datasec[1]
            ysign=1
        #print(ystart,yend,xstart,xend,ydstart,ydend,ysign,xdstart,xdend,xsign)
        if ydend<0 : 
            image[ystart:yend,xstart:xend]=np.copy(image_to_display[i-1][ydstart::ysign,xdstart:xdend:xsign]) 
        else : 
            image[ystart:yend,xstart:xend]=np.copy(image_to_display[i-1][ydstart:ydend:ysign,xdstart:xdend:xsign]) 
    return image
    
def SingleImageIR(image,is_e2v=True):
        # display à la ds9 shape 
        # Display an IR2 image , with amplifiers set at the right place ...there is a DM version which does this better...
        # but here you are in stand alone 
        # the default associated to the image area (pre-overscan excluded) are for e2v IR2 files 
        #
        if (is_e2v) :
            first_col=10
            first_cover=522	   
            first_line=0
            first_lower=2002
        else :
            first_col=3
            first_cover=512	   
            first_line=0
            first_lower=2000
        #
        col_size=first_cover-first_col
        line_size=first_lower-first_line
        #
        spf=np.zeros((line_size*2,col_size*8))
        if(is_e2v):
             for i in range(16) :
                  if i<8 :
                       xx=i*col_size-1
                       for x in range(first_col,first_cover) :  
                            for y in range(first_line,first_lower) :
                                 spf[2*line_size-(y-first_line)-1,xx+col_size-(x-first_col)]=image[i][y,x]
                  else :
                       xx=(15-i)*col_size
                       for y in range(first_line,first_lower) :
                            spf[y-first_line,xx:xx+col_size]=image[i][y,first_col:first_cover]
             return spf
        else:
             for i in range(16) :
                  if i<8 :
                       xx=i*col_size-1
                       for x in range(first_col,first_cover) :  
                            for y in range(first_line,first_lower) :
                                 spf[2*line_size-(y-first_line)-1,xx+col_size-(x-first_col)]=image[i][y,x]
                  else :
                       xx=(15-i)*col_size
                       for x in range(first_col,first_cover) :
                            for y in range(first_line,first_lower) : 
                                 spf[y-first_line,xx+col_size-(x-first_col)-1]=image[i][y,x]
        return spf
