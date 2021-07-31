#    Copyright 2021 SUBHAMOY MAHAJAN 
#    
#    This file is part of InSilicoMicroscopy software.
# 
#    InSilicoMicroscopy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.)
import numpy as np
import os
import multiprocessing as mp
import cv2
import tifffile as tif
from . import plot_image

def psf_dat2tiff(filename,outname,Plmn,dlmn,dtype='uint16'):
    """ Converts PSF.dat to PSF.tiff for its use in ImageJ and similar software.
        Works for 3D PSF.

    Parameters
    ----------
    filename: str
        Filename of the PSF file created in siliscopy format
    outname: str
        Output filename. Should end with .tif or .tiff.
    Plmn: array of floats
        Dimensions of the box withing which PSF is calculated.
    dlmn: array of floats
        Delta l', Delta m', and Delta n'. Voxel size.
    dtype: str
        Data type of the tiff.
    
    Writes
    ------
    [outname]: A Tiff file containing the PSF.
    """
    f=open(filename,'r')
    N=[int(Plmn[i]/dlmn[i]+1E-10)*2-1 for i in range(3)]
    if dtype=='uint16':
        PSF=np.zeros((N[2],N[1],N[0]),dtype=np.uint16)
    elif dtype=='uint8':
        PSF=np.zeros((N[2],N[1],N[0]),dtype=np.uint8)
    else:
        raise exception("Type error")

    x0,y0,z0=[int((a-1)/2) for a in N]
    
    for lines in f:
        if len(lines)==0:
            continue
        if lines[0]=='#':
            continue
        foo=lines.split()
        xyz=[float(foo[i])/dlmn[i] for i in range(3)]
        xyz_int=[int(a+1E-6) for a in xyz]
        for i in range(3):
            if abs(xyz[i]-xyz_int[i])>1E-6:
                continue #not part of the dlmn
        x,y,z=xyz_int
        if dtype=='uint16':
            I=int(float(foo[3])*65535) 
        elif dtype=='uint8':
            I=int(float(foo[3])*255) 

        for i in [-1,1]:
            for j in [-1,1]:
                for k in [-1,1]:
                    PSF[z0+k*z,x0+i*x,y0+j*y]=I
                    PSF[z0+k*z,y0+j*y,x0+i*x]=I
    tif.imsave(outname,PSF,resolution=(1./dlmn[0],1./dlmn[1]),imagej=True,
               metadata={'axes':'CYX', 'unit':'nm' })

def psf_dat2tiff2(filename,outname,Plmn,dlmn,dtype='uint16'):
    """ Converts PSF.dat to PSF.tiff. 0-0.2 is black to red, 0.2-0.8 is red hue 
        to blue hue (0 to 240 degrees), and 0.8 to 1 is blue to white. All other 
        properties are same as ''psf_dat2tiff''
       
    """
    f=open(filename,'r')
    N=[int(Plmn[i]/dlmn[i]+1E-10)*2-1 for i in range(3)]
    if dtype=='uint16':
        PSF=np.zeros((N[2],3,N[1],N[0]),dtype=np.uint16)
    elif dtype=='uint8':
        PSF=np.zeros((N[2],3,N[1],N[0]),dtype=np.uint8)
    else:
        raise exception("Type error")

    x0,y0,z0=[int((a-1)/2) for a in N]
    
    for lines in f:
        if len(lines)==0:
            continue
        if lines[0]=='#':
            continue
        foo=lines.split()
        xyz=[float(foo[i])/dlmn[i] for i in range(3)]
        xyz_int=[int(a+1E-6) for a in xyz]
        for i in range(3):
            if abs(xyz[i]-xyz_int[i])>1E-6:
                continue #not part of the dlmn
        x,y,z=xyz_int
        I=float(foo[3])
        if I<0.2:
            rgb=plot_image.hsv2rgb(0,1,I/0.2)            
        elif I<0.8:
            rgb=plot_image.hsv2rgb(2*(I-0.2)/1.8,1,1)
        else:
            rgb=plot_image.hsv2rgb(2./3.,(1-I)/0.2,1) #2/3 is blue

        if dtype=='uint16':
            col=(np.array(rgb)*65535).astype(np.uint16) 
        elif dtype=='uint8':
            col=(np.array(rgb)*255).astype(np.uint8)

        for i in [-1,1]:
            for j in [-1,1]:
                for k in [-1,1]:
                    PSF[z0+k*z,:,x0+i*x,y0+j*y]=col[:]
                    PSF[z0+k*z,:,y0+j*y,x0+i*x]=col[:]
    tif.imsave(outname,PSF,resolution=(1./dlmn[0],1./dlmn[1]),imagej=True,
               metadata={'spacing':dlmn[2], 'axes':'ZCYX', 'unit':'nm' })

def zstack2tiff(datafile, outname, dn):
    """ Generate a tiff file from multiple z-slice tiffs. 

    Parameter
    ---------
    datafile: str
        File containing n index and corresponding tiff file names
    outname: str
        Outname header for the tiff file. 
    dn: float
        Distance between two consecutive z-slice images.
    Writes
    ------
    [outname]: A volume tiff image
    """ 

    f=open(datafile,'r')
    data=[]
    for lines in f:
        foo=lines.strip()
        foo=foo.split(',')
        data.append(foo)
    nN=len(data)
    #Read Meta data.
    img_data=tif.TiffFile(data[0][1])
    res=(float(img_data.pages[0].tags['XResolution'].value[0]),
         float(img_data.pages[0].tags['YResolution'].value[0]))
    metaD={'spacing': dn, 'unit':'nm'}
    foo=img_data.pages[0].tags['ImageDescription']
    foo=foo.split('\n')
    foo=[ x.split() for x in foo]
    for i in range(len(foo)):
        if foo[i][0]=='finterval':
            metaD['finterval']=int(foo[i][1])

    axes=img_data.pages[0].axes
    sha=img_data.asarray().shape
    if 'T' in axes:
        metaD['axes']='TZ'+axes[1:]
        consts=tuple([sha[0],nN]+list(sha[1:]))
        ### Create image   
        IMG=np.zeros(consts,dtype=img0.dtype)
        if len(sha)==4: #TCYX
            for n in range(nN):
                IMG[:,int(data[n][0]),:,:,:]=tif.imread(data[n][1])
        elif len(sha)==3: #TYX
            for n in range(nN):
                IMG[:,int(data[n][0]),:,:]=tif.imread(data[n][1])

    else:
        metaD['axes']='Z'+axes
        consts=tuple([nN]+list(sha[1:]))
        ### Create image   
        IMG=np.zeros(consts,dtype=img0.dtype)
        for n in range(nN):
            IMG[int(data[n][0])]=tif.imread(data[n][1])

    tif.imwrite(outname, IMG, resolution=res, imagej=True, metadata=metaD)
        
def tstack2tiff(datafile, outname):
    """ Generate a tiff file from multiple time frames.
    
    Parameter
    ---------
    datafile: str
        File containing n index and corresponding tiff file names
    outname: str
        Outname header for the tiff file. 
        
    Writes
    ------
    [outname]: A timelapse tiff image.
    """
    f=open(datafile,'r')
    data=[]
    for lines in f:
        foo=lines.strip()
        foo=foo.split(',')
        data.append(foo)
    tN=len(data)


    img_data=tif.TiffFile(data[0][1])
    res=(float(img_data.pages[0].tags['XResolution'].value[0]),
         float(img_data.pages[0].tags['YResolution'].value[0]))
    metaD={'unit':'nm'}
    foo=img_data.pages[0].tags['ImageDescription']
    foo=foo.split('\n')
    foo=[ x.split() for x in foo]
    for i in range(len(foo)):
        if foo[i][0]=='spacing':
            metaD['spacing']=int(foo[i][1])
    axes=img_data.pages[0].axes
    metaD['axes']='T'+axes
      
    IMG=np.zeros(tuple([tN]+list(img_data.asarray().shape)),dtype=img0.dtype)
    for t in range(tN):
        IMG[int(data[t][0])]=tif.imread(data[n][1])

    tif.imwrite(outname, IMG, resolution=res, imagej=True, metadata=metaD)

def rgb_mix(Is,cols):
    if type(Is[0])==np.uint16:
        M=65535
    elif type(Is[0])==np.uint8:
        M=255
    else:
        raise exception("Data type of Is ("+str(type(Is[0]))+") not supported")

    foo=np.zeros(3)
    rgb=np.zeros(3,dtype=type(Is[0]))
    for i in range(len(Is)):
        foo+=Is[i]*cols[i,:]
    for i in range(3):
        if int(foo[i])>M:
            rgb[i]=M
        else:
            rgb[i]=int(foo[i])

    return rgb

def mt_mix(Is,lam_hues,small=1.9E-3):
    if type(Is[0])==np.uint16:
        M=65535
    elif type(Is[0])==np.uint8:
        M=255
    else:
        raise exception("Data type of Is ("+str(type(Is[0]))+") not supported")
    newIs=[x/M for x in Is]
    small=0.5/M

    XY=[0,0]
    Inonzero=[]
    ncol=0
    for i in range(len(Is)):
        if newIs[i]>small:
            # V*e^(i(hue)) #I0 was multiplied when determining 
            # grey images.
            XY[0]+=newIs[i]*np.cos(lam_hues[i]*np.pi/180)
            XY[1]+=newIs[i]*np.sin(lam_hues[i]*np.pi/180)
            Inonzero.append(newIs[i])
            ncol+=1
    if ncol==0:
        return np.zeros(3,dtype=type(Is))
    else:
        Inonzero=sorted(Inonzero)
        #arctan2 returns between [-pi,pi]. Dividing it by 2pi yields
        # [-0.5,0.5]
        hres=np.arctan2(XY[1],XY[0])/(2*np.pi)
        #this makes hres [0,1]
        if hres<0:
            hres+=1
        vres=Inonzero[-1] #Largest value
        if vres>1:
            raise Exception("Color's Value is more than 1!")
        sres=1
        if ncol>2: #Color should be saturated
            sres=1-Inonzero[-3]/Inonzero[-1]
        if sres<0 or sres>1:
            raise Exception("Color's saturation is more than 1")

        foo=np.array(plot_image.hsv2rgb(hres,sres,vres))
        foo=[int(x*M) for x in foo]
        rgb=np.zeros(3,dtype=type(Is[0]))
        for i in range(3):
            if int(foo[i])>M:
                rgb[i]=M
            else:
                rgb[i]=int(foo[i])
        return rgb


def imgs2color(datafile,outname,imgtype,mixtype,lam_hues,dpi=1200):

    f=open(datafile,'r')
    data=[]
    for lines in f:
        foo=lines.strip()
        foo=foo.split(',')
        data.append(foo)
    IMGs=[]
    cols=3
    if mix_type=='nomix':
        cols=len(data)

    print("Color mixing can lose precision!")
    if imgtype in ["tif", "tiff"]:
        #TZCYX / TCYX / ZCYX / CYX 
        img_data=tif.TiffFile(data[0][1])
        res=(float(img_data.pages[0].tags['XResolution'].value[0]),
             float(img_data.pages[0].tags['YResolution'].value[0]))
        metaD={'unit':'nm'}
        foo=img_data.pages[0].tags['ImageDescription']
        foo=foo.split('\n')
        foo=[ x.split() for x in foo]
        for i in range(len(foo)):
            if foo[i][0]=='spacing':
                metaD['spacing']=int(foo[i][1])
            if foo[i][0]=='finterval':
                metaD['spacing']=int(foo[i][1])
        axes=img_data.pages[0].axes
        sha=IMGs[0].shape
        if axes[0:2]=='TZ':
            metaD['axes']=axes[0:2]+'C'+axes[2:]
            new_sha=tuple(list(sha[0:2])+[cols]+list(sha[2:]))
        elif axes[0]=='T' or axes[0]=='Z':
            metaD['axes']=axes[0]+'C'+axes[1:]
            new_sha=tuple([sha[0],cols]+list(sha[1:]))
        else:
            metaD['axes']='C'+axes
            new_sha=tuple([cols]+list(sha[:]))

        for i in range(len(data)):
            foo=tif.imread(data[i][1])
            dtype=foo.dtype
            foo=foo.astype('float')
            foo=foo/255.0
            IMGs.append(foo) #Lams (TZYX / TYX / ZYX / YX) 
        IMGs=np.array(IMGs)
        if mixtype=='nomix':
            if axes=='TZYX':
                colIMG=np.transpose(IMGs,(1,2,0,3,4))
            elif axes=='TYX' or axes=='ZYX':
                colIMG=np.transpose(IMGs,(1,0,2,3))
        else:  #mt , rgb
            colIMG=np.zeros(new_sha)
            if axes=='TZYX':
                for t in range(new_sha[0]):
                    for z in range(new_sha[1]):
                        foo=np.transpose(IMGs[:,t,z,:,:],(2,1,0))
                        foo=plot_image.add_color(foo,lam_hues,0,mixtype) #XYC
                        foo=np.transpose(foo,(2,1,0)) #CYX
                        colIMG[t,z,:,:,:]=foo[:,:,:]
            elif axes=='TYX' or axes=='ZYX':
                for z in range(new_sha[1]):
                    foo=np.transpose(IMGs[:,z,:,:],(2,1,0))
                    foo=plot_image.add_color(foo,lam_hues,0,mixtype) #XYC
                    foo=np.transpose(foo,(2,1,0)) #CYX
                    colIMG[z,:,:,:]=foo[:,:,:]
            elif axes=='YX':
                foo=np.transpose(IMGs[:,:,:],(2,1,0))
                foo=plot_image.add_color(foo,lam_hues,0,mixtype) #XYC
                foo=np.transpose(foo,(2,1,0)) #CYX
                colIMG[:,:,:]=foo[:,:,:]
        tif.imwrite(outname,colIMG,resolution=res, imagej=True, metadata=metaD) 
    else:
        for i in range(len(data)):
            foo=cv2.imread(data[i][1])
            dtype=foo.dtype
            foo=foo.astype('float')
            foo=foo/255.0
            IMGs.append(foo)  #L XY
        IMGs=np.array(IMGs)
        IMGs=np.transpose(IMGs,(1,2,0)) #XYL
        consts=IMGs.shape
        print('White frame will not be calculated properly!')
        colIMG=plot_image.add_color(IMGs,lam_hues,frame_color=1.0,mix_type=mixtype)
        img_width=consts[1]*3/consts[0]

        fig,ax=plt.subplots(1, 1, figsize=(img_width, img_hei))
        ax.imshow(colIMG)
        ax.set_xticks([])
        ax.set_yticks([])
        plt.axis('off')
        plt.tight_layout(pad=0)
        print('Writing: '+outname)
        plt.savefig(outname,dpi=dpi)
        plt.close()
