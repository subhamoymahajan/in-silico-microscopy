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
import colorsys
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
    tif.imsave(outname,PSF)

def zstack_dat2tiff(filename, outname, lams, I0s, fs, ti, T, MaxBox, zmax, \
    opt_axis, dtype='uint16', add_z=1):
    """ Plots Generate a tiff file from multiple zstack images (volume image). 
        currently only works for monocolors.
    
    Looks for image data file 
    [filename][ti to ti+T]_[x/y/z][nidx]_lam[lam]_fs[fs].dat for ti>0
     or [filename]_[x/y/z][nidx]_lam[lam]_fs[fs].dat if ti < 0.

    Parameter
    ---------
    filename: str
        Filename header for image data file.
    outname: str
        Outname header for the tiff file. 
    lams: array of int
        The wavelength of all flurophore types
    I0s: array of floats
        The maximum image intensity of all fluorophore types
    fs: int
        Scaling factor for wave vector or MS position coordinates.
    ti: int
        timestep of the image data file. -1 if there is no sense of time
        (static system).
    T: int
        Number of timesteps to perform an average.
    MaxBox: array of ints
        Contains the number of pixels in the image in l and m directions.
    zmax: int
        Max number of zstacks. zmax is ideally length of the system in z (or n)
        direction, divided by the voxel lenght in the same direction.
        int(maxlen[opt_axis]/dlmn[2]+1E-3)
    opt_axis: int
        Optical axis 0, 1, and 2 for x, y, and z axis.
    dtype: str, optional
        Datatype of output tiff. (Default is 'uint16')
    add_z: int, optional    
        Adds 1 zstack for every [add_z] zstacks. (default is 1)
        
    Writes
    ------
    [outname]: A volume tiff image (mono color)
    """ 
    zN=int(zmax/add_z)
    xyz='xyz'
    for l in range(len(lams)):
        if dtype=='uint16':
            img=np.zeros((zN,MaxBox[1],MaxBox[0]),dtype=np.uint16)
            M=65535
        if dtype=='uint8':
            img=np.zeros((zN,MaxBox[1],MaxBox[0]),dtype=np.uint8)
            M=255
        Arguments=[]
        for z in range(zN):
            Arguments.append((filename,I0s[l],lams[l],T,ti,fs,MaxBox,\
                True, opt_axis,z*add_z))

        pool=mp.Pool(mp.cpu_count())
        results=pool.starmap(plot_image.get_grey_img,Arguments)
        for z in range(zN):
            foo=results[z]*M
            foo[foo>M]=M
            foo[foo<0]=0
            foo=np.transpose(foo)
            img[zN-z-1,:,:]=foo[:,:] 
        tif.imsave(outname+str(ti) + '_'+xyz[opt_axis] +'_lam'+str(lams[l])+ \
            '_fs'+str(fs) +'_T'+str(T) +'.tiff',img)
        
def tstack_dat2tiff(filename, outname, lams, I0s, fs, t0, tmax, tdiff, T, \
    MaxBox, dtype='uint16'):
    """ Plots Generate a tiff file from multiple time frames for the same 
        zstack (similar to a video). Currently only works for monocolors.
    
    Looks for image data file 

    [filename][ti to ti+T]_lam[lam]_fs[fs].dat for ti>0
    ti belongs to range(t0,tmax,tdiff)

    Parameter
    ---------
    filename: str
        Filename header for image data file.
    outname: str
        Outname header for the tiff file. 
    lams: array of int
        The wavelength of all flurophore types
    I0s: array of floats
        The maximum image intensity of all fluorophore types
    fs: int
        Scaling factor for wave vector or MS position coordinates.
    t0: int
        First timestep
    tmax: int
        Last timestep
    tdiff: int
        Timestep interval
    T: int
        Number of timesteps to perform an average.
    MaxBox: array of ints
        Contains the number of pixels in the image in l and m directions.
    dtype: str, optional
        Datatype of output tiff. (Default is 'uint16')
        
    Writes
    ------
    [outname]: A volume tiff image (mono color)
    """
    ts=np.arange(t0,tmax,tdiff) 
    for l in range(len(lams)):
        if dtype=='uint16':
            img=np.zeros((len(ts),MaxBox[1],MaxBox[0]),dtype=np.uint16)
            M=65535
        if dtype=='uint8':
            img=np.zeros((len(ts),MaxBox[1],MaxBox[0]),dtype=np.uint8)
            M=255
        Arguments=[]
        for t in ts:
            Arguments.append((filename,I0s[l],lams[l],T,t,fs,MaxBox,\
                True, None, None))
        pool=mp.Pool(mp.cpu_count())
        results=pool.starmap(plot_image.get_grey_img,Arguments)
        for i in range(len(ts)):
            foo=results[i]*M
            foo[foo>M]=M
            foo[foo<0]=0
            foo=np.transpose(foo)
            img[i,:,:]=foo[:,:] 
        tif.imsave(outname+'_tstack' + '_lam'+str(lams[l])+ \
            '_fs'+str(fs) +'_T'+str(T) +'.tiff',img)

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

        # for colorsys module hres should be belong to [0,1]
        foo=np.array(colorsys.hsv_to_rgb(hres,sres,vres))
        foo=[int(x*M) for x in foo]
        rgb=np.zeros(3,dtype=type(Is[0]))
        for i in range(3):
            if int(foo[i])>M:
                rgb[i]=M
            else:
                rgb[i]=int(foo[i])
        return rgb


def imgs2color(filename,outname,imgtype,mixtype,lam_hues,I0s=None):
    fnames=[]
    if os.path.exists(filename):
        if imgtype in ["tif1", "tiff1"]:
            IMGs=tif.imread(filename)
            dims=IMGs[0,:,:].shape #Assumes all shapes are the same
        elif imgtype in ["tif2", "tiff2"]:
            IMGs=tif.imread(filename)
            dims=IMGs[:,:,0].shape #Assumes all shapes are the same
        else:
            IMGs=cv2.imread(filename)
            dims=IMGs[:,:,0].shape #Assumes all shapes are the same
    else:
        i=0
        while True:
            if os.path.exists(filename+str(i)+'.'+imgtype):
                fnames.append(filename+str(i)+'.'+imgtype)
                i+=1
            else:
                break
        if imgtype in ["tiff1", "tif1"]:
            foo=tif.imread(fnames[0])
            dims=foo.shape #Assumes all shapes are the same
            IMGs=np.zeros((len(fnames),dims[0],dims[1]),dtype=type(foo[0,0]))
            for i in range(len(fnames)):
                IMGs[i,:,:]=tif.imread(fnames[i])
        elif imgtype in ["tiff2", "tif2"]:
            foo=tif.imread(fnames[0])
            dims=foo.shape #Assumes all shapes are the same
            IMGs=np.zeros((dims[0],dims[1],len(fnames)),dtype=type(foo[0,0]))
            for i in range(len(fnames)):
                IMGs[:,:,i]=tif.imread(fnames[i])
        else:
            foo=cv2.imread(fnames[0],0)
            dims=foo.shape #Assumes all shapes are the same
            IMGs=np.zeros((dims[0],dims[1],len(fnames)),dtype=type(foo[0,0]))
            for i in range(len(fnames)):
                IMGs[:,:,i]=cv2.imread(fnames[i],0)

    if mixtype=='rgb':
        cols=np.zeros((len(lam_hues),3))
        for i in range(len(lam_hues)):
            rgb=colorsys.hsv_to_rgb(lam_hues[i]/360,1,1)
            cols[i,:]=rgb[:]
            print('Hue '+str(lam_hues[i])+' changed to RGB '+str(cols[i,:]))

    print(dims)
    mixIMG=np.zeros((dims[0],dims[1],3)) #output is going to be png
    for i in range(dims[0]):
        for j in range(dims[1]):
            if imgtype in ["tiff1", "tif1"]:
                Is=IMGs[:,i,j]
            elif imgtype in ["tiff2", "tif2"]:
                Is=IMGs[i,j,:]
            else:
                Is=IMGs[i,j,:]
                Is=Is[::-1]
            if I0s!=None:
                Is=np.multiply(Is,I0s)
            if mixtype=='rgb':
                col=rgb_mix(Is,cols)
            elif mixtype=='mt':
                col=mt_mix(Is,lam_hues)
            mixIMG[i,j,:]=col[::-1]
    cv2.imwrite(outname+'.png',mixIMG)

           
#def monotiff_2color(filename,lams,hues):        
        
        

