#    Copyright 2020,2021 SUBHAMOY MAHAJAN 
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

# Creates a colored in-silico microscopy image from in-silico monochrome image 
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import sys
import copy
import tifffile as tif
import os
small=1E-10

def get_grey_img(filename, I0, lam, T, ti, fs, MaxBox, frame=False, \
    opt_axis=None, nidx=None, frame_col=1.0, noise=False, poi_a=None, 
    gauss_b=None, psf_type=0, tsO=None):
    """ Calculates greyscale image
     
    Parameters
    ----------
    filename: str
        Filename header for image data file
    I0: float
        The maximum image intensity 
    lam: int
        The wavelength
    T: int
        Number of timesteps to perform an average.
    ti: int
        timestep of the image data file. -1 if there is no sense of time
        (static system).
    fs: int
        Scaling factor for wave vector or MS position coordinates.
    MaxBox: array of ints
        Contains the number of pixels in the image in l and m directions.
    frame: Bool, optional
        True keeps the image intensity of white frame as -1. False converts
        the image intensities of -1 to 1. (default is False).
    opt_axis: int, optional
        Optical axis 0, 1, and 2 for x, y, and z axis. (default is None)
    nidx: int, optional    
        index of n' axis. (default is None)

    Returns
    -------
    IMG: 2D ndarray
        Image intensities between 0 and 1. Image intensity of -1 implies white 
        frame (absence of molecular simulation system).
    """
    if ti<0 and T>1:
        raise Exception('More than one timesteps cannot be used for a' + \
            ' simulation without time.')

    IMG=np.zeros((MaxBox[0],MaxBox[1]))
    Cnt=np.zeros((MaxBox[0],MaxBox[1]),dtype=int)
    nstr=''
    if nidx != None:
        xyz='xyz'
        nstr='_'+xyz[opt_axis]+str(nidx)
    tstr=''
    if psf_type==1:
        tstr='_tsO'+"%g"%tsO
    for i in range(T):
        if ti>=0:
            fname=filename + str(i+ti) + nstr + tstr + '_lam' + str(lam) + '_fs' + \
                  str(fs) + '.dat'
        else:
            fname=filename + nstr + tstr + '_lam' + str(lam) + '_fs' + str(fs) + '.dat'
        print('Reading: '+fname) 
        f=open(fname,'r')
        j=0
        for lines in f:
            foo=lines.split()
            if foo[0][0]=='#':
                continue
            for k in range(len(IMG[0])):
                I=float(foo[k])*I0
                if I>1:
                    IMG[j,k]+=1.0
                    Cnt[j,k]+=1
                elif I>= -small:
                    IMG[j,k]+=I
                    Cnt[j,k]+=1
                # else its white frame - do nothing
            j+=1                            
    for j in range(len(IMG)):
        for k in range(len(IMG[0])):
            if Cnt[j,k]>0:
                IMG[j,k]=IMG[j,k]/float(Cnt[j,k])
            else:
                if frame:
                    IMG[j,k]=-1
                else:
                    IMG[j,k]=frame_col
    if noise:
        IMG=add_noise(IMG,poi_a,gauss_b) 
    return IMG

def add_scale(IMG, scale, Bm):
    """ Adds a scale bar to the in-silico image
       
    Parameters
    ----------
    IMG: ndarray
        Image intensities
    scale: float
        Size of scale bar in nm.
    Bm:
        Size of maximum box length in m direction
 
    Returns
    -------
    IMG: ndarray
        Image intensities with a scale bar.
    """
    L=int(scale/Bm*len(IMG[0]))
    Llast=int(0.9*len(IMG[0]))
    Hlast=int(0.9*len(IMG))
    wid=int(len(IMG)*0.005)
    for i in range(Llast-L,Llast+1):
        for j in range(Hlast-wid,Hlast+wid+1):
            IMG[j][i]=1.0
    return IMG

def plot_ism(IMG, lam_I0, lam, T, ti, fs, img_hei=3.0, filename=None, 
    show=False, gcolmap='gray', dpi=600, otype='jpeg', psf_type=0, tsO=None):
    """ Plots in-silico microscopy image.
    
    Parameters
    ----------
    IMG: ndarray
        Image intensities
    lam_I0: array of float #change to lam_I0s
        An array containing maximum intensities of fluorophores
    lam: array of integers #change to lams
        An array of wavelength of fluorophores
    T: int
        Number of timesteps to perform an average.
    ti: int
        timestep of the image data file. -1 if there is no sense of time
        (static system).
    fs: int
        Scaling factor for wave vector or MS position coordinates.
    img_hei: float
        Physical image height in inches
    filename: str, optional
        Filename to save the image (default value is None)
    show: bool, optional
        Shows the image instead of saving if it is True. (default is False)
    gcolmap: string, optional
        Color map for monochrome image. (default is 'gray')
    dpi: int, optional
        Dots per square inch for the image if saved to a file.
    otype: str, optional
        Output type of image, 'jpeg','png' or 'tiff'. (default is 'jpeg')

    Writes
    ------
    [filename]: JPEG Image file
        Writes the in-silico microscopy image if show is False and filename 
        is not None
    """
    if filename!=None:
        Istring=''
        for i in range(len(lam_I0)):
            Istring+='_'+str(round(lam_I0[i],4))
        tstr=''
        if psf_type==1:
            tstr='_tsO'+"%g"%tsO
        if ti>=0:#dynamic gro structure 
            if len(lam_I0)==1: #mono
                fname=filename + str(ti) + tstr + '_lam' + str(lam[0]) + '_fs' + \
                      str(fs) + '_T' + str(T) + '_I' + str(lam_I0[0]) +'.'+otype
            else:#color
                fname=filename + str(ti) + tstr + '_fs' + str(fs) + '_T' + str(T) + \
                      '_I' + Istring +'.'+otype
        else:# static gro structure
            if len(lam_I0)==1:#mono
                fname=filename + tstr + '_lam' + str(lam[0]) + '_fs' + str(fs) + '_I' + \
                      str(lam_I0[0]) +'.'+otype
            else:#color
                fname=filename + tstr + '_fs' + str(fs) + '_T' + str(T) + '_I' + \
                      Istring + '.'+otype

    if otype in ["jpeg", "png", "JPEG", "PNG"]:
        consts=IMG.shape
        img_width=consts[1]*img_hei/consts[0]
        
        fig,ax=plt.subplots(1, 1, figsize=(img_width, img_hei))
        if len(consts)==2: #Greyscale
            ax.imshow(IMG, vmin=0, vmax=1, cmap=gcolmap)
        else: #Color
            ax.imshow(IMG)

        ax.set_xticks([])
        ax.set_yticks([])
        plt.axis('off')
        plt.tight_layout(pad=0)
        if filename==None or show==True:
            plt.show()
            return
        print('Writing: '+fname)
        plt.savefig(fname,dpi=dpi)
        plt.close()

    elif otype in ["tiff8", "tif8"]:
        IMG=IMG*255
        IMG=IMG.astype('uint8')
        if len(consts)==2: #Greyscale
            tif.imwrite(fname[:-1],IMG,photometric='minisblack')
        else: #Color
            tif.imwrite(fname[:-1],IMG,photometric='rgb',)

    elif otype in ["tiff16", "tif16"]:
        IMG=IMG*65535
        IMG=IMG.astype('uint16')
        if len(consts)==2: #Greyscale
            tif.imwrite(fname[:-1],IMG,photometric='minisblack')
        else: #Color
            tif.imwrite(fname[:-1],IMG,photometric='rgb')
  
    return

def get_col_img(filename, lam_I0s, lams, lam_hues, T, ti, fs, MaxBox, 
    frame_col=1.0, mix_type='mt', opt_axis=None, nidx=None, noise=False, 
    poi_a=None, gauss_b=None, psf_type=0, tsO=None):
    """ Calculates the image intensity for a color image,

    Parameters
    ----------
    filename: str
        Filename header for image data file
    lam_I0s: array of floats
        The maximum image intensity of all fluorophore types
    lams: array of int
        The wavelength of all flurophore types
    lam_hues: array of floats
        The hue in degree of all fluorophore types
    T: int
        Number of timesteps to perform an average.
    ti: int
        timestep of the image data file. -1 if there is no sense of time
        (static system).
    fs: int
        Scaling factor for wave vector or MS position coordinates.
    MaxBox: array of ints
        Contains the number of pixels in the image in l and m directions.
    frame_col= int, optional
        Color of the frame. (Default is 1 (white)).
    mix_type: str
        Algorithm to mix colors. 'mt' for Mahajan and Tang, 'rgb' for Red-
        Green-Blue. Addition is CMYK is not supported.

    Returns
    ------- 
    IMG: 3D ndarray
        Axis 2 corresponds to red, green and blue channels. Image intensities 
        between 0 and 1.
    """
    IMGs=[]
    for i in range(len(lams)):
        IMGs.append(get_grey_img(filename, lam_I0s[i], lams[i], T, ti, fs, MaxBox, 
            frame=True, opt_axis=opt_axis, nidx=nidx, noise=noise, poi_a=poi_a, 
            gauss_b=gauss_b, psf_type=psf_type, tsO=tsO))
    consts=IMGs[0].shape
    col_IMG=add_color(IMGs,lam_hues, frame_col, mix_type)

    return col_IMG  

def hsv2rgb(h, s, v):
    #h,s,v is between 0 to 1
    if h>1 or s>1 or v>1:
        raise Exception("h, s, or v > 1")

    A1=v #C+m
    A2=v*(1 - s*abs((h*6)%2-1))  #X+m
    A3=v*(1 - s) #m
    if h<1/6.:
        return A1,A2,A3
    elif h<2/6.:
        return A2,A1,A2
    elif h<3/6.:
        return A3,A1,A2
    elif h<4/6.:
        return A3,A2,A1
    elif h<5/6.:
        return A2,A3,A1
    else:
        return A1,A3,A2

def rgb2hsv(r, g, b):
    #r, g, b < 1
    if r>1 or g>1 or b>1:
        raise Exception("r, g, or b > 1")

    Cmax=max(r,g,b)
    Cmin=min(r,g,b)
    delta=Cmax-Cmin
    if delta<1E-6:
        h=0
    elif Cmax==r: 
        h=(((g-b)/delta)%6)/6.
    elif Cmax==g:
        h=(((b-r)/delta)%6)/6.
    elif Cmax==b:
        h=(((r-g)/delta)%6)/6.


    if Cmax==0:
        s=0
    else:
        s=delta/Cmax
   
    return h,s, Cmax 

def add_color(IMGs, lam_hues, frame_col=1.0, mix_type='mt'):
    """ Calculates the image intensity for a color image,

    Parameters
    ----------
    IMGs:list of 2D array
        Monochrome image intensities
    lam_hues: array of floats
        The hue in degree of all fluorophore types
    frame_col: float, optional
        Color of the frame, 0 is black, 1 is white, anything in between is 
        gray. (Default is 1.0)
    mix_type: str
        Algorithm to mix colors. 'mt' for Mahajan and Tang, 'rgb' for Red-
        Green-Blue. Addition is CMYK is not supported.

    Returns
    ------- 
    IMG: 3D ndarray
        Axis 2 corresponds to red, green and blue channels. Image intensities 
        between 0 and 1.
    frame_col= int, optional
        Color of the frame. (Default is 1 (white)).
    """
    consts=IMGs[0].shape
    col_IMG=np.zeros((consts[0],consts[1],3))
    if mix_type=='rgb':
        cols=np.zeros((len(lam_hues),3))
        for lam_id in range(len(lam_hues)):
            rgb=hsv2rgb(lam_hues[lam_id]/360,1,1)
            cols[lam_id,:]=rgb[:]

    for i in range(consts[0]):
        for j in range(consts[1]):
            foo=0
            for lam_id in range(len(lam_hues)):
                if IMGs[lam_id][i,j]>-small:
                    foo=+1
            if foo==0:
                #if all Img_dat are -1 then res is -1
                col_IMG[i,j,:]=-1
            elif mix_type=='mt':
                xres=0
                yres=0
                Is=[]
                ncol=0
                for lam_id in range(len(lam_hues)):
                    if IMGs[lam_id][i,j]>small:
                        # V*e^(i(hue)) #I0 was multiplied when determining 
                        # grey images.
                        xres+=IMGs[lam_id][i,j]*np.cos(lam_hues[lam_id]* \
                                                       2*np.pi/360)   
                        yres+=IMGs[lam_id][i,j]*np.sin(lam_hues[lam_id]* \
                                                       2*np.pi/360)
                        Is.append(IMGs[lam_id][i,j])
                        ncol+=1

                if ncol==0: #No fluorescence; black background (by default)
                    continue
                Is=sorted(Is) 
                #arctan2 returns between [-pi,pi]. Dividing it by 2pi yields
                # [-0.5,0.5]
                hres=np.arctan2(yres,xres)/(2*np.pi)
                #this makes hres [0,1]
                if hres<0:
                    hres+=1
                vres=Is[-1] #Largest value
                if vres>1:
                    raise Exception("Color's Value is more than 1!")
                sres=1
                if ncol>2: #Color should be saturated
                    sres=1-Is[-3]/Is[-1]
                if sres<0 or sres>1:
                    raise Exception("Color's saturation is more than 1")

                rgb=hsv2rgb(hres,sres,vres)
                col_IMG[i,j,:]=rgb[:]
            elif mix_type=='rgb':
                rgb=np.zeros(3)
                for lam_id in range(len(lam_hues)):
                    rgb+=IMGs[lam_id][i,j]*cols[lam_id,:]
                for k in range(3):
                    if rgb[k]>1:
                        col_IMG[i,j,k]=1.0 
                    else:
                        col_IMG[i,j,k]=rgb[k]

    for i in range(consts[0]):
        for j in range(consts[1]):
            for k in range(3):
                if col_IMG[i,j,k]<-small:
                    col_IMG[i,j,k]=frame_col

    return col_IMG  

def add_noise(I,poi_a,gauss_b):#I is single channel
    """ Adds a Poisson-Gaussian noise to a monochrome image

    Parameters
    ----------
    I: 2D array (image)
        Monochrome image
    poi_a: float
        Mean of the Poisson distribution. Should be a positive number.
    gauss_b: float
        Vairance of Gaussian distribution. Should be a positive number. Mean
        of Gaussian distribution is assumed to be 0.

    Returns
    -------
    I: Monochrome image with noise
    """
    dims=I.shape
    x0,y0,xL,yL=[0,0,0,0]
    for i in range(dims[0]):
        if I[i,int(dims[1]/2)-1]==-1 and xL==0:
            x0+=1
        elif I[i,int(dims[1]/2)-1]>-0.5:
            xL+=1       
    for j in range(dims[1]):
        if I[int(dims[0]/2)-1,j]==-1 and yL==0:
            y0+=1
        elif I[int(dims[0]/2)-1,j]>-0.5:
            yL+=1   
    for i in range(x0,x0+xL):
        for j in range(y0,y0+yL):
            P=np.random.poisson(lam=255*I[i,j])
            G=np.random.normal(loc=0.0,scale=gauss_b)
            I[i,j]=G+poi_a*P/255.0+I[i,j]*(1-poi_a)
            if I[i,j]>1:
                I[i,j]=1.0
            if I[i,j]<0: #Note that the space of i,j here is not the backgroumd
                I[i,j]=0.0
    #frame color will be taken care by add_color or separately for monochrome images

    return I

def plot_lumin(filename, lam_I0s, lams, lam_hues, T, ti, fs, MaxBox):
    """ Interactively plot small portions of colored the in-silico microscopy
        image, the relative luminescence, and hue as a funciton in pixel 
        distance from the central pixel.

    Parameters
    ----------
    filename: str
        Filename header for image data file
    lam_I0s: array of floats
        The maximum image intensity of all fluorophore types
    lams: array of int
        The wavelength of all flurophore types
    lam_hues: array of floats
        The hue in degree of all fluorophore types
    T: int
        Number of timesteps to perform an average.
    ti: int
        timestep of the image data file. -1 if there is no sense of time
        (static system).
    fs: int
        Scaling factor for wave vector or MS position coordinates.
    MaxBox: array of ints
        Contains the number of pixels in the image in l and m directions.
   
    Writes
    ------
        Ineteractively write the image generated in a custom filename.
    """
    col_IMGs=get_col_img(filename, lam_I0s, lams, lam_hues, T, ti, fs, MaxBox)
    consts=col_IMGs.shape
    img_hei=3.0
    img_width=consts[1]*img_hei/consts[0]
    fig,ax=plt.subplots(1, 1, figsize=(img_width, img_hei))
    ax.set_xticks([])
    ax.set_yticks([])
    plt.axis('off')
    plt.tight_layout(pad=0)
    ax.imshow(col_IMGs)
    plt.show(block=False)

    while True:
        print("Image dimensions are "+str(consts[0])+','+str(consts[1]))
        foo = input("Enter the pixel of interest (x,y): ")
        x,y=[int(x) for x in foo.split()]  
        width= int(input("Enter width in pixels: "))
        vals=[x,consts[0]-x,y,consts[1]-y]
        if min(vals)<width:
            width=min(vals)
        d=[]
        Lumin=[]
        Hue=[]
        for i in range(-width,width+1):
            for j in range(-width,width+1):
                d.append(np.sqrt(i**2+j**2)) 
                rgb=copy.deepcopy(col_IMGs[i+x,j+y,:])
                for k in range(3):
                   if rgb[k]<=0.03928:
                       rgb[k]=rgb[k]/12.92
                   else:
                       rgb[k]=((rgb[k]+0.055)/1.055)**2.4
                #Reference: https://www.w3.org/Graphics/Color/sRGB.html 
                #Check the published work for the accurate reference.
                hsv=rgb2hsv(rgb[0],rgb[1],rgb[2])
                Hue.append(hsv[0]*360)
                Lumin.append(0.2126*rgb[0]+0.7152*rgb[1]+0.0722*rgb[2])
        fig,ax=plt.subplots(1, 3, figsize=(img_width*3, img_hei))
        ax[0].set_xticks([])
        ax[0].set_yticks([])
        ax[0].axis('off')
        ax[0].imshow(col_IMGs[x-width:x+width+1,y-width:y+width+1])
        ax[1].scatter(d,Lumin,c='k')
        ax[1].set_ylabel('Relative Luminance')
        ax[1].set_xlabel('Pixel distance')
        ax[2].scatter(d,Hue,c='k')
        ax[2].set_ylabel('Hue (deg)')
        ax[2].set_xlabel('Pixel distance')
        plt.tight_layout()
        plt.show(block=False)
        save=input("Save figure, contiune or quit (s/c/q)?")
        if save=='s':
            outname=input("Save as: ")
            foo=input("ylim for Lumin: ")
            ylim1,ylim2 = [float(x) for x in foo.split()]  
            foo=input("ylim for Hue: ")
            ylim3,ylim4 = [float(x) for x in foo.split()]  
            ax[0].set_xticks([])
            ax[0].set_yticks([])
            ax[0].axis('off')
            ax[0].imshow(col_IMGs[x-width:x+width+1,y-width:y+width+1])
            ax[1].scatter(d,Lumin,c='k')
            ax[1].set_ylabel('Relative Luminance')
            ax[1].set_xlabel('Pixel distance')
            ax[1].set_ylim(ylim1,ylim2)
            ax[2].scatter(d,Hue,c='k')
            ax[2].set_ylabel('Hue (deg)')
            ax[2].set_xlabel('Pixel distance')
            ax[2].set_ylim(ylim3,ylim4)
            plt.tight_layout()
            plt.savefig(outname,dpi=600)
            return
        elif save=='q':
            return


def get_region(filename, lam_I0s, lams, lam_hues, T, ti, fs, MaxBox, \
    frame_col=1.0, mix_type='mt'):
    """ Calculates the image intensity for a region image,

    RGB color image is read. Colors are converted to HSV. Hues are changed to be 
    multiples of 10. Values are increased to 1, and the color is converted to 
    RGB.
     
    Parameters
    ----------
    filename: str
        Filename header for image data file
    lam_I0s: array of floats
        The maximum image intensity of all fluorophore types
    lams: array of int
        The wavelength of all flurophore types
    lam_hues: array of floats
        The hue in degree of all fluorophore types
    T: int
        Number of timesteps to perform an average.
    ti: int
        timestep of the image data file. -1 if there is no sense of time
        (static system).
    fs: int
        Scaling factor for wave vector or MS position coordinates.
    MaxBox: array of ints
        Contains the number of pixels in the image in l and m directions.
    frame_col= int, optional
        Color of the frame. (Default is 1 (white)).
    mix_type: str
        Algorithm to mix colors. 'mt' for Mahajan and Tang, 'rgb' for Red-
        Green-Blue. Addition is CMYK is not supported.
   
    Returns
    ------- 
    IMG: 3D ndarray
        Axis 2 corresponds to red, green and blue channels. Image intensities 
        between 0 and 1.
    """
    col_IMGs=get_col_img(filename, lam_I0s, lams, lam_hues, T, ti, fs, MaxBox, \
        frame_col=frame_col, mix_type=mix_type)
    consts=col_IMGs.shape
    for i in range(consts[0]):
        for j in range(consts[1]):
            hsv=rgb2hsv(col_IMGs[i,j,0],col_IMGs[i,j,1],col_IMGs[i,j,2])
            hsv[0]=int(hsv[0]*36)/36
            rgb=hsv2rgb(hsv[0],hsv[1],1)
            col_IMGs[i,j,:]=rgb[:]
    return col_IMGs   

def plot_region(filename, lam_I0s, lams, lam_hues, T, ti, fs, MaxBox, Bm,
    scale, dpi=600, outfile=None):
    """ Plots color image by assigning a specific color for a region.
   
    Parameters
    ----------
    filename: str
        Filename header for image data file
    lam_I0s: array of floats
        The maximum image intensity of all fluorophore types
    lams: array of int
        The wavelength of all flurophore types
    lam_hues: array of floats
        The hue in degree of all fluorophore types
    T: int
        Number of timesteps to perform an average.
    ti: int
        timestep of the image data file. -1 if there is no sense of time
        (static system).
    fs: int
        Scaling factor for wave vector or MS position coordinates.
    MaxBox: array of ints
        Contains the number of pixels in the image in l and m directions.
    Bm: float
        Width of the image in nm.
    scale: float
        Length of scale bar in nm
    dpi: int
        Dots per square inch of output image.
    outfile: str
        Output filename string
  
    Writes
    ------
    Writes a "region" image of in-silico microscopy image.
    """
    IMG=get_region(filename, lam_I0s, lams, lam_hues, T, ti, fs, MaxBox)
    IMG=add_scale(IMG, scale, Bm)
    if outfile!=None:
        plot_ism(IMG, lam_I0s, lams, T, ti, fs, filename='Reg_'+outfile, \
                 dpi=dpi)
    else:
        plot_ism(IMG, lam_I0s, lams, T, ti, fs, filename=None, dpi=dpi)


def plot_grey_img(filename, lam_I0s, lams, T, ti, fs, MaxBox, Bm, scale, 
    dpi=600, outfile=None, frame_col=1.0, noise=False, poi_a=None, gauss_b=None,
    otype='jpeg', psf_type=0, tsO=None):
    """ Plots monochrome image with a scale.
    
    See functions get_grey_img, add_scale and plot_ism for more details.
    """
    for i in range(len(lams)):
        IMG=get_grey_img(filename, lam_I0s[i], lams[i], T, ti, fs, MaxBox, 
            frame_col=frame_col, noise=noise, poi_a=poi_a, gauss_b=gauss_b, 
            psf_type=psf_type, tsO=tsO)
        IMG=add_scale(IMG, scale, Bm)
        plot_ism(IMG, [lam_I0s[i]], [lams[i]], T, ti, fs, filename=outfile,
                 dpi=dpi, otype=otype, psf_type=psf_type, tsO=tsO)

def plot_col_img(filename, lam_I0s, lams, lam_hues, T, ti, fs, MaxBox, Bm, 
    scale, dpi=600, outfile=None, frame_col=1.0, mix_type='mt', noise=False,
    poi_a=None, gauss_b=None, otype='jpeg',psf_type=0, tsO=None):
    """ Plots coloured image with a scale.
    
    See functions get_col_img, add_scale and plot_ism for more details.
    """
    IMG=get_col_img(filename, lam_I0s, lams, lam_hues, T, ti, fs, MaxBox, 
        frame_col, mix_type, noise=noise, poi_a=poi_a, gauss_b=gauss_b,
        psf_type=psf_type, tsO=tsO)
    IMG=add_scale(IMG, scale, Bm)
    plot_ism(IMG,  lam_I0s,  lams,  T,  ti,  fs,  filename=outfile,  dpi=dpi, 
        otype=otype, psf_type=psf_type, tsO=tsO)

def get_grey_3dimg(filename, lam_I0, lam, T, ti, fs, MaxBox, dlmn, nmax,
    opt_axis, add_n=1, outfile=None, frame_col=1.0, otype='tiff8',
    mprocess=True, noise=False, poi_a=None, gauss_b=None, psf_type=0, 
    tsO=None):

    nN=int(nmax/add_n)
    xyz='xyz'
    if otype in ["tiff16","tif16"]:
        img3d=np.zeros((nN,MaxBox[1],MaxBox[0]),dtype=np.uint16)
        M=65535
    elif otype in ["tiff8", "tif8"]:
        img3d=np.zeros((nN,MaxBox[1],MaxBox[0]),dtype=np.uint8)
        M=255
    if mprocess==True: #Multiprocess
        Arguments=[]
        for n in range(nN):
            Arguments.append((filename, lam_I0, lam, T, ti, fs, MaxBox, True, \
                opt_axis, n*add_n, -1, noise, poi_a, gauss_b, psf_type, tsO))
        pool=mp.Pool(mp.cpu_count())
        results=pool.starmap(get_grey_img,Arguments)
    else: #Serial
        results=[]
        for n in range(nN):
            results.append(get_grey_img(filename, lam_I0, lam, T, ti, fs,
                MaxBox, True, opt_axis, n*add_n, -1, noise, poi_a, gauss_b,
                psf_type, tsO))

    for z in range(zN):
        foo=results[z]*M
        foo[foo>M]=M
        foo[foo<0]=int(frame_col*M)
        foo=np.transpose(foo)
        img3d[zN-z-1,:,:]=foo[:,:]

    return img3d

def plot_grey_3dimg(filename, lam_I0s, lams, T, ti, fs, MaxBox, dlmn, nmax,
    opt_axis, add_n=1, outfile=None, frame_col=1.0, otype='tiff8', 
    mprocess=True, noise=False, poi_a=None, gauss_b=None, psf_type=0, tsO=None):

    tstr=''
    if psf_type==1:
        tstr='_tsO'+"%g"%tsO
    for l in range(len(lams)):
        img3d=get_grey_3dimg(filename, lam_I0s[l], lams[l], T, ti, fs, MaxBox,
            dlmn, nmax, opt_axis, add_n=add_n, outfile=outfile, otype=otype,
            frame_col=frame_col, mprocess=mprocess, noise=noise, poi_a=poi_a,
            gauss_b=gauss_b, psf_type=psf_type, tsO=tsO)
        tif.imwrite(outname+str(ti) + '_'+xyz[opt_axis] + tstr + '_lam' + \
            str(lams[l]) + '_fs'+str(fs) +'_T'+str(T) +'I_'+str(lam_I0) + \
            '.tiff', img3d, imagej=True, resolution=(1./dlmn[0],1./dlmn[1]), \
             metadata={'spacing': dlmn[2], 'unit': 'nm', 'axes': 'ZYX'})
    
def plot_grey_3dtimg(filename, lam_I0s, lams, T, tbegin, tmax, tdiff, fs, MaxBox, 
    dlmn, nmax, opt_axis, fps, add_n=1, outfile=None, frame_col=1.0, 
    otype='tiff8', mprocess=True, noise=False, poi_a=None, gauss_b=None, 
    psf_type=0, tsO=None):
    
    nN=int(nmax/add_n)
    tN=int((tmax-1-tbegin)/tdiff)+1
    xyz='xyz'
    if otype in ["tiff16","tif16"]:
        img3dt=np.zeros((tN,nN,MaxBox[1],MaxBox[0]),dtype=np.uint16)
        M=65535
    elif otype in ["tiff8", "tif8"]:
        img3dt=np.zeros((tN,nN,MaxBox[1],MaxBox[0]),dtype=np.uint8)
        M=255

    tstr=''
    if psf_type==1:
        tstr='_tsO'+"%g"%tsO
    #Check if all tiff exists
    for l in range(len(lams)):
        for i in range(tN):
            ti=tbegin+tdiff*i
            fname=filename+str(ti) + '_'+xyz[opt_axis] +tstr +'_lam' + \
                str(lams[l]) +'_fs'+str(fs) +'_T'+str(T) +'I_'+str(lam_I0) + \
                '.tiff'
            if os.path.exists(fname):
                img3d=tif.imread(fname)
            else:
                img3d=get_grey_3dimg(filename, lam_I0s[l], lams[l], T, ti, fs, 
                    MaxBox, dlmn, nmax, opt_axis, add_n=add_n, outfile=outfile, 
                    otype=otype, frame_col=frame_col, mprocess=mprocess, 
                    noise=noise, poi_a=poi_a, gauss_b=gauss_b, 
                    psf_type=psf_type, tsO=tsO)
            img3dt[i,:,:,:]=img3d[:,:,:]
        tif.imwrite(outname+'_' +xyz[opt_axis] +tstr +'_lam'+str(lams[l]) + \
            '_fs'+str(fs) +'_T'+str(T) +'I_'+str(lam_I0) +'.tiff', img3dt, 
            imagej=True, resolution=(1./dlmn[0],1./dlmn[1]), metadata={ 
            'spacing': dlmn[2], 'unit': 'nm', 'finterval': fps, 'axes': 'TZYX'})

def plot_grey_2dtimg(filename, lam_I0s, lams, T, tbegin, tmax, tdiff, fs, MaxBox, 
    dlmn, nmax, opt_axis, fps, add_n=1, outfile=None, frame_col=1.0, 
    otype='tiff8', mprocess=True, noise=False, poi_a=None, gauss_b=None, 
    psf_type=0, tsO=None):
    
    tN=int((tmax-1-tbegin)/tdiff)+1
    if otype in ["tiff16","tif16"]:
        img3dt=np.zeros((tN,MaxBox[1],MaxBox[0]),dtype=np.uint16)
        M=65535
    elif otype in ["tiff8", "tif8"]:
        img3dt=np.zeros((tN,MaxBox[1],MaxBox[0]),dtype=np.uint8)
        M=255
    tstr=''
    if psf_type==1:
        tstr='_tsO'+"%g"%tsO
    #Check if all tiff exists
    for l in range(len(lams)):
        for i in range(tN):
            ti=tbegin+tdiff*i
            fname=filename+str(ti) + '_lam'+str(lams[l]) +tstr +'_fs'+str(fs) +\
                '_T'+str(T) +'I_'+str(lam_I0) +'.tiff'
            if os.path.exists(fname):
                img2d=tif.imread(fname)
            else:
                img2d=get_grey_img(filename, lam_I0s[l], lams[l], T, ti, fs, 
                    MaxBox, frame=False, opt_axis=None, nidx=None, 
                    frame_col=frame_col, noise=noise, poi_a=poi_a, 
                    gauss_b=gauss_b, psf_type=psf_type, tsO=tsO)
                img2d=np.transpose(img2d)
            img2dt[i,:,:]=img2d[:,:]
        tif.imwrite(outname+ '_lam'+str(lams[l]) +tstr +'_fs'+str(fs) +'_T' + \
            str(T) +'I_'+str(lam_I0) +'.tiff', img2dt, imagej=True, 
            resolution=(1./dlmn[0],1./dlmn[1]), metadata={'unit': 'nm', 
            'finterval': fps, 'axes': 'TYX'})

def get_col_3dimg(filename, lam_I0s, lams, lam_hues, T, ti, fs, MaxBox, dlmn, 
    nmax, opt_axis, add_n=1, outfile=None, frame_col=1.0, otype='tiff8', 
    mprocess=True, noise=False, poi_a=None, gauss_b=None, mix_type='mt',
    mix=True, psf_type=0, tsO=None):

    nN=int(nmax/add_n)
    xyz='xyz'
    C=3 #3 color channels for 'mt' image
    if mix==False: #Colors are not mixed.
        C=len(lams)
    if otype in ["tiff16","tif16"]:
        img3d=np.zeros((nN,C,MaxBox[1],MaxBox[0]),dtype=np.uint16)
        M=65535
    elif otype in ["tiff8", "tif8"]:
        img3d=np.zeros((nN,C,MaxBox[1],MaxBox[0]),dtype=np.uint8)
        M=255

    if mix:
        if mprocess==True: #Multiprocess
            Arguments=[]
            for n in range(nN):
                Arguments.append((filename, lam_I0s, lams, lam_hues, T, ti, fs, 
                    MaxBox, frame_col, mix_type, opt_axis, n*add_n, noise, poi_a, 
                    gauss_b, psf_type, tsO))
            pool=mp.Pool(mp.cpu_count())
            results=pool.starmap(get_col_img,Arguments)
        else: #Serial
            results=[]
            for n in range(nN):
                results.append(get_col_img(filename, lam_I0s, lams, lam_hues, 
                    T, ti, fs, MaxBox, frame_col, mix_type, opt_axis, n*add_n, 
                    noise, poi_a, gauss_b, psf_type, tsO))
        
        for z in range(zN):
            foo=results[z]*M
            foo[foo>M]=M
            foo[foo<0]=int(frame_col*M) #foo has x, y, c
            foo=np.transpose(foo,(2,1,0)) # not foo has c, y, x
            img3d[zN-z-1,:,:,:]=foo[:,:,:] #TODO: Check if order is correct.

    else:
        for l in range(len(lams)):
            if mprocess==True: #Multiprocess
                Arguments=[]
                for n in range(nN):
                    Arguments.append((filename, lam_I0s[l], lams[l], T, ti, fs, 
                        MaxBox, True, opt_axis, n*add_n, -1, noise, poi_a, 
                        gauss_b, psf_type, tsO))
                pool=mp.Pool(mp.cpu_count())
                results=pool.starmap(get_grey_img,Arguments)
            else: #Serial
                results=[]
                for n in range(nN):
                    results.append(get_grey_img(filename, lam_I0s[l], lams[l], 
                        T, ti, fs, MaxBox, True, opt_axis, n*add_n, -1, noise,
                        poi_a, gauss_b, psf_type, tsO))
        
            for z in range(zN):
                foo=results[z]*M
                foo[foo>M]=M
                foo[foo<0]=int(frame_col*M)
                foo=np.transpose(foo)
                img3d[zN-z-1,l,:,:]=foo[:,:] #TODO: Check if order is correct.

    return img3d

def plot_col_3dimg(filename, lam_I0s, lams, lam_hues, T, ti, fs, MaxBox, dlmn, 
    nmax, opt_axis, add_n=1, outfile=None, frame_col=1.0, otype='tiff8', 
    mprocess=True, noise=False, poi_a=None, gauss_b=None, mix_type='mt',
    mix=True, psf_type=0, tsO=None):
    
    Istring=''
    for i in range(len(lams)):
        Istring+='_'+lam_I0s[i]
    tstr=''
    if psf_type==1:
        tstr='_tsO'+"%g"%tsO
    img3d=get_col_3dimg(filename, lam_I0s, lams, lam_hues, T, ti, fs, MaxBox,
        dlmn, nmax, opt_axis, add_n=add_n, outfile=outfile, frame_col=frame_col,
        otype=otype, mprocess=mprocess, noise=noise, poi_a=poi_a, 
        gauss_b=gauss_b, mix_type=mix_type, mix=mix, psf_type=psf_type, tsO=tsO)
    tif.imwrite(outname+str(ti) + '_'+xyz[opt_axis] +tstr +'_fs'+str(fs) +'_T'+\
        str(T) +'I'+Istring+'.tiff', img3d, imagej=True, 
        resolution=(1./dlmn[0],1./dlmn[1]), metadata={'spacing': dlmn[2], 
        'unit': 'nm', 'axes': 'ZCYX'})

def plot_col_3dtimg(filename, lam_I0s, lams, lam_hues, T, tbegin, tmax, tdiff, fs, 
    MaxBox, dlmn, nmax, opt_axis, fps, add_n=1, outfile=None, frame_col=1.0, 
    otype='tiff8', mprocess=True, noise=False, poi_a=None, gauss_b=None, 
    mix_type='mt', mix=True, psf_type=0, tsO=None):
    
    nN=int(nmax/add_n)
    tN=int((tmax-1-tbegin)/tdiff)+1
    C=3 #3 color channels for 'mt' image
    if mix==False: #Colors are not mixed.
        C=len(lams)
    
    xyz='xyz'
    if otype in ["tiff16","tif16"]:
        img3dt=np.zeros((tN,nN,C,MaxBox[1],MaxBox[0]),dtype=np.uint16)
        M=65535
    elif otype in ["tiff8", "tif8"]:
        img3dt=np.zeros((tN,nN,C,MaxBox[1],MaxBox[0]),dtype=np.uint8)
        M=255
    Istring=''
    for i in range(len(lams)):
        Istring+='_'+str(lam_hues[i])
    tstr=''
    if psf_type==1:  
        tstr='_tsO'+"%g"%tsO
    #Check if all tiff exists
    for i in range(tN):
        ti=tbegin+tdiff*i
        fname=filename+str(ti) + '_'+xyz[opt_axis] +tstr +'_fs'+str(fs) + \
            '_T'+str(T) +'I'+Istring +'.tiff'
        if os.path.exists(fname):
            img3d=tif.imread(fname)
        else:
            img3d=get_col_3dimg(filename, lam_I0s, lams, lam_hues, T, ti, fs,
                MaxBox, dlmn, nmax, opt_axis, add_n=add_n, outfile=outfile, 
                frame_col=frame_col, otype=otype, mprocess=mprocess, mix=mix,
                noise=noise, poi_a=poi_a, gauss_b=gauss_b, mix_type=mix_type, 
                psf_type=psf_type, tsO=tsO)
        img3dt[i,:,:,:,:]=img3d[:,:,:,:]
    tif.imwrite(outname+'_' +xyz[opt_axis] +tstr +'_fs'+str(fs) +'_T'+str(T)+ \
        'I'+Istring +'.tiff', img3dt, resolution=(1./dlmn[0],1./dlmn[1]), 
        imagej=True, metadata={'spacing': dlmn[2], 'unit': 'nm', 
        'finterval': fps, 'axes': 'TZCYX'})

def plot_col_2dtimg(filename, lam_I0s, lams, lam_hues, T, tbegin, tmax, tdiff, fs, 
    MaxBox, fps, outfile=None, frame_col=1.0, otype='tiff8', mprocess=True, 
    noise=False, poi_a=None, gauss_b=None, mix_type='mt', mix=True, psf_type=0,
    tsO=None):
    
    tN=int((tmax-1-tbegin)/tdiff)+1
    C=3 #3 color channels for 'mt' image
    if mix==False: #Colors are not mixed.
        C=len(lams)
    
    xyz='xyz'
    if otype in ["tiff16","tif16"]:
        img2dt=np.zeros((tN,C,MaxBox[1],MaxBox[0]),dtype=np.uint16)
        M=65535
    elif otype in ["tiff8", "tif8"]:
        img2dt=np.zeros((tN,C,MaxBox[1],MaxBox[0]),dtype=np.uint8)
        M=255
    Istring=''
    for i in range(len(lams)):
        Istring+='_'+str(lam_hues[i])
    tstr=''
    if psf_type==1:
        tstr='_tsO'+"%g"%tsO
    #Check if all tiff exists
    for i in range(tN):
        ti=tbegin+tdiff*i
        fname=filename+str(ti) + tstr + '_fs'+str(fs) +'_T'+str(T)+\
            'I'+Istring +'.tiff'
        if os.path.exists(fname):
            img2d=tif.imread(fname)
            img3dt[i,:,:,:]=img2d[:,:,:]
        elif mix:
            img2d=get_col_img(filename, lam_I0s, lams, lam_hues, T, ti, fs,
                MaxBox, frame_col=frame_col, mix_type=mix_type, opt_axis=None,
                nidx=None, noise=noise, poi_a=poi_a, gauss_b=gauss_b,
                psf_type=psf_type, tsO=tsO)
            img2d=np.transpose(img2d,(2,1,0))    
            img3dt[i,:,:,:]=img2d[:,:,:]
        else:
            for j in range(C):
                img2d=get_grey_img(filename, lam_I0s[j], lams[j], T, ti, fs, 
                    MaxBox, frame=False, opt_axis=None, nidx=None, 
                    frame_col=frame_col, noise=noise, poi_a=poi_a, 
                    gauss_b=gauss_b, psf_type=psf_type, tsO=tsO)
                img2d=np.tranpose(img2d)
                img3dt[i,j,:,:]=img2d[:,:]
    tif.imwrite(outname+'_' +xyz[opt_axis] + tstr + '_fs'+str(fs) +'_T' + \
        str(T) +'I'+Istring +'.tiff', img3dt, 
        resolution=(1./dlmn[0],1./dlmn[1]), imagej=True, 
        metadata={'spacing': dlmn[2], 'unit': 'nm', 'finterval': fps, 
        'axes': 'TZCYX'})

def plot_grey_serial(filename, lam_I0s, lams, T, tbegin, tmax, tdiff, fs, MaxBox,
    Bm, scale, dpi, outname, frame_col, noise, poi_a, gauss_b, otype, psf_type, 
    tsO):
    """ Plots several monochrome images serially.

    See plot_grey_img for more details.
    """
    for i in range(tbegin,tmax,tdiff):
        plot_grey_img(filename, lam_I0s, lams, T, i, fs, MaxBox, Bm, scale, dpi, 
                      outname, frame_col, noise, poi_a, gauss_b, otype, psf_type,
                      tsO)

    
def plot_col_serial(filename, lam_I0s, lams, lam_hues, T, tbegin, tmax, tdiff, fs, 
    MaxBox, Bm, scale, dpi, outname, frame_col, mix_type, noise, poi_a, gauss_b, 
    otype, psf_type, tsO):
    """ Plots several coloured images serially.

    See plot_col_img for more details.
    """
    for i in range(tbegin,tmax,tdiff):
        plot_col_img(filename, lam_I0s, lams, lam_hues, T, i, fs, MaxBox, Bm, 
            scale, dpi, outname, frame_col, mix_type, noise, poi_a, gauss_b, 
            otype, psf_type, tsO)

def plot_grey_mp(filename, lam_I0s, lams, T, tbegin, tmax, tdiff, fs, MaxBox, Bm, 
    scale, dpi, output, frame_col, noise, poi_a, gauss_b, otype, psf_type, tsO):
    """ Plots several monochrome images parallelly.

    See plot_grey_img for more details.
    """

    Arguments=[]
    for i in range(tbegin,tmax,tdiff): 
        Arguments.append([filename, lam_I0s, lams, T, i, fs, MaxBox, Bm, scale, 
                          dpi, output, frame_col, noise, poi_a, gauss_b, otype,
                          psf_type, tsO])
    cpus=mp.cpu_count()
    if len(Arguments)<cpus:
        cpus=len(Arguments)
    pool=mp.Pool(cpus)
    results=pool.starmap(plot_grey_img,Arguments)

def plot_col_mp(filename, lam_I0s, lams, lam_hues, T, tbegin, tmax, tdiff, fs, 
    MaxBox, Bm, scale, dpi, output, frame_col, mix_type, noise, poi_a, gauss_b,
    otype, psf_type, tsO):
    """ Plots several coloured images parallely.

    See plot_col_img for more details.
    """
    Arguments=[]
    for i in range(tbegin,tmax,tdiff):
        Arguments.append([filename, lam_I0s, lams, lam_hues, T, i, fs, MaxBox, 
            Bm, scale, dpi, output,frame_col,mix_type, noise, poi_a, gauss_b, 
            otype,psf_type,tsO])
    cpus=mp.cpu_count()
    if len(Arguments)<cpus:
        cpus=len(Arguments)
    pool=mp.Pool(cpus)
    results=pool.starmap(plot_col_img,Arguments)
