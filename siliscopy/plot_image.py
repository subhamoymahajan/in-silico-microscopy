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
import colorsys
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import sys
small=1E-10

def get_grey_img(filename, I0, lam, T, ti, fs, MaxBox, whiteframe=False):
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
    whiteframe: Bool, optional
        True keeps the image intensity of white frame as -1. False converts
        the image intensities of -1 to 1. (default is False).
   
    Returns
    -------
    IMG: 2D ndarray
        Image intensities between 0 and 1. Image intensity of -1 implies white 
        frame (absence of molecular simulation system).
    """
    if ti<0 and T>1:
        raise Exception('More than one timesteps cannot be accomodate '+ \
                        'for a simulation without time.')

    IMG=np.zeros((MaxBox[0],MaxBox[1]))
    Cnt=np.zeros((MaxBox[0],MaxBox[1]),dtype=int)
    for i in range(T):
        if ti>=0:
            fname=filename + str(i+ti) + '_lam' + str(lam) + '_fs' + \
                  str(fs) + '.dat'
        else:
            fname=filename+'_lam'+str(lam)+'_fs'+str(fs)+'.dat'
        
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
                if whiteframe:
                    IMG[j,k]=-1
                else:
                    IMG[j,k]=1 #white background image
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
    show=False, gcolmap='gray', dpi=600):
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
    gcolmap: string
        Color map for monochrome image. (default is 'gray')
    dpi: int
        Dots per square inch for the image if saved to a file.

    Writes
    ------
    [filename]: JPEG Image file
        Writes the in-silico microscopy image if show is False and filename 
        is not None
    """
    consts=IMG.shape
    img_width=consts[1]*img_hei/consts[0]
    fontS=12
    if img_width<3:
        fontS=img_width*72/18.0
    
    fig,ax=plt.subplots(1, 1, figsize=(img_width, img_hei))
    if len(consts)==2: #Greyscale
        ax.imshow(IMG, vmin=0, vmax=1, cmap=gcolmap)
    else: #Color
        ax.imshow(IMG)
        Istring=''
        for i in range(len(lam_I0)):
            Istring+='_'+str(round(lam_I0[i],4))

    ax.set_xticks([])
    ax.set_yticks([])
    plt.axis('off')
    plt.tight_layout(pad=0)
    if filename==None or show==True:
        plt.show()
        return
    
    if ti>=0:#dynamic gro structure 
        if len(lam_I0)==1: #mono
            fname=filename + str(ti) + '_lam' + str(lam[0]) + '_fs' + \
                  str(fs) + '_T' + str(T) + '_I' + str(lam_I0[0]) + '.jpeg'
        else:#color
            fname=filename + str(ti) + '_fs' + str(fs) + '_T' + str(T) + \
                  '_I' + Istring + '.jpeg'
    else:# static gro structure
        if len(lam_I0)==1:#mono
            fname=filename + '_lam' + str(lam[0]) + '_fs' + str(fs) + '_I' + \
                  str(lam_I0[0]) + '.jpeg'
        else:#color
            fname=filename + '_fs' + str(fs) + '_T' + str(T) + '_I' + \
                  Istring + '.jpeg'

    print('Writing: '+fname)
    plt.savefig(fname,dpi=dpi,quality=100)
    plt.close()
    return

def get_col_img(filename, lam_I0s, lams, lam_hues, T, ti, fs, MaxBox):
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
   
    Returns
    ------- 
    IMG: 3D ndarray
        Axis 2 corresponds to red, green and blue channels. Image intensities 
        between 0 and 1.
    """
    IMGs=[]
    for i in range(len(lams)):
        IMGs.append(get_grey_img(filename, lam_I0s[i], lams[i], T, ti, fs,
                    MaxBox, whiteframe=True))
    consts=IMGs[0].shape
    col_IMG=np.zeros((consts[0],consts[1],3))
    for i in range(consts[0]):
        for j in range(consts[1]):
            foo=0
            for lam_id in range(len(lams)):
                if IMGs[lam_id][i,j]>-small:
                    foo=+1
            if foo==0:
                #if all Img_dat are -1 then res is -1
                col_IMG[i,j,:]=-1
            else:
                xres=0
                yres=0
                Is=[]
                ncol=0
                for lam_id in range(len(lams)):
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

                # for colorsys module hres should be belong to [0,1]
                rgb=list(colorsys.hsv_to_rgb(hres,sres,vres))
                col_IMG[i,j,:]=rgb[:]

    for i in range(consts[0]):
        for j in range(consts[1]):
            for k in range(3):
                if col_IMG[i,j,k]<-small:
                    col_IMG[i,j,k]=1.0
    return col_IMG    

def plot_grey_img(filename, lam_I0s, lams, T, ti, fs, MaxBox, Bm, scale, 
    dpi=600, outfile=None):
    """ Plots greyscale or monochrome image with a scale.
    
    See functions get_grey_img, add_scale and plot_ism for more details.
    """
    for i in range(len(lams)):
        IMG=get_grey_img(filename, lam_I0s[i], lams[i], T, ti, fs, MaxBox)
        IMG=add_scale(IMG, scale, Bm)
        plot_ism(IMG, [lam_I0s[i]], [lams[i]], T, ti, fs, filename=outfile,
                 dpi=dpi)

def plot_col_img(filename, lam_I0s, lams, lam_hues, T, ti, fs, MaxBox, Bm, 
    scale, dpi=600, outfile=None):
    """ Plots coloured image with a scale.
    
    See functions get_col_img, add_scale and plot_ism for more details.
    """
    IMG=get_col_img(filename, lam_I0s, lams, lam_hues, T, ti, fs, MaxBox)
    IMG=add_scale(IMG, scale, Bm)
    plot_ism(IMG,  lam_I0s,  lams,  T,  ti,  fs,  filename=outfile,  dpi=dpi)

def plot_grey_serial(filename, lam_I0s, lams, T, t0, tmax, tdiff, fs, MaxBox,
    Bm, scale, dpi, outname):
    """ Plots several greyscale or monochrome images serially.

    See plot_grey_img for more details.
    """
    for i in range(t0,tmax,tdiff):
        worker_grey([filename, lam_I0s, lams, T, i, fs, MaxBox, Bm, scale, dpi, 
                     outname])
    
def plot_col_serial(filename, lam_I0s, lams, lam_hues, T, t0, tmax, tdiff, fs, 
    MaxBox, Bm, scale, dpi, outname):
    """ Plots several coloured images serially.

    See plot_col_img for more details.
    """
    for i in range(t0,tmax,tdiff):
        worker_col([filename, lam_I0s, lams, lam_hues, T, i, fs, MaxBox, Bm, 
                    scale, dpi, outname])

def plot_grey_mp(filename, lam_I0s, lams, T, t0, tmax, tdiff, fs, MaxBox, Bm, 
    scale, dpi, output):
    """ Plots several greyscale or monochrome images parallelly.

    See plot_grey_img for more details.
    """

    Arguments=[]
    for i in range(t0,tmax,tdiff): 
        Arguments.append([filename, lam_I0s, lams, T, i, fs, MaxBox, Bm, scale, 
                          dpi, output])
    cpus=mp.cpu_count()
    if len(Arguments)<cpus:
        cpus=len(Arguments)
    pool=mp.Pool(cpus)
    results=pool.map(worker_grey, Arguments)

def plot_col_mp(filename, lam_I0s, lams, lam_hues, T, t0, tmax, tdiff, fs, 
    MaxBox, Bm, scale, dpi, output):
    """ Plots several coloured images parallely.

    See plot_col_img for more details.
    """
    Arguments=[]
    for i in range(t0,tmax,tdiff):
        Arguments.append([filename, lam_I0s, lams, lam_hues, T, i, fs, MaxBox, 
                          Bm, scale, dpi, output])
    cpus=mp.cpu_count()
    if len(Arguments)<cpus:
        cpus=len(Arguments)
    pool=mp.Pool(cpus)
    results=pool.map(worker_col, Arguments)

def worker_grey(Args):
    """ Runs plot_grey_img for a list of arguments.

    See plot_grey_img for more details
    """
    #             Args[0]  Args[1]  Args[2]  Args[3]  Args[4]  Args[5] 
    #            filename  lam_I0s    lams      T        ti       fs     
    plot_grey_img(Args[0], Args[1], Args[2], Args[3], Args[4], Args[5], 
                  Args[6], Args[7], Args[8], dpi=Args[9], outfile=Args[10])
               #  Args[6]  Args[7]  Args[8]    Args[9]      Args[10]
               #  MaxBox      Bm     scale      dpi          output
def worker_col(Args):
    """ Runs plot_col_img for a list of arguments.

    See plot_col_img for more details
    """
    #             Args[0] Args[1]  Args[2]  Args[3]  Args[4]  Args[5]  Args[6]
    #            filename lam_I0s    lams   lam_hues    T       ti       fs
    plot_col_img(Args[0], Args[1], Args[2], Args[3], Args[4], Args[5], Args[6],
                  Args[7], Args[8], Args[9], dpi=Args[10], outfile=Args[11])
               #  Args[7]  Args[8]  Args[9]    Args[10]       Args[11]
               #  MaxBox     Bm      scale       dpi           output

