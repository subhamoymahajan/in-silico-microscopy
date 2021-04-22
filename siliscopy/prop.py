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
import numpy as np
from .plot_image import get_grey_img
def get_maxI(filename, lams, fs):
    """ Get maximum intensity from I(l',m'). 

    Parameters
    ----------
    filename: str
        Starting string of filename of image intensity files.
    lams: array of int
        Wavelength of all fluorophore types.
    fs: int
        Full-Width-at-Half-Maximum (FWHM) scaling factor.

    Returns 
    -------
    Imax: array of float
        Maximum intensity for each fluorophore type.
    """
    Imax=np.zeros(len(lams))
    for i in range(len(lams)):
        f=open(filename+'_lam'+str(lams[i])+'_fs'+str(fs)+'.dat')
        maxI=0
        for lines in f:
            if lines[0]=='#' or lines[0]=='@':
                continue
            if len(lines)==0:
                continue
            foo=lines.split()
            foo=[float(x) for x in foo ]
            foo_max=max(foo)
            if foo_max>maxI:
                maxI=foo_max
        Imax[i]=maxI
       # print('lam '+str(lams[i])+' first I0:'+str(round(1/Imax[i],2)))
    return Imax

def get_I0s(filename, lams, fs, iterations=10):
    """ Print I0s to try in different iterations

    filename: str
        Starting string of filename of image intensity files.
    lams: array of int
        Wavelength of all fluorophore types.
    fs: int
        Full-Width-at-Half-Maximum (FWHM) scaling factor.
    iterations: int, optional
        number of iterations of I0 prediction
    """
    Imax=get_maxI(filename,lams,fs)
    Imax=np.array(Imax)
    for i in range(iterations):
        I0=np.divide(1,Imax)
        string=''
        for j in range(len(I0)):
            string+=str(round(I0[j],2))+' '
        print(string)
        Imax=Imax-1
        for j in range(len(Imax)):
            if Imax[j]<0:
                Imax[j]=Imax[j]+1

def get_hist(filename, lams, fs, outname, maxI=20, dI=0.1,norm=False):
    """ Get histogram from I(l',m')
 
    Parameters
    ----------
    filename: str
        Starting string of filename of image intensity files.
    lams: array of int
        Wavelength of all fluorophore types.
    fs: int
        Full-Width-at-Half-Maximum (FWHM) scaling factor.
    outname: str
        Output filename
    maxI: float, optional
        Maximum to consider while generating the histogram. 
        (default value is 20)
    dI: float, optional
        Width of the bin in the histogram. (default value is 0.1)
    norm: Bool, optional
        If true, normalizes the histogram count such that maximum count
        is one.
  
    Writes
    -------
    outname: comma seprated data
        Histogram of all fluorophore types.
    """
    hist=np.zeros((int(maxI/dI)+1,len(lams)),dtype=int)
    for i in range(len(lams)):
        f=open(filename+'_lam'+str(lams[i])+'_fs'+str(fs)+'.dat')
        for lines in f:
            if lines[0]=='#' or lines[0]=='@':
                continue
            if len(lines)==0:
                continue
            foo=lines.split()
            foo=[float(x) for x in foo ]
            for j in range(len(foo)):
                I=int(foo[j]/dI)
                if I>1:
                    hist[I,i]+=1

    if norm==True: #Normalize count
        maxCnt=np.amax(hist,axis=0)
        hist=np.divide(hist,maxCnt)
     
    
    w=open(outname,'w')
    w.write('#iteration')
    for i in range(len(lams)):
        w.write(','+str(lams[i])+' nm')
    w.write('\n')
    for j in range(len(hist)):
        w.write(str(round((j+0.5)*dI,4)))
        for i in range(len(lams)):
            w.write(','+str(round(hist[j,i],4)))
        w.write('\n')
    w.close()

def get_num_area(filename, I0, lam, T, ti, fs, threshold, MaxBox, dlmn, pbc, \
    outname, show=False, white=False):
    """ Calculates number of particles and their area

    Parameters
    ----------
    filename: str
        Starting string of filename of image intensity files.
    I0: float
        The maximum image intensity
    lam: array of int
        Wavelength of the fluorophore type to threshold.
    T: int
        Number of timesteps to perform an average.
    ti: int
        timestep of the image data file. -1 if there is no sense of time
        (static system).
    fs: int
        Full-Width-at-Half-Maximum (FWHM) scaling factor.
    threshold: int
        An integer between 0-255. Image intensity above thresold/255 implies
        the presence of the particle. Periodic boundary conditions are applied.
    MaxBox: array of ints
        Contains the number of pixels in the image in l and m directions.
    dlmn: array of floats
        Delta l', Delta m', and Delta n'. Voxel size.
    pbc: array of ints
        1 if pbc condition is applied and 0 otherwise. First element for l', 
       and second element for m' 
    outname: str
        Output file name.
    show: Boolean, optional
        If true the Particles are shown. (default is False)
    white: Boolean, optional
        If true particles are shown as white or else shown as gray. (default is
        False)
    Writes
    ------
    [outname].dat:  
        Each line containing particle ID and area.
    [outname].png:
        The thresholded image (if show=True).
    """
    IMG_wf=get_grey_img(filename, I0, lam, T, ti, fs, MaxBox, whiteframe=True)
    IMG=[]
    for i in range(len(IMG_wf)):
        foo=IMG_wf[i][IMG_wf[i]>-1E-6]
        if len(foo)>0:
            foo=[int(x+1-threshold/255) for x in foo]
            IMG.append(foo)
    ly=len(IMG)
    lx=len(IMG[0])
    particles=[]
    #axis 0 is particle ID, axis 2 contains particle position i*lx+j
    for i in range(ly):
        for j in range(lx):
            if IMG[i][j]==0:
                continue
            idx=i*lx+j
            surr=[]
            if pbc[1]==1: #pbc in m'
                surr+=[i*lx+(j+1)%lx, i*lx+(j-1)%lx]
            else:
                surr+=[i*lx+j+1, i*lx+j-1]
            if pbc[0]==1: #pbc in l'
                surr+=[((i+1)%ly)*lx+j, ((i-1)%ly)*lx+j]
            else:
                surr+=[(i-1)*lx+j, (i+1)*lx+j]
            surr=set(surr)
            add_in=[]
            for n in range(len(particles)):
                common=particles[n].intersection(surr)
                if len(common)>0:
                    add_in.append(n)
            if len(add_in)>0:
                particles[add_in[0]].add(idx)
                for k in range(len(add_in)-1,0,-1):
                    particles[add_in[0]] |= particles[add_in[k]]
                    particles.pop(add_in[k])
            else:
                particles.append(set([idx]))

    tot_area=0
    w=open(outname+'.dat','w')
    for i in range(len(particles)):
        area=len(particles[i])*dlmn[0]*dlmn[1]
        tot_area+=area
        print('Particle: '+str(i)+' Area: '+str(area))
        w.write(str(i)+','+str(round(area,6))+'\n')
    print('Total Particle: '+str(len(particles))+' Avg area: '+ \
          str(tot_area/len(particles)))
    
    import matplotlib.pyplot as plt
    
    #Plot the particles as grey. If the algorithm fails, white color  
    #will be visible.
    for n1 in range(len(particles)):
        A=list(particles[n1])
        for n2 in range(len(A)):
            a=A[n2]
            i=int(a/lx)
            j=a-i*lx
            if white:
                IMG[i][j]=255
            else:
                IMG[i][j]=128

    #Add the white frame.
    dh=len(IMG_wf)-len(IMG)
    dw=len(IMG_wf[0])-len(IMG[0])

    for i in range(len(IMG_wf)):
        if IMG_wf[i][int(dw/2)+5]>0:
            break
    for j in range(len(IMG_wf)):
        if IMG_wf[int(dh/2)+5][j]>0:
            break
    IMG_new=np.ones((len(IMG_wf),len(IMG_wf[0])),dtype=int)
    IMG_new=IMG_new*255
    IMG_new[i:len(IMG)+i,j:len(IMG[0])+j]=IMG
     
    
    img_width=len(IMG_new[1])*3.0/len(IMG_new[0])
    fig,ax=plt.subplots(1,1,figsize=(img_width,3))
    plt.imshow(IMG_new,vmin=0,vmax=255,cmap='gray')
    ax.set_xticks([])
    ax.set_yticks([])
    plt.axis('off')
    plt.tight_layout(pad=0)
    if show:
        plt.show()
    else:
        plt.savefig(outname+'.png',dpi=1200)
     

