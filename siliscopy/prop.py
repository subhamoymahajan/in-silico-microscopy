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
import tifffile as tif
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
    lams: 
        See ''get_maxI''
    fs: 
        See ''get_maxI''
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
    lams: 
        See ''get_maxI''
    fs: 
        See ''get_maxI''
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

def get_num_area(filename, threshold, dlmn, pbc, outname, min_pix=0):
    """ Calculates number of particles and their area

    Parameters
    ----------
    filename: str
        Filename of 2D or 2DT image intensity tiff files.
    threshold: int
        An integer between 0-255. Image intensity above thresold/255 implies
        the presence of the particle. Periodic boundary conditions are applied.
    pbc: array of ints
        1 if pbc condition is applied and 0 otherwise. First element for l', 
       and second element for m' 
    outname: str
        Output file name.
    min_pix: int, optional
        Any group of pixles containing less than min_pix pixels are not 
        considered as particles.

    Writes
    ------
    bin_[filename]:
        Binary image of the input image, showing particles in white, background
        is black. 
    [outname]:
        Each line contains, timestep (based on tiff), number of particles and 
        total cross-sectional area.
  
    """
    #####TIFF 
    IMG=tif.imread(filename) #Should be a TYX tiff file i.e. 2DT monochrome
    fps=1 
    consts=IMG.shape
    if len(consts)==2:
        foo=np.zeros((1,consts[0],consts[1]),dtype=IMG.dtype)
        foo[0,:,:]=IMG[:,:]
        IMG=foo
        consts=IMG.shape
    ###Binary image
    if IMG.dtype==np.uint8:
        M=255
    elif IMG.dtype==np.uint16:
        M=65535

    IMG[IMG<threshold]=0
    IMG[IMG>threshold]=M
    Bin=np.zeros(consts,dtype=IMG.dtype) 
    w=open(outname,'w')
    lx,ly=[consts[2],consts[1]]
    for t in range(consts[0]):
        max_label=1
        boundary=[]
        equiv={}
        IMG2=np.zeros(consts[1:],dtype=int)
        for j in range(consts[1]):
            for i in range(consts[2]):
                if IMG[t,j,i]==0:
                    continue
                labels=[0,0] #i, j
                if i>0:
                    labels[0]=IMG2[j,i-1]
                if j>0:
                    labels[1]=IMG2[j-1,i]
                #PBC
                if i==consts[2]-1 and pbc[0]==1:
                    labels[1]=IMG2[0,i]
                if j==consts[1]-1 and pbc[1]==1:
                    labels[0]=IMG2[j,0]
                
                if labels[0]==0 and labels[1]!=0:
                    IMG2[j,i]=labels[1]
                elif labels[1]==0 and labels[0]!=0:
                    IMG2[j,i]=labels[0]
                elif labels[0]==labels[1] and labels[0]!=0:
                    IMG2[j,i]=labels[0]
                elif labels[0]*labels[1]>0 and labels[0]!=labels[1]:
                    IMG2[j,i]=labels[0]
                    #label1 = label0
                    foo1=min(labels[0],labels[1])
                    foo2=max(labels[0],labels[1])
                    equiv=add_equiv(equiv,foo1,foo2)
                else:
                    IMG2[j,i]=max_label
                    max_label+=1
                if i==0 or i==consts[2]-1:
                    if pbc[0]==0:
                        if IMG2[j,i] not in boundary:
                            boundary.append(IMG2[j,i])
                if j==0 or j==consts[1]-1:
                    if pbc[1]==0:
                        if IMG2[j,i] not in boundary:
                            boundary.append(IMG2[j,i])
        for k in range(max_label):
            if k not in equiv:
                continue
            IMG2[IMG2==k]=equiv[k]
        bound_new=[]
        for k in range(len(boundary)):
            if boundary[k] in equiv:
                if equiv[boundary[k]] not in bound_new:
                    bound_new.append(equiv[boundary[k]])
        tot_area=0
        n=0 
        for k in range(1,max_label):
            if k in equiv:
                continue
            if k in bound_new:
                IMG2[IMG2==k]=0
                continue
            A=np.where(IMG2==k)
            if len(A[0])>min_pix:
                tot_area+=len(A[0])
                n+=1
            else:
                IMG2[A]=0
        IMG2[IMG2>0]=255
        IMG2=IMG2.astype(IMG.dtype)
        Bin[t,:,:]=IMG2[:,:] 
        w.write(str(t)+','+str(n)+','+str(float(tot_area)*dlmn[0]*dlmn[1])+'\n')
        print('t = '+str(t)+'      ',end='\r')  
    tif.imwrite('bin_'+filename,Bin,imagej=True, resolution=(1./dlmn[0], 
            1./dlmn[1]), metadata={'unit': 'nm', 'finterval': fps, 
            'axes': 'TYX'}) 

def add_equiv(equiv,foo1,foo2):
    """ Add new equivalence relation.
    
    Parameters
    ----------
    equiv: dictionary
        Contains all equivalence relations a->b such that a > b.
    foo1: int
    foo2: int
        foo2->foo1 and foo2 > foo1.

    Returns
    -------
    equiv: New simplified equivalence relation dictionary.
    """
    if foo1 in equiv and foo2 in equiv: #If foo2->a2, foo1->a1 exists and 
        # foo2->foo1 is to be added. Let a1 be the smallest of a1,a2,foo1,foo2.
        # foo2->a1, foo1->a1, a2->a1.
        a=min(foo1,foo2,equiv[foo1],equiv[foo2])
        for b in [foo1,foo2,equiv[foo1],equiv[foo2]]:
            if b!=a:
                equiv[b]=a
        
    elif foo2 in equiv: #foo2->a exists and foo2->foo1 is to be added.
        if foo1<equiv[foo2]: # if a is larger than foo1. 
        #foo2->foo1 and a->foo1
            equiv[equiv[foo2]]=foo1
            equiv[foo2]=foo1
        elif foo1>equiv[foo2]: # if a is smaller than foo1. 
        #foo2->a and foo1->a.
            equiv[foo1]=equiv[foo2]

    elif foo1 in equiv: # foo1->a exists and foo2->foo1. So foo2->a.
        equiv[foo2]=equiv[foo1]
    else: #neither foo1, foo2 is present so simply add foo1->foo2
        equiv[foo2]=foo1
    equiv=simp_equiv(equiv)
    return equiv

def simp_equiv(equiv):
    """ Simplifies an equivalnce dictionary

    Maxes sure a->b is such that b is the lowest possible number.
    """
    for key in equiv:
        foo=equiv[key]
        while True: #iteratively, if a->b and b->c, then a->c
            if foo in equiv:
                foo=equiv[foo]
            else:
                break
        equiv[key]=foo
    return equiv

def get_num_vol(filename, threshold, dlmn, pbc, outname, min_pix=0):
    """ Calculates number of particles and their area

    Parameters
    ----------
    filename: str
        Filename of 3D or 3DT monochrome image intensity tiff file. Or it should
        be a .dat file containing file names of different 2D or 2DT monochrome
        image intensity tiff files.
    threshold: int
        An integer between 0-255. Image intensity above thresold/255 implies
        the presence of the particle. Periodic boundary conditions are applied.
    pbc: array of ints
        1 if pbc condition is applied and 0 otherwise. First element for l', 
        , second element for m' and third for n'.
    outname: str
        Output file name.
    min_pix: int, optional
        Any group of pixles containing less than min_pix pixels are not 
        considered as particles.

    Writes
    ------
    bin_[filename]:
        Binary image of the input image, showing particles in white, background
        is black. 
    [outname]:
        Each line contains, timestep (based on tiff), number of particles and 
        total cross-sectional area.
  
    """
    #####TIFF 
    if filename[-4:]=='.dat':
        f=open(filename,'r')
        IMG=[]
        for lines in f:
            IMG.append(tif.imread(lines.strip()))
        IMG=np.array(IMG)
        consts=IMG.shape
        if len(consts)==3: #it is in ZYX ; many 2D
            IMG=np.array([IMG])
        if len(consts)==4: #ZTYX  ; many 2DT
            IMG=np.transpose(IMG,(1,0,2,3))#TZYX 
    else:
        IMG=tif.imread(filename)
        consts=IMG.shape
        if len(consts)==3: #ZYX ; 3D
            IMG=np.array([IMG])
        #if len(consts)==4: #TZYK ; 3DT

    fps=1 
    consts=IMG.shape

    IMG[IMG<threshold]=0
    IMG[IMG>threshold]=np.iinfo(IMG.dtype).max
    Bin=np.zeros(consts,dtype=IMG.dtype) 
    w=open(outname,'w')
    lx,ly,lz=consts[3],consts[2],consts[1]
    for t in range(consts[0]):
        max_label=1
        boundary=[]
        equiv={}
        IMG2=np.zeros(consts[1:],dtype=int)
        for k in range(consts[1]):
            for j in range(consts[2]):
                for i in range(consts[3]):
                    if IMG[t,k,j,i]==0:
                        continue
                    labels=[0,0,0]
                    if i>0:
                        labels[0]=IMG2[k,j,i-1]
                    if j>0:
                        labels[1]=IMG2[k,j-1,i]
                    if k>0:
                        labels[2]=IMG2[k-1,j,i] 

                    #PBC
                    if i==consts[3]-1 and pbc[0]==1:
                        labels[0]=IMG2[k,0,i]
                    if j==consts[2]-1 and pbc[1]==1:
                        labels[1]=IMG2[k,j,0]
                    if k==consts[1]-1 and pbc[2]==1:
                        labels[2]=IMG2[0,j,i]


                    idx=np.where(labels==0)
                    if len(idx[0])>0:#connected to at least one of voxel
                        IMG2[k,j,i]=labels[idx[0][0]]
                    else:
                        IMG2[k,j,i]=max_label
                        max_label+=1

                    if len(idx[0])>1: #connected to more than one voxels
                        #Need to add equivalence relations
                        #Either two are equal or all are equal.
                        foo=[[k,j,i-1],[k,j-1,i],[k-1,j,i]]                       
                        eq_lab=sorted([IMG2[foo[i][0],foo[i][1],foo[i][2]] \
                                      for i in idx[0]])
                        #Largest -> Smallest
                        equiv=add_equiv(equiv,eq_lab[0],eq_lab[-1])
                        #Mid -> Smalles or same as prev
                        equiv=add_equiv(equiv,eq_lab[0],eq_lan[1])

                    if (i==0 or i==consts[3]-1) and pbc[0]==0:
                        if IMG2[k,j,i] not in boundary:
                            boundary.append(IMG2[k,j,i])
                    elif (j==0 or j==consts[2]-1) and pbc[1]==0:
                        if IMG2[k,j,i] not in boundary:
                            boundary.append(IMG2[k,j,i])
                    elif (k==0 or k==consts[1]-1) and pbc[2]==0:
                        if IMG2[k,j,i] not in boundary:
                            boundary.append(IMG2[k,j,i])

        for l in range(max_label):
            if l not in equiv:
                continue
            IMG2[IMG2==l]=equiv[l]
        bound_new=[]
        for l in range(len(boundary)):
            if boundary[l] in equiv:
                if equiv[boundary[l]] not in bound_new:
                    bound_new.append(equiv[boundary[l]])
        tot_vol=0
        n=0 
        for l in range(1,max_label):
            if l in equiv: 
                continue
            if l in bound_new:
                IMG2[IMG2==l]=0
                continue
            A=np.where(IMG2==l)
            if len(A[0])>min_pix:
                tot_vol+=len(A[0])
                n+=1
            else:
                IMG2[A]=0

        IMG2[IMG2>0]=255
        IMG2=IMG2.astype(IMG.dtype)
        Bin[t,:,:]=IMG2[:,:] 
        w.write(str(t)+','+str(n)+','+str(float(tot_vol)*dlmn[0]*dlmn[1]*dn)+'\n')
        print('t = '+str(t)+'      ',end='\r')  
    tif.imwrite('bin_'+filename,Bin,imagej=True, resolution=(1./dlmn[0], 
            1./dlmn[1]), metadata={'spacing':dn, 'unit': 'nm', 'finterval': fps, 
            'axes': 'TZYX'}) 

