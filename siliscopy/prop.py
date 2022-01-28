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
import networkx as nx
import tifffile as tif
import copy
from .plot_image import get_grey_img
from numba import njit, prange
import multiprocessing as mp

def get_maxI(filename, lams, fs, tstep=None):
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
    if tstep!=None:
        filename=filename+str(tstep)
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
        f.close()
        Imax[i]=maxI
    return Imax

def get_I0s(filename, lams, fs, iterations=10, tstep=None):
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
    if tstep!=None:
        filename+=str(tstep)
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

def get_hist(filename, lams, fs, outname, maxI=20, dI=0.1,norm=False, \
    tstep=None):
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
    if tstep!=None:
        filename+=str(tstep)
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
        f.close()
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

def str2array(string,dtype,width=None):
    """ Convert string to array.
    
    Parameters
    ----------
    string: str
        string of text
    dtype: str
        dtype of values in the string (other than [,])
    width: int
        number of elements present in a row of 2D array. default (None)

    Return
    ------
    arr: 1D or 2D array.

    """
    string=string.replace('[','')    
    string=string.replace(']','')
    if ',' in string:    #bounds
        arr=np.fromstring(string,dtype=dtype,sep=',')
    else: #pbc
        arr=np.fromstring(string,dtype=dtype,sep=' ')

    if width!=None:
        n=len(arr)
        m=int(n/width+0.5)
        arr=arr.reshape(m,width)
    if width==None:
        arr=list(arr)
    else:
        arr=list(arr)
        for i in range(len(arr)):
            arr[i]=list(arr[i])
    return arr

def get_num_area(filename, threshold, outname='test', min_pix=1, 
    col_channel=0, ncoor=0, write=False):
    """ Calculates number of particles and their area

    Parameters
    ----------
    filename: str
        Filename of 2D or 2DT image intensity tiff files.
    threshold: int
        An integer between 0-1. Image intensity above thresold/255 implies
        the presence of the particle. Periodic boundary conditions are applied.
    pbc: array of ints
        1 if pbc condition is applied and 0 otherwise. First element for l', 
       and second element for m' 
    outname: str
        Output file name.
    min_pix: int, optional
        Any group of pixles containing less than min_pix pixels are not 
        considered as particles.

    Returns
    ------
    bin:
        Binary image of the input image, showing particles in white, background
        is black. 
    areas:
  
    """
    #####TIFF 
    IMG=tif.imread(filename) 
    IMG_dat=tif.TiffFile(filename)
    dtype=IMG.dtype
    axes=IMG_dat.series[0].axes
    xres=IMG_dat.pages[0].tags['XResolution'].value
    yres=IMG_dat.pages[0].tags['YResolution'].value
    darea=xres[1]*yres[1]/float(xres[0]*yres[0])
    bounds=IMG_dat.imagej_metadata['bounds'] #In str format.
    pbc=IMG_dat.imagej_metadata['pbc'] #In str format.
    pbc=str2array(pbc,'int')
    threshold=threshold*np.iinfo(IMG.dtype).max
    if col_channel==None:
        col_channel=0

    if 'Z' in axes:
        nres=IMG_dat.imagej_metadata['spacing']
        nidx=int(ncoor/nres+0.5)
        if ncoor<0: raise Exception('ncoor must be positive')
        bounds=str2array(bounds,'int',6)
    else:
        bounds=str2array(bounds,'int',4)
    if 'T' in axes:
        fpns=IMG_dat.imagej_metadata['finterval']
        funit=IMG_dat.imagej_metadata['funit']
    else:
        fpns=1
        funit='ns'


    if axes=='CYX':
        IMG=IMG[col_channel,:,:]
    elif axes=='ZYX':
        IMG=IMG[bounds[0][0]+nidx,:,:]
    elif axes=='ZCYX':
        IMG=IMG[bounds[0][0]+nidx,col_channel,:,:]
    elif axes=='TCYX':
        IMG=IMG[:,col_channel,:,:]
    elif axes=='TZYX':
        foo=copy.deepcopy(IMG)
        IMG=IMG[:,bounds[0][0]+nidx,:,:]
        sha0=foo.shape
        for t in range(1,sha0[0]):
            IMG[t,:,:]=foo[t,bounds[t][0]+nidx,:,:]
        
    elif axes=='TZCYX':
        foo=copy.deepcopy(IMG)
        IMG=IMG[:,bounds[0][0]+nidx,col_channel,:,:]
        sha0=foo.shape
        for t in range(1,sha0[0]):
            IMG[t,:,:]=foo[t,bounds[t][0]+nidx,col_channel,:,:]
  
    sha=IMG.shape #YX or TYX
    if len(sha)==2: #if YX
        IMG=IMG.reshape(1,sha[0],sha[1])
    Bin=np.zeros(IMG.shape,dtype='int')   
    tN,yN,xN=Bin.shape
    areas_t=[]
    for t in range(tN):
        max_label=1
        equiv={}
        for j in range(bounds[t][-2],bounds[t][-1]): #Over X (right)
            for i in range(bounds[t][-4],bounds[t][-3]): #Over Y (dowm)
                if IMG[t,j,i]<threshold:
                    continue
                labels=np.array([0,0,0,0]) #i, j, ipbc, jpbc
                if i>bounds[t][-4]: # if not on top edge, current label0 = 
                                    # label of pixel to the top.
                    labels[0]=Bin[t,j,i-1]
                if j>bounds[t][-2]: # if not on the left edge, current label1 = 
                                    # label of pixel to the ;eft.
                    labels[1]=Bin[t,j-1,i]
                #PBC
                if i==bounds[t][-3]-1 and pbc[0]==1: # if on the bottom edge, 
                            # label2 = label of the top pixel in the row.
                    labels[2]=Bin[t,j,bounds[t][-4]]
                if j==bounds[t][-1]-1 and pbc[1]==1: # if on the right edge, 
                            # label2 = label of the leftmost pixel in the row.
                    labels[3]=Bin[t,bounds[t][-2],i]
                
                foo=np.where(labels>0)
                if len(foo[0])==0:# No connecting pixel has a label. Add new.
                    Bin[t,j,i]=max_label
                    max_label+=1
                elif len(foo[0])==1:# Only one connecting pixel has a lable. 
                                    # Add the label to current
                    Bin[t,j,i]=labels[foo[0][0]]
                elif len(foo[0])>1:# Add the first label to current and 
                                   # add equivalence relations.
                    Bin[t,j,i]=labels[foo[0][0]]
                    for x1 in range(len(foo[0])-1):
                        for x2 in range(x1+1,len(foo[0])):
                            equiv=add_equiv(equiv, labels[foo[0][x1]], \
                                            labels[foo[0][x2]])

        for k in range(1,max_label):
            if k not in equiv:
                continue
            Bin[t,Bin[t]==k]=equiv[k] 
        areas=[]
        foo_idx=[]
        for k in range(1,max_label):
            if k in equiv:
                continue
            A=np.where(Bin[t]==k)
            if len(A[0])>=min_pix:
                areas.append(len(A[0]))
                foo_idx.append(k)
            else:
                Bin[t][A]=0
        areas=[float(x)*darea for x in areas]
        ele=sorted(np.unique(Bin[t]))
        ele.remove(0)
        N_ele=len(ele)
        if N_ele>0:
            foo=int(0.8*np.iinfo(dtype).max/N_ele+0.5)
        where=[]
        for i in range(len(ele)):
            where.append(np.where(Bin[t]==ele[i]))
        for i in range(len(ele)):
            foo1=(i+1)*foo+int(0.2*np.iinfo(dtype).max)
            Bin[t][where[i]]=min(foo1,np.iinfo(dtype).max)
        areas_t.append(areas)
    Bin=Bin.astype(dtype)
    if write:
        w=open(outname,'w')
        w.write('#t,Areas\n')
        for t in range(tN):
            w.write(str(round(t*fpns,5)))
            for i in range(len(areas_t[t])):
                w.write(','+str(round(areas_t[t][i],4)))
            w.write('\n')
        w.close()
        tif.imwrite('bin_'+filename,Bin,imagej=True, 
            resolution=(xres[0]/float(xres[1]), yres[0]/float(yres[0])), 
            metadata={'unit': 'nm', 'finterval':fpns , 'funit':funit,
                'axes': 'TYX', 'bounds':bounds, 'pbc':pbc})
    else:
        return Bin, areas_t 

def add_equiv(equiv,x1,x2):
    """ Add new equivalence relation.
    
    Parameters
    ----------
    equiv: dictionary
        Contains all equivalence relations a->b such that a > b.
    x1: int
    x2: int
        foo2->foo1 and foo2 > foo1.

    Returns
    -------
    equiv: New simplified equivalence relation dictionary.
    """
    if x1==x2: return equiv
    foo1=min(x1,x2)
    foo2=max(x1,x2)
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

def get_num_vol(filename, threshold, outname='test', min_pix=1, 
    col_channel=0, write=False):
    """ Calculates number of particles and their area

    Parameters
    ----------
    filename: str
        Filename of 2D or 2DT image intensity tiff files.
    threshold: int
        An integer between 0-1. Image intensity above thresold/255 implies
        the presence of the particle. Periodic boundary conditions are applied.
    pbc: array of ints
        1 if pbc condition is applied and 0 otherwise. First element for l', 
       and second element for m' 
    outname: str
        Output file name.
    min_pix: int, optional
        Any group of pixles containing less than min_pix pixels are not 
        considered as particles.

    Returns
    ------
    bin:
        Binary image of the input image, showing particles in white, background
        is black. 
    vol_t:
  
    """
    #####TIFF 
    IMG=tif.imread(filename) 
    dtype=IMG.dtype
    IMG_dat=tif.TiffFile(filename)
    axes=IMG_dat.series[0].axes
    xres=IMG_dat.pages[0].tags['XResolution'].value
    yres=IMG_dat.pages[0].tags['YResolution'].value
    bounds=IMG_dat.imagej_metadata['bounds']
    bounds=str2array(bounds,'int',6)
    pbc=IMG_dat.imagej_metadata['pbc']
    pbc=str2array(pbc,'int')
    threshold=threshold*np.iinfo(IMG.dtype).max
    nres=IMG_dat.imagej_metadata['spacing']
    dvol=xres[1]*yres[1]/float(xres[0]*yres[0])*nres
    if 'T' in axes:
        fpns=IMG_dat.imagej_metadata['finterval']
        funit=IMG_dat.imagej_metadata['funit']
    else:
        fpns=1
        funit='ns'

    if col_channel==None:
        col_channel=0

    if axes=='ZCYX':
        IMG=IMG[:,col_channel,:,:]
    elif axes=='TZCYX':
        IMG=IMG[:,:,col_channel,:,:]
    sha=IMG.shape
    if len(sha)==3: 
        IMG=IMG.reshape(1,sha[0],sha[1],sha[2])
    Bin=np.zeros(IMG.shape,dtype='int')    
    tN,zN,yN,xN=Bin.shape
    vols_t=[]
    for t in range(tN):
        max_label=1
        equiv={}
        for k in range(bounds[t][0],bounds[t][1]): #Over Z 
            for j in range(bounds[t][4],bounds[t][5]): #Over X (right)
                for i in range(bounds[t][2],bounds[t][3]): #Over Y (down)
                    if IMG[t,k,j,i]<threshold:
                        continue
                    labels=np.array([0,0,0,0,0,0]) #i, j, k, ipbc, jpbc, kpbc
                    if i>bounds[t][2]:
                        labels[0]=Bin[t,k,j,i-1]
                    if j>bounds[t][4]:
                        labels[1]=Bin[t,k,j-1,i]
                    if k>bounds[t][0]:
                        labels[2]=Bin[t,k-1,j,i]
                    #PBC
                    if i==bounds[t][3]-1 and pbc[0]==1:
                        labels[3]=Bin[t,k,j,bounds[t][2]]
                    if j==bounds[t][5]-1 and pbc[1]==1:
                        labels[4]=Bin[t,k,bounds[t][4],i]
                    if k==bounds[t][1]-1 and pbc[2]==1:
                        labels[5]=Bin[t,bounds[t][0],j,i]
                    
                    foo=np.where(labels>0)
                    if len(foo[0])==0:# Add a new label
                        Bin[t,k,j,i]=max_label
                        max_label+=1
                    elif len(foo[0])==1:# Add the label to current
                        Bin[t,k,j,i]=labels[foo[0][0]]
                    elif len(foo[0])>1:# Add the first label to current 
                                       # and add equiv
                        Bin[t,k,j,i]=labels[foo[0][0]]
                        for x1 in range(len(foo[0])-1):
                            for x2 in range(x1+1,len(foo[0])):
                                equiv=add_equiv(equiv,labels[foo[0][x1]], \
                                    labels[foo[0][x2]])
            print('t = '+str(round((t+1)/tN*100,2))+ '%  '+ 
               str(round((k+1-bounds[t][0])/(bounds[t][1]-bounds[t][0])*100,2))+
                "% done    ",end='\r')
        for l in range(1,max_label):
            if l not in equiv:
                continue
            Bin[t,Bin[t]==l]=equiv[l] #Check how to do
        vols=[]
        for l in range(1,max_label):
            if l in equiv:
                continue
            A=np.where(Bin[t]==l)
            if len(A[0])>=min_pix:
                vols.append(len(A[0]))
            else:
                Bin[t][A]=0
        vols=[float(x)*dvol for x in vols]
        ele=sorted(np.unique(Bin[t]))
        ele.remove(0)
        foo=int(0.8*np.iinfo(dtype).max/len(ele)+0.5)
        where=[]
        for i in range(len(ele)):
            where.append(np.where(Bin[t]==ele[i]))
        for i in range(len(ele)):
            foo1=(i+1)*foo+int(0.2*np.iinfo(dtype).max)
            Bin[t][where[i]]=min(foo1,np.iinfo(dtype).max)
        vols_t.append(vols)

    Bin=Bin.astype('uint8')
    if write:
        w=open(outname,'w')
        w.write('#t,vols\n')
        for t in range(tN):
            w.write(str(t*fpns))
            for i in range(len(vols_t[t])):
                w.write(','+str(round(vols_t[t][i],4)))
            w.write('\n')
        w.close()
        tif.imwrite('bin_vol_'+filename,Bin,imagej=True,
            resolution=(xres[0]/float(xres[1]), yres[0]/float(yres[1])), 
            metadata={'unit': 'nm', 'finterval': fpns, 'funit': funit, 
            'axes': 'TZYX', 'bounds':bounds, 'pbc':pbc,'spacing':nres})
    else:
        return Bin, vols_t 

def max_Nf_t(filename):
    """ Get max number of radiative emissions for a timestep.
    
    Parameters
    ----------
    filename: str
        filename of a specimen format .spm 

    Returns
    -------
    max_Nf: maximum radiative emissions 

    """
    f=open(filename,'r')
    print('filename = '+filename)
    cnt=0
    max_Nf=0
    for lines in f:
        if cnt==0:
            foo=lines.split()
            Natoms=int(foo[1])
        elif cnt>Natoms:
            break
        elif cnt>0 :
            num_Nf=int(lines[70:75])
            for i in range(num_Nf):
                Nf=int(lines[75+15*i:85+15*i])
                if Nf>max_Nf:
                    max_Nf=Nf
        cnt+=1
    f.close()
    return max_Nf


def read_ppfile(filename):

    # pp_file format
    # Nfluor_type dt
    # [ lam1 ]
    # Nstate Nfluor_state
    # states
    # si sj kij
    # ...
    # si sj -1 lam1_1
    # ...
    f=open

def max_Nf(filename, tbegin, tmax, tdiff, pp_file, mprocess=False):
    """ Get max number of radiative emissions
    
    Parameters
    ----------
    filename: str
        filename of a specimen format .spm 
    tbegin: int
        index of first timestep, which is included.
    tmax: int
        index of maximum timestep, tmax is not included.
    tdiff: int 
        difference between timesteps to generate the 3DT image.
    mprocess: bool
        If true multiprocessing is used. (default False)

    Returns
    -------
    max_Nf: maximum radiative emissions 

    """
        

    Arguments=[]
    for i in range(tbegin,tmax,tdiff):
        Arguments.append(filename+str(i)+'.spm')
    if mprocess==True:
        cpus=mp.cpu_count()
        if len(Arguments)<cpus:
            cpus=len(Arguments)
        pool=mp.Pool(cpus)
        results=pool.map(max_Nf_t,Arguments)
        pool.close()
    else:
        results=[]
        for i in range(len(Arguments)):
            results.append(max_Nf_t(Arguments[i]))

    max_Nf=max(results)
    print('max_Nf = '+str(max_Nf)) 



def get_fcs(IMGname,outname,fcs_tmax,nidx=None):
    """ Get Fluorescence Correlation Spectroscopy data from an tiff image.

    Parameters
    ----------
    IMGname: str
        Filename of tiff image.
    fcs_tmax: float?
        Maximum time for the correlation function
    nidx: int (Optional)
        Index of the z coordinate in the tiff file. (default None)

    Writes
    ------
    [outname].pickle: Writes a pickle file storing the fcs.
 
    """
    IMG_info=tif.TiffFile(IMGname)
    xres=IMG_info.pages[0].tags['XResolution'].value
    yres=IMG_info.pages[0].tags['YResolution'].value
    axes=IMG_info.series[0].axes
    dtype=IMG_info.series[0].dtype
    if 'T' not in axes:
        raise Exception('Time should be in Axes of the IMG')
    if 'Z' in axes:
        bounds=str2array(IMG_info.imagej_metadata['bounds'],dtype,width=6)
        if nidx==None:
            raise Exception('Not applicable for 3D images')
    else:
        bounds=str2array(IMG_info.imagej_metadata['bounds'],dtype,width=4)

    dt=1.0/float(IMG_info.imagej_metadata['finterval'])
    IMG=tif.imread(IMGname)
    sha=IMG.shape
    Nmax=int(fcs_tmax/dt)
    if Nmax>sha[0]:
        Nmax=sha[0]
    #find strict bounds
    xbnds=[0,sha[-1]]
    ybnds=[0,sha[-2]]
    for t in range(sha[0]):
        if bounds[t][-2]>xbnds[0]:
            xbnds[0]=bounds[t][-2]
        if bounds[t][-1]<xbnds[1]:
            xbnds[1]=bounds[t][-1]
        if bounds[t][-4]>ybnds[0]:
            ybnds[0]=bounds[t][-4]
        if bounds[t][-3]<ybnds[1]:
            ybnds[1]=bounds[t][-3]
    #TZCYX/ TZYX / TCYX / TYX
    if axes == 'TZCYX': #TZCYX -> TCYX
        newIMG=np.zeros((sha[0],sha[2],ybnds[1]-ybnds[0],xbnds[1]-xbnds[0]),dtype=dtype)
        for t in range(sha[0]):
            newIMG[t,:,:,:]=IMG[t,nidx+bounds[t][0],:,ybnds[0]:ybnds[1],xbnds[0]:xbnds[1]]
    elif axes == 'TZYX': #TZYX -> TYX
        newIMG=np.zeros((sha[0],sha[2],sha[3]),dtype=dtype)
        for t in range(sha[0]):
            newIMG[t,:,:,]=IMG[t,nidx+bounds[t][0],ybnds[0]:ybnds[1],xbnds[0]:xbnds[1]]
    elif axes == 'TYX': 
        newIMG=IMG[:,ybnds[0]:ybnds[1],xbnds[0]:xbnds[1]]
    elif axes == 'TCYX': 
        newIMG=IMG[:,:,ybnds[0]:ybnds[1],xbnds[0]:xbnds[1]]

    newIMG=newIMG.astype('float')
    #now image is either TCYX or TYX
    if 'C' in axes: #TCYX
        foo=axes.find('C')
        t=np.linspace(0,Nmax-1,Nmax)*dt
        data=np.zeros((Nmax,sha[foo]+1)) #time, acf1, err, acf2, err,...
        data[:,0]=t[:]
        for i in range(sha[foo]):#color channels
            ACF=get_acf(newIMG[:,i,:,:],Nmax)
            data[:,1+i]=ACF[:]
    else:
        ACF=get_acf(newIMG,Nmax)
        t=np.linspace(0,Nmax-1,Nmax)*dt
        data=np.zeros((Nmax,3))
        data[:,0]=t[:]
        data[:,1]=ACF[:]

    if len(outname)>7:
        if outname[-7:]=='.pickle':
            nx.write_gpickle(data,outname)
        else:
            nx.write_gpickle(data,outname+'.pickle')
    else:
        nx.write_gpickle(data,outname+'.pickle')

@njit(parallel=True)
def get_acf(I,Nmax):
    """ Get auto correlation of fluorescence intensity

    Parameters
    ----------
    I: numpy array of unsigned integers.
        Tiff image -> fluorescence intensity. Axes is TYX.
    Nmax: int
        Maximum index to find autocorrelation

    Returns
    -------
    ACF: Autocorrelation function of fluorescnce intensity 
    """
    sha=I.shape
    Iavg=np.sum(I,axis=0)
    Iavg=Iavg/sha[0]
    Iavg2=np.square(Iavg)
    ACF=np.zeros((Nmax,sha[1],sha[2]))
    for i in prange(Nmax):
        for j in prange(sha[0]-i):
            ACF[i,:,:]+=np.multiply(I[j,:,:],I[j+i,:,:])
        ACF[i,:,:]=np.divide(ACF[i,:,:],Iavg2*(sha[0]-i))
    #average
    ACF=ACF-1
    ACF_avg=np.sum(ACF,axis=1)
    ACF_avg=np.sum(ACF_avg,axis=1)
    ACF_avg=ACF_avg/(sha[1]*sha[2])
    return ACF_avg

def get_fccs(IMGname, outname, tmax, c1, c2, nidx=None):
    """ Get Fluorescence Cross-Correlation Spectroscopy data from an tiff image.

    Parameters
    ----------
    tmax: float
        Maximum time for the correlation function
    c1: int
        Color channel 1
    c2: int
        Color channel 2
    nidx: int (Optional)
        Index of the z coordinate in the tiff file. (default None)

    Writes
    ------
    [outname].pickle: Writes a pickle file storing the fccs.
 
    """
    IMG_info=tif.TiffFile(IMGname)
    xres=IMG_info.pages[0].tags['XResolution'].value
    yres=IMG_info.pages[0].tags['YResolution'].value
    axes=IMG_info.series[0].axes
    dtype=IMG_info.series[0].dtype
    if 'T' not in axes:
        raise Exception('Time should be in Axes of the IMG')
    if 'C' not in axes:
        raise Exception('Color should be in Axes of the IMG')
    if 'Z' in axes:
        bounds=str2array(IMG_info.imagej_metadata['bounds'],dtype,width=6)
        if nidx==None:
            raise Exception('Not applicable for 3D images')
    else:
        bounds=str2array(IMG_info.imagej_metadata['bounds'],dtype,width=4)

    dt=1.0/float(IMG_info.imagej_metadata['finterval'])
    IMG=tif.imread(IMGname)
    sha=IMG.shape
    Nmax=int(tmax/dt)
    if Nmax>sha[0]:
        Nmax=sha[0]

    #find strict bounds
    xbnds=[0,sha[-1]]
    ybnds=[0,sha[-2]]
    for t in range(sha[0]):
        if bounds[t][-2]>xbnds[0]:
            xbnds[0]=bounds[t][-2]
        if bounds[t][-1]<xbnds[1]:
            xbnds[1]=bounds[t][-1]
        if bounds[t][-4]>ybnds[0]:
            ybnds[0]=bounds[t][-4]
        if bounds[t][-3]<ybnds[1]:
            ybnds[1]=bounds[t][-3]
    #TZCYX/ TCYX 
    if axes == 'TZCYX': #TZCYX -> TCYX
        newIMG=np.zeros((sha[0],sha[2],sha[3],sha[4]),dtype=dtype)
        for t in range(sha[0]):
            newIMG[t,:,:,:]=IMG[t,nidx+bounds[t][0],:,ybnds[0]:ybnds[1],xbnds[0]:xbnds[1]]
    else:
        newIMG=IMG[:,:,ybnds[0]:ybnds[1],xbnds[0]:xbnds[1]]
    newIMG=newIMG.astype(float)
    t=np.linspace(0,Nmax-1,Nmax)*dt
    data=np.zeros((Nmax,2))
    data[:,0]=t[:]
    CCF=get_ccf(newIMG,Nmax,c1,c2)
    data[:,1]=CCF[:]
    if len(outname)>7:
        if outname[-7:]=='.pickle':
            nx.write_gpickle(data,outname)
        else:
            nx.write_gpickle(data,outname+'.pickle')
    else:
        nx.write_gpickle(data,outname+'.pickle')


@njit(parallel=True)
def get_ccf(I,Nmax,c1,c2):
    """ Get cross-correlation of fluorescence intensity

    Parameters
    ----------
    I: numpy array of unsigned integers.
        Tiff image -> fluorescence intensity. Aces is TCXY.
    Nmax: int
        Maximum index to find autocorrelation

    Returns
    -------
    CCF: Cross-correlation function of fluorescnce intensity 
    """
    sha=I.shape
    Iavg1=np.sum(I[:,c1,:,:],axis=0)
    Iavg1=Iavg1/sha[0]
    Iavg2=np.sum(I[:,c2,:,:],axis=0)
    Iavg2=Iavg2/sha[0]
    Iavg12=np.multiply(Iavg1,Iavg2)
    CCF=np.zeros((Nmax,sha[2],sha[3]))
    for i in prange(Nmax):
        for j in prange(sha[0]-i):
            CCF[i,:,:]+=np.multiply(I[j,c1,:,:],I[j+i,c2,:,:])
        CCF[i,:,:]=np.divide(CCF[i,:,:],Iavg12*(sha[0]-i))
    CCF=CCF-1
    CCF_avg=np.sum(CCF,axis=1)
    CCF_avg=np.sum(CCF_avg,axis=1)
    CCF_avg=CCF_avg/(sha[1]*sha[2])
    return CCF_avg

