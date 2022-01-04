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
import cv2
import tifffile as tif
import multiprocessing as mp
from . import plot_image
from . import prop

def psf_dat2tiff(filename,outname,Plmn,dlmn,dtype='uint8',psf_type=0):
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
    psf_type: int
        0 implies depth invariant with circular symmetry, and 1 implies depth 
        variant with circular symmetry.
    
    Writes
    ------
    [outname]: A Tiff file containing the PSF.
    """
    f=open(filename,'r')
    N=[int(Plmn[i]/dlmn[i]+0.5)*2+1 for i in range(3)] #XYZ
    PSF=np.zeros((N[2],N[1],N[0]),dtype=dtype) #ZYX
    x0,y0,z0=[int((a-1)/2) for a in N]
    for lines in f:
        if len(lines)==0:
            continue
        if lines[0]=='#':
            continue
        foo=lines.split()
        foo=[float(x) for x in foo]
        xyz=np.divide(foo[:3],dlmn)
        xyz_int=np.zeros(3,dtype='int')

        for i in range(3):
            if foo[i]>=0:
                xyz_int[i]=int(foo[i]/dlmn[i]+0.5)
            else:
                xyz_int[i]=int(foo[i]/dlmn[i]-0.5)

        for i in range(3):
            if abs(xyz[i]-xyz_int[i])>1E-6:
                continue #not part of the dlmn
        x,y,z=xyz_int
        I=int(float(foo[3])*np.iinfo(dtype).max) 
        for i in [-1,1]:
            for j in [-1,1]:
                if psf_type==0:
                    for k in [-1,1]:
                        PSF[z0+k*z,x0+i*x,y0+j*y]=I
                        PSF[z0+k*z,y0+j*y,x0+i*x]=I
                elif psf_type==1:
                    PSF[z0+z,x0+i*x,y0+j*y]=I
                    PSF[z0+z,y0+i*y,x0+j*x]=I
    f.close()
    tif.imsave(outname,PSF,resolution=(1./dlmn[0],1./dlmn[1]),imagej=True,
               metadata={'axes':'ZYX', 'unit':'nm', 'spacing':dlmn[2] })

def psf_dat2tiff2(filename,outname,Plmn,dlmn,dtype='uint16',psf_type=0):
    """ Converts PSF.dat to PSF.tiff. 0-0.2 is black to red, 0.2-0.8 is red hue 
        to blue hue (0 to 240 degrees), and 0.8 to 1 is blue to white. All other 
        properties are same as ''psf_dat2tiff''
       
    """
    f=open(filename,'r')
    N=[int(Plmn[i]/dlmn[i]+0.5)*2-1 for i in range(3)] #XYZ
    PSF=np.zeros((N[2],3,N[1],N[0]),dtype=dtype) #ZYC
    x0,y0,z0=[int((a-1)/2) for a in N]
    
    for lines in f:
        if len(lines)==0:
            continue
        if lines[0]=='#':
            continue
        foo=lines.split()
        foo=[float(x) for x in foo]
        xyz=np.divide(foo[:3],dlmn)
        xyz_int=np.zeros(3,dtype='int')

        for i in range(3):
            if foo[i]>=0:
                xyz_int[i]=int(foo[i]/dlmn[i]+0.5)
            else:
                xyz_int[i]=int(foo[i]/dlmn[i]-0.5)

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

        col=(np.array(rgb)*np.iinfo(dtype).max).astype(np.iinfo(dtype)) 

        for i in [-1,1]:
            for j in [-1,1]:
                if psf_type==0:
                    for k in [-1,1]:
                        PSF[z0+k*z,:,x0+i*x,y0+j*y]=col[:]
                        PSF[z0+k*z,:,y0+j*y,x0+i*x]=col[:]
                elif psf_type==1:
                    PSF[z0+z,:,x0+i*x,y0+j*y]=col[:]
                    PSF[z0+z,:,y0+j*y,x0+i*x]=col[:]
    f.close()
    tif.imsave(outname,PSF,resolution=(1./dlmn[0],1./dlmn[1]),imagej=True,
               metadata={'spacing':dlmn[2], 'axes':'ZCYX', 'unit':'nm' })

def nstack2tiff(datafile, outname, spacing):
    """ Generate a tiff file from multiple n-slice tiffs. 

    Parameter
    ---------
    datafile: str
        File containing tiff file names in order of n coordinate
    outname: str
        Outname header for the tiff file. 
    Writes
    ------
    [outname]: A volume tiff image
    """ 

    f=open(datafile,'r')
    data=[]
    for lines in f:
        data.append(lines.strip())
    f.close()
    nN=len(data)
    #Read Meta data.
    img_data=tif.TiffFile(data[0])
    res=(img_data.pages[0].tags['XResolution'].value,
         img_data.pages[0].tags['YResolution'].value)
    metaD=img_data.imagej_metadata
    axes=img_data.series[0].axes
    dtyp=img_data.series[0].dtype
    sha=img_data.asarray().shape
    metaD2={}
    bounds=prop.str2array(metaD['bounds'],int,width=4)

    if 'T' in axes:
        metaD2['axes']='TZ'+axes[1:]
        consts=tuple([sha[0],nN]+list(sha[1:]))
        ### Create image   
        IMG=np.zeros(consts,dtype=dtyp)
        bounds_n=[]
        for n in range(nN):
            imgdata_n=tif.TiffFile(data[n])
            bounds_n.append(prop.str2array(imgdata_n.imagej_metadata['bounds'],int,width=4))
        for t in range(sha[0]):
            z0,zN=0,nN
            for i in range(nN):
                if bounds_n[i][t][0]>bounds_n[i][t][1] and bounds_n[i][t][2]>bounds_n[i][t][3]:
                    z0+=1
                else:
                    break
            for i in range(nN):
                if bounds_n[nN-i-1][t][0]>bounds_n[nN-1-i][t][1] and bounds_n[nN-1-i][t][2]>bounds_n[nN-1-i][t][3]:
                    zN-=1
                else:
                    break
            bounds[t]=[z0,zN]+bounds[t]

        if len(sha)==4: #TCYX
            for n in range(nN):
                IMG[:,n,:,:,:]=tif.imread(data[n])
        elif len(sha)==3: #TYX
            for n in range(nN):
                IMG[:,n,:,:]=tif.imread(data[n])
                
    else: #No time axis
        metaD2['axes']='Z'+axes
        consts=tuple([nN]+list(sha))
        ### Create image   
        IMG=np.zeros(consts,dtype=dtyp)
        z0,zN=0,nN
        bounds_n=[]
        for n in range(nN):
            imgdata_n=tif.TiffFile(data[n])
            bounds_n.append(prop.str2array(imgdata_n.imagej_metadata['bounds'],int,width=4))
        for i in range(nN):
            if bounds_n[i][0][0]>bounds_n[i][0][1] and bounds_n[i][0][2]>bounds_n[i][0][3]:
                z0+=1
            else:
                break
        for i in range(nN):
            if bounds_n[nN-i-1][0][0]>bounds_n[nN-1-i][0][1] and bounds_n[nN-1-i][0][2]>bounds_n[nN-1-i][0][3]:
                zN-=1
            else:
                break
        bounds[0]=[z0,zN]+bounds[0]
        for n in range(nN):
            IMG[n]=tif.imread(data[n])
            

    for a in ['finterval', 'unit', 'funit', 'pbc']:
        if a in metaD:
            metaD2[a]=metaD[a]
    metaD2['spacing']=spacing
    metaD2['bounds']=bounds
    tif.imwrite(outname, IMG, imagej=True, metadata=metaD2, 
            resolution=(float(res[0][0]/res[0][1]),float(res[1][0]/res[1][1])))
        
def tstack2tiff(datafile, fpns, outname):
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
        data.append(lines.strip())
    f.close()
    tN=len(data)
    img_data=tif.TiffFile(data[0])
    res=(img_data.pages[0].tags['XResolution'].value,
         img_data.pages[0].tags['YResolution'].value)
    metaD=img_data.imagej_metadata
    axes=img_data.series[0].axes
    bounds='['  
    IMG=np.zeros(tuple([tN]+list(img_data.series[0].shape)),
                 dtype=img_data.series[0].dtype)
    for t in range(tN):
        IMG[t]=tif.imread(data[t])
        foo=tif.TiffFile(data[t])
        bounds+=foo.imagej_metadata['bounds'][1:-1]
        if t<tN-1:
            bounds+=', '
        if t==tN-1:
            bounds+=']'
    if 'Z' in axes:
        tif.imwrite(outname, IMG, imagej=True,
            resolution=(res[0][0]/float(res[0][1]),res[1][0]/float(res[1][1])), 
            metadata={'unit': 'nm', 'finterval':fpns, 'funit': 'ns', 
                      'bounds': bounds, 'pbc':metaD['pbc'], 'axes':'T'+axes, 
                      'spacing': metaD['spacing']})
    else:
        tif.imwrite(outname, IMG, imagej=True,
            resolution=(res[0][0]/float(res[0][1]),res[1][0]/float(res[1][1])), 
            metadata={'unit':'nm', 'finterval':fpns, 'bounds':bounds,
                      'funit':'ns', 'pbc':metaD['pbc'], 'axes':'T'+axes})

def tiff2float(f):
    foo=tif.imread(f)
    dtype=foo.dtype
    foo=foo.astype('float')
    foo=foo/np.iinfo(dtype).max
    return foo

def imgs2color(datafile,outname,mix_type,lam_hues):

    f=open(datafile,'r')
    data=[]
    for lines in f:
        foo=lines.strip()
        data.append(foo)
    f.close()
    IMGs=[]
    cols=3
    if mix_type=='nomix':
        cols=len(data)
    img_data=tif.TiffFile(data[0])
    dtype=img_data.series[0].dtype
    print("Color mixing can lose precision!")
    img_data=tif.TiffFile(data[0])
    res=[img_data.pages[0].tags['XResolution'].value,
         img_data.pages[0].tags['YResolution'].value]
    metaD=img_data.imagej_metadata
    metaD2={}
    for a in ['unit', 'spacing', 'finterval', 'funit', 'pbc', 'bounds']:
        if a in metaD:
            metaD2[a]=metaD[a]
    axes=img_data.series[0].axes

    #Convert TIFF to float type.

    cpus=mp.cpu_count()
    if len(data)<cpus:
        cpus=len(data)
    pool=mp.Pool(cpus)
    IMGs=pool.map(tiff2float,data) # Lams (TZYX/ TYX/ ZYX or YX)
    pool.close()
    #Determine new shape and axes of the combined image.
    sha=IMGs[0].shape
    IMGs=np.array(IMGs) #C, Axes.
    if axes[0:2]=='TZ':
        metaD2['axes']=axes[0:2]+'C'+axes[2:]
        new_sha=tuple(list(sha[0:2])+[cols]+list(sha[2:]))
    elif axes[0]=='T' or axes[0]=='Z':
        metaD2['axes']=axes[0]+'C'+axes[1:]
        new_sha=tuple([sha[0],cols]+list(sha[1:]))
    else:
        metaD2['axes']='C'+axes
        new_sha=tuple([cols]+list(sha[:]))

    if mix_type=='nomix':
        if axes=='TZYX':
            colIMG=np.transpose(IMGs,(1,2,0,3,4))
        elif axes=='TYX' or axes=='ZYX':
            colIMG=np.transpose(IMGs,(1,0,2,3))
    else:  #mt , rgb
        colIMG=np.zeros(new_sha)
        if axes=='TZYX':
            Arguments=[]                
            for t in range(new_sha[0]):
                for z in range(new_sha[1]):
                    Arguments.append([np.transpose(IMGs[:,t,z,:,:],(2,1,0)),lam_hues,0,mix_type]) #XYC
            cpus=mp.cpu_count()
            if len(Arguments)<cpus:
                cpus=len(Arguments)
            pool=mp.Pool(cpus)
            results=pool.starmap(plot_image.add_color,Arguments) #XYC
            pool.close()

            for t in range(new_sha[0]):
                for z in range(new_sha[1]):
                    i=t*new_sha[1]+z
                    foo=np.transpose(results[i],(2,1,0)) #CYX
                    colIMG[t,z,:,:,:]=foo[:,:,:]
        elif axes=='TYX' or axes=='ZYX':
            Arguments=[]                
            for z in range(new_sha[1]):
                Arguments.append([np.transpose(IMGs[:,z,:,:],(2,1,0)),lam_hues,0,mix_type]) #XYC
            cpus=mp.cpu_count()
            if len(Arguments)<cpus:
                cpus=len(Arguments)
            pool=mp.Pool(cpus)
            results=pool.starmap(plot_image.add_color,Arguments) #XYC
            pool.close()

            for z in range(new_sha[1]):
                foo=np.transpose(results[z],(2,1,0)) #Converts to CYX
                colIMG[z,:,:,:]=foo[:,:,:]

        elif axes=='YX':
            foo=np.transpose(IMGs[:,:,:],(2,1,0)) #Converts to XYC
            foo=plot_image.add_color(foo,lam_hues,0,mix_type) 
            foo=np.transpose(foo,(2,1,0)) #Converts to CYX
            colIMG[:,:,:]=foo[:,:,:]
    if 'Z' in axes:
        metaD2['spacing']=metaD['spacing']
    if 'T' in axes:
        metaD2['finterval']=metaD['finterval']
    colIMG,bounds=plot_image.intensity2image(colIMG,dtype,metaD2['axes'])
    tif.imwrite(outname, colIMG, imagej=True, metadata=metaD2, \
        resolution=(res[0][0]/float(res[0][1]),res[1][0]/float(res[1][1])))
