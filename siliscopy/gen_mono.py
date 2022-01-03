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

import os
import multiprocessing as mp
import numpy as np
from . import convert
def gen_mono_c(data,silent=False,photophys=False):
    """ Runs the C-binary to calculate monochome intensities
 
    Parameters
    ----------
    data: array for strings
        data[0] is Gro filename, data[1] is parameter filename, data[2] is PSF 
        file name header, and data[3] is output file name header.
    silent: bool
        True would suppress the default output of gen_mono C binary file.
        (default False)

    Writes
    ------
    [data[3]]_lam[lam[i]]_fs[fs].dat: custom data format.
        [lam[i]] and [fs] is read from paramter file. [i] is an integers 
        greater than 1. Each line contains intensities separated by spaces.
        Intensity of -1 represents the white frame (absence of molecular
        simulation system).
    """
    mono="gen_mono"
    if photophys==True:
        mono="gen_mono_pp"
    if silent:
        os.system(os.path.dirname(__file__) + '/'+mono+' -f ' + data[0] + ' -p ' + 
              data[1]+' -psf '+data[2]+' -o '+data[3] +' -silent')
        print('Writing: '+data[3]+'         ',end='\r')
    else:
        os.system(os.path.dirname(__file__) + '/'+mono+' -f ' + data[0] + ' -p ' + 
              data[1]+' -psf '+data[2]+' -o '+data[3])
    

def read_data(datafile):
    """ Reads comma separated data from a file.
  
    Parameters
    ----------
    datafile: str
        Data file name
    """
    Arguments=[]
    f=open(datafile,'r')
    for lines in f: #Read data 
        foo=lines.split(',')
        Arguments.append(foo)
    f.close()
    return Arguments 

def gen_mono_c_mp(datafile,silent,photophys):
    """ Runs gen_mono_c using multiprocessing. Reads the required filenames 
        from a file.

    Parameters
    ----------
    datafile: str
        Data file name

    Writes
    ------
    multiple files similar to gen_mono_c
    """
    Arguments=read_data(datafile)
    Arg_slice=[]
    Arg_vol=[]
    for i in range(len(Arguments)):
        if Arguments[i][-1].strip()=='volume':
            Arg_vol.append([Arguments[i],silent,photophys])
        else: #slice
            Arg_slice.append([Arguments[i],silent,photophys])
    pool=mp.Pool(mp.cpu_count())
    results=pool.starmap(gen_mono_c,Arg_slice)
    
    for i in range(len(Arg_vol)):
        gen_mono_c_vol(*Arg_vol[i])

def gen_mono_c_serial(datafile, silent, photophys):
    """ Runs gen_mon_c serially multiple times, using the filenames acquired
        from a file.

    Parameters
    ----------
    datafile: str
        Data file name

    Writes
    ------
    multiple files similar to gen_mono_c
    """
    Arguments=read_data(datafile)
    for i in range(len(Arguments)):
        if Arguments[i][-1]=='volume':
            gen_mono_c_vol(Arguments[i], silent, photophys)
        else:#slice
            gen_mono_c(Arguments[i], silent, photophys)

def gen_mono_c_vol(data, silent, photophys, maxlen=None, opt_axis=None, dlmn=None, add_n=1,
    mprocess=True):
    """ Calculates image intensity for multiple slices, which is equivalent
        to a 3D image.

    Parameters
    ----------
    data: list of str
        A list of string containing input filename, parameter filename, psf 
        filename header, and output filename header.
    maxlen: list of floats
        Maximum box length of the system. (default None)
    opt_axis: int
        Optical axis for in-silico microscope (default None)
    dlmn: list of float
        Voxel dimensions (default None)
    add_n: int
        There are maxlen[opt_axis]/dlmn[2] volume slices. Every `add_n` slices 
        is calculated. (default 1)
    mprocess: bool
        If true multiprocessing is used to calculate the 3D image intensity.
        (default True)
    """

    xyz='xyz'
    if maxlen is None:
        lams=np.zeros(10,dtype='int')
        tsO=None
        f=open(data[1],'r')
        for lines in f:
            foo=lines.split('=')
            varname=foo[0].strip()
            val_string=foo[1].split('/')[0].strip()
            val_string=val_string.split()
            if varname == 'maxlen':
                maxlen=[float(x) for x in val_string]
            if varname == 'dlmn':
                dlmn=[float(x) for x in val_string]
            if varname == 'opt_axis':
                opt_axis=int(val_string[0])
            if varname == 'add_n':
                add_n=int(val_string[0])
            if varname == 'psf_type':
                psf_type=int(val_string[0])
            if varname[0:3] == 'lam':
                if varname[3].isdigit():
                    num=int(varname[3:])
                    lams[num-1]=int(val_string[0])
            if varname == 'tsO':
                tsO=float(val_string[0])
            if varname == 'fs':
                fs=int(val_string[0])
        f.close()

    f1=open(data[0],'r')
    for lines in f1:
        foo=lines.split()
    f1.close()
    box=[float(x) for x in foo]
 
    N=int(maxlen[opt_axis]/(dlmn[2]*add_n)+0.5)
    N1=int(box[opt_axis]/(dlmn[2]*add_n)+0.5)
    for n in range(0,int((N-N1)/2)):
        #-1 image
        for i in range(len(lams)):
            if lams[i]==0:
                break
            fend="_lam"+str(lams[i])+"_fs"+str(fs)+".dat"
            if psf_type==1:
                fend="_tsO%g"%tsO+fend
            white_image(data[3]+'_'+xyz[opt_axis]+str(n*add_n)+fend,maxlen,dlmn,opt_axis)
    for n in range(N1+int((N-N1)/2),N):
        #-1 image
        for i in range(len(lams)):
            if lams[i]==0:
                break
            fend="_lam"+str(lams[i])+"_fs"+str(fs)+".dat"
            if psf_type==1:
                fend="_tsO%g"%tsO+fend
            white_image(data[3]+'_'+xyz[opt_axis]+str(n*add_n)+fend,maxlen,dlmn,opt_axis)

 
    w=open(data[0]+'datalist.dat','w')
    for i in range(N1):
        n=round(i*add_n*dlmn[2],4)
        os.system("sed 's/focus_cor.*/focus_cor = "+str(n)+"/g' "+data[1]+" > " + \
            data[0]+"foo_param"+str(i*add_n) + ".dat")
        w.write(data[0] + ','+str(data[0])+'foo_param' + str(i*add_n) + '.dat,' + data[2] + ',' + \
            data[3]+ '_'+ xyz[opt_axis] + str((i+int((N-N1)/2))*add_n) + ',slice\n')
    w.close()
    if mprocess==True:
        gen_mono_c_mp(data[0]+'datalist.dat', silent, photophys)
    else:
        gen_mono_c_serial(data[0]+'datalist.dat', silent, photophys)
 
    #Cleanup
    os.system('rm '+str(data[0])+'foo_param*.dat')
    os.system('rm '+str(data[0])+'datalist.dat')


def white_image(outname,maxlen,dlmn,opt_axis):
    w=open(outname,'w')
    w.write('# White Frame for N axis\n')
    xN=int(maxlen[(opt_axis+1)%3]/dlmn[0]+0.5)
    yN=int(maxlen[(opt_axis+2)%3]/dlmn[1]+0.5)    
    for i in range(xN):
        for j in range(yN):
           w.write('-1.0 ')
        w.write('\n')
    w.close()



