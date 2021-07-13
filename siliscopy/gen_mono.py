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
from . import convert
def gen_mono_c(data,silent=False):
    """ Runs the C-binary to calculate monochome intensities
 
    Parameters
    ----------
    data: array for strings
        data[0] is Gro filename, data[1] is parameter filename, data[2] is PSF 
        file name header, and data[3] is output file name header.

    Writes
    ------
    [data[3]]_lam[lam[i]]_fs[fs].dat: custom data format.
        [lam[i]] and [fs] is read from paramter file. [i] is an integers 
        greater than 1. Each line contains intensities separated by spaces.
        Intensity of -1 represents the white frame (absence of molecular
        simulation system).
    """
    if silent:
        os.system(os.path.dirname(__file__) + '/gen_mono -f ' + data[0] + ' -p ' + 
              data[1]+' -psf '+data[2]+' -o '+data[3] +' 1> /dev/null')
        print('Writing: '+data[3]+'         ',end='\r')
    else:
        os.system(os.path.dirname(__file__) + '/gen_mono -f ' + data[0] + ' -p ' + 
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
    return Arguments 

def gen_mono_c_mp(datafile):
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
            Arg_vol.append(Arguments[i])
        else: #slice
            Arg_slice.append([Arguments[i],True])

    pool=mp.Pool(mp.cpu_count())
    results=pool.starmap(gen_mono_c,Arg_slice)

    for i in range(len(Arg_vol)):
        gen_mono_c_vol(Arg_vol[i])

def gen_mono_c_serial(datafile):
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
            gen_mono_c_vol(Arguments[i])
        else:#slice
            gen_mono_c(Arguments[i])

def gen_mono_c_vol(data, maxlen=None, opt_axis=None, dlmn=None, add_n=1,
    mprocess=True):
    xyz='xyz'
    if maxlen==None:
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
        f.close()

    N=int(maxlen[opt_axis]/dlmn[2]+1E-3)
    N1=int(N/add_n +1E-3)
    w=open(data[0]+'datalist.dat','w')
    for i in range(N1):
        n=round(i*add_n*dlmn[2],4)
        os.system("sed 's/focus_cor.*/focus_cor = "+str(n)+"/g' "+data[1]+" > " + \
            data[0]+"foo_param"+str(i*add_n) + ".dat")
        w.write(data[0] + ','+str(data[0])+'foo_param' + str(i*add_n) + '.dat,' + data[2] + ',' + \
            data[3]+ '_'+ xyz[opt_axis] + str(i*add_n) + ',slice\n')
    w.close()
    if mprocess==True:
        gen_mono_c_mp(data[0]+'datalist.dat')
    else:
        gen_mono_c_serial(data[0]+'datalist.dat')
 
    #Cleanup
    os.system('rm '+str(data[0])+'foo_param*.dat')
    os.system('rm '+str(data[0])+'datalist.dat')

