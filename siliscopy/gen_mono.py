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
def gen_mono_c(data):
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
    pool=mp.Pool(mp.cpu_count())
    results=pool.map(gen_mono_c,Arguments)

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
        gen_mono_c(Arguments[i])
