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
