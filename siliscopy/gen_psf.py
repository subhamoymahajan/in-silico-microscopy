#    Copyright 2020 SUBHAMOY MAHAJAN 
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

# This code generates a point spread function for a in-silico microscope at 3D 
# grid points

import scipy.integrate as integrate
import scipy.special as special
import numpy as np
import multiprocessing as mp
import os 

def psf_gandy(beta, lambd, dlmn, Plmn, fs, outname):
    """ Calculates Point Spread Function based on R. O. Gandy, 1954 Proc. Phys.
        Soc. B 67 825-831. 
    
        Calculates PSF for 0 < l' < Pl'/2 , 0 < m' <= l', 0 < n' < Pn'/2.
    
        Parameters
        ----------
        beta: float
            Maximum half angle as seen from the immersion oil.
        lambd: int
            Wavelenght of fluorescence light
        dlmn: array of floats
            Delta l', Delta m', and Delta n'. Voxel size.
        Plmn: array of floats
            Dimensions of the box withing which PSF is calculated.
        fs: int
            Scaling factor for wavenumber or position coordinates
        outname: str
            Output filename header.
      
        Writes
        ------
        [outname]_lam[lambd]_fs[fs].dat: custom data file.
            Contains l', m', n' and PSF(l',m',n') in each line.
    """
    Nn=int((int(Plmn[2]/dlmn[2])+1)/2)+1
    outname = outname + '_lam' + str(lambd) + '_fs' + str(fs) + '.dat'
    for k in range(Nn):
        psf_gandy_sep(beta, lambd, dlmn, Plmn, fs, outname, 'a', k)       

def psf_gandy_sep(beta, lambd, dlmn, Plmn, fs, outname, otype, nidx):
    """ Calculates Point Spread Function based on R. O. Gandy, 1954 Proc. Phys.
        Soc. B 67 825-831 for a specific n' coordinate.
      
        Parameters
        ----------
        beta: float
            Maximum half angle as seen from the immersion oil.
        lambd: int
            Wavelenght of fluorescence light
        dlmn: array of floats
            Delta l', Delta m', and Delta n'. Voxel size.
        Plmn: array of floats
            Dimensions of the box withing which PSF is calculated.
        fs: int
            Scaling factor for wavenumber or position coordinates
        outname: str
            Output filename header.
        otype: char
            'w' to create a new file. 'a' to append to an existing file
        nidx: int
            Index for the n' coordinate. nidx=int(n'/dn').
    
        Writes
        ------
        [outname]:
            Each line contains l', m', n' and PSF(l',m',n').
    """
    def integrand_real(x, r, d):
         a=np.sin(x)
         b=np.cos(x)
         return a*(b**0.5)*np.cos(d*b)*special.jv(0,r*a)
     
    def integrand_img(x, r, d):
         a=np.sin(x)
         b=np.cos(x)
         return a*(b**0.5)*np.sin(d*b)*special.jv(0,r*a)
     
    A=3/(2*(1-(np.cos(beta)**1.5)))
    K=2*np.pi*fs/lambd 
    
    w=open(outname, otype)
    if otype=='w':
        w.write('# Beta= ' + str(beta) + ' lambda= ' + str(lambd) + \
                ' dlmn= ' + str(dlmn[0]) + ' ' + str(dlmn[1]) + ' ' + \
                str(dlmn[2]) + ' Plmn= ' + str(Plmn[0]) + ' ' + \
                str(Plmn[1]) + ' ' + str(Plmn) + ' fs= ' + str(fs) + '\n')
 
    Nl=int((int(Plmn[0]/dlmn[0])+1)/2)+1
    Nm=int((int(Plmn[0]/dlmn[1])+1)/2)+1
    n=round(nidx*dlmn[2],6)
    Kn=K*n
    for i in range(Nl):
        l=round(i*dlmn[0],6)
        for j in range(i+1):
            m=round(j*dlmn[1],6)
            r=np.sqrt(l*l+m*m)
            Kr=K*r
            H1=integrate.quad(integrand_real, 0, beta, args=(Kr, Kn))
            H2=integrate.quad(integrand_img, 0, beta, args=(Kr, Kn))
            I=round((H1[0]*H1[0]+H2[0]*H2[0])*A*A,6)
            w.write(str(l).rjust(10) + str(m).rjust(10) + str(n).rjust(10) + \
                    str(I).rjust(10) + '\n')
    print('n = '+str(n)+' done')

def worker_gandy(data):
    """Runs psf_gandy_sep. It is used by psf_gandy_mp.
 
    See psf_gandy_sep for more details.
    """
    #              beta     lambd    dlmn     Plmn      fs     outname   otype       
    psf_gandy_sep(data[0], data[1], data[2], data[3], data[4], data[5], data[6], 
                  data[7])
                  # nidx

def psf_gandy_mp(beta, lambd, dlmn, Plmn, fs, outname):
    """ Generates PSF using multiprocssing.
    
    Parameters
    ----------
    beta: float
        Maximum half angle as seen from the immersion oil.
    lambd: int
        Wavelenght of fluorescence light
    dlmn: array of floats
        Delta l', Delta m', and Delta n'. Voxel size.
    Plmn: array of floats
        Dimensions of the box withing which PSF is calculated.
    fs: int
        Scaling factor for wavenumber or position coordinates
    outname: str
        Output filename header.
    
    Writes
    ------
    [outname]_lam[lambd]_fs[fs].dat
        Each line contains l', m', n', and PSF(l',m',n')
    """
    Nn=int((int(Plmn[2]/dlmn[2])+1)/2)+1
    Arguments=[]
    for k in range(Nn):
        Arguments.append([beta, lambd, dlmn, Plmn, fs, outname+'n'+str(k), 'w', 
                          k])
    pool=mp.Pool(mp.cpu_count())
    results=pool.map(worker_gandy,Arguments)
    os.system('mv ' + outname + 'n0 ' + outname + '_lam' + str(lambd) + \
              '_fs' + str(fs) + '.dat')
    for k in range(1,Nn):
        os.system('tail -n +2 ' + outname + 'n' + str(k) + ' >> ' + outname + \
                  '_lam' + str(lambd) + '_fs' + str(fs) + '.dat')
        os.system('rm '+ outname + 'n' + str(k))

#Ref: https://www.micro-shop.zeiss.com/en/us/shop/objectives/420493-9900-000/Objective-EC-Plan-Neofluar-100x-1.3-Oil-Pol-M27#/
#Lens: Plan-Neofluar 100x/1.3
#NA=1.3
#Refractive index on immersion oil = 1.51

