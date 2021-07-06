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

# This code generates a point spread function for a in-silico microscope at 3D 
# grid points

import scipy.integrate as integrate
import scipy.special as special
import numpy as np
import multiprocessing as mp
import os 


#Ref: https://www.micro-shop.zeiss.com/en/us/shop/objectives/420493-9900-000/Objective-EC-Plan-Neofluar-100x-1.3-Oil-Pol-M27#/
#Lens: Plan-Neofluar 100x/1.3
#NA=1.3
#Refractive index on immersion oil = 1.51

def psf_gandy(NA, meu, lambd, dlmn, Plmn, fs, outname):
    """ Calculates Point Spread Function based on R. O. Gandy, 1954 Proc. Phys.
        Soc. B 67 825-831. 
    
        Calculates PSF for 0 < l' < Pl'/2 , 0 < m' <= l', 0 < n' < Pn'/2.
    
        Parameters
        ----------
            See psf_gandy_sep() 

        Writes
        ------
        [outname]_lam[lambd]_fs[fs].dat: custom data file.
            Contains l', m', n' and PSF(l',m',n') in each line.
    """
    Nn=int((int(Plmn[2]/dlmn[2])+1)/2)+1
    outname = outname + '_lam' + str(lambd) + '_fs' + str(fs) + '.dat'
    for k in range(Nn):
        psf_gandy_sep(NA, meu, lambd, dlmn, Plmn, fs, outname, 'a', k)       

def psf_gandy_sep(NA, meu, lambd, dlmn, Plmn, fs, outname, otype, nidx):
    """ Calculates Point Spread Function based on R. O. Gandy, 1954 Proc. Phys.
        Soc. B 67 825-831 for a specific n' coordinate.
      
        Parameters
        ----------
        NA: float
            Numerical aperture.
        meu: float
            Refractive index of immersion oil.
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
    #d= k'\mu n' or -k'\mu n'. The negative sign is inconsequential.
    def integrand_real(x, r, d):
         a=np.sin(x)
         b=np.cos(x)
         return a*(b**0.5)*np.cos(d*b)*special.jv(0,r*a)
     
    def integrand_img(x, r, d):
         a=np.sin(x)
         b=np.cos(x)
         return a*(b**0.5)*np.sin(d*b)*special.jv(0,r*a)
    sinb=NA/meu
    beta=np.arcsin(sinb)
    cos2b=1-sinb**2
    A=3/(2*(1-cos2b**0.75))
    K=2*np.pi*fs*meu/lambd 
    
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
    w.close()

def worker_gandy(data):
    """Runs psf_gandy_sep. It is used by psf_gandy_mp.
 
    See psf_gandy_sep for more details.
    """
    #              NA       meu      lambd    dlmn     Plmn      fs     outname          
    psf_gandy_sep(data[0], data[1], data[2], data[3], data[4], data[5], data[6], 
                  data[7], data[8])
                  # otype  nidx

def psf_gandy_mp(NA, meu, lambd, dlmn, Plmn, fs, outname):
    """ Generates PSF using multiprocssing.
    
    Parameters
    ----------
        See psf_gandy_sep() 
    Writes
    ------
    [outname]_lam[lambd]_fs[fs].dat
        Each line contains l', m', n', and PSF(l',m',n')
    """
    Nn=int((int(Plmn[2]/dlmn[2])+1)/2)+1
    Arguments=[]
    for k in range(Nn):
        Arguments.append([NA, meu, lambd, dlmn, Plmn, fs, outname+'n'+str(k), 'w', 
                          k])
    pool=mp.Pool(mp.cpu_count())
    results=pool.map(worker_gandy,Arguments)
    os.system('mv ' + outname + 'n0 ' + outname + '_lam' + str(lambd) + \
              '_fs' + str(fs) + '.dat')
    for k in range(1,Nn):
        os.system('tail -n +2 ' + outname + 'n' + str(k) + ' >> ' + outname + \
                  '_lam' + str(lambd) + '_fs' + str(fs) + '.dat')
        os.system('rm '+ outname + 'n' + str(k))


def psf_GL1991(NA, meu, meu0, t0, tsO, meus, tg, tg0, meug, meug0, lambd, dlmn,\
    Plmn, fs, outname):
    """ See function psf_GL1991_sep()
    
        Parameters
        ----------
        See function psf_GL1991_sep()
 
        Writes
        ------
        [outname]_lam[lambd]_fs[fs].dat: custom data file.
            Contains l', m', n' and PSF(l',m',n') in each line.
    """
    Nn=int((int(Plmn[2]/dlmn[2])+1)/2)
    outname = outname + '_tsO' + str(tsO) + '_lam' + str(lambd) + '_fs' + \
        str(fs) + '.dat'
    print(Nn) 
    for k in range(-Nn,Nn+1):
        psf_GL1991_sep(NA, meu, meu0, t0, tsO, meus, tg, tg0, meug, meug0, \
            lambd, dlmn, Plmn, fs, outname, 'a', k)       

def psf_GL1991_sep(NA, meu, meu0, t0, tsO, meus, tg, tg0, meug, meug0, lambd,\
    dlmn, Plmn, fs, outname, otype, nidx):
    """ Calculates Point Spread Function based on Gibson et al., J. Opt. Soc. 
        Am. A 1991, 8, 1601-1613. The same publication with better quality 
        image can be found in Gibson et al., J. Opt. Soc. Am. A. 1992, 9, 154-
        166.  

        The microscope is assumed to be in best geometric focus, based on the 
        2nd order approximation.
    
        Calculates PSF for 0 < l' < Pl'/2 , 0 < m' <= l', -Pn'/2 < n' < Pn'/2.
        Unlike Gandy Gibson and Lanni PSF is not symmetric about the object 
        focal plane
    
        Parameters
        ----------
        NA: float
            Numerical Aperture of objective lens
        meu: float
            Refractive index of immersion medium
        meu0: float
            Refractive index of immersion medium in design condition
        t0: float
            Thickness of immersion medium in design condition (in micrometers)
        tsO: float 
            Depth of specimen at which the microscoped is focused 
            (in micrometers)
        meus: float 
            Effective refractive index of specimen
        tg: float
            Thickness of coverslip (in micrometers)
        tg0: float
            Thickness of coverslip in design condition (in micrometers)
        meug: float 
            Refractive index of coverslip 
        meug0: float 
            Refractive index of coverslip in design condition
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
    def integrand_r(rho, kr, z, k, NA, meu, meu0, t0, tsO, meus, tg, tg0, \
        meug, meug0): 
        #kr and k are both provided to reduce computational cost.
        #tsO is focal plane distance from coverslip.
        NArho=NA*rho
        NArho2=NArho**2  #To reduce computational cost
        cmeu=np.sqrt(meu**2-NArho2) # cos(\theta_meu)
        OPD=-z*np.sqrt(meus**2-NArho2) 
        #OPD is in micrometers. z is in nanometers
        if tsO>1E-6:
            OPD+=tsO*(np.sqrt(meus**2-NArho2) - meu/meus*cmeu)
        if abs(tg-tg0)>1E-6 or abs(meug-meug0)>1E-6:
            OPD+=tg*(np.sqrt(meug**2-NArho2) - meu/meug*cmeu)
            OPD-=tg0*(np.sqrt(meug0**2-NArho2) - meu/meug0*cmeu)
        if abs(meu-meu0)>1E-6:
            OPD-=t0*(np.sqrt(meu0**2-NArho2) - meu/meu0*cmeu)
        return special.jv(0,kr*NArho)*np.cos(k*OPD)*rho

    def integrand_i(rho, kr, z, k, NA, meu, meu0, t0, tsO, meus, tg, tg0, \
        meug, meug0): 
        #kr and k are both provided to reduce computational cost.
        #tsO is focal plane distance from coverslip. 
        NArho=NA*rho
        NArho2=NArho**2  #To reduce computational cost
        cmeu=np.sqrt(meu**2-NArho2) # cos(\theta_meu)
        OPD=-z*np.sqrt(meus**2-NArho2)
        #OPD is in micrometers. z is in nanometers
        if tsO>1E-6:
            OPD+=tsO*(np.sqrt(meus**2-NArho2) - meu/meus*cmeu)
        if abs(tg-tg0)>1E-6 or abs(meug-meug0)>1E-6:
            OPD+=tg*(np.sqrt(meug**2-NArho2) - meu/meug*cmeu)
            OPD-=tg0*(np.sqrt(meug0**2-NArho2) - meu/meug0*cmeu)
        if abs(meu-meu0)>1E-6:
            OPD-=t0*(np.sqrt(meu0**2-NArho2) - meu/meu0*cmeu)
        return special.jv(0,kr*NArho)*np.sin(k*OPD)*rho
     
     
    K=2*np.pi*fs/lambd #in 1/nanometers
    
    w=open(outname, otype)
    if otype=='w':
        w.write('# NA= '+str(NA) + ', lambda= '+str(lambd)+' nm' + ', dlmn= ('+\
                str(dlmn[0])+','+str(dlmn[1])+','+str(dlmn[2])+') nm' + \
                ', Plmn= ('+str(Plmn[0])+','+str(Plmn[1])+','+str(Plmn[2])+') nm'+\
                ', fs= '+str(fs) + ', meu= '+str(meu) + ', meu0= '+str(meu0) +\
                ', t0= '+str(t0)+' nm' + ', tsO= '+str(tsO)+' nm' + ', meus= '+\
                str(meus) + ', tg= '+str(tg)+' nm' + ', tg0= '+str(tg0)+' nm' +\
                ', meug= '+str(meug) + ', meug0= ' + str(meug0) + '\n')
    
    Nl=int((int(Plmn[0]/dlmn[0])+1)/2)+1
    Nm=int((int(Plmn[1]/dlmn[1])+1)/2)+1
    n=round(nidx*dlmn[2],6)

    for i in range(Nl):
        l=round(i*dlmn[0],6)
        for j in range(i+1):
            m=round(j*dlmn[1],6)
            r=np.sqrt(l*l+m*m) #l,m,n is in nanometers
            Kr=K*r #dimensionless.
            H1=integrate.quad(integrand_r, 0, 1, args=(Kr, n, K, NA, meu, meu0, \
                t0, tsO, meus, tg, tg0, meug, meug0), limit=150)
            H2=integrate.quad(integrand_i, 0, 1, args=(Kr, n, K, NA, meu, meu0, \
                t0, tsO, meus, tg, tg0, meug, meug0), limit=150)
            I=H1[0]*H1[0]+H2[0]*H2[0]
            I=I*4.0
            w.write(str(l).rjust(10) + str(m).rjust(10) + str(n).rjust(10) + \
                    ' ' + str(I) + '\n')
    print('n = '+str(n)+' done')
    w.close()

def psf_GL1991_mp(NA, meu, meu0, t0, tsO, meus, tg, tg0, meug, meug0, lambd,\
    dlmn, Plmn, fs, outname):
    """ Generates Gibson and Lanni PSF using multiprocssing. 
        See psf_GL1991_sep() for more details
    Parameters
    ----------
        See psf_GL1991_sep() for more details
    
    Writes
    ------
    [outname]_lam[lambd]_fs[fs].dat
        Each line contains l', m', n', and PSF(l',m',n')
    """
    Nn=int((int(Plmn[2]/dlmn[2])+1)/2)
    Arguments=[]
    for k in range(-Nn,Nn+1):
        Arguments.append([NA, meu, meu0, t0, tsO, meus, tg, tg0, meug, meug0,
            lambd, dlmn, Plmn, fs, outname+'n'+str(k), 'w', k])
    pool=mp.Pool(mp.cpu_count())
    results=pool.starmap(psf_GL1991_sep,Arguments)
    os.system('mv ' + outname+'n'+str(-Nn)+' ' + outname+'_tsO'+str(tsO)+'_lam'+str(lambd)+ \
              '_fs'+str(fs)+'.dat')
    for k in range(-Nn+1,Nn+1):
        os.system('tail -n +2 ' + outname+'n'+str(k) + ' >> ' + outname+'_tsO'+ \
                  str(tsO)+'_lam'+str(lambd)+'_fs'+str(fs)+'.dat')
        os.system('rm '+ outname+'n'+str(k))

def psf_Mod_Gandy(NA, meu, meu0, t0, tsO, meus, tg, tg0, meug, meug0, lambd, \
    dlmn, Plmn, fs, outname):
    """ See function psf_Mod_Gandy_sep()
    
        Parameters
        ----------
        See function psf_Mod_Gandy_sep()
 
        Writes
        ------
        [outname]_lam[lambd]_fs[fs].dat: custom data file.
            Contains l', m', n' and PSF(l',m',n') in each line.
    """
    Nn=int((int(Plmn[2]/dlmn[2])+1)/2)
    outname = outname + '_tsO' + str(tsO) + '_lam' + str(lambd) + '_fs' + \
        str(fs) + '.dat'
    print(Nn) 
    for k in range(-Nn,Nn+1):
        psf_Mod_Gandy_sep(NA, meu, meu0, t0, tsO, meus, tg, tg0, meug, meug0,\
             lambd, dlmn, Plmn, fs, outname, 'a', k)       

def psf_Mod_Gandy_sep(NA, meu, meu0, t0, tsO, meus, tg, tg0, meug, meug0, \
    lambd, dlmn, Plmn, fs, outname, otype, nidx):
    """ Calculates Point Spread Function based on Gandy 
        modified using Gibson et al., J. Opt. Soc. Am. A 1991, 8, 1601-1613. The
        same publication with better quality image can be found in Gibson et al.,
         J. Opt. Soc. Am. A. 1992, 9, 154-166.  

        The microscope is assumed to be in best geometric focus, based on the 
        2nd order approximation.
    
        Calculates PSF for 0 < l' < Pl'/2 , 0 < m' <= l', -Pn'/2 < n' < Pn'/2.
        Unlike Gandy Gibson and Lanni PSF is not symmetric about the object 
        focal plane
    
        Parameters
        ----------
        NA: float
            Numerical Aperture of objective lens
        meu: float
            Refractive index of immersion medium
        meu0: float
            Refractive index of immersion medium in design condition
        t0: float
            Thickness of immersion medium in design condition (in micrometers)
        tsO: float 
            Depth of specimen at which the microscoped is focused 
            (in micrometers)
        meus: float 
            Effective refractive index of specimen
        tg: float
            Thickness of coverslip (in micrometers)
        tg0: float
            Thickness of coverslip in design condition (in micrometers)
        meug: float 
            Refractive index of coverslip 
        meug0: float 
            Refractive index of coverslip in design condition
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
    def integrand_r(rho, kr, z, k, NA, meu, meu0, t0, tsO, meus, tg, tg0, \
        meug, meug0): 
        #kr and k are both provided to reduce computational cost.
        #tsO is focal plane distance from coverslip.
        NArho=NA*rho
        NArho2=NArho**2  #To reduce computational cost
        cmeu=np.sqrt(meu**2-NArho2) # cos(\theta_meu)
        OPD=-z*np.sqrt(meus**2-NArho2)
        #OPD is in micrometers. z is in nanometers
        if tsO>1E-6:
            OPD+=tsO*(np.sqrt(meus**2-NArho2) - meu/meus*cmeu)
        if abs(tg-tg0)>1E-6 or abs(meug-meug0)>1E-6:
            OPD+=tg*(np.sqrt(meug**2-NArho2) - meu/meug*cmeu)
            OPD-=tg0*(np.sqrt(meug0**2-NArho2) - meu/meug0*cmeu)
        if abs(meu-meu0)>1E-6:
            OPD-=t0*(np.sqrt(meu0**2-NArho2) - meu/meu0*cmeu)
        #meu^0.5 will be multiplied after integration.
        return special.jv(0,kr*NArho)*np.cos(k*OPD)*rho/np.sqrt(cmeu)

    def integrand_i(rho, kr, z, k, NA, meu, meu0, t0, tsO, meus, tg, tg0, \
        meug, meug0): 
        #kr and k are both provided to reduce computational cost.
        #tsO is focal plane distance from coverslip. 
        NArho=NA*rho
        NArho2=NArho**2  #To reduce computational cost
        cmeu=np.sqrt(meu**2-NArho2) # cos(\theta_meu)
        OPD=-z*np.sqrt(meus**2-NArho2)
        #OPD is in micrometers. z is in nanometers
        if tsO>1E-6:
            OPD+=tsO*(np.sqrt(meus**2-NArho2) - meu/meus*cmeu)
        if abs(tg-tg0)>1E-6 or abs(meug-meug0)>1E-6:
            OPD+=tg*(np.sqrt(meug**2-NArho2) - meu/meug*cmeu)
            OPD-=tg0*(np.sqrt(meug0**2-NArho2) - meu/meug0*cmeu)
        if abs(meu-meu0)>1E-6:
            OPD-=t0*(np.sqrt(meu0**2-NArho2) - meu/meu0*cmeu)
        #meu^0.5 will be multiplied after integration.
        return special.jv(0,kr*NArho)*np.sin(k*OPD)*rho/np.sqrt(cmeu)
     
     
    K=2*np.pi*fs/lambd #in 1/nanometers
    
    w=open(outname, otype)
    if otype=='w':
        w.write('# NA= '+str(NA) + ', lambda= '+str(lambd)+' nm' + ', dlmn= ('+ \
                str(dlmn[0])+','+str(dlmn[1])+','+str(dlmn[2])+ ') nm' + \
                ', Plmn= ('+str(Plmn[0])+','+str(Plmn[1])+','+str(Plmn) + ') nm' +\
                ', fs= '+str(fs) + ', meu= '+str(meu) + ', meu0= '+str(meu0) + \
                ', t0= '+str(t0)+' nm' + ', tsO= '+str(tsO)+' nm' + ', meus= '+ \
                str(meus) + ', tg= '+str(tg)+' nm,' + ' tg0= '+str(tg0)+' nm' + \
                ', meug= '+str(meug) + ', meug0= '+str(meug0) + '\n')
    
    Nl=int((int(Plmn[0]/dlmn[0])+1)/2)+1
    Nm=int((int(Plmn[1]/dlmn[1])+1)/2)+1
    n=round(nidx*dlmn[2],6)
    sinb=NA/meu
    sinb2=sinb**2
    cosb2=1-sinb2
    fac=(1.5*sinb2/(1-cosb2**0.75))**2*meu
    for i in range(Nl):
        l=round(i*dlmn[0],6)
        for j in range(i+1):
            m=round(j*dlmn[1],6)
            r=np.sqrt(l*l+m*m) #l,m,n is in nanometers
            Kr=K*r #dimensionless.
            H1=integrate.quad(integrand_r, 0, 1, args=(Kr, n, K, NA, meu, meu0, \
                t0, tsO, meus, tg, tg0, meug, meug0), limit=150)
            H2=integrate.quad(integrand_i, 0, 1, args=(Kr, n, K, NA, meu, meu0, \
                t0, tsO, meus, tg, tg0, meug, meug0), limit=150)
            I=H1[0]*H1[0]+H2[0]*H2[0]
            I=I*fac
            w.write(str(l).rjust(10) + str(m).rjust(10) + str(n).rjust(10) + \
                    ' ' + str(I) + '\n')
    print('n = '+str(n)+' done')
    w.close()

def psf_Mod_Gandy_mp(NA, meu, meu0, t0, tsO, meus, tg, tg0, meug, meug0, \
    lambd, dlmn, Plmn, fs, outname):
    """ Generates Gibson and Lanni PSF using multiprocssing. 
        See psf_Mod_Gandy_sep() for more details
    Parameters
    ----------
        See psf_Mod_Gandy_sep() for more details
    
    Writes
    ------
    [outname]_lam[lambd]_fs[fs].dat
        Each line contains l', m', n', and PSF(l',m',n')
    """
    Nn=int((int(Plmn[2]/dlmn[2])+1)/2)
    Arguments=[]
    for k in range(-Nn,Nn+1):
        Arguments.append([NA, meu, meu0, t0, tsO, meus, tg, tg0, meug, meug0, 
            lambd, dlmn, Plmn, fs, outname+'n'+str(k), 'w', k])
    pool=mp.Pool(mp.cpu_count())
    results=pool.starmap(psf_Mod_Gandy_sep,Arguments)
    os.system('mv ' + outname+'n'+str(-Nn)+' ' + outname+'_tsO'+str(tsO)+'_lam'+str(lambd)+ \
              '_fs'+str(fs)+'.dat')
    for k in range(-Nn+1,Nn+1):
        os.system('tail -n +2 ' + outname+'n'+str(k) + ' >> ' + outname+'_tsO'+ \
                  str(tsO)+'_lam'+str(lambd)+'_fs'+str(fs)+'.dat')
        os.system('rm '+ outname+'n'+str(k))

