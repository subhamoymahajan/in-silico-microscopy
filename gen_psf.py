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

# This code generates a point spread function for a in-silico microscope at 3D grid points

import scipy.integrate as integrate
import scipy.special as special
import numpy as np

def psf_gandy(beta,lambd,dl,dm,dn,Ll,Lm,Ln,fs,outname):
   Nn=int((int(Ln/dn)+1)/2)+1
   for k in range(Nn):
       psf_gandy_sep(beta,lambd,dl,dm,dn,Ll,Lm,Ln,fs,outname,'a',k)       

def psf_gandy_sep(beta,lambd,dl,dm,dn,Ll,Lm,Ln,fs,outname,otype,nidx):
# Based on  R O Gandy 1954 Proc. Phys. Soc. B 67 825-831
# NA = Numerical aperture
# meu = Refractive index of the immersion oil
# lambd = Wavelength of light
# dl, dm, and dn = Distance between consecutive grid points 
# Ln, Lm, and Ln = Dimensions of the PSF box
# fs = full-width-at-half-maximum (FWHM) scaling factor; scaling factor of wavenumber (k=2pi/lambda); scaling factor of "gro" coordinates. (New Paper Cite) 
# outname = the output file name
# otype = open a file outname to append or write ('a' or  'w')
# nidx = a counter for n coordinate;  int(n/dn)
   def integrand_real(x,r,d):
        a=np.sin(x)
        b=np.cos(x)
        return a*(b**0.5)*np.cos(d*b)*special.jv(0,r*a)
    
   def integrand_img(x,r,d):
        a=np.sin(x)
        b=np.cos(x)
        return a*(b**0.5)*np.sin(d*b)*special.jv(0,r*a)
    
   A=3/(2*(1-(np.cos(beta)**1.5)))
   K=2*np.pi*fs/lambd 
   
   w=open(outname,otype)
   if otype=='w':
       w.write('# NA= '+str(NA)+' meu= '+str(meu)+' lambda= '+str(lambd)+' dl= '+str(dl)+' dm= '+str(dm)+' dn= '+str(dn)+' Ll= '+str(Ll)+' Lm= '+str(Lm)+' Ln= '+str(Ln)+' fs= '+str(fs)+'\n')
   Nl=int((int(Ll/dl)+1)/2)+1
   Nm=int((int(Lm/dm)+1)/2)+1
   n=round(nidx*dn,6)
   Kn=K*n
   for i in range(Nl):
       l=round(i*dl,6)
       for j in range(i+1):
           m=round(j*dm,6)
           r=np.sqrt(l*l+m*m)
           Kr=K*r
           H1=integrate.quad(integrand_real,0,beta,args=(Kr,Kn))
           H2=integrate.quad(integrand_img,0,beta,args=(Kr,Kn))
           I=round((H1[0]*H1[0]+H2[0]*H2[0])*A*A,6)
           w.write(str(l).rjust(10)+str(m).rjust(10)+str(n).rjust(10)+str(I).rjust(10)+'\n')
   print('n = '+str(n)+' done')

#Ref: https://www.micro-shop.zeiss.com/en/us/shop/objectives/420493-9900-000/Objective-EC-Plan-Neofluar-100x-1.3-Oil-Pol-M27#/
#Lens: Plan-Neofluar 100x/1.3
#NA=1.3
#Refractive index on immersion oil = 1.51





