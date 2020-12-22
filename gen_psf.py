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

# Generates point spread function for a in-silico microscope

import scipy.integrate as integrate
import scipy.special as special
import numpy as np

def psf_gandy(NA,meu,lambd,dx,dy,dz,Lx,Ly,Lz,fs,outname):
   Nz=int((int(Lz/dz)+1)/2)+1
   for k in range(Nz):
       psf_gandy_sep(NA,meu,lambd,dx,dy,dz,Lx,Ly,Lz,fs,outname,'a',k)       

def psf_gandy_sep(NA,meu,lambd,dx,dy,dz,Lx,Ly,Lz,fs,outname,otype,zidx):
# Based on  R O Gandy 1954 Proc. Phys. Soc. B 67 825
# NA = numerical aperture
# meu = refractive index of the immersion oil
# lambd = wavelength of light
# dx, dy pixel dimensions. dz is the distance in the sectioned object plane (gro file)
# Lz, Ly, Lz is the dimension of the psf box
# fs = FWHM scaling factor; scaling factor of wavenumber (2pi/lambda); scaling factor of "gro" coordinates. (New Paper Cite) 
# outname = the output file name
# otype = open a file outname with type type ('a' or  'w')
# zidx = a counter for z coordinate
   def integrand_real(x,r,d):
        a=np.sin(x)
        ans=a*(1-a*a)**0.25
        ans=ans*np.cos(d*a)*special.jv(0,r*a)
        return ans
    
   def integrand_img(x,r,d):
        a=np.sin(x)
        ans=a*(1-a*a)**0.25
        ans=ans*np.sin(d*a)*special.jv(0,r*a)
        return ans
    
   beta=np.arcsin(NA/meu)
   A=3/(2*(1-(np.cos(beta)**1.5)))
   K=2*np.pi*fs/lambd 
   
   w=open(outname,otype)
   if otype=='w':
       w.write('# NA= '+str(NA)+' meu= '+str(meu)+' lambda= '+str(lambd)+' dx= '+str(dx)+' dy= '+str(dy)+' dz= '+str(dz)+' Lx= '+str(Lx)+' Ly= '+str(Ly)+' Lz= '+str(Lz)+' fs= '+str(fs)+'\n')
   Nx=int((int(Lx/dx)+1)/2)+1
   Ny=int((int(Ly/dy)+1)/2)+1
   z=round(zidx*dz,6)
   Kz=K*z
   for i in range(Nx):
       x=round(i*dx,6)
       for j in range(i+1):
           y=round(j*dy,6)
           r=np.sqrt(x*x+y*y)
           Kr=K*r
           H1=integrate.quad(integrand_real,0,beta,args=(Kr,Kz))
           H2=integrate.quad(integrand_img,0,beta,args=(Kr,Kz))
           I=round((H1[0]*H1[0]+H2[0]*H2[0])*A*A,6)
           w.write(str(x).rjust(10)+str(y).rjust(10)+str(z).rjust(10)+str(I).rjust(10)+'\n')
   print('z = '+str(z)+' done')

#Ref: https://www.micro-shop.zeiss.com/en/us/shop/objectives/420493-9900-000/Objective-EC-Plan-Neofluar-100x-1.3-Oil-Pol-M27#/
#Lens: Plan-Neofluar 100x/1.3
#NA=1.3
#Refractive index on immersion oil = 1.51





