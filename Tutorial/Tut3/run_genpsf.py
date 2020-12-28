import sys
sys.path.insert(1,'../../')
import gen_psf as gpsf
import numpy as np
NA=1.3
meu=1.51
beta=np.arcsin(NA/meu)
lambd=518
dx=0.1
dy=0.1
dz=0.2
Lx=15
Ly=15
Lz=50
fs=600
gpsf.psf_gandy(beta,lambd,dx,dy,dz,Lx,Ly,Lz,fs,'PSF_gandy_lam'+str(lambd)+'_fs'+str(fs)+'.dat')
lambd=670
gpsf.psf_gandy(beta,lambd,dx,dy,dz,Lx,Ly,Lz,fs,'PSF_gandy_lam'+str(lambd)+'_fs'+str(fs)+'.dat')
