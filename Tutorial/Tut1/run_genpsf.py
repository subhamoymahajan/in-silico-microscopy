import sys
sys.path.insert(1,'../../')
import gen_psf as gpsf
NA=1.3
meu=1.51
lambd=518
dx=0.1
dy=0.1
dz=0.2
Lx=15
Ly=15
Lz=50
fs=800
#gpsf.psf_gandy(NA,meu,lambd,dx,dy,dz,Lx,Ly,Lz,fs,'PSF_gandy_lam'+str(lambd)+'_fs'+str(fs)+'.dat')
lambd=670
#gpsf.psf_gandy(NA,meu,lambd,dx,dy,dz,Lx,Ly,Lz,fs,'PSF_gandy_lam'+str(lambd)+'_fs'+str(fs)+'.dat')
for k in range(83,126):
    gpsf.psf_gandy_sep(NA,meu,lambd,dx,dy,dz,Lx,Ly,Lz,fs,'PSF_gandy_lam'+str(lambd)+'_fs'+str(fs)+'.dat','a',k)
