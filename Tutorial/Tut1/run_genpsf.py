import sys
sys.path.insert(1,'../../')
import gen_psf as gpsf
import numpy as np
NA=1.3
meu=1.51
beta=np.arcsin(NA/meu)
lambd=518
dl=0.1
dm=0.1
dn=0.2
Ll=15
Lm=15
Ln=50
fs=800
gpsf.psf_gandy(beta,lambd,dl,dm,dn,Ll,Lm,Ln,fs,'PSF_gandy_lam'+str(lambd)+'_fs'+str(fs)+'.dat')
lambd=670
gpsf.psf_gandy(beta,lambd,dl,dm,dn,Ll,Lm,Ln,fs,'PSF_gandy_lam'+str(lambd)+'_fs'+str(fs)+'.dat')
