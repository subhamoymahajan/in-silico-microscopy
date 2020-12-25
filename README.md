# In-silico-widefield-microscopy

This open-source software package allows *in-silico* microscopy of molecular simulations such as molecular
dynamics simulations, monte carlo simulations, etc. The current version can perform widefield microscopy 
and approximate forms of confocal, light-sheet, two-photon, and super-resolution microscopy (Cite). 
The programs are developed for GROMACS (.gro) files in mind. Please look at the tutorials for examples.

If you are using these set of programs for your research publications, don't forget to cite (Cite). This software 
package is free to use as long as you comply to the license (GPLv3).

##Manual

###1. gen_psf.py 

Generates a point spread function file where each line represents the l, m, n, and PSF intensity.
The PSF intensity is only printed for $m \leq l$.

#### a. psf_gandy or psf_gandy_sep
Usage: 
import gen_psf
gen_psf.psf_gandy(beta,dl,dm,dn,Ll,Lm,Ln,fs,outname)
gen_psf.psf_gandy(beta,dl,dm,dn,Ll,Lm,Ln,fs,outname,otype,nidx)


This PSF is based on  R. O. Gandy, **1954**, Proc. Phys. Soc. B, 67, 825-831. psf_gandy_sep calculates PSF for a fixed z coordinate, whereas
psf_gandy calculates PSF for all z coordinates.
- beta = maximum half angle as seen from immersion oil; $sin^{-1}(NA/\mu), where NA is numerical aperture
  and μ is the refractive index of of the immersion oil. 
- dl, dm, and dn = are the grid spacing over which PSF is calculated. 
- Ll, Lm, and Ln = are the the dimension of the box over which PSF is calculated. PSF is calculated for (-Ll/2 to Ll/2, -Lm/2 to Lm/2, -Ln/2 to Ln/2). However, because PSF
  by Gandy is radially and axially symmetrical, it is calculated for (0 to Ll/2, 0 to Ll/2, 0 to Ln/2) and Ly = Lx.
- fs = full-width-at-half-maximum (FWHM) scaling factor; scaling factor of wavenumber (k=2pi/lambda); scaling factor of "gro" coordinates. (Cite) 
- outname = the output file name
- otype = open the file "outname" to write or append ('w' or 'a')
- nidx = a counter for n coordinate. This is an integer int(n/dn).

$PSF(r,n')=PSF(l',m',n')=I_0 \left\vert \frac{3}{2(1-\cos^{3/2}\beta)} \int_0^\beta e^{ik'n'\cos\theta}J_0(k'r\sin\theta)\sin\theta\cos^{1/2}\theta d\theta \right\vert^2$

where, $I_0$ is taken as 1.   

**output file**
Creates PSF data file with name [**outname**]

#### b. gen_mono or gen_mono_c
Usage:
gen_mono -f structure.gro -p parameter.dat -o imageheader

This evaluates the convolution between PSF and particle number density (ρ), which is evaluated as, 

$$I(l',m')=\sum_{j=1}^N PSF(l'-l_j,m'-m_j,n_O-n_j)$$

The output file consisted of intensity values for (l',m') coordinates. While evaluating $I$, periodic boundary condition was applied. 
To compare different simulation times, a frame is created which is larger than the molecular simulation box size (over the simulation time). 
The resultant intensity $I$ box size is scaled with respect to the frame and placed centered in it. The intensity $I$ was between 0 and 1 (both included). Intensity of -1 was used in 
to represent absence of molecular simulation system (and the frame).

**parameter file**
- f: (int). full-width-at-half-maximum (FWHM) scaling factor; scaling factor of wavenumber (k=2pi/lambda); scaling factor of "gro" coordinates. ($f_s$ in Cite).
- maxlen: (float float float). Largest dimension of molecular simulation box in x, y, and z directions ($B_l^*$, $B_m^*$, $B_n^*$ in Cite).
- focus_cor: (float). The n-coordinate (in gro file) to which the *in-silico* microscope is focused ($n_O$ in Cite). 
- opt_axis: (int). Optical axis direction 0 for x, 1 for y, and 2 for z ($n$ in Cite).
- lam1: (int). The wavelength of light emitted by fluorophore of type 1. Similar syntax for lam2, ..., lam10. ($\lambda$ in Cite).
- lam_names1: (string separated by spaces). Bead or atom names of fluorophore of type 1. Similar syntax for lam_names2, ..., lam_names10. Upto 100 atom names are supported for each fluorophore type.
- NA: (float). Numerical aperture. 
- dx: (float float float). The grid spacing for the PSF. ($\Delta l'$, $\Delta m'$, $\Delta n'$ in Cite)
- Lpsf: (float float float). The PSF box dimension. 
- psfheader: (string). Starting characters with which PSF was saved.
- pbc: (string). Directions in which PBC was applied. None, x, y, z, xy, yz, xz, or xyz should be used.

The C-code will search for the PSF file named, [**psfheader**]_lam[**lam[i]**]_fs[**fs**].dat, where [psfheader], [fs] refer to the values in the parameter file, and lam[i] refers to value of lam1, ..., lam10 in parameters file.

**structure.gro**
The structure file of the specimen should be in GROMACS structure format. (link)

**output file**
monochrome image data files are saved with names, [**imageheader**]_lam[**lam[i]**]_fs[**fs**].dat.

#### b. render_mono.py

Usage: python render_mono.py -f imageheader -p param.dat -t timestep

This generates an *in-silico* monochrome microscopy image (PNG with 1200 dpi). The frame is generated as white color. 
The monochrome intensity is generated as with grey colormap.

**param.dat**
**output file**

#### c. mono2color.py
Usage:  python mono2color.py -f imageheader -p param.dat -t timestep

This generates a colored *in-silico* microsocpy image (PNG with 1200 dpi).
**param.dat**
**output file**

#### d. create_vid.py
Usage:  python create_vid.py -f imageheader -p param.dat -tmax maxtime -tdiff delta_time

This creates a .AVI video from multiple microscopy images (monochrome or colored).

**param.dat**
**output file**
