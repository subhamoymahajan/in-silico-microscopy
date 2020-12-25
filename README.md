# In-silico-widefield-microscopy

This open-source software package allows *in-silico* microscopy of molecular simulations such as molecular
dynamics simulations, monte carlo simulations, etc. The current version can perform widefield microscopy 
and approximate forms of confocal, light-sheet, two-photon, and super-resolution microscopy (Cite). 
The programs are developed for GROMACS (.gro) files in mind. Please look at the tutorials for examples.

If you are using these set of programs for your research publications, don't forget to cite (Cite). This software 
package is free to use as long as you comply to the license (GPLv3).


## Manual

### 1. gen_psf.py 

Generates a point spread function file where each line represents the <img src="https://render.githubusercontent.com/render/math?math=l^'">,  <img src="https://render.githubusercontent.com/render/math?math=m^'">,  <img src="https://render.githubusercontent.com/render/math?math=n^'">, and PSF intensity.
The PSF intensity is only printed for <img src="https://render.githubusercontent.com/render/math?math=m^' \leq l^'">.

#### psf_gandy or psf_gandy_sep
*Usage:* 

```python
import gen_psf
gen_psf.psf_gandy(beta, dl, dm, dn, Pl, Pm, Pn, fs, outname)
gen_psf.psf_gandy(beta, dl, dm, dn, Pl, Pm, Pn, fs, outname, otype, nidx)
```


This PSF is based on  R. O. Gandy, **1954**, Proc. Phys. Soc. B, 67, 825-831. psf_gandy_sep calculates PSF for a fixed  <img src="https://render.githubusercontent.com/render/math?math=n^'"> coordinate, whereas
psf_gandy calculates PSF for all  <img src="https://render.githubusercontent.com/render/math?math=n^'"> coordinates.

**arguments**

- beta = maximum half angle as seen from immersion oil <img src="https://render.githubusercontent.com/render/math?math=\beta=\sin^{-1}(NA/\mu)">, where NA is numerical aperture and  <img src="https://render.githubusercontent.com/render/math?math=\mu"> is the refractive index of of the immersion oil. 
- dl, dm, and dn = are the grid spacing over which PSF is calculated. ( <img src="https://render.githubusercontent.com/render/math?math=\Delta l^', \Delta m^', \Delta n^'"> )
- Pl, Pm, and Pn = are the the dimension of the box over which PSF is calculated (<img src="https://render.githubusercontent.com/render/math?math=P_{l^'}, P_{m^'}, P_{n^'}">). PSF is calculated for (-Pl/2 to Pl/2, -Pm/2 to Pm/2, -Pn/2 to Pn/2). However, because PSF 
  by Gandy is radially and axially symmetrical, it is calculated for (0 to Pl/2, 0 to Pl/2, 0 to Pn/2) and Pl = Pm.
- fs = full-width-at-half-maximum (FWHM) scaling factor; scaling factor of wavenumber (<img src="https://render.githubusercontent.com/render/math?math=k=2\pi/\lambda">); scaling factor of "gro" coordinates. (Cite) 
- outname = the output file name
- otype = open the file "outname" to write or append ('w' or 'a')
- nidx = a counter for (<img src="https://render.githubusercontent.com/render/math?math=n^'">) coordinate. This is an integer int(n/dn).
<img src="https://render.githubusercontent.com/render/math?math=PSF(r,n^')=PSF(l^',m^',n^')=I_0 \left\vert \frac{3}{2(1-\cos^{3/2}\beta)} \int_0^\beta e^{ik^'n^'\cos\theta}J_0(k^'r\sin\theta)\sin\theta\cos^{1/2}\theta d\theta \right\vert^2">

where,<img src="https://render.githubusercontent.com/render/math?math=k^'=\frac{2\pi f_s}{\lambda}">, and <img src="https://render.githubusercontent.com/render/math?math=I_0 = 1">.   

**output file**

Creates PSF data file with name **outname**

### 2. gen_mono or gen_mono_c

Usage:
```bash
gen_mono -f structure.gro -p parameter.dat -o imageheader
```
This evaluates the convolution between PSF and particle number density (ρ), which is evaluated as, 

<img src="https://render.githubusercontent.com/render/math?math=I(l^',m^')=\sum_{j=1}^N PSF(l^'-l_j,m^'-m_j,n_O-n_j)">

The output file consisted of intensity values for <img src="https://render.githubusercontent.com/render/math?math=(l^',m^')"> coordinates. While evaluating <img src="https://render.githubusercontent.com/render/math?math=I">, periodic boundary condition was applied. 
To compare different simulation times, a frame is created which is larger than the molecular simulation box size (over the simulation time). 
The resultant intensity <img src="https://render.githubusercontent.com/render/math?math=I"> box size is scaled with respect to the frame and placed centered in it. The intensity <img src="https://render.githubusercontent.com/render/math?math=I"> was between 0 and 1 (both included). Intensity of -1 was used in 
to represent absence of molecular simulation system (and the frame).

**arguments**

- f:
- p:
- o:

**parameter file**

- f: (int). full-width-at-half-maximum (FWHM) scaling factor; scaling factor of wavenumber <img src="https://render.githubusercontent.com/render/math?math=\left( k=2\pi/\lambda\right)">; scaling factor of "gro" coordinates. (<img src="https://render.githubusercontent.com/render/math?math=f_s"> in Cite).
- maxlen: (float float float). Largest dimension of molecular simulation box in x, y, and z directions (<img src="https://render.githubusercontent.com/render/math?math=B_l^*, B_m^*, B_n^*"> in Cite).
- focus_cor: (float). The n-coordinate (in gro file) to which the *in-silico* microscope is focused (<img src="https://render.githubusercontent.com/render/math?math=n_O"> in Cite). 
- opt_axis: (int). Optical axis direction 0 for x, 1 for y, and 2 for z (<img src="https://render.githubusercontent.com/render/math?math=n"> in Cite).
- lam1: (int). The wavelength of light emitted by fluorophore of type 1. Similar syntax for lam2, ..., lam10. (<img src="https://render.githubusercontent.com/render/math?math=\lambda"> in Cite).
- lam_names1: (string separated by spaces). Bead or atom names of fluorophore of type 1. Similar syntax for lam_names2, ..., lam_names10. Upto 100 atom names are supported for each fluorophore type.
- NA: (float). Numerical aperture. 
- dx: (float float float). The grid spacing for the PSF. (<img src="https://render.githubusercontent.com/render/math?math=\Delta l^', \Delta m^', \Delta n^'">,  in Cite)
- Lpsf: (float float float). The PSF box dimension. (<img src="https://render.githubusercontent.com/render/math?math=P_{l^'}, P_{m^'}, P_{n^'}">)
- psfheader: (string). Starting characters with which PSF was saved.
- pbc: (string). Directions in which PBC was applied. None, x, y, z, xy, yz, xz, or xyz should be used.

The C-code will search for the PSF file named, **psfheader**_lam**lam[i]**_fs**fs**.dat, where **psfheader**, **lam[i]**, and **fs** refer to the values in the parameter file.

**structure.gro**

The structure file of the specimen should be in GROMACS structure format. 

**output file**

monochrome image data files are saved with names, **imageheader**_lam**lam[i]**_fs**fs**.dat, where **psfheader**, **lam[i]**, and **fs** refer to the values in the parameter file.

### 3. render_mono.py

Usage:
```bash
python render_mono.py -f imageheader -p param.dat -t timestep
```

This generates an *in-silico* monochrome microscopy image (PNG with 1200 dpi). The frame is generated as white color. 
The monochrome intensity is generated as with grey colormap.

**arguments**

- f:
- p:
- t:

**param.dat**

- T: (int). Number of timesteps to generate a time-averaged image. Use the value of 1 to avoid time-averaging. 
- fs: (int). full-width-at-half-maximum (FWHM) scaling factor; scaling factor of wavenumber <img src="https://render.githubusercontent.com/render/math?math=\left( k=2\pi/\lambda\right)">; scaling factor of "gro" coordinates. (<img src="https://render.githubusercontent.com/render/math?math=f_s"> in Cite).
- lam1: (int). The wavelength of light emitted by fluorophore of type 1. Similar syntax for lam2, ..., lam10. (<img src="https://render.githubusercontent.com/render/math?math=\lambda"> in Cite).
- lam1_I0: (float). The maxmium intensity of of light emitted by fluorophore of type 1. Similar syntax for lam2_I0, ..., lam10_I0. (<img src="https://render.githubusercontent.com/render/math?math=I_0"> in Cite).
- size: (float). The largest dimension of molecular simulation box in m direction (<img src="https://render.githubusercontent.com/render/math?math=B_m^*"> in Cite).
- scale: (float). Length of the scale bar to be drawn.



**output file**

### 4. mono2color.py

Usage:  
```bash
python mono2color.py -f imageheader -p param.dat -t timestep
```

This generates a colored *in-silico* microsocpy image (PNG with 1200 dpi).

**arguments**

- f:
- p:
- t:

**param.dat**
- T: (int). Same as **render_mono.py**
- fs: (int). Same as **render_mono.py**
- lam1: (int). Same as **render_mono.py**
- lam1_I0: (float). Same as **render_mono.py**
- size: (float). Same as **render_mono.py**
- scale: (float). Same as **render_mono.py**
- lam1_hue: (int). The artificial hue (in degrees) assigned to fluorophore of type 1. Similar syntax for lam2_hue, ..., lam10_hue. 

**output file**

### 5. create_vid.py

Usage:  
```bash
python create_vid.py -f imageheader -p param.dat -tmax maxtime -tdiff delta_time
```
This creates a .AVI video from multiple microscopy images (monochrome or colored).

**arguments**

- f:
- p:
- tmax:
- tdiff:

**param.dat**
- T: (int). Same as **render_mono.py**
- fs: (int). Same as **render_mono.py**
- lam1: (int). Same as **render_mono.py**
- lam1_I0: (float). Same as **render_mono.py**


**output file**
