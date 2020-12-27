# In-silico-widefield-microscopy

This open-source software package allows *in-silico* microscopy of molecular simulations such as molecular
dynamics, monte carlo, etc. The current version can perform widefield microscopy 
and approximate forms of confocal, light-sheet, two-photon, and super-resolution microscopy (Cite). 
The programs are developed for GROMACS (.gro) files in mind. Please look at the tutorials for examples.

If you are using these set of programs for your research publications, don't forget to cite (Cite). This software 
package is free to use as long as you comply to the license (GPLv3).


## Manual

### 1. gen_psf.py 

Generates a point spread function file where each line represents the <img src="https://render.githubusercontent.com/render/math?math=l^'">,  <img src="https://render.githubusercontent.com/render/math?math=m^'">,  <img src="https://render.githubusercontent.com/render/math?math=n^'">, and PSF intensity.
The PSF intensity is only printed for <img src="https://render.githubusercontent.com/render/math?math=m^' \leq l^'">.

#### psf_gandy or psf_gandy_sep

**Usage:**

```python
import gen_psf
gen_psf.psf_gandy(beta, dl, dm, dn, Pl, Pm, Pn, fs, outname)
gen_psf.psf_gandy(beta, dl, dm, dn, Pl, Pm, Pn, fs, outname, otype, nidx)
```


This PSF is based on  R. O. Gandy, **1954**, Proc. Phys. Soc. B, 67, 825-831. psf_gandy_sep calculates PSF for a fixed  <img src="https://render.githubusercontent.com/render/math?math=n^'"> coordinate, whereas
psf_gandy calculates PSF for all  <img src="https://render.githubusercontent.com/render/math?math=n^'"> coordinates.

**Arguments:**

- beta = maximum half angle as seen from immersion oil <img src="https://render.githubusercontent.com/render/math?math=\beta=\sin^{-1}(NA/\mu)">, where NA is numerical aperture and  <img src="https://render.githubusercontent.com/render/math?math=\mu"> is the refractive index of of the immersion oil. 
- dl, dm, and dn = are the grid spacing for which PSF is calculated. ( <img src="https://render.githubusercontent.com/render/math?math=\Delta l^', \Delta m^', \Delta n^'"> )
- Pl, Pm, and Pn = are the the dimension of the box for which PSF is calculated (<img src="https://render.githubusercontent.com/render/math?math=P_{l^'}, P_{m^'}, P_{n^'}">). PSF is calculated for (-Pl/2 to Pl/2, -Pm/2 to Pm/2, -Pn/2 to Pn/2). However, because PSF 
  by Gandy is radially and axially symmetrical, it is calculated for (0 to Pl/2, 0 to Pl/2, 0 to Pn/2) and Pl = Pm.
- fs = full-width-at-half-maximum (FWHM) scaling factor; scaling factor of wavenumber (<img src="https://render.githubusercontent.com/render/math?math=k=2\pi/\lambda">); scaling factor of "gro" coordinates. (Cite) 
- outname = the output file name
- otype = open the file "outname" to write or append ('w' or 'a')
- nidx = a counter for (<img src="https://render.githubusercontent.com/render/math?math=n^'">) coordinate. This is an integer int(n/dn).
<img src="https://render.githubusercontent.com/render/math?math=PSF(r,n^')=PSF(l^',m^',n^')=I_0 \left\vert \frac{3}{2(1-\cos^{3/2}\beta)} \int_0^\beta e^{ik^'n^'\cos\theta}J_0(k^'r\sin\theta)\sin\theta\cos^{1/2}\theta d\theta \right\vert^2">

where,<img src="https://render.githubusercontent.com/render/math?math=k^'=\frac{2\pi f_s}{\lambda}">, and <img src="https://render.githubusercontent.com/render/math?math=I_0 = 1">.   

**Output:**

Creates PSF data file with name **outname**

### 2. gen_mono or gen_mono_c

**Usage:**
```bash
gen_mono -f structure.gro -p parameter.dat -o imageheader
```
This evaluates the convolution between PSF and particle number density (ρ), which is evaluated as, 

<img src="https://render.githubusercontent.com/render/math?math=I(l^',m^')=\sum_{j=1}^N PSF(l^'-l_j,m^'-m_j,n_O-n_j)">

The output file consisted of intensity values for <img src="https://render.githubusercontent.com/render/math?math=(l^',m^')"> coordinates. While evaluating <img src="https://render.githubusercontent.com/render/math?math=I">, periodic boundary condition was applied. 
To compare different simulation times, a frame is created which is larger than the molecular simulation box size (over the simulation time). 
The resultant intensity <img src="https://render.githubusercontent.com/render/math?math=I"> box size is scaled with respect to the frame and placed centered in it. The intensity <img src="https://render.githubusercontent.com/render/math?math=I"> was between 0 and 1 (both included). Intensity of -1 was used in 
to represent absence of molecular simulation system and the frame (See Cite for more details).

**Arguments:**

- f: After using "-f" option enter the file name of the GROMACS structure file. **It should contain only one time step**.  
- p: After using "-p" option enter the parameter file name. 
- o: After using "-o" option enter the starting strings of the output filename.

**parameter.dat file:**

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

**File requirements:**

The C-code will search for the PSF file named, **psfheader**\_lam**lam[i]**\_fs**fs**.dat, where **psfheader**, **lam[i]**, and **fs** refer to the values in the parameter file.

For example if parameter.dat contains,
```Note
fs = 20
lam1 = 200
lam2 = 300
psf_header = PSF_gandy
```
it will look for the files ```"PSF_gandy_lam200_fs20.dat"``` and ```"PSF_gandy_lam300_fs20.dat"```.

**structure.gro file:**

The structure file of the specimen should be in GROMACS structure format. 

**Output:**

monochrome image data files are saved with names, **imageheader**_lam**lam[i]**_fs**fs**.dat, where **psfheader**, **lam[i]**, and **fs** refer to the values in the parameter file.

For example if the following command is used,
```Note
gen_mono -f structure.gro -p parameter.dat -o ABC
```

and parameter.dat contains,
```Note
fs = 20
lam1 = 200
lam2 = 300
```

then the output files are ```"ABC_lam200_fs20.dat"```, ```"ABC_lam300_fs20.dat"```.


### 3. render_mono.py

**Usage:**
```bash
python render_mono.py -f imageheader -p param.dat -t timestep
```

This generates an *in-silico* monochrome microscopy image (PNG with 1200 dpi). The frame is generated as white color. 
The monochrome intensity is generated as with grey colormap.

**Arguments:**

- f: After using "-f" option enter the starting characters of the image data files created by gen_mono.
- p: After using "-p" option enter the image parameters filename. 
- t: After using "-t" option enter the timestep associated with the image.

**param.dat file:**

- T: (int). Number of timesteps to generate a time-averaged image. Use the value of 1 to avoid time-averaging. 
- fs: (int). full-width-at-half-maximum (FWHM) scaling factor; scaling factor of wavenumber <img src="https://render.githubusercontent.com/render/math?math=\left( k=2\pi/\lambda\right)">; scaling factor of "gro" coordinates. (<img src="https://render.githubusercontent.com/render/math?math=f_s"> in Cite).
- lam1: (int). The wavelength of light emitted by fluorophore of type 1. Similar syntax for lam2, ..., lam10. (<img src="https://render.githubusercontent.com/render/math?math=\lambda"> in Cite).
- lam1_I0: (float). The maxmium intensity of of light emitted by fluorophore of type 1. Similar syntax for lam2_I0, ..., lam10_I0. (<img src="https://render.githubusercontent.com/render/math?math=I_0"> in Cite).
- size: (float). The largest dimension of molecular simulation box in m direction (<img src="https://render.githubusercontent.com/render/math?math=B_m^*"> in Cite).
- scale: (float). Length of the scale bar to be drawn.

**File requirements:**

When the argument after "-t" is greater than or equal to zero, it searches for the files **imageheaderTimestep**\_lam**lam[i]**\_fs**fs**.dat. 

For example if the following command is used,
```Note
python render_mono.py -f ABC -p param.dat -t 10
```

and parameter.dat contains,
```Note
fs = 20
lam1 = 200
lam2 = 300
lam1_I0 = 0.1
lam2_I0 = 0.05
```
then image data files ```"ABC10_lam200_fs20.dat"```, ```"ABC10_lam300_fs20.dat"``` will be used to render the monochrome image.

When the argument after "-t" is less than zero, it searches for the files **imageheader**\_lam**lam[i]**\_fs**fs**.dat. 

For example if the following command is used,
```Note
python render_mono.py -f ABC -p param.dat -t -1
```

and parameter.dat contains,
```Note
fs = 20
lam1 = 200
lam2 = 300
lam1_I0 = 0.1
lam2_I0 = 0.05
```

then image data files ```"ABC_lam200_fs20.dat"```, ```"ABC_lam300_fs20.dat"``` will be used to render the monochrome image.

**Output:**

When the argument after "-t" is greater than or equal to zero, files mono_**imageheaderTimestep**\_lam**lam[i]**\_fs**fs**\_I**lam[i]\_I0**.png will be created. For the above example, ```"mono_ABC10_lam200_fs20_I0.1.png"``` and ```"mono_ABC10_lam300_fs20_I0.05.png"```.

When the argument after "-t" is less than zero, files mono_**imageheader**\_lam**lam[i]**\_fs**fs**\_I**lam[i]\_I0**.png will be created. For the above example, ```"mono_ABC_lam200_fs20_I0.1.png"``` and ```"mono_ABC_lam300_fs20_I0.05.png"```.

### 4. mono2color.py

**Usage:**  
```bash
python mono2color.py -f imageheader -p param.dat -t timestep
```

This generates a colored *in-silico* microsocpy image (PNG with 1200 dpi).

**Arguments:**

- f: Same as **render_mono.py**
- p: Same as **render_mono.py**
- t: Same as **render_mono.py**

**File requirements:**

Same as **render_mono.py**

**param.dat file:**
- T: (int). Same as **render_mono.py**
- fs: (int). Same as **render_mono.py**
- lam1: (int). Same as **render_mono.py**
- lam1_I0: (float). Same as **render_mono.py**
- size: (float). Same as **render_mono.py**
- scale: (float). Same as **render_mono.py**
- lam1_hue: (int). The artificial hue (in degrees) assigned to fluorophore of type 1. Similar syntax for lam2_hue, ..., lam10_hue. 

**Output:**

When the argument after "-t" is greater than or equal to zero, files **imageheaderTimestep**\_fs**fs**\_T**T**\_I_**lam_I0s**.png will be created, where **lam_I0s** is a string of all lam[i]\_I0 separated by \_.

For example if the following command is used,
```Note
python mono2color.py -f ABC -p param.dat -t 10
```

and parameter.dat contains,
```Note
fs = 20
T = 1
lam1 = 200
lam2 = 300
lam1_I0 = 0.1
lam2_I0 = 0.05
```
then image file ```"ABC10_fs20_T1_I_0.1_0.05.png"``` will be created.


When the argument after "-t" is less than zero, files **imageheader**\_fs**fs**\_I_**lam_I0s**.png will be created. 

For example if the following command is used,
```Note
python mono2color.py -f ABC -p param.dat -t -1
```

and parameter.dat contains,
```Note
fs = 20
T = 1
lam1 = 200
lam2 = 300
lam1_I0 = 0.1
lam2_I0 = 0.05
```
then image file ```"ABC_fs20_T1_I_0.1_0.05.png"``` will be created.

### 5. create_vid.py

**Usage:**  
```bash
python create_vid.py -f imageheader -p param.dat -tmax maxtime -tdiff delta_time
```
This creates a .AVI video from multiple microscopy images (monochrome or colored).

**Arguments:**

- f: (str). After using "-f" option enter the starting characters of the image files to combine and form videos.
- p: (str). After using "-p" option enter the image parameters filename. 
- t0: (int). After using "-t0" enter the first timestep to consider for creating the video.
- tmax: (int). After using "-tmax" enter the maximum timestep to consider for creating the video.
- tdiff: (int). After using the "-tdiff" enter the increaments of timesteps to use for creating the video.
 
**param.dat file:**
- T: (int). Same as **render_mono.py**
- fs: (int). Same as **render_mono.py**
- lam1: (int). Same as **render_mono.py**
- lam1_I0: (float). Same as **render_mono.py**

**File requirements:**

The code searches for timesteps 0, **delta_time**, 2**delta_time**, ..., N**delta_time** (less than or equal to **maxtime**). The specific files searched are **imageheaderTimestep**\_fs**fs**\_T**T**\_I_**lam_I0s**.dat, where **imageheader**, **timestep** is determined from the arguments, **fs**, **T**, and **lamI0s** are determined from the param.dat file (see below). **lamI0s** is a string **lam1_I0**\_**lam2_I0**_ ... **lamN_I0**. 

For example if param.dat contains,

```Note
T = 1
fs = 20
lam1_I0 = 0.1
lam2_I0 = 0.3
lam3_I0 = 0.05
```
and the command,

```bash
python create_vid.py -f ABC -p param.dat -tmax 1000 -tdiff 10
```
then **lamI0s** is "0.1_0.3_0.05", and it will look for files ```"ABC0_fs20_T1_I_0.1_0.3_0.05.png"```, ```"ABC10_fs20_T1_I_0.1_0.3_0.05.png"```, ..., ```"ABC1000_fs20_T1_I_0.1_0.3_0.05.png"```. 

**Output:**

Creates a video file **imageheader**\_fs**fs**\_T**T**\_I_**lam_I0s**.avi. For the above example, the file ```"ABC_fs20_T1_I_0.1_0.3_0.05.avi"``` is created.


### Same param.dat file can be used for render_mono.py, mono2color.py, and create_vid.py

## Contribute

More functionalities can be added in each part of the code. 

### 1. gen_psf.py

- To add new PSF, please add it in gen_psf.py
- Add GPU accelerations, improve the algorithm of PSF calculation (It is currently slow).

### 2. gen_mono

- Adding a generic (l,m,n) coordinates rather than (x,y,z), (y,z,x) and (z,x,y).
- Handle molecular simulation format other than GROMACS.
- Handle non-molecular simulations
- Handle continuum simulations.

### 3. render_mono.py

### 4. mono2color.py

### 5. create_vid.py

### 6. Additional functionalities

- Read Cite.
