# In-silico-widefield-microscopy

This open-source software package allows *in-silico* microscopy of molecular simulations such as molecular
dynamics, monte carlo, etc. The current version can perform widefield microscopy 
and approximate forms of confocal, light-sheet, two-photon, and super-resolution microscopy (Cite). 
The programs are developed for GROMACS (.gro) files in mind. Please look at the tutorials for examples.

If you are using these set of programs for your research publications, don't forget to cite (Cite). This software 
package is free to use as long as you comply to the license (GPLv3).

## Installation

```bash
python setup.py install --user
```

## Quickstart

### 1. Generate Point Spread Function

Generates a point spread function file where each line represents the <img src="https://render.githubusercontent.com/render/math?math=l^'">,  <img src="https://render.githubusercontent.com/render/math?math=m^'">,  <img src="https://render.githubusercontent.com/render/math?math=n^'">, and PSF intensity <img src="https://render.githubusercontent.com/render/math?math=PSF(l^',m^',n^')">. The PSF intensity is only printed for <img src="https://render.githubusercontent.com/render/math?math=m^' \leq l^'">.

Currently the PSF is only evaluated based on R. O. Gandy, **1954**, Proc. Phys. Soc. B, 67, 825-831. 

**Usage: parallel processing**

```bash
siliscopy gen_psf --method gandy \
                 --paramfile [parameter file] \
                 --calc all \
                 --output [output file header] \
                 --multiprocess 
```
**Note1:** Details of parameters that should be included in the parameters file can be found at the end of this Readme file.

**Note2:** This can be run serially but is very slow. Serial code is not extensively tested.

**Note3:** Creates file `output file header`\_lam`lam[i]`_fs`fs`.dat, where `lam[i]` and `fs` are read from `parameter file`.

### 2. Generate Monochrome Image Intensity 

Calculate the monochrome image intensity using the convolution between PSF and particle number density (ρ), 

<img src="https://render.githubusercontent.com/render/math?math=I(l^',m^')=\sum_{j=1}^N PSF(l^'-l_j,m^'-m_j,n_O-n_j)">

The output file consisted of intensity values for <img src="https://render.githubusercontent.com/render/math?math=(l^',m^')"> coordinates. 
While evaluating <img src="https://render.githubusercontent.com/render/math?math=I">, periodic boundary condition was applied.
The intensity <img src="https://render.githubusercontent.com/render/math?math=I"> is between 0 and 1 (both included). Intensity of -1 was used in 
to represent the frame where molecular simulation system is absence (See Cite for more details).

**Single file**

```bash
siliscopy gen_mono --file [GRO file with 1 timestep] \
                   --paramfile [Parameter file] \
                   --output [Output file header]
```

**Multiple file with parallel processing**
```bash
siliscopy gen_mono --data [datafile] \
                   --multiprocess
```
To calculate it serially, remove '--multiprocess'.

Template of data file:
Store all filenames, parameter file (can be the same), output file header separated by ','. For each new monochrome image intensity evaluation, write the arguments in a new line. For example,

```data
file1.gro,param1.dat,out1
file2.gro,param2.dat,out2
```

**Note1:** Index each output file with an index, which can be the timestep of the simulation.
**Note2:** The file names [output file header]_lam[lam[i]]_fs[fs].dat will be created. The value [fs] and [lam[i]] are read from [parameter file].
**Note3:** This runs a C-binary code to make the calculation faster.  
**Note4:** Use parallel processing when the number of files is high. Parallelism is not implemented for files.

### 3. plot images

#### Plot Monochrome Image

**Show specific file**
```bash
siliscopy plot --file [filename header] \
               --paramfile [Parameter file] \
               --method mono \
               --timestep [timestep] \
               --calc show
```

**Save a specific JPEG file**
```bash
siliscopy plot --file [filename header] \
               --paramfile [parameter file] \
               --method mono \
               --timestep [timestep] \
               --calc specific \
               --output [output file header] 
```
**Note1:** If outname is not provided its taken to be the same as [filename header]
**Note2:** The file names [output file header]\[timestep]_lam[lam[i]]_fs[fs]_T[T]_I[lam\_I0[i]].jpeg will be created, wher [lam[i]], [lam\_I0[i]], [fs], [T] are read from [parameter file].
**Note3:** timestep can be negative. This would imply the image does not have a timestep. So the output filename is [output file header]_lam[lam[i]]_fs[fs]_T[T]_I[lam\_I0[i]].jpeg
**Note4:** [filename header] is [output header] used in by gen\_mono.
**Note5:** The calculation method 'specific' and 'spec' yeild the same results. 

**Save multiple JPEG files with parallel processing**
```bash
siliscopy plot --file [filename header] \
               --paramfile [parameter file] \
               --method mono \
               --calc all \
               --output [output file header]
               --multiprocess
```
To calculate it serially remove '--multiprocess'
**Note6:** Use parallel processing when the number of images is high. Parallelism is not implemented for individual images.
**Note7:** The method name 'gray', 'grey' and mono can be used interchangably.


#### Plot Color Image

The commands remain the same as monochrome image. Replace the method 'mono' with 'color'.
 
**Note1:** If outname is not provided its taken to be the same as [filename header]
**Note2:** The file names [output file header]\[timestep]_fs[fs]_T[T]_I[lam\_I0s].jpeg will be created, where [fs], [T] are read from [parameter file]. [lam\_I0s] is a all [lam\_I0[i]] in [parameter file] appended together with a '_' separator. Example, '_[lam\_I0[0]]_[lam\_I0[1]]'.
**Note3:** timestep can be negative. This would imply the image does not have a timestep. So the output filename is [output file header]_fs[fs]_T[T]_I[lam\_I0s].jpeg
**Note4:** [filename header] is [output header] used in by gen\_mono.
**Note5:** The method name 'col' and 'color' can be used interchangably.

### 5. Generate Video

**Monochrome video**
```bash
siliscopy video
                --file [filename header] \
                --paramfile [parameter file] \
                --method mono \
```

**Color video**

Use the method 'color' instead of 'mono'

**Custom video from a list of image filenames**
```bash
siliscopy video --method data \
                --data [data file] \
                --paramfile [parameter file] \
                --output [output filename with extension] 
```
