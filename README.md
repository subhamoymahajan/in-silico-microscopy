# *In-silico* microscopy

This open-source software package allows *in-silico* microscopy of molecular simulations such as molecular
dynamics, monte carlo, etc. The current version can perform widefield microscopy 
and approximate forms of confocal, light-sheet, two-photon, and super-resolution microscopy [1]. 
The programs are developed for GROMACS (.gro) files in mind. Please look at the tutorials for examples.

If you are using these set of programs for your research publications, don't forget to cite [1]. This software 
package is free to use as long as you comply to the license (GPLv3).

## Published Work

This toolbox is based on published work.

[1] Mahajan S., and Tang T., Meeting experiments at the diffraction barrier: an in-silico widefield fluorescence microscopy, bioRxiv 2021.03.02.433395;  DOI: [10.1101/2021.03.02.433395](https://doi.org/10.1101/2021.03.02.433395) (Currently in peer-review)

## Installation

```bash
python setup.py install --user
```

## Quickstart

### 1. Generate Point Spread Function (PSF)

Generates a PSF file where each line represents the <img src="https://render.githubusercontent.com/render/math?math=l^'">,  <img src="https://render.githubusercontent.com/render/math?math=m^'">,  <img src="https://render.githubusercontent.com/render/math?math=n^'">, and PSF intensity <img src="https://render.githubusercontent.com/render/math?math=PSF(l^',m^',n^')">. The PSF intensity is only printed for <img src="https://render.githubusercontent.com/render/math?math=m^' \leq l^'">.


```bash
siliscopy gen_psf --method [Gandy/GL1991/Mod_Gandy] \
                 --paramfile [parameter file] \
                 --calc all \
                 --output [output file header] \
                 --multiprocess 
```
> **_Note1:_** Details of parameters that should be included in the parameters file can be found at the end of this Readme file. 
>   
> **_Note2:_** This can be run serially but is very slow. Serial code is not extensively tested.  
>  
> **_Note3:_** Creates file `output file header`\_lam`lam[i]`\_fs`fs`.dat, where `lam[i]` and `fs` are read from `parameter file`.  
>  

Current version 1.2.2 can generate three PSFs. Depth-invariant PSF (`psf_type = 0`) Gandy (`Gandy`), and depth-variant (`psf_type = 1`) Gibson and Lanni (`GL1991`), and Modified Gandy (`Mod_Gandy`). 

[2] Gandy PSF: R. O. Gandy, Out-of-Focus Diffraction Patterns for Microscope Objectives, *Proc. Phys. Soc. B* **1954**, 67, 825-831.
[3] Gibson and Lanni PSF: S. F. Gibson, and F. Lanni, Experimental test of an analytical model of aberration in an oil-immersion objective lens used in three-dimensional light microscopy, *J. Opt. Soc. Am. A* **1991**, 8, 1601-1613.
[4] Gibson and Lanni PSF: S. F. Gibson, and F. Lanni, Experimental test of an analytical model of aberration in an oil-immersion objective lens used in three-dimensional light microscopy, *J. Opt. Soc. Am. A* **1992**, 9, 154-166. (Same as [3] but better images)
[5] Modified Gandy PSF: Upcoming publication.

> **_Note4:_** The functional form of Gibson and Lanni PSF is changed to be compatible with `siliscopy` and will be reported in the upcoming journal.
>  


### 2. Generate Monochrome Image Intensity 

Calculate the monochrome image intensity using the convolution between PSF and particle number density (ρ), 

<img src="https://render.githubusercontent.com/render/math?math=I(l^',m^')=\sum_{j=1}^N PSF(l^'-l_j,m^'-m_j,n_O-n_j)">

The output file consisted of intensity values for <img src="https://render.githubusercontent.com/render/math?math=(l^',m^')"> coordinates. 
While evaluating <img src="https://render.githubusercontent.com/render/math?math=I">, periodic boundary condition was applied.
The intensity <img src="https://render.githubusercontent.com/render/math?math=I"> is between 0 and 1 (both included). Intensity of -1 was used 
in to represent the white frame where molecular simulation system is absent (See Ref [1] for more details). 
The image can be generated for a specific object focal plane (`slice`) or multiple object focal planes (3D image; `volume`) using the method option.
The default method is `slice`. 


Single file:

```bash
siliscopy gen_mono --file [GRO file] \
                   --paramfile [parameter file] \
                   --psf [PSF file header]
                   --output [output file header]
                   --method [slice/volume]
```

Multiple file with parallel processing:
```bash
siliscopy gen_mono --data [datafile] \
                   --multiprocess
```
To calculate it serially, remove `--multiprocess`.

Template of data file:

Store all filenames, parameter file (can be the same), PSF header file name (can be the same), output file header separated by ','. For each new monochrome image intensity evaluation, write the arguments in a new line. For example,

```data
file1.gro,param1.dat,psf.dat,out1,method
file2.gro,param2.dat,psf.dat,out2,method
```

> **_Note1:_** Index each output file with an index, which can be the timestep of the simulation.  
>  
> **_Note2:_** The file names `output file header`\_lam`lam[i]`\_fs`fs`.dat will be created. The value `fs` and `lam[i]` are read from the parameter file.  
>  
> **_Note3:_** This runs a C-binary code for faster calculation.    
>  
> **_Note4:_** Use parallel processing when the number of files is high. Parallelism is not implemented for individual files.  
>  
> **_Note5:_** `GRO file` must contain only one timestep.  
>  
> **_Note6:_** `PSF file header` is the same as `--output` in `gen_psf`.
>

### 3. Plot Images

Generate *In-silico* microscopy images. If saved to a file, all images are saved with `JPEG` file format.

#### Plot monochrome image

Show specific image
```bash
siliscopy plot --file [filename header] \
               --paramfile [Parameter file] \
               --method [mono/mono3d/mono3dt/noise_mono/noise_mono3d/noise_mono3dt] \
               --timestep [timestep] \
               --calc show
```

Save a specific image:
```bash
siliscopy plot --file [filename header] \
               --paramfile [parameter file] \
               --method [mono/mono3d/mono3dt/noise_mono/noise_mono3d/noise_mono3dt] \
               --timestep [timestep] \
               --calc specific \
               --output [output file header] 
               --type [jpeg/png/tiff8/tiff16]
```
> **_Note1:_** If `--output` is not provided, then `output file header` its taken to be the same as `filename header`  
>  
> **_Note2:_** The file names `[output file header][timestep]`\_lam`lam[i]`\_fs`fs`\_T`T`\_I`lam_I0[i]`.jpeg will be created, where `lam[i]`, `lam_I0[i]`, `fs`, `T` are read from `parameter file`.  
>  
> **_Note3:_** The timestep can be negative. This would imply the image does not have a timestep. Then, the output filename will be `output file header`\_lam`lam[i]`\_fs`fs`\_T`T`\_I`lam_I0[i]`.jpeg
>  
> **_Note4:_** `filename header` should be the same as `output header` used in by `gen_mono`.  
>   
> **_Note5:_** The calculation method `specific` and `spec` yeild the same results.   

Save multiple images with parallel processing:
```bash
siliscopy plot --file [filename header] \
               --paramfile [parameter file] \
               --method [mono/mono3d/mono3dt/noise_mono/noise_mono3d/noise_mono3dt] \
               --calc all \
               --output [output file header]
               --type [jpeg/png/tiff8/tiff16]
               --multiprocess
```
To calculate it serially remove `--multiprocess`

> **_Note6:_** Use parallel processing when the number of images is high. Parallelism is not implemented for individual images.  
>  

#### Plot Color Image

The commands remain the same as monochrome image. Replace the method `mono` with `color`.
 
> **_Note1:_** If outname is not provided its taken to be the same as `filename header`  
>  
> **_Note2:_** The file names `[output file header][timestep]`\_fs`fs`\_T`T`\_I`lam_I0s`.jpeg will be created, where `fs`, `T` are read from `parameter file`. `lam_I0s` is a all `lam_I0[i]` in `parameter file` appended together with a '\_' separator. Example, '\_[lam\_I0[0]]\_[lam\_I0[1]]'.  
>  
> **_Note3:_** timestep can be negative. This would imply the image does not have a timestep. So the output filename is `output file header`\_fs`fs`\_T`T`\_I`lam_I0s`.jpeg  
>  
> **_Note4:_** `filename header` is `output header` used in by `gen_mono`.  
>  
#### Other Plots

Use the method `lumin` to plot lumination and hue of pixels about a central pixel interactively. 

Use the method `region` to plot a region plot where hues occur in factors of 10 degrees and values are made one. 

These methods will be demonstrated in future tutorial. For now read the documentation.

  
### 5. Generate Video

Generates a movie from a set of `JPEG` images.

Monochrome video:  
```bash
siliscopy video
                --file [filename header] \
                --paramfile [parameter file] \
                --method mono \
```

Color video:  
Use the method 'color' instead of 'mono'

Custom video from a list of image filenames:

```bash
siliscopy video --method data \
                --data [data file] \
                --paramfile [parameter file] \
                --output [output filename with extension] 
```
### 6. Calculate Properties 

Calculate the maximum intensity use the command:

```bash
siliscopy prop --method maxI  \
               --file [filename header] \
               --paramfile [parameter file] \
```

Predict the maximum intensity to generate images:
```bash
siliscopy prop --method predI0  \
               --file [filename header] \
               --paramfile [parameter file] \
               --timestep [iterations]
```
> **_Note1:_** The parameter for timestep is used as maximum number of iterations to predict I0.   
>  

Calculate histogram of image intensities:
```bash
siliscopy prop --method hist  \
               --file [filename header] \
               --paramfile [parameter file] \
               --output [output filename]
```

To create a histogram with normalized count use the options `--calc norm`.

Calculate number of particles, area, and binary images:
```bash
siliscopy prop --method num_area \
               --file [filename header] \
               --paramfile [parameter file] \
               --calc [show/show-test/test] \ 
               --threshold [intensity] \
               --lambdaID [index of fluorophore] \
               --timestep [timestep] \
               --output [output filename header]
```
Use `calc` value of 'show' to show the figure. This option won't save the file shown. 'show-test' 
is used for check the area calculation algorithm is working correctly. If the pixels used to 
calculate particle's area is shown in grey. In this option if this option shows particles in white,
report the issue. Option 'test' will save the image shown in 'show-test'. Any other option will save
generated binary image. 

`threshold` is used to set the threshold intensity for particle detection. If the value is below 1 
the maximum intensity is assumed to be 1 otherwise maximum intensity is assumed to be 255. 

`lambdaID` is used to chose a specific fluorophore. 

`[output filename header].dat` is created which stores area of different particles. Each line contains,
the particle ID (starting from 0) and its area. `[output filename header].png` is the binary image used 
to calculate the area. It is created based on the value of `calc`.

### 7. Parameter File

The parameter file should contain the following parameters (Not all are used in every step). 

* `fs`: (int). Full-Width-at-Half-Maximum (FWHM) scaling factor. Equivalently scales all molecular simulation coordinates or the wave vector.  
* `maxlen`: (float float float). Maximum molecular simulation box dimensions in nm.
* `focus_cor`:  (float). The n-coordinate (in gro file) at whic the *in-silico* microscope is focused in nm.
* `opt_axis`: (int). The direction of the optical axis. 0 => x, 1 => y, 2 => z.
* `lam[i]`: (int). Wavelenght of [i]th fluorophore type in nm. Replace [i] with integers starting from 1. Currently 10 wavelengths are supported. 
* `lam_names[i]`: (str str ... str). Atom names of [i]th fluorophore type. Replace [i] with integers strating from 1. Currently 200 names are supported.
* `dlmn`: (float float float). The voxel dimensions <img src="https://render.githubusercontent.com/render/math?math=\Delta l^', \Delta m^', \Delta n^'"> in nm.  
* `Plmn`: (float float float). The dimensions of the box in nm within which PSF is calculated; <img src="https://render.githubusercontent.com/render/math?math=P_{l^'}, P_{m^'}, P_{n^'}">
* `pbc`: (str). Directions in which periodic boundary condition is active. `None`, `x`, `y`, `z`, `xy`, `yz`, `xz`, or `xyz`.
* `NA`: (float). Numerical aperture.
* `meu`: (float). Refractive index of immersion oil.

* `T`: (int). Number of consecutive timesteps to average the the *in-silico* microscopy image.
* `lam_hue[i]`: (float). Hue in degrees of [i]th fluorophore type. Replace [i] with integers starting from 1.
* `lam_I0_[i]`:  (float). Maximum image intensity of of [i]th fluorophore type. Replace [i] with integers starting from 1.
* `scale`: (float). Length of scale bar in nm. 
* `dpi`: (int). Dots per square inch of the output image. The image resolution is also dependent on `dlmn`. 
* `tbegin`: (int). Index associated with first timestep (will be included).
* `tmax`: (int). Index associated with last timestep (will not be included).
* `tdiff`: (int). Difference in index associated with timestep.
* `fps`: (int). Frames per second of the output video. (default value is 1)
* `vid_ext`: (str). The extension of the output video. (default value is .mov)
* `fourcc`: ('str'). The four character code of the encoder with which video will be created. (default value is 'mp4v') 

Parameters for Gibson-Lanni and Modified-Gandy PSF.

* `meu0`: (float). Refractive index of immersion oil in design condition. 
* `t0`: (float). Thickness of immersion oil in design condition in nanometer. 
* `meug`: (float). Refractive index of coverslip. 
* `meug0`: (float). Refractive index of coverslip in design condition. 
* `tg`: (float). Thickness of immersion oil in design condition in nanometer. 
* `tg0`: (float). Thickness of immersion oil in design condition in nanometer. 
* `meus`: (float). Refractive index of specimen. 
* `tsO`: (float). Location of object focal plane below the coverslip in nanometers. 
* `psf_type`: (int). Default value is 0 is for Gandy PSF (Circular symmetry and depth-invariant). Value of 1 is for Gibson-Lanni and Modified-Gandy (Circular symmetry and depth-variant).  

# Contributors
1. Subhamoy Mahajan, PhD student, Mechanical Engineering, University of Alberta, Edmonton, Canada (will graduate by Sep, 2021)
2. Tian Tang, Professor, Mechanical Engineering, University of Alberta, Edmonton, Canada 

# Funding Sources
This research was funded by Natural Science and Engineering Research Council of Canada through Tian Tang. My research was also funded by several scholarships, Mitacs Globalink Graduate Fellowship (2016-2019), RR Gilpin Memorial Scholarship (2019), Alberta Graduate Excellence Scholarship (2020), Sadler Graduate Scholarship (2020).
