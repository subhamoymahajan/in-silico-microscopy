# Tutorial 10: Using Siliscopy with Fiji ImageJ

This tutorial will cover how to generate images that are compatible with the [Fiji ImageJ](https://imagej.net/software/fiji/) software.

## Converting Files for ImageJ


### Converting PSF data file to TIFF

Convert PSF to 3D black and white tiff: 
```bash
siliscopy convert --method psf2tiff --file PSF_gandy_lam518_fs530.dat --paramfile parameters.dat --calc uint8 --output PSF_518.tiff
```

To view the PSF in Fiji ImageJ:

1. Open PSF\_518.tiff file
2. Select Plugins > 3D Viewer
3. Make sure the correct `Image` is chosen. 
4. Change the desired color or use default settings and click OK
5. Select Edit > Show bounding box
6. With the cursor find the desired view. Use `View > Center > ` and `View > Fit view to > ` settings if necessary.

![PSF 3D view](https://github.com/subhamoymahajan/in-silico-microscopy/Tutorial/Tut10/PSF_view.png?raw=true)

The PSF generated using psf2tiff can be used for deconvolution in Fiji ImageJ (Shown later).


Convert PSF to 3D color tiff: 
```bash
siliscopy convert --method psf2tiff2 --file PSF_gandy_lam518_fs530.dat --paramfile parameters.dat --calc uint8 --output PSF2_518.tiff
```

The PSF generated using psf2tiff2 is good for visualization and cannot be used for deconvolution in Fiji ImageJ. 
When this PSF is opened it would appear black and white.  To displa color select `Image > Color > Make Composite`. 
Folow the steps above to see the 3D view. Additionally, select `View > Change background color` and change Red, Green and Blue to 255.

![Color PSF 3D view](https://github.com/subhamoymahajan/in-silico-microscopy/Tutorial/Tut10/PSF2_view.png?raw=true)

0-0.2 Changes the color black to fully saturated red. (Changing the background color to white changes it from transparent white to red).

0.2-0.8 changes the color from fully saturated red to fully saturated blue.

0.8-1.0 Chnages the color from fully saturated blue to white.

### Converting a stack of images to a volume TIFF

The spacing for volume TIFF is determined from `dlmn` and `add_n`.

(i) Convert YX image to ZYX.
```bash
siliscopy convert --method nstack2tiff --data nlist1.dat --output IMG_n1.tiff --paramfile parameters.dat 
```

(ii) Convert CYX image to ZCYX.
```bash
siliscopy convert --method nstack2tiff --data nlist2.dat --output IMG_n2.tiff --paramfile parameters.dat
```

(iii) Convert TYX image to TZYX.
```bash
siliscopy convert --method nstack2tiff --data nlist3.dat --output IMG_n3.tiff  --paramfile parameters.dat
```

(iv) Convert TCYX image to TZCYX.
```bash
siliscopy convert --method nstack2tiff --data nlist4.dat --output IMG_n4.tiff  --paramfile parameters.dat
```

### Converting a stack of TIFF images into a t-stacked (time) TIFF

Interval between images are determiend from `fpns` and unit of time is nanoseconds.

(i) Convert YX image to TYX.
```bash
siliscopy convert --method tstack2tiff --data tlist1.dat --output IMG_t1.tiff  --paramfile parameters.dat
```

(ii) Convert CYX image to TCYX.
```bash
siliscopy convert --method tstack2tiff --data tlist2.dat --output IMG_t2.tiff  --paramfile parameters.dat
```

(iii) Convert ZYX image to TZYX.
```bash
siliscopy convert --method tstack2tiff --data tlist3.dat --output IMG_t3.tiff  --paramfile parameters.dat
```

(iv) Convert ZCYX image to TZCYX.
```bash
siliscopy convert --method tstack2tiff --data tlist4.dat --output IMG_t4.tiff  --paramfile parameters.dat
```

Note: It is best to create a 3D TIFF directly using 'mono2dt' or 'color2dt' method.

### Combining multiple greyscale TIFF to a color TIFF

`lams`, `lam_hues`, `mix_type` are read from the parameter file.

(i) Convert YX image to CYX.
```bash
siliscopy convert --method img2color --data clist1.dat --output IMG_c1.tiff --paramfile parameters.dat
```

(ii) Convert ZYX image to ZCYX.
```bash
siliscopy convert --method img2color --data clist2.dat --output IMG_c2.tiff --paramfile parameters.dat
```

(iii) Convert TYX image to TCYX.
```bash
siliscopy convert --method img2color --data clist3.dat --output IMG_c3.tiff --paramfile parameters.dat
```

(iv) Convert TZYX image to TZCYX.
```bash
siliscopy convert --method img2color --data clist4.dat --output IMG_c4.tiff --paramfile parameters.dat
```

## Calculating number of particles

### Using ImageJ (no PBC!)

Number of particles and area can be calculated using the in-built  `Analyze particles` using the following steps.

1. Open an image.
2. Theshold the image and generate a binary image using `Image > Adjust > Threshold` then choose the upper and lower limit and click `Apply`.
3. Uncheck `Calculate threshold for each image` to ensure the same threshold is applied to each image and click `OK`.
4. Select `Analyze > Analyze particle`. Choose minimum and maximum size, range of circularity, and `show` options then click `OK`.
5. For multidimensional images Fiji ImageJ creates a single list. You can choose not to calculate area for all images. Choose `show = oultines` to get identify area of particle with the ouline.

However, Fiji imageJ does not natively apply periodic boundary condition while calculating particles and their area.

To calculate number of paricles and volume,
1. Select `Plugin > 3D > 3D Simple Segmentation` to threshold the image. The plugin has other thresholding methods.
2. Analyze properties such as volume and surface area using `Plugin > 3D > 3D Geometric measure`. The plugin offers many more options to calculate differnet properties.


## Colocalization

Colocalization analysis can be performed using several Fiji ImageJ plugins, such as, Just Another Colocalization Plugin (JACoP).

1. Open Image.
2. Select `Image > Color > Split Channel` to split the channels in the image.
3. Select `Plugins > JACop`.
4. Choose interested images in `Image A` and `Image B`.
5. Tick the interested analysis under `Analysis to perform`.
6. Near the bottom, required panels would turn red, such as `Threshold`.
7. For example, in threshold panel, chooese the threshold intensity for each image and click `Analyze`.

If you want to perform the analysis on a portion of multidimensional images, such as a particular Z slice in a ZYX image, perform the following steps before analysis.

1. Open Image.
2. Select `Image > Duplicate`.
3. Choose the required ranges and click `OK`.

To crop the image,
1. Open Image.
2. Make sure the rectangle tool is selected, and then choose a portion of the image using a mouse.
3. Select `Image > Crop`.


## Deconvolution

Deconvolution can be performed using Fiji ImageJ plugin DeconvolutionLab2 `Plugin > DeconvolutionLab2 > DeconvolutionLab2 Lab`. The Plugin has easy to use intructions.


## Single-Particle Tracking

Single particle tracking can be performed using Fiji ImageJ plugin TrackMate `Plugin > Tracking > TrackMate`. The plugin has easy to use instructions.
