# Tutorial 9: Calculating Properties from Images


## Determining Maximum Intensity

```bash
siliscopy prop --method maxI --file img --paramfile parameters.dat --timestep 100
```

This opens the two image files and reports their maximum intensity. The maximum intensity for 670 and 518 nm light is 12.71 and 4.047 respectively.

## Predicting Ideal Intensities

```bash
siliscopy prop --method predI0 --file img --paramfile parameters.dat --timestep 100 --iterations 10
```
This predicts maximum intensities for next 10 iterations. These maximum intensity values are suggestions, and any value other than the predictions can be used.

For more details refer to the section "Choosing maximum fluorescence intensity and FWHM scaling factor" in Supporting information. ([here](https://www.biorxiv.org/content/10.1101/2021.03.02.433395v1))

## Determining histogram of Intensities

```bash
siliscopy prop --method hist --file img --paramfile parameters.dat --timestep 100 --output hist.dat
```

This outputs the histogram in the file hist.dat, where bin width is 0.1. To use a variable bin width, calculate it using the function siliscopy.prop.geti\_hist(). See the documentation for more details.

## Calculating Number of particles and area

Number and area from a monochrome 2D TIFF.
 
```bash
siliscopy prop --method num_area --file img10_lam670_fs530_T1_I0.13.tiff --paramfile parameters.dat --threshold 0.4 --output area_mono2d.dat
```

Number and area from a monochrome 2DT TIFF.
 
```bash
siliscopy prop --method num_area --file img10-111_lam670_fs530_T1_I0.13.tiff --paramfile parameters.dat --threshold 0.4 --output area_mono2dt.dat
```

Number and area from a monochrome 3D TIFF.
 
```bash
siliscopy prop --method num_area --file img10_z_lam670_fs530_T1_I0.13.tiff --paramfile parameters.dat --threshold 0.4 --output area_mono3d.dat
```

Number and area from a monochrome 3DT TIFF.
 
```bash
siliscopy prop --method num_area --file img10-111_z_lam670_fs530_T1_I0.13.tiff --paramfile parameters.dat --threshold 0.4 --output area_mono3dt.dat
```

Note: To use color TIFF images, use `rgb` or `nomix` color mixing scheme. Code not supported for PNG images.

## Calculating Number of particles and volume


Number and volume from a monochrome 3D TIFF.
 
```bash
siliscopy prop --method num_vol --file img10_z_lam670_fs530_T1_I0.13.tiff --paramfile parameters.dat --threshold 0.6 --output vol_mono3d.dat
```

Number and area from a monochrome 3DT TIFF.
 
```bash
siliscopy prop --method num_vol --file img10-111_z_lam670_fs530_T1_I0.13.tiff --paramfile parameters.dat --threshold 0.6 --output vol_mono3dt.dat
```

Note: To use color TIFF images, use `rgb` or `nomix` color mixing scheme. Code not supported for PNG images.

