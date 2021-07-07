# Tutorial 1: First in-silico microscopy image

## 1. Generating PSF

```bash
siliscopy gen_psf --method gandy --paramfile parameters.dat --calc all --output PSF_gandy --multiprocess
```

## 2. Calculate Monochrome Intensities

```bash
siliscopy gen_mono --file dp100.gro --paramfile parameters.dat --psf PSF_gandy --output img100 --method slice
siliscopy gen_mono --file dp2000.gro --paramfile parameters.dat --psf PSF_gandy --output img2000 --method slice
```

## 3. Plot Monochome In-silico Images

```bash
siliscopy plot --file img --paramfile parameters.dat --method mono --timestep 100 --calc specific --type jpeg
siliscopy plot --file img --paramfile parameters.dat --method mono --timestep 2000 --calc specific --type jpeg 
```

## 4. Plot Colored In-silico Images
```bash
siliscopy plot --file img --paramfile parameters.dat --method color --timestep 100 --calc specific --type jpeg
siliscopy plot --file img --paramfile parameters.dat --method color --timestep 2000 --calc specific --type jpeg
```
