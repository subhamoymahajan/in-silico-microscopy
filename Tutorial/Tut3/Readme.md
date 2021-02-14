# Tutorial 3: Geneerating *in-silico* microscopy image with different resolution and brightness

## 1. Generate the the PSF
```bash
siliscopy gen_psf --method gandy --paramfile param_600.dat --calc all --output PSF_gandy --multiprocess
```

## 2. Generate in-silico monochrome image intensities.
```bash
siliscopy gen_mono --file dp100.gro --paramfile param_800.dat --psf PSF_gandy --output img100
siliscopy gen_mono --file dp100.gro --paramfile param_600.dat --psf PSF_gandy --output img100
```
## 3. Generate colored in-silico microscopy image.

Use the scrip below to generate images with different maximum intensity I0. 
```bash
bash gen_I0_fs.sh
```


