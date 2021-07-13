# Tutorial 8: Generating images with different data-type

## Generate 2D image

PNG:

```bash
siliscopy plot --file img --paramfile parameters.dat --method mono \
               --timestep 10 --calc specific --output img --type png 
siliscopy plot --file img --paramfile parameters.dat --method color \
               --timestep 10 --calc specific --output img --type png 
```

TIFF:

```bash
siliscopy plot --file img --paramfile parameters.dat --method mono \
               --timestep 10 --calc specific --output img --type tiff8 
siliscopy plot --file img --paramfile parameters.dat --method color \
               --timestep 10 --calc specific --output img --type tiff8 
```

TIFF without color mixing:

```bash
siliscopy plot --file img --paramfile parameters2.dat --method color \
               --timestep 10 --calc specific --output nomix_img --type tiff8 
```

## Generate 2Dt image

```bash
siliscopy plot --file img --paramfile parameters.dat --method mono2dt \
               --calc specific --output img --type tiff8 
siliscopy plot --file img --paramfile parameters.dat --method color2dt \
               --calc specific --output img --type tiff8 
```

## Generate 3D image

```bash
siliscopy plot --file img --paramfile parameters.dat --method mono3d \
               --timestep 10 --calc specific --output img --type tiff8 \
               --multiprocess
siliscopy plot --file img --paramfile parameters.dat --method color3d \
               --timestep 10 --calc specific --output img --type tiff8 \
               --multiprocess
```

## Generate 3Dt image

```bash
siliscopy plot --file img --paramfile parameters.dat --method mono3dt \
               --calc specific --output img --type tiff8 --multiprocess
siliscopy plot --file img --paramfile parameters.dat --method color3dt \
              --calc specific --output img --type tiff8 --multiprocess
```


## Images with noise

2D PNG:

```bash
siliscopy plot --file img --paramfile param_noise.dat --method noise_mono \
               --timestep 10 --calc specific --output noise_img --type png 
siliscopy plot --file img --paramfile param_noise.dat --method noise_color \
               --timestep 10 --calc specific --output noise_img --type png 
```

2D TIFF:

```bash
siliscopy plot --file img --paramfile param_noise.dat --method noise_mono \
               --timestep 10 --calc specific --output noise_img --type tiff8 
siliscopy plot --file img --paramfile param_noise.dat --method noise_color \
               --timestep 10 --calc specific --output noise_img --type tiff8 
```

2D TIFF without color mixing

```bash
siliscopy plot --file img --paramfile param_noise2.dat --method noise_color \
               --timestep 10 --calc specific --output noise_nomix_img --type tiff8 
```

2DT image:

```bash
siliscopy plot --file img --paramfile param_noise.dat --method noise_mono2dt \
               --calc specific --output noise_img --type tiff8 
siliscopy plot --file img --paramfile param_noise.dat --method noise_color2dt \
               --calc specific --output noise_img --type tiff8 
```

3D image:

```bash
siliscopy plot --file img --paramfile param_noise.dat --method noise_mono3d \
               --timestep 10 --calc specific --output noise_img --type tiff8 \
               --multiprocess
siliscopy plot --file img --paramfile param_noise.dat --method noise_color3d \
               --timestep 10 --calc specific --output noise_img --type tiff8 \
               --multiprocess
```

3Dt image:

```bash
siliscopy plot --file img --paramfile param_noise.dat --method noise_mono3dt \
               --calc specific --output noise_img --type tiff8 --multiprocess
siliscopy plot --file img --paramfile param_noise.dat --method noise_color3dt \
               --calc specific --output noise_img --type tiff8 --multiprocess
```
