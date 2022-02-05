# Tutorial 5: Generating in silico microscopy images with different hues.

## Red-green
```bash
siliscopy plot --file img --paramfile param_rg.dat --method color --timestep 100 --calc specific --output img_rg_ --type jpeg
```
## Orange-Violet
```bash
siliscopy plot --file img --paramfile param_ov.dat --method color --timestep 100 --calc specific --output img_ov_ --type jpeg
```
## Cyan-Magenta
```bash
siliscopy plot --file img --paramfile param_cm.dat --method color --timestep 100 --calc specific --output img_cm_ --type jpeg
```



