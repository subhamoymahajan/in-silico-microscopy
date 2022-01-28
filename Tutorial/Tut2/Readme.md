# Tutorial 2: Generating Video with different Exposure time

## 1. Extract GRO files
```bash
tar -xzf Struct.tar.gz
```

## 2. Generate Image dat
```bash
siliscopy gen_mono --data imggen.dat --multiprocess
```

## 3. Render color images
```bash
siliscopy plot --file img --paramfile parameters.dat --method color --calc all --multiprocess --type jpeg
siliscopy plot --file img --paramfile parameters2.dat --method color --calc all --multiprocess --type jpeg
```

## 4. Create color videos
```bash
siliscopy video --file img --paramfile parameters.dat --method color 
siliscopy video --file img --paramfile parameters2.dat --method color 
```

Since .mov movies are created `fpns` is essentially the FPS of the video.
