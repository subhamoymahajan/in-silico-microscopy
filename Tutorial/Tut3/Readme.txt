Generating my first in-silico microscopy image.
1) Generate the the PSF
term$ python run_genpsf.py

2) Generate in-silico monochrome images.
 Image data files
term$ ../../gen_mono -p param_800.dat -f dp100.gro -o img100
term$ ../../gen_mono -p param_600.dat -f dp100.gro -o img100

3) Generate colored in-silico microscopy image.
term$ bash gen_I0_fs.sh


