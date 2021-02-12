Generating my first in-silico microscopy image.
1) Generate the the PSF
term$ python run_genpsf.py

It will create two PSF for wavelength 670 nm and 518 nm. The code is currently slow. I will work on GPU accelerations (or hopefully someone else can help me with that).

2) Generate in-silico monochrome images.
a) Image data files
term$ ../../gen_mono -p parameters.dat -f dp100.gro -o img100
term$ ../../gen_mono -p parameters.dat -f dp2000.gro -o img2000

b) Render grey-scale images
term$ python ../../render_mono.py -f img -p png_param.dat -t 100
term$ python ../../render_mono.py -f img -p png_param.dat -t 2000

3) Generate colored in-silico monochrome image.
term$ python ../../mono2color.py -f img -p png_param.dat -t 100
term$ python ../../mono2color.py -f img -p png_param.dat -t 2000



