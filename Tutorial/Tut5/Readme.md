Generating my first in-silico microscopy image.

3) Generate colored in-silico monochrome image.
#Red-green
term$ python ../../mono2color.py -f img100 -p png_param_rg.dat -t -1
term$ mv img1000.png img_rg.png
#Orange-Violet
term$ python ../../mono2color.py -f img100 -p png_param_ov.dat -t -1
term$ mv img1000.png img_ov.png
#Cyan-Magenta
term$ python ../../mono2color.py -f img100 -p png_param_cm.dat -t -1
term$ mv img1000.png img_cm.png



