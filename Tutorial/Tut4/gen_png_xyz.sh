#!/bin/bash
dir_name=(x y z)

for dir in 0 1 2
do
    for val in 3 6 9 
    do
         python ../../mono2color.py -f img_${dir_name[$dir]}${val}_ -p png_param.dat -t 100
    done
done

