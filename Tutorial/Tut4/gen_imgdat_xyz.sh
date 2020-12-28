#!/bin/bash
dir_name=(x y z)

for dir in 0 1 2
do
    sed "s/opt_axis =.*/opt_axis = $dir/g" parameters.dat > foo.dat 
    for val in 3 6 9 
    do
         sed -i "s/focus_cor =.*/focus_cor = $val/g" foo.dat
         ../../gen_mono -f dp100.gro -p foo.dat -o img_${dir_name[$dir]}${val}_100
    done
done

