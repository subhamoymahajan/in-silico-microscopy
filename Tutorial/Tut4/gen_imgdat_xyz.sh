#!/bin/bash
dir_name=(x y z)

for dir in 0 1 2
do
    sed "s/opt_axis\s*=.*/opt_axis = $dir/g" parameters.dat > foo.dat 
    for cor in 3 6 9 
    do
         sed -i "s/focus_cor\s*=.*/focus_cor = $cor/g" foo.dat
         siliscopy gen_mono --file dp100.gro --paramfile foo.dat \
                            --psf PSF_gandy --method slice \
                            --output img_${dir_name[$dir]}${cor}_100  
    done
done
rm foo.dat
