#!/bin/bash
set -e
for fs in 400 530
do
    for I0 in 0.1 0.2 0.3
    do
         # change the value of lam_I0_1 and lam_I0_2
         sed "s/lam_I0_1\s*=.*/lam_I0_1=$I0/g" param_${fs}.dat > foo.dat 
         sed -i "s/lam_I0_2\s*=.*/lam_I0_2=$I0/g" foo.dat 
         siliscopy plot --file img --paramfile foo.dat --method color \
                        --timestep 100 --calc specific --type jpeg
    done
done 
rm foo.dat
