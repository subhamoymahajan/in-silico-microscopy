#!/bin/bash

sed 's/fs=.*/fs=600/g' png_param.dat > foo.dat # change fs to 600
for I0 in 0.1 0.2 0.3
do
     sed -i "s/lam1_I0=.*/lam1_I0=$I0/g" foo.dat # change the value of I0
     sed -i "s/lam2_I0=.*/lam2_I0=$I0/g" foo.dat # change the value of I0
     python ../../mono2color.py -f img -p foo.dat -t 100
done

sed 's/fs=.*/fs=800/g' png_param.dat > foo.dat # change fs to 800
for I0 in 0.1 0.2 0.3
do
     sed -i "s/lam1_I0=.*/lam1_I0=$I0/g" foo.dat # change the value of I0
     sed -i "s/lam2_I0=.*/lam2_I0=$I0/g" foo.dat # change the value of I0
     python ../../mono2color.py -f img -p foo.dat -t 100
done

rm foo.dat
