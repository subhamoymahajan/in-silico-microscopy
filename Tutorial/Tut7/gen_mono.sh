#!/bin/bash

gen_mono (){
    echo "##################################"
    siliscopy gen_mono --file dp100.gro --paramfile foo.dat --psf PSF_GL$1\
                       --output GL$1_img100 --method slice
}

#Design condition
cp parameters.dat foo.dat
gen_mono

#Changing tsO
sed "s/tsO\s*=.*/tsO = 10/g" parameters.dat  > foo.dat
gen_mono
sed "s/tsO\s*=.*/tsO = 20/g" parameters.dat > foo.dat
gen_mono

#Difference between meu and meu0
sed "s/meu\s*=.*/meu = 1.6/g" parameters.dat > foo.dat
gen_mono "_meu1.6"
sed "s/meu\s*=.*/meu = 1.8/g" parameters.dat > foo.dat
gen_mono "_meu1.8"

#Difference between meug and meug0
sed "s/meug\s*=.*/meug = 1.6/g" parameters.dat > foo.dat
gen_mono "_meug1.6"
sed "s/meug\s*=.*/meug = 1.8/g" parameters.dat > foo.dat
gen_mono "_meug1.8"

#Changing meus
sed "s/meus\s*=.*/meus = 1.45/g" parameters.dat > foo.dat 
gen_mono "_meus1.45"
sed "s/meus\s*=.*/meus = 1.55/g" parameters.dat > foo.dat
gen_mono "_meus1.55"

#Difference between tg and tg0
sed "s/tg\s*=.*/tg = 220/g" parameters.dat > foo.dat
gen_mono "_tg220"
sed "s/tg\s*=.*/tg = 420/g" parameters.dat > foo.dat
gen_mono "_tg420"

rm foo.dat
