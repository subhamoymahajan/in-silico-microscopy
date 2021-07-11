#!/bin/bash

gen_psf (){
    echo "###################################"
    siliscopy gen_psf --method GL1991 --paramfile foo.dat --calc all\
                      --output PSF_GL$1 --multiprocess
}

#Design condition
cp parameters.dat foo.dat
gen_psf

#Changing tsO
sed "s/tsO\s*=.*/tsO = 10/g" parameters.dat  > foo.dat
gen_psf
sed "s/tsO\s*=.*/tsO = 20/g" parameters.dat > foo.dat
gen_psf

#Difference between meu and meu0
sed "s/meu\s*=.*/meu = 1.6/g" parameters.dat > foo.dat
gen_psf "_meu1.6"
sed "s/meu\s*=.*/meu = 1.8/g" parameters.dat > foo.dat
gen_psf "_meu1.8"

#Difference between meug and meug0
sed "s/meug\s*=.*/meug = 1.6/g" parameters.dat > foo.dat
gen_psf "_meug1.6"
sed "s/meug\s*=.*/meug = 1.8/g" parameters.dat > foo.dat
gen_psf "_meug1.8"

#Changing meus
sed "s/meus\s*=.*/meus = 1.45/g" parameters.dat > foo.dat 
gen_psf "_meus1.45"
sed "s/meus\s*=.*/meus = 1.55/g" parameters.dat > foo.dat
gen_psf "_meus1.55"

#Difference between tg and tg0
sed "s/tg\s*=.*/tg = 220/g" parameters.dat > foo.dat
gen_psf "_tg220"
sed "s/tg\s*=.*/tg = 420/g" parameters.dat > foo.dat
gen_psf "_tg420"

rm foo.dat
