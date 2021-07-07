#!/bin/bash

for PN in 1 5 25
do
     sed "s/Plmn\s*=.*/Plmn = 15\.0 15\.0 ${PN}\.0/g" parameters.dat > foo.dat 
     siliscopy gen_mono --file dp100.gro --paramfile foo.dat --psf PSF_gandy \
                        --output Pn${PN}_img100 --method slice
done
rm foo.dat
