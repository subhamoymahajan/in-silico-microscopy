#!/bin/bash

for PN in 1 10 20 30 40 50
do
     sed "s/Lpsf =.*/Lpsf = 15\.0 15\.0 ${PN}\.0/g" parameters.dat > foo.dat 
     ../../gen_mono -f dp100.gro -p foo.dat -o Pn${PN}_img100
done
