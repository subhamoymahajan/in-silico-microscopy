#!/bin/bash 

for PN in 1 10 20 30 40 50
do
      python ../../mono2color.py -f Pn${PN}_img -p png_param.dat -t 100
done
