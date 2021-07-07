#!/bin/bash 

for PN in 1 5 25
do
      siliscopy plot --file Pn${PN}_img --paramfile parameters.dat \
                     --method color --timestep 100 --calc specific --type jpeg 
done
