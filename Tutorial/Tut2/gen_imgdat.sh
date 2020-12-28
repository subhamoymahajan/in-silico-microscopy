#!/bin/bash
# Roughly takes 5 mins
for i in {10..110}
do
    ../../gen_mono -p parameters.dat -f dp${i}.gro -o img${i}
done

