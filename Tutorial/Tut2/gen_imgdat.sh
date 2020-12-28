#!/bin/bash

for i in {10..110}
do
    ../../gen_mono -p parameters.dat -f dp${i}.gro -o img${i}
done

