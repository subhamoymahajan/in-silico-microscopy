#!/bin/bash

for i in {10..100..10}
do
    python ../../mono2color.py -f img -p png_param.dat -t ${i}
done

for i in {10..100..10}
do
    python ../../mono2color.py -f img -p png_param_T10.dat -t ${i}
done
