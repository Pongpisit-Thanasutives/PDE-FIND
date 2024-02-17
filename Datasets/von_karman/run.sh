#!/usr/bin/env bash

# Run at Re=100
./ibpm/build/ibpm -ic von_karman/ibpm15000_start.bin -Re 100 -nx 450 -ny 200 -ngrid 4 -length 9 -xoffset -1 -yoffset -2 -xshift 0.75 -nsteps 1500 -restart 100 -tecplot 10 -outdir von_karman/ -geom von_karman/cylinder.geom -dt .02 -ubf 0 -tecplotallgrids 0
