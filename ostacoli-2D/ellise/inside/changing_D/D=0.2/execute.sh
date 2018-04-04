#!/bin/bash 
cd /home/giuseppe/Documents/myCUDA/ostacoli-2D/ellise/inside/D=0.2/tau=0.03
nvcc -o inside inside.cu  -maxrregcount 32 -lineinfo --use_fast_math
./inside
cd /home/giuseppe/Documents/myCUDA/ostacoli-2D/ellise/inside/D=0.2/tau=0.06
nvcc -o inside inside.cu  -maxrregcount 32 -lineinfo --use_fast_math
./inside
cd /home/giuseppe/Documents/myCUDA/ostacoli-2D/ellise/inside/D=0.2/tau=0.12
nvcc -o inside inside.cu  -maxrregcount 32 -lineinfo --use_fast_math
./inside
cd /home/giuseppe/Documents/myCUDA/ostacoli-2D/ellise/inside/D=0.2/tau=0.24
nvcc -o inside inside.cu  -maxrregcount 32 -lineinfo --use_fast_math
./inside
cd /home/giuseppe/Documents/myCUDA/ostacoli-2D/ellise/inside/D=0.2/tau=0.48
nvcc -o inside inside.cu  -maxrregcount 32 -lineinfo --use_fast_math
./inside
cd /home/giuseppe/Documents/myCUDA/ostacoli-2D/ellise/inside/D=0.2/tau=0.96
nvcc -o inside inside.cu  -maxrregcount 32 -lineinfo --use_fast_math
./inside

