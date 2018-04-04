#!/bin/bash
if false
then
cd ./R=1
pwd
nvcc -o inside inside.cu --use_fast_math --maxrregcount=32
./inside
python -c 'from analisi import*;pressione(0,9)'&
cd ..
cd ./R=1.5
pwd
nvcc -o inside inside.cu --use_fast_math --maxrregcount=32
./inside
python -c 'from analisi import*;pressione(0,9)'&
cd ..

cd ./R=2
pwd
nvcc -o inside inside.cu --use_fast_math --maxrregcount=32
./inside
python -c 'from analisi import*;pressione(0,9)'&
cd ..

cd ./R=3
pwd
nvcc -o inside inside.cu --use_fast_math --maxrregcount=32
./inside
python -c 'from analisi import*;pressione(0,9)'&
cd ..

cd ./R=4
pwd
nvcc -o inside inside.cu --use_fast_math --maxrregcount=32
./inside
python -c 'from analisi import*;pressione(0,9)'&
cd ..
then
cd ./R=6
pwd
nvcc -o inside inside.cu --use_fast_math --maxrregcount=32
./inside

python -c 'from analisi import*;pressione(0,9)'&
cd ..

if false 
then
cd ./R=8
pwd
nvcc -o inside inside.cu --use_fast_math --maxrregcount=32
./inside

python -c 'from analisi import*;pressione(0,9)'&
cd ..
fi
cd ./R=12
pwd
nvcc -o inside inside.cu --use_fast_math --maxrregcount=32
./inside

python -c 'from analisi import*;pressione(0,9)'&
cd ..




