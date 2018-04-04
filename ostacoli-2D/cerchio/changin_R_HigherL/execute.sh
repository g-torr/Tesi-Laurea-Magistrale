#!/bin/bash
if false
then
cd ./R=0.5
python -c 'from analisi import*;pressione(0,9)'&
cd ..
cd ./R=0.75
python -c 'from analisi import*;pressione(0,9)'&
cd ..

cd ./R=1
python -c 'from analisi import*;pressione(0,9)'&
cd ..

cd ./R=1.5
python -c 'from analisi import*;pressione(0,9)'&
cd ..

cd ./R=2
python -c 'from analisi import*;pressione(0,9)'&
cd ..
fi



cd ./R=3
pwd
nvcc -o inside inside.cu
./inside
python -c 'from analisi import*;pressione(0,9)'&
cd ..

cd ./R=4
pwd
nvcc -o inside inside.cu
./inside
python -c 'from analisi import*;pressione(0,9)'&
cd ..

cd ./R=6
pwd
nvcc -o inside inside.cu
./inside
python -c 'from analisi import*;pressione(0,9)'&
cd ..

cd ./R=8
pwd
nvcc -o inside inside.cu
./inside
python -c 'from analisi import*;pressione(0,9)'&
cd ..


cd ./R=12
pwd
nvcc -o inside inside.cu
./inside
python -c 'from analisi import*;pressione(0,9)'&
cd ..



