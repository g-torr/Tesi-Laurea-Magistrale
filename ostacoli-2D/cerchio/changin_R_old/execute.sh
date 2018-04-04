#!/bin/bash
if false
then
cd ./R=0.5
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

cd ./R=3
python -c 'from analisi import*;pressione(0,9)'&
cd ..
fi
cd ./R=4
python -c 'from analisi import*;pressione(0,9)'&
cd ..

cd ./R=6
python -c 'from analisi import*;pressione(0,9)'&
cd ..

cd ./R=8
python -c 'from analisi import*;pressione(0,9)'&
cd ..

cd ./R=12
python -c 'from analisi import*;pressione(0,9)'&
cd ..




