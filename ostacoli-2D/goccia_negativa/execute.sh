#!/bin/bash
nvcc -o ./pi:8/inside ./pi:8/inside.cu 
cd pi:8/
pwd
./inside
cd ..
nvcc -o ./pi3:5/inside ./pi3:5/inside.cu 
cd pi3:5
pwd
./inside
cd ..
nvcc -o ./pi:4/inside ./pi:4/inside.cu 
cd pi:4
pwd
./inside

