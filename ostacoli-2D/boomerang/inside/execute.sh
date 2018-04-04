#!/bin/bash
cd /home/giuseppe/Documents/myCUDA/ostacoli-2D/boomerang/inside/tau=0.3
pwd
python -c 'from mkframes import*;media_modulo_corrente()'>>result&

cd /home/giuseppe/Documents/myCUDA/ostacoli-2D/boomerang/inside/tau=0.5
pwd
python -c 'from mkframes import*;media_modulo_corrente()'>>result&

cd /home/giuseppe/Documents/myCUDA/ostacoli-2D/boomerang/inside/tau=0.06
pwd
python -c 'from mkframes import*;media_modulo_corrente()'>>result&

cd /home/giuseppe/Documents/myCUDA/ostacoli-2D/boomerang/inside/tau=1
pwd
python -c 'from mkframes import*;media_modulo_corrente()'>>result&

cd /home/giuseppe/Documents/myCUDA/ostacoli-2D/boomerang/inside/tau=2
pwd
python -c 'from mkframes import*;media_modulo_corrente()'>>result&

cd /home/giuseppe/Documents/myCUDA/ostacoli-2D/boomerang/inside/tau=3
pwd
python -c 'from mkframes import*;media_modulo_corrente()'>>result&

cd /home/giuseppe/Documents/myCUDA/ostacoli-2D/boomerang/inside/tau=4
pwd
python -c 'from mkframes import*;media_modulo_corrente()'>>result&

cd /home/giuseppe/Documents/myCUDA/ostacoli-2D/boomerang/inside/tau=10
pwd
python -c 'from mkframes import*;media_modulo_corrente()'>>result&

cd /home/giuseppe/Documents/myCUDA/ostacoli-2D/boomerang/inside/tau=20
pwd
python -c 'from mkframes import*;media_modulo_corrente()'>>result&




