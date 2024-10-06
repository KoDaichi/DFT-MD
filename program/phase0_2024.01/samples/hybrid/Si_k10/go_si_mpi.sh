#!/bin/sh

phase=$HOME/phase_v1001/bin/phase

# (1) PBE DFT
cd PBE
mpirun -np 47 $phase ne=1 nk=47 >& log

# (2) HSE06 hybrid DFT
cd ../HSE06
cp ../PBE/nfchgt.data .
cp ../PBE/zaj.data .
mpirun -np 47 $phase ne=1 nk=47 >& log

