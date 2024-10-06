#!/bin/bash

#PBS -q S
#PBS -l ncpus=1
#PBS -l nodes=1:ppn=1

cd $PBS_O_WORKDIR
#phase=$HOME/phase_v900/bin/phase
phase=/s0/home2/jkoga/phase_dev/phase_beta/phase_v900/bin/phase

# (1) DFT
cd PBE
$phase >& log

# (2) Hybrid DFT
cd ../PBE0
cp ../PBE/nfchgt.data .
cp ../PBE/zaj.data .
$phase >& log

