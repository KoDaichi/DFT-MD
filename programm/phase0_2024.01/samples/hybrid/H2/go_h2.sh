#!/bin/sh

phase=$HOME/phase_v1001/bin/phase

# (1) PBE DFT
cd PBE
$phase >& log

# (2) HSE06 Hybrid DFT
cd ../HSE06
cp ../PBE/nfchgt.data .
cp ../PBE/zaj.data .
$phase >& log

# (3) PBE0 Hybrid DFT
cd ../PBE0
cp ../PBE/nfchgt.data .
cp ../PBE/zaj.data .
$phase >& log

# (4) Hartree-Fock
cd ../HF
cp ../PBE/nfchgt.data .
cp ../PBE/zaj.data .
$phase >& log

