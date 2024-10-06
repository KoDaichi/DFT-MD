#!/bin/sh
# Shell script for checking PHASE and EKCAL
# Copyright (C) T.Yamamoto

echo "Checking total-energy calculation."
cd scf
rm -f cont*
rm -f output*
../../bin/phase
etot=`grep TH output000 | tail -1 | awk '{print $7}'`
eref="-7.897015064593"
echo " Total energy : ${etot} Hartree/cell"
echo " Reference    : ${eref} Hartree/cell"
rm -f cont* nfdynm* nfefn* nfstop* zaj*
mv -f nfchgt.data ../dos/

echo "Checking band-energy calculation."
cd ../dos
rm -f output*
../../bin/ekcal
vbme=`grep Valence nfenergy.data | awk '{print $5}'`
eref="0.233847"
echo " Valence band maximum : ${vbme} Hartree"
echo " Reference            : ${eref} Hartree"
rm -f cont* dos* nfchgt* nfenergy* nfstop* zaj*
