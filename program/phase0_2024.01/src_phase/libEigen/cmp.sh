#!/bin/sh -xe

cd eigen_h
make
cd ../eigen_s32
make
cd ../eigen_sx1
make
cd ../eigen_ex1
make
cd ..
ar r libEigen.a eigen_h/*.o eigen_s32/*.o eigen_sx1/*.o eigen_ex1/*.o
ranlib libEigen.a
mv eigen_*/lib* .
