#!/bin/sh
nc=50
ic=5.0
dc=0.25
for v in `seq 1 $nc`; do
  c=$( echo \($v-1\)*$dc+$ic | bc -l )
  mkdir -p c$c
  cp -r file_names.data c$c/
  sed "s/__C__/${c}/g" nfinp.data > c$c/nfinp.data
  cd c$c/
  mpiexec -n 1 ../../../../../bin/phase
  cd ../
  echo -n $c >> nfefn.data ; tail -1 c$c/nfefn.data >> nfefn.data
done

