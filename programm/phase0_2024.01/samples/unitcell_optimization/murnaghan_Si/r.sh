#!/bin/sh
NP=1
NE=1
NK=1
nvol=20
ivol=990
dvol=15
rm -f nfefn.data
for v in `seq 1 $nvol`;do
vol=$( echo \($v-1\)*$dvol+$ivol | bc -l )
cd vol$vol
mpiexec -n $NP ../../../../bin/phase ne=$NE nk=$NK < /dev/null
echo -n $vol >> ../nfefn.data; echo -n; tail -1 nfefn.data >> ../nfefn.data
cd ../
done

