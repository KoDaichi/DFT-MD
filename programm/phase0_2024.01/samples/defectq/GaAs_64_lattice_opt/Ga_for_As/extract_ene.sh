#!/bin/sh
for n in q_*; 
do
#    echo $n
    n1=`echo $n | sed -e 's/q_//g'`
    ene=`tail -1 $n/nfefn.data | awk '{print $3}'`
    echo $n1 $ene
done
