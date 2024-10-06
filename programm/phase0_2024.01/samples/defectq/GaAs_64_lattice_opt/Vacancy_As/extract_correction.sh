#!/bin/sh
for n in q_*; 
do
#    echo $n
    n1=`echo $n | sed -e 's/q_//g'`
    ene=0.0
    if [ -f $n/defect*atoms ]; then
	ene=`tail -1 $n/defect*atoms | awk '{print $NF}'`
    fi
    echo $n1 $ene
done
