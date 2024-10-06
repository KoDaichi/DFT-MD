#!/bin/sh
#
E0=`tail -1 q_0/nfefn.data | awk '{print $3}'`
#
echo "# q  vbm [eV}"
for n in q*
do
    E1=`tail -1 $n/nfefn.data | awk '{print $3}'`
    dq=`echo $n | sed -e 's/q_//g'`
    if [ $dq != 0 ];then
	vbm=`echo $E0 $E1 $dq | awk '{print ($1-$2)/$3*27.2116}'`
	echo $dq $vbm
    fi
done

    
