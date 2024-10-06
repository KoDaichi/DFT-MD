#!/bin/bash
erange=" -erange -13,9" #energy range
einc="-einc 1.0" #scale of energy tics
imrfont="-imrfont Helvetica,18"#font of irreducible representation
ticsfont="-ticsfont Helvetica,10"
#color of irreducible representation (two methods)
#imrcolor="-imrcolor #545422" #1st way
imrcolor2="-imrcolor blue" #2nd way
#color of line (two methods)
#linecolor="-linecolor #545422" # 1st way
linecolor2="-linecolor red" #2nd way
imrtype="-imrtype mulliken" #type of representation
#nonecross="-nonecross" #turn off band crossing method 
#kgrouptype="-kgrouptype Schoenflies" #kgroup type
withfermi="-with_fermi nfefermi.data" #using fermi energy
numimr="-numimr 0.15" 
offset="-offset 0.5"
perl band_symm.pl reduce.data ${erange} ${einc} ${imrfont} ${ticsfont} ${imrcolor} ${imrcolor2} ${linecolor} ${linecolor2} ${imrtype} ${nonecross} ${kgrouptype} ${withfermi} ${numimr} ${offset}
#perl ${HOME}/band_symm/perl/band_symm.pl reduce.data ${erange} ${einc} ${imrfont} ${ticsfont} ${imrcolor} ${imrcolor2} ${linecolor} ${linecolor2} ${imrtype} ${noncross} ${kgrouptype} ${withfermi} ${numimr}
