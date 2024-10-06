#imrcolor,lincolorは2通りの指定方法(linecolor,linecolor2)のうち片方を使う
erange=" -erange -13,9"
einc="-einc 1.0"
imrfont="-imrfont Helvetica,14"
#ticsfont="-ticsfont Helvetica,12"
#imrcolor="-imrcolor #545422"
imrcolor2="-imrcolor blue"
#linecolor="-linecolor #545422"
linecolor2="-linecolor red"
imrtype="-imrtype mulliken"
#nonecross="-nonecross"
#kgrouptype="-kgrouptype HermannMauguin"
withfermi="-with_fermi ./nfefermi.data"
numimr="-numimr 0.15"
perl ~/phase/band_symm.09.2016test/perl/band_symm.pl reduce.data ${erange} ${einc} ${imrfont} ${ticsfont} ${imrcolor} ${imrcolor2} ${linecolor} ${linecolor2} ${imrtype} ${nonecross} ${kgrouptype} ${withfermi} ${numimr}
#echo "${erange} ${einc} ${imrfont} ${ticsfont} ${imrcolor} ${imrcolor2} ${linecolor} ${linecolor2} ${imrtype} ${nonecross} ${kgrouptype} ${withfermi} ${numimr}"
