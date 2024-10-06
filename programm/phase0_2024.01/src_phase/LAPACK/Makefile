.SUFFIXES:
.SUFFIXES: .o .f

COMPLEX = zhetrd.o  zlaed0.o  zlarfb.o  zlasr.o   zsteqr.o  zungtr.o  zunmtr.o \
zdrot.o   zlacgv.o  zlaed7.o  zlarfg.o  zlassq.o  zung2l.o  zunm2l.o \
zheevd.o  zlacpy.o  zlaed8.o  zlarft.o  zlatrd.o  zung2r.o  zunm2r.o \
zheevx.o  zlacrm.o  zlanhe.o  zlascl.o  zstedc.o  zungql.o  zunmql.o \
zhetd2.o  zladiv.o  zlarf.o   zlaset.o  zstein.o  zungqr.o  zunmqr.o \
zhpevx.o  zhpgvx.o  zupgtr.o  zhpgst.o  zhptrd.o \
zlanhp.o  zpptrf.o  zupmtr.o  zheev.o \
zgesv.o   zgetrf.o  zgetri.o  zgetrs.o  zgetf2.o  zlaswp.o  ztrtri.o  ztrti2.o \
zhegvx.o  zpotrf.o  zhegst.o  zhegs2.o  zpotf2.o

DOUBLE = dlaebz.o  dlaed7.o  dlamrg.o  dlarft.o  dlassq.o  dorm2l.o  dsteqr.o \
dgesv.o   dlaed0.o  dlaed8.o  dlanst.o  dlarnv.o  dlaswp.o  dorm2r.o  dsterf.o \
dgetf2.o  dlaed1.o  dlaed9.o  dlansy.o  dlartg.o  dlatrd.o  dormql.o  dsyevd.o \
dgetrf.o  dlaed2.o  dlaeda.o  dlapy2.o  dlaruv.o  dorg2l.o  dormqr.o  dsyevx.o \
dgetrs.o  dlaed3.o  dlaev2.o  dlapy3.o  dlascl.o  dorg2r.o  dormtr.o  dsytd2.o \
dlacpy.o  dlaed4.o  dlagtf.o  dlarf.o   dlaset.o  dorgql.o  dstebz.o  dsytrd.o \
dladiv.o  dlaed5.o  dlagts.o  dlarfb.o  dlasr.o   dorgqr.o  dstedc.o \
dlae2.o   dlaed6.o  dlamch.o  dlarfg.o  dlasrt.o  dorgtr.o  dstein.o \
dlansp.o  dopmtr.o  dpptrf.o  dspgvx.o  dspevx.o  dsptrd.o \
dspgst.o  dopgtr.o  dsyev.o   dpotf2.o  dpotrf.o  dpotri.o \
dlauu2.o  dlauum.o  dspev.o   dsyrk.o   dtrtri.o  dtrti2.o \
dposv.o   dpotrs.o

UTIL =  ieeeck.o  ilaenv.o lsame.o   disnan.o  dlaisnan.o

../liblapack.a: $(COMPLEX) $(DOUBLE) $(UTIL)
	$(AR) $@ $?

clean:
	rm -f *.o

.f.o:
	$(F77) -c $(F77FLAGS) $*.f
