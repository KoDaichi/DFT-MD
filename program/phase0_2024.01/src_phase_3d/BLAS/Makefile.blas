.SUFFIXES:
.SUFFIXES: .o .f

BLAS = dasum.o dgemm.o drot.o dsyr2k.o idamax.o zdscal.o zher2.o ztrmv.o \
daxpy.o   dgemv.o   dscal.o  dtrmm.o   xerbla.o  zgemm.o   zher2k.o \
dcabs1.o  dger.o    dswap.o  dtrmv.o   zaxpy.o   zgemv.o   zscal.o \
dcopy.o   dlamch.o  dsymv.o  dtrsm.o   zcopy.o   zgerc.o   zswap.o \
ddot.o    dnrm2.o   dsyr2.o  dznrm2.o  zdotc.o   zhemv.o   ztrmm.o \
dspr.o    dspr2.o   dtpmv.o  dtpsv.o   dspmv.o   zhpr.o    zhpr2.o \
zhpmv.o   ztpmv.o   ztpsv.o \
ztrsm.o   zgeru.o   izamax.o \
zhemm.o   zherk.o   ztrsv.o

../libblas.a: $(BLAS)
	$(AR) $@ $?

clean:
	rm -f *.o

.f.o:
	$(F77) -c $(F77FLAGS) $*.f
