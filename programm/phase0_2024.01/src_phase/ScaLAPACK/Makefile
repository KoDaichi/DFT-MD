#F77 = ifort
#F77FLAGS = -c
#AR = ar vq

OBJ = dlapst.o   pdlaed1.o  pdlaed3.o  pdlasrt.o  pdsyevd.o  pzheevd.o \
      pdlaed0.o  pdlaed2.o  pdlaedz.o  pdstedc.o  pzheev.o   zsteqr2.o

../libscalapack.a: $(OBJ)
	$(AR) $@ $?

clean:
	rm -f *.o

.f.o:
	$(F77) $(F77FLAGS) $*.f
