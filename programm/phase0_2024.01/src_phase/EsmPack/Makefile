.SUFFIXES: .o .f90 .F90

FFLAG   = $(INCLUDE) $(OMPFLAG) $(MPIFLAG)

OBJ = qe_erf.o fft.o vector.o Ewald.o  Esm.o EsmInterface.o
OBJ_P = qe_erf.o fft.o vector.o Ewald.o  Esm.o EsmPack.o
LIBFLAG = -L/usr/local/lib -lfftw3 -lblas -llapack

.f90.o:
	$(FORTRAN) -c $(FFLAG) $<
.F90.o:
	$(FORTRAN) -c $(FFLAG) $<

all: ../libesm.a

../libesm.a: $(OBJ)
	$(AR) $@ $?

exe: $(OBJ_P) 
	$(FORTRAN) $(FFLAG) $(LIBFLAG) $(OBJ_P) -o EsmPack.x  

clean:
	rm -f $(OBJ) $(OBJ_P) ../libesm.a EsmPack.x  *.mod
