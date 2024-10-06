###########################################################################
###<< PLEASE CHANGE THE VARIABLES BELOW ACCORDING TO YOUR ENVIRONMENT >>###
###########################################################################

FC = gfortran -std=legacy
FFLAG = -ffree-line-length-none

###########################################################################
###<< PLEASE CHANGE THE VARIABLES ABOVE ACCORDING TO YOUR ENVIRONMENT >>###
###########################################################################

#_USE_OPENMP=1
ifdef _USE_OPENMP
OMPFLG=-fopenmp
OMPLIB=-lgomp
endif

CPP = 
F77FLAGS = $(FFLAG)
F90FLAGS = $(FFLAG)
LINK=$(FC)
INC=
LIB=

.SUFFIXES:
.SUFFIXES: .o .F .f .F90 .f90 .c .mod

all:
.f.o:
	$(FC) -c $(F77FLAGS) $(OMPFLG) $*.f $(INC)

.f90.o:
	$(FC) -c $(F90FLAGS) $(OMPFLG) $*.f90 $(INC)

.F.o:
	$(FC) -c $(F77FLAGS) $(OMPFLG) -cpp $(CPP) $*.F $(INC)

.F90.o:
	$(FC) -c $(F90FLAGS) $(OMPFLG) -cpp $(CPP) $*.F90 $(INC)
.o:
	$(LINK) -o $@ $*.o $(LIB) $(OMPLIB)

clean:
	rm -f *.o *.mod *.a
erase:
	make clean
	rm -f *.out *.exe
