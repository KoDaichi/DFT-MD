
#FC = gfortran
#FFLAGS = -O

FC = ifort
FFLAGS = -O

CPPFLAGS =

LIBS = 

OBJS = spectrum.o

.SUFFIXES: .o .f90 .F90

.f90.o:
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $< -o $@
.F90.o:
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $< -o $@

all: spectrum

spectrum: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	-\rm -f *.o spectrum

install:
	mv spectrum ../bin/
