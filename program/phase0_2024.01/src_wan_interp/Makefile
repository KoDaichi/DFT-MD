
#FC = gfortran
#FFLAGS = -O

FC = ifort
FFLAGS = -O

CPPFLAGS =

# LAPACK
#LAPACK_LIB = -L/usr/local/lapack -llapack -lblas

# Intel MKL
LAPACK_LIB = -L/opt/intel/composerxe/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

LIBS = $(LAPACK_LIB)

OBJS = wan_interp_phase.o

.SUFFIXES: .o .f90 .F90

.f90.o:
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $< -o $@
.F90.o:
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $< -o $@

all: wan_interp

wan_interp: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	-\rm -f *.o wan_interp

install:
	mv wan_interp ../bin/
