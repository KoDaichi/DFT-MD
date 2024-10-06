.SUFFIXES:
.SUFFIXES: .o .F .f .F90 .f90 .c .mod

ifndef F90
F90 = gfortran
endif
ifndef CC
CC  = gcc 
endif
LINK = $(F90)
F90FLAGS = -I./
CFLAGS = -O -DINTEL
CPPFLAGS = -DLinux

OBJ_INPUTPARSE = input_parse.o
OBJS_BM    =  m_Const_Parameters.o input_interface.o bluemoon.o 

all : bluemoon 

bluemoon : $(OBJS_BM) $(OBJ_INPUTPARSE)
	$(LINK) $(OBJS_BM) $(OBJ_INPUTPARSE) -o $@

$(OBJ_INPUTPARSE):$(@:.o=.c) $(@:.o=.h)
	$(CC) -o $*.o -c $(CFLAGS) $(@:.o=.c)

.f.o:
	$(F90) -o $*.o -c $(F77FLAGS) $*.f

.f90.o:
	$(F90) -o $*.o -c $(F90FLAGS) $*.f90

.F.o:
	$(F90) -o $*.o -c $(F77FLAGS) $(CPPFLAGS) $*.F

.F90.o:
	$(F90) -o $*.o -c $(F90FLAGS) $(CPPFLAGS) $*.F90

clean:
	\rm -f *.o *.mod *.a *.lib *.L *.list bluemoon

install: bluemoon
	\mv -f bluemoon ../bin/
