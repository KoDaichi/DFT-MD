.SUFFIXES:
.SUFFIXES: .o .F .f .F90 .f90 .c .mod

CPPFLAGS = -DLinux

LOWER_MODULES = m_Const_Parameters.o m_ArraySize_Parameters.o

UPPER_MODULES =  m_Control_Parameters.o  m_Charge_File.o  m_Files.o   m_Timing.o\
m_Crystal_Structure.o       m_FFT.o              m_Kpoints.o\
m_PlaneWaveBasisSet.o  m_Electronic_Structure.o 

F_SUBROUTINES = Initialization.o      InputData_Analysis.o \
Preparation.o  rd_WFs_doFFT_and_solve_eq.o  Finalization.o \
miscellaneous.o 

F_OTHERSUBS   = stm_main.o bottom_Subroutines.o

FFT           = jrcat_ft_stm.o

OBJECTSF77 =  decfft_ent.o

TARGET   = stm

COMP=ifort

include makefile.$(COMP)

all : ${TARGET}

stm : $(LOWER_MODULES) $(UPPER_MODULES) $(F_SUBROUTINES) $(F_OTHERSUBS) $(FFT)
	$(LINK) $(LOWER_MODULES) $(UPPER_MODULES) $(F_SUBROUTINES) $(F_OTHERSUBS) $(FFT) -o $@

clean:
	\rm -f *.o  *.mod ${TARGET}

install:
	mv ${TARGET} ../bin/

.f.o:
	$(F90) -c $(F77FLAGS) $*.f

.f90.o:
	$(F90) -c $(F90FLAGS) $*.f90

.F.o:
	$(F90) -c $(F77FLAGS) $(CPPFLAGS) $*.F

.F90.o:
	$(F90) -c $(F90FLAGS) $(CPPFLAGS) $*.F90

