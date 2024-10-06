
#FC = gfortran
#FFLAGS = -O

ifdef F90
FC=${F90}
else
FC = gfortran
endif
FFLAGS = -O

CPPFLAGS =

LIBS =

LOWER_MODULES = m_Const_Parameters.o m_ArraySize_Parameters.o

UPPER_MODULES = m_Files.o m_Kpoints.o  m_Control_Parameters.o  m_Charge_File.o\
	m_Crystal_Structure.o       m_FFT.o \
	m_PlaneWaveBasisSet.o  m_Electronic_Structure.o

F_SUBROUTINES = Initialization.o      InputData_Analysis.o \
	Preparation.o  rd_WFs_doFFT_and_solve_eq.o  Finalization.o

F_OTHERSUBS   = stm_main.o bottom_Subroutines.o

FFT           = jrcat_ft_stm.o

#OBJECTSF77 =  decfft_ent.o

OBJS = $(LOWER_MODULES) $(UPPER_MODULES) $(F_SUBROUTINES) $(F_OTHERSUBS) $(FFT)

TARGET = workfunc

all : $(TARGET)

workfunc : $(OBJS)
	$(FC) -o $@ $(OBJS) $(LIBS)

.SUFFIXES: .o .F .f .F90 .f90
.f.o:
	$(FC) -c $(FFLAGS) $*.f
.f90.o:
	$(FC) -c $(FFLAGS) $*.f90
.F.o:
	$(FC) -c $(FFLAGS) $(CPPFLAGS) $*.F
.F90.o:
	$(FC) -c $(FFLAGS) $(CPPFLAGS) $*.F90

clean:
	\rm -f *.o *.mod ${TARGET}

install: $(TARGET)
	mv ${TARGET} ../bin/

