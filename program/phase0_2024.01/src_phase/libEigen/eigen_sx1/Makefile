MPIR_HOME = /usr/local/mpi
FF = mpifrtpx
FC = mpifrtpx
CC = mpifccpx

# Optimization flags for the Intel compiler
#CCFLAG_BASE = -g -axW -openmp
#CCFLAG  = $(CCFLAG_BASE) -vec-report3 # -ftrapuv -CU -CB -WB
CCFLAG = -Kfast,parallel,ocl,array_private,auto,openmp,preex,simd=2 -X9 -Cpp -DTIMER -DPHASE #-Haefos#u#x

#T2CCFLAG  = #-mP2OPT_hlo_pref_hint=3
#TXCCFLAG  = #-mP2OPT_hlo_pref_hint=3

#SCALAPACK = -L$(HOME)/lib -lscalapack-intel $(BLAS) -lblacsC-intel-a -lblacsF77-intel-a -lblacsBASE-intel-a -L/opt/intel/Compiler/11.0/074/mkl/lib/em64t -lmkl -lmkl_em64t -lblacsCinit -lblacsF77init -lblacs
#BLAS = -L/opt/intel/Compiler/11.0/074/mkl/lib/em64t -lmkl -lmkl_em64t
SCALAPACK = -SSL2BLAMP -SCALAPACK
BLAS =

INCPATH = -I/usr/local/include

INCFLAG = $(INCPATH)

LIB = Eigen_sx
AR     = ar
RANLIB = ranlib
RM     = rm -f

LDFLAG  = $(CCFLAG) $(LDPATH) $(SCALAPACK) $(SCALAPACK)
#LDFLAG  = -g $(LDPATH) -Kopenmp $(SCALAPACK) $(SCALAPACK)

OBJS = getrank26d.o eigen.o pdlaed2.o pdlaed3.o lapack_eigen.o my_pdlaedz.o my_pdlaed0.o my_pdsxedc.o my_pdlaed1.o trbakwy4.o trbakwy4_body.o machine_epsilon.o tred1.o t1.o t2.o t3.o t4.o t5.o t6.o t7.o t8.o CSTAB.o CSTAB_mem.o comm.o dcx.o
OBJSM = main2.o mat_set.o ev_test_2D.o

all : lib$(LIB).a eigen_sx

lib$(LIB).a: $(OBJS)
	$(RM) $@
	$(AR) r $@ $(OBJS)
	$(RANLIB) $@

eigen_sx : $(OBJSM)
	$(FF) -o eigen_sx $(LDFLAG) $(OBJSM) -L. -lEigen_sx

main2.o : main2.f
	$(FF) -c main2.f $(CCFLAG) $(INCFLAG)
tred1.o : tred1.f
	$(FF) -c tred1.f $(CCFLAG) $(INCFLAG)
t1.o : t1.f
	$(FF) -c t1.f $(CCFLAG) $(INCFLAG)
t2.o : t2.f
	$(FC) -c t2.f $(CCFLAG) $(T2CCFLAG) $(INCFLAG)
t3.o : t3.f
	$(FC) -c t3.f $(CCFLAG) $(INCFLAG)
t4.o : t4.f
	$(FC) -c t4.f $(CCFLAG) $(T2CCFLAG) $(INCFLAG)
t5.o : t5.f
	$(FC) -c t5.f $(CCFLAG) $(INCFLAG)
t6.o : t6.f
	$(FC) -c t6.f $(CCFLAG) $(INCFLAG)
t7.o : t7.f
	$(FC) -c t7.f $(CCFLAG) $(INCFLAG)
t8.o : t8.f
	$(FC) -c t8.f $(CCFLAG) $(INCFLAG)
trbakwy4.o : trbakwy4.f
	$(FF) -c trbakwy4.f $(CCFLAG)  $(INCFLAG)
trbakwy4_body.o : trbakwy4_body.f
	$(FF) -c trbakwy4_body.f $(CCFLAG)  $(INCFLAG)
machine_epsilon.o : machine_epsilon.f
	$(FF) -c machine_epsilon.f -O0 $(CCFLAG) -O0 $(INCFLAG)
mat_set.o: mat_set.f
	$(FF) -c mat_set.f $(CCFLAG) $(INCFLAG)
ev_test.o: ev_test.f
	$(FF) -c ev_test.f $(CCFLAG) $(INCFLAG)
dlaed3.o: dlaed3.f
	$(FF) -c dlaed3.f $(CCFLAG) $(INCFLAG)
dlaed4.o: dlaed4.f
	$(FF) -c dlaed4.f $(CCFLAG) $(INCFLAG)
dlaed6.o: dlaed6.f
	$(FF) -c dlaed6.f $(CCFLAG) $(INCFLAG)
dcx.o: dcx.f
	$(FF) -c dcx.f $(CCFLAG) $(INCFLAG)
comm.o : comm.f
	$(FF) -c comm.f  $(CCFLAG) $(INCFLAG)
eigen.o : eigen.f
	$(FF) -c eigen.f  $(CCFLAG) $(INCFLAG)

clean :
	-\rm *.o eigen_sx  lib$(LIB).a

trbakwy4.o: CSTAB.h
htrd1.o: CSTAB.h
tred1.o: CSTAB.h

lapack_eigen.o : lapack_eigen.f
	$(FF) -c lapack_eigen.f $(CCFLAG) $(INCFLAG)
my_pdlaedz.o : my_pdlaedz.f
	$(FF) -c my_pdlaedz.f $(CCFLAG) $(INCFLAG)
my_pdlaed0.o : my_pdlaed0.f
	$(FF) -c my_pdlaed0.f $(CCFLAG) $(INCFLAG)
my_pdsxedc.o : my_pdsxedc.f
	$(FF) -c my_pdsxedc.f $(CCFLAG) $(INCFLAG)
my_pdlaed1.o : my_pdlaed1.f
	$(FF) -c my_pdlaed1.f $(CCFLAG) $(INCFLAG)
pdlaed2.o : pdlaed2.f
	$(FF) -c pdlaed2.f $(CCFLAG) $(INCFLAG)
pdlaed3.o : pdlaed3.f
	$(FF) -c pdlaed3.f $(CCFLAG) $(INCFLAG)



