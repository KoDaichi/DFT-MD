MPIR_HOME = /usr/local/mpi
FF = mpifrtpx
FC = mpifrtpx
CC = mpifccpx

# Optimization flags for the Intel compiler
#CCFLAG_BASE = -g -axW -openmp
#CCFLAG  = $(CCFLAG_BASE) -vec-report3 # -ftrapuv -CU -CB -WB
CCFLAG = -Kfast,parallel,ocl,array_private,auto,openmp,preex,simd=2 -X9 -Cpp -DTIMER -DPHASE

#SCALAPACK = -L$(HOME)/lib -lscalapack -lblacsC -lblacsF77 -lblacsBASE -lmkl_lapack 
#BLAS = -L/opt/intel/Compiler/11.0/074/mkl/lib/em64t -lmkl -lmkl_em64t
SCALAPACK = -SSL2BLAMP -SCALAPACK
BLAS =

INCPATH = -I/usr/local/include

INCFLAG = $(INCPATH)

LIB = Eigen_s
AR     = ar
RANLIB = ranlib
RM     = rm -f

# FOR PAPI AVAILABLE
#LDFLAG  = -g $(LDPATH) -Kopenmp $(SCALAPACK) $(BLAS) $(SCALAPACK) -lpapi
#OBJS = main2-papi.o eigen-papi.o trbakwy4.o trbakwy4_body.o machine_epsilon.o mat_set.o ev_test_2D.o tred1.o t1.o t2.o t3.o t4.o t5.o t5x.o t6.o t7.o t8.o CSTAB.o CSTAB_mem.o comm.o dc.o pdlaed2.o pdlaed3.o dlaed3.o papi-stat.o

LDFLAG  = $(CCFLAG) $(LDPATH) $(SCALAPACK) $(BLAS) $(SCALAPACK)
#LDFLAG  = -g $(LDPATH) -Kopenmp $(SCALAPACK) $(BLAS) $(SCALAPACK)
#LDFLAG  = -g -ftrapuv -CU -CB $(LDPATH) -Kopenmp $(SCALAPACK) $(BLAS) $(SCALAPACK)
OBJS = eigen_real_symmetric.o eigen.o trbakwy4.o trbakwy4_body.o machine_epsilon.o tred1.o t1.o t2.o t3.o t4.o t5.o t5x.o t6.o t7.o t8.o CSTAB.o CSTAB_mem.o comm.o dc.o dc_redist1.o dc_redist2.o pdlaed0.o pdlaed1.o pdlaed2.o pdlaed3.o pdlaedz.o pdlasrt.o dlaed3.o
OBJSM = main2.o ev_test_2D.o mat_set.o 

all : lib$(LIB).a eigen_s

lib$(LIB).a: $(OBJS)
	$(RM) $@
	$(AR) r $@ $(OBJS)
	$(RANLIB) $@

eigen_s : $(OBJSM)
	$(FF) -o eigen_s $(LDFLAG) $(OBJSM) -L. -lEigen_s

main2.o : main2.f
	$(FF) -c main2.f $(CCFLAG) $(INCFLAG)
eigen_real_symmetric.o : eigen_real_symmetric.f
	$(FF) -c eigen_real_symmetric.f $(CCFLAG) $(INCFLAG)
eigen.o : eigen.f
	$(FF) -c eigen.f $(CCFLAG) $(INCFLAG)
tred1.o : tred1.f
	$(FF) -c tred1.f $(CCFLAG) $(INCFLAG)
t1.o : t1.f
	$(FF) -c t1.f $(CCFLAG) $(INCFLAG)
t2.o : t2.f
	$(FC) -c t2.f $(CCFLAG) $(INCFLAG)
t3.o : t3.f
	$(FC) -c t3.f $(CCFLAG) $(INCFLAG)
t4.o : t4.f
	$(FC) -c t4.f $(CCFLAG) $(INCFLAG)
t5.o : t5.f
	$(FC) -c t5.f $(CCFLAG) $(INCFLAG)
t5x.o : t5x.f
	$(FC) -c t5x.f $(CCFLAG) $(INCFLAG)
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
ev_test_2D.o: ev_test_2D.f
	$(FF) -c ev_test_2D.f $(CCFLAG) $(INCFLAG)
comm.o : comm.f
	$(FF) -c comm.f

dc.o : dc.f
	$(FC) -c dc.f -Cpp $(CCFLAG) $(INCFLAG)
dc_redist1.o : dc_redist1.f
	$(FC) -c dc_redist1.f $(CCFLAG) $(INCFLAG)
dc_redist2.o : dc_redist2.f
	$(FC) -c dc_redist2.f $(CCFLAG) $(INCFLAG)
dlaed3.o : dlaed3.f
	$(FC) -c dlaed3.f $(CCFLAG) $(INCFLAG)
pdlaed0.o : pdlaed0.f
	$(FC) -c pdlaed0.f $(CCFLAG) $(INCFLAG)
pdlaed1.o : pdlaed1.f
	$(FC) -c pdlaed1.f $(CCFLAG) $(INCFLAG)
pdlaed2.o : pdlaed2.f
	$(FC) -c pdlaed2.f $(CCFLAG) $(INCFLAG)
pdlaed3.o : pdlaed3.f
	$(FC) -c pdlaed3.f $(CCFLAG) $(INCFLAG)
pdlaed4.o : pdlaed4.f
	$(FC) -c pdlaed4.f $(CCFLAG) $(INCFLAG)
pdlaedz.o : pdlaedz.f
	$(FC) -c pdlaedz.f $(CCFLAG) $(INCFLAG)
pdlasrt.o : pdlasrt.f
	$(FC) -c pdlasrt.f $(CCFLAG) $(INCFLAG)


clean :
	-\rm *.o eigen_s lib$(LIB).a

trbakwy4.o: CSTAB.h trd.h
trbakwy4_body.o: CSTAB.h trd.h
tred1.o: CSTAB.h trd.h
comm.o: trd.h
dc.o: trd.h
ev_test_2D.o: trd.h
main2.o: trd.h
mat_set.o: trd.h
t1.o: trd.h
t2.o: trd.h trd_au.h
t3.o: trd.h
t4.o: trd.h
t5.o: trd.h
t6.o: trd.h
t7.o: trd.h
t8.o: trd.h trd_au.h

