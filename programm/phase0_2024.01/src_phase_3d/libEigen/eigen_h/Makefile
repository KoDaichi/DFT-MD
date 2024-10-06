FF = mpifrtpx
FC = mpifrtpx
CC = mpifccpx

#CCFLAG    = -Kfast,parallel,ocl,preex,array_private,auto,simd=2,openmp -X9 -Cpp
CCFLAG    = -Kfast,parallel,ocl,array_private,auto,openmp,preex,simd=2 -X9 -Cpp

SCALAPACK = -SSL2BLAMP -SCALAPACK
BLAS = 

INCPATH   = 
INCFLAG   = $(INCPATH)
LDFLAG    = $(CCFLAG) $(SCALAPACK) $(BLAS)


OBJS = communication_h.o eigen_hermite.o eigen_dch.o machine_epsilon.o matrix_adjust_h.o dlaed1x.o dlaed3x.o dlaed4x.o dlaed6x.o pdlaed3x.o pdlaed2x.o eigen_hrd.o eigen_hrd_2update.o eigen_hrd_Au.o eigen_hrd_v.o eigen_hrd_u.o eigen_hbk.o

OBJSM = ev_test2.o matrix_set.o main.o

LIB   = Eigen_h
AR    = ar
RANLIB= ranlib
RM    = -rm -f

all : lib$(LIB).a

lib$(LIB).a: $(OBJS)
	$(RM) $@
	$(AR) r $@ $(OBJS)
	$(RANLIB) $@

a.out : main.o $(OBJSM)
	$(FF) -o a.out $(LDFLAG) $(OBJSM) -L.. -lEigen

main.o : main.f
	$(FF) -c main.f $(CCFLAG) -DCHECK $(INCFLAG)
eigen_hermite.o : eigen_hermite.f
	$(FF) -c eigen_hermite.f $(CCFLAG) $(INCFLAG)
eigen_hrd.o : eigen_hrd.f
	$(FF) -c eigen_hrd.f $(CCFLAG) $(INCFLAG)
eigen_hrd_2update.o : eigen_hrd_2update.f
	$(FF) -c eigen_hrd_2update.f $(CCFLAG) $(INCFLAG)
eigen_hrd_Au.o : eigen_hrd_Au.f
	$(FF) -c eigen_hrd_Au.f $(CCFLAG) $(INCFLAG)
eigen_hrd_v.o : eigen_hrd_v.f
	$(FF) -c eigen_hrd_v.f $(CCFLAG) $(INCFLAG)
eigen_hrd_u.o : eigen_hrd_u.f
	$(FF) -c eigen_hrd_u.f $(CCFLAG) $(INCFLAG)
eigen_hbk.o : eigen_hbk.f
	$(FF) -c eigen_hbk.f $(CCFLAG) $(INCFLAG)
machine_epsilon.o : machine_epsilon.f
	$(FF) -c machine_epsilon.f -O0 $(INCFLAG)
matrix_set.o: matrix_set.f
	$(FF) -c matrix_set.f $(CCFLAG) $(INCFLAG)
matrix_adjust_h.o: matrix_adjust_h.f
	$(FF) -c matrix_adjust_h.f $(CCFLAG) $(INCFLAG)
ev_test2.o: ev_test2.f
	$(FF) -c ev_test2.f $(CCFLAG) $(INCFLAG)
communication_h.o : communication_h.f
	$(FF) -c communication_h.f $(CCFLAG) $(INCFLAG)
eigen_dch.o : eigen_dch.f
	$(FF) -c eigen_dch.f $(CCFLAG) $(INCFLAG)
dlaed1x.o : dlaed1x.f
	$(FF) -c dlaed1x.f $(CCFLAG) $(INCFLAG)
dlaed3x.o : dlaed3x.f
	$(FF) -c dlaed3x.f $(CCFLAG) $(INCFLAG)
dlaed4x.o : dlaed4x.f
	$(FF) -c dlaed4x.f $(CCFLAG) $(INCFLAG)
dlaed6x.o : dlaed6x.f
	$(FF) -c dlaed6x.f $(CCFLAG) $(INCFLAG)
pdlaed2x.o : pdlaed2x.f
	$(FF) -c pdlaed2x.f $(CCFLAG) $(INCFLAG)
pdlaed3x.o : pdlaed3x.f
	$(FF) -c pdlaed3x.f $(CCFLAG) $(INCFLAG)

clean :
	$(RM) a.out lib$(LIB).a
	$(RM) *.o *.mod

