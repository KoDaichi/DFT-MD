MPIR_HOME = /usr/local/mpi
#INSTALL_DIR = $(HOME)/lib
INSTALL_DIR =


# select an appropriate make_inc file
#include make_inc.gcc
#include make_inc.Intel
#include make_inc.BX900
include make_inc.K_FX10


OBJS_TRD  = \
	eigen_prd.o \
	eigen_prd_t1.o eigen_prd_t2.o eigen_prd_t4x.o eigen_prd_t5.o \
	eigen_prd_t6_3.o eigen_prd_t7.o eigen_prd_t8.o
OBJS_DCX  = \
	dcx.o \
	my_pdsxedc.o my_pdlaedz.o my_pdlaed0.o my_pdlaed1.o \
	my_pdlaed2.o my_pdlaed3.o my_pdlasrt.o lapack_eigen.o \
	dlaed3.o dlaed6.o dlaed6_init.o \
	dc_redist1.o dc_redist2.o
OBJS_TBK  = \
	trbakwy4.o trbakwy4_body.o
OBJS_ABC  = \
	eigen_sx.o \
	comm.o machine_epsilon.o matrix_adjustment.o \
	CSTAB.o CSTAB_mem.o
OBJS_MOD  = \
	eigen_libs.o eigen_blacs.o eigen_devel.o \
	eigen_house.o eigen_dc.o eigen_trbak.o \
	dlaed6_common.o

OBJS      = $(OBJS_TRD) $(OBJS_DCX) $(OBJS_TBK) $(OBJS_ABC) $(OBJS_MOD)

OBJS_TEST = main2.o mat_set.o ev_test.o

MODEIGEN_PUB  = \
	eigen_blacs.mod eigen_libs.mod
MODEIGEN_SYS  = \
	eigen_devel.mod \
	eigen_house.mod eigen_dc.mod eigen_trbak.mod \
	dlaed6_common.mod

MODEIGEN  = $(MODEIGEN_PUB) $(MODEIGEN_SYS)

LIBEIGEN  = libEigenExa.a


# MACROS
AT_SYNC   = -DUSE_SYNC=1 -D_DEBUG_=0 -D_DEBUG2_=0
AT_BCAST_OVERLAP = -DAT_BCAST_OVERLAP=1 -DOVERLAP_DECISION_TYPE=3
DEBUGFLAG = -DTIMER_PRINT
#DEBUGFLAG =

all : $(LIBEIGEN) $(MODEIGEN)
test : a.out


.o.mod:
	@true
.F.o:
	$(FF) -c $< $(CCFLAG) $(INCFLAG) $(MODFLAG) $(CPPFLAG) $(DEBUGFLAG)
.f.o:
	$(FF) -c $< $(CCFLAG) $(INCFLAG) $(MODFLAG) $(DEBUGFLAG)
.c.o:
	$(CC) -c $< $(CCFLAG) $(INCFLAG) $(DEBUGFLAG)


$(LIBEIGEN) : $(OBJS)
	$(AR) r $(LIBEIGEN) $(OBJS)

a.out : $(OBJS_TEST) $(LIBEIGEN) 
	$(FF) -o $@  $(OBJS_TEST) -L./ -lEigenExa $(LDFLAG) -lEigenExa $(LDFLAG)


eigen_prd.o : eigen_prd.F 
	$(FF) -c $< $(CCFLAG) $(INCFLAG) $(MODFLAG) $(CPPFLAG) $(AT_SYNC) $(DEBUGFLAG)
trbakwy4.o : trbakwy4.F 
	$(FF) -c $< $(CCFLAG) $(INCFLAG) $(MODFLAG) $(CPPFLAG) $(AT_BCAST_OVERLAP) $(DEBUGFLAG)
trbakwy4_body.o : trbakwy4_body.F
	$(FF) -c $< $(CCFLAG) $(INCFLAG) $(MODFLAG) $(CPPFLAG) $(AT_BCAST_OVERLAP) $(DEBUGFLAG)


$(OBJS_TRD) $(OBJS_DCX) $(OBJS_TBK) $(OBJS_ABC) : $(MODEIGEN)
$(OBJS_TEST) : $(MODEIGEN)
eigen_libs.mod : eigen_blacs.mod eigen_devel.mod 


install : all
	@if [ '$(INSTALL_DIR)X' == 'X' ]; then \
		echo 'Please setup INSTALL_DIR in Makefile'; \
	else \
		if [ ! -d $(INSTALL_DIR) ]; then \
			mkdir $(INSTALL_DIR); \
		fi; \
		for f in $(LIBEIGEN) $(MODEIGEN_PUB) ; do \
			echo 'install ./'$$f' $(INSTALL_DIR)/' ; \
			install ./$$f $(INSTALL_DIR)/; \
		done; \
	fi

clean :
	@echo 'rm *.o *.mod '$(LIBEIGEN)' a.out'
	@for f in *.o *.mod $(LIBEIGEN) a.out ; do \
		if [ -e $$f ]; then \
		if [ -f $$f ]; then \
			\rm $$f; \
		fi; \
		fi; \
	done

