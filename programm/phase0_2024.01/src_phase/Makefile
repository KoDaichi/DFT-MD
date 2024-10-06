.SUFFIXES:
.SUFFIXES: .o .F .f .F90 .f90 .c .mod

# Platform    : GNU Linux (EM64T/AMD64)
# Prog. model : Serial
# Compiler    : Intel Fortran compiler
# BLAS/LAPACK : Netlib BLAS/LAPACK
# FFT         : Built-in FFT subroutnes
###########################################################################
###<< PLEASE CHANGE THE VARIABLES BELOW ACCORDING TO YOUR ENVIRONMENT >>###
###########################################################################
F90 = ifort
CC  = gcc -m64
CPP = 
AR  = ar -vq
LINK = ifort
F90FLAGS = -traceback -extend-source -I./no_mpi/include
F77FLAGS = -traceback -extend-source -I./no_mpi/include
CFLAGS = -O -DINTEL

ESM = 
ifdef ESM
CPPESM=-DENABLE_ESM_PACK
LESM=-lesm
else
CPPESM=
endif
NO_MPI=YES
CPPFLAGS = -DLinux -DJRCATFFT_WS -DCD_JRCATFFT_WS -D_NO_MPI_ -D_NO_ARG_DUMMY_ -DDISABLE_VDWDF -D_POT_SMOOTHING_ -DTRANSPOSE -DGGA_ATOMIC_WITH_NEW_GNCPP -DREMOVE_PC_FROM_FORCE -D_POSITRON_ -D_FAST_WAY_ -D_USE_DATE_AND_TIME_ -D__EDA__ -DLIBRARY_BUILD -DUSE_NONBLK_COMM -DRMM_NONLOCAL_NEW ${CPPESM}
LFLAGS = 
F90FLAGS_FIXED = -extend_source -Fl -fixed
F90FLAGS_FREE = -extend_source -Fl
MKLHOME=./
INCLUDE=-I
LIBS =  -L./ -llapack -lblas  
LIBXC=

ifdef LIBXC
  LIBXC_PATH=
  LIBXC_INC=-I$(LIBXC_PATH)/include
  LIBXC_LIB=-L$(LIBXC_PATH)/lib -lxcf03 -lxc
  CPPFLAGS += -DLIBXC
  F90FLAGS += $(LIBXC_INC)
  F77FLAGS += $(LIBXC_INC)
  LIBS += $(LIBXC_LIB)
endif

###########################################################################
###<< PLEASE CHANGE THE VARIABLES ABOVE ACCORDING TO YOUR ENVIRONMENT >>###
###########################################################################

ifdef ESM
ESM_LIB = libesm.a
else
ESM_LIB =
endif
LAPACK = liblapack.a libblas.a
FFTOBJECT = jrcat_fft.o
OBJ_INPUTPARSE = input_parse.o


PHASE_LOWER_MODULES = mpi_dummy.o m_Const_Parameters.o m_ErrorMessages.o m_Parallelization.o \
m_IterationNumbers.o

PHASE_UPPER_MODULES = \
m_Spg_plus_Tetra_Common_Blocks.o \
m_Control_Parameters.o   m_Files.o                m_Timing.o \
m_Crystal_Structure.o    m_FFT.o                  m_Ionic_System.o \
m_CS_SpaceGroup.o        m_CS_Magnetic.o \
m_Orbital_QuantumNum.o  \
m_Kpoints.o              m_PseudoPotential.o      m_PlaneWaveBasisSet.o \
m_Realspace.o \
m_SpinOrbit_Potential.o  m_SpinOrbit_FromFile.o  \
m_ES_NonCollinear.o \
m_NonLocal_Potential.o   m_Electronic_Structure.o \
m_ES_nonlocal.o	 	 m_ES_ortho.o  m_ES_wf_extrpl.o m_ES_initialWF.o 	\
m_ES_occup.o \
m_ES_RSB.o \
m_FiniteElectricField.o \
m_ES_ExactExchange.o \
m_ES_WF_by_SDorCG.o      m_ES_WF_by_Davidson.o \
m_ES_WF_by_ModifiedDavidson.o m_ES_WF_by_RMM.o \
m_Charge_Density.o  m_CD_Mag_Moment.o m_Orbital_Population.o      m_Dipole.o  m_CD_mixing.o \
m_PS_opencore.o \
m_epc_potential.o        \
m_Positron_Wave_Functions.o \
m_ES_WF_mixing.o \
m_Screening_FFT.o        m_Screening.o      m_External_Potential.o \
m_ES_LHXC.o              m_ES_Intgr_VlhxcQlm.o    m_ES_IO.o \
m_Representation.o       m_OP_rotation.o          m_OP_decomp_band.o \
m_ES_WF_by_MatDiagon.o   m_ES_dos.o               m_Hubbard.o \
m_KineticEnergy_Density.o  \
m_KE_mixing.o \
m_vdWDF.o \
m_Ldos.o                 m_XC_Potential.o         \
string.o m_db.o \
m_PAW_Tecplot.o m_PAW_ChargeDensity.o m_PAW_Hartree.o  \
m_PAW_XC_Potential.o \
m_SpinOrbit_RadInt.o \
m_OP_Moment.o \
m_ES_Mag_Constraint.o  \
m_Total_Energy.o \
m_Force.o                m_Stress.o               m_ES_WF_by_submat.o \
m_UnitCell.o \
m_constraints.o \
m_ELF.o \
m_BerryPhase.o           m_BP_Properties.o \
m_Raman.o            m_Phonon.o \
m_Wannier.o m_Wannier90.o m_Replica.o Renewal_of_ChgCtrlParam.o \
m_dimer.o \
m_LinearResponse_Control.o  \
m_ValenceBand_Spectrum.o   m_CoreLevel_Spectrum.o   \
m_Excitation.o \
m_CLS_dipquad.o m_ES_occup_EPS.o         m_Epsilon_ek.o \
m_LinearResponse_Qpt.o  \
m_ES_ChargeState.o \
m_rttddft.o \
m_Potential_Mixing.o       m_ThomasFermiW_Potential.o  \
m_SpinOrbit_SecondVariation.o m_Crystal_Field.o \
m_PlaneWaveEDA.o \
m_Fcp.o m_Band_Unfolding.o  m_Potential_Average.o

PHASE_F_SUBROUTINES = Preparation.o \
InputData_Analysis.o           PseudoPotential_Construction.o \

PHASE_F_OTHERSUBS  = mdmain0.o    constraint_main.o meta_dynamics.o  NEB.o Preparation_for_mpi.o \
Preparation_for_ESM.o \
scf_routines.o \
UnitaryTransform_WF.o \
bottom_Subroutines.o              \
spline.o \
b_Crystal_Structure.o             b_Electronic_Structure.o \
b_Words.o                         b_PseudoPotential.o \
b_PseudoPotential_EXX.o \
b_Fermi.o                         b_Kpoints.o \
b_PlaneWaveBasisSet.o             b_Ionic_System.o \
b_XC_Potential.o                  b_PAW_XC_Potential.o \
b_XC_OmegaPBE.o  \
b_XC_metagga.o \
b_Ldos_f77.o\
input_interface.o                 Initialization.o \
WriteDownData_onto_Files.o        Ewald_and_Structure_Factor.o \
Initial_Electronic_Structure.o \
Renewal_of_WaveFunctions.o \
Renewal_of_pWaveFunctions.o \
Uramping.o \
IterationNumbers_Setting.o        ChargeDensity_Construction.o \
ChargeDensity_Mixing.o            Renewal_of_Potential.o \
Renewal_of_Hubbard_Potential.o    Renewal_of_OccMat.o \
Renewal_of_pPotential.o           Renewal_of_Hubbard_Parameters.o \
Convergence_Check.o               Forces.o \
Move_Ions.o                       Initial_MD_Condition.o \
Stress.o                          Postprocessing.o \
Finalization_of_mpi.o             miscellaneous.o\
b_BerryPhase.o \
heap_sort.o                       real_spherical_harmonics.o \
Real_space_integ.o                crotylm.o \
screening_correction.o \
Initialization_Epsilon.o               Shift_Kpoint.o \
Reset_Kpoint.o                         Preparation_for_Calc_Epsilon.o \
Transition_moment_Epsilon.o            Calc_Epsilon.o \
Nonlinear_Optics_Epsilon.o             WriteDownData_onto_Files_Epsilon.o \
PseudoPotential_ek_Epsilon.o           Dealloc_Radr_and_Wos_Epsilon.o \
                                       WriteDownData_onto_Files_ek.o \
GaussLeg.o lib_int_deri_add.o  \
rttddft_main.o \
Potential_Construction.o        Potential_Mixing.o             ThomasFermiWeiz.o \
Epsilon_postscf.o vdW.o Epsilon_Paramset.o EPS_Prep_Check.o memsize.o scdft_routines.o

PHASE_INTERFACE = phase0_interface.o

# for vc_nl
NLOBJ = vc_nl.o

#
ifdef LIBXC
EXT_XC_OBJ = b_XC_libxc.o
else
EXT_XC_OBJ =
endif

ifndef SX_DGEMM
PHASE_OBJECTSF77 = \
b_PseudoPotential_f77.o \
b_Force_f77.o \
rmmsubs.o \
gncpp_xc_gga_rad.o  \
decfft_ent.o \
spg+tetra.o
else
PHASE_OBJECTSF77 = \
b_PseudoPotential_f77.o \
b_Force_f77.o \
rmmsubs.o \
decfft_ent.o \
gncpp_xc_gga_rad.o  \
spg+tetra.o \
dgemm__.o 
endif

EKCAL_LOWER_MODULES = m_Const_Parameters.o   m_ErrorMessages.o m_Parallelization.o \
m_IterationNumbers.o

EKCAL_UPPER_MODULES = \
m_Spg_plus_Tetra_Common_Blocks.o \
m_Control_Parameters.o   m_Files.o                m_Timing.o \
m_Crystal_Structure.o    m_FFT.o                  m_Ionic_System.o \
m_CS_SpaceGroup.o        m_CS_Magnetic.o \
m_Orbital_QuantumNum.o  \
m_Kpoints.o              m_PseudoPotential.o      m_PlaneWaveBasisSet.o \
m_Realspace.o \
m_SpinOrbit_Potential.o  m_SpinOrbit_FromFile.o  \
m_ES_NonCollinear.o \
m_NonLocal_Potential.o   m_Electronic_Structure.o \
m_ES_nonlocal.o	 	 m_ES_ortho.o  m_ES_wf_extrpl.o m_ES_initialWF.o 	\
m_FiniteElectricField.o \
m_ES_ExactExchange.o \
m_ES_occup.o             m_ES_WF_by_SDorCG.o      m_ES_WF_by_Davidson.o \
m_ES_RSB.o \
m_ES_WF_by_ModifiedDavidson.o m_ES_WF_by_RMM.o \
m_epc_potential.o \
m_Positron_Wave_Functions.o \
m_Charge_Density.o m_CD_Mag_Moment.o m_Orbital_Population.o      m_Dipole.o  m_CD_mixing.o \
m_PS_opencore.o \
m_ES_WF_mixing.o \
m_Screening_FFT.o        m_Screening.o      m_External_Potential.o \
m_ES_LHXC.o              m_ES_Intgr_VlhxcQlm.o    m_ES_IO.o \
m_Representation.o       m_OP_rotation.o          m_OP_decomp_band.o \
m_ES_WF_by_MatDiagon.o   m_ES_dos.o               m_Hubbard.o \
m_KineticEnergy_Density.o  \
m_KE_mixing.o \
m_vdWDF.o \
m_Ldos.o                 m_XC_Potential.o         \
string.o m_db.o \
m_PAW_Tecplot.o m_PAW_ChargeDensity.o m_PAW_Hartree.o  \
m_PAW_XC_Potential.o \
m_SpinOrbit_RadInt.o \
m_ES_Mag_Constraint.o  \
m_constraints.o \
m_Total_Energy.o \
m_OP_Moment.o \
m_Force.o                m_Stress.o m_ES_WF_by_submat.o\
m_UnitCell.o \
m_ELF.o                  m_Wannier.o \
m_Wannier90.o \
m_BerryPhase.o           m_BP_Properties.o \
m_Raman.o               m_Phonon.o \
m_LinearResponse_Control.o  \
m_ValenceBand_Spectrum.o   m_CoreLevel_Spectrum.o \
m_Excitation.o \
m_ES_occup_EPS.o         m_Epsilon_ek.o \
m_LinearResponse_Qpt.o \
m_ES_ChargeState.o \
m_Potential_Mixing.o       m_ThomasFermiW_Potential.o \
m_SpinOrbit_SecondVariation.o m_Crystal_Field.o \
m_PlaneWaveEDA.o \
m_Fcp.o m_CLS_dipquad.o m_Band_Unfolding.o m_Potential_Average.o

EKCAL_F_SUBROUTINES = Preparation.o \
InputData_Analysis.o           PseudoPotential_Construction.o

EKCAL_F_OTHERSUBS  = ekmain.o     constraint_main.o          Preparation_for_mpi.o \
Preparation_for_ESM.o \
scf_routines.o \
UnitaryTransform_WF.o \
bottom_Subroutines.o                   \
spline.o \
b_Crystal_Structure.o                  b_Electronic_Structure.o \
b_Words.o                              b_PseudoPotential.o \
b_PseudoPotential_EXX.o \
b_Fermi.o                              b_Kpoints.o \
b_PlaneWaveBasisSet.o                  b_Ionic_System.o \
b_XC_Potential.o                       b_PAW_XC_Potential.o \
b_XC_OmegaPBE.o  \
b_XC_metagga.o \
b_Ldos_f77.o \
input_interface.o                      Initialization.o \
Ewald_and_Structure_Factor.o \
Initial_Electronic_Structure.o         Renewal_of_WaveFunctions.o \
IterationNumbers_Setting.o \
ChargeDensity_Construction.o                 Renewal_of_Potential.o \
ChargeDensity_Mixing.o \
Renewal_of_Hubbard_Potential.o         Renewal_of_OccMat.o \
Convergence_Check.o                    Forces.o \
WriteDownData_onto_Files.o \
Initial_MD_Condition.o                 WriteDownData_onto_Files_ek.o\
Stress.o                               Postprocessing.o \
Finalization_of_mpi.o \
miscellaneous.o \
b_BerryPhase.o \
heap_sort.o                            real_spherical_harmonics.o \
Real_space_integ.o                     crotylm.o \
screening_correction.o \
mpi_dummy.o \
GaussLeg.o lib_int_deri_add.o \
Uramping.o \
ThomasFermiWeiz.o       Renewal_of_ChgCtrlParam.o vdW.o scdft_routines.o memsize.o

ifndef SX_DGEMM
EKCAL_OBJECTSF77 = \
b_PseudoPotential_f77.o \
b_Force_f77.o \
rmmsubs.o                              spg+tetra.o \
gncpp_xc_gga_rad.o  \
decfft_ent.o
else
EKCAL_OBJECTSF77 = \
b_PseudoPotential_f77.o \
b_Force_f77.o \
rmmsubs.o                              spg+tetra.o \
decfft_ent.o \
gncpp_xc_gga_rad.o  \
dgemm__.o 
endif


EPS_LOWER_MODULES = m_Const_Parameters.o   m_ErrorMessages.o m_Parallelization.o \
m_IterationNumbers.o

EPS_UPPER_MODULES = \
m_Spg_plus_Tetra_Common_Blocks.o \
m_Control_Parameters.o   m_Files.o                m_Timing.o \
m_Crystal_Structure.o    m_FFT.o                  m_Ionic_System.o \
m_CS_SpaceGroup.o        m_CS_Magnetic.o \
m_Orbital_QuantumNum.o  \
m_Kpoints.o              m_PseudoPotential.o      m_PlaneWaveBasisSet.o \
m_Realspace.o \
m_SpinOrbit_Potential.o  m_SpinOrbit_FromFile.o  \
m_ES_NonCollinear.o \
m_NonLocal_Potential.o   m_Electronic_Structure.o \
m_ES_nonlocal.o	 	 m_ES_ortho.o  m_ES_wf_extrpl.o m_ES_initialWF.o 	\
m_FiniteElectricField.o \
m_Positron_Wave_Functions.o \
m_ES_ExactExchange.o \
m_ES_occup.o             m_ES_WF_by_SDorCG.o      m_ES_WF_by_Davidson.o \
m_ES_RSB.o \
m_ES_WF_by_ModifiedDavidson.o m_ES_WF_by_RMM.o \
m_ES_LHXC.o              m_ES_Intgr_VlhxcQlm.o    m_ES_IO.o \
m_Representation.o       m_OP_rotation.o          m_OP_decomp_band.o \
m_ES_dos.o               m_Hubbard.o \
m_epc_potential.o \
m_vdWDF.o \
m_Charge_Density.o m_CD_Mag_Moment.o  m_Orbital_Population.o     m_Dipole.o m_CD_mixing.o \
m_XC_Potential.o \
m_PS_opencore.o \
m_ES_WF_mixing.o \
string.o m_db.o \
m_PAW_Tecplot.o m_PAW_ChargeDensity.o m_PAW_Hartree.o  \
m_PAW_XC_Potential.o \
m_SpinOrbit_RadInt.o \
m_ES_Mag_Constraint.o  \
m_Total_Energy.o \
m_OP_Moment.o \
m_Force.o          m_Stress.o      m_ES_WF_by_submat.o      m_ES_WF_by_MatDiagon.o \
m_UnitCell.o \
m_KineticEnergy_Density.o  \
m_KE_mixing.o \
m_Ldos.o                 m_ELF.o \
m_constraints.o \
m_BerryPhase.o           m_BP_Properties.o \
m_Raman.o       m_Phonon.o               \
m_Screening_FFT.o        m_Screening.o      m_External_Potential.o \
m_Wannier.o \
m_Wannier90.o \
m_ValenceBand_Spectrum.o   m_CoreLevel_Spectrum.o   \
m_Excitation.o \
m_CLS_dipquad.o m_ES_occup_EPS.o         m_Epsilon_ek.o \
m_LinearResponse_Control.o  \
m_LinearResponse_Qpt.o \
m_ES_ChargeState.o \
m_Potential_Mixing.o       m_ThomasFermiW_Potential.o \
m_SpinOrbit_SecondVariation.o m_Crystal_Field.o \
m_PlaneWaveEDA.o \
m_Fcp.o m_Band_Unfolding.o m_Potential_Average.o

EPS_F_SUBROUTINES = Preparation.o \
InputData_Analysis.o           PseudoPotential_Construction.o

EPS_F_OTHERSUBS  = epsmain.o      constraint_main.o     Preparation_for_mpi.o \
Preparation_for_ESM.o \
scf_routines.o \
UnitaryTransform_WF.o \
bottom_Subroutines.o                   \
spline.o \
b_Crystal_Structure.o                  b_Electronic_Structure.o \
b_Words.o                              b_PseudoPotential.o \
b_PseudoPotential_EXX.o \
b_Fermi.o                              b_Kpoints.o \
b_PlaneWaveBasisSet.o                  b_Ionic_System.o \
b_XC_Potential.o                       b_PAW_XC_Potential.o \
b_XC_OmegaPBE.o  \
b_XC_metagga.o \
b_Ldos_f77.o \
input_interface.o                      Initialization.o \
Ewald_and_Structure_Factor.o \
Initial_Electronic_Structure.o         Renewal_of_WaveFunctions.o \
IterationNumbers_Setting.o \
ChargeDensity_Construction.o                 Renewal_of_Potential.o \
ChargeDensity_Mixing.o \
Renewal_of_Hubbard_Potential.o         Renewal_of_OccMat.o \
Convergence_Check.o                    Forces.o \
WriteDownData_onto_Files.o \
Initial_MD_Condition.o                 WriteDownData_onto_Files_ek.o\
Stress.o                               Postprocessing.o \
Finalization_of_mpi.o \
miscellaneous.o \
b_BerryPhase.o \
heap_sort.o                            real_spherical_harmonics.o \
Real_space_integ.o                     crotylm.o \
screening_correction.o \
Initialization_Epsilon.o               Shift_Kpoint.o \
Reset_Kpoint.o                         Preparation_for_Calc_Epsilon.o \
Transition_moment_Epsilon.o            Calc_Epsilon.o \
Nonlinear_Optics_Epsilon.o             WriteDownData_onto_Files_Epsilon.o \
PseudoPotential_ek_Epsilon.o           Dealloc_Radr_and_Wos_Epsilon.o \
EPS_Prep_Check.o \
ThomasFermiWeiz.o         Renewal_of_ChgCtrlParam.o  \
mpi_dummy.o \
Uramping.o \
GaussLeg.o lib_int_deri_add.o vdW.o memsize.o scdft_routines.o

ifndef SX_DGEMM
EPS_OBJECTSF77 = \
b_PseudoPotential_f77.o \
b_Force_f77.o \
rmmsubs.o                              spg+tetra.o \
decfft_ent.o \
gncpp_xc_gga_rad.o  
else
EPS_OBJECTSF77 = \
b_PseudoPotential_f77.o \
b_Force_f77.o \
rmmsubs.o                              spg+tetra.o \
decfft_ent.o \
gncpp_xc_gga_rad.o  \
dgemm__.o
endif


TDLR_LOWER_MODULES = m_Const_Parameters.o   m_ErrorMessages.o m_Parallelization.o \
m_IterationNumbers.o

TDLR_UPPER_MODULES = \
m_Spg_plus_Tetra_Common_Blocks.o \
m_Control_Parameters.o   m_Files.o                m_Timing.o \
m_Crystal_Structure.o    m_FFT.o                  m_Ionic_System.o \
m_CS_SpaceGroup.o        m_CS_Magnetic.o \
m_Orbital_QuantumNum.o  \
m_Kpoints.o              m_PseudoPotential.o      m_PlaneWaveBasisSet.o \
m_Realspace.o \
m_SpinOrbit_Potential.o  m_SpinOrbit_FromFile.o  \
m_ES_NonCollinear.o \
m_NonLocal_Potential.o   m_Electronic_Structure.o \
m_ES_nonlocal.o	 	 m_ES_ortho.o  m_ES_wf_extrpl.o m_ES_initialWF.o 	\
m_FiniteElectricField.o \
m_Positron_Wave_Functions.o \
m_ES_ExactExchange.o \
m_ES_occup.o             m_ES_WF_by_SDorCG.o      m_ES_WF_by_Davidson.o \
m_ES_RSB.o \
m_ES_WF_by_ModifiedDavidson.o m_ES_WF_by_RMM.o \
m_ES_LHXC.o              m_ES_Intgr_VlhxcQlm.o    m_ES_IO.o \
m_Representation.o       m_OP_rotation.o          m_OP_decomp_band.o \
m_ES_dos.o               m_Hubbard.o \
m_epc_potential.o \
m_vdWDF.o \
m_Charge_Density.o m_CD_Mag_Moment.o m_Orbital_Population.o m_Dipole.o m_CD_mixing.o \
m_XC_Potential.o \
m_PS_opencore.o \
m_ES_WF_mixing.o \
string.o m_db.o \
m_PAW_Tecplot.o m_PAW_ChargeDensity.o m_PAW_Hartree.o  \
m_PAW_XC_Potential.o \
m_SpinOrbit_RadInt.o \
m_ES_Mag_Constraint.o  \
m_Total_Energy.o \
m_OP_Moment.o \
m_Force.o                m_Stress.o m_ES_WF_by_submat.o      m_ES_WF_by_MatDiagon.o \
m_UnitCell.o \
m_KineticEnergy_Density.o  \
m_KE_mixing.o \
m_Ldos.o                 m_ELF.o \
m_constraints.o \
m_BerryPhase.o           m_BP_Properties.o \
m_Raman.o          m_Phonon.o               \
m_Wannier.o \
m_Wannier90.o \
m_Screening_FFT.o        m_Screening.o      m_External_Potential.o \
m_Excitation.o \
m_CLS_dipquad.o m_ES_occup_EPS.o m_Epsilon_ek.o       \
m_LinearResponse_Control.o  \
m_ValenceBand_Spectrum.o  m_CoreLevel_Spectrum.o   \
m_Band_Unfolding.o \
m_ES_ChargeState.o \
m_LinearResponse_Qpt.o \
m_LinearResponse_Tools.o     m_LinearResponse_Density.o  \
m_LinearResponse_NonInt.o    m_LinearResponse_ALDA.o \
m_LinearResponse_Kernel.o    m_LinearResponse_BS.o \
m_LinearResponse_Spectrum.o \
m_Potential_Mixing.o       m_ThomasFermiW_Potential.o \
m_SpinOrbit_SecondVariation.o m_Crystal_Field.o \
m_PlaneWaveEDA.o \
m_Potential_Average.o

TDLR_F_SUBROUTINES = Preparation.o \
InputData_Analysis.o           PseudoPotential_Construction.o \

TDLR_F_OTHERSUBS  = tdlrmain.o     constraint_main.o    Preparation_for_mpi.o \
Preparation_for_ESM.o \
scf_routines.o \
UnitaryTransform_WF.o \
bottom_Subroutines.o                   \
spline.o \
b_Crystal_Structure.o                  b_Electronic_Structure.o \
b_Words.o                              b_PseudoPotential.o \
b_PseudoPotential_EXX.o \
b_Fermi.o                              b_Kpoints.o \
b_PlaneWaveBasisSet.o                  b_Ionic_System.o \
b_XC_Potential.o                       b_PAW_XC_Potential.o \
b_XC_OmegaPBE.o  \
b_XC_metagga.o \
b_Ldos_f77.o \
input_interface.o                      Initialization.o \
Ewald_and_Structure_Factor.o \
Initial_Electronic_Structure.o         Renewal_of_WaveFunctions.o \
IterationNumbers_Setting.o \
ChargeDensity_Construction.o                 Renewal_of_Potential.o \
ChargeDensity_Mixing.o \
Renewal_of_Hubbard_Potential.o         Renewal_of_OccMat.o \
Convergence_Check.o                    Forces.o \
WriteDownData_onto_Files.o \
Initial_MD_Condition.o                 WriteDownData_onto_Files_ek.o \
Stress.o                               Postprocessing.o \
Finalization_of_mpi.o \
miscellaneous.o \
b_BerryPhase.o \
heap_sort.o                            real_spherical_harmonics.o \
Real_space_integ.o                     crotylm.o \
screening_correction.o \
PseudoPotential_ek_Epsilon.o           \
mpi_dummy.o \
GaussLeg.o lib_int_deri_add.o \
ThomasFermiWeiz.o       Renewal_of_ChgCtrlParam.o \
b_LinearResponse_Kernel.o  b_LinearResponse_exc.o \
Uramping.o \
LinearResponse_Proc.o     LinearResponse_Spec.o vdW.o scdft_routines.o memsize.o

ifndef SX_DGEMM
TDLR_OBJECTSF77 = \
b_PseudoPotential_f77.o \
b_Force_f77.o \
rmmsubs.o                              spg+tetra.o \
decfft_ent.o \
gncpp_xc_gga_rad.o  
else
TDLR_OBJECTSF77 = \
b_PseudoPotential_f77.o \
b_Force_f77.o \
rmmsubs.o                              spg+tetra.o \
decfft_ent.o \
gncpp_xc_gga_rad.o  \
dgemm__.o
endif

PHASE_OBJECTS = $(PHASE_LOWER_MODULES) $(FFTOBJECT) $(PHASE_UPPER_MODULES) $(PHASE_F_SUBROUTINES) $(PHASE_F_OTHERSUBS) $(PHASE_OBJECTSF77) $(OBJ_INPUTPARSE) $(EXT_XC_OBJ)

EKCAL_OBJECTS = $(EKCAL_LOWER_MODULES) $(FFTOBJECT) $(EKCAL_UPPER_MODULES) $(EKCAL_F_SUBROUTINES) $(EKCAL_F_OTHERSUBS) $(EKCAL_OBJECTSF77) $(OBJ_INPUTPARSE) $(EXT_XC_OBJ)

EPS_OBJECTS = $(EPS_LOWER_MODULES) $(FFTOBJECT) $(EPS_UPPER_MODULES) $(EPS_F_SUBROUTINES) $(EPS_F_OTHERSUBS) $(EPS_OBJECTSF77) $(OBJ_INPUTPARSE) $(EXT_XC_OBJ)

TDLR_OBJECTS = $(TDLR_LOWER_MODULES) $(FFTOBJECT) $(TDLR_UPPER_MODULES) $(TDLR_F_SUBROUTINES) $(TDLR_F_OTHERSUBS) $(TDLR_OBJECTSF77) $(OBJ_INPUTPARSE) $(EXT_XC_OBJ)

PHASE_LIB_OBJECTS = $(PHASE_OBJECTS) $(PHASE_INTERFACE)

all : phase ekcal epsmain tdlrmain

libphase.a : $(ESM_LIB) $(LAPACK) $(PHASE_LIB_OBJECTS)
	$(AR) $@ $?

ifdef ESM
phase : $(ESM_LIB) $(LAPACK) $(PHASE_OBJECTS) $(NLOBJ)
	$(LINK) $(PHASE_OBJECTS) $(NLOBJ) $(LFLAGS) $(LIBS) -o $@
else
phase : $(LAPACK) $(PHASE_OBJECTS) $(NLOBJ)
	$(LINK) $(PHASE_OBJECTS) $(NLOBJ) $(LFLAGS) $(LIBS) -o $@
endif

ekcal : $(LAPACK) $(EKCAL_OBJECTS) $(NLOBJ)
	$(LINK) $(EKCAL_OBJECTS) $(NLOBJ) $(LFLAGS) $(LIBS) -o $@

epsmain : $(LAPACK) $(EPS_OBJECTS) $(NLOBJ)
	$(LINK) $(EPS_OBJECTS) $(NLOBJ) $(LFLAGS) $(LIBS) -o $@

tdlrmain : $(LAPACK) $(TDLR_OBJECTS) $(NLOBJ)
	$(LINK) $(TDLR_OBJECTS) $(NLOBJ) $(LFLAGS) $(LIBS) -o $@

ifdef NO_MPI
libesm.a:
	cd EsmPack; make INCLUDE="$(INCLUDE)" FORTRAN="$(F90)" LIBFLAG="$(LIBS)" MPIFLAG="" AR="$(AR)"
else
libesm.a:
	cd EsmPack; make INCLUDE="$(INCLUDE)" FORTRAN="$(F90)" LIBFLAG="$(LIBS)" MPIFLAG="-D__MPI__" AR="$(AR)"
endif

liblapack.a:
	cd LAPACK; make F77="$(F90)" F77FLAGS="$(F77FLAGS)" AR="$(AR)"

libblas.a:
	cd BLAS; make F77="$(F90)" F77FLAGS="$(F77FLAGS)" AR="$(AR)"

$(OBJ_INPUTPARSE):$(@:.o=.c) $(@:.o=.h)
	$(CC) -c $(CFLAGS) $(@:.o=.c)

.f.o:
	$(F90) -c $(F77FLAGS) $*.f

.f90.o:
	$(F90) -c $(F90FLAGS) $*.f90

.F.o:
	$(F90) -c $(F77FLAGS) $(CPPFLAGS) $*.F

.F90.o:
	$(F90) -c $(F90FLAGS) $(CPPFLAGS) $*.F90

clean:
	\rm -f *.o *.mod *.a *.lib *.L *.list phase ekcal epsmain tdlrmain
	\cd LAPACK; make clean
	\cd BLAS; make clean
	\cd EsmPack; make clean

install: phase ekcal epsmain tdlrmain
	\mv -f phase ../bin/
	\mv -f ekcal ../bin/
	\mv -f epsmain ../bin/
	\mv -f tdlrmain ../bin/
