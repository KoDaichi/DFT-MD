!#define _DEBUG_WRITE_DFTU_MPI_PROCESSES_
!!$#define DEBUG_ITERATION_WRITE
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  MODULE:  m_Total_Energy
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
!=======================================================================
!
!     The original version of this set of the computer programs "PHASE"
!  was developed by the members of the Theory Group of Joint Research
!  Center for Atom Technology (JRCAT), based in Tsukuba, in the period
!  1993-2001.
!
!     Since 2002, this set has been tuned and new functions have been
!  added to it as a part of the national project "Frontier Simulation
!  Software for Industrial Science (FSIS)",  which is supported by
!  the IT program of the Ministry of Education, Culture, Sports,
!  Science and Technology (MEXT) of Japan.
!     Since 2006, this program set has been developed as a part of the
!  national project "Revolutionary Simulation Software (RSS21)", which
!  is supported by the next-generation IT program of MEXT of Japan.
!   Since 2013, this program set has been further developed centering on PHASE System
!  Consortium.
!   The activity of development of this program set has been supervised by Takahisa Ohno.
!
!
module m_Total_Energy
!    ( m_TE )
! $Id: m_Total_Energy.F90 633 2020-12-01 05:11:03Z jkoga $
  use m_Charge_Density,       only : chgq_l, chgqo_l, hsr, hsro
  use m_XC_Potential,         only : vxc_l, exc,eex,ecor
  use m_XC_Potential,         only : m_XC_cal_potential
  use m_Electronic_Structure, only : nrvf_ordr, occup_l, eko_l, totch &
      &                            , band_entropy, dhub &
      &                            , m_ES_eekdif_cond &
      &                            , zaj_l, vloc_esm, m_ES_what_is_evdff_now
  use m_PseudoPotential,      only : psc_l,etot1,ilmt,ltp,mtp,epc,dion, ival &
                                    , flg_paw,ipaw &
                                    , dion_hartree,dion_vxc,dion_paw &
                                    , dion_hartree_now, dion_kin_ion &
                                    , flg_symmtry,ia2ia_symmtry_op, epc_paw
  use m_PlaneWaveBasisSet,    only : gr_l,m_pwBS_kinetic_energies,iba,kg1,igfp_l
  use m_FFT,                  only : nfftps
  use m_Crystal_Structure,    only : univol,nopr,op
  use m_Ionic_System,         only : ntyp,natm,natm2,iwei,ityp,zfm3_l &
      &                            , eewald, ihubbard &
      &                            , evdw, ntyp_vdw &
      &                            , num_regions, regions
  use m_Timing,               only : tstatc0_begin, tstatc0_end
!!$  use m_Control_Parameters,   only : af, kimg, ipri, nspin, edelta, neg &
  use m_Control_Parameters,   only : af, kimg, ipri, nspin, neg &
       &                           , num_extra_bands, printable, xctype, oneshot &
       &                           , width, way_of_smearing, icond &
       &                           , m_CtrlP_ntcnvg_incre &
       &                           , m_CtrlP_ntcnvg_clear &
       &                           , m_CtrlP_ntcnvg_reset &
       &                           , m_CtrlP_get_edelta &
       &                           , m_CtrlP_sub_ntcnvg_incre &
       &                           , m_CtrlP_solver_for_WFs_now &
       &                           , m_CtrlP_get_isolver_now &
       &                           , sub_delta_factor_is_given &
       &                           , sub_delta_factor &
       &                           , sw_dipole_correction, sw_screening_correction &
       &                           , sw_hubbard, proj_attribute &
       &                           , critical_ehub, delta_ehub &
       &                           , num_conduction_bands_lmm &
       &                           , sw_hybrid_functional, sw_eval_vexx, sw_retard_eigval_evaluation, sw_fef &
       &                           , in_line_minimization &
       &                           , sw_external_potential, ekmode, icond, sw_rsb, sw_output_xc_seperately &
       &                           , number_of_cdmixing_applied, cdmixing_names_applied &
       &                           , number_of_solvers_applied, solver_names_applied, len_solvername &
       &                           , sw_potential_mixing &
#ifdef ENABLE_ESM_PACK
       &                           , sw_esm &
#endif
       &                           , m_CtrlP_push_CDMixingNameApplied,sw_harris_functional &
       &                           , nsamp, convergence_criteria &
       &                           , sw_vdw_correction, vdw_method &
       &                           , order_mp
  use m_Kpoints,              only : qwgt,kv3,vkxyz,k_symmetry
  use m_IterationNumbers,     only : iteration, iteration_electronic,iteration_ionic
  use m_Parallelization,      only : MPI_CommGroup,map_k,myrank_k,ierr,np_e &
       &                           , ista_kngp,iend_kngp,npes,mype &
       &                           , myrank_e, map_e, map_z   &
       &                           , ista_atm, iend_atm, myrank_g, nrank_g, nrank_e &
       &                           , ista_e, iend_e, istep_e &
       &                           , mpi_spin_group, ista_spin, iend_spin
  use m_Const_Parameters,     only : Valence_plus_PC_Charge,PAI4,UP,DOWN&
       &                           , EXC_ONLY, VXC_AND_EXC, IINCRE_CRITICAL,DP,CMPLDP &
       &                           , len_tag_total_energy, tag_total_energy &
       &                           , PARABOLIC, COLD, ON, OFF &
       &                           , FIXED_CHARGE, FIXED_CHARGE_CONTINUATION, GAMMA, DELTA &
       &                           , MDKOSUGI, MDDAVIDSON, INITIAL, CONTINUATION, COORDINATE_CONTINUATION &
       &                           , FIXED_CHARGE, FIXED_CHARGE_CONTINUATION, MATRIXDIAGON &
       &                           , DELTA_MOVING_AVERAGE, SLOPE, DELTA_V &
       &                           , VDW_DFTD3 &
       &                           , MP
  use m_Dipole,               only : edip_ion, eext_ion, vdip_l, vext_l
  use m_Screening,            only : screening
  use m_Hubbard,              only : m_Hubbard_energy
  use m_FiniteElectricField,  only : m_FEF_polarization, pmac, pmac_old
  use m_ES_ExactExchange,     only : m_ES_EXX_gather_valence_states,m_ES_EXX_energy &
       &                           , m_ES_EXX_energy2
#if 0
       &  , m_ES_EXX_gather_valence_states_k
#endif
  use m_PAW_XC_Potential,     only : m_PAW_XC_cal_potential,exc_ae,exc_ps &
       &                           , m_PAW_XC_cal_potential_sphex2
!!$                                    , m_PAW_XC_cal_potential_sym
  use m_PAW_Hartree,          only : m_PAWH_get_dion_hartree_now
  use m_External_Potential,   only : espot_g
  use m_PAW_ChargeDensity,    only : calcGaussLegendreIntegration &
       &                            , calcSphericalHarmonicsExpansion


! ====================================== added by K. Tagami ================ 11.0
  use m_Control_Parameters,    only : noncol, ndim_spinor, ndim_magmom, ndim_chgpot
  use m_PseudoPotential,       only : dion0_noncl, nlmt
  use m_Charge_Density,         only : hsi
  use m_ES_NonCollinear,       only : m_ES_MagMom_To_DensMat_hsr, &
       &                              m_ES_MagMom_to_DensMat_Dhub
!
  use m_Crystal_Structure,     only :  sw_magnetic_constraint
  use m_ES_Mag_Constraint,      only : m_ES_calc_MagConstraint_Energy
  use m_Hubbard,              only : m_Hubbard_energy_noncl, &
       &                             m_Hubbard_energy2_noncl, &
       &                             m_Hubbard_energy3_noncl
! ========================================================================== 11.0

! ====================================== added by K. Tagami ================ 11.0
  use m_Control_Parameters,    only : SpinOrbit_Mode
  use m_Const_Parameters,      only : ByPawPot, Neglected, CMPLDP, BuiltIn, EXECUT, &
       &                              ByProjector, ZeffApprox, ReadFromPP
  use m_SpinOrbit_Potential,   only : dsoc, m_SO_set_Dsoc_potential2
  use m_SpinOrbit_RadInt,    only : m_SO_calc_SOC_strength_pawpot
  use m_PseudoPotential,      only : lmta,  m_PP_include_vanderbilt_pot, dion_scr_noncl
  use m_Electronic_Structure,   only : fsr_l, fsi_l, dhub_aimag
! ========================================================================== 11.0


! ======================= KT_add ================== 13.0E
  use m_Const_Parameters,   only : Fermi_Dirac
! ================================================= 13.0E

! =========== KT_add ========== 13.0U2
  use m_Control_Parameters,  only : sw_potential_mixing, use_metagga, vtau_exists
  use m_XC_Potential,        only : vtau_l
  use m_KineticEnergy_Density,  only : ekins_l, ekins_old
! ============================= 13.0U2

#ifndef DISABLE_VDWDF
  use m_vdWDF,  only : ecnl_vdwdf
#endif

  use m_Control_Parameters,  only : m_CtrlP_get_isolver_now

! === Positron SCF ==== 2015/11/28
  use m_Control_Parameters,  only : sw_positron, npeg, positron_method
  use m_Const_Parameters,   only :  positron_GGGC
  use m_epc_potential,  only : ecorr_pztr => epc, m_epc_cal_potential, vepc_l
  use m_PlaneWaveBasisSet, only : kg1_pwf
  use m_Positron_Wave_Functions,  only : pzaj, nprvf_ordr, pchg_l, pchgo_l, pev
  use m_PlaneWaveBasisSet,    only : m_pwBS_pstrn_kinetic_energies
! ==================== 2015/11/28
  use mpi

  implicit none
!  include 'mpif.h'

  real(kind=DP),private          :: eband ! band energy
  real(kind=DP),private          :: eband_extendedrange
  real(kind=DP),private          :: eohxc ! exchange correlation and &
  !                                         Hartree with old charge

! ==== KT_mod === 13.0U2
!!  real(kind=DP),private          :: elocal ! local potential
  real(kind=DP)            :: elocal ! local potential
! =============== 13.0U2

  real(kind=DP),private          :: eloca1
  real(kind=DP),private          :: enonlc ! non-local
!!$  real(kind=DP),private          :: ehartr ! hartree energy with new charge
  real(kind=DP)                  :: ehartr ! hartree energy with new charge
  real(kind=DP),private          :: ekinet ! kinetic energy
  real(kind=DP),private          :: eentropy ! entropic term of free energy (-TS)

  real(kind=DP),private          :: edip   ! dipole energy
  real(kind=DP),private          :: evdip  ! dipole potential
  real(kind=DP),private          :: evext  ! external potential

  real(kind=DP),private          :: etoold ! previous total energy
  real(kind=DP),private          :: edeltb = 1.d99 ! etotal - etoold
  real(kind=DP)                  :: etotal ! total free energy: F(simga)
  real(kind=DP)                  :: etotal0 ! total energy: E(sigma=0)

  real(kind=DP)                  :: eespot ! electrostatic potential

  real(kind=DP), pointer, dimension(:,:,:) :: chgqt,chgqto
  real(kind=DP), pointer, dimension(:,:,:,:) :: hsrt,hsit
  !
  ! Free energy:
  !F(sigma) = E(sigma) - sigma*S(sigma)
  !
  ! 1. E(sigma=0) exprapolation for the cold smearing method
  !E(sigma=0) ~ (2*F(sigma)+E(sigma))/3
  !
  ! 2. E(sigma=0) exprapolation for the parabolic broadening method
  !F(sigma) = E(sigma) - sigma*S(sigma)
  !E(sigma=0) ~ (F(sigma)+E(sigma))/2
  !

  ! DFT+U
  real(kind=DP),private          :: ehub0 ! Hubbard energy
  real(kind=DP),private          :: ehub0_old = 1.d+20 ! Old Hubbard energy
  real(kind=DP),private          :: ehub1 ! Hubbard potential energy

  ! PAW
  real(kind=DP),private          :: eohxc_paw, ehartr_paw, ekin_ion_paw

  ! Hybrid functional method
  real(kind=DP),private          :: vexx,eexx

  ! Finite electric field method
  real(kind=DP),private          :: eplr  ! -EP term

! ============================== added by K. Tagami ================== 11.0
! Magnetic constraint
!
  real(kind=DP),private          :: emag0      ! magnetic constraint energy
  real(kind=DP),private          :: emag1      ! double counting energy
! ==================================================================== 11.0

! ====================================== added by K. Tagami ================ 11.0
! Spin Orbit
!
  real(kind=DP),private          :: espinorb_old, espinorb_now
! ========================================================================== 11.0

! === positron
  real(kind=DP) :: ekin_pztr, elocal_pztr, ehartr_ep, eohxc_pztr
! ===

  integer,private, parameter     :: len_str = 132
  character(len=len_str),private    ::  str

  real(kind=DP), allocatable, dimension(:) :: ehist
  integer :: ihist
  real(kind=DP) :: emova, emovaold
contains
  real(DP) function m_TE_what_is_edeltb_now()
     m_TE_what_is_edeltb_now = edeltb
  end function m_TE_what_is_edeltb_now

  subroutine m_TE_wd_total_energy_with_solvers(nfout)
    integer, intent(in) :: nfout
    logical :: display_on
    display_on = .true.
    call sumup_all_energies(nfout,exc,display_on)

  end subroutine m_TE_wd_total_energy_with_solvers

  subroutine m_TE_total_energy(nfout,display_on,kv3)
#if 0
    use m_Control_Parameters, only : sw_distribute_wf, force_exx_energy1
#else
    use m_Control_Parameters, only : force_exx_energy1
#endif
    integer, intent(in) :: nfout
    logical, intent(in) :: display_on
    integer, intent(in) :: kv3
    integer :: isolver,sw_submat
    real(kind=DP) :: edeltb_now
    integer             :: id_sname = -1
    call tstatc0_begin('Total_Energy(including xc_pot) ',id_sname)
    if(sw_harris_functional == OFF)then
      chgqt => chgq_l
      chgqto => chgqo_l
      hsrt => hsr
      if(noncol) hsit => hsi
    else
      chgqt => chgqo_l
      chgqto => chgqo_l
      hsrt => hsro
      if(noncol) hsit => hsi
    endif
    call get_band_energy(kv3)
    call get_xc_and_HE_of_old_CD(vxc_l)
    call get_local_potential_energy
    call get_hubbard_energy(nfout)
    call get_nonlocal_potential_energy
    if(in_line_minimization) then
      call m_XC_cal_potential(nfout,Valence_plus_PC_Charge, chgqt, EXC_ONLY)
    else
      call m_XC_cal_potential(nfout,Valence_plus_PC_Charge,chgqt,VXC_AND_EXC)
    endif
#ifdef ENABLE_ESM_PACK
    if (sw_esm==OFF) then
       call get_hartree_energy
    endif
#else
    call get_hartree_energy
#endif
    if(sw_dipole_correction ==  ON) call get_dipole_energy(vdip_l,vext_l)
    if(sw_fef == ON) then
       call m_FEF_polarization(nfout,eplr)
    endif
    if(sw_hybrid_functional == ON) then
       if(sw_eval_vexx==ON) call m_ES_EXX_energy(vexx)
       if(icond == INITIAL .or. icond == CONTINUATION .or. icond == COORDINATE_CONTINUATION .or. &
       & ((icond==FIXED_CHARGE.or.icond==FIXED_CHARGE_CONTINUATION)&
       &                             .and. ekmode == OFF)   ) then
          edeltb_now = m_TE_what_is_edeltb_now()/natm2
!!$     write(nfout,'(" edeltb_per_atom = ", d20.12)') edeltb_per_atom
       else
          edeltb_now = m_ES_what_is_evdff_now()
       end if
       !isolver = m_CtrlP_solver_for_WFs_now(iteration_electronic,iteration_ionic,0 &
       !&                              ,edeltb_now,sw_submat)
       isolver = m_CtrlP_get_isolver_now()
       if(sw_retard_eigval_evaluation==OFF.or.isolver==MDDAVIDSON.or.isolver==MDKOSUGI.or. &
       &  force_exx_energy1.or.isolver==MATRIXDIAGON) then
#if 0
          if(sw_distribute_wf==ON)then
             call m_ES_EXX_gather_valence_states_k(nfout)
          else
             call m_ES_EXX_gather_valence_states(nfout)
          endif
#else
          call m_ES_EXX_gather_valence_states(nfout)
#endif

          call m_ES_EXX_energy(eexx)
       else
          call m_ES_EXX_energy2(eexx)
       endif
    end if
    force_exx_energy1 = .false.

    if ( sw_magnetic_constraint == ON ) call get_magnetic_constraint_energy

    call get_entropic_term(nfout)
    if(flg_paw) then
#if 0
        call get_xc_and_HE_of_old_CD_paw(nfout)
!       call get_xc_and_HE_of_old_CD_paw_sym
#endif
        if(calcSphericalHarmonicsExpansion) then
          call m_PAW_XC_cal_potential_sphex2(nfout,EXC_ONLY)
!    call m_PAW_XC_cal_potential_sphex3(nfout,VXC_AND_EXC)
        end if
        if(calcGaussLegendreIntegration)then
           call m_PAW_XC_cal_potential(nfout,EXC_ONLY,flg_symmtry)
        endif
!!$        if(.not.flg_symmtry) then
!!$            call m_PAW_XC_cal_potential(nfout,EXC_ONLY)
!!$        else
!!$            call m_PAW_XC_cal_potential_sym(nfout,EXC_ONLY)
!!$        end if
        call m_PAWH_get_dion_hartree_now(nfout)
        call get_hartree_energy_paw(nfout)
!       call get_hartree_energy_paw_sym
!
#if 1
        call get_kinetic_local_energy_paw(nfout)
#endif
    end if

    if ( sw_positron /= OFF ) then
       if ( positron_method == Positron_GGGC ) then
          call get_local_pot_energy_pztr
          call get_xc_and_HE_of_old_CD_pztr( vepc_l )
          call get_hartree_energy_inter_ep
          call m_epc_cal_potential( nfout, chgqt )

         if (iteration_electronic==1) then
            call get_kinetic_energy_pztr_direct(nfout)
         else
            call get_kinetic_energy_pztr(nfout)
         endif
       endif
    endif

#if 0
    if ( use_metagga .and. vtau_exists ) then
       call get_mgga_ekindens_pot_energy
    endif
#endif

#ifdef ENABLE_ESM_PACK
! ============== KT_mod =============================== 13.0U2
!    if(sw_esm==ON.or.sw_hybrid_functional==ON)then
    if ( sw_esm==ON .or. sw_hybrid_functional==ON &
         &          .or. sw_potential_mixing==ON  &
         &          .or. sw_rsb==ON &
         &          .or. xctype=='vdwdf' &
         &          .or. use_metagga ) then
! ===================================================== 13.0U2
       call get_kinetic_energy_directly(nfout)
    else
       if(iteration_electronic==1) then
          call get_kinetic_energy_directly(nfout)
       else
          call get_kinetic_energy(nfout)
       endif
    endif
#else
! ============== KT_mod =============================== 13.0U2
!    if(sw_hybrid_functional==ON.or.iteration_electronic==1)then
    if ( sw_hybrid_functional==ON .or. iteration_electronic==1 &
         &                        .or. sw_potential_mixing==ON &
         &                        .or. sw_rsb==ON &
         &                        .or. xctype=='vdwdf' &
         &                        .or. use_metagga ) then
! ===================================================== 13.0U2
       call get_kinetic_energy_directly(nfout)
    else
       call get_kinetic_energy(nfout)
    endif
#endif

    call sumup_all_energies(nfout,exc,display_on)

    call tstatc0_end(id_sname)
  end subroutine m_TE_total_energy

! =============================== added by K. Tagami ================= 11.0
  subroutine m_TE_total_energy_noncl(nfout,display_on,kv3)
#if 0
    use m_Control_Parameters, only : sw_distribute_wf
#endif

    integer, intent(in) :: nfout
    logical, intent(in) :: display_on
    integer, intent(in) :: kv3
    integer             :: id_sname = -1
    call tstatc0_begin('Total_Energy(including xc_pot) ',id_sname)
    if(sw_harris_functional == OFF)then
      chgqt => chgq_l
      chgqto => chgqo_l
      hsrt => hsr
      if(noncol) hsit => hsi
    else
      chgqt => chgqo_l
      chgqto => chgqo_l
      hsrt => hsro
      if(noncol) hsit => hsi
    endif

    if(sw_harris_functional == OFF)then
      chgqt => chgq_l
      chgqto => chgqo_l
      hsrt => hsr
      if(noncol) hsit => hsi
    else
      chgqt => chgqo_l
      chgqto => chgqo_l
      hsrt => hsro
      if(noncol) hsit => hsi
    endif

    call get_band_energy_noncl(kv3)
    call get_xc_and_HE_of_old_CD_noncl(vxc_l)
    call get_local_pot_energy_noncl
    call get_hubbard_energy_noncl(nfout)

#if 0
    call get_nonlocal_pot_energy_noncl
#else
    call get_nonlocal_pot_energy_noncl2
#endif

    call m_XC_cal_potential( nfout,Valence_plus_PC_Charge, chgqt, EXC_ONLY)

#ifdef ENABLE_ESM_PACK
    if (sw_esm==OFF) then
       call get_hartree_energy_noncl
    endif
#else
    call get_hartree_energy_noncl
#endif

    if(sw_dipole_correction ==  ON) then
       call get_dipole_energy_noncl(vdip_l,vext_l)
    end if

    if(sw_fef == ON) then
       write(*,*) 'Not supported : FEF '
       call phase_error_with_msg(nfout,'FEF not supported',__LINE__,__FILE__)
       call m_FEF_polarization(nfout,eplr)
    endif

    if(sw_hybrid_functional == ON) then
       write(*,*) 'Not supported : hybrid '
       call phase_error_with_msg(nfout,'hybrid not supported',__LINE__,__FILE__)
       if(sw_eval_vexx==ON) call m_ES_EXX_energy(vexx)
#if 0
       if(sw_distribute_wf==ON)then
          call m_ES_EXX_gather_valence_states_k(nfout)
       else
          call m_ES_EXX_gather_valence_states(nfout)
       endif
#else
       call m_ES_EXX_gather_valence_states(nfout)
#endif

       call m_ES_EXX_energy(eexx)
    end if

    if ( sw_magnetic_constraint == ON ) call get_magnetic_constraint_energy

    call get_entropic_term(nfout)

    if (flg_paw) then
!       stop 'Not supported : PAW'
#if 0
       call get_xc_and_HE_old_CD_paw_noncl(nfout)
#endif
       if (calcSphericalHarmonicsExpansion) then
          call m_PAW_XC_cal_potential_sphex2(nfout,EXC_ONLY)

       end if
       if (calcGaussLegendreIntegration)then
          call phase_error_with_msg(nfout,'Not supported : Gauss Legendre in paw + noncol',__LINE__,__FILE__)
          call m_PAW_XC_cal_potential(nfout,EXC_ONLY,flg_symmtry)
       endif

       call m_PAWH_get_dion_hartree_now(nfout)
       call get_hartree_energy_paw_noncl(nfout)
#if 1
       call get_kinetic_local_energy_paw(nfout)
#endif
    end if
!
!!!#ifdef USE_ESPINORB
    if ( flg_paw ) then
       if ( SpinOrbit_Mode == ByPawPot ) then
#if 1
          call get_spinorbit_energy_noncl2(espinorb_old)
#else
          call get_spinorbit_energy_noncl3(espinorb_old)
#endif
          call m_SO_calc_SOC_strength_pawpot
          call m_SO_set_Dsoc_potential2
!
#if 1
          call get_spinorbit_energy_noncl2(espinorb_now)
#else
          call get_spinorbit_energy_noncl3(espinorb_now)
#endif
      endif
    endif

    if ( SpinOrbit_Mode == ByProjector .or. &
         &  SpinOrbit_Mode == ZeffApprox .or. &
         &  SpinOrbit_Mode == ReadFromPP ) then
#if 1
       call get_spinorbit_energy_noncl2(espinorb_now)
#else
       call get_spinorbit_energy_noncl3(espinorb_now)
#endif
    endif
!!!#endif

#ifdef ENABLE_ESM_PACK
! ============== KT_mod =============================== 13.0U2
!    if(sw_esm==ON.or.sw_hybrid_functional==ON)then
    if ( sw_esm==ON .or. sw_hybrid_functional==ON &
         &          .or. sw_potential_mixing==ON  &
         &          .or. sw_rsb==ON &
         &          .or. xctype=='vdwdf') then
! ===================================================== 13.0U2
       call get_kinetic_energy_directly(nfout)
    else
       if(iteration_electronic==1) then
          call get_kinetic_energy_directly(nfout)
       else
          call get_kinetic_energy(nfout)
       endif
    endif
#else
! ============== KT_mod =============================== 13.0U2
!    if(sw_hybrid_functional==ON.or.iteration_electronic==1)then
    if ( sw_hybrid_functional==ON .or. iteration_electronic==1 &
         &                        .or. sw_potential_mixing==ON &
         &                        .or. sw_rsb==ON &
         &                        .or. xctype=='vdwdf') then
! ===================================================== 13.0U2
       call get_kinetic_energy_directly(nfout)
    else
       call get_kinetic_energy(nfout)
    endif
#endif
!!!    call get_kinetic_energy(nfout)

    call sumup_all_energies(nfout,exc,display_on)

    call tstatc0_end(id_sname)

  end subroutine m_TE_total_energy_noncl
! ==================================================================== 11.0

  subroutine get_band_energy(kv3)
    integer, intent(in) :: kv3
    integer             :: ib, ik, is
    real(kind=DP)       :: eband_mpi   ! MPI
    eband = 0.d0
    do is = ista_spin, iend_spin
    do ik = is, kv3-nspin+is, nspin
       if(map_k(ik) /= myrank_k) cycle ! MPI
       do ib = 1, np_e                 ! MPI
          eband = eband + occup_l(ib,ik)*eko_l(ib,ik)
       end do
    end do
    enddo

    if(npes > 1) then
       call mpi_allreduce(eband,eband_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr) !MPI
       eband = eband_mpi
    end if
!!$    eband = eband_mpi                  ! MPI
    eband = 2*eband/kv3 * (af+1)
  end subroutine get_band_energy

! ======================================== added by K. Tagami ============== 11.0
  subroutine get_band_energy_noncl(kv3)
    integer, intent(in) :: kv3
    integer             :: ib, ik, is
    real(kind=DP)       :: eband_mpi   ! MPI

    eband = 0.d0

    do is = ista_spin, iend_spin
    do ik = is, kv3-nspin+is, ndim_spinor
!    do ik =1, kv3, ndim_spinor
       if(map_k(ik) /= myrank_k) cycle ! MPI
       do ib = 1, np_e                 ! MPI
          eband = eband + occup_l(ib,ik)*eko_l(ib,ik)
       end do
    end do
    enddo

    if (npes > 1) then
       call mpi_allreduce( eband, eband_mpi, 1, mpi_double_precision, &
	&                  mpi_sum, MPI_CommGroup, ierr ) !MPI
       eband = eband_mpi
    end if
!!$    eband = eband_mpi                  ! MPI
    eband = eband/ (kv3 /ndim_spinor )

  end subroutine get_band_energy_noncl
! ===================================================================== 11.0

  subroutine get_xc_and_HE_of_old_CD(vxc_l)
    real(kind=DP), intent(in) :: vxc_l(ista_kngp:iend_kngp,kimg,nspin)
    integer ik, i, ispin
    integer ist !mpi
    real(kind=DP) :: eohxc_mpi

    eohxc = 0.d0
    do ik = 1, kimg
       do ispin = 1, nspin
          if(mype==0) then
             eohxc = eohxc + vxc_l(1,ik,ispin)*chgqt(1,ik,ispin)
             if(sw_external_potential==ON) then
                eohxc = eohxc + espot_g(1,ik,ispin)*chgqt(1,ik,ispin)
             end if
          endif
       end do
       ist = ista_kngp
       if(ist == 1) ist = 2

       if(nspin == 1) then
          do i = ist, iend_kngp  !for mpi
             eohxc = eohxc &
                  & + (vxc_l(i,ik,1)+PAI4*chgqto(i,ik,1)/gr_l(i)**2) &
                  &   * chgqt(i,ik,1)
             if(sw_external_potential==ON) then
                eohxc = eohxc + espot_g(i,ik,1) * chgqt(i,ik,1)
             end if
             if(sw_screening_correction==ON) then
                 eohxc = eohxc &
                  & + chgqto(i,ik,1)*chgqt(i,ik,1)*screening%phik(i)
             end if
             if ( sw_positron /= OFF ) then
                if ( positron_method == positron_GGGC ) then
                   eohxc = eohxc &
                        & - (PAI4 *pchgo_l(i,ik)/gr_l(i)**2) &
                        &   * chgqt(i,ik,1)
                endif
             endif
          end do
       else if(nspin == 2) then
          do i = ist, iend_kngp  !for mpi
             eohxc = eohxc &
                  & +  (vxc_l(i,ik,UP)*chgqt(i,ik,UP)&
                  &   + vxc_l(i,ik,DOWN)*chgqt(i,ik,DOWN)) &
                  & + PAI4*(chgqto(i,ik,UP)+chgqto(i,ik,DOWN))&
                  & /gr_l(i)**2 &
                  &  *(chgqt(i,ik,UP)+chgqt(i,ik,DOWN))
             if(sw_external_potential==ON) then
                eohxc = eohxc + (espot_g(i,ik,UP)*chgqt(i,ik,UP) &
                              + espot_g(i,ik,DOWN)*chgqt(i,ik,DOWN))
             end if
             if(sw_screening_correction==ON) then
                 eohxc = eohxc &
                  & + (chgqto(i,ik,UP)+chgqto(i,ik,DOWN))* &
                  &   (chgqt(i,ik,UP)+chgqt(i,ik,DOWN))*screening%phik(i)
             end if
             if ( sw_positron /= OFF ) then
                if ( positron_method == positron_GGGC ) then
                   eohxc = eohxc &
                        & - (PAI4 *pchgo_l(i,ik)/gr_l(i)**2) &
                        &    *(chgqt(i,ik,UP)+chgqt(i,ik,DOWN))
                end if
             endif
          end do
       end if
    end do
    if(npes > 1) then
       call mpi_allreduce(eohxc,eohxc_mpi,1 &
            &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       eohxc = eohxc_mpi
    end if

    eohxc = univol*eohxc
  end subroutine get_xc_and_HE_of_old_CD

! -----------------
!   meta-gga
! -----------------
#if 0
  subroutine get_mgga_ekindens_pot_energy
    integer :: is, ri, ig
    real(kind=DP) :: ene_wk

    if ( .not.allocated(ekins_old)) write(*,*) "Ekins_old aa"
    if ( .not.allocated(vtau_l)) write(*,*) "vatu_l aa"

    ene_mgga_tau = 0.0d0;     ene_mgga_tau_old = 0.0d0
    do is = 1, ndim_magmom
       do ri = 1, kimg
          do ig = ista_kngp, iend_kngp
             ene_mgga_tau_old = ene_mgga_tau_old +vtau_l(ig,ri,is) *ekins_old(ig,ri,is)
             ene_mgga_tau     = ene_mgga_tau     +vtau_l(ig,ri,is) *ekins_l(ig,ri,is)
          end do
       end do
    end do
    if(npes > 1) then
       call mpi_allreduce( ene_mgga_tau, ene_wk,1 &
            &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       ene_mgga_tau = ene_wk
       call mpi_allreduce( ene_mgga_tau_old, ene_wk,1 &
            &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       ene_mgga_tau_old = ene_wk
    end if
  end subroutine get_mgga_ekindens_pot_energy
#endif

! ============================== added by K. Tagami ====================== 11.0
  subroutine get_xc_and_HE_of_old_CD_noncl(vxc_l)
    real(kind=DP), intent(in) :: vxc_l(ista_kngp:iend_kngp,kimg,ndim_magmom)
    integer ik, is
    integer ist !mpi

    integer :: ri, ig, ni

    real(kind=DP) :: eohxc_mpi
    real(kind=DP) :: csum1, csum2

    eohxc = 0.d0

    if ( sw_external_potential == ON ) then
       call phase_error_with_msg(6,'kt : external_pot not supported ',__LINE__,__FILE__)
    endif

    do is = 1, ndim_magmom
       do ri = 1, kimg
          do ig = ista_kngp, iend_kngp
             eohxc = eohxc + vxc_l(ig,ri,is)*chgqt(ig,ri,is)
          end do
       end do
    end do

    do ri=1, kimg
       ist = ista_kngp
       if(ist == 1) ist = 2

       do ig=ist, iend_kngp
         ni = 1
         csum1 = chgqt( ig, ri, ni )
         csum2 = chgqto( ig, ri, ni )
         eohxc = eohxc + PAI4 *csum2 /gr_l(ig)**2 *csum1

       end do
    end do

    if ( sw_screening_correction == ON ) then
       do ri=1, kimg
          ist = ista_kngp
          if(ist == 1) ist = 2

          do ig=ist, iend_kngp
             ni = 1
             csum1 = chgqt( ig, ri, ni )
             csum2 = chgqto( ig, ri, ni )
             eohxc = eohxc + csum2 *csum1 *screening%phik(ig)
          end do
       end do
    endif

    if (npes > 1) then
       call mpi_allreduce( eohxc, eohxc_mpi, 1 &
            &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr )
       eohxc = eohxc_mpi
    end if

    eohxc = univol*eohxc
  end subroutine get_xc_and_HE_of_old_CD_noncl
! ====================================================================== 11.0

  subroutine get_local_potential_energy

    integer       :: ik, it, i, ig
    real(kind=DP) :: eloca1_mpi
    real(kind=DP) :: eespot_mpi
    real(kind=DP) :: eloclr,eloclr_mpi

    eloca1 = 0.d0
    do ik = 1, kimg
       if(nspin == 1) then
          do it = 1,ntyp
             do i = ista_kngp, iend_kngp !for mpi
                eloca1 = eloca1 &
                     & + psc_l(i,it)*zfm3_l(i,it,ik)*chgqt(i,ik,1)
                if(sw_screening_correction==ON) then
                  eloca1 = eloca1 &
                      & - ival(it)*screening%phik(i)/univol*zfm3_l(i,it,ik)*chgqt(i,ik,1)
                end if
             end do
          end do
       else
          do it = 1, ntyp
             do i = ista_kngp, iend_kngp  !for mpi
                eloca1 = eloca1 &
                     &   + psc_l(i,it)*zfm3_l(i,it,ik)&
                     &             *(chgqt(i,ik,UP)+chgqt(i,ik,DOWN))
                if(sw_screening_correction==ON) then
                   eloca1 = eloca1 &
                        &   - ival(it)*Screening%phik(i)/univol*zfm3_l(i,it,ik)&
                        &             *(chgqt(i,ik,UP)+chgqt(i,ik,DOWN))
                end if
             end do
          end do
       endif
    end do
    if(npes > 1) then
       call mpi_allreduce(eloca1,eloca1_mpi,1 &
            &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       eloca1 = eloca1_mpi
    end if

    eloca1    = univol*eloca1
#ifdef ENABLE_ESM_PACK
    if(sw_esm==OFF) then
       elocal = eloca1 + etot1*totch
    else
       elocal = eloca1
    endif
#else
    elocal    = eloca1 + etot1*totch
#endif

#ifdef ENABLE_ESM_PACK
    if(sw_esm==ON)then
!add the long-range part
      eloclr=0.d0
      if(nspin.eq.1)then
        if(kimg==2)then
          do ig=ista_kngp,iend_kngp
            eloclr=eloclr+dble(vloc_esm(igfp_l(ig)))*chgqt(ig,1,1)+aimag(vloc_esm(igfp_l(ig)))*chgqt(ig,2,1)
          enddo
        else if(kimg==1)then
          do ig=ista_kngp,iend_kngp
            eloclr=eloclr+dble(vloc_esm(igfp_l(ig)))*chgqt(ig,1,1)
          enddo
        endif
      else
        if(kimg==2)then
          do ig=ista_kngp,iend_kngp
            eloclr=eloclr+dble(vloc_esm(igfp_l(ig))) &
   &                     *(chgqt(ig,1,1)+chgqt(ig,1,2)) &
   &                     +aimag(vloc_esm(igfp_l(ig))) &
   &                     *(chgqt(ig,2,1)+chgqt(ig,2,2))
          enddo
        else if(kimg==1)then
          do ig=ista_kngp,iend_kngp
            eloclr=eloclr+dble(vloc_esm(igfp_l(ig))) &
   &                     *(chgqt(ig,1,1)+chgqt(ig,1,2))
          enddo
        endif
      endif
      eloclr_mpi=0.d0
      call mpi_allreduce(eloclr,eloclr_mpi,1,&
      & mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
      eloclr = eloclr_mpi
      eloclr = univol*eloclr
      elocal = elocal + eloclr
    endif
#endif

    eespot = 0.0d0
    if(sw_external_potential==ON) then
    do ik = 1, kimg
!       do ispin = 1, nspin
!          eespot = eespot + espot_g(1,ik,ispin)*chgq_l(1,ik,ispin)
!       end do
       if(nspin == 1) then
!          do i = 2, iend_kngp
          do i = ista_kngp, iend_kngp
            eespot = eespot + espot_g(i,ik,1)*chgqt(i,ik,1)
          end do
       else if( nspin == 2 )then
!          do i = 2, iend_kngp
          do i = ista_kngp, iend_kngp
            eespot = eespot + ( espot_g(i,ik,UP)*chgqt(i,ik,UP) + espot_g(i,ik,DOWN)*chgqt(i,ik,DOWN ) )
          end do
       end if
    end do
    if(npes > 1) then
       call mpi_allreduce(eespot,eespot_mpi,1 &
            &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       eespot = eespot_mpi
    end if

    eespot = univol*eespot
    elocal = elocal + eespot
    end if

  end subroutine get_local_potential_energy

! ====================================== added by K. Tagami ================= 11.0
  subroutine get_local_pot_energy_noncl

    integer       :: ri, it, ig
    real(kind=DP) :: eloca1_mpi
    real(kind=DP) :: csum1
    real(kind=DP) :: eespot_mpi
    real(kind=DP) :: eloclr,eloclr_mpi

    if ( sw_external_potential == ON ) then
       call phase_error_with_msg(6,'kt : external_pot not supported 2 ',__LINE__,__FILE__)
    endif

    eloca1 = 0.d0
    do ri = 1, kimg
      do it = 1, ntyp
         do ig = ista_kngp, iend_kngp  !for mpi
           csum1 = chgqt( ig, ri, 1 )
           eloca1 = eloca1  + psc_l(ig,it)*zfm3_l(ig,it,ri) *csum1
         end do
      end do
    end do

    if ( sw_screening_correction == ON ) then
       do ri = 1, kimg
          do it = 1, ntyp
             do ig = ista_kngp, iend_kngp  !for mpi
                csum1 = chgqt( ig, ri, 1 )
                eloca1 = eloca1 -ival(it) *Screening%phik(ig) /univol &
                     &           *zfm3_l(ig,it,ri) *csum1
             end do
          end do
       end do
    endif

    if (npes > 1) then
       call mpi_allreduce(eloca1,eloca1_mpi,1 &
            &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       eloca1 = eloca1_mpi
    end if

    eloca1    = univol*eloca1

#ifdef ENABLE_ESM_PACK
    if(sw_esm==OFF) then
       elocal = eloca1 + etot1*totch
    else
       elocal = eloca1
    endif
#else
    elocal    = eloca1 + etot1*totch
#endif

#ifdef ENABLE_ESM_PACK
    if(sw_esm==ON)then
!add the long-range part
      eloclr=0.d0
      if(kimg==2)then
         do ig=ista_kngp,iend_kngp
            eloclr = eloclr +dble(vloc_esm(igfp_l(ig))) *chgqt(ig,1,1) &
                 &          +aimag(vloc_esm(igfp_l(ig)))*chgqt(ig,2,1)
         enddo
      else if(kimg==1)then
         do ig=ista_kngp,iend_kngp
            eloclr=eloclr+dble(vloc_esm(igfp_l(ig)))*chgqt(ig,1,1)
         enddo
      endif
      eloclr_mpi=0.d0
      call mpi_allreduce(eloclr,eloclr_mpi,1,&
      & mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
      eloclr = eloclr_mpi
      eloclr = univol*eloclr
      elocal = elocal + eloclr
    endif
#endif

    eespot = 0.0d0
    if(sw_external_potential==ON) then
       call phase_error_with_msg(6,"Not implenented for sw_external_potential = ON",__LINE__,__FILE__)
    endif

  end subroutine get_local_pot_energy_noncl
! ============================================================================ 11.0

  subroutine get_nonlocal_potential_energy
    integer :: ispin, it,lmt1, lmt2, il1, im1, il2, im2, ia
    real(kind=DP) :: fac
    enonlc = 0.d0
    do ispin = 1, nspin, af+1
       do it = 1, ntyp
          if(ipaw(it)==0) then
              do lmt1 = 1, ilmt(it)
                 il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
                 do lmt2 = lmt1, ilmt(it)
                    il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
                    if( il1 /= il2 .or. im1 /= im2) cycle
    !xocl spread do/ind_katm
                    do ia = 1, natm
                       if(ityp(ia) /= it) cycle
                       fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0
                       enonlc = enonlc &
                            & + fac*dion(lmt1,lmt2,it)*hsrt(ia,lmt1,lmt2,ispin)
                    end do
    !xocl end spread
                 end do
              end do
          else
              do lmt1 = 1, ilmt(it)
                 do lmt2 = lmt1, ilmt(it)
    !xocl spread do/ind_katm
                    do ia = 1, natm
                       if(ityp(ia) /= it) cycle
                       fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0
                       enonlc = enonlc &
                            & + fac*dion_paw(lmt1,lmt2,ispin,ia)*hsrt(ia,lmt1,lmt2,ispin)
                    end do
    !xocl end spread
                 end do
              end do
          end if
!!$          do lmt1 = 1, ilmt(it)
!!$             il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
!!$             do lmt2 = lmt1, ilmt(it)
!!$                il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
!!$                if( il1 /= il2 .or. im1 /= im2) cycle
!!$!xocl spread do/ind_katm
!!$                do ia = 1, natm
!!$                   if(ityp(ia) /= it) cycle
!!$                   fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0
!!$                   enonlc = enonlc &
!!$                        & + fac*dion(lmt1,lmt2,it)*hsr(ia,lmt1,lmt2,ispin)
!!$                end do
!!$!xocl end spread
!!$             end do
!!$          end do
       end do
    end do
    enonlc = enonlc*(af+1)
    !!$if(sw_hubbard==ON) enonlc = enonlc - ehub1
  end subroutine get_nonlocal_potential_energy

! ================================== added by K. Tagami =================== 11.0
  subroutine get_nonlocal_pot_energy_noncl

    integer :: it,lmt1, lmt2, il1, im1, il2, im2, ia
    integer :: is
    real(kind=DP) :: fac
    real(kind=DP) :: c1, c2

    real(kind=DP), allocatable :: hsr_densmat( :,:,:,: )
    real(kind=DP), allocatable :: hsi_densmat( :,:,:,: )

    enonlc = 0.d0

! --
!   hsr (magomom) --> hsr, hs ( chgopt)
!

    allocate( hsr_densmat( natm, nlmt, nlmt, ndim_chgpot ) )
    allocate( hsi_densmat( natm, nlmt, nlmt, ndim_chgpot ) )
    hsr_densmat = 0.0d0;  hsi_densmat = 0.0d0

    call m_ES_MagMom_To_DensMat_hsr( natm, nlmt, hsrt, hsit, hsr_densmat, hsi_densmat )
!
    do it = 1, ntyp
       do lmt1 = 1, ilmt(it)
          il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)

          do lmt2 = 1, ilmt(it)
             il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)

             if ( .not. flg_paw ) then
                if ( SpinOrbit_mode == Neglected ) then
                   if ( il1 /= il2 .or. im1 /= im2 ) cycle
                else
                   if( il1 /= il2 ) cycle
                endif
             endif

             do ia = 1, natm
                if(ityp(ia) /= it) cycle

                fac = iwei(ia)

                do is= 1, ndim_chgpot
                   c1 = real( dion0_noncl( lmt1,lmt2,is,ia ) -dsoc(lmt1,lmt2,ia,is) ) &
                        &  *hsr_densmat( ia,lmt1,lmt2,is )

! ---------------
            ! ( The results seem unchanged even if 1 is changed to 0 )
! --------
#if 1
                   c2 = aimag( dion0_noncl( lmt1,lmt2,is,ia ) ) &
                        &  *hsi_densmat( ia,lmt1,lmt2,is )
#else
!                   c2 = -aimag( dion0_noncl( lmt1,lmt2,is,ia ) ) &
!                        &  *hsi_densmat( ia,lmt1,lmt2,is )
#endif
                   c2 = 0.0d0          !!! ASMS 2015/02/20,  just in case
                   enonlc = enonlc + fac*( c1 +c2 )
                end do
             end do
          end do
       end do
    end do
!
    deallocate( hsr_densmat, hsi_densmat )

  end subroutine get_nonlocal_pot_energy_noncl

  subroutine get_nonlocal_pot_energy_noncl2  ! without SOC
    integer :: it,lmt1, lmt2, il1, im1, il2, im2, ia
    integer :: is
    real(kind=DP) :: fac
    real(kind=DP) :: c1, c2

    real(kind=DP), allocatable :: hsr_densmat( :,:,:,: )
    real(kind=DP), allocatable :: hsi_densmat( :,:,:,: )

    enonlc = 0.d0

    do it = 1, ntyp
       do lmt1 = 1, ilmt(it)
          il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
          do lmt2 = lmt1, ilmt(it)
             il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)

             if (ipaw(it)==0 ) then
                if( il1 /= il2 .or. im1 /= im2) cycle

                do ia = 1, natm
                   if(ityp(ia) /= it) cycle
                   fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0
                   enonlc = enonlc &
                        & + fac*dion(lmt1,lmt2,it)*hsrt(ia,lmt1,lmt2,1)
                end do
             else
                do ia = 1, natm
                   if(ityp(ia) /= it) cycle
                   fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0
                   Do is=1, ndim_magmom
                      enonlc = enonlc &
                           & + fac*dion_paw(lmt1,lmt2,is,ia)*hsrt(ia,lmt1,lmt2,is)
                   end do
                end do
             end if
          end do
       end do
    end do
  end subroutine get_nonlocal_pot_energy_noncl2
! ====================================================================== 11.0

  subroutine get_hartree_energy
    integer ik, i
    integer :: ist !mpi
    real(kind=DP) :: ehartr_mpi

    ehartr_mpi = 0.d0
    ehartr = 0.d0
    do ik = 1, kimg

       if(mype==0) then
          if(sw_screening_correction==ON) then
              ehartr_mpi  = ehartr_mpi + chgqt(1,ik,1)**2 * screening%phik(1)
          end if
       end if

       ist = ista_kngp
       if(ist == 1) ist = 2

       if(nspin == 1) then
          do i = ist, iend_kngp !for mpi
!!!!             ehartr_mpi  = ehartr_mpi + (chgq_l(i,ik,1)/gr_l(i))**2
             ehartr_mpi  = ehartr_mpi + PAI4 * (chgqt(i,ik,1)/gr_l(i))**2
             if(sw_screening_correction==ON) then
                 ehartr_mpi  = ehartr_mpi + chgqt(i,ik,1)**2 * screening%phik(i)
             end if
          end do
       else if(nspin == 2) then
          do i = ist, iend_kngp  !for mpi
!!!!             ehartr_mpi  = ehartr_mpi &
!!!!                  &  +((chgq_l(i,ik,UP)+chgq_l(i,ik,DOWN))/gr_l(i))**2
             ehartr_mpi  = ehartr_mpi &
                  &  + PAI4 * ((chgqt(i,ik,UP)+chgqt(i,ik,DOWN))/gr_l(i))**2
             if(sw_screening_correction==ON) then
                 ehartr_mpi  = ehartr_mpi &
                  &  + (chgqt(i,ik,UP)+chgqt(i,ik,DOWN))**2 * screening%phik(i)
             end if
          end do
       end if
    end do
    call mpi_allreduce(ehartr_mpi,ehartr,1 &
         &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
!!!!    ehartr    = univol*PAI4*ehartr*0.5d0
    ehartr    = univol*ehartr*0.5d0
  end subroutine get_hartree_energy

! ================================= added by K. Tagami ==================== 11.0
  subroutine get_hartree_energy_noncl
    integer ik, i
    integer :: ist !mpi
    real(kind=DP) :: ehartr_mpi

    ehartr_mpi = 0.d0
    ehartr = 0.d0

    do ik = 1, kimg
       ist = ista_kngp
       if(ist == 1) ist = 2

       do i = ist, iend_kngp !for mpi
           ehartr_mpi  = ehartr_mpi + (chgqt(i,ik,1)/gr_l(i))**2
       end do
    end do

    ehartr_mpi = ehartr_mpi *PAI4

    if ( sw_screening_correction == ON ) then
       do ik = 1, kimg
          ist = ista_kngp
          if(ist == 1) ist = 2

          do i = ist, iend_kngp !for mpi
             ehartr_mpi  = ehartr_mpi + chgqt(i,ik,1)**2 * screening%phik(i)
          end do
       end do
    endif

    call mpi_allreduce(ehartr_mpi,ehartr,1 &
         &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
    ehartr    = univol *ehartr *0.5d0

  end subroutine get_hartree_energy_noncl
! ========================================================================= 11.0

  subroutine get_dipole_energy(vdip_l,vext_l)
    real(kind=DP), intent(in) :: vdip_l(ista_kngp:iend_kngp,kimg)
    real(kind=DP), intent(in) :: vext_l(ista_kngp:iend_kngp,kimg)
    integer ik, i, ispin
    integer ist !mpi
    real(kind=DP) :: evdip_mpi
    real(kind=DP) :: evext_mpi

    evdip = 0.d0
    evext = 0.d0
    do ik = 1, kimg
       do ispin = 1, nspin
          if(mype==0) then
             evdip = evdip + vdip_l(1,ik)*chgqt(1,ik,ispin)
             evext = evext + vext_l(1,ik)*chgqt(1,ik,ispin)
          endif
       end do
       ist = ista_kngp
       if(ist == 1) ist = 2

       if(nspin == 1) then
          do i = ist, iend_kngp  !for mpi
             evdip = evdip + vdip_l(i,ik)*chgqt(i,ik,1)
             evext = evext + vext_l(i,ik)*chgqt(i,ik,1)
          end do
       else if(nspin == 2) then
          do i = ist, iend_kngp  !for mpi
             evdip = evdip + vdip_l(i,ik)*(chgqt(i,ik,UP)+chgqt(i,ik,DOWN))
             evext = evext + vext_l(i,ik)*(chgqt(i,ik,UP)+chgqt(i,ik,DOWN))
          end do
       end if
    end do
    if(npes > 1) then
       call mpi_allreduce(evdip,evdip_mpi,1 &
            &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       evdip = evdip_mpi
    end if
    evdip = evdip * univol
    evext = evext * univol
    edip = edip_ion + eext_ion + evdip*0.5d0 + evext
  end subroutine get_dipole_energy

! =================================== added by K. Tagami ==================== 11.0
  subroutine get_dipole_energy_noncl(vdip_l,vext_l)
    real(kind=DP), intent(in) :: vdip_l(ista_kngp:iend_kngp,kimg)
    real(kind=DP), intent(in) :: vext_l(ista_kngp:iend_kngp,kimg)
    integer ik, i, ispin
    integer ist !mpi
    real(kind=DP) :: evdip_mpi
    real(kind=DP) :: evext_mpi

    evdip = 0.d0;     evext = 0.d0
    do ik = 1, kimg
      do i = ist, iend_kngp  !for mpi
         evdip = evdip + vdip_l(i,ik)*chgqt(i,ik,1)
         evext = evext + vext_l(i,ik)*chgqt(i,ik,1)
      end do
    end do
    if(npes > 1) then
       call mpi_allreduce(evdip,evdip_mpi,1 &
            &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       evdip = evdip_mpi
    end if
    evdip = evdip * univol
    evext = evext * univol
    edip = edip_ion + eext_ion + evdip*0.5d0 + evext

  end subroutine get_dipole_energy_noncl
! ======================================================================== 11.0

! ==== KT_add === 2015/01/30
  subroutine get_kinetic_local_energy_paw(nfout)
    integer, intent(in) :: nfout
    integer :: ispin, it,lmt1, lmt2, il1, im1, il2, im2, ia
    integer :: ismax
    real(kind=DP) :: fac

    ekin_ion_paw = 0.0d0
    if ( noncol ) then
       ismax = 1
    else
       ismax = nspin
    endif
    do ispin = 1, ismax, af+1
       do it = 1, ntyp
          if(ipaw(it)==0) cycle
          do lmt1 = 1, ilmt(it)
             il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
             do lmt2 = lmt1, ilmt(it)
                il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
                if( il1 /= il2 .or. im1 /= im2) cycle
                do ia = 1, natm
                   if(ityp(ia) /= it) cycle
                   fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0
                   ekin_ion_paw = ekin_ion_paw &
                        & + fac*dion_kin_ion(lmt1,lmt2,it)*hsrt(ia,lmt1,lmt2,ispin)
                end do
             end do
          end do
       end do
    end do
    ekin_ion_paw = ekin_ion_paw*(af+1)
  end subroutine get_kinetic_local_energy_paw
! ============ 2015/01/30

  subroutine get_kinetic_energy(nfout)
    integer, intent(in) :: nfout
#ifdef __TIMER_SUB__
    call timer_sta(747)
#endif
    ekinet = eband - eohxc - eloca1 - enonlc
    if(sw_dipole_correction == ON) then
       ekinet = ekinet - evdip - evext
    end if
    if(sw_hubbard==ON) then
       ekinet = ekinet - ehub1
    end if
    if(sw_hybrid_functional == ON) then
       ekinet = ekinet - vexx
    end if
! ========================= added by K. Tagami ====================== 11.0
    if ( sw_magnetic_constraint==ON ) then
       ekinet = ekinet - emag1
    endif
! =================================================================== 11.0

! === KT_add === 2014/08/01
!!!#ifdef USE_ESPINORB
    if ( noncol ) then
       if ( flg_paw ) then
          if ( SpinOrbit_Mode == ByPawPot ) then
             ekinet = ekinet -espinorb_old
          endif
       endif
       if ( SpinOrbit_Mode == ByProjector .or. &
            &  SpinOrbit_Mode == ZeffApprox .or. &
            &  SpinOrbit_Mode == ReadFromPP ) then
          ekinet = ekinet -espinorb_now
       endif
    endif
!!!#endif
! ============== 2014/08/01

    if(ipri >= 2) then
       write(nfout,'(" !D EBAND  = ",f12.5)') eband
       write(nfout,'(" !D EOHXC  = ",f12.5)') eohxc
       write(nfout,'(" !D ELOCA1 = ",f12.5)') eloca1
       write(nfout,'(" !D ENONLC = ",f12.5)') enonlc
       if(sw_dipole_correction == ON) then
          write(nfout,'(" !D EVDIP  = ",f12.5)') evdip
          write(nfout,'(" !D EVEXT  = ",f12.5)') evext
       end if
       if(sw_hubbard==ON) write(nfout,'(" !D EHUB1  = ",f12.5)') ehub1
       if(sw_hybrid_functional == ON) write(nfout,'(" !D VEXX   = ",f12.5)') vexx
    endif
#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    if(sw_hubbard==ON) then
       write(nfout,'(" !D EBAND= ",F20.12," EOHXC= ",F20.12," ELOCA1= ",F20.12," ENONLC= ",F20.12, " EHUB1= ",F20.12)') &
            & EBAND, EOHXC, ELOCA1, ENONLC, EHUB1
    else       
       write(nfout,'(" !D EBAND= ",F20.12," EOHXC= ",F20.12," ELOCA1= ",F20.12," ENONLC= ",F20.12)') &
            & EBAND, EOHXC, ELOCA1, ENONLC
    end if
#endif
    
#ifdef __TIMER_SUB__
    call timer_end(747)
#endif
  end subroutine get_kinetic_energy

  subroutine get_kinetic_energy_directly(nfout)
     integer,intent(in) :: nfout
     integer :: ik,ib,ig, is
     real(kind=DP), allocatable, dimension(:)       :: ekin
     real(kind=DP) :: ekinet_tmp,ekinet_mpi

     allocate(ekin(kg1));ekin=0.d0
     ekinet = 0.d0

     if ( noncol ) then
        do ik=1, kv3, ndim_spinor
           call m_pwBS_kinetic_energies(ik,vkxyz,ekin)
           if(map_k(ik)/=myrank_k)cycle
           do ib=1,np_e
              ekinet_tmp=0.d0
              do ig=1,iba(ik)
                 ekinet_tmp = ekinet_tmp +2.d0*ekin(ig) &
                      &  *(zaj_l(ig,ib,ik,1)**2 +zaj_l(ig,ib,ik,2)**2 &
                      &   +zaj_l(ig,ib,ik+1,1)**2 +zaj_l(ig,ib,ik+1,2)**2 )
              enddo
!              ekinet = ekinet+occup_l(ib,ik)/dble(kv3/ndim_spinor)*ekinet_tmp
              ekinet = ekinet+occup_l(ib,ik)/dble(kv3)*ekinet_tmp
           enddo
        enddo
     else
        do is = ista_spin, iend_spin
        !do ik=1,kv3,af+1
         do ik = is, kv3-nspin+is, nspin
           call m_pwBS_kinetic_energies(ik,vkxyz,ekin)
           if(map_k(ik)/=myrank_k)cycle
           do ib=1,np_e
              ekinet_tmp=0.d0
              if (kimg==1) then
                 do ig=1,iba(ik)
                    ekinet_tmp = ekinet_tmp+2.d0*ekin(ig)*zaj_l(ig,ib,ik,1)**2
                 enddo
              else if (kimg==2) then
                 do ig=1,iba(ik)
                    ekinet_tmp = ekinet_tmp+2.d0*ekin(ig)*(zaj_l(ig,ib,ik,1)**2 &
                         & +zaj_l(ig,ib,ik,2)**2)
                 enddo
              endif
              if(k_symmetry(ik)==GAMMA) ekinet_tmp = ekinet_tmp*2.d0
              ekinet = ekinet+occup_l(ib,ik)/dble(kv3)*ekinet_tmp
           enddo
         enddo
        enddo
        if(af==1) ekinet = ekinet*2.d0
     endif
     ekinet_mpi = 0.d0
     call mpi_allreduce(ekinet,ekinet_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
     ekinet = ekinet_mpi
     deallocate(ekin)
  end subroutine get_kinetic_energy_directly

! ==== positron ===
 subroutine get_xc_and_HE_of_old_CD_pztr(vepc_l)
    real(kind=DP), intent(in) :: vepc_l(ista_kngp:iend_kngp,kimg,nspin)
    integer ik, i, ispin
    integer ist !mpi
    real(kind=DP) :: eohxc_pztr_mpi

    eohxc_pztr = 0.d0
    do ik = 1, kimg
       do ispin = 1, nspin
          if(mype==0) then
             eohxc_pztr = eohxc_pztr + vepc_l(1,ik,ispin)*pchg_l(1,ik)
          endif
       end do
       ist = ista_kngp
       if(ist == 1) ist = 2

       if(nspin == 1) then
          do i = ist, iend_kngp  !for mpi
             eohxc_pztr = eohxc_pztr &
                  & + (vepc_l(i,ik,1) -PAI4*chgqto(i,ik,1)/gr_l(i)**2) &
                  &   * pchg_l(i,ik)
          end do
       else if(nspin == 2) then
          call phase_error_with_msg(6,"YYY",__LINE__,__FILE__)
       end if
    end do
    if(npes > 1) then
       call mpi_allreduce(eohxc_pztr,eohxc_pztr_mpi,1 &
            &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       eohxc_pztr = eohxc_pztr_mpi
    end if
!
    eohxc_pztr = univol *eohxc_pztr
  end subroutine get_xc_and_HE_of_old_CD_pztr

  subroutine get_kinetic_energy_pztr(nfout)
    integer,intent(in) :: nfout
    integer :: ik,ib,ig
    real(kind=DP) :: ekinet_tmp, ekinet_mpi

    ekin_pztr = 0.0
    Do ib=1, npeg
       if ( nprvf_ordr(ib) /= 1 ) cycle
       ekin_pztr = ekin_pztr + pev(ib)
    End do
    ekin_pztr = ekin_pztr  -eohxc_pztr -elocal_pztr
  end subroutine get_kinetic_energy_pztr

  subroutine get_kinetic_energy_pztr_direct(nfout)
    integer,intent(in) :: nfout
    integer :: ik,ib,ig
    real(kind=DP), allocatable, dimension(:)       :: ekin
    real(kind=DP) :: ekinet_tmp, ekinet_mpi, e1

    allocate(ekin(kg1_pwf));ekin=0.d0

    ekin_pztr = 0.d0

    call m_pwBS_pstrn_kinetic_energies(ekin)

    do ib=1, npeg
       if ( nprvf_ordr(ib) /= 1 ) cycle
       ekinet_tmp=0.d0
       if (kimg==1) then
          do ig=1, kg1_pwf
             ekinet_tmp = ekinet_tmp + ekin(ig)*pzaj(ig,ib,1)**2
          enddo
       else if (kimg==2) then
          do ig=1, kg1_pwf
             ekinet_tmp = ekinet_tmp + ekin(ig)*( pzaj(ig,ib,1)**2 &
                  &                              +pzaj(ig,ib,2)**2)
          enddo
       endif
       ekin_pztr = ekin_pztr +ekinet_tmp
    enddo
    deallocate(ekin)

  end subroutine get_kinetic_energy_pztr_direct

  subroutine get_hartree_energy_inter_ep
    integer ik, i
    integer :: ist !mpi
    real(kind=DP) :: ehartr_mpi

    ehartr_mpi = 0.d0
    ehartr_ep = 0.d0

    do ik = 1, kimg

       ist = ista_kngp
       if(ist == 1) ist = 2

       if(nspin == 1) then
          do i = ist, iend_kngp !for mpi
             ehartr_mpi  = ehartr_mpi + PAI4 *chgqt(i,ik,1) *pchg_l(i,ik) &
                  &                          /gr_l(i)**2
          end do
       else if(nspin == 2) then
          call phase_error_with_msg(6,"PPP",__LINE__,__FILE__)
       endif
    end do
    call mpi_allreduce(ehartr_mpi,ehartr_ep,1 &
         &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)

    ehartr_ep    = -univol *ehartr_ep

  end subroutine get_hartree_energy_inter_ep

  subroutine get_local_pot_energy_pztr
    integer       :: ik, it, i, ig
    real(kind=DP) :: elocal_mpi

    elocal_pztr = 0.d0
    do ik = 1, kimg
       if(nspin == 1) then
          do it = 1,ntyp
             do i = ista_kngp, iend_kngp !for mpi
                elocal_pztr = elocal_pztr &
                     & + psc_l(i,it)*zfm3_l(i,it,ik)*pchg_l(i,ik)
             end do
          end do
       else
          call phase_error_with_msg(6,"PPP",__LINE__,__FILE__)
       endif
    end do
    if(npes > 1) then
       call mpi_allreduce(elocal_pztr,elocal_mpi,1 &
            &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       elocal_pztr = elocal_mpi
    end if

    elocal_pztr    = -univol*elocal_pztr
  end subroutine get_local_pot_energy_pztr
! =====

  subroutine get_xc_and_HE_of_old_CD_paw(nfout)
    integer, intent(in) :: nfout
    integer :: ispin, it,lmt1, lmt2, il1, im1, il2, im2, ia
    integer, parameter :: DEBUGPRINTLEVEL = 2
    real(kind=DP) :: fac
    eohxc_paw = 0.d0
    if(ipri>= DEBUGPRINTLEVEL) write(nfout,'(" ipaw = ",10i8)') (ipaw(it),it=1,ntyp)
    do ispin = 1, nspin, af+1
       do it = 1, ntyp
          if(ipri>= DEBUGPRINTLEVEL) then
             write(nfout,'(" it = ",i8," <<get_xc_and_HE_of_old_CD_paw>>")') it
             write(nfout,'(" ilmt(",i3,") = ",i8)') it, ilmt(it)
             do ia = 1, natm
                if(ityp(ia) /= it) cycle
                write(nfout,'(" -- dion_hartree, dion_vxc, hsr --, ispin = ",i3, " ia = ", i8," ipaw(",i3,") = ",i2)') &
                     & ispin,ia, it,ipaw(it)
                do lmt1 = 1, ilmt(it)
                   il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
                   do lmt2 = lmt1, ilmt(it)
                      il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
!!!!                      if( il1 /= il2 .or. im1 /= im2) cycle
                      write(nfout,'(3f20.8)') dion_hartree(lmt1,lmt2,ia), dion_vxc(lmt1,lmt2,ispin,ia) &
                     &, hsrt(ia,lmt1,lmt2,ispin)
                   end do
                end do
             end do
          end if

          if(ipaw(it)==0) cycle
          do lmt1 = 1, ilmt(it)
             il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
             do lmt2 = lmt1, ilmt(it)
                il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
!!!                if( il1 /= il2 .or. im1 /= im2) cycle     !! ASMS 2015/02/06
!xocl spread do/ind_katm
                do ia = 1, natm
                   if(ityp(ia) /= it) cycle
                   fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0
                   eohxc_paw = eohxc_paw &
                        & + fac*(dion_hartree(lmt1,lmt2,ia)+ &
                                dion_vxc(lmt1,lmt2,ispin,ia))*hsrt(ia,lmt1,lmt2,ispin)
                end do
!xocl end spread
             end do
          end do
       end do
    end do
    eohxc_paw = eohxc_paw*(af+1)
    if(ipri>=DEBUGPRINTLEVEL) write(nfout,'(" eohxc_paw = ",f20.8)') eohxc_paw
  end subroutine get_xc_and_HE_of_old_CD_paw

! ======================= added by K. Tagami ================ 11.0
  subroutine get_xc_and_HE_old_CD_paw_noncl(nfout)
    integer, intent(in) :: nfout
    integer :: ispin, it,lmt1, lmt2, il1, im1, il2, im2, ia
    integer, parameter :: DEBUGPRINTLEVEL = 2
    real(kind=DP) :: fac

    eohxc_paw = 0.d0
    if(ipri>= DEBUGPRINTLEVEL) write(nfout,'(" ipaw = ",10i8)') (ipaw(it),it=1,ntyp)

    do ispin = 1, ndim_magmom

       do it = 1, ntyp

          if(ipri>= DEBUGPRINTLEVEL) then
             write(nfout,'(" it = ",i8," <<get_xc_and_HE_of_old_CD_paw>>")') it
             write(nfout,'(" ilmt(",i3,") = ",i8)') it, ilmt(it)
             do ia = 1, natm
                if(ityp(ia) /= it) cycle
                write(nfout,'(" -- dion_hartree, dion_vxc, hsr --, ispin = ",i3, " ia = ", i8," ipaw(",i3,") = ",i2)') &
                     & ispin,ia, it,ipaw(it)
                do lmt1 = 1, ilmt(it)
                   il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
                   do lmt2 = lmt1, ilmt(it)
                      il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
!!!!                      if( il1 /= il2 .or. im1 /= im2) cycle
                      write(nfout,'(3f20.8)') dion_hartree(lmt1,lmt2,ia), dion_vxc(lmt1,lmt2,ispin,ia) &
                     & , hsrt(ia,lmt1,lmt2,ispin)
                   end do
                end do
             end do
          end if

          if (ipaw(it)==0) cycle
          do lmt1 = 1, ilmt(it)
             il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
             do lmt2 = lmt1, ilmt(it)
                il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
!!!                if( il1 /= il2 .or. im1 /= im2) cycle       !! ASMS 2015/02/06
!xocl spread do/ind_katm
                do ia = 1, natm
                   if(ityp(ia) /= it) cycle
                   fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0

                   if ( ispin == 1 ) then
                      eohxc_paw = eohxc_paw &
                           &     + fac *( dion_hartree(lmt1,lmt2,ia) &
                           &             +dion_vxc(lmt1,lmt2,ispin,ia ) ) &
                           &           *hsrt(ia,lmt1,lmt2,ispin)
                   else
                      eohxc_paw = eohxc_paw &
                           &     + fac * dion_vxc(lmt1,lmt2,ispin,ia ) &
                           &           *hsrt(ia,lmt1,lmt2,ispin)

                   endif

                end do
!xocl end spread
             end do
          end do
       end do
    end do

    if(ipri>=DEBUGPRINTLEVEL) write(nfout,'(" eohxc_paw = ",f20.8)') eohxc_paw
  end subroutine get_xc_and_HE_old_CD_paw_noncl
! ====================================================================== 11.0

  subroutine get_xc_and_HE_of_old_CD_paw_sym
    integer :: ispin, it,lmt1, lmt2, il1, im1, il2, im2, ia
    real(kind=DP) :: fac
    integer :: iopr,ja
    eohxc_paw = 0.d0

    do iopr=1,nopr

    do ispin = 1, nspin, af+1
       do it = 1, ntyp
          if(ipaw(it)==0) cycle
          do lmt1 = 1, ilmt(it)
             il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
             do lmt2 = lmt1, ilmt(it)
                il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
!!!                if( il1 /= il2 .or. im1 /= im2) cycle     ! ASMS 2015/02/06
!xocl spread do/ind_katm
                do ia = 1, natm
                   if(ityp(ia) /= it) cycle
                   ja=abs(ia2ia_symmtry_op(ia,iopr))
                   fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0
                   eohxc_paw = eohxc_paw &
                        & + fac*(dion_hartree(lmt1,lmt2,ia)+ &
                                dion_vxc(lmt1,lmt2,ispin,ja))*hsrt(ja,lmt1,lmt2,ispin)
                end do
!xocl end spread
             end do
          end do
       end do
    end do

    end do
    eohxc_paw = eohxc_paw*(af+1)/dble(nopr)
  end subroutine get_xc_and_HE_of_old_CD_paw_sym

! ================================= added by K. Tagami ============= 11.0
  subroutine get_xc_HE_old_CD_paw_sym_noncl
    integer :: ispin, it,lmt1, lmt2, il1, im1, il2, im2, ia
    real(kind=DP) :: fac
    integer :: iopr,ja
    eohxc_paw = 0.d0

    do iopr=1,nopr

       do ispin = 1, ndim_magmom
          do it = 1, ntyp
             if(ipaw(it)==0) cycle

             do lmt1 = 1, ilmt(it)
                il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
                do lmt2 = lmt1, ilmt(it)
                   il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
!!!                   if( il1 /= il2 .or. im1 /= im2) cycle     !! ASMS 2015/02/06
!xocl spread do/ind_katm
                   do ia = 1, natm
                      if(ityp(ia) /= it) cycle
                      ja=abs(ia2ia_symmtry_op(ia,iopr))
                      fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0

                      if ( ispin==1 ) then
                         eohxc_paw = eohxc_paw &
                              &     + fac *( dion_hartree(lmt1,lmt2,ia) &
                              &             +dion_vxc(lmt1,lmt2,ispin,ja) ) &
                              &           *hsrt(ja,lmt1,lmt2,ispin)
                      else
                         eohxc_paw = eohxc_paw &
                              &     + fac *dion_vxc(lmt1,lmt2,ispin,ja)  &
                              &           *hsrt(ja,lmt1,lmt2,ispin)
                      endif
                   end do
!xocl end spread
                end do
             end do
          end do
       end do

    end do
    eohxc_paw = eohxc_paw /dble(nopr)

  end subroutine get_xc_HE_old_CD_paw_sym_noncl
! ============================================================= 11.0

  subroutine get_hartree_energy_paw(nfout)
    integer, intent(in) :: nfout
    integer :: ispin, it,lmt1, lmt2, il1, im1, il2, im2, ia
    integer, parameter :: DEBUGPRINTLEVEL = 2
    real(kind=DP) :: fac

    if(ipri>=DEBUGPRINTLEVEL) then
       write(nfout,'(" -- iwei -- <<get_hartree_energy_paw>>")')
       write(nfout,'(10i8)') (iwei(ia),ia=1,natm)
       write(nfout,'(" -- ia, lmt1, lmt2, dion_hartree_now, -- hsr --")')
       do ia = 1, natm
          it = ityp(ia)
          do lmt1 = 1, ilmt(it)
             il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
             do lmt2 = lmt1, ilmt(it)
                il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
!!!                if( il1 /= il2 .or. im1 /= im2) cycle
                write(nfout,'(3i6,3d20.8)') ia, lmt1, lmt2, dion_hartree_now(lmt1,lmt2,ia) &
                     & , (hsrt(ia,lmt1,lmt2,ispin),ispin = 1, nspin/(af+1))
             end do
          end do
       end do
    end if

    ehartr_paw = 0.d0
    do ispin = 1, nspin, af+1
       do it = 1, ntyp
          if(ipaw(it)==0) cycle
          do lmt1 = 1, ilmt(it)
             il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
             do lmt2 = lmt1, ilmt(it)
                il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
!!!                if( il1 /= il2 .or. im1 /= im2) cycle    !! ASMS 2015/02/06
!xocl spread do/ind_katm
                do ia = 1, natm
                   if(ityp(ia) /= it) cycle
                   fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0
                   ehartr_paw = ehartr_paw &
                        & + fac*dion_hartree_now(lmt1,lmt2,ia)*hsrt(ia,lmt1,lmt2,ispin)
                end do
!xocl end spread
             end do
          end do
       end do
    end do
    ehartr_paw = ehartr_paw*(af+1)*0.5d0
  end subroutine get_hartree_energy_paw

! =========================== added by K. Tagami =================== 11.0
  subroutine get_hartree_energy_paw_noncl(nfout)
    integer, intent(in) :: nfout
    integer :: ispin, it,lmt1, lmt2, il1, im1, il2, im2, ia
    integer, parameter :: DEBUGPRINTLEVEL = 2
    real(kind=DP) :: fac

    if (ipri>=DEBUGPRINTLEVEL) then
       write(nfout,'(" -- iwei -- <<get_hartree_energy_paw_noncl>>")')
       write(nfout,'(10i8)') (iwei(ia),ia=1,natm)
       write(nfout,'(" -- ia, lmt1, lmt2, dion_hartree_now, -- hsr --")')
       do ia = 1, natm
          it = ityp(ia)
          do lmt1 = 1, ilmt(it)
             il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
             do lmt2 = lmt1, ilmt(it)
                il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
!!!                if( il1 /= il2 .or. im1 /= im2) cycle
                write(nfout,*) ia, lmt1, lmt2, dion_hartree_now(lmt1,lmt2,ia) &
                     & , (hsrt(ia,lmt1,lmt2,ispin),ispin = 1, ndim_magmom)
             end do
          end do
       end do
    end if

    ehartr_paw = 0.d0
    do ispin = 1, 1
       do it = 1, ntyp
          if(ipaw(it)==0) cycle
          do lmt1 = 1, ilmt(it)
             il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
             do lmt2 = lmt1, ilmt(it)
                il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
!!!                if( il1 /= il2 .or. im1 /= im2) cycle     !! ASMS 2015/02/06
!xocl spread do/ind_katm
                do ia = 1, natm
                   if(ityp(ia) /= it) cycle
                   fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0
                   ehartr_paw = ehartr_paw &
                        & + fac*dion_hartree_now(lmt1,lmt2,ia)*hsrt(ia,lmt1,lmt2,ispin)
                end do
!xocl end spread
             end do
          end do
       end do
    end do
    ehartr_paw = ehartr_paw *0.5d0

  end subroutine get_hartree_energy_paw_noncl
! ================================================================= 11.0

  subroutine get_hartree_energy_paw_sym
    integer :: ispin, it,lmt1, lmt2, il1, im1, il2, im2, ia
    real(kind=DP) :: fac
    integer :: iopr,ja

    ehartr_paw = 0.d0

    do iopr=1,nopr

    do ispin = 1, nspin, af+1
       do it = 1, ntyp
          if(ipaw(it)==0) cycle
          do lmt1 = 1, ilmt(it)
             il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
             do lmt2 = lmt1, ilmt(it)
                il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
!!!                if( il1 /= il2 .or. im1 /= im2) cycle     !! ASMS 2015/02/06
!xocl spread do/ind_katm
                do ia = 1, natm
                   if(ityp(ia) /= it) cycle
                   ja=abs(ia2ia_symmtry_op(ia,iopr))
                   fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0
                   ehartr_paw = ehartr_paw &
                        & + fac*dion_hartree_now(lmt1,lmt2,ja)*hsrt(ja,lmt1,lmt2,ispin)
!print *,fac*dion_hartree_now(lmt1,lmt2,ja)*hsr(ja,lmt1,lmt2,ispin)/det,det
                end do
!xocl end spread
             end do
          end do
       end do
    end do

    end do
    ehartr_paw = ehartr_paw*(af+1)*0.5d0/dble(nopr)
  end subroutine get_hartree_energy_paw_sym

! ============================== added by K. Tagami ============== 11.0
  subroutine get_hartree_ene_paw_sym_noncl
    integer :: ispin, it,lmt1, lmt2, il1, im1, il2, im2, ia
    real(kind=DP) :: fac
    integer :: iopr,ja

    ehartr_paw = 0.d0

    do iopr=1,nopr

       do ispin = 1, 1
          do it = 1, ntyp
             if(ipaw(it)==0) cycle
             do lmt1 = 1, ilmt(it)
                il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
                do lmt2 = lmt1, ilmt(it)
                   il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
!!!!                   if( il1 /= il2 .or. im1 /= im2) cycle      !! ASMS 2015/02/06
!xocl spread do/ind_katm
                   do ia = 1, natm
                      if(ityp(ia) /= it) cycle
                      ja=abs(ia2ia_symmtry_op(ia,iopr))
                      fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0
                      ehartr_paw = ehartr_paw &
                           & + fac*dion_hartree_now(lmt1,lmt2,ja)*hsrt(ja,lmt1,lmt2,ispin)
                   end do
!xocl end spread
                end do
             end do
          end do
       end do

    end do
    ehartr_paw = ehartr_paw*0.5d0/dble(nopr)

  end subroutine get_hartree_ene_paw_sym_noncl
! ================================================================= 11.0

  subroutine get_entropic_term(nfout)
    integer, intent(in) :: nfout
#ifdef __TIMER_SUB__
    call timer_sta(748)
#endif
    eentropy = -width*band_entropy
    if(ipri >= 2) then
       write(nfout,'(" !D EENTROPY = ",f20.14)') eentropy
    endif
#ifdef __TIMER_SUB__
    call timer_end(748)
#endif
  end subroutine get_entropic_term

  subroutine get_hubbard_energy(nfout)
    integer, intent(in) :: nfout
    ehub0 = 0.d0
    ehub1 = 0.d0
    if(sw_hubbard==ON) then
      call m_Hubbard_energy(ehub0)
      call get_hubbard_potential_energy(ehub1)
    end if
    if(ipri >= 2) then
       write(nfout,'(" !D EHUB0 = ",f20.14)') ehub0
       write(nfout,'(" !D EHUB1 = ",f20.14)') ehub1
    endif
  end subroutine get_hubbard_energy

! ============================= added by K. Tagami =============== 11.0
  subroutine get_hubbard_energy_noncl(nfout)
    integer, intent(in) :: nfout

    if ( sw_hubbard == on ) then
!!       write(*,*) 'Not supported : DFT-U '

       call m_Hubbard_energy_noncl(ehub0)

!       call m_Hubbard_energy2_noncl(ehub1)
!
!       write(*,*) 'ehub0 ehu1 comp = ', ehub0, ehub1
!
!       call m_Hubbard_energy3_noncl(ehub1)
!       write(*,*) 'ehub3 = ', ehub1
!       stop

       call get_hubbard_pot_energy_noncl3(ehub1)
    end if
    if(ipri >= 2) then
       write(nfout,'(" !D EHUB0 = ",f20.14)') ehub0
       write(nfout,'(" !D EHUB1 = ",f20.14)') ehub1
    endif
  end subroutine get_hubbard_energy_noncl
! ================================================================ 11.0

  subroutine get_hubbard_potential_energy(ehub1)
    real(kind=DP), intent(out) :: ehub1
    integer :: ispin, it,lmt1, lmt2, il1, im1, il2, im2, ia
    integer :: ih, l1p
    real(kind=DP) :: fac
    ehub1 = 0.d0
    do ispin = 1, nspin, af+1
       do ia = 1, natm
          ih = ihubbard(ia)
          if(ih == 0) cycle
          it = ityp(ia)
          l1p = proj_attribute(ih)%l+1
          do lmt1 = 1, ilmt(it)
             il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
             if(il1 /= l1p) cycle
             do lmt2 = lmt1, ilmt(it)
                il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
                if(il2 /= l1p) cycle
                !!$if(im1 /= im2) cycle
                fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0
                ehub1 = ehub1 &
                      & + fac*dhub(lmt1,lmt2,ia,ispin)*hsrt(ia,lmt1,lmt2,ispin)
             end do
          end do
       end do
    end do
    ehub1 = ehub1*(af+1)
  end subroutine get_hubbard_potential_energy

! ================================ added by K. Tagami ================= 11.0
  subroutine get_hubbard_pot_energy_noncl(ehub1)
    real(kind=DP), intent(out) :: ehub1
    integer :: ispin, it,lmt1, lmt2, il1, im1, il2, im2, ia
    integer :: ih, l1p

    integer :: is
    real(kind=DP) :: fac

    ehub1 = 0.d0

    do is=1, ndim_magmom
       do ia = 1, natm
          ih = ihubbard(ia)
          if(ih == 0) cycle
          it = ityp(ia)
          l1p = proj_attribute(ih)%l+1
          do lmt1 = 1, ilmt(it)
             il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
             if(il1 /= l1p) cycle
             do lmt2 = lmt1, ilmt(it)
                il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
                if(il2 /= l1p) cycle
                !!$if(im1 /= im2) cycle
                fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0
                ehub1 = ehub1 &
                      & + fac*dhub(lmt1,lmt2,ia,is)*hsrt(ia,lmt1,lmt2,is)
             end do
          end do
       end do
    end do

  end subroutine get_hubbard_pot_energy_noncl

  subroutine get_hubbard_pot_energy_noncl2(ehub1)
    real(kind=DP), intent(out) :: ehub1
    integer :: ispin, it,lmt1, lmt2, il1, im1, il2, im2, ia
    integer :: ih, l1p

    integer :: is
    real(kind=DP) :: fac, c1, c2

    ehub1 = 0.d0

    do is=1, ndim_magmom
       do ia = 1, natm
          ih = ihubbard(ia)
          if(ih == 0) cycle
          it = ityp(ia)
          l1p = proj_attribute(ih)%l+1
          do lmt1 = 1, ilmt(it)
             il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
             if(il1 /= l1p) cycle

             do lmt2 = 1, ilmt(it)
                il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
                if(il2 /= l1p) cycle
                !!$if(im1 /= im2) cycle
                fac = iwei(ia)
!
                c1 = dhub( lmt1,lmt2,ia,is ) &
                     &  *hsrt( ia,lmt1,lmt2,is )
                c2 = -dhub_aimag( lmt1,lmt2,ia,is ) &
                     &  *hsit( ia,lmt1,lmt2,is )

                ehub1 = ehub1 + fac*( c1 +c2 )
             end do
          end do
       end do
    end do

  end subroutine get_hubbard_pot_energy_noncl2

  subroutine get_hubbard_pot_energy_noncl3(ehub1)
    real(kind=DP), intent(out) :: ehub1
    integer :: ispin, it,lmt1, lmt2, il1, im1, il2, im2, ia
    integer :: ih, l1p

    integer :: is
    real(kind=DP) :: fac, c1, c2
    real(kind=DP), allocatable :: hsr_densmat( :,:,:,: )
    real(kind=DP), allocatable :: hsi_densmat( :,:,:,: )
    complex(kind=CMPLDP), allocatable :: dhub_ssrep(:,:,:,:)

    ehub1 = 0.d0

    allocate( dhub_ssrep(nlmt,nlmt,natm,ndim_magmom) )
    dhub_ssrep = 0.0d0

    allocate( hsr_densmat( natm, nlmt, nlmt, ndim_chgpot ) )
    allocate( hsi_densmat( natm, nlmt, nlmt, ndim_chgpot ) )
    hsr_densmat = 0.0d0;  hsi_densmat = 0.0d0

    call m_ES_MagMom_to_DensMat_Dhub( dhub, dhub_aimag, dhub_ssrep )
    call m_ES_MagMom_To_DensMat_hsr( natm, nlmt, hsrt, hsit, hsr_densmat, hsi_densmat )

    do ia = 1, natm
       ih = ihubbard(ia)
       if(ih == 0) cycle
       it = ityp(ia)
       l1p = proj_attribute(ih)%l+1
       do lmt1 = 1, ilmt(it)
          il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
          if(il1 /= l1p) cycle

          do lmt2 = 1, ilmt(it)
             il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
             if(il2 /= l1p) cycle
                !!$if(im1 /= im2) cycle
             fac = iwei(ia)
!
             Do is=1, ndim_chgpot
                c1 = real(  dhub_ssrep( lmt1,lmt2,ia,is ) ) &
                     &  *hsr_densmat( ia,lmt1,lmt2,is )
                c2 = aimag( dhub_ssrep( lmt1,lmt2,ia,is ) ) &
                     &  *hsi_densmat( ia,lmt1,lmt2,is )
!                c2 = -aimag( dhub_ssrep( lmt1,lmt2,ia,is ) ) &
!                     &  *hsi_densmat( ia,lmt1,lmt2,is )
                ehub1 = ehub1 + fac*( c1 +c2 )
             end do
          end do
       end do
    end do

    deallocate( hsr_densmat, hsi_densmat )
    deallocate( dhub_ssrep )

  end subroutine get_hubbard_pot_energy_noncl3
! ========================================================================== 11.0

! ============================ added by K. Tagami ========================= 11.0
  subroutine get_magnetic_constraint_energy
    if ( sw_magnetic_constraint == ON ) then
       call m_ES_calc_MagConstraint_Energy( emag1, emag0 )
    endif
  end subroutine get_magnetic_constraint_energy
! ========================================================================= 11.0

! ======================== added by K. Tagami =============== 11.0
  subroutine get_spinorbit_energy_noncl2( ene_out )
    real(kind=DP) :: ene_out

    integer :: lmt1, lmt2, il1, im1, il2, im2, it, ia
    integer :: p, q, is1, is2, istmp, ib, ik
    integer :: ierr

    real(kind=DP) :: fac
    real(kind=DP) :: c1, c2, d_factor, csum, w_n

    complex(kind=CMPLDP) :: z1, z2, ztmp, zsum

! ------------------------
    ene_out = 0.d0

    d_factor = 1.d0/ dble( kv3 /ndim_spinor )
    csum = 0.0d0
    zsum = 0.0d0
! ----------------------------------

    do ia = 1, natm
       it = ityp(ia)

       do lmt1 = 1, ilmt(it)
          il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
          p = lmta(lmt1,ia)

          do lmt2 = 1, ilmt(it)
             il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
             q = lmta(lmt2,ia)

             Do is1=1, ndim_spinor
                Do is2=1, ndim_spinor
                   istmp = ( is1 -1 )*ndim_spinor +is2

                   Do ik=1, kv3, ndim_spinor
                      if ( map_k(ik) /= myrank_k ) cycle

                      Do ib=1, np_e
                         w_n = occup_l(ib,ik) *d_factor

                         z1 = cmplx( fsr_l(ib,p,ik+is1-1), fsi_l(ib,p,ik+is1-1) )
                         z2 = cmplx( fsr_l(ib,q,ik+is2-1), fsi_l(ib,q,ik+is2-1) )

#if 1
! -- org
                         csum = csum + w_n *conjg(z1) *dsoc(lmt1,lmt2,ia,istmp) *z2
#else
! -- dbug
!                         csum = csum + w_n *z1 *dsoc(lmt1,lmt2,ia,istmp) *conjg(z2)
#endif
                      End do
                   End do
                end do
             end do
          end do

       End do

    End do

! ---------------
    if ( npes > 1 ) then
       call mpi_allreduce( csum, ene_out, 1, mpi_double_precision, &
            &              mpi_sum, MPI_CommGroup, ierr )
    else
       ene_out = csum
    endif
! -------------

  end subroutine get_spinorbit_energy_noncl2

  subroutine get_spinorbit_energy_noncl3( ene_out )
    real(kind=DP) :: ene_out

    integer :: is, it,lmt1, lmt2, il1, im1, il2, im2, ia
    real(kind=DP) :: fac
    real(kind=DP) :: c1, c2

    real(kind=DP), allocatable :: hsr_densmat( :,:,:,: )
    real(kind=DP), allocatable :: hsi_densmat( :,:,:,: )

    ene_out = 0.d0

! --
!   hsr (magomom) --> hsr, hs ( chgopt)

    allocate( hsr_densmat( natm, nlmt, nlmt, ndim_chgpot ) )
    allocate( hsi_densmat( natm, nlmt, nlmt, ndim_chgpot ) )
    hsr_densmat = 0.0d0;  hsi_densmat = 0.0d0

    call m_ES_MagMom_To_DensMat_hsr( natm, nlmt, hsrt, hsit, hsr_densmat, hsi_densmat )
! ----------
    do ia = 1, natm
       it = ityp(ia)

       do lmt1 = 1, ilmt(it)
          il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
!          p = lmta(lmt1,ia)

          do lmt2 = 1, ilmt(it)
             il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
!             q = lmta(lmt2,ia)

             fac = iwei(ia)

             do is= 1, ndim_chgpot
                c1 = real(  dsoc( lmt1,lmt2,ia,is ) ) &
                     &  *hsr_densmat( ia,lmt1,lmt2,is )
                c2 = aimag( dsoc( lmt1,lmt2,ia,is ) ) &
                     &  *hsi_densmat( ia,lmt1,lmt2,is )
                ene_out = ene_out + fac*( c1 +c2 )
             end do
          end do
       end do
    end do

    deallocate( hsr_densmat, hsi_densmat )

  end subroutine get_spinorbit_energy_noncl3
! =========================================================== 11.0

  subroutine sumup_all_energies(nfout,exc,display_on)
    integer, intent(in) :: nfout
    logical, intent(in) :: display_on
    real(kind=DP), intent(in) :: exc
    integer :: i

    real(kind=DP) :: edel, EPC_t
    character(len("potential_mixing")) :: tag_mixing
    character(len=2) ::ndecimals    ! T. Yamasaki, 2024/03/22
    integer :: ij                   ! T. Yamasaki, 2024/03/22
    character(256) :: fmt,fmt2      ! T. Yamasaki, 2024/03/22
    integer :: iterdigit            ! T. Yamasaki, 2024/03/22
    character(len=2) :: niterdigits ! T. Yamasaki, 2024/03/22
#ifdef DEBUG_ITERATION_WRITE
    integer :: iteration_tmp        ! T. Yamasaki, 2024/03/22
#endif    

#ifdef __TIMER_SUB__
    call timer_sta(749)
#endif

! ===== KT_mod ==== 13.0S
!    etotal0 = ekinet+ehartr+exc+elocal+enonlc+eewald-epc
!    if(flg_paw) etotal0=etotal0-eohxc_paw+ehartr_paw+exc_ae-exc_ps
!
    if ( flg_paw ) then
#if 0
       ekin_ion_paw = enonlc -eohxc_paw      !! ??
#else
       eohxc_paw = enonlc -ekin_ion_paw      ! meaningless, this is just for output
#endif
       etotal0 = ekinet +ehartr +exc +elocal +ekin_ion_paw +eewald -epc_paw
       etotal0 = etotal0 +ehartr_paw +exc_ae -exc_ps
    else
       etotal0 = ekinet+ehartr+exc+elocal+enonlc+eewald-epc
    endif
! ================= 13.0S

! ==== POSITRON SCF == 2015/11/28
    if ( sw_positron /= OFF ) then
       if ( positron_method == Positron_GGGC ) then
          etotal0 = etotal0 +ekin_pztr +elocal_pztr +ecorr_pztr +ehartr_ep
       endif
    endif
! =================== 2015/11/28

    if(sw_dipole_correction == ON) etotal0 = etotal0 + edip
    if(sw_hubbard == ON) etotal0 = etotal0 + ehub0
    if(sw_hybrid_functional == ON) etotal0 = etotal0 + eexx
    if(sw_fef == ON) etotal0 = etotal0 + eplr
    if(ntyp_vdw>0 .or. sw_vdw_correction == ON .and. vdw_method == VDW_DFTD3) etotal0 = etotal0 + evdw
#ifndef DISABLE_VDWDF
    if ( xctype == "vdwdf" ) etotal0 = etotal0 +ecnl_vdwdf
#endif

! =============================== added by K. Tagami =================== 11.0
    if ( sw_magnetic_constraint == ON ) etotal0 = etotal0 + emag0
! ====================================================================== 11.0

! =============================== added by K. Tagami =================== 11.0
!!!#ifdef USE_ESPINORB
    if ( noncol ) then
       if ( SpinOrbit_Mode /= Neglected ) then
          etotal0 = etotal0 + espinorb_now
       endif
    endif
!!!#endif
! ====================================================================== 11.0

    if(num_regions>0)then
       do i=1,num_regions
          if(regions(i)%tally) etotal0 = etotal0+regions(i)%energy
       enddo
    endif

    etotal = etotal0+eentropy
    if(way_of_smearing == COLD) etotal0 = (2.d0*etotal+etotal0)/3.d0
    if(way_of_smearing == PARABOLIC) etotal0 = (etotal+etotal0)*0.5d0

! =============== KT_add ============================================= 13.0E
    if(way_of_smearing == Fermi_Dirac) etotal0 = (etotal+etotal0)*0.5d0
! ==================================================================== 13.0E

    if(way_of_smearing == MP)  &
    &  etotal0 = (1.d0/(dble(order_mp+2)))*(dble(order_mp+1)*etotal+etotal0)

#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    if(display_on) then
#else
    if(display_on .and. ipri >= 1) then
#endif
       edel = cal_edeltb()
       if(number_of_cdmixing_applied<=0)  call m_CtrlP_push_CDMixingNameApplied(" ",1)
       tag_mixing = "Charge-Mixing" ;      if(sw_potential_mixing == ON) tag_mixing = "Potential-Mixing"
       
       if(number_of_solvers_applied == 1 .or. number_of_solvers_applied == 2) then
          if(dabs(etotal) < 1.e13) then
             ij = int(log10(dabs(etotal)))              ! T. Yamasaki, 2023/07/09
             if(dabs(etotal) < 1.e5) ij = 5             ! T. Yamasaki, 2023/07/09
             if(17-ij>9)  write(ndecimals,'(i2)') 17-ij ! T. Yamasaki, 2023/07/09
             if(17-ij<=9) write(ndecimals,'(i1)') 17-ij ! T. Yamasaki, 2023/07/09
          end if
       end if
#ifdef DEBUG_ITERATION_WRITE
       iteration_tmp = iteration
       iteration = iteration_tmp + 100000000
#endif       
       iterdigit = max(nint(log10(amax0(iteration,1)))+2,6)   ! T. Yamasaki, 2024/03/22
       if(iterdigit>=10) write(niterdigits,'(i2)') iterdigit  ! T. Yamasaki, 2024/03/22
       if(iterdigit<10)  write(niterdigits,'(i1)') iterdigit  ! T. Yamasaki, 2024/03/22
       
       if(number_of_solvers_applied == 1 .or. number_of_solvers_applied == 2) then
          if(number_of_solvers_applied == 1) then
             fmt2 = ",' EDEL = ',d14.6, ' : SOLVER = ',A, ' : ',A,' = ',A)"
          else  ! number_of_solvers_applied == 2
             fmt2 = ",' EDEL = ',d14.6, ' : SOLVER = ',A,' + ',A,' : ',A,' = ',A)"
          end if
          if(dabs(etotal) <1.e13) then
             fmt = "(' TOTAL ENERGY FOR',i"//trim(niterdigits)//",' -TH ITER=',f20."//trim(ndecimals)//trim(fmt2)  
          else
             fmt = "(' TOTAL ENERGY FOR',i"//trim(niterdigits)//",' -TH ITER=',d20.12"//trim(fmt2)
          end if
          if(number_of_solvers_applied == 1) then
             write(nfout,fmt) iteration,etotal,edel,trim(solver_names_applied(1)),trim(tag_mixing),trim(cdmixing_names_applied(1))
          !                                               T. Yamasaki, 2023/07/09, 2024/03/22
          else  ! number_of_solvers_applied == 2
             write(nfout,fmt) iteration, etotal, edel, trim(solver_names_applied(1)),trim(solver_names_applied(2)) &
                  & , trim(tag_mixing), trim(cdmixing_names_applied(1))
          !                                               T. Yamasaki, 2023/07/09, 2024/03/22
          end if
       else
          fmt = "(' TOTAL ENERGY FOR',I"//trim(niterdigits)//",' -TH ITER=',F0.12,2x,' edel = ',D14.6)"
          write(nfout,fmt) iteration,etotal,edel
          !                                               T. Yamasaki, 2024/03/22
       end if
       
       EPC_t = epc;  if ( flg_paw ) EPC_t = epc_paw                                          ! T.Y., 2023/07/09
       if(sw_output_xc_seperately==OFF)then
          write(nfout,610) ekinet,ehartr,exc,elocal,enonlc,eewald,-EPC_t,eentropy            ! T.Y., 2023/07/09
       else
          write(nfout,615) ekinet,ehartr,eex,ecor,exc,elocal,enonlc,eewald,-EPC_t,eentropy   ! T.Y., 2023/07/09
       endif
#ifdef DEBUG_ITERATION_WRITE
       iteration = iteration_tmp
#endif       
       if(sw_dipole_correction == ON) write(nfout,630) evdip,evext,edip
       if(sw_hubbard == ON) write(nfout,640) ehub0,ehub1
       if(sw_hybrid_functional == ON.and.sw_eval_vexx==ON) write(nfout,650) vexx,eexx,eexx-vexx
       if(sw_hybrid_functional == ON.and.sw_eval_vexx==OFF) write(nfout,655) eexx
       if(sw_external_potential==ON) write(nfout,670) eespot
       if(sw_fef == ON) write(nfout,680) eplr
       if(ntyp_vdw>0 .or. sw_vdw_correction == ON .and. vdw_method == VDW_DFTD3) write(nfout,690) evdw
#ifndef DISABLE_VDWDF
       if ( xctype == "vdwdf" .and. ( .not. oneshot ) ) write(nfout,690) ecnl_vdwdf
#endif

       if(way_of_smearing == COLD .or. way_of_smearing == PARABOLIC &
            & .and. abs(eentropy)>=DELTA ) write(nfout,620) etotal0

! =============== KT_add =========================================== 13.0E
       if (way_of_smearing == Fermi_Dirac ) write(nfout,620) etotal0
! ================================================================== 13.0E

       if (way_of_smearing == MP) write(nfout,620) etotal0

! === DEBUG by tkato 2011/10/07 ================================================
!      if(flg_paw) &
!      write(nfout,650) eohxc_paw,ehartr_paw,exc_ae,exc_ps,exc_ae-exc_ps
       if(flg_paw) then
          if(dabs(exc_ae) >1.d6 .or. dabs(exc_ps)>1.d6 )then
             write(nfout,661) eohxc_paw,ehartr_paw,exc_ae,exc_ps,exc_ae-exc_ps
          else
             write(nfout,660) eohxc_paw,ehartr_paw,exc_ae,exc_ps,exc_ae-exc_ps
          end if
       end if
! ==============================================================================

! ====================================== added by K. Tagami =============== 11.0
       if ( sw_magnetic_constraint == ON ) write(nfout,700) emag0, emag1
! ========================================================================= 11.0

! ====================================== added by K. Tagami =============== 11.0
!!#ifdef USE_ESPINORB
       if ( flg_paw ) then
          if ( SpinOrbit_Mode == ByPawPot ) then
             write(nfout,710) espinorb_old, espinorb_now
          endif
       endif
       if ( SpinOrbit_Mode == ByProjector .or. SpinOrbit_Mode == ZeffApprox &
            &                             .or. SpinOrbit_Mode == ReadFromPP ) then
          write(nfout,720) espinorb_now
       endif
!!!#endif
! ========================================================================= 11.0

       call flush(nfout)
    end if
#ifdef __TIMER_SUB__
    call timer_end(749)
#endif
!!$600 FORMAT(' ','TOTAL ENERGY FOR',I6,' -TH ITER=',F0.12,2x,' edel = ',D14.6)
610 FORMAT(' KI=',F20.12,' HA=',F20.12,' XC=',F20.12,' LO=',F20.12,/ &
    &      ,' NL=',F20.12,' EW=',F20.12,' PC=',d20.12,' EN=',F20.12)
615 FORMAT(' KI=',F20.12,' HA=',F20.12,' EX=',d20.12,' CR=',F20.12,' XC=',F20.12,' LO=',F20.12,/ &
    &      ,' NL=',F20.12,' EW=',F20.12,' PC=',d20.12,' EN=',F20.12)
620 FORMAT(' ','PHYSICALLY CORRECT ENERGY = ',d20.12)
630 FORMAT(' VD=',F15.7,' VE=',F15.7,' ED=',F15.7)
640 FORMAT(" HE=",F15.7," HP=",F15.7)
650 FORMAT(" VEXX=",F20.12," EEXX=",F20.12," EXX-VEXX=",F20.12)
655 FORMAT(" EXX=",F20.12)
660 FORMAT(' EOHXC_PAW=',F15.7,' HA_PAW=',F15.7,/ &
            ' XC_PAW_AE=',F15.7,' XC_PAW_PS=',F15.7,/ &
            '!XC_PAW_AE-XC_PAW_PS=',F15.7)
661 FORMAT(' EOHXC_PAW=',F15.7,' HA_PAW=',F15.7,/ &
            ' XC_PAW_AE=',d20.8,' XC_PAW_PS=',d20.8,/ &
            '!XC_PAW_AE-XC_PAW_PS=',F15.7)
670 FORMAT(' ES=',F20.12)
680 FORMAT(" EP=",F20.12)
690 FORMAT(" VDW=",F20.12)

! ================================== added by K. Tagami ================== 11.0
700 FORMAT(" Emag0=",F15.7," Emag1=",F15.7)
710 FORMAT(" ESpinOrb_old=",F15.7," ESpinOrb_now=",F15.7)
720 FORMAT(" ESpinOrb=",F15.7)
! ======================================================================== 11.0

  end subroutine sumup_all_energies

  function cal_edeltb() result (res)
    real(kind=DP) :: res
    res = etotal - etoold
    return
  end function cal_edeltb

  logical function m_TE_is_Divergent_core(nfout)
    integer, intent(in) :: nfout
    integer :: iincre = 0

    if(iteration_electronic == 1) iincre = 0
    !edeltb = etotal - etoold
    edeltb = cal_edeltb()
    if(ipri >= 2) write(nfout,'(" ! edeltb = ",d14.6 &
	& ," hr (= ",d14.6," hr/atom ) ( iter = ",i7," )")') edeltb,edeltb/natm2,iteration

    if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then
       m_TE_is_Divergent_core = .false.
       etoold = etotal
       return
    end if
    if(sw_fef == ON)then
       if(printable) write(nfout,'(a,d14.6)') ' !** delta pmac =',sum(pmac-pmac_old)/3.d0
    endif

!!$    if(edeltb > 1.d-1) then
!!$       iincre = iincre + 100
!!$    else if(edeltb > 1.d-2) then
!!$       iincre = iincre + 100
!!$    else if(edeltb >1.d-3) then
!!$       iincre = iincre + 70
!!$    else if(edeltb > 1.d-4) then
!!$       iincre = iincre + 20
!!$    else if(edeltb > 1.d-5) then
!!$       iincre = iincre + 10
!!$    else if(edeltb > 1.d-6) then
!!$       iincre = iincre + 5
!!$    else if(iincre > 0) then
!!$       iincre = iincre - 1
!!$    end if
    if(edeltb > 1.d-7 .and. iincre /= 0 .and. printable) &
         &    write(nfout,*) ' !W IINCRE is increasing as ', iincre
    if(iincre > IINCRE_CRITICAL) then
       if(printable) write(nfout,*)  ' !S IINCRE exceeds ',IINCRE_CRITICAL,'!'
       m_TE_is_Divergent_core = .true.
    else
       m_TE_is_Divergent_core = .false.
    end if
    etoold = etotal
  end function m_TE_is_Divergent_core

  subroutine m_TE_set_etotal_old
    etoold = 0.d0
  end subroutine m_TE_set_etotal_old

! --> T. Yamasaki 08 Aug. 2009
  real(kind=DP) function m_TE_edeltb()
    m_TE_edeltb = edeltb
  end function m_TE_edeltb
! <--
  logical function m_TE_is_converged(nfout)
    integer, intent(in) :: nfout
    integer   :: ncnv
    real(kind=DP) :: edelta

    call m_CtrlP_get_edelta(edelta)

    if(dabs(edeltb) < edelta*natm2 .and. m_ES_eekdif_cond()) then
       ncnv = m_CtrlP_ntcnvg_incre()
       if(printable) write(nfout,'(" edeltb = ",d12.4, " edelta = ",d12.4 &
            & , " ntcnvg = ",i7)') edeltb, edelta, ncnv
! --> T. Yamasaki, 25 July 2008
    else if(sub_delta_factor_is_given .and. dabs(edeltb) <edelta*sub_delta_factor*natm2) then
       ncnv = m_CtrlP_sub_ntcnvg_incre()
       if(printable) write(nfout,'(" edeltb = ",d12.4, " sub_edelta = ",d12.4 &
            & , " sub_ntcnvg = ",i7)') edeltb, edelta*sub_delta_factor, ncnv
    else
! <--
       call m_CtrlP_ntcnvg_reset()! k .Mae 030808
    end if
    m_TE_is_converged = m_CtrlP_ntcnvg_clear()
    if(m_TE_is_converged) call m_CtrlP_ntcnvg_reset()
  end function m_TE_is_converged

  function m_TE_tell_total_energy()
    real(kind=DP)::  m_TE_tell_total_energy
    m_TE_tell_total_energy = etotal
  end function m_TE_tell_total_energy

  function m_TE_tell_total_energy0()
    real(kind=DP)::  m_TE_tell_total_energy0
    m_TE_tell_total_energy0 = etotal0
  end function m_TE_tell_total_energy0

!!$! --> T. Yamasaki 17th Aug. 2009
!!$  function m_TE_tell_total_energy_mdfy()
!!$    real(kind=DP)::  m_TE_tell_total_energy_mdfy
!!$    if(sw_hubbard == ON) then
!!$       m_TE_tell_total_energy_mdfy = etotal - ehub0*0.5d0
!!$    else
!!$       m_TE_tell_total_energy_mdfy = etotal
!!$    end if
!!$  end function m_TE_tell_total_energy_mdfy
!!$! <--

  function m_TE_tell_band_energy(nfout,kv3)
    real(kind=DP)::  m_TE_tell_band_energy
    integer, intent(in) :: nfout,kv3
    integer :: ib, ik

! ============================== modified by K. Tagami ================== 11.0
!    call sum_eigenvalues(kv3)
!
    if ( noncol ) then
      call sum_eigenvalues_noncl(kv3)
    else
      call sum_eigenvalues(kv3)
    endif
! ======================================================================= 11.0

    m_TE_tell_band_energy = eband
    if(ipri >= 2) then
       if(map_k(ik) == myrank_k) then
          write(nfout,'(" -- sum_eigenvalues (m_TE_tell_band_energy) --")')
          write(nfout,'("   -- num_extra_bands = ",i8)') num_extra_bands

! ========================================= modified by K. Tagami ======= 11.0
!          do ik = 1, kv3
          do ik = 1, kv3, ndim_spinor
! ====================================================================== 11.0

             write(nfout,'("  ik = ",i5)') ik
             write(nfout,'(8f8.4)') (eko_l(ib,ik),ib=1,np_e) ! MPI
          end do
       end if
    endif
  end function m_TE_tell_band_energy

  function m_TE_tell_extended_band_energy(nfout,kv3)
    real(kind=DP)::  m_TE_tell_extended_band_energy
    integer, intent(in) :: nfout,kv3
    integer :: ib, ik

! ============================== modified by K. Tagami ================== 11.0
!    call sum_wholeeigenvalues(nfout,kv3)
!
    if ( noncol ) then
      call sum_wholeeigenvalues_noncl(nfout,kv3)
    else
      call sum_wholeeigenvalues(nfout,kv3)
    endif
! ======================================================================= 11.0

    m_TE_tell_extended_band_energy = eband_extendedrange
    if(ipri >= 2) then
       if(map_k(ik) == myrank_k) then
          write(nfout,'(" -- sum_eigenvalues (m_TE_tell_band_energy) --")')
          write(nfout,'("   -- num_conduction_bands_lmm = ",i8)') num_conduction_bands_lmm
! ========================================= modified by K. Tagami ======= 11.0
!          do ik = 1, kv3
          do ik = 1, kv3, ndim_spinor
! ======================================================================= 11.0

             write(nfout,'("  ik = ",i5)') ik
             write(nfout,'(8f8.4)') (eko_l(ib,ik),ib=1,np_e) ! MPI
          end do
       end if
    endif
  end function m_TE_tell_extended_band_energy

  subroutine sum_eigenvalues(kv3)
    integer, intent(in) :: kv3
    integer             :: ib, ik, is
    eband = 0.d0
    do is = ista_spin, iend_spin
    !do ik=1,kv3,af+1
    do ik = is, kv3-nspin+is, nspin
!    do ik =1, kv3, af+1
       if(map_k(ik) /= myrank_k) cycle ! MPI
       do ib = ista_e, iend_e, istep_e                 ! MPI
          if(nrvf_ordr(ib,ik) > neg - num_extra_bands) cycle
          eband = eband + eko_l(map_z(ib),ik)*occup_l(map_z(ib),ik)
       end do
    end do
    enddo
    call mpi_allreduce(MPI_IN_PLACE,eband,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
    eband = 2*eband/kv3 * (af+1)
  end subroutine sum_eigenvalues

! ================================- added by K. Tagami ====================== 11.0
  subroutine sum_eigenvalues_noncl(kv3)
    integer, intent(in) :: kv3
    integer             :: ib, ik
    real(kind=DP)       :: eband_mpi   ! MPI

    eband = 0.d0
    do ik =1, kv3, ndim_spinor
       if(map_k(ik) /= myrank_k) cycle ! MPI
       do ib = 1, np_e                 ! MPI
          if(nrvf_ordr(ib,ik) > neg - num_extra_bands) cycle
          eband = eband + eko_l(ib,ik)*occup_l(ib,ik)
       end do
    end do
!    call mpi_allreduce(eband,eband_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr) !MPI
    call mpi_allreduce(eband,eband_mpi,1,mpi_double_precision,mpi_sum,mpi_spin_group,ierr) !MPI
    eband = eband_mpi                  ! MPI
    eband = eband/ (kv3/ndim_spinor)
  end subroutine sum_eigenvalues_noncl
! ===================================================================== 11.0

  subroutine sum_wholeeigenvalues(nfout,kv3)
    integer, intent(in) :: nfout,kv3
    integer             :: ib, ik, is
    real(kind=DP)       :: eband_mpi   ! MPI
    integer             :: iband
    real(kind=DP)       :: weight, weight0

    if(num_conduction_bands_lmm < 0) then
       iband = neg
    else
       iband = ceiling(totch/2.0) + num_conduction_bands_lmm
    end if
    if(iband > neg - num_extra_bands) iband = neg - num_extra_bands
    eband_extendedrange = 0.d0
    weight = 0.d0; weight0 = 0.d0
    do is = ista_spin, iend_spin
    do ik = is, kv3-nspin+is, nspin
    !do ik =1, kv3, af+1
       weight = weight + qwgt(ik)
       if(map_k(ik) /= myrank_k) cycle ! MPI
       do ib = ista_e, iend_e, istep_e                 ! MPI
          if(ik==1 .and. ib == 1) weight0 = weight0 + occup_l(1,ik)
!!$          if(nrvf_ordr(ib,ik) > neg - num_extra_bands) cycle
          if(nrvf_ordr(ib,ik) > iband) cycle
          eband_extendedrange = eband_extendedrange + eko_l(map_z(ib),ik)*qwgt(ik)
       end do
    end do
    enddo
    if(npes > 1) then
       call mpi_allreduce(eband_extendedrange,eband_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr) !MPI
       eband_extendedrange = eband_mpi                  ! MPI
       call mpi_allreduce(weight0, eband_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       weight0 = eband_mpi
    end if
    eband_extendedrange = (3.0-nspin) * eband_extendedrange * (af+1)
    if(ipri >= 2) then
       weight = weight * (af+1)
       write(nfout,'("   -- num_conduction_bands_lmm = ",i8)') num_conduction_bands_lmm
       write(nfout,'(" weight (sum_wholeeigenvalues) = ", f8.4, " iband = ",i10)') weight,iband
       write(nfout,'(" weight0 (occup_l)             = ", f8.4)') weight0
       write(nfout,'(" eband_extendedrange = ", f16.8, "   eband = ",f16.8)') eband_extendedrange, eband
    end if
  end subroutine sum_wholeeigenvalues

! ================================ added by K. Tagami =============== 11.0
  subroutine sum_wholeeigenvalues_noncl(nfout,kv3)
    integer, intent(in) :: nfout,kv3
    integer             :: ib, ik
    real(kind=DP)       :: eband_mpi   ! MPI
    integer             :: iband
    real(kind=DP)       :: weight, weight0

    if(num_conduction_bands_lmm < 0) then
       iband = neg
    else
       iband = ceiling(totch) + num_conduction_bands_lmm
    end if
    if(iband > neg - num_extra_bands) iband = neg - num_extra_bands
    eband_extendedrange = 0.d0
    weight = 0.d0; weight0 = 0.d0

    do ik =1, kv3, ndim_spinor
       weight = weight + qwgt(ik)
       if(map_k(ik) /= myrank_k) cycle ! MPI
       do ib = 1, np_e                 ! MPI
          if(ik==1 .and. ib == 1) weight0 = weight0 + occup_l(1,ik)
!!$          if(nrvf_ordr(ib,ik) > neg - num_extra_bands) cycle
          if(nrvf_ordr(ib,ik) > iband) cycle
          eband_extendedrange = eband_extendedrange + eko_l(ib,ik)*qwgt(ik)
       end do
    end do
    if(npes > 1) then
!       call mpi_allreduce(eband_extendedrange,eband_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr) !MPI
       call mpi_allreduce(eband_extendedrange,eband_mpi,1,mpi_double_precision,mpi_sum,mpi_spin_group,ierr) !MPI
       eband_extendedrange = eband_mpi                  ! MPI
!       call mpi_allreduce(weight0, eband_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       call mpi_allreduce(weight0, eband_mpi,1,mpi_double_precision,mpi_sum,mpi_spin_group,ierr)
       weight0 = eband_mpi
    end if

    if(ipri >= 2) then
       weight = weight * (af+1)
       write(nfout,'("   -- num_conduction_bands_lmm = ",i8)') num_conduction_bands_lmm
       write(nfout,'(" weight (sum_wholeeigenvalues) = ", f8.4, " iband = ",i10)') weight,iband
       write(nfout,'(" weight0 (occup_l)             = ", f8.4)') weight0
       write(nfout,'(" eband_extendedrange = ", f16.8, "   eband = ",f16.8)') eband_extendedrange, eband
    end if
  end subroutine sum_wholeeigenvalues_noncl
! ============================================================== 11.0

  subroutine m_TE_wd_total_energy(nfcntn)
    integer, intent(in) :: nfcntn
    !!$ print *,' tag_total_energy'
    if(mype==0) then
       write(nfcntn,*) tag_total_energy
       write(nfcntn,'(2d24.16)') etotal, etoold
       !edeltb = etotal - etoold
       !edeltb = cal_edeltb()
    endif
  end subroutine m_TE_wd_total_energy

  subroutine m_TE_rd_total_energy(nfcntn)
    integer, intent(in) :: nfcntn
    logical             :: EOF_reach, tag_is_found

    if(mype==0) then
       call rewind_to_tag0(nfcntn,len_tag_total_energy,tag_total_energy &
            &, EOF_reach, tag_is_found,str,len_str)
       if(.not.tag_is_found) then
          call phase_error_with_msg(6,' tag_total_energy is not found',__LINE__,__FILE__)
       else
          read(nfcntn,*) etotal, etoold
          etoold = etotal
       end if
    endif
    call mpi_bcast(etotal,1 &
         & ,mpi_double_precision,0,MPI_CommGroup,ierr)
    call mpi_bcast(etoold,1 &
         & ,mpi_double_precision,0,MPI_CommGroup,ierr)
  end subroutine m_TE_rd_total_energy

  logical function m_TE_Converged_Hubbard_Energy(nfout)
    integer, intent(in) :: nfout

    m_TE_Converged_Hubbard_Energy = .false.
    !!$if(ehub0 <= critical_ehub) then
    !!$   m_TE_Converged_Hubbard_Energy = .true.
    !!$   write(nfout,'(" Interation will be stoped &
    !!$   &because Hubbard energy is lower than ",f10.5,".")') critical_ehub
    !!$end if

    if(abs(ehub0-ehub0_old) <= delta_ehub) then
       m_TE_Converged_Hubbard_Energy = .true.
       write(nfout,'(" Interation will be stoped &
       &because Hubbard energy is converged within ",f10.5,".")') delta_ehub
    else
       ehub0_old = ehub0
    end if

  end function m_TE_Converged_Hubbard_Energy


end module m_Total_Energy
