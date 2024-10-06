!#define _DEBUG_WRITE_DFTU_MPI_PROCESSES_
!!$#define DEBUG_WRITE
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
! === DEBUG by tkato 2012/11/05 ================================================
!!$  use m_PlaneWaveBasisSet,    only : gr_l,m_pwBS_kinetic_energies_3D,iba,kg1
  use m_Parallelization,      only : np_g1k, ista_g1k, iend_g1k, ista_spin, iend_spin, &
                                     nrank_s,mpi_keg_world
! ==============================================================================
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
       &                           , sw_dipole_correction &
       &                           , sw_hubbard, proj_attribute &
       &                           , critical_ehub, delta_ehub &
       &                           , num_conduction_bands_lmm &
       &                           , sw_hybrid_functional, sw_eval_vexx, sw_retard_eigval_evaluation, sw_fef &
       &                           , in_line_minimization &
       &                           , sw_external_potential, ekmode, icond, sw_rsb, sw_output_xc_seperately &
       &                           , number_of_cdmixing_applied, cdmixing_names_applied &
       &                           , number_of_solvers_applied, solver_names_applied, len_solvername &
       &                           , sw_potential_mixing, sw_communicator_for_chg &
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
       &                           , ista_atm, iend_atm, myrank_g, nrank_g, nrank_e, nrank_chg &
       &                           , ista_e, iend_e, istep_e
  use m_Const_Parameters,     only : Valence_plus_PC_Charge,PAI4,UP,DOWN&
       &                           , EXC_ONLY, VXC_AND_EXC, IINCRE_CRITICAL,DP,CMPLDP &
       &                           , len_tag_total_energy, tag_total_energy &
       &                           , PARABOLIC, COLD, ON, OFF &
       &                           , FIXED_CHARGE, FIXED_CHARGE_CONTINUATION, GAMMA, DELTA &
       &                           , MDKOSUGI, MDDAVIDSON, INITIAL, CONTINUATION, COORDINATE_CONTINUATION &
       &                           , FIXED_CHARGE, FIXED_CHARGE_CONTINUATION, MATRIXDIAGON &
       &                           , DELTA_MOVING_AVERAGE, SLOPE, DELTA_V &
       &                           , VDW_DFTD3 &
       &                           , MP, Hartree
  use m_Dipole,               only : edip_ion, eext_ion, vdip_l, vext_l
  use m_Hubbard,              only : m_Hubbard_energy
  use m_FiniteElectricField,  only : m_FEF_polarization, pmac, pmac_old
  use m_ES_ExactExchange,     only : m_ES_EXX_gather_valence_states,m_ES_EXX_energy &
       &                           , m_ES_EXX_energy2
#if 0
       &  , m_ES_EXX_gather_valence_states_k
#endif
  use m_PAW_XC_Potential,     only : m_PAW_XC_cal_potential,exc_ae,exc_ps &
       &                           , m_PAW_XC_cal_potential_sphex2,m_PAW_XC_allocation, m_PAW_XC_deallocation
!!$                                    , m_PAW_XC_cal_potential_sym
  use m_PAW_Hartree,          only : m_PAWH_get_dion_hartree_now
  use m_PAW_ChargeDensity,    only : calcGaussLegendreIntegration &
       &                            , calcSphericalHarmonicsExpansion

  use m_PlaneWaveBasisSet,    only : kgp, kg

! ====================================== added by K. Tagami ================ 11.0
  use m_Control_Parameters,    only : noncol, ndim_spinor, ndim_magmom, ndim_chgpot
  use m_PseudoPotential,       only : dion0_noncl, nlmt
  use m_Charge_Density,        only : hsi
  use m_ES_NonCollinear,       only : m_ES_MagMom_To_DensMat_hsr, &
       &                              m_ES_MagMom_to_DensMat_Dhub
!
  use m_Crystal_Structure,     only : sw_magnetic_constraint
  use m_ES_Mag_Constraint,     only : m_ES_calc_MagConstraint_Energy
  use m_Hubbard,              only : m_Hubbard_energy_noncl, &
       &                             m_Hubbard_energy2_noncl, &
       &                             m_Hubbard_energy3_noncl
! ========================================================================== 11.0

! ======================= KT_add ================== 13.0E
  use m_Const_Parameters,   only : Fermi_Dirac
! ================================================= 13.0E

! =========== KT_add ========== 13.0U2
  use m_Control_Parameters,  only : sw_potential_mixing, use_metagga
! ============================= 13.0U2

#ifndef DISABLE_VDWDF
  use m_vdWDF,  only : ecnl_vdwdf
#endif

  use m_Control_Parameters,  only : m_CtrlP_get_isolver_now

! === Positron SCF ==== 2015/11/28
  use m_Control_Parameters,  only : sw_positron, positron_method
  use m_Const_Parameters,   only :  positron_GGGC
  use m_epc_potential,  only : ecorr_pztr => epc
! ==================== 2015/11/28
  use m_ES_LHXC, only : m_ESlhxc_delta_vmax
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

! === positron
  real(kind=DP) :: ekin_pztr, elocal_pztr, ehartr_ep, eohxc_pztr
! ===

  integer,private, parameter     :: len_str = 132
  character(len=len_str),private    ::  str

  real(kind=DP), allocatable, dimension(:) :: ehist
  integer :: ihist
  real(kind=DP) :: emova, emovaold

#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
  integer, parameter :: DEBUGPRINTLEVEL = 1
#else
  integer, parameter :: DEBUGPRINTLEVEL = 2
#endif

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

! === DEBUG by tkato 2012/11/05 ================================================
  subroutine get_kinetic_energy_directly(nfout)
     integer,intent(in) :: nfout
     integer :: ik,ib,ig,iadd,ispin
     real(kind=DP), allocatable, dimension(:)       :: ekin
     real(kind=DP) :: ekinet_tmp,ekinet_mpi
     allocate(ekin(1:maxval(np_g1k(:))));ekin=0.d0
     ekinet = 0.d0
     do ispin = ista_spin, iend_spin, (af+1)
     !do ik=1,kv3,af+1
     do ik = ispin, kv3-nspin+ispin, nspin
        call m_pwBS_kinetic_energies(ik,vkxyz,ekin)
        if(map_k(ik)/=myrank_k)cycle
        do ib=1,np_e
           ekinet_tmp=0.d0
           if (kimg==1) then
              do ig=ista_g1k(ik), iend_g1k(ik)
                 iadd = ig - ista_g1k(ik) + 1
                 ekinet_tmp = ekinet_tmp+2.d0*ekin(iadd)*zaj_l(iadd,ib,ik,1)**2
              enddo
           else if (kimg==2) then
              do ig=ista_g1k(ik), iend_g1k(ik)
                 iadd = ig - ista_g1k(ik) + 1
                 ekinet_tmp = ekinet_tmp+2.d0*ekin(iadd)*(zaj_l(iadd,ib,ik,1)**2 &
              & +zaj_l(iadd,ib,ik,2)**2)
              enddo
           endif
           if(k_symmetry(ik)==GAMMA) ekinet_tmp = ekinet_tmp*2.d0
           ekinet = ekinet+occup_l(ib,ik)/dble(kv3)*ekinet_tmp
        enddo
     enddo
     enddo
     if(af==1) ekinet = ekinet*2.d0
     ekinet_mpi = 0.d0
     call mpi_allreduce(ekinet,ekinet_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
     ekinet = ekinet_mpi
     deallocate(ekin)
  end subroutine get_kinetic_energy_directly
! ==============================================================================



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


! ============================ added by K. Tagami ========================= 11.0
  subroutine get_magnetic_constraint_energy
    if ( sw_magnetic_constraint == ON ) then
       call m_ES_calc_MagConstraint_Energy( emag1, emag0 )
    endif
  end subroutine get_magnetic_constraint_energy
! ========================================================================= 11.0


  subroutine sumup_all_energies(nfout,exc,display_on)
    integer, intent(in) :: nfout
    logical, intent(in) :: display_on
    real(kind=DP), intent(in) :: exc
    integer :: i

    real(kind=DP) :: edel, EPC_t
    character(len("potential_mixing")) :: tag_mixing
    character(len=2) ::ndecimals    ! T. Yamasaki, 2023/07/09
    integer :: ij                   ! T. Yamasaki, 2023/07/09
    character(256) :: fmt,fmt2      ! T. Yamasaki, 2023/07/09
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

    if(num_regions>0)then
       do i=1,num_regions
          if(regions(i)%tally) etotal0 = etotal0+regions(i)%energy
       enddo
    endif

    etotal = etotal0+eentropy
    call mpi_bcast(etotal,1 &
         & ,mpi_double_precision,0,MPI_CommGroup,ierr)
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
          else
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
          else
             write(nfout,fmt) iteration, etotal, edel, trim(solver_names_applied(1)),trim(solver_names_applied(2)) &
                  & , trim(tag_mixing), trim(cdmixing_names_applied(1))
          !                                               T. Yamasaki, 2023/07/09, 2024/03/22
          end if
       else
          fmt = "(' TOTAL ENERGY FOR',I"//trim(niterdigits)//",' -TH ITER=',F0.12,2x,' edel = ',D14.6)"
          write(nfout,fmt) iteration,etotal,edel
          !                                               T. Yamasaki, 2024/03/22
       end if
       
#ifdef DEBUG_ITERATION_WRITE
       iteration = iteration_tmp
#endif

       EPC_t = epc;  if ( flg_paw ) EPC_t = epc_paw                                          ! T.Y., 2023/07/09
       if(sw_output_xc_seperately==OFF)then
          write(nfout,610) ekinet,ehartr,exc,elocal,enonlc,eewald,-EPC_t,eentropy            ! T.Y., 2023/07/09
       else
          write(nfout,615) ekinet,ehartr,eex,ecor,exc,elocal,enonlc,eewald,-EPC_t,eentropy   ! T.Y., 2023/07/09
       endif

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

       if (way_of_smearing == Fermi_Dirac ) write(nfout,620) etotal0  ! === KT_add === 13.0E

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
! ======================================================================== 11.0

  end subroutine sumup_all_energies

  function cal_edeltb() result (res)
    real(kind=DP) :: res
    res = 1e+30
    if(convergence_criteria == DELTA_MOVING_AVERAGE) then
       res = emova-emovaold
    else if (convergence_criteria == DELTA_V )then
       res = m_ESlhxc_delta_vmax() * natm2
    else if (convergence_criteria == SLOPE )then
       res = m_TE_get_slope_of_ehist()
    else
       res = etotal - etoold
    endif
    return
  end function cal_edeltb

  logical function m_TE_is_Divergent_core(nfout)
    integer, intent(in) :: nfout
    integer :: iincre = 0

    if(iteration_electronic == 1) iincre = 0
    !edeltb = etotal - etoold
    edeltb = cal_edeltb()
    if(ipri >= DEBUGPRINTLEVEL) write(nfout,'(" ! edeltb = ",d14.6 &
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
!!$    if(iteration_electronic >= 10) then
!!$       if(edeltb/natm2 >= 1.0/Hartree) then
!!$          iincre = iincre + (edeltb/natm2)/(1.0/Hartree)
!!$       else
!!$          iincre = iincre - 1
!!$       end if
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
       if(printable .and. ipri>0) write(nfout,'(" edeltb = ",d12.4, " edelta = ",d12.4 &
            & , " ntcnvg = ",i7)') edeltb, edelta, ncnv
! --> T. Yamasaki, 25 July 2008
    else if(sub_delta_factor_is_given .and. dabs(edeltb) <edelta*sub_delta_factor*natm2) then
       ncnv = m_CtrlP_sub_ntcnvg_incre()
       if(printable .and. ipri>0) write(nfout,'(" edeltb = ",d12.4, " sub_edelta = ",d12.4 &
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

!===============================================================================
!!$  subroutine m_TE_total_energy_3D(nfout,display_on,kv3)
  subroutine m_TE_total_energy(nfout,display_on,kv3)
    use m_XC_Potential,         only : m_XC_cal_potential_3D
    use m_XC_Potential_2D,         only : m_XC_cal_potential
#if 0
    use m_Control_Parameters,   only : sw_distribute_wf,force_exx_energy1
#else
    use m_Control_Parameters,   only : force_exx_energy1
#endif
    integer, intent(in) :: nfout
    logical, intent(in) :: display_on
    integer, intent(in) :: kv3
    real(kind=DP) :: edeltb_now
    real(kind=DP),allocatable,dimension(:) :: ehistmp
    integer :: isolver,sw_submat
    integer :: ih,it
    integer             :: id_sname = -1   , i
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
#ifdef __TIMER_SUB__
    call timer_sta(739)
#endif
    call tstatc0_begin('Total_Energy(including xc_pot) ',id_sname)
    call get_band_energy_3D(kv3)
    call get_xc_and_HE_of_old_CD_3D(vxc_l)
    call get_local_potential_energy_3D
    call get_hubbard_energy_3D(nfout)
    call get_nonlocal_potential_energy_3D
! === DEBUG by tkato 2012/11/07 ================================================
#if 1
    if(in_line_minimization) then
       call m_XC_cal_potential_3D(nfout,Valence_plus_PC_Charge, chgqt, EXC_ONLY)
    else
       call m_XC_cal_potential_3D(nfout,Valence_plus_PC_Charge, chgqt, VXC_AND_EXC)
    endif
#else
    if(in_line_minimization) then
       call m_XC_cal_potential(nfout,Valence_plus_PC_Charge, chgqt, EXC_ONLY)
    else
       call m_XC_cal_potential(nfout,Valence_plus_PC_Charge, chgqt, VXC_AND_EXC)
    endif
#endif
! ==============================================================================
#ifdef ENABLE_ESM_PACK
    if(sw_esm==OFF)then
       call get_hartree_energy_3D
    endif
#else
    call get_hartree_energy_3D
#endif
    if(sw_dipole_correction ==  ON) then
       call get_dipole_energy_3D(vdip_l,vext_l)
    end if
    if(sw_fef == ON) then
       call m_FEF_polarization(nfout,eplr)
    endif
! === Support Hybrid on 3D_Parallel by tkato 2013/02/10 ========================
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
       if(sw_retard_eigval_evaluation==OFF.or.isolver==MDDAVIDSON.or.isolver==MDKOSUGI.or.force_exx_energy1) then
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
       call m_PAW_XC_allocation(nfout)
#if 0
        call get_xc_and_HE_of_old_CD_paw_3D
#endif
! === DEBUG by tkato 2011/10/04 ================================================
!       call m_PAW_XC_cal_potential(nfout,EXC_ONLY,flg_symmtry)
        if(calcSphericalHarmonicsExpansion) then
           call m_PAW_XC_cal_potential_sphex2(nfout,EXC_ONLY)
        end if
        if(calcGaussLegendreIntegration)then
           call m_PAW_XC_cal_potential(nfout,EXC_ONLY,flg_symmtry)
        endif
! ==============================================================================
        call m_PAWH_get_dion_hartree_now(nfout)
        call get_hartree_energy_paw_3D
#if 1
        call get_kinetic_local_energy_paw(nfout)
#endif
       call m_PAW_XC_deallocation(nfout)
    end if
! ==============================================================================
! === DEBUG by tkato 2012/11/05 ================================================
!   call get_kinetic_energy(nfout)
#ifdef ENABLE_ESM_PACK
    if ( sw_esm==ON .or. sw_hybrid_functional==ON .or. sw_communicator_for_chg==ON &
         &          .or. xctype=='vdwdf' .or. use_metagga ) then
       call get_kinetic_energy_directly(nfout)
    else
       if(iteration_electronic == 1) then
          call get_kinetic_energy_directly(nfout)
       else
          call get_kinetic_energy(nfout)
       end if
    endif
#else
    if ( sw_hybrid_functional==ON .or. iteration_electronic==1 &
         &       .or. xctype=='vdwdf' .or. use_metagga ) then
       call get_kinetic_energy_directly(nfout)
    else
       call get_kinetic_energy(nfout)
    endif
#endif
! ==============================================================================
    call sumup_all_energies(nfout,exc,display_on)
    if(convergence_criteria == DELTA_MOVING_AVERAGE .or. convergence_criteria == SLOPE)then
      if(iteration_electronic == 1)then
        if(.not.allocated(ehist)) then
          allocate(ehist(nsamp))
        endif
        ehist = 0.d0
        emova = 0.d0
        ihist = 0
      endif
      if(iteration_electronic<=nsamp) then
        ihist = ihist + 1
        ehist(ihist) = etotal
      else
        allocate(ehistmp(nsamp));ehistmp=ehist
        do ih = 2,nsamp
           ehist(ih-1) = ehistmp(ih)
        enddo
        ehist(nsamp) = etotal
        deallocate(ehistmp)
      endif
      call cal_moving_average()
      if(printable .and. ipri>=2)then
        write(nfout,'(a,3f20.10,i10)') '!** moving average (curr step, last step, delta/natm, ihist)',&
           & emova,emovaold,(emova-emovaold)/real(natm2),ihist
        write(nfout,'(a,10f20.10)') '!** ehist ',(ehist(it),it=1,10)
        write(nfout,'(a,f20.10)')   '!** slope ',m_TE_get_slope_of_ehist()
      endif
    endif
    call tstatc0_end(id_sname)
#ifdef __TIMER_SUB__
    call timer_end(739)
#endif
  end subroutine m_TE_total_energy

  subroutine cal_moving_average()
    real(kind=DP) :: sume
    integer :: ii
    sume = 0.d0
    do ii=1,ihist
       sume = sume + ehist(ii)
    enddo
    emovaold = emova
    emova = sume/real(ihist)
  end subroutine cal_moving_average

  function m_TE_get_slope_of_ehist() result(res)
    integer :: ii
    real(kind=DP) :: res
    real(kind=DP) :: sx,sy,sxx,sxy,syy
    sx=0.d0;sy=0.d0;sxx=0.d0;sxy=0.d0;syy=0.d0
    if(ihist == 1) then
       res = ehist(1)
       return
    endif
    do ii=1,ihist
       sx = sx + real(ii)
       sy = sy + ehist(ii)
       sxx = sxx + real(ii)*real(ii)
       sxy = sxy + ehist(ii)*real(ii)
       syy = ehist(ii)*ehist(ii)
    enddo
    res = (ihist*sxy-sx*sy)/(ihist*sxx-sx*sx)
  end function m_TE_get_slope_of_ehist

  function m_TE_get_moving_average() result(res)
    real(kind=DP) :: res
    res = emova
    return
  end function m_TE_get_moving_average

  function m_TE_get_moving_average_old() result(res)
    real(kind=DP) :: res
    res = emovaold
    return
  end function m_TE_get_moving_average_old

!===============================================================================
  subroutine get_band_energy_3D(kv3)
    use m_Parallelization,      only : mpi_kg_world
! === DEBUG by tkato 2012/04/04 ================================================
    use m_Parallelization,      only : mpi_ge_world, nrank_k
! ==============================================================================

    integer, intent(in) :: kv3
    integer             :: ib, ik, is
    real(kind=DP)       :: eband_mpi   ! MPI
#ifdef __TIMER_SUB__
    call timer_sta(740)
#endif
    eband = 0.d0
#ifdef __TIMER_DO__
  call timer_sta(863)
#endif
    do is = ista_spin, iend_spin
    do ik = is, kv3-nspin+is, nspin
       if(map_k(ik) /= myrank_k) cycle ! MPI
       do ib = 1, np_e                 ! MPI
          eband = eband + occup_l(ib,ik)*eko_l(ib,ik)
       end do
    end do
    end do
#ifdef __TIMER_DO__
  call timer_end(863)
#endif

    if(nrank_e > 1) then
#ifdef __TIMER_COMM__
  call timer_barrier(mpi_kg_world)
  call timer_sta(864)
#endif
       call mpi_allreduce(eband,eband_mpi,1,mpi_double_precision,mpi_sum,mpi_kg_world,ierr) !MPI
#ifdef __TIMER_COMM__
  call timer_end(864)
#endif
       eband = eband_mpi
    end if
! === DEBUG by tkato 2012/04/04 ================================================
    if(nrank_k > 1) then
       call mpi_allreduce(eband,eband_mpi,1,mpi_double_precision,mpi_sum,mpi_ge_world,ierr) !MPI
       eband = eband_mpi
    end if
! ==============================================================================
    if(nrank_s > 1 .and. nrank_k==1) then
       call mpi_allreduce(eband,eband_mpi,1,mpi_double_precision,mpi_sum,mpi_keg_world,ierr) !MPI
       eband = eband_mpi
    endif
    eband = 2*eband/kv3 * (af+1)
#ifdef __TIMER_SUB__
    call timer_end(740)
#endif
  end subroutine get_band_energy_3D

!===============================================================================
  subroutine get_xc_and_HE_of_old_CD_3D(vxc_l)
    use m_Parallelization,     only : mpi_ke_world, mpi_kg_world, mpi_chg_world

    real(kind=DP), intent(in) :: vxc_l(ista_kngp:iend_kngp,kimg,nspin)
    integer ik, i, ispin
    integer ist !mpi
    real(kind=DP) :: eohxc_mpi

#ifdef __TIMER_SUB__
    call timer_sta(741)
#endif
    eohxc = 0.d0
    do ik = 1, kimg
#ifdef __TIMER_DO__
  call timer_sta(865)
#endif
       do ispin = 1, nspin
          if(sw_communicator_for_chg == OFF)then
             if(myrank_g==0) then
                eohxc = eohxc + vxc_l(1,ik,ispin)*chgqt(1,ik,ispin)
             endif
          else
             if(mype==0) then
                eohxc = eohxc + vxc_l(1,ik,ispin)*chgqt(1,ik,ispin)
             endif
          endif
       end do
#ifdef __TIMER_DO__
  call timer_end(865)
#endif
       ist = ista_kngp
       if(ist == 1) ist = 2

#ifdef __TIMER_DO__
  call timer_sta(866)
#endif
       if(nspin == 1) then
          do i = ist, iend_kngp  !for mpi
             eohxc = eohxc &
                  & + (vxc_l(i,ik,1)+PAI4*chgqto(i,ik,1)/gr_l(i)**2) &
                  &   * chgqt(i,ik,1)
          end do
       else if(nspin == 2) then
          do i = ist, iend_kngp  !for mpi
             eohxc = eohxc &
                  & +  (vxc_l(i,ik,UP)*chgqt(i,ik,UP)&
                  &   + vxc_l(i,ik,DOWN)*chgqt(i,ik,DOWN)) &
                  & + PAI4*(chgqto(i,ik,UP)+chgqto(i,ik,DOWN))&
                  & /gr_l(i)**2 &
                  &  *(chgqt(i,ik,UP)+chgqt(i,ik,DOWN))
          end do
       end if
#ifdef __TIMER_DO__
  call timer_end(866)
#endif
    end do
    if(nrank_chg > 1) then
#ifdef __TIMER_COMM__
  call timer_barrier(mpi_chg_world)
  call timer_sta(867)
#endif
       call mpi_allreduce(eohxc,eohxc_mpi,1,mpi_double_precision,mpi_sum,mpi_chg_world,ierr)
#ifdef __TIMER_COMM__
  call timer_end(867)
#endif
       eohxc = eohxc_mpi
    end if

    eohxc = univol*eohxc
#ifdef __TIMER_SUB__
    call timer_end(741)
#endif
  end subroutine get_xc_and_HE_of_old_CD_3D

!===============================================================================
  subroutine get_local_potential_energy_3D
    use m_Parallelization,     only : mpi_ke_world, mpi_chg_world

    integer       :: ik, it, i, ig
    real(kind=DP) :: eloca1_mpi
    real(kind=DP) :: eloclr,eloclr_mpi

#ifdef __TIMER_SUB__
    call timer_sta(742)
#endif
    eloca1 = 0.d0
#ifdef __TIMER_DO__
  call timer_sta(868)
#endif
    do ik = 1, kimg
       if(nspin == 1) then
!OCL NOFLTLD
          do it = 1,ntyp
             do i = ista_kngp, iend_kngp !for mpi
                eloca1 = eloca1 &
                     & + psc_l(i,it)*zfm3_l(i,it,ik)*chgqt(i,ik,1)
             end do
          end do
       else
!OCL NOFLTLD
          do it = 1, ntyp
             do i = ista_kngp, iend_kngp  !for mpi
                eloca1 = eloca1 &
                     &   + psc_l(i,it)*zfm3_l(i,it,ik)&
                     &             *(chgqt(i,ik,UP)+chgqt(i,ik,DOWN))
             end do
          end do
       endif
    end do
#ifdef __TIMER_DO__
  call timer_end(868)
#endif
    if(nrank_chg > 1) then
#ifdef __TIMER_COMM__
  call timer_barrier(mpi_chg_world)
  call timer_sta(869)
#endif
       call mpi_allreduce(eloca1,eloca1_mpi,1,mpi_double_precision,mpi_sum,mpi_chg_world,ierr)
#ifdef __TIMER_COMM__
  call timer_end(869)
#endif
       eloca1 = eloca1_mpi
    end if

    eloca1    = univol*eloca1
#ifdef ENABLE_ESM_PACK
    if(sw_esm==OFF) then
       elocal    = eloca1 + etot1*totch
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
      & mpi_double_precision,mpi_sum,mpi_chg_world,ierr)
      eloclr = eloclr_mpi
      eloclr = univol*eloclr
      elocal = elocal + eloclr
    endif
#endif

#ifdef __TIMER_SUB__
    call timer_end(742)
#endif
  end subroutine get_local_potential_energy_3D

!===============================================================================
  subroutine get_hubbard_energy_3D(nfout)
    integer, intent(in) :: nfout
#ifdef __TIMER_SUB__
    call timer_sta(743)
#endif
    ehub0 = 0.d0
    ehub1 = 0.d0
    if(sw_hubbard==ON) then
      call m_Hubbard_energy(ehub0)
      call get_hubbard_potential_energy_3D(ehub1)
    end if
    if(ipri >= 2) then
       write(nfout,'(" !D EHUB0 = ",f20.14)') ehub0
       write(nfout,'(" !D EHUB1 = ",f20.14)') ehub1
    endif
#ifdef __TIMER_SUB__
    call timer_end(743)
#endif
  end subroutine get_hubbard_energy_3D

!===============================================================================
  subroutine get_hubbard_potential_energy_3D(ehub1)
    use m_Parallelization,     only : mpi_ke_world, mpi_chg_world
    real(kind=DP), intent(out) :: ehub1
    integer :: ispin, it,lmt1, lmt2, il1, im1, il2, im2, ia
    integer :: ih, l1p, ierr
    real(kind=DP) :: fac, mpi
    ehub1 = 0.d0
    do ispin = 1, nspin, af+1
       do ia = ista_atm, iend_atm
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
    if(nrank_chg > 1) then
       call mpi_allreduce(ehub1,mpi,1,mpi_double_precision,mpi_sum,mpi_chg_world,ierr)
       ehub1 = mpi
    end if
    ehub1 = ehub1*(af+1)
  end subroutine get_hubbard_potential_energy_3D

!===============================================================================
  subroutine get_nonlocal_potential_energy_3D
    use m_Parallelization,     only : mpi_ke_world, mpi_chg_world

    integer :: ispin, it,lmt1, lmt2, il1, im1, il2, im2, ia, ierr
    real(kind=DP) :: fac, mpi
#ifdef __TIMER_SUB__
    call timer_sta(744)
#endif
    enonlc = 0.d0
#ifdef __TIMER_DO__
  call timer_sta(870)
#endif
    do ispin = 1, nspin, af+1
       do it = 1, ntyp
          if(ipaw(it)==0) then
              do lmt1 = 1, ilmt(it)
                 il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
                 do lmt2 = lmt1, ilmt(it)
                    il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
                    if( il1 /= il2 .or. im1 /= im2) cycle
    !xocl spread do/ind_katm
                    do ia = ista_atm, iend_atm
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
                    do ia = ista_atm, iend_atm
                       if(ityp(ia) /= it) cycle
                       fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0
                       enonlc = enonlc &
                            & + fac*dion_paw(lmt1,lmt2,ispin,ia)*hsrt(ia,lmt1,lmt2,ispin)
                    end do
    !xocl end spread
                 end do
              end do
          end if
       end do
    end do
#ifdef __TIMER_DO__
  call timer_end(870)
#endif
    if(nrank_g > 1) then
#ifdef __TIMER_COMM__
  call timer_barrier(mpi_ke_world)
  call timer_sta(871)
#endif
       call mpi_allreduce(enonlc,mpi,1,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
#ifdef __TIMER_COMM__
  call timer_end(871)
#endif
       enonlc = mpi
    end if
    enonlc = enonlc*(af+1)
    !!$if(sw_hubbard==ON) enonlc = enonlc - ehub1
#ifdef __TIMER_SUB__
    call timer_end(744)
#endif
  end subroutine get_nonlocal_potential_energy_3D
!===============================================================================

  subroutine get_hartree_energy_3D
    use m_Parallelization,     only : mpi_ke_world, mpi_chg_world

    integer ik, i
    integer :: ist !mpi
    real(kind=DP) :: ehartr_mpi

#ifdef __TIMER_SUB__
    call timer_sta(745)
#endif
    ehartr_mpi = 0.d0
    ehartr = 0.d0
#ifdef __TIMER_DO__
  call timer_sta(872)
#endif
    do ik = 1, kimg
       ist = ista_kngp
       if(ist == 1) ist = 2

       if(nspin == 1) then
          do i = ist, iend_kngp !for mpi
             ehartr_mpi  = ehartr_mpi + (chgqt(i,ik,1)/gr_l(i))**2
          end do
       else if(nspin == 2) then
          do i = ist, iend_kngp  !for mpi
             ehartr_mpi  = ehartr_mpi &
                  &  +((chgqt(i,ik,UP)+chgqt(i,ik,DOWN))/gr_l(i))**2
          end do
       end if
    end do
#ifdef __TIMER_DO__
  call timer_end(872)
#endif
#ifdef __TIMER_COMM__
  call timer_barrier(mpi_chg_world)
  call timer_sta(873)
#endif
    call mpi_allreduce(ehartr_mpi,ehartr,1,mpi_double_precision,mpi_sum,mpi_chg_world,ierr)
#ifdef __TIMER_COMM__
  call timer_end(873)
#endif
    ehartr    = univol*PAI4*ehartr*0.5d0
#ifdef __TIMER_SUB__
    call timer_end(745)
#endif
  end subroutine get_hartree_energy_3D
!===============================================================================

  subroutine get_dipole_energy_3D(vdip_l,vext_l)
    use m_Parallelization,     only : mpi_ke_world, mpi_chg_world

    real(kind=DP), intent(in) :: vdip_l(ista_kngp:iend_kngp,kimg)
    real(kind=DP), intent(in) :: vext_l(ista_kngp:iend_kngp,kimg)
    integer ik, i, ispin, ierr
    integer ist !mpi
    real(kind=DP) :: evdip_mpi
    real(kind=DP) :: evext_mpi
#ifdef __TIMER_SUB__
    call timer_sta(746)
#endif

    evdip = 0.d0
    evext = 0.d0
#ifdef __TIMER_DO__
  call timer_sta(874)
#endif
    do ik = 1, kimg
       do ispin = 1, nspin
          if(sw_communicator_for_chg == OFF)then
             if(myrank_g==0) then
                evdip = evdip + vdip_l(1,ik)*chgqt(1,ik,ispin)
                evext = evext + vext_l(1,ik)*chgqt(1,ik,ispin)
             endif
          else
             if(mype==0) then
                evdip = evdip + vdip_l(1,ik)*chgqt(1,ik,ispin)
                evext = evext + vext_l(1,ik)*chgqt(1,ik,ispin)
             endif
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
#ifdef __TIMER_DO__
  call timer_end(874)
#endif
    if(nrank_chg > 1) then
#ifdef __TIMER_COMM__
  call timer_barrier(mpi_chg_world)
  call timer_sta(875)
#endif
       call mpi_allreduce(evdip,evdip_mpi,1,mpi_double_precision,mpi_sum,mpi_chg_world,ierr)
#ifdef __TIMER_COMM__
  call timer_end(875)
#endif
       evdip =evdip_mpi
    end if
    evdip = evdip * univol
    evext = evext * univol
    edip = edip_ion + eext_ion + evdip*0.5d0 + evext
#ifdef __TIMER_SUB__
    call timer_end(746)
#endif
  end subroutine get_dipole_energy_3D
!===============================================================================

  subroutine get_xc_and_HE_of_old_CD_paw_3D
    use m_Parallelization,     only :mpi_ke_world

    integer :: ispin, it,lmt1, lmt2, il1, im1, il2, im2, ia, ierr
    real(kind=DP) :: fac,mpi
    eohxc_paw = 0.d0
    do ispin = 1, nspin, af+1
       do it = 1, ntyp
          if(ipaw(it)==0) cycle
          do lmt1 = 1, ilmt(it)
! === DEBUG by tkato 2012/11/08 ================================================
             il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
! ==============================================================================
             do lmt2 = lmt1, ilmt(it)
! === DEBUG by tkato 2012/11/08 ================================================
                il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
!!!                if( il1 /= il2 .or. im1 /= im2) cycle      !! ASMS 2015/02/06
! ==============================================================================
!xocl spread do/ind_katm
                do ia = ista_atm,iend_atm
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
    if(nrank_g > 1) then
       call mpi_allreduce(eohxc_paw,mpi,1,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       eohxc_paw = mpi
    end if
    eohxc_paw = eohxc_paw*(af+1)
  end subroutine get_xc_and_HE_of_old_CD_paw_3D
!===============================================================================

  subroutine get_hartree_energy_paw_3D
    use m_Parallelization,     only :mpi_ke_world
    integer :: ispin, it,lmt1, lmt2, il1, im1, il2, im2, ia, ierr
    real(kind=DP) :: fac, mpi

    ehartr_paw = 0.d0
    do ispin = 1, nspin, af+1
       do it = 1, ntyp
          if(ipaw(it)==0) cycle
          do lmt1 = 1, ilmt(it)
! === DEBUG by tkato 2012/11/07 ================================================
             il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
! ==============================================================================
             do lmt2 = lmt1, ilmt(it)
! === DEBUG by tkato 2012/11/07 ================================================
                il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
!!!                if( il1 /= il2 .or. im1 /= im2) cycle    !! ASMS 2015/02/06
! ==============================================================================
!xocl spread do/ind_katm
                do ia = ista_atm, iend_atm
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
    if(nrank_g > 1) then
       call mpi_allreduce(ehartr_paw,mpi,1,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       ehartr_paw = mpi
    end if
    ehartr_paw = ehartr_paw*(af+1)*0.5d0
  end subroutine get_hartree_energy_paw_3D
!===============================================================================

!!$  function m_TE_tell_band_energy_3D(nfout,kv3)
!!$    real(kind=DP)::  m_TE_tell_band_energy_3D
!!$    integer, intent(in) :: nfout,kv3
!!$    integer :: ib, ik
!!$    call sum_eigenvalues_3D(kv3)
!!$    m_TE_tell_band_energy_3D = eband
!!$    if(ipri >= 2) then
!!$       if(map_k(ik) == myrank_k) then
!!$          write(nfout,'(" -- sum_eigenvalues (m_TE_tell_band_energy) --")')
!!$          write(nfout,'("   -- num_extra_bands = ",i8)') num_extra_bands
!!$          do ik = 1, kv3
!!$             write(nfout,'("  ik = ",i5)') ik
!!$             write(nfout,'(8f8.4)') (eko_l(ib,ik),ib=1,np_e) ! MPI
!!$          end do
!!$       end if
!!$    endif
!!$  end function m_TE_tell_band_energy_3D

!!$  function m_TE_tell_extended_band_energy_3D(nfout,kv3)
!!$    real(kind=DP)::  m_TE_tell_extended_band_energy_3D
!!$    integer, intent(in) :: nfout,kv3
!!$    integer :: ib, ik
!!$    call sum_wholeeigenvalues_3D(nfout,kv3)
!!$    m_TE_tell_extended_band_energy_3D = eband_extendedrange
!!$    if(ipri >= 2) then
!!$       if(map_k(ik) == myrank_k) then
!!$          write(nfout,'(" -- sum_eigenvalues (m_TE_tell_band_energy) --")')
!!$          write(nfout,'("   -- num_conduction_bands_lmm = ",i8)') num_conduction_bands_lmm
!!$          do ik = 1, kv3
!!$             write(nfout,'("  ik = ",i5)') ik
!!$             write(nfout,'(8f8.4)') (eko_l(ib,ik),ib=1,np_e) ! MPI
!!$          end do
!!$       end if
!!$    endif
!!$  end function m_TE_tell_extended_band_energy_3D

  subroutine sum_eigenvalues(kv3)
   use m_Parallelization,     only : neg_g                   &
  &                                , mpi_kg_world
    integer, intent(in) :: kv3
    integer             :: ib, ik
    eband = 0.d0
    do ik =1, kv3, af+1
       if(map_k(ik) /= myrank_k) cycle
       do ib = 1, np_e
          if(neg_g(ib) > neg - num_extra_bands) cycle
          eband = eband + eko_l(ib,ik)*occup_l(ib,ik)
       end do
    end do
    call mpi_allreduce(MPI_IN_PLACE,eband,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
    eband = 2*eband/kv3 * (af+1)
  end subroutine sum_eigenvalues

  subroutine sum_wholeeigenvalues(nfout,kv3)
   use m_Parallelization,     only : neg_g                   &
  &                                , mpi_kg_world
    integer, intent(in) :: nfout,kv3
    integer             :: ib, ik
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
    do ik =1, kv3, af+1
       weight = weight + qwgt(ik)
       if(map_k(ik) /= myrank_k) cycle ! MPI
       do ib = 1, np_e                 ! MPI
          if(ik==1 .and. ib == 1) weight0 = weight0 + occup_l(1,ik)
!!$          if(nrvf_ordr(ib,ik) > neg - num_extra_bands) cycle
!!$       if(nrvf_ordr(ib,ik) > iband) cycle
          if(neg_g(ib) > iband) cycle
          eband_extendedrange = eband_extendedrange + eko_l(ib,ik)*qwgt(ik)
       end do
    end do
    if(nrank_e > 1) then
       call mpi_allreduce(eband_extendedrange,eband_mpi,1,mpi_double_precision, &
      &                   mpi_sum,mpi_kg_world,ierr) !MPI
       eband_extendedrange = eband_mpi                  ! MPI
       call mpi_allreduce(weight0, eband_mpi,1,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
       weight0 = eband_mpi
    end if
    eband_extendedrange = (3.0-nspin) * eband_extendedrange * (af+1)
    if(ipri >= 1) then
       weight = weight * (af+1)
       write(nfout,'("   -- num_conduction_bands_lmm = ",i8)') num_conduction_bands_lmm
       write(nfout,'(" weight (sum_wholeeigenvalues) = ", f8.4, " iband = ",i10)') weight,iband
       write(nfout,'(" weight0 (occup_l)             = ", f8.4)') weight0
       write(nfout,'(" eband_extendedrange = ", f16.8, "   eband = ",f16.8)') eband_extendedrange, eband
    end if
  end subroutine sum_wholeeigenvalues

end module m_Total_Energy
