!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 603 $)
!
!  SUBROUINE:  Preparation_for_mpi
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
subroutine Preparation_for_mpi(prepare_communicators)
! $Id: Preparation_for_mpi.F90 603 2020-04-07 03:25:30Z jkoga $
!                           @(#)Preparation_for_mpi.F90 1.10 03/02/19 00:49:14
  use m_Const_Parameters,     only : INITIAL, CONTINUATION &
       &                           , FIXED_CHARGE, FIXED_CHARGE_CONTINUATION,OFF &
       &                           , ALL_AT_ONCE, ON, DP, COORDINATE_CONTINUATION &
       &                           , ONE_BY_ONE
  use m_Kpoints,              only : kv3, kv3_ek
  use m_Files,                only : nfout
  use m_Control_Parameters,   only : ipriparallel,nspin,neg,printable,ngnode_nbmx &
       &                           , flag_mpi_g_dot_r,flag_mpi_g_dot_r_k &
       &                           , icond, ekmode, fixed_charge_k_parallel, sw_rsb, neg_is_given &
       &                           , m_CtrlP_flag_mpi_G_dot_R
  use m_PlaneWaveBasisSet,    only : kg1, kgpm, nbmx, iba, kgp, kgp_reduced
#ifndef PARAMSET
  use m_Charge_Density,       only : m_CD_alloc_chgq
  use m_XC_Potential,         only : m_XC_alloc_vxc
  use m_Parallelization,      only : m_Parallel_init_mpi_mix, m_Parallel_init_mpi_snl &
       &                           , m_Parallel_init_mpi_atm, m_Parallel_init_mpi_atm2 &
       &                           , m_Parallel_init_mpi_iba, m_Parallel_init_mpi_elec &
       &                           , m_Parallel_init_mpi_kv3_ek
  use m_Parallelization,      only : m_Parallel_init_mpi_gga &
       &                           , m_Parallel_init_mpi_nbmx &
       &                           , m_Parallel_init_mpi_ffth
  use m_FFT,                  only : nfftp, nfftps, nfft
  use m_Ionic_System,         only : natm, natm2, m_IS_alloc_fxyzew
  use m_Ionic_System,         only : m_IS_alloc_zfm3
  use m_Force,                only : m_Force_alloc
#ifdef __EDA__
  use m_XC_Potential,         only : m_XC_alloc_vxc, m_XC_alloc_exc_on_a_grid
  use m_Ionic_System,         only : m_IS_alloc_eewald_per_atom, m_IS_alloc_zfm3_EDA
  use m_PseudoPotential,      only : m_PP_alloc_PP_per_atom_etc
#endif
#endif

! ============================== added by K. Tagami ================== 11.0&13.0XX
  use m_Control_Parameters,   only : noncol, ndim_spinor, sw_calc_ekin_density
  use m_KineticEnergy_Density,  only : m_KE_alloc_ekin_density
! ==================================================================== 11.0&13.0XX

#ifdef __EDA__
  use m_Control_Parameters,   only : sw_eda
#endif

  implicit none

  integer, intent(in) :: prepare_communicators
  integer             :: lsize
  real(kind=DP), allocatable, dimension(:,:) :: dfft_l

#ifndef PARAMSET
  call m_IS_alloc_fxyzew()
  call m_Force_alloc()
  call m_IS_alloc_zfm3()
  call m_CD_alloc_chgq()
  call m_XC_alloc_vxc()

! === KT_add ==== 13.0XX
  if ( sw_calc_ekin_density == ON ) call m_KE_alloc_ekin_density
! =============== 13.0XX

  if(prepare_communicators==ON)then
! ================================ modified by K. Tagami =================== 11.0
!     call m_Parallel_init_mpi_elec(nfout,ipriparallel,printable,neg,kv3,nspin,kg1)
!     call m_Parallel_init_mpi_snl(nfout,ipriparallel,printable,nspin)
! -
     if ( noncol ) then
!!$        if(neg_is_given) &
             call m_Parallel_init_mpi_elec( nfout, ipriparallel, printable, neg, &
             &                             kv3, ndim_spinor, kg1 )
        call m_Parallel_init_mpi_snl( nfout, ipriparallel, printable, ndim_spinor )
     else
!!$        if(neg_is_given) &
             call m_Parallel_init_mpi_elec( nfout, ipriparallel, printable, neg, &
             &                             kv3, nspin, kg1 )
        call m_Parallel_init_mpi_snl( nfout, ipriparallel, printable, nspin )
     endif
! ========================================================================== 11.0 
     call m_Parallel_init_mpi_mix(nfout,ipriparallel,printable,kgpm)
     call m_Parallel_init_mpi_atm(nfout,ipriparallel,printable,natm)
     call m_Parallel_init_mpi_atm2(nfout,ipriparallel,printable,natm2)

!!!  call m_CtrlP_flag_mpi_G_dot_R(nfout,nbmx) ! -> flag_mpi_g_dot_r
     call m_Parallel_init_mpi_nbmx(nfout,ipriparallel,printable,nbmx,kg1,ngnode_nbmx,flag_mpi_g_dot_r,flag_mpi_g_dot_r_k)
     call m_Parallel_init_mpi_gga(nfout,ipriparallel,printable,nfftp,nfftps)
     if(sw_rsb==ON) call m_Parallel_init_mpi_ffth(nfout,ipriparallel,printable,nfft)

  endif

#ifdef __EDA__
  if (sw_eda == ON) then
! -----  ascat starts modifying  -----
  call m_XC_alloc_exc_on_a_grid()
  call m_IS_alloc_zfm3_EDA(natm)
  call m_IS_alloc_eewald_per_atom(natm)
  call m_PP_alloc_PP_per_atom_etc
! -----  ascat ceases modifying  -----
  endif
#endif

  if((icond == INITIAL .or. icond == CONTINUATION .or.  icond==COORDINATE_CONTINUATION) &
       & .or.((icond==FIXED_CHARGE.or.icond==FIXED_CHARGE_CONTINUATION).and.ekmode==OFF &
       & .and. fixed_charge_k_parallel == ALL_AT_ONCE)) then
     call m_Parallel_init_mpi_iba(nfout,ipriparallel,printable,kv3,iba) ! -> np_g1k, mp_g1k
  end if
!!!!$!BRANCH_P ORG_Parallel
  if((icond==FIXED_CHARGE .or. icond==FIXED_CHARGE_CONTINUATION) .and. &
    & fixed_charge_k_parallel==ONE_BY_ONE) &
    & call m_Parallel_init_mpi_kv3_ek(nfout,ipriparallel,printable,kv3_ek,nspin)
!!!!$!BRANCH_P_END ORG_Parallel



#endif

end subroutine Preparation_for_mpi

subroutine Preparation_for_mpi_ek
! $Id: Preparation_for_mpi.F90 603 2020-04-07 03:25:30Z jkoga $
!                           @(#)Preparation_for_mpi.F90 1.10 03/02/19 00:49:14
  use m_Kpoints,              only : kv3
  use m_Files,                only : nfout
  use m_Control_Parameters,   only : ipriparallel, printable
  use m_PlaneWaveBasisSet,    only : iba
  use m_Parallelization,      only : m_Parallel_init_mpi_iba

  call m_Parallel_init_mpi_iba(nfout,ipriparallel,printable,kv3,iba) !  -> np_g1k, mp_g1k

end subroutine Preparation_for_mpi_ek

subroutine Preparation_for_mpi_PAW()
  use m_Files,                only : nfout
  use m_Ionic_System,         only : natm
  use m_PseudoPotential,      only : mmesh, flg_paw
end subroutine Preparation_for_mpi_PAW

