!#undef NEC_TIMER
!#ifdef NEC_TIMER
!#  define START_TIMER(a) call start_timer(a)
!#  define STOP_TIMER(a)  call stop_timer(a)
#ifdef FJ_TIMER_RMM
#  define START_TIMER(a) call timer_sta(a)
#  define STOP_TIMER(a)  call timer_end(a)
#else
#  define START_TIMER(a)
#  define STOP_TIMER(a)
#endif
#ifdef __TIMER_COMM__
#   define __TIMER_COMM_START_w_BARRIER(str,a)   call timer_barrier(str) ;   call timer_sta(a)
#   define __TIMER_COMM_START(a)       call timer_sta(a)
#   define __TIMER_COMM_STOP(a)        call timer_end(a)
#else
#   define __TIMER_COMM_START_w_BARRIER(str,a)
#   define __TIMER_COMM_START(a)
#   define __TIMER_COMM_STOP(a)
#endif
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 581 $)
!
!  MODULE: m_ES_WF_by_RMM
!
!  AUTHOR(S): T. Yamasaki, T. Uda   August/20/2003
!
!  FURTHER MODIFICATION: T. Yamasaki, January/13/2004, June/04/2005
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
#ifndef SX
#define DGEMM__       DGEMM
#endif

#ifndef NO_NONLOCAL_RMM_DGEMM
#define NONLOCAL_RMM_DGEMM
#endif

module m_ES_WF_by_RMM
! This module was coded by T. Yamasaki and T. Uda, 1998-2003
!
! previous Id: m_ES_WF_by_RMM.f90,v 1.6 2003/09/16 03:15:32 yamasaki Exp
! $Id: m_ES_WF_by_RMM.F90 581 2018-08-01 08:38:42Z jkoga $
!
  use m_Electronic_Structure,only : zaj_l, afft,bfft, eko_l,vnlph_l,vlhxcQ &
#ifdef SAVE_FFT_TIMES
 &                                , status_saved_phifftr &
#endif
       &                          , m_ES_alloc_fft_related &
       &                          , m_ES_dealloc_fft_related &
       &                          , m_ES_Vlocal_in_Rspace &
       &                          , m_ES_eigen_values_for_each_k &
       &                          , m_ES_WF_in_Rspace, m_ES_decide_precon_factor &
       &                          , m_ES_wd_zaj_small_portion,m_ES_init_chgsoft, m_ES_dealloc_chgsoft &
       &                          , fsr_l,fsi_l, hlocphi_l, vtau_l
  use m_ES_nonlocal,         only : m_ES_Vnonlocal_W,m_ES_betar_dot_WFs_4_each_k
  use m_ES_ortho            ,only : m_ES_MGS_4_each_k
  use m_ES_WF_by_SDorCG,     only : m_ESsd_diff_WFs_to_zaj_old
  use m_ES_WF_by_SDorCG,     only : m_ESsd_copy_zaj_old_to_zaj, m_ESsd_copy_phi_to_zaj_old
  use m_NonLocal_Potential,  only : snl
  use m_PlaneWaveBasisSet,   only : kg1,kgp,nbmx,iba, igf, nbase, ngabc &
       &                          , m_pwBS_find_min_max_G, m_pwBS_kinetic_energies
  use m_PseudoPotential,     only : ilmt, lmtt, lmta, ltp, mtp &
       &                          , q, dion, nlmta1, nlmta2, modnrm, nac, fqwei, nlmta &
       &                          , m_PP_include_vanderbilt_pot &
       &                          , ipaw,dion_paw
  use m_Kpoints,             only : kv3,vkxyz, k_symmetry
  use m_Ionic_System,        only : ntyp,natm,ityp,iwei,pos
  use m_FFT,                 only : fft_box_size_WF
  use m_FFT,                 only : m_FFT_Vlocal_W, m_FFT_WF
  use m_Timing,              only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters,  only : nspin,iprirmm,imGSrmm, rr_Critical_Value, printable &
       &                          , rmm_save_memory_mode, rmm_precal_phase_matm,kimg,neg,af &
       &                          , sw_retard_eigval_evaluation,sw_precalculate, potential_update &
       &                          , sw_serial_fft &
       &                          , sw_gep &
       &                          , sw_keep_hloc_phi &
       &                          , sw_reduce_fft_for_charge &
#ifdef SAVE_FFT_TIMES
       &                          , sw_hybrid_functional, sw_fef, sw_save_fft
#else
       &                          , sw_hybrid_functional, sw_fef
#endif
  use m_IterationNumbers,    only : iteration_electronic, iteration_rmm_start, iteration, nk_in_the_process
  use m_Const_Parameters,    only : DP,CMPLDP,DIRECT,ON,OFF,SKIP,EXECUT,RMM2,RMM2P,RMM3,PAI2 &
       &                          , ORTHONORMALIZATION, NORMALIZATION, OVER,UNDER, GAMMA &
       &                          , SmallestPositiveNumber, ELECTRON, OLD
  use m_Parallelization,     only : MPI_CommGroup &
       &                          , myrank_k,map_k,ista_e,iend_e,istep_e,np_e,map_z,ierr,npes &
       &                          , ista_k,iend_k, mpi_spin_group, ista_spin, iend_spin, nrank_s
#ifdef NEC_TIMER
  use nec_timer
#endif
  use m_ES_ExactExchange,    only : m_ES_EXX_eigenvalue_for_each_k &
       &                          , m_ES_EXX_gather_valence_states, m_ES_EXX_cp_eigenvalue &
       &                          , m_ES_Vexx_add_vexx &
       &                          , m_ES_Vexx_W
  use m_FiniteElectricField, only : m_FEF_build_grad, m_FEF_add_grad_to_vnlph &
       &                          , m_FEF_add_grad_to_vnlph_RMM

! ====================================== added by K. Tagami ============= 11.0
  use m_Control_Parameters,  only : noncol, ndim_spinor, ndim_chgpot, SpinOrbit_mode, &
       &                            sw_hubbard
  use m_Const_Parameters,    only : Neglected, BuiltIn
  use m_FFT,                 only : nfft, m_FFT_Vlocal_W_noncl, &
       &                            m_FFT_alloc_WF_work, m_FFT_dealloc_WF_work

  use m_PseudoPotential,     only : q_noncl, dion_scr_noncl, fqwei_noncl
  use m_ES_ortho            ,only : m_ES_MGS_4_each_k_noncl
  use m_Electronic_Structure, only :  m_ES_eigen_vals_each_k_noncl, &
       &                              m_ES_Vlocal_in_Rspace_noncl, &
       &                              m_ES_sort_eigen_vals_noncl
! ======================================================================= 11.0

  use m_Control_Parameters,  only : use_metagga, vtau_exists
  use m_ES_WF_by_SDorCG,     only : m_ES_contrib_kindens_to_vnlph, &
       &                            m_ES_kindens_to_vnlph_ib, m_ES_kindens_to_vnlph_ib2
  use mpi

  implicit none

  real(kind=DP), private,allocatable,dimension(:,:)     :: zfc1,zfc2,zfc3
  real(kind=DP), private,allocatable,dimension(:,:)     :: zfs1,zfs2,zfs3
  real(kind=DP), private,allocatable,dimension(:,:)     :: phasec,phases
  integer,       private,allocatable,dimension(:,:,:)   :: nglist
  integer,       private,pointer,dimension(:,:)     :: nngabc
  integer,       private,pointer,dimension(:)       :: newp
  real(kind=DP), private,allocatable,dimension(:,:,:) :: bWr, bWi
  real(kind=DP), private,pointer,dimension(:)       :: ekin0
  real(kind=DP), private,pointer,dimension(:,:,:,:) ::      phi  ! d(kg1,kimg,n,0:nrmm-2)
  real(kind=DP), private,allocatable,dimension(:,:,:,:) ::  Rphi ! d(kg1,kimg,n,0:nrmm-2)
  real(kind=DP), private,allocatable,dimension(:,:,:,:) ::  rr_e, psp_e  ! d(0:nrmm-1,0:nrmm-1,n,kimg),
  !                                                      ( n={1(rmm_save_memory_mode)|np_e(otherwise)})

! ============================== added by K. Tagami ===================== 11.0
  real(kind=DP), private,pointer,dimension(:,:,:,:,:) ::      phi_noncl
                                                 ! d(kg1,kimg,n,0:nrmm-2,ndim_spinor)
  real(kind=DP), private,allocatable,dimension(:,:,:,:,:) ::  Rphi_noncl
                                                 ! d(kg1,kimg,n,0:nrmm-2,ndim_spinor)
  real(kind=DP), private,allocatable,dimension(:,:,:,:) :: bWr_noncl, bWi_noncl
! ======================================================================= 11.0

  integer,       private,        dimension(2)       :: pm = (/1,-1/)

  real(kind=DP), private                            :: rr_avr = 0.d0
  integer, private,allocatable,dimension(:,:)     :: rr_is_over_or_under

  integer,       private,save                       :: ng_max = 0, ng_min = 0 !

!   1. m_ESrmm_reset_ng_maxmin     <-(Renewal_of_WaveFunctions)
!   2. m_ESrmm_reset_r_norm_flag    <-(Renewal_of_WaveFunctions)
!   3. alloc_zfc_zfs     <-(9)
!   4. dealloc_zfc_zfs   <-(9)
!   5. alloc_phasecs_bW_phi_Rphi   <-(9)
!   6. dealloc_phasecs_bW_phi_Rphi <-(9)
!   7. alloc_bW_phi_Rphi           <-(9)
!      - sub_error1, - sub_error2
!   8. dealloc_bW_phi_Rphi         <-(9)
!      - sub_error1
!   9. m_ESrmm_renew_WF             <-(Renewal_of_WaveFunctions)
!      - orthonorm_or_norm, - zajold2zaj_phi2zaj_old_all, - zajold2zaj_phi2zaj_old
!      - rr_avr_final,
!      - what_is_the_dimension_of_rmm, - rmm1, - rmm_n, - evolve_WF_using_Residuals,
!      - wd_rr_and_psp, - wd_alpha, - rmm_n_uda, - rg_or_cg,
!      - solveAYB, - LUdecomp, - rr_and_psp_from_Rphi_and_phi, - phsphn,
!      - normalization_of_rr_and_psp, - set_phase_pointer, - zf_listing
!  10. Vnonlocal_W_RMM       <-(9)
!      - alloc_zfsincos_arai, - dealloc_zfsincos_arai, - alloc_scssqcqs,
!      - dealloc_scssqcqs,    - calc_phase_RMM,  - Vnonlocal_W_part_sum_over_lmt1
!      - sumset_rmm,   - add_vnlph_l_without_eko_part1, - add_vnlph_l_with_eko_part1
!  11. Vnonlocal_W_RMMn      <-(9)
!      - alloc_scssqcqs_x, - dealloc_scssqcqs_x, - find_max_n_ialist0, - sumset_rmm_all3
!      - add_vnlph_l_without_eko_part3, -add_vnlph_l_with_eko_part3
!      - Vnonlocal_W_part_sum_over_lmt1b

!  include 'mpif.h'

contains

  subroutine m_ESrmm_reset_ng_maxmin()
    ! this subroutine is called when
    ! (icond == FIXED_CHARGE or icond == FIXED_CHARGE_CONTINUATION)
    ! and iteration_rmm_start is set.
    ng_max = 0; ng_min = 0
  end subroutine m_ESrmm_reset_ng_maxmin

  subroutine m_ESrmm_dealloc_r_norm_flag
     if(allocated(rr_is_over_or_under)) deallocate(rr_is_over_or_under)
  end subroutine m_ESrmm_dealloc_r_norm_flag

  subroutine m_ESrmm_reset_r_norm_flag
    if(.not.allocated(rr_is_over_or_under)) allocate(rr_is_over_or_under(np_e,ista_k:iend_k))
    rr_is_over_or_under = OVER
  end subroutine m_ESrmm_reset_r_norm_flag

  subroutine alloc_zfc_zfs(nfout,n_min,n_max,matm_t,natm,nrmm,matm)
    integer, intent(inout), dimension(3) :: n_min, n_max
    integer, intent(in)                  :: nfout, matm_t, natm,nrmm
    integer, intent(out) :: matm

    matm = matm_t
    if(matm > natm) matm = natm
    if(matm < 0) matm = 0

! ========================================== modified by K. Tagami =========== 11.0
!    call m_pwBS_find_min_max_G(nfout,nspin,ng_max,ng_min,n_max,n_min)
!
    if ( noncol ) then
       call m_pwBS_find_min_max_G(nfout,ndim_spinor,ng_max,ng_min,n_max,n_min)
    else
       call m_pwBS_find_min_max_G(nfout,nspin,ng_max,ng_min,n_max,n_min)
    endif
! ============================================================================ 11.0

    if(iprirmm >= 2) then
       write(nfout,'(" !!rmm kg1, matm = ",2i6)') kg1,matm
       write(nfout,'(" !! matm+1, natm = ",2i8," <<alloc_zfc_zfs>>")') matm+1,natm
    end if

    if(matm+1 <= natm) then
       allocate(zfc1(n_min(1):n_max(1),matm+1:natm))
       allocate(zfc2(n_min(2):n_max(2),matm+1:natm))
       allocate(zfc3(n_min(3):n_max(3),matm+1:natm))
       allocate(zfs1(n_min(1):n_max(1),matm+1:natm))
       allocate(zfs2(n_min(2):n_max(2),matm+1:natm))
       allocate(zfs3(n_min(3):n_max(3),matm+1:natm))
    end if
  end subroutine alloc_zfc_zfs

  subroutine dealloc_zfc_zfs()
    if(iprirmm >= 2) write(6,'(" rmm_save_memory_mode = ",i6," <<dealloc_zfc_zfs>>")') rmm_save_memory_mode
    if(allocated(zfs1)) then
       deallocate(zfs1);    deallocate(zfs2);    deallocate(zfs3)
       deallocate(zfc1);    deallocate(zfc2);    deallocate(zfc3)
    end if
  end subroutine dealloc_zfc_zfs

  subroutine alloc_phasecs_bW_phi_Rphi(nrmm,n_min,n_max,matm)
    integer, intent(in) :: nrmm
    integer, intent(in) :: n_min(3),n_max(3),matm

    allocate(nngabc(kg1,3))
    allocate(newp(kg1))
    allocate(phi(kg1,kimg,1,0:nrmm-2));  phi  = 0.d0
    allocate(Rphi(kg1,kimg,1,0:nrmm-2)); Rphi = 0.d0
    allocate(rr_e(0:nrmm-1,0:nrmm-1,1,kimg))
    allocate(psp_e(0:nrmm-1,0:nrmm-1,1,kimg))
    allocate(phasec(kg1,matm)); allocate(phases(kg1,matm))
    allocate(nglist(n_min(1):n_max(1),n_min(2):n_max(2),n_min(3):n_max(3)))
    allocate(bWr(nlmta,0:nrmm-1,1)); allocate(bWi(nlmta,0:nrmm-1,1))
  end subroutine alloc_phasecs_bW_phi_Rphi

! ====================================== added by K. Tagami =============== 11.0
  subroutine alloc_phasecs_bW_phietc_noncl( nrmm,n_min,n_max,matm )
    integer, intent(in) :: nrmm
    integer, intent(in) :: n_min(3),n_max(3),matm

    allocate(nngabc(kg1,3))
    allocate(newp(kg1))
    allocate(phi_noncl(kg1,kimg,1,0:nrmm-2,ndim_spinor));  phi_noncl  = 0.d0
    allocate(Rphi_noncl(kg1,kimg,1,0:nrmm-2,ndim_spinor)); Rphi_noncl = 0.d0

    allocate(rr_e(0:nrmm-1,0:nrmm-1,1,kimg))
    allocate(psp_e(0:nrmm-1,0:nrmm-1,1,kimg))
    allocate(phasec(kg1,matm)); allocate(phases(kg1,matm))
    allocate(nglist(n_min(1):n_max(1),n_min(2):n_max(2),n_min(3):n_max(3)))

    allocate(bWr_noncl(nlmta,0:nrmm-1,1,ndim_spinor))
    allocate(bWi_noncl(nlmta,0:nrmm-1,1,ndim_spinor))
  end subroutine alloc_phasecs_bW_phietc_noncl
! ============================================================================ 11.0

  subroutine dealloc_phasecs_bW_phi_Rphi()
    deallocate(bWi); deallocate(bWr)
    deallocate(nglist)
    deallocate(phasec); deallocate(phases)
    deallocate(psp_e)
    deallocate(rr_e)
    deallocate(Rphi)
    deallocate(phi)
    deallocate(newp)
    deallocate(nngabc)
  end subroutine dealloc_phasecs_bW_phi_Rphi

! ================================= added by K. Tagami ================= 11.0
  subroutine dealloc_phasecs_bW_phietc_noncl()
    deallocate(bWi_noncl); deallocate(bWr_noncl)
    deallocate(nglist)
    deallocate(phasec); deallocate(phases)
    deallocate(psp_e)
    deallocate(rr_e)
    deallocate(Rphi_noncl)
    deallocate(phi_noncl)
    deallocate(newp)
    deallocate(nngabc)
  end subroutine dealloc_phasecs_bW_phietc_noncl
! ====================================================================== 11.0

  subroutine alloc_bW_phi_Rphi(nrmm)
    integer, intent(in) :: nrmm
    integer :: istat
    character*4 :: name
    if(iprirmm >= 3) then
       write(6,'(" !! nlmta, nrmm, np_e = ",3i10," <<alloc_bW_phi_Rphi>>")') nlmta,nrmm,np_e
       write(6,'(" !! kg1,   kimg, np_e = ",3i10," <<alloc_bW_phi_Rphi>>")') kg1,kimg,np_e
    end if

    allocate(bWr(nlmta,0:nrmm-1,np_e),stat=istat)
    if(istat /= 0) then
       name = "bWr "
       call sub_error1
    end if
    bWr = 0.d0

    allocate(bWi(nlmta,0:nrmm-1,np_e),stat=istat)
    if(istat /= 0) then
       name = "bWi"
       call sub_error1
    end if
    bWi = 0.d0

    allocate(Rphi(kg1,kimg,np_e,0:nrmm-2),stat=istat)
    if(istat /= 0) then
       name = "Rphi"
       call sub_error2
       call sub_error1
    else
       if(iprirmm >= 2) write(6,*) 'Allocation success for Rphi <<alloc_bW_phi_Rphi>>'
    end if
    Rphi = 0.d0

    allocate(phi(kg1,kimg,np_e,0:nrmm-2),stat=istat)
    if(istat /= 0) then
       name = "phi "
       call sub_error2
       call sub_error1
    else
       if(iprirmm >= 2) write(6,*) 'Allocation success for Rphi <<alloc_bW_phi_Rphi>>'
    end if
    phi = 0.d0

    allocate(rr_e(0:nrmm-1,0:nrmm-1,np_e,kimg));    rr_e = 0.d0
    allocate(psp_e(0:nrmm-1,0:nrmm-1,np_e,kimg));   psp_e = 0.d0
  contains
    subroutine sub_error1
      write(6,'(" Allocation error for ",a4," <<m_ES_WF_by_RMM.alloc_bW_phi_Rphi>>")') name
      write(6,*) 'stat = ',istat
      call phase_error_with_msg(6, 'Allocation error <<m_ES_WF_by_RMM.alloc_bW_phi_Rphi>>',__LINE__,__FILE__)
    end subroutine sub_error1
    subroutine sub_error2
      write(6,'(" required allocation size = ",i10," kg1,kimg,np_e,nrmm = ",4i8)') kg1*kimg*np_e*nrmm,kg1,kimg,np_e,nrmm
    end subroutine sub_error2

  end subroutine alloc_bW_phi_Rphi

! ============================================= added by K. Tagami =========== 11.0
  subroutine alloc_bW_phi_Rphi_noncl(nrmm)
    integer, intent(in) :: nrmm
    integer :: istat
    character*4 :: name

    if (iprirmm >= 3) then
       write(6,'(" !! nlmta, nrmm, np_e = ",3i10," <<alloc_bW_phi_Rphi_noncl>>")') &
	&       nlmta,nrmm,np_e
       write(6,'(" !! kg1,   kimg, np_e = ",3i10," <<alloc_bW_phi_Rphi_noncl>>")') &
	&       kg1,kimg,np_e
    end if

    allocate(bWr_noncl(nlmta,0:nrmm-1,np_e,ndim_spinor),stat=istat)
    if(istat /= 0) then
       name = "bWr_noncl ";   call sub_error1
    end if
    bWr_noncl = 0.d0

    allocate(bWi_noncl(nlmta,0:nrmm-1,np_e,ndim_spinor),stat=istat)
    if(istat /= 0) then
       name = "bWi_noncl ";   call sub_error1
    end if
    bWi_noncl = 0.d0

    allocate( Rphi_noncl(kg1,kimg,np_e,0:nrmm-2,ndim_spinor),stat=istat )
    if (istat /= 0) then
       name = "Rphi_noncl";    call sub_error2;    call sub_error1
    else
       if (iprirmm >= 2) then
	  write(6,*) 'Allocation success for Rphi_noncl <<alloc_bW_phi_Rphi_noncl>>'
       endif
    end if
    Rphi_noncl = 0.d0

    allocate( phi_noncl(kg1,kimg,np_e,0:nrmm-2,ndim_spinor),stat=istat )
    if (istat /= 0) then
       name = "phi_noncl ";    call sub_error2;     call sub_error1
    else
       if(iprirmm >= 2) then
	  write(6,*) 'Allocation success for Rphi_noncl <<alloc_bW_phi_Rphi_noncl>>'
       endif
    end if
    phi_noncl = 0.d0

    allocate(rr_e(0:nrmm-1,0:nrmm-1,np_e,kimg));    rr_e = 0.d0
    allocate(psp_e(0:nrmm-1,0:nrmm-1,np_e,kimg));   psp_e = 0.d0

  contains

    subroutine sub_error1
      write(6,'(" Allocation error for ",a10, &
	&    " <<m_ES_WF_by_RMM.alloc_bW_phi_Rphi_noncl>>")') name
      write(6,*) 'stat = ',istat
      call phase_error_with_msg(6,  'Allocation error <<m_ES_WF_by_RMM.alloc_bW_phi_Rphi_noncl>>',__LINE__,__FILE__)
    end subroutine sub_error1

    subroutine sub_error2
      write(6,'(" required allocation size = ",i10," kg1,kimg,np_e,nrmm = ",4i8)') kg1*kimg*np_e*nrmm,kg1,kimg,np_e,nrmm
    end subroutine sub_error2

  end subroutine alloc_bW_phi_Rphi_noncl
! ===================================================================  11.0

  subroutine dealloc_bW_phi_Rphi()
    integer :: istat
    character*4 :: name
    deallocate(psp_e)
    deallocate(rr_e)
    deallocate(phi, stat=istat)
    if(istat /= 0) then
       name = "phi "
       call sub_error1
    end if
    deallocate(Rphi, stat=istat)
    if(istat /= 0) then
       name = "Rphi"
       call sub_error1
    end if
    deallocate(bWi, stat=istat)
    if(istat /= 0) then
       name = "bWi "
       call sub_error1
    end if
    deallocate(bWr, stat=istat)
    if(istat /= 0) then
       name = "bWr "
       call sub_error1
    end if
  contains
    subroutine sub_error1
      write(6,'(" Deallocation error for ",a4,"<<m_ES_WF_by_RMM.dealloc_bW_phi_Rphi>>")') name
      write(6,*) 'stat =', istat
      call phase_error_with_msg(6, 'Deallocation error <<m_ES_WF_by_RMM.dealloc_bW_phi_Rphi>>',__LINE__,__FILE__)
    end subroutine sub_error1
  end subroutine dealloc_bW_phi_Rphi

!  ========================================= added by K. Tagami ============== 11.0
  subroutine dealloc_bW_phi_Rphi_noncl()
    integer :: istat
    character*4 :: name
    deallocate(psp_e)
    deallocate(rr_e)

    deallocate( phi_noncl, stat=istat )
    if (istat /= 0) then
      name = "phi_noncl ";         call sub_error1
    end if
    deallocate( Rphi_noncl, stat=istat )
    if (istat /= 0) then
      name = "Rphi_noncl";         call sub_error1
    end if

    deallocate( bWi_noncl, stat=istat)
    if(istat /= 0) then
       name = "bWi_noncl "
       call sub_error1
    end if
    deallocate( bWr_noncl, stat=istat)
    if(istat /= 0) then
       name = "bWr_noncl "
       call sub_error1
    end if

  contains

    subroutine sub_error1
      write(6,'(" Deallocation error for ",a10, &
	&        "<<m_ES_WF_by_RMM.dealloc_bW_phi_Rphi_noncl>>")') name
      write(6,*) 'stat =', istat
      call phase_error_with_msg(6, 'Deallocation error <<m_ES_WF_by_RMM.dealloc_bW_phi_Rphi_noncl>>',__LINE__,__FILE__)
    end subroutine sub_error1

  end subroutine dealloc_bW_phi_Rphi_noncl
! ===================================================================== 11.0

  subroutine m_ESrmm_renew_WF(nfout,isolver,precon,dtim)
    integer,       intent(in) :: nfout,isolver,precon
    real(kind=DP), intent(in) :: dtim

    !  R   := -(H-eS)
    !  zaj_l(:,i,k,1:kimg) :=  | \Phi_{k,i} >
    !  phi (:,1:kimg,0)    :=  | \Phi_{k,i} >
    !  Rphi(:,1:kimg,0)    := R|\Phi_{k,i}>
    !  phi (:,1:kimg,1)    := KR |\Phi_{k,i}>
    !  Rphi(:,1:kimg,1)    := RKR|\Phi_{k,i}>
    ! when isolver == RMM3,
    !  phi (:,1:kimg,2)    := KRKR |\Phi_{k,i}>
    !  Rphi(:,1:kimg,2)    := RKRKR|\Phi_{k,i}>
    ! when isolver == RMM2P,
    !  phi (:,1:kimg,2)    := |\Phi^{-1}_{k,i}> = |\Phi^{m}_{k,i}> - |\Phi^{m-1}_{k,i}>
    !  Rphi(:,1:kimg,2)    := R|\Phi^{-1}_{k,i}>
    !
    ! output
    !  zaj_l := |\Psi^{m+1}_{k,i}>
    !         = \alpha_0 phi(:,:,0) + \alpha_1 phi(:,:,1) + \alpha_2 phi(:,:,2)
    ! where,
    !  K is preconditioning operator,
    !  alpha_0,alpha_1, and alpha_2 are decided from <Rphi_i|Rphi_j> and <phi_i|S|phi_j>
    !
    ! For details, refer to G. Kresse and J. Furthmuller, PRB54,11169,(1996)
    ! The procedure in a case of isolver == RMM2P has been developped originally by T. Uda.
    !

! === Debug!!! These arrays should be saved for alloc_zfc_zfs call!!! by T.Kato ==========
!   integer :: n_max(3), n_min(3)
    integer, save :: n_max(3) = 0, n_min(3) = 0
! ========================================================================================
    integer :: ispin,ik,ib,iter_rmm, nrmm, matm,iilmt
    real(kind=DP), pointer, dimension(:)   :: ekin
    real(kind=DP), pointer, dimension(:,:) :: phi0
    real(kind=DP), allocatable,dimension(:,:,:) :: vxw_exx
    real(kind=DP), allocatable, dimension(:,:) :: afft_spin
    integer :: id_sname = -1
    real(kind=DP) :: exx
    integer :: ng
    integer :: iup
    logical :: store_e
    integer :: ig
    real(kind=DP) :: c1
    real(kind=DP), allocatable :: vtau_phl(:,:), cfft(:)

!!$    call tstatc0_begin('m_ESrmm_renew_WF ', id_sname,1)

    rr_avr = 0.d0
    call what_is_the_dimension_of_rmm() ! -(contained here) ->nrmm(= 2 or 3)

    if(potential_update>0) sw_retard_eigval_evaluation = OFF

! === DEBUG by tkato 2011/07/12 ================================================
!   allocate(rr_is_over_or_under(np_e,ista_k:iend_k))
    if(.not.allocated(rr_is_over_or_under)) then
       allocate(rr_is_over_or_under(np_e,ista_k:iend_k))
       rr_is_over_or_under = OVER
    endif
! ==============================================================================

    call m_ES_alloc_fft_related()  ! afft, bfft, m_FFT_alloc_WF_work
    if(isolver==RMM2P.and.nrmm==3) call m_ESsd_diff_WFs_to_zaj_old(1.d0) !->zaj_old

    if(rmm_save_memory_mode == ON) then
       call alloc_zfc_zfs(nfout,n_min,n_max,rmm_precal_phase_matm,natm,nrmm,matm) ! -(m_ES_WF_by_RMM)
       call zf_listing()                     ! -(contained here) -> zfc1,zfc2,zfc3,zfs1,zfs2,zfs3
    end if

    allocate(ekin0(kg1))
    ekin => ekin0(1:kg1)

    if ( use_metagga .and. vtau_exists ) allocate( cfft(nfft) )

    if(sw_fef==ON) call m_FEF_build_grad
    call m_ES_init_chgsoft()
    if(nrank_s==2) then
      allocate(afft_spin(nfft,2))
      do ispin = 1, nspin, (af+1)
        call m_ES_Vlocal_in_Rspace(ispin,afft_spin(:,ispin))
      enddo
    endif
!    Loop_spin: do ispin = 1, nspin, af+1
    Loop_spin: do ispin = ista_spin, iend_spin, af+1
       if(nrank_s==1) then
         call m_ES_Vlocal_in_Rspace(ispin,afft)  ! (ptfft1) ->afft
         if ( use_metagga .and. vtau_exists ) then
            call m_ES_Vlocal_in_Rspace( ispin, cfft, vtau_l )    ! r space
         endif
       else
         afft(:) = afft_spin(:,ispin)
       endif

       Loop_kpoints: do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) /= myrank_k) cycle              ! MPI

          if(iprirmm >= 3 .and. printable) call m_ES_wd_zaj_small_portion(nfout,ik," -- rmm --",10)

          call m_pwBS_kinetic_energies(ik,vkxyz,ekin)! (diakin) ->ekin

          if(rmm_save_memory_mode == ON) then
             call alloc_phasecs_bW_phi_Rphi(nrmm,n_min,n_max,matm) ! allocations of bWr, bWi, phi and Rphi. (bWr,bWi) = <\beta |WF=zaj_l>
             call pre_calc_phase !-(contained here) -> phasec,phases,newp

             Loop_band: do ib = ista_e, iend_e, istep_e      ! MPI
                if(rr_is_over_or_under(map_z(ib),ik) == UNDER.and.sw_precalculate==OFF) cycle
                Loop_res_dim: do iter_rmm = 0, nrmm-1
                   if(isolver == RMM2P .and. iter_rmm == nrmm-1) &
                        & call zajold2zaj_phi2zaj_old !-(contained here) zaj_old->zaj, phi0 ->zaj_old

                   call Vnonlocal_W_RMM(ik,ib,ispin,rmm_precal_phase_matm,iter_rmm)! ->vnlph_l,bWr,bWi
!                   if(sw_hybrid_functional==ON) call m_ES_Vexx_W_RMM(ik,ib)
                   if(sw_fef==ON) call m_FEF_add_grad_to_vnlph_RMM(ik,ib)
                   if ( use_metagga .and. vtau_exists ) then
                      call m_ES_kindens_to_vnlph_ib( ispin, ik, map_z(ib), cfft )
                   endif
                   if(iter_rmm>0 .or. sw_keep_hloc_phi==OFF) then
                     call m_ES_WF_in_Rspace(ik,ib,bfft) ! ->bfft(=|WF(R)>)
                     call m_FFT_Vlocal_W(afft,bfft) ! (afft,bfft)->bfft
                     call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON) ! bfft(R-space)->bfft(G-space)
                   endif
                   call rmm1(ib,iter_rmm,nrmm)     ! -(contained here)
                   !          zaj_l,vnlph_l,ekin,eko_l,bfft -> zaj_l,phi,Rphi
                end do Loop_res_dim
                call rmm_n_uda(ik,ib,nrmm)       ! -(contained here) phi,Rphi ->zaj_l, nrmm = 2 or 3
                !                       Residual norm is also checked (->rr_is_oever_or_under)
             end do Loop_band

             call dealloc_phasecs_bW_phi_Rphi()
          else
             call alloc_bW_phi_Rphi(nrmm) ! allocations of bWr, bWi, phi and Rphi. (bWr,bWi) = <\beta |WF=zaj_l>
             Loop_res_dim_b: do iter_rmm = 0, nrmm-1
                iup = 2
                !if(iter_rmm.gt.0) iup = 0
                if(isolver == RMM2P .and. iter_rmm == nrmm-1) &
                     & call zajold2zaj_phi2zaj_old_all !-(contained here) zaj_old->zaj, phi0 ->zaj_old, if rr_is_over_or_under /= UNDER

#ifdef RMM_NONLOCAL_NEW
                call m_ES_betar_dot_WFs_4_each_k(nfout,ik)
                do ib=1,np_e
                  do iilmt=1,nlmta
                    bWr(iilmt,iter_rmm,ib) = fsr_l(ib,iilmt,ik)
                    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) &
                    &   bWi(iilmt,iter_rmm,ib) = fsi_l(ib,iilmt,ik)
                  enddo
                enddo
                call m_ES_Vnonlocal_W(ik,(ik-1)/nspin+1,ispin,ON)
#else
                call Vnonlocal_W_RMMn(nfout,ik,ispin,iter_rmm)! ->vnlph_l,bWr,bWi
#endif
                if(sw_fef==ON) call m_FEF_add_grad_to_vnlph(ik)
                if(sw_hybrid_functional==ON) &
                &  call m_ES_Vexx_W(ik,2,store_exxp=iter_rmm==0)

                if ( use_metagga .and. vtau_exists ) then
                   call m_ES_contrib_kindens_to_vnlph( ispin, ik, cfft )
                endif

                do ib = ista_e, iend_e, istep_e
                   if(rr_is_over_or_under(map_z(ib),ik) == UNDER.and.sw_precalculate==OFF) cycle
                   if(iter_rmm>0 .or. sw_keep_hloc_phi==OFF)then
                     call m_ES_WF_in_Rspace(ik,ib,bfft) ! ->bfft(=|WF(R)>)
                     call m_FFT_Vlocal_W(afft,bfft) ! (afft,bfft)->bfft
                     call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON) ! bfft(R-space)->bfft(G-space)
                   endif
                   call rmm1(ib,iter_rmm,nrmm)     ! -(contained here)
                      !          zaj_l,vnlph_l,ekin,eko_l,bfft -> zaj_l,phi,Rphi
                end do
             end do Loop_res_dim_b

             Loop_band_b: do ib = ista_e, iend_e, istep_e
                !!$if(rr_is_over_or_under(map_z(ib),ik) == UNDER) cycle
                !! rewrote for ifort on em64t 2006/4/19 T.Yamamoto
                if(rr_is_over_or_under(map_z(ib),ik) /= UNDER) &
                   & call rmm_n_uda(ik,ib,nrmm)       ! -(contained here) phi,Rphi ->zaj_l, nrmm = 2 or 3
                !                       Residual norm is also checked (->rr_is_oever_or_under)
             end do Loop_band_b

             call dealloc_bW_phi_Rphi()
          end if

          call orthonorm_or_norm(imGSrmm)  !-(contained here)
          if(sw_hybrid_functional==ON) then
             if(sw_retard_eigval_evaluation==ON) then
                call m_ES_EXX_cp_eigenvalue(ik)
             else
                call m_ES_EXX_eigenvalue_for_each_k(ispin,ik,iupdate=0)
             endif
          endif
          call m_ES_eigen_values_for_each_k(ispin,ik,ekin,afft,eval_charge=sw_reduce_fft_for_charge==ON)

          if ( use_metagga .and. vtau_exists ) then
             allocate( vtau_phl(kg1,kimg) )
             Do ib=1, np_e
                call m_ES_kindens_to_vnlph_ib2( ispin, ik, ib, cfft, vtau_phl )
                c1 = 0.0d0
                Do ig=1, iba(ik)
                   c1 = c1 +vtau_phl(ig,1) *zaj_l(ig,ib,ik,1) &
                        &  +vtau_phl(ig,2) *zaj_l(ig,ib,ik,2)
                End Do
                eko_l(ib,ik) = eko_l(ib,ik) +c1
             ENd Do
             deallocate( vtau_phl )
          endif

       end do Loop_kpoints
    end do Loop_spin
    if(sw_hybrid_functional==ON.and.sw_retard_eigval_evaluation==ON)then
       call m_ES_EXX_gather_valence_states(nfout)
       do ispin=1,nspin,af+1
          do ik=ispin, kv3+ispin-nspin, nspin
             if(map_k(ik) /= myrank_k) cycle ! MPI
             call m_ES_EXX_eigenvalue_for_each_k(ispin,ik,update_eko=.false.,iupdate=2)
          enddo
       enddo
    endif
    if ( allocated(cfft) ) deallocate( cfft )

    if(rmm_save_memory_mode == ON) call dealloc_zfc_zfs()

    call rr_avr_final !-(c.h.)  rr_avr = rr_avr/(kv3*neg)

    deallocate(ekin0)

    call m_ES_dealloc_fft_related()
    if(nrank_s==2) deallocate(afft_spin)

!!$    call tstatc0_end(id_sname)
#ifdef fj_dump
  call dump_zaj(zaj_l)
#endif
  contains
    subroutine orthonorm_or_norm(imGSrmm_t)
      integer, intent(in) :: imGSrmm_t
      if(mod(iteration_electronic - iteration_rmm_start, imGSrmm_t) == 0) then
         if(sw_gep == ON)then
            if(iprirmm >= 2 .and. ik == 1) write(nfout,'(" -- NORMALIZATION --")')
            call m_ES_MGS_4_each_k(nfout,ik,mode=NORMALIZATION)
         else
            if(iprirmm >= 2 .and. ik == 1) write(nfout,'(" -- ORTHONORMALIZATION --")')
            call m_ES_MGS_4_each_k(nfout,ik,mode=ORTHONORMALIZATION)
         endif
      else
         if(iprirmm >= 2 .and. ik == 1) write(nfout,'(" -- NORMALIZATION --")')
         call m_ES_MGS_4_each_k(nfout,ik,mode=NORMALIZATION)
      end if
#ifdef SAVE_FFT_TIMES
      if(sw_save_fft == ON) status_saved_phifftr(:,ik) = OLD
#endif
    end subroutine orthonorm_or_norm
!!$!BRANCH_P_END ORG_Parallel

    subroutine zajold2zaj_phi2zaj_old_all
      integer :: ib, ibto
      do ib = ista_e, iend_e, istep_e
         if(rr_is_over_or_under(map_z(ib),ik) == UNDER.and.sw_precalculate==OFF) cycle
         if(nrmm == 3) call m_ESsd_copy_zaj_old_to_zaj(ik,ib)
         ibto = map_z(ib)
         phi0 => phi(:,:,ibto,0)
         call m_ESsd_copy_phi_to_zaj_old(ik,ib,phi0)  ! phi->zaj_old
      end do
    end subroutine zajold2zaj_phi2zaj_old_all

    subroutine zajold2zaj_phi2zaj_old
      if(nrmm == 3) call m_ESsd_copy_zaj_old_to_zaj(ik,ib)
      phi0 => phi(:,:,1,0)
      call m_ESsd_copy_phi_to_zaj_old(ik,ib,phi0)  ! phi->zaj_old
    end subroutine zajold2zaj_phi2zaj_old

    subroutine rr_avr_final
      integer       :: i,j, n_rr_under, n_mpi
      real(kind=DP) :: rr_avr_mpi

      if(npes >= 2) then
         call mpi_allreduce(rr_avr, rr_avr_mpi,1,mpi_double_precision &
              & , mpi_sum,MPI_CommGroup,ierr)
      else
         rr_avr_mpi = rr_avr
      end if
      rr_avr = rr_avr_mpi/(kv3*neg)
      if(iprirmm >= 2) &
           & write(nfout,'(" Residual Norm (sum_i <R_i|R_i>/sum_i) ( " &
           & ,i6," -th iter (electronic)) = ",d19.7)') iteration_electronic, rr_avr

      n_rr_under = 0
      do j = ista_k, iend_k
         do i = 1, np_e
            if(rr_is_over_or_under(map_z(i),j) == UNDER) n_rr_under = n_rr_under + 1
         end do
      end do
      if(npes >= 2) then
         call mpi_allreduce(n_rr_under,n_mpi, 1, mpi_integer,mpi_sum,MPI_CommGroup,ierr)
      else
         n_mpi = n_rr_under
      end if
      if(n_mpi > 0 .and. iprirmm >= 1) write(nfout,'(" Number of rr_under = ",i5)') n_mpi
    end subroutine rr_avr_final

    subroutine what_is_the_dimension_of_rmm
      if(isolver == RMM2 ) then
         nrmm = 2
      else if(isolver == RMM2p) then
         if(iteration_electronic == iteration_rmm_start) then
            nrmm = 2
         else
            nrmm = 3
         end if
      else if(isolver == RMM3) then
         nrmm = 3
      else
         if(printable) write(nfout,'(" !isolver",i3," is illegal ( m_ESrmm_renew_WF )")') isolver
         call phase_error_with_msg(nfout,'solver is illegal ( m_ESrmm_renew_WF )',__LINE__,__FILE__)
      endif
    end subroutine what_is_the_dimension_of_rmm


    subroutine rmm1(ibo,iter_rmm,nrmm)
      ! Revised by T. Yamasaki, 18th Sep. 2004
      !    <R_i|R_j> and <Phi_i|S|Phi_j> are calculated, then are substituted to
      !   rr_e and psp_e, respectively with not using
      integer, intent(in) :: ibo,iter_rmm,nrmm  ! nrmm = {2|3}
      !               iter_rmm = {0|1} (when nrmm == 2)}, {0|1|2} (when nrmm == 3)
      real(kind=DP) :: evr,devr,dnm,evi,e1,devi
      real(kind=DP) :: rrr, rri, pspr, pspi

      real(kind=DP), pointer, dimension(:,:) :: phi_t, Rphi_t !d(iba(ik),kimg)
      real(kind=DP), allocatable, dimension(:) :: p
      integer       :: i, i1, ib, ibt, ir, ii, j

      allocate(p(kg1))

      if(iter_rmm == nrmm-1) then
         allocate(phi_t(iba(ik),kimg))
         allocate(Rphi_t(iba(ik),kimg))
      end if

      ib = map_z(ibo)                  ! MPI
      if(rmm_save_memory_mode == ON) then
         ibt = 1
      else
         ibt = ib
      end if
      dnm = 1.d0/product(fft_box_size_WF(1:3,1))
      call m_ES_decide_precon_factor(precon,ik,ibo,ekin,p) ! ->p(1:iba(ik))

      if(kimg == 1) then
         if(iter_rmm <= nrmm-2) then
            do i = 1, iba(ik)
               i1    = igf(nbase(i,ik))
               evr   = zaj_l(i,ib,ik,1)
               if(iter_rmm==0 .and. sw_keep_hloc_phi==ON) then
                 devr  = hlocphi_l(i,ib,ik,1)-eko_l(ib,ik)*evr+vnlph_l(i,ib,1)
               else
                 devr  = (ekin(i)-eko_l(ib,ik))*evr + bfft(i1)*dnm+vnlph_l(i,ib,1)
               endif
               zaj_l(i,ib,ik,1)    = -p(i)*devr
               phi(i,1,ibt,iter_rmm)   =  evr
               Rphi(i,1,ibt,iter_rmm)  = -devr
            end do
#ifdef SAVE_FFT_TIMES
            if(sw_save_fft == ON) status_saved_phifftr(ib,ik) = OLD
#endif
            if(k_symmetry(ik) == GAMMA) then
               do i1 = 0, iter_rmm
                  rr_e(i1,iter_rmm,ibt,1) = Rphi(1,1,ibt,i1)*Rphi(1,1,ibt,iter_rmm) + &
                       & 2.d0*dot_product(Rphi(2:iba(ik),1,ibt,i1),Rphi(2:iba(ik),1,ibt,iter_rmm))
                  psp_e(i1,iter_rmm,ibt,1) = phi(1,1,ibt,i1)*phi(1,1,ibt,iter_rmm) + &
                       & 2.d0*dot_product(phi(2:iba(ik),1,ibt,i1),phi(2:iba(ik),1,ibt,iter_rmm))
               end do
            else
               do i1 = 0, iter_rmm
                  rr_e(i1,iter_rmm,ibt,1) = dot_product(Rphi(1:iba(ik),1,ibt,i1),Rphi(1:iba(ik),1,ibt,iter_rmm))
                  psp_e(i1,iter_rmm,ibt,1) = dot_product(phi(1:iba(ik),1,ibt,i1),phi(1:iba(ik),1,ibt,iter_rmm))
               end do
            end if
         else if(iter_rmm == nrmm-1) then
            do i = 1, iba(ik)
               i1    = igf(nbase(i,ik))
               phi_t(i,1)   = zaj_l(i,ib,ik,1)
               Rphi_t(i,1)  = - ((ekin(i)-eko_l(ib,ik))*phi_t(i,1) + bfft(i1)*dnm+vnlph_l(i,ib,1))
!!$               zaj_l(i,ib,ik,1) = evr
               !   when iter_rmm == nrmm-1, zaj_l (to which KR|Psi> is stored) is not
               !   going to be used recursively, so we do not need to store evr to zaj_l,
               !   instead we ues zaj_l as |phi_{iter_rmm}>.
            end do
            if(k_symmetry(ik) == GAMMA) then
               do i1 = 0, iter_rmm-1
                  rr_e(i1,iter_rmm,ibt,1) = Rphi(1,1,ibt,i1)*Rphi_t(1,1) &
                       & + 2.d0*dot_product(Rphi(2:iba(ik),1,ibt,i1),Rphi_t(2:iba(ik),1))
                  psp_e(i1,iter_rmm,ibt,1) = phi(1,1,ibt,i1)*phi_t(1,1) &
                       & + 2.d0*dot_product(phi(2:iba(ik),1,ibt,i1),phi_t(2:iba(ik),1))
               end do
               rr_e(iter_rmm,iter_rmm,ibt,1) = Rphi_t(1,1)*Rphi_t(1,1) &
                    & + 2.d0*dot_product(Rphi_t(2:iba(ik),1),Rphi_t(2:iba(ik),1))
               psp_e(iter_rmm,iter_rmm,ibt,1) = phi_t(1,1)*phi_t(1,1) &
                    & + 2.d0*dot_product(phi_t(2:iba(ik),1),phi_t(2:iba(ik),1))
            else
               do i1 = 0, iter_rmm-1
                  rr_e(i1,iter_rmm,ibt,1) = dot_product(Rphi(1:iba(ik),1,ibt,i1),Rphi_t(1:iba(ik),1))
                  psp_e(i1,iter_rmm,ibt,1) = dot_product(phi(1:iba(ik),1,ibt,i1),phi_t(1:iba(ik),1))
               end do
               rr_e(iter_rmm,iter_rmm,ibt,1) = dot_product(Rphi_t(1:iba(ik),1),Rphi_t(1:iba(ik),1))
               psp_e(iter_rmm,iter_rmm,ibt,1) = dot_product(phi_t(1:iba(ik),1),phi_t(1:iba(ik),1))
            end if
         else
            if(printable) write(nfout,'("  iter_rmm > nrmm-1 ( illegal relation) <<m_ES_WF_by_RMM.rmm1>>")')
            call phase_error_with_msg(nfout, '   iter_rmm > nrmm-1 ( illegal relation) <<m_ES_WF_by_RMM.rmm1>>', &
            __LINE__,__FILE__)
         end if
      else if(kimg == 2) then
         ir = iter_rmm
         if(iter_rmm <= nrmm-2) then
            do i = 1, iba(ik)
               i1    = igf(nbase(i,ik))
               evr   = zaj_l(i,ib,ik,1);    evi   = zaj_l(i,ib,ik,2)
               if(iter_rmm==0 .and. sw_keep_hloc_phi==ON) then
                 e1    = - eko_l(ib,ik)
                 devr  = hlocphi_l(i,ib,ik,1)+e1*evr+vnlph_l(i,ib,1)
                 devi  = hlocphi_l(i,ib,ik,2)+e1*evi+vnlph_l(i,ib,2)
               else
                 e1    = ekin(i) - eko_l(ib,ik)
                 devr  = e1*evr+bfft(2*i1-1)*dnm+vnlph_l(i,ib,1)
                 devi  = e1*evi+bfft(2*i1  )*dnm+vnlph_l(i,ib,2)
               endif
               zaj_l(i,ib,ik,1)   = -p(i)*devr
               zaj_l(i,ib,ik,2)   = -p(i)*devi
               phi(i,   1,ibt,iter_rmm)  = evr;    phi(i, 2,ibt,iter_rmm) = evi
               Rphi(i,  1,ibt,iter_rmm)  = -devr; Rphi(i, 2,ibt,iter_rmm) = -devi
            end do
#ifdef SAVE_FFT_TIMES
            if(sw_save_fft == ON) status_saved_phifftr(ib,ik) = OLD
#endif
            do i1 = 0, iter_rmm
               rrr = 0.d0; rri = 0.d0; pspr = 0.d0; pspi = 0.d0
               if(k_symmetry(ik) == GAMMA) then
                  do i = 2, iba(ik)
                     rrr = rrr + Rphi(i,1,ibt,i1)*Rphi(i,1,ibt,iter_rmm) + Rphi(i,2,ibt,i1)*Rphi(i,2,ibt,iter_rmm)
                     pspr = pspr + phi(i,1,ibt,i1)*phi(i,1,ibt,iter_rmm) + phi(i,2,ibt,i1)*phi(i,2,ibt,iter_rmm)
                  end do
                  rrr = 2.d0*rrr + Rphi(1,1,ibt,i1)*Rphi(1,1,ibt,iter_rmm)
                  pspr = 2.d0*pspr + phi(1,1,ibt,i1)*phi(1,1,ibt,iter_rmm)
               else
                  do i = 1, iba(ik)
                     rrr = rrr + Rphi(i,1,ibt,i1)*Rphi(i,1,ibt,iter_rmm) + Rphi(i,2,ibt,i1)*Rphi(i,2,ibt,iter_rmm)
                     rri = rri + Rphi(i,1,ibt,i1)*Rphi(i,2,ibt,iter_rmm) - Rphi(i,2,ibt,i1)*Rphi(i,1,ibt,iter_rmm)
                     pspr = pspr + phi(i,1,ibt,i1)*phi(i,1,ibt,iter_rmm) + phi(i,2,ibt,i1)*phi(i,2,ibt,iter_rmm)
                     pspi = pspi + phi(i,1,ibt,i1)*phi(i,2,ibt,iter_rmm) - phi(i,2,ibt,i1)*phi(i,1,ibt,iter_rmm)
                  end do
               end if
               rr_e(i1,iter_rmm,ibt,1) = rrr;   rr_e(i1,iter_rmm,ibt,2) = rri
               psp_e(i1,iter_rmm,ibt,1) = pspr; psp_e(i1,iter_rmm,ibt,2) = pspi
            end do
         else if(iter_rmm == nrmm-1) then
            do i = 1, iba(ik)
               i1    = igf(nbase(i,ik))
               phi_t(i,1) = zaj_l(i,ib,ik,1);    phi_t(i,2) = zaj_l(i,ib,ik,2)
               e1    = ekin(i) - eko_l(ib,ik)
               Rphi_t(i,1) = - (e1*phi_t(i,1)+bfft(2*i1-1)*dnm+vnlph_l(i,ib,1))
               Rphi_t(i,2) = - (e1*phi_t(i,2)+bfft(2*i1  )*dnm+vnlph_l(i,ib,2))
            end do
            do i1 = 0, iter_rmm-1
               rrr = 0.d0; rri = 0.d0; pspr = 0.d0; pspi = 0.d0
               if(k_symmetry(ik) == GAMMA) then
                  do i = 2, iba(ik)
                     rrr = rrr + Rphi(i,1,ibt,i1)*Rphi_t(i,1) + Rphi(i,2,ibt,i1)*Rphi_t(i,2)
                     pspr = pspr + phi(i,1,ibt,i1)*phi_t(i,1) + phi(i,2,ibt,i1)*phi_t(i,2)
                  end do
                  rrr = 2.d0*rrr + Rphi(1,1,ibt,i1)*Rphi_t(1,1)
                  pspr = 2.d0*pspr + phi(1,1,ibt,i1)*phi_t(1,1)
               else
                  do i = 1, iba(ik)
                     rrr = rrr + Rphi(i,1,ibt,i1)*Rphi_t(i,1) + Rphi(i,2,ibt,i1)*Rphi_t(i,2)
                     rri = rri + Rphi(i,1,ibt,i1)*Rphi_t(i,2) - Rphi(i,2,ibt,i1)*Rphi_t(i,1)
                     pspr = pspr + phi(i,1,ibt,i1)*phi_t(i,1) + phi(i,2,ibt,i1)*phi_t(i,2)
                     pspi = pspi + phi(i,1,ibt,i1)*phi_t(i,2) - phi(i,2,ibt,i1)*phi_t(i,1)
                  end do
               end if
               rr_e(i1,iter_rmm,ibt,1) = rrr; rr_e(i1,iter_rmm,ibt,2) = rri
               psp_e(i1,iter_rmm,ibt,1) = pspr; psp_e(i1,iter_rmm,ibt,2) = pspi
            end do
            rrr = 0.d0; rri = 0.d0; pspr = 0.d0; pspi = 0.d0
            if(k_symmetry(ik) == GAMMA) then
               do i = 2, iba(ik)
                  rrr = rrr + Rphi_t(i,1)*Rphi_t(i,1) + Rphi_t(i,2)*Rphi_t(i,2)
                  pspr = pspr + phi_t(i,1)*phi_t(i,1) + phi_t(i,2)*phi_t(i,2)
               end do
               rrr = 2.d0*rrr + Rphi_t(1,1)*Rphi_t(1,1)
               pspr = 2.d0*pspr + phi_t(1,1)*phi_t(1,1)
            else
               do i = 1, iba(ik)
                  rrr = rrr + Rphi_t(i,1)*Rphi_t(i,1) + Rphi_t(i,2)*Rphi_t(i,2)
                  rri = rri + Rphi_t(i,1)*Rphi_t(i,2) - Rphi_t(i,2)*Rphi_t(i,1)
                  pspr = pspr + phi_t(i,1)*phi_t(i,1) + phi_t(i,2)*phi_t(i,2)
                  pspi = pspi + phi_t(i,1)*phi_t(i,2) - phi_t(i,2)*phi_t(i,1)
               end do
            end if
            rr_e(iter_rmm,iter_rmm,ibt,1) = rrr; rr_e(iter_rmm,iter_rmm,ibt,2) = rri
            psp_e(iter_rmm,iter_rmm,ibt,1) = pspr; psp_e(iter_rmm,iter_rmm,ibt,2) = pspi
         else
            if(printable) write(nfout,'("  iter_rmm > nrmm-1 ( illegal relation) <<m_ES_WF_by_RMM.rmm1>>")')
            call phase_error_with_msg(nfout, '   iter_rmm > nrmm-1 ( illegal relation) <<m_ES_WF_by_RMM.rmm1>>'&
                                     ,__LINE__,__FILE__)
         end if
      end if

      if(iter_rmm == nrmm-1) then
         deallocate(Rphi_t)
         deallocate(phi_t)
      end if
      deallocate(p)
    end subroutine rmm1

    recursive subroutine rmm_n(ibo,nrmm,dtim)
!
!    * Rewritten by T. Uda and T. Yamasaki, 19th Mar. 2003
!      + alpha(nrmm-1) -> alpha(nrmm-1,kimg)
!
      integer,       intent(in) :: ibo,nrmm
      real(kind=DP), intent(in) :: dtim

      integer                   :: ip(nrmm)
      real(DP) :: ymat(0:nrmm-1,0:nrmm-1,kimg)   &
           &     ,vector(0:nrmm-1,0:nrmm-1,kimg) &
           &     ,rr(0:nrmm-1,0:nrmm-1,kimg),psp(0:nrmm-1,0:nrmm-1,kimg)
      real(DP) :: eigenr(0:nrmm-1), eigeni(0:nrmm-1)
      real(DP) :: alpha(nrmm-1,kimg),ww1(nrmm*kimg),ww2(nrmm*kimg),ww3(nrmm*kimg)
      integer  :: ierrdecomp, ib

      ib = map_z(ibo)                 ! MPI
      call rr_and_psp_from_Rphi_and_phi(ib,nrmm,rr,psp) !-(m_ES_WF_by_RMM)
      !   rr(i,j) = <R_i|R_j>; psp(i,j) = <Phi_i|S|Phi_j>
      !                            S = 1 + \sum_{ij}q_{ij}|b_i><b_j|
      call normalization_of_rr_and_psp(nrmm,rr,psp)
      !   rr(i,j) = rr(i,j)/rr(0,0), psp(i,j) = psp(i,j)/psp(0,0)

!* Solve the Eigen value problem of Bx = epsilon*Ax
!* Here B is the residual matrix: B_ij = <R^i|R^j>,
!* and A is <phi^i|S|phi^j>.
!* It is equivallent to eigen value problem of Yx = epsilon*x
!* where AY = B.
      call LUdecomp(nrmm,psp,ww1,ip,ierrdecomp)
      if(ierrdecomp /= 0) then
         if(printable) write(nfout,*) 'LU decomposition is impossible.';  alpha = 0.0d0
         alpha(1,1) = dtim
      else
         call solveAYB(nrmm,psp,rr,ymat,ip)      ! Solve AY = B
         call rg_or_cg(nrmm,ymat,eigenr,eigeni,vector,ww1,ww2,ww3)
         call decide_alpha(kimg,dtim,nrmm,eigenr,eigeni,vector,alpha)
         !                    -(rmmsubs) ->alpha
      end if

      if(nrmm == 3 .and. ierrdecomp /= 0) then
         call rmm_n(ibo,2,dtim)
      else
         if(iprirmm >= 2) call wd_alpha(nrmm,alpha)
         call evolve_WF_using_Residuals(ib,nrmm,alpha) ! -(contained subr. m_ESrmm_renew_WF)
         if(iprirmm>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')
         !    phi,alpha ->zaj_l(:,ib,ik,:)
         !   |WF> = |phi_0> + alpha(1)|phi_1> + alpha(2)|phi_2> , when nrmm == 3
      end if
    end subroutine rmm_n

    subroutine evolve_WF_using_Residuals(ib,nrmm,alpha)
!
!  * Rewritten by T. Uda and T. Yamasaki, 19th Mar. 2003
!      + alpha(nrmm-1) -> alpha(nrmm-1,kimg)
!
!  * Revised by T. Yamasaki, 18th Sep. 2004
!      + being changed not to use phi(:,:,:,nrmm-1) for saving the memory allocation area.
!      + phi(:,:,:,nrmm-1) is substituted to zaj_l before the program reached to this
!       subroutine.
!
      integer, intent(in) ::                             ib, nrmm
      real(kind=DP),intent(in),dimension(nrmm-1,kimg) :: alpha
      integer ::                                         i, ibt, nr, ir, ii, j, ig
      real(kind=DP) ::                                   zr,zi

      if(rmm_save_memory_mode == ON) then
         ibt = 1
      else
         ibt = ib
      end if

      if(kimg==1) then
         ! previously in sub. rmm1, zaj_l(:,ib,ik,:) <= phi(:,:,ibt,nrmm-1)
         ! Following calculation is equivalent to
         !  zaj_l = phi(:,:,ibt,0) + alpha(nrmm-1)*phi(:,:,ibt,nrmm-1)
         !    |WF> = |phi_0> + alpha(nrmm-1)|phi_{nrmm-1}>
         do i = 1, iba(ik)
            zaj_l(i,ib,ik,1) = phi(i,1,ibt,0) + alpha(nrmm-1,1)*zaj_l(i,ib,ik,1)
         end do
         do ir = 1, nrmm-2
            do i = 1, iba(ik)
               zaj_l(i,ib,ik,1) = zaj_l(i,ib,ik,1) + alpha(ir,1)*phi(i,1,ibt,ir)
            end do
         end do
      else if(kimg==2) then
         ! previously in <sub. rmm1>, zaj_l(:,ib,ik,:) <-- phi(:,:,ibt,nrmm-1)
         nr = nrmm-1
         do i = 1, iba(ik)
            !   Calculation in this do-loop is equivalent to
            !  zaj_l(:,ib,ik,:) = phi(:,:,ibt,0) + alpha(nrmm-1)*phi(:,:,ibt,nrmm-1),
            !  namely, |WF> = |phi_0> + alpha(nrmm-1)|phi_{nrmm-1}>
            zr = zaj_l(i,ib,ik,1)
            zi = zaj_l(i,ib,ik,2)
            zaj_l(i,ib,ik,1) = phi(i,1,ibt,0) + alpha(nr,1)*zr - alpha(nr,2)*zi
            zaj_l(i,ib,ik,2) = phi(i,2,ibt,0) + alpha(nr,1)*zi + alpha(nr,2)*zr
         end do

         do i = 1, nrmm-2
            do ig = 1, iba(ik)
               zaj_l(ig,ib,ik,1) = zaj_l(ig,ib,ik,1) + alpha(i,1)*phi(ig,1,ibt,i) - alpha(i,2)*phi(ig,2,ibt,i)
               zaj_l(ig,ib,ik,2) = zaj_l(ig,ib,ik,2) + alpha(i,1)*phi(ig,2,ibt,i) + alpha(i,2)*phi(ig,1,ibt,i)
            end do
         end do
      else
         call phase_error_with_msg(nfout, ' kimg is illegal (evolve_WF_using_Residuals)',__LINE__,__FILE__)
      end if
#ifdef SAVE_FFT_TIMES
      if(sw_save_fft == ON) status_saved_phifftr(ib,ik) = OLD
#endif

      if(iprirmm >= 3) then
         ibt = 100
         if(ibt > iba(ik)) ibt = iba(ik)
         if(ik <= 2 .or. ik >= kv3-1) then
            write(nfout,'(" !rmm  zaj_l in <<evolve_WF_using_Residuals>>")')
            write(nfout,'(" !rmm  ik, ib = ",2i8)') ik, ib
            if(kimg == 2) write(nfout,'(" !rmm -- real part --")')
            write(nfout,'(" !rmm  ",8f8.4)') (zaj_l(i,ib,ik,1),i=1,ibt)
            if(kimg == 2) write(nfout,'(" !rmm -- imaginary part --")')
            write(nfout,'(" !rmm  ",8f8.4)') (zaj_l(i,ib,ik,2),i=1,ibt)
         end if
      end if

    end subroutine evolve_WF_using_Residuals

    subroutine wd_rr_and_psp(nrmm,rr,psp)
      integer, intent(in)                                        :: nrmm
      real(kind=DP),intent(in),dimension(0:nrmm-1,0:nrmm-1,kimg) :: rr,psp
      integer :: i, j
      write(nfout,'("!rr =",6d15.7)') ((rr(i,j,1),j=i,nrmm-1),i=0,nrmm-1)
      write(nfout,'("!psp=",6d15.7)') ((psp(i,j,1),j=i,nrmm-1),i=0,nrmm-1)
    end subroutine wd_rr_and_psp

    subroutine wd_alpha(nrmm,alpha)
      integer, intent(in)                          :: nrmm
      real(kind=DP), intent(in), dimension(nrmm-1,kimg) :: alpha

      write(nfout,'(" ! alpha(1) = ", f20.10)') alpha(1,1)
      if(nrmm == 3) write(nfout,'(" ! alpha(2) = ", f20.10)') alpha(2,1)
    end subroutine wd_alpha

    subroutine rmm_n_uda(ik,ibo,nrmm)
      integer,       intent(in) :: ik,ibo, nrmm

      real(DP) :: rr(0:nrmm-1,0:nrmm-1,kimg),psp(0:nrmm-1,0:nrmm-1,kimg)
      real(DP) :: alpha(nrmm-1,kimg)
      integer  :: ib
      integer :: id_sname = -1
      call tstatc0_begin('rmm_n_uda ', id_sname)

      ib = map_z(ibo)
!!$      call rr_and_psp_from_Rphi_and_phi(ib,nrmm,rr,psp) !-(m_ES_WF_by_RMM)
      call rr_and_psp_from_rr_e_and_psp_e(ib,nrmm,rr,psp) ! -(m_ES_WF_by_RMM)
      !   rr(i,j) = <R_i|R_j>; psp(i,j) = <Phi_i|S|Phi_j>
      !                            S = 1 + \sum_{ij}q_{ij}|b_i><b_j|
      rr_avr = rr_avr + rr(0,0,1)
      if(rr(0,0,1) <= rr_Critical_Value) rr_is_over_or_under(ib,ik) = UNDER

    !!$  call normalization_of_rr_and_psp(nrmm,rr,psp)
      !   rr(i,j) = rr(i,j)/rr(0,0), psp(i,j) = psp(i,j)/psp(0,0)

      if(iprirmm >= 2) call wd_rr_and_psp(nrmm,rr,psp)

      if(kimg == 1) then
         if(nrmm == 2) then
            call rmm2_uda(iprirmm,nrmm,rr,psp,alpha)
         else
            call rmm3_uda(iprirmm,nrmm,rr,psp,alpha)
         end if
      else
         if(nrmm == 2) then
            call crmm2_uda(iprirmm,nrmm,rr,psp,alpha)
         else
            call crmm3_uda(iprirmm,nrmm,rr,psp,alpha)
         end if
      end if

      if(iprirmm >= 2) call wd_alpha(nrmm,alpha)
      call evolve_WF_using_Residuals(ib,nrmm,alpha) ! -(contained subr. m_ESrmm_renew_WF)
      if(iprirmm>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')
      !   phi,alpha ->zaj_l(:,ib,ik,:)
      !   |WF> = |phi_0> + alpha(1)|phi_1> + alpha(2)|phi_2> , when nrmm == 3
      call tstatc0_end(id_sname)
    end subroutine rmm_n_uda

    subroutine rg_or_cg(nrmm,ymat,eigenr,eigeni,vector,ww1,ww2,ww3)
      integer, intent(in) :: nrmm
      real(DP), dimension(0:nrmm-1,0:nrmm-1,kimg):: ymat
      real(DP), intent(out),dimension(0:nrmm-1)  :: eigenr, eigeni
      real(DP), intent(out),dimension(0:nrmm-1,0:nrmm-1,kimg):: vector
      real(DP), intent(out),dimension(nrmm*kimg) :: ww1,ww2,ww3
      integer                                    :: ier
      if(kimg == 1) then
         call rg_eispack(nrmm,nrmm,ymat,eigenr,eigeni,1,vector,ww1,ww2,ier)
      else if(kimg == 2) then
         call cg_eispack(nrmm,nrmm,ymat(0,0,1),ymat(0,0,kimg),eigenr,eigeni,1 &
              &, vector(0,0,1),vector(0,0,kimg),ww1,ww2,ww3,ier)
      end if
      if(ier /= 0 .and. printable) write(nfout,*) 'ier(rg_or_cg) = ',ier
    end subroutine rg_or_cg

    subroutine solveAYB(nrmm,psp,rr,ymat,ip)
      integer, intent(in) :: nrmm
      real(DP),intent(in), dimension(0:nrmm-1,0:nrmm-1,kimg):: psp,rr
      real(DP),intent(out),dimension(0:nrmm-1,0:nrmm-1,kimg):: ymat
      integer, intent(out),dimension(nrmm) :: ip(nrmm)

      if(kimg == 1) then
         call rsolve(nrmm,nrmm,psp,rr,ymat,ip)
      else
         call csolve2(nrmm,nrmm,psp,rr,ymat,ip)
      end if
    end subroutine solveAYB

    subroutine LUdecomp(nrmm,psp,ww1,ip,ier)
      integer, intent(in) :: nrmm
      real(DP),intent(inout), dimension(0:nrmm-1,0:nrmm-1,kimg):: psp
      real(DP),intent(out),dimension(nrmm*kimg)                :: ww1
      integer, intent(out)                 :: ip(nrmm), ier

      complex(CMPLDP)                      :: cpsp(nrmm*nrmm)
      integer i, j, n
      if(kimg == 1) then
         call rdecomp(nrmm,psp,ww1,ip,ier)
      else
         n = 0
         do j = 0, nrmm-1
            do i = 0, nrmm-1
               n = n + 1
               cpsp(n) = cmplx(psp(i,j,1),psp(i,j,2))
            end do
         end do
         call cdecomp(nrmm,cpsp,ww1,ip,ier)
      end if
    end subroutine LUdecomp

    subroutine rr_and_psp_from_Rphi_and_phi(ib_t,nrmm,rr,psp)
      integer, intent(in)       :: ib_t,nrmm
      real(DP),intent(out),dimension(0:nrmm-1,0:nrmm-1,kimg)::rr,psp
      real(DP), dimension(kimg) :: sumqff
      real(DP) ::                  rrr,rri,pspr,pspi
      integer  ::                  j,i,m,ri, ib

      if(rmm_save_memory_mode == ON) then
         ib = 1
      else
         ib = ib_t
      end if

      rr = 0.d0; psp = 0.d0
      do j = 0, nrmm-1
         do i = 0,j
            call phsphn(i,j,ib,kimg,sumqff)  ! ->sumqff (= \sum_{n1,n2}q_{n1,n2}<W|b_{nl}><b_{n2}|W>
            if(kimg == 1) then
               rr(i,j,1)  = dot_product(Rphi(1:iba(ik),1,ib,i),Rphi(1:iba(ik),1,ib,j))
               psp(i,j,1) = dot_product(phi(1:iba(ik),1,ib,i),phi(1:iba(ik),1,ib,j))
            else if(kimg == 2) then
               rrr = 0.d0;   rri = 0.d0
               pspr = 0.d0;  pspi = 0.d0
               do m = 1,iba(ik)
                  rrr = rrr + Rphi(m,1,ib,i)*Rphi(m,1,ib,j) + Rphi(m,2,ib,i)*Rphi(m,2,ib,j)
                  rri = rri + Rphi(m,1,ib,i)*Rphi(m,2,ib,j) - Rphi(m,2,ib,i)*Rphi(m,1,ib,j)
                  pspr = pspr + phi(m,1,ib,i)* phi(m,1,ib,j) + phi(m,2,ib,i)* phi(m,2,ib,j)
                  pspi = pspi + phi(m,1,ib,i)* phi(m,2,ib,j) - phi(m,2,ib,i)* phi(m,1,ib,j)
               end do
               rr(i,j,1) = rrr;               rr(i,j,2) = rri
               psp(i,j,1) = pspr;             psp(i,j,2) = pspi
            end if
            do ri = 1, kimg
               rr( j,i,ri) = pm(ri) * rr(i,j,ri)      ! pm=1,-1 (ri==1,2)
               psp(i,j,ri) = psp(i,j,ri) + sumqff(ri)
               psp(j,i,ri) = pm(ri) * psp(i,j,ri)
            end do
         end do
      end do

    end subroutine rr_and_psp_from_Rphi_and_phi

    subroutine rr_and_psp_from_rr_e_and_psp_e(ib_t,nrmm,rr,psp)
      integer, intent(in)       :: ib_t,nrmm
      real(DP),intent(out),dimension(0:nrmm-1,0:nrmm-1,kimg)::rr,psp
      real(DP), dimension(kimg) :: sumqff
      integer  ::                  j,i,m,ri, ib

      if(rmm_save_memory_mode == ON) then
         ib = 1
      else
         ib = ib_t
      end if

      rr = 0.d0; psp = 0.d0
      do j = 0, nrmm-1
         do i = 0,j
            call phsphn(i,j,ib,kimg,sumqff)  ! ->sumqff (= \sum_{n1,n2}q_{n1,n2}<W|b_{nl}><b_{n2}|W>
            do ri = 1, kimg
               rr(i,j,ri)  = rr_e(i,j,ib,ri)
               psp(i,j,ri) = psp_e(i,j,ib,ri)
               rr( j,i,ri) = pm(ri) * rr(i,j,ri)      ! pm=1,-1 (ri==1,2)
               psp(i,j,ri) = psp(i,j,ri) + sumqff(ri)
               psp(j,i,ri) = pm(ri) * psp(i,j,ri)
            end do
         end do
      end do

    end subroutine rr_and_psp_from_rr_e_and_psp_e

    subroutine phsphn(i,j,ib,nimg,sumqff)
      integer, intent(in)        :: i, j, ib, nimg  ! ib = 1 when save_memory_mode == ON
      real(kind=DP), intent(out) :: sumqff(nimg)
      integer         :: ia,m,n

      sumqff = 0.d0
      if(modnrm == EXECUT) then
         do ia = 1, nac
            m = nlmta1(ia); n = nlmta2(ia)
            sumqff(1) = sumqff(1) + fqwei(ia)*(bWr(m,i,ib)*bWr(n,j,ib)+bWi(m,i,ib)*bWi(n,j,ib))
         end do
         if(nimg == 2) then
            do ia = 1, nac
               m = nlmta1(ia); n = nlmta2(ia)
               sumqff(2) = sumqff(2) + fqwei(ia)*(bWr(m,i,ib)*bWi(n,j,ib)-bWi(m,i,ib)*bWr(n,j,ib))
            end do
         end if
      end if
    end subroutine phsphn

    subroutine normalization_of_rr_and_psp(nrmm,rr,psp)
      integer, intent(in) :: nrmm
      real(DP),intent(out),dimension(0:nrmm-1,0:nrmm-1,kimg)::rr,psp
      real(DP) :: divrr,divpsp

      divrr = 1.d0/rr(0,0,1); divpsp = 1.d0/psp(0,0,1)
      rr    =  rr*divrr
      psp   = psp*divpsp

      if(iprirmm >= 2) call wd_rr_and_psp(nrmm,rr,psp)
    end subroutine normalization_of_rr_and_psp

    subroutine pre_calc_phase
      call mltpha4(natm,rmm_precal_phase_matm,pos,PAI2,kgp,ngabc,kg1,nbase(1,ik),iba(ik) &
           &     ,phasec,phases)
      if(rmm_precal_phase_matm < natm) then
         call setglist4(n_min(1),n_max(1),n_min(2),n_max(2) &  ! -(b_E.S.)
              &, n_min(3),n_max(3),nbase(1,ik),ngabc,kgp,iba(ik),nglist)
         call crngabc4(n_min(1),n_max(1),n_min(2),n_max(2) &   ! -(b_E.S.)
              &, n_min(3),n_max(3),nglist,kg1,nngabc,newp)
      end if

      if(iprirmm >= 2) then
         if(allocated(nglist)) then
            write(nfout,'(" !!rmm nglist is allocated <<pre_calc_phase>>")')
         else
            write(nfout,'(" !!rmm nglist is not allocated <<pre_calc_phase>>")')
         end if
      end if

    end subroutine pre_calc_phase

    subroutine zf_listing
      integer       :: ia
      real(kind=DP) :: f(3)

      integer :: id_sname = -1
      call tstatc0_begin('zf_listing ',id_sname)
      do ia = rmm_precal_phase_matm+1, natm
         f = pos(ia,1:3)*PAI2
         call zf_list_s(n_min(1),n_max(1),rmm_precal_phase_matm,natm,f(1),ia,zfc1,zfs1)
         call zf_list_s(n_min(2),n_max(2),rmm_precal_phase_matm,natm,f(2),ia,zfc2,zfs2)
         call zf_list_s(n_min(3),n_max(3),rmm_precal_phase_matm,natm,f(3),ia,zfc3,zfs3)
         !                 -(b_Electronic_Structure)
      end do
      call tstatc0_end(id_sname)

    end subroutine zf_listing
  end subroutine m_ESrmm_renew_WF

! ======================================= added by K. Tagami ============== 11.0
  subroutine m_ESrmm_renew_WF_noncl( nfout, isolver, precon, dtim )
    integer,       intent(in) :: nfout,isolver,precon
    real(kind=DP), intent(in) :: dtim

    !  R   := -(H-eS)
    !  zaj_l(:,i,k,1:kimg) :=  | \Phi_{k,i} >
    !  phi (:,1:kimg,0)    :=  | \Phi_{k,i} >
    !  Rphi(:,1:kimg,0)    := R|\Phi_{k,i}>
    !  phi (:,1:kimg,1)    := KR |\Phi_{k,i}>
    !  Rphi(:,1:kimg,1)    := RKR|\Phi_{k,i}>
    ! when isolver == RMM3,
    !  phi (:,1:kimg,2)    := KRKR |\Phi_{k,i}>
    !  Rphi(:,1:kimg,2)    := RKRKR|\Phi_{k,i}>
    ! when isolver == RMM2P,
    !  phi (:,1:kimg,2)    := |\Phi^{-1}_{k,i}> = |\Phi^{m}_{k,i}> - |\Phi^{m-1}_{k,i}>
    !  Rphi(:,1:kimg,2)    := R|\Phi^{-1}_{k,i}>
    !
    ! output
    !  zaj_l := |\Psi^{m+1}_{k,i}>
    !         = \alpha_0 phi(:,:,0) + \alpha_1 phi(:,:,1) + \alpha_2 phi(:,:,2)
    ! where,
    !  K is preconditioning operator,
    !  alpha_0,alpha_1, and alpha_2 are decided from <Rphi_i|Rphi_j> and <phi_i|S|phi_j>
    !
    ! For details, refer to G. Kresse and J. Furthmuller, PRB54,11169,(1996)
    ! The procedure in a case of isolver == RMM2P has been developped originally by T. Uda.
    !

! === Debug!!! These arrays should be saved for alloc_zfc_zfs call!!! by T.Kato ==========
!   integer :: n_max(3), n_min(3)
    integer, save :: n_max(3) = 0, n_min(3) = 0
! ========================================================================================
    integer :: ispin,ik,ib,iter_rmm, nrmm, matm
    integer :: is, is1, is2, istmp, k2

    integer :: ib1

    real(kind=DP), pointer, dimension(:)   :: ekin
    real(kind=DP), pointer, dimension(:,:) :: phi0

    real(kind=DP), allocatable ::  afft_kt(:,:)
    real(kind=DP), allocatable ::  bfft_kt(:,:)
    real(kind=DP), allocatable :: vnlph_noncl(:,:,:,:)

    integer :: id_sname = -1

!!$    call tstatc0_begin('m_ESrmm_renew_WF ', id_sname,1)

    rr_avr = 0.d0
    call what_is_the_dimension_of_rmm() ! -(contained here) ->nrmm(= 2 or 3)

    call m_FFT_alloc_WF_work()  ! allocate(ftw)
    allocate( afft_kt(nfft,ndim_chgpot ) ); afft_kt = 0.0d0
    allocate( bfft_kt(nfft,ndim_spinor ) ); bfft_kt = 0.0d0
    allocate(vnlph_noncl(kg1,np_e,kimg,ndim_spinor)); vnlph_noncl = 0.0d0


    if (isolver==RMM2P.and.nrmm==3) call m_ESsd_diff_WFs_to_zaj_old(1.d0) !->zaj_old
    if (rmm_save_memory_mode == ON) then
       call alloc_zfc_zfs(nfout,n_min,n_max,rmm_precal_phase_matm,natm,nrmm,matm)
                                  ! -(m_ES_WF_by_RMM)
       call zf_listing()          ! -(contained here) -> zfc1,zfc2,zfc3,zfs1,zfs2,zfs3
    end if

    allocate(ekin0(kg1))
    ekin => ekin0(1:kg1)

! -----
    if(sw_fef==ON) then
       call phase_error_with_msg(nfout, 'kt : not supported ... fef ',__LINE__,__FILE__)
       call m_FEF_build_grad
    endif
! ----

    call m_ES_Vlocal_in_Rspace_noncl( afft_kt )
! --
    Loop_Kpoints : Do ik=1, kv3, ndim_spinor
       if ( map_k(ik) /= myrank_k ) cycle              ! MPI

       if ( iprirmm >= 3 .and. printable ) then
	  Do is=1, ndim_spinor
             call m_ES_wd_zaj_small_portion( nfout, ik+is-1," -- rmm --",10 )
          End do
       endif
       call m_pwBS_kinetic_energies( ik, vkxyz, ekin )     ! (diakin) ->ekin

       if ( rmm_save_memory_mode == ON ) then
          call alloc_phasecs_bW_phietc_noncl( nrmm, n_min, n_max, matm )
          ! allocations of bWr, bWi, phi and Rphi. (bWr,bWi) = <\beta |WF=zaj_l>
          call pre_calc_phase !-(contained here) -> phasec,phases,newp

          Loop_band: do ib1 = ista_e, iend_e, istep_e      ! MPI
             if ( rr_is_over_or_under(map_z(ib1),ik) == UNDER .and.sw_precalculate==OFF) cycle

             ib = map_z(ib1)

             Loop_res_dim: do iter_rmm = 0, nrmm-1
                if ( isolver == RMM2P .and. iter_rmm == nrmm-1 ) then
                   call zajold2zaj_phi2zaj_old_noncl( ib1 )
                   !-(contained here) zaj_old->zaj, phi0 ->zaj_old
                endif

                vnlph_noncl(:,ib,:,:) = 0.0d0
                Do is1=1, ndim_spinor
                   Do is2=1, ndim_spinor
                      istmp = ( is1 -1 )*ndim_spinor + is2
                      k2 = ik +is2 -1

                      call Vnonlocal_W_RMM_noncl( k2, ib1, istmp, &
                           &                      rmm_precal_phase_matm, iter_rmm )
                      !  ->vnlph_l, bWri_noncl

                      vnlph_noncl(:,ib,:,is1) = vnlph_noncl(:,ib,:,is1) &
                           &                        + vnlph_l(:,ib,:)
                   End do
                End do

! --------------------------------- under construction ---
!                 if(sw_hybrid_functional==ON) call m_ES_Vexx_W_RMM(ik,ib)
!                   if(sw_fef==ON) call m_FEF_add_grad_to_vnlph_RMM(ik,ib)
! ----------------------------------------------------

                Do is=1, ndim_spinor
                   call m_ES_WF_in_Rspace( ik+is-1, ib1, bfft_kt(:,is) ) ! ->bfft(=|WF(R)>)
                End do
                call m_FFT_Vlocal_W_noncl( afft_kt, bfft_kt, ndim_chgpot, ndim_spinor )
                ! (afft,bfft)->bfft
                Do is=1, ndim_spinor
                   call m_FFT_WF( ELECTRON, nfout,bfft_kt(:,is), DIRECT, ON )
                   ! bfft(R-space)->bfft(G-space)
                End do

                call rmm1_noncl( ib1, iter_rmm, nrmm )     ! -(contained here)
                !            zaj_l,vnlph_l,ekin,eko_l,bfft -> zaj_l,phi,Rphi

             end do Loop_res_dim
             call rmm_n_uda_noncl( ik, ib1, nrmm )
             ! -(contained here) phi,Rphi ->zaj_l, nrmm = 2 or 3
             !         Residual norm is also checked (->rr_is_oever_or_under)

          end do Loop_band
          call dealloc_phasecs_bW_phietc_noncl()

       else

          call alloc_bW_phi_Rphi_noncl(nrmm)
              ! allocations of bWr, bWi, phi and Rphi. (bWr,bWi) = <\beta |WF=zaj_l>

          Loop_res_dim_b: do iter_rmm = 0, nrmm-1
             if (isolver == RMM2P .and. iter_rmm == nrmm-1) then
                call zajold2zaj_phi2zaj_old_allnoncl
                !-(contained here) zaj_old->zaj, phi0 ->zaj_old,
                ! if rr_is_over_or_under /= UNDER
             endif

             vnlph_noncl(:,:,:,:) = 0.0d0
             Do is1=1, ndim_spinor
                Do is2=1, ndim_spinor
                   istmp = ( is1 -1 )*ndim_spinor + is2
                   k2 = ik +is2 -1
                   call Vnonlocal_W_RMMn_noncl( nfout, k2, istmp, iter_rmm )
                   ! -> vnlph_l,bWr,bWi

                   vnlph_noncl(:,:,:,is1) = vnlph_noncl(:,:,:,is1) &
                  &                     + vnlph_l(:,:,:)
                End do
             End do

! ------------------ under construction --
!             if(sw_fef==ON) call m_FEF_add_grad_to_vnlph(ik)
! ---------------------------------------

             do ib1 = ista_e, iend_e, istep_e
                if (rr_is_over_or_under(map_z(ib1),ik) == UNDER.and.sw_precalculate==OFF) cycle

! ------------------ under construction --
!                if(sw_hybrid_functional==ON) call m_ES_Vexx_W_RMM(ik,ib)
! ---------------------------------------
                Do is=1, ndim_spinor
                   call m_ES_WF_in_Rspace( ik+is-1, ib1, bfft_kt(:,is))
                                           ! ->bfft(=|WF(R)>)
                End do
                call m_FFT_Vlocal_W_noncl( afft_kt, bfft_kt, ndim_chgpot, ndim_spinor )

                Do is=1, ndim_spinor
                   call m_FFT_WF( ELECTRON, nfout,bfft_kt(:,is), DIRECT, ON )
                   ! bfft(R-space)->bfft(G-space)
                End do

                call rmm1_noncl( ib1, iter_rmm, nrmm )     ! -(contained here)
                !          zaj_l,vnlph_l,ekin,eko_l,bfft -> zaj_l,phi,Rphi
             end do
          end do Loop_res_dim_b

          Loop_band_b: do ib1 = ista_e, iend_e, istep_e
             if (rr_is_over_or_under(map_z(ib1),ik) /= UNDER) then
                call rmm_n_uda_noncl( ik, ib1, nrmm )
                ! -(contained here) phi,Rphi ->zaj_l, nrmm = 2 or 3
                !   Residual norm is also checked (->rr_is_oever_or_under)
             endif
          end do Loop_band_b

          call dealloc_bW_phi_Rphi_noncl()
       end if

!!$          call m_ES_betar_dot_WFs_4_each_k(nfout,ik)   !  -> fsr_l,fsi_l

       call orthonorm_or_norm_noncl( imGSrmm )  !-(contained here)

! ------------------- under construction --
!       if(sw_hybrid_functional==ON) call m_ES_EXX_eigenvalue_for_each_k(ispin,ik)
! --------------------------------------

       call m_ES_eigen_vals_each_k_noncl( ik, ekin, afft_kt, bfft_kt )

    end do Loop_kpoints

!$$#endif

    if (rmm_save_memory_mode == ON) call dealloc_zfc_zfs()

    call rr_avr_final_noncl        !-(c.h.)  rr_avr = rr_avr/(kv3/ndim_spinor*neg)

! ****************************************** This is important ******
    call m_ES_sort_eigen_vals_noncl()
! *******************************************************************

    deallocate(ekin0)
    call m_FFT_dealloc_WF_work()

    deallocate( afft_kt, bfft_kt )
    deallocate( vnlph_noncl )

!!$    call tstatc0_end(id_sname)
  contains

    subroutine orthonorm_or_norm_noncl(imGSrmm_t)
      integer, intent(in) :: imGSrmm_t
      if (mod(iteration_electronic - iteration_rmm_start, imGSrmm_t) == 0 ) then
         if (iprirmm >= 2 .and. ik == 1) write(nfout,'(" -- ORTHONORMALIZATION --")')
         call m_ES_MGS_4_each_k_noncl( nfout, ik, mode=ORTHONORMALIZATION )
      else
         if (iprirmm >= 2 .and. ik == 1) write(nfout,'(" -- NORMALIZATION --")')
         call m_ES_MGS_4_each_k_noncl( nfout, ik, mode=NORMALIZATION )
      end if
    end subroutine orthonorm_or_norm_noncl

    subroutine zajold2zaj_phi2zaj_old_allnoncl
      integer :: ib, ibto, is

      do ib = ista_e, iend_e, istep_e
         if (rr_is_over_or_under( map_z(ib),ik ) == UNDER.and.sw_precalculate==OFF) cycle

         Do is=1, ndim_spinor
           if (nrmm == 3) call m_ESsd_copy_zaj_old_to_zaj( ik+is-1,ib )
           ibto = map_z(ib)
           phi0 => phi_noncl(:,:,ibto,0,is)
           call m_ESsd_copy_phi_to_zaj_old( ik+is-1, ib, phi0 )  ! phi->zaj_old
         End do
      end do
    end subroutine zajold2zaj_phi2zaj_old_allnoncl

    subroutine zajold2zaj_phi2zaj_old_noncl( ib1 )
      integer, intent(in) :: ib1

      integer :: is

      Do is=1, ndim_spinor
        if (nrmm == 3) call m_ESsd_copy_zaj_old_to_zaj( ik+is-1,ib1 )
        phi0 => phi_noncl(:,:,1,0,is)
        call m_ESsd_copy_phi_to_zaj_old( ik+is-1, ib1, phi0 )  ! phi->zaj_old
      End do
    end subroutine zajold2zaj_phi2zaj_old_noncl

    subroutine rr_avr_final_noncl
      integer       :: i,j, n_rr_under, n_mpi
      real(kind=DP) :: rr_avr_mpi

      if(npes >= 2) then
         call mpi_allreduce(rr_avr, rr_avr_mpi,1,mpi_double_precision &
              & , mpi_sum,MPI_CommGroup,ierr)
!         call mpi_allreduce(rr_avr, rr_avr_mpi,1,mpi_double_precision &
!              & , mpi_sum,mpi_spin_group,ierr)
      else
         rr_avr_mpi = rr_avr
      end if

      rr_avr = rr_avr_mpi/( kv3/ndim_spinor *neg)

      if(iprirmm >= 1) &
           & write(nfout,'(" Residual Norm (sum_i <R_i|R_i>/sum_i) ( " &
           & ,i6," -th iter (electronic)) = ",d19.7)') iteration_electronic, rr_avr

      n_rr_under = 0
      do j = ista_k, iend_k
         do i = 1, np_e
            if(rr_is_over_or_under(map_z(i),j) == UNDER) n_rr_under = n_rr_under + 1
         end do
      end do
      if(npes >= 2) then
         call mpi_allreduce(n_rr_under,n_mpi, 1, mpi_integer,mpi_sum,MPI_CommGroup,ierr)
!         call mpi_allreduce(n_rr_under,n_mpi, 1, mpi_integer,mpi_sum,mpi_spin_group,ierr)
      else
         n_mpi = n_rr_under
      end if
      if(n_mpi > 0 .and. iprirmm >= 1) write(nfout,'(" Number of rr_under = ",i5)') n_mpi
    end subroutine rr_avr_final_noncl

    subroutine what_is_the_dimension_of_rmm
      if(isolver == RMM2 ) then
         nrmm = 2
      else if(isolver == RMM2p) then
         if(iteration_electronic == iteration_rmm_start) then
            nrmm = 2
         else
            nrmm = 3
         end if
      else if(isolver == RMM3) then
         nrmm = 3
      else
         if(printable) write(nfout,'(" !isolver",i3," is illegal ( m_ESrmm_renew_WF )")') isolver
         call phase_error_with_msg(nfout,'solver is illegal ( m_ESrmm_renew_WF )',__LINE__,__FILE__)
      endif
    end subroutine what_is_the_dimension_of_rmm

    subroutine rmm1_noncl( ibo, iter_rmm, nrmm )
      ! Revised by T. Yamasaki, 18th Sep. 2004
      !    <R_i|R_j> and <Phi_i|S|Phi_j> are calculated, then are substituted to
      !   rr_e and psp_e, respectively with not using
      integer, intent(in) :: ibo,iter_rmm,nrmm  ! nrmm = {2|3}
      !               iter_rmm = {0|1} (when nrmm == 2)}, {0|1|2} (when nrmm == 3)
      real(kind=DP) :: evr,devr,dnm,evi,e1,devi
      real(kind=DP) :: rrr, rri, pspr, pspi

      real(kind=DP), pointer, dimension(:,:,:) :: phi_t_noncl
      real(kind=DP), pointer, dimension(:,:,:) :: Rphi_t_noncl
                                                       !d(iba(ik),kimg,ndim_spinor)
      real(kind=DP), allocatable, dimension(:) :: p
      integer       :: i, i1, ib, ibt, ir, ii, j
      integer :: k1

      real(kind=DP) :: ctmp1, ctmp2

      allocate(p(kg1))

      if(iter_rmm == nrmm-1) then
         allocate( phi_t_noncl(iba(ik),kimg,ndim_spinor))
         allocate( Rphi_t_noncl(iba(ik),kimg,ndim_spinor))
      end if

      ib = map_z(ibo)                  ! MPI
      if ( rmm_save_memory_mode == ON ) then
         ibt = 1
      else
         ibt = ib
      end if
      dnm = 1.d0/product(fft_box_size_WF(1:3,1))
      call m_ES_decide_precon_factor( precon, ik, ibo, ekin, p ) ! ->p(1:iba(ik))

      if (kimg == 1) then
         if (iter_rmm <= nrmm-2) then
            Do is=1, ndim_spinor
              k1 = ik +is -1
              do i = 1, iba(ik)
                 i1    = igf(nbase(i,ik))
                 evr   = zaj_l( i, ib, k1, 1 )
                 devr  = ( ekin(i) - eko_l(ib,ik) )*evr + bfft_kt(i1,is)*dnm &
	&               +vnlph_noncl( i, ib, 1, is )
                 zaj_l( i, ib, k1, 1 )    = -p(i)*devr
                 phi_noncl( i, 1, ibt, iter_rmm, is )  =  evr
                 Rphi_noncl( i, 1, ibt,iter_rmm, is )  = -devr
              end do
#ifdef SAVE_FFT_TIMES
              if(sw_save_fft == ON) status_saved_phifftr(ib,k1) = OLD
#endif
            End do
            if ( k_symmetry(ik) == GAMMA ) then
               do i1 = 0, iter_rmm
                 ctmp1 = 0.0d0;  ctmp2 = 0.0d0
                 Do is=1, ndim_spinor
                    ctmp1 = ctmp1 + Rphi_noncl(1,1,ibt,i1,is) &
	&                        *Rphi_noncl(1,1,ibt,iter_rmm,is) &
	&                 +2.d0 *dot_product( Rphi_noncl(2:iba(ik),1,ibt,i1,is),&
	&                                     Rphi_noncl(2:iba(ik),1,ibt,iter_rmm,is) )
!
                    ctmp2 = ctmp2 + phi_noncl(1,1,ibt,i1,is) &
	&                          *phi_noncl(1,1,ibt,iter_rmm,is) &
	&                 + 2.d0 *dot_product( phi_noncl(2:iba(ik),1,ibt,i1,is), &
	&                                      phi_noncl(2:iba(ik),1,ibt,iter_rmm,is) )
                 End do
                 rr_e( i1, iter_rmm, ibt, 1 ) = ctmp1
                 psp_e(i1, iter_rmm, ibt, 1 ) = ctmp2
               end do
            else
               do i1 = 0, iter_rmm
                 ctmp1 = 0.0d0;  ctmp2 = 0.0d0
                 Do is=1, ndim_spinor
	           ctmp1 = ctmp1 +dot_product( Rphi_noncl(1:iba(ik),1,ibt,i1,is),&
	&                                      Rphi_noncl(1:iba(ik),1,ibt,iter_rmm,is))
                   ctmp2 = ctmp2 +dot_product( phi_noncl(1:iba(ik),1,ibt,i1,is),&
	&                                      phi_noncl(1:iba(ik),1,ibt,iter_rmm,is) )
                 End do
                 rr_e( i1, iter_rmm, ibt,1 ) = ctmp1
                 psp_e(i1, iter_rmm, ibt,1 ) = ctmp2
               end do
            end if
         else if ( iter_rmm == nrmm-1 ) then
            Do is=1, ndim_spinor
              k1 = ik + is -1
              do i = 1, iba(ik)
                 i1 = igf(nbase(i,ik))
                 phi_t_noncl(i,1,is)  = zaj_l( i, ib, k1, 1 )
                 Rphi_t_noncl(i,1,is) = -( (ekin(i)-eko_l(ib,ik))*phi_t_noncl(i,1,is)&
	&                               + bfft_kt(i1,is)*dnm +vnlph_noncl(i,ib,1,is) )
               !   when iter_rmm == nrmm-1, zaj_l (to which KR|Psi> is stored) is not
               !   going to be used recursively, so we do not need to store evr to zaj_l,
               !   instead we ues zaj_l as |phi_{iter_rmm}>.
              end do
            End do

            if ( k_symmetry(ik) == GAMMA ) then
               do i1 = 0, iter_rmm-1
                 ctmp1 = 0.0d0;  ctmp2 = 0.0d0
                 DO is=1, ndim_spinor
	           ctmp1 = ctmp1 + Rphi_noncl(1,1,ibt,i1,is) &
	&                         *Rphi_t_noncl(1,1,is) &
        &                 + 2.d0 *dot_product( Rphi_noncl(2:iba(ik),1,ibt,i1,is), &
	&                                      Rphi_t_noncl(2:iba(ik),1,is) )
                   ctmp2 = ctmp2 + phi_noncl(1,1,ibt,i1,is) &
	&                         *phi_t_noncl(1,1,is) &
        &                 + 2.d0 *dot_product( phi_noncl(2:iba(ik),1,ibt,i1,is),&
	&                                      phi_t_noncl(2:iba(ik),1,is) )
                 End do
                 rr_e(i1,iter_rmm,ibt,1)  = ctmp1
                 psp_e(i1,iter_rmm,ibt,1) = ctmp2
               end do

               ctmp1 = 0.0d0;  ctmp2 = 0.0d0
               Do is=1, ndim_spinor
                 ctmp1 = ctmp1 + Rphi_t_noncl(1,1,is) *Rphi_t_noncl(1,1,is) &
        &                      + 2.d0 *dot_product( Rphi_t_noncl(2:iba(ik),1,is),&
	&                                           Rphi_t_noncl(2:iba(ik),1,is) )
                 ctmp2 = ctmp2 + phi_t_noncl(1,1,is) *phi_t_noncl(1,1,is) &
        &                      + 2.d0 *dot_product( phi_t_noncl(2:iba(ik),1,is), &
	&                                           phi_t_noncl(2:iba(ik),1,is) )
               End do
               rr_e(iter_rmm,iter_rmm,ibt,1)  = ctmp1
               psp_e(iter_rmm,iter_rmm,ibt,1) = ctmp2

            else
               do i1 = 0, iter_rmm-1
                 ctmp1 = 0.0d0;  ctmp2 = 0.0d0
                 DO is=1, ndim_spinor
	           ctmp1 = ctmp1 + dot_product( Rphi_noncl(1:iba(ik),1,ibt,i1,is),&
	&                                       Rphi_t_noncl(1:iba(ik),1,is) )
                   ctmp2 = ctmp2 + dot_product( phi_noncl(1:iba(ik),1,ibt,i1,is),&
	&                                       phi_t_noncl(1:iba(ik),1,is) )
                 End do
                 rr_e(i1,iter_rmm,ibt,1)  = ctmp1
                 psp_e(i1,iter_rmm,ibt,1) = ctmp2
               end do

               ctmp1 = 0.0d0;  ctmp2 = 0.0d0
               Do is=1, ndim_spinor
	         ctmp1 = ctmp1 + dot_product( Rphi_t_noncl(1:iba(ik),1,is),&
	&                                     Rphi_t_noncl(1:iba(ik),1,is) )
                 ctmp2 = ctmp2 + dot_product( phi_t_noncl(1:iba(ik),1,is), &
	&                                     phi_t_noncl(1:iba(ik),1,is) )
               End do
               rr_e(iter_rmm,iter_rmm,ibt,1)  = ctmp1
               psp_e(iter_rmm,iter_rmm,ibt,1) = ctmp2
            end if
         else
            if (printable) then
	       write(nfout,*) &
                    & "  iter_rmm > nrmm-1 ( illegal relation) <<m_ES_WF_by_RMM.rmm1_noncl>>"
               call phase_error_with_msg(nfout, '   iter_rmm > nrmm-1 ( illegal relation) <<m_ES_WF_by_RMM.rmm1_noncl>>' &
                                        , __LINE__,__FILE__)
            end if
         endif

      else if (kimg == 2) then
         ir = iter_rmm

         if (iter_rmm <= nrmm-2) then
            Do is=1, ndim_spinor
               k1 = ik + is -1
               do i = 1, iba(ik)
                 i1    = igf(nbase(i,ik))

                 evr   = zaj_l( i, ib, k1, 1 );
	         evi   = zaj_l( i, ib, k1, 2 )
                 e1    = ekin(i) - eko_l(ib,ik)
                 devr  = e1 *evr +bfft_kt(2*i1-1,is)*dnm +vnlph_noncl(i,ib,1,is)
                 devi  = e1 *evi +bfft_kt(2*i1,  is)*dnm +vnlph_noncl(i,ib,2,is)

                 zaj_l( i, ib, k1, 1 ) = -p(i)*devr
                 zaj_l( i, ib, k1, 2 ) = -p(i)*devi
                 phi_noncl( i, 1, ibt, iter_rmm, is ) =  evr;
	         phi_noncl( i, 2, ibt, iter_rmm, is ) =  evi
                 Rphi_noncl( i,1, ibt, iter_rmm, is ) = -devr
                 Rphi_noncl( i,2, ibt, iter_rmm, is ) = -devi
               end do
#ifdef SAVE_FFT_TIMES
               if(sw_save_fft == ON) status_saved_phifftr(ib,k1) = OLD
#endif
            End do

            do i1 = 0, iter_rmm
               rrr = 0.d0; rri = 0.d0; pspr = 0.d0; pspi = 0.d0

               if (k_symmetry(ik) == GAMMA) then
                  Do is=1, ndim_spinor
                     do i = 2, iba(ik)
                       rrr = rrr + Rphi_noncl(i,1,ibt,i1,is) &
	&                         *Rphi_noncl(i,1,ibt,iter_rmm,is) &
	&                        + Rphi_noncl(i,2,ibt,i1,is) &
	&                         *Rphi_noncl(i,2,ibt,iter_rmm,is)
                       pspr = pspr + phi_noncl(i,1,ibt,i1,is) &
	&                           *phi_noncl(i,1,ibt,iter_rmm,is) &
	&                          + phi_noncl(i,2,ibt,i1,is) &
	&                           *phi_noncl(i,2,ibt,iter_rmm,is)
                    end do
                  End do

                  ctmp1 = 0.0d0;  ctmp2 = 0.0d0
                  Do is=1, ndim_spinor
                     ctmp1 = ctmp1 + Rphi_noncl(1,1,ibt,i1,is)&
	&                           *Rphi_noncl(1,1,ibt,iter_rmm,is)
	             ctmp2 = ctmp2 + phi_noncl(1,1,ibt,i1,is) &
	&                           *phi_noncl(1,1,ibt,iter_rmm,is)
	          End do
                  rrr  = 2.d0*rrr  + ctmp1
                  pspr = 2.d0*pspr + ctmp2

               else
                  Do is=1, ndim_spinor
                     do i = 1, iba(ik)
                       rrr = rrr + Rphi_noncl(i,1,ibt,i1,is) &
	&                         *Rphi_noncl(i,1,ibt,iter_rmm,is) &
	&                        + Rphi_noncl(i,2,ibt,i1,is) &
	&                         *Rphi_noncl(i,2,ibt,iter_rmm,is)
                       rri = rri + Rphi_noncl(i,1,ibt,i1,is) &
	&                         *Rphi_noncl(i,2,ibt,iter_rmm,is) &
	&                        - Rphi_noncl(i,2,ibt,i1,is) &
	&                         *Rphi_noncl(i,1,ibt,iter_rmm,is)
                       pspr = pspr + phi_noncl(i,1,ibt,i1,is) &
	&                           *phi_noncl(i,1,ibt,iter_rmm,is) &
	&                          + phi_noncl(i,2,ibt,i1,is) &
	&                           *phi_noncl(i,2,ibt,iter_rmm,is)
                       pspi = pspi + phi_noncl(i,1,ibt,i1,is) &
	&                           *phi_noncl(i,2,ibt,iter_rmm,is) &
	&                          - phi_noncl(i,2,ibt,i1,is) &
	&                           *phi_noncl(i,1,ibt,iter_rmm,is)
                     end do
                  End do
               end if
!
                rr_e(i1,iter_rmm,ibt,1) = rrr;   rr_e(i1,iter_rmm,ibt,2) = rri
               psp_e(i1,iter_rmm,ibt,1) = pspr; psp_e(i1,iter_rmm,ibt,2) = pspi
            end do

         else if(iter_rmm == nrmm-1) then
            Do is=1, ndim_spinor
               k1 = ik + is -1
               do i = 1, iba(ik)
                 i1    = igf(nbase(i,ik))
                 phi_t_noncl(i,1,is) = zaj_l( i, ib, k1, 1 )
	         phi_t_noncl(i,2,is) = zaj_l( i, ib, k1, 2 )
                 e1    = ekin(i) - eko_l(ib,ik)
                 Rphi_t_noncl(i,1,is) = - ( e1*phi_t_noncl(i,1,is) &
	&                                  +bfft_kt(2*i1-1,is)*dnm &
	&                                  +vnlph_noncl(i,ib,1,is) )
                 Rphi_t_noncl(i,2,is) = - ( e1*phi_t_noncl(i,2,is) &
	&                                  +bfft_kt(2*i1,  is)*dnm &
	&                                  +vnlph_noncl(i,ib,2,is) )
              end do
            End do

            do i1 = 0, iter_rmm-1
               rrr = 0.d0; rri = 0.d0; pspr = 0.d0; pspi = 0.d0

               if (k_symmetry(ik) == GAMMA) then
                  Do is=1, ndim_spinor
                     do i = 2, iba(ik)
                       rrr = rrr + Rphi_noncl(i,1,ibt,i1,is) &
	&                         *Rphi_t_noncl(i,1,is) &
	&                        + Rphi_noncl(i,2,ibt,i1,is) &
	&                         *Rphi_t_noncl(i,2,is)
                       pspr = pspr + phi_noncl(i,1,ibt,i1,is) &
	&                           *phi_t_noncl(i,1,is) &
	&                          + phi_noncl(i,2,ibt,i1,is) &
	&                           *phi_t_noncl(i,2,is)
                    end do
                  End do

                  ctmp1 = 0.0d0;  ctmp2 = 0.0d0
	          Do is=1, ndim_spinor
                    ctmp1 = ctmp1 + Rphi_noncl(1,1,ibt,i1,is) &
	&                          *Rphi_t_noncl(1,1,is)
                    ctmp2 = ctmp2 + phi_noncl(1,1,ibt,i1,is) &
	&                          *phi_t_noncl(1,1,is)
                  End do
                  rrr = 2.d0*rrr  + ctmp1
                  pspr= 2.d0*pspr + ctmp2

               else
                  Do is=1, ndim_spinor
                     do i = 1, iba(ik)
                        rrr = rrr + Rphi_noncl(i,1,ibt,i1,is) &
	&                          *Rphi_t_noncl(i,1,is) &
	&                         + Rphi_noncl(i,2,ibt,i1,is) &
	&                          *Rphi_t_noncl(i,2,is)
                        rri = rri + Rphi_noncl(i,1,ibt,i1,is) &
        &                          *Rphi_t_noncl(i,2,is) &
	&                         - Rphi_noncl(i,2,ibt,i1,is) &
	&                          *Rphi_t_noncl(i,1,is)
                        pspr = pspr + phi_noncl(i,1,ibt,i1,is) &
	&                            *phi_t_noncl(i,1,is) &
	&                           + phi_noncl(i,2,ibt,i1,is) &
	&                            *phi_t_noncl(i,2,is)
                        pspi = pspi + phi_noncl(i,1,ibt,i1,is) &
	&                            *phi_t_noncl(i,2,is) &
	&                           - phi_noncl(i,2,ibt,i1,is) &
	&                            *phi_t_noncl(i,1,is)
                    end do
                  End do
               end if

               rr_e(i1,iter_rmm,ibt,1) = rrr;   rr_e(i1,iter_rmm,ibt,2) = rri
               psp_e(i1,iter_rmm,ibt,1) = pspr; psp_e(i1,iter_rmm,ibt,2) = pspi
            end do

            rrr = 0.d0; rri = 0.d0; pspr = 0.d0; pspi = 0.d0

! ---oo
            if (k_symmetry(ik) == GAMMA) then
               Do is=1, ndim_spinor
                  do i = 2, iba(ik)
                     rrr = rrr + Rphi_t_noncl(i,1,is) &
	&                       *Rphi_t_noncl(i,1,is) &
	&                      + Rphi_t_noncl(i,2,is) &
	&                       *Rphi_t_noncl(i,2,is)
                     pspr = pspr + phi_t_noncl(i,1,is) &
	&                         *phi_t_noncl(i,1,is) &
	&                        + phi_t_noncl(i,2,is) &
	&                         *phi_t_noncl(i,2,is)
                 end do
               End do

	       ctmp1 = 0.0d0;  ctmp2 = 0.0d0
               Do is=1, ndim_spinor
	          ctmp1 = ctmp1 + Rphi_t_noncl(1,1,is)*Rphi_t_noncl(1,1,is)
                  ctmp2 = ctmp2 + phi_t_noncl(1,1,is) *phi_t_noncl(1,1,is)
               End do
                rrr = 2.d0*rrr  +ctmp1
               pspr = 2.d0*pspr +ctmp2

            else
               Do is=1, ndim_spinor
                  do i = 1, iba(ik)
                     rrr = rrr + Rphi_t_noncl(i,1,is) &
	&                       *Rphi_t_noncl(i,1,is) &
	&                      + Rphi_t_noncl(i,2,is) &
	&                       *Rphi_t_noncl(i,2,is)
                     rri = rri + Rphi_t_noncl(i,1,is) &
	&                       *Rphi_t_noncl(i,2,is) &
	&                      - Rphi_t_noncl(i,2,is) &
	&                       *Rphi_t_noncl(i,1,is)
                     pspr = pspr + phi_t_noncl(i,1,is) &
	&                         *phi_t_noncl(i,1,is) &
	&                        + phi_t_noncl(i,2,is) &
	&                         *phi_t_noncl(i,2,is)
                     pspi = pspi + phi_t_noncl(i,1,is) &
	&                         *phi_t_noncl(i,2,is) &
	&                        - phi_t_noncl(i,2,is) &
 	&                         *phi_t_noncl(i,1,is)
                 end do
               End do
            end if
! ---oo
             rr_e(iter_rmm,iter_rmm,ibt,1) = rrr;   rr_e(iter_rmm,iter_rmm,ibt,2) = rri
            psp_e(iter_rmm,iter_rmm,ibt,1) = pspr; psp_e(iter_rmm,iter_rmm,ibt,2) = pspi

         else
            if (printable) then
	       write(nfout,*) &
                    &        "  iter_rmm > nrmm-1 ( illegal relation) <<m_ES_WF_by_RMM.rmm1_noncl>>"
               call phase_error_with_msg(nfout, '   iter_rmm > nrmm-1 ( illegal relation) <<m_ES_WF_by_RMM.rmm1_noncl>>'&
                                        , __LINE__, __FILE__)
            end if
         endif
      end if

      if (iter_rmm == nrmm-1) then
         deallocate(Rphi_t_noncl)
         deallocate(phi_t_noncl)
      end if
      deallocate(p)
    end subroutine rmm1_noncl

    recursive subroutine rmm_n_noncl( ibo, nrmm, dtim )
!
!    * Rewritten by T. Uda and T. Yamasaki, 19th Mar. 2003
!      + alpha(nrmm-1) -> alpha(nrmm-1,kimg)
!
      integer,       intent(in) :: ibo,nrmm
      real(kind=DP), intent(in) :: dtim

      integer                   :: ip(nrmm)
      real(DP) :: ymat(0:nrmm-1,0:nrmm-1,kimg)   &
           &     ,vector(0:nrmm-1,0:nrmm-1,kimg) &
           &     ,rr(0:nrmm-1,0:nrmm-1,kimg),psp(0:nrmm-1,0:nrmm-1,kimg)
      real(DP) :: eigenr(0:nrmm-1), eigeni(0:nrmm-1)
      real(DP) :: alpha(nrmm-1,kimg),ww1(nrmm*kimg),ww2(nrmm*kimg),ww3(nrmm*kimg)
      integer  :: ierrdecomp, ib

      ib = map_z(ibo)                 ! MPI
      call rr_and_psp_from_Rphi_etc_noncl( ib, nrmm, rr, psp )  !-(m_ES_WF_by_RMM)
      !   rr(i,j) = <R_i|R_j>; psp(i,j) = <Phi_i|S|Phi_j>
      !                            S = 1 + \sum_{ij}q_{ij}|b_i><b_j|
      call normalization_of_rr_and_psp( nrmm, rr, psp )
      !   rr(i,j) = rr(i,j)/rr(0,0), psp(i,j) = psp(i,j)/psp(0,0)

!* Solve the Eigen value problem of Bx = epsilon*Ax
!* Here B is the residual matrix: B_ij = <R^i|R^j>,
!* and A is <phi^i|S|phi^j>.
!* It is equivallent to eigen value problem of Yx = epsilon*x
!* where AY = B.
      call LUdecomp(nrmm,psp,ww1,ip,ierrdecomp)
      if (ierrdecomp /= 0) then
         if(printable) write(nfout,*) 'LU decomposition is impossible.';  alpha = 0.0d0
         alpha(1,1) = dtim
      else
         call solveAYB(nrmm,psp,rr,ymat,ip)      ! Solve AY = B
         call rg_or_cg(nrmm,ymat,eigenr,eigeni,vector,ww1,ww2,ww3)
         call decide_alpha(kimg,dtim,nrmm,eigenr,eigeni,vector,alpha)
         !                    -(rmmsubs) ->alpha
      end if

      if (nrmm == 3 .and. ierrdecomp /= 0) then
         call rmm_n_noncl(ibo,2,dtim)
      else
         if(iprirmm >= 2) call wd_alpha(nrmm,alpha)
         call evolve_WF_using_Residuals_noncl( ib,nrmm,alpha) ! -(contained subr. m_ESrmm_renew_WF)
         if(iprirmm>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')
         !    phi,alpha ->zaj_l(:,ib,ik,:)
         !   |WF> = |phi_0> + alpha(1)|phi_1> + alpha(2)|phi_2> , when nrmm == 3
      end if
    end subroutine rmm_n_noncl

    subroutine evolve_WF_using_Residuals_noncl( ib, nrmm, alpha )
!
!  * Rewritten by T. Uda and T. Yamasaki, 19th Mar. 2003
!      + alpha(nrmm-1) -> alpha(nrmm-1,kimg)
!
!  * Revised by T. Yamasaki, 18th Sep. 2004
!      + being changed not to use phi(:,:,:,nrmm-1) for saving the memory allocation area.
!      + phi(:,:,:,nrmm-1) is substituted to zaj_l before the program reached to this
!       subroutine.
!
      integer, intent(in) ::                             ib, nrmm
      real(kind=DP),intent(in),dimension(nrmm-1,kimg) :: alpha
      integer ::                                         i, ibt, nr, ir, ii, j, ig
      integer :: is, k1

      real(kind=DP) ::                                   zr,zi

      if (rmm_save_memory_mode == ON) then
         ibt = 1
      else
         ibt = ib
      end if

      if (kimg==1) then
         ! previously in sub. rmm1, zaj_l(:,ib,ik,:) <= phi(:,:,ibt,nrmm-1)
         ! Following calculation is equivalent to
         !  zaj_l = phi(:,:,ibt,0) + alpha(nrmm-1)*phi(:,:,ibt,nrmm-1)
         !    |WF> = |phi_0> + alpha(nrmm-1)|phi_{nrmm-1}>

         Do is=1, ndim_spinor
            k1 = ik +is -1
            do i = 1, iba(ik)
               zaj_l( i, ib, k1, 1 ) = phi_noncl( i, 1, ibt, 0, is ) &
	&                            + alpha(nrmm-1,1) *zaj_l( i, ib, k1, 1 )
            end do
#ifdef SAVE_FFT_TIMES
            if(sw_save_fft == ON) status_saved_phifftr(ib,k1) = OLD
#endif
         End do
         do ir = 1, nrmm-2
            Do is=1, ndim_spinor
               k1 = ik +is -1
               do i = 1, iba(ik)
                  zaj_l( i, ib, k1, 1 ) = zaj_l( i, ib, k1, 1 ) &
	&                               + alpha(ir,1) *phi_noncl( i, 1, ibt, ir, is )
               end do
            End do
         end do

      else if(kimg==2) then
         ! previously in <sub. rmm1>, zaj_l(:,ib,ik,:) <-- phi(:,:,ibt,nrmm-1)
         nr = nrmm-1
         Do is=1, ndim_spinor
            k1 = ik +is -1
            do i = 1, iba(ik)
               !   Calculation in this do-loop is equivalent to
               !  zaj_l(:,ib,ik,:) = phi(:,:,ibt,0) + alpha(nrmm-1)*phi(:,:,ibt,nrmm-1),
               !  namely, |WF> = |phi_0> + alpha(nrmm-1)|phi_{nrmm-1}>

               zr = zaj_l( i, ib, k1, 1 )
               zi = zaj_l( i, ib, k1, 2 )
               zaj_l( i, ib, k1, 1 ) = phi_noncl( i, 1, ibt, 0, is ) &
	&                            + alpha(nr,1)*zr - alpha(nr,2)*zi
               zaj_l( i, ib, k1, 2 ) = phi_noncl( i, 2, ibt, 0, is ) &
	&                            + alpha(nr,1)*zi + alpha(nr,2)*zr
            end do
#ifdef SAVE_FFT_TIMES
            if(sw_save_fft == ON) status_saved_phifftr(ib,k1) = OLD
#endif
         End do

         do i = 1, nrmm-2
            Do is=1, ndim_spinor
               k1 = ik +is -1
               do ig = 1, iba(ik)
                  zaj_l( ig, ib, k1, 1 ) = zaj_l( ig, ib,k1, 1 ) &
	&                                + alpha(i,1) *phi_noncl( ig, 1, ibt, i, is )&
	&                                - alpha(i,2) *phi_noncl( ig, 2, ibt, i, is )
                  zaj_l( ig, ib, k1, 2 ) = zaj_l( ig, ib, k1, 2 )&
	&                                + alpha(i,1) *phi_noncl( ig, 2, ibt, i, is )&
	&                                + alpha(i,2) *phi_noncl( ig, 1, ibt, i, is )
               end do
            End do
         end do
      else
         call phase_error_with_msg(nfout, ' kimg is illegal (evolve_WF_using_Residuals_noncl)',__LINE__,__FILE__)
      end if

      if (iprirmm >= 3) then
         ibt = 100
         if (ibt > iba(ik)) ibt = iba(ik)
         if (ik <= 2 .or. ik >= kv3-1 ) then
            write(nfout,'(" !rmm  zaj_l in <<evolve_WF_using_Residuals_noncl>>")')
            Do is=1, ndim_spinor
              write(nfout,'(" !rmm  ik, ib = ",2i8)') ik +is -1, ib
              if (kimg == 2) write(nfout,'(" !rmm -- real part --")')
              write(nfout,'(" !rmm  ",8f8.4)') (zaj_l(i,ib,ik+is-1,1),i=1,ibt)
              if (kimg == 2) write(nfout,'(" !rmm -- imaginary part --")')
              write(nfout,'(" !rmm  ",8f8.4)') (zaj_l(i,ib,ik+is-1,2),i=1,ibt)
            End do
         end if
      end if

    end subroutine evolve_WF_using_Residuals_noncl

    subroutine wd_rr_and_psp(nrmm,rr,psp)
      integer, intent(in)                                        :: nrmm
      real(kind=DP),intent(in),dimension(0:nrmm-1,0:nrmm-1,kimg) :: rr,psp
      integer :: i, j
      write(nfout,'("!rr =",6d15.7)') ((rr(i,j,1),j=i,nrmm-1),i=0,nrmm-1)
      write(nfout,'("!psp=",6d15.7)') ((psp(i,j,1),j=i,nrmm-1),i=0,nrmm-1)
    end subroutine wd_rr_and_psp

    subroutine wd_alpha(nrmm,alpha)
      integer, intent(in)                          :: nrmm
      real(kind=DP), intent(in), dimension(nrmm-1,kimg) :: alpha

      write(nfout,'(" ! alpha(1) = ", f20.10)') alpha(1,1)
      if(nrmm == 3) write(nfout,'(" ! alpha(2) = ", f20.10)') alpha(2,1)
    end subroutine wd_alpha

    subroutine rmm_n_uda_noncl( ik, ibo, nrmm )
      integer,       intent(in) :: ik,ibo, nrmm

      real(DP) :: rr(0:nrmm-1,0:nrmm-1,kimg),psp(0:nrmm-1,0:nrmm-1,kimg)
      real(DP) :: alpha(nrmm-1,kimg)
      integer  :: ib
      integer :: id_sname = -1
      call tstatc0_begin('rmm_n_uda ', id_sname)

      ib = map_z(ibo)
      call rr_and_psp_from_rr_e_etc_noncl( ib, nrmm, rr, psp )
                                             ! -(m_ES_WF_by_RMM)
      !   rr(i,j) = <R_i|R_j>; psp(i,j) = <Phi_i|S|Phi_j>
      !                            S = 1 + \sum_{ij}q_{ij}|b_i><b_j|

      rr_avr = rr_avr + rr(0,0,1)
      if (rr(0,0,1) <= rr_Critical_Value) rr_is_over_or_under(ib,ik) = UNDER

    !!$  call normalization_of_rr_and_psp(nrmm,rr,psp)
      !   rr(i,j) = rr(i,j)/rr(0,0), psp(i,j) = psp(i,j)/psp(0,0)

      if (iprirmm >= 2) call wd_rr_and_psp(nrmm,rr,psp)

      if (kimg == 1) then
         if(nrmm == 2) then
            call rmm2_uda(iprirmm,nrmm,rr,psp,alpha)
         else
            call rmm3_uda(iprirmm,nrmm,rr,psp,alpha)
         end if
      else
         if(nrmm == 2) then
            call crmm2_uda(iprirmm,nrmm,rr,psp,alpha)
         else
            call crmm3_uda(iprirmm,nrmm,rr,psp,alpha)
         end if
      end if

      if(iprirmm >= 2) call wd_alpha(nrmm,alpha)
      call evolve_WF_using_Residuals_noncl(ib,nrmm,alpha)
                  ! -(contained subr. m_ESrmm_renew_WF)
      !   phi,alpha ->zaj_l(:,ib,ik,:)
      !   |WF> = |phi_0> + alpha(1)|phi_1> + alpha(2)|phi_2> , when nrmm == 3
      call tstatc0_end(id_sname)
    end subroutine rmm_n_uda_noncl

    subroutine rg_or_cg(nrmm,ymat,eigenr,eigeni,vector,ww1,ww2,ww3)
      integer, intent(in) :: nrmm
      real(DP), dimension(0:nrmm-1,0:nrmm-1,kimg):: ymat
      real(DP), intent(out),dimension(0:nrmm-1)  :: eigenr, eigeni
      real(DP), intent(out),dimension(0:nrmm-1,0:nrmm-1,kimg):: vector
      real(DP), intent(out),dimension(nrmm*kimg) :: ww1,ww2,ww3
      integer                                    :: ier
      if(kimg == 1) then
         call rg_eispack(nrmm,nrmm,ymat,eigenr,eigeni,1,vector,ww1,ww2,ier)
      else if(kimg == 2) then
         call cg_eispack(nrmm,nrmm,ymat(0,0,1),ymat(0,0,kimg),eigenr,eigeni,1 &
              &, vector(0,0,1),vector(0,0,kimg),ww1,ww2,ww3,ier)
      end if
      if(ier /= 0 .and. printable) write(nfout,*) 'ier(rg_or_cg) = ',ier
    end subroutine rg_or_cg

    subroutine solveAYB(nrmm,psp,rr,ymat,ip)
      integer, intent(in) :: nrmm
      real(DP),intent(in), dimension(0:nrmm-1,0:nrmm-1,kimg):: psp,rr
      real(DP),intent(out),dimension(0:nrmm-1,0:nrmm-1,kimg):: ymat
      integer, intent(out),dimension(nrmm) :: ip(nrmm)

      if(kimg == 1) then
         call rsolve(nrmm,nrmm,psp,rr,ymat,ip)
      else
         call csolve2(nrmm,nrmm,psp,rr,ymat,ip)
      end if
    end subroutine solveAYB

    subroutine LUdecomp(nrmm,psp,ww1,ip,ier)
      integer, intent(in) :: nrmm
      real(DP),intent(inout), dimension(0:nrmm-1,0:nrmm-1,kimg):: psp
      real(DP),intent(out),dimension(nrmm*kimg)                :: ww1
      integer, intent(out)                 :: ip(nrmm), ier

      complex(CMPLDP)                      :: cpsp(nrmm*nrmm)
      integer i, j, n
      if(kimg == 1) then
         call rdecomp(nrmm,psp,ww1,ip,ier)
      else
         n = 0
         do j = 0, nrmm-1
            do i = 0, nrmm-1
               n = n + 1
               cpsp(n) = cmplx(psp(i,j,1),psp(i,j,2))
            end do
         end do
         call cdecomp(nrmm,cpsp,ww1,ip,ier)
      end if
    end subroutine LUdecomp

    subroutine rr_and_psp_from_Rphi_etc_noncl(ib_t,nrmm,rr,psp)
      integer, intent(in)       :: ib_t,nrmm
      real(DP),intent(out),dimension(0:nrmm-1,0:nrmm-1,kimg)::rr,psp
      real(DP), dimension(kimg) :: sumqff

      real(DP) :: ctmp1, ctmp2
      real(DP) ::                  rrr,rri,pspr,pspi
      integer  ::                  j,i,m,ri, ib
      integer :: is_tmp

      if (rmm_save_memory_mode == ON) then
         ib = 1
      else
         ib = ib_t
      end if

      rr = 0.d0; psp = 0.d0
      do j = 0, nrmm-1
         do i = 0,j
            call phsphn_noncl( i, j, ib, kimg, sumqff )
                    ! ->sumqff (= \sum_{n1,n2}q_{n1,n2}<W|b_{nl}><b_{n2}|W>
            if (kimg == 1) then
               ctmp1 = 0.0d0;   ctmp2 = 0.0d0
               Do is=1, ndim_spinor
	          ctmp1 = ctmp1 + dot_product( Rphi_noncl(1:iba(ik),1,ib,i,is),&
	&                                      Rphi_noncl(1:iba(ik),1,ib,j,is) )
                  ctmp2 = ctmp2 + dot_product( phi_noncl(1:iba(ik),1,ib,i,is), &
	&                                      phi_noncl(1:iba(ik),1,ib,j,is) )
               End do
               rr(i,j,1)  = ctmp1;     psp(i,j,1) = ctmp2

            else if(kimg == 2) then
               rrr = 0.d0;   rri = 0.d0;     pspr = 0.d0;  pspi = 0.d0

               Do is=1, ndim_spinor
                  do m = 1,iba(ik)
                     rrr = rrr + Rphi_noncl(m,1,ib,i,is) &
	&                       *Rphi_noncl(m,1,ib,j,is) &
	&                      + Rphi_noncl(m,2,ib,i,is) &
	&                       *Rphi_noncl(m,2,ib,j,is)
                     rri = rri + Rphi_noncl(m,1,ib,i,is) &
	&                       *Rphi_noncl(m,2,ib,j,is) &
	&                      - Rphi_noncl(m,2,ib,i,is) &
	&                       *Rphi_noncl(m,1,ib,j,is)
                     pspr = pspr + phi_noncl(m,1,ib,i,is) &
	&                         *phi_noncl(m,1,ib,j,is) &
	&                        + phi_noncl(m,2,ib,i,is) &
	&                         *phi_noncl(m,2,ib,j,is)
                     pspi = pspi + phi_noncl(m,1,ib,i,is) &
	&                         *phi_noncl(m,2,ib,j,is) &
	&                        - phi_noncl(m,2,ib,i,is) &
	&                         *phi_noncl(m,1,ib,j,is)
                  end do
               End do
                rr(i,j,1) = rrr;               rr(i,j,2) = rri
               psp(i,j,1) = pspr;             psp(i,j,2) = pspi
            end if

            do ri = 1, kimg
               rr( j,i,ri) = pm(ri) *rr(i,j,ri)      ! pm=1,-1 (ri==1,2)
               psp(i,j,ri) = psp(i,j,ri) + sumqff(ri)
               psp(j,i,ri) = pm(ri) *psp(i,j,ri)
            end do
         end do
      end do

    end subroutine rr_and_psp_from_Rphi_etc_noncl

    subroutine rr_and_psp_from_rr_e_etc_noncl(ib_t,nrmm,rr,psp)
      integer, intent(in)       :: ib_t,nrmm
      real(DP),intent(out),dimension(0:nrmm-1,0:nrmm-1,kimg)::rr,psp
      real(DP), dimension(kimg) :: sumqff
      integer  ::                  j,i,m,ri, ib

      if(rmm_save_memory_mode == ON) then
         ib = 1
      else
         ib = ib_t
      end if

      rr = 0.d0; psp = 0.d0
      do j = 0, nrmm-1
         do i = 0,j
            call phsphn_noncl( i, j, ib, kimg, sumqff )
                      ! ->sumqff (= \sum_{n1,n2}q_{n1,n2}<W|b_{nl}><b_{n2}|W>
            do ri = 1, kimg
               rr(i,j,ri)  = rr_e(i,j,ib,ri)
               psp(i,j,ri) = psp_e(i,j,ib,ri)
               rr( j,i,ri) = pm(ri) * rr(i,j,ri)      ! pm=1,-1 (ri==1,2)
               psp(i,j,ri) = psp(i,j,ri) + sumqff(ri)
               psp(j,i,ri) = pm(ri) * psp(i,j,ri)
            end do
         end do
      end do

    end subroutine rr_and_psp_from_rr_e_etc_noncl

    subroutine phsphn_noncl(i,j,ib,nimg,sumqff)
      integer, intent(in)        :: i, j, ib, nimg  ! ib = 1 when save_memory_mode == ON
      real(kind=DP), intent(out) :: sumqff(nimg)
      integer         :: ia,m,n
      integer :: is1, is2, is_tmp
      real(kind=DP) :: c1, c2, c3, c4

      sumqff = 0.d0
      if (modnrm == EXECUT) then
         do ia = 1, nac
            m = nlmta1(ia); n = nlmta2(ia)
            Do is1=1, ndim_spinor
              Do is2=1, ndim_spinor
                is_tmp = ( is1 -1 )*ndim_spinor +is2

                c1 = real( fqwei_noncl(ia,is_tmp) ) &
                     &   *( bWr_noncl(m,i,ib,is1) *bWr_noncl(n,j,ib,is2) &
                     &    + bWi_noncl(m,i,ib,is1) *bWi_noncl(n,j,ib,is2) )
                c2 =-aimag( fqwei_noncl(ia,is_tmp) ) &
                     &   *( bWr_noncl(m,i,ib,is1) *bWi_noncl(n,j,ib,is2) &
                     &    + bWi_noncl(m,i,ib,is1) *bWr_noncl(n,j,ib,is2) )

                sumqff(1) = sumqff(1) + c1 + c2

              End do
           End do
         end do
         if (nimg == 2) then
            do ia = 1, nac
               m = nlmta1(ia); n = nlmta2(ia)
               Do is1=1, ndim_spinor
                 Do is2=1, ndim_spinor
                   is_tmp = ( is1 -1 )*ndim_spinor +is2

                   c3 = real( fqwei_noncl(ia,is_tmp) ) &
                        &   *( bWr_noncl(m,i,ib,is1) *bWi_noncl(n,j,ib,is2) &
                        &    - bWi_noncl(m,i,ib,is1) *bWr_noncl(n,j,ib,is2) )
                   c4 = aimag( fqwei_noncl(ia,is_tmp) ) &
                        &   *( bWr_noncl(m,i,ib,is1) *bWr_noncl(n,j,ib,is2) &
                        &    + bWi_noncl(m,i,ib,is1) *bWi_noncl(n,j,ib,is2) )

                   sumqff(2)=sumqff(2) + c3 + c4

                 End do
              End do
            end do
         end if
      end if
    end subroutine phsphn_noncl

    subroutine normalization_of_rr_and_psp(nrmm,rr,psp)
      integer, intent(in) :: nrmm
      real(DP),intent(out),dimension(0:nrmm-1,0:nrmm-1,kimg)::rr,psp
      real(DP) :: divrr,divpsp

      divrr = 1.d0/rr(0,0,1); divpsp = 1.d0/psp(0,0,1)
      rr    =  rr*divrr
      psp   = psp*divpsp

      if(iprirmm >= 2) call wd_rr_and_psp(nrmm,rr,psp)
    end subroutine normalization_of_rr_and_psp

    subroutine pre_calc_phase
      call mltpha4( natm, rmm_precal_phase_matm, pos, PAI2, kgp, ngabc, &
	&           kg1, nbase(1,ik), iba(ik), phasec, phases )
      if(rmm_precal_phase_matm < natm) then
         call setglist4(n_min(1),n_max(1),n_min(2),n_max(2) &  ! -(b_E.S.)
              &, n_min(3),n_max(3),nbase(1,ik),ngabc,kgp,iba(ik),nglist)
         call crngabc4(n_min(1),n_max(1),n_min(2),n_max(2) &   ! -(b_E.S.)
              &, n_min(3),n_max(3),nglist,kg1,nngabc,newp)
      end if

      if(iprirmm >= 2) then
         if(allocated(nglist)) then
            write(nfout,'(" !!rmm nglist is allocated <<pre_calc_phase>>")')
         else
            write(nfout,'(" !!rmm nglist is not allocated <<pre_calc_phase>>")')
         end if
      end if

    end subroutine pre_calc_phase

    subroutine zf_listing
      integer       :: ia
      real(kind=DP) :: f(3)

      integer :: id_sname = -1
      call tstatc0_begin('zf_listing ',id_sname)
      do ia = rmm_precal_phase_matm+1, natm
         f = pos(ia,1:3)*PAI2
         call zf_list_s(n_min(1),n_max(1),rmm_precal_phase_matm,natm,f(1),ia,zfc1,zfs1)
         call zf_list_s(n_min(2),n_max(2),rmm_precal_phase_matm,natm,f(2),ia,zfc2,zfs2)
         call zf_list_s(n_min(3),n_max(3),rmm_precal_phase_matm,natm,f(3),ia,zfc3,zfs3)
         !                 -(b_Electronic_Structure)
      end do
      call tstatc0_end(id_sname)

    end subroutine zf_listing
  end subroutine m_ESrmm_renew_WF_noncl
! =================================================================== 11.0

  subroutine Vnonlocal_W_RMM(ik,ibo,ispin,matm,iter_rmm)
    integer, intent(in) :: ik,ibo,ispin,matm,iter_rmm
    integer :: mdvdb, it, ia, lmt2, lmta2, ib
    real(kind=DP), allocatable,dimension(:)   ::   ar, ai             ! d(kg1)
    real(kind=DP), allocatable,dimension(:)   ::   zfcos, zfsin       ! d(nbmx)
    real(kind=DP), allocatable,target,dimension(:)   :: sc,ss,qc,qs   ! d(kg1)

    integer :: id_sname = -1
    call tstatc0_begin('Vnonlocal_W_RMM ',id_sname)

    ib = map_z(ibo)        ! MPI
    vnlph_l(:,ib,:) = 0.d0
    call alloc_zfsincos_arai()
    call alloc_scssqcqs()

    Loop_ntyp: do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       Loop_natm : do ia = 1, natm
          if(ityp(ia) /= it) cycle
          call calc_phase_RMM(ik,ia,matm)  ! (zfc1,2,3,zfs1,2,3)-> zfcos,zfsin
          call sumset_rmm(ik,ib,iter_rmm,it,ia) ! zfcos,zfsin,sc,ss,zaj,snl -> bWr,bWi
          do lmt2 = 1, ilmt(it)
             lmta2 = lmta(lmt2,ia)
             call Vnonlocal_W_part_sum_over_lmt1(ispin,ik,it,ia,lmt2,mdvdb)
             if(mdvdb == SKIP) then
                call add_vnlph_l_without_eko_part1(ik,ib,lmta2,iter_rmm)
             else if(mdvdb == EXECUT) then
                call add_vnlph_l_with_eko_part1(ik,ib,lmta2,iter_rmm)
             endif
          end do
       end do Loop_natm
    end do Loop_ntyp
    call dealloc_scssqcqs()
    call dealloc_zfsincos_arai()
    call tstatc0_end(id_sname)
  contains
    subroutine alloc_zfsincos_arai()
      allocate(zfsin(nbmx)); allocate(zfcos(nbmx))
      allocate(ar(nbmx)); allocate(ai(nbmx))
    end subroutine alloc_zfsincos_arai

    subroutine dealloc_zfsincos_arai()
      deallocate(ai,ar,zfcos,zfsin)
    end subroutine dealloc_zfsincos_arai

    subroutine alloc_scssqcqs()
      allocate(sc(kg1))
      allocate(ss(kg1))
      allocate(qc(kg1))
      allocate(qs(kg1))
    end subroutine alloc_scssqcqs

    subroutine dealloc_scssqcqs()
      deallocate(qs,qc,ss,sc)
    end subroutine dealloc_scssqcqs

    subroutine calc_phase_RMM(ik,ia,matm)
      integer, intent(in) :: ik,ia,matm

      integer  :: i, nb, n1,n2,n3
      real(DP) :: z12r,z12i
      integer :: id_sname = -1
      call tstatc0_begin('calc_phase_RMM ',id_sname)

      if(ia <= matm) then
         do i = 1, iba(ik)
            zfcos(i) = phasec(i,ia)
            zfsin(i) = phases(i,ia)
         end do
      else
         do i = 1, iba(ik)
            nb = newp(i)
            n1=nngabc(i,1); n2=nngabc(i,2); n3=nngabc(i,3)
            z12r = zfc1(n1,ia)*zfc2(n2,ia)-zfs1(n1,ia)*zfs2(n2,ia)
            z12i = zfs1(n1,ia)*zfc2(n2,ia)+zfc1(n1,ia)*zfs2(n2,ia)
            zfcos(nb) = zfc3(n3,ia)*z12r - zfs3(n3,ia)*z12i
            zfsin(nb) = zfs3(n3,ia)*z12r + zfc3(n3,ia)*z12i
         end do
      end if
      call tstatc0_end(id_sname)
    end subroutine calc_phase_RMM

    subroutine Vnonlocal_W_part_sum_over_lmt1(ispin,ik,it,ia,lmt2,mdvdb)
      integer, intent(in) :: ispin,ik,it,ia,lmt2,mdvdb
      integer             :: lmt1, lmtt1, il1, im1, il11, mdl, iksnl,il2,im2, ii, i
      real(kind=DP)       :: tmp

      integer :: id_sname = -1
      call tstatc0_begin('Vnonlocal_W_part_sum_over_lmt1 ',id_sname)

      iksnl = (ik-1)/nspin + 1
      il2   = ltp(lmt2,it)
      im2   = mtp(lmt2,it)

      sc = 0.d0; ss = 0.d0
      if(mdvdb == EXECUT) then
         qc = 0.d0; qs = 0.d0
      endif

      do lmt1 = 1,ilmt(it)
         lmtt1 = lmtt(lmt1,it)
         il1   = ltp(lmt1,it)
         im1   = mtp(lmt1,it)
         il11  = il1 - 1
         mdl   = mod(il11,4)
         if(il1 == il2 .and. im1 == im2) then
!!$            tmp = dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)
          if(ipaw(it)==0) then
              tmp = dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)
          else
              tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
          end if
         else
!!$            tmp = vlhxcQ(lmt1,lmt2,ia,ispin)
          if(ipaw(it)==0) then
              tmp = vlhxcQ(lmt1,lmt2,ia,ispin)
          else
              tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
          end if
         endif
         tmp = tmp * iwei(ia)
         if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
         if(mdl == 0 .or. mdl == 2) then
            do ii = 1, iba(ik)
               sc(ii) = sc(ii) + tmp*zfcos(ii)*snl(ii,lmtt1,iksnl)
               ss(ii) = ss(ii) - tmp*zfsin(ii)*snl(ii,lmtt1,iksnl)
            end do
         else if(mdl == 1 .or. mdl == 3) then
            do ii = 1, iba(ik)
               sc(ii) = sc(ii) - tmp*zfsin(ii)*snl(ii,lmtt1,iksnl)
               ss(ii) = ss(ii) - tmp*zfcos(ii)*snl(ii,lmtt1,iksnl)
            end do
         end if
         if(mdvdb == EXECUT .and. il1 == il2 .and. im1 == im2) then
            tmp = q(lmt1,lmt2,it)*iwei(ia)
            if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
            if(mdl == 0 .or. mdl == 2) then
               do ii = 1, iba(ik)
                  qc(ii) = qc(ii) + tmp*zfcos(ii)*snl(ii,lmtt1,iksnl)
                  qs(ii) = qs(ii) - tmp*zfsin(ii)*snl(ii,lmtt1,iksnl)
               end do
            else if(mdl == 1 .or. mdl == 3) then
               do ii = 1, iba(ik)
                  qc(ii) = qc(ii) - tmp*zfsin(ii)*snl(ii,lmtt1,iksnl)
                  qs(ii) = qs(ii) - tmp*zfcos(ii)*snl(ii,lmtt1,iksnl)
               end do
            end if
         end if
      end do
      call tstatc0_end(id_sname)
    end subroutine Vnonlocal_W_part_sum_over_lmt1

    subroutine sumset_rmm(ik,ib,iter_rmm,it,ia)
      integer, intent(in) :: ik,ib,iter_rmm,it,ia

      integer  :: lmt1, lmtt1, lmta1, il1, i, iksnl, il2
      real(DP) :: temp, fsrt,fsit,ar,ai,crt1,cit1

      iksnl = (ik-1)/nspin + 1
      do lmt1 = 1, ilmt(it)
         lmtt1 = lmtt(lmt1,it)
         lmta1 = lmta(lmt1,ia)
         il1   = ltp(lmt1,it)
         fsrt = 0.d0; fsit = 0.d0
         if(kimg == 1) then
            do i = 1, iba(ik)
               temp = snl(i,lmtt1,iksnl)*zaj_l(i,ib,ik,1)
               fsrt = fsrt + zfcos(i)*temp
               fsit = fsit + zfsin(i)*temp
            end do
         else if(kimg == 2) then
            if(k_symmetry(ik) == GAMMA) then
               il2 = mod(il1,2) ! il2 = 1 -> l=0,2,4,..,: il2 = 0 -> l=1,3,5,7,...
               if(il2 == 1) then
                  do i = 2, iba(ik)
                     ar = zfcos(i)*snl(i,lmtt1,iksnl)
                     ai = zfsin(i)*snl(i,lmtt1,iksnl)
                     crt1 = zaj_l(i,ib,ik,1)
                     cit1 = zaj_l(i,ib,ik,2)
                     fsrt = fsrt + ar*crt1 - ai*cit1
                  end do
                  fsrt = 2.d0*fsrt + zaj_l(1,ib,ik,1)*snl(1,lmtt1,iksnl)
               else
                  do i = 2, iba(ik)
                     ar = zfcos(i)*snl(i,lmtt1,iksnl)
                     ai = zfsin(i)*snl(i,lmtt1,iksnl)
                     crt1 = zaj_l(i,ib,ik,1)
                     cit1 = zaj_l(i,ib,ik,2)
                     fsit = fsit + ai*crt1 + ar*cit1
                  end do
                  fsit = 2.d0*fsit
               end if
            else
               do i    = 1,iba(ik)
                  ar = zfcos(i)*snl(i,lmtt1,iksnl)
                  ai = zfsin(i)*snl(i,lmtt1,iksnl)
                  crt1     = zaj_l(i,ib,ik,1)
                  cit1     = zaj_l(i,ib,ik,2)
                  fsrt = fsrt + ar*crt1-ai*cit1
                  fsit = fsit + ai*crt1+ar*cit1
               end do
            end if
         end if

         i = mod(il1,4)
         if(i == 2) then
            temp =  fsit;         fsit =  fsrt;         fsrt = -temp
         else if(i == 3) then
            fsrt = -fsrt;         fsit = -fsit;
         else if(i == 0) then
            temp =  fsit;         fsit = -fsrt;         fsrt = temp
         end if
         bWr(lmta1,iter_rmm,1) = fsrt; bWi(lmta1,iter_rmm,1) = fsit
      end do
    end subroutine sumset_rmm

    subroutine add_vnlph_l_without_eko_part1(ik,ib,lmta2,iter_rmm)
      integer, intent(in) :: ik,ib,lmta2,iter_rmm
      real(kind=DP) :: bWFr, bWFi
      integer       :: i

      bWFr = bWr(lmta2,iter_rmm,1); bWFi = bWi(lmta2,iter_rmm,1)

      if(kimg == 1) then
         do i = 1, iba(ik)
            vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*sc(i)-bWFi*ss(i)
         enddo
      else if(kimg == 2) then
         do i = 1, iba(ik)
            vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*sc(i)-bWFi*ss(i)
            vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi*sc(i)+bWFr*ss(i)
         enddo
      endif
    end subroutine add_vnlph_l_without_eko_part1

    subroutine add_vnlph_l_with_eko_part1(ik,ib,lmta2,iter_rmm)
      integer, intent(in) :: ik,ib,lmta2,iter_rmm
      real(kind=DP) :: bWFr, bWFi, e
      integer       :: i

      integer :: id_sname = -1
      call tstatc0_begin('add_vnlph_l_with_eko_part1 ',id_sname)

      bWFr = bWr(lmta2,iter_rmm,1); bWFi = bWi(lmta2,iter_rmm,1); e = eko_l(ib,ik)

      if(kimg == 1) then
         do i = 1, iba(ik)
            vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*(sc(i)-e*qc(i))-bWFi*(ss(i)-e*qs(i))
         enddo
      else if(kimg == 2) then
         do i = 1, iba(ik)
            vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*(sc(i)-e*qc(i))-bWFi*(ss(i)-e*qs(i))
            vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi*(sc(i)-e*qc(i))+bWFr*(ss(i)-e*qs(i))
         enddo
      endif
      call tstatc0_end(id_sname)
    end subroutine add_vnlph_l_with_eko_part1
  end subroutine Vnonlocal_W_RMM

! ======================================= added by K. Tagami ================ 11.0
  subroutine Vnonlocal_W_RMM_noncl( ik, ibo, ispin, matm, iter_rmm )
    integer, intent(in) :: ik,ibo,ispin,matm,iter_rmm
    integer :: mdvdb, it, ia, lmt2, lmta2, ib
    real(kind=DP), allocatable,dimension(:)   ::   ar, ai             ! d(kg1)
    real(kind=DP), allocatable,dimension(:)   ::   zfcos, zfsin       ! d(nbmx)
    real(kind=DP), allocatable,target,dimension(:)   :: sc,ss,qc,qs   ! d(kg1)

    integer :: id_sname = -1
    call tstatc0_begin('Vnonlocal_W_RMM_noncl ',id_sname)

    ib = map_z(ibo)        ! MPI
    vnlph_l(:,ib,:) = 0.d0

    call alloc_zfsincos_arai()
    call alloc_scssqcqs()

    Loop_ntyp: do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       Loop_natm : do ia = 1, natm
          if(ityp(ia) /= it) cycle
          call calc_phase_RMM(ik,ia,matm)  ! (zfc1,2,3,zfs1,2,3)-> zfcos,zfsin

          call sumset_rmm_noncl( ik, ib, iter_rmm, it, ia )
                             ! zfcos,zfsin,sc,ss,zaj,snl -> bWr_noncl,bWi_noncl

          do lmt2 = 1, ilmt(it)
             lmta2 = lmta(lmt2,ia)

             call Vnonlocal_W_part_sum_lmt1_noncl( ispin,ik,it,ia,lmt2,mdvdb )

             if(mdvdb == SKIP) then
                call add_vnlph_l_no_eko_pt1_noncl(ik,ib,lmta2,iter_rmm)
             else if(mdvdb == EXECUT) then
                call add_vnlph_l_with_eko_pt1_noncl(ik,ib,lmta2,iter_rmm)
             endif
          end do
       end do Loop_natm
    end do Loop_ntyp
    call dealloc_scssqcqs()
    call dealloc_zfsincos_arai()
    call tstatc0_end(id_sname)
  contains
    subroutine alloc_zfsincos_arai()
      allocate(zfsin(nbmx)); allocate(zfcos(nbmx))
      allocate(ar(nbmx)); allocate(ai(nbmx))
    end subroutine alloc_zfsincos_arai

    subroutine dealloc_zfsincos_arai()
      deallocate(ai,ar,zfcos,zfsin)
    end subroutine dealloc_zfsincos_arai

    subroutine alloc_scssqcqs()
      allocate(sc(kg1))
      allocate(ss(kg1))
      allocate(qc(kg1))
      allocate(qs(kg1))
    end subroutine alloc_scssqcqs

    subroutine dealloc_scssqcqs()
      deallocate(qs,qc,ss,sc)
    end subroutine dealloc_scssqcqs

    subroutine calc_phase_RMM(ik,ia,matm)
      integer, intent(in) :: ik,ia,matm

      integer  :: i, nb, n1,n2,n3
      real(DP) :: z12r,z12i
      integer :: id_sname = -1
      call tstatc0_begin('calc_phase_RMM ',id_sname)

      if(ia <= matm) then
         do i = 1, iba(ik)
            zfcos(i) = phasec(i,ia)
            zfsin(i) = phases(i,ia)
         end do
      else
         do i = 1, iba(ik)
            nb = newp(i)
            n1=nngabc(i,1); n2=nngabc(i,2); n3=nngabc(i,3)
            z12r = zfc1(n1,ia)*zfc2(n2,ia)-zfs1(n1,ia)*zfs2(n2,ia)
            z12i = zfs1(n1,ia)*zfc2(n2,ia)+zfc1(n1,ia)*zfs2(n2,ia)
            zfcos(nb) = zfc3(n3,ia)*z12r - zfs3(n3,ia)*z12i
            zfsin(nb) = zfs3(n3,ia)*z12r + zfc3(n3,ia)*z12i
         end do
      end if
      call tstatc0_end(id_sname)
    end subroutine calc_phase_RMM

    subroutine Vnonlocal_W_part_sum_lmt1_noncl(ispin,ik,it,ia,lmt2,mdvdb)
      integer, intent(in) :: ispin,ik,it,ia,lmt2,mdvdb
      integer             :: lmt1, lmtt1, il1, im1, il11, mdl, iksnl,il2,im2, ii, i
      complex(kind=CMPLDP)       :: tmp
      real(kind=DP) :: c1, c2

      integer :: id_sname = -1
      call tstatc0_begin('Vnonlocal_W_part_sum_lmt1_noncl ',id_sname)

      iksnl = (ik-1)/ndim_spinor + 1
      il2   = ltp(lmt2,it)
      im2   = mtp(lmt2,it)

      sc = 0.d0; ss = 0.d0
      if(mdvdb == EXECUT) then
         qc = 0.d0; qs = 0.d0
      endif

      do lmt1 = 1,ilmt(it)
         lmtt1 = lmtt(lmt1,it)
         il1   = ltp(lmt1,it)
         im1   = mtp(lmt1,it)
         il11  = il1 - 1
         mdl   = mod(il11,4)

! ---------------------------------------------------- 11.0S
#ifdef SKIP_TEST
         if ( mdvdb == SKIP ) then
            if ( il1 /= il2 ) cycle
            if ( SpinOrbit_mode == Neglected .and. sw_hubbard == OFF ) then
               if ( im1 /= im2 ) cycle
            endif
         endif
#endif
! ---------------------------------------------------- 11.0S

         tmp = dion_scr_noncl( lmt1,lmt2,ispin,ia )
         tmp = tmp * iwei(ia)

         if (mdl == 2 .or. mdl == 3) tmp = -1*tmp

         if (mdl == 0 .or. mdl == 2) then
            c1 = real(tmp);  c2 = aimag(tmp)
            do ii = 1, iba(ik)
               sc(ii) = sc(ii) + ( c1*zfcos(ii)+c2*zfsin(ii) )*snl(ii,lmtt1,iksnl)
               ss(ii) = ss(ii) + (-c1*zfsin(ii)+c2*zfcos(ii) )*snl(ii,lmtt1,iksnl)
            end do
         else if(mdl == 1 .or. mdl == 3) then
            c1 = real(tmp);  c2 = aimag(tmp)
            do ii = 1, iba(ik)
               sc(ii) = sc(ii) + (-c1*zfsin(ii)+c2*zfcos(ii) )*snl(ii,lmtt1,iksnl)
               ss(ii) = ss(ii) + (-c1*zfcos(ii)-c2*zfsin(ii) )*snl(ii,lmtt1,iksnl)
            end do
         end if

! ----------------------------------------------------- 11.0S
!          if (mdvdb == EXECUT .and. il1 == il2 .and. im1 == im2) then
         if ( mdvdb == EXECUT .and. il1 == il2 ) then
            if ( SpinOrbit_mode /= BuiltIn ) then
               if ( im1 /= im2 ) cycle
            endif
! ------------------------------------------------------ 11.0S

            tmp = q_noncl(lmt1,lmt2,ispin,it)*iwei(ia)

            if (mdl == 2 .or. mdl == 3) tmp = -1*tmp
            if (mdl == 0 .or. mdl == 2) then
               c1 = real(tmp);  c2 = aimag(tmp)

               do ii = 1, iba(ik)
                  qc(ii) = qc(ii) + ( c1*zfcos(ii)+c2*zfsin(ii) )*snl(ii,lmtt1,iksnl)
                  qs(ii) = qs(ii) + (-c1*zfsin(ii)+c2*zfcos(ii) )*snl(ii,lmtt1,iksnl)
               end do
            else if(mdl == 1 .or. mdl == 3) then
               c1 = real(tmp);  c2 = aimag(tmp)

               do ii = 1, iba(ik)
                  qc(ii) = qc(ii) + (-c1*zfsin(ii)+c2*zfcos(ii) )*snl(ii,lmtt1,iksnl)
                  qs(ii) = qs(ii) + (-c1*zfcos(ii)-c2*zfsin(ii) )*snl(ii,lmtt1,iksnl)
               end do
            end if
         end if
      end do
      call tstatc0_end(id_sname)
    end subroutine Vnonlocal_W_part_sum_lmt1_noncl

    subroutine sumset_rmm_noncl( ik, ib, iter_rmm, it, ia )
      integer, intent(in) :: ik, ib,iter_rmm,it,ia

      integer  :: lmt1, lmtt1, lmta1, il1, i, iksnl, il2
      real(DP) :: temp, fsrt,fsit,ar,ai,crt1,cit1

      integer :: my_is, is,  k1

      iksnl = (ik-1)/ndim_spinor + 1
      my_is = mod( ik-1, ndim_spinor ) +1

!!!      k1 =  ik + is -1

      do lmt1 = 1, ilmt(it)
         lmtt1 = lmtt(lmt1,it)
         lmta1 = lmta(lmt1,ia)
         il1   = ltp(lmt1,it)
         fsrt = 0.d0; fsit = 0.d0
         if (kimg == 1) then
            do i = 1, iba(ik)
               temp = snl(i,lmtt1,iksnl) *zaj_l( i, ib, ik, 1 )
               fsrt = fsrt + zfcos(i)*temp
               fsit = fsit + zfsin(i)*temp
            end do
         else if (kimg == 2) then
            if (k_symmetry(ik) == GAMMA) then
               il2 = mod(il1,2) ! il2 = 1 -> l=0,2,4,..,: il2 = 0 -> l=1,3,5,7,...
               if (il2 == 1) then
                  do i = 2, iba(ik)
                     ar = zfcos(i)*snl(i,lmtt1,iksnl)
                     ai = zfsin(i)*snl(i,lmtt1,iksnl)
                     crt1 = zaj_l( i, ib, ik, 1 )
                     cit1 = zaj_l( i, ib, ik, 2 )
                     fsrt = fsrt + ar*crt1 - ai*cit1
                  end do
                  temp = zaj_l( 1, ib, ik,1 )*snl(1,lmtt1,iksnl)
                  fsrt = 2.d0*fsrt + temp

               else
                  do i = 2, iba(ik)
                     ar = zfcos(i)*snl(i,lmtt1,iksnl)
                     ai = zfsin(i)*snl(i,lmtt1,iksnl)
                     crt1 = zaj_l( i, ib, ik, 1 )
                     cit1 = zaj_l( i, ib, ik, 2 )
                     fsit = fsit + ai*crt1 + ar*cit1
                  end do
                  fsit = 2.d0*fsit
               end if
            else
               do i = 1,iba(ik)
                  ar = zfcos(i)*snl(i,lmtt1,iksnl)
                  ai = zfsin(i)*snl(i,lmtt1,iksnl)
                  crt1 = zaj_l( i, ib, ik, 1 )
                  cit1 = zaj_l( i, ib, ik, 2 )
                  fsrt = fsrt + ar*crt1-ai*cit1
                  fsit = fsit + ai*crt1+ar*cit1
               end do
            end if
         end if

         i = mod(il1,4)
         if(i == 2) then
            temp =  fsit;         fsit =  fsrt;         fsrt = -temp
         else if(i == 3) then
            fsrt = -fsrt;         fsit = -fsit;
         else if(i == 0) then
            temp =  fsit;         fsit = -fsrt;         fsrt = temp
         end if
         bWr_noncl(lmta1,iter_rmm,1,my_is) = fsrt
	 bWi_noncl(lmta1,iter_rmm,1,my_is) = fsit

      end do
    end subroutine sumset_rmm_noncl

    subroutine  add_vnlph_l_no_eko_pt1_noncl(ik,ib,lmta2,iter_rmm)
      integer, intent(in) :: ik,ib,lmta2,iter_rmm
      real(kind=DP) :: bWFr, bWFi
      integer       :: i
      integer :: my_is

      my_is = mod( ik-1, ndim_spinor ) +1
      bWFr = bWr_noncl( lmta2,iter_rmm,1,my_is )
      bWFi = bWi_noncl( lmta2,iter_rmm,1,my_is )

      if(kimg == 1) then
         do i = 1, iba(ik)
            vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*sc(i)-bWFi*ss(i)
         enddo
      else if(kimg == 2) then
         do i = 1, iba(ik)
            vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*sc(i)-bWFi*ss(i)
            vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi*sc(i)+bWFr*ss(i)
         enddo
      endif
    end subroutine add_vnlph_l_no_eko_pt1_noncl

    subroutine add_vnlph_l_with_eko_pt1_noncl(ik,ib,lmta2,iter_rmm)
      integer, intent(in) :: ik,ib,lmta2,iter_rmm
      real(kind=DP) :: bWFr, bWFi, e
      integer       :: i
      integer :: iksnl, my_is, ik_for_pointing_eko

      integer :: id_sname = -1
      call tstatc0_begin('add_vnlph_l_with_eko_part1 ',id_sname)

      iksnl = (ik-1)/ndim_spinor + 1
      my_is = mod( ik-1, ndim_spinor ) +1
      ik_for_pointing_eko = ( iksnl -1 )*ndim_spinor +1

      bWFr = bWr_noncl(lmta2,iter_rmm,1,my_is)
      bWFi = bWi_noncl(lmta2,iter_rmm,1,my_is)

      e = eko_l(ib,ik_for_pointing_eko)

      if(kimg == 1) then
         do i = 1, iba(ik)
            vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*(sc(i)-e*qc(i))-bWFi*(ss(i)-e*qs(i))
         enddo
      else if(kimg == 2) then
         do i = 1, iba(ik)
            vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*(sc(i)-e*qc(i))-bWFi*(ss(i)-e*qs(i))
            vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi*(sc(i)-e*qc(i))+bWFr*(ss(i)-e*qs(i))
         enddo
      endif
      call tstatc0_end(id_sname)
    end subroutine add_vnlph_l_with_eko_pt1_noncl

  end subroutine Vnonlocal_W_RMM_noncl
! ======================================================================== 11.0

  subroutine Vnonlocal_W_RMMn(nfout,ik,ispin,iter_rmm)
    integer, intent(in) :: nfout,ik,ispin,iter_rmm
    integer :: mdvdb, it, ia, lmt2, lmta2, ia_p,iap,lmta1
    integer :: n_ialist, n_ialist0, ia_start, n_iagroup, n_ia, ia_g, max_n_ialist0 ,nbmx_adj, kg1_adj
    real(kind=DP), allocatable,dimension(:,:) ::  ar_x, ai_x                 ! d(kg1_adj,n_ialist0)
    real(kind=DP), allocatable,dimension(:,:) ::  zfcos_x, zfsin_x           ! d(nbmx,n_ialist0)
    real(kind=DP), allocatable,target,dimension(:,:) :: sc_x,ss_x,qc_x,qs_x  ! d(kg1+,n_ialist))
    real(kind=DP), allocatable,target,dimension(:,:) :: sc_lmt,ss_lmt,qc_lmt,qs_lmt  ! d(kg1+,n_ialist))
    integer, allocatable, dimension(:) :: ia_list
    real(kind=DP), allocatable, dimension(:,:) :: fq
#ifdef NONLOCAL_RMM_DGEMM
    real(kind=DP), allocatable, dimension(:,:) :: bWr_lmt, bWi_lmt
#endif
    integer :: id_sname = -1, id_sname2 = -1
    call tstatc0_begin('Vnonlocal_W_RMMn ',id_sname,1)

#ifdef NONLOCAL_RMM_DGEMM
    n_ialist = 16
#else
    n_ialist = 1
#ifdef HIUX
    n_ialist = 16
#endif
#ifdef VPP
    n_ialist = 1
#endif
#ifdef SX
    n_ialist = 8
#endif
#endif
    if(n_ialist <=0) call phase_error_with_msg(nfout, 'n_ialist is illegal <<Vnonlocal_W_RMMn>>',__LINE__,__FILE__)
    call find_max_n_ialist0(max_n_ialist0)
    if(iprirmm >= 2) write(nfout,'(" !Vnonlocal_W_RMMn: max_n_ialist0 = ",i8)') max_n_ialist0
    nbmx_adj = nbmx
    if(mod(nbmx_adj,2) == 0) nbmx_adj = nbmx_adj+1
    kg1_adj = kg1
    if(mod(kg1_adj,2) == 0) kg1_adj = kg1_adj+1

    allocate(zfsin_x(nbmx_adj,max_n_ialist0)); allocate(zfcos_x(nbmx_adj,max_n_ialist0))
    allocate(ar_x(kg1_adj,max_n_ialist0)); allocate(ai_x(kg1_adj,max_n_ialist0))
    allocate(ia_list(n_ialist)); ia_list = 0
    allocate(fq(n_ialist,2)); fq = 0.d0

    vnlph_l(:,:,:) = 0.d0
    Loop_ntyp: do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
#ifdef NONLOCAL_RMM_DGEMM
       call alloc_scssqcqs_lmt(kg1_adj,max_n_ialist0*ilmt(it),mdvdb)
       call alloc_scssqcqs_x(kg1_adj,max_n_ialist0,mdvdb)
#else
       call alloc_scssqcqs_x(kg1_adj,max_n_ialist0,mdvdb)
#endif

       n_ia = 0
       do ia = 1, natm
          if(ityp(ia) == it) n_ia = n_ia + 1
       end do

       n_iagroup = n_ia/n_ialist + 1
       ia_start = 1
       if(iprirmm >= 2) write(nfout,'(" !Vnonlocal_W_RMMn: n_iagroup = ",i8, " ityp = ",i8)') n_iagroup,it
       Loop_ia_group: do ia_g = 1, n_iagroup
          n_ialist0 = 0
          ia_list = 0
          AtomcountLoop: do ia = ia_start, natm
             if(ityp(ia) == it) then
                n_ialist0 = n_ialist0 + 1
                ia_list(n_ialist0) = ia
             end if
              if(n_ialist0 >= n_ialist) exit AtomcountLoop
          end do AtomcountLoop
          ia_start = ia+1
          if(n_ialist0 >= 1 )then
             if(iprirmm >= 2) write(nfout,'(" !m_ES_WF_by_RMM.Vnonlocal_W_RMMn: ia_list = ",8i8)') (ia_list(ia),ia=1,n_ialist0)

                  call tstatc0_begin('calc_phasek_b ',id_sname2)
             call calc_phasek_b(natm,pos,n_ialist0,ia_list,kgp,iba(ik) &
                  & ,ngabc,kg1,nbase(1,ik),nbmx_adj,zfcos_x,zfsin_x) ! b_E.S.
                  call tstatc0_end(id_sname2)

#ifdef NONLOCAL_RMM_DGEMM
             call alloc_bW_lmt(it)  ! alloc bWr_lmt, bWi_lmp (np_e, n_ialist0*ilmt(it))
             call sumset_rmm_all4()                         ! zfcos,zfsin,sc,ss,zaj,snl -> bWr_tmp,bWi_tmp, bWr, bWi

             call Vnonlocal_W_part_sum_ovr_lmt4()           ! -> sc_x,ss_x,qc_x,qs_x
             call add_vnlph_l_part4(mdvdb)

             call dealloc_bW_lmt  ! alloc bWr_lmt, bWi_lmp (np_e, n_ialist0*ilmt(it))
#else
             call sumset_rmm_all3()                         ! zfcos,zfsin,sc,ss,zaj,snl -> bWr,bWi

             do lmt2 = 1, ilmt(it)
                call Vnonlocal_W_part_sum_over_lmt1b()      ! -> sc_x,ss_x,qc_x,qs_x
                if(mdvdb == SKIP) then
                   call add_vnlph_l_without_eko_part3()     ! -> vnlph_l
                else if(mdvdb == EXECUT) then
                   call add_vnlph_l_with_eko_part3()        ! -> vnlph_l
                endif
             end do
#endif
          end if
       end do Loop_ia_group
       call dealloc_scssqcqs_x()
#ifdef NONLOCAL_RMM_DGEMM
       call dealloc_scssqcqs_lmt()
#endif
    end do Loop_ntyp

    deallocate(fq)
    deallocate(ai_x,ar_x)
    deallocate(ia_list)
    deallocate(zfsin_x,zfcos_x)
    call tstatc0_end(id_sname)
  contains
    subroutine alloc_scssqcqs_lmt(kg1_adj,m,mdvdb)
      integer, intent(in) :: kg1_adj,m,mdvdb
      allocate(sc_lmt(kg1_adj,m))
      allocate(ss_lmt(kg1_adj,m))
      if(mdvdb == EXECUT) then
         allocate(qc_lmt(kg1_adj,m))
         allocate(qs_lmt(kg1_adj,m))
      end if
    end subroutine alloc_scssqcqs_lmt

    subroutine dealloc_scssqcqs_lmt()
      deallocate(sc_lmt,ss_lmt)
      if(mdvdb == EXECUT) deallocate(qc_lmt,qs_lmt)
    end subroutine dealloc_scssqcqs_lmt

    subroutine alloc_scssqcqs_x(kg1_adj,m,mdvdb)
      integer, intent(in) :: kg1_adj,m,mdvdb
      allocate(sc_x(kg1_adj,m))
      allocate(ss_x(kg1_adj,m))
      if(mdvdb == EXECUT) then
         allocate(qc_x(kg1_adj,m))
         allocate(qs_x(kg1_adj,m))
      end if
    end subroutine alloc_scssqcqs_x

    subroutine dealloc_scssqcqs_x()
      if(allocated(qs_x)) deallocate(qs_x)
      if(allocated(qc_x)) deallocate(qc_x)
      deallocate(ss_x,sc_x)
    end subroutine dealloc_scssqcqs_x

    subroutine find_max_n_ialist0(n)
      integer, intent(out) :: n

      n = 0
      do it = 1, ntyp

         n_ia = 0
         do ia = 1, natm
            if(ityp(ia) == it) n_ia = n_ia + 1
         end do

         n_iagroup = n_ia/n_ialist + 1
         ia_start = 1
         if(iprirmm >= 2) write(nfout,'(" !m_ES_WF_by_RMM.Vnonlocal_W_RMMn: n_iagroup = ",i8, " ityp = ",i8)') n_iagroup,it
         do ia_g = 1, n_iagroup
            n_ialist0 = 0
            AtomcountLoop: do ia = ia_start, natm
               if(ityp(ia) == it)  n_ialist0 = n_ialist0 + 1
               if(n_ialist0 >= n_ialist) exit AtomcountLoop
            end do AtomcountLoop
            if(n < n_ialist0) n = n_ialist0
            ia_start = ia+1
         end do
      end do

     end subroutine find_max_n_ialist0

#ifdef NONLOCAL_RMM_DGEMM
    subroutine alloc_bW_lmt(it)
      integer, intent(in) :: it
      integer :: nsize1
      nsize1 = n_ialist0 * ilmt(it)
      allocate(bWr_lmt(nsize1,np_e))
      allocate(bWi_lmt(nsize1,np_e))
    end subroutine alloc_bW_lmt

    subroutine dealloc_bW_lmt()
      deallocate(bWr_lmt)
      deallocate(bWi_lmt)
    end subroutine dealloc_bW_lmt

    subroutine sumset_rmm_all4()
      integer  :: lmt1, lmtt1, lmta1, i, iksnl, ib, il1, ia, iap, il2
      integer :: M, N, K
      real(kind=DP) :: fsrt, fsit
      real(kind=DP), allocatable, dimension(:,:) :: bWr_tmp, bWi_tmp
      integer :: id_sname = -1
      call tstatc0_begin('sumset_rmm_all4 ',id_sname)

      allocate(bWr_tmp(n_ialist0,np_e)); bWr_tmp = 0.d0
      allocate(bWi_tmp(n_ialist0,np_e)); bWi_tmp = 0.d0

      iksnl = (ik-1)/nspin + 1
      do lmt1 = 1, ilmt(it)
         lmtt1 = lmtt(lmt1,it)
         il1   = ltp(lmt1,it)
         do iap = 1, n_ialist0
            do i = 1, iba(ik)
               ar_x(i,iap) = zfcos_x(i,iap)*snl(i,lmtt1,iksnl)
               ai_x(i,iap) = zfsin_x(i,iap)*snl(i,lmtt1,iksnl)
            end do
         end do


         M = n_ialist0; N = np_e; K = iba(ik)
         if(kimg == 1) then
            call DGEMM__('T','N',M,N,K,1.d0,ar_x,kg1_adj,zaj_l(1,1,ik,1),kg1,0.d0,bWr_tmp,M)
            call DGEMM__('T','N',M,N,K,1.d0,ai_x,kg1_adj,zaj_l(1,1,ik,1),kg1,0.d0,bWi_tmp,M)

         else if(kimg == 2) then
            if(k_symmetry(ik) == GAMMA) then
               il2 = mod(il1,2) ! il2 = 1 -> l=0,2,4,6,...: il2 = 0 -> l=1,3,5,7,...
               if(il2 == 1) then
                  ai_x(1,1:n_ialist0) = 0.d0;  ar_x(1,1:n_ialist0) = ar_x(1,1:n_ialist0)*0.5d0
                  call DGEMM__('T','N', M,N,K, 1.d0,ar_x,kg1_adj,zaj_l(1,1,ik,1),kg1,0.d0,bWr_tmp,M)
                  call DGEMM__('T','N', M,N,K,-1.d0,ai_x,kg1_adj,zaj_l(1,1,ik,2),kg1,1.d0,bWr_tmp,M)
                  bWr_tmp = bWr_tmp*2.d0;   bWi_tmp = 0.d0
               else
                  ar_x(1,1:n_ialist0) = 0.d0; ai_x(1,1:n_ialist0) = 0.d0
                  call DGEMM__('T','N', M,N,K, 1.d0,ai_x,kg1_adj,zaj_l(1,1,ik,1),kg1,0.d0,bWi_tmp,M)
                  call DGEMM__('T','N', M,N,K, 1.d0,ar_x,kg1_adj,zaj_l(1,1,ik,2),kg1,1.d0,bWi_tmp,M)
                  bWr_tmp = 0.d0;           bWi_tmp = bWi_tmp*2.d0
               end if
            else
               call DGEMM__('T','N', M,N,K, 1.d0,ar_x,kg1_adj,zaj_l(1,1,ik,1),kg1,0.d0,bWr_tmp,M)
               call DGEMM__('T','N', M,N,K,-1.d0,ai_x,kg1_adj,zaj_l(1,1,ik,2),kg1,1.d0,bWr_tmp,M)
               call DGEMM__('T','N', M,N,K, 1.d0,ai_x,kg1_adj,zaj_l(1,1,ik,1),kg1,0.d0,bWi_tmp,M)
               call DGEMM__('T','N', M,N,K, 1.d0,ar_x,kg1_adj,zaj_l(1,1,ik,2),kg1,1.d0,bWi_tmp,M)
            end if
         end if

         i = mod(il1,4)
         do iap = 1, n_ialist0
            if(i==2) then
               do ib = 1, np_e
                  fsrt = bWr_tmp(iap,ib); fsit = bWi_tmp(iap,ib)
                  bWr_tmp(iap,ib) = -fsit; bWi_tmp(iap,ib) = fsrt
               end do
            else if(i==3) then
               do ib = 1, np_e
                  fsrt = bWr_tmp(iap,ib); fsit = bWi_tmp(iap,ib)
                  bWr_tmp(iap,ib) = -fsrt; bWi_tmp(iap,ib) = -fsit
               end do
            else if(i==0) then
               do ib = 1, np_e
                  fsrt = bWr_tmp(iap,ib); fsit = bWi_tmp(iap,ib)
                  bWr_tmp(iap,ib) = fsit; bWi_tmp(iap,ib) = -fsrt
               end do
            end if
         end do

         do ib = 1, np_e
            do iap = 1, n_ialist0
               ia = n_ialist0*(lmt1-1) + iap
               bWr_lmt(ia,ib) = bWr_tmp(iap,ib)
               bWi_lmt(ia,ib) = bWi_tmp(iap,ib)
            end do
         end do

         do iap = 1, n_ialist0
            ia = ia_list(iap)
            lmta1 = lmta(lmt1,ia)
            do ib = 1, np_e
               bWr(lmta1,iter_rmm,ib) =  bWr_tmp(iap,ib)
               bWi(lmta1,iter_rmm,ib) =  bWi_tmp(iap,ib)
            end do
         end do

      end do

      deallocate(bWr_tmp, bWi_tmp)

      call tstatc0_end(id_sname)

    end subroutine sumset_rmm_all4

    subroutine Vnonlocal_W_part_sum_ovr_lmt4()
      integer             :: lmt1, lmtt1, il1, im1, il11, mdl, iksnl,il2,im2, ia_p, ia, i, lmt2,iap
      real(kind=DP)       :: tmp

      integer :: id_sname = -1
      call tstatc0_begin('Vnonlocal_W_part_sum_ovr_lmt4 ',id_sname)

      sc_lmt = 0.d0; ss_lmt = 0.d0
      if(mdvdb == EXECUT) then
         qc_lmt = 0.d0; qs_lmt = 0.d0
      endif

      iksnl = (ik-1)/nspin + 1
      do lmt2 = 1, ilmt(it)
         il2   = ltp(lmt2,it)
         im2   = mtp(lmt2,it)


         do lmt1 = 1,ilmt(it)
            lmtt1 = lmtt(lmt1,it)
            il1   = ltp(lmt1,it)
            im1   = mtp(lmt1,it)
            il11  = il1 - 1
            mdl   = mod(il11,4)

            do ia_p = 1, n_ialist0
               ia = ia_list(ia_p)
               if(il1 == il2 .and. im1 == im2) then
                  if(ipaw(it)==0) then
                     tmp = (dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)) * iwei(ia)
                  else
                     tmp = (dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)) * iwei(ia)
                  end if
               else
                  if(ipaw(it)==0) then
                     tmp = vlhxcQ(lmt1,lmt2,ia,ispin) * iwei(ia)
                  else
                     tmp = (dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)) * iwei(ia)
                  end if
               endif
               if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
               fq(ia_p,1) = tmp
            end do

            if(mdvdb == EXECUT .and. il1 == il2 .and. im1 == im2) then
               do ia_p = 1, n_ialist0
                  ia = ia_list(ia_p)
                  tmp = q(lmt1,lmt2,it)*iwei(ia)
                  if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
                  fq(ia_p,2) = tmp
               end do

               if(mdl == 0 .or. mdl == 2) then
                  do ia_p = 1, n_ialist0
                     iap = n_ialist0*(lmt2-1) + ia_p
                     do i = 1, iba(ik)
                        sc_lmt(i,iap) = sc_lmt(i,iap) + fq(ia_p,1)*zfcos_x(i,ia_p)*snl(i,lmtt1,iksnl)
                        ss_lmt(i,iap) = ss_lmt(i,iap) - fq(ia_p,1)*zfsin_x(i,ia_p)*snl(i,lmtt1,iksnl)
                        qc_lmt(i,iap) = qc_lmt(i,iap) + fq(ia_p,2)*zfcos_x(i,ia_p)*snl(i,lmtt1,iksnl)
                        qs_lmt(i,iap) = qs_lmt(i,iap) - fq(ia_p,2)*zfsin_x(i,ia_p)*snl(i,lmtt1,iksnl)
                     end do
                  end do
               else if(mdl == 1 .or. mdl == 3) then
                  do ia_p = 1, n_ialist0
                     iap = n_ialist0*(lmt2-1) + ia_p
                     do i = 1, iba(ik)
                        sc_lmt(i,iap) = sc_lmt(i,iap) - fq(ia_p,1)*zfsin_x(i,ia_p)*snl(i,lmtt1,iksnl)
                        ss_lmt(i,iap) = ss_lmt(i,iap) - fq(ia_p,1)*zfcos_x(i,ia_p)*snl(i,lmtt1,iksnl)
                        qc_lmt(i,iap) = qc_lmt(i,iap) - fq(ia_p,2)*zfsin_x(i,ia_p)*snl(i,lmtt1,iksnl)
                        qs_lmt(i,iap) = qs_lmt(i,iap) - fq(ia_p,2)*zfcos_x(i,ia_p)*snl(i,lmtt1,iksnl)
                     end do
                  end do
               end if
            else
               if(mdl == 0 .or. mdl == 2) then
                  do ia_p = 1, n_ialist0
                     iap = n_ialist0*(lmt2-1) + ia_p
                     do i = 1, iba(ik)
                        sc_lmt(i,iap) = sc_lmt(i,iap) + fq(ia_p,1)*zfcos_x(i,ia_p)*snl(i,lmtt1,iksnl)
                        ss_lmt(i,iap) = ss_lmt(i,iap) - fq(ia_p,1)*zfsin_x(i,ia_p)*snl(i,lmtt1,iksnl)
                     end do
                  end do
               else if(mdl == 1 .or. mdl == 3) then
                  do ia_p = 1, n_ialist0
                     iap = n_ialist0*(lmt2-1) + ia_p
                     do i = 1, iba(ik)
                        sc_lmt(i,iap) = sc_lmt(i,iap) - fq(ia_p,1)*zfsin_x(i,ia_p)*snl(i,lmtt1,iksnl)
                        ss_lmt(i,iap) = ss_lmt(i,iap) - fq(ia_p,1)*zfcos_x(i,ia_p)*snl(i,lmtt1,iksnl)
                     end do
                  end do
               end if
            end if
         end do
      end do

      call tstatc0_end(id_sname)
    end subroutine Vnonlocal_W_part_sum_ovr_lmt4

    subroutine add_vnlph_l_part4(mdvdb)
       integer, intent(in) :: mdvdb
       integer       :: ib, iap,  M, N, K
       real(kind=DP), allocatable, dimension(:,:) :: bWer, bWei

       integer :: id_sname = -1
       call tstatc0_begin('add_vnlph_l_part4 ',id_sname)

       if(mdvdb == EXECUT) then
          allocate(bWer(ilmt(it)*n_ialist0,np_e))
          allocate(bWei(ilmt(it)*n_ialist0,np_e))
          do iap = 1, ilmt(it)*n_ialist0
             do ib = 1, np_e
                bWer(iap,ib) = bWr_lmt(iap,ib)*eko_l(ib,ik)
                bWei(iap,ib) = bWi_lmt(iap,ib)*eko_l(ib,ik)
             end do
          end do
       end if

       M = iba(ik); N = np_e ; K = ilmt(it)*n_ialist0
       if(kimg == 1) then
          call DGEMM__('N','N',M,N,K,  1.d0, sc_lmt,kg1_adj, bWr_lmt,K, 1.d0, vnlph_l(1,1,1),kg1)
          call DGEMM__('N','N',M,N,K, -1.d0, ss_lmt,kg1_adj, bWi_lmt,K, 1.d0, vnlph_l(1,1,1),kg1)

          if(mdvdb == EXECUT) then
          call DGEMM__('N','N',M,N,K, -1.d0, qc_lmt,kg1_adj, bWer,   K, 1.d0, vnlph_l(1,1,1),kg1)
          call DGEMM__('N','N',M,N,K,  1.d0, qs_lmt,kg1_adj, bWei,   K, 1.d0, vnlph_l(1,1,1),kg1)
          end if
       else if(kimg == 2) then
          call DGEMM__('N','N',M,N,K,  1.d0, sc_lmt,kg1_adj, bWr_lmt,K, 1.d0, vnlph_l(1,1,1),kg1)
          call DGEMM__('N','N',M,N,K, -1.d0, ss_lmt,kg1_adj, bWi_lmt,K, 1.d0, vnlph_l(1,1,1),kg1)
          call DGEMM__('N','N',M,N,K,  1.d0, sc_lmt,kg1_adj, bWi_lmt,K, 1.d0, vnlph_l(1,1,2),kg1)
          call DGEMM__('N','N',M,N,K,  1.d0, ss_lmt,kg1_adj, bWr_lmt,K, 1.d0, vnlph_l(1,1,2),kg1)

          if(mdvdb == EXECUT) then
          call DGEMM__('N','N',M,N,K, -1.d0, qc_lmt,kg1_adj, bWer,   K, 1.d0, vnlph_l(1,1,1),kg1)
          call DGEMM__('N','N',M,N,K,  1.d0, qs_lmt,kg1_adj, bWei,   K, 1.d0, vnlph_l(1,1,1),kg1)
          call DGEMM__('N','N',M,N,K, -1.d0, qc_lmt,kg1_adj, bWei,   K, 1.d0, vnlph_l(1,1,2),kg1)
          call DGEMM__('N','N',M,N,K, -1.d0, qs_lmt,kg1_adj, bWer,   K, 1.d0, vnlph_l(1,1,2),kg1)
          end if
       end if

       if(mdvdb == EXECUT) deallocate(bWer,bWei)
       call tstatc0_end(id_sname)
    end subroutine add_vnlph_l_part4
#else

     subroutine sumset_rmm_all3()

       integer  :: lmt1, lmtt1, lmta1, i, iksnl, ib, il1, ia, iap, il2
       real(kind=DP) :: fsrt,fsit
       integer :: id_sname = -1
       call tstatc0_begin('sumset_rmm_all3 ',id_sname)

       iksnl = (ik-1)/nspin + 1
       do lmt1 = 1, ilmt(it)
          lmtt1 = lmtt(lmt1,it)
          il1   = ltp(lmt1,it)
          if(kimg == 1) then

#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do iap = 1, n_ialist0
                ia = ia_list(iap)
!               lmta1 = lmta(lmt1,ia)
                do i = 1, iba(ik)
                   ar_x(i,iap) = zfcos_x(i,iap)*snl(i,lmtt1,iksnl)
                   ai_x(i,iap) = zfsin_x(i,iap)*snl(i,lmtt1,iksnl)
                end do
             end do
          else if(kimg == 2) then
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do iap = 1,n_ialist0
                ia = ia_list(iap)
!               lmta1 = lmta(lmt1,ia)
                do i = 1, iba(ik)
                   ar_x(i,iap) = zfcos_x(i,iap)*snl(i,lmtt1,iksnl)
                   ai_x(i,iap) = zfsin_x(i,iap)*snl(i,lmtt1,iksnl)
                end do
             end do
          end if

          if(kimg == 1) then
#ifdef HIUX
*poption indep
*poption parallel
#endif
             do iap = 1, n_ialist0
                ia = ia_list(iap)
                lmta1 = lmta(lmt1,ia)
#ifdef VPP
*vocl loop, unroll(4)
#endif
                do ib = 1, np_e
                   fsrt = 0.d0; fsit = 0.d0
                   do i = 1, iba(ik)
                      fsrt = fsrt + ar_x(i,iap)*zaj_l(i,ib,ik,1)
                      fsit = fsit + ai_x(i,iap)*zaj_l(i,ib,ik,1)
                   end do
                   bWr(lmta1,iter_rmm,ib) = fsrt; bWi(lmta1,iter_rmm,ib) = fsit
                end do
             end do
          else if(kimg == 2) then
             if(k_symmetry(ik) == GAMMA) then
                il2 = mod(il1,2) ! il2 = 1 -> l=0,2,4,6,...: il2 = 0 -> l=1,3,5,7,...
#ifdef HIUX
*poption indep
*poption parallel
#endif
                do iap = 1, n_ialist0
                   ia = ia_list(iap)
                   lmta1 = lmta(lmt1,ia)
#ifdef VPP
*vocl loop, unroll(4)
#endif
                   do ib = 1, np_e
                      fsrt = 0.d0; fsit = 0.d0
                      if(il2 == 1) then
                         do i = 2, iba(ik)
                            fsrt = fsrt + ar_x(i,iap)*zaj_l(i,ib,ik,1)-ai_x(i,iap)*zaj_l(i,ib,ik,2)
                         end do
                         fsrt = 2.d0*fsrt + zaj_l(1,ib,ik,1)*ar_x(1,iap)
                      else
                         do i = 2, iba(ik)
                            fsit = fsit + ai_x(i,iap)*zaj_l(i,ib,ik,1)+ar_x(i,iap)*zaj_l(i,ib,ik,2)
                         end do
                         fsit = 2.d0*fsit
                      end if
                      bWr(lmta1,iter_rmm,ib) = fsrt; bWi(lmta1,iter_rmm,ib) = fsit
                   end do
                end do
             else
#ifdef HIUX
*poption indep
*poption parallel
#endif
                do iap = 1, n_ialist0
                   ia = ia_list(iap)
                   lmta1 = lmta(lmt1,ia)
#ifdef VPP
*vocl loop, unroll(4)
#endif
                   do ib = 1, np_e
                      fsrt = 0.d0; fsit = 0.d0
                      do i    = 1,iba(ik)
                         fsrt = fsrt + ar_x(i,iap)*zaj_l(i,ib,ik,1)-ai_x(i,iap)*zaj_l(i,ib,ik,2)
                         fsit = fsit + ai_x(i,iap)*zaj_l(i,ib,ik,1)+ar_x(i,iap)*zaj_l(i,ib,ik,2)
                      end do
                      bWr(lmta1,iter_rmm,ib) = fsrt; bWi(lmta1,iter_rmm,ib) = fsit
                   end do
                end do
             end if
          end if

          i = mod(il1,4)
          do iap = 1, n_ialist0
             ia = ia_list(iap)
             lmta1 = lmta(lmt1,ia)
             if(i == 2) then
                do ib = 1, np_e
                   fsrt = bWr(lmta1,iter_rmm,ib); fsit = bWi(lmta1,iter_rmm,ib)
                   bWr(lmta1,iter_rmm,ib) = -fsit; bWi(lmta1,iter_rmm,ib) = fsrt
                end do
             else if(i == 3) then
                do ib = 1, np_e
                   fsrt = bWr(lmta1,iter_rmm,ib); fsit = bWi(lmta1,iter_rmm,ib)
                   bWr(lmta1,iter_rmm,ib) = -fsrt; bWi(lmta1,iter_rmm,ib) = -fsit
                end do
             else if(i == 0) then
                do ib = 1, np_e
                   fsrt = bWr(lmta1,iter_rmm,ib); fsit = bWi(lmta1,iter_rmm,ib)
                   bWr(lmta1,iter_rmm,ib) =  fsit; bWi(lmta1,iter_rmm,ib) = -fsrt
                end do
             end if
          end do
       end do
       call tstatc0_end(id_sname)
     end subroutine sumset_rmm_all3

     subroutine add_vnlph_l_without_eko_part3()
        real(kind=DP) :: bWFr, bWFi
        integer       :: i, ib, ia, iap, lmta2
        integer :: id_sname = -1
        call tstatc0_begin('add_vnlph_l_without_eko_part3 ',id_sname)

        if(kimg == 1) then
!!$#ifdef HIUX
!!$*poption parallel
!!$#endif
           do iap = 1, n_ialist0
              ia = ia_list(iap)
              lmta2 = lmta(lmt2,ia)
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr = bWr(lmta2,iter_rmm,ib); bWFi = bWi(lmta2,iter_rmm,ib)
                do i = 1, iba(ik)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*sc_x(i,iap)-bWFi*ss_x(i,iap)
                enddo
             end do
          end do
       else if(kimg == 2) then
!!$#ifdef HIUX
!!$*poption parallel
!!$#endif
          do iap = 1, n_ialist0
             ia = ia_list(iap)
             lmta2 = lmta(lmt2,ia)
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr = bWr(lmta2,iter_rmm,ib); bWFi = bWi(lmta2,iter_rmm,ib)
                do i = 1, iba(ik)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*sc_x(i,iap)-bWFi*ss_x(i,iap)
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi*sc_x(i,iap)+bWFr*ss_x(i,iap)
                enddo
             end do
          end do
       endif
       call tstatc0_end(id_sname)
     end subroutine add_vnlph_l_without_eko_part3

     subroutine add_vnlph_l_with_eko_part3()
        real(kind=DP) :: bWFr, bWFi, e, sceqc, sseqs
        integer       :: i, ib, ia, iap, lmta2
#ifdef HIUX
        integer       :: lmta2_1, lmta2_2, lmta2_3, lmta2_4, lmta2_5, lmta2_6, lmta2_7, lmta2_8
        real(kind=DP) :: bWFr_1, bWFr_2, bWFr_3, bWFr_4, bWFr_5, bWFr_6, bWFr_7, bWFr_8
        real(kind=DP) :: bWFi_1, bWFi_2, bWFi_3, bWFi_4, bWFi_5, bWFi_6, bWFi_7, bWFi_8
        real(kind=DP) :: sceqc_1, sseqs_1, sceqc_2, sseqs_2, sceqc_3, sseqs_3, sceqc_4, sseqs_4 &
            &         , sceqc_5, sseqs_5, sceqc_6, sseqs_6, sceqc_7, sseqs_7, sceqc_8, sseqs_8
#endif

       integer :: id_sname = -1
       call tstatc0_begin('add_vnlph_l_with_eko_part3 ',id_sname)

#ifdef HIUX
       if(kimg == 1) then
          if(n_ialist0 == 1) then
             lmta2 = lmta(lmt2,ia_list(1))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr = bWr(lmta2,iter_rmm,ib); bWFi = bWi(lmta2,iter_rmm,ib); e = eko_l(ib,ik)
                do i = 1, iba(ik)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*(sc_x(i,1)-e*qc_x(i,1))-bWFi*(ss_x(i,1)-e*qs_x(i,1))
                enddo
             end do
          else if(n_ialist0 == 2) then
             lmta2_1 = lmta(lmt2,ia_list(1)); lmta2_2 = lmta(lmt2,ia_list(2))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib =1, np_e
                bWFr_1 = bWr(lmta2_1,iter_rmm,ib); bWFi_1 = bWi(lmta2_1,iter_rmm,ib)
                bWFr_2 = bWr(lmta2_2,iter_rmm,ib); bWFi_2 = bWi(lmta2_2,iter_rmm,ib); e = eko_l(ib,ik)
                do i = 1, iba(ik)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*(sc_x(i,1)-e*qc_x(i,1))-bWFi_1*(ss_x(i,1)-e*qs_x(i,1)) &
                        &                            + bWFr_2*(sc_x(i,2)-e*qc_x(i,2))-bWFi_2*(ss_x(i,2)-e*qs_x(i,2))
                end do
             end do
          else if(n_ialist0 == 3) then
             lmta2_1 = lmta(lmt2,ia_list(1)); lmta2_2 = lmta(lmt2,ia_list(2)); lmta2_3 = lmta(lmt2,ia_list(3))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib =1, np_e
                bWFr_1 = bWr(lmta2_1,iter_rmm,ib); bWFi_1 = bWi(lmta2_1,iter_rmm,ib)
                bWFr_2 = bWr(lmta2_2,iter_rmm,ib); bWFi_2 = bWi(lmta2_2,iter_rmm,ib)
                bWFr_3 = bWr(lmta2_3,iter_rmm,ib); bWFi_3 = bWi(lmta2_3,iter_rmm,ib); e = eko_l(ib,ik)
                do i = 1, iba(ik)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*(sc_x(i,1)-e*qc_x(i,1))-bWFi_1*(ss_x(i,1)-e*qs_x(i,1)) &
                        &                            + bWFr_2*(sc_x(i,2)-e*qc_x(i,2))-bWFi_2*(ss_x(i,2)-e*qs_x(i,2)) &
                        &                            + bWFr_3*(sc_x(i,3)-e*qc_x(i,3))-bWFi_3*(ss_x(i,3)-e*qs_x(i,3))
                end do
             end do
          else if(n_ialist0 >= 4) then
             lmta2_1=lmta(lmt2,ia_list(1)); lmta2_2=lmta(lmt2,ia_list(2)); lmta2_3=lmta(lmt2,ia_list(3)); lmta2_4=lmta(lmt2,ia_list(4))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib =1, np_e
                bWFr_1 = bWr(lmta2_1,iter_rmm,ib); bWFi_1 = bWi(lmta2_1,iter_rmm,ib)
                bWFr_2 = bWr(lmta2_2,iter_rmm,ib); bWFi_2 = bWi(lmta2_2,iter_rmm,ib)
                bWFr_3 = bWr(lmta2_3,iter_rmm,ib); bWFi_3 = bWi(lmta2_3,iter_rmm,ib)
                bWFr_4 = bWr(lmta2_4,iter_rmm,ib); bWFi_4 = bWi(lmta2_4,iter_rmm,ib); e = eko_l(ib,ik)
                do i = 1, iba(ik)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*(sc_x(i,1)-e*qc_x(i,1))-bWFi_1*(ss_x(i,1)-e*qs_x(i,1)) &
                        &                            + bWFr_2*(sc_x(i,2)-e*qc_x(i,2))-bWFi_2*(ss_x(i,2)-e*qs_x(i,2)) &
                        &                            + bWFr_3*(sc_x(i,3)-e*qc_x(i,3))-bWFi_3*(ss_x(i,3)-e*qs_x(i,3)) &
                        &                            + bWFr_4*(sc_x(i,4)-e*qc_x(i,4))-bWFi_4*(ss_x(i,4)-e*qs_x(i,4))
                end do
             end do
          end if
          if(n_ialist0 > 4) then
             do iap = 5, n_ialist0
                ia = ia_list(iap)
                lmta2 = lmta(lmt2,ia)
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
                do ib = 1, np_e
                   bWFr = bWr(lmta2,iter_rmm,ib); bWFi = bWi(lmta2,iter_rmm,ib); e = eko_l(ib,ik)
                   do i = 1, iba(ik)
                      vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*(sc_x(i,iap)-e*qc_x(i,iap))-bWFi*(ss_x(i,iap)-e*qs_x(i,iap))
                   enddo
                end do
             end do
          end if
       else if(kimg == 2) then
         if(n_ialist0 == 1) then
             lmta2 = lmta(lmt2,ia_list(1))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr = bWr(lmta2,iter_rmm,ib); bWFi = bWi(lmta2,iter_rmm,ib); e = eko_l(ib,ik)
                do i = 1, iba(ik)
                   sceqc = sc_x(i,1)-e*qc_x(i,1);  sseqs = ss_x(i,1)-e*qs_x(i,1)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*sceqc-bWFi*sseqs
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi*sceqc+bWFr*sseqs
                end do
             end do
          else if(n_ialist0 == 2) then
             lmta2_1 = lmta(lmt2,ia_list(1)); lmta2_2 = lmta(lmt2,ia_list(2))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr_1 = bWr(lmta2_1,iter_rmm,ib); bWFi_1 = bWi(lmta2_1,iter_rmm,ib)
                bWFr_2 = bWr(lmta2_2,iter_rmm,ib); bWFi_2 = bWi(lmta2_2,iter_rmm,ib); e = eko_l(ib,ik)
                do i = 1, iba(ik)
                   sceqc_1 = sc_x(i,1)-e*qc_x(i,1);  sseqs_1 = ss_x(i,1)-e*qs_x(i,1)
                   sceqc_2 = sc_x(i,2)-e*qc_x(i,2);  sseqs_2 = ss_x(i,2)-e*qs_x(i,2)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*sceqc_1-bWFi_1*sseqs_1 + bWFr_2*sceqc_2-bWFi_2*sseqs_2
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi_1*sceqc_1+bWFr_1*sseqs_1 + bWFi_2*sceqc_2+bWFr_2*sseqs_2
                end do
             end do
          else if(n_ialist0 == 3) then
             lmta2_1 = lmta(lmt2,ia_list(1)); lmta2_2 = lmta(lmt2,ia_list(2)); lmta2_3 = lmta(lmt2,ia_list(3))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr_1 = bWr(lmta2_1,iter_rmm,ib); bWFi_1 = bWi(lmta2_1,iter_rmm,ib)
                bWFr_2 = bWr(lmta2_2,iter_rmm,ib); bWFi_2 = bWi(lmta2_2,iter_rmm,ib)
                bWFr_3 = bWr(lmta2_3,iter_rmm,ib); bWFi_3 = bWi(lmta2_3,iter_rmm,ib); e = eko_l(ib,ik)
                do i = 1, iba(ik)
                   sceqc_1 = sc_x(i,1)-e*qc_x(i,1);  sseqs_1 = ss_x(i,1)-e*qs_x(i,1)
                   sceqc_2 = sc_x(i,2)-e*qc_x(i,2);  sseqs_2 = ss_x(i,2)-e*qs_x(i,2)
                   sceqc_3 = sc_x(i,3)-e*qc_x(i,3);  sseqs_3 = ss_x(i,3)-e*qs_x(i,3)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*sceqc_1-bWFi_1*sseqs_1 + bWFr_2*sceqc_2-bWFi_2*sseqs_2 &
                        &                            + bWFr_3*sceqc_3-bWFi_3*sseqs_3
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi_1*sceqc_1+bWFr_1*sseqs_1 + bWFi_2*sceqc_2+bWFr_2*sseqs_2 &
                        &                            + bWFi_3*sceqc_3+bWFr_3*sseqs_3
                end do
             end do
          else if(n_ialist0 == 4) then
             lmta2_1=lmta(lmt2,ia_list(1)); lmta2_2=lmta(lmt2,ia_list(2)); lmta2_3=lmta(lmt2,ia_list(3));lmta2_4=lmta(lmt2,ia_list(4))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr_1 = bWr(lmta2_1,iter_rmm,ib); bWFi_1 = bWi(lmta2_1,iter_rmm,ib)
                bWFr_2 = bWr(lmta2_2,iter_rmm,ib); bWFi_2 = bWi(lmta2_2,iter_rmm,ib)
                bWFr_3 = bWr(lmta2_3,iter_rmm,ib); bWFi_3 = bWi(lmta2_3,iter_rmm,ib)
                bWFr_4 = bWr(lmta2_4,iter_rmm,ib); bWFi_4 = bWi(lmta2_4,iter_rmm,ib); e = eko_l(ib,ik)
                do i = 1, iba(ik)
                   sceqc_1 = sc_x(i,1)-e*qc_x(i,1);  sseqs_1 = ss_x(i,1)-e*qs_x(i,1)
                   sceqc_2 = sc_x(i,2)-e*qc_x(i,2);  sseqs_2 = ss_x(i,2)-e*qs_x(i,2)
                   sceqc_3 = sc_x(i,3)-e*qc_x(i,3);  sseqs_3 = ss_x(i,3)-e*qs_x(i,3)
                   sceqc_4 = sc_x(i,4)-e*qc_x(i,4);  sseqs_4 = ss_x(i,4)-e*qs_x(i,4)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*sceqc_1-bWFi_1*sseqs_1 + bWFr_2*sceqc_2-bWFi_2*sseqs_2 &
                        &                            + bWFr_3*sceqc_3-bWFi_3*sseqs_3 + bWFr_4*sceqc_4-bWFi_4*sseqs_4
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi_1*sceqc_1+bWFr_1*sseqs_1 + bWFi_2*sceqc_2+bWFr_2*sseqs_2 &
                        &                            + bWFi_3*sceqc_3+bWFr_3*sseqs_3 + bWFi_4*sceqc_4+bWFr_4*sseqs_4
                end do
             end do
          else if(n_ialist0 == 5) then
             lmta2_1=lmta(lmt2,ia_list(1)); lmta2_2=lmta(lmt2,ia_list(2)); lmta2_3=lmta(lmt2,ia_list(3));lmta2_4=lmta(lmt2,ia_list(4))
             lmta2_5=lmta(lmt2,ia_list(5))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr_1 = bWr(lmta2_1,iter_rmm,ib); bWFi_1 = bWi(lmta2_1,iter_rmm,ib)
                bWFr_2 = bWr(lmta2_2,iter_rmm,ib); bWFi_2 = bWi(lmta2_2,iter_rmm,ib)
                bWFr_3 = bWr(lmta2_3,iter_rmm,ib); bWFi_3 = bWi(lmta2_3,iter_rmm,ib)
                bWFr_4 = bWr(lmta2_4,iter_rmm,ib); bWFi_4 = bWi(lmta2_4,iter_rmm,ib)
                bWFr_5 = bWr(lmta2_5,iter_rmm,ib); bWFi_5 = bWi(lmta2_5,iter_rmm,ib); e = eko_l(ib,ik)
                do i = 1, iba(ik)
                   sceqc_1 = sc_x(i,1)-e*qc_x(i,1);  sseqs_1 = ss_x(i,1)-e*qs_x(i,1)
                   sceqc_2 = sc_x(i,2)-e*qc_x(i,2);  sseqs_2 = ss_x(i,2)-e*qs_x(i,2)
                   sceqc_3 = sc_x(i,3)-e*qc_x(i,3);  sseqs_3 = ss_x(i,3)-e*qs_x(i,3)
                   sceqc_4 = sc_x(i,4)-e*qc_x(i,4);  sseqs_4 = ss_x(i,4)-e*qs_x(i,4)
                   sceqc_5 = sc_x(i,5)-e*qc_x(i,5);  sseqs_5 = ss_x(i,5)-e*qs_x(i,5)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*sceqc_1-bWFi_1*sseqs_1 + bWFr_2*sceqc_2-bWFi_2*sseqs_2 &
                       &                            + bWFr_3*sceqc_3-bWFi_3*sseqs_3 + bWFr_4*sceqc_4-bWFi_4*sseqs_4 &
                       &                            + bWFr_5*sceqc_5-bWFi_5*sseqs_5
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi_1*sceqc_1+bWFr_1*sseqs_1 + bWFi_2*sceqc_2+bWFr_2*sseqs_2 &
                       &                            + bWFi_3*sceqc_3+bWFr_3*sseqs_3 + bWFi_4*sceqc_4+bWFr_4*sseqs_4 &
                       &                            + bWFi_5*sceqc_5+bWFr_3*sseqs_5
                end do
             end do
          else if(n_ialist0 == 6) then
             lmta2_1=lmta(lmt2,ia_list(1)); lmta2_2=lmta(lmt2,ia_list(2)); lmta2_3=lmta(lmt2,ia_list(3));lmta2_4=lmta(lmt2,ia_list(4))
             lmta2_5=lmta(lmt2,ia_list(5)); lmta2_6=lmta(lmt2,ia_list(6))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr_1 = bWr(lmta2_1,iter_rmm,ib); bWFi_1 = bWi(lmta2_1,iter_rmm,ib)
                bWFr_2 = bWr(lmta2_2,iter_rmm,ib); bWFi_2 = bWi(lmta2_2,iter_rmm,ib)
                bWFr_3 = bWr(lmta2_3,iter_rmm,ib); bWFi_3 = bWi(lmta2_3,iter_rmm,ib)
                bWFr_4 = bWr(lmta2_4,iter_rmm,ib); bWFi_4 = bWi(lmta2_4,iter_rmm,ib)
                bWFr_5 = bWr(lmta2_5,iter_rmm,ib); bWFi_5 = bWi(lmta2_5,iter_rmm,ib)
                bWFr_6 = bWr(lmta2_6,iter_rmm,ib); bWFi_6 = bWi(lmta2_6,iter_rmm,ib); e = eko_l(ib,ik)
                do i = 1, iba(ik)
                   sceqc_1 = sc_x(i,1)-e*qc_x(i,1);  sseqs_1 = ss_x(i,1)-e*qs_x(i,1)
                   sceqc_2 = sc_x(i,2)-e*qc_x(i,2);  sseqs_2 = ss_x(i,2)-e*qs_x(i,2)
                   sceqc_3 = sc_x(i,3)-e*qc_x(i,3);  sseqs_3 = ss_x(i,3)-e*qs_x(i,3)
                   sceqc_4 = sc_x(i,4)-e*qc_x(i,4);  sseqs_4 = ss_x(i,4)-e*qs_x(i,4)
                   sceqc_5 = sc_x(i,5)-e*qc_x(i,5);  sseqs_5 = ss_x(i,5)-e*qs_x(i,5)
                   sceqc_6 = sc_x(i,6)-e*qc_x(i,6);  sseqs_6 = ss_x(i,6)-e*qs_x(i,6)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*sceqc_1-bWFi_1*sseqs_1 + bWFr_2*sceqc_2-bWFi_2*sseqs_2 &
                        &                            + bWFr_3*sceqc_3-bWFi_3*sseqs_3 + bWFr_4*sceqc_4-bWFi_4*sseqs_4 &
                        &                            + bWFr_5*sceqc_5-bWFi_5*sseqs_5 + bWFr_6*sceqc_6-bWFi_6*sseqs_6
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi_1*sceqc_1+bWFr_1*sseqs_1 + bWFi_2*sceqc_2+bWFr_2*sseqs_2 &
                        &                            + bWFi_3*sceqc_3+bWFr_3*sseqs_3 + bWFi_4*sceqc_4+bWFr_4*sseqs_4 &
                        &                            + bWFi_5*sceqc_5+bWFr_5*sseqs_5 + bWFi_6*sceqc_6+bWFr_6*sseqs_6
                end do
             end do
          else if(n_ialist0 == 7) then
             lmta2_1=lmta(lmt2,ia_list(1)); lmta2_2=lmta(lmt2,ia_list(2)); lmta2_3=lmta(lmt2,ia_list(3));lmta2_4=lmta(lmt2,ia_list(4))
             lmta2_5=lmta(lmt2,ia_list(5)); lmta2_6=lmta(lmt2,ia_list(6)); lmta2_7=lmta(lmt2,ia_list(7))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr_1 = bWr(lmta2_1,iter_rmm,ib); bWFi_1 = bWi(lmta2_1,iter_rmm,ib)
                bWFr_2 = bWr(lmta2_2,iter_rmm,ib); bWFi_2 = bWi(lmta2_2,iter_rmm,ib)
                bWFr_3 = bWr(lmta2_3,iter_rmm,ib); bWFi_3 = bWi(lmta2_3,iter_rmm,ib)
                bWFr_4 = bWr(lmta2_4,iter_rmm,ib); bWFi_4 = bWi(lmta2_4,iter_rmm,ib)
                bWFr_5 = bWr(lmta2_5,iter_rmm,ib); bWFi_5 = bWi(lmta2_5,iter_rmm,ib)
                bWFr_6 = bWr(lmta2_6,iter_rmm,ib); bWFi_6 = bWi(lmta2_6,iter_rmm,ib)
                bWFr_7 = bWr(lmta2_7,iter_rmm,ib); bWFi_7 = bWi(lmta2_7,iter_rmm,ib); e = eko_l(ib,ik)
                do i = 1, iba(ik)
                   sceqc_1 = sc_x(i,1)-e*qc_x(i,1);  sseqs_1 = ss_x(i,1)-e*qs_x(i,1)
                   sceqc_2 = sc_x(i,2)-e*qc_x(i,2);  sseqs_2 = ss_x(i,2)-e*qs_x(i,2)
                   sceqc_3 = sc_x(i,3)-e*qc_x(i,3);  sseqs_3 = ss_x(i,3)-e*qs_x(i,3)
                   sceqc_4 = sc_x(i,4)-e*qc_x(i,4);  sseqs_4 = ss_x(i,4)-e*qs_x(i,4)
                   sceqc_5 = sc_x(i,5)-e*qc_x(i,5);  sseqs_5 = ss_x(i,5)-e*qs_x(i,5)
                   sceqc_6 = sc_x(i,6)-e*qc_x(i,6);  sseqs_6 = ss_x(i,6)-e*qs_x(i,6)
                   sceqc_7 = sc_x(i,7)-e*qc_x(i,7);  sseqs_7 = ss_x(i,7)-e*qs_x(i,7)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*sceqc_1-bWFi_1*sseqs_1 + bWFr_2*sceqc_2-bWFi_2*sseqs_2 &
                        &                            + bWFr_3*sceqc_3-bWFi_3*sseqs_3 + bWFr_4*sceqc_4-bWFi_4*sseqs_4 &
                        &                            + bWFr_5*sceqc_5-bWFi_5*sseqs_5 + bWFr_6*sceqc_6-bWFi_6*sseqs_6 &
                        &                            + bWFr_7*sceqc_7-bWFi_7*sseqs_7
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi_1*sceqc_1+bWFr_1*sseqs_1 + bWFi_2*sceqc_2+bWFr_2*sseqs_2 &
                        &                            + bWFi_3*sceqc_3+bWFr_3*sseqs_3 + bWFi_4*sceqc_4+bWFr_4*sseqs_4 &
                        &                            + bWFi_5*sceqc_5+bWFr_5*sseqs_5 + bWFi_6*sceqc_6+bWFr_6*sseqs_6 &
                        &                            + bWFi_7*sceqc_7+bWFr_7*sseqs_7
                end do
             end do
          else if(n_ialist0 >= 8) then
             lmta2_1=lmta(lmt2,ia_list(1)); lmta2_2=lmta(lmt2,ia_list(2)); lmta2_3=lmta(lmt2,ia_list(3));lmta2_4=lmta(lmt2,ia_list(4))
             lmta2_5=lmta(lmt2,ia_list(5)); lmta2_6=lmta(lmt2,ia_list(6)); lmta2_7=lmta(lmt2,ia_list(7));lmta2_8=lmta(lmt2,ia_list(8))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr_1 = bWr(lmta2_1,iter_rmm,ib); bWFi_1 = bWi(lmta2_1,iter_rmm,ib)
                bWFr_2 = bWr(lmta2_2,iter_rmm,ib); bWFi_2 = bWi(lmta2_2,iter_rmm,ib)
                bWFr_3 = bWr(lmta2_3,iter_rmm,ib); bWFi_3 = bWi(lmta2_3,iter_rmm,ib)
                bWFr_4 = bWr(lmta2_4,iter_rmm,ib); bWFi_4 = bWi(lmta2_4,iter_rmm,ib)
                bWFr_5 = bWr(lmta2_5,iter_rmm,ib); bWFi_5 = bWi(lmta2_5,iter_rmm,ib)
                bWFr_6 = bWr(lmta2_6,iter_rmm,ib); bWFi_6 = bWi(lmta2_6,iter_rmm,ib)
                bWFr_7 = bWr(lmta2_7,iter_rmm,ib); bWFi_7 = bWi(lmta2_7,iter_rmm,ib)
                bWFr_8 = bWr(lmta2_8,iter_rmm,ib); bWFi_8 = bWi(lmta2_8,iter_rmm,ib); e = eko_l(ib,ik)
                do i = 1, iba(ik)
                   sceqc_1 = sc_x(i,1)-e*qc_x(i,1);  sseqs_1 = ss_x(i,1)-e*qs_x(i,1)
                   sceqc_2 = sc_x(i,2)-e*qc_x(i,2);  sseqs_2 = ss_x(i,2)-e*qs_x(i,2)
                   sceqc_3 = sc_x(i,3)-e*qc_x(i,3);  sseqs_3 = ss_x(i,3)-e*qs_x(i,3)
                   sceqc_4 = sc_x(i,4)-e*qc_x(i,4);  sseqs_4 = ss_x(i,4)-e*qs_x(i,4)
                   sceqc_5 = sc_x(i,5)-e*qc_x(i,5);  sseqs_5 = ss_x(i,5)-e*qs_x(i,5)
                   sceqc_6 = sc_x(i,6)-e*qc_x(i,6);  sseqs_6 = ss_x(i,6)-e*qs_x(i,6)
                   sceqc_7 = sc_x(i,7)-e*qc_x(i,7);  sseqs_7 = ss_x(i,7)-e*qs_x(i,7)
                   sceqc_8 = sc_x(i,8)-e*qc_x(i,8);  sseqs_8 = ss_x(i,8)-e*qs_x(i,8)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*sceqc_1-bWFi_1*sseqs_1 + bWFr_2*sceqc_2-bWFi_2*sseqs_2 &
                        &                            + bWFr_3*sceqc_3-bWFi_3*sseqs_3 + bWFr_4*sceqc_4-bWFi_4*sseqs_4 &
                        &                            + bWFr_5*sceqc_5-bWFi_5*sseqs_5 + bWFr_6*sceqc_6-bWFi_6*sseqs_6 &
                        &                            + bWFr_7*sceqc_7-bWFi_7*sseqs_7 + bWFr_8*sceqc_8-bWFi_8*sseqs_8
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi_1*sceqc_1+bWFr_1*sseqs_1 + bWFi_2*sceqc_2+bWFr_2*sseqs_2 &
                        &                            + bWFi_3*sceqc_3+bWFr_3*sseqs_3 + bWFi_4*sceqc_4+bWFr_4*sseqs_4 &
                        &                            + bWFi_5*sceqc_5+bWFr_5*sseqs_5 + bWFi_6*sceqc_6+bWFr_6*sseqs_6 &
                        &                            + bWFi_7*sceqc_7+bWFr_7*sseqs_7 + bWFi_8*sceqc_8+bWFr_8*sseqs_8
                end do
             end do
          end if
          if(n_ialist0 > 8) then
             do iap = 9, n_ialist0
                ia = ia_list(iap)
                lmta2 = lmta(lmt2,ia)
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
                do ib = 1, np_e
                   bWFr = bWr(lmta2,iter_rmm,ib); bWFi = bWi(lmta2,iter_rmm,ib); e = eko_l(ib,ik)
                   do i = 1, iba(ik)
                      sceqc = sc_x(i,iap)-e*qc_x(i,iap)
                      sseqs = ss_x(i,iap)-e*qs_x(i,iap)
                      vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*sceqc-bWFi*sseqs
                      vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi*sceqc+bWFr*sseqs
                   end do
                enddo
            end do
         end if
       end if

#else
!  --> for ordinary scaler and vector machines but HIUX
!  = #ifndef HIUX

       if(kimg == 1) then
          do iap = 1, n_ialist0
             ia = ia_list(iap)
             lmta2 = lmta(lmt2,ia)
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr = bWr(lmta2,iter_rmm,ib); bWFi = bWi(lmta2,iter_rmm,ib); e = eko_l(ib,ik)
                do i = 1, iba(ik)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*(sc_x(i,iap)-e*qc_x(i,iap))-bWFi*(ss_x(i,iap)-e*qs_x(i,iap))
                enddo
             end do
          end do
       else if(kimg == 2) then
          do iap = 1, n_ialist0
             ia = ia_list(iap)
             lmta2 = lmta(lmt2,ia)
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr = bWr(lmta2,iter_rmm,ib); bWFi = bWi(lmta2,iter_rmm,ib); e = eko_l(ib,ik)
                do i = 1, iba(ik)
                   sceqc = sc_x(i,iap)-e*qc_x(i,iap)
                   sseqs = ss_x(i,iap)-e*qs_x(i,iap)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*sceqc-bWFi*sseqs
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi*sceqc+bWFr*sseqs
                end do
             enddo
          end do
       endif
#endif
       call tstatc0_end(id_sname)
     end subroutine add_vnlph_l_with_eko_part3

     subroutine Vnonlocal_W_part_sum_over_lmt1b()
       integer             :: lmt1, lmtt1, il1, im1, il11, mdl, iksnl,il2,im2, ia_p, ia, i, ii
       real(kind=DP)       :: tmp

       integer :: id_sname = -1
       call tstatc0_begin('Vnonlocal_W_part_sum_over_lmt1 ',id_sname)

       iksnl = (ik-1)/nspin + 1
       il2   = ltp(lmt2,it)
       im2   = mtp(lmt2,it)

       sc_x = 0.d0; ss_x = 0.d0
       if(mdvdb == EXECUT) then
          qc_x = 0.d0; qs_x = 0.d0
       endif

       do lmt1 = 1,ilmt(it)
          lmtt1 = lmtt(lmt1,it)
          il1   = ltp(lmt1,it)
          im1   = mtp(lmt1,it)
          il11  = il1 - 1
          mdl   = mod(il11,4)

          do ia_p = 1, n_ialist0
             ia = ia_list(ia_p)
             if(il1 == il2 .and. im1 == im2) then
!!$                tmp = (dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)) * iwei(ia)
                if(ipaw(it)==0) then
                    tmp = (dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)) * iwei(ia)
                else
                    tmp = (dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)) * iwei(ia)
                end if
             else
!!$                tmp = vlhxcQ(lmt1,lmt2,ia,ispin) * iwei(ia)
                if(ipaw(it)==0) then
                    tmp = vlhxcQ(lmt1,lmt2,ia,ispin) * iwei(ia)
                else
                    tmp = (dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)) * iwei(ia)
                end if
             endif
             if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
             fq(ia_p,1) = tmp
          end do

          if(mdvdb == EXECUT .and. il1 == il2 .and. im1 == im2) then
             do ia_p = 1, n_ialist0
                ia = ia_list(ia_p)
                tmp = q(lmt1,lmt2,it)*iwei(ia)
                if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
                fq(ia_p,2) = tmp
             end do

             if(mdl == 0 .or. mdl == 2) then
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
                do ia_p = 1, n_ialist0
                   do i = 1, iba(ik)
                      sc_x(i,ia_p) = sc_x(i,ia_p) + fq(ia_p,1)*zfcos_x(i,ia_p)*snl(i,lmtt1,iksnl)
                      ss_x(i,ia_p) = ss_x(i,ia_p) - fq(ia_p,1)*zfsin_x(i,ia_p)*snl(i,lmtt1,iksnl)
                      qc_x(i,ia_p) = qc_x(i,ia_p) + fq(ia_p,2)*zfcos_x(i,ia_p)*snl(i,lmtt1,iksnl)
                      qs_x(i,ia_p) = qs_x(i,ia_p) - fq(ia_p,2)*zfsin_x(i,ia_p)*snl(i,lmtt1,iksnl)
                   end do
                end do
             else if(mdl == 1 .or. mdl == 3) then
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
                do ia_p = 1, n_ialist0
                   do i = 1, iba(ik)
                      sc_x(i,ia_p) = sc_x(i,ia_p) - fq(ia_p,1)*zfsin_x(i,ia_p)*snl(i,lmtt1,iksnl)
                      ss_x(i,ia_p) = ss_x(i,ia_p) - fq(ia_p,1)*zfcos_x(i,ia_p)*snl(i,lmtt1,iksnl)
                      qc_x(i,ia_p) = qc_x(i,ia_p) - fq(ia_p,2)*zfsin_x(i,ia_p)*snl(i,lmtt1,iksnl)
                      qs_x(i,ia_p) = qs_x(i,ia_p) - fq(ia_p,2)*zfcos_x(i,ia_p)*snl(i,lmtt1,iksnl)
                   end do
                end do
             end if
          else
             if(mdl == 0 .or. mdl == 2) then
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
                do ia_p = 1, n_ialist0
                   do i = 1, iba(ik)
                      sc_x(i,ia_p) = sc_x(i,ia_p) + fq(ia_p,1)*zfcos_x(i,ia_p)*snl(i,lmtt1,iksnl)
                      ss_x(i,ia_p) = ss_x(i,ia_p) - fq(ia_p,1)*zfsin_x(i,ia_p)*snl(i,lmtt1,iksnl)
                   end do
                end do
             else if(mdl == 1 .or. mdl == 3) then
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
                do ia_p = 1, n_ialist0
                   do i = 1, iba(ik)
                      sc_x(i,ia_p) = sc_x(i,ia_p) - fq(ia_p,1)*zfsin_x(i,ia_p)*snl(i,lmtt1,iksnl)
                      ss_x(i,ia_p) = ss_x(i,ia_p) - fq(ia_p,1)*zfcos_x(i,ia_p)*snl(i,lmtt1,iksnl)
                   end do
                end do
             end if
          end if
       end do

       call tstatc0_end(id_sname)
     end subroutine Vnonlocal_W_part_sum_over_lmt1b
#endif
  end subroutine Vnonlocal_W_RMMn

! ================================= added by K. Tagami ======================== 11.0
  subroutine Vnonlocal_W_RMMn_noncl( nfout, ik, ispin, iter_rmm )
    integer, intent(in) :: nfout,ik,ispin, iter_rmm
    integer :: mdvdb, it, ia, lmt2, lmta2, ia_p,iap,lmta1
    integer :: n_ialist, n_ialist0, ia_start, n_iagroup, n_ia, ia_g, max_n_ialist0 ,nbmx_adj, kg1_adj
    real(kind=DP), allocatable,dimension(:,:) ::  ar_x, ai_x                 ! d(kg1_adj,n_ialist0)
    real(kind=DP), allocatable,dimension(:,:) ::  zfcos_x, zfsin_x           ! d(nbmx,n_ialist0)
    real(kind=DP), allocatable,target,dimension(:,:) :: sc_x,ss_x,qc_x,qs_x  ! d(kg1+,n_ialist))
    real(kind=DP), allocatable,target,dimension(:,:) :: sc_lmt,ss_lmt,qc_lmt,qs_lmt  ! d(kg1+,n_ialist))
    integer, allocatable, dimension(:) :: ia_list

    complex(kind=CMPLDP), allocatable, dimension(:,:) :: fq

#ifdef NONLOCAL_RMM_DGEMM
    real(kind=DP), allocatable, dimension(:,:) :: bWr_lmt, bWi_lmt
#endif
    integer :: id_sname = -1, id_sname2 = -1
    call tstatc0_begin('Vnonlocal_W_RMMn_noncl ',id_sname,1)

#ifdef NONLOCAL_RMM_DGEMM
    n_ialist = 16
#else
    n_ialist = 1
#ifdef HIUX
    n_ialist = 16
#endif
#ifdef VPP
    n_ialist = 1
#endif
#ifdef SX
    n_ialist = 8
#endif
#endif
    if(n_ialist <=0) call phase_error_with_msg(nfout, 'n_ialist is illegal <<Vnonlocal_W_RMMn>>',__LINE__,__FILE__)
    call find_max_n_ialist0(max_n_ialist0)
    if(iprirmm >= 2) write(nfout,'(" !Vnonlocal_W_RMMn: max_n_ialist0 = ",i8)') max_n_ialist0
    nbmx_adj = nbmx
    if(mod(nbmx_adj,2) == 0) nbmx_adj = nbmx_adj+1
    kg1_adj = kg1
    if(mod(kg1_adj,2) == 0) kg1_adj = kg1_adj+1

    allocate(zfsin_x(nbmx_adj,max_n_ialist0)); allocate(zfcos_x(nbmx_adj,max_n_ialist0))
    allocate(ar_x(kg1_adj,max_n_ialist0)); allocate(ai_x(kg1_adj,max_n_ialist0))
    allocate(ia_list(n_ialist)); ia_list = 0
    allocate(fq(n_ialist,2)); fq = 0.d0

    vnlph_l(:,:,:) = 0.d0
    Loop_ntyp: do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
#ifdef NONLOCAL_RMM_DGEMM
       call alloc_scssqcqs_lmt(kg1_adj,max_n_ialist0*ilmt(it),mdvdb)
       call alloc_scssqcqs_x(kg1_adj,max_n_ialist0,mdvdb)
#else
       call alloc_scssqcqs_x(kg1_adj,max_n_ialist0,mdvdb)
#endif

       n_ia = 0
       do ia = 1, natm
          if(ityp(ia) == it) n_ia = n_ia + 1
       end do

       n_iagroup = n_ia/n_ialist + 1
       ia_start = 1
       if (iprirmm >= 2) write(nfout,'(" !Vnonlocal_W_RMMn: n_iagroup = ",i8, " ityp = ",i8)') n_iagroup,it
       Loop_ia_group: do ia_g = 1, n_iagroup
          n_ialist0 = 0
          ia_list = 0
          AtomcountLoop: do ia = ia_start, natm
             if(ityp(ia) == it) then
                n_ialist0 = n_ialist0 + 1
                ia_list(n_ialist0) = ia
             end if
              if(n_ialist0 >= n_ialist) exit AtomcountLoop
          end do AtomcountLoop
          ia_start = ia+1
          if(n_ialist0 >= 1 )then
             if(iprirmm >= 2) write(nfout,'(" !m_ES_WF_by_RMM.Vnonlocal_W_RMMn: ia_list = ",8i8)') (ia_list(ia),ia=1,n_ialist0)

                  call tstatc0_begin('calc_phasek_b ',id_sname2)
             call calc_phasek_b(natm,pos,n_ialist0,ia_list,kgp,iba(ik) &
                  & ,ngabc,kg1,nbase(1,ik),nbmx_adj,zfcos_x,zfsin_x) ! b_E.S.
                  call tstatc0_end(id_sname2)

#ifdef NONLOCAL_RMM_DGEMM
             call alloc_bW_lmt(it)  ! alloc bWr_lmt, bWi_lmp (np_e, n_ialist0*ilmt(it))
             call sumset_rmm_all4_noncl()
                            ! zfcos,zfsin,sc,ss,zaj,snl -> bWr_tmp,bWi_tmp, bWr, bWi

             call Vnonlocal_W_part_sum_lmt4_noncl()           ! -> sc_x,ss_x,qc_x,qs_x
             call add_vnlph_l_part4_noncl(mdvdb)

             call dealloc_bW_lmt  ! alloc bWr_lmt, bWi_lmp (np_e, n_ialist0*ilmt(it))
#else
             call sumset_rmm_all3_noncl()
                                           ! zfcos,zfsin,sc,ss,zaj,snl -> bWr,bWi

             do lmt2 = 1, ilmt(it)
                call Vnonlocal_W_partsum_lmt1b_noncl()      ! -> sc_x,ss_x,qc_x,qs_x
                if(mdvdb == SKIP) then
                   call add_vnlph_l_no_eko_pt3_noncl()     ! -> vnlph_l
                else if(mdvdb == EXECUT) then
                   call add_vnlph_l_with_eko_pt3_noncl()        ! -> vnlph_l
                endif
             end do
#endif
          end if
       end do Loop_ia_group
       call dealloc_scssqcqs_x()
#ifdef NONLOCAL_RMM_DGEMM
       call dealloc_scssqcqs_lmt()
#endif
    end do Loop_ntyp

    deallocate(fq)
    deallocate(ai_x,ar_x)
    deallocate(ia_list)
    deallocate(zfsin_x,zfcos_x)
    call tstatc0_end(id_sname)

  contains

    subroutine alloc_scssqcqs_lmt(kg1_adj,m,mdvdb)
      integer, intent(in) :: kg1_adj,m,mdvdb
      allocate(sc_lmt(kg1_adj,m))
      allocate(ss_lmt(kg1_adj,m))
      if(mdvdb == EXECUT) then
         allocate(qc_lmt(kg1_adj,m))
         allocate(qs_lmt(kg1_adj,m))
      end if
    end subroutine alloc_scssqcqs_lmt

    subroutine dealloc_scssqcqs_lmt()
      deallocate(sc_lmt,ss_lmt)
      if(mdvdb == EXECUT) deallocate(qc_lmt,qs_lmt)
    end subroutine dealloc_scssqcqs_lmt

    subroutine alloc_scssqcqs_x(kg1_adj,m,mdvdb)
      integer, intent(in) :: kg1_adj,m,mdvdb
      allocate(sc_x(kg1_adj,m))
      allocate(ss_x(kg1_adj,m))
      if(mdvdb == EXECUT) then
         allocate(qc_x(kg1_adj,m))
         allocate(qs_x(kg1_adj,m))
      end if
    end subroutine alloc_scssqcqs_x

    subroutine dealloc_scssqcqs_x()
      if(allocated(qs_x)) deallocate(qs_x)
      if(allocated(qc_x)) deallocate(qc_x)
      deallocate(ss_x,sc_x)
    end subroutine dealloc_scssqcqs_x

    subroutine find_max_n_ialist0(n)
      integer, intent(out) :: n

      n = 0
      do it = 1, ntyp

         n_ia = 0
         do ia = 1, natm
            if(ityp(ia) == it) n_ia = n_ia + 1
         end do

         n_iagroup = n_ia/n_ialist + 1
         ia_start = 1
         if(iprirmm >= 2) write(nfout,'(" !m_ES_WF_by_RMM.Vnonlocal_W_RMMn: n_iagroup = ",i8, " ityp = ",i8)') n_iagroup,it
         do ia_g = 1, n_iagroup
            n_ialist0 = 0
            AtomcountLoop: do ia = ia_start, natm
               if(ityp(ia) == it)  n_ialist0 = n_ialist0 + 1
               if(n_ialist0 >= n_ialist) exit AtomcountLoop
            end do AtomcountLoop
            if(n < n_ialist0) n = n_ialist0
            ia_start = ia+1
         end do
      end do

     end subroutine find_max_n_ialist0

#ifdef NONLOCAL_RMM_DGEMM
    subroutine alloc_bW_lmt(it)
      integer, intent(in) :: it
      integer :: nsize1
      nsize1 = n_ialist0 * ilmt(it)
      allocate(bWr_lmt(nsize1,np_e))
      allocate(bWi_lmt(nsize1,np_e))
    end subroutine alloc_bW_lmt

    subroutine dealloc_bW_lmt()
      deallocate(bWr_lmt)
      deallocate(bWi_lmt)
    end subroutine dealloc_bW_lmt

    subroutine sumset_rmm_all4_noncl()
      integer  :: lmt1, lmtt1, lmta1, i, iksnl, ib, il1, ia, iap, il2
      integer :: M, N, K
      real(kind=DP) :: fsrt, fsit

      integer :: my_is

      real(kind=DP), allocatable, dimension(:,:) :: bWr_tmp, bWi_tmp
      integer :: id_sname = -1

      call tstatc0_begin('sumset_rmm_all4_noncl ',id_sname)

      allocate(bWr_tmp(n_ialist0,np_e)); bWr_tmp = 0.d0
      allocate(bWi_tmp(n_ialist0,np_e)); bWi_tmp = 0.d0

      iksnl = (ik-1)/ndim_spinor + 1
      my_is = mod( ik-1, ndim_spinor ) +1

      do lmt1 = 1, ilmt(it)
         lmtt1 = lmtt(lmt1,it)
         il1   = ltp(lmt1,it)
         do iap = 1, n_ialist0
            do i = 1, iba(ik)
               ar_x(i,iap) = zfcos_x(i,iap)*snl(i,lmtt1,iksnl)
               ai_x(i,iap) = zfsin_x(i,iap)*snl(i,lmtt1,iksnl)
            end do
         end do


         M = n_ialist0; N = np_e; K = iba(ik)
         if(kimg == 1) then
            call DGEMM__('T','N',M,N,K,1.d0,ar_x,kg1_adj,zaj_l(1,1,ik,1),kg1,0.d0,bWr_tmp,M)
            call DGEMM__('T','N',M,N,K,1.d0,ai_x,kg1_adj,zaj_l(1,1,ik,1),kg1,0.d0,bWi_tmp,M)

         else if(kimg == 2) then
            if(k_symmetry(ik) == GAMMA) then
               il2 = mod(il1,2) ! il2 = 1 -> l=0,2,4,6,...: il2 = 0 -> l=1,3,5,7,...
               if(il2 == 1) then
                  ai_x(1,1:n_ialist0) = 0.d0;  ar_x(1,1:n_ialist0) = ar_x(1,1:n_ialist0)*0.5d0
                  call DGEMM__('T','N', M,N,K, 1.d0,ar_x,kg1_adj,zaj_l(1,1,ik,1),kg1,0.d0,bWr_tmp,M)
                  call DGEMM__('T','N', M,N,K,-1.d0,ai_x,kg1_adj,zaj_l(1,1,ik,2),kg1,1.d0,bWr_tmp,M)
                  bWr_tmp = bWr_tmp*2.d0;   bWi_tmp = 0.d0
               else
                  ar_x(1,1:n_ialist0) = 0.d0; ai_x(1,1:n_ialist0) = 0.d0
                  call DGEMM__('T','N', M,N,K, 1.d0,ai_x,kg1_adj,zaj_l(1,1,ik,1),kg1,0.d0,bWi_tmp,M)
                  call DGEMM__('T','N', M,N,K, 1.d0,ar_x,kg1_adj,zaj_l(1,1,ik,2),kg1,1.d0,bWi_tmp,M)
                  bWr_tmp = 0.d0;           bWi_tmp = bWi_tmp*2.d0
               end if
            else
               call DGEMM__('T','N', M,N,K, 1.d0,ar_x,kg1_adj,zaj_l(1,1,ik,1),kg1,0.d0,bWr_tmp,M)
               call DGEMM__('T','N', M,N,K,-1.d0,ai_x,kg1_adj,zaj_l(1,1,ik,2),kg1,1.d0,bWr_tmp,M)
               call DGEMM__('T','N', M,N,K, 1.d0,ai_x,kg1_adj,zaj_l(1,1,ik,1),kg1,0.d0,bWi_tmp,M)
               call DGEMM__('T','N', M,N,K, 1.d0,ar_x,kg1_adj,zaj_l(1,1,ik,2),kg1,1.d0,bWi_tmp,M)
            end if
         end if

         i = mod(il1,4)
         do iap = 1, n_ialist0
            if(i==2) then
               do ib = 1, np_e
                  fsrt = bWr_tmp(iap,ib); fsit = bWi_tmp(iap,ib)
                  bWr_tmp(iap,ib) = -fsit; bWi_tmp(iap,ib) = fsrt
               end do
            else if(i==3) then
               do ib = 1, np_e
                  fsrt = bWr_tmp(iap,ib); fsit = bWi_tmp(iap,ib)
                  bWr_tmp(iap,ib) = -fsrt; bWi_tmp(iap,ib) = -fsit
               end do
            else if(i==0) then
               do ib = 1, np_e
                  fsrt = bWr_tmp(iap,ib); fsit = bWi_tmp(iap,ib)
                  bWr_tmp(iap,ib) = fsit; bWi_tmp(iap,ib) = -fsrt
               end do
            end if
         end do

         do ib = 1, np_e
            do iap = 1, n_ialist0
               ia = n_ialist0*(lmt1-1) + iap
               bWr_lmt(ia,ib) = bWr_tmp(iap,ib)
               bWi_lmt(ia,ib) = bWi_tmp(iap,ib)
            end do
         end do

         do iap = 1, n_ialist0
            ia = ia_list(iap)
            lmta1 = lmta(lmt1,ia)
            do ib = 1, np_e
               bWr_noncl(lmta1,iter_rmm,ib,my_is) =  bWr_tmp(iap,ib)
               bWi_noncl(lmta1,iter_rmm,ib,my_is) =  bWi_tmp(iap,ib)
            end do
         end do

      end do

      deallocate(bWr_tmp, bWi_tmp)
      call tstatc0_end(id_sname)

    end subroutine sumset_rmm_all4_noncl

    subroutine Vnonlocal_W_part_sum_lmt4_noncl()
      integer             :: lmt1, lmtt1, il1, im1, il11, mdl, iksnl,il2,im2, &
           &                 ia_p, ia, i, lmt2,iap
      complex(kind=CMPLDP)       :: tmp

      real(kind=DP) :: rtmp1, rtmp2, atmp1, atmp2
      logical :: enter_flag

      integer :: my_is
      integer :: id_sname = -1

      call tstatc0_begin('Vnonlocal_W_part_sum_lmt4_noncl ',id_sname)

      sc_lmt = 0.d0; ss_lmt = 0.d0
      if(mdvdb == EXECUT) then
         qc_lmt = 0.d0; qs_lmt = 0.d0
      endif

      iksnl = (ik-1)/ndim_spinor + 1

      do lmt2 = 1, ilmt(it)
         il2   = ltp(lmt2,it)
         im2   = mtp(lmt2,it)

         do lmt1 = 1,ilmt(it)
            lmtt1 = lmtt(lmt1,it)
            il1   = ltp(lmt1,it)
            im1   = mtp(lmt1,it)
            il11  = il1 - 1
            mdl   = mod(il11,4)

! ---------------------------------------------------- 11.0S
#ifdef SKIP_TEST
            if ( mdvdb == SKIP ) then
               if ( il1 /= il2 ) then
                  fq(1:n_ialist0,1) = 0.0d0
                  goto 100
               endif
               if ( SpinOrbit_mode == Neglected .and. sw_hubbard == OFF ) then
                  if ( im1 /= im2 ) then
                     fq(1:n_ialist0,1) = 0.0d0
                     goto 100
                  endif
               endif
            endif
#endif
! ---------------------------------------------------- 11.0S

            do ia_p = 1, n_ialist0
               ia = ia_list(ia_p)

               tmp = dion_scr_noncl( lmt1, lmt2, ispin, ia )
               tmp = tmp *iwei(ia)

               if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
               fq(ia_p,1) = tmp
            end do

100         continue

! ----------------------------------------------------- 11.0S
!            if (mdvdb == EXECUT .and. il1 == il2 .and. im1 == im2) then
!
            enter_flag = .false.
            if ( mdvdb == EXECUT ) then
               if ( SpinOrbit_mode /= BuiltIn ) then
                  if ( il1 == il2 .and. im1 == im2 ) enter_flag = .true.
               else
                  if ( il1 == il2 ) enter_flag = .true.
               endif
            endif
!
            if ( enter_flag ) then
! ------------------------------------------------------ 11.0S

               do ia_p = 1, n_ialist0
                  ia = ia_list(ia_p)
                  tmp = q_noncl(lmt1,lmt2,ispin,it) *iwei(ia)
                  if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
                  fq(ia_p,2) = tmp
               end do

               if (mdl == 0 .or. mdl == 2) then

                  do ia_p = 1, n_ialist0
                     iap = n_ialist0*(lmt2-1) + ia_p

                     rtmp1 = real( fq(ia_p,1) );   atmp1 = aimag( fq(ia_p,1) )
                     rtmp2 = real( fq(ia_p,2) );   atmp2 = aimag( fq(ia_p,2) )

                     do i = 1, iba(ik)
                        sc_lmt(i,iap) = sc_lmt(i,iap) &
                             &   + (  rtmp1 *zfcos_x(i,ia_p) &
                             &       +atmp1 *zfsin_x(i,ia_p) ) *snl(i,lmtt1,iksnl)

                        ss_lmt(i,iap) = ss_lmt(i,iap) &
                             &   + ( -rtmp1 *zfsin_x(i,ia_p) &
                             &       +atmp1 *zfcos_x(i,ia_p) ) *snl(i,lmtt1,iksnl)


                        qc_lmt(i,iap) = qc_lmt(i,iap) &
                             &   + (  rtmp2 *zfcos_x(i,ia_p) &
                             &       +atmp2 *zfsin_x(i,ia_p) ) *snl(i,lmtt1,iksnl)


                        qs_lmt(i,iap) = qs_lmt(i,iap) &
                             &   + ( -rtmp2 *zfsin_x(i,ia_p) &
                             &       +atmp2 *zfcos_x(i,ia_p) ) *snl(i,lmtt1,iksnl)

                     end do
                  end do
               else if(mdl == 1 .or. mdl == 3) then
                  do ia_p = 1, n_ialist0
                     iap = n_ialist0*(lmt2-1) + ia_p

                     rtmp1 = real( fq(ia_p,1) );   atmp1 = aimag( fq(ia_p,1) )
                     rtmp2 = real( fq(ia_p,2) );   atmp2 = aimag( fq(ia_p,2) )

                     do i = 1, iba(ik)
                        sc_lmt(i,iap) = sc_lmt(i,iap) &
                             &   + ( -rtmp1 *zfsin_x(i,ia_p) &
                             &       +atmp1 *zfcos_x(i,ia_p) )*snl(i,lmtt1,iksnl)
                        ss_lmt(i,iap) = ss_lmt(i,iap) &
                             &   + ( -rtmp1 *zfcos_x(i,ia_p) &
                             &       -atmp1 *zfsin_x(i,ia_p) )*snl(i,lmtt1,iksnl)

                        qc_lmt(i,iap) = qc_lmt(i,iap) &
                             &   + ( -rtmp2 *zfsin_x(i,ia_p) &
                             &       +atmp2 *zfcos_x(i,ia_p) )*snl(i,lmtt1,iksnl)

                        qs_lmt(i,iap) = qs_lmt(i,iap) &
                             &   + ( -rtmp2 *zfcos_x(i,ia_p) &
                             &       -atmp2 *zfsin_x(i,ia_p) )*snl(i,lmtt1,iksnl)
                     end do
                  end do
               end if
            else
               if(mdl == 0 .or. mdl == 2) then
                  do ia_p = 1, n_ialist0
                     iap = n_ialist0*(lmt2-1) + ia_p
                     rtmp1 = real( fq(ia_p,1) );   atmp1 = aimag( fq(ia_p,1) )

                     do i = 1, iba(ik)
                        sc_lmt(i,iap) = sc_lmt(i,iap) &
                             &   + (  rtmp1 *zfcos_x(i,ia_p) &
                             &       +atmp1 *zfsin_x(i,ia_p) )*snl(i,lmtt1,iksnl)
                        ss_lmt(i,iap) = ss_lmt(i,iap) &
                             &   + ( -rtmp1 *zfsin_x(i,ia_p) &
                             &       +atmp1 *zfcos_x(i,ia_p) )*snl(i,lmtt1,iksnl)

                     end do
                  end do
               else if(mdl == 1 .or. mdl == 3) then
                  do ia_p = 1, n_ialist0
                     iap = n_ialist0*(lmt2-1) + ia_p
                     rtmp1 = real( fq(ia_p,1) );   atmp1 = aimag( fq(ia_p,1) )

                     do i = 1, iba(ik)
                        sc_lmt(i,iap) = sc_lmt(i,iap) &
                             &   + ( -rtmp1 *zfsin_x(i,ia_p) &
                             &       +atmp1 *zfcos_x(i,ia_p) )*snl(i,lmtt1,iksnl)
                        ss_lmt(i,iap) = ss_lmt(i,iap) &
                             &   + ( -rtmp1 *zfcos_x(i,ia_p) &
                             &       -atmp1 *zfsin_x(i,ia_p) )*snl(i,lmtt1,iksnl)

                     end do
                  end do
               end if
            end if
         end do
      end do

      call tstatc0_end(id_sname)
    end subroutine Vnonlocal_W_part_sum_lmt4_noncl

    subroutine add_vnlph_l_part4_noncl(mdvdb)
      integer, intent(in) :: mdvdb
      integer       :: ib, iap,  M, N, K
      real(kind=DP), allocatable, dimension(:,:) :: bWer, bWei

! ---------- ktDEBUG ------------------------------------ 20121001 --
      integer :: iksnl, my_is, ik_for_pointing_eko
! ---------- ktDEBUG ------------------------------------ 20121001 --

      integer :: id_sname = -1

      call tstatc0_begin('add_vnlph_l_part4_noncl ',id_sname)

      if(mdvdb == EXECUT) then
         allocate(bWer(ilmt(it)*n_ialist0,np_e))
         allocate(bWei(ilmt(it)*n_ialist0,np_e))

! ---------- ktDEBUG ------------------------------------ 20121001 --
!         do iap = 1, ilmt(it)*n_ialist0
!            do ib = 1, np_e
!               bWer(iap,ib) = bWr_lmt(iap,ib)*eko_l(ib,ik)
!               bWei(iap,ib) = bWi_lmt(iap,ib)*eko_l(ib,ik)
!            end do
!         end do
!
         iksnl = (ik-1)/ndim_spinor + 1
         my_is = mod( ik-1, ndim_spinor ) +1
         ik_for_pointing_eko = ( iksnl -1 )*ndim_spinor +1

         do iap = 1, ilmt(it)*n_ialist0
            do ib = 1, np_e
               bWer(iap,ib) = bWr_lmt(iap,ib)*eko_l(ib,ik_for_pointing_eko)
               bWei(iap,ib) = bWi_lmt(iap,ib)*eko_l(ib,ik_for_pointing_eko)
            end do
         end do
! ---------- ktDEBUG ------------------------------------ 20121001 --

      end if

      M = iba(ik); N = np_e ; K = ilmt(it)*n_ialist0
      if(kimg == 1) then
         call DGEMM__('N','N',M,N,K,  1.d0, sc_lmt,kg1_adj, bWr_lmt,K, 1.d0, vnlph_l(1,1,1),kg1)
         call DGEMM__('N','N',M,N,K, -1.d0, ss_lmt,kg1_adj, bWi_lmt,K, 1.d0, vnlph_l(1,1,1),kg1)

         if(mdvdb == EXECUT) then
            call DGEMM__('N','N',M,N,K, -1.d0, qc_lmt,kg1_adj, bWer,   K, 1.d0, vnlph_l(1,1,1),kg1)
            call DGEMM__('N','N',M,N,K,  1.d0, qs_lmt,kg1_adj, bWei,   K, 1.d0, vnlph_l(1,1,1),kg1)
         end if
      else if(kimg == 2) then
         call DGEMM__('N','N',M,N,K,  1.d0, sc_lmt,kg1_adj, bWr_lmt,K, 1.d0, vnlph_l(1,1,1),kg1)
         call DGEMM__('N','N',M,N,K, -1.d0, ss_lmt,kg1_adj, bWi_lmt,K, 1.d0, vnlph_l(1,1,1),kg1)
         call DGEMM__('N','N',M,N,K,  1.d0, sc_lmt,kg1_adj, bWi_lmt,K, 1.d0, vnlph_l(1,1,2),kg1)
         call DGEMM__('N','N',M,N,K,  1.d0, ss_lmt,kg1_adj, bWr_lmt,K, 1.d0, vnlph_l(1,1,2),kg1)

         if(mdvdb == EXECUT) then
            call DGEMM__('N','N',M,N,K, -1.d0, qc_lmt,kg1_adj, bWer,   K, 1.d0, vnlph_l(1,1,1),kg1)
            call DGEMM__('N','N',M,N,K,  1.d0, qs_lmt,kg1_adj, bWei,   K, 1.d0, vnlph_l(1,1,1),kg1)
            call DGEMM__('N','N',M,N,K, -1.d0, qc_lmt,kg1_adj, bWei,   K, 1.d0, vnlph_l(1,1,2),kg1)
            call DGEMM__('N','N',M,N,K, -1.d0, qs_lmt,kg1_adj, bWer,   K, 1.d0, vnlph_l(1,1,2),kg1)
         end if
      end if

      if(mdvdb == EXECUT) deallocate(bWer,bWei)
      call tstatc0_end(id_sname)
    end subroutine add_vnlph_l_part4_noncl
#else

    subroutine sumset_rmm_all3_noncl()

      integer  :: lmt1, lmtt1, lmta1, i, iksnl, ib, il1, ia, iap, il2
      real(kind=DP) :: fsrt,fsit

      integer :: my_is
      integer :: id_sname = -1

      call tstatc0_begin('sumset_rmm_all3_noncl ',id_sname)

      iksnl = (ik-1)/ndim_spinor + 1
      my_is = mod( ik-1, ndim_spinor ) +1

      do lmt1 = 1, ilmt(it)
         lmtt1 = lmtt(lmt1,it)
         il1   = ltp(lmt1,it)
         if(kimg == 1) then

#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
          do iap = 1, n_ialist0
                ia = ia_list(iap)
!               lmta1 = lmta(lmt1,ia)
                do i = 1, iba(ik)
                   ar_x(i,iap) = zfcos_x(i,iap)*snl(i,lmtt1,iksnl)
                   ai_x(i,iap) = zfsin_x(i,iap)*snl(i,lmtt1,iksnl)
                end do
             end do
          else if(kimg == 2) then
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do iap = 1,n_ialist0
                ia = ia_list(iap)
!               lmta1 = lmta(lmt1,ia)
                do i = 1, iba(ik)
                   ar_x(i,iap) = zfcos_x(i,iap)*snl(i,lmtt1,iksnl)
                   ai_x(i,iap) = zfsin_x(i,iap)*snl(i,lmtt1,iksnl)
                end do
             end do
          end if

          if(kimg == 1) then
#ifdef HIUX
*poption indep
*poption parallel
#endif
             do iap = 1, n_ialist0
                ia = ia_list(iap)
                lmta1 = lmta(lmt1,ia)
#ifdef VPP
*vocl loop, unroll(4)
#endif
                do ib = 1, np_e
                   fsrt = 0.d0; fsit = 0.d0
                   do i = 1, iba(ik)
                      fsrt = fsrt + ar_x(i,iap)*zaj_l(i,ib,ik,1)
                      fsit = fsit + ai_x(i,iap)*zaj_l(i,ib,ik,1)
                   end do
                   bWr_noncl(lmta1,iter_rmm,ib,my_is) = fsrt;
                   bWi_noncl(lmta1,iter_rmm,ib,my_is) = fsit
                end do
             end do
          else if(kimg == 2) then
             if(k_symmetry(ik) == GAMMA) then
                il2 = mod(il1,2) ! il2 = 1 -> l=0,2,4,6,...: il2 = 0 -> l=1,3,5,7,...
#ifdef HIUX
*poption indep
*poption parallel
#endif
                do iap = 1, n_ialist0
                   ia = ia_list(iap)
                   lmta1 = lmta(lmt1,ia)
#ifdef VPP
*vocl loop, unroll(4)
#endif
                   do ib = 1, np_e
                      fsrt = 0.d0; fsit = 0.d0
                      if(il2 == 1) then
                         do i = 2, iba(ik)
                            fsrt = fsrt + ar_x(i,iap)*zaj_l(i,ib,ik,1)-ai_x(i,iap)*zaj_l(i,ib,ik,2)
                         end do
                         fsrt = 2.d0*fsrt + zaj_l(1,ib,ik,1)*ar_x(1,iap)
                      else
                         do i = 2, iba(ik)
                            fsit = fsit + ai_x(i,iap)*zaj_l(i,ib,ik,1)+ar_x(i,iap)*zaj_l(i,ib,ik,2)
                         end do
                         fsit = 2.d0*fsit
                      end if
                      bWr_noncl(lmta1,iter_rmm,ib,my_is) = fsrt;
                      bWi_noncl(lmta1,iter_rmm,ib,my_is) = fsit
                   end do
                end do
             else
#ifdef HIUX
*poption indep
*poption parallel
#endif
                do iap = 1, n_ialist0
                   ia = ia_list(iap)
                   lmta1 = lmta(lmt1,ia)
#ifdef VPP
*vocl loop, unroll(4)
#endif
                   do ib = 1, np_e
                      fsrt = 0.d0; fsit = 0.d0
                      do i    = 1,iba(ik)
                         fsrt = fsrt + ar_x(i,iap)*zaj_l(i,ib,ik,1)-ai_x(i,iap)*zaj_l(i,ib,ik,2)
                         fsit = fsit + ai_x(i,iap)*zaj_l(i,ib,ik,1)+ar_x(i,iap)*zaj_l(i,ib,ik,2)
                      end do
                      bWr_noncl(lmta1,iter_rmm,ib,my_is) = fsrt;
                      bWi_noncl(lmta1,iter_rmm,ib,my_is) = fsit
                   end do
                end do
             end if
          end if

          i = mod(il1,4)
          do iap = 1, n_ialist0
             ia = ia_list(iap)
             lmta1 = lmta(lmt1,ia)
             if(i == 2) then
                do ib = 1, np_e
                   fsrt = bWr_noncl(lmta1,iter_rmm,ib,my_is);
                   fsit = bWi_noncl(lmta1,iter_rmm,ib,my_is)
                   bWr_noncl(lmta1,iter_rmm,ib,my_is) = -fsit;
                   bWi_noncl(lmta1,iter_rmm,ib,my_is) = fsrt
                end do
             else if(i == 3) then
                do ib = 1, np_e
                   fsrt = bWr_noncl(lmta1,iter_rmm,ib,my_is);
                   fsit = bWi_noncl(lmta1,iter_rmm,ib,my_is)
                   bWr_noncl(lmta1,iter_rmm,ib,my_is) = -fsrt;
                   bWi_noncl(lmta1,iter_rmm,ib,my_is) = -fsit
                end do
             else if(i == 0) then
                do ib = 1, np_e
                   fsrt = bWr_noncl(lmta1,iter_rmm,ib,my_is);
                   fsit = bWi_noncl(lmta1,iter_rmm,ib,my_is)
                   bWr_noncl(lmta1,iter_rmm,ib,my_is) =  fsit;
                   bWi_noncl(lmta1,iter_rmm,ib,my_is) = -fsrt
                end do
             end if
          end do
       end do
       call tstatc0_end(id_sname)
     end subroutine sumset_rmm_all3_noncl

     subroutine add_vnlph_l_no_eko_pt3_noncl()
        real(kind=DP) :: bWFr, bWFi
        integer       :: i, ib, ia, iap, lmta2

        integer :: my_is

        integer :: id_sname = -1

        my_is = mod( ik-1, ndim_spinor ) +1

        call tstatc0_begin('add_vnlph_l_no_eko_pt3_noncl ',id_sname)

        if(kimg == 1) then
!!$#ifdef HIUX
!!$*poption parallel
!!$#endif
           do iap = 1, n_ialist0
              ia = ia_list(iap)
              lmta2 = lmta(lmt2,ia)
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr = bWr_noncl(lmta2,iter_rmm,ib,my_is);
                bWFi = bWi_noncl(lmta2,iter_rmm,ib,my_is)
                do i = 1, iba(ik)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*sc_x(i,iap)-bWFi*ss_x(i,iap)
                enddo
             end do
          end do
       else if(kimg == 2) then
!!$#ifdef HIUX
!!$*poption parallel
!!$#endif
          do iap = 1, n_ialist0
             ia = ia_list(iap)
             lmta2 = lmta(lmt2,ia)
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr = bWr_noncl(lmta2,iter_rmm,ib,my_is);
                bWFi = bWi_noncl(lmta2,iter_rmm,ib,my_is)
                do i = 1, iba(ik)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*sc_x(i,iap)-bWFi*ss_x(i,iap)
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi*sc_x(i,iap)+bWFr*ss_x(i,iap)
                enddo
             end do
          end do
       endif
       call tstatc0_end(id_sname)
     end subroutine add_vnlph_l_no_eko_pt3_noncl

     subroutine add_vnlph_l_with_eko_pt3_noncl()
       real(kind=DP) :: bWFr, bWFi, e, sceqc, sseqs
       integer       :: i, ib, ia, iap, lmta2
#ifdef HIUX
       integer       :: lmta2_1, lmta2_2, lmta2_3, lmta2_4, lmta2_5, lmta2_6, lmta2_7, lmta2_8
       real(kind=DP) :: bWFr_1, bWFr_2, bWFr_3, bWFr_4, bWFr_5, bWFr_6, bWFr_7, bWFr_8
       real(kind=DP) :: bWFi_1, bWFi_2, bWFi_3, bWFi_4, bWFi_5, bWFi_6, bWFi_7, bWFi_8
       real(kind=DP) :: sceqc_1, sseqs_1, sceqc_2, sseqs_2, sceqc_3, sseqs_3, sceqc_4, sseqs_4 &
            &         , sceqc_5, sseqs_5, sceqc_6, sseqs_6, sceqc_7, sseqs_7, sceqc_8, sseqs_8
#endif

       integer :: my_is, iksnl, ik_for_pointing_eko
       integer :: id_sname = -1

       call tstatc0_begin('add_vnlph_l_with_eko_pt3_noncl ',id_sname)

       iksnl = (ik-1)/ndim_spinor + 1
       my_is = mod( ik-1, ndim_spinor ) +1
       ik_for_pointing_eko = ( iksnl -1 )*ndim_spinor +1

#ifdef HIUX
       if(kimg == 1) then
          if(n_ialist0 == 1) then
             lmta2 = lmta(lmt2,ia_list(1))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr = bWr_noncl(lmta2,iter_rmm,ib,my_is);
                bWFi = bWi_noncl(lmta2,iter_rmm,ib,my_is);
                e = eko_l(ib,ik_for_pointing_eko)

                do i = 1, iba(ik)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*(sc_x(i,1)-e*qc_x(i,1))-bWFi*(ss_x(i,1)-e*qs_x(i,1))
                enddo
             end do
          else if(n_ialist0 == 2) then
             lmta2_1 = lmta(lmt2,ia_list(1)); lmta2_2 = lmta(lmt2,ia_list(2))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib =1, np_e
                bWFr_1 = bWr_noncl(lmta2_1,iter_rmm,ib,my_is);
                bWFi_1 = bWi_noncl(lmta2_1,iter_rmm,ib.my_is)
                bWFr_2 = bWr_noncl(lmta2_2,iter_rmm,ib,my_is);
                bWFi_2 = bWi_noncl(lmta2_2,iter_rmm,ib,my_is);
                e = eko_l(ib,ik_for_pointing_eko)

                do i = 1, iba(ik)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*(sc_x(i,1)-e*qc_x(i,1))-bWFi_1*(ss_x(i,1)-e*qs_x(i,1)) &
                        &                            + bWFr_2*(sc_x(i,2)-e*qc_x(i,2))-bWFi_2*(ss_x(i,2)-e*qs_x(i,2))
                end do
             end do
          else if(n_ialist0 == 3) then
             lmta2_1 = lmta(lmt2,ia_list(1)); lmta2_2 = lmta(lmt2,ia_list(2)); lmta2_3 = lmta(lmt2,ia_list(3))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib =1, np_e
                bWFr_1 = bWr_noncl(lmta2_1,iter_rmm,ib,my_is);
                bWFi_1 = bWi_noncl(lmta2_1,iter_rmm,ib,my_is)
                bWFr_2 = bWr_noncl(lmta2_2,iter_rmm,ib,my_is);
                bWFi_2 = bWi_noncl(lmta2_2,iter_rmm,ib,my_is)
                bWFr_3 = bWr_noncl(lmta2_3,iter_rmm,ib,my_is);
                bWFi_3 = bWi_noncl(lmta2_3,iter_rmm,ib,my_is);
                e = eko_l(ib,ik_for_pointing_eko)

                do i = 1, iba(ik)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*(sc_x(i,1)-e*qc_x(i,1))-bWFi_1*(ss_x(i,1)-e*qs_x(i,1)) &
                        &                            + bWFr_2*(sc_x(i,2)-e*qc_x(i,2))-bWFi_2*(ss_x(i,2)-e*qs_x(i,2)) &
                        &                            + bWFr_3*(sc_x(i,3)-e*qc_x(i,3))-bWFi_3*(ss_x(i,3)-e*qs_x(i,3))
                end do
             end do
          else if(n_ialist0 >= 4) then
             lmta2_1=lmta(lmt2,ia_list(1)); lmta2_2=lmta(lmt2,ia_list(2)); lmta2_3=lmta(lmt2,ia_list(3)); lmta2_4=lmta(lmt2,ia_list(4))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib =1, np_e
                bWFr_1 = bWr_noncl(lmta2_1,iter_rmm,ib,my_is);
                bWFi_1 = bWi_noncl(lmta2_1,iter_rmm,ib,my_is)
                bWFr_2 = bWr_noncl(lmta2_2,iter_rmm,ib,my_is);
                bWFi_2 = bWi_noncl(lmta2_2,iter_rmm,ib,my_is)
                bWFr_3 = bWr_noncl(lmta2_3,iter_rmm,ib,my_is);
                bWFi_3 = bWi_noncl(lmta2_3,iter_rmm,ib,my_is)
                bWFr_4 = bWr_noncl(lmta2_4,iter_rmm,ib,my_is);
                bWFi_4 = bWi_noncl(lmta2_4,iter_rmm,ib,my_is);
                e = eko_l(ib,ik_for_pointing_eko)

                do i = 1, iba(ik)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*(sc_x(i,1)-e*qc_x(i,1))-bWFi_1*(ss_x(i,1)-e*qs_x(i,1)) &
                        &                            + bWFr_2*(sc_x(i,2)-e*qc_x(i,2))-bWFi_2*(ss_x(i,2)-e*qs_x(i,2)) &
                        &                            + bWFr_3*(sc_x(i,3)-e*qc_x(i,3))-bWFi_3*(ss_x(i,3)-e*qs_x(i,3)) &
                        &                            + bWFr_4*(sc_x(i,4)-e*qc_x(i,4))-bWFi_4*(ss_x(i,4)-e*qs_x(i,4))
                end do
             end do
          end if
          if(n_ialist0 > 4) then
             do iap = 5, n_ialist0
                ia = ia_list(iap)
                lmta2 = lmta(lmt2,ia)
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
                do ib = 1, np_e
                   bWFr = bWr_noncl(lmta2,iter_rmm,ib,my_is);
                   bWFi = bWi_noncl(lmta2,iter_rmm,ib,my_is);
                   e = eko_l(ib,ik_for_pointing_eko)
                   do i = 1, iba(ik)
                      vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*(sc_x(i,iap)-e*qc_x(i,iap))-bWFi*(ss_x(i,iap)-e*qs_x(i,iap))
                   enddo
                end do
             end do
          end if
       else if(kimg == 2) then
         if(n_ialist0 == 1) then
             lmta2 = lmta(lmt2,ia_list(1))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr = bWr_noncl(lmta2,iter_rmm,ib,my_is);
                bWFi = bWi_noncl(lmta2,iter_rmm,ib,my_is);
                e = eko_l(ib,ik_for_pointing_eko)
                do i = 1, iba(ik)
                   sceqc = sc_x(i,1)-e*qc_x(i,1);  sseqs = ss_x(i,1)-e*qs_x(i,1)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*sceqc-bWFi*sseqs
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi*sceqc+bWFr*sseqs
                end do
             end do
          else if(n_ialist0 == 2) then
             lmta2_1 = lmta(lmt2,ia_list(1)); lmta2_2 = lmta(lmt2,ia_list(2))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr_1 = bWr_noncl(lmta2_1,iter_rmm,ib,my_is);
                bWFi_1 = bWi_noncl(lmta2_1,iter_rmm,ib,my_is)
                bWFr_2 = bWr_noncl(lmta2_2,iter_rmm,ib,my_is);
                bWFi_2 = bWi_noncl(lmta2_2,iter_rmm,ib,my_is);
                e = eko_l(ib,ik_for_pointing_eko)
                do i = 1, iba(ik)
                   sceqc_1 = sc_x(i,1)-e*qc_x(i,1);  sseqs_1 = ss_x(i,1)-e*qs_x(i,1)
                   sceqc_2 = sc_x(i,2)-e*qc_x(i,2);  sseqs_2 = ss_x(i,2)-e*qs_x(i,2)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*sceqc_1-bWFi_1*sseqs_1 + bWFr_2*sceqc_2-bWFi_2*sseqs_2
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi_1*sceqc_1+bWFr_1*sseqs_1 + bWFi_2*sceqc_2+bWFr_2*sseqs_2
                end do
             end do
          else if(n_ialist0 == 3) then
             lmta2_1 = lmta(lmt2,ia_list(1)); lmta2_2 = lmta(lmt2,ia_list(2)); lmta2_3 = lmta(lmt2,ia_list(3))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr_1 = bWr_noncl(lmta2_1,iter_rmm,ib,my_is);
                bWFi_1 = bWi_noncl(lmta2_1,iter_rmm,ib,my_is)
                bWFr_2 = bWr_noncl(lmta2_2,iter_rmm,ib,my_is);
                bWFi_2 = bWi_noncl(lmta2_2,iter_rmm,ib,my_is)
                bWFr_3 = bWr_noncl(lmta2_3,iter_rmm,ib,my_is);
                bWFi_3 = bWi_noncl(lmta2_3,iter_rmm,ib,my_is);
                e = eko_l(ib,ik_for_pointing_eko)
                do i = 1, iba(ik)
                   sceqc_1 = sc_x(i,1)-e*qc_x(i,1);  sseqs_1 = ss_x(i,1)-e*qs_x(i,1)
                   sceqc_2 = sc_x(i,2)-e*qc_x(i,2);  sseqs_2 = ss_x(i,2)-e*qs_x(i,2)
                   sceqc_3 = sc_x(i,3)-e*qc_x(i,3);  sseqs_3 = ss_x(i,3)-e*qs_x(i,3)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*sceqc_1-bWFi_1*sseqs_1 + bWFr_2*sceqc_2-bWFi_2*sseqs_2 &
                        &                            + bWFr_3*sceqc_3-bWFi_3*sseqs_3
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi_1*sceqc_1+bWFr_1*sseqs_1 + bWFi_2*sceqc_2+bWFr_2*sseqs_2 &
                        &                            + bWFi_3*sceqc_3+bWFr_3*sseqs_3
                end do
             end do
          else if(n_ialist0 == 4) then
             lmta2_1=lmta(lmt2,ia_list(1)); lmta2_2=lmta(lmt2,ia_list(2)); lmta2_3=lmta(lmt2,ia_list(3));lmta2_4=lmta(lmt2,ia_list(4))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr_1 = bWr_noncl(lmta2_1,iter_rmm,ib,my_is);
                bWFi_1 = bWi_noncl(lmta2_1,iter_rmm,ib,my_is)
                bWFr_2 = bWr_noncl(lmta2_2,iter_rmm,ib,my_is);
                bWFi_2 = bWi_noncl(lmta2_2,iter_rmm,ib,my_is)
                bWFr_3 = bWr_noncl(lmta2_3,iter_rmm,ib,my_is);
                bWFi_3 = bWi_noncl(lmta2_3,iter_rmm,ib,my_is)
                bWFr_4 = bWr_noncl(lmta2_4,iter_rmm,ib,my_is);
                bWFi_4 = bWi_noncl(lmta2_4,iter_rmm,ib,my_is);
                e = eko_l(ib,ik_for_pointing_eko)
                do i = 1, iba(ik)
                   sceqc_1 = sc_x(i,1)-e*qc_x(i,1);  sseqs_1 = ss_x(i,1)-e*qs_x(i,1)
                   sceqc_2 = sc_x(i,2)-e*qc_x(i,2);  sseqs_2 = ss_x(i,2)-e*qs_x(i,2)
                   sceqc_3 = sc_x(i,3)-e*qc_x(i,3);  sseqs_3 = ss_x(i,3)-e*qs_x(i,3)
                   sceqc_4 = sc_x(i,4)-e*qc_x(i,4);  sseqs_4 = ss_x(i,4)-e*qs_x(i,4)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*sceqc_1-bWFi_1*sseqs_1 + bWFr_2*sceqc_2-bWFi_2*sseqs_2 &
                        &                            + bWFr_3*sceqc_3-bWFi_3*sseqs_3 + bWFr_4*sceqc_4-bWFi_4*sseqs_4
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi_1*sceqc_1+bWFr_1*sseqs_1 + bWFi_2*sceqc_2+bWFr_2*sseqs_2 &
                        &                            + bWFi_3*sceqc_3+bWFr_3*sseqs_3 + bWFi_4*sceqc_4+bWFr_4*sseqs_4
                end do
             end do
          else if(n_ialist0 == 5) then
             lmta2_1=lmta(lmt2,ia_list(1)); lmta2_2=lmta(lmt2,ia_list(2)); lmta2_3=lmta(lmt2,ia_list(3));lmta2_4=lmta(lmt2,ia_list(4))
             lmta2_5=lmta(lmt2,ia_list(5))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr_1 = bWr_noncl(lmta2_1,iter_rmm,ib,my_is);
                bWFi_1 = bWi_noncl(lmta2_1,iter_rmm,ib,my_is)
                bWFr_2 = bWr_noncl(lmta2_2,iter_rmm,ib,my_is);
                bWFi_2 = bWi_noncl(lmta2_2,iter_rmm,ib,my_is)
                bWFr_3 = bWr_noncl(lmta2_3,iter_rmm,ib,my_is);
                bWFi_3 = bWi_noncl(lmta2_3,iter_rmm,ib,my_is)
                bWFr_4 = bWr_noncl(lmta2_4,iter_rmm,ib,my_is);
                bWFi_4 = bWi_noncl(lmta2_4,iter_rmm,ib,my_is)
                bWFr_5 = bWr_noncl(lmta2_5,iter_rmm,ib,my_is);
                bWFi_5 = bWi_noncl(lmta2_5,iter_rmm,ib,my_is);
                e = eko_l(ib,ik_for_pointing_eko)

                do i = 1, iba(ik)
                   sceqc_1 = sc_x(i,1)-e*qc_x(i,1);  sseqs_1 = ss_x(i,1)-e*qs_x(i,1)
                   sceqc_2 = sc_x(i,2)-e*qc_x(i,2);  sseqs_2 = ss_x(i,2)-e*qs_x(i,2)
                   sceqc_3 = sc_x(i,3)-e*qc_x(i,3);  sseqs_3 = ss_x(i,3)-e*qs_x(i,3)
                   sceqc_4 = sc_x(i,4)-e*qc_x(i,4);  sseqs_4 = ss_x(i,4)-e*qs_x(i,4)
                   sceqc_5 = sc_x(i,5)-e*qc_x(i,5);  sseqs_5 = ss_x(i,5)-e*qs_x(i,5)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*sceqc_1-bWFi_1*sseqs_1 + bWFr_2*sceqc_2-bWFi_2*sseqs_2 &
                       &                            + bWFr_3*sceqc_3-bWFi_3*sseqs_3 + bWFr_4*sceqc_4-bWFi_4*sseqs_4 &
                       &                            + bWFr_5*sceqc_5-bWFi_5*sseqs_5
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi_1*sceqc_1+bWFr_1*sseqs_1 + bWFi_2*sceqc_2+bWFr_2*sseqs_2 &
                       &                            + bWFi_3*sceqc_3+bWFr_3*sseqs_3 + bWFi_4*sceqc_4+bWFr_4*sseqs_4 &
                       &                            + bWFi_5*sceqc_5+bWFr_3*sseqs_5
                end do
             end do
          else if(n_ialist0 == 6) then
             lmta2_1=lmta(lmt2,ia_list(1)); lmta2_2=lmta(lmt2,ia_list(2)); lmta2_3=lmta(lmt2,ia_list(3));lmta2_4=lmta(lmt2,ia_list(4))
             lmta2_5=lmta(lmt2,ia_list(5)); lmta2_6=lmta(lmt2,ia_list(6))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr_1 = bWr_noncl(lmta2_1,iter_rmm,ib,my_is);
                bWFi_1 = bWi_noncl(lmta2_1,iter_rmm,ib,my_is)
                bWFr_2 = bWr_noncl(lmta2_2,iter_rmm,ib,my_is);
                bWFi_2 = bWi_noncl(lmta2_2,iter_rmm,ib,my_is)
                bWFr_3 = bWr_noncl(lmta2_3,iter_rmm,ib,my_is);
                bWFi_3 = bWi_noncl(lmta2_3,iter_rmm,ib,my_is)
                bWFr_4 = bWr_noncl(lmta2_4,iter_rmm,ib,my_is);
                bWFi_4 = bWi_noncl(lmta2_4,iter_rmm,ib,my_is)
                bWFr_5 = bWr_noncl(lmta2_5,iter_rmm,ib,my_is);
                bWFi_5 = bWi_noncl(lmta2_5,iter_rmm,ib,my_is)
                bWFr_6 = bWr_noncl(lmta2_6,iter_rmm,ib,my_is);
                bWFi_6 = bWi_noncl(lmta2_6,iter_rmm,ib,my_is);
                e = eko_l(ib,ik_for_pointing_eko)

                do i = 1, iba(ik)
                   sceqc_1 = sc_x(i,1)-e*qc_x(i,1);  sseqs_1 = ss_x(i,1)-e*qs_x(i,1)
                   sceqc_2 = sc_x(i,2)-e*qc_x(i,2);  sseqs_2 = ss_x(i,2)-e*qs_x(i,2)
                   sceqc_3 = sc_x(i,3)-e*qc_x(i,3);  sseqs_3 = ss_x(i,3)-e*qs_x(i,3)
                   sceqc_4 = sc_x(i,4)-e*qc_x(i,4);  sseqs_4 = ss_x(i,4)-e*qs_x(i,4)
                   sceqc_5 = sc_x(i,5)-e*qc_x(i,5);  sseqs_5 = ss_x(i,5)-e*qs_x(i,5)
                   sceqc_6 = sc_x(i,6)-e*qc_x(i,6);  sseqs_6 = ss_x(i,6)-e*qs_x(i,6)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*sceqc_1-bWFi_1*sseqs_1 + bWFr_2*sceqc_2-bWFi_2*sseqs_2 &
                        &                            + bWFr_3*sceqc_3-bWFi_3*sseqs_3 + bWFr_4*sceqc_4-bWFi_4*sseqs_4 &
                        &                            + bWFr_5*sceqc_5-bWFi_5*sseqs_5 + bWFr_6*sceqc_6-bWFi_6*sseqs_6
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi_1*sceqc_1+bWFr_1*sseqs_1 + bWFi_2*sceqc_2+bWFr_2*sseqs_2 &
                        &                            + bWFi_3*sceqc_3+bWFr_3*sseqs_3 + bWFi_4*sceqc_4+bWFr_4*sseqs_4 &
                        &                            + bWFi_5*sceqc_5+bWFr_5*sseqs_5 + bWFi_6*sceqc_6+bWFr_6*sseqs_6
                end do
             end do
          else if(n_ialist0 == 7) then
             lmta2_1=lmta(lmt2,ia_list(1)); lmta2_2=lmta(lmt2,ia_list(2)); lmta2_3=lmta(lmt2,ia_list(3));lmta2_4=lmta(lmt2,ia_list(4))
             lmta2_5=lmta(lmt2,ia_list(5)); lmta2_6=lmta(lmt2,ia_list(6)); lmta2_7=lmta(lmt2,ia_list(7))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr_1 = bWr_noncl(lmta2_1,iter_rmm,ib,my_is);
                bWFi_1 = bWi_noncl(lmta2_1,iter_rmm,ib,my_is)
                bWFr_2 = bWr_noncl(lmta2_2,iter_rmm,ib,my_is);
                bWFi_2 = bWi_noncl(lmta2_2,iter_rmm,ib,my_is)
                bWFr_3 = bWr_noncl(lmta2_3,iter_rmm,ib,my_is);
                bWFi_3 = bWi_noncl(lmta2_3,iter_rmm,ib,my_is)
                bWFr_4 = bWr_noncl(lmta2_4,iter_rmm,ib,my_is);
                bWFi_4 = bWi_noncl(lmta2_4,iter_rmm,ib,my_is)
                bWFr_5 = bWr_noncl(lmta2_5,iter_rmm,ib,my_is);
                bWFi_5 = bWi_noncl(lmta2_5,iter_rmm,ib,my_is)
                bWFr_6 = bWr_noncl(lmta2_6,iter_rmm,ib,my_is);
                bWFi_6 = bWi_noncl(lmta2_6,iter_rmm,ib,my_is)
                bWFr_7 = bWr_noncl(lmta2_7,iter_rmm,ib,my_is);
                bWFi_7 = bWi_noncl(lmta2_7,iter_rmm,ib,my_is);
                e = eko_l(ib,ik_for_pointing_eko)

                do i = 1, iba(ik)
                   sceqc_1 = sc_x(i,1)-e*qc_x(i,1);  sseqs_1 = ss_x(i,1)-e*qs_x(i,1)
                   sceqc_2 = sc_x(i,2)-e*qc_x(i,2);  sseqs_2 = ss_x(i,2)-e*qs_x(i,2)
                   sceqc_3 = sc_x(i,3)-e*qc_x(i,3);  sseqs_3 = ss_x(i,3)-e*qs_x(i,3)
                   sceqc_4 = sc_x(i,4)-e*qc_x(i,4);  sseqs_4 = ss_x(i,4)-e*qs_x(i,4)
                   sceqc_5 = sc_x(i,5)-e*qc_x(i,5);  sseqs_5 = ss_x(i,5)-e*qs_x(i,5)
                   sceqc_6 = sc_x(i,6)-e*qc_x(i,6);  sseqs_6 = ss_x(i,6)-e*qs_x(i,6)
                   sceqc_7 = sc_x(i,7)-e*qc_x(i,7);  sseqs_7 = ss_x(i,7)-e*qs_x(i,7)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*sceqc_1-bWFi_1*sseqs_1 + bWFr_2*sceqc_2-bWFi_2*sseqs_2 &
                        &                            + bWFr_3*sceqc_3-bWFi_3*sseqs_3 + bWFr_4*sceqc_4-bWFi_4*sseqs_4 &
                        &                            + bWFr_5*sceqc_5-bWFi_5*sseqs_5 + bWFr_6*sceqc_6-bWFi_6*sseqs_6 &
                        &                            + bWFr_7*sceqc_7-bWFi_7*sseqs_7
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi_1*sceqc_1+bWFr_1*sseqs_1 + bWFi_2*sceqc_2+bWFr_2*sseqs_2 &
                        &                            + bWFi_3*sceqc_3+bWFr_3*sseqs_3 + bWFi_4*sceqc_4+bWFr_4*sseqs_4 &
                        &                            + bWFi_5*sceqc_5+bWFr_5*sseqs_5 + bWFi_6*sceqc_6+bWFr_6*sseqs_6 &
                        &                            + bWFi_7*sceqc_7+bWFr_7*sseqs_7
                end do
             end do
          else if(n_ialist0 >= 8) then
             lmta2_1=lmta(lmt2,ia_list(1)); lmta2_2=lmta(lmt2,ia_list(2)); lmta2_3=lmta(lmt2,ia_list(3));lmta2_4=lmta(lmt2,ia_list(4))
             lmta2_5=lmta(lmt2,ia_list(5)); lmta2_6=lmta(lmt2,ia_list(6)); lmta2_7=lmta(lmt2,ia_list(7));lmta2_8=lmta(lmt2,ia_list(8))
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr_1 = bWr_noncl(lmta2_1,iter_rmm,ib,my_is);
                bWFi_1 = bWi_noncl(lmta2_1,iter_rmm,ib,my_is)
                bWFr_2 = bWr_noncl(lmta2_2,iter_rmm,ib,my_is);
                bWFi_2 = bWi_noncl(lmta2_2,iter_rmm,ib,my_is)
                bWFr_3 = bWr_noncl(lmta2_3,iter_rmm,ib,my_is);
                bWFi_3 = bWi_noncl(lmta2_3,iter_rmm,ib,my_is)
                bWFr_4 = bWr_noncl(lmta2_4,iter_rmm,ib,my_is);
                bWFi_4 = bWi_noncl(lmta2_4,iter_rmm,ib,my_is)
                bWFr_5 = bWr_noncl(lmta2_5,iter_rmm,ib,my_is);
                bWFi_5 = bWi_noncl(lmta2_5,iter_rmm,ib,my_is)
                bWFr_6 = bWr_noncl(lmta2_6,iter_rmm,ib,my_is);
                bWFi_6 = bWi_noncl(lmta2_6,iter_rmm,ib,my_is)
                bWFr_7 = bWr_noncl(lmta2_7,iter_rmm,ib,my_is);
                bWFi_7 = bWi_noncl(lmta2_7,iter_rmm,ib,my_is)
                bWFr_8 = bWr_noncl(lmta2_8,iter_rmm,ib,my_is);
                bWFi_8 = bWi_noncl(lmta2_8,iter_rmm,ib,my_is);
                e = eko_l(ib,ik_for_pointing_eko)

                do i = 1, iba(ik)
                   sceqc_1 = sc_x(i,1)-e*qc_x(i,1);  sseqs_1 = ss_x(i,1)-e*qs_x(i,1)
                   sceqc_2 = sc_x(i,2)-e*qc_x(i,2);  sseqs_2 = ss_x(i,2)-e*qs_x(i,2)
                   sceqc_3 = sc_x(i,3)-e*qc_x(i,3);  sseqs_3 = ss_x(i,3)-e*qs_x(i,3)
                   sceqc_4 = sc_x(i,4)-e*qc_x(i,4);  sseqs_4 = ss_x(i,4)-e*qs_x(i,4)
                   sceqc_5 = sc_x(i,5)-e*qc_x(i,5);  sseqs_5 = ss_x(i,5)-e*qs_x(i,5)
                   sceqc_6 = sc_x(i,6)-e*qc_x(i,6);  sseqs_6 = ss_x(i,6)-e*qs_x(i,6)
                   sceqc_7 = sc_x(i,7)-e*qc_x(i,7);  sseqs_7 = ss_x(i,7)-e*qs_x(i,7)
                   sceqc_8 = sc_x(i,8)-e*qc_x(i,8);  sseqs_8 = ss_x(i,8)-e*qs_x(i,8)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr_1*sceqc_1-bWFi_1*sseqs_1 + bWFr_2*sceqc_2-bWFi_2*sseqs_2 &
                        &                            + bWFr_3*sceqc_3-bWFi_3*sseqs_3 + bWFr_4*sceqc_4-bWFi_4*sseqs_4 &
                        &                            + bWFr_5*sceqc_5-bWFi_5*sseqs_5 + bWFr_6*sceqc_6-bWFi_6*sseqs_6 &
                        &                            + bWFr_7*sceqc_7-bWFi_7*sseqs_7 + bWFr_8*sceqc_8-bWFi_8*sseqs_8
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi_1*sceqc_1+bWFr_1*sseqs_1 + bWFi_2*sceqc_2+bWFr_2*sseqs_2 &
                        &                            + bWFi_3*sceqc_3+bWFr_3*sseqs_3 + bWFi_4*sceqc_4+bWFr_4*sseqs_4 &
                        &                            + bWFi_5*sceqc_5+bWFr_5*sseqs_5 + bWFi_6*sceqc_6+bWFr_6*sseqs_6 &
                        &                            + bWFi_7*sceqc_7+bWFr_7*sseqs_7 + bWFi_8*sceqc_8+bWFr_8*sseqs_8
                end do
             end do
          end if
          if(n_ialist0 > 8) then
             do iap = 9, n_ialist0
                ia = ia_list(iap)
                lmta2 = lmta(lmt2,ia)
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
                do ib = 1, np_e
                   bWFr = bWr_noncl(lmta2,iter_rmm,ib,my_is);
                   bWFi = bWi_noncl(lmta2,iter_rmm,ib,my_is);
                   e = eko_l(ib,ik_for_pointing_eko)

                   do i = 1, iba(ik)
                      sceqc = sc_x(i,iap)-e*qc_x(i,iap)
                      sseqs = ss_x(i,iap)-e*qs_x(i,iap)
                      vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*sceqc-bWFi*sseqs
                      vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi*sceqc+bWFr*sseqs
                   end do
                enddo
            end do
         end if
       end if

#else
!  --> for ordinary scaler and vector machines but HIUX
!  = #ifndef HIUX

       if(kimg == 1) then
          do iap = 1, n_ialist0
             ia = ia_list(iap)
             lmta2 = lmta(lmt2,ia)
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr = bWr_noncl(lmta2,iter_rmm,ib,my_is);
                bWFi = bWi_noncl(lmta2,iter_rmm,ib,my_is);
                e = eko_l(ib,ik_for_pointing_eko)

                do i = 1, iba(ik)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*(sc_x(i,iap)-e*qc_x(i,iap))-bWFi*(ss_x(i,iap)-e*qs_x(i,iap))
                enddo
             end do
          end do
       else if(kimg == 2) then
          do iap = 1, n_ialist0
             ia = ia_list(iap)
             lmta2 = lmta(lmt2,ia)
#ifdef VPP
*vocl loop, unroll(4)
#endif
             do ib = 1, np_e
                bWFr = bWr_noncl(lmta2,iter_rmm,ib,my_is);
                bWFi = bWi_noncl(lmta2,iter_rmm,ib,my_is);
                e = eko_l(ib,ik_for_pointing_eko)

                do i = 1, iba(ik)
                   sceqc = sc_x(i,iap)-e*qc_x(i,iap)
                   sseqs = ss_x(i,iap)-e*qs_x(i,iap)
                   vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + bWFr*sceqc-bWFi*sseqs
                   vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + bWFi*sceqc+bWFr*sseqs
                end do
             enddo
          end do
       endif
#endif
       call tstatc0_end(id_sname)
     end subroutine add_vnlph_l_with_eko_pt3_noncl

     subroutine Vnonlocal_W_partsum_lmt1b_noncl()
       integer             :: lmt1, lmtt1, il1, im1, il11, mdl, iksnl,il2,im2, &
            &                 ia_p, ia, i, ii
       complex(kind=CMPLDP)       :: tmp
       real(kind=DP) :: rtmp1, rtmp2, atmp1, atmp2
       logical :: enter_flag

       integer :: id_sname = -1
       call tstatc0_begin('Vnonlocal_W_partsum_lmt1b_noncl ',id_sname)

       iksnl = (ik-1)/ndim_spinor + 1
       il2   = ltp(lmt2,it)
       im2   = mtp(lmt2,it)

       sc_x = 0.d0; ss_x = 0.d0
       if(mdvdb == EXECUT) then
          qc_x = 0.d0; qs_x = 0.d0
       endif

       do lmt1 = 1,ilmt(it)
          lmtt1 = lmtt(lmt1,it)
          il1   = ltp(lmt1,it)
          im1   = mtp(lmt1,it)
          il11  = il1 - 1
          mdl   = mod(il11,4)


! ---------------------------------------------------- 11.0S
#ifdef SKIP_TEST
          if ( mdvdb == SKIP ) then
             if ( il1 /= il2 ) then
                fq(1:n_ialist0,1) = 0.0d0
                goto 100
             endif
             if ( SpinOrbit_mode == Neglected .and. sw_hubbard == OFF ) then
                if ( im1 /= im2 ) then
                   fq(1:n_ialist0,1) = 0.0d0
                   goto 100
                endif
             endif
          endif
#endif
! ---------------------------------------------------- 11.0S

          do ia_p = 1, n_ialist0
             ia = ia_list(ia_p)
             tmp = dion_scr_noncl( lmt1, lmt2, ispin, ia )
             if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
             fq(ia_p,1) = tmp
          end do

100       continue

! ----------------------------------------------------- 11.0S
!          if (mdvdb == EXECUT .and. il1 == il2 .and. im1 == im2) then
!
          enter_flag = .false.
          if ( mdvdb == EXECUT ) then
             if ( SpinOrbit_mode /= BuiltIn ) then
                if ( il1 == il2 .and. im1 == im2 ) enter_flag = .true.
             else
                if ( il1 == il2 ) enter_flag = .true.
             endif
          endif
!
          if ( enter_flag ) then
! ------------------------------------------------------ 11.0S

             do ia_p = 1, n_ialist0
                ia = ia_list(ia_p)
                tmp = q_noncl(lmt1,lmt2,ispin,it)*iwei(ia)
                if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
                fq(ia_p,2) = tmp
             end do

             if(mdl == 0 .or. mdl == 2) then
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
                do ia_p = 1, n_ialist0
                   rtmp1 = real( fq(ia_p,1) );   atmp1 = aimag( fq(ia_p,1) )
                   rtmp2 = real( fq(ia_p,2) );   atmp2 = aimag( fq(ia_p,2) )

                   do i = 1, iba(ik)
                      sc_x(i,ia_p) = sc_x(i,ia_p) &
                           &     + (  rtmp1 *zfcos_x(i,ia_p) &
                           &         +atmp1 *zfsin_x(i,ia_p) ) *snl(i,lmtt1,iksnl)
                      ss_x(i,ia_p) = ss_x(i,ia_p) &
                           &     + ( -rtmp1 *zfsin_x(i,ia_p) &
                           &         +atmp1 *zfcos_x(i,ia_p) ) *snl(i,lmtt1,iksnl)
                      qc_x(i,ia_p) = qc_x(i,ia_p) &
                           &     + (  rtmp2 *zfcos_x(i,ia_p) &
                           &         +atmp2 *zfsin_x(i,ia_p) ) *snl(i,lmtt1,iksnl)
                      qs_x(i,ia_p) = qs_x(i,ia_p) &
                           &     + ( -rtmp2 *zfsin_x(i,ia_p) &
                           &         +atmp2 *zfcos_x(i,ia_p) ) *snl(i,lmtt1,iksnl)
                   end do
                end do
             else if(mdl == 1 .or. mdl == 3) then
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
                do ia_p = 1, n_ialist0
                   rtmp1 = real( fq(ia_p,1) );   atmp1 = aimag( fq(ia_p,1) )
                   rtmp2 = real( fq(ia_p,2) );   atmp2 = aimag( fq(ia_p,2) )

                   do i = 1, iba(ik)
                      sc_x(i,ia_p) = sc_x(i,ia_p) &
                           &     + ( -rtmp1 *zfsin_x(i,ia_p) &
                           &         +atmp1 *zfcos_x(i,ia_p) )*snl(i,lmtt1,iksnl)
                      ss_x(i,ia_p) = ss_x(i,ia_p) &
                           &     + ( -rtmp1 *zfcos_x(i,ia_p) &
                           &         -atmp1 *zfsin_x(i,ia_p) )*snl(i,lmtt1,iksnl)
                      qc_x(i,ia_p) = qc_x(i,ia_p) &
                           &     + ( -rtmp2 *zfsin_x(i,ia_p) &
                           &         +atmp2 *zfcos_x(i,ia_p) )*snl(i,lmtt1,iksnl)
                      qs_x(i,ia_p) = qs_x(i,ia_p) &
                           &     + ( -rtmp2 *zfcos_x(i,ia_p) &
                           &         -atmp2 *zfsin_x(i,ia_p) )*snl(i,lmtt1,iksnl)
                   end do
                end do
             end if
          else
             if(mdl == 0 .or. mdl == 2) then
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
                do ia_p = 1, n_ialist0
                   rtmp1 = real( fq(ia_p,1) );   atmp1 = aimag( fq(ia_p,1) )

                   do i = 1, iba(ik)
                      sc_x(i,ia_p) = sc_x(i,ia_p) &
                           &     + (  rtmp1 *zfcos_x(i,ia_p) &
                           &         +atmp1 *zfsin_x(i,ia_p) )*snl(i,lmtt1,iksnl)
                      ss_x(i,ia_p) = ss_x(i,ia_p) &
                           &     + ( -rtmp1 *zfsin_x(i,ia_p) &
                           &         +atmp1 *zfcos_x(i,ia_p) )*snl(i,lmtt1,iksnl)
                   end do
                end do
             else if(mdl == 1 .or. mdl == 3) then
#ifdef HIUX
*poption parallel
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
                do ia_p = 1, n_ialist0
                   rtmp1 = real( fq(ia_p,1) );   atmp1 = aimag( fq(ia_p,1) )

                   do i = 1, iba(ik)
                      sc_x(i,ia_p) = sc_x(i,ia_p) &
                           &     + ( -rtmp1 *zfsin_x(i,ia_p) &
                           &         +atmp1 *zfcos_x(i,ia_p) )*snl(i,lmtt1,iksnl)
                      ss_x(i,ia_p) = ss_x(i,ia_p) &
                           &     + ( -rtmp1 *zfcos_x(i,ia_p) &
                           &         -atmp1 *zfsin_x(i,ia_p) )*snl(i,lmtt1,iksnl)

                   end do
                end do
             end if
          end if
       end do

       call tstatc0_end(id_sname)
     end subroutine Vnonlocal_W_partsum_lmt1b_noncl
#endif
  end subroutine Vnonlocal_W_RMMn_noncl
! =========================================================================== 11.0

end module m_ES_WF_by_RMM
