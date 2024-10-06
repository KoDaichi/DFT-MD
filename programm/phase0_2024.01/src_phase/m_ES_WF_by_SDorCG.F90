!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 593 $)
!
!  MODULE:  m_ES_WF_by_SDorCG
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!  FURTHER MODIFICATION: T. Yamasaki, January/13/2004
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
#undef NEC_TIMER
#ifdef NEC_TIMER
#  define START_TIMER(a) call start_timer(a)
#  define STOP_TIMER(a)  call stop_timer(a)
#else
#  define START_TIMER(a)
#  define STOP_TIMER(a)
#endif
!
#ifdef __TIMER__
#   define __TIMER_START_w_BARRIER(str,a)   call mpi_barrier(str,ierr) ;   call timer_sta(a)
#   define __TIMER_START(a)       call timer_sta(a)
#   define __TIMER_STOP(a)        call timer_end(a)
#else
#   define __TIMER_START_w_BARRIER(str,a)
#   define __TIMER_START(a)
#   define __TIMER_STOP(a)
#endif
#ifdef __TIMER_SUB__
#   define __TIMER_SUB_START(a)  call timer_sta(a)
#   define __TIMER_SUB_STOP(a)   call timer_end(a)
#else
#   define __TIMER_SUB_START(a)
#   define __TIMER_SUB_STOP(a)
#endif
#ifdef __TIMER_DO__
#   define __TIMER_DO_START(a)   call timer_sta(a)
#   define __TIMER_DO_STOP(a)    call timer_end(a)
#else
#   define __TIMER_DO_START(a)
#   define __TIMER_DO_STOP(a)
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
module m_ES_WF_by_SDorCG
!     (m_ESsd or m_ESSD)
! $Id: m_ES_WF_by_SDorCG.F90 593 2019-06-20 03:47:31Z jkoga $
! This module evolves wave functions.
! MSD: modified steepest descent algorithm (Original subroutine name was
! "msdv").
! Following lines are comments in the "msdv".
!
!!$C ref.: Williams and Soler, Bull.Am.Phys.Soc.32 ('87),p562.
!!$c
!!$c      revised by T.Yamasaki
!!$c         #1) 2 folding inline expansion at an do loop of NEG is
!!$c            rewrited as before as a simple do loop
!!$c                                               5th May 1994
!!$c         #2) Parallelization                   9th May 1994
!!$c         #3) Parallelization transient
!!$c                         by T.Yamasaki        15th May 1994
!!$c
!!$c         #4)  Following arrays are removed from the arguments.
!!$c                                  ----> 31th May '94 Y.M
!!$cc   >        VXC,CHGQ,PSC,ZFM3,
!!$c                                  ----> 3rd June '94, Y.M.
!!$cc   W        W31,
!!$c                                  ----> 3rd June '94, Y.M.
!!$cc   W        W310,
!!$c
!!$c         #5)  data array size is reduced by T.Yamasaki
!!$c             snlfmc(kng1,klmt),snlfms(kng1,klmt)
!!$c          -->  snlfmc(kng1), snlfms(kng1)       on 17th Jun. 1994
!!$c
!!$c         #6) Data parallelization               on 24th Jun. 1994
!!$c                 eko --> eko_l             by T.Yamasaki
!!$c         #7) indx0 is removed                      26th Jun. 1994
!!$c         #8) Interface has changed                 27th Jun. 1994
!!$c         #9) imsd is added in the argument list     4th Aug. 1994
!!$cc      imsd = 1
!!$cc    for preconditioning IMSD=2   '94.6.14 Tsuyoshi Miyazaki
!!$cc      imsd = 2
!!$c         #10) ivanl is inserted in the arugment list  13th Sep. 1994
!!$c         #11) A do loop of 450 is tuned up.        15th Sep. 1994
!!$c         #12) Spin polarization is introduced.     11th Dec. 1994
!!$c         #13) vlhxc --> vlhxc_l                    31st Jan. 1995
!!$c                                          by T.Yamasaki
!!$c         #14) gx,gy,gz --> ngabc, vx,vy,vz --> vkxyz
!!$c                   by T.Yamasaki                on 15th Feb. 1995
!!$c         #15) antiferromagnetic calculation is added on 9th Jul. 1996
!!$c                                          by H.Sawada
!
! The "msdv" subroutine was translated into a subroutine,
! "evolve_WFs_in_MSD_direction" by T. Yamasaki in 1999.
!
  use m_Control_Parameters,  only : nspin,ipri,iprisolver,kimg,neg,af,sw_hybrid_functional, sw_fef &
#ifdef SAVE_FFT_TIMES
 &                                , sw_modified_cg_formula, sw_save_fft, sw_retard_eigval_evaluation
#else
 &                                , sw_modified_cg_formula, sw_retard_eigval_evaluation
#endif
  use m_Const_Parameters,    only : DP,SP,DIRECT,ON,SD,MSD,eazyCG,CG,SKIP,EXECUT,zi &
 &                                , ORTHONORMALIZATION, NORMALIZATION, ELECTRON &
 &                                , OTHER_BANDS, SAME_BAND, ALL_BANDS, SCF, OFF, GAMMA, OLD
  use m_Crystal_Structure,   only : rltv
  use m_Electronic_Structure,only : zaj_l, afft,bfft &
 &                                , eko_l,vnlph_l,vlhxcQ &
 &                                , vlhxc_l, vtau_l &
#ifdef SAVE_FFT_TIMES
 &                                , status_saved_phifftr &
#endif
 &                                , m_ES_WF_in_Rspace &
 &                                , m_ES_sort_eigen_values &
 &                                , m_ES_Vlocal_in_Rspace  &
 &                                , m_ES_eigen_values_for_each_k  &
 &                                , m_ES_wd_zaj_small_portion &
 &                                , m_ES_wd_eko &
 &                                , m_ES_decide_precon_factor     &
 &                                , m_ES_sum_of_LocalPart  &
 &                                , m_ES_sum_of_LocalPart2 &
 &                                , m_ES_sum_of_LocalPart3

  use m_ES_ortho            ,only : m_ES_MGS_4_each_k
  use m_ES_ortho            ,only : m_ES_orthogonalize_SD_to_WFs
  use m_ES_nonlocal         ,only : sc,ss,qc                                          &
 &                                , m_ES_alloc_afft_scss_etc                          &
 &                                , m_ES_dealloc_afft_scss_etc
  use m_ES_nonlocal         ,only : m_ES_Vnonlocal_W                                  &
 &                                , m_ES_betar_dot_Psi_4_each_k
  use m_FFT,                 only : nfft, fft_box_size_WF
  use m_FFT,                 only : m_FFT_Vlocal_W, m_FFT_WF
  use m_Ionic_System,        only : ntyp, natm, ityp, iwei
  use m_Kpoints,             only : kv3,vkxyz,k_symmetry
  use m_NonLocal_Potential,  only : snl
  use m_Files,               only : nfout
  use m_Parallelization,     only : MPI_CommGroup,is_kngp,ie_kngp,npes,mype &
       &, nrank_e, myrank_e, map_e, ista_e, iend_e, istep_e, map_z, np_e &
       &, nrank_k, myrank_k, map_k, ista_k, iend_k, mpi_k_world, nis_k, nel_k, ierr, ista_fs, iend_fs &
       &, nrank_s, ista_spin, iend_spin
  use m_PlaneWaveBasisSet,   only : kg1,iba, igf, nbase, m_pwBS_kinetic_energies
  use m_PseudoPotential,     only : ilmt, lmtt, ltp, mtp, q, dion, modnrm, nlmta &
       &                           ,m_PP_include_vanderbilt_pot,ipaw,dion_paw
  use m_Timing,              only : tstatc0_begin, tstatc0_end
  use m_ES_ExactExchange,    only : m_ES_Vexx_W, m_ES_EXX_eigenvalue_for_each_k &
       &                          , m_ES_EXX_Diagonal_part, m_ES_EXX_add_Diagonal_part &
                                  , m_ES_EXX_cp_eigenvalue, m_ES_EXX_gather_valence_states
#ifdef NEC_TIMER
  use nec_timer
#endif
  use m_FiniteElectricField, only : m_FEF_add_grad_to_vnlph, m_FEF_build_grad

! ==================================== added by K. Tagami ============= 11.0
  use m_Control_Parameters,    only : noncol, ndim_spinor, ndim_chgpot, ndim_magmom, &
       &                              SpinOrbit_mode, sw_hubbard
  use m_Const_Parameters,    only : CMPLDP, Neglected
  use m_PseudoPotential,       only : dion_scr_noncl, q_noncl
  use m_ES_ortho,              only : m_ES_orthogonl_SD_to_WFs_noncl, &
       &                              m_ES_MGS_4_each_k_noncl

  use m_Electronic_Structure,   only : m_ES_sum_of_LocalPart2_noncl, &
       &                               m_ES_Vlocal_in_Rspace_noncl, &
       &                               m_ES_sort_eigen_vals_noncl, &
       &                               m_ES_eigen_vals_each_k_noncl, &
       &                               m_ES_sum_of_LocalPart_noncl, &
       &                               m_ES_sum_of_LocalPart3_noncl

  use m_ES_nonlocal,            only : m_ES_alloc_scss_etc, &
       &                               m_ES_dealloc_scss_etc
  use m_FFT,                    only : m_FFT_Vlocal_W_noncl
!! for debug
  use m_Electronic_Structure,     only : fsr_l, fsi_l
! ===================================================================== 11.0

! === 0 divide occurs if abs(wdi) < epsilon. by tkato 2012/12/18 ===============
  use m_Const_Parameters, only : SmallestPositiveNumber
! ==============================================================================

  use m_Control_Parameters,  only : use_metagga, vtau_exists

  use m_IterationNumbers, only : iteration_electronic
  use mpi

  implicit none
  real(kind=DP),private,allocatable,dimension(:,:,:,:):: zaj_old !d(kg1,np_e,ista_k:iend_k,kimg)
  real(kind=DP),private,allocatable,dimension(:,:,:,:):: wfsd_l  !d(kg1,np_e,ik:ik,kimg)
! ----------- Added by T. Yamasaki, 28 June 2008 -----
  real(kind=DP),private,allocatable,dimension(:,:,:,:):: wfred_l !d(kg1,np_e,ista_k:iend_k,kimg)
  logical, private                                    :: wfred_is_allocated = .false.
! ----------------------------------------------------
  real(kind=DP),public ,allocatable,dimension(:,:,:,:):: wfsd_old!d(kg1,np_e,ista_k:iend_k,kimg)
  real(kind=DP),private,allocatable,dimension(:,:,:,:):: wfsd_np !d(kg1,np_e,1,kimg)
  real(kind=DP),private,allocatable,dimension(:,:,:)  :: bsdr_l, bsdi_l !d(np_e,nlmta,1)
  real(kind=DP),private,allocatable,dimension(:,:,:)  :: bsdr_old, bsdi_old !d(np_e,nlmta,1)
  real(kind=DP),private,allocatable,dimension(:,:,:)  :: bsdr_np, bsdi_np !d(np_e,nlmta,1)

  real(kind=DP),private,allocatable,dimension(:)          :: dzajn2  !d(kv3)
  real(kind=DP),private                               :: betacg

!   1. m_ESsd_alloc_dzajn2            <-(Initial_Electronic_Structure)
!   2. m_ESsd_dealloc_dzajn2          <-(Finalization_of_mpi)
!   3. m_ESsd_alloc_zaj_old           <-(Initial_Electronic_Structure)
!   4. m_ESsd_evolve_WFs_again        <-(Renewal_of_WaveFunctions) ->(5)
!   5. evolve_each_WF_again           <-(4)
!   6. m_ESsd_decide_CG_direction     <-(Renewal_of_WaveFunctions) ->(7)
!   7. decide_CG_direction_core       <-(6,) ->(16,29)
!   8. m_ESsd_renew_WF_by_SDorCG      <-(Renewal_of_WaveFunctions) ->(9,10,11,12,13,14)
!   9. evolve_WFs_in_MSD_direction    <-(8) ->(15,19)
!  10. evolve_WFs_in_SD_direction     <-(8)
!  11. evolve_WFs_in_eazyCG_direction <-(8) ->(16,17,18,26,28,30,33)
!  12. evolve_WFs_in_CG_direction     <-(8) ->(18,26,28,29)
!  13. evolve_WFs_in_CG_direction0    <-(8) ->(26,27)
!  14. vlhxc_l_zero_term              <-(8)
!  15. modified_steepest_descent      <-(9)
!  16. square_of_SD_direction(ik,ibo,dz) <-(7,11)
!  17. square_of_SD_direction2(ik,ibo,dz) <-(11)
!  18. SD_direction                   <-(11,12)
!  19. Vnonlocal_Diagonal_part        <-(9)
!  20.  - Vnonlocal_D_vanderbilt_case, - Vnonlocal_D_norm_conserve_case
!  21. m_ESsd_reset_dzajn2            <-(Renewal_of_WaveFunctions)
!  22. m_ESsd_copy_zaj_to_zaj_old     <-(Renewal_of_WaveFunctions)
!  23. m_ESsd_copy_zaj_old_to_zaj     <-(m_ES_WF_by_RMM.zajold2zaj_phi2zaj_old)
!  24. m_ESsd_copy_phi_to_zaj_old     <-(m_ES_WF_by_RMM.zajold2zaj_phi2zaj_old)
!  25. m_ESsd_diff_WFs_to_zaj_old     <-(m_ES_WF_by_RMM.m_ESrmm_renew_WF)
!  26. WF_conjugate_gradient          <-(11,12,13)
!  27. WF_conjugate_gradient0         <-(13)
!  28. make_CG_direction              <-(11,12)
!  29. orthogonalize_SD_drctns        <-(7,12)
!  30. Precon_dot_SD                  <-(11)  ->(31)
!  31. decide_precon_factor_wfsd      <-(30)  ->(32)
!  32. kinetic_energy_wfsd            <-(31)
!  33. cp_wfsd_to_wfsd_old            <-(11)

!  include 'mpif.h'

contains
  subroutine m_ESsd_alloc_dzajn2
    allocate(dzajn2(ista_k:iend_k)); dzajn2 = 0.d0
  end subroutine m_ESsd_alloc_dzajn2

  subroutine m_ESsd_dealloc_dzajn2
    if(allocated(zaj_old)) deallocate(zaj_old)
    if(allocated(dzajn2)) deallocate(dzajn2)
  end subroutine m_ESsd_dealloc_dzajn2

  subroutine m_ESsd_alloc_zaj_old
    allocate(zaj_old(kg1,np_e,ista_k:iend_k,kimg))
  end subroutine m_ESsd_alloc_zaj_old

  subroutine m_ESsd_alloc_wfsd
    if(.not.allocated(wfsd_old)) then
       allocate(wfsd_old(kg1,np_e,ista_k:iend_k,kimg))
       wfsd_old = 0.d0
    end if
    if(.not.allocated(bsdr_old)) then
       allocate(bsdr_old(np_e,nlmta,ista_k:iend_k));bsdr_old=0.0d0
    endif
    if(.not.allocated(bsdi_old)) then
       allocate(bsdi_old(np_e,nlmta,ista_k:iend_k));bsdi_old=0.0d0
    endif
  end subroutine m_ESsd_alloc_wfsd

  subroutine m_ESsd_dealloc_wfsd
    if(allocated(wfsd_old)) deallocate(wfsd_old)
    if(allocated(bsdr_old)) deallocate(bsdr_old)
    if(allocated(bsdi_old)) deallocate(bsdi_old)
  end subroutine m_ESsd_dealloc_wfsd

! -------------------- Added by T. Yamasaki, 28 June 2008 --->>
  subroutine m_ESsd_alloc_wfred()
    allocate(wfred_l(kg1,np_e,ista_k:iend_k,kimg)); wfred_l = 0.d0
    wfred_is_allocated = .true.
  end subroutine m_ESsd_alloc_wfred

  subroutine m_ESsd_dealloc_wfred()
    deallocate(wfred_l)
    wfred_is_allocated = .false.
  end subroutine m_ESsd_dealloc_wfred

  subroutine cp_wfred(ik)
    integer, intent(in) :: ik
    integer :: ir, ib, ig
    do ir = 1, kimg
       do ib = 1, np_e
          do ig = 1, iba(ik)
             wfred_l(ig,ib,ik,ir) = zaj_l(ig,ib,ik,ir) - zaj_old(ig,ib,ik,ir)
          end do
       end do
    end do
  end subroutine cp_wfred
! ----------------------------------<<

  subroutine alloc_wfsd_bsdri(ik)
    integer, intent(in) :: ik
    allocate(wfsd_l(kg1,np_e,ik:ik,kimg)); wfsd_l = 0.d0
    allocate(bsdr_l(np_e,nlmta,ik:ik))
    allocate(bsdi_l(np_e,nlmta,ik:ik)); bsdi_l = 0.0d0
!!$    write(nfout,'(" -- bsdr_l, bsdi_l have been allocated --  np_e, nlmta, ik = ",3i5)') np_e,nlmta,ik
  end subroutine alloc_wfsd_bsdri

! ==================================== added by K. Tagami ============== 11.0
  subroutine alloc_wfsd_bsdri_noncl( ik1, ik2 )
    integer, intent(in) :: ik1, ik2
    allocate(wfsd_l(kg1,np_e,ik1:ik2,kimg)); wfsd_l = 0.d0
    allocate(bsdr_l(np_e,nlmta,ik1:ik2))
    allocate(bsdi_l(np_e,nlmta,ik1:ik2)); bsdi_l = 0.0d0
  end subroutine alloc_wfsd_bsdri_noncl
! ====================================================================== 11.0

  subroutine dealloc_wfsd_bsdri()
    deallocate(wfsd_l)
    deallocate(bsdr_l)
    deallocate(bsdi_l)
  end subroutine dealloc_wfsd_bsdri

  subroutine alloc_wfsdnp_bsdrinp(ik)
    integer, intent(in) :: ik
    allocate(wfsd_np(kg1,np_e,ik:ik,kimg))
    allocate(bsdr_np(np_e,nlmta,ik:ik))
    allocate(bsdi_np(np_e,nlmta,ik:ik))
  end subroutine alloc_wfsdnp_bsdrinp

! ========================================= added by K. Tagami =========== 11.0
  subroutine alloc_wfsdnp_bsdrinp_noncl(ik1, ik2)
    integer, intent(in) :: ik1, ik2
    allocate(wfsd_np(kg1,np_e,ik1:ik2,kimg))
    allocate(bsdr_np(np_e,nlmta,ik1:ik2))
    allocate(bsdi_np(np_e,nlmta,ik1:ik2))
  end subroutine alloc_wfsdnp_bsdrinp_noncl
! ====================================================================== 11.0

  subroutine dealloc_wfsdnp_bsdrinp()
    deallocate(wfsd_np)
    deallocate(bsdr_np)
    deallocate(bsdi_np)
  end subroutine dealloc_wfsdnp_bsdrinp

  subroutine m_ESsd_dealloc_zaj_old
    if(allocated(zaj_old)) deallocate(zaj_old)
!!$    deallocate(wfsd_l)
!!$    deallocate(wfsd_np)
!!$    deallocate(bsdr_l)
!!$    deallocate(bsdi_l)
!!$    deallocate(bsdr_np)
!!$    deallocate(bsdi_np)
  end subroutine m_ESsd_dealloc_zaj_old

  subroutine m_ESsd_evolve_WFs_again(nfout,isolver2,mode,dtim_old,dtim_new,iupdate)
#ifdef __TIMER__
    use m_Const_Parameters,   only : VDB, NORMCONSERVATION
#endif
    integer,       intent(in) :: nfout,isolver2 &
         &                     , mode   ! mode = {ORTHONORMALIZATION | NORMALIZATION}
    real(kind=DP), intent(in) :: dtim_old, dtim_new
    integer, intent(in), optional :: iupdate
    integer       :: ispin, ik, iksnl
    real(kind=DP), pointer, dimension(:) :: ekin
    real(kind=DP), allocatable, dimension(:,:) :: afft_spin
#ifdef __TIMER__
    integer       :: ierr
#endif

    call m_ES_alloc_afft_scss_etc()
    ekin => sc(1:kg1)
    if(nrank_s==2) then
      allocate(afft_spin(nfft,2))
      do ispin=1, nspin
       call m_ES_Vlocal_in_Rspace(ispin,afft_spin(:,ispin))                 ! (ptfft1)
      enddo
    endif
!    do ispin = 1, nspin, af+1
    do ispin = ista_spin, iend_spin, af+1
       if(nrank_s==1) then
         call m_ES_Vlocal_in_Rspace(ispin,afft)                 ! (ptfft1)
       else
         afft(:) = afft_spin(:,ispin)
       endif

       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) == myrank_k ) then                     ! MPI
             call m_pwBS_kinetic_energies(ik,vkxyz,ekin) !-(PWBS) (diakin) ->ekin
             call evolve_each_WF_again(ik,isolver2,dtim_new,dtim_old)  !-(m_ES_WF_by_SDorCG)
             if(ipri>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')
!!$             if(modnrm == EXECUT) call m_ES_betar_dot_WFs_4_each_k(nfout,ik)       ! -> fsr_l, fsi_l
             call m_ES_MGS_4_each_k(nfout,ik,mode)
             if(sw_hybrid_functional==ON) then
                if(sw_retard_eigval_evaluation==OFF)then
                   call m_ES_EXX_eigenvalue_for_each_k(ispin,ik,iupdate=iupdate)
                else
                   call m_ES_EXX_cp_eigenvalue(ik)
                endif
             endif
             call m_ES_eigen_values_for_each_k(ispin,ik,ekin,afft)
          end if
       end do
    end do
    if(sw_hybrid_functional==ON.and.sw_retard_eigval_evaluation==ON)then
       call m_ES_EXX_gather_valence_states(nfout)
       do ispin=1,nspin,af+1
          do ik=ispin, kv3+ispin-nspin, nspin
             if(map_k(ik) /= myrank_k) cycle ! MPI
             call m_ES_EXX_eigenvalue_for_each_k(ispin,ik,update_eko=.false.,iupdate=iupdate)
          enddo
       enddo
    endif
#ifdef LMM_PREVIOUS
    call m_ES_sort_eigen_values    !-(m_Elec.)
! --> T. Yamasaki, 21th Jul. 2009
!!$    call m_ES_sort_eigen_values    !-(m_Elec.)
! <--
#endif
    call m_ES_dealloc_afft_scss_etc()
    if(nrank_s==2) deallocate(afft_spin)

  end subroutine m_ESsd_evolve_WFs_again

! ============================= added by K. Tagami ======================== 11.0
  subroutine m_ESsd_evolve_WFs_again_noncl(nfout,isolver2,mode,dtim_old,dtim_new,iupdate)
    integer,       intent(in) :: nfout,isolver2 &
         &                     , mode   ! mode = {ORTHONORMALIZATION | NORMALIZATION}
    real(kind=DP), intent(in) :: dtim_old, dtim_new
    integer, intent(in), optional :: iupdate

    integer       :: ispin, ik, iksnl
    integer :: is

    real(kind=DP), pointer, dimension(:) :: ekin

    real(kind=DP), allocatable :: afft_kt(:,:)
    real(kind=DP), allocatable :: bfft_kt(:,:)

!    return

    call m_ES_alloc_scss_etc()
    allocate(afft_kt(nfft,ndim_chgpot)); afft_kt = 0.0d0      ! pot
    allocate(bfft_kt(nfft,ndim_spinor)); bfft_kt = 0.0d0      ! wfn

    ekin => sc(1:kg1)

    call m_ES_Vlocal_in_Rspace_noncl( afft_kt )   ! vlhxc_ss -> afft_kt

    Do ik=1, kv3, ndim_spinor
       if (map_k(ik) /= myrank_k ) cycle
       iksnl = (ik-1)/ndim_spinor + 1

       call m_pwBS_kinetic_energies( ik, vkxyz, ekin )     !-(PWBS) (diakin) ->ekin

       Do is=1, ndim_spinor
         call evolve_each_WF_again( ik+is-1, isolver2, dtim_new, dtim_old )
                          !-(m_ES_WF_by_SDorCG)
       End do

       call m_ES_MGS_4_each_k_noncl( nfout, ik, mode )

       if (sw_hybrid_functional==ON) call m_ES_EXX_eigenvalue_for_each_k(ispin,ik,iupdate=iupdate)

       call m_ES_eigen_vals_each_k_noncl( ik, ekin, afft_kt, bfft_kt )

    End do

#ifdef LMM_PREVIOUS
    call m_ES_sort_eigen_vals_noncl           !-(m_Elec.)

#endif
    call m_ES_dealloc_scss_etc()
    deallocate( afft_kt, bfft_kt )

  end subroutine m_ESsd_evolve_WFs_again_noncl
! ========================================================================= 11.0

  subroutine evolve_each_WF_again(ik,isolver2,dt_new,dt_old)
    integer,       intent(in) :: ik,isolver2
    real(kind=DP), intent(in) :: dt_new,dt_old
    integer ::        ir, ib, ig
    real(kind=DP) ::  dtt
    integer   :: id_sname = -1
    call tstatc0_begin('evolve_each_WF_again ', id_sname,1)

!!$    if(isolver2 == CG .or. isolver2 == eazyCG) then
#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) then
       do ib = 1, np_e
          status_saved_phifftr(ib,ik) = OLD
       end do
    end if
#endif

  if(isolver2 == eazyCG) then
       do ir = 1, kimg
          do ib = 1, np_e
             do ig = 1, iba(ik)
                zaj_l(ig,ib,ik,ir) = zaj_old(ig,ib,ik,ir) + dt_new*wfsd_l(ig,ib,ik,ir)
             end do
          end do
       end do
    else
       dtt = dt_new/dt_old
!!$       if(ipri >= 3) write(nfout,'(" dtt = ",f8.4)') dtt
! ------------------- Revised by T. Yamasaki,  31 Oct 2008 --->>
       if(wfred_is_allocated) then
          do ir = 1, kimg
             do ib = 1, np_e
                do ig = 1, iba(ik)
! ------------------- Revised by T. Yamasaki,  03 July 2008 --->>
                   zaj_l(ig,ib,ik,ir) = zaj_old(ig,ib,ik,ir) + dtt*wfred_l(ig,ib,ik,ir)
! ------------------------------------<<
                end do
             end do
          end do
       else if(.not.wfred_is_allocated) then
          do ir = 1, kimg
             do ib = 1, np_e
                do ig = 1, iba(ik)
                   zaj_l(ig,ib,ik,ir) = (1-dtt)*zaj_old(ig,ib,ik,ir) + dtt*zaj_l(ig,ib,ik,ir)
                end do
             end do
          end do
       end if
! <----------------
    end if
!!$    do ir = 1, kimg
!!$       do ib = 1, np_e                         ! MPI
!!$          do i = 1, iba(ik)
!!$             dzz = zaj_l(i,ib,ik,ir) - zaj_old(i,ib,ik,ir)
!!$             zaj_l(i,ib,ik,ir) = zaj_old(i,ib,ik,ir) + dtt*dzz
!!$          end do
!!$       end do
!!$    end do
    call tstatc0_end(id_sname)
  end subroutine evolve_each_WF_again

  subroutine m_ESsd_decide_CG_direction(precon)
    integer, intent(in) :: precon

    integer ispin, iksnl, ik
    real(kind=DP), parameter             :: Delta = 1.d-10
    real(kind=DP), pointer, dimension(:) :: ekin, p
    real(kind=DP)                        :: gmgm, gmmgmm, sumdz2, sumdz
    real(kind=DP), allocatable :: cfft(:)
    real(kind=DP), allocatable, dimension(:,:) :: afft_spin

    if(sw_fef==ON) call m_FEF_build_grad() !->grad_fef
    call m_ES_alloc_afft_scss_etc()
    ekin => sc(1:kg1); p => ss(1:kg1)
    gmgm = 0.d0; gmmgmm = 0.d0

    if ( use_metagga .and. vtau_exists ) allocate( cfft(nfft) )
    if(nrank_s==2) then
      allocate(afft_spin(nfft,2))
      do ispin=1,nspin
        call m_ES_Vlocal_in_Rspace(ispin,afft_spin(:,ispin))      ! (ptfft1) -> afft
      enddo
    endif
    do ispin = 1, nspin, (af+1)
       if(nrank_s==1) then
         call m_ES_Vlocal_in_Rspace(ispin,afft)      ! (ptfft1) -> afft
       else
         afft(:) = afft_spin(:,ispin)
       endif
       if ( use_metagga .and. vtau_exists ) then
          call m_ES_Vlocal_in_Rspace( ispin, cfft, vtau_l )    ! r space
       endif

       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) == myrank_k) then           ! MPI
             iksnl = (ik-1)/nspin + 1
             call alloc_wfsd_bsdri(ik)
             call m_ES_Vnonlocal_W(ik,iksnl,ispin,ON)    ! (nonloc) ->(vnlph_l)
             if(sw_hybrid_functional==ON) call m_ES_Vexx_W(ik) ! (exx) ->(vnlph_l)
             if(sw_fef==ON) call m_FEF_add_grad_to_vnlph(ik) !->vnlph_l

             if ( use_metagga .and. vtau_exists ) then
                call m_ES_contrib_kindens_to_vnlph( ispin, ik, cfft )
             endif

             call m_pwBS_kinetic_energies(ik,vkxyz,ekin) ! -(m_PWBS) (diakin) ->ekin
             call decide_CG_direction_core &     ! -(m_ES_WF_by_SDorCG) ->wfsd_l
                  &(precon,ik,ekin,afft,bfft,p,sumdz2,sumdz)
             gmgm   = gmgm + sumdz
             gmmgmm = gmmgmm + dzajn2(ik)
             dzajn2(ik) = sumdz2
             if(ipri >= 2) write(nfout,'(" dzajn2(",i3,") = ",d20.8)') ik,dzajn2(ik)
             call dealloc_wfsd_bsdri()
          end if
       enddo
    enddo
    if ( allocated(cfft) ) deallocate( cfft )

    if(npes > 1) then
       call mpi_allreduce(MPI_IN_PLACE,gmgm,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)    ! MPI
       call mpi_allreduce(MPI_IN_PLACE,gmmgmm,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)    ! MPI
    end if

    if(gmmgmm > Delta) then
       betacg = gmgm/gmmgmm
    else
       betacg = 0.d0
    endif
    if(ipri >= 1) write(nfout,'(" ! betacg = ",d20.8," gmgm = ",d20.8," gmmgmm = ",d20.8)') betacg,gmgm,gmmgmm
    if(betacg > 1.d0) then
       betacg = 0.d0
       if(ipri >= 1) write(nfout,*) ' beta for CG is too large, and it sets to be zero'
    else if(betacg < 0.d0) then
       betacg = 0.d0
       if(ipri >= 1) write(nfout,*) ' beta for CG is less than zero, and it sets to be zero'
    end if

    call m_ES_dealloc_afft_scss_etc()
!!$  contains
!!$    subroutine bcast_dzajn2                     ! MPI
!!$      integer       :: i
!!$
!!$      do i = 0, nrank_k-1
!!$         call mpi_bcast(dzajn2(nis_k(i)),nel_k(i),mpi_double_precision &
!!$              & ,i*nrank_e,MPI_CommGroup,ierr)
!!$      end do
!!$    end subroutine bcast_dzajn2
  end subroutine m_ESsd_decide_CG_direction

! ====================================== added by K. Tagami ============== 11.0
  subroutine m_ESsd_decide_CG_direc_noncl(precon)
    integer, intent(in) :: precon

    integer ispin, iksnl, ik
    integer :: is1, is2, k2, istmp

    real(kind=DP), parameter             :: Delta = 1.d-10
    real(kind=DP), pointer, dimension(:) :: ekin, p
    real(kind=DP)                        :: gmgm, gmmgmm, sumdz2, sumdz

    real(kind=DP), allocatable :: vnlph_noncl(:,:,:,:)
    real(kind=DP), allocatable :: afft_kt(:,:)
    real(kind=DP), allocatable :: bfft_kt(:,:)

    call m_ES_alloc_scss_etc()
    ekin => sc(1:kg1); p => ss(1:kg1)
    allocate(afft_kt(nfft,ndim_chgpot)); afft_kt = 0.0d0      ! pot
    allocate(bfft_kt(nfft,ndim_spinor)); bfft_kt = 0.0d0      ! wfn

    allocate(vnlph_noncl(kg1,np_e,kimg,ndim_spinor)); vnlph_noncl = 0.0d0

    gmgm = 0.d0; gmmgmm = 0.d0

    call m_ES_Vlocal_in_Rspace_noncl( afft_kt )

    do ik = 1, kv3, ndim_spinor

      if (map_k(ik) /= myrank_k) cycle
      iksnl = (ik-1)/ndim_spinor + 1

      call alloc_wfsd_bsdri_noncl( ik, ik+ndim_spinor-1 )

      vnlph_noncl = 0.0d0

      Do is1=1, ndim_spinor
          Do is2=1, ndim_spinor
            istmp = ( is1-1 )*ndim_spinor +is2
            k2 = ik +is2 -1

            call m_ES_Vnonlocal_W( k2, iksnl, istmp, ON )  ! (nonloc) ->(vnlph_l)

            vnlph_noncl(:,:,:,is1) = vnlph_noncl(:,:,:,is1) &
             &                     + vnlph_l(:,:,:)
          End do
       End do

      if (sw_hybrid_functional==ON) call m_ES_Vexx_W(ik) ! (exx) ->(vnlph_l)

      call m_pwBS_kinetic_energies( ik, vkxyz, ekin ) ! -(m_PWBS) (diakin) ->ekin

      call decide_CG_direction_core_noncl( precon, ik, ekin, &
              &                               afft_kt, bfft_kt, vnlph_noncl, &
              &                               p, sumdz2, sumdz )

!!!!!!!!!!!!!      gmgm   = gmgm + sumdz2
      gmgm   = gmgm + sumdz
      gmmgmm = gmmgmm + dzajn2(ik)
      dzajn2(ik) = sumdz2

      if(ipri >= 2) write(nfout,'(" dzajn2(",i3,") = ",d20.8)') ik,dzajn2(ik)
      call dealloc_wfsd_bsdri()

    enddo

    if(npes > 1) then
       call mpi_allreduce(MPI_IN_PLACE, gmgm,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)    ! MPI
       call mpi_allreduce(MPI_IN_PLACE, gmmgmm,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)    ! MPI
    end if

    if(gmmgmm > Delta) then
       betacg = gmgm/gmmgmm
    else
       betacg = 0.d0
    endif
    if(ipri >= 1) write(nfout,'(" ! betacg = ",d20.8," gmgm = ",d20.8," gmmgmm = ",d20.8)') betacg,gmgm,gmmgmm
    if(betacg > 1.d0) then
       betacg = 0.d0
       if(ipri >= 1) write(nfout,*) ' beta for CG is too large, and it sets to be zero'
    else if(betacg < 0.d0) then
       betacg = 0.d0
       if(ipri >= 1) write(nfout,*) ' beta for CG is less than zero, and it sets to be zero'
    end if

    call m_ES_dealloc_scss_etc()
    deallocate( afft_kt, bfft_kt )
    deallocate( vnlph_noncl )

  end subroutine m_ESsd_decide_CG_direc_noncl
! ================================================================ 11.0

  subroutine decide_CG_direction_core(precon,ik,ekin,afft,bfft,p,sumdz2,sumdz)
    integer, intent(in)        :: precon, ik
    real(kind=DP), intent(in)  :: ekin(kg1)
    real(kind=DP), intent(in)  :: afft(nfft)
    real(kind=DP), intent(out) :: bfft(nfft)
    real(kind=DP)              :: p(kg1)
    real(kind=DP), intent(out) :: sumdz2,sumdz

    integer :: ib
    integer :: id_sname = -1
    call tstatc0_begin('decide_CG_direction_core ', id_sname,1)

    do ib = ista_e, iend_e, istep_e      ! MPI
       call m_ES_WF_in_Rspace(ik,ib,bfft)! (swffft)
       call m_FFT_Vlocal_W(afft,bfft)    ! (afft, bfft) -> (bfft)
       call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON) ! bfft: G-space repres.
       call SD_direction(precon,ik,ib,ekin,bfft,p) !-(m_ES_WF_by_SDorCG)
    end do

    call orthogonalize_SD_drctns(ik,to=OTHER_BANDS)  ! -(m_ES_WF_by_SDorCG) ->(wfsd_l, bsd(ri)_l)

    call WFSD_dot_WFSD                   ! -(contained here)
    !    ~~~~~~~~~~~~~
    call tstatc0_end(id_sname)
  contains
    subroutine WFSD_dot_WFSD
      real(kind=DP) :: dz,dz2
      integer :: ib

      sumdz2 = 0.d0
      sumdz = 0.0d0
      do ib = ista_e, iend_e, istep_e      ! MPI
         call square_of_SD_direction(ik,ib,dz)
         dz2 = 0.d0
         if(sw_modified_cg_formula==ON) call square_of_SD_direction3(ik,ib,dz2)
         sumdz2 = sumdz2 + dz
         sumdz = sumdz + (dz-dz2)
!!$         sumdz = sumdz + dz
      end do
      if(npes > 1) then
         call mpi_allreduce(MPI_IN_PLACE, sumdz2,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
         call mpi_allreduce(MPI_IN_PLACE, sumdz ,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
      end if
    end subroutine WFSD_dot_WFSD
  end subroutine decide_CG_direction_core

! ================================== added by K. Tagami ================== 11.0
  subroutine decide_CG_direction_core_noncl( precon, ik, ekin, afft_kt, bfft_kt, &
       &                                    vnlph_noncl, p, sumdz2, sumdz )
    integer, intent(in)        :: precon, ik
    real(kind=DP), intent(in)  :: ekin(kg1)

    real(kind=DP), intent(inout)  :: afft_kt(nfft,ndim_chgpot)
    real(kind=DP), intent(out) :: bfft_kt(nfft,ndim_spinor)
    real(kind=DP), intent(in)  :: vnlph_noncl(kg1,np_e,kimg,ndim_spinor)

    real(kind=DP)              :: p(kg1)
    real(kind=DP), intent(out) :: sumdz2, sumdz

    integer :: ib, is
    integer :: id_sname = -1
    call tstatc0_begin('decide_CG_direction_core_noncl ', id_sname,1)

    do ib = ista_e, iend_e, istep_e      ! MPI

       Do is=1, ndim_spinor
         call m_ES_WF_in_Rspace( ik+is-1, ib, bfft_kt(:,is) )! (swffft)
       End do

       call m_FFT_Vlocal_W_noncl( afft_kt, bfft_kt, ndim_chgpot, ndim_spinor )
                                           ! (afft, bfft)-> (bfft)
       Do is=1, ndim_spinor
         call m_FFT_WF( ELECTRON, nfout, bfft_kt(:,is), DIRECT, ON )
       End do

       call SD_direction_noncl( precon, ik, ib, ekin, bfft_kt, vnlph_noncl, p )
    end do

    call orthogonalize_SD_drctns_noncl( ik,to=OTHER_BANDS )
                               ! -(m_ES_WF_by_SDorCG) ->(wfsd_l, bsd(ri)_l)

    call WFSD_dot_WFSD_noncl                   ! -(contained here)

    call tstatc0_end(id_sname)

  contains

    subroutine WFSD_dot_WFSD_noncl
      real(kind=DP) :: dz, dz2
      integer :: ib

      sumdz2 = 0.d0
      sumdz = 0.0
      do ib = ista_e, iend_e, istep_e
         call square_of_SD_direction_noncl( ik, ib, dz )
         dz2 = 0.0d0

         if(sw_modified_cg_formula==ON) call square_of_SD_direction3_noncl(ik,ib,dz2)

         sumdz2 = sumdz2 + dz
         sumdz = sumdz + ( dz - dz2 )
      end do
      if(npes > 1) then
         call mpi_allreduce(MPI_IN_PLACE,sumdz2,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
         call mpi_allreduce(MPI_IN_PLACE,sumdz, 1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
      end if
    end subroutine WFSD_dot_WFSD_noncl

  end subroutine decide_CG_direction_core_noncl
! ====================================================================== 11.0

  subroutine m_ESsd_renew_WF_by_SDorCG(nfout,isolver,precon,dtim,iupdate)
    use m_IterationNumbers,     only : iteration
    use m_Control_Parameters,   only : sw_gep
#ifdef __TIMER__
    use m_Const_Parameters,   only : VDB, NORMCONSERVATION
#endif

    integer, intent(in) :: nfout,isolver, precon
    integer, intent(in), optional :: iupdate
    real(kind=DP)       :: dtim

    integer                              :: ispin, iksnl, ik, mode, ipri0
    real(kind=DP), pointer, dimension(:) :: ekin, p, vnldi
    real(kind=DP)                        :: vlhxc0
    real(kind=DP), allocatable :: cfft(:)
    real(kind=DP), allocatable, dimension(:,:) :: afft_spin

    if(sw_gep == ON.and.iteration>=2) then
      mode = NORMALIZATION
    else
      mode = ORTHONORMALIZATION
    endif
    call m_ES_alloc_afft_scss_etc()
    ekin => sc(1:kg1); p => ss(1:kg1); vnldi => qc(1:kg1)

    if(iprisolver >= 2) then
       if(isolver == SD) then
          write(nfout,'(" !! isolver_core = SD")')
       else if(isolver == MSD) then
          write(nfout,'(" !! isolver_core = MSD")')
       else if(isolver == CG) then
          write(nfout,'(" !! isolver_core = CG")')
       end if
    end if

    if ( use_metagga .and. vtau_exists ) allocate( cfft(nfft) )

    if(sw_fef==ON) call m_FEF_build_grad() !->grad_fef
    if(nrank_s==2) then
      allocate(afft_spin(nfft, 2))
      do ispin = 1, nspin
        call m_ES_Vlocal_in_Rspace(ispin,afft_spin(:,ispin))      ! (ptfft1) vlhxc_l->afft
      enddo
    endif
    do ispin = 1, nspin, (af+1)
       if(isolver == MSD) call vlhxc_l_zero_term(vlhxc0,ispin) ! -(m_ES_WF_by_SDorCG) ->vlhxc0
       if(nrank_s==1) then
         call m_ES_Vlocal_in_Rspace(ispin,afft)      ! (ptfft1) vlhxc_l->afft
       else
         afft(:) = afft_spin(:,ispin)
       endif
       if ( use_metagga .and. vtau_exists ) then
          call m_ES_Vlocal_in_Rspace( ispin, cfft, vtau_l )    ! r space
       endif

       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) == myrank_k) then           ! MPI
             iksnl = (ik-1)/nspin + 1
             call m_ES_Vnonlocal_W(ik,iksnl,ispin,ON)    ! (nonloc) ->(vnlph_l)
             if(sw_hybrid_functional==ON) call m_ES_Vexx_W(ik) ! (exx) ->(vnlph_l)
             if(sw_fef==ON) call m_FEF_add_grad_to_vnlph(ik) !->vnlph_l

             if ( use_metagga .and. vtau_exists ) then
                call m_ES_contrib_kindens_to_vnlph( ispin, ik, cfft )
             endif

             call m_pwBS_kinetic_energies(ik,vkxyz,ekin) ! (diakin) ->ekin
             if(isolver == SD) then
                call evolve_WFs_in_SD_direction&      ! -(m_ES_WF_by_SDorCG)
                     &(precon,ik,dtim,ekin, afft,bfft,p)
             else if(isolver == MSD ) then
                call evolve_WFs_in_MSD_direction&     !-(m_ES_WF_by_SDorCG)
                     &(precon,ik,iksnl,ispin,dtim,ekin, afft,bfft,p,vnldi,vlhxc0)
                if(ipri>=3.and.ik==1) &
                     & call m_ES_wd_zaj_small_portion(nfout,ik," -- after MSD --",16)
             else if(isolver == eazyCG) then
                call evolve_WFs_in_eazyCG_direction(precon,ik,dtim,ekin,afft,bfft,p)
             else if(isolver == CG) then
                call evolve_WFs_in_CG_direction(precon,ik,dtim,ekin,afft,bfft,p)
!!$                call evolve_WFs_in_CG_direction0(precon,ik,dtim,ekin,afft,bfft,p)
             endif
!!$             if(modnrm == EXECUT) call m_ES_betar_dot_WFs_4_each_k(nfout,ik) !-> fsr_l,fsi_l
! --------------- Added by T. Yamasaki, 28 June 2008 ---
             if(wfred_is_allocated) call cp_wfred(ik)
! ------------------------------------------------------<<
             call m_ES_MGS_4_each_k(nfout,ik,mode)
             if(sw_hybrid_functional==ON) then
                if(sw_retard_eigval_evaluation==OFF)then
                   call m_ES_EXX_eigenvalue_for_each_k(ispin,ik,iupdate=iupdate)
                else
                   call m_ES_EXX_cp_eigenvalue(ik)
                endif
             endif
             call m_ES_eigen_values_for_each_k(ispin,ik,ekin,afft)
          end if
       enddo
    enddo
    if ( allocated(cfft) ) deallocate(cfft)

    if(sw_hybrid_functional==ON.and.sw_retard_eigval_evaluation==ON)then
       call m_ES_EXX_gather_valence_states(nfout)
       do ispin=1,nspin,af+1
          do ik=ispin, kv3+ispin-nspin, nspin
             if(map_k(ik) /= myrank_k) cycle ! MPI
             call m_ES_EXX_eigenvalue_for_each_k(ispin,ik,update_eko=.false.,iupdate=iupdate)
          enddo
       enddo
    endif

!!$    if(nspin == 2 .and. sw_so == ON) then
!!$       do ik = 1, kv3, nspin

#ifdef LMM_PREVIOUS
    call m_ES_sort_eigen_values()   ! -> neordr, nrvf_ordr
! --> T. Yamasaki, 21th Jul. 2009
!!$    call m_ES_sort_eigen_values()   ! -> neordr, nrvf_ordr
! <--
#endif
    call get_ipri0(ipri, ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko(nfout,mode=SCF)
    call m_ES_dealloc_afft_scss_etc()
  contains
    subroutine get_ipri0(ipri_in, ipri_out)
      integer, intent(in)  :: ipri_in
      integer, intent(out) :: ipri_out
      if(npes > 1) then
         if(mype == 0) ipri_out = ipri_in
         call mpi_bcast(ipri_out,1,mpi_integer,0,MPI_CommGroup,ierr)
      else
         ipri_out = ipri_in
      end if
    end subroutine get_ipri0
  end subroutine m_ESsd_renew_WF_by_SDorCG

! =============================== added by K. Tagami ================== 11.0
  subroutine m_ESsd_renew_WF_by_SDorCG_noncl(nfout,isolver,precon,dtim,iupdate)
    integer, intent(in) :: nfout,isolver, precon
    integer, intent(in), optional :: iupdate
    real(kind=DP)       :: dtim

    integer                              :: ispin, iksnl, ik, mode, ipri0
    integer :: is1, is2, istmp, k2

    real(kind=DP), pointer, dimension(:) :: ekin, p, vnldi

    real(kind=DP)                        :: vlhxc0( ndim_chgpot )
    real(kind=DP), allocatable :: afft_kt(:,:)
    real(kind=DP), allocatable :: bfft_kt(:,:)
    real(kind=DP), allocatable :: vnlph_noncl(:,:,:,:)

    call m_ES_alloc_scss_etc()

! --
!   Note : The followig is done in order to prevent allocating new arrays.
!
    ekin => sc(1:kg1); p => ss(1:kg1); vnldi => qc(1:kg1)
! --
    allocate(afft_kt(nfft,ndim_chgpot)); afft_kt = 0.0d0      ! pot
    allocate(bfft_kt(nfft,ndim_spinor)); bfft_kt = 0.0d0      ! wfn

    allocate(vnlph_noncl(kg1,np_e,kimg,ndim_spinor)); vnlph_noncl = 0.0d0

    mode = ORTHONORMALIZATION

    if(iprisolver >= 2) then
       if(isolver == SD) then
          write(nfout,'(" !! isolver_core = SD")')
       else if(isolver == MSD) then
          write(nfout,'(" !! isolver_core = MSD")')
       else if(isolver == CG) then
          write(nfout,'(" !! isolver_core = CG")')
       end if
    end if

    call m_ES_Vlocal_in_Rspace_noncl( afft_kt )      ! (ptfft1) vlhxc_l->afft
    if ( isolver == MSD ) then
       call vlhxc_l_zero_term_noncl( vlhxc0 )
    endif

    Do ik=1, kv3, ndim_spinor
       if (map_k(ik) /= myrank_k) cycle
       iksnl = (ik-1)/ndim_spinor + 1

       vnlph_noncl = 0.0d0

       Do is1=1, ndim_spinor
          Do is2=1, ndim_spinor
             istmp = ( is1-1 )*ndim_spinor +is2
             k2 = ik +is2 -1

             call m_ES_Vnonlocal_W( k2, iksnl, istmp, ON )
             ! (nonloc) ->(vnlph_l)

             vnlph_noncl(:,:,:,is1) = vnlph_noncl(:,:,:,is1) &
                  &                     + vnlph_l(:,:,:)
          End do
       End do

       if (sw_hybrid_functional==ON) call m_ES_Vexx_W(ik) ! (exx) ->(vnlph_l)

       call m_pwBS_kinetic_energies(ik,vkxyz,ekin) ! (diakin) ->ekin

       select case (isolver)
       case (SD)
          call evolve_WFs_in_SD_direc_noncl( precon, ik, dtim, ekin, &
               &                             afft_kt, bfft_kt, vnlph_noncl, p )
       case (MSD)
          call evolve_WFs_in_MSD_direc_noncl( precon, ik, iksnl, &
               &                              dtim, ekin, afft_kt, bfft_kt,&
               &                              vnlph_noncl, p, vnldi, vlhxc0 )
          if ( ipri>=3.and.ik==1 ) then
             Do is1=1, ndim_spinor
                write(nfout,*) 'is = ', is1
                call m_ES_wd_zaj_small_portion( nfout, ik+is1-1, " -- after MSD  --",21)
             End do
          endif
       case (eazyCG)
          call evolve_WFs_in_ezCG_direc_noncl( precon, ik, dtim, ekin, &
               &                               afft_kt, bfft_kt, vnlph_noncl, p )

       case (CG)
          call evolve_WFs_in_CG_direc_noncl( precon, ik, dtim, ekin, &
               &                             afft_kt, bfft_kt, vnlph_noncl, p )
       end select


       if (wfred_is_allocated) then
          Do is1=1, ndim_spinor
             call cp_wfred(ik+is1-1)
          End do
       endif

       call m_ES_MGS_4_each_k_noncl( nfout, ik, mode )

       if (sw_hybrid_functional==ON) call m_ES_EXX_eigenvalue_for_each_k(ispin,ik,iupdate=iupdate)

       call m_ES_eigen_vals_each_k_noncl( ik, ekin, afft_kt, bfft_kt )

    End do


#ifdef LMM_PREVIOUS
    call m_ES_sort_eigen_vals_noncl
#endif
    call get_ipri0(ipri, ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko(nfout,mode=SCF)
    call m_ES_dealloc_scss_etc()
    deallocate( afft_kt, bfft_kt )

    deallocate( vnlph_noncl )

  contains

    subroutine get_ipri0(ipri_in, ipri_out)
      integer, intent(in)  :: ipri_in
      integer, intent(out) :: ipri_out
      if(npes > 1) then
         if(mype == 0) ipri_out = ipri_in
         call mpi_bcast(ipri_out,1,mpi_integer,0,MPI_CommGroup,ierr)
      else
         ipri_out = ipri_in
      end if
    end subroutine get_ipri0
  end subroutine m_ESsd_renew_WF_by_SDorCG_noncl
! ======================================================================== 11.0

  subroutine evolve_WFs_in_MSD_direction&
       &(precon,ik,iksnl,ispin,dtim,ekin,afft,bfft,p,vnldi,vlhxc0)
    integer, intent(in)        :: precon, ik, iksnl, ispin
    real(kind=DP), intent(in)  :: dtim,ekin(kg1)
    real(kind=DP), intent(in)  :: afft(nfft)
    real(kind=DP), intent(out) :: bfft(nfft)
    real(kind=DP)              :: p(kg1), vnldi(kg1)
    real(kind=DP), intent(in)  :: vlhxc0

    real(kind=DP), allocatable :: vxdi(:) !d(kg1)

    integer ib, ig
    integer :: id_sname = -1
    call tstatc0_begin('evolve_WFs_in_MSD_direction ', id_sname,1)

    if(sw_hybrid_functional==ON) then
       allocate(vxdi(kg1))
       call m_ES_EXX_Diagonal_part(ispin,ik,vxdi)
    end if

#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(bfft,vnldi,p)
#endif
    do ib = ista_e, iend_e, istep_e      ! MPI
       call m_ES_WF_in_Rspace(ik,ib,bfft)! (swffft)

       if(ipri >= 2) then
          write(nfout,'(" --- afft and bfft <<evolve_WFs_in_MSD_direction>>,  ib = ",i5)') ib
          write(nfout,'(" afft: ",5d16.8)') (afft(ig),ig=1,min(5,nfft))
          write(nfout,'(" bfft: ",5d16.8)') (bfft(ig),ig=1,min(5,nfft))
       end if

       call m_FFT_Vlocal_W(afft,bfft)    ! (afft, bfft)-> (bfft)
       call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
       call Vnonlocal_Diagonal_part(ispin,ik,iksnl,ib,vnldi)
                                         ! -(m_ES_WF_by_SDorCG) (nonldj)
       if(sw_hybrid_functional==ON) call m_ES_EXX_add_Diagonal_part(ik,ib,vxdi,vnldi)
       if(ipri >= 2) then
          write(nfout,'(" --- vnlph_l <<evolve_WFs_in_MSD_direction>>,  ib = ",i5)') ib
                       write(nfout,'(" vnlph_l: Re ",5d16.8)') (vnlph_l(ig,ib,1),ig=1,min(5,iba(ik)))
          if(kimg== 2) write(nfout,'(" vnlph_l: Im ",5d16.8)') (vnlph_l(ig,ib,2),ig=1,min(5,iba(ik)))

          write(nfout,'(" vnldi ",5d16.8)') (vnldi(ig),ig=1,min(5,iba(ik)))
          write(nfout,'(" VlocalW: Re ",5d16.8)') (bfft(2*igf(nbase(ig,ik))-1),ig=1,min(5,iba(ik)))
          write(nfout,'(" VlocalW: Im ",5d16.8)') (bfft(2*igf(nbase(ig,ik))),ig=1,min(5,iba(ik)))
          write(nfout,'(" vlhxc0 = ",5d16.8)') vlhxc0
       end if
       call modified_steepest_descent&   ! -(m_ES_WF_by_SDorCG)
            &(precon,ik,ib,dtim,vnldi,vlhxc0,ekin,bfft,p)
    end do
    if(ipri>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')

    if(sw_hybrid_functional==ON) deallocate(vxdi)

    call tstatc0_end(id_sname)
  end subroutine evolve_WFs_in_MSD_direction

! ========================================= added by K. Tagami =============== 11.0
 subroutine evolve_WFs_in_MSD_direc_noncl( precon, ik, iksnl, &
       &                                    dtim, ekin, afft_kt, bfft_kt, &
       &                                    vnlph_noncl, p, vnldi, vlhxc0 )

    integer, intent(in)        :: precon, ik, iksnl
    real(kind=DP), intent(in)  :: dtim,ekin(kg1)
    real(kind=DP), intent(inout)  :: afft_kt(nfft,ndim_chgpot)
    real(kind=DP), intent(out) :: bfft_kt(nfft,ndim_spinor)
    real(kind=DP), intent(in)  :: vnlph_noncl(kg1,np_e,kimg,ndim_spinor)

    real(kind=DP)              :: p(kg1), vnldi(kg1)
    real(kind=DP), intent(in)  :: vlhxc0(ndim_chgpot)
!
#ifdef USE_COMPLEX_VNLDI
    complex(kind=CMPLDP), allocatable :: vnldi_noncl(:,:)
#else
    real(kind=DP), allocatable :: vnldi_noncl(:,:)
#endif
!
    integer ib, ig, is, istmp
    integer :: is1, is2

    integer :: id_sname = -1

    call tstatc0_begin('evolve_WFs_in_MSD_direc_noncl ', id_sname,1)

! ---------------------
    allocate( vnldi_noncl( kg1, ndim_spinor)); vnldi_noncl = 0.0d0
!    Do is=1, ndim_spinor
!      vnldi_noncl(:,is) = vnldi(:)
!    End do
! -----------------------

#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(bfft,vnldi,p)
#endif
    do ib = ista_e, iend_e, istep_e      ! MPI

       Do is=1, ndim_spinor
          call m_ES_WF_in_Rspace( ik +is -1, ib, bfft_kt(:,is) ) ! (swffft)
       End do

       call m_FFT_Vlocal_W_noncl( afft_kt, bfft_kt, ndim_chgpot, ndim_spinor )
                                                       ! (afft, bfft)-> (bfft)
       Do is=1, ndim_spinor
         call m_FFT_WF( ELECTRON, nfout, bfft_kt(:,is), DIRECT, ON )
       End do

       vnldi_noncl = 0.0d0
       Do is1=1, ndim_spinor
          Do is2=1, ndim_spinor
             istmp = ( is1 -1 )*ndim_spinor +is2

             call Vnonlocal_Diagonal_part_noncl( istmp, ik, iksnl, ib, &
                  &                              vnldi_noncl(:,is1) )
!             call Vnonlocal_Diagonal_part_noncl( istmp, ik, iksnl, ib, &
!                  &                              vnldi_noncl(:,is2) )
                                                    ! -(m_ES_WF_by_SDorCG) (nonldj)
          End do
       End do

       call modified_steepest_desc_noncl( precon, ik, ib, dtim, &
            &                             vnldi_noncl, vlhxc0, ekin,&
            &                             bfft_kt, vnlph_noncl, p )
    end do

    call tstatc0_end(id_sname)
    deallocate( vnldi_noncl)

  end subroutine evolve_WFs_in_MSD_direc_noncl
! ========================================================================= 11.0

  subroutine evolve_WFs_in_SD_direction&
       &(precon,ik,dtim,ekin,afft,bfft,p)
    integer, intent(in)        :: precon, ik
    real(kind=DP), intent(in)  :: dtim,ekin(kg1)
    real(kind=DP), intent(in)  :: afft(nfft)
    real(kind=DP), intent(out) :: bfft(nfft)
    real(kind=DP)              :: p(kg1)

    integer :: ib
    integer :: id_sname = -1
    call tstatc0_begin('evolve_WFs_in_SD_direction ', id_sname,1)

    do ib = ista_e, iend_e, istep_e      ! MPI
       call m_ES_WF_in_Rspace(ik,ib,bfft)! (swffft)
       call m_FFT_Vlocal_W(afft,bfft)    ! (afft, bfft)-> (bfft)
       call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
       call steepest_descent(precon,ik,ib,dtim,ekin,bfft,p)
                                         ! -(m_ES_WF_by_SDorCG)
    end do
    if(ipri>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')
    call tstatc0_end(id_sname)
  end subroutine evolve_WFs_in_SD_direction

! =============================== added by K. Tagami =================== 11.0
  subroutine evolve_WFs_in_SD_direc_noncl( precon, ik, dtim, ekin, &
       &                                   afft_kt, bfft_kt, vnlph_noncl, p )
    integer, intent(in)        :: precon, ik
    real(kind=DP), intent(in)  :: dtim,ekin(kg1)
    real(kind=DP), intent(inout)  :: afft_kt( nfft,ndim_chgpot )
    real(kind=DP), intent(out) :: bfft_kt(nfft,ndim_spinor)
    real(kind=DP), intent(in)  :: vnlph_noncl(kg1,np_e,kimg,ndim_spinor)

    real(kind=DP)              :: p(kg1)

    integer :: ib, is
    integer :: id_sname = -1
    call tstatc0_begin('evolve_WFs_in_SD_direc_noncl ', id_sname,1)

    do ib = ista_e, iend_e, istep_e      ! MPI
       Do is=1, ndim_spinor
         call m_ES_WF_in_Rspace( ik+is-1, ib, bfft_kt(:,is) )  ! (swffft)
       End do

       call m_FFT_Vlocal_W_noncl( afft_kt, bfft_kt, nspin, ndim_spinor )

       Do is=1, ndim_spinor
         call m_FFT_WF( ELECTRON, nfout, bfft_kt(:,is), DIRECT, ON )
       End do

       call steepest_descent_noncl( precon, ik, ib, dtim, ekin, bfft_kt, &
            &                       vnlph_noncl, p )
    end do

    call tstatc0_end(id_sname)

  end subroutine evolve_WFs_in_SD_direc_noncl
! =========================================================================== 11.0

  subroutine evolve_WFs_in_eazyCG_direction(precon,ik,dtim,ekin,afft,bfft,p)
    integer, intent(in)        :: precon, ik
    real(kind=DP), intent(in)  :: dtim,ekin(kg1)
    real(kind=DP), intent(in)  :: afft(nfft)
    real(kind=DP), intent(out) :: bfft(nfft)
    real(kind=DP)              :: p(kg1)

    real(kind=DP), parameter   :: Delta = 1.d-10
!!$    real(kind=DP)   :: dz, sumdz2
    integer :: ib
    integer :: id_sname = -1
    call tstatc0_begin('evolve_WFs_in_eazyCG_direction ', id_sname,1)

    if(modnrm == EXECUT) call alloc_wfsd_bsdri(ik)
    if(modnrm == EXECUT .and. precon == ON) call alloc_wfsdnp_bsdrinp(ik)

    do ib = ista_e, iend_e, istep_e      ! MPI
       call m_ES_WF_in_Rspace(ik,ib,bfft)! (swffft)
       call m_FFT_Vlocal_W(afft,bfft)    ! (afft, bfft) -> (bfft)
       call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON) ! bfft: G-space repres.
       call SD_direction(precon,ik,ib,ekin,bfft,p)! -(m_ES_WF_by_SDorCG) zaj_l, etc -> wfsd_l
       !    ~~~~~~~~~~~~
    end do

!!$    call orthogonalize_SD_drctns(ik,to=OTHER_BANDS)   ! -(m_ES_WF_by_SDorCG) ->wfsd_l
    if(precon == ON) call Precon_dot_SD(ik,ekin,p)  ! -(m_ES_WF_by_SDorCG) -> wfsd_np, wfsd_l
!!$    call orthogonalize_SD_drctns(ik,to=ALL_BANDS)     ! -(m_ES_WF_by_SDorCG) ->wfsd_l
    !                                      wfsd_l is orthogonalized to zaj_l's
    if(modnrm == EXECUT) call m_ES_betar_dot_Psi_4_each_k(nfout,wfsd_l,ik,ik,ik,bsdr_l,bsdi_l)
    if(modnrm == EXECUT .and. precon == ON ) then
       write(6,'(" m_ES_betar_dot_Psi_4_each_k is being called (evolve_WFs_in_eazyCG_direction)")')
       call m_ES_betar_dot_Psi_4_each_k(nfout,wfsd_np,ik,ik,ik,bsdr_np,bsdi_np)
    end if
    call get_betacg                      ! -(contained here) wfsd_l ->betacg, dzajn2
    !    ~~~~~~~~~~
!!$    do ib = ista_e, iend_e, istep_e      ! MPI
!!$       call m_ES_WF_in_Rspace(ik,ib,bfft)! (swffft)
!!$       call m_FFT_Vlocal_W(afft,bfft)    ! (afft, bfft) -> (bfft)
!!$       call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON) ! bfft: G-space repres.
!!$       call SD_direction(precon,ik,ib,ekin,bfft,p)! -(m_ES_WF_by_SDorCG) zaj_l, etc -> wfsd_l
!!$       !    ~~~~~~~~~~~~
!!$    end do
!!$
!!$    sumdz2 = 0.d0
!!$    do ib = ista_e, iend_e, istep_e      ! MPI
!!$       call square_of_SD_direction(ik,ib,dz)  ! wfsd_l dot wfsd_l
!!$       write(nfout,'(" dz = ",d20.8)') dz
!!$       sumdz2 = sumdz2 + dz
!!$    end do

    call make_CG_direction(ik)
!!$    write(6,'(" -- before orthonormalization -- ")')
!!$    call check_of_orthonormality(ik)  ! z = <wfsd_l(i)|zaj_l(j)>
!!$    call orthogonalize_SD_drctns(ik,to=OTHER_BANDS)     ! -(m_ES_WF_by_SDorCG) ->wfsd_l
    call cp_wfsd_to_wfsd_old(ik)                        ! -(m_ES_WF_by_SDorCG) wfsd_l ->wfsd_old
!!$    write(6,'(" -- after orthonormalization -- ")')
!!$    call check_of_orthonormality(ik)  ! z = <wfsd_l(i)|zaj_l(j)>

!!$    write(6,'(" -- before evolution -- ")')
!!$    call check_of_orthonormality2(ik)  ! <zaj_l(i)|zaj_l(j)>
    do ib = ista_e, iend_e, istep_e      ! MPI
       call WF_conjugate_gradient(ik,ib,dtim) !-(m_ES_WF_by_SDorCG)
    end do
    if(ipri>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')
!!$    write(6,'(" -- after evolution -- ")')
!!$    call check_of_orthonormality2(ik)  ! <zaj_l(i)|zaj_l(j)>

!!$    do ib = ista_e, iend_e, istep_e      ! MPI
!!$       call m_ES_WF_in_Rspace(ik,ib,bfft)! (swffft)
!!$       call m_FFT_Vlocal_W(afft,bfft)    ! (afft, bfft) -> (bfft)
!!$       call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
!!$       call WF_conjugate_gradient0(precon,ik,ib,dtim,ekin,bfft,p)
!!$       !    ~~~~~~~~~~~~~~~~~~~~~~         -(m_ES_WF_by_SDorCG)
!!$    end do
    if(modnrm == EXECUT) call dealloc_wfsd_bsdri()
    if(modnrm == EXECUT .and. precon == ON) call dealloc_wfsdnp_bsdrinp()
    call tstatc0_end(id_sname)
  contains
    subroutine get_betacg
      real(kind=DP) :: sumdz2, dz, avbetacg, betacg_min, betacg_max, avnewbetacg
      integer       :: ib

      if(ik == ista_k) then
         betacg_min = 1.d+10; betacg_max = -1.d+10; avbetacg = 0.d0
         avnewbetacg = 0.d0
      end if

      sumdz2 = 0.d0
      do ib = ista_e, iend_e, istep_e
         if(ipri >=1 ) write(6,'(" !! square_of_SD_direction(2?) is being called (evolve_WFs_in_eazyCG_direction.betacg)")')
         if(Precon == ON) then
            call square_of_SD_direction2(ik,ib,dz)
         else
            call square_of_SD_direction(ik,ib,dz)
         end if
         sumdz2 = sumdz2 + dz
      end do
      if(npes>1) call mpi_allreduce(MPI_IN_PLACE,sumdz2,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)

      if(dabs(dzajn2(ik)) > Delta) then
         betacg = sumdz2/dzajn2(ik)
      else
         betacg = 0.d0
      endif
      dzajn2(ik) = sumdz2

      if(ipri >= 1) write(nfout,'(" ik = ",i5," sumdz2 = ", d20.8, " betacg = ",d20.8)') ik,sumdz2,betacg
      if(betacg_min > betacg) betacg_min = betacg
      if(betacg_max < betacg) betacg_max = betacg
      avbetacg = avbetacg + betacg/kv3

      if(betacg > 1.d0) betacg = 1.d0
      if(betacg < 0.d0) betacg = 0.d0

      if(kv3 > 1 .and. ik == iend_k) then              ! MPI
         avnewbetacg = avnewbetacg + betacg/kv3
         if(ipri >= 1) write(nfout, '(" ! average of betacg = ",2f10.6," range ",f10.6,"-" &
              &,       f10.6)') avbetacg, avnewbetacg, betacg_min, betacg_max
      end if
    end subroutine get_betacg

    subroutine check_of_orthonormality(ik)
      integer, intent(in) :: ik
      integer             :: i, j, ig
      real(kind=DP)       :: wz

      if(k_symmetry(ik) == GAMMA) then
         do j = ista_e, iend_e, istep_e
            do i = ista_e, iend_e, istep_e
               wz = 0.d0
               do ig = 2, iba(ik)
                  wz = wz + wfsd_l(ig,i,ik,1)*zaj_l(ig,j,ik,1)
               end do
               wz = wz*2.d0 + wfsd_l(1,i,ik,1)*zaj_l(1,j,ik,1)
               if(ipri >= 1) write(6,'(" (i,j) = ",2i5, " wz = ",d20.8)') i,j, wz
            end do
         end do
      else
         do j = ista_e, iend_e, istep_e
            do i = ista_e, iend_e, istep_e
               wz = 0.d0
               do ig = 1, iba(ik)
                  wz = wz + wfsd_l(ig,i,ik,1)*zaj_l(ig,j,ik,1)
               end do
               if(ipri >= 1) write(6,'(" (i,j) = ",2i5, " wz = ",d20.8)') i,j, wz
            end do
         end do
      end if
    end subroutine check_of_orthonormality

    subroutine check_of_orthonormality2(ik)
      integer, intent(in) :: ik
      integer             :: i, j, ig
      real(kind=DP)       :: wz

      if(k_symmetry(ik) == GAMMA) then
         do j = ista_e, iend_e, istep_e
            do i = ista_e, iend_e, istep_e
               wz = 0.d0
               do ig = 2, iba(ik)
                  wz = wz + zaj_l(ig,i,ik,1)*zaj_l(ig,j,ik,1)
               end do
               wz = wz*2.d0 + zaj_l(1,i,ik,1)*zaj_l(1,j,ik,1)
               if(ipri >=1 ) write(6,'(" (i,j) = ",2i5, " zz = ",d20.8)') i,j, wz
            end do
         end do
      else
         do j = ista_e, iend_e, istep_e
            do i = ista_e, iend_e, istep_e
               wz = 0.d0
               do ig = 1, iba(ik)
                  wz = wz + zaj_l(ig,i,ik,1)*zaj_l(ig,j,ik,1)
               end do
               wz = wz*2.d0
               if(ipri >=1 ) write(6,'(" (i,j) = ",2i5, " zz = ",d20.8)') i,j, wz
            end do
         end do
      end if
    end subroutine check_of_orthonormality2
  end subroutine evolve_WFs_in_eazyCG_direction

! ============================ added by K. Tagami ========================= 11.0
  subroutine evolve_WFs_in_ezCG_direc_noncl( precon, ik, dtim, ekin,&
       &                                     afft_kt, bfft_kt, vnlph_noncl,  p )
    integer, intent(in)        :: precon, ik
    real(kind=DP), intent(in)  :: dtim,ekin(kg1)
    real(kind=DP), intent(inout)  :: afft_kt(nfft,ndim_chgpot)
    real(kind=DP), intent(out) :: bfft_kt(nfft,ndim_spinor)
    real(kind=DP), intent(in)  :: vnlph_noncl(kg1,np_e,kimg,ndim_spinor)
    real(kind=DP)              :: p(kg1)

    real(kind=DP), parameter   :: Delta = 1.d-10

    integer :: ib, is
    integer :: id_sname = -1
    call tstatc0_begin('evolve_WFs_in_eazyCG_dire_noncl ', id_sname,1)

    if (modnrm==EXECUT) call alloc_wfsd_bsdri_noncl(ik,ik+ndim_spinor-1)

    if (modnrm==EXECUT .and. precon==ON) then
	call alloc_wfsdnp_bsdrinp_noncl( ik,ik+ndim_spinor-1 )
    endif

    do ib = ista_e, iend_e, istep_e      ! MPI
       Do is=1, ndim_spinor
         call m_ES_WF_in_Rspace( ik+is-1, ib, bfft_kt(:,is) )  ! (swffft)
       End do

       call m_FFT_Vlocal_W_noncl( afft_kt, bfft_kt, nspin, ndim_spinor )

       Do is=1, ndim_spinor
         call m_FFT_WF( ELECTRON, nfout, bfft_kt(:,is), DIRECT, ON )
                                           ! bfft: G-space repres.
       End do
!
       call SD_direction_noncl( precon, ik, ib, ekin, bfft_kt, vnlph_noncl, p )

    end do

    if (precon == ON) then
       call Precon_dot_SD_noncl( ik, ekin, p )
                   ! -(m_ES_WF_by_SDorCG) -> wfsd_np, wfsd_l
    endif

    if (modnrm == EXECUT) then
       Do is=1, ndim_spinor
          call m_ES_betar_dot_Psi_4_each_k(nfout, wfsd_l, ik, ik+ndim_spinor-1, &
	&                                   ik+is-1, bsdr_l, bsdi_l)
       End do
    endif
    if (modnrm == EXECUT .and. precon == ON ) then
       Do is=1, ndim_spinor
          call m_ES_betar_dot_Psi_4_each_k(nfout, wfsd_np, ik, ik+ndim_spinor-1, &
	&                                   ik+is-1, bsdr_np, bsdi_np )
       End do
    end if

    call get_betacg_noncl
                                     ! -(contained here) wfsd_l ->betacg, dzajn2
    Do is=1, ndim_spinor
      call cp_wfsd_to_wfsd_old(ik+is-1)
                                   ! -(m_ES_WF_by_SDorCG) wfsd_l ->wfsd_old
      call make_CG_direction(ik+is-1)
    End do

    do ib = ista_e, iend_e, istep_e      ! MPI
       Do is=1, ndim_spinor
          call WF_conjugate_gradient( ik+is-1, ib, dtim ) !-(m_ES_WF_by_SDorCG)
       End do
    end do

    if(modnrm == EXECUT) call dealloc_wfsd_bsdri()
    if(modnrm == EXECUT .and. precon == ON) call dealloc_wfsdnp_bsdrinp()

    call tstatc0_end(id_sname)

  contains

    subroutine get_betacg_noncl
      real(kind=DP) :: sumdz2, dz, avbetacg, betacg_min, betacg_max, avnewbetacg
      integer       :: ib

      if(ik == ista_k) then                ! MPI
         betacg_min = 1.d+10; betacg_max = -1.d+10; avbetacg = 0.d0
         avnewbetacg = 0.d0
      end if

      sumdz2 = 0.d0
      do ib = ista_e, iend_e, istep_e
         if(ipri >=1 ) write(6,'(" !! square_of_SD_direction(2?) is being called (evolve_WFs_in_eazyCG_direction.betacg)")')
         if(Precon == ON) then
            call square_of_SD_direction2_noncl(ik,ib,dz)
         else
            call square_of_SD_direction_noncl(ik,ib,dz)
         end if
         sumdz2 = sumdz2 + dz
      end do

      if(npes>1) call mpi_allreduce(MPI_IN_PLACE,sumdz2,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)

      if(dabs(dzajn2(ik)) > Delta) then
         betacg = sumdz2/dzajn2(ik)
      else
         betacg = 0.d0
      endif
      dzajn2(ik) = sumdz2

      if(ipri >= 1) write(nfout,'(" ik = ",i5," sumdz2 = ", d20.8, " betacg = ",d20.8)') ik,sumdz2,betacg
      if(betacg_min > betacg) betacg_min = betacg
      if(betacg_max < betacg) betacg_max = betacg
      avbetacg = avbetacg + betacg / ( kv3/ndim_spinor )

      if(betacg > 1.d0) betacg = 1.d0
      if(betacg < 0.d0) betacg = 0.d0

      if(kv3 > 1 .and. ik == iend_k-ndim_spinor+1 )  then              ! MPI
         avnewbetacg = avnewbetacg + betacg/ ( kv3 /ndim_spinor )
         if(ipri >= 1) write(nfout, '(" ! average of betacg = ",2f10.6," range ",f10.6,"-" &
              &,       f10.6)') avbetacg, avnewbetacg, betacg_min, betacg_max
      end if
    end subroutine get_betacg_noncl

    subroutine check_of_orthonormality_noncl(ik)
      integer, intent(in) :: ik
      integer             :: i, j, ig
      real(kind=DP)       :: wz
      real(kind=DP)       :: wz1, wz2

      if(k_symmetry(ik) == GAMMA) then
         do j = ista_e, iend_e, istep_e
            do i = ista_e, iend_e, istep_e
               wz1 = 0.d0; wz2 = 0.0d0
               Do is=1, ndim_spinor
                 do ig = 2, iba(ik)
                    wz1 = wz1 + wfsd_l(ig,i,ik+is-1,1)*zaj_l(ig,j,ik+is-1,1)
                 end do
               End do
               Do is=1, ndim_spinor
                  wz2 = wz2 + wfsd_l(1,i,ik+is-1,1)*zaj_l(1,j,ik+is-1,1)
               End do
               wz = wz1*2.d0 + wz2
               if(ipri >= 1) write(6,'(" (i,j) = ",2i5, " wz = ",d20.8)') i,j, wz
            end do
         end do
      else
         do j = ista_e, iend_e, istep_e
            do i = ista_e, iend_e, istep_e
               wz = 0.d0
               Do is=1, ndim_spinor
                  do ig = 1, iba(ik)
                     wz = wz + wfsd_l(ig,i,ik+is-1,1)*zaj_l(ig,j,ik+is-1,1)
                  end do
               End do
               if(ipri >= 1) write(6,'(" (i,j) = ",2i5, " wz = ",d20.8)') i,j, wz
            end do
         end do
      end if
    end subroutine check_of_orthonormality_noncl

    subroutine check_of_orthonormality2_noncl(ik)
      integer, intent(in) :: ik
      integer             :: i, j, ig
      real(kind=DP)       :: wz
      real(kind=DP)       :: wz1, wz2

      if(k_symmetry(ik) == GAMMA) then
         do j = ista_e, iend_e, istep_e
            do i = ista_e, iend_e, istep_e
               wz1 = 0.d0; wz2 = 0.0d0
               Do is=1, ndim_spinor
                 do ig = 2, iba(ik)
                    wz1 = wz1 + zaj_l(ig,i,ik+is-1,1)*zaj_l(ig,j,ik+is-1,1)
                 end do
               End do
               Do is=1, ndim_spinor
                  wz2 = wz2 + zaj_l(1,i,ik+is-1,1)*zaj_l(1,j,ik+is-1,1)
               End do
               wz = wz1*2.d0 + wz2

               if(ipri >=1 ) write(6,'(" (i,j) = ",2i5, " zz = ",d20.8)') i,j, wz
            end do
         end do
      else
         do j = ista_e, iend_e, istep_e
            do i = ista_e, iend_e, istep_e
               wz = 0.d0
               Do is=1, ndim_spinor
                 do ig = 1, iba(ik)
                    wz = wz + zaj_l(ig,i,ik+is-1,1)*zaj_l(ig,j,ik+is-1,1)
                 end do
               End do
               wz = wz*2.d0
               if(ipri >=1 ) write(6,'(" (i,j) = ",2i5, " zz = ",d20.8)') i,j, wz
            end do
         end do
      end if
    end subroutine check_of_orthonormality2_noncl

  end subroutine evolve_WFs_in_ezCG_direc_noncl
! ===================================================================== 11.0

  subroutine evolve_WFs_in_CG_direction(precon,ik,dtim,ekin,afft,bfft,p)
    integer, intent(in)        :: precon, ik
    real(kind=DP), intent(in)  :: dtim,ekin(kg1)
    real(kind=DP), intent(in)  :: afft(nfft)
    real(kind=DP), intent(out) :: bfft(nfft)
    real(kind=DP)              :: p(kg1)

    integer :: ib
    integer :: id_sname = -1
    call tstatc0_begin('evolve_WFs_in_CG_direction ', id_sname,1)

    call alloc_wfsd_bsdri(ik)

    do ib = ista_e, iend_e, istep_e
       call m_ES_WF_in_Rspace(ik,ib,bfft)! (swffft)
       call m_FFT_Vlocal_W(afft,bfft)    ! (afft, bfft) -> (bfft)
       call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON) ! bfft: G-space repres.
       call SD_direction(precon,ik,ib,ekin,bfft,p)     ! -(m_ES_WF_by_SDorCG) -> wfsd_l
       !                wfsd_l : -(H-e_{k\mu}^m)\Psi_{k\mu}^m
    end do

    call orthogonalize_SD_drctns(ik,to=SAME_BAND) ! -(m_ES_WF_by_SDorCG) ->(wfsd_l,bsd(ri)_l)
    call make_CG_direction(ik)       ! -(m_ES_WF_by_SDorCG) -> wfsd_l + betacg*wfsd_old -> wfsd_l
!!$    call orthogonalize_SD_drctns(ik,to=SAME_BAND) ! Actually, CG direction ->(wfsd_l,bsd(ri)_l)

    do ib = ista_e, iend_e, istep_e
       call WF_conjugate_gradient(ik,ib,dtim) !-(m_ES_WF_by_SDorCG)
    end do

    call cp_wfsd_to_wfsd_old(ik)    ! -(m_ES_WF_by_SDorCG) wfsd_l ->wfsd_old

    call dealloc_wfsd_bsdri()
    call tstatc0_end(id_sname)
  end subroutine evolve_WFs_in_CG_direction

! ======================================= added by K. Tagami ================= 11.0
  subroutine evolve_WFs_in_CG_direc_noncl( precon, ik, dtim, ekin, &
       &                                   afft_kt, bfft_kt, vnlph_noncl, p)
    integer, intent(in)        :: precon, ik
    real(kind=DP), intent(in)  :: dtim,ekin(kg1)
    real(kind=DP), intent(inout)  :: afft_kt(nfft,ndim_chgpot)
    real(kind=DP), intent(out) :: bfft_kt(nfft,ndim_spinor)
    real(kind=DP), intent(in)  :: vnlph_noncl(kg1,np_e,kimg,ndim_spinor)
    real(kind=DP)              :: p(kg1)

    integer :: ib, is
    integer :: id_sname = -1
    call tstatc0_begin('evolve_WFs_in_CG_direc_noncl ', id_sname,1)

    call alloc_wfsd_bsdri_noncl( ik,ik+ndim_spinor-1 )

    do ib = ista_e, iend_e, istep_e

       Do is=1, ndim_spinor
          call m_ES_WF_in_Rspace( ik+is-1, ib, bfft_kt(:,is) )   ! (swffft)
       End do

       call m_FFT_Vlocal_W_noncl( afft_kt, bfft_kt, ndim_chgpot, ndim_spinor )

       Do is=1, ndim_spinor
         call m_FFT_WF( ELECTRON, nfout, bfft_kt(:,is), DIRECT, ON )
       End do

       call SD_direction_noncl( precon, ik, ib, ekin, bfft_kt, vnlph_noncl, p )
                                 ! -(m_ES_WF_by_SDorCG) -> wfsd_l
                                  !        wfsd_l : -(H-e_{k\mu}^m)\Psi_{k\mu}^m

    end do

    call orthogonalize_SD_drctns_noncl( ik,to=SAME_BAND )
                       ! -(m_ES_WF_by_SDorCG) ->(wfsd_l,bsd(ri)_l)

    Do is=1, ndim_spinor
      call make_CG_direction(ik+is-1)
                  ! -(m_ES_WF_by_SDorCG) -> wfsd_l + betacg*wfsd_old -> wfsd_l
    End do
    do ib = ista_e, iend_e, istep_e           ! MPI
       Do is=1, ndim_spinor
          call WF_conjugate_gradient( ik+is-1, ib, dtim )   !-(m_ES_WF_by_SDorCG)
       End do
    end do

    Do is=1, ndim_spinor
       call cp_wfsd_to_wfsd_old( ik+is-1 )
                                       ! -(m_ES_WF_by_SDorCG) wfsd_l ->wfsd_old
    End do

    call dealloc_wfsd_bsdri()
    call tstatc0_end(id_sname)
  end subroutine evolve_WFs_in_CG_direc_noncl
! ================================================================= 11.0

  subroutine evolve_WFs_in_CG_direction0 &
       &(precon,ik,dtim,ekin,afft,bfft,p)
    integer, intent(in)        :: precon, ik
    real(kind=DP), intent(in)  :: dtim,ekin(kg1)
    real(kind=DP), intent(in)  :: afft(nfft)
    real(kind=DP), intent(out) :: bfft(nfft)
    real(kind=DP)              :: p(kg1)

    integer :: ib
    integer :: id_sname = -1
    call tstatc0_begin('evolve_WFs_in_CG_direction0 ', id_sname,1)

    do ib = ista_e, iend_e, istep_e
       call m_ES_WF_in_Rspace(ik,ib,bfft)! (swffft)
       call m_FFT_Vlocal_W(afft,bfft)    ! (afft, bfft) -> (bfft)
       call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
       call WF_conjugate_gradient0(precon,ik,ib,dtim,ekin,bfft,p) !-(m_ES_WF_by_SDorCG)
    end do
    if(ipri>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')
    call tstatc0_end(id_sname)
  end subroutine evolve_WFs_in_CG_direction0

  subroutine vlhxc_l_zero_term(vlhxc0,ispin)
    real(kind=DP), intent(out) :: vlhxc0
    integer, intent(in)        :: ispin

    if(mype == 0) vlhxc0 = vlhxc_l(1,1,ispin)
    call mpi_bcast(vlhxc0,1,mpi_double_precision,0,MPI_CommGroup,ierr)
  end subroutine vlhxc_l_zero_term

! ===================================== added by K. Tagami ============== 11.0
  subroutine vlhxc_l_zero_term_noncl( vlhxc0 )
    real(kind=DP), intent(out) :: vlhxc0(ndim_chgpot)

    vlhxc0 = 0.0d0
!
    if (mype == 0) then
! ------------------------------ uncertain ----
!        vlhxc0(1)           = vlhxc_l(1,1,1) + vlhxc_l(1,1,ndim_magmom)
!        vlhxc0(ndim_chgpot) = vlhxc_l(1,1,1) - vlhxc_l(1,1,ndim_magmom)
!
       vlhxc0 = vlhxc_l(1,1,1)
! --------------------------
    endif
    call mpi_bcast( vlhxc0, ndim_chgpot, mpi_double_precision,0,MPI_CommGroup,ierr)
  end subroutine vlhxc_l_zero_term_noncl
! ========================================================================= 11.0

  subroutine modified_steepest_descent&
       &(precon,ik,ibo,dtim,vnldi,vlhxc0,ekin,VlocalW,p)
    integer      , intent(in)                  :: precon,ik,ibo
    real(kind=DP), intent(in)                  :: dtim
    real(kind=DP), intent(in), dimension(kg1)  :: vnldi
    real(kind=DP), intent(in)                  :: vlhxc0
    real(kind=DP), intent(in), dimension(kg1)  :: ekin
    real(kind=DP), intent(in), dimension(nfft) :: VlocalW
    real(kind=DP)            , dimension(kg1)  :: p

    integer       :: i, i1, ib
    real(kind=DP) :: evr,devr,denom, wdi, evi,e1, devi, fdexp
!!$    integer :: id_sname = -1
!!$    call tstatc0_begin('modified_steepest_descent ', id_sname)

    ib = map_z(ibo)                            ! MPI
    denom = 1.d0/product(fft_box_size_WF(1:3,1))
    call m_ES_decide_precon_factor(precon,ik,ibo,ekin,p)  ! -> p(1:iba(ik))
    if(ipri >= 3) then
       write(nfout,'(" ekin : ",5d16.8)') (ekin(i),i=1,iba(ik))
       write(nfout,'(" p    : ",5d16.8)') (p(i),i=1,iba(ik))
    end if

#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) status_saved_phifftr(ib,ik) = OLD
#endif
    if(kimg == 1) then
       do i = 1, iba(ik)
          i1    = igf(nbase(i,ik))
          evr   = zaj_l(i,ib,ik,1)
          devr  = (ekin(i)-eko_l(ib,ik))*evr&
               & + VlocalW(i1)*denom + vnlph_l(i,ib,1)
          wdi   = ekin(i) + vlhxc0 + vnldi(i) - eko_l(ib,ik)
! === 0 divide occurs if abs(wdi) < epsilon. by tkato 2012/12/18 ===============
!         fdexp = dexp( -p(i) * wdi * dtim)
!         zaj_l(i,ib,ik,1) = (fdexp - 1)*devr/wdi + evr
          if(dabs(wdi) < SmallestPositiveNumber) then
             zaj_l(i,ib,ik,1) = -p(i)*devr*dtim + evr
          else
             fdexp = dexp( -p(i) * wdi * dtim)
             zaj_l(i,ib,ik,1) = (fdexp - 1)*devr/wdi + evr
          end if
! ==============================================================================
       end do
    else if(kimg == 2) then
       do i = 1, iba(ik)
          i1    = igf(nbase(i,ik))
          evr   = zaj_l(i,ib,ik,1);    evi   = zaj_l(i,ib,ik,2)
          e1    = ekin(i) - eko_l(ib,ik)
          devr  = e1*evr+VlocalW(2*i1-1)*denom+vnlph_l(i,ib,1)
          devi  = e1*evi+VlocalW(2*i1  )*denom+vnlph_l(i,ib,2)
          wdi   = ekin(i) + vlhxc0 + vnldi(i) - eko_l(ib,ik)
! === 0 divide occurs if abs(wdi) < epsilon. by tkato 2012/12/18 ===============
!         fdexp = dexp( -p(i) * wdi * dtim)
!         zaj_l(i,ib,ik,1) = (fdexp - 1)*devr/wdi + evr
!         zaj_l(i,ib,ik,2) = (fdexp - 1)*devi/wdi + evi
          if(dabs(wdi) < SmallestPositiveNumber) then
             zaj_l(i,ib,ik,1) = -p(i)*devr*dtim + evr
             zaj_l(i,ib,ik,2) = -p(i)*devi*dtim + evi
          else
             fdexp = dexp( -p(i) * wdi * dtim)
             zaj_l(i,ib,ik,1) = (fdexp - 1)*devr/wdi + evr
             zaj_l(i,ib,ik,2) = (fdexp - 1)*devi/wdi + evi
          end if
! ==============================================================================
       end do
    end if
!!$    call tstatc0_end(id_sname)
  end subroutine modified_steepest_descent

! ================================ added by K. Tagami =================== 11.0
  subroutine modified_steepest_desc_noncl( precon, ik, ibo, dtim, &
       &                                   vnldi_noncl, vlhxc0, ekin, &
       &                                   VlocalW_noncl, vnlph_noncl, p )
    integer      , intent(in)                  :: precon,ik,ibo
    real(kind=DP), intent(in)                  :: dtim

#ifdef USE_COMPLEX_VNLDI
    complex(kind=CMPLDP), intent(in), dimension(kg1,ndim_spinor)  :: vnldi_noncl
#else
    real(kind=DP), intent(in), dimension(kg1,ndim_spinor)  :: vnldi_noncl
#endif

    real(kind=DP), intent(in)                  :: vlhxc0(ndim_chgpot)
    real(kind=DP), intent(in), dimension(kg1)  :: ekin
    real(kind=DP), intent(in), dimension(nfft,ndim_spinor) :: VlocalW_noncl
    real(kind=DP), intent(in), dimension(kg1,np_e,kimg,ndim_spinor) :: vnlph_noncl

    real(kind=DP)            , dimension(kg1)  :: p

    integer       :: i, i1, ib
    integer  :: is, k1, istmp

#ifdef USE_COMPLEX_VNLDI
    real(kind=DP) :: evr,devr,denom, evi, e1, devi
    complex(kind=CMPLDP) :: zwdi, zfdexp, ztmp1
    real(kind=DP) :: ctmp1, stmp1
#else
    real(kind=DP) :: evr,devr,denom, wdi, evi,e1, devi, fdexp
#endif

    ib = map_z(ibo)                            ! MPI
    denom = 1.d0/product(fft_box_size_WF(1:3,1))

    call m_ES_decide_precon_factor(precon,ik,ibo,ekin,p)  ! -> p(1:iba(ik))

    if(ipri >= 3) then
       write(nfout,'(" ekin : ",5d16.8)') (ekin(i),i=1,iba(ik))
       write(nfout,'(" p    : ",5d16.8)') (p(i),i=1,iba(ik))
    end if

#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) then
       do is = 1, ndim_spinor
          k1 = ik+is-1
          status_saved_phifftr(ib,k1) = OLD
       end do
    end if
#endif
#ifdef USE_COMPLEX_VNLDI
    if(kimg == 1) then

    else if(kimg == 2) then

       Do is=1, ndim_spinor
         k1 = ik +is -1
         istmp = ( is -1 )*ndim_spinor + is

         do i = 1, iba(ik)
            i1    = igf(nbase(i,ik))
            evr   = zaj_l(i,ib,k1,1);    evi   = zaj_l(i,ib,k1,2)
            e1    = ekin(i) - eko_l(ib,ik)

            devr  = e1*evr+VlocalW_noncl(2*i1-1,is)*denom +vnlph_noncl(i,ib,1,is)
            devi  = e1*evi+VlocalW_noncl(2*i1,  is)*denom +vnlph_noncl(i,ib,2,is)

            zwdi   = ekin(i) + vlhxc0(istmp) + vnldi_noncl(i,is) - eko_l(ib,ik)

            ztmp1 = -p(i) * zwdi * dtim
            zfdexp = exp( ztmp1 )

            ztmp1 = ( zfdexp - 1.0d0 ) / zwdi
            ztmp1 = ztmp1 *cmplx( devr, devi )

            ctmp1 = real(ztmp1);  stmp1 = aimag(ztmp1)

            zaj_l(i,ib,k1,1) = ctmp1 + evr
            zaj_l(i,ib,k1,2) = stmp1 + evi
          End do
       end do
    end if

#else
    if(kimg == 1) then
       Do is=1, ndim_spinor
         k1 = ik +is -1
         istmp = ( is -1 )*ndim_spinor + is

         do i = 1, iba(ik)
            i1    = igf(nbase(i,ik))
            evr   = zaj_l(i,ib,k1,1)
            devr  = (ekin(i)-eko_l(ib,ik))*evr&
               & + VlocalW_noncl(i1,is)*denom + vnlph_noncl(i,ib,1,is)
            wdi   = ekin(i) + vlhxc0(istmp) + vnldi_noncl(i,is) - eko_l(ib,ik)

            fdexp = dexp( -p(i) * wdi * dtim)
            zaj_l(i,ib,k1,1) = (fdexp - 1)*devr/wdi + evr
          End do
       end do
    else if(kimg == 2) then
       Do is=1, ndim_spinor
         k1 = ik +is -1
         istmp = ( is -1 )*ndim_spinor + is

         do i = 1, iba(ik)
            i1    = igf(nbase(i,ik))
            evr   = zaj_l(i,ib,k1,1);    evi   = zaj_l(i,ib,k1,2)
            e1    = ekin(i) - eko_l(ib,ik)

            devr  = e1*evr+VlocalW_noncl(2*i1-1,is)*denom +vnlph_noncl(i,ib,1,is)
            devi  = e1*evi+VlocalW_noncl(2*i1,  is)*denom +vnlph_noncl(i,ib,2,is)

            wdi   = ekin(i) + vlhxc0(istmp) + vnldi_noncl(i,is) - eko_l(ib,ik)

            fdexp = dexp( -p(i) * wdi * dtim)

            zaj_l(i,ib,k1,1) = (fdexp - 1)*devr/wdi + evr
            zaj_l(i,ib,k1,2) = (fdexp - 1)*devi/wdi + evi
          End do
       end do
    end if

#endif

  end subroutine modified_steepest_desc_noncl
! ================================================================== 11.0

  subroutine steepest_descent(precon,ik,ibo,dtim,ekin,VlocalW,p)
    integer      , intent(in)                  :: precon,ik,ibo
    real(kind=DP), intent(in)                  :: dtim
    real(kind=DP), intent(in), dimension(kg1)  :: ekin
    real(kind=DP), intent(in), dimension(nfft) :: VlocalW
    real(kind=DP)            , dimension(kg1)  :: p

    integer       :: i, i1, ib
    real(kind=DP) :: evr,devr,denom, evi,e1, devi

    ib = map_z(ibo)                            ! MPI
    denom = 1.d0/product(fft_box_size_WF(1:3,1))
    call m_ES_decide_precon_factor(precon,ik,ibo,ekin,p)  ! -> p(1:iba(ik))

#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) status_saved_phifftr(ib,ik) = OLD
#endif
    if(kimg == 1) then
       do i = 1, iba(ik)
          i1    = igf(nbase(i,ik))
          evr   = zaj_l(i,ib,ik,1)
          devr  = (ekin(i)-eko_l(ib,ik))*evr + VlocalW(i1)*denom + vnlph_l(i,ib,1)
          zaj_l(i,ib,ik,1) = evr - p(i)*dtim*devr
       end do
    else if(kimg == 2) then
       do i = 1, iba(ik)
          i1    = igf(nbase(i,ik))
          evr   = zaj_l(i,ib,ik,1);    evi   = zaj_l(i,ib,ik,2)
          e1    = ekin(i) - eko_l(ib,ik)
          devr  = e1*evr+VlocalW(2*i1-1)*denom+vnlph_l(i,ib,1)
          devi  = e1*evi+VlocalW(2*i1  )*denom+vnlph_l(i,ib,2)
          zaj_l(i,ib,ik,1) = evr - p(i)*dtim*devr
          zaj_l(i,ib,ik,2) = evi - p(i)*dtim*devi
       end do
    end if
  end subroutine steepest_descent

! ================================== added by K. Tagami ======================= 11.0
  subroutine steepest_descent_noncl( precon, ik, ibo, dtim, ekin, VlocalW_noncl, &
       &                             vnlph_noncl, p )
    integer      , intent(in)                  :: precon,ik,ibo
    real(kind=DP), intent(in)                  :: dtim
    real(kind=DP), intent(in), dimension(kg1)  :: ekin
    real(kind=DP), intent(in), dimension(nfft,ndim_spinor) :: VlocalW_noncl
    real(kind=DP), intent(in), dimension(kg1,np_e,kimg,ndim_spinor) :: vnlph_noncl

    real(kind=DP)            , dimension(kg1)  :: p

    integer       :: i, i1, ib
    integer :: is, istmp, k1
    real(kind=DP) :: evr,devr,denom, evi,e1, devi

    ib = map_z(ibo)                            ! MPI
    denom = 1.d0/product(fft_box_size_WF(1:3,1))
    call m_ES_decide_precon_factor(precon,ik,ibo,ekin,p)  ! -> p(1:iba(ik))
#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) then
       do is = 1, ndim_spinor
          k1 = ik+is-1
          status_saved_phifftr(ib,k1) = OLD
       end do
    end if
#endif
    if(kimg == 1) then
       Do is=1, ndim_spinor
          k1 = ik +is -1
          istmp = ( is -1 )*ndim_spinor + is

          do i = 1, iba(ik)
             i1    = igf(nbase(i,ik))
             evr   = zaj_l(i,ib,k1,1)
             devr  = ( ekin(i)-eko_l(ib,ik))*evr + VlocalW_noncl(i1,is)*denom &
	&           + vnlph_noncl(i,ib,1,is)
             zaj_l(i,ib,k1,1) = evr - p(i)*dtim*devr
          End do
       end do
    else if(kimg == 2) then
       Do is=1, ndim_spinor
          k1 = ik +is -1
          istmp = ( is -1 )*ndim_spinor + is

          do i = 1, iba(ik)
             i1    = igf(nbase(i,ik))
             evr   = zaj_l(i,ib,k1,1);    evi   = zaj_l(i,ib,k1,2)
             e1    = ekin(i) - eko_l(ib,ik)
             devr  = e1*evr+VlocalW_noncl(2*i1-1,is)*denom+vnlph_noncl(i,ib,1,is)
             devi  = e1*evi+VlocalW_noncl(2*i1,  is)*denom+vnlph_noncl(i,ib,2,is)
             zaj_l(i,ib,k1,1) = evr - p(i)*dtim*devr
             zaj_l(i,ib,k1,2) = evi - p(i)*dtim*devi
          End do
       end do
    end if
  end subroutine steepest_descent_noncl
! ============================================================================ 11.0

  subroutine square_of_SD_direction(ik,ibo,dz)
    integer,       intent(in)       :: ik,ibo
    real(kind=DP), intent(out)      :: dz  ! dz = <wfsd|w|wfsd>
    integer                         :: i,ib

    dz = 0.d0
    if(modnrm == EXECUT) then
       call m_ES_sum_of_LocalPart(ik,ibo,bsdr_l,bsdi_l,dz)
    end if
    ib = map_z(ibo)                     ! MPI
    if(k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          dz = dz + wfsd_l(1,ib,ik,1)*wfsd_l(1,ib,ik,1)
          do i = 2, iba(ik)
             dz = dz + 2.d0*wfsd_l(i,ib,ik,1)*wfsd_l(i,ib,ik,1)
          end do
       else
          dz = dz + wfsd_l(1,ib,ik,1)**2 + wfsd_l(1,ib,ik,2)**2
          do i = 2, iba(ik)
             dz = dz + 2.d0*(wfsd_l(i,ib,ik,1)**2 + wfsd_l(i,ib,ik,2)**2)
          end do
       end if
    else
       if(kimg == 1) then
          do i = 1, iba(ik)
             dz = dz + wfsd_l(i,ib,ik,1)*wfsd_l(i,ib,ik,1)
          end do
       else
          do i = 1, iba(ik)
             dz = dz + wfsd_l(i,ib,ik,1)**2 + wfsd_l(i,ib,ik,2)**2
          end do
       end if
    end if
  end subroutine square_of_SD_direction

! ==================================== added by K. Tagami ================== 11.0
  subroutine square_of_SD_direction_noncl( ik,ibo,dz )
    integer,       intent(in)       :: ik,ibo
    real(kind=DP), intent(out)      :: dz  ! dz = <wfsd|w|wfsd>
    integer                         :: i,ib
    integer :: is, k1

    dz = 0.d0
    if(modnrm == EXECUT) then
       call m_ES_sum_of_LocalPart_noncl( ik, ibo, ik, ik+ndim_spinor-1, &
     &                                   bsdr_l, bsdi_l, dz )
    end if

    ib = map_z(ibo)                     ! MPI
    if(k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          Do is=1, ndim_spinor
            k1 = ik + is -1
            dz = dz + wfsd_l(1,ib,k1,1)*wfsd_l(1,ib,k1,1)
          End do
          Do is=1, ndim_spinor
            k1 = ik + is -1
            do i = 2, iba(ik)
              dz = dz + 2.d0*wfsd_l(i,ib,k1,1)*wfsd_l(i,ib,k1,1)
            end do
          End do
       else
          Do is=1, ndim_spinor
            k1 = ik + is -1
            dz = dz + wfsd_l(1,ib,k1,1)**2 + wfsd_l(1,ib,k1,2)**2
          End do
          Do is=1, ndim_spinor
            do i = 2, iba(ik)
               k1 = ik + is -1
               dz = dz + 2.d0*(wfsd_l(i,ib,k1,1)**2 + wfsd_l(i,ib,k1,2)**2)
            end do
          End do
       end if
    else
       if(kimg == 1) then
          Do is=1, ndim_spinor
            do i = 1, iba(ik)
               k1 = ik + is -1
               dz = dz + wfsd_l(i,ib,k1,1)*wfsd_l(i,ib,k1,1)
            end do
          End do
       else
          Do is=1, ndim_spinor
            do i = 1, iba(ik)
               k1 = ik + is -1
               dz = dz + wfsd_l(i,ib,k1,1)**2 + wfsd_l(i,ib,k1,2)**2
            end do
          End do
       end if
    end if
  end subroutine square_of_SD_direction_noncl
! ======================================================================== 11.0

  subroutine square_of_SD_direction2(ik,ibo,dz)
    integer,       intent(in)       :: ik,ibo
    real(kind=DP), intent(out)      :: dz  ! dz = <wfsd|w|wfsd>
    integer                         :: i,ib

    ib = map_z(ibo)                     ! MPI
    dz = 0.d0
    if(modnrm == EXECUT) then
       call m_ES_sum_of_LocalPart2(ik,ibo,bsdr_l,bsdi_l,bsdr_np,bsdi_np,dz)
    end if
    if(k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          dz = dz + wfsd_l(1,ib,ik,1)*wfsd_np(1,ib,ik,1)
          do i = 2, iba(ik)
             dz = dz + 2.d0*wfsd_l(i,ib,ik,1)*wfsd_np(i,ib,ik,1)
          end do
       else
          dz = dz + wfsd_l(1,ib,ik,1)*wfsd_np(1,ib,ik,1)
          do i = 2, iba(ik)
             dz = dz + 2.d0*( wfsd_l(i,ib,ik,1)*wfsd_np(i,ib,ik,1) &
                  &         + wfsd_l(i,ib,ik,2)*wfsd_np(i,ib,ik,2))
          end do
       end if
    else
       if(kimg == 1) then
          do i = 1, iba(ik)
             dz = dz + wfsd_l(i,ib,ik,1)*wfsd_np(i,ib,ik,1)
          end do
       else
          do i = 1, iba(ik)
             dz = dz + wfsd_l(i,ib,ik,1)*wfsd_np(i,ib,ik,1) &
                  &  + wfsd_l(i,ib,ik,2)*wfsd_np(i,ib,ik,2)
          end do
       end if
    end if
  end subroutine square_of_SD_direction2

! =================================== added by K. Tagami =================== 11.0
  subroutine square_of_SD_direction2_noncl( ik,ibo,dz )
    integer,       intent(in)       :: ik,ibo
    real(kind=DP), intent(out)      :: dz  ! dz = <wfsd|w|wfsd>
    integer                         :: i,ib, is, k1

    ib = map_z(ibo)                     ! MPI
    dz = 0.d0
    if(modnrm == EXECUT) then
       call m_ES_sum_of_LocalPart2_noncl( ik, ibo, ik, ik+ndim_spinor-1, &
	&                                 bsdr_l, bsdi_l, bsdr_np, bsdi_np, dz )
    end if
    if(k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          Do is=1, ndim_spinor
            k1 = ik + is -1
            dz = dz + wfsd_l(1,ib,k1,1)*wfsd_np(1,ib,k1,1)
          End do
          Do is=1, ndim_spinor
            k1 = ik + is -1
            do i = 2, iba(ik)
               dz = dz + 2.d0*wfsd_l(i,ib,k1,1)*wfsd_np(i,ib,k1,1)
            end do
          End do
       else
          Do is=1, ndim_spinor
            k1 = ik + is -1
            dz = dz + wfsd_l(1,ib,k1,1)*wfsd_np(1,ib,k1,1) &
                 &  + wfsd_l(1,ib,k1,2)*wfsd_np(1,ib,k1,2)
          End do
          Do is=1, ndim_spinor
            k1 = ik + is -1
            do i = 2, iba(ik)
               dz = dz + 2.d0*( wfsd_l(i,ib,k1,1)*wfsd_np(i,ib,k1,1) &
                    &         + wfsd_l(i,ib,k1,2)*wfsd_np(i,ib,k1,2))
            end do
          End do
       end if
    else
       if(kimg == 1) then
          Do is=1, ndim_spinor
            k1 = ik + is -1
            do i = 1, iba(ik)
               dz = dz + wfsd_l(i,ib,k1,1)*wfsd_np(i,ib,k1,1)
            end do
          End do
       else
          Do is=1, ndim_spinor
            k1 = ik + is -1
            do i = 1, iba(ik)
               dz = dz + wfsd_l(i,ib,k1,1)*wfsd_np(i,ib,k1,1) &
                    &  + wfsd_l(i,ib,k1,2)*wfsd_np(i,ib,k1,2)
            end do
          End do
       end if
    end if
  end subroutine square_of_SD_direction2_noncl
! ========================================================================== 11.0
!fj$$#endif

  subroutine square_of_SD_direction3(ik,ibo,dz)
    integer,       intent(in)       :: ik,ibo
    real(kind=DP), intent(out)      :: dz  ! dz = <wfsd|w|wfsd>
    integer                         :: i,ib

    dz = 0.d0
    if(modnrm == EXECUT) then
       call m_ES_sum_of_LocalPart3(ik,ibo,bsdr_l,bsdi_l,bsdr_old,bsdi_old,dz)
    end if
    ib = map_z(ibo)                     ! MPI
    if(k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          dz = dz + wfsd_l(1,ib,ik,1)*wfsd_old(1,ib,ik,1)
          do i = 2, iba(ik)
             dz = dz + 2.d0*wfsd_l(i,ib,ik,1)*wfsd_old(i,ib,ik,1)
          end do
       else
          dz = dz + wfsd_l(1,ib,ik,1)*wfsd_old(1,ib,ik,1) + wfsd_l(1,ib,ik,2)*wfsd_old(1,ib,ik,2)
          do i = 2, iba(ik)
             dz = dz + 2.d0*(wfsd_l(i,ib,ik,1)*wfsd_old(i,ib,ik,1) + wfsd_l(i,ib,ik,2)*wfsd_old(i,ib,ik,2))
          end do
       end if
    else
       if(kimg == 1) then
          do i = 1, iba(ik)
             dz = dz + wfsd_l(i,ib,ik,1)*wfsd_old(i,ib,ik,1)
          end do
       else
          do i = 1, iba(ik)
             dz = dz + wfsd_l(i,ib,ik,1)*wfsd_old(i,ib,ik,1) + wfsd_l(i,ib,ik,2)*wfsd_old(i,ib,ik,2)
          end do
       end if
    end if
  end subroutine square_of_SD_direction3

! =================================== added by K. Tagami =================== 11.0
  subroutine square_of_SD_direction3_noncl( ik,ibo,dz )
    integer,       intent(in)       :: ik,ibo
    real(kind=DP), intent(out)      :: dz  ! dz = <wfsd|w|wfsd>
    integer                         :: i,ib, is, k1

    dz = 0.d0
    if (modnrm == EXECUT) then
       call m_ES_sum_of_LocalPart3_noncl( ik, ibo, ik, ik+ndim_spinor-1, &
	&                                 bsdr_l, bsdi_l, bsdr_old, bsdi_old, dz )
    end if
    ib = map_z(ibo)                     ! MPI

    if(k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          Do is=1, ndim_spinor
             k1 = ik + is -1
             dz = dz + wfsd_l(1,ib,k1,1) *wfsd_old(1,ib,k1,1)
          End do
          Do is=1, ndim_spinor
             k1 = ik + is -1
             do i = 2, iba(ik)
                dz = dz + 2.d0*wfsd_l(i,ib,k1,1)*wfsd_old(i,ib,k1,1)
             end do
          End do
       else
          Do is=1, ndim_spinor
             k1 = ik + is -1
             dz = dz + wfsd_l(1,ib,k1,1)*wfsd_old(1,ib,k1,1) &
                  &  + wfsd_l(1,ib,k1,2)*wfsd_old(1,ib,k1,2)
          End do
          Do is=1, ndim_spinor
             k1 = ik + is -1
             do i = 2, iba(ik)
                dz = dz + 2.d0*( wfsd_l(i,ib,k1,1)*wfsd_old(i,ib,k1,1) &
                     &         + wfsd_l(i,ib,k1,2)*wfsd_old(i,ib,k1,2))
             end do
          End do
       end if
    else
       if(kimg == 1) then
          Do is=1, ndim_spinor
             k1 = ik + is -1
             do i = 1, iba(ik)
                dz = dz + wfsd_l(i,ib,k1,1)*wfsd_old(i,ib,k1,1)
             end do
          End do
       else
          Do is=1, ndim_spinor
             k1 = ik + is -1
             do i = 1, iba(ik)
                dz = dz + wfsd_l(i,ib,k1,1)*wfsd_old(i,ib,k1,1) &
                     &  + wfsd_l(i,ib,k1,2)*wfsd_old(i,ib,k1,2)
             end do
          End do
       end if
    end if
  end subroutine square_of_SD_direction3_noncl
! =========================================================================== 11.0
  subroutine SD_direction(precon,ik,ibo,ekin,VlocalW,p)
    integer      , intent(in)                  :: precon,ik,ibo
    real(kind=DP), intent(in), dimension(kg1)  :: ekin
    real(kind=DP), intent(in), dimension(nfft) :: VlocalW
    real(kind=DP)            , dimension(kg1)  :: p

    integer       :: i, i1, ib
    real(kind=DP) :: devr,denom, e1, devi

    ib = map_z(ibo)                                  ! MPI
    denom = 1.d0/product(fft_box_size_WF(1:3,1))
    call m_ES_decide_precon_factor(precon,ik,ibo,ekin,p)  ! -> p(1:iba(ik))

    if(kimg == 1) then
       do i = 1, iba(ik)
          i1    = igf(nbase(i,ik))
          devr  = (ekin(i)-eko_l(ib,ik))*zaj_l(i,ib,ik,1)&
               & + VlocalW(i1)*denom + vnlph_l(i,ib,1)
          wfsd_l(i,ib,ik,1) = - p(i)*devr
       end do
    else if(kimg == 2) then
       do i = 1, iba(ik)
          i1    = igf(nbase(i,ik))
          e1    = ekin(i) - eko_l(ib,ik)
          devr  = e1*zaj_l(i,ib,ik,1) + VlocalW(2*i1-1)*denom+vnlph_l(i,ib,1)
          devi  = e1*zaj_l(i,ib,ik,2) + VlocalW(2*i1  )*denom+vnlph_l(i,ib,2)
          wfsd_l(i,ib,ik,1) = - p(i)*devr
          wfsd_l(i,ib,ik,2) = - p(i)*devi
       end do
    end if
  end subroutine SD_direction

! =================================== added by K. Tagami ================= 11.0
  subroutine SD_direction_noncl( precon, ik, ibo, ekin, VlocalW_noncl, &
       &                         vnlph_noncl, p )
    integer      , intent(in)                  :: precon,ik,ibo
    real(kind=DP), intent(in), dimension(kg1)  :: ekin
    real(kind=DP), intent(in), dimension(nfft,ndim_spinor) :: VlocalW_noncl
    real(kind=DP), intent(in), dimension(kg1,np_e,kimg,ndim_spinor) :: vnlph_noncl

    real(kind=DP)            , dimension(kg1)  :: p

    integer       :: i, i1, ib
    integer :: is, k1
    real(kind=DP) :: devr,denom, e1, devi

    ib = map_z(ibo)                                  ! MPI
    denom = 1.d0/product(fft_box_size_WF(1:3,1))
    call m_ES_decide_precon_factor(precon,ik,ibo,ekin,p)  ! -> p(1:iba(ik))

    if(kimg == 1) then
       Do is=1, ndim_spinor
         k1 = ik + is -1
         do i = 1, iba(ik)
            i1    = igf(nbase(i,ik))
            devr  = (ekin(i)-eko_l(ib,ik))*zaj_l(i,ib,k1,1)&
               &     + VlocalW_noncl(i1,is)*denom + vnlph_noncl(i,ib,1,is)
            wfsd_l(i,ib,k1,1) = - p(i)*devr
         end do
       End do
    else if(kimg == 2) then
       Do is=1, ndim_spinor
         k1 = ik + is -1
         do i = 1, iba(ik)
            i1    = igf(nbase(i,ik))
            e1    = ekin(i) - eko_l(ib,ik)
            devr  = e1*zaj_l(i,ib,k1,1) + VlocalW_noncl(2*i1-1,is)*denom &
	&           +vnlph_noncl(i,ib,1,is)
            devi  = e1*zaj_l(i,ib,k1,2) + VlocalW_noncl(2*i1,  is)*denom &
	&           +vnlph_noncl(i,ib,2,is)
            wfsd_l(i,ib,k1,1) = - p(i)*devr
            wfsd_l(i,ib,k1,2) = - p(i)*devi
         end do
       End do
    end if
  end subroutine SD_direction_noncl
! ===================================================================== 11.0

  subroutine Vnonlocal_Diagonal_part(ispin,ik,iksnl,ibo,vnldi)
    integer, intent(in)                        :: ispin, ik, iksnl,ibo
    real(kind=DP), intent(out), dimension(kg1) :: vnldi

    integer :: it,mdvdb,ib

!!$    integer :: id_sname = -1
!!$    call tstatc0_begin('Vnonlocal_Diagonal_part ', id_sname)

    ib = map_z(ibo)                ! MPI

    vnldi = 0.d0

    do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(mdvdb == SKIP) then
          call Vnonlocal_D_norm_conserve_case
       else if(mdvdb == EXECUT) then
          call Vnonlocal_D_vanderbilt_case
       end if
    end do

!!$    call tstatc0_end(id_sname)
  contains

    subroutine Vnonlocal_D_vanderbilt_case
      integer       :: ia, p1,p2,lmtt1,il1,lmtt2,il2,i
      real(kind=DP) :: ph,fac, dion_eq

      do p1 = 1, ilmt(it)
         lmtt1 = lmtt(p1,it); il1 = ltp(p1,it)
         do p2 = p1, ilmt(it)
            lmtt2 = lmtt(p2,it); il2 = ltp(p2,it)
            if( p1 /= p2) then
               ph = 2.d0*real(zi**(il2-il1))
            else
               ph = 1.d0
            endif
            if(mod(il1+il2,2) == 1) cycle
            dion_eq = dion(p1,p2,it)-eko_l(ib,ik)*q(p1,p2,it)
            fac = 0.d0
            do ia = 1, natm
               if(ityp(ia) /= it) cycle
!!$               fac = fac + ph*iwei(ia) * (dion_eq+vlhxcQ(p1,p2,ia,ispin))
               if(ipaw(it)==0) then
                   fac = fac + ph*iwei(ia) * (dion_eq+vlhxcQ(p1,p2,ia,ispin))
               else
! ======================================== modified by K. Tagami ============= 11.0
!                   fac = ph*iwei(ia) * &
!                        & (dion_paw(p1,p2,ispin,ia)-eko_l(ib,ik)*q(p1,p2,it)+vlhxcQ(p1,p2,ia,ispin))
!
                   fac = fac + ph*iwei(ia) &
        &                   * (dion_paw(p1,p2,ispin,ia)-eko_l(ib,ik)*q(p1,p2,it) &
	&                     +vlhxcQ(p1,p2,ia,ispin) )
! ============================================================================ 11.0
               end if
            end do

            do i = 1, iba(ik)
               vnldi(i) = vnldi(i)+fac*snl(i,lmtt1,iksnl)*snl(i,lmtt2,iksnl)
            end do
         end do
      end do

!!$      do p1 = 1, ilmt(it)
!!$         lmtt1 = lmtt(p1,it); il1 = ltp(p1,it)
!!$         do p2 = p1, ilmt(it)
!!$            lmtt2 = lmtt(p2,it); il2 = ltp(p2,it)
!!$            if( p1 /= p2) then
!!$               ph = 2.d0*real(zi**(il2-il1))
!!$            else
!!$               ph = 1.d0
!!$            endif
!!$            if(mod(il1+il2,2) == 1) cycle
!!$            dion_eq = dion(p1,p2,it)-eko_l(ib,ik)*q(p1,p2,it)
!!$            do ia = 1, natm
!!$               if(ityp(ia) /= it) cycle
!!$!!$               fac = ph*iwei(ia) * &
!!$!!$                    & (dion(p1,p2,it)-eko_l(ib,ik)*q(p1,p2,it)+vlhxcQ(p1,p2,ia,ispin))
!!$               fac = ph*iwei(ia) * (dion_eq+vlhxcQ(p1,p2,ia,ispin))
!!$               do i = 1, iba(ik)
!!$                  vnldi(i) = vnldi(i)+fac*snl(i,lmtt1,iksnl)*snl(i,lmtt2,iksnl)
!!$               end do
!!$            end do
!!$         end do
!!$      end do
    end subroutine Vnonlocal_D_vanderbilt_case

    subroutine Vnonlocal_D_norm_conserve_case
      integer       :: ia, lmt1,lmt2,lmtt1,il1,im1,il2,im2,i
      real(kind=DP) :: ph,fac

      ph = 0.d0
      do ia = 1, natm
         if(ityp(ia) /= it) cycle
         ph = ph + iwei(ia)
      end do
      do lmt1 = 1, ilmt(it)
         lmtt1 = lmtt(lmt1,it); il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
         do lmt2 = lmt1, ilmt(it)
            il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
            if(il1 /= il2 .or. im1 /= im2) cycle
            if(mod(il1+il2,2) == 1) cycle

! ====================================== modified by K. Tagami ============= 11.0
!!            fac = ph * dion(lmt1,lmt2,it)
!
            fac = 0.d0
            Do ia=1, natm
              if(ityp(ia) /= it) cycle
              if(ipaw(it)==0)then
                 fac = fac + iwei(ia) * dion(lmt1,lmt2,it)
              else
                 fac = fac + iwei(ia) * dion_paw(lmt1,lmt2,ispin,ia)
              endif
            End do
! ========================================================================== 11.0

            do i = 1, iba(ik)
               vnldi(i)  = vnldi(i) + fac * snl(i,lmtt1,iksnl)**2
            end do
         end do
      end do
    end subroutine Vnonlocal_D_norm_conserve_case

  end subroutine Vnonlocal_Diagonal_part

! ====================================== added by K. Tagami ================ 11.0
  subroutine Vnonlocal_Diagonal_part_noncl(ispin,ik,iksnl,ibo,vnldi)
    integer, intent(in)                        :: ispin, ik, iksnl,ibo

#ifdef USE_COMPLEX_VNLDI
    complex(kind=CMPLDP), intent(out), dimension(kg1) :: vnldi
#else
    real(kind=DP), intent(out), dimension(kg1) :: vnldi
#endif

    integer :: it,mdvdb,ib

!!$    integer :: id_sname = -1
!!$    call tstatc0_begin('Vnonlocal_Diagonal_part ', id_sname)

    ib = map_z(ibo)                ! MPI

!    vnldi = 0.d0

    do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(mdvdb == SKIP) then
          call Vnonlocal_D_norm_consv_noncl
       else if(mdvdb == EXECUT) then
         call Vnonlocal_D_vanderbilt_noncl
       end if
    end do

!!$    call tstatc0_end(id_sname)
  contains

    subroutine Vnonlocal_D_vanderbilt_noncl
      integer       :: ia, p1,p2,lmtt1,il1,lmtt2,il2,i
      integer :: im1, im2
      complex(kind=CMPLDP) :: ph, fac

      do p1 = 1, ilmt(it)
         lmtt1 = lmtt(p1,it); il1 = ltp(p1,it);  im1 = mtp(p1,it)

         do p2 = 1, ilmt(it)
            lmtt2 = lmtt(p2,it); il2 = ltp(p2,it);  im2 = mtp(p2,it)


! ---------------------------------------------------- 11.0S
#ifdef SKIP_TEST
            if ( il1 /= il2 ) cycle
            if ( SpinOrbit_mode == Neglected .and. sw_hubbard == OFF ) then
               if ( im1 /= im2 ) cycle
            endif
#endif
! ---------------------------------------------------- 11.0S

            if( p1 /= p2 ) then
               ph = zi**(il2-il1)
           else
               ph = 1.d0
            endif

!            if(mod(il1+il2,2) == 1) cycle

            fac = 0.d0

            do ia = 1, natm
               if(ityp(ia) /= it) cycle
               fac = fac + ph*iwei(ia)* ( dion_scr_noncl( p1,p2,ispin,ia ) &
        &                                -eko_l(ib,ik) *q_noncl( p1, p2, ispin, it ) )
            end do

            do i = 1, iba(ik)
               vnldi(i) = vnldi(i)+fac*snl(i,lmtt1,iksnl)*snl(i,lmtt2,iksnl)
            end do
         end do
      end do

    end subroutine Vnonlocal_D_vanderbilt_noncl

    subroutine Vnonlocal_D_norm_consv_noncl
      integer       :: ia, p1,p2,lmtt1,il1,lmtt2,il2,i
      integer :: im1, im2
      complex(kind=CMPLDP) :: ph, fac

      do p1 = 1, ilmt(it)
         lmtt1 = lmtt(p1,it); il1 = ltp(p1,it);  im1 = mtp(p1,it)

         do p2 = 1, ilmt(it)
            lmtt2 = lmtt(p2,it); il2 = ltp(p2,it);  im2 = mtp(p2,it)

! ---------------------------------------------------- 11.0S
            if ( il1 /= il2 ) cycle
            if ( SpinOrbit_mode == Neglected .and. sw_hubbard == OFF ) then
               if ( im1 /= im2 ) cycle
            endif
! ---------------------------------------------------- 11.0S

            if( p1 /= p2 ) then
               ph = zi**(il2-il1)
            else
               ph = 1.d0
            endif

!            if(mod(il1+il2,2) == 1) cycle

            fac = 0.d0

            do ia = 1, natm
               if(ityp(ia) /= it) cycle
               fac = fac + ph*iwei(ia)* dion_scr_noncl( p1,p2,ispin,ia )

            end do

            do i = 1, iba(ik)
               vnldi(i) = vnldi(i)+fac*snl(i,lmtt1,iksnl)*snl(i,lmtt2,iksnl)
            end do
         end do
      end do

    end subroutine Vnonlocal_D_norm_consv_noncl

  end subroutine Vnonlocal_Diagonal_part_noncl
! ==================================================================== 11.0


  subroutine m_ESsd_reset_dzajn2
    dzajn2 = 0.d0
  end subroutine m_ESsd_reset_dzajn2

  subroutine m_ESsd_copy_zaj_to_zaj_old
    integer :: ik, ir
    do ir = 1, kimg
       do ik = 1, kv3, af+1
          if(map_k(ik) /= myrank_k) cycle                       ! MPI
          zaj_old(:,:,ik,ir) = zaj_l(:,:,ik,ir)
       end do
    end do
  end subroutine m_ESsd_copy_zaj_to_zaj_old

  subroutine m_ESsd_copy_zaj_old_to_zaj(ik,ibo)
    integer, intent(in)                           :: ik, ibo

    integer    :: ir, ib

    ib = map_z(ibo)
     do ir = 1, kimg
       zaj_l(:,ib,ik,ir) = zaj_old(:,ib,ik,ir)
    end do
#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) status_saved_phifftr(ib,ik) = OLD
#endif
  end subroutine m_ESsd_copy_zaj_old_to_zaj

  subroutine m_ESsd_copy_phi_to_zaj_old(ik,ibo,phi)
    integer, intent(in)                           :: ik, ibo
    real(kind=DP),intent(in),dimension(kg1,kimg)  :: phi

    integer    :: ir, ib

    ib = map_z(ibo)
     do ir = 1, kimg
       zaj_old(:,ib,ik,ir) = phi(:,ir)
    end do
  end subroutine m_ESsd_copy_phi_to_zaj_old
!!$!BRANCH_P 3D_Parallel
!!$  subroutine m_ESsd_copy_phi_to_zaj_old_3D(ik,ibo,phi)
!!$    integer, intent(in)                           :: ik, ibo
!!$    real(kind=DP),intent(in),dimension(maxval(np_g1k),kimg)  :: phi
!!$
!!$    integer    :: ir, ib
!!$
!!$     do ir = 1, kimg
!!$       zaj_old(:,ib,ik,ir) = phi(:,ir)
!!$    end do
!!$  end subroutine m_ESsd_copy_phi_to_zaj_old_3D
!!$!BRANCH_P_END 3D_Parallel

  subroutine m_ESsd_diff_WFs_to_zaj_old(dt)
    real(kind=DP), intent(in) :: dt
    integer                   :: ik, ir
    real(kind=DP), parameter  :: Delta = 1.d-10
    real(kind=DP)             :: rdt

    if(dt < Delta ) then
       if(ipri >=1 ) write(nfout,'(" !! dt is smaller than ",d20.8)') Delta
       rdt = 0.d0
    else
       rdt = 1.d0/dt
    end if

    do ir = 1, kimg
       do ik = 1, kv3, af+1
          if(map_k(ik) /= myrank_k) cycle                ! MPI
          zaj_old(:,:,ik,ir) = rdt * (zaj_l(:,:,ik,ir) - zaj_old(:,:,ik,ir))
       end do
    end do
  end subroutine m_ESsd_diff_WFs_to_zaj_old

  subroutine WF_conjugate_gradient(ik,ibo,dtim)
! \Psi_{k\mu}^{m+1} = \Psi_{k\mu}^{m} + \Delta t_{opt}^{m} d^{m}
!   d^{m} = g^{m} + \beta_CG d^{m-1}
!   g^{m} = - \Phi_{k\mu}^{m}
!   d^{m-1} = (\Psi_{k\mu}^{m} - \Psi_{k\mu}^{m-1})\frac{1}{\Delta t_{opt}^{m-1}}
!
    integer      , intent(in)                  :: ik,ibo
    real(kind=DP), intent(in)                  :: dtim

    integer       :: i,ib
    real(kind=DP) :: evr,evi

    ib = map_z(ibo)                    ! MPI
#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) status_saved_phifftr(ib,ik) = OLD
#endif
    if(kimg == 1) then
       do i = 1, iba(ik)
          evr   = zaj_l(i,ib,ik,1)
          zaj_l(i,ib,ik,1) = evr + dtim*wfsd_l(i,ib,ik,1)
       end do
    else if(kimg == 2) then
       do i = 1, iba(ik)
          evr   = zaj_l(i,ib,ik,1);  evi   = zaj_l(i,ib,ik,kimg)
          zaj_l(i,ib,ik,1) = evr + dtim*wfsd_l(i,ib,ik,1)
          zaj_l(i,ib,ik,2) = evi + dtim*wfsd_l(i,ib,ik,2)
       end do
    end if
  end subroutine WF_conjugate_gradient

  subroutine WF_conjugate_gradient0 &
       &(precon,ik,ibo,dtim,ekin,VlocalW,p)
! \Psi_{k\mu}^{m+1} = \Psi_{k\mu}^{m} + \Delta t_{opt}^{m} d^{m}
!   d^{m} = g^{m} + \beta_CG d^{m-1}
!   g^{m} = - \Phi_{k\mu}^{m}
!   d^{m-1} = (\Psi_{k\mu}^{m} - \Psi_{k\mu}^{m-1})\frac{1}{\Delta t_{opt}^{m-1}}
!
    integer      , intent(in)                  :: precon,ik,ibo
    real(kind=DP), intent(in)                  :: dtim
    real(kind=DP), intent(in), dimension(kg1)  :: ekin
    real(kind=DP), intent(in), dimension(nfft) :: VlocalW
    real(kind=DP)            , dimension(kg1)  :: p

    integer       :: i, i1, ib
    real(kind=DP) :: evr,devr,denom, evi,e1, devi

    ib = map_z(ibo)
    denom = 1.d0/product(fft_box_size_WF(1:3,1))
    call m_ES_decide_precon_factor(precon,ik,ibo,ekin,p)  ! ->p(1:iba(ik))

#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) status_saved_phifftr(ib,ik) = OLD
#endif
    if(kimg == 1) then
       do i = 1, iba(ik)
          i1    = igf(nbase(i,ik))
          evr   = zaj_l(i,ib,ik,1)
          devr  = (ekin(i)-eko_l(ib,ik))*evr&
               & + VlocalW(i1)*denom + vnlph_l(i,ib,1)
          zaj_l(i,ib,ik,1) = evr + dtim*(-p(i)*devr+betacg*zaj_old(i,ib,ik,1))
          zaj_old(i,ib,ik,1) = evr
       end do
    else if(kimg == 2) then
       do i = 1, iba(ik)
          i1    = igf(nbase(i,ik))
          evr   = zaj_l(i,ib,ik,1)
          evi   = zaj_l(i,ib,ik,kimg)
          e1    = ekin(i) - eko_l(ib,ik)
          devr  = e1*evr+VlocalW(2*i1-1)*denom+vnlph_l(i,ib,1)
          devi  = e1*evi+VlocalW(2*i1  )*denom+vnlph_l(i,ib,2)
          zaj_l(i,ib,ik,1) = evr +dtim*(-p(i)*devr+betacg*zaj_old(i,ib,ik,1))
          zaj_l(i,ib,ik,2) = evi +dtim*(-p(i)*devi+betacg*zaj_old(i,ib,ik,2))
          zaj_old(i,ib,ik,1) = evr; zaj_old(i,ib,ik,2) = evi
       end do
    end if
  end subroutine WF_conjugate_gradient0

  subroutine make_CG_direction(ik)
    integer, intent(in) :: ik

    real(kind=DP), parameter             :: Delta = 1.d-20
    integer :: i, ib
    integer :: id_sname = -1
    call tstatc0_begin('make_CG_direction', id_sname)

    if(betacg > Delta) then
       if(kimg == 1) then
          do ib = 1, np_e
             do i = 1, iba(ik)
                wfsd_l(i,ib,ik,1) = wfsd_l(i,ib,ik,1) + betacg*wfsd_old(i,ib,ik,1)
             end do
          end do
       else
          do ib = 1, np_e
             do i = 1, iba(ik)
                wfsd_l(i,ib,ik,1) = wfsd_l(i,ib,ik,1) + betacg*wfsd_old(i,ib,ik,1)
                wfsd_l(i,ib,ik,2) = wfsd_l(i,ib,ik,2) + betacg*wfsd_old(i,ib,ik,2)
             end do
          end do
       end if
    end if
    call tstatc0_end(id_sname)
  end subroutine make_CG_direction

  subroutine orthogonalize_SD_drctns(ik,to)
    integer, intent(in) :: ik,to

    integer :: id_sname = -1
    call tstatc0_begin('orthogonalize_SD_drctns ', id_sname)

!!$    if(modnrm == EXECUT) call m_ES_betar_dot_Psi_4_each_k(wfsd_l,ista_k,iend_k,ik,bsdr_l,bsdi_l)
    if(modnrm == EXECUT) call m_ES_betar_dot_Psi_4_each_k(nfout,wfsd_l,ik,ik,ik,bsdr_l,bsdi_l)
    !                                                          ->bsd(ri)_l
    call m_ES_orthogonalize_SD_to_WFs(ik,to,wfsd_l,bsdr_l,bsdi_l) ! ->(wfsd_l,bsd(ri)_l)

!!$    call orthogonalize_SD(ik,wfsd_l,bsdr_l,bsdi_l)           !-(m_E.S.) ->(wfsd_l)
    call tstatc0_end(id_sname)
  end subroutine orthogonalize_SD_drctns

! =============================== added by K. Tagami =================== 11.0
  subroutine orthogonalize_SD_drctns_noncl( ik,to )
    integer, intent(in) :: ik,to

    integer :: is
    integer :: id_sname = -1

    call tstatc0_begin('orthogonalize_SD_drctns_noncl ', id_sname)

    if (modnrm == EXECUT) then
       Do is=1, ndim_spinor
          call m_ES_betar_dot_Psi_4_each_k(nfout, wfsd_l, ik, ik+ndim_spinor-1, &
	&                                   ik +is -1, bsdr_l, bsdi_l )
                                               !           ->bsd(ri)_l
       End do
    endif

    call m_ES_orthogonl_SD_to_WFs_noncl( ik, to, ik, ik+ndim_spinor-1, &
	&                                wfsd_l, bsdr_l, bsdi_l)
                                             ! ->(wfsd_l,bsd(ri)_l)

    call tstatc0_end(id_sname)

  end subroutine orthogonalize_SD_drctns_noncl
! ========================================================================= 11.0

  subroutine Precon_dot_SD(ik,ekin,p)
    integer,      intent(in)                 :: ik
    real(kind=DP),intent(in),dimension(kg1)  :: ekin
    real(kind=DP),           dimension(kg1)  :: p

    integer :: ib, ri, i

    do ri = 1, kimg
       do ib = ista_e, iend_e, istep_e
          do i = 1, iba(ik)
             wfsd_np(i,map_z(ib),ik,ri) = wfsd_l(i,map_z(ib),ik,ri)
          end do
       end do
    end do

    do ib = ista_e, iend_e, istep_e   ! MPI
       call decide_precon_factor_wfsd(ik,ib,ekin,p)
       if(kimg == 1) then
          do i = 1, iba(ik)
             wfsd_l(i,map_z(ib),ik,1) = p(i)*wfsd_l(i,map_z(ib),ik,1)
          end do
       else if(kimg == 2) then
          do i = 1, iba(ik)
             wfsd_l(i,map_z(ib),ik,1) = p(i)*wfsd_l(i,map_z(ib),ik,1)
             wfsd_l(i,map_z(ib),ik,2) = p(i)*wfsd_l(i,map_z(ib),ik,2)
          end do
       end if
    end do

  end subroutine Precon_dot_SD

! ======================================== added by K. Tagami =============== 11.0
  subroutine Precon_dot_SD_noncl( ik,ekin,p )
    integer,      intent(in)                 :: ik
    real(kind=DP),intent(in),dimension(kg1)  :: ekin
    real(kind=DP),           dimension(kg1)  :: p

    integer :: ib, ri, i
    integer :: is, k1

    do ri = 1, kimg
       Do is=1, ndim_spinor
         k1 = ik + is -1
         do ib = ista_e, iend_e, istep_e
            do i = 1, iba(ik)
               wfsd_np(i,map_z(ib),k1,ri) = wfsd_l(i,map_z(ib),k1,ri)
            end do
         end do
       End do
    end do

    do ib = ista_e, iend_e, istep_e   ! MPI
       call decide_precon_factor_wfsd( ik, ib, ekin, p )

       if (kimg == 1) then
          Do is=1, ndim_spinor
            k1 = ik + is -1
            do i = 1, iba(ik)
               wfsd_l(i,map_z(ib),k1,1) = p(i)*wfsd_l(i,map_z(ib),k1,1)
            end do
          End do
       else if(kimg == 2) then
          Do is=1, ndim_spinor
            k1 = ik + is -1
            do i = 1, iba(ik)
               wfsd_l(i,map_z(ib),k1,1) = p(i)*wfsd_l(i,map_z(ib),k1,1)
               wfsd_l(i,map_z(ib),k1,2) = p(i)*wfsd_l(i,map_z(ib),k1,2)
            end do
          End do
       end if
    end do

  end subroutine Precon_dot_SD_noncl
! ===================================================================== 11.0

  subroutine decide_precon_factor_wfsd(ik,ib,ekin,p)
    integer, intent(in)                         :: ik,ib
    real(kind=DP), intent(in),  dimension(kg1)  :: ekin
    real(kind=DP), intent(out), dimension(kg1)  :: p

    integer       :: i
    real(kind=DP) :: ektot, x, x1, x2, d_ektot

! ===================================== modified by K. Tagami ========== 11.0
!    call kinetic_energy_wfsd(ik,ib,ekin,ektot)   ! -(m_ES_WF_by_SDorCG)
!
    if ( noncol ) then
       call kinetic_energy_wfsd_noncl(ik,ib,ekin,ektot)
    else
       call kinetic_energy_wfsd(ik,ib,ekin,ektot)   ! -(m_ES_WF_by_SDorCG)
    endif
! ====================================================================== 11.0

    d_ektot = 1.d0/ektot
    p = 0.d0
    do i = 1, iba(ik)
       x = ekin(i)*d_ektot
       x1 = 27 + ( 18 + (12 + 8*x) *x) *x
       x2 = 16*(x*x)*(x*x)
       p(i)  = x1/(x1 + x2 )
    end do
  end subroutine decide_precon_factor_wfsd

  subroutine kinetic_energy_wfsd(ik,ib,dekin,ektot)
    integer, intent(in) :: ik, ib
    real(kind=DP), intent(in), dimension(kg1)  :: dekin
    real(kind=DP), intent(out)                 :: ektot
    integer  :: i, ri
    ektot = 0.d0

    do ri = 1, kimg
       do i = 1, iba(ik)
          ektot = ektot + dekin(i)*wfsd_l(i,map_z(ib),ik,ri)**2   ! MPI
       end do
    end do

    if(k_symmetry(ik) == GAMMA) ektot = ektot*2.d0

  end subroutine kinetic_energy_wfsd

! ====================================== added by K. Tagami ================ 11.0
  subroutine kinetic_energy_wfsd_noncl( ik,ib,dekin,ektot )
    integer, intent(in) :: ik, ib
    real(kind=DP), intent(in), dimension(kg1)  :: dekin
    real(kind=DP), intent(out)                 :: ektot
    integer  :: i, ri
    integer :: is, k1

    ektot = 0.d0

    do ri = 1, kimg
       Do is=1, ndim_spinor
         k1 = ik + is -1
         do i = 1, iba(ik)
            ektot = ektot + dekin(i)*wfsd_l(i,map_z(ib),k1,ri)**2   ! MPI
         end do
       End do
    end do

    if(k_symmetry(ik) == GAMMA) ektot = ektot*2.d0

  end subroutine kinetic_energy_wfsd_noncl
! ========================================================================= 11.0


  subroutine cp_wfsd_to_wfsd_old(ik)
    integer, intent(in) :: ik
    wfsd_old(:,:,ik,:) = wfsd_l(:,:,ik,:)
    bsdr_old(:,:,ik) = bsdr_l(:,:,ik)
    bsdi_old(:,:,ik) = bsdi_l(:,:,ik)
  end subroutine cp_wfsd_to_wfsd_old

  subroutine m_ESsd_dealloc
    if(allocated(zaj_old)) deallocate(zaj_old)
    if(allocated(dzajn2)) deallocate(dzajn2)


  end subroutine m_ESsd_dealloc

! =============
!   meta-gga
!
  subroutine m_ES_contrib_kindens_to_vnlph( is, ik, afft )
    use m_Crystal_Structure,  only : rltv
    use m_PlaneWaveBasisSet,   only : ngabc, kgp
    use m_Electronic_Structure,   only : vtau_l, m_ES_WF_in_Rspace_kt, &
         &                               m_ES_add_it_to_vnlph, &
         &                               m_ES_Vlocal_in_Rspace

    use m_FFT,                 only : m_FFT_Vlocal_W

    integer, intent(in) :: is, ik
    real(kind=DP), intent(in) :: afft(nfft)

    integer :: ib, in, ri, i1, ig
    real(kind=DP) :: denom
    real(kind=DP), allocatable :: G_vec(:,:), k_plus_G(:,:), vlength(:)
    real(kind=DP), allocatable :: bfft(:)
    real(kind=DP), allocatable :: zwk1(:,:,:,:), zwk2(:,:)

    allocate( vlength(kg1) )
    allocate( k_plus_G(kg1,3) );
    allocate( bfft(nfft) )

    call k_plus_G_vectors( ik, kgp, kg1, kv3, iba, nbase, vkxyz, &
         &                 ngabc, rltv, &
         &                 k_plus_G(1,1), k_plus_G(1,2), k_plus_G(1,3), vlength )

    denom = 1.d0 /product(fft_box_size_WF(1:3,1))

    do ib = 1, np_e
       allocate( zwk1(kg1,1,ik:ik,kimg ) );  zwk1 = 0.0d0
       allocate( zwk2(kg1,kimg) );  zwk2 = 0.0d0
#if 0
       Do in=1, 3
          Do ig=1, iba(ik)
             zwk1(ig,1,ik,2) = zaj_l(ig,ib,ik,1) *k_plus_G(ig,in)
             zwk1(ig,1,ik,1) =-zaj_l(ig,ib,ik,2) *k_plus_G(ig,in)
          End do
          call m_ES_WF_in_Rspace_kt( ik ,ik, ik, zwk1, bfft )
          call m_FFT_Vlocal_W( afft, bfft )  ! --> bfft  in R space
          call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)  ! G space

          do ig=1,iba(ik)
             i1 = kimg*igf(nbase(ig,ik)) -1
             zwk2(ig,2) = zwk2(ig,2) +bfft(i1) *k_plus_G(ig,in) *denom
             zwk2(ig,1) = zwk2(ig,1) -bfft(i1+1) *k_plus_G(ig,in) *denom
          enddo
       end Do
       zwk2 = -zwk2 /2.0d0       ! -> vnlph_l   zwk2( kg1, kimg )
#else
       Do in=1, 3
          do ri=1,kimg
             Do ig=1, iba(ik)
                zwk1(ig,1,ik,ri) = zaj_l(ig,ib,ik,ri) *k_plus_G(ig,in)
             End Do
          End do

          call m_ES_WF_in_Rspace_kt( ik ,ik, ik, zwk1, bfft )
          call m_FFT_Vlocal_W( afft, bfft )  ! --> bfft  in R space
          call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)  ! G space

          do ri=1,kimg
             do ig=1,iba(ik)
                i1 = kimg*igf(nbase(ig,ik)) + (ri - kimg)
                zwk2(ig,ri) = zwk2(ig,ri) +bfft(i1) *k_plus_G(ig,in) *denom
             enddo
          enddo
       End Do

       zwk2 = zwk2 /2.0d0       ! -> vnlph_l   zwk2( kg1, kimg )
!       zwk2 = zwk2 /4.0d0       ! -> vnlph_l   zwk2( kg1, kimg )

#endif
       call m_ES_add_it_to_vnlph( ik, ib, zwk2 )
       deallocate( zwk1 );   deallocate( zwk2 )
    end do
    deallocate( k_plus_G );  deallocate( vlength )
    deallocate( bfft )

  end subroutine m_ES_contrib_kindens_to_vnlph

  subroutine m_ES_kindens_to_vnlph_ib( is, ik, ib, afft )
    use m_Crystal_Structure,  only : rltv
    use m_PlaneWaveBasisSet,   only : ngabc, kgp
    use m_FFT,                 only : m_FFT_Vlocal_W
    use m_Electronic_Structure,   only : vtau_l, m_ES_WF_in_Rspace_kt, &
         &                               m_ES_add_it_to_vnlph, &
         &                               m_ES_Vlocal_in_Rspace

    integer, intent(in) :: is, ik, ib
    real(kind=DP), intent(in) :: afft(nfft)

    integer :: in, ri, i1, ig
    real(kind=DP) :: denom
    real(kind=DP), allocatable :: G_vec(:,:), k_plus_G(:,:), vlength(:)
    real(kind=DP), allocatable :: bfft(:)
    real(kind=DP), allocatable :: zwk1(:,:,:,:), zwk2(:,:)

    allocate( vlength(kg1) );   allocate( k_plus_G(kg1,3) );
    allocate( bfft(nfft) )

    call k_plus_G_vectors( ik, kgp, kg1, kv3, iba, nbase, vkxyz, &
         &                 ngabc, rltv, &
         &                 k_plus_G(1,1), k_plus_G(1,2), k_plus_G(1,3), vlength )

    denom = 1.d0 /product(fft_box_size_WF(1:3,1))

    allocate( zwk1(kg1,1,ik:ik,kimg ) );  zwk1 = 0.0d0
    allocate( zwk2(kg1,kimg) );  zwk2 = 0.0d0
    Do in=1, 3
       do ri=1,kimg
          Do ig=1, iba(ik)
             zwk1(ig,1,ik,ri) = zaj_l(ig,ib,ik,ri) *k_plus_G(ig,in)
          End Do
       End do
       call m_ES_WF_in_Rspace_kt( ik ,ik, ik, zwk1, bfft )
       call m_FFT_Vlocal_W( afft, bfft )  ! --> bfft  in R space
       call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)  ! G space

       do ri=1,kimg
          do ig=1,iba(ik)
             i1 = kimg*igf(nbase(ig,ik)) + (ri - kimg)
             zwk2(ig,ri) = zwk2(ig,ri) +bfft(i1) *k_plus_G(ig,in) *denom
          enddo
       enddo
    End Do
    zwk2 = zwk2 /2.0d0       ! -> vnlph_l   zwk2( kg1, kimg )

    if(kimg==1) then
       do ig=1,iba(ik)
          vnlph_l(ig,ib,1) = vnlph_l(ig,ib,1) + zwk2(ig,1)
       end do
    else
       do ig=1,iba(ik)
          vnlph_l(ig,ib,1) = vnlph_l(ig,ib,1) + zwk2(ig,1)
          vnlph_l(ig,ib,2) = vnlph_l(ig,ib,2) + zwk2(ig,2)
       end do
    end if

    deallocate( zwk1 );   deallocate( zwk2 )
    deallocate( k_plus_G );  deallocate( vlength )
    deallocate( bfft )

  end subroutine m_ES_kindens_to_vnlph_ib

  subroutine m_ES_kindens_to_vnlph_ib2( is, ik, ib, afft, vtau_phl )
    use m_Crystal_Structure,  only : rltv
    use m_PlaneWaveBasisSet,   only : ngabc, kgp
    use m_FFT,                 only : m_FFT_Vlocal_W
    use m_Electronic_Structure,   only : vtau_l, m_ES_WF_in_Rspace_kt, &
         &                               m_ES_add_it_to_vnlph, &
         &                               m_ES_Vlocal_in_Rspace

    integer, intent(in) :: is, ik, ib
    real(kind=DP), intent(in) :: afft(nfft)
    real(kind=DP), intent(out) :: vtau_phl(kg1,kimg)

    integer :: in, ri, i1, ig
    real(kind=DP) :: denom
    real(kind=DP), allocatable :: G_vec(:,:), k_plus_G(:,:), vlength(:)
    real(kind=DP), allocatable :: bfft(:)
    real(kind=DP), allocatable :: zwk1(:,:,:,:), zwk2(:,:)

    allocate( vlength(kg1) );   allocate( k_plus_G(kg1,3) );
    allocate( bfft(nfft) )

    call k_plus_G_vectors( ik, kgp, kg1, kv3, iba, nbase, vkxyz, &
         &                 ngabc, rltv, &
         &                 k_plus_G(1,1), k_plus_G(1,2), k_plus_G(1,3), vlength )

    denom = 1.d0 /product(fft_box_size_WF(1:3,1))

    allocate( zwk1(kg1,1,ik:ik,kimg ) );  zwk1 = 0.0d0
    allocate( zwk2(kg1,kimg) );  zwk2 = 0.0d0
    Do in=1, 3
       do ri=1,kimg
          Do ig=1, iba(ik)
             zwk1(ig,1,ik,ri) = zaj_l(ig,ib,ik,ri) *k_plus_G(ig,in)
          End Do
       End do
       call m_ES_WF_in_Rspace_kt( ik ,ik, ik, zwk1, bfft )
       call m_FFT_Vlocal_W( afft, bfft )  ! --> bfft  in R space
       call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)  ! G space

       do ri=1,kimg
          do ig=1,iba(ik)
             i1 = kimg*igf(nbase(ig,ik)) + (ri - kimg)
             zwk2(ig,ri) = zwk2(ig,ri) +bfft(i1) *k_plus_G(ig,in) *denom
          enddo
       enddo
    End Do
    zwk2 = zwk2 /2.0d0       ! -> vnlph_l   zwk2( kg1, kimg )

    if(kimg==1) then
       do ig=1,iba(ik)
          vtau_phl(ig,1) = zwk2(ig,1)
       end do
    else
       do ig=1,iba(ik)
          vtau_phl(ig,1) = zwk2(ig,1)
          vtau_phl(ig,2) = zwk2(ig,2)
       end do
    end if

    deallocate( zwk1 );   deallocate( zwk2 )
    deallocate( k_plus_G );  deallocate( vlength );  deallocate( bfft )

  end subroutine m_ES_kindens_to_vnlph_ib2


end module m_ES_WF_by_SDorCG
