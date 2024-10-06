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
#ifdef MPI_FFTW
  use, intrinsic :: iso_c_binding
#endif
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
 &                                , vlhxc_l &
#ifdef SAVE_FFT_TIMES
 &                                , status_saved_phifftr &
#endif
 &                                , m_ES_sum_of_LocalPart  &
 &                                , m_ES_sum_of_LocalPart2 &
 &                                , m_ES_sum_of_LocalPart3 &
 &                                , m_ES_WF_2D

  use m_ES_ortho            ,only : m_ES_MGS_4_each_k
  use m_ES_nonlocal         ,only : sc,ss,qc                                          &
 &                                , m_ES_dealloc_afft_scss_etc
  use m_FFT,                 only : nfft, fft_box_size_WF
  use m_Ionic_System,        only : ntyp, natm, ityp, iwei
  use m_Kpoints,             only : kv3,vkxyz,k_symmetry
  use m_NonLocal_Potential,  only : snl
  use m_Files,               only : nfout
  use m_Parallelization,     only : MPI_CommGroup,is_kngp,ie_kngp,npes,mype &
       &, nrank_e, myrank_e, map_e, ista_e, iend_e, istep_e, map_z, np_e &
       &, nrank_k, myrank_k, map_k, ista_k, iend_k, mpi_k_world, nis_k, nel_k, ierr, ista_fs, iend_fs
  use m_PlaneWaveBasisSet,   only : kg1,iba, igf, nbase, m_pwBS_kinetic_energies
  use m_PseudoPotential,     only : ilmt, lmtt, ltp, mtp, q, dion, modnrm, nlmta &
       &                           ,m_PP_include_vanderbilt_pot,ipaw,dion_paw
  use m_Timing,              only : tstatc0_begin, tstatc0_end
  use m_Electronic_Structure,only : nrvf_ordr,neordr,neordr_old &
       &                          , fsr_l_2d,fsi_l_2d,fsr_l,fsi_l  &
       &                          , vlhxc_l, vtau_l                                &
       &                          , m_ES_Vlocal_in_Rspace_3D                          &
       &                          , m_ES_WF_in_Rspace_3D                              &
       &                          , m_ES_decide_precon_factor_3D                    &
       &                          , m_ES_eigen_values_for_each_k_3D                   &
       &                          , m_ES_sort_eigen_values_3D                         &
       &                          , m_ES_wd_eko_3D                                    &
       &                          , m_ES_wd_zaj_small_portion_3D                      &
       &                          , nblocksize_mgs_default
!!$  use m_ES_ortho            ,only : m_ES_MGS_4_each_k_3D                              &
!!$ &                                , m_ES_orthogonalize_SD_to_WFs_3D
!!$  use m_ES_ortho            ,only : m_ES_orthogonalize_SD_to_WFs_3D
  use m_ES_ortho            ,only : m_ES_orthogonal_phi_to_WFs
  use m_ES_nonlocal         ,only : sc_l                                              &
 &                                , m_ES_alloc_afft_scss_etc_3D                       &
 &                                , m_ES_Vnonlocal_W_3D                               &
 &                                , m_ES_betar_dot_WFs_4_each_k_3D                    &
 &                                , m_ES_betar_dot_Psi_4_each_k_3D
  use m_FFT                 ,only : m_FFT_Direct_3D, m_FFT_Vlocal_W_3D                &
#ifdef FFT_3D_DIVISION
       &                          , m_FFT_Direct_3DIV_3D, m_FFT_Vlocal_W_3DIV_3D      &
#endif
       &                          , m_FFT_Direct_XYZ_3D
  use m_Parallelization,     only : nel_fft_z, nel_fft_y, nel_fft_x &
       &                          , fft_X_x_nel, fft_X_y_nel, fft_X_z_nel &
       &                          , mp_fft_x, xyz_fft_x &
       &                          , myrank_g, myrank_chg, nrank_g, nrank_chg   &
       &                          , mp_e &
       &                          , nis_e, nie_e, nel_e &
       &                          , ista_g1k, iend_g1k, np_g1k, mp_g1k, nel_g1k  &
       &                          , mpi_ke_world, mpi_chg_world,  mpi_kg_world, ista_g1k, iend_g1k &
       &                          , nis_g1k, nie_g1k, neg_g, neg_g_all &
       &                          , neg_gg, neg_gg_all &
       &                          , ista_atm, iend_atm, np_fs
  use m_PseudoPotential,     only : nlmtt
  use m_Control_Parameters,  only : nblocksize_fftw, nblocksize_fftw_is_given, sw_fft_xzy,sw_serial_fft
!  use z_interface_3D
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
!  use m_ES_ortho,              only : m_ES_orthogonl_SD_to_WFs_noncl, &
!       &                              m_ES_MGS_4_each_k_noncl

!  use m_Electronic_Structure,   only : m_ES_sum_of_LocalPart2_noncl, &
!       &                               m_ES_Vlocal_in_Rspace_noncl, &
!       &                               m_ES_sort_eigen_vals_noncl, &
!       &                               m_ES_eigen_vals_each_k_noncl, &
!       &                               m_ES_sum_of_LocalPart_noncl, &
!       &                               m_ES_sum_of_LocalPart3_noncl

!  use m_ES_nonlocal,            only : m_ES_alloc_scss_etc, &
!       &                               m_ES_dealloc_scss_etc
!!  use m_FFT,                    only : m_FFT_Vlocal_W_noncl
!! for debug
  use m_Electronic_Structure,     only : fsr_l, fsi_l
! ===================================================================== 11.0
! === 0 divide occurs if abs(wdi) < epsilon. by tkato 2012/12/18 ===============
  use m_Const_Parameters, only : SmallestPositiveNumber
! ==============================================================================

  use m_Control_Parameters,  only : use_metagga, vtau_exists

  use m_IterationNumbers, only : iteration_electronic

#ifdef MPI_FFTW
  use m_Electronic_Structure,only : m_ES_WF_in_Rspace_mpifftw, m_ES_Vlocal_in_Rspace_mpifftw, m_ES_Vlocal_in_Rspace_mpifftw3d
  use m_FFT,                 only : m_FFT_Vlocal_W_mpifftw, m_FFT_Direct_MPI_FFTW, afft_mpifftw, bfft_mpifftw &
                                  , m_FFT_Vlocal_W_mpifftw3d, afft_mpifftw_vlocal
#endif
  use mpi

  implicit none
! === This should be modified after zaj_old is 3D-decomposed!!! by T.Kato ======
  real(kind=DP),public,allocatable,dimension(:,:,:,:):: zaj_old    !d(kg1,np_e,ista_k:iend_k,kimg)
!  real(kind=DP),public,allocatable,dimension(:,:,:,:):: zaj_old !d(kg1,np_e,ista_k:iend_k,kimg)
! ==============================================================================
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

#ifdef MPI_FFTW
  integer,  allocatable, dimension(:) :: mmx,mmy,mmz,iiadd,wfdist
  integer :: maxwfdist
  logical :: mpifftw_map_built = .false.
#endif

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

#ifdef MPI_FFTW
  include 'fftw3-mpi.f03'
#endif
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
! === This should be modified after zaj_old is 3D-decomposed!!! by T.Kato ======
    allocate(zaj_old(maxval(np_g1k),np_e,ista_k:iend_k,kimg))
! ==============================================================================
  end subroutine m_ESsd_alloc_zaj_old

  subroutine m_ESsd_alloc_wfsd
    if(.not.allocated(wfsd_old)) then
! ==============================================================================
       allocate(wfsd_old(maxval(np_g1k),np_e,ista_k:iend_k,kimg))
       wfsd_old = 0.d0
! ==============================================================================
    end if
    if(.not.allocated(bsdr_old)) then
       allocate(bsdr_old(np_e,np_fs,ista_k:iend_k));bsdr_old=0.0d0
    endif
    if(.not.allocated(bsdi_old)) then
       allocate(bsdi_old(np_e,np_fs,ista_k:iend_k));bsdi_old=0.0d0
    endif
  end subroutine m_ESsd_alloc_wfsd

  subroutine m_ESsd_dealloc_wfsd
    if(allocated(wfsd_old)) deallocate(wfsd_old)
    if(allocated(bsdr_old)) deallocate(bsdr_old)
    if(allocated(bsdi_old)) deallocate(bsdi_old)
  end subroutine m_ESsd_dealloc_wfsd

! -------------------- Added by T. Yamasaki, 28 June 2008 --->>
  subroutine m_ESsd_alloc_wfred()
    allocate(wfred_l(maxval(np_g1k),np_e,ista_k:iend_k,kimg)); wfred_l = 0.d0
    wfred_is_allocated = .true.
  end subroutine m_ESsd_alloc_wfred

  subroutine m_ESsd_dealloc_wfred()
    deallocate(wfred_l)
    wfred_is_allocated = .false.
  end subroutine m_ESsd_dealloc_wfred

  subroutine cp_wfred_3D(ik)
    integer, intent(in) :: ik
    integer :: ir, ib, ig
    do ir = 1, kimg
       do ib = 1, np_e
          do ig = 1, np_g1k(ik)
             wfred_l(ig,ib,ik,ir) = zaj_l(ig,ib,ik,ir) - zaj_old(ig,ib,ik,ir)
          end do
       end do
    end do
  end subroutine cp_wfred_3D
! ----------------------------------<<

  subroutine alloc_wfsd_bsdri_3D(ik)
    integer, intent(in) :: ik
    allocate(wfsd_l(maxval(np_g1k),np_e,ik:ik,kimg)); wfsd_l = 0.d0
    allocate(bsdr_l(np_e,np_fs,ik:ik))
    allocate(bsdi_l(np_e,np_fs,ik:ik)); bsdi_l = 0.0d0
!!$    write(nfout,'(" -- bsdr_l, bsdi_l have been allocated --  np_e, nlmta, ik = ",3i5)') np_e,nlmta,ik
  end subroutine alloc_wfsd_bsdri_3D

  subroutine dealloc_wfsd_bsdri_3D()
    deallocate(wfsd_l)
    deallocate(bsdr_l)
    deallocate(bsdi_l)
  end subroutine dealloc_wfsd_bsdri_3D

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
! === This should be modified after zaj_old is 3D-decomposed!!! by T.Kato ======
    if(allocated(zaj_old)) deallocate(zaj_old)
! ==============================================================================
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
    use m_ES_ortho,           only : mgs_4_each_k_G_3D
#endif
    integer,       intent(in) :: nfout,isolver2 &
         &                     , mode   ! mode = {ORTHONORMALIZATION | NORMALIZATION}
    real(kind=DP), intent(in) :: dtim_old, dtim_new
    integer, intent(in), optional :: iupdate
    integer       :: ispin, ik, iksnl
    real(kind=DP), pointer, dimension(:) :: ekin
#ifdef __TIMER__
    integer       :: ierr
#endif

    real(kind=DP), pointer, dimension(:) :: ekin_l
    real(kind=DP), allocatable, dimension(:) ::  afft_l
    integer             :: lsize
#ifdef FFT_3D_DIVISION
    lsize = fft_X_x_nel*fft_X_y_nel*fft_X_z_nel
    allocate(afft_l(lsize*2), stat=ierr)
#else
    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
    allocate(afft_l(lsize*kimg), stat=ierr)
#endif
     if(ierr /= 0) then
        write(nfout,*)' m_ESsubmat_Renew_WF : Not allocated afft_l array'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 999, ierr)
     endif
    call m_ES_alloc_afft_scss_etc_3D()
    ekin_l => sc_l
    do ispin = 1, nspin, af+1
                                                  __TIMER_START_w_BARRIER(MPI_CommGroup,1)
       call m_ES_Vlocal_in_Rspace_3D(ispin,afft_l,lsize,1,OFF)   ! (ptfft1)
                                                  __TIMER_STOP(1)
       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) == myrank_k ) then                     ! MPI
             call m_pwBS_kinetic_energies(ik,vkxyz,ekin_l) !-(PWBS) (diakin) ->ekin
             call evolve_each_WF_again_3D(ik,isolver2,dtim_new,dtim_old)  !-(m_ES_WF_by_SDorCG)
             if(ipri>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')
#ifdef __TIMER__
             if(modnrm == EXECUT) then
                                                  __TIMER_START_w_BARRIER(mpi_k_world(myrank_k),4)
                call m_ES_betar_dot_WFs_4_each_k_3D(nfout,ik)   ! -> fsr_l,fsi_l
                                                  __TIMER_STOP(4)
                                                  __TIMER_START_w_BARRIER(mpi_k_world(myrank_k),5)
                call mgs_4_each_k_G_3D(ista_k,iend_k,ik,zaj_l,mode,fsr_l,fsi_l,mod_pot=VDB)!-(m_E.S.)
                                                  __TIMER_STOP(5)
             else
                                                  __TIMER_START_w_BARRIER(mpi_k_world(myrank_k),5)
                call mgs_4_each_k_G_3D(ista_k,iend_k,ik,zaj_l,mode,mod_pot=NORMCONSERVATION)!-(m_E.S.)
                                                  __TIMER_STOP(5)
                                                  __TIMER_START_w_BARRIER(mpi_k_world(myrank_k),4)
                call m_ES_betar_dot_WFs_4_each_k_3D(nfout,ik)   ! -> fsr_l,fsi_l
                                                  __TIMER_STOP(4)
             end if
#else
!!$             call m_ES_MGS_4_each_k_3D(nfout,ik,mode)
             call m_ES_MGS_4_each_k(nfout,ik,mode)
#endif
! === Support Hybrid on 3D_Parallel by tkato 2013/02/10 ========================
             if(sw_hybrid_functional==ON) then
                if(sw_retard_eigval_evaluation==OFF)then
                   call m_ES_EXX_eigenvalue_for_each_k(ispin,ik,iupdate=iupdate)
                else
                   call m_ES_EXX_cp_eigenvalue(ik)
                endif
             endif
! ==============================================================================
                                                  __TIMER_START_w_BARRIER(MPI_CommGroup,6)
             call m_ES_eigen_values_for_each_k_3D(ispin,ik,ekin_l,afft_l,lsize)
                                                  __TIMER_STOP(6)
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
    call m_ES_sort_eigen_values_3D !-(m_Elec.)
    if(ipri>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')
    call m_ESsd_sort_zaj_old_3D()
#endif

    deallocate(afft_l)

    call m_ES_dealloc_afft_scss_etc()

  end subroutine m_ESsd_evolve_WFs_again


  subroutine evolve_each_WF_again_3D(ik,isolver2,dt_new,dt_old)
  use m_Parallelization     ,only : np_e => np_e

    integer,       intent(in) :: ik,isolver2
    real(kind=DP), intent(in) :: dt_new,dt_old
    integer ::        ir, ib, ig
    real(kind=DP) ::  dtt
    integer   :: id_sname = -1
    call tstatc0_begin('evolve_each_WF_again_3D ', id_sname,1)

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
             do ig = 1, np_g1k(ik)
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
                do ig = 1, np_g1k(ik)
! ------------------- Revised by T. Yamasaki,  03 July 2008 --->>
                   zaj_l(ig,ib,ik,ir) = zaj_old(ig,ib,ik,ir) + dtt*wfred_l(ig,ib,ik,ir)
! ------------------------------------<<
                end do
             end do
          end do
       else if(.not.wfred_is_allocated) then
          do ir = 1, kimg
             do ib = 1, np_e
                do ig = 1, np_g1k(ik)
                   zaj_l(ig,ib,ik,ir) = (1-dtt)*zaj_old(ig,ib,ik,ir) + dtt*zaj_l(ig,ib,ik,ir)
                end do
             end do
          end do
       end if
! <----------------
    end if
    call tstatc0_end(id_sname)
  end subroutine evolve_each_WF_again_3D

  subroutine m_ESsd_decide_CG_direction(precon)
    integer, intent(in) :: precon

    integer ispin, iksnl, ik
    real(kind=DP), parameter             :: Delta = 1.d-10
    real(kind=DP), pointer, dimension(:) :: ekin, p
    real(kind=DP)                        :: gmgm, gmmgmm, sumdz2, sumdz
    real(kind=DP), pointer, dimension(:) :: ekin_l
    real(kind=DP), allocatable, dimension(:) ::  afft_l, cfft_l
    integer             :: lsize, ierr
#ifdef FFT_3D_DIVISION
    lsize = fft_X_x_nel*fft_X_y_nel*fft_X_z_nel
    allocate(afft_l(lsize*2), stat=ierr)
    if ( use_metagga .and. vtau_exists ) allocate( cfft_l(lsize*2), stat=ierr )
#else
    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
    allocate(afft_l(lsize*kimg), stat=ierr)
    if ( use_metagga .and. vtau_exists ) allocate( cfft_l(lsize*kimg), stat=ierr )
#endif
     if(ierr /= 0) then
        write(nfout,*)' m_ESsubmat_Renew_WF : Not allocated afft_l array'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 999, ierr)
     endif
    if(sw_fef==ON) call m_FEF_build_grad() !->grad_fef
    call m_ES_alloc_afft_scss_etc_3D()
    ekin_l => sc_l


    gmgm = 0.d0; gmmgmm = 0.d0

    do ispin = 1, nspin, (af+1)
       call m_ES_Vlocal_in_Rspace_3D(ispin,afft_l,lsize,1,OFF)   ! (ptfft1) -> afft
       if ( use_metagga .and. vtau_exists ) then
          call m_ES_Vlocal_in_Rspace_3D( ispin, cfft_l, lsize, 1, OFF, vtau_l ) ! r space
       endif

       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) == myrank_k) then           ! MPI
             iksnl = (ik-1)/nspin + 1
             call alloc_wfsd_bsdri_3D(ik)
!FJFJ             call decomp_eko_l_3D(eko_l,eko_l_3D(:,ik),ik)
!FJFJ             call decomp_snl_l_3D(snl,snl_l_3D,ik,nlmtt,iksnl)
!FJFJ             call decomp_fsr_l_3D(fsr_l,fsr_l_3D,ik,nrvf_ordr,'    ')
!FJFJ             if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
!FJFJ                call decomp_fsr_l_3D(fsi_l,fsi_l_3D,ik,nrvf_ordr,'    ')
!FJFJ             endif

             call m_ES_Vnonlocal_W_3D(ik,iksnl,ispin,ON)    ! (nonloc) ->(vnlph_l)

!FJFJ             call decomp_vnlph_l_r_3D(vnlph_l_3D,vnlph_l,ik,1)
!FJFJ             if (kimg == 2) then
!FJFJ                call decomp_vnlph_l_r_3D(vnlph_l_3D,vnlph_l,ik,2)
!FJFJ             endif

! === Support Hybrid on 3D_Parallel by tkato 2013/02/10 ========================
             if(sw_hybrid_functional==ON) call m_ES_Vexx_W(ik) ! (exx) ->(vnlph_l)
! ==============================================================================
             if(sw_fef==ON) call m_FEF_add_grad_to_vnlph(ik) !->vnlph_l

             if ( use_metagga .and. vtau_exists ) then
                call m_ES_contrib_kindens_to_vnlph( ispin, ik, lsize, cfft_l )
             endif

             call m_pwBS_kinetic_energies(ik,vkxyz,ekin_l) ! -(m_PWBS) (diakin) ->ekin
             call decide_CG_direction_core_3D &     ! -(m_ES_WF_by_SDorCG) ->wfsd_l
                  &(precon,ik,ekin_l,afft_l,lsize,sumdz2,sumdz)
!             gmgm   = gmgm + sumdz2
             gmgm   = gmgm + sumdz
             gmmgmm = gmmgmm + dzajn2(ik)
             dzajn2(ik) = sumdz2
             if(ipri >= 2) write(nfout,'(" dzajn2(",i3,") = ",d20.8)') ik,dzajn2(ik)
             call dealloc_wfsd_bsdri_3D()
          end if
       enddo
    enddo

    if(npes > 1) then
       call mpi_allreduce(MPI_IN_PLACE,gmgm,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)    ! MPI
       call mpi_allreduce(MPI_IN_PLACE,gmmgmm,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)    ! MPI
    end if

    if ( allocated(cfft_l) ) deallocate( cfft_l )

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

    deallocate(afft_l)
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



  subroutine m_ESsd_renew_WF_by_SDorCG(nfout,isolver,precon,dtim,iupdate)
    use m_IterationNumbers,     only : iteration
    use m_Control_Parameters,   only : sw_gep
#ifdef __TIMER__
    use m_Const_Parameters,   only : VDB, NORMCONSERVATION
    use m_ES_ortho,           only : mgs_4_each_k_G_3D
#endif
    integer, intent(in) :: nfout,isolver, precon
    integer, intent(in), optional :: iupdate
    real(kind=DP)       :: dtim

    integer                              :: ispin, iksnl, ik, mode, ipri0
    real(kind=DP), pointer, dimension(:) :: ekin, p, vnldi
    real(kind=DP)                        :: vlhxc0
    real(kind=DP), pointer, dimension(:) :: ekin_l
    real(kind=DP), allocatable, dimension(:) ::  afft_l, cfft_l
    integer             :: lsize, ierr, ii
! ==============================================================================
    real(kind=DP) :: p_l(maxval(np_g1k))
    real(kind=DP) :: vnldi_l(maxval(np_g1k))
! ==============================================================================

    call m_ES_alloc_afft_scss_etc_3D()
#ifdef FFT_3D_DIVISION
    lsize = fft_X_x_nel*fft_X_y_nel*fft_X_z_nel
    allocate(afft_l(lsize*2), stat=ierr)
    if ( use_metagga .and. vtau_exists ) allocate(cfft_l(lsize*2), stat=ierr)
#else
    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
    allocate(afft_l(lsize*kimg), stat=ierr)
    if ( use_metagga .and. vtau_exists ) allocate(cfft_l(lsize*kimg), stat=ierr)
#endif
     if(ierr /= 0) then
        write(nfout,*)' m_ESsubmat_Renew_WF : Not allocated afft_l array'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 201, ierr)
     endif
    ekin_l => sc_l
    if(sw_gep == ON.and.iteration>=2) then
      mode = NORMALIZATION
    else
      mode = ORTHONORMALIZATION
    endif

    if(iprisolver >= 2) then
       if(isolver == SD) then
          write(nfout,'(" !! isolver_core = SD")')
       else if(isolver == MSD) then
          write(nfout,'(" !! isolver_core = MSD")')
       else if(isolver == CG) then
          write(nfout,'(" !! isolver_core = CG")')
       end if
    end if

    if(sw_fef==ON) call m_FEF_build_grad() !->grad_fef
    do ispin = 1, nspin, (af+1)
       if(isolver == MSD) then
          call vlhxc_l_zero_term(vlhxc0,ispin) ! -(m_ES_WF_by_SDorCG) ->vlhxc0
       end if
                                                  __TIMER_START_w_BARRIER(MPI_CommGroup,1)
       call m_ES_Vlocal_in_Rspace_3D(ispin,afft_l,lsize,1,OFF)      ! (ptfft1) vlhxc_l->afft
       if ( use_metagga .and. vtau_exists ) then
          call m_ES_Vlocal_in_Rspace_3D( ispin, cfft_l, lsize, 1, OFF, vtau_l ) ! r space
       endif
                                                  __TIMER_STOP(1)
       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) == myrank_k) then           ! MPI
             iksnl = (ik-1)/nspin + 1
                                                  __TIMER_START_w_BARRIER(mpi_k_world(myrank_k),2)
             call m_ES_Vnonlocal_W_3D(ik,iksnl,ispin,ON)    ! (nonloc) ->(vnlph_l)
                                                  __TIMER_STOP(2)
! === Support Hybrid on 3D_Parallel by tkato 2013/02/10 ========================
             if(sw_hybrid_functional==ON) call m_ES_Vexx_W(ik) ! (exx) ->(vnlph_l)
! ==============================================================================
             if(sw_fef==ON) call m_FEF_add_grad_to_vnlph(ik) !->vnlph_l

             if ( use_metagga .and. vtau_exists ) then
                call m_ES_contrib_kindens_to_vnlph( ispin, ik, lsize, cfft_l )
             endif

             call m_pwBS_kinetic_energies(ik,vkxyz,ekin_l) ! (diakin) ->ekin
             if(isolver == SD) then
                call evolve_WFs_in_SD_direction_3D&      ! -(m_ES_WF_by_SDorCG)
                     &(precon,ik,dtim,ekin_l, afft_l,lsize)
             else if(isolver == MSD ) then
                                                  __TIMER_START_w_BARRIER(mpi_k_world(myrank_k),3)
                call evolve_WFs_in_MSD_direction_3D&     !-(m_ES_WF_by_SDorCG)
                     &(precon,ik,iksnl,ispin,dtim,ekin_l,afft_l,lsize,vlhxc0)
                                                  __TIMER_STOP(3)
                if(ipri>=3.and.ik==1) &
                     & call m_ES_wd_zaj_small_portion_3D(nfout,ik," -- after MSD --",16)
             else if(isolver == eazyCG) then
                write(nfout,'("Not support Multi Parallel for eazyCG")')
                call flush(nfout)
                write(1002,'("Not support Multi Parallel for eazyCG")')
                call flush(1002)
!               call evolve_WFs_in_eazyCG_direction(precon,ik,dtim,ekin,afft,bfft,p)
             else if(isolver == CG) then
! ==============================================================================
START_TIMER('CG_Initialize')
!FJFJ           call decomp_wfsd_l(wfsd_old(:,:,ik:ik,:), wfsd_old(:,:,ik:ik,:), ik)
!!              call decomp_fsr_l_2D(fsr_l,fsr_l_2D,ik)
!!              if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
!!                 call decomp_fsr_l_2D(fsi_l,fsi_l_2D,ik)
!!              endif
!===============================================================================
! The reordering "wfsd_old" with "nrvf_ordr/neordr" is necessary
! to obtain the correct results.
! But, such a reordering needs much memory and expends overhead...
! by tkato 2011/08/19
!===============================================================================
                call sort_wfsd_old_in(wfsd_old, ik)
!===============================================================================
STOP_TIMER('CG_Initialize')
! ==============================================================================
!               call evolve_WFs_in_CG_direction_3D(precon,ik,dtim,ekin,afft,bfft,p_l)
                call evolve_WFs_in_CG_direction_3D(precon,ik,dtim,ekin_l,afft_l,lsize)
!!$                call evolve_WFs_in_CG_direction0(precon,ik,dtim,ekin,afft,bfft,p)
! ==============================================================================
START_TIMER('CG_Finalize')
!===============================================================================
                call sort_wfsd_old_out(wfsd_old, ik)
!===============================================================================
!FJFJ           call decomp_wfsd_l_r_3D(wfsd_old(:,:,ik:ik,:), wfsd_old(:,:,ik:ik,:), ik)
STOP_TIMER('CG_Finalize')
! ==============================================================================
             endif
! --------------- Added by T. Yamasaki, 28 June 2008 ---
             if(wfred_is_allocated) then
                call cp_wfred_3D(ik)
             end if
! ------------------------------------------------------<<
#ifdef __TIMER__
             if(modnrm == EXECUT) then
                                                  __TIMER_START_w_BARRIER(mpi_k_world(myrank_k),4)
                call m_ES_betar_dot_WFs_4_each_k_3D(nfout,ik)   ! -> fsr_l,fsi_l
                                                  __TIMER_STOP(4)
                                                  __TIMER_START_w_BARRIER(mpi_k_world(myrank_k),5)
                call mgs_4_each_k_G_3D(ista_k,iend_k,ik,zaj_l,mode,fsr_l,fsi_l,mod_pot=VDB)!-(m_E.S.)
                                                  __TIMER_STOP(5)
             else
                                                  __TIMER_START_w_BARRIER(mpi_k_world(myrank_k),5)
                call mgs_4_each_k_G_3D(ista_k,iend_k,ik,zaj_l,mode,mod_pot=NORMCONSERVATION)!-(m_E.S.)
                                                  __TIMER_STOP(5)
                                                  __TIMER_START_w_BARRIER(mpi_k_world(myrank_k),4)
                call m_ES_betar_dot_WFs_4_each_k_3D(nfout,ik)   ! -> fsr_l,fsi_l
                                                  __TIMER_STOP(4)
             end if
#else
!!$             call m_ES_MGS_4_each_k_3D(nfout,ik,mode)
             call m_ES_MGS_4_each_k(nfout,ik,mode)
#endif
! === Support Hybrid on 3D_Parallel by tkato 2013/02/10 ========================
             if(sw_hybrid_functional==ON) then
                if(sw_retard_eigval_evaluation==OFF)then
                   call m_ES_EXX_eigenvalue_for_each_k(ispin,ik,iupdate=iupdate)
                else
                   call m_ES_EXX_cp_eigenvalue(ik)
                endif
             endif
! ==============================================================================

                                                  __TIMER_START_w_BARRIER(mpi_k_world(myrank_k),6)
             call m_ES_eigen_values_for_each_k_3D(ispin,ik,ekin_l,afft_l,lsize)
                                                  __TIMER_STOP(6)
          end if
       enddo
    enddo

    if(sw_hybrid_functional==ON.and.sw_retard_eigval_evaluation==ON)then
       call m_ES_EXX_gather_valence_states(nfout)
       do ispin=1,nspin,af+1
          do ik=ispin, kv3+ispin-nspin, nspin
             if(map_k(ik) /= myrank_k) cycle ! MPI
             call m_ES_EXX_eigenvalue_for_each_k(ispin,ik,update_eko=.false.,iupdate=iupdate)
          enddo
       enddo
    endif

    if ( allocated(cfft_l) ) deallocate( cfft_l )

#ifdef LMM_PREVIOUS
    call m_ES_sort_eigen_values_3D()   ! -> neordr, nrvf_ordr
    if(ipri>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')
    call m_ESsd_sort_zaj_old_3D()
#endif

    deallocate(afft_l)
    call get_ipri0(ipri, ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko_3D(nfout,mode=SCF)
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





  subroutine evolve_WFs_in_CG_direction_3D(precon,ik,dtim,ekin_l,afft_l,lsize)
    integer, intent(in)        :: precon, ik, lsize
    real(kind=DP), intent(in)  :: dtim,ekin_l(np_g1k(ik))
#ifdef FFT_3D_DIVISION
    real(kind=DP), intent(in)  :: afft_l(lsize*2   )
#else
    real(kind=DP), intent(in)  :: afft_l(lsize*kimg)
#endif
    real(kind=DP)              :: p(maxval(np_g1k))

    integer :: ib
    integer :: id_sname = -1
! === FFT Marge. by T.Kato ===============================================================
    integer :: ib1, ib2, ibsize, ibesize
    integer :: isrsize, fft_l_size
    real(kind=DP), allocatable, dimension(:,:) :: wk_bfft_l
    real(kind=DP), allocatable, dimension(:,:) :: VlocalW
    real(kind=DP), allocatable, dimension(:,:) :: p_l
! ========================================================================================
    call tstatc0_begin('evolve_WFs_in_CG_direction_3D ', id_sname,1)
! === FFT Marge. by T.Kato ===============================================================
    ibsize = 1
    if (nblocksize_fftw_is_given) then
       ibsize = nblocksize_fftw
       if (ibsize < 1) ibsize = 1
    endif
    isrsize = min(lsize,mp_g1k(ik))
    fft_l_size  = nel_fft_x(myrank_g)
#ifdef FFT_3D_DIVISION
    allocate(wk_bfft_l(lsize*2   ,ibsize) ,stat=ierr)
    allocate(VlocalW(lsize*2   ,ibsize) ,stat=ierr)
#else
    allocate(wk_bfft_l(lsize*kimg,ibsize) ,stat=ierr)
    allocate(VlocalW(lsize*kimg,ibsize) ,stat=ierr)
#endif
!P! allocate(p_l(mp_g1k(ik),ibsize) ,stat=ierr)
    allocate(p_l(mp_g1k(ik),np_e) ,stat=ierr)
    if (ierr /= 0) then
       write(nfout,*)' evolve_WFs_in_CG_direction_3D :  Not allocate '
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 205, ierr)
    endif
! ========================================================================================

    call alloc_wfsd_bsdri_3D(ik)

    call m_ES_decide_precon_factor_3D(precon,ik,1,np_e,np_e,ekin_l,p_l)  ! -> p(1:iba(ik))

    do ib1 = 1, np_e, ibsize
       ib2 = min(ib1+ibsize-1,np_e)
       ibesize = ib2 - ib1 + 1
START_TIMER('CG_temporary_FFT')
#ifdef __TIMER_COMM__
       call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,wk_bfft_l, 0)
#else
       call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,wk_bfft_l)
#endif
#ifdef FFT_3D_DIVISION
       call m_FFT_Vlocal_W_3DIV_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_x(myrank_g))
       call m_FFT_Direct_3DIV_3D(nfout,wk_bfft_l,lsize,ibsize)
#else
       if (sw_fft_xzy > 0) then
          call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_y(myrank_g))
          call m_FFT_Direct_3D(nfout,wk_bfft_l,lsize,ibsize)
       else
          call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_z(myrank_g))
          if(sw_serial_fft == ON) then
             call m_ES_WF_2D(ik,wk_bfft_l,ib2,ib1,ibsize,lsize,DIRECT)
          else
             call m_FFT_Direct_XYZ_3D (nfout, wk_bfft_l, lsize, ibsize)
          endif
       end if
#endif
!!$       call map_fft_to_WF_3D(ik,lsize,ibsize,wk_bfft_l,VlocalW,isrsize,fft_l_size)
       call map_fft_to_WF_3D(ik,lsize,ibesize,wk_bfft_l,VlocalW,isrsize,fft_l_size)
STOP_TIMER('CG_temporary_FFT')
START_TIMER('CG_SD_direction_3D')
       call SD_direction_3D(precon,ik,ib1,ib2,ibesize,ekin_l,VlocalW,lsize,p_l)
                                                          ! -(m_ES_WF_by_SDorCG) -> wfsd_l
STOP_TIMER('CG_SD_direction_3D')
       !                wfsd_l : -(H-e_{k\mu}^m)\Psi_{k\mu}^m
    end do

START_TIMER('CG_orthogonalize_SD_drctns_3D')
    call orthogonalize_SD_drctns_3D(ik,to=SAME_BAND) ! -(m_ES_WF_by_SDorCG) ->(wfsd_l,bsd(ri)_l)
STOP_TIMER('CG_orthogonalize_SD_drctns_3D')
START_TIMER('CG_make_CG_direction_3D')
    call make_CG_direction_3D(ik)       ! -(m_ES_WF_by_SDorCG) -> wfsd_l + betacg*wfsd_old -> wfsd_l
STOP_TIMER('CG_make_CG_direction_3D')
!!$    call orthogonalize_SD_drctns(ik,to=SAME_BAND) ! Actually, CG direction ->(wfsd_l,bsd(ri)_l)

START_TIMER('CG_WF_conjugate_gradient_3D')
    do ib = 1, np_e
       call WF_conjugate_gradient_3D(ik,ib,dtim) !-(m_ES_WF_by_SDorCG)
    end do
STOP_TIMER('CG_WF_conjugate_gradient_3D')
    if(ipri>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')

START_TIMER('CG_cp_wfsd_to_wfsd_old')
    call cp_wfsd_to_wfsd_old(ik)    ! -(m_ES_WF_by_SDorCG) wfsd_l ->wfsd_old
STOP_TIMER('CG_cp_wfsd_to_wfsd_old')

    call dealloc_wfsd_bsdri_3D()
! === FFT Marge. by T.Kato ===============================================================
    deallocate(wk_bfft_l)
    deallocate(VlocalW)
    deallocate(p_l)
! ========================================================================================
    call tstatc0_end(id_sname)
  end subroutine evolve_WFs_in_CG_direction_3D


  subroutine vlhxc_l_zero_term(vlhxc0,ispin)
    real(kind=DP), intent(out) :: vlhxc0
    integer, intent(in)        :: ispin

    if(myrank_chg == 0) vlhxc0 = vlhxc_l(1,1,ispin)
    call mpi_bcast(vlhxc0,1,mpi_double_precision,0,mpi_chg_world,ierr)
  end subroutine vlhxc_l_zero_term


  subroutine square_of_SD_direction_3D(ik,ib,dz)
    integer,       intent(in)       :: ik,ib
    real(kind=DP), intent(out)      :: dz  ! dz = <wfsd|w|wfsd>
    integer                         :: i,iadd,ierr
    real(kind=DP),allocatable,dimension(:,:,:) :: bsdrt,bsdit
    dz = 0.d0
    if(modnrm == EXECUT) then
       call m_ES_sum_of_LocalPart(ik,ib,bsdr_l,bsdi_l,dz)
    end if
    if(k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          do i = ista_g1k(ik), iend_g1k(ik)
             if (i == 1) then
                dz = dz + wfsd_l(1,ib,ik,1)*wfsd_l(1,ib,ik,1)
             else
                iadd = i-ista_g1k(ik)+1
                dz = dz + 2.d0*wfsd_l(iadd,ib,ik,1)*wfsd_l(iadd,ib,ik,1)
             end if
          end do
       else
          do i = ista_g1k(ik), iend_g1k(ik)
             if (i == 1) then
                dz = dz + wfsd_l(1,ib,ik,1)**2 + wfsd_l(1,ib,ik,2)**2
             else
                iadd = i-ista_g1k(ik)+1
                dz = dz + 2.d0*(wfsd_l(iadd,ib,ik,1)**2 + wfsd_l(iadd,ib,ik,2)**2)
             end if
          end do
       end if
    else
       if(kimg == 1) then
          do i = ista_g1k(ik), iend_g1k(ik)
             iadd = i-ista_g1k(ik)+1
             dz = dz + wfsd_l(iadd,ib,ik,1)*wfsd_l(iadd,ib,ik,1)
          end do
       else
          do i = ista_g1k(ik), iend_g1k(ik)
             iadd = i-ista_g1k(ik)+1
             dz = dz + wfsd_l(iadd,ib,ik,1)**2 + wfsd_l(iadd,ib,ik,2)**2
          end do
       end if
    end if
  end subroutine square_of_SD_direction_3D

  subroutine square_of_SD_direction3_3D(ik,ibo,dz)
    integer,       intent(in)       :: ik,ibo
    real(kind=DP), intent(out)      :: dz  ! dz = <wfsd|w|wfsd>
    integer                         :: i,ib,iadd
    real(kind=DP),allocatable,dimension(:,:,:) :: bsdrt,bsdit,bsdrto,bsdito

    dz = 0.d0
    if(modnrm == EXECUT) then
       call m_ES_sum_of_LocalPart3(ik,ibo,bsdr_l,bsdi_l,bsdr_old,bsdi_old,dz)
    end if
    ib = map_z(ibo)                     ! MPI
    if(k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          do i = ista_g1k(ik),iend_g1k(ik)
             if(i == 1)then
                dz = dz + wfsd_l(1,ib,ik,1)*wfsd_old(1,ib,ik,1)
             else
                iadd = i-ista_g1k(ik)+1
                dz = dz + 2.d0*wfsd_l(iadd,ib,ik,1)*wfsd_old(iadd,ib,ik,1)
             endif
          end do
       else
          do i = ista_g1k(ik),iend_g1k(ik)
             if(i == 1)then
                dz = dz + wfsd_l(1,ib,ik,1)*wfsd_old(1,ib,ik,1) + wfsd_l(1,ib,ik,2)*wfsd_old(1,ib,ik,2)
             else
                iadd = i-ista_g1k(ik)+1
                dz = dz + 2.d0*(wfsd_l(iadd,ib,ik,1)*wfsd_old(iadd,ib,ik,1) &
                &  + wfsd_l(iadd,ib,ik,2)*wfsd_old(iadd,ib,ik,2))
             endif
          end do
       end if
    else
       if(kimg == 1) then
          do i = ista_g1k(ik),iend_g1k(ik)
             iadd = i-ista_g1k(ik)+1
             dz = dz + wfsd_l(iadd,ib,ik,1)*wfsd_old(iadd,ib,ik,1)
          end do
       else
          do i = ista_g1k(ik),iend_g1k(ik)
             iadd = i-ista_g1k(ik)+1
             dz = dz + wfsd_l(iadd,ib,ik,1)*wfsd_old(iadd,ib,ik,1) + wfsd_l(iadd,ib,ik,2)*wfsd_old(iadd,ib,ik,2)
          end do
       end if
    end if
  end subroutine square_of_SD_direction3_3D


  subroutine SD_direction_3D(precon,ik,ib1,ib2,ibesize,ekin_l,VlocalW,lsize,p_l)
    integer      , intent(in)                  :: precon,ik,ib1,ib2,ibesize,lsize
    real(kind=DP), intent(in), dimension(np_g1k(ik))  :: ekin_l
#ifdef FFT_3D_DIVISION
    real(kind=DP), intent(in), dimension(lsize*2   ,ibesize) :: VlocalW
#else
    real(kind=DP), intent(in), dimension(lsize*kimg,ibesize) :: VlocalW
#endif
    real(kind=DP)            , dimension(mp_g1k(ik),np_e) :: p_l

    integer       :: i, i1, ib, iadd
    real(kind=DP) :: devr,denom, e1, devi

    denom = 1.d0/product(fft_box_size_WF(1:3,1))

    if(kimg == 1) then
       do i = ista_g1k(ik), iend_g1k(ik)
          iadd = i - ista_g1k(ik) + 1
!xx       i1    = igf(nbase(i,ik))
          do ib = ib1, ib2
             devr  = (ekin_l(iadd)-eko_l(ib,ik))*zaj_l(iadd,ib,ik,1)&
                  & + VlocalW(iadd,ib-ib1+1)*denom + vnlph_l(iadd,ib,1)
             wfsd_l(iadd,ib,ik,1) = - p_l(iadd,ib)*devr
          end do
       end do
    else if(kimg == 2) then
       do i = ista_g1k(ik), iend_g1k(ik)
          iadd = i - ista_g1k(ik) + 1
!xx       i1    = igf(nbase(i,ik))
          do ib = ib1, ib2
             e1    = ekin_l(iadd) - eko_l(ib,ik)
             devr  = e1*zaj_l(iadd,ib,ik,1) + VlocalW(2*iadd-1,ib-ib1+1)*denom+vnlph_l(iadd,ib,1)
             devi  = e1*zaj_l(iadd,ib,ik,2) + VlocalW(2*iadd  ,ib-ib1+1)*denom+vnlph_l(iadd,ib,2)
             wfsd_l(iadd,ib,ik,1) = - p_l(iadd,ib)*devr
             wfsd_l(iadd,ib,ik,2) = - p_l(iadd,ib)*devi
          end do
       end do
    end if
  end subroutine SD_direction_3D

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

    ib = ibo
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

    ib = ibo
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

  subroutine WF_conjugate_gradient_3D(ik,ibo,dtim)
! \Psi_{k\mu}^{m+1} = \Psi_{k\mu}^{m} + \Delta t_{opt}^{m} d^{m}
!   d^{m} = g^{m} + \beta_CG d^{m-1}
!   g^{m} = - \Phi_{k\mu}^{m}
!   d^{m-1} = (\Psi_{k\mu}^{m} - \Psi_{k\mu}^{m-1})\frac{1}{\Delta t_{opt}^{m-1}}
!
    integer      , intent(in)                  :: ik,ibo
    real(kind=DP), intent(in)                  :: dtim

    integer       :: i,ib,iadd
    real(kind=DP) :: evr,evi

    ib = ibo                    ! MPI
#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) status_saved_phifftr(ib,ik) = OLD
#endif
    if(kimg == 1) then
       do i = ista_g1k(ik), iend_g1k(ik)
          iadd = i - ista_g1k(ik) + 1
          evr   = zaj_l(iadd,ib,ik,1)
          zaj_l(iadd,ib,ik,1) = evr + dtim*wfsd_l(iadd,ib,ik,1)
       end do
    else if(kimg == 2) then
       do i = ista_g1k(ik), iend_g1k(ik)
          iadd = i - ista_g1k(ik) + 1
          evr   = zaj_l(iadd,ib,ik,1);  evi   = zaj_l(iadd,ib,ik,kimg)
          zaj_l(iadd,ib,ik,1) = evr + dtim*wfsd_l(iadd,ib,ik,1)
          zaj_l(iadd,ib,ik,2) = evi + dtim*wfsd_l(iadd,ib,ik,2)
       end do
    end if
  end subroutine WF_conjugate_gradient_3D


  subroutine make_CG_direction_3D(ik)
    integer, intent(in) :: ik

    real(kind=DP), parameter             :: Delta = 1.d-20
    integer :: i, ib, iadd
    integer :: id_sname = -1
    call tstatc0_begin('make_CG_direction_3D', id_sname)

    if(betacg > Delta) then
       if(kimg == 1) then
          do ib = 1, np_e
             do i = ista_g1k(ik), iend_g1k(ik)
                iadd = i - ista_g1k(ik) + 1
                wfsd_l(iadd,ib,ik,1) = wfsd_l(iadd,ib,ik,1) + betacg*wfsd_old(iadd,ib,ik,1)
             end do
          end do
       else
          do ib = 1, np_e
             do i = ista_g1k(ik), iend_g1k(ik)
                iadd = i - ista_g1k(ik) + 1
                wfsd_l(iadd,ib,ik,1) = wfsd_l(iadd,ib,ik,1) + betacg*wfsd_old(iadd,ib,ik,1)
                wfsd_l(iadd,ib,ik,2) = wfsd_l(iadd,ib,ik,2) + betacg*wfsd_old(iadd,ib,ik,2)
             end do
          end do
       end if
    end if
    call tstatc0_end(id_sname)
  end subroutine make_CG_direction_3D

  subroutine orthogonalize_SD_drctns_3D(ik,to)
    integer, intent(in) :: ik,to

    integer :: id_sname = -1
    call tstatc0_begin('orthogonalize_SD_drctns_3D ', id_sname)

!!$    if(modnrm == EXECUT) call m_ES_betar_dot_Psi_4_each_k(wfsd_l,ista_k,iend_k,ik,bsdr_l,bsdi_l)
START_TIMER('CG2_m_ES_betar_dot_Psi_4_each_k_3D')
    if(modnrm == EXECUT) call m_ES_betar_dot_Psi_4_each_k_3D(nfout,wfsd_l,ik,ik,ik,bsdr_l,bsdi_l)
STOP_TIMER('CG2_m_ES_betar_dot_Psi_4_each_k_3D')
    !                                                          ->bsd(ri)_l
START_TIMER('CG2_m_ES_orthogonal_phi_to_WFs')
!!$    call m_ES_orthogonalize_SD_to_WFs_3D(ik,to,wfsd_l,bsdr_l,bsdi_l) ! ->(wfsd_l,bsd(ri)_l)
    call m_ES_orthogonal_phi_to_WFs(ik,wfsd_l,bsdr_l,bsdi_l) ! ->(wfsd_l,bsd(ri)_l)
STOP_TIMER('CG2_m_ES_orthogonal_phi_to_WFs')

!!$    call orthogonalize_SD(ik,wfsd_l,bsdr_l,bsdi_l)           !-(m_E.S.) ->(wfsd_l)
    call tstatc0_end(id_sname)
  end subroutine orthogonalize_SD_drctns_3D

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

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine evolve_WFs_in_SD_direction_3D&
       &(precon,ik,dtim,ekin_l,afft_l,lsize)
    integer, intent(in)        :: precon, ik, lsize
#ifdef FFT_3D_DIVISION
    real(kind=DP), intent(in)  :: dtim, ekin_l(np_g1k(ik)), afft_l(lsize*2   )
#else
    real(kind=DP), intent(in)  :: dtim, ekin_l(np_g1k(ik)), afft_l(lsize*kimg)
#endif

    real(kind=DP), allocatable,dimension(:,:) :: wk_bfft_l
    real(kind=DP), allocatable,dimension(:,:) :: bfft_l
    real(kind=DP), allocatable,dimension(:,:) :: p_l

    integer :: ib1, ib2, ibsize, ibesize
    integer :: isrsize, fft_l_size

    integer :: id_sname = -1
    call tstatc0_begin('evolve_WFs_in_SD_direction_3D ', id_sname,1)

    ibsize = 1
    if (nblocksize_fftw_is_given) then
       ibsize = nblocksize_fftw
       if (ibsize < 1) ibsize = 1
    endif
    isrsize = min(lsize,mp_g1k(ik))
    fft_l_size  = nel_fft_x(myrank_g)

#ifdef FFT_3D_DIVISION
    allocate(wk_bfft_l(lsize*2   ,ibsize) ,stat=ierr)
    allocate(bfft_l(lsize*2   ,ibsize) ,stat=ierr)
#else
    allocate(wk_bfft_l(lsize*kimg,ibsize) ,stat=ierr)
    allocate(bfft_l(lsize*kimg,ibsize) ,stat=ierr)
#endif
!P! allocate(p_l(mp_g1k(ik),ibsize) ,stat=ierr)
    allocate(p_l(mp_g1k(ik),np_e) ,stat=ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_ES_WF_in_Rspace_3D :  Not allocate '
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 203, ierr)
     endif
!P!
    call m_ES_decide_precon_factor_3D(precon,ik,1,np_e,np_e,ekin_l,p_l)! -> p(1:iba(ik))
!P!
    do ib1 = 1, np_e, ibsize
       ib2 = min(ib1+ibsize-1,np_e)
       ibesize = ib2 - ib1 + 1

#ifdef __TIMER_COMM__
       call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,wk_bfft_l, 0)
#else
       call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,wk_bfft_l)
#endif
#ifdef FFT_3D_DIVISION
       call m_FFT_Vlocal_W_3DIV_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_x(myrank_g))
       call m_FFT_Direct_3DIV_3D(nfout,wk_bfft_l,lsize,ibsize)
#else
       if (sw_fft_xzy > 0) then
          call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_y(myrank_g))
          call m_FFT_Direct_3D(nfout,wk_bfft_l,lsize,ibsize)
       else
          call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_z(myrank_g))
          if(sw_serial_fft == ON) then
             call m_ES_WF_2D(ik,wk_bfft_l,ib2,ib1,ibsize,lsize,DIRECT)
          else
             call m_FFT_Direct_XYZ_3D (nfout, wk_bfft_l, lsize, ibsize)
          endif
       end if
#endif

       call map_fft_to_WF_3D(ik, lsize, ibesize, wk_bfft_l, bfft_l, isrsize, fft_l_size)

       call steepest_descent_3D(precon,ik,ib1,ib2,ibesize,dtim,ekin_l,bfft_l,lsize,p_l)
                                         ! -(m_ES_WF_by_SDorCG)
    end do
    if(ipri>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')

    deallocate(wk_bfft_l)
    deallocate(bfft_l)
    deallocate(p_l)

    call tstatc0_end(id_sname)
  end subroutine evolve_WFs_in_SD_direction_3D

!===============================================================================
  subroutine steepest_descent_3D(precon,ik,ib1,ib2,ibesize,dtim,ekin_l,VlocalW,lsize,p_l)
    integer      , intent(in)                  :: precon, ik, ib1, ib2, ibesize, lsize
    real(kind=DP), intent(in)                  :: dtim
    real(kind=DP), intent(in), dimension(np_g1k(ik))  :: ekin_l
#ifdef FFT_3D_DIVISION
    real(kind=DP), intent(in), dimension(lsize*2   ,ibesize) :: VlocalW
#else
    real(kind=DP), intent(in), dimension(lsize*kimg,ibesize) :: VlocalW
#endif
!P! real(kind=DP)            , dimension(mp_g1k(ik),ibesize)  :: p_l
    real(kind=DP)            , dimension(mp_g1k(ik),np_e)  :: p_l

    integer       :: iadd, i, ib
    real(kind=DP) :: evr,devr,denom, evi,e1, devi

    denom = 1.d0/product(fft_box_size_WF(1:3,1))

#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) status_saved_phifftr(ib1:ib2,ik) = OLD
#endif
    if(kimg == 1) then
       do i=ista_g1k(ik), iend_g1k(ik)
          iadd = i-ista_g1k(ik)+1
          do ib = ib1, ib2
             evr   = zaj_l(iadd,ib,ik,1)
             e1    = ekin_l(iadd)-eko_l(ib,ik)
             devr  = e1*evr + VlocalW(iadd,ib-ib1+1)*denom + vnlph_l(iadd,ib,1)
!P!          zaj_l(iadd,ib,ik,1) = evr - p_l(iadd,ib-ib1+1)*dtim*devr
             zaj_l(iadd,ib,ik,1) = evr - p_l(iadd,ib)*dtim*devr
          enddo
       end do
    else if(kimg == 2) then
       do i=ista_g1k(ik), iend_g1k(ik)
          iadd = i-ista_g1k(ik)+1
          do ib = ib1, ib2
             evr   = zaj_l(iadd,ib,ik,1)
             evi   = zaj_l(iadd,ib,ik,2)
             e1    = ekin_l(iadd) - eko_l(ib,ik)
             devr  = e1*evr + VlocalW(2*iadd-1,ib-ib1+1)*denom + vnlph_l(iadd,ib,1)
             devi  = e1*evi + VlocalW(2*iadd  ,ib-ib1+1)*denom + vnlph_l(iadd,ib,2)
!P!          zaj_l(iadd,ib,ik,1) = evr - p_l(iadd,ib-ib1+1)*dtim*devr
!P!          zaj_l(iadd,ib,ik,2) = evi - p_l(iadd,ib-ib1+1)*dtim*devi
             zaj_l(iadd,ib,ik,1) = evr - p_l(iadd,ib)*dtim*devr
             zaj_l(iadd,ib,ik,2) = evi - p_l(iadd,ib)*dtim*devi
         end do
       end do
    end if
  end subroutine steepest_descent_3D

!===============================================================================
  subroutine evolve_WFs_in_MSD_direction_3D  &
 &           (precon,ik,iksnl,ispin,dtim,ekin_l,afft_l,lsize,vlhxc0)
    integer, intent(in)        :: precon, ik, iksnl, ispin, lsize
#ifdef FFT_3D_DIVISION
    real(kind=DP), intent(in)  :: dtim, ekin_l(np_g1k(ik)), afft_l(lsize*2), vlhxc0
#else
    real(kind=DP), intent(in)  :: dtim, ekin_l(np_g1k(ik)), afft_l(lsize*kimg), vlhxc0
#endif

    real(kind=DP), allocatable,dimension(:,:) :: p_l, vnldi_l

    real(kind=DP), allocatable,dimension(:,:) :: wk_bfft_l
    real(kind=DP), allocatable,dimension(:,:) :: bfft_l

    integer :: isrsize, fft_l_size
    integer ib1, ib2, ibsize, ibesize, ib, ig
    integer :: id_sname = -1
    integer, save :: print_ibsize = ON
! === Support Hybrid on 3D_Parallel by tkato 2013/02/10 ========================
    real(kind=DP), allocatable :: vxdi(:) !d(kg1)
! ==============================================================================
                                                  __TIMER_SUB_START(301)

    call tstatc0_begin('evolve_WFs_in_MSD_direction_3D ', id_sname,1)
! === Support Hybrid on 3D_Parallel by tkato 2013/02/10 ========================
    if(sw_hybrid_functional==ON) then
       allocate(vxdi(maxval(np_g1k)))
       call m_ES_EXX_Diagonal_part(ispin,ik,vxdi)
    end if
! ==============================================================================

    ibsize = 1
    if (nblocksize_fftw_is_given) then
       ibsize = nblocksize_fftw
       if (ibsize < 1) ibsize = 1
    endif
    if(ipri>=1 .and. print_ibsize == ON) then
       write(nfout,'(" !! ibsize in evolve_WFs_in_MSD_direction_3D = ",i8)') ibsize
       print_ibsize = OFF
    end if
    isrsize = min(lsize,mp_g1k(ik))
    fft_l_size  = nel_fft_x(myrank_g)

#ifdef FFT_3D_DIVISION
    allocate(wk_bfft_l(lsize*2   ,ibsize) ,stat=ierr)
    allocate(bfft_l(lsize*2   ,ibsize) ,stat=ierr)
#else
    allocate(wk_bfft_l(lsize*kimg,ibsize) ,stat=ierr)
    allocate(bfft_l(lsize*kimg,ibsize) ,stat=ierr)
#endif
!P! allocate(p_l(mp_g1k(ik),ibsize), stat=ierr)
    allocate(p_l(mp_g1k(ik),np_e), stat=ierr)
    allocate(vnldi_l(mp_g1k(ik),ibsize), stat=ierr)
     if (ierr /= 0) then
        write(nfout,*)' evolve_WFs_in_MSD_direction_3D :  Not allocate '
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 205, ierr)
     endif
!P!
    call m_ES_decide_precon_factor_3D(precon,ik,1,np_e,np_e,ekin_l,p_l)! -> p(1:iba(ik))
!P!
    do ib1 = 1, np_e, ibsize
       ib2 = min(ib1+ibsize-1,np_e)
       ibesize = ib2 - ib1 + 1

#ifdef __TIMER_COMM__
       call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,wk_bfft_l, 3)
#else
       call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,wk_bfft_l)
#endif

       if(ipri >= 2) then
          write(nfout,'(" --- afft and bfft <<evolve_WFs_in_MSD_direction_3D>>,  ib1 = ",i5)') ib1
          write(nfout,'(" afft: ",5d16.8)') (afft_l(ig),ig=1,min(5,lsize))
          write(nfout,'(" bfft: ",5d16.8)') (wk_bfft_l(ig,1),ig=1,min(5,lsize))
       end if

#ifdef FFT_3D_DIVISION
       call m_FFT_Vlocal_W_3DIV_3D(afft_l,wk_bfft_l,lsize,ibesize,nel_fft_x(myrank_g))
#else
       if (sw_fft_xzy > 0) then
          call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibesize,nel_fft_y(myrank_g))
       else
          call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibesize,nel_fft_z(myrank_g))
       end if
#endif
                                                  __TIMER_COMM_START(328)
#ifdef FFT_3D_DIVISION
       call m_FFT_Direct_3DIV_3D (nfout, wk_bfft_l, lsize, ibesize)
#else
       if (sw_fft_xzy > 0) then
          call m_FFT_Direct_3D (nfout, wk_bfft_l, lsize, ibesize)
       else
          if(sw_serial_fft == ON) then
             call m_ES_WF_2D(ik,wk_bfft_l,ib2,ib1,ibsize,lsize,DIRECT)
          else
             call m_FFT_Direct_XYZ_3D (nfout, wk_bfft_l, lsize, ibsize)
          endif
       end if
#endif
                                                  __TIMER_COMM_STOP(328)

       call Vnonlocal_Diagonal_part_3D(ispin,ik,iksnl,ib1,ib2,ibesize,vnldi_l)
                                         ! -(m_ES_WF_by_SDorCG) (nonldj)

       call map_fft_to_WF_3D(ik, lsize, ibesize, wk_bfft_l, bfft_l, isrsize, fft_l_size)

! === Support Hybrid on 3D_Parallel by tkato 2013/02/10 ========================
       if(sw_hybrid_functional==ON) call m_ES_EXX_add_Diagonal_part(ik,ib,vxdi,vnldi_l)
! ==============================================================================
       if(ipri >= 2) then
        do ib = ib1, ib2
          write(nfout,'(" --- vnlph_l <<evolve_WFs_in_MSD_direction>>,  ib = ",i5)') ib
                      write(nfout,'(" vnlph_l: Re ",5d16.8)') (vnlph_l(ig,ib,1),ig=1,min(5,iba(ik)))
          if(kimg== 2) then
             write(nfout,'(" vnlph_l: Im ",5d16.8)') (vnlph_l(ig,ib,2),ig=1,min(5,iba(ik)))
          end if
          write(nfout,'(" vnldi ",5d16.8)') (vnldi_l(ig,1),ig=1,min(5,iba(ik)))
          write(nfout,'(" VlocalW: Re ",5d16.8)') (bfft_l(ig,1),ig=1,min(5,iba(ik)))
          write(nfout,'(" VlocalW: Im ",5d16.8)') (bfft_l(ig,1),ig=1,min(5,iba(ik)))
          write(nfout,'(" vlhxc0 = ",5d16.8)') vlhxc0
        end do
       end if

       call modified_steepest_descent_3D&   ! -(m_ES_WF_by_SDorCG)
            &(precon,ik,ib1,ib2,ibesize,dtim,vnldi_l,vlhxc0,ekin_l,bfft_l,lsize,p_l)
    end do
    if(ipri>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')

    deallocate(wk_bfft_l)
    deallocate(bfft_l)
    deallocate(p_l)
    deallocate(vnldi_l)
! === Support Hybrid on 3D_Parallel by tkato 2013/02/10 ========================
    if(sw_hybrid_functional==ON) deallocate(vxdi)
! ==============================================================================

    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(301)
  end subroutine evolve_WFs_in_MSD_direction_3D

!===============================================================================
  subroutine Vnonlocal_Diagonal_part_3D(ispin,ik,iksnl,ib1,ib2,ibesize,vnldi_l)
    integer, intent(in)                        :: ispin, ik, iksnl, ib1, ib2, ibesize
    real(kind=DP), intent(out), dimension(mp_g1k(ik),ibesize) :: vnldi_l
    integer :: it,mdvdb,ib,    i,j

#ifdef VNONLOCAL_SPEED
    integer :: p1, p2, lmtt1, lmtt2, il1, il2, icnt
    integer, save :: icnt_lmt = 0
    logical, save :: first_call = .true.
#endif
                                                  __TIMER_SUB_START(305)

#ifdef VNONLOCAL_SPEED
                                                  __TIMER_DO_START(316)
    if (first_call) then
      do it = 1, ntyp
         icnt = 0
         do p1 = 1, ilmt(it)
            lmtt1 = lmtt(p1,it)
            il1 = ltp(p1,it)
            do p2 = p1, ilmt(it)
               lmtt2 = lmtt(p2,it)
               il2 = ltp(p2,it)
               if(mod(il1+il2,2) == 1) cycle
               icnt = icnt + 1
            end do
         end do
         if (icnt_lmt < icnt) then
            icnt_lmt = icnt
         end if
      end do
      first_call = .false.
      write(nfout,'("Vnonlocal_Diagonal_part_3D ntyp=",i8)') ntyp
      write(nfout,'("Vnonlocal_Diagonal_part_3D icnt_lmt=",i8)') icnt_lmt
    end if
                                                  __TIMER_DO_STOP(316)
#endif

    vnldi_l = 0.d0
    do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(mdvdb == SKIP) then
          call Vnonlocal_D_norm_conserve_case
       else if(mdvdb == EXECUT) then
          call Vnonlocal_D_vanderbilt_case
       end if
    end do
                                                  __TIMER_SUB_STOP(305)
  contains

#ifdef VNONLOCAL_SPEED
    subroutine Vnonlocal_D_vanderbilt_case
      integer       :: ia, p1,p2,lmtt1,il1,lmtt2,il2,i,iadd,icnt
      real(kind=DP) :: ph, dion_eq
      real(kind=DP) :: fac_l(ib1:ib2,icnt_lmt), fac
                                                  __TIMER_SUB_START(307)

      if(ipaw(it)==0) then
         icnt = 0
         fac_l = 0.0d0
                                                  __TIMER_DO_START(319)
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
               icnt = icnt + 1

               do ib = ib1, ib2
                  do ia = ista_atm, iend_atm
                     if(ityp(ia) /= it) cycle
                     dion_eq = dion(p1,p2,it)-eko_l(ib,ik)*q(p1,p2,it)
                     fac_l(ib,icnt) = fac_l(ib,icnt) + ph*iwei(ia) * (dion_eq+vlhxcQ(p1,p2,ia,ispin))
                  end do
               end do
            end do
         end do
                                                  __TIMER_DO_STOP(319)
                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,320)
         if (nrank_g > 1) &
              & call mpi_allreduce(MPI_IN_PLACE,fac_l,(ib2-ib1+1)*icnt_lmt,mpi_double_precision &
              &                   ,mpi_sum,mpi_ke_world,ierr)
                                                  __TIMER_COMM_STOP(320)
         icnt = 0
                                                  __TIMER_DO_START(321)
         do p1 = 1, ilmt(it)
            lmtt1 = lmtt(p1,it)
            il1 = ltp(p1,it)
            do p2 = p1, ilmt(it)
               lmtt2 = lmtt(p2,it)
               il2 = ltp(p2,it)
               if(mod(il1+il2,2) == 1) cycle
               icnt = icnt + 1

               do ib = ib1, ib2
                  do i = ista_g1k(ik), iend_g1k(ik)
                     iadd = i-ista_g1k(ik)+1
                     vnldi_l(iadd,ib-ib1+1) = vnldi_l(iadd,ib-ib1+1) + fac_l(ib,icnt) &
                    &                         *snl(iadd,lmtt1,iksnl)*snl(iadd,lmtt2,iksnl)
                  end do
               end do
            end do
         end do
                                                  __TIMER_DO_STOP(321)
      else
                                                  __TIMER_DO_START(322)
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
               do ib = ib1, ib2
                  fac = 0.d0
                  do ia = natm, natm
                     if(ityp(ia) /= it) cycle
                     fac = ph*iwei(ia) * (dion_paw(p1,p2,ispin,ia) &
                    &                   -eko_l(ib,ik)*q(p1,p2,it)+vlhxcQ(p1,p2,ia,ispin))
                  end do
                  do i = ista_g1k(ik), iend_g1k(ik)
                        iadd = i-ista_g1k(ik)+1
                        vnldi_l(iadd,ib-ib1+1) = vnldi_l(iadd,ib-ib1+1) &
                       &                + fac*snl(iadd,lmtt1,iksnl)*snl(iadd,lmtt2,iksnl)
                  end do
               end do
            end do
         end do
                                                  __TIMER_DO_STOP(322)
      end if
                                                  __TIMER_SUB_STOP(307)
    end subroutine Vnonlocal_D_vanderbilt_case
#else
    subroutine Vnonlocal_D_vanderbilt_case
      integer       :: ia, p1,p2,lmtt1,il1,lmtt2,il2,i,iadd
      real(kind=DP) :: ph,fac, dion_eq
                                                  __TIMER_SUB_START(307)

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

            do ib = ib1, ib2
               dion_eq = dion(p1,p2,it)-eko_l(ib,ik)*q(p1,p2,it)
               fac = 0.d0
               do ia = 1, natm
                  if(ityp(ia) /= it) cycle
                  if(ipaw(it)==0) then
                      fac = fac + ph*iwei(ia) * (dion_eq+vlhxcQ(p1,p2,ia,ispin))
                  else
! ==== DEBUG by tkato 2012/11/07 ===============================================
!                     fac = ph*iwei(ia) * (dion_paw(p1,p2,ispin,ia) &
                      fac = fac + ph*iwei(ia) * (dion_paw(p1,p2,ispin,ia) &
! ==============================================================================
                     &                     -eko_l(ib,ik)*q(p1,p2,it)+vlhxcQ(p1,p2,ia,ispin))
                  end if
               end do

               do i = ista_g1k(ik), iend_g1k(ik)
                  iadd = i-ista_g1k(ik)+1
                  vnldi_l(iadd,ib-ib1+1) = vnldi_l(iadd,ib-ib1+1) &
                 &                       + fac*snl(iadd,lmtt1,iksnl)*snl(iadd,lmtt2,iksnl)
               end do
            end do

         end do
      end do
                                                  __TIMER_SUB_STOP(307)
    end subroutine Vnonlocal_D_vanderbilt_case
#endif
!! #endif VNONLOCAL_SPEED

    subroutine Vnonlocal_D_norm_conserve_case
      integer       :: ia, lmt1,lmt2,lmtt1,il1,im1,il2,im2,i,iadd
      real(kind=DP) :: ph,fac
                                                  __TIMER_SUB_START(306)
      ph = 0.d0
                                                  __TIMER_DO_START(317)
      do ia = 1, natm
         if(ityp(ia) /= it) cycle
         ph = ph + iwei(ia)
      end do
                                                  __TIMER_DO_STOP(317)
                                                  __TIMER_DO_START(318)
      do lmt1 = 1, ilmt(it)
         lmtt1 = lmtt(lmt1,it); il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
         do lmt2 = lmt1, ilmt(it)
            il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
            if(il1 /= il2 .or. im1 /= im2) cycle
            if(mod(il1+il2,2) == 1) cycle
! ==== DEBUG by tkato 2012/11/07 ===============================================
!           fac = ph * dion(lmt1,lmt2,it)
            fac = 0.d0
            Do ia=1, natm
              if(ityp(ia) /= it) cycle
              if(ipaw(it)==0)then
                 fac = fac + iwei(ia) * dion(lmt1,lmt2,it)
              else
                 fac = fac + iwei(ia) * dion_paw(lmt1,lmt2,ispin,ia)
              endif
            End do
! ==============================================================================
            do ib = ib1, ib2
               do i = ista_g1k(ik), iend_g1k(ik)
                  iadd = i-ista_g1k(ik)+1
                  vnldi_l(iadd,ib-ib1+1)  = vnldi_l(iadd,ib-ib1+1)  &
                 &                           + fac * snl(iadd,lmtt1,iksnl)**2
               end do
            end do
         end do
      end do
                                                  __TIMER_DO_STOP(318)
                                                  __TIMER_SUB_STOP(306)
    end subroutine Vnonlocal_D_norm_conserve_case
  end subroutine Vnonlocal_Diagonal_part_3D

!===============================================================================
  subroutine modified_steepest_descent_3D&
       &(precon,ik,ib1,ib2,ibesize,dtim,vnldi_l,vlhxc0,ekin_l,VlocalW,lsize,p_l)
    integer      , intent(in)                  :: precon,ik,ib1,ib2,ibesize,lsize
    real(kind=DP), intent(in)                  :: dtim
    real(kind=DP), intent(in)                  :: vlhxc0
    real(kind=DP), intent(in), dimension(np_g1k(ik))  :: ekin_l
#ifdef FFT_3D_DIVISION
    real(kind=DP), intent(in), dimension(lsize*2   ,ibesize) :: VlocalW
#else
    real(kind=DP), intent(in), dimension(lsize*kimg,ibesize) :: VlocalW
#endif
    real(kind=DP), intent(in), dimension(mp_g1k(ik),ibesize)  :: vnldi_l
!P! real(kind=DP)            , dimension(mp_g1k(ik),ibesize)  :: p_l
    real(kind=DP)            , dimension(mp_g1k(ik),np_e)  :: p_l

    integer       :: i, ib, iadd,j
    real(kind=DP) :: evr,devr,denom, wdi, evi,e1, devi, fdexp
                                                  __TIMER_SUB_START(309)

    denom = 1.d0/product(fft_box_size_WF(1:3,1))
    if(ipri >= 3) then
       write(nfout,'(" ekin : ",5d16.8)') (ekin_l(i),i=1,np_g1k(ik))
       write(nfout,'(" p    : ",5d16.8)') (p_l(i,1),i=1,iba(ik))
    end if

                                                  __TIMER_DO_START(326)
#ifdef SAVE_FFT_TIMES
  if(sw_save_fft == ON) status_saved_phifftr(ib1:ib2,ik) = OLD
#endif
    if(kimg == 1) then
!OCL NORECURRENCE
!OCL NOFLTLD
       do i = ista_g1k(ik), iend_g1k(ik)
          iadd = i-ista_g1k(ik)+1
          do ib = ib1, ib2
             evr   = zaj_l(iadd,ib,ik,1)
             e1    = ekin_l(iadd)-eko_l(ib,ik)
             devr  = e1*evr + VlocalW(iadd,ib-ib1+1)*denom + vnlph_l(iadd,ib,1)
             wdi   = ekin_l(iadd) + vlhxc0 + vnldi_l(iadd,ib-ib1+1) - eko_l(ib,ik)
!P!          fdexp = dexp( -p_l(iadd,ib-ib1+1) * wdi * dtim)
! === 0 divide occurs if abs(wdi) < epsilon. by tkato 2012/12/18 ===============
!            fdexp = dexp( -p_l(iadd,ib) * wdi * dtim)
!            zaj_l(iadd,ib,ik,1) = (fdexp - 1)*devr/wdi + evr
             if(dabs(wdi) < SmallestPositiveNumber) then
                zaj_l(iadd,ib,ik,1) = -p_l(iadd,ib)*devr*dtim + evr
             else
                fdexp = dexp( -p_l(iadd,ib) * wdi * dtim)
                zaj_l(iadd,ib,ik,1) = (fdexp - 1)*devr/wdi + evr
             end if
! ==============================================================================
          enddo
       end do
    else if(kimg == 2) then
!OCL NORECURRENCE
!OCL NOFLTLD
       do i = ista_g1k(ik), iend_g1k(ik)
          iadd = i-ista_g1k(ik)+1
          do ib = ib1, ib2
             evr   = zaj_l(iadd,ib,ik,1)
             evi   = zaj_l(iadd,ib,ik,2)
             e1    = ekin_l(iadd) - eko_l(ib,ik)
             devr  = e1*evr + VlocalW(2*iadd-1,ib-ib1+1)*denom + vnlph_l(iadd,ib,1)
             devi  = e1*evi + VlocalW(2*iadd,  ib-ib1+1)*denom + vnlph_l(iadd,ib,2)
             wdi   = ekin_l(iadd) + vlhxc0 + vnldi_l(iadd,ib-ib1+1) - eko_l(ib,ik)
!P!          fdexp = dexp( -p_l(iadd,ib-ib1+1) * wdi * dtim)
! === 0 divide occurs if abs(wdi) < epsilon. by tkato 2012/12/18 ===============
!            fdexp = dexp( -p_l(iadd,ib) * wdi * dtim)
!            zaj_l(iadd,ib,ik,1) = (fdexp - 1)*devr/wdi + evr
!            zaj_l(iadd,ib,ik,2) = (fdexp - 1)*devi/wdi + evi
             if(dabs(wdi) < SmallestPositiveNumber) then
                zaj_l(iadd,ib,ik,1) = -p_l(iadd,ib)*devr*dtim + evr
                zaj_l(iadd,ib,ik,2) = -p_l(iadd,ib)*devi*dtim + evi
             else
                fdexp = dexp( -p_l(iadd,ib) * wdi * dtim)
                zaj_l(iadd,ib,ik,1) = (fdexp - 1)*devr/wdi + evr
                zaj_l(iadd,ib,ik,2) = (fdexp - 1)*devi/wdi + evi
             end if
! ==============================================================================
          end do
       end do
    end if
                                                  __TIMER_DO_STOP(326)
!!$    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(309)
  end subroutine modified_steepest_descent_3D

!===============================================================================
  subroutine map_fft_to_WF_3D(ik,lsize, ibesize, wk_bfft_l, bfft_l, isrsize, fftsize)
    use m_Parallelization,     only : fft_wf_scnt, fft_wf_rcnt &
   &                                , fft_wf_recv &
   &                                , fft_wf_index, fft_wf_dist &
   &                                , fft_wf_maxrecv, fft_wf_maxsend

    integer, intent(in)  :: ik, lsize, ibesize, isrsize, fftsize
#ifdef FFT_3D_DIVISION
    real(kind=DP), dimension(lsize*2,ibesize), intent(in)  :: wk_bfft_l
    real(kind=DP), dimension(lsize*2,ibesize), intent(out) :: bfft_l
#else
    real(kind=DP), dimension(lsize*kimg,ibesize), intent(in)  :: wk_bfft_l
    real(kind=DP), dimension(lsize*kimg,ibesize), intent(out) :: bfft_l
#endif
    integer, dimension(0:nrank_g-1)                       ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)       ::sta_r, sta_s
    integer, parameter :: itag = 11

    real(kind=DP), allocatable, dimension(:,:) :: sendbuf, recvbuf
    integer :: icnt_send, icnt_recv, ierr, lrank, i, j, k, iadd
! === DEBUG by tkato 2012/06/05 ================================================
#ifdef USE_ALLTOALLV
    integer, allocatable, dimension(:) :: sdsp, rdsp
#endif
! ==============================================================================
                                                  __TIMER_SUB_START(308)

    if (fft_wf_maxsend(ik) /= 0) then
       allocate(sendbuf(fft_wf_maxsend(ik)*kimg*ibesize,0:nrank_g-1), stat=ierr)
!       sendbuf = 0.0d0
! === DEBUG by tkato 2012/06/05 ================================================
    else
       allocate(sendbuf(1,1), stat=ierr)
! ==============================================================================
    endif
    if (fft_wf_maxrecv(ik) /= 0) then
       allocate(recvbuf(fft_wf_maxrecv(ik)*kimg*ibesize,0:nrank_g-1), stat=ierr)
!       recvbuf = 0.0d0
! === DEBUG by tkato 2012/06/05 ================================================
    else
       allocate(recvbuf(1,1), stat=ierr)
! ==============================================================================
    endif
     if (ierr /= 0) then
        write(nfout,*)' map_info_fft_to_WF_3D :  Not allocate '
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 211, ierr)
     endif

#ifndef USE_ALLTOALLV
                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,323)
    if (fft_wf_maxrecv(ik) /= 0) then
       icnt_recv = 0
       lrank = mod(myrank_g,nrank_g)
       do i = 0, nrank_g - 1
          lrank = lrank + 1
          if (lrank > nrank_g -1) lrank = 0
          if (fft_wf_rcnt(lrank,ik) /= 0) then
             call mpi_irecv(recvbuf(1,lrank), fft_wf_rcnt(lrank,ik)*kimg*ibesize, &
            &               mpi_double_precision, lrank, itag, mpi_ke_world, req_r(icnt_recv), ierr)
              if (ierr /= 0) then
                 write(nfout,*)' map_info_fft_to_WF_3D :  mpi_irecv error'
                 call flush(nfout)
                 call mpi_abort(mpi_comm_world, 212, ierr)
              endif
             icnt_recv = icnt_recv + 1
          endif
       enddo
    endif

    if (fft_wf_maxsend(ik) /= 0) then
#endif
                                                  __TIMER_DO_START(324)
       if (kimg == 1) then
#ifdef FFT_3D_DIVISION
       integer :: i1, lx, ly, lz, mx, my, mz, mm, kx1p, kx2p, kx3p, jadd
         lx = fft_box_size_WF(1,0)
         ly = fft_box_size_WF(2,0)
         lz = fft_box_size_WF(3,0)
         kx1p = fft_X_x_nel
         kx2p = fft_X_y_nel
         kx3p = fft_X_z_nel
!OCL NORECURRENCE
          do k = 1, nel_fft_x(myrank_g)
             if(fft_wf_index(k,ik) == 0) cycle
             i1 = mp_fft_x(k)
             mz = (i1-1)/(lx*ly)+1
             mm = mod(i1,(lx*ly))
             if (mm==0) mm=lx*ly
             my = (mm-1)/lx+1
             mx = mod(mm,lx)
             if (mx==0) mx = lx
             jadd = mx-xyz_fft_x(1,1)+1+kx1p*(my-xyz_fft_x(1,2))+kx1p*kx2p*(mz-xyz_fft_x(1,3))
             do i = 1, ibesize
                sendbuf(ibesize*(fft_wf_index(k,ik)-1)+i,fft_wf_dist(k,ik)) = wk_bfft_l(jadd*2-1,i)
             enddo
          end do
#else
!OCL NORECURRENCE
          do k = 1, nel_fft_x(myrank_g)
             if(fft_wf_index(k,ik) == 0) cycle
             do i = 1, ibesize
                sendbuf(ibesize*(fft_wf_index(k,ik)-1)+i,fft_wf_dist(k,ik)) = wk_bfft_l(k,i)
             enddo
          end do
#endif
       else
#ifdef FFT_3D_DIVISION
         lx = fft_box_size_WF(1,0)
         ly = fft_box_size_WF(2,0)
         lz = fft_box_size_WF(3,0)
         kx1p = fft_X_x_nel
         kx2p = fft_X_y_nel
         kx3p = fft_X_z_nel
!OCL NORECURRENCE
          do k = 1, nel_fft_x(myrank_g)
             if(fft_wf_index(k,ik) == 0) cycle
             i1 = mp_fft_x(k)
             mz = (i1-1)/(lx*ly)+1
             mm = mod(i1,(lx*ly))
             if (mm==0) mm=lx*ly
             my = (mm-1)/lx+1
             mx = mod(mm,lx)
             if (mx==0) mx = lx
             jadd = mx-xyz_fft_x(1,1)+1+kx1p*(my-xyz_fft_x(1,2))+kx1p*kx2p*(mz-xyz_fft_x(1,3))
             do i = 1, ibesize
                iadd = ibesize*2*(fft_wf_index(k,ik)-1)+i*2
                sendbuf(iadd-1,fft_wf_dist(k,ik)) = wk_bfft_l(jadd*2-1,i)
                sendbuf(iadd,  fft_wf_dist(k,ik)) = wk_bfft_l(jadd*2  ,i)
             enddo
          end do
#else
!OCL NORECURRENCE
          do k = 1, nel_fft_x(myrank_g)
             if(fft_wf_index(k,ik) == 0) cycle
             do i = 1, ibesize
                iadd = ibesize*2*(fft_wf_index(k,ik)-1)+i*2
                sendbuf(iadd-1,fft_wf_dist(k,ik)) = wk_bfft_l(k*2-1,i)
                sendbuf(iadd,  fft_wf_dist(k,ik)) = wk_bfft_l(k*2  ,i)
             enddo
          end do
#endif
       end if
                                                  __TIMER_DO_STOP(324)

#ifndef USE_ALLTOALLV
       icnt_send = 0
       lrank = mod((myrank_g+1),nrank_g)
       do i = 0, nrank_g - 1
          lrank = lrank + 1
          if (lrank > (nrank_g - 1)) lrank = 0
          if (fft_wf_scnt(lrank,ik) /= 0) then
             call mpi_isend(sendbuf(1,lrank), fft_wf_scnt(lrank,ik)*kimg*ibesize, &
            &               mpi_double_precision, lrank, itag, mpi_ke_world, req_s(icnt_send), ierr)
              if (ierr /= 0) then
                 write(nfout,*)' map_info_fft_to_WF_3D :  mpi_isend error'
                 call flush(nfout)
                 call mpi_abort(mpi_comm_world, 213, ierr)
              endif
             icnt_send = icnt_send + 1
          endif
       enddo
    endif

    if (fft_wf_maxrecv(ik) /= 0) then
       call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
        if (ierr /= 0) then
           write(nfout,*)' map_info_fft_to_WF_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 214, ierr)
        endif
    endif

    if (fft_wf_maxsend(ik) /= 0) then
       call mpi_waitall(icnt_send, req_s, sta_s, ierr)
        if (ierr /= 0) then
           write(nfout,*)' map_info_fft_to_WF_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 215, ierr)
        endif
    endif
                                                  __TIMER_COMM_STOP(323)
#else
                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,330)
! === DEBUG by tkato 2012/06/05 ================================================
!   integer, allocatable, dimension(:) :: sdsp, rdsp
! ==============================================================================
    allocate(sdsp(0:nrank_g-1), stat=ierr)
    allocate(rdsp(0:nrank_g-1), stat=ierr)
    do i = 0, nrank_g - 1
       sdsp(i)=fft_wf_maxsend(ik)*kimg*ibesize*i
       rdsp(i)=fft_wf_maxrecv(ik)*kimg*ibesize*i
    enddo
    call MPI_ALLTOALLV(      sendbuf, fft_wf_scnt(:,ik)*kimg*ibesize, sdsp, &
   &   mpi_double_precision, recvbuf, fft_wf_rcnt(:,ik)*kimg*ibesize, rdsp, &
   &   mpi_double_precision, mpi_ke_world, ierr )
    if (ierr /= 0) then
       write(nfout,*)' map_info_fft_to_WF_3D :  mpi_alltoallv error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 216, ierr)
    endif
    deallocate(sdsp)
    deallocate(rdsp)
                                                  __TIMER_COMM_STOP(330)
#endif

!    bfft_l = 0.0d0
                                                  __TIMER_DO_START(325)
    if (kimg == 1) then
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          if (fft_wf_rcnt(i,ik) /= 0) then
             do k = 1, fft_wf_rcnt(i,ik)
                do j = 1, ibesize
                   bfft_l(fft_wf_recv(k,ik,i),j) = recvbuf(ibesize*(k-1)+j,i)
                enddo
             end do
          end if
       end do
    else
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          if (fft_wf_rcnt(i,ik) /= 0) then
             do k = 1, fft_wf_rcnt(i,ik)
                do j = 1, ibesize
                   iadd = ibesize*2*(k-1)+j*2
                   bfft_l(fft_wf_recv(k,ik,i)*2-1,j) = recvbuf(iadd-1,i)
                   bfft_l(fft_wf_recv(k,ik,i)*2  ,j) = recvbuf(iadd,  i)
                enddo
             end do
          end if
       end do
    end if
                                                  __TIMER_DO_STOP(325)

    if (allocated(sendbuf)) deallocate(sendbuf)
    if (allocated(recvbuf)) deallocate(recvbuf)
                                                  __TIMER_SUB_STOP(308)

  end subroutine map_fft_to_WF_3D

#ifdef MPI_FFTW
!===============================================================================
  subroutine gen_fft_to_WF_map(ik)
    use m_Parallelization,     only : fft_wf_scnt_mfftw, fft_wf_rcnt_mfftw &
   &                                , fft_wf_recv_mfftw &
   &                                , fft_wf_index_mfftw, fft_wf_dist_mfftw &
   &                                , fft_wf_maxrecv_mfftw, fft_wf_maxsend_mfftw
    integer, intent(in) :: ik
    integer :: nsize
    integer :: k, mm,i1,icount
    integer(C_INTPTR_T) :: lz,ly,lx,mx,my,mz,local_n, local_n_offset, alloc_local

    lx = fft_box_size_WF(1,0)
    ly = fft_box_size_WF(2,0)
    lz = fft_box_size_WF(3,0)
    alloc_local = fftw_mpi_local_size_3d(lz,ly,lx,mpi_ke_world,local_n,local_n_offset)
    if(allocated(mmx)) deallocate(mmx)
    if(allocated(mmy)) deallocate(mmy)
    if(allocated(mmz)) deallocate(mmz)
    if(allocated(iiadd)) deallocate(iiadd)
    if(allocated(wfdist)) deallocate(wfdist)
    nsize = local_n*lx*ly
    allocate(mmx(nsize))
    allocate(mmy(nsize))
    allocate(mmz(nsize))
    allocate(iiadd(nsize))
    allocate(wfdist(nsize))
    icount = 0
    if(kimg==2) then
      do k = 1, local_n*lx*ly
         if(fft_wf_index_mfftw(k,ik) == 0) cycle
         icount = icount+1
         i1 = k
         mz = (i1-1)/(lx*ly)+1
         mm = mod(i1,(lx*ly))
         if (mm==0) mm=lx*ly
         my = (mm-1)/lx+1
         mx = mod(mm,lx)
         if (mx==0) mx = lx
         iiadd(icount) = 2*(fft_wf_index_mfftw(k,ik))
         wfdist(icount) = fft_wf_dist_mfftw(k,ik)
         mmx(icount) = mx
         mmy(icount) = my
         mmz(icount) = mz
      end do
      maxwfdist = icount
    else
      do k = 1, local_n*lx*ly
         if(fft_wf_index_mfftw(k,ik) == 0) cycle
         icount = icount+1
         i1 = k
         mz = (i1-1)/(lx*ly)+1
         mm = mod(i1,lx*ly)
         if (mm==0) then
            my = ly
         else
            my = (mm-1)/lx+1
         end if
         mx = mod(mm,lx)
         if(mx==0) mx = lx
         iiadd(icount)  = fft_wf_index_mfftw(k,ik)
         wfdist(icount) = fft_wf_dist_mfftw(k,ik)
         mmx(icount) = mx
         mmy(icount) = my
         mmz(icount) = mz
      end do
      maxwfdist = icount
    endif
  end subroutine gen_fft_to_WF_map

  subroutine map_fft_to_WF_mpifftw(ik,lsize, ibesize, bfft_l, isrsize, fftsize)
    use m_Parallelization,     only : fft_wf_scnt_mfftw, fft_wf_rcnt_mfftw &
   &                                , fft_wf_recv_mfftw &
   &                                , fft_wf_index_mfftw, fft_wf_dist_mfftw &
   &                                , fft_wf_maxrecv_mfftw, fft_wf_maxsend_mfftw
    use m_FFT,                 only : bfft_mpifftw, bfft_mpifftw_kimg1

    integer, intent(in)  :: ik, lsize, ibesize, isrsize, fftsize
    real(kind=DP), dimension(lsize*kimg,ibesize), intent(out) :: bfft_l
    integer, dimension(0:nrank_g-1)                       ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)       ::sta_r, sta_s
    integer, parameter :: itag = 11

    real(kind=DP), allocatable, dimension(:,:) :: sendbuf, recvbuf
    integer :: icnt_send, icnt_recv, ierr, lrank, i, j, k, iadd, iwf
    integer :: i1,mm
    integer(C_INTPTR_T) :: lz,ly,lx,mx,my,mz,local_n, local_n_offset, alloc_local
    integer :: id_sname = -1

    call tstatc0_begin('map_fft_to_WF_mpifftw ', id_sname)
    lx = fft_box_size_WF(1,0)
    ly = fft_box_size_WF(2,0)
    lz = fft_box_size_WF(3,0)
    alloc_local = fftw_mpi_local_size_3d(lz,ly,lx,mpi_ke_world,local_n,local_n_offset)

    if (fft_wf_maxsend_mfftw(ik) /= 0) then
       allocate(sendbuf(fft_wf_maxsend_mfftw(ik)*kimg,0:nrank_g-1), stat=ierr)
!       sendbuf = 0.0d0
! === DEBUG by tkato 2012/06/05 ================================================
    else
       allocate(sendbuf(1,1), stat=ierr)
! ==============================================================================
    endif
    if (fft_wf_maxrecv_mfftw(ik) /= 0) then
       allocate(recvbuf(fft_wf_maxrecv_mfftw(ik)*kimg,0:nrank_g-1), stat=ierr)
!       recvbuf = 0.0d0
! === DEBUG by tkato 2012/06/05 ================================================
    else
       allocate(recvbuf(1,1), stat=ierr)
! ==============================================================================
    endif
     if (ierr /= 0) then
        write(nfout,*)' map_info_fft_to_WF_3D :  Not allocate '
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 211, ierr)
     endif

    if (fft_wf_maxrecv_mfftw(ik) /= 0) then
       icnt_recv = 0
       lrank = mod(myrank_g,nrank_g)
       do i = 0, nrank_g - 1
          lrank = lrank + 1
          if (lrank > nrank_g -1) lrank = 0
          if (fft_wf_rcnt_mfftw(lrank,ik) /= 0) then
             call mpi_irecv(recvbuf(1,lrank), fft_wf_rcnt_mfftw(lrank,ik)*kimg, &
            &               mpi_double_precision, lrank, itag, mpi_ke_world, req_r(icnt_recv), ierr)
              if (ierr /= 0) then
                 write(nfout,*)' map_info_fft_to_WF_3D :  mpi_irecv error'
                 call flush(nfout)
                 call mpi_abort(mpi_comm_world, 212, ierr)
              endif
             icnt_recv = icnt_recv + 1
          endif
       enddo
    endif

    if (fft_wf_maxsend_mfftw(ik) /= 0) then
       if (kimg == 1) then
!OCL NORECURRENCE

!          do k = 1, local_n*lx*ly
!             if(fft_wf_index_mfftw(k,ik) == 0) cycle
!             i1 = k
!             mz = (i1-1)/(lx*ly)+1
!             mm = mod(i1,lx*ly)
!             if (mm==0) then
!                my = ly
!             else
!                my = (mm-1)/lx+1
!             end if
!             mx = mod(mm,lx)
!             if(mx==0) mx = lx
!             sendbuf(fft_wf_index_mfftw(k,ik),fft_wf_dist_mfftw(k,ik)) = bfft_mpifftw_kimg1(mx,my,mz)
!          end do
          do k=1, maxwfdist
            iadd = iiadd(k)
            iwf  = wfdist(k)
            mx   = mmx(k)
            my   = mmy(k)
            mz   = mmz(k)
            sendbuf(iadd,iwf) = bfft_mpifftw_kimg1(mx,my,mz)
          enddo
       else
!OCL NORECURRENCE
!          do k = 1, local_n*lx*ly
!             if(fft_wf_index_mfftw(k,ik) == 0) cycle
!             i1 = k
!             mz = (i1-1)/(lx*ly)+1
!             mm = mod(i1,(lx*ly))
!             if (mm==0) mm=lx*ly
!             my = (mm-1)/lx+1
!             mx = mod(mm,lx)
!             if (mx==0) mx = lx
!             iadd = 2*(fft_wf_index_mfftw(k,ik))
!             sendbuf(iadd-1,fft_wf_dist_mfftw(k,ik)) = real (bfft_mpifftw(mx,my,mz))
!             sendbuf(iadd,  fft_wf_dist_mfftw(k,ik)) = aimag(bfft_mpifftw(mx,my,mz))
!          end do
          do k=1, maxwfdist
            iadd = iiadd(k)
            iwf  = wfdist(k)
            mx   = mmx(k)
            my   = mmy(k)
            mz   = mmz(k)
            sendbuf(iadd-1,iwf) = real (bfft_mpifftw(mx,my,mz))
            sendbuf(iadd,  iwf) = aimag(bfft_mpifftw(mx,my,mz))
          enddo
       end if

       icnt_send = 0
       lrank = mod((myrank_g+1),nrank_g)
       do i = 0, nrank_g - 1
          lrank = lrank + 1
          if (lrank > (nrank_g - 1)) lrank = 0
          if (fft_wf_scnt_mfftw(lrank,ik) /= 0) then
             call mpi_isend(sendbuf(1,lrank), fft_wf_scnt_mfftw(lrank,ik)*kimg, &
            &               mpi_double_precision, lrank, itag, mpi_ke_world, req_s(icnt_send), ierr)
              if (ierr /= 0) then
                 write(nfout,*)' map_info_fft_to_WF_3D :  mpi_isend error'
                 call flush(nfout)
                 call mpi_abort(mpi_comm_world, 213, ierr)
              endif
             icnt_send = icnt_send + 1
          endif
       enddo
    endif

    if (fft_wf_maxrecv_mfftw(ik) /= 0) then
       call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
        if (ierr /= 0) then
           write(nfout,*)' map_info_fft_to_WF_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 214, ierr)
        endif
    endif

    if (fft_wf_maxsend_mfftw(ik) /= 0) then
       call mpi_waitall(icnt_send, req_s, sta_s, ierr)
        if (ierr /= 0) then
           write(nfout,*)' map_info_fft_to_WF_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 215, ierr)
        endif
    endif

!    bfft_l = 0.0d0
    if (kimg == 1) then
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          if (fft_wf_rcnt_mfftw(i,ik) /= 0) then
             do k = 1, fft_wf_rcnt_mfftw(i,ik)
                bfft_l(fft_wf_recv_mfftw(k,ik,i),1) = recvbuf(k,i)
             end do
          end if
       end do
    else
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          if (fft_wf_rcnt_mfftw(i,ik) /= 0) then
             do k = 1, fft_wf_rcnt_mfftw(i,ik)
                iadd = 2*k
                bfft_l(fft_wf_recv_mfftw(k,ik,i)*2-1,1) = recvbuf(iadd-1,i)
                bfft_l(fft_wf_recv_mfftw(k,ik,i)*2  ,1) = recvbuf(iadd,  i)
             end do
          end if
       end do
    end if

    if (allocated(sendbuf)) deallocate(sendbuf)
    if (allocated(recvbuf)) deallocate(recvbuf)

    call tstatc0_end(id_sname)
  end subroutine map_fft_to_WF_mpifftw
#endif

  subroutine m_ESsd_sort_zaj_old_3D()
    real(kind=DP), allocatable, dimension(:,:,:)   :: wk_zaj,wk_ball0,wk_ball1
    real(kind=DP), allocatable, dimension(:,:,:,:) :: wk_mpi
    real(kind=DP), allocatable, dimension(:,:,:)     :: wk_ball00,wk_ball01
    integer,                    dimension(neg)     :: t_ordr, k_ordr, tr_neg
!!  integer(kind=4),            dimension(0:nrank_g-1) :: recvcnt, recvdsp
    integer               :: max_g1k, ib, jb, kb, ik, nb, i

    max_g1k = maxval(np_g1k(:))

    do ib = 1, neg
       tr_neg(neg_g_all(ib)) = ib
    end do

!!  do i = 0, nrank_g - 1
!!     recvcnt(i) = nel_e_3D(i) * max_g1k * kimg
!!  end do
!!  recvdsp(0) = 0
!!  do i = 1, nrank_g - 1
!!     recvdsp(i) = recvdsp(i-1) + recvcnt(i-1)
!!  end do

    do ik = 1, kv3, af+1
       if(map_k(ik) /= myrank_k) cycle    ! MPI

       allocate(wk_zaj(max_g1k,mp_e,kimg), stat=ierr)
       wk_zaj = 0.0d0
       allocate(wk_mpi(max_g1k,mp_e,kimg,0:nrank_e-1), stat=ierr)

       do kb = 1, kimg
          do jb = 1, np_e
             do ib = 1, np_g1k(ik)
                wk_zaj(ib,jb,kb) = zaj_old(ib,jb,ik,kb)
             end do
          end do
       end do

       call mpi_allgather(wk_zaj, max_g1k*mp_e*kimg, MPI_DOUBLE_PRECISION, &
      &                   wk_mpi, max_g1k*mp_e*kimg, MPI_DOUBLE_PRECISION, mpi_kg_world, ierr )
!!     call mpi_allgatherv(wk_zaj, max_g1k*np_e*kimg, MPI_DOUBLE_PRECISION, &
!!    &                    wk_mpi, recvcnt, recvdsp, MPI_DOUBLE_PRECISION, mpi_kg_world, ierr )

       deallocate(wk_zaj)
       allocate(wk_ball0(max_g1k,neg,kimg), stat=ierr)

       do nb = 0, nrank_e-1
          do kb = 1, kimg
             do jb = nis_e(nb), nie_e(nb)
                do ib = 1, np_g1k(ik)
                   wk_ball0(ib,neg_g_all(jb),kb) = wk_mpi(ib,jb-nis_e(nb)+1,kb,nb)
                end do
             end do
          end do
       end do

       deallocate(wk_mpi)
       allocate(wk_ball1(max_g1k,neg,kimg), stat=ierr)

       do kb = 1, kimg
          do jb = 1, neg
             do ib = 1, np_g1k(ik)
                wk_ball1(ib,tr_neg(jb),kb) = wk_ball0(ib,neordr_old(jb,ik),kb)
             end do
          end do
       end do

       do kb = 1, kimg
          do jb = ista_e, iend_e
             do ib = 1, np_g1k(ik)
                zaj_old(ib,jb-ista_e+1,ik,kb) = wk_ball1(ib,jb,kb)
             end do
          end do
       end do

       deallocate(wk_ball0)
       deallocate(wk_ball1)
    end do

  end subroutine m_ESsd_sort_zaj_old_3D

!===============================================================================
  subroutine decide_CG_direction_core_3D(precon,ik,ekin_l,afft_l,lsize,sumdz2,sumdz)
    integer, intent(in)        :: precon, ik, lsize
    real(kind=DP), intent(in)  :: ekin_l(np_g1k(ik))
#ifdef FFT_3D_DIVISION
    real(kind=DP), intent(in)  :: afft_l(lsize*2   )
#else
    real(kind=DP), intent(in)  :: afft_l(lsize*kimg)
#endif
    real(kind=DP), intent(out) :: sumdz2,sumdz

    real(kind=DP), allocatable, dimension(:,:) :: wk_bfft_l
    real(kind=DP), allocatable,dimension(:,:) :: VlocalW
    real(kind=DP), allocatable,dimension(:,:) :: p_l

    integer :: ibsize, flag, iesize
    integer :: ib1, ib2, ibesize
    integer :: isrsize, fft_l_size

    integer :: ib
    integer :: id_sname = -1
    call tstatc0_begin('decide_CG_direction_core ', id_sname,1)

    ibsize = 1
    if (nblocksize_fftw_is_given) then
       ibsize = nblocksize_fftw
       if (ibsize < 1) ibsize = 1
    endif
    isrsize = min(lsize,mp_g1k(ik))
    fft_l_size  = nel_fft_x(myrank_g)

#ifdef FFT_3D_DIVISION
    allocate(wk_bfft_l(lsize*2   ,ibsize) ,stat=ierr)
    allocate(VlocalW(lsize*2   ,ibsize) ,stat=ierr)
#else
    allocate(wk_bfft_l(lsize*kimg,ibsize) ,stat=ierr)
    allocate(VlocalW(lsize*kimg,ibsize) ,stat=ierr)
#endif
!P! allocate(p_l(mp_g1k(ik),ibsize) ,stat=ierr)
    allocate(p_l(mp_g1k(ik),np_e) ,stat=ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_ES_WF_in_Rspace_3D :  Not allocate '
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 203, ierr)
     endif
!P!
    call m_ES_decide_precon_factor_3D(precon,ik,1,np_e,np_e,ekin_l,p_l)! -> p(1:iba(ik))
!P!

    do ib1 = 1, np_e, ibsize           ! MPI
       ib2 = min(ib1+ibsize-1,np_e)
       ibesize = ib2 - ib1 + 1

#ifdef __TIMER_COMM__
       call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l, 0)
#else
       call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,wk_bfft_l)
#endif
#ifdef FFT_3D_DIVISION
       call m_FFT_Vlocal_W_3DIV_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_x(myrank_g))
       call m_FFT_Direct_3DIV_3D(nfout,wk_bfft_l,lsize,ibsize)
#else
       if (sw_fft_xzy > 0) then
          call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_y(myrank_g))
          call m_FFT_Direct_3D(nfout,wk_bfft_l,lsize,ibsize)
       else
          call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_z(myrank_g))
          if(sw_serial_fft == ON) then
             call m_ES_WF_2D(ik,wk_bfft_l,ib2,ib1,ibsize,lsize,DIRECT)
          else
             call m_FFT_Direct_XYZ_3D (nfout, wk_bfft_l, lsize, ibsize)
          endif
       end if
#endif
       call map_fft_to_WF_3D(ik, lsize, ibesize, wk_bfft_l, VlocalW, isrsize, fft_l_size)

       call SD_direction_3D(precon,ik,ib1,ib2,ibesize,ekin_l,VlocalW,lsize,p_l)
    end do

!fj$    call orthogonalize_SD_drctns(ik,to=OTHER_BANDS)  ! -(m_ES_WF_by_SDorCG) ->(wfsd_l, bsd(ri)_l)
    call orthogonalize_SD_drctns_3D(ik,to=OTHER_BANDS)  ! -(m_ES_WF_by_SDorCG) ->(wfsd_l, bsd(ri)_l)

    call WFSD_dot_WFSD                   ! -(contained here)
    !    ~~~~~~~~~~~~~

    deallocate(wk_bfft_l)
    deallocate(VlocalW)
    deallocate(p_l)

    call tstatc0_end(id_sname)
  contains
    subroutine WFSD_dot_WFSD
      real(kind=DP) :: dz,dz2
      integer :: ib

      sumdz2 = 0.d0
      sumdz = 0.0d0
      do ib = 1, np_e
         call square_of_SD_direction_3D(ik,ib,dz)
         dz2 = 0.d0
         if(sw_modified_cg_formula==ON) call square_of_SD_direction3_3D(ik,ib,dz2)
         sumdz2 = sumdz2 + dz
         sumdz = sumdz + (dz-dz2)
      end do
      if(npes > 1) then
        call mpi_allreduce(MPI_IN_PLACE,sumdz2,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
        call mpi_allreduce(MPI_IN_PLACE,sumdz ,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
      endif
    end subroutine WFSD_dot_WFSD
  end subroutine decide_CG_direction_core_3D
!===============================================================================
!fj$$#endif
!===============================================================================
! The reordering "wfsd_old" with "nrvf_ordr/neordr" is necessary
! to obtain the correct results.
! But, such a reordering needs much memory and expends overhead...
! by tkato 2011/08/19
!===============================================================================
  subroutine sort_wfsd_old_in(array, ik_)
  real(kind=DP), intent(out) :: array(maxval(np_g1k),np_e,ista_k:iend_k,kimg)
  integer, intent(in) :: ik_
  !!real(kind=DP) :: work(maxval(np_g1k),neg,kimg)
  real(kind=DP), allocatable, dimension(:,:,:) :: work
  integer :: ib, ncnt, ierr
  allocate(work(maxval(np_g1k),neg,kimg))
  ncnt = maxval(np_g1k)*neg*kimg
  work = 0.0d0
  do ib = 1, np_e
     work(:,nrvf_ordr(neg_g(ib),ik_),:) = array(:,ib,ik_,:)
  enddo
  call mpi_allreduce(MPI_IN_PLACE,work,ncnt,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
  do ib = 1, np_e
    array(:,ib,ik_,:) = work(:,neg_g(ib),:)
  enddo
  deallocate(work)
  return
  end subroutine sort_wfsd_old_in
!===============================================================================
  subroutine sort_wfsd_old_out(array, ik_)
  real(kind=DP), intent(out) :: array(maxval(np_g1k),np_e,ista_k:iend_k,kimg)
  integer, intent(in) :: ik_
  !!real(kind=DP) :: work(maxval(np_g1k),neg,kimg)
  real(kind=DP), allocatable, dimension(:,:,:) :: work
  integer :: ib, ncnt, ierr
  allocate(work(maxval(np_g1k),neg,kimg))
  ncnt = maxval(np_g1k)*neg*kimg
  work = 0.0d0
  do ib = 1, np_e
     work(:,neordr(neg_g(ib),ik_),:) = array(:,ib,ik_,:)
  enddo
  call mpi_allreduce(MPI_IN_PLACE,work,ncnt,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
  do ib = 1, np_e
    array(:,ib,ik_,:) = work(:,neg_g(ib),:)
  enddo
  deallocate(work)
  return
  end subroutine sort_wfsd_old_out
!===============================================================================
! === DEBUG by tkato 2014/03/20 ================================================
  subroutine m_ES_WF_by_SDorCG_epsmain_reallocate()
     if(allocated(zaj_old)) deallocate(zaj_old)
     allocate(zaj_old(maxval(np_g1k),np_e,ista_k:iend_k,kimg)); zaj_old = 0.0d0
  end subroutine m_ES_WF_by_SDorCG_epsmain_reallocate
! ==============================================================================

! =============
!   meta-gga
!
  subroutine m_ES_contrib_kindens_to_vnlph( is, ik, lsize, afft_l )
    use m_Crystal_Structure,  only : rltv
    use m_PlaneWaveBasisSet,   only : ngabc, kg
    use m_Electronic_Structure,   only : vtau_l, m_ES_WF_in_Rspace_3D1, &
         &                               m_ES_add_it_to_vnlph, &
         &                               m_ES_Vlocal_in_Rspace_3D

    use m_FFT,      only : nfft,fft_box_size_WF, m_FFT_Vlocal_W_3D, m_FFT_Direct_3D &
#ifdef FFT_3D_DIVISION
 &                        , m_FFT_Vlocal_W_3DIV_3D, m_FFT_Direct_3DIV_3D             &
#endif
 &                        , m_FFT_Direct_XYZ_3D

    integer, intent(in) :: is, ik, lsize
    real(kind=DP), intent(in) :: afft_l(lsize)

    integer :: ib, in, ri, i1, ig, ig2, ibsize, isrsize, fft_l_size
    real(kind=DP) :: denom
    real(kind=DP), allocatable :: G_vec(:,:), k_plus_G(:,:)
    real(kind=DP), allocatable :: bfft_l(:), wkbfft_l(:)
    real(kind=DP), allocatable :: zwk1(:,:,:,:), zwk2(:,:)

    if ( kimg == 1 ) stop "kimg = 1 does not work"

    ibsize = 1
    if (nblocksize_fftw_is_given) then
       ibsize = nblocksize_fftw
       if (ibsize < 1) ibsize = 1
    endif
    isrsize = min(lsize,mp_g1k(ik))
    fft_l_size  = nel_fft_x(myrank_g)

    allocate( k_plus_G(maxval(np_g1k),3) );

#ifdef FFT_3D_DIVISION
    allocate(bfft_l(lsize*2   ), stat=ierr)
    allocate(wkbfft_l(lsize*2   ), stat=ierr)
#else
    allocate(bfft_l(lsize*kimg), stat=ierr)
    allocate(wkbfft_l(lsize*kimg), stat=ierr)
#endif

    call k_plus_G_vectors_m_3D( ik, kg, kg1, kv3, iba, nbase, vkxyz, &
         &                      ngabc, rltv, &
         &                      k_plus_G(1,1), k_plus_G(1,2), k_plus_G(1,3) )

    denom = 1.d0 /product(fft_box_size_WF(1:3,1))

    allocate( zwk1(maxval(np_g1k),np_e,ik:ik,kimg) );  zwk1 = 0.0d0

    Do in=1, 3

       zwk1 = 0.0d0
       do ib = 1, np_e
          do ri=1,kimg
             Do ig=1, np_g1k(ik)
                zwk1(ig,ib,ik,ri) = zaj_l(ig,ib,ik,ri) *k_plus_G(ig,in)
             End Do
          End do
       end do

       do ib = 1, np_e
          allocate( zwk2(maxval(np_g1k),kimg) );  zwk2 = 0.0d0

          call m_ES_WF_in_Rspace_3D1( ik, ik, ik, ib, ib, ibsize, lsize, zwk1, wkbfft_l )
#ifdef FFT_3D_DIVISION
          call m_FFT_Vlocal_W_3DIV_3D(afft_l,wkbfft_l,lsize,ibsize,nel_fft_x(myrank_g))
          call m_FFT_Direct_3DIV_3D(nfout,wkbfft_l,lsize,ibsize)
#else
          if (sw_fft_xzy > 0) then
             call m_FFT_Vlocal_W_3D(afft_l,wkbfft_l,lsize,ibsize,nel_fft_y(myrank_g))
             call m_FFT_Direct_3D(nfout,wkbfft_l,lsize,ibsize)
          else
             call m_FFT_Vlocal_W_3D(afft_l,wkbfft_l,lsize,ibsize,nel_fft_z(myrank_g))
             if(sw_serial_fft == ON) then
!!                call m_ES_WF_2D(ik,bfft_l,ib2,ib1,ibsize,lsize,DIRECT)
                stop "serial fft does not work"
             else
                call m_FFT_Direct_XYZ_3D(nfout,wkbfft_l,lsize,ibsize)
             endif
          end if
#endif
          call map_fft_to_WF_3D( ik, lsize, ibsize, wkbfft_l, bfft_l, &
               &                 isrsize, fft_l_size )

          if ( kimg == 2 ) then
             do ig=1, np_g1k(ik)
                zwk2(ig,1) = zwk2(ig,1) +bfft_l(2*ig-1) *k_plus_G(ig,in) *denom
                zwk2(ig,2) = zwk2(ig,2) +bfft_l(2*ig  ) *k_plus_G(ig,in) *denom
             end do
          endif
          zwk2 = zwk2 /2.0d0
          call m_ES_add_it_to_vnlph( ik, ib, zwk2 )
          deallocate( zwk2 )
       end do
    end Do
    deallocate( zwk1 )
    deallocate( k_plus_G )
    deallocate( bfft_l );     deallocate( wkbfft_l )
!    stop

  end subroutine m_ES_contrib_kindens_to_vnlph

  subroutine m_ES_kindens_to_vnlph_ib( is, ik, ib, lsize, afft_l )
    use m_Crystal_Structure,  only : rltv
    use m_PlaneWaveBasisSet,   only : ngabc, kg
    use m_Electronic_Structure,   only : vtau_l, m_ES_WF_in_Rspace_kt_3D, &
         &                               m_ES_Vlocal_in_Rspace_3D

    use m_FFT,      only : nfft,fft_box_size_WF, m_FFT_Vlocal_W_3D, m_FFT_Direct_3D &
#ifdef FFT_3D_DIVISION
 &                        , m_FFT_Vlocal_W_3DIV_3D, m_FFT_Direct_3DIV_3D             &
#endif
 &                        , m_FFT_Direct_XYZ_3D

    integer, intent(in) :: is, ik, ib, lsize
    real(kind=DP), intent(in) :: afft_l(lsize)

    integer :: in, ri, i1, ig, ig2, ibsize, isrsize, fft_l_size
    real(kind=DP) :: denom
    real(kind=DP), allocatable :: G_vec(:,:), k_plus_G(:,:)
    real(kind=DP), allocatable :: bfft_l(:), wkbfft_l(:)
    real(kind=DP), allocatable :: zwk1(:,:,:,:), zwk2(:,:)

    if ( kimg == 1 ) stop "kimg = 1 does not work"

    ibsize = 1
    if (nblocksize_fftw_is_given) then
       ibsize = nblocksize_fftw
       if (ibsize < 1) ibsize = 1
    endif
    isrsize = min(lsize,mp_g1k(ik))
    fft_l_size  = nel_fft_x(myrank_g)

    allocate( k_plus_G( maxval(np_g1k),3 ) );

#ifdef FFT_3D_DIVISION
    allocate(bfft_l(lsize*2   ), stat=ierr)
    allocate(wkbfft_l(lsize*2   ), stat=ierr)
#else
    allocate(bfft_l(lsize*kimg), stat=ierr)
    allocate(wkbfft_l(lsize*kimg), stat=ierr)
#endif

    call k_plus_G_vectors_m_3D( ik, kg, kg1, kv3, iba, nbase, vkxyz, &
         &                      ngabc, rltv, &
         &                      k_plus_G(1,1), k_plus_G(1,2), k_plus_G(1,3) )

    denom = 1.d0 /product(fft_box_size_WF(1:3,1))

    allocate( zwk1(maxval(np_g1k),1,ik:ik,kimg) );  zwk1 = 0.0d0
    allocate( zwk2(maxval(np_g1k),kimg) );             zwk2 = 0.0d0
    Do in=1, 3
       do ri=1,kimg
          Do ig=1, np_g1k(ik)
             zwk1(ig,1,ik,ri) = zaj_l(ig,ib,ik,ri) *k_plus_G(ig,in)
         End Do
       End do
       call m_ES_WF_in_Rspace_kt_3D( ik ,ik, ik, ibsize, lsize, zwk1, wkbfft_l )

#ifdef FFT_3D_DIVISION
       call m_FFT_Vlocal_W_3DIV_3D(afft_l,wkbfft_l,lsize,ibsize,nel_fft_x(myrank_g))
       call m_FFT_Direct_3DIV_3D(nfout,wkbfft_l,lsize,ibsize)
#else
       if (sw_fft_xzy > 0) then
          call m_FFT_Vlocal_W_3D(afft_l,wkbfft_l,lsize,ibsize,nel_fft_y(myrank_g))
          call m_FFT_Direct_3D(nfout,wkbfft_l,lsize,ibsize)
       else
          call m_FFT_Vlocal_W_3D(afft_l,wkbfft_l,lsize,ibsize,nel_fft_z(myrank_g))
          if(sw_serial_fft == ON) then
             stop "serial_fft does not work"
          else
             call m_FFT_Direct_XYZ_3D(nfout,wkbfft_l,lsize,ibsize)
          endif
       end if
#endif
       call map_fft_to_WF_3D( ik, lsize, ibsize, wkbfft_l, bfft_l, &
            &                 isrsize, fft_l_size )

       if ( kimg == 2 ) then
          do ig=1, np_g1k(ik)
             zwk2(ig,1) = zwk2(ig,1) +bfft_l(2*ig-1) *k_plus_G(ig,in) *denom
             zwk2(ig,2) = zwk2(ig,2) +bfft_l(2*ig  ) *k_plus_G(ig,in) *denom
          end do
       endif
    End Do
    zwk2 = zwk2 /2.0d0       ! -> vnlph_l   zwk2( kg1, kimg )

    if(kimg==1) then
       do ig=1,np_g1k(ik)
          vnlph_l(ig,ib,1) = vnlph_l(ig,ib,1) + zwk2(ig,1)
       end do
    else
       do ig=1,np_g1k(ik)
          vnlph_l(ig,ib,1) = vnlph_l(ig,ib,1) + zwk2(ig,1)
          vnlph_l(ig,ib,2) = vnlph_l(ig,ib,2) + zwk2(ig,2)
       end do
    end if

    deallocate( zwk1 );   deallocate( zwk2 )
    deallocate( k_plus_G )
    deallocate( bfft_l );    deallocate( wkbfft_l )

  end subroutine m_ES_kindens_to_vnlph_ib

  subroutine m_ES_kindens_to_vnlph_ib2( is, ik, ib, lsize, afft_l, vtau_phl )
    use m_Crystal_Structure,  only : rltv
    use m_PlaneWaveBasisSet,   only : ngabc, kg
    use m_Electronic_Structure,   only : vtau_l, m_ES_WF_in_Rspace_kt_3D, &
         &                               m_ES_Vlocal_in_Rspace_3D

    use m_FFT,      only : nfft,fft_box_size_WF, m_FFT_Vlocal_W_3D, m_FFT_Direct_3D &
#ifdef FFT_3D_DIVISION
 &                        , m_FFT_Vlocal_W_3DIV_3D, m_FFT_Direct_3DIV_3D             &
#endif
 &                        , m_FFT_Direct_XYZ_3D

    integer, intent(in) :: is, ik, ib, lsize
    real(kind=DP), intent(out) :: vtau_phl(maxval(np_g1k),kimg)
    real(kind=DP), intent(in) :: afft_l(lsize)

    integer :: in, ri, i1, ig, ig2, ibsize, isrsize, fft_l_size
    real(kind=DP) :: denom
    real(kind=DP), allocatable :: G_vec(:,:), k_plus_G(:,:), vlength(:)
    real(kind=DP), allocatable :: bfft_l(:), wkbfft_l(:)
    real(kind=DP), allocatable :: zwk1(:,:,:,:), zwk2(:,:)

    if ( kimg == 1 ) stop "kimg = 1 does not work"

    ibsize = 1
    if (nblocksize_fftw_is_given) then
       ibsize = nblocksize_fftw
       if (ibsize < 1) ibsize = 1
    endif
    isrsize = min(lsize,mp_g1k(ik))
    fft_l_size  = nel_fft_x(myrank_g)

    allocate( k_plus_G(maxval(np_g1k),3) );

#ifdef FFT_3D_DIVISION
!    lsize = fft_X_x_nel*fft_X_y_nel*fft_X_z_nel
    allocate(bfft_l(lsize*2   ), stat=ierr)
    allocate(wkbfft_l(lsize*2   ), stat=ierr)
#else
!    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
    allocate(bfft_l(lsize*kimg), stat=ierr)
    allocate(wkbfft_l(lsize*kimg), stat=ierr)
#endif

    call k_plus_G_vectors_m_3D( ik, kg, kg1, kv3, iba, nbase, vkxyz, &
         &                      ngabc, rltv, &
         &                      k_plus_G(1,1), k_plus_G(1,2), k_plus_G(1,3) )

    denom = 1.d0 /product(fft_box_size_WF(1:3,1))

    allocate( zwk1(maxval(np_g1k),1,ik:ik,kimg) );  zwk1 = 0.0d0
    allocate( zwk2(maxval(np_g1k),kimg) );             zwk2 = 0.0d0
    Do in=1, 3
       do ri=1,kimg
          Do ig=1, np_g1k(ik)
             zwk1(ig,1,ik,ri) = zaj_l(ig,ib,ik,ri) *k_plus_G(ig,in)
         End Do
       End do
       call m_ES_WF_in_Rspace_kt_3D( ik ,ik, ik, ibsize, lsize, zwk1, wkbfft_l )

#ifdef FFT_3D_DIVISION
       call m_FFT_Vlocal_W_3DIV_3D(afft_l,wkbfft_l,lsize,ibsize,nel_fft_x(myrank_g))
       call m_FFT_Direct_3DIV_3D(nfout,wkbfft_l,lsize,ibsize)
#else
       if (sw_fft_xzy > 0) then
          call m_FFT_Vlocal_W_3D(afft_l,wkbfft_l,lsize,ibsize,nel_fft_y(myrank_g))
          call m_FFT_Direct_3D(nfout,wkbfft_l,lsize,ibsize)
       else
          call m_FFT_Vlocal_W_3D(afft_l,wkbfft_l,lsize,ibsize,nel_fft_z(myrank_g))
          if(sw_serial_fft == ON) then
             stop "serial_fft does not work"
          else
             call m_FFT_Direct_XYZ_3D(nfout,wkbfft_l,lsize,ibsize)
          endif
       end if
#endif
       call map_fft_to_WF_3D( ik, lsize, ibsize, wkbfft_l, bfft_l, &
            &                 isrsize, fft_l_size )

       if ( kimg == 2 ) then
          do ig=1, np_g1k(ik)
             zwk2(ig,1) = zwk2(ig,1) +bfft_l(2*ig-1) *k_plus_G(ig,in) *denom
             zwk2(ig,2) = zwk2(ig,2) +bfft_l(2*ig  ) *k_plus_G(ig,in) *denom
          end do
       endif
    End Do
    zwk2 = zwk2 /2.0d0       ! -> vnlph_l   zwk2( kg1, kimg

    if(kimg==1) then
       do ig=1,np_g1k(ik)
          vtau_phl(ig,1) = zwk2(ig,1)
       end do
    else
       do ig=1,np_g1k(ik)
          vtau_phl(ig,1) = zwk2(ig,1)
          vtau_phl(ig,2) = zwk2(ig,2)
       end do
    end if

    deallocate( zwk1 );   deallocate( zwk2 )
    deallocate( k_plus_G )
    deallocate( bfft_l );    deallocate( wkbfft_l )

  end subroutine m_ES_kindens_to_vnlph_ib2

#ifdef MPI_FFTW
  subroutine m_ES_con_kindens_to_vnlph_mpfw( is, ik, lsize, &
       &                         lx, ly, lz, local_n, afft_mpifftw_vlocal )
    use m_Crystal_Structure,  only : rltv
    use m_PlaneWaveBasisSet,   only : ngabc, kg
    use m_Electronic_Structure,   only : vtau_l, m_ES_WF_in_Rspace_3D1, &
         &                               m_ES_add_it_to_vnlph, &
         &                               m_ES_Vlocal_in_Rspace_3D, m_ES_fftbox_map

    use m_FFT,      only : nfft,fft_box_size_WF, m_FFT_Vlocal_W_3D, m_FFT_Direct_3D &
 &                        , m_FFT_Direct_XYZ_3D

    integer, intent(in) :: is, ik, lsize
    integer(C_INTPTR_T), intent(in) :: local_n, lx, ly, lz

    real(kind=DP), intent(in) :: afft_mpifftw_vlocal(lx,lz,local_n)

    integer :: ib, in, ri, i1, ig, ig2, ibsize, isrsize, fft_l_size, nsize
    real(kind=DP) :: denom
    real(kind=DP), allocatable :: G_vec(:,:), k_plus_G(:,:)
    real(kind=DP), allocatable :: bfft_l(:), wkbfft_l(:)
    real(kind=DP), allocatable :: zwk1(:,:,:,:), zwk2(:,:)

    if ( kimg == 1 ) stop "kimg = 1 does not work"

    ibsize = 1
    if (nblocksize_fftw_is_given) then
       ibsize = nblocksize_fftw
       if (ibsize < 1) ibsize = 1
    endif
    isrsize = min(lsize,mp_g1k(ik))
    fft_l_size  = lsize

    allocate( k_plus_G(maxval(np_g1k),3) );
    allocate(bfft_l(lsize*kimg), stat=ierr)

    call k_plus_G_vectors_m_3D( ik, kg, kg1, kv3, iba, nbase, vkxyz, &
         &                      ngabc, rltv, &
         &                      k_plus_G(1,1), k_plus_G(1,2), k_plus_G(1,3) )

    denom = 1.d0 /product(fft_box_size_WF(1:3,1))

    allocate( zwk1(maxval(np_g1k),np_e,ik:ik,kimg) );  zwk1 = 0.0d0

    call gen_fft_to_WF_map(ik)
    call m_ES_fftbox_map(ik)

    Do in=1, 3

       zwk1 = 0.0d0
       do ib = 1, np_e
          do ri=1,kimg
             Do ig=1, np_g1k(ik)
                zwk1(ig,ib,ik,ri) = zaj_l(ig,ib,ik,ri) *k_plus_G(ig,in)
             End Do
          End do
       end do

       do ib = 1, np_e
          allocate( zwk2(maxval(np_g1k),kimg) );  zwk2 = 0.0d0
          bfft_l = 0.0d0

          call m_ES_WF_in_Rspace_mpifftw( ik, ik, ik, ib, zwk1 )

          nsize = local_n*lx*ly
          call m_FFT_Vlocal_W_mpifftw3d(afft_mpifftw_vlocal,lx,local_n,lz)
          call m_FFT_Direct_MPI_FFTW(nfout)

          call map_fft_to_WF_mpifftw(ik,lsize,1,bfft_l,isrsize,fft_l_size)

          if ( kimg == 2 ) then
             do ig=1, np_g1k(ik)
                zwk2(ig,1) = zwk2(ig,1) +bfft_l(2*ig-1) *k_plus_G(ig,in) *denom
                zwk2(ig,2) = zwk2(ig,2) +bfft_l(2*ig  ) *k_plus_G(ig,in) *denom
             end do
          endif
          zwk2 = zwk2 /2.0d0
          call m_ES_add_it_to_vnlph( ik, ib, zwk2 )
          deallocate( zwk2 )
       end do
    end Do

    deallocate( zwk1 )
    deallocate( k_plus_G )
    deallocate( bfft_l )

  end subroutine m_ES_con_kindens_to_vnlph_mpfw

  subroutine m_ES_kindens_to_vnlph_ib_mpfw( is, ik, ib, lsize, &
       &                         lx, ly, lz, local_n, afft_mpifftw_vlocal )
    use m_Crystal_Structure,  only : rltv
    use m_PlaneWaveBasisSet,   only : ngabc, kg
    use m_Electronic_Structure,   only : vtau_l, m_ES_WF_in_Rspace_kt_3D, &
         &                               m_ES_Vlocal_in_Rspace_3D, m_ES_fftbox_map

    use m_FFT,      only : nfft,fft_box_size_WF, m_FFT_Vlocal_W_3D, m_FFT_Direct_3D &
 &                        , m_FFT_Direct_XYZ_3D

    integer, intent(in) :: is, ik, ib, lsize
    integer(C_INTPTR_T), intent(in) :: local_n, lx, ly, lz

    real(kind=DP), intent(in) :: afft_mpifftw_vlocal(lx,lz,local_n)

    integer :: in, ri, i1, ig, ig2, ibsize, isrsize, fft_l_size, nsize
    real(kind=DP) :: denom
    real(kind=DP), allocatable :: G_vec(:,:), k_plus_G(:,:)
    real(kind=DP), allocatable :: bfft_l(:), wkbfft_l(:)
    real(kind=DP), allocatable :: zwk1(:,:,:,:), zwk2(:,:)

    if ( kimg == 1 ) stop "kimg = 1 does not work"

    ibsize = 1
    if (nblocksize_fftw_is_given) then
       ibsize = nblocksize_fftw
       if (ibsize < 1) ibsize = 1
    endif
    isrsize = min(lsize,mp_g1k(ik))
    fft_l_size  = lsize

    allocate( k_plus_G( maxval(np_g1k),3 ) );
    allocate(bfft_l(lsize*kimg), stat=ierr)

    call k_plus_G_vectors_m_3D( ik, kg, kg1, kv3, iba, nbase, vkxyz, &
         &                      ngabc, rltv, &
         &                      k_plus_G(1,1), k_plus_G(1,2), k_plus_G(1,3) )

    denom = 1.d0 /product(fft_box_size_WF(1:3,1))

!    allocate( zwk1(maxval(np_g1k),1,ik:ik,kimg) );  zwk1 = 0.0d0
    allocate( zwk1(maxval(np_g1k),np_e,ik:ik,kimg) );  zwk1 = 0.0d0
    allocate( zwk2(maxval(np_g1k),kimg) );             zwk2 = 0.0d0

    call gen_fft_to_WF_map(ik)
    call m_ES_fftbox_map(ik)

    Do in=1, 3

       zwk1 = 0.0d0
       bfft_l = 0.0d0
       do ri=1,kimg
          Do ig=1, np_g1k(ik)
!             zwk1(ig,1,ik,ri) = zaj_l(ig,ib,ik,ri) *k_plus_G(ig,in)
             zwk1(ig,ib,ik,ri) = zaj_l(ig,ib,ik,ri) *k_plus_G(ig,in)
         End Do
       End do

!       call m_ES_WF_in_Rspace_mpifftw( ik, ik, ik, 1, zwk1 )
       call m_ES_WF_in_Rspace_mpifftw( ik, ik, ik, ib, zwk1 )

       nsize = local_n*lx*ly
       call m_FFT_Vlocal_W_mpifftw3d(afft_mpifftw_vlocal,lx,local_n,lz)
       call m_FFT_Direct_MPI_FFTW(nfout)

       call map_fft_to_WF_mpifftw(ik,lsize,1,bfft_l,isrsize,fft_l_size)

       if ( kimg == 2 ) then
          do ig=1, np_g1k(ik)
             zwk2(ig,1) = zwk2(ig,1) +bfft_l(2*ig-1) *k_plus_G(ig,in) *denom
             zwk2(ig,2) = zwk2(ig,2) +bfft_l(2*ig  ) *k_plus_G(ig,in) *denom
          end do
       endif
    End Do
    zwk2 = zwk2 /2.0d0       ! -> vnlph_l   zwk2( kg1, kimg )

    if(kimg==1) then
       do ig=1,np_g1k(ik)
!          vnlph_l(ig,ib,1) = vnlph_l(ig,ib,1) + zwk2(ig,1)
          vnlph_l(ig,ib,1) = zwk2(ig,1)
       end do
    else
       do ig=1,np_g1k(ik)
!          vnlph_l(ig,ib,1) = vnlph_l(ig,ib,1) + zwk2(ig,1)
!          vnlph_l(ig,ib,2) = vnlph_l(ig,ib,2) + zwk2(ig,2)
          vnlph_l(ig,ib,1) = zwk2(ig,1)
          vnlph_l(ig,ib,2) = zwk2(ig,2)
       end do
    end if

    deallocate( zwk1 );   deallocate( zwk2 )
    deallocate( k_plus_G )
    deallocate( bfft_l )

  end subroutine m_ES_kindens_to_vnlph_ib_mpfw

  subroutine m_ES_kindens_to_vnlph_ib2_mpfw( is, ik, ib, lsize, &
       &                         lx, ly, lz, local_n, afft_mpifftw_vlocal, vtau_phl )
    use m_Crystal_Structure,  only : rltv
    use m_PlaneWaveBasisSet,   only : ngabc, kg
    use m_Electronic_Structure,   only : vtau_l, m_ES_WF_in_Rspace_kt_3D, &
         &                               m_ES_Vlocal_in_Rspace_3D, m_ES_fftbox_map

    use m_FFT,      only : nfft,fft_box_size_WF, m_FFT_Vlocal_W_3D, m_FFT_Direct_3D &
 &                        , m_FFT_Direct_XYZ_3D

    integer, intent(in) :: is, ik, ib, lsize
    integer(C_INTPTR_T), intent(in) :: local_n, lx, ly, lz

    real(kind=DP), intent(in) :: afft_mpifftw_vlocal(lx,lz,local_n)
    real(kind=DP), intent(out) :: vtau_phl(maxval(np_g1k),kimg)

    integer :: in, ri, i1, ig, ig2, ibsize, isrsize, fft_l_size, nsize
    real(kind=DP) :: denom
    real(kind=DP), allocatable :: G_vec(:,:), k_plus_G(:,:), vlength(:)
    real(kind=DP), allocatable :: bfft_l(:), wkbfft_l(:)
    real(kind=DP), allocatable :: zwk1(:,:,:,:), zwk2(:,:)

    if ( kimg == 1 ) stop "kimg = 1 does not work"

    ibsize = 1
    if (nblocksize_fftw_is_given) then
       ibsize = nblocksize_fftw
       if (ibsize < 1) ibsize = 1
    endif
    isrsize = min(lsize,mp_g1k(ik))
    fft_l_size  = lsize

    allocate( k_plus_G(maxval(np_g1k),3) );
    allocate(bfft_l(lsize*kimg), stat=ierr)

    call gen_fft_to_WF_map(ik)
    call m_ES_fftbox_map(ik)

    call k_plus_G_vectors_m_3D( ik, kg, kg1, kv3, iba, nbase, vkxyz, &
         &                      ngabc, rltv, &
         &                      k_plus_G(1,1), k_plus_G(1,2), k_plus_G(1,3) )

    denom = 1.d0 /product(fft_box_size_WF(1:3,1))

!    allocate( zwk1(maxval(np_g1k),1,ik:ik,kimg) );  zwk1 = 0.0d0
    allocate( zwk1(maxval(np_g1k),np_e,ik:ik,kimg) );  zwk1 = 0.0d0
    allocate( zwk2(maxval(np_g1k),kimg) );             zwk2 = 0.0d0

    Do in=1, 3
       zwk1 = 0.0d0
       bfft_l = 0.0d0
       do ri=1,kimg
          Do ig=1, np_g1k(ik)
!             zwk1(ig,1,ik,ri) = zaj_l(ig,ib,ik,ri) *k_plus_G(ig,in)
             zwk1(ig,ib,ik,ri) = zaj_l(ig,ib,ik,ri) *k_plus_G(ig,in)
         End Do
       End do

!       call m_ES_WF_in_Rspace_mpifftw( ik, ik, ik, 1, zwk1 )
       call m_ES_WF_in_Rspace_mpifftw( ik, ik, ik, ib, zwk1 )
       nsize = local_n*lx*ly
       call m_FFT_Vlocal_W_mpifftw3d(afft_mpifftw_vlocal,lx,local_n,lz)
       call m_FFT_Direct_MPI_FFTW(nfout)

       call map_fft_to_WF_mpifftw(ik,lsize,1,bfft_l,isrsize,fft_l_size)

       if ( kimg == 2 ) then
          do ig=1, np_g1k(ik)
             zwk2(ig,1) = zwk2(ig,1) +bfft_l(2*ig-1) *k_plus_G(ig,in) *denom
             zwk2(ig,2) = zwk2(ig,2) +bfft_l(2*ig  ) *k_plus_G(ig,in) *denom
          end do
       endif
    End Do
    zwk2 = zwk2 /2.0d0

    if(kimg==1) then
       do ig=1,np_g1k(ik)
          vtau_phl(ig,1) = zwk2(ig,1)
       end do
    else
       do ig=1,np_g1k(ik)
          vtau_phl(ig,1) = zwk2(ig,1)
          vtau_phl(ig,2) = zwk2(ig,2)
       end do
    end if

    deallocate( zwk1 );   deallocate( zwk2 )
    deallocate( k_plus_G )
    deallocate( bfft_l )

  end subroutine m_ES_kindens_to_vnlph_ib2_mpfw

#endif

end module m_ES_WF_by_SDorCG
