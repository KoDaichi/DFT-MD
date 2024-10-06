!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 603 $)
!
!  MODULE: m_Electronic_Structure
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!  FURTHER MODIFICATION: T. Yamasaki, January/13/2004, April/15/2006, September/02/2008
!  FURTHER MODIFICATION: T. Yamasaki, T. Uda and T. Ohno, September 2009 (MGS_DGEMM)
!  FURTHER MODIFICATION: T. Yamasaki and T. Yamamoto,   October 2009  (NONLOCAL_DGEMM)
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

!   This module has been revised for the GAMMA point (k=(0,0,0)) by T. Yamasaki
!  in April 2006. Number of operations for the Gamma point have been tremendously
!  reduced in subroutines of m_ES_betar_dot_Wfs, m_ES_Vnonlocal_W, and
!  m_ES_modified_gram_schmidt.
!
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
#ifdef __TIMER_DGEMM__
#   define __TIMER_DGEMM_START(a)  call timer_sta(a)
#   define __TIMER_DGEMM_STOP(a)   call timer_end(a)
#else
#   define __TIMER_DGEMM_START(a)
#   define __TIMER_DGEMM_STOP(a)
#endif

!  2010.03
!  the performance of the Gram-Schmidt orthonormalization procedure
!  has been tuned through a joint effort with RIKEN.

#ifndef SX
#define DGEMM__       DGEMM
#endif

#ifndef NO_MGS_DGEMM
#define MGS_DGEMM
#endif

module m_ES_ortho
  use m_Control_Parameters, only : nspin,ipri,kimg, neg, printable, af &
       &                         , m_CtrlP_cachesize &
       &                         , sw_gep &
#ifdef SAVE_FFT_TIMES
       &                         , nblocksize_mgs, nblocksize_mgs_is_given, sw_save_fft
#else
       &                         , nblocksize_mgs, nblocksize_mgs_is_given
#endif
  use m_Const_Parameters,   only : DP, EXECUT, ON, OFF, DELTA, PAI2 &
       &                         , ORTHOGONALIZATION, ORTHONORMALIZATION, NORMALIZATION &
       &                         , NORMCONSERVATION, VDB, VANDERBILT_TYPE &
       &                         , OTHER_BANDS, SAME_BAND, GAMMA, OLD
  use m_Electronic_Structure,only: zaj_l, neordr, nrvf_ordr, fsr_l, fsi_l &
#ifdef SAVE_FFT_TIMES
 &                                , status_saved_phifftr &
#endif
       &                         , m_ES_wd_zaj_small_portion              &
       &                         , nblocksize_mgs_default
  use m_ES_nonlocal,        only : m_ES_betar_dot_WFs_4_each_k
  use m_Files,              only : nfout
  use m_Kpoints,            only : kv3,k_symmetry
  use m_Parallelization,    only : MPI_CommGroup,npes,mype &
       &                         , nrank_e,myrank_e,map_e,ista_e,iend_e,istep_e &
       &                         , map_z,np_e,mpi_k_world,myrank_k,ista_k,iend_k &
       &                         , ierr,mp_e,nis_e,nie_e,nel_e &
       &                         , ista_g1k,iend_g1k,np_g1k,mp_g1k,nis_g1k,nel_g1k &
       &                         , np_fs, mp_fs, ista_fs, iend_fs, nis_fs, nie_fs, nel_fs  &
       &                         , myrank_g, nrank_g &
       &                         , is_ffth, ista_ffth, nel_ffth, np_ffth, mp_ffth
  use m_PlaneWaveBasisSet,  only : kg1,kg,iba
  use m_PseudoPotential,    only : nlmt,nlmta,lmta,lmtt,ltp,q &
       &                         , modnrm,nac,fqwei,nlmta1,nlmta2 &
       &                         , nac_p, fqwei_p, nlmta1_p, nlmta2_p
  use m_Timing,             only : tstatc0_begin, tstatc0_end

! ================================ added by K. Tagami ================ 11.0
  use m_Control_Parameters,   only : noncol, ndim_spinor
  use m_PseudoPotential,      only : fqwei_noncl, fqwei_p_noncl
! ==================================================================== 11.0
  use mpi

  implicit none

  integer, private, parameter                         :: sw_timing_2ndlevel = ON

  integer                                             :: np_g1k_x     ! np_g1k_x = np_g1k(ik) | np_g1k(ik)+1
  real(kind=DP),private,allocatable,target,dimension(:,:,:)  :: psi_t  ! d(np_g1k_x,neg,kimg) a work_array for mgs
  real(kind=DP),private,allocatable,dimension(:,:)    :: p1Sp2
  integer                                             :: np_fs_x      ! np_fs_x = np_fs | np_fs+1
  real(kind=DP),private,allocatable,dimension(:,:)    :: bpr_t, bpi_t ! d(np_fs_x,neg) work arrays for mgs
#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
  real(kind=DP),private,allocatable,dimension(:)      :: psi_ii,psi_ir
  real(kind=DP),private,allocatable,dimension(:)      :: bp_ii,bp_ir
#endif
#ifdef MGS_DGEMM
!-$ tune FUJITSU for block ----->>
  real(kind=DP),private,allocatable,dimension(:,:)    :: bpr_tw1, bpr_tw2, bpi_tw1, bpi_tw2 ! d(nac_p,neg) work arrays for mgs
  integer,private                                     :: NB
#ifdef SX
  integer,private                                     :: MB, LMB
#endif
! <--
!-$ tune FUJITSU for block <<-----
#endif
! <--

!  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

! ============================= added by K. Tagami ================ 11.0
  real(kind=DP), private, allocatable, target :: psi_t_noncl(:,:,:,:)
  real(kind=DP), private, allocatable :: bpr_t_noncl(:,:,:)
  real(kind=DP), private, allocatable :: bpi_t_noncl(:,:,:)

#ifdef MGS_DGEMM
  real(kind=DP), private, allocatable :: bpr_tw1_noncl(:,:,:)
  real(kind=DP), private, allocatable :: bpr_tw2_noncl(:,:,:)
  real(kind=DP), private, allocatable :: bpi_tw1_noncl(:,:,:)
  real(kind=DP), private, allocatable :: bpi_tw2_noncl(:,:,:)
#endif
! ================================================================= 11.0

contains
!  1-18. m_ESortho_mgs_alloc
!        - mgs_vdb_alloc, - mgs_nrc_alloc
!  1-19. m_ESortho_mgs_dealloc
!        - mgs_vdb_dealloc, - mgs_nrc_dealloc
!  2-28. m_ES_modified_gram_schmidt
!  2-29. m_ES_orthogonalize_SD_to_WFs -> (2-30)
!  2-30. mgs_sd2wf_each_k_G           <- (2-29)  -> (2-34), (2-44)
!        - broadcast_fs, - Psi1SPhi2_t, - modify_bsd_and_phi_t, - alloc_phi_w_and_brd_w
!        - dealloc_phi_w_and_brd_w, - WSW, - normalize_bsd_and_phi, - Psi1SPhi2
!        - modify_bsd_and_phi
!  2-31. orthogonalize_SD             -> (2-33)
!  2-32. m_ES_MGS_4_each_k            -> (2-17), (2-33)
!        - wd_title_of_the_operation
!  2-33. mgs_4_each_k_G               <- (2-31), (2-32)  -> (2-34), (2-44)
!        - WSW_t, - normalize_bp_and_psi_t, - W1SW2_t_r
!        - modify_bp_and_psi_t_r, - substitute_jto_ib2back, - W1SW2_t
!        - modify_bp_and_psi_t, - alloc_and_brd_w, - dealloc_and_brd_w
!        - WSW, - normalize_bp_and_psi, - W1SW2, - modify_bp_and_psi
!  2-34. set_npzri                    <- (2-30), (2-33)
!  2-37. m_ES_W_transpose_back2
!  2-38. m_ES_W_transpose_r
!  2-39. m_ES_W_transpose_back_r
!  2-42. m_ES_F_transpose_r
!  2-43. m_ES_F_transpose_back_r

!  2-44. cp_bp_and_psi_2_brd_w2      <- (2-30), (2-33)

#ifdef VPP
#define _ODD_BOUNDARY_
#endif
#ifdef SX
#define _ODD_BOUNDARY_
#endif
#ifdef NEC_TUNE1
#define _ODD_BOUNDARY_
#endif


  subroutine m_ESortho_set_np_g1k_x(ik)
    integer, intent(in) :: ik
!   OUTPUT : np_g1k_x
!
#ifdef _ODD_BOUNDARY_
    if(mod(np_g1k(ik),2) == 0) then
       np_g1k_x = np_g1k(ik) + 1
    else
       np_g1k_x = np_g1k(ik)
    end if
#else
    np_g1k_x = np_g1k(ik) + 1
#endif
  end subroutine m_ESortho_set_np_g1k_x

  subroutine m_ESortho_set_np_fs_x()
#ifdef _ODD_BOUNDARY_
    if(mod(np_fs,2) == 0) then
       np_fs_x = np_fs+1
    else
       np_fs_x = np_fs
    end if
#else
    np_fs_x = np_fs + 1
#endif
  end subroutine m_ESortho_set_np_fs_x

  subroutine m_ESortho_mgs_alloc(ik)
    integer, intent(in) :: ik
    call m_ESortho_set_np_g1k_x(ik) ! -> np_g1k_x
!!$    call set_np_g1k_x() ! -> np_g1k_x

    if(modnrm == EXECUT) then
       call m_ESortho_set_np_fs_x() ! -> np_fs_x
!!$       call set_np_fs_x() ! -> np_fs_x
! ============================== modified by K. Tagami ============= 11.0
!       call mgs_vdb_alloc
       if ( noncol ) then
          call mgs_vdb_alloc_noncl(neg)
       else
          call mgs_vdb_alloc(neg)
       end if
! =================================================================== 11.0
    else
! ============================== modified by K. Tagami ============= 11.0
!       call mgs_nrc_alloc
       if ( noncol ) then
          call mgs_nrc_alloc_noncl(neg)
       else
          call mgs_nrc_alloc(neg)
       end if
! ===================================================================== 11.0
    end if
  contains
!!$    subroutine set_np_g1k_x()
!!$!   OUTPUT : np_g1k_x
!!$!
!!$!BRANCH_P ORG_Parallel
!!$#ifdef _ODD_BOUNDARY_
!!$      if(mod(np_g1k(ik),2) == 0) then
!!$         np_g1k_x = np_g1k(ik) + 1
!!$      else
!!$         np_g1k_x = np_g1k(ik)
!!$      end if
!!$#else
!!$!BRANCH_P_END ORG_Parallel
!!$      np_g1k_x = np_g1k(ik) + 1
!!$!BRANCH_P ORG_Parallel
!!$#endif
!!$!BRANCH_P_END ORG_Parallel
!!$    end subroutine set_np_g1k_x
!!$
!!$    subroutine set_np_fs_x()
!!$
!!$!BRANCH_P ORG_Parallel
!!$#ifdef _ODD_BOUNDARY_
!!$      if(mod(np_fs,2) == 0) then
!!$         np_fs_x = np_fs+1
!!$      else
!!$         np_fs_x = np_fs
!!$      end if
!!$#else
!!$!BRANCH_P_END ORG_Parallel
!!$      np_fs_x = np_fs + 1
!!$!BRANCH_P ORG_Parallel
!!$#endif
!!$!BRANCH_P_END ORG_Parallel
!!$    end subroutine set_np_fs_x

! ============================== added by K. Tagami ================== 11.0
    subroutine mgs_vdb_alloc_noncl(neg)
      integer, intent(in) ::neg
      integer :: kimg_t

      allocate(psi_t_noncl(np_g1k_x,neg,kimg,ndim_spinor))

      if((k_symmetry(ik) == GAMMA .and. kimg == 2) .or. kimg==1) then
         kimg_t = 1
      else
         kimg_t = 2
      end if
      allocate(p1Sp2(neg,kimg_t))

      if((k_symmetry(ik) == GAMMA .and. kimg == 2)) then
         kimg_t = 1
      else
         kimg_t = 2
      end if

      if(kimg_t == 1) then
         allocate(bpr_t_noncl(np_fs_x,neg,ndim_spinor))
#ifdef MGS_DGEMM
         allocate(bpr_tw1_noncl(nac_p,neg,ndim_spinor))
         allocate(bpr_tw2_noncl(nac_p,neg,ndim_spinor))
         bpr_tw1_noncl = 0.0
         bpr_tw2_noncl = 0.0
#endif
      else
         allocate(bpr_t_noncl(np_fs_x,neg,ndim_spinor))
         allocate(bpi_t_noncl(np_fs_x,neg,ndim_spinor))
#ifdef MGS_DGEMM
         allocate(bpr_tw1_noncl(nac_p,neg,ndim_spinor))
         allocate(bpr_tw2_noncl(nac_p,neg,ndim_spinor))
         allocate(bpi_tw1_noncl(nac_p,neg,ndim_spinor))
         allocate(bpi_tw2_noncl(nac_p,neg,ndim_spinor))
         bpr_tw1_noncl = 0.0
         bpr_tw2_noncl = 0.0
         bpi_tw1_noncl = 0.0
         bpi_tw2_noncl = 0.0
#endif
      end if
    end subroutine mgs_vdb_alloc_noncl
! ===================================================================== 11.0
    subroutine mgs_vdb_alloc(neg)
      integer,intent(in) :: neg
      integer :: kimg_t

      allocate(psi_t(np_g1k_x,neg,kimg))

#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
      allocate(psi_ir(np_g1k_x))
      if(kimg == 2) allocate(psi_ii(np_g1k_x))
      if(np_g1k(ik) < np_g1k_x) then
         psi_ir(np_g1k_x) = 0.d0
         if(kimg==2) psi_ii(np_g1k_x) = 0.d0
      end if
#endif

      if((k_symmetry(ik) == GAMMA .and. kimg == 2) .or. kimg==1) then
         kimg_t = 1
      else
         kimg_t = 2
      end if
      allocate(p1Sp2(neg,kimg_t))

      if((k_symmetry(ik) == GAMMA .and. kimg == 2)) then
         kimg_t = 1
      else
         kimg_t = 2
      end if


      if(kimg_t == 1) then
         allocate(bpr_t(np_fs_x,neg))
         allocate(bpi_t(1,1))
#ifdef MGS_DGEMM
! === Debug by Intel "-check all" option! by T.Kato 2013/07/02 =================
         if(nac_p > 0) then
!  =============================================================================
         allocate(bpr_tw1(nac_p,neg))
         allocate(bpr_tw2(nac_p,neg))
         allocate(bpi_tw1(nac_p,neg))
         allocate(bpi_tw2(nac_p,neg))
!         allocate(bpi_tw1(1,1))
!         allocate(bpi_tw2(1,1))
         bpr_tw1 = 0.0
         bpr_tw2 = 0.0
! === Debug by Intel "-check all" option! by T.Kato 2013/07/02 =================
         end if
!  =============================================================================
#endif
      else
         allocate(bpr_t(np_fs_x,neg))
         allocate(bpi_t(np_fs_x,neg))
#ifdef MGS_DGEMM
! === Debug by Intel "-check all" option! by T.Kato 2013/07/02 =================
         if(nac_p > 0) then
!  =============================================================================
         allocate(bpr_tw1(nac_p,neg))
         allocate(bpr_tw2(nac_p,neg))
         allocate(bpi_tw1(nac_p,neg))
         allocate(bpi_tw2(nac_p,neg))
         bpr_tw1 = 0.0
         bpr_tw2 = 0.0
         bpi_tw1 = 0.0
         bpi_tw2 = 0.0
! === Debug by Intel "-check all" option! by T.Kato 2013/07/02 =================
         end if
!  =============================================================================
#endif
      end if

#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
      allocate(bp_ir(max(np_fs,np_g1k_x)));bp_ir=0.d0
      if(kimg_t == 2) then
         allocate(bp_ii(max(np_fs,np_g1k_x)));bp_ii=0.d0
      endif
!!$      allocate(bp_ii(np_g1k_x))
      if(np_g1k(ik) < np_g1k_x) then
         bp_ir(np_g1k(ik)+1:np_g1k_x) = 0.d0
!!$         if(kimg==2) bp_ii(np_g1k_x) = 0.d0
         if(kimg_t == 2) bp_ii(np_g1k(ik)+1:np_g1k_x) = 0.d0
!!$         bp_ii(np_g1k(ik)+1:np_g1k_x) = 0.d0
      end if
#endif
    end subroutine mgs_vdb_alloc


! ============================== added by K. Tagami ================== 11.0
    subroutine mgs_nrc_alloc_noncl(neg)
      integer,intent(in) :: neg
      integer :: kimg_t
      allocate(psi_t_noncl(np_g1k_x,neg,kimg,ndim_spinor))

      if((k_symmetry(ik) == GAMMA .and. kimg == 2).or.kimg==1) then
         kimg_t = 1
      else
         kimg_t = 2
      end if
      allocate(p1Sp2(neg,kimg_t))
    end subroutine mgs_nrc_alloc_noncl
! ================================================================= 11.0
    subroutine mgs_nrc_alloc(neg)
      integer,intent(in) :: neg
      integer :: kimg_t
      allocate(psi_t(np_g1k_x,neg,kimg))
#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
      allocate(psi_ir(np_g1k_x))
      if(kimg==2) allocate(psi_ii(np_g1k_x))
      if(np_g1k(ik) < np_g1k_x) then
         psi_ir(np_g1k_x) = 0.d0
         if(kimg==2) psi_ii(np_g1k_x) = 0.d0
      end if
#endif

      if((k_symmetry(ik) == GAMMA .and. kimg == 2).or.kimg==1) then
         kimg_t = 1
      else
         kimg_t = 2
      end if
      allocate(p1Sp2(neg,kimg_t))
    end subroutine mgs_nrc_alloc
  end subroutine m_ESortho_mgs_alloc


  subroutine m_ESortho_mgs_dealloc()
    if(modnrm == EXECUT) then
! ================================ modified by K. Tagami ========== 11.0
!!       call mgs_vdb_dealloc()
!
       if ( noncol ) then
         call mgs_vdb_dealloc_noncl()
       else
         call mgs_vdb_dealloc()
       endif
! ================================================================= 11.0
    else
! ================================ modified by K. Tagami ========== 11.0
!       call mgs_nrc_dealloc()

       if ( noncol ) then
          call mgs_nrc_dealloc_noncl()
       else
          call mgs_nrc_dealloc()
       endif
! ==================================================================== 11.0
    end if
  contains
    subroutine mgs_vdb_dealloc()
      deallocate(bpr_t)
#ifdef MGS_DGEMM
! === Debug by Intel "-check all" option! by T.Kato 2013/07/02 =================
!     deallocate(bpr_tw1)
!     deallocate(bpr_tw2)
      if(allocated(bpr_tw1)) deallocate(bpr_tw1)
      if(allocated(bpr_tw2)) deallocate(bpr_tw2)
!  =============================================================================
#endif
      if(allocated(bpi_t)) deallocate(bpi_t)
#ifdef MGS_DGEMM
      if(allocated(bpi_tw1)) deallocate(bpi_tw1)
      if(allocated(bpi_tw2)) deallocate(bpi_tw2)
#endif
#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
      deallocate(psi_ir)
      if(kimg == 2) deallocate(psi_ii)
#endif
      deallocate(psi_t)
      deallocate(p1Sp2)
#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
      if(allocated(bp_ir)) deallocate(bp_ir)
      if(allocated(bp_ii)) deallocate(bp_ii)
#endif
    end subroutine mgs_vdb_dealloc

! ============================= added by K. Tagami ================= 11.0
    subroutine mgs_vdb_dealloc_noncl()
      deallocate(bpr_t_noncl)
#ifdef MGS_DGEMM
      deallocate(bpr_tw1_noncl)
      deallocate(bpr_tw2_noncl)
#endif
      if(allocated(bpi_t_noncl)) deallocate(bpi_t_noncl)
#ifdef MGS_DGEMM
      if(allocated(bpi_tw1_noncl)) deallocate(bpi_tw1_noncl)
      if(allocated(bpi_tw2_noncl)) deallocate(bpi_tw2_noncl)
#endif
      deallocate(psi_t_noncl)
      deallocate(p1Sp2)
!!$#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
!!$      if(mod_pot == VANDERBILT_TYPE) then
!!$         if(allocated(bp_ir)) deallocate(bp_ir)
!!$         if(allocated(bp_ii)) deallocate(bp_ii)
!!$      end if
!!$#endif
    end subroutine mgs_vdb_dealloc_noncl
! ================================================================= 11.0

    subroutine mgs_nrc_dealloc()
#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
      deallocate(psi_ir)
      if(kimg == 2) deallocate(psi_ii)
#endif
      deallocate(psi_t)
      deallocate(p1Sp2)
    end subroutine mgs_nrc_dealloc

! ============================== added by K. Tagami =============== 11.0
    subroutine mgs_nrc_dealloc_noncl()
      deallocate(psi_t_noncl)
      deallocate(p1Sp2)
    end subroutine mgs_nrc_dealloc_noncl
! ================================================================== 11.0

  end subroutine m_ESortho_mgs_dealloc

  subroutine m_ES_modified_gram_schmidt(nfout)
  use m_IterationNumbers,     only : iteration
    integer, intent(in) :: nfout
    integer :: ik

    do ik = ista_k, iend_k, af+1                            ! MPI
       if(ipri>=2) call m_ES_wd_zaj_small_portion(nfout,ik," -- before GS --",16)

       if(sw_gep == OFF .or. iteration<2)then
          call m_ES_MGS_4_each_k(nfout,ik,mode=ORTHONORMALIZATION)
       else
          call m_ES_MGS_4_each_k(nfout,ik,mode=NORMALIZATION)
       endif
       !    ~~~~~~~~~~~~~~~~~~~~~~
       if(ipri>=2) call m_ES_wd_zaj_small_portion(nfout,ik," -- after GS --",15)

    end do

  end subroutine m_ES_modified_gram_schmidt

! ====================================== added by K. Tagami ============= 11.0
  subroutine m_ES_modified_gramschmidt_noncl(nfout)
    integer, intent(in) :: nfout

    integer :: ik, is

    do ik = ista_k, iend_k, ndim_spinor
      do is=1, ndim_spinor
        if(ipri>=2) call m_ES_wd_zaj_small_portion( &
                  nfout, ik+is-1,  " -- before GS --",16)
      end do
      call m_ES_MGS_4_each_k_noncl( nfout,ik,mode=ORTHONORMALIZATION )
      do is=1, ndim_spinor
        if(ipri>=2) call m_ES_wd_zaj_small_portion( &
     &              nfout, ik+is-1,  " -- after GS --",15)
      end do
   end do
  end subroutine m_ES_modified_gramschmidt_noncl
! ========================================================================11.0

!!$!BRANCH_P 3D_Parallel
!!$  subroutine m_ES_modified_gram_schmidt_3D(nfout)
!!$    integer, intent(in) :: nfout
!!$    integer :: ik
!!$
!!$    do ik = ista_k, iend_k, af+1                            ! MPI
!!$       if(ipri>=2) call m_ES_wd_zaj_small_portion_3D(nfout,ik," -- before GS --",16)
!!$
!!$       call m_ES_MGS_4_each_k(nfout,ik,mode=ORTHONORMALIZATION)
!!$       !    ~~~~~~~~~~~~~~~~~~~~~~
!!$       if(ipri>=2) call m_ES_wd_zaj_small_portion_3D(nfout,ik," -- after GS --",15)
!!$    end do
!!$
!!$  end subroutine m_ES_modified_gram_schmidt_3D
!!$!BRANCH_P_END 3D_Parallel

  subroutine m_ES_orthogonal_phi_to_WFs(ik,wfsd_l,bsdr_l,bsdi_l)
    integer, intent(in)                                     :: ik
    real(kind=DP),intent(inout),dimension(kg1,np_e,ik:ik,kimg) :: wfsd_l
    real(kind=DP),intent(inout),dimension(np_e,nlmta,ik:ik)    :: bsdr_l,bsdi_l

    integer  :: id_sname = -1, id_sname2 = -1

    call tstatc0_begin('m_ES_orthogonal_phi_to_WFs ',id_sname,1)

    call tstatc0_begin('mpi_barrier(ortho_phi_to_WFs) ',id_sname2)
    call mpi_barrier(mpi_k_world(myrank_k),ierr)
    call tstatc0_end(id_sname2)

    if(modnrm == EXECUT) then
       call mgs_phi2wf_each_k_G(ik,wfsd_l,ORTHOGONALIZATION,bsdr_l,bsdi_l,mod_pot=VDB)
    else
       call mgs_phi2wf_each_k_G(ik,wfsd_l,ORTHOGONALIZATION,mod_pot=NORMCONSERVATION)
    end if

    call tstatc0_end(id_sname)
  end subroutine m_ES_orthogonal_phi_to_WFs

  subroutine m_ES_orthogonalize_SD_to_WFs(ik,to_which_band,wfsd_l,bsdr_l,bsdi_l)
    integer, intent(in)                                     :: ik,to_which_band
    real(kind=DP),intent(inout),dimension(kg1,np_e,ik:ik,kimg) :: wfsd_l
    real(kind=DP),intent(inout),dimension(np_e,nlmta,ik:ik)    :: bsdr_l,bsdi_l

    integer  :: id_sname = -1
    call tstatc0_begin('m_ES_orthogonalize_SD_to_WFs ',id_sname)

!!$    call phi_alloc(ik)
!!$    call m_ES_mgs_alloc(ik)
    if(modnrm == EXECUT) then
       call mgs_sd2wf_each_k_G(ik,to_which_band,wfsd_l,bsdr_l,bsdi_l,mod_pot=VDB)
    else
       call mgs_sd2wf_each_k_G(ik,to_which_band,wfsd_l,mod_pot=NORMCONSERVATION)
    end if
    ! -(m_E.S.)   ->wfsd_l,bsdr_l,bsdi_l!!$    if(ik == iend_k) call phi_dealloc()    !-(m_E.S.)
!!$    call m_ES_mgs_dealloc

!!$    if(ik == iend_k) call phi_dealloc()    !-(m_E.S.)

    call tstatc0_end(id_sname)
  end subroutine m_ES_orthogonalize_SD_to_WFs

! ========================= added by K. Tagami ====================== 11.0
  subroutine m_ES_orthogonl_SD_to_WFs_noncl( ik, to_which_band, ikst, iken, &
	&                                    wfsd_l, bsdr_l, bsdi_l)
    integer, intent(in)              :: ik,to_which_band
    integer, intent(in) :: ikst, iken

    real(kind=DP),intent(inout),dimension(kg1,np_e,ikst:iken,kimg) :: wfsd_l
    real(kind=DP),intent(inout),dimension(np_e,nlmta,ikst:iken)    :: bsdr_l
    real(kind=DP),intent(inout),dimension(np_e,nlmta,ikst:iken)    :: bsdi_l

    integer  :: id_sname = -1
    call tstatc0_begin('m_ES_orthogonl_SD_to_WFs_noncl ',id_sname)

    if(modnrm == EXECUT) then
       call mgs_sd2wf_each_k_G_noncl( ik, to_which_band, wfsd_l, ikst, iken, &
	&                             bsdr_l, bsdi_l, mod_pot=VDB )
    else
       call mgs_sd2wf_each_k_G_noncl( ik, to_which_band, wfsd_l, ikst, iken, &
	&                             mod_pot=NORMCONSERVATION )
    end if
    call tstatc0_end(id_sname)

  end subroutine m_ES_orthogonl_SD_to_WFs_noncl

! ================================================================== 11.0

  subroutine mgs_sd2wf_each_k_G(ik,to,phi_l,bsdr_l,bsdi_l,mod_pot)
    integer, intent(in)             :: ik,to
    real(kind=DP), intent(inout)    :: phi_l(kg1,np_e,ik:ik,kimg)
    real(kind=DP), optional,intent(inout), dimension(np_e,nlmta,ik:ik) :: bsdr_l,bsdi_l
    integer, intent(in)             :: mod_pot
    integer               :: i

    integer :: kimg_t
    real(kind=DP),allocatable,dimension(:,:)    :: bp_w
    real(kind=DP),allocatable, dimension(:,:,:) :: phi_t

!!  call m_ES_mgs_alloc(ik)
    call m_ESortho_mgs_alloc(ik)
    allocate(phi_t(np_g1k_x,neg,kimg))
    if(k_symmetry(ik) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = 2
    end if
    allocate(bp_w(nlmta,kimg_t))

    call m_ES_W_transpose_r(.false.,ista_k,iend_k,ik,zaj_l,psi_t) ! zaj_l(ig,ie,ik,ri) -> psi_t(ig,ie,ri) (transp.)
    call m_ES_W_transpose_r(.false.,ik,ik,ik,phi_l,phi_t) ! phi_l(ig,ie,ik,ri) -> phi_t(ig,ie,ri) (transposed)

    do i = 1, neg
       if(mod_pot == VDB) call broadcast_fs(i)  ! (fs(ri)_l(i,*,ik)-> bp_w
       call Psi1SPhi2_t(ik,i)                      ! -(contained here) -> p1Sp2
       call modify_bsd_and_phi_t(i)                ! -(c.h.) p1Sp2,phi_t -> bsd(ri),phi_t
    end do
    call m_ES_W_transpose_back_r(.false.,ik,ik,ik,phi_l,phi_t)       ! phi_t -> phi_l
    deallocate(phi_t)
    deallocate(bp_w)
!xx call m_ES_mgs_dealloc
    call m_ESortho_mgs_dealloc
  contains
    subroutine broadcast_fs(io)
      integer, intent(in) :: io
      integer  :: ia, i, kimg_t

      i = map_z(io)                                 ! MPI
      if(k_symmetry(ik) == GAMMA) then
         kimg_t = 1
      else
         kimg_t = 2
      end if

! ==== DEBUG by tkato 2011/09/01 ===============================================
!     if(map_e(i) == myrank_e) then                 ! MPI
      if(map_e(io) == myrank_e) then                 ! MPI
! ==============================================================================
         if(k_symmetry(ik) == GAMMA) then
            do ia = 1, nlmta
               bp_w(ia,1) = fsr_l(i,ia,ik)
            end do
         else
            do ia = 1, nlmta
               bp_w(ia,1) = fsr_l(i,ia,ik)
               bp_w(ia,2) = fsi_l(i,ia,ik)
            end do
         end if
      end if                                        ! MPI
! ==== DEBUG by tkato 2011/09/01 ===============================================
!     call mpi_bcast(bp_w,nlmta*kimg_t,mpi_double_precision,map_e(i) &
      call mpi_bcast(bp_w,nlmta*kimg_t,mpi_double_precision,map_e(io) &
! ==============================================================================
           & ,mpi_k_world(myrank_k),ierr)           ! MPI
    end subroutine broadcast_fs

    subroutine Psi1SPhi2_t(ik,i)
      integer, intent(in) :: ik,i

      real(kind=DP),allocatable,dimension(:,:)  :: p1Sp2_w, p1Sp2_w2

      integer   :: j,ia,p,q, kimg_t
      real(DP)  :: ar,ai

!!    DEBUG asms
!!    if((k_symmetry(1) == GAMMA .and. kimg == 2).or.kimg==1) then
      if((k_symmetry(ik) == GAMMA .and. kimg == 2).or.kimg==1) then
         kimg_t = 1
      else
         kimg_t = 2
      end if

      if(npes >= 2) then
         if(mod_pot == VANDERBILT_TYPE) then
            allocate(p1Sp2_w(neg,kimg_t))
         end if
         allocate(p1Sp2_w2(neg,kimg_t))
      end if

      p1Sp2 = 0.d0
      if(mod_pot == VANDERBILT_TYPE) then
         do j = 1, np_e                             ! MPI
            if(to == OTHER_BANDS .and. j+ista_e-1 == i) cycle               ! MPI
            if(to == SAME_BAND .and. j+ista_e-1 /= i) cycle
            ar = 0.d0; if(kimg == 2) ai = 0.d0
            if(kimg == 1) then
               do ia = 1, nac
                  p = nlmta1(ia);      q = nlmta2(ia)
                  ar = ar + fqwei(ia)*(bp_w(p,1)*bsdr_l(j,q,ik) + bp_w(p,2)*bsdi_l(j,q,ik))
               end do
! ==== DEBUG by tkato 2011/09/01 ===============================================
!              p1Sp2(j,1) = ar
               p1Sp2(j+ista_e-1,1) = ar
! ==============================================================================
            else if(kimg == 2) then
               if(k_symmetry(ik) == GAMMA) then
                  do ia = 1, nac
                     p = nlmta1(ia);      q = nlmta2(ia)
                     ar = ar + fqwei(ia)*(bp_w(p,1)*bsdr_l(j,q,ik))
                  end do
! ==== DEBUG by tkato 2011/09/01 ===============================================
!                 p1Sp2(j,1) = ar
                  p1Sp2(j+ista_e-1,1) = ar
! ==============================================================================
               else
                  do ia = 1, nac
                     p = nlmta1(ia);      q = nlmta2(ia)
                     ar = ar + fqwei(ia)*(bp_w(p,1)*bsdr_l(j,q,ik) + bp_w(p,2)*bsdi_l(j,q,ik))
                     ai = ai + fqwei(ia)*(bp_w(p,1)*bsdi_l(j,q,ik) - bp_w(p,2)*bsdr_l(j,q,ik))
                  end do
! ==== DEBUG by tkato 2011/09/01 ===============================================
!                 p1Sp2(j,1) = ar; p1Sp2(j,2) = ai
                  p1Sp2(j+ista_e-1,1) = ar; p1Sp2(j+ista_e-1,2) = ai
! ==============================================================================
               end if
            end if
         end do
         if(npes >= 2) then
            call mpi_allreduce(p1Sp2,p1Sp2_w,neg*kimg_t,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)  ! MPI
            p1Sp2 = 0.d0
         end if
      end if

      do j = 1, neg
         if(to == OTHER_BANDS .and. j == i) cycle
         if(to == SAME_BAND .and. j /= i) cycle
         if(kimg == 1) then
            do ia = 1, np_g1k(ik)                                           ! MPI
               p1Sp2(j,1) = p1Sp2(j,1) + psi_t(ia,i,1)*phi_t(ia,j,1)
            end do
         else if(kimg == 2) then
            if(k_symmetry(ik) == GAMMA) then
               do ia = 2, np_g1k(ik)                                        ! MPI
                  ar = psi_t(ia,i,1)
                  ai = psi_t(ia,i,2)
                  p1Sp2(j,1) = p1Sp2(j,1) + (ar*phi_t(ia,j,1) + ai*phi_t(ia,j,2))*2.d0
               end do
!!$               if(mype /= 0) then
               if(myrank_e /= 0) then
                  ar = psi_t(1,i,1)
                  ai = psi_t(1,i,2)
                  p1Sp2(j,1) = p1Sp2(j,1) + (ar*phi_t(1,j,1) + ai*phi_t(1,j,2))*2.d0
               else
! ==== DEBUG by tkato 2011/09/03 ===============================================
!                 p1Sp2(j,1) = p1Sp2(j,1) + psi_t(1,i,1)*psi_t(1,j,1)
                  p1Sp2(j,1) = p1Sp2(j,1) + psi_t(1,i,1)*phi_t(1,j,1)
! ==============================================================================
               end if
            else
               do ia = 1, np_g1k(ik)                                         ! MPI
                  ar = psi_t(ia,i,1)
                  ai = psi_t(ia,i,2)
                  p1Sp2(j,1) = p1Sp2(j,1) + ar*phi_t(ia,j,1) + ai*phi_t(ia,j,2)
                  p1Sp2(j,2) = p1Sp2(j,2) + ar*phi_t(ia,j,2) - ai*phi_t(ia,j,1)
               end do
            end if
         end if
      end do
      if(npes >= 2) then
         call mpi_allreduce(p1Sp2,p1Sp2_w2,neg*kimg_t,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)  ! MPI
         if(mod_pot == VANDERBILT_TYPE) then                                ! MPI
            p1Sp2 = p1Sp2_w + p1Sp2_w2                                      ! MPI
         else                                                               ! MPI
            p1Sp2 = p1Sp2_w2                                                ! MPI
         end if                                                             ! MPI
      end if

      if(npes >= 2) then
         deallocate(p1Sp2_w2)
         if(mod_pot == VANDERBILT_TYPE) then
            deallocate(p1Sp2_w)
         end if
      end if

    end subroutine Psi1SPhi2_t

    subroutine modify_bsd_and_phi_t(i)
      integer, intent(in) :: i

      integer             :: j, ia
      real(DP)            :: sr, si

      if(mod_pot == VANDERBILT_TYPE) then
         do j = 1, np_e         ! MPI
            if(to == OTHER_BANDS .and. j+ista_e-1 == i) cycle
            if(to == SAME_BAND .and. j+ista_e-1 /= i) cycle
            if(kimg == 1) then
               do ia = 1, nlmta
! ==== DEBUG by tkato 2011/09/02 ===============================================
!                 bsdr_l(j,ia,ik) = bsdr_l(j,ia,ik)-p1Sp2(j,1)*bp_w(ia,1)
!                 bsdi_l(j,ia,ik) = bsdi_l(j,ia,ik)-p1Sp2(j,1)*bp_w(ia,2)
                  bsdr_l(j,ia,ik) = bsdr_l(j,ia,ik)-p1Sp2(j+ista_e-1,1)*bp_w(ia,1)
                  bsdi_l(j,ia,ik) = bsdi_l(j,ia,ik)-p1Sp2(j+ista_e-1,1)*bp_w(ia,2)
! ==============================================================================
               end do
            else if(kimg == 2) then
               if(k_symmetry(ik) == GAMMA) then
                  do ia = 1, nlmta
                     sr = bp_w(ia,1)
! ==== DEBUG by tkato 2011/09/02 ===============================================
!                    bsdr_l(j,ia,ik) = bsdr_l(j,ia,ik)-p1Sp2(j,1)*sr
                     bsdr_l(j,ia,ik) = bsdr_l(j,ia,ik)-p1Sp2(j+ista_e-1,1)*sr
! ==============================================================================
                  end do
               else
                  do ia = 1, nlmta
                     sr = bp_w(ia,1);  si = bp_w(ia,2)
! ==== DEBUG by tkato 2011/09/02 ===============================================
!                    bsdr_l(j,ia,ik) = bsdr_l(j,ia,ik)-p1Sp2(j,1)*sr+p1Sp2(j,2)*si
!                    bsdi_l(j,ia,ik) = bsdi_l(j,ia,ik)-p1Sp2(j,1)*si-p1Sp2(j,2)*sr
                     bsdr_l(j,ia,ik) = bsdr_l(j,ia,ik)-p1Sp2(j+ista_e-1,1)*sr+p1Sp2(j+ista_e-1,2)*si
                     bsdi_l(j,ia,ik) = bsdi_l(j,ia,ik)-p1Sp2(j+ista_e-1,1)*si-p1Sp2(j+ista_e-1,2)*sr
! ==============================================================================
                  end do
               end if
            end if
         end do
      end if

      do j = 1, neg
         if(to == OTHER_BANDS .and. j == i) cycle
         if(to == SAME_BAND .and. j /= i) cycle
         if(kimg == 1) then
            do ia = 1, np_g1k(ik)                                        ! MPI
               phi_t(ia,j,1) = phi_t(ia,j,1) - p1Sp2(j,1)*psi_t(ia,i,1)
            end do
         else if(kimg == 2) then
            if(k_symmetry(ik) == GAMMA) then
               do ia = 1, np_g1k(ik)                                         ! MPI
                  sr = psi_t(ia,i,1) ;  si = psi_t(ia,i,2)
                  phi_t(ia,j,1) = phi_t(ia,j,1) - p1Sp2(j,1)*sr
                  phi_t(ia,j,2) = phi_t(ia,j,2) - p1Sp2(j,1)*si
               end do
            else
               do ia = 1, np_g1k(ik)                                       ! MPI
                  sr = psi_t(ia,i,1) ;  si = psi_t(ia,i,2)
                  phi_t(ia,j,1) = phi_t(ia,j,1) - p1Sp2(j,1)*sr+p1Sp2(j,2)*si
                  phi_t(ia,j,2) = phi_t(ia,j,2) - p1Sp2(j,1)*si-p1Sp2(j,2)*sr
               end do
            end if
         end if
      end do
    end subroutine modify_bsd_and_phi_t

  end subroutine mgs_sd2wf_each_k_G


! =============================== added by K. Tagami ================ 11.0
  subroutine mgs_sd2wf_each_k_G_noncl( ik, to, phi_l, ikst, iken, &
	&                              bsdr_l, bsdi_l, mod_pot )
    integer, intent(in)             :: ik, to, ikst, iken
    real(kind=DP), intent(inout)    :: phi_l(kg1,np_e,ikst:iken,kimg)
    real(kind=DP), optional,intent(inout), dimension(np_e,nlmta,ikst:iken) :: bsdr_l
    real(kind=DP), optional,intent(inout), dimension(np_e,nlmta,ikst:iken) :: bsdi_l

    integer, intent(in)             :: mod_pot
    integer               :: i, is

    integer :: kimg_t
    real(kind=DP),allocatable,dimension(:,:,:)    :: bp_w_noncl
    real(kind=DP),allocatable, dimension(:,:,:,:) :: phi_t_noncl

    call m_ESortho_mgs_alloc(ik)
    allocate(phi_t_noncl(np_g1k_x,neg,kimg,ndim_spinor))

    if(k_symmetry(ik) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = 2
    end if
    allocate(bp_w_noncl(nlmta,kimg_t,ndim_spinor))

    Do is=1, ndim_spinor
      call m_ES_W_transpose_r(.false., ista_k, iend_k, ik+is-1, zaj_l, psi_t_noncl(:,:,:,is) )
                  ! zaj_l(ig,ie,ik,ri) -> psi_t(ig,ie,ri) (transposed)
    End do

    Do is=1, ndim_spinor
      call m_ES_W_transpose_r(.false., ikst, iken, ik+is-1, phi_l, phi_t_noncl(:,:,:,is) )
                   ! phi_l(ig,ie,ik,ri) -> phi_t(ig,ie,ri) (transposed)
    End do
!
    do i = 1, neg
       if(mod_pot == VDB) call broadcast_fs_noncl(i)  ! (fs(ri)_l(i,*,ik)-> bp_w
       call Psi1SPhi2_t_noncl(ik,i)                      ! -(contained here) -> p1Sp2
       call modify_bsd_and_phi_t_noncl(i)                ! -(c.h.) p1Sp2,phi_t -> bsd(ri),phi_t
    end do

    Do is=1, ndim_spinor
      call m_ES_W_transpose_back_r(.false., ikst, iken, ik+is-1, phi_l, phi_t_noncl(:,:,:,is) )
                                             ! phi_t -> phi_l
    End do
    deallocate(phi_t_noncl)
    deallocate(bp_w_noncl)

!xx call m_ES_mgs_dealloc
    call m_ESortho_mgs_dealloc

  contains
    subroutine broadcast_fs_noncl(io)
      integer, intent(in) :: io
      integer  :: ia, i, kimg_t

      i = map_z(io)                                 ! MPI
      if(k_symmetry(ik) == GAMMA) then
         kimg_t = 1
      else
         kimg_t = 2
      end if

! ======= katDEBUG ===================================================== 11.0
!      if(map_e(i) == myrank_e) then                 ! MPI
      if (map_e(io) == myrank_e) then                 ! MPI
! ====================================================================== 11.0
         if(k_symmetry(ik) == GAMMA) then
            do ia = 1, nlmta
	      Do is=1, ndim_spinor
                 bp_w_noncl(ia,1,is) = fsr_l(i,ia,ik+is-1)
              End do
            end do
         else
            do ia = 1, nlmta
	      Do is=1, ndim_spinor
                 bp_w_noncl(ia,1,is) = fsr_l(i,ia,ik+is-1)
                 bp_w_noncl(ia,2,is) = fsi_l(i,ia,ik+is-1)
              End do
            end do
         end if
      end if                                        ! MPI

! ======== katDEBUG ====================================================== 11.0
!      call mpi_bcast( bp_w_noncl, nlmta*kimg_t*ndim_spinor, mpi_double_precision, &
!	   &          map_e(i), mpi_k_world(myrank_k), ierr )           ! MPI
      call mpi_bcast( bp_w_noncl, nlmta*kimg_t*ndim_spinor, mpi_double_precision, &
	   &          map_e(io), mpi_k_world(myrank_k), ierr )           ! MPI
! =================================================================== 11.0

    end subroutine broadcast_fs_noncl

    subroutine Psi1SPhi2_t_noncl(ik,i)
      integer, intent(in) :: ik,i

      integer   :: j,ia,p,q, kimg_t
      integer :: is1, is2, k1, k2
      integer :: is_tmp

      real(DP)  :: ar,ai, c1, c2, c3, c4

      real(kind=DP),allocatable,dimension(:,:)  :: p1Sp2_w, p1Sp2_w2

      if((k_symmetry(1) == GAMMA .and. kimg == 2).or.kimg==1) then
         kimg_t = 1
      else
         kimg_t = 2
      end if

      if(npes >= 2) then
         if(mod_pot == VANDERBILT_TYPE) then
            allocate(p1Sp2_w(neg,kimg_t))
         end if
         allocate(p1Sp2_w2(neg,kimg_t))
      end if

      p1Sp2 = 0.d0
      if ( mod_pot == VANDERBILT_TYPE ) then
         do j = 1, np_e                             ! MPI
            if (to == OTHER_BANDS .and. j+ista_e-1 == i) cycle               ! MPI
            if (to == SAME_BAND .and. j+ista_e-1 /= i) cycle
            ar = 0.d0; if(kimg == 2) ai = 0.d0

            if (kimg == 1) then
               do ia = 1, nac
                  p = nlmta1(ia);      q = nlmta2(ia)
                  Do is1=1, ndim_spinor
                     Do is2=1, ndim_spinor
                        is_tmp = ( is1 -1 ) *ndim_spinor + is2
                        k2 = ik + is2 -1

                        c1 = real( fqwei_noncl(ia,is_tmp) ) &
                             &   *( bp_w_noncl(p,1,is1) *bsdr_l(j,q,k2) &
                             &    + bp_w_noncl(p,2,is1) *bsdi_l(j,q,k2) )

                        c2 =-aimag( fqwei_noncl(ia,is_tmp) ) &
                             &   *( bp_w_noncl(p,1,is1) *bsdi_l(j,q,k2) &
                             &    + bp_w_noncl(p,2,is1) *bsdr_l(j,q,k2) )

                        ar = ar + c1 + c2

                     End do
                  End do
               end do
! ======== katDEBUG ========================================== 11.0
!               p1Sp2( j,1 ) = ar
               p1Sp2( j +ista_e -1,1 ) = ar
! ============================================================ 11.0
            else if(kimg == 2) then
               if(k_symmetry(ik) == GAMMA) then
                  do ia = 1, nac
                     p = nlmta1(ia);      q = nlmta2(ia)
                     Do is1=1, ndim_spinor
                        Do is2=1, ndim_spinor
                           is_tmp = ( is1 -1 ) *ndim_spinor + is2
                           k2 = ik + is2 -1

                           c1 = real( fqwei_noncl(ia,is_tmp) ) &
                                &   * bp_w_noncl(p,1,is1) *bsdr_l(j,q,k2)
                           c2 =-aimag( fqwei_noncl(ia,is_tmp) ) &
                                &   *bp_w_noncl(p,2,is1) *bsdr_l(j,q,k2)

                           ar = ar + c1 + c2

                        End Do
                     End do
                  end do

! ======= katDEBUG ============================================== 11.0
!               p1Sp2( j,1 ) = ar
                  p1Sp2( j+ista_e -1,1) = ar
! =============================================================== 11.0
               else
                  do ia = 1, nac
                     p = nlmta1(ia);      q = nlmta2(ia)
                     Do is1=1, ndim_spinor
                        Do is2=1, ndim_spinor
                           is_tmp = ( is1 -1 ) *ndim_spinor + is2
                           k2 = ik + is2 -1

                           c1 = real( fqwei_noncl(ia,is_tmp) ) &
                                &   *( bp_w_noncl(p,1,is1) *bsdr_l(j,q,k2) &
                                &    + bp_w_noncl(p,2,is1) *bsdi_l(j,q,k2) )
                           c2 =-aimag( fqwei_noncl(ia,is_tmp) ) &
                                &   *( bp_w_noncl(p,1,is1) *bsdi_l(j,q,k2) &
                                &    + bp_w_noncl(p,2,is1) *bsdr_l(j,q,k2) )

                           ar = ar + c1 + c2

                           c3 = real( fqwei_noncl(ia,is_tmp) ) &
                                &   *( bp_w_noncl(p,1,is1) *bsdi_l(j,q,k2) &
                                &    - bp_w_noncl(p,2,is1) *bsdr_l(j,q,k2) )
                           c4 = aimag( fqwei_noncl(ia,is_tmp) ) &
                                &   *( bp_w_noncl(p,1,is1) *bsdr_l(j,q,k2) &
                                &    + bp_w_noncl(p,2,is1) *bsdi_l(j,q,k2) )

                           ai = ai + c3 + c4

                        End Do
                     End do
                  end do
! ====== katDEBUG ============================================ 11.0
!                  p1Sp2(j,1) = ar; p1Sp2(j,2) = ai
                  p1Sp2(j+ista_e-1,1) = ar; p1Sp2(j+ista_e-1,2) = ai
! ============================================================= 11.0
               end if
            end if
         end do
!
         if (npes >= 2) then
            call mpi_allreduce( p1Sp2, p1Sp2_w, neg*kimg_t, mpi_double_precision, &
                 &              mpi_sum, mpi_k_world(myrank_k), ierr )  ! MPI
            p1Sp2 = 0.d0
         end if
      end if

      do j = 1, neg
         if (to == OTHER_BANDS .and. j == i) cycle
         if (to == SAME_BAND .and. j /= i) cycle
         if (kimg == 1) then
            do ia = 1, np_g1k(ik)                                           ! MPI
               Do is1=1, ndim_spinor
                  p1Sp2(j,1) = p1Sp2(j,1) + psi_t_noncl(ia,i,1,is1) &
                       &                   *phi_t_noncl(ia,j,1,is1)
               End do
            end do
         else if(kimg == 2) then
            if (k_symmetry(ik) == GAMMA) then
               do ia = 2, np_g1k(ik)                                        ! MPI
                  Do is1=1, ndim_spinor
                     ar = psi_t_noncl(ia,i,1,is1)
                     ai = psi_t_noncl(ia,i,2,is1)
                     p1Sp2(j,1) = p1Sp2(j,1) + ( ar *phi_t_noncl(ia,j,1,is1) &
                          &                    + ai *phi_t_noncl(ia,j,2,is1) )*2.d0
                  End do
               end do
!!$               if(mype /= 0) then
               if (myrank_e /= 0) then
                  Do is1=1, ndim_spinor
                     ar = psi_t_noncl(1,i,1,is1)
                     ai = psi_t_noncl(1,i,2,is1)
                     p1Sp2(j,1) = p1Sp2(j,1) + ( ar *phi_t_noncl(1,j,1,is1) &
                          &                    + ai *phi_t_noncl(1,j,2,is1) )*2.d0
                  End do
               else
                  Do is1=1, ndim_spinor
                     p1Sp2(j,1) = p1Sp2(j,1) + psi_t_noncl(1,i,1,is1) &
! ===== katDEBUG ======================================================== 11.0
!                          &                   *psi_t_noncl(1,j,1,is1)
                          &                   *phi_t_noncl(1,j,1,is1)
! ======================================================================= 11.0
                  End do
               end if
            else
               do ia = 1, np_g1k(ik)                                         ! MPI
                  Do is1=1, ndim_spinor
                     ar = psi_t_noncl(ia,i,1,is1)
                     ai = psi_t_noncl(ia,i,2,is1)
                     p1Sp2(j,1) = p1Sp2(j,1) + ar *phi_t_noncl(ia,j,1,is1) &
                          &                  + ai *phi_t_noncl(ia,j,2,is1)
                     p1Sp2(j,2) = p1Sp2(j,2) + ar *phi_t_noncl(ia,j,2,is1) &
                          &                  - ai *phi_t_noncl(ia,j,1,is1)
                  End do
               end do
            end if
         end if
      end do
!
      if ( npes >= 2 ) then
         call mpi_allreduce( p1Sp2, p1Sp2_w2, neg*kimg_t, mpi_double_precision, &
              &              mpi_sum, mpi_k_world(myrank_k), ierr )  ! MPI
         if (mod_pot == VANDERBILT_TYPE) then                                ! MPI
            p1Sp2 = p1Sp2_w + p1Sp2_w2                                      ! MPI
         else                                                               ! MPI
            p1Sp2 = p1Sp2_w2                                                ! MPI
         end if                                                             ! MPI
      end if

      if (npes >= 2) then
         deallocate(p1Sp2_w2)
         if (mod_pot == VANDERBILT_TYPE) then
            deallocate(p1Sp2_w)
         end if
      end if

    end subroutine Psi1SPhi2_t_noncl

   subroutine modify_bsd_and_phi_t_noncl(i)
      integer, intent(in) :: i

      integer :: is1, k1
      integer             :: j, ia
      real(DP)            :: sr, si

      if (mod_pot == VANDERBILT_TYPE) then
         do j = 1, np_e         ! MPI
            if (to == OTHER_BANDS .and. j+ista_e-1 == i) cycle
            if (to == SAME_BAND .and. j+ista_e-1 /= i) cycle
            if (kimg == 1) then
               do ia = 1, nlmta
                  Do is1=1, ndim_spinor
                     k1 = ik + is1 -1
                     bsdr_l(j,ia,k1) = bsdr_l(j,ia,k1) &
! ============ katDEBUG =============================================== 11.0
!!                          &          -p1Sp2(j,1) *bp_w_noncl(ia,1,is1)
                          &          -p1Sp2(j+ista_e-1,1) *bp_w_noncl(ia,1,is1)
! ==================================================================== 11.0

                     bsdi_l(j,ia,k1) = bsdi_l(j,ia,k1) &
! =========== katDEBUG ============================================= 11.0
!!                          &          -p1Sp2(j,1) *bp_w_noncl(ia,2,is1)
                          &          -p1Sp2(j+ista_e-1,1) *bp_w_noncl(ia,2,is1)
! ================================================================== 11.0
                  End do
               end do
            else if(kimg == 2) then
               if(k_symmetry(ik) == GAMMA) then
                  do ia = 1, nlmta
                     Do is1=1, ndim_spinor
                        k1 = ik + is1 -1
                        sr = bp_w_noncl(ia,1,is1)
! ===== katDEBUG ==================================================== 11.0
!                        bsdr_l(j,ia,k1) = bsdr_l(j,ia,k1) -p1Sp2(j,1)*sr
                        bsdr_l(j,ia,k1) = bsdr_l(j,ia,k1) -p1Sp2(j+ista_e-1,1)*sr
! =================================================================== 11.0
                     End do
                  end do
               else
                  do ia = 1, nlmta
                     Do is1=1, ndim_spinor
                        k1 = ik + is1 -1
                        sr = bp_w_noncl(ia,1,is1)
                        si = bp_w_noncl(ia,2,is1)
                        bsdr_l(j,ia,k1) = bsdr_l(j,ia,k1) &
! ======= katDEBUG ======================================================= 11.0
!                             &           -p1Sp2(j,1)*sr +p1Sp2(j,2)*si
                             &           -p1Sp2(j+ista_e-1,1)*sr +p1Sp2(j+ista_e-1,2)*si
!========================================================================= 11.0

                        bsdi_l(j,ia,k1) = bsdi_l(j,ia,k1) &
! =========katDEBUG ====================================================== 11.0
!                             &           -p1Sp2(j,1)*si -p1Sp2(j,2)*sr
                             &           -p1Sp2(j+ista_e-1,1)*si -p1Sp2(j+ista_e-1,2)*sr
! ======================================================================== 11.0
                     End do
                  end do
               end if
            end if
         end do
      end if

      do j = 1, neg
         if (to == OTHER_BANDS .and. j == i) cycle
         if (to == SAME_BAND .and. j /= i) cycle
         if (kimg == 1) then
            do ia = 1, np_g1k(ik)                                        ! MPI
               Do is1=1, ndim_spinor
                  k1 = ik + is1 -1
                  phi_t_noncl(ia,j,1,is1) = phi_t_noncl(ia,j,1,is1) &
                       &                   -p1Sp2(j,1) *psi_t_noncl(ia,i,1,is1)
               End do
            end do
         else if(kimg == 2) then
            if(k_symmetry(ik) == GAMMA) then
               do ia = 1, np_g1k(ik)                                         ! MPI
                  Do is1=1, ndim_spinor
                     k1 = ik + is1 -1
                     sr = psi_t_noncl(ia,i,1,is1)
                     si = psi_t_noncl(ia,i,2,is1)

                     phi_t_noncl(ia,j,1,is1) = phi_t_noncl(ia,j,1,is1)&
                          &                   -p1Sp2(j,1) *sr
                     phi_t_noncl(ia,j,2,is1) = phi_t_noncl(ia,j,2,is1) &
                          &                   -p1Sp2(j,1) *si
                  End do
               end do
            else
               do ia = 1, np_g1k(ik)                                       ! MPI
                  Do is1=1, ndim_spinor
                     k1 = ik + is1 -1
                     sr = psi_t_noncl(ia,i,1,is1)
                     si = psi_t_noncl(ia,i,2,is1)

                     phi_t_noncl(ia,j,1,is1) = phi_t_noncl(ia,j,1,is1)&
                          &                   -p1Sp2(j,1)*sr +p1Sp2(j,2)*si
                     phi_t_noncl(ia,j,2,is1) = phi_t_noncl(ia,j,2,is1)&
                          &                   -p1Sp2(j,1)*si -p1Sp2(j,2)*sr
                  End do
               end do
            end if
         end if
      end do
    end subroutine modify_bsd_and_phi_t_noncl

  end subroutine mgs_sd2wf_each_k_G_noncl
! ========================================================================== 11.0



  subroutine orthogonalize_SD(ik,wfsd_l,bsdr_l,bsdi_l)
    integer, intent(in)                                     :: ik
    real(kind=DP),intent(inout),dimension(kg1,np_e,ista_k:iend_k,kimg):: wfsd_l ! MPI
    real(kind=DP),intent(inout),dimension(np_e,nlmta,1)    :: bsdr_l,bsdi_l
    integer  :: id_sname = -1
    call tstatc0_begin('orthogonalize_SD ',id_sname)

!!$    call m_ES_mgs_alloc(ik)
    if(modnrm == EXECUT) then
       call mgs_4_each_k_G(ista_k,iend_k,ik,wfsd_l,ORTHONORMALIZATION &
            & ,bsdr_l,bsdi_l,mod_pot=VDB)!-(m_E.S.)
    else
       call mgs_4_each_k_G(ista_k,iend_k,ik,wfsd_l,ORTHONORMALIZATION &
            & ,mod_pot=NORMCONSERVATION)!-(m_E.S.)
    end if
!!$    call m_ES_mgs_dealloc()

    call tstatc0_end(id_sname)
  end subroutine orthogonalize_SD

  subroutine m_ES_MGS_4_each_k(nfout,ik,mode)
    integer, intent(in) :: nfout,mode,ik ! mode={ORTHONORMALIZATION | NORMALIZATION}
    integer :: id_sname = -1,i
                                                  __TIMER_SUB_START(501)
    call wd_title_of_the_operation()                        !-(c.h.)

!!$    call m_ES_mgs_alloc(ik)
    if(mode == ORTHONORMALIZATION) then
       call tstatc0_begin('modified_gram_schmidt ',id_sname,1)
    else
       call tstatc0_begin('modified_gram_schmidt (normalization only)',id_sname,1)
    endif
    if(modnrm == EXECUT) then
       call m_ES_betar_dot_WFs_4_each_k(nfout,ik)   ! -> fsr_l,fsi_l
       call mgs_4_each_k_G(ista_k,iend_k,ik,zaj_l,mode,fsr_l,fsi_l,mod_pot=VDB)!-(m_E.S.)
       if(ipri>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')
       call tstatc0_end(id_sname)
!!$       if(mode==NORMALIZATION) call m_ES_betar_dot_WFs_4_each_k(nfout,ik)   ! -> fsr_l,fsi_l
    else
       call mgs_4_each_k_G(ista_k,iend_k,ik,zaj_l,mode,mod_pot=NORMCONSERVATION)!-(m_E.S.)
       if(ipri>=2 .and. ik==1) write(nfout,'(" !### zaj_l is new,  bfft is old")')
       call tstatc0_end(id_sname)
       call m_ES_betar_dot_WFs_4_each_k(nfout,ik)   ! -> fsr_l,fsi_l
    end if
#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) status_saved_phifftr(:,ik) = OLD
#endif
!!$    call m_ES_mgs_dealloc
                                                  __TIMER_SUB_STOP(501)
  contains
    subroutine wd_title_of_the_operation
      integer, save :: iflag = 0
       if(iflag == 0 .and. ipri >= 2 .and. ik == 1) then
          if(modnrm == EXECUT) then
             write(nfout,*) ' <<< modified_gram_schmidt_vanderbilt_type >>>'
          else
             write(nfout,*) ' <<< modified_gram_schmidt_norm_conserve >>>'
          end if
       end if
       if(iflag == 0) iflag = 1
     end subroutine wd_title_of_the_operation
  end subroutine m_ES_MGS_4_each_k

! =================================== added by K. Tagami ============ 11.0
  subroutine m_ES_MGS_4_each_k_noncl( nfout,ik,mode )
    integer, intent(in) :: nfout,mode,ik ! mode={ORTHONORMALIZATION | NORMALIZATION}

    integer :: id_sname = -1,i
    integer :: is

    call wd_title_of_the_operation()                        !-(c.h.)

    if (modnrm == EXECUT) then
       call tstatc0_begin('modified_gram_schmidt_noncl ',id_sname,1)

       Do is=1, ndim_spinor
          call m_ES_betar_dot_WFs_4_each_k(nfout,ik+is-1)   ! -> fsr_l,fsi_l
       End do

       call mgs_4_each_k_G_noncl( ista_k, iend_k, ik, zaj_l, mode, &
                 &                fsr_l, fsi_l, mod_pot=VDB )

       call tstatc0_end(id_sname)

    else
       call tstatc0_begin('modified_gram_schmidt_noncl ',id_sname,1)
       call mgs_4_each_k_G_noncl( ista_k, iend_k, ik, zaj_l, mode, &
      &                           mod_pot=NORMCONSERVATION )

       call tstatc0_end(id_sname)

       Do is=1, ndim_spinor
         call m_ES_betar_dot_WFs_4_each_k(nfout,ik+is-1)   ! -> fsr_l,fsi_l
       End do
    end if
#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) then
       do is=1,ndim_spinor
          status_saved_phifftr(:,ik+is-1) = OLD
       end do
    end if
#endif

  contains
    subroutine wd_title_of_the_operation
      integer, save :: iflag = 0
       if(iflag == 0 .and. ipri >= 1 .and. ik == 1) then
          if(modnrm == EXECUT) then
             write(nfout,*) ' <<< modified_gram_schmidt_vanderbilt_type_noncl >>>'
          else
             write(nfout,*) ' <<< modified_gram_schmidt_norm_conserve_noncl >>>'
          end if
       end if
       if(iflag == 0) iflag = 1
     end subroutine wd_title_of_the_operation

  end subroutine m_ES_MGS_4_each_k_noncl
! ============================================================================ 11.0

!!$!BRANCH_P 3D_Parallel
!!$  subroutine m_ES_MGS_4_each_k_3D(nfout,ik,mode)
!!$    integer, intent(in) :: nfout,mode,ik ! mode={ORTHONORMALIZATION | NORMALIZATION}
!!$    integer :: id_sname = -1,i
!!$#ifdef __TIMER_SUB__
!!$  call timer_sta(502)
!!$#endif
!!$
!!$    call wd_title_of_the_operation()                        !-(c.h.)
!!$
!!$       call tstatc0_begin('modified_gram_schmidt ',id_sname,1)
!!$    if(modnrm == EXECUT) then
!!$       call m_ES_betar_dot_WFs_4_each_k_3D(nfout,ik)   ! -> fsr_l,fsi_l
!!$       call mgs_4_each_k_G_3D(ista_k,iend_k,ik,zaj_l,mode &
!!$      &                      ,fsr_l,fsi_l,mod_pot=VDB)!-(m_E.S.)
!!$       call tstatc0_end(id_sname)
!!$    else
!!$       call mgs_4_each_k_G_3D(ista_k,iend_k,ik,zaj_l,mode &
!!$      &                      ,mod_pot=NORMCONSERVATION)!-(m_E.S.)
!!$       call tstatc0_end(id_sname)
!!$       call m_ES_betar_dot_WFs_4_each_k_3D(nfout,ik)   ! -> fsr_l,fsi_l
!!$    end if
!!$
!!$#ifdef __TIMER_SUB__
!!$  call timer_end(502)
!!$#endif
!!$  contains
!!$    subroutine wd_title_of_the_operation
!!$      integer, save :: iflag = 0
!!$       if(iflag == 0 .and. ipri >= 1 .and. ik == 1) then
!!$          if(modnrm == EXECUT) then
!!$             write(nfout,*) ' <<< modified_gram_schmidt_vanderbilt_type >>>'
!!$          else
!!$             write(nfout,*) ' <<< modified_gram_schmidt_norm_conserve >>>'
!!$          end if
!!$       end if
!!$       if(iflag == 0) iflag = 1
!!$     end subroutine wd_title_of_the_operation
!!$  end subroutine m_ES_MGS_4_each_k_3D
!!$!BRANCH_P_END 3D_Parallel

  subroutine mgs_4_each_k_G(k1,k2,ik,psi_l,mode,bpr_l,bpi_l,mod_pot)
! Revised by T. Yamasaki in April 2006
    integer, intent(in)                                :: k1,k2,ik
    real(kind=DP),dimension(kg1,np_e,k1:k2,kimg)      :: psi_l
    integer, intent(in)                                :: mode,mod_pot
    real(kind=DP),intent(inout),optional,dimension(np_e,nlmta,k1:k2) :: bpr_l,bpi_l
!!$    real(kind=DP),optional,dimension(np_e,nlmta,k1:k2) :: bpr_l,bpi_l

    integer ::       i
    real(kind=DP) :: fr
    integer  :: kimg_t_wk
#ifdef MGS_DGEMM
    real(kind=DP),allocatable,dimension(:,:,:) :: psi_t_dia
    real(kind=DP),allocatable,dimension(:,:)   :: bpr_t_dia, bpi_t_dia
    integer ::       NB_END, NB_END2, i1, i2, ri, iy, ix
!!$    real(kind=DP), allocatable, dimension(:,:) :: bpr_tw1_BLAS
    real(kind=DP), allocatable, dimension(:,:,:) :: p1Sp2_NB
    integer, save :: ibsize_print = OFF
#endif
!!$    real(kind=DP), allocatable, dimension(:,:) :: psi_iri

#ifndef SX
! NEC tune ------------------------------->
    integer :: ibl1,ibl2,ibsize,ncache
    ncache = (m_CtrlP_cachesize()*1024)*3/4
! NEC tune <-------------------------------
#endif

#ifdef MGS_DGEMM
    if(nblocksize_mgs_is_given) then
       NB = nblocksize_mgs
    else
       NB = nblocksize_mgs_default
    end if
    if(ipri >= 2) then
       if(ibsize_print == OFF) then
          if(nblocksize_mgs_is_given) then
             write(nfout,'(" ! nblocksize_mgs_is_given")')
          else
             write(nfout,'(" ! nblocksize_mgs_is_given is false")')
          end if
          write(nfout,'( "! NB(=nblocksize_mgs) (mgs_4_each_k_G) = ",i8)') NB
          ibsize_print = ON
       end if
    end if
#endif

#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
!xx    call m_ES_mgs_alloc(ik)
    kimg_t_wk = 2
    if(k_symmetry(ik) == GAMMA) kimg_t_wk = 1
    call m_ESortho_mgs_alloc(ik)
!!$    allocate(psi_iri(np_g1k_x,kimg))
!!$    kimg_t_wk = kimg
!!$    if(kimg==2 .and. k_symmetry(ik) == GAMMA) kimg_t_wk = 1
#ifdef MGS_DGEMM
    allocate(p1Sp2_NB(NB,NB,kimg_t_wk))
    allocate(psi_t_dia(np_g1k_x,NB,kimg))
#endif
    if(mod_pot == VANDERBILT_TYPE) then
       if(kimg_t_wk == 1) then
          call m_ES_F_transpose_r(.true.,k1,k2,ik,bpr_l,bpr_t)     ! bpr_l -> bpr_t
#ifdef MGS_DGEMM
          allocate(bpr_t_dia(np_fs_x,NB))
! === Debug by Intel "-check all" option! by T.Kato 2013/07/02 =================
          allocate(bpi_t_dia(1,1))
! ==============================================================================
#endif
       else
          call m_ES_F_transpose_r(.true.,k1,k2,ik,bpr_l,bpr_t,bpi_l,bpi_t)  ! bp[ri]_l -> bp[ri]_t
#ifdef MGS_DGEMM
          allocate(bpr_t_dia(np_fs_x,NB))
          allocate(bpi_t_dia(np_fs_x,NB))
#endif
       end if
    end if

    call m_ES_W_transpose_r(.true.,k1,k2,ik,psi_l,psi_t)    !-(m_E.S.) psi_ l-> psi_t
#ifdef MGS_DGEMM
    do i = 1, neg,NB
       NB_END = i + NB -1
       if( NB_END > neg ) NB_END = neg
!diagonal
    do i1 = i, NB_END
       if(mode == ORTHONORMALIZATION .or. mode == NORMALIZATION) then
          if(mod_pot == VANDERBILT_TYPE) then
             call WSW_t_g(ik,i1,mod_pot,fr,psi_t,np_g1k_x,neg,kimg,bpr_t,bpi_t)
          else
             call WSW_t_g(ik,i1,mod_pot,fr,psi_t,np_g1k_x,neg,kimg)
          endif

!!$                                    !   fr = 1/dsqrt(<Psi(i)|S|Psi(i)>)
          if(dabs(fr-1.d0) > DELTA) then
             if(mod_pot == VANDERBILT_TYPE) then
                call normalize_bp_and_psi_t_g(ik,i1,fr,mod_pot, &
                     & psi_t,np_g1k_x,neg,kimg,bpr_t,bpi_t)
             else
                call normalize_bp_and_psi_t_g(ik,i1,fr,mod_pot, &
                     & psi_t,np_g1k_x,neg,kimg)
             endif
          endif
          !   |Psi(i)> = |Psi(i)> * fr,  <beta|Psi(i)> = <beta|Psi(i)> * fr

          if(mod_pot == VANDERBILT_TYPE) call cp_bpr2bprtw(i1,neg)  ! bpr_t -> bpr_tw1, bpr_tw2
       end if
       if(mode /= NORMALIZATION) then
          if(i1 == neg) cycle
          call cp_psi2psii_g(ik,i1) ! psi_t(:,i,:) -> psi_ir,psi_ii
          call W1SW2_t_r_g(ik,i1,NB_END,mod_pot,psi_t,np_g1k_x,neg,kimg) ! -> p1Sp2
          if(mod_pot == VANDERBILT_TYPE) &
               & call cp_bpr2bpi_g(kimg_t_wk,i1,bpr_t,bpi_t) ! bpr_t(:,i1)->bp_ir,bp_ii

          if(mod_pot == VANDERBILT_TYPE) then
             call modify_bp_and_psi_t_r_g(ik,i1,NB_END,mod_pot &
                  & ,psi_t,np_g1k_x,neg,kimg,bpr_t,bpi_t)
             ! psi_t, bpr_t, pbi_t, p1Sp2 -> psi_t, bpr_t, bpi_t
          else
             call modify_bp_and_psi_t_r_g(ik,i1,NB_END,mod_pot &
                  & ,psi_t,np_g1k_x,neg,kimg )
          endif

       end if
    end do   ! i1-loop
    if(mode /= NORMALIZATION) then
!lower
       if(mod_pot == VANDERBILT_TYPE) then
          call cp_psi_bpri2dias_g(ik,i,kimg_t_wk,NB,psi_t_dia,np_g1k_x,NB,kimg,np_fs_x,bpr_t_dia,bpi_t_dia)
       else
          call cp_psi_bpri2dias_g(ik,i,kimg_t_wk,NB,psi_t_dia,np_g1k_x,NB,kimg,np_fs_x)
       end if

    do i2 = i+NB, neg,NB
       if(NB_END == neg) cycle
       NB_END2 = i2 + NB -1
       if( NB_END2 > neg ) NB_END2 = neg
       if(mode /= NORMALIZATION) then
          if(allocated(bpr_tw1) .and. allocated(bpi_tw1)) then
             call W1SW2_t_r_block_g(ik,i,i2,p1Sp2_NB,NB_END,NB_END2,kimg_t_wk & ! -> p1Sp2
                  & , mod_pot,psi_t,psi_t_dia,np_g1k_x,neg,kimg,bpr_tw1(1,i),bpi_tw1(1,i))
          else if(allocated(bpr_tw1)) then
             call W1SW2_t_r_block_g(ik,i,i2,p1Sp2_NB,NB_END,NB_END2,kimg_t_wk & ! -> p1Sp2
                  & , mod_pot,psi_t,psi_t_dia,np_g1k_x,neg,kimg,bpr_tw1(1,i))
          else
             call W1SW2_t_r_block_g(ik,i,i2,p1Sp2_NB,NB_END,NB_END2,kimg_t_wk & ! -> p1Sp2
                  & , mod_pot,psi_t,psi_t_dia,np_g1k_x,neg,kimg)
          end if
          if(nrank_e > 1) call mpi_allreduce(MPI_IN_PLACE, p1Sp2_NB, NB*NB*kimg_t_wk &
               & ,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)

          if(mod_pot == VANDERBILT_TYPE) then
             call modify_bp_and_psi_t_r_blk_g(ik,i,i2,p1Sp2_NB,NB_END,NB_END2,kimg_t_wk &
                  & , mod_pot,psi_t,psi_t_dia,np_g1k_x,neg,kimg &
                  & ,  bpr_t,bpi_t,bpr_t_dia,bpi_t_dia)
             !               psi_t, bpr_t, pbi_t, p1Sp2 -> psi_t, bpr_t, bpi_t
          else
             call modify_bp_and_psi_t_r_blk_g(ik,i,i2,p1Sp2_NB,NB_END,NB_END2,kimg_t_wk &
                  & , mod_pot,psi_t,psi_t_dia,np_g1k_x,neg,kimg )
          endif

       end if
    end do   ! i2-loop
    end if   ! if(mode /= NORMALIZATION)
    end do   ! i-loop
#else
    do i = 1, neg
       if(mode == ORTHONORMALIZATION .or. mode == NORMALIZATION) then
          call WSW_t_g(ik,i,mod_pot,fr,psi_t,np_g1k_x,neg,kimg,bpr_t,bpi_t)
!!$                                  !   fr = 1/dsqrt(<Psi(i)|S|Psi(i)>)
          if(dabs(fr-1.d0) > DELTA) &
               & call normalize_bp_and_psi_t_g(ik,i,fr,mod_pot &
               & ,psi_t,np_g1k_x,neg,kimg, bpr_t,bpi_t)
          !   |Psi(i)> = |Psi(i)> * fr,  <beta|Psi(i)> = <beta|Psi(i)> * fr
       end if
       if(mode /= NORMALIZATION) then
          if(i == neg) cycle
          call cp_psi2psii_g(ik,i) ! psi_t(:,i,:) -> psi_ir,psi_ii
          if(mod_pot == VANDERBILT_TYPE) &
               & call cp_bpr2bpi_g(kimg_t_wk, i,bpr_t,bpi_t)  ! -> bp_ir, bp_ii
          call W1SW2_t_r_g(ik,i,neg,mod_pot,psi_t,np_g1k_x,neg,kimg & ! -> p1Sp2
               & ,bpr_t,bpi_t)
          call modify_bp_and_psi_t_r_g(ik,i,neg,mod_pot &
               & ,psi_t,np_g1k_x,neg,kimg,bpr_t,bpi_t)
                         ! psi_t, bpr_t, pbi_t, p1Sp2 -> psi_t, bpr_t, bpi_t
       end if
    end do
#endif


    call m_ES_W_transpose_back_r(.true.,k1,k2,ik,psi_l,psi_t)
    if(mod_pot == VANDERBILT_TYPE) then
       if(kimg_t_wk == 1) then
          call m_ES_F_transpose_back_r(.true.,k1,k2,ik,bpr_l,bpr_t)
#ifdef MGS_DGEMM
          deallocate(bpr_t_dia)
! === Debug by Intel "-check all" option! by T.Kato 2013/07/02 =================
          deallocate(bpi_t_dia)
! ==============================================================================
#endif
       else
          call m_ES_F_transpose_back_r(.true.,k1,k2,ik,bpr_l,bpr_t,bpi_l,bpi_t)
#ifdef MGS_DGEMM
          deallocate(bpr_t_dia,bpi_t_dia)
#endif
       end if
    end if
#ifdef MGS_DGEMM
    deallocate(p1Sp2_NB)
    deallocate(psi_t_dia)
#endif
#else
    integer ::       nmax, ito
    integer, allocatable, dimension(:)         :: ib2to_a, ib2back_a
    allocate(ib2to_a(neg)); allocate(ib2back_a(neg))

!xx call m_ES_mgs_alloc(ik)
    call m_ESortho_mgs_alloc(ik)
!!$    if(mod_pot == VANDERBILT_TYPE) call m_ES_F_transpose(k1,k2,ik,bpr_l,bpi_l,bpr_t,bpi_t)
    if(mod_pot == VANDERBILT_TYPE) &
         & call m_ES_F_transpose_r(.false.,k1,k2,ik,bpr_l,bpr_t,bpi_l,bpi_t)
    call m_ES_W_transpose_r(.false.,k1,k2,ik,psi_l,psi_t)    !-(m_E.S.) psi_ l-> psi_t
    do i = 1, neg
       ito = neordr(i,ik)
       if(mode == ORTHONORMALIZATION .or. mode == NORMALIZATION) then
          call WSW_t_g(ik,ito,mod_pot,fr,psi_t,np_g1k_x,neg,kimg,bpr_t,bpi_t)  !fr = 1/dsqrt(<Psi(ito)|S|Psi(ito)>)
          if(dabs(fr-1.d0) > DELTA) &
               & call normalize_bp_and_psi_t_g(ik,ito,fr,mod_pot &
               & ,psi_t,np_g1k_x,neg,kimg, bpr_t,bpi_t)
          !   |Psi(ito)> = |Psi(ito)> * fr,  <beta|Psi(ito)> = <beta|Psi(ito)> * fr
       end if
       if(mode /= NORMALIZATION) then
          call substitute_jto_ib2back(i,nmax) !-(c.h.) ->ib2to_a,ib2back_a,
          if(nmax == 0) cycle
          call W1SW2_t(i,ito) ! -> p1Sp2
          call modify_bp_and_psi_t(i,ito) ! psi_t, bpr_t, pbi_t, p1Sp2 -> psi_t, bpr_t, bpi_t
       end if
    end do

    call m_ES_W_transpose_back_r(.false.,k1,k2,ik,psi_l,psi_t)
!!$    if(mod_pot == VANDERBILT_TYPE) call m_ES_F_transpose_back(k1,k2,ik,bpr_l,bpi_l,bpr_t,bpi_t)
    if(mod_pot == VANDERBILT_TYPE) &
         & call m_ES_F_transpose_back_r(.false.,k1,k2,ik,bpr_l,bpr_t,bpi_l,bpi_t)
    deallocate(ib2back_a,ib2to_a)
#endif
!xx call m_ES_mgs_dealloc()
    call m_ESortho_mgs_dealloc()

  contains

#ifdef MGS_DGEMM
    subroutine cp_bpr2bprtw(i1,i2)
      integer, intent(in) :: i1,i2
      integer :: i, ia, p, q, i_last
!!$      if(i1==1) then
!!$         i_last = i2
!!$      else
!!$         i_last = i1
!!$      end if
#ifdef MGS_DGEMM_DEBUG
      i_last = i1
#else
      i_last = i2
#endif
      if(k_symmetry(ik) == GAMMA) then
         do i = i1, i_last
            do ia = 1, nac_p
               p = nlmta1_p(ia);     q = nlmta2_p(ia)
               bpr_tw1(ia,i) = bpr_t(p,i)
               bpr_tw2(ia,i) = bpr_t(q,i)
            end do
         end do
      else
         do i = i1, i_last
            do ia = 1, nac_p
               p = nlmta1_p(ia);     q = nlmta2_p(ia)
               bpr_tw1(ia,i) = bpr_t(p,i)
               bpr_tw2(ia,i) = bpr_t(q,i)
               bpi_tw1(ia,i) = bpi_t(p,i)
               bpi_tw2(ia,i) = bpi_t(q,i)
            end do
         end do
      end if
    end subroutine cp_bpr2bprtw
#endif

#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT

#else

    subroutine substitute_jto_ib2back(i,nmax)
      integer, intent(in)  :: i
      integer, intent(out) :: nmax

      integer              ::jto, j
      jto = 0
      do j = 1, neg
         if(nrvf_ordr(j,ik) <= i) cycle
         jto = jto + 1
         ib2to_a(jto)  = j
         ib2back_a(j)  = jto
      end do
      nmax = jto
      if(ipri >= 2 .and. ik == 1) then
         write(nfout,'(" <<substitute_jto_ib2back>>")')
         do j = 1, neg
            write(nfout,'(" ib2back_a(",i3,") = ",i3)') j, ib2back_a(j)
         end do
         do j = 1, nmax
            write(nfout,'(" ib2to_a(",i3,") = ",i3)') j, ib2to_a(j)
         end do
      end if
    end subroutine substitute_jto_ib2back

    subroutine W1SW2_t(i,ito)
      integer, intent(in)        :: i,ito

      integer       :: j, ia, jto, p, q, kimg_t
      real(kind=DP) :: ar, ai
      real(kind=DP), allocatable, dimension(:,:)  :: p1Sp2_t1, p1Sp2_t2             ! MPI

!!$      call alloc_p1Sp2(nmax)
      if(npes > 1) then
         if((k_symmetry(ik) == GAMMA .and. kimg==2) .or. kimg == 1) then
            kimg_t = 1
         else
            kimg_t = 2
         end if
         allocate(p1Sp2_t2(nmax,kimg_t)); p1Sp2_t2 = 0.d0                  ! MPI
         allocate(p1Sp2_t1(nmax,kimg_t)); p1Sp2_t1 = 0.d0                  ! MPI
      end if

      p1Sp2 = 0.d0

      if(mod_pot == VANDERBILT_TYPE) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(jto,ia,ar,ai,p,q)
#endif
         do jto = 1, nmax                                                 ! MPI
            j = ib2to_a(jto)                                              ! MPI
            if(kimg == 1) then
               ar = 0.d0
               do ia = 1, nac_p
                  p = nlmta1_p(ia);            q = nlmta2_p(ia)
                  ar = ar + fqwei_p(ia)*(bpr_t(p,ito)*bpr_t(q,j)+bpi_t(p,ito)*bpi_t(q,j))
               end do
               p1Sp2(jto,1) = ar
            else if(kimg == 2) then
               if(k_symmetry(ik) == GAMMA) then
                  ar = 0.d0
                  do ia = 1, nac_p
                     p = nlmta1_p(ia);         q = nlmta2_p(ia)
                     ar = ar + fqwei_p(ia)*(bpr_t(p,ito)*bpr_t(q,j))
                  end do
                  p1Sp2(jto,1) = ar
               else
                  ar = 0.d0; ai = 0.d0
                  do ia = 1, nac_p
                     p = nlmta1_p(ia);         q = nlmta2_p(ia)
                     ar = ar + fqwei_p(ia)*(bpr_t(p,ito)*bpr_t(q,j)+bpi_t(p,ito)*bpi_t(q,j))
                     ai = ai + fqwei_p(ia)*(bpr_t(p,ito)*bpi_t(q,j)-bpi_t(p,ito)*bpr_t(q,j))
                  end do
                  p1Sp2(jto,1) = ar;  p1Sp2(jto,2) = ai
               end if
            end if
         end do
      end if

      if(ipri >= 2 .and. ik == 1) then
         if(nmax > 1) then
            write(nfout,'(" <<W1SW2_t>> ik = ",i9,"  i = ",i3, " ito = ", i3)') ik,i,ito
            write(nfout,'(" (real) ",6d11.3, 99(/8x,6d11.3))') (p1Sp2(jto,1),jto=1, nmax)
            if(kimg == 2 .and. k_symmetry(ik) /= GAMMA) &
                 & write(nfout,'(" (imag) ",6d11.3, 99(/8x,6d11.3))') (p1Sp2(jto,kimg),jto=1, nmax)
         end if
      end if

      if(kimg == 1) then
         do jto = 1, nmax
            j = ib2to_a(jto)
            do ia = 1, np_g1k(ik)                                           ! MPI
               p1Sp2(jto,1) = p1Sp2(jto,1) + psi_t(ia,ito,1)*psi_t(ia,j,1 ) ! MPI
            end do
         end do
      else if(kimg == 2) then
         if(k_symmetry(ik) == GAMMA) then
            do jto = 1, nmax
               j = ib2to_a(jto)
               do ia = 2, np_g1k(ik)
                  ar  = psi_t(ia,ito,1)
                  ai  = psi_t(ia,ito,2)
                  p1Sp2(jto,1) = p1Sp2(jto,1)+(ar*psi_t(ia,j,1)+ai*psi_t(ia,j,2))*2.d0 ! MPI
               end do
!!$               if(mype /= 0) then
               if(myrank_e /= 0) then
                  ar = psi_t(1,ito,1)
                  ai = psi_t(1,ito,2)
                  p1Sp2(jto,1) = p1Sp2(jto,1) + (ar*psi_t(1,j,1) + ai*psi_t(1,j,2))*2.d0
               else
                  p1Sp2(jto,1) = p1Sp2(jto,1) + psi_t(1,ito,1)*psi_t(1,j,1)
               end if
            end do
         else
            do jto = 1, nmax
               j = ib2to_a(jto)
               do ia = 1, np_g1k(ik)                                          ! MPI
                  ar  = psi_t(ia,ito,1)
                  ai  = psi_t(ia,ito,2)
                  p1Sp2(jto,1) = p1Sp2(jto,1)+ar*psi_t(ia,j,1)+ai*psi_t(ia,j,2) ! MPI
                  p1Sp2(jto,2) = p1Sp2(jto,2)+ar*psi_t(ia,j,2)-ai*psi_t(ia,j,1) ! MPI
               end do
            end do
         end if
      end if
      if(npes > 1 ) then
         do q = 1, kimg_t
            do jto = 1, nmax
               p1Sp2_t1(jto,q) = p1Sp2(jto,q)
            end do
         end do
         call mpi_allreduce(p1Sp2_t1, p1Sp2_t2,nmax*kimg_t,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
         do q = 1, kimg_t
            do jto = 1, nmax
               p1Sp2(jto,q) = p1Sp2_t2(jto,q)
            end do
         end do
      end if

      if(ipri >= 2 .and. ik == 1) then
         if(nmax > 1) then
            write(nfout,'(" <<W1SW2_t>> ik = ",i9,"  i = ",i3, " ito = ", i3)') ik,i,ito
            write(nfout,'(" (real) ",6d11.3, 99(/8x,6d11.3))') (p1Sp2(jto,1),jto=1, nmax)
            if(kimg == 2 .and. k_symmetry(ik) /= GAMMA) &
                 & write(nfout,'(" (imag) ",6d11.3, 99(/8x,6d11.3))') (p1Sp2(jto,kimg),jto=1, nmax)
         end if
      end if

      if(npes > 1) deallocate(p1Sp2_t1,p1Sp2_t2)
    end subroutine W1SW2_t

    subroutine modify_bp_and_psi_t(i,ito)
      integer, intent(in) :: i, ito
      integer             :: j,ia,jto
      real(kind=DP) :: sr, si

      if(mod_pot == VANDERBILT_TYPE) then
         do j = 1, neg
            if(nrvf_ordr(j,ik) <= i) cycle
            jto = ib2back_a(j)
            if(kimg == 1) then
               do ia = 1, np_fs
                  bpr_t(ia,j) = bpr_t(ia,j) - p1Sp2(jto,1)*bpr_t(ia,ito)
                  bpi_t(ia,j) = bpi_t(ia,j) - p1Sp2(jto,1)*bpi_t(ia,ito)
               end do
            else if(kimg == 2) then
               if(k_symmetry(ik) == GAMMA) then
                  do ia = 1, np_fs
                     bpr_t(ia,j) = bpr_t(ia,j) - p1Sp2(jto,1)*bpr_t(ia,ito)
                  end do
               else
                  do ia = 1, np_fs
                     sr  =  bpr_t(ia,ito);     si  =  bpi_t(ia,ito)
                     bpr_t(ia,j) = bpr_t(ia,j) - p1Sp2(jto,1)*sr+p1Sp2(jto,2)*si
                     bpi_t(ia,j) = bpi_t(ia,j) - p1Sp2(jto,1)*si-p1Sp2(jto,2)*sr
                  end do
               end if
            end if
         end do
      end if

#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!CDIR NOSYNC
!CDIR CONCUR(BY=1)
#endif
      do j = 1,neg
         if(nrvf_ordr(j,ik) <= i) cycle
         jto = ib2back_a(j)         ! <-- 99Jan19 T. Y.
         if(kimg == 1) then
            do ia = 1, np_g1k(ik)                 ! MPI
               psi_t(ia,j,1) = psi_t(ia,j,1) - p1Sp2(jto,1)*psi_t(ia,ito,1 )
            end do
         else if(kimg == 2) then
            if(k_symmetry(ik) == GAMMA) then
               do ia = 1, np_g1k(ik)
                  sr  =  psi_t(ia,ito,1 );   si  =  psi_t(ia,ito,2 )
                  psi_t(ia,j,1) = psi_t(ia,j,1) - p1Sp2(jto,1)*sr
                  psi_t(ia,j,2) = psi_t(ia,j,2) - p1Sp2(jto,1)*si
               end do
            else
               do ia = 1, np_g1k(ik)                ! MPI
                  sr  =  psi_t(ia,ito,1 );   si  =  psi_t(ia,ito,2 )
                  psi_t(ia,j,1) = psi_t(ia,j,1) - p1Sp2(jto,1)*sr+p1Sp2(jto,2)*si
                  psi_t(ia,j,2) = psi_t(ia,j,2) - p1Sp2(jto,1)*si-p1Sp2(jto,2)*sr
               end do
            end if
         end if
      end do
!!$      call dealloc_p1Sp2()
    end subroutine modify_bp_and_psi_t
#endif

  end subroutine mgs_4_each_k_G

! ======================================= added by K. Tagami ======== 11.0
  subroutine mgs_4_each_k_G_noncl(k1,k2,ik,psi_l,mode,bpr_l,bpi_l,mod_pot)
    integer, intent(in)                                :: k1,k2,ik
    integer, intent(in)                                :: mode,mod_pot

    integer :: is
    real(kind=DP),dimension(kg1,np_e,k1:k2,kimg)      :: psi_l
    real(kind=DP),intent(inout),optional,dimension(np_e,nlmta,k1:k2) :: bpr_l,bpi_l

    integer ::       i
    real(kind=DP) :: fr
    integer  :: kimg_t_wk
#ifdef MGS_DGEMM
    integer ::       NB_END, NB_END2, i1, i2
!!$    real(kind=DP), allocatable, dimension(:,:) :: bpr_tw1_BLAS
    real(kind=DP), allocatable, dimension(:,:,:) :: p1Sp2_NB
    integer, save :: ibsize_print = OFF
#endif
#ifndef SX
! NEC tune ------------------------------->
    integer :: ibl1,ibl2,ibsize,ncache
    ncache = (m_CtrlP_cachesize()*1024)*3/4
! NEC tune <-------------------------------
#endif

#ifdef MGS_DGEMM
    if(nblocksize_mgs_is_given) then
       NB = nblocksize_mgs
    else
       NB = nblocksize_mgs_default
    end if
    if(ipri >= 2) then
       if(ibsize_print == OFF) then
          if(nblocksize_mgs_is_given) then
             write(nfout,'(" ! nblocksize_mgs_is_given")')
          else
             write(nfout,'(" ! nblocksize_mgs_is_given is false")')
          end if
          write(nfout,'( "! NB(=nblocksize_mgs) (mgs_4_each_k_G_noncl) = ",i8)') NB
          ibsize_print = ON
       end if
    end if
#endif


#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
!xx    call m_ES_mgs_alloc(ik)

    call m_ESortho_mgs_alloc(ik)

!!$    kimg_t_wk = kimg
!!$    if(kimg==2 .and. k_symmetry(ik) == GAMMA) kimg_t_wk = 1
    kimg_t_wk = 2
    if(k_symmetry(ik) == GAMMA) kimg_t_wk = 1
#ifdef MGS_DGEMM
    allocate(p1Sp2_NB(NB,NB,kimg_t_wk))
#endif

    if(mod_pot == VANDERBILT_TYPE) then
       if(kimg_t_wk == 1) then
	  Do is=1, ndim_spinor
            call m_ES_F_transpose_r(.true.,k1,k2,ik+is-1,bpr_l,bpr_t_noncl(:,:,is) )
                                         ! bpr_l -> bpr_t
          End do
       else
	  Do is=1, ndim_spinor
            call m_ES_F_transpose_r(.true.,k1,k2,ik+is-1,bpr_l,bpr_t_noncl(:,:,is),&
	&                                             bpi_l,bpi_t_noncl(:,:,is) )
                                                  ! bp[ri]_l -> bp[ri]_t
          End do
       end if
    end if

    Do is=1, ndim_spinor
       call m_ES_W_transpose_r(.true., k1, k2, ik+is-1, psi_l, psi_t_noncl(:,:,:,is) )
                                           !-(m_E.S.) psi_ l-> psi_t
    End do

#ifdef MGS_DGEMM
    do i = 1, neg,NB
       NB_END = i + NB -1
       if( NB_END > neg ) NB_END = neg
!diagonal
       do i1 = i, NB_END
          if ( mode == ORTHONORMALIZATION .or. mode == NORMALIZATION ) then
            call WSW_t_noncl(i1,fr)            !-(c.h.) fr = 1/dsqrt(<Psi(i)|S|Psi(i)>)

            if ( dabs(fr-1.d0) > DELTA ) then
                call normalize_bp_and_psi_t_noncl(i1,fr) !-(c.h.)
               !   |Psi(i)> = |Psi(i)> * fr,  <beta|Psi(i)> = <beta|Psi(i)> * fr
            endif
            if ( mod_pot == VANDERBILT_TYPE ) then
                call cp_bpr2bprtw_noncl(i1,neg)  ! bpr_t -> bpr_tw1, bpr_tw2
            endif
          end if

          if ( mode /= NORMALIZATION ) then
             if (i1 == neg) cycle
             call W1SW2_t_r_noncl( i1,NB_END ) ! -> p1Sp2
             call modify_bp_and_psi_t_r_noncl(i1,NB_END)
                                ! psi_t, bpr_t, pbi_t, p1Sp2 -> psi_t, bpr_t, bpi_t
          end if
       end do   ! i1-loop

!lower
       do i2 = i+NB, neg,NB
          if(NB_END == neg) cycle
          NB_END2 = i2 + NB -1
          if( NB_END2 > neg ) NB_END2 = neg
          if ( mode /= NORMALIZATION ) then
             call W1SW2_t_r_block_noncl( i, i2, p1Sp2_NB, NB_END, NB_END2, &
                  &                      kimg_t_wk ) ! -> p1Sp2
             call modify_bp_psi_t_r_block_noncl( i, i2, p1Sp2_NB, NB_END,NB_END2, &
                  &                      kimg_t_wk )
                      ! psi_t, bpr_t, pbi_t, p1Sp2 -> psi_t, bpr_t, bpi_t
          end if
       end do   ! i2-loop
    end do   ! i-loop


#else
    do i = 1, neg
       if ( mode == ORTHONORMALIZATION .or. mode == NORMALIZATION ) then
          call WSW_t_noncl(i,fr)            !-(c.h.) fr = 1/dsqrt(<Psi(i)|S|Psi(i)>)
          if ( dabs(fr-1.d0) > DELTA ) then
             call normalize_bp_and_psi_t_noncl(i,fr) !-(c.h.)
          !   |Psi(i)> = |Psi(i)> * fr,  <beta|Psi(i)> = <beta|Psi(i)> * fr
          endif
       end if
       if ( mode /= NORMALIZATION ) then
          if(i == neg) cycle
          call W1SW2_t_r_noncl(i,neg) ! -> p1Sp2
          call modify_bp_and_psi_t_r_noncl(i,neg)
                      ! psi_t, bpr_t, pbi_t, p1Sp2 -> psi_t, bpr_t, bpi_t
       end if
    end do
#endif

    Do is=1, ndim_spinor
      call m_ES_W_transpose_back_r(.true., k1, k2, ik+is-1, psi_l, psi_t_noncl(:,:,:,is) )
    End do

    if ( mod_pot == VANDERBILT_TYPE ) then
       if ( kimg_t_wk == 1 ) then
          Do is=1, ndim_spinor
            call m_ES_F_transpose_back_r(.true.,k1,k2,ik+is-1, bpr_l, bpr_t_noncl(:,:,is) )
          End do
       else
          Do is=1, ndim_spinor
            call m_ES_F_transpose_back_r(.true.,k1,k2,ik+is-1, bpr_l, bpr_t_noncl(:,:,is),&
	                             &                   bpi_l, bpi_t_noncl(:,:,is) )
          End do
       end if
    end if
#ifdef MGS_DGEMM
    deallocate(p1Sp2_NB)
#endif
#else
    integer ::       nmax, ito
    integer, allocatable, dimension(:)         :: ib2to_a, ib2back_a
    allocate(ib2to_a(neg)); allocate(ib2back_a(neg))

!xx call m_ES_mgs_alloc(ik)
    call m_ESortho_mgs_alloc(ik)

    if (mod_pot == VANDERBILT_TYPE) then
       Do is=1, ndim_spinor
!!$         call m_ES_F_transpose( k1, k2, ik+is-1,  bpr_l, bpi_l, &
!!$            &                 bpr_t_noncl(:,:,:,is), bpi_t_noncl(:,:,is) )
         call m_ES_F_transpose_r(.false.,k1,k2,ik+is-1,bpr_l,bpr_t_noncl(:,:,is) &
              &                                       ,bpi_l,bpi_t_noncl(:,:,is) )
       End do
    endif

    Do is=1, ndim_spinor
      call m_ES_W_transpose_r(.false., k1, k2, ik+is-1, psi_l, psi_t_noncl(:,:,:,is) )
    End do

    do i = 1, neg
       ito = neordr(i,ik)
       if ( mode == ORTHONORMALIZATION .or. mode == NORMALIZATION ) then
          call WSW_t_noncl(ito,fr)        !-(c.h.) fr = 1/dsqrt(<Psi(ito)|S|Psi(ito)>)
          if ( dabs(fr-1.d0) > DELTA ) then
             call normalize_bp_psi_t_noncl(ito,fr) !-(c.h.)
            !   |Psi(ito)> = |Psi(ito)> * fr,  <beta|Psi(ito)> = <beta|Psi(ito)> * fr
          endif
       end if
       if ( mode /= NORMALIZATION ) then
          call substitute_jto_ib2back(i,nmax) !-(c.h.) ->ib2to_a,ib2back_a,
          if ( nmax == 0 ) cycle
          call W1SW2_t_noncl(i,ito) ! -> p1Sp2
          call modify_bp_and_psi_t_noncl(i,ito)
                    ! psi_t, bpr_t, pbi_t, p1Sp2 -> psi_t, bpr_t, bpi_t
       end if
    end do

    Do is=1, ndim_spinor
      call m_ES_W_transpose_back_r(.false., k1, k2, ik+is-1, psi_l, psi_t_noncl(:,:,:,is) )
    End do

    if ( mod_pot == VANDERBILT_TYPE ) then
       Do is=1, ndim_spinor
          call m_ES_F_transpose_back_r(.false.,k1,k2,ik+is-1,bpr_l,bpr_t_noncl(:,:,is) &
               & bpi_l, bpi_t_noncl(:,:,is) )
       End do
    endif
    deallocate(ib2back_a,ib2to_a)


#endif
!xx call m_ES_mgs_dealloc()
    call m_ESortho_mgs_dealloc()

  contains

    subroutine WSW_t_noncl(j,fr)
      integer, intent(in)        :: j
      real(kind=DP), intent(out) :: fr

      integer :: ia,p,q, i, ig1
      integer :: is1, is2, is_tmp
      real(kind=DP) :: c1, c2, c3, c4

      real(kind=DP) :: fr1, fr2

      fr = 0.d0
      if ( mod_pot == VANDERBILT_TYPE ) then
         if ( k_symmetry(ik) == GAMMA ) then
            do ia = 1, nac_p
               p = nlmta1_p(ia);     q = nlmta2_p(ia)
               Do is1=1, ndim_spinor
                  Do is2=1, ndim_spinor
                     is_tmp = ( is1 -1 ) *ndim_spinor + is2

                     c1 = real( fqwei_p_noncl(ia,is_tmp) ) &
                     &   *bpr_t_noncl(p,j,is1) *bpr_t_noncl(q,j,is2)
                     c2 =-aimag( fqwei_p_noncl(ia,is_tmp) ) &
                          &   *bpr_t_noncl(p,j,is1) *bpr_t_noncl(q,j,is2)

                     fr = fr +c1 +c2

                  End do
               End do
            end do
         else
            do ia = 1, nac_p
               p = nlmta1_p(ia);     q = nlmta2_p(ia)
               Do is1=1, ndim_spinor
                  Do is2=1, ndim_spinor
                     is_tmp = ( is1 -1 )*ndim_spinor + is2

                     c1 = real( fqwei_p_noncl(ia,is_tmp) ) &
                          &   *( bpr_t_noncl(p,j,is1) *bpr_t_noncl(q,j,is2) &
                          &    + bpi_t_noncl(p,j,is1) *bpi_t_noncl(q,j,is2) )

                     c2 =-aimag( fqwei_p_noncl(ia,is_tmp) ) &
                          &   *( bpr_t_noncl(p,j,is1) *bpi_t_noncl(q,j,is2) &
                          &    + bpi_t_noncl(p,j,is1) *bpr_t_noncl(q,j,is2) )

                     fr = fr + c1 + c2

                  End do
               End do
            end do
         end if
      end if

      fr1 = 0.d0
      if ( kimg == 1 ) then

         do i = 1, np_g1k(ik)                       ! MPI
            Do is1=1, ndim_spinor
               fr1 = fr1 + psi_t_noncl(i,j,1,is1) *psi_t_noncl(i,j,1,is1)
            End Do
         end do

      else if(kimg == 2) then
!!$         ig1 = 1; if(k_symmetry(ik) == GAMMA .and. mype == 0) ig1 = 2
         ig1 = 1; if(k_symmetry(ik) == GAMMA .and. myrank_e == 0) ig1 = 2
         do i = ig1, np_g1k(ik)
            Do is1=1, ndim_spinor
               fr1 = fr1 + psi_t_noncl(i,j,1,is1) *psi_t_noncl(i,j,1,is1) &
                    &    + psi_t_noncl(i,j,2,is1) *psi_t_noncl(i,j,2,is1)
            End do
         end do
         if (k_symmetry(ik) == GAMMA) fr1 = fr1*2.d0
         if (ig1 == 2) then
            Do is1=1, ndim_spinor
               fr1 = fr1 + (psi_t_noncl(1,j,1,is1) *psi_t_noncl(1,j,1,is1) )
            End Do
         end if
      endif
      fr = fr+fr1

      if ( npes > 1 ) then
         call mpi_allreduce( fr, fr2, 1, mpi_double_precision, mpi_sum, &
	          &          mpi_k_world(myrank_k),ierr )
         fr = fr2
      end if
      fr = 1.d0/dsqrt(fr)
      if(ipri >= 2 .and. (ik == 1 .or. ik == 2)) then
         write(nfout,'(" ((WSW_t)) ik = ",i8," fr = ",f21.15)') ik, fr
      end if
    end subroutine WSW_t_noncl

#ifdef MGS_DGEMM
    subroutine cp_bpr2bprtw_noncl(i1,i2)
      integer, intent(in) :: i1,i2
      integer :: i, ia, p, q, i_last
!!$      if(i1==1) then
!!$         i_last = i2
!!$      else
!!$         i_last = i1
!!$      end if
#ifdef MGS_DGEMM_DEBUG
      i_last = i1
#else
      i_last = i2
#endif
      if(k_symmetry(ik) == GAMMA) then
         do i = i1, i_last
            do ia = 1, nac_p
               p = nlmta1_p(ia);     q = nlmta2_p(ia)
               bpr_tw1_noncl(ia,i,:) = bpr_t_noncl(p,i,:)
               bpr_tw2_noncl(ia,i,:) = bpr_t_noncl(q,i,:)
            end do
         end do
      else
         do i = i1, i_last
            do ia = 1, nac_p
               p = nlmta1_p(ia);     q = nlmta2_p(ia)
               bpr_tw1_noncl(ia,i,:) = bpr_t_noncl(p,i,:)
               bpr_tw2_noncl(ia,i,:) = bpr_t_noncl(q,i,:)
               bpi_tw1_noncl(ia,i,:) = bpi_t_noncl(p,i,:)
               bpi_tw2_noncl(ia,i,:) = bpi_t_noncl(q,i,:)
            end do
         end do
      end if
    end subroutine cp_bpr2bprtw_noncl

#endif

    subroutine normalize_bp_and_psi_t_noncl(ibo,fr)
      integer, intent(in)        :: ibo
      real(kind=DP),intent(in)   :: fr

      integer                    :: ia,ri

      if(mod_pot == VANDERBILT_TYPE) then
         if(k_symmetry(ik) == GAMMA) then
            do ia = 1, np_fs
               bpr_t_noncl(ia,ibo,:) = fr*bpr_t_noncl(ia,ibo,:)
            end do
         else
            do ia = 1, np_fs
               bpr_t_noncl(ia,ibo,:) = fr*bpr_t_noncl(ia,ibo,:)
               bpi_t_noncl(ia,ibo,:) = fr*bpi_t_noncl(ia,ibo,:)
            end do
         end if
      end if

      do ri = 1, kimg
         psi_t_noncl(1:np_g1k(ik),ibo,ri,:) &
              &  = fr * psi_t_noncl(1:np_g1k(ik),ibo,ri,:)       ! MPI
      end do

     if(ipri >= 2 .and. ik == 1 .and. ndim_spinor ==2  ) then
         write(nfout,'(" ((normalize_bp_and_psi_t_noncl)) up   ik = ",i8," ibo = ",i8," fr = ",f21.15)') ik,ibo,fr

         write(nfout,'(6d14.6)') (psi_t_noncl(ia,ibo,1,1),ia=1, 6)
         if ( kimg == 2 ) write(nfout,'(3d20.10)') (psi_t_noncl(ia,ibo,kimg,1),ia=1, 6)

         write(nfout,'(" ((normalize_bp_and_psi_t_noncl)) down ik = ",i8," ibo = ",i8," fr = ",f21.15)') ik,ibo,fr
         write(nfout,'(6d14.6)') (psi_t_noncl(ia,ibo,1,2),ia=1, 6)

      end if
    end subroutine normalize_bp_and_psi_t_noncl


!!$    subroutine alloc_p1Sp2(n)
!!$      integer, intent(in) :: n
!!$      if(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2) then
!!$         kimg_t = 1
!!$      else
!!$         kimg_t = 2
!!$      end if
!!$      allocate(p1Sp2(n,kimg_t))
!!$    end subroutine alloc_p1Sp2
!!$
!!$    subroutine dealloc_p1Sp2
!!$      deallocate(p1Sp2)
!!$    end subroutine dealloc_p1Sp2

#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT

    subroutine W1SW2_t_r_noncl(i,NB_END)
      integer, intent(in)        :: i, NB_END
! Coded by T. Yamasaki in April 2006
! Revised according to the RIKEN phase tuning project 2009, 13 Sep 2009

      integer       :: j, ia, jto, p, q,  kimg_t
      real(kind=DP) :: ar, ai
      real(kind=DP) :: c1, c2, c3, c4

      real(kind=DP), allocatable, dimension(:,:)  :: p1Sp2_t2, p1Sp2_t1  ! MPI
      character*4 F_RSVTASK
      integer       :: nt, mpant, mmdnt, ipar, ist, ied, mm
#ifndef SX
      real(kind=DP), allocatable, dimension(:,:) :: psi_ir, psi_ii
      integer :: n_unroll, jmax, ia_start
#endif

      integer :: is1, is2, is_tmp
      integer :: id_sname = -1, id_sname2 = -1

      if (sw_timing_2ndlevel == ON) call tstatc0_begin('W1SW2_t_r_noncl ',id_sname)

      if (nrank_e > 1) then
         if((k_symmetry(ik) == GAMMA .and. kimg==2) .or. kimg == 1) then
            kimg_t = 1
         else
            kimg_t = 2
         end if
         p = NB_END-i
         allocate(p1Sp2_t2(p,kimg_t)); p1Sp2_t2 = 0.d0                  ! MPI
         allocate(p1Sp2_t1(p,kimg_t)); p1Sp2_t1 = 0.d0                  ! MPI
      end if

      p1Sp2 = 0.d0

      if(mod_pot == VANDERBILT_TYPE) then
         if(kimg == 1) then
#ifndef SX
! NEC tune ------------------------------->
           if(ncache.eq.0) then
             ibsize=nac_p
           else
             ibsize=ncache/(8*(NB_END*2+3))
           endif
           do ibl1=1,nac_p,ibsize
            ibl2=ibl1+ibsize-1
            if(ibl2.gt.nac_p) ibl2=nac_p
! NEC tune <-------------------------------
#endif
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(ia,ar,ai,p,q)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
#ifndef ES
!CDIR OUTERUNROLL=4
#endif
#endif
            do j = i+1, NB_END
#ifdef SX
! NEC tune ------------------------------->
               ar = 0.d0
               do ia = 1, nac_p
                  p = nlmta1_p(ia);            q = nlmta2_p(ia)
                  Do is1=1, ndim_spinor
                     Do is2=1, ndim_spinor
                        is_tmp = ( is1 -1 )*ndim_spinor + is2

                        c1 = real( fqwei_p_noncl(ia,is_tmp) ) &
                             &   *( bpr_t_noncl(p,i,is1) *bpr_t_noncl(q,j,is2) &
                             &    + bpi_t_noncl(p,i,is1) *bpi_t_noncl(q,j,is2) )

                        c2 =-aimag( fqwei_p_noncl(ia,is_tmp) ) &
                             &   *( bpr_t_noncl(p,i,is1) *bpi_t_noncl(q,j,is2) &
                             &    + bpi_t_noncl(p,i,is1) *bpr_t_noncl(q,j,is2) )

                        ar = ar + c1 + c2

                     End Do
                  End do
               end do
               p1Sp2(j,1) = ar
#else
! NEC tune <-------------------------------
               ar = p1Sp2(j,1)
               do ia = ibl1, ibl2
                  p = nlmta1_p(ia);            q = nlmta2_p(ia)
                  Do is1=1, ndim_spinor
                     Do is2=1, ndim_spinor
                        is_tmp = ( is1 -1 )*ndim_spinor + is2

                        c1 = real( fqwei_p_noncl(ia,is_tmp) ) &
                             &   *( bpr_t_noncl(p,i,is1) *bpr_t_noncl(q,j,is2) &
                             &    + bpi_t_noncl(p,i,is1) *bpi_t_noncl(q,j,is2) )

                        c2 =-aimag( fqwei_p_noncl(ia,is_tmp) ) &
                             &   *( bpr_t_noncl(p,i,is1) *bpi_t_noncl(q,j,is2) &
                             &    + bpi_t_noncl(p,i,is1) *bpr_t_noncl(q,j,is2) )

                        ar = ar + c1 + c2
                     End do
                  End do
               end do
               p1Sp2(j,1) = ar

#endif
            end do
#ifndef SX
! NEC tune
           end do
#endif
         else if(kimg == 2) then
            if(k_symmetry(ik) == GAMMA) then
!!$#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO private(ia,ar,p,q)
!!$#endif
!!$#ifdef VPP
!!$*vocl loop, unroll(4)
!!$#endif
!!$!CDIR OUTERUNROLL=4
!!$               do j = i+1, NB_END                                                ! MPI
!!$                  ar = 0.d0
!!$                  do ia = 1, nac_p
!!$                     p = nlmta1_p(ia);         q = nlmta2_p(ia)
!!$                     ar = ar + fqwei_p(ia)*(bpr_t(p,i)*bpr_t(q,j))
!!$                  end do
!!$                  p1Sp2(j,1) = ar
!!$               end do
#ifdef NEC_TUNE_SMP
               call getenv('F_RSVTASK',F_RSVTASK)
               read (F_RSVTASK,'(i4)') nt
#else
               nt = 1
#endif
               mpant = (NB_END-i)/nt
               mmdnt = mod(NB_END-i,nt)

#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(mm,ist,ied,ia,j,p,q)
#endif
! NEC tune
!!#ifdef NEC_TUNE_SMP
#ifdef SX
               do ipar = 1, min(nt,NB_END-i)                                            ! SX
                  IF (IPAR.LE.MMDNT) THEN                                               ! SX
                     MM = MPANT+1                                                       ! SX
                  ELSE                                                                  ! SX
                     MM = MPANT                                                         ! SX
                  ENDIF                                                                 ! SX
                  IST = (IPAR-1)*MPANT + MIN(MMDNT+1,IPAR) + i                          ! SX
                  IED = IST + mm - 1                                                    ! SX
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                  do ia = 1, nac_p                                                      ! SX
                     do j = IST, IED
		     ! SX
#ifdef MGS_DGEMM


                        Do is1=1, ndim_spinor
                           Do is2=1, ndim_spinor
                              is_tmp = ( is1 -1 )*ndim_spinor + is2

                              c1 = real( fqwei_p_noncl(ia,is_tmp) ) &
                                   &   *bpr_tw1_noncl(ia,i,is1) *bpr_tw2_noncl(ia,j,is2)
                              c2 =-aimag( fqwei_p_noncl(ia,is_tmp) ) &
                                   &   *bpr_tw1_noncl(ia,i,is1) *bpr_tw2_noncl(ia,j,is2)

                              p1Sp2(j,1) = p1Sp2(j,1) + c1 +c2


                  ! SX
                           End Do
                        End Do
#else
                        p = nlmta1_p(ia);    q = nlmta2_p(ia)         ! SX
                        Do is1=1, ndim_spinor
                           Do is2=1, ndim_spinor
                              is_tmp = ( is1 -1 )*ndim_spinor + is2

                              c1 = real( fqwei_p_noncl(ia,is_tmp) ) &
                                   &   *bpr_t_noncl(p,i,is1) *bpr_t_noncl(q,j,is2)
                              c2 =-aimag( fqwei_p_noncl(ia,is_tmp) ) &
                                   &   *bpr_t_noncl(p,i,is1) *bpr_t_noncl(q,j,is2)

                              p1Sp2(j,1) = p1Sp2(j,1) + c1 +c2

                           End do
                        End do
#endif
                     end do                                       ! SX
                  end do                                          ! SX
               end do

#else
! NEC tune ------------------------------->
              if(ncache.eq.0) then
                ibsize=nac_p
              else
                ibsize=ncache/(8*(NB_END*1+2))
              endif
              do ibl1=1,nac_p,ibsize
                ibl2=ibl1+ibsize-1
                if(ibl2.gt.nac_p) ibl2=nac_p
                do j = i+1, NB_END
                  if(ibl1.eq.1)then
                    ar=0.0d0
                  else
                    ar=p1Sp2(j,1)
                  endif
                  do ia = ibl1, ibl2
#ifdef MGS_DGEMM
                    Do is1=1, ndim_spinor
                        Do is2=1, ndim_spinor
                           is_tmp = ( is1 -1 )*ndim_spinor + is2

                           c1 = real( fqwei_p_noncl(ia,is_tmp) ) &
                                &   *bpr_tw1_noncl(ia,i,is1) *bpr_tw2_noncl(ia,j,is2)
                           c2 =-aimag( fqwei_p_noncl(ia,is_tmp) ) &
                                &   *bpr_tw1_noncl(ia,i,is1) *bpr_tw2_noncl(ia,j,is2)

                           ar = ar + c1 +c2

                        End Do
                     End do
#else
                     p = nlmta1_p(ia);         q = nlmta2_p(ia)
                     Do is1=1, ndim_spinor
                        Do is2=1, ndim_spinor
                           is_tmp = ( is1 -1 )*ndim_spinor + is2

                           c1 = real( fqwei_p_noncl(ia,is_tmp) ) &
                                &   *bpr_t_noncl(p,i,is1) *bpr_t_noncl(q,j,is2)
                           c2 =-aimag( fqwei_p_noncl(ia,is_tmp) ) &
                                &   *bpr_t_noncl(p,i,is1) *bpr_t_noncl(q,j,is2)

                           ar = ar + c1 +c2

                        End Do
                     End do
#endif
                  end do
                  p1Sp2(j,1)=ar
                end do
             end do

! NEC tune <-------------------------------
#endif
            else      ! kimg==2 .and. k_symmetry(ik) /= GAMMA
#ifndef SX
! NEC tune ------------------------------->
              if(ncache.eq.0) then
                ibsize=nac_p
              else
                ibsize=ncache/(8*(NB_END*4+3))
              endif
#ifdef DEBUG_MGS
              if(ipri >= 1) then
                 write(nfout,'(" ibsize = ",i8," <<W1SW2_t_r>>")') ibsize
              end if
#endif
              do ibl1=1,nac_p,ibsize
               ibl2=ibl1+ibsize-1
               if(ibl2.gt.nac_p) ibl2=nac_p
! NEC tune <-------------------------------
#endif

#ifdef NEC_TUNE_SMP
#ifdef MGS_DGEMM
!CDIR PARALLEL DO private(ia,ar,ai)
#else
!CDIR PARALLEL DO private(ia,ar,ai,p,q)
#endif
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
#ifndef ES
!CDIR OUTERUNROLL=4
#endif
#endif
               do j = i+1, NB_END

#ifdef SX
                  ar = 0.d0; ai = 0.d0
                  do ia = 1, nac_p
#else
                  if(ibl1.eq.1) then
                     ar = 0.d0; ai = 0.d0
                  else
                    ar = p1Sp2(j,1)
                    ai = p1Sp2(j,2)

                  end if
                  do ia = ibl1, ibl2
#endif
#ifdef MGS_DGEMM
                     Do is1=1, ndim_spinor
                        Do is2=1, ndim_spinor
                           is_tmp = ( is1 -1 )*ndim_spinor +is2

                           c1 = real( fqwei_p_noncl(ia,is_tmp) ) &
                                &   *( bpr_tw1_noncl(ia,i,is1) *bpr_tw2_noncl(ia,j,is2) &
                                &    + bpi_tw1_noncl(ia,i,is1) *bpi_tw2_noncl(ia,j,is2) )
                           c2 =-aimag( fqwei_p_noncl(ia,is_tmp) ) &
                                &   *( bpr_tw1_noncl(ia,i,is1) *bpi_tw2_noncl(ia,j,is2) &
                                &    + bpi_tw1_noncl(ia,i,is1) *bpr_tw2_noncl(ia,j,is2) )

                           ar = ar + c1 + c2

                           c3 = real( fqwei_p_noncl(ia,is_tmp) ) &
                                &   *( bpr_tw1_noncl(ia,i,is1) *bpi_tw2_noncl(ia,j,is2) &
                                &    - bpi_tw1_noncl(ia,i,is1) *bpr_tw2_noncl(ia,j,is2) )
                           c4 = aimag( fqwei_p_noncl(ia,is_tmp) ) &
                                &   *( bpr_tw1_noncl(ia,i,is1) *bpr_tw2_noncl(ia,j,is2) &
                                &    + bpi_tw1_noncl(ia,i,is1) *bpi_tw2_noncl(ia,j,is2) )

                           ai = ai + c3 + c4

                        End do
                     End do
#else
                     p = nlmta1_p(ia);         q = nlmta2_p(ia)
                     Do is1=1, ndim_spinor
                        Do is2=1, ndim_spinor
                           is_tmp = ( is1 -1 )*ndim_spinor +is2

                           c1 = real( fqwei_p_noncl(ia,is_tmp) ) &
                                &   *( bpr_t_noncl(p,i,is1) *bpr_t_noncl(q,j,is2) &
                                &    + bpi_t_noncl(p,i,is1) *bpi_t_noncl(q,j,is2) )
                           c2 =-aimag( fqwei_p_noncl(ia,is_tmp) ) &
                                &   *( bpr_t_noncl(p,i,is1) *bpi_t_noncl(q,j,is2) &
                                &    + bpi_t_noncl(p,i,is1) *bpr_t_noncl(q,j,is2) )

                           ar = ar + c1 + c2

                           c3 = real( fqwei_p_noncl(ia,is_tmp) ) &
                                &   *( bpr_t_noncl(p,i,is1) *bpi_t_noncl(q,j,is2) &
                                &    - bpi_t_noncl(p,i,is1) *bsr_t_noncl(q,j,is2) )
                           c4 = aimag( fqwei_p_noncl(ia,is_tmp) ) &
                                &   *( bpr_t_noncl(p,i,is1) *bsr_t_noncl(q,j,is2) &
                                &    + bpi_t_noncl(p,i,is1) *bsi_t_noncl(q,j,is2) )

                           ai = ai + c1 + c2

                        End do
                     End do
#endif
                  end do
                  p1Sp2(j,1) = ar;  p1Sp2(j,2) = ai
               end do
#ifndef SX
! NEC tune
              end do
#endif
            end if
         end if
      end if
! --- <Psi_i|Psi_j> ---
      if(kimg == 1) then
#ifndef SX
! NEC tune ------------------------------------------------------------->
        if(ncache.eq.0) then
          ibsize=np_g1k(ik)
        else
          ibsize=ncache/(8*(NB_END*1+1))
        endif
        do ibl1=1,np_g1k(ik),ibsize
         ibl2=ibl1+ibsize-1
         if(ibl2.gt.np_g1k(ik)) ibl2=np_g1k(ik)
! NEC tune <-------------------------------------------------------------
#endif
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(ia,ar,ai,p,q)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
#ifndef ES
!CDIR OUTERUNROLL=4
#endif
#endif
         do j = i+1, NB_END
#ifdef SX
! NEC tune ------------------------------------------------------------->
            do ia = 1, np_g1k(ik)                                           ! MPI
               Do is1=1, ndim_spinor
                  p1Sp2(j,1) = p1Sp2(j,1) &
                       & + psi_t_noncl(ia,i,1,is1)*psi_t_noncl(ia,j,1,is1 ) ! MPI
               End do
#else
            ai = p1Sp2(j,1)
            do ia = ibl1, ibl2
               Do is1=1, ndim_spinor
                  ai = ai + psi_t_noncl(ia,i,1,is1)*psi_t_noncl(ia,j,1,is1) ! MPI
               End do
! NEC tune <-------------------------------------------------------------
#endif
            end do

#ifndef SX
! NEC tune
            p1Sp2(j,1)=ai
#endif

         end do

#ifndef SX
! NEC tune
        end do
#endif

      else if(kimg == 2) then
         if(k_symmetry(ik) == GAMMA) then
!!$            if(mype /= 0) then
!!$#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO private(ia,ar,ai,p,q)
!!$#endif
!!$#ifdef VPP
!!$*vocl loop, unroll(4)
!!$#endif
!!$!CDIR OUTERUNROLL=4
!!$               do j = i+1, NB_END
!!$                  do ia = 1, np_g1k(ik)
!!$                     ar  = psi_t(ia,i,1)
!!$                     ai  = psi_t(ia,i,2)
!!$                     p1Sp2(j,1) = p1Sp2(j,1)+(ar*psi_t(ia,j,1)+ai*psi_t(ia,j,2))*2.d0 ! MPI
!!$                  end do
!!$               end do
#ifdef NEC_TUNE_SMP
            call getenv('F_RSVTASK',F_RSVTASK)
            read (F_RSVTASK,'(i4)') nt
#else
            nt = 1
#endif

! NEC tune
!!#ifdef NEC_TUNE_SMP
#ifdef SX
            if( (NB_END-i)/nt .lt. 256 ) then

               if(myrank_e /= 0) then
                  do j = i+1, NB_END
                     do ia = 1, np_g1k(ik)
                        Do is1=1, ndim_spinor
                           ar  = psi_t_noncl(ia,i,1,is1)
                           ai  = psi_t_noncl(ia,i,2,is1)
                           p1Sp2(j,1) = p1Sp2(j,1)+( ar*psi_t_noncl(ia,j,1,is1)&
                                &                   +ai*psi_t_noncl(ia,j,2,is1) )*2.d0
                        End do
                     end do
                  end do
               else if(myrank_e == 0) then
                  do j = i+1, NB_END
                     do ia = 2, np_g1k(ik)
                        Do is1=1, ndim_spinor
                           ar  = psi_t_noncl(ia,i,1,is1)
                           ai  = psi_t_noncl(ia,i,2,is1)
                           p1Sp2(j,1) = p1Sp2(j,1)+( ar*psi_t_noncl(ia,j,1,is1) &
                                &                   +ai*psi_t_noncl(ia,j,2,is1))*2.d0
                        End do
                     end do
                     Do is1=1, ndim_spinor
                        p1Sp2(j,1) = p1Sp2(j,1) + psi_t_noncl(1,i,1,is1)*psi_t_noncl(1,j,1,is1)
                     End do
                  end do
               end if

            else

               mpant = (NB_END-i)/nt
               mmdnt = mod(NB_END-i,nt)
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(mm,ist,ied,ia,j,ar,ai)
#endif
               do ipar = 1, min(nt,NB_END-i)
                  IF (IPAR.LE.MMDNT) THEN
                     MM = MPANT+1
                  ELSE
                     MM = MPANT
                  ENDIF
                  IST = (IPAR-1)*MPANT + MIN(MMDNT+1,IPAR) + i
                  IED = IST + mm - 1
                  if(myrank_e /= 0) then
                     !     write(*,*)'ipar, IST, IED=', ipar, IST, IED
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                     do ia = 1, np_g1k(ik)
                        do j = IST, IED
                           Do is1=1, ndim_spinor
                              ar  = psi_t_noncl(ia,i,1,is1)
                              ai  = psi_t_noncl(ia,i,2,is1)
                              p1Sp2(j,1) = p1Sp2(j,1) &
                                   &      +( ar *psi_t_noncl(ia,j,1,is1) &
                                   &        +ai *psi_t_noncl(ia,j,2,is1) )*2.d0
                           End do
                        end do
                     end do
                  else if(myrank_e == 0) then
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                     do ia = 2, np_g1k(ik)
                        do j = IST, IED
                           Do is1=1, ndim_spinor
                              ar  = psi_t_noncl(ia,i,1,is1)
                              ai  = psi_t_noncl(ia,i,2,is1)
                              p1Sp2(j,1) = p1Sp2(j,1) &
                                   &      +( ar *psi_t_noncl(ia,j,1,is1) &
                                   &        +ai *psi_t_noncl(ia,j,2,is1) )*2.d0
                           End do
                        end do
                     enddo
                     do j = IST, IED
                        Do is1=1, ndim_spinor
                           p1Sp2(j,1) = p1Sp2(j,1) + psi_t_noncl(1,i,1,is1) &
                                &                   *psi_t_noncl(1,j,1,is1)
                        End do
                     end do
                  end if
               end do
            end if
! NEC tune ------------------------------------>
#else
!!$            if(ncache.eq.0) then
!!$               ibsize=np_g1k(ik)
!!$            else
!!$               ibsize=ncache/(8*(NB_END*2+2))
!!$            endif

!!$            allocate(psi_ir(np_g1k(ik)))
!!$            allocate(psi_ii(np_g1k(ik)))
!!$            do ia = 1, np_g1k(ik)
!!$               psi_ir(ia) = psi_t(ia,i,1)
!!$               psi_ii(ia) = psi_t(ia,i,2)
!!$            end do
!!$
!!$            n_unroll = 4
!!$            jto = (NB_END-i)/n_unroll
!!$            jmax = i + n_unroll*jto

            if(sw_timing_2ndlevel == ON) call tstatc0_begin('W1SW2_t_r(core) ',id_sname2)
            if(myrank_e /= 0) then
               ia_start = 1
            else if(myrank_e == 0) then
               ia_start = 2
            end if

            do j = i+1, NB_END
               ar=p1Sp2(j,1)
               do ia = ia_start, np_g1k(ik)
                  Do is1=1, ndim_spinor
                     ar = ar + ( psi_t_noncl(ia,i,1,is1)*psi_t_noncl(ia,j,1,is1) &
                          &     +psi_t_noncl(ia,i,2,is1)*psi_t_noncl(ia,j,2,is1) )*2.d0
                  End do
               end do
               if(myrank_e == 0) then
                 Do is1=1, ndim_spinor
                     ar = ar + psi_t_noncl(1,i,1,is1)*psi_t_noncl(1,j,1,is1)
                  End do
               endif
               p1Sp2(j,1)=ar
            end do

!!$            do j = i+1, jmax, n_unroll
!!$               do ia = ia_start, np_g1k(ik)
!!$                  ar = psi_ir(ia)*2.d0; ai = psi_ii(ia)*2.d0
!!$                  p1Sp2(j  ,1) = p1Sp2(j  ,1) + (ar*psi_t(ia,j  ,1)+ai*psi_t(ia,j  ,2))
!!$                  p1Sp2(j+1,1) = p1Sp2(j+1,1) + (ar*psi_t(ia,j+1,1)+ai*psi_t(ia,j+1,2))
!!$                  p1Sp2(j+2,1) = p1Sp2(j+2,1) + (ar*psi_t(ia,j+2,1)+ai*psi_t(ia,j+2,2))
!!$                  p1Sp2(j+3,1) = p1Sp2(j+3,1) + (ar*psi_t(ia,j+3,1)+ai*psi_t(ia,j+3,2))
!!$               end do
!!$            end do
!!$            do j = jmax+1, neg
!!$               do ia = ia_start, np_g1k(ik)
!!$                  ar = psi_ir(ia); ai = psi_ii(ia)
!!$                  p1Sp2(j  ,1) = p1Sp2(j  ,1) + (ar*psi_t(ia,j,  1)+ai*psi_t(ia,j,  2))*2.d0
!!$               end do
!!$            end do
!!$            if(mype == 0) then
!!$               do j = i+1, neg
!!$                  p1Sp2(j,1) = p1Sp2(j,1) + psi_ir(1)*psi_t(1,j, 1)
!!$               end do
!!$            end if

            if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname2)
#endif
! NEC tune <------------------------------------
         else    ! kimg==2 .and. k_symmetry(ik) /= GAMMA
#ifdef SX
! NEC tune ------------------------------------>
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(ia,ar,ai,p,q)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
#ifndef ES
!CDIR OUTERUNROLL=4
#endif
#endif
            do j = i+1, NB_END
               do ia = 1, np_g1k(ik)                                          ! MPI
                  Do is1=1, ndim_spinor
                     ar  = psi_t_noncl(ia,i,1,is1)
                     ai  = psi_t_noncl(ia,i,2,is1)
                     p1Sp2(j,1) = p1Sp2(j,1) &
                          &      + ar *psi_t_noncl(ia,j,1,is1) &
                          &      + ai *psi_t_noncl(ia,j,2,is1) ! MPI
                     p1Sp2(j,2) = p1Sp2(j,2)&
                          &      + ar *psi_t_noncl(ia,j,2,is1) &
                          &      - ai *psi_t_noncl(ia,j,1,is1) ! MPI
                  End do
               end do
            end do
#else
!!$            if(ncache.eq.0) then
               ibsize=np_g1k(ik)
!!$            else
!!$!!$               ibsize=ncache/(8*(neg*2+2))
!!$               ibsize=ncache/(8*(NB_END*2+2))
!!$            endif
            allocate( psi_ir( np_g1k(ik),ndim_spinor ) )
            allocate( psi_ii( np_g1k(ik),ndim_spinor ) )

            do ia = 1, np_g1k(ik)
               Do is1=1, ndim_spinor
                  psi_ir(ia,is1) = psi_t_noncl(ia,i,1,is1)
                  psi_ii(ia,is1) = psi_t_noncl(ia,i,2,is1)
               End do
            end do

            do ibl1=1,np_g1k(ik),ibsize
               ibl2=ibl1+ibsize-1
               if(ibl2.gt.np_g1k(ik)) ibl2=np_g1k(ik)
               if(sw_timing_2ndlevel == ON) call tstatc0_begin('W1SW2_t_r(core) ',id_sname2)
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(ia,ar,ai,p,q)
#endif
#ifdef NEC_TUNE2
#ifndef ES
!CDIR OUTERUNROLL=4
#define _DIRECTIVE_UNROLLING_
#endif
#endif
#ifdef _DIRECTIVE_UNROLLING_
               do j = i+1, NB_END
                  ar=p1Sp2(j,1)
                  ai=p1Sp2(j,2)
                 do ia = ibl1, ibl2
                     Do is1=1, ndim_spinor
                        ar = ar + psi_ir(ia,is1) *psi_t_noncl(ia,j,1,is1)&
                             &  + psi_ii(ia,is1) *psi_t_noncl(ia,j,2,is1)  ! MPI
                        ai = ai + psi_ir(ia,is1) *psi_t_noncl(ia,j,2,is1) &
                             &  - psi_ii(ia,is1) *psi_t_noncl(ia,j,1,is1) ! MPI
                     End Do
                  end do
                  p1Sp2(j,1)=ar
                  p1Sp2(j,2)=ai
               end do
#else
               n_unroll = 4
               jto = (NB_END-i)/n_unroll
               jmax = i + n_unroll*jto
               do j = i+1, jmax, n_unroll
                  do ia = ibl1, ibl2
                     Do is1=1, ndim_spinor
                        ar = psi_ir(ia,is1); ai = psi_ii(ia,is1)
                        p1Sp2(j,1)   = p1Sp2(j,1)  &
                             &        +ar *psi_t_noncl(ia,j,1,is1) &
                             &        +ai *psi_t_noncl(ia,j,2,is1) ! MPI
                        p1Sp2(j,2)   = p1Sp2(j,2)  &
                             &        +ar *psi_t_noncl(ia,j,2,is1) &
                             &        -ai *psi_t_noncl(ia,j,1,is1) ! MPI

                        p1Sp2(j+1,1) = p1Sp2(j+1,1) &
                             &        +ar *psi_t_noncl(ia,j+1,1,is1) &
                             &        +ai *psi_t_noncl(ia,j+1,2,is1) ! MPI
                        p1Sp2(j+1,2) = p1Sp2(j+1,2) &
                             &        +ar *psi_t_noncl(ia,j+1,2,is1) &
                             &        -ai *psi_t_noncl(ia,j+1,1,is1) ! MPI
                        p1Sp2(j+2,1) = p1Sp2(j+2,1) &
                             &        +ar *psi_t_noncl(ia,j+2,1,is1) &
                             &        +ai *psi_t_noncl(ia,j+2,2,is1) ! MPI
                        p1Sp2(j+2,2) = p1Sp2(j+2,2) &
                             &        +ar *psi_t_noncl(ia,j+2,2,is1) &
                             &        -ai *psi_t_noncl(ia,j+2,1,is1) ! MPI
                        p1Sp2(j+3,1) = p1Sp2(j+3,1) &
                             &        +ar *psi_t_noncl(ia,j+3,1,is1) &
                             &        +ai *psi_t_noncl(ia,j+3,2,is1) ! MPI
                        p1Sp2(j+3,2) = p1Sp2(j+3,2) &
                             &        +ar *psi_t_noncl(ia,j+3,2,is1) &
                             &        -ai *psi_t_noncl(ia,j+3,1,is1) ! MPI
                     End do
                  end do
               end do
               do j = jmax+1, NB_END
                  do ia = 1, ibl2-ibl1+1
                     Do is1=1, ndim_spinor
                        ar = psi_ir(ia,is1); ai=psi_ii(ia,is1)
                        p1Sp2(j,1)   = p1Sp2(j,1) &
                             &        +ar *psi_t_noncl(ia,j,1,is1) &
                             &        +ai *psi_t_noncl(ia,j,2,is1)
                        p1Sp2(j,2)   = p1Sp2(j,2) &
                             &        +ar *psi_t_noncl(ia,j,2,is1) &
                             &        -ai *psi_t_noncl(ia,j,1,is1)
                     End Do
                  end do
               end do
#endif
               if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname2)
            end do
           deallocate(psi_ii, psi_ir)
! NEC tune <------------------------------------
#endif
         end if
      end if

      if(nrank_e > 1 .and. NB_END > i ) then
         do q = 1, kimg_t
            do j = 1, NB_END-i
               p = j+i
               p1Sp2_t1(j,q) = p1Sp2(p,q)
            end do
         end do
!!$         p = (neg-i)*kimg_t
         call mpi_allreduce( p1Sp2_t1, p1Sp2_t2, (NB_END-i)*kimg_t, &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k),ierr )

         do q = 1, kimg_t
            do j = 1, NB_END-i
               p = j+i
               p1Sp2(p,q) = p1Sp2_t2(j,q)
            end do
         end do
!!$         p1Sp2(:,1:kimg_t) = p1Sp2_t2(:,1:kimg_t)
      end if

      if(ipri >= 2 .and. ik == 1) then
         write(nfout,'(" <<W1SW2_t>> ik = ",i9,"  i = ",i3)') ik,i
         write(nfout,'(" (real) ",6d11.3, 99(/8x,6d11.3))') (p1Sp2(j,1),j=i+1,NB_END)
         if(kimg == 2) &
              & write(nfout,'(" (imag) ",6d11.3, 99(/8x,6d11.3))') (p1Sp2(j,kimg),j=i+1,NB_END)
      end if

     if(nrank_e > 1) deallocate(p1Sp2_t1,p1Sp2_t2)
     if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
    end subroutine W1SW2_t_r_noncl

    subroutine modify_bp_and_psi_t_r_noncl(i,NB_END)
! Coded by T. Yamasaki in April 2006
! Revised according to the RIKEN phase tuning project 2009, 13 Sep 2009
      integer, intent(in) :: i, NB_END
      integer             :: j,ia

      integer :: is1

      real(kind=DP) :: sr, si
#ifndef SX
      real(kind=DP) :: ar, ai
#endif
#ifdef MGS_DGEMM
      integer :: p, q
#endif
      integer :: id_sname = -1

      if(i == neg) return
      if(sw_timing_2ndlevel == ON) call tstatc0_begin('modify_bp_and_psi_t_r ',id_sname)

      if(mod_pot == VANDERBILT_TYPE) then
         if(kimg == 1) then
#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!!$!CDIR NOSYNC
!!$!CDIR CONCUR(BY=1)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
            do j = i+1, NB_END

#ifdef SX
               do ia = 1, np_fs
                  Do is1=1, ndim_spinor
                     bpr_t_noncl(ia,j,is1) = bpr_t_noncl(ia,j,is1) &
                          &                 -p1Sp2(j,1) *bpr_t_noncl(ia,i,is1)
                     bpi_t_noncl(ia,j,is1) = bpi_t_noncl(ia,j,is1) &
                          &                 -p1Sp2(j,1) *bpi_t_noncl(ia,i,is1)
                  End Do
               end do
#else
               ar=p1Sp2(j,1)
               do ia = 1, np_fs
                  Do is1=1, ndim_spinor
                     bpr_t_noncl(ia,j,is1) = bpr_t_noncl(ia,j,is1) &
                          &                 -ar *bpr_t_noncl(ia,i,is1)
                     bpi_t_noncl(ia,j,is1) = bpi_t_noncl(ia,j,is1) &
                          &                 -ar *bpi_t_noncl(ia,i,is1)
                  End do
               end do
#endif

            end do  ! j-loop
         else if(kimg == 2) then
            if(k_symmetry(ik) == GAMMA) then
#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!!$!CDIR NOSYNC
!!$!CDIR CONCUR(BY=1)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
               do j = i+1, NB_END
#ifdef SX
! NEC tune ---------------------------------------------------------->
                  do ia = 1, np_fs
                     Do is1=1, ndim_spinor
                        bpr_t_noncl(ia,j,is1) = bpr_t_noncl(ia,j,is1) &
                             &                 -p1Sp2(j,1) *bpr_t_noncl(ia,i,is1)
                     End do
                  end do
#else
                  ar=p1Sp2(j,1)
                  do ia = 1, np_fs
                     Do is1=1, ndim_spinor
                        bpr_t_noncl(ia,j,is1) = bpr_t_noncl(ia,j,is1) &
                             &                 -ar *bpr_t_noncl(ia,i,is1)
                     End Do
                  end do
! NEC tune <----------------------------------------------------------
#endif
               end do
            else  ! kimg==2, k_symmetry(ik) /= GAMMA
#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!!$!CDIR NOSYNC
!!$!CDIR CONCUR(BY=1)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
               do j = i+1, NB_END
#ifdef SX
! NEC tune ---------------------------------------------------------->
                 do ia = 1, np_fs
                     Do is1=1, ndim_spinor
                        sr  =  bpr_t_noncl(ia,i,is1);
                        si  =  bpi_t_noncl(ia,i,is1)
                        bpr_t_noncl(ia,j,is1) = bpr_t_noncl(ia,j,is1) &
                             &                 -p1Sp2(j,1) *sr +p1Sp2(j,2) *si
                        bpi_t_noncl(ia,j,is1) = bpi_t_noncl(ia,j,is1) &
                             &                 -p1Sp2(j,1) *si -p1Sp2(j,2) *sr
                     End do
                  end do
#else
                  ar=p1Sp2(j,1)
                  ai=p1Sp2(j,2)
                  do ia = 1, np_fs
                     Do is1=1, ndim_spinor
                        sr  =  bpr_t_noncl(ia,i,is1);
                        si  =  bpi_t_noncl(ia,i,is1)
                        bpr_t_noncl(ia,j,is1) = bpr_t_noncl(ia,j,is1) -ar*sr +ai*si
                        bpi_t_noncl(ia,j,is1) = bpi_t_noncl(ia,j,is1) -ar*si -ai*sr
                     End do
                  end do
#endif
               end do
! NEC tune <----------------------------------------------------------
            end if
         end if
      end if

      if(kimg == 1) then
#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!!$!CDIR NOSYNC
!!$!CDIR CONCUR(BY=1)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
         do j = i+1,NB_END
#ifdef SX
! NEC tune ---------------------------------------------------------->
            do ia = 1, np_g1k(ik)                 ! MPI
               Do is1=1, ndim_spinor
                  psi_t_noncl(ia,j,1,is1) = psi_t_noncl(ia,j,1,is1)&
                       &                   -p1Sp2(j,1) *psi_t_noncl(ia,i,1,is1)
               End do
            end do
#else
            ar=p1Sp2(j,1)
            do ia = 1, np_g1k(ik)
               Do is1=1, ndim_spinor
                  psi_t_noncl(ia,j,1,is1) = psi_t_noncl(ia,j,1,is1) &
                       &                   -ar *psi_t_noncl(ia,i,1,is1)
               End do
            end do
! NEC tune <----------------------------------------------------------
#endif
         end do
      else if(kimg == 2) then
         if(k_symmetry(ik) == GAMMA) then
#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!!$!CDIR NOSYNC
!!$!CDIR CONCUR(BY=1)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
            do j = i+1,NB_END
#ifdef SX
! NEC tune ---------------------------------------------------------->
              do ia = 1, np_g1k(ik)
                  Do is1=1, ndim_spinor
                     psi_t_noncl(ia,j,1,is1) = psi_t_noncl(ia,j,1,is1) &
                          &                   -p1Sp2(j,1) *psi_t_noncl(ia,i,1,is1)
                     psi_t_noncl(ia,j,2,is1) = psi_t_noncl(ia,j,2,is1) &
                          &                   -p1Sp2(j,1) *psi_t_noncl(ia,i,2,is1)
                  End do
               end do
#else
               ar=p1Sp2(j,1)
               do ia = 1, np_g1k(ik)
                  Do is1=1, ndim_spinor
                     psi_t_noncl(ia,j,1,is1) = psi_t_noncl(ia,j,1,is1) &
                          &                   -ar *psi_t_noncl(ia,i,1,is1 )
                     psi_t_noncl(ia,j,2,is1) = psi_t_noncl(ia,j,2,is1) &
                          &                   -ar *psi_t_noncl(ia,i,2,is1 )
                  End Do
               end do
! NEC tune <----------------------------------------------------------
#endif
            end do
         else   ! kimg==2 .and. k_symmetry(ik) /= GAMMA
#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!!$!CDIR NOSYNC
!!$!CDIR CONCUR(BY=1)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
            do j = i+1,NB_END
#ifdef SX
! NEC tune ---------------------------------------------------------->
               do ia = 1, np_g1k(ik)                ! MPI
                  Do is1=1, ndim_spinor
                     sr  =  psi_t_noncl(ia,i,1,is1 )
                     si  =  psi_t_noncl(ia,i,2,is1 )
                     psi_t_noncl(ia,j,1,is1) = psi_t_noncl(ia,j,1,is1) &
                          &                   -p1Sp2(j,1) *sr +p1Sp2(j,2) *si
                     psi_t_noncl(ia,j,2,is1) = psi_t_noncl(ia,j,2,is1) &
                          &                   -p1Sp2(j,1) *si -p1Sp2(j,2) *sr
                  End do
               end do
#else
               ar=p1Sp2(j,1)
               ai=p1Sp2(j,2)
               do ia = 1, np_g1k(ik)
                  Do is1=1, ndim_spinor
                     sr  =  psi_t_noncl(ia,i,1,is1)
                     si  =  psi_t_noncl(ia,i,2,is1)
                     psi_t_noncl(ia,j,1,is1) = psi_t_noncl(ia,j,1,is1) -ar*sr +ai*si
                     psi_t_noncl(ia,j,2,is1) = psi_t_noncl(ia,j,2,is1) -ar*si -ai*sr
                  End do
               end do
#endif
! NEC tune <----------------------------------------------------------
            end do
         end if
      end if
      if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
    end subroutine modify_bp_and_psi_t_r_noncl

#ifdef MGS_DGEMM
!-$ tune FUJITSU for block ----->>
!     13th Sep 2009

    subroutine W1SW2_t_r_block_noncl(i,i2,p1Sp2_NB,NB_END,NB_END2,kimg_t_wk)
      integer, intent(in)        :: i, i2, NB_END, NB_END2, kimg_t_wk

      real(kind=DP), intent(out), dimension(NB,NB,kimg_t_wk) :: p1Sp2_NB

      real(kind=DP), allocatable :: bpr_tw1_noncl_BLAS(:,:,:)
      real(kind=DP), allocatable :: bpi_tw1_noncl_BLAS(:,:,:)
!
      real(kind=DP), allocatable, dimension(:,:,:)  ::  p1Sp2_t1_NB  ! MPI
      integer       :: iblk, jblk, i_NB, j_NB
      integer       :: ii, ia,  M, N, LDA, LDB, p
      integer       :: LDC

      real(kind=DP) :: c1, c2, c3, c4

      integer  :: is1, is2, is_tmp

      integer :: id_sname = -1, id_sname2 = -1
      if (sw_timing_2ndlevel == ON) then
         call tstatc0_begin('W1SW2_t_r_block_noncl ',id_sname)
      endif

      p1Sp2_NB = 0.0d0
      M = NB_END2-i2+1; N = NB_END-i+1

      LDC = NB

      if(mod_pot == VANDERBILT_TYPE .and. nac_p > 0) then
         allocate( bpr_tw1_noncl_BLAS(nac_p, NB, ndim_spinor ) )
         bpr_tw1_noncl_BLAS = 0.0

         do iblk = i, NB_END
            ii = iblk - i + 1
            do ia = 1, nac_p
               Do is1=1, ndim_spinor
                  Do is2=1, ndim_spinor
                     is_tmp = ( is1 -1 ) *ndim_spinor + is2

                     c1 = real( fqwei_p_noncl(ia,is_tmp) ) &
                          &   * bpr_tw1_noncl(ia,iblk,is2)

                     if ( allocated( bpi_tw1_noncl ) ) then
                        c2 =-aimag( fqwei_p_noncl(ia,is_tmp) ) &
                             &   * bpi_tw1_noncl(ia,iblk,is2)
                     else
                        c2 = 0.0d0
                     endif

                     bpr_tw1_noncl_BLAS(ia,ii,is1) = bpr_tw1_noncl_BLAS(ia,ii,is1) &
                          &                         + c1 + c2

                  End do
               End do
            end do
         end do
         if (k_symmetry(ik) /= GAMMA) then
            allocate(bpi_tw1_noncl_BLAS(nac_p,NB, ndim_spinor ))
            bpi_tw1_noncl_BLAS = 0.d0
            do iblk = i, NB_END
               ii = iblk - i + 1
               do ia = 1, nac_p

                  Do is1=1, ndim_spinor
                     Do is2=1, ndim_spinor
                        is_tmp = ( is1 -1 ) *ndim_spinor + is2

                        c3 = real( fqwei_p_noncl(ia,is_tmp) ) &
                             &   * bpi_tw1_noncl(ia,iblk,is2)

                        c4 = aimag( fqwei_p_noncl(ia,is_tmp) ) &
                                &   * bpr_tw1_noncl(ia,iblk,is2)

                        bpi_tw1_noncl_BLAS(ia,ii,is1) = bpi_tw1_noncl_BLAS(ia,ii,is1) &
                             &                         + c3 + c4
                     End Do
                  End do
               end do

            end do
         end if

         if (kimg==1) then
            Do is1=1, ndim_spinor
               call DGEMM__( 'T', 'N', M, N, nac_p, 1.0d0, &
                    &         bpr_tw2_noncl(1,i2,is1), nac_p, &
                    &         bpr_tw1_noncl_BLAS(:,:,is1), nac_p, &
                    &         1.0d0, p1Sp2_NB(1,1,1), LDC )
               call DGEMM__( 'T', 'N', M, N, nac_p, 1.0d0, &
                    &         bpi_tw2_noncl(1,i2,is1), nac_p, &
                    &         bpi_tw1_noncl_BLAS(:,:,is1), nac_p, &
                    &         1.0d0, p1Sp2_NB(1,1,1), LDC )
            End do
!!
!!$             ar = ar + fqwei_p(ia)*(bpr_t(p,i)*bpr_t(q,j)+bpi_t(p,i)*bpi_t(q,j))
!!
         else if ( kimg==2 .and. k_symmetry(ik) == GAMMA ) then
            Do is1=1, ndim_spinor
               call DGEMM__( 'T', 'N', M, N, nac_p, 1.0d0, &
                    &         bpr_tw2_noncl(1,i2,is1), nac_p, &
                    &         bpr_tw1_noncl_BLAS(:,:,is1), nac_p, &
                    &         1.0d0, p1Sp2_NB(1,1,1), LDC )
            End do
!!
!!$                     ar = ar + fqwei_p(ia)*(bpr_tw1(ia,i)*bpr_tw2(ia,j))
!!
         else if ( kimg==2 .and. k_symmetry(ik) /= GAMMA ) then
            Do is1=1, ndim_spinor
               call DGEMM__( 'T', 'N', M, N, nac_p, 1.0d0, &
                    &         bpr_tw2_noncl(1,i2,is1), nac_p, &
                    &         bpr_tw1_noncl_BLAS(:,:,is1), nac_p, &
                    &         1.0d0, p1Sp2_NB(1,1,1), LDC )
               call DGEMM__( 'T', 'N', M, N, nac_p, 1.0d0, &
                    &         bpi_tw2_noncl(1,i2,is1), nac_p, &
                    &         bpr_tw1_noncl_BLAS(:,:,is1), nac_p, &
                    &         1.0d0, p1Sp2_NB(1,1,2), LDC )

               call DGEMM__( 'T', 'N', M, N, nac_p, 1.0d0, &
                    &         bpi_tw2_noncl(1,i2,is1), nac_p, &
                    &         bpi_tw1_noncl_BLAS(:,:,is1), nac_p, &
                    &         1.0d0, p1Sp2_NB(1,1,1), LDC )
               call DGEMM__( 'T', 'N', M, N, nac_p, -1.0d0, &
                    &         bpr_tw2_noncl(1,i2,is1), nac_p, &
                    &         bpi_tw1_noncl_BLAS(:,:,is1) ,nac_p, &
                    &         1.0d0, p1Sp2_NB(1,1,2), LDC )
!!
!!$                     ar = ar + fqwei_p(ia)*(bpr_tw1(ia,i)*bpr_tw2(ia,j)+bpi_tw1(ia,i)*bpi_tw2(ia,j))
!!$                     ai = ai + fqwei_p(ia)*(bpr_tw1(ia,i)*bpi_tw2(ia,j)-bpi_tw1(ia,i)*bpr_tw2(ia,j))
!!
            End Do
         end if

         if (k_symmetry(ik) /= GAMMA) deallocate(bpi_tw1_noncl_BLAS)
         deallocate(bpr_tw1_noncl_BLAS)
      end if
! --------------------------------------
      if (sw_timing_2ndlevel == ON) then
         call tstatc0_begin('W1SW2_t_r_block_noncl(core) ',id_sname2)
      endif
      LDA = np_g1k_x; LDB = np_g1k_x
      if ( kimg == 1 ) then
         Do is1=1, ndim_spinor
            call DGEMM__( 'T', 'N', M, N, np_g1k(ik), 1.0d0, &
                 &        psi_t_noncl(1,i2,1,is1), LDA, &
                 &        psi_t_noncl(1, i,1,is1), LDB, &
                 &        1.0d0, p1Sp2_NB(1,1,1), LDC )
         End do
      else if ( kimg == 2 .and. k_symmetry(ik) == GAMMA ) then
         Do is1=1, ndim_spinor
            call DGEMM__( 'T', 'N', M, N, np_g1k(ik), 2.0d0, &
                 &        psi_t_noncl(1,i2,1,is1), LDA, &
                 &        psi_t_noncl(1, i,1,is1) ,LDB, &
                 &        1.0d0, p1Sp2_NB(1,1,1), LDC )
            call DGEMM__( 'T', 'N', M, N, np_g1k(ik), 2.0d0, &
                 &        psi_t_noncl(1,i2,2,is1), LDA, &
                 &        psi_t_noncl(1, i,2,is1), LDB, &
                 &        1.0d0, p1Sp2_NB(1,1,1), LDC )
         End do

         if ( myrank_e == 0 ) then
            do iblk = i, NB_END                !for BLAS3
               i_NB = iblk - i + 1
               do jblk = i2, NB_END2
                  j_NB = jblk - i2 + 1
                  Do is1=1, ndim_spinor
                     p1Sp2_NB(j_NB,i_NB,1) = p1Sp2_NB(j_NB,i_NB,1) &
                          & - ( psi_t_noncl(1,iblk,1,is1)*psi_t_noncl(1,jblk,1,is1) &
                          &    +psi_t_noncl(1,iblk,2,is1)*psi_t_noncl(1,jblk,2,is1)*2.d0)
                  End Do
               end do
            end do
         end if

      else  ! kimg==2 .and. k_symmetry(ik) /= GAMMA
         Do is1=1, ndim_spinor
            call DGEMM__( 'T', 'N', M, N, np_g1k(ik), 1.0d0, &
                 &        psi_t_noncl(1,i2,1,is1), LDA, &
                 &        psi_t_noncl(1, i,1,is1), LDB, &
                 &        1.0d0, p1Sp2_NB(1,1,1), LDC )
            call DGEMM__( 'T', 'N', M, N, np_g1k(ik), 1.0d0, &
                 &        psi_t_noncl(1,i2,2,is1), LDA, &
                 &        psi_t_noncl(1, i,2,is1), LDB, &
                 &        1.0d0, p1Sp2_NB(1,1,1), LDC )
            call DGEMM__( 'T', 'N', M, N, np_g1k(ik), 1.0d0, &
                 &        psi_t_noncl(1,i2,2,is1), LDA, &
                 &        psi_t_noncl(1, i,1,is1), LDB, &
                 &        1.0d0, p1Sp2_NB(1,1,2), LDC )
            call DGEMM__( 'T', 'N', M, N, np_g1k(ik), -1.0d0, &
                 &        psi_t_noncl(1,i2,1,is1), LDA, &
                 &        psi_t_noncl(1, i,2,is1), LDB, &
                 &        1.0d0, p1Sp2_NB(1,1,2), LDC )
         End Do
      end if

      if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname2)

      if (nrank_e > 1 ) then
         if(ipri>=2) then
            write(nfout,'(" NB_END, NB_END2, NB_END-i+1, NB_END2-i2+1, NB, kimg_t_wk = ",6i8)') &
                 & NB_END, NB_END2,  NB_END-i+1, NB_END2-i2+1, NB, kimg_t_wk
         end if
#ifdef SX
         allocate(p1Sp2_t1_NB(LMB,NB,kimg_t_wk)); p1Sp2_t1_NB = 0.d0           ! MPI
         call mpi_allreduce(p1Sp2_NB, p1Sp2_t1_NB,LMB*NB*kimg_t_wk,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
#else
         allocate(p1Sp2_t1_NB(NB,NB,kimg_t_wk)); p1Sp2_t1_NB = 0.d0            ! MPI
         call mpi_allreduce(p1Sp2_NB, p1Sp2_t1_NB,NB*NB*kimg_t_wk,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
#endif
         p1Sp2_NB = p1Sp2_t1_NB
         deallocate(p1Sp2_t1_NB)
      end if

     if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)

    end subroutine W1SW2_t_r_block_noncl

!-$ tune FUJITSU for block <<-----

    subroutine modify_bp_psi_t_r_block_noncl(i,i2,p1Sp2_NB,NB_END,NB_END2,kimg_t_wk)
      integer, intent(in) :: i,i2
      real(kind=DP), intent(in),  dimension(NB,NB,kimg_t_wk) :: p1Sp2_NB
      integer, intent(in) :: NB_END, NB_END2, kimg_t_wk
      integer       ::  N, K, LDA, LDB, LDC
      integer :: is1

      integer :: id_sname = -1
      if(i == neg) return

      if(sw_timing_2ndlevel == ON) then
         call tstatc0_begin('modify_bp_psi_t_r_block_noncl ',id_sname)
      endif

      N = NB_END2-i2+1
      K = NB_END-i+1

      if(mod_pot == VANDERBILT_TYPE .and. np_fs.gt.0) then
         LDA = np_fs_x; LDB = NB; LDC = np_fs_x

         Do is1=1, ndim_spinor
            call DGEMM__( 'N', 'T', np_fs, N, K, -1.0d0, &
                 &        bpr_t_noncl(1,i,is1), LDA, p1Sp2_NB(1,1,1), LDB, 1.0d0, &
                 &        bpr_t_noncl(1,i2,is1),LDC )
         End do

         if(kimg == 1 .or. (kimg==2 .and. k_symmetry(ik) /= GAMMA)) then
            Do is1=1, ndim_spinor
               call DGEMM__( 'N', 'T', np_fs, N, K, -1.0d0, &
                    &        bpi_t_noncl(1,i,is1), LDA, p1Sp2_NB(1,1,1), LDB, 1.0d0, &
                    &        bpi_t_noncl(1,i2,is1),LDC)
            End do
         end if

         if(kimg==2 .and. k_symmetry(ik) /= GAMMA) then
            Do is1=1, ndim_spinor
               call DGEMM__( 'N', 'T', np_fs, N, K, 1.0d0, &
                    &         bpi_t_noncl(1,i,is1), LDA, p1Sp2_NB(1,1,2), LDB, 1.0d0,&
                    &         bpr_t_noncl(1,i2,is1), LDC )
               call DGEMM__( 'N', 'T', np_fs, N, K, -1.0d0,&
                    &         bpr_t_noncl(1,i,is1), LDA, p1Sp2_NB(1,1,2), LDB, 1.0d0,&
                    &         bpi_t_noncl(1,i2,is1),LDC )
            End do
         end if
      end if

      LDA = np_g1k_x; LDB = NB; LDC = np_g1k_x
      if(kimg == 1) then
         Do is1=1, ndim_spinor
            call DGEMM__( 'N', 'T', np_g1k(ik), N, K, -1.0d0, &
                 &        psi_t_noncl(1,i,1,is1), LDA, p1Sp2_NB(1,1,1), LDB, 1.0d0,&
                 &        psi_t_noncl(1,i2,1,is1),LDC )
         End do
      else if(kimg==2 .and. k_symmetry(ik) == GAMMA) then
         Do is1=1, ndim_spinor
            call DGEMM__( 'N', 'T', np_g1k(ik), N, K, -1.0d0, &
                 &        psi_t_noncl(1,i,1,is1), LDA, p1Sp2_NB(1,1,1), LDB,1.0d0, &
                 &        psi_t_noncl(1,i2,1,is1),LDC )
            call DGEMM__( 'N', 'T', np_g1k(ik), N, K, -1.0d0, &
                 &        psi_t_noncl(1,i,2,is1), LDA, p1Sp2_NB(1,1,1), LDB, 1.0d0,&
                 &        psi_t_noncl(1,i2,2,is1),LDC )
         End do
      else  ! kimg == 2 .and. k_symmetry(ik) /= GAMMA
         Do is1=1, ndim_spinor
            call DGEMM__( 'N', 'T', np_g1k(ik), N, K, -1.0d0, &
                 &        psi_t_noncl(1,i,1,is1), LDA, p1Sp2_NB(1,1,1),LDB,1.0d0, &
                 &        psi_t_noncl(1,i2,1,is1),LDC )
            call DGEMM__( 'N', 'T', np_g1k(ik), N, K, -1.0d0,&
                 &         psi_t_noncl(1,i,2,is1), LDA, p1Sp2_NB(1,1,1),LDB,1.0d0,&
                 &         psi_t_noncl(1,i2,2,is1),LDC )
            call DGEMM__( 'N', 'T', np_g1k(ik), N, K, 1.0d0, &
                 &         psi_t_noncl(1,i,2,is1), LDA, p1Sp2_NB(1,1,2),LDB,1.0d0,&
                 &         psi_t_noncl(1,i2,1,is1),LDC )
            call DGEMM__( 'N', 'T', np_g1k(ik), N, K, -1.0d0,&
                 &         psi_t_noncl(1,i,1,is1), LDA, p1Sp2_NB(1,1,2),LDB,1.0d0,&
                 &         psi_t_noncl(1,i2,2,is1),LDC )
         End do
      end if

      if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)

    end subroutine modify_bp_psi_t_r_block_noncl



#endif

#else

    subroutine substitute_jto_ib2back(i,nmax)
      integer, intent(in)  :: i
      integer, intent(out) :: nmax

      integer              ::jto, j
      jto = 0
      do j = 1, neg
         if(nrvf_ordr(j,ik) <= i) cycle
         jto = jto + 1
         ib2to_a(jto)  = j
         ib2back_a(j)  = jto
      end do
      nmax = jto
      if(ipri >= 2 .and. ik == 1) then
         write(nfout,'(" <<substitute_jto_ib2back>>")')
         do j = 1, neg
            write(nfout,'(" ib2back_a(",i3,") = ",i3)') j, ib2back_a(j)
         end do
         do j = 1, nmax
            write(nfout,'(" ib2to_a(",i3,") = ",i3)') j, ib2to_a(j)
         end do
      end if
    end subroutine substitute_jto_ib2back

    subroutine W1SW2_t_noncl(i,ito)
      integer, intent(in)        :: i,ito

      integer       :: j, ia, jto, p, q, kimg_t
      integer :: is1, is2, is_tmp

      real(kind=DP) :: ar, ai
      real(kind=DP), allocatable, dimension(:,:)  :: p1Sp2_t1, p1Sp2_t2             ! MPI

!!!$      call alloc_p1Sp2(nmax)
      if(npes > 1) then
         if((k_symmetry(ik) == GAMMA .and. kimg==2) .or. kimg == 1) then
            kimg_t = 1
         else
            kimg_t = 2
         end if
         allocate(p1Sp2_t2(nmax,kimg_t)); p1Sp2_t2 = 0.d0                  ! MPI
         allocate(p1Sp2_t1(nmax,kimg_t)); p1Sp2_t1 = 0.d0                  ! MPI
      end if

      p1Sp2 = 0.d0

      if(mod_pot == VANDERBILT_TYPE) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(jto,ia,ar,ai,p,q)
#endif
         do jto = 1, nmax                                                 ! MPI
            j = ib2to_a(jto)                                              ! MPI
            if(kimg == 1) then
               ar = 0.d0
               do ia = 1, nac_p
                  p = nlmta1_p(ia);            q = nlmta2_p(ia)
                  Do is1=1, ndim_spinor
                     Do is2=1, ndim_spinor
                        is_tmp =  ( is1 - 1)*ndim_spinor + is2
                        ar = ar + fqwei_p_eff(ia,is_tmp) &
                                 & *( bpr_t_noncl(p,ito,is1) *bpr_t_noncl(q,j,is2) &
                                 &   +bpi_t_noncl(p,ito,is1) *bpi_t_noncl(q,j,is2) )
                     End Do
                  End Do
               end do
               p1Sp2(jto,1) = ar
            else if(kimg == 2) then
               if(k_symmetry(ik) == GAMMA) then
                  ar = 0.d0
                  do ia = 1, nac_p
                     p = nlmta1_p(ia);         q = nlmta2_p(ia)
                     Do is1=1, ndim_spinor
                        Do is2=1, ndim_spinor
                           is_tmp =  ( is1 - 1)*ndim_spinor + is2
                           ar = ar + fqwei_eff(ia,is_tmp) &
                                &   *( bpr_t_noncl(p,ito,is1) *bpr_t_noncl(q,j,is2))
                        End do
                     End Do
                  end do
                  p1Sp2(jto,1) = ar
               else
                  ar = 0.d0; ai = 0.d0
                  do ia = 1, nac_p
                     p = nlmta1_p(ia);         q = nlmta2_p(ia)
                     Do is1=1, ndim_spinor
                        Do is2=1, ndim_spinor
                           is_tmp =  ( is1 - 1)*ndim_spinor + is2
                           ar = ar + fqwei_eff(ia,is_tmp) &
                                &   *( bpr_t_noncl(p,ito,is1) *bpr_t_noncl(q,j,is2) &
                                &     +bpi_t_noncl(p,ito,is1) *bpi_t_noncl(q,j,is2) )
                           ai = ai + fqwei_eff(ia,is_tmp) &
                                &   *( bpr_t_noncl(p,ito,is1) *bpi_t_noncl(q,j,is2) &
                                &     -bpi_t_noncl(p,ito,is1) *bpr_t_noncl(q,j,is2) )
                        End Do
                     End Do
                  end do
                  p1Sp2(jto,1) = ar;  p1Sp2(jto,2) = ai
               end if
            end if
         end do
      end if

      if(ipri >= 2 .and. ik == 1) then
         if(nmax > 1) then
            write(nfout,'(" <<W1SW2_t>> ik = ",i9,"  i = ",i3, " ito = ", i3)') ik,i,ito
            write(nfout,'(" (real) ",6d11.3, 99(/8x,6d11.3))') (p1Sp2(jto,1),jto=1, nmax)
            if(kimg == 2 .and. k_symmetry(ik) /= GAMMA) &
                 & write(nfout,'(" (imag) ",6d11.3, 99(/8x,6d11.3))') (p1Sp2(jto,kimg),jto=1, nmax)
         end if
      end if

     if(kimg == 1) then
         do jto = 1, nmax
            j = ib2to_a(jto)
            do ia = 1, np_g1k(ik)                                           ! MPI
               Do is1=1, ndim_spinor
                  p1Sp2(jto,1) = p1Sp2(jto,1) &
                       &        + psi_t_noncl(ia,ito,1,is1) &
                       &         *psi_t_noncl(ia,  j,1,is1) ! MPI
               End do
            end do
         end do
      else if(kimg == 2) then
         if(k_symmetry(ik) == GAMMA) then
            do jto = 1, nmax
               j = ib2to_a(jto)
               do ia = 2, np_g1k(ik)
                  Do is1=1, ndim_spinor
                     ar  = psi_t_noncl(ia,ito,1,is1)
                     ai  = psi_t_noncl(ia,ito,2,is1)
                     p1Sp2(jto,1) = p1Sp2(jto,1) &
                          &        + ( ar *psi_t_noncl(ia,j,1,is1) &
                          &           +ai *psi_t_noncl(ia,j,2,is1) )*2.d0 ! MPI
                  End do
               end do
!!!$               if(mype /= 0) then
               if(myrank_e /= 0) then
                  Do is1=1, ndim_spinor
                     ar = psi_t_noncl(1,ito,1,is1)
                     ai = psi_t_noncl(1,ito,2,is1)
                     p1Sp2(jto,1) = p1Sp2(jto,1) &
                          &        + ( ar *psi_t_noncl(1,j,1,is1) &
                          &           +ai *psi_t_noncl(1,j,2,is1) )*2.d0
                  End  do
               else
                  Do is1=1, ndim_spinor
                     p1Sp2(jto,1) = p1Sp2(jto,1) &
                          &        + psi_t_noncl(1,ito,1,is1)*psi_t_noncl(1,j,1,is1)
                  End do
               end if
            end do
         else
            do jto = 1, nmax
               j = ib2to_a(jto)
               do ia = 1, np_g1k(ik)                                          ! MPI
                  Do is1=1, ndim_spinor
                     ar  = psi_t(ia,ito,1,is1)
                     ai  = psi_t(ia,ito,2,is1)
                     p1Sp2(jto,1) = p1Sp2(jto,1) &
                          &        +ar *psi_t_noncl(ia,j,1,is1) &
                          &        +ai *psi_t_noncl(ia,j,2,is1) ! MPI
                     p1Sp2(jto,2) = p1Sp2(jto,2) &
                          &        +ar *psi_t_noncl(ia,j,2,is1) &
                          &        -ai *psi_t_noncl(ia,j,1,is1) ! MPI
                  End Do
               end do
            end do
         end if
      end if

      if(npes > 1 ) then
         do q = 1, kimg_t
            do jto = 1, nmax
               p1Sp2_t1(jto,q) = p1Sp2(jto,q)
            end do
         end do
         call mpi_allreduce( p1Sp2_t1, p1Sp2_t2, nmax*kimg_t, mpi_double_precision, &
              &              mpi_sum, mpi_k_world(myrank_k), ierr )
         do q = 1, kimg_t
            do jto = 1, nmax
               p1Sp2(jto,q) = p1Sp2_t2(jto,q)
            end do
         end do
      end if

      if(ipri >= 2 .and. ik == 1) then
         if(nmax > 1) then
            write(nfout,'(" <<W1SW2_t>> ik = ",i9,"  i = ",i3, " ito = ", i3)') ik,i,ito
            write(nfout,'(" (real) ",6d11.3, 99(/8x,6d11.3))') (p1Sp2(jto,1),jto=1, nmax)
            if(kimg == 2 .and. k_symmetry(ik) /= GAMMA) &
                 & write(nfout,'(" (imag) ",6d11.3, 99(/8x,6d11.3))') &
                 &                     (p1Sp2(jto,kimg),jto=1, nmax)
         end if
      end if

      if(npes > 1) deallocate(p1Sp2_t1,p1Sp2_t2)
    end subroutine W1SW2_t_noncl

    subroutine modify_bp_and_psi_t_noncl(i,ito)
      integer, intent(in) :: i, ito
      integer             :: j,ia,jto
      real(kind=DP) :: sr, si

      if(mod_pot == VANDERBILT_TYPE) then
         do j = 1, neg
            if(nrvf_ordr(j,ik) <= i) cycle
            jto = ib2back_a(j)
            if(kimg == 1) then
               do ia = 1, np_fs
                  Do is1=1, ndim_spinor
                     bpr_t_noncl(ia,j,is1) = bpr_t_noncl(ia,j,is1) &
                          &                 -p1Sp2(jto,1) *bpr_t_noncl(ia,ito,is1)
                     bpi_t_noncl(ia,j,is1) = bpi_t_noncl(ia,j,is1) &
                          &                 -p1Sp2(jto,1) *bpi_t_noncl(ia,ito,is1)
                  Enddo
               end do
            else if(kimg == 2) then
               if(k_symmetry(ik) == GAMMA) then
                  do ia = 1, np_fs
                     Do is1=1, ndim_spinor
                        bpr_t_noncl(ia,j,is1) = bpr_t_noncl(ia,j,is1) &
                             &                 -p1Sp2(jto,1) *bpr_t_noncl(ia,ito,is1)
                     End do
                  end do
               else
                  do ia = 1, np_fs
                     Do is1=1, ndim_spinor
                        sr  =  bpr_t_noncl(ia,ito,is1)
                        si  =  bpi_t_noncl(ia,ito,is1)
                        bpr_t_noncl(ia,j,is1) = bpr_t_noncl(ia,j,is1) &
                             &                 -p1Sp2(jto,1)*sr +p1Sp2(jto,2)*si
                        bpi_t_noncl(ia,j,is1) = bpi_t_noncl(ia,j,is1) &
                             &                 -p1Sp2(jto,1)*si -p1Sp2(jto,2)*sr
                     End do
                  end do
               end if
            end if
         end do
      end if

#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!CDIR NOSYNC
!CDIR CONCUR(BY=1)
#endif
      do j = 1,neg
         if(nrvf_ordr(j,ik) <= i) cycle
         jto = ib2back_a(j)         ! <-- 99Jan19 T. Y.
         if(kimg == 1) then
            do ia = 1, np_g1k(ik)                 ! MPI
               Do is1=1, ndim_spinor
                  psi_t_noncl(ia,j,1,is1) = psi_t_noncl(ia,j,1,is1) &
                       &                   -p1Sp2(jto,1) *psi_t_noncl(ia,ito,1,is1 )
               End do
            end do
         else if(kimg == 2) then
            if(k_symmetry(ik) == GAMMA) then
               do ia = 1, np_g1k(ik)
                  Do is1=1, ndim_spinor
                     sr  =  psi_t_noncl(ia,ito,1,is1)
                     si  =  psi_t_noncl(ia,ito,2,is1)
                     psi_t_noncl(ia,j,1,is1) = psi_t_noncl(ia,j,1,is1) &
                          &                   -p1Sp2(jto,1) *sr
                     psi_t_noncl(ia,j,2,is1) = psi_t_noncl(ia,j,2,is1) &
                          &                   -p1Sp2(jto,1) *si
                  End do
               end do
            else
               do ia = 1, np_g1k(ik)                ! MPI
                  Do is1=1, ndim_spinor
                     sr  =  psi_t_noncl(ia,ito,1,is1)
                     si  =  psi_t_noncl(ia,ito,2,is1)
                     psi_t_noncl(ia,j,1,is1) = psi_t_noncl(ia,j,1,is1) &
                          &                   -p1Sp2(jto,1) *sr +p1Sp2(jto,2) *si
                     psi_t_noncl(ia,j,2,is1) = psi_t_noncl(ia,j,2,is1) &
                          &                   -p1Sp2(jto,1) *si -p1Sp2(jto,2) *sr
                  End do
               end do
            end if
         end if
      end do
!!$      call dealloc_p1Sp2()
    end subroutine modify_bp_and_psi_t_noncl
#endif

  end subroutine mgs_4_each_k_G_noncl
! ================================================================= 11.0



!!$  subroutine set_npzri(mod_pot,ik,npzr_o,npzi_o,npzri_o)
!!$    integer, intent(in)  :: mod_pot,ik
!!$    integer, intent(out) :: npzr_o, npzi_o, npzri_o(2)
!!$
!!$    if(mod_pot == NORMCONSERVATION) then
!!$       npzr_o = 0
!!$    else
!!$       if(k_symmetry(ik) == GAMMA) then
!!$          npzr_o = nlmta
!!$       else
!!$          npzr_o = nlmta * 2
!!$       end if
!!$    end if
!!$    npzi_o = npzr_o + iba(ik)
!!$    npzri_o(1) = npzr_o; npzri_o(2) = npzi_o
!!$  end subroutine set_npzri

#ifdef TRANSPOSE
!!$#ifndef TRANSPOSE_ORIGINAL

!!$  subroutine m_ES_W_transpose(k1,k2,ik,psi_l,psi_t)
!!$  end subroutine m_ES_W_transpose

!!$  subroutine m_ES_W_transpose_back(k1,k2,ik,psi_l,psi_t)
!!$  end subroutine m_ES_W_transpose_back

  subroutine m_ES_W_transpose_back2(k1,k2,ik,meg,neordr,psi_l,psi_t)
    integer, intent(in)                                       :: k1,k2,ik,meg
    integer, intent(in), dimension(neg,k1:k2)                 :: neordr
    real(kind=DP), intent(inout),dimension(kg1,np_e,k1:k2,kimg) :: psi_l
    real(kind=DP), intent(in) ,dimension(np_g1k_x,neg,kimg) :: psi_t

    real(kind=DP), allocatable, dimension(:,:,:,:) :: tmp_s,tmp_r
    integer, allocatable, dimension(:) :: req_s,req_r

    integer :: ib, ri, datasize, nb_proc, j, ix,iy
    integer :: pe_s,pe_r
    integer :: mp_g1k_x

    if(npes == 1) then
       do ri = 1, kimg
          do iy = meg+1, neg
             do ix = 1, iba(ik)
                psi_l(ix,iy,ik,ri) = psi_t(ix,neordr(iy,ik),ri)
             end do
          end do
       end do
       return
    end if

#ifdef _ODD_BOUNDARY_
    if(mod(mp_g1k(ik),2) == 0) then
       mp_g1k_x = mp_g1k(ik) + 1
    else
       mp_g1k_x = mp_g1k(ik)
    end if
#else
    mp_g1k_x = mp_g1k(ik)
#endif
    datasize = mp_g1k_x*mp_e*kimg                       ! MPI
    nb_proc  = nrank_e

    allocate(tmp_s(mp_g1k_x,mp_e,kimg,0:nb_proc-1))
    allocate(tmp_r(mp_g1k_x,mp_e,kimg,0:nb_proc-1))
    allocate(req_s(nb_proc-1)) ! MPI
    allocate(req_r(nb_proc-1)) ! MPI

! sbuf + isr
#ifdef USE_NONBLK_COMM
                                                  __TIMER_COMM_START_w_BARRIER(mpi_k_world,582)
#endif
    do j = 1, nb_proc-1
       pe_s=mod(myrank_e+j,        nb_proc)
       pe_r=mod(myrank_e-j+nb_proc,nb_proc)

          do ri = 1, kimg
             do iy = 1, nel_e(pe_s)                   ! MPI
                ib = iy + nis_e(pe_s)-1               ! MPI
                if(ib <= meg) cycle
                do ix = 1, nel_g1k(myrank_e,ik)
                   tmp_s(ix,iy,ri,pe_s) = psi_t(ix,neordr(ib,ik),ri)
                end do
             end do
          end do

#ifdef USE_NONBLK_COMM
       call mpi_irecv(tmp_r(1,1,1,pe_r),datasize,mpi_double_precision, &
                    & pe_r,pe_r    ,mpi_k_world(myrank_k),req_r(j),ierr)
       call mpi_isend(tmp_s(1,1,1,pe_s),datasize,mpi_double_precision, &
                    & pe_s,myrank_e,mpi_k_world(myrank_k),req_s(j),ierr)
#endif

    end do

! diagonal blocks
    do ri = 1, kimg
       do iy = 1, np_e                                              ! MPI
          ib = iy + ista_e-1                                        ! MPI
          if(ib <= meg) cycle                                       ! MPI
          do ix = 1, np_g1k(ik)                                     ! MPI
             psi_l(ista_g1k(ik)-1+ix,iy,ik,ri) = psi_t(ix,neordr(ib,ik),ri)! MPI
          end do
       end do
    end do

! wait
#ifdef USE_NONBLK_COMM
    do j = 1, nb_proc-1
      call mpi_wait(req_r(j),istatus,ierr)
      call mpi_wait(req_s(j),istatus,ierr)
    end do
                                                  __TIMER_COMM_STOP(582)
#else
                                                  __TIMER_COMM_START_w_BARRIER(mpi_k_world,583)
     call MPI_ALLTOALL( tmp_s,datasize,mpi_double_precision, &
    &                   tmp_r,datasize,mpi_double_precision, &
    &                                  mpi_k_world,ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_ES_F_transpose_back2 :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 261, ierr)
     endif
                                                  __TIMER_COMM_STOP(583)
#endif

! rbuf
    do j = 1, nb_proc-1
       pe_s=mod(myrank_e+j,        nb_proc)
       pe_r=mod(myrank_e-j+nb_proc,nb_proc)

          do ri = 1, kimg
             do iy = 1, np_e                               ! MPI
                ib = iy + ista_e - 1                       ! MPI
                if(ib <= meg) cycle
                do ix = 1, nel_g1k(pe_r,ik)
                   psi_l(nis_g1k(pe_r,ik)-1+ix,iy,ik,ri) = tmp_r(ix,iy,ri,pe_r)
                end do
             end do
          end do
    end do
!!$    call mpi_barrier(mpi_k_world(myrank_k),ierr)

    deallocate(tmp_s,tmp_r)
    deallocate(req_s,req_r)
  end subroutine m_ES_W_transpose_back2
!!$#endif !TRANSPOSE_ORIGINAL

!!$#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
!!$#ifndef TRANSPOSE_ORIGINAL
  subroutine m_ES_W_transpose_r(rearrangement,k1,k2,ik,psi_l,psi_t)
    logical, intent(in)                                       :: rearrangement
    integer, intent(in)                                       :: k1,k2,ik
    real(kind=DP), intent(in), dimension(kg1,np_e,k1:k2,kimg) :: psi_l  ! MPI
    real(kind=DP), intent(out),dimension(np_g1k_x,neg,kimg)   :: psi_t  ! MPI

    real(kind=DP), allocatable, dimension(:,:,:,:) :: tmp_s,tmp_r
    integer, allocatable, dimension(:) :: req_s,req_r
    integer, allocatable, dimension(:) :: indexn

    integer :: ri, datasize, nb_proc, j, ix,iy, ifrom, ito
    integer :: pe_s,pe_r
    integer :: mp_g1k_x
    integer :: id_sname = -1

    allocate(indexn(neg))
    if(npes == 1) then
       if(np_g1k(ik) > kg1 .or. neg > np_e) psi_t = 0.d0
       if(rearrangement) then
          indexn(:) = neordr(:,ik)
       else
          indexn(1:neg) = (/(j,j=1,neg)/)
       end if
       if(k_symmetry(ik) == GAMMA) then ! assuming kimg == 2 when k_symmetry(ik) == GAMMA.
          do iy = 1, neg
             ifrom = indexn(iy)
             do ix = 2, iba(ik)
                psi_t(ix,iy,1) = psi_l(ix,ifrom,ik,1)
                psi_t(ix,iy,2) = psi_l(ix,ifrom,ik,2)
             end do
             psi_t(1,iy,1) = psi_l(1,ifrom,ik,1)
             psi_t(1,iy,2) = 0.d0
          end do
       else
          do ri = 1, kimg
             do iy = 1, neg
                ifrom = indexn(iy)
                do ix = 1, iba(ik)
                   psi_t(ix,iy,ri) = psi_l(ix,ifrom,ik,ri)
                end do
             end do
          end do
       end if
       deallocate(indexn)
       return
    end if

    if(sw_timing_2ndlevel == ON) call tstatc0_begin('m_ES_W_transpose_r ',id_sname)

#ifdef _ODD_BOUNDARY_
    if(mod(mp_g1k(ik),2) == 0) then
       mp_g1k_x = mp_g1k(ik) + 1
    else
       mp_g1k_x = mp_g1k(ik)
    end if
#else
    mp_g1k_x = mp_g1k(ik)
#endif
    datasize = mp_g1k_x*mp_e*kimg                       ! MPI
    nb_proc  = nrank_e                               ! MPI

    allocate(tmp_s(mp_g1k_x,mp_e,kimg,0:nb_proc-1)); tmp_s = 0.d0 ! MPI
    allocate(tmp_r(mp_g1k_x,mp_e,kimg,0:nb_proc-1)); tmp_r = 0.d0 ! MPI
    allocate(req_s(nb_proc-1)) ! MPI
    allocate(req_r(nb_proc-1)) ! MPI

    psi_t = 0.d0

    if(rearrangement) then
       indexn(1:neg) = nrvf_ordr(1:neg,ik)
    else
       indexn(1:neg) = (/(j,j=1,neg)/)
    end if

! sbuf + isr
#ifdef USE_NONBLK_COMM
                                                  __TIMER_COMM_START_w_BARRIER(mpi_k_world,582)
#endif
    do j = 1, nb_proc-1
       pe_s=mod(myrank_e+j,        nb_proc)
       pe_r=mod(myrank_e-j+nb_proc,nb_proc)

          do ri = 1, kimg
             do iy = 1, np_e                            ! MPI
                do ix = 1, nel_g1k(pe_s,ik)       ! MPI
                   tmp_s(ix,iy,ri,pe_s) = psi_l(nis_g1k(pe_s,ik)-1+ix,iy,ik,ri) ! MPI
                end do
             end do
          end do

#ifdef USE_NONBLK_COMM
       call mpi_irecv(tmp_r(1,1,1,pe_r),datasize,mpi_double_precision, &
                    & pe_r,pe_r    ,mpi_k_world(myrank_k),req_r(j),ierr)
       call mpi_isend(tmp_s(1,1,1,pe_s),datasize,mpi_double_precision, &
                    & pe_s,myrank_e,mpi_k_world(myrank_k),req_s(j),ierr)
#endif

    end do

! diagonal blocks
    do ri = 1, kimg
       do iy = 1, np_e                               ! MPI
          ito = indexn(ista_e-1+iy)
          do ix = 1, np_g1k(ik)                      ! MPI
             psi_t(ix,ito,ri) = psi_l(ista_g1k(ik)-1+ix,iy,ik,ri)  ! MPI
          end do
       end do
    end do

! wait
#ifdef USE_NONBLK_COMM
    do j = 1, nb_proc-1
      call mpi_wait(req_r(j),istatus,ierr)
      call mpi_wait(req_s(j),istatus,ierr)
    end do
                                                  __TIMER_COMM_STOP(582)
#else
                                                  __TIMER_COMM_START_w_BARRIER(mpi_k_world,583)
     call MPI_ALLTOALL( tmp_s,datasize,mpi_double_precision, &
    &                   tmp_r,datasize,mpi_double_precision, &
    &                                  mpi_k_world,ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_ES_W_transpose_r :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 261, ierr)
     endif
                                                  __TIMER_COMM_STOP(583)
#endif

! rbuf
    do j = 1, nb_proc-1
!!$       pe_s=mod(myrank_e+j,        nb_proc)
       pe_r=mod(myrank_e-j+nb_proc,nb_proc)
          do ri = 1, kimg
             do iy = 1, nel_e(pe_r)
                ito = indexn(nis_e(pe_r)-1+iy)
                do ix = 1, nel_g1k(myrank_e,ik)
                   psi_t(ix,ito,ri) = tmp_r(ix,iy,ri,pe_r)
                end do
             end do
          end do
    end do

    if(k_symmetry(ik) == GAMMA .and. myrank_e == 0 .and. kimg == 2) then
       do iy = 1, neg
          psi_t(1,iy,2) = 0.d0
       end do
    end if

!!$    call mpi_barrier(mpi_k_world(myrank_k),ierr)
    deallocate(tmp_s,tmp_r)
    deallocate(req_s,req_r)

    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
  end subroutine m_ES_W_transpose_r

  subroutine m_ES_W_transpose_back_r(rearrangement,k1,k2,ik,psi_l,psi_t)
    logical, intent(in)  :: rearrangement
    integer, intent(in)                                       :: k1,k2,ik
    real(kind=DP), intent(inout),dimension(kg1,np_e,k1:k2,kimg) :: psi_l
    real(kind=DP), intent(in) ,dimension(np_g1k_x,neg,kimg)   :: psi_t

    real(kind=DP), allocatable, dimension(:,:,:,:) :: tmp_s,tmp_r
    integer, allocatable, dimension(:) :: req_s,req_r,indexn

    integer :: ri, datasize, nb_proc, j, ix,iy, ifrom,ito
    integer :: pe_s,pe_r
    integer :: mp_g1k_x
    integer :: id_sname = -1

    allocate(indexn(neg))
    if(npes == 1) then
       if(rearrangement) then
          indexn(1:neg) = neordr(1:neg,ik)
       else
          indexn(1:neg) = (/(j,j=1,neg)/)
       end if

       do ri = 1, kimg
          do iy = 1, neg
!!$             ito = neordr(iy,ik)
             ito = indexn(iy)
             do ix = 1, iba(ik)
                psi_l(ix,ito,ik,ri) = psi_t(ix,iy,ri)
             end do
          end do
       end do
       deallocate(indexn)
       return
    end if

    if(sw_timing_2ndlevel == ON) call tstatc0_begin('m_ES_W_transpose_back_r ',id_sname)

#ifdef _ODD_BOUNDARY_
    if(mod(mp_g1k(ik),2) == 0) then
       mp_g1k_x = mp_g1k(ik) + 1
    else
       mp_g1k_x = mp_g1k(ik)
    end if
#else
    mp_g1k_x = mp_g1k(ik)
#endif
    datasize = mp_g1k_x*mp_e*kimg                       ! MPI
    nb_proc  = nrank_e

    allocate(tmp_s(mp_g1k_x,mp_e,kimg,0:nb_proc-1)); tmp_s = 0.d0 ! MPI
    allocate(tmp_r(mp_g1k_x,mp_e,kimg,0:nb_proc-1)); tmp_r = 0.d0 ! MPI
    allocate(req_s(nb_proc-1)) ! MPI
    allocate(req_r(nb_proc-1)) ! MPI

    if(rearrangement) then
       indexn(1:neg) = nrvf_ordr(1:neg,ik)
    else
       indexn(1:neg) = (/(j,j=1,neg)/)
    end if

! sbuf + isr
#ifdef USE_NONBLK_COMM
                                                  __TIMER_COMM_START_w_BARRIER(mpi_k_world,582)
#endif
    do j = 1, nb_proc-1
       pe_s=mod(myrank_e+j,        nb_proc)
       pe_r=mod(myrank_e-j+nb_proc,nb_proc)

          do ri = 1, kimg
             do iy = 1, nel_e(pe_s)                   ! MPI
!!$                ifrom = nrvf_ordr(nis_e(pe_s)-1+iy,ik)
                ifrom = indexn(nis_e(pe_s)-1+iy)
                do ix = 1, nel_g1k(myrank_e,ik)             ! MPI
                   tmp_s(ix,iy,ri,pe_s) = psi_t(ix,ifrom,ri) ! MPI
                end do
             end do
          end do

#ifdef USE_NONBLK_COMM
       call mpi_irecv(tmp_r(1,1,1,pe_r),datasize,mpi_double_precision, &
                    & pe_r,pe_r    ,mpi_k_world(myrank_k),req_r(j),ierr)
       call mpi_isend(tmp_s(1,1,1,pe_s),datasize,mpi_double_precision, &
                    & pe_s,myrank_e,mpi_k_world(myrank_k),req_s(j),ierr)
#endif

    end do

! diagonal blocks
    do ri = 1, kimg
       do iy = 1, np_e                                               ! MPI
!!$          ifrom = nrvf_ordr(ista_e-1+iy,ik)
          ifrom = indexn(ista_e-1+iy)
          do ix = 1, np_g1k(ik)                                      ! MPI
             psi_l(ista_g1k(ik)-1+ix,iy,ik,ri) = psi_t(ix,ifrom,ri) ! MPI
          end do
       end do
    end do

! wait
#ifdef USE_NONBLK_COMM
    do j = 1, nb_proc-1
      call mpi_wait(req_r(j),istatus,ierr)
      call mpi_wait(req_s(j),istatus,ierr)
    end do
                                                  __TIMER_COMM_STOP(582)
#else
                                                  __TIMER_COMM_START_w_BARRIER(mpi_k_world,583)
     call MPI_ALLTOALL( tmp_s,datasize,mpi_double_precision, &
    &                   tmp_r,datasize,mpi_double_precision, &
    &                                  mpi_k_world,ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_ES_W_transpose_back_r :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 261, ierr)
     endif
                                                  __TIMER_COMM_STOP(583)
#endif

! rbuf
    do j = 1, nb_proc-1
!!$       pe_s=mod(myrank_e+j,        nb_proc)
       pe_r=mod(myrank_e-j+nb_proc,nb_proc)

          do ri = 1, kimg
             do iy = 1, np_e                                ! MPI
                do ix = 1, nel_g1k(pe_r,ik)
                   psi_l(nis_g1k(pe_r,ik)-1+ix,iy,ik,ri) = tmp_r(ix,iy,ri,pe_r)
                end do
             end do
          end do

    end do

!!$    call mpi_barrier(mpi_k_world(myrank_k),ierr)

    deallocate(tmp_s,tmp_r)
    deallocate(req_s,req_r)

    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)

    deallocate(indexn)

  end subroutine m_ES_W_transpose_back_r

  subroutine m_ES_W_transpose_back_r_cmplx(rearrangement,k1,k2,ik,psi_l,psi_t)
    logical, intent(in)  :: rearrangement
    integer, intent(in)                                       :: k1,k2,ik
    real(kind=DP), intent(inout),dimension(kg1,np_e,k1:k2,kimg) :: psi_l
    real(kind=DP), intent(in) ,dimension(kimg*np_g1k_x,neg,1)   :: psi_t

    real(kind=DP), allocatable, dimension(:,:,:,:) :: tmp_s,tmp_r
    integer, allocatable, dimension(:) :: req_s,req_r,indexn

    integer :: ri, datasize, nb_proc, j, ix,iy, ifrom,ito
    integer :: pe_s,pe_r
    integer :: mp_g1k_x
    integer :: id_sname = -1

    allocate(indexn(neg))
    if(npes == 1) then
       if(rearrangement) then
          indexn(1:neg) = neordr(1:neg,ik)
       else
          indexn(1:neg) = (/(j,j=1,neg)/)
       end if

       do iy = 1, neg
!!$          ito = neordr(iy,ik)
          ito = indexn(iy)
          do ix = 1, iba(ik)
             psi_l(ix,ito,ik,1) = psi_t(kimg*ix-1,iy,1)
             if(kimg==2) psi_l(ix,ito,ik,2) = psi_t(kimg*ix,iy,1)
          end do
       end do
       deallocate(indexn)
       return
    end if

    if(sw_timing_2ndlevel == ON) call tstatc0_begin('m_ES_W_transpose_back_r ',id_sname)

#ifdef _ODD_BOUNDARY_
    if(mod(mp_g1k(ik),2) == 0) then
       mp_g1k_x = mp_g1k(ik) + 1
    else
       mp_g1k_x = mp_g1k(ik)
    end if
#else
    mp_g1k_x = mp_g1k(ik)
#endif
    datasize = mp_g1k_x*mp_e*kimg                       ! MPI
    nb_proc  = nrank_e

    allocate(tmp_s(mp_g1k_x,mp_e,kimg,0:nb_proc-1)); tmp_s = 0.d0 ! MPI
    allocate(tmp_r(mp_g1k_x,mp_e,kimg,0:nb_proc-1)); tmp_r = 0.d0 ! MPI
    allocate(req_s(nb_proc-1)) ! MPI
    allocate(req_r(nb_proc-1)) ! MPI

    if(rearrangement) then
       indexn(1:neg) = nrvf_ordr(1:neg,ik)
    else
       indexn(1:neg) = (/(j,j=1,neg)/)
    end if

! sbuf + isr
#ifdef USE_NONBLK_COMM
                                                  __TIMER_COMM_START_w_BARRIER(mpi_k_world,582)
#endif
    do j = 1, nb_proc-1
       pe_s=mod(myrank_e+j,        nb_proc)
       pe_r=mod(myrank_e-j+nb_proc,nb_proc)

!          do ri = 1, kimg
             do iy = 1, nel_e(pe_s)                   ! MPI
!!$                ifrom = nrvf_ordr(nis_e(pe_s)-1+iy,ik)
                ifrom = indexn(nis_e(pe_s)-1+iy)
                do ix = 1, nel_g1k(myrank_e,ik)             ! MPI
                   tmp_s(ix,iy,1,pe_s) = psi_t(kimg*ix-1,ifrom,1) ! MPI
                   if(kimg==2) tmp_s(ix,iy,2,pe_s) = psi_t(kimg*ix,ifrom,1) ! MPI
                end do
             end do
!          end do

#ifdef USE_NONBLK_COMM
       call mpi_irecv(tmp_r(1,1,1,pe_r),datasize,mpi_double_precision, &
                    & pe_r,pe_r    ,mpi_k_world(myrank_k),req_r(j),ierr)
       call mpi_isend(tmp_s(1,1,1,pe_s),datasize,mpi_double_precision, &
                    & pe_s,myrank_e,mpi_k_world(myrank_k),req_s(j),ierr)
#endif

    end do

! diagonal blocks
!    do ri = 1, kimg
     do iy = 1, np_e                                               ! MPI
!!$        ifrom = nrvf_ordr(ista_e-1+iy,ik)
        ifrom = indexn(ista_e-1+iy)
        do ix = 1, np_g1k(ik)                                      ! MPI
           psi_l(ista_g1k(ik)-1+ix,iy,ik,1) = psi_t(kimg*ix-1,ifrom,1) ! MPI
           if(kimg==2) psi_l(ista_g1k(ik)-1+ix,iy,ik,2) = psi_t(kimg*ix,ifrom,1) ! MPI
        end do
     end do
!    end do

! wait
#ifdef USE_NONBLK_COMM
    do j = 1, nb_proc-1
      call mpi_wait(req_r(j),istatus,ierr)
      call mpi_wait(req_s(j),istatus,ierr)
    end do
                                                  __TIMER_COMM_STOP(582)
#else
                                                  __TIMER_COMM_START_w_BARRIER(mpi_k_world,583)
     call MPI_ALLTOALL( tmp_s,datasize,mpi_double_precision, &
    &                   tmp_r,datasize,mpi_double_precision, &
    &                                  mpi_k_world,ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_ES_W_transpose_back_r :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 261, ierr)
     endif
                                                  __TIMER_COMM_STOP(583)
#endif

! rbuf
    do j = 1, nb_proc-1
!!$       pe_s=mod(myrank_e+j,        nb_proc)
       pe_r=mod(myrank_e-j+nb_proc,nb_proc)

          do ri = 1, kimg
             do iy = 1, np_e                                ! MPI
                do ix = 1, nel_g1k(pe_r,ik)
                   psi_l(nis_g1k(pe_r,ik)-1+ix,iy,ik,ri) = tmp_r(ix,iy,ri,pe_r)
                end do
             end do
          end do

    end do

!!$    call mpi_barrier(mpi_k_world(myrank_k),ierr)

    deallocate(tmp_s,tmp_r)
    deallocate(req_s,req_r)

    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)

    deallocate(indexn)

  end subroutine m_ES_W_transpose_back_r_cmplx


  subroutine m_ES_W_transpose_fft_cmplx(nffth,np_ffth,psi_l,psi_t)
    integer, intent(in) :: nffth,np_ffth
    complex(kind=DP), intent(in), dimension(nffth,np_e) :: psi_l  ! MPI
    complex(kind=DP), intent(out),dimension(np_ffth,neg)   :: psi_t  ! MPI

    complex(kind=DP), allocatable, dimension(:,:,:) :: tmp_s,tmp_r
    integer, allocatable, dimension(:) :: req_s,req_r
    integer, allocatable, dimension(:) :: indexn
    complex(kind=DP) :: zero=(0.d0,0.d0)
    integer :: datasize, nb_proc, j, ix, iy, ito
    integer :: pe_s,pe_r
    integer :: id_sname = -1

    if(npes == 1) then
       psi_t = zero
       do iy = 1, neg
          do ix = 1, nffth
             psi_t(ix,iy) = psi_l(ix,iy)
          end do
       end do
       return
    end if

    if(sw_timing_2ndlevel == ON) call tstatc0_begin('m_ES_W_transpose_fft_cmplx ',id_sname)

    datasize = mp_ffth*mp_e
    nb_proc  = nrank_e                               ! MPI

    allocate(tmp_s(mp_ffth,mp_e,0:nb_proc-1)); tmp_s = zero ! MPI
    allocate(tmp_r(mp_ffth,mp_e,0:nb_proc-1)); tmp_r = zero ! MPI
    allocate(req_s(nb_proc-1)) ! MPI
    allocate(req_r(nb_proc-1)) ! MPI

    psi_t = zero

! sbuf + isr
    do j = 1, nb_proc-1
       pe_s=mod(myrank_e+j,        nb_proc)
       pe_r=mod(myrank_e-j+nb_proc,nb_proc)

       do iy = 1, np_e                            ! MPI
          do ix = 1, nel_ffth(pe_s)       ! MPI
             tmp_s(ix,iy,pe_s) = psi_l(is_ffth(pe_s)-1+ix,iy) ! MPI
          end do
       end do

#ifdef USE_NONBLK_COMM
       call mpi_irecv(tmp_r(1,1,pe_r),datasize,mpi_double_complex, &
                    & pe_r,pe_r    ,mpi_k_world(myrank_k),req_r(j),ierr)
       call mpi_isend(tmp_s(1,1,pe_s),datasize,mpi_double_complex, &
                    & pe_s,myrank_e,mpi_k_world(myrank_k),req_s(j),ierr)
#endif

    end do

! diagonal blocks
    do iy = 1, np_e                               ! MPI
       ito = ista_e-1+iy
       do ix = 1, np_ffth                      ! MPI
          psi_t(ix,ito) = psi_l(ista_ffth-1+ix,iy)  ! MPI
       end do
    end do

! wait
#ifdef USE_NONBLK_COMM
    do j = 1, nb_proc-1
      call mpi_wait(req_r(j),istatus,ierr)
      call mpi_wait(req_s(j),istatus,ierr)
    end do
#else
    call MPI_ALLTOALL( tmp_s,datasize,mpi_double_complex, &
   &                   tmp_r,datasize,mpi_double_complex, &
   &                                  mpi_k_world,ierr)
    if (ierr /= 0) then
       write(nfout,*)' m_ES_W_transpose_fft_cmplx :  mpi_alltoall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 261, ierr)
    endif
#endif

! rbuf
    do j = 1, nb_proc-1
!!$       pe_s=mod(myrank_e+j,        nb_proc)
       pe_r=mod(myrank_e-j+nb_proc,nb_proc)

       do iy = 1, nel_e(pe_r)
          ito = nis_e(pe_r)-1+iy
          do ix = 1, nel_ffth(myrank_e)
             psi_t(ix,ito) = tmp_r(ix,iy,pe_r)
          end do
       end do

    end do

!!$    call mpi_barrier(mpi_k_world(myrank_k),ierr)
    deallocate(tmp_s,tmp_r)
    deallocate(req_s,req_r)

    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
  end subroutine m_ES_W_transpose_fft_cmplx

  subroutine m_ES_W_transpose_back_fft_cmplx(nffth,nffth_el,psi_l,psi_t)
    integer, intent(in) :: nffth,nffth_el
    complex(kind=DP), intent(out), dimension(nffth,np_e) :: psi_l  ! MPI
    complex(kind=DP), intent(in),dimension(nffth_el,neg)   :: psi_t  ! MPI

    complex(kind=DP), allocatable, dimension(:,:,:) :: tmp_s,tmp_r
    integer, allocatable, dimension(:) :: req_s,req_r

    integer :: datasize, nb_proc, j, ix, iy, ito, ifrom
    integer :: pe_s,pe_r
    complex(kind=DP) :: zero=(0.d0,0.d0)

    integer :: id_sname = -1

    if(npes == 1) then
       do iy = 1, neg
          ito = iy
          do ix = 1, nffth
             psi_l(ix,ito) = psi_t(ix,iy)
          end do
       end do
       return
    end if

    if(sw_timing_2ndlevel == ON) call tstatc0_begin('m_ES_W_transpose_back_fft_cmplx ',id_sname)

    datasize = mp_ffth*mp_e
    nb_proc  = nrank_e

    allocate(tmp_s(mp_ffth,mp_e,0:nb_proc-1)); tmp_s = zero ! MPI
    allocate(tmp_r(mp_ffth,mp_e,0:nb_proc-1)); tmp_r = zero ! MPI
    allocate(req_s(nb_proc-1)) ! MPI
    allocate(req_r(nb_proc-1)) ! MPI

! sbuf + isr
    do j = 1, nb_proc-1
       pe_s=mod(myrank_e+j,        nb_proc)
       pe_r=mod(myrank_e-j+nb_proc,nb_proc)

       do iy = 1, nel_e(pe_s)                   ! MPI
          ifrom = nis_e(pe_s)-1+iy
          do ix = 1, nel_ffth(myrank_e)             ! MPI
             tmp_s(ix,iy,pe_s) = psi_t(ix,ifrom) ! MPI
          end do
       end do

#ifdef USE_NONBLK_COMM
       call mpi_irecv(tmp_r(1,1,pe_r),datasize,mpi_double_complex, &
                    & pe_r,pe_r    ,mpi_k_world(myrank_k),req_r(j),ierr)
       call mpi_isend(tmp_s(1,1,pe_s),datasize,mpi_double_complex, &
                    & pe_s,myrank_e,mpi_k_world(myrank_k),req_s(j),ierr)
#endif

    end do

! diagonal blocks
    do iy = 1, np_e                                               ! MPI
       ifrom = ista_e-1+iy
       do ix = 1, nffth_el
          psi_l(is_ffth-1+ix,iy) = psi_t(ix,ifrom) ! MPI
       end do
    end do

! wait
#ifdef USE_NONBLK_COMM
    do j = 1, nb_proc-1
      call mpi_wait(req_r(j),istatus,ierr)
      call mpi_wait(req_s(j),istatus,ierr)
    end do
#else
    call MPI_ALLTOALL( tmp_s,datasize,mpi_double_complex, &
   &                   tmp_r,datasize,mpi_double_complex, &
   &                                  mpi_k_world,ierr)
    if (ierr /= 0) then
       write(nfout,*)' m_ES_W_transpose_back_fft_cmplx :  mpi_alltoall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 261, ierr)
    endif
#endif

! rbuf
    do j = 1, nb_proc-1
       pe_r=mod(myrank_e-j+nb_proc,nb_proc)
       do iy = 1, np_e                                ! MPI
          do ix = 1, nel_ffth(pe_r)
             psi_l(is_ffth(pe_r)-1+ix,iy) = tmp_r(ix,iy,pe_r)
          end do
       end do
    end do

    deallocate(tmp_s,tmp_r)
    deallocate(req_s,req_r)

    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)

  end subroutine m_ES_W_transpose_back_fft_cmplx

#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
!!$#ifndef TRANSPOSE_ORIGINAL
  subroutine m_ES_F_transpose_r(rearrangement,k1,k2,ik,bpr_l,bpr_t,bpi_l,bpi_t)
    logical, intent(in) :: rearrangement
    integer, intent(in)                                                :: k1,k2,ik
    real(kind=DP), intent(in), dimension(np_e,nlmta,k1:k2)             :: bpr_l
    real(kind=DP), intent(out),dimension(np_fs_x,neg)                  :: bpr_t
    real(kind=DP), intent(in), optional,dimension(np_e,nlmta,k1:k2)    :: bpi_l
    real(kind=DP), intent(out),optional,dimension(np_fs_x,neg)         :: bpi_t

    real(kind=DP), allocatable, dimension(:,:,:,:) :: tmp_s,tmp_r
    integer, allocatable, dimension(:) :: req_s,req_r,indexn

    integer :: datasize, nb_proc, j, ix,iy, kimg_t, ifrom, ito
    integer :: pe_s,pe_r
    integer :: id_sname = -1
#ifdef MGS_DGEMM_DEBUG
    integer :: ia,p,q
#endif

    if(k_symmetry(ik) == GAMMA) then ! assuming kimg == 2 when k_symmetry(ik) == GAMMA.
       kimg_t = 1
    else
       kimg_t = 2
    end if

    if(npes == 1) then
       allocate(indexn(neg))
       if(rearrangement) then
          indexn(1:neg) = neordr(1:neg,ik)
       else
          indexn(1:neg) = (/(j,j=1,neg)/)
       end if
       if(kimg_t == 1) then
          if(np_fs_x > nlmta .or. neg > np_e) bpr_t = 0.d0
          do iy = 1, neg
!!$             ifrom = neordr(iy,ik)
             ifrom = indexn(iy)
             do ix = 1, nlmta
                bpr_t(ix,iy) = bpr_l(ifrom,ix,ik)
             end do
          end do
       else if(kimg_t == 2) then
          if(np_fs_x > nlmta .or. neg > np_e)then
             bpr_t = 0.d0; bpi_t = 0.d0
          end if

          do iy = 1, neg
!!$             ifrom = neordr(iy,ik)
             ifrom = indexn(iy)
             do ix = 1, nlmta
                bpr_t(ix,iy) = bpr_l(ifrom,ix,ik)
                bpi_t(ix,iy) = bpi_l(ifrom,ix,ik)
             end do
          end do
       end if

       deallocate(indexn)
#ifdef MGS_DGEMM_DEBUG
       goto 9999
#else
       return
#endif
    end if

    if(sw_timing_2ndlevel == ON) call tstatc0_begin('m_ES_F_transpose_r ',id_sname)

    datasize = mp_fs*mp_e*kimg_t
    nb_proc  = nrank_e

    if(ipri >= 3 )then
       write(nfout,'(" !m_ES_F_transpose  datasize = ",i5)') datasize
       write(nfout,'(" !m_ES_F_transpose  nb_proc  = ",i5)') nb_proc
       write(nfout,'(" !m_ES_F_transpose  kimg_t   = ",i5)') kimg_t
    end if

    allocate(tmp_s(mp_fs,mp_e,kimg_t,0:nb_proc-1)); tmp_s = 0.d0
    allocate(tmp_r(mp_fs,mp_e,kimg_t,0:nb_proc-1)); tmp_r = 0.d0
    allocate(req_s(nb_proc-1))
    allocate(req_r(nb_proc-1))

    bpr_t = 0.d0
    if(kimg_t == 2) bpi_t = 0.d0

    allocate(indexn(neg))
    if(rearrangement) then
       indexn(1:neg) = nrvf_ordr(1:neg,ik)
    else
       indexn(1:neg) = (/(j,j=1,neg)/)
    end if
! sbuf + istr
    do j = 1, nb_proc-1
       pe_s=mod(myrank_e+j,        nb_proc)
       pe_r=mod(myrank_e-j+nb_proc,nb_proc)

          if(kimg_t == 1) then
             do iy = 1, np_e
                do ix = 1, nel_fs(pe_s)
                   tmp_s(ix,iy,1,pe_s) = bpr_l(iy,nis_fs(pe_s)-1+ix,ik)
                end do
             end do
          else
             do iy = 1, np_e
                do ix = 1, nel_fs(pe_s)
                   tmp_s(ix,iy,1,pe_s) = bpr_l(iy,nis_fs(pe_s)-1+ix,ik)
                   tmp_s(ix,iy,2,pe_s) = bpi_l(iy,nis_fs(pe_s)-1+ix,ik)
                end do
             end do
          end if

#ifdef USE_NONBLK_COMM
       call mpi_irecv(tmp_r(1,1,1,pe_r),datasize,mpi_double_precision, &
                    & pe_r,pe_r    ,mpi_k_world(myrank_k),req_r(j),ierr)
       call mpi_isend(tmp_s(1,1,1,pe_s),datasize,mpi_double_precision, &
                    & pe_s,myrank_e,mpi_k_world(myrank_k),req_s(j),ierr)
#endif

    end do

    if(ipri >= 3 )then
       write(nfout,'(" !m_ES_F_transpose sbuf + istr = OK")')
    end if

! diagonal blocks
    if(kimg_t == 1) then
       do iy = 1, np_e
!!$          ito = nrvf_ordr(ista_e-1+iy,ik)
          ito = indexn(ista_e-1+iy)
          do ix = 1, np_fs
             bpr_t(ix,ito) = bpr_l(iy,ista_fs-1+ix,ik)
          end do
       end do
    else
       do iy = 1, np_e
!!$          ito = nrvf_ordr(ista_e-1+iy,ik)
          ito = indexn(ista_e-1+iy)
          do ix = 1, np_fs
             bpr_t(ix,ito) = bpr_l(iy,ista_fs-1+ix,ik)
             bpi_t(ix,ito) = bpi_l(iy,ista_fs-1+ix,ik)
          end do
       end do
    end if

    if(ipri >= 3 )then
       write(nfout,'(" !m_ES_F_transpose  diagonal block = OK")')
    end if

! wait
#ifdef USE_NONBLK_COMM
    do j = 1, nb_proc-1
      call mpi_wait(req_r(j),istatus,ierr)
      call mpi_wait(req_s(j),istatus,ierr)
    end do
                                                  __TIMER_COMM_STOP(582)
#else
                                                  __TIMER_COMM_START_w_BARRIER(mpi_k_world,583)
     call MPI_ALLTOALL( tmp_s,datasize,mpi_double_precision, &
    &                   tmp_r,datasize,mpi_double_precision, &
    &                                  mpi_k_world,ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_ES_F_transpose :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 261, ierr)
     endif
                                                  __TIMER_COMM_STOP(583)
#endif

! rbuf
    do j = 1, nb_proc-1
       pe_r=mod(myrank_e-j+nb_proc,nb_proc)

          if(kimg_t == 1) then
             do iy = 1, nel_e(pe_r)
!!$                ito = nrvf_ordr(nis_e(pe_r)-1+iy,ik)
                ito = indexn(nis_e(pe_r)-1+iy)
                do ix = 1, nel_fs(myrank_e)
                   bpr_t(ix,ito) = tmp_r(ix,iy,1,pe_r)
                end do
             end do
          else if(kimg_t == 2) then
             do iy = 1, nel_e(pe_r)
!!$                ito = nrvf_ordr(nis_e(pe_r)-1+iy,ik)
                ito = indexn(nis_e(pe_r)-1+iy)
                do ix = 1, nel_fs(myrank_e)
                   bpr_t(ix,ito) = tmp_r(ix,iy,1,pe_r)
                   bpi_t(ix,ito) = tmp_r(ix,iy,2,pe_r)
                end do
             end do
          end if
    end do

    if(ipri >= 3 )then
       write(nfout,'(" !m_ES_F_transpose  rbuf part = OK")')
    end if

!!$    call mpi_barrier(mpi_k_world(myrank_k),ierr)
    deallocate(tmp_s,tmp_r)
    deallocate(req_s,req_r)
    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
    deallocate(indexn)

#ifdef MGS_DGEMM_DEBUG
 9999 continue
    if(kimg_t == 1) then
       do j = 1, neg
         do ia = 1, nac_p
            p = nlmta1_p(ia);     q = nlmta2_p(ia)
            bpr_tw1(ia,j) = bpr_t(p,j)
            bpr_tw2(ia,j) = bpr_t(q,j)
         end do
       end do
    else if(kimg_t == 2) then
       do j = 1, neg
         do ia = 1, nac_p
            p = nlmta1_p(ia);     q = nlmta2_p(ia)
            bpr_tw1(ia,j) = bpr_t(p,j)
            bpr_tw2(ia,j) = bpr_t(q,j)
            bpi_tw1(ia,j) = bpi_t(p,j)
            bpi_tw2(ia,j) = bpi_t(q,j)
         end do
       end do
    end if
#ifdef DEBUG_MGS
    if(ipri >= 1 )then
       if(kimg_t == 1) then
          write(nfout,'(" bpr_tw1, bpr_tw2")')
          do j = 1, 5
             write(nfout,'(" j = ",i8)') j
             do ia = 1, nac_p
                write(nfout,'(" ia = ",i8," : ",2f8.4)') ia,bpr_tw1(ia,j),bpr_tw2(ia,j)
             end do
          end do
       else if(kimg_t == 2) then
          write(nfout,'(" bpr_tw1, bpr_tw2, bpi_tw1, bpi_tw2")')
          do j = 1, 5
             write(nfout,'(" j = ",i8)') j
             do ia = 1, nac_p
                write(nfout,'(" ia = ",i8," : ",4f8.4)') ia, bpr_tw1(ia,j),bpr_tw2(ia,j),bpi_tw1(ia,j),bpi_tw2(ia,j)
             end do
          end do
       end if
    end if
#endif
#endif

  end subroutine m_ES_F_transpose_r

  subroutine m_ES_F_transpose_back_r(rearrangement,k1,k2,ik,bpr_l,bpr_t,bpi_l,bpi_t)
    logical, intent(in) :: rearrangement
    integer, intent(in)                                                :: k1,k2,ik
    real(kind=DP), intent(inout),dimension(np_e,nlmta,k1:k2)             :: bpr_l
    real(kind=DP), intent(in) ,dimension(np_fs_x,neg)                  :: bpr_t
    real(kind=DP), intent(inout),optional,dimension(np_e,nlmta,k1:k2)    :: bpi_l
    real(kind=DP), intent(in) ,optional,dimension(np_fs_x,neg)         :: bpi_t

    real(kind=DP), allocatable, dimension(:,:,:,:) :: tmp_s,tmp_r
    integer, allocatable, dimension(:) :: req_s,req_r,indexn

    integer :: datasize, nb_proc, j, ix,iy, kimg_t,ito,ifrom
    integer :: pe_s,pe_r
    integer :: id_sname = -1

    if(k_symmetry(ik) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = 2
    end if

    if(npes == 1) then
       allocate(indexn(1:neg))
       if(rearrangement) then
          indexn(1:neg) = neordr(1:neg,ik)
       else
          indexn(1:neg) = (/(j,j=1,neg)/)
       end if
       if(kimg_t == 1) then
          do iy = 1, neg
             ito = neordr(iy,ik)
             do ix = 1, nlmta
                bpr_l(ito,ix,ik) = bpr_t(ix,iy)
             end do
          end do
       else
          do iy = 1, neg
             ito = neordr(iy,ik)
             do ix = 1, nlmta
                bpr_l(ito,ix,ik) = bpr_t(ix,iy)
                bpi_l(ito,ix,ik) = bpi_t(ix,iy)
             end do
          end do
       end if
       deallocate(indexn)
       return
    end if

    if(sw_timing_2ndlevel == ON) call tstatc0_begin('m_ES_F_transpose_back_r ',id_sname)

    datasize = mp_fs*mp_e*kimg_t
    nb_proc  = nrank_e

    if(ipri >= 3 )then
       write(nfout,'(" !m_ES_F_transpose_back  datasize = ",i5)') datasize
       write(nfout,'(" !m_ES_F_transpose_back  nb_proc  = ",i5)') nb_proc
       write(nfout,'(" !m_ES_F_transpose_back  kimg_t   = ",i5)') kimg_t
    end if

    allocate(tmp_s(mp_fs,mp_e,kimg_t,0:nb_proc-1)); tmp_s = 0.d0
    allocate(tmp_r(mp_fs,mp_e,kimg_t,0:nb_proc-1)); tmp_r = 0.d0
    allocate(req_s(nb_proc-1))
    allocate(req_r(nb_proc-1))

    allocate(indexn(1:neg))
    if(rearrangement) then
       indexn(1:neg) = nrvf_ordr(1:neg,ik)
    else
       indexn(1:neg) = (/(j,j=1,neg)/)
    end if
! sbuf + isr
#ifdef USE_NONBLK_COMM
                                                  __TIMER_COMM_START_w_BARRIER(mpi_k_world,582)
#endif
    do j = 1, nb_proc-1
       pe_s=mod(myrank_e+j,        nb_proc)
       pe_r=mod(myrank_e-j+nb_proc,nb_proc)

          if(kimg_t == 1) then
             do iy = 1, nel_e(pe_s)
!!$                ifrom = nrvf_ordr(nis_e(pe_s)-1+iy,ik)
                ifrom = indexn(nis_e(pe_s)-1+iy)
                do ix = 1, nel_fs(myrank_e)
                   tmp_s(ix,iy,1,pe_s) = bpr_t(ix,ifrom)
                end do
             end do
          else
             do iy = 1, nel_e(pe_s)
!!$                ifrom = nrvf_ordr(nis_e(pe_s)-1+iy,ik)
                ifrom = indexn(nis_e(pe_s)-1+iy)
                do ix = 1, nel_fs(myrank_e)
                   tmp_s(ix,iy,1,pe_s) = bpr_t(ix,ifrom)
                   tmp_s(ix,iy,2,pe_s) = bpi_t(ix,ifrom)
                end do
             end do
          end if

#ifdef USE_NONBLK_COMM
       call mpi_irecv(tmp_r(1,1,1,pe_r),datasize,mpi_double_precision, &
                    & pe_r,pe_r    ,mpi_k_world(myrank_k),req_r(j),ierr)
       call mpi_isend(tmp_s(1,1,1,pe_s),datasize,mpi_double_precision, &
                    & pe_s,myrank_e,mpi_k_world(myrank_k),req_s(j),ierr)
#endif

    end do

! diagonal blocks
    if(kimg_t == 1) then
       do ix = 1, np_fs
          do iy = 1, np_e
!!$             ifrom = nrvf_ordr(ista_e-1+iy,ik)
             ifrom = indexn(ista_e-1+iy)
             bpr_l(iy,ista_fs-1+ix,ik) = bpr_t(ix,ifrom)
          end do
       end do
    else
       do ix = 1, np_fs
          do iy = 1, np_e
!!$             ifrom = nrvf_ordr(ista_e-1+iy,ik)
             ifrom = indexn(ista_e-1+iy)
             bpr_l(iy,ista_fs-1+ix,ik) = bpr_t(ix,ifrom)
             bpi_l(iy,ista_fs-1+ix,ik) = bpi_t(ix,ifrom)
          end do
       end do
    end if

! wait
#ifdef USE_NONBLK_COMM
    do j = 1, nb_proc-1
      call mpi_wait(req_r(j),istatus,ierr)
      call mpi_wait(req_s(j),istatus,ierr)
    end do
                                                  __TIMER_COMM_STOP(582)
#else
                                                  __TIMER_COMM_START_w_BARRIER(mpi_k_world,583)
     call MPI_ALLTOALL( tmp_s,datasize,mpi_double_precision, &
    &                   tmp_r,datasize,mpi_double_precision, &
    &                                  mpi_k_world,ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_ES_F_transpose :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 261, ierr)
     endif
                                                  __TIMER_COMM_STOP(583)
#endif

! rbuf
    do j = 1, nb_proc-1
       pe_r=mod(myrank_e-j+nb_proc,nb_proc)

          if(kimg_t == 1) then
             do iy = 1, np_e
                do ix = 1, nel_fs(pe_r)
                   bpr_l(iy,nis_fs(pe_r)-1+ix,ik) = tmp_r(ix,iy,1,pe_r)
                end do
             end do
          else
             do iy = 1, np_e
                do ix = 1, nel_fs(pe_r)
                   bpr_l(iy,nis_fs(pe_r)-1+ix,ik) = tmp_r(ix,iy,1,pe_r)
                   bpi_l(iy,nis_fs(pe_r)-1+ix,ik) = tmp_r(ix,iy,2,pe_r)
                end do
             end do
          end if

    end do

    deallocate(indexn)
!!$    call mpi_barrier(mpi_k_world(myrank_k),ierr)

    deallocate(tmp_s,tmp_r)
    deallocate(req_s,req_r)
    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
  end subroutine m_ES_F_transpose_back_r
#endif

#else
!! #ifdef TRANSPOSE else

#endif
!! #ifdef TRANSPOSE end


! -- GLOBAL --
  subroutine WSW_t_g(ik,j,mod_pot,fr,psi_t,LA,LB,LC,bpr_t,bpi_t)
    integer, intent(in)        :: ik,j,mod_pot
    real(kind=DP), intent(out) :: fr
    integer, intent(in)        :: LA,LB,LC
    real(kind=DP),intent(in)   :: psi_t(LA,LB,LC)
    real(kind=DP),intent(in),optional :: bpr_t(:,:),bpi_t(:,:)

    integer              :: ia,p,q, i, ig1
    real(kind=DP)        :: fr1
    integer :: id_sname = -1
    if(sw_timing_2ndlevel == ON) call tstatc0_begin('WSW_t_g ',id_sname)

                                                  __TIMER_SUB_START(504)
    fr = 0.d0
    if(mod_pot == VANDERBILT_TYPE) then
       if(k_symmetry(ik) == GAMMA) then
                                                  __TIMER_DO_START(548)
          do ia = 1, nac_p
             p = nlmta1_p(ia);     q = nlmta2_p(ia)
             fr = fr+fqwei_p(ia)*(bpr_t(p,j)*bpr_t(q,j))
          end do
                                                  __TIMER_DO_STOP(548)
       else
                                                  __TIMER_DO_START(549)
          do ia = 1, nac_p
             p = nlmta1_p(ia);     q = nlmta2_p(ia)
             fr = fr+fqwei_p(ia)*(bpr_t(p,j)*bpr_t(q,j)+bpi_t(p,j)*bpi_t(q,j))
          end do
                                                  __TIMER_DO_STOP(549)
       end if
    end if

    fr1 = 0.d0
    if(kimg == 1) then
                                                  __TIMER_DO_START(550)
       do i = 1, np_g1k(ik)                       ! MPI
          fr1 = fr1 + psi_t(i,j,1)*psi_t(i,j,1)
       end do
                                                  __TIMER_DO_STOP(550)
    else if(kimg == 2) then
       ig1 = 1; if(k_symmetry(ik) == GAMMA .and. myrank_e == 0) ig1 = 2
                                                  __TIMER_DO_START(551)
       do i = ig1, np_g1k(ik)
          fr1 = fr1 + psi_t(i,j,1)*psi_t(i,j,1) + psi_t(i,j,2)*psi_t(i,j,2)
       end do
                                                  __TIMER_DO_STOP(551)
       if(k_symmetry(ik) == GAMMA) fr1 = fr1*2.d0
       if(ig1 == 2) fr1 = fr1 + (psi_t(1,j,1)*psi_t(1,j,1))
    end if
    fr = fr+fr1
    if(npes > 1) then

                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,536)
       call mpi_allreduce(MPI_IN_PLACE,fr,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
                                                  __TIMER_COMM_STOP(536)
    end if
    fr = 1.d0/dsqrt(fr)
    if(ipri >= 2 .and. (ik == 1 .or. ik == 2)) then
       write(nfout,'(" ((WSW_t_g)) ik = ",i8," fr = ",f21.15)') ik, fr
    end if
                                                  __TIMER_SUB_STOP(504)
  if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
  end subroutine WSW_t_g

! -- GLOBAL --
  subroutine normalize_bp_and_psi_t_g(ik,ibo,fr,mod_pot,psi_t,LA,LB,LC &
       & ,bpr_t,bpi_t)
    integer, intent(in)        :: ik,ibo
    real(kind=DP),intent(in)   :: fr
    integer, intent(in)        :: mod_pot,LA,LB,LC
    real(kind=DP),intent(inout),dimension(LA,LB,LC) :: psi_t
    real(kind=DP),intent(inout),optional :: bpr_t(:,:),bpi_t(:,:)

    integer                    :: ia,ri
    integer :: id_sname = -1

    if(sw_timing_2ndlevel == ON) call tstatc0_begin('WSW_t_g ',id_sname)
                                                  __TIMER_SUB_START(506)
                                                  __TIMER_DO_START(554)
    if(mod_pot == VANDERBILT_TYPE) then
       if(k_symmetry(ik) == GAMMA) then
          do ia = 1, np_fs
             bpr_t(ia,ibo) = fr*bpr_t(ia,ibo)
          end do
       else
          do ia = 1, np_fs
             bpr_t(ia,ibo) = fr*bpr_t(ia,ibo)
             bpi_t(ia,ibo) = fr*bpi_t(ia,ibo)
          end do
       end if
    end if
                                                  __TIMER_DO_STOP(554)
                                                  __TIMER_DO_START(555)
    do ri = 1, kimg
!-F         psi_t(1:np_g1k(ik),ibo,ri) = fr * psi_t(1:np_g1k(ik),ibo,ri)
       psi_t(1:np_g1k(ik),ibo,ri) = fr * psi_t(1:np_g1k(ik),ibo,ri)
    end do
                                                  __TIMER_DO_STOP(555)
    if(ipri >= 2 .and. (ik == 1 .or. ik == 2)) then
       write(nfout,'(" ((normalize_bp_and_psi_t_g)) ik = ",i8," ibo = ",i8," fr = ",f21.15)') ik,ibo,fr
       write(nfout,'(6d14.6)') (psi_t(ia,ibo,1),ia=1, 6)
       if(kimg == 2) write(nfout,'(3d20.10)') (psi_t(ia,ibo,kimg),ia=1, 6)
    end if
                                                  __TIMER_SUB_STOP(506)
    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
  end subroutine normalize_bp_and_psi_t_g

! -- GLOBAL --
  subroutine cp_bpr2bpi_g(kimg_t,i,bpr_t,bpi_t)
    integer, intent(in) :: kimg_t,i
    real(kind=DP),intent(in)           :: bpr_t(:,:)
    real(kind=DP),intent(in),optional  :: bpi_t(:,:)
    integer ::  nel
    integer :: id_sname = -1
    if(sw_timing_2ndlevel == ON) call tstatc0_begin('cp_bpr2bpi_g ',id_sname)
    nel = min(np_fs,np_g1k_x)
    if(kimg_t == 1) then
       bp_ir(1:nel) = bpr_t(1:nel,i)
    else
       bp_ir(1:nel) = bpr_t(1:nel,i)
       bp_ii(1:nel) = bpi_t(1:nel,i)
    end if
    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
  end subroutine cp_bpr2bpi_g

! -- GLOBAL --
  subroutine W1SW2_t_r_g(ik,i,NB_END,mod_pot,phi_t,LA,LB,LC &
!!$#ifndef MGS_DGEMM
       & ,phifr_t,phifi_t &
!!$#endif
       & )
    integer, intent(in)        :: ik,i, NB_END,mod_pot
    integer, intent(in)        :: LA,LB,LC
    real(kind=DP),intent(in)          :: phi_t(LA,LB,LC)
!!$#ifndef MGS_DGEMM
    real(kind=DP),intent(in),optional :: phifr_t(:,:),phifi_t(:,:)
!!$#endif

! Coded by T. Yamasaki in April 2006
! Revised according to the RIKEN phase tuning project 2009, 13 Sep 2009
! This subroutine is moved out from subroutine mgs_4_each_k_G() in a contained state.

    integer       :: j, ia, jto, p, q,  kimg_t
    real(kind=DP) :: ar, ai
    real(kind=DP), allocatable, dimension(:,:)  :: p1Sp2_t2, p1Sp2_t1  ! MPI
    character*4 F_RSVTASK
    integer       :: nt, mpant, mmdnt, ipar, ist, ied, mm

    integer :: id_sname = -1, id_sname2 = -1
    integer :: myrank_g_common
#ifndef SX
    integer :: n_unroll, jmax, ia_start
    integer :: ibsize,ibl1,ibl2
    integer :: ncache
                                                  __TIMER_SUB_START(507)
    ncache = (m_CtrlP_cachesize()*1024)*3/4
#endif

    if(sw_timing_2ndlevel == ON) call tstatc0_begin('W1SW2_t_r_g ',id_sname)

    myrank_g_common = myrank_e


    if(nrank_e > 1) then
       if((k_symmetry(ik) == GAMMA .and. kimg==2) .or. kimg == 1) then
          kimg_t = 1
       else
          kimg_t = 2
       end if
       p = NB_END-i
       allocate(p1Sp2_t2(p,kimg_t)); p1Sp2_t2 = 0.d0                  ! MPI
       allocate(p1Sp2_t1(p,kimg_t)); p1Sp2_t1 = 0.d0                  ! MPI
    end if

    p1Sp2 = 0.d0

    if(mod_pot == VANDERBILT_TYPE) then
       if(kimg == 1) then
#ifndef SX
! NEC tune ------------------------------->
          if(ncache.eq.0) then
             ibsize=nac_p
!f
             if(ibsize == 0) ibsize = 1
          else
             ibsize=ncache/(8*(NB_END*2+3))
          endif
                                                  __TIMER_DO_START(556)
          do ibl1=1,nac_p,ibsize
             ibl2=ibl1+ibsize-1
             if(ibl2.gt.nac_p) ibl2=nac_p
! NEC tune <-------------------------------
#endif
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(ia,ar,ai,p,q)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
#ifndef ES
!CDIR OUTERUNROLL=4
#endif
#endif
!OCL NOFLTLD
             do j = i+1, NB_END
#ifdef SX
! NEC tune ------------------------------->
                ar = 0.d0
                do ia = 1, nac_p
#ifdef MGS_DGEMM
                   ar = ar + fqwei_p(ia)*(bpr_tw1(ia,i)*bpr_tw2(ia,j)+bpi_tw1(ia,i)*bpi_tw2(ia,j))
#else
                   p = nlmta1_p(ia);            q = nlmta2_p(ia)
                   ar = ar + fqwei_p(ia)*(bp_ir(p)*phifr_t(q,j)+bp_ii(p)*phifi_t(q,j))
#endif
                end do
                p1Sp2(j,1) = ar
#else
! NEC tune <-------------------------------
                ar = p1Sp2(j,1)
                do ia = ibl1, ibl2
#ifdef MGS_DGEMM
                   ar = ar + fqwei_p(ia)*(bpr_tw1(ia,i)*bpr_tw2(ia,j)+bpi_tw1(ia,i)*bpi_tw2(ia,j))
#else
                   p = nlmta1_p(ia);            q = nlmta2_p(ia)
                   ar = ar + fqwei_p(ia)*(bp_ir(p)*phifr_t(q,j)+bp_ii(p)*phifi_t(q,j))
#endif
                end do
                p1Sp2(j,1) = ar
#endif
            end do
#ifndef SX
! NEC tune
         end do
                                                  __TIMER_DO_STOP(556)
#endif
      else if(kimg == 2) then
         if(k_symmetry(ik) == GAMMA) then
!!!$#ifdef NEC_TUNE_SMP
!!!$!CDIR PARALLEL DO private(ia,ar,p,q)
!!!$#endif
!!!$#ifdef VPP
!!!$*vocl loop, unroll(4)
!!!$#endif
!!!$!CDIR OUTERUNROLL=4
!!!$               do j = i+1, NB_END                                                ! MPI
!!!$                  ar = 0.d0
!!!$                  do ia = 1, nac_p
!!!$                     p = nlmta1_p(ia);         q = nlmta2_p(ia)
!!!$                     ar = ar + fqwei_p(ia)*(bpr_t(p,i)*bpr_t(q,j))
!!!$                  end do
!!!$                  p1Sp2(j,1) = ar
!!!$               end do
#ifdef NEC_TUNE_SMP
             call getenv('F_RSVTASK',F_RSVTASK)
             read (F_RSVTASK,'(i4)') nt
#else
             nt = 1
#endif
             mpant = (NB_END-i)/nt
             mmdnt = mod(NB_END-i,nt)

#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(mm,ist,ied,ia,j,p,q)
#endif
! NEC tune
!!#ifdef NEC_TUNE_SMP
#ifdef SX
             do ipar = 1, min(nt,NB_END-i)                                            ! SX
                IF (IPAR.LE.MMDNT) THEN                                               ! SX
                   MM = MPANT+1                                                       ! SX
                ELSE                                                                  ! SX
                   MM = MPANT                                                         ! SX
                ENDIF                                                                 ! SX
                IST = (IPAR-1)*MPANT + MIN(MMDNT+1,IPAR) + i                          ! SX
                IED = IST + mm - 1                                                    ! SX
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                do ia = 1, nac_p                                                      ! SX
                   do j = IST, IED
		     ! SX
#ifdef MGS_DGEMM
                      p1Sp2(j,1) = p1Sp2(j,1) + fqwei_p(ia)*(bpr_tw1(ia,i)*bpr_tw2(ia,j))   ! SX
#else
                      p = nlmta1_p(ia);         q = nlmta2_p(ia)                      ! SX
                      p1Sp2(j,1) = p1Sp2(j,1) + fqwei_p(ia)*(bp_ir(p)*phifr_t(q,j))   ! SX
#endif
                   end do                                                             ! SX
                end do                                                                ! SX
             end do                                                                   ! SX
#else
! NEC tune ------------------------------->
             if(ncache.eq.0) then
                ibsize=nac_p
!f
                if(ibsize == 0) ibsize = 1
             else
                ibsize=ncache/(8*(NB_END*1+2))
             endif
                                                  __TIMER_DO_START(557)
!OCL NOFLTLD
             do ibl1=1,nac_p,ibsize
                ibl2=ibl1+ibsize-1
                if(ibl2.gt.nac_p) ibl2=nac_p
                do j = i+1, NB_END
                   if(ibl1.eq.1)then
                      ar=0.0d0
                   else
                      ar=p1Sp2(j,1)
                   endif
                   do ia = ibl1, ibl2
#ifdef MGS_DGEMM
                      ar = ar + fqwei_p(ia)*(bpr_tw1(ia,i)*bpr_tw2(ia,j))
#else
                      p = nlmta1_p(ia);         q = nlmta2_p(ia)
                      ar = ar + fqwei_p(ia)*(bp_ir(p)*phifr_t(q,j))
#endif
                   end do
                   p1Sp2(j,1)=ar
                end do
             end do
                                                  __TIMER_DO_STOP(557)

! NEC tune <-------------------------------
#endif
          else      ! kimg==2 .and. k_symmetry(ik) /= GAMMA
#ifndef SX
! NEC tune ------------------------------->
             if(ncache.eq.0) then
                ibsize=nac_p
!f
                if(ibsize == 0) ibsize = 1
             else
                ibsize=ncache/(8*(NB_END*4+3))
             endif
#ifdef DEBUG_MGS
             if(ipri >= 1) then
                write(nfout,'(" ibsize = ",i8," <<W1SW2_t_r>>")') ibsize
             end if
#endif
                                                  __TIMER_DO_START(558)
             do ibl1=1,nac_p,ibsize
                ibl2=ibl1+ibsize-1
                if(ibl2.gt.nac_p) ibl2=nac_p
! NEC tune <-------------------------------
#endif

#ifdef NEC_TUNE_SMP
#ifdef MGS_DGEMM
!CDIR PARALLEL DO private(ia,ar,ai)
#else
!CDIR PARALLEL DO private(ia,ar,ai,p,q)
#endif
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
#ifndef ES
!CDIR OUTERUNROLL=4
#endif
#endif
                do j = i+1, NB_END

#ifdef SX
                   ar = 0.d0; ai = 0.d0
                   do ia = 1, nac_p
#else
                   if(ibl1.eq.1) then
                      ar = 0.d0; ai = 0.d0
                   else
                      ar = p1Sp2(j,1)
                      ai = p1Sp2(j,2)

                   end if
                   do ia = ibl1, ibl2
#endif
#ifdef MGS_DGEMM
                      ar = ar + fqwei_p(ia)*(bpr_tw1(ia,i)*bpr_tw2(ia,j)+bpi_tw1(ia,i)*bpi_tw2(ia,j))
                      ai = ai + fqwei_p(ia)*(bpr_tw1(ia,i)*bpi_tw2(ia,j)-bpi_tw1(ia,i)*bpr_tw2(ia,j))
#else
                      p = nlmta1_p(ia);         q = nlmta2_p(ia)
                      ar = ar + fqwei_p(ia)*(bp_ir(p)*phifr_t(q,j)+bp_ii(p)*phifi_t(q,j))
                      ai = ai + fqwei_p(ia)*(bp_ir(p)*phifi_t(q,j)-bp_ii(p)*phifr_t(q,j))
#endif
                   end do
                   p1Sp2(j,1) = ar;  p1Sp2(j,2) = ai
                end do
#ifndef SX
! NEC tune
             end do
#endif
                                                  __TIMER_DO_STOP(558)
          end if
       end if
    end if
! --- <Psi_i|Psi_j> ---
    if(kimg == 1) then
#ifndef SX
! NEC tune ------------------------------------------------------------->
       if(ncache.eq.0) then
          ibsize=np_g1k(ik)
       else
          ibsize=ncache/(8*(NB_END*1+1))
       endif
                                                  __TIMER_DO_START(559)
       do ibl1=1,np_g1k(ik),ibsize
          ibl2=ibl1+ibsize-1
          if(ibl2.gt.np_g1k(ik)) ibl2=np_g1k(ik)
! NEC tune <-------------------------------------------------------------
#endif
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(ia,ar,ai,p,q)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
#ifndef ES
!CDIR OUTERUNROLL=4
#endif
#endif
          do j = i+1, NB_END
#ifdef SX
! NEC tune ------------------------------------------------------------->
             do ia = 1, np_g1k(ik)                                           ! MPI
                p1Sp2(j,1) = p1Sp2(j,1) + psi_ir(ia)*phi_t(ia,j,1 ) ! MPI
#else
             ai = p1Sp2(j,1)
             do ia = ibl1, ibl2
                ai = ai + psi_ir(ia)*phi_t(ia,j,1 ) ! MPI
! NEC tune <-------------------------------------------------------------
#endif
             end do

#ifndef SX
! NEC tune
             p1Sp2(j,1)=ai
#endif

          end do

#ifndef SX
! NEC tune
       end do
#endif
                                                  __TIMER_DO_STOP(559)

    else if(kimg == 2) then
       if(k_symmetry(ik) == GAMMA) then
!!!$            if(mype /= 0) then
!!!$#ifdef NEC_TUNE_SMP
!!!$!CDIR PARALLEL DO private(ia,ar,ai,p,q)
!!!$#endif
!!!$#ifdef VPP
!!!$*vocl loop, unroll(4)
!!!$#endif
!!!$!CDIR OUTERUNROLL=4
!!!$               do j = i+1, NB_END
!!!$                  do ia = 1, np_g1k(ik)
!!!$                     ar  = psi_t(ia,i,1)
!!!$                     ai  = psi_t(ia,i,2)
!!!$                     p1Sp2(j,1) = p1Sp2(j,1)+(ar*psi_t(ia,j,1)+ai*psi_t(ia,j,2))*2.d0 ! MPI
!!!$                  end do
!!!$               end do
#ifdef NEC_TUNE_SMP
          call getenv('F_RSVTASK',F_RSVTASK)
          read (F_RSVTASK,'(i4)') nt
#else
          nt = 1
#endif

! NEC tune
!!#ifdef NEC_TUNE_SMP
#ifdef SX
          if( (NB_END-i)/nt .lt. 256 ) then

             if(myrank_g_common /= 0) then

                do j = i+1, NB_END
                   do ia = 1, np_g1k(ik)
                      ar  = psi_ir(ia)
                      ai  = psi_ii(ia)
                      p1Sp2(j,1) = p1Sp2(j,1)+(ar*phi_t(ia,j,1)+ai*phi_t(ia,j,2))*2.d0
                   end do
                end do
             else if(myrank_g_common == 0) then

                do j = i+1, NB_END
                   do ia = 2, np_g1k(ik)
                      ar  = psi_ir(ia)
                      ai  = psi_ii(ia)
                      p1Sp2(j,1) = p1Sp2(j,1)+(ar*phi_t(ia,j,1)+ai*phi_t(ia,j,2))*2.d0
                   end do
                   p1Sp2(j,1) = p1Sp2(j,1) + psi_ir(1)*phi_t(1,j,1)
                end do
             end if

          else

             mpant = (NB_END-i)/nt
             mmdnt = mod(NB_END-i,nt)
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(mm,ist,ied,ia,j,ar,ai)
#endif
                                                  __TIMER_DO_START(560)
             do ipar = 1, min(nt,NB_END-i)
                IF (IPAR.LE.MMDNT) THEN
                   MM = MPANT+1
                ELSE
                   MM = MPANT
                ENDIF
                IST = (IPAR-1)*MPANT + MIN(MMDNT+1,IPAR) + i
                IED = IST + mm - 1
                if(myrank_g_common /= 0) then
                     !     write(*,*)'ipar, IST, IED=', ipar, IST, IED
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                   do ia = 1, np_g1k(ik)
                      do j = IST, IED
                         ar  = psi_ir(ia)
                         ai  = psi_ii(ia)
                         p1Sp2(j,1) = p1Sp2(j,1)+(ar*phi_t(ia,j,1)+ai*phi_t(ia,j,2))*2.d0
                      end do
                   end do
                else if(myrank_g_common == 0) then
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                   do ia = 2, np_g1k(ik)
                      do j = IST, IED
                         ar  = psi_ir(ia)
                         ai  = psi_ii(ia)
                         p1Sp2(j,1) = p1Sp2(j,1)+(ar*phi_t(ia,j,1)+ai*phi_t(ia,j,2))*2.d0
                      end do
                   enddo
                   do j = IST, IED
                      p1Sp2(j,1) = p1Sp2(j,1) + psi_ir(1)*phi_t(1,j,1)
                   end do
                end if
             end do
                                                  __TIMER_DO_STOP(560)
          end if
! NEC tune ------------------------------------>
#else
!!$            if(ncache.eq.0) then
!!$               ibsize=np_g1k(ik)
!!$            else
!!$               ibsize=ncache/(8*(NB_END*2+2))
!!$            endif

!!$            allocate(psi_ir(np_g1k(ik)))
!!$            allocate(psi_ii(np_g1k(ik)))
!!$            do ia = 1, np_g1k(ik)
!!$               psi_ir(ia) = psi_t(ia,i,1)
!!$               psi_ii(ia) = psi_t(ia,i,2)
!!$            end do
!!$
!!$            n_unroll = 4
!!$            jto = (NB_END-i)/n_unroll
!!$            jmax = i + n_unroll*jto

          if(sw_timing_2ndlevel == ON) call tstatc0_begin('W1SW2_t_r(core) ',id_sname2)
          if(myrank_g_common /= 0) then
             ia_start = 1
          else if(myrank_g == 0) then
             ia_start = 2
          end if

                                                  __TIMER_DO_START(561)
          do j = i+1, NB_END
             ar=p1Sp2(j,1)
             do ia = ia_start, np_g1k(ik)
                ar = ar+(psi_ir(ia)*phi_t(ia,j,1)+psi_ii(ia)*phi_t(ia,j,2))*2.d0
             end do
             if(myrank_g_common == 0) then
                ar = ar + psi_ir(1)*phi_t(1,j,1)
             endif
             p1Sp2(j,1)=ar
          end do
                                                  __TIMER_DO_STOP(561)

!!$            do j = i+1, jmax, n_unroll
!!$               do ia = ia_start, np_g1k(ik)
!!$                  ar = psi_ir(ia)*2.d0; ai = psi_ii(ia)*2.d0
!!$                  p1Sp2(j  ,1) = p1Sp2(j  ,1) + (ar*psi_t(ia,j  ,1)+ai*psi_t(ia,j  ,2))
!!$                  p1Sp2(j+1,1) = p1Sp2(j+1,1) + (ar*psi_t(ia,j+1,1)+ai*psi_t(ia,j+1,2))
!!$                  p1Sp2(j+2,1) = p1Sp2(j+2,1) + (ar*psi_t(ia,j+2,1)+ai*psi_t(ia,j+2,2))
!!$                  p1Sp2(j+3,1) = p1Sp2(j+3,1) + (ar*psi_t(ia,j+3,1)+ai*psi_t(ia,j+3,2))
!!$               end do
!!$            end do
!!$            do j = jmax+1, neg
!!$               do ia = ia_start, np_g1k(ik)
!!$                  ar = psi_ir(ia); ai = psi_ii(ia)
!!$                  p1Sp2(j  ,1) = p1Sp2(j  ,1) + (ar*psi_t(ia,j,  1)+ai*psi_t(ia,j,  2))*2.d0
!!$               end do
!!$            end do
!!$            if(mype == 0) then
!!$               do j = i+1, neg
!!$                  p1Sp2(j,1) = p1Sp2(j,1) + psi_ir(1)*psi_t(1,j, 1)
!!$               end do
!!$            end if

          if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname2)
#endif
! NEC tune <------------------------------------
       else    ! kimg==2 .and. k_symmetry(ik) /= GAMMA
#ifdef SX
! NEC tune ------------------------------------>
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(ia,ar,ai,p,q)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
#ifndef ES
!CDIR OUTERUNROLL=4
#endif
#endif
                                                  __TIMER_DO_START(562)
          do j = i+1, NB_END
             do ia = 1, np_g1k(ik)                                          ! MPI
                ar  = psi_ir(ia)
                ai  = psi_ii(ia)
                p1Sp2(j,1) = p1Sp2(j,1)+ar*phi_t(ia,j,1)+ai*phi_t(ia,j,2) ! MPI
                p1Sp2(j,2) = p1Sp2(j,2)+ar*phi_t(ia,j,2)-ai*phi_t(ia,j,1) ! MPI
             end do
          end do
                                                  __TIMER_DO_STOP(562)
#else
!!$            if(ncache.eq.0) then
          ibsize=np_g1k(ik)
!!$            else
!!$!!$               ibsize=ncache/(8*(neg*2+2))
!!$               ibsize=ncache/(8*(NB_END*2+2))
!!$            endif
!!$            allocate(psi_ir(np_g1k(ik)))
!!$            allocate(psi_ii(np_g1k(ik)))
!!$            do ia = 1, np_g1k(ik)
!!$               psi_ir(ia) = psi_iri(ia,1)
!!$               psi_ii(ia) = psi_iri(ia,2)
!!$            end do

                                                  __TIMER_DO_START(563)
          do ibl1=1,np_g1k(ik),ibsize
             ibl2=ibl1+ibsize-1
             if(ibl2.gt.np_g1k(ik)) ibl2=np_g1k(ik)
             if(sw_timing_2ndlevel == ON) call tstatc0_begin('W1SW2_t_r(core) ',id_sname2)
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO private(ia,ar,ai,p,q)
#endif
#ifdef NEC_TUNE2
#ifndef ES
!CDIR OUTERUNROLL=4
#define _DIRECTIVE_UNROLLING_
#endif
#endif
#ifdef _DIRECTIVE_UNROLLING_
             do j = i+1, NB_END
                ar=p1Sp2(j,1)
                ai=p1Sp2(j,2)
                do ia = ibl1, ibl2
                   ar = ar+psi_ir(ia)*phi_t(ia,j,1)+psi_ii(ia)*phi_t(ia,j,2)
                   ai = ai+psi_ir(ia)*phi_t(ia,j,2)-psi_ii(ia)*phi_t(ia,j,1)
                end do
                p1Sp2(j,1)=ar
                p1Sp2(j,2)=ai
             end do
#else
             n_unroll = 4
             jto = (NB_END-i)/n_unroll
             jmax = i + n_unroll*jto
             do j = i+1, jmax, n_unroll
                do ia = ibl1, ibl2
                   ar = psi_ir(ia); ai = psi_ii(ia)
                   p1Sp2(j,1)   = p1Sp2(j,1)  +ar*phi_t(ia,j,1)  +ai*phi_t(ia,j,2)
                   p1Sp2(j,2)   = p1Sp2(j,2)  +ar*phi_t(ia,j,2)  -ai*phi_t(ia,j,1)
                   p1Sp2(j+1,1) = p1Sp2(j+1,1)+ar*phi_t(ia,j+1,1)+ai*phi_t(ia,j+1,2)
                   p1Sp2(j+1,2) = p1Sp2(j+1,2)+ar*phi_t(ia,j+1,2)-ai*phi_t(ia,j+1,1)
                   p1Sp2(j+2,1) = p1Sp2(j+2,1)+ar*phi_t(ia,j+2,1)+ai*phi_t(ia,j+2,2)
                   p1Sp2(j+2,2) = p1Sp2(j+2,2)+ar*phi_t(ia,j+2,2)-ai*phi_t(ia,j+2,1)
                   p1Sp2(j+3,1) = p1Sp2(j+3,1)+ar*phi_t(ia,j+3,1)+ai*phi_t(ia,j+3,2)
                   p1Sp2(j+3,2) = p1Sp2(j+3,2)+ar*phi_t(ia,j+3,2)-ai*phi_t(ia,j+3,1)
                end do
             end do
             do j = jmax+1, NB_END
                do ia = 1, ibl2-ibl1+1
                   ar = psi_ir(ia); ai=psi_ii(ia)
                   p1Sp2(j,1)   = p1Sp2(j,1)   + ar*phi_t(ia,j,1)+ai*phi_t(ia,j,2)
                   p1Sp2(j,2)   = p1Sp2(j,2)   + ar*phi_t(ia,j,2)-ai*phi_t(ia,j,1)
                end do
             end do
#endif
             if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname2)
          end do
                                                  __TIMER_DO_STOP(563)
!!$          deallocate(psi_ii, psi_ir)
! NEC tune <------------------------------------
#endif
       end if
    end if

    if(nrank_e > 1 .and. NB_END > i ) then
                                                  __TIMER_COMM_START(537)
       do q = 1, kimg_t
          do j = 1, NB_END-i
             p = j+i
             p1Sp2_t1(j,q) = p1Sp2(p,q)
          end do
       end do
                                                  __TIMER_COMM_STOP(537)
                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,538)
!!$         p = (neg-i)*kimg_t
       call mpi_allreduce(p1Sp2_t1, p1Sp2_t2,(NB_END-i)*kimg_t,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
!!$         call mpi_allreduce(p1Sp2_t1, p1Sp2_t2,(neg-i)*kimg_t,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
                                                  __TIMER_COMM_STOP(538)
                                                  __TIMER_COMM_START(539)
       do q = 1, kimg_t
          do j = 1, NB_END-i
             p = j+i
             p1Sp2(p,q) = p1Sp2_t2(j,q)
          end do
       end do
!!$         p1Sp2(:,1:kimg_t) = p1Sp2_t2(:,1:kimg_t)
                                                  __TIMER_COMM_STOP(539)
    end if

    if(ipri >= 2 .and. ik == 1) then
       write(nfout,'(" <<W1SW2_t>> ik = ",i9,"  i = ",i3)') ik,i
       write(nfout,'(" (real) ",6d11.3, 99(/8x,6d11.3))') (p1Sp2(j,1),j=i+1,NB_END)
! === DEBUG by tkato 2013/11/05 ================================================
       if(nrank_e > 1) then
! ==============================================================================
       if(kimg_t == 2) &
            & write(nfout,'(" (imag) ",6d11.3, 99(/8x,6d11.3))') (p1Sp2(j,kimg_t),j=i+1,NB_END)
! === DEBUG by tkato 2013/11/05 ================================================
       end if
! ==============================================================================
    end if
    if(nrank_e > 1) deallocate(p1Sp2_t1,p1Sp2_t2)
    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(507)
  end subroutine W1SW2_t_r_g

! -- GLOBAL --
!-$ tune FUJITSU for block ----->>
!     13th Sep 2009
#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
#ifdef MGS_DGEMM
  subroutine W1SW2_t_r_block_g(ik,i,i2,p1Sp2_NB,NB_END,NB_END2,kimg_t_wk &
       &  ,   mod_pot,psi_t,psi_t_dia,LA,LB,LC,bpr_tw1_dia,bpi_tw1_dia)
    integer, intent(in)        :: ik,i, i2, NB_END, NB_END2, kimg_t_wk
    integer, intent(in)        :: mod_pot,LA,LB,LC
    real(kind=DP), intent(out), dimension(NB,NB,kimg_t_wk) :: p1Sp2_NB
    real(kind=DP), intent(in), dimension(LA,LB,LC)      :: psi_t
    real(kind=DP), intent(in), dimension(LA,NB,LC)      :: psi_t_dia
    real(kind=DP), intent(in), dimension(nac_p,*),optional :: bpr_tw1_dia,bpi_tw1_dia

    real(kind=DP), allocatable, dimension(:,:) :: bpr_tw1_BLAS, bpi_tw1_BLAS ! work

    integer       :: iblk, jblk, i_NB, j_NB
    integer       :: ii, ia,  M, N, LDA, LDB, p

    integer :: id_sname = -1, id_sname2 = -1
                                                  __TIMER_SUB_START(509)
    if(sw_timing_2ndlevel == ON) call tstatc0_begin('W1SW2_t_r_block_g ',id_sname)

    p1Sp2_NB = 0.0d0
    M = NB_END2-i2+1; N = NB_END-i+1
    if(mod_pot == VANDERBILT_TYPE .and. nac_p > 0) then
       allocate(bpr_tw1_BLAS(nac_p, NB))
       bpr_tw1_BLAS = 0.0
                                                  __TIMER_DO_START(570)
!OCL NOFLTLD
       do iblk = i, NB_END
          ii = iblk - i + 1
          do ia = 1, nac_p
             bpr_tw1_BLAS(ia,ii) = fqwei_p(ia)*bpr_tw1_dia(ia,ii)
          end do
       end do
                                                  __TIMER_DO_STOP(570)
       if(k_symmetry(ik) /= GAMMA) then
          allocate(bpi_tw1_BLAS(nac_p,NB))
          bpi_tw1_BLAS = 0.d0
                                                  __TIMER_DO_START(571)
          do iblk = i, NB_END
             ii = iblk - i + 1
             do ia = 1, nac_p
                bpi_tw1_BLAS(ia,ii) = fqwei_p(ia)*bpi_tw1_dia(ia,ii)
             end do
          end do
                                                  __TIMER_DO_STOP(571)
       end if

       if(kimg==1) then
                                                  __TIMER_DGEMM_START(515)
          call DGEMM__('T','N',M,N,nac_p,1.0d0,bpr_tw2(1,i2),nac_p,bpr_tw1_BLAS,nac_p,1.0d0,p1Sp2_NB(1,1,1),NB)
          call DGEMM__('T','N',M,N,nac_p,1.0d0,bpi_tw2(1,i2),nac_p,bpi_tw1_BLAS,nac_p,1.0d0,p1Sp2_NB(1,1,1),NB)
                                                  __TIMER_DGEMM_STOP(515)
!!$                     ar = ar + fqwei_p(ia)*(bpr_t(p,i)*bpr_t(q,j)+bpi_t(p,i)*bpi_t(q,j))
       else if(kimg==2 .and. k_symmetry(ik) == GAMMA) then
                                                  __TIMER_DGEMM_START(516)
          call DGEMM__('T','N',M,N,nac_p,1.0d0,bpr_tw2(1,i2),nac_p,bpr_tw1_BLAS,nac_p,1.0d0,p1Sp2_NB(1,1,1),NB)
!!$                     ar = ar + fqwei_p(ia)*(bpr_tw1(ia,i)*bpr_tw2(ia,j))
                                                  __TIMER_DGEMM_STOP(516)
       else if(kimg==2 .and. k_symmetry(ik) /= GAMMA) then
                                                  __TIMER_DGEMM_START(517)
          call DGEMM__('T','N',M,N,nac_p, 1.0d0,bpr_tw2(1,i2),nac_p,bpr_tw1_BLAS,nac_p,1.0d0,p1Sp2_NB(1,1,1),NB)
          call DGEMM__('T','N',M,N,nac_p, 1.0d0,bpi_tw2(1,i2),nac_p,bpr_tw1_BLAS,nac_p,1.0d0,p1Sp2_NB(1,1,2),NB)
          call DGEMM__('T','N',M,N,nac_p, 1.0d0,bpi_tw2(1,i2),nac_p,bpi_tw1_BLAS,nac_p,1.0d0,p1Sp2_NB(1,1,1),NB)
          call DGEMM__('T','N',M,N,nac_p,-1.0d0,bpr_tw2(1,i2),nac_p,bpi_tw1_BLAS,nac_p,1.0d0,p1Sp2_NB(1,1,2),NB)
!!$                     ar = ar + fqwei_p(ia)*(bpr_tw1(ia,i)*bpr_tw2(ia,j)+bpi_tw1(ia,i)*bpi_tw2(ia,j))
!!$                     ai = ai + fqwei_p(ia)*(bpr_tw1(ia,i)*bpi_tw2(ia,j)-bpi_tw1(ia,i)*bpr_tw2(ia,j))
                                                  __TIMER_DGEMM_STOP(517)
       end if

       if(k_symmetry(ik) /= GAMMA) deallocate(bpi_tw1_BLAS)
       deallocate(bpr_tw1_BLAS)
    end if
! --------------------------------------
    if(sw_timing_2ndlevel == ON) call tstatc0_begin('W1SW2_t_r_block_g(core) ',id_sname2)
    LDA = np_g1k_x; LDB = np_g1k_x
    if(kimg == 1) then
                                                  __TIMER_DGEMM_START(518)
       call DGEMM__('T','N',M,N,np_g1k(ik),1.0d0,psi_t(1,i2,1),LDA,psi_t_dia(1,1,1), LDB,1.0d0,p1Sp2_NB(1,1,1),NB)
                                                  __TIMER_DGEMM_STOP(518)
    else if(kimg == 2 .and. k_symmetry(ik) == GAMMA ) then
                                                  __TIMER_DGEMM_START(519)
       call DGEMM__('T','N',M,N,np_g1k(ik),2.0d0,psi_t(1,i2,1),LDA,psi_t_dia(1,1,1),LDB,1.0d0,p1Sp2_NB(1,1,1),NB)
       call DGEMM__('T','N',M,N,np_g1k(ik),2.0d0,psi_t(1,i2,2),LDA,psi_t_dia(1,1,2),LDB,1.0d0,p1Sp2_NB(1,1,1),NB)
                                                  __TIMER_DGEMM_STOP(519)
       if (myrank_e == 0 ) then
                                                  __TIMER_DO_START(572)
          do iblk = i, NB_END                !for BLAS3
             i_NB = iblk - i + 1
             do jblk = i2, NB_END2
                j_NB = jblk - i2 + 1
!f                  p1Sp2_NB(j_NB,i_NB,1) = p1Sp2_NB(j_NB,i_NB,1) &
!f                       & - (psi_t(1,iblk,1)*psi_t(1,jblk,1)+psi_t(1,iblk,2)*psi_t(1,jblk,2)*2.d0)
                p1Sp2_NB(j_NB,i_NB,1) = p1Sp2_NB(j_NB,i_NB,1) &
                     & - (psi_t_dia(1,i_NB,1)*psi_t(1,jblk,1)+psi_t_dia(1,i_NB,2)*psi_t(1,jblk,2)*2.d0)
             end do
          end do
                                                  __TIMER_DO_STOP(572)
       end if
    else  ! kimg==2 .and. k_symmetry(ik) /= GAMMA
                                                  __TIMER_DGEMM_START(520)
       call DGEMM__('T','N',M,N,np_g1k(ik), 1.0d0,psi_t(1,i2,1),LDA,psi_t_dia(1,1,1),LDB,1.0d0,p1Sp2_NB(1,1,1),NB)
       call DGEMM__('T','N',M,N,np_g1k(ik), 1.0d0,psi_t(1,i2,2),LDA,psi_t_dia(1,1,2),LDB,1.0d0,p1Sp2_NB(1,1,1),NB)
       call DGEMM__('T','N',M,N,np_g1k(ik), 1.0d0,psi_t(1,i2,2),LDA,psi_t_dia(1,1,1),LDB,1.0d0,p1Sp2_NB(1,1,2),NB)
       call DGEMM__('T','N',M,N,np_g1k(ik),-1.0d0,psi_t(1,i2,1),LDA,psi_t_dia(1,1,2),LDB,1.0d0,p1Sp2_NB(1,1,2),NB)
                                                  __TIMER_DGEMM_STOP(570)
    end if
    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname2)

!-F    if(nrank_e > 1 ) then
!-F       if(ipri>=2) then
!-F          write(nfout,'(" NB_END, NB_END2, NB_END-i+1, NB_END2-i2+1, NB, kimg_t_wk = ",6i8)') &
!-F               & NB_END, NB_END2,  NB_END-i+1, NB_END2-i2+1, NB, kimg_t_wk
!-F       end if
!-F       call mpi_allreduce(MPI_IN_PLACE,p1Sp2_NB,NB*NB*kimg_t_wk,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
!-F    end if

    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(509)
  end subroutine W1SW2_t_r_block_g
#endif
!-$ tune FUJITSU for block <<-----
#endif

! -- GLOBAL --
#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
#ifdef MGS_DGEMM
    subroutine modify_bp_and_psi_t_r_blk_g(ik,i,i2,p1Sp2_NB,NB_END,NB_END2,kimg_t_wk &
	&   , mod_pot,psi_t,psi_t_dia,LA,LB,LC,bpr_t,bpi_t,bpr_t_dia,bpi_t_dia)
      integer, intent(in) :: ik,i,i2,mod_pot,LA,LB,LC
      real(kind=DP), intent(in),  dimension(NB,NB,kimg_t_wk) :: p1Sp2_NB
      real(kind=DP), intent(inout), dimension(LA,LB,LC) :: psi_t
      real(kind=DP), intent(in), dimension(LA,NB,LC)    :: psi_t_dia
      real(kind=DP), intent(inout), optional :: bpr_t(:,:),bpi_t(:,:)
      real(kind=DP), intent(in), optional :: bpr_t_dia(:,:),bpi_t_dia(:,:)

      integer, intent(in) :: NB_END, NB_END2, kimg_t_wk
      integer       ::  N, K, LDA, LDB, LDC
      integer :: id_sname = -1
                                                  __TIMER_SUB_START(510)
      if(i == neg) return

      if(sw_timing_2ndlevel == ON) call tstatc0_begin('modify_bp_and_psi_t_r_blk_g ',id_sname)

      N = NB_END2-i2+1
      K = NB_END-i+1

      if(mod_pot == VANDERBILT_TYPE .and. np_fs.gt.0) then
         LDA = np_fs_x; LDB = NB; LDC = np_fs_x
                                                  __TIMER_DGEMM_START(521)
         call DGEMM__('N','T',np_fs,N,K,-1.0d0,bpr_t_dia(1,1),LDA,p1Sp2_NB(1,1,1),LDB,1.0d0,bpr_t(1,i2),LDC)
                                                  __TIMER_DGEMM_STOP(521)

         if(kimg == 1 .or. (kimg==2 .and. k_symmetry(ik) /= GAMMA)) then
                                                  __TIMER_DGEMM_START(522)
            call DGEMM__('N','T',np_fs,N,K,-1.0d0,bpi_t_dia(1,1),LDA,p1Sp2_NB(1,1,1),LDB, 1.0d0,bpi_t(1,i2),LDC)
                                                  __TIMER_DGEMM_STOP(522)
         end if

         if(kimg==2 .and. k_symmetry(ik) /= GAMMA) then
                                                  __TIMER_DGEMM_START(523)
            call DGEMM__('N','T',np_fs,N,K, 1.0d0,bpi_t_dia(1,1),LDA,p1Sp2_NB(1,1,2),LDB,1.0d0,bpr_t(1,i2),LDC)
            call DGEMM__('N','T',np_fs,N,K,-1.0d0,bpr_t_dia(1,1),LDA,p1Sp2_NB(1,1,2),LDB,1.0d0,bpi_t(1,i2),LDC)
                                                  __TIMER_DGEMM_STOP(523)
         end if
      end if

      LDA = np_g1k_x; LDB = NB; LDC = np_g1k_x
      if(kimg == 1) then
                                                  __TIMER_DGEMM_START(524)
         call DGEMM__('N','T',np_g1k(ik),N,K,-1.d0,psi_t_dia(1,1,1),LDA,p1Sp2_NB(1,1,1),LDB,1.d0,psi_t(1,i2,1),LDC)
                                                  __TIMER_DGEMM_STOP(524)
      else if(kimg==2 .and. k_symmetry(ik) == GAMMA) then
                                                  __TIMER_DGEMM_START(525)
         call DGEMM__('N','T',np_g1k(ik),N,K,-1.d0,psi_t_dia(1,1,1),LDA,p1Sp2_NB(1,1,1),LDB,1.d0,psi_t(1,i2,1),LDC)
         call DGEMM__('N','T',np_g1k(ik),N,K,-1.d0,psi_t_dia(1,1,2),LDA,p1Sp2_NB(1,1,1),LDB,1.d0,psi_t(1,i2,2),LDC)
                                                  __TIMER_DGEMM_STOP(525)
      else  ! kimg == 2 .and. k_symmetry(ik) /= GAMMA
                                                  __TIMER_DGEMM_START(526)
         call DGEMM__('N','T',np_g1k(ik),N,K,-1.d0,psi_t_dia(1,1,1),LDA,p1Sp2_NB(1,1,1),LDB,1.d0,psi_t(1,i2,1),LDC)
         call DGEMM__('N','T',np_g1k(ik),N,K,-1.d0,psi_t_dia(1,1,2),LDA,p1Sp2_NB(1,1,1),LDB,1.d0,psi_t(1,i2,2),LDC)
         call DGEMM__('N','T',np_g1k(ik),N,K, 1.d0,psi_t_dia(1,1,2),LDA,p1Sp2_NB(1,1,2),LDB,1.d0,psi_t(1,i2,1),LDC)
         call DGEMM__('N','T',np_g1k(ik),N,K,-1.d0,psi_t_dia(1,1,1),LDA,p1Sp2_NB(1,1,2),LDB,1.d0,psi_t(1,i2,2),LDC)
                                                  __TIMER_DGEMM_STOP(526)
      end if
      if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(510)
    end subroutine modify_bp_and_psi_t_r_blk_g
!-$ tune FUJITSU for block <<-----
#endif
#endif

#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT
! -- GLOBAL --
  subroutine modify_bp_and_psi_t_r_g(ik,i,NB_END,mod_pot &
       & ,psi_t,LA,LB,LC,bpr_t,bpi_t)
! Coded by T. Yamasaki in April 2006
! Revised according to the RIKEN phase tuning project 2009, 13 Sep 2009
    integer, intent(in) :: ik,i, NB_END,mod_pot
    integer, intent(in) :: LA,LB,LC
    real(kind=DP),intent(inout) :: psi_t(LA,LB,LC)
    real(kind=DP),intent(inout),optional :: bpr_t(:,:),bpi_t(:,:)

    integer             :: j,ia
    real(kind=DP) :: sr, si
#ifndef SX
    real(kind=DP) :: ar, ai
#endif
#ifdef MGS_DGEMM
    integer :: p, q
#endif
    integer :: id_sname = -1
                                                  __TIMER_SUB_START(508)
    if(i == neg) return
    if(sw_timing_2ndlevel == ON) call tstatc0_begin('modify_bp_and_psi_t_r_g ',id_sname)

    if(mod_pot == VANDERBILT_TYPE) then
       if(kimg == 1) then
#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!!$!CDIR NOSYNC
!!$!CDIR CONCUR(BY=1)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                                                  __TIMER_DO_START(564)
          do j = i+1, NB_END
#ifdef SX
             do ia = 1, np_fs
                bpr_t(ia,j) = bpr_t(ia,j) - p1Sp2(j,1)*bp_ir(ia)
                bpi_t(ia,j) = bpi_t(ia,j) - p1Sp2(j,1)*bp_ii(ia)
             end do
#else
             ar=p1Sp2(j,1)
             do ia = 1, np_fs
                bpr_t(ia,j) = bpr_t(ia,j) - ar*bp_ir(ia)
                bpi_t(ia,j) = bpi_t(ia,j) - ar*bp_ii(ia)
             end do
#endif

          end do  ! j-loop
                                                  __TIMER_DO_STOP(564)
       else if(kimg == 2) then
          if(k_symmetry(ik) == GAMMA) then
#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!!$!CDIR NOSYNC
!!$!CDIR CONCUR(BY=1)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                                                  __TIMER_DO_START(565)
             do j = i+1, NB_END
#ifdef SX
! NEC tune ---------------------------------------------------------->
                do ia = 1, np_fs
                   bpr_t(ia,j) = bpr_t(ia,j) - p1Sp2(j,1)*bp_ir(ia)
                end do
#else
                ar=p1Sp2(j,1)
                do ia = 1, np_fs
                   bpr_t(ia,j) = bpr_t(ia,j) - ar*bp_ir(ia)
                end do
! NEC tune <----------------------------------------------------------
#endif
             end do
                                                  __TIMER_DO_STOP(565)
          else  ! kimg==2, k_symmetry(ik) /= GAMMA
#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!!$!CDIR NOSYNC
!!$!CDIR CONCUR(BY=1)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                                                  __TIMER_DO_START(566)
             do j = i+1, NB_END
#ifdef SX
! NEC tune ---------------------------------------------------------->
                do ia = 1, np_fs
                   sr  =  bp_ir(ia);     si  =  bp_ii(ia)
                   bpr_t(ia,j) = bpr_t(ia,j) - p1Sp2(j,1)*sr+p1Sp2(j,2)*si
                   bpi_t(ia,j) = bpi_t(ia,j) - p1Sp2(j,1)*si-p1Sp2(j,2)*sr
                end do
#else
                ar=p1Sp2(j,1)
                ai=p1Sp2(j,2)
                do ia = 1, np_fs
                   sr  =  bp_ir(ia);     si  =  bp_ii(ia)
                   bpr_t(ia,j) = bpr_t(ia,j) - ar*sr+ai*si
                   bpi_t(ia,j) = bpi_t(ia,j) - ar*si-ai*sr
                end do
#endif
             end do
                                                  __TIMER_DO_STOP(566)
! NEC tune <----------------------------------------------------------
          end if
       end if
    end if

    if(kimg == 1) then
#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!!$!CDIR NOSYNC
!!$!CDIR CONCUR(BY=1)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                                                  __TIMER_DO_START(567)
       do j = i+1,NB_END
#ifdef SX
! NEC tune ---------------------------------------------------------->
          do ia = 1, np_g1k(ik)                 ! MPI
!!$             psi_t(ia,j,1) = psi_t(ia,j,1) - p1Sp2(j,1)*psi_t(ia,i,1 )
             psi_t(ia,j,1) = psi_t(ia,j,1) - p1Sp2(j,1)*psi_ir(ia )
          end do
#else
          ar=p1Sp2(j,1)
          do ia = 1, np_g1k(ik)
!!$             psi_t(ia,j,1) = psi_t(ia,j,1) - ar*psi_t(ia,i,1 )
             psi_t(ia,j,1) = psi_t(ia,j,1) - ar*psi_ir(ia)
          end do
! NEC tune <----------------------------------------------------------
#endif
       end do
                                                  __TIMER_DO_STOP(567)
    else if(kimg == 2) then
       if(k_symmetry(ik) == GAMMA) then
#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!!$!CDIR NOSYNC
!!$!CDIR CONCUR(BY=1)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                                                  __TIMER_DO_START(568)
          do j = i+1,NB_END
#ifdef SX
! NEC tune ---------------------------------------------------------->
             do ia = 1, np_g1k(ik)
!!$                psi_t(ia,j,1) = psi_t(ia,j,1) - p1Sp2(j,1)*psi_t(ia,i,1 )
!!$                psi_t(ia,j,2) = psi_t(ia,j,2) - p1Sp2(j,1)*psi_t(ia,i,2 )
                psi_t(ia,j,1) = psi_t(ia,j,1) - p1Sp2(j,1)*psi_ir(ia)
                psi_t(ia,j,2) = psi_t(ia,j,2) - p1Sp2(j,1)*psi_ii(ia)
             end do
#else
             ar=p1Sp2(j,1)
             do ia = 1, np_g1k(ik)
!!$                psi_t(ia,j,1) = psi_t(ia,j,1) - ar*psi_t(ia,i,1 )
!!$                psi_t(ia,j,2) = psi_t(ia,j,2) - ar*psi_t(ia,i,2 )
                psi_t(ia,j,1) = psi_t(ia,j,1) - ar*psi_ir(ia)
                psi_t(ia,j,2) = psi_t(ia,j,2) - ar*psi_ii(ia)
             end do
! NEC tune <----------------------------------------------------------
#endif
          end do
                                                  __TIMER_DO_STOP(568)
       else   ! kimg==2 .and. k_symmetry(ik) /= GAMMA
#ifdef NEC_TUNE_SMP
!!$!CDIR PARALLEL DO PRIVATE(jto,ia,sr,si)
!!$!CDIR NOSYNC
!!$!CDIR CONCUR(BY=1)
#endif
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
                                                  __TIMER_DO_START(569)
          do j = i+1,NB_END
#ifdef SX
! NEC tune ---------------------------------------------------------->
             do ia = 1, np_g1k(ik)                ! MPI
!!$                sr  =  psi_t(ia,i,1 );   si  =  psi_t(ia,i,2 )
                sr  =  psi_ir(ia);   si  =  psi_ii(ia)
                psi_t(ia,j,1) = psi_t(ia,j,1) - p1Sp2(j,1)*sr+p1Sp2(j,2)*si
                psi_t(ia,j,2) = psi_t(ia,j,2) - p1Sp2(j,1)*si-p1Sp2(j,2)*sr
             end do
#else
             ar=p1Sp2(j,1)
             ai=p1Sp2(j,2)
             do ia = 1, np_g1k(ik)
!!$                sr  =  psi_t(ia,i,1 );   si  =  psi_t(ia,i,2 )
                sr  =  psi_ir(ia);   si  =  psi_ii(ia)
                psi_t(ia,j,1) = psi_t(ia,j,1) - ar*sr+ai*si
                psi_t(ia,j,2) = psi_t(ia,j,2) - ar*si-ai*sr
             end do
#endif
! NEC tune <----------------------------------------------------------
          end do
                                                  __TIMER_DO_STOP(569)
       end if
    end if
    if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(508)
  end subroutine modify_bp_and_psi_t_r_g
#endif

#ifndef TRANSPOSE_WITHOUT_REARRANGEMENT

! -- GLOBAL --
  subroutine cp_psi2psii_g(ik,i)
    integer, intent(in) :: ik,i
    psi_ir(1:np_g1k(ik)) = psi_t(1:np_g1k(ik),i,1)
    if(kimg==2) psi_ii(1:np_g1k(ik)) = psi_t(1:np_g1k(ik),i,kimg)
  end subroutine cp_psi2psii_g
#endif

  ! -- GLOBAL --
#ifdef MGS_DGEMM
  subroutine cp_psi_bpri2dias_g(ik,i,kimg_t_wk,NB,psi_t_dia,LA,LB,LC,MA &
       &                           ,bpr_t_dia,bpi_t_dia)
    integer, intent(in) :: ik,i,kimg_t_wk,NB,LA,LB,LC,MA
    real(kind=DP),intent(out),dimension(LA,LB,LC) :: psi_t_dia
    real(kind=DP),intent(out),dimension(MA,NB),optional  :: bpr_t_dia,bpi_t_dia

    integer :: ri, iy, ix

    do ri = 1, kimg
       do iy = 1, NB
          if(iy+i-1 > neg) cycle
          do ix = 1, np_g1k(ik)
             psi_t_dia(ix,iy,ri) = psi_t(ix,iy+i-1,ri)
    end do; end do; end do;
    if(np_g1k(ik) < LA) psi_t_dia(LA,:,:) = 0.d0

    if(kimg_t_wk == 1) then
       if(present(bpr_t_dia)) then
          do iy = 1, NB
             if(iy+i-1 > neg) cycle
             do ix = 1, np_fs
                bpr_t_dia(ix,iy) = bpr_t(ix,iy+i-1)
          end do;     end do;
          if(np_fs_x > MA) bpr_t_dia(MA,:) = 0.d0
       end if
    else
       if(present(bpr_t_dia) .and. present(bpi_t_dia)) then
          do iy = 1, NB
             if(iy+i-1 > neg) cycle
             do ix = 1, np_fs
                bpr_t_dia(ix,iy) = bpr_t(ix,iy+i-1)
                bpi_t_dia(ix,iy) = bpi_t(ix,iy+i-1)
          end do;     end do;
          if(np_fs_x > MA) then
             bpr_t_dia(MA,:) = 0.d0
             bpi_t_dia(MA,:) = 0.d0
          end if
       end if
    end if
  end subroutine cp_psi_bpri2dias_g
#endif

  subroutine mgs_phi2wf_each_k_G(ik,phi_l,mode,bsdr_l,bsdi_l,mod_pot)
    integer, intent(in)             :: ik
    real(kind=DP), intent(inout)    :: phi_l(kg1,np_e,ik:ik,kimg)
    integer, intent(in)             :: mode
    real(kind=DP), optional,intent(inout), dimension(np_e,nlmta,ik:ik) :: bsdr_l,bsdi_l
    integer, intent(in)             :: mod_pot
    integer               :: i

    real(kind=DP) :: fr
    integer :: kimg_t_wk

    real(kind=DP),allocatable,dimension(:,:)    :: phifr_t, phifi_t
    real(kind=DP),allocatable, dimension(:,:,:) :: phi_t
#ifdef MGS_DGEMM
    real(kind=DP),allocatable,dimension(:,:,:) :: psi_t_dia
    real(kind=DP),allocatable,dimension(:,:)   :: bpr_t_dia, bpi_t_dia
    integer ::       NB_END, NB_END2, i1, i2, ri, iy, ix
    real(kind=DP), allocatable, dimension(:,:) :: bpr_tw1_BLAS
    real(kind=DP), allocatable, dimension(:,:,:) :: p1Sp2_NB
    integer, save :: ibsize_print = OFF
#endif

#ifdef TRANSPOSE_WITHOUT_REARRANGEMENT
    integer :: ierror
#endif
!x!$#ifndef SX
    integer :: ibl1,ibl2,ibsize,ncache
!!$    integer  :: id_sname = -1, id_sname2 = -1, id_sname3 = -1, id_sname4 = -1, id_sname1 = -1
    integer  :: id_sname2 = -1, id_sname3 = -1, id_sname4 = -1, id_sname1 = -1

!!$    if (sw_timing_2ndlevel == ON) call tstatc0_begin('mgs_phi2wf_each_k_G ',id_sname)
    ncache = (m_CtrlP_cachesize()*1024)*3/4
!x!$#endif

    if (sw_timing_2ndlevel == ON) call tstatc0_begin('mgs_phi2wf_each_k_G(1) ',id_sname1)
#ifdef TRANSPOSE_WITHOUT_REARRANGEMENT
    if(ipri>=1) write(nfout,'("A CPP definition of " &
         & ,"TRANSPOSE_WITHOUT_REARRANGEMENT can not be set for MDDAVIDSON")')
    ierror = CPP_DEFINE_ERROR
#ifdef DEBUG_ERRORS
    call phase_error_wo_filename(ierror,nfout,line=__LINE__,modulefile=__FILE__)
#else
    call phase_error_wo_filename(ierror,nfout)
#endif
#endif

#ifdef MGS_DGEMM
    if(nblocksize_mgs_is_given) then
       NB = nblocksize_mgs
    else
       NB = nblocksize_mgs_default
    end if
#endif

    call m_ESortho_mgs_alloc(ik)
    kimg_t_wk = 2
    if(k_symmetry(ik) == GAMMA) kimg_t_wk = 1
#ifdef MGS_DGEMM
    allocate(p1Sp2_NB(NB,NB,kimg_t_wk))
    allocate(psi_t_dia(np_g1k_x,NB,kimg))
#endif

    if(mod_pot == VANDERBILT_TYPE) then
       if(kimg_t_wk == 1) then
          call m_ES_F_transpose_r(.true.,ista_k,iend_k,ik,fsr_l,bpr_t)              ! fsr_l -> bpr_t
          allocate(phifr_t(np_fs_x,neg))
          call m_ES_F_transpose_r(.true.,ik,    ik,    ik,bsdr_l,phifr_t)
#ifdef MGS_DGEMM
          allocate(bpr_t_dia(np_fs_x,NB))
! === Debug by Intel "-check all" option! by T.Kato 2013/07/02 =================
          allocate(bpi_t_dia(1,1))
! ==============================================================================
#endif
       else
          call m_ES_F_transpose_r(.true.,ista_k,iend_k,ik,fsr_l,bpr_t,fsi_l,bpi_t)  ! fs[ri]_l -> bp[ri]_t
          allocate(phifr_t(np_fs_x,neg),phifi_t(np_fs_x,neg))
          call m_ES_F_transpose_r(.true.,ik,    ik,    ik,bsdr_l,phifr_t,bsdi_l,phifi_t)
#ifdef MGS_DGEMM
          allocate(bpr_t_dia(np_fs_x,NB))
          allocate(bpi_t_dia(np_fs_x,NB))
#endif
       end if
    end if

    allocate(phi_t(np_g1k_x,neg,kimg))

    call m_ES_W_transpose_r(.true.,ista_k,iend_k,ik,zaj_l,psi_t) ! zaj_l(ig,ie,ik,ri) -> psi_t(ig,ie,ri) (transposed)
    call m_ES_W_transpose_r(.true.,ik,ik,ik,phi_l,phi_t) ! phi_l(ig,ie,ik,ri) -> phi_t(ig,ie,ri) (transposed)

    if (sw_timing_2ndlevel == ON) call tstatc0_end(id_sname1)
#ifdef MGS_DGEMM
    do i = 1, neg,NB
       NB_END = i + NB -1
       if( NB_END > neg ) NB_END = neg
!diagonal
    if (sw_timing_2ndlevel == ON) call tstatc0_begin('mgs_phi2wf_each_k_G(2) ',id_sname2)
    do i1 = i, NB_END
       if(mode == ORTHONORMALIZATION .or. mode == NORMALIZATION) then
          call WSW_t_g(ik,i1,mod_pot,fr,phi_t,np_g1k_x,neg,kimg,phifr_t,phifi_t)
          if(dabs(fr-1.d0) > DELTA) &
               & call normalize_bp_and_psi_t_g(ik,i1,fr,mod_pot &
               &   , phi_t,np_g1k_x,neg,kimg,phifr_t,phifi_t)
       end if
       if(mod_pot == VANDERBILT_TYPE) call cp_bpr_bsdr2bprtw(i1,neg)  ! phif[ri]_t, bp[ri]_t -> bpr_tw1, bpr_tw2
       if(i1 == neg) cycle
       call cp_psi2psii_g(ik,i1) ! psi_t(:,i,:) -> psi_ir,psi_ii
       if(mod_pot == VANDERBILT_TYPE) &
            & call cp_bpr2bpi_g(kimg_t_wk,i1,bpr_t,bpi_t) ! -> bp_ir, bp_ii
       call W1SW2_t_r_g(ik,i1,NB_END,mod_pot,phi_t,np_g1k_x,neg,kimg) ! -> p1Sp2
       call modify_bp_and_psi_t_r_g(ik,i1,NB_END,mod_pot &
            & ,phi_t,np_g1k_x,neg,kimg,phifr_t,phifi_t) ! psi_t, bpr_t, pbi_t, p1Sp2 -> psi_t, bpr_t, bpi_t
    end do   ! i1-loop
    if (sw_timing_2ndlevel == ON) call tstatc0_end(id_sname2)
!lower
    if (sw_timing_2ndlevel == ON) call tstatc0_begin('mgs_phi2wf_each_k_G(3) ',id_sname3)

    if(mod_pot == VANDERBILT_TYPE) then
       call cp_psi_bpri2dias_g(ik,i,kimg_t_wk,NB,psi_t_dia,np_g1k_x,NB,kimg,np_fs_x,bpr_t_dia,bpi_t_dia)
    else
       call cp_psi_bpri2dias_g(ik,i,kimg_t_wk,NB,psi_t_dia,np_g1k_x,NB,kimg,np_fs_x)
    end if

    do i2 = i+NB, neg,NB
       if(NB_END == neg) cycle
       NB_END2 = i2 + NB -1
       if( NB_END2 > neg ) NB_END2 = neg
! === Debug by Intel "-check all" option! by T.Kato 2013/07/02 =================
!      call W1SW2_t_r_block_g(ik,i,i2,p1Sp2_NB,NB_END,NB_END2,kimg_t_wk &
!           &    , mod_pot,phi_t,psi_t_dia,np_g1k_x,neg,kimg,bpr_tw1(1,i),bpi_tw1(1,i)) ! -> p1Sp2
       if(allocated(bpr_tw1) .and. allocated(bpi_tw1)) then
          call W1SW2_t_r_block_g(ik,i,i2,p1Sp2_NB,NB_END,NB_END2,kimg_t_wk &
               &    , mod_pot,phi_t,psi_t_dia,np_g1k_x,neg,kimg,bpr_tw1(1,i),bpi_tw1(1,i)) ! -> p1Sp2
       else if(allocated(bpr_tw1)) then
          call W1SW2_t_r_block_g(ik,i,i2,p1Sp2_NB,NB_END,NB_END2,kimg_t_wk &
               &    , mod_pot,phi_t,psi_t_dia,np_g1k_x,neg,kimg,bpr_tw1(1,i)) ! -> p1Sp2
       else
          call W1SW2_t_r_block_g(ik,i,i2,p1Sp2_NB,NB_END,NB_END2,kimg_t_wk &
               &    , mod_pot,phi_t,psi_t_dia,np_g1k_x,neg,kimg) ! -> p1Sp2
       end if
! ==============================================================================
       if(nrank_e > 1) call mpi_allreduce(MPI_IN_PLACE, p1Sp2_NB, NB*NB*kimg_t_wk &
            &              ,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
       call modify_bp_and_psi_t_r_blk_g(ik,i,i2,p1Sp2_NB,NB_END,NB_END2,kimg_t_wk &
            &    , mod_pot,phi_t,psi_t_dia,np_g1k_x,neg,kimg,phifr_t,phifi_t,bpr_t_dia,bpi_t_dia)

    end do   ! i2-loop
    if (sw_timing_2ndlevel == ON) call tstatc0_end(id_sname3)
    end do   ! i-loop
#else
    do i = 1, neg
       if(mode == ORTHONORMALIZATION .or. mode == NORMALIZATION) then
          call WSW_t_g(ik,i,mod_pot,fr,phi_t,np_g1k_x,neg,kimg,phifr_t,phifi_t)
          if(dabs(fr-1.d0) > DELTA)  &
               & call normalize_bp_and_psi_t_g(ik,i,fr,mod_pot &
               &   , phi_t,np_g1k_x,neg,kimg,phifr_t,phifi_t)
       end if
!x!$       if(mod_pot == VANDERBILT_TYPE) call cp_bpr_bsdr2bprtw(i,neg)
!x!$       if(i == neg) cycle
          call cp_psi2psii_g(ik,i) ! psi_t(:,i,:) -> psi_ir,psi_ii
       if(mod_pot == VANDERBILT_TYPE) &
            & call cp_bpr2bpi_g(kimg_t_wk,i,bpr_t,bpi_t) ! -> bp_ir, bp_ii
       call W1SW2_t_r_g(ik,i,neg,mod_pot,phi_t,np_g1k_x,neg,kimg & ! ->p1Sp2
            & ,phifr_t,phifi_t)
       call modify_bp_and_psi_t_r_g(ik,i,neg,mod_pot &
            & ,phi_t,np_g1k_x,neg,kimg,phifr_t,phifi_t)
       ! phi_t, psi_t, phifr_t,phifi_t,bpr_t, pbi_t, p1Sp2 -> phi_t, phifr_t,phifi_t

    end do
#endif

    if (sw_timing_2ndlevel == ON) call tstatc0_begin('mgs_phi2wf_each_k_G(4) ',id_sname4)
    call m_ES_W_transpose_back_r(.true.,ik,ik,ik,phi_l,phi_t)       ! phi_t -> phi_l
    if(mod_pot == VANDERBILT_TYPE) then
       if(kimg_t_wk == 1) then
          call m_ES_F_transpose_back_r(.true.,ik,ik,ik,bsdr_l,phifr_t)
#ifdef MGS_DGEMM
          deallocate(bpr_t_dia)
! === Debug by Intel "-check all" option! by T.Kato 2013/07/02 =================
          deallocate(bpi_t_dia)
! ==============================================================================
#endif
       else
          call m_ES_F_transpose_back_r(.true.,ik,ik,ik,bsdr_l,phifr_t,bsdi_l,phifi_t)
#ifdef MGS_DGEMM
          deallocate(bpi_t_dia,bpr_t_dia)
#endif
       end if
    end if

    deallocate(phi_t)

#ifdef MGS_DGEMM
    deallocate(p1Sp2_NB)
    deallocate(psi_t_dia)
#endif

    if(mod_pot == VANDERBILT_TYPE) then
       if(kimg_t_wk == 1) then
          deallocate(phifr_t)
       else
          deallocate(phifr_t,phifi_t)
       end if
    end if

    call m_ESortho_mgs_dealloc
    if (sw_timing_2ndlevel == ON) call tstatc0_end(id_sname4)
!!$    if (sw_timing_2ndlevel == ON) call tstatc0_end(id_sname)

#ifdef MGS_DGEMM
  contains
    subroutine cp_bpr_bsdr2bprtw(i1,i2)
      integer, intent(in) :: i1,i2
      integer :: i, ia, p, q, i_last
#ifdef MGS_DGEMM_DEBUG
      i_last = i1
#else
      i_last = i2
#endif
      if(k_symmetry(ik) == GAMMA) then
         do i = i1, i_last
            do ia = 1, nac_p
               p = nlmta1_p(ia);     q = nlmta2_p(ia)
               bpr_tw1(ia,i) = bpr_t(p,i)
               bpr_tw2(ia,i) = phifr_t(q,i)
            end do
         end do
      else
         do i = i1, i_last
            do ia = 1, nac_p
               p = nlmta1_p(ia);     q = nlmta2_p(ia)
               bpr_tw1(ia,i) = bpr_t(p,i)
               bpr_tw2(ia,i) = phifr_t(q,i)
               bpi_tw1(ia,i) = bpi_t(p,i)
               bpi_tw2(ia,i) = phifi_t(q,i)
!x!$               bpr_tw1(ia,i) = phifr_t(p,i)
!x!$               bpr_tw2(ia,i) = bpr_t(q,i)
!x!$               bpi_tw1(ia,i) = phifi_t(p,i)
!x!$               bpi_tw2(ia,i) = bpi_t(q,i)
            end do
         end do
      end if
    end subroutine cp_bpr_bsdr2bprtw
#endif
  end subroutine mgs_phi2wf_each_k_G

end module m_ES_ortho
