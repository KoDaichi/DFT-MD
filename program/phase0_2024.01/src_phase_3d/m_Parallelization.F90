#ifdef HIUX
*option MP(P(0))
#endif
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 635 $)
!
!  MODULE: m_Parallelization
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!  FURTHER MODIFICATION: T. Yamasaki, August/23/2006
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
! ***************************************************************
!
! This module "m_Parallelization" is for mpi-parallelization and was
! coded by T. Yamasaki (JRCAT-ATP, Fujitsu Laboratories Ltd.) in 1999.
!
#ifdef __TIMER__
#   define __TIMER_START(a) call timer_init() ; call timer_sta(a)
#   define __TIMER_STOP     call timer_fin()
#else
#   define __TIMER_START(a)
#   define __TIMER_STOP
#endif
#ifdef __TIMER_SUB__
#   define __TIMER_SUB_START(a)  call timer_sta(a)
#   define __TIMER_SUB_STOP(a)   call timer_end(a)
#else
#   define __TIMER_SUB_START(a)
#   define __TIMER_SUB_STOP(a)
#endif
#ifdef __TIMER_INICOMM__
#   define __TIMER_INICOMM_START_w_BARRIER(str,a)   call timer_barrier(str) ; call timer_sta(a)
#   define __TIMER_INICOMM_STOP(a)         call timer_end(a)
#else
#   define __TIMER_INICOMM_START_w_BARRIER(str,a)
#   define __TIMER_INICOMM_STOP(a)
#endif

module m_Parallelization
!   (m_Parallel)
! $Id: m_Parallelization.F90 635 2021-02-26 07:16:10Z jkoga $
#ifdef MPI_FFTW
  use, intrinsic :: iso_c_binding
#endif
  use m_Const_Parameters, only       : ON, OFF, tag_npes_etc
  use m_ErrorMessages
! === For nrc decomposion. by takto 2012/12/07 =================================
  use m_Const_Parameters, only       : paw_nrc_min_elements
! ==============================================================================
  use mpi

  implicit none
  integer :: npes, mype, nrank_e, nrank_k, nrank_g, nrank_s
  integer, allocatable, dimension(:)   :: lrank   ! lrank(nbs) rank_e number of nbs
  integer, allocatable, dimension(:)   :: nbsn    ! nbsn(nbs) local block number of nbs
  integer, allocatable, dimension(:)   :: nbsn_sta, nbsn_end    ! nbsn_sta(nbsn) loop of local block
  integer, allocatable, dimension(:)   :: nbs_sta, nbs_end    ! nbs_sta(nbs) loop of global block
  integer, allocatable, dimension(:)   :: neg_g    ! local index
  integer, allocatable, dimension(:)   :: neg_g_all  ! order of block-cyclic for all neg
  integer, allocatable, dimension(:)   :: ball_buff,ball_addr
  integer :: nbs      ! block number
  integer :: nbs_num      ! number of block
  integer :: nbsn_num      ! number of local block
  integer, allocatable, dimension(:) :: map_fft_z, map_fft_x, map_fft_y
  integer, allocatable, dimension(:) :: nel_fft_z, nel_fft_x, nel_fft_y
  integer, allocatable, dimension(:) ::  mp_fft_z,  mp_fft_x,  mp_fft_y
  integer                            ::  np_fft_z,  np_fft_x,  np_fft_y
  integer, dimension(2,3)            :: xyz_fft_z, xyz_fft_x, xyz_fft_y
  integer                            :: fft_X_z_dim, fft_X_y_dim
  integer                            :: fft_X_x_dim
  integer                            :: fft_X_x_nel, fft_X_y_nel, fft_X_z_nel
  integer                            :: fft_Y_x_dim, fft_Y_z_dim
  integer                            :: fft_Z_x_dim, fft_Z_y_dim
  integer, allocatable, dimension(:) :: nis_fft_X_z, nis_fft_X_y, nie_fft_X_z, nie_fft_X_y
  integer, allocatable, dimension(:) :: nis_fft_Y_x, nis_fft_Y_z, nie_fft_Y_x, nie_fft_Y_z
  integer, allocatable, dimension(:) :: nis_fft_Z_x, nis_fft_Z_y, nie_fft_Z_x, nie_fft_Z_y
  integer, allocatable, dimension(:) :: nis_fft_X_x, nie_fft_X_x
  integer, allocatable, dimension(:) :: map_fftcd_z, map_fftcd_x, map_fftcd_y
  integer, allocatable, dimension(:) :: nel_fftcd_z, nel_fftcd_x, nel_fftcd_y
  integer, allocatable, dimension(:) ::  mp_fftcd_z,  mp_fftcd_x,  mp_fftcd_y
  integer                            ::  np_fftcd_z,  np_fftcd_x,  np_fftcd_y
  integer, dimension(2,3)            :: xyz_fftcd_z, xyz_fftcd_x, xyz_fftcd_y
  integer                            :: fftcd_X_z_dim, fftcd_X_y_dim
  integer                            :: fftcd_X_x_dim
  integer                            :: fftcd_X_x_nel, fftcd_X_y_nel, fftcd_X_z_nel
  integer                            :: fftcd_Y_x_dim, fftcd_Y_z_dim
  integer                            :: fftcd_Z_x_dim, fftcd_Z_y_dim
  integer, allocatable, dimension(:) :: nis_fftcd_X_z, nis_fftcd_X_y, nie_fftcd_X_z, nie_fftcd_X_y
  integer, allocatable, dimension(:) :: nis_fftcd_Y_x, nis_fftcd_Y_z, nie_fftcd_Y_x, nie_fftcd_Y_z
  integer, allocatable, dimension(:) :: nis_fftcd_Z_x, nis_fftcd_Z_y, nie_fftcd_Z_x, nie_fftcd_Z_y
  integer, allocatable, dimension(:) :: nis_fftcd_X_x, nie_fftcd_X_x
  integer                            :: ista_kngp_B, iend_kngp_B, np_kngp_B, mp_kngp_B
  integer, allocatable, dimension(:) :: is_kngp_B, ie_kngp_B, nel_kngp_B

  integer, allocatable, dimension(:,:)   :: wf_fft_scnt,    wf_fft_rcnt
  integer, allocatable, dimension(:,:,:) :: wf_fft_send,    wf_fft_recv
  integer, allocatable, dimension(:,:)   :: wf_fft_dist,    wf_fft_index
  integer, allocatable, dimension(:)   :: wf_fft_maxsend, wf_fft_maxrecv
  integer, allocatable, dimension(:,:)   :: fft_wf_scnt,    fft_wf_rcnt
  integer, allocatable, dimension(:,:,:) :: fft_wf_send,    fft_wf_recv
  integer, allocatable, dimension(:,:)   :: fft_wf_dist,    fft_wf_index
  integer, allocatable, dimension(:)   :: fft_wf_maxsend, fft_wf_maxrecv
  integer, allocatable, dimension(:)   :: fft_chgq_scnt,    fft_chgq_rcnt
  integer, allocatable, dimension(:,:) :: fft_chgq_send,    fft_chgq_recv
  integer, allocatable, dimension(:)   :: fft_chgq_dist,    fft_chgq_index
  integer                              :: fft_chgq_maxsend, fft_chgq_maxrecv
#ifdef MPI_FFTW
  integer, allocatable, dimension(:)   :: fft_chgq_scnt_mpifftw,    fft_chgq_rcnt_mpifftw
  integer, allocatable, dimension(:,:) :: fft_chgq_send_mpifftw,    fft_chgq_recv_mpifftw
  integer, allocatable, dimension(:)   :: fft_chgq_dist_mpifftw,    fft_chgq_index_mpifftw
  integer                              :: fft_chgq_maxsend_mpifftw, fft_chgq_maxrecv_mpifftw
#endif
  integer, allocatable, dimension(:)   :: chgq_fftcd_scnt,    chgq_fftcd_rcnt
  integer, allocatable, dimension(:,:) :: chgq_fftcd_send,    chgq_fftcd_recv
  integer, allocatable, dimension(:)   :: chgq_fftcd_index,   chgq_fftcd_dist
  integer                              :: chgq_fftcd_maxrecv, chgq_fftcd_maxsend
  integer, allocatable, dimension(:)   :: fftcd_chgq_scnt,    fftcd_chgq_rcnt
  integer, allocatable, dimension(:,:) :: fftcd_chgq_send,    fftcd_chgq_recv
  integer, allocatable, dimension(:)   :: fftcd_chgq_dist,    fftcd_chgq_index
  integer                              :: fftcd_chgq_maxsend, fftcd_chgq_maxrecv
  integer                              :: mpi_fft_xz_world,   myrank_fft_xz,   nrank_fft_xz
  integer                              :: mpi_fft_zy_world,   myrank_fft_zy,   nrank_fft_zy
  integer                              :: mpi_fftcd_xz_world, myrank_fftcd_xz, nrank_fftcd_xz
  integer                              :: mpi_fftcd_zy_world, myrank_fftcd_zy, nrank_fftcd_zy
  integer                              :: mpi_fft_xy_world,   myrank_fft_xy,   nrank_fft_xy
  integer                              :: mpi_fft_yz_world,   myrank_fft_yz,   nrank_fft_yz
  integer                              :: mpi_fftcd_xy_world,   myrank_fftcd_xy,   nrank_fftcd_xy
  integer                              :: mpi_fftcd_yz_world,   myrank_fftcd_yz,   nrank_fftcd_yz
  integer, allocatable, dimension(:)   :: igfp_full
  logical                              :: init_zaj_para

#ifdef MPI_FFTW
  integer, allocatable, dimension(:,:)   :: wf_fft_scnt_mfftw,    wf_fft_rcnt_mfftw
  integer, allocatable, dimension(:,:,:) :: wf_fft_send_mfftw,    wf_fft_recv_mfftw
  integer, allocatable, dimension(:,:)   :: wf_fft_dist_mfftw,    wf_fft_index_mfftw
  integer, allocatable, dimension(:)     :: wf_fft_maxsend_mfftw, wf_fft_maxrecv_mfftw
  integer, allocatable, dimension(:,:)   :: fft_wf_scnt_mfftw,    fft_wf_rcnt_mfftw
  integer, allocatable, dimension(:,:,:) :: fft_wf_send_mfftw,    fft_wf_recv_mfftw
  integer, allocatable, dimension(:,:)   :: fft_wf_dist_mfftw,    fft_wf_index_mfftw
  integer, allocatable, dimension(:)     :: fft_wf_maxsend_mfftw, fft_wf_maxrecv_mfftw
#endif

  integer :: nrank_g1
  integer, allocatable, dimension(:) :: map_e, map_z              ! d(neg)  --> expl.1 and 5
  integer                            :: myrank_e,ista_e,iend_e,istep_e,np_e,mp_e ! --> expl.2 to 4

#ifdef NEC_TUNE_SOFT
  integer, allocatable, dimension(:) :: ista_e_smp, iend_e_smp
  integer                            :: itask
#elif NEC_TUNE_FFT
  integer                            :: itask
#endif
  integer, allocatable, dimension(:) :: nis_e, nie_e, nel_e, idisp_e   ! d(nrank_e)
  integer, allocatable, dimension(:) :: map_k                     ! d(kv3)
  integer, allocatable, dimension(:,:):: map_ek                   ! d(neg,kv3)
  integer, allocatable, dimension(:) :: map_s                   ! d(neg,kv3)
  integer                            :: myrank_k,ista_k,iend_k
!  integer, allocatable, dimension(:,:,:) :: map_rank_gek  ! d(nrank_g,nrank_e,nrank_k)
  integer, allocatable, dimension(:,:,:,:) :: map_rank_gek  ! d(nrank_g,nrank_e,nrank_k,nrak_s)
  integer, allocatable, dimension(:) :: nis_k,  nie_k,  nel_k
  integer, allocatable, dimension(:,:,:):: map_ekg                  ! d(neg,kg1,kv3)
  integer, allocatable, dimension(:) :: nis_eg, nie_eg, nel_eg, neg_gg, neg_gg_all ! d(nrank_g)
!#ifdef _USE_SCALAPACK_
  integer :: scl_row, scl_col, my_row, my_col
  integer :: scl_comm_max, scl_comm_rank, scl_comm_max_r, scl_comm_rank_r
  integer , allocatable, dimension(:)   :: scl_comm_rno, scl_comm_rno_r
  integer , allocatable, dimension(:)   :: neg_row, neg_col
  integer , allocatable, dimension(:)   :: nis_row, nis_col
  integer , allocatable, dimension(:)   :: nie_row, nie_col
  integer , allocatable, dimension(:)   :: nel_row, nel_col
  integer , allocatable, dimension(:,:) :: nrm_rank, scl_rank
!#endif
! add for matrixdiagon >
  integer :: scl_row_md, scl_col_md
  integer :: scl_md_comm_max, scl_md_comm_rank, scl_md_comm_max_r, scl_md_comm_rank_r
  integer , allocatable, dimension(:)   :: scl_md2_comm_rno, scl_md2_comm_rno_r
  integer :: scl_md2_comm_max, scl_md2_comm_rank, scl_md2_comm_max_r, scl_md2_comm_rank_r
  integer , allocatable, dimension(:)   :: scl_md_comm_rno, scl_md_comm_rno_r
  integer , allocatable, dimension(:)   :: nmatsz_row, nmatsz_col
  integer , allocatable, dimension(:)   :: nis_row_md, nis_col_md
  integer , allocatable, dimension(:)   :: nie_row_md, nie_col_md
  integer , allocatable, dimension(:)   :: nel_row_md, nel_col_md
! < add for matrixdiagon


!!$#ifdef TRANSPOSE
  integer                            :: ista_g1, iend_g1, np_g1, mp_g1
  integer, allocatable, dimension(:) :: nis_g1, nie_g1, nel_g1      !d(0:nrank_g1-1)
  integer, allocatable, dimension(:) :: ista_g1k, iend_g1k, np_g1k, mp_g1k, np_g1k_prev !d(kv3)
  integer, allocatable, dimension(:) :: ista_g1k_prev, iend_g1k_prev
  integer, allocatable, dimension(:,:):: nis_g1k, nie_g1k, nel_g1k  !d(0:nrank_g1-1,kv3)

  integer                            :: ista_fs, iend_fs, np_fs, mp_fs
  integer, allocatable, dimension(:) :: nis_fs, nie_fs, nel_fs      !d(0:nrank_g-1)
  integer                            :: ista_fs_atm, iend_fs_atm, np_fs_atm, mp_fs_atm
  integer, allocatable, dimension(:) :: nis_fs_atm, nie_fs_atm, nel_fs_atm !d(0:nrank_e-1)
  integer                            :: myrank_g
  integer                            :: ista_lmta,iend_lmta
  integer, allocatable, dimension(:) :: lmta_atm
! add for matrixdiagon >
  integer, allocatable, dimension(:) :: ista_G_g1k, iend_G_g1k, np_G_g1k, mp_G_g1k
  integer, allocatable, dimension(:,:):: nis_G_g1k, nie_G_g1k, nel_G_g1k, map_G_g1k
  integer, allocatable, dimension(:) :: ista_B_g1k, iend_B_g1k, np_B_g1k, mp_B_g1k
  integer, allocatable, dimension(:,:):: nis_B_g1k, nie_B_g1k, nel_B_g1k, map_B_g1k
! < add for matrixdiagon

!!$#endif
  integer, allocatable, dimension(:) :: mpi_k_world  ! kd(0:nrank_k-1)
  integer, allocatable, dimension(:) :: mpi_e_world  ! kd(0:nrank_k-1)
  integer                            :: mpi_kg_world
  integer                            :: mpi_ke_world
  integer                            :: mpi_kr_world

  integer                            :: mpi_g_world  ! kd(0:nrank_k-1)
  integer                            :: myrank_ke
  integer                            :: nrank_ke

  integer                            :: mpi_chg_world
  integer                            :: myrank_chg
  integer                            :: nrank_chg

! === DEBUG by tkato 2012/04/04 ================================================
  integer                            :: mpi_ge_world
! ==============================================================================
  integer                            :: ista_kngp, iend_kngp, np_kngp, mp_kngp
  integer                            :: ista_kngp_prev, iend_kngp_prev, np_kngp_prev
  integer, allocatable, dimension(:) :: is_kngp, ie_kngp, nel_kngp

  integer                            :: ista_kngp_gw, iend_kngp_gw, np_kngp_gw, mp_kngp_gw
  integer, allocatable, dimension(:) :: is_kngp_gw, ie_kngp_gw, nel_kngp_gw
  integer                            :: ista_kngp_exx, iend_kngp_exx, np_kngp_exx, mp_kngp_exx

  integer, allocatable, dimension(:) :: is_kngp_exx, ie_kngp_exx, nel_kngp_exx
! natm for fxyzew_l
  integer                            :: ista_atm, iend_atm, np_atm, mp_atm
  integer, allocatable, dimension(:) :: is_atm, ie_atm, nel_atm

  integer                            :: ista_atm_ke, iend_atm_ke, np_atm_ke, mp_atm_ke
  integer, allocatable, dimension(:) :: is_atm_ke, ie_atm_ke, nel_atm_ke,map_atm_ke

  integer                            :: ista_atm_f, iend_atm_f, np_atm_f, mp_atm_f
  integer, allocatable, dimension(:) :: is_atm_f, ie_atm_f, nel_atm_f

  integer                            :: ista_rspace_aug, iend_rspace_aug, np_rspace_aug, mp_rspace_aug
  integer, allocatable, dimension(:) :: ista_rspace_aug_atm, iend_rspace_aug_atm
  integer, allocatable, dimension(:) :: is_rspace_aug, ie_rspace_aug, nel_rspace_aug

  integer                            :: ista_nn,iend_nn,np_nn,mp_nn
  integer, allocatable, dimension(:) :: is_nn, ie_nn, nel_nn

  integer                            :: ista_nq, iend_nq, np_nq, mp_nq
  integer, allocatable, dimension(:) :: is_nq, ie_nq, nel_nq, map_nq, map_z_nq

!  integer                            :: ista_nval,iend_nval,np_nval,mp_nval
!  integer, allocatable, dimension(:) :: is_nval, ie_nval, nel_nval

  integer                            :: ista_atm_B, iend_atm_B, np_atm_B, mp_atm_B
  integer, allocatable, dimension(:) :: is_atm_B, ie_atm_B, nel_atm_B
  integer, allocatable, dimension(:) :: mem_atm_B ! (1:natm)
! natm2 for s_ew in m_Ionic_System
  integer                            :: ista_atm2, iend_atm2, np_atm2, mp_atm2
  integer, allocatable, dimension(:) :: is_atm2, ie_atm2, nel_atm2

! BROYDEN or DFP MIXING METHOD
  integer                            :: ista_kgpm,iend_kgpm, np_kgpm, mp_kgpm
  integer, allocatable, dimension(:) :: is_kgpm,  ie_kgpm, nel_kgpm
  integer                            :: ista_urec_hsr, iend_urec_hsr
  logical, save ::                      ista_and_iend_urec_hsr_set = .false.

! SPIN
  integer                            :: ista_spin, iend_spin, np_spin, mp_spin
  integer, allocatable, dimension(:) :: nis_spin, nie_spin, nel_spin
  integer                            :: myrank_spin
  integer                            :: mpi_spin_group, mpi_keg_world
  integer                            :: mpi_skg_world, mpi_sge_world, mpi_ske_world

! FFT BOX
  integer                            :: npes_cdfft, nrank_ggacmp, nrest_cdfft
!!$  integer                            :: myrank_cdfft, myrank_ggacmp, max_ggacmp
  integer                            :: myrank_cdfft, myrank_ggacmp
  integer, allocatable, dimension(:) :: map_ggacmp ! d(3)
  integer, allocatable, dimension(:) :: mpi_ggacmp_cross_world ! d(0:npes_cdfft-1)
  integer, allocatable, dimension(:) :: mpi_cdfft_world        ! d(0:nrank_ggacmp-1)
  integer, allocatable, dimension(:) :: map_pe2ggacmp, map_pe2cdfft    ! d(0:npes-1)

  integer                            :: ista_sfftp, iend_sfftp, np_sfftp, mp_sfftp
  integer                            :: ista_sfftph,iend_sfftph
  integer, allocatable, dimension(:) :: nis_sfftp,  nie_sfftp,  nel_sfftp, idisp_sfftp
  integer, allocatable, dimension(:) :: nis_sfftph, nie_sfftph, nel_sfftph
  integer                            :: ista_fftp, iend_fftp, np_fftp, mp_fftp
  integer                            :: ista_fftph,iend_fftph
  integer, allocatable, dimension(:) :: nis_fftp,  nie_fftp,  nel_fftp,  idisp_fftp
  integer, allocatable, dimension(:) :: nis_fftph, nie_fftph, nel_fftph, idisp_fftph


! NonLocalPotential snl
  integer                            :: ista_snl, iend_snl

! nbmx
  integer  :: ng_nbmx, myrank_nbmx, nrank_nbmx, np_nbmx, mp_nbmx,ista_nbmx, iend_nbmx
  integer  :: nbmx_ext
  integer, allocatable, dimension(:) :: nis_nbmx,nie_nbmx,nel_nbmx, idisp_nbmx
  integer, allocatable, dimension(:) :: mpi_nbmx_world  ! kd(0:nrank_nbmx-1)
  integer  :: ng_nbmx_k, myrank_nbmx_k, nrank_nbmx_k, np_nbmx_k, mp_nbmx_k &
       & ,ista_nbmx_k, iend_nbmx_k
  integer, allocatable, dimension(:) :: nis_nbmx_k, nie_nbmx_k, nel_nbmx_k, idisp_nbmx_k
! kg1
  integer, allocatable, dimension(:) :: mpi_nbmx_world_k  ! kd(0:nrank_nbmx_k-1)
  integer  :: np_kg1_k, mp_kg1_k, ista_kg1_k, iend_kg1_k
  integer  :: kg1_ext
  integer, allocatable, dimension(:) :: nis_kg1_k, nie_kg1_k, nel_kg1_k, idisp_kg1_k

  integer                            :: ierr
  character(len=50)                  :: workdir
  integer                            :: sw_wdir = 0

  integer                            :: MPI_CommGroup
  integer                            :: color,key

  logical,save                       :: conf_para = .false.
  integer,save                       :: mype_conf = 0
  integer                            :: nrank_conf
!!$#ifdef PAW3D
  integer :: nrank_natm, nrank_nrc
  integer :: myrank_natm, myrank_nrc
  integer :: mpi_natm_world, mpi_nrc_world
  integer, allocatable, dimension(:) :: is_natm, ie_natm, nel_natm
  integer, allocatable, dimension(:) :: is_nrc, ie_nrc, nel_nrc
  integer :: ista_natm, iend_natm, ne_natm
  integer :: ista_nrc,  iend_nrc,  ne_nrc
! === For nrc decomposion. by takto 2012/12/07 =================================
  logical :: paw_no_work_to_do     = .false., &
             paw_last_rank_on_nrc  = .false.
! ==============================================================================
!!$#endif

  integer :: nrank_nvale,myrank_nvale
  integer :: ista_nvale,iend_nvale,np_nvale,mp_nvale
  integer, allocatable, dimension(:) :: nis_nvale,nie_nvale,nel_nvale
  integer, allocatable, dimension(:) :: map_nvale,map_z_nvale
  logical,private :: mpi_nvale_enabled = .false.

  integer :: nrank_kv3_ek,myrank_kv3_ek
  integer :: ista_kv3_ek,iend_kv3_ek,np_kv3_ek,mp_kv3_ek
  integer, allocatable, dimension(:) :: nis_kv3_ek,nie_kv3_ek,nel_kv3_ek
  integer, allocatable, dimension(:) :: map_kv3_ek,map_z_kv3_ek

  integer :: nrank_nval,myrank_nval
  integer :: ista_nval,iend_nval,np_nval,mp_nval
  integer, allocatable, dimension(:) :: nis_nval,nie_nval,nel_nval
  integer, allocatable, dimension(:) :: map_nval,map_z_nval
  logical,private :: mpi_nval_enabled = .false.

! for nfft
  integer                            :: ista_ffth, iend_ffth, np_ffth, mp_ffth
  integer, allocatable, dimension(:) :: is_ffth, ie_ffth, nel_ffth

! whether parallelization was configured from the command line or not
  logical :: read_from_args = .false.

#ifdef KMATH_FFT3D
  integer :: kmath3d_handle_wf_direct, kmath3d_handle_wf_inverse
  integer, dimension(3) :: nproc_fft3d_inverse, nproc_fft3d_direct,kmath3d_box_size
#endif
! 1. map_e(neg)
! rank# = map_e(ieg) : rank# is the rank number which operates ieg-th data
!  -- examples --
!  (case 1) cyclic partitioning case (assuming nrank_e == 3 and neg == 8)
!   ieg :   1   2   3   4   5   6   7   8
! rank# :   0   1   2   0   1   2   0   1
!
!  (case 2) block  partitioning case (also assuming nrank_e == 3 and neg == 8)
!   ieg :   1   2   3   4   5   6   7   8
! rank# :   0   0   0   1   1   1   2   2
!
! 2. myrank_e
!  This lets each node know its rank number according to an axis-e.
!  -- an example --
!  Assuming nrank_e == 4 and nrank_k == 2
!   #rank (or #pe)  0  1  2  3  4  5  6  7
!   myrank_e        0  1  2  3  0  1  2  3
!
! 3. ista_e, iend_e, istep_e, nis_e, nie_e
!  ista_e, iend_e : Bounds of a window which is a range 'e' runs in each nrank_e.
!  nis_e,  nie_e  : Bounds of windows. This shows all windows ranges.
!  istep_e        : = nrank_e (the cyclic partitioning case)
!                   = 1       (the block partitioning case)
!  -- examples --
!  Assuming nrank_e == 3, nrank_k == 2 and neg ==  8
!  (case 1) cyclic partitioning case
!   #rank (or #pe)  0  1  2  3  4  5
!   ista_e          1  2  3  1  2  3
!   iend_e          7  8  6  7  8  6
!   istep_e         3  3  3  3  3  3
!   nis_e(#rank)    1  2  3  1  2  3
!   nie_e(#rank)    7  8  6  7  8  6
!
!  (case 2) block partitioning case
!   #rank (or #pe)  0  1  2  3  4  5
!   ista_e          1  4  7  1  4  7
!   iend_e          3  6  8  3  6  8
!   istep_e         1  1  1  1  1  1
!   nis_e(#rank)    1  4  7  1  4  7
!   nie_e(#rank)    3  6  8  3  6  8
!
! 4. np_e, mp_e, nel_e
!  np_e is the number of elements of e which each rank_e operates
!  mp_e is the maximum of np_e's over all nodes
!  -- examples --
!  Assuming nrank_e == 3, nrank_k == 2 and neg ==  8
!   #rank (or #pe)  0  1  2  3  4  5
!   np_e            3  3  2  3  3  2
!   mp_e            3  3  3  3  3  3
!   nel_e(#rank)    3  3  2  3  3  2
!
! 5. map_z(neg)
! j = map_z(ieg) : ieg-th data is mapped on j-th position on the
!                 map_e(i)-th array. Common on all ranks
!  -- examples --
! Assuming nrank_e = 3 and neg = 8
!  (case 1) cyclic partitioning case
!   ieg    :   1   2   3   4   5   6   7   8
!  rank_e# :   0   1   2   0   1   2   0   1
!    j     :   1   1   1   2   2   2   3   3
!  (case 2) block partitioning case
!   ieg    :   1   2   3   4   5   6   7   8
!  rank_e# :   0   0   0   1   1   1   2   2
!    j     :   1   2   3   1   2   3   1   2
!
! 6. map_k(kv3)
! rank_k# = map_k(ik) : rank_k# is the rank number which operates ik-th data
!  -- an examples --
!  Assuming nrank_k == 2 and kv3 == 9)
!    ik   :   1   2   3   4   5   6   7   8   9
! rank_k# :   0   0   0   0   0   1   1   1   1
!
! 7. map_ek(neg,kv3)
! rank# = map_ek(ie,ik) : rank# is,  in the whole using nodes, the rank number which operates
!                         ie-th and -k-th data
!  -- an example --
!  Assuming nrank_e == 2, nrank_k = 2, npes == 4, neg == 8 and kv3 = 5
!  and taking block partitioning for an axis-'e'.
!   ik   : 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5
!   ie   : 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8
! rank_k#: 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
! rank_e#: 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1
! rank#  : 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3 2 2 2 2 3 3 3 3
!
! 8. myrank_k
!  This lets each node know its rank number according to an axis-k.
!  -- an example --
!  Assuming nrank_e == 4 and nrank_k == 2
!   #rank (or #pe)  0  1  2  3  4  5  6  7
!   myrank_k        0  0  0  0  1  1  1  1
!   myrank_e        0  1  2  3  0  1  2  3  (for reference, assuming block partitioning)
!
! 9. ista_k, iend_k
!  ista_k, iend_k : bounds of a window which is a range each nrank_k node operates
!
!10. ista_g1, iend_g1, nis_g1, nie_g1, nel_g1
!  ista_g1,iend_g1: bounds of a window which is a range 'g' runs in each nrank_e after transpose.
!  -- An example --
!  Assuming nrank_e == 3, nrank_k == 2 and kg1 == 101
!   #rank (or #pe)  0   1   2   3   4   5
!   ista_e          1  35  69   1  35  69
!   iend_e         34  68 101  34  68 101
!
!11. np_g1, mp_g1
!  np_g1 is the number of elements of 'g' which each rank_e operates
!  mp_g1 is the maximum of np_g1's over all nodes
!  -- examples --
!  Assuming nrank_e == 3, nrank_k == 2 and kg1 ==  101
!   #rank (or #pe)  0  1  2  3  4  5
!   np_g1          34 34 33 34 34 33
!   mp_g1          34 34 34 34 34 34
!
!12. mpi_k_world(0:nrank_k-1)
! Sub-group of rank_k
!
!13. ng_nbmx, myrank_nbmx, nrank_nbmx,np_nbmx, mp_nbmx, nel_nbmx(#rank_nbmx),
!    ista_nbmx, and iend_nbmx
!  np_nbmx is the number of elements of G which each rank_nbmx operates
!  mp_nbmx is the maximum of np_nbmx's over all nodes in a mpi_nbmx_world
!  -- example --
!  Assuming nbmx = 152, nrank_nbmx = 3 or 2, npe = 5, and ngnode_nbmx = 2
!  #rank (or #pe)               0   1   2   3   4
!  ng_nbmx                      0   1   0   1   0
!  myrank_nbmx                  0   0   1   1   2
!  nrank_nbmx                   3   2   3   2   3
!  np_nbmx                     51  76  51  76  50
!  mp_nbmx                     51  76  51  76  51
!  nel_nbmx(#rank_nbmx)        51  76  51  76  50
!  ista_nbmx                    1   1  52  77 103
!  iend_nbmx                   51  76 102 152 152
!
!                          ||
!
!  ng_nbmx            (      0       ) (    1    )
!  #rank ( or #node)    0    2    4       1   3
!  nrank_nbmx           3    3    3       2   2
!  myrank_nbmx          0    1    2       0   1
!    np_nbmx           51   51   50      76  76
!    mp_nbmx           51   51   51      76  76
! nel_nbmx(#rank_nbmx) 51   51   50      76  76
!    ista_nbmx          1   52  103       1  77
!    iend_nbmx         51  102  152      76 152
!
!14. mpi_nbmx_world(0:nrank_nbmx-1)
! Sub-group of rank_nbmx
!
#ifdef MPI_FFTW
  include 'fftw3-mpi.f03'
#endif
!  include 'mpif.h'             ! MPI

  integer                            :: istatus(MPI_STATUS_SIZE)
contains

  subroutine m_Parallel_init_comm_world(init_mpi)
!!$    character :: F_OUT*13, name_mype*3
!!$    character :: name_jobstep*3
!!$    integer ::   i, js
!!$    logical ::   existence
    integer, optional, intent(in) :: init_mpi
    integer :: inimpi
    integer :: i,lenc,lencom,ntmp,numdir,np,hispe
    logical :: exist_dirlist
    character(len=50) :: buf
    character(len=50), allocatable :: workdir_array(:)

    integer, allocatable :: color(:),key(:),nproc(:)
    integer :: my_color, my_key
    integer :: maxproc, count, idir, npes1, numdir1
    logical, save :: firstcall = .true.

    workdir = ''

    inimpi=ON
    if(present(init_mpi)) inimpi = init_mpi

    if(inimpi==ON) then

    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world, npes, ierr)
    call mpi_comm_rank(mpi_comm_world, mype, ierr)
                                    __TIMER_START(17)
#ifdef FJ_TIMER
                 call timer_sta(21)
#endif

!! Group creation
    inquire(file='dirlist',exist=exist_dirlist)
    if(.not.exist_dirlist) then

       !if(firstcall) call mpi_comm_dup(mpi_comm_world,MPI_CommGroup,ierr)
      if(.not.firstcall) call mpi_comm_free(MPI_CommGroup, ierr)
      call mpi_comm_dup(mpi_comm_world,MPI_CommGroup,ierr)

       !!write(*,*) 'mpi_comm_world=',mpi_comm_world
       !!write(*,*) 'MPI_CommGroup=',MPI_CommGroup

    else
       !!write(*,*) 'dirlist exists.'
       !! switching working directory
       sw_wdir = 1

       if(mype == 0) then
          open(unit=100,file="dirlist",form='formatted')
          read(100,*) numdir
          !! debug
          !!write(*,*) 'numdir=',numdir
          allocate(workdir_array(numdir))
          ntmp=0
          do i=1,numdir
             buf(1:50) = ' '
             read(100,*,err=500,end=500) buf
             ntmp=ntmp+1
             lenc = len_trim(buf)
             buf(lenc+1:lenc+1) = '/'
             workdir_array(i) = buf
          end do
500       continue
          numdir = ntmp

          allocate(color(npes))
          allocate(key(npes))
          allocate(nproc(npes))

          npes1 = npes
          numdir1 = numdir
          count = 0
          my_color = 0
          do while(npes1 > 0)
             maxproc = npes1/numdir1
             if(npes1-maxproc*numdir1>0) maxproc=maxproc+1
             my_color = my_color + 1
             do i=1,maxproc
                nproc(i+count) = maxproc
                color(i+count) = my_color
                key(i+count) = i-1
             end do
             count = count + maxproc
             npes1 = npes1 - maxproc
             numdir1 = numdir1 -1
          end do

          !!do i=1,npes
          !!   key(i) = mod(i-1,maxproc)
          !!   color(i) = (i-key(i))/maxproc
          !!end do
          !!write(*,'("nproc:",50(1x,i2))') nproc
          !!write(*,'("color:",50(1x,i2))') color
          !!write(*,'("key:",50(1x,i2))') key

          do i=1,npes
             idir = color(i)
             buf = workdir_array(idir)
             !!$write(*,*) buf(1:lenc)
             hispe = i-1
             if(mype /= hispe) then
                call MPI_Send(workdir_array(idir),50,MPI_CHARACTER &
                           & ,hispe,1,MPI_COMM_WORLD,ierr)
                call MPI_Send(color(i),1,MPI_INTEGER &
                           & ,hispe,1,MPI_COMM_WORLD,ierr)
                call MPI_Send(key(i),1,MPI_INTEGER &
                           & ,hispe,1,MPI_COMM_WORLD,ierr)
             else
                workdir = workdir_array(idir)
                my_color = color(i)
                my_key = key(i)
             end if
          end do

          deallocate(workdir_array)
       else
          call MPI_Recv(workdir,50,MPI_CHARACTER, &
                   & MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,istatus,ierr)
          call MPI_Recv(my_color,1,MPI_INTEGER, &
                   & MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,istatus,ierr)
          call MPI_Recv(my_key,1,MPI_INTEGER, &
                   & MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,istatus,ierr)
       end if
       call MPI_barrier(MPI_COMM_WORLD,ierr)
       !!write(*,*) 'workdir=',trim(workdir)
       !!my_color = mype
       !!my_key = 0
       !if(firstcall) call mpi_comm_split(MPI_COMM_WORLD,my_color,my_key,MPI_CommGroup,ierr)
       if(.not.firstcall) call mpi_comm_free(MPI_CommGroup, ierr)
       call mpi_comm_split(MPI_COMM_WORLD,my_color,my_key,MPI_CommGroup,ierr)
       call mpi_comm_size(MPI_CommGroup, npes, ierr)
       call mpi_comm_rank(MPI_CommGroup, mype, ierr)
       !!write(*,*) 'MPI_CommGroup=',MPI_CommGroup
       !!write(*,*) 'npes=',npes
       !!write(*,*) 'myes=',mype
       !!write(*,*) 'workdir=',trim(workdir)

       !! debug
       !!call mpi_Finalize(ierr)
       !!stop 'group_creation end'

    end if
!! Group creation

    else

    call mpi_comm_size(MPI_CommGroup, npes, ierr)
    call mpi_comm_rank(MPI_CommGroup, mype, ierr)

    endif

    firstcall = .false.

!!$    if(npes > 1) then
!!$       write(name_mype,'(i3)') mype
!!$       do i = 1, 3
!!$          if(name_mype(i:i) == ' ') name_mype(i:i)='0'
!!$       enddo
!!$
!!$       if(mype==0) then
!!$          do js = 0, 999
!!$             write(name_jobstep,'(i3)') js
!!$             do i = 1, 3
!!$                if(name_jobstep(i:i) == ' ') name_jobstep(i:i)='0'
!!$             enddo
!!$             F_OUT = 'output'//name_jobstep
!!$             inquire(file=F_OUT, exist=existence)
!!$             if(.not.existence) goto 1001
!!$          end do
!!$1001      if(js >=1000) js = 999
!!$       end if
!!$       if(npes > 1) call mpi_bcast(js,1,mpi_integer,0,mpi_comm_world,ierr)
!!$       write(name_jobstep,'(i3)') js
!!$       do i = 1, 3
!!$          if(name_jobstep(i:i) == ' ') name_jobstep(i:i)='0'
!!$       enddo
!!$
!!$       if(mype==0) then
!!$          F_OUT = 'output'//name_jobstep
!!$       else
!!$          F_OUT = 'output'//name_jobstep//'_'//name_mype
!!$       end if
!!$    else if(npes == 1) then
!!$       do js = 0, 999
!!$          write(name_jobstep,'(i3)') js
!!$          do i = 1, 3
!!$             if(name_jobstep(i:i) == ' ') name_jobstep(i:i)='0'
!!$          enddo
!!$          F_OUT = 'output'//name_jobstep
!!$          inquire(file=F_OUT, exist=existence)
!!$          if(.not.existence) goto 1000
!!$       end do
!!$1000   continue
!!$    end if
!!$    open(6, file=F_OUT, status='unknown', form='formatted')

  end subroutine m_Parallel_init_comm_world

  subroutine m_Parallel_init_comm_group(mpi_comm_root)
!!$    character :: F_OUT*13, name_mype*3
!!$    character :: name_jobstep*3
!!$    integer ::   i, js
!!$    logical ::   existence
    integer,  intent(in) :: mpi_comm_root
    integer :: inimpi

    logical, save :: firstcall = .true.

    workdir = ''

    inimpi=ON

    if(inimpi==ON) then

      if(.not.firstcall) call mpi_comm_free(MPI_CommGroup, ierr)
      call mpi_comm_dup(mpi_comm_root,MPI_CommGroup,ierr)
      call MPI_barrier(MPI_CommGroup,ierr)
      call mpi_comm_size(MPI_CommGroup, npes, ierr)
      call mpi_comm_rank(MPI_CommGroup, mype, ierr)
                                    __TIMER_START(17)
#ifdef FJ_TIMER
                 call timer_sta(21)
#endif


    else

      call mpi_comm_size(MPI_CommGroup, npes, ierr)
      call mpi_comm_rank(MPI_CommGroup, mype, ierr)

    endif

    firstcall = .false.

  end subroutine m_Parallel_init_comm_group

  subroutine m_Parallel_init_mpi(printable)
    logical, intent(in) :: printable
    integer   :: i
    integer, dimension(8) :: ipresent_time

!!$    call mpi_init(ierr)
!!$    call mpi_comm_size(mpi_comm_world, npes, ierr)
!!$    call mpi_comm_rank(mpi_comm_world, mype, ierr)

    call date_and_time(values=ipresent_time)

    if(printable) then
       write(6,'(" !|| program start    ",i2,":",i2,":",i2,"  ",i2,"/",i2,"/",i4)') &
            & ipresent_time(5:7),ipresent_time(3),ipresent_time(2),ipresent_time(1)

       write(6,'(" npes = ",i6)') npes
       write(6,'(" mype = ",i3," mpi_comm_world = ",i3)') mype,MPI_CommGroup
       write(6,'(" mype = ",i3," mpi_comm_world = ",i3)') mype,MPI_CommGroup
    end if

    call mpi_barrier(MPI_CommGroup,ierr)
    if(printable) write(6,'(" mpi_comm_world = ",i5," mype = ",i5)') MPI_CommGroup, mype
    call mpi_barrier(MPI_CommGroup,ierr)

  end subroutine m_Parallel_init_mpi




  subroutine m_Parallel_init_mpi_nbmx(nfout,ipri,printable,nbmx,kg1,ngnode_nbmx,flag_mpi,flag_mpi_k)
    integer, intent(in) :: nfout,ipri,nbmx, kg1, ngnode_nbmx
    logical, intent(in) :: printable,flag_mpi,flag_mpi_k
    integer             :: i, j, ip, icolor, key, np
    integer             :: newpes, newmype
    logical, save       :: firstcall=.true.
                                                  __TIMER_SUB_START(1235)

    if(.not.flag_mpi) then
! (( ng_nbmx, myrank_nbmx, nrank_nbmx, np_nbmx, mp_nbmx, ista_nbmx, iend_nbmx
       ng_nbmx = mype
       myrank_nbmx = 0
       nrank_nbmx  = 1
       np_nbmx = nbmx
       mp_nbmx = nbmx
       ista_nbmx = 1
       iend_nbmx = nbmx
       nbmx_ext = nbmx
       if(ipri >= 2 .and. printable) then
          write(nfout,'(" !|| -- parallelization parameters for nbmx --")')
          write(nfout,'(" !||  - ng_nbmx     = ",i12)') ng_nbmx
          write(nfout,'(" !||  - myrank_nbmx = ",i12)') myrank_nbmx
          write(nfout,'(" !||  - nrank_nbmx  = ",i12)') nrank_nbmx
          write(nfout,'(" !||  - nbmx, nbmx_ext = ",2i12)') nbmx, nbmx_ext
          write(nfout,'(" !||  - np_nbmx     = ",i12)') np_nbmx
          write(nfout,'(" !||  - mp_nbmx     = ",i12)') mp_nbmx
          write(nfout,'(" !||  - ista_nbmx   = ",i12)') ista_nbmx
          write(nfout,'(" !||  - iend_nbmx   = ",i12)') iend_nbmx
       end if
    else
! (( ng_nbmx, myrank_nbmx, nrank_nbmx, np_nbmx, mp_nbmx, ista_nbmx, iend_nbmx
!      , nis_nbmx, nie_nbmx, nel_nbmx, idisp_nbmx))
       ng_nbmx     = floor(mod(mype+1, ngnode_nbmx) - 0.5)
       if(ng_nbmx < 0) ng_nbmx = ng_nbmx + ngnode_nbmx
       myrank_nbmx = floor((mype+0.5)/ngnode_nbmx)
       nrank_nbmx  = (npes - ng_nbmx-1)/ngnode_nbmx + 1

       if(.not.allocated(mpi_nbmx_world)) allocate(mpi_nbmx_world(0:ngnode_nbmx-1))

       allocate(nis_nbmx(0:nrank_nbmx-1))
       allocate(nie_nbmx(0:nrank_nbmx-1))
       allocate(nel_nbmx(0:nrank_nbmx-1))
       allocate(idisp_nbmx(0:nrank_nbmx-1))

       call set_block_range4allgather(nbmx,nrank_nbmx,nel_nbmx,nis_nbmx,nie_nbmx,idisp_nbmx)
       ista_nbmx = nis_nbmx(myrank_nbmx)
       iend_nbmx = nie_nbmx(myrank_nbmx)
       np_nbmx = nel_nbmx(myrank_nbmx)
       mp_nbmx = maxval(nel_nbmx)
       nbmx_ext = mp_nbmx*nrank_nbmx

       if(ipri > 1 .and. printable) then
          write(nfout,'(" !|| -- parallelization parameters for nbmx --")')
          write(nfout,'(" !||  - ng_nbmx     = ",i12)') ng_nbmx
          write(nfout,'(" !||  - myrank_nbmx = ",i12)') myrank_nbmx
          write(nfout,'(" !||  - nrank_nbmx  = ",i12)') nrank_nbmx
          write(nfout,'(" !||  - nbmx, nbmx_ext = ",2i12)') nbmx, nbmx_ext
          write(nfout,'(" !||  - np_nbmx     = ",i12)') np_nbmx
          write(nfout,'(" !||  - mp_nbmx     = ",i12)') mp_nbmx
          write(nfout,'(" !||  - ista_nbmx   = ",i12)') ista_nbmx
          write(nfout,'(" !||  - iend_nbmx   = ",i12)') iend_nbmx
          write(nfout,'(" !||   ( nis_nbmx )",8i10)') (nis_nbmx(i),i=0,nrank_nbmx-1)
          write(nfout,'(" !||   ( nie_nbmx )",8i10)') (nie_nbmx(i),i=0,nrank_nbmx-1)
          write(nfout,'(" !||   ( nel_nbmx )",8i10)') (nel_nbmx(i),i=0,nrank_nbmx-1)
          write(nfout,'(" !||   (idisp_nbmx)",8i10)') (idisp_nbmx(i),i=0,nrank_nbmx-1)
       end if

! (( mpi_nbmx_world ))
       do j = 0, ngnode_nbmx-1
          icolor = 0
          key = 0
          if(ng_nbmx == j) then
             icolor = 1
             key = myrank_nbmx
          end if
!          if(firstcall) call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_nbmx_world(j),ierr)
          if(.not.firstcall) call mpi_comm_free(mpi_nbmx_world(j), ierr)
          call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_nbmx_world(j),ierr)
          call mpi_comm_size(mpi_nbmx_world(j), newpes, ierr)
          call mpi_comm_rank(mpi_nbmx_world(j), newmype,ierr)
       end do
    end if

    if(.not.flag_mpi_k) then
! (( ng_nbmx_k, myrank_nbmx_k, nrank_nbmx_k, np_nbmx_k, mp_nbmx_k, ista_nbmx_k, iend_nbmx_k,
!    ista_kg1_k, iend_kg1_k ))
       ng_nbmx_k = myrank_k
       myrank_nbmx_k = 0
       nrank_nbmx_k  = 1
       np_nbmx_k = nbmx
       mp_nbmx_k = nbmx
       ista_nbmx_k = 1
       iend_nbmx_k = nbmx

       ista_kg1_k = 1
       iend_kg1_k = kg1
       kg1_ext = kg1
    else
! (( ng_nbmx_k, myrank_nbmx_k, nrank_nbmx_k, np_nbmx_k, mp_nbmx_k, nel_nbmx_k,
!  , ista_nbmx_k, iend_nbmx_k,  nis_kg1_k, nie_kg1_k, nel_kg1_k, ista_kg1_k, iend_kg1_k ))
       ng_nbmx_k     = floor(mod(myrank_e+1, ngnode_nbmx) - 0.5)
       if(ng_nbmx_k < 0) ng_nbmx_k = ng_nbmx_k + ngnode_nbmx
       myrank_nbmx_k = floor((myrank_e+0.5)/ngnode_nbmx)
       nrank_nbmx_k  = (nrank_e - ng_nbmx-1)/ngnode_nbmx + 1

       allocate(nel_nbmx_k(0:nrank_nbmx_k-1))
       allocate(nis_nbmx_k(0:nrank_nbmx_k-1))
       allocate(nie_nbmx_k(0:nrank_nbmx_k-1))
       allocate(idisp_nbmx_k(0:nrank_nbmx_k-1))
       if(.not.allocated(mpi_nbmx_world_k)) then
       allocate(mpi_nbmx_world_k(0:ngnode_nbmx-1))
       endif

       call set_block_range4allgather(nbmx,nrank_nbmx_k,nel_nbmx_k,nis_nbmx_k,nie_nbmx_k,idisp_nbmx_k)
       ista_nbmx_k = nis_nbmx_k(myrank_nbmx_k)
       iend_nbmx_k = nie_nbmx_k(myrank_nbmx_k)
       np_nbmx_k = nel_nbmx_k(myrank_nbmx_k)
       mp_nbmx_k = maxval(nel_nbmx_k)
       np = mp_nbmx_k*nrank_nbmx_k
       if(np .gt. nbmx_ext) nbmx_ext = np

! (( np_kg1_k, mp_kg1_k, nel_kg1_k, ista_kg1_k, iend_kg1_k, idisp_kg1_k ))
       allocate(nel_kg1_k(0:nrank_nbmx_k-1))
       allocate(nis_kg1_k(0:nrank_nbmx_k-1))
       allocate(nie_kg1_k(0:nrank_nbmx_k-1))
       allocate(idisp_kg1_k(0:nrank_nbmx_k-1))

       call set_block_range4allgather(kg1,nrank_nbmx_k,nel_kg1_k,nis_kg1_k,nie_kg1_k,idisp_kg1_k)
       np_kg1_k = nel_kg1_k(myrank_nbmx_k)
       mp_kg1_k = maxval(nel_kg1_k)
       ista_kg1_k = nis_kg1_k(myrank_nbmx_k)
       iend_kg1_k = nie_kg1_k(myrank_nbmx_k)
       kg1_ext = mp_kg1_k*nrank_nbmx_k

       if(ipri > 1 .and. printable) then
          write(nfout,'(" !|| -- parallelization parameters for nbmx in mpi_k_world--")')
          write(nfout,'(" !||  - ng_nbmx_k     = ",i12)') ng_nbmx_k
          write(nfout,'(" !||  - myrank_nbmx_k = ",i12)') myrank_nbmx_k
          write(nfout,'(" !||  - nrank_nbmx_k  = ",i12)') nrank_nbmx_k
          write(nfout,'(" !||  - nbmx, nbmx_ext = ",2i12)') nbmx, nbmx_ext
          write(nfout,'(" !||  - np_nbmx_k     = ",i12)') np_nbmx_k
          write(nfout,'(" !||  - mp_nbmx_k     = ",i12)') mp_nbmx_k
          write(nfout,'(" !||  - ista_nbmx_k   = ",i12)') ista_nbmx_k
          write(nfout,'(" !||  - iend_nbmx_k   = ",i12)') iend_nbmx_k
          write(nfout,'(" !||   ( nis_nbmx_k )",8i10)') (nis_nbmx_k(i),i=0,nrank_nbmx_k-1)
          write(nfout,'(" !||   ( nie_nbmx_k )",8i10)') (nie_nbmx_k(i),i=0,nrank_nbmx_k-1)
          write(nfout,'(" !||   ( nel_nbmx_k )",8i10)') (nel_nbmx_k(i),i=0,nrank_nbmx_k-1)
          write(nfout,'(" !||   (idisp_nbmx_k)",8i10)') (idisp_nbmx_k(i),i=0,nrank_nbmx_k-1)
          write(nfout,'(" !||  - kg1, kg1_ext  = ",2i12)') kg1, kg1_ext
          write(nfout,'(" !||  - np_kg1_k      = ",i12)') np_kg1_k
          write(nfout,'(" !||  - mp_kg1_k      = ",i12)') mp_kg1_k
          write(nfout,'(" !||  - ista_kg1_k    = ",i12)') ista_kg1_k
          write(nfout,'(" !||  - iend_kg1_k    = ",i12)') iend_kg1_k
          write(nfout,'(" !||   ( nis_kg1_k )",8i10)') (nis_kg1_k(i),i=0,nrank_nbmx_k-1)
          write(nfout,'(" !||   ( nie_kg1_k )",8i10)') (nie_kg1_k(i),i=0,nrank_nbmx_k-1)
          write(nfout,'(" !||   ( nel_kg1_k )",8i10)') (nel_kg1_k(i),i=0,nrank_nbmx_k-1)
          write(nfout,'(" !||   (idisp_kg1_k)",8i10)') (idisp_kg1_k(i),i=0,nrank_nbmx_k-1)
       end if

! (( mpi_nbmx_world ))
       do i = 0, nrank_k-1
          if(i == myrank_k) then
             do j = 0, ngnode_nbmx-1
                icolor = 0
                key = 0
                if(ng_nbmx_k == j) then
                   icolor = 1
                   key = myrank_nbmx_k
                end if
!                if(firstcall) call mpi_comm_split(mpi_k_world(myrank_k), icolor, key, mpi_nbmx_world_k(j),ierr)
                if(.not.firstcall) call mpi_comm_free(mpi_nbmx_world_k(j), ierr)
                call mpi_comm_split(mpi_k_world(myrank_k), icolor, key, mpi_nbmx_world_k(j),ierr)
                call mpi_comm_size(mpi_nbmx_world_k(j), newpes, ierr)
                call mpi_comm_rank(mpi_nbmx_world_k(j), newmype,ierr)
             end do
          end if
       end do
    end if
    firstcall = .false.
                                                  __TIMER_SUB_STOP(1235)
  end subroutine m_Parallel_init_mpi_nbmx

  subroutine m_Parallel_init_mpi_urec_hsr(nfout,ipri,nsize_rho_hsr)
    ! Coded by T. Ymasaki, 2023/07/07
    integer, intent(in) :: nfout, ipri, nsize_rho_hsr
    integer, allocatable, dimension(:) :: is_hsr, ie_hsr
    integer :: iwork, k
    allocate(is_hsr(0:npes-1),ie_hsr(0:npes-1))
    iwork = ( nsize_rho_hsr - 1)/nrank_chg + 1
    do k = 0, nrank_chg-1
       is_hsr(k) = min(k*iwork+1,nsize_rho_hsr+1)
       ie_hsr(k) = min(is_hsr(k)+iwork-1,nsize_rho_hsr)
    end do
    ista_urec_hsr = is_hsr(myrank_chg)
    iend_urec_hsr = ie_hsr(myrank_chg)
    if(ipri >= 1) then
       write(nfout,'(" !|| ista_urec_hsr = ",i8, " iend_urec_hsr = ",i8, " myrank_chg = ",i4)') &
            & ista_urec_hsr, iend_urec_hsr, myrank_chg
       call flush(nfout)
    end if
    deallocate(is_hsr, ie_hsr)
    ista_and_iend_urec_hsr_set = .true.
  end subroutine m_Parallel_init_mpi_urec_hsr

  subroutine set_block_range4allgather(ne,np,nel_p,nis_p,nie_p,idisp_p)
    integer, intent(in)                     :: ne ! number of total elements
    integer, intent(in)                     :: np ! number of ranks (or processors)
    integer, intent(out), dimension(0:np-1) :: nel_p,nis_p,nie_p,idisp_p

    integer :: i, j, npf

    nel_p = 0
    i = ne/np
    j = i*np
    if(j == ne) then
       nel_p = i
    else
       npf = ne/(i+1)
       nel_p(0:npf-1) = i+1
       nel_p(npf)     = ne - (i+1)*npf
    end if
    nis_p(0) = 1
    do i = 1, np-1
       nis_p(i)   = nis_p(i-1) + nel_p(i-1)
       nie_p(i-1) = nis_p(i)   - 1
    end do
    nie_p(np-1) = ne
    do i = 0, np-1
       idisp_p(i) = nis_p(i)-1
    end do
  end subroutine set_block_range4allgather

! === necessary to make 3D_Parallel, too!!! by tkato ===========================

! ==============================================================================
  subroutine m_Parallel_dealloc_mpi_fft_box
    if(allocated(nis_fft_X_z)) deallocate(nis_fft_X_z)
    if(allocated(nis_fft_X_y)) deallocate(nis_fft_X_y)
    if(allocated(nis_fft_Z_x)) deallocate(nis_fft_Z_x)
    if(allocated(nis_fft_Z_y)) deallocate(nis_fft_Z_y)
    if(allocated(nis_fft_Y_x)) deallocate(nis_fft_Y_x)
    if(allocated(nis_fft_Y_z)) deallocate(nis_fft_Y_z)
    if(allocated(nie_fft_X_z)) deallocate(nie_fft_X_z)
    if(allocated(nie_fft_X_y)) deallocate(nie_fft_X_y)
    if(allocated(nie_fft_Z_x)) deallocate(nie_fft_Z_x)
    if(allocated(nie_fft_Z_y)) deallocate(nie_fft_Z_y)
    if(allocated(nie_fft_Y_x)) deallocate(nie_fft_Y_x)
    if(allocated(nie_fft_Y_z)) deallocate(nie_fft_Y_z)

    if(allocated(map_fft_x)) deallocate(map_fft_x)
    if(allocated(map_fft_y)) deallocate(map_fft_y)
    if(allocated(map_fft_z)) deallocate(map_fft_z)
    if(allocated(nel_fft_x)) deallocate(nel_fft_x)
    if(allocated(nel_fft_y)) deallocate(nel_fft_y)
    if(allocated(nel_fft_z)) deallocate(nel_fft_z)
    if(allocated(mp_fft_x)) deallocate(mp_fft_x)
    if(allocated(mp_fft_y)) deallocate(mp_fft_y)
    if(allocated(mp_fft_z)) deallocate(mp_fft_z)

    if(allocated(nel_fftcd_x)) deallocate(nel_fftcd_x)
    if(allocated(nel_fftcd_y)) deallocate(nel_fftcd_y)
    if(allocated(nel_fftcd_z)) deallocate(nel_fftcd_z)
    if(allocated(map_fftcd_x)) deallocate(map_fftcd_x)
    if(allocated(map_fftcd_y)) deallocate(map_fftcd_y)
    if(allocated(map_fftcd_z)) deallocate(map_fftcd_z)

    if(allocated(mp_fftcd_x)) deallocate(mp_fftcd_x)
    if(allocated(mp_fftcd_y)) deallocate(mp_fftcd_y)
    if(allocated(mp_fftcd_z)) deallocate(mp_fftcd_z)

    if(allocated(nis_fftcd_X_z)) deallocate(nis_fftcd_X_z)
    if(allocated(nis_fftcd_X_y)) deallocate(nis_fftcd_X_y)
    if(allocated(nie_fftcd_X_z)) deallocate(nie_fftcd_X_z)
    if(allocated(nie_fftcd_X_y)) deallocate(nie_fftcd_X_y)

    if(allocated(nis_fftcd_Y_x)) deallocate(nis_fftcd_Y_x)
    if(allocated(nis_fftcd_Y_z)) deallocate(nis_fftcd_Y_z)
    if(allocated(nie_fftcd_Y_x)) deallocate(nie_fftcd_Y_x)
    if(allocated(nie_fftcd_Y_z)) deallocate(nie_fftcd_Y_z)

    if(allocated(nis_fftcd_Z_x)) deallocate(nis_fftcd_Z_x)
    if(allocated(nis_fftcd_Z_y)) deallocate(nis_fftcd_Z_y)
    if(allocated(nie_fftcd_Z_x)) deallocate(nie_fftcd_Z_x)
    if(allocated(nie_fftcd_Z_y)) deallocate(nie_fftcd_Z_y)

    if(allocated(nis_fftcd_X_x)) deallocate(nis_fftcd_X_x)
    if(allocated(nie_fftcd_X_x)) deallocate(nie_fftcd_X_x)
  end subroutine m_Parallel_dealloc_mpi_fft_box

  subroutine m_Parallel_dealloc_mpi_kngp_B()
    if(allocated(is_kngp_B))  deallocate(is_kngp_B)
    if(allocated(ie_kngp_B))  deallocate(ie_kngp_B)
    if(allocated(nel_kngp_B)) deallocate(nel_kngp_B)
  end subroutine m_Parallel_dealloc_mpi_kngp_B
  subroutine m_Parallel_dealloc_mpi_elec()

    if(allocated(map_e)) deallocate(map_e)
    if(allocated(map_z)) deallocate(map_z)
    if(allocated(map_k)) deallocate(map_k)
    if(allocated(map_ek)) deallocate(map_ek)
    if(allocated(map_s))  deallocate(map_s)
!fj    deallocate(mpi_k_world)
    if(allocated(nis_e)) deallocate(nis_e)
    if(allocated(nie_e)) deallocate(nie_e)
    if(allocated(nel_e)) deallocate(nel_e)
    if(allocated(idisp_e)) deallocate(idisp_e)
    if(allocated(nis_g1)) deallocate(nis_g1)
    if(allocated(nie_g1)) deallocate(nie_g1)
    if(allocated(nel_g1)) deallocate(nel_g1)
    if(allocated(nis_k)) deallocate(nis_k)
    if(allocated(nie_k)) deallocate(nie_k)
    if(allocated(nel_k)) deallocate(nel_k)

#ifdef NEC_TUNE_SOFT
    if(allocated(ista_e_smp)) deallocate(ista_e_smp)
    if(allocated(iend_e_smp)) deallocate(iend_e_smp)
#endif

!for 3D
    if(allocated(lrank)) deallocate(lrank)
    if(allocated(nbsn)) deallocate(nbsn)
    if(allocated(nbsn_sta)) deallocate(nbsn_sta)
    if(allocated(nbsn_end)) deallocate(nbsn_end)
    if(allocated(nbs_sta)) deallocate(nbs_sta)
    if(allocated(nbs_end)) deallocate(nbs_end)
    if(allocated(neg_g)) deallocate(neg_g)
    if(allocated(neg_g_all)) deallocate(neg_g_all)

    if(allocated(neg_gg)) deallocate(neg_gg)
    if(allocated(neg_gg_all)) deallocate(neg_gg_all)
    if(allocated(nis_eg)) deallocate(nis_eg)
    if(allocated(nie_eg)) deallocate(nie_eg)
    if(allocated(nel_eg)) deallocate(nel_eg)

!    if(allocated(ball_buff)) deallocate(ball_buff)
!    if(allocated(ball_addr)) deallocate(ball_addr)

    if(allocated(nis_spin)) deallocate(nis_spin)
    if(allocated(nie_spin)) deallocate(nie_spin)
    if(allocated(nel_spin)) deallocate(nel_spin)
  end subroutine m_Parallel_dealloc_mpi_elec


  subroutine m_Parallel_end_mpi
                                                  __TIMER_STOP
    call mpi_finalize(ierr)
  end subroutine m_Parallel_end_mpi

  subroutine set_block_range(ne,np,nel_p,nis_p,nie_p,set_mapping_func,map_p)
    integer, intent(in)                     :: ne ! number of total elements
    integer, intent(in)                     :: np ! number of ranks (or processors)
    integer, intent(out), dimension(0:np-1) :: nel_p,nis_p,nie_p
    logical, intent(in)                     :: set_mapping_func
    integer, intent(out), optional, dimension(ne)  :: map_p

    integer :: j,i

    if(np == 0) call phase_error_with_msg(6,' np == 0',__LINE__,__FILE__)
    nel_p = ne/np
    j = mod(ne,np)
    do i = 0, j-1
       nel_p(i) = nel_p(i) + 1
    end do

    nis_p(0) = 1
    do i = 1, np-1
       nis_p(i)   = nis_p(i-1) + nel_p(i-1)
       nie_p(i-1) = nis_p(i) - 1
    end do
    nie_p(np-1) = ne

    if(set_mapping_func) then
       j = 0
       do i = 1, ne
          if(nie_p(j) < i) j = j + 1
          map_p(i) = j
       end do
    end if

  end subroutine set_block_range

  subroutine m_Parallel_init_mpi_kngp_prev_3D(kngp,comm_for_chg)
    integer, intent(in) :: kngp
    logical, intent(in) :: comm_for_chg
    integer :: iwork, i
    integer, allocatable, dimension(:) :: is_kngp,ie_kngp,nel_kngp
    if(comm_for_chg)then
      mpi_chg_world = MPI_CommGroup
      nrank_chg = npes
      myrank_chg = mype
    else
      mpi_chg_world = mpi_ke_world
      nrank_chg = nrank_g
      myrank_chg = myrank_g
    endif
    allocate(is_kngp(0:nrank_chg-1))
    allocate(ie_kngp(0:nrank_chg-1))
    allocate(nel_kngp(0:nrank_chg-1))
    iwork = ( kngp - 1 ) / nrank_chg + 1
    do i = 0, nrank_chg-1
       is_kngp(i) = min(i*iwork+1, kngp+1)
       ie_kngp(i) = min(is_kngp(i)+iwork-1, kngp)
       nel_kngp(i) = ie_kngp(i) - is_kngp(i) + 1
    enddo
    ista_kngp_prev = is_kngp(myrank_chg)
    iend_kngp_prev = ie_kngp(myrank_chg)
    np_kngp_prev   = nel_kngp(myrank_chg)
    deallocate(is_kngp)
    deallocate(ie_kngp)
    deallocate(nel_kngp)
  end subroutine m_Parallel_init_mpi_kngp_prev_3D

  subroutine m_Parallel_init_mpi_cdfft(nfout,ipri,ggacmp_parallel)
    integer, intent(in)  :: nfout,ipri,ggacmp_parallel
    integer :: i, iprilevel

    if(ggacmp_parallel == ON) then
       call split_into_ggablock_cdfft(npes,nrank_ggacmp,npes_cdfft,nrest_cdfft)
    else
       nrank_ggacmp = 1
       npes_cdfft   = npes
       nrest_cdfft  = 0
    end if
    if(ipri > 1) write(nfout,'(" !|| ggacmp_parallel = ",i3)') ggacmp_parallel
    if(ggacmp_parallel==ON) then
       iprilevel=1
    else
       iprilevel=2
    end if
    if(ipri > iprilevel) write(nfout,'(" !|| nrank_ggacmp, npes_cdfft, nrest_cdfft",/ &
         & ," !||   = ",3i10," <<m_Parallel_init_mpi_cdfft>>")') &
         & nrank_ggacmp, npes_cdfft, nrest_cdfft
    myrank_ggacmp = mype/npes_cdfft
    myrank_cdfft  = mype - myrank_ggacmp*npes_cdfft
    if(ipri > iprilevel) write(nfout,'(" !|| myrank_ggacmp, myrank_cdfft",/ &
         & ," !|| = ",2i10," <<m_Parallel_init_mpi_cdfft>>")') &
         & myrank_ggacmp, myrank_cdfft
    if(allocated(map_pe2ggacmp)) deallocate(map_pe2ggacmp)
    if(allocated(map_pe2cdfft)) deallocate(map_pe2cdfft)
    allocate(map_pe2ggacmp(0:npes-1))
    allocate(map_pe2cdfft(0:npes-1))
    do i = 0, npes-1
       map_pe2ggacmp(i) = i/npes_cdfft
       map_pe2cdfft(i) = i - map_pe2ggacmp(i)*npes_cdfft
    end do
    if(ipri > iprilevel) write(nfout,'(" !|| map_pe2ggacmp",/,99( " !|| ",10i5,/))') &
         & (map_pe2ggacmp(i),i=0,npes-1)
    if(ipri > iprilevel) write(nfout,'(" !|| map_pe2cdfft",/,99( " !|| ",10i5,/))') &
         & (map_pe2cdfft(i),i=0,npes-1)
  end subroutine m_Parallel_init_mpi_cdfft

  subroutine m_Parallel_dealloc_mpi_gga()
    if(allocated(map_ggacmp))  deallocate(map_ggacmp)
    if(allocated(nis_fftp))    deallocate(nis_fftp)
    if(allocated(nie_fftp))    deallocate(nie_fftp)
    if(allocated(nel_fftp))    deallocate(nel_fftp)
    if(allocated(idisp_fftp))  deallocate(idisp_fftp)
    if(allocated(nis_fftph))   deallocate(nis_fftph)
    if(allocated(nie_fftph))   deallocate(nie_fftph)
    if(allocated(nel_fftph))   deallocate(nel_fftph)
    if(allocated(idisp_fftph)) deallocate(idisp_fftph)
    if(allocated(nis_sfftp))   deallocate(nis_sfftp)
    if(allocated(nie_sfftp))   deallocate(nie_sfftp)
    if(allocated(nel_sfftp))   deallocate(nel_sfftp)
    if(allocated(idisp_sfftp)) deallocate(idisp_sfftp)
    if(allocated(nis_sfftph))  deallocate(nis_sfftph)
    if(allocated(nie_sfftph))  deallocate(nie_sfftph)
    if(allocated(nel_sfftph))  deallocate(nel_sfftph)
  end subroutine m_Parallel_dealloc_mpi_gga

  subroutine m_Parallel_init_mpi_gga(nfout,ipri,printable,nfftp,nfftps)
    integer, intent(in)  :: nfout,ipri,nfftp,nfftps
    logical, intent(in)  :: printable
    integer              :: nfftph, i
    logical              :: set_mapping_func
    integer              :: icolor, key, j, newpes, newmype
    logical, save        :: firstcall = .true.
                                                  __TIMER_SUB_START(1236)
                                                  __TIMER_SUB_STOP(1236)
!    return

    allocate(map_ggacmp(3)); map_ggacmp = 0
    allocate(nis_fftp(0:npes_cdfft-1)); nis_fftp = 0
    allocate(nie_fftp(0:npes_cdfft-1)); nie_fftp = 0
    allocate(nel_fftp(0:npes_cdfft-1)); nel_fftp = 0
    allocate(idisp_fftp(0:npes_cdfft-1)); idisp_fftp = 0
    allocate(nis_fftph(0:npes_cdfft-1));nis_fftph = 0
    allocate(nie_fftph(0:npes_cdfft-1)); nie_fftph = 0
    allocate(nel_fftph(0:npes_cdfft-1)); nel_fftph = 0
    allocate(idisp_fftph(0:npes_cdfft-1)); idisp_fftph = 0
    if(.not.allocated(mpi_cdfft_world))then
    allocate(mpi_cdfft_world(0:nrank_ggacmp-1)); mpi_cdfft_world = 0
    endif
    if(.not.allocated(mpi_ggacmp_cross_world)) then
    allocate(mpi_ggacmp_cross_world(0:npes_cdfft-1)); mpi_ggacmp_cross_world = 0
    endif

    if(nrank_ggacmp > 1) then
       do i = 1, 3
          map_ggacmp(i) = i-1
       end do
    else
       map_ggacmp(1:3) = 0
    end if

    set_mapping_func = .false.
    nfftph = nfftp/2
    if(ipri > 1) write(nfout,'(" !|| nfftp, nfftph = ",2i10," <<m_Parallel_init_mpi_gga>>")') nfftp, nfftph
    call set_block_range(nfftph,npes_cdfft,nel_fftph,nis_fftph,nie_fftph,set_mapping_func)

    ista_fftph = nis_fftph(myrank_cdfft)
    iend_fftph = nie_fftph(myrank_cdfft)

    if(ipri > 1) write(nfout,'(" !|| ista_fftph, iend_fftph = ",2i10 &
         & ," <<m_Parallel_init_mpi_gga>>")') ista_fftph, iend_fftph

    do i = 0, npes_cdfft -1
       nis_fftp(i) = nis_fftph(i)*2 - 1
       nie_fftp(i) = nie_fftph(i)*2
       nel_fftp(i) = nel_fftph(i)*2
       idisp_fftp(i) = nis_fftp(i) - 1
       idisp_fftph(i) = nis_fftph(i) - 1
    end do
    ista_fftp = nis_fftp(myrank_cdfft)
    iend_fftp = nie_fftp(myrank_cdfft)
    mp_fftp = maxval(nel_fftp)
    np_fftp = nel_fftp(myrank_cdfft)

    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| --- nis_fftp,nie_fftp,nel_fftp,idisp_fftp ---")')
       write(nfout,'(" !||  npes_cdfft = ",i12)') npes_cdfft
       write(nfout,'(" !|| ( npe_cdfft )",10i12)') (i,i=1,npes_cdfft)
       write(nfout,'(" !|| ( nis_fftp )",10i12)') (nis_fftp(i),i=0,npes_cdfft-1)
       write(nfout,'(" !|| ( nie_fftp )",10i12)') (nie_fftp(i),i=0,npes_cdfft-1)
       write(nfout,'(" !|| ( nel_fftp )",10i12)') (nel_fftp(i),i=0,npes_cdfft-1)
       write(nfout,'(" !|| (idisp_fftp)",10i12)') (idisp_fftp(i),i=0,npes_cdfft-1)
       write(nfout,'(" !|| ( ista_fftp,  iend_fftp )",2i12)') ista_fftp, iend_fftp
       write(nfout,'(" !|| --- nis_fftph,nie_fftph,nel_fftph ---")')
       write(nfout,'(" !|| ( npe_cdfft )",10i12)') (i,i=1,npes_cdfft)
       write(nfout,'(" !|| ( nis_fftph)",10i12)') (nis_fftph(i),i=0,npes_cdfft-1)
       write(nfout,'(" !|| ( nie_fftph)",10i12)') (nie_fftph(i),i=0,npes_cdfft-1)
       write(nfout,'(" !|| ( nel_fftph)",10i12)') (nel_fftph(i),i=0,npes_cdfft-1)
       write(nfout,'(" !|| ( ista_fftph, iend_fftph)",2i12)') ista_fftph, iend_fftph
    end if

    allocate(nis_sfftp(0:npes-1)); nis_sfftp = 0
    allocate(nie_sfftp(0:npes-1)); nie_sfftp = 0
    allocate(nel_sfftp(0:npes-1)); nel_sfftp = 0
    allocate(idisp_sfftp(0:npes-1)); idisp_sfftp = 0
    allocate(nis_sfftph(0:npes-1)); nis_sfftph = 0
    allocate(nie_sfftph(0:npes-1)); nie_sfftph = 0
    allocate(nel_sfftph(0:npes-1)); nel_sfftph = 0
    set_mapping_func = .false.
    if(ipri > 1 .and. printable) write(nfout,'(" !|| nfftp = ",i10," <<m_Parallel_init_mpi_gga>>")') nfftp

!!$    call set_block_range(nfftp,npes,nel_sfftp,nis_sfftp,nie_sfftp,set_mapping_func)
!!$    ista_sfftp = nis_sfftp(mype)
!!$    iend_sfftp = nie_sfftp(mype)
!!$    idisp_sfftp(0) = 0
!!$    do i = 1, npes-1
!!$       idisp_sfftp(i) = idisp_sfftp(i-1) + nel_sfftp(i-1)
!!$    end do

    nfftph = nfftps/2
    call set_block_range(nfftph,npes,nel_sfftph,nis_sfftph,nie_sfftph,set_mapping_func)
    ista_sfftph = nis_sfftph(mype)
    iend_sfftph = nie_sfftph(mype)

    do i = 0, npes -1
       nis_sfftp(i) = nis_sfftph(i)*2 - 1
       nie_sfftp(i) = nie_sfftph(i)*2
       nel_sfftp(i) = nel_sfftph(i)*2
       idisp_sfftp(i) = nis_sfftp(i) - 1
    end do
    ista_sfftp = nis_sfftp(mype)
    iend_sfftp = nie_sfftp(mype)
    mp_sfftp = maxval(nel_sfftp)
    np_sfftp = nel_sfftp(mype)

    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| --- nis_sfftp,nie_sfftp,nel_sfftp,idisp_sfftp ---")')
       write(nfout,'(" !|| ( npe      )",10i12)') (i,i=1,npes)
       write(nfout,'(" !|| ( nis_sfftp )",10i12)') (nis_sfftp(i),i=0,npes-1)
       write(nfout,'(" !|| ( nie_sfftp )",10i12)') (nie_sfftp(i),i=0,npes-1)
       write(nfout,'(" !|| ( nel_sfftp )",10i12)') (nel_sfftp(i),i=0,npes-1)
       write(nfout,'(" !|| (idisp_sfftp)",10i12)') (idisp_sfftp(i),i=0,npes-1)
       write(nfout,'(" !|| ( ista_sfftp,  iend_sfftp )",2i12)') ista_sfftp, iend_sfftp
       write(nfout,'(" !|| --- nis_sfftph,nie_sfftph,nel_sfftph ---")')
       write(nfout,'(" !|| ( npe       )",10i12)') (i,i=1,npes)
       write(nfout,'(" !|| ( nis_sfftph )",10i12)') (nis_sfftph(i),i=0,npes-1)
       write(nfout,'(" !|| ( nie_sfftph )",10i12)') (nie_sfftph(i),i=0,npes-1)
       write(nfout,'(" !|| ( nel_sfftph )",10i12)') (nel_sfftph(i),i=0,npes-1)
       write(nfout,'(" !|| ( ista_sfftph, iend_sfftph)",2i12)') ista_sfftph, iend_sfftph
    end if

    do j = 0, nrank_ggacmp-1
       icolor = 0
       key = 0
       if(myrank_ggacmp == j) then
          icolor = 1
          key = myrank_cdfft
       end if
!!$       do i = 0, npes_cdfft-1
!!$          if(mype == i + npes_cdfft*j) then
!!$             icolor = 1
!!$             key = i
!!$          end if
!!$       end do
!       if(firstcall) call mpi_comm_split(MPI_CommGroup,icolor,key,mpi_cdfft_world(j),ierr)
       if(.not.firstcall) call mpi_comm_free(mpi_cdfft_world(j), ierr)
       call mpi_comm_split(MPI_CommGroup,icolor,key,mpi_cdfft_world(j),ierr)
       call mpi_comm_size(mpi_cdfft_world(j), newpes, ierr)
       call mpi_comm_rank(mpi_cdfft_world(j), newmype, ierr)
    end do

    do j = 0, npes_cdfft-1
       icolor = 0
       key = 0
       if(myrank_cdfft == j) then
          icolor = 1
          key = myrank_ggacmp
       end if
!!$       do i = 0, nrank_ggacmp-1
!!$          if(mype == i + nrank_ggacmp*j) then
!!$             icolor = 1
!!$             key = i
!!$          end if
!!$       end do
!       if(firstcall) call mpi_comm_split(MPI_CommGroup,icolor,key,mpi_ggacmp_cross_world(j),ierr)
       if(.not.firstcall) call mpi_comm_free(mpi_ggacmp_cross_world(j), ierr)
       call mpi_comm_split(MPI_CommGroup,icolor,key,mpi_ggacmp_cross_world(j),ierr)
       call mpi_comm_size(mpi_ggacmp_cross_world(j),newpes,ierr)
       call mpi_comm_rank(mpi_ggacmp_cross_world(j),newmype,ierr)
       if(ipri >= 2) then
! === DEBUG by tkato 2011/07/12 ================================================
!         write(nfout,'(" !|| mype, newmype, newpes, mpi_ggacmp_cross_world = ",4i8)') &
!              & mype, newmype,newpes,mpi_ggacmp_cross_world(j)
          write(nfout,'(" !|| mype, newmype, newpes, mpi_ggacmp_cross_world = ",4i12)') &
               & mype, newmype,newpes,mpi_ggacmp_cross_world(j)
! ==============================================================================
       end if
    end do
    firstcall = .false.
  end subroutine m_Parallel_init_mpi_gga


  subroutine m_Parallel_dealloc_mpi_mix()
    if(allocated(is_kgpm))  deallocate(is_kgpm)
    if(allocated(ie_kgpm))  deallocate(ie_kgpm)
    if(allocated(nel_kgpm)) deallocate(nel_kgpm)
  end subroutine m_Parallel_dealloc_mpi_mix

  subroutine m_Parallel_init_mpi_mix(nfout,ipri,printable,kgpm,comm_for_chg)
    integer, intent(in) :: nfout,ipri,kgpm
    logical, intent(in) :: printable
    logical, intent(in) :: comm_for_chg
    integer :: iwork, i
    integer :: npes_, mype_
                                                  __TIMER_SUB_START(1241)
    if(comm_for_chg)then
      npes_ = npes
      mype_ = mype
    else
      npes_ = nrank_g
      mype_ = myrank_g
    endif
    allocate(is_kgpm(0:npes_-1))
    allocate(ie_kgpm(0:npes_-1))
    allocate(nel_kgpm(0:npes_-1))
    iwork = ( kgpm - 1 ) / npes_ + 1
    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| << init_mpi_mix >>")')
       write(nfout,'(" !|| kgpm = ",i12)') kgpm
       write(nfout,'(" !|| -- is_kgpm, ie_kgpm --")')
    end if
    do i = 0, npes_-1
       is_kgpm(i) = min(i*iwork+1, kgpm+1)
       ie_kgpm(i) = min(is_kgpm(i)+iwork-1, kgpm)
       nel_kgpm(i) = ie_kgpm(i) - is_kgpm(i) + 1
       if(ipri >= 2 .and. printable) write(nfout,'(" !|| ",2i12)') is_kgpm(i),ie_kgpm(i)
    enddo
    ista_kgpm = is_kgpm(mype_)
    iend_kgpm = ie_kgpm(mype_)
    np_kgpm   = nel_kgpm(mype_)
    mp_kgpm   = maxval(nel_kgpm)
    if(ipri > 1 .and. printable) then
       if(iend_kgpm <1000000) then
          write(nfout,'(" !|| ista_kgpm, iend_kgpm = ",i6,",",i6,", mp_kgpm = ",i6,", kgpm = ",i12 &
               & ,"  << m_Parallel_init_mpi_mix>>")') ista_kgpm,iend_kgpm,mp_kgpm,kgpm
       elseif(iend_kgpm <100000000)then
          write(nfout,'(" !|| ista_kgpm, iend_kgpm = ",i8,",",i8,", mp_kgpm = ",i8,", kgpm = ",i12 &
               & ,"  << m_Parallel_init_mpi_mix>>")') ista_kgpm,iend_kgpm,mp_kgpm,kgpm
       else
          write(nfout,'(" !|| ista_kgpm, iend_kgpm = ",i12,",",i12,", mp_kgpm = ",i12,", kgpm = ",i0 &
               & ,"  << m_Parallel_init_mpi_mix>>")') ista_kgpm,iend_kgpm,mp_kgpm,kgpm
       end if
    end if
                                                  __TIMER_SUB_START(1241)
  end subroutine m_Parallel_init_mpi_mix



  subroutine m_Parallel_init_mpi_ffth(nfout,ipri,printable,nfft)
    integer, intent(in) :: nfout,ipri,nfft
    logical, intent(in) :: printable
    integer :: iwork, i, nffth
    nffth = nfft/2

    if ( .not. allocated(is_ffth) ) allocate(is_ffth(0:npes-1))
    if ( .not. allocated(ie_ffth) ) allocate(ie_ffth(0:npes-1))
    if ( .not. allocated(nel_ffth) ) allocate(nel_ffth(0:npes-1))

    iwork = ( nffth - 1 ) / npes + 1
    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| << init_mpi_ffth >>")')
       write(nfout,'(" !|| nffth = ",i12)') nffth
       write(nfout,'(" !|| -- is_ffth, ie_ffth --")')
    end if
    do i = 0, npes-1
       is_ffth(i) = min(i*iwork+1, nffth+1)
       ie_ffth(i) = min(is_ffth(i)+iwork-1, nffth)
       nel_ffth(i) = ie_ffth(i) - is_ffth(i) + 1
       if(ipri >= 2 .and. printable) write(nfout,'(" !|| ",2i12)') is_ffth(i),ie_ffth(i)
    enddo
    ista_ffth = is_ffth(mype)
    iend_ffth = ie_ffth(mype)
    np_ffth   = nel_ffth(mype)
    mp_ffth   = maxval(nel_ffth)
    if(ipri > 1 .and. printable) then
       if(iend_ffth <1000) then
          write(nfout,'(" !|| ista_ffth, iend_ffth = ",i3,",",i3,", mp_ffth = ",i3,", nffth = ",i5 &
               & ,"  << m_Parallel_init_mpi_ffth >>")') ista_ffth,iend_ffth,mp_ffth,nffth
       elseif(iend_ffth <100000)then
          write(nfout,'(" !|| ista_ffth, iend_ffth = ",i5,",",i5,", mp_ffth = ",i5,", ffth = ",i6 &
               & ,"  << m_Parallel_init_mpi_ffth >>")') ista_ffth,iend_ffth,mp_ffth,nffth
       else
          write(nfout,'(" !|| ista_ffth, iend_ffth = ",i0,",",i0,", mp_atm = ",i0,", natm = ",i0 &
               & ,"  << m_Parallel_init_mpi_ffth >>")') ista_ffth,iend_ffth,mp_ffth,nffth
       end if
    end if
  end subroutine m_Parallel_init_mpi_ffth

  subroutine m_Parallel_init_mpi_snl_3D(nfout,ipri,printable,nspin)
    integer, intent(in) :: nfout,ipri,nspin
    logical, intent(in) :: printable
                                                  __TIMER_SUB_START(1237)
    ista_snl = (ista_k + nspin - 1)/nspin
    iend_snl = iend_k/nspin
    if(ipri > 1 .and. printable) then
       write(nfout,'(" !|| ista_snl, iend_snl = ",i5,",",i5," << m_Parallel_init_mpi_snl_3D >>")') &
            & ista_snl,iend_snl
    end if
                                                  __TIMER_SUB_STOP(1237)
  end subroutine m_Parallel_init_mpi_snl_3D

  subroutine m_Parallel_init_mpi_nn(nfout,ipri,printable,nn)
    integer, intent(in) :: nfout,ipri,nn
    logical, intent(in) :: printable
    integer :: iwork, i,ia

    allocate(is_nn(0:nrank_e-1))
    allocate(ie_nn(0:nrank_e-1))
    allocate(nel_nn(0:nrank_e-1))
    iwork = ( nn - 1 ) / nrank_e + 1
    do i = 0, nrank_e-1
       is_nn(i) = min(i*iwork+1, nn+1)
       ie_nn(i) = min(is_nn(i)+iwork-1, nn)
       nel_nn(i) = ie_nn(i) - is_nn(i) + 1
       if(ipri >= 2 .and. printable) write(nfout,'(" !|| ",2i12)') is_nn(i),ie_nn(i)
    enddo
    ista_nn = is_nn(myrank_e)
    iend_nn = ie_nn(myrank_e)
                                                  __TIMER_SUB_STOP(1238)
  end subroutine m_Parallel_init_mpi_nn

  subroutine m_Parallel_init_mpi_rspace_aug(nfout,ipri,printable,natm,nmesh_rs_aug)
    integer, intent(in) :: nfout,ipri,natm
    logical, intent(in) :: printable
    integer, dimension(natm), intent(in) :: nmesh_rs_aug
    integer :: iwork, i,ia
    integer :: npes, mype
                                                  __TIMER_SUB_START(1238)
    npes = nrank_g
    mype = myrank_g

    if(allocated(is_rspace_aug))deallocate(is_rspace_aug)
    if(allocated(ie_rspace_aug))deallocate(ie_rspace_aug)
    if(allocated(nel_rspace_aug))deallocate(nel_rspace_aug)
    if(allocated(ista_rspace_aug_atm))deallocate(ista_rspace_aug_atm)
    if(allocated(iend_rspace_aug_atm))deallocate(iend_rspace_aug_atm)
    allocate(is_rspace_aug(0:npes-1))
    allocate(ie_rspace_aug(0:npes-1))
    allocate(nel_rspace_aug(0:npes-1))
    allocate(ista_rspace_aug_atm(natm));ista_rspace_aug_atm=0
    allocate(iend_rspace_aug_atm(natm));iend_rspace_aug_atm=0
    do ia=1,natm
    iwork = ( nmesh_rs_aug(ia) - 1 ) / npes + 1
    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| << init_mpi_rspace_aug >>")')
       write(nfout,'(" !|| nmesh_rs_aug = ",i12)') nmesh_rs_aug(ia)
       write(nfout,'(" !|| -- is_rspace_aug, ie_rspace_aug --")')
    end if
    do i = 0, npes-1
       is_rspace_aug(i) = min(i*iwork+1, nmesh_rs_aug(ia)+1)
       ie_rspace_aug(i) = min(is_rspace_aug(i)+iwork-1, nmesh_rs_aug(ia))
       nel_rspace_aug(i) = ie_rspace_aug(i) - is_rspace_aug(i) + 1
       if(ipri >= 2 .and. printable) write(nfout,'(" !|| ",2i12)') is_rspace_aug(i),ie_rspace_aug(i)
    enddo
    ista_rspace_aug_atm(ia) = is_rspace_aug(mype)
    iend_rspace_aug_atm(ia) = ie_rspace_aug(mype)
    if(ipri >= 2 .and. printable) then
       if(iend_rspace_aug_atm(ia) <1000) then
          write(nfout,'(" !|| ista_rspace_aug, iend_rspace_aug = ",i3,",",i3,", mp_rspace_aug = ",i3,", nmesh_rs_aug_max = ",i5 &
          & ,"  << m_Parallel_init_mpi_rspace_aug >>")') &
          & ista_rspace_aug_atm(ia),iend_rspace_aug_atm(ia),mp_rspace_aug,nmesh_rs_aug(ia)
       else
          write(nfout,'(" !|| ista_rspace_aug, iend_rspace_aug = ",i0,",",i0,", mp_rspace_aug = ",i0,", nmesh_rs_aug_max = ",i0 &
          & ,"  << m_Parallel_init_mpi_atm >>")') ista_rspace_aug_atm(ia),iend_rspace_aug_atm(ia),mp_rspace_aug,nmesh_rs_aug(ia)
       end if
    end if
    enddo
                                                  __TIMER_SUB_STOP(1238)
  end subroutine m_Parallel_init_mpi_rspace_aug

  subroutine m_Parallel_init_mpi_atm_ke(nfout,ipri,printable,natm)
    integer, intent(in) :: nfout,ipri,natm
    logical, intent(in) :: printable
    integer :: iwork, i
    integer :: npes, mype
                                                  __TIMER_SUB_START(1238)
    npes = nrank_ke
    mype = myrank_ke

    allocate(is_atm_ke(0:npes-1))
    allocate(ie_atm_ke(0:npes-1))
    allocate(nel_atm_ke(0:npes-1))
    iwork = ( natm - 1 ) / npes + 1
    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| << init_mpi_atm >>")')
       write(nfout,'(" !|| natm = ",i12)') natm
       write(nfout,'(" !|| -- is_natm, ie_natm --")')
    end if
    do i = 0, npes-1
       is_atm_ke(i) = min(i*iwork+1, natm+1)
       ie_atm_ke(i) = min(is_atm_ke(i)+iwork-1, natm)
       nel_atm_ke(i) = ie_atm_ke(i) - is_atm_ke(i) + 1
       if(ipri >= 2 .and. printable) write(nfout,'(" !|| ",2i12)') is_atm_ke(i),ie_atm_ke(i)
    enddo
    ista_atm_ke = is_atm_ke(mype)
    iend_atm_ke = ie_atm_ke(mype)
    np_atm_ke   = nel_atm_ke(mype)
    mp_atm_ke   = maxval(nel_atm_ke)
!!$    if(ipri == 1 .and. printable) then
    if(ipri>1) then
      write(nfout,*) ' ==== m_Parallel_init_mpi_atm_ke ===='
      write(nfout,*) 'npes : ', npes, ' = nrank_ke : ', nrank_ke
      write(nfout,*) ' ----- natm_ke Decomposed as ... ----'
      write(nfout,*) ' i, is_atm_ke, ie_atm_ke, nel_atm_ke'
      do i = 0, npes-1
         write(nfout,*) i, is_atm_ke(i), ie_atm_ke(i), nel_atm_ke(i)
      end do
      if(iend_atm <1000) then
         write(nfout,'(" !|| ista_atm_ke, iend_atm_ke = ",i3,",",i3,", mp_atm = ",i3,", natm = ",i5 &
              & ,"  << m_Parallel_init_mpi_atm >>")') ista_atm_ke,iend_atm_ke,mp_atm_ke,natm
      elseif(iend_kgpm <100000)then
         write(nfout,'(" !|| ista_atm_ke, iend_atm_ke = ",i5,",",i5,", mp_atm = ",i5,", natm = ",i6 &
              & ,"  << m_Parallel_init_mpi_atm >>")') ista_atm_ke,iend_atm_ke,mp_atm_ke,natm
      else
         write(nfout,'(" !|| ista_atm_ke, iend_atm_ke = ",i0,",",i0,", mp_atm = ",i0,", natm = ",i0 &
              & ,"  << m_Parallel_init_mpi_atm >>")') ista_atm_ke,iend_atm_ke,mp_atm_ke,natm
      end if
    end if
  end subroutine m_Parallel_init_mpi_atm_ke

  subroutine m_Parallel_init_mpi_nq(nfout,ipri,printable,nq)
    integer, intent(in) :: nfout, ipri,nq
    logical, intent(in) :: printable
    integer :: i,j,ip

    allocate(is_nq(0:nrank_ke-1));ista_nq=nq+1
    allocate(ie_nq(0:nrank_ke-1));iend_nq=0
    allocate(nel_nq(0:nrank_ke-1));nel_nq=0
    allocate(map_nq(nq));map_nq=0
    allocate(map_z_nq(nq));map_z_nq=0

    call set_block_range(nq,nrank_ke,nel_nq,is_nq,ie_nq,.true.,map_nq)
    ista_nq = is_nq(myrank_ke)
    iend_nq = ie_nq(myrank_ke)
    np_nq = nel_nq(myrank_ke)
    mp_nq = maxval(nel_nq)
    j = 0
    do ip = 1, nrank_ke
       do i = 1, nel_nq(ip-1)
          j = j + 1
          map_z_nq(j) = i
       end do
    end do
  end subroutine m_Parallel_init_mpi_nq

  subroutine m_Parallel_init_mpi_atm_f(nfout,ipri,printable,natm)
    integer, intent(in) :: nfout,ipri,natm
    logical, intent(in) :: printable
    integer :: iwork, i
    !!integer :: npes, mype
                                                  __TIMER_SUB_START(1238)
    !!npes = nrank_g
    !!mype = myrank_g

    if(allocated(is_atm_f)) deallocate(is_atm_f)
    if(allocated(ie_atm_f)) deallocate(ie_atm_f)
    if(allocated(nel_atm_f)) deallocate(nel_atm_f)
    allocate(is_atm_f(0:npes-1))
    allocate(ie_atm_f(0:npes-1))
    allocate(nel_atm_f(0:npes-1))
    iwork = ( natm - 1 ) / npes + 1
    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| << init_mpi_atm >>")')
       write(nfout,'(" !|| natm = ",i12)') natm
       write(nfout,'(" !|| -- is_natm, ie_natm --")')
    end if
    do i = 0, npes-1
       is_atm_f(i) = min(i*iwork+1, natm+1)
       ie_atm_f(i) = min(is_atm_f(i)+iwork-1, natm)
       nel_atm_f(i) = ie_atm_f(i) - is_atm_f(i) + 1
       if(ipri >= 2 .and. printable) write(nfout,'(" !|| ",2i12)') is_atm_f(i),ie_atm_f(i)
    enddo
    ista_atm_f = is_atm_f(mype)
    iend_atm_f = ie_atm_f(mype)
    np_atm_f   = nel_atm_f(mype)
    mp_atm_f   = maxval(nel_atm_f)
    if(ipri > 1 .and. printable) then
       if(iend_atm <1000) then
          write(nfout,'(" !|| ista_atm_f, iend_atm_f = ",i3,",",i3,", mp_atm = ",i3,", natm = ",i5 &
               & ,"  << m_Parallel_init_mpi_atm_f >>")') ista_atm_f,iend_atm_f,mp_atm_f,natm
       elseif(iend_kgpm <100000)then
          write(nfout,'(" !|| ista_atm_f, iend_atm_f = ",i5,",",i5,", mp_atm = ",i5,", natm = ",i6 &
               & ,"  << m_Parallel_init_mpi_atm_f >>")') ista_atm_f,iend_atm_f,mp_atm_f,natm
       else
          write(nfout,'(" !|| ista_atm_f, iend_atm_f = ",i0,",",i0,", mp_atm = ",i0,", natm = ",i0 &
               & ,"  << m_Parallel_init_mpi_atm_f >>")') ista_atm_f,iend_atm_f,mp_atm_f,natm
       end if
    end if
                                                  __TIMER_SUB_STOP(1238)
  end subroutine m_Parallel_init_mpi_atm_f

  subroutine m_Parallel_init_mpi_atm(nfout,ipri,printable,natm)
    integer, intent(in) :: nfout,ipri,natm
    logical, intent(in) :: printable
    integer :: iwork, i
    integer :: npes, mype
                                                  __TIMER_SUB_START(1238)
    npes = nrank_g
    mype = myrank_g

    allocate(is_atm(0:npes-1))
    allocate(ie_atm(0:npes-1))
    allocate(nel_atm(0:npes-1))
    iwork = ( natm - 1 ) / npes + 1
    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| << init_mpi_atm >>")')
       write(nfout,'(" !|| natm = ",i12)') natm
       write(nfout,'(" !|| -- is_natm, ie_natm --")')
    end if
    do i = 0, npes-1
       is_atm(i) = min(i*iwork+1, natm+1)
       ie_atm(i) = min(is_atm(i)+iwork-1, natm)
       nel_atm(i) = ie_atm(i) - is_atm(i) + 1
       if(ipri >= 2 .and. printable) write(nfout,'(" !|| ",2i12)') is_atm(i),ie_atm(i)
    enddo
    ista_atm = is_atm(mype)
    iend_atm = ie_atm(mype)
    np_atm   = nel_atm(mype)
    mp_atm   = maxval(nel_atm)
! ====================================================================
    if(ipri>1) then
       write(nfout,*) '===== m_Parallel_init_mpi_atm ====='
       write(nfout,*) 'npes: ', npes, ' = nrank_g : ',nrank_g
       write(nfout,*) 'mype: ', mype
       write(nfout,*) 'ista_atm, iend_atm, np_natm, mp_atm: ', ista_atm, iend_atm, np_atm, mp_atm
       write(nfout,*) '----- natm Decomposed as... -----'
       write(nfout,*) 'natm:  ', natm
       write(nfout,*) ' i, is_atm, ie_atm, nel_atm'
       do i = 0, npes-1
         write(nfout,*) i, is_atm(i), ie_atm(i), nel_atm(i)
       enddo
! =====================================================================
       if(iend_atm <1000) then
          write(nfout,'(" !|| ista_atm, iend_atm = ",i3,",",i3,", mp_atm = ",i3,", natm = ",i5 &
               & ,"  << m_Parallel_init_mpi_atm >>")') ista_atm,iend_atm,mp_atm,natm
       elseif(iend_kgpm <100000)then
          write(nfout,'(" !|| ista_atm, iend_atm = ",i5,",",i5,", mp_atm = ",i5,", natm = ",i6 &
               & ,"  << m_Parallel_init_mpi_atm >>")') ista_atm,iend_atm,mp_atm,natm
       else
          write(nfout,'(" !|| ista_atm, iend_atm = ",i0,",",i0,", mp_atm = ",i0,", natm = ",i0 &
               & ,"  << m_Parallel_init_mpi_atm >>")') ista_atm,iend_atm,mp_atm,natm
       end if
    end if
                                                  __TIMER_SUB_STOP(1238)
  end subroutine m_Parallel_init_mpi_atm

  subroutine m_Parallel_init_mpi_atm2(nfout,ipri,printable,natm2)
    integer, intent(in) :: nfout,ipri,natm2
    logical, intent(in) :: printable
    integer :: iwork, i
    integer :: npes, mype
                                                  __TIMER_SUB_START(1239)
    npes = nrank_g
    mype = myrank_g
    if(allocated(is_atm2)) deallocate(is_atm2)
    if(allocated(ie_atm2)) deallocate(ie_atm2)
    if(allocated(nel_atm2)) deallocate(nel_atm2)
    allocate(is_atm2(0:npes-1))
    allocate(ie_atm2(0:npes-1))
    allocate(nel_atm2(0:npes-1))
    iwork = ( natm2 - 1 ) / npes + 1
    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| << init_mpi_atm2 >>")')
       write(nfout,'(" !|| natm2 = ",i12)') natm2
       write(nfout,'(" !|| -- is_natm2, ie_natm2 --")')
    end if
    do i = 0, npes-1
       is_atm2(i) = min(i*iwork+1, natm2+1)
       ie_atm2(i) = min(is_atm2(i)+iwork-1, natm2)
       nel_atm2(i) = ie_atm2(i) - is_atm2(i) + 1
       if(ipri >= 2 .and. printable) write(nfout,'(" !|| ",2i12)') is_atm2(i),ie_atm2(i)
    enddo
    ista_atm2 = is_atm2(mype)
    iend_atm2 = ie_atm2(mype)
    np_atm2   = nel_atm2(mype)
    mp_atm2   = maxval(nel_atm2)

    if(ipri > 1 .and. printable) then
       if(iend_atm2 <1000) then
          write(nfout,'(" !|| ista_atm2, iend_atm2 = ",i3,",",i3,", mp_atm2 = ",i3,", natm2 = ",i5 &
               & ,"  << m_Parallel_init_mpi_atm2 >>")') ista_atm2,iend_atm2,mp_atm2,natm2
       elseif(iend_kgpm <100000)then
          write(nfout,'(" !|| ista_atm2, iend_atm2 = ",i5,",",i5,", mp_atm2 = ",i5,", natm2 = ",i6 &
               & ,"  << m_Parallel_init_mpi_atm2 >>")') ista_atm2,iend_atm2,mp_atm2,natm2
       else
          write(nfout,'(" !|| ista_atm2, iend_atm2 = ",i0,",",i0,", mp_atm2 = ",i0,", natm2 = ",i0 &
               & ,"  << m_Parallel_init_mpi_atm2 >>")') ista_atm2,iend_atm2,mp_atm2,natm2
       end if
    end if

                                                  __TIMER_SUB_STOP(1239)
  end subroutine m_Parallel_init_mpi_atm2


  subroutine split_into_ggablock_cdfft(npes,nrank_ggacmp,npes_cdfft,nrest_cdfft)
    integer, intent(in)  :: npes
    integer, intent(out) :: nrank_ggacmp, npes_cdfft, nrest_cdfft

#ifdef _NO_GGA_SPLIT_
    nrank_ggacmp = 1
    npes_cdfft  = npes
    nrest_cdfft = 0
#else
    if(npes <= 2) then
       nrank_ggacmp = 1
       npes_cdfft  = npes
       nrest_cdfft = 0
    else
       nrank_ggacmp = 3
       npes_cdfft  = npes/nrank_ggacmp
       nrest_cdfft = npes - npes_cdfft*nrank_ggacmp
    end if
#endif
  end subroutine split_into_ggablock_cdfft

  subroutine m_Parallel_dealloc_mpi_nlmta()
    if(allocated(nis_fs)) deallocate(nis_fs)
    if(allocated(nie_fs)) deallocate(nie_fs)
    if(allocated(nel_fs)) deallocate(nel_fs)
    if(allocated(nis_fs_atm)) deallocate(nis_fs_atm)
    if(allocated(nie_fs_atm)) deallocate(nie_fs_atm)
    if(allocated(nel_fs_atm)) deallocate(nel_fs_atm)
    if(allocated(ball_buff)) deallocate(ball_buff)
    if(allocated(ball_addr)) deallocate(ball_addr)
  end subroutine m_Parallel_dealloc_mpi_nlmta

  subroutine m_Parallel_dealloc_mpi_exx()
    if(allocated(is_kngp_exx)) deallocate(is_kngp_exx)
    if(allocated(ie_kngp_exx)) deallocate(ie_kngp_exx)
    if(allocated(nel_kngp_exx)) deallocate(nel_kngp_exx)
  end subroutine m_Parallel_dealloc_mpi_exx

  subroutine m_Parallel_cp_g1k(kv3)
    integer, intent(in) :: kv3
    if(.not.allocated(ista_g1k_prev)) allocate(ista_g1k_prev(kv3))
    if(.not.allocated(iend_g1k_prev)) allocate(iend_g1k_prev(kv3))
    ista_g1k_prev = ista_g1k
    iend_g1k_prev = iend_g1k
    ista_kngp_prev = ista_kngp
    iend_kngp_prev = iend_kngp
  end subroutine m_Parallel_cp_g1k

  subroutine m_Parallel_dealloc(neb_mode)
    logical, intent(in), optional :: neb_mode
    logical :: neb
    neb = .false.
    if(present(neb_mode)) neb = neb_mode

    if(allocated(is_kngp)) deallocate(is_kngp)
    if(allocated(ie_kngp)) deallocate(ie_kngp)
    if(allocated(nel_kngp)) deallocate(nel_kngp)

!    if(.not.neb)then
    if(allocated(map_e)) deallocate(map_e)
    if(allocated(map_z)) deallocate(map_z)
    if(allocated(map_k)) deallocate(map_k)
    if(allocated(map_ek)) deallocate(map_ek)
    if(allocated(map_s)) deallocate(map_s)
!    if(allocated(mpi_k_world)) deallocate(mpi_k_world)
    if(allocated(nis_e)) deallocate(nis_e)
    if(allocated(nie_e)) deallocate(nie_e)
    if(allocated(nel_e)) deallocate(nel_e)
    if(allocated(idisp_e)) deallocate(idisp_e)
#ifdef NEC_TUNE_SOFT
    if(allocated(ista_e_smp)) deallocate(ista_e_smp)
    if(allocated(iend_e_smp)) deallocate(iend_e_smp)
#endif
    if(allocated(nis_k)) deallocate(nis_k)
    if(allocated(nie_k)) deallocate(nie_k)
    if(allocated(nel_k)) deallocate(nel_k)

    if(allocated(ista_g1k)) deallocate(ista_g1k)
    if(allocated(iend_g1k)) deallocate(iend_g1k)
    if(allocated(np_g1k)) deallocate(np_g1k)
    if(allocated(np_g1k_prev)) deallocate(np_g1k_prev)
    if(allocated(mp_g1k)) deallocate(mp_g1k)

    if(allocated(nis_fftp)) deallocate(nis_fftp)
    if(allocated(nie_fftp)) deallocate(nie_fftp)
    if(allocated(nel_fftp)) deallocate(nel_fftp)
    if(allocated(idisp_fftp)) deallocate(idisp_fftp)
    if(allocated(nis_fftph)) deallocate(nis_fftph)
    if(allocated(nie_fftph)) deallocate(nie_fftph)
    if(allocated(nel_fftph)) deallocate(nel_fftph)

    if(allocated(is_kgpm)) deallocate(is_kgpm)
    if(allocated(ie_kgpm)) deallocate(ie_kgpm)
    if(allocated(nel_kgpm)) deallocate(nel_kgpm)

    if(allocated(is_atm)) deallocate(is_atm)
    if(allocated(ie_atm)) deallocate(ie_atm)
    if(allocated(nel_atm)) deallocate(nel_atm)

    if(allocated(is_atm2)) deallocate(is_atm2)
    if(allocated(ie_atm2)) deallocate(ie_atm2)
    if(allocated(nel_atm2)) deallocate(nel_atm2)
!!$#ifndef TRANSPOSE
!!$    endif
!!$#endif
!!$#ifdef TRANSPOSE
    if(allocated(nis_g1)) deallocate(nis_g1)
    if(allocated(nie_g1)) deallocate(nie_g1)
    if(allocated(nel_g1)) deallocate(nel_g1)
    if(.not.neb)then
    if(allocated(nis_fs)) deallocate(nis_fs)
    if(allocated(nie_fs)) deallocate(nie_fs)
    if(allocated(nel_fs)) deallocate(nel_fs)
    if(allocated(nis_fs_atm)) deallocate(nis_fs_atm)
    if(allocated(nie_fs_atm)) deallocate(nie_fs_atm)
    if(allocated(nel_fs_atm)) deallocate(nel_fs_atm)
    endif

    if(allocated(nis_g1k)) deallocate(nis_g1k)
    if(allocated(nie_g1k)) deallocate(nie_g1k)
    if(allocated(nel_g1k)) deallocate(nel_g1k)
!!$#endif

    if(allocated(map_pe2ggacmp)) deallocate(map_pe2ggacmp)
    if(allocated(map_pe2cdfft)) deallocate(map_pe2cdfft)

!    if(.not.neb)then
    if(allocated(map_ggacmp)) deallocate(map_ggacmp)
    if(allocated(idisp_fftph)) deallocate(idisp_fftph)
!    if(allocated(mpi_cdfft_world)) deallocate(mpi_cdfft_world)
    !if(allocated(mpi_ggacmp_cross_world)) deallocate(mpi_ggacmp_cross_world)

    if(allocated(nis_sfftp)) deallocate(nis_sfftp)
    if(allocated(nie_sfftp)) deallocate(nie_sfftp)
    if(allocated(nel_sfftp)) deallocate(nel_sfftp)
    if(allocated(idisp_sfftp)) deallocate(idisp_sfftp)
    if(allocated(nis_sfftph)) deallocate(nis_sfftph)
    if(allocated(nie_sfftph)) deallocate(nie_sfftph)
    if(allocated(nel_sfftph)) deallocate(nel_sfftph)

!    if(allocated(mpi_nbmx_world)) deallocate(mpi_nbmx_world)
    if(allocated(nel_nbmx)) deallocate(nel_nbmx)
    if(allocated(nis_nbmx)) deallocate(nis_nbmx)
    if(allocated(nie_nbmx)) deallocate(nie_nbmx)
    if(allocated(idisp_nbmx)) deallocate(idisp_nbmx)

    if(allocated(nel_nbmx_k)) deallocate(nel_nbmx_k)
    if(allocated(nis_nbmx_k)) deallocate(nis_nbmx_k)
    if(allocated(nie_nbmx_k)) deallocate(nie_nbmx_k)
    if(allocated(idisp_nbmx_k)) deallocate(idisp_nbmx_k)
!    if(allocated(mpi_nbmx_world_k)) deallocate(mpi_nbmx_world_k)
    if(allocated(nel_kg1_k)) deallocate(nel_kg1_k)
    if(allocated(nis_kg1_k)) deallocate(nis_kg1_k)
    if(allocated(nie_kg1_k)) deallocate(nie_kg1_k)
    if(allocated(idisp_kg1_k)) deallocate(idisp_kg1_k)
!    endif
    if(allocated(wf_fft_scnt)) deallocate(wf_fft_scnt)
    if(allocated(wf_fft_rcnt)) deallocate(wf_fft_rcnt)
    if(allocated(wf_fft_index)) deallocate(wf_fft_index)
    if(allocated(wf_fft_dist)) deallocate(wf_fft_dist)
    if(allocated(wf_fft_send)) deallocate(wf_fft_send)
    if(allocated(wf_fft_recv)) deallocate(wf_fft_recv)
    if(allocated(wf_fft_maxsend)) deallocate(wf_fft_maxsend)
    if(allocated(wf_fft_maxrecv)) deallocate(wf_fft_maxrecv)

    if(allocated(fft_wf_scnt)) deallocate(fft_wf_scnt)
    if(allocated(fft_wf_rcnt)) deallocate(fft_wf_rcnt)
    if(allocated(fft_wf_index)) deallocate(fft_wf_index)
    if(allocated(fft_wf_dist)) deallocate(fft_wf_dist)
    if(allocated(fft_wf_send)) deallocate(fft_wf_send)
    if(allocated(fft_wf_recv)) deallocate(fft_wf_recv)
    if(allocated(fft_wf_maxsend)) deallocate(fft_wf_maxsend)
    if(allocated(fft_wf_maxrecv)) deallocate(fft_wf_maxrecv)

    if(allocated(fft_chgq_scnt)) deallocate(fft_chgq_scnt)
    if(allocated(fft_chgq_rcnt)) deallocate(fft_chgq_rcnt)
    if(allocated(fft_chgq_index)) deallocate(fft_chgq_index)
    if(allocated(fft_chgq_dist)) deallocate(fft_chgq_dist)
    if(allocated(fft_chgq_send)) deallocate(fft_chgq_send)
    if(allocated(fft_chgq_recv)) deallocate(fft_chgq_recv)

#ifdef MPI_FFTW
    if(allocated(fft_chgq_scnt_mpifftw)) deallocate(fft_chgq_scnt_mpifftw)
    if(allocated(fft_chgq_rcnt_mpifftw)) deallocate(fft_chgq_rcnt_mpifftw)
    if(allocated(fft_chgq_index_mpifftw)) deallocate(fft_chgq_index_mpifftw)
    if(allocated(fft_chgq_dist_mpifftw)) deallocate(fft_chgq_dist_mpifftw)
    if(allocated(fft_chgq_send_mpifftw)) deallocate(fft_chgq_send_mpifftw)
    if(allocated(fft_chgq_recv_mpifftw)) deallocate(fft_chgq_recv_mpifftw)
#endif

    if(allocated(chgq_fftcd_scnt)) deallocate(chgq_fftcd_scnt)
    if(allocated(chgq_fftcd_rcnt)) deallocate(chgq_fftcd_rcnt)
    if(allocated(chgq_fftcd_index)) deallocate(chgq_fftcd_index)
    if(allocated(chgq_fftcd_dist)) deallocate(chgq_fftcd_dist)
    if(allocated(chgq_fftcd_send)) deallocate(chgq_fftcd_send)
    if(allocated(chgq_fftcd_recv)) deallocate(chgq_fftcd_recv)

    if(allocated(fftcd_chgq_scnt)) deallocate(fftcd_chgq_scnt)
    if(allocated(fftcd_chgq_rcnt)) deallocate(fftcd_chgq_rcnt)
    if(allocated(fftcd_chgq_index)) deallocate(fftcd_chgq_index)
    if(allocated(fftcd_chgq_dist)) deallocate(fftcd_chgq_dist)
    if(allocated(fftcd_chgq_send)) deallocate(fftcd_chgq_send)
    if(allocated(fftcd_chgq_recv)) deallocate(fftcd_chgq_recv)

    if(allocated(is_atm_B)) deallocate(is_atm_B)
    if(allocated(ie_atm_B)) deallocate(ie_atm_B)
    if(allocated(nel_atm_B)) deallocate(nel_atm_B)
    if(allocated(mem_atm_B)) deallocate(mem_atm_B)

    if(allocated(is_natm)) deallocate(is_natm)
    if(allocated(ie_natm)) deallocate(ie_natm)
    if(allocated(nel_natm)) deallocate(nel_natm)
    if(allocated(is_nrc)) deallocate(is_nrc)
    if(allocated(ie_nrc)) deallocate(ie_nrc)
    if(allocated(nel_nrc)) deallocate(nel_nrc)

    if(allocated(is_atm_ke)) deallocate(is_atm_ke)
    if(allocated(ie_atm_ke)) deallocate(ie_atm_ke)
    if(allocated(nel_atm_ke)) deallocate(nel_atm_ke)

    if(.not.neb)then
    if(allocated(ball_buff)) deallocate(ball_buff)
    if(allocated(ball_addr)) deallocate(ball_addr)
    endif

#ifdef MPI_FFTW
    if(allocated(wf_fft_scnt_mfftw)) deallocate(wf_fft_scnt_mfftw)
    if(allocated(wf_fft_rcnt_mfftw)) deallocate(wf_fft_rcnt_mfftw)
    if(allocated(wf_fft_index_mfftw)) deallocate(wf_fft_index_mfftw)
    if(allocated(wf_fft_dist_mfftw)) deallocate(wf_fft_dist_mfftw)
    if(allocated(wf_fft_send_mfftw)) deallocate(wf_fft_send_mfftw)
    if(allocated(wf_fft_recv_mfftw)) deallocate(wf_fft_recv_mfftw)
    if(allocated(wf_fft_maxsend_mfftw)) deallocate(wf_fft_maxsend_mfftw)
    if(allocated(wf_fft_maxrecv_mfftw)) deallocate(wf_fft_maxrecv_mfftw)

    if(allocated(fft_wf_scnt_mfftw)) deallocate(fft_wf_scnt_mfftw)
    if(allocated(fft_wf_rcnt_mfftw)) deallocate(fft_wf_rcnt_mfftw)
    if(allocated(fft_wf_index_mfftw)) deallocate(fft_wf_index_mfftw)
    if(allocated(fft_wf_dist_mfftw)) deallocate(fft_wf_dist_mfftw)
    if(allocated(fft_wf_send_mfftw)) deallocate(fft_wf_send_mfftw)
    if(allocated(fft_wf_recv_mfftw)) deallocate(fft_wf_recv_mfftw)
    if(allocated(fft_wf_maxsend_mfftw)) deallocate(fft_wf_maxsend_mfftw)
    if(allocated(fft_wf_maxrecv_mfftw)) deallocate(fft_wf_maxrecv_mfftw)
#endif
  end subroutine m_Parallel_dealloc

 logical function m_Parallel_resolve_conf_para()
   integer :: i,nn
   integer :: n1
   integer iargc,narg
   character(100) :: q1
   character(100) arg
   nrank_conf=-1
   m_Parallel_resolve_conf_para = .false.
   narg = iargc()
   n1 = -1
   do nn=1,narg
      call getarg(nn,q1)

!!!!!!!!!!!!!!!!! modified by mizouchi@adv 2003.03.20 !!!!!!!!!!!!!!
      q1 = trim(adjustl(q1))
!!!!!!!!!!!!!!!!! modified by mizouchi@adv 2003.03.20 !!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!! modified by mizouchi@adv 2003.03.20 !!!!!!!!!!!!!!
!!          if(index(trim(adjustl(q1)),"ne").ne.0) then

      if(q1(1:1) == "n".and. (q1(2:2) == "r".or.q1(2:2) == "c" .or. q1(2:2) == "i")) then
!!!!!!!!!!!!!!!!! modified by mizouchi@adv 2003.03.20 !!!!!!!!!!!!!!
         if((q1(3:3) == "=".or.q1(3:3) == ":") .and. len_trim(q1(4:)).ne.0) then
             q1 = q1(4:)
             if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                read(q1,*) n1
             else
                call phase_error_with_msg(6,"wrong nr ",__LINE__,__FILE__)
             end if
         else if(q1(3:3) == "=".or. q1(3:3) == ":") then
             call getarg(nn+1,q1)
             if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                read(q1,*) n1
             else
                call phase_error_with_msg(6,"wrong nr ",__LINE__,__FILE__)
             end if
         else
              call getarg(nn+1,q1)
              if(q1(1:1) == ":".and. len_trim(q1(2:)).ne.0) then
                  q1 = q1(2:)
                  if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                     read(q1,*) n1
                  else
                     call phase_error_with_msg(6,"wrong nr ",__LINE__,__FILE__)
                  end if
              else if((q1(1:1) == "=".or.q1(1:1) == ":").and. len_trim(q1(2:)).eq.0) then
                  call getarg(nn+2,q1)
                  if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                     read(q1,*) n1
                  else
                     call phase_error_with_msg(6,"wrong nr ",__LINE__,__FILE__)
                  end if
               else
                  call phase_error_with_msg(6,"wrong nr ",__LINE__,__FILE__)
               end if
         end if
      end if
   enddo

   nrank_conf = n1

   if(nrank_conf>=1) m_Parallel_resolve_conf_para = .true.
   conf_para = m_Parallel_resolve_conf_para
   if(conf_para)then
     call config_parallel_init()
   endif

   contains

   subroutine config_parallel_init()
     integer :: i,j
     integer npestmp,mypetmp,icolor,ikey,mpi_err
     integer mype_e, npes_e
     integer, allocatable ::  new_comm_world(:)
     logical, save        ::  firstcall=.true.
     call mpi_init(mpi_err)
     call mpi_comm_size(mpi_comm_world,npestmp,mpi_err)
     call mpi_comm_rank(mpi_comm_world,mypetmp,mpi_err)
     allocate(new_comm_world(0:nrank_conf-1))
     do i=0,nrank_conf-1
       icolor = 0
       ikey = 0
       do j=0, npestmp/nrank_conf-1
         if(mypetmp == j + (npestmp/nrank_conf)*i ) then
           icolor = 1
           ikey = i
         end if
       end do
!       if(firstcall) call mpi_comm_split( mpi_comm_world, icolor, ikey, new_comm_world(i), mpi_err )
       if(.not.firstcall) call mpi_comm_free(new_comm_world(i), mpi_err)
       call mpi_comm_split( mpi_comm_world, icolor, ikey, new_comm_world(i), mpi_err )
     end do

     mype_conf = mypetmp/(npestmp/nrank_conf)
     MPI_CommGroup = new_comm_world(mype_conf)
     firstcall = .false.
   end subroutine config_parallel_init

 end function m_Parallel_resolve_conf_para

  subroutine m_Parallel_init_mpi_nvale(nfout,ipri,printable,nvale)
    integer, intent(in) :: nfout, ipri,nvale
    logical, intent(in) :: printable
    integer :: i,j,ip
    nrank_nvale = nrank_e
    myrank_nvale = myrank_e
    if(mpi_nvale_enabled) call m_Parallel_dealloc_mpi_nvale()
    call m_Parallel_alloc_mpi_nvale(nvale)
    call set_block_range(nvale,nrank_nvale,nel_nvale,nis_nvale,nie_nvale,.true.,map_nvale)
    mpi_nvale_enabled = .true.
    if(ipri>=2) call wd_val_range_3D()
    ista_nvale = nis_nvale(myrank_nvale)
    iend_nvale = nie_nvale(myrank_nvale)
    np_nvale = nel_nvale(myrank_nvale)
    mp_nvale = maxval(nel_nvale)

    j = 0
    do ip = 1, nrank_nvale
       do i = 1, nel_nvale(ip-1)
          j = j + 1
          map_z_nvale(j) = i
       end do
    end do

    contains

    subroutine wd_val_range_3D()
      if(printable) then
         write(nfout,'(" !|| nrank_nvale")')
         write(nfout,'(" !||    i    : ",20i4)')(i,i=0,nrank_nvale-1)
         write(nfout,'(" !|| ista_nvale   : ",20i4)')(nis_nvale(i),i=0,nrank_nvale-1)
         write(nfout,'(" !|| iend_nvale   : ",20i4)')(nie_nvale(i),i=0,nrank_nvale-1)
         write(nfout,'(" !|| nel_nvale    : ",20i4)')(nel_nvale(i),i=0,nrank_nvale-1)
      end if
    end subroutine wd_val_range_3D
  end subroutine m_Parallel_init_mpi_nvale

  subroutine m_Parallel_alloc_mpi_nvale(nvale)
    integer, intent(in) :: nvale
    allocate(nis_nvale(0:nrank_nvale-1));ista_nvale=nvale+1
    allocate(nie_nvale(0:nrank_nvale-1));iend_nvale=0
    allocate(nel_nvale(0:nrank_nvale-1));nel_nvale=0
    allocate(map_nvale(nvale));map_nvale=0
    allocate(map_z_nvale(nvale));map_z_nvale=0
  end subroutine m_Parallel_alloc_mpi_nvale

  subroutine m_Parallel_dealloc_mpi_nvale()
    if(allocated(nis_nvale)) deallocate(nis_nvale)
    if(allocated(nie_nvale)) deallocate(nie_nvale)
    if(allocated(nel_nvale)) deallocate(nel_nvale)
    if(allocated(map_nvale)) deallocate(map_nvale)
    if(allocated(map_z_nvale)) deallocate(map_z_nvale)
    mpi_nvale_enabled = .false.
  end subroutine m_Parallel_dealloc_mpi_nvale

  subroutine m_Parallel_init_mpi_kv3_ek(nfout,ipri,printable,kv3_ek,nspin)
    integer, intent(in) :: nfout, ipri,kv3_ek,nspin
    logical, intent(in) :: printable
    integer :: i,j,ip
    integer :: siz
    siz = kv3_ek/nspin
    nrank_kv3_ek = nrank_k
    myrank_kv3_ek = myrank_k
    call m_Parallel_alloc_mpi_kv3_ek(siz)
    call set_block_range(siz,nrank_kv3_ek,nel_kv3_ek,nis_kv3_ek,nie_kv3_ek,.true.,map_kv3_ek)
    ista_kv3_ek = nis_kv3_ek(myrank_kv3_ek)
    iend_kv3_ek = nie_kv3_ek(myrank_kv3_ek)
    np_kv3_ek = nel_kv3_ek(myrank_kv3_ek)
    mp_kv3_ek = maxval(nel_kv3_ek)

    j = 0
    do ip = 1, nrank_kv3_ek
       do i = 1, nel_kv3_ek(ip-1)
          j = j + 1
          map_z_kv3_ek(j) = i
       end do
    end do
  end subroutine m_Parallel_init_mpi_kv3_ek

  subroutine m_Parallel_alloc_mpi_kv3_ek(kv3_ek)
    integer, intent(in) :: kv3_ek
    allocate(nis_kv3_ek(0:nrank_kv3_ek-1));ista_kv3_ek=kv3_ek+1
    allocate(nie_kv3_ek(0:nrank_kv3_ek-1));iend_kv3_ek=0
    allocate(nel_kv3_ek(0:nrank_kv3_ek-1));nel_kv3_ek=0
    allocate(map_kv3_ek(kv3_ek));map_kv3_ek=0
    allocate(map_z_kv3_ek(kv3_ek));map_z_kv3_ek=0
  end subroutine m_Parallel_alloc_mpi_kv3_ek

!add fujitsu
!===============================================================================

  subroutine make_index_band_3D(nfout,ipri,printable,kv3,neg, nblocksize_mgs,nblocksize_mgs_is_given,nblocksize_mgs_default)
    integer, intent(in) :: nfout, ipri
    logical, intent(in) :: printable
    integer :: NB, neg, kv3
    integer :: i, j, k, neg_tmp

    integer            :: nblocksize_mgs
    integer            :: nblocksize_mgs_default
    logical            :: nblocksize_mgs_is_given
    integer, allocatable, dimension(:)   :: wk_sta, wk_sta1, wk_end, wk_end1, ir_buff
! === Change nblocksize_mgs if it's too large. by tkato 2014/===================
    integer            :: max_block_size
! ==============================================================================
                                                  __TIMER_SUB_START(1230)

    if(nblocksize_mgs_is_given) then
       NB = nblocksize_mgs
    else
       NB = nblocksize_mgs_default
    end if
    if (.not.nblocksize_mgs_is_given.and. &
    &   mod(int(real(neg)/real(nrank_e)),NB)/=0.and. &
!     if (mod(int(real(neg)/real(nrank_e)),NB)/=0.and. &
    &   mod(neg,nrank_e)==0) then
       NB = int(real(neg)/real(nrank_e))
       if(NB<1) NB=1
       nblocksize_mgs = NB
       nblocksize_mgs_is_given = .true.
    endif
! === Change nblocksize_mgs if it's too large. by tkato 2014/===================
    if(nrank_e /= 1) then
       max_block_size = int(real(neg)/real(nrank_e))
       if(NB > max_block_size) then
          if(mype == 0) then
             write(0,'(a)') '=== WARNING!!! =============================================='
             write(0,'(a)') 'Block size for block-cyclic division on band is too large!'
             write(0,'(a,i5,a)') 'Block size should not be greater than ', max_block_size, '!'
             write(0,'(a,i5,a)') 'So, nblocksize_mgs is changed into ', max_block_size, '!'
             write(0,'(a)') 'FYI: '
             write(0,'(a,i8)') '   neg: ', neg
             write(0,'(a,i5)') '   Specified block Size: ', NB
             write(0,'(a,i5)') '   Num. of procs. for band: ', nrank_e
             write(0,'(a,i5,i5)') '   nblocksize_mgs, nblocksize_mgs_default: ', &
                nblocksize_mgs, nblocksize_mgs_default
             write(0,'(a)') '=== WARNING!!! =============================================='
          end if
          NB = max_block_size
          nblocksize_mgs = max_block_size
          nblocksize_mgs_is_given = .true.
       end if
    end if
! ==============================================================================

    nbs_num  = (neg - 1)/NB + 1
    nbsn_num = (nbs_num - 1)/nrank_e + 1
    if(allocated(lrank))     deallocate(lrank)
    if(allocated(nbsn))      deallocate(nbsn)
    if(allocated(nbsn_sta))  deallocate(nbsn_sta)
    if(allocated(nbsn_end))  deallocate(nbsn_end)
    if(allocated(nbs_sta))   deallocate(nbs_sta)
    if(allocated(nbs_end))   deallocate(nbs_end)
    if(allocated(neg_g))     deallocate(neg_g)
!    if(allocated(neg_g_all)) deallocate(neg_g_all)

    allocate( lrank(nbs_num) )    ; lrank = 0
    allocate( nbsn(nbs_num) )     ; nbsn  = 0
    allocate( nbsn_sta(nbsn_num) )  ; nbsn_sta = 0
    allocate( nbsn_end(nbsn_num) )  ; nbsn_end = 0
    allocate( nbs_sta(nbs_num) )  ; nbs_sta = 0
    allocate( nbs_end(nbs_num) )  ; nbs_end = 0
    allocate( neg_g(NB*nbsn_num) ); neg_g = 0
    if(.not.allocated(neg_g_all)) allocate( neg_g_all(neg) ); neg_g_all = 0

    allocate( wk_sta1(nbs_num) )  ; wk_sta1 = 0
    allocate( wk_end1(nbs_num) )  ; wk_end1 = 0
    allocate( wk_sta(0:nrank_e) )  ; wk_sta = 0
    allocate( wk_end(0:nrank_e) )  ; wk_end = 0

    do i = 1,neg,NB
      nbs        = (i - 1)/NB + 1
      lrank(nbs) = mod(nbs-1,nrank_e)
      nbsn(nbs)  = (nbs - 1)/nrank_e + 1
      if( myrank_e == lrank(nbs)) then
        j         = nbsn(nbs)
        nbsn_sta(j) = (j-1)*NB+1
        nbsn_end(j) = j*NB
        if ( i+NB-1 > neg ) nbsn_end(j) = nbsn_sta(j) + mod(neg,NB) - 1
        neg_tmp   = i
        do k=nbsn_sta(j),nbsn_end(j)
          neg_g(k) = neg_tmp
          neg_tmp  = neg_tmp + 1
        end do
      end if
    end do
    np_e = maxval(nbsn_end)

  allocate(ir_buff(0:nrank_e-1))
  call mpi_allgather(np_e,1,mpi_integer,ir_buff(0),1,mpi_integer,mpi_kg_world,ierr)
! === DEBUG by tkato 2011/09/23 ================================================
  do i = 0, nrank_e-1
     if(ir_buff(i) == 0) then
        if(printable) write(nfout,'(" !|| ",a,i3,a)') 'np_e = 0 at nrank_e = ', i, '!!!'
        call flush(nfout)
        call mpi_barrier(MPI_CommGroup, ierr)
        call mpi_abort(MPI_CommGroup, -1, ierr)
        call mpi_finalize(ierr)
        call phase_error_with_msg(nfout,'nrank_e is inconsistent',__LINE__,__FILE__)
     endif
  enddo
! ==============================================================================
  if(ipri >=2 .and. printable) then
     write(nfout,'(" !|| np_e = ",i8)') np_e
     write(nfout,'(" !|| ir_buff = ",8i5)') ir_buff(0:nrank_e-1)
  end if

  wk_sta(0) = 1
  do i = 1,nrank_e-1
     wk_sta(i) = wk_sta(i-1) + ir_buff(i-1)
  enddo

    do i = 1,nbs_num
      if( myrank_e == lrank(i)) then
        wk_sta1(i) = wk_sta( lrank(i) ) +nbsn_sta(nbsn(i)) -1
!        wk_end1(i) = min(neg,wk_sta1(i) + NB-1)
        wk_end1(i) = min(neg,wk_sta1(i) + nbsn_end(nbsn(i)) - nbsn_sta(nbsn(i)) )
      end if
    end do
   call mpi_allreduce(wk_sta1,nbs_sta,nbs_num,mpi_integer,mpi_sum, mpi_kg_world,ierr)
   call mpi_allreduce(wk_end1,nbs_end,nbs_num,mpi_integer,mpi_sum, mpi_kg_world,ierr)

!!!!
  nel_e = ir_buff
  nis_e(0) = 1
  nie_e(0) = nel_e(0)
  do i = 1,nrank_e-1
     nis_e(i) = nis_e(i-1) + nel_e(i-1)
     nie_e(i) = nie_e(i-1) + nel_e(i)
  enddo
  ista_e = nis_e(myrank_e)
  iend_e = nie_e(myrank_e)
  mp_e = maxval(nel_e)
! wk_sta(:) = nis_e(:) - 1
  wk_sta(0:nrank_e-1) = nis_e(0:nrank_e-1) - 1
  call mpi_allgatherv(neg_g, nel_e(myrank_e), mpi_integer, &
 &                    neg_g_all, nel_e, wk_sta, mpi_integer, mpi_kg_world, ierr)

  do i = 0, nrank_e - 1
    k = 0
    do j = nis_e(i), nie_e(i)
      k = k + 1
      map_e(neg_g_all(j)) = i
      map_z(neg_g_all(j)) = k
    enddo
  enddo

    do j = 1, kv3
       do i = 1, neg
          map_ek(i,j) = map_e(i) + map_k(j)*(nrank_e)
       end do
    end do

    if(ipri>=2 .and. printable) then
       write(nfout,'(" !|| <<make_index_band_3D>>")')
       write(nfout,'(" !|| neg, NB, nbs_num, nbsn_num = ",4i8)') neg, NB, nbs_num, nbsn_num
       write(nfout,'(" !||    nbs, lrank,  nbsn, nbsn_sta, nbsn_end")')
       do i = 1, nbsn_num
          write(nfout,'(" !||", 3i7,2i10)') i, lrank(i),nbsn(i),nbsn_sta(i),nbsn_end(i)
       end do
    end if
!!!!
    deallocate( wk_sta1 )
    deallocate( wk_end1 )
    deallocate( wk_sta )
    deallocate( wk_end )
    deallocate( ir_buff )

                                                  __TIMER_SUB_STOP(1230)
#ifdef _FJ_DBG0_
     character(len=3) :: p
     character(len=15) :: fdbg1002
     character(len=15) :: fdbg1003
     character(len=15) :: fdbg1004
     character(len=15) :: fdbg1005
     write(p,'(i3.3)') mype
     fdbg1002='dbginfo1002.'//trim(p)
     fdbg1003='dbginfo1003.'//trim(p)
     fdbg1004='dbginfo1004.'//trim(p)
     fdbg1005='dbginfo1005.'//trim(p)
     open(1002, file=trim(fdbg1002), form='formatted')
!     open(1003, file=trim(fdbg1003), form='formatted')
!     open(1004, file=trim(fdbg1004), form='formatted')
!     open(1005, file=trim(fdbg1005), form='formatted')
!    write(1002,'("neg_g_all")')
!    write(1002,'(8i7)')(neg_g_all(i),i=1,neg)
#endif
  end subroutine make_index_band_3D

!===============================================================================

  subroutine make_index_band_for_Gdiv_3D(neg, nblocksize_mgs,nblocksize_mgs_is_given,nblocksize_mgs_default)
    integer :: NB, neg
    integer :: i, j, k, neg_tmp

    integer            :: nblocksize_mgs
    integer            :: nblocksize_mgs_default
    logical            :: nblocksize_mgs_is_given

    integer            :: nbs_num, nbsn_num, ista_eg, iend_eg, np_eg
    integer, allocatable, dimension(:) :: lrank, nbsn, nbsn_sta, nbsn_end, nbs_sta, nbs_end
    integer, allocatable, dimension(:) :: wk_sta, wk_sta1, wk_end, wk_end1, ir_buff
                                                  __TIMER_SUB_START(1231)

    if(nblocksize_mgs_is_given) then
       NB = nblocksize_mgs
    else
       NB = nblocksize_mgs_default
    end if

    nbs_num  = (neg - 1)/NB + 1
!   nbsn_num = (nbs_num - 1)/nrank_e + 1
    nbsn_num = (nbs_num - 1)/nrank_g + 1
    if(allocated(neg_gg))     deallocate(neg_gg)
!    if(allocated(neg_gg_all)) deallocate(neg_gg_all)
!    if(allocated(nis_eg))     deallocate(nis_eg)
!    if(allocated(nie_eg))     deallocate(nie_eg)
!    if(allocated(nel_eg))     deallocate(nel_eg)
    allocate( neg_gg(NB*nbsn_num) ); neg_gg = 0
    if(.not.allocated(neg_gg_all)) allocate( neg_gg_all(neg) ); neg_gg_all = 0
    if(.not.allocated(nis_eg)) allocate( nis_eg(0:nrank_g-1)); nis_eg = 0
    if(.not.allocated(nie_eg)) allocate( nie_eg(0:nrank_g-1)); nie_eg = 0
    if(.not.allocated(nel_eg)) allocate( nel_eg(0:nrank_g-1)); nel_eg = 0

    allocate( lrank(nbs_num) )    ; lrank = 0
    allocate( nbsn(nbs_num) )     ; nbsn  = 0
    allocate( nbsn_sta(nbsn_num) )  ; nbsn_sta = 0
    allocate( nbsn_end(nbsn_num) )  ; nbsn_end = 0
    allocate( nbs_sta(nbs_num) )  ; nbs_sta = 0
    allocate( nbs_end(nbs_num) )  ; nbs_end = 0

    allocate( wk_sta1(nbs_num) )  ; wk_sta1 = 0
    allocate( wk_end1(nbs_num) )  ; wk_end1 = 0
    allocate( wk_sta(0:nrank_g) )  ; wk_sta = 0
    allocate( wk_end(0:nrank_g) )  ; wk_end = 0

    do i = 1,neg,NB
      nbs        = (i - 1)/NB + 1
      lrank(nbs) = mod(nbs-1,nrank_g)
      nbsn(nbs)  = (nbs - 1)/nrank_g + 1
      if( myrank_g == lrank(nbs)) then
        j         = nbsn(nbs)
        nbsn_sta(j) = (j-1)*NB+1
        nbsn_end(j) = j*NB
        if ( i+NB-1 > neg ) nbsn_end(j) = nbsn_sta(j) + mod(neg,NB) - 1
        neg_tmp   = i
        do k=nbsn_sta(j),nbsn_end(j)
          neg_gg(k) = neg_tmp
          neg_tmp  = neg_tmp + 1
        end do
      end if
    end do
    np_eg = maxval(nbsn_end)

  allocate(ir_buff(0:nrank_g-1))
  call mpi_allgather(np_eg,1,mpi_integer,ir_buff(0),1,mpi_integer,mpi_ke_world,ierr)
  wk_sta(0) = 1
  do i = 1,nrank_g-1
     wk_sta(i) = wk_sta(i-1) + ir_buff(i-1)
  enddo

    do i = 1,nbs_num
      if( myrank_g == lrank(i)) then
        wk_sta1(i) = wk_sta( lrank(i) ) +nbsn_sta(nbsn(i)) -1
!        wk_end1(i) = min(neg,wk_sta1(i) + NB-1)
        wk_end1(i) = min(neg,wk_sta1(i) + nbsn_end(nbsn(i)) - nbsn_sta(nbsn(i)) )
      end if
    end do

   call mpi_allreduce(wk_sta1,nbs_sta,nbs_num,mpi_integer,mpi_sum, mpi_ke_world,ierr)
   call mpi_allreduce(wk_end1,nbs_end,nbs_num,mpi_integer,mpi_sum, mpi_ke_world,ierr)

!!!!
  nel_eg = ir_buff
  nis_eg(0) = 1
  nie_eg(0) = nel_eg(0)
  do i = 1,nrank_g-1
     nis_eg(i) = nis_eg(i-1) + nel_eg(i-1)
     nie_eg(i) = nie_eg(i-1) + nel_eg(i)
  enddo
  ista_eg = nis_eg(myrank_g)
  iend_eg = nie_eg(myrank_g)

! wk_sta(:) = nis_eg(:) - 1
  wk_sta(0:nrank_g-1) = nis_eg(0:nrank_g-1) - 1
  call mpi_allgatherv(neg_gg, nel_eg(myrank_g), mpi_integer, &
 &                    neg_gg_all, nel_eg, wk_sta, mpi_integer, mpi_ke_world, ierr)

! do i = 0, nrank_g - 1
!   k = 0
!   do j = nis_eg(i), nie_eg(i)
!     k = k + 1
!     map_eg(neg_gg_all(j)) = i
!     map_zg_3D(neg_gg_all(j)) = k
!   enddo
! enddo
!!!!
    deallocate( lrank )
    deallocate( nbsn )
    deallocate( nbsn_sta )
    deallocate( nbsn_end )
    deallocate( nbs_sta )
    deallocate( nbs_end )
    deallocate( wk_sta1 )
    deallocate( wk_end1 )
    deallocate( wk_sta )
    deallocate( wk_end )
    deallocate( ir_buff )
                                                  __TIMER_SUB_STOP(1231)
  end subroutine make_index_band_for_Gdiv_3D

!===============================================================================
#ifdef _USE_SCALAPACK_
#ifdef EIGEN_6D
  subroutine make_index_band_for_scalapack(neg, meg, mgs_bs, scl_bs, M, N, usermap, irank_c, irank_r)
#else
  subroutine make_index_band_for_scalapack(neg, meg, mgs_bs, scl_bs, M, N)
#endif

    integer, intent(in) :: neg, meg, mgs_bs, scl_bs, M, N

    integer :: nbs_num, nbsn_num
    integer , allocatable, dimension(:,:) :: wk_row, wk_col
    integer , allocatable, dimension(:) :: wk_rk
    integer :: i, k, j, l, ll, bno, bno_c, bno_r, mm, mmmax(2), mmmax_mpi(2)
    integer :: lrk, lad, lno, row, col, ng, ng_c, ng_r, lrk_r, lrk_c, rank_r, rank_c
#ifdef EIGEN_6D
    integer, intent(in) :: usermap(M,N)
    integer, intent(in) :: irank_c(nrank_e*nrank_g*nrank_k-1)
    integer, intent(in) :: irank_r(nrank_e*nrank_g*nrank_k-1)
#endif

    scl_row = M
    scl_col = N
    my_row = mod((mype-nrank_e*nrank_g*myrank_k),scl_col)
    my_col = (mype-nrank_e*nrank_g*myrank_k)/scl_col

    nbs_num  = (neg - 1)/scl_bs + 1
    nbsn_num = (nbs_num - 1)/scl_row + 1

    allocate( wk_row(scl_bs*nbsn_num,0:scl_row-1) ); wk_row = 0
    allocate( neg_row(neg) ); neg_row = 0
    allocate( nis_row(0:scl_row-1) )
    allocate( nie_row(0:scl_row-1) )
    allocate( nel_row(0:scl_row-1) )

    do i = 1, neg
      if (mod(i,scl_bs) > 0) then
        bno = i / scl_bs + 1
      else
        bno = i / scl_bs
      end if
      lrk = mod((bno-1),scl_row)
      lad = mod(i,scl_bs)
      if (lad == 0) lad = scl_bs
      mm  = ((i-1)/scl_bs)/scl_row
      lno = mm * scl_bs + lad
      wk_row(lno,lrk) = i
    end do

    k = 0
    do i = 0, scl_row-1
      l = 0
      do j = 1, scl_bs*nbsn_num
        if (wk_row(j,i) /= 0) then
          k = k + 1
          l = l + 1
          neg_row(k) = wk_row(j,i)
        end if
      end do
      nel_row(i) = l
    end do

    nis_row(0) = 1
    nie_row(0) = nis_row(0) + nel_row(0) - 1
    do i = 1, scl_row-1
      nis_row(i) = nis_row(i-1) + nel_row(i-1)
      nie_row(i) = nis_row(i) + nel_row(i) - 1
    end do

    nbs_num  = (neg - 1)/scl_bs + 1
    nbsn_num = (nbs_num - 1)/scl_col + 1

    allocate( wk_col(scl_bs*nbsn_num,0:scl_col-1) ); wk_col = 0
    allocate( neg_col(neg) ); neg_col = 0
    allocate( nis_col(0:scl_col-1) )
    allocate( nie_col(0:scl_col-1) )
    allocate( nel_col(0:scl_col-1) )

    do i = 1, neg
      if (mod(i,scl_bs) > 0) then
        bno = i / scl_bs + 1
      else
        bno = i / scl_bs
      end if
      lrk = mod((bno-1),scl_col)
      lad = mod(i,scl_bs)
      if (lad == 0) lad = scl_bs
      mm  = ((i-1)/scl_bs)/scl_col
      lno = mm * scl_bs + lad
      wk_col(lno,lrk) = i
    end do

    k = 0
    do i = 0, scl_col-1
      l = 0
      do j = 1, scl_bs*nbsn_num
        if (wk_col(j,i) /= 0) then
          k = k + 1
          l = l + 1
          neg_col(k) = wk_col(j,i)
        end if
      end do
      nel_col(i) = l
    end do

    nis_col(0) = 1
    nie_col(0) = nis_col(0) + nel_col(0) - 1
    do i = 1, scl_col-1
      nis_col(i) = nis_col(i-1) + nel_col(i-1)
      nie_col(i) = nis_col(i) + nel_col(i) - 1
    end do

!!  write(nfout,'("mype=",i4,", myrank_e=",i4,", myrank_g=",i4)') mype, myrank_e, myrank_g
!!  write(nfout,'("scl_col=",i4,", scl_row=",i4,", scl_bs=",i4)') scl_col, scl_row, scl_bs
!!  write(nfout,'("nrank_g=",i4,", nrank_e=",i4,", mgs_bs=",i4)') nrank_g, nrank_e, mgs_bs
!!  call flush(nfout)

!x  allocate(nrm_rank(0:nrank_g-1,0:nrank_e-1))
!x  allocate(scl_rank(0:scl_row-1 ,0:scl_col-1))
!x  do i = 0, nrank_e*nrank_g - 1
!x    col = mod(i,nrank_e)
!x    row = i/nrank_e
!x    nrm_rank(row,col) = i
!x  end do
!x  do i = 0, scl_col*scl_row - 1
!x    col = mod(i,scl_col)
!x    row = i/scl_col
!x    scl_rank(row,col) = i
!x  end do
!x  write(nfout,'("scl_rank :")')
!x  do i = 0, scl_row - 1
!x    write(nfout,'(20(i4,","))') scl_rank(i,:)
!x  end do
!x  call flush(nfout)
!x  write(nfout,'("nrm_rank :")')
!x  do i = 0, nrank_g - 1
!x    write(nfout,'(20(i4,","))') nrm_rank(i,:)
!x  end do
!x  call flush(nfout)
!x
!x  deallocate(nrm_rank)
!x  deallocate(scl_rank)

    allocate(wk_rk(0:nrank_e*nrank_g*nrank_k-1))
    mmmax = 0
    mmmax_mpi = 0

    scl_comm_rank = 0
!   do k = 0, npes - 1
    do ll = mype, mype
      l = ll-(nrank_e*nrank_g*myrank_k)
      rank_c = mod(l,nrank_e)
      rank_r = l/nrank_e
      wk_rk(:) = 0
      do k = nis_e(rank_c), nie_e(rank_c)
        ng_c = neg_g_all(k)
        if(ng_c > meg) cycle
        if (mod(ng_c,scl_bs) > 0) then
          bno_c = ng_c / scl_bs + 1
        else
          bno_c = ng_c / scl_bs
        end if
        lrk_c = mod((bno_c-1),scl_col)
        do j = nis_eg(rank_r), nie_eg(rank_r)
          ng_r = neg_gg_all(j)
          if(ng_r > meg) cycle
          if (mod(ng_r,scl_bs) > 0) then
            bno_r = ng_r / scl_bs + 1
          else
            bno_r = ng_r / scl_bs
          end if
          lrk_r = mod((bno_r-1),scl_row)
!!!!      lrk = scl_rank(lrk_r,lrk_c)
#ifdef EIGEN_6D
          lrk = usermap(lrk_r+1,lrk_c+1)
#else
#if defined(_ASSIGN_ROW_) || defined(USE_EIGENLIB)
          lrk = scl_row*lrk_c+lrk_r
#else
          lrk = scl_col*lrk_r+lrk_c
#endif
#endif
          wk_rk(lrk) = wk_rk(lrk) + 1
        end do
      end do
!xxxxxxxx
      wk_rk(l) = 0
!xxxxxxxx
      mmmax(1) = maxval(wk_rk(:))

      do k = 0, nrank_e*nrank_g*nrank_k-1
         if(wk_rk(k) /= 0) then
            scl_comm_rank = scl_comm_rank + 1
         end if
      end do
      allocate(scl_comm_rno(scl_comm_rank))
      scl_comm_rno = 0
      scl_comm_rank = 0
      do k = 0, nrank_e*nrank_g*nrank_k-1
         if(wk_rk(k) /= 0) then
            scl_comm_rank = scl_comm_rank + 1
            scl_comm_rno(scl_comm_rank) = k
         end if
      end do

    end do

    scl_comm_rank_r = 0
    do ll = mype, mype
      l = ll-(nrank_e*nrank_g*myrank_k)

#ifdef EIGEN_6D
      rank_c = irank_c(l)
      rank_r = irank_r(l)
      if (rank_c<0 .or. rank_r<0 ) exit
#else
      if (l > (scl_col*scl_row-1)) exit
#if defined(_ASSIGN_ROW_) || defined(USE_EIGENLIB)
      rank_c = l/scl_row
      rank_r = mod(l,scl_row)
#else
      rank_c = mod(l,scl_col)
      rank_r = l/scl_col
#endif
#endif

      wk_rk(:) = 0
      do k = nis_col(rank_c), nie_col(rank_c)
        ng_c = neg_col(k)
        if(ng_c > meg) cycle
        if (mod(ng_c,mgs_bs) > 0) then
          bno_c = ng_c / mgs_bs + 1
        else
          bno_c = ng_c / mgs_bs
        end if
        lrk_c = mod((bno_c-1),nrank_e)
        do j = nis_row(rank_r), nie_row(rank_r)
          ng_r = neg_row(j)
          if(ng_r > meg) cycle
          if (mod(ng_r,mgs_bs) > 0) then
            bno_r = ng_r / mgs_bs + 1
          else
            bno_r = ng_r / mgs_bs
          end if
          lrk_r = mod((bno_r-1),nrank_g)
!!!!      lrk = nrm_rank(lrk_r,lrk_c)
#ifdef ASSIGN_G_PREVIOUS
          lrk = nrank_g*lrk_c+lrk_r
#else
          lrk = nrank_e*lrk_r+lrk_c
#endif
          wk_rk(lrk) = wk_rk(lrk) + 1
        end do
      end do
!xxxxxxxx
      wk_rk(l) = 0
!xxxxxxxx
      mmmax(2) = maxval(wk_rk(:))

      do k = 0, nrank_e*nrank_g*nrank_k-1
         if(wk_rk(k) /= 0) then
            scl_comm_rank_r = scl_comm_rank_r + 1
         end if
      end do
      allocate(scl_comm_rno_r(scl_comm_rank_r))
      scl_comm_rno_r = 0
      scl_comm_rank_r = 0
      do k = 0, nrank_e*nrank_g*nrank_k-1
         if(wk_rk(k) /= 0) then
            scl_comm_rank_r = scl_comm_rank_r + 1
            scl_comm_rno_r(scl_comm_rank_r) = k
         end if
      end do
    end do

    call mpi_allreduce(mmmax,mmmax_mpi,2,mpi_integer,mpi_max,mpi_k_world(myrank_k),ierr)
    scl_comm_max   = mmmax_mpi(1)
    scl_comm_max_r = mmmax_mpi(2)

    deallocate(wk_row, wk_col)
    deallocate(wk_rk)

!!  write(nfout,*) 'neg_row='
!!  write(nfout,'(20(i4,","))') neg_row(:)
!!  write(nfout,*) 'nis_row='
!!  write(nfout,'(20(i4,","))') nis_row(:)
!!  write(nfout,*) 'nie_row='
!!  write(nfout,'(20(i4,","))') nie_row(:)
!!  write(nfout,*) 'nel_row='
!!  write(nfout,'(20(i4,","))') nel_row(:)

!!  write(nfout,*) 'neg_col='
!!  write(nfout,'(20(i4,","))') neg_col(:)
!!  write(nfout,*) 'nis_col='
!!  write(nfout,'(20(i4,","))') nis_col(:)
!!  write(nfout,*) 'nie_col='
!!  write(nfout,'(20(i4,","))') nie_col(:)
!!  write(nfout,*) 'nel_col='
!!  write(nfout,'(20(i4,","))') nel_col(:)

  end subroutine make_index_band_for_scalapack


  subroutine make_index_band_for_scalapack_md(neg, nmatsz, ik, mgs_bs, scl_bs, M, N)
    integer, intent(in) :: neg, nmatsz, ik, mgs_bs, scl_bs, M, N

    integer :: nbs_num, nbsn_num
    integer , allocatable, dimension(:,:) :: wk_row, wk_col
    integer , allocatable, dimension(:) :: wk_rk
    integer :: i, k, j, l, ll, bno, bno_c, bno_r, mm, mmmax(2), mmmax_mpi(2)
    integer :: lrk, lad, lno, row, col, ng, ng_c, ng_r, lrk_r, lrk_c, rank_r, rank_c
    integer :: ib, nn, lrk_c1, lrk_c2, lno_c, lrk_r1, lrk_r2, lno_r, lrk1, lrk2
    integer :: k_ng, j_ng
                                                  __TIMER_SUB_START(1473)

    scl_row_md = M
    scl_col_md = N
    my_row = mod((mype-nrank_e*nrank_g*myrank_k),scl_col_md)
    my_col = (mype-nrank_e*nrank_g*myrank_k)/scl_col_md

    nbs_num  = (nmatsz - 1)/scl_bs + 1
    nbsn_num = (nbs_num - 1)/scl_row_md + 1

    allocate( wk_row(scl_bs*nbsn_num,0:scl_row_md-1) ); wk_row = 0
    allocate( nmatsz_row(nmatsz) ); nmatsz_row = 0
    allocate( nis_row_md(0:scl_row_md-1) )
    allocate( nie_row_md(0:scl_row_md-1) )
    allocate( nel_row_md(0:scl_row_md-1) )

    do i = 1, nmatsz
      if (mod(i,scl_bs) > 0) then
        bno = i / scl_bs + 1
      else
        bno = i / scl_bs
      end if
      lrk = mod((bno-1),scl_row_md)
      lad = mod(i,scl_bs)
      if (lad == 0) lad = scl_bs
      mm  = ((i-1)/scl_bs)/scl_row_md
      lno = mm * scl_bs + lad
      wk_row(lno,lrk) = i
    end do

    k = 0
    do i = 0, scl_row_md-1
      l = 0
      do j = 1, scl_bs*nbsn_num
        if (wk_row(j,i) /= 0) then
          k = k + 1
          l = l + 1
          nmatsz_row(k) = wk_row(j,i)
        end if
      end do
      nel_row_md(i) = l
    end do

    nis_row_md(0) = 1
    nie_row_md(0) = nis_row_md(0) + nel_row_md(0) - 1
    do i = 1, scl_row_md-1
      nis_row_md(i) = nis_row_md(i-1) + nel_row_md(i-1)
      nie_row_md(i) = nis_row_md(i) + nel_row_md(i) - 1
    end do

    nbs_num  = (nmatsz - 1)/scl_bs + 1
    nbsn_num = (nbs_num - 1)/scl_col_md + 1

    allocate( wk_col(scl_bs*nbsn_num,0:scl_col_md-1) ); wk_col = 0
    allocate( nmatsz_col(nmatsz) ); nmatsz_col = 0
    allocate( nis_col_md(0:scl_col_md-1) )
    allocate( nie_col_md(0:scl_col_md-1) )
    allocate( nel_col_md(0:scl_col_md-1) )

    do i = 1, nmatsz
      if (mod(i,scl_bs) > 0) then
        bno = i / scl_bs + 1
      else
        bno = i / scl_bs
      end if
      lrk = mod((bno-1),scl_col_md)
      lad = mod(i,scl_bs)
      if (lad == 0) lad = scl_bs
      mm  = ((i-1)/scl_bs)/scl_col_md
      lno = mm * scl_bs + lad
      wk_col(lno,lrk) = i
    end do

    k = 0
    do i = 0, scl_col_md-1
      l = 0
      do j = 1, scl_bs*nbsn_num
        if (wk_col(j,i) /= 0) then
          k = k + 1
          l = l + 1
          nmatsz_col(k) = wk_col(j,i)
        end if
      end do
      nel_col_md(i) = l
    end do

    nis_col_md(0) = 1
    nie_col_md(0) = nis_col_md(0) + nel_col_md(0) - 1
    do i = 1, scl_col_md-1
      nis_col_md(i) = nis_col_md(i-1) + nel_col_md(i-1)
      nie_col_md(i) = nis_col_md(i) + nel_col_md(i) - 1
    end do

    allocate(wk_rk(0:nrank_e*nrank_g-1))
    mmmax = 0
    mmmax_mpi = 0

    scl_md_comm_rank = 0
    do ll = mype, mype
      l = ll-(nrank_e*nrank_g*myrank_k)
      rank_r = mod(l,nrank_e)
      rank_c = l/nrank_e
      wk_rk(:) = 0
      do k = nis_G_g1k(rank_c,ik), nie_G_g1k(rank_c,ik)
        if (mod(k,scl_bs) > 0) then
          bno_c = k / scl_bs + 1
        else
          bno_c = k / scl_bs
        end if
        lrk_c = mod((bno_c-1),scl_col_md)
        do j = nis_B_g1k(rank_r,ik), nie_B_g1k(rank_r,ik)
          if (mod(j,scl_bs) > 0) then
            bno_r = j / scl_bs + 1
          else
            bno_r = j / scl_bs
          end if
          lrk_r = mod((bno_r-1),scl_row_md)
          lrk = scl_row_md*lrk_c+lrk_r
          wk_rk(lrk) = wk_rk(lrk) + 1
        end do
      end do
      wk_rk(l) = 0
      mmmax(1) = maxval(wk_rk(:))

      do k = 0, scl_col_md*scl_row_md - 1
         if(wk_rk(k) /= 0) then
            scl_md_comm_rank = scl_md_comm_rank + 1
         end if
      end do
      allocate(scl_md_comm_rno(scl_md_comm_rank))
      scl_md_comm_rno = 0
      scl_md_comm_rank = 0
      do k = 0, scl_col_md*scl_row_md - 1
         if(wk_rk(k) /= 0) then
            scl_md_comm_rank = scl_md_comm_rank + 1
            scl_md_comm_rno(scl_md_comm_rank) = k
         end if
      end do

    end do

    scl_md_comm_rank_r = 0
    do ll = mype, mype
      l = ll-(nrank_e*nrank_g*myrank_k)

      if (l > (scl_col_md*scl_row_md-1)) exit

      rank_r = mod(l,scl_row_md)
      rank_c = l/scl_row_md
      wk_rk(:) = 0
      do k = nis_col_md(rank_c), nie_col_md(rank_c)
        ng_c = nmatsz_col(k)
        lrk_c = map_G_g1k(ng_c,ik)
        do j = nis_row_md(rank_r), nie_row_md(rank_r)
          ng_r = nmatsz_row(j)
          lrk_r = map_B_g1k(ng_r,ik)
          lrk = nrank_e*lrk_c+lrk_r
          wk_rk(lrk) = wk_rk(lrk) + 1
        end do
      end do
      wk_rk(l) = 0
      mmmax(2) = maxval(wk_rk(:))

      do k = 0, nrank_e*nrank_g - 1
         if(wk_rk(k) /= 0) then
            scl_md_comm_rank_r = scl_md_comm_rank_r + 1
         end if
      end do
      allocate(scl_md_comm_rno_r(scl_md_comm_rank_r))
      scl_md_comm_rno_r = 0
      scl_md_comm_rank_r = 0
      do k = 0, nrank_e*nrank_g - 1
         if(wk_rk(k) /= 0) then
            scl_md_comm_rank_r = scl_md_comm_rank_r + 1
            scl_md_comm_rno_r(scl_md_comm_rank_r) = k
         end if
      end do

    end do

    call mpi_allreduce(mmmax,mmmax_mpi,2,mpi_integer,mpi_max,mpi_k_world(myrank_k),ierr)
    scl_md_comm_max   = mmmax_mpi(1)
    scl_md_comm_max_r = mmmax_mpi(2)

!-- for trans_scalapack_r
    mmmax = 0
    mmmax_mpi = 0
    scl_md2_comm_rank = 0

    do ll = mype, mype
      l = ll-(nrank_e*nrank_g*myrank_k)

      if (l > (scl_col_md*scl_row_md-1)) exit

      rank_r = mod(l,scl_row_md)
      rank_c = l/scl_row_md

      wk_rk(:) = 0
      do k = 1, neg
        ib = neg_g_all(k)
        if (mod(ib,mgs_bs) > 0) then
           nn = ib / mgs_bs + 1
        else
           nn = ib / mgs_bs
        end if
        lrk_c1 = mod((nn-1),nrank_e)
        lrk_c2 = mod((ib-1)/scl_bs,scl_col_md)
        do j = nis_row_md(rank_r), nie_row_md(rank_r)
          j_ng = nmatsz_row(j)
          lrk_r1 = map_G_g1k(j_ng,ik)
          lrk_r2 = mod((j_ng-1)/scl_bs,scl_row_md)
          lrk1 = nrank_e*lrk_r1+lrk_c1
          lrk2 = scl_row_md*lrk_c2+lrk_r2
          if (lrk2 == mype) then
                   wk_rk(lrk1) = wk_rk(lrk1) + 1
          end if
        end do
      end do

      wk_rk(l) = 0
      mmmax(1) = maxval(wk_rk(:))

      do k = 0, nrank_e*nrank_g - 1
        if(wk_rk(k) /= 0) then
           scl_md2_comm_rank = scl_md2_comm_rank + 1
        end if
      end do
      allocate(scl_md2_comm_rno(scl_md2_comm_rank))
      scl_md2_comm_rno = 0
      scl_md2_comm_rank = 0
      do k = 0, nrank_e*nrank_g - 1
        if(wk_rk(k) /= 0) then
           scl_md2_comm_rank = scl_md2_comm_rank + 1
           scl_md2_comm_rno(scl_md2_comm_rank) = k
        end if
      end do
    end do

    scl_md2_comm_rank_r = 0
    do ll = mype, mype
      l = ll-(nrank_e*nrank_g*myrank_k)
      wk_rk(:) = 0
      rank_c = myrank_e
      rank_r = myrank_g

      do k = nis_e(rank_c), nie_e(rank_c)
        ib = neg_g_all(k)
        lrk_c = mod((ib-1)/scl_bs,scl_col_md)
        do j = nis_G_g1k(rank_r,ik), nie_G_g1k(rank_r,ik)
          lrk_r = mod((j-1)/scl_bs,scl_row_md)
          lrk = scl_row_md*lrk_c+lrk_r
                   wk_rk(lrk) = wk_rk(lrk) + 1
        end do
      end do

      wk_rk(l) = 0
      mmmax(2) = maxval(wk_rk(:))

      do k = 0, scl_col_md*scl_row_md - 1
        if(wk_rk(k) /= 0) then
           scl_md2_comm_rank_r = scl_md2_comm_rank_r + 1
        end if
      end do
      allocate(scl_md2_comm_rno_r(scl_md2_comm_rank_r))
      scl_md2_comm_rno_r = 0
      scl_md2_comm_rank_r = 0
      do k = 0, scl_col_md*scl_row_md - 1
        if(wk_rk(k) /= 0) then
           scl_md2_comm_rank_r = scl_md2_comm_rank_r + 1
           scl_md2_comm_rno_r(scl_md2_comm_rank_r) = k
        end if
      end do

    end do

    call mpi_allreduce(mmmax,mmmax_mpi,2,mpi_integer,mpi_max,mpi_k_world(myrank_k),ierr)
    scl_md2_comm_max   = mmmax_mpi(1)
    scl_md2_comm_max_r = mmmax_mpi(2)
!-- for trans_scalapack_r

    deallocate(wk_row, wk_col)
    deallocate(wk_rk)
                                                  __TIMER_SUB_STOP(1473)

  end subroutine make_index_band_for_scalapack_md

  subroutine dealloate_index_band_for_scalapack_md

    if (allocated(nmatsz_row)) deallocate( nmatsz_row )
    if (allocated(nis_row_md)) deallocate( nis_row_md )
    if (allocated(nie_row_md)) deallocate( nie_row_md )
    if (allocated(nel_row_md)) deallocate( nel_row_md )
    if (allocated(nmatsz_col)) deallocate( nmatsz_col )
    if (allocated(nis_col_md)) deallocate( nis_col_md )
    if (allocated(nie_col_md)) deallocate( nie_col_md )
    if (allocated(nel_col_md)) deallocate( nel_col_md )
    if (allocated(scl_md_comm_rno)) deallocate( scl_md_comm_rno    )
    if (allocated(scl_md_comm_rno_r)) deallocate( scl_md_comm_rno_r  )
    if (allocated(scl_md2_comm_rno)) deallocate( scl_md2_comm_rno   )
    if (allocated(scl_md2_comm_rno_r)) deallocate( scl_md2_comm_rno_r )

  end subroutine dealloate_index_band_for_scalapack_md

#endif
!===============================================================================

  subroutine m_Parallel_get_nproc_from_arg_3D(printable,nfout,ne_in,nk_in,ng_in)
    logical, intent(in) :: printable
    integer, intent(in) :: nfout
    integer, intent(in), optional :: ne_in,nk_in,ng_in
!f
    integer, pointer, dimension(:) :: nproc2rank_k, nproc2rank_e, nproc2rank_g, nproc2rank_s
    character(100) :: q1
    integer     :: narg, iargc
    integer     :: nn , neflag, nkflag, ngflag, nsflag, flag_from_nmlfile
    integer     ::  naflag, nrflag, flag_from_nmlfile_for_nnatm
    integer             :: i, j, k, ip, is
    integer             :: icolor,key
    integer             :: newpes, newmype
!f
    integer     :: err
    integer     :: ne, nk, ng, ns
    logical     :: zaj_para
    namelist / decomp3d /  ne, nk, ng, ns, zaj_para
!f    namelist / decomp3d /  ne, nk, ng
    integer,save :: na,nr
! === For nrc decomposion. by takto 2012/12/05 =================================
    integer,save :: nnatm, nnrc
    namelist / decomppaw / nnatm, nnrc
! ==============================================================================
#ifdef FJ_DBG
    character(len=5) :: p
    character(len=77) :: fdbg
    character(len=64) :: dir_name
#endif
    integer :: nein,nkin,ngin
    logical :: nenkng_present,exi
    logical :: firstcall = .true.
    nenkng_present = .false.
    if(present(ne_in) .and. present(nk_in) .and. present(ng_in)) then
      nenkng_present = .true.
      ne = ne_in
      nk = nk_in
      ng = ng_in
      ns = 1
      nrank_e = ne_in
      nrank_k = nk_in
      nrank_g = ng_in
      nrank_s = 1
    endif

    if(printable) then
      if(conf_para) then
        write(nfout,'(a)') ''
        write(nfout,'(a)') '-- configuration-parallelization scheme enabled --'
        write(nfout,'(a,i8)')  '   nrank_config : ',nrank_conf
        write(nfout,'(a,i8)')  '   mype_config  : ',mype_conf
        write(nfout,'(a,i20)') '   communicator : ',MPI_CommGroup
        write(nfout,'(a)') ''
      endif
    end if

    if(printable) then
       write(nfout,'(" npes = ",i6," << m_Parallel_get_nproc_from_arg_3D>>")') npes
    end if
    if(.not. nenkng_present)then
    ne = npes
    nk = 1
    ng = 1
    ns = 1
    endif
    zaj_para = .true.

    read_from_args = nenkng_present
    if(mype==0) then
       neflag=0
       nkflag=0
       ngflag=0
       nsflag=0
       flag_from_nmlfile = 0
       naflag=0
       nrflag=0
       flag_from_nmlfile_for_nnatm = 0
       if(.not. nenkng_present)then
       narg = iargc()

       ne = 0; ng = 0; nk = 0; ns = 0
       na = 0; nr = 0
       nn = 1
       do while(nn<=narg)
!!$       do nn = 1, narg
          call getarg(nn,q1)
          q1 = trim(adjustl(q1))
          if(index(q1,"ne")/=0) then
             neflag = neflag+1
             ne = getint_from_q1(nn)
!!$             if(printable) write(6,'(" ne = ",i5)') ne
          else if(index(q1,"ng")/=0) then
             ngflag = ngflag+1
             ng = getint_from_q1(nn)
!!$             if(printable) write(6,'(" ng = ",i5)') ng
          else if(index(q1,"nk")/=0) then
             nkflag = nkflag+1
             nk = getint_from_q1(nn)
!!$             if(printable) write(6,'(" nk = ",i5)') nk
          else if(index(q1,"ns")/=0) then
             nsflag = nsflag+1
             ns = getint_from_q1(nn)
          else if(index(q1,"na")/=0 .or. index(q1,"nnatm")/=0 ) then
             naflag = naflag+1
             na = getint_from_q1(nn)
             if(printable) write(nfout,'(" na = ",i5)') na
          else if(index(q1,"nnrc")/=0) then
             nrflag = nrflag+1
             nr = getint_from_q1(nn)
             if(printable) write(nfout,'(" nr = ",i5)') nr
          end if
          nn=nn+1
       end do
       if(nkflag==1 .and. neflag==1 .and. ngflag==1 .and. nsflag==1) then
          read_from_args = .true.
       else if(nkflag==1 .and. neflag==1 .and. ngflag==1 .and. nsflag==0) then
          ns = 1
          read_from_args = .true.
       else if(nkflag==1 .and. neflag==1 .and. ngflag==0) then
          ng = npes/(ne*nk)
          ns = 1
          read_from_args = .true.
       else if(nkflag==1 .and. neflag==0 .and. ngflag==1) then
          ne = npes/(ng*nk)
          ns = 1
          read_from_args = .true.
       else if(nkflag==0 .and. neflag==1 .and. ngflag==1) then
          nk = npes/(ng*ne)
          ns = 1
          read_from_args = .true.
       else if(nkflag==0 .and. neflag==0 .and. ngflag==1) then
          nk = 1
          ne = npes/(nk*ng)
          ns = 1
          read_from_args = .true.
       else if(nkflag==0 .and. neflag==1 .and. ngflag==0) then
          nk = 1
          ng = npes/(nk*ne)
          ns = 1
          read_from_args = .true.
       else if(nkflag==1 .and. neflag==0 .and. ngflag==0) then
          ne = 1
          ng = npes/(nk*ne)
          ns = 1
          read_from_args = .true.
       else if(nkflag==0 .and. neflag==0 .and. ngflag==0) then
          flag_from_nmlfile = 1
       end if
       endif

       if(flag_from_nmlfile==1) then

          ne = npes
          nk = 1
          ng = 1
          ns = 1
          inquire(file='nml.lst',exist=exi)
          if(exi)then
            open(1001, file="nml.lst",status="old", err=998)
            goto 111
998         write(nfout,*) '[ERROR] File open errr. nml.lst in m_Parallel_get_nproc_from_arg_3D'
111         read(1001, nml=decomp3d,iostat=err)
            if(err/=0) then
               write(nfout,*) '[ERROR] File read errr. nml.lst in m_Parallel_get_nproc_from_arg_3D'
            else
               read_from_args = .true.
            end if
            close(1001)
          endif
       end if

       if(naflag==1 .and. nrflag==1) then
       else if(naflag==1 .and. nrflag==0) then
          nr = npes/na
       else if(naflag==0 .and. nrflag==1) then
          na = npes/nr
       else if(naflag==0 .and. nrflag==0) then
          flag_from_nmlfile_for_nnatm = 1
       end if


       if(flag_from_nmlfile_for_nnatm==1) then
          nnatm = ng
          nnrc  = npes/ng
          inquire(file='nml.lst',exist=exi)
          if(exi)then
             open(1001, file="nml.lst",status="old", err=997)
             goto 112
997          write(nfout,*) '[ERROR] File open errr(2). nml.lst in m_Parallel_get_nproc_from_arg_3D'
112          read(1001, nml=decomppaw,iostat=err)
             if(err/=0) then
                write(nfout,'(" nnatm = ",i3)') nnatm
                write(nfout,'(" nnrc  = ",i3)') nnrc
                write(nfout,*) '[ERROR] File read errr(2). nml.lst in m_Parallel_get_nproc_from_arg_3D'
             end if
             close(1001)
!!$          if(nnatm*nnrc /= npes) then
!!$             write(0,'(" npes(=",i2,") /= nnatm * nnrc (=",i2,"*",i2")")') npes, nnatm, nnrc
!!$             call mpi_barrier(MPI_CommGroup,err)
!!$             call mpi_abort(MPI_CommGroup,err)
!!$          else
!!$             nrank_natm = nnatm
!!$             nrank_nrc  = nnrc
!!$          end if
          endif
       else
          nnatm = na
          nnrc  = nr
       end if
    end if

    if(npes>1) call mpi_bcast(nk,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(npes>1) call mpi_bcast(ne,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(npes>1) call mpi_bcast(ng,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(npes>1) call mpi_bcast(ns,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(npes>1) call mpi_bcast(read_from_args,1,mpi_logical,0,MPI_CommGroup,ierr)
    if(npes>1) call mpi_bcast(nkflag,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(npes>1) call mpi_bcast(neflag,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(npes>1) call mpi_bcast(ngflag,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(npes>1) call mpi_bcast(flag_from_nmlfile,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(npes>1) call mpi_bcast(naflag,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(npes>1) call mpi_bcast(nrflag,1,mpi_integer,0,MPI_CommGroup,ierr)

    if(nkflag>=2 .or. neflag>=2 .or. ngflag>=2) then
       if(printable) write(nfout,'(" wrong nk, ne, and ng in m_Parallel_get_nproc_from_arg_3D")')
       call phase_error_with_msg(nfout, &
       ' Wrong nk, ne, and ng in m_Parallel_get_nproc_from_arg_3D',__LINE__,__FILE__)
    end if

    if(flag_from_nmlfile==0) then
       if(nk*ne*ng*ns/=npes) then
          if(printable) write(nfout,'(" nk*ne*ng*ns /= npes")')
          if(printable) write(nfout,'(" nk = ",i4, " ne = ", i4, " ng = ",i4, " ns = ",i4, " npes = ",i4)') nk,ne,ng,ns,npes
          call phase_error_with_msg(nfout,'nk*ne*ng*ns /= npes',__LINE__,__FILE__)
       end if
    end if

    if(naflag>=2 .or. nrflag>=2) then
       if(printable) write(nfout,'(" wrong na and nr in m_Parallel_get_nproc_from_arg_3D")')
       call phase_error_with_msg(nfout, &
       ' Wrong na and nr in m_Parallel_get_nproc_from_arg_3D',__LINE__,__FILE__)
    end if

!!$    if ( err/=0 ) goto 999

    if ( ne*nk*ng*ns/=npes ) then
      if(mype == 0) then
         if(printable) write(nfout,'(" npes(=",i2,") /= ne * nk * ng (=",i2,"*",i2,"*",i2,")")') npes,ne,nk,ng
         call mpi_abort(mpi_comm_world, 1, err)
         call phase_error_with_msg(nfout,' ! invalid argument lists',__LINE__,__FILE__)
      endif
    else
      nrank_e = ne
      nrank_k = nk
      nrank_g = ng
      nrank_s = ns
!      if( mype == 0 .and. nk/=1 ) then
!        write(nfout,*) '[ERROR] nk =', nk
!        call mpi_abort(mpi_comm_world, 2, err)
!      endif
    endif

    if ( zaj_para ) then
      init_zaj_para = .true.
    else
      init_zaj_para = .false.
    endif

    if(printable)then
      if(read_from_args .or. nenkng_present) write(nfout,'(" nrank_e = ",i6)') nrank_e
      if(read_from_args .or. nenkng_present) write(nfout,'(" nrank_k = ",i6)') nrank_k
      if(read_from_args .or. nenkng_present) write(nfout,'(" nrank_g = ",i6)') nrank_g
      if(read_from_args .or. nenkng_present) write(nfout,'(" nrank_s = ",i6)') nrank_s
    endif

! (( myrank_e, myrank_k ))
    allocate(nproc2rank_e(0:npes-1))
    allocate(nproc2rank_k(0:npes-1))
    allocate(nproc2rank_g(0:npes-1))
    allocate(nproc2rank_s(0:npes-1))
    ip = 0
    do is= 0, nrank_s-1
    do i = 0, nrank_k-1
#ifdef ASSIGN_G_PREVIOUS
      do k = 0, nrank_e-1
       do j = 0, nrank_g-1
          nproc2rank_s(ip) = is
          nproc2rank_k(ip) = i
          nproc2rank_e(ip) = k
          nproc2rank_g(ip) = j
          ip = ip + 1
       end do
      end do
#else
      do k = 0, nrank_g-1
       do j = 0, nrank_e-1
          nproc2rank_s(ip) = is
          nproc2rank_k(ip) = i
          nproc2rank_e(ip) = j
          nproc2rank_g(ip) = k
          ip = ip + 1
       end do
      end do
#endif
    end do
    end do
    myrank_k = nproc2rank_k(mype)
    myrank_e = nproc2rank_e(mype)
    myrank_g = nproc2rank_g(mype)
    myrank_spin = nproc2rank_s(mype)
    deallocate(nproc2rank_k)
    deallocate(nproc2rank_e)
    deallocate(nproc2rank_g)
    deallocate(nproc2rank_s)
    icolor = myrank_g + nrank_g*myrank_k + nrank_g*nrank_k*myrank_e
    key    = myrank_spin
    if(.not.firstcall) call mpi_comm_free(mpi_keg_world, ierr)
    call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_keg_world,ierr)
    call mpi_comm_size(mpi_keg_world, newpes,ierr)
    call mpi_comm_rank(mpi_keg_world, newmype,ierr)

    icolor = myrank_spin
    key    = mype
    if(.not.firstcall) call mpi_comm_free(mpi_spin_group,ierr)
    call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_spin_group,ierr)
    call mpi_comm_size(mpi_spin_group, newpes,ierr)
    call mpi_comm_rank(mpi_spin_group, newmype,ierr)

    icolor = myrank_g + nrank_g*myrank_k
    key    = mype
    if(.not.firstcall) call mpi_comm_free(mpi_skg_world,ierr)
    call mpi_comm_split(mpi_spin_group, icolor, key, mpi_skg_world,ierr)
    call mpi_comm_size(mpi_skg_world, newpes,ierr)
    call mpi_comm_rank(mpi_skg_world, newmype,ierr)

    icolor = myrank_g + nrank_g*myrank_e
    key    = mype
    if(.not.firstcall) call mpi_comm_free(mpi_sge_world,ierr)
    call mpi_comm_split(mpi_spin_group, icolor, key, mpi_sge_world,ierr)
    call mpi_comm_size(mpi_sge_world, newpes,ierr)
    call mpi_comm_rank(mpi_sge_world, newmype,ierr)

    icolor = myrank_e + nrank_e*myrank_k
    key    = mype
    if(.not.firstcall) call mpi_comm_free(mpi_ske_world,ierr)
    call mpi_comm_split(mpi_spin_group, icolor, key, mpi_ske_world,ierr)
    call mpi_comm_size(mpi_ske_world, newpes,ierr)
    call mpi_comm_rank(mpi_ske_world, newmype,ierr)

! (( mpi_kg_world ))
!    icolor = myrank_g + nrank_g*myrank_k
    icolor = myrank_g + nrank_g*myrank_k + nrank_g*nrank_k*myrank_spin
    key    = mype
    if(.not.firstcall) call mpi_comm_free(mpi_kg_world,ierr)
    call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_kg_world,ierr)
    call mpi_comm_size(mpi_kg_world, newpes,ierr)
    call mpi_comm_rank(mpi_kg_world, newmype,ierr)

! (( mpi_ke_world ))
    icolor = myrank_e + nrank_e*myrank_k + nrank_e*nrank_k*myrank_spin
    key    = mype
    if(.not.firstcall) call mpi_comm_free(mpi_ke_world,ierr)
    call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_ke_world,ierr)
!    call mpi_comm_split(mpi_spin_group, icolor, key, mpi_ke_world,ierr)
    call mpi_comm_size(mpi_ke_world, newpes,ierr)
    call mpi_comm_rank(mpi_ke_world, newmype,ierr)

! (( mpi_k_world ))
    if(allocated(mpi_k_world)) deallocate(mpi_k_world)
    allocate(mpi_k_world(0:nrank_k-1))
    do j = 0, nrank_k-1
       icolor = 0
       key = 0
       do i = 0, nrank_e*nrank_g-1
          if(mype == i + nrank_e*nrank_g*j + myrank_spin*nrank_e*nrank_k*nrank_g) then
             icolor = 1
             key = i
          end if
       end do
       if(.not.firstcall) call mpi_comm_free(mpi_k_world(j),ierr)
!       call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_k_world(j),ierr)
       call mpi_comm_split(mpi_spin_group, icolor, key, mpi_k_world(j),ierr)
       call mpi_comm_size(mpi_k_world(j), newpes,ierr)
       call mpi_comm_rank(mpi_k_world(j), newmype,ierr)
    end do
    if(nrank_conf>1) then
      icolor = myrank_k+mype_conf*nrank_k
      key = myrank_e+myrank_g*nrank_e
      call mpi_comm_split(mpi_comm_world, icolor, key, mpi_kr_world,ierr)
      call mpi_comm_size(mpi_kr_world, newpes,ierr)
      call mpi_comm_rank(mpi_kr_world, newmype,ierr)
    endif
! === DEBUG by tkato 2012/04/04 ================================================
! (( mpi_ge_world ))
    icolor = myrank_g + nrank_g*myrank_e
    key    = mype
    if(.not.firstcall) call mpi_comm_free(mpi_ge_world,ierr)
    call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_ge_world, ierr)
!    call mpi_comm_split(mpi_spin_group, icolor, key, mpi_ge_world, ierr)
    call mpi_comm_size(mpi_ge_world, newpes, ierr)
    call mpi_comm_rank(mpi_ge_world, newmype, ierr)
! ==============================================================================
! === For nrc decomposion. by takto 2012/12/05 =================================
    if(npes>1) call mpi_bcast(nnatm,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(npes>1) call mpi_bcast(nnrc, 1,mpi_integer,0,MPI_CommGroup,ierr)

!!$    nrank_natm = npes
!!$    nrank_nrc  = 1
!!$
!!$    read(1001, nml=decomppaw, iostat=err)
!!$
!!$    if(err /= 0) goto 998

    if(nnatm*nnrc /= npes) then
       if(mype == 0) then
          write(0,'(" npes(=",i2,") /= nnatm * nnrc (=",i2,"*",i2")")') npes, nnatm, nnrc
       endif
       call mpi_barrier(MPI_CommGroup,err)
       call mpi_abort(MPI_CommGroup, 999, err)
    else
       nrank_natm = nnatm
       nrank_nrc  = nnrc
    endif

!!$998 continue

    if(printable)then
      if( read_from_args .or. nenkng_present) write(nfout,'(" nrank_natm = ",i3)') nrank_natm
      if( read_from_args .or. nenkng_present) write(nfout,'(" nrank_nrc  = ",i3)') nrank_nrc
    endif
! ==============================================================================
#ifdef FJ_DBG
    call getenv('FJ_DBG_DIR', dir_name)
    if (len_trim(dir_name) .eq. 0) then
       dir_name = '.'
    else
       call system('if [ ! -d '//dir_name//' ]; then mkdir -p '//dir_name//'; fi')
    endif

    write(p,'(i5.5)') mype
    fdbg=trim(dir_name)//'/dbginfo.'//trim(p)
    open(1002, file=trim(fdbg), form='formatted')
    write(1002,'(" !|| << m_Parallel_get_nproc_from_arg_3D >>")')
    write(1002,'(" npes(=",i2,") /= ne * nk * ng (=",i2,"*",i2,"*",i2,")")') npes,ne,nk,ng
#endif

    if(read_from_args .or. nenkng_present) firstcall = .false.

!!$    stop ' test stop'
    return

!!$999 write(nfout,*) '[ERROR] File open errr. nml.lst in m_Parallel_get_nproc_from_arg_3D'

  contains
    integer function getint_from_q1(nn)
      integer,intent(inout) ::nn
      integer :: i, n
      i = scan(q1,'=;:')
      if(i/=0) then
         if(len_trim(q1(i+1:)).ne.0) then
            q1 = q1(i+1:)
         else
            nn=nn+1
            if(nn>narg) goto 300
            call getarg(nn,q1)
         end if
      else
         nn=nn+1
         if(nn>narg) goto 300
         call getarg(nn,q1)
         i = scan(q1,'=;:')
         if(i/=0) then
            if(len_trim(q1(i+1:)).ne.0) then
               q1 = q1(i+1:)
            else
               nn = nn+1
               if(nn>narg) goto 300
               call getarg(nn,q1)
            end if
         end if
      end if

      read(q1,*,err=200) n
      if(n<=0) call phase_error_with_msg(nfout,' a negative value or zero is read for ne, nk or ng in an argument list'&
                                        ,__LINE__,__FILE__)
      getint_from_q1=n
      return
200   call phase_error_with_msg(nfout,' wrong argument for ne,nk,or ng',__LINE__,__FILE__)
300   call phase_error_with_msg(nfout,' something is lacking in the argument list for ne,nk and ng',__LINE__,__FILE__)
    end function getint_from_q1

  end subroutine m_Parallel_get_nproc_from_arg_3D

!===============================================================================

  subroutine m_Parallel_init_mpi_elec_3D(nfout,ipri,printable,neg,kv3,nspin,kg1,iba)
    integer, intent(in) :: nfout,ipri,neg, kv3, nspin, kg1, iba(kv3)
    logical, intent(in) :: printable
    logical, save :: firstcall=.true.

#ifdef NEC_TUNE_SOFT
    character*4 F_RSVTASK
    integer :: tmp_a1,tmp_a2,tmp_b1,tmp_b2
#elif NEC_TUNE_FFT
    character*4 F_RSVTASK
#endif

    integer             :: i, j, k, ip, icolor, key, kv3_half
    integer             :: ik
    integer             :: newpes, newmype
    logical             :: set_mapping_func
    integer, pointer, dimension(:) :: nproc2rank_k, nproc2rank_e, nproc2rank_g

    integer, parameter  :: nc = 24
    integer, parameter  :: nc1 = 16
!!$    integer, parameter  :: nc2 = 18
    integer             :: i0, i1, sw_title
    character*7         :: strmap
                                                  __TIMER_SUB_START(1229)

    sw_title = 1
    if(ipri > 1 .and. printable) then
       write(nfout,'(" !|| << m_Parallel_init_mpi_elec_3D >>")')
       write(nfout,'(" !|| neg, kv3 = ",2i5)') neg, kv3
       do ik = 1, kv3
         write(nfout,'(" !|| iba(",i5,") = ",i6)') ik, iba(ik)
       enddo
    end if

    if(neg <= 0 ) call phase_error_with_msg(nfout," neg is not positive value ",__LINE__,__FILE__)
    if(kv3 <= 0 ) call phase_error_with_msg(nfout," kv3 is not positive value ",__LINE__,__FILE__)

!!!!!!!!!!!!!!!!! added by mizouchi@adv 2003.02.26 !!!!!!!!!!!!!!
#ifdef IRIX64
    if((kv3 == 1.and.nspin ==1).or.(kv3 == 2.and.nspin ==2)) then
       nrank_e  = npes
       nrank_k  = 1
       if(ipri >= 1 .and. printable) then
          write(nfout,'(" !| modified nrank_e and nrank_k due to 1 kpoint calculation  ")')
          write(nfout,'(" !| modified nrank_e = ", i4)') nrank_e
          write(nfout,'(" !| modified nrank_k = ", i4)') nrank_k
       end if
#endif
!!!!!!!!!!!!!!!!! added by mizouchi@adv 2003.02.26 !!!!!!!!!!!!!!

    if(nrank_k > kv3) call phase_error_with_msg(nfout,' nrank_k > kv3 (m_Parallel_init_mpi_elec)',__LINE__,__FILE__)
    if(nrank_s>nspin)  &
    call phase_error_with_msg(nfout,"ns>nspin",__LINE__,__FILE__)

    allocate(map_e(neg))
    allocate(map_z(neg))
    allocate(map_k(kv3))
    allocate(map_ek(neg,kv3))
    allocate(map_s(nspin))
!fj    allocate(map_ek(neg,kv3))
!fj    allocate(mpi_k_world(0:nrank_k-1))
    allocate(nis_e(0:nrank_e-1)); nis_e = neg+1
    allocate(nie_e(0:nrank_e-1)); nie_e = 0
    allocate(nel_e(0:nrank_e-1)); nel_e = 0
    allocate(idisp_e(0:nrank_e-1)); idisp_e = 0
    allocate(nis_g1(0:nrank_g-1)); nis_g1 = kg1 + 1
    allocate(nie_g1(0:nrank_g-1)); nie_g1 = 0
    allocate(nel_g1(0:nrank_g-1)); nel_g1 = 0
    allocate(nis_k(0:nrank_k-1)); nis_k = 0
    allocate(nie_k(0:nrank_k-1))
    allocate(nel_k(0:nrank_k-1))

    allocate(nis_spin(0:nrank_s-1)); nis_spin = 0
    allocate(nie_spin(0:nrank_s-1))
    allocate(nel_spin(0:nrank_s-1))

! (( myrank_e, myrank_k ))
!f    allocate(nproc2rank_e(0:npes-1))
!f    allocate(nproc2rank_k(0:npes-1))
!f    allocate(nproc2rank_g(0:npes-1))
!f    ip = 0
!f    do i = 0, nrank_k-1
!f      do k = 0, nrank_g-1
!f       do j = 0, nrank_e-1
!f          nproc2rank_k(ip) = i
!f          nproc2rank_e(ip) = j
!f          nproc2rank_g(ip) = k
!f          ip = ip + 1
!f       end do
!f      end do
!f    end do
!f    myrank_k = nproc2rank_k(mype)
!f    myrank_e = nproc2rank_e(mype)
!f    myrank_g = nproc2rank_g(mype)
!f    deallocate(nproc2rank_k)
!f    deallocate(nproc2rank_e)
!f    deallocate(nproc2rank_g)
! (( map_e, nis_e, nie_e, nel_e, ista_e, iend_e, istep_e, np_e, mp_e  ))
!!$#ifdef TRANSPOSE
    set_mapping_func = .true.
    call set_block_range(neg,nrank_e,nel_e,nis_e,nie_e,set_mapping_func, map_e)
    do i = 0, nrank_e-1
       idisp_e(i) = nis_e(i)-1
    end do
    if(ipri >= 1) call wd_e_range_3D
    istep_e = 1
!!$#else
!!$!Fj not modified
!!$    do i = 1, neg
!!$       map_e(i) = mod(i-1,nrank_e)
!!$    end do
!!$    do i = 1, neg
!!$       ip = map_e(i)
!!$       nel_e(ip) = nel_e(ip) + 1
!!$       if(nis_e(ip) > i) nis_e(ip) = i
!!$       if(nie_e(ip) < i) nie_e(ip) = i
!!$    end do
!!$    istep_e = nrank_e
!!$
!!$#endif
    if(ipri >= 2 .and. printable) then
       strmap = "   i   "
       call wd_maparray_3D(neg,map_e,strmap," map_e ",sw_title)
    end if

    ista_e = nis_e(myrank_e)
    iend_e = nie_e(myrank_e)
    np_e = nel_e(myrank_e)
    mp_e = maxval(nel_e)

#ifdef NEC_TUNE_SOFT
!Fj not modified
    call getenv('F_RSVTASK',F_RSVTASK)
    read (F_RSVTASK,'(i4)') itask
    allocate(ista_e_smp(itask))
    allocate(iend_e_smp(itask))

    tmp_a1 = (iend_e-ista_e+1)/istep_e
    tmp_a2 = mod(iend_e-ista_e+1,istep_e)
    if (tmp_a2 .ne. 0) tmp_a1=tmp_a1+1
    if (itask == 0) itask = 1
    tmp_b1 = tmp_a1/itask
    tmp_b2 = mod(tmp_a1,itask)
    if (tmp_b2 .ne. 0) tmp_b1=tmp_b1+1

    do i=1,itask
      ista_e_smp(i) = ista_e + tmp_b1 * istep_e * (i -1)
      iend_e_smp(i) = min((ista_e + tmp_b1 * istep_e * i -1),iend_e)
    end do
#elif NEC_TUNE_FFT
!Fj not modified
    call getenv('F_RSVTASK',F_RSVTASK)
    read (F_RSVTASK,'(i4)') itask
#endif

    if(ipri >= 2 .and. printable) &
         & write(nfout,'(" !|| m_P_init_mpi_elec_3D -- ista_e, iend_e, istep_e, np_e, mp_e = ",5i6)') &
         &  ista_e,iend_e,istep_e, np_e, mp_e


! (( map_z ))
!!$#ifdef TRANSPOSE
    j = 0
    do ip = 1, nrank_e
       do i = 1, nel_e(ip-1)
          j = j + 1
          map_z(j) = i
       end do
    end do
!!$#else
!!$!Fj not modified
!!$    ip = 1
!!$    do i = 1, neg
!!$       map_z(i) = ip
!!$       if(mod(i,nrank_e) == 0) ip = ip + 1
!!$    end do
!!$#endif
    if(ipri >= 2 .and. printable) then
       strmap = "   i   "
       call wd_maparray_3D(neg,map_z,strmap," map_z ",sw_title)
    end if
! (( map_k, nis_k, nie_k, nel_k, ista_k, iend_k ))
    if(nspin == 1) then
       set_mapping_func = .true.
       call set_block_range(kv3,nrank_k,nel_k,nis_k,nie_k,set_mapping_func,map_k)
       ista_k = nis_k(myrank_k)
       iend_k = nie_k(myrank_k)
    else if(nspin == 2) then
       set_mapping_func = .false.
       kv3_half = kv3/nspin
       call set_block_range(kv3_half,nrank_k,nel_k,nis_k,nie_k, set_mapping_func)
       nel_k = nel_k*nspin
       nie_k = nie_k*nspin
       nis_k = nis_k*nspin-1
       j = 0
       do i = 1, kv3
          if(nie_k(j) < i) j = j + 1
          map_k(i) = j
       end do
       ista_k = nis_k(myrank_k)
       iend_k = nie_k(myrank_k)
    end if

    if(ipri >= 2 .and. printable) then
       strmap = "   i   "
       call wd_maparray_3D(kv3,map_k,strmap," map_k ",sw_title)
       write(nfout,'(" !|| -- myrank_k = ",i8)') myrank_k
       write(nfout,'(" !|| -- ista_k, iend_k = ",2i8)') ista_k,iend_k
    end if

    set_mapping_func = .true.
    call set_block_range(nspin,nrank_s,nel_spin,nis_spin,nie_spin, set_mapping_func, map_s)
    ista_spin = nis_spin(myrank_spin)
    iend_spin = nie_spin(myrank_spin)
    np_spin   = nel_spin(myrank_spin)
    mp_spin   = np_spin

   if(ista_k > iend_k) call phase_error_with_msg(nfout,' !! illegal combination of ista_k, and iend_k',__LINE__,__FILE__)

   icolor = myrank_g
   key = myrank_e+myrank_k*nrank_e
!   if(firstcall) call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_g_world,ierr)
   if(.not.firstcall) call mpi_comm_free(mpi_g_world, ierr)
   call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_g_world,ierr)
   call mpi_comm_size(mpi_g_world, nrank_ke,  ierr)
   call mpi_comm_rank(mpi_g_world, myrank_ke, ierr)

!! coded by T. Yamasaki, 2021.02.25 -->
   if(allocated(map_rank_gek)) deallocate(map_rank_gek)
   allocate(map_rank_gek(0:nrank_g-1,0:nrank_e-1,0:nrank_k-1,0:nrank_s-1)); map_rank_gek = 0
   map_rank_gek(myrank_g,myrank_e,myrank_k,myrank_spin) = mype
   call mpi_allreduce(MPI_IN_PLACE,map_rank_gek,npes,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
!!$
!!$   ip = 0
!!$   do k = 0, nrank_k-1
!!$      do i = 0, nrank_g-1
!!$         do j = 0, nrank_e-1
!!$            map_rank_gek(i,j,k) = ip
!!$            ip = ip + 1
!!$         end do
!!$      end do
!!$   end do
   if(ipri>=2 .and. printable) then
      ip = 0
      write(nfout,'(" !||  No.  rank_k,rank_e,rank_g,   map_rank_gek")')
      do k = 0, nrank_k-1
         do j = 0, nrank_e-1
            do i = 0, nrank_g-1
               write(nfout,'(" !|| ",4i6," ",i6)') ip,k,j,i, map_rank_gek(i,j,k,1)
               ip = ip + 1
            end do
         end do
      end do
      write(nfout,'(" myrank_k, myrank_e, myrank_g, mype")')
      write(nfout,'(4i8)') myrank_k,myrank_e,myrank_g,mype
   end if
   if(map_rank_gek(myrank_g,myrank_e,myrank_k,myrank_spin) /= mype) then
      write(nfout,'(" map_rank_gek is not right")')
      call phase_error_with_msg(nfout,' map_rank_gek is not right',__LINE__,__FILE__)
   end if
!! <--
#if 0
! (( map_ek ))
    do j = 1, kv3
       do i = 1, neg
          map_ek(i,j) = map_e(i) + map_k(j)*(nrank_e)
       end do
    end do

    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| --- map_ek ---")')
       sw_title = 0
       do j = 1, kv3
          write(strmap,'("ik=",i4)') j
          call wd_maparray_3D(neg,map_ek(1,j),strmap,"map_ek ",sw_title)
       end do
    end if

! (( map_ekg ))
    do j = 1, kv3
       do i = 1, neg
         do k = 1, iba(ik)
          map_ekg(k,i,j) =
       end do
       end do
    end do

    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| --- map_ek ---")')
       sw_title = 0
       do j = 1, kv3
         do k = 1, iba(ik)
           write(strmap,'("ik=",i4)') j
           call wd_maparray_3D(neg,map_ekg(1,j,k),strmap,"map_ekg",sw_title)
         end do
       end do
    end if
#endif

!!$#ifdef TRANSPOSE
! (( ista_g1, iend_g1 ))
    set_mapping_func = .false.
    if(ipri >= 2 .and. printable) write(nfout,'(" !|| kg1 = ",i10)') kg1
    call set_block_range(kg1,nrank_g,nel_g1,nis_g1,nie_g1, set_mapping_func)
    ista_g1 = nis_g1(myrank_g)
    iend_g1 = nie_g1(myrank_g)
    np_g1   = nel_g1(myrank_g)
    mp_g1   = maxval(nel_g1)
    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| --- mis_g1,nie_g1,nel_g1 ---")')
       write(nfout,'(" !|| ( rank_e  )",15i7)')(i,i=1,nrank_g)
       write(nfout,'(" !|| ( nis_g1  )",15i7)')(nis_g1(i),i=0,nrank_g-1)
       write(nfout,'(" !|| ( nie_g1  )",15i7)')(nie_g1(i),i=0,nrank_g-1)
       write(nfout,'(" !|| ( nel_g1  )",15i7)')(nel_g1(i),i=0,nrank_g-1)
       write(nfout,'(" !|| (myrank_e = ",i3," kg1,ista_g1,iend_g1,np_g1,mp_g1 = ",5i7)') &
            & myrank_e,kg1,ista_g1,iend_g1,np_g1,mp_g1
    end if
    firstcall=.false.
!!$#endif

! (( mpi_kg_world ))
!f    icolor = myrank_g + nrank_g*myrank_k
!f    key    = mype
!f    call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_kg_world,ierr)
!f    call mpi_comm_size(mpi_kg_world, newpes,ierr)
!f    call mpi_comm_rank(mpi_kg_world, newmype,ierr)
!!$       print '(" mype = ",i3," newmype = ", i3," newpes = ",i3,"mpi_kg_world = ", i3 &
!!$            &, " icolor = ",i3," key = ",i3)' &
!!$            & ,mype,newmype, newpes,mpi_kg_world,icolor, key

! (( mpi_ke_world ))
!f    icolor = myrank_e + nrank_e*myrank_k
!f    key    = mype
!f    call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_ke_world,ierr)
!f    call mpi_comm_size(mpi_ke_world, newpes,ierr)
!f    call mpi_comm_rank(mpi_ke_world, newmype,ierr)
!!$       print '(" mype = ",i3," newmype = ", i3," newpes = ",i3,"mpi_ke_world = ", i3 &
!!$            &, " icolor = ",i3," key = ",i3)' &
!!$            & ,mype,newmype, newpes,mpi_ke_world,icolor, key
                                                  __TIMER_SUB_STOP(1229)

  contains
    subroutine wd_maparray_3D(nelement,map_a,map_title0,map_title,sw_title)
      integer, intent(in) ::                 nelement
      integer, dimension(nelement), intent(in) :: map_a
      character*7, intent(in) ::             map_title0,map_title
      integer, intent(in) ::                 sw_title
      integer :: nca, j, i0, i1, i, iloop

      if(sw_title==1)  write(nfout,'(" !|| --- ",a7," ---")') map_title
      if(nelement < 100) then
         nca = nc      ! = 24
      else if(nelement < 1000) then
         nca = nc*3/4  ! = 18
      else if(nelement < 10000) then
         nca = nc*3/5  ! = 14
      else if(nelement < 100000) then
         nca = nc*3/6  ! = 12
      else
         nca = nc*3/9  ! = 8
      end if
      iloop = ceiling(dble(nelement)/nca)
      do j = 1, iloop
         i0 = (j-1)*nca + 1
         i1 = i0 + nca -1
         if(i1 > nelement) i1 = nelement
         if(nelement < 100) then
            write(nfout,'(" !|| (",a7,")",24i3)') map_title0,(i,i=i0,i1)
            write(nfout,'(" !|| (",a7,")",24i3)') map_title, (map_a(i),i=i0,i1)
         else if(nelement < 1000) then
            write(nfout,'(" !|| (",a7,")",18i4)') map_title0,(i,i=i0,i1)
            write(nfout,'(" !|| (",a7,")",18i4)') map_title, (map_a(i),i=i0,i1)
         else if(nelement < 10000) then
            write(nfout,'(" !|| (",a7,")",14i5)') map_title0,(i,i=i0,i1)
            write(nfout,'(" !|| (",a7,")",14i5)') map_title, (map_a(i),i=i0,i1)
         else if(nelement < 100000) then
            write(nfout,'(" !|| (",a7,")",12i5)') map_title0,(i,i=i0,i1)
            write(nfout,'(" !|| (",a7,")",12i5)') map_title, (map_a(i),i=i0,i1)
         else
            write(nfout,'(" !|| (",a7,")",8i9)') map_title0,(i,i=i0,i1)
            write(nfout,'(" !|| (",a7,")",8i9)') map_title, (map_a(i),i=i0,i1)
         end if
      end do
    end subroutine wd_maparray_3D

    subroutine wd_e_range_3D
      if(printable) then
         write(nfout,'(" !|| nrank_e")')
         write(nfout,'(" !||    i    : ",20i4)')(i,i=0,nrank_e-1)
         write(nfout,'(" !|| nis_e   : ",20i4)')(nis_e(i),i=0,nrank_e-1)
         write(nfout,'(" !|| nie_e   : ",20i4)')(nie_e(i),i=0,nrank_e-1)
         write(nfout,'(" !|| nel_e   : ",20i4)')(nel_e(i),i=0,nrank_e-1)
         write(nfout,'(" !||    i : ",20i4)')(i,i=0,nrank_e-1)
         write(nfout,'(" !|| nis_e : ",20i4)')(nis_e(i),i=0,nrank_e-1)
         write(nfout,'(" !|| nie_e : ",20i4)')(nie_e(i),i=0,nrank_e-1)
         write(nfout,'(" !|| nel_e : ",20i4)')(nel_e(i),i=0,nrank_e-1)
      end if
!!$#ifdef TRANSPOSE
         write(nfout,'(" !|| idisp_e : ",20i4)')(idisp_e(i),i=0,nrank_e-1)
!!$#endif
    end subroutine wd_e_range_3D

  end subroutine m_Parallel_init_mpi_elec_3D

!===============================================================================

  subroutine m_Parallel_init_mpi_nval(nfout,ipri,printable,nval)
    integer, intent(in) :: nfout, ipri,nval
    logical, intent(in) :: printable
    integer :: i,j,ip
    nrank_nval = nrank_g
    myrank_nval = myrank_g
    if(mpi_nval_enabled) call m_Parallel_dealloc_mpi_nval()
    call m_Parallel_alloc_mpi_nval(nval)
    call set_block_range(nval,nrank_nval,nel_nval,nis_nval,nie_nval,.true.,map_nval)
    mpi_nval_enabled = .true.
    if(ipri>=2) call wd_val_range_3D()
    ista_nval = nis_nval(myrank_nval)
    iend_nval = nie_nval(myrank_nval)
    np_nval = nel_nval(myrank_nval)
    mp_nval = maxval(nel_nval)

    j = 0
    do ip = 1, nrank_nval
       do i = 1, nel_nval(ip-1)
          j = j + 1
          map_z_nval(j) = i
       end do
    end do

    contains

    subroutine wd_val_range_3D()
      if(printable) then
         write(nfout,'(" !|| nrank_nval")')
         write(nfout,'(" !||    i    : ",20i4)')(i,i=0,nrank_nval-1)
         write(nfout,'(" !|| ista_nval   : ",20i4)')(nis_nval(i),i=0,nrank_nval-1)
         write(nfout,'(" !|| iend_nval   : ",20i4)')(nie_nval(i),i=0,nrank_nval-1)
         write(nfout,'(" !|| nel_nval    : ",20i4)')(nel_nval(i),i=0,nrank_nval-1)
      end if
    end subroutine wd_val_range_3D
  end subroutine m_Parallel_init_mpi_nval

  subroutine m_Parallel_alloc_mpi_nval(nval)
    integer, intent(in) :: nval
    allocate(nis_nval(0:nrank_nval-1));ista_nval=nval+1
    allocate(nie_nval(0:nrank_nval-1));iend_nval=0
    allocate(nel_nval(0:nrank_nval-1));nel_nval=0
    allocate(map_nval(nval));map_nval=0
    allocate(map_z_nval(nval));map_z_nval=0
  end subroutine m_Parallel_alloc_mpi_nval

  subroutine m_Parallel_dealloc_mpi_nval()
    if(allocated(nis_nval)) deallocate(nis_nval)
    if(allocated(nie_nval)) deallocate(nie_nval)
    if(allocated(nel_nval)) deallocate(nel_nval)
    if(allocated(map_nval)) deallocate(map_nval)
    if(allocated(map_z_nval)) deallocate(map_z_nval)
    mpi_nval_enabled = .false.
  end subroutine m_Parallel_dealloc_mpi_nval

  subroutine m_Parallel_init_mpi_iba_3D_prev(kv3,iba_prev)
    integer, intent(in) :: kv3, iba_prev(kv3)

    logical             :: set_mapping_func
    integer             :: i, ik
    integer, allocatable, dimension(:,:) :: nis_g1k,nie_g1k,nel_g1k

    if(allocated(np_g1k_prev)) deallocate(np_g1k_prev);  allocate(np_g1k_prev(kv3))
    if(allocated(ista_g1k_prev)) deallocate(ista_g1k_prev); allocate(ista_g1k_prev(kv3))
    if(allocated(iend_g1k_prev)) deallocate(iend_g1k_prev); allocate(iend_g1k_prev(kv3))

    allocate(nis_g1k(0:nrank_g-1,kv3))
    do ik = 1, kv3
       nis_g1k(:,ik) = iba_prev(ik)
    end do

    allocate(nie_g1k(0:nrank_g-1,kv3))
    nie_g1k = 0

    allocate(nel_g1k(0:nrank_g-1,kv3))
    nel_g1k = 0

    do ik = 1, kv3
       set_mapping_func = .false.
       call set_block_range(iba_prev(ik),nrank_g,nel_g1k(0,ik),nis_g1k(0,ik) &
            & ,nie_g1k(0,ik), set_mapping_func)
       ista_g1k_prev(ik) = nis_g1k(myrank_g,ik)
       iend_g1k_prev(ik) = nie_g1k(myrank_g,ik)
       np_g1k_prev(ik)   = nel_g1k(myrank_g,ik)
    end do
    deallocate(nis_g1k)
    deallocate(nie_g1k)
    deallocate(nel_g1k)
  end subroutine m_Parallel_init_mpi_iba_3D_prev

  subroutine m_Parallel_init_mpi_iba_3D(nfout,ipri,printable,kv3,iba)
    integer, intent(in) :: nfout,ipri, kv3, iba(kv3)
    logical, intent(in) :: printable

    logical             :: set_mapping_func
    integer             :: i, ik
    integer, parameter  :: nc = 24
    integer, parameter  :: nc1 = 16
                                                  __TIMER_SUB_START(1242)

!!$#ifdef TRANSPOSE
    if(nrank_k > kv3) call phase_error_with_msg(nfout,' nrank_k > kv3 (m_Parallel_init_mpi_elec)',__LINE__,__FILE__)

! (( ista_g1k, iend_g1k ))
    if(allocated(np_g1k)) deallocate(np_g1k);  allocate(np_g1k(kv3))
    if(allocated(mp_g1k)) deallocate(mp_g1k);  allocate(mp_g1k(kv3))
    if(allocated(ista_g1k)) deallocate(ista_g1k); allocate(ista_g1k(kv3))
    if(allocated(iend_g1k)) deallocate(iend_g1k); allocate(iend_g1k(kv3))

    if(allocated(nis_g1k)) deallocate(nis_g1k);  allocate(nis_g1k(0:nrank_g-1,kv3))
    do ik = 1, kv3
       nis_g1k(:,ik) = iba(ik)
    end do

    if(allocated(nie_g1k)) deallocate(nie_g1k);  allocate(nie_g1k(0:nrank_g-1,kv3))
    nie_g1k = 0

    if(allocated(nel_g1k)) deallocate(nel_g1k);  allocate(nel_g1k(0:nrank_g-1,kv3))
    nel_g1k = 0

    do ik = 1, kv3
       set_mapping_func = .false.
       if(ipri >= 2 .and. printable) write(nfout,'(" !|| ik = ",i10 &
            & ," iba(",i3,") = ",i10," <<m_Parallel_init_mpi_iba_3D>>")') ik, ik, iba(ik)
       call set_block_range(iba(ik),nrank_g,nel_g1k(0,ik),nis_g1k(0,ik) &
            & ,nie_g1k(0,ik), set_mapping_func)
       ista_g1k(ik) = nis_g1k(myrank_g,ik)
       iend_g1k(ik) = nie_g1k(myrank_g,ik)
       np_g1k(ik)   = nel_g1k(myrank_g,ik)
       mp_g1k(ik)   = maxval(nel_g1k(0:nrank_g-1,ik))
       if(ipri >= 2 .and. printable) then
          write(nfout,'(" !|| --- mis_g1k,nie_g1k,nel_g1k ---")')
          write(nfout,'(" !|| ( rank_g  )",15i7)')(i,i=1,nrank_g)
          write(nfout,'(" !|| ( nis_g1k )",15i7)')(nis_g1k(i,ik),i=0,nrank_g-1)
          write(nfout,'(" !|| ( nie_g1k )",15i7)')(nie_g1k(i,ik),i=0,nrank_g-1)
          write(nfout,'(" !|| ( nel_g1k )",15i7)')(nel_g1k(i,ik),i=0,nrank_g-1)
          write(nfout,'(" !|| (myrank_g = ",i3," iba(ik),ista_g1k,iend_g1k,np_g1k,mp_g1k = ",5i7)') &
               & myrank_g,iba(ik),ista_g1k(ik),iend_g1k(ik),np_g1k(ik),mp_g1k(ik)
       end if

    end do
!!$#endif
                                                  __TIMER_SUB_STOP(1242)
  end subroutine m_Parallel_init_mpi_iba_3D

!===============================================================================

!!$#ifdef TRANSPOSE
  subroutine m_Parallel_init_mpi_nlmta_3D(nfout,ipri,nlmta,nlmt,natm,lmta,ntyp,ilmt,ityp)
    integer, intent(in) :: nfout, ipri, nlmta, nlmt,natm,ntyp
!!$ASASASASAS
!!$    integer, intent(in), dimension(nlmt,ntyp) :: lmta
    integer, intent(in), dimension(nlmt,natm) :: lmta
!!$ASASASASAS
    integer, intent(in), dimension(ntyp)      :: ilmt
    integer, intent(in), dimension(natm)      :: ityp

    integer :: ip, ia, it, it0, n1, i, nsum, lmt1
!!$ASASASASAS
    integer :: itmp
!!$ASASASASAS
    real(kind=8), allocatable, dimension(:) :: np
    real(kind=8) :: n_l, n_r
    integer :: ibuf_temp
    allocate(np(0:nrank_g-1)); np = 0.d0

    if(allocated(nis_fs)) deallocate(nis_fs)
    if(allocated(nie_fs)) deallocate(nie_fs)
    if(allocated(nel_fs)) deallocate(nel_fs)
    if(allocated(nis_fs_atm)) deallocate(nis_fs_atm)
    if(allocated(nie_fs_atm)) deallocate(nie_fs_atm)
    if(allocated(nel_fs_atm)) deallocate(nel_fs_atm)
    allocate(nis_fs(0:nrank_g-1))
    allocate(nie_fs(0:nrank_g-1))
    allocate(nel_fs(0:nrank_g-1))
    allocate(nis_fs_atm(0:nrank_g-1))
    allocate(nie_fs_atm(0:nrank_g-1))
    allocate(nel_fs_atm(0:nrank_g-1))

    do ip = 1, nrank_g
       np(ip-1) = dble(nlmta)/nrank_g *ip
       if(ipri >= 2) write(nfout,'(" !|| ip = ",i5," np = ",f12.4)') ip,np(ip-1)
    end do

    ip = 0
    do ia = 1, natm
       it = ityp(ia)
!!$ASASASASAS
       itmp=ilmt(it)
       if(itmp==0)cycle
!!$ASASASASAS
       n1 = lmta(ilmt(it),ia)
       if(abs(n1-np(ip)) < 1.d-10) then
          ip = ip + 1
          nie_fs(ip-1) = n1
          nie_fs_atm(ip-1) = ia
       else if(n1 > np(ip)) then
          ip = ip+1
          if(ia == 1) then
             n_l = np(ip-1)+1
          else if(ia >= 2) then
             it0 = ityp(ia-1)
             itmp=ilmt(it0)
!!$ASASASASAS
             if(itmp/=0) n_l = np(ip-1)-lmta(ilmt(it0),ia-1)
!!$             n_l = np(ip-1)-lmta(ilmt(it0),ia-1)
!!$ASASASASAS
          end if
          n_r = n1 - np(ip-1)
          if(n_l >= n_r) then
             nie_fs(ip-1) = n1
             nie_fs_atm(ip-1) = ia
          else
             nie_fs(ip-1) = np(ip-1) - n_l
!!$ASASASASAS
!!$             if(nie_fs(ip-1) < 0) nie_fs(ip-1) = n1
             if(nie_fs(ip-1) < 0) nie_fs(ip-1) = 0
!!$ASASASASAS
             nie_fs_atm(ip-1) = ia-1
          end if
       end if
    end do
    do i = ip+1, nrank_g
       nie_fs(i-1) = nlmta
       nie_fs_atm(i-1) = natm
    end do
    nis_fs_atm(0) = 1
    nis_fs(0) = 1
    nel_fs(0) = nie_fs(0)
    nel_fs_atm(0) = nie_fs_atm(0)
    do ip = 1, nrank_g-1
       nis_fs_atm(ip) = nie_fs_atm(ip-1)+1
       nis_fs(ip)  = nie_fs(ip-1)+1
       nel_fs(ip) = nie_fs(ip) - nie_fs(ip-1)
       nel_fs_atm(ip) = nie_fs_atm(ip) - nie_fs_atm(ip-1)
    end do

    ista_fs = nis_fs(myrank_g)
    iend_fs = nie_fs(myrank_g)
    np_fs   = nel_fs(myrank_g)
! == DEBUG by tkato 2011/09/23 =================================================
!   do i = 0, nrank_g-1
!      if(nel_fs(i) == 0) then
!         write(6,'(a,i3,a)') 'np_fs = 0 at nrank_g = ', i, '!!!'
!         call mpi_barrier(MPI_CommGroup, ierr)
!         call mpi_abort(MPI_CommGroup, -1, ierr)
!         stop
!      endif
!   enddo
! ==============================================================================
    mp_fs   = maxval(nel_fs)
    ista_fs_atm = nis_fs_atm(myrank_g)
    iend_fs_atm = nie_fs_atm(myrank_g)
    np_fs_atm = nel_fs_atm(myrank_g)
    mp_fs_atm = maxval(nel_fs_atm)

    if(ipri >= 2) then
       nsum = 0
       write(nfout,'(" !||",9x,"ip",4x,"nis_fs",4x,"nie_fs",7x,"nel",11x &
         & ,"atm(nis,nie)",5x,"lmta")')
       do ip = 1, nrank_g
          it = ityp(nie_fs_atm(ip-1))
          if(nie_fs_atm(ip-1) > 0 .and. nie_fs_atm(ip-1) <= natm ) then
             write(nfout,'(" !|| ",4i10,2x,3i10)') &
                  &   ip,nis_fs(ip-1),nie_fs(ip-1),nel_fs(ip-1),nis_fs_atm(ip-1), nie_fs_atm(ip-1) &
                  & , lmta(ilmt(it),nie_fs_atm(ip-1))
          else
             write(nfout,'(" !|| ",4i10,2x,2i10)') &
                  &   ip,nis_fs(ip-1),nie_fs(ip-1),nel_fs(ip-1),nis_fs_atm(ip-1), nie_fs_atm(ip-1)
          end if

          nsum = nsum + nel_fs(ip-1)
       end do
       write(nfout,'(" !|| nsum = ", i5)') nsum

       write(nfout,'(" !|| --- nis_fs,nie_fs,nel_fs ---")')
       write(nfout,'(" !|| ( rank_g  )",15i7)')(i,i=1,nrank_g)
       write(nfout,'(" !|| ( nis_fs  )",15i7)')(nis_fs(i),i=0,nrank_g-1)
       write(nfout,'(" !|| ( nie_fs  )",15i7)')(nie_fs(i),i=0,nrank_g-1)
       write(nfout,'(" !|| ( nel_fs  )",15i7)')(nel_fs(i),i=0,nrank_g-1)
       write(nfout,'(" !|| myrank_g = ",i3 &
            & ," nlmta,ista_fs,iend_fs,np_fs,mp_fs = ",/,10x,5i7)') &
            & myrank_g,nlmta,ista_fs,iend_fs,np_fs,mp_fs

       write(nfout,'(" !|| --- nis_fs_atm,nie_fs_atm,nel_fs_atm ---")')
       write(nfout,'(" !|| ( rank_g  )",15i7)')(i,i=1,nrank_g)
       write(nfout,'(" !|| ( nis_fs_atm )",15i7)')(nis_fs_atm(i),i=0,nrank_g-1)
       write(nfout,'(" !|| ( nie_fs_atm )",15i7)')(nie_fs_atm(i),i=0,nrank_g-1)
       write(nfout,'(" !|| ( nel_fs_atm )",15i7)')(nel_fs_atm(i),i=0,nrank_g-1)
       write(nfout,'(" !|| myrank_g = ",i3 &
            & ," nlmta,ista_fs_atm,iend_fs_atm,np_fs_atm,mp_fs_atm = ",/,10x,5i7)') &
            & myrank_g,nlmta,ista_fs_atm,iend_fs_atm,np_fs_atm,mp_fs_atm
    end if

! for zaj_l_ball fs(ri)_l_ball
  if(allocated(ball_buff)) deallocate(ball_buff)
  if(allocated(ball_addr)) deallocate(ball_addr)
  allocate(ball_buff(0:nrank_e-1))
  allocate(ball_addr(0:nrank_e-1))
  ibuf_temp = np_e*np_fs
  call mpi_allgather(ibuf_temp,1,mpi_integer,ball_buff(0),1,mpi_integer,mpi_kg_world,ierr)
  ball_addr(0) = 0
  do i = 1,nrank_e-1
     ball_addr(i) = ball_addr(i-1) + ball_buff(i-1)
  enddo


  deallocate(np)
  ista_lmta = HUGE(0)
  iend_lmta = -HUGE(0)
  do ia=ista_atm,iend_atm
     do lmt1=1,ilmt(ityp(ia))
        if(lmta(lmt1,ia)>iend_lmta) iend_lmta = lmta(lmt1,ia)
        if(lmta(lmt1,ia)<ista_lmta) ista_lmta = lmta(lmt1,ia)
     enddo
  enddo
  if(allocated(lmta_atm)) deallocate(lmta_atm)
  allocate(lmta_atm(ista_lmta:iend_lmta));lmta_atm=0
  do ia=ista_atm,iend_atm
     do lmt1=1,ilmt(ityp(ia))
         lmta_atm(lmta(lmt1,ia)) = ia
     enddo
  enddo
  end subroutine m_Parallel_init_mpi_nlmta_3D
!!$#endif

!===============================================================================

  subroutine m_Parallel_init_mpi_kngp_3D(nfout,ipri,kngp,comm_for_chg)
    integer, intent(in) :: nfout,ipri,kngp
    logical, intent(in) :: comm_for_chg
    integer :: iwork, i

    if(comm_for_chg)then
      mpi_chg_world = MPI_CommGroup
      nrank_chg = npes
      myrank_chg = mype
    else
      mpi_chg_world = mpi_ke_world
      nrank_chg = nrank_g
      myrank_chg = myrank_g
    endif
    if(allocated(is_kngp)) deallocate(is_kngp)
    if(allocated(ie_kngp)) deallocate(ie_kngp)
    if(allocated(nel_kngp)) deallocate(nel_kngp)
    allocate(is_kngp(0:nrank_chg-1))
    allocate(ie_kngp(0:nrank_chg-1))
    allocate(nel_kngp(0:nrank_chg-1))
    if(.not.comm_for_chg) then
      if(allocated(is_kngp_gw)) deallocate(is_kngp_gw)
      if(allocated(ie_kngp_gw)) deallocate(ie_kngp_gw)
      if(allocated(nel_kngp_gw)) deallocate(nel_kngp_gw)
      allocate(is_kngp_gw(0:nrank_chg-1))
      allocate(ie_kngp_gw(0:nrank_chg-1))
      allocate(nel_kngp_gw(0:nrank_chg-1))
    endif
    iwork = ( kngp - 1 ) / nrank_chg + 1
    if(ipri >= 2) then
       write(nfout,'(" << m_Parallel_init_mpi_kngp_3D >>")')
       write(nfout,'(" !|| kngp = ",i12)') kngp
       write(nfout,'(" !|| -- is_kngp(ista_kngp),  ie_kngp(iend_kngp) --")')
       write(nfout,'("mpi_chg_world=",i12)') mpi_chg_world
    end if
    do i = 0, nrank_chg-1
       is_kngp(i) = min(i*iwork+1, kngp+1)
       ie_kngp(i) = min(is_kngp(i)+iwork-1, kngp)
       nel_kngp(i) = ie_kngp(i) - is_kngp(i) + 1
       if(ipri >= 2) write(nfout,'(" !|| ",2i10)') is_kngp(i),ie_kngp(i)
    enddo
    ista_kngp = is_kngp(myrank_chg)
    iend_kngp = ie_kngp(myrank_chg)
    np_kngp   = nel_kngp(myrank_chg)
    mp_kngp   = maxval(nel_kngp)
    if(.not.comm_for_chg)then
      ista_kngp_gw = ista_kngp
      iend_kngp_gw = iend_kngp
      np_kngp_gw = np_kngp
      mp_kngp_gw = mp_kngp
      is_kngp_gw = is_kngp
      ie_kngp_gw = ie_kngp
      nel_kngp_gw = nel_kngp
    endif
    if(comm_for_chg) call mpi_kngp_3D_gw()

    contains

    subroutine mpi_kngp_3D_gw()
      integer :: nproc_per_g
      nproc_per_g = npes/nrank_g
      if(.not.allocated(is_kngp_gw)) allocate(is_kngp_gw(0:nrank_g-1))
      if(.not.allocated(ie_kngp_gw)) allocate(ie_kngp_gw(0:nrank_g-1))
      if(.not.allocated(nel_kngp_gw)) allocate(nel_kngp_gw(0:nrank_g-1))
      iwork = ( kngp - 1 ) / nrank_g + 1
      if(ipri >= 2) then
         write(nfout,'(" << m_Parallel_init_mpi_kngp_3D >>")')
         write(nfout,'(" !|| kngp = ",i12)') kngp
         write(nfout,'(" !|| -- is_kngp_gw(ista_kngp),  ie_kngp_gw(iend_kngp) --")')
      end if
      do i = 0, nrank_g-1
         is_kngp_gw(i) = is_kngp(nproc_per_g*i)
         ie_kngp_gw(i) = ie_kngp(nproc_per_g*i+(nproc_per_g-1))
         !is_kngp_gw(i) = min(i*iwork+1, kngp+1)
         !ie_kngp_gw(i) = min(is_kngp_gw(i)+iwork-1, kngp)
         nel_kngp_gw(i) = ie_kngp_gw(i) - is_kngp_gw(i) + 1
         if(ipri >= 2) write(nfout,'(" !|| ",2i10)') is_kngp_gw(i),ie_kngp_gw(i)
      enddo
      ista_kngp_gw = is_kngp_gw(myrank_g)
      iend_kngp_gw = ie_kngp_gw(myrank_g)
      np_kngp_gw   = nel_kngp_gw(myrank_g)
      mp_kngp_gw   = maxval(nel_kngp_gw)
    end subroutine mpi_kngp_3D_gw

  end subroutine m_Parallel_init_mpi_kngp_3D

!===============================================================================

  subroutine m_Parallel_init_mpi_kngp_exx(nfout,ipri,nmax_G_hyb)
    integer, intent(in) :: nfout,ipri,nmax_G_hyb
    integer :: iwork, i
    integer :: npes, mype

    npes = nrank_g
    mype = myrank_g

    if(.not.allocated(is_kngp_exx)) allocate(is_kngp_exx(0:npes-1))
    if(.not.allocated(ie_kngp_exx)) allocate(ie_kngp_exx(0:npes-1))
    if(.not.allocated(nel_kngp_exx)) allocate(nel_kngp_exx(0:npes-1))
    iwork = ( nmax_G_hyb - 1 ) / npes + 1
    if(ipri >= 2) then
       write(nfout,'(" << m_Parallel_init_mpi_kngp_exx >>")')
       write(nfout,'(" !|| nmax_G_hyb = ",i12)') nmax_G_hyb
       write(nfout,'(" !|| -- is_kngp_exx(ista_kngp_exx),  ie_kngp_exx(iend_kngp_exx) --")')
    end if
    do i = 0, npes-1
       is_kngp_exx(i) = min(i*iwork+1, nmax_G_hyb+1)
       ie_kngp_exx(i) = min(is_kngp_exx(i)+iwork-1, nmax_G_hyb)
       nel_kngp_exx(i) = ie_kngp_exx(i) - is_kngp_exx(i) + 1
       if(ipri >= 2) write(nfout,'(" !|| ",2i10)') is_kngp_exx(i),ie_kngp_exx(i)
    enddo
    ista_kngp_exx = is_kngp_exx(mype)
    iend_kngp_exx = ie_kngp_exx(mype)
    np_kngp_exx   = nel_kngp_exx(mype)
    mp_kngp_exx   = maxval(nel_kngp_exx)
  end subroutine m_Parallel_init_mpi_kngp_exx

  subroutine m_Parallel_init_mpi_kngp_B_3D(nfout,ipri,kngp)
    integer, intent(in) :: nfout,ipri,kngp
    integer :: iwork, i
    integer :: npes, mype
                                                  __TIMER_SUB_START(1234)
    npes = nrank_e
    mype = myrank_e

    allocate(is_kngp_B(0:npes-1))
    allocate(ie_kngp_B(0:npes-1))
    allocate(nel_kngp_B(0:npes-1))
    iwork = ( np_kngp - 1 ) / npes + 1
    if(ipri >= 2) then
       write(nfout,'(" << m_Parallel_init_mpi_kngp_3D >>")')
       write(nfout,'(" !|| kngp = ",i12)') kngp
       write(nfout,'(" !|| -- is_kngp_B(ista_kngp_B),  ie_kngp_B(iend_kngp_B) --")')
    end if
    do i = 0, npes-1
       is_kngp_B(i) = min(i*iwork+1, np_kngp+1)
       ie_kngp_B(i) = min(is_kngp_B(i)+iwork-1, np_kngp)
       nel_kngp_B(i) = ie_kngp_B(i) - is_kngp_B(i) + 1
       if(ipri >= 2) write(nfout,'(" !|| ",2i10)') is_kngp_B(i),ie_kngp_B(i)
    enddo
    is_kngp_B(:) = is_kngp_B(:) + ista_kngp - 1
    ie_kngp_B(:) = ie_kngp_B(:) + ista_kngp - 1
    ista_kngp_B = is_kngp_B(mype)
    iend_kngp_B = ie_kngp_B(mype)
    np_kngp_B   = nel_kngp_B(mype)
    mp_kngp_B   = maxval(nel_kngp_B)
                                                  __TIMER_SUB_STOP(1234)
  end subroutine m_Parallel_init_mpi_kngp_B_3D

!===============================================================================

!!$  subroutine m_Parallel_init_mpi_mix_3D(nfout,ipri,printable,kgpm)
!!$    integer, intent(in) :: nfout,ipri,kgpm
!!$    logical, intent(in) :: printable
!!$    integer :: iwork, i, npes, mype
!!$#ifdef __TIMER_SUB__
!!$  call timer_sta(1241)
!!$#endif
!!$
!!$    npes = nrank_g
!!$    mype = myrank_g
!!$    allocate(is_kgpm(0:npes-1))
!!$    allocate(ie_kgpm(0:npes-1))
!!$    allocate(nel_kgpm(0:npes-1))
!!$    iwork = ( kgpm - 1 ) / npes + 1
!!$    if(ipri >= 1 .and. printable) then
!!$       write(nfout,'(" !|| << init_mpi_mix_3D >>")')
!!$       write(nfout,'(" !|| kgpm = ",i12)') kgpm
!!$       write(nfout,'(" !|| -- is_kgpm, ie_kgpm --")')
!!$    end if
!!$    do i = 0, npes-1
!!$       is_kgpm(i) = min(i*iwork+1, kgpm+1)
!!$       ie_kgpm(i) = min(is_kgpm(i)+iwork-1, kgpm)
!!$       nel_kgpm(i) = ie_kgpm(i) - is_kgpm(i) + 1
!!$       if(ipri >= 1 .and. printable) write(nfout,'(" !|| ",2i12)') is_kgpm(i),ie_kgpm(i)
!!$    enddo
!!$    ista_kgpm = is_kgpm(mype)
!!$    iend_kgpm = ie_kgpm(mype)
!!$    np_kgpm   = nel_kgpm(mype)
!!$    mp_kgpm   = maxval(nel_kgpm)
!!$#ifdef __TIMER_SUB__
!!$  call timer_end(1241)
!!$#endif
!!$  end subroutine m_Parallel_init_mpi_mix_3D
!!$!===============================================================================

!===============================================================================
  subroutine m_Parallel_mpi_fft_box_3div(nfout,ipri,printable,fft_box_size_WF,kimg, &
 &                                  fftbox_div_1,fftbox_div_2,fftbox_div_3)
    implicit none

    integer, intent(in) :: nfout, ipri, kimg, fftbox_div_1 ,fftbox_div_2, fftbox_div_3
    logical, intent(in) :: printable
    integer, intent(in) :: fft_box_size_WF(3,0:1)
    integer, allocatable, dimension(:,:) :: wk_x
    integer, allocatable, dimension(:,:,:) :: wk_fft_box
    integer :: div1, div2, div3, igf, igf1, igf2, igf3
    integer :: i, j, id, id1, id2, id3, ierr, ichkalloc
    integer :: ii, jj, kk, ix, iy, iz
    integer :: nfft, n, nmrank, myrank, nfftm2
    integer :: rank1, rank2, rank3, nfft1, nfft2, nfft3, rank
    integer :: kx1p, kx2p, kx3p, n1s, n2s, n3s, n1e, n2e, n3e, n1p, n2p, n3p
    integer :: key, color
    character(len=12) :: strng,strngmax
                                                  __TIMER_SUB_START(1232)
    logical,save :: firstcall=.true.
    nmrank = nrank_g
    myrank = myrank_g

    if ((fftbox_div_1 < 0) .or. (fftbox_div_2 < 0) .or. (fftbox_div_3 < 0)) then
      write(nfout,'("Either of the division numbers of FFT-BOX is negative.")')
      call flush(nfout)
      call mpi_abort(mpi_comm_world, 1234567 , ierr)
      call phase_error_with_msg(nfout,'Either of the division numbers of FFT-BOX is negative.',__LINE__,__FILE__)
    endif
    if ((fftbox_div_1*fftbox_div_2*fftbox_div_3) > nmrank ) then
      write(nfout,'("The three-dimensional division number of FFT-BOX is over the number of processes.")')
      call flush(nfout)
      call mpi_abort(mpi_comm_world, 1234568 , ierr)
      call phase_error_with_msg(nfout,'The three-dimensional division number of FFT-BOX is over the number of processes.'&
                               ,__LINE__,__FILE__)
    end if
    if ((fftbox_div_1*fftbox_div_2*fftbox_div_3) == 0 ) then
! === DEBUG by tkato 2013/09/18 ================================================
!     div1 = int(cbrt(real(nmrank)))
#ifdef __FUJITSU
      div1 = int(cbrt(real(nmrank)))
#else
      div1 = int((real(nmrank))**(1.0d0/3.0d0))
#endif
! ==============================================================================
      div2 = div1
      div3 = div1
    else
      div1 = fftbox_div_1
      div2 = fftbox_div_2
      div3 = fftbox_div_3
    end if

    nfft   =  fft_box_size_WF(1,0)*fft_box_size_WF(2,0)*fft_box_size_WF(3,0)
    nfftm2 =  fft_box_size_WF(1,1)*fft_box_size_WF(2,1)*fft_box_size_WF(3,1)
    nfft1 = fft_box_size_WF(1,1)
    nfft2 = fft_box_size_WF(2,1)
    nfft3 = fft_box_size_WF(3,1)

    if(mod(nfft1,div1) > 0) then
      kx1p = nfft1/div1+1
    else
      kx1p = nfft1/div1
    end if
    if(mod(nfft2,div2) > 0) then
      kx2p = nfft2/div2+1
    else
      kx2p = nfft2/div2
    end if
    if(mod(nfft3,div3) > 0) then
      kx3p = nfft3/div3+1
    else
      kx3p = nfft3/div3
    end if

    allocate(wk_x( 6,nmrank ), stat=ichkalloc)
    if(ichkalloc /= 0) then
       write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 1'
       call mpi_abort(mpi_comm_world, 3 , ichkalloc)
    endif
    wk_x(:,:) = 0

    fft_X_x_dim = div1
    fft_X_y_dim = div2
    fft_X_z_dim = div3
    fft_X_x_nel = kx1p
    fft_X_y_nel = kx2p
    fft_X_z_nel = kx3p

    allocate(nis_fft_X_z(fft_X_z_dim), stat=ichkalloc)
    allocate(nie_fft_X_z(fft_X_z_dim), stat=ichkalloc)
    allocate(nis_fft_X_y(fft_X_y_dim), stat=ichkalloc)
    allocate(nie_fft_X_y(fft_X_y_dim), stat=ichkalloc)
    allocate(nis_fft_X_x(fft_X_x_dim), stat=ichkalloc)
    allocate(nie_fft_X_x(fft_X_x_dim), stat=ichkalloc)

    do rank = 0, div1*div2*div3-1
      rank1 = mod(rank,div1)
      rank2 = mod(rank/div1,div2)
      rank3 = rank/(div1*div2)

      n1s = kx1p*rank1+1
      n1e = min(nfft1,kx1p*(rank1+1))
      n1p = max(0,n1e-n1s+1)
      n2s = kx2p*rank2+1
      n2e = min(nfft2,kx2p*(rank2+1))
      n2p = max(0,n2e-n2s+1)
      n3s = kx3p*rank3+1
      n3e = min(nfft3,kx3p*(rank3+1))
      n3p = max(0,n3e-n3s+1)

      wk_x(1,rank+1) = n1s
      wk_x(2,rank+1) = n1e
      wk_x(3,rank+1) = n2s
      wk_x(4,rank+1) = n2e
      wk_x(5,rank+1) = n3s
      wk_x(6,rank+1) = n3e

      nis_fft_X_x(rank1+1) = n1s
      nie_fft_X_x(rank1+1) = n1e
      nis_fft_X_y(rank2+1) = n2s
      nie_fft_X_y(rank2+1) = n2e
      nis_fft_X_z(rank3+1) = n3s
      nie_fft_X_z(rank3+1) = n3e
    end do

!!!!
!!  n = 1
!!!!
    n = 0
    if (wk_x(1,myrank+1) /= 0) then
      xyz_fft_x(1,1) = wk_x(1,myrank+1)
      xyz_fft_x(2,1) = wk_x(2,myrank+1)
      xyz_fft_x(1,2) = wk_x(3,myrank+1)
      xyz_fft_x(2,2) = wk_x(4,myrank+1)
      xyz_fft_x(1,3) = wk_x(5,myrank+1)
      xyz_fft_x(2,3) = wk_x(6,myrank+1)
    endif

    if(allocated(map_fft_x)) deallocate(map_fft_x)
    allocate(map_fft_x( nfft ), stat=ichkalloc)
    if(ichkalloc /= 0) then
       write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 2'
       call mpi_abort(mpi_comm_world, 4 , ichkalloc)
    endif
    map_fft_x(:) = 0
    np_fft_x = 0

    id1 = fft_box_size_WF(1,0)
    id2 = fft_box_size_WF(2,0)
    id3 = fft_box_size_WF(3,0)

    allocate( wk_fft_box(id1, id2, id3) )
    do igf3 = 1, fft_box_size_WF(3,n)
       do igf2 = 1, fft_box_size_WF(2,n)
          do igf1 = 1, fft_box_size_WF(1,n)
             igf = igf1 + (igf2-1)*id1 + (igf3-1)*id1*id2

             wk_fft_box(igf1,igf2,igf3) = igf

          enddo
       enddo
    enddo

    do rank = 1, div1*div2*div3
       if (wk_x(1,rank) == 0) cycle
       do iz = wk_x(5,rank) , wk_x(6,rank)
          do iy = wk_x(3,rank) , wk_x(4,rank)
             do ix = wk_x(1,rank) , wk_x(2,rank)
                map_fft_x(id1*id2*(iz-1)+id1*(iy-1)+ix) = rank
             end do
          end do
       end do
    end do

    if ((div1*div2*div3) > myrank) then
       if (wk_x(1,myrank+1) /= 0) then
          np_fft_x = (wk_x(2,myrank+1)-wk_x(1,myrank+1)+1)* &
         &           (wk_x(4,myrank+1)-wk_x(3,myrank+1)+1)* &
         &           (wk_x(6,myrank+1)-wk_x(5,myrank+1)+1)
          if(allocated(mp_fft_x)) deallocate(mp_fft_x)
          allocate(mp_fft_x(np_fft_x), stat=ichkalloc)
       end if
       if(ichkalloc /= 0) then
          write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 3'
          call mpi_abort(mpi_comm_world, 5 , ichkalloc)
       endif
    else
       np_fft_x = 0
    end if

    deallocate(wk_x)

       if (np_fft_x /= 0) then
          id = 0
          do kk = xyz_fft_x(1,3), xyz_fft_x(2,3)
            do jj = xyz_fft_x(1,2), xyz_fft_x(2,2)
              do ii = xyz_fft_x(1,1), xyz_fft_x(2,1)
                id = id + 1
                mp_fft_x(id) = wk_fft_box(ii,jj,kk)
              enddo
            enddo
          enddo
       end if

    deallocate(wk_fft_box)

    if(allocated(nel_fft_x)) deallocate(nel_fft_x)
    allocate(nel_fft_x( 0:nmrank-1 ), stat=ichkalloc)
    if(ichkalloc /= 0) then
       write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 4'
       call mpi_abort(mpi_comm_world, 6 , ichkalloc)
    endif
    nel_fft_x(:) = 0

    do i = 1, nfft
      if (map_fft_x(i) .gt. 0) then
         nel_fft_x( map_fft_x(i)-1 ) = nel_fft_x( map_fft_x(i)-1 ) + 1
      endif
    enddo

    if (ipri > 1) then
      write(nfout,'("|||| m_Parallel_mpi_fft_box")')
      write(nfout,'("|||| div1 = ",i5," , div2 = ",i5," , div3 = ",i5)')   div1,div2,div3
      write(nfout,'("|||| fft-box-size (",i4," ,",i4," ,",i4,")")') &
     &          ((fft_box_size_WF(i,j),i=1,3),j=0,0)
      write(nfout,'("|||| fft-box-div[X-axis] (",i4," :",i4," ,",i4," :",i4," ,",i4," :",i4,")")') &
     &          ((xyz_fft_x(i,j),i=1,2),j=1,3)
!     write(nfout,'("|||| ir_x2z=",i4," , ir_z2y=",i4," , ir_y2z=",i4," , ir_z2x=",i4)') &
!    &                    ir_x2z,ir_z2y,ir_y2z,ir_z2x
!     write(nfout,'("|||| is_x2z=",i4," , is_z2y=",i4," , is_y2z=",i4," , is_z2x=",i4)') &
!    &                    is_x2z,is_z2y,is_y2z,is_z2x
      write(nfout,'("|||| np_fft_x=",i8)')          np_fft_x
      write(nfout,'("|||| nel_fft_x:",/,"||||",10(i7))') (nel_fft_x(i),i=0,nmrank-1)
      write(nfout,'("|||| fft_X_x_dim=",i4,",  fft_X_y_dim=",i4,",  fft_X_z_dim=",i4)') &
     &                    fft_X_x_dim,fft_X_y_dim, fft_X_z_dim
      write(nfout, '("|||| nis_fft_X_x :")')
      write(nfout,'(("||||",10(i8)))') (nis_fft_X_x(i),i=1,fft_X_x_dim)
      write(nfout, '("|||| nie_fft_X_x :")')
      write(nfout,'(("||||",10(i8)))') (nie_fft_X_x(i),i=1,fft_X_x_dim)
      write(nfout, '("|||| nis_fft_X_y :")')
      write(nfout,'(("||||",10(i8)))') (nis_fft_X_y(i),i=1,fft_X_y_dim)
      write(nfout, '("|||| nie_fft_X_y :")')
      write(nfout,'(("||||",10(i8)))') (nie_fft_X_y(i),i=1,fft_X_y_dim)
      write(nfout, '("|||| nis_fft_X_z :")')
      write(nfout,'(("||||",10(i8)))') (nis_fft_X_z(i),i=1,fft_X_z_dim)
      write(nfout, '("|||| nie_fft_X_z :")')
      write(nfout,'(("||||",10(i8)))') (nie_fft_X_z(i),i=1,fft_X_z_dim)
    endif

    color = nrank_g
!    if (myrank_g > (div1*div2*div3-1)) then
!       color = MPI_UNDEFINED
!    end if
    if (nrank_g>(div1*div2*div3)) then
      write(strng,'(i10)') nrank_g
      write(strngmax,'(i10)') div1*div2*div3
      call phase_error_with_msg(nfout,'ng = '//trim(strng)//', max ng ='//trim(strngmax), &
      __LINE__,__FILE__)
    endif
    key = myrank_g
!    if(firstcall) then
    if(.not. firstcall) call mpi_comm_free(mpi_fft_zy_world, ierr)
    call mpi_comm_split(mpi_ke_world,color,key,mpi_fft_zy_world,ierr)
    if (ierr /= 0) then
       write(nfout,*)' m_Parallel_mpi_fft_box :  mpi_comm_split error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 100003, ierr)
    endif
!    endif
    firstcall = .false.
                                                  __TIMER_SUB_STOP(1232)
  end subroutine m_Parallel_mpi_fft_box_3div
!===============================================================================
!===============================================================================
  subroutine m_Parallel_mpi_fft_box_xyz(nfout,ipri,printable,fft_box_size_WF,kimg,fft_div1,fft_div2)
    implicit none

    integer, intent(in) :: nfout,ipri, kimg, fft_div1, fft_div2
    logical, intent(in) :: printable
    integer, intent(in) :: fft_box_size_WF(3,0:1)
    integer, allocatable, dimension(:,:) :: wk_z, wk_x, wk_y, wk1, wk2
    integer, allocatable, dimension(:,:,:) :: wk_fft_box
    integer :: dim1, dim2, igf, igf1, igf2, igf3
    integer :: i, j, id, id1, id2, id3, ierr, ichkalloc
    integer :: ii, jj, kk, ix, iy, iz
    integer :: nfft, n, nmrank, myrank
    integer :: key, color
    integer :: mindim1, mindim2
    logical,save :: firstcall = .true.
    character(len=12) :: strng,strngmax
                                                  __TIMER_SUB_START(1232)
    nmrank = nrank_g
    myrank = myrank_g

    mindim1 = 0
    mindim2 = 0
    if (((fft_div1*fft_div2) /= 0) .and. ((fft_div1*fft_div2) <= nmrank)) then
       dim1 = fft_div1
       dim2 = fft_div2
    else
       if (kimg == 1) then
          mindim1 = min(fft_box_size_WF(1,0)/2,fft_box_size_WF(2,0))
       else
          mindim1 = min(fft_box_size_WF(1,0),fft_box_size_WF(2,0))
       end if
       mindim2 = min(fft_box_size_WF(2,0),fft_box_size_WF(3,0))

       call decide_div_fft(nmrank, mindim1, mindim2, dim1, dim2)
    end if

    nfft = fft_box_size_WF(1,0)*fft_box_size_WF(2,0)*fft_box_size_WF(3,0)

    allocate(wk_z( 4,nmrank ), stat=ichkalloc)
    allocate(wk_x( 4,nmrank ), stat=ichkalloc)
    allocate(wk_y( 4,nmrank ), stat=ichkalloc)
    if(ichkalloc /= 0) then
       write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 1'
       call mpi_abort(mpi_comm_world, 3 , ichkalloc)
    endif
    wk_z(:,:) = 0
    wk_x(:,:) = 0
    wk_y(:,:) = 0

    allocate(wk1(2,dim1))
    allocate(wk2(2,dim2))
    wk1 = 0
    wk2 = 0

    fft_X_y_dim = dim1
    fft_X_z_dim = dim2
    fft_Y_x_dim = dim1
    fft_Y_z_dim = dim2
    fft_Z_x_dim = dim1
    fft_Z_y_dim = dim2
    allocate(nis_fft_X_z(fft_X_z_dim), stat=ichkalloc)
    allocate(nie_fft_X_z(fft_X_z_dim), stat=ichkalloc)
    allocate(nis_fft_X_y(fft_X_y_dim), stat=ichkalloc)
    allocate(nie_fft_X_y(fft_X_y_dim), stat=ichkalloc)
    allocate(nis_fft_Z_x(fft_Z_x_dim), stat=ichkalloc)
    allocate(nie_fft_Z_x(fft_Z_x_dim), stat=ichkalloc)
    allocate(nis_fft_Z_y(fft_Z_y_dim), stat=ichkalloc)
    allocate(nie_fft_Z_y(fft_Z_y_dim), stat=ichkalloc)
    allocate(nis_fft_Y_x(fft_Y_x_dim), stat=ichkalloc)
    allocate(nie_fft_Y_x(fft_Y_x_dim), stat=ichkalloc)
    allocate(nis_fft_Y_z(fft_Y_z_dim), stat=ichkalloc)
    allocate(nie_fft_Y_z(fft_Y_z_dim), stat=ichkalloc)

    call index_div_fft(kimg, wk1, dim1, fft_box_size_WF(2,0), 0)
    call index_div_fft(kimg, wk2, dim2, fft_box_size_WF(3,0), 0)
    nis_fft_X_y(:) = wk1(1,:)
    nie_fft_X_y(:) = wk1(2,:)
    nis_fft_X_z(:) = wk2(1,:)
    nie_fft_X_z(:) = wk2(2,:)
    ii = 1
    jj = 1
    do i = 1, dim1*dim2
       wk_x(1,i) = wk1(1,ii)
       wk_x(2,i) = wk1(2,ii)
       wk_x(3,i) = wk2(1,jj)
       wk_x(4,i) = wk2(2,jj)
       ii = ii + 1
       if (ii > dim1) then
          ii = 1
          jj = jj + 1
       end if
    end do

     call index_div_fft(kimg, wk1, dim1, fft_box_size_WF(1,0), 1)
    nis_fft_Y_x(:) = wk1(1,:)
    nie_fft_Y_x(:) = wk1(2,:)
    nis_fft_Y_z(:) = wk2(1,:)
    nie_fft_Y_z(:) = wk2(2,:)
    ii = 1
    jj = 1
    do i = 1, dim1*dim2
       wk_y(1,i) = wk1(1,ii)
       wk_y(2,i) = wk1(2,ii)
       wk_y(3,i) = wk2(1,jj)
       wk_y(4,i) = wk2(2,jj)
       ii = ii + 1
       if (ii > dim1) then
          ii = 1
          jj = jj + 1
       end if
    end do

    call index_div_fft(kimg, wk2, dim2, fft_box_size_WF(2,0), 0)
    nis_fft_Z_x(:) = wk1(1,:)
    nie_fft_Z_x(:) = wk1(2,:)
    nis_fft_Z_y(:) = wk2(1,:)
    nie_fft_Z_y(:) = wk2(2,:)
    ii = 1
    jj = 1
    do i = 1, dim1*dim2
       wk_z(1,i) = wk1(1,ii)
       wk_z(2,i) = wk1(2,ii)
       wk_z(3,i) = wk2(1,jj)
       wk_z(4,i) = wk2(2,jj)
       ii = ii + 1
       if (ii > dim1) then
          ii = 1
          jj = jj + 1
       end if
    end do

!!!!
!!  n = 1
!!!!
    n = 0
    if (wk_x(1,myrank+1) /= 0) then
      xyz_fft_x(1,1) = 1
      xyz_fft_x(2,1) = fft_box_size_WF(1,n)
      xyz_fft_x(1,2) = wk_x(1,myrank+1)
      xyz_fft_x(2,2) = wk_x(2,myrank+1)
      xyz_fft_x(1,3) = wk_x(3,myrank+1)
      xyz_fft_x(2,3) = wk_x(4,myrank+1)
    endif

    if (wk_y(1,myrank+1) /= 0) then
      xyz_fft_y(1,1) = wk_y(1,myrank+1)
      xyz_fft_y(2,1) = wk_y(2,myrank+1)
      xyz_fft_y(1,2) = 1
      xyz_fft_y(2,2) = fft_box_size_WF(2,n)
      xyz_fft_y(1,3) = wk_y(3,myrank+1)
      xyz_fft_y(2,3) = wk_y(4,myrank+1)
    endif

    if (wk_z(1,myrank+1) /= 0) then
      xyz_fft_z(1,1) = wk_z(1,myrank+1)
      xyz_fft_z(2,1) = wk_z(2,myrank+1)
      xyz_fft_z(1,2) = wk_z(3,myrank+1)
      xyz_fft_z(2,2) = wk_z(4,myrank+1)
      xyz_fft_z(1,3) = 1
      xyz_fft_z(2,3) = fft_box_size_WF(3,n)
    endif

    if(allocated(map_fft_x)) deallocate(map_fft_x)
    if(allocated(map_fft_z)) deallocate(map_fft_z)
    if(allocated(map_fft_y)) deallocate(map_fft_y)
    allocate(map_fft_x( nfft ), stat=ichkalloc)
    allocate(map_fft_z( nfft ), stat=ichkalloc)
    allocate(map_fft_y( nfft ), stat=ichkalloc)
    if(ichkalloc /= 0) then
       write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 2'
       call mpi_abort(mpi_comm_world, 4 , ichkalloc)
    endif
    map_fft_x(:) = 0
    map_fft_z(:) = 0
    map_fft_y(:) = 0

    np_fft_x = 0
    np_fft_z = 0
    np_fft_y = 0

    id1 = fft_box_size_WF(1,0)
    id2 = fft_box_size_WF(2,0)
    id3 = fft_box_size_WF(3,0)

    allocate( wk_fft_box(fft_box_size_WF(1,0), fft_box_size_WF(2,0), fft_box_size_WF(3,0)))
    do igf3 = 1, fft_box_size_WF(3,n)
       do igf2 = 1, fft_box_size_WF(2,n)
          do igf1 = 1, fft_box_size_WF(1,n)
             igf = igf1 + (igf2-1)*id1 + (igf3-1)*id1*id2
             wk_fft_box(igf1,igf2,igf3) = igf
          enddo
       enddo
    enddo

    do ii = 1, dim1*dim2
       if (wk_x(1,ii) == 0) cycle
       do iz = wk_x(3,ii) , wk_x(4,ii)
          do iy = wk_x(1,ii) , wk_x(2,ii)
             do ix = 1, id1
                map_fft_x(id1*id2*(iz-1)+id1*(iy-1)+ix) = ii
             end do
          end do
       end do
    end do

    do ii = 1, dim1*dim2
       if (wk_z(1,ii) == 0) cycle
       do iz = 1, id3
          do iy = wk_z(3,ii) , wk_z(4,ii)
             do ix = wk_z(1,ii), wk_z(2,ii)
                map_fft_z(id1*id2*(iz-1)+id1*(iy-1)+ix) = ii
             end do
          end do
       end do
    end do

    do ii = 1, dim1*dim2
       if (wk_y(1,ii) == 0) cycle
       do iz = wk_y(3,ii) , wk_y(4,ii)
          do iy = 1, id2
             do ix = wk_y(1,ii), wk_y(2,ii)
                map_fft_y(id1*id2*(iz-1)+id1*(iy-1)+ix) = ii
             end do
          end do
       end do
    end do

    if ((dim1*dim2) > myrank) then
       if (wk_x(1,myrank+1) /= 0) then
          np_fft_x = id1*(wk_x(2,myrank+1)-wk_x(1,myrank+1)+1)*(wk_x(4,myrank+1)-wk_x(3,myrank+1)+1)
          if(allocated(mp_fft_x)) deallocate(mp_fft_x)
          allocate(mp_fft_x(np_fft_x), stat=ichkalloc)
       end if
       if (wk_z(1,myrank+1) /= 0) then
          np_fft_z = (wk_z(4,myrank+1)-wk_z(3,myrank+1)+1)*(wk_z(2,myrank+1)-wk_z(1,myrank+1)+1)*id3
          if(allocated(mp_fft_z)) deallocate(mp_fft_z)
          allocate(mp_fft_z(np_fft_z), stat=ichkalloc)
       end if
       if (wk_y(1,myrank+1) /= 0) then
          np_fft_y = (wk_y(4,myrank+1)-wk_y(3,myrank+1)+1)*id2*(wk_y(2,myrank+1)-wk_y(1,myrank+1)+1)
          if(allocated(mp_fft_y)) deallocate(mp_fft_y)
          allocate(mp_fft_y(np_fft_y), stat=ichkalloc)
       end if
       if(ichkalloc /= 0) then
          write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 3'
          call mpi_abort(mpi_comm_world, 5 , ichkalloc)
       endif
    else
       np_fft_x = 0
       np_fft_z = 0
       np_fft_y = 0
    end if

    deallocate(wk_x)
    deallocate(wk_z)
    deallocate(wk_y)

    if (kimg == 1) then
       if (np_fft_x /= 0) then
          id = 0
          do kk = xyz_fft_x(1,3), xyz_fft_x(2,3)
            do jj = xyz_fft_x(1,2), xyz_fft_x(2,2)
              do ii = xyz_fft_x(1,1), xyz_fft_x(2,1) , 2
                id = id + 1
                mp_fft_x(id) = wk_fft_box(ii,jj,kk)
                id = id + 1
                mp_fft_x(id) = wk_fft_box(ii+1,jj,kk)
              enddo
            enddo
          enddo
       end if
       if (np_fft_y /= 0) then
          id = 0
            do jj = xyz_fft_y(1,2), xyz_fft_y(2,2)
          do kk = xyz_fft_y(1,3), xyz_fft_y(2,3)
              do ii = xyz_fft_y(1,1), xyz_fft_y(2,1) , 2
                id = id + 1
                mp_fft_y(id) = wk_fft_box(ii,jj,kk)
                id = id + 1
                mp_fft_y(id) = wk_fft_box(ii+1,jj,kk)
              enddo
            enddo
          enddo
       end if
       if (np_fft_z /= 0) then
          id = 0
          do kk = xyz_fft_z(1,3), xyz_fft_z(2,3)
            do jj = xyz_fft_z(1,2), xyz_fft_z(2,2)
              do ii = xyz_fft_z(1,1), xyz_fft_z(2,1) , 2
                id = id + 1
                mp_fft_z(id) = wk_fft_box(ii,jj,kk)
                id = id + 1
                mp_fft_z(id) = wk_fft_box(ii+1,jj,kk)
              enddo
            enddo
          enddo
       end if
    else
       if (np_fft_x /= 0) then
          id = 0
          do kk = xyz_fft_x(1,3), xyz_fft_x(2,3)
            do jj = xyz_fft_x(1,2), xyz_fft_x(2,2)
              do ii = xyz_fft_x(1,1), xyz_fft_x(2,1)
                id = id + 1
                mp_fft_x(id) = wk_fft_box(ii,jj,kk)
              enddo
            enddo
          enddo
       end if
       if (np_fft_y /= 0) then
          id = 0
            do jj = xyz_fft_y(1,2), xyz_fft_y(2,2)
          do kk = xyz_fft_y(1,3), xyz_fft_y(2,3)
              do ii = xyz_fft_y(1,1), xyz_fft_y(2,1)
                id = id + 1
                mp_fft_y(id) = wk_fft_box(ii,jj,kk)
              enddo
            enddo
          enddo
       end if
       if (np_fft_z /= 0) then
          id = 0
          do kk = xyz_fft_z(1,3), xyz_fft_z(2,3)
            do jj = xyz_fft_z(1,2), xyz_fft_z(2,2)
              do ii = xyz_fft_z(1,1), xyz_fft_z(2,1)
                id = id + 1
                mp_fft_z(id) = wk_fft_box(ii,jj,kk)
              enddo
            enddo
          enddo
       end if
    endif

    deallocate(wk_fft_box)

    if(allocated(nel_fft_x)) deallocate(nel_fft_x)
    if(allocated(nel_fft_z)) deallocate(nel_fft_z)
    if(allocated(nel_fft_y)) deallocate(nel_fft_y)
    allocate(nel_fft_x( 0:nmrank-1 ), stat=ichkalloc)
    allocate(nel_fft_z( 0:nmrank-1 ), stat=ichkalloc)
    allocate(nel_fft_y( 0:nmrank-1 ), stat=ichkalloc)
    if(ichkalloc /= 0) then
       write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 4'
       call mpi_abort(mpi_comm_world, 6 , ichkalloc)
    endif
    nel_fft_z(:) = 0
    nel_fft_x(:) = 0
    nel_fft_y(:) = 0

    do i = 1, nfft
      if (map_fft_x(i) .gt. 0) then
         nel_fft_x( map_fft_x(i)-1 ) = nel_fft_x( map_fft_x(i)-1 ) + 1
      endif
      if (map_fft_z(i) .gt. 0) then
         nel_fft_z( map_fft_z(i)-1 ) = nel_fft_z( map_fft_z(i)-1 ) + 1
      endif
      if (map_fft_y(i) .gt. 0) then
         nel_fft_y( map_fft_y(i)-1 ) = nel_fft_y( map_fft_y(i)-1 ) + 1
      endif
    enddo

    if (ipri > 1) then
      write(nfout,'("|||| m_Parallel_mpi_fft_box")')
      write(nfout,'("|||| dim1 = ",i5," , dim2 = ",i5)')   dim1,dim2
      write(nfout,'("|||| fft-box-size (",i4," ,",i4," ,",i4,")")') &
     &          ((fft_box_size_WF(i,j),i=1,3),j=0,0)
      write(nfout,'("|||| fft-box-div[X-axis] (",i4," :",i4," ,",i4," :",i4," ,",i4," :",i4,")")') &
     &          ((xyz_fft_x(i,j),i=1,2),j=1,3)
      write(nfout,'("|||| fft-box-div[Y-axis] (",i4," :",i4," ,",i4," :",i4," ,",i4," :",i4,")")') &
     &          ((xyz_fft_y(i,j),i=1,2),j=1,3)
      write(nfout,'("|||| fft-box-div[Z-axis] (",i4," :",i4," ,",i4," :",i4," ,",i4," :",i4,")")') &
     &          ((xyz_fft_z(i,j),i=1,2),j=1,3)
!     write(nfout,'("|||| ir_x2z=",i4," , ir_z2y=",i4," , ir_y2z=",i4," , ir_z2x=",i4)') &
!    &                    ir_x2z,ir_z2y,ir_y2z,ir_z2x
!     write(nfout,'("|||| is_x2z=",i4," , is_z2y=",i4," , is_y2z=",i4," , is_z2x=",i4)') &
!    &                    is_x2z,is_z2y,is_y2z,is_z2x
      write(nfout,'("|||| np_fft_x=",i8)')          np_fft_x
      write(nfout,'("|||| nel_fft_x:",/,"||||",10(i7))') (nel_fft_x(i),i=0,nmrank-1)
      write(nfout,'("|||| np_fft_y=",i8)')          np_fft_y
      write(nfout,'("|||| nel_fft_y:",/,"||||",10(i7))') (nel_fft_y(i),i=0,nmrank-1)
      write(nfout,'("|||| np_fft_z=",i8)')          np_fft_z
      write(nfout,'("|||| nel_fft_z:",/,"||||",10(i7))') (nel_fft_z(i),i=0,nmrank-1)
      write(nfout, '("|||| fft_X_y_dim=",i4,",  fft_X_z_dim=",i4)') fft_X_y_dim, fft_X_z_dim
      write(nfout, '("|||| nis_fft_X_y :")')
      write(nfout,'(("||||",10(i8)))') (nis_fft_X_y(i),i=1,fft_X_y_dim)
      write(nfout, '("|||| nie_fft_X_y :")')
      write(nfout,'(("||||",10(i8)))') (nie_fft_X_y(i),i=1,fft_X_y_dim)
      write(nfout, '("|||| nis_fft_X_z :")')
      write(nfout,'(("||||",10(i8)))') (nis_fft_X_z(i),i=1,fft_X_z_dim)
      write(nfout, '("|||| nie_fft_X_z :")')
      write(nfout,'(("||||",10(i8)))') (nie_fft_X_z(i),i=1,fft_X_z_dim)
      write(nfout, '("|||| fft_Z_x_dim=",i4,",  fft_Z_y_dim=",i4)') fft_Z_x_dim, fft_Z_y_dim
      write(nfout, '("|||| nis_fft_Z_x :")')
      write(nfout,'(("||||",10(i8)))') (nis_fft_Z_x(i),i=1,fft_Z_x_dim)
      write(nfout, '("|||| nie_fft_Z_x :")')
      write(nfout,'(("||||",10(i8)))') (nie_fft_Z_x(i),i=1,fft_Z_x_dim)
      write(nfout, '("|||| nis_fft_Z_y :")')
      write(nfout,'(("||||",10(i8)))') (nis_fft_Z_y(i),i=1,fft_Z_y_dim)
      write(nfout, '("|||| nie_fft_Z_y :")')
      write(nfout,'(("||||",10(i8)))') (nie_fft_Z_y(i),i=1,fft_Z_y_dim)
      write(nfout, '("|||| fft_Y_x_dim=",i4,",  fft_Y_z_dim=",i4)') fft_Y_x_dim, fft_Y_z_dim
      write(nfout, '("|||| nis_fft_Y_x :")')
      write(nfout,'(("||||",10(i8)))') (nis_fft_Y_x(i),i=1,fft_Y_x_dim)
      write(nfout, '("|||| nie_fft_Y_x :")')
      write(nfout,'(("||||",10(i8)))') (nie_fft_Y_x(i),i=1,fft_Y_x_dim)
      write(nfout, '("|||| nis_fft_Y_z :")')
      write(nfout,'(("||||",10(i8)))') (nis_fft_Y_z(i),i=1,fft_Y_z_dim)
      write(nfout, '("|||| nie_fft_Y_z :")')
      write(nfout,'(("||||",10(i8)))') (nie_fft_Y_z(i),i=1,fft_Y_z_dim)
    endif

#ifdef FFT_ALLTOALL
    if (ipri > 1) then
      write(nfout,'("nrank_e=",i4,", myrank_e=",i4)') nrank_e, myrank_e
      write(nfout,'("nrank_g=",i4,", myrank_g=",i4)') nrank_g, myrank_g
! === DEBUG by tkato 2012/06/06 ================================================
!     write(nfout,'("mpi_ke_world=",i4)') mpi_ke_world
!     write(nfout,'("mpi_kg_world=",i4)') mpi_kg_world
      write(nfout,'("mpi_ke_world=",i12)') mpi_ke_world
      write(nfout,'("mpi_kg_world=",i12)') mpi_kg_world
! ==============================================================================
      write(nfout,'("dim1=",i4,", dim2=",i4)') dim1, dim2
    endif
!   color = mod(myrank_g,dim2)
!   color = myrank_g/dim2
!   color = mod(myrank_g,dim1)
    color = myrank_g/dim1
!    if (myrank_g > (dim1*dim2-1)) then
!       color = MPI_UNDEFINED
!    end if
    if( nrank_g>dim1*dim2 ) then
      write(strng,'(i10)') nrank_g
      write(strngmax,'(i10)') dim1*dim2
      call phase_error_with_msg(nfout,'ng = '//trim(strng)//', max ng ='//trim(strngmax), &
      __LINE__,__FILE__)
    endif
    key = myrank_g
!    if(firstcall) call mpi_comm_split(mpi_ke_world,color,key,mpi_fft_xy_world,ierr)
    if(.not.firstcall) call mpi_comm_free(mpi_fft_xy_world, ierr)
    call mpi_comm_split(mpi_ke_world,color,key,mpi_fft_xy_world,ierr)
    call mpi_comm_size(mpi_fft_xy_world, nrank_fft_xy, ierr)
    call mpi_comm_rank(mpi_fft_xy_world, myrank_fft_xy, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_mpi_fft_box :  mpi_comm_split error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 100002, ierr)
     endif

     call mpi_barrier(mpi_fft_xy_world, ierr)

    if (ipri > 1) then
      write(nfout,'(" mpi_fft_xy_world=",i12)') mpi_fft_xy_world
      write(nfout,'("  color=",i4,", key=",i4)') color, key
      write(nfout,'("  nrank_fft_xy=",i4,", myrank_fft_xy=",i4)') nrank_fft_xy, myrank_fft_xy
      call flush(nfout)
    end if

!   color = myrank_g/dim2
!   color = mod(myrank_g,dim2)
!   color = myrank_g/dim1
    color = mod(myrank_g,dim1)
    !if (myrank_g > (dim1*dim2-1)) then
    !   color = MPI_UNDEFINED
    !end if
    if( nrank_g>dim1*dim2 ) then
      write(strng,*) nrank_g
      write(strngmax,*) dim1*dim2
      call phase_error_with_msg(nfout,'ng = '//trim(strng)//', max ng ='//trim(strngmax), &
      __LINE__,__FILE__)
    endif
    key = myrank_g
!    if(firstcall) call mpi_comm_split(mpi_ke_world,color,key,mpi_fft_yz_world,ierr)
    if(.not.firstcall) call mpi_comm_free(mpi_fft_yz_world, ierr)
    call mpi_comm_split(mpi_ke_world,color,key,mpi_fft_yz_world,ierr)
    call mpi_comm_size(mpi_fft_yz_world, nrank_fft_yz, ierr)
    call mpi_comm_rank(mpi_fft_yz_world, myrank_fft_yz, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_mpi_fft_box :  mpi_comm_split error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 100002, ierr)
     endif

    if (ipri > 1) then
      write(nfout,'(" mpi_fft_yz_world=",i12)') mpi_fft_yz_world
      write(nfout,'("  color=",i4,", key=",i4)') color, key
      write(nfout,'("  nrank_fft_yz=",i4,", myrank_fft_yz=",i4)') nrank_fft_yz, myrank_fft_yz
      call flush(nfout)
    endif
#endif
    firstcall=.false.

                                                  __TIMER_SUB_STOP(1232)
  end subroutine m_Parallel_mpi_fft_box_xyz
!===============================================================================
!===============================================================================
  subroutine m_Parallel_mpi_fft_box(nfout,ipri,printable,fft_box_size_WF,kimg,divide)
    implicit none

    integer, intent(in) :: nfout,ipri, kimg, divide
    logical, intent(in) :: printable
    integer, intent(in) :: fft_box_size_WF(3,0:1)
    integer, allocatable, dimension(:,:) :: wk_z, wk_x, wk_y, wk1, wk2
    integer, allocatable, dimension(:,:,:) :: wk_fft_box
    integer :: dim1, dim2, igf, igf1, igf2, igf3, kdim1, kdim2
    integer :: i, j, id, id1, id2, id3, ierr, ichkalloc
    integer :: ii, jj, kk, ix, iy, iz, mindim1, mindim2
    integer :: nfft, n, nmrank, myrank
    integer :: key, color
    logical, save :: firstcall=.true.
    character(len=12) :: strng,strngmax
                                                  __TIMER_SUB_START(1232)

!!  nmrank = nrank_e
!!  myrank = myrank_e
    nmrank = nrank_g
    myrank = myrank_g

    if ((divide > 0) .and. (nmrank > 3)) then
       dim1 = int(sqrt(real(nmrank)))
       dim2 = dim1
    else
       mindim1 = min(fft_box_size_WF(2,0),fft_box_size_WF(3,0))
       if (kimg == 1) then
          mindim2 = min(fft_box_size_WF(1,0)/2,fft_box_size_WF(3,0))
!!$       if (nmrank < 8) then
!!$          sdim = 1
       else
          mindim2 = min(fft_box_size_WF(1,0),fft_box_size_WF(3,0))
!!$          sdim = 2
       end if

       call decide_div_fft(nmrank, mindim1, mindim2, dim2, dim1)
    end if

    nfft = fft_box_size_WF(1,0)*fft_box_size_WF(2,0)*fft_box_size_WF(3,0)

    allocate(wk_z( 4,nmrank ), stat=ichkalloc)
    allocate(wk_x( 4,nmrank ), stat=ichkalloc)
    allocate(wk_y( 4,nmrank ), stat=ichkalloc)
    if(ichkalloc /= 0) then
       write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 1'
       call mpi_abort(mpi_comm_world, 3 , ichkalloc)
    endif
    wk_z(:,:) = 0
    wk_x(:,:) = 0
    wk_y(:,:) = 0

    allocate(wk1(2,dim1))
    allocate(wk2(2,dim2))
    wk1 = 0
    wk2 = 0

    call index_div_fft(kimg, wk2, dim2, fft_box_size_WF(2,0), 0)
    ii = 0
    do i = 1, dim2
       if (wk2(1,i) == 0) cycle
       ii = ii + 1
    end do
    kdim2 = ii
    if (kdim2 /= dim2) then
       deallocate(wk2)
       dim2 = kdim2
       allocate(wk2(2,dim2))
    end if
    wk2 = 0

    call index_div_fft(kimg, wk1, dim1, fft_box_size_WF(3,0), 0)
    ii = 0
    do i = 1, dim1
       if (wk1(1,i) == 0) cycle
       ii = ii + 1
    end do
    kdim1 = ii

    if (kdim1 /= dim1) then
       deallocate(wk1)
       dim1 = kdim1
       allocate(wk1(2,dim1))
    end if
    wk1 = 0

    fft_X_z_dim = dim1
    fft_X_y_dim = dim2
    fft_Z_x_dim = dim1
    fft_Z_y_dim = dim2
    fft_Y_x_dim = dim1
    fft_Y_z_dim = dim2

    allocate(nis_fft_X_z(fft_X_z_dim), stat=ichkalloc)
    allocate(nie_fft_X_z(fft_X_z_dim), stat=ichkalloc)
    allocate(nis_fft_X_y(fft_X_y_dim), stat=ichkalloc)
    allocate(nie_fft_X_y(fft_X_y_dim), stat=ichkalloc)
    allocate(nis_fft_Z_x(fft_Z_x_dim), stat=ichkalloc)
    allocate(nie_fft_Z_x(fft_Z_x_dim), stat=ichkalloc)
    allocate(nis_fft_Z_y(fft_Z_y_dim), stat=ichkalloc)
    allocate(nie_fft_Z_y(fft_Z_y_dim), stat=ichkalloc)
    allocate(nis_fft_Y_x(fft_Y_x_dim), stat=ichkalloc)
    allocate(nie_fft_Y_x(fft_Y_x_dim), stat=ichkalloc)
    allocate(nis_fft_Y_z(fft_Y_z_dim), stat=ichkalloc)
    allocate(nie_fft_Y_z(fft_Y_z_dim), stat=ichkalloc)

    call index_div_fft(kimg, wk2, dim2, fft_box_size_WF(2,0), 0)
    call index_div_fft(kimg, wk1, dim1, fft_box_size_WF(3,0), 0)
    nis_fft_X_y(:) = wk2(1,:)
    nie_fft_X_y(:) = wk2(2,:)
    nis_fft_X_z(:) = wk1(1,:)
    nie_fft_X_z(:) = wk1(2,:)
    ii = 1
    jj = 1
    do i = 1, dim1*dim2
       wk_x(1,i) = wk2(1,ii)
       wk_x(2,i) = wk2(2,ii)
       wk_x(3,i) = wk1(1,jj)
       wk_x(4,i) = wk1(2,jj)
       ii = ii + 1
       if (ii > dim2) then
          ii = 1
          jj = jj + 1
       end if
    end do
    call index_div_fft(kimg, wk1, dim1, fft_box_size_WF(1,0), 1)
!   call index_div(wk2, dim2, fft_box_size_WF(2,0))
    nis_fft_Z_x(:) = wk1(1,:)
    nie_fft_Z_x(:) = wk1(2,:)
    nis_fft_Z_y(:) = wk2(1,:)
    nie_fft_Z_y(:) = wk2(2,:)
    ii = 1
    jj = 1
    do i = 1, dim1*dim2
       wk_z(1,i) = wk2(1,ii)
       wk_z(2,i) = wk2(2,ii)
       wk_z(3,i) = wk1(1,jj)
       wk_z(4,i) = wk1(2,jj)
       ii = ii + 1
       if (ii > dim2) then
          ii = 1
          jj = jj + 1
       end if
    end do
!   call index_div(wk1, dim1, fft_box_size_WF(1,0))
    call index_div_fft(kimg, wk2, dim2, fft_box_size_WF(3,0), 0)
    nis_fft_Y_x(:) = wk1(1,:)
    nie_fft_Y_x(:) = wk1(2,:)
    nis_fft_Y_z(:) = wk2(1,:)
    nie_fft_Y_z(:) = wk2(2,:)
    ii = 1
    jj = 1
    do i = 1, dim1*dim2
       wk_y(1,i) = wk2(1,ii)
       wk_y(2,i) = wk2(2,ii)
       wk_y(3,i) = wk1(1,jj)
       wk_y(4,i) = wk1(2,jj)
       ii = ii + 1
       if (ii > dim2) then
          ii = 1
          jj = jj + 1
       end if
    end do

!!!!
!!  n = 1
!!!!
    n = 0
!!  call div_index(fft_box_size_WF(2,n), fft_box_size_WF(3,n), dim1, dim2, wk_x, nmrank)
    if (wk_x(1,myrank+1) /= 0) then
      xyz_fft_x(1,1) = 1
      xyz_fft_x(2,1) = fft_box_size_WF(1,n)
      xyz_fft_x(1,2) = wk_x(1,myrank+1)
      xyz_fft_x(2,2) = wk_x(2,myrank+1)
      xyz_fft_x(1,3) = wk_x(3,myrank+1)
      xyz_fft_x(2,3) = wk_x(4,myrank+1)
    endif

!!  call div_index(fft_box_size_WF(2,n), fft_box_size_WF(1,n), dim1, dim2, wk_z, nmrank)
    if (wk_z(1,myrank+1) /= 0) then
      xyz_fft_z(1,1) = wk_z(3,myrank+1)
      xyz_fft_z(2,1) = wk_z(4,myrank+1)
      xyz_fft_z(1,2) = wk_z(1,myrank+1)
      xyz_fft_z(2,2) = wk_z(2,myrank+1)
      xyz_fft_z(1,3) = 1
      xyz_fft_z(2,3) = fft_box_size_WF(3,n)
    endif

!!  call div_index(fft_box_size_WF(3,n), fft_box_size_WF(1,n), dim1, dim2, wk_y, nmrank)
    if (wk_y(1,myrank+1) /= 0) then
      xyz_fft_y(1,1) = wk_y(3,myrank+1)
      xyz_fft_y(2,1) = wk_y(4,myrank+1)
      xyz_fft_y(1,2) = 1
      xyz_fft_y(2,2) = fft_box_size_WF(2,n)
      xyz_fft_y(1,3) = wk_y(1,myrank+1)
      xyz_fft_y(2,3) = wk_y(2,myrank+1)
    endif

    if(allocated(map_fft_x)) deallocate(map_fft_x)
    if(allocated(map_fft_z)) deallocate(map_fft_z)
    if(allocated(map_fft_y)) deallocate(map_fft_y)
    allocate(map_fft_x( nfft ), stat=ichkalloc)
    allocate(map_fft_z( nfft ), stat=ichkalloc)
    allocate(map_fft_y( nfft ), stat=ichkalloc)
    if(ichkalloc /= 0) then
       write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 2'
       call mpi_abort(mpi_comm_world, 4 , ichkalloc)
    endif
    map_fft_x(:) = 0
    map_fft_z(:) = 0
    map_fft_y(:) = 0

    np_fft_x = 0
    np_fft_z = 0
    np_fft_y = 0

    id1 = fft_box_size_WF(1,0)
    id2 = fft_box_size_WF(2,0)
    id3 = fft_box_size_WF(3,0)

    allocate( wk_fft_box(fft_box_size_WF(1,0), fft_box_size_WF(2,0), fft_box_size_WF(3,0)))
    do igf3 = 1, fft_box_size_WF(3,n)
       do igf2 = 1, fft_box_size_WF(2,n)
          do igf1 = 1, fft_box_size_WF(1,n)
             igf = igf1 + (igf2-1)*id1 + (igf3-1)*id1*id2

             wk_fft_box(igf1,igf2,igf3) = igf

!            do j = 1, nmrank
!              if(((igf2.ge.wk_x(1,j)).and.(igf2.le.wk_x(2,j))) .and. &
!             &   ((igf3.ge.wk_x(3,j)).and.(igf3.le.wk_x(4,j)))) then
!                map_fft_x(igf) = j
!                if (myrank .eq. (j-1)) np_fft_x = np_fft_x + 1
!              endif
!              if(((igf1.ge.wk_z(3,j)).and.(igf1.le.wk_z(4,j))) .and. &
!             &   ((igf2.ge.wk_z(1,j)).and.(igf2.le.wk_z(2,j)))) then
!                map_fft_z(igf) = j
!                if (myrank .eq. (j-1)) np_fft_z = np_fft_z + 1
!              endif
!              if(((igf1.ge.wk_y(3,j)).and.(igf1.le.wk_y(4,j))) .and. &
!             &   ((igf3.ge.wk_y(1,j)).and.(igf3.le.wk_y(2,j)))) then
!                map_fft_y(igf) = j
!                if (myrank .eq. (j-1)) np_fft_y = np_fft_y + 1
!              endif
!            enddo

          enddo
       enddo
    enddo

!x  do ii = 1, nmrank
    do ii = 1, dim1*dim2
       if (wk_x(1,ii) == 0) cycle
       do iz = wk_x(3,ii) , wk_x(4,ii)
          do iy = wk_x(1,ii) , wk_x(2,ii)
             do ix = 1, id1
                map_fft_x(id1*id2*(iz-1)+id1*(iy-1)+ix) = ii
             end do
          end do
       end do
    end do

!x  do ii = 1, nmrank
    do ii = 1, dim1*dim2
       if (wk_z(1,ii) == 0) cycle
       do iz = 1, id3
          do iy = wk_z(1,ii) , wk_z(2,ii)
             do ix = wk_z(3,ii), wk_z(4,ii)
                map_fft_z(id1*id2*(iz-1)+id1*(iy-1)+ix) = ii
             end do
          end do
       end do
    end do

!x  do ii = 1, nmrank
    do ii = 1, dim1*dim2
       if (wk_y(1,ii) == 0) cycle
       do iz = wk_y(1,ii) , wk_y(2,ii)
          do iy = 1, id2
             do ix = wk_y(3,ii), wk_y(4,ii)
                map_fft_y(id1*id2*(iz-1)+id1*(iy-1)+ix) = ii
             end do
          end do
       end do
    end do

    if ((dim1*dim2) > myrank) then
       if (wk_x(1,myrank+1) /= 0) then
          np_fft_x = id1*(wk_x(2,myrank+1)-wk_x(1,myrank+1)+1)*(wk_x(4,myrank+1)-wk_x(3,myrank+1)+1)
          if(allocated(mp_fft_x)) deallocate(mp_fft_x)
          allocate(mp_fft_x(np_fft_x), stat=ichkalloc)
       end if
       if (wk_z(1,myrank+1) /= 0) then
          np_fft_z = (wk_z(4,myrank+1)-wk_z(3,myrank+1)+1)*(wk_z(2,myrank+1)-wk_z(1,myrank+1)+1)*id3
          if(allocated(mp_fft_z)) deallocate(mp_fft_z)
          allocate(mp_fft_z(np_fft_z), stat=ichkalloc)
       end if
       if (wk_y(1,myrank+1) /= 0) then
          np_fft_y = (wk_y(4,myrank+1)-wk_y(3,myrank+1)+1)*id2*(wk_y(2,myrank+1)-wk_y(1,myrank+1)+1)
          if(allocated(mp_fft_y)) deallocate(mp_fft_y)
          allocate(mp_fft_y(np_fft_y), stat=ichkalloc)
       end if
       if(ichkalloc /= 0) then
          write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 3'
          call mpi_abort(mpi_comm_world, 5 , ichkalloc)
       endif
    else
       np_fft_x = 0
       np_fft_z = 0
       np_fft_y = 0
    end if

    deallocate(wk_x)
    deallocate(wk_z)
    deallocate(wk_y)

    if (kimg == 1) then
       if (np_fft_x /= 0) then
          id = 0
          do kk = xyz_fft_x(1,3), xyz_fft_x(2,3)
            do jj = xyz_fft_x(1,2), xyz_fft_x(2,2)
              do ii = xyz_fft_x(1,1), xyz_fft_x(2,1) , 2
                id = id + 1
                mp_fft_x(id) = wk_fft_box(ii,jj,kk)
                id = id + 1
                mp_fft_x(id) = wk_fft_box(ii+1,jj,kk)
              enddo
            enddo
          enddo
       end if
       if (np_fft_y /= 0) then
          id = 0
            do jj = xyz_fft_y(1,2), xyz_fft_y(2,2)
          do kk = xyz_fft_y(1,3), xyz_fft_y(2,3)
              do ii = xyz_fft_y(1,1), xyz_fft_y(2,1) , 2
                id = id + 1
                mp_fft_y(id) = wk_fft_box(ii,jj,kk)
                id = id + 1
                mp_fft_y(id) = wk_fft_box(ii+1,jj,kk)
              enddo
            enddo
          enddo
       end if
       if (np_fft_z /= 0) then
          id = 0
          do kk = xyz_fft_z(1,3), xyz_fft_z(2,3)
            do jj = xyz_fft_z(1,2), xyz_fft_z(2,2)
              do ii = xyz_fft_z(1,1), xyz_fft_z(2,1) , 2
                id = id + 1
                mp_fft_z(id) = wk_fft_box(ii,jj,kk)
                id = id + 1
                mp_fft_z(id) = wk_fft_box(ii+1,jj,kk)
              enddo
            enddo
          enddo
       end if
    else
       if (np_fft_x /= 0) then
          id = 0
          do kk = xyz_fft_x(1,3), xyz_fft_x(2,3)
            do jj = xyz_fft_x(1,2), xyz_fft_x(2,2)
              do ii = xyz_fft_x(1,1), xyz_fft_x(2,1)
                id = id + 1
                mp_fft_x(id) = wk_fft_box(ii,jj,kk)
              enddo
            enddo
          enddo
       end if
       if (np_fft_y /= 0) then
          id = 0
            do jj = xyz_fft_y(1,2), xyz_fft_y(2,2)
          do kk = xyz_fft_y(1,3), xyz_fft_y(2,3)
              do ii = xyz_fft_y(1,1), xyz_fft_y(2,1)
                id = id + 1
                mp_fft_y(id) = wk_fft_box(ii,jj,kk)
              enddo
            enddo
          enddo
       end if
       if (np_fft_z /= 0) then
          id = 0
          do kk = xyz_fft_z(1,3), xyz_fft_z(2,3)
            do jj = xyz_fft_z(1,2), xyz_fft_z(2,2)
              do ii = xyz_fft_z(1,1), xyz_fft_z(2,1)
                id = id + 1
                mp_fft_z(id) = wk_fft_box(ii,jj,kk)
              enddo
            enddo
          enddo
       end if
    endif

    deallocate(wk_fft_box)

    if(allocated(nel_fft_x)) deallocate(nel_fft_x)
    if(allocated(nel_fft_z)) deallocate(nel_fft_z)
    if(allocated(nel_fft_y)) deallocate(nel_fft_y)
    allocate(nel_fft_x( 0:nmrank-1 ), stat=ichkalloc)
    allocate(nel_fft_z( 0:nmrank-1 ), stat=ichkalloc)
    allocate(nel_fft_y( 0:nmrank-1 ), stat=ichkalloc)
    if(ichkalloc /= 0) then
       write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 4'
       call mpi_abort(mpi_comm_world, 6 , ichkalloc)
    endif
    nel_fft_z(:) = 0
    nel_fft_x(:) = 0
    nel_fft_y(:) = 0

    do i = 1, nfft
      if (map_fft_x(i) .gt. 0) then
         nel_fft_x( map_fft_x(i)-1 ) = nel_fft_x( map_fft_x(i)-1 ) + 1
      endif
      if (map_fft_z(i) .gt. 0) then
         nel_fft_z( map_fft_z(i)-1 ) = nel_fft_z( map_fft_z(i)-1 ) + 1
      endif
      if (map_fft_y(i) .gt. 0) then
         nel_fft_y( map_fft_y(i)-1 ) = nel_fft_y( map_fft_y(i)-1 ) + 1
      endif
    enddo

    if (ipri > 1) then
      write(nfout,'("|||| m_Parallel_mpi_fft_box")')
      write(nfout,'("|||| dim1 = ",i5," , dim2 = ",i5)')   dim1,dim2
      write(nfout,'("|||| fft-box-size (",i4," ,",i4," ,",i4,")")') &
     &          ((fft_box_size_WF(i,j),i=1,3),j=0,0)
      write(nfout,'("|||| fft-box-div[X-axis] (",i4," :",i4," ,",i4," :",i4," ,",i4," :",i4,")")') &
     &          ((xyz_fft_x(i,j),i=1,2),j=1,3)
      write(nfout,'("|||| fft-box-div[Z-axis] (",i4," :",i4," ,",i4," :",i4," ,",i4," :",i4,")")') &
     &          ((xyz_fft_z(i,j),i=1,2),j=1,3)
      write(nfout,'("|||| fft-box-div[Y-axis] (",i4," :",i4," ,",i4," :",i4," ,",i4," :",i4,")")') &
     &          ((xyz_fft_y(i,j),i=1,2),j=1,3)
!     write(nfout,'("|||| ir_x2z=",i4," , ir_z2y=",i4," , ir_y2z=",i4," , ir_z2x=",i4)') &
!    &                    ir_x2z,ir_z2y,ir_y2z,ir_z2x
!     write(nfout,'("|||| is_x2z=",i4," , is_z2y=",i4," , is_y2z=",i4," , is_z2x=",i4)') &
!    &                    is_x2z,is_z2y,is_y2z,is_z2x
      write(nfout,'("|||| np_fft_x=",i8)')          np_fft_x
      write(nfout,'("|||| nel_fft_x:",/,"||||",10(i7))') (nel_fft_x(i),i=0,nmrank-1)
      write(nfout,'("|||| np_fft_y=",i8)')          np_fft_y
      write(nfout,'("|||| nel_fft_y:",/,"||||",10(i7))') (nel_fft_y(i),i=0,nmrank-1)
      write(nfout,'("|||| np_fft_z=",i8)')          np_fft_z
      write(nfout,'("|||| nel_fft_z:",/,"||||",10(i7))') (nel_fft_z(i),i=0,nmrank-1)
      write(nfout, '("|||| fft_X_y_dim=",i4,",  fft_X_z_dim=",i4)') fft_X_y_dim, fft_X_z_dim
      write(nfout, '("|||| nis_fft_X_y :")')
      write(nfout,'(("||||",10(i8)))') (nis_fft_X_y(i),i=1,fft_X_y_dim)
      write(nfout, '("|||| nie_fft_X_y :")')
      write(nfout,'(("||||",10(i8)))') (nie_fft_X_y(i),i=1,fft_X_y_dim)
      write(nfout, '("|||| nis_fft_X_z :")')
      write(nfout,'(("||||",10(i8)))') (nis_fft_X_z(i),i=1,fft_X_z_dim)
      write(nfout, '("|||| nie_fft_X_z :")')
      write(nfout,'(("||||",10(i8)))') (nie_fft_X_z(i),i=1,fft_X_z_dim)
      write(nfout, '("|||| fft_Z_x_dim=",i4,",  fft_Z_y_dim=",i4)') fft_Z_x_dim, fft_Z_y_dim
      write(nfout, '("|||| nis_fft_Z_x :")')
      write(nfout,'(("||||",10(i8)))') (nis_fft_Z_x(i),i=1,fft_Z_x_dim)
      write(nfout, '("|||| nie_fft_Z_x :")')
      write(nfout,'(("||||",10(i8)))') (nie_fft_Z_x(i),i=1,fft_Z_x_dim)
      write(nfout, '("|||| nis_fft_Z_y :")')
      write(nfout,'(("||||",10(i8)))') (nis_fft_Z_y(i),i=1,fft_Z_y_dim)
      write(nfout, '("|||| nie_fft_Z_y :")')
      write(nfout,'(("||||",10(i8)))') (nie_fft_Z_y(i),i=1,fft_Z_y_dim)
      write(nfout, '("|||| fft_Y_x_dim=",i4,",  fft_Y_z_dim=",i4)') fft_Y_x_dim, fft_Y_z_dim
      write(nfout, '("|||| nis_fft_Y_x :")')
      write(nfout,'(("||||",10(i8)))') (nis_fft_Y_x(i),i=1,fft_Y_x_dim)
      write(nfout, '("|||| nie_fft_Y_x :")')
      write(nfout,'(("||||",10(i8)))') (nie_fft_Y_x(i),i=1,fft_Y_x_dim)
      write(nfout, '("|||| nis_fft_Y_z :")')
      write(nfout,'(("||||",10(i8)))') (nis_fft_Y_z(i),i=1,fft_Y_z_dim)
      write(nfout, '("|||| nie_fft_Y_z :")')
      write(nfout,'(("||||",10(i8)))') (nie_fft_Y_z(i),i=1,fft_Y_z_dim)
    endif

#ifdef FFT_ALLTOALL
    if (ipri > 1) then
      write(nfout,'("nrank_e=",i4,", myrank_e=",i4)') nrank_e, myrank_e
      write(nfout,'("nrank_g=",i4,", myrank_g=",i4)') nrank_g, myrank_g
! === DEBUG by tkato 2012/06/06 ================================================
!     write(nfout,'("mpi_ke_world=",i4)') mpi_ke_world
!     write(nfout,'("mpi_kg_world=",i4)') mpi_kg_world
      write(nfout,'("mpi_ke_world=",i12)') mpi_ke_world
      write(nfout,'("mpi_kg_world=",i12)') mpi_kg_world
! ==============================================================================
      write(nfout,'("dim1=",i4,", dim2=",i4)') dim1, dim2
    endif

    color = mod(myrank_g,dim2)
    !if (myrank_g > (dim1*dim2-1)) then
    !   color = MPI_UNDEFINED
    !end if
    if( nrank_g>dim1*dim2 ) then
      write(strng,'(i10)') nrank_g
      write(strngmax,'(i10)') dim1*dim2
      call phase_error_with_msg(nfout,'ng = '//trim(strng)//', max ng ='//trim(strngmax), &
      __LINE__,__FILE__)
    endif
    key = myrank_g
!    if(firstcall) call mpi_comm_split(mpi_ke_world,color,key,mpi_fft_xz_world,ierr)
    if(.not.firstcall) call mpi_comm_free(mpi_fft_xz_world, ierr)
    call mpi_comm_split(mpi_ke_world,color,key,mpi_fft_xz_world,ierr)
    call mpi_comm_size(mpi_fft_xz_world, nrank_fft_xz, ierr)
    call mpi_comm_rank(mpi_fft_xz_world, myrank_fft_xz, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_mpi_fft_box :  mpi_comm_split error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 100002, ierr)
     endif

     call mpi_barrier(mpi_fft_xz_world, ierr)

    if (ipri > 1) then
! === DEBUG by tkato 2012/06/06 ================================================
!     write(nfout,'(" mpi_fft_xz_world=",i8)') mpi_fft_xz_world
      write(nfout,'(" mpi_fft_xz_world=",i12)') mpi_fft_xz_world
! ==============================================================================
      write(nfout,'("  color=",i4,", key=",i4)') color, key
      write(nfout,'("  nrank_fft_xz=",i4,", myrank_fft_xz=",i4)') nrank_fft_xz, myrank_fft_xz
      call flush(nfout)
    end if

    color = myrank_g/dim2
    !if (myrank_g > (dim1*dim2-1)) then
    !   color = MPI_UNDEFINED
    !end if
    if( nrank_g>dim1*dim2 ) then
      write(strng,*) nrank_g
      write(strngmax,*) dim1*dim2
      call phase_error_with_msg(nfout,'ng = '//trim(strng)//', max ng ='//trim(strngmax), &
      __LINE__,__FILE__)
    endif
    key = myrank_g
!    if(firstcall) call mpi_comm_split(mpi_ke_world,color,key,mpi_fft_zy_world,ierr)
    if(.not.firstcall) call mpi_comm_free(mpi_fft_zy_world, ierr)
    call mpi_comm_split(mpi_ke_world,color,key,mpi_fft_zy_world,ierr)
    call mpi_comm_size(mpi_fft_zy_world, nrank_fft_zy, ierr)
    call mpi_comm_rank(mpi_fft_zy_world, myrank_fft_zy, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_mpi_fft_box :  mpi_comm_split error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 100002, ierr)
     endif

    if (ipri > 1) then
! === DEBUG by tkato 2012/06/06 ================================================
!     write(nfout,'(" mpi_fft_zy_world=",i8)') mpi_fft_zy_world
      write(nfout,'(" mpi_fft_zy_world=",i12)') mpi_fft_zy_world
! ==============================================================================
      write(nfout,'("  color=",i4,", key=",i4)') color, key
      write(nfout,'("  nrank_fft_zy=",i4,", myrank_fft_zy=",i4)') nrank_fft_zy, myrank_fft_zy
      call flush(nfout)
    endif

#endif
    firstcall=.false.

                                                  __TIMER_SUB_STOP(1232)
  end subroutine m_Parallel_mpi_fft_box
!===============================================================================
!===============================================================================
  subroutine m_Parallel_mpi_fft_box_cd_3div(nfout,ipri,printable,fft_box_size_CD,kimg, &
 &                                  fftbox_div_1,fftbox_div_2,fftbox_div_3)
    implicit none

    integer, intent(in) :: nfout, ipri, kimg, fftbox_div_1 ,fftbox_div_2, fftbox_div_3
    logical, intent(in) :: printable
    integer, intent(in) :: fft_box_size_CD(3,0:1)
    integer, allocatable, dimension(:,:) :: wk_x
    integer, allocatable, dimension(:,:,:) :: wk_fft_box
    integer :: div1, div2, div3, igf, igf1, igf2, igf3
    integer :: i, j, id, id1, id2, id3, ierr, ichkalloc
    integer :: ii, jj, kk, ix, iy, iz
    integer :: nfft, n, nmrank, myrank, nfftm2
    integer :: rank1, rank2, rank3, nfft1, nfft2, nfft3, rank
    integer :: kx1p, kx2p, kx3p, n1s, n2s, n3s, n1e, n2e, n3e, n1p, n2p, n3p
    integer :: key, color
    logical, save :: firstcall=.true.
    character(len=12) :: strng,strngmax
#ifdef __TIMER_SUB__
  call timer_sta(1233)
#endif


    nmrank = nrank_g
    myrank = myrank_g

    if ((fftbox_div_1 < 0) .or. (fftbox_div_2 < 0) .or. (fftbox_div_3 < 0)) then
       write(nfout,'("Either of the division numbers of FFT-BOX is negative.")')
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 1234567 , ierr)
    end if
    if ((fftbox_div_1*fftbox_div_2*fftbox_div_3) > nmrank ) then
       write(nfout,'("The three-dimensional division number of FFT-BOX is over the number of processes.")')
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 1234568 , ierr)
    end if

    if ((fftbox_div_1*fftbox_div_2*fftbox_div_3) == 0 ) then
! === DEBUG by tkato 2013/09/18 ================================================
!     div1 = int(cbrt(real(nmrank)))
#ifdef __FUJITSU
      div1 = int(cbrt(real(nmrank)))
#else
      div1 = int((real(nmrank))**(1.0d0/3.0d0))
#endif
! ==============================================================================
      div2 = div1
      div3 = div1
    else
      div1 = fftbox_div_1
      div2 = fftbox_div_2
      div3 = fftbox_div_3
    end if

    nfft   = fft_box_size_CD(1,0)*fft_box_size_CD(2,0)*fft_box_size_CD(3,0)
    nfftm2 = fft_box_size_CD(1,1)*fft_box_size_CD(2,1)*fft_box_size_CD(3,1)
    nfft1 = fft_box_size_CD(1,1)
    nfft2 = fft_box_size_CD(2,1)
    nfft3 = fft_box_size_CD(3,1)

    if(mod(nfft1,div1) > 0) then
      kx1p = nfft1/div1+1
    else
      kx1p = nfft1/div1
    end if
    if(mod(nfft2,div2) > 0) then
      kx2p = nfft2/div2+1
    else
      kx2p = nfft2/div2
    end if
    if(mod(nfft3,div3) > 0) then
      kx3p = nfft3/div3+1
    else
      kx3p = nfft3/div3
    end if

    allocate(wk_x( 6,nmrank ), stat=ichkalloc)
    if(ichkalloc /= 0) then
       write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 1'
       call mpi_abort(mpi_comm_world, 3 , ichkalloc)
    endif
    wk_x(:,:) = 0

    fftcd_X_x_dim = div1
    fftcd_X_y_dim = div2
    fftcd_X_z_dim = div3
    fftcd_X_x_nel = kx1p
    fftcd_X_y_nel = kx2p
    fftcd_X_z_nel = kx3p

    allocate(nis_fftcd_X_z(fftcd_X_z_dim), stat=ichkalloc)
    allocate(nie_fftcd_X_z(fftcd_X_z_dim), stat=ichkalloc)
    allocate(nis_fftcd_X_y(fftcd_X_y_dim), stat=ichkalloc)
    allocate(nie_fftcd_X_y(fftcd_X_y_dim), stat=ichkalloc)
    allocate(nis_fftcd_X_x(fftcd_X_x_dim), stat=ichkalloc)
    allocate(nie_fftcd_X_x(fftcd_X_x_dim), stat=ichkalloc)

    do rank = 0, div1*div2*div3-1
      rank1 = mod(rank,div1)
      rank2 = mod(rank/div1,div2)
      rank3 = rank/(div1*div2)

      n1s = kx1p*rank1+1
      n1e = min(nfft1,kx1p*(rank1+1))
      n1p = max(0,n1e-n1s+1)
      n2s = kx2p*rank2+1
      n2e = min(nfft2,kx2p*(rank2+1))
      n2p = max(0,n2e-n2s+1)
      n3s = kx3p*rank3+1
      n3e = min(nfft3,kx3p*(rank3+1))
      n3p = max(0,n3e-n3s+1)

      wk_x(1,rank+1) = n1s
      wk_x(2,rank+1) = n1e
      wk_x(3,rank+1) = n2s
      wk_x(4,rank+1) = n2e
      wk_x(5,rank+1) = n3s
      wk_x(6,rank+1) = n3e

      nis_fftcd_X_x(rank1+1) = n1s
      nie_fftcd_X_x(rank1+1) = n1e
      nis_fftcd_X_y(rank2+1) = n2s
      nie_fftcd_X_y(rank2+1) = n2e
      nis_fftcd_X_z(rank3+1) = n3s
      nie_fftcd_X_z(rank3+1) = n3e
    end do

!!!!
!!  n = 1
!!!!
    n = 0
    if (wk_x(1,myrank+1) /= 0) then
      xyz_fftcd_x(1,1) = wk_x(1,myrank+1)
      xyz_fftcd_x(2,1) = wk_x(2,myrank+1)
      xyz_fftcd_x(1,2) = wk_x(3,myrank+1)
      xyz_fftcd_x(2,2) = wk_x(4,myrank+1)
      xyz_fftcd_x(1,3) = wk_x(5,myrank+1)
      xyz_fftcd_x(2,3) = wk_x(6,myrank+1)
    endif

    if(allocated(map_fftcd_x)) deallocate(map_fftcd_x)
    allocate(map_fftcd_x( nfft ), stat=ichkalloc)
    if(ichkalloc /= 0) then
       write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 2'
       call mpi_abort(mpi_comm_world, 4 , ichkalloc)
    endif
    map_fftcd_x(:) = 0

    np_fftcd_x = 0

    id1 = fft_box_size_CD(1,0)
    id2 = fft_box_size_CD(2,0)
    id3 = fft_box_size_CD(3,0)

    allocate( wk_fft_box(id1, id2, id3) )
    do igf3 = 1, fft_box_size_CD(3,n)
       do igf2 = 1, fft_box_size_CD(2,n)
          do igf1 = 1, fft_box_size_CD(1,n)
             igf = igf1 + (igf2-1)*id1 + (igf3-1)*id1*id2

             wk_fft_box(igf1,igf2,igf3) = igf

          enddo
       enddo
    enddo

    do rank = 1, div1*div2*div3
       if (wk_x(1,rank) == 0) cycle
       do iz = wk_x(5,rank) , wk_x(6,rank)
          do iy = wk_x(3,rank) , wk_x(4,rank)
             do ix = wk_x(1,rank) , wk_x(2,rank)
                map_fftcd_x(id1*id2*(iz-1)+id1*(iy-1)+ix) = rank
             end do
          end do
       end do
    end do

    if ((div1*div2*div3) > myrank) then
       if (wk_x(1,myrank+1) /= 0) then
          np_fftcd_x = (wk_x(2,myrank+1)-wk_x(1,myrank+1)+1)* &
         &           (wk_x(4,myrank+1)-wk_x(3,myrank+1)+1)* &
         &           (wk_x(6,myrank+1)-wk_x(5,myrank+1)+1)
          allocate(mp_fftcd_x(np_fftcd_x), stat=ichkalloc)
       end if
       if(ichkalloc /= 0) then
          write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box_cd 3'
          call mpi_abort(mpi_comm_world, 5 , ichkalloc)
       endif
    else
       np_fftcd_x = 0
    end if

    deallocate(wk_x)

       if (np_fftcd_x /= 0) then
          id = 0
          do kk = xyz_fftcd_x(1,3), xyz_fftcd_x(2,3)
            do jj = xyz_fftcd_x(1,2), xyz_fftcd_x(2,2)
              do ii = xyz_fftcd_x(1,1), xyz_fftcd_x(2,1)
                id = id + 1
                mp_fftcd_x(id) = wk_fft_box(ii,jj,kk)
              enddo
            enddo
          enddo
       end if

    deallocate(wk_fft_box)

    allocate(nel_fftcd_x( 0:nmrank-1 ), stat=ichkalloc)
    if(ichkalloc /= 0) then
       write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box_cd 4'
       call mpi_abort(mpi_comm_world, 6 , ichkalloc)
    endif
    nel_fftcd_x(:) = 0

    do i = 1, nfft
      if (map_fftcd_x(i) .gt. 0) then
         nel_fftcd_x( map_fftcd_x(i)-1 ) = nel_fftcd_x( map_fftcd_x(i)-1 ) + 1
      endif
    enddo

    if (ipri > 1) then
      write(nfout,'("|||| m_Parallel_mpi_fft_box_cd")')
      write(nfout,'("|||| div1 = ",i5," , div2 = ",i5," , div3 = ",i5)')   div1,div2,div3
      write(nfout,'("|||| fft-box-size (",i4," ,",i4," ,",i4,")")') &
     &          ((fft_box_size_CD(i,j),i=1,3),j=0,0)
      write(nfout,'("|||| fft-box-div[X-axis] (",i4," :",i4," ,",i4," :",i4," ,",i4," :",i4,")")') &
     &          ((xyz_fftcd_x(i,j),i=1,2),j=1,3)
      write(nfout,'("|||| np_fftcd_x=",i8)')          np_fftcd_x
      write(nfout,'("|||| nel_fftcd_x:",/,"||||",10(i7))') (nel_fftcd_x(i),i=0,nmrank-1)
      write(nfout,'("|||| fftcd_X_x_dim=",i4,",  fftcd_X_y_dim=",i4,",  fftcd_X_z_dim=",i4)') &
     &                    fftcd_X_x_dim,fftcd_X_y_dim, fftcd_X_z_dim
      write(nfout, '("|||| nis_fftcd_X_x :")')
      write(nfout,'(("||||",10(i8)))') (nis_fftcd_X_x(i),i=1,fftcd_X_x_dim)
      write(nfout, '("|||| nie_ffcdt_X_x :")')
      write(nfout,'(("||||",10(i8)))') (nie_fftcd_X_x(i),i=1,fftcd_X_x_dim)
      write(nfout, '("|||| nis_fftcd_X_y :")')
      write(nfout,'(("||||",10(i8)))') (nis_fftcd_X_y(i),i=1,fftcd_X_y_dim)
      write(nfout, '("|||| nie_fftcd_X_y :")')
      write(nfout,'(("||||",10(i8)))') (nie_fftcd_X_y(i),i=1,fftcd_X_y_dim)
      write(nfout, '("|||| nis_fftcd_X_z :")')
      write(nfout,'(("||||",10(i8)))') (nis_fftcd_X_z(i),i=1,fftcd_X_z_dim)
      write(nfout, '("|||| nie_fftcd_X_z :")')
      write(nfout,'(("||||",10(i8)))') (nie_fftcd_X_z(i),i=1,fftcd_X_z_dim)
    endif

    color = nrank_g
!    if (myrank_g > (div1*div2*div3-1)) then
!       color = MPI_UNDEFINED
!    end if
    if( nrank_g>div1*div2*div3 ) then
      write(strng,'(i10)') nrank_g
      write(strngmax,'(i10)') div1*div2*div3
      call phase_error_with_msg(nfout,'ng = '//trim(strng)//', max ng ='//trim(strngmax), &
      __LINE__,__FILE__)
    endif
    key = myrank_g
!    if(firstcall) call mpi_comm_split(mpi_ke_world,color,key,mpi_fftcd_zy_world,ierr)
    if(.not.firstcall) call mpi_comm_free(mpi_fftcd_zy_world, ierr)
    call mpi_comm_split(mpi_ke_world,color,key,mpi_fftcd_zy_world,ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_mpi_fft_box_cd :  mpi_comm_split error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 100002, ierr)
     endif
     firstcall=.false.

#ifdef __TIMER_SUB__
  call timer_end(1233)
#endif

  end subroutine m_Parallel_mpi_fft_box_cd_3div
!===============================================================================
!===============================================================================
  subroutine m_Parallel_mpi_fft_box_cd_xyz(nfout,ipri,printable,fft_box_size_CD,kimg,fft_div1,fft_div2)
    implicit none

    integer, intent(in) :: nfout,ipri, kimg, fft_div1, fft_div2
    logical, intent(in) :: printable
    integer, intent(in) :: fft_box_size_CD(3,0:1)
    integer, allocatable, dimension(:,:) :: wk_z, wk_x, wk_y, wk1, wk2
    integer, allocatable, dimension(:,:,:) :: wk_fft_box
    integer :: dim1, dim2, igf, igf1, igf2, igf3
    integer :: i, j, id, id1, id2, id3, ierr, ichkalloc
    integer :: ii, jj, kk, ix, iy, iz
    integer :: nfft, n, nmrank, myrank
    integer :: mindim1, mindim2
    integer :: key, color
    logical, save :: firstcall=.true.
    character(len=12) :: strng,strngmax
                                                  __TIMER_SUB_START(1232)

    nmrank = nrank_g
    myrank = myrank_g

    mindim1 = 0
    mindim2 = 0
    if (((fft_div1*fft_div2) /= 0) .and. ((fft_div1*fft_div2) <= nmrank)) then
       dim1 = fft_div1
       dim2 = fft_div2
    else
       if (kimg == 1) then
          mindim1 = min(fft_box_size_CD(1,0)/2,fft_box_size_CD(2,0))
       else
          mindim1 = min(fft_box_size_CD(1,0),fft_box_size_CD(2,0))
       end if
       mindim2 = min(fft_box_size_CD(2,0),fft_box_size_CD(3,0))

       call decide_div_fft(nmrank, mindim1, mindim2, dim1, dim2)
    end if

    nfft = fft_box_size_CD(1,0)*fft_box_size_CD(2,0)*fft_box_size_CD(3,0)

    allocate(wk_z( 4,nmrank ), stat=ichkalloc)
    allocate(wk_x( 4,nmrank ), stat=ichkalloc)
    allocate(wk_y( 4,nmrank ), stat=ichkalloc)
    if(ichkalloc /= 0) then
       write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 1'
       call mpi_abort(mpi_comm_world, 3 , ichkalloc)
    endif
    wk_z(:,:) = 0
    wk_x(:,:) = 0
    wk_y(:,:) = 0

    allocate(wk1(2,dim1))
    allocate(wk2(2,dim2))
    wk1 = 0
    wk2 = 0

    fftcd_X_y_dim = dim1
    fftcd_X_z_dim = dim2
    fftcd_Y_x_dim = dim1
    fftcd_Y_z_dim = dim2
    fftcd_Z_x_dim = dim1
    fftcd_Z_y_dim = dim2

    allocate(nis_fftcd_X_z(fftcd_X_z_dim), stat=ichkalloc)
    allocate(nie_fftcd_X_z(fftcd_X_z_dim), stat=ichkalloc)
    allocate(nis_fftcd_X_y(fftcd_X_y_dim), stat=ichkalloc)
    allocate(nie_fftcd_X_y(fftcd_X_y_dim), stat=ichkalloc)
    allocate(nis_fftcd_Z_x(fftcd_Z_x_dim), stat=ichkalloc)
    allocate(nie_fftcd_Z_x(fftcd_Z_x_dim), stat=ichkalloc)
    allocate(nis_fftcd_Z_y(fftcd_Z_y_dim), stat=ichkalloc)
    allocate(nie_fftcd_Z_y(fftcd_Z_y_dim), stat=ichkalloc)
    allocate(nis_fftcd_Y_x(fftcd_Y_x_dim), stat=ichkalloc)
    allocate(nie_fftcd_Y_x(fftcd_Y_x_dim), stat=ichkalloc)
    allocate(nis_fftcd_Y_z(fftcd_Y_z_dim), stat=ichkalloc)
    allocate(nie_fftcd_Y_z(fftcd_Y_z_dim), stat=ichkalloc)

    call index_div_fft(kimg, wk1, dim1, fft_box_size_CD(2,0), 0)
    call index_div_fft(kimg, wk2, dim2, fft_box_size_CD(3,0), 0)
    nis_fftcd_X_y(:) = wk1(1,:)
    nie_fftcd_X_y(:) = wk1(2,:)
    nis_fftcd_X_z(:) = wk2(1,:)
    nie_fftcd_X_z(:) = wk2(2,:)
    ii = 1
    jj = 1
    do i = 1, dim1*dim2
       wk_x(1,i) = wk1(1,ii)
       wk_x(2,i) = wk1(2,ii)
       wk_x(3,i) = wk2(1,jj)
       wk_x(4,i) = wk2(2,jj)
       ii = ii + 1
       if (ii > dim1) then
          ii = 1
          jj = jj + 1
       end if
    end do

     call index_div_fft(kimg, wk1, dim1, fft_box_size_CD(1,0), 1)
    nis_fftcd_Y_x(:) = wk1(1,:)
    nie_fftcd_Y_x(:) = wk1(2,:)
    nis_fftcd_Y_z(:) = wk2(1,:)
    nie_fftcd_Y_z(:) = wk2(2,:)
    ii = 1
    jj = 1
    do i = 1, dim1*dim2
       wk_y(1,i) = wk1(1,ii)
       wk_y(2,i) = wk1(2,ii)
       wk_y(3,i) = wk2(1,jj)
       wk_y(4,i) = wk2(2,jj)
       ii = ii + 1
       if (ii > dim1) then
          ii = 1
          jj = jj + 1
       end if
    end do

    call index_div_fft(kimg, wk2, dim2, fft_box_size_CD(2,0), 0)
    nis_fftcd_Z_x(:) = wk1(1,:)
    nie_fftcd_Z_x(:) = wk1(2,:)
    nis_fftcd_Z_y(:) = wk2(1,:)
    nie_fftcd_Z_y(:) = wk2(2,:)
    ii = 1
    jj = 1
    do i = 1, dim1*dim2
       wk_z(1,i) = wk1(1,ii)
       wk_z(2,i) = wk1(2,ii)
       wk_z(3,i) = wk2(1,jj)
       wk_z(4,i) = wk2(2,jj)
       ii = ii + 1
       if (ii > dim1) then
          ii = 1
          jj = jj + 1
       end if
    end do

!!!!
!!  n = 1
!!!!
    n = 0
    if (wk_x(1,myrank+1) /= 0) then
      xyz_fftcd_x(1,1) = 1
      xyz_fftcd_x(2,1) = fft_box_size_CD(1,n)
      xyz_fftcd_x(1,2) = wk_x(1,myrank+1)
      xyz_fftcd_x(2,2) = wk_x(2,myrank+1)
      xyz_fftcd_x(1,3) = wk_x(3,myrank+1)
      xyz_fftcd_x(2,3) = wk_x(4,myrank+1)
    endif

    if (wk_y(1,myrank+1) /= 0) then
      xyz_fftcd_y(1,1) = wk_y(1,myrank+1)
      xyz_fftcd_y(2,1) = wk_y(2,myrank+1)
      xyz_fftcd_y(1,2) = 1
      xyz_fftcd_y(2,2) = fft_box_size_CD(2,n)
      xyz_fftcd_y(1,3) = wk_y(3,myrank+1)
      xyz_fftcd_y(2,3) = wk_y(4,myrank+1)
    endif

    if (wk_z(1,myrank+1) /= 0) then
      xyz_fftcd_z(1,1) = wk_z(1,myrank+1)
      xyz_fftcd_z(2,1) = wk_z(2,myrank+1)
      xyz_fftcd_z(1,2) = wk_z(3,myrank+1)
      xyz_fftcd_z(2,2) = wk_z(4,myrank+1)
      xyz_fftcd_z(1,3) = 1
      xyz_fftcd_z(2,3) = fft_box_size_CD(3,n)
    endif

    if(allocated(map_fftcd_x)) deallocate(map_fftcd_x)
    if(allocated(map_fftcd_z)) deallocate(map_fftcd_z)
    if(allocated(map_fftcd_y)) deallocate(map_fftcd_y)
    allocate(map_fftcd_x( nfft ), stat=ichkalloc)
    allocate(map_fftcd_z( nfft ), stat=ichkalloc)
    allocate(map_fftcd_y( nfft ), stat=ichkalloc)
    if(ichkalloc /= 0) then
       write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 2'
       call mpi_abort(mpi_comm_world, 4 , ichkalloc)
    endif
    map_fftcd_x(:) = 0
    map_fftcd_z(:) = 0
    map_fftcd_y(:) = 0

    np_fftcd_x = 0
    np_fftcd_z = 0
    np_fftcd_y = 0

    id1 = fft_box_size_CD(1,0)
    id2 = fft_box_size_CD(2,0)
    id3 = fft_box_size_CD(3,0)

    allocate( wk_fft_box(fft_box_size_CD(1,0), fft_box_size_CD(2,0), fft_box_size_CD(3,0)))
    do igf3 = 1, fft_box_size_CD(3,n)
       do igf2 = 1, fft_box_size_CD(2,n)
          do igf1 = 1, fft_box_size_CD(1,n)
             igf = igf1 + (igf2-1)*id1 + (igf3-1)*id1*id2
             wk_fft_box(igf1,igf2,igf3) = igf
          enddo
       enddo
    enddo

    do ii = 1, dim1*dim2
       if (wk_x(1,ii) == 0) cycle
       do iz = wk_x(3,ii) , wk_x(4,ii)
          do iy = wk_x(1,ii) , wk_x(2,ii)
             do ix = 1, id1
                map_fftcd_x(id1*id2*(iz-1)+id1*(iy-1)+ix) = ii
             end do
          end do
       end do
    end do

    do ii = 1, dim1*dim2
       if (wk_z(1,ii) == 0) cycle
       do iz = 1, id3
          do iy = wk_z(3,ii) , wk_z(4,ii)
             do ix = wk_z(1,ii), wk_z(2,ii)
                map_fftcd_z(id1*id2*(iz-1)+id1*(iy-1)+ix) = ii
             end do
          end do
       end do
    end do

    do ii = 1, dim1*dim2
       if (wk_y(1,ii) == 0) cycle
       do iz = wk_y(3,ii) , wk_y(4,ii)
          do iy = 1, id2
             do ix = wk_y(1,ii), wk_y(2,ii)
                map_fftcd_y(id1*id2*(iz-1)+id1*(iy-1)+ix) = ii
             end do
          end do
       end do
    end do

    if ((dim1*dim2) > myrank) then
       if (wk_x(1,myrank+1) /= 0) then
          np_fftcd_x = id1*(wk_x(2,myrank+1)-wk_x(1,myrank+1)+1)*(wk_x(4,myrank+1)-wk_x(3,myrank+1)+1)
          if(allocated(mp_fftcd_x)) deallocate(mp_fftcd_x)
          allocate(mp_fftcd_x(np_fftcd_x), stat=ichkalloc)
       end if
       if (wk_z(1,myrank+1) /= 0) then
          np_fftcd_z = (wk_z(4,myrank+1)-wk_z(3,myrank+1)+1)*(wk_z(2,myrank+1)-wk_z(1,myrank+1)+1)*id3
          if(allocated(mp_fftcd_z)) deallocate(mp_fftcd_z)
          allocate(mp_fftcd_z(np_fftcd_z), stat=ichkalloc)
       end if
       if (wk_y(1,myrank+1) /= 0) then
          np_fftcd_y = (wk_y(4,myrank+1)-wk_y(3,myrank+1)+1)*id2*(wk_y(2,myrank+1)-wk_y(1,myrank+1)+1)
          if(allocated(mp_fftcd_y)) deallocate(mp_fftcd_y)
          allocate(mp_fftcd_y(np_fftcd_y), stat=ichkalloc)
       end if
       if(ichkalloc /= 0) then
          write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 3'
          call mpi_abort(mpi_comm_world, 5 , ichkalloc)
       endif
    else
       np_fftcd_x = 0
       np_fftcd_z = 0
       np_fftcd_y = 0
    end if

    deallocate(wk_x)
    deallocate(wk_z)
    deallocate(wk_y)

    if (kimg == 1) then
       if (np_fftcd_x /= 0) then
          id = 0
          do kk = xyz_fftcd_x(1,3), xyz_fftcd_x(2,3)
            do jj = xyz_fftcd_x(1,2), xyz_fftcd_x(2,2)
              do ii = xyz_fftcd_x(1,1), xyz_fftcd_x(2,1) , 2
                id = id + 1
                mp_fftcd_x(id) = wk_fft_box(ii,jj,kk)
                id = id + 1
                mp_fftcd_x(id) = wk_fft_box(ii+1,jj,kk)
              enddo
            enddo
          enddo
       end if
       if (np_fftcd_y /= 0) then
          id = 0
            do jj = xyz_fftcd_y(1,2), xyz_fftcd_y(2,2)
          do kk = xyz_fftcd_y(1,3), xyz_fftcd_y(2,3)
              do ii = xyz_fftcd_y(1,1), xyz_fftcd_y(2,1) , 2
                id = id + 1
                mp_fftcd_y(id) = wk_fft_box(ii,jj,kk)
                id = id + 1
                mp_fftcd_y(id) = wk_fft_box(ii+1,jj,kk)
              enddo
            enddo
          enddo
       end if
       if (np_fftcd_z /= 0) then
          id = 0
#if 0
          do kk = xyz_fftcd_z(1,3), xyz_fftcd_z(2,3)
            do jj = xyz_fftcd_z(1,2), xyz_fftcd_z(2,2)
              do ii = xyz_fftcd_z(1,1), xyz_fftcd_z(2,1) , 2
#else
            do jj = xyz_fftcd_z(1,2), xyz_fftcd_z(2,2)
              do ii = xyz_fftcd_z(1,1), xyz_fftcd_z(2,1) , 2
          do kk = xyz_fftcd_z(1,3), xyz_fftcd_z(2,3)
#endif
                id = id + 1
                mp_fftcd_z(id) = wk_fft_box(ii,jj,kk)
                id = id + 1
                mp_fftcd_z(id) = wk_fft_box(ii+1,jj,kk)
              enddo
            enddo
          enddo
       end if
    else
       if (np_fftcd_x /= 0) then
          id = 0
          do kk = xyz_fftcd_x(1,3), xyz_fftcd_x(2,3)
            do jj = xyz_fftcd_x(1,2), xyz_fftcd_x(2,2)
              do ii = xyz_fftcd_x(1,1), xyz_fftcd_x(2,1)
                id = id + 1
                mp_fftcd_x(id) = wk_fft_box(ii,jj,kk)
              enddo
            enddo
          enddo
       end if
       if (np_fftcd_y /= 0) then
          id = 0
            do jj = xyz_fftcd_y(1,2), xyz_fftcd_y(2,2)
          do kk = xyz_fftcd_y(1,3), xyz_fftcd_y(2,3)
              do ii = xyz_fftcd_y(1,1), xyz_fftcd_y(2,1)
                id = id + 1
                mp_fftcd_y(id) = wk_fft_box(ii,jj,kk)
              enddo
            enddo
          enddo
       end if
       if (np_fftcd_z /= 0) then
          id = 0
#if 0
          do kk = xyz_fftcd_z(1,3), xyz_fftcd_z(2,3)
            do jj = xyz_fftcd_z(1,2), xyz_fftcd_z(2,2)
              do ii = xyz_fftcd_z(1,1), xyz_fftcd_z(2,1)
#else
            do jj = xyz_fftcd_z(1,2), xyz_fftcd_z(2,2)
              do ii = xyz_fftcd_z(1,1), xyz_fftcd_z(2,1)
          do kk = xyz_fftcd_z(1,3), xyz_fftcd_z(2,3)
#endif
                id = id + 1
                mp_fftcd_z(id) = wk_fft_box(ii,jj,kk)
              enddo
            enddo
          enddo
       end if
    endif

    deallocate(wk_fft_box)

    if(allocated(nel_fftcd_x)) deallocate(nel_fftcd_x)
    if(allocated(nel_fftcd_z)) deallocate(nel_fftcd_z)
    if(allocated(nel_fftcd_y)) deallocate(nel_fftcd_y)
    allocate(nel_fftcd_x( 0:nmrank-1 ), stat=ichkalloc)
    allocate(nel_fftcd_z( 0:nmrank-1 ), stat=ichkalloc)
    allocate(nel_fftcd_y( 0:nmrank-1 ), stat=ichkalloc)
    if(ichkalloc /= 0) then
       write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 4'
       call mpi_abort(mpi_comm_world, 6 , ichkalloc)
    endif
    nel_fftcd_z(:) = 0
    nel_fftcd_x(:) = 0
    nel_fftcd_y(:) = 0

    do i = 1, nfft
      if (map_fftcd_x(i) .gt. 0) then
         nel_fftcd_x( map_fftcd_x(i)-1 ) = nel_fftcd_x( map_fftcd_x(i)-1 ) + 1
      endif
      if (map_fftcd_z(i) .gt. 0) then
         nel_fftcd_z( map_fftcd_z(i)-1 ) = nel_fftcd_z( map_fftcd_z(i)-1 ) + 1
      endif
      if (map_fftcd_y(i) .gt. 0) then
         nel_fftcd_y( map_fftcd_y(i)-1 ) = nel_fftcd_y( map_fftcd_y(i)-1 ) + 1
      endif
    enddo

    if (ipri > 1) then
      write(nfout,'("|||| m_Parallel_mpi_fft_box")')
      write(nfout,'("|||| dim1 = ",i5," , dim2 = ",i5)')   dim1,dim2
      write(nfout,'("|||| fft-box-size (",i4," ,",i4," ,",i4,")")') &
     &          ((fft_box_size_CD(i,j),i=1,3),j=0,0)
      write(nfout,'("|||| fft-box-div[X-axis] (",i4," :",i4," ,",i4," :",i4," ,",i4," :",i4,")")') &
     &          ((xyz_fftcd_x(i,j),i=1,2),j=1,3)
      write(nfout,'("|||| fft-box-div[Y-axis] (",i4," :",i4," ,",i4," :",i4," ,",i4," :",i4,")")') &
     &          ((xyz_fftcd_y(i,j),i=1,2),j=1,3)
      write(nfout,'("|||| fft-box-div[Z-axis] (",i4," :",i4," ,",i4," :",i4," ,",i4," :",i4,")")') &
     &          ((xyz_fftcd_z(i,j),i=1,2),j=1,3)
!     write(nfout,'("|||| ir_x2z=",i4," , ir_z2y=",i4," , ir_y2z=",i4," , ir_z2x=",i4)') &
!    &                    ir_x2z,ir_z2y,ir_y2z,ir_z2x
!     write(nfout,'("|||| is_x2z=",i4," , is_z2y=",i4," , is_y2z=",i4," , is_z2x=",i4)') &
!    &                    is_x2z,is_z2y,is_y2z,is_z2x
      write(nfout,'("|||| np_fftcd_x=",i8)')          np_fftcd_x
      write(nfout,'("|||| nel_fftcd_x:",/,"||||",10(i7))') (nel_fftcd_x(i),i=0,nmrank-1)
      write(nfout,'("|||| np_fftcd_y=",i8)')          np_fftcd_y
      write(nfout,'("|||| nel_fftcd_y:",/,"||||",10(i7))') (nel_fftcd_y(i),i=0,nmrank-1)
      write(nfout,'("|||| np_fftcd_z=",i8)')          np_fftcd_z
      write(nfout,'("|||| nel_fftcd_z:",/,"||||",10(i7))') (nel_fftcd_z(i),i=0,nmrank-1)
      write(nfout, '("|||| fftcd_X_y_dim=",i4,",  fftcd_X_z_dim=",i4)') fftcd_X_y_dim, fftcd_X_z_dim
      write(nfout, '("|||| nis_fftcd_X_y :")')
      write(nfout,'(("||||",10(i8)))') (nis_fftcd_X_y(i),i=1,fftcd_X_y_dim)
      write(nfout, '("|||| nie_fftcd_X_y :")')
      write(nfout,'(("||||",10(i8)))') (nie_fftcd_X_y(i),i=1,fftcd_X_y_dim)
      write(nfout, '("|||| nis_fftcd_X_z :")')
      write(nfout,'(("||||",10(i8)))') (nis_fftcd_X_z(i),i=1,fftcd_X_z_dim)
      write(nfout, '("|||| nie_fftcd_X_z :")')
      write(nfout,'(("||||",10(i8)))') (nie_fftcd_X_z(i),i=1,fftcd_X_z_dim)
      write(nfout, '("|||| fftcd_Z_x_dim=",i4,",  fftcd_Z_y_dim=",i4)') fftcd_Z_x_dim, fftcd_Z_y_dim
      write(nfout, '("|||| nis_fftcd_Z_x :")')
      write(nfout,'(("||||",10(i8)))') (nis_fftcd_Z_x(i),i=1,fftcd_Z_x_dim)
      write(nfout, '("|||| nie_fftcd_Z_x :")')
      write(nfout,'(("||||",10(i8)))') (nie_fftcd_Z_x(i),i=1,fftcd_Z_x_dim)
      write(nfout, '("|||| nis_fftcd_Z_y :")')
      write(nfout,'(("||||",10(i8)))') (nis_fftcd_Z_y(i),i=1,fftcd_Z_y_dim)
      write(nfout, '("|||| nie_fftcd_Z_y :")')
      write(nfout,'(("||||",10(i8)))') (nie_fftcd_Z_y(i),i=1,fftcd_Z_y_dim)
      write(nfout, '("|||| fftcd_Y_x_dim=",i4,",  fftcd_Y_z_dim=",i4)') fftcd_Y_x_dim, fftcd_Y_z_dim
      write(nfout, '("|||| nis_fftcd_Y_x :")')
      write(nfout,'(("||||",10(i8)))') (nis_fftcd_Y_x(i),i=1,fftcd_Y_x_dim)
      write(nfout, '("|||| nie_fftcd_Y_x :")')
      write(nfout,'(("||||",10(i8)))') (nie_fftcd_Y_x(i),i=1,fftcd_Y_x_dim)
      write(nfout, '("|||| nis_fftcd_Y_z :")')
      write(nfout,'(("||||",10(i8)))') (nis_fftcd_Y_z(i),i=1,fftcd_Y_z_dim)
      write(nfout, '("|||| nie_fftcd_Y_z :")')
      write(nfout,'(("||||",10(i8)))') (nie_fftcd_Y_z(i),i=1,fftcd_Y_z_dim)
    endif

#ifdef FFT_ALLTOALL
    if (ipri > 1) then
      write(nfout,'("nrank_e=",i4,", myrank_e=",i4)') nrank_e, myrank_e
      write(nfout,'("nrank_g=",i4,", myrank_g=",i4)') nrank_g, myrank_g
      write(nfout,'("mpi_ke_world=",i12)') mpi_ke_world
      write(nfout,'("mpi_kg_world=",i12)') mpi_kg_world
      write(nfout,'("dim1=",i4,", dim2=",i4)') dim1, dim2
    endif

!   color = mod(myrank_g,dim2)
!   color = myrank_g/dim2
!   color = mod(myrank_g,dim1)
    color = myrank_g/dim1
!    if (myrank_g > (dim1*dim2-1)) then
!       color = MPI_UNDEFINED
!    end if
    if( nrank_g>dim1*dim2 ) then
      write(strng,'(i10)') nrank_g
      write(strngmax,'(i10)') dim1*dim2
      call phase_error_with_msg(nfout,'ng = '//trim(strng)//', max ng ='//trim(strngmax), &
      __LINE__,__FILE__)
    endif
    key = myrank_g
!    if(firstcall) call mpi_comm_split(mpi_ke_world,color,key,mpi_fftcd_xy_world,ierr)
    if(.not.firstcall) call mpi_comm_free(mpi_fftcd_xy_world, ierr)
    call mpi_comm_split(mpi_ke_world,color,key,mpi_fftcd_xy_world,ierr)
    call mpi_comm_size(mpi_fftcd_xy_world, nrank_fftcd_xy, ierr)
    call mpi_comm_rank(mpi_fftcd_xy_world, myrank_fftcd_xy, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_mpi_fft_box_cd :  mpi_comm_split error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 100002, ierr)
     endif

     call mpi_barrier(mpi_fftcd_xy_world, ierr)

    if (ipri > 1) then
      write(nfout,'(" mpi_fftcd_xy_world=",i12)') mpi_fftcd_xy_world
      write(nfout,'("  color=",i4,", key=",i4)') color, key
      write(nfout,'("  nrank_fftcd_xy=",i4,", myrank_fftcd_xy=",i4)') nrank_fftcd_xy, myrank_fftcd_xy
      call flush(nfout)
    end if

!   color = myrank_g/dim2
!   color = mod(myrank_g,dim2)
!   color = myrank_g/dim1
    color = mod(myrank_g,dim1)
!    if (myrank_g > (dim1*dim2-1)) then
!       color = MPI_UNDEFINED
!    end if
    if( nrank_g>dim1*dim2 ) then
      write(strng,*) nrank_g
      write(strngmax,*) dim1*dim2
      call phase_error_with_msg(nfout,'ng = '//trim(strng)//', max ng ='//trim(strngmax), &
      __LINE__,__FILE__)
    endif
    key = myrank_g
!    if(firstcall) call mpi_comm_split(mpi_ke_world,color,key,mpi_fftcd_yz_world,ierr)
    if(.not.firstcall) call mpi_comm_free(mpi_fftcd_yz_world, ierr)
    call mpi_comm_split(mpi_ke_world,color,key,mpi_fftcd_yz_world,ierr)
    call mpi_comm_size(mpi_fftcd_yz_world, nrank_fftcd_yz, ierr)
    call mpi_comm_rank(mpi_fftcd_yz_world, myrank_fftcd_yz, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_mpi_fft_boxcd :  mpi_comm_split error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 100002, ierr)
     endif

    if (ipri > 1) then
      write(nfout,'(" mpi_fftcd_yz_world=",i12)') mpi_fftcd_yz_world
      write(nfout,'("  color=",i4,", key=",i4)') color, key
      write(nfout,'("  nrank_fftcd_yz=",i4,", myrank_fftcd_yz=",i4)') nrank_fftcd_yz, myrank_fftcd_yz
      call flush(nfout)
    endif
#endif
    firstcall=.false.

                                                  __TIMER_SUB_STOP(1232)
  end subroutine m_Parallel_mpi_fft_box_cd_xyz
!===============================================================================
!===============================================================================
  subroutine m_Parallel_mpi_fft_box_cd(nfout,ipri,printable,fft_box_size_CD,kimg,divide)
    implicit none

    integer, intent(in) :: nfout,ipri, kimg, divide
    logical, intent(in) :: printable
    integer, intent(in) :: fft_box_size_CD(3,0:1)
    integer, allocatable, dimension(:,:,:) :: wk_fft_box
    integer, allocatable, dimension(:,:) :: wk_z, wk_x, wk_y, wk1, wk2
    integer :: dim1, dim2, igf, igf1, igf2, igf3, kdim1, kdim2
    integer :: i, j, id, id1, id2, id3, ierr, ichkalloc
    integer :: ii, jj, kk, ix, iy, iz, mindim1, mindim2
    integer :: nfft, n, nmrank, myrank, mpicom
    character(len=12) :: strng,strngmax
    logical, save :: firstcall=.true.
!   integer :: ir_x2z,ir_z2y,ir_y2z,ir_z2x
!   integer :: is_x2z,is_z2y,is_y2z,is_z2x
                                                  __TIMER_SUB_START(1233)
#ifdef CD_FFT_ALL
!   nmrank  = npes
!   myrank  = mype
!   mpicom  = mpi_comm_world
#else
    nmrank  = nrank_g
    myrank  = myrank_g
    mpicom  = mpi_ke_world
#endif

    if ((divide > 0) .and. (nmrank > 3)) then
       dim1 = int(sqrt(real(nmrank)))
       dim2 = dim1
    else
       mindim1 = min(fft_box_size_CD(2,0),fft_box_size_CD(3,0))
       if (kimg == 1) then
          mindim2 = min(fft_box_size_CD(1,0)/2,fft_box_size_CD(3,0))
       else
          mindim2 = min(fft_box_size_CD(1,0),fft_box_size_CD(3,0))
       end if

       call decide_div_fft(nmrank, mindim1, mindim2, dim2, dim1)
    end if

    nfft = fft_box_size_CD(1,0)*fft_box_size_CD(2,0)*fft_box_size_CD(3,0)

    allocate(wk_z( 4,nmrank ), stat=ichkalloc)
    allocate(wk_x( 4,nmrank ), stat=ichkalloc)
    allocate(wk_y( 4,nmrank ), stat=ichkalloc)
    if(ichkalloc /= 0) then
       write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 1'
       call mpi_abort(mpi_comm_world, 3 , ichkalloc)
    endif
    wk_z(:,:) = 0
    wk_x(:,:) = 0
    wk_y(:,:) = 0

    allocate(wk1(2,dim1))
    allocate(wk2(2,dim2))
    wk1 = 0
    wk2 = 0

    call index_div_fft(kimg, wk2, dim2, fft_box_size_CD(2,0), 0)
    ii = 0
    do i = 1, dim2
       if (wk2(1,i) == 0) cycle
       ii = ii + 1
    end do
    kdim2 = ii
    if (kdim2 /= dim2) then
       deallocate(wk2)
       dim2 = kdim2
       allocate(wk2(2,dim2))
    end if
    wk2 = 0

    call index_div_fft(kimg, wk1, dim1, fft_box_size_CD(3,0), 0)
    ii = 0
    do i = 1, dim1
       if (wk1(1,i) == 0) cycle
       ii = ii + 1
    end do
    kdim1 = ii

    if (kdim1 /= dim1) then
       deallocate(wk1)
       dim1 = kdim1
       allocate(wk1(2,dim1))
    end if
    wk1 = 0

    fftcd_X_z_dim = dim1
    fftcd_X_y_dim = dim2
    fftcd_Z_x_dim = dim1
    fftcd_Z_y_dim = dim2
    fftcd_Y_x_dim = dim1
    fftcd_Y_z_dim = dim2

    allocate(nis_fftcd_X_z(fftcd_X_z_dim), stat=ichkalloc)
    allocate(nie_fftcd_X_z(fftcd_X_z_dim), stat=ichkalloc)
    allocate(nis_fftcd_X_y(fftcd_X_y_dim), stat=ichkalloc)
    allocate(nie_fftcd_X_y(fftcd_X_y_dim), stat=ichkalloc)
    allocate(nis_fftcd_Z_x(fftcd_Z_x_dim), stat=ichkalloc)
    allocate(nie_fftcd_Z_x(fftcd_Z_x_dim), stat=ichkalloc)
    allocate(nis_fftcd_Z_y(fftcd_Z_y_dim), stat=ichkalloc)
    allocate(nie_fftcd_Z_y(fftcd_Z_y_dim), stat=ichkalloc)
    allocate(nis_fftcd_Y_x(fftcd_Y_x_dim), stat=ichkalloc)
    allocate(nie_fftcd_Y_x(fftcd_Y_x_dim), stat=ichkalloc)
    allocate(nis_fftcd_Y_z(fftcd_Y_z_dim), stat=ichkalloc)
    allocate(nie_fftcd_Y_z(fftcd_Y_z_dim), stat=ichkalloc)

    call index_div_fft(kimg, wk2, dim2, fft_box_size_CD(2,0), 0)
    call index_div_fft(kimg, wk1, dim1, fft_box_size_CD(3,0), 0)
    nis_fftcd_X_y(:) = wk2(1,:)
    nie_fftcd_X_y(:) = wk2(2,:)
    nis_fftcd_X_z(:) = wk1(1,:)
    nie_fftcd_X_z(:) = wk1(2,:)
    ii = 1
    jj = 1
    do i = 1, dim1*dim2
       wk_x(1,i) = wk2(1,ii)
       wk_x(2,i) = wk2(2,ii)
       wk_x(3,i) = wk1(1,jj)
       wk_x(4,i) = wk1(2,jj)
       ii = ii + 1
       if (ii > dim2) then
          ii = 1
          jj = jj + 1
       end if
    end do
    call index_div_fft(kimg, wk1, dim1, fft_box_size_CD(1,0), 1)
!   call index_div(wk2, dim2, fft_box_size_CD(2,0))
    nis_fftcd_Z_x(:) = wk1(1,:)
    nie_fftcd_Z_x(:) = wk1(2,:)
    nis_fftcd_Z_y(:) = wk2(1,:)
    nie_fftcd_Z_y(:) = wk2(2,:)
    ii = 1
    jj = 1
    do i = 1, dim1*dim2
       wk_z(1,i) = wk2(1,ii)
       wk_z(2,i) = wk2(2,ii)
       wk_z(3,i) = wk1(1,jj)
       wk_z(4,i) = wk1(2,jj)
       ii = ii + 1
       if (ii > dim2) then
          ii = 1
          jj = jj + 1
       end if
    end do
!   call index_div(wk1, dim1, fft_box_size_CD(1,0))
    call index_div_fft(kimg, wk2, dim2, fft_box_size_CD(3,0), 0)
    nis_fftcd_Y_x(:) = wk1(1,:)
    nie_fftcd_Y_x(:) = wk1(2,:)
    nis_fftcd_Y_z(:) = wk2(1,:)
    nie_fftcd_Y_z(:) = wk2(2,:)
    ii = 1
    jj = 1
    do i = 1, dim1*dim2
       wk_y(1,i) = wk2(1,ii)
       wk_y(2,i) = wk2(2,ii)
       wk_y(3,i) = wk1(1,jj)
       wk_y(4,i) = wk1(2,jj)
       ii = ii + 1
       if (ii > dim2) then
          ii = 1
          jj = jj + 1
       end if
    end do

!!!!
!!  n = 1
!!!!
    n = 0
!!  call div_index(fft_box_size_CD(2,n), fft_box_size_CD(3,n), dim1, dim2, wk_x, nmrank)
    if (wk_x(1,myrank+1) /= 0) then
      xyz_fftcd_x(1,1) = 1
      xyz_fftcd_x(2,1) = fft_box_size_CD(1,n)
      xyz_fftcd_x(1,2) = wk_x(1,myrank+1)
      xyz_fftcd_x(2,2) = wk_x(2,myrank+1)
      xyz_fftcd_x(1,3) = wk_x(3,myrank+1)
      xyz_fftcd_x(2,3) = wk_x(4,myrank+1)
    endif

!!  call div_index(fft_box_size_CD(2,n), fft_box_size_CD(1,n), dim1, dim2, wk_z, nmrank)
    if (wk_z(1,myrank+1) /= 0) then
      xyz_fftcd_z(1,1) = wk_z(3,myrank+1)
      xyz_fftcd_z(2,1) = wk_z(4,myrank+1)
      xyz_fftcd_z(1,2) = wk_z(1,myrank+1)
      xyz_fftcd_z(2,2) = wk_z(2,myrank+1)
      xyz_fftcd_z(1,3) = 1
      xyz_fftcd_z(2,3) = fft_box_size_CD(3,n)
    endif

!!  call div_index(fft_box_size_CD(3,n), fft_box_size_CD(1,n), dim1, dim2, wk_y, nmrank)
    if (wk_y(1,myrank+1) /= 0) then
      xyz_fftcd_y(1,1) = wk_y(3,myrank+1)
      xyz_fftcd_y(2,1) = wk_y(4,myrank+1)
      xyz_fftcd_y(1,2) = 1
      xyz_fftcd_y(2,2) = fft_box_size_CD(2,n)
      xyz_fftcd_y(1,3) = wk_y(1,myrank+1)
      xyz_fftcd_y(2,3) = wk_y(2,myrank+1)
    endif

    if(allocated(map_fftcd_x)) deallocate(map_fftcd_x)
    if(allocated(map_fftcd_z)) deallocate(map_fftcd_z)
    if(allocated(map_fftcd_y)) deallocate(map_fftcd_y)
    allocate(map_fftcd_x( nfft ), stat=ichkalloc)
    allocate(map_fftcd_z( nfft ), stat=ichkalloc)
    allocate(map_fftcd_y( nfft ), stat=ichkalloc)
    if(ichkalloc /= 0) then
       write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 2'
       call mpi_abort(mpi_comm_world, 4 , ichkalloc)
    endif
    map_fftcd_x(:) = 0
    map_fftcd_z(:) = 0
    map_fftcd_y(:) = 0

    np_fftcd_x = 0
    np_fftcd_z = 0
    np_fftcd_y = 0

    id1 = fft_box_size_CD(1,0)
    id2 = fft_box_size_CD(2,0)
    id3 = fft_box_size_CD(3,0)

    allocate( wk_fft_box(fft_box_size_CD(1,0), fft_box_size_CD(2,0), fft_box_size_CD(3,0)))
    do igf3 = 1, fft_box_size_CD(3,n)
       do igf2 = 1, fft_box_size_CD(2,n)
          do igf1 = 1, fft_box_size_CD(1,n)
             igf = igf1 + (igf2-1)*id1 + (igf3-1)*id1*id2

             wk_fft_box(igf1,igf2,igf3) = igf

!            do j = 1, nmrank
!              if(((igf2.ge.wk_x(1,j)).and.(igf2.le.wk_x(2,j))) .and. &
!             &   ((igf3.ge.wk_x(3,j)).and.(igf3.le.wk_x(4,j)))) then
!                map_fftcd_x(igf) = j
!                if (myrank .eq. (j-1)) np_fftcd_x = np_fftcd_x + 1
!              endif
!              if(((igf1.ge.wk_z(3,j)).and.(igf1.le.wk_z(4,j))) .and. &
!             &   ((igf2.ge.wk_z(1,j)).and.(igf2.le.wk_z(2,j)))) then
!                map_fftcd_z(igf) = j
!                if (myrank .eq. (j-1)) np_fftcd_z = np_fftcd_z + 1
!              endif
!              if(((igf1.ge.wk_y(3,j)).and.(igf1.le.wk_y(4,j))) .and. &
!             &   ((igf3.ge.wk_y(1,j)).and.(igf3.le.wk_y(2,j)))) then
!                map_fftcd_y(igf) = j
!                if (myrank .eq. (j-1)) np_fftcd_y = np_fftcd_y + 1
!              endif
!            enddo

          enddo
       enddo
    enddo

!x  do ii = 1, nmrank
    do ii = 1, dim1*dim2
       if (wk_x(1,ii) == 0) cycle
       do iz = wk_x(3,ii) , wk_x(4,ii)
          do iy = wk_x(1,ii) , wk_x(2,ii)
             do ix = 1, id1
                map_fftcd_x(id1*id2*(iz-1)+id1*(iy-1)+ix) = ii
             end do
          end do
       end do
    end do

!x  do ii = 1, nmrank
    do ii = 1, dim1*dim2
       if (wk_z(1,ii) == 0) cycle
       do iz = 1, id3
          do iy = wk_z(1,ii) , wk_z(2,ii)
             do ix = wk_z(3,ii), wk_z(4,ii)
                map_fftcd_z(id1*id2*(iz-1)+id1*(iy-1)+ix) = ii
             end do
          end do
       end do
    end do

!x  do ii = 1, nmrank
    do ii = 1, dim1*dim2
       if (wk_y(1,ii) == 0) cycle
       do iz = wk_y(1,ii) , wk_y(2,ii)
          do iy = 1, id2
             do ix = wk_y(3,ii), wk_y(4,ii)
                map_fftcd_y(id1*id2*(iz-1)+id1*(iy-1)+ix) = ii
             end do
          end do
       end do
    end do
    if ((dim1*dim2) > myrank) then
       if (wk_x(1,myrank+1) /= 0) then
          np_fftcd_x = id1*(wk_x(2,myrank+1)-wk_x(1,myrank+1)+1)*(wk_x(4,myrank+1)-wk_x(3,myrank+1)+1)
          if(allocated(mp_fftcd_x)) deallocate(mp_fftcd_x)
          allocate(mp_fftcd_x(np_fftcd_x), stat=ichkalloc)
       end if
       if (wk_z(1,myrank+1) /= 0) then
          np_fftcd_z = (wk_z(4,myrank+1)-wk_z(3,myrank+1)+1)*(wk_z(2,myrank+1)-wk_z(1,myrank+1)+1)*id3
          if(allocated(mp_fftcd_z)) deallocate(mp_fftcd_z)
          allocate(mp_fftcd_z(np_fftcd_z), stat=ichkalloc)
       end if
       if (wk_y(1,myrank+1) /= 0) then
          np_fftcd_y = (wk_y(4,myrank+1)-wk_y(3,myrank+1)+1)*id2*(wk_y(2,myrank+1)-wk_y(1,myrank+1)+1)
          if(allocated(mp_fftcd_y)) deallocate(mp_fftcd_y)
          allocate(mp_fftcd_y(np_fftcd_y), stat=ichkalloc)
       end if
       if(ichkalloc /= 0) then
          write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 3'
          call mpi_abort(mpi_comm_world, 5 , ichkalloc)
       endif
    else
       np_fftcd_x = 0
       np_fftcd_z = 0
       np_fftcd_y = 0
    end if

    deallocate(wk_x)
    deallocate(wk_z)
    deallocate(wk_y)

    if (kimg == 1) then
       if (np_fftcd_x /= 0) then
          id = 0
          do kk = xyz_fftcd_x(1,3), xyz_fftcd_x(2,3)
            do jj = xyz_fftcd_x(1,2), xyz_fftcd_x(2,2)
              do ii = xyz_fftcd_x(1,1), xyz_fftcd_x(2,1) , 2
                id = id + 1
                mp_fftcd_x(id) = wk_fft_box(ii,jj,kk)
                id = id + 1
                mp_fftcd_x(id) = wk_fft_box(ii+1,jj,kk)
              enddo
            enddo
          enddo
       end if
       if (np_fftcd_y /= 0) then
          id = 0
            do jj = xyz_fftcd_y(1,2), xyz_fftcd_y(2,2)
          do kk = xyz_fftcd_y(1,3), xyz_fftcd_y(2,3)
              do ii = xyz_fftcd_y(1,1), xyz_fftcd_y(2,1) , 2
                id = id + 1
                mp_fftcd_y(id) = wk_fft_box(ii,jj,kk)
                id = id + 1
                mp_fftcd_y(id) = wk_fft_box(ii+1,jj,kk)
              enddo
            enddo
          enddo
       end if
       if (np_fftcd_z /= 0) then
          id = 0
          do kk = xyz_fftcd_z(1,3), xyz_fftcd_z(2,3)
            do jj = xyz_fftcd_z(1,2), xyz_fftcd_z(2,2)
              do ii = xyz_fftcd_z(1,1), xyz_fftcd_z(2,1) , 2
                id = id + 1
                mp_fftcd_z(id) = wk_fft_box(ii,jj,kk)
                id = id + 1
                mp_fftcd_z(id) = wk_fft_box(ii+1,jj,kk)
              enddo
            enddo
          enddo
       end if
    else
       if (np_fftcd_x /= 0) then
          id = 0
          do kk = xyz_fftcd_x(1,3), xyz_fftcd_x(2,3)
            do jj = xyz_fftcd_x(1,2), xyz_fftcd_x(2,2)
              do ii = xyz_fftcd_x(1,1), xyz_fftcd_x(2,1)
                id = id + 1
                mp_fftcd_x(id) = wk_fft_box(ii,jj,kk)
              enddo
            enddo
          enddo
       end if
       if (np_fftcd_y /= 0) then
          id = 0
            do jj = xyz_fftcd_y(1,2), xyz_fftcd_y(2,2)
          do kk = xyz_fftcd_y(1,3), xyz_fftcd_y(2,3)
              do ii = xyz_fftcd_y(1,1), xyz_fftcd_y(2,1)
                id = id + 1
                mp_fftcd_y(id) = wk_fft_box(ii,jj,kk)
              enddo
            enddo
          enddo
       end if
       if (np_fftcd_z /= 0) then
          id = 0
          do kk = xyz_fftcd_z(1,3), xyz_fftcd_z(2,3)
            do jj = xyz_fftcd_z(1,2), xyz_fftcd_z(2,2)
              do ii = xyz_fftcd_z(1,1), xyz_fftcd_z(2,1)
                id = id + 1
                mp_fftcd_z(id) = wk_fft_box(ii,jj,kk)
              enddo
            enddo
          enddo
       end if
    endif

    deallocate(wk_fft_box)

    if(allocated(nel_fftcd_x)) deallocate(nel_fftcd_x)
    if(allocated(nel_fftcd_z)) deallocate(nel_fftcd_z)
    if(allocated(nel_fftcd_y)) deallocate(nel_fftcd_y)
    allocate(nel_fftcd_x( 0:nmrank-1 ), stat=ichkalloc)
    allocate(nel_fftcd_z( 0:nmrank-1 ), stat=ichkalloc)
    allocate(nel_fftcd_y( 0:nmrank-1 ), stat=ichkalloc)
    if(ichkalloc /= 0) then
       write(nfout,*)'could not allocate in m_Parallel_mpi_fft_box 4'
       call mpi_abort(mpi_comm_world, 6 , ichkalloc)
    endif
    nel_fftcd_z(:) = 0
    nel_fftcd_x(:) = 0
    nel_fftcd_y(:) = 0

    do i = 1, nfft
      if (map_fftcd_x(i) .gt. 0) then
         nel_fftcd_x( map_fftcd_x(i)-1 ) = nel_fftcd_x( map_fftcd_x(i)-1 ) + 1
      endif
      if (map_fftcd_z(i) .gt. 0) then
         nel_fftcd_z( map_fftcd_z(i)-1 ) = nel_fftcd_z( map_fftcd_z(i)-1 ) + 1
      endif
      if (map_fftcd_y(i) .gt. 0) then
         nel_fftcd_y( map_fftcd_y(i)-1 ) = nel_fftcd_y( map_fftcd_y(i)-1 ) + 1
      endif
    enddo

    if (ipri > 1) then
      write(nfout,'("|||| m_Parallel_mpi_fft_boxcd")')
      write(nfout,'("|||| dim1 = ",i5," , dim2 = ",i5)')   dim1,dim2
      write(nfout,'("|||| fft-box-size (",i4," ,",i4," ,",i4,")")') &
     &          ((fft_box_size_CD(i,j),i=1,3),j=0,0)
      write(nfout,'("|||| fft-box-div[X-axis] (",i4," :",i4," ,",i4," :",i4," ,",i4," :",i4,")")') &
     &          ((xyz_fftcd_x(i,j),i=1,2),j=1,3)
      write(nfout,'("|||| fft-box-div[Z-axis] (",i4," :",i4," ,",i4," :",i4," ,",i4," :",i4,")")') &
     &          ((xyz_fftcd_z(i,j),i=1,2),j=1,3)
      write(nfout,'("|||| fft-box-div[Y-axis] (",i4," :",i4," ,",i4," :",i4," ,",i4," :",i4,")")') &
     &          ((xyz_fftcd_y(i,j),i=1,2),j=1,3)
!     write(nfout,'("|||| ir_x2z=",i4," , ir_z2y=",i4," , ir_y2z=",i4," , ir_z2x=",i4)') &
!    &                    ir_x2z,ir_z2y,ir_y2z,ir_z2x
!     write(nfout,'("|||| is_x2z=",i4," , is_z2y=",i4," , is_y2z=",i4," , is_z2x=",i4)') &
!    &                    is_x2z,is_z2y,is_y2z,is_z2x
      write(nfout,'("|||| np_fftcd_x=",i8)')          np_fftcd_x
      write(nfout,'("|||| nel_fftcd_x:",/,"||||",10(i7))') (nel_fftcd_x(i),i=0,nmrank-1)
      write(nfout,'("|||| np_fftcd_y=",i8)')          np_fftcd_y
      write(nfout,'("|||| nel_fftcd_y:",/,"||||",10(i7))') (nel_fftcd_y(i),i=0,nmrank-1)
      write(nfout,'("|||| np_fftcd_z=",i8)')          np_fftcd_z
      write(nfout,'("|||| nel_fftcd_z:",/,"||||",10(i7))') (nel_fftcd_z(i),i=0,nmrank-1)
    endif

#ifdef FFT_ALLTOALL
    color = mod(myrank_g,dim2)
!    if (myrank_g > (dim1*dim2-1)) then
!       color = MPI_UNDEFINED
!    end if
    if( nrank_g>dim1*dim2 ) then
      write(strng,'(i10)') nrank_g
      write(strngmax,'(i10)') dim1*dim2
      call phase_error_with_msg(nfout,'ng = '//trim(strng)//', max ng ='//trim(strngmax), &
      __LINE__,__FILE__)
    endif
    key = myrank_g
!    if(firstcall) call mpi_comm_split(mpi_ke_world,color,key,mpi_fftcd_xz_world,ierr)
    if(.not.firstcall) call mpi_comm_free(mpi_fftcd_xz_world, ierr)
    call mpi_comm_split(mpi_ke_world,color,key,mpi_fftcd_xz_world,ierr)
    call mpi_comm_size(mpi_fftcd_xz_world, nrank_fftcd_xz, ierr)
    call mpi_comm_rank(mpi_fftcd_xz_world, myrank_fftcd_xz, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_mpi_fftcd_box :  mpi_comm_split error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 100002, ierr)
     endif

     call mpi_barrier(mpi_fftcd_xz_world, ierr)

    if (ipri > 1) then
! === DEBUG by tkato 2012/06/06 ================================================
!     write(nfout,'(" mpi_fftcd_xz_world=",i8)') mpi_fftcd_xz_world
      write(nfout,'(" mpi_fftcd_xz_world=",i12)') mpi_fftcd_xz_world
! ==============================================================================
      write(nfout,'("  color=",i4,", key=",i4)') color, key
      write(nfout,'("  nrank_fftcd_xz=",i4,", myrank_fftcd_xz=",i4)') nrank_fftcd_xz, myrank_fftcd_xz
      call flush(nfout)
    endif

    color = myrank_g/dim2
!    if (myrank_g > (dim1*dim2-1)) then
!       color = MPI_UNDEFINED
!    end if
    if( nrank_g>dim1*dim2 ) then
      write(strng,*) nrank_g
      write(strngmax,*) dim1*dim2
      call phase_error_with_msg(nfout,'ng = '//trim(strng)//', max ng ='//trim(strngmax), &
      __LINE__,__FILE__)
    endif
    key = myrank_g
!    if(firstcall) call mpi_comm_split(mpi_ke_world,color,key,mpi_fftcd_zy_world,ierr)
    if(.not.firstcall) call mpi_comm_free(mpi_fftcd_zy_world, ierr)
    call mpi_comm_split(mpi_ke_world,color,key,mpi_fftcd_zy_world,ierr)
    call mpi_comm_size(mpi_fftcd_zy_world, nrank_fftcd_zy, ierr)
    call mpi_comm_rank(mpi_fftcd_zy_world, myrank_fftcd_zy, ierr)
     if (ierr /= 0) then
! === DEBUG by tkato 2012/06/06 ================================================
!     write(nfout,'(" mpi_fftcd_zy_world=",i8)') mpi_fftcd_zy_world
      write(nfout,'(" mpi_fftcd_zy_world=",i12)') mpi_fftcd_zy_world
! ==============================================================================
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 100002, ierr)
     endif

    if (ipri > 1) then
      write(nfout,'(" mpi_fftcd_zy_world=",i12)') mpi_fftcd_zy_world
      write(nfout,'("  color=",i4,", key=",i4)') color, key
      write(nfout,'("  nrank_fftcd_zy=",i4,", myrank_fftcd_zy=",i4)') nrank_fftcd_zy, myrank_fftcd_zy
      call flush(nfout)
    endif
#endif
    firstcall=.false.

                                                  __TIMER_SUB_STOP(1233)
  end subroutine m_Parallel_mpi_fft_box_cd
!===============================================================================

!===============================================================================
!!$
!!$  subroutine m_Parallel_init_mpi_atm_3D(nfout,ipri,printable,natm)
!!$    integer, intent(in) :: nfout,ipri,natm
!!$    logical, intent(in) :: printable
!!$    integer :: iwork, i, npes, mype
!!$#ifdef __TIMER_SUB__
!!$  call timer_sta(1238)
!!$#endif
!!$
!!$    npes = nrank_g
!!$    mype = myrank_g
!!$    allocate(is_atm(0:npes-1))
!!$    allocate(ie_atm(0:npes-1))
!!$    allocate(nel_atm(0:npes-1))
!!$    iwork = ( natm - 1 ) / npes + 1
!!$    if(ipri >= 1 .and. printable) then
!!$       write(nfout,'(" !|| << init_mpi_atm >>")')
!!$       write(nfout,'(" !|| natm = ",i12)') natm
!!$       write(nfout,'(" !|| -- is_natm, ie_natm --")')
!!$    end if
!!$    do i = 0, npes-1
!!$       is_atm(i) = min(i*iwork+1, natm+1)
!!$       ie_atm(i) = min(is_atm(i)+iwork-1, natm)
!!$       nel_atm(i) = ie_atm(i) - is_atm(i) + 1
!!$       if(ipri >= 1 .and. printable) write(nfout,'(" !|| ",2i12)') is_atm(i),ie_atm(i)
!!$    enddo
!!$    ista_atm = is_atm(mype)
!!$    iend_atm = ie_atm(mype)
!!$    np_atm   = nel_atm(mype)
!!$    mp_atm   = maxval(nel_atm)
!!$#ifdef __TIMER_SUB__
!!$  call timer_end(1238)
!!$#endif
!!$  end subroutine m_Parallel_init_mpi_atm_3D

!===============================================================================

  subroutine m_Parallel_init_mpi_atm_B_3D(nfout,ipri,printable,natm,comm_for_chg)
    integer, intent(in) :: nfout,ipri,natm
    logical, intent(in) :: printable
    logical, intent(in) :: comm_for_chg
    integer :: iwork, i, npes, mype, iadd, iblock, ishift
                                                  __TIMER_SUB_START(1240)
    integer :: npes_,mype_
    if(comm_for_chg)then
    npes_ = 1
    mype_ = 0
    else
    npes_ = nrank_e
    mype_ = myrank_e
    endif
    allocate(is_atm_B(0:npes_-1))
    allocate(ie_atm_B(0:npes_-1))
    allocate(nel_atm_B(0:npes_-1))
    iwork = ( natm - 1 ) / npes_ + 1
    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| << init_mpi_atm >>")')
       write(nfout,'(" !|| natm = ",i12)') natm
       write(nfout,'(" !|| -- is_natm, ie_natm --")')
    end if
    do i = 0, npes_-1
       is_atm_B(i) = min(i*iwork+1, natm+1)
       ie_atm_B(i) = min(is_atm_B(i)+iwork-1, natm)
       nel_atm_B(i) = ie_atm_B(i) - is_atm_B(i) + 1
       if(ipri >= 1 .and. printable) write(nfout,'(" !|| ",2i12)') is_atm_B(i),ie_atm_B(i)
    enddo
    ista_atm_B = is_atm_B(mype_)
    iend_atm_B = ie_atm_B(mype_)
    np_atm_B   = nel_atm_B(mype_)
    mp_atm_B   = maxval(nel_atm_B)
    allocate(mem_atm_B(1:natm))
    iblock = 0
    ishift = 1
    do i = 1, natm
       iadd = ishift + iwork * iblock
       if (iadd > natm) then
          iblock = 0
          ishift = ishift + 1
          iadd = ishift + iwork * iblock
       end if
       iblock = iblock + 1
       mem_atm_B(iadd) = i
    end do
                                                  __TIMER_SUB_STOP(1240)
  end subroutine m_Parallel_init_mpi_atm_B_3D
!===============================================================================

#ifdef FFT_3D_DIVISION
!===============================================================================
  subroutine m_Parallel_wf_onto_fft_3D(nfout,fft_box_size_WF,igf,nbase,nbase_gamma, &
 &                                     k_symmetry,GAMMA,kg,kg_gamma,kv3,kk)
    integer, intent(in)  :: nfout, kg, kg_gamma, kv3, kk, GAMMA
    integer, intent(in)  :: fft_box_size_WF(1:3,0:1)
    integer, intent(in)  :: igf(kg)
    integer, intent(in)  :: nbase(kg1_ext,kv3)
    integer, intent(in)  :: nbase_gamma(kg_gamma,2)
    integer, intent(in)  :: k_symmetry(kv3)
    integer, dimension(0:nrank_g-1)                       ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)       ::sta_r, sta_s

    integer, allocatable, dimension(:,:,:) :: xyz
    integer, allocatable, dimension(:,:,:) :: work
    integer :: max_fft_x, lx, ly, lz, mx, my, mz, mm, i, j, i1, klen, ik, max_np_g1k, max_mp_g1k
    integer :: iadd, ladd, ista, len, itag = 10
    integer :: kx1p, kx2p, kx3p
                                                  __TIMER_SUB_START(1243)

     max_fft_x = maxval(nel_fft_x(:))
     allocate(xyz(2,3,0:nrank_g-1))
#ifdef USE_NONBLK_COMM
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1352)
     req_r = 0
     req_s = 0
     sta_r = 0
     sta_s = 0
     do i = 0, nrank_g - 1
        call mpi_irecv(xyz(1,1,i), 6, mpi_integer, &
       &               i, itag, mpi_ke_world, req_r(i), ierr)
         if (ierr /= 0) then
            write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_irecv error'
            call flush(nfout)
            call mpi_abort(mpi_comm_world,170,ierr)
         endif
     enddo
     do i = 0, nrank_g - 1
        call mpi_isend(xyz_fft_x, 6, mpi_integer, &
       &               i, itag, mpi_ke_world, req_s(i), ierr)
         if (ierr /= 0) then
            write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_isend error'
            call flush(nfout)
            call mpi_abort(mpi_comm_world,171,ierr)
         endif
     enddo
     call mpi_waitall(nrank_g, req_r, sta_r, ierr)
      if (ierr /= 0) then
         write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_waitall error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world,172,ierr)
      endif
     call mpi_waitall(nrank_g, req_s, sta_s, ierr)
      if (ierr /= 0) then
         write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_waitall error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world,173,ierr)
      endif
                                                  __TIMER_INICOMM_STOP(1352)
#else
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1353)
     integer, allocatable, dimension(:,:) :: rbuf
     allocate(rbuf(6,0:nrank_g-1))
     call MPI_ALLGATHER(xyz_fft_x, 6, mpi_integer, &
    &                   rbuf,6, mpi_integer, mpi_ke_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_allgather error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 174, ierr)
     endif
     do i = 1,2
       do j = 1,3
         xyz(i,j,:)=rbuf(i+(j-1)*2,:)
       enddo
     enddo
     deallocate(rbuf)
                                                  __TIMER_INICOMM_STOP(1353)
#endif

     lx = fft_box_size_WF(1,0)
     ly = fft_box_size_WF(2,0)
     lz = fft_box_size_WF(3,0)

     kx1p = fft_X_x_nel
     kx2p = fft_X_y_nel
     kx3p = fft_X_z_nel

     allocate(wf_fft_scnt(0:nrank_g-1,ista_k:iend_k))
     allocate(wf_fft_rcnt(0:nrank_g-1,ista_k:iend_k))
     len = 1
     do ik = ista_k, iend_k
        if(k_symmetry(ik) == GAMMA) then
           len = 2
        end if
     end do
     max_np_g1k = maxval(np_g1k(:))
     max_mp_g1k = maxval(mp_g1k(:))
     allocate(wf_fft_index(max_np_g1k*len,ista_k:iend_k))
     allocate(wf_fft_dist (max_np_g1k*len,ista_k:iend_k))
     allocate(wf_fft_send (max_mp_g1k*len,ista_k:iend_k,0:nrank_g-1))
     allocate(wf_fft_recv (max_mp_g1k*len,ista_k:iend_k,0:nrank_g-1))
     allocate(wf_fft_maxsend(ista_k:iend_k))
     allocate(wf_fft_maxrecv(ista_k:iend_k))
     wf_fft_dist(:,:) = -1
     wf_fft_send(:,:,:) = 0
     wf_fft_recv(:,:,:) = 0
     wf_fft_scnt(:,:) = 0
     wf_fft_rcnt(:,:) = 0
     wf_fft_maxsend(:) = 0
     wf_fft_maxrecv(:) = 0

     klen = iend_k - ista_k + 1

     do ik = ista_k, iend_k

        if(k_symmetry(ik) == GAMMA) then
           ista = ista_g1k(ik)
           if (ista == 1) then
              iadd = 1
              i1 = igf(1)
              mz = (i1-1)/(lx*ly)+1
              mm = mod(i1,(lx*ly))
              if (mm==0) mm=lx*ly
!!!           my = (mm-1)/ly+1
!!!           mx = mod(mm,ly)
              my = (mm-1)/lx+1
              mx = mod(mm,lx)
              if (mx==0) mx = ly
              B_4 : do i = 0, nrank_g-1
                 if ((xyz(1,1,i)<=mx).and.(mx<=xyz(2,1,i)).and.    &
                &    (xyz(1,2,i)<=my).and.(my<=xyz(2,2,i)).and.    &
                &    (xyz(1,3,i)<=mz).and.(mz<=xyz(2,3,i))) then
                    ladd = mx-xyz(1,1,i)+1+kx1p*(my-xyz(1,2,i))+kx1p*kx2p*(mz-xyz(1,3,i))
                    wf_fft_scnt(i,ik) = wf_fft_scnt(i,ik) + 1
                    wf_fft_index(iadd*2-1,ik) = wf_fft_scnt(i,ik)
                    wf_fft_dist (iadd*2-1,ik) = i
                    wf_fft_send(wf_fft_scnt(i,ik),ik,i) = ladd
                    wf_fft_scnt(i,ik) = wf_fft_scnt(i,ik) + 1
                    wf_fft_index(iadd*2  ,ik) = wf_fft_scnt(i,ik)
                    wf_fft_dist (iadd*2  ,ik) = i
                    wf_fft_send(wf_fft_scnt(i,ik),ik,i) = ladd
                    exit
                 endif
              enddo B_4
              ista = 2
           endif
           B_1 : do j = ista, iend_g1k(ik)
              iadd = j-ista_g1k(ik)+1
              i1 = igf(nbase(j,ik))
              mz = (i1-1)/(lx*ly)+1
              mm = mod(i1,(lx*ly))
              if (mm==0) mm=lx*ly
              my = (mm-1)/lx+1
              mx = mod(mm,lx)
              if (mx==0) mx = lx
              B_2 : do i = 0, nrank_g-1
                 if ((xyz(1,1,i)<=mx).and.(mx<=xyz(2,1,i)).and.    &
                &    (xyz(1,2,i)<=my).and.(my<=xyz(2,2,i)).and.    &
                &    (xyz(1,3,i)<=mz).and.(mz<=xyz(2,3,i))) then
                    ladd = mx-xyz(1,1,i)+1+kx1p*(my-xyz(1,2,i))+kx1p*kx2p*(mz-xyz(1,3,i))
                       wf_fft_scnt(i,ik) = wf_fft_scnt(i,ik) + 1
                    wf_fft_index(iadd*2-1,ik) = wf_fft_scnt(i,ik)
                    wf_fft_dist (iadd*2-1,ik) = i
                    wf_fft_send(wf_fft_scnt(i,ik),ik,i) = ladd
                    exit
                 endif
              enddo B_2
              i1 = igf(nbase_gamma(j,2))
              mz = (i1-1)/(lx*ly)+1
              mm = mod(i1,(lx*ly))
              if (mm==0) mm=lx*ly
              my = (mm-1)/lx+1
              mx = mod(mm,lx)
              if (mx==0) mx = lx
              B_3 : do i = 0, nrank_g-1
                 if ((xyz(1,1,i)<=mx).and.(mx<=xyz(2,1,i)).and.    &
                   & (xyz(1,2,i)<=my).and.(my<=xyz(2,2,i)).and.    &
                   & (xyz(1,3,i)<=mz).and.(mz<=xyz(2,3,i))) then
                    ladd = mx-xyz(1,1,i)+1+kx1p*(my-xyz(1,2,i))+kx1p*kx2p*(mz-xyz(1,3,i))
                    wf_fft_scnt(i,ik) = wf_fft_scnt(i,ik) + 1
                    wf_fft_index(iadd*2  ,ik) = wf_fft_scnt(i,ik)
                    wf_fft_dist (iadd*2  ,ik) = i
                    wf_fft_send(wf_fft_scnt(i,ik),ik,i) = ladd
                    exit
                 endif
              enddo B_3
           enddo B_1
        else
           B_11 : do j = ista_g1k(ik), iend_g1k(ik)
              iadd = j-ista_g1k(ik)+1
              i1 = igf(nbase(j,ik))
              mz = (i1-1)/(lx*ly)+1
                 mm = mod(i1,(lx*ly))
              if (mm==0) mm=lx*ly
              my = (mm-1)/lx+1
              mx = mod(mm,lx)
              if (mx==0) mx = lx
              B_12 : do i = 0, nrank_g-1
                 if ((xyz(1,1,i)<=mx).and.(mx<=xyz(2,1,i)).and.    &
                &    (xyz(1,2,i)<=my).and.(my<=xyz(2,2,i)).and.    &
                &    (xyz(1,3,i)<=mz).and.(mz<=xyz(2,3,i))) then
                    ladd = mx-xyz(1,1,i)+1+kx1p*(my-xyz(1,2,i))+kx1p*kx2p*(mz-xyz(1,3,i))
                    wf_fft_scnt(i,ik) = wf_fft_scnt(i,ik) + 1
                    wf_fft_index(iadd,ik) = wf_fft_scnt(i,ik)
                    wf_fft_dist (iadd,ik) = i
                    wf_fft_send(wf_fft_scnt(i,ik),ik,i) = ladd
                    exit
                 endif
              enddo B_12
           enddo B_11
        endif

     end do
     deallocate(xyz)

#ifdef USE_NONBLK_COMM
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1354)
     req_r = 0
     req_s = 0
     sta_r = 0
     sta_s = 0
     do i = 0, nrank_g - 1
        call mpi_irecv(wf_fft_recv(1,ista_k,i), max_mp_g1k*len*klen, mpi_integer, &
       &               i, itag, mpi_ke_world, req_r(i), ierr)
         if (ierr /= 0) then
            write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_irecv error'
            call flush(nfout)
            call mpi_abort(mpi_comm_world,174,ierr)
         endif
     enddo
     do i = 0, nrank_g - 1
        call mpi_isend(wf_fft_send(1,ista_k,i), max_mp_g1k*len*klen, mpi_integer, &
       &               i, itag, mpi_ke_world, req_s(i), ierr)
         if (ierr /= 0) then
            write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_isend error'
            call flush(nfout)
            call mpi_abort(mpi_comm_world,175,ierr)
         endif
     enddo
     call mpi_waitall(nrank_g, req_r, sta_r, ierr)
      if (ierr /= 0) then
         write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_waitall error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world,176,ierr)
    endif
     call mpi_waitall(nrank_g, req_s, sta_s, ierr)
      if (ierr /= 0) then
         write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_waitall error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world,177,ierr)
      endif
                                                  __TIMER_INICOMM_STOP(1354)
#else
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1355)
     call MPI_ALLTOALL(wf_fft_send(:,ista_k,:), max_mp_g1k*len*klen, mpi_integer, &
    &                  wf_fft_recv(:,ista_k,:), max_mp_g1k*len*klen, mpi_integer, &
    &                                          mpi_ke_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 178, ierr)
     endif
                                                  __TIMER_INICOMM_STOP(1355)
#endif

     do i = 0, nrank_g - 1
        do ik = ista_k, iend_k
           do j = 1, mp_g1k(ik)*len
              if (wf_fft_recv(j,ik,i) == 0) then
                 exit
              end if
              wf_fft_rcnt(i,ik) = wf_fft_rcnt(i,ik) + 1
           enddo
        enddo
     end do

     do ik = ista_k, iend_k
        wf_fft_maxsend(ik) = maxval(wf_fft_scnt(0:nrank_g-1,ik))
        wf_fft_maxrecv(ik) = maxval(wf_fft_rcnt(0:nrank_g-1,ik))
     end do

     allocate(work(maxval(wf_fft_maxrecv(:)),ista_k:iend_k,0:nrank_g-1))
     do i = 0, nrank_g - 1
        do ik = ista_k, iend_k
           do j = 1, wf_fft_rcnt(i,ik)
              work(j,ik,i) = wf_fft_recv(j,ik,i)
           end do
        end do
     end do
     deallocate(wf_fft_recv)
     allocate(wf_fft_recv(maxval(wf_fft_maxrecv(:)),ista_k:iend_k,0:nrank_g-1))
     do i = 0, nrank_g - 1
        do ik = ista_k, iend_k
           do j = 1, wf_fft_rcnt(i,ik)
              wf_fft_recv(j,ik,i) = work(j,ik,i)
           end do
        end do
     end do
     deallocate(work)

     deallocate(wf_fft_send)
                                                  __TIMER_SUB_START(1243)
  end subroutine m_Parallel_wf_onto_fft_3D
!===============================================================================
#else
  subroutine m_Parallel_wf_onto_fft_3D(nfout,fft_box_size_WF,igf,nbase,nbase_gamma, &
 &                                     k_symmetry,GAMMA,kg,kg_gamma,kv3,kk)
    integer, intent(in)  :: nfout, kg, kg_gamma, kv3, kk, GAMMA
    integer, intent(in)  :: fft_box_size_WF(1:3,0:1)
    integer, intent(in)  :: igf(kg)
    integer, intent(in)  :: nbase(kg1_ext,kv3)
    integer, intent(in)  :: nbase_gamma(kg_gamma,2)
    integer, intent(in)  :: k_symmetry(kv3)
    integer, dimension(0:nrank_g-1)                       ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)       ::sta_r, sta_s

    integer, allocatable, dimension(:,:,:) :: xyz
    integer, allocatable, dimension(:,:,:) :: work
    integer :: max_fft_x, lx, ly, lz, mx, my, mz, mm, i, j, i1, klen, ik, max_np_g1k, max_mp_g1k
    integer :: iadd, ladd, ista, len, itag = 10
! === DEBUG by tkato 2012/06/04 ================================================
#ifndef USE_NONBLK_COMM
     integer, allocatable, dimension(:,:) :: rbuf
#endif
! ==============================================================================
                                                  __TIMER_SUB_STOP(1243)

     max_fft_x = maxval(nel_fft_x(:))
     allocate(xyz(2,3,0:nrank_g-1))
#ifdef USE_NONBLK_COMM
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1352)
     req_r = 0
     req_s = 0
     sta_r = 0
     sta_s = 0
     do i = 0, nrank_g - 1
        call mpi_irecv(xyz(1,1,i), 6, mpi_integer, &
       &               i, itag, mpi_ke_world, req_r(i), ierr)
         if (ierr /= 0) then
            write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_irecv error'
            call flush(nfout)
            call mpi_abort(mpi_comm_world,170,ierr)
         endif
     enddo
     do i = 0, nrank_g - 1
        call mpi_isend(xyz_fft_x, 6, mpi_integer, &
       &               i, itag, mpi_ke_world, req_s(i), ierr)
         if (ierr /= 0) then
            write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_isend error'
            call flush(nfout)
            call mpi_abort(mpi_comm_world,171,ierr)
         endif
     enddo
     call mpi_waitall(nrank_g, req_r, sta_r, ierr)
      if (ierr /= 0) then
         write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_waitall error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world,172,ierr)
      endif
     call mpi_waitall(nrank_g, req_s, sta_s, ierr)
      if (ierr /= 0) then
         write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_waitall error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world,173,ierr)
      endif
                                                  __TIMER_INICOMM_STOP(1352)
#else
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1353)
! === DEBUG by tkato 2012/06/04 ================================================
!    integer, allocatable, dimension(:,:) :: rbuf
! ==============================================================================
     allocate(rbuf(6,0:nrank_g-1))
     call MPI_ALLGATHER(xyz_fft_x, 6, mpi_integer, &
    &                   rbuf,6, mpi_integer, mpi_ke_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_allgather error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 174, ierr)
     endif
     do i = 1,2
       do j = 1,3
         xyz(i,j,:)=rbuf(i+(j-1)*2,:)
       enddo
     enddo
     deallocate(rbuf)
                                                  __TIMER_INICOMM_STOP(1353)
#endif

     lx = fft_box_size_WF(1,0)
     ly = fft_box_size_WF(2,0)
     lz = fft_box_size_WF(3,0)
     allocate(wf_fft_scnt(0:nrank_g-1,ista_k:iend_k))
     allocate(wf_fft_rcnt(0:nrank_g-1,ista_k:iend_k))
     len = 1
     do ik = ista_k, iend_k
        if(k_symmetry(ik) == GAMMA) then
           len = 2
        end if
     end do
     max_np_g1k = maxval(np_g1k(:))
     max_mp_g1k = maxval(mp_g1k(:))
     allocate(wf_fft_index(max_np_g1k*len,ista_k:iend_k))
     allocate(wf_fft_dist (max_np_g1k*len,ista_k:iend_k))
     allocate(wf_fft_send (max_mp_g1k*len,ista_k:iend_k,0:nrank_g-1))
     allocate(wf_fft_recv (max_mp_g1k*len,ista_k:iend_k,0:nrank_g-1))
     allocate(wf_fft_maxsend(ista_k:iend_k))
     allocate(wf_fft_maxrecv(ista_k:iend_k))
     wf_fft_dist(:,:) = -1
     wf_fft_send(:,:,:) = 0
     wf_fft_recv(:,:,:) = 0
     wf_fft_scnt(:,:) = 0
     wf_fft_rcnt(:,:) = 0
     wf_fft_maxsend(:) = 0
     wf_fft_maxrecv(:) = 0

     klen = iend_k - ista_k + 1

     do ik = ista_k, iend_k

        if(k_symmetry(ik) == GAMMA) then
           ista = ista_g1k(ik)
           if (ista == 1) then
              iadd = 1
              i1 = igf(1)
              mz = (i1-1)/(lx*ly)+1
              mm = mod(i1,(lx*ly))
              if (mm==0) mm=lx*ly
!!!           my = (mm-1)/ly+1
!!!           mx = mod(mm,ly)
              my = (mm-1)/lx+1
              mx = mod(mm,lx)
              if (mx==0) mx = ly
              B_4 : do i = 0, nrank_g-1
                 if ((xyz(1,2,i)<=my).and.(my<=xyz(2,2,i)).and.    &
                &    (xyz(1,3,i)<=mz).and.(mz<=xyz(2,3,i))) then
                    ladd = mx+lx*(my-xyz(1,2,i))+lx*(xyz(2,2,i)-xyz(1,2,i)+1)*(mz-xyz(1,3,i))
                    wf_fft_scnt(i,ik) = wf_fft_scnt(i,ik) + 1
                    wf_fft_index(iadd*2-1,ik) = wf_fft_scnt(i,ik)
                    wf_fft_dist (iadd*2-1,ik) = i
                    wf_fft_send(wf_fft_scnt(i,ik),ik,i) = ladd
                    wf_fft_scnt(i,ik) = wf_fft_scnt(i,ik) + 1
                    wf_fft_index(iadd*2  ,ik) = wf_fft_scnt(i,ik)
                    wf_fft_dist (iadd*2  ,ik) = i
                    wf_fft_send(wf_fft_scnt(i,ik),ik,i) = ladd
                    exit
                 endif
              enddo B_4
              ista = 2
           endif
           B_1 : do j = ista, iend_g1k(ik)
              iadd = j-ista_g1k(ik)+1
              i1 = igf(nbase(j,ik))
              mz = (i1-1)/(lx*ly)+1
              mm = mod(i1,(lx*ly))
              if (mm==0) mm=lx*ly
              my = (mm-1)/lx+1
              mx = mod(mm,lx)
              if (mx==0) mx = lx
              B_2 : do i = 0, nrank_g-1
                 if ((xyz(1,2,i)<=my).and.(my<=xyz(2,2,i)).and.    &
                &    (xyz(1,3,i)<=mz).and.(mz<=xyz(2,3,i))) then
                    ladd = mx+lx*(my-xyz(1,2,i))+lx*(xyz(2,2,i)-xyz(1,2,i)+1)*(mz-xyz(1,3,i))
                       wf_fft_scnt(i,ik) = wf_fft_scnt(i,ik) + 1
                    wf_fft_index(iadd*2-1,ik) = wf_fft_scnt(i,ik)
                    wf_fft_dist (iadd*2-1,ik) = i
                    wf_fft_send(wf_fft_scnt(i,ik),ik,i) = ladd
                    exit
                 endif
              enddo B_2
              i1 = igf(nbase_gamma(j,2))
              mz = (i1-1)/(lx*ly)+1
              mm = mod(i1,(lx*ly))
              if (mm==0) mm=lx*ly
              my = (mm-1)/lx+1
              mx = mod(mm,lx)
              if (mx==0) mx = lx
              B_3 : do i = 0, nrank_g-1
                 if ((xyz(1,2,i)<=my).and.(my<=xyz(2,2,i)).and.    &
                   &    (xyz(1,3,i)<=mz).and.(mz<=xyz(2,3,i))) then
                    ladd = mx+lx*(my-xyz(1,2,i))+lx*(xyz(2,2,i)-xyz(1,2,i)+1)*(mz-xyz(1,3,i))
                    wf_fft_scnt(i,ik) = wf_fft_scnt(i,ik) + 1
                    wf_fft_index(iadd*2  ,ik) = wf_fft_scnt(i,ik)
                    wf_fft_dist (iadd*2  ,ik) = i
                    wf_fft_send(wf_fft_scnt(i,ik),ik,i) = ladd
                    exit
                 endif
              enddo B_3
           enddo B_1
        else
           B_11 : do j = ista_g1k(ik), iend_g1k(ik)
              iadd = j-ista_g1k(ik)+1
              i1 = igf(nbase(j,ik))
              mz = (i1-1)/(lx*ly)+1
                 mm = mod(i1,(lx*ly))
              if (mm==0) mm=lx*ly
              my = (mm-1)/lx+1
              mx = mod(mm,lx)
              if (mx==0) mx = lx
              B_12 : do i = 0, nrank_g-1
                 if ((xyz(1,2,i)<=my).and.(my<=xyz(2,2,i)).and.    &
                &    (xyz(1,3,i)<=mz).and.(mz<=xyz(2,3,i))) then
                    ladd = mx+lx*(my-xyz(1,2,i))+lx*(xyz(2,2,i)-xyz(1,2,i)+1)*(mz-xyz(1,3,i))
                    wf_fft_scnt(i,ik) = wf_fft_scnt(i,ik) + 1
                    wf_fft_index(iadd,ik) = wf_fft_scnt(i,ik)
                    wf_fft_dist (iadd,ik) = i
                    wf_fft_send(wf_fft_scnt(i,ik),ik,i) = ladd
                    exit
                 endif
              enddo B_12
           enddo B_11
        endif

     end do
     deallocate(xyz)

#ifdef USE_NONBLK_COMM
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1354)
     req_r = 0
     req_s = 0
     sta_r = 0
     sta_s = 0
     do i = 0, nrank_g - 1
        call mpi_irecv(wf_fft_recv(1,ista_k,i), max_mp_g1k*len*klen, mpi_integer, &
       &               i, itag, mpi_ke_world, req_r(i), ierr)
         if (ierr /= 0) then
            write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_irecv error'
            call flush(nfout)
            call mpi_abort(mpi_comm_world,174,ierr)
         endif
     enddo
     do i = 0, nrank_g - 1
        call mpi_isend(wf_fft_send(1,ista_k,i), max_mp_g1k*len*klen, mpi_integer, &
       &               i, itag, mpi_ke_world, req_s(i), ierr)
         if (ierr /= 0) then
            write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_isend error'
            call flush(nfout)
            call mpi_abort(mpi_comm_world,175,ierr)
         endif
     enddo
     call mpi_waitall(nrank_g, req_r, sta_r, ierr)
      if (ierr /= 0) then
         write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_waitall error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world,176,ierr)
    endif
     call mpi_waitall(nrank_g, req_s, sta_s, ierr)
      if (ierr /= 0) then
         write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_waitall error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world,177,ierr)
      endif
                                                  __TIMER_INICOMM_STOP(1354)
#else
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1355)
! === DEBUG by tkato 2012/06/05 ================================================
!    call MPI_ALLTOALL(wf_fft_send(:,ista_k,:), max_mp_g1k*len*klen, mpi_integer, &
!   &                  wf_fft_recv(:,ista_k,:), max_mp_g1k*len*klen, mpi_integer, &
!   &                                          mpi_ke_world, ierr )
     call MPI_ALLTOALL(wf_fft_send, max_mp_g1k*len*klen, mpi_integer, &
    &                  wf_fft_recv, max_mp_g1k*len*klen, mpi_integer, &
    &                                          mpi_ke_world, ierr )
! ==============================================================================
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 178, ierr)
     endif
                                                  __TIMER_INICOMM_STOP(1355)
#endif

     do i = 0, nrank_g - 1
        do ik = ista_k, iend_k
           do j = 1, mp_g1k(ik)*len
              if (wf_fft_recv(j,ik,i) == 0) then
                 exit
              end if
              wf_fft_rcnt(i,ik) = wf_fft_rcnt(i,ik) + 1
           enddo
        enddo
     end do

     do ik = ista_k, iend_k
        wf_fft_maxsend(ik) = maxval(wf_fft_scnt(0:nrank_g-1,ik))
        wf_fft_maxrecv(ik) = maxval(wf_fft_rcnt(0:nrank_g-1,ik))
     end do

     allocate(work(maxval(wf_fft_maxrecv(:)),ista_k:iend_k,0:nrank_g-1))
     do i = 0, nrank_g - 1
        do ik = ista_k, iend_k
           do j = 1, wf_fft_rcnt(i,ik)
              work(j,ik,i) = wf_fft_recv(j,ik,i)
           end do
        end do
     end do
     deallocate(wf_fft_recv)
     allocate(wf_fft_recv(maxval(wf_fft_maxrecv(:)),ista_k:iend_k,0:nrank_g-1))
     do i = 0, nrank_g - 1
        do ik = ista_k, iend_k
           do j = 1, wf_fft_rcnt(i,ik)
              wf_fft_recv(j,ik,i) = work(j,ik,i)
           end do
        end do
     end do
     deallocate(work)

     deallocate(wf_fft_send)
                                                  __TIMER_SUB_STOP(1243)
  end subroutine m_Parallel_wf_onto_fft_3D
!===============================================================================
#endif

!===============================================================================
  subroutine m_Parallel_wf_onto_fft_dealloc_3D
     deallocate(wf_fft_rcnt)
     deallocate(wf_fft_scnt)
     deallocate(wf_fft_recv)
!    deallocate(wf_fft_send)
     deallocate(wf_fft_index)
     deallocate(wf_fft_dist)
     deallocate(wf_fft_maxrecv)
     deallocate(wf_fft_maxsend)
  end subroutine m_Parallel_wf_onto_fft_dealloc_3D

!===============================================================================
  subroutine m_Parallel_fft_onto_wf_3D(nfout,fft_box_size_WF,igf,nbase,kg,kv3,nfft,kk)
    integer, intent(in)  :: nfout, kg, kv3, nfft, kk
    integer, intent(in)  :: igf(kg)
    integer, intent(in)  :: nbase(kg1_ext,kv3)
    integer, intent(in)  :: fft_box_size_WF(1:3,0:1)
    integer, dimension(0:nrank_g-1)                       ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)       ::sta_r, sta_s

    integer, allocatable,dimension(:,:,:) :: fftigf
    integer, allocatable, dimension(:,:,:) :: work
    integer :: i1, lrank, i, j, k, lsize, isrsize,fft_l_size, klen, ik
    integer, parameter :: itag = 10
                                                  __TIMER_SUB_START(1244)
!   lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
    lsize = maxval(nel_fft_x(:))
    isrsize = min(lsize,maxval(mp_g1k(:)))
    fft_l_size  = nel_fft_x(myrank_g)

    allocate(fft_wf_scnt(0:nrank_g-1,ista_k:iend_k), stat=ierr)
    allocate(fft_wf_rcnt(0:nrank_g-1,ista_k:iend_k), stat=ierr)
    allocate(fft_wf_send(isrsize,ista_k:iend_k,0:nrank_g-1), stat=ierr)
    allocate(fft_wf_recv(isrsize,ista_k:iend_k,0:nrank_g-1), stat=ierr)
    allocate(fft_wf_dist(fft_l_size,ista_k:iend_k), stat=ierr)
    allocate(fft_wf_index(fft_l_size,ista_k:iend_k), stat=ierr)
    allocate(fft_wf_maxsend(ista_k:iend_k), stat=ierr)
    allocate(fft_wf_maxrecv(ista_k:iend_k), stat=ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fft_onto_wf_3D :  Not allocate '
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 204, ierr)
     endif
    fft_wf_scnt(:,:) = 0
    fft_wf_rcnt(:,:) = 0
    fft_wf_send(:,:,:) = 0
    fft_wf_recv(:,:,:) = 0
    fft_wf_dist(:,:) = -1
    fft_wf_index(:,:) = 0

    klen = iend_k - ista_k + 1

    allocate(fftigf(nfft,2,ista_k:iend_k), stat=ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fft_onto_wf_3D :  Not allocate '
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 206, ierr)
     endif

    do ik = ista_k, iend_k

       fftigf(:,1,ik) = -1

       do i = 0, nrank_g - 1
          do j = nis_g1k(i,ik), nie_g1k(i,ik)
             i1 = igf(nbase(j,ik))
             if (i1 > nfft) cycle
             fftigf(i1,1,ik) = i
             fftigf(i1,2,ik) = j - nis_g1k(i,ik) + 1
          enddo
       enddo
       do k = 1, nel_fft_x(myrank_g)

          i1 = mp_fft_x(k)
          if (fftigf(i1,1,ik) < 0) cycle
          lrank = fftigf(i1,1,ik)
          fft_wf_scnt(lrank,ik) = fft_wf_scnt(lrank,ik) + 1
          fft_wf_send(fft_wf_scnt(lrank,ik),ik,lrank) = fftigf(i1,2,ik)
          fft_wf_dist(k,ik) = lrank
          fft_wf_index(k,ik) = fft_wf_scnt(lrank,ik)

       end do
    end do
    deallocate(fftigf)

#ifdef USE_NONBLK_COMM
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1356)
    lrank = mod(myrank_g,nrank_g)
    do i = 0, nrank_g - 1
       lrank = lrank + 1
       if (lrank > (nrank_g - 1)) lrank = 0
       call mpi_irecv(fft_wf_recv(1,ista_k,lrank), isrsize*klen, mpi_integer, &
      &               lrank, itag, mpi_ke_world, req_r(lrank), ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Direct_3D :  mpi_irecv error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 207, ierr)
        endif
    enddo

    lrank = mod((myrank_g+1),nrank_g)
    do i = 0, nrank_g - 1
       lrank = lrank + 1
       if (lrank > (nrank_g - 1)) lrank = 0
       call mpi_isend(fft_wf_send(1,ista_k,lrank), isrsize*klen, mpi_integer, &
      &               lrank, itag, mpi_ke_world, req_s(lrank), ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_Parallel_fft_onto_wf_3D :  mpi_isend error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 208, ierr)
        endif
    enddo

    call mpi_waitall(nrank_g, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fft_onto_wf_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 209, ierr)
     endif
    call mpi_waitall(nrank_g, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fft_onto_wf_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 210, ierr)
     endif
                                                  __TIMER_INICOMM_STOP(1356)
#else
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1357)
! === DEBUG by tkato 2012/06/05 ================================================
!    call MPI_ALLTOALL(fft_wf_send(:,ista_k,:), isrsize*klen, mpi_integer, &
!   &                  fft_wf_recv(:,ista_k,:), isrsize*klen, mpi_integer, &
!   &                                             mpi_ke_world, ierr )
     call MPI_ALLTOALL(fft_wf_send, isrsize*klen, mpi_integer, &
    &                  fft_wf_recv, isrsize*klen, mpi_integer, &
    &                                             mpi_ke_world, ierr )
! ==============================================================================
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fft_onto_wf_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 211, ierr)
     endif
                                                  __TIMER_INICOMM_STOP(1357)
#endif

    do i = 0, nrank_g - 1
       do ik = ista_k, iend_k
          do j = 1, isrsize
             if (fft_wf_recv(j,ik,i) == 0) then
                exit
             end if
             fft_wf_rcnt(i,ik) = fft_wf_rcnt(i,ik) + 1
          enddo
       enddo
    end do

    do ik = ista_k, iend_k
       fft_wf_maxsend(ik) = maxval(fft_wf_scnt(:,ik))
       fft_wf_maxrecv(ik) = maxval(fft_wf_rcnt(:,ik))
    end do

     allocate(work(maxval(fft_wf_maxrecv(:)),ista_k:iend_k,0:nrank_g-1))
     do i = 0, nrank_g - 1
        do ik = ista_k, iend_k
           do j = 1, fft_wf_rcnt(i,ik)
              work(j,ik,i) = fft_wf_recv(j,ik,i)
           end do
        end do
     end do
     deallocate(fft_wf_recv)
     allocate(fft_wf_recv(maxval(fft_wf_maxrecv(:)),ista_k:iend_k,0:nrank_g-1))
     do i = 0, nrank_g - 1
        do ik = ista_k, iend_k
           do j = 1, fft_wf_rcnt(i,ik)
              fft_wf_recv(j,ik,i) = work(j,ik,i)
           end do
        end do
     end do
     deallocate(work)
     deallocate(fft_wf_send)

                                                  __TIMER_SUB_STOP(1244)
  end subroutine m_Parallel_fft_onto_wf_3D

!===============================================================================

  subroutine m_Parallel_fft_onto_wf_dealloc_3D
     if(allocated(fft_wf_rcnt)) deallocate(fft_wf_rcnt)
     if(allocated(fft_wf_scnt)) deallocate(fft_wf_scnt)
     if(allocated(fft_wf_recv)) deallocate(fft_wf_recv)
!    deallocate(fft_wf_send)
     if(allocated(fft_wf_index)) deallocate(fft_wf_index)
     if(allocated(fft_wf_dist)) deallocate(fft_wf_dist)
     if(allocated(fft_wf_maxsend)) deallocate(fft_wf_maxsend)
     if(allocated(fft_wf_maxrecv)) deallocate(fft_wf_maxrecv)
  end subroutine m_Parallel_fft_onto_wf_dealloc_3D

!===============================================================================

#ifdef MPI_FFTW
  subroutine m_Parallel_wf_onto_fft_mpifftw(nfout,fft_box_size_WF,igf,nbase,nbase_gamma, &
 &                                     k_symmetry,GAMMA,kg,kg_gamma,kv3,kk)
     integer, intent(in)  :: nfout, kg, kg_gamma, kv3, kk, GAMMA
     integer, intent(in)  :: fft_box_size_WF(1:3,0:1)
     integer, intent(in)  :: igf(kg)
     integer, intent(in)  :: nbase(kg1_ext,kv3)
     integer, intent(in)  :: nbase_gamma(kg_gamma,2)
     integer, intent(in)  :: k_symmetry(kv3)
     integer, dimension(0:nrank_g-1)                       ::req_r,req_s
     integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)       ::sta_r, sta_s

     integer, allocatable, dimension(:,:,:) :: work
     integer :: mx, my, mz, mm, i, j, i1, klen, ik, max_np_g1k, max_mp_g1k
     integer :: iadd, ladd, ista, len, itag = 10

     integer(C_INTPTR_T)  :: local_n, local_n_offset, alloc_local, lx, ly, lz
     integer, allocatable, dimension(:) :: local_ns, local_n_offsets

     call fftw_mpi_init()

     lx = fft_box_size_WF(1,0)
     ly = fft_box_size_WF(2,0)
     lz = fft_box_size_WF(3,0)

     alloc_local = fftw_mpi_local_size_3d(lz,ly,lx,mpi_ke_world,local_n,local_n_offset)

     allocate(local_ns(0:nrank_g-1));local_ns=0
     allocate(local_n_offsets(0:nrank_g-1));local_n_offsets=0
     local_ns(myrank_g) = local_n
     local_n_offsets(myrank_g) = local_n_offset
     call mpi_allreduce(mpi_in_place, local_ns, nrank_g, mpi_integer, mpi_sum, mpi_ke_world, ierr)
     call mpi_allreduce(mpi_in_place, local_n_offsets, nrank_g, mpi_integer, mpi_sum, mpi_ke_world, ierr)
     allocate(wf_fft_scnt_mfftw(0:nrank_g-1,ista_k:iend_k))
     allocate(wf_fft_rcnt_mfftw(0:nrank_g-1,ista_k:iend_k))
     len = 1
     do ik = ista_k, iend_k
        if(k_symmetry(ik) == GAMMA) then
           len = 2
        end if
     end do
     max_np_g1k = maxval(np_g1k(:))
     max_mp_g1k = maxval(mp_g1k(:))
     allocate(wf_fft_index_mfftw(max_np_g1k*len,ista_k:iend_k))
     allocate(wf_fft_dist_mfftw (max_np_g1k*len,ista_k:iend_k))
     allocate(wf_fft_send_mfftw (max_mp_g1k*len,ista_k:iend_k,0:nrank_g-1))
     allocate(wf_fft_recv_mfftw (max_mp_g1k*len,ista_k:iend_k,0:nrank_g-1))
     allocate(wf_fft_maxsend_mfftw(ista_k:iend_k))
     allocate(wf_fft_maxrecv_mfftw(ista_k:iend_k))
     wf_fft_dist_mfftw(:,:) = -1
     wf_fft_send_mfftw(:,:,:) = 0
     wf_fft_recv_mfftw(:,:,:) = 0
     wf_fft_scnt_mfftw(:,:) = 0
     wf_fft_rcnt_mfftw(:,:) = 0
     wf_fft_maxsend_mfftw(:) = 0
     wf_fft_maxrecv_mfftw(:) = 0

     klen = iend_k - ista_k + 1

     do ik = ista_k, iend_k

        if(k_symmetry(ik) == GAMMA) then
           ista = ista_g1k(ik)
           if (ista == 1) then
              iadd = 1
              i1 = igf(1)
              mz = (i1-1)/(lx*ly)+1
              mm = mod(i1,(lx*ly))
              if (mm==0) mm=lx*ly
!!!           my = (mm-1)/ly+1
!!!           mx = mod(mm,ly)
              my = (mm-1)/lx+1
              mx = mod(mm,lx)
              if (mx==0) mx = ly
              B_4 : do i = 0, nrank_g-1
                 if (mz>local_n_offsets(i) .and. mz<=local_n_offsets(i)+local_ns(i)) then
                    ladd = mx+lx*(my-1)+lx*ly*(mz-local_n_offsets(i)-1)
                    wf_fft_scnt_mfftw(i,ik) = wf_fft_scnt_mfftw(i,ik) + 1
                    wf_fft_index_mfftw(iadd*2-1,ik) = wf_fft_scnt_mfftw(i,ik)
                    wf_fft_dist_mfftw (iadd*2-1,ik) = i
                    wf_fft_send_mfftw(wf_fft_scnt_mfftw(i,ik),ik,i) = ladd
                    wf_fft_scnt_mfftw(i,ik) = wf_fft_scnt_mfftw(i,ik) + 1
                    wf_fft_index_mfftw(iadd*2  ,ik) = wf_fft_scnt_mfftw(i,ik)
                    wf_fft_dist_mfftw (iadd*2  ,ik) = i
                    wf_fft_send_mfftw(wf_fft_scnt_mfftw(i,ik),ik,i) = ladd
                    exit
                 endif
              enddo B_4
              ista = 2
           endif
           B_1 : do j = ista, iend_g1k(ik)
              iadd = j-ista_g1k(ik)+1
              i1 = igf(nbase(j,ik))
              mz = (i1-1)/(lx*ly)+1
              mm = mod(i1,(lx*ly))
              if (mm==0) mm=lx*ly
              my = (mm-1)/lx+1
              mx = mod(mm,lx)
              if (mx==0) mx = lx
              B_2 : do i = 0, nrank_g-1
                 if (mz>local_n_offsets(i) .and. mz<=local_n_offsets(i)+local_ns(i)) then
                    ladd = mx+lx*(my-1)+lx*ly*(mz-local_n_offsets(i)-1)
                    wf_fft_scnt_mfftw(i,ik) = wf_fft_scnt_mfftw(i,ik) + 1
                    wf_fft_index_mfftw(iadd*2-1,ik) = wf_fft_scnt_mfftw(i,ik)
                    wf_fft_dist_mfftw (iadd*2-1,ik) = i
                    wf_fft_send_mfftw(wf_fft_scnt_mfftw(i,ik),ik,i) = ladd
                    exit
                 endif
              enddo B_2
              i1 = igf(nbase_gamma(j,2))
              mz = (i1-1)/(lx*ly)+1
              mm = mod(i1,(lx*ly))
              if (mm==0) mm=lx*ly
              my = (mm-1)/lx+1
              mx = mod(mm,lx)
              if (mx==0) mx = lx
              B_3 : do i = 0, nrank_g-1
                 if (mz>local_n_offsets(i) .and. mz<=local_n_offsets(i)+local_ns(i)) then
                    ladd = mx+lx*(my-1)+lx*ly*(mz-local_n_offsets(i)-1)
                    wf_fft_scnt_mfftw(i,ik) = wf_fft_scnt_mfftw(i,ik) + 1
                    wf_fft_index_mfftw(iadd*2  ,ik) = wf_fft_scnt_mfftw(i,ik)
                    wf_fft_dist_mfftw (iadd*2  ,ik) = i
                    wf_fft_send_mfftw(wf_fft_scnt_mfftw(i,ik),ik,i) = ladd
                    exit
                 endif
              enddo B_3
           enddo B_1
        else
           B_11 : do j = ista_g1k(ik), iend_g1k(ik)
              iadd = j-ista_g1k(ik)+1
              i1 = igf(nbase(j,ik))
              mz = (i1-1)/(lx*ly)+1
              mm = mod(i1,(lx*ly))
              if (mm==0) mm=lx*ly
              my = (mm-1)/lx+1
              mx = mod(mm,lx)
              if (mx==0) mx = lx
              B_12 : do i = 0, nrank_g-1
                 if (mz>local_n_offsets(i) .and. mz<=local_n_offsets(i)+local_ns(i)) then
                    ladd = mx+lx*(my-1)+lx*ly*(mz-local_n_offsets(i)-1)
                    wf_fft_scnt_mfftw(i,ik) = wf_fft_scnt_mfftw(i,ik) + 1
                    wf_fft_index_mfftw(iadd,ik) = wf_fft_scnt_mfftw(i,ik)
                    wf_fft_dist_mfftw (iadd,ik) = i
                    wf_fft_send_mfftw(wf_fft_scnt_mfftw(i,ik),ik,i) = ladd
                    exit
                 endif
              enddo B_12
           enddo B_11
        endif

     end do
#ifdef USE_NONBLK_COMM
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1354)
     req_r = 0
     req_s = 0
     sta_r = 0
     sta_s = 0
     do i = 0, nrank_g - 1
        call mpi_irecv(wf_fft_recv_mfftw(1,ista_k,i), max_mp_g1k*len*klen, mpi_integer, &
       &               i, itag, mpi_ke_world, req_r(i), ierr)
         if (ierr /= 0) then
            write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_irecv error'
            call flush(nfout)
            call mpi_abort(mpi_comm_world,174,ierr)
         endif
     enddo
     do i = 0, nrank_g - 1
        call mpi_isend(wf_fft_send_mfftw(1,ista_k,i), max_mp_g1k*len*klen, mpi_integer, &
       &               i, itag, mpi_ke_world, req_s(i), ierr)
         if (ierr /= 0) then
            write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_isend error'
            call flush(nfout)
            call mpi_abort(mpi_comm_world,175,ierr)
         endif
     enddo
     call mpi_waitall(nrank_g, req_r, sta_r, ierr)
      if (ierr /= 0) then
         write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_waitall error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world,176,ierr)
    endif
     call mpi_waitall(nrank_g, req_s, sta_s, ierr)
      if (ierr /= 0) then
         write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_waitall error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world,177,ierr)
      endif
                                                  __TIMER_INICOMM_STOP(1354)
#else
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1355)
! === DEBUG by tkato 2012/06/05 ================================================
!    call MPI_ALLTOALL(wf_fft_send(:,ista_k,:), max_mp_g1k*len*klen, mpi_integer, &
!   &                  wf_fft_recv(:,ista_k,:), max_mp_g1k*len*klen, mpi_integer, &
!   &                                          mpi_ke_world, ierr )
     call MPI_ALLTOALL(wf_fft_send_mfftw, max_mp_g1k*len*klen, mpi_integer, &
    &                  wf_fft_recv_mfftw, max_mp_g1k*len*klen, mpi_integer, &
    &                                          mpi_ke_world, ierr )
! ==============================================================================
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_wf_onto_fft_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 178, ierr)
     endif
                                                  __TIMER_INICOMM_STOP(1355)
#endif

     do i = 0, nrank_g - 1
        do ik = ista_k, iend_k
           do j = 1, mp_g1k(ik)*len
              if (wf_fft_recv_mfftw(j,ik,i) == 0) then
                 exit
              end if
              wf_fft_rcnt_mfftw(i,ik) = wf_fft_rcnt_mfftw(i,ik) + 1
           enddo
        enddo
     end do

     do ik = ista_k, iend_k
        wf_fft_maxsend_mfftw(ik) = maxval(wf_fft_scnt_mfftw(0:nrank_g-1,ik))
        wf_fft_maxrecv_mfftw(ik) = maxval(wf_fft_rcnt_mfftw(0:nrank_g-1,ik))
     end do

     allocate(work(maxval(wf_fft_maxrecv_mfftw(:)),ista_k:iend_k,0:nrank_g-1))
     do i = 0, nrank_g - 1
        do ik = ista_k, iend_k
           do j = 1, wf_fft_rcnt_mfftw(i,ik)
              work(j,ik,i) = wf_fft_recv_mfftw(j,ik,i)
           end do
        end do
     end do
     deallocate(wf_fft_recv_mfftw)
     allocate(wf_fft_recv_mfftw(maxval(wf_fft_maxrecv_mfftw(:)),ista_k:iend_k,0:nrank_g-1))
     do i = 0, nrank_g - 1
        do ik = ista_k, iend_k
           do j = 1, wf_fft_rcnt_mfftw(i,ik)
              wf_fft_recv_mfftw(j,ik,i) = work(j,ik,i)
           end do
        end do
     end do
     deallocate(work)

     deallocate(wf_fft_send_mfftw)
     deallocate(local_ns)
     deallocate(local_n_offsets)
                                                  __TIMER_SUB_STOP(1243)
  end subroutine m_Parallel_wf_onto_fft_mpifftw

!===============================================================================
  subroutine m_Parallel_wf_onto_fft_dealloc_mpifftw
     deallocate(wf_fft_rcnt_mfftw)
     deallocate(wf_fft_scnt_mfftw)
     deallocate(wf_fft_recv_mfftw)
!    deallocate(wf_fft_send_mfftw)
     deallocate(wf_fft_index_mfftw)
     deallocate(wf_fft_dist_mfftw)
     deallocate(wf_fft_maxrecv_mfftw)
     deallocate(wf_fft_maxsend_mfftw)
  end subroutine m_Parallel_wf_onto_fft_dealloc_mpifftw

!===============================================================================
  subroutine m_Parallel_fft_onto_wf_mpifftw(nfout,fft_box_size_WF,igf,nbase,kg,kv3,nfft,kk)
    integer, intent(in)  :: nfout, kg, kv3, nfft, kk
    integer, intent(in)  :: igf(kg)
    integer, intent(in)  :: nbase(kg1_ext,kv3)
    integer, intent(in)  :: fft_box_size_WF(1:3,0:1)
    integer, dimension(0:nrank_g-1)                       ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)       ::sta_r, sta_s

    integer, allocatable,dimension(:,:,:) :: fftigf
    integer, allocatable, dimension(:,:,:) :: work
    integer :: i1, lrank, i, j, k, isrsize,fft_l_size, klen, ik, i2
    integer, parameter :: itag = 10
    integer(C_INTPTR_T) :: alloc_local, lx, ly, lz, local_n, local_n_offset
    integer :: lsize, icount
                                                  __TIMER_SUB_START(1244)
!   lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
    lx = fft_box_size_WF(1,0)
    ly = fft_box_size_WF(2,0)
    lz = fft_box_size_WF(3,0)

    alloc_local = fftw_mpi_local_size_3d(lz,ly,lx,mpi_ke_world,local_n,local_n_offset)
    lsize = local_n
    call mpi_allreduce(mpi_in_place,lsize,1,mpi_integer,mpi_max,mpi_ke_world,ierr)

    fft_l_size  = lx*ly*lsize
    isrsize     = lx*ly*lsize

    allocate(fft_wf_scnt_mfftw(0:nrank_g-1,ista_k:iend_k), stat=ierr)
    allocate(fft_wf_rcnt_mfftw(0:nrank_g-1,ista_k:iend_k), stat=ierr)
    allocate(fft_wf_send_mfftw(isrsize,ista_k:iend_k,0:nrank_g-1), stat=ierr)
    allocate(fft_wf_recv_mfftw(isrsize,ista_k:iend_k,0:nrank_g-1), stat=ierr)
    allocate(fft_wf_dist_mfftw(fft_l_size,ista_k:iend_k), stat=ierr)
    allocate(fft_wf_index_mfftw(fft_l_size,ista_k:iend_k), stat=ierr)
    allocate(fft_wf_maxsend_mfftw(ista_k:iend_k), stat=ierr)
    allocate(fft_wf_maxrecv_mfftw(ista_k:iend_k), stat=ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fft_onto_wf_3D :  Not allocate '
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 204, ierr)
     endif
    fft_wf_scnt_mfftw(:,:) = 0
    fft_wf_rcnt_mfftw(:,:) = 0
    fft_wf_send_mfftw(:,:,:) = 0
    fft_wf_recv_mfftw(:,:,:) = 0
    fft_wf_dist_mfftw(:,:) = -1
    fft_wf_index_mfftw(:,:) = 0

    klen = iend_k - ista_k + 1

    allocate(fftigf(nfft,2,ista_k:iend_k), stat=ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fft_onto_wf_3D :  Not allocate '
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 206, ierr)
     endif

    do ik = ista_k, iend_k

       fftigf(:,1,ik) = -1

       do i = 0, nrank_g - 1
          do j = nis_g1k(i,ik), nie_g1k(i,ik)
             i1 = igf(nbase(j,ik))
             if (i1 > nfft) cycle
             fftigf(i1,1,ik) = i
             fftigf(i1,2,ik) = j - nis_g1k(i,ik) + 1
          enddo
       enddo
       icount = 0
       do k=local_n_offset+1,local_n_offset+local_n
       do j=1,ly
       do i=1,lx
          i1 = i + (j-1)*lx + (k-1)*lx*ly
          icount = icount+1
          if (fftigf(i1,1,ik) < 0) cycle
          i2 = i+(j-1)*lx+(k-local_n_offset-1)*lx*ly
          lrank = fftigf(i1,1,ik)
          fft_wf_scnt_mfftw(lrank,ik) = fft_wf_scnt_mfftw(lrank,ik) + 1
          fft_wf_send_mfftw(fft_wf_scnt_mfftw(lrank,ik),ik,lrank) = fftigf(i1,2,ik)
          fft_wf_dist_mfftw(icount,ik) = lrank
          fft_wf_index_mfftw(icount,ik) = fft_wf_scnt_mfftw(lrank,ik)
       end do
       end do
       end do
    end do
    deallocate(fftigf)

#ifdef USE_NONBLK_COMM
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1356)
    lrank = mod(myrank_g,nrank_g)
    do i = 0, nrank_g - 1
       lrank = lrank + 1
       if (lrank > (nrank_g - 1)) lrank = 0
       call mpi_irecv(fft_wf_recv_mfftw(1,ista_k,lrank), isrsize*klen, mpi_integer, &
      &               lrank, itag, mpi_ke_world, req_r(lrank), ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_FFT_Direct_3D :  mpi_irecv error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 207, ierr)
        endif
    enddo

    lrank = mod((myrank_g+1),nrank_g)
    do i = 0, nrank_g - 1
       lrank = lrank + 1
       if (lrank > (nrank_g - 1)) lrank = 0
       call mpi_isend(fft_wf_send_mfftw(1,ista_k,lrank), isrsize*klen, mpi_integer, &
      &               lrank, itag, mpi_ke_world, req_s(lrank), ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_Parallel_fft_onto_wf_3D :  mpi_isend error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 208, ierr)
        endif
    enddo

    call mpi_waitall(nrank_g, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fft_onto_wf_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 209, ierr)
     endif
    call mpi_waitall(nrank_g, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fft_onto_wf_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 210, ierr)
     endif
                                                  __TIMER_INICOMM_STOP(1356)
#else
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1357)
! === DEBUG by tkato 2012/06/05 ================================================
!    call MPI_ALLTOALL(fft_wf_send(:,ista_k,:), isrsize*klen, mpi_integer, &
!   &                  fft_wf_recv(:,ista_k,:), isrsize*klen, mpi_integer, &
!   &                                             mpi_ke_world, ierr )
     call MPI_ALLTOALL(fft_wf_send_mfftw, isrsize*klen, mpi_integer, &
    &                  fft_wf_recv_mfftw, isrsize*klen, mpi_integer, &
    &                                             mpi_ke_world, ierr )
! ==============================================================================
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fft_onto_wf_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 211, ierr)
     endif
                                                  __TIMER_INICOMM_STOP(1357)
#endif

    do i = 0, nrank_g - 1
       do ik = ista_k, iend_k
          do j = 1, isrsize
             if (fft_wf_recv_mfftw(j,ik,i) == 0) then
                exit
             end if
             fft_wf_rcnt_mfftw(i,ik) = fft_wf_rcnt_mfftw(i,ik) + 1
          enddo
       enddo
    end do

    do ik = ista_k, iend_k
       fft_wf_maxsend_mfftw(ik) = maxval(fft_wf_scnt_mfftw(:,ik))
       fft_wf_maxrecv_mfftw(ik) = maxval(fft_wf_rcnt_mfftw(:,ik))
    end do

     allocate(work(maxval(fft_wf_maxrecv_mfftw(:)),ista_k:iend_k,0:nrank_g-1))
     do i = 0, nrank_g - 1
        do ik = ista_k, iend_k
           do j = 1, fft_wf_rcnt_mfftw(i,ik)
              work(j,ik,i) = fft_wf_recv_mfftw(j,ik,i)
           end do
        end do
     end do
     deallocate(fft_wf_recv_mfftw)
     allocate(fft_wf_recv_mfftw(maxval(fft_wf_maxrecv_mfftw(:)),ista_k:iend_k,0:nrank_g-1))
     do i = 0, nrank_g - 1
        do ik = ista_k, iend_k
           do j = 1, fft_wf_rcnt_mfftw(i,ik)
              fft_wf_recv_mfftw(j,ik,i) = work(j,ik,i)
           end do
        end do
     end do
     deallocate(work)
     deallocate(fft_wf_send_mfftw)

                                                  __TIMER_SUB_STOP(1244)
  end subroutine m_Parallel_fft_onto_wf_mpifftw

  subroutine m_Parallel_fft_onto_wf_dealloc_mpifftw
     if(allocated(fft_wf_rcnt_mfftw)) deallocate(fft_wf_rcnt_mfftw)
     if(allocated(fft_wf_scnt_mfftw)) deallocate(fft_wf_scnt_mfftw)
     if(allocated(fft_wf_recv_mfftw)) deallocate(fft_wf_recv_mfftw)
!    deallocate(fft_wf_send)
     if(allocated(fft_wf_index_mfftw)) deallocate(fft_wf_index_mfftw)
     if(allocated(fft_wf_dist_mfftw)) deallocate(fft_wf_dist_mfftw)
     if(allocated(fft_wf_maxsend_mfftw)) deallocate(fft_wf_maxsend)
     if(allocated(fft_wf_maxrecv_mfftw)) deallocate(fft_wf_maxrecv_mfftw)
  end subroutine m_Parallel_fft_onto_wf_dealloc_mpifftw

!===============================================================================

#endif


  subroutine m_Parallel_fft_onto_chgq_3D(nfout,fft_box_size_WF,igf,kg,nfft)

    integer, intent(in)  :: nfout, kg, nfft
    integer, intent(in)  :: igf(kg)
    integer, intent(in)  :: fft_box_size_WF(1:3,0:1)
    integer, dimension(0:nrank_g-1)                       ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)       ::sta_r, sta_s

    integer, allocatable, dimension(:,:) :: work
    integer,allocatable,dimension(:,:) :: fftigf
    integer :: i1, lrank, i, j, k, is, ie, lsize, isrsize, fft_l_size
    integer, parameter :: itag = 10

                                                  __TIMER_SUB_START(1245)
!fj.2012s
!   lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
    lsize = maxval(nel_fft_x(:))
!fj.2012e
    isrsize = min(lsize,mp_kngp_gw)
    fft_l_size  = nel_fft_x(myrank_g)
!   nfft = fft_box_size_WF(1,0)*fft_box_size_WF(2,0)*fft_box_size_WF(3,0)

    allocate(fft_chgq_scnt(0:nrank_g-1), stat=ierr)
    allocate(fft_chgq_rcnt(0:nrank_g-1), stat=ierr)
    allocate(fft_chgq_send(isrsize,0:nrank_g-1), stat=ierr)
    allocate(fft_chgq_recv(isrsize,0:nrank_g-1), stat=ierr)
    allocate(fft_chgq_dist(fft_l_size), stat=ierr)
    allocate(fft_chgq_index(fft_l_size), stat=ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fft_onto_chgq_3D :  Not allocate '
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 250, ierr)
     endif
    fft_chgq_recv(:,:) = 0

    allocate(fftigf(nfft,2), stat=ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fft_onto_chgq_3D :  Not allocate '
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 256, ierr)
     endif
    fftigf(:,1) = -1

#ifdef USE_NONBLK_COMM
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1358)
    lrank = myrank_g
    do i = 0, nrank_g - 1
       lrank = lrank + 1
       if (lrank > (nrank_g-1)) lrank = 0
       call mpi_irecv(fft_chgq_recv(1,lrank), isrsize, mpi_integer, &
      &               lrank, itag, mpi_ske_world, req_r(i), ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_Parallel_fft_onto_chgq_3D :  mpi_irecv error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 257, ierr)
        endif
    enddo
#endif

    do i = 0, nrank_g - 1
       if (is_kngp_gw(i) <= kg) then
          is = is_kngp_gw(i)
          ie = ie_kngp_gw(i)
          if (ie > kg) ie = kg
          do j = is, ie
             i1 = igf(j)
             if (i1 > nfft) cycle
             fftigf(i1,1) = i
             fftigf(i1,2) = j - is + 1
          end do
       end if
    enddo
    fft_chgq_scnt(:) = 0
    fft_chgq_rcnt(:) = 0
    fft_chgq_send(:,:) = 0
    if (fft_l_size /= 0) then
       fft_chgq_dist(:) = -1
       fft_chgq_index(:) = 0
    end if

    do k = 1, nel_fft_x(myrank_g)
       i1 = mp_fft_x(k)
       if (fftigf(i1,1) < 0) cycle
       lrank = fftigf(i1,1)
       fft_chgq_scnt(lrank) = fft_chgq_scnt(lrank) + 1
       fft_chgq_send(fft_chgq_scnt(lrank),lrank) = fftigf(i1,2)
       fft_chgq_dist(k) = lrank
       fft_chgq_index(k) = fft_chgq_scnt(lrank)
    end do
    deallocate(fftigf)

#ifdef USE_NONBLK_COMM
!   do i = 0, nrank_g - 1
!      call mpi_irecv(fft_chgq_recv(1,i), isrsize, mpi_integer, &
!     &               i, itag, mpi_ke_world, req_r(i), ierr)
!       if (ierr /= 0) then
!          write(nfout,*)' m_Parallel_fft_onto_chgq_3D :  mpi_irecv error'
!          call flush(nfout)
!          call mpi_abort(mpi_comm_world, 257, ierr)
!       endif
!   enddo

!   do i = 0, nrank_g - 1
!      call mpi_isend(fft_chgq_send(1,i), isrsize, mpi_integer, &
!     &               i, itag, mpi_ke_world, req_s(i), ierr)
!       if (ierr /= 0) then
!          write(nfout,*)' m_Parallel_fft_onto_chgq_3D :  mpi_isend error'
!          call flush(nfout)
!          call mpi_abort(mpi_comm_world, 258, ierr)
!       endif
!   enddo

    lrank = myrank_g
    do i = 0, nrank_g - 1
       lrank = lrank - 1
       if (lrank < 0) lrank = nrank_g - 1
       call mpi_isend(fft_chgq_send(1,lrank), isrsize, mpi_integer, &
      &               lrank, itag, mpi_ske_world, req_s(i), ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_Parallel_fft_onto_chgq_3D :  mpi_isend error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 258, ierr)
        endif
    enddo

    call mpi_waitall(nrank_g, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fft_onto_chgq_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 259, ierr)
     endif
    call mpi_waitall(nrank_g, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fft_onto_chgq_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 260, ierr)
     endif
                                                  __TIMER_INICOMM_STOP(1358)
#else
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1359)
     call MPI_ALLTOALL(fft_chgq_send, isrsize, mpi_integer, &
    &                  fft_chgq_recv, isrsize, mpi_integer, &
    &                                          mpi_ske_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fft_onto_chgq_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 261, ierr)
     endif
                                                  __TIMER_INICOMM_STOP(1359)
#endif

    do i = 0, nrank_g - 1
       do j = 1, isrsize
          if (fft_chgq_recv(j,i) == 0) then
             exit
          end if
          fft_chgq_rcnt(i) = fft_chgq_rcnt(i) + 1
       enddo
    end do
    fft_chgq_maxsend = maxval(fft_chgq_scnt(:))
    fft_chgq_maxrecv = maxval(fft_chgq_rcnt(:))

    allocate(work(1:fft_chgq_maxrecv,0:nrank_g-1))
    do i = 0, nrank_g - 1
       do j = 1, fft_chgq_rcnt(i)
          work(j,i) = fft_chgq_recv(j,i)
       end do
    end do
    deallocate(fft_chgq_recv)
    allocate(fft_chgq_recv(fft_chgq_maxrecv,0:nrank_g-1))
    do i = 0, nrank_g - 1
       do j = 1, fft_chgq_rcnt(i)
          fft_chgq_recv(j,i) = work(j,i)
       end do
    end do
    deallocate(work)

     deallocate(fft_chgq_send)
                                                  __TIMER_SUB_STOP(1245)
  end subroutine m_Parallel_fft_onto_chgq_3D
!===============================================================================

  subroutine m_Parallel_fft_onto_chgq_dealloc()
    deallocate(fft_chgq_scnt)
    deallocate(fft_chgq_rcnt)
    deallocate(fft_chgq_dist)
    deallocate(fft_chgq_index)
    deallocate(fft_chgq_recv)
  end subroutine m_Parallel_fft_onto_chgq_dealloc

#ifdef MPI_FFTW
  subroutine m_Parallel_fft_onto_chgq_mpifftw(nfout,fft_box_size_WF,igf,kg,nfft)

    integer, intent(in)  :: nfout, kg, nfft
    integer, intent(in)  :: igf(kg)
    integer, intent(in)  :: fft_box_size_WF(1:3,0:1)
    integer, dimension(0:nrank_g-1)                       ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)       ::sta_r, sta_s

    integer, allocatable, dimension(:,:) :: work
    integer,allocatable,dimension(:,:) :: fftigf
    integer :: i1, i2, lrank, i, j, k, is, ie, lsize, isrsize, fft_l_size, icount
    integer, parameter :: itag = 10
    integer(C_INTPTR_T) :: alloc_local, lx, ly, lz, local_n, local_n_offset

    lx = fft_box_size_WF(1,0)
    ly = fft_box_size_WF(2,0)
    lz = fft_box_size_WF(3,0)

    alloc_local = fftw_mpi_local_size_3d(lz,ly,lx,mpi_ke_world,local_n,local_n_offset)
    lsize = local_n
    call mpi_allreduce(mpi_in_place,lsize,1,mpi_integer,mpi_max,mpi_ke_world,ierr)
    fft_l_size  = lx*ly*lsize
    isrsize     = lx*ly*lsize

    allocate(fft_chgq_scnt_mpifftw(0:nrank_g-1), stat=ierr)
    allocate(fft_chgq_rcnt_mpifftw(0:nrank_g-1), stat=ierr)
    allocate(fft_chgq_send_mpifftw(isrsize,0:nrank_g-1), stat=ierr)
    allocate(fft_chgq_recv_mpifftw(isrsize,0:nrank_g-1), stat=ierr)
    allocate(fft_chgq_dist_mpifftw(fft_l_size), stat=ierr)
    allocate(fft_chgq_index_mpifftw(fft_l_size), stat=ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fft_onto_chgq_3D :  Not allocate '
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 250, ierr)
     endif
    fft_chgq_recv_mpifftw(:,:) = 0

    allocate(fftigf(nfft,2), stat=ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fft_onto_chgq_3D :  Not allocate '
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 256, ierr)
     endif
    fftigf(:,1) = -1

    lrank = myrank_g
    do i = 0, nrank_g - 1
       lrank = lrank + 1
       if (lrank > (nrank_g-1)) lrank = 0
       call mpi_irecv(fft_chgq_recv_mpifftw(1,lrank), isrsize, mpi_integer, &
      &               lrank, itag, mpi_ke_world, req_r(i), ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_Parallel_fft_onto_chgq_3D :  mpi_irecv error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 257, ierr)
        endif
    enddo

    do i = 0, nrank_g - 1
       if (is_kngp_gw(i) <= kg) then
          is = is_kngp_gw(i)
          ie = ie_kngp_gw(i)
          if (ie > kg) ie = kg
          do j = is, ie
             i1 = igf(j)
             if (i1 > nfft) cycle
             fftigf(i1,1) = i
             fftigf(i1,2) = j - is + 1
          end do
       end if
    enddo
    fft_chgq_scnt_mpifftw(:) = 0
    fft_chgq_rcnt_mpifftw(:) = 0
    fft_chgq_send_mpifftw(:,:) = 0
    if (fft_l_size /= 0) then
       fft_chgq_dist_mpifftw(:) = -1
       fft_chgq_index_mpifftw(:) = 0
    end if
    icount = 0
    do k=local_n_offset+1,local_n_offset+local_n
    do j=1,ly
    do i=1,lx
       i1 = i + (j-1)*lx + (k-1)*lx*ly
       icount = icount+1
       if (fftigf(i1,1) < 0) cycle
       i2 = i+(j-1)*lx+(k-local_n_offset-1)*lx*ly
       lrank = fftigf(i1,1)
       fft_chgq_scnt_mpifftw(lrank) = fft_chgq_scnt_mpifftw(lrank) + 1
       fft_chgq_send_mpifftw(fft_chgq_scnt_mpifftw(lrank),lrank) = fftigf(i1,2)
       fft_chgq_dist_mpifftw(icount) = lrank
       fft_chgq_index_mpifftw(icount) = fft_chgq_scnt_mpifftw(lrank)
    enddo
    enddo
    end do
    deallocate(fftigf)

    lrank = myrank_g
    do i = 0, nrank_g - 1
       lrank = lrank - 1
       if (lrank < 0) lrank = nrank_g - 1
       call mpi_isend(fft_chgq_send_mpifftw(1,lrank), isrsize, mpi_integer, &
      &               lrank, itag, mpi_ke_world, req_s(i), ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_Parallel_fft_onto_chgq_3D :  mpi_isend error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 258, ierr)
        endif
    enddo

    call mpi_waitall(nrank_g, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fft_onto_chgq_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 259, ierr)
     endif
    call mpi_waitall(nrank_g, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fft_onto_chgq_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 260, ierr)
     endif

    do i = 0, nrank_g - 1
       do j = 1, isrsize
          if (fft_chgq_recv_mpifftw(j,i) == 0) then
             exit
          end if
          fft_chgq_rcnt_mpifftw(i) = fft_chgq_rcnt_mpifftw(i) + 1
       enddo
    end do
    fft_chgq_maxsend_mpifftw = maxval(fft_chgq_scnt_mpifftw(:))
    fft_chgq_maxrecv_mpifftw = maxval(fft_chgq_rcnt_mpifftw(:))

    allocate(work(1:fft_chgq_maxrecv_mpifftw,0:nrank_g-1))
    do i = 0, nrank_g - 1
       do j = 1, fft_chgq_rcnt_mpifftw(i)
          work(j,i) = fft_chgq_recv_mpifftw(j,i)
       end do
    end do
    deallocate(fft_chgq_recv_mpifftw)
    allocate(fft_chgq_recv_mpifftw(fft_chgq_maxrecv_mpifftw,0:nrank_g-1))
    do i = 0, nrank_g - 1
       do j = 1, fft_chgq_rcnt_mpifftw(i)
          fft_chgq_recv_mpifftw(j,i) = work(j,i)
       end do
    end do
    deallocate(work)

     deallocate(fft_chgq_send_mpifftw)
                                                  __TIMER_SUB_STOP(1245)
  end subroutine m_Parallel_fft_onto_chgq_mpifftw
!===============================================================================

  subroutine m_Parallel_fft_onto_chgq_mpifftw_dealloc()
    deallocate(fft_chgq_scnt_mpifftw)
    deallocate(fft_chgq_rcnt_mpifftw)
    deallocate(fft_chgq_dist_mpifftw)
    deallocate(fft_chgq_index_mpifftw)
    deallocate(fft_chgq_recv_mpifftw)
  end subroutine m_Parallel_fft_onto_chgq_mpifftw_dealloc
#endif

#ifdef FFT_3D_DIVISION_CD
!===============================================================================
  subroutine m_Parallel_chgq_onto_fftcd_3D(nfout,kgp,fft_box_size_CD_3D,igfp_l)

    integer, intent(in) :: nfout, kgp
    integer, dimension(3,0:1),intent(in) :: fft_box_size_CD_3D
    integer, dimension(ista_kngp:iend_kngp),intent(in) :: igfp_l

    integer, allocatable, dimension(:,:,:) :: xyz
    integer, allocatable, dimension(:)   :: wk_recvdsp
    integer, allocatable, dimension(:,:) :: work

    integer, dimension(0:npes-1)                       :: req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:npes-1)       :: sta_r, sta_s

    integer :: max_fft_x, lx, ly, lz, mx, my, mz, mm, i, j, i1
    integer :: iadd, ladd, lrank, iend
    integer, parameter :: itag = 18
    integer :: kx1p, kx2p, kx3p
                                                  __TIMER_SUB_START(1246)

    if(allocated(igfp_full)) deallocate(igfp_full)
    allocate(igfp_full(1:kgp), stat=ierr)

    allocate(wk_recvdsp(0:nrank_g-1))
    wk_recvdsp(0) = 0
    do i = 1, nrank_g - 1
       wk_recvdsp(i) = wk_recvdsp(i-1) + nel_kngp(i-1)
    end do
    call mpi_allgatherv(igfp_l, np_kngp, mpi_integer, igfp_full, &
   &                    nel_kngp, wk_recvdsp, mpi_integer, mpi_ke_world, ierr)
    deallocate(wk_recvdsp)

     max_fft_x = maxval(nel_fftcd_x(:))
     allocate(xyz(2,3,0:nrank_g-1))
#ifdef USE_NONBLK_COMM
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1360)
     do i = 0, nrank_g - 1
        call mpi_irecv(xyz(1,1,i), 6, mpi_integer, &
       &               i, itag, mpi_ke_world, req_r(i), ierr)
         if (ierr /= 0) then
            write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_irecv error'
            call flush(nfout)
            call mpi_abort(mpi_comm_world,154,ierr)
         endif
     enddo
     do i = 0, nrank_g - 1
        call mpi_isend(xyz_fftcd_x, 6, mpi_integer, &
       &               i, itag, mpi_ke_world, req_s(i), ierr)
         if (ierr /= 0) then
            write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_isend error'
            call flush(nfout)
            call mpi_abort(mpi_comm_world,155,ierr)
         endif
     enddo
     call mpi_waitall(nrank_g, req_r, sta_r, ierr)
      if (ierr /= 0) then
         write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_waitall error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world,156,ierr)
      endif
     call mpi_waitall(nrank_g, req_s, sta_s, ierr)
      if (ierr /= 0) then
         write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_waitall error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world,157,ierr)
      endif
                                                  __TIMER_INICOMM_STOP(1360)
#else
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1701)
     integer, allocatable, dimension(:,:) :: rbuf
     allocate(rbuf(6,0:nrank_g-1))
     call MPI_ALLGATHER(xyz_fftcd_x, 6, mpi_integer, &
    &                   rbuf,6, mpi_integer, &
    &                                mpi_ke_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_allgather error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 174, ierr)
     endif
     do i = 1,2
       do j = 1,3
         xyz(i,j,:)=rbuf(i+(j-1)*2,:)
       enddo
     enddo
     deallocate(rbuf)
                                                  __TIMER_INICOMM_STOP(1701)
#endif

     lx = fft_box_size_CD_3D(1,0)
     ly = fft_box_size_CD_3D(2,0)
     lz = fft_box_size_CD_3D(3,0)

     kx1p = fftcd_X_x_nel
     kx2p = fftcd_X_y_nel
     kx3p = fftcd_X_z_nel

     allocate(chgq_fftcd_scnt(0:nrank_g-1))
     allocate(chgq_fftcd_rcnt(0:nrank_g-1))
     chgq_fftcd_scnt(:) = 0
     chgq_fftcd_rcnt(:) = 0

     allocate(chgq_fftcd_index(np_kngp_gw))
     allocate(chgq_fftcd_dist(np_kngp_gw))
     allocate(chgq_fftcd_send(mp_kngp_gw,0:nrank_g-1))
     allocate(chgq_fftcd_recv(mp_kngp_gw,0:nrank_g-1))
     chgq_fftcd_index(:) = 0
     chgq_fftcd_dist(:) = -1
     chgq_fftcd_send(:,:) = 0
     chgq_fftcd_recv(:,:) = 0
     iend = iend_kngp_gw
     if (iend > kgp) iend = kgp
     do j = ista_kngp_gw, iend
        iadd = j-ista_kngp_gw+1
        i1 = igfp_l(j)
        mz = (i1-1)/(lx*ly)+1
        mm = mod(i1,(lx*ly))
        if (mm==0) mm=lx*ly
        my = (mm-1)/lx+1
        mx = mod(mm,lx)
        if (mx==0) mx = lx
        do i = 0, nrank_g-1
           if ((xyz(1,1,i)<=mx).and.(mx<=xyz(2,1,i)).and.    &
          &    (xyz(1,2,i)<=my).and.(my<=xyz(2,2,i)).and.    &
          &    (xyz(1,3,i)<=mz).and.(mz<=xyz(2,3,i))) then
              ladd = mx-xyz(1,1,i)+1+kx1p*(my-xyz(1,2,i))+kx1p*kx2p*(mz-xyz(1,3,i))
              chgq_fftcd_scnt(i) = chgq_fftcd_scnt(i) + 1
              chgq_fftcd_index(iadd) = chgq_fftcd_scnt(i)
              chgq_fftcd_dist (iadd) = i
              chgq_fftcd_send(chgq_fftcd_scnt(i),i) = ladd
              exit
           endif
        enddo
     enddo
     deallocate(xyz)

#ifdef USE_NONBLK_COMM
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1702)
     do lrank = 0, nrank_g-1
        call mpi_irecv(chgq_fftcd_recv(1,lrank), mp_kngp_gw, mpi_integer, &
       &               lrank, itag, mpi_ke_world, req_r(lrank), ierr)
         if (ierr /= 0) then
            write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_irecv error'
            call flush(nfout)
            call mpi_abort(mpi_comm_world,158,ierr)
         endif
     enddo
     do lrank = 0, nrank_g - 1
        call mpi_isend(chgq_fftcd_send(1,lrank), mp_kngp_gw, mpi_integer, &
       &               lrank, itag, mpi_ke_world, req_s(lrank), ierr)
         if (ierr /= 0) then
            write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_isend error'
            call flush(nfout)
            call mpi_abort(mpi_comm_world,159,ierr)
         endif
     enddo
     call mpi_waitall(nrank_g, req_r, sta_r, ierr)
      if (ierr /= 0) then
         write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_waitall error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world,160,ierr)
      endif
     call mpi_waitall(nrank_g, req_s, sta_s, ierr)
      if (ierr /= 0) then
         write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_waitall error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world,161,ierr)
      endif
                                                  __TIMER_INICOMM_STOP(1702)
#else
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1703)
     call MPI_ALLTOALL(chgq_fftcd_send, mp_kngp_gw, mpi_integer, &
    &                  chgq_fftcd_recv, mp_kngp_gw, mpi_integer, &
    &                                            mpi_ke_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 162, ierr)
     endif
                                                  __TIMER_INICOMM_STOP(1703)
#endif

     do i = 0, nrank_g - 1
        do j = 1, nel_kngp_gw(i)
           if (chgq_fftcd_recv(j,i) == 0) then
              exit
           end if
           chgq_fftcd_rcnt(i) = chgq_fftcd_rcnt(i) + 1
        enddo
     end do
     chgq_fftcd_maxsend = maxval(chgq_fftcd_scnt)
     chgq_fftcd_maxrecv = maxval(chgq_fftcd_rcnt)
     allocate(work(1:chgq_fftcd_maxrecv,0:nrank_g-1))
     do i = 0, nrank_g - 1
        do j = 1, chgq_fftcd_rcnt(i)
           work(j,i) = chgq_fftcd_recv(j,i)
        end do
     end do
     deallocate(chgq_fftcd_recv)
     allocate(chgq_fftcd_recv(chgq_fftcd_maxrecv,0:nrank_g-1))
     do i = 0, nrank_g - 1
        do j = 1, chgq_fftcd_rcnt(i)
           chgq_fftcd_recv(j,i) = work(j,i)
        end do
     end do
     deallocate(work)

     deallocate(chgq_fftcd_send)
                                                  __TIMER_SUB_START(1246)
  end subroutine m_Parallel_chgq_onto_fftcd_3D
!===============================================================================
#else
!ifdef FFT_3D_DIVISION_CD
!===============================================================================
  subroutine m_Parallel_chgq_onto_fftcd_3D(nfout,kgp,fft_box_size_CD_3D,igfp_l)

    integer, intent(in) :: nfout, kgp
    integer, dimension(3,0:1),intent(in) :: fft_box_size_CD_3D
    integer, dimension(ista_kngp:iend_kngp),intent(in) :: igfp_l

    integer, allocatable, dimension(:,:,:) :: xyz
    integer, allocatable, dimension(:)   :: wk_recvdsp
    integer, allocatable, dimension(:,:) :: work

    integer, dimension(0:npes-1)                       :: req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:npes-1)       :: sta_r, sta_s

    integer :: max_fft_x, lx, ly, lz, mx, my, mz, mm, i, j, i1
    integer :: iadd, ladd, lrank, iend
    integer, parameter :: itag = 18
! === DEBUG by tkato 2012/06/04 ================================================
#ifndef USE_NONBLK_COMM
     integer, allocatable, dimension(:,:) :: rbuf
#endif
! ==============================================================================

                                                  __TIMER_SUB_START(1246)

    if(allocated(igfp_full)) deallocate(igfp_full)
    allocate(igfp_full(1:kgp), stat=ierr)

    allocate(wk_recvdsp(0:nrank_chg-1))
    wk_recvdsp(0) = 0
    do i = 1, nrank_chg - 1
       wk_recvdsp(i) = wk_recvdsp(i-1) + nel_kngp(i-1)
    end do
    call mpi_allgatherv(igfp_l, np_kngp, mpi_integer, igfp_full, &
   &                    nel_kngp, wk_recvdsp, mpi_integer, mpi_chg_world, ierr)
    deallocate(wk_recvdsp)

     max_fft_x = maxval(nel_fftcd_x(:))
     allocate(xyz(2,3,0:nrank_g-1))
#ifdef USE_NONBLK_COMM
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1360)
     do i = 0, nrank_g - 1
        call mpi_irecv(xyz(1,1,i), 6, mpi_integer, &
       &               i, itag, mpi_ke_world, req_r(i), ierr)
         if (ierr /= 0) then
            write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_irecv error'
            call flush(nfout)
            call mpi_abort(mpi_comm_world,154,ierr)
         endif
     enddo
     do i = 0, nrank_g - 1
        call mpi_isend(xyz_fftcd_x, 6, mpi_integer, &
       &               i, itag, mpi_ke_world, req_s(i), ierr)
         if (ierr /= 0) then
            write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_isend error'
            call flush(nfout)
            call mpi_abort(mpi_comm_world,155,ierr)
         endif
     enddo
     call mpi_waitall(nrank_g, req_r, sta_r, ierr)
      if (ierr /= 0) then
         write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_waitall error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world,156,ierr)
      endif
     call mpi_waitall(nrank_g, req_s, sta_s, ierr)
      if (ierr /= 0) then
         write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_waitall error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world,157,ierr)
      endif
                                                  __TIMER_INICOMM_STOP(1360)
#else
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1701)
! === DEBUG by tkato 2012/06/04 ================================================
!    integer, allocatable, dimension(:,:) :: rbuf
! ==============================================================================
     allocate(rbuf(6,0:nrank_g-1))
     call MPI_ALLGATHER(xyz_fftcd_x, 6, mpi_integer, &
    &                   rbuf,6, mpi_integer, &
    &                                mpi_ke_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_allgather error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 174, ierr)
     endif
     do i = 1,2
       do j = 1,3
         xyz(i,j,:)=rbuf(i+(j-1)*2,:)
       enddo
     enddo
     deallocate(rbuf)
                                                  __TIMER_INICOMM_STOP(1701)
#endif

     lx = fft_box_size_CD_3D(1,0)
     ly = fft_box_size_CD_3D(2,0)
     lz = fft_box_size_CD_3D(3,0)
     allocate(chgq_fftcd_scnt(0:nrank_g-1))
     allocate(chgq_fftcd_rcnt(0:nrank_g-1))
     chgq_fftcd_scnt(:) = 0
     chgq_fftcd_rcnt(:) = 0

     allocate(chgq_fftcd_index(np_kngp_gw))
     allocate(chgq_fftcd_dist(np_kngp_gw))
     allocate(chgq_fftcd_send(mp_kngp_gw,0:nrank_g-1))
     allocate(chgq_fftcd_recv(mp_kngp_gw,0:nrank_g-1))
     chgq_fftcd_index(:) = 0
     chgq_fftcd_dist(:) = -1
     chgq_fftcd_send(:,:) = 0
     chgq_fftcd_recv(:,:) = 0
     iend = iend_kngp_gw
     if (iend > kgp) iend = kgp
     do j = ista_kngp_gw, iend
        iadd = j-ista_kngp_gw+1
!        i1 = igfp_l(j)
        i1 = igfp_full(j)
        mz = (i1-1)/(lx*ly)+1
        mm = mod(i1,(lx*ly))
        if (mm==0) mm=lx*ly
        my = (mm-1)/lx+1
        mx = mod(mm,lx)
        if (mx==0) mx = lx
        do i = 0, nrank_g-1
           if ((xyz(1,2,i)<=my).and.(my<=xyz(2,2,i)).and.    &
          &    (xyz(1,3,i)<=mz).and.(mz<=xyz(2,3,i))) then
              ladd = mx+lx*(my-xyz(1,2,i))+lx*(xyz(2,2,i)-xyz(1,2,i)+1)*(mz-xyz(1,3,i))
              chgq_fftcd_scnt(i) = chgq_fftcd_scnt(i) + 1
              chgq_fftcd_index(iadd) = chgq_fftcd_scnt(i)
              chgq_fftcd_dist (iadd) = i
              chgq_fftcd_send(chgq_fftcd_scnt(i),i) = ladd
              exit
           endif
        enddo
     enddo
     deallocate(xyz)

#ifdef USE_NONBLK_COMM
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1702)
     do lrank = 0, nrank_g-1
        call mpi_irecv(chgq_fftcd_recv(1,lrank), mp_kngp_gw, mpi_integer, &
       &               lrank, itag, mpi_ke_world, req_r(lrank), ierr)
         if (ierr /= 0) then
            write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_irecv error'
            call flush(nfout)
            call mpi_abort(mpi_comm_world,158,ierr)
         endif
     enddo
     do lrank = 0, nrank_g - 1
        call mpi_isend(chgq_fftcd_send(1,lrank), mp_kngp_gw, mpi_integer, &
       &               lrank, itag, mpi_ke_world, req_s(lrank), ierr)
         if (ierr /= 0) then
            write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_isend error'
            call flush(nfout)
            call mpi_abort(mpi_comm_world,159,ierr)
         endif
     enddo
     call mpi_waitall(nrank_g, req_r, sta_r, ierr)
      if (ierr /= 0) then
         write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_waitall error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world,160,ierr)
      endif
     call mpi_waitall(nrank_g, req_s, sta_s, ierr)
      if (ierr /= 0) then
         write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_waitall error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world,161,ierr)
      endif
                                                  __TIMER_INICOMM_STOP(1702)
#else
                                                  __TIMER_INICOMM_START_w_BARRIER(mpi_ke_world,1703)
     call MPI_ALLTOALL(chgq_fftcd_send, mp_kngp_gw, mpi_integer, &
    &                  chgq_fftcd_recv, mp_kngp_gw, mpi_integer, &
    &                                            mpi_ke_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_chgq_onto_fftcd_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 162, ierr)
     endif
                                                  __TIMER_INICOMM_STOP(1703)
#endif

     do i = 0, nrank_g - 1
        do j = 1, nel_kngp_gw(i)
           if (chgq_fftcd_recv(j,i) == 0) then
              exit
           end if
           chgq_fftcd_rcnt(i) = chgq_fftcd_rcnt(i) + 1
        enddo
     end do
     chgq_fftcd_maxsend = maxval(chgq_fftcd_scnt)
     chgq_fftcd_maxrecv = maxval(chgq_fftcd_rcnt)
     allocate(work(1:chgq_fftcd_maxrecv,0:nrank_g-1))
     do i = 0, nrank_g - 1
        do j = 1, chgq_fftcd_rcnt(i)
           work(j,i) = chgq_fftcd_recv(j,i)
        end do
     end do
     deallocate(chgq_fftcd_recv)
     allocate(chgq_fftcd_recv(chgq_fftcd_maxrecv,0:nrank_g-1))
     do i = 0, nrank_g - 1
        do j = 1, chgq_fftcd_rcnt(i)
           chgq_fftcd_recv(j,i) = work(j,i)
        end do
     end do
     deallocate(work)

     deallocate(chgq_fftcd_send)
                                                  __TIMER_SUB_STOP(1246)
  end subroutine m_Parallel_chgq_onto_fftcd_3D
!===============================================================================
#endif
!ifdef FFT_3D_DIVISION_CD
  subroutine m_Parallel_chgq_onto_fftcd_dealloc()
      deallocate(chgq_fftcd_scnt)
      deallocate(chgq_fftcd_rcnt)
      deallocate(chgq_fftcd_index)
      deallocate(chgq_fftcd_dist)
      deallocate(chgq_fftcd_recv)
  end subroutine m_Parallel_chgq_onto_fftcd_dealloc

!===============================================================================
  subroutine m_Parallel_fftcd_onto_chgq_3D(nfout,kgp,fft_box_size_CD_3D,igfp_l)

    integer,intent(in) :: nfout,kgp
    integer, dimension(3,0:1),intent(in) :: fft_box_size_CD_3D
    integer, dimension(ista_kngp:iend_kngp),intent(in) :: igfp_l
    integer, allocatable, dimension(:,:) :: fftigf

#ifdef CD_FFT_ALL
    integer, dimension(0:npes-1)                       :: req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:npes-1)       :: sta_r, sta_s
#else
    integer, dimension(0:nrank_g-1)                 :: req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1) :: sta_r, sta_s
#endif
    integer, allocatable, dimension(:,:) :: work

    integer :: i1, lrank, i, j, k, is, ie, jrank_e, jrank_g
    integer :: nfftps, isrsize, fftsize,ltmp, nmrank, myrank, mpicom
    integer, parameter :: itag = 10
                                                  __TIMER_SUB_START(1247)

#ifdef CD_FFT_ALL
    mpicom = MPI_CommGroup
#else
    mpicom = mpi_ke_world
#endif
    call mpi_comm_size(mpicom, nmrank, ierr)
    call mpi_comm_rank(mpicom, myrank, ierr)

    nfftps = fft_box_size_CD_3D(1,0) * fft_box_size_CD_3D(2,0) * fft_box_size_CD_3D(3,0)
#ifdef FFT_3D_DIVISION_CD
    ltmp = maxval(nel_fftcd_x(:))
#else
    ltmp = max(maxval(nel_fftcd_x(:)),maxval(nel_fftcd_y(:)),maxval(nel_fftcd_z(:)))
#endif
    isrsize = min(ltmp,mp_kngp_gw)
    fftsize = nel_fftcd_x(myrank)

    allocate(fftigf(nfftps,2), stat=ierr)
    if (ierr /= 0) then
       write(nfout,*)' m_Parallel_fftcd_onto_chgq_3D :  Not allocate '
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 256, ierr)
    endif
   fftigf(:,1) = -1
   allocate(fftcd_chgq_scnt(0:nmrank-1), stat=ierr)
   allocate(fftcd_chgq_rcnt(0:nmrank-1), stat=ierr)
   allocate(fftcd_chgq_send(isrsize,0:nmrank-1), stat=ierr)
   allocate(fftcd_chgq_recv(isrsize,0:nmrank-1), stat=ierr)
   allocate(fftcd_chgq_dist(fftsize), stat=ierr)
   allocate(fftcd_chgq_index(fftsize), stat=ierr)
    if (ierr /= 0) then
       write(nfout,*)' m_Parallel_fftcd_onto_chgq_3D :  Not allocate '
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 250, ierr)
    endif
   fftcd_chgq_recv(:,:) = 0


   do i = 0, nrank_g - 1
!     if (is_kngp(i) <= kg) then
      if (is_kngp_gw(i) <= kgp) then
         is = is_kngp_gw(i)
         ie = ie_kngp_gw(i)
!        if (ie > kg) ie = kg
         if (ie > kgp) ie = kgp
         do j = is, ie
            i1 = igfp_full(j)
!           if (i1 > nfftp) cycle
            if (i1 > nfftps) cycle
            fftigf(i1,1) = i
            fftigf(i1,2) = j - is + 1
         end do
      end if
   enddo
   fftcd_chgq_scnt(:) = 0
   fftcd_chgq_rcnt(:) = 0
   fftcd_chgq_send(:,:) = 0
   if (fftsize /= 0) then
      fftcd_chgq_dist(:) = -1
      fftcd_chgq_index(:) = 0
   end if

   do k = 1, nel_fftcd_x(myrank)
      i1 = mp_fftcd_x(k)
      if (fftigf(i1,1) < 0) cycle

      lrank = fftigf(i1,1)
      fftcd_chgq_scnt(lrank) = fftcd_chgq_scnt(lrank) + 1
      fftcd_chgq_send(fftcd_chgq_scnt(lrank),lrank) = fftigf(i1,2)
      fftcd_chgq_dist(k) = lrank
      fftcd_chgq_index(k) = fftcd_chgq_scnt(lrank)
   end do
   deallocate(fftigf)

#ifdef USE_NONBLK_COMM
                                                  __TIMER_INICOMM_START_w_BARRIER(mpicom,1704)
   do i = 0, nmrank - 1
      call mpi_irecv(fftcd_chgq_recv(1,i), isrsize, mpi_integer, &
     &               i, itag, mpicom, req_r(i), ierr)
       if (ierr /= 0) then
          write(nfout,*)' m_Parallel_fftcd_onto_chgq_3D :  mpi_irecv error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 257, ierr)
       endif
   enddo

#ifdef CD_FFT_ALL
   jrank_e = 0
   jrank_g = 0
   do i = 0, nmrank - 1
      call mpi_isend(fftcd_chgq_send(1,jrank_g), isrsize, mpi_integer, &
     &               i, itag, mpicom, req_s(i), ierr)
       if (ierr /= 0) then
            write(nfout,*)' m_Parallel_fftcd_onto_chgq_3D :  mpi_isend error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 258, ierr)
       endif

      jrank_e = jrank_e + 1
      if (jrank_e > (nrank_e-1)) then
         jrank_e = 0
         jrank_g = jrank_g + 1
      end if
   end do
#else
   do i = 0, nmrank - 1
      call mpi_isend(fftcd_chgq_send(1,i), isrsize, mpi_integer, &
     &               i, itag, mpicom, req_s(i), ierr)
       if (ierr /= 0) then
            write(nfout,*)' m_Parallel_fftcd_onto_chgq_3D :  mpi_isend error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 258, ierr)
       endif
   enddo
#endif

   call mpi_waitall(nmrank, req_r, sta_r, ierr)
    if (ierr /= 0) then
         write(nfout,*)' m_Parallel_fftcd_onto_chgq_3D :  mpi_waitall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 259, ierr)
    endif
   call mpi_waitall(nmrank, req_s, sta_s, ierr)
    if (ierr /= 0) then
       write(nfout,*)' m_Parallel_fftcd_onto_chgq_3D :  mpi_waitall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 260, ierr)
    endif
                                                  __TIMER_INICOMM_STOP(1704)
#else
                                                  __TIMER_INICOMM_START_w_BARRIER(mpicom,1705)
     call MPI_ALLTOALL(fftcd_chgq_send, isrsize, mpi_integer, &
    &                  fftcd_chgq_recv, isrsize, mpi_integer, &
    &                                            mpicom, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_Parallel_fftcd_onto_chgq_3D :  mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 261, ierr)
     endif
                                                  __TIMER_INICOMM_STOP(1705)
#endif

   do i = 0, nmrank - 1
      do j = 1, isrsize
         if (fftcd_chgq_recv(j,i) == 0) then
            exit
         end if
         fftcd_chgq_rcnt(i) = fftcd_chgq_rcnt(i) + 1
      enddo
   end do
   fftcd_chgq_maxsend = maxval(fftcd_chgq_scnt(:))
   fftcd_chgq_maxrecv = maxval(fftcd_chgq_rcnt(:))

   allocate(work(1:fftcd_chgq_maxrecv,0:nmrank-1))
   do i = 0, nmrank - 1
      do j = 1, fftcd_chgq_rcnt(i)
         work(j,i) = fftcd_chgq_recv(j,i)
      end do
   end do
   deallocate(fftcd_chgq_recv)
   allocate(fftcd_chgq_recv(fftcd_chgq_maxrecv,0:nmrank-1))
   do i = 0, nmrank - 1
      do j = 1, fftcd_chgq_rcnt(i)
         fftcd_chgq_recv(j,i) = work(j,i)
      end do
   end do
   deallocate(work)

   deallocate(fftcd_chgq_send)

                                                  __TIMER_SUB_STOP(1247)
  end subroutine m_Parallel_fftcd_onto_chgq_3D

  subroutine m_Parallel_fftcd_onto_chgq_dealloc()
    deallocate(fftcd_chgq_scnt)
    deallocate(fftcd_chgq_rcnt)
    deallocate(fftcd_chgq_dist)
    deallocate(fftcd_chgq_index)
    deallocate(fftcd_chgq_recv)
  end subroutine m_Parallel_fftcd_onto_chgq_dealloc

!===============================================================================

  subroutine m_Parallel_init_mpi_G_iba2_3D(nfout,ipri,printable,kv3,iba2)
    integer, intent(in) :: nfout,ipri, kv3, iba2(kv3)
    logical, intent(in) :: printable

    logical             :: set_mapping_func
    integer             :: i, ik
    integer, parameter  :: nc = 24
    integer, parameter  :: nc1 = 16

    if(nrank_k > kv3) call phase_error_with_msg(nfout,' nrank_k > kv3 (m_Parallel_init_mpi_elec)',__LINE__,__FILE__)

! (( ista_G_g1k, iend_G_g1k ))
    if(allocated(np_G_g1k)) deallocate(np_G_g1k);  allocate(np_G_g1k(kv3))
    if(allocated(mp_G_g1k)) deallocate(mp_G_g1k);  allocate(mp_G_g1k(kv3))
    if(allocated(ista_G_g1k)) deallocate(ista_G_g1k); allocate(ista_G_g1k(kv3))
    if(allocated(iend_G_g1k)) deallocate(iend_G_g1k); allocate(iend_G_g1k(kv3))

    if(allocated(nis_G_g1k)) deallocate(nis_G_g1k);  allocate(nis_G_g1k(0:nrank_g-1,kv3))
    do ik = 1, kv3
       nis_G_g1k(:,ik) = iba2(ik)
    end do

    if(allocated(nie_G_g1k)) deallocate(nie_G_g1k);  allocate(nie_G_g1k(0:nrank_g-1,kv3))
    nie_G_g1k = 0

    if(allocated(nel_G_g1k)) deallocate(nel_G_g1k);  allocate(nel_G_g1k(0:nrank_g-1,kv3))
    nel_G_g1k = 0

    if(allocated(map_G_g1k)) deallocate(map_G_g1k);  allocate(map_G_g1k(maxval(iba2),kv3))
    map_G_g1k = 0

    do ik = 1, kv3
       set_mapping_func = .true.
       if(ipri > 1 .and. printable) write(nfout,'(" !|| ik = ",i10 &
            & ," iba2(",i3,") = ",i10," <<m_Parallel_init_mpi_G_iba2_3D>>")') ik, ik, iba2(ik)
       call set_block_range(iba2(ik),nrank_g,nel_G_g1k(0,ik),nis_G_g1k(0,ik) &
            & ,nie_G_g1k(0,ik), set_mapping_func,map_G_g1k(1,ik))
       ista_G_g1k(ik) = nis_G_g1k(myrank_g,ik)
       iend_G_g1k(ik) = nie_G_g1k(myrank_g,ik)
       np_G_g1k(ik)   = nel_G_g1k(myrank_g,ik)
       mp_G_g1k(ik)   = maxval(nel_G_g1k(0:nrank_g-1,ik))
       if(ipri > 1 .and. printable) then
          write(nfout,'(" !|| --- mis_G_g1k,nie_G_g1k,nel_G_g1k ---")')
          write(nfout,'(" !|| ( rank_g  )",15i7)')(i,i=1,nrank_g)
          write(nfout,'(" !|| ( nis_G_g1k )",15i7)')(nis_G_g1k(i,ik),i=0,nrank_g-1)
          write(nfout,'(" !|| ( nie_G_g1k )",15i7)')(nie_G_g1k(i,ik),i=0,nrank_g-1)
          write(nfout,'(" !|| ( nel_G_g1k )",15i7)')(nel_G_g1k(i,ik),i=0,nrank_g-1)
          write(nfout,'(" !|| (myrank_g = ",i3," iba2(ik),ista_G_g1k,iend_G_g1k,np_G_g1k,mp_G_g1k = ",5i7)') &
               & myrank_g,iba2(ik),ista_G_g1k(ik),iend_G_g1k(ik),np_G_g1k(ik),mp_G_g1k(ik)
          write(nfout,'(" !|| ( map_G_g1k )",15i7)')(map_G_g1k(i,ik),i=nis_G_g1k(0,ik),nis_G_g1k(0,ik)+50)
! === DEBUG by tkato 2012/11/06 ================================================
          if(nrank_g > 1) then
! ==============================================================================
          write(nfout,'(" !|| ( map_G_g1k )",15i7)')(map_G_g1k(i,ik),i=nis_G_g1k(1,ik),nis_G_g1k(1,ik)+50)
! === DEBUG by tkato 2012/11/06 ================================================
          end if
! ==============================================================================
! === DEBUG by tkato 2012/12/05 ================================================
!         write(nfout,'(" !|| ( map_G_g1k )",15i7)')(map_G_g1k(i,ik),i=nis_G_g1k(nrank_g-1,ik),nis_G_g1k(nrank_g-1,ik)+50)
          write(nfout,'(" !|| ( map_G_g1k )",15i7)')(map_G_g1k(i,ik),i=nis_G_g1k(nrank_g-1,ik), &
        &                                                              min(nis_G_g1k(nrank_g-1,ik)+50,nie_G_g1k(nrank_g-1,ik)))
! ==============================================================================
       end if

    end do

  end subroutine m_Parallel_init_mpi_G_iba2_3D

  subroutine m_Parallel_init_mpi_B_iba2_3D(nfout,ipri,printable,kv3,iba2)
    integer, intent(in) :: nfout,ipri, kv3, iba2(kv3)
    logical, intent(in) :: printable

    logical             :: set_mapping_func
    integer             :: i, ik
    integer, parameter  :: nc = 24
    integer, parameter  :: nc1 = 16

    if(nrank_k > kv3) call phase_error_with_msg(nfout,' nrank_k > kv3 (m_Parallel_init_mpi_elec)',__LINE__,__FILE__)

! (( ista_B_g1k, iend_B_g1k ))
    if(allocated(np_B_g1k)) deallocate(np_B_g1k);  allocate(np_B_g1k(kv3))
    if(allocated(mp_B_g1k)) deallocate(mp_B_g1k);  allocate(mp_B_g1k(kv3))
    if(allocated(ista_B_g1k)) deallocate(ista_B_g1k); allocate(ista_B_g1k(kv3))
    if(allocated(iend_B_g1k)) deallocate(iend_B_g1k); allocate(iend_B_g1k(kv3))

    if(allocated(nis_B_g1k)) deallocate(nis_B_g1k);  allocate(nis_B_g1k(0:nrank_e-1,kv3))
    do ik = 1, kv3
       nis_B_g1k(:,ik) = iba2(ik)
    end do

    if(allocated(nie_B_g1k)) deallocate(nie_B_g1k);  allocate(nie_B_g1k(0:nrank_e-1,kv3))
    nie_B_g1k = 0

    if(allocated(nel_B_g1k)) deallocate(nel_B_g1k);  allocate(nel_B_g1k(0:nrank_e-1,kv3))
    nel_B_g1k = 0

    if(allocated(map_B_g1k)) deallocate(map_B_g1k);  allocate(map_B_g1k(maxval(iba2),kv3))
    map_B_g1k = 0

    do ik = 1, kv3
       set_mapping_func = .true.
       if(ipri > 1 .and. printable) write(nfout,'(" !|| ik = ",i10 &
            & ," iba2(",i3,") = ",i10," <<m_Parallel_init_mpi_B_iba2_3D>>")') ik, ik, iba2(ik)
       call set_block_range(iba2(ik),nrank_e,nel_B_g1k(0,ik),nis_B_g1k(0,ik) &
            & ,nie_B_g1k(0,ik), set_mapping_func,map_B_g1k(1,ik))
       ista_B_g1k(ik) = nis_B_g1k(myrank_e,ik)
       iend_B_g1k(ik) = nie_B_g1k(myrank_e,ik)
       np_B_g1k(ik)   = nel_B_g1k(myrank_e,ik)
       mp_B_g1k(ik)   = maxval(nel_B_g1k(0:nrank_e-1,ik))
       if(ipri > 1 .and. printable) then
          write(nfout,'(" !|| --- mis_B_g1k,nie_B_g1k,nel_B_g1k ---")')
          write(nfout,'(" !|| ( rank_e  )",15i7)')(i,i=1,nrank_e)
          write(nfout,'(" !|| ( nis_B_g1k )",15i7)')(nis_B_g1k(i,ik),i=0,nrank_e-1)
          write(nfout,'(" !|| ( nie_B_g1k )",15i7)')(nie_B_g1k(i,ik),i=0,nrank_e-1)
          write(nfout,'(" !|| ( nel_B_g1k )",15i7)')(nel_B_g1k(i,ik),i=0,nrank_e-1)
          write(nfout,'(" !|| (myrank_e = ",i3," iba2(ik),ista_B_g1k,iend_B_g1k,np_B_g1k,mp_B_g1k = ",5i7)') &
               & myrank_e,iba2(ik),ista_B_g1k(ik),iend_B_g1k(ik),np_B_g1k(ik),mp_B_g1k(ik)
          write(nfout,'(" !|| ( map_B_g1k )",15i7)')(map_B_g1k(i,ik),i=nis_B_g1k(0,ik),nis_B_g1k(0,ik)+50)
! === DEBUG by tkato 2012/11/06 ================================================
          if(nrank_e > 1) then
! ==============================================================================
          write(nfout,'(" !|| ( map_B_g1k )",15i7)')(map_B_g1k(i,ik),i=nis_B_g1k(1,ik),nis_B_g1k(1,ik)+50)
! === DEBUG by tkato 2012/11/06 ================================================
          end if
! ==============================================================================
          write(nfout,'(" !|| ( map_B_g1k )",15i7)')(map_B_g1k(i,ik),i=nis_B_g1k(nrank_e-1,ik),nis_B_g1k(nrank_e-1,ik)+50)
       end if

    end do

  end subroutine m_Parallel_init_mpi_B_iba2_3D

!!$  subroutine m_Parallel_init_mpi_atm2_3D(nfout,ipri,printable,natm2)
!!$    integer, intent(in) :: nfout,ipri,natm2
!!$    logical, intent(in) :: printable
!!$    integer :: iwork, i, npes, mype
!!$#ifdef __TIMER_SUB__
!!$  call timer_sta(1239)
!!$#endif
!!$
!!$    npes = nrank_g
!!$    mype = myrank_g
!!$    allocate(is_atm2(0:npes-1))
!!$    allocate(ie_atm2(0:npes-1))
!!$    allocate(nel_atm2(0:npes-1))
!!$    iwork = ( natm2 - 1 ) / npes + 1
!!$    if(ipri >= 1 .and. printable) then
!!$       write(nfout,'(" !|| << init_mpi_atm2 >>")')
!!$       write(nfout,'(" !|| natm2 = ",i12)') natm2
!!$       write(nfout,'(" !|| -- is_natm2, ie_natm2 --")')
!!$    end if
!!$    do i = 0, npes-1
!!$       is_atm2(i) = min(i*iwork+1, natm2+1)
!!$       ie_atm2(i) = min(is_atm2(i)+iwork-1, natm2)
!!$       nel_atm2(i) = ie_atm2(i) - is_atm2(i) + 1
!!$       if(ipri >= 1 .and. printable) write(nfout,'(" !|| ",2i12)') is_atm2(i),ie_atm2(i)
!!$    enddo
!!$    ista_atm2 = is_atm2(mype)
!!$    iend_atm2 = ie_atm2(mype)
!!$    np_atm2   = nel_atm2(mype)
!!$    mp_atm2   = maxval(nel_atm2)
!!$#ifdef __TIMER_SUB__
!!$  call timer_end(1239)
!!$#endif
!!$  end subroutine m_Parallel_init_mpi_atm2_3D

  subroutine m_Parallel_wd_npes_etc_3D(nfcntn)
    integer, intent(in) :: nfcntn
    if(mype==0) then
       write(nfcntn,*) tag_npes_etc
       write(nfcntn,'(4i10)') npes,nrank_e,nrank_g,nrank_k
    end if
  end subroutine m_Parallel_wd_npes_etc_3D

  subroutine make_ball_buff
    integer :: ibuf_temp, i

! for zaj_l_ball fs(ri)_l_ball
    if(allocated(ball_buff)) deallocate(ball_buff)
    allocate(ball_buff(0:nrank_e-1))
    if(allocated(ball_addr)) deallocate(ball_addr)
    allocate(ball_addr(0:nrank_e-1))
    ibuf_temp = np_e*np_fs
    call mpi_allgather(ibuf_temp,1,mpi_integer,ball_buff(0),1,mpi_integer,mpi_kg_world,ierr)
    ball_addr(0) = 0
    do i = 1,nrank_e-1
       ball_addr(i) = ball_addr(i-1) + ball_buff(i-1)
    enddo

  end subroutine make_ball_buff

  subroutine m_Parallel_init_mpi_paw_3D(nfout, ipri, natm, mmesh)
    integer, intent(in) :: nfout, ipri, natm, mmesh
    integer :: width, reminder
    integer :: i
! === For nrc decomposion. by takto 2012/12/07 =================================
    integer :: num_avail
! ==============================================================================
    logical, save :: firstcall=.true.

    myrank_natm = mod(mype,nrank_natm)
    myrank_nrc  = mype/nrank_natm

    if(.not.firstcall) then
      call mpi_comm_free(mpi_natm_world, ierr)
      call mpi_comm_free(mpi_nrc_world,  ierr)
    endif
!    if(firstcall)then
    call mpi_comm_split(MPI_CommGroup,myrank_natm,myrank_nrc, mpi_natm_world,ierr)
    call mpi_comm_split(MPI_CommGroup,myrank_nrc, myrank_natm,mpi_nrc_world, ierr)
!    endif

    allocate(is_natm (0:nrank_natm-1))
    allocate(ie_natm (0:nrank_natm-1))
    allocate(nel_natm(0:nrank_natm-1))

    allocate(is_nrc (0:nrank_nrc-1))
    allocate(ie_nrc (0:nrank_nrc-1))
    allocate(nel_nrc(0:nrank_nrc-1))

    width = natm/nrank_natm
    reminder = mod(natm,nrank_natm)
    do i = 0, nrank_natm-1
       nel_natm(i) = width
       if(reminder > i) nel_natm(i) = nel_natm(i) + 1
       is_natm(i) = width*i + min(i,reminder) + 1
       ie_natm(i) = is_natm(i) + nel_natm(i) - 1
    enddo
    ista_natm = is_natm(myrank_natm)
    iend_natm = ie_natm(myrank_natm)
    ne_natm   = nel_natm(myrank_natm)

! === For nrc decomposion. by takto 2012/12/07 =================================
    if(mmesh/nrank_nrc >= paw_nrc_min_elements) then
! ==============================================================================
    do i = 0, nrank_nrc-1
       width = mmesh/nrank_nrc
       nel_nrc(i) = width
       reminder = mod(mmesh,nrank_nrc)
       if(reminder > i) nel_nrc(i) = nel_nrc(i) + 1
       is_nrc(i) = width*i + min(i,reminder) + 1
       ie_nrc(i) = is_nrc(i) + nel_nrc(i) - 1
    enddo
! === For nrc decomposion. by takto 2012/12/07 =================================
       if(myrank_nrc == nrank_nrc-1) paw_last_rank_on_nrc = .true.
    else
       num_avail = mmesh/paw_nrc_min_elements
       if(num_avail == 0) num_avail = 1
       do i = 0, num_avail-1
          width = mmesh/num_avail
          nel_nrc(i) = width
          reminder = mod(mmesh,num_avail)
          if(reminder > i) nel_nrc(i) = nel_nrc(i) + 1
          is_nrc(i) = width*i + min(i,reminder) + 1
          ie_nrc(i) = is_nrc(i) + nel_nrc(i) - 1
       enddo
       do i = num_avail, nrank_nrc-1
          nel_nrc(i) = 0
          is_nrc(i) = 0
          ie_nrc(i) = 0
       end do
       if(myrank_nrc == num_avail-1) paw_last_rank_on_nrc = .true.
       if(myrank_nrc > num_avail-1) paw_no_work_to_do = .true.
       write(nfout,*) 'NOTE: Only ', num_avail, ' procs. are available for nrc decomposition!!!'
    end if
! ==============================================================================
    ista_nrc = is_nrc(myrank_nrc)
    iend_nrc = ie_nrc(myrank_nrc)
    ne_nrc   = nel_nrc(myrank_nrc)

! ====================================================================
    if(ipri>=2) then
    write(nfout,*) '===== PAW Decomposition!!! ====='
    write(nfout,*) 'npes: ', npes
    write(nfout,*) 'mype: ', mype
    write(nfout,*) 'nrank_natm, nrank_nrc:   ', nrank_natm,  nrank_nrc
    write(nfout,*) 'ista_natm, iend_natm, ne_natm: ', ista_natm, iend_natm, ne_natm
    write(nfout,*) 'myrank_natm, myrank_nrc: ', myrank_natm, myrank_nrc
    write(nfout,*) 'ista_nrc, iend_nrc, nel_nrc:   ', ista_nrc, iend_nrc, ne_nrc
    write(nfout,*) '----- natm Decomposed as... -----'
    write(nfout,*) 'natm:  ', natm
    do i = 0, nrank_natm-1
      write(nfout,*) is_natm(i), ie_natm(i), nel_natm(i)
    enddo
    write(nfout,*) '----- mmesh Decomposed as... -----'
    write(nfout,*) 'mmesh: ', mmesh
    do i = 0, nrank_nrc-1
      write(nfout,*) is_nrc(i), ie_nrc(i), nel_nrc(i)
    enddo
    write(nfout,*) '===== PAW Decomposition!!! ====='
    endif
! ====================================================================
    firstcall=.false.
  end subroutine m_Parallel_init_mpi_paw_3D

  subroutine index_div_fft(kimg, wk, ndim, nelm, iflg)
    integer, intent(in) :: kimg, ndim, nelm, iflg
    integer, intent(out) :: wk(2,ndim)
    integer :: nel, n, i
    integer :: tmp(ndim)

    wk(:,:) = 0
    tmp(:) = 0
    nel = nelm
    n = nel / ndim
    if ((kimg == 1) .and. (iflg == 1)) then
      if (mod(n,2) /=0) then
         n = n + 1
      end if
    end if
    tmp = 0
    do i = 1, ndim
       if (nel > 0) then
          if (nel >= n) then
             tmp(i) = n
             nel = nel - n
          else
             tmp(i) = nel
             nel = 0
          end if
       end if
    end do
    if (nel > 0) then
       n = ((nel-1) / ndim ) + 1
       if ((kimg == 1) .and. (iflg == 1)) then
         if (mod(n,2) /=0) then
            n = n + 1
         end if
       end if
       do i = 1, ndim
          if (nel > 0) then
             tmp(i) = tmp(i) + n
             nel = nel - n
          end if
       end do
    end if

    wk(1,1) = 1
    wk(2,1) = tmp(1)
    do i = 2, ndim
       if (tmp(i) == 0) cycle
       wk(1,i) = wk(2,i-1) + 1
       wk(2,i) = wk(2,i-1) + tmp(i)
    end do
  end subroutine index_div_fft

  subroutine decide_div_fft(nmrank, mindim1, mindim2, dim1, dim2)
    integer, intent(in)  :: nmrank, mindim1, mindim2
    integer, intent(out) :: dim1, dim2
    integer :: jsdim, jldim, ksdim, kldim, krem, is, il

       if ((mindim1*mindim2) <= nmrank) then
          dim1 = mindim1
          dim2 = mindim2
       else
          if (mindim1<mindim2) then
             jsdim = mindim1
             jldim = mindim2
          else
             jsdim = mindim2
             jldim = mindim1
          end if
          ksdim = 0
          kldim = 0
          krem  = 999999
          do is = 1, jsdim
             do il = is, jldim
                if ((is*il) <= nmrank) then
                   if ((nmrank-(is*il)) <= krem) then
                      ksdim = is
                      kldim = il
                      krem  = nmrank-(is*il)
                   end if
                end if
             end do
          end do
          if (mindim1<mindim2) then
             dim1 = ksdim
             dim2 = kldim
          else
             dim1 = kldim
             dim2 = ksdim
          end if
       end if
  end subroutine
! === For epsmain by tkato 2013/11/14 ==========================================
  subroutine m_Parallel_epsmain_reallocate()
     if(allocated(wf_fft_scnt)) deallocate(wf_fft_scnt)
     if(allocated(wf_fft_rcnt)) deallocate(wf_fft_rcnt)
     if(allocated(wf_fft_index)) deallocate(wf_fft_index)
     if(allocated(wf_fft_dist)) deallocate(wf_fft_dist)
     if(allocated(wf_fft_send)) deallocate(wf_fft_send)
     if(allocated(wf_fft_recv)) deallocate(wf_fft_recv)
     if(allocated(wf_fft_maxsend)) deallocate(wf_fft_maxsend)
     if(allocated(wf_fft_maxrecv)) deallocate(wf_fft_maxrecv)
     if(allocated(fft_wf_scnt)) deallocate(fft_wf_scnt)
     if(allocated(fft_wf_rcnt)) deallocate(fft_wf_rcnt)
     if(allocated(fft_wf_send)) deallocate(fft_wf_send)
     if(allocated(fft_wf_recv)) deallocate(fft_wf_recv)
     if(allocated(fft_wf_dist)) deallocate(fft_wf_dist)
     if(allocated(fft_wf_index)) deallocate(fft_wf_index)
     if(allocated(fft_wf_maxsend)) deallocate(fft_wf_maxsend)
     if(allocated(fft_wf_maxrecv)) deallocate(fft_wf_maxrecv)
  end subroutine m_Parallel_epsmain_reallocate
! ==============================================================================

  subroutine m_Parallel_dealloc_mpi_paw_3D()
    if(allocated(is_natm)) deallocate(is_natm)
    if(allocated(ie_natm)) deallocate(ie_natm)
    if(allocated(nel_natm)) deallocate(nel_natm)
    if(allocated(is_nrc)) deallocate(is_nrc)
    if(allocated(ie_nrc)) deallocate(ie_nrc)
    if(allocated(nel_nrc)) deallocate(nel_nrc)
  end subroutine m_Parallel_dealloc_mpi_paw_3D

!! TY 2019.06.25 -->
  subroutine m_Parallel_store_prev_np_g1k(kv3)
    integer, intent(in) :: kv3
    call m_Parallel_alloc_np_g1k_prev(kv3)
    if (.not.allocated(np_g1k)) return
    np_g1k_prev = np_g1k
  end subroutine m_Parallel_store_prev_np_g1k

  subroutine m_Parallel_alloc_np_g1k_prev(kv3)
    integer, intent(in) :: kv3
    if(allocated(np_g1k_prev)) deallocate(np_g1k_prev);  allocate(np_g1k_prev(kv3))
  end subroutine m_Parallel_alloc_np_g1k_prev
!! <--


  subroutine m_Parallel_resolve_decomp(nfout,kv3,done_something)
    integer, intent(in) :: nfout,kv3
    logical, intent(out) :: done_something
    integer :: i,neng,icount,itarget,nk,ne
    integer, allocatable, dimension(:) :: cands1,cands2
    real(kind=8) :: closest,fac
    done_something = .false.
    if(read_from_args) return
    done_something = .true.
    neng = 1;nk=npes
    if (npes>kv3) then
       do i=kv3,1,-1
          if(mod(npes,i)==0 .and. mod(kv3,i)==0) then
            neng = npes/i
            nk = npes/neng
            exit
          endif
       enddo
    else if(npes<kv3) then
       do i=npes,1,-1
          if(mod(npes,i)==0 .and. mod(kv3,i)==0) then
            neng = npes/i
            nk = npes/neng
            exit
          endif
       enddo
    endif
    ne = neng

    nrank_k = nk
    allocate(cands1(npes));allocate(cands2(npes))
    icount=0
    do i=1,max(neng/2,1)
       if(mod(neng,i)==0)then
          icount = icount+1
          cands1(icount) = i
          cands2(icount) = neng/i
       endif
    enddo
    closest = 1.d+30
    itarget = 1
    do i=1,icount
       fac = real(cands1(i))/real(cands2(i))
       if (abs(fac-0.5d0)<closest) then
          closest = abs(fac-0.5d0)
          itarget = i
       endif
    enddo

    nrank_e = cands1(itarget)
    nrank_g = cands2(itarget)

!    if(mype==0) then
!       write(nfout,'(a,3i5)') ' !** default ne nk ng : ',nrank_e,nrank_k,nrank_g
!    endif
    deallocate(cands1);deallocate(cands2)
    call m_Parallel_get_nproc_from_arg_3D(.true.,nfout,nrank_e,nrank_k,nrank_g)
  end subroutine m_Parallel_resolve_decomp

#ifdef KMATH_FFT3D
  subroutine m_Parallel_kmath3d_init(fft_box_size_WF, nstage_fft3d)
    use kmath_fft3d_mod
    use kmath_time_mod
    integer, intent(in) :: fft_box_size_WF(3,0:1)
    integer, dimension(3), intent(in) :: nstage_fft3d
    integer, dimension(3) :: bsize,nproc
    integer :: nx,ny,nz,ierr,tmp,tmpx,tmpy

    call KMATH_Time_Init()

    nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
    ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
    nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1
    call mpi_allreduce(mpi_in_place, nx, 1, mpi_integer, mpi_max, mpi_ke_world,ierr)
    call mpi_allreduce(mpi_in_place, ny, 1, mpi_integer, mpi_max, mpi_ke_world,ierr)
    call mpi_allreduce(mpi_in_place, nz, 1, mpi_integer, mpi_max, mpi_ke_world,ierr)

    bsize(:) = fft_box_size_WF(:,1)
    nproc(1) = nint(bsize(1)/real(nx))
    nproc(2) = nint(bsize(2)/real(ny))
    nproc(3) = nint(bsize(3)/real(nz))
    tmp = nint(bsize(1)/real(ny))

    kmath3d_box_size(1) = tmp      * ny
    kmath3d_box_size(2) = nproc(2) * ny
    kmath3d_box_size(3) = nproc(3) * nz
    nproc_fft3d_inverse = nproc

    call KMATH_FFT3D_Init(kmath3d_handle_wf_inverse, mpi_ke_world, fft_box_size_WF(:,1), &
   &     nproc, nstage_fft3d)

!    nx = xyz_fft_x(2,1)-xyz_fft_x(1,1)+1
!    ny = xyz_fft_x(2,2)-xyz_fft_x(1,2)+1
!    nz = xyz_fft_x(2,3)-xyz_fft_x(1,3)+1

!    call mpi_allreduce(mpi_in_place, nx, 1, mpi_integer, mpi_max, mpi_ke_world,ierr)
!    call mpi_allreduce(mpi_in_place, ny, 1, mpi_integer, mpi_max, mpi_ke_world,ierr)
!    call mpi_allreduce(mpi_in_place, nz, 1, mpi_integer, mpi_max, mpi_ke_world,ierr)

    nproc(1) = nint(bsize(1)/real(nx))
    nproc(2) = nint(bsize(2)/real(ny))
    nproc(3) = nint(bsize(3)/real(nz))
    !write(0,*) 'Nproc',nproc
!    tmpx = nproc(1)
!    tmpy = nproc(2)
    nproc_fft3d_direct = nproc
!    nproc(1) = tmpy;nproc(2)=tmpx
!    nproc = 1
!    nproc(3) = 2
!    nproc(2) = 3
    call KMATH_FFT3D_Init(kmath3d_handle_wf_direct, mpi_ke_world, fft_box_size_WF(:,1), &
   &     nproc, nstage_fft3d)
  end subroutine m_Parallel_kmath3d_init

  subroutine m_Parallel_kmath3d_finalize()
    use kmath_fft3d_mod
    use kmath_time_mod
    call KMATH_FFT3D_Finalize(kmath3d_handle_wf_inverse)
    call KMATH_FFT3D_Finalize(kmath3d_handle_wf_direct)
    call KMATH_Time_Finalize
  end subroutine m_Parallel_kmath3d_finalize
#endif

end module m_Parallelization
