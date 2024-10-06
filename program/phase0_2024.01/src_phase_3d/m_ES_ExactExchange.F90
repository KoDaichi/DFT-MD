#define HYBRID_DGEMM

#define DERIVATIVE_IN_CARTS
#define USE_DERIVATIVE_QRSPS_RSPACE
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 616 $)
!
!  MODULE: m_ES_ExactExchange
!
!  AUTHOR(S): T. Yamamoto   November/22/2009
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

module m_ES_ExactExchange
  use m_Electronic_Structure,only: totch,zaj_l,occup_l,neordr,eko_l,vnlph_l,fsr_l,fsi_l
  use m_NonLocal_Potential, only : snl
  use m_PlaneWaveBasisSet,  only : ngabc,igf,kg1,kg,kgp,nbase,nbmx,iba,m_pwBS_kinetic_energies &
       &                         , nbase_gamma,igfp_l,ngpt_l,m_pwBS_sphrp_exx,kgp_exx, igfp_exx, n_rGpv &
       &                         , igf_exx,kg_exx, kg_gamma, get_kg1_ext
  use m_Kpoints,            only : kv3,vkxyz,k_symmetry,qwgt,mp_index &
       &                         , m_Kp_sample_mesh, kshift &
       &                         , kv3_ek, vkxyz_ek, qwgt_ek
  use m_FFT,                only : nfft,fft_box_size_WF, fft_box_size_EXX &
       &                         , m_FFT_alloc_WF_work &
       &                         , m_FFT_dealloc_WF_work &
       &                         , m_FFT_WF, fft_box_size_CD, m_FFT_CD0 &
       &                         , m_FFT_CD0_exx &
#ifdef FFTW3
       &                         , nfftp_exx_nonpara, fft_box_size_CD_exx, m_FFT_exx, nfft_exx
#else
       &                         , nfftp_exx_nonpara, fft_box_size_CD_exx, nfft_exx
#endif
  use m_Ionic_System,       only : zeta1, qex, iatom, ntyp, natm, ityp, pos, iwei, napt
  use m_PseudoPotential,    only : qitg_l, iqitg, nqitg, m_PP_include_vanderbilt_pot &
       &                         , ilmt,iqitg,nloc,ntau,ival,isph,il2p,dl2p,modnrm &
       &                         , nlmta,lmta,lmtt,ltp,mtp,taup,nloc,m_PP_find_maximum_l &
       &                         , m_PP_tell_lmtt_l_m_tau &
       &                         , m_PP_set_index_arrays1,m_PP_set_index_arrays2 &
       &                         , radr,xh,rmax,nmesh,mmesh, qrspspw, nlmtt &
       &                         , m_PP_qitgft_qmk

  use m_Files,              only : nfout, m_Files_open_nfhybridinfo, m_Files_close_nfhybridinfo, nfhybridinfo &
       &                         , m_Files_open_nfscfzaj, m_Files_close_nfscfzaj, nfscfzaj
  use m_Timing,             only : tstatc0_begin,tstatc0_end
  use m_Control_Parameters, only : nspin,ipri,kimg,neg,printable,af,ekmode &
       &                         , alpha_exx,omega_exx,sw_screened_exchange  &
       &                         , sw_singular_correction,reduction_factor_exx &
       &                         , way_ksample,m_CtrlP_cachesize,sw_rspace_hyb &
       &                         , sw_precalculate,sw_rsb,potential_update &
       &                         , sw_rspace_hyb_dgm &
       &                         , icond, sw_hybrid_functional, ipriexx &
       &                         , precision_WFfile, fixed_charge_k_parallel
  use m_Crystal_Structure,  only : univol,altv,rltv,nopr,tau,op
  ! nopr : # of operations
  ! tau  : nonprimitive translation vector
  ! op   : rotation matrix in real space
  use m_Const_Parameters,   only : DP, CMPLDP, SKIP, EXECUT, ON, OFF, DIRECT, INVERSE, DELTA, DELTA10 &
       &                         , DELTA07, PAI2, PAI4, PAI, ELECTRON, GAMMA, BUCS &
       &                         , MONKHORST_PACK, MESH &
       &                         , INITIAL, CONTINUATION, FIXED_CHARGE_CONTINUATION, FIXED_CHARGE &
       &                         , CRDTYP, SP, ONE_BY_ONE, ALL_AT_ONCE
  use m_Parallelization,    only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype &
       &                         , nrank_e,nrank_k,myrank_e,map_e,ista_e,iend_e,istep_e,idisp_e &
       &                         , map_z,np_e,mpi_k_world,mpi_e_world,myrank_k,map_k,ista_k,iend_k &
       &                         , ista_kg1_k, np_kg1_k, mp_kg1_k &
!       &                         , m_Parallel_mpi_nval,np_nval, mp_nval,ista_kngp_exx,iend_kngp_exx &
       &                         , ista_kngp_exx,iend_kngp_exx &
#ifdef TRANSPOSE
       &                         , ierr,mp_e,nis_e,nie_e,nel_e
#else
       &                         , ierr,mp_e,nel_e
#endif

! ============================ KT_Test ================================ 12.5Exp
  use m_Control_Parameters,    only : nmax_G_hyb, hybrid_functional_type
  use m_PlaneWaveBasisSet,     only : kgp, kg1, kg
! ===================================================================== 12.5Exp
  use m_PlaneWaveBasisSet, only : ngabc_kngp_l, kg_gamma
  use m_Parallelization, only : mpi_ke_world, mpi_kg_world, mpi_ge_world, np_g1k, ista_g1k, iend_g1k &
                              , ista_fs, iend_fs, np_fs, nrank_g, nel_fft_x, nel_fft_y, nel_fft_z &
                              , xyz_fft_x, xyz_fft_y, xyz_fft_z, mp_g1k, kg1_ext, myrank_g, map_fft_x, map_fft_y &
       &                          , fft_X_x_nel, fft_X_y_nel, fft_X_z_nel &
                              , is_kngp, ie_kngp, mp_fft_x, nel_kngp, np_kngp &
       &                      , m_Parallel_init_mpi_nval &
       &                      , ista_nval,iend_nval,np_nval,mp_nval,map_nval,map_z_nval,myrank_nval
  use m_FFT, only : m_FFT_Inverse_3D, m_FFT_Direct_3D
! === FFT xzy ==================================================================
  use m_FFT, only : m_FFT_Inverse_XYZ_3D, m_FFT_Direct_XYZ_3D
  use m_Control_Parameters,  only : sw_fft_xzy,sw_change_axis
! ==============================================================================

  use m_Realspace, only : nmesh_rs_aug_max,nmesh_rs_aug,meshx_rs_aug,meshy_rs_aug,meshz_rs_aug,meshxyz_rs_aug &
  & , qr_clm_ylm,dqr_clm_ylm,plmt1,plmt2,nlmtpair
! ==== ASMS ====
  use m_Kpoints,  only : kv3_fbz,  vkxyz_fbz
! ==== ASMS ====

  use m_IterationNumbers, only : nk_in_the_process
  use mpi

  implicit none

  integer :: sw_update_wfv = ON
  integer :: nval
  integer :: nval_old
  integer :: nfftwf
  integer :: ntrs ! ntrs = 1, if TRS is used. Otherwise, ntrs=0
  integer, allocatable :: ngpt_exx(:,:,:) ! d(kg,nopr,0:ntrs)
  real(kind=DP), allocatable :: qitg_exx(:,:,:) ! d(kgp,nqitg,nqmk)
  real(kind=DP), allocatable :: ylm_exx(:,:,:) ! d(kgp,maxylm,nqmk)
  real(kind=DP), allocatable :: vc(:,:) ! d(kgp,nqmk)
  real(kind=DP), allocatable :: wfv(:,:,:,:) ! d(kg1,nval,kv3,kimg)
  real(kind=DP), allocatable :: wfv_old(:,:,:,:) ! d(kg1,nval_old,kv3,kimg)
  real(kind=DP), allocatable :: fsr_exx(:,:,:) ! d(nval,nlmta,kv3)
  real(kind=DP), allocatable :: fsi_exx(:,:,:) ! d(nval,nlmta,kv3)
  real(kind=DP), allocatable :: fsr_exx_old(:,:,:) ! d(nval,nlmta,kv3)
  real(kind=DP), allocatable :: fsi_exx_old(:,:,:) ! d(nval,nlmta,kv3)
  real(kind=DP), allocatable :: occup_val(:,:) ! d(nval,kv3)
  real(kind=DP), allocatable :: occup_val_old(:,:) ! d(nval_old,kv3)
  real(kind=DP), allocatable :: dfsr_exx(:,:,:,:) ! d(nval,nlmta,kv3,3)
  real(kind=DP), allocatable :: dfsi_exx(:,:,:,:) ! d(nval,nlmta,kv3,3)
  real(kind=DP), allocatable :: dfsr_l(:,:,:,:) ! d(np_e,nlmta,ista_k:iend_k,3)
  real(kind=DP), allocatable :: dfsi_l(:,:,:,:) ! d(np_e,nlmta,ista_k:iend_k,3)

  integer :: ikcenter
  real(kind=DP) :: kshift_norm2
  logical, allocatable :: qmk_zero(:,:) ! d(kv3bz,kv3)
  logical, allocatable :: qmk_zero2(:) ! d(kv3bz*kv3)
  real(kind=DP), allocatable :: qmk(:,:) !d(kv3bz*kv3,3)
  integer, allocatable :: iqmk(:,:) !d(kv3bz,kv3))
  integer :: nqmk
  integer :: kv3bz

  real(kind=DP), allocatable :: eexx_kb(:,:)
  !
  ! q(BZ) = S*k(IBZ)
  !
  ! Symmetry operations in space group:
  !    Os={S|u}
  !    Or={R|t}
  !
  ! Wavefunction coefficient:
  !    c(q,SG) = c(k,G) exp(i(k+G)*t)
  !            = c(k,G) exp(-i(q+SG)*u)
  !
  !    If time reversal symmetry is present,
  !         c(-q,-G) = c*(q,G)
  !
  integer, allocatable :: ikp(:) ! d(kv3bz)
  integer, allocatable :: iops(:) ! d(kv3bz) No. of S_inv = R :  Or={R|t}
  integer, allocatable :: iopt(:) ! d(kv3bz) No. of S         :  Os={S|u}
  integer, allocatable :: itrs(:) ! d(kv3bz) Time reversal symmetry, 0 or 1
  real(kind=DP), allocatable :: vkbz(:,:) ! d(kv3bz,3) k-points in Full BZ
  real(kind=DP), allocatable :: wbz(:) ! d(kv3bz) weight for k-points in Full BZ

  real(kind=DP), allocatable :: oprec(:,:,:) ! d(3,3,nopr) Reciprocal space rotation matrix
  real(kind=DP), allocatable :: opdir(:,:,:) ! d(3,3,nopr) Direct space rotation matrix

  integer :: nface
  integer, allocatable :: iface_bz(:,:) ! d(3,nface)
  integer :: n_bz_mesh(3)

  real(kind=DP), allocatable :: crotylm(:,:,:) ! d(mmax,nsph,nopr)
  integer, allocatable :: iylm(:,:,:) ! d(mmax,nsph,nopr)
  integer, allocatable :: nylm(:,:) ! d(nsph,nopr)
  integer, allocatable :: isph_lmta(:) ! d(nlmta))
  integer, allocatable :: lmta_rot(:,:,:) ! d(mmax,nlmta,nopr))
  integer, allocatable :: ia_lmta(:) !(nlmta)

  logical, save :: preproc_done=.false.
  integer :: n,ibsize,iwidth,maxm
  integer, parameter :: mcritical = 4*2+1
  integer, allocatable :: il3(:)
  integer, allocatable, dimension(:) :: nqitg_sp,nqitg_sp0
  integer, allocatable, dimension(:) :: iq2l3 ! d(nqitg)
  integer, allocatable, dimension(:,:) :: nc  ! d(maxm,nqitg)
  integer :: mc ! maxval(nc)
  integer, allocatable, dimension(:,:,:) :: nc2lmt1, nc2lmt2, nc2n ! d(mc,maxm,nqitg)
  integer, allocatable, dimension(:,:,:,:)   :: wf_fft_scnt_exx,    wf_fft_rcnt_exx
  integer, allocatable, dimension(:,:,:,:,:) :: wf_fft_send_exx,    wf_fft_recv_exx
  integer, allocatable, dimension(:,:,:,:)   :: wf_fft_dist_exx,    wf_fft_index_exx
  integer, allocatable, dimension(:,:,:)     :: wf_fft_maxsend_exx, wf_fft_maxrecv_exx

  integer, allocatable, dimension(:)   :: fft_wf_scnt_rhog,    fft_wf_rcnt_rhog
  integer, allocatable, dimension(:,:) :: fft_wf_send_rhog,    fft_wf_recv_rhog
  integer, allocatable, dimension(:)   :: fft_wf_dist_rhog,    fft_wf_index_rhog
  integer                              :: fft_wf_maxsend_rhog, fft_wf_maxrecv_rhog

  integer, allocatable, dimension(:)   :: wf_fft_scnt_rhog,    wf_fft_rcnt_rhog
  integer, allocatable, dimension(:,:) :: wf_fft_send_rhog,    wf_fft_recv_rhog
  integer, allocatable, dimension(:)   :: wf_fft_dist_rhog,    wf_fft_index_rhog
  integer                              :: wf_fft_maxsend_rhog, wf_fft_maxrecv_rhog

  integer, allocatable, dimension(:) :: mapg2lx, mapg2ly


  real(kind=DP), allocatable, dimension(:,:) :: eko_hyb
  real(kind=DP), allocatable, dimension(:,:,:,:,:) :: fsrqm,fsiqm

  real(kind=DP), allocatable, dimension(:,:,:,:) :: exx_potential
  integer, save :: id_sname_cdfft = -1
  integer, allocatable, dimension(:,:) :: ngabc_red

  integer, pointer :: kv3_for_innerloop, kv3_for_outerloop
  integer, pointer :: kg1_for_innerloop
  real(kind=DP), dimension(:), pointer :: qwgt_for_innerloop,qwgt_for_outerloop
  real(kind=DP), dimension(:,:,:), pointer :: vkxyz_for_innerloop, vkxyz_for_outerloop
  integer, dimension(:), pointer :: iba_for_innerloop
  integer, dimension(:,:), pointer :: nbase_for_innerloop,nbase_gamma_for_innerloop

!  include 'mpif.h'                                      ! MPI

contains

  function k_index(ik) result(res)
    integer, intent(in) :: ik
    integer :: res
    res = ik
    if((icond==FIXED_CHARGE .or. icond==FIXED_CHARGE_CONTINUATION) .and. &
    & fixed_charge_k_parallel == ONE_BY_ONE) then
      res = nk_in_the_process + ik-1
      if (res>kv3_ek) res = kv3_ek
    endif
    return
  end function k_index

  subroutine m_ES_EXX_wd_hybridinfo()
    if(mype == 0)then
      write(nfout,'(a)') '!** output hybrid info to F_HYBRIDINFO'
    endif
    if(sw_hybrid_functional == OFF)then
      call m_ES_EXX_init()
    endif
    call m_Files_open_nfhybridinfo()
    if(mype==0)then
      write(nfhybridinfo) kv3,kv3bz
      write(nfhybridinfo) get_kg1_ext()
      write(nfhybridinfo) kshift
      write(nfhybridinfo) vkxyz
      write(nfhybridinfo) iba
      write(nfhybridinfo) nbase
      write(nfhybridinfo) allocated(nbase_gamma)
      if (allocated(nbase_gamma)) then
      write(nfhybridinfo) nbase_gamma
      endif
      write(nfhybridinfo) qwgt
      write(nfhybridinfo) vkbz
      write(nfhybridinfo) wbz
      write(nfhybridinfo) n_bz_mesh
      write(nfhybridinfo) ikp
      write(nfhybridinfo) iops
      write(nfhybridinfo) iopt
      write(nfhybridinfo) itrs
    endif
    call m_ES_EXX_gather_valence_states(nfout)
    if(mype==0)then
      write(nfhybridinfo) nval
      write(nfhybridinfo) neg
      write(nfhybridinfo) modnrm
      write(nfhybridinfo) nlmta
      write(nfhybridinfo) kg1
      write(nfhybridinfo) occup_val
      !!$write(nfhybridinfo) wfv
      if(modnrm == EXECUT) then
        write(nfhybridinfo) fsr_exx
        write(nfhybridinfo) fsi_exx
      endif
    endif
    call m_Files_close_nfhybridinfo()
  end subroutine m_ES_EXX_wd_hybridinfo

  subroutine m_ES_EXX_rd_hybridinfo()
    logical :: logi
    integer :: ierr,isize,mdn
    integer :: i,j,negp
    real(kind=DP), allocatable, dimension(:,:) :: wfd
    real(kind=SP), allocatable, dimension(:,:) :: wfs
    integer :: ip0,ip1
    real(kind=SP), allocatable, dimension(:,:,:,:) :: wfvt
    real(kind=DP), allocatable, dimension(:,:,:,:) :: wfvtd
    call m_Files_open_nfhybridinfo()
    allocate(kv3_for_innerloop)
    allocate(kg1_for_innerloop)
    if(mype==0)then
      read(nfhybridinfo) kv3_for_innerloop,kv3bz
      read(nfhybridinfo) isize
    endif
    call mpi_bcast(kv3_for_innerloop,1,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(kv3bz,1,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(isize,1,mpi_integer,0,MPI_CommGroup,ierr)
    kg1_for_innerloop = isize
    allocate(vkxyz_for_innerloop(kv3_for_innerloop,3,CRDTYP))
    allocate(iba_for_innerloop(kv3_for_innerloop))
    allocate(nbase_for_innerloop(isize,kv3_for_innerloop))
    logi = .false.
    if(mype == 0)then
      read(nfhybridinfo) kshift
      read(nfhybridinfo) vkxyz_for_innerloop
      read(nfhybridinfo) iba_for_innerloop
      read(nfhybridinfo) nbase_for_innerloop
      read(nfhybridinfo) logi
    endif
    call mpi_bcast(kshift,3,mpi_double_precision,0,MPI_CommGroup,ierr)
    call mpi_bcast(vkxyz_for_innerloop,kv3_for_innerloop*3*CRDTYP,mpi_double_precision,0,MPI_CommGroup,ierr)
    call mpi_bcast(iba_for_innerloop,kv3_for_innerloop,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(nbase_for_innerloop,isize*kv3_for_innerloop,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(logi,1,mpi_logical,0,MPI_CommGroup,ierr)
    if(logi)then
      allocate(nbase_gamma_for_innerloop(kg_gamma,2))
      if(mype == 0)then
        read(nfhybridinfo) nbase_gamma_for_innerloop
      endif
      call mpi_bcast(nbase_gamma_for_innerloop,kg_gamma*2,mpi_integer,0,MPI_CommGroup,ierr)
    endif
    allocate(qwgt_for_innerloop(kv3_for_innerloop))
    allocate(vkbz(kv3bz,3))
    allocate(wbz(kv3bz))
    allocate(ikp(kv3bz))
    allocate(iops(kv3bz))
    allocate(iopt(kv3bz))
    allocate(itrs(kv3bz))
    if(mype==0) then
      read(nfhybridinfo) qwgt_for_innerloop
      read(nfhybridinfo) vkbz
      read(nfhybridinfo) wbz
      read(nfhybridinfo) n_bz_mesh
      read(nfhybridinfo) ikp
      read(nfhybridinfo) iops
      read(nfhybridinfo) iopt
      read(nfhybridinfo) itrs
      read(nfhybridinfo) nval
      read(nfhybridinfo) negp
      read(nfhybridinfo) mdn
      read(nfhybridinfo) nlmta
      read(nfhybridinfo) isize
    endif
    call mpi_bcast(qwgt_for_innerloop,kv3_for_innerloop,mpi_double_precision,0,MPI_CommGroup,ierr)
    call mpi_bcast(vkbz,kv3bz*3,mpi_double_precision,0,MPI_CommGroup,ierr)
    call mpi_bcast(wbz,kv3bz,mpi_double_precision,0,MPI_CommGroup,ierr)
    call mpi_bcast(n_bz_mesh,3,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(ikp,kv3bz,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(iops,kv3bz,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(iopt,kv3bz,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(itrs,kv3bz,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(nval,1,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(negp,1,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(mdn,1,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(nlmta,1,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(isize,1,mpi_integer,0,MPI_CommGroup,ierr)
    ntrs = maxval(itrs)
    allocate(occup_val(nval,kv3_for_innerloop))
    if(mdn == EXECUT) then
      allocate(fsr_exx(nval,nlmta,kv3_for_innerloop))
      allocate(fsi_exx(nval,nlmta,kv3_for_innerloop))
    endif
    if(mype==0)then
      read(nfhybridinfo) occup_val
!!      read(nfhybridinfo) wfv
      if(mdn == EXECUT) then
        read(nfhybridinfo) fsr_exx
        read(nfhybridinfo) fsi_exx
      endif
    endif
    call mpi_bcast(occup_val,nval*kv3_for_innerloop,mpi_double_precision,0,MPI_CommGroup,ierr)
    if(mdn==EXECUT)then
      call mpi_bcast(fsr_exx,nval*nlmta*kv3_for_innerloop,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(fsi_exx,nval*nlmta*kv3_for_innerloop,mpi_double_precision,0,MPI_CommGroup,ierr)
    endif
    call m_Files_close_nfhybridinfo()

    if(sw_change_axis==ON)then
       call m_Parallel_init_mpi_nval(nfout,ipri,printable,nval)
    else
       call phase_error_with_msg(nfout,'fixed-charge-hybrid calculation is unsupported when sw_change_axis = off'&
                                ,__LINE__,__FILE__)
    endif

    call m_Files_open_nfscfzaj()
    allocate(wfv(isize,mp_nval,kv3_for_innerloop,kimg));wfv=0.d0
    if(precision_WFfile == SP) then
      allocate(wfvt(isize,mp_e,kv3_for_innerloop,kimg));wfvt=0.d0
    else
      allocate(wfvtd(isize,mp_e,kv3_for_innerloop,kimg));wfvt=0.d0
    endif
    if(precision_WFfile == SP) then
      allocate(wfs(isize,kimg));wfs=0.d0
      do i=1,kv3_for_innerloop
         do j=1,negp
            if(mype==0) read(nfscfzaj) wfs
            call mpi_bcast(wfs,isize*kimg,mpi_real,0,MPI_CommGroup,ierr)
            if(map_e(j) == myrank_e) wfvt(:,map_z(j),i,:) = wfs(:,:)
         enddo
      enddo
      deallocate(wfs)
    else
      allocate(wfd(isize,kimg));wfd=0.d0
      do i=1,kv3_for_innerloop
         do j=1,negp
            read(nfscfzaj) wfd
            if(map_z_nval(j).le.nval) wfv(:,map_z_nval(j),i,:) = wfd(:,:)
         enddo
      enddo
      deallocate(wfd)
    endif

    do i=1,kv3_for_innerloop
       do j=1,nval
          ip0=map_e(j)
          ip1=map_nval(j)
          if(ip0==myrank_e .and. ip1==myrank_nval) &
         & wfv(:,map_z_nval(j),i,:) = wfvt(:,map_z(j),i,:)
       enddo
    enddo
    if(precision_WFfile == SP) then
      deallocate(wfvt)
    else
      deallocate(wfvtd)
    endif
    call mpi_allreduce(MPI_IN_PLACE,wfv,isize*mp_nval*kv3_for_innerloop*kimg,&
       & mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
    call m_Files_close_nfscfzaj()
  end subroutine m_ES_EXX_rd_hybridinfo

! ====================== KT_Test ==================================== 12.5
  subroutine m_ES_EXX_set_nmax_G_hyb
    if ( nmax_G_hyb == 0 ) then
       nmax_G_hyb = kgp
    else if ( nmax_G_hyb == -1 ) then
       nmax_G_hyb = kg
    else if ( nmax_G_hyb == -2 ) then
       nmax_G_hyb = kg /2
    else if ( nmax_G_hyb == -3 ) then
       nmax_G_hyb = kg1
    endif
    if(sw_rspace_hyb==ON) nmax_G_hyb = kgp_exx
    write(nfout,*) 'nmax_G_hyb is set to ', nmax_G_hyb

  end subroutine m_ES_EXX_set_nmax_G_hyb
! ================================================================== 12.5

  subroutine m_ES_EXX_update(on_or_off)
    implicit none
    integer, intent(in) :: on_or_off

    sw_update_wfv = on_or_off
    if(printable) then
       if(sw_update_wfv==ON) then
          write(nfout,'(" EXX: wfv will be updated.")')
       else
          write(nfout,'(" EXX: wfv will be conserved.")')
       end if
    end if
  end subroutine m_ES_EXX_update

  subroutine m_ES_EXX_kbz
    implicit none
    integer :: ik, jk, iopr, i, j, ikbz, jkbz, jtrs
    integer :: i1, i2, i3
    real(kind=DP) :: rkxyz(3), kxyz(3), df(3), ko(3)
    logical :: ltrue
    integer :: g0(3)
    real(kind=DP) :: x(3), y(3)
    real(kind=DP) :: wtot
    real(kind=DP) :: ss(3,3)
    real(kind=DP) :: dk(3)

    allocate(oprec(3,3,nopr))
    do iopr=1,nopr
       do j=1,3
          x = matmul(op(:,:,iopr),rltv(:,j))
          do i=1,3
             y = altv(:,i)
             oprec(i,j,iopr) = dot_product(y,x)/PAI2
          end do
       end do
    end do

    if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION ) return

    write(nfout,'("k-points in IBZ")')
    write(nfout,'("ik, k1, k2, k3")')
    do ik=1,kv3_for_innerloop,nspin
       write(nfout,'(i5,3(1x,f10.5))') (ik-1)/nspin+1, vkxyz_for_innerloop(ik,1:3,BUCS)
    end do

    call m_ES_EXX_fbz_faces ! -> iface_bz

    if(kv3_for_innerloop/nspin==1) then
       kv3bz = 1
       n_bz_mesh(1:3) = 1
    else
       if(way_ksample == MONKHORST_PACK) then
          n_bz_mesh(1:3) = mp_index(1:3)
       else
          call m_Kp_sample_mesh(n_bz_mesh)
       end if
       kv3bz = n_bz_mesh(1)*n_bz_mesh(2)*n_bz_mesh(3)
    end if

    allocate(vkbz(kv3bz,3))
    allocate(wbz(kv3bz))
    wbz(1:kv3bz) = product(reduction_factor_exx(1:3))/dble(kv3bz)

    ikbz = 0
    do i1=0,n_bz_mesh(1)-1
       do i2=0,n_bz_mesh(2)-1
          do i3=0,n_bz_mesh(3)-1
             rkxyz(1) = (dble(i1) + kshift(1)) / dble(n_bz_mesh(1))  ! USAMI
             rkxyz(2) = (dble(i2) + kshift(2)) / dble(n_bz_mesh(2))
             rkxyz(3) = (dble(i3) + kshift(3)) / dble(n_bz_mesh(3))
             call m_ES_EXX_move_k_into_fbz(rkxyz,ko)
             do ik=1,kv3_for_innerloop
                dk(1:3) = vkxyz_for_innerloop(ik,1:3,BUCS)-ko(1:3)
                if( abs(dk(1)-nint(dk(1)))<DELTA07 .and. &
                  & abs(dk(2)-nint(dk(2)))<DELTA07 .and. &
                  & abs(dk(3)-nint(dk(3)))<DELTA07 ) then
                  ko(1:3) = vkxyz_for_innerloop(ik,1:3,BUCS)
                  exit
                end if
             end do
             ikbz = ikbz + 1
             vkbz(ikbz,1:3) = ko(1:3)
          end do
       end do
    end do

    write(nfout,'("k-points in BZ")')
    write(nfout,'("ikbz k1 k2 k3 wbz")')
    do ikbz=1,kv3bz
       write(nfout,'(i7,4(1x,f20.5))') ikbz, vkbz(ikbz,1:3), wbz(ikbz)
    end do

    allocate(ikp(kv3bz))
    allocate(iops(kv3bz))
    allocate(iopt(kv3bz)); iopt = 0
    allocate(itrs(kv3bz)); itrs = 0

    ! q(bz) = S*k'(ibz)
    do ikbz=1,kv3bz
       ikp(ikbz) = 0
       LOOP_S: do jtrs=0,1
          do iopr=1,nopr
             rkxyz = (1-2*jtrs) * matmul(oprec(:,:,iopr),vkbz(ikbz,1:3))
             do jk=1,kv3_for_innerloop,nspin
                dk(1:3) = vkxyz_for_innerloop(jk,1:3,BUCS)-rkxyz(1:3)
                if( abs(dk(1)-nint(dk(1)))<DELTA07 .and. &
                  & abs(dk(2)-nint(dk(2)))<DELTA07 .and. &
                  & abs(dk(3)-nint(dk(3)))<DELTA07 ) then
                   g0(1:3) = nint(dk(1:3))
                   ikp(ikbz) = jk
                   iops(ikbz) = iopr
                   itrs(ikbz) = jtrs
                  exit LOOP_S
                end if
             end do
          end do
       end do LOOP_S
       if(ikp(ikbz) == 0) then
          write(nfout,'("ikbz=",i5)') ikbz
          write(nfout,'(" ikp=",i5)') ikp(ikbz)
          write(nfout,'("rkxyz=",3(1x,f20.5))') rkxyz(1:3)
          call phase_error_with_msg(nfout,"Did not find the index of k' in IBZ.",__LINE__,__FILE__)
       end if

       do iopr=1,nopr
          ss = matmul(op(:,:,iopr),op(:,:,iops(ikbz)))
          if(abs(ss(1,1)-1.d0)<DELTA07 .and. &
           & abs(ss(2,2)-1.d0)<DELTA07 .and. &
           & abs(ss(3,3)-1.d0)<DELTA07 .and. &
           & abs(ss(1,2))<DELTA07 .and. &
           & abs(ss(1,3))<DELTA07 .and. &
           & abs(ss(2,3))<DELTA07) then
             iopt(ikbz) = iopr
             exit
          end if
       end do

       !   k'(bz) is transferable to k in IBZ by a symmetry operation.
       !   g0 = k(ibz) - R k(bz)
       !    0 = k(ibz) - R k'(bz)
       !    k'(bz) = k(bz) + R^inv g0

       rkxyz = (1-2*itrs(ikbz)) * matmul(oprec(:,:,iopt(ikbz)),g0(1:3))
       vkbz(ikbz,1:3) = vkbz(ikbz,1:3) + rkxyz(1:3)

    end do

    ntrs = maxval(itrs)

    write(nfout,'("k-points in BZ")')
    write(nfout,'("ikbz k1 k2 k3 wbz")')
    do ikbz=1,kv3bz
       write(nfout,'(i7,4(1x,f20.5))') ikbz, vkbz(ikbz,1:3), wbz(ikbz)
    end do

    write(nfout,'(2x,"ikbz",3x,"ikp",2x,"iops",2x,"iopt",2x,"itrs")')
    do ikbz=1,kv3bz
       write(nfout,'(i5,4(1x,i5))') ikbz, ikp(ikbz), iops(ikbz), iopt(ikbz), itrs(ikbz)
    end do
  end subroutine m_ES_EXX_kbz

  subroutine m_ES_EXX_kbz2
    implicit none
    integer :: ik, jk, iopr, i, j, ikbz, jkbz, jtrs
    integer :: i1, i2, i3
    real(kind=DP) :: rkxyz(3), kxyz(3), df(3), ko(3)
    logical :: ltrue
    integer :: g0(3), jtrs_max
    real(kind=DP) :: x(3), y(3)
    real(kind=DP) :: wtot
    real(kind=DP) :: ss(3,3)
    real(kind=DP) :: dk(3)

    if(.not.allocated(oprec)) allocate(oprec(3,3,nopr))
    do iopr=1,nopr
       do j=1,3
          x = matmul(op(:,:,iopr),rltv(:,j))
          do i=1,3
             y = altv(:,i)
             oprec(i,j,iopr) = dot_product(y,x)/PAI2
          end do
       end do
    end do

    if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION ) return

    if(kv3_for_innerloop/nspin==1) then
       n_bz_mesh(1:3) = 1
    else
       if(way_ksample == MONKHORST_PACK) then
          n_bz_mesh(1:3) = mp_index(1:3)
       else
          call m_Kp_sample_mesh(n_bz_mesh)
       end if
    end if

    kv3bz = kv3_fbz /nspin
    if(.not.allocated(vkbz)) allocate(vkbz(kv3bz,3));  if(.not.allocated(wbz)) allocate(wbz(kv3bz))

    wbz(1:kv3bz) = product(reduction_factor_exx(1:3))/dble(kv3bz)
    Do ik=1, kv3bz
       vkbz(ik,1:3) = vkxyz_fbz(nspin*(ik-1)+1,1:3,BUCS)
    End do

    write(nfout,'("k-points in IBZ")')
    write(nfout,'("ik, k1, k2, k3")')
    do ik=1,kv3_for_innerloop,nspin
       write(nfout,'(i5,3(1x,f10.5))') (ik-1)/nspin+1, vkxyz_for_innerloop(ik,1:3,BUCS)
    end do

    write(nfout,'("k-points in BZ")')
    write(nfout,'("ikbz k1 k2 k3 wbz")')
    do ikbz=1,kv3bz
       write(nfout,'(i7,4(1x,f20.5))') ikbz, vkbz(ikbz,1:3), wbz(ikbz)
    end do

    if(.not.allocated(ikp)) allocate(ikp(kv3bz))
    if(.not.allocated(iops)) allocate(iops(kv3bz))
    if(.not.allocated(iopt)) allocate(iopt(kv3bz)); iopt = 0
    if(.not.allocated(itrs)) allocate(itrs(kv3bz)); itrs = 0

    ! q(bz) = S*k'(ibz)
    do ikbz=1,kv3bz
       ikp(ikbz) = 0

       LOOP_S: do jtrs=0, 1
          do iopr=1,nopr
             rkxyz = (1-2*jtrs) * matmul(oprec(:,:,iopr),vkbz(ikbz,1:3))
             do jk=1,kv3_for_innerloop,nspin
                dk(1:3) = vkxyz_for_innerloop(jk,1:3,BUCS)-rkxyz(1:3)
                if( abs(dk(1)-nint(dk(1)))<DELTA07 .and. &
                  & abs(dk(2)-nint(dk(2)))<DELTA07 .and. &
                  & abs(dk(3)-nint(dk(3)))<DELTA07 ) then
                   g0(1:3) = nint(dk(1:3))
                   ikp(ikbz) = jk
                   iops(ikbz) = iopr
                   itrs(ikbz) = jtrs
                  exit LOOP_S
                end if
             end do
          end do
       end do LOOP_S
       if(ikp(ikbz) == 0) then
          write(nfout,'("ikbz=",i5)') ikbz
          write(nfout,'(" ikp=",i5)') ikp(ikbz)
          write(nfout,'("rkxyz=",3(1x,f20.5))') rkxyz(1:3)
          call phase_error_with_msg(nfout,"Did not find the index of k' in IBZ.",__LINE__,__FILE__)
       end if

       do iopr=1,nopr
          ss = matmul(op(:,:,iopr),op(:,:,iops(ikbz)))
          if(abs(ss(1,1)-1.d0)<DELTA07 .and. &
           & abs(ss(2,2)-1.d0)<DELTA07 .and. &
           & abs(ss(3,3)-1.d0)<DELTA07 .and. &
           & abs(ss(1,2))<DELTA07 .and. &
           & abs(ss(1,3))<DELTA07 .and. &
           & abs(ss(2,3))<DELTA07) then
             iopt(ikbz) = iopr
             exit
          end if
       end do

       !   k'(bz) is transferable to k in IBZ by a symmetry operation.
       !   g0 = k(ibz) - R k(bz)
       !    0 = k(ibz) - R k'(bz)
       !    k'(bz) = k(bz) + R^inv g0

       rkxyz = (1-2*itrs(ikbz)) * matmul(oprec(:,:,iopt(ikbz)),g0(1:3))
       vkbz(ikbz,1:3) = vkbz(ikbz,1:3) + rkxyz(1:3)
    end do

    ntrs = maxval(itrs)

    write(nfout,'("k-points in BZ")')
    write(nfout,'("ikbz k1 k2 k3 wbz")')
    do ikbz=1,kv3bz
       write(nfout,'(i7,4(1x,f20.5))') ikbz, vkbz(ikbz,1:3), wbz(ikbz)
    end do

    write(nfout,'(2x,"ikbz",3x,"ikp",2x,"iops",2x,"iopt",2x,"itrs")')
    do ikbz=1,kv3bz
       write(nfout,'(i5,4(1x,i5))') ikbz, ikp(ikbz), iops(ikbz), iopt(ikbz), itrs(ikbz)
    end do
  end subroutine m_ES_EXX_kbz2

  subroutine m_ES_EXX_fbz_faces
    implicit none
    integer :: n
    integer :: i,j,k
    real(kind=DP) :: mtr(3,3)
    integer, parameter :: nface_max = 125
    logical :: exists_in_fbz(nface_max)
    integer :: iface(3,nface_max)
    real(kind=DP), parameter :: lambda = 1.d0+1.d-8
    real(kind=DP) :: gg, gg_ref

    do i=1,3
       do j=1,3
          mtr(i,j) = dot_product(rltv(:,i),rltv(:,j))
       end do
    end do

    n=0
    do i=-2,2
       do j=-2,2
          do k=-2,2
             if(i == 0 .and. j == 0 .and. k == 0) cycle
             n=n+1
             iface(1,n)=i
             iface(2,n)=j
             iface(3,n)=k
          end do
       end do
    end do

    exists_in_fbz(1:nface_max) = .true.
    do i=1,n
       if(exists_in_fbz(i)) then
          gg_ref = dot_product(iface(:,i),matmul(mtr,iface(:,i)))
          do j=1,n
             if(i /= j) then
                gg = dot_product(iface(:,i),matmul(mtr,iface(:,j)))
                if(gg*lambda > gg_ref) exists_in_fbz(j) = .false.
             end if
          end do
       end if
    end do

    nface = 0
    do i=1,n
       if(exists_in_fbz(i)) nface = nface + 1
    end do
    allocate(iface_bz(3,nface))
    nface = 0
    do i=1,n
       if(exists_in_fbz(i)) then
          nface = nface + 1
          iface_bz(1:3,nface) = iface(1:3,i)
       end if
    end do

    ! Write FBZ faces
    write(nfout,'("FBZ faces: iface_bz(3,nface)")')
    write(nfout,'("nface=",i5)') nface
    write(nfout,'("i, n1, n2, n3")')
    do i=1,nface
       write(nfout,'(i5,3(1x,i5))') i, iface_bz(1:3,i)
    end do
  end subroutine m_ES_EXX_fbz_faces

  subroutine m_ES_EXX_move_k_into_fbz(ki,ko)
    implicit none
    real(kind=DP), intent(in) :: ki(3)
    real(kind=DP), intent(out) :: ko(3)

    integer :: i, j
    real(kind=DP) :: mtr(3,3)
    real(kind=DP), parameter :: lambda = 1.d0+1.d-8
    real(kind=DP) :: gg, gg_ref

    do i=1,3
       do j=1,3
          mtr(i,j) = dot_product(rltv(:,i),rltv(:,j))
       end do
    end do

    ko(1:3) = ki(1:3)
    do i=1,nface
       gg_ref = dot_product(iface_bz(:,i),matmul(mtr,iface_bz(:,i))) * 0.5d0
       do
          gg = dot_product(iface_bz(:,i),matmul(mtr,ko))
          if(gg > gg_ref*lambda) then
             ko(1:3) = ko(1:3) - dble(iface_bz(1:3,i))
          else
             exit
          end if
       end do
    end do
  end subroutine m_ES_EXX_move_k_into_fbz

  subroutine m_ES_EXX_init0
    integer :: ii,ie
    if(potential_update>0) then
       if(.not.allocated(exx_potential)) &
       & allocate(exx_potential(maxval(np_g1k),np_e,ista_k:iend_k,kimg))
       exx_potential = 0.d0
    endif
    if(sw_change_axis==ON.and.modnrm==EXECUT)then
       if(.not.allocated(ngabc_red)) allocate(ngabc_red(nmax_G_hyb,3))
       ngabc_red=0

       ie=iend_kngp
       if(iend_kngp>nmax_G_hyb) ie=nmax_G_hyb
       do ii=ista_kngp,ie
          ngabc_red(ii,1:3) = ngabc_kngp_l(ii,1:3)
       enddo
       call mpi_allreduce(MPI_IN_PLACE,ngabc_red,nmax_G_hyb*3,mpi_integer,mpi_sum,mpi_ke_world,ierr)
    endif
  end subroutine m_ES_EXX_init0

  subroutine m_ES_EXX_init
    implicit none
    integer :: ik, ikbz, ii,ierr,ie
    real(kind=DP) :: dk(3), v0(3)
    integer,save  :: id_sname = -1
    call tstatc0_begin('m_ES_EXX_init ',id_sname,level=1)

    if(icond == INITIAL .or. icond == CONTINUATION)then
        kg1_for_innerloop   => kg1
        kv3_for_innerloop   => kv3
        qwgt_for_innerloop  => qwgt
        vkxyz_for_innerloop => vkxyz
        iba_for_innerloop   => iba
        nbase_for_innerloop => nbase
        if (allocated(nbase_gamma)) nbase_gamma_for_innerloop => nbase_gamma
        kv3_for_outerloop   => kv3
        qwgt_for_outerloop  => qwgt
        vkxyz_for_outerloop => vkxyz
    else if (icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION)then
        call m_ES_EXX_rd_hybridinfo()
        if(fixed_charge_k_parallel == ONE_BY_ONE)then
          kv3_for_outerloop   => kv3_ek
          qwgt_for_outerloop  => qwgt_ek
          vkxyz_for_outerloop => vkxyz_ek
        else if (fixed_charge_k_parallel == ALL_AT_ONCE)then
          kv3_for_outerloop   => kv3
          qwgt_for_outerloop  => qwgt
          vkxyz_for_outerloop => vkxyz
        else
          call phase_error_with_msg(nfout,'invalid fixed_charge_k_parallel',__LINE__,__FILE__)
        endif
    else
        call phase_error_with_msg(nfout,'invalid icond',__LINE__,__FILE__)
    endif

! === ASMS ===
!    call m_ES_EXX_kbz ! -> ntrs, etc.
    call m_ES_EXX_kbz2    ! -> ntrs, etc.
! === ASMS ===

! === ASMS ===
!    ikgamma = 1
!    do ik=1,kv3,nspin
!       if( abs(vkxyz(ik,1,BUCS))<DELTA07 .and. &
!         & abs(vkxyz(ik,2,BUCS))<DELTA07 .and. &
!         & abs(vkxyz(ik,3,BUCS))<DELTA07 ) then
!          ikgamma = ik
!          exit
!       end if
!    end do

    kshift_norm2 = kshift(1)**2 +kshift(2)**2 +kshift(3)**2
    ikcenter = 1
    v0(1) = dble(kshift(1)) / dble(n_bz_mesh(1))
    v0(2) = dble(kshift(2)) / dble(n_bz_mesh(2))
    v0(3) = dble(kshift(3)) / dble(n_bz_mesh(3))

    do ik=1,kv3_for_outerloop,nspin
       if( abs(vkxyz_for_outerloop(ik,1,BUCS)-v0(1))<DELTA07 .and. &
            & abs(vkxyz_for_outerloop(ik,2,BUCS)-v0(2))<DELTA07 .and. &
            & abs(vkxyz_for_outerloop(ik,3,BUCS)-v0(3))<DELTA07 ) then
          ikcenter = ik
          exit
       end if
    end do
! === ASMS ===

    if(allocated(qmk_zero)) deallocate(qmk_zero)
    allocate(qmk_zero(kv3bz,kv3_for_outerloop)); qmk_zero = .false.
    if(allocated(qmk)) deallocate(qmk)
    allocate(qmk(kv3bz*kv3_for_outerloop,3))
    if(allocated(qmk_zero2)) deallocate(qmk_zero2)
    allocate(qmk_zero2(kv3bz*kv3_for_outerloop))
    if(allocated(iqmk)) deallocate(iqmk)
    allocate(iqmk(kv3bz,kv3_for_outerloop))
    nqmk = 0
    do ik=1,kv3_for_outerloop
       do ikbz=1,kv3bz
          if(.not.q_on_k_centered_mesh(ikbz,ik)) cycle
          if( abs(vkbz(ikbz,1)-vkxyz_for_outerloop(ik,1,BUCS))<DELTA07 .and. &
            & abs(vkbz(ikbz,2)-vkxyz_for_outerloop(ik,2,BUCS))<DELTA07 .and. &
            & abs(vkbz(ikbz,3)-vkxyz_for_outerloop(ik,3,BUCS))<DELTA07 ) then
             qmk_zero(ikbz,ik) = .true.
          end if
          if(nqmk==0) then
             nqmk = nqmk + 1
             qmk(nqmk,1:3) = vkbz(ikbz,1:3)-vkxyz_for_outerloop(ik,1:3,BUCS)
             iqmk(ikbz,ik) = nqmk
             qmk_zero2(nqmk) = qmk_zero(ikbz,ik)
          else
             dk(1:3) = vkbz(ikbz,1:3)-vkxyz_for_outerloop(ik,1:3,BUCS)
             do ii=1,nqmk
                if( abs(qmk(ii,1)-dk(1))<DELTA07 .and. &
                  & abs(qmk(ii,2)-dk(2))<DELTA07 .and. &
                  & abs(qmk(ii,3)-dk(3))<DELTA07 ) then
                  iqmk(ikbz,ik) = ii
                  go to 100
                end if
             end do
             nqmk = nqmk + 1
             qmk(nqmk,1:3) = dk(1:3)
             iqmk(ikbz,ik) = nqmk
             qmk_zero2(nqmk) = qmk_zero(ikbz,ik)
         100 continue
          end if
       end do
    end do

    if(sw_hybrid_functional==ON) call m_PP_qitgft_qmk()

#ifdef FFTW3
    nfftwf = product(fft_box_size_EXX(1:3,1))
#else
    nfftwf = product(fft_box_size_WF(1:3,1))
    nfft_exx = product(fft_box_size_WF(1:3,0)) * kimg
#endif

    call tstatc0_end(id_sname)
  end subroutine m_ES_EXX_init

  subroutine m_ES_EXX_crotylm
    implicit none
    integer :: l1max, mmax, nsph
    integer :: ia1,it1,il1,im1,tau1,lmt1,lmtt1,isph1
    integer :: ia2,it2,il2,im2,tau2,lmt2,lmtt2,isph2
    integer :: iy,ilmta1,iopr
    integer :: iopr1,iopr2
    integer :: invop(nopr)
    integer :: i,j
    real(kind=DP) :: ss(3,3)
    real(kind=DP) :: x(3), y(3)
    integer,save  :: id_sname = -1
    call tstatc0_begin('m_ES_EXX_crotylm ',id_sname,level=1)

    do iopr1=1,nopr
       do iopr2=1,nopr
          ss = matmul(op(:,:,iopr1),op(:,:,iopr2))
          if(abs(ss(1,1)-1.d0)<DELTA07 .and. &
          & abs(ss(2,2)-1.d0)<DELTA07 .and. &
          & abs(ss(3,3)-1.d0)<DELTA07 .and. &
          & abs(ss(1,2))<DELTA07 .and. &
          & abs(ss(1,3))<DELTA07 .and. &
          & abs(ss(2,3))<DELTA07) then
             invop(iopr1) = iopr2
             exit
          end if
       end do
    end do

    l1max = nloc
    mmax  = 2*l1max-1
    nsph = l1max**2
    if(.not.allocated(crotylm)) allocate(crotylm(mmax,nsph,nopr))
    crotylm = 0.0d0
    if(.not.allocated(iylm))    allocate(iylm(mmax,nsph,nopr))
    if(.not.allocated(nylm))    allocate(nylm(nsph,nopr))
    call get_crotylm(l1max,mmax,nsph,nopr,crotylm,iylm,nylm,op)

    if(.not.allocated(isph_lmta)) allocate(isph_lmta(nlmta))
    if(.not.allocated(lmta_rot))  allocate(lmta_rot(mmax,nlmta,nopr))
    lmta_rot = 0
    if(.not.allocated(ia_lmta))   allocate(ia_lmta(nlmta))
    do ia1=1,natm
       it1=ityp(ia1)
       do lmt1=1,ilmt(it1)

          il1  = ltp( lmt1, it1)
          im1  = mtp( lmt1, it1)
          tau1 = taup(lmt1, it1)
          lmtt1 = lmtt(lmt1,it1)
          isph1 = (il1-1)**2 + im1

          ilmta1 = lmta(lmt1,ia1)
          ia_lmta(ilmta1) = ia1
          isph_lmta(ilmta1) = isph1
          do iopr=1,nopr
             ia2=napt(ia1,invop(iopr))
             it2=ityp(ia2)
             if(it1/=it2) call phase_error_with_msg(nfout,'m_ES_EXX_crotylm: it1 and it2 differ.',__LINE__,__FILE__)
             do lmt2=1,ilmt(it2)
                il2  = ltp( lmt2, it2)
                tau2 = taup(lmt2, it2)
                if(il1==il2.and.tau1==tau2) then
                   im2  = mtp( lmt2, it2)
                   lmtt2 = lmtt(lmt2,it2)
                   isph2 = (il2-1)**2 + im2
                   do iy=1,nylm(isph1,iopr)
                      if(isph2 == iylm(iy,isph1,iopr)) then
                         lmta_rot(iy,ilmta1,iopr) = lmta(lmt2,ia2)
                      end if
                   end do
                end if
             end do
          end do
       end do
    end do

    if(.not.allocated(opdir)) allocate(opdir(3,3,nopr))
    do iopr=1,nopr
       do j=1,3
          x = matmul(op(:,:,iopr),altv(:,j))
          do i=1,3
             y = rltv(:,i)
             opdir(i,j,iopr) = dot_product(y,x)/PAI2
          end do
       end do
    end do

    call tstatc0_end(id_sname)
  end subroutine m_ES_EXX_crotylm

  subroutine m_ES_EXX_dealloc
    implicit none
    deallocate(vc)
    deallocate(wfv)
    deallocate(occup_val)
    deallocate(ngpt_exx)
    if(modnrm == EXECUT) then
       if(sw_rspace_hyb==OFF) then
         deallocate(ylm_exx)
         deallocate(qitg_exx)
       endif
       deallocate(fsr_exx)
       deallocate(fsi_exx)
    end if
  end subroutine m_ES_EXX_dealloc

  subroutine m_ES_EXX_gather_valence_states(nfout,transform)
    implicit none
    integer, intent(in) :: nfout
    logical, optional, intent(in) :: transform

    integer :: ik,ib,ib1,ibm,irev
    real(kind=DP), allocatable :: wfv_mpi(:,:,:,:) ! d(kg1,nval,kv3,kimg)
    real(kind=DP), allocatable :: occup_val_mpi(:,:) ! d(nval,kv3)
    real(kind=DP), allocatable :: fsr_mpi(:,:,:) ! d(nval,nlmta,kv3)
    real(kind=DP), allocatable :: fsi_mpi(:,:,:) ! d(nval,nlmta,kv3)
    real(kind=DP), allocatable :: efsr_l(:,:) ! d(np_e,nlmta)
    real(kind=DP), allocatable :: efsi_l(:,:) ! d(np_e,nlmta)

    real(kind=DP), allocatable :: zaj_buf(:,:),zaj_buf2(:,:),wfvv(:,:,:,:)
    integer :: ig,iadd
    integer :: ip0,ip1

    logical :: trans
    integer, allocatable, dimension(:) :: ista
    integer,save  :: id_sname = -1

    if(sw_update_wfv==OFF) return
    if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) return

    allocate(ista(MPI_STATUS_SIZE))

    trans = .true.
    if(present(transform)) trans = transform
    call tstatc0_begin('m_ES_EXX_gather_valence_states ',id_sname,level=1)



    ibm = 0
    do ik=1,kv3,af+1
       if(map_k(ik) /= myrank_k) cycle
       do ib=1,neg
          ib1 = neordr(ib,ik)
          if(map_e(ib1) == myrank_e) then
             if(occup_l(map_z(ib1),ik) > DELTA) then
                ibm = max(ib,ibm)
             end if
          end if
       end do
    end do
    if(npes>1) then
       call mpi_allreduce(ibm,irev,1,MPI_INTEGER,MPI_MAX,mpi_kg_world,ierr)
! === DEBUG by tkato 2014/04/16 ================================================
!      nval = irev
       ibm = irev
! ==============================================================================
       call mpi_allreduce(ibm,irev,1,MPI_INTEGER,MPI_MAX,mpi_ge_world,ierr)
       nval = irev
    else
       nval = ibm
    end if

    if(sw_change_axis==ON)then
       call m_Parallel_init_mpi_nval(nfout,ipri,printable,nval)
    endif

    if(allocated(wfv)) deallocate(wfv)
    if(allocated(occup_val)) deallocate(occup_val)
    if(sw_change_axis==ON)then
       allocate(wfv(kg1,mp_nval,kv3,kimg));wfv=0.d0
       allocate(wfvv(kg1,mp_e,kv3,kimg));wfvv=0.d0
    else
       allocate(wfv(maxval(np_g1k),nval,kv3,kimg))
    endif
    allocate(occup_val(nval,kv3))
    wfv = 0.d0
    occup_val = 0.d0
    if(modnrm == EXECUT) then
       if(allocated(fsr_exx)) deallocate(fsr_exx)
       if(allocated(fsi_exx)) deallocate(fsi_exx)
       allocate(fsr_exx(nval,nlmta,kv3))
       allocate(fsi_exx(nval,nlmta,kv3))
       allocate(efsr_l(np_e,nlmta));efsr_l=0.d0
       allocate(efsi_l(np_e,nlmta));efsi_l=0.d0
       fsr_exx = 0.d0
       fsi_exx = 0.d0
    end if
    do ik=1,kv3,af+1
       if(map_k(ik) /= myrank_k) cycle
       if(modnrm == EXECUT) call get_expkt_fs(ik,fsr_l,fsi_l,efsr_l,efsi_l)
       do ib=1,nval
          ib1 = neordr(ib,ik)
          if(map_e(ib1) == myrank_e) then
             if(sw_change_axis==ON)then
                do ig=ista_g1k(ik),iend_g1k(ik)
                   iadd = ig-ista_g1k(ik)+1
                   wfvv(ig,map_z(ib),ik,1:kimg) = zaj_l(iadd,map_z(ib),ik,1:kimg)
                enddo
             else
                wfv(:,ib,ik,1:kimg) = zaj_l(:,map_z(ib1),ik,1:kimg)
             endif
             occup_val(ib,ik) = occup_l(map_z(ib1),ik)
             if(modnrm == EXECUT) then
                fsr_exx(ib,1:nlmta,ik) = efsr_l(map_z(ib1),1:nlmta)
                fsi_exx(ib,1:nlmta,ik) = efsi_l(map_z(ib1),1:nlmta)
             end if
          end if
       end do
    end do
    if(modnrm == EXECUT) then
       deallocate(efsr_l)
       deallocate(efsi_l)
    end if

    allocate(occup_val_mpi(nval,kv3))
    if(sw_change_axis/=ON)then
       allocate(wfv_mpi(maxval(np_g1k),nval,kv3,kimg))
       call mpi_allreduce(wfv,wfv_mpi,maxval(np_g1k)*nval*kv3*kimg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
       wfv = wfv_mpi
       call mpi_allreduce(wfv,wfv_mpi,maxval(np_g1k)*nval*kv3*kimg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ge_world,ierr)
       wfv = wfv_mpi
    else
       allocate(wfv_mpi(kg1,mp_e,kv3,kimg));wfv_mpi=0.d0
       call mpi_allreduce(wfvv,wfv_mpi,kg1*mp_e*kv3*kimg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
       wfvv = wfv_mpi
       call mpi_allreduce(wfvv,wfv_mpi,kg1*mp_e*kv3*kimg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ge_world,ierr)
       wfvv = wfv_mpi
       do ik=1,kv3
          do ib=1,nval
             ip0=map_e(ib)
             ip1=map_nval(ib)
             if(ip0==myrank_e .and. ip1==myrank_nval) &
            & wfv(1:kg1,map_z_nval(ib),ik,1:kimg) = wfvv(1:kg1,map_z(ib),ik,1:kimg)
          enddo
       enddo
       call mpi_allreduce(MPI_IN_PLACE,wfv,kg1*mp_nval*kv3*kimg,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
       deallocate(wfvv)
    endif
    call mpi_allreduce(occup_val,occup_val_mpi,nval*kv3,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
    occup_val = occup_val_mpi
    call mpi_allreduce(occup_val,occup_val_mpi,nval*kv3,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ge_world,ierr)
    occup_val = occup_val_mpi
    deallocate(wfv_mpi)
    deallocate(occup_val_mpi)
    if(modnrm == EXECUT) then
       allocate(fsr_mpi(nval,nlmta,kv3))
       allocate(fsi_mpi(nval,nlmta,kv3))
       call mpi_allreduce(fsr_exx,fsr_mpi,nval*nlmta*kv3,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
       fsr_exx = fsr_mpi
       call mpi_allreduce(fsr_exx,fsr_mpi,nval*nlmta*kv3,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ge_world,ierr)
       call mpi_allreduce(fsi_exx,fsi_mpi,nval*nlmta*kv3,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
       fsi_exx = fsi_mpi
       call mpi_allreduce(fsi_exx,fsi_mpi,nval*nlmta*kv3,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ge_world,ierr)
       fsr_exx = fsr_mpi
       fsi_exx = fsi_mpi
       deallocate(fsr_mpi)
       deallocate(fsi_mpi)
    end if

    if(sw_rspace_hyb==ON .and. modnrm==EXECUT.and.sw_precalculate==ON)then
       if(allocated(fsrqm)) deallocate(fsrqm)
       if(allocated(fsiqm)) deallocate(fsiqm)
       allocate(fsrqm(nmesh_rs_aug_max,nlmta,nval,kv3bz,nspin));fsrqm=0.d0
       allocate(fsiqm(nmesh_rs_aug_max,nlmta,nval,kv3bz,nspin));fsiqm=0.d0
       call build_fsqm()
    endif


    deallocate(ista)
    call tstatc0_end(id_sname)
  end subroutine m_ES_EXX_gather_valence_states

  subroutine build_fsqm()
    integer :: jkbz,kk,kop,jop,jtrs,m
    integer :: ia,it,imesh,lmt1,lmt2,nma,ilmta1,ilmta2,ispin,lmtp
    real(kind=DP), allocatable, dimension(:) :: fsr, fsi !d(nlmta)
    real(kind=DP), allocatable, dimension(:) :: cosqmkr,sinqmkr
    real(kind=DP) :: fr,fi,fr0,fi0
    real(kind=DP) :: qm
    allocate(fsr(nlmta));fsr=0.d0
    allocate(fsi(nlmta));fsi=0.d0
    fsrqm = 0.d0
    fsiqm = 0.d0
    do ispin=1,nspin,af+1
       do jkbz=1,kv3bz
          ! q(bz) = S*k(ibz)
          kk = ikp(jkbz) + ispin - 1
          kop = iops(jkbz)
          jop = iopt(jkbz)
          jtrs = itrs(jkbz)
          do m=1,nval
             call get_Rot_betar_dot_WFs(kk,m,jop,kop,jtrs,fsr,fsi)
             do ia=1,natm
                it = ityp(ia)
                if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
                nma = nmesh_rs_aug(ia)
                do lmtp=1,nlmtpair(ia)
                   lmt1 = plmt1(lmtp,ia)
                   lmt2 = plmt2(lmtp,ia)
                   ilmta1 = lmta(lmt1,ia)
                   ilmta2 = lmta(lmt2,ia)
                   fr = fsr(ilmta2)
                   fi = fsi(ilmta2)
                   fr0 = fsr(ilmta1)
                   fi0 = fsi(ilmta1)
                   do imesh=1,nma
                      qm = qr_clm_ylm(imesh,ia,lmtp)
                      fsrqm(imesh,ilmta1,m,jkbz,ispin) = fsrqm(imesh,ilmta1,m,jkbz,ispin)+qm*fr
                      fsiqm(imesh,ilmta1,m,jkbz,ispin) = fsiqm(imesh,ilmta1,m,jkbz,ispin)+qm*fi
                      if(lmt1.ne.lmt2)then
                         fsrqm(imesh,ilmta2,m,jkbz,ispin) = fsrqm(imesh,ilmta2,m,jkbz,ispin)+qm*fr0
                         fsiqm(imesh,ilmta2,m,jkbz,ispin) = fsiqm(imesh,ilmta2,m,jkbz,ispin)+qm*fi0
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    deallocate(fsr)
    deallocate(fsi)
  end subroutine build_fsqm

  subroutine m_ES_EXX_kernel(nfout)
    implicit none
    integer, intent(in) :: nfout

    integer :: ik,ikbz,ig,kgs
    real(kind=DP) :: fac, wi, kg(3), vzero, g2
    real(kind=DP), dimension(6) :: ttr
    integer :: igs,ige
    integer,save  :: id_sname = -1
    call tstatc0_begin('m_ES_EXX_kernel ',id_sname,level=1)

    if(sw_change_axis==ON)then
    if(.not.allocated(vc)) allocate(vc(nmax_G_hyb,nqmk))
    else
    if(.not.allocated(vc)) allocate(vc(ista_kngp:iend_kngp,nqmk))
    endif

    vc = 0.d0

    call getttr(rltv,ttr)

    fac = PAI4/univol

    if(sw_singular_correction==ON) then
       call singular_correction(vzero)
    else
       vzero = 0.d0
    end if
    vzero = vzero*dble(kv3bz)/product(reduction_factor_exx(1:3))

    if ( hybrid_functional_type == 'gaupbe' ) then
       call case_gaupbe
    else
       if(sw_screened_exchange==OFF) then
          do ik=1,nqmk
             kgs=1
             if(qmk_zero2(ik)) then
                kgs=2
                if(ista_kngp == 1) vc(1,ik) = vzero
             end if
             igs = max(kgs,ista_kngp);ige=iend_kngp
             if(ige.gt.nmax_G_hyb) ige=nmax_G_hyb
             do ig=igs,ige
                kg(1:3) = qmk(ik,1:3) + ngabc_kngp_l(ig,1:3)
                g2          = ttr(1)*kg(1)*kg(1) &
                     &           + ttr(2)*kg(2)*kg(2) &
                     &           + ttr(3)*kg(3)*kg(3) &
                     &           + ttr(4)*kg(1)*kg(2) &
                     &           + ttr(5)*kg(2)*kg(3) &
                     &           + ttr(6)*kg(3)*kg(1)
                if(g2<1e-12) then
                   vc(ig,ik) = vzero
                else
                   vc(ig,ik) = fac / g2
                endif
             end do
          end do
       else
          wi = 1.d0/(4.d0*omega_exx*omega_exx)
          do ik=1,nqmk
             kgs=1
             if(qmk_zero2(ik)) then
                kgs=2
                if(ista_kngp == 1) vc(1,ik) = vzero
             end if
             igs = max(kgs,ista_kngp);ige=iend_kngp
             if(ige.gt.nmax_G_hyb) ige=nmax_G_hyb
             do ig=igs,ige
                kg(1:3) = qmk(ik,1:3) + ngabc_kngp_l(ig,1:3)
                g2          = ttr(1)*kg(1)*kg(1) &
                     &           + ttr(2)*kg(2)*kg(2) &
                     &           + ttr(3)*kg(3)*kg(3) &
                     &           + ttr(4)*kg(1)*kg(2) &
                     &           + ttr(5)*kg(2)*kg(3) &
                     &           + ttr(6)*kg(3)*kg(1)
                if(g2<1e-12) then
                   vc(ig,ik) = vzero
                else
                   vc(ig,ik) = fac / g2 * (1.d0-exp(-g2*wi))
                endif
             end do
          end do
       end if
       if(sw_change_axis == ON) &
            & call mpi_allreduce(MPI_IN_PLACE,vc,nmax_G_hyb*nqmk,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)

       vc = -alpha_exx * vc
    endif

    call tstatc0_end(id_sname)

    !!stop 'm_ES_EXX_kernel: Singular correction'

  contains

    subroutine case_gaupbe        ! J. Chem. Phys. 138 (2013) 241101.
      integer :: ik, ig, igs, ige
      real(kind=DP) :: fac, c1

      real(kind=DP), parameter :: alpha_gaupbe = 0.15d0
      real(kind=DP), parameter :: beta_gaupbe  = 0.24d0

      fac = ( PAI /alpha_gaupbe )**(3.d0/2.d0) /univol *beta_gaupbe

      do ik=1,nqmk
         kgs=1
         igs = max(kgs,ista_kngp);ige=iend_kngp
         if(ige.gt.nmax_G_hyb) ige=nmax_G_hyb
         do ig=igs,ige
            kg(1:3) = qmk(ik,1:3) + ngabc_kngp_l(ig,1:3)
            g2          = ttr(1)*kg(1)*kg(1) &
                 &           + ttr(2)*kg(2)*kg(2) &
                 &           + ttr(3)*kg(3)*kg(3) &
                 &           + ttr(4)*kg(1)*kg(2) &
                 &           + ttr(5)*kg(2)*kg(3) &
                 &           + ttr(6)*kg(3)*kg(1)
            c1 = exp( -g2 /4.0d0 /alpha_gaupbe )
            vc(ig,ik) = -fac *c1
         end do
      end do
      if(sw_change_axis == ON) then
         call mpi_allreduce(MPI_IN_PLACE,vc,nmax_G_hyb*nqmk,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
      endif
    end subroutine case_gaupbe

    subroutine singular_correction(chi)
      implicit none
      real(kind=DP), intent(out) :: chi

      real(kind=DP), parameter :: gam  = 1.d0 ! Bohr^-2
      real(kind=DP) :: chig1, chig2, chig3, gam3
      real(kind=DP) :: chi0, chi1, chi2

      !chig1 = chig(gam1)
      !chig2 = chig(gam2)
      !chi_app = chig1 + (chig1-chig2)/(gam2-gam1) * gam1

      call chig_deriv(gam,chi0,chi1,chi2)
      chi = chi0 - chi1*gam + 0.5d0*chi2*gam**2

      if(sw_screened_exchange==ON) then
         gam3 = 1.d0/(4.d0*omega_exx**2)
         chig3 = chig(gam3)
         chi = chi - chig3
      !   chi_app = chi_app - chig3
      end if

      if(printable) then
         write(nfout,'("EXX: Singular correction: chi(0)=",f10.5)') chi
         write(nfout,'("EXX: chi0(",f10.5,") =",f10.5)') gam, chi0
         write(nfout,'("EXX: chi1(",f10.5,") =",f10.5)') gam, chi1
         write(nfout,'("EXX: chi2(",f10.5,") =",f10.5)') gam, chi2
         if(sw_screened_exchange==ON) &
          & write(nfout,'("EXX: chi0(",f10.5,") =",f10.5)') gam3, chig3
      end if
    end subroutine singular_correction

    function chig(gam)
      implicit none
      real(kind=DP) :: chig
      real(kind=DP), intent(in) :: gam

      integer :: ig, ikbz, kgs, igs, ige
      real(kind=DP) :: sumg, sumgk, v0(3)

      sumg = 0.d0
      do ikbz=1,kv3bz
! ==== ASMS ===
!         if(.not.q_on_k_centered_mesh(ikbz,ikgamma)) cycle
!         kgs=1
!         if(qmk_zero(ikbz,ikgamma)) kgs=2
!
         if(.not.q_on_k_centered_mesh(ikbz,ikcenter)) cycle
         kgs=1
#if 1
         if ( qmk_zero(ikbz,ikcenter) ) kgs=2
         v0(:) = vkxyz_for_outerloop(ikcenter,:,BUCS)
#else
         if ( kshift_norm2 < DELTA07 .and. qmk_zero(ikbz,ikcenter) ) kgs=2
         v0(:) = 0.0d0
#endif
! ==== ASMS ===

         sumgk = 0.d0
         do ig=max(kgs,ista_kngp),iend_kngp
            kg(1:3) = vkbz(ikbz,1:3) -v0(1:3) + ngabc_kngp_l(ig,1:3)
            g2          = ttr(1)*kg(1)*kg(1) &
            &           + ttr(2)*kg(2)*kg(2) &
            &           + ttr(3)*kg(3)*kg(3) &
            &           + ttr(4)*kg(1)*kg(2) &
            &           + ttr(5)*kg(2)*kg(3) &
            &           + ttr(6)*kg(3)*kg(1)
            sumgk = sumgk + exp(-gam*g2) / g2
         end do
         sumg = sumg + wbz(ikbz) * sumgk
      end do
      sumg = fac * sumg
      call mpi_allreduce(MPI_IN_PLACE,sumg,1,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)

      chig = 1.d0/sqrt(PAI*gam) - sumg
    end function chig

    subroutine chig_deriv(gam,chi0,chi1,chi2)
      implicit none
      real(kind=DP), intent(in) :: gam
      real(kind=DP), intent(out) :: chi0, chi1, chi2

      integer :: ig, ikbz, kgs
      real(kind=DP) :: sumg0, sumgk0
      real(kind=DP) :: sumg1, sumgk1
      real(kind=DP) :: sumg2, sumgk2
      real(kind=DP) :: expg2, v0(3)

      sumg0 = 0.d0
      sumg1 = 0.d0
      sumg2 = 0.d0
      do ikbz=1,kv3bz
! ==== ASMS ===
!         if(.not.q_on_k_centered_mesh(ikbz,ikgamma)) cycle
!         kgs=1
!         if(qmk_zero(ikbz,ikgamma)) kgs=2
!
         if(.not.q_on_k_centered_mesh(ikbz,ikcenter)) cycle
         kgs=1
#if 1
         if ( qmk_zero(ikbz,ikcenter) ) kgs=2
         v0(:) = vkxyz_for_outerloop(ikcenter,:,BUCS)
#else
         if ( kshift_norm2 < DELTA07 .and. qmk_zero(ikbz,ikcenter) ) kgs=2
         v0(:) = 0.0d0
#endif
! ==== ASMS ===
         sumgk0 = 0.d0
         sumgk1 = 0.d0
         sumgk2 = 0.d0
         do ig=max(kgs,ista_kngp),iend_kngp
            kg(1:3) = vkbz(ikbz,1:3) -v0(1:3) + ngabc_kngp_l(ig,1:3)
            g2          = ttr(1)*kg(1)*kg(1) &
            &           + ttr(2)*kg(2)*kg(2) &
            &           + ttr(3)*kg(3)*kg(3) &
            &           + ttr(4)*kg(1)*kg(2) &
            &           + ttr(5)*kg(2)*kg(3) &
            &           + ttr(6)*kg(3)*kg(1)
            expg2 = exp(-gam*g2)
            sumgk0 = sumgk0 + expg2 / g2
            sumgk1 = sumgk1 + expg2
            sumgk2 = sumgk2 + expg2 * g2
         end do
         sumg0 = sumg0 + wbz(ikbz) * sumgk0
         sumg1 = sumg1 + wbz(ikbz) * sumgk1
         sumg2 = sumg2 + wbz(ikbz) * sumgk2
      end do
      call mpi_allreduce(MPI_IN_PLACE,sumg0,1,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
      call mpi_allreduce(MPI_IN_PLACE,sumg1,1,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
      call mpi_allreduce(MPI_IN_PLACE,sumg2,1,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
      chi0 =  1.d0/sqrt(PAI*gam) - fac*sumg0
      chi1 = -0.5d0/sqrt(PAI*gam**3) + fac*sumg1
      chi2 =  0.75d0/sqrt(PAI*gam**5) - fac*sumg2
    end subroutine chig_deriv

  end subroutine m_ES_EXX_kernel

  subroutine m_ES_EXX_potential(nfout,ispin,ib,ik,nff,fsr,fsi,vxw,iupdate,store_potential,exx,exx_only)
    implicit none
    integer, intent(in) :: nfout
    integer, intent(in) :: ispin,ib,ik
    integer, intent(in) :: nff
    integer, intent(in), optional :: iupdate
    logical, intent(in), optional :: store_potential
    integer :: iup
    real(kind=DP), intent(in)  :: fsr(nlmta),fsi(nlmta)
    real(kind=DP), intent(out) :: vxw(nff,kimg)
    real(kind=DP), intent(out), optional :: exx
    logical, intent(in), optional :: exx_only

    real(kind=DP), allocatable :: efsr_l(:)
    real(kind=DP), allocatable :: efsi_l(:)
    logical :: eo
    real(kind=DP) :: ene
    integer :: ig,iadd
    logical :: store_p
    real(kind=DP), allocatable, dimension(:) :: zajbuf_r,zajbuf_i
    integer :: kgw,kgv
    integer :: id_sname=-1
    call tstatc0_begin('m_ES_EXX_potential ',id_sname,level=1)
    store_p = .true.
    if(present(store_potential)) store_p = store_potential

    eo = .false.
    if(present(exx_only)) eo = exx_only
    vxw = 0.0d0

    iup = 2
    if(present(iupdate)) iup = iupdate
    if(iup.lt.potential_update)then
       vxw(1:np_g1k(ik),1:kimg) = exx_potential(1:np_g1k(ik),ib,ik,1:kimg)
       if(present(exx))then
         exx=0.d0
         if(kimg==1)then
            do ig = ista_g1k(ik), iend_g1k(ik)
               iadd = ig - ista_g1k(ik) + 1
               exx = exx + zaj_l(iadd,ib,ik,1)*vxw(iadd,1)
            enddo
         else
            if(k_symmetry(ik) == GAMMA) then
               do ig = max(2,ista_g1k(ik)), iend_g1k(ik)
                  iadd = ig - ista_g1k(ik) + 1
                  exx = exx + zaj_l(iadd,ib,ik,1) * vxw(iadd,1) &
                      &     + zaj_l(iadd,ib,ik,2) * vxw(iadd,2)
               end do
! === DEBUG by tkato 2014/04/14 ================================================
!            if(ista_g1k(ik) == 1) eko = eko * 2.d0 + zaj_l(1,ib,ik,1) * vxw(1,1)
               if(ista_g1k(ik) == 1) then
                  exx = exx * 2.d0 + zaj_l(1,ib,ik,1) * vxw(1,1)
               else
                  exx = exx * 2.d0
               end if
! ==============================================================================
            else
               do ig = ista_g1k(ik), iend_g1k(ik)
                  iadd = ig - ista_g1k(ik) + 1
                  exx = exx + zaj_l(iadd,ib,ik,1)*vxw(iadd,1) &
                      &     + zaj_l(iadd,ib,ik,2)*vxw(iadd,2)
               end do
            endif
         endif
       endif
       return
    endif
    allocate(efsr_l(nlmta));efsr_l=0.d0
    allocate(efsi_l(nlmta));efsi_l=0.d0
    if(modnrm == EXECUT) call get_expkt_fs_b(ik,ib,fsr,fsi,efsr_l,efsi_l)
    if(sw_change_axis==ON)then
       allocate(zajbuf_r(kg1));zajbuf_r=0.d0
       allocate(zajbuf_i(kg1));zajbuf_i=0.0d0
       do ig=ista_g1k(ik),iend_g1k(ik)
          iadd = ig-ista_g1k(ik)+1
          zajbuf_r(ig) = zaj_l(iadd,ib,ik,1)
          zajbuf_i(ig) = zaj_l(iadd,ib,ik,kimg)
       enddo
       call mpi_allreduce(MPI_IN_PLACE,zajbuf_r,kg1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
       call mpi_allreduce(MPI_IN_PLACE,zajbuf_i,kg1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
       kgw = kg1;kgv=maxval(np_g1k)
    else
       allocate(zajbuf_r(maxval(np_g1k)));zajbuf_r(1:np_g1k(ik))=zaj_l(1:np_g1k(ik),ib,ik,1)
       allocate(zajbuf_i(maxval(np_g1k)));zajbuf_i(1:np_g1k(ik))=zaj_l(1:np_g1k(ik),ib,ik,kimg)
       kgw = maxval(np_g1k);kgv=maxval(np_g1k)
    endif

    if(sw_change_axis==ON)then
       if(sw_rspace_hyb == ON .and. sw_rspace_hyb_dgm == ON.and.modnrm==EXECUT) then
          if(present(exx))then
             call apply_Vx_to_WF_2D_rs_dgm( ispin, ib, ik, kgw, kgv, zajbuf_r, zajbuf_i, efsr_l, efsi_l, &
                   &  vxw, ene, eo )
          else
             call apply_Vx_to_WF_2D_rs_dgm( ispin, ib, ik, kgw, kgv, zajbuf_r, zajbuf_i, efsr_l, efsi_l, vxw )
          endif
       else
          if(present(exx))then
             call apply_Vx_to_WF_2D( ispin, ib, ik, kgw, kgv, zajbuf_r, zajbuf_i, efsr_l, efsi_l, vxw, ene, eo )
          else
             call apply_Vx_to_WF_2D( ispin, ib, ik, kgw, kgv, zajbuf_r, zajbuf_i, efsr_l, efsi_l, vxw )
          endif
       endif
    else
    if(present(exx))then
       call apply_Vx_to_WF( ispin, ib, ik, kgw, kgv, zajbuf_r, zajbuf_i, efsr_l, efsi_l, vxw, ene, eo )
    else
       call apply_Vx_to_WF( ispin, ib, ik, kgw, kgv, zajbuf_r, zajbuf_i, efsr_l, efsi_l, vxw )
    endif
    endif

    deallocate(efsr_l)
    deallocate(efsi_l)
    deallocate(zajbuf_r)
    deallocate(zajbuf_i)
    if(present(exx)) exx = ene
    if(.not.eo.and.potential_update>0.and.store_p) &
    &  exx_potential(1:np_g1k(ik),ib,ik,1:kimg) = vxw(1:np_g1k(ik),1:kimg)
    call tstatc0_end(id_sname)
  end subroutine m_ES_EXX_potential

  subroutine m_ES_EXX_eigenvalue_for_each_k(ispin,ik,update_eko,iupdate)
    implicit none
    integer, intent(in) :: ispin, ik
    logical, intent(in), optional :: update_eko
    integer, intent(in), optional :: iupdate

    integer :: ib,ig
    integer :: iup
    real(kind=DP), allocatable, dimension(:,:,:,:) :: vxw !d(kg1,kimg)
    integer :: iadd
    real(kind=DP), allocatable, dimension(:,:,:) :: vexx
    real(kind=DP) :: exx
    integer :: ng
    logical :: u_eko
    integer,save  :: id_sname = -1
    call tstatc0_begin('m_ES_EXX_eigenvalue_for_each_k ',id_sname,level=1)

    iup = 2
    if(present(iupdate)) iup = iupdate

    if(.not.allocated(eexx_kb)) allocate(eexx_kb(np_e,ista_k:iend_k))

    if(.not.allocated(eko_hyb))then
       allocate(eko_hyb(np_e,ista_k:iend_k));eko_hyb=0.d0
    endif
    eko_hyb(:,ik) = 0.0d0
    u_eko = .true.
    if(.not.present(update_eko))then
       eko_l(:,ik) = 0.d0
    else
       if(update_eko) eko_l(:,ik) = 0.d0
       u_eko = update_eko
    endif
    ng=maxval(np_g1k)
    allocate(vxw(ng,np_e,ista_k:iend_k,kimg))
       do ib=1,np_e   ! MPI
          if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
             call m_ES_EXX_potential(nfout,ispin,ib,ik,ng,fsr_l(ib,:,ik),fsi_l(ib,:,ik),&
             & vxw(1:ng,ib,ik,1:kimg),iup,exx=exx,exx_only=.true.,store_potential=.false.)
          else
             call m_ES_EXX_potential(nfout,ispin,ib,ik,ng,fsr_l(ib,:,ik),fsr_l(ib,:,ik),&
             & vxw(1:ng,ib,ik,1:kimg),iup,exx=exx,exx_only=.true.,store_potential=.false.)
          endif
          call mpi_allreduce(MPI_IN_PLACE,exx,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
          eexx_kb(ib,ik) = occup_l(ib,ik)*exx
          eko_hyb(ib,ik) = exx
          if(u_eko) then
             eko_l(ib,ik) = eko_l(ib,ik) + exx
          endif
       end do

    deallocate(vxw)

    call tstatc0_end(id_sname)
  end subroutine m_ES_EXX_eigenvalue_for_each_k

  subroutine m_ES_EXX_cp_eigenvalue(ik)
     integer, intent(in) :: ik
     if(.not.allocated(eko_hyb)) return
     eko_l(1:np_e,ik) = eko_hyb(1:np_e,ik)
  end subroutine m_ES_EXX_cp_eigenvalue

  subroutine m_ES_EXX_Diagonal_part(ispin,ik,vxdi)
    implicit none
    integer, intent(in) :: ispin, ik
    real(kind=DP), intent(out), dimension(maxval(np_g1k)) :: vxdi

    integer :: ib,ig
    integer :: ii,m,jkbz,kk,jk,jop,jtrs
    integer :: iadd
    real(kind=DP) :: fac, occupation
    real(kind=DP), allocatable, dimension(:,:) :: rho, pot, vx, vxq !d(nfft)
    real(kind=DP), allocatable, dimension(:) :: rhor, rhoi !d(kg1)
    real(kind=DP), allocatable, dimension(:) :: wfr, wfi !d(kg1)
    integer :: lsize, ibsize

    real(kind=DP), allocatable, dimension(:,:,:,:) :: vxdi_t

    integer,save  :: id_sname = -1

    if(sw_change_axis==ON)then
       call m_ES_EXX_Diagonal_part_(ispin,ik,vxdi)
       return
    endif

    call tstatc0_begin('m_ES_EXX_Diagonal_part ',id_sname,level=1)

#ifdef FFT_3D_DIVISION
    lsize = fft_X_x_nel*fft_X_y_nel*fft_X_z_nel
#else
    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
#endif
    ibsize = 1

    allocate(rho(lsize*kimg,ibsize))
    allocate(pot(lsize*kimg,ibsize))
    allocate(vx(lsize*kimg,ibsize)); vx = 0.d0
    allocate(vxq(lsize*kimg,ibsize))
    allocate(rhor(ista_kngp:iend_kngp)); rhor = 0.d0
    allocate(rhoi(ista_kngp:iend_kngp)); rhoi = 0.d0
    allocate(wfr(maxval(np_g1k))); wfr = 0.d0
    allocate(wfi(maxval(np_g1k))); wfi = 0.d0

    ! rho(G) = sum_q sum_m |phi_m(q+G)|^2
    do jkbz=1,kv3bz
       if(.not.q_on_k_centered_mesh(jkbz,ik)) cycle
       ! q(bz) = S*k(ibz)
       kk = ikp(jkbz) + ispin - 1
       jop = iopt(jkbz)
       jtrs = itrs(jkbz)
       wfr = 0.0d0
       wfi = 0.0d0
       do m=1,nval
          if(occup_val(m,kk) < DELTA) cycle
          occupation = wbz(jkbz)*occup_val(m,kk)/qwgt_for_innerloop(kk)/dble(kv3_for_innerloop)
          if(kimg==1) then
             do ii = ista_g1k(kk), iend_g1k(kk)
                iadd = ii - ista_g1k(kk) + 1
                wfr(iadd) = wfr(iadd) + occupation * ( wfv(iadd,m,kk,1) * wfv(iadd,m,kk,1) )
             end do
          else
             do ii = ista_g1k(kk), iend_g1k(kk)
                iadd = ii - ista_g1k(kk) + 1
                wfr(iadd) = wfr(iadd) + occupation * ( wfv(iadd,m,kk,1) * wfv(iadd,m,kk,1) &
                                                   & + wfv(iadd,m,kk,2) * wfv(iadd,m,kk,2) )
             end do
          end if
       end do
       call map_Rot_WFG_on_FFT_box_3D(kk,1,1,ibsize,lsize,rho,wfr,wfi,jop,jtrs)

       do ii=ista_kngp,min(kg,iend_kngp)
          rhor(ii) = vc(ii,iqmk(jkbz,ik))
          rhoi(ii) = 0.d0
       end do
       call map_RHOG_on_FFT_box_3D(ik,1,1,ibsize,lsize,pot,rhor,rhoi)

! === FFT xzy ==================================================================
       if(sw_fft_xzy > 0) then
! ==============================================================================
       call product_on_FFT_box_3D(pot,rho,vxq,lsize,ibsize,nel_fft_y(myrank_g))
! === FFT xzy ==================================================================
       else
       call product_on_FFT_box_3D(pot,rho,vxq,lsize,ibsize,nel_fft_z(myrank_g))
       end if
! ==============================================================================
       vx = vx + vxq
    end do

! === FFT xzy ==================================================================
    if(sw_fft_xzy > 0) then
! ==============================================================================
    call m_FFT_Direct_3D (nfout, vx, lsize, ibsize)
! === FFT xzy ==================================================================
    else
    call m_FFT_Direct_XYZ_3D (nfout, vx, lsize, ibsize)
    end if
! ==============================================================================
    call map_FFT_box_on_WFG_3D(ik,lsize,ibsize,vx,1,1,wfr,wfi)
    fac = 1.d0/dble(nfftwf)
    do ii = ista_g1k(ik), iend_g1k(ik)
       iadd = ii - ista_g1k(ik) + 1
       vxdi(iadd) = fac * wfr(iadd)
    end do

    deallocate(rho)
    deallocate(pot)
    deallocate(vx)
    deallocate(rhor)
    deallocate(rhoi)
    deallocate(wfr)
    deallocate(wfi)

    call tstatc0_end(id_sname)
  end subroutine m_ES_EXX_Diagonal_part

  subroutine m_ES_EXX_Diagonal_part_(ispin,ik,vxdi)
    implicit none
    integer, intent(in) :: ispin, ik
    real(kind=DP), intent(out), dimension(maxval(np_g1k)) :: vxdi

    integer :: ib,ig
    integer :: ii,m,jkbz,kk,jk,jop,jtrs,mm
    real(kind=DP) :: fac, occupation
    real(kind=DP), allocatable, dimension(:) :: rho, pot, vx, vxq !d(nfft)
    real(kind=DP), allocatable, dimension(:) :: rhor, rhoi !d(kg1)

    real(kind=DP), allocatable, dimension(:,:,:,:) :: vxdi_t

    integer,save  :: id_sname = -1
    call tstatc0_begin('m_ES_EXX_Diagonal_part ',id_sname,level=1)


    allocate(rho(nfft))
    allocate(pot(nfft))
    allocate(vx(nfft)); vx = 0.d0
    allocate(vxq(nfft))
    allocate(rhor(nmax_G_hyb)); rhor = 0.d0
    allocate(rhoi(nmax_G_hyb)); rhoi = 0.d0

    ! rho(G) = sum_q sum_m |phi_m(q+G)|^2
    do jkbz=1,kv3bz
       if(.not.q_on_k_centered_mesh(jkbz,ik)) cycle
       ! q(bz) = S*k(ibz)
       kk = ikp(jkbz) + ispin - 1
       jop = iopt(jkbz)
       jtrs = itrs(jkbz)
! === DEBUB by tkato 2013/02/21 ================================================
       rhor = 0.0d0
       rhoi = 0.0d0
! ==============================================================================
       do m=ista_nval,iend_nval
          if(occup_val(m,kk) < DELTA) cycle
          mm = map_z_nval(m)
          occupation = wbz(jkbz)*occup_val(m,kk)/qwgt_for_innerloop(kk)/dble(kv3_for_innerloop)
          if(kimg==1) then
             do ii=1,iba_for_innerloop(kk)
                rhor(ii) = rhor(ii) + occupation * ( wfv(ii,mm,kk,1) * wfv(ii,mm,kk,1) )
             end do
          else
             do ii=1, min( iba_for_innerloop(kk), nmax_G_hyb )
                rhor(ii) = rhor(ii) + occupation * ( wfv(ii,mm,kk,1) * wfv(ii,mm,kk,1) &
                                                 & + wfv(ii,mm,kk,2) * wfv(ii,mm,kk,2) )
             end do
          end if
       end do
       call map_Rot_WFG_on_FFT_box(kk,jop,jtrs,rhor,rhoi,rho)
#ifdef FFTW3
       call m_FFT_exx(nfout,rho,INVERSE,OFF) ! rho(R)
#else
       call m_FFT_WF(ELECTRON,nfout,rho,INVERSE,OFF) ! rho(R)
#endif

       do ii=1, min( kg, nmax_G_hyb )
          rhor(ii) = vc(ii,iqmk(jkbz,ik))
          rhoi(ii) = 0.d0
       end do
       call map_RHOG_on_FFT_box(rhor,rhoi,pot)
#ifdef FFTW3
       call m_FFT_exx(nfout,pot,INVERSE,OFF) ! pot(R)
#else
       call m_FFT_WF(ELECTRON,nfout,pot,INVERSE,OFF) ! pot(R)
#endif

       call product_on_FFT_box(pot,rho,vxq) ! vx(R)
       vx = vx + vxq
    end do

#ifdef FFTW3
    call m_FFT_exx(nfout,vx,DIRECT,OFF) ! vx(G)
#else
    call m_FFT_WF(ELECTRON,nfout,vx,DIRECT,OFF) ! vx(G)
#endif
    call map_FFT_box_on_RHOG(rhor,rhoi,vx)
    fac = 1.d0/dble(nfftwf)
    vxdi = 0.0d0
!    do ii=1,iba(ik)
    do ii=ista_g1k(ik),iend_g1k(ik)
       ig = nbase(ii,ik)
       if ( ig <= nmax_G_hyb ) then
          vxdi(ii-ista_g1k(ik)+1) = fac * rhor(ig)
       endif
    end do

    deallocate(rho)
    deallocate(pot)
    deallocate(vx)
    deallocate(rhor)
    deallocate(rhoi)

    call tstatc0_end(id_sname)
  end subroutine m_ES_EXX_Diagonal_part_

  subroutine m_ES_EXX_add_Diagonal_part(ik,ibo,vxdi,vnldi)
    implicit none
    integer, intent(in) :: ik, ibo
    real(kind=DP), intent(in), dimension(maxval(np_g1k)) :: vxdi
    real(kind=DP), intent(inout), dimension(maxval(np_g1k)) :: vnldi

    integer :: ib
    integer :: i, iadd

    do i = ista_g1k(ik), iend_g1k(ik)
       iadd = i - ista_g1k(ik) + 1
       vnldi(iadd) = vnldi(iadd) + vxdi(iadd)
    end do
  end subroutine m_ES_EXX_add_Diagonal_part

  subroutine m_ES_Vexx_W(ik,iupdate,store_exxp)
    implicit none
    integer, intent(in) :: ik
    integer, intent(in), optional :: iupdate
    logical, intent(in), optional :: store_exxp

    logical :: store_e

    integer :: ib,ig,ispin,iup
    real(kind=DP), allocatable, dimension(:,:,:,:) :: vxw !d(kg1,kimg)
    integer :: iadd
    real(kind=DP), allocatable, dimension(:,:,:) :: vexx
    real(kind=DP) :: exx
    integer :: ng
    integer :: npprev
    real(kind=DP), allocatable, dimension(:,:,:,:) :: exx_potential_buf
    integer,save  :: id_sname = -1
    iup = 2
    if(present(iupdate)) iup = iupdate
    store_e = .true.
    if(present(store_exxp)) store_e = store_exxp
    call tstatc0_begin('m_ES_Vexx_W ',id_sname,level=1)
    if(allocated(exx_potential)) then
      npprev = size(exx_potential, 1)
      if (npprev<maxval(np_g1k)) then
        allocate(exx_potential_buf(npprev, np_e, ista_k:iend_k, kimg))
        exx_potential_buf = exx_potential
        deallocate(exx_potential)
        allocate(exx_potential(maxval(np_g1k),np_e,ista_k:iend_k,kimg))
        exx_potential(1:npprev, np_e, ista_k:iend_k, kimg) = exx_potential_buf(1:npprev, np_e, ista_k:iend_k, kimg)
      endif
    endif

    ispin = mod(ik-1,nspin)+1

    ng = maxval(np_g1k)
       allocate(vxw(ng,np_e,ista_k:iend_k,kimg))
       do ib=1,np_e   ! MPI
          if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
             call m_ES_EXX_potential(nfout,ispin,ib,ik,ng,fsr_l(ib,:,ik),fsi_l(ib,:,ik),&
             &    vxw(1:ng,ib,ik,1:kimg),iup,store_potential=store_e)
          else
             call m_ES_EXX_potential(nfout,ispin,ib,ik,ng,fsr_l(ib,:,ik),fsr_l(ib,:,ik),&
             &    vxw(1:ng,ib,ik,1:kimg),iup,store_potential=store_e)
          endif
       end do
       do ib=1,np_e
          call m_ES_Vexx_add_vexx(ik,ib,vxw(1:ng,ib,ik,1:kimg))
       enddo
       deallocate(vxw)
    call tstatc0_end(id_sname)
  end subroutine m_ES_Vexx_W

  subroutine m_ES_Vexx_add_vexx(ik,ib,vxw_exx)
    integer, intent(in) :: ib,ik
    real(kind=DP),dimension(maxval(np_g1k),kimg),intent(in) :: vxw_exx
    integer :: ig,iadd
    logical ma
    if(kimg==1) then
       do ig = ista_g1k(ik), iend_g1k(ik)
          iadd = ig - ista_g1k(ik) + 1
          vnlph_l(iadd,ib,1) = vnlph_l(iadd,ib,1) + vxw_exx(iadd,1)
       end do
    else
       do ig = ista_g1k(ik), iend_g1k(ik)
          iadd = ig - ista_g1k(ik) + 1
          vnlph_l(iadd,ib,1) = vnlph_l(iadd,ib,1) + vxw_exx(iadd,1)
          vnlph_l(iadd,ib,2) = vnlph_l(iadd,ib,2) + vxw_exx(iadd,2)
       end do
    end if
  end subroutine m_ES_Vexx_add_vexx

  subroutine apply_Vx_to_WF(ispin,ib,ik,kgw,kgv,wfr,wfi,bdwr,bdwi,vxw,eexx,eonly,force_l,dbdwr,dbdwi)
    implicit none
    integer, intent(in)                       :: ispin, ib,ik,kgw,kgv
    real(kind=DP), intent(in), dimension(kgw) :: wfr, wfi
    real(kind=DP), intent(in), dimension(nlmta) :: bdwr, bdwi
    real(kind=DP), intent(out), optional, dimension(kgv,kimg) :: vxw
    real(kind=DP), intent(out), optional :: eexx
    logical, intent(in), optional :: eonly

    real(kind=DP), intent(in), dimension(nlmta,3), optional :: dbdwr, dbdwi
    real(kind=DP), intent(inout), dimension(natm,3), optional :: force_l
    logical :: force_mode=.false.
    real(kind=DP) :: exx
    integer :: kk,m,ig,ii,jop,kop,jtrs
    integer :: jkbz
    integer :: iadd
    real(kind=DP) :: fac, ph, occupation, dnorm
    real(kind=DP) :: skg(3)
    real(kind=DP), allocatable, dimension(:,:) :: rho,phi !d(nfft)
    real(kind=DP), allocatable, dimension(:,:) :: wfn,wfm !d(nfft)
    real(kind=DP), allocatable, dimension(:,:) :: sumdel !d(nfft)
    real(kind=DP), allocatable, dimension(:) :: rhogr, rhogi !d(kgp)
    real(kind=DP), allocatable, dimension(:) :: wfsr,wfsi !d(kg1)
    real(kind=DP), allocatable, dimension(:) :: cosgt, singt !d(kg1)
    real(kind=DP), allocatable, dimension(:) :: fsr, fsi !d(nlmta)
    real(kind=DP), allocatable, dimension(:,:) :: dfsr, dfsi !d(nlmta,3)
    real(kind=DP), allocatable, dimension(:) :: qvr, qvi !d(nlmta)
#ifdef HYBRID_DGEMM
    real(kind=DP), allocatable, dimension(:,:,:) :: tqvr, tqvi !d(nlmta,nval,kv3bz)
#endif
    real(kind=DP), allocatable, dimension(:,:) :: dqvr, dqvi !d(nlmta,3)
    real(kind=DP), allocatable, dimension(:,:) :: gqvr, gqvi !d(nlmta,3)
    real(kind=DP), allocatable, dimension(:) :: sumqvr, sumqvi !d(nlmta)
    real(kind=DP), allocatable, dimension(:,:) :: sumdqvr, sumdqvi !d(nlmta)
    real(kind=DP), allocatable, dimension(:) :: rhor,rhoi
    real(kind=DP), allocatable, dimension(:) :: afft
    real(kind=DP) :: ifac
    logical :: eo
    integer :: ikk
    integer,save  :: id_sname = -1,id_sname1=-2,id_sname2=-3,id_sname3=-4
    integer, allocatable, dimension(:,:) :: ngabc_tmp
    integer :: lsize, ibsize
#ifdef HYBRID_DGEMM
    real(kind=DP), allocatable, dimension(:,:) :: cosgt_, singt_
    real(kind=DP), allocatable, dimension(:,:,:) :: rhogr_, rhogi_
    real(kind=DP), allocatable, dimension(:,:,:,:) :: wfm_
    real(kind=DP), allocatable, dimension(:,:,:) :: fsr_, fsi_
    real(kind=DP), allocatable, dimension(:,:,:,:) :: tdqvr, tdqvi !d(nlmta,nval,kv3bz,3)
    real(kind=DP), allocatable, dimension(:,:,:,:) :: tgqvr, tgqvi !d(nlmta,nval,kv3bz,3)
#endif
    call tstatc0_begin('apply_Vx_to_WF ',id_sname)

    ikk = k_index(ik)

    eo=.false.
    if(present(eonly))then
      eo = eonly
    endif
    if(present(eexx))then
      eexx = 0.d0
    endif
#ifdef FFT_3D_DIVISION
    lsize = fft_X_x_nel*fft_X_y_nel*fft_X_z_nel
#else
    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
#endif
    ibsize = 1

    force_mode = present(force_l).and.present(dbdwr).and.present(dbdwi)

    allocate(rho(lsize*kimg,ibsize))
    allocate(phi(lsize*kimg,ibsize))
    allocate(wfn(lsize*kimg,ibsize))
    allocate(wfm(lsize*kimg,ibsize))
    allocate(sumdel(lsize*kimg,ibsize))
    allocate(rhogr(ista_kngp:iend_kngp))
    if(kimg==2) allocate(rhogi(ista_kngp:iend_kngp))
    allocate(wfsr(maxval(np_g1k)))
    if(kimg==2) allocate(wfsi(maxval(np_g1k)))
    allocate(cosgt(maxval(np_g1k)))
    if(kimg==2) allocate(singt(maxval(np_g1k)))
    if(modnrm == EXECUT) then
       allocate(fsr(nlmta))
       allocate(fsi(nlmta))
       allocate(qvr(nlmta))
       allocate(qvi(nlmta))
#ifdef HYBRID_DGEMM
       allocate(tqvr(nlmta,nval,kv3bz))
       allocate(tqvi(nlmta,nval,kv3bz))
#endif
       allocate(sumqvr(nlmta))
       allocate(sumqvi(nlmta))
       fac = 1.d0/dble(nfftwf)
       sumqvr = 0.d0
       sumqvi = 0.d0
       if(force_mode)then
          allocate(dfsr(nlmta,3))
          allocate(dfsi(nlmta,3))
          allocate(dqvr(nlmta,3))
          allocate(dqvi(nlmta,3))
          allocate(gqvr(nlmta,3))
          allocate(gqvi(nlmta,3))
          allocate(sumdqvr(nlmta,3))
          allocate(sumdqvi(nlmta,3))
          sumdqvr=0.d0;sumdqvi=0.d0
       endif
    end if
    allocate(ngabc_tmp(kg1,3))

    call map_WFG_on_FFT_box_3D(ik,1,1,ibsize,lsize,wfn,wfr,wfi)

    sumdel = 0.d0
#ifdef HYBRID_DGEMM
!#########################################################################################
    if(kimg == 2 .and. modnrm == EXECUT .and. sw_fft_xzy <= 0) then
    allocate(cosgt_(maxval(np_g1k),kv3bz))
    allocate(singt_(maxval(np_g1k),kv3bz))
    allocate(rhogr_(ista_kngp:iend_kngp,nval,kv3bz))
    allocate(rhogi_(ista_kngp:iend_kngp,nval,kv3bz))
    allocate(wfm_(lsize*kimg,ibsize,nval,kv3bz))
    allocate(fsr_(nlmta,nval,kv3bz))
    allocate(fsi_(nlmta,nval,kv3bz))

    if(force_mode)then
       allocate(tdqvr(nlmta,nval,kv3bz,3))
       allocate(tdqvi(nlmta,nval,kv3bz,3))
       allocate(tgqvr(nlmta,nval,kv3bz,3))
       allocate(tgqvi(nlmta,nval,kv3bz,3))
    end if

    do jkbz=1,kv3bz
       if(.not.q_on_k_centered_mesh(jkbz,ikk)) cycle
       kk = ikp(jkbz) + ispin - 1
       kop = iops(jkbz)
       jop = iopt(jkbz)
       jtrs = itrs(jkbz)
       ngabc_tmp = 0
       do ig=1,iba_for_innerloop(kk)
          ii = nbase_for_innerloop(ig,kk)
          if(ii >= ista_kngp .and. ii <= iend_kngp) then
             ngabc_tmp(ig,:) = ngabc_kngp_l(ii,:)
          end if
       end do
       call mpi_allreduce(MPI_IN_PLACE,ngabc_tmp,kg1*3,MPI_INTEGER,MPI_SUM,mpi_ke_world,ierr)
       do ig = ista_g1k(kk), iend_g1k(kk)
          iadd = ig - ista_g1k(kk) + 1
          ii = nbase_for_innerloop(ig,kk)
          skg(1:3) = vkxyz_for_innerloop(kk,1:3,BUCS) + ngabc_tmp(ig,1:3)
          ph = PAI2 * dot_product(skg,tau(:,kop,BUCS))
          cosgt_(iadd,jkbz) = cos(ph)
          singt_(iadd,jkbz) = sin(ph)
       end do
    end do

    do jkbz=1,kv3bz
       if(.not.q_on_k_centered_mesh(jkbz,ikk)) cycle
       kk = ikp(jkbz) + ispin - 1
       kop = iops(jkbz)
       jop = iopt(jkbz)
       jtrs = itrs(jkbz)
       do m=1,nval
          if(occup_val(m,kk) < DELTA) cycle
!         occupation = wbz(jkbz) * occup_val(m,kk)/qwgt(kk)/dble(kv3)
          if(jtrs==0) then
             do ig = ista_g1k(kk), iend_g1k(kk)
                iadd = ig - ista_g1k(kk) + 1
                wfsr(iadd) = wfv(iadd,m,kk,1) * cosgt_(iadd,jkbz) - wfv(iadd,m,kk,2) * singt_(iadd,jkbz)
                wfsi(iadd) = wfv(iadd,m,kk,2) * cosgt_(iadd,jkbz) + wfv(iadd,m,kk,1) * singt_(iadd,jkbz)
             end do
          else
             do ig = ista_g1k(kk), iend_g1k(kk)
                iadd = ig - ista_g1k(kk) + 1
                wfsr(iadd) = wfv(iadd,m,kk,1) * cosgt_(iadd,jkbz) - wfv(iadd,m,kk,2) * singt_(iadd,jkbz)
                wfsi(iadd) = -wfv(iadd,m,kk,2) * cosgt_(iadd,jkbz) - wfv(iadd,m,kk,1) * singt_(iadd,jkbz)
             end do
          end if
          call map_Rot_WFG_on_FFT_box_3D(kk,1,1,ibsize,lsize,wfm_(1,1,m,jkbz),wfsr,wfsi,jop,jtrs)
          call product_on_FFT_box_3D(wfn,wfm_(1,1,m,jkbz),rho,lsize,ibsize,nel_fft_z(myrank_g))
          call m_FFT_Direct_XYZ_3D (nfout, rho, lsize, ibsize)
          call map_FFT_box_on_RHOG_3D(ik, lsize, ibsize, rho, 1, 1,   &
                                  rhogr_(ista_kngp,m,jkbz), rhogi_(ista_kngp,m,jkbz))
          do ii=ista_kngp,min(kg,iend_kngp)
             rhogr_(ii,m,jkbz) = fac * rhogr_(ii,m,jkbz)
             rhogi_(ii,m,jkbz) = fac * rhogi_(ii,m,jkbz)
          end do
          do ii=max(kg+1,ista_kngp),iend_kngp
             rhogr_(ii,m,jkbz) = 0.d0
             rhogi_(ii,m,jkbz) = 0.d0
          end do
       end do
    end do

    do jkbz=1,kv3bz
       if(.not.q_on_k_centered_mesh(jkbz,ikk)) cycle
       kk = ikp(jkbz) + ispin - 1
       kop = iops(jkbz)
       jop = iopt(jkbz)
       jtrs = itrs(jkbz)
       do m=1,nval
          if(occup_val(m,kk) < DELTA) cycle
          if(force_mode)then
             call get_Rot_betar_dot_WFs(kk,m,jop,kop,jtrs,fsr_(1,m,jkbz),fsi_(1,m,jkbz),dfsr,dfsi)
          else
             call get_Rot_betar_dot_WFs(kk,m,jop,kop,jtrs,fsr_(1,m,jkbz),fsi_(1,m,jkbz))
          endif
       end do
    end do

    call add_RHOG_hard_part_(iqmk,rhogr_,rhogi_,bdwr,bdwi,fsr_,fsi_,ikk,ispin)

    do jkbz=1,kv3bz
       if(.not.q_on_k_centered_mesh(jkbz,ikk)) cycle
       kk = ikp(jkbz) + ispin - 1
       do m=1,nval
          if(occup_val(m,kk) < DELTA) cycle
          occupation = wbz(jkbz) * occup_val(m,kk)/qwgt_for_innerloop(kk)/dble(kv3_for_innerloop)
          if(present(eexx))then
             exx = 0.d0
             call sum_rho_vc_rho(rhogr_(ista_kngp,m,jkbz),rhogi_(ista_kngp,m,jkbz),   &
                         vc(ista_kngp,iqmk(jkbz,ikk)),exx)
             eexx = eexx + occupation * exx
             if(eo) cycle
          endif
          do ii=ista_kngp,iend_kngp
             rhogr_(ii,m,jkbz) = vc(ii,iqmk(jkbz,ikk)) * rhogr_(ii,m,jkbz) ! phi(G) = Vc(G,q) * rho(G)
             rhogi_(ii,m,jkbz) = vc(ii,iqmk(jkbz,ikk)) * rhogi_(ii,m,jkbz) ! phi(G) = Vc(G,q) * rho(G)
          end do
       end do
    end do

!    do jkbz=1,kv3bz
!       if(.not.q_on_k_centered_mesh(jkbz,ik)) cycle
!       kk = ikp(jkbz) + ispin - 1
!       kop = iops(jkbz)
!       jop = iopt(jkbz)
!       jtrs = itrs(jkbz)
!       do m=1,nval
!          if(occup_val(m,kk) < DELTA) cycle
          if(force_mode)then
             call integrate_QijVnm_(iqmk,  &
                      rhogr_(ista_kngp,1,1),rhogi_(ista_kngp,1,1),  &
                      fsr_, fsi_,   &
                      tqvr,tqvi,ikk,ispin,dfsr,dfsi,tdqvr,tdqvi,tgqvr,tgqvi)
          else
             call integrate_QijVnm_(iqmk,  &
                      rhogr_(ista_kngp,1,1),rhogi_(ista_kngp,1,1),  &
                      fsr_, fsi_,    &
                      tqvr,tqvi,ikk,ispin)
          endif
    do jkbz=1,kv3bz
       if(.not.q_on_k_centered_mesh(jkbz,ikk)) cycle
       kk = ikp(jkbz) + ispin - 1
       if(.not.force_mode)then
          do m=1,nval
             if(occup_val(m,kk) < DELTA) cycle
             occupation = wbz(jkbz) * occup_val(m,kk)/qwgt_for_innerloop(kk)/dble(kv3_for_innerloop)
                call mpi_allreduce(MPI_IN_PLACE,tqvr(1,m,jkbz),nlmta,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
                call mpi_allreduce(MPI_IN_PLACE,tqvi(1,m,jkbz),nlmta,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
                sumqvr = sumqvr + occupation * tqvr(1:nlmta,m,jkbz) ! sum_m qvr_m
                sumqvi = sumqvi + occupation * tqvi(1:nlmta,m,jkbz) ! sum_m qvi_m
                call map_RHOG_on_FFT_box_3D(ik,1,1,ibsize,lsize,phi,  &
                         rhogr_(ista_kngp,m,jkbz),rhogi_(ista_kngp,m,jkbz))
                call product_on_FFT_box_3D(phi,wfm_(1,1,m,jkbz),rho,lsize,ibsize,nel_fft_z(myrank_g))
                sumdel = sumdel + occupation * rho ! sum delta(R)
          end do
       else
          do m=1,nval
             if(occup_val(m,kk) < DELTA) cycle
             occupation = wbz(jkbz) * occup_val(m,kk)/qwgt_for_innerloop(kk)/dble(kv3_for_innerloop)
             sumqvr = sumqvr + occupation * tqvr(1:nlmta,m,jkbz) ! sum_m qvr_m
             sumqvi = sumqvi + occupation * tqvi(1:nlmta,m,jkbz) ! sum_m qvi_m
             sumdqvr = sumdqvr + occupation * (tdqvr(1:nlmta,m,jkbz,1:3)+tgqvr(1:nlmta,m,jkbz,1:3)) ! sum_m (dqvr_m + gqvr_m)
             sumdqvi = sumdqvi + occupation * (tdqvi(1:nlmta,m,jkbz,1:3)+tgqvi(1:nlmta,m,jkbz,1:3)) ! sum_m (dqvi_m + gqvi_m)
          end do
       endif
    end do

    deallocate(cosgt_)
    deallocate(singt_)
    deallocate(rhogr_)
    deallocate(rhogi_)
    deallocate(wfm_)
    deallocate(fsr_)
    deallocate(fsi_)

    if(force_mode)then
       deallocate(tdqvr)
       deallocate(tdqvi)
       deallocate(tgqvr)
       deallocate(tgqvi)
    endif

    else
!#########################################################################################
#endif
    do jkbz=1,kv3bz
       if(.not.q_on_k_centered_mesh(jkbz,ikk)) cycle
       ! q(bz) = S*k(ibz)
       kk = ikp(jkbz) + ispin - 1
       kop = iops(jkbz)
       jop = iopt(jkbz)
       jtrs = itrs(jkbz)
       ngabc_tmp = 0
       do ig=1,iba_for_innerloop(kk)
          ii = nbase_for_innerloop(ig,kk)
          if(ii >= ista_kngp .and. ii <= iend_kngp) then
             ngabc_tmp(ig,:) = ngabc_kngp_l(ii,:)
          end if
       end do
       call mpi_allreduce(MPI_IN_PLACE,ngabc_tmp,kg1*3,MPI_INTEGER,MPI_SUM,mpi_ke_world,ierr)
       if(kimg==1) then
          do ig = ista_g1k(kk), iend_g1k(kk)
             iadd = ig - ista_g1k(kk) + 1
             ii = nbase_for_innerloop(ig,kk)
             skg(1:3) = vkxyz_for_innerloop(kk,1:3,BUCS) + ngabc_tmp(ig,1:3)
             ph = PAI2 * dot_product(skg,tau(:,kop,BUCS))
             cosgt(iadd) = cos(ph)
          end do
       else
          do ig = ista_g1k(kk), iend_g1k(kk)
             iadd = ig - ista_g1k(kk) + 1
             ii = nbase_for_innerloop(ig,kk)
             skg(1:3) = vkxyz_for_innerloop(kk,1:3,BUCS) + ngabc_tmp(ig,1:3)
             ph = PAI2 * dot_product(skg,tau(:,kop,BUCS))
             cosgt(iadd) = cos(ph)
             singt(iadd) = sin(ph)
          end do
       end if
       do m=1,nval
          if(occup_val(m,kk) < DELTA) cycle
          occupation = wbz(jkbz) * occup_val(m,kk)/qwgt_for_innerloop(kk)/dble(kv3_for_innerloop)
          if(kimg==1) then
             dnorm = 0.d0
             do ig = ista_g1k(kk), iend_g1k(kk)
                iadd = ig - ista_g1k(kk) + 1
                wfsr(iadd) = wfv(iadd,m,kk,1) * cosgt(iadd)
                dnorm = dnorm + wfsr(iadd)*wfsr(iadd)
             end do
             call mpi_allreduce(MPI_IN_PLACE,dnorm,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
             dnorm = 1.d0/sqrt(dnorm)
             wfsr = wfsr * dnorm
             call map_Rot_WFG_on_FFT_box_3D(kk,1,1,ibsize,lsize,wfm,wfsr,wfsr,jop,jtrs)
          else
             if(jtrs==0) then
                do ig = ista_g1k(kk), iend_g1k(kk)
                   iadd = ig - ista_g1k(kk) + 1
                   wfsr(iadd) = wfv(iadd,m,kk,1) * cosgt(iadd) - wfv(iadd,m,kk,2) * singt(iadd)
                   wfsi(iadd) = wfv(iadd,m,kk,2) * cosgt(iadd) + wfv(iadd,m,kk,1) * singt(iadd)
                end do
             else
                do ig = ista_g1k(kk), iend_g1k(kk)
                   iadd = ig - ista_g1k(kk) + 1
                   wfsr(iadd) = wfv(iadd,m,kk,1) * cosgt(iadd) - wfv(iadd,m,kk,2) * singt(iadd)
                   wfsi(iadd) = -wfv(iadd,m,kk,2) * cosgt(iadd) - wfv(iadd,m,kk,1) * singt(iadd)
                end do
             end if
             call map_Rot_WFG_on_FFT_box_3D(kk,1,1,ibsize,lsize,wfm,wfsr,wfsi,jop,jtrs)
          end if
! === FFT xzy ==================================================================
          if(sw_fft_xzy > 0) then
! ==============================================================================
          call product_on_FFT_box_3D(wfn,wfm,rho,lsize,ibsize,nel_fft_y(myrank_g))
          call m_FFT_Direct_3D (nfout, rho, lsize, ibsize)
! === FFT xzy ==================================================================
          else
          call product_on_FFT_box_3D(wfn,wfm,rho,lsize,ibsize,nel_fft_z(myrank_g))
          call m_FFT_Direct_XYZ_3D (nfout, rho, lsize, ibsize)
          end if
! ==============================================================================

          if(kimg==1) then
             call map_FFT_box_on_RHOG_3D(ikk,lsize, ibsize, rho, 1, 1, rhogr, rhogi)
             do ii=ista_kngp,iend_kngp
                rhogr(ii) = vc(ii,iqmk(jkbz,ikk)) * rhogr(ii) ! phi(G) = Vc(G,q) * rho(G)
             end do
             call map_RHOG_on_FFT_box_3D(ik,1,1,ibsize,lsize,phi,rhogr,rhogr)
          else
             call map_FFT_box_on_RHOG_3D(ikk, lsize, ibsize, rho, 1, 1, rhogr, rhogi)
             !! rho(G) = rho_soft(G) + rho_hard(G)
             if(modnrm == EXECUT) then
                do ii=ista_kngp,min(kg,iend_kngp)
                   rhogr(ii) = fac * rhogr(ii)
                   rhogi(ii) = fac * rhogi(ii)
                end do
                do ii=max(kg+1,ista_kngp),iend_kngp
                   rhogr(ii) = 0.d0
                   rhogi(ii) = 0.d0
                end do
                if(force_mode)then
                   call get_Rot_betar_dot_WFs(kk,m,jop,kop,jtrs,fsr,fsi,dfsr,dfsi)
                else
                   call get_Rot_betar_dot_WFs(kk,m,jop,kop,jtrs,fsr,fsi)
                endif
                call add_RHOG_hard_part(iqmk(jkbz,ikk),rhogr,rhogi,bdwr,bdwi,fsr,fsi)
             end if
             if(present(eexx))then
                exx = 0.d0
                call sum_rho_vc_rho(rhogr,rhogi,vc(ista_kngp,iqmk(jkbz,ikk)),exx)
                eexx = eexx + occupation * exx
                if(eo) cycle
             endif
             do ii=ista_kngp,iend_kngp
                rhogr(ii) = vc(ii,iqmk(jkbz,ikk)) * rhogr(ii) ! phi(G) = Vc(G,q) * rho(G)
                rhogi(ii) = vc(ii,iqmk(jkbz,ikk)) * rhogi(ii) ! phi(G) = Vc(G,q) * rho(G)
             end do
             !! Dij[Vnm] = Sum_G Qij(G) * Vnm(G)
             !! Hi[Vnm] = Sum_j Dij[Vnm] <beta_j|psi_m>
             if (modnrm == EXECUT) then
                if(force_mode)then
                   call integrate_QijVnm(iqmk(jkbz,ikk),rhogr,rhogi,fsr,fsi,qvr,qvi,dfsr,dfsi,dqvr,dqvi,gqvr,gqvi)
                else
                   call integrate_QijVnm(iqmk(jkbz,ikk),rhogr,rhogi,fsr,fsi,qvr,qvi)
                endif
             endif
! =============================================================== 12.5Exp
             if(.not.force_mode)then
! === DEBUG by tkato 2014/04/11 ================================================
                if(modnrm == EXECUT) then
! ==============================================================================
                   call mpi_allreduce(MPI_IN_PLACE,qvr,nlmta,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
                   call mpi_allreduce(MPI_IN_PLACE,qvi,nlmta,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
! === DEBUG by tkato 2014/04/11 ================================================
                end if
! ==============================================================================
                call map_RHOG_on_FFT_box_3D(ik,1,1,ibsize,lsize,phi,rhogr,rhogi)
             endif
          end if

          if(force_mode)then
             sumqvr = sumqvr + occupation * qvr ! sum_m qvr_m
             sumqvi = sumqvi + occupation * qvi ! sum_m qvi_m
             sumdqvr = sumdqvr + occupation * (dqvr+gqvr) ! sum_m (dqvr_m + gqvr_m)
             sumdqvi = sumdqvi + occupation * (dqvi+gqvi) ! sum_m (dqvi_m + gqvi_m)
          else
! === FFT xzy ==================================================================
             if(sw_fft_xzy > 0) then
! ==============================================================================
                call product_on_FFT_box_3D(phi,wfm,rho,lsize,ibsize,nel_fft_y(myrank_g))
! === FFT xzy ==================================================================
             else
                call product_on_FFT_box_3D(phi,wfm,rho,lsize,ibsize,nel_fft_z(myrank_g))
             end if
! ==============================================================================
             sumdel = sumdel + occupation * rho ! sum delta(R)
             !! Bin = Sum_m fm * Hi[Vnm]
             if(modnrm == EXECUT) then
                sumqvr = sumqvr + occupation * qvr ! sum_m qvr_m
                sumqvi = sumqvi + occupation * qvi ! sum_m qvi_m
             end if
          endif
       end do
    end do
#ifdef HYBRID_DGEMM
    end if
#endif
    if(force_mode)then
       !! dEXX/dR(n) = Sum_i Re[ <psi_n|dbeta_i/dR> * Bin +  <psi_n|beta_i> * Cin ]
       call sum_EXX_force_terms(force_l,bdwr,bdwi,dbdwr,dbdwi,sumqvr,sumqvi,sumdqvr,sumdqvi)
    else
! === FFT xzy ==================================================================
       if(sw_fft_xzy > 0) then
! ==============================================================================
          call m_FFT_Direct_3D (nfout, sumdel, lsize, ibsize)
! === FFT xzy ==================================================================
       else
          call m_FFT_Direct_XYZ_3D (nfout, sumdel, lsize, ibsize)
       end if
! ==============================================================================
       if(kimg==1) then
          call map_FFT_box_on_WFG_3D(ik,lsize,ibsize,sumdel,1,1,vxw(1,1),vxw(1,1))
       else
          call map_FFT_box_on_WFG_3D(ik,lsize,ibsize,sumdel,1,1,vxw(1,1),vxw(1,2))
       end if

       if(modnrm == EXECUT)  then
          fac = 1.d0/dble(nfftwf)
       else
          fac = 1.d0/dble(nfftwf)**2
          if(present(eexx)) eexx = eexx*fac
       end if
       vxw = fac*vxw
       if(modnrm == EXECUT) then
          call add_Vx_hard_part(ik,vxw,sumqvr,sumqvi)
       endif
    endif
    deallocate(rho)
    deallocate(phi)
    deallocate(wfn)
    deallocate(wfm)
    deallocate(sumdel)
    deallocate(rhogr)
    if(kimg==2) deallocate(rhogi)
    deallocate(wfsr)
    if(kimg==2) deallocate(wfsi)
    deallocate(cosgt)
    if(kimg==2) deallocate(singt)
    if(modnrm == EXECUT) then
       deallocate(fsr)
       deallocate(fsi)
       deallocate(qvr)
       deallocate(qvi)
       deallocate(sumqvr)
       deallocate(sumqvi)
       if(force_mode)then
          deallocate(dfsr)
          deallocate(dfsi)
          deallocate(dqvr)
          deallocate(dqvi)
          deallocate(gqvr)
          deallocate(gqvi)
          deallocate(sumdqvr)
          deallocate(sumdqvi)
       endif
    end if
    deallocate(ngabc_tmp)


    call tstatc0_end(id_sname)
  end subroutine apply_Vx_to_WF

  subroutine apply_Vx_to_WF_2D(ispin,ib,ik,kgw,kgv,wfr,wfi,bdwr,bdwi,vxw,eexx,eonly,force_l,dbdwr,dbdwi)
    implicit none
    integer, intent(in)                       :: ispin, ib,ik,kgw,kgv
    real(kind=DP), intent(in), dimension(kgw) :: wfr, wfi
    real(kind=DP), intent(in), dimension(nlmta) :: bdwr, bdwi
    real(kind=DP), intent(out), optional, dimension(kgv,kimg) :: vxw
    real(kind=DP), intent(out), optional :: eexx
    logical, intent(in), optional :: eonly
    real(kind=DP), allocatable, dimension(:,:) :: vxw_t
    real(kind=DP), intent(in), dimension(nlmta,3), optional :: dbdwr, dbdwi
    real(kind=DP), intent(inout), dimension(natm,3), optional :: force_l
    logical :: force_mode=.false.
    real(kind=DP) :: exx
    integer :: kk,m,ig,ii,jop,kop,jtrs,mm,ib1
    integer :: jkbz
    real(kind=DP) :: fac, ph, occupation, dnorm
    real(kind=DP) :: skg(3)
    real(kind=DP), allocatable, dimension(:) :: rho,phi !d(nfft)
    real(kind=DP), allocatable, dimension(:) :: wfn,wfm !d(nfft)
    real(kind=DP), allocatable, dimension(:) :: sumdel !d(nfft)
    real(kind=DP), allocatable, dimension(:) :: rhogr, rhogi !d(kgp)
    real(kind=DP), allocatable, dimension(:) :: wfsr,wfsi !d(kg1)
    real(kind=DP), allocatable, dimension(:) :: cosgt, singt !d(kg1)
    real(kind=DP), allocatable, dimension(:) :: fsr, fsi !d(nlmta)
    real(kind=DP), allocatable, dimension(:,:) :: dfsr, dfsi !d(nlmta,3)
    real(kind=DP), allocatable, dimension(:) :: qvr, qvi !d(nlmta)
    real(kind=DP), allocatable, dimension(:,:) :: dqvr, dqvi !d(nlmta,3)
    real(kind=DP), allocatable, dimension(:,:) :: gqvr, gqvi !d(nlmta,3)
    real(kind=DP), allocatable, dimension(:) :: sumqvr, sumqvi !d(nlmta)
    real(kind=DP), allocatable, dimension(:,:) :: sumdqvr, sumdqvi !d(nlmta)
    real(kind=DP), allocatable, dimension(:) :: rhor,rhoi
    real(kind=DP), allocatable, dimension(:) :: afft
    real(kind=DP) :: ifac
    logical :: eo
    integer :: iadd
    integer :: ikk
    integer,save  :: id_sname = -1,id_sname1=-2,id_sname2=-3,id_sname3=-4
    integer,save  :: id_sname_cdfft=-1
    call tstatc0_begin('apply_Vx_to_WF ',id_sname)
    eo=.false.
    if(present(eonly))then
      eo = eonly
    endif
    if(present(eexx))then
      eexx = 0.d0
    endif

    force_mode = present(force_l).and.present(dbdwr).and.present(dbdwi)

    allocate(rho(nfft));rho=0.d0
    allocate(phi(nfft));phi=0.d0
    allocate(wfn(nfft));wfn=0.d0
    allocate(wfm(nfft));wfm=0.d0
    allocate(sumdel(nfft));sumdel=0.d0
    allocate(rhogr(nmax_G_hyb));rhogr=0.d0
    if(sw_rspace_hyb==ON.and.modnrm==EXECUT)then
       allocate(rhor(nfftp_exx_nonpara/2));rhor=0.d0
       allocate(rhoi(nfftp_exx_nonpara/2));rhoi=0.d0
       allocate(afft(nfftp_exx_nonpara));afft=0.d0
    endif
    if(kimg==2) then
        allocate(rhogi(nmax_G_hyb))
        rhogi = 0.d0
    endif
    allocate(wfsr(kg1));wfsr=0.d0
    if(kimg==2) then
    allocate(wfsi(kg1));wfsi=0.d0
    endif
    allocate(cosgt(kg1))
    if(kimg==2) allocate(singt(kg1))
    allocate(vxw_t(kg1,kimg));vxw_t=0.d0
    if(modnrm == EXECUT) then
       allocate(fsr(nlmta))
       allocate(fsi(nlmta))
       allocate(qvr(nlmta))
       allocate(qvi(nlmta))
       allocate(sumqvr(nlmta))
       allocate(sumqvi(nlmta))
       fac = 1.d0/dble(nfftwf)
       sumqvr = 0.d0
       sumqvi = 0.d0
       if(force_mode)then
          allocate(dfsr(nlmta,3))
          allocate(dfsi(nlmta,3))
          allocate(dqvr(nlmta,3))
          allocate(dqvi(nlmta,3))
          allocate(gqvr(nlmta,3))
          allocate(gqvi(nlmta,3))
          allocate(sumdqvr(nlmta,3))
          allocate(sumdqvi(nlmta,3))
          sumdqvr=0.d0;sumdqvi=0.d0
       endif
    end if

    call map_WFG_on_FFT_box(ik,wfr,wfi,wfn)
#ifdef FFTW3
    call m_FFT_exx(nfout,wfn,INVERSE,OFF) ! wfn(R)
#else
    call m_FFT_WF(ELECTRON,nfout,wfn,INVERSE,OFF) ! wfn(R)
#endif

    ikk = k_index(ik)
    sumdel = 0.d0
    do jkbz=1,kv3bz
       if(.not.q_on_k_centered_mesh(jkbz,ikk)) cycle
#ifdef MEMORY_SAVE_EXX
       call m_ES_EXX_ylm_each_k(iqmk(jkbz,ikk))
#ifdef MEMORY_SAVE_MORE_EXX
       call qitgft_qmk_each_k(iqmk(jkbz,ikk))
#endif
#endif
       ! q(bz) = S*k(ibz)
       kk = ikp(jkbz) + ispin - 1
       kop = iops(jkbz)
       jop = iopt(jkbz)
       jtrs = itrs(jkbz)
       if(kimg==1) then
          do ig=1,iba_for_innerloop(kk)
             ii = nbase_for_innerloop(ig,kk)
             skg(1:3) = vkxyz_for_innerloop(kk,1:3,BUCS) + ngabc(ii,1:3)
             ph = PAI2 * dot_product(skg,tau(:,kop,BUCS))
             cosgt(ig) = cos(ph)
          end do
       else
          do ig=1,iba_for_innerloop(kk)
             ii = nbase_for_innerloop(ig,kk)
             skg(1:3) = vkxyz_for_innerloop(kk,1:3,BUCS) + ngabc(ii,1:3)
             ph = PAI2 * dot_product(skg,tau(:,kop,BUCS))
             cosgt(ig) = cos(ph)
             singt(ig) = sin(ph)
          end do
       end if
       do m=ista_nval,iend_nval
          if(occup_val(m,kk) < DELTA) cycle
!          if(sw_rsb==ON)then
!             if(.not.overlap(ib,m))cycle
!          endif
          mm = map_z_nval(m)
          occupation = wbz(jkbz) * occup_val(m,kk)/qwgt_for_innerloop(kk)/dble(kv3_for_innerloop)
          if(kimg==1) then
             dnorm = 0.d0
             do ig=1,iba_for_innerloop(kk)
                wfsr(ig) = wfv(ig,mm,kk,1) * cosgt(ig)
                dnorm = dnorm + wfsr(ig)*wfsr(ig)
             end do
             dnorm = 1.d0/sqrt(dnorm)
             wfsr = wfsr * dnorm
             call map_Rot_WFG_on_FFT_box(kk,jop,jtrs,wfsr,wfsr,wfm)
          else
             if(jtrs==0) then
                do ig=1,iba_for_innerloop(kk)
                   wfsr(ig) = wfv(ig,mm,kk,1) * cosgt(ig) - wfv(ig,mm,kk,2) * singt(ig)
                   wfsi(ig) = wfv(ig,mm,kk,2) * cosgt(ig) + wfv(ig,mm,kk,1) * singt(ig)
                end do
             else
                do ig=1,iba_for_innerloop(kk)
                   wfsr(ig) = wfv(ig,mm,kk,1) * cosgt(ig) - wfv(ig,mm,kk,2) * singt(ig)
                   wfsi(ig) = -wfv(ig,mm,kk,2) * cosgt(ig) - wfv(ig,mm,kk,1) * singt(ig)
                end do
             end if
             call map_Rot_WFG_on_FFT_box(kk,jop,jtrs,wfsr,wfsi,wfm)
          end if
#ifdef FFTW3
          call m_FFT_exx(nfout,wfm,INVERSE,OFF) ! wfm(R)
#else
          call m_FFT_WF(ELECTRON,nfout,wfm,INVERSE,OFF) ! wfm(R)
#endif
          call product_on_FFT_box(wfn,wfm,rho) ! rho_soft(R)
#ifdef FFTW3
          call m_FFT_exx(nfout,rho,DIRECT,OFF) ! rho_soft(G)
#else
          call m_FFT_WF(ELECTRON,nfout,rho,DIRECT,OFF) ! rho_soft(G)
#endif

          if(kimg==1) then
             call map_FFT_box_on_RHOG(rhogr,rhogr,rho)
             do ii=1,nmax_G_hyb
                rhogr(ii) = vc(ii,iqmk(jkbz,ikk)) * rhogr(ii) ! phi(G) = Vc(G,q) * rho(G)
             end do
             call map_RHOG_on_FFT_box(rhogr,rhogr,phi)
          else
             call map_FFT_box_on_RHOG(rhogr,rhogi,rho)
             !! rho(G) = rho_soft(G) + rho_hard(G)
             if(modnrm == EXECUT) then
                do ii=1, min( kg, nmax_G_hyb )
                   rhogr(ii) = fac * rhogr(ii)
                   rhogi(ii) = fac * rhogi(ii)
                end do
                do ii=kg+1, nmax_G_hyb
                   rhogr(ii) = 0.d0
                   rhogi(ii) = 0.d0
                end do
                if(force_mode)then
                   call get_Rot_betar_dot_WFs(kk,m,jop,kop,jtrs,fsr,fsi,dfsr,dfsi)
                else
                   if(sw_rspace_hyb==OFF .or. (sw_rspace_hyb==ON .and. sw_precalculate==OFF)) &
                  &   call get_Rot_betar_dot_WFs(kk,m,jop,kop,jtrs,fsr,fsi)
                endif
                if(sw_rspace_hyb==ON)then
                   if(sw_precalculate==OFF)then
                      call add_RHOG_hard_part_rs(iqmk(jkbz,ikk),rhor,rhoi,bdwr,bdwi,fsr,fsi)
                   else
                      call add_RHOG_hard_part_rs2(iqmk(jkbz,ikk),rhor,rhoi,&
                     & fsrqm(:,:,m,jkbz,ispin),fsiqm(:,:,m,jkbz,ispin),bdwr,bdwi)
                   endif
                   call map_RHOG_on_FFT_box_hard(rhor,rhoi,afft)
                   call m_FFT_CD0_exx(nfout,afft,DIRECT)
                   call map_FFT_box_on_RHOG_hard(rhor,rhoi,afft)
                   do ii=1,nmax_G_hyb
                      rhogr(ii) = rhogr(ii) + rhor(ii)
                      rhogi(ii) = rhogi(ii) + rhoi(ii)
                   end do
                else
                   call add_RHOG_hard_part_2D(iqmk(jkbz,ikk),rhogr,rhogi,bdwr,bdwi,fsr,fsi)
                endif
             end if
             if(present(eexx))then
                exx = 0.d0
                call sum_rho_vc_rho_2D(rhogr,rhogi,vc(1,iqmk(jkbz,ikk)),exx)
                eexx = eexx + occupation * exx
                if(eo) cycle
             endif
             do ii=1,nmax_G_hyb
                rhogr(ii) = vc(ii,iqmk(jkbz,ikk)) * rhogr(ii) ! phi(G) = Vc(G,q) * rho(G)
                rhogi(ii) = vc(ii,iqmk(jkbz,ikk)) * rhogi(ii) ! phi(G) = Vc(G,q) * rho(G)
             end do
             !! Dij[Vnm] = Sum_G Qij(G) * Vnm(G)
             !! Hi[Vnm] = Sum_j Dij[Vnm] <beta_j|psi_m>
!
             if (modnrm == EXECUT) then
                if(sw_rspace_hyb==ON)then
                   call map_RHOG_on_FFT_box_hard_inv(rhogr,rhogi,afft)
                   call m_FFT_CD0_exx(nfout,afft,INVERSE)
                   call map_FFT_box_on_RHOG_hard_inv(rhor,rhoi,afft)
                   if(force_mode)then
                      call integrate_QijVnm_rs(iqmk(jkbz,ikk),rhor,rhoi,fsr,fsi,qvr,qvi,dfsr,dfsi,dqvr,dqvi,gqvr,gqvi)
                   else
                   if(sw_precalculate==OFF)then
                      call integrate_QijVnm_rs(iqmk(jkbz,ikk),rhor,rhoi,fsr,fsi,qvr,qvi)
                   else
                      call integrate_QijVnm_rs2(iqmk(jkbz,ikk),rhor,rhoi, &
                    & fsrqm(:,:,m,jkbz,ispin),fsiqm(:,:,m,jkbz,ispin),qvr,qvi)
                   endif
                   endif
                else
                   if(force_mode)then
                      call integrate_QijVnm_2D(iqmk(jkbz,ikk),rhogr,rhogi,fsr,fsi,qvr,qvi,dfsr,dfsi,dqvr,dqvi,gqvr,gqvi)
                   else
                      call integrate_QijVnm_2D(iqmk(jkbz,ikk),rhogr,rhogi,fsr,fsi,qvr,qvi)
                   endif
                endif
             endif
! =============================================================== 12.5Exp
             if(.not.force_mode)then
                call map_RHOG_on_FFT_box(rhogr,rhogi,phi)
             endif
          end if

          if(force_mode)then
             sumqvr = sumqvr + occupation * qvr ! sum_m qvr_m
             sumqvi = sumqvi + occupation * qvi ! sum_m qvi_m
             sumdqvr = sumdqvr + occupation * (dqvr+gqvr) ! sum_m (dqvr_m + gqvr_m)
             sumdqvi = sumdqvi + occupation * (dqvi+gqvi) ! sum_m (dqvi_m + gqvi_m)
          else
#ifdef FFTW3
             call m_FFT_exx(nfout,phi,INVERSE,OFF) ! phi(R)
#else
             call m_FFT_WF(ELECTRON,nfout,phi,INVERSE,OFF) ! phi(R)
#endif
             call product_on_FFT_box(phi,wfm,rho) ! delta(R)
             sumdel = sumdel + occupation * rho ! sum delta(R)
             !! Bin = Sum_m fm * Hi[Vnm]
             if(modnrm == EXECUT) then
                sumqvr = sumqvr + occupation * qvr ! sum_m qvr_m
                sumqvi = sumqvi + occupation * qvi ! sum_m qvi_m
             end if
          endif
       end do
    end do
    if(force_mode)then
       !! dEXX/dR(n) = Sum_i Re[ <psi_n|dbeta_i/dR> * Bin +  <psi_n|beta_i> * Cin ]
       call mpi_allreduce(MPI_IN_PLACE,sumqvr,nlmta,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       call mpi_allreduce(MPI_IN_PLACE,sumqvi,nlmta,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       call mpi_allreduce(MPI_IN_PLACE,sumdqvr,nlmta*3,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       call mpi_allreduce(MPI_IN_PLACE,sumdqvi,nlmta*3,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       call sum_EXX_force_terms(force_l,bdwr,bdwi,dbdwr,dbdwi,sumqvr,sumqvi,sumdqvr,sumdqvi)
    else
       call mpi_allreduce(MPI_IN_PLACE,sumdel,nfft,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
#ifdef FFTW3
       call m_FFT_exx(nfout,sumdel,DIRECT,OFF) ! sum delta(G)
#else
       call m_FFT_WF(ELECTRON,nfout,sumdel,DIRECT,OFF) ! sum delta(G)
#endif
       if(kimg==1) then
          call map_FFT_box_on_WFG(ik,vxw_t(1,1),vxw_t(1,1),sumdel)
       else
          call map_FFT_box_on_WFG(ik,vxw_t(1,1),vxw_t(1,2),sumdel)
       end if

       if(modnrm == EXECUT)  then
          fac = 1.d0/dble(nfftwf)
       else
          fac = 1.d0/dble(nfftwf)**2
          if(present(eexx)) eexx = eexx*fac
       end if
       vxw_t= fac*vxw_t
       if(modnrm == EXECUT) then
          call mpi_allreduce(MPI_IN_PLACE,sumqvr,nlmta,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
          call mpi_allreduce(MPI_IN_PLACE,sumqvi,nlmta,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
          call add_Vx_hard_part_2D(ik,vxw_t,sumqvr,sumqvi)
       endif
       do ig=ista_g1k(ik),iend_g1k(ik)
          iadd = ig-ista_g1k(ik)+1
          vxw(iadd,1:kimg) = vxw_t(ig,1:kimg)
       enddo

    endif
    deallocate(rho)
    deallocate(phi)
    deallocate(wfn)
    deallocate(wfm)
    deallocate(sumdel)
    deallocate(rhogr)
    if(kimg==2) deallocate(rhogi)
    deallocate(wfsr)
    if(kimg==2) deallocate(wfsi)
    deallocate(cosgt)
    if(kimg==2) deallocate(singt)
    deallocate(vxw_t)
    if(modnrm == EXECUT) then
       deallocate(fsr)
       deallocate(fsi)
       deallocate(qvr)
       deallocate(qvi)
       deallocate(sumqvr)
       deallocate(sumqvi)
       if(force_mode)then
          deallocate(dfsr)
          deallocate(dfsi)
          deallocate(dqvr)
          deallocate(dqvi)
          deallocate(gqvr)
          deallocate(gqvi)
          deallocate(sumdqvr)
          deallocate(sumdqvi)
       endif
    end if

    if(sw_rspace_hyb==ON.and.modnrm==EXECUT)then
       deallocate(rhor)
       deallocate(rhoi)
       deallocate(afft)
    endif

    call tstatc0_end(id_sname)
  end subroutine apply_Vx_to_WF_2D

  subroutine apply_Vx_to_WF_2D_rs_dgm(ispin,ib,ik,kgw,kgv,wfr,wfi,bdwr,bdwi,vxw,eexx,eonly,force_l,dbdwr,dbdwi)
    implicit none
    integer, intent(in)                       :: ispin, ib,ik,kgw,kgv
    real(kind=DP), intent(in), dimension(kgw) :: wfr, wfi
    real(kind=DP), intent(in), dimension(nlmta) :: bdwr, bdwi
    real(kind=DP), intent(out), optional, dimension(kgv,kimg) :: vxw
    real(kind=DP), intent(out), optional :: eexx
    logical, intent(in), optional :: eonly
    real(kind=DP), allocatable, dimension(:,:) :: vxw_t
    real(kind=DP), intent(in), dimension(nlmta,3), optional :: dbdwr, dbdwi
    real(kind=DP), intent(inout), dimension(natm,3), optional :: force_l
    logical :: force_mode=.false.
    real(kind=DP) :: exx
    integer :: kk,m,ig,ii,jop,kop,jtrs,mm,ib1
    integer :: jkbz
    real(kind=DP) :: fac, ph, occupation, dnorm
    real(kind=DP) :: skg(3)
    real(kind=DP), allocatable, dimension(:) :: rho,phi !d(nfft)
    real(kind=DP), allocatable, dimension(:) :: wfn !d(nfft)
    real(kind=DP), allocatable, dimension(:,:) :: wfm !d(nfft)
    real(kind=DP), allocatable, dimension(:) :: sumdel !d(nfft)
    real(kind=DP), allocatable, dimension(:) :: rhogr, rhogi !d(kgp)
    real(kind=DP), allocatable, dimension(:) :: wfsr,wfsi !d(kg1)
    real(kind=DP), allocatable, dimension(:) :: cosgt, singt !d(kg1)
    real(kind=DP), allocatable, dimension(:) :: fsr, fsi !d(nlmta)
    real(kind=DP), allocatable, dimension(:,:) :: dfsr, dfsi !d(nlmta,3)
    real(kind=DP), allocatable, dimension(:,:,:) :: dfsr2, dfsi2 !d(nlmta,3)
    real(kind=DP), allocatable, dimension(:) :: qvr, qvi !d(nlmta)
    real(kind=DP), allocatable, dimension(:,:) :: qvr2, qvi2 !d(nlmta)
    real(kind=DP), allocatable, dimension(:,:) :: dqvr, dqvi !d(nlmta,3)
    real(kind=DP), allocatable, dimension(:,:) :: gqvr, gqvi !d(nlmta,3)
    real(kind=DP), allocatable, dimension(:,:,:) :: dqvr2, dqvi2 !d(nlmta,3)
    real(kind=DP), allocatable, dimension(:,:,:) :: gqvr2, gqvi2 !d(nlmta,3)
    real(kind=DP), allocatable, dimension(:) :: sumqvr, sumqvi !d(nlmta)
    real(kind=DP), allocatable, dimension(:,:) :: sumdqvr, sumdqvi !d(nlmta)
    real(kind=DP), allocatable, dimension(:) :: rhor,rhoi
    real(kind=DP), allocatable, dimension(:) :: afft
    real(kind=DP), allocatable, dimension(:,:) :: rhor2,rhoi2
    real(kind=DP), allocatable, dimension(:,:) :: fsr2,fsi2
    real(kind=DP), allocatable, dimension(:,:) :: rhogr2,rhogi2
    real(kind=DP) :: ifac
    logical :: eo
    integer :: iadd
    integer :: ikk
    integer,save  :: id_sname = -1,id_sname1=-2,id_sname2=-3,id_sname3=-4
    integer,save  :: id_sname_cdfft=-1,id_sname_cdfft2=-1,id_sname_cdfft3=-1
    call tstatc0_begin('apply_Vx_to_WF ',id_sname)
    eo=.false.
    if(present(eonly))then
      eo = eonly
    endif
    if(present(eexx))then
      eexx = 0.d0
    endif

    force_mode = present(force_l).and.present(dbdwr).and.present(dbdwi)

    allocate(rho(nfft))
    allocate(phi(nfft))
    allocate(wfn(nfft))
    allocate(wfm(nfft,ista_nval:iend_nval))
    allocate(sumdel(nfft))
    allocate(rhogr(nmax_G_hyb));rhogr=0.d0
    allocate(rhogr2(nmax_G_hyb,ista_nval:iend_nval));rhogr2=0.d0
    allocate(rhor(nfftp_exx_nonpara/2));rhor=0.d0
    allocate(rhoi(nfftp_exx_nonpara/2));rhoi=0.d0
    allocate(afft(nfftp_exx_nonpara));afft=0.d0
    allocate(rhor2(nfftp_exx_nonpara/2,ista_nval:iend_nval));rhor2=0.d0
    allocate(rhoi2(nfftp_exx_nonpara/2,ista_nval:iend_nval));rhoi2=0.d0
    allocate(fsr2(ista_nval:iend_nval,nlmta))
    allocate(fsi2(ista_nval:iend_nval,nlmta))
    if(kimg==2) then
        allocate(rhogi(nmax_G_hyb))
        allocate(rhogi2(nmax_G_hyb,ista_nval:iend_nval));rhogi2=0.d0
        rhogi = 0.d0
    endif
    allocate(wfsr(kg1))
    if(kimg==2) allocate(wfsi(kg1))
    allocate(cosgt(kg1))
    if(kimg==2) allocate(singt(kg1))
    allocate(vxw_t(kg1,kimg));vxw_t=0.d0
    allocate(fsr(nlmta))
    allocate(fsi(nlmta))
    allocate(qvr(nlmta))
    allocate(qvi(nlmta))
    allocate(qvr2(nlmta,ista_nval:iend_nval))
    allocate(qvi2(nlmta,ista_nval:iend_nval))
    allocate(sumqvr(nlmta))
    allocate(sumqvi(nlmta))
    fac = 1.d0/dble(nfftwf)
    sumqvr = 0.d0
    sumqvi = 0.d0
    if(force_mode)then
       allocate(dfsr(nlmta,3))
       allocate(dfsi(nlmta,3))
       allocate(dfsr2(nlmta,ista_nval:iend_nval,3))
       allocate(dfsi2(nlmta,ista_nval:iend_nval,3))
       allocate(dqvr2(nlmta,ista_nval:iend_nval,3))
       allocate(dqvi2(nlmta,ista_nval:iend_nval,3))
       allocate(gqvr2(nlmta,ista_nval:iend_nval,3))
       allocate(gqvi2(nlmta,ista_nval:iend_nval,3))
       allocate(sumdqvr(nlmta,3))
       allocate(sumdqvi(nlmta,3))
       sumdqvr=0.d0;sumdqvi=0.d0
    endif

    call map_WFG_on_FFT_box(ik,wfr,wfi,wfn)
#ifdef FFTW3
    call m_FFT_exx(nfout,wfn,INVERSE,OFF) ! wfn(R)
#else
    call m_FFT_WF(ELECTRON,nfout,wfn,INVERSE,OFF) ! wfn(R)
#endif

    sumdel = 0.d0
    ikk = k_index(ik)
    do jkbz=1,kv3bz
       if(.not.q_on_k_centered_mesh(jkbz,ikk)) cycle
#ifdef MEMORY_SAVE_EXX
       call m_ES_EXX_ylm_each_k(iqmk(jkbz,ikk))
#ifdef MEMORY_SAVE_MORE_EXX
       call qitgft_qmk_each_k(iqmk(jkbz,ikk))
#endif
#endif
       ! q(bz) = S*k(ibz)
       kk = ikp(jkbz) + ispin - 1
       kop = iops(jkbz)
       jop = iopt(jkbz)
       jtrs = itrs(jkbz)
       if(kimg==1) then
          do ig=1,iba_for_innerloop(kk)
             ii = nbase_for_innerloop(ig,kk)
             skg(1:3) = vkxyz_for_innerloop(kk,1:3,BUCS) + ngabc(ii,1:3)
             ph = PAI2 * dot_product(skg,tau(:,kop,BUCS))
             cosgt(ig) = cos(ph)
          end do
       else
          do ig=1,iba_for_innerloop(kk)
             ii = nbase_for_innerloop(ig,kk)
             skg(1:3) = vkxyz_for_innerloop(kk,1:3,BUCS) + ngabc(ii,1:3)
             ph = PAI2 * dot_product(skg,tau(:,kop,BUCS))
             cosgt(ig) = cos(ph)
             singt(ig) = sin(ph)
          end do
       end if

       do m=ista_nval,iend_nval
          if(occup_val(m,kk) < DELTA) cycle
          if(force_mode)then
             call get_Rot_betar_dot_WFs(kk,m,jop,kop,jtrs,fsr2(m,:),fsi2(m,:),dfsr2(:,m,:),dfsi2(:,m,:))
          else
             call get_Rot_betar_dot_WFs(kk,m,jop,kop,jtrs,fsr2(m,:),fsi2(m,:))
          endif
       enddo
       call add_RHOG_hard_part_rs3(iqmk(jkbz,ikk),kk,rhor2,rhoi2,bdwr,bdwi,fsr2,fsi2)
       do m=ista_nval,iend_nval
          if(occup_val(m,kk) < DELTA) cycle
!          if(sw_rsb==ON)then
!             if(.not.overlap(ib,m))cycle
!          endif
          mm = map_z_nval(m)
          occupation = wbz(jkbz) * occup_val(m,kk)/qwgt_for_innerloop(kk)/dble(kv3_for_innerloop)
          if(kimg==1) then
             dnorm = 0.d0
             do ig=1,iba_for_innerloop(kk)
                wfsr(ig) = wfv(ig,mm,kk,1) * cosgt(ig)
                dnorm = dnorm + wfsr(ig)*wfsr(ig)
             end do
             dnorm = 1.d0/sqrt(dnorm)
             wfsr = wfsr * dnorm
             call map_Rot_WFG_on_FFT_box(kk,jop,jtrs,wfsr,wfsr,wfm(:,m))
          else
             if(jtrs==0) then
                do ig=1,iba_for_innerloop(kk)
                   wfsr(ig) = wfv(ig,mm,kk,1) * cosgt(ig) - wfv(ig,mm,kk,2) * singt(ig)
                   wfsi(ig) = wfv(ig,mm,kk,2) * cosgt(ig) + wfv(ig,mm,kk,1) * singt(ig)
                end do
             else
                do ig=1,iba_for_innerloop(kk)
                   wfsr(ig) = wfv(ig,mm,kk,1) * cosgt(ig) - wfv(ig,mm,kk,2) * singt(ig)
                   wfsi(ig) = -wfv(ig,mm,kk,2) * cosgt(ig) - wfv(ig,mm,kk,1) * singt(ig)
                end do
             end if
             call map_Rot_WFG_on_FFT_box(kk,jop,jtrs,wfsr,wfsi,wfm(:,m))
          end if
#ifdef FFTW3
          call m_FFT_exx(nfout,wfm(:,m),INVERSE,OFF) ! wfm(R)
#else
          call m_FFT_WF(ELECTRON,nfout,wfm(:,m),INVERSE,OFF) ! wfm(R)
#endif
          call product_on_FFT_box(wfn,wfm(:,m),rho) ! rho_soft(R)
#ifdef FFTW3
          call m_FFT_exx(nfout,rho,DIRECT,OFF) ! rho_soft(G)
#else
          call m_FFT_WF(ELECTRON,nfout,rho,DIRECT,OFF) ! rho_soft(G)
#endif

          call map_FFT_box_on_RHOG(rhogr2(:,m),rhogi2(:,m),rho)

          !! rho(G) = rho_soft(G) + rho_hard(G)
          !do ii=1, min( kg, nmax_G_hyb )
          rhogr2(:,m) = fac * rhogr2(:,m)
          rhogi2(:,m) = fac * rhogi2(:,m)
          !end do
          do ii=kg+1, nmax_G_hyb
             rhogr2(ii,m) = 0.d0
             rhogi2(ii,m) = 0.d0
          end do
          !fsr(:) = fsr2(m,:)
          !fsi(:) = fsi2(m,:)

          !rhor(:) = rhor2(:,m)
          !rhoi(:) = rhoi2(:,m)
          call map_RHOG_on_FFT_box_hard(rhor2(:,m),rhoi2(:,m),afft)
          call m_FFT_CD0_exx(nfout,afft,DIRECT)
          call map_FFT_box_on_RHOG_hard(rhor2(:,m),rhoi2(:,m),afft)
          do ii=1,nmax_G_hyb
             rhogr2(ii,m) = rhogr2(ii,m) + rhor2(ii,m)
             rhogi2(ii,m) = rhogi2(ii,m) + rhoi2(ii,m)
          end do

          if(present(eexx))then
             exx = 0.d0
             call sum_rho_vc_rho_2D(rhogr2(:,m),rhogi2(:,m),vc(1,iqmk(jkbz,ikk)),exx)
             eexx = eexx + occupation * exx
             if(eo) cycle
          endif
          !do ii=1,nmax_G_hyb
          rhogr2(:,m) = vc(:,iqmk(jkbz,ikk)) * rhogr2(:,m) ! phi(G) = Vc(G,q) * rho(G)
          rhogi2(:,m) = vc(:,iqmk(jkbz,ikk)) * rhogi2(:,m) ! phi(G) = Vc(G,q) * rho(G)
          !end do
          call map_RHOG_on_FFT_box_hard_inv(rhogr2(:,m),rhogi2(:,m),afft)
          call m_FFT_CD0_exx(nfout,afft,INVERSE)
          call map_FFT_box_on_RHOG_hard_inv(rhor2(:,m),rhoi2(:,m),afft)
          !rhor2(:,m) = rhor(:)
          !rhoi2(:,m) = rhoi(:)
       enddo
       if(.not.eo)then
          if(force_mode)then
             call integrate_QijVnm_rs3(iqmk(jkbz,ikk),rhor2,rhoi2,fsr2,fsi2,qvr2,qvi2,dfsr2,dfsi2,dqvr2,dqvi2,gqvr2,gqvi2)
          else
             call integrate_QijVnm_rs3(iqmk(jkbz,ikk),rhor2,rhoi2,fsr2,fsi2,qvr2,qvi2)
          endif
! =============================================================== 12.5Exp
          do m=ista_nval,iend_nval
             if(occup_val(m,kk) < DELTA) cycle
             occupation = wbz(jkbz) * occup_val(m,kk)/qwgt_for_innerloop(kk)/dble(kv3_for_innerloop)
             !rhogr(:) = rhogr2(:,m)
             !rhogi(:) = rhogi2(:,m)
             !if(.not.force_mode)then
             !endif
             !do ii=1,nlmta
             !   write(0,*) m,ii,qvr2(m,ii),qvi2(m,ii)
             !enddo
             if(force_mode)then
                sumqvr = sumqvr + occupation * qvr2(:,m) ! sum_m qvr_m
                sumqvi = sumqvi + occupation * qvi2(:,m) ! sum_m qvi_m
                sumdqvr = sumdqvr + occupation * (dqvr2(:,m,:)+gqvr2(:,m,:)) ! sum_m (dqvr_m + gqvr_m)
                sumdqvi = sumdqvi + occupation * (dqvi2(:,m,:)+gqvi2(:,m,:)) ! sum_m (dqvi_m + gqvi_m)
             else
                call map_RHOG_on_FFT_box(rhogr2(:,m),rhogi2(:,m),phi)
#ifdef FFTW3
                call m_FFT_exx(nfout,phi,INVERSE,OFF) ! phi(R)
#else
                call m_FFT_WF(ELECTRON,nfout,phi,INVERSE,OFF) ! phi(R)
#endif
                call product_on_FFT_box(phi,wfm(:,m),rho) ! delta(R)
                sumdel = sumdel + occupation * rho ! sum delta(R)
                !! Bin = Sum_m fm * Hi[Vnm]
                sumqvr = sumqvr + occupation * qvr2(:,m) ! sum_m qvr_m
                sumqvi = sumqvi + occupation * qvi2(:,m) ! sum_m qvi_m
             endif
          end do
       endif
    end do
    if(force_mode)then
       !! dEXX/dR(n) = Sum_i Re[ <psi_n|dbeta_i/dR> * Bin +  <psi_n|beta_i> * Cin ]
       call mpi_allreduce(MPI_IN_PLACE,sumqvr,nlmta,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       call mpi_allreduce(MPI_IN_PLACE,sumqvi,nlmta,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       call mpi_allreduce(MPI_IN_PLACE,sumdqvr,nlmta*3,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       call mpi_allreduce(MPI_IN_PLACE,sumdqvi,nlmta*3,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       call sum_EXX_force_terms(force_l,bdwr,bdwi,dbdwr,dbdwi,sumqvr,sumqvi,sumdqvr,sumdqvi)
    else
       call mpi_allreduce(MPI_IN_PLACE,sumdel,nfft,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
#ifdef FFTW3
       call m_FFT_exx(nfout,sumdel,DIRECT,OFF) ! sum delta(G)
#else
       call m_FFT_WF(ELECTRON,nfout,sumdel,DIRECT,OFF) ! sum delta(G)
#endif
       if(kimg==1) then
          call map_FFT_box_on_WFG(ik,vxw_t(1,1),vxw_t(1,1),sumdel)
       else
          call map_FFT_box_on_WFG(ik,vxw_t(1,1),vxw_t(1,2),sumdel)
       end if

       fac = 1.d0/dble(nfftwf)
       vxw_t= fac*vxw_t

       call mpi_allreduce(MPI_IN_PLACE,sumqvr,nlmta,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       call mpi_allreduce(MPI_IN_PLACE,sumqvi,nlmta,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       call add_Vx_hard_part_2D(ik,vxw_t,sumqvr,sumqvi)

       do ig=ista_g1k(ik),iend_g1k(ik)
          iadd = ig-ista_g1k(ik)+1
          vxw(iadd,1:kimg) = vxw_t(ig,1:kimg)
       enddo

    endif
    deallocate(rho)
    deallocate(phi)
    deallocate(wfn)
    deallocate(wfm)
    deallocate(sumdel)
    deallocate(rhogr)
    if(kimg==2) deallocate(rhogi)
    deallocate(wfsr)
    if(kimg==2) deallocate(wfsi)
    deallocate(cosgt)
    if(kimg==2) deallocate(singt)
    deallocate(vxw_t)
    deallocate(fsr)
    deallocate(fsi)
    deallocate(qvr)
    deallocate(qvi)
    deallocate(sumqvr)
    deallocate(sumqvi)
    if(force_mode)then
       deallocate(dfsr)
       deallocate(dfsi)
       deallocate(dfsr2)
       deallocate(dfsi2)
       deallocate(dqvr2)
       deallocate(dqvi2)
       deallocate(gqvr2)
       deallocate(gqvi2)
       deallocate(sumdqvr)
       deallocate(sumdqvi)
    endif

    deallocate(rhor)
    deallocate(rhoi)
    deallocate(afft)
    deallocate(rhor2)
    deallocate(rhoi2)
    deallocate(fsr2)
    deallocate(fsi2)

    call tstatc0_end(id_sname)
  end subroutine apply_Vx_to_WF_2D_rs_dgm

  subroutine m_ES_EXX_energy2(eexx)
    implicit none
    real(kind=DP), intent(out) :: eexx
    real(kind=DP) :: eexx_mpi
    integer :: ispin,ik,ib,ib1,ierr
    eexx = 0.0d0
    do ispin=1,nspin,af+1
       do ik=ispin, kv3+ispin-nspin, nspin
          if(map_k(ik) /= myrank_k) cycle ! MPI
          do ib=1,np_e   ! MPI
             if(occup_l(ib,ik) < DELTA) cycle
             eexx = eexx + eexx_kb(ib,ik)
          enddo
       enddo
    enddo
    eexx = 0.5d0*eexx/dble(kv3/nspin)
    if(nspin==1) eexx = 2.d0 * eexx

    if(npes>1) then
       call mpi_allreduce(eexx,eexx_mpi,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
       eexx = eexx_mpi
       call mpi_allreduce(eexx,eexx_mpi,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ge_world,ierr)
       eexx = eexx_mpi
    end if

  end subroutine m_ES_EXX_energy2

  subroutine m_ES_EXX_energy(eexx)
    implicit none
    real(kind=DP), intent(out) :: eexx
    integer :: ispin,ik
    integer :: id_sname=-1
    call tstatc0_begin('m_ES_EXX_energy ',id_sname,level=1)

    do ispin=1,nspin,af+1
       do ik=ispin, kv3+ispin-nspin, nspin
          if(map_k(ik) /= myrank_k) cycle ! MPI
          call m_ES_EXX_eigenvalue_for_each_k(ispin,ik,.false.,iupdate=2)
       enddo
    enddo
    call m_ES_EXX_energy2(eexx)

    call tstatc0_end(id_sname)
  end subroutine m_ES_EXX_energy

  subroutine product_on_FFT_box(wfn,wfm,rho)
    implicit none
    real(kind=DP), intent(in), dimension(nfft) :: wfn,wfm
    real(kind=DP), intent(out), dimension(nfft) :: rho

    integer :: i, ire, iim, ik
    rho=0.d0
    if(kv3_for_innerloop/nspin==1 .and. k_symmetry(1) == GAMMA) then
       do i=1,nfft_exx,2
          ire = i
          iim = i+1
          rho(ire) = wfn(ire)*wfm(ire)
          rho(iim) = 0.d0
       end do
    else
       do i=1,nfft_exx,2
          ire = i
          iim = i+1
          rho(ire) = wfn(ire)*wfm(ire) + wfn(iim)*wfm(iim)
          rho(iim) = wfn(ire)*wfm(iim) - wfn(iim)*wfm(ire)
       end do
    end if
  end subroutine product_on_FFT_box
  subroutine product_on_FFT_box_3D(wfn,wfm,rho,lsize,ibsize,nfft_l)
    implicit none
    integer, intent(in) :: lsize, ibsize, nfft_l
    real(kind=DP), intent(in), dimension(lsize*kimg,ibsize) :: wfn,wfm
    real(kind=DP), intent(out), dimension(lsize*kimg) :: rho

    integer :: i, ire, iim, ik, j

    if(kv3_for_innerloop/nspin==1 .and. k_symmetry(1) == GAMMA) then
       do j = 1, ibsize
! === DEBUG by tkato 2014/04/14 ================================================
!      do i=1,nfft_l,2
       do i=1,nfft_l*kimg,2
! ==============================================================================
          ire = i
          iim = i+1
          rho(ire) = wfn(ire,j)*wfm(ire,j)
          rho(iim) = 0.d0
       end do
       end do
    else
       do j = 1, ibsize
#ifdef FFT_USE_SSL2_3D
          do i = 1, nel_fft_x(myrank_g)
#else
          do i = 1, nfft_l
#endif
          ire = 2*i-1
          iim = 2*i
          rho(ire) = wfn(ire,j)*wfm(ire,j) + wfn(iim,j)*wfm(iim,j)
          rho(iim) = wfn(ire,j)*wfm(iim,j) - wfn(iim,j)*wfm(ire,j)
       end do
       end do
    end if
  end subroutine product_on_FFT_box_3D

  subroutine sum_rho_vc_rho_2D(rhor,rhoi,vc,exx)
    implicit none
    real(kind=DP), intent(in), dimension(nmax_G_hyb) :: rhor,rhoi
    real(kind=DP), intent(in), dimension(nmax_G_hyb) :: vc
    real(kind=DP), intent(out) :: exx

    integer :: ig

    exx = 0.d0

    if(kimg==1) then
       do ig=1,nmax_G_hyb
          exx = exx + vc(ig) * rhor(ig)*rhor(ig)
       end do
    else
       do ig=1,nmax_G_hyb
          exx = exx + vc(ig) * (rhor(ig)*rhor(ig)+rhoi(ig)*rhoi(ig))
       end do
    end if
  end subroutine sum_rho_vc_rho_2D

  subroutine sum_rho_vc_rho(rhor,rhoi,vc,exx)
    implicit none
    real(kind=DP), intent(in), dimension(ista_kngp:iend_kngp) :: rhor,rhoi
    real(kind=DP), intent(in), dimension(ista_kngp:iend_kngp) :: vc
    real(kind=DP), intent(out) :: exx

    integer :: ig

    exx = 0.d0

    if(kimg==1) then
       do ig=ista_kngp,iend_kngp
          exx = exx + vc(ig) * rhor(ig)*rhor(ig)
       end do
    else
       do ig=ista_kngp,iend_kngp
          exx = exx + vc(ig) * (rhor(ig)*rhor(ig)+rhoi(ig)*rhoi(ig))
       end do
    end if
  end subroutine sum_rho_vc_rho

  subroutine map_WFG_on_FFT_box(ik,wfr,wfi,bfft)
    implicit none
    integer, intent(in)                         :: ik
    real(kind=DP), intent(in), dimension(kg1)   :: wfr, wfi
    real(kind=DP), intent(out), dimension(nfft) :: bfft

    integer :: i,i1,ri, j, i2, ii

    bfft = 0.d0
    if(k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          i1 = igf_exx(1)
          bfft(i1) = wfr(1)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba(ik)
             i = nbase(ii,ik)
             if(i<=kg_exx)then
               i1 = igf_exx(i)
               bfft(i1) = wfr(ii)
             endif
             j = nbase_gamma(ii,2)
             if(j<=kg_exx)then
               i2 = igf_exx(j)
               bfft(i2) = wfr(ii)
             endif
          end do
       else ! kimg == 2
          i1 = 2*igf_exx(1) - 1
          bfft(i1)   = wfr(1)
          bfft(i1+1) = wfi(1)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba(ik)
             i = nbase(ii,ik)
             if(i<=kg_exx)then
               i1 = 2*igf_exx(i)-1
               bfft(i1  ) = wfr(ii)
               bfft(i1+1) = wfi(ii)
             endif
             j = nbase_gamma(ii,2)
             if(j<=kg_exx)then
               i2 = 2*igf_exx(j)-1
               bfft(i2  ) = wfr(ii)
               bfft(i2+1) = -wfi(ii)
             endif
          end do
       end if
    else
      if(kimg == 1) then
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
        do i = 1, iba(ik)
           if(nbase(i,ik)<=kg_exx)then
             i1 = igf_exx(nbase(i,ik))
             bfft(i1) = wfr(i)
           endif
        end do
      else ! kimg == 2
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
        do i = 1, iba(ik)
           if(nbase(i,ik)<=kg_exx)then
             i1 = 2*igf_exx(nbase(i,ik)) - 1
             bfft(i1) = wfr(i)
             i2 = 2*igf_exx(nbase(i,ik))
             bfft(i2) = wfi(i)
           endif
        end do
      end if
    end if
  end subroutine map_WFG_on_FFT_box

  subroutine map_Rot_WFG_on_FFT_box(ik,nop,jtrs,wfr,wfi,bfft)
    implicit none
    integer, intent(in)                         :: ik,nop,jtrs
    real(kind=DP), intent(in), dimension(kg1_for_innerloop)   :: wfr, wfi
    real(kind=DP), intent(out), dimension(nfft) :: bfft

    integer :: i,i1,ri, j, i2, ii, ip

    bfft = 0.d0
    if(k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          i1 = igf_exx(1)
          bfft(i1) = wfr(1)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba_for_innerloop(ik)
             i = nbase_for_innerloop(ii,ik)
             if(ngpt_exx(i,nop,jtrs)<=kg_exx)then
                i1 = igf_exx(ngpt_exx(i,nop,jtrs))
                bfft(i1) = wfr(ii)
             endif
             j = nbase_gamma_for_innerloop(ii,2)
             if(ngpt_exx(j,nop,jtrs)<=kg_exx)then
                i2 = igf_exx(ngpt_exx(j,nop,jtrs))
                bfft(i2) = wfr(ii)
             endif
          end do
       else ! kimg == 2
          i1 = 2*igf_exx(1) - 1
          bfft(i1)   = wfr(1)
          bfft(i1+1) = wfi(1)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba_for_innerloop(ik)
             i = nbase_for_innerloop(ii,ik)
             if(ngpt_exx(i,nop,jtrs)<=kg_exx)then
               i1 = 2*igf_exx(ngpt_exx(i,nop,jtrs))-1
               bfft(i1  ) = wfr(ii)
               bfft(i1+1) = wfi(ii)
             endif
             j = nbase_gamma_for_innerloop(ii,2)
             if(ngpt_exx(j,nop,jtrs)<=kg_exx)then
               i2 = 2*igf_exx(ngpt_exx(j,nop,jtrs))-1
               bfft(i2  ) = wfr(ii)
               bfft(i2+1) = -wfi(ii)
             endif
          end do
       end if
    else
      if(kimg == 1) then
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
        do i = 1, iba_for_innerloop(ik)
           if(ngpt_exx(nbase_for_innerloop(i,ik),nop,jtrs)<=kg_exx)then
             i1 = igf_exx(ngpt_exx(nbase_for_innerloop(i,ik),nop,jtrs))
             bfft(i1) = wfr(i)
           endif
        end do
      else ! kimg == 2
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
        do i = 1, iba_for_innerloop(ik)
           if(ngpt_exx(nbase_for_innerloop(i,ik),nop,jtrs)<=kg_exx)then
             ip = igf_exx(ngpt_exx(nbase_for_innerloop(i,ik),nop,jtrs))
             i1 = 2*ip - 1
             bfft(i1) = wfr(i)
             i2 = 2*ip
             bfft(i2) = wfi(i)
           endif
        end do
      end if
    end if
  end subroutine map_Rot_WFG_on_FFT_box

  subroutine map_FFT_box_on_WFG(ik,wfr,wfi,bfft)
    implicit none
    integer, intent(in)                        :: ik
    real(kind=DP), intent(out), dimension(kg1) :: wfr, wfi
    real(kind=DP), intent(in), dimension(nfft) :: bfft

    integer :: i,i1,ri, j, i2, ii

    if(k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          wfr = 0.d0

          i1 = igf_exx(1)
          wfr(1) = bfft(i1)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba(ik)
             i = nbase(ii,ik)
             if(i<=kg_exx)then
               i1 = igf_exx(i)
               wfr(ii) = bfft(i1)
             endif
          end do
       else ! kimg == 2
          wfr = 0.d0
          wfi = 0.d0

          i1 = 2*igf_exx(1) - 1
          wfr(1) = bfft(i1)
          wfi(1) = bfft(i1+1)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba(ik)
             i = nbase(ii,ik)
             if(i<=kg_exx)then
               i1 = 2*igf_exx(i)-1
               wfr(ii) = bfft(i1  )
               wfi(ii) = bfft(i1+1)
             endif
          end do
       end if
    else
       if(kimg == 1) then
          wfr = 0.d0
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do i = 1, iba(ik)
             if(nbase(i,ik)<=kg_exx)then
               i1 = igf_exx(nbase(i,ik))
               wfr(i) = bfft(i1)
             endif
          end do
       else ! kimg == 2
          wfr = 0.d0
          wfi = 0.d0
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do i = 1, iba(ik)
             if(nbase(i,ik)<=kg_exx)then
               i1 = 2*igf_exx(nbase(i,ik)) - 1
               wfr(i) = bfft(i1)
               i2 = 2*igf_exx(nbase(i,ik))
               wfi(i) = bfft(i2)
             endif
          end do
       end if
    end if
  end subroutine map_FFT_box_on_WFG

  subroutine map_FFT_box_on_RHOG(rhor,rhoi,afft)
    implicit none
    real(kind=DP), intent(out), dimension(nmax_G_hyb) :: rhor, rhoi
    real(kind=DP), intent(in), dimension(nfft) :: afft

    integer :: i,i1,i2

    if(kimg == 1) then
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
       do i = 1, min( kg_exx, nmax_G_hyb )
          i1 = igf_exx(i)
          rhor(i) = afft(i1)
       end do
    else ! kimg==2
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
       do i = 1, min( kg_exx, nmax_G_hyb )
          i1 = 2*igf_exx(i) - 1
          rhor(i) = afft(i1)
          i2 = 2*igf_exx(i)
          rhoi(i) = afft(i2)
       end do
    end if
  end subroutine map_FFT_box_on_RHOG

  subroutine map_RHOG_on_FFT_box(rhor,rhoi,afft)
    implicit none
    real(kind=DP), intent(in), dimension(nmax_G_hyb) :: rhor, rhoi
    real(kind=DP), intent(out), dimension(nfft) :: afft

    integer :: i,i1,i2

    afft = 0.d0

    if(kimg == 1) then
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
       do i = 1, min( kg_exx, nmax_G_hyb )
          i1 = igf_exx(i)
          afft(i1) = rhor(i)
       end do
    else ! kimg==2
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
       do i = 1, min( kg_exx, nmax_G_hyb )
          i1 = 2*igf_exx(i) - 1
          afft(i1) = rhor(i)
          i2 = 2*igf_exx(i)
          afft(i2) = rhoi(i)
       end do
    end if
  end subroutine map_RHOG_on_FFT_box

  subroutine map_RHOG_on_FFT_box_hard(rhor,rhoi,afft)
    implicit none
    real(kind=DP), intent(in), dimension(nfftp_exx_nonpara/2) :: rhor, rhoi
    real(kind=DP), intent(out), dimension(nfftp_exx_nonpara) :: afft

    integer :: i,i1,i2

!    afft(:) = 0.d0
    do i = 1, nfftp_exx_nonpara/2
!       i1 = (igfp_l(i)-1)*kimg+1
       i1 = (i-1)*kimg+1
       afft(i1)   = rhor(i)
       afft(i1+1) = rhoi(i)
    end do
  end subroutine map_RHOG_on_FFT_box_hard

  subroutine map_RHOG_on_FFT_box_hard_inv(rhor,rhoi,afft)
    implicit none
    real(kind=DP), intent(in), dimension(nfftp_exx_nonpara/2) :: rhor, rhoi
    real(kind=DP), intent(out), dimension(nfftp_exx_nonpara) :: afft
    integer :: i,i1,i2

    afft(:) = 0.d0

    do i = 1, kgp_exx
       i1 = (igfp_exx(i)-1)*kimg+1
       afft(i1) = rhor(i)
       afft(i1+1) = rhoi(i)
    end do
  end subroutine map_RHOG_on_FFT_box_hard_inv

  subroutine map_FFT_box_on_RHOG_hard(rhor,rhoi,afft)
    implicit none
    real(kind=DP), intent(in), dimension(nfftp_exx_nonpara) :: afft
    real(kind=DP), intent(out), dimension(nfftp_exx_nonpara/2) :: rhor, rhoi
    real(kind=DP) :: rinplw
    integer :: i,i1,i2
    integer :: id_sname = -1

    rinplw = 1.d0/product(fft_box_size_CD_exx(1:3,1))
    rhor(:)=0.d0
    rhoi(:)=0.d0
    do i = 1, kgp_exx
       i1 = (igfp_exx(i)-1)*kimg+1
       rhor(i) = afft(i1)*rinplw
       rhoi(i) = afft(i1+1)*rinplw
    end do
  end subroutine map_FFT_box_on_RHOG_hard

  subroutine map_FFT_box_on_RHOG_hard_inv(rhor,rhoi,afft)
    implicit none
    real(kind=DP), intent(in), dimension(nfftp_exx_nonpara) :: afft
    real(kind=DP), intent(out), dimension(nfftp_exx_nonpara/2) :: rhor, rhoi
    real(kind=DP) :: rinplw
    integer :: i,i1,i2

    rinplw = 1.d0/product(fft_box_size_CD_exx(1:3,1))
!    rhor(:)=0.d0
!    rhoi(:)=0.d0
    do i = 1, nfftp_exx_nonpara/2
       i1 = (i-1)*kimg+1
       rhor(i) = afft(i1)*rinplw
       rhoi(i) = afft(i1+1)*rinplw
    end do
  end subroutine map_FFT_box_on_RHOG_hard_inv

  subroutine m_ES_EXX_ngpt()
    implicit none

    integer :: i,j,iopr,ii
    integer :: ia,ib,ic
    integer :: namin,namax,nbmin,nbmax,ncmin,ncmax
    integer, allocatable, dimension(:,:,:) :: g_list
    integer, allocatable, dimension(:) :: ngpt_exx_tmp, ngpt_exx0_tmp
    integer,save  :: id_sname = -1
    call tstatc0_begin('m_ES_EXX_ngpt ',id_sname,level=1)

!!$    allocate(ngpt_exx(kgp,nopr,0:ntrs)); ngpt_exx = 0
    allocate(ngpt_exx(kg,nopr,0:ntrs)); ngpt_exx = 0

    !! Time reversal symmetry

    if(ntrs>0) then

       namax = n_rGpv(1); nbmax = n_rGpv(2); ncmax = n_rGpv(3)
       namin = -namax   ; nbmin = -nbmax   ; ncmin = -ncmax
       allocate(g_list(namin:namax,nbmin:nbmax,ncmin:ncmax)); g_list = 0

       do i = ista_kngp, iend_kngp
          ia = ngabc_kngp_l(i,1)
          ib = ngabc_kngp_l(i,2)
          ic = ngabc_kngp_l(i,3)
          g_list(ia,ib,ic) = i
       end do
       call mpi_allreduce(MPI_IN_PLACE, g_list,(namax-namin+1)*(nbmax-nbmin+1)*(ncmax-ncmin+1), &
       &  mpi_integer, mpi_sum,mpi_ke_world,ierr)

    end if

    if(npes > 1) then
       allocate(ngpt_exx0_tmp(kgp))
       do iopr=1,nopr
          ngpt_exx0_tmp = 0
          do i = ista_kngp, iend_kngp
             ngpt_exx0_tmp(i) = ngpt_l(i,iopr)
          end do
          call mpi_allreduce(MPI_IN_PLACE,ngpt_exx0_tmp,kgp,mpi_integer,mpi_sum,mpi_ke_world,ierr)
          ngpt_exx(1:kg,iopr,0) = ngpt_exx0_tmp(1:kg)

          if(ntrs>0) then
             allocate(ngpt_exx_tmp(kg))
             ngpt_exx_tmp = 0
             do i=1,kg
                ii = ngpt_exx0_tmp(i)
                if(ista_kngp<=ii .and. ii<=iend_kngp) then
                   ia = -ngabc_kngp_l(ii,1)
                   ib = -ngabc_kngp_l(ii,2)
                   ic = -ngabc_kngp_l(ii,3)
                   ngpt_exx_tmp(i) = g_list(ia,ib,ic)
                end if
             end do
             call mpi_allreduce(MPI_IN_PLACE, ngpt_exx_tmp, kg, mpi_integer, mpi_sum,mpi_ke_world,ierr)
             ngpt_exx(1:kg,iopr,1) = ngpt_exx_tmp(1:kg)
! === ASMS ===
             deallocate(ngpt_exx_tmp)
! === ASMS ===
          end if
       end do
       deallocate(ngpt_exx0_tmp)
    else
       do iopr=1,nopr
          if(ntrs==0) then
             do i = 1, kg
                ngpt_exx(i,iopr,0) = ngpt_l(i,iopr)
             end do
          else
             do i = 1, kg
                ngpt_exx(i,iopr,0) = ngpt_l(i,iopr)
                ii = ngpt_l(i,iopr)
                ia = -ngabc_kngp_l(ii,1)
                ib = -ngabc_kngp_l(ii,2)
                ic = -ngabc_kngp_l(ii,3)
                ngpt_exx(i,iopr,1) = g_list(ia,ib,ic)
             end do
          end if
       end do
    end if

!!$    allocate(ngpt_exx_tmp(kgp))
!!$    do iopr=1,nopr
!!$       ngpt_exx_tmp = 0
!!$       do i=1,kgp
!!$          ii = ngpt_exx(i,iopr,0)
!!$          if(ista_kngp<=ii .and. ii<=iend_kngp) then
!!$             ia = -ngabc_kngp_l(ii,1)
!!$             ib = -ngabc_kngp_l(ii,2)
!!$             ic = -ngabc_kngp_l(ii,3)
!!$             ngpt_exx_tmp(i) = g_list(ia,ib,ic)
!!$          end if
!!$       end do
!!$       call mpi_allreduce(MPI_IN_PLACE, ngpt_exx_tmp, kgp, mpi_integer, mpi_sum,mpi_ke_world,ierr)
!!$       ngpt_exx(:,iopr,1) = ngpt_exx_tmp(:)
!!$    end do
!!$    deallocate(ngpt_exx_tmp)
!!$
!!$          ia = -ngabc(ii,1)
!!$          ib = -ngabc(ii,2)
!!$          ic = -ngabc(ii,3)
!!$          ngpt_exx(i,iopr,1) = g_list(ia,ib,ic)
!!$       end do
!!$    end do
    if(ntrs>0) deallocate(g_list)

    j = 0
    do ii = 0, ntrs
       do iopr = 1, nopr
          do i = 1, kg
             if(ngpt_exx(i,iopr,ii) <= 0) j = j + 1
          end do
       end do
    end do
    if(j >= 1) then
       write(nfout,'(" !! check of ngpt_exx")')
       do ii = 0, ntrs
          do iopr = 1, nopr
             do i = 1, kg
                if(ngpt_exx(i,iopr,ii) <= 0) write(nfout,'(" ngpt_exx(",i8,",",i8,",",i8,") = ",i20)') &
                & i,iopr,ii, ngpt_exx(i,iopr,ii)
             end do
          end do
       end do
       write(nfout,'(" !! total number of negative values for ngpt_exx = ",i8)') j
       write(nfout,'(" !! out of check of ngpt_exx")')
    end if

    if(sw_change_axis /= ON) then
! === Make FFT box index arrays. ===============================================
    call Parallelize_wf_onto_fft_exx_3D(nfout,fft_box_size_WF,igf,nbase,nbase_gamma, &
   &                                   k_symmetry,GAMMA,kg,kg_gamma,kv3)
    call Parallelize_fft_onto_wf_rhog_3D(nfout,igf,kg,nfft)
    call Parallelize_wf_onto_fft_rhog_3D(nfout,fft_box_size_WF,igf,kg)
! ==============================================================================
    endif
    call tstatc0_end(id_sname)
    !!!stop 'Check: G_list'
  end subroutine m_ES_EXX_ngpt

  subroutine m_ES_EXX_occup(nfout)
    implicit none
    integer, intent(in) :: nfout

    real(kind=DP) :: chgt, chg(nspin), zet
    integer :: ik,ib,is,it
    integer,save  :: id_sname = -1
    call tstatc0_begin('m_ES_EXX_occup ',id_sname,level=1)

    chg(1:nspin) = 0.d0
    do it=1,ntyp
       chgt = ival(it) * iatom(it) + qex(it)
       if(nspin==1) then
          chg(1) = chg(1) + chgt
       else
          zet = zeta1(it)
          chg(1) = chg(1) + chgt * (1.d0+zet)/2.d0
          chg(2) = chg(2) + chgt * (1.d0-zet)/2.d0
       end if
    end do
    if(nspin==1) chg(1) = chg(1)/2.d0

    nval = nint(maxval(chg)+0.1d0)

    occup_l(1:np_e,ista_k:iend_k) = 0.d0
    do ik=1,kv3,af+1
       if(map_k(ik) /= myrank_k) cycle ! MPI
       do ib = 1, neg
          if(neordr(ib,ik)>nval) cycle
          if(map_e(ib) == myrank_e) then
!!$       do ib=ista_e,iend_e,istep_e
!!$          !!if(ib>nval) cycle
!!$          if(neordr(ib,ik)>nval) cycle
             occup_l(map_z(ib),ik) = kv3*qwgt(ik)
          end if
       end do
    end do
    call tstatc0_end(id_sname)
  end subroutine m_ES_EXX_occup

  subroutine m_ES_EXX_store_wfv()
    implicit none
    if(.not.allocated(wfv_old)) allocate(wfv_old(kg1,nval,kv3,kimg))
    if(.not.allocated(occup_val_old)) allocate(occup_val_old(nval,kv3))
    if(modnrm == EXECUT) then
       if(.not.allocated(fsr_exx_old)) allocate(fsr_exx_old(nval,nlmta,kv3))
       if(.not.allocated(fsi_exx_old)) allocate(fsi_exx_old(nval,nlmta,kv3))
    end if

    wfv_old  = wfv
    occup_val_old = occup_val
    nval_old = nval
    if(modnrm == EXECUT) then
       fsr_exx_old = fsr_exx
       fsi_exx_old = fsi_exx
    end if
  end subroutine m_ES_EXX_store_wfv

  subroutine m_ES_EXX_restore_wfv()
    implicit none
    if(allocated(wfv)) deallocate(wfv)
    if(allocated(occup_val)) deallocate(occup_val)
    if(modnrm == EXECUT) then
       if(allocated(fsr_exx)) deallocate(fsr_exx)
       if(allocated(fsi_exx)) deallocate(fsi_exx)
    end if
    nval = nval_old
    allocate(wfv(kg1,nval,kv3,kimg))
    allocate(occup_val(nval,kv3))
    if(modnrm == EXECUT) then
       allocate(fsr_exx(nval,nlmta,kv3))
       allocate(fsi_exx(nval,nlmta,kv3))
    end if
    wfv = wfv_old
    occup_val = occup_val_old
    if(modnrm == EXECUT) then
       fsr_exx = fsr_exx_old
       fsi_exx = fsi_exx_old
    end if
    deallocate(wfv_old)
    deallocate(occup_val_old)
    if(modnrm == EXECUT) then
       deallocate(fsr_exx_old)
       deallocate(fsi_exx_old)
    end if
  end subroutine m_ES_EXX_restore_wfv

  function q_on_k_centered_mesh(iq,ik) result(ltf)
    implicit none
    integer, intent(in) :: iq, ik

    logical :: ltf
    real(kind=DP) :: qv(3), kv(3)
    integer :: dv(3)

    qv(1:3) = vkbz(iq,1:3)
    kv(1:3) = vkxyz_for_outerloop(ik,1:3,BUCS)

    dv(1:3) = nint( (qv(1:3) - kv(1:3)) * n_bz_mesh(1:3) )

    ltf = .false.
    if(mod(dv(1),reduction_factor_exx(1)) == 0 .and. &
     & mod(dv(2),reduction_factor_exx(2)) == 0 .and. &
     & mod(dv(3),reduction_factor_exx(3)) == 0 ) then
       ltf = .true.
    end if
  end function q_on_k_centered_mesh

!! For debug
  function check_ibz_k(ik) result(ltf)
    implicit none
    integer, intent(in) :: ik

    logical :: ltf
    real(kind=DP) :: vkx, vky, vkz

    open(1000,file="ibz.data")
    read(1000,*) vkx, vky, vkz
    close(1000)

    ltf = .true.
    if(abs(vkxyz(ik,1,BUCS) - vkx) < DELTA07 .and. &
     & abs(vkxyz(ik,2,BUCS) - vky) < DELTA07 .and. &
     & abs(vkxyz(ik,3,BUCS) - vkz) < DELTA07 ) then
       ltf = .false.
       write(nfout,'("IBZ k: ",3f20.5)') vkxyz(ik,1:3,BUCS)
    end if
  end function check_ibz_k

  function check_bz_k(ik) result(ltf)
    implicit none
    integer, intent(in) :: ik

    logical :: ltf
    real(kind=DP) :: vkx, vky, vkz

    open(1100,file="bz.data")
    read(1100,*) vkx, vky, vkz
    close(1100)

    ltf = .true.
    if(abs(vkbz(ik,1) - vkx) < DELTA07 .and. &
     & abs(vkbz(ik,2) - vky) < DELTA07 .and. &
     & abs(vkbz(ik,3) - vkz) < DELTA07 ) then
       write(nfout,'(" BZ k: ",3f20.5)') vkbz(ik,1:3)
       ltf = .false.
    end if
  end function check_bz_k

!! For ultrasoft pseudopotentials
  subroutine m_ES_EXX_ylm()
    implicit none

    integer :: i,n,ik,ig
    integer :: iend_kngp0
    real(kind=DP), allocatable, dimension(:) :: ylm_t ! d(kngp,n)
    real(kind=DP), allocatable, dimension(:,:) :: gqmk ! d(ista_kngp:iend_kngp,3)
    real(kind=DP), allocatable, dimension(:) :: gqmkr ! d(ista_kngp:iend_kngp)

    real(kind=DP) :: kg(3), ttr(6), g2
    integer,save  :: id_sname = -1
    call tstatc0_begin('m_ES_EXX_ylm ',id_sname,level=1)

    call getttr(rltv,ttr)

    call m_PP_find_maximum_l(n)    ! n-1: maximum l
    n = (n-1) + (n-1) + 1
    n = n*n
    if(.not.allocated(ylm_exx)) then
      if(sw_change_axis==ON)then
        allocate(ylm_exx(nmax_G_hyb,n,nqmk))
      else
        allocate(ylm_exx(ista_kngp:iend_kngp,n,nqmk))
      endif
    endif

! === ASMS ====
    iend_kngp0 = iend_kngp
    if(iend_kngp0.gt.nmax_G_hyb) iend_kngp0 = nmax_G_hyb
! === ASMS ====

    if( ista_kngp <= nmax_G_hyb ) then
       allocate(gqmk(ista_kngp:iend_kngp0,3))
       allocate(gqmkr(ista_kngp:iend_kngp0))
    endif
    do ik=1,nqmk
       do ig=ista_kngp,iend_kngp0
          kg(1:3) = qmk(ik,1:3) + ngabc_kngp_l(ig,1:3)
          g2          = ttr(1)*kg(1)*kg(1) &
          &           + ttr(2)*kg(2)*kg(2) &
          &           + ttr(3)*kg(3)*kg(3) &
          &           + ttr(4)*kg(1)*kg(2) &
          &           + ttr(5)*kg(2)*kg(3) &
          &           + ttr(6)*kg(3)*kg(1)
          gqmk(ig,1:3) = kg(1:3)
          gqmkr(ig) = sqrt(g2)
       end do
       if(sw_change_axis==ON) then
       if(npes > 1) then
          allocate(ylm_t(nmax_G_hyb))
          do i=1,n
             ylm_t = 0.d0
             if ( ista_kngp <= nmax_G_hyb ) then
                call m_pwBS_sphrp_exx( i, rltv, ista_kngp, iend_kngp0, &
                     &                 gqmk, gqmkr, ylm_t(ista_kngp:iend_kngp0) )
             endif
             call mpi_allreduce(ylm_t,ylm_exx(1,i,ik),nmax_G_hyb,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
          end do
          deallocate(ylm_t)
       else
          do i=1,n
             call m_pwBS_sphrp_exx(i,rltv,1,nmax_G_hyb,gqmk,gqmkr,ylm_exx(1,i,ik))
          end do
       end if
       else
          do i=1,n
             if ( ista_kngp <= nmax_G_hyb ) then
                call m_pwBS_sphrp_exx(i,rltv,ista_kngp,nmax_G_hyb,gqmk,gqmkr,ylm_exx(ista_kngp,i,ik))
             endif
          end do
       endif

    end do
    !!stop 'm_ES_EXX_ylm'
    if ( allocated(gqmk) )  deallocate(gqmk)
    if ( allocated(gqmkr) ) deallocate(gqmkr)

    call tstatc0_end(id_sname)
  end subroutine m_ES_EXX_ylm

  subroutine check_qitg()
    implicit none
    integer :: iq, ips , ipe

    ips = 1
    ipe = nmax_G_hyb
    ips = ista_kngp
    ipe = iend_kngp

    do iq=1,nqitg
       write(nfout,'("iq=",i5,1x,"qitg_l=",f20.5,1x,"qitg_exx=",f20.5,1x)') iq, qitg_l(ips,iq), qitg_exx(ips,iq,1)/univol
    end do

    do iq=1,nqitg
       write(nfout,'("iq=",i5,1x,"diff=",f20.5)') iq, sum(qitg_l(ips:ipe,iq)-qitg_exx(ips:ipe,iq,1)/univol)
    end do

    call phase_error_with_msg(nfout,'check_qitg',__LINE__,__FILE__)
  end subroutine check_qitg

  subroutine check_qitg_qmk()
    implicit none
    integer :: iq, ik, ips, ipe

    ips = 1
    ipe = nmax_G_hyb
    ips = ista_kngp
    ipe = iend_kngp

    do ik=1,nqmk
    do iq=1,nqitg
       write(nfout,'("iq=",i5,1x,"qitg_exx=",f20.5,1x)') iq, qitg_exx(ips,iq,ik)/univol
    end do
    end do

    do ik=1,nqmk
    do iq=1,nqitg
       write(nfout,'("iq=",i5,1x,"sum=",f20.5)') iq, sum(qitg_exx(ips:ipe,iq,ik)/univol)
    end do
    end do

    call phase_error_with_msg(nfout,'check_qitg_qmk',__LINE__,__FILE__)
  end subroutine check_qitg_qmk

  subroutine check_ylm_exx()
    implicit none
    integer :: i, n, ik, ips, ipe

    call m_PP_find_maximum_l(n)    ! n-1: maximum l
    n = (n-1) + (n-1) + 1
    n = n*n

    ips = 1
    ipe = nmax_G_hyb
    ips = ista_kngp
    ipe = iend_kngp

    do ik=1,nqmk
    do i=1,n
       write(nfout,'("i=",i5,1x,"ylm_exx=",f20.5,1x)') i, ylm_exx(ips,i,ik)
    end do
    end do

    do ik=1,nqmk
    do i=1,n
       write(nfout,'("i=",i5,1x,"sum=",f20.5)')  i, sum(qitg_exx(ips:ipe,i,ik))
    end do
    end do

    call phase_error_with_msg(nfout,'check_ylm_exx',__LINE__,__FILE__)
  end subroutine check_ylm_exx

  subroutine hard_part_preproc()
    integer :: it,ncache
    call m_PP_find_maximum_l(n)    ! n-1: maximum l
    n = (n-1) + (n-1) + 1
    allocate(il3(n**2))
    call substitute_il3(n**2,il3) ! -(b_Elec..)

    allocate(nqitg_sp(ntyp)); allocate(nqitg_sp0(ntyp))
    allocate(iq2l3(nqitg))
    allocate(nc(mcritical,nqitg));nc=0
! =================================== Added by K. Tagami =======
    nqitg_sp = 0; nqitg_sp0 = 0; iq2l3 = 0
! ==============================================================
    call m_PP_set_index_arrays1(nfout,ntyp,nqitg,mcritical,n**2,il3 &
    & ,maxm,mc,nqitg_sp,nqitg_sp0,iq2l3,nc,.true.)
    allocate(nc2lmt1(mc,maxm,nqitg))
    allocate(nc2lmt2(mc,maxm,nqitg))
    allocate(nc2n(mc,maxm,nqitg))
! ==================================== Added by K. Tagami ======
    nc2lmt1 = 0;  nc2lmt2 = 0; nc2n = 0
! =============================================================
    call m_PP_set_index_arrays2(nfout,mc,maxm,nqitg,mcritical,n**2,il3,iq2l3 &
    & ,nc2lmt1,nc2lmt2,nc2n,nc,.true.) ! -> nc2lmt1, nc2lmt2, nc2n, nc

    ncache = (m_CtrlP_cachesize()*1024)*3/4
    if(ncache == 0) then
       ibsize =  nmax_G_hyb
    else
       iwidth = 0
       do it = 1, ntyp
          iwidth = max(iwidth,nqitg_sp(it)-nqitg_sp0(it)+1)
       end do
       ibsize=ncache/(8*(iwidth + 2))
    end if
    if(ibsize.gt.nmax_G_hyb) ibsize = nmax_G_hyb
#ifdef HYB_DISABLE_BLOCKING
    ibsize = nmax_G_hyb
#endif
    preproc_done = .true.
  end subroutine hard_part_preproc

  subroutine add_RHOG_hard_part_2D(iqmk,rhogr,rhogi,fnr,fni,fmr,fmi)
    implicit none
    integer, intent(in) :: iqmk
    real(kind=DP), intent(inout) :: rhogr(nmax_G_hyb)
    real(kind=DP), intent(inout) :: rhogi(nmax_G_hyb)
    real(kind=DP), intent(in) :: fnr(nlmta)
    real(kind=DP), intent(in) :: fni(nlmta)
    real(kind=DP), intent(in) :: fmr(nlmta)
    real(kind=DP), intent(in) :: fmi(nlmta)

    integer :: it,ia,i,lmt1,lmt2,np,il1,il2,tau1,tau2,ilmta1,ilmta2
    integer :: l3,ilm3,iil
    real(kind=DP) :: ph, hr, hi, fac, f0
    real(kind=DP), allocatable :: cosgt(:), singt(:)

    real(kind=DP), allocatable, dimension(:) :: ylm_sum,ylm_sum2

    real(kind=DP), allocatable :: qitg_red(:,:),ylm_red(:,:)
    real(kind=DP), allocatable :: rhogr_red(:),rhogi_red(:)
    real(kind=DP) :: yr,yi
    integer :: ibl1,ibl2,iq,inn,ip,m

    integer,save  :: id_sname = -1

    call tstatc0_begin('add_RHOG_hard_part ',id_sname)

    if(.not.preproc_done) call hard_part_preproc()

    allocate(cosgt(ibsize))
    allocate(singt(ibsize))
    allocate(ylm_sum(ibsize));ylm_sum=0.d0
    allocate(ylm_sum2(ibsize));ylm_sum2=0.d0

    allocate(qitg_red(ibsize,nqitg))
    allocate(ylm_red(ibsize,n*n))
    allocate(rhogr_red(ibsize))
    allocate(rhogi_red(ibsize))
! --
    do ibl1=1, nmax_G_hyb, ibsize
       rhogr_red=0.d0;  rhogi_red=0.d0

       ibl2=ibl1+ibsize-1

       if (ibl2.gt.nmax_G_hyb) ibl2 = nmax_G_hyb

       do iq=1,nqitg
          do i=1,ibl2-ibl1+1
             qitg_red(i,iq) = qitg_exx(i+ibl1-1,iq,iqmk)
          enddo
       enddo
       do inn=1,n*n
          do i=1,ibl2-ibl1+1
             ylm_red(i,inn) = ylm_exx(i+ibl1-1,inn,iqmk)
          enddo
       enddo

       do ia=1,natm
          it = ityp(ia)
          if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
          do i=1,ibl2-ibl1+1
             ph = PAI2 * dot_product(pos(ia,1:3),ngabc_red(i+ibl1-1,1:3)+qmk(iqmk,1:3))
             cosgt(i) = cos(ph);  singt(i) = sin(ph)
          end do
          do iq = nqitg_sp0(it), nqitg_sp(it)
             l3 = iq2l3(iq)
             iil = mod(l3,4)
             ylm_sum  = 0.d0;  ylm_sum2 = 0.d0

             do m = 1, 2*l3+1
                ilm3 = l3*l3+m
                yr=0.d0;  yi=0.d0

                do ip=1,nc(m,iq)
                   lmt1 = nc2lmt1(ip,m,iq);  lmt2 = nc2lmt2(ip,m,iq)
                   ilmta1 = lmta(lmt1,ia);   ilmta2 = lmta(lmt2,ia)
                   np = nc2n(ip,m,iq)
                   hr = fnr(ilmta1)*fmr(ilmta2) + fni(ilmta1)*fmi(ilmta2)
                   hi = fnr(ilmta1)*fmi(ilmta2) - fni(ilmta1)*fmr(ilmta2)
                   f0 = iwei(ia) * dl2p(lmt1,lmt2,np,it)
                   yr = yr + f0*hr;  yi = yi + f0*hi
                   !do i=1,ibl2-ibl1+1
                   !   ylm_sum(i)  = ylm_sum(i)  + f0*ylm_red(i,ilm3)*hr
                   !   ylm_sum2(i) = ylm_sum2(i) + f0*ylm_red(i,ilm3)*hi
                   !enddo
                enddo
                do i=1,ibl2-ibl1+1
                   ylm_sum(i)  = ylm_sum(i)  + yr*ylm_red(i,ilm3)
                   ylm_sum2(i) = ylm_sum2(i) + yi*ylm_red(i,ilm3)
                enddo
             enddo

             if(iil==0) then
                do i=1,ibl2-ibl1+1
                   fac =  qitg_red(i,iq)
                   rhogr_red(i) = rhogr_red(i) + fac * (cosgt(i)*ylm_sum(i) +singt(i)*ylm_sum2(i))
                   rhogi_red(i) = rhogi_red(i) + fac * (cosgt(i)*ylm_sum2(i)-singt(i)*ylm_sum(i))
                end do
             else if(iil==1) then
                do i=1,ibl2-ibl1+1
                   fac =  qitg_red(i,iq)
                   rhogr_red(i) = rhogr_red(i) + fac * (cosgt(i)*ylm_sum2(i)-singt(i)*ylm_sum(i))
                   rhogi_red(i) = rhogi_red(i) - fac * (cosgt(i)*ylm_sum(i) +singt(i)*ylm_sum2(i))
                end do
             else if(iil==2) then
                do i=1,ibl2-ibl1+1
                   fac =  qitg_red(i,iq)
                   rhogr_red(i) = rhogr_red(i) - fac * (cosgt(i)*ylm_sum(i) +singt(i)*ylm_sum2(i))
                   rhogi_red(i) = rhogi_red(i) - fac * (cosgt(i)*ylm_sum2(i)-singt(i)*ylm_sum(i))
                end do
             else if(iil==3) then
                do i=1,ibl2-ibl1+1
                   fac =  qitg_red(i,iq)
                   rhogr_red(i) = rhogr_red(i) - fac * (cosgt(i)*ylm_sum2(i)-singt(i)*ylm_sum(i))
                   rhogi_red(i) = rhogi_red(i) + fac * (cosgt(i)*ylm_sum(i) +singt(i)*ylm_sum2(i))
                end do
             end if
          enddo
       enddo

       do i=ibl1,ibl2
          rhogr(i) = rhogr(i) + rhogr_red(i-ibl1+1)
          rhogi(i) = rhogi(i) + rhogi_red(i-ibl1+1)
       enddo

    enddo

    deallocate(cosgt); deallocate(singt)
    deallocate(ylm_sum); deallocate(ylm_sum2)
    deallocate(qitg_red); deallocate(ylm_red)
    deallocate(rhogr_red); deallocate(rhogi_red)

    call tstatc0_end(id_sname)
  end subroutine add_RHOG_hard_part_2D

  subroutine add_RHOG_hard_part(iqmk,rhogr,rhogi,fnr,fni,fmr,fmi)
    implicit none
    integer, intent(in) :: iqmk
    real(kind=DP), intent(inout) :: rhogr(ista_kngp:iend_kngp)
    real(kind=DP), intent(inout) :: rhogi(ista_kngp:iend_kngp)
    real(kind=DP), intent(in) :: fnr(nlmta)
    real(kind=DP), intent(in) :: fni(nlmta)
    real(kind=DP), intent(in) :: fmr(nlmta)
    real(kind=DP), intent(in) :: fmi(nlmta)

    integer :: it,ia,i,lmt1,lmt2,np,il1,il2,tau1,tau2,ilmta1,ilmta2
    integer :: l3,ilm3,iil
    real(kind=DP) :: ph, hr, hi, fac, f0
    real(kind=DP), allocatable :: cosgt(:), singt(:)

    real(kind=DP), allocatable, dimension(:) :: ylm_sum,ylm_sum2

    real(kind=DP), allocatable :: qitg_red(:,:),ylm_red(:,:)
    real(kind=DP), allocatable :: rhogr_red(:),rhogi_red(:)
    real(kind=DP) :: yr,yi
    integer :: ibl1,ibl2,iq,inn,ip,m, ips, ipe

    integer,save  :: id_sname = -1

    call tstatc0_begin('add_RHOG_hard_part ',id_sname)

    if(.not.preproc_done) call hard_part_preproc()

    allocate(cosgt(ibsize))
    allocate(singt(ibsize))
    allocate(ylm_sum(ibsize));ylm_sum=0.d0
    allocate(ylm_sum2(ibsize));ylm_sum2=0.d0

    allocate(qitg_red(ibsize,nqitg))
    allocate(ylm_red(ibsize,n*n))
    allocate(rhogr_red(ibsize))
    allocate(rhogi_red(ibsize))
! --
    ips = 1
    ipe = nmax_G_hyb
    ips = ista_kngp
    ipe = iend_kngp

    do ibl1=ips, ipe, ibsize
       rhogr_red=0.d0;  rhogi_red=0.d0

       ibl2=min(ipe,ibl1+ibsize-1)
       do iq=1,nqitg
          do i=1,ibl2-ibl1+1
             qitg_red(i,iq) = qitg_exx(i+ibl1-1,iq,iqmk)
          enddo
       enddo
       do inn=1,n*n
          do i=1,ibl2-ibl1+1
             ylm_red(i,inn) = ylm_exx(i+ibl1-1,inn,iqmk)
          enddo
       enddo

       do ia=1,natm
          it = ityp(ia)
          if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
          do i=1,ibl2-ibl1+1
! === DEBUG by tkato 2015/03/19 ================================================
!            ph = PAI2 * dot_product(pos(ia,1:3),ngabc(i+ibl1-1,1:3)+qmk(iqmk,1:3))
             ph = PAI2 * dot_product(pos(ia,1:3),ngabc_kngp_l(i+ibl1-1,1:3)+qmk(iqmk,1:3))
! ==============================================================================
             cosgt(i) = cos(ph);  singt(i) = sin(ph)
          end do
          do iq = nqitg_sp0(it), nqitg_sp(it)
             l3 = iq2l3(iq)
             iil = mod(l3,4)
             ylm_sum  = 0.d0;  ylm_sum2 = 0.d0

             do m = 1, 2*l3+1
                ilm3 = l3*l3+m
                yr=0.d0;  yi=0.d0

                do ip=1,nc(m,iq)
                   lmt1 = nc2lmt1(ip,m,iq);  lmt2 = nc2lmt2(ip,m,iq)
                   ilmta1 = lmta(lmt1,ia);   ilmta2 = lmta(lmt2,ia)
                   np = nc2n(ip,m,iq)
                   hr = fnr(ilmta1)*fmr(ilmta2) + fni(ilmta1)*fmi(ilmta2)
                   hi = fnr(ilmta1)*fmi(ilmta2) - fni(ilmta1)*fmr(ilmta2)
                   f0 = iwei(ia) * dl2p(lmt1,lmt2,np,it)
                   yr = yr + f0*hr;  yi = yi + f0*hi
                   !do i=1,ibl2-ibl1+1
                   !   ylm_sum(i)  = ylm_sum(i)  + f0*ylm_red(i,ilm3)*hr
                   !   ylm_sum2(i) = ylm_sum2(i) + f0*ylm_red(i,ilm3)*hi
                   !enddo
                enddo
                do i=1,ibl2-ibl1+1
                   ylm_sum(i)  = ylm_sum(i)  + yr*ylm_red(i,ilm3)
                   ylm_sum2(i) = ylm_sum2(i) + yi*ylm_red(i,ilm3)
                enddo
             enddo

             if(iil==0) then
                do i=1,ibl2-ibl1+1
                   fac =  qitg_red(i,iq)
                   rhogr_red(i) = rhogr_red(i) + fac * (cosgt(i)*ylm_sum(i) +singt(i)*ylm_sum2(i))
                   rhogi_red(i) = rhogi_red(i) + fac * (cosgt(i)*ylm_sum2(i)-singt(i)*ylm_sum(i))
                end do
             else if(iil==1) then
                do i=1,ibl2-ibl1+1
                   fac =  qitg_red(i,iq)
                   rhogr_red(i) = rhogr_red(i) + fac * (cosgt(i)*ylm_sum2(i)-singt(i)*ylm_sum(i))
                   rhogi_red(i) = rhogi_red(i) - fac * (cosgt(i)*ylm_sum(i) +singt(i)*ylm_sum2(i))
                end do
             else if(iil==2) then
                do i=1,ibl2-ibl1+1
                   fac =  qitg_red(i,iq)
                   rhogr_red(i) = rhogr_red(i) - fac * (cosgt(i)*ylm_sum(i) +singt(i)*ylm_sum2(i))
                   rhogi_red(i) = rhogi_red(i) - fac * (cosgt(i)*ylm_sum2(i)-singt(i)*ylm_sum(i))
                end do
             else if(iil==3) then
                do i=1,ibl2-ibl1+1
                   fac =  qitg_red(i,iq)
                   rhogr_red(i) = rhogr_red(i) - fac * (cosgt(i)*ylm_sum2(i)-singt(i)*ylm_sum(i))
                   rhogi_red(i) = rhogi_red(i) + fac * (cosgt(i)*ylm_sum(i) +singt(i)*ylm_sum2(i))
                end do
             end if
          enddo
       enddo

       do i=ibl1,ibl2
          rhogr(i) = rhogr(i) + rhogr_red(i-ibl1+1)
          rhogi(i) = rhogi(i) + rhogi_red(i-ibl1+1)
       enddo

    enddo

    deallocate(cosgt); deallocate(singt)
    deallocate(ylm_sum); deallocate(ylm_sum2)
    deallocate(qitg_red); deallocate(ylm_red)
    deallocate(rhogr_red); deallocate(rhogi_red)

    call tstatc0_end(id_sname)
  end subroutine add_RHOG_hard_part

#ifdef HYBRID_DGEMM

  subroutine add_RHOG_hard_part_(iqmk_,rhogr,rhogi,fnr,fni,fmr,fmi,ik,ispin)
    implicit none
    integer, intent(in) :: iqmk_(kv3bz,kv3), ik, ispin
    real(kind=DP), intent(inout) :: rhogr(ista_kngp:iend_kngp,nval,kv3bz)
    real(kind=DP), intent(inout) :: rhogi(ista_kngp:iend_kngp,nval,kv3bz)
    real(kind=DP), intent(in) :: fnr(nlmta)
    real(kind=DP), intent(in) :: fni(nlmta)
    real(kind=DP), intent(in) :: fmr(nlmta,nval,kv3bz)
    real(kind=DP), intent(in) :: fmi(nlmta,nval,kv3bz)

    integer :: it,ia,i,lmt1,lmt2,np,il1,il2,tau1,tau2,ilmta1,ilmta2
    integer :: l3,ilm3,iil
    real(kind=DP) :: ph, hr, hi, fac, f0
    real(kind=DP), allocatable :: cosgt(:,:), singt(:,:)

    real(kind=DP), allocatable :: qitg_red(:,:),ylm_red(:,:)
    real(kind=DP), allocatable :: rhogr_red(:,:),rhogi_red(:,:)
    real(kind=DP), allocatable :: yr(:,:),yi(:,:)
    real(kind=DP), allocatable :: tmp(:,:,:)
    integer :: ibl1,ibl2,iq,inn,ip,m

    integer,save  :: id_sname = -1
    integer :: jkbz, kk, jtrs, iqmk, m_, max_iaiqm, iaiqm

    call tstatc0_begin('add_RHOG_hard_part_ ',id_sname)

    if(.not.preproc_done) call hard_part_preproc()

    allocate(cosgt(ibsize,natm))
    allocate(singt(ibsize,natm))

    allocate(qitg_red(ibsize,nqitg))
    allocate(ylm_red(ibsize,n*n))
    allocate(rhogr_red(ibsize,nval))
    allocate(rhogi_red(ibsize,nval))

    max_iaiqm = 0
    do ia = 1, natm
       it = ityp(ia)
       if(m_PP_include_vanderbilt_pot(it) == SKIP) cycle
       do iq = nqitg_sp0(it), nqitg_sp(it)
          l3 = iq2l3(iq)
          do m = 1, 2*l3+1
             max_iaiqm = max_iaiqm + 1
          end do
       end do
    end do
    allocate(yr(nval,max_iaiqm))
    allocate(yi(nval,max_iaiqm))
    allocate(tmp(ibsize,max_iaiqm,4))

    do ibl1 = ista_kngp, iend_kngp, ibsize
       ibl2 = ibl1 + ibsize - 1
       if(ibl2 .gt. iend_kngp) ibl2 = iend_kngp
       do jkbz = 1, kv3bz
          if(.not. q_on_k_centered_mesh(jkbz,ik)) cycle
          kk = ikp(jkbz) + ispin - 1
          iqmk = iqmk_(jkbz,ik)
          do iq = 1, nqitg
             do i = 1, ibl2-ibl1+1
                qitg_red(i,iq) = qitg_exx(i+ibl1-1,iq,iqmk)
             end do
          end do
          do inn = 1, n*n
             do i = 1, ibl2-ibl1+1
                ylm_red(i,inn) = ylm_exx(i+ibl1-1,inn,iqmk)
             end do
          end do
          do ia = 1, natm
             it = ityp(ia)
             if(m_PP_include_vanderbilt_pot(it) == SKIP) cycle
             do i = 1, ibl2-ibl1+1
! === DEBUG by tkato 2015/03/19 ================================================
!               ph = PAI2 * dot_product(pos(ia,1:3),ngabc(i+ibl1-1,1:3)+qmk(iqmk,1:3))
                ph = PAI2 * dot_product(pos(ia,1:3),ngabc_kngp_l(i+ibl1-1,1:3)+qmk(iqmk,1:3))
! ==============================================================================
                cosgt(i,ia) = cos(ph)
                singt(i,ia) = sin(ph)
             end do
          end do
          yr = 0.0d0
          yi = 0.0d0
          iaiqm = 0
          do ia = 1, natm
             it = ityp(ia)
             if(m_PP_include_vanderbilt_pot(it) == SKIP) cycle
             do iq = nqitg_sp0(it), nqitg_sp(it)
                l3 = iq2l3(iq)
                iil = mod(l3,4)
                do m = 1, 2*l3+1
                   iaiqm = iaiqm + 1
                   do m_ = 1, nval
                      if(occup_val(m_,kk) < DELTA) cycle
                      do ip = 1, nc(m,iq)
                         lmt1 = nc2lmt1(ip,m,iq)
                         lmt2 = nc2lmt2(ip,m,iq)
                         ilmta1 = lmta(lmt1,ia)
                         ilmta2 = lmta(lmt2,ia)
                         np = nc2n(ip,m,iq)
                         hr = fnr(ilmta1)*fmr(ilmta2,m_,jkbz) + fni(ilmta1)*fmi(ilmta2,m_,jkbz)
                         hi = fnr(ilmta1)*fmi(ilmta2,m_,jkbz) - fni(ilmta1)*fmr(ilmta2,m_,jkbz)
                         f0 = iwei(ia) * dl2p(lmt1,lmt2,np,it)
                         yr(m_,iaiqm) = yr(m_,iaiqm) + f0*hr
                         yi(m_,iaiqm) = yi(m_,iaiqm) + f0*hi
                      end do
                   end do
                   if(iil == 0) then
                      do i = 1, ibl2-ibl1+1
                         fac =  qitg_red(i,iq)
                         tmp(i,iaiqm,1) =    fac*cosgt(i,ia)*ylm_red(i,l3*l3+m)
                         tmp(i,iaiqm,2) =    fac*singt(i,ia)*ylm_red(i,l3*l3+m)
                         tmp(i,iaiqm,3) =  - fac*singt(i,ia)*ylm_red(i,l3*l3+m)
                         tmp(i,iaiqm,4) =    fac*cosgt(i,ia)*ylm_red(i,l3*l3+m)
                      end do
                   else if(iil == 1) then
                      do i = 1, ibl2-ibl1+1
                         fac =  qitg_red(i,iq)
                         tmp(i,iaiqm,1) =  - fac*singt(i,ia)*ylm_red(i,l3*l3+m)
                         tmp(i,iaiqm,2) =    fac*cosgt(i,ia)*ylm_red(i,l3*l3+m)
                         tmp(i,iaiqm,3) =  - fac*cosgt(i,ia)*ylm_red(i,l3*l3+m)
                         tmp(i,iaiqm,4) =  - fac*singt(i,ia)*ylm_red(i,l3*l3+m)
                      end do
                   else if(iil == 2) then
                      do i = 1, ibl2-ibl1+1
                         fac =  qitg_red(i,iq)
                         tmp(i,iaiqm,1) =  - fac*cosgt(i,ia)*ylm_red(i,l3*l3+m)
                         tmp(i,iaiqm,2) =  - fac*singt(i,ia)*ylm_red(i,l3*l3+m)
                         tmp(i,iaiqm,3) =    fac*singt(i,ia)*ylm_red(i,l3*l3+m)
                         tmp(i,iaiqm,4) =  - fac*cosgt(i,ia)*ylm_red(i,l3*l3+m)
                      end do
                   else if(iil == 3) then
                      do i = 1, ibl2-ibl1+1
                         fac =  qitg_red(i,iq)
                         tmp(i,iaiqm,1) =    fac*singt(i,ia)*ylm_red(i,l3*l3+m)
                         tmp(i,iaiqm,2) =  - fac*cosgt(i,ia)*ylm_red(i,l3*l3+m)
                         tmp(i,iaiqm,3) =    fac*cosgt(i,ia)*ylm_red(i,l3*l3+m)
                         tmp(i,iaiqm,4) =    fac*singt(i,ia)*ylm_red(i,l3*l3+m)
                      end do
                   end if
                end do
             end do
          end do
          call dgemm('N','T',ibl2-ibl1+1,nval,max_iaiqm,1.0d0,tmp(1,1,1),ibsize,yr,nval,0.0d0,rhogr_red,ibsize)
          call dgemm('N','T',ibl2-ibl1+1,nval,max_iaiqm,1.0d0,tmp(1,1,2),ibsize,yi,nval,1.0d0,rhogr_red,ibsize)
          call dgemm('N','T',ibl2-ibl1+1,nval,max_iaiqm,1.0d0,tmp(1,1,3),ibsize,yr,nval,0.0d0,rhogi_red,ibsize)
          call dgemm('N','T',ibl2-ibl1+1,nval,max_iaiqm,1.0d0,tmp(1,1,4),ibsize,yi,nval,1.0d0,rhogi_red,ibsize)
          do m_ = 1, nval
             if(occup_val(m_,kk) < DELTA) cycle
             do i = ibl1, ibl2
                rhogr(i,m_,jkbz) = rhogr(i,m_,jkbz) + rhogr_red(i-ibl1+1,m_)
                rhogi(i,m_,jkbz) = rhogi(i,m_,jkbz) + rhogi_red(i-ibl1+1,m_)
             end do
          end do
       end do
    end do

    deallocate(cosgt)
    deallocate(singt)
    deallocate(qitg_red)
    deallocate(ylm_red)
    deallocate(rhogr_red)
    deallocate(rhogi_red)
    deallocate(yr)
    deallocate(yi)
    deallocate(tmp)

    call tstatc0_end(id_sname)
  end subroutine add_RHOG_hard_part_
#endif
  subroutine add_RHOG_hard_part_rs(iqmk,rhor,rhoi,fnr,fni,fmr,fmi)
    implicit none
    integer, intent(in) :: iqmk
    real(kind=DP), intent(out) :: rhor(nfftp_exx_nonpara/2)
    real(kind=DP), intent(out) :: rhoi(nfftp_exx_nonpara/2)
    real(kind=DP), intent(in) :: fnr(nlmta)
    real(kind=DP), intent(in) :: fni(nlmta)
    real(kind=DP), intent(in) :: fmr(nlmta)
    real(kind=DP), intent(in) :: fmi(nlmta)
    integer :: n,it,ia,imesh,lmt1,lmt2,il1,il2,tau1,tau2,ilmta1,ilmta2,lmtp
    integer :: nma
    real(kind=DP), dimension(:), allocatable :: cosqmkr,sinqmkr
    real(kind=DP) :: rr,ii,si,co,qm,rr0,ii0
    integer :: ind
    real(kind=DP), allocatable, dimension(:) :: rtmp,itmp
    real(kind=DP) :: fac
    integer,save  :: id_sname = -1
    call tstatc0_begin('add_RHOG_hard_part_rs ',id_sname)
    allocate(cosqmkr(nmesh_rs_aug_max));cosqmkr=0.d0
    allocate(sinqmkr(nmesh_rs_aug_max));sinqmkr=0.d0
    allocate(rtmp(nmesh_rs_aug_max));rtmp=0.d0
    allocate(itmp(nmesh_rs_aug_max));itmp=0.d0
    rhor(:) = 0.d0;rhoi(:) = 0.d0
    do ia=1,natm
       it = ityp(ia)
       if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
       call qmk_dot_r(iqmk,ia,cosqmkr,sinqmkr)
       nma = nmesh_rs_aug(ia)
       rtmp=0.d0;itmp=0.d0
       do lmtp = 1,nlmtpair(ia)
          lmt1 = plmt1(lmtp,ia)
          lmt2 = plmt2(lmtp,ia)
          ilmta1 = lmta(lmt1,ia)
          ilmta2 = lmta(lmt2,ia)
          rr = fnr(ilmta1)*fmr(ilmta2) + fni(ilmta1)*fmi(ilmta2)
          ii = fnr(ilmta1)*fmi(ilmta2) - fni(ilmta1)*fmr(ilmta2)
          if(lmt1.ne.lmt2)then
            rr = rr+fnr(ilmta2)*fmr(ilmta1) + fni(ilmta2)*fmi(ilmta1)
            ii = ii+fnr(ilmta2)*fmi(ilmta1) - fni(ilmta2)*fmr(ilmta1)
          endif
          do imesh=1,nma
             co = cosqmkr(imesh)
             si = sinqmkr(imesh)
             qm = qr_clm_ylm(imesh,ia,lmtp)
             rtmp(imesh) = rtmp(imesh) + qm*(rr*co+ii*si)
             itmp(imesh) = itmp(imesh) + qm*(ii*co-rr*si)
          enddo
       enddo
       do imesh=1,nma
          ind = meshxyz_rs_aug(imesh,ia)
          rhor(ind) = rhor(ind) + rtmp(imesh)
          rhoi(ind) = rhoi(ind) + itmp(imesh)
       enddo
    enddo
    deallocate(cosqmkr)
    deallocate(sinqmkr)
    deallocate(rtmp)
    deallocate(itmp)
    call tstatc0_end(id_sname)
  end subroutine add_RHOG_hard_part_rs

  subroutine add_RHOG_hard_part_rs2(iqmk,rhor,rhoi,qmfnr,qmfni,fmr,fmi)
    integer, intent(in) :: iqmk
    real(kind=DP), intent(out) :: rhor(nfftp_exx_nonpara/2)
    real(kind=DP), intent(out) :: rhoi(nfftp_exx_nonpara/2)
    real(kind=DP), intent(in) :: qmfnr(nmesh_rs_aug_max,nlmta)
    real(kind=DP), intent(in) :: qmfni(nmesh_rs_aug_max,nlmta)
    real(kind=DP), intent(in) :: fmr(nlmta)
    real(kind=DP), intent(in) :: fmi(nlmta)
    real(kind=DP), dimension(:), allocatable :: cosqmkr,sinqmkr
    real(kind=DP), allocatable, dimension(:) :: rtmp,itmp
    real(kind=DP) :: fr,fi,co,si,qmr,qmi,rr,ii
    integer :: ia,it,lmt1,ilmta1,ind,imesh,nma
    integer,save  :: id_sname = -1
    call tstatc0_begin('add_RHOG_hard_part_rs2 ',id_sname)
    allocate(cosqmkr(nmesh_rs_aug_max));cosqmkr=0.d0
    allocate(sinqmkr(nmesh_rs_aug_max));sinqmkr=0.d0
    allocate(rtmp(nmesh_rs_aug_max));rtmp=0.d0
    allocate(itmp(nmesh_rs_aug_max));itmp=0.d0
    rhor(:) = 0.d0;rhoi(:) = 0.d0
    do ia=1,natm
       it = ityp(ia)
       if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
       call qmk_dot_r(iqmk,ia,cosqmkr,sinqmkr)
       nma = nmesh_rs_aug(ia)
       rtmp=0.d0;itmp=0.d0
       do lmt1=1,ilmt(it)
          ilmta1 = lmta(lmt1,ia)
          fr =  fmr(ilmta1)
          fi = -fmi(ilmta1) ! c.c. of fsr
          do imesh=1,nma
             co =  cosqmkr(imesh)
             si = -sinqmkr(imesh) ! c.c. of exp(i(G+q-k))
             qmr = qmfnr(imesh,ilmta1)
             qmi = qmfni(imesh,ilmta1)
             rr = co*qmr-si*qmi
             ii = co*qmi+si*qmr
             rtmp(imesh) = rtmp(imesh) + rr*fr-ii*fi
             itmp(imesh) = itmp(imesh) + ii*fr+rr*fi
          enddo
       enddo
       do imesh=1,nma
          ind = meshxyz_rs_aug(imesh,ia)
          rhor(ind) = rhor(ind) + rtmp(imesh)
          rhoi(ind) = rhoi(ind) + itmp(imesh)
       enddo
    enddo
    deallocate(cosqmkr)
    deallocate(sinqmkr)
    deallocate(rtmp)
    deallocate(itmp)
    call tstatc0_end(id_sname)
  end subroutine add_RHOG_hard_part_rs2

  subroutine add_RHOG_hard_part_rs3(iqmk,kk,rhor,rhoi,fnr,fni,fmr,fmi)
    implicit none
    integer, intent(in) :: iqmk,kk
    real(kind=DP), intent(out) :: rhor(nfftp_exx_nonpara/2,ista_nval:iend_nval)
    real(kind=DP), intent(out) :: rhoi(nfftp_exx_nonpara/2,ista_nval:iend_nval)
    real(kind=DP), intent(in) :: fnr(nlmta)
    real(kind=DP), intent(in) :: fni(nlmta)
    real(kind=DP), intent(in) :: fmr(ista_nval:iend_nval,nlmta)
    real(kind=DP), intent(in) :: fmi(ista_nval:iend_nval,nlmta)
    integer :: n,it,ia,imesh,lmt1,lmt2,il1,il2,tau1,tau2,ilmta1,ilmta2,lmtp,m
    integer :: nma
    real(kind=DP), dimension(:), allocatable :: cosqmkr,sinqmkr
    real(kind=DP) :: rr,ii,si,co,qm,rr0,ii0,ialmtpm
    real(kind=DP), allocatable, dimension(:,:) :: rrarray,iiarray
    integer :: ind
    real(kind=DP), allocatable, dimension(:) :: rtmp,itmp
    real(kind=DP) :: fac
    real(kind=DP), allocatable, dimension(:,:) :: rhor2,rhoi2
    real(kind=DP), allocatable, dimension(:,:) :: tmp1,tmp2
    integer,save  :: id_sname = -1,id_sname_d=-1
    call tstatc0_begin('add_RHOG_hard_part_rs ',id_sname)
    allocate(cosqmkr(nmesh_rs_aug_max));cosqmkr=0.d0
    allocate(sinqmkr(nmesh_rs_aug_max));sinqmkr=0.d0
    allocate(rhor2(nmesh_rs_aug_max,np_nval))
    allocate(rhoi2(nmesh_rs_aug_max,np_nval))
    allocate(rrarray(np_nval,maxval(nlmtpair)))
    allocate(iiarray(np_nval,maxval(nlmtpair)))
    allocate(tmp1(nmesh_rs_aug_max,maxval(nlmtpair)))
    allocate(tmp2(nmesh_rs_aug_max,maxval(nlmtpair)))
    rhor(:,:) = 0.d0;rhoi(:,:) = 0.d0
    do ia=1,natm
       it = ityp(ia)
       if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
       call qmk_dot_r(iqmk,ia,cosqmkr,sinqmkr)
       nma = nmesh_rs_aug(ia)
       rhor2 = 0.d0;rhoi2 = 0.d0
       rrarray = 0.d0;iiarray = 0.d0
       tmp1 = 0.d0;tmp2 = 0.d0
       do lmtp = 1,nlmtpair(ia)
          do m=ista_nval,iend_nval
             if(occup_val(m,kk) < DELTA) cycle
             lmt1 = plmt1(lmtp,ia)
             lmt2 = plmt2(lmtp,ia)
             ilmta1 = lmta(lmt1,ia)
             ilmta2 = lmta(lmt2,ia)
             rr = fnr(ilmta1)*fmr(m,ilmta2) + fni(ilmta1)*fmi(m,ilmta2)
             ii = fnr(ilmta1)*fmi(m,ilmta2) - fni(ilmta1)*fmr(m,ilmta2)
             if(lmt1.ne.lmt2)then
               rr = rr+fnr(ilmta2)*fmr(m,ilmta1) + fni(ilmta2)*fmi(m,ilmta1)
               ii = ii+fnr(ilmta2)*fmi(m,ilmta1) - fni(ilmta2)*fmr(m,ilmta1)
             endif
             rrarray(m-ista_nval+1,lmtp) = rr
             iiarray(m-ista_nval+1,lmtp) = ii
          enddo
          do imesh=1,nma
             co = cosqmkr(imesh)
             si = sinqmkr(imesh)
             qm = qr_clm_ylm(imesh,ia,lmtp)
             tmp1(imesh,lmtp) = qm*co
             tmp2(imesh,lmtp) = qm*si
          enddo
       enddo
       call dgemm('N','T',nmesh_rs_aug(ia),np_nval,nlmtpair(ia),1.d0,&
                 & tmp1,nmesh_rs_aug_max,rrarray,np_nval,0.d0,rhor2,nmesh_rs_aug_max)
       call dgemm('N','T',nmesh_rs_aug(ia),np_nval,nlmtpair(ia),1.d0,&
                 & tmp2,nmesh_rs_aug_max,iiarray,np_nval,1.d0,rhor2,nmesh_rs_aug_max)
       call dgemm('N','T',nmesh_rs_aug(ia),np_nval,nlmtpair(ia),1.d0,&
                 & tmp2,nmesh_rs_aug_max,rrarray,np_nval,0.d0,rhoi2,nmesh_rs_aug_max)
       call dgemm('N','T',nmesh_rs_aug(ia),np_nval,nlmtpair(ia),1.d0,&
                 & tmp1,nmesh_rs_aug_max,iiarray,np_nval,-1.d0,rhoi2,nmesh_rs_aug_max)
       do m=ista_nval,iend_nval
          if(occup_val(m,kk) < DELTA) cycle
          do imesh=1,nma
             ind = meshxyz_rs_aug(imesh,ia)
             rhor(ind,m) = rhor(ind,m) + rhor2(imesh,m-ista_nval+1)
             rhoi(ind,m) = rhoi(ind,m) + rhoi2(imesh,m-ista_nval+1)
          enddo
       enddo
    enddo

    deallocate(cosqmkr)
    deallocate(sinqmkr)
    deallocate(rhor2)
    deallocate(rhoi2)
    deallocate(rrarray)
    deallocate(iiarray)
    deallocate(tmp1)
    deallocate(tmp2)
    call tstatc0_end(id_sname)
  end subroutine add_RHOG_hard_part_rs3


  subroutine qmk_dot_r(iqmk,ia,zc_ar,zs_ar)
    integer, intent(in) :: iqmk,ia
    real(kind=DP),dimension(nmesh_rs_aug_max),intent(out) :: zc_ar,zs_ar
    integer :: i
    real(kind=DP) :: inl,inm,inn
    real(kind=DP) :: rx,ry,rz,kdr
    integer :: id_sname = -1
    inl = 1.d0/dble(fft_box_size_CD_exx(1,1))
    inm = 1.d0/dble(fft_box_size_CD_exx(2,1))
    inn = 1.d0/dble(fft_box_size_CD_exx(3,1))
    do i=1,nmesh_rs_aug(ia)
       rx = dble(meshx_rs_aug(i,ia))*inl
       ry = dble(meshy_rs_aug(i,ia))*inm
       rz = dble(meshz_rs_aug(i,ia))*inn
       kdr = (rx*qmk(iqmk,1)+ry*qmk(iqmk,2)+rz*qmk(iqmk,3))*PAI2
       zc_ar(i) = dcos(kdr)
       zs_ar(i) = dsin(kdr)
    enddo
  end subroutine qmk_dot_r

  subroutine integrate_QijVnm_2D(iqmk,potr,poti,fmr,fmi,qvr,qvi,dfmr,dfmi,dqvr,dqvi,gqvr,gqvi)
    implicit none
    integer, intent(in) :: iqmk
    real(kind=DP), intent(in) :: potr(nmax_G_hyb)
    real(kind=DP), intent(in) :: poti(nmax_G_hyb)
    real(kind=DP), intent(in) :: fmr(nlmta)
    real(kind=DP), intent(in) :: fmi(nlmta)
    real(kind=DP), intent(out) :: qvr(nlmta)
    real(kind=DP), intent(out) :: qvi(nlmta)
    real(kind=DP), optional, intent(in) :: dfmr(nlmta,3)
    real(kind=DP), optional, intent(in) :: dfmi(nlmta,3)
    real(kind=DP), optional, intent(out) :: dqvr(nlmta,3)
    real(kind=DP), optional, intent(out) :: dqvi(nlmta,3)
    real(kind=DP), optional, intent(out) :: gqvr(nlmta,3)
    real(kind=DP), optional, intent(out) :: gqvi(nlmta,3)

    integer :: it,ia,i,lmt1,lmt2,np,il1,il2,tau1,tau2,ilmta1,ilmta2
    integer :: l3,ilm3,iiqitg,iil
    real(kind=DP) :: ph, hr, hi, fac, f0
    real(kind=DP) :: cosgt, singt
    real(kind=DP) :: dr, di
    real(kind=DP), allocatable :: zsr(:), zsi(:)

    real(kind=DP), allocatable :: qitg_red(:,:),ylm_red(:,:)
    integer :: ibl1,ibl2,iq,inn,ip,m
    logical :: force_mode = .false.
    real(kind=DP), allocatable :: gvec(:,:)
    real(kind=DP) :: er(3), ei(3), vec1(3)
    integer,save  :: id_sname = -1

    call tstatc0_begin('integrate_QijVnm_2D ',id_sname)

    force_mode = present(dfmr).and.present(dfmi).and.present(dqvr).and.present(dqvi).and. &
               & present(gqvr).and.present(gqvi)
    if(.not.preproc_done) call hard_part_preproc()

    qvr(1:nlmta) = 0.d0;  qvi(1:nlmta) = 0.d0
    if(force_mode)then
       dqvr(1:nlmta,1:3) = 0.d0
       dqvi(1:nlmta,1:3) = 0.d0
       gqvr(1:nlmta,1:3) = 0.d0
       gqvi(1:nlmta,1:3) = 0.d0
    endif
    allocate(zsr(ibsize)); allocate(zsi(ibsize))
    allocate(qitg_red(ibsize,nqitg))
    allocate(ylm_red(ibsize,n*n))
    if(force_mode) allocate(gvec(ibsize,3))
    do ibl1=1, nmax_G_hyb, ibsize
       ibl2=ibl1+ibsize-1
       if(ibl2.gt.nmax_G_hyb) ibl2=nmax_G_hyb
       do iq=1,nqitg
          do i=1,ibl2-ibl1+1
             qitg_red(i,iq) = qitg_exx(i+ibl1-1,iq,iqmk)
          enddo
       enddo
       do inn=1,n*n
          do i=1,ibl2-ibl1+1
             ylm_red(i,inn) = ylm_exx(i+ibl1-1,inn,iqmk)
          enddo
       enddo

       if(force_mode)then
          do i=1,ibl2-ibl1+1
             gvec(i,1:3) = ngabc_red(i+ibl1-1,1:3) + qmk(iqmk,1:3)
          end do
#ifdef DERIVATIVE_IN_CARTS
          do i=1,ibl2-ibl1+1
             vec1(1) = rltv(1,1)*gvec(i,1) +rltv(1,2)*gvec(i,2) +rltv(1,3)*gvec(i,3)
             vec1(2) = rltv(2,1)*gvec(i,1) +rltv(2,2)*gvec(i,2) +rltv(2,3)*gvec(i,3)
             vec1(3) = rltv(3,1)*gvec(i,1) +rltv(3,2)*gvec(i,2) +rltv(3,3)*gvec(i,3)
             gvec(i,1:3) = vec1(1:3)
          end do
#endif
       endif
       do ia=1,natm
          it = ityp(ia)
          if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
          do i=1,ibl2-ibl1+1
             ph = PAI2 * dot_product(pos(ia,1:3),ngabc_red(i+ibl1-1,1:3)+qmk(iqmk,1:3))
             cosgt = cos(ph); singt = sin(ph)
             zsr(i) = cosgt * potr(i+ibl1-1) - singt * poti(i+ibl1-1)
             zsi(i) =-cosgt * poti(i+ibl1-1) - singt * potr(i+ibl1-1)
          end do
          do iq=nqitg_sp0(it), nqitg_sp(it)
             l3 = iq2l3(iq)
             iil = mod(l3,4)         ! (-i)^L

             do m = 1, 2*l3+1
                ilm3 = l3*l3+m

                dr = 0.d0; di = 0.d0
                if(force_mode)then
                   er(1:3) = 0.d0;ei(1:3)=0.d0
                endif
                if(iil==0) then
                   do i=1,ibl2-ibl1+1
                      fac =  qitg_red(i,iq) * ylm_red(i,ilm3)
                      dr = dr + fac * zsr(i)
                      di = di + fac * zsi(i)
                   end do
                   if(force_mode)then
                      do i=1,ibl2-ibl1+1
                         fac =  qitg_red(i,iq) * ylm_red(i,ilm3)
                         er = er + (fac * zsr(i))*gvec(i,1:3)
                         ei = ei + (fac * zsi(i))*gvec(i,1:3)
                      end do
                   endif
                else if(iil==1) then
                   do i=1,ibl2-ibl1+1
                      fac =  qitg_red(i,iq) * ylm_red(i,ilm3)
! ==== ASMS ====
                      dr = dr + fac * zsi(i)
                      di = di - fac * zsr(i)
! ==== ASMS ====
                   end do
                   if(force_mode)then
                      do i=1,ibl2-ibl1+1
                         fac =  qitg_red(i,iq) * ylm_red(i,ilm3)
!                         er = er - (fac * zsi(i))*gvec(i,1:3)
!                         ei = ei + (fac * zsr(i))*gvec(i,1:3)
                         er = er + (fac * zsi(i))*gvec(i,1:3)
                         ei = ei - (fac * zsr(i))*gvec(i,1:3)
                      end do
                   endif
                else if(iil==2) then
                   do i=1,ibl2-ibl1+1
                      fac =  qitg_red(i,iq) * ylm_red(i,ilm3)
                      dr = dr - fac * zsr(i)
                      di = di - fac * zsi(i)
                   end do
                   if(force_mode)then
                      do i=1,ibl2-ibl1+1
                         fac =  qitg_red(i,iq) * ylm_red(i,ilm3)
!                         er = er + (fac * zsi(i))*gvec(i,1:3)
!                         ei = ei - (fac * zsr(i))*gvec(i,1:3)
                         er = er - (fac * zsr(i))*gvec(i,1:3)
                         ei = ei - (fac * zsi(i))*gvec(i,1:3)
                      end do
                   endif
                else if(iil==3) then
                   do i=1,ibl2-ibl1+1
                      fac =  qitg_red(i,iq) * ylm_red(i,ilm3)
! ==== ASMS ====
                      dr = dr - fac * zsi(i)
                      di = di + fac * zsr(i)
! ==== ASMS ====
                   end do
                   if(force_mode)then
                      do i=1,ibl2-ibl1+1
                         fac =  qitg_red(i,iq) * ylm_red(i,ilm3)
! ==== ASMS ====
                         er = er - (fac * zsi(i))*gvec(i,1:3)
                         ei = ei + (fac * zsr(i))*gvec(i,1:3)
! ==== ASMS ====
                      end do
                   endif
                end if

! ==== ASMS ==== 2017/09/22
!                di = -di
!                if(force_mode) ei=-ei
! ==== ASMS ==== 2017/09/22

                do ip=1,nc(m,iq)
                   lmt1 = nc2lmt1(ip,m,iq); lmt2 = nc2lmt2(ip,m,iq)
                   ilmta1 = lmta(lmt1,ia);  ilmta2 = lmta(lmt2,ia)
                   np = nc2n(ip,m,iq)
                   f0 = iwei(ia) * dl2p(lmt1,lmt2,np,it)
                   qvr(ilmta1) = qvr(ilmta1) + f0 * (dr * fmr(ilmta2) - di * fmi(ilmta2))
                   qvi(ilmta1) = qvi(ilmta1) + f0 * (di * fmr(ilmta2) + dr * fmi(ilmta2))
                end do
                if(force_mode)then
                   do ip=1,nc(m,iq)
                      lmt1 = nc2lmt1(ip,m,iq); lmt2 = nc2lmt2(ip,m,iq)
                      ilmta1 = lmta(lmt1,ia);  ilmta2 = lmta(lmt2,ia)
                      np = nc2n(ip,m,iq)
                      f0 = iwei(ia) * dl2p(lmt1,lmt2,np,it)
                      dqvr(ilmta1,1:3) = dqvr(ilmta1,1:3) + f0 * (dr * dfmr(ilmta2,1:3) - di * dfmi(ilmta2,1:3))
                      dqvi(ilmta1,1:3) = dqvi(ilmta1,1:3) + f0 * (di * dfmr(ilmta2,1:3) + dr * dfmi(ilmta2,1:3))
                      gqvr(ilmta1,1:3) = gqvr(ilmta1,1:3) + f0 * (ei(1:3) * fmr(ilmta2) + er(1:3) * fmi(ilmta2))
                      gqvi(ilmta1,1:3) = gqvi(ilmta1,1:3) + f0 * (-er(1:3) * fmr(ilmta2) + ei(1:3) * fmi(ilmta2))
                   enddo
                endif
             end do
          end do
       end do

    enddo

    deallocate(zsr);  deallocate(zsi)
    deallocate(qitg_red); deallocate(ylm_red)
    if(force_mode) deallocate(gvec)
    call tstatc0_end(id_sname)
  end subroutine integrate_QijVnm_2D

  subroutine add_Vx_hard_part_2D(ik,vxw,sumqvr,sumqvi)
  !!subroutine add_Vx_hard_part(ik,vxw1,sumqvr,sumqvi)
    implicit none
    integer, intent(in) :: ik
    real(kind=DP), intent(inout) :: vxw(kg1,kimg)
    !!real(kind=DP), intent(inout) :: vxw1(kg1,kimg)
    real(kind=DP), intent(in) :: sumqvr(nlmta)
    real(kind=DP), intent(in) :: sumqvi(nlmta)

    integer :: it,ia,i,il1,ilmtt1,ilmta1,iil,iksnl,ikk
    integer :: lmt1,lmtt1
    real(kind=DP) :: fac, ph
    real(kind=DP), allocatable :: cosgt(:) ! d(kg1)
    real(kind=DP), allocatable :: singt(:) ! d(kg1)
    real(kind=DP), allocatable, dimension(:,:) :: snl_t
    integer,save  :: id_sname = -1
    call tstatc0_begin('add_Vx_hard_part ',id_sname)

    allocate(cosgt(kg1))
    allocate(singt(kg1))

    iksnl = (ik-1)/nspin + 1
    allocate(snl_t(kg1,nlmtt));snl_t=0.d0
    do i=ista_g1k(ik),iend_g1k(iK)
       snl_t(i,1:nlmtt) = snl(i-ista_g1k(ik)+1,1:nlmtt,iksnl)
    enddo
    call mpi_allreduce(MPI_IN_PLACE,snl_t,kg1*nlmtt,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
    ikk = k_index(ik)
    do ia=1,natm
       it = ityp(ia)
       if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
       do i=1,iba(ik)
          ph = PAI2 * dot_product(pos(ia,1:3),vkxyz_for_outerloop(ikk,1:3,BUCS)+ngabc(nbase(i,ik),1:3))
          cosgt(i) = cos(ph)
          singt(i) = sin(ph)
       end do
       do lmt1 = 1, ilmt(it)
          il1  = ltp(lmt1,it)
          lmtt1 = lmtt(lmt1,it)
          ilmta1 = lmta(lmt1,ia)
          iil=mod(il1-1,4)
          if(iil==0) then
             do i=1,iba(ik)
                fac = iwei(ia)*snl_t(i,lmtt1)
                vxw(i,1) = vxw(i,1) + fac * (cosgt(i)*sumqvr(ilmta1)+singt(i)*sumqvi(ilmta1))
                if(kimg==2) vxw(i,2) = vxw(i,2) + fac * (cosgt(i)*sumqvi(ilmta1)-singt(i)*sumqvr(ilmta1))
             end do
          else if(iil==1) then
             do i=1,iba(ik)
                fac = iwei(ia)*snl_t(i,lmtt1)
                vxw(i,1) = vxw(i,1) + fac * (cosgt(i)*sumqvi(ilmta1)-singt(i)*sumqvr(ilmta1))
                if(kimg==2) vxw(i,2) = vxw(i,2) - fac * (cosgt(i)*sumqvr(ilmta1)+singt(i)*sumqvi(ilmta1))
             end do
          else if(iil==2) then
             do i=1,iba(ik)
                fac = -iwei(ia)*snl_t(i,lmtt1)
                vxw(i,1) = vxw(i,1) + fac * (cosgt(i)*sumqvr(ilmta1)+singt(i)*sumqvi(ilmta1))
                if(kimg==2) vxw(i,2) = vxw(i,2) + fac * (cosgt(i)*sumqvi(ilmta1)-singt(i)*sumqvr(ilmta1))
             end do
          else if(iil==3) then
             do i=1,iba(ik)
                fac = -iwei(ia)*snl_t(i,lmtt1)
                vxw(i,1) = vxw(i,1) + fac * (cosgt(i)*sumqvi(ilmta1)-singt(i)*sumqvr(ilmta1))
                if(kimg==2) vxw(i,2) = vxw(i,2) - fac * (cosgt(i)*sumqvr(ilmta1)+singt(i)*sumqvi(ilmta1))
             end do
          end if
       end do
    end do

    deallocate(cosgt)
    deallocate(singt)

    deallocate(snl_t)
    call tstatc0_end(id_sname)
  end subroutine add_Vx_hard_part_2D

! ================================== KT_Test ========================= 12.5Exp
  subroutine integrate_QijVnm(iqmk,potr,poti,fmr,fmi,qvr,qvi,dfmr,dfmi,dqvr,dqvi,gqvr,gqvi)
    implicit none
    integer, intent(in) :: iqmk
    real(kind=DP), intent(in) :: potr(ista_kngp:iend_kngp)
    real(kind=DP), intent(in) :: poti(ista_kngp:iend_kngp)
    real(kind=DP), intent(in) :: fmr(nlmta)
    real(kind=DP), intent(in) :: fmi(nlmta)
    real(kind=DP), intent(out) :: qvr(nlmta)
    real(kind=DP), intent(out) :: qvi(nlmta)
    real(kind=DP), optional, intent(in) :: dfmr(nlmta,3)
    real(kind=DP), optional, intent(in) :: dfmi(nlmta,3)
    real(kind=DP), optional, intent(out) :: dqvr(nlmta,3)
    real(kind=DP), optional, intent(out) :: dqvi(nlmta,3)
    real(kind=DP), optional, intent(out) :: gqvr(nlmta,3)
    real(kind=DP), optional, intent(out) :: gqvi(nlmta,3)

    integer :: it,ia,i,lmt1,lmt2,np,il1,il2,tau1,tau2,ilmta1,ilmta2
    integer :: l3,ilm3,iiqitg,iil
    real(kind=DP) :: ph, hr, hi, fac, f0
    real(kind=DP) :: cosgt, singt
    real(kind=DP) :: dr, di
    real(kind=DP), allocatable :: zsr(:), zsi(:)

    real(kind=DP), allocatable :: qitg_red(:,:),ylm_red(:,:)
    integer :: ibl1,ibl2,iq,inn,ip,m, ips, ipe
    logical :: force_mode = .false.
    real(kind=DP), allocatable :: gvec(:,:)
    real(kind=DP) :: er(3), ei(3), vec1(3)
    integer,save  :: id_sname = -1

    call tstatc0_begin('integrate_QijVnm ',id_sname)
    force_mode = present(dfmr).and.present(dfmi).and.present(dqvr).and.present(dqvi).and. &
               & present(gqvr).and.present(gqvi)
    if(.not.preproc_done) call hard_part_preproc()

    qvr(1:nlmta) = 0.d0;  qvi(1:nlmta) = 0.d0
    if(force_mode)then
       dqvr(1:nlmta,1:3) = 0.d0
       dqvi(1:nlmta,1:3) = 0.d0
       gqvr(1:nlmta,1:3) = 0.d0
       gqvi(1:nlmta,1:3) = 0.d0
    endif
    allocate(zsr(ibsize)); allocate(zsi(ibsize))
    allocate(qitg_red(ibsize,nqitg))
    allocate(ylm_red(ibsize,n*n))
    if(force_mode) allocate(gvec(ibsize,3))

    ips = 1
    ipe = nmax_G_hyb
    ips = ista_kngp
    ipe = iend_kngp

    do ibl1=ips, ipe, ibsize
       ibl2=min(ipe,ibl1+ibsize-1)
       do iq=1,nqitg
          do i=1,ibl2-ibl1+1
             qitg_red(i,iq) = qitg_exx(i+ibl1-1,iq,iqmk)
          enddo
       enddo
       do inn=1,n*n
          do i=1,ibl2-ibl1+1
             ylm_red(i,inn) = ylm_exx(i+ibl1-1,inn,iqmk)
          enddo
       enddo

       if(force_mode)then
          do i=1,ibl2-ibl1+1
! === DEBUG by tkato 2015/03/19 ================================================
!            gvec(i,1:3) = ngabc(i+ibl1-1,1:3) + qmk(iqmk,1:3)
             gvec(i,1:3) = ngabc_kngp_l(i+ibl1-1,1:3) + qmk(iqmk,1:3)
! ==============================================================================
          end do
#ifdef DERIVATIVE_IN_CARTS
          do i=1,ibl2-ibl1+1
             vec1(1) = rltv(1,1)*gvec(i,1) +rltv(1,2)*gvec(i,2) +rltv(1,3)*gvec(i,3)
             vec1(2) = rltv(2,1)*gvec(i,1) +rltv(2,2)*gvec(i,2) +rltv(2,3)*gvec(i,3)
             vec1(3) = rltv(3,1)*gvec(i,1) +rltv(3,2)*gvec(i,2) +rltv(3,3)*gvec(i,3)
             gvec(i,1:3) = vec1(1:3)
          end do
#endif
       endif
       do ia=1,natm
          it = ityp(ia)
          if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
          do i=1,ibl2-ibl1+1
! === DEBUG by tkato 2015/03/19 ================================================
!            ph = PAI2 * dot_product(pos(ia,1:3),ngabc(i+ibl1-1,1:3)+qmk(iqmk,1:3))
             ph = PAI2 * dot_product(pos(ia,1:3),ngabc_kngp_l(i+ibl1-1,1:3)+qmk(iqmk,1:3))
! ==============================================================================
             cosgt = cos(ph); singt = sin(ph)
! ==== ASMS ====
             zsr(i) = cosgt * potr(i+ibl1-1) - singt * poti(i+ibl1-1)
             zsi(i) =-cosgt * poti(i+ibl1-1) - singt * potr(i+ibl1-1)
! ==== ASMS ====
          end do

          do iq=nqitg_sp0(it), nqitg_sp(it)
             l3 = iq2l3(iq)
             iil = mod(l3,4)         ! (-i)^L

             do m = 1, 2*l3+1
                ilm3 = l3*l3+m

                dr = 0.d0; di = 0.d0
                if(force_mode)then
                   er(1:3) = 0.d0;ei(1:3)=0.d0
                endif
                if(iil==0) then
                   do i=1,ibl2-ibl1+1
                      fac =  qitg_red(i,iq) * ylm_red(i,ilm3)
                      dr = dr + fac * zsr(i)
                      di = di + fac * zsi(i)
                   end do
                   if(force_mode)then
                      do i=1,ibl2-ibl1+1
                         fac =  qitg_red(i,iq) * ylm_red(i,ilm3)
                         er = er + (fac * zsr(i))*gvec(i,1:3)
                         ei = ei + (fac * zsi(i))*gvec(i,1:3)
                      end do
                   endif
                else if(iil==1) then
                   do i=1,ibl2-ibl1+1
                      fac =  qitg_red(i,iq) * ylm_red(i,ilm3)
! ==== ASMS ====
                      dr = dr + fac * zsi(i)
                      di = di - fac * zsr(i)
! ==== ASMS ====
                   end do
                   if(force_mode)then
                      do i=1,ibl2-ibl1+1
                         fac =  qitg_red(i,iq) * ylm_red(i,ilm3)
! ==== ASMS ====
                         er = er + (fac * zsi(i))*gvec(i,1:3)
                         ei = ei - (fac * zsr(i))*gvec(i,1:3)
! ==== ASMS ====
                      end do
                   endif
                else if(iil==2) then
                   do i=1,ibl2-ibl1+1
                      fac =  qitg_red(i,iq) * ylm_red(i,ilm3)
                      dr = dr - fac * zsr(i)
                      di = di - fac * zsi(i)
                   end do
                   if(force_mode)then
                      do i=1,ibl2-ibl1+1
                         fac =  qitg_red(i,iq) * ylm_red(i,ilm3)
                         er = er - (fac * zsr(i))*gvec(i,1:3)
                         ei = ei - (fac * zsi(i))*gvec(i,1:3)
                      end do
                   endif
                else if(iil==3) then
                   do i=1,ibl2-ibl1+1
                      fac =  qitg_red(i,iq) * ylm_red(i,ilm3)
! ==== ASMS ====
                      dr = dr - fac * zsi(i)
                      di = di + fac * zsr(i)
! ==== ASMS ====
                   end do
                   if(force_mode)then
                      do i=1,ibl2-ibl1+1
                         fac =  qitg_red(i,iq) * ylm_red(i,ilm3)
! ==== ASMS ====
                         er = er - (fac * zsi(i))*gvec(i,1:3)
                         ei = ei + (fac * zsr(i))*gvec(i,1:3)
! ==== ASMS ====
                      end do
                   endif
                end if

! ==== ASMS ==== 2017/09/22
!                di = -di
!                if(force_mode) ei=-ei
! ==== ASMS ==== 2017/09/22

                do ip=1,nc(m,iq)
                   lmt1 = nc2lmt1(ip,m,iq); lmt2 = nc2lmt2(ip,m,iq)
                   ilmta1 = lmta(lmt1,ia);  ilmta2 = lmta(lmt2,ia)
                   np = nc2n(ip,m,iq)
                   f0 = iwei(ia) * dl2p(lmt1,lmt2,np,it)
                   qvr(ilmta1) = qvr(ilmta1) + f0 * (dr * fmr(ilmta2) - di * fmi(ilmta2))
                   qvi(ilmta1) = qvi(ilmta1) + f0 * (di * fmr(ilmta2) + dr * fmi(ilmta2))
                end do
                if(force_mode)then
                   do ip=1,nc(m,iq)
                      lmt1 = nc2lmt1(ip,m,iq); lmt2 = nc2lmt2(ip,m,iq)
                      ilmta1 = lmta(lmt1,ia);  ilmta2 = lmta(lmt2,ia)
                      np = nc2n(ip,m,iq)
                      f0 = iwei(ia) * dl2p(lmt1,lmt2,np,it)
                      dqvr(ilmta1,1:3) = dqvr(ilmta1,1:3) + f0 * (dr * dfmr(ilmta2,1:3) - di * dfmi(ilmta2,1:3))
                      dqvi(ilmta1,1:3) = dqvi(ilmta1,1:3) + f0 * (di * dfmr(ilmta2,1:3) + dr * dfmi(ilmta2,1:3))
                      gqvr(ilmta1,1:3) = gqvr(ilmta1,1:3) + f0 * (ei(1:3) * fmr(ilmta2) + er(1:3) * fmi(ilmta2))
                      gqvi(ilmta1,1:3) = gqvi(ilmta1,1:3) + f0 * (-er(1:3) * fmr(ilmta2) + ei(1:3) * fmi(ilmta2))
                   enddo
                endif
             end do
          end do
       end do

    enddo

    deallocate(zsr);  deallocate(zsi)
    deallocate(qitg_red); deallocate(ylm_red)
    if(force_mode) deallocate(gvec)
    call tstatc0_end(id_sname)
  end subroutine integrate_QijVnm
#ifdef HYBRID_DGEMM
  subroutine integrate_QijVnm_(iqmk_,potr,poti,fmr,fmi,qvr,qvi,ik,ispin,dfmr,dfmi,dqvr,dqvi,gqvr,gqvi)
    implicit none
    integer, intent(in) :: iqmk_(kv3bz,kv3), ik, ispin
    real(kind=DP), intent(in) :: potr(ista_kngp:iend_kngp, nval, kv3bz)
    real(kind=DP), intent(in) :: poti(ista_kngp:iend_kngp, nval, kv3bz)
    real(kind=DP), intent(in) :: fmr(nlmta, nval, kv3bz)
    real(kind=DP), intent(in) :: fmi(nlmta, nval, kv3bz)
    real(kind=DP), intent(out) :: qvr(nlmta, nval, kv3bz)
    real(kind=DP), intent(out) :: qvi(nlmta, nval, kv3bz)
    real(kind=DP), optional, intent(in) :: dfmr(nlmta,3)
    real(kind=DP), optional, intent(in) :: dfmi(nlmta,3)
    real(kind=DP), optional, intent(out) :: dqvr(nlmta,nval,kv3bz,3)
    real(kind=DP), optional, intent(out) :: dqvi(nlmta,nval,kv3bz,3)
    real(kind=DP), optional, intent(out) :: gqvr(nlmta,nval,kv3bz,3)
    real(kind=DP), optional, intent(out) :: gqvi(nlmta,nval,kv3bz,3)

    integer :: it,ia,i,lmt1,lmt2,np,il1,il2,tau1,tau2,ilmta1,ilmta2
    integer :: l3,ilm3,iiqitg,iil
    real(kind=DP) :: ph, hr, hi,  f0, vec1(3)
    real(kind=DP), allocatable :: cosgt(:,:), singt(:,:)
    real(kind=DP), allocatable :: dr(:,:), di(:,:)
    real(kind=DP), allocatable :: zsr(:), zsi(:)
    real(kind=DP), allocatable :: tmp(:,:,:)

    real(kind=DP), allocatable :: qitg_red(:,:),ylm_red(:,:),fac(:)
    integer :: ibl1,ibl2,iq,inn,ip,m
    logical :: force_mode = .false.
    real(kind=DP), allocatable :: tmp2(:,:,:,:)
    real(kind=DP), allocatable :: gvec(:,:)
    real(kind=DP), allocatable :: er(:,:,:), ei(:,:,:)
    integer,save  :: id_sname = -1
    integer :: jkbz, jtrs, m_, max_iaiqm, iaiqm, kngpsize, kk, iqmk
    kngpsize = iend_kngp - ista_kngp +1

    call tstatc0_begin('integrate_QijVnm_',id_sname)
    force_mode = present(dfmr).and.present(dfmi).and.present(dqvr).and.present(dqvi).and. &
               & present(gqvr).and.present(gqvi)
    if(.not.preproc_done) call hard_part_preproc()

    qvr = 0.d0;  qvi = 0.d0
    if(force_mode)then
       dqvr = 0.d0
       dqvi = 0.d0
       gqvr = 0.d0
       gqvi = 0.d0
    endif
    allocate(cosgt(ibsize,natm))
    allocate(singt(ibsize,natm))
    allocate(zsr(ibsize)); allocate(zsi(ibsize))
    allocate(qitg_red(ibsize,nqitg))
    allocate(ylm_red(ibsize,n*n))
    allocate(fac(ibsize))

    max_iaiqm = 0
    do ia = 1, natm
       it = ityp(ia)
       if(m_PP_include_vanderbilt_pot(it) == SKIP) cycle
       do iq = nqitg_sp0(it), nqitg_sp(it)
          l3 = iq2l3(iq)
          do m = 1, 2*l3+1
             max_iaiqm = max_iaiqm + 1
          end do
       end do
    end do
    allocate(dr(max_iaiqm,nval))
    allocate(di(max_iaiqm,nval))
    allocate(tmp(max_iaiqm,ibsize,4))

    if(force_mode) then
       allocate(tmp2(max_iaiqm,ibsize,3,4))
       allocate(er(max_iaiqm,nval,3))
       allocate(ei(max_iaiqm,nval,3))
       allocate(gvec(ibsize,3))
    endif

    do ibl1=ista_kngp,iend_kngp,ibsize
       ibl2=ibl1+ibsize-1
       if(ibl2.gt.iend_kngp) ibl2=iend_kngp
       do jkbz=1,kv3bz
          if(.not.q_on_k_centered_mesh(jkbz,ik)) cycle
          kk = ikp(jkbz) + ispin - 1
          iqmk = iqmk_(jkbz,ik)
          do iq=1,nqitg
             do i=1,ibl2-ibl1+1
                qitg_red(i,iq) = qitg_exx(i+ibl1-1,iq,iqmk)
             enddo
          enddo
          do inn=1,n*n
             do i=1,ibl2-ibl1+1
                ylm_red(i,inn) = ylm_exx(i+ibl1-1,inn,iqmk)
             enddo
          enddo

          if(force_mode)then
             do i=1,ibl2-ibl1+1
! === DEBUG by tkato 2015/03/19 ================================================
!               gvec(i,1:3) = ngabc(i+ibl1-1,1:3) + qmk(iqmk,1:3)
                gvec(i,1:3) = ngabc_kngp_l(i+ibl1-1,1:3) + qmk(iqmk,1:3)
!===============================================================================
             end do
#ifdef DERIVATIVE_IN_CARTS
             do i=1,ibl2-ibl1+1
                vec1(1) = rltv(1,1)*gvec(i,1) +rltv(1,2)*gvec(i,2) +rltv(1,3)*gvec(i,3)
                vec1(2) = rltv(2,1)*gvec(i,1) +rltv(2,2)*gvec(i,2) +rltv(2,3)*gvec(i,3)
                vec1(3) = rltv(3,1)*gvec(i,1) +rltv(3,2)*gvec(i,2) +rltv(3,3)*gvec(i,3)
                gvec(i,1:3) = vec1(1:3)
             end do
#endif
          endif

          do ia=1,natm
             it = ityp(ia)
             if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
             do i=1,ibl2-ibl1+1
! === DEBUG by tkato 2015/03/19 ================================================
!               ph = PAI2 * dot_product(pos(ia,1:3),ngabc(i+ibl1-1,1:3)+qmk(iqmk,1:3))
                ph = PAI2 * dot_product(pos(ia,1:3),ngabc_kngp_l(i+ibl1-1,1:3)+qmk(iqmk,1:3))
!===============================================================================
                cosgt(i,ia) = cos(ph); singt(i,ia) = sin(ph)
             end do
          end do
          iaiqm = 0
          dr = 0.d0; di = 0.d0
          if(force_mode)then
             er = 0.d0; ei = 0.d0
          endif

          do ia=1,natm
             it = ityp(ia)
             if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
             do iq=nqitg_sp0(it), nqitg_sp(it)
                l3 = iq2l3(iq)
                iil = mod(l3,4)

                do m = 1, 2*l3+1
                   iaiqm = iaiqm + 1
                   ilm3 = l3*l3+m
                   do i=1,ibl2-ibl1+1
                      fac(i) =  qitg_red(i,iq) * ylm_red(i,ilm3)
                   enddo

                   if(iil==0) then
                      do i=1,ibl2-ibl1+1
                         tmp(iaiqm,i,1) =  fac(i) * cosgt(i,ia)
                         tmp(iaiqm,i,2) = -fac(i) * singt(i,ia)
                         tmp(iaiqm,i,3) = -fac(i) * singt(i,ia)
                         tmp(iaiqm,i,4) = -fac(i) * cosgt(i,ia)
                      end do
                      if(force_mode)then
                         do i=1,ibl2-ibl1+1
                            tmp2(iaiqm,i,1:3,1) =  fac(i) * cosgt(i,ia) * gvec(i,1:3)
                            tmp2(iaiqm,i,1:3,2) = -fac(i) * singt(i,ia) * gvec(i,1:3)
                            tmp2(iaiqm,i,1:3,3) = -fac(i) * singt(i,ia) * gvec(i,1:3)
                            tmp2(iaiqm,i,1:3,4) = -fac(i) * cosgt(i,ia) * gvec(i,1:3)
                         end do
                      endif
                  else if(iil==1) then
                      do i=1,ibl2-ibl1+1
                         tmp(iaiqm,i,1) = -fac(i) * singt(i,ia)
                         tmp(iaiqm,i,2) = -fac(i) * cosgt(i,ia)
                         tmp(iaiqm,i,3) = -fac(i) * cosgt(i,ia)
                         tmp(iaiqm,i,4) =  fac(i) * singt(i,ia)
                      end do
                      if(force_mode)then
                         do i=1,ibl2-ibl1+1
                            tmp2(iaiqm,i,1:3,1) = -fac(i) * singt(i,ia) * gvec(i,1:3)
                            tmp2(iaiqm,i,1:3,2) = -fac(i) * cosgt(i,ia) * gvec(i,1:3)
                            tmp2(iaiqm,i,1:3,3) = -fac(i) * cosgt(i,ia) * gvec(i,1:3)
                            tmp2(iaiqm,i,1:3,4) =  fac(i) * singt(i,ia) * gvec(i,1:3)
                         end do
                      endif
                   else if(iil==2) then
                      do i=1,ibl2-ibl1+1
                         tmp(iaiqm,i,1) = -fac(i) * cosgt(i,ia)
                         tmp(iaiqm,i,2) =  fac(i) * singt(i,ia)
                         tmp(iaiqm,i,3) =  fac(i) * singt(i,ia)
                         tmp(iaiqm,i,4) =  fac(i) * cosgt(i,ia)
                      end do
                      if(force_mode)then
                         do i=1,ibl2-ibl1+1
                            tmp2(iaiqm,i,1:3,1) = -fac(i) * cosgt(i,ia) * gvec(i,1:3)
                            tmp2(iaiqm,i,1:3,2) =  fac(i) * singt(i,ia) * gvec(i,1:3)
                            tmp2(iaiqm,i,1:3,3) =  fac(i) * singt(i,ia) * gvec(i,1:3)
                            tmp2(iaiqm,i,1:3,4) =  fac(i) * cosgt(i,ia) * gvec(i,1:3)
                         end do
                      endif
                   else if(iil==3) then
                      do i=1,ibl2-ibl1+1
                         tmp(iaiqm,i,1) =  fac(i) * singt(i,ia)
                         tmp(iaiqm,i,2) =  fac(i) * cosgt(i,ia)
                         tmp(iaiqm,i,3) = -fac(i) * cosgt(i,ia)
                         tmp(iaiqm,i,4) =  fac(i) * singt(i,ia)
                      end do
                      if(force_mode)then
                         do i=1,ibl2-ibl1+1
                            tmp2(iaiqm,i,1:3,1) =  fac(i) * singt(i,ia) * gvec(i,1:3)
                            tmp2(iaiqm,i,1:3,2) =  fac(i) * cosgt(i,ia) * gvec(i,1:3)
                            tmp2(iaiqm,i,1:3,3) = -fac(i) * cosgt(i,ia) * gvec(i,1:3)
                            tmp2(iaiqm,i,1:3,4) =  fac(i) * singt(i,ia) * gvec(i,1:3)
                         end do
                      endif
                   end if
                end do
             end do
          end do
          call dgemm('N','N',max_iaiqm,nval,ibl2-ibl1+1,1.0d0,tmp(1,1,1),max_iaiqm,potr(ibl1,1,jkbz),&
             & kngpsize,0.0d0,dr(1,1),max_iaiqm)
          call dgemm('N','N',max_iaiqm,nval,ibl2-ibl1+1,1.0d0,tmp(1,1,2),max_iaiqm,poti(ibl1,1,jkbz),&
             & kngpsize,1.0d0,dr(1,1),max_iaiqm)
          call dgemm('N','N',max_iaiqm,nval,ibl2-ibl1+1,1.0d0,tmp(1,1,3),max_iaiqm,potr(ibl1,1,jkbz),&
             & kngpsize,0.0d0,di(1,1),max_iaiqm)
          call dgemm('N','N',max_iaiqm,nval,ibl2-ibl1+1,1.0d0,tmp(1,1,4),max_iaiqm,poti(ibl1,1,jkbz),&
             & kngpsize,1.0d0,di(1,1),max_iaiqm)
          if(force_mode)then
             call dgemm('N','N',max_iaiqm,nval*3,ibl2-ibl1+1,1.0d0,tmp2(1,1,1,1),max_iaiqm,potr(ibl1,1,jkbz),&
             & kngpsize,0.0d0,er(1,1,1),max_iaiqm)
             call dgemm('N','N',max_iaiqm,nval*3,ibl2-ibl1+1,1.0d0,tmp2(1,1,1,2),max_iaiqm,poti(ibl1,1,jkbz),&
             & kngpsize,1.0d0,er(1,1,1),max_iaiqm)
             call dgemm('N','N',max_iaiqm,nval*3,ibl2-ibl1+1,1.0d0,tmp2(1,1,1,3),max_iaiqm,potr(ibl1,1,jkbz),&
             & kngpsize,0.0d0,ei(1,1,1),max_iaiqm)
             call dgemm('N','N',max_iaiqm,nval*3,ibl2-ibl1+1,1.0d0,tmp2(1,1,1,4),max_iaiqm,poti(ibl1,1,jkbz),&
             & kngpsize,1.0d0,ei(1,1,1),max_iaiqm)
          endif

!          di = -di
!          if(force_mode) ei = -ei

          iaiqm = 0
          do ia=1,natm
             it = ityp(ia)
             if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
             do iq=nqitg_sp0(it), nqitg_sp(it)
                l3 = iq2l3(iq)
                iil = mod(l3,4)

                do m = 1, 2*l3+1
                   iaiqm = iaiqm + 1
                   ilm3 = l3*l3+m
                   do m_ = 1, nval
                      if(occup_val(m_,kk) < DELTA) cycle
                      do ip=1,nc(m,iq)
                         lmt1 = nc2lmt1(ip,m,iq); lmt2 = nc2lmt2(ip,m,iq)
                         ilmta1 = lmta(lmt1,ia);  ilmta2 = lmta(lmt2,ia)
                         np = nc2n(ip,m,iq)
                         f0 = iwei(ia) * dl2p(lmt1,lmt2,np,it)
                         qvr(ilmta1,m_,jkbz) = qvr(ilmta1,m_,jkbz) + f0 * (dr(iaiqm,m_) * fmr(ilmta2,m_,jkbz) &
                                             & - di(iaiqm,m_) * fmi(ilmta2,m_,jkbz))
                         qvi(ilmta1,m_,jkbz) = qvi(ilmta1,m_,jkbz) + f0 * (di(iaiqm,m_) * fmr(ilmta2,m_,jkbz) &
                                             & + dr(iaiqm,m_) * fmi(ilmta2,m_,jkbz))
                      end do

                      if(force_mode)then
                         do ip=1,nc(m,iq)
                            lmt1 = nc2lmt1(ip,m,iq); lmt2 = nc2lmt2(ip,m,iq)
                            ilmta1 = lmta(lmt1,ia);  ilmta2 = lmta(lmt2,ia)
                            np = nc2n(ip,m,iq)
                            f0 = iwei(ia) * dl2p(lmt1,lmt2,np,it)
                            dqvr(ilmta1,m_,jkbz,1:3) = dqvr(ilmta1,m_,jkbz,1:3) &
                            & + f0 * (dr(iaiqm,m_) * dfmr(ilmta2,1:3) - di(iaiqm,m_) * dfmi(ilmta2,1:3))
                            dqvi(ilmta1,m_,jkbz,1:3) = dqvi(ilmta1,m_,jkbz,1:3) &
                            & + f0 * (di(iaiqm,m_) * dfmr(ilmta2,1:3) + dr(iaiqm,m_) * dfmi(ilmta2,1:3))
                            gqvr(ilmta1,m_,jkbz,1:3) = gqvr(ilmta1,m_,jkbz,1:3) &
                            & + f0 * ( ei(iaiqm,m_,1:3) * fmr(ilmta2,m_,jkbz) + er(iaiqm,m_,1:3) * fmi(ilmta2,m_,jkbz))
                            gqvi(ilmta1,m_,jkbz,1:3) = gqvi(ilmta1,m_,jkbz,1:3) &
                            & + f0 * (-er(iaiqm,m_,1:3) * fmr(ilmta2,m_,jkbz) + ei(iaiqm,m_,1:3) * fmi(ilmta2,m_,jkbz))
                         enddo
                      endif

                   end do
                end do
             end do
          end do
       end do
    enddo

    deallocate(zsr);  deallocate(zsi)
    deallocate(qitg_red); deallocate(ylm_red)
    if(force_mode)then
       deallocate(tmp2)
       deallocate(er)
       deallocate(ei)
       deallocate(gvec)
    endif
    call tstatc0_end(id_sname)
  end subroutine integrate_QijVnm_
#endif
! ========================================================================= 125.Exp

  subroutine integrate_QijVnm_rs2(iqmk,potr,poti,fmrq,fmiq,qvr,qvi)
    integer, intent(in) :: iqmk
    real(kind=DP), intent(in) :: potr(nfftp_exx_nonpara/2)
    real(kind=DP), intent(in) :: poti(nfftp_exx_nonpara/2)
    real(kind=DP), intent(in) :: fmrq(nmesh_rs_aug_max,nlmta)
    real(kind=DP), intent(in) :: fmiq(nmesh_rs_aug_max,nlmta)
    real(kind=DP), intent(out) :: qvr(nlmta)
    real(kind=DP), intent(out) :: qvi(nlmta)
    real(kind=DP), dimension(:), allocatable :: cosqmkr,sinqmkr
    real(kind=DP), allocatable, dimension(:) :: rra,iia
    integer :: ia,it,lmt1,ilmta1
    integer :: imesh,nma,ind
    real(kind=DP) :: co,si
    real(kind=DP) :: qr,qi,qmr,qmi,rr,ii
    integer :: id_sname=-1
    call tstatc0_begin('integrate_QijVnm_rs2 ',id_sname)
    allocate(cosqmkr(nmesh_rs_aug_max));cosqmkr=0.d0
    allocate(sinqmkr(nmesh_rs_aug_max));sinqmkr=0.d0
    allocate(rra(nmesh_rs_aug_max));rra=0.d0
    allocate(iia(nmesh_rs_aug_max));iia=0.d0
    qvr(:) = 0.d0
    qvi(:) = 0.d0
    do ia=1,natm
       it = ityp(ia)
       if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
       call qmk_dot_r(iqmk,ia,cosqmkr,sinqmkr)
       nma = nmesh_rs_aug(ia)
       do imesh=1,nma
          ind = meshxyz_rs_aug(imesh,ia)
          co  =  cosqmkr(imesh)
          si  = -sinqmkr(imesh)
          rra(imesh) =  potr(ind)*co+poti(ind)*si
          iia(imesh) = -poti(ind)*co+potr(ind)*si
       enddo
       do lmt1=1,ilmt(it)
          ilmta1 = lmta(lmt1,ia)
          qr=0.d0;qi=0.d0
          do imesh=1,nma
             qmr = fmrq(imesh,ilmta1)
             qmi = fmiq(imesh,ilmta1)
             rr = rra(imesh)
             ii = iia(imesh)
             qr = qr+qmr*rr-qmi*ii
             qi = qi+qmr*ii+qmi*rr
          enddo
          qvr(ilmta1) = qr
          qvi(ilmta1) = qi
       enddo
    enddo
    deallocate(cosqmkr)
    deallocate(sinqmkr)
    deallocate(rra,iia)
    call tstatc0_end(id_sname)
  end subroutine integrate_QijVnm_rs2

  subroutine integrate_QijVnm_rs(iqmk,potr,poti,fmr,fmi,qvr,qvi,dfmr,dfmi,dqvr,dqvi,gqvr,gqvi)
    integer, intent(in) :: iqmk
    real(kind=DP), intent(in) :: potr(nfftp_exx_nonpara/2)
    real(kind=DP), intent(in) :: poti(nfftp_exx_nonpara/2)
    real(kind=DP), intent(in) :: fmr(nlmta)
    real(kind=DP), intent(in) :: fmi(nlmta)
    real(kind=DP), intent(out) :: qvr(nlmta)
    real(kind=DP), intent(out) :: qvi(nlmta)
    real(kind=DP), intent(in),optional  :: dfmr(nlmta,3)
    real(kind=DP), intent(in),optional  :: dfmi(nlmta,3)
    real(kind=DP), intent(out),optional :: dqvr(nlmta,3)
    real(kind=DP), intent(out),optional :: dqvi(nlmta,3)
    real(kind=DP), intent(out),optional :: gqvr(nlmta,3)
    real(kind=DP), intent(out),optional :: gqvi(nlmta,3)
    integer :: n,it,ia,imesh,lmt1,lmt2,il1,il2,tau1,tau2,ilmta1,ilmta2,lmtp
    real(kind=DP), dimension(:), allocatable :: cosqmkr,sinqmkr
    real(kind=DP) :: co,si,rr,ii,rrr,iii,qr,qi,qqr,qqi,qm,rr0,ii0
    real(kind=DP), allocatable, dimension(:) :: rra,iia
    integer :: nma,ind
    real(kind=DP) :: er(3), ei(3),dqm(3),dqqr(3),dqqi(3),der(3),dei(3),drr(3),dii(3),drr0(3),dii0(3)
    real(kind=DP) :: dqqr0(3),dqqi0(3),der0(3),dei0(3)
    logical :: force_mode = .false.
    integer :: id_sname=-1
    call tstatc0_begin('integrate_QijVnm_rs ',id_sname)
    force_mode = present(dfmr).and.present(dfmi).and. &
               & present(dqvr).and.present(dqvi).and. &
               & present(gqvr).and.present(gqvi)
    allocate(cosqmkr(nmesh_rs_aug_max));cosqmkr=0.d0
    allocate(sinqmkr(nmesh_rs_aug_max));sinqmkr=0.d0
    allocate(rra(nmesh_rs_aug_max));rra=0.d0
    allocate(iia(nmesh_rs_aug_max));iia=0.d0
    qvr(:) = 0.d0
    qvi(:) = 0.d0
    if(force_mode)then
       dqvr(1:nlmta,1:3) = 0.d0
       dqvi(1:nlmta,1:3) = 0.d0
       gqvr(1:nlmta,1:3) = 0.d0
       gqvi(1:nlmta,1:3) = 0.d0
    endif
    do ia=1,natm
       it = ityp(ia)
       if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
       call qmk_dot_r(iqmk,ia,cosqmkr,sinqmkr)
       nma = nmesh_rs_aug(ia)
       do imesh=1,nma
          ind = meshxyz_rs_aug(imesh,ia)
          co  =  cosqmkr(imesh)
          si  = -sinqmkr(imesh)
          rra(imesh) =  potr(ind)*co+poti(ind)*si
          iia(imesh) = -poti(ind)*co+potr(ind)*si
       enddo
       do lmtp=1,nlmtpair(ia)
          lmt1 = plmt1(lmtp,ia)
          lmt2 = plmt2(lmtp,ia)
          ilmta1 = lmta(lmt1,ia)
          ilmta2 = lmta(lmt2,ia)
          rr0 =  fmr(ilmta1)
          ii0 =  fmi(ilmta1)
          qqr=0.d0;qqi=0.d0
          qr=0.d0;qi=0.d0
          rr =  fmr(ilmta2)
          ii =  fmi(ilmta2)
          if(force_mode)then
            dqqr=0.d0;dqqi=0.d0
            dqqr0=0.d0;dqqi0=0.d0
            der=0.d0;dei=0.d0
            der0=0.d0;dei0=0.d0
            drr(1:3) =  dfmr(ilmta2,1:3)
            dii(1:3) =  dfmi(ilmta2,1:3)
            er=0.d0;ei=0.d0
          endif
          do imesh=1,nma
             qm = qr_clm_ylm(imesh,ia,lmtp)
             qr = qr+qm*rra(imesh)
             qi = qi+qm*iia(imesh)
          enddo
          if(force_mode)then
             do imesh=1,nma
                dqm(1:3) = dqr_clm_ylm(imesh,ia,lmtp,1:3)
                er(1:3) = er(1:3)+dqm(1:3)*rra(imesh)
                ei(1:3) = ei(1:3)+dqm(1:3)*iia(imesh)
             enddo
          endif
          qqr = qr*rr-qi*ii
          qqi = qr*ii+qi*rr
          qvr(ilmta1) = qvr(ilmta1)+qqr
          qvi(ilmta1) = qvi(ilmta1)+qqi
          if(force_mode)then
             dqqr(1:3)  = qr*drr(1:3)-qi*dii(1:3)
             dqqi(1:3)  = qr*dii(1:3)+qi*drr(1:3)
             der(1:3)   = er(1:3)*rr-ei(1:3)*ii
             dei(1:3)   = er(1:3)*ii+ei(1:3)*rr
             dqvr(ilmta1,1:3) = dqvr(ilmta1,1:3)+dqqr(1:3)
             dqvi(ilmta1,1:3) = dqvi(ilmta1,1:3)+dqqi(1:3)
#ifdef USE_DERIVATIVE_QRSPS_RSPACE
             gqvr(ilmta1,1:3) = gqvr(ilmta1,1:3)+der(1:3)
             gqvi(ilmta1,1:3) = gqvi(ilmta1,1:3)+dei(1:3)
#endif
          endif
          if(lmt1.ne.lmt2)then
             qvr(ilmta2) = qvr(ilmta2)+(qr*rr0-qi*ii0)
             qvi(ilmta2) = qvi(ilmta2)+(qr*ii0+qi*rr0)
          endif
          if(force_mode.and.lmt1.ne.lmt2)then
             drr0(1:3) =  dfmr(ilmta1,1:3)
             dii0(1:3) =  dfmi(ilmta1,1:3)
             dqqr0(1:3) = qr*drr0(1:3)-qi*dii0(1:3)
             dqqi0(1:3) = qr*dii0(1:3)+qi*drr0(1:3)
             der0(1:3)   = er(1:3)*rr0-ei(1:3)*ii0
             dei0(1:3)   = er(1:3)*ii0+ei(1:3)*rr0
             dqvr(ilmta2,1:3) = dqvr(ilmta2,1:3)+dqqr0(1:3)
             dqvi(ilmta2,1:3) = dqvi(ilmta2,1:3)+dqqi0(1:3)
#ifdef USE_DERIVATIVE_QRSPS_RSPACE
             gqvr(ilmta2,1:3) = gqvr(ilmta2,1:3)+der0(1:3)
             gqvi(ilmta2,1:3) = gqvi(ilmta2,1:3)+dei0(1:3)
#endif
          endif
       enddo
    enddo
    deallocate(cosqmkr)
    deallocate(sinqmkr)
    deallocate(rra,iia)
    call tstatc0_end(id_sname)
  end subroutine integrate_QijVnm_rs

  subroutine integrate_QijVnm_rs3(iqmk,potr,poti,fmr,fmi,qvr,qvi,dfmr,dfmi,dqvr,dqvi,gqvr,gqvi)
    integer, intent(in) :: iqmk
    real(kind=DP), intent(in) :: potr(nfftp_exx_nonpara/2,ista_nval:iend_nval)
    real(kind=DP), intent(in) :: poti(nfftp_exx_nonpara/2,ista_nval:iend_nval)
    real(kind=DP), intent(in) :: fmr(ista_nval:iend_nval,nlmta)
    real(kind=DP), intent(in) :: fmi(ista_nval:iend_nval,nlmta)
    real(kind=DP), intent(out) :: qvr(nlmta,ista_nval:iend_nval)
    real(kind=DP), intent(out) :: qvi(nlmta,ista_nval:iend_nval)
    real(kind=DP), intent(in),optional  :: dfmr(nlmta,ista_nval:iend_nval,3)
    real(kind=DP), intent(in),optional  :: dfmi(nlmta,ista_nval:iend_nval,3)
    real(kind=DP), intent(out),optional :: dqvr(nlmta,ista_nval:iend_nval,3)
    real(kind=DP), intent(out),optional :: dqvi(nlmta,ista_nval:iend_nval,3)
    real(kind=DP), intent(out),optional :: gqvr(nlmta,ista_nval:iend_nval,3)
    real(kind=DP), intent(out),optional :: gqvi(nlmta,ista_nval:iend_nval,3)
    integer :: n,it,ia,imesh,lmt1,lmt2,il1,il2,tau1,tau2,ilmta1,ilmta2,lmtp,m,i
    real(kind=DP), dimension(:), allocatable :: cosqmkr,sinqmkr
    real(kind=DP) :: co,si,rr,ii,rrr,iii,qr,qi,qqr,qqi,qm,rr0,ii0
    real(kind=DP), allocatable, dimension(:,:) :: dres1,dres2,dres3,dres4
    real(kind=DP), allocatable, dimension(:,:,:) :: dres5,dres6,dres7,dres8
    real(kind=DP), allocatable, dimension(:,:) :: tmp1,tmp2
    real(kind=DP), allocatable, dimension(:,:,:) :: tmp3,tmp4
    real(kind=DP), allocatable, dimension(:,:) :: potrt,potit
    integer :: nma,ind
    real(kind=DP) :: er(3), ei(3),dqm(3),dqqr(3),dqqi(3),der(3),dei(3),drr(3),dii(3),drr0(3),dii0(3)
    real(kind=DP) :: dqqr0(3),dqqi0(3),der0(3),dei0(3)
    logical :: force_mode = .false.
    integer :: id_sname=-1
    call tstatc0_begin('integrate_QijVnm_rs ',id_sname)
    force_mode = present(dfmr).and.present(dfmi).and. &
               & present(dqvr).and.present(dqvi).and. &
               & present(gqvr).and.present(gqvi)
    allocate(cosqmkr(nmesh_rs_aug_max));cosqmkr=0.d0
    allocate(sinqmkr(nmesh_rs_aug_max));sinqmkr=0.d0
    allocate(tmp1(nmesh_rs_aug_max,maxval(nlmtpair)))
    allocate(tmp2(nmesh_rs_aug_max,maxval(nlmtpair)))
    allocate(dres1(maxval(nlmtpair),np_nval))
    allocate(dres2(maxval(nlmtpair),np_nval))
    allocate(dres3(maxval(nlmtpair),np_nval))
    allocate(dres4(maxval(nlmtpair),np_nval))
    if(force_mode)then
       allocate(tmp3(nmesh_rs_aug_max,maxval(nlmtpair),3))
       allocate(tmp4(nmesh_rs_aug_max,maxval(nlmtpair),3))
       allocate(dres5(maxval(nlmtpair),np_nval,3))
       allocate(dres6(maxval(nlmtpair),np_nval,3))
       allocate(dres7(maxval(nlmtpair),np_nval,3))
       allocate(dres8(maxval(nlmtpair),np_nval,3))
    endif
    allocate(potrt(nmesh_rs_aug_max,np_nval))
    allocate(potit(nmesh_rs_aug_max,np_nval))
    qvr(:,:) = 0.d0
    qvi(:,:) = 0.d0
    if(force_mode)then
       dqvr(1:nlmta,ista_nval:iend_nval,1:3) = 0.d0
       dqvi(1:nlmta,ista_nval:iend_nval,1:3) = 0.d0
       gqvr(1:nlmta,ista_nval:iend_nval,1:3) = 0.d0
       gqvi(1:nlmta,ista_nval:iend_nval,1:3) = 0.d0
    endif
    do ia=1,natm
       it = ityp(ia)
       if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
       call qmk_dot_r(iqmk,ia,cosqmkr,sinqmkr)
       nma = nmesh_rs_aug(ia)
       potrt=0.d0;potit=0.d0
       do m=ista_nval,iend_nval
          do imesh=1,nma
             ind = meshxyz_rs_aug(imesh,ia)
             potrt(imesh,m-ista_nval+1) = potr(ind,m)
             potit(imesh,m-ista_nval+1) = poti(ind,m)
          enddo
       enddo
       tmp1=0.d0;tmp2=0.d0
       do lmtp=1,nlmtpair(ia)
          do imesh=1,nma
             co  =  cosqmkr(imesh)
             si  = -sinqmkr(imesh)
             qm = qr_clm_ylm(imesh,ia,lmtp)
             tmp1(imesh,lmtp) = qm*co
             tmp2(imesh,lmtp) = qm*si
             if(force_mode)then
                dqm(1:3) = dqr_clm_ylm(imesh,ia,lmtp,1:3)
                tmp3(imesh,lmtp,:) = co * dqm(:)
                tmp4(imesh,lmtp,:) = si * dqm(:)
!!                dqm(1:3) = dqr_clm_ylm(imesh,ia,lmtp,1:3)
             endif
          enddo
       enddo

       call dgemm('T','N',nlmtpair(ia),np_nval,nma,1.d0,&
                 & tmp1,nmesh_rs_aug_max,potrt,nmesh_rs_aug_max,0.0d0,dres1,maxval(nlmtpair))
       call dgemm('T','N',nlmtpair(ia),np_nval,nma,1.d0,&
                 & tmp2,nmesh_rs_aug_max,potit,nmesh_rs_aug_max,0.0d0,dres2,maxval(nlmtpair))
       call dgemm('T','N',nlmtpair(ia),np_nval,nma,1.d0,&
                 & tmp1,nmesh_rs_aug_max,potit,nmesh_rs_aug_max,0.0d0,dres3,maxval(nlmtpair))
       call dgemm('T','N',nlmtpair(ia),np_nval,nma,1.d0,&
                 & tmp2,nmesh_rs_aug_max,potrt,nmesh_rs_aug_max,0.0d0,dres4,maxval(nlmtpair))
       if (force_mode) then
          do i=1,3
             call dgemm('T','N',nlmtpair(ia),np_nval,nma,1.d0,&
                    & tmp3(1,1,i),nmesh_rs_aug_max,potrt,nmesh_rs_aug_max,0.0d0,dres5(1,1,i),maxval(nlmtpair))
             call dgemm('T','N',nlmtpair(ia),np_nval,nma,1.d0,&
                    & tmp4(1,1,i),nmesh_rs_aug_max,potit,nmesh_rs_aug_max,0.0d0,dres6(1,1,i),maxval(nlmtpair))
             call dgemm('T','N',nlmtpair(ia),np_nval,nma,1.d0,&
                    & tmp3(1,1,i),nmesh_rs_aug_max,potit,nmesh_rs_aug_max,0.0d0,dres7(1,1,i),maxval(nlmtpair))
             call dgemm('T','N',nlmtpair(ia),np_nval,nma,1.d0,&
                    & tmp4(1,1,i),nmesh_rs_aug_max,potrt,nmesh_rs_aug_max,0.0d0,dres8(1,1,i),maxval(nlmtpair))
          enddo
       endif

       do m=ista_nval,iend_nval
       do lmtp=1,nlmtpair(ia)
          lmt1 = plmt1(lmtp,ia)
          lmt2 = plmt2(lmtp,ia)
          ilmta1 = lmta(lmt1,ia)
          ilmta2 = lmta(lmt2,ia)
          rr0 =  fmr(m,ilmta1)
          ii0 =  fmi(m,ilmta1)
          qqr=0.d0;qqi=0.d0
          qr=0.d0;qi=0.d0
          rr =  fmr(m,ilmta2)
          ii =  fmi(m,ilmta2)
          if(force_mode)then
            dqqr=0.d0;dqqi=0.d0
            dqqr0=0.d0;dqqi0=0.d0
            der=0.d0;dei=0.d0
            der0=0.d0;dei0=0.d0
            drr(1:3) =  dfmr(ilmta2,m,1:3)
            dii(1:3) =  dfmi(ilmta2,m,1:3)
            er=0.d0;ei=0.d0
          endif

          qr =  dres1(lmtp,m-ista_nval+1) + dres2(lmtp,m-ista_nval+1)
          qi = -dres3(lmtp,m-ista_nval+1) + dres4(lmtp,m-ista_nval+1)
          qqr = qr*rr-qi*ii
          qqi = qr*ii+qi*rr
          qvr(ilmta1,m) = qvr(ilmta1,m)+qqr
          qvi(ilmta1,m) = qvi(ilmta1,m)+qqi
          if(force_mode)then
             er(:) =  dres5(lmtp,m-ista_nval+1,:) + dres6(lmtp,m-ista_nval+1,:)
             ei(:) = -dres7(lmtp,m-ista_nval+1,:) + dres8(lmtp,m-ista_nval+1,:)
             dqqr(1:3)  = qr*drr(1:3)-qi*dii(1:3)
             dqqi(1:3)  = qr*dii(1:3)+qi*drr(1:3)
             der(1:3)   = er(1:3)*rr-ei(1:3)*ii
             dei(1:3)   = er(1:3)*ii+ei(1:3)*rr
             dqvr(ilmta1,m,1:3) = dqvr(ilmta1,m,1:3)+dqqr(1:3)
             dqvi(ilmta1,m,1:3) = dqvi(ilmta1,m,1:3)+dqqi(1:3)
#ifdef USE_DERIVATIVE_QRSPS_RSPACE
             gqvr(ilmta1,m,1:3) = gqvr(ilmta1,m,1:3)+der(1:3)
             gqvi(ilmta1,m,1:3) = gqvi(ilmta1,m,1:3)+dei(1:3)
#endif
          endif
          if(lmt1.ne.lmt2)then
             qvr(ilmta2,m) = qvr(ilmta2,m)+(qr*rr0-qi*ii0)
             qvi(ilmta2,m) = qvi(ilmta2,m)+(qr*ii0+qi*rr0)
          endif
          if(force_mode.and.lmt1.ne.lmt2)then
             drr0(1:3) =  dfmr(ilmta1,m,1:3)
             dii0(1:3) =  dfmi(ilmta1,m,1:3)
             dqqr0(1:3) = qr*drr0(1:3)-qi*dii0(1:3)
             dqqi0(1:3) = qr*dii0(1:3)+qi*drr0(1:3)
             der0(1:3)   = er(1:3)*rr0-ei(1:3)*ii0
             dei0(1:3)   = er(1:3)*ii0+ei(1:3)*rr0
             dqvr(ilmta2,m,1:3) = dqvr(ilmta2,m,1:3)+dqqr0(1:3)
             dqvi(ilmta2,m,1:3) = dqvi(ilmta2,m,1:3)+dqqi0(1:3)
#ifdef USE_DERIVATIVE_QRSPS_RSPACE
             gqvr(ilmta2,m,1:3) = gqvr(ilmta2,m,1:3)+der0(1:3)
             gqvi(ilmta2,m,1:3) = gqvi(ilmta2,m,1:3)+dei0(1:3)
#endif
          endif
       enddo
       enddo
    enddo
    deallocate(cosqmkr)
    deallocate(sinqmkr)
    deallocate(tmp1)
    deallocate(tmp2)
    deallocate(dres1)
    deallocate(dres2)
    deallocate(dres3)
    deallocate(dres4)
    if(force_mode)then
      deallocate(tmp3)
      deallocate(tmp4)
      deallocate(dres5)
      deallocate(dres6)
      deallocate(dres7)
      deallocate(dres8)
    endif
    deallocate(potrt)
    deallocate(potit)
    call tstatc0_end(id_sname)
  end subroutine integrate_QijVnm_rs3

  subroutine add_Vx_hard_part(ik,vxw,sumqvr,sumqvi)
  !!subroutine add_Vx_hard_part(ik,vxw1,sumqvr,sumqvi)
    implicit none
    integer, intent(in) :: ik
    real(kind=DP), intent(inout) :: vxw(maxval(np_g1k),kimg)
    !!real(kind=DP), intent(inout) :: vxw1(kg1,kimg)
    real(kind=DP), intent(in) :: sumqvr(nlmta)
    real(kind=DP), intent(in) :: sumqvi(nlmta)

    integer :: it,ia,i,il1,ilmtt1,ilmta1,iil,iksnl,ikk
    integer :: lmt1,lmtt1
    real(kind=DP) :: fac, ph
    real(kind=DP), allocatable :: cosgt(:) ! d(kg1)
    real(kind=DP), allocatable :: singt(:) ! d(kg1)
    integer,save  :: id_sname = -1
    integer :: iadd
    call tstatc0_begin('add_Vx_hard_part ',id_sname)

    allocate(cosgt(kg1))
    allocate(singt(kg1))
    ikk = k_index(ik)
    iksnl = (ik-1)/nspin + 1

    do ia=1,natm
       it = ityp(ia)
       if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
       do i=1,iba(ik)
          ph = PAI2 * dot_product(pos(ia,1:3),vkxyz_for_outerloop(ikk,1:3,BUCS)+ngabc(nbase(i,ik),1:3))
          cosgt(i) = cos(ph)
          singt(i) = sin(ph)
       end do
       do lmt1 = 1, ilmt(it)
          il1  = ltp(lmt1,it)
          lmtt1 = lmtt(lmt1,it)
          ilmta1 = lmta(lmt1,ia)
          iil=mod(il1-1,4)
          if(iil==0) then
             do i = ista_g1k(ik), iend_g1k(ik)
                iadd = i - ista_g1k(ik) + 1
                fac = iwei(ia)*snl(iadd,lmtt1,iksnl)
                vxw(iadd,1) = vxw(iadd,1) + fac * (cosgt(i)*sumqvr(ilmta1)+singt(i)*sumqvi(ilmta1))
                vxw(iadd,2) = vxw(iadd,2) + fac * (cosgt(i)*sumqvi(ilmta1)-singt(i)*sumqvr(ilmta1))
             end do
          else if(iil==1) then
             do i = ista_g1k(ik), iend_g1k(ik)
                iadd = i - ista_g1k(ik) + 1
                fac = iwei(ia)*snl(iadd,lmtt1,iksnl)
                vxw(iadd,1) = vxw(iadd,1) + fac * (cosgt(i)*sumqvi(ilmta1)-singt(i)*sumqvr(ilmta1))
                vxw(iadd,2) = vxw(iadd,2) - fac * (cosgt(i)*sumqvr(ilmta1)+singt(i)*sumqvi(ilmta1))
             end do
          else if(iil==2) then
             do i = ista_g1k(ik), iend_g1k(ik)
                iadd = i - ista_g1k(ik) + 1
                fac = -iwei(ia)*snl(iadd,lmtt1,iksnl)
                vxw(iadd,1) = vxw(iadd,1) + fac * (cosgt(i)*sumqvr(ilmta1)+singt(i)*sumqvi(ilmta1))
                vxw(iadd,2) = vxw(iadd,2) + fac * (cosgt(i)*sumqvi(ilmta1)-singt(i)*sumqvr(ilmta1))
             end do
          else if(iil==3) then
             do i = ista_g1k(ik), iend_g1k(ik)
                iadd = i - ista_g1k(ik) + 1
                fac = -iwei(ia)*snl(iadd,lmtt1,iksnl)
                vxw(iadd,1) = vxw(iadd,1) + fac * (cosgt(i)*sumqvi(ilmta1)-singt(i)*sumqvr(ilmta1))
                vxw(iadd,2) = vxw(iadd,2) - fac * (cosgt(i)*sumqvr(ilmta1)+singt(i)*sumqvi(ilmta1))
             end do
          end if
       end do
    end do

    deallocate(cosgt)
    deallocate(singt)

    call tstatc0_end(id_sname)
  end subroutine add_Vx_hard_part

  subroutine get_Rot_betar_dot_WFs(ik,ib,jop,kop,jtrs,fsr,fsi,dfsr,dfsi)
    implicit none
    integer, intent(in)                          :: ik,ib,jop,kop,jtrs
    real(kind=DP), intent(out), dimension(nlmta) :: fsr,fsi
    real(kind=DP), intent(out), optional, dimension(nlmta,3) :: dfsr,dfsi

    integer :: ilmta1,ilmta2,isph1,isph2,iy
    integer :: ia,i
    real(kind=DP) :: prim(3)
    real(kind=DP) :: coskt(natm), sinkt(natm)
    real(kind=DP) :: ph, tmpr, tmpi
    logical :: force_mode=.false.
    integer,save  :: id_sname = -1

! === ASMS === 2019/10/28
    integer :: ixyz1, ixyz2
    real(kind=DP) :: c1, cvec1(3), cvec2(3)
! === ASMS === 2019/10/28

    call tstatc0_begin('get_Rot_betar_dot_WFs ',id_sname)
    force_mode = present(dfsr).and.present(dfsi)
    if(jop==1) then
       fsr(1:nlmta) = fsr_exx(ib,1:nlmta,ik)
       fsi(1:nlmta) = fsi_exx(ib,1:nlmta,ik)
       if(force_mode)then
          dfsr(1:nlmta,1:3) = dfsr_exx(ib,1:nlmta,ik,1:3)
          dfsi(1:nlmta,1:3) = dfsi_exx(ib,1:nlmta,ik,1:3)
       endif
    else
       fsr = 0.d0
       fsi = 0.d0
       if(force_mode)then
          dfsr = 0.d0
          dfsi = 0.d0
       endif
       do ilmta1=1,nlmta
          isph1=isph_lmta(ilmta1)
          do iy=1,nylm(isph1,jop)
             isph2 = iylm(iy,isph1,jop)
             ilmta2=lmta_rot(iy,ilmta1,jop)
             fsr(ilmta1) = fsr(ilmta1) +crotylm(iy,isph1,jop) *fsr_exx(ib,ilmta2,ik)
             fsi(ilmta1) = fsi(ilmta1) +crotylm(iy,isph1,jop) *fsi_exx(ib,ilmta2,ik)
             if(force_mode)then
                dfsr(ilmta1,1:3) = dfsr(ilmta1,1:3) &
                     &            +crotylm(iy,isph1,jop) *dfsr_exx(ib,ilmta2,ik,1:3)
                dfsi(ilmta1,1:3) = dfsi(ilmta1,1:3) &
                     &            +crotylm(iy,isph1,jop) *dfsi_exx(ib,ilmta2,ik,1:3)
             endif
          end do
       end do
    end if

    if(jop>1) then
       do ia=1,natm
          prim(1:3) = matmul(opdir(:,:,kop),pos(ia,1:3)) + tau(1:3,kop,BUCS) - pos(napt(ia,kop),1:3)
          ph = PAI2 * dot_product(vkxyz_for_innerloop(ik,1:3,BUCS),prim)
          coskt(ia) = cos(ph)
          sinkt(ia) = sin(ph)
       end do
       do ilmta1=1,nlmta
          ia = ia_lmta(ilmta1)
          tmpr = fsr(ilmta1)
          tmpi = fsi(ilmta1)
          fsr(ilmta1) = coskt(ia) * tmpr - sinkt(ia) * tmpi
          fsi(ilmta1) = coskt(ia) * tmpi + sinkt(ia) * tmpr
          if(force_mode)then
             do i=1,3
                tmpr = dfsr(ilmta1,i)
                tmpi = dfsi(ilmta1,i)
                dfsr(ilmta1,i) = coskt(ia) * tmpr - sinkt(ia) * tmpi
                dfsi(ilmta1,i) = coskt(ia) * tmpi + sinkt(ia) * tmpr
             end do
          endif
       end do
    end if

    if(jop>1 .and. force_mode )then
       do ilmta1=1,nlmta
          cvec1(:) = dfsr(ilmta1,:)
          cvec2(:) = dfsi(ilmta1,:)
          Do ixyz1=1, 3
             tmpr = 0.d0
             tmpi = 0.d0
             Do ixyz2=1, 3
#ifdef DERIVATIVE_IN_CARTS
                c1 = op(ixyz2,ixyz1,kop)
#else
                c1 = opdir(ixyz2,ixyz1,kop)
#endif
                tmpr = tmpr +c1 *cvec1(ixyz2)
                tmpi = tmpi +c1 *cvec2(ixyz2)
             End Do
             dfsr(ilmta1,ixyz1) = tmpr
             dfsi(ilmta1,ixyz1) = tmpi
          End Do
       End Do
    endif

    if(jtrs==1) fsi = -fsi
    if(force_mode.and.jtrs==1) dfsi = -dfsi
    call tstatc0_end(id_sname)
  end subroutine get_Rot_betar_dot_WFs

  subroutine get_expkt_fs_b(ik,ib,fsr_l,fsi_l,efsr_l,efsi_l)
    implicit none
    integer, intent(in) :: ik,ib
    real(kind=DP), intent(in), dimension(np_fs) :: fsr_l,fsi_l
    real(kind=DP), intent(out), dimension(nlmta) :: efsr_l,efsi_l

    integer :: it, ia, lmt, ilmta
    integer :: iadd
    real(kind=DP) :: ph, coskt, sinkt
    efsr_l(1:nlmta) = 0.d0
    efsi_l(1:nlmta) = 0.d0
    if(k_symmetry(ik)==GAMMA)then
       do ia=1,natm
          it = ityp(ia)
          if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
          do lmt = 1, ilmt(it)
             ilmta = lmta(lmt,ia)
             if(ilmta >= ista_fs .and. ilmta <= iend_fs) then
                iadd = ilmta - ista_fs + 1
                efsr_l(ilmta) = fsr_l(iadd)
             end if
          end do
       end do
    else
       do ia=1,natm
          it = ityp(ia)
          if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
          ph = PAI2 * dot_product(pos(ia,1:3),vkxyz(ik,1:3,BUCS))
          coskt = cos(ph)
          sinkt = sin(ph)
          do lmt = 1, ilmt(it)
             ilmta = lmta(lmt,ia)
             if(ilmta >= ista_fs .and. ilmta <= iend_fs) then
                iadd = ilmta - ista_fs + 1
                efsr_l(ilmta) = coskt * fsr_l(iadd) - sinkt * fsi_l(iadd)
                efsi_l(ilmta) = sinkt * fsr_l(iadd) + coskt * fsi_l(iadd)
             end if
          end do
       end do
    endif
    call mpi_allreduce(MPI_IN_PLACE,efsr_l,nlmta,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
    call mpi_allreduce(MPI_IN_PLACE,efsi_l,nlmta,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)

  end subroutine get_expkt_fs_b

  subroutine get_expkt_fs(ik,fsr_l,fsi_l,efsr_l,efsi_l)
    implicit none
    integer, intent(in) :: ik
    real(kind=DP), intent(in), dimension(np_e,np_fs,ista_k:iend_k) :: fsr_l,fsi_l
    real(kind=DP), intent(out), dimension(np_e,nlmta) :: efsr_l,efsi_l

    integer :: it, ia, lmt, ilmta
    integer :: iadd
    real(kind=DP) :: ph, coskt, sinkt
    integer :: ib
    efsr_l = 0.d0
    efsi_l = 0.d0
    do ib=1,np_e
       call get_expkt_fs_b(ik,ib,fsr_l(ib,1:np_fs,ik),fsi_l(ib,1:np_fs,ik),efsr_l(ib,1:nlmta),efsi_l(ib,1:nlmta))
    enddo

  end subroutine get_expkt_fs

  subroutine m_ES_EXX_Force(force)
    implicit none
    real(kind=DP), intent(out) :: force(natm,3)

    integer :: ik,ib,ig,ispin,ia,iadd
    real(kind=DP), allocatable :: force_l(:,:) ! d(natm,3)
    real(kind=DP), allocatable :: force_mpi(:,:) ! d(natm,3)
    real(kind=DP), allocatable :: efsr_l(:,:) ! d(np_e,nlmta)
    real(kind=DP), allocatable :: efsi_l(:,:) ! d(np_e,nlmta)
    real(kind=DP), allocatable :: defsr_l(:,:,:) ! d(np_e,nlmta,3)
    real(kind=DP), allocatable :: defsi_l(:,:,:) ! d(np_e,nlmta,3)
    real(kind=DP), allocatable, dimension(:) :: zajbuf_r,zajbuf_i
    integer :: kgw,kgv
    integer,save  :: id_sname = -1

    real(kind=DP), allocatable :: force_kdep(:,:,:), force_kdep_mpi(:,:,:)

    if(modnrm /= EXECUT) then
       force = 0.d0
       return
    end if

    call tstatc0_begin('m_ES_EXX_Force ',id_sname,level=1)

    call m_FFT_alloc_WF_work()

    allocate(force_l(natm,3))
    allocate(efsr_l(np_e,nlmta));efsr_l=0.d0
    allocate(efsi_l(np_e,nlmta));efsi_l=0.d0

    !!write(nfout,'("m_ES_EXX_force")')

    allocate(dfsr_l(np_e,nlmta,ista_k:iend_k,3));dfsr_l=0.d0
    allocate(dfsi_l(np_e,nlmta,ista_k:iend_k,3));dfsi_l=0.d0
    call drv_betar_dot_WFs_exx() !-> dfsr_l, dfsi_l

    allocate(dfsr_exx(nval,nlmta,kv3,3));dfsr_exx=0.d0
    allocate(dfsi_exx(nval,nlmta,kv3,3));dfsi_exx=0.d0
    call gather_drv_bdw_exx() !-> dfsr_exx, dfsi_exx

    if ( ipriexx >= 2 ) then
       allocate( force_kdep(natm,3,kv3) ); force_kdep = 0.0d0
    endif

    force = 0.d0
    do ispin=1,nspin,af+1
       do ik=ispin, kv3+ispin-nspin, nspin
          if(map_k(ik) /= myrank_k) cycle ! MPI
          call get_expkt_fs(ik,fsr_l,fsi_l,efsr_l,efsi_l)
          do ib=1,np_e   ! MPI
             if(occup_l(ib,ik) < DELTA) cycle
             force_l = 0.d0
             if(sw_change_axis==ON)then
                allocate(zajbuf_r(kg1));zajbuf_r=0.d0
                allocate(zajbuf_i(kg1));zajbuf_i=0.0d0
                do ig=ista_g1k(ik),iend_g1k(ik)
                   iadd = ig-ista_g1k(ik)+1
                   zajbuf_r(ig) = zaj_l(iadd,ib,ik,1)
                   zajbuf_i(ig) = zaj_l(iadd,ib,ik,kimg)
                enddo
                call mpi_allreduce(MPI_IN_PLACE,zajbuf_r,kg1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
                call mpi_allreduce(MPI_IN_PLACE,zajbuf_i,kg1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
                kgw = kg1;kgv=maxval(np_g1k)
             else
                allocate(zajbuf_r(maxval(np_g1k)));zajbuf_r(1:np_g1k(ik))=zaj_l(1:np_g1k(ik),ib,ik,1)
                allocate(zajbuf_i(maxval(np_g1k)));zajbuf_i(1:np_g1k(ik))=zaj_l(1:np_g1k(ik),ib,ik,kimg)
                kgw = maxval(np_g1k);kgv=maxval(np_g1k)
             endif

             if(sw_change_axis==ON)then
               if(sw_rspace_hyb == ON .and. sw_rspace_hyb_dgm == ON.and.modnrm==EXECUT) then
                 if(kimg==1) then
                    call apply_Vx_to_WF_2D_rs_dgm(ispin,ib,ik,kgw,kgv,zajbuf_r,zajbuf_r &
                        & ,efsr_l(ib,1:nlmta),efsr_l(ib,1:nlmta) &
                        & ,dbdwr=dfsr_l(ib,1:nlmta,ik,1:3),dbdwi=dfsr_l(ib,1:nlmta,ik,1:3),force_l=force_l)
                 else
                    call apply_Vx_to_WF_2D_rs_dgm(ispin,ib,ik,kgw,kgv,zajbuf_r,zajbuf_i &
                        & ,efsr_l(ib,1:nlmta),efsi_l(ib,1:nlmta) &
                        & ,dbdwr=dfsr_l(ib,1:nlmta,ik,1:3),dbdwi=dfsi_l(ib,1:nlmta,ik,1:3),force_l=force_l)
                 end if
               else
                 if(kimg==1) then
                    call apply_Vx_to_WF_2D(ispin,ib,ik,kgw,kgv,zajbuf_r,zajbuf_r &
                        & ,efsr_l(ib,1:nlmta),efsr_l(ib,1:nlmta) &
                        & ,dbdwr=dfsr_l(ib,1:nlmta,ik,1:3),dbdwi=dfsr_l(ib,1:nlmta,ik,1:3),force_l=force_l)
                 else
                    call apply_Vx_to_WF_2D(ispin,ib,ik,kgw,kgv,zajbuf_r,zajbuf_i &
                        & ,efsr_l(ib,1:nlmta),efsi_l(ib,1:nlmta) &
                        & ,dbdwr=dfsr_l(ib,1:nlmta,ik,1:3),dbdwi=dfsi_l(ib,1:nlmta,ik,1:3),force_l=force_l)
                 end if
               endif
             else
               if(kimg==1) then
                  call apply_Vx_to_WF(ispin,ib,ik,kgw,kgv,zajbuf_r,zajbuf_r &
                      & ,efsr_l(ib,1:nlmta),efsr_l(ib,1:nlmta) &
                      & ,dbdwr=dfsr_l(ib,1:nlmta,ik,1:3),dbdwi=dfsr_l(ib,1:nlmta,ik,1:3),force_l=force_l)
               else
                  call apply_Vx_to_WF(ispin,ib,ik,kgw,kgv,zajbuf_r,zajbuf_i &
                      & ,efsr_l(ib,1:nlmta),efsi_l(ib,1:nlmta) &
                      & ,dbdwr=dfsr_l(ib,1:nlmta,ik,1:3),dbdwi=dfsi_l(ib,1:nlmta,ik,1:3),force_l=force_l)
               end if
             endif
             force = force + occup_l(ib,ik) * force_l
             if ( ipriexx >= 2 ) then
                force_kdep(:,:,ik) = force_kdep(:,:,ik) +occup_l(ib,ik) * force_l
             endif

             deallocate(zajbuf_r)
             deallocate(zajbuf_i)
          end do
       end do
    end do

    force = -force / dble(kv3/nspin)
    if(nspin==1) force = 2.d0 * force

    deallocate(force_l)
    deallocate(efsr_l)
    deallocate(efsi_l)
    deallocate(dfsr_l)
    deallocate(dfsi_l)
    deallocate(dfsr_exx)
    deallocate(dfsi_exx)

    call m_FFT_dealloc_WF_work()

    if(npes>1) then
       allocate(force_mpi(natm,3))
!       call mpi_allreduce(force,force_mpi,natm*3,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_wor!d,ierr)
!       force = force_mpi
!       call mpi_allreduce(force,force_mpi,natm*3,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ge_world,ierr)
       call mpi_allreduce(force,force_mpi,natm*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_CommGroup,ierr)
       force_mpi = force_mpi /dble(nrank_g)
       force = force_mpi
       deallocate(force_mpi)
    end if

    !! Reduced coordinate -> Cartesian coordinate

! ==== ASMS ====
!    do ia=1,natm
!       force(ia,1:3) = matmul(transpose(rltv),force(ia,1:3))
!    end do
!
#ifndef DERIVATIVE_IN_CARTS
    do ia=1,natm
       call pucv2cart( rltv, force(ia,1), force(ia,2), force(ia,3) )
    End do
#endif
! ==== ASMS ====

    if ( ipriexx >= 2 ) then
       force_kdep = -force_kdep / dble(kv3/nspin)
       if(nspin==1) force_kdep = 2.d0 * force_kdep
       if(npes>1) then
          allocate(force_kdep_mpi(natm,3,kv3))
          call mpi_allreduce( force_kdep,force_kdep_mpi,natm*3*kv3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_CommGroup,ierr)
          force_kdep_mpi = force_kdep_mpi /dble(nrank_g)
          force_kdep = force_kdep_mpi
          deallocate(force_kdep_mpi)
       end if
#ifndef DERIVATIVE_IN_CARTS
       Do ik=1, kv3
          do ia=1,natm
             call pucv2cart( rltv, force_kdep(ia,1,ik), &
                  &          force_kdep(ia,2,ik), force_kdep(ia,3,ik) )
          End do
       End Do
#endif
       if ( mype == 0 ) then
          write(nfout,*)
#ifdef DERIVATIVE_IN_CARTS
          write(nfout,'(A)') "#DERIVATIVE_IN_CARTS is defined"
#endif
#ifdef USE_DERIVATIVE_QRSPS_RSPACE
          write(nfout,'(A)') "#USE_DERIVATIVE_QRSPS_RSPACE is defined"
#endif
          Do ik=1, kv3
             write(nfout,'(A,I5,3F20.10)') "# ik = ", ik, vkxyz(ik,1:3,BUCS)
             write(nfout,'(A)')            "# ia,   fx,      fy,      fz"
             Do ia=1, natm
                write(nfout,'(I5,3F20.10)') ia, force_kdep(ia,1:3,ik)
             End Do
          End Do
       endif
       deallocate( force_kdep )
    endif

    call tstatc0_end(id_sname)
  end subroutine m_ES_EXX_Force

  subroutine sum_EXX_force_terms(force_l,bdwr,bdwi,dbdwr,dbdwi,sumqvr,sumqvi,sumdqvr,sumdqvi)
    implicit none
    real(kind=DP), intent(inout) :: force_l(natm,3)
    real(kind=DP), intent(in) :: bdwr(nlmta), bdwi(nlmta)
    real(kind=DP), intent(in) :: dbdwr(nlmta,3), dbdwi(nlmta,3)
    real(kind=DP), intent(in) :: sumqvr(nlmta), sumqvi(nlmta)
    real(kind=DP), intent(in) :: sumdqvr(nlmta,3), sumdqvi(nlmta,3)

    integer :: ia,it,lmt1,ip

    do ia=1,natm
       it=ityp(ia)
       if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
       do lmt1 = 1, ilmt(it)
          ip = lmta(lmt1,ia)
          force_l(ia,1:3) = force_l(ia,1:3) &
                        & + dbdwr(ip,1:3)*sumqvr(ip) + dbdwi(ip,1:3)*sumqvi(ip) &
                        & + bdwr(ip)*sumdqvr(ip,1:3) + bdwi(ip)*sumdqvi(ip,1:3)
       end do
    end do
  end subroutine sum_EXX_force_terms

  subroutine gather_drv_bdw_exx()
    implicit none

    integer :: ik,ib,ib1,ibm,irev
    real(kind=DP), allocatable :: dfsr_mpi(:,:,:,:) ! d(nval,nlmta,kv3,3)
    real(kind=DP), allocatable :: dfsi_mpi(:,:,:,:) ! d(nval,nlmta,kv3,3)

    integer,save  :: id_sname = -1

    call tstatc0_begin('gather_drv_bdw_exx ',id_sname)

    dfsr_exx = 0.d0
    dfsi_exx = 0.d0

    do ik=1,kv3,af+1
       if(map_k(ik) /= myrank_k) cycle
       do ib=1,nval
          ib1 = neordr(ib,ik)
          if(map_e(ib1) == myrank_e) then
             dfsr_exx(ib,1:nlmta,ik,1:3) = dfsr_l(map_z(ib1),1:nlmta,ik,1:3)
             dfsi_exx(ib,1:nlmta,ik,1:3) = dfsi_l(map_z(ib1),1:nlmta,ik,1:3)
          end if
       end do
    end do

    if(npes>1) then
       allocate(dfsr_mpi(nval,nlmta,kv3,3))
       allocate(dfsi_mpi(nval,nlmta,kv3,3))
       call mpi_allreduce(dfsr_exx,dfsr_mpi,nval*nlmta*kv3*3,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
       dfsr_exx = dfsr_mpi
       call mpi_allreduce(dfsr_exx,dfsr_mpi,nval*nlmta*kv3*3,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ge_world,ierr)
       call mpi_allreduce(dfsi_exx,dfsi_mpi,nval*nlmta*kv3*3,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
       dfsi_exx = dfsi_mpi
       call mpi_allreduce(dfsi_exx,dfsi_mpi,nval*nlmta*kv3*3,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ge_world,ierr)
       dfsr_exx = dfsr_mpi
       dfsi_exx = dfsi_mpi
       deallocate(dfsr_mpi)
       deallocate(dfsi_mpi)
    end if

    call tstatc0_end(id_sname)
  end subroutine gather_drv_bdw_exx

  subroutine drv_betar_dot_WFs_exx()
    implicit none

    integer :: ik,ib,it,ia,i,il1,ilmtt1,ilmta1,iil,iksnl
    integer :: lmt1,lmtt1
    integer :: iadd, commsize
    real(kind=DP) :: fac, ph, fr, fi, fac2, vec1(3)
    real(kind=DP), allocatable :: cosgt(:) ! d(kg1)
    real(kind=DP), allocatable :: singt(:) ! d(kg1)
    real(kind=DP), allocatable :: gvec(:,:) ! d(kg1,3)

    integer,save  :: id_sname = -1
    call tstatc0_begin('drv_betar_dot_WFs_exx ',id_sname)

    allocate(cosgt(kg1))
    allocate(singt(kg1))
    allocate(gvec(kg1,3))

    dfsr_l = 0.d0
    dfsi_l = 0.d0

    do ik=1,kv3,af+1
       if(map_k(ik) /= myrank_k) cycle ! MPI
       iksnl = (ik-1)/nspin + 1
       fac2 = 1.d0
       if(k_symmetry(ik)==GAMMA) fac2 = 2.d0
       do ia=1,natm
          it = ityp(ia)
          if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
          do i=1,iba(ik)
             gvec(i,1:3) = vkxyz_for_outerloop(ik,1:3,BUCS)+ngabc(nbase(i,ik),1:3)
             ph = PAI2 * dot_product(pos(ia,1:3),gvec(i,1:3))
             cosgt(i) = cos(ph)
             singt(i) = sin(ph)
          end do
#ifdef DERIVATIVE_IN_CARTS
          do i=1,iba(ik)
             vec1(1) = rltv(1,1)*gvec(i,1) +rltv(1,2)*gvec(i,2) +rltv(1,3)*gvec(i,3)
             vec1(2) = rltv(2,1)*gvec(i,1) +rltv(2,2)*gvec(i,2) +rltv(2,3)*gvec(i,3)
             vec1(3) = rltv(3,1)*gvec(i,1) +rltv(3,2)*gvec(i,2) +rltv(3,3)*gvec(i,3)
             gvec(i,1:3) = vec1(1:3)
          end do
#endif
          do lmt1 = 1, ilmt(it)
             il1  = ltp(lmt1,it)
             lmtt1 = lmtt(lmt1,it)
             ilmta1 = lmta(lmt1,ia)
             iil=mod(il1,4) !! i^(L+1)
             if(iil==0) then
                do ib=1,np_e
                   if(occup_l(ib,ik) < DELTA) cycle
                   do i = ista_g1k(ik), iend_g1k(ik)
                      iadd = i - ista_g1k(ik) + 1
                      fac = snl(iadd,lmtt1,iksnl)*fac2
                      fr = fac * ( cosgt(i)*zaj_l(iadd,ib,ik,1)-singt(i)*zaj_l(iadd,ib,ik,2) )
                      fi = fac * ( cosgt(i)*zaj_l(iadd,ib,ik,2)+singt(i)*zaj_l(iadd,ib,ik,1) )
                      dfsr_l(ib,ilmta1,ik,1:3) = dfsr_l(ib,ilmta1,ik,1:3) + gvec(i,1:3) * fr
                      dfsi_l(ib,ilmta1,ik,1:3) = dfsi_l(ib,ilmta1,ik,1:3) + gvec(i,1:3) * fi
                   end do
                end do
             else if(iil==1) then
                do ib=1,np_e
                   if(occup_l(ib,ik) < DELTA) cycle
                   do i = ista_g1k(ik), iend_g1k(ik)
                      iadd = i - ista_g1k(ik) + 1
                      fac = snl(iadd,lmtt1,iksnl)*fac2
                      fr = -fac * ( cosgt(i)*zaj_l(iadd,ib,ik,2)+singt(i)*zaj_l(iadd,ib,ik,1) )
                      fi =  fac * ( cosgt(i)*zaj_l(iadd,ib,ik,1)-singt(i)*zaj_l(iadd,ib,ik,2) )
                      dfsr_l(ib,ilmta1,ik,1:3) = dfsr_l(ib,ilmta1,ik,1:3) + gvec(i,1:3) * fr
                      dfsi_l(ib,ilmta1,ik,1:3) = dfsi_l(ib,ilmta1,ik,1:3) + gvec(i,1:3) * fi
                   end do
                end do
             else if(iil==2) then
                do ib=1,np_e
                   if(occup_l(ib,ik) < DELTA) cycle
                   do i = ista_g1k(ik), iend_g1k(ik)
                      iadd = i - ista_g1k(ik) + 1
                      fac = snl(iadd,lmtt1,iksnl)*fac2
                      fr = -fac * ( cosgt(i)*zaj_l(iadd,ib,ik,1)-singt(i)*zaj_l(iadd,ib,ik,2) )
                      fi = -fac * ( cosgt(i)*zaj_l(iadd,ib,ik,2)+singt(i)*zaj_l(iadd,ib,ik,1) )
                      dfsr_l(ib,ilmta1,ik,1:3) = dfsr_l(ib,ilmta1,ik,1:3) + gvec(i,1:3) * fr
                      dfsi_l(ib,ilmta1,ik,1:3) = dfsi_l(ib,ilmta1,ik,1:3) + gvec(i,1:3) * fi
                   end do
                end do
             else if(iil==3) then
                do ib=1,np_e
                   if(occup_l(ib,ik) < DELTA) cycle
                   do i = ista_g1k(ik), iend_g1k(ik)
                      iadd = i - ista_g1k(ik) + 1
                      fac = snl(iadd,lmtt1,iksnl)*fac2
                      fr =  fac * ( cosgt(i)*zaj_l(iadd,ib,ik,2)+singt(i)*zaj_l(iadd,ib,ik,1) )
                      fi = -fac * ( cosgt(i)*zaj_l(iadd,ib,ik,1)-singt(i)*zaj_l(iadd,ib,ik,2) )
                      dfsr_l(ib,ilmta1,ik,1:3) = dfsr_l(ib,ilmta1,ik,1:3) + gvec(i,1:3) * fr
                      dfsi_l(ib,ilmta1,ik,1:3) = dfsi_l(ib,ilmta1,ik,1:3) + gvec(i,1:3) * fi
                   end do
                end do
             end if
          end do
       end do
    end do
    commsize = np_e*nlmta*(iend_k-ista_k+1)*3
    call mpi_allreduce(MPI_IN_PLACE,dfsr_l,commsize,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
    call mpi_allreduce(MPI_IN_PLACE,dfsi_l,commsize,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)

    deallocate(cosgt)
    deallocate(singt)
    deallocate(gvec)

    call tstatc0_end(id_sname)
  end subroutine drv_betar_dot_WFs_exx

! ==============================================================================
! === Make FFT box index arrays. ===============================================
! ==============================================================================
  subroutine Parallelize_wf_onto_fft_exx_3D(nfout,fft_box_size_WF,igf,nbase,nbase_gamma, &
 &                                         k_symmetry,GAMMA,kg,kg_gamma,kv3)
    integer, intent(in)  :: nfout, kg, kg_gamma, kv3, GAMMA
    integer, intent(in)  :: fft_box_size_WF(1:3,0:1)
    integer, intent(in)  :: igf(kg)
    integer, intent(in)  :: nbase(kg1_ext,kv3)
    integer, intent(in)  :: nbase_gamma(kg_gamma,2)
    integer, intent(in)  :: k_symmetry(kv3)
    integer, dimension(0:nrank_g-1)                       ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)       ::sta_r, sta_s

    integer, allocatable, dimension(:,:,:) :: xyz
    integer, allocatable, dimension(:,:,:,:,:) :: work
    integer :: max_fft_x, lx, ly, lz, mx, my, mz, mm, i, j, i1, klen, ik, max_np_g1k, max_mp_g1k
    integer :: iadd, ladd, ista, len, itag = 10
#ifndef USE_NONBLK_COMM
    integer, allocatable, dimension(:,:) :: rbuf
#endif
    integer :: itrs, iopr

    max_fft_x = maxval(nel_fft_x(:))
    allocate(xyz(2,3,0:nrank_g-1))
#ifdef USE_NONBLK_COMM
    req_r = 0
    req_s = 0
    sta_r = 0
    sta_s = 0
    do i = 0, nrank_g - 1
       call mpi_irecv(xyz(1,1,i), 6, mpi_integer, &
      &               i, itag, mpi_ke_world, req_r(i), ierr)
       if(ierr /= 0) then
          write(nfout,*)' Parallelize_wf_onto_fft_exx_3D :  mpi_irecv error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world,170,ierr)
       endif
    enddo
    do i = 0, nrank_g - 1
       call mpi_isend(xyz_fft_x, 6, mpi_integer, &
      &               i, itag, mpi_ke_world, req_s(i), ierr)
       if(ierr /= 0) then
          write(nfout,*)' Parallelize_wf_onto_fft_exx_3D :  mpi_isend error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world,171,ierr)
       endif
    enddo
    call mpi_waitall(nrank_g, req_r, sta_r, ierr)
    if(ierr /= 0) then
       write(nfout,*)' Parallelize_wf_onto_fft_exx_3D :  mpi_waitall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world,172,ierr)
    endif
    call mpi_waitall(nrank_g, req_s, sta_s, ierr)
    if(ierr /= 0) then
       write(nfout,*)' Parallelize_wf_onto_fft_exx_3D :  mpi_waitall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world,173,ierr)
    endif
#else
    allocate(rbuf(6,0:nrank_g-1))
    call MPI_ALLGATHER(xyz_fft_x, 6, mpi_integer, &
   &                   rbuf,6, mpi_integer, mpi_ke_world, ierr )
    if(ierr /= 0) then
       write(nfout,*)' Parallelize_wf_onto_fft_exx_3D :  mpi_allgather error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 174, ierr)
    endif
    do i = 1,2
       do j = 1,3
          xyz(i,j,:)=rbuf(i+(j-1)*2,:)
       enddo
    enddo
    deallocate(rbuf)
#endif

    lx = fft_box_size_WF(1,0)
    ly = fft_box_size_WF(2,0)
    lz = fft_box_size_WF(3,0)
    allocate(wf_fft_scnt_exx(0:nrank_g-1,kv3,nopr,0:ntrs))
    allocate(wf_fft_rcnt_exx(0:nrank_g-1,kv3,nopr,0:ntrs))
    len = 1
    do ik = 1, kv3
       if(k_symmetry(ik) == GAMMA) then
          len = 2
       end if
    end do
    max_np_g1k = maxval(np_g1k(:))
    max_mp_g1k = maxval(mp_g1k(:))
    allocate(wf_fft_index_exx(max_np_g1k*len,kv3,nopr,0:ntrs))
    allocate(wf_fft_dist_exx (max_np_g1k*len,kv3,nopr,0:ntrs))
    allocate(wf_fft_send_exx (max_mp_g1k*len,kv3,0:nrank_g-1,nopr,0:ntrs))
    allocate(wf_fft_recv_exx (max_mp_g1k*len,kv3,0:nrank_g-1,nopr,0:ntrs))
    allocate(wf_fft_maxsend_exx(kv3,nopr,0:ntrs))
    allocate(wf_fft_maxrecv_exx(kv3,nopr,0:ntrs))
    wf_fft_dist_exx = -1
    wf_fft_send_exx = 0
    wf_fft_recv_exx = 0
    wf_fft_scnt_exx = 0
    wf_fft_rcnt_exx = 0
    wf_fft_maxsend_exx = 0
    wf_fft_maxrecv_exx = 0

    klen = kv3

! ==============================================================================
    do itrs = 0, ntrs
    do iopr = 1, nopr
! ==============================================================================
       do ik = 1, kv3
          if(k_symmetry(ik) == GAMMA) then
             ista = ista_g1k(ik)
             if(ista == 1) then
                iadd = 1
                i1 = igf(1)
                mz = (i1-1)/(lx*ly)+1
                mm = mod(i1,(lx*ly))
                if(mm==0) mm=lx*ly
                my = (mm-1)/lx+1
                mx = mod(mm,lx)
                if(mx==0) mx = ly
          B_4 : do i = 0, nrank_g-1
                   if((xyz(1,2,i)<=my).and.(my<=xyz(2,2,i)) .and. &
                  &   (xyz(1,3,i)<=mz).and.(mz<=xyz(2,3,i))) then
                      ladd = mx+lx*(my-xyz(1,2,i))+lx*(xyz(2,2,i)-xyz(1,2,i)+1)*(mz-xyz(1,3,i))
                      wf_fft_scnt_exx(i,ik,iopr,itrs) = wf_fft_scnt_exx(i,ik,iopr,itrs) + 1
                      wf_fft_index_exx(iadd*2-1,ik,iopr,itrs) = wf_fft_scnt_exx(i,ik,iopr,itrs)
                      wf_fft_dist_exx (iadd*2-1,ik,iopr,itrs) = i
                      wf_fft_send_exx(wf_fft_scnt_exx(i,ik,iopr,itrs),ik,i,iopr,itrs) = ladd
                      wf_fft_scnt_exx(i,ik,iopr,itrs) = wf_fft_scnt_exx(i,ik,iopr,itrs) + 1
                      wf_fft_index_exx(iadd*2  ,ik,iopr,itrs) = wf_fft_scnt_exx(i,ik,iopr,itrs)
                      wf_fft_dist_exx (iadd*2  ,ik,iopr,itrs) = i
                      wf_fft_send_exx(wf_fft_scnt_exx(i,ik,iopr,itrs),ik,i,iopr,itrs) = ladd
                      exit
                   endif
                enddo B_4
                ista = 2
             endif
       B_1 : do j = ista, iend_g1k(ik)
                iadd = j-ista_g1k(ik)+1
                i1 = igf(ngpt_exx(nbase(j,ik),iopr,itrs))
                mz = (i1-1)/(lx*ly)+1
                mm = mod(i1,(lx*ly))
                if(mm==0) mm=lx*ly
                my = (mm-1)/lx+1
                mx = mod(mm,lx)
                if(mx==0) mx = lx
          B_2 : do i = 0, nrank_g-1
                   if((xyz(1,2,i)<=my).and.(my<=xyz(2,2,i)) .and. &
                  &   (xyz(1,3,i)<=mz).and.(mz<=xyz(2,3,i))) then
                      ladd = mx+lx*(my-xyz(1,2,i))+lx*(xyz(2,2,i)-xyz(1,2,i)+1)*(mz-xyz(1,3,i))
                      wf_fft_scnt_exx(i,ik,iopr,itrs) = wf_fft_scnt_exx(i,ik,iopr,itrs) + 1
                      wf_fft_index_exx(iadd*2-1,ik,iopr,itrs) = wf_fft_scnt_exx(i,ik,iopr,itrs)
                      wf_fft_dist_exx (iadd*2-1,ik,iopr,itrs) = i
                      wf_fft_send_exx(wf_fft_scnt_exx(i,ik,iopr,itrs),ik,i,iopr,itrs) = ladd
                      exit
                   endif
                enddo B_2
                i1 = igf(ngpt_exx(nbase_gamma(j,2),iopr,itrs))
                mz = (i1-1)/(lx*ly)+1
                mm = mod(i1,(lx*ly))
                if(mm==0) mm=lx*ly
                my = (mm-1)/lx+1
                mx = mod(mm,lx)
                if(mx==0) mx = lx
          B_3 : do i = 0, nrank_g-1
                   if((xyz(1,2,i)<=my).and.(my<=xyz(2,2,i)) .and. &
                  &   (xyz(1,3,i)<=mz).and.(mz<=xyz(2,3,i))) then
                      ladd = mx+lx*(my-xyz(1,2,i))+lx*(xyz(2,2,i)-xyz(1,2,i)+1)*(mz-xyz(1,3,i))
                      wf_fft_scnt_exx(i,ik,iopr,itrs) = wf_fft_scnt_exx(i,ik,iopr,itrs) + 1
                      wf_fft_index_exx(iadd*2  ,ik,iopr,itrs) = wf_fft_scnt_exx(i,ik,iopr,itrs)
                      wf_fft_dist_exx (iadd*2  ,ik,iopr,itrs) = i
                      wf_fft_send_exx(wf_fft_scnt_exx(i,ik,iopr,itrs),ik,i,iopr,itrs) = ladd
                      exit
                   endif
                enddo B_3
             enddo B_1
          else
!!$             write(nfout,'(" ista_g1k, iend_g1k for " ,i8," = ",2i8)') ik, ista_g1k(ik), iend_g1k(ik)
!!$             write(nfout,'(" iopr, itrs = ",2i8)') iopr, itrs
             call flush(6)
      B_11 : do j = ista_g1k(ik), iend_g1k(ik)
                iadd = j-ista_g1k(ik)+1
                i1 = igf(ngpt_exx(nbase(j,ik),iopr,itrs))
                mz = (i1-1)/(lx*ly)+1
                mm = mod(i1,(lx*ly))
                if(mm==0) mm=lx*ly
                my = (mm-1)/lx+1
                mx = mod(mm,lx)
                if(mx==0) mx = lx
         B_12 : do i = 0, nrank_g-1
                   if((xyz(1,2,i)<=my).and.(my<=xyz(2,2,i)) .and. &
                  &   (xyz(1,3,i)<=mz).and.(mz<=xyz(2,3,i))) then
                      ladd = mx+lx*(my-xyz(1,2,i))+lx*(xyz(2,2,i)-xyz(1,2,i)+1)*(mz-xyz(1,3,i))
                      wf_fft_scnt_exx(i,ik,iopr,itrs) = wf_fft_scnt_exx(i,ik,iopr,itrs) + 1
                      wf_fft_index_exx(iadd,ik,iopr,itrs) = wf_fft_scnt_exx(i,ik,iopr,itrs)
                      wf_fft_dist_exx (iadd,ik,iopr,itrs) = i
                      wf_fft_send_exx(wf_fft_scnt_exx(i,ik,iopr,itrs),ik,i,iopr,itrs) = ladd
                      exit
                   endif
                enddo B_12
             enddo B_11
          endif
       end do
! ==============================================================================
    end do
    end do
! ==============================================================================
    deallocate(xyz)

! ==============================================================================
    do itrs = 0, ntrs
    do iopr = 1, nopr
! ==============================================================================
#ifdef USE_NONBLK_COMM
       req_r = 0
       req_s = 0
       sta_r = 0
       sta_s = 0
       do i = 0, nrank_g - 1
          call mpi_irecv(wf_fft_recv_exx(1,1,i,iopr,itrs), max_mp_g1k*len*klen, mpi_integer, &
         &               i, itag, mpi_ke_world, req_r(i), ierr)
          if(ierr /= 0) then
             write(nfout,*)' Parallelize_wf_onto_fft_exx_3D :  mpi_irecv error'
             call flush(nfout)
             call mpi_abort(mpi_comm_world,174,ierr)
          endif
       enddo
       do i = 0, nrank_g - 1
          call mpi_isend(wf_fft_send_exx(1,1,i,iopr,itrs), max_mp_g1k*len*klen, mpi_integer, &
         &               i, itag, mpi_ke_world, req_s(i), ierr)
          if(ierr /= 0) then
             write(nfout,*)' Parallelize_wf_onto_fft_exx_3D :  mpi_isend error'
             call flush(nfout)
             call mpi_abort(mpi_comm_world,175,ierr)
          endif
       enddo
       call mpi_waitall(nrank_g, req_r, sta_r, ierr)
       if(ierr /= 0) then
          write(nfout,*)' Parallelize_wf_onto_fft_exx_3D :  mpi_waitall error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world,176,ierr)
       endif
       call mpi_waitall(nrank_g, req_s, sta_s, ierr)
       if(ierr /= 0) then
          write(nfout,*)' Parallelize_wf_onto_fft_exx_3D :  mpi_waitall error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world,177,ierr)
       endif
#else
       call MPI_ALLTOALL(wf_fft_send_exx(1,1,0,iopr,itrs), max_mp_g1k*len*klen, mpi_integer, &
      &                  wf_fft_recv_exx(1,1,0,iopr,itrs), max_mp_g1k*len*klen, mpi_integer, &
      &                  mpi_ke_world, ierr)
       if(ierr /= 0) then
          write(nfout,*)' Parallelize_wf_onto_fft_exx_3D :  mpi_alltoall error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 178, ierr)
       endif
#endif

       do i = 0, nrank_g - 1
          do ik = 1, kv3
             do j = 1, mp_g1k(ik)*len
                if(wf_fft_recv_exx(j,ik,i,iopr,itrs) == 0) then
                   exit
                end if
                wf_fft_rcnt_exx(i,ik,iopr,itrs) = wf_fft_rcnt_exx(i,ik,iopr,itrs) + 1
             enddo
          enddo
       enddo

       do ik = 1, kv3
          wf_fft_maxsend_exx(ik,iopr,itrs) = maxval(wf_fft_scnt_exx(0:nrank_g-1,ik,iopr,itrs))
          wf_fft_maxrecv_exx(ik,iopr,itrs) = maxval(wf_fft_rcnt_exx(0:nrank_g-1,ik,iopr,itrs))
       end do
! ==============================================================================
    end do
    end do
! ==============================================================================

    allocate(work(maxval(wf_fft_maxrecv_exx(:,:,:)),kv3,0:nrank_g-1,nopr,0:ntrs))
! ==============================================================================
    do itrs = 0, ntrs
    do iopr = 1, nopr
! ==============================================================================
       do i = 0, nrank_g - 1
          do ik = 1, kv3
             do j = 1, wf_fft_rcnt_exx(i,ik,iopr,itrs)
                work(j,ik,i,iopr,itrs) = wf_fft_recv_exx(j,ik,i,iopr,itrs)
             end do
          end do
       end do
! ==============================================================================
    end do
    end do
! ==============================================================================
    deallocate(wf_fft_recv_exx)
    allocate(wf_fft_recv_exx(maxval(wf_fft_maxrecv_exx(:,:,:)),kv3,0:nrank_g-1,nopr,0:ntrs))
! ==============================================================================
    do itrs = 0, ntrs
    do iopr = 1, nopr
! ==============================================================================
       do i = 0, nrank_g - 1
          do ik = 1, kv3
             do j = 1, wf_fft_rcnt_exx(i,ik,iopr,itrs)
                wf_fft_recv_exx(j,ik,i,iopr,itrs) = work(j,ik,i,iopr,itrs)
             end do
          end do
       end do
! ==============================================================================
    end do
    end do
! ==============================================================================
    deallocate(work)

    deallocate(wf_fft_send_exx)
  end subroutine Parallelize_wf_onto_fft_exx_3D

  subroutine m_Parallel_wf_onto_fft_dealloc_exx_3D
    deallocate(wf_fft_rcnt_exx)
    deallocate(wf_fft_scnt_exx)
    deallocate(wf_fft_recv_exx)
    deallocate(wf_fft_index_exx)
    deallocate(wf_fft_dist_exx)
    deallocate(wf_fft_maxrecv_exx)
    deallocate(wf_fft_maxsend_exx)
  end subroutine m_Parallel_wf_onto_fft_dealloc_exx_3D

  subroutine Parallelize_fft_onto_wf_rhog_3D(nfout,igf,kg,nfft)
    integer, intent(in)  :: nfout, kg, nfft
    integer, intent(in)  :: igf(kg)
    integer, dimension(0:nrank_g-1)                       ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)       ::sta_r, sta_s

    integer, allocatable,dimension(:,:) :: fftigf
    integer, allocatable, dimension(:,:) :: work
    integer :: i1, lrank, i, j, k, lsize, isrsize,fft_l_size, klen
    integer, parameter :: itag = 10

    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
    isrsize = min(lsize,maxval(nel_kngp(:)))
    fft_l_size  = nel_fft_x(myrank_g)

    allocate(fft_wf_scnt_rhog(0:nrank_g-1), stat=ierr)
    allocate(fft_wf_rcnt_rhog(0:nrank_g-1), stat=ierr)
    allocate(fft_wf_send_rhog(isrsize,0:nrank_g-1), stat=ierr)
    allocate(fft_wf_recv_rhog(isrsize,0:nrank_g-1), stat=ierr)
    allocate(fft_wf_dist_rhog(fft_l_size), stat=ierr)
    allocate(fft_wf_index_rhog(fft_l_size), stat=ierr)
    if(ierr /= 0) then
       write(nfout,*)' Parallelize_fft_onto_wf_rhog_3D :  Not allocate '
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 204, ierr)
    endif
    fft_wf_scnt_rhog = 0
    fft_wf_rcnt_rhog = 0
    fft_wf_send_rhog = 0
    fft_wf_recv_rhog = 0
    fft_wf_dist_rhog = -1
    fft_wf_index_rhog = 0

    allocate(fftigf(nfft,2), stat=ierr)
    if(ierr /= 0) then
       write(nfout,*)' Parallelize_fft_onto_wf_rhog_3D :  Not allocate '
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 206, ierr)
    endif

    fftigf(:,1) = -1

    do i = 0, nrank_g - 1
       do j = is_kngp(i), min(kg,ie_kngp(i))
         i1 = igf(j)
         if(i1 > nfft) cycle
         fftigf(i1,1) = i
         fftigf(i1,2) = j - is_kngp(i) + 1
       enddo
    enddo
    do k = 1, nel_fft_x(myrank_g)
       i1 = mp_fft_x(k)
       if(fftigf(i1,1) < 0) cycle
       lrank = fftigf(i1,1)
       fft_wf_scnt_rhog(lrank) = fft_wf_scnt_rhog(lrank) + 1
       fft_wf_send_rhog(fft_wf_scnt_rhog(lrank),lrank) = fftigf(i1,2)
       fft_wf_dist_rhog(k) = lrank
       fft_wf_index_rhog(k) = fft_wf_scnt_rhog(lrank)
    end do

    deallocate(fftigf)

#ifdef USE_NONBLK_COMM
    lrank = mod(myrank_g,nrank_g)
    do i = 0, nrank_g - 1
       lrank = lrank + 1
       if(lrank > (nrank_g - 1)) lrank = 0
       call mpi_irecv(fft_wf_recv_rhog(1,lrank), isrsize, mpi_integer, &
      &               lrank, itag, mpi_ke_world, req_r(lrank), ierr)
       if(ierr /= 0) then
          write(nfout,*)' Parallelize_fft_onto_wf_rhog_3D :  mpi_irecv error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 207, ierr)
       endif
    enddo

    lrank = mod((myrank_g+1),nrank_g)
    do i = 0, nrank_g - 1
       lrank = lrank + 1
       if(lrank > (nrank_g - 1)) lrank = 0
       call mpi_isend(fft_wf_send_rhog(1,lrank), isrsize, mpi_integer, &
      &               lrank, itag, mpi_ke_world, req_s(lrank), ierr)
       if(ierr /= 0) then
          write(nfout,*)' Parallelize_fft_onto_wf_rhog_3D :  mpi_isend error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 208, ierr)
       endif
    enddo

    call mpi_waitall(nrank_g, req_r, sta_r, ierr)
    if(ierr /= 0) then
       write(nfout,*)' Parallelize_fft_onto_wf_rhog_3D :  mpi_waitall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 209, ierr)
    endif
    call mpi_waitall(nrank_g, req_s, sta_s, ierr)
    if(ierr /= 0) then
       write(nfout,*)' Parallelize_fft_onto_wf_rhog_3D :  mpi_waitall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 210, ierr)
    endif
#else
    call MPI_ALLTOALL(fft_wf_send_rhog, isrsize, mpi_integer, &
   &                  fft_wf_recv_rhog, isrsize, mpi_integer, &
   &                  mpi_ke_world, ierr)
    if(ierr /= 0) then
       write(nfout,*)' Parallelize_fft_onto_wf_rhog_3D :  mpi_alltoall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 211, ierr)
    endif
#endif

    do i = 0, nrank_g - 1
       do j = 1, isrsize
          if(fft_wf_recv_rhog(j,i) == 0) then
             exit
          end if
          fft_wf_rcnt_rhog(i) = fft_wf_rcnt_rhog(i) + 1
       enddo
    enddo

    fft_wf_maxsend_rhog = maxval(fft_wf_scnt_rhog(:))
    fft_wf_maxrecv_rhog = maxval(fft_wf_rcnt_rhog(:))

    allocate(work(fft_wf_maxrecv_rhog,0:nrank_g-1))
    do i = 0, nrank_g - 1
       do j = 1, fft_wf_rcnt_rhog(i)
          work(j,i) = fft_wf_recv_rhog(j,i)
       end do
    end do
    deallocate(fft_wf_recv_rhog)
    allocate(fft_wf_recv_rhog(fft_wf_maxrecv_rhog,0:nrank_g-1))
    do i = 0, nrank_g - 1
       do j = 1, fft_wf_rcnt_rhog(i)
          fft_wf_recv_rhog(j,i) = work(j,i)
       end do
    end do
    deallocate(work)
    deallocate(fft_wf_send_rhog)
  end subroutine Parallelize_fft_onto_wf_rhog_3D

  subroutine m_Parallel_fft_onto_wf_dealloc_rhog_3D
    deallocate(fft_wf_rcnt_rhog)
    deallocate(fft_wf_scnt_rhog)
    deallocate(fft_wf_recv_rhog)
    deallocate(fft_wf_index_rhog)
    deallocate(fft_wf_dist_rhog)
  end subroutine m_Parallel_fft_onto_wf_dealloc_rhog_3D

  subroutine Parallelize_wf_onto_fft_rhog_3D(nfout,fft_box_size_WF,igf,kg)
    integer, intent(in)  :: nfout, kg
    integer, intent(in)  :: fft_box_size_WF(1:3,0:1)
    integer, intent(in)  :: igf(kg)
    integer, dimension(0:nrank_g-1)                       ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)       ::sta_r, sta_s

    integer, allocatable, dimension(:,:,:) :: xyz
    integer, allocatable, dimension(:,:) :: work
    integer :: max_fft_x, lx, ly, lz, mx, my, mz, mm, i, j, i1, klen, max_np_g1k, max_mp_g1k
    integer :: iadd, ladd, ista, len, itag = 10
#ifndef USE_NONBLK_COMM
    integer, allocatable, dimension(:,:) :: rbuf
#endif

    max_fft_x = maxval(nel_fft_x(:))
    allocate(xyz(2,3,0:nrank_g-1))
#ifdef USE_NONBLK_COMM
    req_r = 0
    req_s = 0
    sta_r = 0
    sta_s = 0
    do i = 0, nrank_g - 1
       call mpi_irecv(xyz(1,1,i), 6, mpi_integer, &
      &               i, itag, mpi_ke_world, req_r(i), ierr)
       if(ierr /= 0) then
          write(nfout,*)' Parallelize_wf_onto_fft_rhog_3D :  mpi_irecv error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world,170,ierr)
       endif
    enddo
    do i = 0, nrank_g - 1
       call mpi_isend(xyz_fft_x, 6, mpi_integer, &
      &               i, itag, mpi_ke_world, req_s(i), ierr)
       if(ierr /= 0) then
          write(nfout,*)' Parallelize_wf_onto_fft_rhog_3D :  mpi_isend error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world,171,ierr)
       endif
    enddo
    call mpi_waitall(nrank_g, req_r, sta_r, ierr)
    if(ierr /= 0) then
       write(nfout,*)' Parallelize_wf_onto_fft_rhog_3D :  mpi_waitall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world,172,ierr)
    endif
    call mpi_waitall(nrank_g, req_s, sta_s, ierr)
    if(ierr /= 0) then
       write(nfout,*)' Parallelize_wf_onto_fft_rhog_3D :  mpi_waitall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world,173,ierr)
    endif
#else
    allocate(rbuf(6,0:nrank_g-1))
    call MPI_ALLGATHER(xyz_fft_x, 6, mpi_integer, &
   &                   rbuf, 6, mpi_integer, mpi_ke_world, ierr )
    if(ierr /= 0) then
       write(nfout,*)' Parallelize_wf_onto_fft_rhog_3D :  mpi_allgather error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 174, ierr)
    endif
    do i = 1,2
       do j = 1,3
          xyz(i,j,:)=rbuf(i+(j-1)*2,:)
       enddo
    enddo
    deallocate(rbuf)
#endif

    lx = fft_box_size_WF(1,0)
    ly = fft_box_size_WF(2,0)
    lz = fft_box_size_WF(3,0)
    allocate(wf_fft_scnt_rhog(0:nrank_g-1))
    allocate(wf_fft_rcnt_rhog(0:nrank_g-1))
    len = 1
    max_np_g1k = np_kngp
    max_mp_g1k = maxval(nel_kngp(:))
    allocate(wf_fft_index_rhog(max_np_g1k*len))
    allocate(wf_fft_dist_rhog (max_np_g1k*len))
    allocate(wf_fft_send_rhog (max_mp_g1k*len,0:nrank_g-1))
    allocate(wf_fft_recv_rhog (max_mp_g1k*len,0:nrank_g-1))
    wf_fft_dist_rhog = -1
    wf_fft_send_rhog = 0
    wf_fft_recv_rhog = 0
    wf_fft_scnt_rhog = 0
    wf_fft_rcnt_rhog = 0
    wf_fft_maxsend_rhog = 0
    wf_fft_maxrecv_rhog = 0

B_11 : do j = ista_kngp, min(kg,iend_kngp)
          iadd = j-ista_kngp+1
          i1 = igf(j)
          mz = (i1-1)/(lx*ly)+1
          mm = mod(i1,(lx*ly))
          if(mm==0) mm=lx*ly
          my = (mm-1)/lx+1
          mx = mod(mm,lx)
          if(mx==0) mx = lx
   B_12 : do i = 0, nrank_g-1
             if((xyz(1,2,i)<=my).and.(my<=xyz(2,2,i)) .and. &
            &   (xyz(1,3,i)<=mz).and.(mz<=xyz(2,3,i))) then
                ladd = mx+lx*(my-xyz(1,2,i))+lx*(xyz(2,2,i)-xyz(1,2,i)+1)*(mz-xyz(1,3,i))
                wf_fft_scnt_rhog(i) = wf_fft_scnt_rhog(i) + 1
                wf_fft_index_rhog(iadd) = wf_fft_scnt_rhog(i)
                wf_fft_dist_rhog (iadd) = i
                wf_fft_send_rhog(wf_fft_scnt_rhog(i),i) = ladd
                exit
             endif
          enddo B_12
       enddo B_11
    deallocate(xyz)

#ifdef USE_NONBLK_COMM
    req_r = 0
    req_s = 0
    sta_r = 0
    sta_s = 0
    do i = 0, nrank_g - 1
       call mpi_irecv(wf_fft_recv_rhog(1,i), max_mp_g1k*len, mpi_integer, &
      &               i, itag, mpi_ke_world, req_r(i), ierr)
       if(ierr /= 0) then
          write(nfout,*)' Parallelize_wf_onto_fft_rhog_3D :  mpi_irecv error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world,174,ierr)
       endif
    enddo
    do i = 0, nrank_g - 1
       call mpi_isend(wf_fft_send_rhog(1,i), max_mp_g1k*len, mpi_integer, &
      &               i, itag, mpi_ke_world, req_s(i), ierr)
       if(ierr /= 0) then
          write(nfout,*)' Parallelize_wf_onto_fft_rhog_3D :  mpi_isend error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world,175,ierr)
       endif
    enddo
    call mpi_waitall(nrank_g, req_r, sta_r, ierr)
    if(ierr /= 0) then
       write(nfout,*)' Parallelize_wf_onto_fft_rhog_3D :  mpi_waitall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world,176,ierr)
    endif
    call mpi_waitall(nrank_g, req_s, sta_s, ierr)
    if(ierr /= 0) then
       write(nfout,*)' Parallelize_wf_onto_fft_rhog_3D :  mpi_waitall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world,177,ierr)
    endif
#else
    call MPI_ALLTOALL(wf_fft_send_rhog, max_mp_g1k*len, mpi_integer, &
   &                  wf_fft_recv_rhog, max_mp_g1k*len, mpi_integer, &
   &                  mpi_ke_world, ierr)
    if(ierr /= 0) then
       write(nfout,*)' Parallelize_wf_onto_fft_rhog_3D :  mpi_alltoall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 178, ierr)
    endif
#endif

    do i = 0, nrank_g - 1
       do j = 1, np_kngp*len
          if(wf_fft_recv_rhog(j,i) == 0) then
             exit
          end if
          wf_fft_rcnt_rhog(i) = wf_fft_rcnt_rhog(i) + 1
       enddo
    enddo

    wf_fft_maxsend_rhog = maxval(wf_fft_scnt_rhog(0:nrank_g-1))
    wf_fft_maxrecv_rhog = maxval(wf_fft_rcnt_rhog(0:nrank_g-1))

    allocate(work(wf_fft_maxrecv_rhog,0:nrank_g-1))
    do i = 0, nrank_g - 1
       do j = 1, wf_fft_rcnt_rhog(i)
          work(j,i) = wf_fft_recv_rhog(j,i)
       end do
    end do
    deallocate(wf_fft_recv_rhog)
    allocate(wf_fft_recv_rhog(wf_fft_maxrecv_rhog,0:nrank_g-1))
    do i = 0, nrank_g - 1
       do j = 1, wf_fft_rcnt_rhog(i)
          wf_fft_recv_rhog(j,i) = work(j,i)
       end do
    end do
    deallocate(work)

    deallocate(wf_fft_send_rhog)
  end subroutine Parallelize_wf_onto_fft_rhog_3D

  subroutine m_Parallel_wf_onto_fft_dealloc_rhog_3D
    deallocate(wf_fft_rcnt_rhog)
    deallocate(wf_fft_scnt_rhog)
    deallocate(wf_fft_recv_rhog)
    deallocate(wf_fft_index_rhog)
    deallocate(wf_fft_dist_rhog)
  end subroutine m_Parallel_wf_onto_fft_dealloc_rhog_3D

! ==============================================================================
! === Convert FFT box <-> kg/kg1 ===============================================
! ==============================================================================
  subroutine map_WFG_on_FFT_box_3D(ik,ib1,ib2,ibsize,lsize,bfft_l,wfr,wfi)
    use m_Parallelization,     only : wf_fft_scnt, wf_fft_rcnt &
   &                                , wf_fft_recv &
   &                                , wf_fft_index, wf_fft_dist &
   &                                , wf_fft_maxrecv, wf_fft_maxsend
    integer, intent(in) :: ik, ib1,ib2, ibsize,lsize
    real(kind=DP), intent(inout), dimension(lsize*kimg,ibsize) :: bfft_l
    real(kind=DP), intent(in), dimension(maxval(np_g1k)) :: wfr, wfi
    integer :: i, j, k, ii, jj, iesize, iadd, ierr
    integer, dimension(0:nrank_g-1) ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)::sta_r, sta_s
    real(kind=DP), allocatable, dimension(:,:),save :: sendbuf, recvbuf
    integer :: icnt_send, icnt_recv, lrank
    integer, parameter :: itag = 21
    integer, save :: savesize = 0, savesend=0, saverecv=0
    integer :: nx, ny, nz, nxp, nn, ix, iy, iz
    integer,save  :: id_sname = -1
#ifdef USE_ALLTOALLV
    integer, allocatable, dimension(:) :: sdsp, rdsp
#endif
    call tstatc0_begin('map_WFG_on_FFT_box_3D ',id_sname)

    iesize = ib2 - ib1 + 1

    if(allocated(sendbuf)) deallocate(sendbuf)
    allocate(sendbuf(wf_fft_maxsend(ik)*kimg*iesize,0:nrank_g-1))
    if(allocated(recvbuf)) deallocate(recvbuf)
    allocate(recvbuf(wf_fft_maxrecv(ik)*kimg*iesize,0:nrank_g-1))

    savesize = iesize
    sendbuf = 0.0d0
    recvbuf = 0.0d0

    if(k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                sendbuf(iesize*(wf_fft_index(ii*2-1,ik)-1)+1,wf_fft_dist(ii*2-1,ik)) = wfr(ii)
                sendbuf(iesize*(wf_fft_index(ii*2  ,ik)-1)+1,wf_fft_dist(ii*2  ,ik)) = wfr(ii)
             enddo
          enddo
       else
!OCL NORECURRENCE
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                iadd = iesize*2*(wf_fft_index(ii*2-1,ik)-1)+jj*2
                sendbuf(iadd-1,wf_fft_dist(ii*2-1,ik)) = wfr(ii)
                sendbuf(iadd,  wf_fft_dist(ii*2-1,ik)) = wfi(ii)
                iadd = iesize*2*(wf_fft_index(ii*2,ik)-1)+jj*2
                sendbuf(iadd-1,wf_fft_dist(ii*2  ,ik)) =  wfr(ii)
                sendbuf(iadd,  wf_fft_dist(ii*2  ,ik)) = -wfi(ii)
             enddo
          enddo
       endif
    else
       if(kimg == 1) then
!OCL NORECURRENCE
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                iadd = iesize*(wf_fft_index(ii,ik)-1)+jj
                sendbuf(iadd,wf_fft_dist(ii,ik)) = wfr(ii)
             enddo
          enddo
       else
!OCL NORECURRENCE
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                iadd = iesize*2*(wf_fft_index(ii,ik)-1)+jj*2
                sendbuf(iadd-1,wf_fft_dist(ii,ik)) = wfr(ii)
                sendbuf(iadd  ,wf_fft_dist(ii,ik)) = wfi(ii)
             enddo
          enddo
       endif
    endif

#ifndef USE_ALLTOALLV
    icnt_recv = 0
    do i = 0, nrank_g - 1
       if(wf_fft_rcnt(i,ik) /= 0) then
          call mpi_irecv(recvbuf(1,i), wf_fft_rcnt(i,ik)*kimg*iesize, mpi_double_precision, &
         &               i, itag, mpi_ke_world, req_r(icnt_recv), ierr)
          if(ierr /= 0) then
             write(nfout,*)' map_WFG_on_FFT_box_3D :  mpi_irecv error'
             call flush(nfout)
             call mpi_abort(mpi_comm_world,166,ierr)
          endif
          icnt_recv = icnt_recv + 1
       endif
    enddo
    icnt_send = 0
    do i = 0, nrank_g - 1
       if(wf_fft_scnt(i,ik) /= 0) then
          call mpi_isend(sendbuf(1,i), wf_fft_scnt(i,ik)*kimg*iesize, mpi_double_precision, &
         &               i, itag, mpi_ke_world, req_s(icnt_send), ierr)
          if(ierr /= 0) then
             write(nfout,*)' map_WFG_on_FFT_box_3D :  mpi_isend error'
             call flush(nfout)
             call mpi_abort(mpi_comm_world,167,ierr)
          endif
          icnt_send = icnt_send + 1
       endif
    enddo
    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
    if(ierr /= 0) then
       write(nfout,*)' map_WFG_on_FFT_box_3D :  mpi_waitall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world,168,ierr)
    endif
    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
    if(ierr /= 0) then
       write(nfout,*)' map_WFG_on_FFT_box_3D :  mpi_waitall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world,169,ierr)
    endif
#else
    allocate(sdsp(0:nrank_g-1), stat=ierr)
    allocate(rdsp(0:nrank_g-1), stat=ierr)
    do i = 0, nrank_g - 1
       sdsp(i)=wf_fft_maxsend(ik)*kimg*iesize*i
       rdsp(i)=wf_fft_maxrecv(ik)*kimg*iesize*i
    enddo
    call MPI_ALLTOALLV(      sendbuf, wf_fft_scnt(:,ik)*kimg*iesize, sdsp,&
   &   mpi_double_precision, recvbuf, wf_fft_rcnt(:,ik)*kimg*iesize, rdsp,&
   &   mpi_double_precision, mpi_ke_world, ierr)
    if(ierr /= 0) then
       write(nfout,*)' m_ES_Vlocal_in_Rspace_3D : mpi_alltoallv error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 170, ierr)
    endif
    deallocate(sdsp)
    deallocate(rdsp)
#endif

    bfft_l = 0.0d0
#ifdef FFT_USE_SSL2
    nx = xyz_fft_x(2,1) - xyz_fft_x(1,1) + 1
    ny = xyz_fft_x(2,2) - xyz_fft_x(1,2) + 1
    nz = xyz_fft_x(2,3) - xyz_fft_x(1,3) + 1
    nxp = nx

    if(kimg == 1) then
!OCL NORECURRENCE
!OCL NOFLTLD
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             iadd = wf_fft_recv(j,ik,i)
             iz = (iadd-1)/(nx*ny)+1
             nn = mod(iadd,nx*ny)
             if(nn==0) then
                iy = ny
             else
                iy = (nn-1)/nx+1
             end if
             ix = mod(nn,nx)
             if(ix==0) ix = nx
             do k = 1, iesize
                bfft_l(ix+(iy-1)*nx+(iz-1)*nx*ny,k) = recvbuf(iesize*(j-1)+k,i)
             enddo
          enddo
       enddo
    else
!OCL NORECURRENCE
!OCL NOFLTLD
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             iadd = wf_fft_recv(j,ik,i)
             iz = (iadd-1)/(nx*ny)+1
             nn = mod(iadd,nx*ny)
             if(nn==0) then
                iy = ny
             else
                iy = (nn-1)/nx+1
             end if
             ix = mod(nn,nx)
             if(ix==0) ix = nx
             do k = 1, iesize
                bfft_l((ix+(iy-1)*nxp+(iz-1)*nxp*ny)*2-1,k) = recvbuf(iesize*2*(j-1)+k*2-1,i)
                bfft_l((ix+(iy-1)*nxp+(iz-1)*nxp*ny)*2  ,k) = recvbuf(iesize*2*(j-1)+k*2  ,i)
             enddo
          enddo
       enddo
    endif
#else
    if(kimg == 1) then
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             do k = 1, iesize
                bfft_l(wf_fft_recv(j,ik,i),k) = recvbuf(iesize*(j-1)+k,i)
             enddo
          enddo
       enddo
    else
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             do k = 1, iesize
                bfft_l(wf_fft_recv(j,ik,i)*2-1,k) = recvbuf(iesize*2*(j-1)+k*2-1,i)
                bfft_l(wf_fft_recv(j,ik,i)*2  ,k) = recvbuf(iesize*2*(j-1)+k*2  ,i)
             enddo
          enddo
       enddo
    endif
#endif

! === FFT xzy ==================================================================
    if(sw_fft_xzy > 0) then
! ==============================================================================
    call m_FFT_Inverse_3D(nfout, bfft_l, lsize, iesize)
! === FFT xzy ==================================================================
    else
    call m_FFT_Inverse_XYZ_3D(nfout, bfft_l, lsize, iesize)
    end if
! ==============================================================================

    call tstatc0_end(id_sname)
  end subroutine map_WFG_on_FFT_box_3D

  subroutine map_Rot_WFG_on_FFT_box_3D(ik,ib1,ib2,ibsize,lsize,bfft_l,wfr,wfi,nop,jtrs)
    integer, intent(in) :: ik, ib1,ib2, ibsize,lsize
    integer, intent(in) :: nop, jtrs
    real(kind=DP), intent(inout), dimension(lsize*kimg,ibsize) :: bfft_l
    real(kind=DP), intent(in), dimension(maxval(np_g1k)) :: wfr, wfi
    integer :: i, j, k, ii, jj, iesize, iadd, ierr
    integer, dimension(0:nrank_g-1) ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)::sta_r, sta_s
    real(kind=DP), allocatable, dimension(:,:),save :: sendbuf, recvbuf
    integer :: icnt_send, icnt_recv, lrank
    integer, parameter :: itag = 21
    integer, save :: savesize = 0, savesend=0, saverecv=0
    integer :: nx, ny, nz, nxp, nn, ix, iy, iz
    integer,save  :: id_sname = -1
#ifdef USE_ALLTOALLV
    integer, allocatable, dimension(:) :: sdsp, rdsp
#endif
    call tstatc0_begin('map_Rot_WFG_on_FFT_box_3D ',id_sname)

    iesize = ib2 - ib1 + 1

    if(allocated(sendbuf)) deallocate(sendbuf)
    allocate(sendbuf(wf_fft_maxsend_exx(ik,nop,jtrs)*kimg*iesize,0:nrank_g-1))
    if(allocated(recvbuf)) deallocate(recvbuf)
    allocate(recvbuf(wf_fft_maxrecv_exx(ik,nop,jtrs)*kimg*iesize,0:nrank_g-1))

    savesize = iesize
    sendbuf = 0.0d0
    recvbuf = 0.0d0

    if(k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                sendbuf(iesize*(wf_fft_index_exx(ii*2-1,ik,nop,jtrs)-1)+1, &
                        wf_fft_dist_exx(ii*2-1,ik,nop,jtrs)) = wfr(ii)
                sendbuf(iesize*(wf_fft_index_exx(ii*2  ,ik,nop,jtrs)-1)+1, &
                        wf_fft_dist_exx(ii*2  ,ik,nop,jtrs)) = wfr(ii)
             enddo
          enddo
       else
!OCL NORECURRENCE
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                iadd = iesize*2*(wf_fft_index_exx(ii*2-1,ik,nop,jtrs)-1)+jj*2
                sendbuf(iadd-1,wf_fft_dist_exx(ii*2-1,ik,nop,jtrs)) = wfr(ii)
                sendbuf(iadd,  wf_fft_dist_exx(ii*2-1,ik,nop,jtrs)) = wfi(ii)
                iadd = iesize*2*(wf_fft_index_exx(ii*2,ik,nop,jtrs)-1)+jj*2
                sendbuf(iadd-1,wf_fft_dist_exx(ii*2  ,ik,nop,jtrs)) =  wfr(ii)
                sendbuf(iadd,  wf_fft_dist_exx(ii*2  ,ik,nop,jtrs)) = -wfi(ii)
             enddo
          enddo
       endif
    else
       if(kimg == 1) then
!OCL NORECURRENCE
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                iadd = iesize*(wf_fft_index_exx(ii,ik,nop,jtrs)-1)+jj
                sendbuf(iadd,wf_fft_dist_exx(ii,ik,nop,jtrs)) = wfr(ii)
             enddo
          enddo
       else
!OCL NORECURRENCE
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                iadd = iesize*2*(wf_fft_index_exx(ii,ik,nop,jtrs)-1)+jj*2
                sendbuf(iadd-1,wf_fft_dist_exx(ii,ik,nop,jtrs)) = wfr(ii)
                sendbuf(iadd  ,wf_fft_dist_exx(ii,ik,nop,jtrs)) = wfi(ii)
             enddo
          enddo
       endif
    endif

#ifndef USE_ALLTOALLV
    icnt_recv = 0
    do i = 0, nrank_g - 1
       if(wf_fft_rcnt_exx(i,ik,nop,jtrs) /= 0) then
          call mpi_irecv(recvbuf(1,i), wf_fft_rcnt_exx(i,ik,nop,jtrs)*kimg*iesize, mpi_double_precision, &
      &                  i, itag, mpi_ke_world, req_r(icnt_recv), ierr)
          if(ierr /= 0) then
             write(nfout,*)' map_Rot_WFG_on_FFT_box_3D :  mpi_irecv error'
             call flush(nfout)
             call mpi_abort(mpi_comm_world,166,ierr)
          endif
          icnt_recv = icnt_recv + 1
       endif
    enddo
    icnt_send = 0
    do i = 0, nrank_g - 1
       if(wf_fft_scnt_exx(i,ik,nop,jtrs) /= 0) then
          call mpi_isend(sendbuf(1,i), wf_fft_scnt_exx(i,ik,nop,jtrs)*kimg*iesize, mpi_double_precision, &
      &                  i, itag, mpi_ke_world, req_s(icnt_send), ierr)
          if(ierr /= 0) then
             write(nfout,*)' map_Rot_WFG_on_FFT_box_3D :  mpi_isend error'
             call flush(nfout)
             call mpi_abort(mpi_comm_world,167,ierr)
          endif
          icnt_send = icnt_send + 1
       endif
    enddo
    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
    if(ierr /= 0) then
       write(nfout,*)' map_Rot_WFG_on_FFT_box_3D :  mpi_waitall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world,168,ierr)
    endif
    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
    if(ierr /= 0) then
       write(nfout,*)' map_Rot_WFG_on_FFT_box_3D :  mpi_waitall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world,169,ierr)
    endif
#else
    allocate(sdsp(0:nrank_g-1), stat=ierr)
    allocate(rdsp(0:nrank_g-1), stat=ierr)
    do i = 0, nrank_g - 1
       sdsp(i)=wf_fft_maxsend_exx(ik,nop,jtrs)*kimg*iesize*i
       rdsp(i)=wf_fft_maxrecv_exx(ik,nop,jtrs)*kimg*iesize*i
    enddo
    call MPI_ALLTOALLV(      sendbuf, wf_fft_scnt_exx(:,ik,nop,jtrs)*kimg*iesize, sdsp,&
   &   mpi_double_precision, recvbuf, wf_fft_rcnt_exx(:,ik,nop,jtrs)*kimg*iesize, rdsp,&
   &   mpi_double_precision, mpi_ke_world, ierr)
    if(ierr /= 0) then
       write(nfout,*)' m_ES_Vlocal_in_Rspace_3D : mpi_alltoallv error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 170, ierr)
    endif
    deallocate(sdsp)
    deallocate(rdsp)
#endif

    bfft_l = 0.0d0
#ifdef FFT_USE_SSL2
    nx = xyz_fft_x(2,1) - xyz_fft_x(1,1) + 1
    ny = xyz_fft_x(2,2) - xyz_fft_x(1,2) + 1
    nz = xyz_fft_x(2,3) - xyz_fft_x(1,3) + 1
    nxp = nx

    if(kimg == 1) then
!OCL NORECURRENCE
!OCL NOFLTLD
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt_exx(i,ik,nop,jtrs)
             iadd = wf_fft_recv_exx(j,ik,i,nop,jtrs)
             iz = (iadd-1)/(nx*ny)+1
             nn = mod(iadd,nx*ny)
             if(nn==0) then
                iy = ny
             else
                iy = (nn-1)/nx+1
             end if
             ix = mod(nn,nx)
             if(ix==0) ix = nx
             do k = 1, iesize
                bfft_l(ix+(iy-1)*nx+(iz-1)*nx*ny,k) = recvbuf(iesize*(j-1)+k,i)
             enddo
          enddo
       enddo
    else
!OCL NORECURRENCE
!OCL NOFLTLD
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt_exx(i,ik,nop,jtrs)
             iadd = wf_fft_recv_exx(j,ik,i,nop,jtrs)
             iz = (iadd-1)/(nx*ny)+1
             nn = mod(iadd,nx*ny)
             if(nn==0) then
                iy = ny
             else
                iy = (nn-1)/nx+1
             end if
             ix = mod(nn,nx)
             if(ix==0) ix = nx
             do k = 1, iesize
                bfft_l((ix+(iy-1)*nxp+(iz-1)*nxp*ny)*2-1,k) = recvbuf(iesize*2*(j-1)+k*2-1,i)
                bfft_l((ix+(iy-1)*nxp+(iz-1)*nxp*ny)*2  ,k) = recvbuf(iesize*2*(j-1)+k*2  ,i)
             enddo
          enddo
       enddo
    endif
#else
    if(kimg == 1) then
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt_exx(i,ik,nop,jtrs)
             do k = 1, iesize
                bfft_l(wf_fft_recv_exx(j,ik,i,nop,jtrs),k) = recvbuf(iesize*(j-1)+k,i)
             enddo
          enddo
       enddo
    else
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt_exx(i,ik,nop,jtrs)
             do k = 1, iesize
                bfft_l(wf_fft_recv_exx(j,ik,i,nop,jtrs)*2-1,k) = recvbuf(iesize*2*(j-1)+k*2-1,i)
                bfft_l(wf_fft_recv_exx(j,ik,i,nop,jtrs)*2  ,k) = recvbuf(iesize*2*(j-1)+k*2  ,i)
             enddo
          enddo
       enddo
    endif
#endif

! === FFT xzy ==================================================================
    if(sw_fft_xzy > 0) then
! ==============================================================================
    call m_FFT_Inverse_3D(nfout, bfft_l, lsize, iesize)
! === FFT xzy ==================================================================
    else
    call m_FFT_Inverse_XYZ_3D(nfout, bfft_l, lsize, iesize)
    end if
! ==============================================================================

    call tstatc0_end(id_sname)
  end subroutine map_Rot_WFG_on_FFT_box_3D

  subroutine map_FFT_box_on_RHOG_3D(ik,lsize, ibesize, wk_bfft_l, isrsize, fftsize, rhogr, rhogi)
    integer, intent(in)  :: ik, lsize, ibesize, isrsize, fftsize
    real(kind=DP), dimension(lsize*kimg,ibesize), intent(in)  :: wk_bfft_l
    real(kind=DP), dimension(ista_kngp:iend_kngp), intent(out) :: rhogr, rhogi
    integer, dimension(0:nrank_g-1)                       ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)       ::sta_r, sta_s
    integer, parameter :: itag = 11

    real(kind=DP), allocatable, dimension(:,:) :: sendbuf, recvbuf
    integer :: icnt_send, icnt_recv, ierr, lrank, i, j, k, iadd
#ifdef USE_ALLTOALLV
    integer, allocatable, dimension(:) :: sdsp, rdsp
#endif

    if(fft_wf_maxsend_rhog /= 0) then
       allocate(sendbuf(fft_wf_maxsend_rhog*kimg*ibesize,0:nrank_g-1), stat=ierr)
       sendbuf = 0.0d0
    else
       allocate(sendbuf(1,1), stat=ierr)
    endif
    if(fft_wf_maxrecv_rhog /= 0) then
       allocate(recvbuf(fft_wf_maxrecv_rhog*kimg*ibesize,0:nrank_g-1), stat=ierr)
       recvbuf = 0.0d0
    else
       allocate(recvbuf(1,1), stat=ierr)
    endif
    if(ierr /= 0) then
       write(nfout,*)' map_FFT_box_on_RHOG_3D :  Not allocate '
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 211, ierr)
    endif

#ifndef USE_ALLTOALLV
    if(fft_wf_maxrecv_rhog /= 0) then
       icnt_recv = 0
       lrank = mod(myrank_g,nrank_g)
       do i = 0, nrank_g - 1
          lrank = lrank + 1
          if(lrank > nrank_g -1) lrank = 0
          if(fft_wf_rcnt_rhog(lrank) /= 0) then
             call mpi_irecv(recvbuf(1,lrank), fft_wf_rcnt_rhog(lrank)*kimg*ibesize, &
            &               mpi_double_precision, lrank, itag, mpi_ke_world, req_r(icnt_recv), ierr)
              if(ierr /= 0) then
                 write(nfout,*)' map_FFT_box_on_RHOG_3D :  mpi_irecv error'
                 call flush(nfout)
                 call mpi_abort(mpi_comm_world, 212, ierr)
              endif
             icnt_recv = icnt_recv + 1
          endif
       enddo
    endif

    if(fft_wf_maxsend_rhog /= 0) then
#endif
       if(kimg == 1) then
!OCL NORECURRENCE
          do k = 1, nel_fft_x(myrank_g)
             if(fft_wf_index_rhog(k) == 0) cycle
             do i = 1, ibesize
                sendbuf(ibesize*(fft_wf_index_rhog(k)-1)+i,fft_wf_dist_rhog(k)) = wk_bfft_l(k,i)
             enddo
          end do
       else
!OCL NORECURRENCE
          do k = 1, nel_fft_x(myrank_g)
             if(fft_wf_index_rhog(k) == 0) cycle
             do i = 1, ibesize
                iadd = ibesize*2*(fft_wf_index_rhog(k)-1)+i*2
                sendbuf(iadd-1,fft_wf_dist_rhog(k)) = wk_bfft_l(k*2-1,i)
                sendbuf(iadd,  fft_wf_dist_rhog(k)) = wk_bfft_l(k*2  ,i)
             enddo
          end do
       end if

#ifndef USE_ALLTOALLV
       icnt_send = 0
       lrank = mod((myrank_g+1),nrank_g)
       do i = 0, nrank_g - 1
          lrank = lrank + 1
          if(lrank > (nrank_g - 1)) lrank = 0
          if(fft_wf_scnt_rhog(lrank) /= 0) then
             call mpi_isend(sendbuf(1,lrank), fft_wf_scnt_rhog(lrank)*kimg*ibesize, &
            &               mpi_double_precision, lrank, itag, mpi_ke_world, req_s(icnt_send), ierr)
              if(ierr /= 0) then
                 write(nfout,*)' map_FFT_box_on_RHOG_3D :  mpi_isend error'
                 call flush(nfout)
                 call mpi_abort(mpi_comm_world, 213, ierr)
              endif
             icnt_send = icnt_send + 1
          endif
       enddo
    endif

    if(fft_wf_maxrecv_rhog /= 0) then
       call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
       if(ierr /= 0) then
          write(nfout,*)' map_FFT_box_on_RHOG_3D :  mpi_waitall error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 214, ierr)
       endif
    endif

    if(fft_wf_maxsend_rhog /= 0) then
       call mpi_waitall(icnt_send, req_s, sta_s, ierr)
       if(ierr /= 0) then
          write(nfout,*)' map_FFT_box_on_RHOG_3D :  mpi_waitall error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 215, ierr)
       endif
    endif
#else
    allocate(sdsp(0:nrank_g-1), stat=ierr)
    allocate(rdsp(0:nrank_g-1), stat=ierr)
    do i = 0, nrank_g - 1
       sdsp(i)=fft_wf_maxsend_rhog*kimg*ibesize*i
       rdsp(i)=fft_wf_maxrecv_rhog*kimg*ibesize*i
    enddo
    call MPI_ALLTOALLV(      sendbuf, fft_wf_scnt_rhog(:)*kimg*ibesize, sdsp, &
   &   mpi_double_precision, recvbuf, fft_wf_rcnt_rhog(:)*kimg*ibesize, rdsp, &
   &   mpi_double_precision, mpi_ke_world, ierr)
    if(ierr /= 0) then
       write(nfout,*)' map_FFT_box_on_RHOG_3D :  mpi_alltoallv error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 216, ierr)
    endif
    deallocate(sdsp)
    deallocate(rdsp)
#endif

    if(kimg == 1) then
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          if(fft_wf_rcnt_rhog(i) /= 0) then
             do k = 1, fft_wf_rcnt_rhog(i)
                do j = 1, ibesize
! === DEBUG by tkato 2014/04/14 ================================================
!                  rhogr(fft_wf_recv_rhog(k,i)) = recvbuf(ibesize*(k-1)+j,i)
                   rhogr(fft_wf_recv_rhog(k,i)+ista_kngp-1) = recvbuf(ibesize*(k-1)+j,i)
! ==============================================================================
                enddo
             end do
          end if
       end do
    else
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          if(fft_wf_rcnt_rhog(i) /= 0) then
             do k = 1, fft_wf_rcnt_rhog(i)
                do j = 1, ibesize
                   iadd = ibesize*2*(k-1)+j*2
! === DEBUG by tkato 2014/04/14 ================================================
!                  rhogr(fft_wf_recv_rhog(k,i)) = recvbuf(iadd-1,i)
!                  rhogi(fft_wf_recv_rhog(k,i)) = recvbuf(iadd,  i)
                   rhogr(fft_wf_recv_rhog(k,i)+ista_kngp-1) = recvbuf(iadd-1,i)
                   rhogi(fft_wf_recv_rhog(k,i)+ista_kngp-1) = recvbuf(iadd,  i)
! ==============================================================================
                enddo
             end do
          end if
       end do
    end if

    if(allocated(sendbuf)) deallocate(sendbuf)
    if(allocated(recvbuf)) deallocate(recvbuf)
  end subroutine map_FFT_box_on_RHOG_3D

  subroutine map_RHOG_on_FFT_box_3D(ik,ib1,ib2,ibsize,lsize,bfft_l,rhogr,rhogi)
    integer, intent(in)                           :: ik, ib1,ib2, ibsize,lsize
    real(kind=DP), intent(inout), dimension(lsize*kimg,ibsize) :: bfft_l
    real(kind=DP), intent(in), dimension(ista_kngp:iend_kngp) :: rhogr, rhogi
    integer :: i, j, k, ii, jj, iesize, iadd, ierr
    integer, dimension(0:nrank_g-1) ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)::sta_r, sta_s
    real(kind=DP), allocatable, dimension(:,:),save :: sendbuf, recvbuf
    integer :: icnt_send, icnt_recv, lrank
    integer, parameter :: itag = 21
    integer, save :: savesize = 0, savesend=0, saverecv=0
    integer :: nx, ny, nz, nxp, nn, ix, iy, iz
    integer,save  :: id_sname = -1
#ifdef USE_ALLTOALLV
    integer, allocatable, dimension(:) :: sdsp, rdsp
#endif
    call tstatc0_begin('map_RHOG_on_FFT_box_3D ',id_sname)

    iesize = ib2 - ib1 + 1

    if(allocated(sendbuf)) deallocate(sendbuf)
    allocate(sendbuf(wf_fft_maxsend_rhog*kimg*iesize,0:nrank_g-1))
    if(allocated(recvbuf)) deallocate(recvbuf)
    allocate(recvbuf(wf_fft_maxrecv_rhog*kimg*iesize,0:nrank_g-1))
    savesize = iesize
    sendbuf = 0.0d0
    recvbuf = 0.0d0

    if(kimg == 1) then
!OCL NORECURRENCE
       do jj = 1, iesize
          do ii = 1, np_kngp
             if(ii+ista_kngp-1 > kg) cycle
             iadd = iesize*(wf_fft_index_rhog(ii)-1)+jj
             sendbuf(iadd,wf_fft_dist_rhog(ii)) = rhogr(ii+ista_kngp-1)
          enddo
       enddo
    else
!OCL NORECURRENCE
       do jj = 1, iesize
          do ii = 1, np_kngp
             if(ii+ista_kngp-1 > kg) cycle
             iadd = iesize*2*(wf_fft_index_rhog(ii)-1)+jj*2
             sendbuf(iadd-1,wf_fft_dist_rhog(ii)) = rhogr(ii+ista_kngp-1)
             sendbuf(iadd  ,wf_fft_dist_rhog(ii)) = rhogi(ii+ista_kngp-1)
          enddo
       enddo
    endif

#ifndef USE_ALLTOALLV
    icnt_recv = 0
    do i = 0, nrank_g - 1
       if(wf_fft_rcnt_rhog(i) /= 0) then
          call mpi_irecv(recvbuf(1,i), wf_fft_rcnt_rhog(i)*kimg*iesize, mpi_double_precision, &
         &               i, itag, mpi_ke_world, req_r(icnt_recv), ierr)
          if(ierr /= 0) then
             write(nfout,*)' map_RHOG_on_FFT_box_3D :  mpi_irecv error'
             call flush(nfout)
             call mpi_abort(mpi_comm_world,166,ierr)
          endif
          icnt_recv = icnt_recv + 1
       endif
    enddo
    icnt_send = 0
    do i = 0, nrank_g - 1
       if(wf_fft_scnt_rhog(i) /= 0) then
          call mpi_isend(sendbuf(1,i), wf_fft_scnt_rhog(i)*kimg*iesize, mpi_double_precision, &
         &               i, itag, mpi_ke_world, req_s(icnt_send), ierr)
          if(ierr /= 0) then
             write(nfout,*)' map_RHOG_on_FFT_box_3D :  mpi_isend error'
             call flush(nfout)
             call mpi_abort(mpi_comm_world,167,ierr)
          endif
          icnt_send = icnt_send + 1
       endif
    enddo
    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
    if(ierr /= 0) then
       write(nfout,*)' map_RHOG_on_FFT_box_3D :  mpi_waitall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world,168,ierr)
    endif
    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
    if(ierr /= 0) then
       write(nfout,*)' map_RHOG_on_FFT_box_3D :  mpi_waitall error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world,169,ierr)
    endif
#else
    allocate(sdsp(0:nrank_g-1), stat=ierr)
    allocate(rdsp(0:nrank_g-1), stat=ierr)
    do i = 0, nrank_g - 1
       sdsp(i)=wf_fft_maxsend_rhog*kimg*iesize*i
       rdsp(i)=wf_fft_maxrecv_rhog*kimg*iesize*i
    enddo
    call MPI_ALLTOALLV(      sendbuf, wf_fft_scnt_rhog(:)*kimg*iesize, sdsp,&
   &   mpi_double_precision, recvbuf, wf_fft_rcnt_rhog(:)*kimg*iesize, rdsp,&
   &   mpi_double_precision, mpi_ke_world, ierr)
    if(ierr /= 0) then
       write(nfout,*)' map_RHOG_on_FFT_box_3D : mpi_alltoallv error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 170, ierr)
    endif
    deallocate(sdsp)
    deallocate(rdsp)
#endif

    bfft_l = 0.0d0
#ifdef FFT_USE_SSL2
    nx = xyz_fft_x(2,1) - xyz_fft_x(1,1) + 1
    ny = xyz_fft_x(2,2) - xyz_fft_x(1,2) + 1
    nz = xyz_fft_x(2,3) - xyz_fft_x(1,3) + 1
    nxp = nx

    if(kimg == 1) then
!OCL NORECURRENCE
!OCL NOFLTLD
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt_rhog(i)
             iadd = wf_fft_recv_rhog(j,i)
             iz = (iadd-1)/(nx*ny)+1
             nn = mod(iadd,nx*ny)
             if(nn==0) then
                iy = ny
             else
                iy = (nn-1)/nx+1
             end if
             ix = mod(nn,nx)
             if(ix==0) ix = nx
             do k = 1, iesize
                bfft_l(ix+(iy-1)*nx+(iz-1)*nx*ny,k) = recvbuf(iesize*(j-1)+k,i)
             enddo
          enddo
       enddo
    else
!OCL NORECURRENCE
!OCL NOFLTLD
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt_rhog(i)
             iadd = wf_fft_recv_rhog(j,i)
             iz = (iadd-1)/(nx*ny)+1
             nn = mod(iadd,nx*ny)
             if(nn==0) then
                iy = ny
             else
                iy = (nn-1)/nx+1
             end if
             ix = mod(nn,nx)
             if(ix==0) ix = nx
             do k = 1, iesize
                bfft_l((ix+(iy-1)*nxp+(iz-1)*nxp*ny)*2-1,k) = recvbuf(iesize*2*(j-1)+k*2-1,i)
                bfft_l((ix+(iy-1)*nxp+(iz-1)*nxp*ny)*2  ,k) = recvbuf(iesize*2*(j-1)+k*2  ,i)
             enddo
          enddo
       enddo
    endif
#else
    if(kimg == 1) then
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt_rhog(i)
             do k = 1, iesize
                bfft_l(wf_fft_recv_rhog(j,i),k) = recvbuf(iesize*(j-1)+k,i)
             enddo
          enddo
       enddo
    else
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt_rhog(i)
             do k = 1, iesize
                bfft_l(wf_fft_recv_rhog(j,i)*2-1,k) = recvbuf(iesize*2*(j-1)+k*2-1,i)
                bfft_l(wf_fft_recv_rhog(j,i)*2  ,k) = recvbuf(iesize*2*(j-1)+k*2  ,i)
             enddo
          enddo
       enddo
    endif
#endif

! === FFT xzy ==================================================================
    if(sw_fft_xzy > 0) then
! ==============================================================================
    call m_FFT_Inverse_3D(nfout, bfft_l, lsize, iesize)
! === FFT xzy ==================================================================
    else
    call m_FFT_Inverse_XYZ_3D(nfout, bfft_l, lsize, iesize)
    end if
! ==============================================================================

    call tstatc0_end(id_sname)
  end subroutine map_RHOG_on_FFT_box_3D

  subroutine map_FFT_box_on_WFG_3D(ik,lsize, ibesize, wk_bfft_l, isrsize, fftsize, vxwr, vxwi)
    use m_Parallelization,     only : fft_wf_scnt, fft_wf_rcnt &
   &                                , fft_wf_recv &
   &                                , fft_wf_index, fft_wf_dist &
   &                                , fft_wf_maxrecv, fft_wf_maxsend
    integer, intent(in)  :: ik, lsize, ibesize, isrsize, fftsize
    real(kind=DP), dimension(lsize*kimg,ibesize), intent(in)  :: wk_bfft_l
    real(kind=DP), dimension(maxval(np_g1k)), intent(out) :: vxwr, vxwi
    integer, dimension(0:nrank_g-1)                       ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)       ::sta_r, sta_s
    integer, parameter :: itag = 11

    real(kind=DP), allocatable, dimension(:,:) :: sendbuf, recvbuf
    integer :: icnt_send, icnt_recv, ierr, lrank, i, j, k, iadd
#ifdef USE_ALLTOALLV
    integer, allocatable, dimension(:) :: sdsp, rdsp
#endif

    if(fft_wf_maxsend(ik) /= 0) then
       allocate(sendbuf(fft_wf_maxsend(ik)*kimg*ibesize,0:nrank_g-1), stat=ierr)
       sendbuf = 0.0d0
    else
       allocate(sendbuf(1,1), stat=ierr)
    endif
    if(fft_wf_maxrecv(ik) /= 0) then
       allocate(recvbuf(fft_wf_maxrecv(ik)*kimg*ibesize,0:nrank_g-1), stat=ierr)
       recvbuf = 0.0d0
    else
       allocate(recvbuf(1,1), stat=ierr)
    endif
    if(ierr /= 0) then
       write(nfout,*)' map_FFT_box_on_WFG_3D :  Not allocate '
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 211, ierr)
    endif

#ifndef USE_ALLTOALLV
    if(fft_wf_maxrecv(ik) /= 0) then
       icnt_recv = 0
       lrank = mod(myrank_g,nrank_g)
       do i = 0, nrank_g - 1
          lrank = lrank + 1
          if(lrank > nrank_g -1) lrank = 0
          if(fft_wf_rcnt(lrank,ik) /= 0) then
             call mpi_irecv(recvbuf(1,lrank), fft_wf_rcnt(lrank,ik)*kimg*ibesize, &
            &               mpi_double_precision, lrank, itag, mpi_ke_world, req_r(icnt_recv), ierr)
             if(ierr /= 0) then
                write(nfout,*)' map_FFT_box_on_WFG_3D :  mpi_irecv error'
                call flush(nfout)
                call mpi_abort(mpi_comm_world, 212, ierr)
             endif
             icnt_recv = icnt_recv + 1
          endif
       enddo
    endif

    if(fft_wf_maxsend(ik) /= 0) then
#endif
       if(kimg == 1) then
!OCL NORECURRENCE
          do k = 1, nel_fft_x(myrank_g)
             if(fft_wf_index(k,ik) == 0) cycle
             do i = 1, ibesize
                sendbuf(ibesize*(fft_wf_index(k,ik)-1)+i,fft_wf_dist(k,ik)) = wk_bfft_l(k,i)
             enddo
          end do
       else
!OCL NORECURRENCE
          do k = 1, nel_fft_x(myrank_g)
             if(fft_wf_index(k,ik) == 0) cycle
             do i = 1, ibesize
                iadd = ibesize*2*(fft_wf_index(k,ik)-1)+i*2
                sendbuf(iadd-1,fft_wf_dist(k,ik)) = wk_bfft_l(k*2-1,i)
                sendbuf(iadd,  fft_wf_dist(k,ik)) = wk_bfft_l(k*2  ,i)
             enddo
          end do
       end if

#ifndef USE_ALLTOALLV
       icnt_send = 0
       lrank = mod((myrank_g+1),nrank_g)
       do i = 0, nrank_g - 1
          lrank = lrank + 1
          if(lrank > (nrank_g - 1)) lrank = 0
          if(fft_wf_scnt(lrank,ik) /= 0) then
             call mpi_isend(sendbuf(1,lrank), fft_wf_scnt(lrank,ik)*kimg*ibesize, &
            &               mpi_double_precision, lrank, itag, mpi_ke_world, req_s(icnt_send), ierr)
             if(ierr /= 0) then
                write(nfout,*)' map_FFT_box_on_WFG_3D :  mpi_isend error'
                call flush(nfout)
                call mpi_abort(mpi_comm_world, 213, ierr)
             endif
             icnt_send = icnt_send + 1
          endif
       enddo
    endif

    if(fft_wf_maxrecv(ik) /= 0) then
       call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
       if(ierr /= 0) then
          write(nfout,*)' map_FFT_box_on_WFG_3D :  mpi_waitall error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 214, ierr)
       endif
    endif

    if(fft_wf_maxsend(ik) /= 0) then
       call mpi_waitall(icnt_send, req_s, sta_s, ierr)
       if(ierr /= 0) then
          write(nfout,*)' map_FFT_box_on_WFG_3D :  mpi_waitall error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 215, ierr)
       endif
    endif
#else
    allocate(sdsp(0:nrank_g-1), stat=ierr)
    allocate(rdsp(0:nrank_g-1), stat=ierr)
    do i = 0, nrank_g - 1
       sdsp(i)=fft_wf_maxsend(ik)*kimg*ibesize*i
       rdsp(i)=fft_wf_maxrecv(ik)*kimg*ibesize*i
    enddo
    call MPI_ALLTOALLV(      sendbuf, fft_wf_scnt(:,ik)*kimg*ibesize, sdsp, &
   &   mpi_double_precision, recvbuf, fft_wf_rcnt(:,ik)*kimg*ibesize, rdsp, &
   &   mpi_double_precision, mpi_ke_world, ierr)
    if(ierr /= 0) then
       write(nfout,*)' map_FFT_box_on_WFG_3D :  mpi_alltoallv error'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 216, ierr)
    endif
    deallocate(sdsp)
    deallocate(rdsp)
#endif

    if(kimg == 1) then
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          if(fft_wf_rcnt(i,ik) /= 0) then
             do k = 1, fft_wf_rcnt(i,ik)
                do j = 1, ibesize
                   vxwr(fft_wf_recv(k,ik,i)) = recvbuf(ibesize*(k-1)+j,i)
                enddo
             end do
          end if
       end do
    else
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          if(fft_wf_rcnt(i,ik) /= 0) then
             do k = 1, fft_wf_rcnt(i,ik)
                do j = 1, ibesize
                   iadd = ibesize*2*(k-1)+j*2
                   vxwr(fft_wf_recv(k,ik,i)) = recvbuf(iadd-1,i)
                   vxwi(fft_wf_recv(k,ik,i)) = recvbuf(iadd,  i)
                enddo
             end do
          end if
       end do
    end if

    if(allocated(sendbuf)) deallocate(sendbuf)
    if(allocated(recvbuf)) deallocate(recvbuf)
  end subroutine map_FFT_box_on_WFG_3D

! ==============================================================================
! === Subroutines for DEBUG use ================================================
! ==============================================================================
    subroutine make_map()
      use m_FFT, only : xyz_fft_y, xyz_fft_x
      integer :: id1, id2, nl, nm, nn, index_l, index_g
      integer :: i, j, k, ri, i_, j_, k_

      id1 = fft_box_size_WF(1,0)
      id2 = fft_box_size_WF(2,0)

      nl = xyz_fft_y(2,1) - xyz_fft_y(1,1) + 1
      nm = xyz_fft_y(2,2) - xyz_fft_y(1,2) + 1
      nn = xyz_fft_y(2,3) - xyz_fft_y(1,3) + 1

      if(kimg == 1) then
         do k = xyz_fft_y(1,3), xyz_fft_y(2,3)
            do j = xyz_fft_y(1,2), xyz_fft_y(2,2)
               do i = xyz_fft_y(1,1), xyz_fft_y(2,1), 2
                  do ri = 0, 1
                     i_ = i - xyz_fft_y(1,1) + 1
                     j_ = j - xyz_fft_y(1,2) + 1
                     k_ = k - xyz_fft_y(1,3) + 1
                     index_g = id1*id2*(k-1)+id1*(j-1)+i+ri
                     index_l = nl*nm*(k_-1)+2*nm*((i_-1)/2)+2*(j_-1)+1+ri
                     mapg2ly(index_g) = index_l
                  end do
               end do
            end do
         end do
      else
         do k = xyz_fft_y(1,3), xyz_fft_y(2,3)
            do j = xyz_fft_y(1,2), xyz_fft_y(2,2)
               do i = xyz_fft_y(1,1), xyz_fft_y(2,1)
                  do ri = 0, 1
                     i_ = i - xyz_fft_y(1,1) + 1
                     j_ = j - xyz_fft_y(1,2) + 1
                     k_ = k - xyz_fft_y(1,3) + 1
                     index_g = id1*id2*2*(k-1)+id1*2*(j-1)+2*(i-1)+ri+1
                     index_l = nl*nm*2*(k_-1)+nm*2*(i_-1)+2*(j_-1)+ri+1
                     mapg2ly(index_g) = index_l
                  end do
               end do
            end do
         end do
      end if

      nl = xyz_fft_x(2,1) - xyz_fft_x(1,1) + 1
      nm = xyz_fft_x(2,2) - xyz_fft_x(1,2) + 1
      nn = xyz_fft_x(2,3) - xyz_fft_x(1,3) + 1

      if(kimg == 1) then
         do k = xyz_fft_x(1,3), xyz_fft_x(2,3)
            do j = xyz_fft_x(1,2), xyz_fft_x(2,2)
               do i = xyz_fft_x(1,1), xyz_fft_x(2,1)
                  i_ = i - xyz_fft_x(1,1) + 1
                  j_ = j - xyz_fft_x(1,2) + 1
                  k_ = k - xyz_fft_x(1,3) + 1
                  index_g = id1*id2*(k-1)+id1*(j-1)+i
                  index_l = nl*nm*(k_-1)+nl*(j_-1)+i_
                  mapg2lx(index_g) = index_l
               end do
            end do
         end do
      else
         do k = xyz_fft_x(1,3), xyz_fft_x(2,3)
            do j = xyz_fft_x(1,2), xyz_fft_x(2,2)
               do i = xyz_fft_x(1,1), xyz_fft_x(2,1)
                  do ri = 0, 1
                     i_ = i - xyz_fft_x(1,1) + 1
                     j_ = j - xyz_fft_x(1,2) + 1
                     k_ = k - xyz_fft_x(1,3) + 1
                     index_g = id1*id2*2*(k-1)+id1*2*(j-1)+2*(i-1)+ri+1
                     index_l = nl*nm*2*(k_-1)+nl*2*(j_-1)+2*(i_-1)+ri+1
                     mapg2lx(index_g) = index_l
                  end do
               end do
            end do
         end do
      end if
    end subroutine make_map
end module m_ES_ExactExchange
