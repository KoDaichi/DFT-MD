
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

  use m_Realspace, only : nmesh_rs_aug_max,nmesh_rs_aug,meshx_rs_aug,meshy_rs_aug,meshz_rs_aug,meshxyz_rs_aug &
  & , qr_clm_ylm,dqr_clm_ylm,plmt1,plmt2,nlmtpair
  use m_ES_RSB, only : overlap,m_ES_RSB_unitary_transform_wf,m_ES_RSB_unitary_transform_vec,m_ES_RSB_alloc,m_ES_RSB_doit
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
#ifndef MEMORY_SAVE_EXX
  real(kind=DP), allocatable :: qitg_exx(:,:,:) ! d(kgp,nqitg,nqmk)
  real(kind=DP), allocatable :: ylm_exx(:,:,:) ! d(kgp,maxylm,nqmk)
#else
#ifndef MEMORY_SAVE_MORE_EXX
  real(kind=DP), allocatable :: qitg_exx(:,:,:) ! d(kgp,nqitg,nqmk)
#else
  real(kind=DP), allocatable :: qitg_exx(:,:) ! d(kgp,nqitg)
#endif
  real(kind=DP), allocatable :: ylm_exx(:,:) ! d(kgp,maxylm)
#endif
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


  real(kind=DP), allocatable, dimension(:,:) :: eko_hyb
  real(kind=DP), allocatable, dimension(:,:,:,:,:) :: fsrqm,fsiqm

  real(kind=DP), allocatable, dimension(:,:,:,:) :: exx_potential
#if defined(MEMORY_SAVE_EXX) && defined(MEMORY_SAVE_MORE_EXX)
  integer :: it, lcmax
  integer, allocatable, dimension(:) :: nmm_il3
  integer, allocatable, dimension(:,:) :: mm_il3
  real(kind=DP), allocatable, dimension(:,:) :: qrsps_mm
  real(kind=DP), allocatable, dimension(:) :: h, radr_kept, wos_kept
#endif
  integer, save :: id_sname_cdfft = -1

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


    call m_Files_open_nfscfzaj()
    allocate(wfv(isize,nval,kv3_for_innerloop,kimg));wfv=0.d0
    if(mype==0)then
      if(precision_WFfile == SP) then
        allocate(wfs(isize,kimg));wfs=0.d0
        do i=1,kv3_for_innerloop
           do j=1,negp
              read(nfscfzaj) wfs
              if(j.le.nval) wfv(:,j,i,:) = wfs(:,:)
           enddo
        enddo
        deallocate(wfs)
      else
        allocate(wfd(isize,kimg));wfd=0.d0
        do i=1,kv3_for_innerloop
           do j=1,negp
              read(nfscfzaj) wfd
              if(j.le.nval) wfv(:,j,i,:) = wfd(:,:)
           enddo
        enddo
        deallocate(wfd)
      endif
    endif
    call mpi_bcast(wfv,isize*nval*kv3_for_innerloop*kimg,mpi_double_precision,0,MPI_CommGroup,ierr)
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
          call phase_error_with_msg(nfout, "Did not find the index of k' in IBZ.",__LINE__,__FILE__)
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
          call phase_error_with_msg(nfout, "Did not find the index of k' in IBZ.",__LINE__,__FILE__)
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
       if(.not.allocated(exx_potential)) allocate(exx_potential(kg1,np_e,ista_k:iend_k,kimg));exx_potential = 0.d0
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
          call phase_error_with_msg(nfout, 'invalid fixed_charge_k_parallel',__LINE__,__FILE__)
        endif
    else
        call phase_error_with_msg(nfout, 'invalid icond',__LINE__,__FILE__)
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
             if(it1/=it2) call phase_error_with_msg(nfout, 'm_ES_EXX_crotylm: it1 and it2 differ.',__LINE__,__FILE__)
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
#ifndef MEMORY_SAVE_EXX
       if(sw_rspace_hyb==OFF) then
         deallocate(ylm_exx)
         deallocate(qitg_exx)
       endif
#else
#ifndef MEMORY_SAVE_MORE_EXX
       if(sw_rspace_hyb==OFF) deallocate(qitg_exx)
#endif
#endif
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


    logical :: trans
    integer, allocatable, dimension(:) :: ista
    integer,save  :: id_sname = -1

    if(sw_update_wfv==OFF) return
    if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) return

    allocate(ista(MPI_STATUS_SIZE))

    trans = .true.
    if(present(transform)) trans = transform
    call tstatc0_begin('m_ES_EXX_gather_valence_states ',id_sname,level=1)


    if(sw_rsb==ON.and.trans) then
       do ik=1,kv3,af+1
          if(map_k(ik) /= myrank_k) cycle
          call UnitaryTransform_WF(INVERSE,ik,zaj_l,.true.,.true.)
       enddo
    endif

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
       call mpi_allreduce(ibm,irev,1,MPI_INTEGER,MPI_MAX,MPI_CommGroup,ierr)
       nval = irev
    else
       nval = ibm
    end if


    if(allocated(wfv)) deallocate(wfv)
    if(allocated(occup_val)) deallocate(occup_val)
    allocate(wfv(kg1,nval,kv3,kimg))
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
             wfv(1:kg1,ib,ik,1:kimg) = zaj_l(1:kg1,map_z(ib1),ik,1:kimg)
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

    allocate(wfv_mpi(kg1,nval,kv3,kimg))
    allocate(occup_val_mpi(nval,kv3))
    call mpi_allreduce(wfv,wfv_mpi,kg1*nval*kv3*kimg,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_CommGroup,ierr)
    call mpi_allreduce(occup_val,occup_val_mpi,nval*kv3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_CommGroup,ierr)
    wfv = wfv_mpi
    occup_val = occup_val_mpi
    deallocate(wfv_mpi)
    deallocate(occup_val_mpi)
    if(modnrm == EXECUT) then
       allocate(fsr_mpi(nval,nlmta,kv3))
       allocate(fsi_mpi(nval,nlmta,kv3))
       call mpi_allreduce(fsr_exx,fsr_mpi,nval*nlmta*kv3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_CommGroup,ierr)
       call mpi_allreduce(fsi_exx,fsi_mpi,nval*nlmta*kv3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_CommGroup,ierr)
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

    if(sw_rsb==ON.and.trans) then
       do ik=1,kv3,af+1
          if(map_k(ik) /= myrank_k) cycle
          call UnitaryTransform_WF(DIRECT,ik,zaj_l,.true.,.true.)
       enddo
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

    if(.not.allocated(vc)) allocate(vc(nmax_G_hyb,nqmk))

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
                vc(1,ik) = vzero
             end if
             igs = kgs;ige=nmax_G_hyb
             do ig=igs,ige
                kg(1:3) = qmk(ik,1:3) + ngabc(ig,1:3)
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
                vc(1,ik) = vzero
             end if
             igs = kgs;ige=nmax_G_hyb
             do ig=igs,ige
                kg(1:3) = qmk(ik,1:3) + ngabc(ig,1:3)
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
         igs = kgs;ige=nmax_G_hyb
         do ig=igs,ige
            kg(1:3) = qmk(ik,1:3) + ngabc(ig,1:3)
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
         do ig=kgs,nmax_G_hyb
            kg(1:3) = vkbz(ikbz,1:3) -v0(1:3) +ngabc(ig,1:3)
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
         do ig=kgs,nmax_G_hyb
            kg(1:3) = vkbz(ikbz,1:3) -v0(1:3) +ngabc(ig,1:3)
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
       vxw(1:iba(ik),1:kimg) = exx_potential(1:iba(ik),ib,ik,1:kimg)
       if(present(exx))then
         exx=0.d0
         if(kimg==1)then
            do ig=1,iba(ik)
               exx = exx + zaj_l(ig,ib,ik,1)*vxw(ig,1)
            enddo
         else
            if(k_symmetry(ik) == GAMMA) then
               do ig=2,iba(ik)
                  exx = exx + zaj_l(ig,ib,ik,1) * vxw(ig,1) &
                      &     + zaj_l(ig,ib,ik,2) * vxw(ig,2)
               end do
               exx = exx * 2.d0 + zaj_l(1,ib,ik,1) * vxw(1,1)
            else
               do ig=1,iba(ik)
                  exx = exx + zaj_l(ig,ib,ik,1)*vxw(ig,1) &
                      &     + zaj_l(ig,ib,ik,2)*vxw(ig,2)
               end do
            endif
         endif
       endif
       return
    endif
    allocate(efsr_l(nlmta));efsr_l=0.d0
    allocate(efsi_l(nlmta));efsi_l=0.d0
    if(modnrm == EXECUT) call get_expkt_fs_b(ik,ib,fsr,fsi,efsr_l,efsi_l)
    allocate(zajbuf_r(kg1));zajbuf_r(:) = zaj_l(:,ib,ik,1)
    allocate(zajbuf_i(kg1));zajbuf_i(:) = zaj_l(:,ib,ik,kimg)
    kgw = kg1;kgv=kg1

    if(present(exx))then
       call apply_Vx_to_WF( ispin, ib, ik, kgw, kgv, zajbuf_r, zajbuf_i, efsr_l, efsi_l, vxw, ene, eo )
    else
       call apply_Vx_to_WF( ispin, ib, ik, kgw, kgv, zajbuf_r, zajbuf_i, efsr_l, efsi_l, vxw )
    endif

    deallocate(efsr_l)
    deallocate(efsi_l)
    deallocate(zajbuf_r)
    deallocate(zajbuf_i)
    if(present(exx)) exx = ene
    if(.not.eo.and.potential_update>0.and.store_p) &
    &  exx_potential(1:iba(ik),ib,ik,1:kimg) = vxw(1:iba(ik),1:kimg)
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
    ng=kg1
    allocate(vxw(ng,np_e,ista_k:iend_k,kimg))
       if(sw_rsb==ON.and.iup.ge.potential_update) call UnitaryTransform_WF(INVERSE,ik,zaj_l,.true.,.true.)
       do ib=1,np_e   ! MPI
          if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
             call m_ES_EXX_potential(nfout,ispin,ib,ik,ng,fsr_l(ib,:,ik),fsi_l(ib,:,ik),&
             & vxw(1:ng,ib,ik,1:kimg),iup,exx=exx,exx_only=.true.,store_potential=.false.)
          else
             call m_ES_EXX_potential(nfout,ispin,ib,ik,ng,fsr_l(ib,:,ik),fsr_l(ib,:,ik),&
             & vxw(1:ng,ib,ik,1:kimg),iup,exx=exx,exx_only=.true.,store_potential=.false.)
          endif
          eexx_kb(ib,ik) = occup_l(ib,ik)*exx
          eko_hyb(ib,ik) = exx
          if(u_eko) then
             eko_l(ib,ik) = eko_l(ib,ik) + exx
          endif
       end do
       if(sw_rsb==ON.and.iup.ge.potential_update) call UnitaryTransform_WF(DIRECT,ik,zaj_l,.true.,.true.)

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
    real(kind=DP), intent(out), dimension(kg1) :: vxdi

    integer :: ib,ig
    integer :: ii,m,jkbz,kk,jk,jop,jtrs
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
       do m=1,nval
          if(occup_val(m,kk) < DELTA) cycle
          occupation = wbz(jkbz)*occup_val(m,kk)/qwgt_for_innerloop(kk)/dble(kv3_for_innerloop)
          if(kimg==1) then
             do ii=1,iba_for_innerloop(kk)
                rhor(ii) = rhor(ii) + occupation * ( wfv(ii,m,kk,1) * wfv(ii,m,kk,1) )
             end do
          else
             do ii=1, min( iba_for_innerloop(kk), nmax_G_hyb )
                rhor(ii) = rhor(ii) + occupation * ( wfv(ii,m,kk,1) * wfv(ii,m,kk,1) &
                                                 & + wfv(ii,m,kk,2) * wfv(ii,m,kk,2) )
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
    do ii=1,iba(ik)
       ig = nbase(ii,ik)
       if ( ig <= nmax_G_hyb ) then
          vxdi(ii) = fac * rhor(ig)
       endif
    end do

    deallocate(rho)
    deallocate(pot)
    deallocate(vx)
    deallocate(rhor)
    deallocate(rhoi)

    call tstatc0_end(id_sname)
  end subroutine m_ES_EXX_Diagonal_part


  subroutine m_ES_EXX_add_Diagonal_part(ik,ibo,vxdi,vnldi)
    implicit none
    integer, intent(in) :: ik, ibo
    real(kind=DP), intent(in), dimension(kg1) :: vxdi
    real(kind=DP), intent(out), dimension(kg1) :: vnldi

    integer :: ib

    ib = map_z(ibo)
    vnldi = vnldi + vxdi
  end subroutine m_ES_EXX_add_Diagonal_part

  subroutine m_ES_Vexx_W(ik,iupdate,store_exxp)
    implicit none
    integer, intent(in) :: ik
    integer, intent(in), optional :: iupdate
    logical, intent(in), optional :: store_exxp

    logical :: store_e

    integer :: ib,ig,ispin,iup
    real(kind=DP), allocatable, dimension(:,:,:,:) :: vxw !d(kg1,kimg)
    real(kind=DP), allocatable, dimension(:,:,:) :: vexx
    real(kind=DP) :: exx
    integer :: ng
    integer,save  :: id_sname = -1
    iup = 2
    if(present(iupdate)) iup = iupdate
    store_e = .true.
    if(present(store_exxp)) store_e = store_exxp
    call tstatc0_begin('m_ES_Vexx_W ',id_sname,level=1)

    ispin = mod(ik-1,nspin)+1

    ng = kg1
       allocate(vxw(ng,np_e,ista_k:iend_k,kimg))
       if(sw_rsb==ON) call UnitaryTransform_WF(INVERSE,ik,zaj_l,.true.,.true.)
       do ib=1,np_e   ! MPI
          if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
             call m_ES_EXX_potential(nfout,ispin,ib,ik,ng,fsr_l(ib,:,ik),fsi_l(ib,:,ik),&
             &    vxw(1:ng,ib,ik,1:kimg),iup,store_potential=store_e)
          else
             call m_ES_EXX_potential(nfout,ispin,ib,ik,ng,fsr_l(ib,:,ik),fsr_l(ib,:,ik),&
             &    vxw(1:ng,ib,ik,1:kimg),iup,store_potential=store_e)
          endif
       end do
       if(sw_rsb==ON)then
          call UnitaryTransform_WF(DIRECT,ik,vxw,.false.,.false.)
          call UnitaryTransform_WF(DIRECT,ik,zaj_l,.true.,.true.)
          if(potential_update>0.and.store_e) &
          &  call UnitaryTransform_WF(DIRECT,ik,exx_potential,.false.,.false.)
       endif
       do ib=1,np_e
          call m_ES_Vexx_add_vexx(ik,ib,vxw(1:ng,ib,ik,1:kimg))
       enddo
       deallocate(vxw)
    call tstatc0_end(id_sname)
  end subroutine m_ES_Vexx_W

  subroutine m_ES_Vexx_add_vexx(ik,ib,vxw_exx)
    integer, intent(in) :: ib,ik
    real(kind=DP),dimension(kg1,kimg),intent(in) :: vxw_exx
    integer :: ig,iadd
    logical ma
    if(kimg==1) then
       do ig=1,iba(ik)
          vnlph_l(ig,ib,1) = vnlph_l(ig,ib,1) + vxw_exx(ig,1)
       end do
    else
       do ig=1,iba(ik)
          vnlph_l(ig,ib,1) = vnlph_l(ig,ib,1) + vxw_exx(ig,1)
          vnlph_l(ig,ib,2) = vnlph_l(ig,ib,2) + vxw_exx(ig,2)
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
    integer :: ikk
    integer,save  :: id_sname = -1,id_sname1=-2,id_sname2=-3,id_sname3=-4
    call tstatc0_begin('apply_Vx_to_WF ',id_sname)

    ikk = k_index(ik)

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
    allocate(wfsr(kg1_for_innerloop));wfsr=0.d0
    if(kimg==2) then
       allocate(wfsi(kg1_for_innerloop));wfsi=0.d0
    endif
    allocate(cosgt(kg1_for_innerloop))
    if(kimg==2) allocate(singt(kg1_for_innerloop))
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
       do m=1,nval
          if(occup_val(m,kk) < DELTA) cycle
          if(sw_rsb==ON)then
             if(.not.overlap(ib,m))cycle
          endif
          occupation = wbz(jkbz) * occup_val(m,kk)/qwgt_for_innerloop(kk)/dble(kv3_for_innerloop)
          if(kimg==1) then
             dnorm = 0.d0
             do ig=1,iba_for_innerloop(kk)
                wfsr(ig) = wfv(ig,m,kk,1) * cosgt(ig)
                dnorm = dnorm + wfsr(ig)*wfsr(ig)
             end do
             dnorm = 1.d0/sqrt(dnorm)
             wfsr = wfsr * dnorm
             call map_Rot_WFG_on_FFT_box(kk,jop,jtrs,wfsr,wfsr,wfm)
          else
             if(jtrs==0) then
                do ig=1,iba_for_innerloop(kk)
                   wfsr(ig) = wfv(ig,m,kk,1) * cosgt(ig) - wfv(ig,m,kk,2) * singt(ig)
                   wfsi(ig) = wfv(ig,m,kk,2) * cosgt(ig) + wfv(ig,m,kk,1) * singt(ig)
                end do
             else
                do ig=1,iba_for_innerloop(kk)
                   wfsr(ig) = wfv(ig,m,kk,1) * cosgt(ig) - wfv(ig,m,kk,2) * singt(ig)
                   wfsi(ig) = -wfv(ig,m,kk,2) * cosgt(ig) - wfv(ig,m,kk,1) * singt(ig)
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
                   call add_RHOG_hard_part(iqmk(jkbz,ikk),rhogr,rhogi,bdwr,bdwi,fsr,fsi)
                endif
             end if
             if(present(eexx))then
                exx = 0.d0
                call sum_rho_vc_rho(rhogr,rhogi,vc(1,iqmk(jkbz,ikk)),exx)
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
                      call integrate_QijVnm(iqmk(jkbz,ikk),rhogr,rhogi,fsr,fsi,qvr,qvi,dfsr,dfsi,dqvr,dqvi,gqvr,gqvi)
                   else
                      call integrate_QijVnm(iqmk(jkbz,ikk),rhogr,rhogi,fsr,fsi,qvr,qvi)
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
       call sum_EXX_force_terms(force_l,bdwr,bdwi,dbdwr,dbdwi,sumqvr,sumqvi,sumdqvr,sumdqvi)
    else
#ifdef FFTW3
       call m_FFT_exx(nfout,sumdel,DIRECT,OFF) ! sum delta(G)
#else
       call m_FFT_WF(ELECTRON,nfout,sumdel,DIRECT,OFF) ! sum delta(G)
#endif
       if(kimg==1) then
          call map_FFT_box_on_WFG(ik,vxw(1,1),vxw(1,1),sumdel)
       else
          call map_FFT_box_on_WFG(ik,vxw(1,1),vxw(1,2),sumdel)
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

    if(sw_rspace_hyb==ON.and.modnrm==EXECUT)then
       deallocate(rhor)
       deallocate(rhoi)
       deallocate(afft)
    endif

    call tstatc0_end(id_sname)
  end subroutine apply_Vx_to_WF


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
       call mpi_allreduce(eexx,eexx_mpi,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_CommGroup,ierr)
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


  subroutine sum_rho_vc_rho(rhor,rhoi,vc,exx)
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

       do i = 1, kgp
          ia = ngabc(i,1)
          ib = ngabc(i,2)
          ic = ngabc(i,3)
          g_list(ia,ib,ic) = i
       end do
    end if

    if(npes > 1) then
       allocate(ngpt_exx0_tmp(kgp))
       do iopr=1,nopr
          ngpt_exx0_tmp = 0
          do i = ista_kngp, iend_kngp
             ngpt_exx0_tmp(i) = ngpt_l(i,iopr)
          end do
          call mpi_allreduce(MPI_IN_PLACE,ngpt_exx0_tmp,kgp,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
          ngpt_exx(1:kg,iopr,0) = ngpt_exx0_tmp(1:kg)

          if(ntrs>0) then
             allocate(ngpt_exx_tmp(kg))
             ngpt_exx_tmp = 0
             do i=1,kg
                ii = ngpt_exx0_tmp(i)
                ia = -ngabc(ii,1)
                ib = -ngabc(ii,2)
                ic = -ngabc(ii,3)
                ngpt_exx_tmp(i) = g_list(ia,ib,ic)
             end do
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
                ia = -ngabc(ii,1)
                ib = -ngabc(ii,2)
                ic = -ngabc(ii,3)
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
#ifndef MEMORY_SAVE_EXX
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
    if(.not.allocated(ylm_exx)) allocate(ylm_exx(nmax_G_hyb,n,nqmk))

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
          kg(1:3) = qmk(ik,1:3) + ngabc(ig,1:3)
          g2          = ttr(1)*kg(1)*kg(1) &
          &           + ttr(2)*kg(2)*kg(2) &
          &           + ttr(3)*kg(3)*kg(3) &
          &           + ttr(4)*kg(1)*kg(2) &
          &           + ttr(5)*kg(2)*kg(3) &
          &           + ttr(6)*kg(3)*kg(1)
          gqmk(ig,1:3) = kg(1:3)
          gqmkr(ig) = sqrt(g2)
       end do
       if(npes > 1) then
          allocate(ylm_t(nmax_G_hyb))
          do i=1,n
             ylm_t = 0.d0
             if ( ista_kngp <= nmax_G_hyb ) then
                call m_pwBS_sphrp_exx( i, rltv, ista_kngp, iend_kngp0, &
                     &                 gqmk, gqmkr, ylm_t(ista_kngp:iend_kngp0) )
             endif
             call mpi_allreduce( ylm_t, ylm_exx(1,i,ik), nmax_G_hyb, &
                  &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
          end do
          deallocate(ylm_t)
       else
          do i=1,n
             call m_pwBS_sphrp_exx(i,rltv,1,nmax_G_hyb,gqmk,gqmkr,ylm_exx(1,i,ik))
          end do
       end if

    end do
    !!stop 'm_ES_EXX_ylm'
    if ( allocated(gqmk) )  deallocate(gqmk)
    if ( allocated(gqmkr) ) deallocate(gqmkr)

    call tstatc0_end(id_sname)
  end subroutine m_ES_EXX_ylm
#else
  subroutine m_ES_EXX_ylm_each_k(ik)
     implicit none
     integer, intent(in) :: ik
     real(kind=DP), allocatable, dimension(:,:) :: gqmk ! d(kgp,3)
     real(kind=DP), allocatable, dimension(:) :: gqmkr ! d(kgp)
     integer :: i, n, ig

     real(kind=DP) :: kg(3), ttr(6), g2
     integer,save  :: id_sname = -1
     call tstatc0_begin('m_ES_EXX_ylm ', id_sname, level=1)

     call getttr(rltv,ttr)

     call m_PP_find_maximum_l(n) ! n-1: maximum l
     n = (n-1) + (n-1) + 1
     n = n*n
     if(.not. allocated(ylm_exx)) allocate(ylm_exx(nmax_G_hyb,n))

     allocate(gqmk(kgp,3))
     allocate(gqmkr(kgp))
     do ig = 1, kgp
          kg(1:3) = qmk(ik,1:3) + ngabc(ig,1:3)
          g2          = ttr(1)*kg(1)*kg(1) &
          &           + ttr(2)*kg(2)*kg(2) &
          &           + ttr(3)*kg(3)*kg(3) &
          &           + ttr(4)*kg(1)*kg(2) &
          &           + ttr(5)*kg(2)*kg(3) &
          &           + ttr(6)*kg(3)*kg(1)
          gqmk(ig,1:3) = kg(1:3)
          gqmkr(ig) = sqrt(g2)
     end do
     do i = 1, n
        call m_pwBS_sphrp_exx(i,rltv,1,nmax_G_hyb,gqmk,gqmkr,ylm_exx(1,i))
     end do
     deallocate(gqmk)
     deallocate(gqmkr)

     call tstatc0_end(id_sname)
  end subroutine m_ES_EXX_ylm_each_k
#endif

#ifndef MEMORY_SAVE_MORE_EXX
  subroutine check_qitg()
    implicit none
    integer :: iq, ips , ipe

    ips = 1
    ipe = nmax_G_hyb

    do iq=1,nqitg
       write(nfout,'("iq=",i5,1x,"qitg_l=",f20.5,1x,"qitg_exx=",f20.5,1x)') iq, qitg_l(ips,iq), qitg_exx(ips,iq,1)/univol
    end do

    do iq=1,nqitg
       write(nfout,'("iq=",i5,1x,"diff=",f20.5)') iq, sum(qitg_l(ips:ipe,iq)-qitg_exx(ips:ipe,iq,1)/univol)
    end do

    call phase_error_with_msg(nfout, 'check_qitg',__LINE__,__FILE__)
  end subroutine check_qitg

  subroutine check_qitg_qmk()
    implicit none
    integer :: iq, ik, ips, ipe

    ips = 1
    ipe = nmax_G_hyb

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

    call phase_error_with_msg(nfout, 'check_qitg_qmk',__LINE__,__FILE__)
  end subroutine check_qitg_qmk
#endif

#ifndef MEMORY_SAVE_EXX
  subroutine check_ylm_exx()
    implicit none
    integer :: i, n, ik, ips, ipe

    call m_PP_find_maximum_l(n)    ! n-1: maximum l
    n = (n-1) + (n-1) + 1
    n = n*n

    ips = 1
    ipe = nmax_G_hyb

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

    call phase_error_with_msg(nfout, 'check_ylm_exx',__LINE__,__FILE__)
  end subroutine check_ylm_exx
#endif

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


  subroutine add_RHOG_hard_part(iqmk,rhogr,rhogi,fnr,fni,fmr,fmi)
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

    do ibl1=ips, ipe, ibsize
       rhogr_red=0.d0;  rhogi_red=0.d0

       ibl2=min(ipe,ibl1+ibsize-1)
       do iq=1,nqitg
          do i=1,ibl2-ibl1+1
#if defined(MEMORY_SAVE_EXX) && defined(MEMORY_SAVE_MORE_EXX)
             qitg_red(i,iq) = qitg_exx(i+ibl1-1,iq)
#else
             qitg_red(i,iq) = qitg_exx(i+ibl1-1,iq,iqmk)
#endif
          enddo
       enddo
       do inn=1,n*n
          do i=1,ibl2-ibl1+1
#ifndef MEMORY_SAVE_EXX
             ylm_red(i,inn) = ylm_exx(i+ibl1-1,inn,iqmk)
#else
             ylm_red(i,inn) = ylm_exx(i+ibl1-1,inn)
#endif
          enddo
       enddo

       do ia=1,natm
          it = ityp(ia)
          if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
          do i=1,ibl2-ibl1+1
             ph = PAI2 * dot_product(pos(ia,1:3),ngabc(i+ibl1-1,1:3)+qmk(iqmk,1:3))
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


! ================================== KT_Test ========================= 12.5Exp
  subroutine integrate_QijVnm(iqmk,potr,poti,fmr,fmi,qvr,qvi,dfmr,dfmi,dqvr,dqvi,gqvr,gqvi)
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

    do ibl1=ips, ipe, ibsize
       ibl2=min(ipe,ibl1+ibsize-1)
       do iq=1,nqitg
          do i=1,ibl2-ibl1+1
#if defined(MEMORY_SAVE_EXX) && defined(MEMORY_SAVE_MORE_EXX)
             qitg_red(i,iq) = qitg_exx(i+ibl1-1,iq)
#else
             qitg_red(i,iq) = qitg_exx(i+ibl1-1,iq,iqmk)
#endif
          enddo
       enddo
       do inn=1,n*n
          do i=1,ibl2-ibl1+1
#ifndef MEMORY_SAVE_EXX
             ylm_red(i,inn) = ylm_exx(i+ibl1-1,inn,iqmk)
#else
             ylm_red(i,inn) = ylm_exx(i+ibl1-1,inn)
#endif
          enddo
       enddo

       if(force_mode)then
          do i=1,ibl2-ibl1+1
             gvec(i,1:3) = ngabc(i+ibl1-1,1:3) + qmk(iqmk,1:3)
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
             ph = PAI2 * dot_product(pos(ia,1:3),ngabc(i+ibl1-1,1:3)+qmk(iqmk,1:3))
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


  subroutine add_Vx_hard_part(ik,vxw,sumqvr,sumqvi)
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
    integer,save  :: id_sname = -1
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
             do i=1,iba(ik)
                fac = iwei(ia)*snl(i,lmtt1,iksnl)
                vxw(i,1) = vxw(i,1) + fac * (cosgt(i)*sumqvr(ilmta1)+singt(i)*sumqvi(ilmta1))
                if(kimg==2) vxw(i,2) = vxw(i,2) + fac * (cosgt(i)*sumqvi(ilmta1)-singt(i)*sumqvr(ilmta1))
             end do
          else if(iil==1) then
             do i=1,iba(ik)
                fac = iwei(ia)*snl(i,lmtt1,iksnl)
                vxw(i,1) = vxw(i,1) + fac * (cosgt(i)*sumqvi(ilmta1)-singt(i)*sumqvr(ilmta1))
                if(kimg==2) vxw(i,2) = vxw(i,2) - fac * (cosgt(i)*sumqvr(ilmta1)+singt(i)*sumqvi(ilmta1))
             end do
          else if(iil==2) then
             do i=1,iba(ik)
                fac = -iwei(ia)*snl(i,lmtt1,iksnl)
                vxw(i,1) = vxw(i,1) + fac * (cosgt(i)*sumqvr(ilmta1)+singt(i)*sumqvi(ilmta1))
                if(kimg==2) vxw(i,2) = vxw(i,2) + fac * (cosgt(i)*sumqvi(ilmta1)-singt(i)*sumqvr(ilmta1))
             end do
          else if(iil==3) then
             do i=1,iba(ik)
                fac = -iwei(ia)*snl(i,lmtt1,iksnl)
                vxw(i,1) = vxw(i,1) + fac * (cosgt(i)*sumqvi(ilmta1)-singt(i)*sumqvr(ilmta1))
                if(kimg==2) vxw(i,2) = vxw(i,2) - fac * (cosgt(i)*sumqvr(ilmta1)+singt(i)*sumqvi(ilmta1))
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
    real(kind=DP), intent(in), dimension(nlmta) :: fsr_l,fsi_l
    real(kind=DP), intent(out), dimension(nlmta) :: efsr_l,efsi_l

    integer :: it, ia, lmt, ilmta
    real(kind=DP) :: ph, coskt, sinkt
    efsr_l(1:nlmta) = 0.d0
    efsi_l(1:nlmta) = 0.d0
    if(k_symmetry(ik)==GAMMA)then
       efsr_l(:) = fsr_l(:)
    else
       do ia=1,natm
          it = ityp(ia)
          if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
          ph = PAI2 * dot_product(pos(ia,1:3),vkxyz(ik,1:3,BUCS))
          coskt = cos(ph)
          sinkt = sin(ph)
          do lmt = 1, ilmt(it)
             ilmta = lmta(lmt,ia)
             efsr_l(ilmta) = coskt * fsr_l(ilmta) - sinkt * fsi_l(ilmta)
             efsi_l(ilmta) = sinkt * fsr_l(ilmta) + coskt * fsi_l(ilmta)
          end do
       end do
    endif

  end subroutine get_expkt_fs_b

  subroutine get_expkt_fs(ik,fsr_l,fsi_l,efsr_l,efsi_l)
    implicit none
    integer, intent(in) :: ik
    real(kind=DP), intent(in), dimension(np_e,nlmta,ista_k:iend_k) :: fsr_l,fsi_l
    real(kind=DP), intent(out), dimension(np_e,nlmta) :: efsr_l,efsi_l

    integer :: it, ia, lmt, ilmta
    real(kind=DP) :: ph, coskt, sinkt
    integer :: ib
    efsr_l = 0.d0
    efsi_l = 0.d0
    do ib=1,np_e
       call get_expkt_fs_b(ik,ib,fsr_l(ib,1:nlmta,ik),fsi_l(ib,1:nlmta,ik),efsr_l(ib,1:nlmta),efsi_l(ib,1:nlmta))
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
             allocate(zajbuf_r(kg1));zajbuf_r(:) = zaj_l(:,ib,ik,1)
             allocate(zajbuf_i(kg1));zajbuf_i(:) = zaj_l(:,ib,ik,kimg)
             kgw = kg1;kgv=kg1

               if(kimg==1) then
                  call apply_Vx_to_WF(ispin,ib,ik,kgw,kgv,zajbuf_r,zajbuf_r &
                      & ,efsr_l(ib,1:nlmta),efsr_l(ib,1:nlmta) &
                      & ,dbdwr=dfsr_l(ib,1:nlmta,ik,1:3),dbdwi=dfsr_l(ib,1:nlmta,ik,1:3),force_l=force_l)
               else
                  call apply_Vx_to_WF(ispin,ib,ik,kgw,kgv,zajbuf_r,zajbuf_i &
                      & ,efsr_l(ib,1:nlmta),efsi_l(ib,1:nlmta) &
                      & ,dbdwr=dfsr_l(ib,1:nlmta,ik,1:3),dbdwi=dfsi_l(ib,1:nlmta,ik,1:3),force_l=force_l)
               end if
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
       call mpi_allreduce(force,force_mpi,natm*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_CommGroup,ierr)
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
       call mpi_allreduce(dfsr_exx,dfsr_mpi,nval*nlmta*kv3*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_CommGroup,ierr)
       call mpi_allreduce(dfsi_exx,dfsi_mpi,nval*nlmta*kv3*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_CommGroup,ierr)
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
                   do i=1,iba(ik)
                      fac = snl(i,lmtt1,iksnl)*fac2
                      fr = fac * ( cosgt(i)*zaj_l(i,ib,ik,1)-singt(i)*zaj_l(i,ib,ik,2) )
                      fi = fac * ( cosgt(i)*zaj_l(i,ib,ik,2)+singt(i)*zaj_l(i,ib,ik,1) )
                      dfsr_l(ib,ilmta1,ik,1:3) = dfsr_l(ib,ilmta1,ik,1:3) + gvec(i,1:3) * fr
                      dfsi_l(ib,ilmta1,ik,1:3) = dfsi_l(ib,ilmta1,ik,1:3) + gvec(i,1:3) * fi
                   end do
                end do
             else if(iil==1) then
                do ib=1,np_e
                   if(occup_l(ib,ik) < DELTA) cycle
                   do i=1,iba(ik)
                      fac = snl(i,lmtt1,iksnl)*fac2
                      fr = -fac * ( cosgt(i)*zaj_l(i,ib,ik,2)+singt(i)*zaj_l(i,ib,ik,1) )
                      fi =  fac * ( cosgt(i)*zaj_l(i,ib,ik,1)-singt(i)*zaj_l(i,ib,ik,2) )
                      dfsr_l(ib,ilmta1,ik,1:3) = dfsr_l(ib,ilmta1,ik,1:3) + gvec(i,1:3) * fr
                      dfsi_l(ib,ilmta1,ik,1:3) = dfsi_l(ib,ilmta1,ik,1:3) + gvec(i,1:3) * fi
                   end do
                end do
             else if(iil==2) then
                do ib=1,np_e
                   if(occup_l(ib,ik) < DELTA) cycle
                   do i=1,iba(ik)
                      fac = snl(i,lmtt1,iksnl)*fac2
                      fr = -fac * ( cosgt(i)*zaj_l(i,ib,ik,1)-singt(i)*zaj_l(i,ib,ik,2) )
                      fi = -fac * ( cosgt(i)*zaj_l(i,ib,ik,2)+singt(i)*zaj_l(i,ib,ik,1) )
                      dfsr_l(ib,ilmta1,ik,1:3) = dfsr_l(ib,ilmta1,ik,1:3) + gvec(i,1:3) * fr
                      dfsi_l(ib,ilmta1,ik,1:3) = dfsi_l(ib,ilmta1,ik,1:3) + gvec(i,1:3) * fi
                   end do
                end do
             else if(iil==3) then
                do ib=1,np_e
                   if(occup_l(ib,ik) < DELTA) cycle
                   do i=1,iba(ik)
                      fac = snl(i,lmtt1,iksnl)*fac2
                      fr =  fac * ( cosgt(i)*zaj_l(i,ib,ik,2)+singt(i)*zaj_l(i,ib,ik,1) )
                      fi = -fac * ( cosgt(i)*zaj_l(i,ib,ik,1)-singt(i)*zaj_l(i,ib,ik,2) )
                      dfsr_l(ib,ilmta1,ik,1:3) = dfsr_l(ib,ilmta1,ik,1:3) + gvec(i,1:3) * fr
                      dfsi_l(ib,ilmta1,ik,1:3) = dfsi_l(ib,ilmta1,ik,1:3) + gvec(i,1:3) * fi
                   end do
                end do
             end if
          end do
       end do
    end do

    deallocate(cosgt)
    deallocate(singt)
    deallocate(gvec)

    call tstatc0_end(id_sname)
  end subroutine drv_betar_dot_WFs_exx
#if defined(MEMORY_SAVE_EXX) && defined(MEMORY_SAVE_MORE_EXX)
  subroutine qitgft_qmk_each_k_setup(it_in,nmm_il3_in,mm_il3_in,qrsps_mm_in,lcmax_in,h_in)
     use m_PseudoPotential, only : wos, nqitg_sp
     implicit none
     integer, intent(in) :: it_in
     integer, intent(in) :: nmm_il3_in(lcmax_in+1)
     integer, intent(in) :: mm_il3_in(nqitg_sp(it_in),lcmax_in+1)
     real(kind=DP), intent(in) :: qrsps_mm_in(mmesh,nqitg_sp(it_in))
     integer, intent(in) :: lcmax_in
     real(kind=DP), intent(in) :: h_in(it_in)

     it = it_in
     lcmax = lcmax_in
     allocate(nmm_il3(lcmax+1)); nmm_il3 = nmm_il3_in
     allocate(mm_il3(nqitg_sp(it),lcmax+1)); mm_il3 = mm_il3_in
     allocate(qrsps_mm(mmesh,nqitg_sp(it))); qrsps_mm = qrsps_mm_in
     allocate(h(it)); h = h_in
     allocate(radr_kept(mmesh)); radr_kept = radr
     allocate(wos_kept(mmesh)); wos_kept = wos

     return
  end subroutine qitgft_qmk_each_k_setup

  subroutine qitgft_qmk_each_k(ik)
     use m_PseudoPotential,    only : nqitg_sp
     implicit none
     integer, intent(in) :: ik

     integer :: mm0, i, il3, mm, mmp, n, idp, iq
     real(kind=DP) :: ttr(6), kg(3)
     real(kind=DP) :: qitg_sh, g2, gabs
     real(kind=DP) :: wkx(mmesh), wky(mmesh)
     real(kind=DP), allocatable :: gqmk_l(:) ! d(kgp)
     real(kind=DP), allocatable :: qitg_exx_l(:,:) ! d(kgp,nqitg_sp(it))
     integer,save  :: id_sname = -1

     if(m_PP_include_vanderbilt_pot(it) == SKIP .or. sw_rspace_hyb == ON) return

     call tstatc0_begin('qitgft_qmk ', id_sname)

     call getttr(rltv,ttr)

     allocate(qitg_exx_l(kgp,nqitg_sp(it))); qitg_exx_l = 0.0d0
     allocate(gqmk_l(kgp))
     mm0 = 0
     do i = 1, it-1
        mm0 = mm0 + nqitg_sp(i)
     end do
     do i = 1, kgp
        kg(1:3) = qmk(ik,1:3) + ngabc(i,1:3)
        g2          = ttr(1)*kg(1)*kg(1) &
        &           + ttr(2)*kg(2)*kg(2) &
        &           + ttr(3)*kg(3)*kg(3) &
        &           + ttr(4)*kg(1)*kg(2) &
        &           + ttr(5)*kg(2)*kg(3) &
        &           + ttr(6)*kg(3)*kg(1)
        gqmk_l(i) = sqrt(g2)
     end do
     do il3 = 1, lcmax+1
        if(nmm_il3(il3) <= 0) cycle
        do i = 1, kgp
           gabs = gqmk_l(i)
           if(gabs < DELTA) then
              idp = nmesh(it)+1
           else
              idp = ceiling(nmesh(it) - dlog(rmax(it)*gabs)/h(it))
           end if
! ==== ASMS ====
           if ( idp > nmesh(it) +1 ) idp = nmesh(it) +1
! ==== ASMS ====

           wkx(1:nmesh(it)) = gabs*radr_kept(1:nmesh(it))
           call dsjnvn(il3-1,nmesh(it),idp,wkx,wky)
           do mm = 1, nmm_il3(il3)
              mmp = mm_il3(mm,il3)
              if ( mmp <= 0 ) cycle

              qitg_sh = 0.d0
              do n = 1, nmesh(it)
                 qitg_sh = qitg_sh + wos_kept(n)*qrsps_mm(n,mmp)*wky(n)
              end do
              qitg_exx_l(i,mmp) = qitg_sh*PAI4
           end do
        end do
     end do
     do iq = 1, nqitg_sp(it)
        do i = 1, kgp
           qitg_exx(i,mm0+iq) = qitg_exx_l(i,iq)
        end do
     end do
     deallocate(qitg_exx_l)
     deallocate(gqmk_l)

     call tstatc0_end(id_sname)
  end subroutine qitgft_qmk_each_k
#endif

end module m_ES_ExactExchange
