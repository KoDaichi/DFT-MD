#define NEC_TUNE
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 606 $)
!
!  MODULE: m_NonLocal_Potential
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
#ifdef __TIMER_SUB__
#   define __TIMER_SUB_START(a)  call timer_sta(a)
#   define __TIMER_SUB_STOP(a)   call timer_end(a)
#else
#   define __TIMER_SUB_START(a)
#   define __TIMER_SUB_STOP(a)
#endif
#ifdef __TIMER_IODO__
#   define __TIMER_IODO_START(a)   call timer_sta(a)
#   define __TIMER_IODO_STOP(a)    call timer_end(a)
#else
#   define __TIMER_IODO_START(a)
#   define __TIMER_IODO_STOP(a)
#endif
#ifdef __TIMER_IOCOMM__
#   define __TIMER_IOCOMM_START_w_BARRIER(str,a)   call timer_barrier(str) ;   call timer_sta(a)
#   define __TIMER_IOCOMM_START(a)       call timer_sta(a)
#   define __TIMER_IOCOMM_STOP(a)        call timer_end(a)
#else
#   define __TIMER_IOCOMM_START_w_BARRIER(str,a)
#   define __TIMER_IOCOMM_START(a)
#   define __TIMER_IOCOMM_STOP(a)
#endif

module m_NonLocal_Potential
!    ( m_NLP )
!  $Id: m_NonLocal_Potential.F90 606 2020-04-15 06:45:49Z ktagami $
!
!  This module contains following major subroutines
!  1.  m_NLP_alloc_snl
!  2.  m_NLP_wd_snl(nfcntn_bin,kv3)
!  3.  m_NLP_rd_snl(nfcntn_bin,kv3)
!  4.  m_NLP_betar_dot_PWs(nfout,kv3,vkxyz)
!  5.  m_NLP_betar_dot_PWs_diff(nfout,kv3,vkxyz)
!
! The subroutine "m_NLP_betar_dot_WFs" originates from a subroutine "kbint".
! The subroutine "kbint" was coded mainly by Y. Morikawa in 1993 or earlier.
! The subroutine "kbint" has following comments
!!!$C Fourier transformation of function beta(r).
!!$c  #1) gx,gy,gz --> ngabc by T.Yamasaki on 15th Feb. 1995
!!$c  #2) vx,vy,vz --> vkxyz by T.Yamasaki on 15th Feb. 1995
!!$c  #3) betar --> betar_l  by T.Yamasaki on 3rd Mar. 1995
!!$c  #4) snl2 is introduced by Y. Morikawa. 14th May. 1996
!!$c  #5) antiferromagnetic calculation is added on 9th Jul. 1996
!!$c                                          by H.Sawada
!
! The subroutine "m_NLP_betar_dot_PWs_diff" originates from a subroutine
! "kbint_diff" coded by H. Sawada at 8th May 1997.
! This was translated into the subroutine "m_NLP_betar_dot_PWs_diff" using
! fortran90+MPI by H. Sawada and T. Yamasaki in 1999.
!
  use m_PlaneWaveBasisSet,    only : kgp,kg1,ngabc,iba,nbase
  use m_PseudoPotential,      only : wos,radr,betar,nmesh,mmesh,ilmt,xh,rmax,ntau &
       &                           , m_PP_tell_lmtt_l_m_tau,nlmtt &
       &                           , m_PP_tell_lmtt_l_m_tau_phi &
       &                           , m_PP_tell_lmtt_l_m_tau_add &
       &                           , m_PP_tell_lmtt_l_m_tau_pao &
       &                           , phirt,ilmt_phi, nlmtt_phi &
       &                           , betar_add,ilmt_add, nlmtt_add &
       &                           , paor,ilmt_pao, nlmtt_pao &
       &                           , nlmta,lmta,ilmt,lmtt, nloc &
       &                           , qrspspw,m_PP_find_maximum_l &
       &                           , ltp,taup,il2p,nqitg,iqitg,isph, modnrm
  use m_Crystal_Structure,    only : rltv, univol, op, altv
  use m_Ionic_System,         only : ntyp,pos,natm,ityp
  use m_Timing,               only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters,   only : nspin, ipri, iprisnl, istress &
       &                           , sw_orb_popu, ipriphig, ipripao, iprirs, printable, r0_factor &
       &                           , dq, gamma_factor,gmax, projector_optimization,gmaxp,ipribetar&
       &                           , kimg,qr_optimization,sw_betar_dot_wfs_exp, sw_precalculate_phase_vnonlocal
  use m_Const_Parameters,     only : DP, PAI2, PAI4, BUCS, CRDTYP, ON, DELTA, PREFITTING, MASK_FUNCTION, PAI, GAMMA
  use m_Parallelization,      only : map_k,myrank_k,ista_k,iend_k,mype &
       &                           , myrank_e,nrank_e,mpi_k_world,ierr &
       &                           , ista_snl, iend_snl, MPI_CommGroup, npes &
       &                           , mpi_spin_group

! ================================ added by K. Tagami ============== 11.0
  use m_Control_Parameters,   only : noncol, ndim_spinor
! ================================================================== 11.0

  use m_FFT,                  only : fft_box_size_WF

  use m_Realspace,            only : nmesh_rs, nmesh_rs_h, nmesh_rs_max, nmesh_rs_max_h &
  &  ,meshx_rs, meshy_rs, meshz_rs, meshxyz_rs, meshxyz_rs_h, map_h, map_h_i, meshxyz_rs_conjg &
  &  ,rcut_betar, rcut_betar_rs, rcut_qr, rcut_qr_rs, m_RS_R_minus_pos

  use m_Kpoints,              only : k_symmetry
  use mpi
  implicit none
!  include 'mpif.h'                              ! MPI
  integer istatus(mpi_status_size)              ! MPI

  real(kind=DP), allocatable,target,dimension(:,:,:)     :: snl  !d(kg1,nlmtt,ista_snl:iend_snl)
  real(kind=DP), allocatable,target,dimension(:,:,:,:,:) :: snld !d(kg1,nlmtt,ista_snl:iend_snl,3)

  real(kind=DP), allocatable,target,dimension(:,:,:,:)   :: i_exp_snl  !d(kg1,nlmtt,ista_snl:iend_snl)
  real(kind=DP), allocatable,target,dimension(:,:,:,:)   :: AtaulmaG,BtaulmaG

  real(kind=DP), allocatable,dimension(:,:,:)     :: snl_add  !d(kg1,nlmtt_add,ista_snl:iend_snl)

  real(kind=DP),private,allocatable,dimension(:) :: qx,qy,qz,vlength,snl2,wka,wkb,ylm
  real(kind=DP),private,allocatable,dimension(:)     :: snl3,wkc,wkd
  real(kind=DP),private,allocatable,dimension(:,:)   :: ylmd,alinvt

  real(kind=DP), allocatable,dimension(:,:,:) :: phig !d(kg1,nlmtt_phi,ista_snl:iend_snl)
  real(kind=DP), allocatable,dimension(:,:) :: norm_phig !d(nlmtt_phi,ista_snl:iend_snl)

  real(kind=DP), allocatable,dimension(:,:,:) :: paog !d(kg1,nlmtt_pao,ista_snl:iend_snl)

  ! rs-related variables
!  integer, allocatable, dimension(:)   :: nmesh_rs
!  integer, allocatable, dimension(:)   :: nmesh_rs_h
!  integer                              :: nmesh_rs_max
!  integer                              :: nmesh_rs_max_h
!  integer, allocatable, dimension(:,:) :: meshx_rs
!  integer, allocatable, dimension(:,:) :: meshy_rs
!  integer, allocatable, dimension(:,:) :: meshz_rs
!  integer, allocatable, dimension(:,:) :: meshxyz_rs
!  integer, allocatable, dimension(:,:) :: meshxyz_rs_h
!  integer, allocatable, dimension(:,:) :: map_h
!  integer, allocatable, dimension(:,:) :: map_h_i
!  integer, allocatable, dimension(:,:) :: meshxyz_rs_conjg
!  integer :: neix=2
!  integer :: neiy=2
!  integer :: neiz=2

  real(kind=DP), allocatable, target, dimension(:,:) :: snl_rs
  real(kind=DP), allocatable, target, dimension(:,:) :: snl_rs_h
  real(kind=DP), allocatable, dimension(:,:,:) :: snld_rs
  real(kind=DP), allocatable, dimension(:,:,:,:) :: betar_optimized
  real(kind=DP), allocatable, dimension(:,:) :: qr_optimized
  logical :: done_optimization = .false.

!  real(kind=DP), allocatable, dimension(:) :: rcut_betar_rs
!  real(kind=DP), allocatable, dimension(:) :: rcut_betar

  real(kind=DP) :: drmask = 0.005d0
  integer, parameter :: nmask=201
  real(kind=DP),dimension(nmask) :: maskr15
  data maskr15 &
  & / 0.10000000E+01, 0.10000000E+01, 0.99948662E+00, 0.99863154E+00, 0.99743557E+00, 0.99589985E+00 &
  &  ,0.99402586E+00, 0.99181538E+00, 0.98927052E+00, 0.98639370E+00, 0.98318766E+00, 0.97965544E+00 &
  &  ,0.97580040E+00, 0.97162618E+00, 0.96713671E+00, 0.96233623E+00, 0.95722924E+00, 0.95182053E+00 &
  &  ,0.94611516E+00, 0.94011842E+00, 0.93383589E+00, 0.92727338E+00, 0.92043693E+00, 0.91333282E+00 &
  &  ,0.90596753E+00, 0.89834777E+00, 0.89048044E+00, 0.88237263E+00, 0.87403161E+00, 0.86546483E+00 &
  &  ,0.85667987E+00, 0.84768450E+00, 0.83848659E+00, 0.82909416E+00, 0.81951535E+00, 0.80975838E+00 &
  &  ,0.79983160E+00, 0.78974340E+00, 0.77950227E+00, 0.76911677E+00, 0.75859548E+00, 0.74794703E+00 &
  &  ,0.73718009E+00, 0.72630334E+00, 0.71532544E+00, 0.70425508E+00, 0.69310092E+00, 0.68187158E+00 &
  &  ,0.67057566E+00, 0.65922170E+00, 0.64781819E+00, 0.63637355E+00, 0.62489612E+00, 0.61339415E+00 &
  &  ,0.60187581E+00, 0.59034914E+00, 0.57882208E+00, 0.56730245E+00, 0.55579794E+00, 0.54431609E+00 &
  &  ,0.53286431E+00, 0.52144984E+00, 0.51007978E+00, 0.49876105E+00, 0.48750040E+00, 0.47630440E+00 &
  &  ,0.46517945E+00, 0.45413176E+00, 0.44316732E+00, 0.43229196E+00, 0.42151128E+00, 0.41083069E+00 &
  &  ,0.40025539E+00, 0.38979038E+00, 0.37944042E+00, 0.36921008E+00, 0.35910371E+00, 0.34912542E+00 &
  &  ,0.33927912E+00, 0.32956851E+00, 0.31999705E+00, 0.31056799E+00, 0.30128436E+00, 0.29214897E+00 &
  &  ,0.28316441E+00, 0.27433307E+00, 0.26565709E+00, 0.25713844E+00, 0.24877886E+00, 0.24057988E+00 &
  &  ,0.23254283E+00, 0.22466884E+00, 0.21695884E+00, 0.20941357E+00, 0.20203357E+00, 0.19481920E+00 &
  &  ,0.18777065E+00, 0.18088790E+00, 0.17417080E+00, 0.16761900E+00, 0.16123200E+00, 0.15500913E+00 &
  &  ,0.14894959E+00, 0.14305240E+00, 0.13731647E+00, 0.13174055E+00, 0.12632327E+00, 0.12106315E+00 &
  &  ,0.11595855E+00, 0.11100775E+00, 0.10620891E+00, 0.10156010E+00, 0.97059268E-01, 0.92704295E-01 &
  &  ,0.88492966E-01, 0.84422989E-01, 0.80492001E-01, 0.76697569E-01, 0.73037197E-01, 0.69508335E-01 &
  &  ,0.66108380E-01, 0.62834685E-01, 0.59684561E-01, 0.56655284E-01, 0.53744102E-01, 0.50948236E-01 &
  &  ,0.48264886E-01, 0.45691239E-01, 0.43224469E-01, 0.40861744E-01, 0.38600231E-01, 0.36437098E-01 &
  &  ,0.34369520E-01, 0.32394681E-01, 0.30509780E-01, 0.28712032E-01, 0.26998673E-01, 0.25366964E-01 &
  &  ,0.23814193E-01, 0.22337676E-01, 0.20934765E-01, 0.19602844E-01, 0.18339338E-01, 0.17141711E-01 &
  &  ,0.16007467E-01, 0.14934157E-01, 0.13919377E-01, 0.12960772E-01, 0.12056034E-01, 0.11202905E-01 &
  &  ,0.10399183E-01, 0.96427132E-02, 0.89313983E-02, 0.82631938E-02, 0.76361106E-02, 0.70482151E-02 &
  &  ,0.64976294E-02, 0.59825322E-02, 0.55011581E-02, 0.50517982E-02, 0.46327998E-02, 0.42425662E-02 &
  &  ,0.38795566E-02, 0.35422853E-02, 0.32293218E-02, 0.29392897E-02, 0.26708663E-02, 0.24227820E-02 &
  &  ,0.21938194E-02, 0.19828122E-02, 0.17886449E-02, 0.16102512E-02, 0.14466132E-02, 0.12967606E-02 &
  &  ,0.11597692E-02, 0.10347601E-02, 0.92089812E-03, 0.81739110E-03, 0.72348823E-03, 0.63847906E-03 &
  &  ,0.56169212E-03, 0.49249371E-03, 0.43028657E-03, 0.37450862E-03, 0.32463165E-03, 0.28016004E-03 &
  &  ,0.24062948E-03, 0.20560566E-03, 0.17468305E-03, 0.14748362E-03, 0.12365560E-03, 0.10287226E-03 &
  &  ,0.84830727E-04, 0.69250769E-04, 0.55873673E-04, 0.44461100E-04, 0.34793983E-04, 0.26671449E-04 &
  &  ,0.19909778E-04, 0.14341381E-04, 0.98138215E-05 /

contains
  subroutine m_NLP_alloc_snl
!!$    ista_snl = (ista_k + nspin - 1)/nspin
!!$    iend_snl = iend_k/nspin
!!$    print '(" ista_snl = ",i3)',ista_snl
!!$    print '(" iend_snl = ",i3)',iend_snl
    if ( allocated(snl) ) deallocate(snl)
    if ( allocated(snld) ) deallocate(snld)
    if ( allocated(i_exp_snl) ) deallocate(i_exp_snl)
!!$#ifndef PARA3D
    allocate(snl(kg1,nlmtt,ista_snl:iend_snl)); snl = 0.d0
!!$#endif
    if(sw_betar_dot_wfs_exp==ON) then
      allocate(i_exp_snl(kg1,nlmta,ista_snl:iend_snl,2)); i_exp_snl = 0.d0
    endif
    if(istress == ON) then
       allocate(snld(kg1,nlmtt,ista_snl:iend_snl,3,3)); snld = 0.d0

    end if
  end subroutine m_NLP_alloc_snl

  subroutine m_NLP_alloc_taulmaG()
    if(sw_precalculate_phase_vnonlocal==ON) then
      allocate(AtaulmaG(kg1,nlmta,ista_k:iend_k,2));AtaulmaG = 0.d0
      if(modnrm==ON) allocate(BtaulmaG(kg1,nlmta,ista_k:iend_k,2));BtaulmaG = 0.d0
    endif
  end subroutine m_NLP_alloc_taulmaG

  subroutine m_NLP_alloc_phig
    if ( allocated(phig) ) deallocate(phig)
    if ( allocated(norm_phig) ) deallocate(norm_phig)
    if(sw_orb_popu == ON) then
       allocate(phig(kg1,nlmtt_phi,ista_snl:iend_snl)); phig = 0.d0
       allocate(norm_phig(nlmtt_phi,ista_snl:iend_snl)); norm_phig = 1.d0
    end if
  end subroutine m_NLP_alloc_phig

  subroutine m_NLP_alloc_snl_add
    if ( allocated(snl_add) ) deallocate(snl_add)
    allocate(snl_add(kg1,nlmtt_add,ista_snl:iend_snl)); snl_add = 0.d0

  end subroutine m_NLP_alloc_snl_add

  subroutine m_NLP_alloc_paog
    if ( allocated(paog) ) deallocate(paog)
    allocate(paog(kg1,nlmtt_pao,ista_snl:iend_snl)); paog = 0.d0
  end subroutine m_NLP_alloc_paog

  subroutine m_NLP_wd_snl(nfout,nfcntn_bin,F_CNTN_BIN_partitioned,kv3)
    integer, intent(in) :: nfout,nfcntn_bin
    logical, intent(in) :: F_CNTN_BIN_partitioned
    integer, intent(in) :: kv3

! ===================================== added by K. Tagami ================ 11.0
    integer :: ikskip
! ========================================================================= 11.0

    integer                                :: i,ik,iksnl,p,q,iend  ! MPI
    real(kind=DP), allocatable, dimension(:,:) :: snl_wk  ! MPI
    integer             :: id_sname = -1

    call tstatc0_begin('m_NLP_wd_snl ',id_sname)

    allocate(snl_wk(kg1,nlmtt))                      ! MPI
!!$ASASASASAS
    snl_wk = 0.d0
!!$ASASASASAS
    if(istress == 0) then
       iend = 0
    else
       iend = 9
    end if
    do i = 0, iend
       if(F_CNTN_BIN_partitioned) then
          do iksnl = ista_snl, iend_snl
             if(i==0) then
                snl_wk = snl(:,:,iksnl)
             else
                p = mod(i,3); if(p == 0) p = 3
                q = (i-p)/3 + 1
                snl_wk = snld(:,:,iksnl,p,q)
             end if
             write(nfcntn_bin) snl_wk
          end do
       else
! ==================================== added by K. Tagami ============== 11.0
          if ( noncol ) then
             ikskip = ndim_spinor
          else
             ikskip = nspin
          endif
!========================================================================== 11.0

! ===================================== modified by K. Tagami ============= 11.0
!          do ik = 1, kv3, nspin                             ! MPI
!
          do ik = 1, kv3, ikskip
! ========================================================================== 11.0

             if(iprisnl >= 2) write(nfout,'(" ! ik = ",i5," <<m_NLP_wd_snl>>")') ik

! ===================================== modified by K. Tagami ============= 11.0
!             iksnl = (ik-1)/nspin + 1                       ! MPI
!
             if ( noncol ) then
                iksnl = (ik-1)/ndim_spinor + 1
             else
                iksnl = (ik-1)/nspin + 1
             endif
! ========================================================================= 11.0

             if(map_k(ik) == myrank_k .and. myrank_e == 0 ) then
                if(i==0) then
                   snl_wk = snl(:,:,iksnl)
                else
                   p = mod(i,3); if(p == 0) p = 3
                   q = (i-p)/3 + 1
                   snl_wk = snld(:,:,iksnl,p,q)
                end if
!                if(map_k(ik) /= 0) call mpi_send(snl_wk,kg1*nlmtt,mpi_double_precision,0,1,MPI_CommGroup,ierr)
                if(map_k(ik) /= 0) call mpi_send(snl_wk,kg1*nlmtt,mpi_double_precision,0,1,mpi_spin_group,ierr)
             else if(mype == 0 .and. map_k(ik) /= 0) then
!                call mpi_recv(snl_wk,kg1*nlmtt,mpi_double_precision,map_k(ik)*nrank_e,1,MPI_CommGroup,istatus,ierr)
                call mpi_recv(snl_wk,kg1*nlmtt,mpi_double_precision,map_k(ik)*nrank_e,1,mpi_spin_group,istatus,ierr)
             end if
             if(mype == 0) write(nfcntn_bin) snl_wk
          end do
       end if
    end do
    if(iprisnl >= 2) then
       write(nfout,'(" ! snl is written (m_NLP_wd_snl)")')
       if(F_CNTN_BIN_partitioned) then
          write(nfout,'(" ! snl size = ",i9)') kg1*nlmtt*(iend_snl-ista_snl+1)*8
       else
          write(nfout,'(" ! snl size = ",i9)') kg1*nlmtt*kv3*8
       end if
    end if
    deallocate(snl_wk)
    call tstatc0_end(id_sname)
  end subroutine m_NLP_wd_snl

  subroutine m_NLP_rd_snl(nfout,nfcntn_bin,F_CNTN_BIN_partitioned,kv3)
    integer, intent(in) :: nfout,nfcntn_bin
    logical, intent(in) :: F_CNTN_BIN_partitioned
    integer, intent(in) :: kv3

! ===================================== added by K. Tagami ================ 11.0
    integer :: ikskip
! ========================================================================= 11.0

    integer                                :: i,ik,iksnl,p,q,iend  ! MPI
    real(kind=DP), allocatable, dimension(:,:) :: snl_wk  ! MPI
    integer             :: id_sname = -1
    call tstatc0_begin('m_NLP_rd_snl ',id_sname)
    allocate(snl_wk(kg1,nlmtt))                      ! MPI
!!$ASASASASAS
    snl_wk = 0.d0
!!$ASASASASAS
    if(iprisnl >= 2) write(nfout,'(" kg1 = ",i6, " nlmtt = ",i6)') kg1,nlmtt
    if(iprisnl >= 2) write(nfout,'(" ! nfcntn_bin = ", i6)') nfcntn_bin
    if(istress == 0) then
       iend = 0
    else
       iend = 9
    end if
    do i = 0, iend                                       ! MPI
       if(F_CNTN_BIN_partitioned) then
          do iksnl = ista_snl, iend_snl
             if(iprisnl >= 2) write(nfout,'(" iksnl = ",i6, " <<m_NLP_rd_snl>>")') iksnl
             read(nfcntn_bin) snl_wk
             if(i==0) then
                snl(:,:,iksnl) = snl_wk
             else
                p = mod(i,3); if(p == 0) p = 3           ! MPI
                q = (i-p)/3 + 1                          ! MPI
                snld(:,:,iksnl,p,q) = snl_wk             ! MPI
             end if
          end do
       else
! ==================================== added by K. Tagami ============== 11.0
          if ( noncol ) then
             ikskip = ndim_spinor
          else
             ikskip = nspin
          endif
!========================================================================== 11.0

! ===================================== modified by K. Tagami ============= 11.0
!          do ik = 1, kv3, nspin                             ! MPI
!
          do ik = 1, kv3, ikskip
! ========================================================================== 11.0

             if(iprisnl >= 2) write(nfout,'(" ik = ",i6, " <<m_NLP_rd_snl>>")') ik

! ===================================== modified by K. Tagami ============= 11.0
!             iksnl = (ik-1)/nspin + 1                       ! MPI
!
             if ( noncol ) then
                iksnl = (ik-1)/ndim_spinor + 1
             else
                iksnl = (ik-1)/nspin + 1
             endif
! ========================================================================= 11.0

             if(mype == 0) read(nfcntn_bin) snl_wk
             if(iprisnl >= 2) write(nfout,'(" snl_wk is read")')
             if(mype == 0 .and. map_k(ik) /= 0) then
!                call mpi_send(snl_wk,kg1*nlmtt,mpi_double_precision,map_k(ik)*nrank_e,1,MPI_CommGroup,ierr)
                call mpi_send(snl_wk,kg1*nlmtt,mpi_double_precision,map_k(ik)*nrank_e,1,mpi_spin_group,ierr)
             else if(map_k(ik) /= 0 .and. map_k(ik) == myrank_k .and. myrank_e == 0) then
!                call mpi_recv(snl_wk,kg1*nlmtt,mpi_double_precision,0,1,MPI_CommGroup,istatus,ierr)
                call mpi_recv(snl_wk,kg1*nlmtt,mpi_double_precision,0,1,mpi_spin_group,istatus,ierr)
             end if
             if(map_k(ik) == myrank_k) then
                call mpi_bcast(snl_wk,kg1*nlmtt,mpi_double_precision,0,mpi_k_world(myrank_k),ierr)
                if(i==0) then
                   snl(:,:,iksnl) = snl_wk
                else
                   p = mod(i,3); if(p == 0) p = 3
                   q = (i-p)/3 + 1
                   snld(:,:,iksnl,p,q) = snl_wk
                end if
             end if
          end do
       end if
    end do
    deallocate(snl_wk)
    call tstatc0_end(id_sname)
  end subroutine m_NLP_rd_snl

  subroutine innerPr_allocate
    allocate(qx(kg1)); qx = 0.d0
    allocate(qy(kg1)); qy = 0.d0
    allocate(qz(kg1)); qz = 0.d0
    allocate(vlength(kg1)); vlength = 0.d0
    allocate(snl2(kg1)); snl2 = 0.d0
    allocate(wka(kg1));  wka  = 0.d0
    allocate(wkb(kg1));  wkb  = 0.d0
    allocate(ylm(kg1));  ylm  = 0.d0
  end subroutine innerPr_allocate

  subroutine innerPr_deallocate
    deallocate(ylm)
    deallocate(wkb)
    deallocate(wka)
    deallocate(snl2)
    deallocate(vlength)
    deallocate(qz)
    deallocate(qy)
    deallocate(qx)
  end subroutine innerPr_deallocate

  subroutine innerPr1_allocate
    allocate(qx(kg1)); qx = 0.d0
    allocate(qy(kg1)); qy = 0.d0
    allocate(qz(kg1)); qz = 0.d0
    allocate(vlength(kg1)); vlength = 0.d0
    allocate(snl2(kg1)); snl2 = 0.d0
    allocate(wka(kg1));  wka  = 0.d0
    allocate(wkb(kg1));  wkb  = 0.d0
    allocate(wkc(kg1));  wkc  = 0.d0
    allocate(wkd(kg1));  wkd  = 0.d0
    allocate(ylm(kg1));  ylm  = 0.d0
    allocate(snl3(kg1)); snl3 = 0.d0
    allocate(ylmd(kg1,3)); ylmd = 0.d0
    allocate(alinvt(3,3)); alinvt = 0.d0
  end subroutine innerPr1_allocate

  subroutine innerPr1_deallocate
    deallocate(alinvt)
    deallocate(ylmd)
    deallocate(snl3)
    deallocate(ylm)
    deallocate(wkd)
    deallocate(wkc)
    deallocate(wkb)
    deallocate(wka)
    deallocate(snl2)
!!$    deallocate(wos)
!!$    deallocate(radr)
    deallocate(vlength)
    deallocate(qz)
    deallocate(qy)
    deallocate(qx)
  end subroutine innerPr1_deallocate

  logical function use_sphr_general( ik, kv3, vkxyz )
    integer,       intent(in)       :: ik, kv3
    real(kind=DP), intent(in)       :: vkxyz(kv3,3,CRDTYP)

    real(kind=DP) :: knorm1, knorm2

    knorm1 = vkxyz(ik,1,BUCS)**2 +vkxyz(ik,2,BUCS)**2 +vkxyz(ik,3,BUCS)**2
    knorm2 = ( vkxyz(ik,1,BUCS) -floor(vkxyz(ik,1,BUCS)) )**2 &
         &  +( vkxyz(ik,2,BUCS) -floor(vkxyz(ik,2,BUCS)) )**2 &
         &  +( vkxyz(ik,3,BUCS) -floor(vkxyz(ik,3,BUCS)) )**2
    knorm1 = sqrt( knorm1 );      knorm2 = sqrt( knorm2 )
    if ( knorm1 > 1.d-20 .and. knorm2 < 1.0D-20 ) then
       use_sphr_general = .true.
    else
       use_sphr_general = .false.
    endif
    return
  end function use_sphr_general

  subroutine m_NLP_betar_dot_PWs(nfout,kv3,vkxyz)
    integer,       intent(in)       :: nfout,kv3
    real(kind=DP), intent(in)       :: vkxyz(kv3,3,CRDTYP)

    real(kind=DP)       :: fac, facr
    integer             :: ik,iksnl,it,n,lmt1,lmtt1,il1,im1,tau1,nspher, ig, n1, n2
    integer             :: id_sname = -1
    real(kind=DP), dimension(:,:), allocatable :: snl2_mpi, snl1

! ===================================== added by K. Tagami ================ 11.0
    integer :: ikskip
! ========================================================================= 11.0

#ifndef _mNLP_no_loop_exchange_
    integer, parameter  :: lcmax = 4
    integer :: ip
    integer :: mil, mp  ! mp = maxval(np)
    integer, dimension(:),     allocatable :: nil  ! d(ntyp)
    integer, dimension(:,:),   allocatable :: nt   ! d(ntyp,mil)
    integer, dimension(:,:,:), allocatable :: tu2p ! d(ntau,ntyp,mil)
    integer, dimension(:),     allocatable :: np   ! d(mil), np=sum_{u=1}^{ntyp}(nt(u,mil))
    integer, dimension(:,:,:), allocatable ::  pm2lmtt1 ! d(mp,2*(mil-1)+1,mil)
    real(kind=DP), dimension(:,:), allocatable :: snl_t, snl_t_mpi ! d(iba(ik),mp)
    real(kind=DP), dimension(:,:), allocatable :: radr_p, wos_p ! d(n1:n2,ntyp)
    integer :: p, mradr_indp
    logical, dimension(:,:), allocatable :: flag_radr !d(ntyp,mil)
    integer, dimension(:), allocatable   :: ip_radr !d(ntyp)
    real(kind=DP) :: r, w
#endif

    call tstatc0_begin('m_NLP_betar_dot_PWs ',id_sname,1)
    call innerPr_allocate()

    fac = PAI4/dsqrt(univol)

#ifndef _mNLP_no_loop_exchange_
!!$    allocate(nil(ntyp)); nil = 0
!!$    allocate(nt(ntyp,lcmax)); allocate(tu2p(ntau,ntyp,lcmax)); allocate(np(lcmax))
!!$    nt = 0; tu2p = 0
!!$    call set_index_arrays1(ntyp,nil,nt,tu2p,mil,mp) ! mil, mp, nil, nt, tu2p, np, pm2lmtt1, contained here
!!$    allocate(pm2lmtt1(mp,2*(mil-1)+1,mil))
!!$    call set_index_arrays2(mp,pm2lmtt1)
    call set_index_arrays() ! mil, mp, nil, nt, tu2p, np, pm2lmtt1, contained here

    n = ceiling(dble(mmesh)/nrank_e)
    if(n == 0) n = 1
    n1 = n*myrank_e + 1
    n2 = n*(myrank_e+1)
    if(n2 > mmesh) n2 = mmesh

    allocate(flag_radr(ntyp,mil))
    allocate(ip_radr(ntyp))
!!$ASASASASAS
    flag_radr = .false.
    ip_radr = 0
!!$ASASASASAS
    call set_flag_radr_and_ip_radr() ! -> mradr_indp,flag_radr, ip_radr

    if(n1 <= n2) then
       allocate(radr_p(n1:n2,mradr_indp))
       allocate(wos_p(n1:n2,mradr_indp))
!!$ASASASASAS
       radr_p = 0; wos_p = 0
!!$ASASASASAS
       call radr_and_wos_p(n1,n2,ntyp,mradr_indp,ip_radr,radr_p,wos_p) ! -> radr_p, wos_p
    end if

! ==================================== added by K. Tagami ============== 11.0
    if ( noncol ) then
       ikskip = ndim_spinor
    else
       ikskip = nspin
    endif
!========================================================================== 11.0

! ===================================== modified by K. Tagami ============= 11.0
!    do ik = 1, kv3, nspin
!
    do ik = 1, kv3, ikskip
! ========================================================================== 11.0

       if(ipri >= 2) write(nfout,'(" ik = ",i8)') ik
       if(map_k(ik) /= myrank_k) cycle                     ! MPI
!!$ASASASASAS
!!$       allocate(snl_t(iba(ik),mp))
!!$       if(nrank_e > 1) allocate(snl_t_mpi(iba(ik),mp))
       allocate(snl_t(iba(ik),mp)); snl_t = 0
       if(nrank_e > 1) then
          allocate(snl_t_mpi(iba(ik),mp)); snl_t_mpi = 0
       endif
!!$ASASASASAS

       if(ipri >= 2) write(nfout,'(" entering k_plus_G_vectors")')
       if(ipri >= 2) then
          write(nfout,'(" iba(ik) = ",i8)') iba(ik)
          write(nfout,'(" vkxyz(ik) = ",3f8.4)') vkxyz(ik,1:3,BUCS)
       end if
       if(.not.allocated(qx)) call phase_error_with_msg(nfout,' qx is not allocated',__LINE__,__FILE__)
       if(.not.allocated(vlength)) call phase_error_with_msg(nfout,' vlength is not allocated',__LINE__,__FILE__)
       if(.not.allocated(nbase)) call phase_error_with_msg(nfout,' nbase is not allocated',__LINE__,__FILE__)
       call k_plus_G_vectors(ik,kgp,kg1,kv3,iba,nbase,vkxyz,ngabc,rltv&
            &,qx,qy,qz,vlength)                            ! ->(bottom_Subr.)

! ===================================== modified by K. Tagami ============= 11.0
!       iksnl = (ik-1)/nspin + 1
!
       if ( noncol ) then
          iksnl = (ik-1)/ndim_spinor + 1
       else
          iksnl = (ik-1)/nspin + 1
       endif
! ========================================================================== 11.0

       do il1 = 1, mil
          snl_t = 0.d0
          do n = n1, n2
             do it = 1, ntyp
                if(il1 > nil(it)) cycle
                if(flag_radr(it,il1)) then
                   r = radr_p(n,ip_radr(it))
                   do ig = 1, iba(ik)
                      wka(ig) = vlength(ig)*r
                   end do
                   call dsjnv(il1-1,iba(ik),wka,wkb)     ! -(bottom_Subr.)
                end if

                w = fac * wos_p(n,ip_radr(it)) * radr_p(n,ip_radr(it))
                do tau1 = 1, nt(it,il1)
                   facr = w * betar(n,il1,tau1,it)
                   ip = tu2p(tau1,it,il1)
                   do ig = 1, iba(ik)
                      snl_t(ig,ip) = snl_t(ig,ip) + facr*wkb(ig)
                   end do
                end do
             end do
          end do
          if(nrank_e > 1) then
             call mpi_allreduce(snl_t,snl_t_mpi,iba(ik)*np(il1),mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
             snl_t = snl_t_mpi
          end if

          do im1 = 1, 2*(il1-1)+1
             nspher = (il1-1)**2 + im1
             if ( use_sphr_general(ik,kv3,vkxyz) ) then
                call sphr_general(iba(ik),nspher,qx,qy,qz,ylm)        ! -(bottom_Subr.)
             else
                call sphr(iba(ik),nspher,qx,qy,qz,ylm)        ! -(bottom_Subr.)
             endif
             do p = 1, np(il1)
                ip = pm2lmtt1(p,im1,il1)
                do ig = 1, iba(ik)
                   snl(ig,ip,iksnl) = snl_t(ig,p)*ylm(ig)
                end do
             end do
          end do
       end do
       if(nrank_e > 1) deallocate(snl_t_mpi)
       deallocate(snl_t)
    end do

    if(n1 <= n2) deallocate(wos_p,radr_p)
    deallocate(ip_radr,flag_radr)
    call dealloc_index_arrays()

#else
    if(nrank_e > 1) then
       allocate(snl1(kg1,nlmtt)); snl1 = 0.d0
       allocate(snl2_mpi(kg1,nlmtt))
!!$ASASASASAS
       snl2_mpi = 0.d0
!!$ASASASASAS
    end if

! ==================================== added by K. Tagami ============== 11.0
    if ( noncol ) then
       ikskip = ndim_spinor
    else
       ikskip = nspin
    endif
!========================================================================== 11.0

! ===================================== modified by K. Tagami ============= 11.0
!    do ik = 1, kv3, nspin
!
    do ik = 1, kv3, ikskip
! ========================================================================== 11.0

       if(map_k(ik) /= myrank_k) cycle                     ! MPI

       call k_plus_G_vectors(ik,kgp,kg1,kv3,iba,nbase,vkxyz,ngabc,rltv&
            &,qx,qy,qz,vlength)                            ! ->(bottom_Subr.)

! ===================================== modified by K. Tagami ============= 11.0
!       iksnl = (ik-1)/nspin + 1
!
       if ( noncol ) then
          iksnl = (ik-1)/ndim_spinor + 1
       else
          iksnl = (ik-1)/nspin + 1
       endif
! ========================================================================== 11.0

       do it=1,ntyp

          n = ceiling(dble(nmesh(it))/nrank_e)
          if(n == 0) n = 1
          n1 = n*myrank_e + 1
          n2 = n*(myrank_e+1)
          if(n2 > nmesh(it)) n2 = nmesh(it)

          call new_radr_and_wos(ik,it)                 ! --> radr, wos
          do lmt1 = 1,ilmt(it)
             call m_PP_tell_lmtt_l_m_tau(lmt1,it,lmtt1,il1,im1,tau1,nspher)
             if ( use_sphr_general(ik,kv3,vkxyz) ) then
                call sphr_general(iba(ik),nspher,qx,qy,qz,ylm)        ! -(bottom_Subr.)
             else
                call sphr(iba(ik),nspher,qx,qy,qz,ylm)        ! -(bottom_Subr.)
             endif
             if(iprisnl >= 2) call wd_lmt_l_m_tau_etc &
                                (nfout,it,lmt1,il1,im1,tau1,lmtt1,nspher)
             snl2 = 0.d0

             do n = n1, n2
                facr = fac*wos(n)*radr(n)*betar(n,il1,tau1,it)
                do ig = 1, iba(ik)
                   wka(ig) = vlength(ig)*radr(n)
                end do
                call dsjnv(il1-1,iba(ik),wka,wkb)     ! -(bottom_Subr.)
                do ig = 1, iba(ik)
                   snl2(ig) = snl2(ig) + facr*wkb(ig)
                end do
             end do

             if(nrank_e == 1) then
                do ig = 1, iba(ik)
                   snl(ig,lmtt1,iksnl) = snl2(ig)*ylm(ig)
                end do
             else
                do ig = 1, iba(ik)
                   snl1(ig,lmtt1) = snl2(ig)*ylm(ig)
                end do
             end if
          end do
       end do

       if(nrank_e > 1) then
          call mpi_allreduce(snl1, snl2_mpi, kg1*nlmtt,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
          do lmtt1 = 1, nlmtt
             do ig = 1, iba(ik)
                snl(ig,lmtt1,iksnl) = snl2_mpi(ig,lmtt1)
             end do
          end do
       end if
    end do
    if(nrank_e > 1) deallocate(snl2_mpi,snl1)
#endif
    if(iprisnl >= 2) call wd_snl

    call innerPr_deallocate
    call tstatc0_end(id_sname)

  contains
#ifndef _mNLP_no_loop_exchange_
    subroutine set_flag_radr_and_ip_radr()
      integer :: it, il, ip, it2, is, it1

      mradr_indp = 1
      ip_radr(1) = 1
      do it = 2, ntyp
         is = 0
         search_loop: do it2 = 1, mradr_indp
            ip = ip_radr(it2)
            if((nmesh(ip) == nmesh(it))) then
               if( dabs(xh(ip)-xh(it)) < DELTA .and. dabs(rmax(ip)-rmax(it)) < DELTA) then
                  is = it2
                  exit search_loop
               end if
            end if
         end do search_loop
         if(is <= 0 .or. is > mradr_indp) then
            mradr_indp = mradr_indp+1
            is = mradr_indp
         end if
         ip_radr(it) = is
      end do

      if(ipri >= 2) then
         write(nfout,'(" !mNLP : mradr_indp = ", i5)') mradr_indp
         do it = 1, ntyp
            write(nfout,'(" !mNLP : ip_radr(",i5,") = ",i5)') it, ip_radr(it)
         end do
      end if

      do il = 1, mil
         flag_radr(1:ntyp,il) = .true.
         it1 = 0
         do it = 1, ntyp
            if(il > nil(it)) cycle
            it1 = it1+1
            if(it1 == 1) then
               flag_radr(it,il) = .true.
            else if(it1 >= 2) then
               if(ip_radr(ip) /= ip_radr(it)) then
                  flag_radr(it,il) = .true.
               else
                  flag_radr(it,il) = .false.
               end if
            end if
            ip = it
         end do
      end do

      if(ipri >= 2) then
         do il = 1, mil
            do it = 1, ntyp
               write(nfout,'(" !mNLP : flag_radr(",i5,",",i5,") = ",l3)') it,il,flag_radr(it,il)
            end do
         end do
      end if
    end subroutine set_flag_radr_and_ip_radr

    subroutine set_index_arrays()
      integer :: it, lmt1, n, il1, im1

      allocate(nil(ntyp)); nil = 0
      allocate(nt(ntyp,lcmax)); allocate(tu2p(ntau,ntyp,lcmax)); allocate(np(lcmax))
      nt = 0; tu2p = 0
      mil = 0
      if(ipri>=2) write(nfout,'(" !mNLP:    it, lmt1: lmtt1, il1, im1, tau1, nspher")')
      do it=1,ntyp
         do lmt1 = 1, ilmt(it)
            call m_PP_tell_lmtt_l_m_tau(lmt1,it,lmtt1,il1,im1,tau1,nspher)
            if(ipri>=2) write(nfout,'(" !mNLP: ",2i5," : ",5i5)') it,lmt1, lmtt1,il1,im1,tau1,nspher
            if(mil < il1) mil = il1
            if(nil(it) < il1) nil(it) = il1
            if(nt(it,il1) < tau1) nt(it,il1) = tau1
         end do
      end do

      do il1=1,mil
         n = 0
         np(il1) = 0
         do it=1,ntyp
            if(il1 > nil(it)) cycle
            np(il1) = np(il1)+nt(it,il1)
            do tau1=1,nt(it,il1)
               n = n+1
               tu2p(tau1,it,il1)= n
            end do
         end do
      end do

      mp = maxval(np(1:mil))
      allocate(pm2lmtt1(mp,2*(mil-1)+1,mil))
!!$ASASASASAS
      pm2lmtt1 = 0
!!$ASASASASAS
      do it = 1, ntyp
         do lmt1=1,ilmt(it)
            call m_PP_tell_lmtt_l_m_tau(lmt1,it,lmtt1,il1,im1,tau1,nspher)
            n = tu2p(tau1,it,il1)
            pm2lmtt1(n,im1,il1) = lmtt1
         end do
      end do

      if(ipri>=2) then
         write(nfout,'(" !mNLP: mil = ",i5)') mil
         do it = 1, ntyp
            write(nfout,'(" !mNLP: nil(",i3,") = ",i5)') it,nil(it)
         end do
         do it = 1, ntyp
            write(nfout,'(" !mNLP: nt(",i3,", : ) = ",5i5)') it,(nt(it,il1),il1=1,nil(it))
         end do

         write(nfout,'(" !mNLP: mp = ",i5)') mp
         do il1 = 1, mil
            write(nfout,'(" !mNLP: np(",i3,")=",i5)') il1, np(il1)
         end do
         do il1 = 1, mil
            do it = 1, ntyp
               if(il1 > nil(it)) cycle
               do tau1 = 1, nt(it,il1)
                  write(nfout,'(" !mNLP: tu2p(",i3,",",i3,",",i3,") = ",i5)') tau1,it,il1,tu2p(tau1,it,il1)
               end do
            end do
         end do
         do il1 = 1, mil
            do im1 = 1, 2*(il1-1)+1
               do n = 1, np(il1)
                  write(nfout,'(" !mNLP: il1 = ",i3," im1 = ",i3," n = ", i3, " pm2lmtt1 = ",i5)') &
                       & il1,im1, n,pm2lmtt1(n,im1,il1)
               end do
            end do
         end do
      end if
    end subroutine set_index_arrays

    subroutine dealloc_index_arrays()
      deallocate(pm2lmtt1,np,tu2p,nt,nil)
    end subroutine dealloc_index_arrays
#endif

    subroutine wd_snl
      integer, parameter :: MSNLSIZE = 20
      integer :: i, ilmtt, ik, j, iksnl, icycle, icolumn, max_elements, istart, iend, ic
      write(nfout,'(" << m_NLP_betar_dot_PWs.wd_snl >>")')
      write(nfout,'(10("(",3i2,")"))') ((ngabc(i,j),j=1,3),i=1,30)

! ==================================== added by K. Tagami ============== 11.0
      if ( noncol ) then
         ikskip = ndim_spinor
      else
         ikskip = nspin
      endif
!========================================================================== 11.0

! ===================================== modified by K. Tagami ============= 11.0
!      do ik = 1, kv3, nspin
!
      do ik = 1, kv3, ikskip
! ========================================================================== 11.0

         if(map_k(ik) /= myrank_k) cycle        ! MPI

! ===================================== modified by K. Tagami ============= 11.0
!         iksnl = (ik-1)/nspin + 1
!
         if ( noncol ) then
            iksnl = (ik-1)/ndim_spinor + 1
         else
            iksnl = (ik-1)/nspin + 1
         endif
! ========================================================================== 11.0

         write(nfout,'(" ik = ",i5)') iksnl
         write(nfout,'(" nbase(1:8,",i5,")",i8,9i12)') ik,(nbase(i,ik),i=1,8)
         if(iprisnl >= 3) then
            max_elements = iba(ik)
         else
            max_elements = min(MSNLSIZE,iba(ik))
         end if
         icolumn = 10
         icycle = ceiling(dble(min(max_elements,kg1))/icolumn)
         do ilmtt = 1, nlmtt
!!$            write(nfout,'("(ilmtt = ",i5,")",8f10.5,99(/15x,8f10.5))') ilmtt,(snl(i,ilmtt,iksnl),i=1,kg1)
            write(nfout,'(" !nlp (ilmtt = ",i5,")")') ilmtt
            istart = 1
            do ic = 1, icycle
               iend = min(istart+icolumn-1,max_elements,kg1)
               write(nfout,'(" !nlp (nx)    ",10i12)') (ngabc(nbase(i,ik),1),i=istart,iend)
               write(nfout,'(" !nlp (ny)    ",10i12)') (ngabc(nbase(i,ik),2),i=istart,iend)
               write(nfout,'(" !nlp (nz)    ",10i12)') (ngabc(nbase(i,ik),3),i=istart,iend)
               write(nfout,'(" !nlp (snl)   ",10d12.4)') (snl(i,ilmtt,iksnl),i=istart,iend)
               istart = iend + 1
            end do
         end do
      end do
    end subroutine wd_snl
  end subroutine m_NLP_betar_dot_PWs

  subroutine m_NLP_betar_dot_PWs_diff(nfout,kv3,vkxyz)
    integer,       intent(in)       :: nfout,kv3
    real(kind=DP), intent(in)       :: vkxyz(kv3,3,CRDTYP)

! ===================================== added by K. Tagami ================ 11.0
    integer :: ikskip
! ========================================================================= 11.0

    real(kind=DP)       :: fac, facr
    integer             :: ik,iksnl,it,n,lmt1,lmtt1,il1,im1,tau1,nspher,j, ig
    integer             :: id_sname = -1

    call tstatc0_begin('m_NLP_betar_dot_PWs_diff ',id_sname,1)
    if(iprisnl >= 2) write(nfout,*) ' <<< m_NLP_betar_dot_PWs_diff >>>'
    call innerPr1_allocate
    if(iprisnl >= 2) write(nfout,*) ' after innerPr1_allocate'

    alinvt = rltv / PAI2
    fac = PAI4/dsqrt(univol)

! ==================================== added by K. Tagami ============== 11.0
    if ( noncol ) then
       ikskip = ndim_spinor
    else
       ikskip = nspin
    endif
!========================================================================== 11.0

! ===================================== modified by K. Tagami ============= 11.0
!    do ik = 1, kv3, nspin
!
    do ik = 1, kv3, ikskip
! ========================================================================== 11.0

       if(map_k(ik) /= myrank_k) cycle             ! MPI

       call k_plus_G_vectors(ik,kgp,kg1,kv3,iba,nbase,vkxyz,ngabc,rltv&
            &,qx,qy,qz,vlength)                            ! ->(bottom_Subr.)

! ===================================== modified by K. Tagami ============= 11.0
!       iksnl = (ik-1)/nspin + 1
!
       if ( noncol ) then
          iksnl = (ik-1)/ndim_spinor + 1
       else
          iksnl = (ik-1)/nspin + 1
       endif
! ========================================================================== 11.0

       do it=1,ntyp
          call new_radr_and_wos(ik,it)                 ! --> radr, wos
          do lmt1 = 1,ilmt(it)
             call m_PP_tell_lmtt_l_m_tau(lmt1,it,lmtt1,il1,im1,tau1,nspher)
             ylm = 0.d0; ylmd = 0.d0
             call sphr(iba(ik),nspher,qx,qy,qz,ylm)        ! -(bottom_Subr.)
             call sphr_diff(kg1,iba(ik),nspher,qx,qy,qz,ylmd) ! -(bottom_Subr.)
             if(ipri >= 2) call wd_lmt_l_m_tau_etc &
                                (nfout,it,lmt1,il1,im1,tau1,lmtt1,nspher)
             snl2 = 0.d0; snl3 = 0.d0
             wkb = 0.d0; wkc = 0.d0; wkd = 0.d0
!xocl spread do/ind_kmesh
             do n = 1,nmesh(it)
                facr = fac*wos(n)*radr(n)*betar(n,il1,tau1,it)
                wka = vlength*radr(n)
                if(il1-1.ne.0) call dsjnv(il1-2,iba(ik),wka,wkc)
                !                                            -(bottom_Subr.)
                call dsjnv(il1-1,iba(ik),wka,wkb)     ! -(bottom_Subr.)
                call dsjnv(il1  ,iba(ik),wka,wkd)     ! -(bottom_Subr.)
                do ig = 1, iba(ik)
                   snl2(ig) = snl2(ig) + facr*wkb(ig)
                end do
!!$                snl2 = snl2 + facr*wkb
                if(il1-1 == 0) then
                   do ig = 1, iba(ik)
                      snl3(ig) = snl3(ig) - facr*radr(n)/(2.d0*il1-1.d0)*il1*wkd(ig)
                   end do
                else
                   do ig =1, iba(ik)
                      snl3(ig) = snl3(ig) + &
                           &  facr*radr(n)/(2.d0*il1-1.d0)*((il1-1)*wkc(ig)-il1*wkd(ig))
                   end do
                endif
             end do
!xocl end spread sum(snl2,snl3)
             snl(1:iba(ik),lmtt1,iksnl) = snl2(1:iba(ik)) * ylm(1:iba(ik))
             where( vlength > 1.d-15 )
                snl3 = - snl3 * ylm / vlength
             elsewhere
                snl3 = 0.d0
             end where

             do j = 1,3
                do ig = 1, iba(ik)
                   snld(ig,lmtt1,iksnl,1,j) = &
                        & snl3(ig)*qx(ig)*(qx(ig)*alinvt(1,j) &
                        &                      +qy(ig)*alinvt(2,j) &
                        &                      +qz(ig)*alinvt(3,j)) &
                        & -snl2(ig)*qx(ig)*(ylmd(ig,1)*alinvt(1,j) &
                        &                       +ylmd(ig,2)*alinvt(2,j) &
                        &                       +ylmd(ig,3)*alinvt(3,j)) &
                        & -snl2(ig)*ylm(ig)/2.d0*alinvt(1,j)
                   snld(ig,lmtt1,iksnl,2,j) = &
                        & snl3(ig)*qy(ig)*(qx(ig)*alinvt(1,j) &
                        &                      +qy(ig)*alinvt(2,j) &
                        &                      +qz(ig)*alinvt(3,j)) &
                        & -snl2(ig)*qy(ig)*(ylmd(ig,1)*alinvt(1,j) &
                        &                       +ylmd(ig,2)*alinvt(2,j) &
                        &                       +ylmd(ig,3)*alinvt(3,j)) &
                        & -snl2(ig)*ylm(ig)/2.d0*alinvt(2,j)
                   snld(ig,lmtt1,iksnl,3,j) = &
                        & snl3(ig)*qz(ig)*(qx(ig)*alinvt(1,j) &
                        &                      +qy(ig)*alinvt(2,j) &
                        &                      +qz(ig)*alinvt(3,j)) &
                        & -snl2(ig)*qz(ig)*(ylmd(ig,1)*alinvt(1,j) &
                        &                       +ylmd(ig,2)*alinvt(2,j) &
                        &                       +ylmd(ig,3)*alinvt(3,j)) &
                        & -snl2(ig)*ylm(ig)/2.d0*alinvt(3,j)
                end do
             end do
          end do
       end do
    end do
    if(ipri >= 2) call wd_snld
    call innerPr1_deallocate
    call tstatc0_end(id_sname)

  contains
    subroutine wd_snld
      integer :: i, ilmtt, ik, j, iksnl, k
      write(nfout,'(" << m_NLP_betar_dot_PWs_diff.wd_snl >>")')
      write(nfout,'(10(''('',3i2,'')''))') ((ngabc(i,j),j=1,3),i=1,30)

! ==================================== added by K. Tagami ============== 11.0
      if ( noncol ) then
         ikskip = ndim_spinor
      else
         ikskip = nspin
      endif
!========================================================================== 11.0

! ===================================== modified by K. Tagami ============= 11.0
!      do ik = 1, kv3, nspin
!
      do ik = 1, kv3, ikskip
! ========================================================================== 11.0

         if(map_k(ik) /= myrank_k) cycle   ! MPI

! ===================================== modified by K. Tagami ============= 11.0
!       iksnl = (ik-1)/nspin + 1
!
         if ( noncol ) then
            iksnl = (ik-1)/ndim_spinor + 1
         else
            iksnl = (ik-1)/nspin + 1
         endif
! ========================================================================== 11.0

         write(nfout,'(" ik = ",i5)') iksnl
         write(nfout,'(8i3)') (nbase(i,ik),i=1,8)
         do ilmtt = 1, nlmtt
	    write(nfout,'('' ilmtt ='',i3)') ilmtt
            write(nfout,'(8f10.5)') (snl(i,ilmtt,iksnl),i=1,8)
              do i = 1,3
              do j = 1,3
              write(nfout,'(8f10.5)') (snld(k,ilmtt,IKSNL,i,j),k=1,8)
              end do
              end do
         end do
      end do
    end subroutine wd_snld
  end subroutine m_NLP_betar_dot_PWs_diff

    subroutine wd_lmt_l_m_tau_etc(nfout,it,lmt1,il1,im1,tau1,lmtt1,nspher)
    integer, intent(in)  :: nfout,it,lmt1,il1,im1,tau1,lmtt1,nspher
      write(nfout,230) it,lmt1,il1,im1,tau1,lmtt1,nspher
230   format(' ',' it=',i2,' lmt1=',i2,' il1,im1,tau1=',3i2, &
           &                   ' lmtt1=',i3,' nspher=',i3)
    end subroutine wd_lmt_l_m_tau_etc

    subroutine new_radr_and_wos(ik,it)
      integer, intent(in)  :: ik,it
      real(kind=DP), parameter :: delta = 1.d-40
      real(kind=DP)      :: hn
      real(kind=DP),save :: xhn = 0.d0 , rmaxn = 0.d0
      integer,save       :: nmeshn = 0
      if((ik == ista_k .and. it == 1) .or. &            ! MPI
           & (nmeshn /= nmesh(it) .or. dabs(xhn-xh(it)) > delta  .or. &
           &  dabs(rmaxn-rmax(it)) > delta) ) then
         call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,hn) ! -(b_PP)
         call coef_simpson_integration(nmesh(it),nmesh(it),xh(it),radr,wos) ! -(b_PP)
         xhn = xh(it); rmaxn = rmax(it); nmeshn = nmesh(it)
      endif
    end subroutine new_radr_and_wos

    subroutine radr_and_wos_p(n1,n2,mtyp,mradr,ip_radr,radr_p,wos_p)
      integer, intent(in) :: n1, n2, mtyp,mradr
      integer, intent(in), dimension(mtyp) :: ip_radr
      real(kind=DP), intent(out), dimension(n1:n2,mradr) :: radr_p, wos_p

      real(kind=DP)      :: hn
      integer :: it, i, ip, ip0
      ip = 0
      do it = 1, ntyp
         if(ntyp > mtyp) cycle
         if(ip_radr(it) /= ip) then
            call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,hn) ! -(b_PP)
            call coef_simpson_integration(nmesh(it),nmesh(it),xh(it),radr,wos) ! -(b_PP)
            ip = ip_radr(it)
            do i = n1, n2
               radr_p(i,ip) = radr(i)
               wos_p(i,ip) = wos(i)
            end do
         end if
      end do
    end subroutine radr_and_wos_p

    subroutine find_critical_point(r,n,x,idp)
      real(kind=DP), intent(in) :: r
      integer,  intent(in)      :: n
      real(kind=DP), intent(in), dimension(n) :: x
      integer,  intent(out)     :: idp
      integer :: i
      do i = 1, n
         idp = i
         if(x(i) > r) exit
      end do
    end subroutine find_critical_point

  subroutine m_NLP_phir_dot_PWs(nfout,kv3,vkxyz)
    integer,       intent(in)       :: nfout,kv3
    real(kind=DP), intent(in)       :: vkxyz(kv3,3,CRDTYP)

! ===================================== added by K. Tagami ================ 11.0
    integer :: ikskip
! ========================================================================= 11.0

    real(kind=DP)       :: fac, facr
    integer             :: ik,ikphig,it,n,lmt1,lmtt1,il1,im1,tau1,nspher
    integer             :: iopr
    integer             :: id_sname = -1

    call tstatc0_begin('m_NLP_phir_dot_PWs ',id_sname,1)
    if(ipriphig >= 2) &
         & write(nfout,*) ' <<< m_NLP_phir_dot_PWs >>> START'
    call innerPr_allocate()

    fac = PAI4/dsqrt(univol)

! ==================================== added by K. Tagami ============== 11.0
    if ( noncol ) then
       ikskip = ndim_spinor
    else
       ikskip = nspin
    endif
!========================================================================== 11.0

! ===================================== modified by K. Tagami ============= 11.0
!    do ik = 1, kv3, nspin
!
    do ik = 1, kv3, ikskip
! ========================================================================== 11.0

       if(map_k(ik) /= myrank_k) cycle                     ! MPI
!!$       write(nfout,*) ' ik = ', ik
!!$       write(nfout,*) ' kgp,kg1,kv3 = ',kgp,kg1,kv3
!!$       write(nfout,*) ' allocated(iba)   = ',allocated(iba)
!!$       write(nfout,*) ' allocated(nbase) = ',allocated(nbase)
!!$       write(nfout,*) ' allocated(ngabc) = ',allocated(ngabc)
!!$       write(nfout,*) ' allocated(qx)    = ',allocated(qx)
!!$       write(nfout,*) ' allocated(vlength)=',allocated(vlength)

       call k_plus_G_vectors(ik,kgp,kg1,kv3,iba,nbase,vkxyz,ngabc,rltv&
               &,qx,qy,qz,vlength)                            ! ->(bottom_Subr.)

! ===================================== modified by K. Tagami ============= 11.0
!     ikphig = (ik-1)/nspin + 1
!
       if ( noncol ) then
          ikphig = (ik-1)/ndim_spinor + 1
       else
          ikphig = (ik-1)/nspin + 1
       endif
! ========================================================================== 11.0

       do it=1,ntyp
          call new_radr_and_wos(ik,it)                 ! --> radr, wos
          do lmt1 = 1,ilmt_phi(it)
             call m_PP_tell_lmtt_l_m_tau_phi(lmt1,it,lmtt1,il1,im1,tau1,nspher)
             if ( use_sphr_general(ik,kv3,vkxyz) ) then
                call sphr_general(iba(ik),nspher,qx,qy,qz,ylm)        ! -(bottom_Subr.)
             else
                call sphr(iba(ik),nspher,qx,qy,qz,ylm)        ! -(bottom_Subr.)
             endif
             if(ipri >= 2) call wd_lmt_l_m_tau_etc &
                                (nfout,it,lmt1,il1,im1,tau1,lmtt1,nspher)
             snl2 = 0.d0

             if(.not.allocated(wos)) then
                write(nfout,'(" wos is not allocated")')
                call phase_error_with_msg(nfout,' wos is not allocated',__LINE__,__FILE__)
             end if
             if(.not.allocated(radr)) then
                write(nfout,'(" radr is not allocated")')
                call phase_error_with_msg(nfout,' wos is not allocated',__LINE__,__FILE__)
             end if
             if(ipri>=2) then
                write(nfout,'(" nmesh(it) = ",i8)') nmesh(it)
             end if
             if(.not.allocated(wka)) call phase_error_with_msg(nfout,' wka is not allocated',__LINE__,__FILE__)
             if(.not.allocated(wkb)) call phase_error_with_msg(nfout,' wkb is not allocated',__LINE__,__FILE__)
             if(.not.allocated(phirt)) call phase_error_with_msg(nfout,' phirt is not allocated',__LINE__,__FILE__)
             do n = 1,nmesh(it)
                facr = fac*wos(n)*radr(n)*phirt(n,il1,tau1,it)
                wka = vlength*radr(n)
!!$                call find_critical_point(1.d0,iba(ik),wka,idp)
!!$                call dsjnvn(il1-1,iba(ik),wka,idp,wkb)     ! -(bottom_Subr.)
                call dsjnv(il1-1,iba(ik),wka,wkb)     ! -(bottom_Subr.)
                snl2 = snl2 + facr*wkb*ylm
             end do
             norm_phig(lmtt1,ikphig) = sum(snl2(1:kg1)*snl2(1:kg1))
             phig(1:kg1,lmtt1,ikphig) = snl2(1:kg1)/sqrt(norm_phig(lmtt1,ikphig))
          end do
       end do
    end do

    if(ipriphig >= 2) call wd_phig

    call innerPr_deallocate
    call tstatc0_end(id_sname)

  contains
    subroutine wd_phig
      integer :: i, ilmtt, ik, j, ikphig, iopr
      write(nfout,'(" << m_NLP_phir_dot_PWs.wd_phig >>")')
      write(nfout,'(" nlmtt_phi=",i3)') nlmtt_phi
      write(nfout,'(10("(",3i2,")"))') ((ngabc(i,j),j=1,3),i=1,30)

! ==================================== added by K. Tagami ============== 11.0
      if ( noncol ) then
         ikskip = ndim_spinor
      else
         ikskip = nspin
      endif
!========================================================================== 11.0

! ===================================== modified by K. Tagami ============= 11.0
!      do ik = 1, kv3, nspin
!
      do ik = 1, kv3, ikskip
! ========================================================================== 11.0

         if(map_k(ik) /= myrank_k) cycle        ! MPI

! ===================================== modified by K. Tagami ============= 11.0
!       ikphig = (ik-1)/nspin + 1
!
         if ( noncol ) then
            ikphig = (ik-1)/ndim_spinor + 1
         else
            ikphig = (ik-1)/nspin + 1
         endif
! ========================================================================== 11.0

         write(nfout,'(" ik = ",i5)') ikphig
         write(nfout,'(" nbase(1:8,",i5,")",8i10)') ik,(nbase(i,ik),i=1,8)
         do ilmtt = 1, nlmtt_phi
            write(nfout,'("(ilmtt = ",i5,")",8f10.5,99(/15x,8f10.5))') ilmtt,(phig(i,ilmtt,ikphig),i=1,kg1)
         end do
         do ilmtt = 1, nlmtt_phi
            write(nfout,'("(ilmtt = ",i5,")",f10.5)') ilmtt,norm_phig(ilmtt,ikphig)
         end do
      end do
      write(nfout,'(" << m_NLP_phir_dot_PWs.wd_phig >> END")')
    end subroutine wd_phig
  end subroutine m_NLP_phir_dot_PWs

  subroutine m_NLP_add_betar_dot_PWs(nfout,kv3,vkxyz)
    integer,       intent(in)       :: nfout,kv3
    real(kind=DP), intent(in)       :: vkxyz(kv3,3,CRDTYP)

! ===================================== added by K. Tagami ================ 11.0
    integer :: ikskip
! ========================================================================= 11.0

    real(kind=DP)       :: fac, facr
    integer             :: ik,iksnl,it,n,lmt1,lmtt1,il1,im1,tau1,nspher
    integer             :: id_sname = -1
#ifdef NEC_TUNE
    integer :: ista, iend, ilen
    real(kind=DP), allocatable, dimension(:,:) :: workarray
! ======== KT_add === 13.0S
    real(kind=DP), allocatable, dimension(:,:) :: workarray2
! =================== 13.0S
#endif

    call tstatc0_begin('m_NLP_add_betar_dot_PWs ',id_sname,1)
!!$    write(nfout,*) ' <<< m_NLP_add_betar_dot_PWs >>>'
    call innerPr_allocate()

    fac = PAI4/dsqrt(univol)

! ==================================== added by K. Tagami ============== 11.0
    if ( noncol ) then
       ikskip = ndim_spinor
    else
       ikskip = nspin
    endif
!========================================================================== 11.0

! ===================================== modified by K. Tagami ============= 11.0
!    do ik = 1, kv3, nspin
!
#ifdef NEC_TUNE
    allocate(workarray(kg1,nlmtt_add))
#endif
    do ik = 1, kv3, ikskip
! ========================================================================== 11.0

       if(map_k(ik) /= myrank_k) cycle                     ! MPI
!!$       write(nfout,*) ' ik = ', ik
!!$       write(nfout,*) ' kgp,kg1,kv3 = ',kgp,kg1,kv3
!!$       write(nfout,*) ' allocated(iba)   = ',allocated(iba)
!!$       write(nfout,*) ' allocated(nbase) = ',allocated(nbase)
!!$       write(nfout,*) ' allocated(ngabc) = ',allocated(ngabc)
!!$       write(nfout,*) ' allocated(qx)    = ',allocated(qx)
!!$       write(nfout,*) ' allocated(vlength)=',allocated(vlength)

       call k_plus_G_vectors(ik,kgp,kg1,kv3,iba,nbase,vkxyz,ngabc,rltv&
            &,qx,qy,qz,vlength)                            ! ->(bottom_Subr.)

! ===================================== modified by K. Tagami ============= 11.0
!     iksnl = (ik-1)/nspin + 1
!
       if ( noncol ) then
          iksnl = (ik-1)/ndim_spinor + 1
       else
          iksnl = (ik-1)/nspin + 1
       endif
! ========================================================================== 11.0

#ifdef NEC_TUNE
       workarray = 0.0d0
#endif
       do it=1,ntyp
#ifdef NEC_TUNE
          ilen = (nmesh(it) - 1)/nrank_e + 1
          ista = ilen*myrank_e + 1
          iend = min(ista + ilen - 1, nmesh(it))
#endif
          call new_radr_and_wos(ik,it)                 ! --> radr, wos
          do lmt1 = 1,ilmt_add(it)
             call m_PP_tell_lmtt_l_m_tau_add(lmt1,it,lmtt1,il1,im1,tau1,nspher)
             if ( use_sphr_general(ik,kv3,vkxyz) ) then
                call sphr_general(iba(ik),nspher,qx,qy,qz,ylm)        ! -(bottom_Subr.)
             else
                call sphr(iba(ik),nspher,qx,qy,qz,ylm)        ! -(bottom_Subr.)
             endif
             if(ipri >= 2) call wd_lmt_l_m_tau_etc &
                                (nfout,it,lmt1,il1,im1,tau1,lmtt1,nspher)
#ifndef NEC_TUNE
             snl2 = 0.d0
             do n = 1,nmesh(it)
#else
             do n = ista, iend
#endif
                facr = fac*wos(n)*radr(n)*betar_add(n,it)
                wka(1:iba(ik)) = vlength(1:iba(ik))*radr(n)
!!$                call find_critical_point(1.d0,iba(ik),wka,idp)
!!$                call dsjnvn(il1-1,iba(ik),wka,idp,wkb)     ! -(bottom_Subr.)
                call dsjnv(il1-1,iba(ik),wka,wkb)     ! -(bottom_Subr.)
#ifndef NEC_TUNE
                snl2(1:iba(ik)) = snl2(1:iba(ik)) + facr*wkb(1:iba(ik))*ylm(1:iba(ik))
#else
                workarray(1:iba(ik),lmtt1) = workarray(1:iba(ik),lmtt1) + facr*wkb(1:iba(ik))*ylm(1:iba(ik))
#endif
             end do
#ifndef NEC_TUNE
             snl_add(1:iba(ik),lmtt1,iksnl) = snl2(1:iba(ik))
#endif
          end do
       end do

#ifdef NEC_TUNE
#ifdef forsafe
       if ( npes > 1 ) then
          allocate( workarray2(kg1,nlmtt_add) )
          call MPI_Allreduce( workarray, workarray2, kg1*nlmtt_add, &
               &              MPI_DOUBLE_PRECISION, MPI_SUM, mpi_k_world(myrank_k),ierr )
          workarray = workarray2
          deallocate( workarray2 )
       endif
#else
       call MPI_Allreduce(MPI_IN_PLACE,workarray,kg1*nlmtt_add,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)
#endif
       do it = 1,ntyp
          do lmt1 = 1, ilmt_add(it)
             call m_PP_tell_lmtt_l_m_tau_add(lmt1,it,lmtt1,il1,im1,tau1,nspher)
             snl_add(1:iba(ik),lmtt1,iksnl) = workarray(1:iba(ik),lmtt1)
          end do
       end do
#endif
    end do

#ifdef NEC_TUNE
    deallocate(workarray)
#endif

    if(iprisnl >= 2) call wd_snl_add
    call innerPr_deallocate
    call tstatc0_end(id_sname)

  contains
    subroutine wd_snl_add
      integer :: i, ilmtt, ik, j, iksnl
      write(nfout,'(" << m_NLP_add_betar_dot_PWs.wd_snl_add >>")')
      write(nfout,'(10("(",3i2,")"))') ((ngabc(i,j),j=1,3),i=1,30)

! ==================================== added by K. Tagami ============== 11.0
      if ( noncol ) then
         ikskip = ndim_spinor
      else
         ikskip = nspin
      endif
!========================================================================== 11.0

! ===================================== modified by K. Tagami ============= 11.0
!      do ik = 1, kv3, nspin
!
      do ik = 1, kv3, ikskip
! ========================================================================== 11.0

         if(map_k(ik) /= myrank_k) cycle        ! MPI

! ===================================== modified by K. Tagami ============= 11.0
!       iksnl = (ik-1)/nspin + 1
!
         if ( noncol ) then
            iksnl = (ik-1)/ndim_spinor + 1
         else
            iksnl = (ik-1)/nspin + 1
         endif
! ========================================================================== 11.0

         write(nfout,'(" ik = ",i5)') iksnl
         write(nfout,'(" nbase(1:8,",i5,")",8i10)') ik,(nbase(i,ik),i=1,8)
         do ilmtt = 1, nlmtt_add
            write(nfout,'("(ilmtt = ",i5,")",8f10.5,99(/15x,8f10.5))') ilmtt,(snl_add(i,ilmtt,iksnl),i=1,kg1)
         end do
      end do
    end subroutine wd_snl_add
  end subroutine m_NLP_add_betar_dot_PWs

  subroutine m_NLP_paor_dot_PWs(nfout,kv3,vkxyz)
    integer,       intent(in)       :: nfout,kv3
    real(kind=DP), intent(in)       :: vkxyz(kv3,3,CRDTYP)

! ===================================== added by K. Tagami ================ 11.0
    integer :: ikskip
! ========================================================================= 11.0

    real(kind=DP)       :: fac, facr
    integer             :: ik,iksnl,it,n,lmt1,lmtt1,il1,im1,tau1,nspher
    integer             :: n1,n2,ig
    integer             :: id_sname = -1
    real(kind=DP), dimension(:,:), allocatable :: snl2_mpi, snl1

#ifndef _mNLP_no_loop_exchange_
    integer, parameter  :: lcmax = 4
    integer :: ip
    integer :: mil, mp  ! mp = maxval(np)
    integer, dimension(:),     allocatable :: nil  ! d(ntyp)
    integer, dimension(:,:),   allocatable :: nt   ! d(ntyp,mil)
    integer, dimension(:,:,:), allocatable :: tu2p ! d(ntau,ntyp,mil)
    integer, dimension(:),     allocatable :: np   ! d(mil), np=sum_{u=1}^{ntyp}(nt(u,mil))
    integer, dimension(:,:,:), allocatable ::  pm2lmtt1 ! d(mp,2*(mil-1)+1,mil)
    real(kind=DP), dimension(:,:), allocatable :: snl_t, snl_t_mpi ! d(iba(ik),mp)
    real(kind=DP), dimension(:,:), allocatable :: radr_p, wos_p ! d(n1:n2,ntyp)
    integer :: p, mradr_indp
    logical, dimension(:,:), allocatable :: flag_radr !d(ntyp,mil)
    integer, dimension(:), allocatable   :: ip_radr !d(ntyp)
    real(kind=DP) :: r, w
#endif

    call tstatc0_begin('m_NLP_paor_dot_PWs ',id_sname,1)
!!$    write(nfout,*) ' <<< m_NLP_paor_dot_PWs >>>'
    call innerPr_allocate()

    fac = PAI4/dsqrt(univol)
#ifndef _mNLP_no_loop_exchange_
    call set_index_arrays() ! mil, mp, nil, nt, tu2p, np, pm2lmtt1, contained here

    n = ceiling(dble(mmesh)/nrank_e)
    if(n == 0) n = 1
    n1 = n*myrank_e + 1
    n2 = n*(myrank_e+1)
    if(n2 > mmesh) n2 = mmesh

    allocate(flag_radr(ntyp,mil))
    allocate(ip_radr(ntyp))
!!$ASASASASAS
    flag_radr = .false. ;  ip_radr = 0
!!$ASASASASAS
    call set_flag_radr_and_ip_radr() ! -> mradr_indp,flag_radr, ip_radr

    if(n1 <= n2) then
       allocate(radr_p(n1:n2,mradr_indp))
       allocate(wos_p(n1:n2,mradr_indp))
!!$ASASASASAS
       radr_p = 0; wos_p = 0
!!$ASASASASAS
       call radr_and_wos_p(n1,n2,ntyp,mradr_indp,ip_radr,radr_p,wos_p) ! -> radr_p, wos_p
    end if

! ==================================== added by K. Tagami ============== 11.0
    if ( noncol ) then
       ikskip = ndim_spinor
    else
        ikskip = nspin
    endif
!========================================================================== 11.0

! ===================================== modified by K. Tagami ============= 11.0
!    do ik = 1, kv3, nspin
!
    do ik = 1, kv3, ikskip
! ========================================================================== 11.0

       if(map_k(ik) /= myrank_k) cycle                     ! MPI
!!$ASASASASAS
!!$       allocate(snl_t(iba(ik),mp))
!!$       if(nrank_e > 1) allocate(snl_t_mpi(iba(ik),mp))
       allocate(snl_t(iba(ik),mp)); snl_t = 0
       if(nrank_e > 1) then
          allocate(snl_t_mpi(iba(ik),mp)); snl_t_mpi = 0
       endif
!!$ASASASASAS

       call k_plus_G_vectors(ik,kgp,kg1,kv3,iba,nbase,vkxyz,ngabc,rltv&
            &,qx,qy,qz,vlength)                            ! ->(bottom_Subr.)

! ===================================== modified by K. Tagami ============= 11.0
!     iksnl = (ik-1)/nspin + 1
!
       if ( noncol ) then
          iksnl = (ik-1)/ndim_spinor + 1
       else
          iksnl = (ik-1)/nspin + 1
       endif
! ========================================================================== 11.0

       do il1 = 1, mil
          snl_t = 0.d0
          do n = n1, n2
             do it = 1, ntyp
                if(il1 > nil(it)) cycle
                if(flag_radr(it,il1)) then
                   r = radr_p(n,ip_radr(it))
                   do ig = 1, iba(ik)
                      wka(ig) = vlength(ig)*r
                   end do
                   call dsjnv(il1-1,iba(ik),wka,wkb)     ! -(bottom_Subr.)
                end if

                w = fac * wos_p(n,ip_radr(it)) * radr_p(n,ip_radr(it))
                do tau1 = 1, nt(it,il1)
                   facr = w * paor(n,il1,tau1,it)
                   ip = tu2p(tau1,it,il1)
                   do ig = 1, iba(ik)
                      snl_t(ig,ip) = snl_t(ig,ip) + facr*wkb(ig)
                   end do
                end do
             end do
          end do
          if(nrank_e > 1) then
             call mpi_allreduce(snl_t,snl_t_mpi,iba(ik)*np(il1),mpi_double_precision &
                  &            ,mpi_sum,mpi_k_world(myrank_k),ierr)
             snl_t = snl_t_mpi
          end if

          do im1 = 1, 2*(il1-1)+1
             nspher = (il1-1)**2 + im1
             if ( use_sphr_general(ik,kv3,vkxyz) ) then
                call sphr_general(iba(ik),nspher,qx,qy,qz,ylm)        ! -(bottom_Subr.)
             else
                call sphr(iba(ik),nspher,qx,qy,qz,ylm)        ! -(bottom_Subr.)
             endif
             do p = 1, np(il1)
                ip = pm2lmtt1(p,im1,il1)
                do ig = 1, iba(ik)
                   paog(ig,ip,iksnl) = snl_t(ig,p)*ylm(ig)
                end do
             end do
          end do
       end do
       if(nrank_e > 1) deallocate(snl_t_mpi)
       deallocate(snl_t)
    end do

    if(n1 <= n2) deallocate(wos_p,radr_p)
    deallocate(ip_radr,flag_radr)
    call dealloc_index_arrays()

#else
    if(nrank_e > 1) then
       allocate(snl1(kg1,nlmtt_pao)); snl1 = 0.d0
       allocate(snl2_mpi(kg1,nlmtt_pao))
!!$ASASASASAS
        snl2_mpi = 0.d0
!!$ASASASASAS
    end if

! ==================================== added by K. Tagami ============== 11.0
    if ( noncol ) then
       ikskip = ndim_spinor
    else
        ikskip = nspin
    endif
!========================================================================== 11.0

! ===================================== modified by K. Tagami ============= 11.0
!    do ik = 1, kv3, nspin
!
    do ik = 1, kv3, ikskip
! ========================================================================== 11.0

       if(map_k(ik) /= myrank_k) cycle                     ! MPI

       call k_plus_G_vectors(ik,kgp,kg1,kv3,iba,nbase,vkxyz,ngabc,rltv&
            &,qx,qy,qz,vlength)                            ! ->(bottom_Subr.)

! ===================================== modified by K. Tagami ============= 11.0
!     iksnl = (ik-1)/nspin + 1
!
       if ( noncol ) then
          iksnl = (ik-1)/ndim_spinor + 1
       else
          iksnl = (ik-1)/nspin + 1
       endif
! ========================================================================== 11.0

       do it=1,ntyp

          n = ceiling(dble(nmesh(it))/nrank_e)
          if(n == 0) n = 1
          n1 = n*myrank_e + 1
          n2 = n*(myrank_e+1)
          if(n2 > nmesh(it)) n2 = nmesh(it)

          call new_radr_and_wos(ik,it)                 ! --> radr, wos
          do lmt1 = 1,ilmt_pao(it)
             call m_PP_tell_lmtt_l_m_tau_pao(lmt1,it,lmtt1,il1,im1,tau1,nspher)
             if ( use_sphr_general(ik,kv3,vkxyz) ) then
                call sphr_general(iba(ik),nspher,qx,qy,qz,ylm)        ! -(bottom_Subr.)
             else
                call sphr(iba(ik),nspher,qx,qy,qz,ylm)        ! -(bottom_Subr.)
             endif
             if(ipri >= 2) call wd_lmt_l_m_tau_etc &
                                (nfout,it,lmt1,il1,im1,tau1,lmtt1,nspher)
             snl2 = 0.d0

             do n = n1, n2
                facr = fac*wos(n)*radr(n)*paor(n,il1,tau1,it)
                do ig = 1, iba(ik)
                   wka(ig) = vlength(ig)*radr(n)
                end do
                call dsjnv(il1-1,iba(ik),wka,wkb)     ! -(bottom_Subr.)
                do ig = 1, iba(ik)
                   snl2(ig) = snl2(ig) + facr*wkb(ig)
                end do
             end do

             if(nrank_e == 1) then
                do ig = 1, iba(ik)
                   paog(ig,lmtt1,iksnl) = snl2(ig)*ylm(ig)
                end do
             else
                do ig = 1, iba(ik)
                   snl1(ig,lmtt1) = snl2(ig)*ylm(ig)
                end do
             end if
          end do
       end do

       if(nrank_e > 1) then
          call mpi_allreduce(snl1, snl2_mpi, kg1*nlmtt_pao,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
          do lmtt1 = 1, nlmtt_pao
             do ig = 1, iba(ik)
                paog(ig,lmtt1,iksnl) = snl2_mpi(ig,lmtt1)
             end do
          end do
       end if
    end do
    if(nrank_e > 1) deallocate(snl2_mpi,snl1)
#endif
    if(ipripao >= 2) call wd_paog
    call innerPr_deallocate
    call tstatc0_end(id_sname)

  contains
#ifndef _mNLP_no_loop_exchange_
    subroutine set_flag_radr_and_ip_radr()
      integer :: it, il, ip, it2, is, it1

      mradr_indp = 1
      ip_radr(1) = 1
      do it = 2, ntyp
         is = 0
         search_loop: do it2 = 1, mradr_indp
            ip = ip_radr(it2)
            if((nmesh(ip) == nmesh(it) .and. dabs(xh(ip)-xh(it)) < DELTA .and. &
              & dabs(rmax(ip)-rmax(it)) < DELTA)) then
               is = it2
               exit search_loop
            end if
         end do search_loop
         if(is <= 0 .or. is > mradr_indp) then
            mradr_indp = mradr_indp+1
            is = mradr_indp
         end if
         ip_radr(it) = is
      end do

      if(ipri >= 2) then
         write(nfout,'(" !mNLP : mradr_indp = ", i5)') mradr_indp
         do it = 1, ntyp
            write(nfout,'(" !mNLP : ip_radr(",i5,") = ",i5)') it, ip_radr(it)
         end do
      end if

      do il = 1, mil
         flag_radr(1:ntyp,il) = .true.
         it1 = 0
         do it = 1, ntyp
            if(il > nil(it)) cycle
            it1 = it1+1
            if(it1 == 1) then
               flag_radr(it,il) = .true.
            else if(it1 >= 2) then
               if(ip_radr(ip) /= ip_radr(it)) then
                  flag_radr(it,il) = .true.
               else
                  flag_radr(it,il) = .false.
               end if
            end if
            ip = it
         end do
      end do

      if(ipri >= 2) then
         do il = 1, mil
            do it = 1, ntyp
               write(nfout,'(" !mNLP : flag_radr(",i5,",",i5,") = ",l3)') it,il,flag_radr(it,il)
            end do
         end do
      end if
    end subroutine set_flag_radr_and_ip_radr

    subroutine set_index_arrays()
      integer :: it, lmt1, n, il1, im1

      allocate(nil(ntyp)); nil = 0
      allocate(nt(ntyp,lcmax)); allocate(tu2p(ntau,ntyp,lcmax)); allocate(np(lcmax))
      nt = 0; tu2p = 0
      mil = 0
      if(ipri>=2) write(nfout,'(" !mNLP:    it, lmt1: lmtt1, il1, im1, tau1, nspher")')
      do it=1,ntyp
         do lmt1 = 1, ilmt_pao(it)
            call m_PP_tell_lmtt_l_m_tau_pao(lmt1,it,lmtt1,il1,im1,tau1,nspher)
            if(ipri>=2) write(nfout,'(" !mNLP: ",2i5," : ",5i5)') it,lmt1, lmtt1,il1,im1,tau1,nspher
            if(mil < il1) mil = il1
            if(nil(it) < il1) nil(it) = il1
            if(nt(it,il1) < tau1) nt(it,il1) = tau1
         end do
      end do

      do il1=1,mil
         n = 0
         np(il1) = 0
         do it=1,ntyp
            if(il1 > nil(it)) cycle
            np(il1) = np(il1)+nt(it,il1)
            do tau1=1,nt(it,il1)
               n = n+1
               tu2p(tau1,it,il1)= n
            end do
         end do
      end do

      mp = maxval(np(1:mil))
      allocate(pm2lmtt1(mp,2*(mil-1)+1,mil))
!!$ASASASASAS
      pm2lmtt1 = 0
!!$ASASASASAS
      do it = 1, ntyp
         do lmt1=1,ilmt_pao(it)
            call m_PP_tell_lmtt_l_m_tau_pao(lmt1,it,lmtt1,il1,im1,tau1,nspher)
            n = tu2p(tau1,it,il1)
            pm2lmtt1(n,im1,il1) = lmtt1
         end do
      end do

      if(ipri>=2) then
         write(nfout,'(" !mNLP: mil = ",i5)') mil
         do it = 1, ntyp
            write(nfout,'(" !mNLP: nil(",i3,") = ",i5)') it,nil(it)
         end do
         do it = 1, ntyp
            write(nfout,'(" !mNLP: nt(",i3,", : ) = ",5i5)') it,(nt(it,il1),il1=1,nil(it))
         end do

         write(nfout,'(" !mNLP: mp = ",i5)') mp
         do il1 = 1, mil
            write(nfout,'(" !mNLP: np(",i3,")=",i5)') il1, np(il1)
         end do
         do il1 = 1, mil
            do it = 1, ntyp
               if(il1 > nil(it)) cycle
               do tau1 = 1, nt(it,il1)
                  write(nfout,'(" !mNLP: tu2p(",i3,",",i3,",",i3,") = ",i5)') tau1,it,il1,tu2p(tau1,it,il1)
               end do
            end do
         end do
         do il1 = 1, mil
            do im1 = 1, 2*(il1-1)+1
               do n = 1, np(il1)
                  write(nfout,'(" !mNLP: il1 = ",i3," im1 = ",i3," n = ", i3, " pm2lmtt1 = ",i5)') &
                       & il1,im1, n,pm2lmtt1(n,im1,il1)
               end do
            end do
         end do
      end if
    end subroutine set_index_arrays

    subroutine dealloc_index_arrays()
      deallocate(pm2lmtt1,np,tu2p,nt,nil)
    end subroutine dealloc_index_arrays
#endif

    subroutine wd_paog
      integer :: i, ilmtt, ik, j, iksnl
      write(nfout,'(" << m_NLP_poor_dot_PWs.wd_paog >>")')
      write(nfout,'(10("(",3i2,")"))') ((ngabc(i,j),j=1,3),i=1,30)

! ==================================== added by K. Tagami ============== 11.0
      if ( noncol ) then
         ikskip = ndim_spinor
      else
          ikskip = nspin
      endif
!========================================================================== 11.0

! ===================================== modified by K. Tagami ============= 11.0
!      do ik = 1, kv3, nspin
!
      do ik = 1, kv3, ikskip
! ========================================================================== 11.0

         if(map_k(ik) /= myrank_k) cycle        ! MPI

! ===================================== modified by K. Tagami ============= 11.0
!       iksnl = (ik-1)/nspin + 1
!
         if ( noncol ) then
            iksnl = (ik-1)/ndim_spinor + 1
         else
            iksnl = (ik-1)/nspin + 1
         endif
! ========================================================================== 11.0

         write(nfout,'(" ik = ",i5)') iksnl
         write(nfout,'(" nbase(1:8,",i5,")",8i10)') ik,(nbase(i,ik),i=1,8)
         do ilmtt = 1, nlmtt_pao
            write(nfout,'("(ilmtt = ",i5,")",8f10.5,99(/15x,8f10.5))') ilmtt,(paog(i,ilmtt,iksnl),i=1,kg1)
         end do
      end do
    end subroutine wd_paog
  end subroutine m_NLP_paor_dot_PWs

  subroutine m_NLP_dealloc
    if(allocated(snl)) deallocate(snl)
    if(allocated(i_exp_snl)) deallocate(i_exp_snl)
    if(allocated(AtaulmaG)) deallocate(AtaulmaG)
    if(allocated(BtaulmaG)) deallocate(BtaulmaG)
    if(allocated(snl_add)) deallocate(snl_add)

    if(allocated(paog)) deallocate(paog)
    if(istress == ON) deallocate(snld)

! =========== KT_add ====== 13.0AS
    if (sw_orb_popu == ON) then
       deallocate(phig);  deallocate(norm_phig)
    end if
! ========================= 13.0AS

  end subroutine m_NLP_dealloc

!===============================================================================

!===============================================================================

  ! snl_rs = betar x Ylm
  ! snld_rs = dbetar x Ylm + betar x dYlm
  ! note: betar is expected to be optimized in some way
  ! this operation must be done every time the atomic coordinate changes
  subroutine m_NLP_build_snl_in_rspace(nfout)
    integer, intent(in) :: nfout
    integer :: ia,it,lmt,lmta1,nma
    integer :: imesh
    real(kind=DP) :: inl,inm,inn
    real(kind=DP) :: bi
    real(kind=DP), allocatable, dimension(:) :: cx,cy,cz,rdiff
    real(kind=DP), allocatable, dimension(:) :: ylm
    real(kind=DP), allocatable, dimension(:,:) :: dylm
    real(kind=DP), allocatable, dimension(:,:,:,:) :: b2
    real(kind=DP), allocatable, dimension(:) :: btmp,b2tmp
    real(kind=DP) :: s,ds,dds,hn
    real(kind=DP) :: rmin
    integer :: ilmtt,il,im,tau,nspher
    integer :: id_sname = -1

    call tstatc0_begin('m_NLP_build_snl_in_rspace ',id_sname,1)

    if(allocated(snl_rs)) deallocate(snl_rs)
    if(allocated(snld_rs)) deallocate(snld_rs)

    allocate(snl_rs(nmesh_rs_max,nlmta))
    allocate(snld_rs(nmesh_rs_max,nlmta,3))

    if(kimg==1)then
       if(allocated(snl_rs_h))then
          deallocate(snl_rs_h)
       endif
    endif
    if(kimg==1) allocate(snl_rs_h(nmesh_rs_max_h,nlmta))

    if(.not.done_optimization)then
       call build_optimized_betar()
       done_optimization = .true.
    endif
    inl = 1.d0/dble(fft_box_size_WF(1,1))
    inm = 1.d0/dble(fft_box_size_WF(2,1))
    inn = 1.d0/dble(fft_box_size_WF(3,1))
    snl_rs = 0.d0
    if(kimg==1) snl_rs_h = 0.d0
    snld_rs = 0.d0
    allocate(cx(nmesh_rs_max));cx=0.d0
    allocate(cy(nmesh_rs_max));cy=0.d0
    allocate(cz(nmesh_rs_max));cz=0.d0
    allocate(rdiff(nmesh_rs_max));rdiff=0.d0
    allocate(ylm(nmesh_rs_max));ylm=0.d0
    allocate(dylm(nmesh_rs_max,3));dylm=0.d0
    allocate(b2(mmesh,nloc,ntau,ntyp));b2=0.d0
    allocate(btmp(mmesh));btmp=0.d0
    allocate(b2tmp(mmesh));b2tmp=0.d0

    do it=1,ntyp
       call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,hn) ! -(b_PP)
       do lmt=1,ilmt(it)
          call m_PP_tell_lmtt_l_m_tau(lmt,it,ilmtt,il,im,tau,nspher)
          call init_cubic_spline( &
        & nmesh(it),radr(1:nmesh(it)),betar_optimized(1:nmesh(it),il,tau,it),b2(1:nmesh(it),il,tau,it))
       enddo
    enddo
!    Loop_atm:do ia=1,natm
!       nma = nmesh_rs(ia)
!       it = ityp(ia)
!       call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,hn) ! -(b_PP)
!       call m_RS_R_minus_pos(pos,ia,nma,inl,inm,inn,cx,cy,cz,rdiff,meshx_rs,meshy_rs,meshz_rs)
!    enddo Loop_atm
    do ia=1,natm
       nma = nmesh_rs(ia)
       it = ityp(ia)
       call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,hn) ! -(b_PP)
       call m_RS_R_minus_pos(pos,ia,nmesh_rs_max,inl,inm,inn,cx,cy,cz,rdiff,meshx_rs,meshy_rs,meshz_rs)
       do lmt=1,ilmt(it)
          call m_PP_tell_lmtt_l_m_tau(lmt,it,ilmtt,il,im,tau,nspher)
          call sphr(nma,nspher,cx(1:nma),cy(1:nma),cz(1:nma),ylm(1:nma))
          call sphr_diff(nma,nma,nspher,cx(1:nma),cy(1:nma),cz(1:nma),dylm(1:nma,1:3))
          lmta1 = lmta(lmt,ia)
          btmp(1:nmesh(it)) = betar_optimized(1:nmesh(it),il,tau,it)
          b2tmp(1:nmesh(it)) = b2(1:nmesh(it),il,tau,it)
          do imesh = 1,nma
             if(rdiff(imesh).gt.radr(1))then
                call cubic_spline(nmesh(it),radr,btmp(1:nmesh(it)),b2tmp(1:nmesh(it)),rdiff(imesh),s,ds)
                snl_rs(imesh,lmta1)    = s*ylm(imesh)
                snld_rs(imesh,lmta1,1) = -(ds*cx(imesh)/rdiff(imesh))*ylm(imesh)-s*dylm(imesh,1)
                snld_rs(imesh,lmta1,2) = -(ds*cy(imesh)/rdiff(imesh))*ylm(imesh)-s*dylm(imesh,2)
                snld_rs(imesh,lmta1,3) = -(ds*cz(imesh)/rdiff(imesh))*ylm(imesh)-s*dylm(imesh,3)
             else
                call cubic_spline(nmesh(it),radr,btmp(1:nmesh(it)),b2tmp(1:nmesh(it)),rdiff(imesh),s,ds)
                snl_rs(imesh,lmta1)    = s*ylm(imesh)
             endif
             if(kimg==1)then
               if(map_h_i(imesh,ia)>0) snl_rs_h(map_h_i(imesh,ia),lmta1) = s*ylm(imesh)
             endif
          enddo
       enddo
    enddo
    deallocate(cx,cy,cz,rdiff,ylm,dylm,b2,btmp,b2tmp)

    call tstatc0_end(id_sname)

  end subroutine m_NLP_build_snl_in_rspace

  subroutine build_optimized_betar()
     real(kind=DP) :: gamm
     integer :: id_sname = -1
     call tstatc0_begin('build_optimized_betar ',id_sname,1)
     if(projector_optimization==MASK_FUNCTION)then
        call optimize_betar_by_maskfunction()
     else if (projector_optimization==PREFITTING) then
        call optimize_betar_by_prefitting()
     else
        betar_optimized = betar ! no optimization of the projector; for debugging purposes
     endif
     call tstatc0_end(id_sname)
  end subroutine build_optimized_betar

  subroutine optimize_betar_by_maskfunction()
     integer :: it,imesh,lmt1,imask
     integer :: lmtt1,il1,im1,tau1,nspher,iq,ier
     integer :: nq
     real(kind=DP) :: hn,r
     real(kind=DP),allocatable,dimension(:) :: betar_tmp
     real(kind=DP),allocatable,dimension(:,:) :: mask
     real(kind=DP),allocatable,dimension(:) :: maskr15_2,rmask
     real(kind=DP) :: fac,topi,m1
     real(kind=DP),allocatable,dimension(:) :: ft,q,waq,wbq,wosq,war,wbr
     nq = int((gamma_factor*gmax)/dq)+1
     allocate(betar_tmp(mmesh));betar_tmp=0.d0
     allocate(ft(nq));ft=0.d0
     allocate(q(nq));q=0.d0
     allocate(waq(nq));waq=0.d0
     allocate(wbq(nq));wbq=0.d0
     allocate(wosq(nq));wosq=0.d0

     do iq=1,nq
        q(iq) = dq*(iq-1)
     enddo
     topi = 2.d0/PAI
     call set_weight_unif(ier,1,nq,q,wosq)

     allocate(mask(mmesh,ntyp));mask=0.d0
     allocate(maskr15_2(nmask));maskr15_2=0.d0
     allocate(rmask(nmask));rmask=0.d0
     do it=1,ntyp
        call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,hn) ! -(b_PP)
        do imask=1,nmask
           rmask(imask) = (imask-1)*rcut_betar_rs(it)*drmask
        enddo
        call init_cubic_spline(nmask,rmask,maskr15,maskr15_2)
        do imesh = 1,nmesh(it)
           r = radr(imesh)
           if(r<=rcut_betar_rs(it)) call cubic_spline(nmask,rmask,maskr15,maskr15_2,r,mask(imesh,it),m1)
        enddo
     enddo
     do it=1,ntyp
        call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,hn) ! -(b_PP)
        call coef_simpson_integration(nmesh(it),nmesh(it),xh(it),radr,wos) ! -(b_PP)
        allocate(war(nmesh(it)));war=0.d0
        allocate(wbr(nmesh(it)));wbr=0.d0
        betar_tmp=0.d0
        do lmt1=1,ilmt(it)
           call m_PP_tell_lmtt_l_m_tau(lmt1,it,lmtt1,il1,im1,tau1,nspher)
           betar_tmp(1:nmesh(it)) = betar(1:nmesh(it),il1,tau1,it)
           ! betar' = betar/maskfunction
           do imesh = 1,nmesh(it)
              r = radr(imesh)
              if(ipribetar>=2) write(100+lmt1,*) r,betar(imesh,il1,tau1,it)
              if (r<=rcut_betar(it)) then
                 if(mask(imesh,it).gt.1.e-7) betar_tmp(imesh) = betar_tmp(imesh)/mask(imesh,it)
              endif
           enddo

           ! FT betar'(r) to q-space
           ft = 0.d0
           waq = 0.0d0
           wbq = 0.0d0
           do imesh=1,nmesh(it)
              fac = wos(imesh)*radr(imesh)*betar_tmp(imesh)
              waq(:) = q(:)*radr(imesh)
              call dsjnv(il1-1,nq,waq,wbq)
              ft(:) = ft(:) + fac*wbq(:)
           enddo

           ! FT betar'(q) back to r-space
           betar_tmp = 0.d0
           war = 0.0d0
           wbr = 0.0d0
           do iq=1,nq
              fac = topi*wosq(iq)*q(iq)*q(iq)*ft(iq)
              war(:) = q(iq)*radr(:)
              call dsjnv(il1-1,nmesh(it),war,wbr)
              betar_tmp(:) = betar_tmp(:) + fac*wbr(:)
           enddo

           ! betar_optimized = betar' * maskfunction
           do imesh = 1,nmesh(it)
              r = radr(imesh)
              betar_optimized(imesh,il1,tau1,it) = betar_tmp(imesh)*mask(imesh,it)
              if(ipribetar>=2) write(200+lmt1,*) r,betar_optimized(imesh,il1,tau1,it)
           enddo

        enddo
        deallocate(war)
        deallocate(wbr)
     enddo
     deallocate(betar_tmp)
     deallocate(ft)
     deallocate(q)
     deallocate(waq)
     deallocate(wbq)
     deallocate(wosq)
     deallocate(mask)
     deallocate(maskr15_2)
     deallocate(rmask)
  end subroutine optimize_betar_by_maskfunction

  subroutine optimize_betar_by_prefitting()
     integer :: i,j,it,lmt1
     integer :: lmtt1,il1,im1,tau1,nspher
     real(kind=DP) :: gamm
     real(kind=DP),allocatable,dimension(:) :: q
     integer :: nq,ngmax
     real(kind=DP) :: hn
     real(kind=DP) :: Aqqp
     real(kind=DP),allocatable,dimension(:) :: wost
     real(kind=DP),allocatable,dimension(:) :: waq,wbq
     real(kind=DP),allocatable,dimension(:) :: war,wbr
     real(kind=DP),allocatable,dimension(:) :: wa,wb,wc,wd
     real(kind=DP),allocatable,dimension(:,:) :: A
     real(kind=DP),allocatable,dimension(:,:) :: B
     real(kind=DP),allocatable,dimension(:) :: betar_q
     real(kind=DP) :: fac
     real(kind=DP) :: topi,pih
     integer,allocatable,dimension(:) :: ipiv
     integer :: info,ier

     gamm = gmax*gamma_factor
     nq = int(gamm/dq)+1
     ngmax = int(gmax/dq)+1
     topi = 2.d0/PAI
     pih  = 0.5d0*PAI
     allocate(q(nq))
     do i=1,nq
        q(i) = (i-1)*dq
     enddo

     allocate(waq(nq));waq=0.d0
     allocate(wbq(nq));wbq=0.d0
     allocate(wa(nq));wa=0.d0
     allocate(wb(nq));wb=0.d0
     allocate(wc(nq));wc=0.d0
     allocate(wd(nq));wd=0.d0
     allocate(wost(nq));wost=0.d0
     allocate(A(nq-ngmax,nq-ngmax));A=0.d0
     allocate(B(nq-ngmax,1));B=0.d0
     allocate(betar_q(nq));betar_q=0.d0
     allocate(ipiv(nq-ngmax));ipiv=0
     allocate(war(mmesh));war=0.d0
     allocate(wbr(mmesh));wbr=0.d0
     do it=1,ntyp
        call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,hn) ! -(b_PP)
        call coef_simpson_integration(nmesh(it),nmesh(it),xh(it),radr,wos) ! -(b_PP)
        wa(:) = rcut_betar_rs(it)*q(:)
        do lmt1=1,ilmt(it)
           call m_PP_tell_lmtt_l_m_tau(lmt1,it,lmtt1,il1,im1,tau1,nspher)

           ! FT to reciprocal space
           betar_q = 0.d0
           do i=1,nmesh(it)
              fac = wos(i)*radr(i)*betar(i,il1,tau1,it)
              waq(:) = q(:)*radr(i)
              call dsjnv(il1-1,nq,waq,wbq)
              betar_q(:) = betar_q(:)+fac*wbq(:)
              if(ipribetar>=2) write(100+lmtt1,'(2f25.15)') radr(i),betar(i,il1,tau1,it)
           enddo

           ! prepare the sph. Bessel functions necessary for the calculation of A(q,q')
           if(il1>=2) then
              call dsjnv(il1-2,nq,wa,wb)
           else
              do i=1,nq
                 wb(i) = dcos(wa(i))/wa(i) ! -n0(x)
              enddo
           endif
           call dsjnv(il1-1,nq,wa,wc)
           call dsjnv(il1,nq,wa,wd)

           call set_weight_unif(ier,1,nq,q,wost)

           ! build the matrix B
           B = 0.d0
           do i=ngmax+1,nq
              do j=1,ngmax
                 call calculate_Aqqp(il1,nq,i,j,q,wb,wc,wd,rcut_betar_rs(it),Aqqp)
                 B(i-ngmax,1) = B(i-ngmax,1)+Aqqp*betar_q(j)*wost(j)
              enddo
           enddo

           ! build the matrix A
           A = 0.d0
           do i=ngmax+1,nq
              do j=ngmax+1,nq
                 call calculate_Aqqp(il1,nq,i,j,q,wb,wc,wd,rcut_betar_rs(it),Aqqp)
                 A(i-ngmax,j-ngmax) = -wost(j)*Aqqp
                 if(i==j) then
                    A(i-ngmax,i-ngmax) = A(i-ngmax,i-ngmax)+pih*q(i)*q(i)
                 endif
              enddo
           enddo

           if(ipribetar>=2)then
              do i=1,nq
                 write(200+lmtt1,'(2f25.15)') q(i),betar_q(i)
              enddo
           endif

           ! solve the linear equation, B = AX;
           call dgesv(nq-ngmax,1,A,nq-ngmax,ipiv,B,nq-ngmax,info)
           do i=ngmax+1,nq
              betar_q(i) = B(i-ngmax,1)
           enddo

           if(ipribetar>=2)then
              do i=1,nq
                 write(300+lmtt1,'(2f25.15)') q(i),betar_q(i)
              enddo
           endif

           ! FT back to Rspace
           betar_optimized(:,il1,tau1,it) = 0.d0
           war = 0.0d0
           wbr = 0.0d0
           do i=1,nq
              fac = topi*wost(i)*q(i)*q(i)*betar_q(i)
              war(:) = q(i)*radr(:)
              call dsjnv(il1-1,nmesh(it),war,wbr)
              betar_optimized(:,il1,tau1,it) = betar_optimized(:,il1,tau1,it) + fac*wbr(:)
           enddo
           if(ipribetar>=2)then
              do i=1,nmesh(it)
                 write(400+lmtt1,'(2f25.15)') radr(i),betar_optimized(i,il1,tau1,it)
              enddo
           endif
        enddo
     enddo
     deallocate(q)
     deallocate(war)
     deallocate(wbr)
     deallocate(waq,wbq)
     deallocate(wa,wb,wc,wd)
     deallocate(wost)
     deallocate(A)
     deallocate(B)
     deallocate(betar_q)
     deallocate(ipiv)
  end subroutine optimize_betar_by_prefitting

  ! calculate A(q,q') from direct integration; for debugging purposes
  subroutine calculate_Aqqp_directly(n,il,radr,wos,q,qp,R0,Aqqp)
     integer,intent(in) :: n
     integer,intent(in) :: il
     real(kind=DP),dimension(mmesh),intent(in) :: radr,wos
     real(kind=DP),intent(in) :: q,qp,R0
     real(kind=DP),intent(out) :: Aqqp
     real(kind=DP),allocatable,dimension(:) :: wa,wb,wc
     integer :: i

     allocate(wa(n))
     allocate(wb(n))
     allocate(wc(n))

     wa(1:n) = q*radr(1:n)
     call dsjnv(il-1,n,wa,wb)
     wa(1:n) = qp*radr(1:n)
     call dsjnv(il-1,n,wa,wc)

     Aqqp = 0.d0
     do i=1,n
        Aqqp = Aqqp+radr(i)*radr(i)*wb(i)*wc(i)*wos(i)
     enddo
     Aqqp = Aqqp*q*q*qp*qp

     deallocate(wa,wb,wc)
  end subroutine calculate_Aqqp_directly

  ! calculate A(q,q') from an analytical expression
  subroutine calculate_Aqqp(il,n,i,ip,q,jlm1,jl0,jl1,R0,Aqqp)
     integer,intent(in) :: il
     integer,intent(in) :: n
     integer,intent(in) :: i,ip
     real(kind=DP),dimension(n),intent(in) :: q
     real(kind=DP),dimension(n) :: jlm1,jl0,jl1
     real(kind=DP),intent(in) :: R0
     real(kind=DP),intent(out) :: Aqqp
     real(kind=DP) :: jl,jlp,djl,djlp
     real(kind=DP) :: qq,qqp
     real(kind=DP) :: q2,fac
     qq  = q(i)
     qqp = q(ip)
     jl   = jl0(i)
     jlp  = jl0(ip)
     djl  = (dble(il)/(R0))*jl0(i)-qq*jl1(i)
     djlp = (dble(il)/(R0))*jl0(ip)-qqp*jl1(ip)
     if(i.ne.ip)then
        Aqqp = (R0*R0/(qqp*qqp-qq*qq))*(djl*jlp-djlp*jl)
     else
        Aqqp = 0.5d0*R0*R0*R0*(jl*jl-jl1(i)*jlm1(i))
     endif
     Aqqp = Aqqp*qq*qq*qqp*qqp
  end subroutine calculate_Aqqp

  subroutine m_NLP_alloc_betar_optimized()
     allocate(betar_optimized(mmesh,nloc,ntau,ntyp));betar_optimized=0.d0
     betar_optimized(:,:,:,:) = betar(:,:,:,:)
  end subroutine m_NLP_alloc_betar_optimized

  subroutine m_NLP_alloc_qr_optimized()
     integer :: id_sname = -1
     integer :: it,i,iq
     real(kind=DP) :: rr,hn
     allocate(qr_optimized(mmesh,nqitg));qr_optimized=0.d0
     do it=1,ntyp
        call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,hn) ! -(b_PP)
        do iq=1,nqitg
           do i=1,nmesh(it)
              rr = radr(i)*radr(i)
              qr_optimized(i,iq) = qrspspw(i,iq)/rr
           enddo
        enddo
     enddo
  end subroutine m_NLP_alloc_qr_optimized

  subroutine m_NLP_cal_i_l_exp_snl()
    complex(kind=DP), allocatable, dimension(:,:) :: wkexp
    integer :: ia,it,ig,iksnl,ik,lmt1,lmtt1,lmta1,il1,mil
    complex(kind=DP) :: ctmp,i1=(0.d0,1.d0)
    real(kind=DP) :: fac
    integer :: id_sname = -1
    call tstatc0_begin('m_NLP_cal_i_l_exp_snl ',id_sname,1)
    allocate(wkexp(maxval(iba),natm));wkexp=0.d0
    do iksnl=ista_snl,iend_snl
      if(noncol)then
        ik = (iksnl-1)*ndim_spinor+1
      else
        ik = (iksnl-1)*nspin+1
      endif
      fac = 1.d0
      if(k_symmetry(ik)==GAMMA) fac = 2.d0
      call cal_phase()
      do ia=1,natm
        it = ityp(ia)
        do lmt1 = 1, ilmt(it)
          do ig=1,iba(ik)
            lmtt1 = lmtt(lmt1,it)
            lmta1 = lmta(lmt1,ia)
            il1   = ltp(lmt1,it)
            mil   = mod(il1-1,4)
            ctmp  = (i1**mil)*wkexp(ig,ia)*snl(ig,lmtt1,iksnl)
            i_exp_snl(ig,lmta1,iksnl,1) = fac *  real(ctmp)
            i_exp_snl(ig,lmta1,iksnl,2) = fac * dimag(ctmp)
            if(k_symmetry(ik)==GAMMA .and. ig==1) then
              if(mil==0 .or. mil==2) then
                i_exp_snl(1,lmta1,iksnl,1) = 0.d0
                i_exp_snl(1,lmta1,iksnl,2) = 0.d0
              else if (mil==1 .or. mil==3) then
                i_exp_snl(1,lmta1,iksnl,1) = 0.5d0*i_exp_snl(1,lmta1,iksnl,1)
                i_exp_snl(1,lmta1,iksnl,2) = 0.d0
              endif
            endif
          enddo
        enddo
      enddo
    enddo

    deallocate(wkexp)
    call tstatc0_end(id_sname)

    contains

    subroutine cal_phase()
      real(kind=DP) :: f1,f2,f3,ph
      integer :: nb
      do ia=1, natm
        f1 = pos(ia,1)*PAI2; f2 = pos(ia,2)*PAI2; f3 = pos(ia,3)*PAI2
        do ig=1,iba(ik)
          nb = nbase(ig,ik)
          ph = ngabc(nb,1)*f1+ngabc(nb,2)*f2+ngabc(nb,3)*f3
          wkexp(ig,ia) = dcmplx(dcos(ph),dsin(ph))
        enddo
      enddo
    end subroutine cal_phase

  end subroutine m_NLP_cal_i_l_exp_snl

!  subroutine m_NLP_dealloc_mesh_for_rspace()
!     deallocate(meshx_rs)
!     deallocate(meshy_rs)
!     deallocate(meshz_rs)
!     deallocate(meshxyz_rs)
!     if(kimg==1) then
!        deallocate(meshxyz_rs_h)
!        deallocate(map_h)
!        deallocate(map_h_i)
!        deallocate(meshxyz_rs_conjg)
!        deallocate(snl_rs_h)
!     endif
!     deallocate(snl_rs)
!     deallocate(snld_rs)
!  end subroutine m_NLP_dealloc_mesh_for_rspace
!
!  subroutine m_NLP_resolve_mesh_for_rspace(nfout)
!     integer, intent(in) :: nfout
!     integer :: nl,nm,nn,mm,id,nlhf
!     real(kind=DP) :: inl,inm,inn
!     real(kind=DP) :: ntot
!     integer :: id_sname = -1
!     logical,save :: firstcall=.true.
!     call tstatc0_begin('m_NLP_resolve_mesh_for_rspace ',id_sname,1)
!     if (firstcall) then
!        call resolve_cutoff_of_betar(nfout)
!        if(.not.allocated(nmesh_rs)) allocate(nmesh_rs(natm))
!        if(kimg==1) then
!           if(.not.allocated(nmesh_rs_h)) allocate(nmesh_rs_h(natm))
!        endif
!     endif
!     id = fft_box_size_WF(1,0)
!     mm = fft_box_size_WF(2,0)
!     nl = fft_box_size_WF(1,1);inl = 1.d0/dble(nl)
!     nm = fft_box_size_WF(2,1);inm = 1.d0/dble(nm)
!     nn = fft_box_size_WF(3,1);inn = 1.d0/dble(nn)
!     if(kimg==1)then
!       nlhf = id/2
!     else
!       nlhf = id
!     endif
!     ntot = dble(nl*nm*nn)
!
!     if(.not.firstcall)then
!       call m_NLP_dealloc_mesh_for_rspace()
!     endif
!     nmesh_rs=0
!     if(kimg==1) nmesh_rs_h=0
!     call resolve_atom_centered_mesh(prealloc=.true.,maxcount=nmesh_rs_max,maxcount_h=nmesh_rs_max_h)
!     allocate(meshx_rs(nmesh_rs_max,natm));meshx_rs = -1
!     allocate(meshy_rs(nmesh_rs_max,natm));meshy_rs = -1
!     allocate(meshz_rs(nmesh_rs_max,natm));meshz_rs = -1
!     allocate(meshxyz_rs(nmesh_rs_max,natm));meshxyz_rs = -1
!
!     if(kimg==1) then
!        allocate(meshxyz_rs_conjg(nmesh_rs_max,natm));meshxyz_rs_conjg = 1
!        allocate(meshxyz_rs_h(nmesh_rs_max_h,natm));meshxyz_rs_h = -1
!        allocate(map_h(nmesh_rs_max_h,natm));map_h = -1
!        allocate(map_h_i(nmesh_rs_max,natm));map_h_i = -1
!     endif
!     allocate(snl_rs(nmesh_rs_max,nlmta));snl_rs=0.d0
!     if(kimg==1) then
!        allocate(snl_rs_h(nmesh_rs_max_h,nlmta));snl_rs_h=0.d0
!     endif
!     allocate(snld_rs(nmesh_rs_max,nlmta,3));snld_rs=0.d0
!     call resolve_atom_centered_mesh(prealloc=.false.)
!
!     if(iprirs>=2.and.printable)then
!         call print_mesh()
!     endif
!     firstcall = .false.
!
!     call tstatc0_end(id_sname)
!
!     contains
!
!     subroutine resolve_atom_centered_mesh(prealloc,maxcount,maxcount_h)
!        logical, intent(in) :: prealloc
!        integer, optional, intent(out) :: maxcount,maxcount_h
!        real(kind=DP) :: rcut_betar_max,rcut,rcut2,rr,rcutmax
!        real(kind=DP) :: ex,ey,ez
!        integer :: nl_per_atm,nm_per_atm,nn_per_atm
!        integer :: nnx,nny,nnz
!        integer :: it,ia
!        integer :: i1,j1,k1,ii1,jj1,kk1,iil,iim,iin,cl,cm,cn
!        integer :: i1min,i1max,j1min,j1max,k1min,k1max
!        integer :: icount,icounth
!        real(kind=DP) :: fl,fm,fn
!        real(kind=DP) :: dx,dy,dz,r2,r,cx,cy,cz
!        integer :: maxc,maxch
!        integer :: pm
!        logical :: smallx,smally,smallz,tof
!
!        rcut_betar_max = 0.d0
!        do it=1,ntyp
!           if(rcut_betar_max<rcut_betar_rs(it)) rcut_betar_max = rcut_betar_rs(it)
!        enddo
!
!        ex = neix * dsqrt(altv(1,1)**2+altv(1,2)**2+altv(1,3)**2)
!        ey = neiy * dsqrt(altv(2,1)**2+altv(2,2)**2+altv(2,3)**2)
!        ez = neiz * dsqrt(altv(3,1)**2+altv(3,2)**2+altv(3,3)**2)
!
!        nnx = floor(ex/rcut_betar_max)
!        nny = floor(ey/rcut_betar_max)
!        nnz = floor(ez/rcut_betar_max)
!        smallx = nnx<2
!        smally = nny<2
!        smallz = nnz<2
!
!        nl_per_atm = floor(dble((2*neix+1)*nl)/dble(nnx))+1
!        nm_per_atm = floor(dble((2*neiy+1)*nm)/dble(nny))+1
!        nn_per_atm = floor(dble((2*neiz+1)*nn)/dble(nnz))+1
!        if(iprirs>=2)then
!          write(nfout,'(a,3i8)') ' !RS number of fft elements surrounding an atom  : ', &
!          & nl_per_atm,nm_per_atm,nn_per_atm
!        endif
!        maxc = 0
!        if(kimg==1) maxch = 0
!        do ia=1,natm
!           if(smallx)then
!              i1min = -nl*neix
!              i1max = +nl*(neix+1)
!           else
!              cl = -neix*nl+floor((2*neix+1)*nl*dble(neix+pos(ia,1))/dble(2*neix+1))
!              i1min = cl-nl_per_atm
!              i1max = cl+nl_per_atm
!           endif
!           if(smally)then
!              j1min = -nm*neiy
!              j1max = +nm*(neiy+1)
!           else
!              cm = -neiy*nm+floor((2*neiy+1)*nm*dble(neiy+pos(ia,2))/dble(2*neiy+1))
!              j1min = cm-nm_per_atm
!              j1max = cm+nm_per_atm
!           endif
!           if(smallz)then
!              k1min = -nn*neiz
!              k1max = +nn*(neiz+1)
!           else
!              cn = -neiz*nn+floor((2*neiz+1)*nn*dble(neiz+pos(ia,3))/dble(2*neiz+1))
!              k1min = cn-nn_per_atm
!              k1max = cn+nn_per_atm
!           endif
!           if(iprirs>=2.and. .not. (smallx.or.smally.or.smallz)) &
!           &  write(nfout,'(a,i8,a,3i8)') ' !RS center of mesh for atom ',ia,':',cl,cm,cn
!           it = ityp(ia)
!           rcut2=rcut_betar_rs(it)*rcut_betar_rs(it)
!           icount = 0
!           if(kimg==1) icounth = 0
!           do i1=i1min,i1max
!           do j1=j1min,j1max
!           do k1=k1min,k1max
!              iil = i1-1
!              iim = j1-1
!              iin = k1-1
!
!              ii1 = mod(iil+(2*neix+1)*nl,nl)+1
!              jj1 = mod(iim+(2*neiy+1)*nm,nm)+1
!              kk1 = mod(iin+(2*neiz+1)*nn,nn)+1
!
!              fl = dble(iil) * inl
!              fm = dble(iim) * inm
!              fn = dble(iin) * inn
!
!              dx = fl-pos(ia,1)
!              dy = fm-pos(ia,2)
!              dz = fn-pos(ia,3)
!
!              cx = altv(1,1)*dx+altv(1,2)*dy+altv(1,3)*dz
!              cy = altv(2,1)*dx+altv(2,2)*dy+altv(2,3)*dz
!              cz = altv(3,1)*dx+altv(3,2)*dy+altv(3,3)*dz
!              r2 = cx*cx+cy*cy+cz*cz
!              if(r2<rcut2)then
!                 tof = .false.
!                 icount = icount+1
!                 if(kimg==1.and.ii1<=nlhf) then
!                    icounth = icounth+1
!                    tof = .true.
!                 endif
!                 if(.not.prealloc)then
!                    meshx_rs(icount,ia) = iil
!                    meshy_rs(icount,ia) = iim
!                    meshz_rs(icount,ia) = iin
!                    if(kimg==1)then
!                       if(ii1>nlhf)then
!                          ii1 = id - ii1
!                          jj1 = nm+2 - jj1
!                          kk1 = nn+2 - kk1
!                          if(jj1>nm) jj1 = jj1-nm
!                          if(jj1.eq.0) jj1=1
!                          if(kk1>nn) kk1 = kk1-nn
!                          if(kk1.eq.0) kk1=1
!                          pm = +1.0d0
!                       else
!                          pm = -1.d0
!                       endif
!                    endif
!                    meshxyz_rs(icount,ia) = nlhf*mm*(kk1-1)+nlhf*(jj1-1)+ii1
!                    if(kimg==1) then
!                       meshxyz_rs_conjg(icount,ia) = pm
!                       if(tof) then
!                          meshxyz_rs_h(icounth,ia) = meshxyz_rs(icount,ia)
!                          map_h(icounth,ia) = icount
!                          map_h_i(icount,ia) = icounth
!                       endif
!                    endif
!                 endif
!              endif
!           enddo
!           enddo
!           enddo
!           if(prealloc)then
!             if(icount>maxc) maxc = icount
!             nmesh_rs(ia) = icount
!             if(kimg==1) then
!                if(icounth>maxch) maxch = icounth
!                nmesh_rs_h(ia) = icounth
!             endif
!           endif
!        enddo
!        if(prealloc.and.present(maxcount))then
!           maxcount = maxc
!        endif
!        if(prealloc.and.kimg==1.and.present(maxcount_h))then
!           maxcount_h = maxch
!        endif
!
!     end subroutine resolve_atom_centered_mesh
!
!     subroutine print_mesh()
!        integer :: ia,nm
!        integer :: totnmesh
!        totnmesh = 0
!        do ia=1,natm
!           write(nfout,'(a,i0,a,i0)') ' !RS number of mesh points associated to atom ',ia,' : ',nmesh_rs(ia)
!           if(iprirs>=3)then
!             write(nfout,'(a)')              ' !RS associated mesh points ... '
!             do nm = 1,nmesh_rs(ia)
!                 write(nfout,'(a,4i8)')      ' !RS ',meshx_rs(nm,ia),meshy_rs(nm,ia),meshz_rs(nm,ia),meshxyz_rs(nm,ia)
!             enddo
!           endif
!           totnmesh = totnmesh+nmesh_rs(ia)
!        enddo
!        totnmesh = totnmesh/natm
!        write(nfout,'(a,i0,a,i0,a,f10.5)') ' !RS average number of mesh points associated to an atom / total FFT mesh = ' &
!        &     ,totnmesh,'/',int(ntot),' = ',dble(totnmesh)/ntot
!     end subroutine print_mesh
!
!  end subroutine m_NLP_resolve_mesh_for_rspace

end module m_NonLocal_Potential
