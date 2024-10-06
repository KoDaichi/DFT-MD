!#define _DEBUG_WRITE_DFTU_MPI_PROCESSES_
!#define _DUPLICATION_HSR_DOTPRODUCT_
#define _PARALLEL_HSR_
!=======================================================================
!
!  SOFTWARE NAME : PHASE/0 2023.01
!
!  MODULE: m_Charge_Density
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!      Further modification by T. Yamasaki   Feb. 2004
!      Further modification by T. Yamasaki   Apr/15/2006
!      Further modification by T. Yamasaki   Aug/31/2007
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
!   Revised for the GAMMA point (k=(0,0,0)) by T. Yamasaki, April 2006.
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
#   define __TIMER_COMM_STOP_w_BARRIER(str,a)    call timer_barrier(str) ;   call timer_end(a)
#   define __TIMER_COMM_START(a)       call timer_sta(a)
#   define __TIMER_COMM_STOP(a)        call timer_end(a)
#else
#   define __TIMER_COMM_START_w_BARRIER(str,a)
#   define __TIMER_COMM_STOP_w_BARRIER(str,a)
#   define __TIMER_COMM_START(a)
#   define __TIMER_COMM_STOP(a)
#endif

module m_CD_mixing
! $Id: m_CD_mixing.F90 593 2019-06-20 03:47:31Z jkoga $
  use m_Const_Parameters,    only : DP, OFF, BOHR,NO,YES,ANTIFERRO &
       &                          , ANEW,RENEW,ON, SIMPLE,BROYD1,BROYD2,DFP,PULAY
  use m_IterationNumbers,    only : iteration,iteration_for_cmix
  use m_Parallelization,     only : m_Parallel_init_mpi_urec_hsr, MPI_CommGroup &
       &                          , ista_kngp,iend_kngp,is_kngp,ie_kngp,np_kngp,mp_kngp &
       &                          , npes,mype,ierr &
       &                          , is_kgpm,ie_kgpm,ista_kgpm,iend_kgpm,mp_kgpm &
       &                          , nis_fftp, nie_fftp, myrank_g, nrank_g, ista_atm, iend_atm &
       &                          , ista_urec_hsr,iend_urec_hsr, ista_and_iend_urec_hsr_set
  use m_Timing,              only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters,  only : nspin,ipri,c_precon &
       &                          , amix,bmix,hownew,nbxmix,istrbr &
       &                          , kimg,af,neg,ipripulay,iprichargemixing &
       &                          , sw_recomposing, spin_density_mixfactor &
       &                          , amin, sw_precon_diff, sw_metric_diff,metric_ratio &
       &                          , sw_force_simple_mixing,printable, sw_control_stepsize, max_stepsize &
       &                          , m_CtrlP_set_rmx, ommix_factor
  use m_Crystal_Structure,   only : univol,nopr
  use m_PlaneWaveBasisSet,   only : kg,kgp,gr_l,kgpm
  use m_Charge_Density,      only : chgq_l, chgqo_l ,symmtrz_of_ff
  use m_Charge_Density,      only : charge_average
  use m_Charge_Density,      only : work ! === DEBUG by tkato 2011/09/09 ===

! === Added by tkato 2011/11/09 ================================================
  use m_Control_Parameters,  only : sw_mix_bothspins_sametime &
                                  , sw_recomposing_hsr, sw_force_simple_mixing_hsr &
                                  , num_proj_elems, proj_group, proj_attribute, num_projectors &
                                  , max_projs
  use m_Ionic_System,        only : ityp, natm,iproj_group
  use m_PseudoPotential,     only : ilmt, nlmt
  use m_Charge_Density,      only : hsr, hsro
!===============================================================================

! =================================== added by K. Tagami ============== 11.0
  use m_Control_Parameters,  only : ndim_magmom, noncol, sw_mix_imaginary_hardpart
  use m_Charge_Density,      only : hsi, hsio
  use m_PseudoPotential,     only : ipaw,ia2ia_symmtry_op_inv
! ===================================================================== 11.0

  use m_Control_Parameters,  only : sw_gradient_simplex, alpha_pulay, alpha_pulay_damp, alpha_pulay_org, alpha_pulay_damp_thres

  use m_Control_Parameters,  only : precon_mode  ! ===== KT_add ===== 13.0U3
  use m_Control_Parameters,  only : sw_mix_charge_hardpart, sw_mix_charge_with_ekindens ! ==== KT_Add === 2014/09/16

  use m_Charge_Density,      only : m_CD_symmtrz_of_ff_noncl_C

  use m_Orbital_Population, only : om, omold, ommix, om_aimag, omold_aimag, ommix_aimag
!!$#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
!!$  use m_Files, only : nfout
!!$#endif
  use mpi

  implicit none

  real(kind=DP), allocatable, dimension(:):: rmxtrc ! d(nspin_m)
  real(kind=DP), allocatable, dimension(:,:,:) :: chgqstore_l, chgqostore_l

  real(kind=DP),private,pointer, dimension(:,:,:) ::       rho,rhoo  ! MPI
  !         rho => chgq_l, rhoo => chgqo_l ( when kgpm == kgp)
  !         rho and rhoo is projected from chgq_l and chgqo_l, respectively ( otherwise)

  real(kind=DP),private,allocatable,dimension(:,:,:,:) :: rhoj,Frhoj
  real(kind=DP),private,allocatable,dimension(:,:,:)   :: rhojo,Frhojo
  real(kind=DP),private,allocatable,dimension(:,:,:)   :: rhoj_hsr,Frhoj_hsr
  real(kind=DP),private,allocatable,dimension(:,:)   :: rhoj_hsro,Frhoj_hsro

! === DEBUG by tkato 2011/09/09 ================================================
! real(kind=DP),private, pointer, dimension(:,:)        :: work
! ==============================================================================

!!$  integer, private,parameter   :: n_ratio_q1 = 20
  integer, private             :: n_ratio_q1 = 20
  real(DP),private,parameter   :: q0_default = 1.5*BOHR
  real(DP),private             :: amix_cprec, bmix_cprec
  logical, private             :: param_cprecon_decided = .false.
  real(DP),private             :: q0 = -1.d0, q1 = -1.d0
  real(DP),private,save        :: fg2, frg2

! --> T. Yamasaki, 03rd Aug. 2009
  real(DP),private,pointer,dimension(:,:)     :: c_p !d(ista_kngp:iend_kngp,nspin/(af+1))
  real(DP),private,pointer,dimension(:,:)     :: c_pm !d(ista_kgpm:iend_kgpm,nspin/(af+1))
! <--

  ! -- For Broyden and DFP mixing method --
  integer, private,parameter                :: iU = 1, iVD = 1, iW = 2, iY = 2, iV = 2
  integer, private                          :: nspin_m
  real(DP),private,allocatable,dimension(:) :: f_p !d(ista_kgpm:iend_kgpm)
  real(DP),private,allocatable,dimension(:,:,:) :: d0_l,u_l,v_l,w_l,dout,dd_l
  real(DP),private,allocatable,dimension(:,:,:,:) :: d0_l_h
  real(DP),private,pointer,dimension(:,:,:)         :: F_l
  real(DP),private,allocatable,target,dimension(:,:,:) :: din
  !                                             d(kgpm,kimg,nspin_m)
  real(DP),private,allocatable,dimension(:,:,:)     :: dF_l
  !     dF_l(deltaF):= \Delta \cal F^{m} = \cal F^{m} - \cal F^{m-1}
  !              = F[\rho^{m}] - F[\rho^{m-1}] - (\rho^{m} - \rho^{m-1})
  real(DP),private,allocatable,target,dimension(:,:,:,:,:) :: urec_l
  real(DP),private,allocatable,dimension(:,:,:)     :: f      !d(nbxmix,nbxmix,nspin)
  real(DP),private,allocatable,dimension(:)         :: g      !d(nbxmix)
  !    f and g are used only when hownew == RENEW
  integer, private,allocatable,dimension(:)         :: ncrspd !d(nbxmix)
  real(DP),private,allocatable,dimension(:,:,:)     :: uuf    !d(nbxmix,nspin,2),
                                                   !only for DFP method
  real(DP),private,allocatable,dimension(:,:)       :: uuf_p
  real(DP),private,allocatable,dimension(:,:)       :: g_p
  real(DP),private,allocatable,dimension(:)         :: prj_wk
#ifdef _CDMIX_USE_POINTER_
! Tsuyoshi Miyazaki tmp
  real(DP),private,pointer,dimension(:,:,:)         :: urec_l_3
  real(DP),private,pointer,dimension(:,:,:)         :: urec_l_3_2
#endif

  logical, save :: is_and_ie_hsr_set = .false.

  logical          :: force_dealloc = .false.
  integer, private :: previous_waymix = 0

  real(kind=DP) ::  rmx_max = 0.95d0

!  include 'mpif.h'
  integer istatus(mpi_status_size)

! ========================== adde by K. Tagami ========================== 5.0
  integer :: nsize_rho_hsr
  integer :: nsize_rho_hsr0
  integer :: nsize_rho_om
  integer, private, allocatable :: imap_hsr(:)    ! d(nsize_rho_hsr)
  integer, private, allocatable :: imap_om(:)    ! d(nsize_rho_hsr)
  real(kind=DP),private,allocatable, dimension(:,:) ::  rho_hsr, rhoo_hsr  ! d(nsize_rho_hsr,nspin)

  real(DP),private,allocatable,dimension(:,:) :: d0_hsr, u_hsr, v_hsr, w_hsr, &
        &                                        dout_hsr, dd_hsr
  real(DP),private,pointer,dimension(:,:)         :: FF_hsr
  real(DP),private,allocatable,target,dimension(:,:) :: din_hsr
  !                                             d( nsize_rho_hsr,nspin_m)
  real(DP),private,allocatable,dimension(:,:)     :: dF_hsr
  real(DP),private,allocatable,target,dimension(:,:,:,:) :: urec_hsr

  real(DP),private,allocatable,dimension(:,:,:) :: d0_hsr_h
  logical,allocatable,dimension(:) :: diag_elem

  logical, save :: first = .true.

  real(kind=DP), allocatable :: hsr_store(:,:,:,:)
  real(kind=DP), allocatable :: hsro_store(:,:,:,:)
!
  real(kind=DP), allocatable :: rho_store(:,:)
  real(kind=DP), allocatable :: rhoo_store(:,:)

! ==============================================================

! ========================== adde by K. Tagami ========================== 11.0
  integer :: nsize_rho_hsr_realpart
  integer :: nsize_rho_om_realpart
!  integer :: sw_mix_imaginary_hardpart = OFF
!  integer :: sw_mix_imaginary_hardpart = ON
! ======================================================================= 11.0

  real(kind=DP) :: step_control_factor = 1.d0

  real(kind=DP), allocatable, dimension(:,:) :: ynorm
  real(kind=DP), allocatable, dimension(:) :: sum12

! --------------------------
! meta-gga ( Kinetic Energy Density )
! --------------------------
  real(kind=DP), allocatable, target:: din_ekinq(:,:,:)
  real(kind=DP), allocatable :: dout_ekinq(:,:,:)
  real(kind=DP), allocatable :: dF_l_ekinq(:,:,:)
  real(kind=DP), allocatable, target :: urec_l_ekinq(:,:,:,:,:)
!
  real(kind=DP), pointer :: c_p_ekinq(:,:)
  real(kind=DP), pointer :: c_pm_ekinq(:,:)
  real(kind=DP), allocatable :: f_p_ekinq(:)
!
  real(kind=DP), allocatable :: d0_l_ekinq(:,:,:)
  real(kind=DP), allocatable :: d0_l_h_ekinq(:,:,:,:)
  real(kind=DP), pointer :: F_l_ekinq(:,:,:)
  real(kind=DP), allocatable :: u_l_ekinq(:,:,:)
  real(kind=DP), allocatable :: v_l_ekinq(:,:,:)
!
  real(kind=DP), pointer :: rho_ekinq(:,:,:),  rhoo_ekinq(:,:,:)

  real(kind=DP), allocatable :: rhoj_ekinq(:,:,:,:),  Frhoj_ekinq(:,:,:,:)
  real(kind=DP), allocatable :: rhojo_ekinq(:,:,:), Frhojo_ekinq(:,:,:)
!
  real(kind=DP), allocatable, dimension(:,:) :: ynorm_ekinq

  real(kind=DP), pointer :: ekinq_l(:,:,:), ekinqo_l(:,:,:)
  real(kind=DP), allocatable, dimension(:,:,:) :: ekinqstore_l, ekinqostore_l
! --------------------------

  integer, public, allocatable :: i2lp(:) ! d(num_projectors)

! ================================ modified by K. Tagami ================ 11.0
!  integer, private :: max2lp ! max. of i2lp
  integer, public :: max2lp ! max. of i2lp
! ======================================================================= 11.0

  integer, private :: l1max ! max. of l1
  integer, private :: nyymax

! --- contained subroutines ---
!   7. m_CD_prepare_precon       <-(ChargeDensity_Mixing)
!  10. m_CD_simple_mixing        <-(ChargeDensity_Mixing)
!  22. precon_4_charge_mix       <-(@10), (@53),(@54),(@55),(@60)
!  23. precon_4_mult             <-(@37),(@45),(48)
!  24. iter_from_reset           <-(@53),(@54),(@55),(@60)
!  25. icrspd_is                 <-(@53),(@54),(@55)
!  26. mult1s                    <-(@49),(@50),(@51),(@53),(@54),(@55),(@60)
!  27. mult1s5                   <-(@53),(@55),(@60)
!  28. mult1s10                  <-(@60)
!  29. subtr_j_th_term           <-(@49),(@50),(@53),(@55)
!  30. store_to_urec2            <-(@53),(@55)
!  31. set_ncrspd_mxiter_etc     <-(@53),(@54),(@55)
!    - rotate_cmix_arrays
!  32. simple_mix1               <-(@53),(@54),(@55),(@60)
!  33. scatter_chg_onto_d        <-(@32),(@40),(@51)
!  34. scatter_cp_onto_cpm       <-(@23),(@40)
!  35. concentrate_d_to_chg      <-(@51),(@55),(@60)
!  36. mix_dealloc_previous      <-(@53),(@54),(@55),(@60)
!  37. mix_broyden_allocate      <-(@51),(@54)
!  38. mix_broyden_deallocate    <-(@36)
!  39. mix_broyden_alloc2        <-(@54)
!  40. alloc_rho_rhoo_and_cpm    <-(@39),(@43),(@46),(@58)
!  41. mix_broyden_dealloc2      <-(@54)
!  42. dealloc_rho_rhoo_and_cpm  <-(@41),(@44),(@47),(@59)
!  43. mix_broyden_alloc3        <-(@53)
!  44. mix_broyden_dealloc3      <-(@53)
!  45. mix_DFP_allocate          <-(@55)
!  46. mix_DFP_alloc2            <-(@55)
!  47. mix_DFP_dealloc2          <-(@55)
!  -x. devide_v_with_vdF
!  49. renew_u_br                <-(@53),(@54)
!  50. renew_d_br                <-(@53),(@54)
!  51. renew_d_last_br           <-(@53),(@54)
!  52. simple_mix_large_Gc       <-(@51),(@55),(@60)
!  53. m_CD_mix_broyden1         <-(ChargeDensity_Mixing)
!    - dF_F_d0_u_v_and_dd - renew_v
!  54. m_CD_mix_broyden2         <-(ChargeDensity_Mixing)
!    - dF_F_d0_u_and_v
!  55. m_CD_mix_DFP              <-(ChargeDensity_Mixing)
!    - dF_F_d0_u_and_w - renew_w - renew_d - renew_d_last
!  56. mix_pulay_allocate        <-(@60)
!  57. mix_pulay_deallocate      <-(@36)
!  58. mix_pulay_alloc2          <-(@60)
!  59. mix_pulay_dealloc2        <-(@60)
!  60. m_CD_mix_pulay            <-(ChargeDensity_Mixing)
!    - mix_pulay_alloc3 - mix_pulay_dealloc3
!    - Resid_and_dd_into_urec - Ri_dot_Rj - Rj_dot_d
!    - get_finv -get_matrix -renew_d_using_g

contains

  subroutine m_CD_mixing_write_DEFINITION(nfout)
    ! coded by T. Yamasaki, 2023/07/08
    integer, intent(in) :: nfout
#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    if(ipri>=1) then
#endif
       write(nfout,'(" !!")')
       write(nfout,'(" !! <<m_CD_mixing_write_DEFINITION>>")')
#ifdef _DUPLICATION_HSR_DOTPRODUCT_
       write(nfout,'(" !! Compiler Defintion in (m_CD_mixing.F90) is _DUPLICATION_HSR_DOTPRODUCT_ ,", &
            & "namely asis in HSR related dotproduction")')
#else
#ifdef _PARALLEL_HSR_
       write(nfout,'(" !! Compiler Defintion in (m_CD_mixing.F90) is _PARALLEL_HSR__ ,",&
            & "namely parallelized HSR related dotproduct")')
#else
       write(nfout,'(" !! Compiler Defintion in (m_CD_mixing.F90) is nothing,  namely mpi_bcast after HSR related dotproduct")')
#endif
#endif
       write(nfout,'(" !!")')
#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    endif
#endif
  end subroutine m_CD_mixing_write_DEFINITION
  
  subroutine alloc_chgqstore_recompose_chgq(rmxt,rmxtrc)
    real(kind=DP),intent(in) :: rmxt
    real(kind=DP),intent(out),dimension(nspin_m) :: rmxtrc
                                                  __TIMER_SUB_START(1104)
    allocate(chgqstore_l(ista_kngp:iend_kngp,kimg,nspin))
    allocate(chgqostore_l(ista_kngp:iend_kngp,kimg,nspin))
                                                  __TIMER_DO_START(1149)
    chgqstore_l = chgq_l
    chgqostore_l = chgqo_l
    chgq_l(:,:,1)  = chgqstore_l(:,:,1)  + chgqstore_l(:,:,2)
    chgq_l(:,:,2)  = chgqstore_l(:,:,1)  - chgqstore_l(:,:,2)
    chgqo_l(:,:,1) = chgqostore_l(:,:,1) + chgqostore_l(:,:,2)
    chgqo_l(:,:,2) = chgqostore_l(:,:,1) - chgqostore_l(:,:,2)
                                                  __TIMER_DO_STOP(1149)
    rmxtrc(1) = rmxt
    rmxtrc(2) = rmxt*spin_density_mixfactor
                                                  __TIMER_SUB_STOP(1104)
  end subroutine alloc_chgqstore_recompose_chgq

  subroutine compose_chgq_dealloc_chgqstore
                                                  __TIMER_SUB_START(1106)
                                                  __TIMER_DO_START(1153)
    chgqstore_l = chgq_l
    chgq_l(:,:,1) = 0.5d0*(chgqstore_l(:,:,1) + chgqstore_l(:,:,2))
    chgq_l(:,:,2) = 0.5d0*(chgqstore_l(:,:,1) - chgqstore_l(:,:,2))
    chgqo_l = chgqostore_l
                                                  __TIMER_DO_STOP(1153)
    deallocate(chgqostore_l, chgqstore_l)
                                                  __TIMER_SUB_STOP(1106)
  end subroutine compose_chgq_dealloc_chgqstore

  subroutine m_CD_prepare_precon(nfout,rmxt)
    integer, intent(in)      :: nfout
    real(kind=DP), intent(in):: rmxt
    real(kind=DP) :: G_longest, G_shortest, x, gg,gg_mpi
    integer       :: i, n
    real(kind=DP) :: G_longest_mpi, G_shortest_mpi
    integer       :: ist !mpi

    if(param_cprecon_decided) return
                                                  __TIMER_SUB_START(1102)

    if(iprichargemixing >= 2) write(nfout,*) ' << d_para_cprec >>'

!mpi    G_longest = maxval(gr_l)
!mpi    G_shortest = minval(gr_l(2:kgp))
    G_longest_mpi = maxval(gr_l(ista_kngp:iend_kngp))

    ist = ista_kngp
    if(ist == 1) ist = 2

    G_shortest_mpi = minval(gr_l(ist:iend_kngp))
    call mpi_allreduce(G_longest_mpi,G_longest,1 &
                   &  ,mpi_double_precision,mpi_max,MPI_CommGroup,ierr)
    call mpi_allreduce(G_shortest_mpi,G_shortest,1 &
                   &  ,mpi_double_precision,mpi_min,MPI_CommGroup,ierr)
    if(iprichargemixing >= 2) then
       write(nfout,*) ' G_longest = ', G_longest
       write(nfout,*) ' G_shortest = ', G_shortest
    end if

    amix_cprec = amix; bmix_cprec = bmix
    if(amix_cprec < 0.d0) then
       amix_cprec = rmxt
       if(iprichargemixing >= 2) write(nfout,*) ' amix_cprec = ', amix_cprec
    end if

! ====== KT_mod ======== 13.0U3
!    if(bmix_cprec < 0.d0) then
!       q0 = q0_default
!    else
!       q0 = bmix_cprec*G_shortest
!    end if
!    q0 = q0*q0
!
    select case (precon_mode)
    case (1)                               ! default
       if(bmix_cprec < 0.d0) then
          q0 = q0_default
       else
          q0 = bmix_cprec*G_shortest        ! Kerker
       end if
       q0 = q0*q0
    case (2)                                ! Kresse
       q0 = bmix_cprec**2
    end select
! ========================= 13.0U3

    if(iprichargemixing >= 2) then
       write(nfout,*) ' bmix_cprec = ', bmix_cprec
       write(nfout,*) ' q0   = ', q0
    end if

    if (metric_ratio>0) n_ratio_q1 = metric_ratio

    x = G_longest**2  - n_ratio_q1*G_shortest**2
    if(x < 0.d0) then
       n = (G_longest/G_shortest)**2 - 1.d0
       if(n <= 0) n = 1
       x = G_longest**2 - n * G_shortest**2
    else
       n = n_ratio_q1
    end if
    if(iprichargemixing >= 2) write(nfout,*) ' n = ', n
    q1 = (n-1) * G_shortest**2 * G_longest**2/x
    q1 = dsqrt(q1)
    if(iprichargemixing >= 2) write(nfout,*) ' q1 = ', q1

    if(nspin == 2 .or. noncol ) then
       i = 2
!!$       gg = gr_l(i)**2
       gg = 0.d0
       if(ista_kngp <= i .and. i <= iend_kngp) gg = gr_l(i)**2
                                                 __TIMER_COMM_START_w_BARRIER(mpi_chg_world,1147)
       if(npes > 1) then
          call mpi_allreduce(gg,gg_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
          gg = gg_mpi
       end if
                                                 __TIMER_COMM_STOP_w_BARRIER(mpi_chg_world,1147)

       if(iprichargemixing >= 2) write(nfout,'(" !! gg = ",d20.8)') gg

       if (amin<=0) then
          amin = 0.d0
       endif
       fg2 = max(gg/(gg+q0),amin)
       frg2 = 1.0d0 + q1**2/gg
       if(iprichargemixing >= 2) then
          write(nfout,*) ' ! amix_cprec = ', amix_cprec, ' fg2 = ', fg2
          write(nfout,*) ' ! amix_cprec*g_S**2/(g_S**2+q0) = ', amix_cprec*fg2
       end if
    else
       fg2 = 0.d0
       frg2 = 0.d0
    endif

    param_cprecon_decided = .true.
                                                  __TIMER_SUB_STOP(1102)
  end subroutine m_CD_prepare_precon

  subroutine m_CD_force_dealloc()
    force_dealloc = .true.
  end subroutine m_CD_force_dealloc

  subroutine m_CD_simple_mixing(nfout,rmxt)
    integer ,intent(in)      :: nfout
    real(kind=DP),intent(in) :: rmxt

    integer       :: is, k
    real(kind=DP) :: rmxtt
    integer       :: id_sname = -1
                                                  __TIMER_SUB_START(1103)
    call tstatc0_begin('m_CD_simple_mixing ',id_sname,1)



    if ( noncol ) then             ! === modified by K. Tagami === 11.0
       nspin_m = ndim_magmom
    else
       nspin_m  = nspin/(af+1)
    endif


    if(previous_waymix /= SIMPLE.or.force_dealloc) then
       call mix_dealloc_previous()
       call mix_dealloc_previous_hsr()  ! ---  ktDEBUG ---------------- 20121030
       force_dealloc = .false.
    end if

    allocate(rmxtrc(nspin_m))
    if ( noncol ) then                  ! === modified by K. Tagami === 11.0
       rmxtrc = rmxt
       rmxtrc(2:nspin_m) = min( rmxt *spin_density_mixfactor, rmx_max )
    else
       if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
          call alloc_chgqstore_recompose_chgq(rmxt,rmxtrc) ! --> chgq_l, chgqo_l, rmxtrc
       else
          rmxtrc = rmxt
       endif
    end if

    if(ipri >= 2) write(nfout,'(" rmxt = ",d20.8)') rmxt

    allocate(c_p(ista_kngp:iend_kngp,nspin_m)); c_p = 0.d0 ! =================== 11.0

! ================================ modified by K. Tagami =============== 11.0
!!    call precon_4_charge_mix(rmxtrc,c_p)
!
    if ( noncol ) then
       call precon_4_charge_mix_noncl(rmxtrc,c_p)
    else
       call precon_4_charge_mix(rmxtrc,c_p)
    endif
! ======================================================================= 11.0

                                                  __TIMER_DO_START(1148)
    do is = 1, ndim_magmom, af+1  !!  do is = 1, nspin, af+1  ! === modified by K. Tagami === 11.0
       do k = 1, kimg
          chgq_l(:,k,is) = c_p(:,is)*chgq_l(:,k,is) + (1.0d0-c_p(:,is))*chgqo_l(:,k,is)
       end do
    end do
                                                  __TIMER_DO_STOP(1148)
    deallocate(c_p)
    if ( .not. noncol ) then     ! === modified by K. Tagami ========= 11.0
       if (sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
          call compose_chgq_dealloc_chgqstore()
       end if
    endif
    deallocate(rmxtrc)

    if(af /= 0)  then
       allocate(work(kgp,kimg))
       call charge_average(ANTIFERRO,chgq_l)
       deallocate(work)
    endif

    previous_waymix = SIMPLE
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(1103)
  end subroutine m_CD_simple_mixing

  subroutine precon_4_charge_mix(pmix,c_p)
    real(DP),intent(in),dimension(nspin_m)                       :: pmix
    real(DP),intent(out), dimension(ista_kngp:iend_kngp,nspin_m) :: c_p
    integer              :: i
    real(DP)             :: gg, agg, tmp
    integer              :: ist  !for mpi
    integer              :: is
                                                  __TIMER_SUB_START(1105)

    if(iprichargemixing >= 2) write(6,'("! pmix(precon_4_charge_mix) = ",f8.4)') pmix
    if(c_precon) then
       ist = ista_kngp
       if(ist == 1) ist = 2
       if(nspin_m == 2) then
                                                  __TIMER_DO_START(1150)
          do i = ist, iend_kngp  !for mpi
             gg = gr_l(i)*gr_l(i)
             tmp = max(gg/(gg+q0),amin)
!!$             agg = amix_cprec*gg/(gg+q0)
             agg = amix_cprec*tmp
             c_p(i,1) = agg * pmix(1)
             if (sw_recomposing==ON .and. sw_precon_diff==NO)then
                c_p(i,2) = amix_cprec * pmix(2)
             else
                c_p(i,2) = agg * pmix(2)
             endif
          enddo
                                                  __TIMER_DO_STOP(1150)
          if(mype == 0) then
             c_p(1,1) = amix_cprec*fg2 * pmix(1)
             if (sw_recomposing==ON .and. sw_precon_diff==NO)then
                c_p(1,2) = amix_cprec * pmix(2)
             else
                c_p(1,2) = amix_cprec*fg2 * pmix(2)
             endif
          end if
       else
                                                  __TIMER_DO_START(1151)
          do i = ist, iend_kngp  !for mpi
             gg = gr_l(i)*gr_l(i)
             tmp = max(gg/(gg+q0),amin)
!!$             c_p(i,1) = amix_cprec*gg/(gg+q0) * pmix(1)
             c_p(i,1) = amix_cprec* tmp * pmix(1)
          end do
                                                  __TIMER_DO_STOP(1151)
          if(mype == 0) c_p(1,1) = amix_cprec*fg2 * pmix(1)
       end if
    else
                                                  __TIMER_DO_START(1152)
       do is = 1, nspin, af+1
          c_p(:,is) = pmix(is)
       end do
                                                  __TIMER_DO_STOP(1152)
    endif
                                                  __TIMER_SUB_STOP(1105)
  end subroutine precon_4_charge_mix

! ===================================== added by K. Tagami ================ 11.0
  subroutine precon_4_charge_mix_noncl(pmix,c_p)
    real(DP),intent(in),dimension(nspin_m)                       :: pmix
    real(DP),intent(out), dimension(ista_kngp:iend_kngp,nspin_m) :: c_p
    integer              :: i
    real(DP)             :: gg, agg, tmp
    integer              :: ist  !for mpi
    integer              :: is

    if(iprichargemixing >= 2) then
       write(6,'("! pmix(precon_4_charge_mix_noncl) = ",f8.4)') pmix
    end if

!
!    write(*,*) " amix = ", amix_cprec
!    write(*,*) 'fg2 ', fg2
!    write(*,*) 'pmix ', pmix
!
    if(c_precon) then
       ist = ista_kngp
       if(ist == 1) ist = 2

       do i = ist, iend_kngp  !for mpi
          gg = gr_l(i)*gr_l(i)
          tmp = max(gg/(gg+q0),amin)
!!$             agg = amix_cprec*gg/(gg+q0)
          agg = amix_cprec*tmp

! === 2015/09/29
!!!        c_p(i,:) = agg * pmix(:)

          c_p(i,1) = agg * pmix(1)
          c_p(i,2:nspin_m) = amix_cprec * pmix(2:nspin_m)
! === 2015/09/29

       end do
       if (mype == 0) then
! === 2015/09/29
!          c_p(1,:) = amix_cprec*fg2 * pmix(:)

          c_p(1,1) = amix_cprec*fg2 * pmix(1)
          c_p(1,2:nspin_m) = amix_cprec * pmix(2:nspin_m)
! === 2015/09/29
       end if
    else
       do is = 1, ndim_magmom
          c_p(:,is) = pmix(is)
       end do
    endif

  end subroutine precon_4_charge_mix_noncl
! ============================================================== 11.0

  subroutine precon_4_mult(f_q)
    real(DP),intent(out), dimension(ista_kgpm:iend_kgpm)   :: f_q
    integer  :: i, ist      !mpi
    real(kind=DP), pointer, dimension(:) :: gr_l_m

    f_q = 0.d0
    if(c_precon ) then
       ist  = ista_kgpm
       if(ist == 1) ist = 2
       if(kgp == kgpm .or. npes == 1) then
          do i = ist, iend_kgpm  !for mpi
             f_q(i) = 1.0d0 + (q1/gr_l(i))**2
          end do
       else
          allocate(gr_l_m(ista_kngp:iend_kngp))
! ============================== by K. Tagami =============
        gr_l_m = 0.0d0
! ========================================================
          call scatter_cp_onto_cpm(gr_l,gr_l_m)
          do i = ist, iend_kgpm  !for mpi
             f_q(i) = 1.0d0 + (q1/gr_l_m(i))**2
          end do
          deallocate(gr_l_m)
       end if
       if(mype==0) f_q(1) = frg2
!mpi       f_q(2:kgpm) = 1 + (q1/gr_l(2:kgpm))**2
!mpi       f_q(1)       = frg2
    else
       f_q          = 1.d0
    end if

  end subroutine precon_4_mult

  function iter_from_reset()
    integer             :: n, nbox
    integer             :: iter_from_reset
    if(hownew ==  ANEW) then
       n = (iteration_for_cmix - istrbr - 1)/(nbxmix-1)
       if(n < 0) n = 0
       nbox = iteration_for_cmix - (n*(nbxmix-1) + istrbr +1) + 2
       iter_from_reset = nbox + istrbr - 1
    else
       iter_from_reset = iteration_for_cmix
    endif
  end function iter_from_reset

  function icrspd_is(iter)
    integer, intent(in) :: iter
    integer             :: icrspd_is
    if(iter-istrbr+1 < nbxmix) then
       icrspd_is = ncrspd(iter-istrbr+1)
    else
       icrspd_is = ncrspd(nbxmix)
    endif
  end function icrspd_is

  subroutine mult_urec_hsr(nfout,u_hsr,v_hsr,fdpsum)
    ! Coded by T. Yamasaki, 2023/07/07
    integer :: nfout
    real(DP), intent(in), dimension(nsize_rho_hsr,nspin_m) :: u_hsr, v_hsr
    real(DP), intent(inout), dimension(nspin_m) :: fdpsum
    real(DP), dimension(nspin_m) :: fmult
    integer :: is, k
#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    real(DP) :: time0, time1
#endif

#ifdef _PARALLEL_HSR_
    if(.not.ista_and_iend_urec_hsr_set) &
         & call m_Parallel_init_mpi_urec_hsr(nfout,nsize_rho_hsr) !-> ista_urec_hsr, iend_urec_hsr
#endif
#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    call MPI_Barrier( MPI_CommGroup,ierr)
    time0 = MPI_Wtime()
#endif
    fmult = 0.d0
#ifdef _PARALLEL_HSR_
    do is = 1, ndim_magmom, af+1
       do k = ista_urec_hsr, iend_urec_hsr
          fmult(is) = fmult(is) + u_hsr(k,is) * v_hsr(k,is)
       end do
    end do
    if(npes>=2) then
       call mpi_allreduce(MPI_IN_PLACE, fmult, nspin_m, mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
    end if
#else
#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    write(nfout,'(" not _PARALLEL_HSR_")')
#endif
    do is = 1, ndim_magmom, af+1
       fmult(is) = fmult(is) + sum( u_hsr(:,is)*v_hsr(:,is) )
    end do
    call mpi_bcast(fmult, nspin_m, mpi_double_precision, 0, MPI_CommGroup,ierr)
#endif
#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    call MPI_Barrier( MPI_CommGroup,ierr)
    time1 = MPI_Wtime()
#ifdef _PARALLEL_HSR_
    write(nfout,'(" time in <<mult_urec_hsr>> = ",f20.8, " (mpi_allreduce)")') time1-time0
#else
    write(nfout,'(" time in <<mult_urec_hsr>> = ",f20.8, " (mpi_bcast)")') time1-time0
#endif
#endif
    fdpsum = fdpsum + fmult
  end subroutine mult_urec_hsr

  subroutine mult_urec_hsr5(nfout,u_hsr,mb,muv,j,iuv,v_hsr,fdpsum)
    ! Coded by T. Yamasaki, 2023/07/07
    integer, intent(in) :: nfout,mb,muv,j,iuv
    real(DP), intent(in), dimension(nsize_rho_hsr,nspin_m,mb,muv) :: u_hsr
    real(DP), intent(in), dimension(nsize_rho_hsr,nspin_m) ::        v_hsr
    real(DP), intent(out), dimension(nspin_m) :: fdpsum
    real(DP),              dimension(nspin_m) :: fmult
    integer :: is, k
#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    real(DP) :: time0, time1
#endif

#ifdef _PARALLEL_HSR_
    if(.not.is_and_ie_hsr_set) then
       call m_Parallel_init_mpi_urec_hsr(nfout,nsize_rho_hsr) !-> ista_urec_hsr, iend_urec_hsr
       is_and_ie_hsr_set = .true.
    end if
#endif
#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    call MPI_Barrier( MPI_CommGroup,ierr)
    time0 = MPI_Wtime()
#endif
    fmult = 0.d0
#ifdef _PARALLEL_HSR_
    do is = 1, ndim_magmom, af+1     !      do is = 1, nspin, af+1  ! === modified by K. Tagami === 11.0
       do k = ista_urec_hsr, iend_urec_hsr
          fmult(is) = fmult(is) + u_hsr(k,is,j,iuv) * v_hsr(k,is)
       end do
    end do
    if(npes>=2) then
       call mpi_allreduce(MPI_IN_PLACE, fmult, nspin_m, mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
    end if
#else
#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    write(nfout,'(" not _PARALLEL_HSR_")')
#endif
    do is = 1, ndim_magmom, af+1    !   do is = 1, nspin, af+1  ! === modified by K. Tagami === 11.0
       fmult(is) = fmult(is) + sum( u_hsr(:,is,j,iuv)*v_hsr(:,is) )
    end do
    call mpi_bcast(fmult, nspin_m, mpi_double_precision, 0, MPI_CommGroup,ierr)
#endif
#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    call MPI_Barrier( MPI_CommGroup,ierr)
    time1 = MPI_Wtime()
#ifdef _PARALLEL_HSR_
    write(nfout,'(" time in <<mult_urec_hsr5>> = ",f20.8, " (mpi_allreduce)")') time1-time0
#else
    write(nfout,'(" time in <<mult_urec_hsr5>> = ",f20.8, " (mpi_bcast)")') time1-time0
#endif
#endif
    fdpsum = fdpsum + fmult
  end subroutine mult_urec_hsr5

  subroutine mult1s(u,v,f_q,fmult)
    real(DP),intent(in), dimension(ista_kgpm:iend_kgpm,kimg,nspin_m) :: u,v
    real(DP),intent(in), dimension(ista_kgpm:iend_kgpm):: f_q
    real(DP),intent(out),dimension(nspin_m)            :: fmult

    real(DP) :: p, p_mpi, fac
    integer  :: is,ik,i
                                                  __TIMER_SUB_START(1114)

    fmult = 0.d0

! ================================ modified by K. Tagami ============== 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ===================================================================== 11.0

       p = 0.d0
       fac=1.0d0
                                                  __TIMER_DO_START(1160)
       do ik = 1,kimg
          do i = ista_kgpm, iend_kgpm   ! mpi
! ========================================== modified by K. Tagami ======== 11.0
!             if (is==1 .or. sw_recomposing==OFF .or. sw_metric_diff==ON) fac=f_q(i)
!
             if ( noncol ) then
                fac=f_q(i)
             else
                if (is==1 .or. sw_recomposing==OFF .or. sw_metric_diff==ON) then
                   fac=f_q(i)
                endif
             endif
! ======================================================================== 11.0

!             p = p + f_q(i)*u(i,ik,is)*v(i,ik,is)
             p = p + fac*u(i,ik,is)*v(i,ik,is)
          end do
       end do
                                                  __TIMER_DO_STOP(1160)
       if( npes >= 2) then
          call mpi_allreduce(p,p_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
          p = p_mpi
       end if
       fmult(is) = p*univol
    enddo
                                                  __TIMER_SUB_STOP(1114)
  end subroutine mult1s

  subroutine mult1s_reduce_spin(u,v,f_q,fmult)
    real(DP),intent(in), dimension(ista_kgpm:iend_kgpm,kimg,nspin_m) :: u,v
    real(DP),intent(in), dimension(ista_kgpm:iend_kgpm):: f_q
!!$    real(DP),intent(out),dimension(nspin_m)            :: fmult
    real(DP),intent(out)            :: fmult

    real(DP) :: p, p_mpi, fac
    integer  :: is,ik,i

    fmult = 0.d0
    p = 0.d0

! ================================ modified by K. Tagami ============== 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ===================================================================== 11.0

!!$       p = 0.d0
       fac=1.0d0
       do ik = 1,kimg
          do i = ista_kgpm, iend_kgpm   ! mpi
! ========================================== modified by K. Tagami ======== 11.0
!             if (is==1 .or. sw_recomposing==OFF .or. sw_metric_diff==ON) fac=f_q(i)
!
             if ( noncol ) then
                fac=f_q(i)
             else
                if (is==1 .or. sw_recomposing==OFF .or. sw_metric_diff==ON) then
                   fac=f_q(i)
                endif
             endif
! ======================================================================== 11.0

!             p = p + f_q(i)*u(i,ik,is)*v(i,ik,is)
             p = p + fac*u(i,ik,is)*v(i,ik,is)
          end do
       end do
       if( npes >= 2) then
          call mpi_allreduce(p,p_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
          p = p_mpi
       end if
    enddo
    fmult = p
  end subroutine mult1s_reduce_spin

  subroutine mult1s5(u,mb,muv,j,iuv,v,f_q,fmult)
    integer,intent(in) :: mb,muv,j,iuv
    real(DP),intent(in), dimension(ista_kgpm:iend_kgpm,kimg,nspin_m,mb,muv) :: u
    real(DP),intent(in), dimension(ista_kgpm:iend_kgpm,kimg,nspin_m) :: v
    real(DP),intent(in), dimension(ista_kgpm:iend_kgpm):: f_q
    real(DP),intent(out),dimension(nspin_m)            :: fmult

    real(DP) :: p, p_mpi, fac
    integer  :: is,ik,i
                                                  __TIMER_SUB_START(1115)
    fmult = 0.d0
! ================================ modified by K. Tagami ============== 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ===================================================================== 11.0

       p = 0.d0
       fac=1.0d0
                                                  __TIMER_DO_START(1162)
       do ik = 1,kimg
          do i = ista_kgpm, iend_kgpm   ! mpi
! ========================================== modified by K. Tagami ======== 11.0
!             if (is==1 .or. sw_recomposing==OFF .or. sw_metric_diff==ON) fac=f_q(i)
!
             if ( noncol ) then
                fac=f_q(i)
             else
                if (is==1 .or. sw_recomposing==OFF .or. sw_metric_diff==ON) then
                   fac=f_q(i)
                endif
             end if
! ========================================================================= 11.0

             p = p + fac*u(i,ik,is,j,iuv)*v(i,ik,is)
          end do
       end do
                                                  __TIMER_DO_STOP(1162)
       if( npes >= 2) then
          call mpi_allreduce(p,p_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
          p = p_mpi
       end if
       fmult(is) = p*univol
    enddo
                                                  __TIMER_SUB_STOP(1115)
  end subroutine mult1s5

  subroutine mult1s5_reduce_spin(u,mb,muv,j,iuv,v,f_q,fmult)
    integer,intent(in) :: mb,muv,j,iuv
    real(DP),intent(in), dimension(ista_kgpm:iend_kgpm,kimg,nspin_m,mb,muv) :: u
    real(DP),intent(in), dimension(ista_kgpm:iend_kgpm,kimg,nspin_m) :: v
    real(DP),intent(in), dimension(ista_kgpm:iend_kgpm):: f_q
!!$    real(DP),intent(out),dimension(nspin_m)            :: fmult
    real(DP),intent(out)            :: fmult

    real(DP) :: p, p_mpi, fac
    integer  :: is,ik,i

    fmult = 0.d0
    p = 0.d0

! ================================ modified by K. Tagami ============== 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ===================================================================== 11.0

!!$       p = 0.d0
       fac = 1.0d0
       do ik = 1,kimg
          do i = ista_kgpm, iend_kgpm   ! mpi
! ========================================== modified by K. Tagami ======== 11.0
!             if (is==1 .or. sw_recomposing==OFF .or. sw_metric_diff==ON) fac=f_q(i)
!
             if ( noncol ) then
                fac=f_q(i)
             else
                if (is==1 .or. sw_recomposing==OFF .or. sw_metric_diff==ON) then
                   fac=f_q(i)
                endif
             end if
! ========================================================================= 11.0

             p = p + fac*u(i,ik,is,j,iuv)*v(i,ik,is)
          end do
       end do
    enddo
    if( npes >= 2) then
       call mpi_allreduce(p,p_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       p = p_mpi
    end if
    fmult = p*univol
  end subroutine mult1s5_reduce_spin

  subroutine mult1s10(u,mb,muv,i,iu,v,j,iv,f_q,fmult)
    integer,intent(in) :: mb,muv,i,iu,j,iv
    real(DP),intent(in), dimension(ista_kgpm:iend_kgpm,kimg,nspin_m,mb,muv) :: u,v
    real(DP),intent(in), dimension(ista_kgpm:iend_kgpm):: f_q
    real(DP),intent(out),dimension(nspin_m)            :: fmult

    real(DP) :: p, p_mpi, fac
    integer  :: is,ik,ig
                                                  __TIMER_SUB_START(1137)
    fmult = 0.d0

! ====================================== modified by K. Tagami =========== 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ======================================================================== 11.0

       p = 0.d0
       fac = 1.0d0
                                                  __TIMER_DO_START(1186)
       do ik = 1,kimg
          do ig = ista_kgpm, iend_kgpm   ! mpi
! ====================================== modified by K. Tagami ============= 11.0
!!             if (is==1 .or. sw_recomposing==OFF .or. sw_metric_diff==ON) fac=f_q(ig)
!
             if ( noncol ) then
                fac=f_q(ig)
             else
                if (is==1 .or. sw_recomposing==OFF .or. sw_metric_diff==ON) then
                   fac=f_q(ig)
                endif
             end if
! ========================================================================== 11.0

             p = p + fac*u(ig,ik,is,i,iu)*v(ig,ik,is,j,iv)
          end do
       end do
                                                  __TIMER_DO_STOP(1186)
       if( npes >= 2) then
          call mpi_allreduce(p,p_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
          p = p_mpi
       end if
       fmult(is) = p*univol
    enddo
                                                  __TIMER_SUB_STOP(1137)
  end subroutine mult1s10

  subroutine mult1s10_reduce_spin(u,mb,muv,i,iu,v,j,iv,f_q,fmult)
    integer,intent(in) :: mb,muv,i,iu,j,iv
    real(DP),intent(in), dimension(ista_kgpm:iend_kgpm,kimg,nspin_m,mb,muv) :: u,v
    real(DP),intent(in), dimension(ista_kgpm:iend_kgpm):: f_q
    real(DP),intent(out)            :: fmult

    real(DP) :: p, p_mpi, fac
    integer  :: is,ik,ig

    fmult = 0.d0
    p = 0.d0

! ====================================== modified by K. Tagami =========== 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ======================================================================== 11.0

!!$       p = 0.d0
       fac = 1.0d0
       do ik = 1,kimg
          do ig = ista_kgpm, iend_kgpm   ! mpi

! ====================================== modified by K. Tagami ============= 11.0
!!             if (is==1 .or. sw_recomposing==OFF .or. sw_metric_diff==ON) fac=f_q(ig)
!
             if ( noncol ) then
                fac=f_q(ig)
             else
                if (is==1 .or. sw_recomposing==OFF .or. sw_metric_diff==ON) then
                   fac=f_q(ig)
                endif
             end if
! ========================================================================== 11.0

             p = p + fac*u(ig,ik,is,i,iu)*v(ig,ik,is,j,iv)
          end do
       end do
    enddo
    if( npes >= 2) then
       call mpi_allreduce(p,p_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       p = p_mpi
    end if
    fmult = p*univol
  end subroutine mult1s10_reduce_spin

  subroutine subtr_j_th_term(f,iuv,j,urec_l,um)
! === DEBUG by tkato 2011/11/24 ================================================
!   real(DP), intent(in),   dimension(nspin) :: f
    real(DP), intent(in),   dimension(nspin_m) :: f
! ==============================================================================
    integer,  intent(in)    :: iuv,j
    real(DP), intent(in)    :: urec_l(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix,2)
    real(DP), intent(inout) :: um(ista_kgpm:iend_kgpm,kimg,nspin_m)

    integer :: is, ik, i, istart
                                                  __TIMER_SUB_START(1116)

    istart = ista_kgpm
    if(istart == 1) istart = 2
                                                  __TIMER_DO_START(1164)

! =============================== modified by K. Tagami =================== 11.0
!!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ========================================================================== 11.0

       do ik = 1, kimg
          do i = istart, iend_kgpm   ! mpi
             um(i,ik,is) = um(i,ik,is) - f(is)*urec_l(i,ik,is,j,iuv)
          end do
       end do
    end do
                                                  __TIMER_DO_STOP(1164)
                                                  __TIMER_SUB_STOP(1116)
  end subroutine subtr_j_th_term

  subroutine store_to_urec2(v,f,j,iuv)
    real(DP), intent(in)                     :: v(ista_kgpm:iend_kgpm,kimg,nspin_m)
! === DEBUG by tkato 2011/11/24 ================================================
!   real(DP), intent(in),  dimension(nspin)  :: f
    real(DP), intent(in),  dimension(nspin_m)  :: f
! ==============================================================================
    integer , intent(in)                     :: j,iuv

    real(DP)          :: dv
    integer           :: is
                                                  __TIMER_SUB_START(1119)
                                                  __TIMER_DO_START(1165)
! ====================================== modified by K. Tagami ============= 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ========================================================================== 11.0

       dv = 1.d0/f(is)
       urec_l(:,:,is,j,iuv) = v(:,:,is)*dv
    end do
                                                  __TIMER_DO_STOP(1165)
                                                  __TIMER_SUB_STOP(1119)
  end subroutine store_to_urec2

  subroutine set_ncrspd_mxiter_etc(iter,iuv,mxiter)
    integer, intent(in)  :: iuv,iter
    integer, intent(out) :: mxiter
                                                  __TIMER_SUB_START(1111)
    if(hownew == RENEW) then
       if((iter-istrbr+1) >= 3) then
          if((iter-istrbr+1) > nbxmix) then   ! When the box overflows
             call rotate_cmix_arrays          !-(contained here) ->mxiter,ncrspd,urec_l,f,g
          else
             mxiter = (iter-istrbr+1) - 1
             ncrspd(iter-istrbr+1) = iter-istrbr+1
          endif
       else
          mxiter = (iter-istrbr+1) - 1
          ncrspd(1) = 1
          ncrspd(2) = 2
       endif
    else ! if(hownew == ANEW)
       mxiter = (iter-istrbr+1) - 1
       ncrspd(iter-istrbr+1) = iter-istrbr+1
    endif
!!$#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
!!$    write(nfout,'(" mxiter = ",i5," <<set_ncrspd_mxiter_etc>>")') mxiter
!!$    write(nfout,'(" ncrspd = ",15i5)') ncrspd(1:mxiter)
!!$    call flush(nfout)
!!$#endif
                                                  __TIMER_SUB_STOP(1111)
  contains
    subroutine rotate_cmix_arrays
      integer :: is,j,i,icr,jcr,iwork
                                                  __TIMER_SUB_START(1112)
                                                  __TIMER_DO_START(1158)
! ======================================= modified by K. Tagami ========= 11.0
!      do is = 1, nspin, af+1
      do is = 1, ndim_magmom, af+1
! ======================================================================= 11.0

         do j = 3, nbxmix
            icr = ncrspd(2)
            jcr = ncrspd(j)
            g(j) = f(icr,jcr,is)
            do i = 3, j-1
               icr = ncrspd(i)
               g(j) = g(j) - f(icr,jcr,is)*g(i)
            enddo
            icr = ncrspd(2)
            urec_l(:,:,is,jcr,iuv) &
                 & = urec_l(:,:,is,jcr,iuv) + g(j)*urec_l(:,:,is,icr,iuv)
         enddo
      enddo
                                                  __TIMER_DO_STOP(1158)
      mxiter = nbxmix-1
      iwork = ncrspd(2)
                                                  __TIMER_DO_START(1159)
      do i = 2, mxiter
         ncrspd(i)= ncrspd(i+1)
      end do
                                                  __TIMER_DO_STOP(1159)
      ncrspd(mxiter+1) = iwork
                                                  __TIMER_SUB_STOP(1112)
    end subroutine rotate_cmix_arrays
  end subroutine set_ncrspd_mxiter_etc

  subroutine simple_mix2(p)
    real(DP), intent(in), dimension(ista_kngp:iend_kngp,nspin_m) :: p
    integer  :: is,k

! ==================================== added by K. Tagami ============ 11.0
    if ( noncol ) return
! ===================================================================11.0

    if(nspin<2 .or. af==1) return

    if(kgpm == kgp .or. npes == 1) then
!!$       write(6,'(" ! kgpm == kgp")')
       do is = 2, 2
          din (ista_kgpm:iend_kgpm,:,is) = chgqo_l(ista_kgpm:iend_kgpm,:,is)
          dout(ista_kgpm:iend_kgpm,:,is) = chgq_l (ista_kgpm:iend_kgpm,:,is)
       end do
    else
       call scatter_chg_onto_d(chgqo_l,din )  ! -(m_C.D.)
       call scatter_chg_onto_d(chgq_l, dout)  ! -(m_C.D.)
    end if

    do is = 2,2
       do k = 1, kimg
          chgq_l(:,k,is) = p(:,is)*chgq_l(:,k,is) + (1.0d0-p(:,is))*chgqo_l(:,k,is)
       end do
    end do
  end subroutine simple_mix2

  subroutine simple_mix1(p)
    real(DP), intent(in), dimension(ista_kngp:iend_kngp,nspin_m) :: p
    integer  :: is,k
                                                  __TIMER_SUB_START(1108)

    if(kgpm == kgp .or. npes == 1) then
                                                  __TIMER_DO_START(1154)
! ========================= modified by K. Tagami ==================== 11.0
!       do is = 1, nspin, af+1
       do is = 1, ndim_magmom, af+1
! ==================================================================== 11.0

          din (ista_kgpm:iend_kgpm,:,is) = chgqo_l(ista_kgpm:iend_kgpm,:,is)
          dout(ista_kgpm:iend_kgpm,:,is) = chgq_l (ista_kgpm:iend_kgpm,:,is)
       end do
                                                  __TIMER_DO_STOP(1154)
    else
       call scatter_chg_onto_d(chgqo_l,din )  ! -(m_C.D.)
       call scatter_chg_onto_d(chgq_l, dout)  ! -(m_C.D.)
    end if
                                                  __TIMER_DO_START(1155)
! ========================= modified by K. Tagami ==================== 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ==================================================================== 11.0
       do k = 1, kimg
          chgq_l(:,k,is) = p(:,is)*chgq_l(:,k,is) + (1.0d0-p(:,is))*chgqo_l(:,k,is)
       end do
    end do
                                                  __TIMER_DO_STOP(1155)
                                                  __TIMER_SUB_STOP(1108)
  end subroutine simple_mix1

  subroutine scatter_chg_onto_d(c,d)
! ================================== modified by K. Tagami ============= 11.0
!    real(DP),intent(in), dimension(ista_kngp:iend_kngp,kimg,nspin) :: c
!    real(DP),intent(out),dimension(ista_kgpm:iend_kgpm,kimg,nspin) :: d
!
    real(DP),intent(in), dimension(ista_kngp:iend_kngp,kimg,ndim_magmom) :: c
    real(DP),intent(out),dimension(ista_kgpm:iend_kgpm,kimg,ndim_magmom) :: d
! ======================================================================= 11.0

    integer :: ip,is,istart,iend,nelmnt,i,ik,ipbase
                                                  __TIMER_SUB_START(1121)
    nelmnt = mp_kngp*kimg*nspin_m
    prj_wk = 0.d0
    do ip = 0, npes - 1
       ! (1)  coping input data onto a work array, and broadcasting
       if(is_kngp(ip) > kgpm) exit
                                                  __TIMER_DO_START(1168)
       if(ip == mype) then

! =============================== modified by K. Tagami ============== 11.0
!          do is = 1, nspin, af+1
          do is = 1, ndim_magmom, af+1
! ==================================================================== 11.0

             do ik = 1, kimg
                ipbase = mp_kngp*(ik-1) + mp_kngp*kimg*(is-1)
                do i = 1, iend_kngp-ista_kngp+1
                   prj_wk(i + ipbase) = c(ista_kngp-1+i,ik,is)
!!$                   prj_wk(i,ik,is) = c(ista_kngp-1+i,ik,is)
                end do
             end do
          end do
       end if
                                                  __TIMER_DO_STOP(1168)
       call mpi_bcast(prj_wk,nelmnt,mpi_double_precision,ip,MPI_CommGroup,ierr)

       ! (2) projection
       istart = ista_kgpm; if(istart < is_kngp(ip)) istart = is_kngp(ip)
       iend   = iend_kgpm; if(iend   > ie_kngp(ip)) iend   = ie_kngp(ip)
       if(iend < istart) cycle
                                                  __TIMER_DO_START(1170)
! ====================================== modified by K. Tagami ============ 11.0
!       do is = 1, nspin, af+1
       do is = 1, ndim_magmom, af+1
! =========================================================================== 11.0

          do ik = 1, kimg
             ipbase = mp_kngp*(ik-1) + mp_kngp*kimg*(is-1)
             do i = istart, iend
                d(i,ik,is) = prj_wk(i + ipbase)
!!$                d(i,ik,is) = prj_wk(i - is_kngp(ip)+1,ik,is)
             end do
          end do
!!$          d(istart:iend,:,is) = prj_wk(istart-is_kngp(ip)+1:iend-is_kngp(ip)+1,:,is)
       end do
                                                  __TIMER_DO_STOP(1170)
    end do
                                                  __TIMER_SUB_STOP(1121)
  end subroutine scatter_chg_onto_d

  subroutine scatter_cp_onto_cpm(cp,cpm)
    real(DP),intent(in), dimension(ista_kngp:iend_kngp,nspin_m) :: cp
    real(DP),intent(out),dimension(ista_kgpm:iend_kgpm,nspin_m) :: cpm

    integer :: ip,istart,iend,i, is, ibase, nelmnt

! --> T. Yamasaki, 03rd Aug. 2009
    nelmnt = mp_kngp*nspin_m
! <--
!!$    print '(" -- scatter_cp_onto_cpm -- ")'
    do ip = 0, npes - 1
       ! (1)  coping input data onto a work array, and broadcasting
       if(is_kngp(ip) > kgpm) exit
       if(ip == mype) then
! --> T. Yamasaki, 03rd Aug. 2009

! ============================= modified by K. Tagami ================= 11.0
!          do is = 1, nspin, af+1
          do is = 1, ndim_magmom, af+1
! ==================================================================== 11.0

             ibase = mp_kngp*(is-1)
             do i = 1, iend_kngp-ista_kngp+1
                prj_wk(i+ibase) = cp(ista_kngp-1+i,is)
             end do
          end do
!!$          do i = 1, iend_kngp-ista_kngp+1
!!$             prj_wk(i)     = cp(ista_kngp-1+i)
!!$!!$             prj_wk(i,1,1) = cp(ista_kngp-1+i)
!!$          end do
       end if
!!$       call mpi_bcast(prj_wk,mp_kngp,mpi_double_precision,ip,MPI_CommGroup,ierr)
       call mpi_bcast(prj_wk,mp_kngp*nspin_m,mpi_double_precision,ip,MPI_CommGroup,ierr)
! <--
       ! (2) projection
       istart = ista_kgpm; if(istart < is_kngp(ip)) istart = is_kngp(ip)
       iend   = iend_kgpm; if(iend   > ie_kngp(ip)) iend   = ie_kngp(ip)
       if(iend < istart) cycle
! --> T. Yamasaki, 03rd Aug. 2009

! ==================================== modified by K. Tagami ============== 11.0
!       do is = 1, nspin, af+1
       do is = 1, ndim_magmom, af+1
! ========================================================================= 11.0
          ibase = mp_kngp*(is-1)
          do i = istart, iend
             cpm(i,is) = prj_wk(i + ibase)
          end do
       end do
!!$       do i = istart, iend
!!$          cpm(i) = prj_wk(i - is_kngp(ip)+1)
!!$!!$          cpm(i) = prj_wk(i -is_kngp(ip)+1,1,1)
!!$       end do
!!$!!$       cpm(istart:iend) = prj_wk(istart-is_kngp(ip)+1:iend-is_kngp(ip)+1,1,1)
! <--
    end do
  end subroutine scatter_cp_onto_cpm

  subroutine concentrate_d_to_chg(d,c)
! ============================= modified by K. Tagami ====================== 11.0
!    real(DP),intent(in),dimension(ista_kgpm:iend_kgpm,kimg,nspin) :: d
!    real(DP),intent(out), dimension(ista_kngp:iend_kngp,kimg,nspin) :: c
    real(DP),intent(in),dimension(ista_kgpm:iend_kgpm,kimg,ndim_magmom) :: d
    real(DP),intent(out), dimension(ista_kngp:iend_kngp,kimg,ndim_magmom) :: c
! =========================================================================== 11.0

    integer :: ip,is,istart,iend,nelmnt,ik,i,ip2,ipbase
!!$    integer :: istart_p, iend_p

                                                  __TIMER_SUB_START(1122)
    if(kgpm < kgp .and. npes /= 1) then
       nelmnt = mp_kgpm*kimg*nspin_m
       do ip = 0, npes - 1
          if(is_kngp(ip) > kgpm) exit
          do ip2 = 0, npes - 1
             istart = is_kgpm(ip2); if(istart < is_kngp(ip)) istart = is_kngp(ip)
             iend   = ie_kgpm(ip2); if(iend   > ie_kngp(ip)) iend   = ie_kngp(ip)
             if(iend < istart) cycle
             if(mype == ip2) then
                                                  __TIMER_DO_START(1171)
! ================================ modified by K. Tagami ================= 11.0
!                do is = 1, nspin, af+1
                do is = 1, ndim_magmom, af+1
! ======================================================================== 11.0
                   do ik = 1, kimg
                      ipbase = mp_kgpm*(ik-1) + mp_kgpm*kimg*(is-1) - istart + 1
                      do i = istart, iend
                         prj_wk(i+ipbase) = d(i,ik,is)
                      end do
                   end do
                end do
                                                  __TIMER_DO_STOP(1171)
                call mpi_send(prj_wk,nelmnt,mpi_double_precision,ip,1,MPI_CommGroup,ierr)
             end if
             if(mype == ip) then
                call mpi_recv(prj_wk,nelmnt,mpi_double_precision,ip2,1,MPI_CommGroup,istatus,ierr)
                                                  __TIMER_DO_START(1173)
! ================================ modified by K. Tagami ================= 11.0
!                do is = 1, nspin, af+1
                do is = 1, ndim_magmom, af+1
! ======================================================================== 11.0
                   do ik = 1, kimg
                      ipbase = mp_kgpm*(ik-1) + mp_kgpm*kimg*(is-1) - istart + 1
                      do i = istart, iend
                         c(i,ik,is) = prj_wk(i+ipbase)
                      end do
                   end do
                end do
                                                  __TIMER_DO_STOP(1173)
             end if
             call mpi_barrier(MPI_CommGroup,ierr)
          end do
       end do
    end if
                                                  __TIMER_SUB_STOP(1122)
  end subroutine concentrate_d_to_chg

  subroutine mix_dealloc_previous()
    if(previous_waymix == BROYD1) then
       call mix_broyden_deallocate()
    else if(previous_waymix == BROYD2) then
       call mix_broyden_deallocate()
    else if(previous_waymix == DFP) then
       call mix_DFP_deallocate()
    else if(previous_waymix == PULAY) then
       call mix_PULAY_deallocate()
    end if
  end subroutine mix_dealloc_previous

  subroutine mix_broyden_allocate
!!$    if(allocated(f_p)) return

! ========================= modified by K. Tagami =================== 11.0
!    nspin_m  = nspin/(af+1)
!
    if ( noncol ) then
       nspin_m  = ndim_magmom
    else
       nspin_m  = nspin/(af+1)
    endif
! ====================================================================== 11.0

! ================================ Modofied by K. Tagami ===========
!    allocate(f_p(ista_kgpm:iend_kgpm)); call precon_4_mult(f_p) !-(m_CD)
    allocate(f_p(ista_kgpm:iend_kgpm)); f_p = 0.0d0
    call precon_4_mult(f_p) !-(m_CD)
! ================================================================
    allocate(din(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(dout(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(dF_l(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(urec_l(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix,2))
    allocate(prj_wk(mp_kngp*kimg*nspin_m))
! ======================================Added by K. Tagami ========
    din = 0.0d0; dout = 0.0d0; dF_l = 0.0d0; urec_l = 0.0d0; prj_wk = 0.0d0
! ==================================================================
    if(hownew == RENEW) then
! ============================= modified by K. Tagami =========== 11.0
!       allocate(f(nbxmix,nbxmix,nspin))
!
       if ( noncol ) then
          allocate(f(nbxmix,nbxmix,ndim_magmom))
       else
          allocate(f(nbxmix,nbxmix,nspin))
       endif
! =============================================================== 11.0
       allocate(g(nbxmix))
! ================================= Added by K. Tagami ==========
        f = 0.0d0; g = 0.0d0
! ==============================================================
    end if
    allocate(ncrspd(nbxmix))
! ================================= Added by K. Tagami ==========
        ncrspd = 0
! ==============================================================
  end subroutine mix_broyden_allocate

  subroutine mix_broyden_deallocate
    if(allocated(f_p)) deallocate(f_p)
    if(allocated(din)) deallocate(din)
    if(allocated(dout)) deallocate(dout)
    if(allocated(dF_l)) deallocate(dF_l)
    if(allocated(urec_l)) deallocate(urec_l)
    if(allocated(prj_wk)) deallocate(prj_wk)
    if(allocated(f)) deallocate(f)
    if(allocated(g)) deallocate(g)
    if(allocated(ncrspd)) deallocate(ncrspd)
  end subroutine mix_broyden_deallocate

  subroutine mix_broyden_alloc2
    allocate(d0_l(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(u_l(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(v_l(ista_kgpm:iend_kgpm,kimg,nspin_m))
! =========================================== Added by K. Tagami =======
    d0_l = 0; u_l = 0; v_l = 0
! =======================================================================
    call alloc_rho_rhoo_and_cpm
  end subroutine mix_broyden_alloc2

  subroutine alloc_rho_rhoo_and_cpm
    if(kgpm < kgp .and. npes /= 1 ) then
       allocate(rho(ista_kgpm:iend_kgpm,kimg,nspin_m))
       allocate(rhoo(ista_kgpm:iend_kgpm,kimg,nspin_m))
       allocate(c_pm(ista_kgpm:iend_kgpm,nspin_m))
! ============================================= Added by K. Tagami ======
       rho = 0.0d0; rhoo = 0.0d0 ; c_pm = 0.0d0
! =======================================================================
       call scatter_chg_onto_d(chgq_l,rho)
       call scatter_chg_onto_d(chgqo_l,rhoo)
       call scatter_cp_onto_cpm(c_p,c_pm)
    else
       rho => chgq_l; rhoo => chgqo_l; c_pm => c_p
    end if
  end subroutine alloc_rho_rhoo_and_cpm

  subroutine mix_broyden_dealloc2
    deallocate(d0_l); deallocate(u_l); deallocate(v_l)
    call dealloc_rho_rhoo_and_cpm
  end subroutine mix_broyden_dealloc2

  subroutine dealloc_rho_rhoo_and_cpm
    if(kgpm < kgp .and. npes /= 1) deallocate(rho)
    if(kgpm < kgp .and. npes /= 1) deallocate(rhoo)
    if(kgpm < kgp .and. npes /= 1) deallocate(c_pm)
  end subroutine dealloc_rho_rhoo_and_cpm

  subroutine mix_broyden_alloc3
                                                  __TIMER_SUB_START(1109)
    allocate(d0_l(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(u_l(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(v_l(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(dd_l(ista_kgpm:iend_kgpm,kimg,nspin_m))
! =========================================== Added by K. Tagami ===
    d0_l = 0.0d0; u_l = 0.0d0; v_l = 0.0d0; dd_l = 0.0d0
! =================================================================
    call alloc_rho_rhoo_and_cpm
                                                  __TIMER_SUB_STOP(1109)
  end subroutine mix_broyden_alloc3

  subroutine mix_broyden_dealloc3
    deallocate(d0_l); deallocate(u_l); deallocate(v_l); deallocate(dd_l)
    call dealloc_rho_rhoo_and_cpm
  end subroutine mix_broyden_dealloc3

  subroutine mix_DFP_allocate
!!$    if(allocated(f_p)) return

! ============================ modified by K. Tagami ================ 11.0
!    nspin_m  = nspin/(af+1)
!
    if ( noncol ) then
       nspin_m  = ndim_magmom
    else
       nspin_m  = nspin/(af+1)
    endif
! ==================================================================== 11.0

! =============================== Modified by K. Tagami ==========
!    allocate(f_p(ista_kgpm:iend_kgpm)); call precon_4_mult(f_p) !-(m_CD)
    allocate(f_p(ista_kgpm:iend_kgpm)); f_p = 0.0d0
    call precon_4_mult(f_p) !-(m_CD)
! ===================================================================
    allocate(din(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(dout(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(dF_l(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(urec_l(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix,2))
    allocate(ncrspd(nbxmix))
    allocate(uuf(nbxmix,nspin_m,2))
    allocate(prj_wk(mp_kngp*kimg*nspin_m))
! =============================Added by K. Tagami =============
    din = 0.0d0; dout = 0.0d0; dF_l = 0.0d0; urec_l = 0.0d0
    ncrspd = 0; uuf = 0.0d0; prj_wk = 0.0d0
! ============================================================
  end subroutine mix_DFP_allocate

  subroutine mix_DFP_deallocate()
    if(allocated(f_p)) deallocate(f_p)
    if(allocated(din)) deallocate(din)
    if(allocated(dout)) deallocate(dout)
    if(allocated(dF_l)) deallocate(dF_l)
    if(allocated(urec_l)) deallocate(urec_l)
    if(allocated(ncrspd)) deallocate(ncrspd)
    if(allocated(uuf)) deallocate(uuf)
    if(allocated(prj_wk)) deallocate(prj_wk)
  end subroutine mix_DFP_deallocate

  subroutine mix_DFP_alloc2
    allocate(d0_l(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(u_l(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(w_l(ista_kgpm:iend_kgpm,kimg,nspin_m))
! ======================================== Added by K. Tagami ======
    d0_l = 0; u_l = 0; w_l = 0
! ==================================================================
    call alloc_rho_rhoo_and_cpm
  end subroutine mix_DFP_alloc2

  subroutine mix_DFP_dealloc2
    deallocate(d0_l); deallocate(u_l); deallocate(w_l)
    call dealloc_rho_rhoo_and_cpm
  end subroutine mix_DFP_dealloc2

!!!
!!!
!!! ========================== added by K. Tagami ================ 5.0
  subroutine alloc_rhostore_recomp( rmxt, rmxtrc )
    real(kind=DP),intent(in) :: rmxt
    real(kind=DP),intent(out),dimension(nspin_m) :: rmxtrc

    allocate( rhoo_store( nsize_rho_hsr,nspin) )
    allocate( rho_store( nsize_rho_hsr,nspin) )

    rho_store = rho_hsr;      rhoo_store = rhoo_hsr

     rho_hsr(:,1) =  rho_store(:,1) +  rho_store(:,2)
     rho_hsr(:,2) =  rho_store(:,1) -  rho_store(:,2)
    rhoo_hsr(:,1) = rhoo_store(:,1) + rhoo_store(:,2)
    rhoo_hsr(:,2) = rhoo_store(:,1) - rhoo_store(:,2)
    rmxtrc(1) = rmxt;     rmxtrc(2) = rmxt*spin_density_mixfactor

  end subroutine alloc_rhostore_recomp

  subroutine compose_rho_dealloc_store
    rho_store = rho_hsr

    rho_hsr(:,1) = 0.5d0*( rho_store(:,1) + rho_store(:,2) )
    rho_hsr(:,2) = 0.5d0*( rho_store(:,1) - rho_store(:,2) )

    rhoo_hsr(:,1) = 0.5d0*( rhoo_store(:,1) + rhoo_store(:,2) )
    rhoo_hsr(:,2) = 0.5d0*( rhoo_store(:,1) - rhoo_store(:,2) )

!    rhoo_hsr = rhoo_store
    deallocate( rho_store, rhoo_store )

  end subroutine compose_rho_dealloc_store

  subroutine alloc_hsrstore_recomp( rmxt, rmxtrc )
    real(kind=DP),intent(in) :: rmxt
    real(kind=DP),intent(out),dimension(nspin_m) :: rmxtrc

    allocate( hsr_store ( natm, nlmt, nlmt, nspin ) ); hsr_store=0.0d0
    allocate( hsro_store( natm, nlmt, nlmt, nspin ) ); hsro_store=0.0d0

    hsr_store = hsr;      hsro_store = hsro

    hsr(:,:,:,1) =   hsr_store(:,:,:,1) +  hsr_store(:,:,:,2)
    hsr(:,:,:,2) =   hsr_store(:,:,:,1) -  hsr_store(:,:,:,2)
    hsro(:,:,:,1) =  hsro_store(:,:,:,1) + hsro_store(:,:,:,2)
    hsro(:,:,:,2) =  hsro_store(:,:,:,1) - hsro_store(:,:,:,2)

    rmxtrc(1) = rmxt;     rmxtrc(2) = rmxt*spin_density_mixfactor
  end subroutine alloc_hsrstore_recomp

  subroutine compose_hsr_dealloc_store
    hsr_store = hsr
    hsro_store = hsro

    hsr(:,:,:,1) = 0.5d0*( hsr_store(:,:,:,1) + hsr_store(:,:,:,2) )
    hsr(:,:,:,2) = 0.5d0*( hsr_store(:,:,:,1) - hsr_store(:,:,:,2) )
    hsro(:,:,:,1) = 0.5d0*( hsro_store(:,:,:,1) + hsro_store(:,:,:,2) )
    hsro(:,:,:,2) = 0.5d0*( hsro_store(:,:,:,1) - hsro_store(:,:,:,2) )

!    hsro = hsro_store
    deallocate( hsr_store, hsro_store )
  end subroutine compose_hsr_dealloc_store

 subroutine simple_mix_kt(rmx_this)
    real(kind=DP), intent(in) :: rmx_this(nspin_m)

    din_hsr   = rhoo_hsr ! chgqo
    dout_hsr  = rho_hsr  ! chgq
    rho_hsr(:,1) = rmx_this(1) *dout_hsr(:,1) + ( 1.0D0 - rmx_this(1) )* din_hsr(:,1)
    if(nspin_m==2) rho_hsr(:,2) = rmx_this(2) *dout_hsr(:,2) + ( 1.0D0 - rmx_this(2) )* din_hsr(:,2)

! ===================== added by K. Tagami ================================ 11.0
    if ( noncol ) then
       rho_hsr(:,2) = rmx_this(2) *dout_hsr(:,2) + ( 1.0D0 -rmx_this(2) )* din_hsr(:,2)
       rho_hsr(:,3) = rmx_this(3) *dout_hsr(:,3) + ( 1.0D0 -rmx_this(3) )* din_hsr(:,3)
       rho_hsr(:,4) = rmx_this(4) *dout_hsr(:,4) + ( 1.0D0 -rmx_this(4) )* din_hsr(:,4)
    endif
! ========================================================================= 11.0

  end subroutine simple_mix_kt

  subroutine simple_mix2_kt(rmx_this)
    real(kind=DP), intent(in) :: rmx_this(nspin_m)

    if(nspin_m==2)then
        din_hsr(:,2)  = rhoo_hsr(:,2)         ! chgqo
       dout_hsr(:,2)  =  rho_hsr(:,2)          ! chgq
       rho_hsr(:,2) = rmx_this(2) *dout_hsr(:,2) + ( 1.0D0 - rmx_this(2) )* din_hsr(:,2)
    endif
! ===================== added by K. Tagami ================================ 11.0
    if ( noncol ) then
       din_hsr(:,2)  = rhoo_hsr(:,2);     dout_hsr(:,2)  =  rho_hsr(:,2)
       rho_hsr(:,2) = rmx_this(2) *dout_hsr(:,2) + ( 1.0D0 -rmx_this(2) )* din_hsr(:,2)
       din_hsr(:,3)  = rhoo_hsr(:,3);     dout_hsr(:,3)  =  rho_hsr(:,3)
       rho_hsr(:,3) = rmx_this(3) *dout_hsr(:,3) + ( 1.0D0 -rmx_this(3) )* din_hsr(:,3)
       din_hsr(:,4)  = rhoo_hsr(:,4);     dout_hsr(:,4)  =  rho_hsr(:,4)
       rho_hsr(:,4) = rmx_this(4) *dout_hsr(:,4) + ( 1.0D0 -rmx_this(4) )* din_hsr(:,4)
    endif
! ========================================================================== 11.0

  end subroutine simple_mix2_kt

  subroutine mix_dealloc_previous_hsr()
    if(previous_waymix == BROYD1) then
       call mix_broyden_deallocate_hsr()
    else if(previous_waymix == BROYD2) then
       call mix_broyden_deallocate_hsr()
    else if(previous_waymix == DFP) then
       call mix_DFP_deallocate_hsr()
    else if(previous_waymix == PULAY) then
       call mix_PULAY_deallocate_hsr()
    end if
  end subroutine mix_dealloc_previous_hsr

  subroutine mix_broyden_allocate_hsr
! ========================= modified by K. Tagami =================== 11.0
!    nspin_m  = nspin/(af+1)
!
    if ( noncol ) then
       nspin_m  = ndim_magmom
    else
       nspin_m  = nspin/(af+1)
    endif
! ====================================================================== 11.0
    allocate( din_hsr(nsize_rho_hsr,nspin_m) )
    allocate( dout_hsr(nsize_rho_hsr,nspin_m))
    allocate( dF_hsr(nsize_rho_hsr,nspin_m))
    allocate( urec_hsr(nsize_rho_hsr,nspin_m,nbxmix,2) )
!
    din_hsr = 0.0d0; dout_hsr = 0.0d0; dF_hsr = 0.0d0; urec_hsr = 0.0d0
!
  end subroutine mix_broyden_allocate_hsr

  subroutine mix_broyden_deallocate_hsr
    if(allocated(din_hsr)) deallocate(din_hsr)
    if(allocated(dout_hsr)) deallocate(dout_hsr)
    if(allocated(dF_hsr)) deallocate(dF_hsr)
    if(allocated(urec_hsr)) deallocate(urec_hsr)
  end subroutine mix_broyden_deallocate_hsr

  subroutine mix_broyden_alloc2_hsr
    allocate( d0_hsr( nsize_rho_hsr,nspin_m ) ); d0_hsr = 0.0d0
    allocate(  u_hsr( nsize_rho_hsr,nspin_m ) );  u_hsr = 0.0d0
    allocate(  v_hsr( nsize_rho_hsr,nspin_m ) );  v_hsr = 0.0d0
  end subroutine mix_broyden_alloc2_hsr

 subroutine mix_broyden_dealloc2_hsr
    deallocate(d0_hsr); deallocate(u_hsr); deallocate(v_hsr)
  end subroutine mix_broyden_dealloc2_hsr

  subroutine mix_broyden_alloc3_hsr
    allocate( d0_hsr( nsize_rho_hsr,nspin_m ) ); d0_hsr = 0.0d0
    allocate(  u_hsr( nsize_rho_hsr,nspin_m ) );  u_hsr = 0.0d0
    allocate(  v_hsr( nsize_rho_hsr,nspin_m ) );  v_hsr = 0.0d0
    allocate( dd_hsr( nsize_rho_hsr,nspin_m ) ); dd_hsr = 0.0d0
  end subroutine mix_broyden_alloc3_hsr

  subroutine mix_broyden_dealloc3_hsr
    deallocate(d0_hsr); deallocate(u_hsr); deallocate(v_hsr); deallocate(dd_hsr)
  end subroutine mix_broyden_dealloc3_hsr

  subroutine mix_DFP_allocate_hsr
! ========================= modified by K. Tagami =================== 11.0
!    nspin_m  = nspin/(af+1)
!
    if ( noncol ) then
       nspin_m  = ndim_magmom
    else
       nspin_m  = nspin/(af+1)
    endif
! ====================================================================== 11.0
    allocate(  din_hsr( nsize_rho_hsr,nspin_m ) );  din_hsr = 0.0d0
    allocate( dout_hsr( nsize_rho_hsr,nspin_m ) ); dout_hsr = 0.0d0
    allocate(   dF_hsr( nsize_rho_hsr,nspin_m ) );   dF_hsr = 0.0d0
    allocate( urec_hsr( nsize_rho_hsr,nspin_m,nbxmix,2) ); urec_hsr = 0.0d0
  end subroutine mix_DFP_allocate_hsr

  subroutine mix_DFP_deallocate_hsr()
    if ( allocated( din_hsr)) deallocate( din_hsr)
    if ( allocated(dout_hsr)) deallocate(dout_hsr)
    if ( allocated(  dF_hsr)) deallocate(  dF_hsr)
    if ( allocated(urec_hsr)) deallocate(urec_hsr)
  end subroutine mix_DFP_deallocate_hsr

  subroutine mix_DFP_alloc2_hsr
    allocate( d0_hsr( nsize_rho_hsr,nspin_m) ); d0_hsr = 0.0d0
    allocate(  u_hsr( nsize_rho_hsr,nspin_m) );  u_hsr = 0.0d0
    allocate(  w_hsr( nsize_rho_hsr,nspin_m) );  w_hsr = 0.0d0
  end subroutine mix_DFP_alloc2_hsr

  subroutine mix_DFP_dealloc2_hsr
    deallocate( d0_hsr, u_hsr, w_hsr )
  end subroutine mix_DFP_dealloc2_hsr

 subroutine mix_pulay_allocate_hsr
! ========================= modified by K. Tagami =================== 11.0
!    nspin_m  = nspin/(af+1)
!
    if ( noncol ) then
       nspin_m  = ndim_magmom
    else
       nspin_m  = nspin/(af+1)
    endif
! ====================================================================== 11.0
    allocate(  din_hsr( nsize_rho_hsr,nspin_m ) );  din_hsr = 0.0d0
    allocate( dout_hsr( nsize_rho_hsr,nspin_m ) ); dout_hsr = 0.0d0
    allocate(   dF_hsr( nsize_rho_hsr,nspin_m ) );   dF_hsr = 0.0d0
    allocate( urec_hsr( nsize_rho_hsr,nspin_m,nbxmix,2) ); urec_hsr = 0.0d0
    if(sw_gradient_simplex==ON)then
       allocate(rhoj_hsr(nsize_rho_hsr,nspin_m,nbxmix));rhoj_hsr=0.d0
       allocate(Frhoj_hsr(nsize_rho_hsr,nspin_m,nbxmix));Frhoj_hsr=0.d0
       allocate(rhoj_hsro(nsize_rho_hsr,nspin_m));rhoj_hsro=0.d0
       allocate(Frhoj_hsro(nsize_rho_hsr,nspin_m));Frhoj_hsro=0.d0
    endif
  end subroutine mix_pulay_allocate_hsr

  subroutine mix_pulay_deallocate_hsr
    if ( allocated( din_hsr)) deallocate( din_hsr)
    if ( allocated(dout_hsr)) deallocate(dout_hsr)
    if ( allocated(  dF_hsr)) deallocate(  dF_hsr)
    if ( allocated(urec_hsr)) deallocate(urec_hsr)
    if ( allocated(ynorm)) deallocate(ynorm)
    if ( allocated(sum12)) deallocate(sum12)
    if ( allocated(rhoj_hsr)) deallocate(rhoj_hsr)
    if ( allocated(Frhoj_hsr)) deallocate(Frhoj_hsr)
    if ( allocated(rhoj_hsro)) deallocate(rhoj_hsro)
    if ( allocated(Frhoj_hsro)) deallocate(Frhoj_hsro)
    if ( allocated(d0_l_h)) deallocate(d0_l_h)
  end subroutine mix_pulay_deallocate_hsr

  subroutine mix_pulay_alloc2_hsr
    allocate( d0_hsr( nsize_rho_hsr,nspin_m) ); d0_hsr = 0.0d0
  end subroutine mix_pulay_alloc2_hsr

  subroutine mix_pulay_dealloc2_hsr
    deallocate(d0_hsr)
  end subroutine mix_pulay_dealloc2_hsr

! ===================================================================== 5.0

  integer function nspin_for_qnewton()
! ================================= modified by K. Tagami ============= 11.0
!     nspin_for_qnewton=nspin
!     if (sw_force_simple_mixing==ON .and. sw_recomposing==ON) nspin_for_qnewton=1
!
    if ( noncol ) then
       nspin_for_qnewton=ndim_magmom
    else
       nspin_for_qnewton=nspin
       if (sw_force_simple_mixing==ON .and. sw_recomposing==ON) nspin_for_qnewton=1
    endif
! ========================================================================== 11.0
  end function nspin_for_qnewton

  subroutine renew_u_br(j,i)
    integer, intent(in) :: j,i

! === DEBUG by tkato 2011/11/24 ================================================
!   real(DP)      :: v_dF(nspin)
    real(DP)      :: v_dF(nspin_m)
! ==============================================================================
                                                  __TIMER_SUB_START(1113)
#ifdef _CDMIX_USE_POINTER_
    urec_l_3 => urec_l(ista_kgpm:iend_kgpm,1:kimg,1:nspin_m,j,iV)
    call mult1s(urec_l_3,dF_l,f_p,v_dF)!-(m_CD);<v|dF> ->v_dF
#else
    call mult1s5(urec_l,nbxmix,2,j,iV,dF_l,f_p,v_dF)
#endif
    if ( nspin_m==2 .and. sw_mix_bothspins_sametime == YES ) then
      v_dF(1) = v_dF(1) + v_dF(2)
      v_df(2) = v_dF(1)
    endif

! ============================== added by K.Tagami ================= 11.0
    if ( noncol ) then
       v_dF(1) = sum( v_dF(:) )
       v_dF(:) = v_dF(1)
    endif
! ================================================================== 11.0

    call subtr_j_th_term(v_dF,iU,j,urec_l,u_l)  !-(m_CD)

    !                        |u(m)> = |u(m)> - <v(j)|dF(m)>|u(j)>
! === DEBUG by tkato 2011/11/24 ================================================
!   if(hownew == RENEW) f(j,i,1:nspin) = v_dF(1:nspin)
    if(hownew == RENEW) f(j,i,1:nspin_m) = v_dF(1:nspin_m)
! ==============================================================================
                                                  __TIMER_SUB_STOP(1113)
  end subroutine renew_u_br

  subroutine renew_d_br(j)
    integer, intent(in) :: j

    real(DP)  :: vF(nspin_m)  !   real(DP)  :: vF(nspin) ! === DEBUG by tkato 2011/11/24 ====
    
                                                  __TIMER_SUB_START(1118)
#ifdef _CDMIX_USE_POINTER_
    urec_l_3 => urec_l(ista_kgpm:iend_kgpm,1:kimg,1:nspin_m,j,iV)
    call mult1s(urec_l_3,F_l,f_p,vF) !-(m_CD);<v|F>  ->vF
#else
    call mult1s5(urec_l,nbxmix,2,j,iV,F_l,f_p,vF) !-(m_CD);<v|F>  ->vF
#endif

    if ( nspin_m==2 .and. sw_mix_bothspins_sametime == YES ) then  ! nspin==2 --> nspin_m==2    ! === DEBUG by tkato 2011/11/24 ===
      vF(1) = vF(1) + vF(2)
      vF(2) = vF(1)
    endif

! ================================= added by K. Tagami =============== 11.0
    if ( noncol ) then
       vF(1) = sum( vF(:) )
       vF(:) = vF(1)
    endif
! ===================================================================== 11.0

    call subtr_j_th_term(vF,iU,j,urec_l,d0_l) !-(m_CD)
    !                        |d(m)> = |d(m)> - <v(j)|F(m)>|u(j)>
                                                  __TIMER_SUB_STOP(1118)
  end subroutine renew_d_br

  subroutine renew_d_last_br(p)
    real(DP), intent(in), dimension(ista_kngp:iend_kngp) :: p
    integer   :: is, ik, i, ns

    real(DP)  :: vF(nspin_m)  !   real(DP)  :: vF(nspin)  ! === DEBUG by tkato 2011/11/24 ===
                                                  __TIMER_SUB_START(1120)

    call mult1s(v_l,F_l,f_p,vF)              !-(m_CD) <v|F> ->vF

    if ( nspin_m==2 .and. sw_mix_bothspins_sametime == YES ) then    !  nspin==2 --> nspin_m==2  ! === DEBUG by tkato 2011/11/24 ===
      vF(1) = vF(1) + vF(2)
      vF(2) = vF(1)
    endif

! ============================ added by K.Tagami ====================== 11.0
    if ( noncol ) then
       vF(1) = sum( vF(:) )
       vF(:) = vF(1)
    endif
! ====================================================================== 11.0

    if(kgpm == kgp .or. npes == 1) then
                                                  __TIMER_DO_START(1166)

       do is = 1, ndim_magmom, af+1  !       do is = 1, nspin, af+1 ! === modified by K. Tagami === 11.0
          din (:,:,is) = chgqo_l(ista_kgpm:iend_kgpm,:,is)
          dout(:,:,is) = chgq_l (ista_kgpm:iend_kgpm,:,is)
       end do
                                                  __TIMER_DO_STOP(1166)
    else
       call scatter_chg_onto_d(chgqo_l,din )  ! -(m_C.D.)
       call scatter_chg_onto_d(chgq_l, dout)  ! -(m_C.D.)
    end if

!!$    do is = 1, nspin, af+1
    ns = nspin_for_qnewton()
                                                  __TIMER_DO_START(1167)
    do is = 1, ns,af+1
       do ik = 1, kimg
          do i = ista_kgpm,iend_kgpm
             rho(i,ik,is) = d0_l(i,ik,is) - vF(is)*u_l(i,ik,is)
          end do
       end do
    end do
                                                  __TIMER_DO_STOP(1167)

! ====================================== modified by K. Tagami ============ 11.0
!    if(sw_force_simple_mixing==ON .and. sw_recomposing==ON) call simple_mix2(c_p)
!
    if ( .not. noncol ) then
       if (sw_force_simple_mixing==ON .and. sw_recomposing==ON) then
          call simple_mix2(c_p)
       endif
    endif

! =========================================================================== 11.0

    if(kgpm < kgp) then
       call concentrate_d_to_chg(rho,chgq_l) !-(m_C.D.)
       call simple_mix_large_Gc(p,chgqo_l,chgq_l)      !-(m_C.D.) chgq,chgqo,p ->chgq
    end if
                                                  __TIMER_SUB_STOP(1120)
  end subroutine renew_d_last_br

#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
  subroutine printvF(nfout,aorb,title,n,vF,j,i)
    integer, intent(in) :: nfout,n
    character(len=1), intent(in) :: aorb
    character(len=n), intent(in) :: title
    real(DP), intent(in) :: vF(nspin_m)
    integer, intent(in), optional :: j,i
!!$    character(len=80) :: fmt = ''
!!$    write(fmt,*) "(a",n,")"
    if(present(i) .and. present(j)) then
       if(nspin_m == 1) then
          write(nfout,'(" (",a1,") (j,i) = (",i2,",",i2,") ",a4,"(1) = ",f20.12)') aorb,j,i, title, vF(1)
       else
          write(nfout,'(" (",a1,") (j,i) = (",i2,",",i2,") ",a4,"(1) = ",f20.12, 2x,a4,"(2) = ",f20.12)') &
               & aorb,j,i, title, vF(1), title, vF(2)
       end if
    else if(present(j)) then
       if(nspin_m == 1) then
          write(nfout,'(" (",a1,") j = ",i2,1x,a2,"(1) = ",f20.12)') aorb, j, title, vF(1)
       else
          write(nfout,'(" (",a1,") j = ",i2,1x,a2,"(1) = ",f20.12, 2x,a4,"(2) = ",f20.12)')  aorb, j, title, vF(1), title, vF(2)
       end if
    else
       if(nspin_m == 1) then
          write(nfout,'(" (",a1,") ",a2,"(1) = ",f20.12)') aorb, title, vF(1)
       else
          write(nfout,'(" (",a1,") ",a2,"(1) = ",f20.12,2x, a4,"(2) = ",f20.12)') aorb, title, vF(1),title, vF(2)
       end if
    end if
  end subroutine printvF
#endif
          
! =========================== added by K. Tagami ================================== 5.0
  subroutine renew_u_br_with_hsr(nfout,j,i)
    integer, intent(in) :: nfout,j,i

    integer       :: is
    real(DP)      :: v_dF(nspin_m) !  v_dF(nspin),  revised by tkato 2011/11/24
    
    v_dF = 0.d0
#ifdef _CDMIX_USE_POINTER_
    urec_l_3 => urec_l(ista_kgpm:iend_kgpm,1:kimg,1:nspin_m,j,iV)
    call mult1s(urec_l_3,dF_l,f_p,v_dF)!-(m_CD);<v|dF> ->v_dF
#else
    call mult1s5(urec_l,nbxmix,2,j,iV,dF_l,f_p,v_dF)
#endif

#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    call printvF(nfout,"a","v_dF",4,v_dF,j,i)
    write(nfout,'(" nsize_rho_hsr = ",i8)') nsize_rho_hsr
#endif
#ifdef _DUPLICATION_HSR_DOTPRODUCT_
    do is = 1, ndim_magmom, af+1 !   do is = 1, nspin, af+1 ! === modified by K. Tagami === 11.0
       v_dF(is) = v_dF(is) + sum( urec_hsr(:,is,j,iV)*dF_hsr(:,is) )
    end do
#else
    call mult_urec_hsr5(nfout,urec_hsr,nbxmix,2,j,iV,dF_hsr,v_dF)  ! -(m_CD_mixing) <v_hsr|FF_hsr> -> v_dF_hsr by T. Yamasaki 2023/07/07
#endif

    if ( nspin_m==2 .and. sw_mix_bothspins_sametime == YES ) then !   nspin==2 --> nspin_m==2, by tkato 2011/11/24
      v_dF(1) = v_dF(1) + v_dF(2)
      v_df(2) = v_dF(1)
    endif

    if ( noncol ) then           ! === added by K. Tagami === 11.0
       v_dF(1) = sum( v_dF(:) )
       v_dF(:) = v_dF(1)
    endif
    
#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    call printvF(nfout,"b","v_dF",4,v_dF,j,i)
#endif
    
    call subtr_j_th_term(v_dF,iU,j,urec_l,u_l)  !-(m_CD)
    !                        |u(m)> = |u(m)> - <v(j)|dF(m)>|u(j)>

    do is = 1, ndim_magmom, af+1 !   do is = 1, nspin, af+1 ! === modified by K. Tagami === 11.0
       u_hsr(:,is) = u_hsr(:,is) - v_dF(is) *urec_hsr(:,is,j,iU)
    end do

    if(hownew == RENEW) f(j,i,1:nspin_m) = v_dF(1:nspin_m) !  nspin --> nspin_m ! === DEBUG by tkato 2011/11/24 ===

  end subroutine renew_u_br_with_hsr

  subroutine renew_d_br_with_hsr(nfout,j)
    integer, intent(in) :: nfout,j
    real(DP)  :: vF(nspin_m) !   real(DP)  :: vF(nspin) ! === DEBUG by tkato 2011/11/24 ===

    integer :: is

    vF = 0.d0

#ifdef _CDMIX_USE_POINTER_
    urec_l_3 => urec_l(ista_kgpm:iend_kgpm,1:kimg,1:nspin_m,j,iV)
    call mult1s(urec_l_3,F_l,f_p,vF) !-(m_CD);<v|F>  ->vF
#else
    call mult1s5(urec_l,nbxmix,2,j,iV,F_l,f_p,vF) !-(m_CD);<v|F>  ->vF
#endif

#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    call printvF(nfout,"a","vF",2,vF,j)
#endif
    
#ifdef _DUPLICATION_HSR_DOTPRODUCT_
    do is = 1, ndim_magmom, af+1 !  do is = 1, nspin, af+1  === modified by K. Tagami === 11.0
       vF(is) = vF(is) + sum( urec_hsr(:,is,j,iV)*FF_hsr(:,is) )
    end do
#else
    call mult_urec_hsr5(nfout,urec_hsr,nbxmix,2,j,iV,FF_hsr, vF) ! by T. Yamasaki 2023/07/07
#endif

    if ( nspin_m==2 .and. sw_mix_bothspins_sametime == YES ) then  ! nspin --> nspin_m, ! === DEBUG by tkato 2011/11/24 ===
      vF(1) = vF(1) + vF(2)
      vF(2) = vF(1)
    endif

! ======================== added by K. Tagami ==================== 11.0
    if ( noncol ) then
       vF(1) = sum( vF(:) )
       vF(:) = vF(1)
    endif
! ================================================================ 11.0
    
#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    call printvF(nfout,"b","vF",2,vF,j)
#endif
    call subtr_j_th_term(vF,iU,j,urec_l,d0_l) !-(m_CD)
    !                        |d(m)> = |d(m)> - <v(j)|F(m)>|u(j)>

    do is = 1, ndim_magmom, af+1 !     do is = 1, nspin, af+1   ! === modified by K. Tagami === 11.0
       d0_hsr(:,is) = d0_hsr(:,is) - vF(is) *urec_hsr(:,is,j,iU)
    end do
  end subroutine renew_d_br_with_hsr

  subroutine renew_d_last_br_with_hsr(nfout, p, rmxtrc_hsr )
    integer, intent(in) :: nfout
    real(DP), intent(in), dimension(ista_kngp:iend_kngp) :: p
    real(DP), intent(in) :: rmxtrc_hsr(nspin_m)

    integer   :: is, ik, i, ns

    real(DP)  :: vF(nspin_m)  !  real(DP)  :: vF(nspin)  ! === DEBUG by tkato 2011/11/24 ===

    vF = 0.0d0
    call mult1s(v_l,F_l,f_p,vF)              !-(m_CD) <v|F> ->vF

#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    call printvF(nfout,"a","vF",2,vF)
#endif

#ifdef _DUPLICATION_HSR_DOTPRODUCT_
    do is = 1, ndim_magmom, af+1
       vF(is) = vF(is) + sum( v_hsr(:,is)*FF_hsr(:,is) )
    End do
#else
    call mult_urec_hsr(nfout,v_hsr,FF_hsr,vF)  ! -(m_CD_mixing) <v_hsr|FF_hsr> -> vF_hsr, by T. Yamasaki 2023/07/07
#endif

    if ( nspin_m==2 .and. sw_mix_bothspins_sametime == YES ) then  ! from nspin to nspin_m !== DEBUG by tkato 2011/11/24 ===
      vF(1) = vF(1) + vF(2)
      vF(2) = vF(1)
    endif

    if ( noncol ) then       ! === added by K. Tagami ========= 11.0
       vF(1) = sum( vF(:) )
       vF(:) = vF(1)
    endif

    if(kgpm == kgp .or. npes == 1) then
       do is = 1, ndim_magmom, af+1   !    do is = 1, nspin, af+1  ! === modified by K. Tagami === 11.0
          din (:,:,is) = chgqo_l(ista_kgpm:iend_kgpm,:,is)
          dout(:,:,is) = chgq_l (ista_kgpm:iend_kgpm,:,is)
       end do
    else
       call scatter_chg_onto_d(chgqo_l,din )  ! -(m_C.D.)
       call scatter_chg_onto_d(chgq_l, dout)  ! -(m_C.D.)
    end if

     do is = 1, ndim_magmom, af+1  !  do is = 1, nspin, af+1  ! === modified by K. Tagami === 11.0
       din_hsr (:,is) = rhoo_hsr(:,is) ! chgqo
       dout_hsr(:,is) = rho_hsr (:,is) ! chgq
    end do

#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    call printvF(nfout,"b","vF",2,vF)
#endif
    ns = nspin_for_qnewton()
    do is = 1, ns,af+1
       do ik = 1, kimg
          do i = ista_kgpm,iend_kgpm
             rho(i,ik,is) = d0_l(i,ik,is) - vF(is)*u_l(i,ik,is)
          end do
       end do
    end do

    do is = 1, ns, af+1
       rho_hsr(:,is) = d0_hsr(:,is) - vF(is) *u_hsr(:,is)
    end do

! ====================================== modified by K. Tagami ============ 11.0
!    if(sw_force_simple_mixing==ON .and. sw_recomposing==ON) call simple_mix2(c_p)
!
!    if ( sw_force_simple_mixing_hsr==ON .and. sw_recomposing_hsr==ON ) then
!       call simple_mix2_kt( rmxtrc_hsr )
!    endif

    if ( .not. noncol ) then
       if (sw_force_simple_mixing==ON .and. sw_recomposing==ON) then
          call simple_mix2(c_p)
       endif
       if ( sw_force_simple_mixing_hsr==ON .and. sw_recomposing_hsr==ON ) then
          call simple_mix2_kt( rmxtrc_hsr )
       endif
    endif
! =========================================================================== 11.0

    if(kgpm < kgp) then
       call concentrate_d_to_chg(rho,chgq_l) !-(m_C.D.)
       call simple_mix_large_Gc(p,chgqo_l,chgq_l)     !-(m_C.D.) chgq,chgqo,p ->chgq
    end if
  end subroutine renew_d_last_br_with_hsr

! ============================================================================ 5.0

  subroutine simple_mix_large_Gc(p,chgqo_l,chgq_l)
    real(kind=DP), intent(in), dimension(ista_kngp:iend_kngp) :: p
    real(kind=DP), intent(in) :: chgqo_l(ista_kngp:iend_kngp,kimg,ndim_magmom)
    real(kind=DP), intent(inout) :: chgq_l(ista_kngp:iend_kngp,kimg,ndim_magmom)

    integer :: istart, iend, is, ik, i
                                                  __TIMER_SUB_START(1123)
    istart = kgpm + 1; if(istart < ista_kngp) istart = ista_kngp
    iend   = kgp;      if(iend   > iend_kngp) iend   = iend_kngp
                                                  __TIMER_DO_START(1174)
    if(iend >= istart) then
! =================================== modified by K. Tagami ============ 11.0
!       do is = 1, nspin, af+1
       do is = 1, ndim_magmom, af+1
! ======================================================================= 11.0
          do ik = 1, kimg
             do i = istart, iend
                chgq_l(i,ik,is) = p(i)*chgq_l(i,ik,is) + (1.0d0-p(i))*chgqo_l(i,ik,is)
             end do
          end do
       end do
    end if
                                                  __TIMER_DO_STOP(1174)
                                                  __TIMER_SUB_STOP(1123)
  end subroutine simple_mix_large_Gc

! <<< Quasi-Newton Methods >>>
!  1. Broyden's 1st method
!  2. Broyden's 2nd method
!  3. DFP method
!
! In Quasi-Newton method, following set of equations are used.
!
!\begin{equation}
!\rho^{(m+1)}  = \rho^{(m)} - \lambda^{(m)} \left[ {\bf J} ^{(m)}\right]^{-1} {\cal F}^{(m)}
!\end{equation}
!
!\begin{equation}
!  \left[ {\bf J}^{(m+1)}\right]^{-1} \Delta {\cal F}^{(m+1)} = \Delta
!  \rho ^{(m+1)}
!\end{equation}
!
! $\lambda^{(m)}$ in the first euqation is a parameter of one
! dimensional search. And ${\bf J}^{(m+1)}$ is an approximate value of
! true Jacobian{\cal J}. In this program we fix this value of
! $\lambda^{(m)}$ to be 1.
! $\Delta {\cal F}^{(m+1)}$ and $\Delta \rho^{(m+1)}$ are definded as
! follows.
!
!\begin{equation}
!  \Delta {\cal F}^{(m+1)} = {\cal F}^{(m+1)} - {\cal F}^{(m)},
!\end{equation}
! where ${\cal F}^{(m)} = \rho^{out,(m)} - \rho^{in,(m)}$.
!
!\begin{equation}
!  \Delta \rho^{(m+1)} = \rho^{(m+1)} - \rho^{(m)}
!\end{equation}
!
! Three kinds of mixing ways implemented in this program,
! i.e. Broyden's 1st, Broyden's 2nd and DFP methods, differ in the
! way of approximation of the $\left[ {\bf J} ^{(m)}\right]^{-1}$. For
! combinience, we simply this description of the inverse of Jacobian
! as to be ${\bf J}^{-1(m)}$, hereafter.
!
! <<1>> Broyden's 1st method
!\begin{equation}
!{\bf J}^{-1(m)} = {\bf J}^{-1(m-1)} + \frac{\left[ | \Delta
!    \rho^{(m)}\rangle - {\bf J}^{-1(m-1)} | \Delta {\cal
!      F}^{(m)}\rangle \right] \otimes \langle \Delta \rho^{(m)} | {\bf
!    J}^{-1(m-1)}}{\langle\rho^{(m)}|{\bf J}^{-1(m-1)}| \Delta {\cal
!    F}^{(m)}\rangle}
!\end{equation}
!
! <<2>> Broyden's 2nd method
!\begin{equation}
!  {\bf J}^{-1(m)} = {\bf J}^{-1(m-1)} + \left[ |\Delta \rho^{(m)}
!    \rangle - {\bf J}^{-1(m-1)} | \Delta {\cal F}^{(m)} \right]
!  \otimes \frac{\langle \Delta {\cal F}^{(m)}|}{\| \Delta {\cal
!      F}^{(m)}\|^2}
!\end{equation}
!
! <<3>> DFP method
!DFP(Davidon-Fletcher-Powell) algorithm
!\begin{equation}
!  {\bf J}^{-1(m)} = {\bf J}^{-1(m-1)} + \frac{| \Delta
!    \rho^{(m)}\rangle \otimes \langle \Delta \rho^{(m)} |}
!   {\langle\rho^{(m)}|\Delta {\cal F}^{(m)}\rangle}
!   - \frac{|{\bf J}^{-1(m-1)}\Delta {\cal F}^{(m)}\rangle \otimes
!   \langle \Delta {\cal F}^{(m)}{\bf J}^{-1(m-1)}|}{\langle \Delta
!   {\cal F}^{(m)}|{\bf J}^{-1(m-1)}|\Delta {\cal F}^{(m)}\rangle}
!\end{equation}
!
! For reduction of memory space, we rewrite the Jacobian as linear
! combination of dyadic products.
!\begin{equation}
!  {\bf J}^{-1(m)} = -\alpha{\bf 1} + \sum_{i=2}^{m}|u^{(i)}\rangle
!  \otimes \langle v^{(i)}|
!\label{eq:JQN}
!\end{equation}
!Here, we used diagonal form for ${\bf J}^{-1(1)}$.
!\begin{equation}
!  {\bf J}^{-1(1)} = -\alpha {\bf 1}
!\end{equation}
!We can easily obtain $u^{(m)}$ and $v^{(m)}$ for the Broyden's 2nd
!method formula (eq.(\ref{eq:broyden2})).
!\begin{equation}
!  u^{(m)} = \Delta \rho^{(m)} + \alpha \Delta {\cal F}^{(m)} -
!  \sum_{i=2}^{m-1}\langle v^{(i)}|\Delta {\cal F}^{(m)}\rangle u^{(i)}
!\end{equation}
!\begin{equation}
!  v^{(m)} = \frac{1}{\parallel \Delta {\cal F}^{(m)}\parallel^2}
!  \Delta {\cal F}^{(m)}
!\end{equation}
!
!    a     : $\alpha$
!    d     : $\rho$
!    dd    : $\Delta \rho$
!    F     : $\cal F$
!    dF    : $\Delta F$
!

! dF_l   = dF  : $\Delta {\cal F}^{(m)}$
! d0_l    = d(m) + aF(m) - sum(i=2,m-1)<u(i)|F(m)>u(i)
! u_l      = u(m)
!          : $ u^{(m)} = \Delta \rho^{(m)} + \alpha \Delta {\cal F}^{(m)} -
!  \sum_{i=2}^{m-1}\langle v^{(i)}|\Delta {\cal F}^{(m)}\rangle u^{(i)}$
! v        = v(m)
!          : $ v^{(m)} = \frac{1}{\parallel \Delta {\cal F}^{(m)}\parallel^2}$
! din    = d(input,(m-1))   : $ \rho^{in, (m-1)}$ (:first and last),
!         or F(m)             : ${\cal F}^{(m)}$ (: as a transient array)
! dout    = d(output,(m-1))  : $ \rho^{out,(m-1)}$
!
!  J(m) = a + sum_{i=2}^m |u(i)><v(i)|
!  u(m) = dd(m) + a dF(m) - sum_{i=2}^{m-1} <v(i)|dF(m)>|u(i)>
!  v(m) = dF(m)/||dF(m)||
!

  subroutine m_CD_mix_broyden1(rmx)
    real(DP),intent(in) :: rmx
    integer   :: iter,j,mxiter,icr,jcr
! === DEBUG by tkato 2011/11/24 ================================================
!   real(DP)  :: vdF(nspin)
    real(DP)  :: vdF(nspin_m)
! ==============================================================================
    integer   :: id_sname = -1
                                                  __TIMER_SUB_START(1107)
    call tstatc0_begin('m_CD_mix_broyden1 ',id_sname,1)

    if(previous_waymix /= BROYD1.or.force_dealloc) then
       call mix_dealloc_previous()
       call mix_broyden_allocate();    F_l => din
       force_dealloc = .false.
    end if

! ================================== modified by K. Tagami =============== 11.0
!! --> T. Yamasaki  03 Aug. 2009
!    allocate(rmxtrc(nspin_m))
!    if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
!       call alloc_chgqstore_recompose_chgq(rmx,rmxtrc) ! --> chgq_l, chgqo_l, rmxtrc
!    else
!       rmxtrc(1:nspin_m) = rmx
!    end if
!! <--
!
    allocate(rmxtrc(nspin_m))
    if ( noncol ) then
       rmxtrc(1:nspin_m) = rmx
       rmxtrc(2:nspin_m) = min( rmx *spin_density_mixfactor, rmx_max )
    else
       if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
          call alloc_chgqstore_recompose_chgq(rmx,rmxtrc) ! --> chgq_l, chgqo_l, rmxtrc
       else
          rmxtrc(1:nspin_m) = rmx
       end if
    end if
! ======================================================================== 11.0

! ================================== Modified by K. Tagami =============
!    allocate(c_p(ista_kngp:iend_kngp)); call precon_4_charge_mix(rmx,c_p)
! --> T. Yamasaki, 03rd Aug. 2009
!!$    allocate(c_p(ista_kngp:iend_kngp)); c_p = 0; call precon_4_charge_mix(rmx,c_p)
    allocate(c_p(ista_kngp:iend_kngp,nspin_m)); c_p = 0

! =======================================================================

! ============================= modified by K. Tagami =============== 11.0
!    call precon_4_charge_mix(rmxtrc,c_p)

    if ( noncol ) then
       call precon_4_charge_mix_noncl(rmxtrc,c_p)
    else
       call precon_4_charge_mix(rmxtrc,c_p)
    endif
! ==================================================================== 11.0


    iter = iter_from_reset()                 !-(m_CD)

    if((iter-istrbr+1) <= 1) then
       call simple_mix1(c_p)                 !-(m_CD)
       !   din=chgqo_l; dout=chgq_l; (din,dout,c_p)->chgq_l
    else
       call mix_broyden_alloc3   !-(m_CD) d0_l,u_l,v_l, and dd_l are allocated
       call dF_F_d0_u_v_and_dd   !-(c.h.)   dF_l, F_l, dd_l, initial u_l,v_l,d0_l

       call set_ncrspd_mxiter_etc(iter,iU,mxiter) !-(m_CD) ->mxiter,ncrspd
       icr = icrspd_is(iter)                 !-(m_CD) function
       do j = 2, mxiter
          jcr = ncrspd(j)
          call renew_u_br(jcr,icr) !-(m_CD) |u(m)> = |u(m)> - <v(j)|dF(m)>|u(j)>
          call renew_v(jcr)        !-(c.h.) |v(m)> = |v(m)> - <u(j)|dd(m)>|v(j)>
          call renew_d_br(jcr)        !-(m_CD) |d(m)> = |d(m)> - <v(j)|F(m)> |u(j)>
       enddo!j-loop

       urec_l(:,:,:,icr,iU) = u_l(:,:,:)
       call mult1s(v_l,dF_l,f_p,vdF)
       call store_to_urec2(v_l,vdF,icr,iV)  !-(m_CD) v_l->urec_l(icr,iV)

       call renew_d_last_br(c_p)            !-(m_CD)
                                            ! chgq_l(|d(m)>) = |d(m)>-<v(m)|F(m)>|u(m)>
       call mix_broyden_dealloc3()          !-(m_CD)
    endif

! =================================== modified by K. Tagami =========== 11.0
!! --> T. Yamasaki  03 Aug. 2009
!    if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
!       call compose_chgq_dealloc_chgqstore()
!    end if
!    deallocate(rmxtrc)
!! <--
!
    if ( .not. noncol ) then
       if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
          call compose_chgq_dealloc_chgqstore()
       end if
    endif
    deallocate(rmxtrc)
! ======================================================================== 11.0

    if(af /= 0)  then
       allocate(work(kgp,kimg))
       work = 0  ! === Added by K. Tagami ===
       call charge_average(ANTIFERRO,chgq_l)
       deallocate(work)
    endif


    deallocate(c_p)
    previous_waymix = BROYD1
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(1107)
  contains
    subroutine dF_F_d0_u_v_and_dd
      !   dF_l(=deltaF) = (rho - dout) - (rhoo - din)
      !   F_l  = rho - rhoo (=\cal F^{m}); u_l  = (rhoo - din) + c_p*dF_l;
      !   dd_l = rhoo - din
      !   d0_l = rhoo+c_p* F_l;              v_l = c_p* dd_l

      integer                      :: is,k,i
                                                  __TIMER_SUB_START(1110)
                                                  __TIMER_DO_START(1156)


      do is = 1, ndim_magmom, af+1  !      do is = 1, nspin, af+1  ! === modified by K. Tagami === 11.0
         do k = 1, kimg
            do i = ista_kgpm,iend_kgpm
!  Revised by T. Yamasaki, 2009/05/28 (Pointed out by Fukata-san (NEC))
!!$               dF_l(i,k,is) = (rho(i,k,is)-rhoo(i,k,is)) - (dout(i,k,is)-din(i,k,is))
!!$               d0_l(i,k,is) = rhoo(i,k,is)+ c_pm(i)*(rho(i,k,is) - rhoo(i,k,is))
!!$               dd_l(i,k,is) = rhoo(i,k,is) - din(i,k,is)
               dF_l(i,k,is) = (rho(i,k,is)-rhoo(i,k,is)) - (dout(i,k,is)-F_l(i,k,is))
               d0_l(i,k,is) = rhoo(i,k,is)+ c_pm(i,is)*(rho(i,k,is) - rhoo(i,k,is))
               dd_l(i,k,is) = rhoo(i,k,is) - F_l(i,k,is)
! ---
               u_l(i,k,is)  = c_pm(i,is)*dF_l(i,k,is) + dd_l(i,k,is)
               F_l(i,k,is)  = rho(i,k,is) - rhoo(i,k,is)
            end do
            if(mype == 0) u_l(1,k,is) = 0.d0
         end do
      end do
                                                  __TIMER_DO_STOP(1156)
                                                  __TIMER_DO_START(1157)

      do is = 1, ndim_magmom, af+1 !   do is = 1, nspin, af+1 ! === modified by K. Tagami === 11.0
         do k = 1, kimg
            v_l(:,k,is) = c_pm(:,is)*dd_l(:,k,is)
         end do
      end do
                                                  __TIMER_DO_STOP(1157)
                                                  __TIMER_SUB_STOP(1110)
    end subroutine dF_F_d0_u_v_and_dd

    subroutine renew_v(j)
      integer, intent(in) :: j

      real(DP)  :: u_dd(nspin_m)  !     real(DP)  :: u_dd(nspin) ! === DEBUG by tkato 2011/11/24 ===
                                                  __TIMER_SUB_START(1117)
#ifdef _CDMIX_USE_POINTER_
      urec_l_3 => urec_l(ista_kgpm:iend_kgpm,1:kimg,1:nspin_m,j,iU)
      call mult1s(urec_l_3,dd_l,f_p,u_dd)!-(m_CD);<u(j)|dd(m)> ->u_dd
#else
      call mult1s5(urec_l,nbxmix,2,j,iU,dd_l,f_p,u_dd)!-(m_CD);<u(j)|dd(m)> ->u_dd
#endif
      call subtr_j_th_term(u_dd,iV,j,urec_l,v_l)  !-(m_CD)
      !                        |v(m)> = |v(m)> - <u(j)|dd(m)>|v(j)>
                                                  __TIMER_SUB_STOP(1117)
    end subroutine renew_v
  end subroutine m_CD_mix_broyden1

  subroutine m_CD_mix_broyden2(nfout,rmx)
    integer, intent(in) :: nfout
    real(DP),intent(in) :: rmx

    integer   :: iter,j,mxiter,icr,jcr
!!$    real(DP)  :: v_dF(nspin),vF(nspin)
    integer   :: id_sname = -1
! --> T. Yamasaki  03 Aug. 2009
    real(kind=DP), allocatable, dimension(:):: rmxtrc ! d(nspin_m)
!   real(kind=DP), allocatable, dimension(:,:,:) :: chgqstore_l, chgqostore_l
! <--
                                                  __TIMER_SUB_START(1124)
    call tstatc0_begin('m_CD_mix_broyden2 ',id_sname,1)

    if(previous_waymix /= BROYD2.or.force_dealloc) then
       call mix_dealloc_previous()
       call mix_broyden_allocate();    F_l => din
       force_dealloc = .false.
    end if

! ===================================== modified by K. Tagami ============ 11.0
!! --> T. Yamasaki  03 Aug. 2009
!    allocate(rmxtrc(nspin_m))
!    if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
!       call alloc_chgqstore_recompose_chgq(rmx,rmxtrc) ! --> chgq_l, chgqo_l, rmxtrc
!    else
!       rmxtrc(1:nspin_m) = rmx
!    end if
!! <--

    allocate(rmxtrc(nspin_m))
    if ( noncol ) then
       rmxtrc(1:nspin_m) = rmx
       rmxtrc(2:nspin_m) = min( rmx *spin_density_mixfactor, rmx_max )
    else
       if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
          call alloc_chgqstore_recompose_chgq(rmx,rmxtrc) ! --> chgq_l, chgqo_l, rmxtrc
       else
          rmxtrc(1:nspin_m) = rmx
       end if
    endif
! ======================================================================== 11.0

! ====================== Modified by K. Tagami =========
!    allocate(c_p(ista_kngp:iend_kngp)); call precon_4_charge_mix(rmx,c_p)
! --> T. Yamasaki, 03rd Aug. 2009
!!$    allocate(c_p(ista_kngp:iend_kngp)); c_p = 0; call precon_4_charge_mix(rmx,c_p)
    allocate(c_p(ista_kngp:iend_kngp,nspin_m)); c_p = 0
! =======================================================


! ============================== modified by K. Tagami ================== 11.0
!    call precon_4_charge_mix(rmxtrc,c_p)
!
    if ( noncol ) then
       call precon_4_charge_mix_noncl(rmxtrc,c_p)
    else
       call precon_4_charge_mix(rmxtrc,c_p)
    endif
! ======================================================================== 11.0

    iter = iter_from_reset()                 !-(m_CD)

    if((iter-istrbr+1) <= 1) then
       call simple_mix1(c_p)                 !-(m_CD)
       !   din=chgqo_l; dout=chgq_l; (din,dout,c_p)->chgq_l
    else
!!$       stop ' -- iter-istrbr+1 > 1 (m_CD_mix_broyden2) --'
       call mix_broyden_alloc2   !-(m_CD) d0_l,u_l, and v_l are allocated
       call dF_F_d0_u_and_v      !-(c.h.)   dF_l, F_l, initial u_l,v_l,d0_l

       call set_ncrspd_mxiter_etc(iter,iU,mxiter) !-(m_CD) ->mxiter,ncrspd
       !                  when hownew == RENEW: f,g,ncrspd, and urec_l are reset.
       icr = icrspd_is(iter)                 !-(m_CD) function
       do j = 2, mxiter
          jcr = ncrspd(j)
          call renew_u_br(jcr,icr) !-(m_CD) |u(m)> = |u(m)> - <v(j)|dF(m)>|u(j)>
          call renew_d_br(jcr)     !-(m_CD) |d(m)> = |d(m)> - <v(j)|F(m)> |u(j)>
       enddo!j-loop

       urec_l(:,:,:,icr,iU) = u_l(:,:,:)  ! storing
       urec_l(:,:,:,icr,iV) = v_l(:,:,:)  ! storing

       call renew_d_last_br(c_p)       !-(m_CD) chgq_l(|d(m)>) = |d(m)>-<v(m)|F(m)>|u(m)>
       call mix_broyden_dealloc2                      !-(m_CD)
    endif

! ============================== modified by K. Tagami ================= 11.0
!! --> T. Yamasaki  03 Aug. 2009
!    if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
!       call compose_chgq_dealloc_chgqstore()
!    end if
!    deallocate(rmxtrc)
!! <--
!
    if ( .not. noncol ) then
       if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
          call compose_chgq_dealloc_chgqstore()
       end if
    endif
    deallocate(rmxtrc)
! =========================================================================== 11.0

    if(af /= 0)  then
       allocate(work(kgp,kimg))
       work = 0 ! === Added by K. Tagami ===
       call charge_average(ANTIFERRO,chgq_l)
       deallocate(work)
    endif

    deallocate(c_p)
    previous_waymix = BROYD2
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(1124)
  contains
    subroutine dF_F_d0_u_and_v
      !   dF_l(=deltaF) = (rho - dout) - (rhoo - din)
      !   F_l = rho - rhoo (=\cal F^{m}); u_l  = (rhoo - din) + c_p*dF_l;
      !   d0_l = rhoo+c_p* F_l;              v_l = dF_l/( |dF_l| )

      integer                      :: is,k,i
      real(DP), dimension(nspin_m) :: fff
                                                  __TIMER_SUB_START(1125)
                                                  __TIMER_DO_START(1175)
      do is = 1, ndim_magmom, af+1  !    do is = 1, nspin, af+1  ==== modified by K. Tagami === 11.0
         do k = 1, kimg
            do i = ista_kgpm,iend_kgpm
!  Revised by T. Yamasaki, 2009/05/28 (Pointed out by Fukata-san (NEC))
!!$               dF_l(i,k,is) = (rho (i,k,is)-rhoo(i,k,is)) - (dout(i,k,is)-din(i,k,is))
!!$               d0_l(i,k,is) = rhoo(i,k,is) + c_pm(i)*(rho(i,k,is) - rhoo(i,k,is))
!!$               u_l(i,k,is)  = c_pm(i)*dF_l(i,k,is) + (rhoo(i,k,is) - din(i,k,is))
               dF_l(i,k,is) = (rho (i,k,is)-rhoo(i,k,is)) - (dout(i,k,is)-F_l(i,k,is))
               d0_l(i,k,is) = rhoo(i,k,is) + c_pm(i,is)*(rho(i,k,is) - rhoo(i,k,is))
               u_l(i,k,is)  = c_pm(i,is)*dF_l(i,k,is) + (rhoo(i,k,is) - F_l(i,k,is))
! ----
               F_l(i,k,is)  = rho(i,k,is) - rhoo(i,k,is)
            end do
            if(mype == 0) u_l(1,k,is) = 0.d0
         end do
      end do
                                                  __TIMER_DO_STOP(1175)

      call mult1s(dF_l,dF_l,f_p,fff)
      if(sum(fff) < 1.d-40)  call phase_error_with_msg(nfout, ' fmult is too small',__LINE__,__FILE__)

      if ( nspin_m == 2 .and. sw_mix_bothspins_sametime == YES ) then  !  if ( nspin == 2 .and. ...)  ! === DEBUG by tkato 2011/11/24 ===
        fff(1) = fff(1) + fff(2)
        fff(2) = fff(1)
      endif

! ========================= added by K. Tagami =========================== 11.0
      if ( noncol ) then
         fff(1) = sum( fff(:) )
         fff(:) = fff(1)
      endif
! ======================================================================== 11.0

                                                  __TIMER_DO_START(1176)
      do is = 1, ndim_magmom, af+1 !   do is = 1, nspin, af+1 ! === modified by K. Tagami === 11.0
         v_l(:,:,is) = dF_l(:,:,is)/fff(is)
      end do
                                                  __TIMER_DO_STOP(1176)
                                                  __TIMER_SUB_STOP(1125)
    end subroutine dF_F_d0_u_and_v

  end subroutine m_CD_mix_broyden2

! ===================== added by K. Tagami ============================== 5.0
  subroutine m_CD_mix_broyden2_with_hsr(nfout,rmx,mixocc)
    integer, intent(in) :: nfout
    real(DP),intent(in) :: rmx
    logical, intent(in) :: mixocc

    integer   :: iter,j,mxiter,icr,jcr
    integer   :: id_sname = -1
! --> T. Yamasaki  03 Aug. 2009
    real(kind=DP), allocatable, dimension(:):: rmxtrc ! d(nspin_m)
    real(kind=DP), allocatable, dimension(:,:,:) :: chgqstore_l, chgqostore_l
! <--
    call tstatc0_begin('m_CD_mix_broyden2_hsr ',id_sname,1)

    if (previous_waymix /= BROYD2.or.force_dealloc) then
       force_dealloc = .false.
       if ( first ) then
          if(mixocc) call set_i2lp_max2lp()
          call create_map_func(.true.,mixocc)
          call alloc_rho_hsr(mixocc)
          call create_map_func(.false.,mixocc)
          first = .false.
       endif
       call mix_dealloc_previous()
       call mix_dealloc_previous_hsr()
       call mix_broyden_allocate();        F_l => din
       call mix_broyden_allocate_hsr();    FF_hsr => din_hsr
    end if

! ==================================== modified by K. Tagami ========== 11.0
!    call map_hsr_to_rho( hsr, rho_hsr )
!    call map_hsr_to_rho( hsro,rhoo_hsr )
!
    if ( noncol ) then
       call map_hsr_to_rho_noncl( hsr, hsi, rho_hsr )
       call map_hsr_to_rho_noncl( hsro,hsio,rhoo_hsr )
       if(mixocc)then
          call map_om_to_rho_noncl( om, om_aimag, rho_hsr )
          call map_om_to_rho_noncl( omold, omold_aimag, rhoo_hsr )
       endif
    else
       call map_hsr_to_rho( hsr, rho_hsr )
       call map_hsr_to_rho( hsro,rhoo_hsr )
       if(mixocc)then
          call map_om_to_rho( om, rho_hsr )
          call map_om_to_rho( omold,rhoo_hsr )
       endif
    endif
! ========================================================================= 11.0

! ==================================== modified by K. Tagami ============= 11.0
!! --> T. Yamasaki  03 Aug. 2009
!    allocate(rmxtrc(nspin_m))
!    if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
!       call alloc_chgqstore_recompose_chgq(rmx,rmxtrc) ! --> chgq_l, chgqo_l, rmxtrc
!    else
!       rmxtrc(1:nspin_m) = rmx
!    end if
!! <--
!
    allocate(rmxtrc(nspin_m))
    if ( noncol ) then
       rmxtrc(1:nspin_m) = rmx
       rmxtrc(2:nspin_m) = min( rmx *spin_density_mixfactor, rmx_max )
    else
       if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
          call alloc_chgqstore_recompose_chgq(rmx,rmxtrc) ! --> chgq_l, chgqo_l, rmxtrc
       else
          rmxtrc(1:nspin_m) = rmx
       end if
    endif
! ========================================================================= 11.0

! ========================= modified by K. Tagami ======================= 11.0
!    if ( sw_recomposing_hsr == YES .and. af == 0 .and. nspin == 2 ) then
!       call alloc_rhostore_recomp( rmx, rmxtrc )
!    else
!       rmxtrc = rmx
!    endif

    if ( noncol ) then
       rmxtrc = rmx
       rmxtrc(2:nspin_m) = min( rmx *spin_density_mixfactor, rmx_max )
    else
       if ( sw_recomposing_hsr == YES .and. af == 0 .and. nspin == 2 ) then
          call alloc_rhostore_recomp( rmx, rmxtrc )
       else
          rmxtrc = rmx
       endif
    endif
! ========================================================================== 11.0

    allocate(c_p(ista_kngp:iend_kngp,nspin_m)); c_p = 0

! ============================= modiifed by K. Tagami =================== 11.0
!    call precon_4_charge_mix(rmxtrc,c_p)
!
    if ( noncol ) then
       call precon_4_charge_mix_noncl(rmxtrc,c_p)
    else
       call precon_4_charge_mix(rmxtrc,c_p)
    endif
! ====================================================================== 11.0

    iter = iter_from_reset()                 !-(m_CD)

    if((iter-istrbr+1) <= 1) then
       call simple_mix1(c_p)                 !-(m_CD)
       !   din=chgqo_l; dout=chgq_l; (din,dout,c_p)->chgq_l
       call simple_mix_kt( rmxtrc )
    else
!!$       stop ' -- iter-istrbr+1 > 1 (m_CD_mix_broyden2) --'

       call mix_broyden_alloc2   !-(m_CD) d0_l,u_l, and v_l are allocated
       call mix_broyden_alloc2_hsr

       call dF_F_d0_u_and_v_with_hsr

       call set_ncrspd_mxiter_etc(iter,iU,mxiter) !-(m_CD) ->mxiter,ncrspd
       !                  when hownew == RENEW: f,g,ncrspd, and urec_l are reset.
       icr = icrspd_is(iter)                 !-(m_CD) function
       do j = 2, mxiter
          jcr = ncrspd(j)
          call renew_u_br_with_hsr(nfout,jcr,icr) !-(m_CD) |u(m)> = |u(m)> - <v(j)|dF(m)>|u(j)>
          call renew_d_br_with_hsr(nfout,jcr)     !-(m_CD) |d(m)> = |d(m)> - <v(j)|F(m)> |u(j)>
       enddo!j-loop

       urec_l(:,:,:,icr,iU) = u_l(:,:,:)  ! storing
       urec_l(:,:,:,icr,iV) = v_l(:,:,:)  ! storing

       urec_hsr(:,:,icr,iU) = u_hsr(:,:)  ! storing
       urec_hsr(:,:,icr,iV) = v_hsr(:,:)  ! storing

       call renew_d_last_br_with_hsr(nfout, c_p, rmxtrc )   ! u_l, v_l, u_hsr, v_hsr --> rho, rho_hsr using vF(nspin_m)
                                !-(m_CD) chgq_l(|d(m)>) = |d(m)>-<v(m)|F(m)>|u(m)>

       call mix_broyden_dealloc2                      !-(m_CD)
       call mix_broyden_dealloc2_hsr
    endif

! ============================== modified by K. Tagami ================= 11.0
!! --> T. Yamasaki  03 Aug. 2009
!    if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
!       call compose_chgq_dealloc_chgqstore()
!    end if
!! <--
!    if ( sw_recomposing_hsr == YES .and. af == 0 .and. nspin == 2 ) then
!       call compose_rho_dealloc_store
!    end if
!    call map_rho_to_hsr( hsr, rho_hsr )
!    deallocate(rmxtrc)
!
    if ( .not. noncol ) then
       if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
          call compose_chgq_dealloc_chgqstore()
       end if
       if ( sw_recomposing_hsr == YES .and. af == 0 .and. nspin == 2 ) then
          call compose_rho_dealloc_store
       end if
    endif

    if ( noncol ) then
       call map_rho_to_hsr_noncl( hsr, hsi, rho_hsr )
       if ( mixocc ) call map_rho_to_om_noncl( om, om_aimag, rho_hsr )
    else
       call map_rho_to_hsr( hsr, rho_hsr )
       if(mixocc) call map_rho_to_om( om, rho_hsr )
    endif

    deallocate(rmxtrc)
! =========================================================================== 11.0

    if(af /= 0)  then
       allocate(work(kgp,kimg)); work = 0.0d0
       call charge_average(ANTIFERRO,chgq_l)
       deallocate(work)
    endif

    deallocate(c_p)

    previous_waymix = BROYD2
    call tstatc0_end(id_sname)
  contains

    subroutine dF_F_d0_u_and_v_with_hsr
      !   dF_l(=deltaF) = (rho - dout) - (rhoo - din)
      !   F_l = rho - rhoo (=\cal F^{m}); u_l  = (rhoo - din) + c_p*dF_l;
      !   d0_l = rhoo+c_p* F_l;              v_l = dF_l/( |dF_l| )

      integer                      :: is,k,i
      real(DP), dimension(nspin_m) :: fff

      do is = 1, ndim_magmom, af+1  !  do is = 1, nspin, af+1 ! === modified by K. Tagami === 11.0
         do k = 1, kimg
            do i = ista_kgpm,iend_kgpm
!  Revised by T. Yamasaki, 2009/05/28 (Pointed out by Fukata-san (NEC))
!!$               dF_l(i,k,is) = (rho (i,k,is)-rhoo(i,k,is)) - (dout(i,k,is)-din(i,k,is))
!!$               d0_l(i,k,is) = rhoo(i,k,is) + c_pm(i)*(rho(i,k,is) - rhoo(i,k,is))
!!$               u_l(i,k,is)  = c_pm(i)*dF_l(i,k,is) + (rhoo(i,k,is) - din(i,k,is))
               dF_l(i,k,is) = (rho (i,k,is)-rhoo(i,k,is)) - (dout(i,k,is)-F_l(i,k,is))
               d0_l(i,k,is) = rhoo(i,k,is) + c_pm(i,is)*(rho(i,k,is) - rhoo(i,k,is))
               u_l(i,k,is)  = c_pm(i,is)*dF_l(i,k,is) + (rhoo(i,k,is) - F_l(i,k,is))
! ----
               F_l(i,k,is)  = rho(i,k,is) - rhoo(i,k,is)
            end do
            if(mype == 0) u_l(1,k,is) = 0.d0
         end do
      end do

      do is = 1, ndim_magmom, af+1  !   do is = 1, nspin, af+1  ! === modified by K. Tagami === 11.0
         dF_hsr(:,is) = ( rho_hsr(:,is)-rhoo_hsr(:,is)) - ( dout_hsr(:,is)-FF_hsr(:,is))
         d0_hsr(:,is) = rhoo_hsr(:,is) + rmxtrc(is) *( rho_hsr(:,is) - rhoo_hsr(:,is))
          u_hsr(:,is) = rmxtrc(is) *dF_hsr(:,is) + ( rhoo_hsr(:,is) - FF_hsr(:,is) )
         FF_hsr(:,is) = rho_hsr(:,is) - rhoo_hsr(:,is)
      end do

      call mult1s(dF_l,dF_l,f_p,fff)

      do is = 1, ndim_magmom, af+1  !      do is = 1, nspin, af+1 ! === modified by K. Tagami === 11.0
         fff(is) = fff(is) + sum( dF_hsr(:,is)*dF_hsr(:,is) )
      end do

      if(sum(fff) < 1.d-40)  call phase_error_with_msg(nfout,' fmult is too small',__LINE__,__FILE__)

      if ( nspin_m == 2 .and. sw_mix_bothspins_sametime == YES ) then !  if ( nspin == 2 .and. ...) then ! === DEBUG by tkato 2011/11/24 ===
	fff(1) = fff(1) + fff(2)
        fff(2) = fff(1)
      endif

      if ( noncol ) then     ! === added by K. Tagami === 11.0
         fff(1) = sum( fff(:) )
         fff(:) = fff(1)
      endif

      do is = 1, ndim_magmom, af+1 !      do is = 1, nspin, af+1 ! === modified by K. Tagami === 11.0
         v_l(:,:,is) = dF_l(:,:,is)/fff(is)
      end do

      do is = 1, ndim_magmom, af+1  !      do is = 1, nspin, af+1  ! === modified by K. Tagami === 11.0
         do i=1,nsize_rho_hsr
            v_hsr(i,is) = dF_hsr(i,is)/fff(is)
         end do
      end do

    end subroutine dF_F_d0_u_and_v_with_hsr

  end subroutine m_CD_mix_broyden2_with_hsr
! ==================================================================== 5.0

  subroutine m_CD_mix_DFP(rmx)
    real(DP),intent(in) :: rmx
    integer   :: iter,j,mxiter,icr,jcr
    real(DP), pointer, dimension(:,:,:) :: F_l

    real(DP)  :: udF(nspin_m),wdF(nspin_m)  ! === DEBUG by tkato 2011/11/24 ===
! --> T. Yamasaki  03 Aug. 2009
    real(kind=DP), allocatable, dimension(:):: rmxtrc ! d(nspin_m)
!   real(kind=DP), allocatable, dimension(:,:,:) :: chgqstore_l, chgqostore_l
! <--
    integer   :: id_sname = -1
                                                  __TIMER_SUB_START(1126)
    call tstatc0_begin('m_CD_mix_DFP ',id_sname,1)

    if(previous_waymix /= DFP.or.force_dealloc) then
       force_dealloc = .false.
       call mix_dealloc_previous()
       call mix_DFP_allocate();    F_l => din
    end if

! ==================================== modified by K. Tagami ============= 11.0
!! --> T. Yamasaki  03 Aug. 2009
!    allocate(rmxtrc(nspin_m))
!    if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
!       call alloc_chgqstore_recompose_chgq(rmx,rmxtrc) ! --> chgq_l, chgqo_l, rmxtrc
!    else
!       rmxtrc(1:nspin_m) = rmx
!    end if
!! <--
!
    allocate(rmxtrc(nspin_m))
    if ( noncol ) then
       rmxtrc(1:nspin_m) = rmx
       rmxtrc(2:nspin_m) = min( rmx *spin_density_mixfactor, rmx_max )
    else
       if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
          call alloc_chgqstore_recompose_chgq(rmx,rmxtrc) ! --> chgq_l, chgqo_l, rmxtrc
       else
          rmxtrc(1:nspin_m) = rmx
       end if
    endif
! ========================================================================= 11.0

! ============================= Modified by K. Tagami ================
!    allocate(c_p(kgp)); call precon_4_charge_mix(rmx,c_p)
! --> T. Yamasaki, 03rd Aug. 2009
!!$    allocate(c_p(ista_kngp:iend_kngp)); c_p = 0; call precon_4_charge_mix(rmx,c_p)
    allocate(c_p(ista_kngp:iend_kngp,nspin_m)); c_p = 0.0d0
!!$    allocate(c_p(kgp)); c_p = 0; call precon_4_charge_mix(rmx,c_p)
! ====================================================================


! ============================= modiifed by K. Tagami =================== 11.0
!    call precon_4_charge_mix(rmxtrc,c_p)
!
    if ( noncol ) then
       call precon_4_charge_mix_noncl(rmxtrc,c_p)
    else
       call precon_4_charge_mix(rmxtrc,c_p)
    endif
! ====================================================================== 11.0

    iter = iter_from_reset()                 !-(m_CD)

    if((iter-istrbr+1) <= 1) then
       call simple_mix1(c_p)                 !-(m_CD)
       !   din=chgqo_l; dout=chgq_l; (din,dout,c_p)->chgq_l
    else
       call mix_DFP_alloc2   !-(m_CD) d0_l,u_l, and w_l are allocated

       call dF_F_d0_u_and_w  !-(c.h.)   dF_l, F_l, initial u_l,w_l,d0_l
       call set_ncrspd_mxiter_etc(iter,iU,mxiter) !-(m_CD) ->mxiter,ncrspd
       icr = icrspd_is(iter)                   !-(m_CD) function
       do j = 2, mxiter
          jcr = ncrspd(j)
          call renew_w(jcr)  !-(c.h.)    m      m      j   m   j      j   m   j
          !                            |w > = |w > - <v |dF >|u > - <y |dF >|w >
          call renew_d(jcr)  !-(c.h.)    m      m      j  m   j      j  m   j
          !                            |d > = |d > - <v |F >|u > - <y |F >|w >
       enddo!j-loop

       call mult1s(dF_l,u_l,f_p,udF)                 !  ->udF = <u(m)|dF(m)>
       call mult1s(dF_l,w_l,f_p,wdF)                 !  ->wdF = <w(m)|dF(m)>
       uuf(icr,1:nspin_m,iU) = udF(1:nspin_m)
       uuf(icr,1:nspin_m,iW) = wdF(1:nspin_m)

       call renew_d_last(udF,wdF) !-(c.h.)    m       m    m  m   m    m  m   m
       !                             chgq_l(|d >) = |d >-<v |F >|u >-<y |F >|w >
       call store_to_urec2(u_l,udF,icr,iVD)!-(m_CD) |v>=|u>/<u|dF> ->urec_l(icr,iVD)
       call store_to_urec2(w_l,wdF,icr,iY) !-(m_CD) |y>=|w>/<w|dF> ->urec_l(icr,iY)

       call mix_DFP_dealloc2                    !-(m_CD)
    endif

! ================================== modified by K. Tagami ============== 11.0
!! --> T. Yamasaki  03 Aug. 2009
!    if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) &
!         & call compose_chgq_dealloc_chgqstore()
!    deallocate(rmxtrc)
!! <--

    if ( .not. noncol ) then
       if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) &
            & call compose_chgq_dealloc_chgqstore()
    endif
    deallocate(rmxtrc)
! ========================================================================= 11.0

    if(af /= 0)  then
       allocate(work(kgp,kimg))
       work = 0   ! === Added by K. Tagami ===
       call charge_average(ANTIFERRO,chgq_l)
       deallocate(work)
    endif

    deallocate(c_p)
    previous_waymix = DFP
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(1126)
  contains
    subroutine dF_F_d0_u_and_w
      !   dF_l(=deltaF) = (rho - dout) - (rhoo - din)
      !   u_l  = rhoo-din;       F_l = rho - rhoo (= \cal F^{m})
      !   d0_l = rhoo+c_pm* F_l;   w_l = c_pm*dF_l;

      integer                      :: is,k,i
                                                  __TIMER_SUB_START(1127)
                                                  __TIMER_DO_START(1177)

      do is = 1, ndim_magmom, af+1  !      do is = 1, nspin, af+1  ! === modified by K. Tagami === 11.0
         do k = 1, kimg
            do i = ista_kgpm,iend_kgpm
!  Revised by T. Yamasaki, 2009/05/28 (Pointed out by Fukata-san (NEC))
               dF_l(i,k,is) = (rho(i,k,is)-rhoo(i,k,is)) - (dout(i,k,is)-F_l(i,k,is))
               d0_l(i,k,is) = rhoo(i,k,is) + c_pm(i,is)*(rho(i,k,is) - rhoo(i,k,is))
               u_l(i,k,is)  = rhoo(i,k,is) - F_l(i,k,is)
! ---
               w_l(i,k,is)  = c_pm(i,is)*dF_l(i,k,is)
               F_l(i,k,is)  = rho(i,k,is) - rhoo(i,k,is)
            end do
            if(mype == 0)  u_l(1,k,is) = 0.d0
         end do
      end do
                                                  __TIMER_DO_STOP(1177)
                                                  __TIMER_SUB_STOP(1127)
    end subroutine dF_F_d0_u_and_w

    subroutine renew_w(j)
      integer   :: j

      real(DP)  :: y_dF(nspin_m),v_dF(nspin_m) !     real(DP)  :: y_dF(nspin),v_dF(nspin)  ! === DEBUG by tkato 2011/11/24 ===
                                                  __TIMER_SUB_START(1128)

#ifdef _CDMIX_USE_POINTER_
      urec_l_3 => urec_l(ista_kgpm:iend_kgpm,1:kimg,1:nspin_m,j,iY)
      call mult1s(urec_l_3,dF_l,f_p,y_dF) !-(m_CD) <y(j)|dF(m)> ->y_dF
#else
      call mult1s5(urec_l,nbxmix,2,j,iY,dF_l,f_p,y_dF)
#endif
      y_dF = uuf(j,1:nspin_m,iW)*y_dF      ! = <w(j)|dF(j)><y(j)|dF(m)>
      !                                         <y(j)|dF(j)>*y(j) = w(j)
      call subtr_j_th_term(y_dF,iY,j,urec_l,w_l)  ! |w(m)> = |w(m)> - <y(j)|dF(m)>|w(j)>

#ifdef _CDMIX_USE_POINTER_
      urec_l_3 => urec_l(ista_kgpm:iend_kgpm,1:kimg,1:nspin_m,j,iVD)
      call mult1s(urec_l_3,dF_l,f_p,v_dF) ! <v(j)|dF(m)> ->v_dF
#else
      call mult1s5(urec_l,nbxmix,2,j,iVD,dF_l,f_p,v_dF)
#endif
      v_dF = uuf(j,1:nspin_m,iU)*v_dF  ! = <u(j)|dF(j)><v(j)|dF(m)>
      !                                         <u(j)|dF(j)>*v(j) = u(j)
      call subtr_j_th_term(v_dF,iVD,j,urec_l,w_l)  ! |w(m)> = |w(m)> - <v(j)|dF(m)>|u(j)>
                                                  __TIMER_SUB_STOP(1128)
    end subroutine renew_w

    subroutine renew_d(j)
      integer   :: j
      real(DP)  :: yF(nspin_m),vF(nspin_m)  !  real(DP)  :: yF(nspin),vF(nspin) ! === DEBUG by tkato 2011/11/24 ===
                                                  __TIMER_SUB_START(1129)

#ifdef _CDMIX_USE_POINTER_
      urec_l_3 => urec_l(ista_kgpm:iend_kgpm,1:kimg,1:nspin_m,j,iY)
      call mult1s(urec_l_3,F_l,f_p,yF) !-(m_CD);<y(j)|F(m)>  ->yF
#else
      call mult1s5(urec_l,nbxmix,2,j,iY,F_l,f_p,yF) !-(m_CD);<y(j)|F(m)>  ->yF
#endif
      yF = uuf(j,1:nspin_m,iW)*yf        ! = (y(j)|F(m)><w(j)|dF(j)>
      call subtr_j_th_term(yF,iY,j,urec_l,d0_l) ! |d(m)> = |d(m)> - <y(j)|F(m)>|w(j)>

#ifdef _CDMIX_USE_POINTER_
      urec_l_3 => urec_l(ista_kgpm:iend_kgpm,1:kimg,1:nspin_m,j,iVD)
      call mult1s(urec_l_3,F_l,f_p,vF)  ! <v(j)|F(m)>
#else
      call mult1s5(urec_l,nbxmix,2,j,iVD,F_l,f_p,vF)  ! <v(j)|F(m)>
#endif

      vF = uuf(j,1:nspin_m,iU)*vF
      call subtr_j_th_term(vF,iVD,j,urec_l,d0_l)! |d(m)> = |d(m)> - <v(j)|F(m)>|u(j)>
                                                  __TIMER_SUB_STOP(1129)
    end subroutine renew_d

    subroutine renew_d_last(udF,wdF)
      real(DP),intent(in)  :: udF(nspin_m),wdF(nspin_m) !     real(DP),intent(in)  :: udF(nspin),wdF(nspin) ! === DEBUG by tkato 2011/11/24 ===

      integer              :: is
      real(DP)             :: uF(nspin_m),wF(nspin_m)  !     real(DP)  :: uF(nspin),wF(nspin) ! === DEBUG by tkato 2011/11/24 ===

                                                  __TIMER_SUB_START(1130)
      call mult1s(F_l,u_l,f_p,uF)                   ! ->uF = <u(m)|F(m)>
      uF = uF/udF

      call mult1s(F_l,w_l,f_p,wF)                   ! ->wF = <w(m)|F(m)>
      wF = wF/wdF

      din  = rhoo
      dout = rho
                                                  __TIMER_DO_START(1178)

      do is = 1, ndim_magmom, af+1  !   do is = 1, nspin, af+1  ! === modified by K. Tagami === 11.0
         rho(:,:,is)  = d0_l(:,:,is)-uF(is)*u_l(:,:,is)-wF(is)*w_l(:,:,is)
      enddo
                                                  __TIMER_DO_STOP(1178)

      if(kgpm < kgp) then
         call concentrate_d_to_chg(rho,chgq_l) ! -(m_C.D.)
         call simple_mix_large_Gc(c_p,chgqo_l,chgq_l)     ! -(m_C.D.) chgq,chgqo,c_p -> chgq
      end if
                                                  __TIMER_SUB_STOP(1130)
    end subroutine renew_d_last

  end subroutine m_CD_mix_DFP

  subroutine mix_pulay_allocate

! =============================== modified by K. Tagami ========== 11.0
!    nspin_m  = nspin/(af+1)
!
    if ( noncol ) then
       nspin_m  = ndim_magmom
    else
       nspin_m  = nspin/(af+1)
    endif
! ================================================================= 11.0

!   allocate(f_p(ista_kgpm:iend_kgpm));          call precon_4_mult(f_p) !-(m_CD)
    allocate(f_p(ista_kgpm:iend_kgpm)); f_p = 0; call precon_4_mult(f_p) !-(m_CD), === Modified by K. Tagami ===

    allocate(din(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(dout(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(urec_l(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix,2))
    allocate(uuf_p(nbxmix,nspin_m))
    allocate(f(nbxmix,nbxmix,nspin_m))
    allocate(g_p(nbxmix,nspin_m))
    allocate(prj_wk(mp_kngp*kimg*nspin_m))
    allocate(ncrspd(nbxmix))

    if(sw_control_stepsize==ON) allocate(d0_l_h(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix))

    if(sw_gradient_simplex==ON)then
       allocate(rhoj(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix));rhoj=0.d0
       allocate(Frhoj(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix));Frhoj=0.d0
       allocate(rhojo(ista_kgpm:iend_kgpm,kimg,nspin_m));rhojo=0.d0
       allocate(Frhojo(ista_kgpm:iend_kgpm,kimg,nspin_m));Frhojo=0.d0
    endif

    allocate(ynorm(nbxmix,nspin_m));ynorm=1.d0
    allocate(sum12(nspin_m))
! ======================================= Added by K. Tagami ===========
    din = 0.0d0; dout = 0.0d0; urec_l = 0.0d0; uuf_p = 0.0d0; f = 0.0d0
    g_p = 0.0d0; prj_wk = 0.0d0; ncrspd = 0
! ======================================================================
  end subroutine mix_pulay_allocate

  subroutine mix_pulay_deallocate
    if(allocated(f_p)) deallocate(f_p)
    if(allocated(din)) deallocate(din)
    if(allocated(dout)) deallocate(dout)
    if(allocated(urec_l)) deallocate(urec_l)
    if(allocated(uuf_p)) deallocate(uuf_p)
    if(allocated(f)) deallocate(f)
    if(allocated(g_p)) deallocate(g_p)
    if(allocated(prj_wk)) deallocate(prj_wk)
    if(allocated(ncrspd)) deallocate(ncrspd)
    if(allocated(rhoj)) deallocate(rhoj)
    if(allocated(Frhoj)) deallocate(Frhoj)
    if(allocated(rhojo)) deallocate(rhojo)
    if(allocated(Frhojo)) deallocate(Frhojo)
    if(allocated(d0_l_h)) deallocate(d0_l_h)
    if (allocated(ynorm)) deallocate(ynorm)
    if(allocated(sum12)) deallocate(sum12)
  end subroutine mix_pulay_deallocate

  subroutine mix_pulay_alloc2
    allocate(d0_l(ista_kgpm:iend_kgpm,kimg,nspin_m))
    d0_l = 0.0d0  ! === Added by K. Tagami ===
    call alloc_rho_rhoo_and_cpm
  end subroutine mix_pulay_alloc2

  subroutine mix_pulay_dealloc2
    deallocate(d0_l)
    call dealloc_rho_rhoo_and_cpm
  end subroutine mix_pulay_dealloc2

  subroutine m_CD_mix_pulay(nfout,rmx)
    integer, parameter  :: iRho = 1, iResid = 2
    integer, intent(in) :: nfout
    real(DP),intent(in) :: rmx
    integer   :: iter, mxiter
    real(DP),pointer,dimension(:)  :: e_wk, f_wk, ww1, finv
    integer, pointer,dimension(:)  :: ip
! --> T. Yamasaki  03 Aug. 2009
    real(kind=DP), allocatable, dimension(:):: rmxtrc ! d(nspin_m)
!   real(kind=DP), allocatable, dimension(:,:,:) :: chgqstore_l, chgqostore_l
! <--
    real(kind=DP) :: rmxtt
    integer   :: id_sname = -1
                                                  __TIMER_SUB_START(1131)
    call tstatc0_begin('m_CD_mix_pulay ',id_sname,1)

    if(previous_waymix /= PULAY.or.force_dealloc) then
       force_dealloc = .false.
       call mix_dealloc_previous()
       call mix_pulay_allocate()
    end if

! ============================= modified by K. Tagami ================== 11.0
!! --> T. Yamasaki  03 Aug. 2009
!    allocate(rmxtrc(nspin_m))
!    if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
!       call alloc_chgqstore_recompose_chgq(rmx,rmxtrc) ! --> chgq_l, chgqo_l, rmxtrc
!    else
!       rmxtrc(1:nspin_m) = rmx
!    end if
!! <--
!
    allocate(rmxtrc(nspin_m))
    if ( noncol ) then
       rmxtrc(1:nspin_m) = rmx
       rmxtrc(2:nspin_m) = min( rmx *spin_density_mixfactor, rmx_max )
    else
       if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
          call alloc_chgqstore_recompose_chgq(rmx,rmxtrc) ! --> chgq_l, chgqo_l, rmxtrc
       else
          rmxtrc(1:nspin_m) = rmx
       end if
    endif
! ========================================================================= 11.0

    if(sw_control_stepsize==ON)then
       rmxtt = rmx*step_control_factor
       if(rmxtt>max_stepsize) rmxtt = max_stepsize
       rmxtrc = rmxtt
       call m_CtrlP_set_rmx(rmxtt)
       if(printable) write(nfout,'(a,f10.5)') 'step size for the current iteration : ',rmxtt
    endif

! ====================================== Modified by K. Tagami =========
!    allocate(c_p(ista_kngp:iend_kngp)); call precon_4_charge_mix(rmx,c_p)
    allocate(c_p(ista_kngp:iend_kngp,nspin_m)); c_p = 0

    if ( noncol ) then
       call precon_4_charge_mix_noncl(rmxtrc,c_p)
    else
       call precon_4_charge_mix(rmxtrc,c_p)
    endif
! ========================================================================= 11.0

    iter = iter_from_reset()                 !-(m_CD)
    if((iter-istrbr+1) <= 1) then
       call simple_mix1(c_p)                 !-(m_CD)
       !   din=chgqo_l; dout=chgq_l; (din,dout,c_p)->chgq_l
    else
       call mix_pulay_alloc2   !-(m_CD) d0_l,u_l, and w_l are allocated
       call set_ncrspd_mxiter(nbxmix,iter-istrbr,mxiter) ! -> ncrspd, mxiter
!!$       call mix_pulay_alloc3(nbxmix,iter-istrbr)   !-(c.h.) e_wk,f_wk,ww1,finv,ip
       call mix_pulay_alloc3(mxiter)   !-(c.h.) e_wk,f_wk,ww1,finv,ip

!!$       call Resid_and_dd_into_urec(iter-istrbr) !-(c.h.)
!!$       !                               dF ->urec_l; dd ->urec_l; d0_l,din,dout
!!$       call Ri_dot_Rj(iter-istrbr)          !-(c.h.) <R(i)|R(j)>->f
!!$       call get_finv(nbxmix,iter-istrbr,f)  !-(c.h.) f -> f^{-1}= <R(i)|R(j)>^{-1}
!!$
!!$       call Rj_dot_d(iter-istrbr)           !-(c.h.) <R(j)|d>,(j=1,iter-istrb) -> uuf_p
!!$
!!$       call get_gmatrix(iter-istrbr)        !-(c.h.) (f,uuf_p)->g
!!$       call renew_d_using_g(iter-istrbr,c_pm)     !-(c.h.)
       call Resid_and_dd_into_urec(mxiter) !-(c.h.)
       !                               dF ->urec_l; dd ->urec_l; d0_l,din,dout
       call Ri_dot_Rj(mxiter)          !-(c.h.) <R(i)|R(j)>->f
!!$       call get_finv(nbxmix,mxiter,f)  !-(c.h.) f -> f^{-1}= <R(i)|R(j)>^{-1}
       call get_finv_lapack(nbxmix,mxiter,f)  !-(c.h.) f -> f^{-1}= <R(i)|R(j)>^{-1}

       call Rj_dot_d(mxiter)           !-(c.h.) <R(j)|d>,(j=1,iter-istrb) -> uuf_p

       call get_gmatrix(mxiter)        !-(c.h.) (f,uuf_p)->g
       call renew_d_using_g(mxiter,c_pm)     !-(c.h.)

       call mix_pulay_dealloc3                    !-(c.h.)
       call mix_pulay_dealloc2                    !-(m_CD)
    endif

! ============================== modified by K. Tagami =================== 11.0
!! --> T. Yamasaki  03 Aug. 2009
!    if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) &
!         & call compose_chgq_dealloc_chgqstore()
!    deallocate(rmxtrc)
!! <--
!
    if ( .not. noncol ) then
       if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) &
            & call compose_chgq_dealloc_chgqstore()
    endif
    deallocate(rmxtrc)
! ========================================================================= 11.0

    if(af /= 0)  then
       allocate(work(kgp,kimg))
       work = 0 ! === Added by K. Tagami ===
       call charge_average(ANTIFERRO,chgq_l)
       deallocate(work)
    endif

    deallocate(c_p)
    previous_waymix = PULAY
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(1131)
  contains
    subroutine mix_pulay_alloc3(m)
      integer, intent(in) :: m
      allocate(e_wk(m*m)); allocate(f_wk(m*m)); allocate(ww1(m)); allocate(finv(m*m))
      allocate(ip(m))
! ===================================== Added by K. Tagami ============
      e_wk = 0; f_wk = 0; ww1 = 0; finv = 0; ip = 0
! =====================================================================
    end subroutine mix_pulay_alloc3

    subroutine set_ncrspd_mxiter(n,iter,m)
      integer, intent(in)  :: n, iter
      integer, intent(out) :: m
      integer :: i, nx
      if(hownew == ANEW) then
         m = iter
!!$         ncrspd(:) = (/(i,i=1,m)/)
         do i=1,iter
            ncrspd(i) = i
         end do
      else ! hownew == RENEW
         if(iter <= n) then
            m = iter
!!$            ncrspd(:) = (/(i,i=1,m)/)
            do i=1,iter
               ncrspd(i) = i
            end do
         else
            m = n
            nx = ncrspd(1)
            do i = 1, m-1
               ncrspd(i) = ncrspd(i+1)
            end do
            ncrspd(m) = nx
         end if
      end if
    end subroutine set_ncrspd_mxiter

    subroutine mix_pulay_dealloc3
      deallocate(e_wk); deallocate(f_wk); deallocate(ww1); deallocate(finv)
      deallocate(ip)
    end subroutine mix_pulay_dealloc3

    subroutine Resid_and_dd_into_urec(iter)
      integer, intent(in) :: iter
      integer             :: itc,itc0,itc1
      integer :: i,j,k,imix
      real(kind=DP) :: sum1,sum2
                                                  __TIMER_SUB_START(1132)
      itc = ncrspd(iter)
      if(sw_gradient_simplex==ON)then
         do imix=2,iter-1
            itc0 = ncrspd(imix)
            itc1 = ncrspd(imix-1)
            do i=1,nspin_m
               do j=1,kimg
                  do k=ista_kgpm,iend_kgpm
                     urec_l(k,j,i,itc0,iResid) = rho(k,j,i)-rhoo(k,j,i)-(Frhoj(k,j,i,itc1)-rhoj(k,j,i,itc1))
                     urec_l(k,j,i,itc0,iRho  ) = rhoo(k,j,i)-rhoj(k,j,i,itc1)
                  enddo
               enddo
            enddo
         enddo
         do i=1,nspin_m
            do j=1,kimg
               do k=ista_kgpm,iend_kgpm
                  urec_l(k,j,i,itc,iResid) = rho(k,j,i)-rhoo(k,j,i)-(dout(k,j,i)-din(k,j,i))
                  urec_l(k,j,i,itc,iRho  ) = rhoo(k,j,i)-din(k,j,i)
                  rhoj(k,j,i,itc) = rhoo(k,j,i)
                  Frhoj(k,j,i,itc) = rho(k,j,i)
                  d0_l(k,j,i) = rho(k,j,i) - rhoo(k,j,i)
                  din(k,j,i)  = rhoo(k,j,i)
                  dout(k,j,i) = rho(k,j,i)
               enddo
            enddo
         enddo
      else
         do i=1,nspin_m
            do j=1,kimg
               do k=ista_kgpm,iend_kgpm
                  urec_l(k,j,i,itc,iResid) = rho(k,j,i) - rhoo(k,j,i) - (dout(k,j,i) - din(k,j,i)) ! =dF(=delta F^i)
                  urec_l(k,j,i,itc,iRho  ) = rhoo(k,j,i) - din(k,j,i)                ! =dd
                  d0_l(k,j,i) = rho(k,j,i) - rhoo(k,j,i)
                  din(k,j,i)  = rhoo(k,j,i)
                  dout(k,j,i) = rho(k,j,i)
               enddo
            enddo
         enddo
      endif
!!$      ynorm(itc,:)=0.d0
      do i=1,nspin_m
         sum12(i) = 0.d0
         do j=1,kimg
            do k=ista_kgpm,iend_kgpm
!!$               ynorm(itc,i) = ynorm(itc,i)+urec_l(k,j,i,itc,iResid)*urec_l(k,j,i,itc,iResid)
               sum12(i) = sum12(i) + urec_l(k,j,i,itc,iResid)*urec_l(k,j,i,itc,iResid)
            enddo
         enddo
      enddo
!!$      call mpi_allreduce(MPI_IN_PLACE,ynorm(itc,1),nspin_m,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
!!$      ynorm(itc,:) = 1.d0/sqrt(univol*ynorm(itc,:))
      call mpi_allreduce(MPI_IN_PLACE,sum12,nspin_m,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
      ynorm(itc,:) = 1.d0/sqrt(univol*sum12(:))
      sum1 = 0.d0
      do i=1,nspin_m
         do j=1,kimg
            do k=ista_kgpm,iend_kgpm
               sum1 = sum1 + d0_l(k,j,i)*d0_l(k,j,i)
            enddo
         enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,sum1,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
      sum1 = dsqrt(sum1)
      if(sum1<alpha_pulay_damp_thres) then
         alpha_pulay = alpha_pulay*alpha_pulay_damp
         if(alpha_pulay>0.and.printable) &
         & write(nfout,'(a,e12.5)') '!** new value for the parameter alpha : ',alpha_pulay
      endif
      if(sw_control_stepsize==ON) d0_l_h(:,:,:,itc) = d0_l(:,:,:)
      if(iter>=2.and.sw_control_stepsize==ON)then
         sum2 = 0.d0
         itc1 = ncrspd(iter-1)
         do i=1,nspin_m
            do j=1,kimg
               do k=ista_kgpm,iend_kgpm
                  sum2 = sum2 + d0_l_h(k,j,i,itc1)*d0_l_h(k,j,i,itc1)
               enddo
            enddo
         enddo
         call mpi_allreduce(MPI_IN_PLACE,sum2,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
         sum2 = dsqrt(sum2)
         step_control_factor = max(0.5d0,min(2.0d0,sum2/sum1))
      endif
                                                  __TIMER_SUB_STOP(1132)
    end subroutine Resid_and_dd_into_urec

    subroutine Ri_dot_Rj(n)
      integer, intent(in) :: n
      integer  :: it,jt,itc,jtc
! === DEBUG by tkato 2011/11/24 ================================================
!     real(DP) :: ff1(nspin),ff1tmp
      real(DP) :: ff1(nspin_m),ff1tmp
! ==============================================================================
                                                  __TIMER_SUB_START(1133)
                                                  __TIMER_DO_START(1179)
      do it = 1, n
         itc = ncrspd(it)
         do jt = it, n
            jtc = ncrspd(jt)
#ifdef _CDMIX_USE_POINTER_
            urec_l_3 => urec_l(ista_kgpm:iend_kgpm,1:kimg,1:nspin_m,itc,iResid)
            urec_l_3_2 => urec_l(ista_kgpm:iend_kgpm,1:kimg,1:nspin_m,jtc,iResid)
            call mult1s(urec_l_3,urec_l_3_2,f_p,ff1)   ! <delta F^i|delta F^j>
#else
            if(sw_mix_bothspins_sametime==ON .and. nspin_m>1)then
               call mult1s10_reduce_spin(urec_l,nbxmix,2,itc,iResid,urec_l,jtc,iResid,f_p,ff1tmp)   ! <delta F^i|delta F^j>
               ff1(1)=ff1tmp;ff1(2)=ff1tmp
            else
               call mult1s10(urec_l,nbxmix,2,itc,iResid,urec_l,jtc,iResid,f_p,ff1)   ! <delta F^i|delta F^j>
            endif

! ============================= added by K. Tagami ======================= 11.0
            if ( noncol ) then
               call mult1s10_reduce_spin( urec_l, nbxmix, 2, itc, iResid, &
                    &                     urec_l, jtc, iResid, f_p, ff1tmp )
                                                        ! <delta F^i|delta F^j>
               ff1(:) = ff1tmp
            endif
! ======================================================================== 11.0

#endif
            f(it,jt,1:nspin_m) = ff1(1:nspin_m)
            if(jt /= it) f(jt,it,1:nspin_m) = f(it,jt,1:nspin_m)
         end do
      end do
                                                  __TIMER_DO_STOP(1179)
                                                  __TIMER_SUB_STOP(1133)
    end subroutine Ri_dot_Rj

    subroutine Rj_dot_d(n)
      integer, intent(in) :: n
      integer  :: jt, jtc
! === DEBUG by tkato 2011/11/24 ================================================
!     real(DP) :: ff1(nspin),ff1tmp
      real(DP) :: ff1(nspin_m),ff1tmp
! ==============================================================================
                                                  __TIMER_SUB_START(1138)
      do jt = 1, n
         jtc = ncrspd(jt)
#ifdef _CDMIX_USE_POINTER_
         urec_l_3 => urec_l(ista_kgpm:iend_kgpm,1:kimg,1:nspin_m,jtc,iResid)
         call mult1s(urec_l_3,d0_l,f_p,ff1)
#else
         if(sw_mix_bothspins_sametime==ON .and. nspin_m>1)then
            call mult1s5_reduce_spin(urec_l,nbxmix,2,jtc,iResid,d0_l,f_p,ff1tmp)
            ff1(1) = ff1tmp;ff1(2)=ff1tmp
         else
            call mult1s5(urec_l,nbxmix,2,jtc,iResid,d0_l,f_p,ff1)
         endif

! ============================= added by K. Tagami ======================= 11.0
         if ( noncol ) then
            call mult1s5_reduce_spin(urec_l,nbxmix,2,jtc,iResid,d0_l,f_p,ff1tmp)
            ff1(:) = ff1tmp
         endif
! ======================================================================== 11.0

#endif
         uuf_p(jt,1:nspin_m) = ff1(1:nspin_m)
      end do
                                                  __TIMER_SUB_STOP(1138)
    end subroutine Rj_dot_d

    subroutine get_finv_lapack(m,n,f)
      integer,intent(in)                             :: m,n
      real(DP),intent(inout),dimension(m,m,nspin_m) :: f
      real(DP), allocatable,dimension(:,:) :: fwork
      integer :: is,inf,it,jt,kt,nnspin
      real(DP) :: div,tmp
      allocate(fwork(n,n))
      nnspin = nspin
      if(sw_mix_bothspins_sametime==ON .or. af==1) nnspin=1

! ======================= added by K. Tagami ============= 11.0
      if ( noncol ) then
         nnspin = 1
      end if
! ======================================================== 11.0

      do is=1,nnspin
         if(ipripulay >= 2) then
            write(nfout,600) n,(('(',it,jt,')',f(it,jt,is),jt=1,n),it=1,n)
600         format(//11x,"**input matrix**"/12x &
                 & ,"horder=",I5/(2x,4(1x,1a,i2,",",i2,1a,e14.6)))
         end if
         fwork=0
         do it=1,n
            do jt=1,n
               fwork(jt,it) = f(jt,it,is)*ynorm(jt,is)*ynorm(it,is)
               if(it==jt) fwork(jt,it)=fwork(jt,it)+alpha_pulay
            enddo
         enddo
         call dpotrf('U',n,fwork,n,inf)
         call dpotri('U',n,fwork,n,inf)
         do it=1,n-1
            do jt=it+1,n
               fwork(jt,it) = fwork(it,jt)
            enddo
         enddo
         do it=1,n
            do jt=1,n
               f(jt,it,is) = fwork(jt,it)*ynorm(jt,is)*ynorm(it,is)
            enddo
         enddo
         if(ipripulay >= 2) then
            write(nfout,630) (('(',it,jt,')',f(it,jt,is),it=1,n),jt=1,n)
630         format(/11x, "**inverse matrix**" &
                 & ,/(2x,4(1x,1a,i2,",",i2,1a,e14.6)))
         end if
      enddo
      if(sw_mix_bothspins_sametime==ON .and. nspin_m>1)then
         do it=1,n
            do jt=1,n
               f(jt,it,2) = f(jt,it,1)
            enddo
         enddo
      endif
! ============================== added by K. Tagami ========== 11.0
      if ( noncol ) then
         do it=1,n
            do jt=1,n
               f(jt,it,:) = f(jt,it,1)
            enddo
         end do
      endif
! ============================================================ 11.0
      deallocate(fwork)

    end subroutine get_finv_lapack

    subroutine get_finv(m,n,f)
      integer,intent(in)                             :: m,n
      real(DP),intent(inout),dimension(m,m,nspin_m) :: f

      integer                        :: icount,is,jt,it,icon
      real(DP)                       :: div
                                                  __TIMER_SUB_START(1134)

      e_wk = 0.d0
      do it = 1, n
         e_wk(it*it) = 1.d0
      end do

      do is = 1, ndim_magmom, af+1  !      do is = 1, nspin, af+1  ! === modified by K. Tagami === 11.0
         div = 1.d0/f(1,1,is)
         icount = 1
                                                  __TIMER_DO_START(1180)
         do jt = 1, n
            do it = 1, n
               f_wk(icount) = f(it,jt,is)*div
               icount = icount + 1
            end do
         end do
                                                  __TIMER_DO_STOP(1180)
         if(ipripulay >= 1) then
            write(nfout,600) n,(('(',it,jt,')',f(it,jt,is)*div,jt=1,n),it=1,n)
600         format(//11x,"**input matrix**"/12x &
                 & ,"horder=",I5/(2x,4(1x,1a,i2,",",i2,1a,e14.6)))
         end if
         call rdecomp(n,f_wk,ww1,ip,icon)
         if(icon /= 0) then
            call phase_error_with_msg(nfout, 'LU decomposition is impossible.',__LINE__,__FILE__)
         else
            call rsolve(n,n,f_wk,e_wk,finv,ip)
         endif

         icount = 1
                                                  __TIMER_DO_START(1181)
         do jt = 1, n
            do it = 1, n
               f(it,jt,is) = finv(icount)
               icount = icount + 1
            end do
         end do
                                                  __TIMER_DO_STOP(1181)
         if(ipripulay >= 1) then
            write(nfout,630) (('(',it,jt,')',f(it,jt,is),it=1,n),jt=1,n)
630         format(/11x, "**inverse matrix**" &
                 & ,/(2x,4(1x,1a,i2,",",i2,1a,e14.6)))
         end if
      end do
                                                  __TIMER_SUB_STOP(1134)
    end subroutine get_finv

    subroutine get_gmatrix(n)
      integer,intent(in) :: n
      integer :: is, it, jt, nnspin
                                                  __TIMER_SUB_START(1139)
      nnspin = nspin
      if(sw_mix_bothspins_sametime==ON .or. af==1) nnspin=1

! ============================ added by K. Tagami ============= 11.0
      if ( noncol ) nnspin = 1
! ============================================================== 11.0

      g_p = 0.d0
      do is = 1, nnspin
                                                  __TIMER_DO_START(1188)
         do it = 1, n
            do jt = 1, n
               g_p(it,is) = g_p(it,is) - f(jt,it,is)*uuf_p(jt,is)
            end do
         end do
                                                  __TIMER_DO_STOP(1188)
         if(ipripulay >= 2) then
            write(nfout,'(" -- g_p(1:",i3,") --")') n
            write(nfout,'(8f20.12)') (g_p(it,is),it=1,n)
         end if
      end do
      if(sw_mix_bothspins_sametime==ON .and. nspin_m>1)then
         do it = 1,n
            g_p(it,2) = g_p(it,1)
         enddo
      endif
! ============================== added by K. Tagami ============ 11.0
      if ( noncol ) then
         do it = 1,n
            g_p(it,:) = g_p(it,1)
         enddo
      endif
! ============================================================== 11.0

                                                  __TIMER_SUB_STOP(1139)
    end subroutine get_gmatrix

    subroutine renew_d_using_g(n,p)
      integer, intent(in)                                :: n
      real(DP),intent(in),dimension(ista_kgpm:iend_kgpm,nspin_m) :: p
      integer    :: is, k, i, it, itc, ns
                                                  __TIMER_SUB_START(1140)

!!$      do is = 1, nspin, af+1
      ns = nspin_for_qnewton()
      do is = 1, ns,af+1
         do k = 1, kimg
                                                  __TIMER_DO_START(1189)
            do i = ista_kngp, iend_kngp
               rho(i,k,is)  = rhoo(i,k,is) + p(i,is)*d0_l(i,k,is)
            end do
                                                  __TIMER_DO_STOP(1189)
                                                  __TIMER_DO_STOP(1190)
            do it = 1, n
               itc = ncrspd(it)
               do i = ista_kngp, iend_kngp
                  rho(i,k,is) = rho(i,k,is) + g_p(it,is)* &
                       &        (urec_l(i,k,is,itc,iRho) + p(i,is)*urec_l(i,k,is,itc,iResid))
               end do
            end do
                                                  __TIMER_DO_STOP(1190)
         end do
      end do

! ============================== modified by K. Tagami ================ 11.0
!      if(sw_force_simple_mixing==ON .and. sw_recomposing==ON) call simple_mix2(c_p)
!
      if ( .not. noncol ) then
         if(sw_force_simple_mixing==ON .and. sw_recomposing==ON) then
            call simple_mix2(c_p)
         endif
      endif
! ===================================================================== 11.0

      if(kgpm < kgp) then
         call concentrate_d_to_chg(rho,chgq_l) !-(m_C.D.)
         call simple_mix_large_Gc(c_p,chgqo_l,chgq_l)     !-(m_C.D.) chgq,chgqo,c_p ->chgq
      end if
                                                  __TIMER_SUB_STOP(1140)
    end subroutine renew_d_using_g
  end subroutine m_CD_mix_pulay

!!$ 11.07 AS Pulay version of 'sw_mix_charge_hardpart'
  subroutine m_CD_mix_pulay_with_hsr(nfout,rmx,mixocc)
    integer, parameter  :: iRho = 1, iResid = 2
    integer, intent(in) :: nfout
    real(DP),intent(in) :: rmx
    logical, intent(in) :: mixocc
    integer   :: iter, mxiter
    real(DP),pointer,dimension(:)  :: e_wk, f_wk, ww1, finv
    integer, pointer,dimension(:)  :: ip
! --> T. Yamasaki  03 Aug. 2009
    real(kind=DP), allocatable, dimension(:):: rmxtrc ! d(nspin_m)
    real(kind=DP), allocatable, dimension(:,:,:) :: chgqstore_l, chgqostore_l
! <--
    real(kind=DP) :: rmxtt
    integer   :: id_sname = -1
    call tstatc0_begin('m_CD_mix_pulay_with_hsr ',id_sname,1)
    if(previous_waymix /= PULAY.or.force_dealloc) then
       force_dealloc = .false.
       if ( first ) then
          if(mixocc) call set_i2lp_max2lp()
          call create_map_func(.true.,mixocc)
          call alloc_rho_hsr(mixocc)
          call create_map_func(.false.,mixocc)
          first = .false.
       endif
       call mix_dealloc_previous()
       call mix_dealloc_previous_hsr()
       call mix_pulay_allocate()
       call mix_pulay_allocate_hsr()
    end if

! ==================================== modified by K. Tagami ========== 11.0
!    call map_hsr_to_rho( hsr, rho_hsr )
!    call map_hsr_to_rho( hsro,rhoo_hsr )
!
    if ( noncol ) then
       call map_hsr_to_rho_noncl( hsr, hsi, rho_hsr )
       call map_hsr_to_rho_noncl( hsro,hsio,rhoo_hsr )
       if(mixocc)then
          call map_om_to_rho_noncl( om, om_aimag, rho_hsr )
          call map_om_to_rho_noncl( omold, omold_aimag, rhoo_hsr )
       endif
    else
       call map_hsr_to_rho( hsr, rho_hsr )
       call map_hsr_to_rho( hsro,rhoo_hsr )
       if(mixocc)then
          call map_om_to_rho( om, rho_hsr )
          call map_om_to_rho( omold,rhoo_hsr )
       endif
    endif
! ========================================================================= 11.0

! ==================================== modified by K. Tagami ============= 11.0
!! --> T. Yamasaki  03 Aug. 2009
!    allocate(rmxtrc(nspin_m))
!    if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
!       call alloc_chgqstore_recompose_chgq(rmx,rmxtrc) ! --> chgq_l, chgqo_l, rmxtrc
!    else
!       rmxtrc(1:nspin_m) = rmx
!    end if
!! <--
!
    allocate(rmxtrc(nspin_m))
    if ( noncol ) then
       rmxtrc(1:nspin_m) = rmx
       rmxtrc(2:nspin_m) = min( rmx *spin_density_mixfactor, rmx_max )
    else
       if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
          call alloc_chgqstore_recompose_chgq(rmx,rmxtrc) ! --> chgq_l, chgqo_l, rmxtrc
       else
          rmxtrc(1:nspin_m) = rmx
       end if
    endif
! ========================================================================= 11.0

! ========================= modified by K. Tagami ======================= 11.0
!    if ( sw_recomposing_hsr == YES .and. af == 0 .and. nspin == 2 ) then
!       call alloc_rhostore_recomp( rmx, rmxtrc )
!    else
!       rmxtrc = rmx
!    endif

    if ( noncol ) then
       rmxtrc(1:nspin_m) = rmx
       rmxtrc(2:nspin_m) = min( rmx *spin_density_mixfactor, rmx_max )
    else
       if ( sw_recomposing_hsr == YES .and. af == 0 .and. nspin == 2 ) then
          call alloc_rhostore_recomp( rmx, rmxtrc )
       else
          rmxtrc(1:nspin_m) = rmx
       endif
    endif
! ========================================================================== 11.0

    if(sw_control_stepsize==ON)then
       rmxtt = rmx*step_control_factor
       if(rmxtt>max_stepsize) rmxtt = max_stepsize
       rmxtrc = rmxtt
       call m_CtrlP_set_rmx(rmxtt)
       if(printable) write(nfout,'(a,f10.5)') 'step size for the current iteration : ',rmxtt
    endif

! ====================================== Modified by K. Tagami =========
!    allocate(c_p(ista_kngp:iend_kngp)); call precon_4_charge_mix(rmx,c_p)
! --> T. Yamasaki, 03rd Aug. 2009
!!$    allocate(c_p(ista_kngp:iend_kngp)); c_p = 0; call precon_4_charge_mix(rmx,c_p)
    allocate(c_p(ista_kngp:iend_kngp,nspin_m)); c_p = 0
! ========================================================================


! ============================= modiifed by K. Tagami =================== 11.0
!    call precon_4_charge_mix(rmxtrc,c_p)
!
    if ( noncol ) then
       call precon_4_charge_mix_noncl(rmxtrc,c_p)
    else
       call precon_4_charge_mix(rmxtrc,c_p)
    endif
! ====================================================================== 11.0

    iter = iter_from_reset()                 !-(m_CD)
    if(iter.eq.1) then
       alpha_pulay = alpha_pulay_org
    endif
    if((iter-istrbr+1) <= 1) then
       call simple_mix1(c_p)                 !-(m_CD)
       call simple_mix_kt( ommix_factor*rmxtrc )
       !   din=chgqo_l; dout=chgq_l; (din,dout,c_p)->chgq_l
    else
       call mix_pulay_alloc2   !-(m_CD) d0_l,u_l, and w_l are allocated
       call mix_pulay_alloc2_hsr()
       call set_ncrspd_mxiter(nbxmix,iter-istrbr,mxiter) ! -> ncrspd, mxiter
!!$       call mix_pulay_alloc3(mxiter)   !-(c.h.) e_wk,f_wk,ww1,finv,ip

       call Resid_and_dd_into_urec_with_hsr(mxiter) !-(c.h.)
       call Ri_dot_Rj_with_hsr(mxiter)          !-(c.h.) <R(i)|R(j)>->f
       call get_finv_lapack_with_hsr(nbxmix,mxiter,f)  !-(c.h.) f -> f^{-1}= <R(i)|R(j)>^{-1}
       call Rj_dot_d_with_hsr(mxiter)           !-(c.h.) <R(j)|d>,(j=1,iter-istrb) -> uuf_p
       call get_gmatrix_with_hsr(mxiter)        !-(c.h.) (f,uuf_p)->g
       call renew_d_using_g_with_hsr(mxiter,c_pm)     !-(c.h.)

!!$       call mix_pulay_dealloc3                    !-(c.h.)
       call mix_pulay_dealloc2                    !-(m_CD)
       call mix_pulay_dealloc2_hsr                    !-(m_CD)
    endif

! ============================== modified by K. Tagami ================= 11.0
!! --> T. Yamasaki  03 Aug. 2009
!    if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
!       call compose_chgq_dealloc_chgqstore()
!    end if
!! <--
!    if ( sw_recomposing_hsr == YES .and. af == 0 .and. nspin == 2 ) then
!       call compose_rho_dealloc_store
!    end if
!    call map_rho_to_hsr( hsr, rho_hsr )
!    deallocate(rmxtrc)
!
    if ( .not. noncol ) then
       if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
          call compose_chgq_dealloc_chgqstore()
       end if
       if ( sw_recomposing_hsr == YES .and. af == 0 .and. nspin == 2 ) then
          call compose_rho_dealloc_store
       end if
    endif

    if ( noncol ) then
       call map_rho_to_hsr_noncl( hsr, hsi, rho_hsr )
       if(mixocc) call map_rho_to_om_noncl( om, om_aimag, rho_hsr )
    else
       call map_rho_to_hsr( hsr, rho_hsr )
       if(mixocc) call map_rho_to_om( om, rho_hsr )
    endif

    deallocate(rmxtrc)
! =========================================================================== 11.0

    if(af /= 0)  then
       allocate(work(kgp,kimg))
       work = 0 ! === Added by K. Tagami ===
       call charge_average(ANTIFERRO,chgq_l)
       deallocate(work)
    endif

    deallocate(c_p)
    previous_waymix = PULAY
    call tstatc0_end(id_sname)
  contains
    subroutine mix_pulay_alloc3(m)
      integer, intent(in) :: m
      allocate(e_wk(m*m)); allocate(f_wk(m*m)); allocate(ww1(m)); allocate(finv(m*m))
      allocate(ip(m))
      e_wk = 0; f_wk = 0; ww1 = 0; finv = 0; ip = 0  ! === Added by K. Tagami ===
    end subroutine mix_pulay_alloc3

    subroutine set_ncrspd_mxiter(n,iter,m)
      integer, intent(in)  :: n, iter
      integer, intent(out) :: m
      integer :: i, nx
      if(hownew == ANEW) then
         m = iter
!         ncrspd(:) = (/(i,i=1,m)/)
         do i=1,m
           ncrspd(i) = i
         enddo
      else ! hownew == RENEW
         if(iter <= n) then
            m = iter
!            ncrspd(:) = (/(i,i=1,m)/)
            do i=1,m
              ncrspd(i) = i
            enddo
         else
            m = n
            nx = ncrspd(1)
            do i = 1, m-1
               ncrspd(i) = ncrspd(i+1)
            end do
            ncrspd(m) = nx
         end if
      end if
    end subroutine set_ncrspd_mxiter

    subroutine mix_pulay_dealloc3
      deallocate(e_wk); deallocate(f_wk); deallocate(ww1); deallocate(finv)
      deallocate(ip)
    end subroutine mix_pulay_dealloc3

    subroutine Resid_and_dd_into_urec_with_hsr(iter)
      integer, intent(in) :: iter
      integer             :: itc,itc0,itc1
      integer :: i,j,k,nmix,imix,ierr
      real(kind=DP) :: sum1,sum2
      itc = ncrspd(iter)
      if(sw_gradient_simplex==ON)then
         do imix=2,iter-1
            itc0 = ncrspd(imix)
            itc1 = ncrspd(imix-1)
            do i=1,nspin_m
               do j=1,kimg
                  do k=ista_kgpm,iend_kgpm
                     urec_l(k,j,i,itc0,iResid) = rho(k,j,i)-rhoo(k,j,i)-(Frhoj(k,j,i,itc1)-rhoj(k,j,i,itc1))
                     urec_l(k,j,i,itc0,iRho  ) = rhoo(k,j,i)-rhoj(k,j,i,itc1)
                  enddo
               enddo
            enddo
            do i=1,nspin_m
               do j=1,nsize_rho_hsr
                  urec_hsr(j,i,itc0,iResid) = rho_hsr(j,i) - rhoo_hsr(j,i) - (Frhoj_hsr(j,i,itc1) - rhoj_hsr(j,i,itc1))
                  urec_hsr(j,i,itc0,iRho  ) = rhoo_hsr(j,i) - rhoj_hsr(j,i,itc1)
               enddo
            enddo
         enddo
         do i=1,nspin_m
            do j=1,kimg
               do k=ista_kgpm,iend_kgpm
                  urec_l(k,j,i,itc,iResid) = rho(k,j,i)-rhoo(k,j,i)-(dout(k,j,i)-din(k,j,i))
                  urec_l(k,j,i,itc,iRho  ) = rhoo(k,j,i)-din(k,j,i)
                  rhoj(k,j,i,itc) = rhoo(k,j,i)
                  Frhoj(k,j,i,itc) = rho(k,j,i)
                  d0_l(k,j,i) = rho(k,j,i) - rhoo(k,j,i)
                  din(k,j,i)  = rhoo(k,j,i)
                  dout(k,j,i) = rho(k,j,i)
               enddo
            enddo
         enddo

         do i=1,nspin_m
            do j=1,nsize_rho_hsr
               urec_hsr(j,i,itc,iResid) = rho_hsr(j,i) - rhoo_hsr(j,i) - (dout_hsr(j,i) - din_hsr(j,i)) ! =dF(=delta F^i)
               urec_hsr(j,i,itc,iRho  ) = rhoo_hsr(j,i) - din_hsr(j,i)                ! =dd
               d0_hsr(j,i) = rho_hsr(j,i) - rhoo_hsr(j,i)
               din_hsr(j,i) = rhoo_hsr(j,i)
               dout_hsr(j,i) = rho_hsr(j,i)
               rhoj_hsr(j,i,itc) = rhoo_hsr(j,i)
               Frhoj_hsr(j,i,itc) = rho_hsr(j,i)
            enddo
         enddo
      else
         do i=1,nspin_m
            do j=1,kimg
               do k=ista_kgpm,iend_kgpm
                  urec_l(k,j,i,itc,iResid) = rho(k,j,i) - rhoo(k,j,i) - (dout(k,j,i) - din(k,j,i)) ! =dF(=delta F^i)
                  urec_l(k,j,i,itc,iRho  ) = rhoo(k,j,i) - din(k,j,i)                ! =dd
                  d0_l(k,j,i) = rho(k,j,i) - rhoo(k,j,i)
                  din(k,j,i)  = rhoo(k,j,i)
                  dout(k,j,i) = rho(k,j,i)
               enddo
            enddo
         enddo
         do i=1,nspin_m
            do j=1,nsize_rho_hsr
               urec_hsr(j,i,itc,iResid) = rho_hsr(j,i) - rhoo_hsr(j,i) - (dout_hsr(j,i) - din_hsr(j,i)) ! =dF(=delta F^i)
               urec_hsr(j,i,itc,iRho  ) = rhoo_hsr(j,i) - din_hsr(j,i)                ! =dd
               d0_hsr(j,i) = rho_hsr(j,i) - rhoo_hsr(j,i)
               din_hsr(j,i) = rhoo_hsr(j,i)
               dout_hsr(j,i) = rho_hsr(j,i)
            enddo
         enddo
      endif
!!$      ynorm(itc,:)=0.d0
      do i=1,nspin_m
         do j=1,kimg
            do k=ista_kgpm,iend_kgpm
!!$               ynorm(itc,i) = ynorm(itc,i)+urec_l(k,j,i,itc,iResid)*urec_l(k,j,i,itc,iResid)
               sum12(i) = sum12(i)+urec_l(k,j,i,itc,iResid)*urec_l(k,j,i,itc,iResid)
            enddo
         enddo
      enddo
!!$      call mpi_allreduce(MPI_IN_PLACE,ynorm(itc,1),nspin_m,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
      call mpi_allreduce(MPI_IN_PLACE,sum12,nspin_m,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
!!$      ynorm(itc,:) = univol*ynorm(itc,:)
      ynorm(itc,:) = univol*sum12(:)
      do i=1,nspin_m
         do j=1,nsize_rho_hsr
            ynorm(itc,i) = ynorm(itc,i)+urec_hsr(j,i,itc,iResid)*urec_hsr(j,i,itc,iResid)
         enddo
      enddo
      ynorm(itc,:) = 1.d0/sqrt(ynorm(itc,:))
      sum1 = 0.d0
      do i=1,nspin_m
         do j=1,kimg
            do k=ista_kgpm,iend_kgpm
               sum1 = sum1 + d0_l(k,j,i)*d0_l(k,j,i)
            enddo
         enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,sum1,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
      sum1 = dsqrt(sum1)
      if(sum1<alpha_pulay_damp_thres) then
         alpha_pulay = alpha_pulay*alpha_pulay_damp
         if(alpha_pulay>0.and.printable) &
         & write(nfout,'(a,e12.5)') '!** new value for the parameter alpha : ',alpha_pulay
      endif
      if(sw_control_stepsize==ON) d0_l_h(:,:,:,itc) = d0_l(:,:,:)
      if(iter>=2.and.sw_control_stepsize==ON)then
         sum2 = 0.d0
         itc1 = ncrspd(iter-1)
         do i=1,nspin_m
            do j=1,kimg
               do k=ista_kgpm,iend_kgpm
                  sum2 = sum2 + d0_l_h(k,j,i,itc1)*d0_l_h(k,j,i,itc1)
               enddo
            enddo
         enddo
         call mpi_allreduce(MPI_IN_PLACE,sum2,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
         sum2 = dsqrt(sum2)
         step_control_factor = max(0.5d0,min(2.0d0,sum2/sum1))
      endif
    end subroutine Resid_and_dd_into_urec_with_hsr

    subroutine Ri_dot_Rj_with_hsr(n)
      integer, intent(in) :: n
      integer  :: it,jt,itc,jtc,is
      real(DP) :: ff1(nspin_m),ff2(nspin_m),ff1tmp  ! real(DP) :: ff1(nspin),ff2(nspin),ff1tmp ! === DEBUG by tkato 2011/11/24 ===

      do it = 1, n
         itc = ncrspd(it)
         do jt = it, n
            jtc = ncrspd(jt)
#ifdef _CDMIX_USE_POINTER_
            urec_l_3 => urec_l(ista_kgpm:iend_kgpm,1:kimg,1:nspin_m,itc,iResid)
            urec_l_3_2 => urec_l(ista_kgpm:iend_kgpm,1:kimg,1:nspin_m,jtc,iResid)
            call mult1s(urec_l_3,urec_l_3_2,f_p,ff1)   ! <delta F^i|delta F^j>
#else
            if(sw_mix_bothspins_sametime==ON .and. nspin_m>1)then
               ff1tmp=0.d0
               call mult1s10_reduce_spin(urec_l,nbxmix,2,itc,iResid,urec_l,jtc,iResid,f_p,ff1tmp)   ! <delta F^i|delta F^j>
               do is=1,nspin_m,af+1
                  ff1tmp = ff1tmp+sum(urec_hsr(:,is,itc,iResid) * urec_hsr(:,is,jtc,iResid))
               enddo
               ff1(1) = ff1tmp
               ff1(2) = ff1tmp
            else
               call mult1s10(urec_l,nbxmix,2,itc,iResid,urec_l,jtc,iResid,f_p,ff1)   ! <delta F^i|delta F^j>
               do is=1,nspin_m,(af+1)
                  ff1(is) = ff1(is)+sum(urec_hsr(:,is,itc,iResid) * urec_hsr(:,is,jtc,iResid))
               enddo
            endif

! ============================= added by K. Tagami ======================= 11.0
            if ( noncol ) then
               ff1tmp=0.d0
               call mult1s10_reduce_spin(urec_l,nbxmix,2,itc,iResid,urec_l,jtc,iResid,f_p,ff1tmp)   ! <delta F^i|delta F^j>
               do is=1,ndim_magmom
                  ff1tmp = ff1tmp+sum(urec_hsr(:,is,itc,iResid) * urec_hsr(:,is,jtc,iResid))
               enddo
               ff1(:) = ff1tmp
            endif
! ======================================================================== 11.0
#endif
            f(it,jt,1:nspin_m) = ff1(1:nspin_m)
            if(jt /= it) f(jt,it,1:nspin_m) = f(it,jt,1:nspin_m)
         end do
      end do
    end subroutine Ri_dot_Rj_with_hsr

    subroutine Rj_dot_d_with_hsr(n)
      integer, intent(in) :: n
      integer  :: jt, jtc, is
      real(DP) :: ff1(nspin_m),ff2(nspin_m) !  real(DP) :: ff1(nspin),ff2(nspin) ! === DEBUG by tkato 2011/11/24 ===
      real(DP) :: ff1tmp
      do jt = 1, n
         jtc = ncrspd(jt)
#ifdef _CDMIX_USE_POINTER_
         urec_l_3 => urec_l(ista_kgpm:iend_kgpm,1:kimg,1:nspin_m,jtc,iResid)
         call mult1s(urec_l_3,d0_l,f_p,ff1)
#else
         if(sw_mix_bothspins_sametime==ON .and. nspin_m>1)then
            ff1tmp=0.d0
            call mult1s5_reduce_spin(urec_l,nbxmix,2,jtc,iResid,d0_l,f_p,ff1tmp)
            do is=1,nspin_m,af+1
               ff1tmp = ff1tmp+sum(urec_hsr(:,is,jtc,iResid) * d0_hsr(:,is))
            enddo
            ff1(1) = ff1tmp;ff1(2) = ff1tmp
         else
            call mult1s5(urec_l,nbxmix,2,jtc,iResid,d0_l,f_p,ff1)
            ff2=0.d0
            do is=1,nspin_m,af+1
               ff2(is) = ff2(is)+sum(urec_hsr(:,is,jtc,iResid) * d0_hsr(:,is))
            enddo
            ff1(:) = ff1(:)+ff2(:)
         endif

! =========================== added by K. Tagami ================== 11.0
         if ( noncol ) then
            ff1tmp=0.d0
            call mult1s5_reduce_spin(urec_l,nbxmix,2,jtc,iResid,d0_l,f_p,ff1tmp)
            do is=1,ndim_magmom
               ff1tmp = ff1tmp+sum(urec_hsr(:,is,jtc,iResid) * d0_hsr(:,is))
            enddo
            ff1(:) = ff1tmp
         endif
! ================================================================ 11.0

#endif
         uuf_p(jt,1:nspin_m) = ff1(1:nspin_m)
      end do
    end subroutine Rj_dot_d_with_hsr

    subroutine get_finv_lapack_with_hsr(m,n,f)
      integer,intent(in)                             :: m,n
      real(DP),intent(inout),dimension(m,m,nspin_m) :: f
      real(DP), allocatable,dimension(:,:) :: fwork
      integer :: is,inf,it,jt,kt,nnspin
      real(DP) :: div,tmp
      allocate(fwork(n,n))
      nnspin = nspin
      if(sw_mix_bothspins_sametime==ON .or. af==1) nnspin=1

! ============================== added by K. Tagami ============== 11.0
      if ( noncol )  nnspin = 1
! ================================================================ 11.0

      do is=1,nnspin
         if(ipripulay >= 2) then
            write(nfout,600) n,(('(',it,jt,')',f(it,jt,is),jt=1,n),it=1,n)
600         format(//11x,"**input matrix**"/12x &
                 & ,"horder=",I5/(2x,4(1x,1a,i2,",",i2,1a,e14.6)))
         end if
         fwork=0
         do it=1,n
            do jt=1,n
               fwork(jt,it) = f(jt,it,is)*ynorm(jt,is)*ynorm(it,is)
               if(it==jt) fwork(jt,it)=fwork(jt,it)+alpha_pulay
            enddo
         enddo
         call dpotrf('U',n,fwork,n,inf)
         call dpotri('U',n,fwork,n,inf)
         do it=1,n-1
            do jt=it+1,n
               fwork(jt,it) = fwork(it,jt)
            enddo
         enddo
         do it=1,n
            do jt=1,n
               f(jt,it,is) = fwork(jt,it)*ynorm(jt,is)*ynorm(it,is)
            enddo
         enddo
         if(ipripulay >= 2) then
            write(nfout,630) (('(',it,jt,')',f(it,jt,is),it=1,n),jt=1,n)
630         format(/11x, "**inverse matrix**" &
                 & ,/(2x,4(1x,1a,i2,",",i2,1a,e14.6)))
         end if
      enddo
      if(sw_mix_bothspins_sametime==ON .and. nspin_m>1)then
         do it=1,n
            do jt=1,n
               f(jt,it,2) = f(jt,it,1)
            enddo
         enddo
      endif

! ================================= added by K. Tagami ================= 11.0
      if ( noncol ) then
         do it=1,n
            do jt=1,n
               f(jt,it,:) = f(jt,it,1)
            enddo
         end do
      endif
! ====================================================================== 11.0
      deallocate(fwork)
    end subroutine get_finv_lapack_with_hsr

    subroutine get_gmatrix_with_hsr(n)
      integer,intent(in) :: n
      integer :: is, it, jt, nnspin
      nnspin = nspin
      if(sw_mix_bothspins_sametime==ON .or. af==1) nnspin=1

! ============================== added by K. Tagami ============== 11.0
      if ( noncol )  nnspin = 1
! ================================================================ 11.0

      g_p = 0.d0
      do is = 1, nnspin
         do it = 1, n
            do jt = 1, n
               g_p(it,is) = g_p(it,is) - f(jt,it,is)*uuf_p(jt,is)
            end do
         end do
         if(ipripulay >= 2) then
            write(nfout,'(" -- g_p(1:",i3,") --")') n
            write(nfout,'(8f20.12)') (g_p(it,is),it=1,n)
         end if
      end do
      if(sw_mix_bothspins_sametime==ON .and. nspin_m>1)then
         do it = 1,n
            g_p(it,2) = g_p(it,1)
         enddo
      endif
! ================================= added by K. Tagami ================= 11.0
      if ( noncol ) then
         do it = 1,n
            g_p(it,:) = g_p(it,1)
         end do
      end if
! ====================================================================== 11.0
    end subroutine get_gmatrix_with_hsr

    subroutine renew_d_using_g_with_hsr(n,p)
      integer, intent(in)                                :: n
      real(DP),intent(in),dimension(ista_kgpm:iend_kgpm,nspin_m) :: p
      integer    :: is, k, i, it, itc, ns

!!$      do is = 1, nspin, af+1
      ns = nspin_for_qnewton()
      do is = 1, ns,af+1
         do k = 1, kimg
            do i = ista_kngp, iend_kngp
               rho(i,k,is)  = rhoo(i,k,is) + p(i,is)*d0_l(i,k,is)
            end do
            do it = 1, n
               itc = ncrspd(it)
               do i = ista_kngp, iend_kngp
                  rho(i,k,is) = rho(i,k,is) + g_p(it,is)* &
                       &        (urec_l(i,k,is,itc,iRho) + p(i,is)*urec_l(i,k,is,itc,iResid))
               end do
            end do
         end do
      end do
      do is=1,ns,af+1
         rho_hsr(:,is) = rhoo_hsr(:,is) + ommix_factor*rmxtrc(is) * d0_hsr(:,is)
         do it = 1, n
            itc = ncrspd(it)
            rho_hsr(:,is) = rho_hsr(:,is) + g_p(it,is) * &
          & (urec_hsr(:,is,itc,iRho) + ommix_factor*rmxtrc(is)*urec_hsr(:,is,itc,iResid))
         enddo
      enddo

! ============================== modified by K. Tagami ================ 11.0
!      if(sw_force_simple_mixing==ON .and. sw_recomposing==ON) call simple_mix2(c_p)
!      if(sw_force_simple_mixing_hsr==ON .and. sw_recomposing_hsr==ON ) call simple_mix2_kt( rmxtrc )
!
      if ( .not. noncol ) then
         if(sw_force_simple_mixing==ON .and. sw_recomposing==ON) then
            call simple_mix2(c_p)
         endif
         if(sw_force_simple_mixing_hsr==ON .and. sw_recomposing_hsr==ON ) then
            call simple_mix2_kt( rmxtrc )
         endif
      endif
! ===================================================================== 11.0

      if(kgpm < kgp) then
         call concentrate_d_to_chg(rho,chgq_l) !-(m_C.D.)
         call simple_mix_large_Gc(c_p,chgqo_l,chgq_l)     !-(m_C.D.) chgq,chgqo,c_p ->chgq
      end if
    end subroutine renew_d_using_g_with_hsr
  end subroutine m_CD_mix_pulay_with_hsr
!!$ 11.07 AS Pulay version of 'sw_mix_charge_hardpart'

! ========================= added by K. Tagami ===================== 5.0
  subroutine create_map_func(paramset,mixocc)
    logical :: paramset,mixocc
    integer :: n, ia, it
    integer :: lmt1, lmt2
    integer :: ig,m2,m1,i,ip

    n=0
    do ia=1,natm
       it = ityp(ia)

       do lmt1=1, ilmt(it)
          do lmt2 = lmt1, ilmt(it)
             n=n+1

             if(.not.paramset) then
                imap_hsr(n) = ia + natm *(lmt1-1) + natm*nlmt*( lmt2 -1 )
                diag_elem(n) = lmt1.eq.lmt2
             endif
          end do
       end do
    end do
    if ( noncol ) then
       nsize_rho_hsr_realpart = n

       if ( sw_mix_imaginary_hardpart == ON ) then
          do ia=1,natm
             it = ityp(ia)
             do lmt1=1, ilmt(it)
                do lmt2 = lmt1 +1, ilmt(it)
                   n=n+1
                   if(.not.paramset) &
                        & imap_hsr(n) = ia + natm *(lmt1-1) + natm*nlmt*( lmt2 -1 )
                end do
             end do
          end do
          nsize_rho_hsr = n
       endif
    endif

! ----
    nsize_rho_hsr0 = n
    nsize_rho_hsr = nsize_rho_hsr0

    if (mixocc) then
       n=0
       do ia=1,natm
          ig = iproj_group(ia)
          if(ig<1) cycle
          do i=1,num_proj_elems(ig)
             ip=proj_group(i,ig)
             it = proj_attribute(ip)%ityp
             do m2=1,i2lp(ip)
                do m1=m2,i2lp(ip)
                   n=n+1
                   if(.not.paramset) &
                        & imap_om(n) = m1 +max2lp *( m2 -1 &
                        &                 +max2lp*( i -1 +max_projs*( ia-1 ) ) )
                end do
             end do
          end do
       enddo
       if ( noncol .and. sw_mix_imaginary_hardpart == ON ) then
          nsize_rho_om_realpart = n

          do ia=1,natm
             ig = iproj_group(ia)
             if(ig<1) cycle
             do i=1,num_proj_elems(ig)
                ip=proj_group(i,ig)
                it = proj_attribute(ip)%ityp
                do m2=1,i2lp(ip)
                   do m1=m2,i2lp(ip)
                      n=n+1
                      if(.not.paramset) &
                           & imap_om(n) = m1 +max2lp *( m2 -1 &
                           &                 +max2lp*( i -1 +max_projs*( ia-1 ) ) )
                   end do
                end do
             end do
          end do
       endif
       nsize_rho_om = n
       nsize_rho_hsr = nsize_rho_hsr0 +nsize_rho_om
    endif

  end subroutine create_map_func

  subroutine alloc_rho_hsr(mixocc)
    logical, intent(in) :: mixocc
    if ( noncol ) then
       nspin_m  = ndim_magmom
    else
       nspin_m = nspin/(af+1)
    endif
    if ( noncol ) then
       nspin_m  = ndim_magmom
       allocate( rho_hsr( nsize_rho_hsr,ndim_magmom)); rho_hsr = 0.0d0
       allocate( rhoo_hsr(nsize_rho_hsr,ndim_magmom));rhoo_hsr = 0.0d0
    else
       nspin_m  = nspin/(af+1)
       allocate( rho_hsr( nsize_rho_hsr,nspin_m)); rho_hsr = 0.0d0
       allocate( rhoo_hsr(nsize_rho_hsr,nspin_m));rhoo_hsr = 0.0d0
    endif
! ================================ modified by K. Tagami ============== 11.0
!    allocate( rho_hsr( nsize_rho_hsr,nspin)); rho_hsr = 0.0d0
!    allocate( rhoo_hsr(nsize_rho_hsr,nspin));rhoo_hsr = 0.0d0
! ====================================================================== 11.0
    allocate( imap_hsr(nsize_rho_hsr0) ); imap_hsr = 0
    allocate(diag_elem(nsize_rho_hsr0));diag_elem=.false.
    if(mixocc)then
        allocate( imap_om(nsize_rho_om) ); imap_om = 0
    endif
  end subroutine alloc_rho_hsr

  subroutine set_i2lp_max2lp()
    integer :: it,ip
    integer, parameter :: ntau0=2
    integer :: nsize ! === added by K. Tagami === 11.0

    allocate(i2lp(num_projectors))
    do ip=1,num_projectors
       i2lp(ip) = 2*proj_attribute(ip)%l+1
    end do
    max2lp = 0
    do ip=1,num_projectors
       if(i2lp(ip) > max2lp) then
          max2lp = i2lp(ip)
          l1max  = proj_attribute(ip)%l+1
       end if
    end do

! =========================== modified by K. Tagami ====================== 11.0
!!
!!    nyymax = ntau0*l1max**2*(l1max**2+1)/2
!
    nsize = ntau0*( 2*( l1max -1 )+1 )
    nyymax = nsize *( nsize +1 ) /2

  end subroutine set_i2lp_max2lp

  subroutine dealloc_rho_hsr
    deallocate(rho_hsr)
    deallocate(rhoo_hsr)
    deallocate(imap_hsr)
    deallocate(diag_elem)
    if(allocated(imap_om)) deallocate(imap_om)
    if(allocated(i2lp))    deallocate(i2lp)
  end subroutine dealloc_rho_hsr

  subroutine map_hsr_to_rho( hsr,rho )
    real(kind=DP), intent(in) :: hsr( natm *nlmt *nlmt, nspin )
    real(kind=DP), intent(out) :: rho( nsize_rho_hsr,nspin)

    integer :: i,is

    do is=1,nspin,(af+1)
       rho(:,is)=0.d0
       do i=1,nsize_rho_hsr0
          rho(i,is) = hsr( imap_hsr(i),is )
       end do
    end do
  end subroutine map_hsr_to_rho

! ============================================ added by K. Tagami ======== 11.0
  subroutine map_hsr_to_rho_noncl( hsr, hsi, rho )
    real(kind=DP), intent(in) :: hsr( natm *nlmt *nlmt, ndim_magmom )
    real(kind=DP), intent(in) :: hsi( natm *nlmt *nlmt, ndim_magmom )
    real(kind=DP), intent(out) :: rho( nsize_rho_hsr,ndim_magmom )

    integer :: i,is

    do is=1,ndim_magmom,(af+1)

       do i=1,nsize_rho_hsr_realpart
          rho(i,is) = hsr( imap_hsr(i),is )
       end do
       if ( sw_mix_imaginary_hardpart == ON ) then
          do i=nsize_rho_hsr_realpart+1, nsize_rho_hsr0
             rho(i,is) = hsi( imap_hsr(i),is )
          end do
       endif
    end do
  end subroutine map_hsr_to_rho_noncl
! ===================================================================== 11.0

  subroutine map_om_to_rho(om,rho)
    real(kind=DP), intent(in) :: om(max2lp*max2lp*max_projs*natm,nspin)
    real(kind=DP), intent(out) :: rho(nsize_rho_hsr,nspin)

    integer :: i,is

    do is=1,nspin,(af+1)
       do i=nsize_rho_hsr0+1,nsize_rho_hsr
          rho(i,is) = om(imap_om(i-nsize_rho_hsr0),is)
       end do
    end do
  end subroutine map_om_to_rho

  subroutine map_om_to_rho_noncl(om_re,om_im,rho)
    real(kind=DP), intent(in) :: om_re(max2lp*max2lp*max_projs*natm,ndim_magmom)
    real(kind=DP), intent(in) :: om_im(max2lp*max2lp*max_projs*natm,ndim_magmom)
    real(kind=DP), intent(out) :: rho(nsize_rho_hsr,ndim_magmom)

    integer :: i,is

    do is=1, ndim_magmom
       do i=nsize_rho_hsr0+1, nsize_rho_hsr0 +nsize_rho_om_realpart
          rho(i,is) = om_re(imap_om(i-nsize_rho_hsr0),is)
       end do
       if ( sw_mix_imaginary_hardpart == ON ) then
          do i=nsize_rho_hsr0 +nsize_rho_om_realpart+1, nsize_rho_hsr
             rho(i,is) = om_im(imap_om(i-nsize_rho_hsr0),is)
          end do
       endif
    end do
  end subroutine map_om_to_rho_noncl

  subroutine map_rho_to_om(om,rho)
    real(kind=DP), intent(out) :: om(max2lp*max2lp*max_projs*natm,nspin)
    real(kind=DP), intent(in) :: rho(nsize_rho_hsr,nspin)
    integer :: i,is,ia,ig,ip,it,m1,m2

    do is=1,nspin,(af+1)
       do i=nsize_rho_hsr0+1,nsize_rho_hsr
          om(imap_om(i-nsize_rho_hsr0),is) = rho(i,is)
       end do
    end do
    call symmetrize(om)

  contains

    subroutine symmetrize(om)
      real(kind=DP), intent(inout) :: om(max2lp,max2lp,max_projs,natm,nspin)

      do is=1,nspin,(af+1)
         do ia=1,natm
            ig = iproj_group(ia)
            if(ig<1) cycle
            do i=1,num_proj_elems(ig)
               ip=proj_group(i,ig)
               it = proj_attribute(ip)%ityp
               do m2=1,i2lp(ip)
                  do m1=m2,i2lp(ip)
                     if(m1/=m2) om(m2,m1,i,ia,is) = om(m1,m2,i,ia,is)
                  end do
               end do
            end do
         end do
      end do
    end subroutine symmetrize

  end subroutine map_rho_to_om

  subroutine map_rho_to_om_noncl(om_re,om_im,rho)
    real(kind=DP), intent(out) :: om_re(max2lp*max2lp*max_projs*natm,ndim_magmom)
    real(kind=DP), intent(out) :: om_im(max2lp*max2lp*max_projs*natm,ndim_magmom)
    real(kind=DP), intent(in) :: rho(nsize_rho_hsr,ndim_magmom)
    integer :: i,is,ia,ig,ip,it,m1,m2

    do is=1,ndim_magmom
       do i=nsize_rho_hsr0+1, nsize_rho_hsr0 +nsize_rho_om_realpart
          om_re(imap_om(i-nsize_rho_hsr0),is) = rho(i,is)
       end do
       if ( sw_mix_imaginary_hardpart == ON ) then
          do i=nsize_rho_hsr0 +nsize_rho_om_realpart+1, nsize_rho_hsr
             om_im(imap_om(i-nsize_rho_hsr0),is) = rho(i,is)
          end do
       endif
    end do
    call symmetrize(om_re,om_im)

  contains

    subroutine symmetrize(om_re, om_im)
      real(kind=DP), intent(inout) :: om_re(max2lp,max2lp,max_projs,natm,ndim_magmom)
      real(kind=DP), intent(inout) :: om_im(max2lp,max2lp,max_projs,natm,ndim_magmom)

      do is=1, ndim_magmom
         do ia=1,natm
            ig = iproj_group(ia)
            if(ig<1) cycle
            do i=1,num_proj_elems(ig)
               ip=proj_group(i,ig)
               it = proj_attribute(ip)%ityp
               do m2=1,i2lp(ip)
                  do m1=m2,i2lp(ip)
                     if (m1/=m2) then
                        om_re(m2,m1,i,ia,is) =  om_re(m1,m2,i,ia,is)
                        om_im(m2,m1,i,ia,is) = -om_im(m1,m2,i,ia,is)
                     endif
                  end do
               end do
            end do
         end do
      end do
    end subroutine symmetrize

  end subroutine map_rho_to_om_noncl

  subroutine map_rho_to_hsr( hsr,rho )
    real(kind=DP), intent(out) :: hsr( natm *nlmt *nlmt, nspin )
    real(kind=DP), intent(in) :: rho( nsize_rho_hsr,nspin )

    integer :: i, is, ia, lmt1, lmt2


    do is=1,nspin,(af+1)
       hsr(1:nsize_rho_hsr0,is) = 0.0d0
       do i=1,nsize_rho_hsr0
          hsr(imap_hsr(i),is) = rho(i,is)
       end do
    end do
!    call symmetrize(hsr)
    call symmtrz_of_ff()
  contains
    subroutine symmetrize(hsr)
      real(kind=DP), intent(inout) :: hsr( natm, nlmt, nlmt, nspin )
      integer :: ia,ja
      do is=1,nspin,(af+1)
         do ia=1,natm
            do lmt1=1, nlmt
               do lmt2=lmt1, nlmt
                  if ( lmt1/=lmt2 ) hsr(ia,lmt2,lmt1,is) = hsr(ia,lmt1,lmt2,is)
               end do
            end do
         end do
      end do
      if(af/=0)then
         do ia = 1, natm
            ja=abs(ia2ia_symmtry_op_inv(ia,nopr+af))
            if(ja <= 0 .or. natm < ia) cycle
            hsr(ia,:,:,nspin) = hsr(ja,:,:,1)
         end do
      endif
    end subroutine symmetrize
  end subroutine map_rho_to_hsr
! ========================================================================= 5.0

  subroutine map_rho_to_hsr_noncl( hsr,hsi,rho )
    real(kind=DP), intent(out) :: hsr( natm *nlmt *nlmt, ndim_magmom )
    real(kind=DP), intent(out) :: hsi( natm *nlmt *nlmt, ndim_magmom )
    real(kind=DP), intent(in) :: rho( nsize_rho_hsr,ndim_magmom )

    integer :: i, is, ia, lmt1, lmt2

    hsr = 0.0d0
    if ( sw_mix_imaginary_hardpart == ON ) hsi = 0.0d0

    do is=1,ndim_magmom,(af+1)
       do i=1,nsize_rho_hsr_realpart
          hsr(imap_hsr(i),is) = rho(i,is)
       end do
       if ( sw_mix_imaginary_hardpart == ON ) then
          do i=nsize_rho_hsr_realpart+1, nsize_rho_hsr0
             hsi(imap_hsr(i),is) = rho(i,is)
          end do
       endif
    end do

!    call symmetrize(hsr,hsi)
    call m_CD_symmtrz_of_ff_noncl_C( hsr, hsi )

  contains
    subroutine symmetrize(hsr,hsi)
      real(kind=DP), intent(inout) :: hsr( natm, nlmt, nlmt, ndim_magmom )
      real(kind=DP), intent(inout) :: hsi( natm, nlmt, nlmt, ndim_magmom )
      integer :: ia,ja
      do is=1,ndim_magmom,(af+1)
         do ia=1,natm
            do lmt1=1, nlmt
               do lmt2=lmt1, nlmt
                  if ( lmt1/=lmt2 ) then
                     hsr(ia,lmt2,lmt1,is) = hsr(ia,lmt1,lmt2,is)
                     if ( sw_mix_imaginary_hardpart == ON ) then
                        hsi(ia,lmt2,lmt1,is) = -hsi(ia,lmt1,lmt2,is)
                     endif
                  endif
               end do
            end do
         end do
      end do
      if(af/=0)then
         do ia = 1, natm
            ja=abs(ia2ia_symmtry_op_inv(ia,nopr+af))
            if(ja <= 0 .or. natm < ia) cycle
            hsr(ia,:,:,nspin) = hsr(ja,:,:,1)
         end do
      endif
    end subroutine symmetrize

  end subroutine map_rho_to_hsr_noncl
! ========================================================================= 11.0

! ========================= added by K. Tagami========================== 5.0
  subroutine m_CD_simple_mixing_hard( nfout,rmxt )
    integer, intent(in) :: nfout
    real(kind=DP), intent(in) :: rmxt

    integer :: is

    if ( noncol ) then   !    nspin_m  = nspin/(af+1)  ! === modified by K.Tagami === 11.0
       nspin_m = ndim_magmom
    else
       nspin_m  = nspin/(af+1)   
    endif

    allocate( rmxtrc(nspin_m) )

    if ( noncol ) then  ! === modified by K. Tagami === 11.0
       rmxtrc = rmxt
       rmxtrc(2:nspin_m) = min( rmxt *spin_density_mixfactor, rmx_max )
    else
       if( sw_recomposing_hsr == YES .and. af == 0 .and. nspin == 2) then
          call alloc_hsrstore_recomp(rmxt,rmxtrc)
       else
          rmxtrc = rmxt
       endif
    end if

    Do is=1, ndim_magmom, af+1  !    Do is=1, nspin, af+1 ! === modified by K. Tagami === 11.0
       hsr(:,:,:,is) = rmxtrc(is) *hsr(:,:,:,is) &
            &             + ( 1.d0-rmxtrc(is) )*hsro(:,:,:,is)
    End do

! ================================= added by K. Tagami ============= 11.0
    if ( noncol .and. sw_mix_imaginary_hardpart == ON ) then
       Do is=1, ndim_magmom
          hsi(:,:,:,is) = rmxtrc(is) *hsi(:,:,:,is) &
               &             + ( 1.d0-rmxtrc(is) )*hsio(:,:,:,is)
       End do
    endif
! ================================================================== 11.0

    if ( .not. noncol ) then  ! === modified by K. Tagami === 11.0
       if ( sw_recomposing_hsr == YES .and. af == 0 .and. nspin == 2) then
          call compose_hsr_dealloc_store
       end if
    endif
    deallocate(rmxtrc)

  end subroutine m_CD_simple_mixing_hard
! ============================================================== 5.0

  subroutine m_CD_simple_mixing_hsr(nfout,rmxt)
    integer ,intent(in)      :: nfout
    real(kind=DP),intent(in) :: rmxt

    integer:: ia,it,is,lmt1,lmt2
                                                  __TIMER_SUB_START(1144)
                                                  __TIMER_DO_START(1192)
    do ia = 1, natm
       it = ityp(ia)
       if(ipaw(it)/=1) cycle

       if ( noncol ) then   ! ==================================== modified by K. Tagami ============ 11.0
          do is = 1, ndim_magmom  !    do is = 1, nspin, af+1
             do lmt1 = 1, ilmt(it)
                do lmt2 = lmt1, ilmt(it)
                   hsr(ia,lmt1,lmt2,is) = rmxt*hsr(ia,lmt1,lmt2,is) + &
                        & (1-rmxt)*hsro(ia,lmt1,lmt2,is)
                end do! lmt2
             end do! lmt1
          end do! is
! --
          if ( sw_mix_imaginary_hardpart == ON ) then
             do is = 1, ndim_magmom !    do is = 1, nspin, af+1
                do lmt1 = 1, ilmt(it)
                   do lmt2 = lmt1, ilmt(it)
                      hsi(ia,lmt1,lmt2,is) = rmxt*hsi(ia,lmt1,lmt2,is) &
                           &               +(1-rmxt)*hsio(ia,lmt1,lmt2,is)
                   end do! lmt2
                end do! lmt1
             end do
          end if
! -----
       else
          do is = 1, nspin, af+1
             do lmt1 = 1, ilmt(it)
                do lmt2 = lmt1, ilmt(it)
                   hsr(ia,lmt1,lmt2,is) = rmxt*hsr(ia,lmt1,lmt2,is) &
                        &               +(1-rmxt)*hsro(ia,lmt1,lmt2,is)
                 end do! lmt2
              end do! lmt1
         end do! is
      endif                 ! ========================================================================== 11.0

    end do! ia
                                                  __TIMER_DO_STOP(1192)
                                                  __TIMER_SUB_STOP(1144)
  end subroutine m_CD_simple_mixing_hsr

! == KT_add ==== 2014/09/19
  subroutine mix_broyden_allocate_intg
    if ( noncol ) then
       nspin_m  = ndim_magmom
    else
       nspin_m  = nspin/(af+1)
    endif

    allocate(f_p(ista_kgpm:iend_kgpm)); f_p = 0.0d0
    call precon_4_mult(f_p) !-(m_CD)

    if(hownew == RENEW) then
       if ( noncol ) then
          allocate(f(nbxmix,nbxmix,ndim_magmom))
       else
          allocate(f(nbxmix,nbxmix,nspin))
       endif
       allocate(g(nbxmix));
       f = 0.0d0;   g = 0.0d0
    end if

    allocate(ncrspd(nbxmix)); ncrspd = 0

! ---
    allocate(din(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(dout(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(dF_l(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(urec_l(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix,2))
    allocate(prj_wk(mp_kngp*kimg*nspin_m))
    din = 0.0d0; dout = 0.0d0; dF_l = 0.0d0; urec_l = 0.0d0; prj_wk = 0.0d0

    if ( sw_mix_charge_hardpart == ON ) then
       allocate( din_hsr(nsize_rho_hsr,nspin_m) )
       allocate( dout_hsr(nsize_rho_hsr,nspin_m))
       allocate( dF_hsr(nsize_rho_hsr,nspin_m))
       allocate( urec_hsr(nsize_rho_hsr,nspin_m,nbxmix,2) )

       din_hsr = 0.0d0; dout_hsr = 0.0d0; dF_hsr = 0.0d0; urec_hsr = 0.0d0
    endif

    if ( sw_mix_charge_with_ekindens == ON ) then
       allocate(din_ekinq(ista_kgpm:iend_kgpm,kimg,nspin_m))
       allocate(dout_ekinq(ista_kgpm:iend_kgpm,kimg,nspin_m))
       allocate(dF_l_ekinq(ista_kgpm:iend_kgpm,kimg,nspin_m))
       allocate(urec_l_ekinq(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix,2))
    endif

    if ( sw_mix_charge_with_ekindens == ON ) then
       allocate(f_p_ekinq(ista_kgpm:iend_kgpm)); f_p_ekinq = 1.0d0
    endif

  end subroutine mix_broyden_allocate_intg

  subroutine mix_broyden_deallocate_intg
    if(allocated(f_p)) deallocate(f_p)
    if(allocated(din)) deallocate(din)
    if(allocated(dout)) deallocate(dout)
    if(allocated(dF_l)) deallocate(dF_l)
    if(allocated(urec_l)) deallocate(urec_l)
    if(allocated(prj_wk)) deallocate(prj_wk)
    if(allocated(f)) deallocate(f)
    if(allocated(g)) deallocate(g)
    if(allocated(ncrspd)) deallocate(ncrspd)
!
    if ( sw_mix_charge_hardpart == ON ) then
       if(allocated(din_hsr)) deallocate(din_hsr)
       if(allocated(dout_hsr)) deallocate(dout_hsr)
       if(allocated(dF_hsr)) deallocate(dF_hsr)
       if(allocated(urec_hsr)) deallocate(urec_hsr)
    endif

    if ( sw_mix_charge_with_ekindens == ON ) then
       if(allocated(din_ekinq)) deallocate(din_ekinq)
       if(allocated(dout_ekinq)) deallocate(dout_ekinq)
       if(allocated(dF_l_ekinq)) deallocate(dF_l_ekinq)
       if(allocated(urec_l_ekinq)) deallocate(urec_l_ekinq)
       if(allocated(f_p_ekinq)) deallocate(f_p_ekinq)
    endif

  end subroutine mix_broyden_deallocate_intg

  subroutine mix_PULAY_deallocate_intg
    if(allocated(f_p)) deallocate(f_p)
    if(allocated(din)) deallocate(din)
    if(allocated(dout)) deallocate(dout)
    if(allocated(urec_l)) deallocate(urec_l)
    if(allocated(uuf_p)) deallocate(uuf_p)
    if(allocated(f)) deallocate(f)
    if(allocated(g_p)) deallocate(g_p)
    if(allocated(prj_wk)) deallocate(prj_wk)
    if(allocated(ncrspd)) deallocate(ncrspd)
    if(allocated(rhoj)) deallocate(rhoj)
    if(allocated(Frhoj)) deallocate(Frhoj)
    if(allocated(rhojo)) deallocate(rhojo)
    if(allocated(Frhojo)) deallocate(Frhojo)
    if(allocated(d0_l_h)) deallocate(d0_l_h)
    if (allocated(ynorm)) deallocate(ynorm)

    if( sw_mix_charge_hardpart == ON ) then
       if ( allocated( din_hsr)) deallocate( din_hsr)
       if ( allocated(dout_hsr)) deallocate(dout_hsr)
       if ( allocated(urec_hsr)) deallocate(urec_hsr)
       if ( allocated(rhoj_hsr)) deallocate(rhoj_hsr)
       if ( allocated(Frhoj_hsr)) deallocate(Frhoj_hsr)
       if ( allocated(rhoj_hsro)) deallocate(rhoj_hsro)
       if ( allocated(Frhoj_hsro)) deallocate(Frhoj_hsro)
    endif

    if ( sw_mix_charge_with_ekindens == ON ) then
       if(allocated(din_ekinq)) deallocate(din_ekinq)
       if(allocated(dout_ekinq)) deallocate(dout_ekinq)
       if(allocated(urec_l_ekinq)) deallocate(urec_l_ekinq)

       if(allocated(rhoj_ekinq)) deallocate(rhoj_ekinq)
       if(allocated(Frhoj_ekinq)) deallocate(Frhoj_ekinq)
       if(allocated(rhojo_ekinq)) deallocate(rhojo_ekinq)
       if(allocated(Frhojo_ekinq)) deallocate(Frhojo_ekinq)
       if(allocated(d0_l_h_ekinq)) deallocate(d0_l_h_ekinq)

       if(allocated(f_p_ekinq)) deallocate(f_p_ekinq)
    endif

  end subroutine mix_PULAY_deallocate_intg

  subroutine mix_dealloc_previous_intg()
    if(previous_waymix == BROYD2) then
       call mix_broyden_deallocate_intg()
    else if(previous_waymix == PULAY) then
       call mix_PULAY_deallocate_intg()
    end if
  end subroutine mix_dealloc_previous_intg

  subroutine mix_broyden_alloc2_intg
    allocate( d0_l(ista_kgpm:iend_kgpm,kimg,nspin_m) ); d0_l = 0.0d0
    allocate(  u_l(ista_kgpm:iend_kgpm,kimg,nspin_m) );  u_l = 0.0d0
    allocate(  v_l(ista_kgpm:iend_kgpm,kimg,nspin_m) );  v_l = 0.0d0

    if ( sw_mix_charge_hardpart == ON ) then
       allocate( d0_hsr( nsize_rho_hsr,nspin_m ) ); d0_hsr = 0.0d0
       allocate(  u_hsr( nsize_rho_hsr,nspin_m ) );  u_hsr = 0.0d0
       allocate(  v_hsr( nsize_rho_hsr,nspin_m ) );  v_hsr = 0.0d0
    endif
    if ( sw_mix_charge_with_ekindens == ON ) then
       allocate( d0_l_ekinq(ista_kgpm:iend_kgpm,kimg,nspin_m ) ); d0_l_ekinq = 0.0d0
       allocate(  u_l_ekinq(ista_kgpm:iend_kgpm,kimg,nspin_m ) );  u_l_ekinq = 0.0d0
       allocate(  v_l_ekinq(ista_kgpm:iend_kgpm,kimg,nspin_m ) );  v_l_ekinq = 0.0d0
    endif

    call alloc_rho_rhoo_and_cpm_intg
  end subroutine mix_broyden_alloc2_intg

  subroutine alloc_rho_rhoo_and_cpm_intg
    if(kgpm < kgp .and. npes /= 1 ) then
       allocate(  rho(ista_kgpm:iend_kgpm,kimg,nspin_m ) )
       allocate( rhoo(ista_kgpm:iend_kgpm,kimg,nspin_m ) )
       allocate( c_pm(ista_kgpm:iend_kgpm,nspin_m) )
       rho = 0.0d0; rhoo = 0.0d0 ; c_pm = 0.0d0

       call scatter_chg_onto_d(chgq_l,rho)
       call scatter_chg_onto_d(chgqo_l,rhoo)
       call scatter_cp_onto_cpm(c_p,c_pm)
    else
       rho => chgq_l; rhoo => chgqo_l; c_pm => c_p
    end if

    if ( sw_mix_charge_with_ekindens == ON )  then
       if(kgpm < kgp .and. npes /= 1 ) then
          allocate(  rho_ekinq(ista_kgpm:iend_kgpm,kimg,nspin_m ) )
          allocate( rhoo_ekinq(ista_kgpm:iend_kgpm,kimg,nspin_m ) )
          allocate( c_pm_ekinq(ista_kgpm:iend_kgpm,nspin_m) )
          rho_ekinq = 0.0d0; rhoo_ekinq = 0.0d0; c_pm_ekinq = 0.0d0

          call scatter_chg_onto_d(ekinq_l, rho_ekinq)
          call scatter_chg_onto_d(ekinqo_l,rhoo_ekinq)
          call scatter_chg_onto_d(c_p_ekinq, c_pm_ekinq)
       else
          rho_ekinq => ekinq_l; rhoo_ekinq => ekinqo_l;  c_pm_ekinq => c_p_ekinq
       end if
    endif
  end subroutine alloc_rho_rhoo_and_cpm_intg

  subroutine simple_mix1_intg( p, rmx_this, p2 )
    real(kind=DP), intent(in) :: p(ista_kngp:iend_kngp,nspin_m)
    real(kind=DP), intent(in) :: rmx_this(nspin_m)
    real(kind=DP), intent(in), optional :: p2(ista_kngp:iend_kngp,nspin_m)

    integer  :: is,k

    if(kgpm == kgp .or. npes == 1) then
       do is = 1, ndim_magmom, af+1
          din (ista_kgpm:iend_kgpm,:,is) = chgqo_l(ista_kgpm:iend_kgpm,:,is)
          dout(ista_kgpm:iend_kgpm,:,is) = chgq_l (ista_kgpm:iend_kgpm,:,is)
       end do
    else
       call scatter_chg_onto_d(chgqo_l,din )  ! -(m_C.D.)
       call scatter_chg_onto_d(chgq_l, dout)  ! -(m_C.D.)
    end if
    do is = 1, ndim_magmom, af+1
       do k = 1, kimg
          chgq_l(:,k,is) = p(:,is)*chgq_l(:,k,is) + (1.0d0-p(:,is))*chgqo_l(:,k,is)
       end do
    end do

    if ( sw_mix_charge_hardpart == ON ) then
       din_hsr   = rhoo_hsr ! chgqo
       dout_hsr  = rho_hsr  ! chgq
       rho_hsr(:,1) = rmx_this(1) *dout_hsr(:,1) + ( 1.0D0 -rmx_this(1) )*din_hsr(:,1)
       if (nspin_m==2) then
          rho_hsr(:,2) = rmx_this(2) *dout_hsr(:,2) + ( 1.0D0 -rmx_this(2) )*din_hsr(:,2)
       endif
       if ( noncol ) then
          rho_hsr(:,2) = rmx_this(2) *dout_hsr(:,2) + ( 1.0D0 -rmx_this(2) )*din_hsr(:,2)
          rho_hsr(:,3) = rmx_this(3) *dout_hsr(:,3) + ( 1.0D0 -rmx_this(3) )*din_hsr(:,3)
          rho_hsr(:,4) = rmx_this(4) *dout_hsr(:,4) + ( 1.0D0 -rmx_this(4) )*din_hsr(:,4)
       endif
    endif

    if ( sw_mix_charge_with_ekindens == ON ) then
      if(kgpm == kgp .or. npes == 1) then
          do is = 1, ndim_magmom, af+1
             din_ekinq (ista_kgpm:iend_kgpm,:,is) = ekinqo_l(ista_kgpm:iend_kgpm,:,is)
             dout_ekinq(ista_kgpm:iend_kgpm,:,is) = ekinq_l (ista_kgpm:iend_kgpm,:,is)
          end do
       else
!          call scatter_chg_onto_d( ekinqo_l,din )  ! -(m_C.D.)
!          call scatter_chg_onto_d( ekinq_l, dout)  ! -(m_C.D.)
          call scatter_chg_onto_d( ekinqo_l,din_ekinq )  ! -(m_C.D.)
          call scatter_chg_onto_d( ekinq_l, dout_ekinq)  ! -(m_C.D.)
       end if
       do is = 1, ndim_magmom, af+1
          do k = 1, kimg
             ekinq_l(:,k,is) = p2(:,is) *ekinq_l(:,k,is) &
                  &          + ( 1.0d0 -p2(:,is) ) *ekinqo_l(:,k,is)
          end do
       end do
    endif

  end subroutine simple_mix1_intg

  subroutine set_ncrspd_mxiter_etc_intg(iter,iuv,mxiter)
    integer, intent(in)  :: iuv,iter
    integer, intent(out) :: mxiter

    if(hownew == RENEW) then
       if((iter-istrbr+1) >= 3) then
          if((iter-istrbr+1) > nbxmix) then   ! When the box overflows
             call rotate_cmix_arrays   !-(contained here) ->mxiter,ncrspd,urec_l,f,g
          else
             mxiter = (iter-istrbr+1) - 1
             ncrspd(iter-istrbr+1) = iter-istrbr+1
          endif
       else
          mxiter = (iter-istrbr+1) - 1
          ncrspd(1) = 1;  ncrspd(2) = 2
       endif
    else ! if(hownew == ANEW)
       mxiter = (iter-istrbr+1) - 1
       ncrspd(iter-istrbr+1) = iter-istrbr+1
    endif

  contains

    subroutine rotate_cmix_arrays
      integer :: is,j,i,icr,jcr,iwork

      do is = 1, ndim_magmom, af+1
         do j = 3, nbxmix
            icr = ncrspd(2); jcr = ncrspd(j)
            g(j) = f(icr,jcr,is)
            do i = 3, j-1
               icr = ncrspd(i)
               g(j) = g(j) - f(icr,jcr,is)*g(i)
            enddo
            icr = ncrspd(2)

            urec_l(:,:,is,jcr,iuv) &
                 & = urec_l(:,:,is,jcr,iuv) + g(j)*urec_l(:,:,is,icr,iuv)

            if ( sw_mix_charge_hardpart == ON ) then
               urec_hsr(:,is,jcr,iuv) &
                    & = urec_hsr(:,is,jcr,iuv) + g(j)*urec_hsr(:,is,icr,iuv)
            endif
            if ( sw_mix_charge_with_ekindens == ON ) then
               urec_l_ekinq(:,:,is,jcr,iuv) &
                    & = urec_l_ekinq(:,:,is,jcr,iuv) + g(j)*urec_l_ekinq(:,:,is,icr,iuv)
            endif
         enddo
      enddo

      mxiter = nbxmix-1
      iwork = ncrspd(2)

      do i = 2, mxiter
         ncrspd(i)= ncrspd(i+1)
      end do
      ncrspd(mxiter+1) = iwork
    end subroutine rotate_cmix_arrays

  end subroutine set_ncrspd_mxiter_etc_intg

  subroutine renew_u_br_intg(j,i)
    integer, intent(in) :: j,i

    integer       :: is
    real(DP)      :: v_dF(nspin_m), v_dF_tmp(nspin_m)

    v_dF = 0.d0;  v_dF_tmp = 0.0d0

    call mult1s5( urec_l, nbxmix,2, j, iV, dF_l, f_p, v_dF )
    if ( sw_mix_charge_hardpart == ON ) then
       do is = 1, ndim_magmom, af+1
          v_dF(is) = v_dF(is) + sum( urec_hsr(:,is,j,iV)*dF_hsr(:,is) )
       End do
    endif
    if ( sw_mix_charge_with_ekindens == ON ) then
       call mult1s5( urec_l_ekinq, nbxmix,2, j, iV, dF_l_ekinq, f_p_ekinq, v_dF_tmp )
       v_dF = v_dF +v_df_tmp
    endif

    if ( nspin_m==2 .and. sw_mix_bothspins_sametime == YES ) then
      v_dF(1) = v_dF(1) + v_dF(2)
      v_df(2) = v_dF(1)
    endif
    if ( noncol ) then
       v_dF(1) = sum( v_dF(:) )
       v_dF(:) = v_dF(1)
    endif

    call subtr_j_th_term(v_dF,iU,j,urec_l,u_l)  !-(m_CD)
                                       !     |u(m)> = |u(m)> - <v(j)|dF(m)>|u(j)>
    if ( sw_mix_charge_hardpart == ON ) then
       do is = 1, ndim_magmom, af+1
          u_hsr(:,is) = u_hsr(:,is) - v_dF(is) *urec_hsr(:,is,j,iU)
       End do
    endif
    if ( sw_mix_charge_with_ekindens == ON ) then
       call subtr_j_th_term( v_dF,iU,j,urec_l_ekinq,u_l_ekinq )
    endif

    if(hownew == RENEW) f(j,i,1:nspin_m) = v_dF(1:nspin_m)

  end subroutine renew_u_br_intg

  subroutine renew_d_last_br_intg( p, p2 )
    real(DP), intent(in) :: p(ista_kngp:iend_kngp)
    real(DP), intent(in), optional :: p2(ista_kngp:iend_kngp)

    integer   :: is, ik, i, ns
    real(DP)  :: vF(nspin_m), vF_tmp(nspin_m)

    vF = 0.0d0; vF_tmp = 0.0d0

    call mult1s(v_l,F_l,f_p,vF)              !-(m_CD) <v|F> ->vF

    if ( sw_mix_charge_hardpart == ON ) then
       do is = 1, ndim_magmom, af+1
          vF(is) = vF(is) + sum( v_hsr(:,is)*FF_hsr(:,is) )
       End do
    endif
    if ( sw_mix_charge_with_ekindens == ON ) then
       call mult1s(v_l_ekinq,F_l_ekinq,f_p_ekinq,vF_tmp)              !-(m_CD) <v|F> ->vF
       vF = vF +vF_tmp
    endif

    if ( nspin_m==2 .and. sw_mix_bothspins_sametime == YES ) then
      vF(1) = vF(1) + vF(2)
      vF(2) = vF(1)
    endif

    if ( noncol ) then
       vF(1) = sum( vF(:) )
       vF(:) = vF(1)
    endif

! ---
    if(kgpm == kgp .or. npes == 1) then
       do is = 1, ndim_magmom, af+1
          din (:,:,is) = chgqo_l(ista_kgpm:iend_kgpm,:,is)
          dout(:,:,is) = chgq_l (ista_kgpm:iend_kgpm,:,is)
       end do
    else
       call scatter_chg_onto_d(chgqo_l,din )  ! -(m_C.D.)
       call scatter_chg_onto_d(chgq_l, dout)  ! -(m_C.D.)
    end if

    if ( sw_mix_charge_hardpart == ON ) then
       do is = 1, ndim_magmom, af+1
          din_hsr (:,is) = rhoo_hsr(:,is) ! chgqo
          dout_hsr(:,is) = rho_hsr (:,is) ! chgq
       end do
    endif

    if ( sw_mix_charge_with_ekindens == ON ) then
      if(kgpm == kgp .or. npes == 1) then
          do is = 1, ndim_magmom, af+1
             din_ekinq (:,:,is) = ekinqo_l(ista_kgpm:iend_kgpm,:,is)
             dout_ekinq(:,:,is) = ekinq_l (ista_kgpm:iend_kgpm,:,is)
          end do
       else
          call scatter_chg_onto_d( ekinqo_l,din )  ! -(m_C.D.)
          call scatter_chg_onto_d( ekinq_l, dout)  ! -(m_C.D.)
       end if
    endif

! ---
    ns = nspin_for_qnewton()
    do is = 1, ns,af+1
       do ik = 1, kimg
          do i = ista_kgpm,iend_kgpm
             rho(i,ik,is) = d0_l(i,ik,is) - vF(is)*u_l(i,ik,is)
          end do
       end do
    end do
    if ( sw_mix_charge_hardpart == ON ) then
       do is = 1, ns, af+1
          rho_hsr(:,is) = d0_hsr(:,is) - vF(is) *u_hsr(:,is)
       end do
    endif
    if ( sw_mix_charge_with_ekindens == ON ) then
       do is = 1, ns,af+1
          do ik = 1, kimg
             do i = ista_kgpm,iend_kgpm
                rho_ekinq(i,ik,is) = d0_l_ekinq(i,ik,is) - vF(is)*u_l_ekinq(i,ik,is)
             end do
          end do
       end do
    endif

    if(kgpm < kgp) then
       call concentrate_d_to_chg(rho,chgq_l) !-(m_C.D.)
       call simple_mix_large_Gc( p, chgqo_l, chgq_l )   !-(m_C.D.) chgq,chgqo,p ->chgq
    endif
    if ( sw_mix_charge_with_ekindens == ON ) then
       if(kgpm < kgp) then
          call concentrate_d_to_chg( rho_ekinq, ekinq_l ) !-(m_C.D.)
          call simple_mix_large_Gc( p2, ekinqo_l,ekinq_l )
       endif
    end if

  end subroutine renew_d_last_br_intg

  subroutine mix_broyden_dealloc2_intg
    deallocate(d0_l); deallocate(u_l); deallocate(v_l)
    if ( sw_mix_charge_hardpart == ON ) then
       deallocate(d0_hsr); deallocate(u_hsr); deallocate(v_hsr)
    endif
    if ( sw_mix_charge_with_ekindens == ON ) then
       deallocate(d0_l_ekinq); deallocate(u_l_ekinq); deallocate(v_l_ekinq)
    endif
    call dealloc_rho_rhoo_and_cpm_intg
  end subroutine mix_broyden_dealloc2_intg

  subroutine dealloc_rho_rhoo_and_cpm_intg
    if(kgpm < kgp .and. npes /= 1) then
       deallocate(rho); deallocate(rhoo); deallocate(c_pm)
    endif
    if ( sw_mix_charge_with_ekindens == ON ) then
      if(kgpm < kgp .and. npes /= 1) then
          deallocate(rho_ekinq); deallocate(rhoo_ekinq); deallocate(c_pm_ekinq)
       endif
    endif
  end subroutine dealloc_rho_rhoo_and_cpm_intg

  subroutine renew_d_br_intg(j)
    integer, intent(in) :: j

    real(DP)  :: vF(nspin_m), vF_tmp(nspin_m)
    integer :: is

    vF = 0.d0; vF_tmp = 0.0d0

    call mult1s5( urec_l, nbxmix, 2, j, iV, F_l, f_p, vF ) !-(m_CD);<v|F>  ->vF
    if ( sw_mix_charge_hardpart == ON ) then
       do is = 1, ndim_magmom, af+1
          vF(is) = vF(is) + sum( urec_hsr(:,is,j,iV)*FF_hsr(:,is) )
       End do
    endif
    if ( sw_mix_charge_with_ekindens == ON ) then
       call mult1s5( urec_l_ekinq, nbxmix,2, j, iV, F_l_ekinq, f_p_ekinq, vF_tmp )
       vF = vF +vF_tmp
    endif

    if ( nspin_m==2 .and. sw_mix_bothspins_sametime == YES ) then
      vF(1) = vF(1) + vF(2)
      vF(2) = vF(1)
    endif

    if ( noncol ) then
       vF(1) = sum( vF(:) )
       vF(:) = vF(1)
    endif

    call subtr_j_th_term(vF,iU,j,urec_l,d0_l) !-(m_CD)
    !                        |d(m)> = |d(m)> - <v(j)|F(m)>|u(j)>
    if ( sw_mix_charge_hardpart == ON ) then
       do is = 1, ndim_magmom, af+1
          d0_hsr(:,is) = d0_hsr(:,is) - vF(is) *urec_hsr(:,is,j,iU)
       end do
    endif
    if ( sw_mix_charge_with_ekindens == ON ) then
       call subtr_j_th_term( vF, iU, j, urec_l_ekinq, d0_l_ekinq )
    endif

  end subroutine renew_d_br_intg

  subroutine m_CD_mix_broyden2_intg(nfout,rmx,mixocc)
    integer, intent(in) :: nfout
    real(DP),intent(in) :: rmx
    logical, intent(in) :: mixocc
    integer   :: iter,j,mxiter,icr,jcr

    integer   :: id_sname = -1

    real(kind=DP), allocatable, dimension(:):: rmxtrc ! d(nspin_m)

    call tstatc0_begin('m_CD_mix_broyden2_intg ',id_sname,1)

! -- init ---
    if (previous_waymix /= BROYD2.or.force_dealloc) then

       if ( sw_mix_charge_hardpart == ON ) then
          if ( first ) then
             call create_map_func(.true.,mixocc)
             call alloc_rho_hsr(mixocc)
             call create_map_func(.false.,mixocc)
             first = .false.
          endif
       endif
       call mix_dealloc_previous_intg()
       call mix_broyden_allocate_intg();

       F_l => din
       if ( sw_mix_charge_hardpart == ON ) FF_hsr => din_hsr
       if ( sw_mix_charge_with_ekindens == ON ) F_l_ekinq => din_ekinq

       force_dealloc = .false.

    end if

    if ( sw_mix_charge_hardpart == ON ) then
       if ( noncol ) then
          call map_hsr_to_rho_noncl( hsr, hsi, rho_hsr )
          call map_hsr_to_rho_noncl( hsro,hsio,rhoo_hsr )
       else
          call map_hsr_to_rho( hsr, rho_hsr )
          call map_hsr_to_rho( hsro,rhoo_hsr )
          if(mixocc) then
            call map_om_to_rho(om, rho_hsr)
            call map_om_to_rho(omold, rho_hsr)
          endif
       endif
    endif

    allocate(rmxtrc(nspin_m))
    rmxtrc(1:nspin_m) = rmx
    if ( noncol ) rmxtrc(2:nspin_m) = min( rmx *spin_density_mixfactor, rmx_max )

    if ( sw_recomposing_hsr == YES .and. af == 0 .and. nspin == 2 ) then
       call alloc_chgqstore_recompose_chgq(rmx,rmxtrc) ! --> chgq_l, chgqo_l, rmxtrc
       if ( sw_mix_charge_hardpart==ON ) then
          call alloc_rhostore_recomp( rmx, rmxtrc )
       endif
       if ( sw_mix_charge_with_ekindens == ON ) then
          call alloc_kinqstore_recompos_kinq( rmx, rmxtrc )
       endif
    endif

    allocate(c_p(ista_kngp:iend_kngp,nspin_m)); c_p = 0
    if ( noncol ) then
       call precon_4_charge_mix_noncl(rmxtrc,c_p)
    else
       call precon_4_charge_mix(rmxtrc,c_p)
    endif

    if ( sw_mix_charge_with_ekindens == ON ) then
       allocate(c_p_ekinq(ista_kngp:iend_kngp,nspin_m)); c_p_ekinq = 0
       c_p_ekinq = rmx
    endif

! -- start ---
    iter = iter_from_reset()                 !-(m_CD)

    if((iter-istrbr+1) <= 1) then
       if ( sw_mix_charge_with_ekindens == ON ) then
          call simple_mix1_intg( c_p, rmxtrc, c_p_ekinq )                 !-(m_CD)
       else
          call simple_mix1_intg( c_p, rmxtrc )                 !-(m_CD)
       endif
    else

       call mix_broyden_alloc2_intg    !-(m_CD) d0_l,u_l, and v_l are allocated
       call dF_F_d0_u_and_v_intg

       call set_ncrspd_mxiter_etc_intg(iter,iU,mxiter) !-(m_CD) ->mxiter,ncrspd
       !                  when hownew == RENEW: f,g,ncrspd, and urec_l are reset.
       icr = icrspd_is(iter)                 !-(m_CD) function
       do j = 2, mxiter
          jcr = ncrspd(j)
          call renew_u_br_intg(jcr,icr) !-(m_CD) |u(m)> = |u(m)> - <v(j)|dF(m)>|u(j)>
          call renew_d_br_intg(jcr)     !-(m_CD) |d(m)> = |d(m)> - <v(j)|F(m)> |u(j)>
       enddo!j-loop

       urec_l(:,:,:,icr,iU) = u_l(:,:,:)  ! storing
       urec_l(:,:,:,icr,iV) = v_l(:,:,:)  ! storing

       if ( sw_mix_charge_hardpart == ON ) then
          urec_hsr(:,:,icr,iU) = u_hsr(:,:)  ! storing
          urec_hsr(:,:,icr,iV) = v_hsr(:,:)  ! storing
       endif
       if ( sw_mix_charge_with_ekindens == ON ) then
          urec_l_ekinq(:,:,:,icr,iU) = u_l_ekinq(:,:,:)  ! storing
          urec_l_ekinq(:,:,:,icr,iV) = v_l_ekinq(:,:,:)  ! storing
       endif

       call renew_d_last_br_intg( c_p, c_p_ekinq )
                                !-(m_CD) chgq_l(|d(m)>) = |d(m)>-<v(m)|F(m)>|u(m)>
       call mix_broyden_dealloc2_intg                      !-(m_CD)

    endif

    if (sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
       call compose_chgq_dealloc_chgqstore()
       if ( sw_mix_charge_hardpart == ON ) then
          call compose_rho_dealloc_store
       end if
       if ( sw_mix_charge_with_ekindens == ON ) then
          call compos_kinq_dealloc_kinqstore()
       endif
    endif

    if ( sw_mix_charge_hardpart == ON ) then
       if ( noncol ) then
          call map_rho_to_hsr_noncl( hsr, hsi, rho_hsr )
       else
          call map_rho_to_hsr( hsr, rho_hsr )
          if(mixocc) call map_rho_to_om(om,rho_hsr)
       endif
    endif

    deallocate(rmxtrc)
! =========================================================================== 11.0

    deallocate(c_p)

    previous_waymix = BROYD2
    call tstatc0_end(id_sname)
  contains

    subroutine dF_F_d0_u_and_v_intg
      !   dF_l(=deltaF) = (rho - dout) - (rhoo - din)
      !   F_l = rho - rhoo (=\cal F^{m}); u_l  = (rhoo - din) + c_p*dF_l;
      !   d0_l = rhoo+c_p* F_l;              v_l = dF_l/( |dF_l| )

      integer                      :: is,k,i
      real(DP), dimension(nspin_m) :: fff, fff_tmp

      do is = 1, ndim_magmom, af+1
         do k = 1, kimg
            do i = ista_kgpm,iend_kgpm
               dF_l(i,k,is) = (rho (i,k,is)-rhoo(i,k,is)) - (dout(i,k,is)-F_l(i,k,is))
               d0_l(i,k,is) = rhoo(i,k,is) + c_pm(i,is)*(rho(i,k,is) - rhoo(i,k,is))
               u_l(i,k,is)  = c_pm(i,is)*dF_l(i,k,is) + (rhoo(i,k,is) - F_l(i,k,is))
! ----
               F_l(i,k,is)  = rho(i,k,is) - rhoo(i,k,is)
            end do
            if(mype == 0) u_l(1,k,is) = 0.d0
        end do
      end do

      if ( sw_mix_charge_hardpart == ON ) then
         do is = 1, ndim_magmom, af+1
            dF_hsr(:,is) = ( rho_hsr(:,is)-rhoo_hsr(:,is)) &
                 &         - ( dout_hsr(:,is)-FF_hsr(:,is))
            d0_hsr(:,is) = rhoo_hsr(:,is) &
                 &        + rmxtrc(is) *( rho_hsr(:,is) - rhoo_hsr(:,is))
            u_hsr(:,is) = rmxtrc(is) *dF_hsr(:,is) + ( rhoo_hsr(:,is) - FF_hsr(:,is) )
            FF_hsr(:,is) = rho_hsr(:,is) - rhoo_hsr(:,is)
         end do
      endif

      if ( sw_mix_charge_with_ekindens == ON ) then
         do is = 1, ndim_magmom, af+1
            do k = 1, kimg
               do i = ista_kgpm,iend_kgpm
                  dF_l_ekinq(i,k,is) = (rho_ekinq (i,k,is)-rhoo_ekinq(i,k,is)) &
                       &               - (dout_ekinq(i,k,is) -F_l_ekinq(i,k,is))
                  d0_l_ekinq(i,k,is) = rhoo_ekinq(i,k,is) &
                       &              + c_pm_ekinq(i,is) &
                       &               *( rho_ekinq(i,k,is) -rhoo_ekinq(i,k,is))
                  u_l_ekinq(i,k,is)  = c_pm_ekinq(i,is)*dF_l_ekinq(i,k,is) &
                       &               + ( rhoo_ekinq(i,k,is) -F_l_ekinq(i,k,is) )
! ----
                  F_l_ekinq(i,k,is)  = rho_ekinq(i,k,is) - rhoo_ekinq(i,k,is)
               end do
               if(mype == 0) u_l_ekinq(1,k,is) = 0.d0
            end do
         end do
      endif

! - calc fff ---
      call mult1s(dF_l,dF_l,f_p,fff)
      if ( sw_mix_charge_hardpart == ON ) then
         do is = 1, ndim_magmom, af+1
            fff(is) = fff(is) + sum( dF_hsr(:,is)*dF_hsr(:,is) )
         end do
      endif
      if ( sw_mix_charge_with_ekindens == ON ) then
         call mult1s(dF_l_ekinq,dF_l_ekinq,f_p,fff_tmp)
      endif
      fff = fff + fff_tmp

      if(sum(fff) < 1.d-40)  call phase_error_with_msg(nfout,' fmult is too small',__LINE__,__FILE__)

      if ( nspin_m == 2 .and. sw_mix_bothspins_sametime == YES ) then
        fff(1) = fff(1) + fff(2)
        fff(2) = fff(1)
      endif
      if ( noncol ) then
         fff(1) = sum( fff(:) )
         fff(:) = fff(1)
      endif

! - calc v_l ---
      do is = 1, ndim_magmom, af+1
         v_l(:,:,is) = dF_l(:,:,is)/fff(is)
      end do

      if ( sw_mix_charge_hardpart == ON ) then
         do is = 1, ndim_magmom, af+1
            do i=1,nsize_rho_hsr
               v_hsr(i,is) = dF_hsr(i,is)/fff(is)
            end do
         end do
      endif
      if ( sw_mix_charge_with_ekindens == ON ) then
         do is = 1, ndim_magmom, af+1
            v_l_ekinq(:,:,is) = dF_l_ekinq(:,:,is)/fff(is)
         end do
      endif
    end subroutine dF_F_d0_u_and_v_intg

  end subroutine m_CD_mix_broyden2_intg

  subroutine mix_pulay_allocate_intg
    if ( noncol ) then
       nspin_m  = ndim_magmom
    else
       nspin_m  = nspin/(af+1)
    endif

    allocate(f_p(ista_kgpm:iend_kgpm)); f_p = 0;
    call precon_4_mult(f_p) !-(m_CD)

    allocate(din(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(dout(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(urec_l(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix,2))

    allocate(uuf_p(nbxmix,nspin_m))
    allocate(f(nbxmix,nbxmix,nspin_m))
    allocate(g_p(nbxmix,nspin_m))
    allocate(prj_wk(mp_kngp*kimg*nspin_m))
    allocate(ncrspd(nbxmix))

    if(sw_control_stepsize==ON) allocate(d0_l_h(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix))

    if(sw_gradient_simplex==ON)then
       allocate(rhoj(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix));rhoj=0.d0
       allocate(Frhoj(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix));Frhoj=0.d0
       allocate(rhojo(ista_kgpm:iend_kgpm,kimg,nspin_m));rhojo=0.d0
       allocate(Frhojo(ista_kgpm:iend_kgpm,kimg,nspin_m));Frhojo=0.d0
    endif

    allocate(ynorm(nbxmix,nspin_m));ynorm=1.d0
    din = 0.0d0; dout = 0.0d0; urec_l = 0.0d0; uuf_p = 0.0d0; f = 0.0d0
    g_p = 0.0d0; prj_wk = 0.0d0; ncrspd = 0

! --
    if ( sw_mix_charge_hardpart == ON ) then
       allocate(  din_hsr( nsize_rho_hsr,nspin_m ) );  din_hsr = 0.0d0
       allocate( dout_hsr( nsize_rho_hsr,nspin_m ) ); dout_hsr = 0.0d0
       allocate(   dF_hsr( nsize_rho_hsr,nspin_m ) );   dF_hsr = 0.0d0
       allocate( urec_hsr( nsize_rho_hsr,nspin_m,nbxmix,2) ); urec_hsr = 0.0d0
       if(sw_gradient_simplex==ON)then
          allocate(rhoj_hsr(nsize_rho_hsr,nspin_m,nbxmix));rhoj_hsr=0.d0
          allocate(Frhoj_hsr(nsize_rho_hsr,nspin_m,nbxmix));Frhoj_hsr=0.d0
          allocate(rhoj_hsro(nsize_rho_hsr,nspin_m));rhoj_hsro=0.d0
          allocate(Frhoj_hsro(nsize_rho_hsr,nspin_m));Frhoj_hsro=0.d0
       endif
       if(sw_control_stepsize==ON) then
          allocate(d0_hsr_h( nsize_rho_hsr,nspin_m,nbxmix) )
          d0_hsr_h = 0.0d0
       endif
    endif

    if ( sw_mix_charge_with_ekindens == ON ) then
       allocate(din_ekinq(ista_kgpm:iend_kgpm,kimg,nspin_m))
       allocate(dout_ekinq(ista_kgpm:iend_kgpm,kimg,nspin_m))
       allocate(urec_l_ekinq(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix,2))
       din_ekinq = 0.0d0;  dout_ekinq = 0.0d0;  urec_l_ekinq = 0.0d0

       if(sw_gradient_simplex==ON)then
          allocate(rhoj_ekinq(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix))
          allocate(Frhoj_ekinq(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix))
          allocate(rhojo_ekinq(ista_kgpm:iend_kgpm,kimg,nspin_m))
          allocate(Frhojo_ekinq(ista_kgpm:iend_kgpm,kimg,nspin_m))
          rhoj_ekinq = 0.0d0; Frhoj_ekinq = 0.0d0;
          rhojo_ekinq = 0.0d0; Frhojo_ekinq = 0.0d0
       endif
       if(sw_control_stepsize==ON) then
          allocate(d0_l_h_ekinq(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix))
          d0_l_h_ekinq = 0.0d0
       endif
       allocate(ynorm_ekinq(nbxmix,nspin_m)); ynorm_ekinq=1.d0

       allocate(f_p_ekinq(ista_kgpm:iend_kgpm)); f_p_ekinq = 1.0d0
    endif

  end subroutine mix_pulay_allocate_intg

  subroutine mix_pulay_alloc2_intg
    allocate(d0_l(ista_kgpm:iend_kgpm,kimg,nspin_m));  d0_l = 0.0d0
    if ( sw_mix_charge_hardpart == ON ) then
       allocate( d0_hsr( nsize_rho_hsr,nspin_m) ); d0_hsr = 0.0d0
    endif
    if ( sw_mix_charge_with_ekindens == ON ) then
       allocate(d0_l_ekinq(ista_kgpm:iend_kgpm,kimg,nspin_m));  d0_l_ekinq = 0.0d0
    endif

    call alloc_rho_rhoo_and_cpm_intg
  end subroutine mix_pulay_alloc2_intg

  subroutine mix_pulay_dealloc2_intg
    deallocate(d0_l)
    if ( sw_mix_charge_hardpart == ON ) deallocate(d0_hsr)
    if ( sw_mix_charge_with_ekindens == ON ) deallocate( d0_l_ekinq )

    call dealloc_rho_rhoo_and_cpm_intg

  end subroutine mix_pulay_dealloc2_intg

  subroutine m_CD_mix_pulay_intg(nfout,rmx,mixocc)
    integer, intent(in) :: nfout
    real(DP),intent(in) :: rmx
    logical, intent(in) :: mixocc

    integer, parameter  :: iRho = 1, iResid = 2
    integer   :: iter, mxiter
    real(DP),pointer,dimension(:)  :: e_wk, f_wk, ww1, finv
    integer, pointer,dimension(:)  :: ip

    real(kind=DP), allocatable, dimension(:):: rmxtrc ! d(nspin_m)
    real(kind=DP) :: rmxtt
    integer   :: id_sname = -1

    call tstatc0_begin('m_CD_mix_pulay_intg ',id_sname,1)

! -- init ---
    if (previous_waymix /= PULAY.or.force_dealloc) then

       if ( sw_mix_charge_hardpart == ON ) then
          if ( first ) then
             call set_i2lp_max2lp()
             call create_map_func(.true.,mixocc)
             call alloc_rho_hsr(mixocc)
             call create_map_func(.false.,mixocc)
             first = .false.
          endif
       endif
       call mix_dealloc_previous_intg()
       call mix_pulay_allocate_intg()

       force_dealloc = .false.
    end if

    if ( sw_mix_charge_hardpart == ON ) then
       if ( noncol ) then
          call map_hsr_to_rho_noncl( hsr, hsi, rho_hsr )
          call map_hsr_to_rho_noncl( hsro,hsio,rhoo_hsr )
       else
          call map_hsr_to_rho( hsr, rho_hsr )
          call map_hsr_to_rho( hsro,rhoo_hsr )
          if(mixocc) then
             call map_om_to_rho(om, rho_hsr)
             call map_om_to_rho(omold, rho_hsr)
          endif
       endif
    endif

    allocate(rmxtrc(nspin_m))
    rmxtrc(1:nspin_m) = rmx
    if ( noncol ) rmxtrc(2:nspin_m) = min( rmx *spin_density_mixfactor, rmx_max )

    if ( sw_recomposing_hsr == YES .and. af == 0 .and. nspin == 2 ) then
       call alloc_chgqstore_recompose_chgq(rmx,rmxtrc) ! --> chgq_l, chgqo_l, rmxtrc
       if ( sw_mix_charge_hardpart==ON ) then
          call alloc_rhostore_recomp( rmx, rmxtrc )
       endif
       if ( sw_mix_charge_with_ekindens == ON ) then
          call alloc_kinqstore_recompos_kinq( rmx, rmxtrc )
       endif
    endif

    if(sw_control_stepsize==ON)then
       rmxtt = rmx*step_control_factor
       if(rmxtt>max_stepsize) rmxtt = max_stepsize
       rmxtrc = rmxtt
       call m_CtrlP_set_rmx(rmxtt)
       if(printable) write(nfout,'(a,f10.5)') 'step size for the current iteration : ',rmxtt
    endif

    allocate(c_p(ista_kngp:iend_kngp,nspin_m)); c_p = 0
    if ( noncol ) then
       call precon_4_charge_mix_noncl(rmxtrc,c_p)
    else
       call precon_4_charge_mix(rmxtrc,c_p)
    endif

    if ( sw_mix_charge_with_ekindens == ON ) then
       allocate(c_p_ekinq(ista_kngp:iend_kngp,nspin_m)); c_p_ekinq = 0
       c_p_ekinq = rmx
    endif

! -- start ---
    iter = iter_from_reset()                 !-(m_CD)

    if(iter.eq.1) then
       alpha_pulay = alpha_pulay_org
    endif

    if((iter-istrbr+1) <= 1) then
       if ( sw_mix_charge_with_ekindens == ON ) then
          call simple_mix1_intg( c_p, ommix_factor*rmxtrc, c_p_ekinq )     !-(m_CD)
       else
          call simple_mix1_intg( c_p, ommix_factor*rmxtrc )                 !-(m_CD)
       endif
       !   din=chgqo_l; dout=chgq_l; (din,dout,c_p)->chgq_l

    else
       call mix_pulay_alloc2_intg   !-(m_CD) d0_l,u_l, and w_l are allocated
       call set_ncrspd_mxiter(nbxmix,iter-istrbr,mxiter) ! -> ncrspd, mxiter

       call Resid_and_dd_into_urec_intg(mxiter) !-(c.h.)
       call Ri_dot_Rj_intg(mxiter)          !-(c.h.) <R(i)|R(j)>->f
       call get_finv_lapack_intg(nbxmix,mxiter,f)  !-(c.h.) f -> f^{-1}= <R(i)|R(j)>^{-1}
       call Rj_dot_d_intg(mxiter)           !-(c.h.) <R(j)|d>,(j=1,iter-istrb) -> uuf_p
       call get_gmatrix_intg(mxiter)        !-(c.h.) (f,uuf_p)->g

       if ( sw_mix_charge_with_ekindens == ON ) then
          call renew_d_using_g_intg( mxiter, c_pm, ommix_factor*rmxtrc, c_pm_ekinq )
       else
          call renew_d_using_g_intg( mxiter, c_pm, ommix_factor*rmxtrc )     !-(c.h.)
       endif

       call mix_pulay_dealloc2_intg                    !-(m_CD)
    endif

    if (sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
       call compose_chgq_dealloc_chgqstore()
       if ( sw_mix_charge_hardpart == ON ) then
          call compose_rho_dealloc_store
       end if
       if ( sw_mix_charge_with_ekindens == ON ) then
          call compos_kinq_dealloc_kinqstore()
       endif
    endif

    if ( sw_mix_charge_hardpart == ON ) then
       if ( noncol ) then
          call map_rho_to_hsr_noncl( hsr, hsi, rho_hsr )
       else
          call map_rho_to_hsr( hsr, rho_hsr )
          if(mixocc) call map_rho_to_om(om,rho_hsr)
       endif
    endif

    deallocate(rmxtrc)
    deallocate(c_p)
    if ( sw_mix_charge_with_ekindens == ON ) then
       deallocate(c_p_ekinq)
    endif

    previous_waymix = PULAY
    call tstatc0_end(id_sname)

  contains

    subroutine set_ncrspd_mxiter(n,iter,m)
      integer, intent(in)  :: n, iter
      integer, intent(out) :: m
      integer :: i, nx
      if(hownew == ANEW) then
         m = iter
!         ncrspd(:) = (/(i,i=1,m)/)
         do i=1,m
           ncrspd(i) = i
         enddo
      else ! hownew == RENEW
         if(iter <= n) then
            m = iter
!            ncrspd(:) = (/(i,i=1,m)/)
            do i=1,m
              ncrspd(i) = i
            enddo
         else
            m = n
            nx = ncrspd(1)
            do i = 1, m-1
               ncrspd(i) = ncrspd(i+1)
            end do
            ncrspd(m) = nx
         end if
      end if
    end subroutine set_ncrspd_mxiter

    subroutine Resid_and_dd_into_urec_intg(iter)
      integer, intent(in) :: iter
      integer             :: itc,itc0,itc1
      integer :: i,j,k,nmix,imix,ierr
      real(kind=DP) :: sum1,sum2, ctmp2, ctmp3

      itc = ncrspd(iter)
      if(sw_gradient_simplex==ON)then
         do imix=2,iter-1
            itc0 = ncrspd(imix)
            itc1 = ncrspd(imix-1)
            do i=1,nspin_m
               do j=1,kimg
                  do k=ista_kgpm,iend_kgpm
                     urec_l(k,j,i,itc0,iResid) &
                          &        = rho(k,j,i) -rhoo(k,j,i) &
                          &        - (Frhoj(k,j,i,itc1) -rhoj(k,j,i,itc1) )
                     urec_l(k,j,i,itc0,iRho  ) = rhoo(k,j,i) -rhoj(k,j,i,itc1)
                  enddo
               enddo
            enddo
            if ( sw_mix_charge_hardpart == ON ) then
               do i=1,nspin_m
                  do j=1,nsize_rho_hsr
                     urec_hsr(j,i,itc0,iResid) &
                          &     = rho_hsr(j,i) -rhoo_hsr(j,i) &
                          &      - ( Frhoj_hsr(j,i,itc1) -rhoj_hsr(j,i,itc1) )
                     urec_hsr(j,i,itc0,iRho  ) = rhoo_hsr(j,i) -rhoj_hsr(j,i,itc1)
                  enddo
               enddo
            end if
            if ( sw_mix_charge_with_ekindens == ON ) then
               do i=1,nspin_m
                  do j=1,kimg
                     do k=ista_kgpm,iend_kgpm
                        urec_l_ekinq(k,j,i,itc0,iResid) &
                             &     = rho_ekinq(k,j,i) -rhoo_ekinq(k,j,i) &
                             &     - ( Frhoj_ekinq(k,j,i,itc1) -rhoj_ekinq(k,j,i,itc1) )
                        urec_l_ekinq(k,j,i,itc0,iRho  ) &
                             &     = rhoo_ekinq(k,j,i) -rhoj_ekinq(k,j,i,itc1)
                     enddo
                  enddo
               enddo
            endif
         end do

         do i=1,nspin_m
            do j=1,kimg
               do k=ista_kgpm,iend_kgpm
                  urec_l(k,j,i,itc,iResid) = rho(k,j,i)-rhoo(k,j,i) &
                       &                   -( dout(k,j,i) -din(k,j,i) )
                  urec_l(k,j,i,itc,iRho  ) = rhoo(k,j,i) -din(k,j,i)
                  rhoj(k,j,i,itc) = rhoo(k,j,i)
                  Frhoj(k,j,i,itc) = rho(k,j,i)
                  d0_l(k,j,i) = rho(k,j,i) - rhoo(k,j,i)
                  din(k,j,i)  = rhoo(k,j,i)
                  dout(k,j,i) = rho(k,j,i)
               enddo
            enddo
         enddo
         if ( sw_mix_charge_hardpart == ON ) then
            do i=1,nspin_m
               do j=1,nsize_rho_hsr
                  urec_hsr(j,i,itc,iResid) &
                       &    = rho_hsr(j,i) -rhoo_hsr(j,i) &
                       &     - (dout_hsr(j,i) -din_hsr(j,i)) ! =dF(=delta F^i)
                  urec_hsr(j,i,itc,iRho  ) = rhoo_hsr(j,i) -din_hsr(j,i)       ! =dd
                  d0_hsr(j,i) = rho_hsr(j,i) -rhoo_hsr(j,i)
                  din_hsr(j,i) = rhoo_hsr(j,i)
                  dout_hsr(j,i) = rho_hsr(j,i)
                  rhoj_hsr(j,i,itc) = rhoo_hsr(j,i)
                  Frhoj_hsr(j,i,itc) = rho_hsr(j,i)
               enddo
            enddo
         endif
         if ( sw_mix_charge_with_ekindens == ON ) then
            do i=1,nspin_m
               do j=1,kimg
                  do k=ista_kgpm,iend_kgpm
                     urec_l_ekinq(k,j,i,itc,iResid) &
                          &     = rho_ekinq(k,j,i) -rhoo_ekinq(k,j,i) &
                          &      -( dout_ekinq(k,j,i) -din_ekinq(k,j,i) )
                     urec_l_ekinq(k,j,i,itc,iRho  ) &
                          &     = rhoo_ekinq(k,j,i) -din_ekinq(k,j,i)
                     rhoj_ekinq(k,j,i,itc) = rhoo_ekinq(k,j,i)
                     Frhoj_ekinq(k,j,i,itc) = rho_ekinq(k,j,i)
                     d0_l_ekinq(k,j,i) = rho_ekinq(k,j,i) - rhoo_ekinq(k,j,i)
                     din_ekinq(k,j,i)  = rhoo_ekinq(k,j,i)
                     dout_ekinq(k,j,i) = rho_ekinq(k,j,i)
                  enddo
               enddo
            enddo
         endif

      else
         do i=1,nspin_m
            do j=1,kimg
               do k=ista_kgpm,iend_kgpm
                  urec_l(k,j,i,itc,iResid) = rho(k,j,i) -rhoo(k,j,i) &
                       &                    - (dout(k,j,i) -din(k,j,i)) ! =dF(=delta F^i)
                  urec_l(k,j,i,itc,iRho  ) = rhoo(k,j,i) -din(k,j,i)            ! =dd
                  d0_l(k,j,i) = rho(k,j,i) -rhoo(k,j,i)
                  din(k,j,i)  = rhoo(k,j,i)
                  dout(k,j,i) = rho(k,j,i)
               enddo
            enddo
         enddo
         if ( sw_mix_charge_hardpart == ON ) then
            do i=1,nspin_m
               do j=1,nsize_rho_hsr
                  urec_hsr(j,i,itc,iResid) &
                       &       = rho_hsr(j,i) -rhoo_hsr(j,i) &
                       &        - (dout_hsr(j,i) -din_hsr(j,i)) ! =dF(=delta F^i)
                  urec_hsr(j,i,itc,iRho  ) = rhoo_hsr(j,i) -din_hsr(j,i)       ! =dd
                  d0_hsr(j,i) = rho_hsr(j,i) -rhoo_hsr(j,i)
                  din_hsr(j,i) = rhoo_hsr(j,i)
                  dout_hsr(j,i) = rho_hsr(j,i)
               enddo
            enddo
         endif
         if ( sw_mix_charge_with_ekindens == ON ) then
            do i=1,nspin_m
               do j=1,kimg
                  do k=ista_kgpm,iend_kgpm
                     urec_l_ekinq(k,j,i,itc,iResid) &
                          &   = rho_ekinq(k,j,i) -rhoo_ekinq(k,j,i) &
                          &   - (dout_ekinq(k,j,i) -din_ekinq(k,j,i)) ! =dF(=delta F^i)
                     urec_l_ekinq(k,j,i,itc,iRho  ) &
                          &   = rhoo_ekinq(k,j,i) -din_ekinq(k,j,i)            ! =dd
                     d0_l_ekinq(k,j,i) = rho_ekinq(k,j,i) -rhoo_ekinq(k,j,i)
                     din_ekinq(k,j,i)  = rhoo_ekinq(k,j,i)
                     dout_ekinq(k,j,i) = rho_ekinq(k,j,i)
                  enddo
               enddo
            enddo
         endif
      endif

!!$      ynorm(itc,:)=0.d0
      do i=1,nspin_m
         sum12(i) = 0.d0
         do j=1,kimg
            do k=ista_kgpm,iend_kgpm
!!$               ynorm(itc,i) = ynorm(itc,i) &
!!$                    &        +urec_l(k,j,i,itc,iResid) *urec_l(k,j,i,itc,iResid)
               sum12(i) = sum12(i) + urec_l(k,j,i,itc,iResid) *urec_l(k,j,i,itc,iResid)
            enddo
         enddo
      enddo
!!$      call mpi_allreduce(MPI_IN_PLACE,ynorm(itc,1),nspin_m,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
      call mpi_allreduce(MPI_IN_PLACE,sum12,nspin_m,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
      ynorm(itc,:) = sum12(:)

      if ( sw_mix_charge_hardpart == ON ) then
         ynorm(itc,:) = univol*ynorm(itc,:)
         do i=1,nspin_m
            do j=1,nsize_rho_hsr
               ynorm(itc,i) = ynorm(itc,i) &
                    &        +urec_hsr(j,i,itc,iResid) *urec_hsr(j,i,itc,iResid)
            enddo
         enddo
      endif
      if ( sw_mix_charge_with_ekindens == ON ) then
!!$         ynorm_ekinq(itc,:)=0.d0
         do i=1,nspin_m
            sum12(i) = 0.d0
            do j=1,kimg
               do k=ista_kgpm,iend_kgpm
!!$                  ynorm_ekinq(itc,i) = ynorm(itc,i) &
!!$                       &            +urec_l_ekinq(k,j,i,itc,iResid) &
!!$                       &            *urec_l_ekinq(k,j,i,itc,iResid)
                  sum12(i) = sum12(i) + urec_l_ekinq(k,j,i,itc,iResid) &
                       &               *urec_l_ekinq(k,j,i,itc,iResid)
               enddo
            enddo
         enddo
!!$         call mpi_allreduce( MPI_IN_PLACE, ynorm_ekinq(itc,1), nspin_m, &
!!$              &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
!!$         ynorm = ynorm +ynorm_ekinq
         call mpi_allreduce( MPI_IN_PLACE, sum12, nspin_m, mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
         ynorm_ekinq(itc,:) = sum12(:)
         ynorm(itc,:) = ynorm(itc,:) +ynorm_ekinq(itc,:)
      endif

      ynorm(itc,:) = 1.d0/sqrt(ynorm(itc,:))

!---
      sum1 = 0.0d0
      do i=1,nspin_m
         do j=1,kimg
            do k=ista_kgpm,iend_kgpm
               sum1 = sum1 + d0_l(k,j,i)*d0_l(k,j,i)
            enddo
         enddo
      enddo
      call mpi_allreduce( MPI_IN_PLACE, sum1, 1, mpi_double_precision, &
           &             mpi_sum, MPI_CommGroup, ierr )

      if ( sw_mix_charge_hardpart == ON ) then
         ctmp2 = 0.0d0
         do i=1,nspin_m
            do j=1,nsize_rho_hsr
               ctmp2 = ctmp2 + d0_hsr(j,i)*d0_hsr(j,i)
            enddo
         enddo
         sum1 = sum1 +ctmp2
      endif
      if ( sw_mix_charge_with_ekindens == ON ) then
         ctmp3 = 0.0d0
         do i=1,nspin_m
            do j=1,kimg
               do k=ista_kgpm,iend_kgpm
                  ctmp3 = ctmp3 + d0_l_ekinq(k,j,i)*d0_l_ekinq(k,j,i)
               enddo
            enddo
         enddo
         call mpi_allreduce( MPI_IN_PLACE, ctmp3, 1, mpi_double_precision, &
              &             mpi_sum, MPI_CommGroup, ierr )
         sum1 = sum1 +ctmp3
      endif

      sum1 = dsqrt(sum1)

      if(sum1<alpha_pulay_damp_thres) then
         alpha_pulay = alpha_pulay*alpha_pulay_damp
         if(alpha_pulay>0.and.printable) &
         &       write(nfout,'(a,e12.5)') '!** new value for the parameter alpha : ', &
         &       alpha_pulay
      endif

      if(sw_control_stepsize==ON ) d0_l_h(:,:,:,itc) = d0_l(:,:,:)
      if ( sw_mix_charge_hardpart == ON ) then
         if (sw_control_stepsize==ON) d0_hsr_h(:,:,itc) = d0_hsr(:,:)
      endif
      if ( sw_mix_charge_with_ekindens == ON ) then
         if (sw_control_stepsize==ON) d0_l_h_ekinq(:,:,:,itc) = d0_l_ekinq(:,:,:)
      endif

      if(iter>=2.and.sw_control_stepsize==ON)then
         sum2 = 0.d0
         itc1 = ncrspd(iter-1)

         do i=1,nspin_m
            do j=1,kimg
               do k=ista_kgpm,iend_kgpm
                  sum2 = sum2 + d0_l_h(k,j,i,itc1)*d0_l_h(k,j,i,itc1)
               enddo
            enddo
         enddo
         call mpi_allreduce( MPI_IN_PLACE, sum2, 1, mpi_double_precision, &
              &              mpi_sum, MPI_CommGroup, ierr )

         if ( sw_mix_charge_hardpart == ON ) then
            ctmp2 = 0.0d0
            do i=1,nspin_m
               do j=1,nsize_rho_hsr
                  ctmp2 = ctmp2 + d0_hsr_h(j,i,itc1)*d0_hsr_h(j,i,itc1)
               enddo
            enddo
            sum2 = sum2 +ctmp2
         endif
         if ( sw_mix_charge_with_ekindens == ON ) then
            ctmp3 = 0.0d0
            do i=1,nspin_m
               do j=1,kimg
                  do k=ista_kgpm,iend_kgpm
                     ctmp3 = ctmp3 + d0_l_h_ekinq(k,j,i,itc1)*d0_l_h_ekinq(k,j,i,itc1)
                  enddo
               enddo
            enddo
            call mpi_allreduce( MPI_IN_PLACE, ctmp3, 1, mpi_double_precision, &
              &              mpi_sum, MPI_CommGroup, ierr )
            sum2 = sum2 +ctmp3
         endif

         sum2 = dsqrt(sum2)
         step_control_factor = max(0.5d0,min(2.0d0,sum2/sum1))
      endif
    end subroutine Resid_and_dd_into_urec_intg

    subroutine Ri_dot_Rj_intg(n)
      integer, intent(in) :: n

      integer  :: it,jt,itc,jtc,is
      real(DP) :: ff1(nspin_m),ff2(nspin_m),ff1tmp, ff2tmp

      do it = 1, n
         itc = ncrspd(it)
         do jt = it, n
            jtc = ncrspd(jt)

            if (sw_mix_bothspins_sametime==ON .and. nspin_m>1)then
               ff1tmp=0.d0
               call mult1s10_reduce_spin( urec_l, nbxmix, 2, itc, iResid, &
                    &                     urec_l, jtc, iResid, f_p, ff1tmp )
                                                 ! <delta F^i|delta F^j>
               if ( sw_mix_charge_hardpart == ON ) then
                  do is=1,nspin_m,af+1
                     ff1tmp = ff1tmp + sum( urec_hsr(:,is,itc,iResid) &
                          &                *urec_hsr(:,is,jtc,iResid) )
                  enddo
               endif
               if ( sw_mix_charge_with_ekindens == ON ) then
                  call mult1s10_reduce_spin( urec_l_ekinq, nbxmix, 2, itc, iResid, &
                       &                     urec_l_ekinq, jtc, iResid, f_p_ekinq, &
                       &                     ff2tmp )
                  ff1tmp = ff1tmp +ff2tmp
               endif

               ff1(1) = ff1tmp;  ff1(2) = ff1tmp

            else
               call mult1s10( urec_l, nbxmix, 2, itc, iResid, &
                    &         urec_l, jtc, iResid, f_p, ff1 )   ! <delta F^i|delta F^j>
               if ( sw_mix_charge_hardpart == ON ) then
                  do is=1,nspin_m,(af+1)
                     ff1(is) = ff1(is) +sum( urec_hsr(:,is,itc,iResid) &
                          &                 * urec_hsr(:,is,jtc,iResid) )
                  enddo
               endif
               if ( sw_mix_charge_with_ekindens == ON ) then
                  call mult1s10( urec_l_ekinq, nbxmix, 2, itc, iResid, &
                       &         urec_l_ekinq, jtc, iResid, f_p_ekinq, ff2 )
                  ff1 = ff1 +ff2
               endif
            endif

            if ( noncol ) then
               ff1tmp=0.d0
               call mult1s10_reduce_spin( urec_l, nbxmix, 2, itc, iResid, &
                    &                     urec_l, jtc, iResid, f_p, ff1tmp )
                                              ! <delta F^i|delta F^j>
               if ( sw_mix_charge_hardpart == ON ) then
                  do is=1,ndim_magmom
                     ff1tmp = ff1tmp +sum( urec_hsr(:,is,itc,iResid) &
                          &               *urec_hsr(:,is,jtc,iResid) )
                  enddo
               endif
               if ( sw_mix_charge_with_ekindens == ON ) then
                  call mult1s10_reduce_spin( urec_l_ekinq, nbxmix, 2, itc, iResid, &
                       &                     urec_l_ekinq, jtc, iResid, f_p_ekinq, &
                       &                     ff2tmp )
                  ff1tmp = ff1tmp +ff2tmp
               endif

               ff1(:) = ff1tmp
            endif

            f(it,jt,1:nspin_m) = ff1(1:nspin_m)
            if(jt /= it) f(jt,it,1:nspin_m) = f(it,jt,1:nspin_m)
         end do
      end do
    end subroutine Ri_dot_Rj_intg

    subroutine Rj_dot_d_intg(n)
      integer, intent(in) :: n

      integer  :: jt, jtc, is
      real(DP) :: ff1(nspin_m),ff2(nspin_m)
      real(DP) :: ff1tmp, ff2tmp

      do jt = 1, n
         jtc = ncrspd(jt)
         if (sw_mix_bothspins_sametime==ON .and. nspin_m>1)then
            ff1tmp=0.d0
            call mult1s5_reduce_spin( urec_l, nbxmix, 2, jtc, iResid, &
                 &                    d0_l, f_p, ff1tmp )
            if ( sw_mix_charge_hardpart == ON ) then
               do is=1,nspin_m,af+1
                  ff1tmp = ff1tmp+sum(urec_hsr(:,is,jtc,iResid) * d0_hsr(:,is))
               enddo
            endif
            if ( sw_mix_charge_with_ekindens == ON ) then
               ff2tmp=0.d0
               call mult1s5_reduce_spin( urec_l_ekinq, nbxmix, 2, jtc, iResid, &
                    &                    d0_l_ekinq, f_p_ekinq, ff2tmp )
               ff1tmp = ff1tmp +ff2tmp
            endif

            ff1(1) = ff1tmp; ff1(2) = ff1tmp

         else
            call mult1s5( urec_l, nbxmix, 2, jtc, iResid, d0_l, f_p, ff1 )
            if ( sw_mix_charge_hardpart == ON ) then
               ff2 = 0.d0
               do is=1,nspin_m,af+1
                  ff2(is) = ff2(is) +sum( urec_hsr(:,is,jtc,iResid) *d0_hsr(:,is) )
               enddo
               ff1(:) = ff1(:) +ff2(:)
            endif
            if ( sw_mix_charge_with_ekindens == ON ) then
               ff2 = 0.0d0
               call mult1s5( urec_l_ekinq, nbxmix, 2, jtc, iResid, &
                    &        d0_l_ekinq, f_p_ekinq, ff2 )
               ff1(:) = ff1(:) +ff2(:)
            endif
         endif

         if ( noncol ) then
            ff1tmp=0.d0
            call mult1s5_reduce_spin( urec_l, nbxmix, 2, jtc, iResid, &
                 &                    d0_l, f_p, ff1tmp )
            if ( sw_mix_charge_hardpart == ON ) then
               do is=1,ndim_magmom
                  ff1tmp = ff1tmp+sum(urec_hsr(:,is,jtc,iResid) * d0_hsr(:,is))
               enddo
            endif
            if ( sw_mix_charge_with_ekindens == ON ) then
               ff2tmp=0.d0
               call mult1s5_reduce_spin( urec_l_ekinq, nbxmix, 2, jtc, iResid, &
                    &                    d0_l_ekinq, f_p_ekinq, ff2tmp )
               ff1tmp = ff1tmp +ff2tmp
            endif

            ff1(:) = ff1tmp
         endif

         uuf_p(jt,1:nspin_m) = ff1(1:nspin_m)
      end do
    end subroutine Rj_dot_d_intg

    subroutine get_finv_lapack_intg(m,n,f)
      integer,intent(in)                             :: m,n
      real(DP),intent(inout),dimension(m,m,nspin_m) :: f

      real(DP), allocatable,dimension(:,:) :: fwork
      integer :: is,inf,it,jt,kt,nnspin
      real(DP) :: div,tmp

      allocate(fwork(n,n))
      nnspin = nspin
      if(sw_mix_bothspins_sametime==ON .or. af==1) nnspin=1

! ============================== added by K. Tagami ============== 11.0
      if ( noncol )  nnspin = 1
! ================================================================ 11.0

      do is=1,nnspin
         if(ipripulay >= 2) then
            write(nfout,600) n,(('(',it,jt,')',f(it,jt,is),jt=1,n),it=1,n)
600         format(//11x,"**input matrix**"/12x &
                 & ,"horder=",I5/(2x,4(1x,1a,i2,",",i2,1a,e14.6)))
         end if
         fwork=0
         do it=1,n
            do jt=1,n
               fwork(jt,it) = f(jt,it,is)*ynorm(jt,is)*ynorm(it,is)
               if(it==jt) fwork(jt,it)=fwork(jt,it)+alpha_pulay
            enddo
         enddo
         call dpotrf('U',n,fwork,n,inf)
         call dpotri('U',n,fwork,n,inf)
         do it=1,n-1
            do jt=it+1,n
               fwork(jt,it) = fwork(it,jt)
            enddo
         enddo
         do it=1,n
            do jt=1,n
               f(jt,it,is) = fwork(jt,it)*ynorm(jt,is)*ynorm(it,is)
            enddo
         enddo
         if(ipripulay >= 2) then
            write(nfout,630) (('(',it,jt,')',f(it,jt,is),it=1,n),jt=1,n)
630         format(/11x, "**inverse matrix**" &
                 & ,/(2x,4(1x,1a,i2,",",i2,1a,e14.6)))
         end if
      enddo
      if(sw_mix_bothspins_sametime==ON .and. nspin_m>1)then
         do it=1,n
            do jt=1,n
               f(jt,it,2) = f(jt,it,1)
            enddo
         enddo
      endif

! ================================= added by K. Tagami ================= 11.0
      if ( noncol ) then
         do it=1,n
            do jt=1,n
               f(jt,it,:) = f(jt,it,1)
            enddo
         end do
      endif
! ====================================================================== 11.0
      deallocate(fwork)
    end subroutine get_finv_lapack_intg

    subroutine get_gmatrix_intg(n)
      integer,intent(in) :: n
      integer :: is, it, jt, nnspin
      nnspin = nspin
      if(sw_mix_bothspins_sametime==ON .or. af==1) nnspin=1

! ============================== added by K. Tagami ============== 11.0
      if ( noncol )  nnspin = 1
! ================================================================ 11.0

      g_p = 0.d0
      do is = 1, nnspin
         do it = 1, n
            do jt = 1, n
               g_p(it,is) = g_p(it,is) - f(jt,it,is)*uuf_p(jt,is)
            end do
         end do
         if(ipripulay >= 2) then
            write(nfout,'(" -- g_p(1:",i3,") --")') n
            write(nfout,'(8f20.12)') (g_p(it,is),it=1,n)
         end if
      end do
      if(sw_mix_bothspins_sametime==ON .and. nspin_m>1)then
         do it = 1,n
            g_p(it,2) = g_p(it,1)
         enddo
      endif
! ================================= added by K. Tagami ================= 11.0
      if ( noncol ) then
         do it = 1,n
            g_p(it,:) = g_p(it,1)
         end do
      end if
! ====================================================================== 11.0
    end subroutine get_gmatrix_intg

    subroutine renew_d_using_g_intg(n,p,rmx_this,p2)
      integer, intent(in)                                :: n
      real(DP),intent(in) :: p(ista_kgpm:iend_kgpm,nspin_m)
      real(DP),intent(in) :: rmx_this(nspin_m)
      real(DP),intent(in), optional :: p2(ista_kgpm:iend_kgpm,nspin_m)

      integer    :: is, k, i, it, itc, ns

!!$      do is = 1, nspin, af+1
      ns = nspin_for_qnewton()
      do is = 1, ns,af+1
         do k = 1, kimg
            do i = ista_kngp, iend_kngp
               rho(i,k,is)  = rhoo(i,k,is) + p(i,is)*d0_l(i,k,is)
            end do
            do it = 1, n
               itc = ncrspd(it)
               do i = ista_kngp, iend_kngp
                  rho(i,k,is) = rho(i,k,is) &
                       &      + g_p(it,is) * ( urec_l(i,k,is,itc,iRho) &
                       &                     + p(i,is)*urec_l(i,k,is,itc,iResid) )
               end do
            end do
         end do
      end do
      if ( sw_mix_charge_hardpart == ON ) then
         do is=1,ns,af+1
            rho_hsr(:,is) = rhoo_hsr(:,is) &
                 &        + rmx_this(is) * d0_hsr(:,is)
            do it = 1, n
               itc = ncrspd(it)
               rho_hsr(:,is) = rho_hsr(:,is) &
                    &        + g_p(it,is) * ( urec_hsr(:,is,itc,iRho) &
                    &                        + rmx_this(is) *urec_hsr(:,is,itc,iResid) )
            enddo
         enddo
      endif
      if ( sw_mix_charge_with_ekindens == ON ) then
         do is = 1, ns,af+1
            do k = 1, kimg
               do i = ista_kngp, iend_kngp
                  rho_ekinq(i,k,is)  = rhoo_ekinq(i,k,is) + p2(i,is)*d0_l_ekinq(i,k,is)
               end do
               do it = 1, n
                  itc = ncrspd(it)
                  do i = ista_kngp, iend_kngp
                     rho_ekinq(i,k,is) = rho_ekinq(i,k,is) &
                          &      + g_p(it,is) &
                          &          * ( urec_l_ekinq(i,k,is,itc,iRho) &
                          &             + p2(i,is)*urec_l_ekinq(i,k,is,itc,iResid) )
                  end do
               end do
            end do
         end do
      endif

      if(kgpm < kgp) then
         call concentrate_d_to_chg(rho,chgq_l) !-(m_C.D.)
         call simple_mix_large_Gc(c_p,chgqo_l,chgq_l)  !-(m_C.D.) chgq,chgqo,c_p ->chgq
      end if
      if ( sw_mix_charge_with_ekindens == ON ) then
         if(kgpm < kgp) then
            call concentrate_d_to_chg( rho_ekinq, ekinq_l )
            call simple_mix_large_Gc( c_p_ekinq, ekinqo_l, ekinq_l )
         end if
      endif
    end subroutine renew_d_using_g_intg

  end subroutine m_CD_mix_pulay_intg

  subroutine m_CD_simple_mixing_intg(nfout,rmxt)
    integer ,intent(in)      :: nfout
    real(kind=DP),intent(in) :: rmxt

    integer       :: is, k, ia, it, lmt1, lmt2
    real(kind=DP) :: rmxtt
    integer       :: id_sname = -1

    if ( noncol ) then
       nspin_m = ndim_magmom
    else
       nspin_m  = nspin/(af+1)
    endif

    call tstatc0_begin('m_CD_simple_mixing_intg ',id_sname,1)

    if(previous_waymix /= SIMPLE.or.force_dealloc) then
       call mix_dealloc_previous_intg()
       force_dealloc = .false.
    end if

    allocate(rmxtrc(nspin_m))
    rmxtrc = rmxt
    if ( noncol ) rmxtrc(2:nspin_m) = min( rmxt *spin_density_mixfactor, rmx_max )

    if ( sw_recomposing_hsr == YES .and. af == 0 .and. nspin == 2 ) then
       call alloc_chgqstore_recompose_chgq(rmxt,rmxtrc) ! --> chgq_l, chgqo_l, rmxtrc
       if ( sw_mix_charge_hardpart==ON ) then
          call alloc_rhostore_recomp( rmxt, rmxtrc )
       endif
       if ( sw_mix_charge_with_ekindens == ON ) then
          call alloc_kinqstore_recompos_kinq( rmxt, rmxtrc )
       endif
    endif

    if(ipri >= 2) write(nfout,'(" rmxt = ",d20.8)') rmxt

    allocate(c_p(ista_kngp:iend_kngp,nspin_m))
    c_p = 0.0d0
    if ( noncol ) then
       call precon_4_charge_mix_noncl(rmxtrc,c_p)
    else
       call precon_4_charge_mix(rmxtrc,c_p)
    endif
    if ( sw_mix_charge_with_ekindens == ON ) then
       allocate(c_p_ekinq(ista_kngp:iend_kngp,nspin_m)); c_p_ekinq = 0
       c_p_ekinq = rmxt
    endif

    do is = 1, ndim_magmom, af+1
       do k = 1, kimg
          chgq_l(:,k,is) = c_p(:,is)*chgq_l(:,k,is) + (1.0d0-c_p(:,is))*chgqo_l(:,k,is)
       end do
    end do

    if ( sw_mix_charge_hardpart == ON ) then
       do ia = 1, natm
          it = ityp(ia)
          if(ipaw(it)/=1) cycle
          if ( noncol ) then
             do is = 1, ndim_magmom
                do lmt1 = 1, ilmt(it)
                   do lmt2 = lmt1, ilmt(it)
                      hsr(ia,lmt1,lmt2,is) = rmxt*hsr(ia,lmt1,lmt2,is) + &
                           (1-rmxt)*hsro(ia,lmt1,lmt2,is)
                   end do! lmt2
                end do! lmt1
             end do! is
             if ( sw_mix_imaginary_hardpart == ON ) then
                do is = 1, ndim_magmom
                   do lmt1 = 1, ilmt(it)
                      do lmt2 = lmt1, ilmt(it)
                         hsi(ia,lmt1,lmt2,is) = rmxt*hsi(ia,lmt1,lmt2,is) &
                              &               +(1-rmxt)*hsio(ia,lmt1,lmt2,is)
                      end do! lmt2
                   end do! lmt1
                end do
             end if
          else
             do is = 1, nspin, af+1
                do lmt1 = 1, ilmt(it)
                   do lmt2 = lmt1, ilmt(it)
                      hsr(ia,lmt1,lmt2,is) = rmxt*hsr(ia,lmt1,lmt2,is) + &
                           (1-rmxt)*hsro(ia,lmt1,lmt2,is)
                   end do! lmt2
                end do! lmt1
             end do! is
          endif
       end do! ia
    endif

    if ( sw_mix_charge_with_ekindens == ON ) then
       do is = 1, ndim_magmom, af+1
          do k = 1, kimg
             ekinq_l(:,k,is) = c_p_ekinq(:,is) *ekinq_l(:,k,is) &
                  &          + (1.0d0-c_p_ekinq(:,is)) *ekinqo_l(:,k,is)
          end do
       end do
    endif

    if (sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
       call compose_chgq_dealloc_chgqstore()
       if ( sw_mix_charge_hardpart == ON ) then
          call compose_rho_dealloc_store
       end if
       if ( sw_mix_charge_with_ekindens == ON ) then
          call compos_kinq_dealloc_kinqstore()
       endif
    endif

    deallocate(c_p)
    if ( sw_mix_charge_with_ekindens == ON ) deallocate(c_p_ekinq)
    deallocate(rmxtrc)

    previous_waymix = SIMPLE
    call tstatc0_end(id_sname)

  end subroutine m_CD_simple_mixing_intg
! ====== 2014/09/19

  subroutine m_CD_hsr_diff(nfout)
    integer, intent(in) :: nfout
    integer :: i,j,ndiag,nnondiag
    real(kind=DP) :: sumhsr_diag,sumhsr_nondiag
    sumhsr_diag = 0.d0
    sumhsr_nondiag = 0.d0
    ndiag = 0
    nnondiag = 0
    do i=1,nspin_m
       do j=1,nsize_rho_hsr
          if(diag_elem(j))then
             sumhsr_diag = sumhsr_diag+abs(rhoo_hsr(j,i)-rho_hsr(j,i))
             ndiag = ndiag+1
          else
             sumhsr_nondiag = sumhsr_nondiag+abs(rhoo_hsr(j,i)-rho_hsr(j,i))
             nnondiag = nnondiag+1
          endif
       enddo
    enddo
    if(printable) write(nfout,'(a,f15.10)') '!** dhsr_diag   ',sumhsr_diag/dble(ndiag)
    if(printable) write(nfout,'(a,f15.10)') '!** dhsr_nondiag',sumhsr_nondiag/dble(nnondiag)
  end subroutine m_CD_hsr_diff

  subroutine alloc_kinqstore_recompos_kinq(rmxt,rmxtrc)
    real(kind=DP),intent(in) :: rmxt
    real(kind=DP),intent(out),dimension(nspin_m) :: rmxtrc

    allocate(ekinqstore_l(ista_kngp:iend_kngp,kimg,nspin))
    allocate(ekinqostore_l(ista_kngp:iend_kngp,kimg,nspin))
    ekinqstore_l = ekinq_l
    ekinqostore_l = ekinqo_l

    ekinq_l(:,:,1)  = ekinqstore_l(:,:,1)  + ekinqstore_l(:,:,2)
    ekinq_l(:,:,2)  = ekinqstore_l(:,:,1)  - ekinqstore_l(:,:,2)
    ekinqo_l(:,:,1) = ekinqostore_l(:,:,1) + ekinqostore_l(:,:,2)
    ekinqo_l(:,:,2) = ekinqostore_l(:,:,1) - ekinqostore_l(:,:,2)

    rmxtrc(1) = rmxt
    rmxtrc(2) = rmxt*spin_density_mixfactor
  end subroutine alloc_kinqstore_recompos_kinq

  subroutine compos_kinq_dealloc_kinqstore
    ekinqstore_l = ekinq_l
    ekinq_l(:,:,1) = 0.5d0*(ekinqstore_l(:,:,1) + ekinqstore_l(:,:,2))
    ekinq_l(:,:,2) = 0.5d0*(ekinqstore_l(:,:,1) - ekinqstore_l(:,:,2))
    ekinqo_l = ekinqostore_l
    deallocate(ekinqostore_l, ekinqstore_l)
  end subroutine compos_kinq_dealloc_kinqstore

end module m_CD_mixing
