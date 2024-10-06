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
#   define __TIMER_COMM_START_w_BARRIER(str,a)   call timer_barrier(str) ;   call timer_st
a(a)
#   define __TIMER_COMM_STOP_w_BARRIER(str,a)    call timer_barrier(str) ;   call timer_en
d(a)
#   define __TIMER_COMM_START(a)       call timer_sta(a)
#   define __TIMER_COMM_STOP(a)        call timer_end(a)
#else
#   define __TIMER_COMM_START_w_BARRIER(str,a)
#   define __TIMER_COMM_STOP_w_BARRIER(str,a)
#   define __TIMER_COMM_START(a)
#   define __TIMER_COMM_STOP(a)
#endif

module m_KE_mixing
! $Id: m_KE_mixing.F90 285 2013-01-01 04:33:41Z ktagami $
!
  use m_Const_Parameters,    only : BUCS, DP, OFF &
       &                          , EXECUT,SIMPLE_CUBIC,BOHR,NO,ANTIFERRO &
       &                          , ANEW,RENEW,ON, SIMPLE,BROYD1,BROYD2,DFP,PULAY &
       &                          , OLD, NEXT, PAI, VTK &
       &                          , DELTA10 &
       &                          , unit_conv_byname, UMICRO, GAMMA, DELTA &
       &                          , ELECTRON, INVERSE, YES, CMPLDP
  use m_IterationNumbers,    only : iteration,iteration_for_cmix
  use m_Parallelization,     only : MPI_CommGroup &
       &                          , ista_kngp,iend_kngp,is_kngp,ie_kngp,np_kngp,mp_kngp &
       &                          , npes,mype,ierr &
       &                          , is_kgpm,ie_kgpm,ista_kgpm,iend_kgpm,mp_kgpm &
       &                          , nis_fftp, nie_fftp, myrank_g, nrank_g

  use m_Timing,              only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters,  only : nspin,ipri,ipriwf, c_precon &
       &                          , amix,bmix,hownew,nbxmix,istrbr &
       &                          , kimg,af,neg,ipripulay &
       &                          , iprichargemixing &
       &                          , sw_recomposing, spin_density_mixfactor &
       &                          , amin, sw_precon_diff, sw_metric_diff,metric_ratio &
       &                          , sw_force_simple_mixing,printable, &
       &                            g0_wf_precon, amin_wf_precon

  use m_Crystal_Structure,   only : univol, sw_magnetic_constraint
  use m_PlaneWaveBasisSet,   only : kg,kgp,ngpt_l,ngabc,gr_l,kgpm


  use m_KineticEnergy_Density,  only : ekins_l, ekins_old, ekina_l, ekina_old

  use m_Control_Parameters,  only : sw_mix_bothspins_sametime, &
                                  & sw_recomposing_hsr, sw_force_simple_mixing_hsr, &
                                  & noncol, ndim_magmom

  use m_Ionic_System,        only : ityp, natm
  use m_PseudoPotential,     only : ilmt, nlmt

  use m_Control_Parameters,  only : sw_gradient_simplex, alpha_pulay, &
       &                            use_symm_ekin_density, use_asymm_ekin_density
  use m_CD_mixing,    only : ekinq_l, ekinqo_l
  use mpi

  implicit none

  real(kind=DP), allocatable, dimension(:):: rmxtrc ! d(nspin_m)
  real(kind=DP), allocatable, dimension(:,:,:) :: kinqstore_l, kinqostore_l

  real(kind=DP),private,pointer, dimension(:,:,:) ::       rho,rhoo  ! MPI
  !         rho => vlhxc_l, rhoo => vlhxco_l ( when kgpm == kgp)
  !         rho and rhoo is projected from vlhxc_l and vlhxco_l, respectively ( otherwise)

  real(kind=DP),private,allocatable,dimension(:,:,:,:) :: rhoj,Frhoj
  real(kind=DP),private,allocatable,dimension(:,:,:)   :: rhojo,Frhojo
  real(kind=DP),private,allocatable,dimension(:,:,:)   :: rhoj_vhsr,Frhoj_vhsr
  real(kind=DP),private,allocatable,dimension(:,:)   :: rhoj_vhsro,Frhoj_vhsro

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
!!  real(DP),private,pointer,dimension(:,:)     :: c_p !d(ista_kngp:iend_kngp,nspin/(af+1))
  real(DP),public,pointer,dimension(:,:)     :: c_p !d(ista_kngp:iend_kngp,nspin/(af+1))

  real(DP),private,pointer,dimension(:,:)     :: c_pm !d(ista_kgpm:iend_kgpm,nspin/(af+1))
! <--

  ! -- For Broyden and DFP mixing method --
  integer, private,parameter                :: iU = 1, iVD = 1, iW = 2, iY = 2, iV = 2
  integer, private                          :: nspin_m
  real(DP),private,allocatable,dimension(:) :: f_p !d(ista_kgpm:iend_kgpm)
  real(DP),private,allocatable,dimension(:,:,:) :: d0_l,u_l,v_l,w_l,dout,dd_l
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

  logical          :: force_dealloc = .false.
  integer, private :: previous_waymix = 0

!  include 'mpif.h'
  integer istatus(mpi_status_size)

! ========================== adde by K. Tagami ========================== 5.0
  integer :: nsize_rho_vhsr
  integer, private, allocatable :: imap_vhsr(:)    ! d(nsize_rho_vhsr)
  real(kind=DP),private,allocatable, dimension(:,:) ::   rho_vhsr, rhoo_vhsr
                                                  ! d(nsize_rho_vhsr,nspin)

  real(DP),private,allocatable,dimension(:,:) :: d0_vhsr, u_vhsr, v_vhsr, w_vhsr, &
        &                                        dout_vhsr, dd_vhsr
  real(DP),private,pointer,dimension(:,:)         :: FF_vhsr
  real(DP),private,allocatable,target,dimension(:,:) :: din_vhsr
  !                                             d( nsize_rho_vhsr,nspin_m)
  real(DP),private,allocatable,dimension(:,:)     :: dF_vhsr
  real(DP),private,allocatable,target,dimension(:,:,:,:) :: urec_vhsr

  logical, save :: first = .true.

  real(kind=DP), allocatable :: rho_store(:,:), vhsr(:,:,:,:)
  real(kind=DP), allocatable :: rhoo_store(:,:), vhsro(:,:,:,:)
! ==============================================================

! ========================== adde by K. Tagami ========================== 11.0
  integer :: nsize_rho_vhsr_unit
!  integer :: sw_mix_imaginary_hardpart = OFF
!  integer :: sw_mix_imaginary_hardpart = ON
! ======================================================================= 11.0

! ====
  real(kind=DP), private, pointer :: kinq_l(:,:,:), kinqo_l(:,:,:)
! ===

! --- contained subroutines ---
!   7. m_KE_prepare_precon       <-(KE_Mixing)
!  10. m_KE_simple_mixing        <-(KE_Mixing)
!  22. precon_4_KE_mix       <-(@10), (@53),(@54),(@55),(@60)
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
!  53. m_KE_mix_broyden1         <-(KE_Mixing)
!    - dF_F_d0_u_v_and_dd - renew_v
!  54. m_KE_mix_broyden2         <-(KE_Mixing)
!    - dF_F_d0_u_and_v
!  55. m_KE_mix_DFP              <-(KE_Mixing)
!    - dF_F_d0_u_and_w - renew_w - renew_d - renew_d_last
!  56. mix_pulay_allocate        <-(@60)
!  57. mix_pulay_deallocate      <-(@36)
!  58. mix_pulay_alloc2          <-(@60)
!  59. mix_pulay_dealloc2        <-(@60)
!  60. m_KE_mix_pulay            <-(KE_Mixing)
!    - mix_pulay_alloc3 - mix_pulay_dealloc3
!    - Resid_and_dd_into_urec - Ri_dot_Rj - Rj_dot_d
!    - get_finv -get_matrix -renew_d_using_g

contains

! ------------------------------------

  subroutine m_KE_set_pointer_ekinq
    if ( use_symm_ekin_density ) then
       ekinq_l => ekins_l
       ekinqo_l => ekins_old
    else if ( use_asymm_ekin_density ) then
       ekinq_l => ekina_l
       ekinqo_l => ekina_old
    endif
  end subroutine m_KE_set_pointer_ekinq

  subroutine m_KE_set_pointer_ekinq_00
    if ( use_symm_ekin_density ) then
       kinq_l => ekins_l
       kinqo_l => ekins_old
    else if ( use_asymm_ekin_density ) then
       kinq_l => ekina_l
       kinqo_l => ekina_old
    endif
  end subroutine m_KE_set_pointer_ekinq_00

  subroutine alloc_kinqstore_recompos_kinq(rmxt,rmxtrc)
    real(kind=DP),intent(in) :: rmxt
    real(kind=DP),intent(out),dimension(nspin_m) :: rmxtrc
#ifdef __TIMER_SUB__
    call timer_sta(1104)
#endif
    allocate(kinqstore_l(ista_kngp:iend_kngp,kimg,nspin))
    allocate(kinqostore_l(ista_kngp:iend_kngp,kimg,nspin))
#ifdef __TIMER_DO__
    call timer_sta(1149)
#endif
    kinqstore_l = kinq_l
    kinqostore_l = kinqo_l

    kinq_l(:,:,1)  = kinqstore_l(:,:,1)  + kinqstore_l(:,:,2)
    kinq_l(:,:,2)  = kinqstore_l(:,:,1)  - kinqstore_l(:,:,2)
    kinqo_l(:,:,1) = kinqostore_l(:,:,1) + kinqostore_l(:,:,2)
    kinqo_l(:,:,2) = kinqostore_l(:,:,1) - kinqostore_l(:,:,2)

#ifdef __TIMER_DO__
  call timer_end(1149)
#endif
    rmxtrc(1) = rmxt
    rmxtrc(2) = rmxt*spin_density_mixfactor
#ifdef __TIMER_SUB__
    call timer_end(1104)
#endif
  end subroutine alloc_kinqstore_recompos_kinq

  subroutine compos_kinq_dealloc_kinqstore
#ifdef __TIMER_SUB__
    call timer_sta(1106)
#endif
#ifdef __TIMER_DO__
  call timer_sta(1153)
#endif
  kinqstore_l = kinq_l
  kinq_l(:,:,1) = 0.5d0*(kinqstore_l(:,:,1) + kinqstore_l(:,:,2))
  kinq_l(:,:,2) = 0.5d0*(kinqstore_l(:,:,1) - kinqstore_l(:,:,2))
  kinqo_l = kinqostore_l

#ifdef __TIMER_DO__
  call timer_end(1153)
#endif
    deallocate(kinqostore_l, kinqstore_l)
#ifdef __TIMER_SUB__
    call timer_end(1106)
#endif
  end subroutine compos_kinq_dealloc_kinqstore

  subroutine m_KE_prepare_precon(nfout,rmxt)
    integer, intent(in)      :: nfout
    real(kind=DP), intent(in):: rmxt
    real(kind=DP) :: G_longest, G_shortest, x, gg,gg_mpi
    integer       :: i, n
    real(kind=DP) :: G_longest_mpi, G_shortest_mpi
    integer       :: ist !mpi

    if(param_cprecon_decided) return
#ifdef __TIMER_SUB__
    call timer_sta(1102)
#endif

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
    if(bmix_cprec < 0.d0) then
       q0 = q0_default
    else
       q0 = bmix_cprec*G_shortest
    end if
    q0 = q0*q0
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
       if(npes > 1) then
          call mpi_allreduce(gg,gg_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
          gg = gg_mpi
       end if

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
#ifdef __TIMER_SUB__
    call timer_end(1102)
#endif
  end subroutine m_KE_prepare_precon

  subroutine m_KE_force_dealloc()
    force_dealloc = .true.
  end subroutine m_KE_force_dealloc

  subroutine m_KE_simple_mixing(nfout,rmxt)
    integer ,intent(in)      :: nfout
    real(kind=DP),intent(in) :: rmxt

    integer       :: is, k
    integer       :: id_sname = -1
#ifdef __TIMER_SUB__
    call timer_sta(1103)
#endif

! ================================ modified by K. Tagami =============== 11.0
!!! --> T. Yamasaki  03 Aug. 2009
!!    nspin_m  = nspin/(af+1)
!!! <--
!
    if ( noncol ) then
       nspin_m = ndim_magmom
    else
       nspin_m  = nspin/(af+1)
    endif
! ====================================================================== 11.0

    call tstatc0_begin('m_KE_simple_mixing ',id_sname,1)

    if(previous_waymix /= SIMPLE.or.force_dealloc) then
       call mix_dealloc_previous()
! ------------------------------  ktDEBUG -------------------- 20121030
       call mix_dealloc_previous_vhsr()
! ------------------------------  ktDEBUG -------------------- 20121030
       force_dealloc = .false.
    end if

! ================================ modified by K. Tagami =============== 11.0
!! --> T. Yamasaki  03 Aug. 2009
!    allocate(rmxtrc(nspin_m))
!
!    if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
!       call alloc_kinqstore_recompos_kinq(rmxt,rmxtrc) ! --> kinq_l, kinqo_l, rmxtrc
!    else
!       rmxtrc = rmxt
!    end if
!    if(ipri >= 2) write(nfout,'(" rmxt = ",d20.8)') rmxt
!! --> T. Yamasaki  03 Aug. 2009
!
!
    allocate(rmxtrc(nspin_m))

    if ( noncol ) then
       rmxtrc = rmxt
    else
       if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
          call alloc_kinqstore_recompos_kinq(rmxt,rmxtrc) ! --> kinq_l, kinqo_l, rmxtrc
       else
          rmxtrc = rmxt
       endif
    end if
    if(ipri >= 2) write(nfout,'(" rmxt = ",d20.8)') rmxt
! ====================================================================== 11.0

!!$    allocate(c_p(ista_kngp:iend_kngp))
    allocate(c_p(ista_kngp:iend_kngp,nspin_m))
    c_p = 0.0d0                 ! ===== Adde by K. Tagami =========


! ================================ modified by K. Tagami =============== 11.0
!!    call precon_4_KE_mix(rmxtrc,c_p)
!
    if ( noncol ) then
       call precon_4_KE_mix_noncl(rmxtrc,c_p)
    else
       call precon_4_KE_mix(rmxtrc,c_p)
    endif
! ======================================================================= 11.0


#ifdef __TIMER_DO__
  call timer_sta(1148)
#endif
! ================================ modified by K. Tagami ================ 11.0
!!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ======================================================================= 11.0
       do k = 1, kimg
          kinq_l(:,k,is) = c_p(:,is)*kinq_l(:,k,is) &
               &  + (1.0d0-c_p(:,is))*kinqo_l(:,k,is)
       end do
    end do

#ifdef __TIMER_DO__
  call timer_end(1148)
#endif
    deallocate(c_p)

!
    if ( .not. noncol ) then
       if (sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
          call compos_kinq_dealloc_kinqstore()
       end if
    endif
    deallocate(rmxtrc)

    previous_waymix = SIMPLE
    call tstatc0_end(id_sname)

#ifdef __TIMER_SUB__
    call timer_end(1103)
#endif
  end subroutine m_KE_simple_mixing

  subroutine precon_4_KE_mix(pmix,c_p)
    real(DP),intent(in),dimension(nspin_m)                       :: pmix
    real(DP),intent(out), dimension(ista_kngp:iend_kngp,nspin_m) :: c_p
    integer              :: i
    real(DP)             :: gg, agg, tmp
    integer              :: ist  !for mpi
    integer              :: is
#ifdef __TIMER_SUB__
    call timer_sta(1105)
#endif

    do is = 1, nspin, af+1
       c_p(:,is) = pmix(is)
    end do
    return


    if(iprichargemixing >= 2) write(6,'("! pmix(precon_4_KE_mix) = ",f8.4)') pmix
    if(c_precon) then
       ist = ista_kngp
       if(ist == 1) ist = 2
       if(nspin_m == 2) then
#ifdef __TIMER_DO__
  call timer_sta(1150)
#endif
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
#ifdef __TIMER_DO__
  call timer_end(1150)
#endif
          if(mype == 0) then
             c_p(1,1) = amix_cprec*fg2 * pmix(1)
             if (sw_recomposing==ON .and. sw_precon_diff==NO)then
                c_p(1,2) = amix_cprec * pmix(2)
             else
                c_p(1,2) = amix_cprec*fg2 * pmix(2)
             endif
          end if
       else
#ifdef __TIMER_DO__
  call timer_sta(1151)
#endif
          do i = ist, iend_kngp  !for mpi
             gg = gr_l(i)*gr_l(i)
             tmp = max(gg/(gg+q0),amin)
!!$             c_p(i,1) = amix_cprec*gg/(gg+q0) * pmix(1)
             c_p(i,1) = amix_cprec* tmp * pmix(1)
          end do
#ifdef __TIMER_DO__
  call timer_end(1151)
#endif
          if(mype == 0) c_p(1,1) = amix_cprec*fg2 * pmix(1)
       end if
    else
#ifdef __TIMER_DO__
  call timer_sta(1152)
#endif
       do is = 1, nspin, af+1
          c_p(:,is) = pmix(is)
       end do
#ifdef __TIMER_DO__
  call timer_end(1152)
#endif
    endif
#ifdef __TIMER_SUB__
    call timer_end(1105)
#endif
  end subroutine precon_4_KE_mix

! ===================================== added by K. Tagami ================ 11.0
  subroutine precon_4_KE_mix_noncl(pmix,c_p)
    real(DP),intent(in),dimension(nspin_m)                       :: pmix
    real(DP),intent(out), dimension(ista_kngp:iend_kngp,nspin_m) :: c_p
    integer              :: i
    real(DP)             :: gg, agg, tmp
    integer              :: ist  !for mpi
    integer              :: is

    if(iprichargemixing >= 2) then
       write(6,'("! pmix(precon_4_KE_mix_noncl) = ",f8.4)') pmix
    end if

    if(c_precon) then
       ist = ista_kngp
       if(ist == 1) ist = 2

       do i = ist, iend_kngp  !for mpi
          gg = gr_l(i)*gr_l(i)
          tmp = max(gg/(gg+q0),amin)
!!$             agg = amix_cprec*gg/(gg+q0)
          agg = amix_cprec*tmp
          c_p(i,:) = agg * pmix(:)
       end do
       if (mype == 0) then
          c_p(1,:) = amix_cprec*fg2 * pmix(:)
       end if
    else
       do is = 1, ndim_magmom
          c_p(:,is) = pmix(is)
       end do
    endif

  end subroutine precon_4_KE_mix_noncl
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

  subroutine mult1s(u,v,f_q,fmult)
    real(DP),intent(in), dimension(ista_kgpm:iend_kgpm,kimg,nspin_m) :: u,v
    real(DP),intent(in), dimension(ista_kgpm:iend_kgpm):: f_q
    real(DP),intent(out),dimension(nspin_m)            :: fmult

    real(DP) :: p, p_mpi, fac
    integer  :: is,ik,i
#ifdef __TIMER_SUB__
    call timer_sta(1114)
#endif

    fmult = 0.d0

! ================================ modified by K. Tagami ============== 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ===================================================================== 11.0

       p = 0.d0
       fac=1.0d0
#ifdef __TIMER_DO__
  call timer_sta(1160)
#endif
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
#ifdef __TIMER_DO__
  call timer_end(1160)
#endif
       if( npes >= 2) then
          call mpi_allreduce(p,p_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
          p = p_mpi
       end if
       fmult(is) = p*univol
    enddo
#ifdef __TIMER_SUB__
    call timer_end(1114)
#endif
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

#ifdef __TIMER_SUB__
    call timer_sta(1115)
#endif
    fmult = 0.d0
! ================================ modified by K. Tagami ============== 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ===================================================================== 11.0

       p = 0.d0
       fac=1.0d0
#ifdef __TIMER_DO__
  call timer_sta(1162)
#endif
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
#ifdef __TIMER_DO__
  call timer_end(1162)
#endif

       if( npes >= 2) then
          call mpi_allreduce(p,p_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
          p = p_mpi
       end if
       fmult(is) = p*univol
    enddo
#ifdef __TIMER_SUB__
    call timer_sta(1115)
#endif
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
#ifdef __TIMER_SUB__
    call timer_sta(1137)
#endif

    fmult = 0.d0

! ====================================== modified by K. Tagami =========== 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ======================================================================== 11.0

       p = 0.d0
       fac = 1.0d0
#ifdef __TIMER_DO__
  call timer_sta(1186)
#endif
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
#ifdef __TIMER_DO__
  call timer_end(1186)
#endif

       if( npes >= 2) then
          call mpi_allreduce(p,p_mpi,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
          p = p_mpi
       end if
       fmult(is) = p*univol
    enddo
#ifdef __TIMER_SUB__
    call timer_end(1137)
#endif
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

  subroutine subtr_j_th_term(f,iuv,j,um)
! === DEBUG by tkato 2011/11/24 ================================================
!   real(DP), intent(in),   dimension(nspin) :: f
    real(DP), intent(in),   dimension(nspin_m) :: f
! ==============================================================================
    integer,  intent(in)                     :: iuv,j
    real(DP), intent(inout)                  :: um(ista_kgpm:iend_kgpm,kimg,nspin_m)
    integer :: is, ik, i, istart
#ifdef __TIMER_SUB__
    call timer_sta(1116)
#endif

    istart = ista_kgpm
    if(istart == 1) istart = 2
#ifdef __TIMER_DO__
  call timer_sta(1164)
#endif

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
#ifdef __TIMER_DO__
  call timer_end(1164)
#endif
#ifdef __TIMER_SUB__
    call timer_end(1116)
#endif
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
#ifdef __TIMER_SUB__
    call timer_sta(1119)
#endif

#ifdef __TIMER_DO__
  call timer_sta(1165)
#endif
! ====================================== modified by K. Tagami ============= 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ========================================================================== 11.0

       dv = 1.d0/f(is)
       urec_l(:,:,is,j,iuv) = v(:,:,is)*dv
    end do
#ifdef __TIMER_DO__
  call timer_end(1165)
#endif
#ifdef __TIMER_SUB__
    call timer_end(1119)
#endif
  end subroutine store_to_urec2

  subroutine set_ncrspd_mxiter_etc(iter,iuv,mxiter)
    integer, intent(in)  :: iuv,iter
    integer, intent(out) :: mxiter
#ifdef __TIMER_SUB__
    call timer_sta(1111)
#endif
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

#ifdef __TIMER_SUB__
    call timer_end(1111)
#endif
  contains
    subroutine rotate_cmix_arrays
      integer :: is,j,i,icr,jcr,iwork
#ifdef __TIMER_SUB__
    call timer_sta(1112)
#endif

#ifdef __TIMER_DO__
  call timer_sta(1158)
#endif
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
#ifdef __TIMER_DO__
  call timer_end(1158)
#endif
      mxiter = nbxmix-1
      iwork = ncrspd(2)
#ifdef __TIMER_DO__
  call timer_sta(1159)
#endif
      do i = 2, mxiter
         ncrspd(i)= ncrspd(i+1)
      end do
#ifdef __TIMER_DO__
  call timer_end(1159)
#endif
      ncrspd(mxiter+1) = iwork
#ifdef __TIMER_SUB__
    call timer_end(1112)
#endif
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
       do is = 2, 2
          din (ista_kgpm:iend_kgpm,:,is) = kinqo_l(ista_kgpm:iend_kgpm,:,is)
          dout(ista_kgpm:iend_kgpm,:,is) = kinq_l (ista_kgpm:iend_kgpm,:,is)
       end do
    else
       call scatter_chg_onto_d(kinqo_l,din )  ! -(m_C.D.)
       call scatter_chg_onto_d(kinq_l, dout)  ! -(m_C.D.)
    end if

    do is = 2,2
       do k = 1, kimg
          kinq_l(:,k,is) = p(:,is)*kinq_l(:,k,is) + (1.0d0-p(:,is))*kinqo_l(:,k,is)
       end do
    end do
  end subroutine simple_mix2

  subroutine simple_mix1(p)
    real(DP), intent(in), dimension(ista_kngp:iend_kngp,nspin_m) :: p
    integer  :: is,k
#ifdef __TIMER_SUB__
    call timer_sta(1108)
#endif

    if(kgpm == kgp .or. npes == 1) then
#ifdef __TIMER_DO__
  call timer_sta(1154)
#endif
! ========================= modified by K. Tagami ==================== 11.0
!       do is = 1, nspin, af+1
       do is = 1, ndim_magmom, af+1
! ==================================================================== 11.0

          din (ista_kgpm:iend_kgpm,:,is) = kinqo_l(ista_kgpm:iend_kgpm,:,is)
          dout(ista_kgpm:iend_kgpm,:,is) = kinq_l (ista_kgpm:iend_kgpm,:,is)
       end do
#ifdef __TIMER_DO__
  call timer_end(1154)
#endif
    else
       call scatter_chg_onto_d(kinqo_l,din )  ! -(m_C.D.)
       call scatter_chg_onto_d(kinq_l, dout)  ! -(m_C.D.)
    end if

#ifdef __TIMER_DO__
  call timer_sta(1155)
#endif
! ========================= modified by K. Tagami ==================== 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ==================================================================== 11.0
       do k = 1, kimg
          kinq_l(:,k,is) = p(:,is)*kinq_l(:,k,is) + (1.0d0-p(:,is))*kinqo_l(:,k,is)
       end do
    end do
#ifdef __TIMER_DO__
  call timer_end(1155)
#endif

#ifdef __TIMER_SUB__
    call timer_end(1108)
#endif
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
#ifdef __TIMER_SUB__
    call timer_sta(1121)
#endif

    nelmnt = mp_kngp*kimg*nspin_m
    prj_wk = 0.d0
    do ip = 0, npes - 1
       ! (1)  coping input data onto a work array, and broadcasting
       if(is_kngp(ip) > kgpm) exit
#ifdef __TIMER_DO__
  call timer_sta(1168)
#endif
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
#ifdef __TIMER_DO__
  call timer_end(1168)
#endif
       call mpi_bcast(prj_wk,nelmnt,mpi_double_precision,ip,MPI_CommGroup,ierr)

       ! (2) projection
       istart = ista_kgpm; if(istart < is_kngp(ip)) istart = is_kngp(ip)
       iend   = iend_kgpm; if(iend   > ie_kngp(ip)) iend   = ie_kngp(ip)
       if(iend < istart) cycle
#ifdef __TIMER_DO__
  call timer_sta(1170)
#endif
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
#ifdef __TIMER_DO__
  call timer_end(1170)
#endif
    end do
#ifdef __TIMER_SUB__
    call timer_end(1121)
#endif
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

#ifdef __TIMER_SUB__
    call timer_sta(1122)
#endif
    if(kgpm < kgp .and. npes /= 1) then
       nelmnt = mp_kgpm*kimg*nspin_m
       do ip = 0, npes - 1
          if(is_kngp(ip) > kgpm) exit
          do ip2 = 0, npes - 1
             istart = is_kgpm(ip2); if(istart < is_kngp(ip)) istart = is_kngp(ip)
             iend   = ie_kgpm(ip2); if(iend   > ie_kngp(ip)) iend   = ie_kngp(ip)
             if(iend < istart) cycle
             if(mype == ip2) then

#ifdef __TIMER_DO__
  call timer_sta(1171)
#endif
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
#ifdef __TIMER_SUB__
    call timer_end(1122)
#endif
  end subroutine concentrate_d_to_chg

  subroutine mix_dealloc_previous()
    if(previous_waymix == BROYD2) then
       call mix_broyden_deallocate()
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
!    allocate(f_p(ista_kgpm:iend_kgpm)); call precon_4_mult(f_p) !-(m_KE)
    allocate(f_p(ista_kgpm:iend_kgpm)); f_p = 0.0d0
    call precon_4_mult(f_p) !-(m_KE)
    f_p = 1.0d0
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
       call scatter_chg_onto_d(kinq_l,rho)
       call scatter_chg_onto_d(kinqo_l,rhoo)
       call scatter_cp_onto_cpm(c_p,c_pm)
    else
       rho => kinq_l; rhoo => kinqo_l; c_pm => c_p
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
#ifdef __TIMER_SUB__
    call timer_sta(1109)
#endif
    allocate(d0_l(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(u_l(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(v_l(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(dd_l(ista_kgpm:iend_kgpm,kimg,nspin_m))
! =========================================== Added by K. Tagami ===
    d0_l = 0.0d0; u_l = 0.0d0; v_l = 0.0d0; dd_l = 0.0d0
! =================================================================
    call alloc_rho_rhoo_and_cpm
#ifdef __TIMER_SUB__
    call timer_end(1109)
#endif
  end subroutine mix_broyden_alloc3

  subroutine mix_broyden_dealloc3
    deallocate(d0_l); deallocate(u_l); deallocate(v_l); deallocate(dd_l)
    call dealloc_rho_rhoo_and_cpm
  end subroutine mix_broyden_dealloc3

!!! ========================== added by K. Tagami ================ 5.0
  subroutine alloc_rhostore_recomp( rmxt, rmxtrc )
    real(kind=DP),intent(in) :: rmxt
    real(kind=DP),intent(out),dimension(nspin_m) :: rmxtrc

    allocate( rhoo_store( nsize_rho_vhsr,nspin) )
    allocate( rho_store( nsize_rho_vhsr,nspin) )

    rho_store = rho_vhsr;      rhoo_store = rhoo_vhsr

     rho_vhsr(:,1) =  rho_store(:,1) +  rho_store(:,2)
     rho_vhsr(:,2) =  rho_store(:,1) -  rho_store(:,2)
    rhoo_vhsr(:,1) = rhoo_store(:,1) + rhoo_store(:,2)
    rhoo_vhsr(:,2) = rhoo_store(:,1) - rhoo_store(:,2)
    rmxtrc(1) = rmxt;     rmxtrc(2) = rmxt*spin_density_mixfactor

  end subroutine alloc_rhostore_recomp

  subroutine compose_rho_dealloc_store
    rho_store = rho_vhsr

    rho_vhsr(:,1) = 0.5d0*( rho_store(:,1) + rho_store(:,2) )
    rho_vhsr(:,2) = 0.5d0*( rho_store(:,1) - rho_store(:,2) )

    rhoo_vhsr(:,1) = 0.5d0*( rhoo_store(:,1) + rhoo_store(:,2) )
    rhoo_vhsr(:,2) = 0.5d0*( rhoo_store(:,1) - rhoo_store(:,2) )

!    rhoo_vhsr = rhoo_store
    deallocate( rho_store, rhoo_store )

  end subroutine compose_rho_dealloc_store

 subroutine simple_mix_kt(rmx_this)
    real(kind=DP), intent(in) :: rmx_this(nspin_m)

    din_vhsr   = rhoo_vhsr ! kinqo
    dout_vhsr  = rho_vhsr  ! kinq
    rho_vhsr(:,1) = rmx_this(1) *dout_vhsr(:,1) + ( 1.0D0 - rmx_this(1) )* din_vhsr(:,1)
    if(nspin_m==2) rho_vhsr(:,2) = rmx_this(2) *dout_vhsr(:,2) + ( 1.0D0 - rmx_this(2) )* din_vhsr(:,2)

! ===================== added by K. Tagami ================================ 11.0
    if ( noncol ) then
       rho_vhsr(:,2) = rmx_this(2) *dout_vhsr(:,2) + ( 1.0D0 -rmx_this(2) )* din_vhsr(:,2)
       rho_vhsr(:,3) = rmx_this(3) *dout_vhsr(:,3) + ( 1.0D0 -rmx_this(3) )* din_vhsr(:,3)
       rho_vhsr(:,4) = rmx_this(4) *dout_vhsr(:,4) + ( 1.0D0 -rmx_this(4) )* din_vhsr(:,4)
    endif
! ========================================================================= 11.0

  end subroutine simple_mix_kt

  subroutine simple_mix2_kt(rmx_this)
    real(kind=DP), intent(in) :: rmx_this(nspin_m)

    if(nspin_m==2)then
        din_vhsr(:,2)  = rhoo_vhsr(:,2)         ! kinqo
       dout_vhsr(:,2)  =  rho_vhsr(:,2)          ! kinq
       rho_vhsr(:,2) = rmx_this(2) *dout_vhsr(:,2) + ( 1.0D0 - rmx_this(2) )* din_vhsr(:,2)
    endif
! ===================== added by K. Tagami ================================ 11.0
    if ( noncol ) then
       din_vhsr(:,2)  = rhoo_vhsr(:,2);     dout_vhsr(:,2)  =  rho_vhsr(:,2)
       rho_vhsr(:,2) = rmx_this(2) *dout_vhsr(:,2) + ( 1.0D0 -rmx_this(2) )* din_vhsr(:,2)
       din_vhsr(:,3)  = rhoo_vhsr(:,3);     dout_vhsr(:,3)  =  rho_vhsr(:,3)
       rho_vhsr(:,3) = rmx_this(3) *dout_vhsr(:,3) + ( 1.0D0 -rmx_this(3) )* din_vhsr(:,3)
       din_vhsr(:,4)  = rhoo_vhsr(:,4);     dout_vhsr(:,4)  =  rho_vhsr(:,4)
       rho_vhsr(:,4) = rmx_this(4) *dout_vhsr(:,4) + ( 1.0D0 -rmx_this(4) )* din_vhsr(:,4)
    endif
! ========================================================================== 11.0

  end subroutine simple_mix2_kt

  subroutine mix_dealloc_previous_vhsr()
    if(previous_waymix == BROYD2) then
       call mix_broyden_deallocate_vhsr()
    else if(previous_waymix == PULAY) then
       call mix_PULAY_deallocate_vhsr()
    end if
  end subroutine mix_dealloc_previous_vhsr

  subroutine mix_broyden_allocate_vhsr
! ========================= modified by K. Tagami =================== 11.0
!    nspin_m  = nspin/(af+1)
!
    if ( noncol ) then
       nspin_m  = ndim_magmom
    else
       nspin_m  = nspin/(af+1)
    endif
! ====================================================================== 11.0
    allocate( din_vhsr(nsize_rho_vhsr,nspin_m) )
    allocate( dout_vhsr(nsize_rho_vhsr,nspin_m))
    allocate( dF_vhsr(nsize_rho_vhsr,nspin_m))
    allocate( urec_vhsr(nsize_rho_vhsr,nspin_m,nbxmix,2) )
!
    din_vhsr = 0.0d0; dout_vhsr = 0.0d0; dF_vhsr = 0.0d0; urec_vhsr = 0.0d0
!
  end subroutine mix_broyden_allocate_vhsr

  subroutine mix_broyden_deallocate_vhsr
    if(allocated(din_vhsr)) deallocate(din_vhsr)
    if(allocated(dout_vhsr)) deallocate(dout_vhsr)
    if(allocated(dF_vhsr)) deallocate(dF_vhsr)
    if(allocated(urec_vhsr)) deallocate(urec_vhsr)
  end subroutine mix_broyden_deallocate_vhsr

  subroutine mix_broyden_alloc2_vhsr
    allocate( d0_vhsr( nsize_rho_vhsr,nspin_m ) ); d0_vhsr = 0.0d0
    allocate(  u_vhsr( nsize_rho_vhsr,nspin_m ) );  u_vhsr = 0.0d0
    allocate(  v_vhsr( nsize_rho_vhsr,nspin_m ) );  v_vhsr = 0.0d0
  end subroutine mix_broyden_alloc2_vhsr

 subroutine mix_broyden_dealloc2_vhsr
    deallocate(d0_vhsr); deallocate(u_vhsr); deallocate(v_vhsr)
  end subroutine mix_broyden_dealloc2_vhsr

  subroutine mix_broyden_alloc3_vhsr
    allocate( d0_vhsr( nsize_rho_vhsr,nspin_m ) ); d0_vhsr = 0.0d0
    allocate(  u_vhsr( nsize_rho_vhsr,nspin_m ) );  u_vhsr = 0.0d0
    allocate(  v_vhsr( nsize_rho_vhsr,nspin_m ) );  v_vhsr = 0.0d0
    allocate( dd_vhsr( nsize_rho_vhsr,nspin_m ) ); dd_vhsr = 0.0d0
  end subroutine mix_broyden_alloc3_vhsr

  subroutine mix_broyden_dealloc3_vhsr
    deallocate(d0_vhsr); deallocate(u_vhsr); deallocate(v_vhsr); deallocate(dd_vhsr)
  end subroutine mix_broyden_dealloc3_vhsr

 subroutine mix_pulay_allocate_vhsr
! ========================= modified by K. Tagami =================== 11.0
!    nspin_m  = nspin/(af+1)
!
    if ( noncol ) then
       nspin_m  = ndim_magmom
    else
       nspin_m  = nspin/(af+1)
    endif
! ====================================================================== 11.0
    allocate(  din_vhsr( nsize_rho_vhsr,nspin_m ) );  din_vhsr = 0.0d0
    allocate( dout_vhsr( nsize_rho_vhsr,nspin_m ) ); dout_vhsr = 0.0d0
    allocate(   dF_vhsr( nsize_rho_vhsr,nspin_m ) );   dF_vhsr = 0.0d0
    allocate( urec_vhsr( nsize_rho_vhsr,nspin_m,nbxmix,2) ); urec_vhsr = 0.0d0
    if(sw_gradient_simplex==ON)then
       allocate(rhoj_vhsr(nsize_rho_vhsr,nspin_m,nbxmix));rhoj_vhsr=0.d0
       allocate(Frhoj_vhsr(nsize_rho_vhsr,nspin_m,nbxmix));Frhoj_vhsr=0.d0
       allocate(rhoj_vhsro(nsize_rho_vhsr,nspin_m));rhoj_vhsro=0.d0
       allocate(Frhoj_vhsro(nsize_rho_vhsr,nspin_m));Frhoj_vhsro=0.d0
    endif
  end subroutine mix_pulay_allocate_vhsr

  subroutine mix_pulay_deallocate_vhsr
    if ( allocated( din_vhsr)) deallocate( din_vhsr)
    if ( allocated(dout_vhsr)) deallocate(dout_vhsr)
    if ( allocated(  dF_vhsr)) deallocate(  dF_vhsr)
    if ( allocated(urec_vhsr)) deallocate(urec_vhsr)
    if ( allocated(rhoj_vhsr)) deallocate(rhoj_vhsr)
    if ( allocated(Frhoj_vhsr)) deallocate(Frhoj_vhsr)
    if ( allocated(rhoj_vhsro)) deallocate(rhoj_vhsro)
    if ( allocated(Frhoj_vhsro)) deallocate(Frhoj_vhsro)
  end subroutine mix_pulay_deallocate_vhsr

  subroutine mix_pulay_alloc2_vhsr
    allocate( d0_vhsr( nsize_rho_vhsr,nspin_m) ); d0_vhsr = 0.0d0
  end subroutine mix_pulay_alloc2_vhsr

  subroutine mix_pulay_dealloc2_vhsr
    deallocate(d0_vhsr)
  end subroutine mix_pulay_dealloc2_vhsr

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
#ifdef __TIMER_SUB__
    call timer_sta(1113)
#endif
#ifdef _CDMIX_USE_POINTER_
    urec_l_3 => urec_l(ista_kgpm:iend_kgpm,1:kimg,1:nspin_m,j,iV)
    call mult1s(urec_l_3,dF_l,f_p,v_dF)!-(m_KE);<v|dF> ->v_dF
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

    call subtr_j_th_term(v_dF,iU,j,u_l)  !-(m_KE)

    !                        |u(m)> = |u(m)> - <v(j)|dF(m)>|u(j)>
! === DEBUG by tkato 2011/11/24 ================================================
!   if(hownew == RENEW) f(j,i,1:nspin) = v_dF(1:nspin)
    if(hownew == RENEW) f(j,i,1:nspin_m) = v_dF(1:nspin_m)
! ==============================================================================
#ifdef __TIMER_SUB__
    call timer_end(1113)
#endif
  end subroutine renew_u_br

  subroutine renew_d_br(j)
    integer, intent(in) :: j
! === DEBUG by tkato 2011/11/24 ================================================
!   real(DP)  :: vF(nspin)
    real(DP)  :: vF(nspin_m)
! ==============================================================================
#ifdef __TIMER_SUB__
    call timer_sta(1118)
#endif
#ifdef _CDMIX_USE_POINTER_
    urec_l_3 => urec_l(ista_kgpm:iend_kgpm,1:kimg,1:nspin_m,j,iV)
    call mult1s(urec_l_3,F_l,f_p,vF) !-(m_KE);<v|F>  ->vF
#else
    call mult1s5(urec_l,nbxmix,2,j,iV,F_l,f_p,vF) !-(m_KE);<v|F>  ->vF
#endif
! === DEBUG by tkato 2011/11/24 ================================================
!   if ( nspin==2 .and. sw_mix_bothspins_sametime == YES ) then
    if ( nspin_m==2 .and. sw_mix_bothspins_sametime == YES ) then
! ==============================================================================
      vF(1) = vF(1) + vF(2)
      vF(2) = vF(1)
    endif

! ================================= added by K. Tagami =============== 11.0
    if ( noncol ) then
       vF(1) = sum( vF(:) )
       vF(:) = vF(1)
    endif
! ===================================================================== 11.0

    call subtr_j_th_term(vF,iU,j,d0_l) !-(m_KE)
    !                        |d(m)> = |d(m)> - <v(j)|F(m)>|u(j)>
#ifdef __TIMER_SUB__
    call timer_end(1118)
#endif
  end subroutine renew_d_br

  subroutine renew_d_last_br(p)
    real(DP), intent(in), dimension(ista_kngp:iend_kngp) :: p
    integer   :: is, ik, i, ns
! === DEBUG by tkato 2011/11/24 ================================================
!   real(DP)  :: vF(nspin)
    real(DP)  :: vF(nspin_m)
! ==============================================================================
#ifdef __TIMER_SUB__
    call timer_sta(1120)
#endif

    call mult1s(v_l,F_l,f_p,vF)              !-(m_KE) <v|F> ->vF

! === DEBUG by tkato 2011/11/24 ================================================
!   if ( nspin==2 .and. sw_mix_bothspins_sametime == YES ) then
    if ( nspin_m==2 .and. sw_mix_bothspins_sametime == YES ) then
! ==============================================================================
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

#ifdef __TIMER_DO__
  call timer_sta(1166)
#endif
! ===================================== modified by K. Tagami ============== 11.0
!       do is = 1, nspin, af+1
       do is = 1, ndim_magmom, af+1
! ========================================================================== 11.0
          din (:,:,is) = kinqo_l(ista_kgpm:iend_kgpm,:,is)
          dout(:,:,is) = kinq_l (ista_kgpm:iend_kgpm,:,is)
       end do
#ifdef __TIMER_DO__
  call timer_end(1166)
#endif
    else
       call scatter_chg_onto_d(kinqo_l,din )  ! -(m_C.D.)
       call scatter_chg_onto_d(kinq_l, dout)  ! -(m_C.D.)
    end if

!!$    do is = 1, nspin, af+1
    ns = nspin_for_qnewton()
#ifdef __TIMER_DO__
  call timer_sta(1167)
#endif
    do is = 1, ns,af+1
       do ik = 1, kimg
          do i = ista_kgpm,iend_kgpm
             rho(i,ik,is) = d0_l(i,ik,is) - vF(is)*u_l(i,ik,is)
          end do
       end do
    end do
#ifdef __TIMER_DO__
  call timer_end(1167)
#endif

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
       call concentrate_d_to_chg(rho,kinq_l) !-(m_C.D.)
       call simple_mix_large_Gc(p)           !-(m_C.D.) kinq,kinqo,p ->kinq
    end if
#ifdef __TIMER_SUB__
    call timer_end(1120)
#endif
  end subroutine renew_d_last_br

! =========================== added by K. Tagami ================================== 5.0
  subroutine renew_u_br_with_vhsr(j,i)
    integer, intent(in) :: j,i

    integer       :: is
! === DEBUG by tkato 2011/11/24 ================================================
!   real(DP)      :: v_dF(nspin)
    real(DP)      :: v_dF(nspin_m)
! ==============================================================================

    v_dF = 0.d0

#ifdef _CDMIX_USE_POINTER_
    urec_l_3 => urec_l(ista_kgpm:iend_kgpm,1:kimg,1:nspin_m,j,iV)
    call mult1s(urec_l_3,dF_l,f_p,v_dF)!-(m_KE);<v|dF> ->v_dF
#else
    call mult1s5(urec_l,nbxmix,2,j,iV,dF_l,f_p,v_dF)
#endif

! =================================== modified by K. Tagami ============ 11.0
!       do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ======================================================================= 11.0
       v_dF(is) = v_dF(is) + sum( urec_vhsr(:,is,j,iV)*dF_vhsr(:,is) )
    End do
!
! === DEBUG by tkato 2011/11/24 ================================================
!   if ( nspin==2 .and. sw_mix_bothspins_sametime == YES ) then
    if ( nspin_m==2 .and. sw_mix_bothspins_sametime == YES ) then
! ==============================================================================
      v_dF(1) = v_dF(1) + v_dF(2)
      v_df(2) = v_dF(1)
    endif
!
! ======================== added by K. Tagami ==================== 11.0
    if ( noncol ) then
       v_dF(1) = sum( v_dF(:) )
       v_dF(:) = v_dF(1)
    endif
! ================================================================ 11.0

    call subtr_j_th_term(v_dF,iU,j,u_l)  !-(m_KE)
    !                        |u(m)> = |u(m)> - <v(j)|dF(m)>|u(j)>

! =================================== modified by K. Tagami ============ 11.0
!       do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ======================================================================= 11.0
       u_vhsr(:,is) = u_vhsr(:,is) - v_dF(is) *urec_vhsr(:,is,j,iU)
    End do
! === DEBUG by tkato 2011/11/24 ================================================
!   if(hownew == RENEW) f(j,i,1:nspin) = v_dF(1:nspin)
    if(hownew == RENEW) f(j,i,1:nspin_m) = v_dF(1:nspin_m)
! ==============================================================================

  end subroutine renew_u_br_with_vhsr

  subroutine renew_d_br_with_vhsr(j)
    integer, intent(in) :: j
! === DEBUG by tkato 2011/11/24 ================================================
!   real(DP)  :: vF(nspin)
    real(DP)  :: vF(nspin_m)
! ==============================================================================
    integer :: is

    vF = 0.d0

#ifdef _CDMIX_USE_POINTER_
    urec_l_3 => urec_l(ista_kgpm:iend_kgpm,1:kimg,1:nspin_m,j,iV)
    call mult1s(urec_l_3,F_l,f_p,vF) !-(m_KE);<v|F>  ->vF
#else
    call mult1s5(urec_l,nbxmix,2,j,iV,F_l,f_p,vF) !-(m_KE);<v|F>  ->vF
#endif

! =================================== modified by K. Tagami ============ 11.0
!       do is = 1, nspin, af+1
       do is = 1, ndim_magmom, af+1
! ======================================================================= 11.0
      vF(is) = vF(is) + sum( urec_vhsr(:,is,j,iV)*FF_vhsr(:,is) )
    End do

! === DEBUG by tkato 2011/11/24 ================================================
!   if ( nspin==2 .and. sw_mix_bothspins_sametime == YES ) then
    if ( nspin_m==2 .and. sw_mix_bothspins_sametime == YES ) then
! ==============================================================================
      vF(1) = vF(1) + vF(2)
      vF(2) = vF(1)
    endif

! ======================== added by K. Tagami ==================== 11.0
    if ( noncol ) then
       vF(1) = sum( vF(:) )
       vF(:) = vF(1)
    endif
! ================================================================ 11.0

    call subtr_j_th_term(vF,iU,j,d0_l) !-(m_KE)
    !                        |d(m)> = |d(m)> - <v(j)|F(m)>|u(j)>

! =================================== modified by K. Tagami ============ 11.0
!       do is = 1, nspin, af+1
       do is = 1, ndim_magmom, af+1
! ======================================================================= 11.0
       d0_vhsr(:,is) = d0_vhsr(:,is) - vF(is) *urec_vhsr(:,is,j,iU)
    end do
  end subroutine renew_d_br_with_vhsr

  subroutine renew_d_last_br_with_vhsr( p, rmxtrc_vhsr )
    real(DP), intent(in), dimension(ista_kngp:iend_kngp) :: p
    real(DP), intent(in) :: rmxtrc_vhsr(nspin_m)

    integer   :: is, ik, i, ns
! === DEBUG by tkato 2011/11/24 ================================================
!   real(DP)  :: vF(nspin)
    real(DP)  :: vF(nspin_m)
! ==============================================================================

    vF = 0.0d0
    call mult1s(v_l,F_l,f_p,vF)              !-(m_KE) <v|F> ->vF

! =================================== modified by K. Tagami ============ 11.0
!       do is = 1, nspin, af+1
       do is = 1, ndim_magmom, af+1
! ======================================================================= 11.0
       vF(is) = vF(is) + sum( v_vhsr(:,is)*FF_vhsr(:,is) )
    End do

! === DEBUG by tkato 2011/11/24 ================================================
!   if ( nspin==2 .and. sw_mix_bothspins_sametime == YES ) then
    if ( nspin_m==2 .and. sw_mix_bothspins_sametime == YES ) then
! ==============================================================================
      vF(1) = vF(1) + vF(2)
      vF(2) = vF(1)
    endif

! ========================== added by K. Tagami =============== 11.0
    if ( noncol ) then
       vF(1) = sum( vF(:) )
       vF(:) = vF(1)
    endif
! ============================================================= 11.0

    if(kgpm == kgp .or. npes == 1) then

! =================================== modified by K. Tagami ============ 11.0
!       do is = 1, nspin, af+1
       do is = 1, ndim_magmom, af+1
! ======================================================================= 11.0
          din (:,:,is) = kinqo_l(ista_kgpm:iend_kgpm,:,is)
          dout(:,:,is) = kinq_l (ista_kgpm:iend_kgpm,:,is)
       end do
    else
       call scatter_chg_onto_d(kinqo_l,din )  ! -(m_C.D.)
       call scatter_chg_onto_d(kinq_l, dout)  ! -(m_C.D.)
    end if

! =================================== modified by K. Tagami ============ 11.0
!     do is = 1, nspin, af+1
     do is = 1, ndim_magmom, af+1
! ======================================================================= 11.0
       din_vhsr (:,is) = rhoo_vhsr(:,is) ! kinqo
       dout_vhsr(:,is) = rho_vhsr (:,is) ! kinq
    end do

!!$    do is = 1, nspin, af+1
    ns = nspin_for_qnewton()
    do is = 1, ns,af+1
       do ik = 1, kimg
          do i = ista_kgpm,iend_kgpm
             rho(i,ik,is) = d0_l(i,ik,is) - vF(is)*u_l(i,ik,is)
          end do
       end do
    end do

    do is = 1, ns, af+1
       rho_vhsr(:,is) = d0_vhsr(:,is) - vF(is) *u_vhsr(:,is)
    end do

! ====================================== modified by K. Tagami ============ 11.0
!    if(sw_force_simple_mixing==ON .and. sw_recomposing==ON) call simple_mix2(c_p)
!
!    if ( sw_force_simple_mixing_hsr==ON .and. sw_recomposing_vhsr==ON ) then
!       call simple_mix2_kt( rmxtrc_vhsr )
!    endif

    if ( .not. noncol ) then
       if (sw_force_simple_mixing==ON .and. sw_recomposing==ON) then
          call simple_mix2(c_p)
       endif
       if ( sw_force_simple_mixing_hsr==ON .and. sw_recomposing_hsr==ON ) then
          call simple_mix2_kt( rmxtrc_vhsr )
       endif
    endif
! =========================================================================== 11.0

    if(kgpm < kgp) then
       call concentrate_d_to_chg(rho,kinq_l) !-(m_C.D.)
       call simple_mix_large_Gc(p)           !-(m_C.D.) kinq,kinqo,p ->kinq
    end if
  end subroutine renew_d_last_br_with_vhsr

! ============================================================================ 5.0

  subroutine simple_mix_large_Gc(p)
    real(kind=DP), intent(in), dimension(ista_kngp:iend_kngp) :: p

    integer :: istart, iend, is, ik, i
#ifdef __TIMER_SUB__
    call timer_sta(1123)
#endif
    istart = kgpm + 1; if(istart < ista_kngp) istart = ista_kngp
    iend   = kgp;      if(iend   > iend_kngp) iend   = iend_kngp
#ifdef __TIMER_DO__
  call timer_sta(1174)
#endif
    if(iend >= istart) then
! =================================== modified by K. Tagami ============ 11.0
!       do is = 1, nspin, af+1
       do is = 1, ndim_magmom, af+1
! ======================================================================= 11.0
          do ik = 1, kimg
             do i = istart, iend
                kinq_l(i,ik,is) = p(i)*kinq_l(i,ik,is) + (1.0d0-p(i))*kinqo_l(i,ik,is)
             end do
          end do
       end do
    end if
#ifdef __TIMER_DO__
  call timer_end(1174)
#endif
#ifdef __TIMER_SUB__
    call timer_end(1123)
#endif
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
  subroutine m_KE_mix_broyden2(nfout,rmx)
    integer, intent(in) :: nfout
    real(DP),intent(in) :: rmx

    integer   :: iter,j,mxiter,icr,jcr
!!$    real(DP)  :: v_dF(nspin),vF(nspin)
    integer   :: id_sname = -1
! --> T. Yamasaki  03 Aug. 2009
    real(kind=DP), allocatable, dimension(:):: rmxtrc ! d(nspin_m)
!   real(kind=DP), allocatable, dimension(:,:,:) :: kinqstore_l, kinqostore_l
! <--
#ifdef __TIMER_SUB__
    call timer_sta(1124)
#endif
    call tstatc0_begin('m_KE_mix_broyden2 ',id_sname,1)

    if(previous_waymix /= BROYD2.or.force_dealloc) then
       call mix_dealloc_previous()
       call mix_broyden_allocate();    F_l => din
       force_dealloc = .false.
    end if

    allocate(rmxtrc(nspin_m))
    if ( noncol ) then
       rmxtrc(1:nspin_m) = rmx
    else
       if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
          call alloc_kinqstore_recompos_kinq(rmx,rmxtrc) ! --> kinq_l, kinqo_l, rmxtrc
       else
          rmxtrc(1:nspin_m) = rmx
       end if
    endif
! ======================================================================== 11.0

! ====================== Modified by K. Tagami =========
!    allocate(c_p(ista_kngp:iend_kngp)); call precon_4_KE_mix(rmx,c_p)
! --> T. Yamasaki, 03rd Aug. 2009
!!$    allocate(c_p(ista_kngp:iend_kngp)); c_p = 0; call precon_4_KE_mix(rmx,c_p)
    allocate(c_p(ista_kngp:iend_kngp,nspin_m)); c_p = 0
! =======================================================


! ============================== modified by K. Tagami ================== 11.0
!    call precon_4_KE_mix(rmxtrc,c_p)
!
    if ( noncol ) then
       call precon_4_KE_mix_noncl(rmxtrc,c_p)
    else
       call precon_4_KE_mix(rmxtrc,c_p)
    endif
! ======================================================================== 11.0

    iter = iter_from_reset()                 !-(m_KE)

    if((iter-istrbr+1) <= 1) then
       call simple_mix1(c_p)                 !-(m_KE)
       !   din=kinqo_l; dout=kinq_l; (din,dout,c_p)->kinq_l
    else
!!$       stop ' -- iter-istrbr+1 > 1 (m_KE_mix_broyden2) --'
       call mix_broyden_alloc2   !-(m_KE) d0_l,u_l, and v_l are allocated
       call dF_F_d0_u_and_v      !-(c.h.)   dF_l, F_l, initial u_l,v_l,d0_l

       call set_ncrspd_mxiter_etc(iter,iU,mxiter) !-(m_KE) ->mxiter,ncrspd
       !                  when hownew == RENEW: f,g,ncrspd, and urec_l are reset.
       icr = icrspd_is(iter)                 !-(m_KE) function
       do j = 2, mxiter
          jcr = ncrspd(j)
          call renew_u_br(jcr,icr) !-(m_KE) |u(m)> = |u(m)> - <v(j)|dF(m)>|u(j)>
          call renew_d_br(jcr)     !-(m_KE) |d(m)> = |d(m)> - <v(j)|F(m)> |u(j)>
       enddo!j-loop

       urec_l(:,:,:,icr,iU) = u_l(:,:,:)  ! storing
       urec_l(:,:,:,icr,iV) = v_l(:,:,:)  ! storing

       call renew_d_last_br(c_p)       !-(m_KE) kinq_l(|d(m)>) = |d(m)>-<v(m)|F(m)>|u(m)>
       call mix_broyden_dealloc2                      !-(m_KE)
    endif

    if ( .not. noncol ) then
       if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
          call compos_kinq_dealloc_kinqstore()
       end if
    endif
    deallocate(rmxtrc)

    deallocate(c_p)
    previous_waymix = BROYD2
    call tstatc0_end(id_sname)
#ifdef __TIMER_SUB__
    call timer_end(1124)
#endif
  contains
    subroutine dF_F_d0_u_and_v
      !   dF_l(=deltaF) = (rho - dout) - (rhoo - din)
      !   F_l = rho - rhoo (=\cal F^{m}); u_l  = (rhoo - din) + c_p*dF_l;
      !   d0_l = rhoo+c_p* F_l;              v_l = dF_l/( |dF_l| )

      integer                      :: is,k,i
      real(DP), dimension(nspin_m) :: fff
#ifdef __TIMER_SUB__
    call timer_sta(1125)
#endif
#ifdef __TIMER_DO__
  call timer_sta(1175)
#endif
! ======================================= modified by K. Tagami =========== 11.0
!      do is = 1, nspin, af+1
      do is = 1, ndim_magmom, af+1
! ========================================================================= 11.0

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
#ifdef __TIMER_DO__
  call timer_end(1175)
#endif

      call mult1s(dF_l,dF_l,f_p,fff)
      if(sum(fff) < 1.d-40)  call phase_error_with_msg(nfout,' fmult is too small',__LINE__,__FILE__)

! === DEBUG by tkato 2011/11/24 ================================================
!     if ( nspin == 2 .and. sw_mix_bothspins_sametime == YES ) then
      if ( nspin_m == 2 .and. sw_mix_bothspins_sametime == YES ) then
! ==============================================================================
        fff(1) = fff(1) + fff(2)
        fff(2) = fff(1)
      endif

! ========================= added by K. Tagami =========================== 11.0
      if ( noncol ) then
         fff(1) = sum( fff(:) )
         fff(:) = fff(1)
      endif
! ======================================================================== 11.0

#ifdef __TIMER_DO__
  call timer_sta(1176)
#endif
! ========================================= modified by K. Tagami ========== 11.0
!      do is = 1, nspin, af+1
      do is = 1, ndim_magmom, af+1
! =========================================================================== 11.0
         v_l(:,:,is) = dF_l(:,:,is)/fff(is)
      end do
#ifdef __TIMER_DO__
  call timer_end(1176)
#endif

#ifdef __TIMER_SUB__
    call timer_end(1125)
#endif
    end subroutine dF_F_d0_u_and_v

  end subroutine m_KE_mix_broyden2

  subroutine mix_pulay_allocate
!!$    if(allocated(f_p)) return

! =============================== modified by K. Tagami ========== 11.0
!    nspin_m  = nspin/(af+1)
!
    if ( noncol ) then
       nspin_m  = ndim_magmom
    else
       nspin_m  = nspin/(af+1)
    endif
! ================================================================= 11.0

! =========================================== Modified by K. Tagami =========
!    allocate(f_p(ista_kgpm:iend_kgpm)); call precon_4_mult(f_p) !-(m_KE)
    allocate(f_p(ista_kgpm:iend_kgpm)); f_p = 0; call precon_4_mult(f_p) !-(m_KE)
! ============================================================================
    f_p = 1.0d0

    allocate(din(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(dout(ista_kgpm:iend_kgpm,kimg,nspin_m))
    allocate(urec_l(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix,2))
    allocate(uuf_p(nbxmix,nspin_m))
    allocate(f(nbxmix,nbxmix,nspin_m))
    allocate(g_p(nbxmix,nspin_m))
    allocate(prj_wk(mp_kngp*kimg*nspin_m))
    allocate(ncrspd(nbxmix))

    if(sw_gradient_simplex==ON)then
       allocate(rhoj(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix));rhoj=0.d0
       allocate(Frhoj(ista_kgpm:iend_kgpm,kimg,nspin_m,nbxmix));Frhoj=0.d0
       allocate(rhojo(ista_kgpm:iend_kgpm,kimg,nspin_m));rhojo=0.d0
       allocate(Frhojo(ista_kgpm:iend_kgpm,kimg,nspin_m));Frhojo=0.d0
    endif

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
  end subroutine mix_pulay_deallocate

  subroutine mix_pulay_alloc2
    allocate(d0_l(ista_kgpm:iend_kgpm,kimg,nspin_m))
! =========================================== Added by K. Tagami ========
    d0_l = 0.0d0
! =======================================================================
    call alloc_rho_rhoo_and_cpm
  end subroutine mix_pulay_alloc2

  subroutine mix_pulay_dealloc2
    deallocate(d0_l)
    call dealloc_rho_rhoo_and_cpm
  end subroutine mix_pulay_dealloc2

  subroutine m_KE_mix_pulay(nfout,rmx)
    integer, parameter  :: iRho = 1, iResid = 2
    integer, intent(in) :: nfout
    real(DP),intent(in) :: rmx
    integer   :: iter, mxiter
    real(DP),pointer,dimension(:)  :: e_wk, f_wk, ww1, finv
    integer, pointer,dimension(:)  :: ip
! --> T. Yamasaki  03 Aug. 2009
    real(kind=DP), allocatable, dimension(:):: rmxtrc ! d(nspin_m)
!   real(kind=DP), allocatable, dimension(:,:,:) :: kinqstore_l, kinqostore_l
! <--
    integer   :: id_sname = -1
#ifdef __TIMER_SUB__
    call timer_sta(1131)
#endif
    call tstatc0_begin('m_KE_mix_pulay ',id_sname,1)

    if(previous_waymix /= PULAY.or.force_dealloc) then
       force_dealloc = .false.
       call mix_dealloc_previous()
       call mix_pulay_allocate()
    end if

    allocate(rmxtrc(nspin_m))
    if ( noncol ) then
       rmxtrc(1:nspin_m) = rmx
    else
       if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) then
          call alloc_kinqstore_recompos_kinq(rmx,rmxtrc) ! --> kinq_l,kinqo_l, rmxtrc
       else
          rmxtrc(1:nspin_m) = rmx
       end if
    endif
! ========================================================================= 11.0

! ====================================== Modified by K. Tagami =========
!    allocate(c_p(ista_kngp:iend_kngp)); call precon_4_KE_mix(rmx,c_p)
! --> T. Yamasaki, 03rd Aug. 2009
!!$    allocate(c_p(ista_kngp:iend_kngp)); c_p = 0; call precon_4_KE_mix(rmx,c_p)
    allocate(c_p(ista_kngp:iend_kngp,nspin_m)); c_p = 0
! ========================================================================


! =================================== modified by K. Tagami =============== 11.0
!    call precon_4_KE_mix(rmxtrc,c_p)
!
    if ( noncol ) then
       call precon_4_KE_mix_noncl(rmxtrc,c_p)
    else
       call precon_4_KE_mix(rmxtrc,c_p)
    endif
! ========================================================================= 11.0

    iter = iter_from_reset()                 !-(m_KE)
    if((iter-istrbr+1) <= 1) then
       call simple_mix1(c_p)                 !-(m_KE)
       !   din=kinqo_l; dout=kinq_l; (din,dout,c_p)->kinq_l
    else
       call mix_pulay_alloc2   !-(m_KE) d0_l,u_l, and w_l are allocated
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
       call mix_pulay_dealloc2                    !-(m_KE)
    endif

    if ( .not. noncol ) then
       if(sw_recomposing == YES .and. af == 0 .and. nspin == 2) &
            & call compos_kinq_dealloc_kinqstore()
    endif
    deallocate(rmxtrc)

    deallocate(c_p)
    previous_waymix = PULAY
    call tstatc0_end(id_sname)
#ifdef __TIMER_SUB__
    call timer_end(1131)
#endif
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
#ifdef __TIMER_SUB__
    call timer_sta(1132)
#endif
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
#ifdef __TIMER_SUB__
    call timer_end(1132)
#endif
    end subroutine Resid_and_dd_into_urec

    subroutine Ri_dot_Rj(n)
      integer, intent(in) :: n
      integer  :: it,jt,itc,jtc
! === DEBUG by tkato 2011/11/24 ================================================
!     real(DP) :: ff1(nspin),ff1tmp
      real(DP) :: ff1(nspin_m),ff1tmp
! ==============================================================================
#ifdef __TIMER_SUB__
    call timer_sta(1133)
#endif
#ifdef __TIMER_DO__
  call timer_sta(1179)
#endif
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
#ifdef __TIMER_DO__
  call timer_end(1179)
#endif
#ifdef __TIMER_SUB__
    call timer_end(1133)
#endif
    end subroutine Ri_dot_Rj

    subroutine Rj_dot_d(n)
      integer, intent(in) :: n
      integer  :: jt, jtc
! === DEBUG by tkato 2011/11/24 ================================================
!     real(DP) :: ff1(nspin),ff1tmp
      real(DP) :: ff1(nspin_m),ff1tmp
! ==============================================================================
#ifdef __TIMER_SUB__
    call timer_sta(1138)
#endif
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
#ifdef __TIMER_SUB__
    call timer_end(1138)
#endif
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
               fwork(jt,it) = f(jt,it,is)
               if(it==jt) fwork(jt,it)=fwork(jt,it)+alpha_pulay
            enddo
         enddo
         call dpotrf('U',n,fwork,n,inf)
         call dpotri('U',n,fwork,n,inf)
         if(inf/=0)then
           if(printable) write(nfout,*) ' !** failed calculation of the inverse matrix for the pulay mixer'
           fwork=0
           do it=1,n
              do jt=1,n
                 fwork(jt,it) = f(jt,it,is)
                 if(it==jt) fwork(jt,it)=fwork(jt,it)+1.d-7
              enddo
           enddo
           call dpotrf('U',n,fwork,n,inf)
           call dpotri('U',n,fwork,n,inf)
         endif
         do it=1,n-1
            do jt=it+1,n
               fwork(jt,it) = fwork(it,jt)
            enddo
         enddo
         do it=1,n
            do jt=1,n
               f(jt,it,is) = fwork(jt,it)
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
#ifdef __TIMER_SUB__
    call timer_sta(1134)
#endif

      e_wk = 0.d0
      do it = 1, n
         e_wk(it*it) = 1.d0
      end do

! ======================================= modified by K. Tagami =========== 11.0
!      do is = 1, nspin, af+1
      do is = 1, ndim_magmom, af+1
! ========================================================================== 11.0
         div = 1.d0/f(1,1,is)
         icount = 1
#ifdef __TIMER_DO__
  call timer_sta(1180)
#endif
         do jt = 1, n
            do it = 1, n
               f_wk(icount) = f(it,jt,is)*div
               icount = icount + 1
            end do
         end do
#ifdef __TIMER_DO__
  call timer_end(1180)
#endif
         if(ipripulay >= 1) then
            write(nfout,600) n,(('(',it,jt,')',f(it,jt,is)*div,jt=1,n),it=1,n)
600         format(//11x,"**input matrix**"/12x &
                 & ,"horder=",I5/(2x,4(1x,1a,i2,",",i2,1a,e14.6)))
         end if
         call rdecomp(n,f_wk,ww1,ip,icon)
         if(icon /= 0) then
            call phase_error_with_msg(nfout,'LU decomposition is impossible.',__LINE__,__FILE__)
         else
            call rsolve(n,n,f_wk,e_wk,finv,ip)
         endif

         icount = 1
#ifdef __TIMER_DO__
  call timer_sta(1181)
#endif
         do jt = 1, n
            do it = 1, n
               f(it,jt,is) = finv(icount)
               icount = icount + 1
            end do
         end do
#ifdef __TIMER_DO__
  call timer_end(1181)
#endif

         if(ipripulay >= 1) then
            write(nfout,630) (('(',it,jt,')',f(it,jt,is),it=1,n),jt=1,n)
630         format(/11x, "**inverse matrix**" &
                 & ,/(2x,4(1x,1a,i2,",",i2,1a,e14.6)))
         end if
      end do
#ifdef __TIMER_SUB__
    call timer_end(1134)
#endif
    end subroutine get_finv

    subroutine get_gmatrix(n)
      integer,intent(in) :: n
      integer :: is, it, jt, nnspin
#ifdef __TIMER_SUB__
    call timer_sta(1139)
#endif
      nnspin = nspin
      if(sw_mix_bothspins_sametime==ON .or. af==1) nnspin=1

! ============================ added by K. Tagami ============= 11.0
      if ( noncol ) nnspin = 1
! ============================================================== 11.0

      g_p = 0.d0
      do is = 1, nnspin
#ifdef __TIMER_DO__
  call timer_sta(1188)
#endif
         do it = 1, n
            do jt = 1, n
               g_p(it,is) = g_p(it,is) - f(jt,it,is)*uuf_p(jt,is)
            end do
         end do
#ifdef __TIMER_DO__
  call timer_end(1188)
#endif
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

#ifdef __TIMER_SUB__
    call timer_end(1139)
#endif
    end subroutine get_gmatrix

    subroutine renew_d_using_g(n,p)
      integer, intent(in)                                :: n
      real(DP),intent(in),dimension(ista_kgpm:iend_kgpm,nspin_m) :: p
      integer    :: is, k, i, it, itc, ns
#ifdef __TIMER_SUB__
    call timer_sta(1140)
#endif

!!$      do is = 1, nspin, af+1
      ns = nspin_for_qnewton()
      do is = 1, ns,af+1
         do k = 1, kimg
#ifdef __TIMER_DO__
  call timer_sta(1189)
#endif
            do i = ista_kngp, iend_kngp
               rho(i,k,is)  = rhoo(i,k,is) + p(i,is)*d0_l(i,k,is)
            end do
#ifdef __TIMER_DO__
  call timer_end(1189)
#endif
#ifdef __TIMER_DO__
  call timer_sta(1190)
#endif
            do it = 1, n
               itc = ncrspd(it)
               do i = ista_kngp, iend_kngp
                  rho(i,k,is) = rho(i,k,is) + g_p(it,is)* &
                       &        (urec_l(i,k,is,itc,iRho) + p(i,is)*urec_l(i,k,is,itc,iResid))
               end do
            end do
#ifdef __TIMER_DO__
  call timer_end(1190)
#endif
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
         call concentrate_d_to_chg(rho,kinq_l) !-(m_C.D.)
         call simple_mix_large_Gc(c_p)         !-(m_C.D.) kinq,kinqo,c_p ->kinq
      end if
#ifdef __TIMER_SUB__
    call timer_end(1140)
#endif
    end subroutine renew_d_using_g

  end subroutine m_KE_mix_pulay

  subroutine m_KE_simple_mixing_00( rmxt )
    real(kind=DP), intent(in) :: rmxt

    if ( use_symm_ekin_density ) then
       ekins_l = ekins_old + rmxt*( ekins_l -ekins_old )
    endif
    if ( use_asymm_ekin_density ) then
       ekina_l = ekina_old + rmxt*( ekina_l -ekina_old )
    endif
  end subroutine m_KE_simple_mixing_00

  subroutine m_KE_kerker_mixing( rmxt )
    real(kind=DP), intent(in) :: rmxt

    integer :: i, is, ri
    real(kind=DP) :: gg, q0, fac_min
    real(kind=DP) :: rmxtrc(nspin)
!
    real(kind=DP), allocatable :: c_p_ekin(:,:), dekin(:,:,:)

    allocate( c_p_ekin( ista_kngp:iend_kngp,nspin) );  c_p_ekin = 0.0d0
    allocate( dekin( ista_kngp:iend_kngp,kimg,nspin) );     dekin = 0.0d0

    q0 = 1.0d0
!    q0 = 5.0d0
!    q0 = 2.0d0

    q0 = g0_wf_precon **2
    fac_min = amin_wf_precon

    Do i=ista_kngp, iend_kngp
       gg = gr_l(i)*gr_l(i)
       Do is=1, nspin
          c_p_ekin(i,is) = max( gg/(gg+q0), fac_min )
!          write(1000+mype,*) i, is, c_p_ekin(i,is)
       End do
    End Do
!    stop

    c_p_ekin = c_p_ekin *rmxt

    if ( use_symm_ekin_density ) then
       dekin = ekins_l - ekins_old
       Do is=1, nspin
          Do ri=1, kimg
             Do i=ista_kngp, iend_kngp
                ekins_l(i,ri,is) = ekins_old(i,ri,is) &
                     &  +c_p_ekin(i,is) *dekin(i,ri,is)
             End do
          End do
       End do
    endif
    if ( use_asymm_ekin_density ) then
       dekin = ekina_l - ekina_old
       Do is=1, nspin
          Do ri=1, kimg
             Do i=ista_kngp, iend_kngp
                ekina_l(i,ri,is) = ekina_old(i,ri,is) &
                     &  +c_p_ekin(i,is) *dekin(i,ri,is)
             End do
          End do
       End do
    endif

    deallocate( c_p_ekin)
    deallocate( dekin )

  end subroutine m_KE_kerker_mixing

end module m_KE_mixing
