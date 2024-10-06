!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  MODULE: m_XC_Potential
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!  FURTHER MODIFICATION: T. Yamasaki, January/13/2004, August/23/2006
!                                     September/19/2006
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
! =========== Contributions ===================================
!
! Through the courtesy of contributors, the following functions are added.
!
! Company:  ASMS Co.,Ltd.
! Functions:  [Identifier: 13.0XX]
!                 modified BJ meta-gga functional is implemented (experimental)
!
! =============================================================
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

#define XC_PACK_FFT
module m_XC_Potential
! $Id: m_XC_Potential.F90 633 2020-12-01 05:11:03Z jkoga $
!
!  Upgraded on 23rd Aug. 2006 by T. Yamasaki
!    Differentials of the charge density function in GGA calculation are
!    parallelized according to the components of x, y, and z.
!  Upgraded on 19th Sep. 2006 by T. Yamasaki
!    FFTs in GGA calculations are parallelized.
!
  use m_PlaneWaveBasisSet,    only : ngabc,gr_l,igfp_l,kg,kgp,kgp_reduced, ylm_l&
       &                           , m_pwBS_sphrp2, igfp_nonpara
#ifdef _MPIFFT_
  use m_PlaneWaveBasisSet,    only : igfp_l_c
#endif
  use m_PseudoPotential,      only : itpcc,ilmt,ltp,taup,il2p,isph,iqitg,rhpcg_l,qitg_l&
       &                           , rhpcg_diff_l,qitg_diff_l,dl2p &
       &                           , m_PP_include_vanderbilt_pot&
       &                           , m_PP_set_index_arrays1, m_PP_set_index_arrays2 &
       &                           , m_PP_find_maximum_l,nlmt, rhcg_l, nqitg &
       &                           , m_PP_set_index_arrays1 &
       &                           , m_PP_set_index_arrays2
  use m_Crystal_Structure,    only : rltv,univol
  use m_Ionic_System,         only : ntyp,natm,iwei,ityp,pos,zfm3_l
  use m_FFT,                  only : fft_box_size_CD, fft_box_size_CD_c, nfftp &
       &                           , m_FFT_check_of_negative_CD, fft_box_size_CD_nonpara, nfftp_nonpara
#ifdef _MPIFFT_
  use m_FFT,                  only : m_FFT_set_cdata  &
       &                           , lx,ly,lz,ly_d,lz_d, ny_d,nz_d
#else
  use m_FFT,                  only : m_FFT_alloc_CD_box &
       &                           , m_FFT_dealloc_CD_box
#endif
  use m_Timing,               only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters,   only : xctype,len_xctype,ipri,iprixc,nspin,kimg,af,nel_Ylm     &
                                   , sw_hybrid_functional,sw_exchange_only,alpha_exx &
       &                           , sw_screened_exchange,omega_exx_pbe &
       &                           , hybrid_functional_type &
       &                           , iprinegativecharge, iprinegativecharge0 &
       &                           , sw_output_xc_seperately, exchange_pot_type &
       &                           , vdwdf_version, m_CtrlP_cachesize &
       &                           , sw_communicator_for_chg &
       &                           , m_CtrlP_cachesize &
#ifdef ENABLE_ESM_PACK
       &                           , max_warnings_negativecharge,sw_esm               &
#else
       &                           , max_warnings_negativecharge                      &
#endif
       &                           , sw_fft_xzy
  use m_Const_Parameters,     only : DP,DOWN,UP,Partial_Core_Charge,Valence_plus_PC_Charge &
       &                           , LDA,GGA,PAI,PAI2,VXC_AND_EXC,STRESS_,zi,SKIP, OFF, ON
  use m_Parallelization,      only : mp_kngp,ista_kngp,iend_kngp,npes,mype, ierr &
       &                           , npes_cdfft,nrank_ggacmp, myrank_cdfft, myrank_ggacmp &
       &                           , ista_fftp, iend_fftp, ista_fftph, iend_fftph &
       &                           , nis_fftp, nie_fftp, nel_fftp, idisp_fftp,np_fftp,mp_fftp &
       &                           , nel_fftph, idisp_fftph, nrest_cdfft &
       &                           , mpi_ggacmp_cross_world,mpi_cdfft_world, map_ggacmp &
       &                           , map_pe2ggacmp, map_pe2cdfft, MPI_CommGroup &
       &                           , myrank_g, myrank_e, nrank_g, nrank_e, ista_atm_ke, iend_atm_ke
!!$       &                           , ista_sfftp,iend_sfftp,ista_sfftph,iend_sfftph &
!!$       &                           , nis_sfftp,nie_sfftp,nel_sfftp,idisp_sfftp,np_sfftp,mp_sfftp &
!!$       &                 , npes_cdfft,nrank_ggacmp, myrank_cdfft, myrank_ggacmp, max_ggacmp &
  use m_FFT,                  only : fft_box_size_CD_3D
  use m_PlaneWaveBasisSet,   only : ngabc_kngp_l, ngabc_kngp_B_l

! ================================ added by K. Tagami ================ 11.0
  use m_Control_Parameters,    only : noncol, ndim_chgpot, ndim_magmom
  use m_Const_Parameters,      only : Valence_Charge_Only, Valence_plus_Core_Charge, &
       &                              EXC_ONLY
  use m_ES_NonCollinear,       only : m_ES_cp_VxcR_to_afft, &
       &                              m_ES_DensMat_To_MagMom_Gspace, &
       &                              m_ES_get_Angles_MagMom_Rspace, &
       &                              m_ES_SpinDens_Along_QuantzAxis, &
       &                              m_ES_DensMat_to_MagMom_vlhxcl, &
       &                              m_ES_get_Angles_MagMom_Rspace2, &
       &                              m_ES_cp_VxcR_to_afft2, &
       &                              m_ES_SpinDens_Along_QuantzAxis2,&
       &                              m_ES_Chgsoft_Along_QuantzAxis, &
       &                          m_ES_QuantzAxis_inv_flg_atm, &
       &                          m_ES_Chghard_Along_QuantzAxis, &
       &                          m_ES_GradCorr_Along_QuantzAxis2
  use m_Crystal_Structure,  only : sw_constrain_on_grad_correction
! ==================================================================== 11.0

! ================================ added by K. Tagami ================ 11.0
!!!  use m_PseudoPotential,  only : flg_paw
!!  use m_Charge_Density,   only : magmom_local_now
! ==================================================================== 11.0

! ===== KT_add ====== 13.0XX, 2014/09/19
  use m_Control_Parameters,     only : use_metagga, ekin_density_is_active,oneshot, &
       &                               use_modeled_ekin_density, &
       &                               use_asymm_ekin_density, use_symm_ekin_density, &
       &                               sw_rspace_ekin_density, &
       &                               sw_calc_ekin_density_hardpart, &
       &                               sw_add_ekin_hardpart_on_Gspace, &
       &                               val_c_tb09, sw_fix_val_c_tb09, val_g_tb09_paw, &
       &                               ekin_density_type, icond
  use m_KineticEnergy_Density,  only : ekins_l, ekina_l, &
       &                               m_KE_set_modeled_ekin_density
! =================== 13.0XX. 2014/09/19

! === KT_add === 2014/08/04
  use m_Crystal_Structure,  only : sw_neglect_low_vorticity, sw_neglect_low_helicity, &
       &                           sw_monitor_magnetic_vorticity
  use m_ES_Noncollinear,   only : m_ES_add_contrib_to_Vorticity, &
       &                          m_ES_neglect_low_Helicity, &
       &                          m_ES_neglect_low_Vorticity
! ============== 2014/08/04

  use m_XC_Variables,      only : vxc_l, vxco_l, vxcpc_l, exc, excpc, eex, ecor
  use m_XC_Potential_2D,   only : m_XC_cal_potential

  use m_Control_Parameters,  only : sw_opencore, sw_xc_opencore_ae_only
  use m_PS_opencore,  only : has_opencore, rmag_opencore_l, mag_opencore_pol

  use m_IterationNumbers,   only : iteration
  use m_Control_Parameters,  only : vtau_exists
  use m_Electronic_Structure,  only : vtau_l
#ifdef LIBXC
  use m_Control_Parameters,  only : xc_family_exch, xc_family_corr, xc_name_exch
  use xc_f03_lib_m
#endif
  use mpi

!!#ifdef __EDA__
!!  use m_Control_Parameters, only : sw_eda
!!#endif
  implicit none

!  include 'mpif.h'
  integer istatus(mpi_status_size)

#ifdef __EDA__
! -----  ascat starts modifying  -----
!!  real(kind=DP),public,target,allocatable,dimension(:) :: exc_on_a_grid
  real(kind=DP),private,allocatable,dimension(:) :: exc_on_a_grid_wk
! -----  ascat ceases modifying  -----
#endif

  real(kind=DP),private,allocatable,dimension(:):: afft        ! d(ista_fftp:iend_fftp)
  real(kind=DP),private,pointer,dimension(:,:)  :: chgrhr_l    ! MPI d(ista_fftph:iend_fftph,ispin)
! GGA arrays
  integer, private,pointer,dimension(:)         :: inx,jnx,knx ! MPI d(ista_fftp:iend_fftp)
  real(kind=DP),private,pointer,dimension(:,:)  :: chden_l         ! MPI d(ista_fftp :iend_fftp)
  real(kind=DP),private,pointer,dimension(:)    :: grad_trho       ! MPI d(ista_fftph:iend_fftph)
  real(kind=DP),private,allocatable,dimension(:,:)  :: grad_rho    ! MPI d(ista_fftph:iend_fftph, nspin)
#ifndef _XC_SAVE_MEMORY_
  real(kind=DP),private,allocatable,dimension(:,:,:):: cgrad_rho ! MPI d(ista_fftph:iend_fftph,3,nspin)
#else
  real(kind=DP),private,allocatable,dimension(:):: cggawk13 ! MPI d(ista_fftph:iend_fftph)
#endif
  real(kind=DP),private,pointer,dimension(:,:)  :: dF_drho,dF_dgradrho ! MPI d(ista_fftph:iend_fftph,nspin)

! === meta-gga
  real(kind=DP),private,allocatable,dimension(:,:)  :: dF_dlaplrho, dF_dtau

  real(kind=DP) :: val_g
  real(kind=DP),private,allocatable,dimension(:,:) :: grad2_rho, ekin_dens
  real(kind=DP),private,allocatable,dimension(:,:) :: rho_with_core, grad_rho_with_core
                                                 ! MPI d(ista_fftph:iend_fftph,nspin)


  real(kind=DP),        allocatable, dimension(:,:) :: s_gga1,s_gga2      ! d(3,3)

  real(kind=DP),private, pointer, dimension(:)      :: grinv,zfcos,zfsin,mzfsin  ! MPI d(ista_kngp:iend_kngp)
  real(kind=DP),private, pointer, dimension(:)      :: zfcos_blk,zfsin_blk,mzfsin_blk
  real(kind=DP),private, pointer, dimension(:)      :: flchgq             ! MPI d(nspin)
  real(kind=DP),private, pointer, dimension(:)      :: chgfft             ! MPI d(ista_fftp:iend_fftp)
  real(kind=DP),private, pointer, dimension(:,:)    :: alinvt,g
  real(kind=DP),private, allocatable, dimension(:,:,:) :: drhodh,flchgqd,ylmd
  real(kind=DP),private, pointer, dimension(:)      :: ylm
  real(kind=DP),private, pointer, dimension(:)      :: f2or1 ! MPI d(ista_fftph,iend_fftph)
  real(kind=DP),private, pointer, dimension(:)      :: chgfft_mpi            ! MPI
  logical,      private, dimension(3)               :: lmn_even

!!$  integer, private, allocatable, dimension(:,:)     :: np_send_map_charge2fftbox  ! d(0:npes-1,0:npes-1)
!!$  integer, private, allocatable, dimension(:)       :: ip_map_charge2fftbox ! d(ista_kngp:iend_kngp)
  logical, private                                  :: np_send_is_set = .false.

  integer, allocatable, dimension(:)  :: ip ! d(ista_kngp:iend_kngp,3)
  integer, allocatable, dimension(:,:):: igfp_ijk ! d(ista_kngp:iend_kngp+(0,1),3)
  integer, allocatable, dimension(:,:):: np_send ! d(0:npes-1,0:npes-1)
  real(kind=DP),private,allocatable,dimension(:,:):: afft_l   !
  logical, private                                  :: np_send_is_set_3D = .false.
  integer, private                             :: nfft_div_size

  real(kind=DP), private :: ecnl = 0.d0

  integer :: nblk = 1000

  real(kind=DP), allocatable, target, dimension(:,:) :: ylm_ext
! =================================== added by K. Tagami ============== 11.0
!  integer :: sw_constrain_on_grad_correction = OFF
!  integer :: sw_constrain_on_grad_correction = ON
! ===================================================================== 11.0

!   1. check_of_xctype (function)
!   2. m_XC_alloc_vxc
!   3. m_XC_dealloc_vxc
!   4. cp_vxc_to_vxco
!   5. xc_allocate
!       set_f2or1
!   6. xc_deallocate
!   7. m_XC_alloc_s_gga
!   8. m_XC_dealloc_s_gga
!   9. m_XC_cal_potential
!       ggaxcp, ggaxcp_diff, ggaxcp0, sum_s_gga12, stress_exchange_part
!      ,stress_correlation_part, rhos_diff, rhopc_diff, rhoh_diff, even_case
!      ,odd_case, real_case, complex_case, map_drhodh1, subst_A_into_A_reduced, map_drhodh2
!      ,dgrhodh1, dgrhodh2, gradrho1, gradrho2, allocate_for_stress, deallocate_for_stress
!      ,abs_grad_rho_up_down_total, dFxc_over_ddgradrho, finally_gga_xc_pot
!      ,add_negative_afft_to_grad_rho,G_xyz_afft, dFxc_dgradrho_dot_gradrho
!      ,dFxc_dgradrho_dot_gradrho2, x_dot_cggawk13_into_afft, cp_afft_to_cggawk13
!      ,g_xyz_total_chden_l,add_sq_afft_to_grad_trho,add_sq_afft_to_grad_rho
!      ,g_xyz_chden_l,check_lmn_even
!      ,boundary_zero_into_afft, initialize_cggawk_arrays
!      ,cp_afft_to_cgrad_rho, cp_afft_to_chgrhr, wd_small_part_of_vxc
!      ,wd_small_part_of_afft, map_charge_onto_a_fft_box, set_ispin,cpafft_to_vxc_or_vxcpc
!      ,cpafft
!  10. xcpotf          <-- (9)
!  11. check_of_pcc
!  12. afft_allgatherv

interface
#ifdef __EDA__
   subroutine ex_ggapw91(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc &
        &                    ,dFx_drho,dFx_dgradrho, exc_on_a_grid_wk,ist,ien)
#else
   subroutine ex_ggapw91(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc &
        &                    ,dFx_drho,dFx_dgradrho,ist,ien)
#endif
  !!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
     use m_Const_Parameters,  only : DP
     implicit none

     integer,intent(in)        :: nspin,ispin,ista_r,iend_r
     integer,intent(in)        :: ist, ien
     real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
     real(kind=DP),intent(in)  :: grad_rho(ista_r:iend_r,nspin)
     real(kind=DP),intent(in)  :: wos(ista_r:iend_r)
     real(kind=DP),intent(out) :: exc
     real(kind=DP),intent(out) :: dFx_drho(ista_r:iend_r,nspin)
     real(kind=DP),intent(out) :: dFx_dgradrho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
     real(kind=DP),intent(out) :: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif
   end subroutine ex_ggapw91

#ifdef __EDA__
  subroutine cr_ggapw91(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_trho,wos,exc,dF_drho,exc_on_a_grid_wk,ist,ien)
#else
  subroutine cr_ggapw91(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_trho,wos,exc,dF_drho,ist,ien)
#endif
!!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : PAI,DP
  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in)        :: ist, ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(inout) :: grad_trho(ista_r:iend_r)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)
  real(kind=DP),intent(inout) :: exc
  real(kind=DP),intent(inout) :: dF_drho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif
  end subroutine cr_ggapw91

#ifndef PREV_EX_GGAPBE
#ifdef __EDA__
  subroutine ex_ggapbe(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dFx_drho,dFx_dgradrho,revPBE,exc_on_a_grid_wk, ist,ien)
#else
  subroutine ex_ggapbe(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dFx_drho,dFx_dgradrho,revPBE,ist,ien)
#endif
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,PAI
  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in)        :: ist, ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)  :: grad_rho(ista_r:iend_r,nspin)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)
  real(kind=DP),intent(out) :: exc
  real(kind=DP),intent(out) :: dFx_drho(ista_r:iend_r,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif
  logical, intent(in),optional :: revPBE
  end subroutine ex_ggapbe
#else
#ifdef __EDA__
  subroutine ex_ggapbe(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dFx_drho,dFx_dgradrho,exc_on_a_grid_wk, ien)
#else
  subroutine ex_ggapbe(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dFx_drho,dFx_dgradrho,ien)
#endif
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,PAI
  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in),optional :: ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)  :: grad_rho(ista_r:iend_r,nspin)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)
  real(kind=DP),intent(out) :: exc
  real(kind=DP),intent(out) :: dFx_drho(ista_r:iend_r,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif
  end subroutine ex_ggapbe
#endif
#ifndef PREV_CR_GGAPBE
#ifdef __EDA__
  subroutine cr_ggapbe(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_trho,wos,exc,dF_drho,ecor,exc_on_a_grid_wk, ist,ien)
#else
  subroutine cr_ggapbe(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_trho,wos,exc,dF_drho,ecor,ien)
#endif
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : PAI,DP
  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
#ifdef __EDA__
  integer,intent(in)        :: ist
#endif
  integer,intent(in)        :: ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(inout) :: grad_trho(ista_r:iend_r)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)
  real(kind=DP),intent(inout) :: exc
  real(kind=DP),intent(inout) :: dF_drho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif
  real(kind=DP),intent(out),optional :: ecor
  end subroutine cr_ggapbe
#else
#ifdef __EDA__
  subroutine cr_ggapbe(nspin,ispin,chgrhr_l,grad_trho,f2or1,exc,dF_drho,exc_on_a_grid_wk)
#else
  subroutine cr_ggapbe(nspin,ispin,chgrhr_l,grad_trho,f2or1,exc,dF_drho)
#endif
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : PAI,DP
  use m_Parallelization,   only : ista_fftph, iend_fftph
  implicit none

  integer,intent(in)        :: nspin,ispin
  real(kind=DP),intent(in)  :: chgrhr_l(ista_fftph:iend_fftph,ispin)
  real(kind=DP),intent(inout) :: grad_trho(ista_fftph:iend_fftph)
  real(kind=DP),intent(in)  :: f2or1(ista_fftph:iend_fftph)
  real(kind=DP),intent(inout) :: exc
  real(kind=DP),intent(inout) :: dF_drho(ista_fftph:iend_fftph,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_fftph:iend_fftph)
! -----  ascat ceases modifying  -----
#endif
  end subroutine cr_ggapbe
#endif
#ifdef __EDA__
  subroutine xclda(nspin,ispin,ista_r,iend_r,chgrhr_l,wos,exc,dF_drho,exc_on_a_grid_wk, ien)
#else
  subroutine xclda(nspin,ispin,ista_r,iend_r,chgrhr_l,wos,exc,dF_drho,ien)
#endif
  use m_Const_Parameters,  only : DP, PAI
  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in),optional :: ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)
  real(kind=DP),intent(out) :: exc
!!$  real(kind=DP),intent(inout) :: dF_drho(ista_r:iend_r,nspin)
  real(kind=DP),intent(out) :: dF_drho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif
  end subroutine xclda

#ifdef __EDA__
  subroutine ggabek(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dF_drho,dF_dgradrho,exc_on_a_grid_wk, ien)
#else
  subroutine ggabek(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dF_drho,dF_dgradrho,ien)
#endif
  use m_Const_Parameters,  only : DP
  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in),optional :: ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)  :: grad_rho(ista_r:iend_r,nspin)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)
  real(kind=DP),intent(out) :: exc
  real(kind=DP),intent(inout) :: dF_drho(ista_r:iend_r,nspin)
  real(kind=DP),intent(out) :: dF_dgradrho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif
  end subroutine ggabek

#ifdef __EDA__
  subroutine ggaprd(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dF_drho,dF_dgradrho,exc_on_a_grid_wk, ien)
#else
  subroutine ggaprd(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dF_drho,dF_dgradrho,ien)
#endif
  use m_Const_Parameters,  only : PAI,DP
  implicit none

  integer,intent(in)         :: nspin,ispin,ista_r,iend_r
  integer,intent(in),optional :: ien
  real(kind=DP),intent(in)   :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)   :: grad_rho(ista_r:iend_r,nspin)
  real(kind=DP),intent(in)   :: wos(ista_r:iend_r)
  real(kind=DP),intent(out)  :: exc
  real(kind=DP),intent(inout):: dF_drho(ista_r:iend_r,nspin)
  real(kind=DP),intent(inout)  :: dF_dgradrho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif
  end subroutine ggaprd

#ifdef __EDA__
  subroutine cr_lda(nspin,ispin,ista_r,iend_r,chgrhr_l,exc,dF_drho,exc_on_a_grid_wk, ien)
#else
  subroutine cr_lda(nspin,ispin,ista_r,iend_r,chgrhr_l,exc,dF_drho,ien)
#endif
  use m_Const_Parameters,  only : PAI,DP
  integer, intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in),optional :: ien
  real(kind=DP), intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP), intent(inout) :: exc
  real(kind=DP), intent(inout) :: dF_drho(ista_r:iend_r,ispin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif
  end subroutine cr_lda

#ifdef __EDA__
  subroutine cr_gga_library( nspin, ispin, ista_r, iend_r, chgrhr_l, grad_trho, wos, &
     &                     exc, dF_drho, exc_on_a_grid_wk, ecor, pot_type, ist, ien )
#else
  subroutine cr_gga_library( nspin, ispin, ista_r, iend_r, chgrhr_l, grad_trho, wos, &
     &                     exc, dF_drho, ecor, pot_type, ist, ien )
#endif
  use m_Const_Parameters,  only : PAI,DP
  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r, pot_type
  integer,intent(in)        :: ist, ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(inout) :: grad_trho(ista_r:iend_r)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)
  real(kind=DP),intent(inout) :: exc
  real(kind=DP),intent(inout) :: dF_drho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif
  real(kind=DP),intent(out),optional :: ecor
  end subroutine cr_gga_library

  subroutine ex_ggapbe_3D(nspin,ispin,chgrhr_l,grad_rho,f2or1,exc,dFx_drho,dFx_dgradrho,nfft_y,iteration,revPBE)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,PAI
  implicit none

  integer,intent(in)        :: nspin,ispin,nfft_y,iteration
  real(kind=DP),intent(in)  :: chgrhr_l(1:nfft_y,ispin)
  real(kind=DP),intent(in)  :: grad_rho(1:nfft_y,nspin)
  real(kind=DP),intent(in)  :: f2or1(1:nfft_y)
  real(kind=DP),intent(out) :: exc
  real(kind=DP),intent(out) :: dFx_drho(1:nfft_y,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(1:nfft_y,nspin)
  logical, intent(in),optional :: revPBE
  end subroutine ex_ggapbe_3D

#if 0
#ifdef LIBXC
#ifdef __EDA__
  subroutine ex_gga_libxc( nspin, ispin, ista_r, iend_r, chgrhr_l, grad_rho, grad_trho, &
       &                   wos, exc, dFx_drho, dFx_dgradrho, exc_on_a_grid_wk, ist, ien )
#else
  subroutine ex_gga_libxc( nspin, ispin, ista_r, iend_r, chgrhr_l, grad_rho, grad_trho, &
     &                   wos, exc, dFx_drho, dFx_dgradrho, ist, ien )
#endif
  use m_Const_Parameters,  only : DP,PAI

  use m_Control_Parameters,  only : xc_func_exch
  use xc_f03_lib_m
  use m_Parallelization,    only : mype

#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif

  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in), optional :: ist, ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)  :: grad_rho(ista_r:iend_r,nspin)
  real(kind=DP),intent(in)  :: grad_trho(ista_r:iend_r)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)
  real(kind=DP),intent(out) :: exc
  real(kind=DP),intent(out) :: dFx_drho(ista_r:iend_r,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif
end subroutine ex_gga_libxc
#endif
#endif

end interface

contains
  integer function check_of_xctype()
    integer :: ii

    if(xctype == 'ggapw91' .or. xctype == 'ldapw91' .or. xctype == 'ggabp  ' &
         & .or. xctype == 'ggapbe ' .or. xctype == 'ldapbe '&
         & .or. xctype == 'ggapbex' .or. xctype == 'vdwdf'  &
         & .or. xctype == 'ggapbey' &
! =========================================== KT_add =============== 13.0A
         & .or. xctype == 'rpbe   ' .or. xctype == 'wc06   ' &!
         & .or. xctype == 'htbs   ' .or. xctype == 'pbesol ' &
         & .or. xctype == 'pbeint ' &
         & .or. xctype == 'ev93   ' .or. xctype == 'evpw91 ' &
         & .or. xctype == 'lb94   ' &
! ================================================================== 13.0A
         & .or. xctype == 'revpbe ' &

         & .or. xctype == 'katopbe' .or. xctype == 'ggapbek' ) then
       check_of_xctype = GGA
    else
       check_of_xctype = LDA
    endif

! ==== KT_add ==== 13.0XX
    if ( xctype == "tb09" ) then
       check_of_xctype = GGA
    endif
#ifdef LIBXC
    if ( xctype == "libxc" ) then
       ii = max( xc_family_exch, xc_family_corr )
       select case (ii)
       case (XC_FAMILY_LDA)
          check_of_xctype = GGA
       case (XC_FAMILY_GGA)
          check_of_xctype = GGA
       case default
          check_of_xctype = GGA
       end select
    endif
#endif
! ================ 13.0XX
  end function check_of_xctype


  subroutine cp_vxc_to_vxco
    vxco_l = vxc_l
  end subroutine cp_vxc_to_vxco



  subroutine m_XC_alloc_s_gga
    allocate(s_gga1(3,3)); s_gga1=0.d0
    allocate(s_gga2(3,3)); s_gga2=0.d0
  end subroutine m_XC_alloc_s_gga

  subroutine m_XC_dealloc_s_gga
    deallocate(s_gga1)
    deallocate(s_gga2)
  end subroutine m_XC_dealloc_s_gga


  subroutine check_of_pcc(nopcc)
    logical, intent(out) :: nopcc
    integer              :: i

    nopcc = .true.

    do i = 1, ntyp
       if(itpcc(i) /= 0) then
          nopcc = .false.
          goto 2201
       endif
    end do

    vxcpc_l = 0.d0
    exc = 0.d0
2201 continue
  end subroutine check_of_pcc

  subroutine afft_allgatherv(afft_mpi0,afft)
!  Upgraded on 19 Sep. 2006 by T. Yamasaki
!    * MPIFFT
    real(kind=DP),intent(in),dimension(:) :: afft_mpi0(ista_fftp:iend_fftp)
    real(kind=DP),intent(out),dimension(:):: afft(nfftp)
#ifdef _MPIFFT_
    real(kind=DP), allocatable, dimension(:) :: afft_tmp
    integer  :: nstart, nend, lx, lxy, llx, lly, llxy, i, j, jj, k, ipin, ipout
    integer  :: nb, datasize, pe_s, pe_r
    integer, allocatable, dimension(:) :: req_s, req_r
    real(kind=DP), allocatable, dimension(:) :: tmp_s, tmp_r
#endif

    integer  :: id_sname = -1
    call tstatc0_begin('afft_allgatherv(in m_XC_Potential) ',id_sname)
    if(npes>=2)  call mpi_barrier(MPI_CommGroup,ierr)

    if(npes_cdfft >= 2) then
       if(myrank_ggacmp < nrank_ggacmp) then
#ifdef _MPIFFT_
          afft = 0.d0
#ifndef _ALLREDUCE_AFFT_ALLGAHTERV_

          datasize = nfftp/npes_cdfft
          allocate(tmp_s(datasize)); tmp_s = 0.d0
          allocate(tmp_r(datasize)); tmp_r = 0.d0
          allocate(req_r(npes-1),req_s(npes-1))

          lx     = fft_box_size_CD_c(1,0)
          lxy    = lx*ly_d
          llx    = fft_box_size_CD_c(1,0)
          lly    = fft_box_size_CD_c(2,0)
          llxy   = llx*lly

          do i = 1, iend_fftp-ista_fftp+1
             tmp_s(i) = afft_mpi0(i+ista_fftp-1)
          end do

          do nb = 1, npes_cdfft-1
             pe_s = mod(myrank_cdfft+nb,           npes_cdfft)
             pe_r = mod(myrank_cdfft-nb+npes_cdfft,npes_cdfft)
             call mpi_irecv(tmp_r,datasize,mpi_double_precision, &
                  & pe_r,pe_r,        mpi_cdfft_world(myrank_ggacmp),req_r(nb),ierr)
             call mpi_isend(tmp_s,datasize,mpi_double_precision, &
                  & pe_s,myrank_cdfft,mpi_cdfft_world(myrank_ggacmp),req_s(nb),ierr)
             call mpi_wait(req_r(nb),istatus,ierr)
             call mpi_wait(req_s(nb),istatus,ierr)

             nstart = pe_r*ny_d+1
             nend   = min(fft_box_size_CD(2,1),(pe_r+1)*ny_d)

             if(kimg == 1) then
                do k = 1, fft_box_size_CD(3,1)
                   do j = nstart, nend
                      jj = j - nstart + 1
                      do i = 1, fft_box_size_CD(1,1)
                         ipin  = i + (jj-1)* lx + (k-1)* lxy
                         ipout = i + ( j-1)*llx + (k-1)*llxy
                         afft(ipout)   = tmp_r(ipin)
                      end do
                   end do
                end do
             else if(kimg == 2) then
                do k = 1, fft_box_size_CD(3,1)
                   do j = nstart, nend
                      jj = j - nstart + 1
                      do i = 1, fft_box_size_CD(1,1)
                         ipin  = 2*(i + (jj-1)* lx + (k-1)* lxy)-1
                         ipout = 2*(i + ( j-1)*llx + (k-1)*llxy)-1
                         afft(ipout)   = tmp_r(ipin)
                         afft(ipout+1) = tmp_r(ipin+1)
                      end do
                   end do
                end do
             end if
          end do
          deallocate(req_s,req_r)
          deallocate(tmp_s)
          deallocate(tmp_r)

! diagonal part
          nstart = myrank_cdfft*ny_d+1
          nend   = min(fft_box_size_CD(2,1),(myrank_cdfft+1)*ny_d)
!!$          allocate(afft_tmp(nfftp))

          if(kimg == 1) then
             do k = 1, fft_box_size_CD(3,1)
                do j = nstart, nend
                   jj = j - nstart + 1
                   do i = 1, fft_box_size_CD(1,1)
                      ipin  = ista_fftp-1 + (i + (jj-1)* lx + (k-1)* lxy)
                      ipout =               (i + ( j-1)*llx + (k-1)*llxy)
                      afft(ipout)   = afft_mpi0(ipin)
                   end do
                end do
             end do
          else if(kimg == 2) then
             do k = 1, fft_box_size_CD(3,1)
                do j = nstart, nend
                   jj = j - nstart + 1
                   do i = 1, fft_box_size_CD(1,1)
                      ipin  = ista_fftp-1 + 2*(i + (jj-1)* lx + (k-1)* lxy)-1
                      ipout =               2*(i + ( j-1)*llx + (k-1)*llxy)-1
                      afft(ipout)   = afft_mpi0(ipin)
                      afft(ipout+1) = afft_mpi0(ipin+1)
                   end do
                end do
             end do
          end if

#else
          nstart = myrank_cdfft*ny_d+1
          nend   = min(fft_box_size_CD(2,1),(myrank_cdfft+1)*ny_d)
          lx     = fft_box_size_CD_c(1,0)
          lxy    = lx*ly_d
          llx    = fft_box_size_CD_c(1,0)
          lly    = fft_box_size_CD_c(2,0)
          llxy   = llx*lly
          allocate(afft_tmp(nfftp))

          if(kimg == 1) then
             do k = 1, fft_box_size_CD(3,1)
                do j = nstart, nend
                   jj = j - nstart + 1
                   do i = 1, fft_box_size_CD(1,1)
                      ipin  = ista_fftp-1 + (i + (jj-1)* lx + (k-1)* lxy)
                      ipout =               (i + ( j-1)*llx + (k-1)*llxy)
                      afft(ipout)   = afft_mpi0(ipin)
                   end do
                end do
             end do
          else if(kimg == 2) then
             do k = 1, fft_box_size_CD(3,1)
                do j = nstart, nend
                   jj = j - nstart + 1
                   do i = 1, fft_box_size_CD(1,1)
                      ipin  = ista_fftp-1 + 2*(i + (jj-1)* lx + (k-1)* lxy)-1
                      ipout =               2*(i + ( j-1)*llx + (k-1)*llxy)-1
                      afft(ipout)   = afft_mpi0(ipin)
                      afft(ipout+1) = afft_mpi0(ipin+1)
                   end do
                end do
             end do
          end if
          call mpi_allreduce(afft, afft_tmp, nfftp, mpi_double_precision, mpi_sum &
               &                                  , mpi_cdfft_world(myrank_ggacmp),ierr)
          afft = afft_tmp
#endif

#else
          call mpi_allgatherv(afft_mpi0,nel_fftp(myrank_cdfft),mpi_double_precision &  ! MPI
               &       ,afft,nel_fftp,idisp_fftp,mpi_double_precision,mpi_cdfft_world(myrank_ggacmp),ierr)
#endif
       end if
    else
       afft = afft_mpi0
    end if

    if(nrest_cdfft >= 1) then
       if(mype > npes_cdfft*nrank_ggacmp-1) then
          call mpi_recv(afft,nfftp,mpi_double_precision &
               & , mype-npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,istatus,ierr)
       end if
       if(mype < nrest_cdfft) then
          call mpi_send(afft,nfftp,mpi_double_precision &
               & , mype+npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,ierr)
       end if
    end if

    call tstatc0_end(id_sname)

  end subroutine afft_allgatherv


!===============================================================================

  subroutine m_XC_alloc_vxc_3D

    allocate(vxc_l(ista_kngp:iend_kngp,kimg,nspin)); vxc_l = 0.d0
    allocate(vxco_l(ista_kngp:iend_kngp,kimg,nspin)); vxco_l = 0.d0
    allocate(vxcpc_l(ista_kngp:iend_kngp,kimg)); vxcpc_l = 0.d0
  end subroutine m_XC_alloc_vxc_3D

!===============================================================================
  subroutine m_XC_dealloc_vxc_3D
    if(allocated(vxc_l)) deallocate(vxc_l)
    if(allocated(vxco_l)) deallocate(vxco_l)
    if(allocated(vxcpc_l)) deallocate(vxcpc_l)
  end subroutine m_XC_dealloc_vxc_3D

!===============================================================================

!!#ifdef __EDA__
!!! -----  ascat starts modifying  -----
!!  subroutine m_XC_alloc_exc_on_a_grid()
!!    use m_FFT, only : nfftp_nonpara
!!    integer :: lsize
!!    allocate(exc_on_a_grid(nfftp_nonpara)); exc_on_a_grid = 0.d0
!!  end subroutine m_XC_alloc_exc_on_a_grid
!!
!!  subroutine m_XC_dealloc_exc_on_a_grid()
!!    deallocate(exc_on_a_grid)
!!  end subroutine m_XC_dealloc_exc_on_a_grid
!!! -----  ascat ceases modifying  -----
!!#endif

  subroutine xc_allocate_3D(ispin,nfout)
    use m_Parallelization,      only : np_fftcd_y, mp_fftcd_y   &
   &                                 , np_fftcd_x, mp_fftcd_x   &
   &                                 , np_fftcd_z, mp_fftcd_z   &
   &                                 , xyz_fftcd_y, xyz_fftcd_x, xyz_fftcd_z

    integer, intent(in)      :: ispin,nfout

    integer                  :: idp,nlp,nmp,nnp,nd2p,nd3p,ip, idph, nlph
    integer                  :: n, nn, i, j, k
    integer                  :: np0, j0, i2, j2, k2, n0
    integer  :: lx,ly,lz,mx,my,mz,i1,mm,nx,ny,nz,ierr
!XX!integer, allocatable, dimension(:) :: wk_mp_fft_y

                                                  __TIMER_SUB_START(751)
    allocate(chgrhr_l(1:nfft_div_size,ispin), stat=ierr)    ! MPI
     if (ierr /= 0) then
        write(nfout,*)' xc_allocate_3D :  Not allocated'
        call flush(nfout)
        call mpi_abort(mpi_comm_world,499,ierr)
     end if
    allocate(   f2or1(1:nfft_div_size), stat=ierr)
     if (ierr /= 0) then
        write(nfout,*)' xc_allocate_3D :  Not allocated'
        call flush(nfout)
        call mpi_abort(mpi_comm_world,500,ierr)
     end if

    nlp  = fft_box_size_CD_3D(1,1)
    nmp  = fft_box_size_CD_3D(2,1)
    nnp  = fft_box_size_CD_3D(3,1)
    idp  = fft_box_size_CD_3D(1,0)
    nd2p = fft_box_size_CD_3D(2,0)
    nd3p = fft_box_size_CD_3D(3,0)

    lx = fft_box_size_CD_3D(1,0)
    ly = fft_box_size_CD_3D(2,0)
    lz = fft_box_size_CD_3D(3,0)

    if (nfft_div_size > 0) then
#ifdef FFT_3D_DIVISION_CD
       nx = xyz_fftcd_x(2,1)-xyz_fftcd_x(1,1)+1
       ny = xyz_fftcd_x(2,2)-xyz_fftcd_x(1,2)+1
       nz = xyz_fftcd_x(2,3)-xyz_fftcd_x(1,3)+1
#else
       if (sw_fft_xzy > 0) then
          nx = xyz_fftcd_y(2,1)-xyz_fftcd_y(1,1)+1
          ny = xyz_fftcd_y(2,2)-xyz_fftcd_y(1,2)+1
          nz = xyz_fftcd_y(2,3)-xyz_fftcd_y(1,3)+1
       else
          nx = xyz_fftcd_z(2,1)-xyz_fftcd_z(1,1)+1
          ny = xyz_fftcd_z(2,2)-xyz_fftcd_z(1,2)+1
          nz = xyz_fftcd_z(2,3)-xyz_fftcd_z(1,3)+1
       end if
#endif

       call set_f2or1_3D(npes_cdfft,ista_fftph,iend_fftph)
    end if

    if(check_of_xctype() == GGA) then
       allocate(inx(1:np_fftcd_x*kimg)); inx = 0
       allocate(jnx(1:np_fftcd_x*kimg)); jnx = 0
       allocate(knx(1:np_fftcd_x*kimg)); knx = 0

       if (np_fftcd_x /= 0) then
          nx = xyz_fftcd_x(2,1)-xyz_fftcd_x(1,1)+1
          ny = xyz_fftcd_x(2,2)-xyz_fftcd_x(1,2)+1
          nz = xyz_fftcd_x(2,3)-xyz_fftcd_x(1,3)+1

                                                  __TIMER_DO_START(876)
!OCL NORECURRENCE
          do n = 1, np_fftcd_x
             i1 = mp_fftcd_x(n)
             mz = (i1-1)/(lx*ly)+1
             mm = mod(i1,(lx*ly))
             if (mm==0) mm=lx*ly
             my = (mm-1)/lx+1
             mx = mod(mm,lx)
             if (mx==0) mx = lx
             inx(n*kimg) = mx - 1
             if(2*inx(n*kimg) > nlp) inx(n*kimg) = inx(n*kimg) - nlp
             jnx(n*kimg) = my - 1
             if(2*jnx(n*kimg) > nmp) jnx(n*kimg) = jnx(n*kimg) - nmp
             knx(n*kimg) = mz - 1
             if(2*knx(n*kimg) > nnp) knx(n*kimg) = knx(n*kimg) - nnp
             if (kimg==2) then
                inx(n*kimg-1) = inx(n*kimg)
                jnx(n*kimg-1) = jnx(n*kimg)
                knx(n*kimg-1) = knx(n*kimg)
             end if
          end do
                                                  __TIMER_DO_STOP(876)
       end if

#ifdef FFT_3D_DIVISION_CD
       allocate(    chden_l(1:np_fftcd_x*2,nspin), stat=ierr)
#else
       allocate(    chden_l(1:np_fftcd_x*kimg,nspin), stat=ierr)
#endif

       allocate(  grad_trho(1:nfft_div_size), stat=ierr)
       allocate(   grad_rho(1:nfft_div_size,nspin), stat=ierr)
       allocate(    dF_drho(1:nfft_div_size,nspin), stat=ierr)      ! MPI
       allocate(dF_dgradrho(1:nfft_div_size,nspin), stat=ierr)  ! MPI
       allocate(  cgrad_rho(1:nfft_div_size,3,nspin), stat=ierr)  ! MPI
        if (ierr /= 0) then
           write(nfout,*)' xc_allocate_3D :  Not allocated array'
           call flush(nfout)
           call mpi_abort(mpi_comm_world,501,ierr)
        end if
       chden_l = 0.d0 ! MPI
       grad_rho = 0.d0  ! MPI
       cgrad_rho = 0.0d0
       dF_drho = 0.0d0
       dF_dgradrho = 0.0d0
    end if
                                                  __TIMER_SUB_STOP(751)

    if ( use_metagga ) then
       allocate( ekin_dens(1:nfft_div_size,nspin) ); ekin_dens = 0.0d0
       allocate( grad2_rho(1:nfft_div_size,nspin) ); grad2_rho = 0.0d0
       if ( vtau_exists ) then
          allocate( dF_dtau(1:nfft_div_size,nspin) ); dF_dtau = 0.0d0
          allocate( dF_dlaplrho(1:nfft_div_size,nspin) ); dF_dlaplrho = 0.0d0
       endif
    endif

  contains
!------------------------------------------------------------------------------
    subroutine set_f2or1_3D(npes,ista,iend)
      integer, intent(in) :: npes,ista, iend
      integer :: idph,nlph,ip,i,j,k,md,inc,loop

                                                  __TIMER_SUB_START(752)
#ifdef FFT_3D_DIVISION_CD
      if(iprixc >= 2 ) write(nfout,'(" ix kimg = 2 <<set_f2or1>>")')
      f2or1 = 1.d0
                                                  __TIMER_DO_START(878)
       do n = 1, np_fftcd_x
          i1 = mp_fftcd_x(n)
          mz = (i1-1)/(lx*ly)+1
          mm = mod(i1,(lx*ly))
          if (mm==0) mm=lx*ly
          my = (mm-1)/lx+1
          mx = mod(mm,lx)
          if (mx==0) mx = lx

          if (nlp < mx) then
             f2or1(n) = 0.d0
             cycle
          end if
          if (nmp < my) then
             f2or1(n) = 0.d0
             cycle
          end if
          if (nnp < mz) then
             f2or1(n) = 0.d0
             cycle
          end if
       end do
                                                  __TIMER_DO_STOP(878)
#else
      if(kimg == 1) then
         idph = idp/2
         nlph = nlp/2
!!       idph = lx

         f2or1 = 0.d0
!        do k = 1, min(nz_d, nnp-nz_d*myrank_cdfft)
!           do j = 1, nmp
!              do i = 1, nlph
!                 ip = i + idph*(j-1) + idph*ly*(k-1) + idph*ly*lz_d*myrank_cdfft
!                 f2or1(ip) = 2.d0
!              end do
!                 ip = 1 + idph*(j-1) + idph*ly*(k-1) + idph*ly*lz_d*myrank_cdfft
!              f2or1(ip) = 1.d0
!              ip = nlph+1 + idph*(j-1) + idph*ly*(k-1) + idph*ly*lz_d*myrank_cdfft
!              f2or1(ip) = 1.d0
!           end do
!        end do
         f2or1 = 0.d0
          inc = 0
                                                  __TIMER_DO_START(877)
          if (sw_fft_xzy > 0) then
             loop = np_fftcd_y
          else
             loop = np_fftcd_z
          end if
          do n = 1, loop , 2
             inc = inc + 1
!XX!         i1 = wk_mp_fft_y(n)
             if (sw_fft_xzy > 0) then
                i1 = mp_fftcd_y(n)
             else
                i1 = mp_fftcd_z(n)
             end if
             mz = (i1-1)/(lx*ly)+1
             mm = mod(i1,(lx*ly))
             if (mm==0) mm=lx*ly
             my = (mm-1)/lx+1
             mx = mod(mm,lx)
             if (mx==0) mx = lx

             md = (mx-1)/2+1
             md = mod(md,idph)
             if (md > 1) then
                f2or1(inc) = 2.0d0
             else
                f2or1(inc) = 1.0d0
             endif
!!           md = mod(mx,idph)
!!           if ((2 < md) .and. ( md < (idph - 1))) then
!!              f2or1(n) = 2.0d0
!!           else
!!              f2or1(n) = 1.0d0
!!           end if
          end do
                                                  __TIMER_DO_STOP(877)
      else
         if(iprixc >= 2 ) write(nfout,'(" ix kimg = 2 <<set_f2or1>>")')
         f2or1 = 1.d0
                                                  __TIMER_DO_START(878)
          if (sw_fft_xzy > 0) then
             loop = np_fftcd_y
          else
             loop = np_fftcd_z
          end if
          do n = 1, loop
             if (sw_fft_xzy > 0) then
                i1 = mp_fftcd_y(n)
             else
                i1 = mp_fftcd_z(n)
             end if
             mz = (i1-1)/(lx*ly)+1
             mm = mod(i1,(lx*ly))
             if (mm==0) mm=lx*ly
             my = (mm-1)/lx+1
             mx = mod(mm,lx)
             if (mx==0) mx = lx

             if (nlp < mx) then
                f2or1(n) = 0.d0
                cycle
             end if
             if (nmp < my) then
                f2or1(n) = 0.d0
                cycle
             end if
             if (nnp < mz) then
                f2or1(n) = 0.d0
                cycle
             end if
          end do
                                                  __TIMER_DO_STOP(878)
      end if
#endif
                                                  __TIMER_SUB_STOP(752)
    end subroutine set_f2or1_3D
!------------------------------------------------------------------------------
  end subroutine xc_allocate_3D

!===============================================================================

  subroutine xc_deallocate_3D
#ifndef _MPIFFT_
    call m_FFT_dealloc_CD_box()
#endif
    deallocate(afft_l)
    deallocate(chgrhr_l)
    deallocate(f2or1)

    if(check_of_xctype() == GGA) then
       deallocate(inx); deallocate(jnx); deallocate(knx)
       deallocate(chden_l)
       deallocate(grad_trho)
       deallocate(cgrad_rho)
       deallocate(grad_rho); deallocate(dF_drho)
       deallocate(dF_dgradrho)
    end if
    if ( use_metagga ) then
       deallocate( ekin_dens );       deallocate( grad2_rho )
       if ( allocated( dF_dtau ) ) deallocate( dF_dtau )
       if ( allocated( dF_dlaplrho ) ) deallocate( dF_dlaplrho )
    endif
  end subroutine xc_deallocate_3D

!===============================================================================

  subroutine m_XC_cal_potential_3D(nfout,input_charge,chgq_l, vflag  &
    &  , chgsoft,hsr,hsrd)
!   Upgraded on 19th Sep. 2006 by T. Yamasaki
!     * MPIFFT
    use m_Parallelization,     only : nel_fftcd_x, nel_fftcd_y, nel_fftcd_z &
   &                                , xyz_fftcd_x, xyz_fftcd_y, xyz_fftcd_z &
   &                                , mp_fftcd_x , mp_fftcd_y , mp_fftcd_z  &
   &                                , np_fftcd_x , np_fftcd_y , np_fftcd_z  &
   &                                , mpi_ke_world, mpi_kg_world   &
   &                                , chgq_fftcd_scnt, chgq_fftcd_rcnt &
   &                                , chgq_fftcd_send, chgq_fftcd_recv &
   &                                , chgq_fftcd_index, chgq_fftcd_dist &
   &                                , chgq_fftcd_maxrecv, chgq_fftcd_maxsend &
   &                                , fftcd_chgq_scnt, fftcd_chgq_rcnt &
   &                                , fftcd_chgq_send, fftcd_chgq_recv &
   &                                , fftcd_chgq_index, fftcd_chgq_dist &
   &                                , fftcd_chgq_maxrecv, fftcd_chgq_maxsend &
   &                                , fftcd_X_x_nel , fftcd_X_y_nel , fftcd_X_z_nel  &
   &                                , igfp_full
    use m_FFT,                 only : m_FFT_CD_Direct_3D, m_FFT_CD_Inverse_3D           &
#ifdef FFT_3D_DIVISION_CD
   &                                , m_FFT_CD_Direct_3DIV_3D, m_FFT_CD_Inverse_3DIV_3D &
#endif
   &                                , m_FFT_CD_Direct_XYZ_3D, m_FFT_CD_Inverse_XYZ_3D,m_FFT_CD_inverse0
    use m_IterationNumbers,     only : iteration

    integer, intent(in) :: nfout, input_charge
    real(DP),intent(in) :: chgq_l(ista_kngp:iend_kngp,kimg,nspin)
    !                                    charge density
    integer, intent(in) :: vflag
    real(DP),intent(in),optional,dimension(ista_kngp:iend_kngp,kimg,nspin) &
         &              :: chgsoft  ! soft charge density
    real(DP),intent(in),optional,dimension(natm,nlmt,nlmt,nspin)     :: hsr
    real(DP),intent(in),optional,dimension(natm,nlmt,nlmt,nspin,3,3) :: hsrd

    integer, dimension(0:nrank_g-1) ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)::sta_r, sta_s
    integer :: nfftsize, lsize, isrsize, fft_l_size, i, ierr
    real(kind=DP), allocatable, dimension(:,:) :: chgq_plus
    integer       :: ispin, iloop
    logical       :: nopcc
    real(kind=DP) :: rinplw

    logical,save                   :: firstcall = .true.
    integer       :: id_sname = -1
#ifdef __FAPP__
    call fapp_start('xc_cal_potential',1,1)
#endif

!    iteration = 0
!    read(2000) chgrhr_l

    if(sw_communicator_for_chg == ON) then
       if(present(chgsoft)) then
         call m_XC_cal_potential(nfout,input_charge,chgq_l,vflag,chgsoft,hsr,hsrd)
       else
         call m_XC_cal_potential(nfout,input_charge,chgq_l,vflag)
       endif
       return
    endif
                                                  __TIMER_SUB_START(750)
    call tstatc0_begin('m_XC_cal_potential_3D ',id_sname,level=1)
    rinplw = 1.d0/product(fft_box_size_CD_3D(1:3,1))

    if(input_charge == Partial_Core_Charge) then
       call check_of_pcc(nopcc)
       if(nopcc) goto 4001
    endif

#ifdef FFT_3D_DIVISION_CD
     nfft_div_size = np_fftcd_x
#else
    if (sw_fft_xzy > 0) then
       if (kimg == 1) then
          nfft_div_size = np_fftcd_y / 2
       else
          nfft_div_size = np_fftcd_y
       end if
    else
       if (kimg == 1) then
          nfft_div_size = np_fftcd_z / 2
       else
          nfft_div_size = np_fftcd_z
       end if
    end if
#endif

    call set_ispin(ispin)
    call xc_allocate_3D(ispin,nfout)

    nfftsize = fft_box_size_CD_3D(1,0) * fft_box_size_CD_3D(2,0) * fft_box_size_CD_3D(3,0)
#ifdef FFT_3D_DIVISION_CD
    lsize = fftcd_X_x_nel*fftcd_X_y_nel*fftcd_X_z_nel
#else
    lsize = max(maxval(nel_fftcd_x(:)),maxval(nel_fftcd_y(:)),maxval(nel_fftcd_z(:)))
#endif
    isrsize = min(lsize,mp_kngp)
!   fft_l_size  = nel_fftcd_y(myrank_g)
    fft_l_size  = nel_fftcd_x(myrank_g)
#ifdef FFT_3D_DIVISION_CD
    allocate(afft_l(lsize*2,1), stat=ierr)
#else
    allocate(afft_l(lsize*kimg,1), stat=ierr)
#endif
     if(ierr /= 0) then
        write(nfout,*)' m_XC_cal_potential_3D : Not allocated afft_l array'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 251, ierr)
     endif

    afft_l = 0.0d0
    call set_chgrhr_case_collinear( ispin )

    if ( use_metagga .and. .not. use_modeled_ekin_density ) then
       if ( noncol ) then
       else
          if ( use_symm_ekin_density ) then
             call set_ekindens_case_collinear_3D( ispin, ekins_l, ekin_dens )
          endif
       endif
    endif

    call check_lmn_even_3D  ! -(contained in subr. m_XC_cal_potential) ->lmn_even

    if(check_of_xctype() == GGA) then
       if(vflag == STRESS_) then
! === This case is NOT supported on 3D_Parallel! by tkato 2014/01/30 ===========
!!          if(mype == 0) then
!!             write(0,*) '======================================================================'
!!             write(0,*) '= ERROR!!!                                                           ='
!!             write(0,*) '= This case ("check_of_xctype() == GGA" and "vflag == STRESS_")      ='
!!             write(0,*) '= is NOT supported on 3D_Parallel because MPI communicators          ='
!!             write(0,*) '= mpi_cdfft_world and mpi_ggacmp_cross_world used in ggaxcp_diff_3D  ='
!!             write(0,*) '= are NOT made.                                                      ='
!!             write(0,*) '======================================================================'
!!          end if
!!          call MPI_Barrier(MPI_CommGroup,ierr)
!!          call MPI_Abort(MPI_CommGroup,999,ierr)
! ==============================================================================
          call ggaxcp_diff_3D           ! -(contained here) ->chgrhr,exc,etc.
          !    ~~~~~~~~~~~
       else
          call ggaxcp_3D()              ! -(contained here) ->chgrhr,exc
          !    ~~~~~~
       endif
    else
       call xcpotf_3D(ispin,input_charge) ! -(m_XC_Potential) ->chgrhr,exc
       !    ~~~~~~
    end if
    if (vflag == VXC_AND_EXC .or. vflag == STRESS_) then
       call set_vxc_case_collinear(ispin)
    endif

    call xc_deallocate_3D
#ifdef __FAPP__
    call fapp_stop('xc_cal_potential',1,1)
#endif
4001 continue

!      write(5000+mype,*) "exc = ", exc
!      write(5100+mype,*) dF_drho(1:5,1)
!      write(5200+mype,*) dF_dtau(1:5,1)
!      write(5300+mype,*) vxc_l(1:20,1,1)
!      write(5400+mype,*) vxc_l(1:20,2,1)

!      stop

    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(750)
  contains

    subroutine set_vxc_case_collinear(ispin)
      integer, intent(in) :: ispin
      integer :: iloop
      real(kind=DP), allocatable :: bfft_l(:,:)

      if ( use_metagga .and. vtau_exists ) then
#ifdef FFT_3D_DIVISION_CD
         allocate(bfft_l(lsize*2,1), stat=ierr)
#else
         allocate(bfft_l(lsize*kimg,1), stat=ierr)
#endif
      endif

      do iloop = 1, ispin
!         call cpafft_3D(iloop)                     ! chgrhr_l -> afft
         call cp_vals_to_afft_rspace_3D( iloop, chgrhr_l, afft_l )
#ifdef FFT_3D_DIVISION_CD
         call m_FFT_CD_Direct_3DIV_3D(nfout,afft_l,lsize,1)
#else
         if (sw_fft_xzy > 0) then
            call m_FFT_CD_Direct_3D(nfout,afft_l,lsize,1)
         else
            call m_FFT_CD_Direct_XYZ_3D(nfout,afft_l,lsize,1)
         end if
#endif
         if ( use_metagga .and. vtau_exists ) then
            call cp_vals_to_afft_rspace_3D( iloop, dF_dtau, bfft_l )
#ifdef FFT_3D_DIVISION_CD
            call m_FFT_CD_Direct_3DIV_3D(nfout,bfft_l,lsize,1)
#else
            if (sw_fft_xzy > 0) then
               call m_FFT_CD_Direct_3D(nfout,bfft_l,lsize,1)
            else
               call m_FFT_CD_Direct_XYZ_3D(nfout,bfft_l,lsize,1)
            end if
#endif
         endif
         call map_fftcd_onto_charge( iloop, lsize, 1, afft_l, isrsize, fft_l_size, &
              &                      nfout, input_charge )
         if ( vtau_exists ) then
            call map_fftcd_onto_charge( iloop, lsize, 1, bfft_l, isrsize, fft_l_size, &
                 &                      nfout, -1 )
         endif
      end do
      if(iprixc >= 2) call wd_small_part_of_vxc
      if ( allocated(bfft_l) ) deallocate( bfft_l )
    end subroutine set_vxc_case_collinear

    subroutine set_chgrhr_case_collinear( ispin )
      integer, intent(in) :: ispin
      integer :: iloop

      do iloop = 1, ispin
!!       call map_charge_onto_a_fftcd_box(iloop)     !-(m_XC_Pot.) -> afft(*) (xcchg2)
         call map_charge_onto_a_fftcd_box( chgq_l, nspin, iloop, input_charge )

         if(check_of_xctype() == GGA) then
#ifdef FFT_3D_DIVISION_CD
            integer :: lx, ly, lz, nxyz, i1, mx, my, mz, mm, jadd, kx1p, kx2p, kx3p
            kx1p = fftcd_X_x_nel
            kx2p = fftcd_X_y_nel
            kx3p = fftcd_X_z_nel
            lx = fft_box_size_CD_3D(1,0)
            ly = fft_box_size_CD_3D(2,0)
            lz = fft_box_size_CD_3D(3,0)
            do nxyz = 1, np_fftcd_x
               i1 = mp_fftcd_x(nxyz)
               mz = (i1-1)/(lx*ly)+1
               mm = mod(i1,(lx*ly))
               if (mm==0) mm=lx*ly
               my = (mm-1)/lx+1
               mx = mod(mm,lx)
               if (mx==0) mx = lx
               jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
               chden_l(nxyz*2-1,iloop) = afft_l(jadd*2-1,1)
               chden_l(nxyz*2  ,iloop) = afft_l(jadd*2  ,1)
            end do
#else
            chden_l(1:np_fftcd_x*kimg,iloop) = afft_l(1:np_fftcd_x*kimg,1)
#endif
         end if
#ifdef FFT_3D_DIVISION_CD
         call m_FFT_CD_Inverse_3DIV_3D(nfout,afft_l,lsize,1)
#else
         if (sw_fft_xzy > 0) then
            call m_FFT_CD_Inverse_3D(nfout,afft_l,lsize,1)
         else
            call m_FFT_CD_Inverse_XYZ_3D(nfout,afft_l,lsize,1)
         end if
#endif
         call check_of_negative_CD_3D(nfout, ispin)
         call cp_afft_to_chgrhr_3D( iloop )      ! -(contained here) afft -> chgrhr_l
      end do

    end subroutine set_chgrhr_case_collinear

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!!    subroutine map_charge_onto_a_fftcd_box(is)
    subroutine map_charge_onto_a_fftcd_box( chgq_l, ndim_dens, iloop, operation_mode )
      integer, intent(in) :: iloop, ndim_dens, operation_mode
      real(kind=DP), intent(in) :: chgq_l(ista_kngp:iend_kngp,kimg,ndim_dens)
      real(kind=DP)      :: fac
      integer            :: mm, it, i, j
      real(kind=DP), allocatable, dimension(:,:) :: chgqplus_l
      real(kind=DP), allocatable, dimension(:,:) :: sendbuf, recvbuf
      integer :: icnt_send, icnt_recv, lrank
      integer :: iend, iadd, mpi_comm_all

      integer :: ia
      real(kind=DP) :: weight, grt

      integer, parameter :: itag = 19
! === DEBUG by tkato 2012/06/05 ================================================
#ifdef USE_ALLTOALLV
      integer, allocatable, dimension(:) :: sdsp, rdsp
#endif
! ==============================================================================
                                                  __TIMER_SUB_START(753)

      allocate(chgqplus_l(ista_kngp:iend_kngp,kimg))

      if(operation_mode == Valence_plus_PC_Charge) then
         fac = 1.d0/nspin
      else if ( operation_mode == Valence_Charge_Only ) then
         fac = 1.d0
      else if(operation_mode == Partial_Core_Charge) then
         fac = 1.d0
      else
        call phase_error_with_msg(nfout,' Error of ixc in <<<map_charge_onto_a_fft_box2>>>',__LINE__,__FILE__)
      endif

      chgqplus_l = 0.d0
      mm = 0
                                                  __TIMER_DO_START(879)
      if ( operation_mode == Valence_Charge_Only ) goto 1200

      do it = 1, ntyp
         if(itpcc(it) == 0) cycle
         mm = mm + 1
         do j = 1, kimg
            !do i = ista_kngp, iend_kngp
            do i = ista_kngp, min(kgp_reduced, iend_kngp)  !for mpi
               chgqplus_l(i,j) = chgqplus_l(i,j) + fac*zfm3_l(i,it,j)*rhpcg_l(i,mm)
            end do
         end do
      end do
                                                  __TIMER_DO_STOP(879)
                                                  __TIMER_DO_START(880)
1200  continue

      if ( operation_mode == Valence_plus_PC_Charge &
           &     .or. operation_mode == Valence_Charge_Only ) then
         do j = 1, kimg
            !do i = ista_kngp, iend_kngp
            do i = ista_kngp, min(kgp_reduced, iend_kngp)  !for mpi
               chgqplus_l(i,j) = chgqplus_l(i,j) + chgq_l(i,j,iloop)
            end do
         end do
      end if
                                                  __TIMER_DO_STOP(880)

!=== OPENCORE ===
      if ( nspin == 1 .or. sw_opencore == OFF ) goto 1500
      if ( sw_xc_opencore_ae_only == ON ) goto 1500

      Do ia=1, natm
         it = ityp(ia)
!         if (itpcc(it) == 0) cycle
         if ( has_opencore(it) == 0 ) cycle

         if ( noncol ) then
            if ( iloop == 1 ) cycle
            if ( iloop == 2 ) then
               weight = mag_opencore_pol(ia,1)
            else if ( iloop == 3 ) then
               weight = mag_opencore_pol(ia,2)
            else if ( iloop == 4 ) then
               weight = mag_opencore_pol(ia,3)
            endif
            do i = ista_kngp, iend_kngp  !for mpi
               grt = (pos(ia,1)*ngabc(i,1) + pos(ia,2)*ngabc(i,2) &
                    &    + pos(ia,3)*ngabc(i,3))*PAI2
               chgqplus_l(i,1) = chgqplus_l(i,1) &
                    &          + dcos(grt)*iwei(ia) *rmag_opencore_l(i,it) *weight
               chgqplus_l(i,2) = chgqplus_l(i,2) &
                    &          - dsin(grt)*iwei(ia) *rmag_opencore_l(i,it) *weight
            end do

         else
            if ( iloop == 1 ) then
               weight = 0.5d0 *mag_opencore_pol(ia,1)
            else
               weight = -0.5d0 *mag_opencore_pol(ia,1)
            endif
            if ( kimg == 1 ) then
               do i = ista_kngp, iend_kngp  !for mpi
                  grt = (pos(ia,1)*ngabc(i,1) + pos(ia,2)*ngabc(i,2) &
                       &    + pos(ia,3)*ngabc(i,3))*PAI2
                  chgqplus_l(i,1) = chgqplus_l(i,1) &
                       &          + dcos(grt)*iwei(ia) *rmag_opencore_l(i,it) *weight
               end do
            else
               do i = ista_kngp, iend_kngp  !for mpi
                  grt = (pos(ia,1)*ngabc(i,1) + pos(ia,2)*ngabc(i,2) &
                       &    + pos(ia,3)*ngabc(i,3))*PAI2
                  chgqplus_l(i,1) = chgqplus_l(i,1) &
                       &          + dcos(grt)*iwei(ia) *rmag_opencore_l(i,it) *weight
                  chgqplus_l(i,2) = chgqplus_l(i,2) &
                       &          - dsin(grt)*iwei(ia) *rmag_opencore_l(i,it) *weight
               end do
            end if
         end if
      End do
1500  continue

       mpi_comm_all = mpi_ke_world

       if (chgq_fftcd_maxsend /= 0) then
          allocate(sendbuf(chgq_fftcd_maxsend*kimg,0:nrank_g-1))
          sendbuf = 0.0d0
                                                  __TIMER_DO_START(881)
          if (kimg == 1) then
             iend = iend_kngp
!!           if (iend > kg) iend = kg
             if (iend > kgp_reduced) iend = kgp_reduced
             do i = ista_kngp, iend
                iadd = i-ista_kngp+1
#ifdef CD_FFT_ALL
                if (chgq_fftcd_dist(iadd) < 0) cycle
#endif
                sendbuf(chgq_fftcd_index(iadd),chgq_fftcd_dist(iadd)) = chgqplus_l(i,1)
             enddo
          else
             iend = iend_kngp
!            if (iend > kg) iend = kg
!             if (iend > kgp) iend = kgp
             if (iend > kgp_reduced) iend = kgp_reduced
             do i = ista_kngp, iend
                iadd = i-ista_kngp+1
#ifdef CD_FFT_ALL
                if (chgq_fftcd_dist(iadd) < 0) cycle
#endif
                sendbuf(chgq_fftcd_index(iadd)*2-1,chgq_fftcd_dist(iadd)) = chgqplus_l(i,1)
                sendbuf(chgq_fftcd_index(iadd)*2,  chgq_fftcd_dist(iadd)) = chgqplus_l(i,2)
             enddo
          endif
                                                  __TIMER_DO_STOP(881)
       endif
       if (chgq_fftcd_maxrecv /= 0) then
          allocate(recvbuf(chgq_fftcd_maxrecv*kimg,0:nrank_g-1))
          recvbuf = 0.0d0
       endif

       deallocate(chgqplus_l)

#ifndef USE_ALLTOALLV
       icnt_recv = 0
                                                 __TIMER_COMM_START_w_BARRIER(mpi_comm_all,882)
       do lrank = 0, nrank_g - 1
          if (chgq_fftcd_rcnt(lrank) /= 0) then
             call mpi_irecv(recvbuf(1,lrank), chgq_fftcd_rcnt(lrank)*kimg, mpi_double_precision, &
         &                  lrank, itag, mpi_comm_all, req_r(icnt_recv), ierr)
              if (ierr /= 0) then
                 write(nfout,*)' map_charge_onto_a_fftcd_box :  mpi_irecv error'
                 call flush(nfout)
                 call mpi_abort(mpi_comm_world,162,ierr)
              endif
             icnt_recv = icnt_recv + 1
          endif
       enddo
       icnt_send = 0
       do lrank = 0, nrank_g - 1
          if (chgq_fftcd_scnt(lrank) /= 0) then
             call mpi_isend(sendbuf(1,lrank), chgq_fftcd_scnt(lrank)*kimg, mpi_double_precision, &
         &                  lrank, itag, mpi_comm_all, req_s(icnt_send), ierr)
              if (ierr /= 0) then
                 write(nfout,*)' map_charge_onto_a_fftcd_box :  mpi_isend error'
                 call flush(nfout)
                 call mpi_abort(mpi_comm_world,163,ierr)
              endif
             icnt_send = icnt_send + 1
          endif
       enddo
       call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
        if (ierr /= 0) then
           write(nfout,*)' map_charge_onto_a_fftcd_box :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world,164,ierr)
        endif
       call mpi_waitall(icnt_send, req_s, sta_s, ierr)
        if (ierr /= 0) then
           write(nfout,*)' map_charge_onto_a_fftcd_box :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world,165,ierr)
        endif
                                                 __TIMER_COMM_STOP(882)
#else
                                                 __TIMER_COMM_START_w_BARRIER(mpi_comm_all,696)
! === DEBUG by tkato 2012/06/05 ================================================
!      integer, allocatable, dimension(:) :: sdsp, rdsp
! ==============================================================================
       allocate(sdsp(0:nrank_g-1), stat=ierr)
       allocate(rdsp(0:nrank_g-1), stat=ierr)
       do i = 0, nrank_g - 1
          sdsp(i)=chgq_fftcd_maxsend*kimg*i
          rdsp(i)=chgq_fftcd_maxrecv*kimg*i
       enddo
       call MPI_ALLTOALLV(      sendbuf, chgq_fftcd_scnt*kimg, sdsp, &
      &   mpi_double_precision, recvbuf, chgq_fftcd_rcnt*kimg, rdsp, &
      &   mpi_double_precision, mpi_comm_all, ierr )
       if (ierr /= 0) then
          write(nfout,*)' map_charge_onto_a_fftcd_box :  mpi_alltoallv error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 166, ierr)
       endif
       deallocate(sdsp)
       deallocate(rdsp)
                                                 __TIMER_COMM_STOP(696)
#endif

       afft_l = 0.0d0
                                                 __TIMER_DO_START(883)
       if (kimg == 1) then
#ifdef FFT_3D_DIVISION_CD
          do i = 0, nrank_g - 1
             do j = 1, chgq_fftcd_rcnt(i)
                if(chgq_fftcd_recv(j,i)<1) cycle
                afft_l(chgq_fftcd_recv(j,i)*2-1,1) = recvbuf(j,i)
                afft_l(chgq_fftcd_recv(j,i)*2  ,1) = 0.0d0
             enddo
          enddo
#else
          do i = 0, nrank_g - 1
             do j = 1, chgq_fftcd_rcnt(i)
                if(chgq_fftcd_recv(j,i)<1) cycle
                afft_l(chgq_fftcd_recv(j,i),1) = recvbuf(j,i)
             enddo
          enddo
#endif
       else
          do i = 0, nrank_g - 1
             do j = 1, chgq_fftcd_rcnt(i)
                if(chgq_fftcd_recv(j,i)<1) cycle
                afft_l(chgq_fftcd_recv(j,i)*2-1,1) = recvbuf(j*2-1,i)
                afft_l(chgq_fftcd_recv(j,i)*2  ,1) = recvbuf(j*2  ,i)
             enddo
          enddo
       endif
                                                 __TIMER_DO_STOP(883)
       if (allocated(sendbuf)) deallocate(sendbuf)
       if (allocated(recvbuf)) deallocate(recvbuf)

                                                  __TIMER_SUB_STOP(753)
    end subroutine map_charge_onto_a_fftcd_box
!-------------------------------------------------------------------------------
    subroutine map_fftcd_onto_charge(is, lsize, ibesize, inn_fft_l, isrsize, fftsize, &
         &                           nfout, mode )
      integer, intent(in)  :: is, lsize, ibesize, isrsize, fftsize, nfout, mode
#ifdef FFT_3D_DIVISION_CD
      real(kind=DP), dimension(lsize*2   ,ibesize), intent(in)    :: inn_fft_l
#else
      real(kind=DP), dimension(lsize*kimg,ibesize), intent(in)    :: inn_fft_l
#endif
#ifdef CD_FFT_ALL
      integer, dimension(0:npes-1)                             ::req_r,req_s
      integer, dimension(MPI_STATUS_SIZE,0:npes-1)             ::sta_r, sta_s
#else
      integer, dimension(0:nrank_g-1)                       ::req_r,req_s
      integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)       ::sta_r, sta_s
#endif
      integer, parameter :: itag = 11

      real(kind=DP), allocatable, dimension(:,:) :: sendbuf, recvbuf
      integer :: icnt_send, icnt_recv, ierr, lrank, i, j, k, iadd
      integer :: nmrank, myrank, mpicom, jrank_e, jrank_g
! === DEBUG by tkato 2012/06/04 ================================================
#ifndef USE_NONBLK_COMM
      integer, allocatable, dimension(:) :: sdsp, rdsp
#endif
! ==============================================================================
                                                  __TIMER_SUB_START(785)

#ifdef CD_FFT_ALL
      nmrank = npes
      myrank = mype
      mpicom = MPI_CommGroup
#else
      nmrank = nrank_g
      myrank = myrank_g
      mpicom = mpi_ke_world
#endif
      ierr = 0
!      if (fftcd_chgq_maxsend /= 0) then
         allocate(sendbuf(fftcd_chgq_maxsend*kimg*ibesize,0:nmrank-1), stat=ierr)
         sendbuf = 0.0d0
!      endif
!      if (fftcd_chgq_maxrecv /= 0) then
         allocate(recvbuf(fftcd_chgq_maxrecv*kimg*ibesize,0:nmrank-1), stat=ierr)
         recvbuf = 0.0d0
!      endif
       if (ierr /= 0) then
          write(nfout,*)' map_fft_onto_charge :  Not allocate '
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 261, ierr)
       endif

#ifdef USE_NONBLK_COMM
                                                 __TIMER_COMM_START_w_BARRIER(mpicom,1021)

      if (fftcd_chgq_maxrecv /= 0) then
         icnt_recv = 0
         lrank = mod(myrank,nmrank)
           do i = 0, nmrank - 1
            lrank = lrank + 1
            if (lrank > (nmrank-1)) lrank = 0
            if (fftcd_chgq_rcnt(lrank) /= 0) then
               call mpi_irecv(recvbuf(1,lrank), fftcd_chgq_rcnt(lrank)*kimg*ibesize, &
              &               mpi_double_precision, lrank, itag, mpicom, req_r(icnt_recv), ierr)
                if (ierr /= 0) then
                   write(nfout,*)' map_fft_onto_charge :  mpi_irecv error'
                   call flush(nfout)
                   call mpi_abort(mpi_comm_world, 262, ierr)
                endif
               icnt_recv = icnt_recv + 1
            endif
         enddo
      endif
#endif

      if (fftcd_chgq_maxsend /= 0) then
                                                 __TIMER_DO_START(1022)
#ifdef FFT_3D_DIVISION_CD
         integer :: i1, lx, ly, lz, mx, my, mz, mm, kx1p, kx2p, kx3p, jadd
         if (kimg == 1) then
           lx = fft_box_size_CD_3D(1,0)
           ly = fft_box_size_CD_3D(2,0)
           lz = fft_box_size_CD_3D(3,0)
           kx1p = fftcd_X_x_nel
           kx2p = fftcd_X_y_nel
           kx3p = fftcd_X_z_nel

            do k = 1, nel_fftcd_x(myrank_g)
               if(fftcd_chgq_index(k) == 0) cycle
               i1 = mp_fftcd_x(k)
               mz = (i1-1)/(lx*ly)+1
               mm = mod(i1,(lx*ly))
               if (mm==0) mm=lx*ly
               my = (mm-1)/lx+1
               mx = mod(mm,lx)
               if (mx==0) mx = lx
               jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
               do i = 1, ibesize
                  sendbuf(ibesize*(fftcd_chgq_index(k)-1)+i,fftcd_chgq_dist(k)) = inn_fft_l(jadd*2-1,i)
               enddo
            end do
         else
           lx = fft_box_size_CD_3D(1,0)
           ly = fft_box_size_CD_3D(2,0)
           lz = fft_box_size_CD_3D(3,0)
           kx1p = fftcd_X_x_nel
           kx2p = fftcd_X_y_nel
           kx3p = fftcd_X_z_nel

            do k = 1, nel_fftcd_x(myrank_g)
               if(fftcd_chgq_index(k) == 0) cycle
               i1 = mp_fftcd_x(k)
               mz = (i1-1)/(lx*ly)+1
               mm = mod(i1,(lx*ly))
               if (mm==0) mm=lx*ly
               my = (mm-1)/lx+1
               mx = mod(mm,lx)
               if (mx==0) mx = lx
               jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
               do i = 1, ibesize
                  iadd = ibesize*2*(fftcd_chgq_index(k)-1)+i*2
                  sendbuf(iadd-1,fftcd_chgq_dist(k)) = inn_fft_l(jadd*2-1,i)
                  sendbuf(iadd,  fftcd_chgq_dist(k)) = inn_fft_l(jadd*2  ,i)
               enddo
            end do
         endif
#else
         if (kimg == 1) then
            do k = 1, nel_fftcd_x(myrank)
               if(fftcd_chgq_index(k) == 0) cycle
               do i = 1, ibesize
                  sendbuf(ibesize*(fftcd_chgq_index(k)-1)+i,fftcd_chgq_dist(k)) = inn_fft_l(k,i)
               enddo
            end do
         else
            do k = 1, nel_fftcd_x(myrank)
               if(fftcd_chgq_index(k) == 0) cycle
               do i = 1, ibesize
                  iadd = ibesize*2*(fftcd_chgq_index(k)-1)+i*2
                  sendbuf(iadd-1,fftcd_chgq_dist(k)) = inn_fft_l(k*2-1,i)
                  sendbuf(iadd,  fftcd_chgq_dist(k)) = inn_fft_l(k*2  ,i)
               enddo
            end do
         end if
#endif
                                                 __TIMER_DO_STOP(1022)

#ifdef USE_NONBLK_COMM
         icnt_send = 0
#ifdef CD_FFT_ALL
         jrank_e = 0
         jrank_g = 0
         do lrank = 0, nmrank - 1
            if (fftcd_chgq_scnt(jrank_g) /= 0) then
               call mpi_isend(sendbuf(1,jrank_g), fftcd_chgq_scnt(jrank_g)*kimg*ibesize, &
              &               mpi_double_precision, lrank, itag, mpicom, req_s(icnt_send), ierr)
                if (ierr /= 0) then
                   write(nfout,*)' map_fft_onto_charge :  mpi_isend error'
                   call flush(nfout)
                   call mpi_abort(mpi_comm_world, 263, ierr)
                endif
               icnt_send = icnt_send + 1
            endif
            jrank_e = jrank_e + 1
            if (jrank_e > (nrank_e-1)) then
               jrank_e = 0
               jrank_g = jrank_g + 1
            end if
         enddo
#else
         lrank = mod((myrank+1),nmrank)
         do i = 0, nmrank - 1
            lrank = lrank + 1
            if (lrank > (nmrank - 1)) lrank = 0
            if (fftcd_chgq_scnt(lrank) /= 0) then
               call mpi_isend(sendbuf(1,lrank), fftcd_chgq_scnt(lrank)*kimg*ibesize, &
              &               mpi_double_precision, lrank, itag, mpicom, req_s(icnt_send), ierr)
                if (ierr /= 0) then
                   write(nfout,*)' map_fft_onto_charge :  mpi_isend error'
                   call flush(nfout)
                   call mpi_abort(mpi_comm_world, 263, ierr)
                endif
               icnt_send = icnt_send + 1
            endif
         enddo
#endif
#endif
      endif

#ifdef USE_NONBLK_COMM
      if (fftcd_chgq_maxrecv /= 0) then
         call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
          if (ierr /= 0) then
             write(nfout,*)' map_fft_onto_charge :  mpi_waitall error'
             call flush(nfout)
             call mpi_abort(mpi_comm_world, 264, ierr)
          endif
      endif

      if (fftcd_chgq_maxsend /= 0) then
         call mpi_waitall(icnt_send, req_s, sta_s, ierr)
          if (ierr /= 0) then
             write(nfout,*)' map_fft_onto_charge :  mpi_waitall error'
             call flush(nfout)
             call mpi_abort(mpi_comm_world, 265, ierr)
          endif
      endif
                                                 __TIMER_COMM_STOP(1021)
#else
                                                 __TIMER_COMM_START_w_BARRIER(mpicom,695)
! === DEBUG by tkato 2012/06/04 ================================================
!     integer, allocatable, dimension(:) :: sdsp, rdsp
! ==============================================================================
      allocate(sdsp(0:nrank_g-1), stat=ierr)
      allocate(rdsp(0:nrank_g-1), stat=ierr)
      do i = 0, nrank_g - 1
         sdsp(i)=fftcd_chgq_maxsend*kimg*ibesize*i
         rdsp(i)=fftcd_chgq_maxrecv*kimg*ibesize*i
      enddo
      call MPI_ALLTOALLV(      sendbuf, fftcd_chgq_scnt*kimg*ibesize, sdsp, &
     &   mpi_double_precision, recvbuf, fftcd_chgq_rcnt*kimg*ibesize, rdsp, &
     &   mpi_double_precision, mpicom, ierr )
      if (ierr /= 0) then
         write(nfout,*)' map_fft_onto_charge :  mpi_alltoallv error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world, 266, ierr)
      endif
      deallocate(sdsp)
      deallocate(rdsp)
                                                 __TIMER_COMM_STOP(695)
#endif


!!      if (input_charge == Partial_Core_Charge) then
      if ( mode == Partial_Core_Charge) then
                                                 __TIMER_DO_START(1023)
         if (kimg == 1) then
            do i = 0, nmrank - 1
               if (fftcd_chgq_rcnt(i) /= 0) then
                  do k = 1, fftcd_chgq_rcnt(i)
                     do j = 1, ibesize
                        vxcpc_l(ista_kngp-1+fftcd_chgq_recv(k,i),1) = recvbuf(ibesize*(k-1)+j,i)*rinplw
                     enddo
                  end do
               end if
            end do
         else
            do i = 0, nmrank - 1
               if (fftcd_chgq_rcnt(i) /= 0) then
                  do k = 1, fftcd_chgq_rcnt(i)
                     do j = 1, ibesize
                        iadd = ibesize*2*(k-1)+j*2
                        vxcpc_l(ista_kngp-1+fftcd_chgq_recv(k,i),1) = recvbuf(iadd-1,i)*rinplw
                        vxcpc_l(ista_kngp-1+fftcd_chgq_recv(k,i),2) = recvbuf(iadd,  i)*rinplw
                     enddo
                  end do
               end if
            end do
         end if
                                                 __TIMER_DO_STOP(1023)
!!      else if (input_charge == Valence_plus_PC_Charge) then
      else if ( mode == Valence_plus_PC_Charge) then
                                                 __TIMER_DO_START(1024)
         if (kimg == 1) then
            do i = 0, nmrank - 1
               if (fftcd_chgq_rcnt(i) /= 0) then
                  do k = 1, fftcd_chgq_rcnt(i)
                     do j = 1, ibesize
                        vxc_l(ista_kngp-1+fftcd_chgq_recv(k,i),1,is) = recvbuf(ibesize*(k-1)+j,i)*rinplw
                     enddo
                  end do
               end if
            end do
         else
            do i = 0, nmrank - 1
               if (fftcd_chgq_rcnt(i) /= 0) then
                  do k = 1, fftcd_chgq_rcnt(i)
                     do j = 1, ibesize
                        iadd = ibesize*2*(k-1)+j*2
                        vxc_l(ista_kngp-1+fftcd_chgq_recv(k,i),1,is) = recvbuf(iadd-1,i)*rinplw
                        vxc_l(ista_kngp-1+fftcd_chgq_recv(k,i),2,is) = recvbuf(iadd,  i)*rinplw
                     enddo
                  end do
               end if
            end do
         end if
                                                 __TIMER_DO_STOP(1024)

      else if ( mode == -1 ) then

         if (kimg == 1) then
            do i = 0, nmrank - 1
               if (fftcd_chgq_rcnt(i) /= 0) then
                  do k = 1, fftcd_chgq_rcnt(i)
                     do j = 1, ibesize
                        vtau_l(ista_kngp-1+fftcd_chgq_recv(k,i),1,is) = recvbuf(ibesize*(k-1)+j,i)*rinplw
                     enddo
                  end do
               end if
            end do
         else
            do i = 0, nmrank - 1
               if (fftcd_chgq_rcnt(i) /= 0) then
                  do k = 1, fftcd_chgq_rcnt(i)
                     do j = 1, ibesize
                        iadd = ibesize*2*(k-1)+j*2
                        vtau_l(ista_kngp-1+fftcd_chgq_recv(k,i),1,is) = recvbuf(iadd-1,i)*rinplw
                        vtau_l(ista_kngp-1+fftcd_chgq_recv(k,i),2,is) = recvbuf(iadd,  i)*rinplw
                     enddo
                  end do
               end if
            end do
         end if

      end if

      if (allocated(sendbuf)) deallocate(sendbuf)
      if (allocated(recvbuf)) deallocate(recvbuf)

                                                 __TIMER_SUB_STOP(785)
    end subroutine map_fftcd_onto_charge
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    subroutine map_fftcd_onto_chargelike_data(is, lsize, ibesize, inn_fft_l, isrsize, fftsize, nfout, chgdata)
      integer, intent(in)  :: is, lsize, ibesize, isrsize, fftsize, nfout
#ifdef FFT_3D_DIVISION_CD
      real(kind=DP), dimension(lsize*2   ,ibesize), intent(in)    :: inn_fft_l
#else
      real(kind=DP), dimension(lsize*kimg,ibesize), intent(in)    :: inn_fft_l
#endif
      real(kind=DP), dimension(ista_kngp:iend_kngp,kimg,nspin), intent(out)  :: chgdata
#ifdef CD_FFT_ALL
      integer, dimension(0:npes-1)                             ::req_r,req_s
      integer, dimension(MPI_STATUS_SIZE,0:npes-1)             ::sta_r, sta_s
#else
      integer, dimension(0:nrank_g-1)                       ::req_r,req_s
      integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)       ::sta_r, sta_s
#endif
      integer, parameter :: itag = 11

      real(kind=DP), allocatable, dimension(:,:) :: sendbuf, recvbuf
      integer :: icnt_send, icnt_recv, ierr, lrank, i, j, k, iadd
      integer :: nmrank, myrank, mpicom, jrank_e, jrank_g
! === DEBUG by tkato 2012/06/04 ================================================
#ifndef USE_NONBLK_COMM
      integer, allocatable, dimension(:) :: sdsp, rdsp
#endif
! ==============================================================================

#ifdef CD_FFT_ALL
      nmrank = npes
      myrank = mype
      mpicom = MPI_CommGroup
#else
      nmrank = nrank_g
      myrank = myrank_g
      mpicom = mpi_ke_world
#endif
      ierr = 0
!      if (fftcd_chgq_maxsend /= 0) then
         allocate(sendbuf(fftcd_chgq_maxsend*kimg*ibesize,0:nmrank-1), stat=ierr)
         sendbuf = 0.0d0
!      endif
!      if (fftcd_chgq_maxrecv /= 0) then
         allocate(recvbuf(fftcd_chgq_maxrecv*kimg*ibesize,0:nmrank-1), stat=ierr)
         recvbuf = 0.0d0
!      endif
       if (ierr /= 0) then
          write(nfout,*)' map_fft_onto_charge :  Not allocate '
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 261, ierr)
       endif

#ifdef USE_NONBLK_COMM

      if (fftcd_chgq_maxrecv /= 0) then
         icnt_recv = 0
         lrank = mod(myrank,nmrank)
           do i = 0, nmrank - 1
            lrank = lrank + 1
            if (lrank > (nmrank-1)) lrank = 0
            if (fftcd_chgq_rcnt(lrank) /= 0) then
               call mpi_irecv(recvbuf(1,lrank), fftcd_chgq_rcnt(lrank)*kimg*ibesize, &
              &               mpi_double_precision, lrank, itag, mpicom, req_r(icnt_recv), ierr)
                if (ierr /= 0) then
                   write(nfout,*)' map_fft_onto_charge :  mpi_irecv error'
                   call flush(nfout)
                   call mpi_abort(mpi_comm_world, 262, ierr)
                endif
               icnt_recv = icnt_recv + 1
            endif
         enddo
      endif
#endif

      if (fftcd_chgq_maxsend /= 0) then
#ifdef FFT_3D_DIVISION_CD
         integer :: i1, lx, ly, lz, mx, my, mz, mm, kx1p, kx2p, kx3p, jadd
         if (kimg == 1) then
           lx = fft_box_size_CD_3D(1,0)
           ly = fft_box_size_CD_3D(2,0)
           lz = fft_box_size_CD_3D(3,0)
           kx1p = fftcd_X_x_nel
           kx2p = fftcd_X_y_nel
           kx3p = fftcd_X_z_nel

            do k = 1, nel_fftcd_x(myrank_g)
               if(fftcd_chgq_index(k) == 0) cycle
               i1 = mp_fftcd_x(k)
               mz = (i1-1)/(lx*ly)+1
               mm = mod(i1,(lx*ly))
               if (mm==0) mm=lx*ly
               my = (mm-1)/lx+1
               mx = mod(mm,lx)
               if (mx==0) mx = lx
               jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
               do i = 1, ibesize
                  sendbuf(ibesize*(fftcd_chgq_index(k)-1)+i,fftcd_chgq_dist(k)) = inn_fft_l(jadd*2-1,i)
               enddo
            end do
         else
           lx = fft_box_size_CD_3D(1,0)
           ly = fft_box_size_CD_3D(2,0)
           lz = fft_box_size_CD_3D(3,0)
           kx1p = fftcd_X_x_nel
           kx2p = fftcd_X_y_nel
           kx3p = fftcd_X_z_nel

            do k = 1, nel_fftcd_x(myrank_g)
               if(fftcd_chgq_index(k) == 0) cycle
               i1 = mp_fftcd_x(k)
               mz = (i1-1)/(lx*ly)+1
               mm = mod(i1,(lx*ly))
               if (mm==0) mm=lx*ly
               my = (mm-1)/lx+1
               mx = mod(mm,lx)
               if (mx==0) mx = lx
               jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
               do i = 1, ibesize
                  iadd = ibesize*2*(fftcd_chgq_index(k)-1)+i*2
                  sendbuf(iadd-1,fftcd_chgq_dist(k)) = inn_fft_l(jadd*2-1,i)
                  sendbuf(iadd,  fftcd_chgq_dist(k)) = inn_fft_l(jadd*2  ,i)
               enddo
            end do
         endif
#else
         if (kimg == 1) then
            do k = 1, nel_fftcd_x(myrank)
               if(fftcd_chgq_index(k) == 0) cycle
               do i = 1, ibesize
                  sendbuf(ibesize*(fftcd_chgq_index(k)-1)+i,fftcd_chgq_dist(k)) = inn_fft_l(k,i)
               enddo
            end do
         else
            do k = 1, nel_fftcd_x(myrank)
               if(fftcd_chgq_index(k) == 0) cycle
               do i = 1, ibesize
                  iadd = ibesize*2*(fftcd_chgq_index(k)-1)+i*2
                  sendbuf(iadd-1,fftcd_chgq_dist(k)) = inn_fft_l(k*2-1,i)
                  sendbuf(iadd,  fftcd_chgq_dist(k)) = inn_fft_l(k*2  ,i)
               enddo
            end do
         end if
#endif

#ifdef USE_NONBLK_COMM
         icnt_send = 0
#ifdef CD_FFT_ALL
         jrank_e = 0
         jrank_g = 0
         do lrank = 0, nmrank - 1
            if (fftcd_chgq_scnt(jrank_g) /= 0) then
               call mpi_isend(sendbuf(1,jrank_g), fftcd_chgq_scnt(jrank_g)*kimg*ibesize, &
              &               mpi_double_precision, lrank, itag, mpicom, req_s(icnt_send), ierr)
                if (ierr /= 0) then
                   write(nfout,*)' map_fft_onto_charge :  mpi_isend error'
                   call flush(nfout)
                   call mpi_abort(mpi_comm_world, 263, ierr)
                endif
               icnt_send = icnt_send + 1
            endif
            jrank_e = jrank_e + 1
            if (jrank_e > (nrank_e-1)) then
               jrank_e = 0
               jrank_g = jrank_g + 1
            end if
         enddo
#else
         lrank = mod((myrank+1),nmrank)
         do i = 0, nmrank - 1
            lrank = lrank + 1
            if (lrank > (nmrank - 1)) lrank = 0
            if (fftcd_chgq_scnt(lrank) /= 0) then
               call mpi_isend(sendbuf(1,lrank), fftcd_chgq_scnt(lrank)*kimg*ibesize, &
              &               mpi_double_precision, lrank, itag, mpicom, req_s(icnt_send), ierr)
                if (ierr /= 0) then
                   write(nfout,*)' map_fft_onto_charge :  mpi_isend error'
                   call flush(nfout)
                   call mpi_abort(mpi_comm_world, 263, ierr)
                endif
               icnt_send = icnt_send + 1
            endif
         enddo
#endif
#endif
      endif

#ifdef USE_NONBLK_COMM
      if (fftcd_chgq_maxrecv /= 0) then
         call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
          if (ierr /= 0) then
             write(nfout,*)' map_fft_onto_charge :  mpi_waitall error'
             call flush(nfout)
             call mpi_abort(mpi_comm_world, 264, ierr)
          endif
      endif

      if (fftcd_chgq_maxsend /= 0) then
         call mpi_waitall(icnt_send, req_s, sta_s, ierr)
          if (ierr /= 0) then
             write(nfout,*)' map_fft_onto_charge :  mpi_waitall error'
             call flush(nfout)
             call mpi_abort(mpi_comm_world, 265, ierr)
          endif
      endif
#else
! === DEBUG by tkato 2012/06/04 ================================================
!     integer, allocatable, dimension(:) :: sdsp, rdsp
! ==============================================================================
      allocate(sdsp(0:nrank_g-1), stat=ierr)
      allocate(rdsp(0:nrank_g-1), stat=ierr)
      do i = 0, nrank_g - 1
         sdsp(i)=fftcd_chgq_maxsend*kimg*ibesize*i
         rdsp(i)=fftcd_chgq_maxrecv*kimg*ibesize*i
      enddo
      call MPI_ALLTOALLV(      sendbuf, fftcd_chgq_scnt*kimg*ibesize, sdsp, &
     &   mpi_double_precision, recvbuf, fftcd_chgq_rcnt*kimg*ibesize, rdsp, &
     &   mpi_double_precision, mpicom, ierr )
      if (ierr /= 0) then
         write(nfout,*)' map_fft_onto_charge :  mpi_alltoallv error'
         call flush(nfout)
         call mpi_abort(mpi_comm_world, 266, ierr)
      endif
      deallocate(sdsp)
      deallocate(rdsp)
#endif

         if (kimg == 1) then
            do i = 0, nmrank - 1
               if (fftcd_chgq_rcnt(i) /= 0) then
                  do k = 1, fftcd_chgq_rcnt(i)
                     do j = 1, ibesize
                        chgdata(ista_kngp-1+fftcd_chgq_recv(k,i),1,is) = recvbuf(ibesize*(k-1)+j,i)*rinplw
                     enddo
                  end do
               end if
            end do
         else
            do i = 0, nmrank - 1
               if (fftcd_chgq_rcnt(i) /= 0) then
                  do k = 1, fftcd_chgq_rcnt(i)
                     do j = 1, ibesize
                        iadd = ibesize*2*(k-1)+j*2
                        chgdata(ista_kngp-1+fftcd_chgq_recv(k,i),1,is) = recvbuf(iadd-1,i)*rinplw
                        chgdata(ista_kngp-1+fftcd_chgq_recv(k,i),2,is) = recvbuf(iadd,  i)*rinplw
                     enddo
                  end do
               end if
            end do
         end if

      if (allocated(sendbuf)) deallocate(sendbuf)
      if (allocated(recvbuf)) deallocate(recvbuf)

    end subroutine map_fftcd_onto_chargelike_data

    subroutine set_ekindens_case_collinear_3D( ispin, ekin_l, ekin_dens )
      integer, intent(in) :: ispin
      real(kind=DP), intent(in)  :: ekin_l(ista_kngp:iend_kngp,kimg,nspin)
      real(kind=DP), intent(out) :: ekin_dens(1:nfft_div_size,nspin)

      integer :: iloop, i
      integer :: lx, ly, lz, i1, mx, my, mz, mm, jadd, kx1p, kx2p, kx3p

      do iloop = 1, ispin
         call map_charge_onto_a_fftcd_box( ekin_l, nspin, iloop, Valence_Charge_Only )

#ifdef FFT_3D_DIVISION_CD
         call m_FFT_CD_Inverse_3DIV_3D(nfout,afft_l,lsize,1)
#else
         if (sw_fft_xzy > 0) then
            call m_FFT_CD_Inverse_3D(nfout,afft_l,lsize,1)
         else
            call m_FFT_CD_Inverse_XYZ_3D(nfout,afft_l,lsize,1)
         end if
#endif

#ifdef FFT_3D_DIVISION_CD
         kx1p = fftcd_X_x_nel
         kx2p = fftcd_X_y_nel
         kx3p = fftcd_X_z_nel
         lx = fft_box_size_CD_3D(1,0)
         ly = fft_box_size_CD_3D(2,0)
         lz = fft_box_size_CD_3D(3,0)
         do i = 1, np_fftcd_x
            i1 = mp_fftcd_x(i)
            mz = (i1-1)/(lx*ly)+1
            mm = mod(i1,(lx*ly))
            if (mm==0) mm=lx*ly
            my = (mm-1)/lx+1
            mx = mod(mm,lx)
            if (mx==0) mx = lx
            jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2)) &
                 & +kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
            ekin_dens(i,iloop) = afft_l(jadd*2-1,1)
         end do
#else
         do i = 1, nfft_div_size
            ekin_dens(i,iloop) = afft_l(i*2-1,1)
         end do
#endif
      end do
    end subroutine set_ekindens_case_collinear_3D

    subroutine ggaxcp_3D()
      !        ~~~~~~
      integer       :: is, id_sname = -1
      real(kind=DP) :: exc_mpi
                                                 __TIMER_SUB_START(771)
      call tstatc0_begin('ggaxcp(in xc_pot.) ',id_sname)

#ifdef LIBXC
      if ( xctype == "libxc" ) then
         call ggaxcp0_3D_libxc()
      else
         call ggaxcp0_3D()
      endif
#else
      call ggaxcp0_3D()
#endif
      !    ~~~~~~~
             !   dF/d|rho(r)| (=vxc) --> dF_drho
             !   dFx/d|grad(rho(r))| --> dF_dgradrho
             !   dFc/d|grad(rho(r))| --> grad_trho
             !   ->exc

      if(vflag == VXC_AND_EXC .or. vflag == STRESS_) then
         if(xctype /= 'ldapw91' .and. xctype /= 'ldapbe ') then
            call dFxc_over_ddgradrho_3D()
            !   sum(q=xyz)(dFxc/dq(d|grad(rho)|)) --> grad_rho
         end if

         if ( use_metagga ) then
            if ( xctype == "tb09" ) then
               grad_rho = 0.0d0
            else
               call calc_lapl_dF_dlapl_3D;
            endif
#ifdef LIBXC
            if ( xc_name_exch == "mgga_x_tb09" ) then
               grad_rho = 0.0d0;    dF_dlaplrho = 0.0d0
            endif
#endif
         endif

         do is = 1, ispin
            call finally_gga_xc_pot_3D(is)
            ! dFxc/d(rho)  dF_drho + grad_rho --> chgrhr
         end do
      end if

      exc = exc*univol*rinplw

#ifdef CD_FFT_ALL
      if(npes >= 2) then
                                                 __TIMER_COMM_START_w_BARRIER(MPI_CommGroup,1005)
         call mpi_allreduce(exc,exc_mpi,1,mpi_double_precision,mpi_sum, MPI_CommGroup,ierr)
                                                 __TIMER_COMM_STOP(1005)
         exc = exc_mpi
      end if
#else
      if(nrank_g >= 2) then
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,1006)
         call mpi_allreduce(exc,exc_mpi,1,mpi_double_precision,mpi_sum, mpi_ke_world,ierr)
                                                 __TIMER_COMM_STOP(1006)
         exc = exc_mpi
      end if
#endif
      call tstatc0_end(id_sname)
                                                 __TIMER_SUB_STOP(771)
    end subroutine ggaxcp_3D

    subroutine ggaxcp_diff_3D
      !        ~~~~~~~~~~~
      real(DP):: exc_mpi               ! MPI
      integer :: is,i1,i2,in
      integer :: nfftx, nffty, nfftz
      real(kind=DP),allocatable,dimension(:) :: chgrhr_in
      real(kind=DP),allocatable,dimension(:) :: grad_rho_in
      real(kind=DP),allocatable,dimension(:,:) :: cgrad_rho_in
                                                 __TIMER_SUB_START(757)

#ifdef LIBXC
      if ( xctype == "libxc" ) then
         call ggaxcp0_3D_libxc()
      else
         call ggaxcp0_3D()
      endif
#else
      call ggaxcp0_3D()
#endif
      !    ~~~~~~~
             !   dF/d|rho(r)| (=vxc) --> dF_drho
             !   dFx/d|grad(rho(r))| --> dF_dgradrho
             !   dFc/d|grad(rho(r))| --> grad_trho
             !   ->exc

      exc = exc*univol*rinplw

#ifdef CD_FFT_ALL
      if(npes >= 2) then
                                                 __TIMER_COMM_START_w_BARRIER(MPI_CommGroup,890)
         call mpi_allreduce(exc,exc_mpi,1,mpi_double_precision,mpi_sum, MPI_CommGroup,ierr)
                                                 __TIMER_COMM_STOP(890)
         exc = exc_mpi
      end if
#else
      if(nrank_g >= 2) then
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,891)
         call mpi_allreduce(exc,exc_mpi,1,mpi_double_precision,mpi_sum, mpi_ke_world,ierr)
                                                 __TIMER_COMM_STOP(891)
         exc = exc_mpi
      end if
#endif

      if ( xctype .eq. 'vdwdf' .and. ( .not. oneshot ) ) then
         nfftx = fft_box_size_CD(1,1);
         nffty = fft_box_size_CD(2,1);
         nfftz = fft_box_size_CD(3,1)
         allocate(chgrhr_in(nfft_div_size))
         allocate(grad_rho_in(nfft_div_size))
         allocate(cgrad_rho_in(nfft_div_size,3))
         if (ispin == 1) then
            chgrhr_in(:)   = chgrhr_l(:,1)
            grad_rho_in(:) = grad_rho(:,1)
            cgrad_rho_in(:,1:3) = cgrad_rho(:,1:3,1)
         else
            chgrhr_in(:)   = chgrhr_l(:,1)+chgrhr_l(:,2)
            grad_rho_in(:) = grad_rho(:,1)+grad_rho(:,2)
            cgrad_rho_in(:,1:3) = cgrad_rho(:,1:3,1)+cgrad_rho(:,1:3,2)
         endif
         call vdW_scf_stress( nspin, ispin, nfft_div_size, nfftx, nffty, nfftz, &
              &               chgrhr_in, grad_rho_in, cgrad_rho_in, &
              &               vdwdf_version, input_charge )
         deallocate(chgrhr_in);deallocate(grad_rho_in);deallocate(cgrad_rho_in)
      endif

      chgrhr_l = dF_drho

! stress tensor comes from gradient correction
!
      call allocate_for_stress_3D
      do i1=1,3
      do i2=1,3
        call rhos_diff(i1,i2)
        call rhopc_diff_3D(i1,i2)
        call rhoh_diff_3D(i1,i2)
        do is=1,nspin-af
           call map_drhodh1_3D(is) ! drhodh -> chgfft
           do in=1,3
              grad_rho(1:nfft_div_size,is)=0.d0
              dF_drho(1:nfft_div_size,is)=0.d0
              call dgrhodh1_3D(i1,i2,in,is)
              call stress_exchange_part_3D(i1,i2,in,is)   ! ->(s_gga1,s_gga2)
           enddo !in
        enddo !is
        if(nspin==2) then
           if(af == 1) then
              is = 1
              call map_drhodh1_3D(is) ! drhodh -> chgfft
              do in = 1, 3
                 call dgrhodh1_3D(i1,i2,in,is)
                 call stress_correlation_part_3D(i1,i2,in)   ! ->(s_gga2)
              end do
           else
              call map_drhodh2_3D
              do in=1,3
                 call dgrhodh2_3D(i1,i2,in)
                 call stress_correlation_part_3D(i1,i2,in)   ! ->(s_gga2)
              end do
           end if
        end if
      enddo
      enddo

      call sum_s_gga12
      call deallocate_for_stress
                                                 __TIMER_SUB_STOP(757)
    end subroutine ggaxcp_diff_3D

    subroutine ggaxcp0_3D
      integer :: i
      integer :: pot_type
      real(kind=DP), allocatable :: dummy1(:,:), dummy2(:)

                                                 __TIMER_SUB_START(758)
      call initialize_cggawk_arrays_3D

      if(xctype == 'ldapw91' .or. xctype == 'ldapbe ') then
         grad_rho = 0.d0; grad_trho = 0.d0
         cgrad_rho = 0.d0
      else
         call abs_grad_rho_up_down_total_3D()
      !   chden_l ->
      !   |grad(rho_up(r))|,|grad(rho_down(r))|  --> grad_rho
      !   |grad(rho_up(r) + rho_down(r))|        --> grad_trho
      end if

      if ( use_metagga ) then
         call grad2_rho_up_down_3D()
         if ( use_modeled_ekin_density ) then
            call m_KE_set_modeled_ekin_density( ista_fftph, iend_fftph, &
                 &                              chgrhr_l, grad_rho, grad2_rho, &
                 &                              ekin_dens )
         else if ( use_asymm_ekin_density ) then
            ekin_dens = ekin_dens + 0.25d0 *grad2_rho       ! ????
         endif
      endif

      if ( use_metagga ) then
         if ( xctype /= 'tb09' ) then
            write(nfout,'(" xctype = ",a7)') xctype
            call phase_error_with_msg(nfout,' xctype is not set properly (ggaxcp0)',__LINE__,__FILE__)
         endif
      else
         if (xctype /= 'ggapw91' .and. xctype /= 'ldapw91' &
                 & .and. xctype /= 'ggapbe' .and. xctype /= 'ldapbe' &
                 & .and. xctype /= 'katopbe' .and. xctype /= 'ggapbek' &
                 & .and. xctype /= 'rpbe' .and. xctype /= 'wc06' &
                 & .and. xctype /= 'htbs' .and. xctype /= 'pbesol' &
                 & .and. xctype /= 'pbeint' &
                 & .and. xctype /= 'ev93' .and. xctype /= 'evpw91' &
                 & .and. xctype /= 'revpbe' &
                 & .and. xctype /= 'ggapbey' &
                 & .and. xctype /= 'vdwdf' ) then
            write(nfout,'(" xctype = ",a7)') xctype
            call phase_error_with_msg(nfout,' xctype is not set properly (ggaxcp0)',__LINE__,__FILE__)
         end if
      endif

#ifdef __EDA__
!      if(sw_eda==ON) then
!      allocate(exc_on_a_grid_wk(nfft_div_size))
!      allocate(exc_on_a_grid_wk(1))
!      else
      allocate(exc_on_a_grid_wk(1))
!      endif
#endif

      if(xctype == 'ggapw91' .or. xctype == 'ldapw91') then
         call ex_ggapw91_3D(nspin,ispin,chgrhr_l,grad_rho,f2or1,exc,dF_drho,dF_dgradrho,nfft_div_size)
! === Support Hybrid on 3D_Parallel by tkato 2013/02/10 ========================
         if(sw_hybrid_functional==ON) call scale_exchange(alpha_exx)
! ==============================================================================
         call cr_ggapw91_3D(nspin,ispin,chgrhr_l,grad_trho,f2or1,exc,dF_drho,nfft_div_size)
      else if(xctype == 'ggapbe ' .or. xctype == 'ldapbe '.or.xctype=='revpbe') then
         call ex_ggapbe_3D (nspin,ispin,chgrhr_l,grad_rho,f2or1,exc,dF_drho,dF_dgradrho,nfft_div_size,iteration,xctype=='revpbe')
! === Support Hybrid on 3D_Parallel by tkato 2013/02/10 ========================
         if(sw_hybrid_functional==ON) then
            if ( hybrid_functional_type == 'gaupbe' ) then
               call ex_gaupbe( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
              &               grad_rho, f2or1, exc, dF_drho, dF_dgradrho )
            else if ( hybrid_functional_type == 'hse06-hjs' ) then
               call ex_omegapbe_library_HJS( -alpha_exx, omega_exx_pbe, &
                    &                        nspin, ispin, 1, nfft_div_size, &
                    &                        chgrhr_l, grad_rho, f2or1, exc, &
                    &                        dF_drho, dF_dgradrho, 1 )
            else
               if(sw_screened_exchange==ON) then
                  call ex_omegapbe(-alpha_exx,omega_exx_pbe,nspin,ispin,ista_fftph,iend_fftph,chgrhr_l, &
                       &      grad_rho,f2or1,exc,dF_drho,dF_dgradrho,nfft_div_size)
               else
                  call scale_exchange(alpha_exx)
               end if
            endif
         end if
! ==============================================================================
         call cr_ggapbe_3D (nspin,ispin,chgrhr_l,grad_trho,f2or1,exc,dF_drho,nfft_div_size)
! =========================== KT_add ========================= 13.0A
      else if ( xctype == 'rpbe' .or. xctype == 'wc06' .or. xctype == 'htbs' &
           &                     .or. xctype == 'pbesol' .or. xctype == 'pbeint' &
           &                     .or. xctype == 'ev93' .or. xctype == 'evpw91' ) then
         if ( xctype == 'ggapbe' )  pot_type = 1
         if ( xctype == 'revpbe' )  pot_type = 2
         if ( xctype == 'rpbe' )    pot_type = 3
         if ( xctype == 'wc06' )    pot_type = 4
         if ( xctype == 'htbs' )    pot_type = 5

         if ( xctype == 'pbesol' )    pot_type = 6
         if ( xctype == 'pbeint' )    pot_type = 7

         if ( xctype == 'ev93' )    pot_type = 20
         if ( xctype == 'evpw91' )  pot_type = 20

#ifdef __EDA__
         call ex_gga_library( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
              &               grad_rho, f2or1, exc, dF_drho, dF_dgradrho, &
              &               exc_on_a_grid_wk, pot_type, 1, nfft_div_size )
#else
         call ex_gga_library( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
              &               grad_rho, f2or1, exc, dF_drho, dF_dgradrho, pot_type, &
              &               1, nfft_div_size )
#endif

         eex = exc
         if ( sw_hybrid_functional == ON ) then
            if ( hybrid_functional_type == 'hsesol-hjs' ) then
               call ex_omegapbe_library_HJS( -alpha_exx, omega_exx_pbe, &
                    &                        nspin, ispin, 1, nfft_div_size, &
                    &                        chgrhr_l, grad_rho, f2or1, exc, &
                    &                        dF_drho, dF_dgradrho, 6 )
            else
               if ( sw_screened_exchange == ON ) then
                  call ex_omegapbe( -alpha_exx, omega_exx_pbe, nspin, ispin, &
                       &            ista_fftph, iend_fftph, chgrhr_l, &
                       &            grad_rho, f2or1, exc, dF_drho, dF_dgradrho,nfft_div_size )
               else
                  call scale_exchange( alpha_exx )
               end if
            endif
         end if
         if ( sw_exchange_only == OFF ) then
            if ( xctype == 'ev93' ) then
               grad_trho=0.d0
#ifdef __EDA__
               call cr_gga_library( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
                    &               grad_trho, f2or1, exc, dF_drho, exc_on_a_grid_wk, &
                    &               ecor, pot_type, 1, nfft_div_size )
#else
               call cr_gga_library( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
                    &               grad_trho, f2or1, exc, dF_drho, &
                    &               ecor, pot_type, 1, nfft_div_size )
#endif
            else if ( xctype == 'evpw91' ) then
#ifdef __EDA__
               call cr_ggapw91( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
                   &            grad_trho, f2or1, exc, dF_drho, exc_on_a_grid_wk, &
                   &            1, nfft_div_size )
#else
               call cr_ggapw91( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
                   &            grad_trho, f2or1, exc, dF_drho, 1, nfft_div_size )
#endif
            else
#ifdef __EDA__
               call cr_gga_library( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
                    &               grad_trho, f2or1, exc, dF_drho, exc_on_a_grid_wk, &
                    &               ecor, pot_type, 1, nfft_div_size )
#else
               call cr_gga_library( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
                    &               grad_trho, f2or1, exc, dF_drho, ecor, pot_type, &
                    &               1, nfft_div_size )
#endif
            endif
         endif
! ============================================================= 13.0A
      else if(xctype == 'vdwdf')then
         if ( exchange_pot_type == 'pbe' ) then

           call ex_ggapbe_3D (nspin,ispin,chgrhr_l,grad_rho,f2or1,exc,dF_drho,dF_dgradrho, &
         & nfft_div_size,iteration,.false.)
         else if (exchange_pot_type == 'revpbe' ) then
           call ex_ggapbe_3D (nspin,ispin,chgrhr_l,grad_rho,f2or1,exc,dF_drho,dF_dgradrho, &
         & nfft_div_size,iteration,.true.)
         else
            if ( exchange_pot_type == 'b86r' )  pot_type = 11
            if ( exchange_pot_type == 'optpbe' )  pot_type = 12
            if ( exchange_pot_type == 'optb86b' )  pot_type = 13
            if ( exchange_pot_type == 'pw86r' )  pot_type = 14
            if ( exchange_pot_type == 'c09x' )  pot_type = 15
            if ( exchange_pot_type == 'lvpw86r' )  pot_type = 16

#ifdef __EDA__
            call ex_gga_library( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
                 &               grad_rho, f2or1, exc, dF_drho, dF_dgradrho, &
                 &               exc_on_a_grid_wk, pot_type, 1, nfft_div_size )
#else
            call ex_gga_library( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
                 &               grad_rho, f2or1, exc, dF_drho, dF_dgradrho, pot_type, &
                 &               1, nfft_div_size )
#endif
         endif
! === Support Hybrid on 3D_Parallel by tkato 2013/02/10 ========================
         if(sw_hybrid_functional==ON) then
            if(sw_screened_exchange==ON) then
               call ex_omegapbe(-alpha_exx,omega_exx_pbe,nspin,ispin,ista_fftph,iend_fftph,chgrhr_l, &
             &      grad_rho,f2or1,exc,dF_drho,dF_dgradrho,nfft_div_size)
            else
               call scale_exchange(alpha_exx)
            end if
         end if
! ==============================================================================
         if(sw_exchange_only==OFF) then
           grad_trho=0.d0
           call cr_ggapbe_3D (nspin,ispin,chgrhr_l,grad_trho,f2or1,exc,dF_drho,nfft_div_size)
         endif
         if (.not. oneshot ) call add_vdwdf_nonlocal_energy( vdwdf_version )
      else if(xctype == 'katopbe' .or. xctype == 'ggapbek' ) then
         call ex_ggapbe_3D (nspin,ispin,chgrhr_l,grad_rho,f2or1,exc,dF_drho,dF_dgradrho,nfft_div_size,iteration)
! === Support Hybrid on 3D_Parallel by tkato 2013/02/10 ========================
         if(sw_hybrid_functional==ON) then
            if(sw_screened_exchange==ON) then
               call ex_omegapbe(-alpha_exx,omega_exx_pbe,nspin,ispin,ista_fftph,iend_fftph,chgrhr_l, &
             & grad_rho,f2or1,exc,dF_drho,dF_dgradrho,nfft_div_size)
            else
               call scale_exchange(alpha_exx)
            end if
         end if
! ==============================================================================
         call cr_ggapbe_3D (nspin,ispin,chgrhr_l,grad_trho,f2or1,exc,dF_drho,nfft_div_size)
      else if(xctype == 'ggabp  ') then
         call xclda_3D(nspin,ispin,chgrhr_l,f2or1,exc,dF_drho,nfft_div_size)
         call ggabek_3D(nspin,ispin,chgrhr_l,grad_rho,f2or1,exc,dF_drho,dF_dgradrho,nfft_div_size)
         call ggaprd_3D(nspin,ispin,chgrhr_l,grad_rho,f2or1,exc,dF_drho,dF_dgradrho,nfft_div_size)
      end if
             !   dF/d|rho(r)| (=vxc) --> dF_drho
             !   dFx/d|grad(rho(r))| --> dF_dgradrho
             !   dFc/d|grad(rho(r))| --> grad_trho


      if ( use_metagga ) then
         if ( xctype == "tb09" ) then
            if ( sw_fix_val_c_tb09 == OFF ) then
               call set_gval_tb09( nspin, 1, nfft_div_size, &
                    &              chgrhr_l, grad_rho, f2or1, val_g )
               val_g = val_g +val_g_tb09_paw
               call set_cval_tb09( val_g, val_c_tb09 )
            endif
            write(nfout,*) "val_c_tb09 is now ", val_c_tb09
         endif

         if ( xctype == "tb09" ) then
            allocate( dummy1(1:nfft_div_size,nspin ) ); dummy1 = 0.0d0
            allocate( dummy2(1:nfft_div_size) );       dummy2 = 0.0d0

            call ex_ggapw91_3D( nspin, ispin, chgrhr_l, dummy1, f2or1, &
                 &              exc, dF_drho, dF_dgradrho, nfft_div_size )
            dF_drho = 0.0d0

            if ( ekin_density_is_active ) then
               call ex_mgga_tb09( nspin, ispin, 1, nfft_div_size, 1, &
                    &              chgrhr_l, grad_rho, grad2_rho, ekin_dens, &
                    &              dF_drho, val_c_tb09 )
            endif
            call cr_ggapw91_3D( nspin, ispin, chgrhr_l, dummy2, f2or1, &
                 &             exc, dF_drho, nfft_div_size )

            deallocate( dummy1 ); deallocate( dummy2 )
            dF_dgradrho = 0.0d0
         endif
      end if
#ifdef __EDA__
!      if(sw_eda==ON) then
!      call cp_exc_for_EDA
!      endif
      deallocate(exc_on_a_grid_wk)
#endif

      if(iprixc >= 2) then
         write(nfout,'(" !XC -- chgrhr_l -- <<ggaxcp0>>")')
         write(nfout,'(" !XC ",6f12.6)') (chgrhr_l(i,1),i=1,18)
         write(nfout,'(" !XC -- grad_rho -- <<ggaxcp0>>")')
         write(nfout,'(" !XC ",6f12.6)') (grad_rho(i,1),i=1,18)
         write(nfout,'(" !XC -- grad_trho -- <<ggaxcp0>>")')
         write(nfout,'(" !XC ",6f12.6)') (grad_trho(i),i=1,18)
         write(nfout,'(" !XC -- dF_drho -- <<ggaxcp0>>")')
         write(nfout,'(" !XC ",6f12.6)') (dF_drho(i,1),i=1,18)
         write(nfout,'(" !XC -- dF_dgradrho -- <<ggaxcp0>>")')
         write(nfout,'(" !XC ",6f12.6)') (dF_dgradrho(i,1),i=1,18)
      end if
                                                 __TIMER_SUB_STOP(758)
    end subroutine ggaxcp0_3D

#ifdef LIBXC
    subroutine ggaxcp0_3D_libxc
      integer :: i, key_save

      call initialize_cggawk_arrays_3D
      exc = 0.0d0

      if ( xc_family_exch == XC_FAMILY_LDA .and. xc_family_corr == XC_FAMILY_LDA  ) then
         grad_rho = 0.d0; grad_trho = 0.d0
         cgrad_rho = 0.d0
      else
         call abs_grad_rho_up_down_total_3D()
      !   chden_l ->
      !   |grad(rho_up(r))|,|grad(rho_down(r))|  --> grad_rho
      !   |grad(rho_up(r) + rho_down(r))|        --> grad_trho
      end if

! ====== KT_add ======= 13.0XX
      if ( use_metagga ) then
         call grad2_rho_up_down_3D()

!         if ( use_modeled_ekin_density ) then
         if ( use_modeled_ekin_density .or. ( icond==0 .and. iteration == 0 ) ) then
            if ( iteration == 0 ) then
               key_save = ekin_density_type;  ekin_density_type = 3
            endif
            call m_KE_set_modeled_ekin_density( 1, nfft_div_size, &
                 &                              chgrhr_l, grad_rho, grad2_rho, &
                 &                              ekin_dens )
            if ( iteration == 0 ) ekin_density_type = key_save
         else if ( use_asymm_ekin_density ) then
            ekin_dens = ekin_dens + 0.25d0 *grad2_rho       ! ????
         endif

      endif

#ifdef __EDA__
      allocate(exc_on_a_grid_wk(1))
#endif

      select case (xc_family_exch)
      case (XC_FAMILY_LDA)
#ifdef __EDA__
         call ex_lda_libxc( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
              &             f2or1, exc, dF_drho, exc_on_a_grid_wk, 1, nfft_div_size )
#else
         call ex_lda_libxc( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
              &             f2or1, exc, dF_drho, 1, nfft_div_size )
#endif
      case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
#ifdef __EDA__
         call ex_gga_libxc( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
              &             grad_rho, grad_trho, f2or1, exc, dF_drho, dF_dgradrho, &
              &             exc_on_a_grid_wk, 1, nfft_div_size )
#else
         call ex_gga_libxc( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
              &             grad_rho, grad_trho, f2or1, exc, dF_drho, dF_dgradrho, &
              &             1, nfft_div_size )
#endif
      case (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
         if ( xc_name_exch == "mgga_x_tb09" ) then
            if ( sw_fix_val_c_tb09 == OFF ) then
               call set_gval_tb09( nspin, 1, nfft_div_size, &
                    &              chgrhr_l, grad_rho, f2or1, val_g )
               val_g = val_g +val_g_tb09_paw
               call set_cval_tb09( val_g, val_c_tb09 )
            endif
            write(nfout,*) "val_c_tb09 is now ", val_c_tb09
         endif
         call ex_mgga_libxc( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
              &              grad_rho, grad_trho, grad2_rho, ekin_dens, f2or1, exc, &
              &              dF_drho, dF_dgradrho, dF_dlaplrho, dF_dtau )
      end select

      select case (xc_family_corr)
      case (XC_FAMILY_LDA)
#ifdef __EDA__
         call cr_lda_libxc( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
              &             f2or1, exc, dF_drho, &
              &             exc_on_a_grid_wk, 1, nfft_div_size )
#else
         call cr_lda_libxc( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
              &             f2or1, exc, dF_drho, 1, nfft_div_size )
#endif
      case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
#ifdef __EDA__
         call cr_gga_libxc( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
              &             grad_rho, grad_trho, f2or1, exc, dF_drho, grad_trho, &
              &             exc_on_a_grid_wk, 1, nfft_div_size )
#else
         call cr_gga_libxc( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
              &             grad_rho, grad_trho, f2or1, exc, dF_drho, grad_trho, &
              &             1, nfft_div_size )
#endif
      case (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
         call cr_mgga_libxc( nspin, ispin, 1, nfft_div_size, chgrhr_l, &
              &              grad_rho, grad_trho, grad2_rho, ekin_dens, f2or1, exc, &
              &              dF_drho, grad_trho, dF_dlaplrho, dF_dtau )
      end select
#ifdef __EDA__
      deallocate(exc_on_a_grid_wk)
#endif

    end subroutine ggaxcp0_3D_libxc
#endif

    subroutine calc_lapl_dF_dlapl_3D
      integer :: iloop, n, i
      real(kind=DP) :: gx, gy, gz, g2
      real(kind=DP), allocatable :: afft_l(:,:)
#ifdef FFT_3D_DIVISION_CD
      integer :: lx, ly, lz, n, ixyz, mx, my, mz, mm, jadd, kx1p, kx2p, kx3p
#endif

#ifdef FFT_3D_DIVISION_CD
      allocate(afft_l(lsize*2,1))
#else
      allocate(afft_l(lsize*kimg,1))
#endif

      Do iloop=1, nspin
         call cp_vals_to_afft_rspace_3D( iloop, dF_dlaplrho, afft_l )
#ifdef FFT_3D_DIVISION_CD
         call m_FFT_CD_Direct_3DIV_3D(nfout,afft_l,lsize,1)
#else
         if (sw_fft_xzy > 0) then
            call m_FFT_CD_Direct_3D(nfout,afft_l,lsize,1)
         else
            call m_FFT_CD_Direct_XYZ_3D(nfout,afft_l,lsize,1)
         end if
#endif

#ifdef FFT_3D_DIVISION_CD
         if (kimg == 1) then
            kx1p = fftcd_X_x_nel
            kx2p = fftcd_X_y_nel
            kx3p = fftcd_X_z_nel
            lx = fft_box_size_CD_3D(1,0)
            ly = fft_box_size_CD_3D(2,0)
            lz = fft_box_size_CD_3D(3,0)
            do n = 1, np_fftcd_x
               ixyz = mp_fftcd_x(n)
               mz = (ixyz-1)/(lx*ly)+1
               mm = mod(ixyz,(lx*ly))
               if (mm==0) mm=lx*ly
               my = (mm-1)/lx+1
               mx = mod(mm,lx)
               if (mx==0) mx = lx
               jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
               gx = rltv(1,1)*inx(n    )+rltv(1,2)*jnx(n    )+rltv(1,3)*knx(n    )
               gy = rltv(2,1)*inx(n    )+rltv(2,2)*jnx(n    )+rltv(2,3)*knx(n    )
               gz = rltv(3,1)*inx(n    )+rltv(3,2)*jnx(n    )+rltv(3,3)*knx(n    )
               g2 = gx**2 +gy**2 +gz**2
               afft_l(jadd*2-1,1) = -g2 *afft_l(jadd*2-1,1)

               gx = rltv(1,1)*inx(n    )+rltv(1,2)*jnx(n    )+rltv(1,3)*knx(n    )
               gy = rltv(2,1)*inx(n    )+rltv(2,2)*jnx(n    )+rltv(2,3)*knx(n    )
               gz = rltv(3,1)*inx(n    )+rltv(3,2)*jnx(n    )+rltv(3,3)*knx(n    )
               g2 = gx**2 +gy**2 +gz**2
               afft_l(jadd*2  ,1) = -g2 *afft_l(jadd*2,1)
            end do
         else
            kx1p = fftcd_X_x_nel
            kx2p = fftcd_X_y_nel
            kx3p = fftcd_X_z_nel
            lx = fft_box_size_CD_3D(1,0)
            ly = fft_box_size_CD_3D(2,0)
            lz = fft_box_size_CD_3D(3,0)
            do n = 1, np_fftcd_x
               ixyz = mp_fftcd_x(n)
               mz = (ixyz-1)/(lx*ly)+1
               mm = mod(ixyz,(lx*ly))
               if (mm==0) mm=lx*ly
               my = (mm-1)/lx+1
               mx = mod(mm,lx)
               if (mx==0) mx = lx
               jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
               gx = rltv(1,1)*inx(n*2-1)+rltv(1,2)*jnx(n*2-1)+rltv(1,3)*knx(n*2-1)
               gy = rltv(2,1)*inx(n*2-1)+rltv(2,2)*jnx(n*2-1)+rltv(2,3)*knx(n*2-1)
               gz = rltv(3,1)*inx(n*2-1)+rltv(3,2)*jnx(n*2-1)+rltv(3,3)*knx(n*2-1)
               g2 = gx**2 +gy**2 +gz**2
               afft_l(jadd*2-1,1) = -g2 *afft_l(jadd*2-1,1)

               gx = rltv(1,1)*inx(n*2  )+rltv(1,2)*jnx(n*2  )+rltv(1,3)*knx(n*2  )
               gy = rltv(2,1)*inx(n*2  )+rltv(2,2)*jnx(n*2  )+rltv(2,3)*knx(n*2  )
               gz = rltv(3,1)*inx(n*2  )+rltv(3,2)*jnx(n*2  )+rltv(3,3)*knx(n*2  )
               g2 = gx**2 +gy**2 +gz**2
               afft_l(jadd*2  ,1) = -g2 *afft_l(jadd*2,1)
            end do
         end if
#else
         do n = 1, np_fftcd_x*kimg  ! MPI
            gx = rltv(1,1)*inx(n)+rltv(1,2)*jnx(n)+rltv(1,3)*knx(n)
            gy = rltv(2,1)*inx(n)+rltv(2,2)*jnx(n)+rltv(2,3)*knx(n)
            gz = rltv(3,1)*inx(n)+rltv(3,2)*jnx(n)+rltv(3,3)*knx(n)
            g2 = gx**2 +gy**2 +gz**2
            afft_l(n,1) = -g2 *afft_l(n,1)
         enddo
#endif

#ifdef FFT_3D_DIVISION_CD
         call m_FFT_CD_Inverse_3DIV_3D(nfout,afft_l,lsize,1)
#else
         if (sw_fft_xzy > 0) then
            call m_FFT_CD_Inverse_3D(nfout,afft_l,lsize,1)
         else
            call m_FFT_CD_Inverse_XYZ_3D(nfout,afft_l,lsize,1)
         end if
#endif

#ifdef FFT_3D_DIVISION_CD
         kx1p = fftcd_X_x_nel
         kx2p = fftcd_X_y_nel
         kx3p = fftcd_X_z_nel
         lx = fft_box_size_CD_3D(1,0)
         ly = fft_box_size_CD_3D(2,0)
         lz = fft_box_size_CD_3D(3,0)
         do i = 1, np_fftcd_x
            i1 = mp_fftcd_x(i)
            mz = (i1-1)/(lx*ly)+1
            mm = mod(i1,(lx*ly))
            if (mm==0) mm=lx*ly
            my = (mm-1)/lx+1
            mx = mod(mm,lx)
            if (mx==0) mx = lx
            jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2)) &
                 & +kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
            dF_dlaplrho(i,iloop) = afft_l(jadd*2-1,1)
         end do
#else
         do i = 1, nfft_div_size
            dF_dlaplrho(i,iloop) = afft_l(i*2-1,1)
         end do
#endif
      End Do
      deallocate( afft_l )
    end subroutine calc_lapl_dF_dlapl_3D

    subroutine add_vdwdf_nonlocal_energy( version_vdwdf )
      integer, intent(in) :: version_vdwdf
      integer :: i, nfftcd,nfftx,nffty,nfftz,ix,iy,iz,ixyz
      integer :: ix2, iy2, iz2, nlphf,idp, mmp

      real(kind=DP),allocatable,dimension(:) :: chgrhr_l_in
      real(kind=DP),allocatable,dimension(:) :: grad_rho_in
      real(kind=DP),allocatable,dimension(:) :: dfdrho_vdw,dfddrho_vdw

      nfftx = fft_box_size_CD_3D(1,1)
      nffty = fft_box_size_CD_3D(2,1)
      nfftz = fft_box_size_CD_3D(3,1)
      nfftcd = nfftx*nffty*nfftz

      allocate(dfdrho_vdw (nfft_div_size));dfdrho_vdw=0.d0
      allocate(dfddrho_vdw(nfft_div_size));dfddrho_vdw=0.d0
      allocate(chgrhr_l_in(nfft_div_size));chgrhr_l_in=0.0d0
      allocate(grad_rho_in(nfft_div_size));grad_rho_in=0.d0

      if(ispin==1)then
         chgrhr_l_in(:) = chgrhr_l(:,1)
         grad_rho_in(:) = grad_rho(:,1)
      else
         chgrhr_l_in(:) = chgrhr_l(:,1)+chgrhr_l(:,2)
         grad_rho_in(:) = grad_rho(:,1)+grad_rho(:,2)
      endif

      call vdW_scf( nspin, ispin, nfft_div_size, nfftx, nffty, nfftz, chgrhr_l_in, grad_rho_in, &
           &        ecnl, dfdrho_vdw, dfddrho_vdw, version_vdwdf, vflag )

      if ( kimg == 1 ) then
         do i=1,nfft_div_size
            dF_drho(i,1:ispin) = dF_drho(i,1:ispin) + dfdrho_vdw(i)
            dF_dgradrho(i,1:ispin) = dF_dgradrho(i,1:ispin) + dfddrho_vdw(i)
         enddo
      else
         do i=1,nfft_div_size
            dF_drho(i,1:ispin) = dF_drho(i,1:ispin) + dfdrho_vdw(i)
            dF_dgradrho(i,1:ispin) = dF_dgradrho(i,1:ispin) + dfddrho_vdw(i)
         enddo
      endif

      deallocate(dfdrho_vdw); deallocate(dfddrho_vdw)
      deallocate(chgrhr_l_in); deallocate(grad_rho_in)
    end subroutine add_vdwdf_nonlocal_energy

! === Support Hybrid on 3D_Parallel by tkato 2013/02/10 ========================
    subroutine scale_exchange(alpha)
      implicit none
      real(kind=DP), intent(in) :: alpha

      real(kind=DP) :: fexx

      fexx = 1.d0-alpha
      exc = fexx * exc
      dF_drho = fexx * dF_drho
      dF_dgradrho = fexx * dF_dgradrho
    end subroutine scale_exchange
! ==============================================================================

    subroutine sum_s_gga12
      real(kind=DP), allocatable, dimension(:,:)    :: s_gga_mpi        ! MPI d(3,3)
                                                 __TIMER_SUB_START(795)
      call mpi_allreduce(MPI_IN_PLACE,s_gga1,9,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
      call mpi_allreduce(MPI_IN_PLACE,s_gga2,9,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)

!      if(npes>=2 .and. (nrank_ggacmp>1.or.npes_cdfft>1)) allocate(s_gga_mpi(3,3))
!      if(nrank_ggacmp > 1) then
!                              __TIMER_COMM_START_w_BARRIER(mpi_ggacmp_cross_world(myrank_cdfft),1034)
!         call mpi_allreduce(s_gga1,s_gga_mpi,9,mpi_double_precision, mpi_sum &
!              &                                          , mpi_ggacmp_cross_world(myrank_cdfft),ierr)
!                                                 __TIMER_COMM_STOP(1034)
!         s_gga1 = s_gga_mpi
!
!                              __TIMER_COMM_START_w_BARRIER(mpi_ggacmp_cross_world(myrank_cdfft),1035)
!         call mpi_allreduce(s_gga2,s_gga_mpi,9,mpi_double_precision,mpi_sum &
!              &                                          , mpi_ggacmp_cross_world(myrank_cdfft),ierr)
!                                                 __TIMER_COMM_STOP(1035)
!         s_gga2 = s_gga_mpi
!      end if
!
!      if(npes_cdfft >= 2)  then
!         if(myrank_ggacmp < nrank_ggacmp) then
!                              __TIMER_COMM_START_w_BARRIER(mpi_cdfft_world(myrank_ggacmp),1036)
!            call mpi_allreduce(s_gga1,s_gga_mpi,9,mpi_double_precision,mpi_sum &
!                 &                                       , mpi_cdfft_world(myrank_ggacmp),ierr)
!                                                 __TIMER_COMM_STOP(1036)
!            s_gga1 = s_gga_mpi
!                              __TIMER_COMM_START_w_BARRIER(mpi_cdfft_world(myrank_ggacmp),1037)
!            call mpi_allreduce(s_gga2,s_gga_mpi,9,mpi_double_precision, mpi_sum &
!                 &                                       , mpi_cdfft_world(myrank_ggacmp),ierr)
!                                                 __TIMER_COMM_STOP(1037)
!            s_gga2 = s_gga_mpi
!         end if
!      end if
!      if(nrest_cdfft >= 1) then
!                                                 __TIMER_COMM_START_w_BARRIER(MPI_CommGroup,1038)
!         if(mype > npes_cdfft*nrank_ggacmp - 1) then
!            call mpi_recv(s_gga1,9,mpi_double_precision &
!                 & ,mype - npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,istatus,ierr)
!         end if
!         if(mype < nrest_cdfft) then
!            call mpi_send(s_gga1,9,mpi_double_precision &
!                 & ,mype + npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,ierr)
!         end if
!                                                 __TIMER_COMM_STOP(1038)
!
!                                                 __TIMER_COMM_START_w_BARRIER(MPI_CommGroup,1039)
!         if(mype > npes_cdfft*nrank_ggacmp - 1) then
!            call mpi_recv(s_gga2,9,mpi_double_precision &
!                 & ,mype - npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,istatus,ierr)
!         end if
!         if(mype < nrest_cdfft) then
!            call mpi_send(s_gga2,9,mpi_double_precision &
!                 & ,mype + npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,ierr)
!         end if
!                                                 __TIMER_COMM_STOP(1039)
!      end if
!
!!$         call mpi_allreduce(s_gga1,s_gga_mpi,9 &
!!$              &, mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)   ! MPI
!!$         s_gga1 = s_gga_mpi                  ! MPI
!!$
!!$         call mpi_allreduce(s_gga2,s_gga_mpi,9 &
!!$              &, mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)   ! MPI
!!$         s_gga2 = s_gga_mpi                  ! MPI

!      if(npes >= 2 .and. (nrank_ggacmp>1.or.npes_cdfft>1)) deallocate(s_gga_mpi)



                                                 __TIMER_SUB_STOP(795)
    end subroutine sum_s_gga12

#ifndef _XC_SAVE_MEMORY_

    subroutine stress_exchange_part_3D(i1,i2,in,is)
      integer, intent(in) :: i1,i2,in,is
      integer :: m
                                                 __TIMER_SUB_START(797)
                                                 __TIMER_DO_START(1042)
      do m = 1, nfft_div_size       ! MPI
         s_gga1(i1,i2)=s_gga1(i1,i2) &
              &  - univol * rinplw * dF_dgradrho(m,is) &
              &    * cgrad_rho(m,in,is) * dF_drho(m,is) * f2or1(m)
      end do
                                                 __TIMER_DO_STOP(1042)

      ! correlation part -->
                                                 __TIMER_DO_START(1043)
      if(nspin.eq.1) then
         do m = 1, nfft_div_size    ! MPI
            s_gga2(i1,i2)=s_gga2(i1,i2) &
                 &  - univol * rinplw * grad_trho(m) &
                 &    * cgrad_rho(m,in,is) * dF_drho(m,1) * f2or1(m)
         end do
      endif
                                                 __TIMER_DO_STOP(1043)
      ! <--
                                                 __TIMER_SUB_STOP(797)
    end subroutine stress_exchange_part_3D

    subroutine stress_correlation_part_3D(i1,i2,in)
      integer, intent(in) :: i1,i2,in
      integer :: m
                                                 __TIMER_SUB_START(794)
                                                 __TIMER_DO_START(1033)
      do m = 1, nfft_div_size         ! MPI
         s_gga2(i1,i2)=s_gga2(i1,i2) &
              &  - univol * rinplw * grad_trho(m) &
              &    * (cgrad_rho(m,in,1)+cgrad_rho(m,in,2)) * dF_drho(m,1) * f2or1(m)
      end do
                                                 __TIMER_DO_STOP(1033)
                                                 __TIMER_SUB_STOP(794)
    end subroutine stress_correlation_part_3D
#else
    subroutine stress_exchange_part_3D(i1,i2,is)
      integer, intent(in) :: i1,i2,is
      integer :: m

      do m = 1, nfft_div_size       ! MPI
         s_gga1(i1,i2)=s_gga1(i1,i2) &
              &  + univol * rinplw * dF_dgradrho(m,is) &
              &    * cggawk13(m) * dF_drho(m,is) * f2or1(m)
      end do

      if(nspin.eq.1) then
         do m = 1, nfft_div_size    ! MPI
            s_gga2(i1,i2)=s_gga2(i1,i2) &
                 &  + univol * rinplw * grad_trho(m) &
                 &    * cggawk13(m) * dF_drho(m,1) * f2or1(m)
         end do
      endif
    end subroutine stress_exchange_part_3D

    subroutine stress_correlation_part_3D(i1,i2)
      integer, intent(in) :: i1,i2
      integer :: m
                                                 __TIMER_SUB_START(794)
      do m = 1, nfft_div_size         ! MPI
         s_gga2(i1,i2)=s_gga2(i1,i2) &
              &  + univol * rinplw * grad_trho(m) &
              &    * cggawk13(m) * dF_drho(m,1) * f2or1(m)
      end do
                                                 __TIMER_SUB_STOP(794)
    end subroutine stress_correlation_part_3D
#endif

    subroutine rhos_diff(i1,i2)
      integer, intent(in) :: i1,i2
      integer :: is,ig
      integer :: iend  !mpi
                                                 __TIMER_SUB_START(786)
      drhodh = 0.d0
                                                 __TIMER_DO_START(1025)
      if(kimg.eq.1) then
         do is = 1, nspin - af
            iend = iend_kngp
            if( iend_kngp > kg ) iend = kg
            if( ista_kngp <= iend ) then
               do ig = ista_kngp, iend  !for mpi
                  drhodh(ig,1,is) = - chgsoft(ig,1,is) * alinvt(i1,i2)
               enddo
            endif
!!xocl end spread sum(drhodh)
         enddo
      else
         do is = 1, nspin - af
            iend = iend_kngp
            if( iend_kngp > kg ) iend = kg
            if( ista_kngp <= iend ) then
               do ig = ista_kngp, iend  !for mpi
                  drhodh(ig,1,is) = - chgsoft(ig,1,is) * alinvt(i1,i2)
                  drhodh(ig,kimg,is) = - chgsoft(ig,kimg,is) * alinvt(i1,i2)
               enddo
            endif
!!xocl end spread sum(drhodh)
         enddo
      endif
                                                 __TIMER_DO_STOP(1025)
                                                 __TIMER_SUB_STOP(786)
    end subroutine rhos_diff

    subroutine rhopc_diff_3D(i1,i2)
      integer, intent(in) :: i1,i2
      integer :: ia,it,ipc,ig
      real(kind=DP) :: pc
                                                 __TIMER_SUB_START(787)
      do ia=1,natm
      it=ityp(ia)
      if(itpcc(it)==0) cycle
!      call calc_phase2(natm,pos,ia,kgp,ngabc,ista_kngp,iend_kngp,zfcos,zfsin)
      call calc_phase2(natm,pos,ia,kgp,ngabc_kngp_l,ista_kngp,iend_kngp,zfcos,zfsin)
                       ! -(b_Elec.)  -> zfcos, zfsin
      ipc=itpcc(it)
                                                 __TIMER_DO_START(1026)
      if(nspin.eq.1 .or. af.eq.1) then
         if(kimg.eq.1) then
            do ig = ista_kngp, iend_kngp  !for mpi
               pc = (rhpcg_l(ig,ipc)*alinvt(i1,i2) &
                 &   +rhpcg_diff_l(ig,ipc)*grinv(ig)*g(ig,i1) &
                 &   *(g(ig,1)*alinvt(1,i2)+g(ig,2)*alinvt(2,i2) &
                 &    +g(ig,3)*alinvt(3,i2))) &
                 &   * dble(iwei(ia))
               drhodh(ig,1,1) = drhodh(ig,1,1) &
                 &   - zfcos(ig) * pc / nspin
            enddo
         elseif(kimg.eq.2) then
            do ig = ista_kngp, iend_kngp  !for mpi
               pc = (rhpcg_l(ig,ipc)*alinvt(i1,i2) &
                 &   +rhpcg_diff_l(ig,ipc)*grinv(ig)*g(ig,i1) &
                 &   *(g(ig,1)*alinvt(1,i2)+g(ig,2)*alinvt(2,i2) &
                 &    +g(ig,3)*alinvt(3,i2))) &
                 &   * dble(iwei(ia))
               drhodh(ig,1,1) = drhodh(ig,1,1) &
                 &   - zfcos(ig) * pc / nspin
               drhodh(ig,kimg,1) = drhodh(ig,kimg,1) &
                 &   + zfsin(ig) * pc / nspin
            enddo
         endif
      else
         if(kimg.eq.1) then
            do ig = ista_kngp, iend_kngp  !for mpi
               pc = (rhpcg_l(ig,ipc)*alinvt(i1,i2) &
                 &   +rhpcg_diff_l(ig,ipc)*grinv(ig)*g(ig,i1) &
                 &   *(g(ig,1)*alinvt(1,i2)+g(ig,2)*alinvt(2,i2) &
                 &    +g(ig,3)*alinvt(3,i2))) &
                 &   * dble(iwei(ia))
               drhodh(ig,1,1) = drhodh(ig,1,1) &
                 &   - zfcos(ig) * pc / nspin
               drhodh(ig,1,nspin) = drhodh(ig,1,nspin) &
                 &   - zfcos(ig) * pc / nspin
            enddo
         elseif(kimg.eq.2) then
            do ig = ista_kngp, iend_kngp  !for mpi
               pc = (rhpcg_l(ig,ipc)*alinvt(i1,i2) &
                 &   +rhpcg_diff_l(ig,ipc)*grinv(ig)*g(ig,i1) &
                 &   *(g(ig,1)*alinvt(1,i2)+g(ig,2)*alinvt(2,i2) &
                 &    +g(ig,3)*alinvt(3,i2))) &
                 &   * dble(iwei(ia))
               drhodh(ig,1,1) = drhodh(ig,1,1) &
                 &   - zfcos(ig) * pc / nspin
               drhodh(ig,1,nspin) = drhodh(ig,1,nspin) &
                 &   - zfcos(ig) * pc / nspin
               drhodh(ig,kimg,1) = drhodh(ig,kimg,1) &
                 &   + zfsin(ig) * pc / nspin
               drhodh(ig,kimg,nspin) = drhodh(ig,kimg,nspin) &
                 &   + zfsin(ig) * pc / nspin
            enddo
         endif
      endif
                                                 __TIMER_DO_STOP(1026)
      enddo
                                                 __TIMER_SUB_STOP(787)
    end subroutine rhopc_diff_3D

    subroutine rhoh_diff_3D(i1,i2)
      use m_PlaneWaveBasisSet,    only : m_pwBS_sphrp2_3D,m_pwBS_sphrp2_diff_3D

      integer, intent(in) :: i1,i2
      integer    :: n,ilm3,i,it,mdvdb,ia
      integer    :: lmt1,lmt2,il1,il2,tau1,tau2,l3,iiqitg,ii
      real(kind=DP) :: fac,dga
      integer, pointer, dimension(:)  :: il3       ! d(n**2)
      real(kind=DP),allocatable,target,dimension(:) :: ylm_t
      real(kind=DP), pointer, dimension(:,:)   :: ylmd_mpi
      integer :: ibl1,ibl2
#ifndef _m_XC_no_loop_exchange_
      integer :: m, maxm, ip, np, iq
      integer, parameter :: mcritical = 4*2+1
      integer, allocatable, dimension(:) :: nqitg_sp, nqitg_sp0 !d(ntyp)
      integer, allocatable, dimension(:) :: iq2l3 ! d(nqitg)
      integer, allocatable, dimension(:,:) :: nc  ! d(maxm,nqitg)
      integer :: mc ! maxval(nc)
      integer, allocatable, dimension(:,:,:) :: nc2lmt1, nc2lmt2, nc2n ! d(mc,maxm,nqitg)
#endif
      real(kind=DP), allocatable, target, dimension(:,:) :: ylm_ext
      real(kind=DP), allocatable, dimension(:,:,:) :: drhodh_tmp,drhodh_l
      integer :: ierr
                                                 __TIMER_SUB_START(788)
      if(iprixc >= 2) write(6,'(" << rhoh_diff >>")')
      call m_PP_find_maximum_l(n)   ! n-1: maximum l
      n = (n-1) + (n-1) + 1
      allocate(il3(n**2)); call substitute_il3(n**2,il3) ! -(b_Elec..)

      allocate(drhodh_l(ista_kngp:iend_kngp,kimg,nspin)); drhodh_l = drhodh
      drhodh = 0.d0
#ifndef _m_XC_no_loop_exchange_
      allocate(nqitg_sp(ntyp)); allocate(nqitg_sp0(ntyp))
      allocate(iq2l3(nqitg))
      allocate(nc(mcritical,nqitg));nc=0
      call m_PP_set_index_arrays1(nfout,ntyp,nqitg,mcritical,n**2,il3 &
           & ,maxm,mc,nqitg_sp,nqitg_sp0,iq2l3,nc)
      allocate(nc2lmt1(mc,maxm,nqitg))
      allocate(nc2lmt2(mc,maxm,nqitg))
      allocate(nc2n(mc,maxm,nqitg))
      call m_PP_set_index_arrays2(nfout,mc,maxm,nqitg,mcritical,n**2,il3,iq2l3 &
           & ,nc2lmt1,nc2lmt2,nc2n,nc) ! -> nc2lmt1, nc2lmt2, nc2n, nc
#endif
!      allocate(ylm_t(ista_kngp:iend_kngp));ylm_t = 0.d0
      allocate(ylmd(ista_kngp:iend_kngp,3,n**2)); ylmd = 0.d0
      allocate(ylmd_mpi(ista_kngp:iend_kngp,3)); ylmd_mpi = 0.d0
      if(n**2 > nel_Ylm) then
         allocate(ylm_ext(ista_kngp:iend_kngp,nel_Ylm+1:n**2));  ylm_ext = 0.d0
         allocate(ylm_t(ista_kngp:iend_kngp)); ylm_t = 0.d0
         do ilm3 = nel_Ylm+1, n**2
            call m_pwBS_sphrp2_3D(ilm3,rltv,ista_kngp,iend_kngp,ylm_t)  ! (ilm3,rltv,ngabc,gr_l)->(ylm)
            ylm_ext(:,ilm3) = ylm_t(:)
         end do
         deallocate(ylm_t)
      end if
      do ilm3 = 1, n**2
         call m_pwBS_sphrp2_diff_3D(ilm3,rltv,ylmd_mpi)
         do i = ista_kngp, iend_kngp  !for mpi
            ylmd(i,1,ilm3) = ylmd_mpi(i,1)
            ylmd(i,2,ilm3) = ylmd_mpi(i,2)
            ylmd(i,3,ilm3) = ylmd_mpi(i,3)
         enddo
      end do
      deallocate(ylmd_mpi)
      call blksize(n)
      do ibl1 = ista_kngp, iend_kngp, nblk
      ibl2 = ibl1+nblk-1
      if(ibl2 .gt. iend_kngp) ibl2 = iend_kngp
#ifndef _m_XC_no_loop_exchange_
      allocate(zfsin_blk(ibl1:ibl2));zfsin_blk=0.d0
      allocate(zfcos_blk(ibl1:ibl2));zfcos_blk=0.d0
      allocate(mzfsin_blk(ibl1:ibl2));mzfsin_blk=0.d0
!      do ia = 1, natm
      do ia = ista_atm_ke, iend_atm_ke
         it = ityp(ia)
         mdvdb = m_PP_include_vanderbilt_pot(it)
         if(mdvdb == SKIP) cycle
         call calc_phase2_blk(natm,pos,ia,kgp,ngabc_kngp_l,ibl1,ibl2,ista_kngp,iend_kngp,zfcos_blk,zfsin_blk)
         mzfsin_blk = -zfsin_blk
            ! -(b_Elec.)  -> zfcos, zfsin
         do iq = nqitg_sp0(it), nqitg_sp(it)
            l3 = iq2l3(iq)
            do m = 1, 2*l3+1
               ilm3 = l3*l3+m
               if(ilm3 <= nel_Ylm) then
                  ylm => ylm_l(ibl1:ibl2,ilm3)
               else
                  ylm => ylm_ext(ibl1:ibl2,ilm3)
               end if
               do ip = 1, nc(m,iq)
                  lmt1 = nc2lmt1(ip,m,iq)
                  lmt2 = nc2lmt2(ip,m,iq)
                  np = nc2n(ip,m,iq)
                  dga = dl2p(lmt1,lmt2,np,it)
                  fac = 2.d0 ; if(lmt1 == lmt2) fac = 1.d0
                  if(mod(l3,2) == 0) then
                     call even_case_3D(ibl1,ibl2,i1,i2,ilm3,iq,l3,ia,lmt1,lmt2,fac,dga)
                  else
                     call odd_case_3D(ibl1,ibl2,i1,i2,ilm3,iq,l3,ia,lmt1,lmt2,fac,dga)
                  endif
               enddo
            enddo
         enddo
      enddo
      deallocate(mzfsin_blk)
      deallocate(zfsin_blk)
      deallocate(zfcos_blk)
#else
      do ia = 1, natm
         it = ityp(ia)
         mdvdb = m_PP_include_vanderbilt_pot(it)
         if(mdvdb == SKIP) cycle
!         call calc_phase2(natm,pos,ia,kgp,ngabc,ista_kngp,iend_kngp,zfcos,zfsin)
         call calc_phase2(natm,pos,ia,kgp,ngabc_kngp_l,ista_kngp,iend_kngp,zfcos,zfsin)
            ! -(b_Elec.)  -> zfcos, zfsin
         do lmt1 = 1,ilmt(it)
            il1 = ltp(lmt1,it); tau1 = taup(lmt1,it)
            do lmt2 = lmt1, ilmt(it)
              il2 = ltp(lmt2,it); tau2 = taup(lmt2,it)
              fac = 2.d0 ; if(lmt1 == lmt2) fac = 1.d0
              do n = 1, il2p(lmt1,lmt2,it)
                ilm3 = isph(lmt1,lmt2,n,it)
                l3   =  il3(ilm3)
                iiqitg = iqitg(il1,tau1,il2,tau2,l3+1,it)
                if(iiqitg == 0) cycle
                if(ilm3 <= nel_Ylm) then
                   do ii = ista_kngp,iend_kngp
                      ylm(ii-ista_kngp+1) =  ylm_l(ii,ilm3)
                   enddo
                else
                   call m_pwBS_sphrp2_3D(ilm3,rltv,ista_kngp,iend_kngp,ylm_t)
                   do ii = ista_kngp,iend_kngp
                      ylm(ii-ista_kngp+1) =  ylm_t(ii)
                   enddo
                end if
                dga = dl2p(lmt1,lmt2,n,it)
                if(mod(il1+il2,2) == 0) then
                   call even_case_3D(i1,i2,ilm3,iiqitg,l3,ia,lmt1,lmt2,fac,dga)
                else
                   call odd_case_3D(i1,i2,ilm3,iiqitg,l3,ia,lmt1,lmt2,fac,dga)
                endif
              enddo
            enddo
         enddo
      enddo
#endif
      enddo

      allocate(drhodh_tmp(kgp,kimg,nspin));drhodh_tmp=0.d0
      drhodh_tmp(ista_kngp:iend_kngp,:,:) = drhodh(ista_kngp:iend_kngp,:,:)
      call mpi_allreduce( &
      &    mpi_in_place,drhodh_tmp,kgp*kimg*nspin,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
      drhodh(ista_kngp:iend_kngp,:,:) = drhodh_l(ista_kngp:iend_kngp,:,:) &
                                    & + drhodh_tmp(ista_kngp:iend_kngp,:,:)
      deallocate(drhodh_tmp);deallocate(drhodh_l)
#ifndef _m_XC_no_loop_exchange_
      deallocate(nc2n,nc2lmt2,nc2lmt1,nc,iq2l3,nqitg_sp,nqitg_sp0)
#endif
      !deallocate(il3);deallocate(ylm_t);deallocate(ylmd)
      deallocate(il3);deallocate(ylmd)
                                                 __TIMER_SUB_STOP(788)
    end subroutine rhoh_diff_3D

    subroutine blksize(n)
      integer, intent(in) :: n
      integer :: ncache
      ncache = (m_CtrlP_cachesize()*1024)*3/4
      if(ncache == 0) then
         nblk = iend_kngp-ista_kngp+1
      else
         nblk=ncache/(8*(n**2))
      end if
      if (nblk<32) nblk = 1000
      !write(nfout,'(a,i8)') ' !** block size : ',nblk
    end subroutine blksize

    subroutine even_case_3D(ibl1,ibl2,i1,i2,ilm3,iiqitg,l3,ia,lmt1,lmt2,fac,dga)
      integer, intent(in) :: ibl1,ibl2,i1,i2,ilm3,iiqitg,l3,ia,lmt1,lmt2
      real(kind=DP), intent(in) :: fac,dga
      integer :: is
                                                 __TIMER_SUB_START(789)
                                                 __TIMER_DO_START(1027)
      do is=1,nspin,af+1
        flchgq(is) = fac*real(zi**(-l3))*dga &
                   & *hsr(ia,lmt1,lmt2,is)
        flchgqd(i1,i2,is) = fac*real(zi**(-l3))*dga &
                          & *hsrd(ia,lmt1,lmt2,is,i1,i2)
      enddo
                                                 __TIMER_DO_STOP(1027)
      if(kimg == 1) then
         call real_case_3D(ibl1,ibl2,i1,i2,ilm3,iiqitg,ia,zfcos_blk)
      else
         call complex_case_3D(ibl1,ibl2,i1,i2,ilm3,iiqitg,zfcos_blk,mzfsin_blk)
      endif
                                                 __TIMER_SUB_STOP(789)
    end subroutine even_case_3D

    subroutine odd_case_3D(ibl1,ibl2,i1,i2,ilm3,iiqitg,l3,ia,lmt1,lmt2,fac,dga)
      integer, intent(in) :: ibl1,ibl2,i1,i2,ilm3,iiqitg,l3,ia,lmt1,lmt2
      real(kind=DP), intent(in) :: fac,dga
      integer :: is
                                                 __TIMER_SUB_START(790)
                                                 __TIMER_DO_START(1028)
      do is=1,nspin,af+1
        flchgq(is) = fac*aimag(zi**(-l3))*dga &
                   & *hsr(ia,lmt1,lmt2,is)
        flchgqd(i1,i2,is) = fac*aimag(zi**(-l3))*dga &
                          & *hsrd(ia,lmt1,lmt2,is,i1,i2)
      enddo
                                                 __TIMER_DO_STOP(1028)
      if(kimg == 1) then
         call real_case_3D(ibl1,ibl2,i1,i2,ilm3,iiqitg,ia,zfsin_blk)
      else
         call complex_case_3D(ibl1,ibl2,i1,i2,ilm3,iiqitg,zfsin_blk,zfcos_blk)
      endif
                                                 __TIMER_SUB_STOP(790)
    end subroutine odd_case_3D

    subroutine real_case_3D(ibl1,ibl2,i1,i2,ilm3,iiqitg,ia,zf)

      integer,       intent(in) :: ibl1,ibl2,i1,i2,ilm3,iiqitg,ia
      real(kind=DP), intent(in), dimension(ibl1:ibl2) :: zf
      integer :: is, i, iy
      real(kind=DP) :: c1,c2,c3up,c3dn
                                                 __TIMER_SUB_START(791)
                                                 __TIMER_DO_START(1029)
      do is=1,nspin,af+1
        flchgq(is) = flchgq(is) * dble(iwei(ia))
        flchgqd(i1,i2,is) = flchgqd(i1,i2,is) * dble(iwei(ia))
      enddo
                                                 __TIMER_DO_STOP(1029)
                                                 __TIMER_DO_START(1030)
      if(nspin==1 .or. af==1) then
         do i = ibl1, ibl2  !for mpi
            iy = i - ibl1 + 1
            c1=(  g(i,1)*alinvt(1,i2) &
             &  + g(i,2)*alinvt(2,i2) &
             &  + g(i,3)*alinvt(3,i2) )  *g(i,i1) * grinv(i)
            c2=(  ylmd(i,1,ilm3)*alinvt(1,i2) &
             &  + ylmd(i,2,ilm3)*alinvt(2,i2) &
             &  + ylmd(i,3,ilm3)*alinvt(3,i2) )  *g(i,i1)
            c3up =  &
             &  - flchgq(1) &
             &  *(ylm(iy)*(qitg_l(i,iiqitg)*alinvt(i1,i2) &
             &  +qitg_diff_l(i,iiqitg)*c1) &
             &  +c2*qitg_l(i,iiqitg)) &
             &  + flchgqd(i1,i2,1)*ylm(iy)*qitg_l(i,iiqitg)
            drhodh(i,1,1) = drhodh(i,1,1) &
             &  + c3up * zf(i)
         enddo
      else
         do i = ibl1, ibl2  !for mpi
            iy = i - ibl1 + 1
            c1=(  g(i,1)*alinvt(1,i2) &
             &  + g(i,2)*alinvt(2,i2) &
             &  + g(i,3)*alinvt(3,i2) )  *g(i,i1) * grinv(i)
            c2=(  ylmd(i,1,ilm3)*alinvt(1,i2) &
             &  + ylmd(i,2,ilm3)*alinvt(2,i2) &
             &  + ylmd(i,3,ilm3)*alinvt(3,i2) )  *g(i,i1)
            c3up =  &
             &  - flchgq(1) &
             &  *(ylm(iy)*(qitg_l(i,iiqitg)*alinvt(i1,i2) &
             &  +qitg_diff_l(i,iiqitg)*c1) &
             &  +c2*qitg_l(i,iiqitg)) &
             &  + flchgqd(i1,i2,1)*ylm(iy)*qitg_l(i,iiqitg)
            c3dn =  &
             &  - flchgq(nspin) &
             &  *(ylm(iy)*(qitg_l(i,iiqitg)*alinvt(i1,i2) &
             &  +qitg_diff_l(i,iiqitg)*c1) &
             &  +c2*qitg_l(i,iiqitg)) &
             &  + flchgqd(i1,i2,nspin)*ylm(iy)*qitg_l(i,iiqitg)
            drhodh(i,1,1) = drhodh(i,1,1) &
             &  + c3up * zf(i)
            drhodh(i,1,nspin) = drhodh(i,1,nspin) &
             &  + c3dn * zf(i)
         enddo
      endif
                                                 __TIMER_DO_STOP(1030)
                                                 __TIMER_SUB_STOP(791)
    end subroutine real_case_3D

    subroutine complex_case_3D(is,ie,i1,i2,ilm3,iiqitg,zf1,zf2)

      integer,       intent(in) :: is,ie,i1,i2,ilm3,iiqitg
      real(kind=DP), intent(in), dimension(is:ie) &
           &        :: zf1,zf2      ! MPI
      real(kind=DP) :: c1,c2,c3up,c3dn
      integer       :: i, iy

                                                 __TIMER_SUB_START(792)
                                                 __TIMER_DO_START(1031)
      if(nspin==1 .or. af==1) then
         do i = is, ie  !for mpi
            iy = i - is+1
            c1=(  g(i,1)*alinvt(1,i2) &
             &  + g(i,2)*alinvt(2,i2) &
             &  + g(i,3)*alinvt(3,i2) )  *g(i,i1) * grinv(i)
            c2=(  ylmd(i,1,ilm3)*alinvt(1,i2) &
             &  + ylmd(i,2,ilm3)*alinvt(2,i2) &
             &  + ylmd(i,3,ilm3)*alinvt(3,i2) )  *g(i,i1)
            c3up =  &
             &  - flchgq(1) &
             &  *(ylm(iy)*(qitg_l(i,iiqitg)*alinvt(i1,i2) &
             &  +qitg_diff_l(i,iiqitg)*c1) &
             &  +c2*qitg_l(i,iiqitg)) &
             &  + flchgqd(i1,i2,1)*ylm(iy)*qitg_l(i,iiqitg)
            drhodh(i,1,1) = drhodh(i,1,1) &
             &  + c3up * zf1(i)
            drhodh(i,kimg,1) = drhodh(i,kimg,1) &
             &  + c3up * zf2(i)
         enddo
      else
         do i = is, ie !for mpi
            iy = i - is+1
            c1=(  g(i,1)*alinvt(1,i2) &
             &  + g(i,2)*alinvt(2,i2) &
             &  + g(i,3)*alinvt(3,i2) )  *g(i,i1) * grinv(i)
            c2=(  ylmd(i,1,ilm3)*alinvt(1,i2) &
             &  + ylmd(i,2,ilm3)*alinvt(2,i2) &
             &  + ylmd(i,3,ilm3)*alinvt(3,i2) )  *g(i,i1)
            c3up =  &
             &  - flchgq(1) &
             &  *(ylm(iy)*(qitg_l(i,iiqitg)*alinvt(i1,i2) &
             &  +qitg_diff_l(i,iiqitg)*c1) &
             &  +c2*qitg_l(i,iiqitg)) &
             &  + flchgqd(i1,i2,1)*ylm(iy)*qitg_l(i,iiqitg)
            c3dn =  &
             &  - flchgq(nspin) &
             &  *(ylm(iy)*(qitg_l(i,iiqitg)*alinvt(i1,i2) &
             &  +qitg_diff_l(i,iiqitg)*c1) &
             &  +c2*qitg_l(i,iiqitg)) &
             &  + flchgqd(i1,i2,nspin)*ylm(iy)*qitg_l(i,iiqitg)
            drhodh(i,1,1) = drhodh(i,1,1) &
                 &  + c3up * zf1(i)
            drhodh(i,kimg,1) = drhodh(i,kimg,1) &
                 &  + c3up * zf2(i)
            drhodh(i,1,nspin) = drhodh(i,1,nspin) &
                 &  + c3dn * zf1(i)
            drhodh(i,kimg,nspin) = drhodh(i,kimg,nspin) &
                 &  + c3dn * zf2(i)
        enddo
      endif
                                                 __TIMER_DO_STOP(1031)
                                                 __TIMER_SUB_STOP(792)
    end subroutine complex_case_3D

    subroutine map_drhodh1_3D(is)
      integer, intent(in) :: is
      integer :: i,ip,ri
      real(kind=DP), allocatable, dimension(:) :: chgfft_mpi4
                                                 __TIMER_SUB_START(798)

      chgfft_mpi = 0.d0
                                                 __TIMER_DO_START(1044)
      do ri = 1, kimg
!         do i = ista_kngp, iend_kngp  !for mpi
         do i = ista_kngp, min(kgp_reduced, iend_kngp)  !for mpi
            ip = igfp_l(i)*kimg - kimg + ri
            chgfft_mpi(ip)  = chgfft_mpi(ip) + drhodh(i,ri,is)
         end do
      end do
                                                 __TIMER_DO_STOP(1044)

#ifdef CD_FFT_ALL
      if(npes >= 2) then
         allocate(chgfft_mpi4(nfftsize*kimg))
                                                 __TIMER_COMM_START_w_BARRIER(MPI_CommGroup,1045)
         call mpi_allreduce(chgfft_mpi,chgfft_mpi4,nfftsize*kimg,mpi_double_precision, mpi_sum, MPI_CommGroup,ierr)
                                                 __TIMER_COMM_STOP(1045)
         chgfft = chgfft_mpi4
         deallocate(chgfft_mpi4)
#else
      if(nrank_g >= 2) then
         allocate(chgfft_mpi4(nfftsize*kimg))
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,1046)
         call mpi_allreduce(chgfft_mpi,chgfft_mpi4,nfftsize*kimg,mpi_double_precision, mpi_sum, mpi_ke_world,ierr)
                                                 __TIMER_COMM_STOP(1046)
         chgfft = chgfft_mpi4
         deallocate(chgfft_mpi4)
#endif
      else
         chgfft = chgfft_mpi
      end if
                                                 __TIMER_SUB_STOP(798)
    end subroutine map_drhodh1_3D

    subroutine map_drhodh2_3D
      integer :: i,ip,ri
      real(kind=DP), allocatable, dimension(:) :: chgfft_mpi4
                                                 __TIMER_SUB_START(799)
      chgfft_mpi = 0.d0
                                                 __TIMER_DO_START(1047)
      do ri = 1, kimg
!         do i = ista_kngp, iend_kngp  ! for mpi
         do i = ista_kngp, min(kgp_reduced, iend_kngp)  ! for mpi
            ip = igfp_l(i)*kimg - kimg + ri
            chgfft_mpi(ip) = chgfft_mpi(ip) + drhodh(i,ri,1) + drhodh(i,ri,nspin)
         end do
      end do
                                                 __TIMER_DO_STOP(1047)

#ifdef CD_FFT_ALL
      if(npes >= 2) then
         allocate(chgfft_mpi4(nfftsize*kimg))
                                                 __TIMER_COMM_START_w_BARRIER(MPI_CommGroup,1048)
         call mpi_allreduce(chgfft_mpi,chgfft_mpi4,nfftsize*kimg,mpi_double_precision, mpi_sum, MPI_CommGroup,ierr)
                                                 __TIMER_COMM_STOP(1048)
         chgfft = chgfft_mpi4
         deallocate(chgfft_mpi4)
#else
      if(nrank_g >= 2) then
         allocate(chgfft_mpi4(nfftsize*kimg))
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,1049)
         call mpi_allreduce(chgfft_mpi,chgfft_mpi4,nfftsize*kimg,mpi_double_precision, mpi_sum, mpi_ke_world,ierr)
                                                 __TIMER_COMM_STOP(1049)
         chgfft = chgfft_mpi4
         deallocate(chgfft_mpi4)
#endif
      else
         chgfft = chgfft_mpi
      end if
                                                 __TIMER_SUB_STOP(799)
    end subroutine map_drhodh2_3D

    subroutine dgrhodh1_3D(i1,i2,in,is)
      integer, intent(in) :: i1,i2,in,is
      integer       :: n
      real(kind=DP) :: g1,g2

                                                 __TIMER_SUB_START(793)
                                                 __TIMER_DO_START(1032)
#ifdef FFT_3D_DIVISION_CD
      integer :: lx, ly, lz, n, ixyz, mx, my, mz, mm, jadd, kx1p, kx2p, kx3p
      if (kimg == 1) then
         kx1p = fftcd_X_x_nel
         kx2p = fftcd_X_y_nel
         kx3p = fftcd_X_z_nel
         lx = fft_box_size_CD_3D(1,0)
         ly = fft_box_size_CD_3D(2,0)
         lz = fft_box_size_CD_3D(3,0)
         do n = 1, np_fftcd_x
            ixyz = mp_fftcd_x(n)
            mz = (ixyz-1)/(lx*ly)+1
            mm = mod(ixyz,(lx*ly))
            if (mm==0) mm=lx*ly
            my = (mm-1)/lx+1
            mx = mod(mm,lx)
            if (mx==0) mx = lx
            jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
            g1 = rltv(in,1)*inx(n    )+rltv(in,2)*jnx(n    )+rltv(in,3)*knx(n    )
            g2 = rltv(i1,1)*inx(n    )+rltv(i1,2)*jnx(n    )+rltv(i1,3)*knx(n    )
            afft_l(jadd*2-1,1) = g1*chgfft(mp_fftcd_x(n)*2-1) - chden_l(n*2-1,is)*g2*alinvt(in,i2)
            g1 = rltv(in,1)*inx(n  )+rltv(in,2)*jnx(n  )+rltv(in,3)*knx(n  )
            g2 = rltv(i1,1)*inx(n  )+rltv(i1,2)*jnx(n  )+rltv(i1,3)*knx(n  )
            afft_l(jadd*2,1) = g1*chgfft(mp_fftcd_x(n)*2) - chden_l(n*2,is)*g2*alinvt(in,i2)
         end do
      else
         kx1p = fftcd_X_x_nel
         kx2p = fftcd_X_y_nel
         kx3p = fftcd_X_z_nel
         lx = fft_box_size_CD_3D(1,0)
         ly = fft_box_size_CD_3D(2,0)
         lz = fft_box_size_CD_3D(3,0)
         do n = 1, np_fftcd_x
            ixyz = mp_fftcd_x(n)
            mz = (ixyz-1)/(lx*ly)+1
            mm = mod(ixyz,(lx*ly))
            if (mm==0) mm=lx*ly
            my = (mm-1)/lx+1
            mx = mod(mm,lx)
            if (mx==0) mx = lx
            jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
            g1 = rltv(in,1)*inx(n*2-1)+rltv(in,2)*jnx(n*2-1)+rltv(in,3)*knx(n*2-1)
            g2 = rltv(i1,1)*inx(n*2-1)+rltv(i1,2)*jnx(n*2-1)+rltv(i1,3)*knx(n*2-1)
            afft_l(jadd*2-1,1) = g1*chgfft(mp_fftcd_x(n)*2-1) - chden_l(n*2-1,is)*g2*alinvt(in,i2)
            g1 = rltv(in,1)*inx(n*2)+rltv(in,2)*jnx(n*2)+rltv(in,3)*knx(n*2)
            g2 = rltv(i1,1)*inx(n*2)+rltv(i1,2)*jnx(n*2)+rltv(i1,3)*knx(n*2)
            afft_l(jadd*2,1) = g1*chgfft(mp_fftcd_x(n)*2) - chden_l(n*2,is)*g2*alinvt(in,i2)
         end do
      end if
#else
      if (kimg == 1) then
         do n = 1, np_fftcd_x                   ! MPI
            g1 = rltv(in,1)*inx(n)+rltv(in,2)*jnx(n)+rltv(in,3)*knx(n)
            g2 = rltv(i1,1)*inx(n)+rltv(i1,2)*jnx(n)+rltv(i1,3)*knx(n)
            afft_l(n,1) = g1*chgfft(mp_fftcd_x(n)) - chden_l(n,is)*g2*alinvt(in,i2)
         enddo
      else
         do n = 1, np_fftcd_x                     ! MPI
            g1 = rltv(in,1)*inx(n*2-1)+rltv(in,2)*jnx(n*2-1)+rltv(in,3)*knx(n*2-1)
            g2 = rltv(i1,1)*inx(n*2-1)+rltv(i1,2)*jnx(n*2-1)+rltv(i1,3)*knx(n*2-1)
            afft_l(n*2-1,1) = g1*chgfft(mp_fftcd_x(n)*2-1) - chden_l(n*2-1,is)*g2*alinvt(in,i2)
            g1 = rltv(in,1)*inx(n*2)+rltv(in,2)*jnx(n*2)+rltv(in,3)*knx(n*2)
            g2 = rltv(i1,1)*inx(n*2)+rltv(i1,2)*jnx(n*2)+rltv(i1,3)*knx(n*2)
            afft_l(n*2,1) = g1*chgfft(mp_fftcd_x(n)*2) - chden_l(n*2,is)*g2*alinvt(in,i2)
         enddo
      end if
#endif
                                                 __TIMER_DO_STOP(1032)
      call boundary_zero_into_afft_3D(in)  ! -(contained in subr. m_XC_cal_potential)

#ifdef FFT_3D_DIVISION_CD
      call m_FFT_CD_Inverse_3DIV_3D(nfout,afft_l,lsize,1)
#else
      if (sw_fft_xzy > 0) then
         call m_FFT_CD_Inverse_3D(nfout,afft_l,lsize,1)
      else
         call m_FFT_CD_Inverse_XYZ_3D(nfout,afft_l,lsize,1)
      end if
#endif

#ifdef FFT_3D_DIVISION_CD
      kx1p = fftcd_X_x_nel
      kx2p = fftcd_X_y_nel
      kx3p = fftcd_X_z_nel
      lx = fft_box_size_CD_3D(1,0)
      ly = fft_box_size_CD_3D(2,0)
      lz = fft_box_size_CD_3D(3,0)
      do n = 1, np_fftcd_x
         ixyz = mp_fftcd_x(n)
         mz = (ixyz-1)/(lx*ly)+1
         mm = mod(ixyz,(lx*ly))
         if (mm==0) mm=lx*ly
         my = (mm-1)/lx+1
         mx = mod(mm,lx)
         if (mx==0) mx = lx
         jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
         dF_drho(n,is) = afft_l(jadd*2,1)
      end do
#else
      do n = 1, nfft_div_size     ! MPI
         dF_drho(n,is) = afft_l(2*n,1)
      enddo
#endif
                                                 __TIMER_SUB_STOP(793)
    end subroutine dgrhodh1_3D

    subroutine dgrhodh2_3D(i1,i2,in)
      integer, intent(in) :: i1,i2,in
      integer       :: n
      real(kind=DP) :: g1,g2

                                                 __TIMER_SUB_START(796)
                                                 __TIMER_DO_START(1040)
#ifdef FFT_3D_DIVISION_CD
      integer :: lx, ly, lz, n, ixyz, mx, my, mz, mm, jadd, kx1p, kx2p, kx3p
      if (kimg == 1) then
         kx1p = fftcd_X_x_nel
         kx2p = fftcd_X_y_nel
         kx3p = fftcd_X_z_nel
         lx = fft_box_size_CD_3D(1,0)
         ly = fft_box_size_CD_3D(2,0)
         lz = fft_box_size_CD_3D(3,0)
         do n = 1, np_fftcd_x
            ixyz = mp_fftcd_x(n)
            mz = (ixyz-1)/(lx*ly)+1
            mm = mod(ixyz,(lx*ly))
            if (mm==0) mm=lx*ly
            my = (mm-1)/lx+1
            mx = mod(mm,lx)
            if (mx==0) mx = lx
            jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))

            g1 = rltv(in,1)*inx(n    )+rltv(in,2)*jnx(n    )+rltv(in,3)*knx(n    )
            g2 = rltv(i1,1)*inx(n    )+rltv(i1,2)*jnx(n    )+rltv(i1,3)*knx(n    )
            afft_l(jadd*2-1,1) = g1*chgfft(mp_fftcd_x(n)*2-1) &
                 &   - (chden_l(n*2-1,1)+chden_l(n*2-1,nspin))*g2*alinvt(in,i2)
            g1 = rltv(in,1)*inx(n    )+rltv(in,2)*jnx(n    )+rltv(in,3)*knx(n    )
            g2 = rltv(i1,1)*inx(n    )+rltv(i1,2)*jnx(n    )+rltv(i1,3)*knx(n    )
            afft_l(jadd*2  ,1) = g1*chgfft(mp_fftcd_x(n)*2  ) &
                 &   - (chden_l(n*2  ,1)+chden_l(n*2  ,nspin))*g2*alinvt(in,i2)

         end do
      else
         kx1p = fftcd_X_x_nel
         kx2p = fftcd_X_y_nel
         kx3p = fftcd_X_z_nel
         lx = fft_box_size_CD_3D(1,0)
         ly = fft_box_size_CD_3D(2,0)
         lz = fft_box_size_CD_3D(3,0)
         do n = 1, np_fftcd_x
            ixyz = mp_fftcd_x(n)
            mz = (ixyz-1)/(lx*ly)+1
            mm = mod(ixyz,(lx*ly))
            if (mm==0) mm=lx*ly
            my = (mm-1)/lx+1
            mx = mod(mm,lx)
            if (mx==0) mx = lx
            jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))

            g1 = rltv(in,1)*inx(n*2-1)+rltv(in,2)*jnx(n*2-1)+rltv(in,3)*knx(n*2-1)
            g2 = rltv(i1,1)*inx(n*2-1)+rltv(i1,2)*jnx(n*2-1)+rltv(i1,3)*knx(n*2-
1)
            afft_l(jadd*2-1,1) = g1*chgfft(mp_fftcd_x(n)*2-1) &
                 &   - (chden_l(n*2-1,1)+chden_l(n*2-1,nspin))*g2*alinvt(in,i2)
            g1 = rltv(in,1)*inx(n*2  )+rltv(in,2)*jnx(n*2  )+rltv(in,3)*knx(n*2  )
            g2 = rltv(i1,1)*inx(n*2  )+rltv(i1,2)*jnx(n*2  )+rltv(i1,3)*knx(n*2
 )
            afft_l(jadd*2  ,1) = g1*chgfft(mp_fftcd_x(n)*2  ) &
                 &   - (chden_l(n*2  ,1)+chden_l(n*2  ,nspin))*g2*alinvt(in,i2)

         end do
      end if
#else
      if (kimg == 1) then
         do n = 1, np_fftcd_x
            g1 = rltv(in,1)*inx(n)+rltv(in,2)*jnx(n)+rltv(in,3)*knx(n)
            g2 = rltv(i1,1)*inx(n)+rltv(i1,2)*jnx(n)+rltv(i1,3)*knx(n)
            afft_l(n,1) = g1*chgfft(mp_fftcd_x(n)) &
                 &   - (chden_l(n,1)+chden_l(n,nspin))*g2*alinvt(in,i2)
         enddo
      else
         do n = 1, np_fftcd_x
            g1 = rltv(in,1)*inx(n*2-1)+rltv(in,2)*jnx(n*2-1)+rltv(in,3)*knx(n*2-1)
            g2 = rltv(i1,1)*inx(n*2-1)+rltv(i1,2)*jnx(n*2-1)+rltv(i1,3)*knx(n*2-1)
            afft_l(n*2-1,1) = g1*chgfft(mp_fftcd_x(n)*2-1) &
                 &   - (chden_l(n*2-1,1)+chden_l(n*2-1,nspin))*g2*alinvt(in,i2)
            g1 = rltv(in,1)*inx(n*2  )+rltv(in,2)*jnx(n*2  )+rltv(in,3)*knx(n*2  )
            g2 = rltv(i1,1)*inx(n*2  )+rltv(i1,2)*jnx(n*2  )+rltv(i1,3)*knx(n*2  )
            afft_l(n*2  ,1) = g1*chgfft(mp_fftcd_x(n)*2  ) &
                 &   - (chden_l(n*2  ,1)+chden_l(n*2  ,nspin))*g2*alinvt(in,i2)
         enddo
      end if
#endif
                                                 __TIMER_DO_STOP(1040)
      call boundary_zero_into_afft_3D(in)   ! -(contained in subr. m_XC_cal_potential)

! === DEBUG by tkato 2013/09/26 ================================================
!     call m_FFT_CD_Inverse_3D(nfout,afft_l,lsize,1)
#ifdef FFT_3D_DIVISION_CD
      call m_FFT_CD_Inverse_3DIV_3D(nfout,afft_l,lsize,1)
#else
      if (sw_fft_xzy > 0) then
         call m_FFT_CD_Inverse_3D(nfout,afft_l,lsize,1)
      else
         call m_FFT_CD_Inverse_XYZ_3D(nfout,afft_l,lsize,1)
      end if
#endif
! ==============================================================================

                                                 __TIMER_DO_START(1041)
#ifdef FFT_3D_DIVISION_CD
         kx1p = fftcd_X_x_nel
         kx2p = fftcd_X_y_nel
         kx3p = fftcd_X_z_nel
         lx = fft_box_size_CD_3D(1,0)
         ly = fft_box_size_CD_3D(2,0)
         lz = fft_box_size_CD_3D(3,0)
         do n = 1, np_fftcd_x
            ixyz = mp_fftcd_x(n)
            mz = (ixyz-1)/(lx*ly)+1
            mm = mod(ixyz,(lx*ly))
            if (mm==0) mm=lx*ly
            my = (mm-1)/lx+1
            mx = mod(mm,lx)
            if (mx==0) mx = lx
            jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
            dF_drho(n,1) = afft_l(jadd*2,1)
         end do
#else
      do n = 1, nfft_div_size     ! MPI
         dF_drho(n,1) = afft_l(2*n,1)
      enddo
#endif
                                                 __TIMER_DO_STOP(1041)
                                                 __TIMER_SUB_STOP(796)
    end subroutine dgrhodh2_3D

#ifdef _XC_SAVE_MEMORY_
    subroutine gradrho1(in,is)
      integer, intent(in) :: in,is
      integer :: m
      call g_xyz_chden_l(in,is)          ! G_xyz * rho(G) --> afft
      call m_FFT_CD_inverse_c(nfout,afft)  ! (-i)*d(rho(r))/d(x|y|z)
      do m = ista_fftph, iend_fftph      ! MPI
         cggawk13(m) = afft(2*m)
      enddo
    end subroutine gradrho1

    subroutine gradrho2(in)
      integer, intent(in) ::in
      integer :: m
      call g_xyz_total_chden_l(in)       ! G_xyz*(rho(G)up+rho(G)down) -> afft
      call m_FFT_CD_inverse_c(nfout,afft)  !(-i)*d(rho_total(r))/d(x|y|z)
      do m = ista_fftph, iend_fftph      ! MPI
         cggawk13(m) = afft(2*m)
      enddo
    end subroutine gradrho2
#endif

    subroutine allocate_for_stress_3D

      integer :: i
      real(kind=DP) :: ga,gb,gc
      allocate(alinvt(3,3)); alinvt = rltv / PAI2
      allocate(drhodh(ista_kngp:iend_kngp,kimg,nspin)); drhodh = 0.d0
      allocate(zfcos(ista_kngp:iend_kngp)); zfcos = 0.d0
      allocate(zfsin(ista_kngp:iend_kngp)); zfsin = 0.d0
      allocate(flchgq(nspin)); flchgq = 0.d0
      allocate(flchgqd(3,3,nspin)); flchgqd = 0.d0
!      allocate(chgfft(ista_fftp:iend_fftp)); chgfft = 0.d0  ! MPI
      allocate(chgfft(nfftsize*kimg)); chgfft = 0.d0  ! MPI
      allocate(grinv(ista_kngp:iend_kngp)); grinv = 0.d0
      allocate(g(ista_kngp:iend_kngp,3)); g = 0.d0

      do i = ista_kngp, iend_kngp  !for mpi
        if(gr_l(i).lt.1.0d-20) then
          grinv(i)=0.d0
        else
          grinv(i)=1.d0/gr_l(i)
        endif
      enddo

      do i = ista_kngp, iend_kngp  !for mpi
        ga = ngabc_kngp_l(i,1)
        gb = ngabc_kngp_l(i,2)
        gc = ngabc_kngp_l(i,3)
        g(i,1) = rltv(1,1)*ga+rltv(1,2)*gb+rltv(1,3)*gc
        g(i,2) = rltv(2,1)*ga+rltv(2,2)*gb+rltv(2,3)*gc
        g(i,3) = rltv(3,1)*ga+rltv(3,2)*gb+rltv(3,3)*gc
      enddo
      allocate(chgfft_mpi(nfftsize*kimg))             ! MPI

!      allocate(ylm(1:(iend_kngp-ista_kngp+1)));ylm=0.d0

    end subroutine allocate_for_stress_3D

    subroutine deallocate_for_stress
      deallocate(alinvt)
      deallocate(drhodh)
      deallocate(zfcos)
      deallocate(zfsin)
      deallocate(flchgq)
      deallocate(flchgqd)
      deallocate(chgfft)
      deallocate(grinv)
      deallocate(g)
      deallocate(chgfft_mpi)
!      deallocate(ylm)
    end subroutine deallocate_for_stress

    subroutine abs_grad_rho_up_down_total_3D()
      real(kind=DP),allocatable,dimension(:)    :: grad_rho_c2  ! MPI d(1:nfft_div_size)

      integer  :: is, in, i
      real(kind=DP) :: x,y,z
#ifdef XC_PACK_FFT
      real(kind=DP) ,allocatable, dimension(:,:) :: bfft_l
#ifdef FFT_3D_DIVISION_CD
      allocate(bfft_l(lsize*2,3))
#else
      allocate(bfft_l(lsize*kimg,3))
#endif
#endif

                                                 __TIMER_SUB_START(759)
      do is = 1, ispin
         grad_rho(:,is) = 0.d0
#ifdef XC_PACK_FFT
         do in = 1, 3
            call g_xyz_chden_l_3D(in,is)        ! G_xyz * rho_{up|down}(G)-->afft
#ifdef FFT_3D_DIVISION_CD
            integer :: lx, ly, lz, ixyz, mx, my, mz, mm, jadd, kx1p, kx2p, kx3

            kx1p = fftcd_X_x_nel
            kx2p = fftcd_X_y_nel
            kx3p = fftcd_X_z_nel
            lx = fft_box_size_CD_3D(1,0)
            ly = fft_box_size_CD_3D(2,0)
            lz = fft_box_size_CD_3D(3,0)
            do i = 1, np_fftcd_x
               ixyz = mp_fftcd_x(i)
               mz = (ixyz-1)/(lx*ly)+1
               mm = mod(ixyz,(lx*ly))
               if (mm==0) mm=lx*ly
               my = (mm-1)/lx+1
               mx = mod(mm,lx)
               if (mx==0) mx = lx
               jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
               bfft_l(jadd*2-1,in) = afft_l(jadd*2-1,1)
               bfft_l(jadd*2  ,in) = afft_l(jadd*2  ,1)
            end do
#else
            do i = 1, np_fftcd_x*kimg
               bfft_l(i,in) = afft_l(i,1)
            end do
#endif
         end do
#ifdef FFT_3D_DIVISION_CD
         call m_FFT_CD_Inverse_3DIV_3D(nfout,bfft_l,lsize,3)
#else
         if (sw_fft_xzy > 0) then
            call m_FFT_CD_Inverse_3D(nfout,bfft_l,lsize,3)
         else
            call m_FFT_CD_Inverse_XYZ_3D(nfout,bfft_l,lsize,3)
         end if
#endif
!                                                 __TIMER_DO_START(892)
                                                 __TIMER_SUB_START(762)
                                                 __TIMER_SUB_START(763)
                                                 __TIMER_DO_START(896)
                                                 __TIMER_DO_START(897)
         do in = 1, 3
#ifdef FFT_3D_DIVISION_CD
            kx1p = fftcd_X_x_nel
            kx2p = fftcd_X_y_nel
            kx3p = fftcd_X_z_nel
            lx = fft_box_size_CD_3D(1,0)
            ly = fft_box_size_CD_3D(2,0)
            lz = fft_box_size_CD_3D(3,0)
            do i = 1, np_fftcd_x
               ixyz = mp_fftcd_x(i)
               mz = (ixyz-1)/(lx*ly)+1
               mm = mod(ixyz,(lx*ly))
               if (mm==0) mm=lx*ly
               my = (mm-1)/lx+1
               mx = mod(mm,lx)
               if (mx==0) mx = lx
               jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
               cgrad_rho(i,in,is) = -bfft_l(jadd*2,in)
               grad_rho(i,is) = grad_rho(i,is) + bfft_l(jadd*2,in)*bfft_l(jadd*2,in)
            end do
#else
!OCL PARALLEL_STRONG
            do i = 1, nfft_div_size
               cgrad_rho(i,in,is) = -bfft_l(i*2,in)
               grad_rho(i,is) = grad_rho(i,is) + bfft_l(i*2,in)*bfft_l(i*2,in)
#endif
            end do
         end do
                                                 __TIMER_DO_STOP(896)
                                                 __TIMER_DO_STOP(897)
                                                 __TIMER_SUB_STOP(762)
                                                 __TIMER_SUB_STOP(763)
!                                                 __TIMER_DO_STOP(892)
#else
         do in = 1, 3
            call g_xyz_chden_l_3D(in,is)        ! G_xyz * rho_{up|down}(G)-->afft
#ifdef FFT_3D_DIVISION_CD
            call m_FFT_CD_Inverse_3DIV_3D(nfout,afft_l,lsize,1)
#else
            if (sw_fft_xzy > 0) then
               call m_FFT_CD_Inverse_3D(nfout,afft_l,lsize,1)
            else
               call m_FFT_CD_Inverse_XYZ_3D(nfout,afft_l,lsize,1)
            end if
#endif
                                                                ! (-i)*d(rho_{up|down}(r))/d(x|y|z)
            call cp_afft_to_cgrad_rho_3D(is,in) ! -> cgrad_rho(i,in,is) = -afft(i*2)
            call add_sq_afft_to_grad_rho_3D(is) ! grad_rho <--  + afft**2
            if(iprixc >= 2) then
               write(nfout,'(" !XC after add_sq_afft_to_grad_rho: is, in = ",2i8)') is, in
               write(nfout,'(" !XC -- grad_rho -- <<abs_grad_rho_up_down_total>>")')
               write(nfout,'(" !XC ",6f12.6)') (grad_rho(i,1),i=1,+17)
               write(nfout,'(" !XC ",6f12.6)') (grad_rho(i,1),i=nfft_div_size-17,nfft_div_size)
               write(nfout,'(" !XC -- cgrad_rho -- <<abs_grad_rho_up_down_total>>")')
               write(nfout,'(" !XC ",6f12.6)') (cgrad_rho(i,in,1),i=1,1+17)
               write(nfout,'(" !XC ",6f12.6)') (cgrad_rho(i,in,1),i=nfft_div_size-17,nfft_div_size)
            end if

         end do
#endif
         if(iprixc >= 2) write(nfout,'(" !XC after add_sq_afft_to_grad_rho")')
         grad_rho(:,is) = dsqrt(grad_rho(:,is))
      enddo
#ifdef XC_PACK_FFT
      deallocate(bfft_l)
#endif
      if(ispin == 2) then
                                                 __TIMER_DO_START(893)
         do in = 1, nfft_div_size
            x = cgrad_rho(in,1,1) + cgrad_rho(in,1,2)
            y = cgrad_rho(in,2,1) + cgrad_rho(in,2,2)
            z = cgrad_rho(in,3,1) + cgrad_rho(in,3,2)
            grad_trho(in) = dsqrt(x*x+y*y+z*z)
         end do
                                                 __TIMER_DO_STOP(893)
      else
         grad_trho(:) = grad_rho(:,1)
      end if
                                                 __TIMER_SUB_STOP(759)
    end subroutine abs_grad_rho_up_down_total_3D

    subroutine grad2_rho_up_down_3D()
      integer  :: is, i
      real(kind=DP) :: x,y,z

      do is = 1, ispin
         grad2_rho(:,is) = 0.d0
         call g2_xyz_chden_l_3D(is)        ! G2 * rho_{up|down}(G)-->afft
#ifdef FFT_3D_DIVISION_CD
         call m_FFT_CD_Inverse_3DIV_3D(nfout,afft_l,lsize,1)
#else
         if (sw_fft_xzy > 0) then
            call m_FFT_CD_Inverse_3D(nfout,afft_l,lsize,1)
         else
            call m_FFT_CD_Inverse_XYZ_3D(nfout,afft_l,lsize,1)
         end if
#endif
         call cp_afft_to_grad2_rho_3D(is)
      enddo

    end subroutine grad2_rho_up_down_3D

    subroutine dFxc_over_ddgradrho_3D()
      real(kind=DP),allocatable,dimension(:)    :: grad_rho_c2  ! MPI
      real(kind=DP),allocatable,dimension(:)    :: grad_rho_c4  ! MPI
      integer  :: is, in
      real(kind=DP),pointer,dimension(:)  :: x
#ifdef XC_PACK_FFT
      real(kind=DP) ,allocatable, dimension(:,:) :: bfft_l
      integer  :: m
      real(kind=DP) :: gxyz
#endif
                                                 __TIMER_SUB_START(772)

#ifdef XC_PACK_FFT
#ifdef FFT_3D_DIVISION_CD
      allocate(bfft_l(lsize*2,3))
#else
      allocate(bfft_l(lsize*kimg,3))
#endif
      bfft_l = 0.d0

      do is = 1, ispin
         grad_rho(:,is) = 0.d0
         do in = 1, 3
                                                 __TIMER_SUB_START(773)
                                                 __TIMER_DO_START(1007)
            if(ispin == 2) then
#ifdef FFT_3D_DIVISION_CD
               integer :: lx, ly, lz, ixyz, mx, my, mz, mm, jadd, kx1p, kx2p, kx3p
               kx1p = fftcd_X_x_nel
               kx2p = fftcd_X_y_nel
               kx3p = fftcd_X_z_nel
               lx = fft_box_size_CD_3D(1,0)
               ly = fft_box_size_CD_3D(2,0)
               lz = fft_box_size_CD_3D(3,0)
               do m = 1, np_fftcd_x
                  ixyz = mp_fftcd_x(m)
                  mz = (ixyz-1)/(lx*ly)+1
                  mm = mod(ixyz,(lx*ly))
                  if (mm==0) mm=lx*ly
                  my = (mm-1)/lx+1
                  mx = mod(mm,lx)
                  if (mx==0) mx = lx
                  jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
                  bfft_l(2*jadd,in) = rinplw*(dF_dgradrho(m,is)*cgrad_rho(m,in,is) &
                 &            +grad_trho(m)*(cgrad_rho(m,in,1)+cgrad_rho(m,in,2)) )
               end do
#else
               do m = 1, nfft_div_size
                  bfft_l(2*m,in) = rinplw*(dF_dgradrho(m,is)*cgrad_rho(m,in,is) &
                 &            +grad_trho(m)*(cgrad_rho(m,in,1)+cgrad_rho(m,in,2)) )
               end do
#endif
            else if(ispin == 1) then
#ifdef FFT_3D_DIVISION_CD
               kx1p = fftcd_X_x_nel
               kx2p = fftcd_X_y_nel
               kx3p = fftcd_X_z_nel
               lx = fft_box_size_CD_3D(1,0)
               ly = fft_box_size_CD_3D(2,0)
               lz = fft_box_size_CD_3D(3,0)
               do m = 1, np_fftcd_x
                  ixyz = mp_fftcd_x(m)
                  mz = (ixyz-1)/(lx*ly)+1
                  mm = mod(ixyz,(lx*ly))
                  if (mm==0) mm=lx*ly
                  my = (mm-1)/lx+1
                  mx = mod(mm,lx)
                  if (mx==0) mx = lx
                  jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
                  bfft_l(2*jadd,in) = rinplw*(dF_dgradrho(m,is)+grad_trho(m))*cgrad_rho(m,in,is)
               end do
#else
               do m = 1, nfft_div_size
                  bfft_l(2*m,in) = rinplw*(dF_dgradrho(m,is)+grad_trho(m))*cgrad_rho(m,in,is)
               end do
#endif
            end if
         end do
                                                 __TIMER_DO_STOP(1007)
                                                 __TIMER_SUB_STOP(773)
#ifdef FFT_3D_DIVISION_CD
         call m_FFT_CD_Direct_3DIV_3D(nfout,bfft_l,lsize,3)
#else
         if (sw_fft_xzy > 0) then
            call m_FFT_CD_Direct_3D(nfout,bfft_l,lsize,3)
         else
            call m_FFT_CD_Direct_XYZ_3D(nfout,bfft_l,lsize,3)
         end if
#endif
                                                 __TIMER_SUB_START(774)
                                                 __TIMER_DO_START(1008)
         do in = 1, 3
#ifdef FFT_3D_DIVISION_CD
            if (kimg == 1) then
               kx1p = fftcd_X_x_nel
               kx2p = fftcd_X_y_nel
               kx3p = fftcd_X_z_nel
               lx = fft_box_size_CD_3D(1,0)
               ly = fft_box_size_CD_3D(2,0)
               lz = fft_box_size_CD_3D(3,0)
               do m = 1, np_fftcd_x
                  ixyz = mp_fftcd_x(m)
                  mz = (ixyz-1)/(lx*ly)+1
                  mm = mod(ixyz,(lx*ly))
                  if (mm==0) mm=lx*ly
                  my = (mm-1)/lx+1
                  mx = mod(mm,lx)
                  if (mx==0) mx = lx
                  jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
                  gxyz = rltv(in,1)*inx(m    )+rltv(in,2)*jnx(m    )+rltv(in,3)*knx(m    )
                  afft_l(jadd*2-1,1) = gxyz*bfft_l(m*2-1,in)
                  gxyz = rltv(in,1)*inx(m    )+rltv(in,2)*jnx(m    )+rltv(in,3)*knx(m    )
                  afft_l(jadd*2  ,1) = gxyz*bfft_l(m*2  ,in)
               end do
            else
               kx1p = fftcd_X_x_nel
               kx2p = fftcd_X_y_nel
               kx3p = fftcd_X_z_nel
               lx = fft_box_size_CD_3D(1,0)
               ly = fft_box_size_CD_3D(2,0)
               lz = fft_box_size_CD_3D(3,0)
               do m = 1, np_fftcd_x
                  ixyz = mp_fftcd_x(m)
                  mz = (ixyz-1)/(lx*ly)+1
                  mm = mod(ixyz,(lx*ly))
                  if (mm==0) mm=lx*ly
                  my = (mm-1)/lx+1
                  mx = mod(mm,lx)
                  if (mx==0) mx = lx
                  jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
                  gxyz = rltv(in,1)*inx(m*2-1)+rltv(in,2)*jnx(m*2-1)+rltv(in,3)*knx(m*2-1)
                  afft_l(jadd*2-1,1) = gxyz*bfft_l(m*2-1,in)
                  gxyz = rltv(in,1)*inx(m*2  )+rltv(in,2)*jnx(m*2  )+rltv(in,3)*knx(m*2  )
                  afft_l(jadd*2  ,1) = gxyz*bfft_l(m*2  ,in)
               end do
            end if
#else
            do m = 1, np_fftcd_x*kimg   ! MPI
               gxyz = rltv(in,1)*inx(m)+rltv(in,2)*jnx(m)+rltv(in,3)*knx(m)
               afft_l(m,1) = gxyz*bfft_l(m,in)
            enddo
#endif
            call boundary_zero_into_afft_3D(in)  ! -(contained in subr. m_XC_cal_potential)
#ifdef FFT_3D_DIVISION_CD
               kx1p = fftcd_X_x_nel
               kx2p = fftcd_X_y_nel
               kx3p = fftcd_X_z_nel
               lx = fft_box_size_CD_3D(1,0)
               ly = fft_box_size_CD_3D(2,0)
               lz = fft_box_size_CD_3D(3,0)
               do m = 1, np_fftcd_x
                  ixyz = mp_fftcd_x(m)
                  mz = (ixyz-1)/(lx*ly)+1
                  mm = mod(ixyz,(lx*ly))
                  if (mm==0) mm=lx*ly
                  my = (mm-1)/lx+1
                  mx = mod(mm,lx)
                  if (mx==0) mx = lx
                  jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
                  bfft_l(jadd*2-1,in) = afft_l(jadd*2-1,1)
                  bfft_l(jadd*2  ,in) = afft_l(jadd*2  ,1)
               end do
#else
            do m = 1, np_fftcd_x*kimg   ! MPI
               bfft_l(m,in) = afft_l(m,1)
            end do
#endif
         end do
                                                 __TIMER_DO_STOP(1008)
                                                 __TIMER_SUB_STOP(774)
#ifdef FFT_3D_DIVISION_CD
         call m_FFT_CD_Inverse_3DIV_3D(nfout,bfft_l,lsize,3)
#else
         if (sw_fft_xzy > 0) then
            call m_FFT_CD_Inverse_3D(nfout,bfft_l,lsize,3)
         else
            call m_FFT_CD_Inverse_XYZ_3D(nfout,bfft_l,lsize,3)
         end if
#endif
                                                 __TIMER_SUB_START(775)
                                                 __TIMER_DO_START(1009)
         do in = 1, 3
#ifdef FFT_3D_DIVISION_CD
               kx1p = fftcd_X_x_nel
               kx2p = fftcd_X_y_nel
               kx3p = fftcd_X_z_nel
               lx = fft_box_size_CD_3D(1,0)
               ly = fft_box_size_CD_3D(2,0)
               lz = fft_box_size_CD_3D(3,0)
               do m = 1, np_fftcd_x
                  ixyz = mp_fftcd_x(m)
                  mz = (ixyz-1)/(lx*ly)+1
                  mm = mod(ixyz,(lx*ly))
                  if (mm==0) mm=lx*ly
                  my = (mm-1)/lx+1
                  mx = mod(mm,lx)
                  if (mx==0) mx = lx
                  jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
                  grad_rho(m,is) = grad_rho(m,is) - bfft_l(2*jadd-1,in)
               end do
#else
            do m = 1, nfft_div_size      ! MPI
               grad_rho(m,is) = grad_rho(m,is) - bfft_l(2*m-1,in)
            end do
#endif
            if(ipri >= 2) then
               write(6,'(" !XC -- afft -- ,is = ",i8," <<add_negative_afft_to_grad_rho>>")') is
               write(6,'(" !XC ",6f12.6)') (bfft_l(m*2-1,in),m=1,18)
               write(6,'(" !XC -- grad_rho -- ,is = ",i8," <<add_negative_afft_to_grad_rho>>")') is
               write(6,'(" !XC ",6f12.6)') (grad_rho(m,is),m=1,18)
            end if
         end do
                                                 __TIMER_DO_STOP(1009)
                                                 __TIMER_SUB_STOP(775)
      end do
      deallocate(bfft_l)

#else
!else #ifdef XC_PACK_FFT

      allocate(grad_rho_c2(1:nfft_div_size),grad_rho_c4(1:nfft_div_size))   ! MPI

      do is = 1, ispin
         grad_rho(:,is) = 0.d0
         do in = 1, 3
            call dFxc_dgradrho_dot_gradrho2_3D(rinplw,in,ispin,is)
             ! grad_{in}(rho_{is})*dFx/d|grad(rho_{is})| + grad_{in}(rho)*dFc/d|grad(rho)|  -> afft
#ifdef FFT_3D_DIVISION_CD
            call m_FFT_CD_Direct_3DIV_3D(nfout,afft_l,lsize,1)
#else
            if (sw_fft_xzy > 0) then
               call m_FFT_CD_Direct_3D(nfout,afft_l,lsize,1)
            else
               call m_FFT_CD_Direct_XYZ_3D(nfout,afft_l,lsize,1)
            end if
#endif
            call G_xyz_afft_3D(in)                                       ! G_{in}q_{in}(G)  -> afft
#ifdef FFT_3D_DIVISION_CD
            call m_FFT_CD_Inverse_3DIV_3D(nfout,afft_l,lsize,1)
#else
            if (sw_fft_xzy > 0) then
               call m_FFT_CD_Inverse_3D(nfout,afft_l,lsize,1)
            else
               call m_FFT_CD_Inverse_XYZ_3D(nfout,afft_l,lsize,1)
            end if
#endif
            call add_negative_afft_to_grad_rho_3D(is)                            ! afft -> grad_rho
         end do
!        if(nrank_g > 1) then
!           grad_rho_c2(:) = grad_rho(:,is)
!           call mpi_allreduce(grad_rho_c2,grad_rho_c4,nfft_div_size &
!                & ,mpi_double_precision, mpi_sum, mpi_ke_world,ierr)
!           grad_rho(:,is) = grad_rho_c4(:)
!        end if
      end do
      deallocate(grad_rho_c4,grad_rho_c2)    ! MPI

#endif
!endif #ifdef XC_PACK_FFT

                                                 __TIMER_SUB_STOP(772)
    end subroutine dFxc_over_ddgradrho_3D

    subroutine finally_gga_xc_pot(is)
      integer, intent(in) :: is
      integer       :: m

      do m = ista_fftph, iend_fftph   ! MPI
         chgrhr_l(m,is) = dF_drho(m,is)+grad_rho(m,is)
      end do
      if ( use_metagga ) then
         if ( xctype /= "tb09" ) then
            do m = ista_fftph, iend_fftph   ! MPI
               chgrhr_l(m,is) = chgrhr_l(m,is) +dF_dlaplrho(m,is)
            end do
         endif
      end if

      if(iprixc >= 2) then
         write(nfout,'(" !XC -- chgrhr_l -- <<finally_gga_xc_pot>>")')
         write(nfout,'(" !XC ",6f12.6)') (chgrhr_l(i,1),i=ista_fftph,ista_fftph+17)
         write(nfout,'(" !XC -- grad_rho -- <<finally_gga_xc_pot>>")')
         write(nfout,'(" !XC ",6f12.6)') (grad_rho(i,1),i=ista_fftph,ista_fftph+17)
         write(nfout,'(" !XC -- dF_drho -- <<finally_gga_xc_pot>>")')
         write(nfout,'(" !XC ",6f12.6)') (dF_drho(i,1),i=ista_fftph,ista_fftph+17)
      end if
    end subroutine finally_gga_xc_pot

    subroutine finally_gga_xc_pot_3D(is)
      integer, intent(in) :: is
      integer       :: m
                                                 __TIMER_SUB_START(776)
                                                 __TIMER_DO_START(1010)
      do m = 1, nfft_div_size   ! MPI
         chgrhr_l(m,is) = dF_drho(m,is)+grad_rho(m,is)
      end do
      if ( use_metagga ) then
         if ( xctype /= "tb09" ) then
            do m = 1, nfft_div_size   ! MPI
               chgrhr_l(m,is) = chgrhr_l(m,is) +dF_dlaplrho(m,is)
            end do
         endif
      end if
                                                 __TIMER_DO_STOP(1010)
      if(iprixc >= 2) then
         write(nfout,'(" !XC -- chgrhr_l -- <<finally_gga_xc_pot>>")')
! === DEBUG by tkato 2013/08/28 ================================================
!        write(nfout,'(" !XC ",6f12.6)') (chgrhr_l(i,1),i=ista_fftph,ista_fftph+17)
         write(nfout,'(" !XC ",6f12.6)') (chgrhr_l(i,1),i=1,18)
! ==============================================================================
         write(nfout,'(" !XC -- grad_rho -- <<finally_gga_xc_pot>>")')
! === DEBUG by tkato 2013/08/28 ================================================
!        write(nfout,'(" !XC ",6f12.6)') (grad_rho(i,1),i=ista_fftph,ista_fftph+17)
         write(nfout,'(" !XC ",6f12.6)') (grad_rho(i,1),i=1,18)
! ==============================================================================
         write(nfout,'(" !XC -- dF_drho -- <<finally_gga_xc_pot>>")')
! === DEBUG by tkato 2013/08/28 ================================================
!        write(nfout,'(" !XC ",6f12.6)') (dF_drho(i,1),i=ista_fftph,ista_fftph+17)
         write(nfout,'(" !XC ",6f12.6)') (dF_drho(i,1),i=1,18)
! ==============================================================================
      end if
                                                 __TIMER_SUB_STOP(776)
    end subroutine finally_gga_xc_pot_3D

    subroutine add_negative_afft_to_grad_rho_3D(is)
      integer, intent(in) :: is
      integer             :: m
                                                 __TIMER_SUB_START(775)
                                                 __TIMER_DO_START(1009)
#ifdef FFT_3D_DIVISION_CD
          integer :: lx, ly, lz, ixyz, mx, my, mz, mm, jadd, kx1p, kx2p, kx3p
          kx1p = fftcd_X_x_nel
          kx2p = fftcd_X_y_nel
          kx3p = fftcd_X_z_nel
          lx = fft_box_size_CD_3D(1,0)
          ly = fft_box_size_CD_3D(2,0)
          lz = fft_box_size_CD_3D(3,0)
          do m = 1, np_fftcd_x
             ixyz = mp_fftcd_x(m)
             mz = (ixyz-1)/(lx*ly)+1
             mm = mod(ixyz,(lx*ly))
             if (mm==0) mm=lx*ly
             my = (mm-1)/lx+1
             mx = mod(mm,lx)
             if (mx==0) mx = lx
             jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
             grad_rho(m,is) = grad_rho(m,is) - afft_l(2*jadd-1,1)
          end do
#else
      do m = 1, nfft_div_size      ! MPI
         grad_rho(m,is) = grad_rho(m,is) - afft_l(2*m-1,1)
      end do
#endif
                                                 __TIMER_DO_STOP(1009)
      if(iprixc >= 2) then
         write(6,'(" !XC -- afft -- ,is = ",i8," <<add_negative_afft_to_grad_rho>>")') is
         write(6,'(" !XC ",6f12.6)') (afft_l(m*2-1,1),m=1,18)

         write(6,'(" !XC -- grad_rho -- ,is = ",i8," <<add_negative_afft_to_grad_rho>>")') is
         write(6,'(" !XC ",6f12.6)') (grad_rho(m,is),m=1,18)
      end if
                                                 __TIMER_SUB_STOP(775)
    end subroutine add_negative_afft_to_grad_rho_3D

    subroutine G_xyz_afft_3D(in)
      integer, intent(in) :: in
      integer       :: n, i,j,k,ip
      real(kind=DP) :: gxyz
                                                 __TIMER_SUB_START(774)
                                                 __TIMER_DO_START(1008)
#ifdef FFT_3D_DIVISION_CD
      integer :: lx, ly, lz, ixyz, mx, my, mz, mm, jadd, kx1p, kx2p, kx3p
      if (kimg == 1) then
         kx1p = fftcd_X_x_nel
         kx2p = fftcd_X_y_nel
         kx3p = fftcd_X_z_nel
         lx = fft_box_size_CD_3D(1,0)
         ly = fft_box_size_CD_3D(2,0)
         lz = fft_box_size_CD_3D(3,0)
         do n = 1, np_fftcd_x
            ixyz = mp_fftcd_x(n)
            mz = (ixyz-1)/(lx*ly)+1
            mm = mod(ixyz,(lx*ly))
            if (mm==0) mm=lx*ly
            my = (mm-1)/lx+1
            mx = mod(mm,lx)
            if (mx==0) mx = lx
            jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
            gxyz = rltv(in,1)*inx(n    )+rltv(in,2)*jnx(n    )+rltv(in,3)*knx(n    )
            afft_l(jadd*2-1,1) = gxyz*afft_l(jadd*2-1,1)
            gxyz = rltv(in,1)*inx(n    )+rltv(in,2)*jnx(n    )+rltv(in,3)*knx(n    )
            afft_l(jadd*2  ,1) = gxyz*afft_l(jadd*2  ,1)
         end do
      else
         kx1p = fftcd_X_x_nel
         kx2p = fftcd_X_y_nel
         kx3p = fftcd_X_z_nel
         lx = fft_box_size_CD_3D(1,0)
         ly = fft_box_size_CD_3D(2,0)
         lz = fft_box_size_CD_3D(3,0)
         do n = 1, np_fftcd_x
            ixyz = mp_fftcd_x(n)
            mz = (ixyz-1)/(lx*ly)+1
            mm = mod(ixyz,(lx*ly))
            if (mm==0) mm=lx*ly
            my = (mm-1)/lx+1
            mx = mod(mm,lx)
            if (mx==0) mx = lx
            jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
            gxyz = rltv(in,1)*inx(n*2-1)+rltv(in,2)*jnx(n*2-1)+rltv(in,3)*knx(n*2-1)
            afft_l(jadd*2-1,1) = gxyz*afft_l(jadd*2-1,1)
            gxyz = rltv(in,1)*inx(n*2  )+rltv(in,2)*jnx(n*2  )+rltv(in,3)*knx(n*2  )
            afft_l(jadd*2  ,1) = gxyz*afft_l(jadd*2  ,1)
         end do
      end if
#else
      do n = 1, np_fftcd_x*kimg   ! MPI
         gxyz = rltv(in,1)*inx(n)+rltv(in,2)*jnx(n)+rltv(in,3)*knx(n)
         afft_l(n,1) = gxyz*afft_l(n,1)
      enddo
#endif
                                                 __TIMER_DO_STOP(1008)
      call boundary_zero_into_afft_3D(in)  ! -(contained in subr. m_XC_cal_potential)
                                                 __TIMER_SUB_STOP(774)
    end subroutine G_xyz_afft_3D

    subroutine dFxc_dgradrho_dot_gradrho2_3D(rinplw,in,ispin,is)
      real(kind=DP), intent(in) :: rinplw
      integer, intent(in)       :: ispin,in,is
      integer                   :: m
                                                 __TIMER_SUB_START(773)
      afft_l(:,1) = 0.d0
                                                 __TIMER_DO_START(1007)
      if(ispin == 2) then
#ifdef FFT_3D_DIVISION_CD
            integer :: lx, ly, lz, ixyz, mx, my, mz, mm, jadd, kx1p, kx2p, kx3p
            kx1p = fftcd_X_x_nel
            kx2p = fftcd_X_y_nel
            kx3p = fftcd_X_z_nel
            lx = fft_box_size_CD_3D(1,0)
            ly = fft_box_size_CD_3D(2,0)
            lz = fft_box_size_CD_3D(3,0)
            do m = 1, np_fftcd_x
               ixyz = mp_fftcd_x(m)
               mz = (ixyz-1)/(lx*ly)+1
               mm = mod(ixyz,(lx*ly))
               if (mm==0) mm=lx*ly
               my = (mm-1)/lx+1
               mx = mod(mm,lx)
               if (mx==0) mx = lx
               jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
               afft_l(2*jadd,1) = rinplw*(dF_dgradrho(m,is)*cgrad_rho(m,in,is) &
                    &            +grad_trho(m)*(cgrad_rho(m,in,1)+cgrad_rho(m,in,2)) )
            end do
#else
         do m = 1, nfft_div_size
            afft_l(2*m,1) = rinplw*(dF_dgradrho(m,is)*cgrad_rho(m,in,is) &
                 &            +grad_trho(m)*(cgrad_rho(m,in,1)+cgrad_rho(m,in,2)) )
         end do
#endif
      else if(ispin == 1) then
#ifdef FFT_3D_DIVISION_CD
            kx1p = fftcd_X_x_nel
            kx2p = fftcd_X_y_nel
            kx3p = fftcd_X_z_nel
            lx = fft_box_size_CD_3D(1,0)
            ly = fft_box_size_CD_3D(2,0)
            lz = fft_box_size_CD_3D(3,0)
            do m = 1, np_fftcd_x
               ixyz = mp_fftcd_x(m)
               mz = (ixyz-1)/(lx*ly)+1
               mm = mod(ixyz,(lx*ly))
               if (mm==0) mm=lx*ly
               my = (mm-1)/lx+1
               mx = mod(mm,lx)
               if (mx==0) mx = lx
               jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
               afft_l(2*jadd,1) = rinplw*(dF_dgradrho(m,is)+grad_trho(m))*cgrad_rho(m,in,is)
            end do
#else
         do m = 1, nfft_div_size
            afft_l(2*m,1) = rinplw*(dF_dgradrho(m,is)+grad_trho(m))*cgrad_rho(m,in,is)
         end do
#endif
      end if
                                                 __TIMER_DO_STOP(1007)
                                                 __TIMER_SUB_STOP(773)
    end subroutine dFxc_dgradrho_dot_gradrho2_3D

    subroutine g_xyz_total_chden_l(in)
      integer, intent(in) ::in
      integer       :: n
      real(kind=DP) :: gxyz

      do n = ista_fftp, iend_fftp   ! MPI
         gxyz = rltv(in,1)*inx(n)+rltv(in,2)*jnx(n)+rltv(in,3)*knx(n)
         afft(n) = gxyz*(chden_l(n, 1)+chden_l(n,ispin))
      enddo
      call boundary_zero_into_afft(in)  ! -(contained in subr. m_XC_cal_potential)

    end subroutine g_xyz_total_chden_l

    subroutine add_sq_afft_to_grad_trho
      integer   :: i
      do i = ista_fftph, iend_fftph     ! MPI
         grad_trho(i) = grad_trho(i) + afft(i*2)**2
      end do
    end subroutine add_sq_afft_to_grad_trho

    subroutine add_sq_afft_to_grad_rho(is)
      integer, intent(in) :: is
      integer :: i
      do i = ista_fftph, iend_fftph     ! MPI
         grad_rho(i,is) = grad_rho(i,is) + afft(i*2)**2
      end do
    end subroutine add_sq_afft_to_grad_rho

    subroutine add_sq_afft_to_grad_rho_3D(is)
      integer, intent(in) :: is
      integer :: i
                                                 __TIMER_SUB_START(763)
                                                 __TIMER_DO_START(897)
#ifdef FFT_3D_DIVISION_CD
      integer :: lx, ly, lz, ixyz, mx, my, mz, mm, jadd, kx1p, kx2p, kx3p
         kx1p = fftcd_X_x_nel
         kx2p = fftcd_X_y_nel
         kx3p = fftcd_X_z_nel
         lx = fft_box_size_CD_3D(1,0)
         ly = fft_box_size_CD_3D(2,0)
         lz = fft_box_size_CD_3D(3,0)
         do i = 1, np_fftcd_x
            ixyz = mp_fftcd_x(i)
            mz = (ixyz-1)/(lx*ly)+1
            mm = mod(ixyz,(lx*ly))
            if (mm==0) mm=lx*ly
            my = (mm-1)/lx+1
            mx = mod(mm,lx)
            if (mx==0) mx = lx
            jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
            grad_rho(i,is) = grad_rho(i,is) + afft_l(jadd*2,1)**2
         end do
#else
      do i = 1, nfft_div_size     ! MPI
         grad_rho(i,is) = grad_rho(i,is) + afft_l(i*2,1)**2
      end do
#endif
                                                 __TIMER_DO_STOP(897)
                                                 __TIMER_SUB_STOP(763)
    end subroutine add_sq_afft_to_grad_rho_3D

    subroutine g_xyz_chden_l(in,is)
      integer, intent(in) :: in,is
      integer       :: n
      real(kind=DP) :: gxyz

!!$      deallocate(afft)
!!$      allocate(afft(ista_fftp:iend_fftp))

      afft = 0.d0
      do n = ista_fftp, iend_fftp  ! MPI
         gxyz = rltv(in,1)*inx(n)+rltv(in,2)*jnx(n)+rltv(in,3)*knx(n)
         afft(n) = gxyz*chden_l(n,is)
      enddo
      call boundary_zero_into_afft(in)  ! -(contained in subr. m_XC_cal_potential)
    end subroutine g_xyz_chden_l

    subroutine g_xyz_chden_l_3D(in,is)
      integer, intent(in) :: in,is
      integer       :: n
      real(kind=DP) :: gxyz
                                                 __TIMER_SUB_START(760)
      afft_l = 0.d0
                                                 __TIMER_DO_START(894)
#ifdef FFT_3D_DIVISION_CD
      integer :: lx, ly, lz, n, ixyz, mx, my, mz, mm, jadd, kx1p, kx2p, kx3p
      if (kimg == 1) then
         kx1p = fftcd_X_x_nel
         kx2p = fftcd_X_y_nel
         kx3p = fftcd_X_z_nel
         lx = fft_box_size_CD_3D(1,0)
         ly = fft_box_size_CD_3D(2,0)
         lz = fft_box_size_CD_3D(3,0)
         do n = 1, np_fftcd_x
            ixyz = mp_fftcd_x(n)
            mz = (ixyz-1)/(lx*ly)+1
            mm = mod(ixyz,(lx*ly))
            if (mm==0) mm=lx*ly
            my = (mm-1)/lx+1
            mx = mod(mm,lx)
            if (mx==0) mx = lx
            jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
            gxyz = rltv(in,1)*inx(n    )+rltv(in,2)*jnx(n    )+rltv(in,3)*knx(n    )
            afft_l(jadd*2-1,1) = gxyz*chden_l(n*2-1,is)
            gxyz = rltv(in,1)*inx(n    )+rltv(in,2)*jnx(n    )+rltv(in,3)*knx(n    )
            afft_l(jadd*2  ,1) = gxyz*chden_l(n*2  ,is)
         end do
      else
         kx1p = fftcd_X_x_nel
         kx2p = fftcd_X_y_nel
         kx3p = fftcd_X_z_nel
         lx = fft_box_size_CD_3D(1,0)
         ly = fft_box_size_CD_3D(2,0)
         lz = fft_box_size_CD_3D(3,0)
         do n = 1, np_fftcd_x
            ixyz = mp_fftcd_x(n)
            mz = (ixyz-1)/(lx*ly)+1
            mm = mod(ixyz,(lx*ly))
            if (mm==0) mm=lx*ly
            my = (mm-1)/lx+1
            mx = mod(mm,lx)
            if (mx==0) mx = lx
            jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
            gxyz = rltv(in,1)*inx(n*2-1)+rltv(in,2)*jnx(n*2-1)+rltv(in,3)*knx(n*2-1)
            afft_l(jadd*2-1,1) = gxyz*chden_l(n*2-1,is)
            gxyz = rltv(in,1)*inx(n*2  )+rltv(in,2)*jnx(n*2  )+rltv(in,3)*knx(n*2  )
            afft_l(jadd*2  ,1) = gxyz*chden_l(n*2  ,is)
         end do
      end if
#else
      do n = 1, np_fftcd_x*kimg  ! MPI
         gxyz = rltv(in,1)*inx(n)+rltv(in,2)*jnx(n)+rltv(in,3)*knx(n)
         afft_l(n,1) = gxyz*chden_l(n,is)
      enddo
#endif
                                                 __TIMER_DO_STOP(894)
      call boundary_zero_into_afft_3D(in)  ! -(contained in subr. m_XC_cal_potential)
                                                 __TIMER_SUB_STOP(760)
    end subroutine g_xyz_chden_l_3D

    subroutine g2_xyz_chden_l_3D(is)
      integer, intent(in) :: is
      integer       :: n
      real(kind=DP) :: gx, gy, gz, g2

      afft_l = 0.d0

#ifdef FFT_3D_DIVISION_CD
      integer :: lx, ly, lz, n, ixyz, mx, my, mz, mm, jadd, kx1p, kx2p, kx3p
      if (kimg == 1) then
         kx1p = fftcd_X_x_nel
         kx2p = fftcd_X_y_nel
         kx3p = fftcd_X_z_nel
         lx = fft_box_size_CD_3D(1,0)
         ly = fft_box_size_CD_3D(2,0)
         lz = fft_box_size_CD_3D(3,0)
         do n = 1, np_fftcd_x
            ixyz = mp_fftcd_x(n)
            mz = (ixyz-1)/(lx*ly)+1
            mm = mod(ixyz,(lx*ly))
            if (mm==0) mm=lx*ly
            my = (mm-1)/lx+1
            mx = mod(mm,lx)
            if (mx==0) mx = lx
            jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
            gx = rltv(1,1)*inx(n    )+rltv(1,2)*jnx(n    )+rltv(1,3)*knx(n    )
            gy = rltv(2,1)*inx(n    )+rltv(2,2)*jnx(n    )+rltv(2,3)*knx(n    )
            gz = rltv(3,1)*inx(n    )+rltv(3,2)*jnx(n    )+rltv(3,3)*knx(n    )
            g2 = gx**2 +gy**2 +gz**2
            afft_l(jadd*2-1,1) = g2 *chden_l(n*2-1,is)

            gx = rltv(1,1)*inx(n    )+rltv(1,2)*jnx(n    )+rltv(1,3)*knx(n    )
            gy = rltv(2,1)*inx(n    )+rltv(2,2)*jnx(n    )+rltv(2,3)*knx(n    )
            gz = rltv(3,1)*inx(n    )+rltv(3,2)*jnx(n    )+rltv(3,3)*knx(n    )
            g2 = gx**2 +gy**2 +gz**2

            afft_l(jadd*2  ,1) = g2 *chden_l(n*2  ,is)
         end do
      else
         kx1p = fftcd_X_x_nel
         kx2p = fftcd_X_y_nel
         kx3p = fftcd_X_z_nel
         lx = fft_box_size_CD_3D(1,0)
         ly = fft_box_size_CD_3D(2,0)
         lz = fft_box_size_CD_3D(3,0)
         do n = 1, np_fftcd_x
            ixyz = mp_fftcd_x(n)
            mz = (ixyz-1)/(lx*ly)+1
            mm = mod(ixyz,(lx*ly))
            if (mm==0) mm=lx*ly
            my = (mm-1)/lx+1
            mx = mod(mm,lx)
            if (mx==0) mx = lx
            jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
            gx = rltv(1,1)*inx(n*2-1)+rltv(1,2)*jnx(n*2-1)+rltv(1,3)*knx(n*2-1)
            gy = rltv(2,1)*inx(n*2-1)+rltv(2,2)*jnx(n*2-1)+rltv(2,3)*knx(n*2-1)
            gz = rltv(3,1)*inx(n*2-1)+rltv(3,2)*jnx(n*2-1)+rltv(3,3)*knx(n*2-1)
            g2 = gx**2 +gy**2 +gz**2
            afft_l(jadd*2-1,1) = g2 *chden_l(n*2-1,is)

            gx = rltv(1,1)*inx(n*2  )+rltv(1,2)*jnx(n*2  )+rltv(1,3)*knx(n*2  )
            gy = rltv(2,1)*inx(n*2  )+rltv(2,2)*jnx(n*2  )+rltv(2,3)*knx(n*2  )
            gz = rltv(3,1)*inx(n*2  )+rltv(3,2)*jnx(n*2  )+rltv(3,3)*knx(n*2  )
            afft_l(jadd*2  ,1) = g**2 *chden_l(n*2  ,is)
         end do
      end if
#else
      do n = 1, np_fftcd_x*kimg  ! MPI
         gx = rltv(1,1)*inx(n)+rltv(1,2)*jnx(n)+rltv(1,3)*knx(n)
         gy = rltv(2,1)*inx(n)+rltv(2,2)*jnx(n)+rltv(2,3)*knx(n)
         gz = rltv(3,1)*inx(n)+rltv(3,2)*jnx(n)+rltv(3,3)*knx(n)
         g2 = gx**2 +gy**2 +gz**2
         afft_l(n,1) = g2 *chden_l(n,is)
      enddo
#endif
    end subroutine g2_xyz_chden_l_3D

    subroutine check_lmn_even_3D
      integer             :: nlmn, i
                                                 __TIMER_SUB_START(756)
                                                 __TIMER_DO_START(889)
      do i = 1, 3
         nlmn = fft_box_size_CD_3D(i,1)/2
         if(2*nlmn == fft_box_size_CD_3D(i,1)) then
            lmn_even(i) = .true.
         else
            lmn_even(i) = .false.
         end if
      end do
                                                 __TIMER_DO_STOP(889)
                                                 __TIMER_SUB_STOP(756)
    end subroutine check_lmn_even_3D

    subroutine boundary_zero_into_afft(in)
      integer, intent(in) :: in
      integer             :: i,j,k,nn,n,idp,nlp,nmp,nnp,nd2p, j0

      nlp = fft_box_size_CD(1,1)
      nmp = fft_box_size_CD(2,1)
      nnp = fft_box_size_CD(3,1)
#ifdef _MPIFFT_
      idp = fft_box_size_CD_c(1,0)
      nd2p = fft_box_size_CD_c(2,0)
#else
      idp = fft_box_size_CD(1,0)
      nd2p = fft_box_size_CD(2,0)
#endif

      if(kimg == 1) then
         if(lmn_even(in)) then
            if( in == 1) then
#ifdef _MPIFFT_
               do j = 1, ly_d*lz
!!$               do j0 = 1, ly_d
!!$                  do k = 1, nnp
!!$                     nn = nlp/2 + 1 + idp*(j0-1) + idp*ly_d*(k-1) + idp*ly_d*lz*myrank_cdfft
                  nn = nlp/2 + 1 + idp*(j-1) + idp*ly_d*lz*myrank_cdfft
                  afft(nn) = 0.d0
!!$                  end do
               end do
#else
               do j = 1, nmp      ! y
                  do k = 1, nnp   ! z
!!$                  nn = nlp/2 + 1 + idp*(j-1) + idp*nmp*(k-1)
                     nn = nlp/2 + 1 + idp*(j-1) + idp*nd2p*(k-1)
                     if(nn >= ista_fftp .and. nn <= iend_fftp) afft(nn) = 0.d0
                  end do
               end do
#endif
            else if( in == 2) then
#ifdef _MPIFFT_
               j0 = nmp/2 + 1 - ny_d*myrank_cdfft
               if(j0 >= 1 .and. j0 <= ny_d) then
                  do i = 1, idp
                     do k = 1, nnp
                        nn = i + idp*(j0-1) + idp*ly_d*(k-1) + idp*ly_d*lz*myrank_cdfft
                        afft(nn) = 0.d0
                     end do
                  end do
               end if
#else
               do i = 1, idp      ! x
                  do k = 1, nnp   ! z
                     nn = i + idp*(nmp/2) + idp*nd2p*(k-1)
                     if(nn >= ista_fftp .and. nn <= iend_fftp) afft(nn) = 0.d0
                  end do
               end do
#endif
            else if(in == 3) then
#ifdef _MPIFFT_
               do i = 1, idp
                  do j = 1+ny_d*myrank_cdfft, min(nmp,ny_d*(myrank_cdfft+1))
                     nn = i + idp*(j-1-ny_d*myrank_cdfft) + idp*ly_d*(nnp/2) &
                          &                          + idp*ly_d*lz*myrank_cdfft
                     afft(nn) = 0.d0
                  end do
               end do
#else
               do i = 1, idp      ! x
                  do j = 1, nmp   ! y
                     nn = i + idp*(j-1) + idp*nd2p*(nnp/2)
                     if(nn >= ista_fftp .and. nn <= iend_fftp) afft(nn) = 0.d0
                  end do
               end do
#endif
            end if
         end if
      else if(kimg == 2) then
         if(lmn_even(in)) then
            if( in == 1) then
#ifdef _MPIFFT_
               do j = 1, ly_d*lz
                  nn = nlp/2 + 1 + idp*(j-1) + idp*ly_d*lz*myrank_cdfft
                  n = nn*2 - 1
                  afft(n) = 0.d0
                  afft(n+1) = 0.d0
               end do
#else
               do j = 1, nmp
                  do k = 1, nnp
!!$                  nn = nlp/2 + 1 + idp*(j-1) + idp*nmp*(k-1)
                     nn = nlp/2 + 1 + idp*(j-1) + idp*nd2p*(k-1)
!!$                     nn = nlp/2 + 1 + idp*(j-1) -idp*ly_d*myrank_cdfft &
!!$                          & +idp*ly_d*(k-1+myrank_cdfft)
                     n  = nn*2 - 1
                     if(n >= ista_fftp .and. n <= iend_fftp) afft(n) = 0.d0
                     n  = nn*2
                     if(n >= ista_fftp .and. n <= iend_fftp) afft(n) = 0.d0
                  end do
               end do
#endif
            else if( in == 2) then
#ifdef _MPIFFT_
               j0 = nmp/2 + 1 - ny_d*myrank_cdfft
               if(j0 >= 1 .and. j0 <= ny_d) then
                  do i = 1, idp
                     do k = 1, nnp
                        nn = i + idp*(j0-1) + idp*ly_d*(k-1) + idp*ly_d*lz*myrank_cdfft
                        n  = nn*2 - 1
                        afft(n) = 0.d0; afft(n+1) = 0.d0
                     end do
                  end do
               end if
#else
               do i = 1, idp
                  do k = 1, nnp
!!$                     nn = i + idp*(nmp/2) + idp*nmp*(k-1)
                     nn = i + idp*(nmp/2) + idp*nd2p*(k-1)
                     n  = nn*2 - 1
                     if(n >= ista_fftp .and. n <= iend_fftp) afft(n) = 0.d0
                     n  = nn*2
                     if(n >= ista_fftp .and. n <= iend_fftp) afft(n) = 0.d0
                  end do
               end do
#endif
            else if(in == 3) then
#ifdef _MPIFFT_
               do i = 1, idp
                  do j = 1+ny_d*myrank_cdfft, min(nmp,ny_d*(myrank_cdfft+1))
                     nn = i + idp*(j-1-ny_d*myrank_cdfft) + idp*ly_d*(nnp/2) &
                          &                          + idp*ly_d*lz*myrank_cdfft
                     n  = nn*2 - 1
                     afft(n) = 0.d0
                     afft(n+1) = 0.d0
                  end do
               end do
#else
               do i = 1, idp
                  do j = 1, nmp
!!$                  nn = i + idp*(j-1) + idp*nmp*(nnp/2)
                     nn = i + idp*(j-1) + idp*nd2p*(nnp/2)
                     n  = nn*2 - 1
                     if(n >= ista_fftp .and. n <= iend_fftp) afft(n) = 0.d0
                     n  = nn*2
                     if(n >= ista_fftp .and. n <= iend_fftp) afft(n) = 0.d0
                  end do
               end do
#endif
            end if
         end if
      end if
    end subroutine boundary_zero_into_afft

    subroutine boundary_zero_into_afft_3D(in)
      integer, intent(in) :: in
      integer             :: lx,ly,lz,nlp,nmp,nnp,mx,my,mz
      integer             :: i,j,k,nn,n,mm,i1,j0
! #ifdef FFT_3D_DIVISION_CD
!       integer :: ixyz, mx, my, mz, mm, jadd, kx1p, kx2p, kx3p
! #endif

                                                 __TIMER_SUB_START(761)
      lx = fft_box_size_CD_3D(1,0)
      ly = fft_box_size_CD_3D(2,0)
      lz = fft_box_size_CD_3D(3,0)

      nlp = fft_box_size_CD_3D(1,1)
      nmp = fft_box_size_CD_3D(2,1)
      nnp = fft_box_size_CD_3D(3,1)
#ifdef FFT_3D_DIVISION_CD
      kx1p = fftcd_X_x_nel
      kx2p = fftcd_X_y_nel
      kx3p = fftcd_X_z_nel
#endif
                                                 __TIMER_DO_START(895)
      if(kimg == 1) then
         if(lmn_even(in)) then
            if( in == 1) then
               do n = 1, np_fftcd_x
                  i1 = mp_fftcd_x(n)
                  mz = (i1-1)/(lx*ly)+1
                  mm = mod(i1,(lx*ly))
                  if (mm==0) mm=lx*ly
                  my = (mm-1)/lx+1
                  mx = mod(mm,lx)
                  if (mx==0) mx = lx
                  if ((nlp/2+1) == mx) then
#ifdef FFT_3D_DIVISION_CD
                     jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
                     afft_l(jadd*2-1,1) = 0.0d0
#else
                     afft_l(n,1) = 0.0d0
#endif
                  end if
               end do
            else if( in == 2) then
               do n = 1, np_fftcd_x
                  i1 = mp_fftcd_x(n)
                  mz = (i1-1)/(lx*ly)+1
                  mm = mod(i1,(lx*ly))
                  if (mm==0) mm=lx*ly
                  my = (mm-1)/lx+1
                  if ((nmp/2+1) == my) then
#ifdef FFT_3D_DIVISION_CD
                     mx = mod(mm,lx)
                     if (mx==0) mx = lx
                     jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
                     afft_l(jadd*2-1,1) = 0.0d0
#else
                     afft_l(n,1) = 0.0d0
#endif
                  end if
               end do
            else if(in == 3) then
               do n = 1, np_fftcd_x
                  i1 = mp_fftcd_x(n)
                  mz = (i1-1)/(lx*ly)+1
                  if ((nnp/2+1) == mz) then
#ifdef FFT_3D_DIVISION_CD
                     mm = mod(i1,(lx*ly))
                     if (mm==0) mm=lx*ly
                     my = (mm-1)/lx+1
                     mx = mod(mm,lx)
                     if (mx==0) mx = lx
                     jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
                     afft_l(jadd*2-1,1) = 0.0d0
#else
                     afft_l(n,1) = 0.0d0
#endif
                  end if
               end do
            end if
         end if
      else if(kimg == 2) then
         if(lmn_even(in)) then
            if( in == 1) then
               do n = 1, np_fftcd_x
                  i1 = mp_fftcd_x(n)
                  mz = (i1-1)/(lx*ly)+1
                  mm = mod(i1,(lx*ly))
                  if (mm==0) mm=lx*ly
                  my = (mm-1)/lx+1
                  mx = mod(mm,lx)
                  if (mx==0) mx = lx
                  if ((nlp/2+1) == mx) then
#ifdef FFT_3D_DIVISION_CD
                     jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
                     afft_l(jadd*2-1,1) = 0.0d0
                     afft_l(jadd*2  ,1) = 0.0d0
#else
                     afft_l(n*2-1,1) = 0.0d0
                     afft_l(n*2  ,1) = 0.0d0
#endif
                  end if
               end do
            else if( in == 2) then
               do n = 1, np_fftcd_x
                  i1 = mp_fftcd_x(n)
                  mz = (i1-1)/(lx*ly)+1
                  mm = mod(i1,(lx*ly))
                  if (mm==0) mm=lx*ly
                  my = (mm-1)/lx+1
                  if ((nmp/2+1) == my) then
#ifdef FFT_3D_DIVISION_CD
                     mx = mod(mm,lx)
                     if (mx==0) mx = lx
                     jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
                     afft_l(jadd*2-1,1) = 0.0d0
                     afft_l(jadd*2  ,1) = 0.0d0
#else
                     afft_l(n*2-1,1) = 0.0d0
                     afft_l(n*2  ,1) = 0.0d0
#endif
                     afft_l(n*2-1,1) = 0.0d0
                     afft_l(n*2  ,1) = 0.0d0
                  end if
               end do
            else if(in == 3) then
               do n = 1, np_fftcd_x
                  i1 = mp_fftcd_x(n)
                  mz = (i1-1)/(lx*ly)+1
                  if ((nnp/2+1) == mz) then
#ifdef FFT_3D_DIVISION_CD
                     mm = mod(i1,(lx*ly))
                     if (mm==0) mm=lx*ly
                     my = (mm-1)/lx+1
                     mx = mod(mm,lx)
                     if (mx==0) mx = lx
                     jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
                     afft_l(jadd*2-1,1) = 0.0d0
                     afft_l(jadd*2  ,1) = 0.0d0
#else
                     afft_l(n*2-1,1) = 0.0d0
                     afft_l(n*2  ,1) = 0.0d0
#endif
                  end if
               end do
            end if
         end if
      end if
                                                 __TIMER_DO_STOP(895)
                                                 __TIMER_SUB_STOP(761)
    end subroutine boundary_zero_into_afft_3D

    subroutine initialize_cggawk_arrays_3D
      grad_trho = 0.d0
#ifndef _XC_SAVE_MEMORY_
      cgrad_rho = 0.d0
#else
      cggawk13 = 0.d0
#endif
      grad_rho = 0.d0
      dF_drho = 0.d0
      dF_dgradrho = 0.d0
    end subroutine initialize_cggawk_arrays_3D

!!$    subroutine cp_chgrhr_onto_afft(is)
!!$      integer, intent(in) :: is
!!$      integer   :: i
!!$
!!$      afft = 0.d0
!!$      do i = ista_fftph, iend_fftph   ! MPI
!!$         afft(2*i-1) = chgrhr_l(i,is)*rinplw
!!$      end do
!!$    end subroutine cp_chgrhr_onto_afft

    subroutine cp_afft_to_cgrad_rho(is,in)
      integer, intent(in) :: is, in
      integer :: i
      do i = ista_fftph, iend_fftph
         cgrad_rho(i,in,is) = -afft(i*2)
      end do
    end subroutine cp_afft_to_cgrad_rho

    subroutine cp_afft_to_cgrad_rho_3D(is,in)
      integer, intent(in) :: is, in
      integer :: i
                                                 __TIMER_SUB_START(762)
                                                 __TIMER_DO_START(896)
#ifdef FFT_3D_DIVISION_CD
         integer :: lx, ly, lz, nxyz, i1, mx, my, mz, mm, jadd, kx1p, kx2p, kx3p
         kx1p = fftcd_X_x_nel
         kx2p = fftcd_X_y_nel
         kx3p = fftcd_X_z_nel
         lx = fft_box_size_CD_3D(1,0)
         ly = fft_box_size_CD_3D(2,0)
         lz = fft_box_size_CD_3D(3,0)
         do i = 1, np_fftcd_x
            i1 = mp_fftcd_x(i)
            mz = (i1-1)/(lx*ly)+1
            mm = mod(i1,(lx*ly))
            if (mm==0) mm=lx*ly
            my = (mm-1)/lx+1
            mx = mod(mm,lx)
            if (mx==0) mx = lx
            jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
            cgrad_rho(i,in,is) = -afft_l(jadd*2,1)
         end do
#else
      do i = 1, nfft_div_size
         cgrad_rho(i,in,is) = -afft_l(i*2,1)
      end do
#endif
                                                 __TIMER_DO_STOP(896)
                                                 __TIMER_SUB_STOP(762)
    end subroutine cp_afft_to_cgrad_rho_3D

    subroutine cp_afft_to_grad2_rho_3D(is)
      integer, intent(in) :: is
      integer :: i
#ifdef FFT_3D_DIVISION_CD
      integer :: lx, ly, lz, nxyz, i1, mx, my, mz, mm, jadd, kx1p, kx2p, kx3p
      kx1p = fftcd_X_x_nel
      kx2p = fftcd_X_y_nel
      kx3p = fftcd_X_z_nel
      lx = fft_box_size_CD_3D(1,0)
      ly = fft_box_size_CD_3D(2,0)
      lz = fft_box_size_CD_3D(3,0)
      do i = 1, np_fftcd_x
         i1 = mp_fftcd_x(i)
         mz = (i1-1)/(lx*ly)+1
         mm = mod(i1,(lx*ly))
         if (mm==0) mm=lx*ly
         my = (mm-1)/lx+1
         mx = mod(mm,lx)
         if (mx==0) mx = lx
         jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
         grad2_rho(i,is) = -afft_l(jadd*2-1,1)
      end do
#else
      do i = 1, nfft_div_size
         grad2_rho(i,is) = -afft_l(i*2-1,1)
      end do
#endif
    end subroutine cp_afft_to_grad2_rho_3D

    subroutine cp_afft_to_chgrhr_3D( iloop )
      integer, intent(in) :: iloop

      integer  :: i
      real(kind=DP) :: tmp
                                                 __TIMER_SUB_START(755)
                                                 __TIMER_DO_START(888)
#ifdef FFT_3D_DIVISION_CD
         integer :: lx, ly, lz, i1, mx, my, mz, mm, jadd, kx1p, kx2p, kx3p
         kx1p = fftcd_X_x_nel
         kx2p = fftcd_X_y_nel
         kx3p = fftcd_X_z_nel
         lx = fft_box_size_CD_3D(1,0)
         ly = fft_box_size_CD_3D(2,0)
         lz = fft_box_size_CD_3D(3,0)
         do i = 1, np_fftcd_x
            i1 = mp_fftcd_x(i)
            mz = (i1-1)/(lx*ly)+1
            mm = mod(i1,(lx*ly))
            if (mm==0) mm=lx*ly
            my = (mm-1)/lx+1
            mx = mod(mm,lx)
            if (mx==0) mx = lx
            jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
            chgrhr_l(i,iloop) = afft_l(jadd*2-1,1)
         end do
#else
      do i = 1, nfft_div_size
         chgrhr_l(i,iloop) = afft_l(i*2-1,1)
      end do
#endif
                                                 __TIMER_DO_STOP(888)
      if(iprixc >= 2) then
         write(nfout,'(" !XC -- chgrhr_l -- <<cp_afft_to_chgrhr>>")')
         write(nfout,'(" !XC ",6f12.6)') (chgrhr_l(i,1),i=1, 18)
         write(nfout,'(" !XC -- afft_l --     <<cp_afft_to_chgrhr>>")')
         write(nfout,'(" !XC ",6f12.6)') (afft_l(i*2-1,1),i=1, 18)
      end if
                                                 __TIMER_SUB_STOP(755)
    end subroutine cp_afft_to_chgrhr_3D

    subroutine wd_small_part_of_vxc
      integer :: i

!!$      if(mype /= 0) return
      if(kimg==1) then
         write(nfout,*) 'vxc'
         if(input_charge == Partial_Core_Charge) then
            do i = ista_kngp, min(ista_kngp+14, iend_kngp)
               write(nfout,200) vxcpc_l(i,1), i
            end do
         else
            if(nspin == 1) then
               do i = ista_kngp, min(ista_kngp+14,iend_kngp)
                  write(nfout,200) vxc_l(i,1,1), i
               end do
            else
               do i = ista_kngp, min(ista_kngp+14, iend_kngp)
                  write(nfout,201) vxc_l(i,1,1),vxc_l(i,1,nspin), i
               end do
            end if
         endif
200      format(' ',d20.12,i8)
201      format(' ','(',d12.4,',',d12.4,')',i8)
         write(nfout,*) 'exc'
         write(nfout,200) exc
      else if(kimg == 2) then
         write(nfout,*) 'vxc'
         if(input_charge == Partial_Core_Charge) then
            do i = ista_kngp, min(ista_kngp+29,iend_kngp)
               write(nfout,203) vxcpc_l(i,1   ),vxcpc_l(i,kimg), i
            end do
         else
            if(ispin == 1) then
               do i = ista_kngp, min(ista_kngp+29,iend_kngp)
                  write(nfout,203) vxc_l(i,1   ,1), vxc_l(i,kimg,1), i
               end do
            else
               do i = ista_kngp, min(ista_kngp+29,iend_kngp)
                  write(nfout,202) vxc_l(i,1   ,1),vxc_l(i,1,   nspin), i,1
                  write(nfout,202) vxc_l(i,kimg,1),vxc_l(i,kimg,nspin), i,2
               end do
            endif
         end if
203      format(' ','(',d12.4,',',d12.4,')',i8)
202      format(' ','(',d12.4,',',d12.4,')',i8, ' kimg= ',i3)
         write(nfout,*) 'exc'
         write(nfout,200) exc
      end if
    end subroutine wd_small_part_of_vxc

    subroutine wd_small_part_of_afft(len_str,str,ne)
      integer, intent(in)        :: len_str
      character(len=len_str),intent(in) :: str
      integer, intent(in)        :: ne

      integer                    :: i
!!$      if(mype == 0) then
         write(nfout, *) 'afft (',str,')'
         if(nspin == 2) write(nfout,*) ' #spin = ', iloop
         if(ne*2 > np_fftp) then
            write(nfout,*)  ' all elements'
            write(nfout,'(6f12.6)') (afft(i),i=ista_fftp, iend_fftp)
         else
            write(nfout,*) ' first ',ne,' elements'
!!$            write(nfout,'(6f12.6)') (afft(i),i=1,ne)
            write(nfout,'(6f12.6)') (afft(ista_fftp-1+i),i=1,ne)
            write(nfout,*) ' last ', ne,' elements'
            write(nfout,'(6f12.6)') (afft(iend_fftp-ne+i),i=1,ne)
         endif
!!$         write(nfout,'(" iend_fftp = ",i8)') iend_fftp
!!$         write(nfout,*) ' last ', lx*2,' elements'
!!$         write(nfout,'(" j, k = ", 2i8)') ly_d+ly_d*myrank_cdfft, lz
!!$         n = 2*((ly_d-1)*lx + (lz-1)*lx*ly_d) + ista_fftp - 1
!!$         write(nfout,'(8d12.4)') (afft( n+i),i=1,lx*2)
!!$         write(nfout,'(" n+lx*2 = ",i8)') n+lx*2
!!$            end do
!!$         end do
!!$         n = 0
!!$         do i = ista_fftp, iend_fftp
!!$            if(abs(afft(i)) > 1.d-8) n = n+1
!!$         end do
!!$         write(nfout,'(" number of non-zero term of afft = ",i8)') n
!!$      end if
    end subroutine wd_small_part_of_afft

#ifdef _MPIFFT_
    subroutine set_cdata(np,aa,a)
      integer, intent(in) :: np
      real(kind=DP), dimension(nfftp), intent(in) :: aa
      real(kind=DP), dimension(mp_fftp), intent(out) :: a
      integer :: i, j, jj, k, nstart, nend, ipin, ipout
      integer :: lx, lxy, llx, lly, llxy

      nstart = np*ny_d+1
      nend   = min(fft_box_size_CD(2,1),(np+1)*ny_d)

      lx   = fft_box_size_CD_c(1,0)
      lxy  = lx*ly_d
      llx  = fft_box_size_CD_c(1,0)
      lly  = fft_box_size_CD_c(2,0)
      llxy = llx*lly

      a = 0.d0
      if(kimg==1)then
         do i = 1, fft_box_size_CD(1,1)
            do j = nstart, nend
               jj = j-nstart+1
               do k = 1, fft_box_size_CD(3,1)
                  ipin  =  i + (j-1)*llx + (k-1)*llxy
                  ipout =  i +(jj-1)* lx + (k-1)*lxy
                  A(ipout) = aa(ipin)
               end do
            end do
         end do
      else if(kimg == 2) then
         do i = 1, fft_box_size_CD(1,1)
            do j = nstart, nend
               jj = j-nstart+1
               do k = 1, fft_box_size_CD(3,1)
                  ipin  = 2*(i + (j-1)*llx + (k-1)*llxy)-1
                  ipout = 2*(i + (jj-1)*lx + (k-1)*lxy)-1
                  A(ipout  ) = aa(ipin)
                  A(ipout+1) = aa(ipin+1)
               end do
            end do
         end do
      end if
    end subroutine set_cdata
#endif

    subroutine set_ispin(ispin)
      integer, intent(out) :: ispin
      if(nspin == 2 .and. input_charge == Valence_plus_PC_Charge) then
         ispin = 2
      else
         ispin = 1
      endif
    end subroutine set_ispin

    subroutine cpafft_to_vxc_or_vxcpc(is)
      integer, intent(in) :: is
      integer :: i, ik, ip
      real(kind=DP),allocatable,dimension(:)    :: afft_mpi1   ! MPI d(nfftp)
      integer :: id_sname = -1
      call tstatc0_begin('cpafft_to_vxc_or_vxcpc ',id_sname)

      allocate(afft_mpi1(nfftp))

      call afft_allgatherv(afft,afft_mpi1)

      do ik = 1, kimg
         if(input_charge == Partial_Core_Charge) then
!            do i = ista_kngp, iend_kngp  ! mpi
            do i = ista_kngp, min(kgp_reduced,iend_kngp)  ! mpi
#ifdef _MPIFFT_
               ip = (igfp_l_c(i)-1)*kimg + ik
#else
               ip = (igfp_l(i)-1)*kimg + ik
#endif
               vxcpc_l(i,ik) = afft_mpi1(ip)*rinplw
            end do
         else if(input_charge == Valence_plus_PC_Charge) then
!            do i = ista_kngp, iend_kngp  ! mpi
            do i = ista_kngp, min(kgp_reduced,iend_kngp)  ! mpi
#ifdef _MPIFFT_
               ip = (igfp_l_c(i)-1)*kimg + ik
#else
               ip = (igfp_l(i)-1)*kimg + ik
#endif
               vxc_l(i,ik,is)   = afft_mpi1(ip)*rinplw
            end do
         end if
      end do
      deallocate(afft_mpi1)
      call tstatc0_end(id_sname)
    end subroutine cpafft_to_vxc_or_vxcpc

    subroutine cpafft(is)
      integer, intent(in) :: is
      integer             :: i

      afft = 0.d0
      do i = ista_fftph, iend_fftph                  ! MPI
         afft(2*i-1) = chgrhr_l(i,is)
      end do

    end subroutine cpafft

    subroutine cpafft_3D(is)
      integer, intent(in) :: is
      integer             :: i
                                                 __TIMER_SUB_START(784)
      afft_l = 0.d0
                                                 __TIMER_DO_START(1020)
#ifdef FFT_3D_DIVISION_CD
         integer :: lx, ly, lz, i1, mx, my, mz, mm, jadd, kx1p, kx2p, kx3p
         kx1p = fftcd_X_x_nel
         kx2p = fftcd_X_y_nel
         kx3p = fftcd_X_z_nel
         lx = fft_box_size_CD_3D(1,0)
         ly = fft_box_size_CD_3D(2,0)
         lz = fft_box_size_CD_3D(3,0)
         do i = 1, np_fftcd_x
            i1 = mp_fftcd_x(i)
            mz = (i1-1)/(lx*ly)+1
            mm = mod(i1,(lx*ly))
            if (mm==0) mm=lx*ly
            my = (mm-1)/lx+1
            mx = mod(mm,lx)
            if (mx==0) mx = lx
            jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
            afft_l(jadd*2-1,1) = chgrhr_l(i,iloop)
         end do
#else
      do i = 1, nfft_div_size                  ! MPI
         afft_l(2*i-1,1) = chgrhr_l(i,is)
      end do
#endif
                                                 __TIMER_DO_STOP(1020)
                                                 __TIMER_SUB_STOP(784)
    end subroutine cpafft_3D

    subroutine cp_vals_to_afft_rspace_3D( is, vals_l, afft_l )
      integer, intent(in) :: is
      real(kind=DP), intent(in) :: vals_l(1:nfft_div_size,ispin)
#ifdef FFT_3D_DIVISION_CD
      real(kind=DP), intent(out) :: afft_l(lsize*2,1)
#else
      real(kind=DP), intent(out) :: afft_l(lsize*kimg,1)
#endif
      integer             :: i
                                                 __TIMER_SUB_START(784)
      afft_l = 0.d0
                                                 __TIMER_DO_START(1020)
#ifdef FFT_3D_DIVISION_CD
         integer :: lx, ly, lz, i1, mx, my, mz, mm, jadd, kx1p, kx2p, kx3p
         kx1p = fftcd_X_x_nel
         kx2p = fftcd_X_y_nel
         kx3p = fftcd_X_z_nel
         lx = fft_box_size_CD_3D(1,0)
         ly = fft_box_size_CD_3D(2,0)
         lz = fft_box_size_CD_3D(3,0)
         do i = 1, np_fftcd_x
            i1 = mp_fftcd_x(i)
            mz = (i1-1)/(lx*ly)+1
            mm = mod(i1,(lx*ly))
            if (mm==0) mm=lx*ly
            my = (mm-1)/lx+1
            mx = mod(mm,lx)
            if (mx==0) mx = lx
            jadd = mx-xyz_fftcd_x(1,1)+1+kx1p*(my-xyz_fftcd_x(1,2))+kx1p*kx2p*(mz-xyz_fftcd_x(1,3))
            afft_l(jadd*2-1,1) = vals_l(i,iloop)
         end do
#else
      do i = 1, nfft_div_size                  ! MPI
         afft_l(2*i-1,1) = vals_l(i,is)
      end do
#endif
                                                 __TIMER_DO_STOP(1020)
                                                 __TIMER_SUB_STOP(784)
    end subroutine cp_vals_to_afft_rspace_3D


!!#ifdef __EDA__
! -----  ascat stars modifying  -----
!!    subroutine cp_exc_for_EDA
!!      use m_Parallelization, only : nel_fftcd_x, nel_fftcd_y, nel_fftcd_z
!!      use m_FFT, only : nfftp_nonpara
!!      integer             :: i, ierr, myoffset, ip, mx, my, mz
!!      real(kind=DP), allocatable, dimension(:,:,:) :: chg,f21
!!      integer, allocatable, dimension(:) :: offset
!!
!!      allocate(offset(0:nrank_g-1)); offset = 0
!!      offset(myrank_g) = nfft_div_size
!!      call mpi_allreduce(mpi_in_place, offset, nrank_g, mpi_integer, mpi_sum, mpi_ke_world,ierr)
!!      myoffset = 0
!!      do i=0,myrank_g-1
!!         myoffset = myoffset+offset(i)
!!      enddo
!!      exc_on_a_grid_wk = univol*exc_on_a_grid_wk
!!      exc_on_a_grid = 0.d0
!!      do i=1,nfft_div_size
!!         call convert_index(i,ip,mx,my,mz)
!!         exc_on_a_grid(2*ip-1) = exc_on_a_grid_wk(i)/f2or1(i)
!!      enddo
!!      deallocate(offset)
!!      call mpi_allreduce(mpi_in_place, exc_on_a_grid, nfftp_nonpara, mpi_double_precision,mpi_sum, &
!!      &                  mpi_ke_world, ierr)
!!
!!    end subroutine cp_exc_for_EDA
! -----  ascat ceases modifying  -----
!!
!!#endif

  end subroutine m_XC_cal_potential_3D

!===============================================================================

  subroutine xcpotf_3D(ispin,input_charge)
    use m_Parallelization,     only : mpi_ke_world
    integer, intent(in) :: ispin,input_charge

! #1) 1994/11/08 by T.Yamasaki
!    Coding for the case of xctype='PERZUN ' and  'XALPHA ' are done.
! #2) Spin-polarization is introduced by T. Yamasaki at 15th Dec. 1994
! #3) f77 -> f90     4th April 1999  by T. Yamasaki

    real(kind=DP) :: rinplw
    real(kind=DP) :: exc_mpi     ! MPI

    real(kind=DP) :: DELTA
    data DELTA/1.d-40/
                                                 __TIMER_SUB_START(777)

    rinplw = 1.d0/product(fft_box_size_CD_3D(1:3,1))

    exc = 0.d0
    if(xctype == 'wign   '.or. xctype == 'wigner ')then
       call xcpotf_wigner_3D(nspin,ispin,input_charge,DELTA,chgrhr_l,f2or1,exc,nfft_div_size)
    else if(  xctype ==  'pzold  ') then
       call xcpotf_pzold_3D(nspin,ispin,input_charge,DELTA,chgrhr_l,f2or1,exc,nfft_div_size)
    else if(  xctype ==  'xalfa  ') then
       call xcpotf_xalfa_3D(nspin,ispin,input_charge,DELTA,chgrhr_l,f2or1,exc,nfft_div_size)
    else if(  xctype == 'perzun '.or. xctype == 'pz     ') then
       call xcpotf_pz_3D(nspin,ispin,input_charge,DELTA,chgrhr_l,f2or1,exc,nfft_div_size)
    else if(  xctype == 'vwn    ') then
       call xcpotf_vwn_3D(nspin,ispin,input_charge,DELTA,chgrhr_l,f2or1,exc,nfft_div_size)
    else if(  xctype=='mjw    '.or. xctype=='bh     ' .or.xctype=='gl     ') then
       call xcpotf_mjw_bh_gl_3D(len_xctype,xctype,nspin,ispin,input_charge,DELTA,chgrhr_l,f2or1,exc,nfft_div_size)
    endif
                        ! all xcpotf_* subroutines are in -(b_XC_Potential) ->chgrhr_l,exc
    exc = exc*univol*rinplw

#ifdef CD_FFT_ALL
    if(npes >= 2) then
                                                 __TIMER_COMM_START_w_BARRIER(MPI_CommGroup,1011)
       call mpi_allreduce(exc,exc_mpi,1,mpi_double_precision,mpi_sum, MPI_CommGroup,ierr)
                                                 __TIMER_COMM_STOP(1011)
       exc = exc_mpi
    end if
#else
    if(nrank_g >= 2) then
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,1012)
       call mpi_allreduce(exc,exc_mpi,1,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
                                                 __TIMER_COMM_STOP(1012)
       exc = exc_mpi
    end if
#endif
                                                 __TIMER_SUB_STOP(777)
  end subroutine xcpotf_3D

!===============================================================================

  subroutine check_of_negative_CD_3D(nfout, ispin)
    use m_Parallelization,      only : np_fftcd_y, mp_fftcd_y   &
   &                                 , np_fftcd_x, mp_fftcd_x   &
   &                                 , np_fftcd_z, mp_fftcd_z   &
   &                                 , xyz_fftcd_y, xyz_fftcd_x, xyz_fftcd_z &
   &                                 , mpi_ke_world
    integer, intent(in) :: nfout, ispin
!XX!integer, allocatable, dimension(:) :: wk_mp_fft_y
    integer  :: nx, ny, nz, i, j, k, ierr
    integer  :: lx, ly, lz, kx, ky, kz, loop
    integer  :: mx, my, mz, i1, mm, n, ip, isw
    integer  :: ipadx, ipady, ipadz
    integer  :: icwarn_NEGA, icwarn_IMAG, icwarn_NEGA_total, icwarn_IMAG_total
    real(kind=DP),parameter   :: chgdel = 8.d-5
    real(kind=DP),parameter   :: D_min  = 1.d-40

                                                 __TIMER_SUB_START(754)
    if (nfft_div_size > 0) then
#ifdef FFT_3D_DIVISION_CD
       nx = xyz_fftcd_x(2,1)-xyz_fftcd_x(1,1)+1
       ny = xyz_fftcd_x(2,2)-xyz_fftcd_x(1,2)+1
       nz = xyz_fftcd_x(2,3)-xyz_fftcd_x(1,3)+1
#else
       if (sw_fft_xzy > 0) then
          nx = xyz_fftcd_y(2,1)-xyz_fftcd_y(1,1)+1
          ny = xyz_fftcd_y(2,2)-xyz_fftcd_y(1,2)+1
          nz = xyz_fftcd_y(2,3)-xyz_fftcd_y(1,3)+1
       else
          nx = xyz_fftcd_z(2,1)-xyz_fftcd_z(1,1)+1
          ny = xyz_fftcd_z(2,2)-xyz_fftcd_z(1,2)+1
          nz = xyz_fftcd_z(2,3)-xyz_fftcd_z(1,3)+1
       end if
#endif

!XX!   allocate(wk_mp_fft_y(nx*ny*nz) ,stat=ierr)
!XX!   do k = 0, nz-1
!XX!      do j = 0, ny-1
!XX!            do i = 0, nx-1
!XX!            wk_mp_fft_y(i+j*nx+k*nx*ny+1) = mp_fftcd_y(i+k*nx+j*nx*nz+1)
!XX!         end do
!XX!      end do
!XX!   end do
    end if

    lx = fft_box_size_CD_3D(1,0)
    ly = fft_box_size_CD_3D(2,0)
    lz = fft_box_size_CD_3D(3,0)
    kx  = fft_box_size_CD_3D(1,1)
    ky  = fft_box_size_CD_3D(2,1)
    kz  = fft_box_size_CD_3D(3,1)

                                                 __TIMER_DO_START(884)
    if (kimg == 1) then
       ipadx = fft_box_size_CD_3D(1,0) - fft_box_size_CD_3D(1,1)
       ipady = fft_box_size_CD_3D(2,0) - fft_box_size_CD_3D(2,1)
       ipadz = fft_box_size_CD_3D(3,0) - fft_box_size_CD_3D(3,1)
       if ((ipadx >= 4) .or. (ipady >= 2) .or. (ipadz >= 2)) then
          if (sw_fft_xzy > 0) then
             loop = np_fftcd_y
          else
             loop = np_fftcd_z
          end if
          do n = 1, loop , 2
             if (sw_fft_xzy > 0) then
                i1 = mp_fftcd_y(n)
             else
                i1 = mp_fftcd_z(n)
             end if
             mz = (i1-1)/(lx*ly)+1
             mm = mod(i1,(lx*ly))
             if (mm==0) mm=lx*ly
             my = (mm-1)/lx+1
             mx = mod(mm,lx)
             if (mx==0) mx = lx

             if (ipadx >= 4) then
                if ((kx+2) < mx) then
                   afft_l(n,1) = 0.d0
                end if
             end if
             if (ipady >= 2) then
                if (ky < my) then
                   afft_l(n,1) = 0.d0
                end if
             end if
             if (ipadz >= 2) then
                if (kz < mz) then
                   afft_l(n,1) = 0.d0
                   cycle
                end if
             end if
          end do
       end if
    else
#ifdef FFT_3D_DIVISION_CD
       do n = 1, np_fftcd_x
          i1 = mp_fftcd_x(n)
          mz = (i1-1)/(lx*ly)+1
          mm = mod(i1,(lx*ly))
          if (mm==0) mm=lx*ly
          my = (mm-1)/lx+1
          mx = mod(mm,lx)
          if (mx==0) mx = lx

          if (kx < mx) then
             afft_l(n*2-1,1) = 0.d0
             afft_l(n*2  ,1) = 0.d0
             cycle
          end if
          if (ky < my) then
             afft_l(n*2-1,1) = 0.d0
             afft_l(n*2  ,1) = 0.d0
            cycle
          end if
          if (kz < mz) then
             afft_l(n*2-1,1) = 0.d0
             afft_l(n*2  ,1) = 0.d0
             cycle
          end if
       end do
#else
       if (sw_fft_xzy > 0) then
          loop = np_fftcd_y
       else
          loop = np_fftcd_z
       end if
       do n = 1, loop
          if (sw_fft_xzy > 0) then
             i1 = mp_fftcd_y(n)
          else
             i1 = mp_fftcd_z(n)
          end if
          mz = (i1-1)/(lx*ly)+1
          mm = mod(i1,(lx*ly))
          if (mm==0) mm=lx*ly
          my = (mm-1)/lx+1
          mx = mod(mm,lx)
          if (mx==0) mx = lx

          if (kx < mx) then
             afft_l(n*2-1,1) = 0.d0
             afft_l(n*2  ,1) = 0.d0
             cycle
          end if
          if (ky < my) then
             afft_l(n*2-1,1) = 0.d0
             afft_l(n*2  ,1) = 0.d0
            cycle
          end if
          if (kz < mz) then
             afft_l(n*2-1,1) = 0.d0
             afft_l(n*2  ,1) = 0.d0
             cycle
          end if
       end do
#endif
    end if
                                                 __TIMER_DO_STOP(884)

    icwarn_NEGA = 0
    icwarn_IMAG = 0
    if(iprinegativecharge0 < 1 .or. max_warnings_negativecharge <= 0) then
                                                 __TIMER_DO_START(885)
       if (kimg == 1) then
          if (sw_fft_xzy > 0) then
             loop = np_fftcd_y
          else
             loop = np_fftcd_z
          end if
          do n = 1, loop
             if (sw_fft_xzy > 0) then
                i1 = mp_fftcd_y(n)
             else
                i1 = mp_fftcd_z(n)
             end if
             mz = (i1-1)/(lx*ly)+1
             mm = mod(i1,(lx*ly))
             if (mm==0) mm=lx*ly
             my = (mm-1)/lx+1
             mx = mod(mm,lx)
             if (mx==0) mx = lx

             if (mod(n,2) == 1) then
                if ((afft_l(n,1) <= -chgdel) .and. (mod(mx+1,lx) /= 0)) then
                   icwarn_NEGA = icwarn_NEGA + 1
                   afft_l(n,1) = D_min
                else if (afft_l(n,1) <= 0.0d0) then
                   afft_l(n,1) = D_min
                end if
             else
                if ((abs(afft_l(n,1)) > chgdel) .and. (mod(mx+1,lx) /= 0)) then
                   icwarn_IMAG = icwarn_IMAG + 1
                   afft_l(n,1) = 0.0d0
                else if (abs(afft_l(n,1)) > 0.0d0) then
                   afft_l(n,1) = 0.0d0
                end if
             end if
          end do
       else
#ifdef FFT_3D_DIVISION_CD
          do n = 1, np_fftcd_x
             i1 = mp_fftcd_x(n)
             mz = (i1-1)/(lx*ly)+1
             mm = mod(i1,(lx*ly))
             if (mm==0) mm=lx*ly
             my = (mm-1)/lx+1
             mx = mod(mm,lx)
             if (mx==0) mx = lx

             ip = n*2-1
             if ((afft_l(ip,1) <= -chgdel) .and. (mod(mx+1,lx) /= 0)) then
                icwarn_NEGA = icwarn_NEGA + 1
                afft_l(ip,1) = D_min
             else if (afft_l(ip,1) <= 0.0d0) then
                afft_l(ip,1) = D_min
             end if
             if ((abs(afft_l(ip+1,1)) > chgdel) .and. (mod(mx+1,lx) /= 0)) then
                icwarn_IMAG = icwarn_IMAG + 1
                afft_l(ip+1,1) = 0.0d0
             else if (abs(afft_l(ip+1,1)) > 0.0d0) then
                afft_l(ip+1,1) = 0.0d0
             end if
          end do
#else
          if (sw_fft_xzy > 0) then
             loop = np_fftcd_y
          else
             loop = np_fftcd_z
          end if
          do n = 1, loop
             if (sw_fft_xzy > 0) then
                i1 = mp_fftcd_y(n)
             else
                i1 = mp_fftcd_z(n)
             end if
             mz = (i1-1)/(lx*ly)+1
             mm = mod(i1,(lx*ly))
             if (mm==0) mm=lx*ly
             my = (mm-1)/lx+1
             mx = mod(mm,lx)
             if (mx==0) mx = lx

             ip = n*2-1
             if ((afft_l(ip,1) <= -chgdel) .and. (mod(mx+1,lx) /= 0)) then
                icwarn_NEGA = icwarn_NEGA + 1
                afft_l(ip,1) = D_min
             else if (afft_l(ip,1) <= 0.0d0) then
                afft_l(ip,1) = D_min
             end if
             if ((abs(afft_l(ip+1,1)) > chgdel) .and. (mod(mx+1,lx) /= 0)) then
                icwarn_IMAG = icwarn_IMAG + 1
                afft_l(ip+1,1) = 0.0d0
             else if (abs(afft_l(ip+1,1)) > 0.0d0) then
                afft_l(ip+1,1) = 0.0d0
             end if
          end do
#endif
       endif
                                                 __TIMER_DO_STOP(885)
    else if(iprinegativecharge0 >= 1 .and. max_warnings_negativecharge > 0) then
                                                 __TIMER_DO_START(886)
       isw = 0
       if (kimg == 1) then
          if (sw_fft_xzy > 0) then
             loop = np_fftcd_y
          else
             loop = np_fftcd_z
          end if
          do n = 1, loop
             if (sw_fft_xzy > 0) then
                i1 = mp_fftcd_y(n)
             else
                i1 = mp_fftcd_z(n)
             end if
             mz = (i1-1)/(lx*ly)+1
             mm = mod(i1,(lx*ly))
             if (mm==0) mm=lx*ly
             my = (mm-1)/lx+1
             mx = mod(mm,lx)
             if (mx==0) mx = lx

             if (mod(n,2) == 1) then
                if ((afft_l(n,1) <= -chgdel) .and. (mod(mx+1,lx) /= 0)) then
                   icwarn_NEGA = icwarn_NEGA + 1
                   if(icwarn_NEGA <= max_warnings_negativecharge) then
                      if(icwarn_NEGA == 1 .and. nspin == 2) then
                         if(iprinegativecharge>=1) then
                            write(nfout,'(" #spin = ",i3)') ispin
                            isw = 1
                         end if
                      endif
                      if(iprinegativecharge >=1) then
! === DEBUG by tkato 2011/11/19 ================================================
!                        write(nfout,'(" *** WARN CHG.DEN = ",d15.7, " < 0.0 AT ",i8," ***")') afft_l(n,1),i
                         write(nfout,'(" *** WARN CHG.DEN = ",d15.7, " < 0.0 AT ",i8," ***")') afft_l(n,1),n
! ==============================================================================
                      end if
                   end if
                   afft_l(n,1) = D_min
                else if (afft_l(n,1) <= 0.0d0) then
                   afft_l(n,1) = D_min
                end if
             else
                if ((abs(afft_l(n,1)) > chgdel) .and. (mod(mx+1,lx) /= 0)) then
                   icwarn_IMAG = icwarn_IMAG + 1
                   if(icwarn_IMAG <= max_warnings_negativecharge) then
                      if(icwarn_IMAG == 1 .and. nspin == 2) then
                         if(iprinegativecharge >=1 .and. isw == 0) then
                            write(nfout,'(" #spin = ",i3)') ispin
                         end if
                      end if
                      if(iprinegativecharge>=1) then
! === DEBUG by tkato 2011/11/19 ================================================
!                        write(nfout,'(" *** WARN IMAG(CHG) = ",d15.7, " > 0.0 AT ",i8," ***")') afft_l(n,1),i
                         write(nfout,'(" *** WARN IMAG(CHG) = ",d15.7, " > 0.0 AT ",i8," ***")') afft_l(n,1),n
! ==============================================================================
                      end if
                   end if
                   afft_l(n,1) = 0.0d0
                else if (abs(afft_l(n,1)) > 0.0d0) then
                   afft_l(n,1) = 0.0d0
                end if
             end if
          end do
       else
#ifdef FFT_3D_DIVISION_CD
          do n = 1, np_fftcd_x
             i1 = mp_fftcd_x(n)
             mz = (i1-1)/(lx*ly)+1
             mm = mod(i1,(lx*ly))
             if (mm==0) mm=lx*ly
             my = (mm-1)/lx+1
             mx = mod(mm,lx)
             if (mx==0) mx = lx

             ip = n*2-1
             if ((afft_l(ip,1) <= -chgdel) .and. (mod(mx+1,lx) /= 0)) then
                icwarn_NEGA = icwarn_NEGA + 1
                if(icwarn_NEGA <= max_warnings_negativecharge) then
                   if(icwarn_NEGA == 1 .and. nspin == 2) then
                      if(iprinegativecharge>=1) then
                         write(nfout,'(" #spin = ",i3)') ispin
                         isw = 1
                      end if
                   endif
                   if(iprinegativecharge >=1) then
! === DEBUG by tkato 2011/11/19 ================================================
!                     write(nfout,'(" *** WARN CHG.DEN = ",d15.7, " < 0.0 AT ",i8," ***")') afft_l(ip,1),i
                      write(nfout,'(" *** WARN CHG.DEN = ",d15.7, " < 0.0 AT ",i8," ***")') afft_l(ip,1),ip
! ==============================================================================
                   end if
                end if
                afft_l(ip,1) = D_min
             else if (afft_l(ip,1) <= 0.0d0) then
                afft_l(ip,1) = D_min
             end if
             if ((abs(afft_l(ip+1,1)) > chgdel) .and. (mod(mx+1,lx) /= 0)) then
                icwarn_IMAG = icwarn_IMAG + 1
                if(icwarn_IMAG <= max_warnings_negativecharge) then
                   if(icwarn_IMAG == 1 .and. nspin == 2) then
                      if(iprinegativecharge >=1 .and. isw == 0) then
                         write(nfout,'(" #spin = ",i3)') ispin
                      end if
                   end if
                   if(iprinegativecharge>=1) then
! === DEBUG by tkato 2011/11/19 ================================================
!                     write(nfout,'(" *** WARN IMAG(CHG) = ",d15.7, " > 0.0 AT ",i8," ***")') afft_l(ip+1,1),i
                      write(nfout,'(" *** WARN IMAG(CHG) = ",d15.7, " > 0.0 AT ",i8," ***")') &
                           & afft_l(ip+1,1),ip+1
! ==============================================================================
                   end if
                end if
                afft_l(ip+1,1) = 0.0d0
             else if (abs(afft_l(ip+1,1)) > 0.0d0) then
                afft_l(ip+1,1) = 0.0d0
             end if
          end do
#else
          if (sw_fft_xzy > 0) then
             loop = np_fftcd_y
          else
             loop = np_fftcd_z
          end if
          do n = 1, loop
             if (sw_fft_xzy > 0) then
                i1 = mp_fftcd_y(n)
             else
                i1 = mp_fftcd_z(n)
             end if
             mz = (i1-1)/(lx*ly)+1
             mm = mod(i1,(lx*ly))
             if (mm==0) mm=lx*ly
             my = (mm-1)/lx+1
             mx = mod(mm,lx)
             if (mx==0) mx = lx

             ip = n*2-1
             if ((afft_l(ip,1) <= -chgdel) .and. (mod(mx+1,lx) /= 0)) then
                icwarn_NEGA = icwarn_NEGA + 1
                if(icwarn_NEGA <= max_warnings_negativecharge) then
                   if(icwarn_NEGA == 1 .and. nspin == 2) then
                      if(iprinegativecharge>=1) then
                         write(nfout,'(" #spin = ",i3)') ispin
                         isw = 1
                      end if
                   endif
                   if(iprinegativecharge >=1) then
! === DEBUG by tkato 2011/11/19 ================================================
!                     write(nfout,'(" *** WARN CHG.DEN = ",d15.7, " < 0.0 AT ",i8," ***")') afft_l(ip,1),i
                      write(nfout,'(" *** WARN CHG.DEN = ",d15.7, " < 0.0 AT ",i8," ***")') afft_l(ip,1),ip
! ==============================================================================
                   end if
                end if
                afft_l(ip,1) = D_min
             else if (afft_l(ip,1) <= 0.0d0) then
                afft_l(ip,1) = D_min
             end if
             if ((abs(afft_l(ip+1,1)) > chgdel) .and. (mod(mx+1,lx) /= 0)) then
                icwarn_IMAG = icwarn_IMAG + 1
                if(icwarn_IMAG <= max_warnings_negativecharge) then
                   if(icwarn_IMAG == 1 .and. nspin == 2) then
                      if(iprinegativecharge >=1 .and. isw == 0) then
                         write(nfout,'(" #spin = ",i3)') ispin
                      end if
                   end if
                   if(iprinegativecharge>=1) then
! === DEBUG by tkato 2011/11/19 ================================================
!                     write(nfout,'(" *** WARN IMAG(CHG) = ",d15.7, " > 0.0 AT ",i8," ***")') afft_l(ip+1,1),i
                      write(nfout,'(" *** WARN IMAG(CHG) = ",d15.7, " > 0.0 AT ",i8," ***")') afft_l(ip+1,1),ip+1
! ==============================================================================
                   end if
                end if
                afft_l(ip+1,1) = 0.0d0
             else if (abs(afft_l(ip+1,1)) > 0.0d0) then
                afft_l(ip+1,1) = 0.0d0
             end if
          end do
#endif
       endif
                                                 __TIMER_DO_STOP(886)
    endif

    if(iprinegativecharge0 >= 1) then
       if(nrank_g > 1) then
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,887)
          call mpi_allreduce(icwarn_NEGA,icwarn_NEGA_total,1,mpi_integer,mpi_sum,mpi_ke_world,ierr)
          call mpi_allreduce(icwarn_IMAG,icwarn_IMAG_total,1,mpi_integer,mpi_sum,mpi_ke_world,ierr)
                                                 __TIMER_COMM_STOP(887)
       else
          icwarn_NEGA_total = icwarn_NEGA
          icwarn_IMAG_total = icwarn_IMAG
       end if
    end if

    if(iprinegativecharge >=1) then
       if(icwarn_NEGA_total >= 1) then
          if(nspin == 2) write(nfout,'(" #spin = ",i3," : 1 = UP, 2 = DOWN")') ispin
          if(nrank_g > 1) then
             write(nfout,'(" *** WARN  # of <<Negative Charge Density>>  = ",i9)') icwarn_NEGA_total
          else
             write(nfout,'(" *** WARN  # of <<Negative Charge Density>>  = ",i9,", (node ",i3,") = ",i9)') &
           & icwarn_NEGA_total, myrank_g,icwarn_NEGA
          end if
       endif

       if(icwarn_IMAG_total >= 1) then
          if(nspin == 2) write(nfout,'(" #spin = ",i3," : 1 = UP, 2 = DOWN")') ispin
          if(nrank_g > 1) then
             write(nfout,'(" *** WARN  # of <<Imaginary Charge Density>> = ",i9)') icwarn_NEGA_total
          else
             write(nfout,'(" *** WARN  # of <<Imaginary Charge Density>> = ",i9,", (node ",i3,") = ",i9)') &
           & icwarn_NEGA_total, myrank_g,icwarn_NEGA
          end if
       end if
    end if

!XX!if (nfft_div_size > 0) then
!XX!   deallocate(wk_mp_fft_y)
!XX!end if

                                                 __TIMER_SUB_STOP(754)
  end subroutine check_of_negative_CD_3D

  subroutine m_XC_rst_npsend()
     np_send_is_set = .false.
     if(allocated(ip)) deallocate(ip)
     if(allocated(igfp_ijk)) deallocate(igfp_ijk)
     if(allocated(np_send)) deallocate(np_send)
  end subroutine m_XC_rst_npsend

end module m_XC_Potential
