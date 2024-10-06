!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 573 $)
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
!module m_XC_Potential

module m_XC_Variables
  use m_Const_Parameters, only : DP
  implicit none
  real(kind=DP), allocatable, dimension(:,:,:) :: vxc_l    ! d(ista_kngp:iend_kngp,kimg,nspin)  MPI
  real(kind=DP), allocatable, dimension(:,:,:) :: vxco_l   ! d(ista_kngp:iend_kngp,kimg,nspin)  MPI
  real(kind=DP), allocatable, dimension(:,:)   :: vxcpc_l  ! d(ista_kngp:iend_kngp,kimg)        MPI
  real(kind=DP)                                :: exc, excpc
  real(kind=DP)                                :: eex, ecor
end module m_XC_Variables

module m_XC_Potential_2D
! $Id: m_XC_Potential.F90 573 2017-05-08 07:17:02Z jkoga $
!
!  Upgraded on 23rd Aug. 2006 by T. Yamasaki
!    Differentials of the charge density function in GGA calculation are
!    parallelized according to the components of x, y, and z.
!  Upgraded on 19th Sep. 2006 by T. Yamasaki
!    FFTs in GGA calculations are parallelized.
!
  use m_PlaneWaveBasisSet,    only : ngabc,gr_l,igfp_l,kg,kgp,kgp_reduced, ylm_l&
       &                           , m_pwBS_sphrp2_3D
  use m_PlaneWaveBasisSet,    only : m_pwBS_sphrp2_diff_3D
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
       &                           , m_FFT_CD_inverse_c &
       &                           , m_FFT_CD_direct_c  &
       &                           , m_FFT_check_of_negative_CD, fft_box_size_CD_nonpara
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
#ifdef ENABLE_ESM_PACK
       &                           , max_warnings_negativecharge,sw_esm
#else
       &                           , max_warnings_negativecharge
#endif
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
       &                               m_KE_set_modeled_ekin_density, &
       &                               m_KE_add_ekindens_hard_rspace
! =================== 13.0XX. 2014/09/19

! === KT_add === 2014/08/04
  use m_Crystal_Structure,  only : sw_neglect_low_vorticity, sw_neglect_low_helicity, &
       &                           sw_monitor_magnetic_vorticity
  use m_ES_Noncollinear,   only : m_ES_add_contrib_to_Vorticity, &
       &                          m_ES_neglect_low_Helicity, &
       &                          m_ES_neglect_low_Vorticity
! ============== 2014/08/04
  use m_XC_Variables,      only : vxc_l, vxco_l, vxcpc_l, exc, excpc, eex, ecor

  use m_Control_Parameters,  only : sw_opencore, sw_xc_opencore_ae_only
  use m_PS_opencore,  only : has_opencore, rmag_opencore_l, mag_opencore_pol

  use m_IterationNumbers,   only : iteration

  use m_Control_Parameters,  only : vtau_exists
  use m_Electronic_Structure,  only : vtau_l
#ifdef LIBXC
  use m_Control_Parameters,  only : xc_family_exch, xc_family_corr, xc_name_exch
  use xc_f03_lib_m
#endif

#ifdef __EDA__
  use m_Control_Parameters, only : sw_eda
#endif
  use mpi
  implicit none

!  include 'mpif.h'
  integer istatus(mpi_status_size)

#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),public,target,allocatable,dimension(:) :: exc_on_a_grid
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
     &                    ,dFx_drho,dFx_dgradrho, exc_on_a_grid_wk, ist,ien)
#else
  subroutine ex_ggapw91(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc &
     &                    ,dFx_drho,dFx_dgradrho, ist,ien)
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
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif
  end subroutine ex_ggapw91

#ifdef __EDA__
  subroutine cr_ggapw91(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_trho,wos,exc,dF_drho,exc_on_a_grid_wk, ist,ien)
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
  subroutine ex_ggapbe(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dFx_drho,dFx_dgradrho,exc_on_a_grid_wk, revPBE,ist,ien)
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
  subroutine ex_ggapbe(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dFx_drho,dFx_dgradrho,exc_on_a_grid_wk,ien)
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
  subroutine cr_ggapbe(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_trho,wos,exc,dF_drho,exc_on_a_grid_wk,ecor,ist,ien)
#else
  subroutine cr_ggapbe(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_trho,wos,exc,dF_drho,ecor,ist,ien)
#endif
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
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
  subroutine xclda(nspin,ispin,ista_r,iend_r,chgrhr_l,wos,exc,dF_drho,exc_on_a_grid_wk,ien)
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
  subroutine ggabek(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dF_drho,dF_dgradrho,exc_on_a_grid_wk,ien)
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
  subroutine ggaprd(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dF_drho,dF_dgradrho,exc_on_a_grid_wk,ien)
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
  subroutine cr_lda(nspin,ispin,ista_r,iend_r,chgrhr_l,exc,dF_drho,exc_on_a_grid_wk,ien)
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
#ifdef __EDA__
  subroutine ex_ggapbe_3D(nspin,ispin,chgrhr_l,grad_rho,f2or1,exc,dFx_drho,dFx_dgradrho,nfft_y,iteration,exc_on_a_grid_wk,revPBE)
#else
  subroutine ex_ggapbe_3D(nspin,ispin,chgrhr_l,grad_rho,f2or1,exc,dFx_drho,dFx_dgradrho,nfft_y,iteration,revPBE)
#endif
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
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(1:nfft_y)
! -----  ascat ceases modifying  -----
#endif
  logical, intent(in),optional :: revPBE
  end subroutine ex_ggapbe_3D
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

  subroutine m_XC_alloc_vxc
! ============================== modified by K. Tagami =========== 11.0
!    allocate(vxc_l(ista_kngp:iend_kngp,kimg,nspin)); vxc_l = 0.d0
!    allocate(vxco_l(ista_kngp:iend_kngp,kimg,nspin)); vxco_l = 0.d0
!    allocate(vxcpc_l(ista_kngp:iend_kngp,kimg)); vxcpc_l = 0.d0

    if ( noncol ) then
       allocate(vxc_l(ista_kngp:iend_kngp,kimg,ndim_chgpot)); vxc_l = 0.d0
       allocate(vxco_l(ista_kngp:iend_kngp,kimg,ndim_chgpot)); vxco_l = 0.d0
       allocate(vxcpc_l(ista_kngp:iend_kngp,kimg)); vxcpc_l = 0.d0
    else
       allocate(vxc_l(ista_kngp:iend_kngp,kimg,nspin)); vxc_l = 0.d0
       allocate(vxco_l(ista_kngp:iend_kngp,kimg,nspin)); vxco_l = 0.d0
       allocate(vxcpc_l(ista_kngp:iend_kngp,kimg)); vxcpc_l = 0.d0
    endif
! ==================================================================== 11.0
  end subroutine m_XC_alloc_vxc

  subroutine m_XC_dealloc_vxc
!f    deallocate(vxc_l)
!f    deallocate(vxco_l)
!f    deallocate(vxcpc_l)
    if(allocated(vxc_l))     deallocate(vxc_l)
    if(allocated(vxco_l))    deallocate(vxco_l)
    if(allocated(vxcpc_l))   deallocate(vxcpc_l)
  end subroutine m_XC_dealloc_vxc

#ifdef __EDA__
! -----  ascat starts modifying  -----
  subroutine m_XC_alloc_exc_on_a_grid()
    use m_FFT, only : nfftp_nonpara
!    allocate(exc_on_a_grid(ista_fftp:iend_fftp)); exc_on_a_grid = 0.d0
!    allocate(exc_on_a_grid(nfftp_nonpara)); exc_on_a_grid = 0.d0
    allocate(exc_on_a_grid(nfftp_nonpara)); exc_on_a_grid = 0.d0
  end subroutine m_XC_alloc_exc_on_a_grid

  subroutine m_XC_dealloc_exc_on_a_grid()
    deallocate(exc_on_a_grid)
  end subroutine m_XC_dealloc_exc_on_a_grid
! -----  ascat ceases modifying  -----
#endif

  subroutine cp_vxc_to_vxco
    vxco_l = vxc_l
  end subroutine cp_vxc_to_vxco

  subroutine xc_allocate(ispin,nfout)
!#ifdef PARA3D
!    use m_Parallelization,      only : np_fftcd_y
!#endif
    integer, intent(in)      :: ispin,nfout

    integer                  :: idp,nlp,nmp,nnp,nd2p,nd3p,ip, idph, nlph
    integer                  :: n, nn, i, j, k
#ifdef _MPIFFT_
    integer                  :: np0, j0, i2, j2, k2, n0
#endif

#ifndef _MPIFFT_
    call m_FFT_alloc_CD_box()
#endif
    allocate(afft(ista_fftp:iend_fftp));afft=0.d0
    allocate(chgrhr_l(ista_fftph:iend_fftph,ispin))    ! MPI
    allocate(f2or1(ista_fftph:iend_fftph))

    nlp  = fft_box_size_CD(1,1)
    nmp  = fft_box_size_CD(2,1)
    nnp  = fft_box_size_CD(3,1)
#ifdef _MPIFFT_
    idp  = fft_box_size_CD_c(1,0)
    nd2p = fft_box_size_CD_c(2,0)
    nd3p = fft_box_size_CD_c(3,0)
#else
    idp  = fft_box_size_CD(1,0)
    nd2p = fft_box_size_CD(2,0)
    nd3p = fft_box_size_CD(3,0)
#endif

    call set_f2or1(npes_cdfft,ista_fftph,iend_fftph,f2or1)

    if(check_of_xctype() == GGA) then
       allocate(inx(ista_fftp:iend_fftp)); inx = 0
       allocate(jnx(ista_fftp:iend_fftp)); jnx = 0
       allocate(knx(ista_fftp:iend_fftp)); knx = 0

#ifdef _MPIFFTTESTTEST_
       if(kimg == 1) then
          do n = ista_fftp, iend_fftp
             nn = n
             if(nn == ista_fftp-1) write(nfout,'(" n == ista_fftp-1 = ",i8)') nn
             np0 = nn - idp*ly_d*lz*myrank_cdfft
             i  = mod(np0-1,idp)+1
             j0 = mod((np0-i)/idp,ly_d) + 1
             j  = j0 + ny_d*myrank_cdfft
             k  = (np0-i-(j0-1)*idp)/(idp*ly_d) + 1
             if(i >= idp) i = i - idp
             inx(n) = i - 1;  if(2*inx(n) > nlp) inx(n) = inx(n) - nlp
             jnx(n) = j - 1;  if(2*jnx(n) > nmp) jnx(n) = jnx(n) - nmp
             knx(n) = k - 1;  if(2*knx(n) > nnp) knx(n) = knx(n) - nnp
          end do

!!$          do k = 0, nnp-1
!!$             k2 = k; if(2*k2 > nnp) k2 = k2 - nnp
!!$             do j0 = 0, min(ny_d, nmp-ny_d*myrank_cdfft)-1
!!$                j = j0 + ny_d*myrank_cdfft
!!$                j2 = j; if(2*j2 > nmp) j2 = j2 - nmp
!!$                n0 = 1+j0*idp + k*idp*ly_d + ista_fftp-1
!!$                do i = 0, nlp-1
!!$                   n = i+n0
!!$                   inx(n) = i; if(2*inx(n) > nlp) inx(n) = inx(n) - nlp
!!$                   if(n == ista_fftp-1) write(nfout,'(" n = ista_fftp-1 = ",i8)') n
!!$                   jnx(n) = j2
!!$                   knx(n) = k2
!!$                end do
!!$             end do
!!$          end do

!!$          do i = 0, nlp-1
!!$             i2 = i; if(2*i2 > nlp) i2 = i2 - nlp
!!$             do j0 = 0, min(ny_d, nmp-ny_d*myrank_cdfft)-1
!!$                j = j0 + ny_d*myrank_cdfft
!!$                j2 = j; if(2*j2 > nmp) j2 = j2 - nmp
!!$                do k = 0, nnp/2
!!$                   n = i+1 + j0*idp + k*idp*ly_d + ista_fftp-1
!!$                   inx(n) = i2
!!$                   jnx(n) = j2
!!$                   knx(n) = k
!!$                end do
!!$                do k = nnp/2+1, nnp-1
!!$                   n = i+1 + j0*idp + k*idp*ly_d + ista_fftp-1
!!$                   inx(n) = i2
!!$                   jnx(n) = j2
!!$                   knx(n) = k - nnp
!!$                end do
!!$             end do
!!$          end do
       else if(kimg == 2) then
          do k = 0, nnp-1
             k2 = k; if(2*k2 > nnp) k2 = k2 - nnp
             do j0 = 0, min(ny_d, nmp-ny_d*myrank_cdfft)-1
                j = j0 + ny_d*myrank_cdfft
                j2 = j; if(2*j2 > nmp) j2 = j2 - nmp
                do i = 0, nlp/2
                   n = i+1 + j0*idp + k*idp*ly_d
                   nn = n*2-1 + ista_fftp-1
                   inx(nn) = i ; inx(nn+1) = i
                   jnx(nn) = j2; jnx(nn+1) = j2
                   knx(nn) = k2; knx(nn+1) = k2
                end do
                do i = nlp/2+1, nlp-1
                   n = i+1 + j0*idp + k*idp*ly_d
                   nn = n*2-1 + ista_fftp-1
                   inx(nn) = i-nlp ; inx(nn+1) = i-nlp
                   jnx(nn) = j2; jnx(nn+1) = j2
                   knx(nn) = k2; knx(nn+1) = k2
                end do
             end do
          end do
       end if
#else
       do n = ista_fftp, iend_fftp
          nn = (n+kimg-1)/kimg    ! nn = n (when kimg=1), nn = (n+1)/2 (when kimg=2)
#ifdef _MPIFFT_
!!$          if(nn == ista_fftp-1) write(nfout,'(" n == ista_fftp-1 = ",i8)') nn
          np0 = nn - idp*ly_d*lz*myrank_cdfft
          i  = mod(np0-1,idp)+1
          j0 = mod((np0-i)/idp,ly_d) + 1
          j  = j0 + ny_d*myrank_cdfft
          k  = (np0-i-(j0-1)*idp)/(idp*ly_d) + 1
          if(i >= idp) i = i - idp
#else
          i  = mod(nn,idp)
          j  = mod((nn-1)/idp,nd2p) + 1
          k  = (nn - (j-1)*idp - i)/(idp*nd2p) + 1
#endif
          inx(n) = i - 1;  if(2*inx(n) > nlp) inx(n) = inx(n) - nlp
          jnx(n) = j - 1;  if(2*jnx(n) > nmp) jnx(n) = jnx(n) - nmp
          knx(n) = k - 1;  if(2*knx(n) > nnp) knx(n) = knx(n) - nnp
       end do
#endif

! =============================== modified by K. Tagami ================ 11.0
!       allocate(chden_l(ista_fftp:iend_fftp,nspin)); chden_l = 0.d0 ! MPI
!       allocate(grad_trho(ista_fftph:iend_fftph))                   ! MPI
!       allocate(grad_rho(ista_fftph:iend_fftph,nspin)); grad_rho = 0.d0  ! MPI

       if ( noncol .and. sw_constrain_on_grad_correction == OFF ) then
          allocate(chden_l(ista_fftp:iend_fftp,ndim_magmom)); chden_l = 0.d0
          allocate(grad_trho(ista_fftph:iend_fftph))
          allocate(grad_rho(ista_fftph:iend_fftph,nspin)); grad_rho = 0.d0
       else
          allocate(chden_l(ista_fftp:iend_fftp,nspin)); chden_l = 0.d0 ! MPI
          allocate(grad_trho(ista_fftph:iend_fftph))                   ! MPI
          allocate(grad_rho(ista_fftph:iend_fftph,nspin)); grad_rho = 0.d0  ! MPI
       endif
! ======================================================================= 11.0

#ifndef _XC_SAVE_MEMORY_
       if(nrank_ggacmp > 1) then
          allocate(cgrad_rho(ista_fftph:iend_fftph,myrank_ggacmp+1:myrank_ggacmp+1,nspin))  ! MPI
       else
          allocate(cgrad_rho(ista_fftph:iend_fftph,3,nspin))  ! MPI
       end if
       cgrad_rho =0.d0
#else
       allocate(cggawk13(ista_fftph:iend_fftph))           ! MPI
#endif
       allocate(dF_drho(ista_fftph:iend_fftph,nspin))      ! MPI
       dF_drho = 0.d0
       allocate(dF_dgradrho(ista_fftph:iend_fftph,nspin))  ! MPI
       dF_dgradrho = 0.d0
    end if

! ====== KT_add ======== 13.0XX
    if ( use_metagga ) then
       allocate( ekin_dens(ista_fftph:iend_fftph,nspin) ); ekin_dens = 0.0d0
       allocate( grad2_rho(ista_fftph:iend_fftph,nspin) ); grad2_rho = 0.0d0
       if ( vtau_exists ) then
          allocate( dF_dtau(ista_fftph:iend_fftph,nspin) ); dF_dtau = 0.0d0
          allocate( dF_dlaplrho(ista_fftph:iend_fftph,nspin) ); dF_dlaplrho = 0.0d0
       endif
    endif
! ====================== 13.0XX

  contains
    subroutine set_f2or1(npes,ista,iend,f2or1)
      integer, intent(in) :: npes,ista, iend
      real(kind=DP), intent(out), dimension(ista:iend) :: f2or1
      integer :: idph,nlph,ip,i,j,k

      if(kimg == 1) then
         idph = idp/2
         nlph = nlp/2
#ifdef _MPIFFT_

         f2or1 = 0.d0
!!$         do j = 1, min(nz_d,nnp-nz_d*myrank_cdfft)*nmp
!!$         do j = 1, min(nz_d,nnp-nz_d*myrank_cdfft)*nd2p
         do k = 1, min(nz_d, nnp-nz_d*myrank_cdfft)
            do j = 1, nmp
               do i = 1, nlph
!!$                  ip = i + idph*(j-1) + idph*ly*lz_d*myrank_cdfft
                  ip = i + idph*(j-1) + idph*ly*(k-1) + idph*ly*lz_d*myrank_cdfft
                  f2or1(ip) = 2.d0
               end do
               ip = 1 + idph*(j-1) + idph*ly*(k-1) + idph*ly*lz_d*myrank_cdfft
               f2or1(ip) = 1.d0
               ip = nlph+1 + idph*(j-1) + idph*ly*(k-1) + idph*ly*lz_d*myrank_cdfft
               f2or1(ip) = 1.d0
            end do
         end do
!!$            ip = idph*(j-1) + 1 + idph*ly*lz_d*myrank_cdfft
!!$            f2or1(ip) = 1.d0
!!$            ip = idph*(j-1)+ nlph + 1 + idph*ly*lz_d*myrank_cdfft
!!$            f2or1(ip) = 1.d0
!!$         end do
#else
         f2or1 = 2.d0
         if(npes >= 2) then
            do i = 1, nd2p*nnp
               ip = idph*(i-1) + 1
               if(ip>= ista .and. ip <= iend) f2or1(ip) = 1.d0
            end do
            do i = 1, nd2p*nnp
               ip = idph*(i-1) + nlph + 1
               if(ip>= ista .and. ip <= iend) f2or1(ip) = 1.d0
            end do
            do j = nlph+2, idph
               do i = 1,nd2p*nnp
                  ip = idph*(i-1)+j
                  if(ip>= ista .and. ip <= iend) f2or1(ip) = 0.d0
               end do
            end do
            do j = nmp+1, nd2p
               do k = 1, nnp
                  do i = 1, nlph
                     ip = i + idph*(j-1) + idph*nd2p*(k-1)
                     if(ip >= ista .and. ip <= iend) f2or1(ip) = 0.d0
                  end do
               end do
            end do
! for SX6 ASL 20040817
            do k = nnp+1, nd3p
               do i = 1, idph*nd2p
                  ip = i + idph*nd2p*(k-1)
                  if(ip >= ista .and. ip <= iend) f2or1(ip) = 0.d0
               end do
            end do
         else
            do i = 1, nd2p*nnp
               ip = idph*(i-1) + 1
               f2or1(ip) = 1.d0
            end do
            do i = 1, nd2p*nnp
               ip = idph*(i-1) + nlph + 1
               f2or1(ip) = 1.d0
            end do
            do j = nlph+2, idph
               do i = 1,nd2p*nnp
                  ip = idph*(i-1)+j
                  f2or1(ip) = 0.d0
               end do
            end do
            do j = nmp+1, nd2p
               do k = 1, nnp
                  do i = 1, nlph
                     ip = i + idph*(j-1) + idph*nd2p*(k-1)
                     f2or1(ip) = 0.d0
                  end do
               end do
            end do
! for SX6 ASL 20040817
            do k = nnp+1, nd3p
               do i = 1, idph*nd2p
                  ip = i + idph*nd2p*(k-1)
                  f2or1(ip) = 0.d0
               end do
            end do
         end if
!!$       do i = ista, iend
!!$          if(mod(i*2,idp) == 2 .or. mod(i*2,idp) == 0) f2or1(i) = 1.d0
!!$       end do
#endif
      else
#ifdef _MPIFFT_
!!$         f2or1 = 0.d0                               ! f2or1 works to the fft data in R space.
!!$         do k = 1, min(nz_d,nnp-nz_d*myrank_cdfft)
!!$            do j = 1, nmp     ! nmp = fft_box_size_CD(2,1)
!!$               do i = 1, nlp  ! nlp = fft_box_size_CD(1,1)
!!$                  ip = i+(j-1)*idp+(k-1)*idp*nd2p+idp*nd2p*lz_d*myrank_cdfft
!!$                  f2or1(ip) = 1.d0
!!$               end do
!!$            end do
!!$         end do
         if(iprixc >= 2 ) write(nfout,'(" ix kimg = 2 <<set_f2or1>>")')
         f2or1 = 1.d0
         do j = nlp+1, idp      ! x
            do i = 1, ly*nz_d
               ip = idp*(i-1)+j+ista-1
               f2or1(ip) = 0.d0
            end do
         end do
         if(iprixc >= 2 ) write(nfout,'(" iy kimg = 2 <<set_f2or1>>")')
         do j = nmp+1,ly         ! y
            do k = 1, nz_d
               do i = 1, nlp
                  ip = i + idp*(j-1) + idp*ly*(k-1) + ista-1
                  f2or1(ip) = 0.d0
               end do
            end do
         end do
         if(iprixc >= 2 ) write(nfout,'(" iz kimg = 2 <<set_f2or1>>")')
         do  k = nz_d+1, lz_d   ! z
            do i = 1, idp*ly
               ip = i + idp*ly*(k-1) + ista-1
               f2or1(ip) = 0.d0
            end do
         end do
#else
         f2or1 = 1.d0
         if(npes >= 2) then
            do j = nlp+1, idp    ! x
               do i = 1, nd2p*nnp
                  ip = idp*(i-1)+j
                  if(ip>= ista .and. ip <= iend) f2or1(ip) = 0.d0
               end do
            end do
            do j = nmp+1, nd2p   ! y
               do k = 1, nnp
                  do i = 1, nlp
                     ip = i + idp*(j-1) + idp*nd2p*(k-1)
                     if(ip >= ista .and. ip <= iend) f2or1(ip) = 0.d0
                  end do
               end do
            end do
! for SX6 ASL 20040817
            do k = nnp+1, nd3p   ! z
               do i = 1, idp*nd2p
                  ip = i + idp*nd2p*(k-1)
                  if(ip >= ista .and. ip <= iend) f2or1(ip) = 0.d0
               end do
            end do
         else
            do j = nlp+1, idp    ! x
               do i = 1, nd2p*nnp
                  ip = idp*(i-1)+j
! ================================ modifed by K. Tagami ====( uncertain )== 11.0
!                  f2or1(ip) = 0.d0
                  if(ip >= ista .and. ip <= iend) f2or1(ip) = 0.d0
! ===================================================================== 11.0
               end do
            end do
            do j = nmp+1, nd2p   ! y
               do k = 1, nnp
                  do i = 1, nlp
                     ip = i + idp*(j-1) + idp*nd2p*(k-1)
! ================================ modifed by K. Tagami ====( uncertain )== 11.0
!                     f2or1(ip) = 0.d0
                     if(ip >= ista .and. ip <= iend) f2or1(ip) = 0.d0
! ===================================================================== 11.0
                  end do
               end do
            end do
! for SX6 ASL 20040817
            do k = nnp+1, nd3p   ! z
               do i = 1, idp*nd2p
                  ip = i + idp*nd2p*(k-1)
! ================================ modifed by K. Tagami ====( uncertain )== 11.0
!                  f2or1(ip) = 0.d0
                  if(ip >= ista .and. ip <= iend) f2or1(ip) = 0.d0
! ===================================================================== 11.0
               end do
            end do
         end if
#endif
      end if
    end subroutine set_f2or1
  end subroutine xc_allocate

  subroutine xc_deallocate
#ifndef _MPIFFT_
    call m_FFT_dealloc_CD_box()
#endif
    deallocate(afft)
    deallocate(chgrhr_l)
    deallocate(f2or1)

    if(check_of_xctype() == GGA) then
       deallocate(inx); deallocate(jnx); deallocate(knx)
       deallocate(chden_l)
       deallocate(grad_trho)
#ifndef _XC_SAVE_MEMORY_
       deallocate(cgrad_rho)
#else
       deallocate(cggawk13)
#endif
       deallocate(grad_rho); deallocate(dF_drho)
       deallocate(dF_dgradrho)
    end if

! ====== KT_add ======== 13.0XX
    if ( use_metagga ) then
       deallocate( ekin_dens );   deallocate( grad2_rho )
       if ( allocated( dF_dtau ) ) deallocate( dF_dtau )
       if ( allocated( dF_dlaplrho ) ) deallocate( dF_dlaplrho )
    endif
! ====================== 13.0XX
  end subroutine xc_deallocate

  subroutine m_XC_alloc_s_gga
    allocate(s_gga1(3,3)); s_gga1=0.d0
    allocate(s_gga2(3,3)); s_gga2=0.d0
  end subroutine m_XC_alloc_s_gga

  subroutine m_XC_dealloc_s_gga
    deallocate(s_gga1)
    deallocate(s_gga2)
  end subroutine m_XC_dealloc_s_gga

  subroutine m_XC_cal_potential(nfout,input_charge,chgq_l, vflag  &
    &  , chgsoft,hsr,hsrd)
!   Upgraded on 19th Sep. 2006 by T. Yamasaki
!     * MPIFFT
    integer, intent(in) :: nfout, input_charge

! =================================== modified by K. Tagami =============== 11.0
!    real(DP),intent(in) :: chgq_l(ista_kngp:iend_kngp,kimg,nspin)
    real(DP),intent(in) :: chgq_l(ista_kngp:iend_kngp,kimg,ndim_magmom )
! ========================================================================= 11.0
    !                                    charge density
    integer, intent(in) :: vflag

! ================================== modified by K. Tagami ================ 11.0
!    real(DP),intent(in),optional,dimension(ista_kngp:iend_kngp,kimg,nspin) &
!         &              :: chgsoft  ! soft charge density
!    real(DP),intent(in),optional,dimension(natm,nlmt,nlmt,nspin)     :: hsr
!    real(DP),intent(in),optional,dimension(natm,nlmt,nlmt,nspin,3,3) :: hsrd
!
    real(DP),intent(in),optional,dimension(ista_kngp:iend_kngp,kimg,ndim_magmom) &
         &              :: chgsoft  ! soft charge density
    real(DP),intent(in),optional,dimension(natm,nlmt,nlmt,ndim_magmom)     :: hsr
    real(DP),intent(in),optional,dimension(natm,nlmt,nlmt,ndim_magmom,3,3) :: hsrd
!
! ========================================================================== 11.0

! ================================== added by K. Tagami ================== 11.0
    Real(kind=DP), allocatable :: RhoMag_G(:,:,:)
    Real(kind=DP), allocatable :: RhoMag_R(:,:)
    Real(kind=DP), allocatable :: Rot_Angles(:,:)
    Real(kind=DP), allocatable :: bfft_kt(:,:)
    Real(kind=DP), allocatable :: bfft(:)
    integer, allocatable :: quantz_axis_inversion_flg_mesh(:)
! ======================================================================== 11.0

! ================================== added by K. Tagami ================== 11.0
!!!!!    integer, allocatable :: quantz_axis_inversion_flg_atm(:)
    real(kind=DP), allocatable :: chgsoft_projected_on_Gspace(:,:,:)
    real(kind=DP), allocatable :: chgsoft_on_Rspace(:,:)
    real(kind=DP), allocatable :: chgsoft_along_QuantzAxis(:,:)
    real(kind=DP), allocatable :: hsr_along_QuantzAxis(:,:,:,:)
    real(kind=DP), allocatable :: hsrd_along_QuantzAxis(:,:,:,:,:,:)
! ======================================================================= 11.0

! ==== KT_add === 2014/08/04
    Real(kind=DP), allocatable :: MagVorticity(:,:)
! =============== 2014/08/04

    integer       :: ispin, iloop, i
    logical       :: nopcc
    real(kind=DP) :: rinplw
    integer       :: id_sname = -1
    call tstatc0_begin('m_XC_cal_potential ',id_sname,level=1)

    rinplw = 1.d0/product(fft_box_size_CD(1:3,1))

    if(input_charge == Partial_Core_Charge) then
       if(iprixc >= 1) write(nfout, '(" input_charge == Partial_Core_Charge")')
       call check_of_pcc(nopcc)
       if(nopcc) goto 4001
    endif

! ============================= added by K. Tagami ================== 11.0
    if ( noncol ) then
      if ( input_charge == Valence_plus_PC_Charge ) then
!        allocate( Rot_Angles(ista_fftph:iend_fftph,2 )) ; Rot_Angles = 0.0d0
        call check_of_pcc(nopcc)
      endif
    endif
! ===================================================================== 11.0
    call set_ispin(ispin)
          if(iprixc >= 2) write(nfout,*) ' ! ispin = ', ispin
    call xc_allocate(ispin,nfout)

! ========================== modified by K. Tagami ==================== 11.0
!    do iloop = 1, ispin
!       if(iprixc >= 2) write(nfout,'(" ! (m_XC_cal_potential) , iloop = ",i5)') iloop
!#ifdef _OLD_MAP_CHARGE_
!       call map_charge_onto_a_fft_box(chgq_l)     !-(m_XC_Pot.) -> afft(*) (xcchg2)
!#else
!       call map_charge_onto_a_fft_box2(nfout,chgq_l)     !-(m_XC_Pot.) -> afft(*) (xcchg2)
!#endif
!             if(iprixc >= 2) write(nfout,*) " just after map_charge_onto_a_fft_box "
!             if(iprixc >= 2) call wd_small_part_of_afft(7,'G space',120)
!             if(iprixc >= 2) then
!                write(nfout,'(" ista_fftp, iend_fftp, np_fftp, mp_fftp = ",4i8)') &
!                     & ista_fftp, iend_fftp, np_fftp, mp_fftp
!             end if
!       if(check_of_xctype() == GGA) then
!          chden_l(ista_fftp:iend_fftp,iloop) = afft(ista_fftp:iend_fftp)
!       end if
!       if(myrank_ggacmp < nrank_ggacmp) then
!          call m_FFT_CD_inverse_c(nfout,afft)  ! afft(G_sp.) -> afft(R_sp.)
!             if(iprixc >= 2) call wd_small_part_of_afft(7,'R space',120)
!          call m_FFT_check_of_negative_CD(npes_cdfft,ista_fftp,iend_fftp&
!               & ,ista_fftph,iend_fftph,afft,nfout,nspin,iloop)
!       else
!          call m_FFT_CD_inverse_c(nfout,afft,mode=0)
!       end if
!!!$       if(npes > 1 .and. max_ggacmp >= nrank_ggacmp) call mpi_barrier(MPI_CommGroup,ierr)
!       if(iprixc >= 2) write(nfout,'(" ! out of m_FFT_check_of_negative_CD <<m_XC_cal_potential>>")')
!
!       call cp_afft_to_chgrhr        ! -(contained here) afft -> chgrhr_l
!    end do

    if ( noncol ) then
       if ( input_charge == Partial_Core_Charge .and. vflag == EXC_ONLY ) then
          call set_chgrhr_case_collinear( ispin )
       else
          if ( sw_constrain_on_grad_correction == ON ) then
             call set_chgrhr_case_noncollinear2
          else
             call set_chgrhr_case_noncollinear3
          endif
       endif
    else
       call set_chgrhr_case_collinear( ispin )
    endif

! --
!    write(nfout,*) 'chden_l'
!    Do i=ista_fftp, iend_fftp
!       write(nfout,*) i, chden_l(i,1), chden_l(i,2)
!   End do
!   stop
! ======================================================================== 11.0
! ======================================================================== 11.0

! ==== KT_add ==== 13.0XX, 2014/09/19
    if ( use_metagga .and. .not. use_modeled_ekin_density ) then
       if ( noncol ) then
       else
          if ( use_asymm_ekin_density ) then
             call set_ekindens_case_collinear( ispin, ekina_l, ekin_dens )
                         ! (asymmetric) ekin_l (G) -> ekin_dens (R)
          else if ( use_symm_ekin_density ) then
             call set_ekindens_case_collinear( ispin, ekins_l, ekin_dens )
                         ! (symmetric) ekin_l (G) -> ekin_dens (R)
             if ( sw_rspace_ekin_density == ON &
                  &        .and. sw_calc_ekin_density_hardpart == ON &
                  &        .and. sw_add_ekin_hardpart_on_Gspace == OFF ) then
                call m_KE_add_ekindens_hard_rspace( ista_fftph, iend_fftph, ekin_dens )
             endif
          endif
       endif
    endif
! ================= 13.0XX, 2014/09/19

    call check_lmn_even  ! -(contained in subr. m_XC_cal_potential) ->lmn_even

    if(check_of_xctype() == GGA) then
       if(vflag == STRESS_) then
          call ggaxcp_diff           ! -(contained here) ->chgrhr,exc,etc.
          !    ~~~~~~~~~~~
       else
          call ggaxcp()              ! -(contained here) ->chgrhr,exc
          !    ~~~~~~
       endif
    else
       call xcpotf(ispin,input_charge) ! -(m_XC_Potential) ->chgrhr,exc
       !    ~~~~~~
    end if
    if(vflag == VXC_AND_EXC .or. vflag == STRESS_) then
! ======================= modified by K. Tagami ====================== 11.0
!       do iloop = 1, ispin
!          if(myrank_ggacmp < nrank_ggacmp) then
!             call cpafft(iloop)                     ! chgrhr_l -> afft
!             call m_FFT_CD_direct_c(nfout,afft)     ! afft(R space) -> afft(G sp.)
!          else
!             call m_FFT_CD_direct_c(nfout,afft,mode=0)
!          end if
!             if(iprixc >= 2) call wd_small_part_of_afft(7,'G space',120)
!             if(npes > 1) call mpi_barrier(MPI_CommGroup,ierr)
!          call cpafft_to_vxc_or_vxcpc(iloop)      ! -> vxc_l, vxcpc_l
!       end do
!       if(iprixc >= 2) call wd_small_part_of_vxc
!
        if ( noncol ) then
!!!!!          call set_vxc_case_noncollinear
          call set_vxc_case_noncollinear2
        else
          call set_vxc_case_collinear(ispin)
        endif
! ===================================================================== 11.0
    end if

    call xc_deallocate

! ============================= added by K. Tagami ================== 11.0
    if ( noncol ) then
!       deallocate( Rot_Angles )
    endif
! ===================================================================== 11.0

! ==== KT_add === 2014/08/04
    if ( noncol .and. sw_monitor_magnetic_vorticity == ON ) then
       if ( allocated( MagVorticity) ) deallocate( MagVorticity)
    endif
! =============== 2014/08/04

4001 continue

    if ( input_charge == Partial_Core_Charge ) excpc = exc
    call tstatc0_end(id_sname)
  contains

! ============================== added by K. Tagami ================ 11.0
    subroutine set_chgrhr_case_collinear( ispin )
      integer, intent(in) :: ispin
      integer :: iloop

      do iloop = 1, ispin
         if(iprixc >= 2) write(nfout,'(" ! (m_XC_cal_potential) , iloop = ",i5)') iloop
#ifdef _OLD_MAP_CHARGE_
         call map_charge_onto_a_fft_box( chgq_l, nspin, iloop, input_charge )
	                                     !-(m_XC_Pot.) -> afft(*) (xcchg2)
#else
         call map_charge_onto_a_fft_box2( nfout, chgq_l, nspin, iloop, &
	&                                 input_charge )
                                              !-(m_XC_Pot.) -> afft(*) (xcchg2)
#endif
         if(iprixc >= 2) write(nfout,*) " just after map_charge_onto_a_fft_box "
         if(iprixc >= 2) call wd_small_part_of_afft(7,'G space',120)
         if(iprixc >= 2) then
            write(nfout,'(" ista_fftp, iend_fftp, np_fftp, mp_fftp = ",4i8)') &
                     & ista_fftp, iend_fftp, np_fftp, mp_fftp
         end if
         if(check_of_xctype() == GGA) then
            chden_l(ista_fftp:iend_fftp,iloop) = afft(ista_fftp:iend_fftp)

         end if

         if(myrank_ggacmp < nrank_ggacmp) then
            call m_FFT_CD_inverse_c(nfout,afft)  ! afft(G_sp.) -> afft(R_sp.)
            if(iprixc >= 2) call wd_small_part_of_afft(7,'R space',120)
            call m_FFT_check_of_negative_CD(npes_cdfft,ista_fftp,iend_fftp&
               & ,ista_fftph,iend_fftph,afft,nfout,nspin,iloop)
         else
            call m_FFT_CD_inverse_c(nfout,afft,mode=0)
         end if
!!$       if(npes > 1 .and. max_ggacmp >= nrank_ggacmp) call mpi_barrier(MPI_CommGroup,ierr)
         if(iprixc >= 2) write(nfout,'(" ! out of m_FFT_check_of_negative_CD <<m_XC_cal_potential>>")')

         call cp_afft_to_chgrhr(iloop)        ! -(contained here) afft -> chgrhr_l
      end do
    end subroutine set_chgrhr_case_collinear

    subroutine set_chgrhr_case_noncollinear3
      integer :: iloop
! --
      allocate( RhoMag_R(ista_fftph:iend_fftph,ndim_magmom)) ; RhoMag_R = 0.0d0
! --------------------------------------------
! 1) calc pcc
! -------------------------------------------
      if ( .not. nopcc ) then
         allocate( bfft(ista_fftp:iend_fftp) ); bfft = 0.0d0
#ifdef _OLD_MAP_CHARGE_
         call map_charge_onto_a_fft_box( chgq_l, ndim_magmom, 1, Partial_Core_Charge )
#else
         call map_charge_onto_a_fft_box2( nfout, chgq_l, ndim_magmom, 1, &
              &                           Partial_Core_Charge )
#endif
         bfft = afft
      endif
! --------------------------------------------
! 2) add valence contribution
! -------------------------------------------
      Do iloop = 1, ndim_magmom
         afft = 0.0d0
#ifdef _OLD_MAP_CHARGE_
         call map_charge_onto_a_fft_box( chgq_l, ndim_magmom, iloop, &
        &                                Valence_Charge_Only )
                                            !-(m_XC_Pot.) -> afft(*) (xcchg2)
#else
         call map_charge_onto_a_fft_box2( nfout, chgq_l, ndim_magmom, iloop, &
        &                                 Valence_Charge_Only )
                                              !-(m_XC_Pot.) -> afft(*) (xcchg2)
#endif
         if(iprixc >= 2) write(nfout,*) " just after map_charge_onto_a_fft_box p1"
         if(iprixc >= 2) call wd_small_part_of_afft(7,'G space',120)
!

         if ( (.not. nopcc) .and. iloop == 1 ) afft = afft + bfft

         if (check_of_xctype() == GGA) then
            chden_l(ista_fftp:iend_fftp,iloop) = afft(ista_fftp:iend_fftp)
         end if

         call m_FFT_CD_inverse_c(nfout,afft)  ! afft(G_sp.) -> afft(R_sp.)
         if (iprixc >= 2) then
            write(nfout,'(" ! out of m_FFT_check_of_negative_CD <<m_XC_cal_potential>>")')
         endif
         call cp_afft_to_MagMom_Rspace(iloop)

      End do

      if (allocated( bfft )) deallocate( bfft )

! -----------------
! 3) project on local quan. axis
! -----------------
      allocate( quantz_axis_inversion_flg_mesh(ista_fftph:iend_fftph ) )
      quantz_axis_inversion_flg_mesh = 0

      call m_ES_SpinDens_Along_QuantzAxis2( RhoMag_R, chgrhr_l, &
           &                                quantz_axis_inversion_flg_mesh, .false. )

    end subroutine set_chgrhr_case_noncollinear3

    subroutine set_chgrhr_case_noncollinear2
      integer :: iloop
! -----------------------------------------------------------------
!!
!!! 1) calculate the up/down spin densities along the local quantization axis
!!                             ( valence electron )
! ----------------------------------------------------------------
      integer :: i

      allocate( RhoMag_R(ista_fftph:iend_fftph,ndim_magmom)) ; RhoMag_R = 0.0d0
      call get_MagMom_on_Rspace                ! chgq_l --> RhoMag_R
!
      allocate( quantz_axis_inversion_flg_mesh(ista_fftph:iend_fftph ) )
      quantz_axis_inversion_flg_mesh = 0

      call m_ES_SpinDens_Along_QuantzAxis2( RhoMag_R, chgrhr_l, &
           &                                quantz_axis_inversion_flg_mesh, .true. )

! ----------------------------------------------------------------
!!
!!! 2) add pcc contribution to chgrhr_l
!!
! ---------------------------------------------------------------
      if ( .not. nopcc ) then
         allocate( bfft(ista_fftp:iend_fftp) ); bfft = 0.0d0
      endif

      Do iloop=1, nspin
! --------------------- set valence charge on G-space ----------
        afft = 0.0d0

!!! ------------------------------- origial ----
!!!        if ( myrank_ggacmp < nrank_ggacmp ) then
!!!	  call cp_MagMom_to_afft_Rspace( iloop )
!!!          call m_FFT_CD_direct_c( nfout, afft )
!!!        else
!!!          call m_FFT_CD_direct_c( nfout, afft, mode=0 )
!!!        end if

! ---- for debug ----
        call cp_MagMom_to_afft_Rspace( iloop )
        call m_FFT_CD_direct_c( nfout, afft )
! -
!
! --
        afft = afft * rinplw

! -- dbug
!        goto 1200
! --------------------- set core charge on G-space ----------
        if ( .not. nopcc ) then

           bfft = afft                ! save afft
           afft = 0.0d0
!
! ----
!                      In the followings, chgq_l is a dummy matrix.
! ----
!
#ifdef _OLD_MAP_CHARGE_
           call map_charge_onto_a_fft_box( chgq_l, ndim_magmom, iloop, &
	&                                  Partial_Core_Charge )
#else
           call map_charge_onto_a_fft_box2( nfout, chgq_l, ndim_magmom, iloop, &
	&                                   Partial_Core_Charge )
#endif

           afft = bfft + afft / dble(nspin)         ! nspin == 2
        end if
! --
1200    continue

        if ( check_of_xctype() == GGA ) then
           chden_l(ista_fftp:iend_fftp,iloop) = afft(ista_fftp:iend_fftp)
        end if

! -----------------------------------------------
!  This if-sentence is not necessary.
!
!!!!!!!!        if ( .not. nopcc ) then
! -----------------------------------------------

           if ( myrank_ggacmp < nrank_ggacmp ) then
              call m_FFT_CD_inverse_c( nfout, afft )    ! afft(G_sp.) -> afft(R_sp.)
              if ( iprixc >= 2 ) call wd_small_part_of_afft( 7,'R space',120 )
              call m_FFT_check_of_negative_CD( npes_cdfft, ista_fftp, iend_fftp, &
                   &                           ista_fftph, iend_fftph, afft, &
                   &                           nfout, nspin, iloop )
           else
              call m_FFT_CD_inverse_c( nfout, afft, mode=0 )
           end if
           if (iprixc >= 2) then
              write(nfout,'(" ! out of m_FFT_check_of_negative_CD <<m_XC_cal_potential>>")')
           endif

           call cp_afft_to_chgrhr(iloop)       ! -(contained here) afft -> chgrhr_l

! ------------------------------
!!!!!!!!!!!!!!        endif
! -----------------------------
      end do

      if ( allocated(bfft) ) deallocate( bfft )

    end subroutine set_chgrhr_case_noncollinear2

    subroutine set_chgrhr_case_noncollinear
      integer :: iloop
! -----------------------------------------------------------------
!!
!!! 1) calculate the up/down spin densities along the local quantization axis
!!                             ( valence electron )
! ----------------------------------------------------------------
      integer :: i

      allocate( RhoMag_R(ista_fftph:iend_fftph,ndim_magmom)) ; RhoMag_R = 0.0d0
      call get_MagMom_on_Rspace                ! chgq_l --> RhoMag_R
!
      call m_ES_get_Angles_MagMom_Rspace2( RhoMag_R, Rot_Angles )
      call m_ES_SpinDens_Along_QuantzAxis( RhoMag_R, Rot_Angles, chgrhr_l )

      deallocate( RhoMag_R )

! ----------------------------------------------------------------
!!
!!! 2) add pcc contribution to chgrhr_l
!!
! ---------------------------------------------------------------
      if ( .not. nopcc ) then
         allocate( bfft(ista_fftp:iend_fftp) ); bfft = 0.0d0
      endif

      Do iloop=1, nspin
! --------------------- set valence charge on G-space ----------
        afft = 0.0d0

!!! ------------------------------- origial ----
!!!        if ( myrank_ggacmp < nrank_ggacmp ) then
!!!	  call cp_MagMom_to_afft_Rspace( iloop )
!!!          call m_FFT_CD_direct_c( nfout, afft )
!!!        else
!!!          call m_FFT_CD_direct_c( nfout, afft, mode=0 )
!!!        end if

! ---- for debug ----
        call cp_MagMom_to_afft_Rspace( iloop )
        call m_FFT_CD_direct_c( nfout, afft )
! -
!
! --
        afft = afft * rinplw

! -- dbug
!        goto 1200
! --------------------- set core charge on G-space ----------
        if ( .not. nopcc ) then

           bfft = afft                ! save afft
           afft = 0.0d0
!
! ----
!                      In the followings, chgq_l is a dummy matrix.
! ----
!
#ifdef _OLD_MAP_CHARGE_
           call map_charge_onto_a_fft_box( chgq_l, ndim_magmom, iloop, &
	&                                  Partial_Core_Charge )
#else
           call map_charge_onto_a_fft_box2( nfout, chgq_l, ndim_magmom, iloop, &
	&                                   Partial_Core_Charge )
#endif

           afft = bfft + afft / dble(nspin)         ! nspin == 2
        end if
! --
1200    continue

        if ( check_of_xctype() == GGA ) then
           chden_l(ista_fftp:iend_fftp,iloop) = afft(ista_fftp:iend_fftp)
        end if

! -----------------------------------------------
!  This if-sentence is not necessary.
!
!!!!!!!!        if ( .not. nopcc ) then
! -----------------------------------------------

           if ( myrank_ggacmp < nrank_ggacmp ) then
              call m_FFT_CD_inverse_c( nfout, afft )    ! afft(G_sp.) -> afft(R_sp.)
              if ( iprixc >= 2 ) call wd_small_part_of_afft( 7,'R space',120 )
              call m_FFT_check_of_negative_CD( npes_cdfft, ista_fftp, iend_fftp, &
                   &                           ista_fftph, iend_fftph, afft, &
                   &                           nfout, nspin, iloop )
           else
              call m_FFT_CD_inverse_c( nfout, afft, mode=0 )
           end if
           if (iprixc >= 2) then
              write(nfout,'(" ! out of m_FFT_check_of_negative_CD <<m_XC_cal_potential>>")')
           endif

           call cp_afft_to_chgrhr(iloop)       ! -(contained here) afft -> chgrhr_l

! ------------------------------
!!!!!!!!!!!!!!        endif
! -----------------------------
      end do

      if ( allocated(bfft) ) deallocate( bfft )

    end subroutine set_chgrhr_case_noncollinear

! ------------------------------------------------------------
!!
!!!  n(G), mx(G), my(G), mz(G) ---> n(R), mx(R), my(R), mz(R)
!!
! -------------------------------------------------------------
    subroutine get_MagMom_on_Rspace
      integer :: iloop

      do iloop = 1, ndim_magmom
         afft = 0.0d0
#ifdef _OLD_MAP_CHARGE_
         call map_charge_onto_a_fft_box( chgq_l, ndim_magmom, iloop, &
	&                                Valence_Charge_Only )
	                                     !-(m_XC_Pot.) -> afft(*) (xcchg2)
#else
         call map_charge_onto_a_fft_box2( nfout, chgq_l, ndim_magmom, iloop, &
	&                                 Valence_Charge_Only )
                                              !-(m_XC_Pot.) -> afft(*) (xcchg2)
#endif
         if(iprixc >= 2) write(nfout,*) " just after map_charge_onto_a_fft_box p1"
         if(iprixc >= 2) call wd_small_part_of_afft(7,'G space',120)

!!! ------------ original ---
!         if(myrank_ggacmp < nrank_ggacmp) then
!            call m_FFT_CD_inverse_c(nfout,afft)  ! afft(G_sp.) -> afft(R_sp.)
!            if(iprixc >= 2) call wd_small_part_of_afft(7,'R space',120)
!            call m_FFT_check_of_negative_CD( npes_cdfft, ista_fftp, iend_fftp, &
!                 &                           ista_fftph, iend_fftph, afft, &
!                 &                           nfout, 100, iloop )
!         else
!            call m_FFT_CD_inverse_c(nfout,afft,mode=0)
!         end if

!----
!         if ( .not. nopcc ) then
!            call m_FFT_CD_inverse_c(nfout,afft)  ! afft(G_sp.) -> afft(R_sp.)
!         else
!            if(myrank_ggacmp < nrank_ggacmp) then
!               call m_FFT_CD_inverse_c(nfout,afft)  ! afft(G_sp.) -> afft(R_sp.)
!!               if(iprixc >= 2) call wd_small_part_of_afft(7,'R space',120)
!               call m_FFT_check_of_negative_CD( npes_cdfft, ista_fftp, iend_fftp, &
!                    &                           ista_fftph, iend_fftph, afft, &
!                    &                           nfout, 100, iloop )
!            else
!               call m_FFT_CD_inverse_c(nfout,afft,mode=0)
!            end if
!         endif

! ------ for debug --
         call m_FFT_CD_inverse_c(nfout,afft)  ! afft(G_sp.) -> afft(R_sp.)
! --
         call cp_afft_to_MagMom_Rspace(iloop)
      end do

    end subroutine get_MagMom_on_Rspace

    subroutine cp_MagMom_to_afft_Rspace( iloop )
      integer, intent(in) :: iloop
      integer :: i

      do i = ista_fftph, iend_fftph    ! MPI
         afft(2*i-1) = chgrhr_l(i,iloop)
      end do
    end subroutine cp_MagMom_to_afft_Rspace

    subroutine cp_afft_to_MagMom_Rspace(iloop)
      integer, intent(in) :: iloop
      integer :: i

! ==================
!      if(myrank_ggacmp < nrank_ggacmp) then
!         do i = ista_fftph, iend_fftph    ! MPI
!            RhoMag_R(i,iloop) = afft(i*2-1)
!         end do
!      else
!         do i = ista_fftph, iend_fftph    ! MPI
!            RhoMag_R(i,iloop) = 0.d0
!         end do
!      end if

! --- for debug -------------- This does not resolve the problem.
!
      do i = ista_fftph, iend_fftph    ! MPI
         RhoMag_R(i,iloop) = afft(i*2-1)
      end do
! ----

! ----------------
!      Do i = ista_fftph, iend_fftph
!         if ( RhoMag_R( i,iloop ) < 0.0 ) then
!            RhoMag_R( i,iloop ) = 0.0d0
!         end if
!      End Do

      if(iprixc >= 2) then
         write(nfout,*) "!XC == "
         write(nfout,*) "!XC        iloop = ", iloop

         write(nfout,'(" !XC -- RhoMag_R -- <<cp_afft_to_MagMom_Rspace>>")')
         write(nfout,'(" !XC ",6f12.6)') (RhoMag_R(i,iloop),i=ista_fftph,ista_fftph+17)
         write(nfout,'(" !XC -- afft --     <<cp_afft_to_MagMom_Rspace>>")')
         write(nfout,'(" !XC ",6f12.6)') (afft(i*2-1),i=ista_fftph,ista_fftph+17)
      end if
    end subroutine cp_afft_to_MagMom_Rspace

    subroutine set_vxc_case_collinear(ispin)
      integer, intent(in) :: ispin
      integer :: iloop
      real(kind=DP), allocatable :: bfft(:)

      if ( use_metagga .and. vtau_exists ) allocate( bfft(ista_fftp:iend_fftp) )
      do iloop = 1, ispin
         if(myrank_ggacmp < nrank_ggacmp) then
            call cpafft(iloop)                     ! chgrhr_l -> afft
            call m_FFT_CD_direct_c(nfout,afft)     ! afft(R space) -> afft(G sp.)
         else
            call m_FFT_CD_direct_c(nfout,afft,mode=0)
         end if
         if(iprixc >= 2) call wd_small_part_of_afft(7,'G space',120)
         if ( use_metagga .and. vtau_exists ) then
            if (myrank_ggacmp < nrank_ggacmp) then
               call cp_vals_to_afft_rspace1( iloop, dF_dtau, bfft )
               call m_FFT_CD_direct_c(nfout,bfft)     ! afft(R space) -> afft(G sp.)
            else
               call m_FFT_CD_direct_c(nfout,bfft,mode=0)
            end if
         endif
         if(npes > 1) call mpi_barrier(MPI_CommGroup,ierr)
         call cpafft_to_vxc_or_vxcpc(iloop)      ! -> vxc_l, vxcpc_l
         if ( vtau_exists ) call cpafft_to_vtau(iloop,bfft)         ! -> vtau_l
      end do
      if(iprixc >= 2) call wd_small_part_of_vxc
    end subroutine set_vxc_case_collinear

    subroutine set_vxc_case_noncollinear
      integer :: iloop
      real(kind=DP), allocatable :: vxc_magmom(:,:,:)

      do iloop = 1, ndim_chgpot
         if ( myrank_ggacmp < nrank_ggacmp ) then
            call m_ES_cp_VxcR_to_afft( iloop, chgrhr_l, Rot_Angles, afft )
            call m_FFT_CD_direct_c(nfout,afft)     ! afft(R space) -> afft(G sp.)
         else
            call m_FFT_CD_direct_c(nfout,afft,mode=0)
         end if
         if(iprixc >= 2) call wd_small_part_of_afft(7,'G space',120)
         if(npes > 1) call mpi_barrier(MPI_CommGroup,ierr)

         call cpafft_to_vxc_or_vxcpc(iloop)      ! -> vxc_l, vxcpc_l
      end do

! ---
      if(iprixc >= 2) call wd_small_part_of_vxc_ssrep
                                                 !ss-representation

! -- convert to magmom-representation

      allocate( vxc_magmom( ista_kngp:iend_kngp,kimg,ndim_magmom ) )
      vxc_magmom = 0.0d0
      call m_ES_DensMat_to_MagMom_vlhxcl( vxc_l, vxc_magmom )
      vxc_l = vxc_magmom
      deallocate( vxc_magmom )

    end subroutine set_vxc_case_noncollinear

    subroutine set_vxc_case_noncollinear2
      integer :: iloop

      do iloop = 1, ndim_magmom
         if ( myrank_ggacmp < nrank_ggacmp ) then
            call m_ES_cp_VxcR_to_afft2( iloop, quantz_axis_inversion_flg_mesh, &
                 &                      RhoMag_R, chgrhr_l, afft )
            call m_FFT_CD_direct_c(nfout,afft)     ! afft(R space) -> afft(G sp.)
         else
            call m_FFT_CD_direct_c(nfout,afft,mode=0)
         end if
         if(iprixc >= 2) call wd_small_part_of_afft(7,'G space',120)
         if(npes > 1) call mpi_barrier(MPI_CommGroup,ierr)

         call cpafft_to_vxc_or_vxcpc(iloop)      ! -> vxc_l, vxcpc_l
      end do

! ---
!!      if(iprixc >= 2) call wd_small_part_of_vxc_ssrep
!!                                                 !ss-representation

      deallocate( RhoMag_R )
      deallocate( quantz_axis_inversion_flg_mesh )

    end subroutine set_vxc_case_noncollinear2
! ====================================================================== 11.0

! ==== KT_add ==== 13.0XX
    subroutine set_ekindens_case_collinear( ispin, ekin_l, ekin_dens )
      integer, intent(in) :: ispin
      real(kind=DP), intent(in)  :: ekin_l(ista_kngp:iend_kngp,kimg,nspin)
      real(kind=DP), intent(out) :: ekin_dens(ista_fftph:iend_fftph,nspin)

      integer :: iloop, i

      do iloop = 1, ispin
         call map_charge_onto_a_fft_box2( nfout, ekin_l, nspin, iloop, &
        &                                 Valence_Charge_Only )
         if(myrank_ggacmp < nrank_ggacmp) then
            call m_FFT_CD_inverse_c(nfout,afft)  ! afft(G_sp.) -> afft(R_sp.)
            do i = ista_fftph, iend_fftph    ! MPI
               ekin_dens(i,iloop) = afft(i*2-1)
            end do
         else
            call m_FFT_CD_inverse_c(nfout,afft,mode=0)
            do i = ista_fftph, iend_fftph    ! MPI
               ekin_dens(i,iloop) = 0.d0
            end do
         end if
      end do
    end subroutine set_ekindens_case_collinear
! ================= 13.0XX

    subroutine ggaxcp()
      !        ~~~~~~
      integer       :: is, id_sname = -1
      real(kind=DP) :: exc_mpi
      integer :: i,ierr
      call tstatc0_begin('ggaxcp(in xc_pot.) ',id_sname)

      if(iprixc >= 2) write(nfout,*) " -- start ggaxcp -- "
!!$      if(xctype /= 'ldapw91') call charge_in_gspace() ! chgrhr --> chden_l
#ifdef LIBXC
      if ( xctype == "libxc" ) then
         call ggaxcp0_libxc()
      else
         call ggaxcp0()
      endif
#else
      call ggaxcp0()
#endif
      !    ~~~~~~~
             !   dF/d|rho(r)| (=vxc) --> dF_drho
             !   dFx/d|grad(rho(r))| --> dF_dgradrho
             !   dFc/d|grad(rho(r))| --> grad_trho
             !   ->exc

      if (iprixc >= 2) write(nfout,'(" !XC after ggaxcp0 <<ggaxcp>>")')

      if (vflag == VXC_AND_EXC .or. vflag == STRESS_) then
         if (xctype /= 'ldapw91' .and. xctype /= 'ldapbe ') then

! ======================================== modified by K. Tagami ============= 11.0
!            call dFxc_over_ddgradrho()
                    !   sum(q=xyz)(dFxc/dq(d|grad(rho)|)) --> grad_rho
!
            if ( noncol .and. sw_constrain_on_grad_correction == OFF ) then
               call dFxc_over_ddgradrho_noncl()
            else
               call dFxc_over_ddgradrho()
            endif
! ============================================================================== 11.0
         end if

! ==== KT_add ====== 13.0XX
         if ( use_metagga ) then
            if ( xctype == "tb09" ) then
               grad_rho = 0.0d0
            else
               call calc_lapl_dF_dlapl
            endif
#ifdef LIBXC
            if ( xc_name_exch == "mgga_x_tb09" ) then
               grad_rho = 0.0d0
            endif
#endif
         endif
! ================== 13.0XX

         do is = 1, ispin
            call finally_gga_xc_pot(is)
            ! dFxc/d(rho)  dF_drho + grad_rho --> chgrhr
         end do
      end if

      exc  = exc*univol*rinplw
      if(sw_output_xc_seperately==ON.or.xctype=='ggapbey')then
         eex  = eex*univol*rinplw
         ecor = ecor*univol*rinplw
      endif

      if(npes_cdfft >= 2) then               ! MPI
         if(myrank_ggacmp < nrank_ggacmp) then
            call mpi_allreduce(exc,exc_mpi,1,mpi_double_precision,mpi_sum &
                 & ,mpi_cdfft_world(myrank_ggacmp),ierr)  ! MPI
            exc = exc_mpi
            if(sw_output_xc_seperately==ON.or.xctype=='ggapbey')then
               call mpi_allreduce(eex,exc_mpi,1,mpi_double_precision,mpi_sum &
                    & ,mpi_cdfft_world(myrank_ggacmp),ierr)  ! MPI
               eex = exc_mpi
               call mpi_allreduce(ecor,exc_mpi,1,mpi_double_precision,mpi_sum &
                    & ,mpi_cdfft_world(myrank_ggacmp),ierr)  ! MPI
               ecor = exc_mpi
            endif
         end if
      end if

      if(nrest_cdfft >=1 ) then
         if(mype > npes_cdfft*nrank_ggacmp-1) then
            call mpi_recv(exc,1,mpi_double_precision &
                 &   , mype-npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,istatus,ierr)
            if(sw_output_xc_seperately==ON.or.xctype=='ggapbey')then
               call mpi_recv(eex,1,mpi_double_precision &
                    &   , mype-npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,istatus,ierr)
               call mpi_recv(ecor,1,mpi_double_precision &
                    &   , mype-npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,istatus,ierr)
            endif
         end if
         if(mype < nrest_cdfft) then
            call mpi_send(exc,1,mpi_double_precision &
                 &   , mype+npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,ierr)
            if(sw_output_xc_seperately==ON.or.xctype=='ggapbey')then
               call mpi_send(eex,1,mpi_double_precision &
                    &   , mype+npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,ierr)
               call mpi_send(ecor,1,mpi_double_precision &
                    &   , mype+npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,ierr)
            endif
         end if
      end if

!      if(xctype.eq.'vdwdf' .and. .not.oneshot) exc = exc + ecnl

      if(iprixc >= 2) write(nfout,'(" !XC exc = ",f20.8)') exc

      call tstatc0_end(id_sname)
    end subroutine ggaxcp

    subroutine ggaxcp_diff
      !        ~~~~~~~~~~~
      integer :: nfftcd, nfftx, nffty, nfftz
      real(kind=DP),allocatable,dimension(:) :: chgrhr_red
      real(kind=DP),allocatable,dimension(:) :: grad_rho_red
      real(kind=DP),allocatable,dimension(:,:) :: cgrad_rho_red

      real(DP):: exc_mpi               ! MPI
      integer :: is,i1,i2,in

!!$      call charge_in_gspace()  ! chgrhr -> chden_l
#ifdef LIBXC
      if ( xctype == "libxc" ) then
         call ggaxcp0_libxc()
      else
         call ggaxcp0()
      endif
#else
      call ggaxcp0()
#endif
      !    ~~~~~~~
             !   dF/d|rho(r)| (=vxc) --> dF_drho
             !   dFx/d|grad(rho(r))| --> dF_dgradrho
             !   dFc/d|grad(rho(r))| --> grad_trho
             !   ->exc

      exc = exc*univol*rinplw

      if(npes_cdfft >= 2) then               ! MPI
         if(myrank_ggacmp < nrank_ggacmp) then
            call mpi_allreduce(exc,exc_mpi,1,mpi_double_precision,mpi_sum &
                 & ,mpi_cdfft_world(myrank_ggacmp),ierr)  ! MPI
            exc = exc_mpi
         end if
      end if
      if(nrest_cdfft >=1 ) then
         if(mype > npes_cdfft*nrank_ggacmp-1) then
            call mpi_recv(exc,1,mpi_double_precision &
                 &   , mype-npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,istatus,ierr)
         end if
         if(mype < nrest_cdfft) then
            call mpi_send(exc,1,mpi_double_precision &
                 &   , mype+npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,ierr)
         end if
      end if
!!$      if(npes_cdfft >= 2) then               ! MPI
!!$         call mpi_allreduce(exc,exc_mpi,1,mpi_double_precision,mpi_sum &
!!$              & ,mpi_cdfft_world(myrank_ggacmp),ierr)  ! MPI
!!$         exc = exc_mpi
!!$      end if

!      if ( xctype .eq. 'vdwdf' ) exc = exc +ecnl

#ifndef DISABLE_VDWDF
      if ( xctype .eq. 'vdwdf' .and. ( .not. oneshot ) ) then
         nfftx = fft_box_size_CD(1,1);
         nffty = fft_box_size_CD(2,1);
         nfftz = fft_box_size_CD(3,1)
         nfftcd = nfftx*nffty*nfftz

         allocate(chgrhr_red(nfftcd));chgrhr_red=0.0d0
         allocate(grad_rho_red(nfftcd));grad_rho_red=0.d0
         allocate(cgrad_rho_red(nfftcd,3));cgrad_rho_red=0.d0

         if(ispin==1)then
            do i=ista_fftph,iend_fftph
               chgrhr_red(i) = chgrhr_l(i,1)
               grad_rho_red(i) = grad_rho(i,1)
               cgrad_rho_red(i,1:3) = cgrad_rho(i,1:3,1)
            enddo
         else
            do i=ista_fftph,iend_fftph
               chgrhr_red(i) = chgrhr_l(i,1)+chgrhr_l(i,2)
               grad_rho_red(i) = grad_rho(i,1)+grad_rho(i,2)
               cgrad_rho_red(i,1:3) = cgrad_rho(i,1:3,1)+cgrad_rho(i,1:3,2)
            enddo
         endif
         call mpi_allreduce(MPI_IN_PLACE,chgrhr_red,nfftcd,mpi_double_precision,mpi_sum,&
              & mpi_cdfft_world(myrank_ggacmp),ierr)
         call mpi_allreduce(MPI_IN_PLACE,grad_rho_red,nfftcd,mpi_double_precision,&
              & mpi_sum, mpi_cdfft_world(myrank_ggacmp),ierr)
         call mpi_allreduce(MPI_IN_PLACE,cgrad_rho_red,nfftcd*3, mpi_double_precision,&
              & mpi_sum, mpi_cdfft_world(myrank_ggacmp),ierr)
         call vdW_scf_stress( nspin, ispin, nfftx, nffty, nfftz, &
              &               chgrhr_red, grad_rho_red, cgrad_rho_red, &
              &               vdwdf_version, input_charge )
         deallocate( chgrhr_red )
         deallocate( grad_rho_red )
         deallocate( cgrad_rho_red )
      endif
#endif

      chgrhr_l = dF_drho

! stress tensor comes from gradient correction
!
! ==================================== added by K. Tagami ================ 11.0
      if ( noncol ) then
         if ( sw_constrain_on_grad_correction == ON ) then
            call calc_contrib_for_stress_nonclA
         else
            call calc_contrib_for_stress_nonclB
         endif
      else
         call calc_contrib_for_stress
      endif
! ========================================================================== 11.0
    end subroutine ggaxcp_diff

    subroutine calc_lapl_dF_dlapl
      integer :: iloop, n, i
      real(kind=DP) :: gx, gy, gz, g2
      real(kind=DP), allocatable :: afft(:)

      allocate(afft(ista_fftp:iend_fftp) )

      Do iloop=1, nspin
         if(myrank_ggacmp < nrank_ggacmp) then
            call cp_vals_to_afft_rspace1( iloop, dF_dlaplrho, afft )
            call m_FFT_CD_direct_c(nfout,afft)     ! afft(R space) -> afft(G sp.)
         else
            call m_FFT_CD_direct_c(nfout,afft,mode=0)
         end if
         do n = ista_fftp, iend_fftp  ! MPI
            gx = rltv(1,1)*inx(n) +rltv(1,2)*jnx(n) +rltv(1,3)*knx(n)
            gy = rltv(2,1)*inx(n) +rltv(2,2)*jnx(n) +rltv(2,3)*knx(n)
            gz = rltv(3,1)*inx(n) +rltv(3,2)*jnx(n) +rltv(3,3)*knx(n)
            g2 = gx**2 +gy**2 +gz**2
            afft(n) = -g2 *afft(n)
         enddo
         if(myrank_ggacmp < nrank_ggacmp) then
            call m_FFT_CD_inverse_c( nfout, afft )
         else
            call m_FFT_CD_inverse_c( nfout, afft, mode=0 )
         end if
         if(myrank_ggacmp < nrank_ggacmp) then
            do i = ista_fftph, iend_fftph    ! MPI
               dF_dlaplrho(i,iloop) = afft(i*2-1)
            end do
         else
            do i = ista_fftph, iend_fftph    ! MPI
               dF_dlaplrho(i,iloop) = 0.d0
            end do
         end if
      End Do
      deallocate( afft )
    end subroutine calc_lapl_dF_dlapl

! ==================================== added by K. Tagami ================ 11.0
    subroutine calc_contrib_for_stress
      integer :: is,i1,i2,in

      call allocate_for_stress
      s_gga1 = 0.0d0; s_gga2 = 0.0d0

      do i1=1,3
         do i2=1,3
            call rhos_diff(i1,i2);  call rhopc_diff(i1,i2);  call rhoh_diff(i1,i2)

            do is=1,nspin-af
               call map_drhodh1(is) ! drhodh -> chgfft
               do in=1,3
                  if(map_ggacmp(in) /= myrank_ggacmp) cycle
#ifdef _XC_SAVE_MEMORY_
                  cggawk13=0.d0
#endif
                  grad_rho(ista_fftph:iend_fftph,is)=0.d0
                  dF_drho(ista_fftph:iend_fftph,is)=0.d0
                  call dgrhodh1(i1,i2,in,is)
#ifndef _XC_SAVE_MEMORY_
                  call stress_exchange_part(i1,i2,in,is)   ! ->(s_gga1,s_gga2)
#else
                  call gradrho1(in,is)
                  call stress_exchange_part(i1,i2,is)   ! ->(s_gga1,s_gga2)
#endif
               enddo !in
            enddo !is

            if(nspin==2) then
               if(af == 1) then
                  is = 1
                  call map_drhodh1(is) ! drhodh -> chgfft
                  do in = 1, 3
                     if(map_ggacmp(in) /= myrank_ggacmp) cycle
                     call dgrhodh1(i1,i2,in,is)
#ifndef _XC_SAVE_MEMORY_
                     call stress_correlation_part(i1,i2,in)   ! ->(s_gga2)
#else
                     call gradrho2(in)
                     call stress_correlation_part(i1,i2)   ! ->(s_gga2)
#endif
                  end do
               else
                  call map_drhodh2
                  do in=1,3
                     if(map_ggacmp(in) /= myrank_ggacmp) cycle
                     call dgrhodh2(i1,i2,in)
#ifndef _XC_SAVE_MEMORY_
                     call stress_correlation_part(i1,i2,in)   ! ->(s_gga2)
#else
                     call gradrho2(in)
                     call stress_correlation_part(i1,i2)   ! ->(s_gga2)
#endif
                  end do
               end if
            end if
         enddo
      enddo

      call sum_s_gga12
      call deallocate_for_stress

    end subroutine calc_contrib_for_stress

    subroutine calc_contrib_for_stress_nonclA
      integer :: is,i1,i2,in

      call allocate_for_stress
      call alloc_chgsoft_project_noncol
      call set_chgsoft_projected_noncol
      call alloc_chghard_project_noncol
      call set_chghard_projected_noncol

      do i1=1,3
         do i2=1,3
            call rhos_diff_noncl(i1,i2);  call rhopc_diff_noncl(i1,i2)
            call rhoh_diff(i1,i2)

            do is=1,nspin
               call map_drhodh1(is) ! drhodh -> chgfft
               do in=1,3
                  if(map_ggacmp(in) /= myrank_ggacmp) cycle
#ifdef _XC_SAVE_MEMORY_
                  cggawk13=0.d0
#endif
                  grad_rho(ista_fftph:iend_fftph,is)=0.d0
                  dF_drho(ista_fftph:iend_fftph,is)=0.d0
                  call dgrhodh1(i1,i2,in,is)
#ifndef _XC_SAVE_MEMORY_
                  call stress_exchange_part(i1,i2,in,is)   ! ->(s_gga1,s_gga2)
#else
                  call gradrho1(in,is)
                  call stress_exchange_part(i1,i2,is)   ! ->(s_gga1,s_gga2)
#endif
               enddo !in
            enddo !is

            call map_drhodh2
            do in=1,3
               if(map_ggacmp(in) /= myrank_ggacmp) cycle
               call dgrhodh2(i1,i2,in)
#ifndef _XC_SAVE_MEMORY_
               call stress_correlation_part(i1,i2,in)   ! ->(s_gga2)
#else
               call gradrho2(in)
               call stress_correlation_part(i1,i2)   ! ->(s_gga2)
#endif
            end do

         enddo
      enddo

      call sum_s_gga12
      call deallocate_for_stress
      call dealloc_chgsoft_project_noncol
      call dealloc_chghard_project_noncol

    end subroutine calc_contrib_for_stress_nonclA

    subroutine calc_contrib_for_stress_nonclB
      integer :: is,i1,i2,in,nn

      call allocate_for_stress
      allocate( bfft_kt( ista_fftp:iend_fftp,ndim_magmom )); bfft_kt = 0.0d0

      do i1=1,3
         do i2=1,3
            call rhos_diff_noncl(i1,i2)
            call rhopc_diff_noncl(i1,i2)
            call rhoh_diff(i1,i2)                      ! drohdh is 1:ndim_magmom

            do in=1,3
               if (map_ggacmp(in) /= myrank_ggacmp) cycle

               bfft_kt = 0.0d0
               do is=1,ndim_magmom
                  call map_drhodh1(is)                  ! drhodh -> chgfft
                  call dgrhodh1_noncl(i1,i2,in,is)      ! chgfft, chden_l --> afft
                  bfft_kt(:,is) = afft(:)
               end do

               call m_ES_GradCorr_Along_QuantzAxis2( RhoMag_R, bfft_kt, &
                    &                                quantz_axis_inversion_flg_mesh )

! ====== KT_DEBUG ========== 2013/10/30
!               Do is=1, nspin
!                  dF_drho(:,is) = bfft_kt(:,is)
!               End do
!
               Do is=1, nspin
                  do nn = ista_fftph, iend_fftph              ! MPI
                     dF_drho(nn,is) = bfft_kt(2*nn,is)
                  enddo
               End Do
! =========================== 2013/10/30

#ifndef _XC_SAVE_MEMORY_
               Do is=1, nspin
                  call stress_exchange_part(i1,i2,in,is)   ! ->(s_gga1,s_gga2)
               End do
#else
               Do is=1, ndim_magmom
                  call gradrho1(in,is)
                  bfft_kt(:,is) = afft(:)
               End do

               call m_ES_GradCorr_Along_QuantzAxis2( RhoMag_R, bfft_kt, &
                    &                                quantz_axis_inversion_flg_mesh )

               Do is=1, nspin
                  cggawk13(:) = bfft_kt(:,is)
                  call stress_exchange_part(i1,i2,is)   ! ->(s_gga1,s_gga2)
               End do
#endif
            end do

            is = 1
            call map_drhodh1(is)                  ! drhodh -> chgfft

            do in=1,3
               if(map_ggacmp(in) /= myrank_ggacmp) cycle

               call dgrhodh1(i1,i2,in,is)
#ifndef _XC_SAVE_MEMORY_
               call stress_correlation_part(i1,i2,in)   ! ->(s_gga2)
#else
               call gradrho1(in,is)
               call stress_correlation_part(i1,i2)   ! ->(s_gga2)
#endif
            end do

         enddo
      enddo

      call sum_s_gga12

      call deallocate_for_stress
      deallocate( bfft_kt )

    end subroutine calc_contrib_for_stress_nonclB
! ======================================================================= 11.0

    subroutine ggaxcp0
      integer :: pot_type

! ==== KT_add ======== 13.0XX
      real(kind=DP), allocatable :: dummy1(:,:), dummy2(:)
! ==================== 13.0XX

      call initialize_cggawk_arrays

!!$      if(myrank_ggacmp >= nrank_ggacmp) goto 1002

      if(xctype == 'ldapw91' .or. xctype == 'ldapbe ') then
         grad_rho = 0.d0; grad_trho = 0.d0
#ifndef _XC_SAVE_MEMORY_
         cgrad_rho = 0.d0
#endif
      else
! ===================================== modiifed by K. Tagami ============= 11.0
!!         call abs_grad_rho_up_down_total()
      !   chden_l ->
      !   |grad(rho_up(r))|,|grad(rho_down(r))|  --> grad_rho
      !   |grad(rho_up(r) + rho_down(r))|        --> grad_trho
         if ( noncol .and. sw_constrain_on_grad_correction == OFF ) then
            if ( input_charge == Partial_Core_Charge .and. vflag == EXC_ONLY ) then
               call abs_grad_rho_up_down_total()
            else
               call abs_grad_rho_updown_total_noncl()
            endif
         else
            call abs_grad_rho_up_down_total()
         endif
! ========================================================================== 11.0
      end if

! ====== KT_add ======= 13.0XX
      if ( use_metagga ) then
         call grad2_rho_up_down()

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
! ==================================  KT_add ================ 13.0A
                 & .and. xctype /= 'rpbe' .and. xctype /= 'wc06' &
                 & .and. xctype /= 'htbs' .and. xctype /= 'pbesol' &
                 & .and. xctype /= 'pbeint' &
                 & .and. xctype /= 'ev93' .and. xctype /= 'evpw91' &
                 & .and. xctype /= 'lb94' &
! =========================================================== 13.0A
                 & .and. xctype /= 'revpbe' &
                 & .and. xctype /= 'ggapbe'.and. xctype /= 'ggapbex' &
#ifdef DISABLE_VDWDF
                 & .and. xctype /= 'ggapbey') then
#else
                 & .and. xctype /= 'ggapbey' &
                 & .and. xctype /= 'vdwdf' ) then
#endif
            write(nfout,'(" xctype = ",a7)') xctype
            call phase_error_with_msg(nfout,' xctype is not set properly (ggaxcp0)',__LINE__,__FILE__)
         end if
      endif
! ===================== 13.0XX

#ifdef __EDA__
      if(sw_eda==ON) then
      allocate(exc_on_a_grid_wk(ista_fftph:iend_fftph))
      else
      allocate(exc_on_a_grid_wk(1))
      endif
#endif

      if(myrank_ggacmp < nrank_ggacmp) then
      if(xctype == 'ggapw91' .or. xctype == 'ldapw91') then
#ifdef __EDA__
         call ex_ggapw91( nspin, ispin, ista_fftph, iend_fftph, &
              &           chgrhr_l, grad_rho, f2or1, exc, dF_drho, dF_dgradrho, &
              &           exc_on_a_grid_wk, ista_fftph, iend_fftph )
#else
         call ex_ggapw91( nspin, ispin, ista_fftph, iend_fftph, &
              &           chgrhr_l, grad_rho, f2or1, exc, dF_drho, dF_dgradrho, &
              &           ista_fftph, iend_fftph )
#endif
         if(sw_hybrid_functional==ON) call scale_exchange(alpha_exx)
#ifdef __EDA__
         if(sw_exchange_only==OFF) &
              & call cr_ggapw91( nspin, ispin, ista_fftph, iend_fftph, &
              &                  chgrhr_l, grad_trho, f2or1, exc, &
              &                  dF_drho, exc_on_a_grid_wk, ista_fftph, iend_fftph )
#else
         if(sw_exchange_only==OFF) &
              & call cr_ggapw91( nspin, ispin, ista_fftph, iend_fftph, &
              &                  chgrhr_l, grad_trho, f2or1, exc, &
              &                  dF_drho, ista_fftph, iend_fftph )
#endif
      else if(xctype == 'ggapbe ' .or. xctype == 'ldapbe ' &
           &                      .or. xctype=='ggapbey'.or.xctype=='revpbe') then
#ifdef __EDA__
         call ex_ggapbe( nspin, ispin, ista_fftph, iend_fftph, &
              &          chgrhr_l, grad_rho, f2or1, exc, dF_drho, dF_dgradrho, &
         &               exc_on_a_grid_wk, xctype=='revpbe', ista_fftph ,iend_fftph )
#else
         call ex_ggapbe( nspin, ispin, ista_fftph, iend_fftph, &
              &          chgrhr_l, grad_rho, f2or1, exc, dF_drho, dF_dgradrho, &
              &          xctype=='revpbe', ista_fftph, iend_fftph )
#endif
         eex = exc

         if(sw_hybrid_functional==ON) then
            if ( hybrid_functional_type == 'gaupbe' ) then
               call ex_gaupbe( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
              &               grad_rho, f2or1, exc, dF_drho, dF_dgradrho )
            else if ( hybrid_functional_type == 'hse06-hjs' ) then
               call ex_omegapbe_library_HJS( -alpha_exx, omega_exx_pbe, &
                    &                        nspin, ispin, ista_fftph, iend_fftph, &
                    &                        chgrhr_l, grad_rho, f2or1, exc, &
                    &                        dF_drho, dF_dgradrho, 1 )
            else
               if(sw_screened_exchange==ON) then
                  call ex_omegapbe( -alpha_exx, omega_exx_pbe, nspin, ispin, &
                       &            ista_fftph, iend_fftph, chgrhr_l, grad_rho, &
                       &            f2or1, exc, dF_drho, dF_dgradrho )
               else
                  call scale_exchange(alpha_exx)
               end if
            endif
         end if
#ifdef __EDA__
         if(sw_exchange_only==OFF) &
              & call cr_ggapbe( nspin, ispin, ista_fftph, iend_fftph, &
              &                 chgrhr_l, grad_trho, f2or1, exc, dF_drho, &
              &                 exc_on_a_grid_wk, ecor, ista_fftph, iend_fftph )
#else
         if(sw_exchange_only==OFF) &
              & call cr_ggapbe( nspin, ispin, ista_fftph, iend_fftph, &
              &                 chgrhr_l, grad_trho, f2or1, exc, dF_drho, &
              &                 ecor, ista_fftph, iend_fftph )
#endif

      else if(xctype == 'ggapbex') then
#ifdef __EDA__
         call ex_ggapbe( nspin, ispin, ista_fftph, iend_fftph, &
              &          chgrhr_l, grad_rho, f2or1, exc, dF_drho, dF_dgradrho, &
         &               exc_on_a_grid_wk, xctype=='revpbe', ista_fftph, iend_fftph )
#else
         call ex_ggapbe( nspin, ispin, ista_fftph, iend_fftph, &
              &          chgrhr_l, grad_rho, f2or1, exc, dF_drho, dF_dgradrho, &
              &          xctype=='revpbe', ista_fftph, iend_fftph )
#endif

! =========================== KT_add ========================= 13.0A
      else if ( xctype == 'rpbe' .or. xctype == 'wc06' .or. xctype == 'htbs' &
           &                     .or. xctype == 'pbesol' .or. xctype == 'pbeint' &
           &                     .or. xctype == 'ev93' .or. xctype == 'evpw91' &
           &                     .or. xctype == 'lb94' ) then
         if ( xctype == 'ggapbe' )  pot_type = 1
         if ( xctype == 'revpbe' )  pot_type = 2
         if ( xctype == 'rpbe' )    pot_type = 3
         if ( xctype == 'wc06' )    pot_type = 4
         if ( xctype == 'htbs' )    pot_type = 5

         if ( xctype == 'pbesol' )    pot_type = 6
         if ( xctype == 'pbeint' )    pot_type = 7

         if ( xctype == 'ev93' )    pot_type = 20
         if ( xctype == 'evpw91' )  pot_type = 20
         if ( xctype == 'lb94' )  pot_type = 25

#ifdef __EDA__
         call ex_gga_library( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
              &               grad_rho, f2or1, exc, dF_drho, dF_dgradrho, &
              &               exc_on_a_grid_wk, pot_type, ista_fftph, iend_fftph )
#else
         call ex_gga_library( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
              &               grad_rho, f2or1, exc, dF_drho, dF_dgradrho, pot_type, &
              &               ista_fftph, iend_fftph )
#endif

         eex = exc
         if ( sw_hybrid_functional == ON ) then        ! pbe-only
            if ( hybrid_functional_type == 'hsesol-hjs' ) then
               call ex_omegapbe_library_HJS( -alpha_exx, omega_exx_pbe, &
                    &                        nspin, ispin, ista_fftph, iend_fftph, &
                    &                        chgrhr_l, grad_rho, f2or1, exc, &
                    &                        dF_drho, dF_dgradrho, 6 )
            else
               if ( sw_screened_exchange == ON ) then
                  call ex_omegapbe( -alpha_exx, omega_exx_pbe, nspin, ispin, &
                       &            ista_fftph, iend_fftph, chgrhr_l, &
                       &            grad_rho, f2or1, exc, dF_drho, dF_dgradrho )
               else
                  call scale_exchange( alpha_exx )
               end if
            endif
         end if
         if ( sw_exchange_only == OFF ) then
#ifdef __EDA__
            if ( xctype == 'ev93' .or. xctype == 'lb94' ) then
               grad_trho=0.d0
               call cr_gga_library( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
                    &               grad_trho, f2or1, exc, dF_drho, exc_on_a_grid_wk, &
                    &               ecor, pot_type, ista_fftph, iend_fftph )
            else if ( xctype == 'evpw91' ) then
               call cr_ggapw91( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
                   &            grad_trho, f2or1, exc, dF_drho, exc_on_a_grid_wk, &
                   &            ista_fftph, iend_fftph )
            else
               call cr_gga_library( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
                    &               grad_trho, f2or1, exc, dF_drho, exc_on_a_grid_wk, &
                    &               ecor, pot_type, ista_fftph, iend_fftph )
            endif
#else
            if ( xctype == 'ev93' .or. xctype == 'lb94' ) then
               grad_trho=0.d0
               call cr_gga_library( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
                    &               grad_trho, f2or1, exc, dF_drho, ecor, pot_type, &
                    &               ista_fftph, iend_fftph )
            else if ( xctype == 'evpw91' ) then
               call cr_ggapw91( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
                   &            grad_trho, f2or1, exc, dF_drho, &
                   &            ista_fftph, iend_fftph )
            else
               call cr_gga_library( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
                    &               grad_trho, f2or1, exc, dF_drho, ecor, pot_type, &
                    &               ista_fftph, iend_fftph )
            endif
#endif
         endif
! ============================================================= 13.0A
#ifndef DISABLE_VDWDF
      else if(xctype == 'vdwdf') then
#ifdef __EDA__
         if ( exchange_pot_type == 'pbe' ) then
            call ex_ggapbe( nspin, ispin, ista_fftph, iend_fftph, &
                 &          chgrhr_l, grad_rho, f2or1, exc, dF_drho, dF_dgradrho, &
                 &          exc_on_a_grid_wk, .false., ista_fftph, iend_fftph )
         else if ( exchange_pot_type == 'revpbe' ) then
            call ex_ggapbe( nspin, ispin, ista_fftph, iend_fftph, &
                 &          chgrhr_l, grad_rho, f2or1, exc, dF_drho, dF_dgradrho, &
                 &          exc_on_a_grid_wk, .true., ista_fftph, iend_fftph )
         else
            if ( exchange_pot_type == 'b86r' )  pot_type = 11
            if ( exchange_pot_type == 'optpbe' )  pot_type = 12
            if ( exchange_pot_type == 'optb86b' )  pot_type = 13
            if ( exchange_pot_type == 'pw86r' )  pot_type = 14
            if ( exchange_pot_type == 'c09x' )  pot_type = 15
            if ( exchange_pot_type == 'lvpw86r' )  pot_type = 16

            call ex_gga_library( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
                 &               grad_rho, f2or1, exc, dF_drho, dF_dgradrho, &
                 &               exc_on_a_grid_wk, pot_type, ista_fftph, iend_fftph )
         endif

         eex = exc
         if(sw_hybrid_functional==ON) then
            if(sw_screened_exchange==ON) then
               call ex_omegapbe( -alpha_exx, omega_exx_pbe, nspin, ispin, &
                    &            ista_fftph, iend_fftph, chgrhr_l, &
                    &            grad_rho, f2or1, exc, dF_drho, dF_dgradrho )
            else
               call scale_exchange(alpha_exx)
            end if
         end if

         if(sw_exchange_only==OFF) then
            grad_trho=0.d0
            call cr_ggapbe( nspin, ispin, ista_fftph, iend_fftph, &
                 &          chgrhr_l, grad_trho, f2or1, exc, dF_drho, &
                 &          exc_on_a_grid_wk, ecor, ista_fftph, iend_fftph )
         endif
         if (.not. oneshot ) call add_vdwdf_nonlocal_energy( vdwdf_version )
#else
         if ( exchange_pot_type == 'pbe' ) then
            call ex_ggapbe( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, grad_rho, &
                 &          f2or1, exc, dF_drho, dF_dgradrho, .false. ,&
                 &          ista_fftph, iend_fftph )
         else if ( exchange_pot_type == 'revpbe' ) then
            call ex_ggapbe( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, grad_rho, &
                 &          f2or1, exc, dF_drho, dF_dgradrho, .true., &
                 &          ista_fftph, iend_fftph )
         else
            if ( exchange_pot_type == 'b86r' )  pot_type = 11
            if ( exchange_pot_type == 'optpbe' )  pot_type = 12
            if ( exchange_pot_type == 'optb86b' )  pot_type = 13
            if ( exchange_pot_type == 'pw86r' )  pot_type = 14
            if ( exchange_pot_type == 'c09x' )  pot_type = 15
            if ( exchange_pot_type == 'lvpw86r' )  pot_type = 16

            call ex_gga_library( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
                 &               grad_rho, f2or1, exc, dF_drho, dF_dgradrho, pot_type, &
                 &               ista_fftph, iend_fftph )
         endif

         eex = exc
         if(sw_hybrid_functional==ON) then
            if(sw_screened_exchange==ON) then
               call ex_omegapbe( -alpha_exx, omega_exx_pbe, nspin, ispin, &
                    &            ista_fftph, iend_fftph, chgrhr_l, &
                    &            grad_rho, f2or1, exc, dF_drho, dF_dgradrho )
            else
               call scale_exchange(alpha_exx)
            end if
         end if

         if(sw_exchange_only==OFF) then
            grad_trho=0.d0
            call cr_ggapbe( nspin, ispin, ista_fftph, iend_fftph, &
                 &          chgrhr_l, grad_trho, f2or1, exc, dF_drho,ecor, &
                 &          ista_fftph, iend_fftph )
         endif
         if (.not. oneshot ) call add_vdwdf_nonlocal_energy( vdwdf_version )
#endif
#endif
      else if(xctype == 'katopbe' .or. xctype == 'ggapbek') then
#ifdef __EDA__
         call ex_ggapbe( nspin, ispin, ista_fftph, iend_fftph, &
              &          chgrhr_l, grad_rho, f2or1, exc, dF_drho, dF_dgradrho, &
         &               exc_on_a_grid_wk, xctype=='revpbe', ista_fftph, iend_fftph )
#else
         call ex_ggapbe( nspin, ispin, ista_fftph, iend_fftph, &
              &          chgrhr_l, grad_rho, f2or1, exc, dF_drho, dF_dgradrho, &
              &          xctype=='revpbe', ista_fftph, iend_fftph )
#endif
         if(sw_hybrid_functional==ON) then
            if(sw_screened_exchange==ON) then
               call ex_omegapbe( -alpha_exx, omega_exx_pbe, nspin, ispin, &
                    &            ista_fftph, iend_fftph, chgrhr_l, grad_rho, &
                    &            f2or1,exc, dF_drho, dF_dgradrho )
            else
               call scale_exchange(alpha_exx)
            end if
         end if
#ifdef __EDA__
         if(sw_exchange_only==OFF) &
              & call cr_ggapbe( nspin, ispin, ista_fftph, iend_fftph, &
              &                 chgrhr_l, grad_trho, f2or1, exc, dF_drho, &
              &                 exc_on_a_grid_wk, ecor, ista_fftph, iend_fftph )
#else
         if(sw_exchange_only==OFF) &
              & call cr_ggapbe( nspin, ispin, ista_fftph, iend_fftph, &
              &                 chgrhr_l, grad_trho, f2or1, exc, dF_drho, &
              &                 ecor, ista_fftph, iend_fftph )
#endif
      else if(xctype == 'ggabp  ') then
#ifdef __EDA__
         call xclda(nspin,ispin,ista_fftph,iend_fftph,chgrhr_l,f2or1,exc,dF_drho, exc_on_a_grid_wk)
         call ggabek(nspin,ispin,ista_fftph,iend_fftph,chgrhr_l,grad_rho,f2or1,exc,dF_drho,dF_dgradrho, exc_on_a_grid_wk)
         call ggaprd(nspin,ispin,ista_fftph,iend_fftph,chgrhr_l,grad_rho,f2or1,exc,dF_drho,dF_dgradrho, exc_on_a_grid_wk)
#else
         call xclda(nspin,ispin,ista_fftph,iend_fftph,chgrhr_l,f2or1,exc,dF_drho)
         call ggabek(nspin,ispin,ista_fftph,iend_fftph,chgrhr_l,grad_rho,f2or1,exc,dF_drho,dF_dgradrho)
         call ggaprd(nspin,ispin,ista_fftph,iend_fftph,chgrhr_l,grad_rho,f2or1,exc,dF_drho,dF_dgradrho)
#endif
      end if
      else
        exc = 0.d0
      end if

#ifdef __EDA__
      if(sw_eda==ON) then
      call cp_exc_for_EDA
      endif
      deallocate(exc_on_a_grid_wk)
#endif

             !   dF/d|rho(r)| (=vxc) --> dF_drho
             !   dFx/d|grad(rho(r))| --> dF_dgradrho
             !   dFc/d|grad(rho(r))| --> grad_trho
!!$1002  continue

! ====== KT_add ===== 13.0XX
      if ( use_metagga ) then
         if ( xctype == "tb09" ) then
            if ( sw_fix_val_c_tb09 == OFF ) then
               call set_gval_tb09( nspin, ista_fftph, iend_fftph, &
                    &              chgrhr_l, grad_rho, f2or1, val_g )
               val_g = val_g +val_g_tb09_paw
               call set_cval_tb09( val_g, val_c_tb09 )
            endif
            write(nfout,*) "val_c_tb09 is now ", val_c_tb09
         endif

         if (myrank_ggacmp < nrank_ggacmp ) then

            if ( xctype == "tb09" ) then
               allocate( dummy1(ista_fftph:iend_fftph,nspin ) ); dummy1 = 0.0d0
               allocate( dummy2(ista_fftph:iend_fftph ) );       dummy2 = 0.0d0

#ifdef __EDA__
               call ex_ggapw91( nspin, ispin, ista_fftph, iend_fftph, &
                    &           chgrhr_l, dummy1, f2or1, exc, dF_drho, dF_dgradrho, exc_on_a_grid_wk, ista_fftph, iend_fftph )
#else
               call ex_ggapw91( nspin, ispin, ista_fftph, iend_fftph, &
                    &           chgrhr_l, dummy1, f2or1, exc, dF_drho, dF_dgradrho, &
                    &           ista_fftph, iend_fftph )
#endif

               if ( ekin_density_is_active ) then
                  call ex_mgga_tb09( nspin, ispin, ista_fftph, iend_fftph, 1, &
                       &              chgrhr_l, grad_rho, grad2_rho, ekin_dens, &
                       &              dF_drho, val_c_tb09 )
               endif
#ifdef __EDA__
               call cr_ggapw91( nspin, ispin, ista_fftph, iend_fftph, &
                    &           chgrhr_l, dummy2,  f2or1, exc, dF_drho, &
                    &           exc_on_a_grid_wk, ista_fftph, iend_fftph )
#else
               call cr_ggapw91( nspin, ispin, ista_fftph, iend_fftph, &
                    &           chgrhr_l, dummy2,  f2or1, exc, dF_drho, &
                    &           ista_fftph, iend_fftph )
#endif

               deallocate( dummy1 ); deallocate( dummy2 )
               dF_dgradrho = 0.0d0
            endif
         end if

      endif
! =================== 13.0XX

      if(iprixc >= 2) then
         write(nfout,'(" !XC -- chgrhr_l -- <<ggaxcp0>>")')
         write(nfout,'(" !XC ",6f12.6)') (chgrhr_l(i,1),i=ista_fftph,ista_fftph+17)
         write(nfout,'(" !XC -- grad_rho -- <<ggaxcp0>>")')
         write(nfout,'(" !XC ",6f12.6)') (grad_rho(i,1),i=ista_fftph,ista_fftph+17)
         write(nfout,'(" !XC -- grad_trho -- <<ggaxcp0>>")')
         write(nfout,'(" !XC ",6f12.6)') (grad_trho(i),i=ista_fftph,ista_fftph+17)
         write(nfout,'(" !XC -- dF_drho -- <<ggaxcp0>>")')
         write(nfout,'(" !XC ",6f12.6)') (dF_drho(i,1),i=ista_fftph,ista_fftph+17)
         write(nfout,'(" !XC -- dF_dgradrho -- <<ggaxcp0>>")')
         write(nfout,'(" !XC ",6f12.6)') (dF_dgradrho(i,1),i=ista_fftph,ista_fftph+17)
      end if
    end subroutine ggaxcp0

#ifdef LIBXC
    subroutine ggaxcp0_libxc
      integer :: i, key_save

! ==== KT_add ======== 13.0XX
      real(kind=DP), allocatable :: dummy1(:,:), dummy2(:)
! ==================== 13.0XX

      call initialize_cggawk_arrays
      exc = 0.0d0

!!$      if(myrank_ggacmp >= nrank_ggacmp) goto 1002

      if ( xc_family_exch == XC_FAMILY_LDA .and. xc_family_corr == XC_FAMILY_LDA  ) then
         grad_rho = 0.d0; grad_trho = 0.d0
#ifndef _XC_SAVE_MEMORY_
         cgrad_rho = 0.d0
#endif
      else
         if ( noncol .and. sw_constrain_on_grad_correction == OFF ) then
            if ( input_charge == Partial_Core_Charge .and. vflag == EXC_ONLY ) then
               call abs_grad_rho_up_down_total()
            else
               call abs_grad_rho_updown_total_noncl()
            endif
         else
            call abs_grad_rho_up_down_total()
         endif
      end if

! ====== KT_add ======= 13.0XX
      if ( use_metagga ) then
         call grad2_rho_up_down()

#if 0
         if ( use_modeled_ekin_density ) then
            call m_KE_set_modeled_ekin_density( ista_fftph, iend_fftph, &
                 &                              chgrhr_l, grad_rho, grad2_rho, &
                 &                              ekin_dens )
         else if ( use_asymm_ekin_density ) then
            ekin_dens = ekin_dens + 0.25d0 *grad2_rho       ! ????
         endif
#else
!         if ( use_modeled_ekin_density ) then
         if ( use_modeled_ekin_density .or. ( icond==0 .and. iteration == 0 ) ) then
            if ( iteration == 0 ) then
               key_save = ekin_density_type;  ekin_density_type = 3
            endif
            call m_KE_set_modeled_ekin_density( ista_fftph, iend_fftph, &
                 &                              chgrhr_l, grad_rho, grad2_rho, &
                 &                              ekin_dens )
            if ( iteration == 0 ) ekin_density_type = key_save
         else if ( use_asymm_ekin_density ) then
            ekin_dens = ekin_dens + 0.25d0 *grad2_rho       ! ????
         endif
#endif
      endif

#ifdef __EDA__
      if(sw_eda==ON) then
         allocate(exc_on_a_grid_wk(ista_fftph:iend_fftph))
      else
         allocate(exc_on_a_grid_wk(1))
      endif
#endif

      if (myrank_ggacmp < nrank_ggacmp) then
         select case (xc_family_exch)
         case (XC_FAMILY_LDA)
#ifdef __EDA__
            call ex_lda_libxc( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
                 &             f2or1, exc, dF_drho, exc_on_a_grid_wk, &
                 &             ista_fftph, iend_fftph )
#else
            call ex_lda_libxc( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
                 &             f2or1, exc, dF_drho, ista_fftph, iend_fftph )
#endif
         case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
#ifdef __EDA__
            call ex_gga_libxc( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
                 &             grad_rho, grad_trho, f2or1, exc, dF_drho, dF_dgradrho, &
                 &             exc_on_a_grid_wk, ista_fftph, iend_fftph )
#else
            call ex_gga_libxc( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
                 &             grad_rho, grad_trho, f2or1, exc, dF_drho, dF_dgradrho, &
                 &             ista_fftph, iend_fftph )
#endif
         case (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
            if ( xc_name_exch == "mgga_x_tb09" ) then
               if ( sw_fix_val_c_tb09 == OFF ) then
                  call set_gval_tb09( nspin, ista_fftph, iend_fftph, &
                       &              chgrhr_l, grad_rho, f2or1, val_g )
                  val_g = val_g +val_g_tb09_paw
                  call set_cval_tb09( val_g, val_c_tb09 )
               endif
               write(nfout,*) "val_c_tb09 is now ", val_c_tb09
            endif
            call ex_mgga_libxc( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
                 &              grad_rho, grad_trho, grad2_rho, ekin_dens, f2or1, exc, &
                 &              dF_drho, dF_dgradrho, dF_dlaplrho, dF_dtau )
         end select

         select case (xc_family_corr)
         case (XC_FAMILY_LDA)
#ifdef __EDA__
            call cr_lda_libxc( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
                 &             f2or1, exc, dF_drho, &
                 &             exc_on_a_grid_wk, ista_fftph, iend_fftph )
#else
            call cr_lda_libxc( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
                 &             f2or1, exc, dF_drho, ista_fftph, iend_fftph )
#endif
         case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
#ifdef __EDA__
            call cr_gga_libxc( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
                 &             grad_rho, grad_trho, f2or1, exc, dF_drho, grad_trho, &
                 &             exc_on_a_grid_wk, ista_fftph, iend_fftph )
#else
            call cr_gga_libxc( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
                 &             grad_rho, grad_trho, f2or1, exc, dF_drho, grad_trho, &
                 &             ista_fftph, iend_fftph )
#endif
         case (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
            call cr_mgga_libxc( nspin, ispin, ista_fftph, iend_fftph, chgrhr_l, &
                 &              grad_rho, grad_trho, grad2_rho, ekin_dens, f2or1, exc, &
                 &              dF_drho, grad_trho, dF_dlaplrho, dF_dtau )
         end select
      end if
#ifdef __EDA__
      if(sw_eda==ON) then
         call cp_exc_for_EDA
      endif
      deallocate(exc_on_a_grid_wk)
#endif

    end subroutine ggaxcp0_libxc
#endif


#ifndef DISABLE_VDWDF
    subroutine add_vdwdf_nonlocal_energy( version_vdwdf )
      integer, intent(in) :: version_vdwdf
      integer :: i, nfftcd,nfftx,nffty,nfftz,ix,iy,iz,ixyz
      integer :: ix2, iy2, iz2, nlphf,idp, mmp

      real(kind=DP),allocatable,dimension(:) :: chgrhr_red
      real(kind=DP),allocatable,dimension(:) :: grad_rho_red
      real(kind=DP),allocatable,dimension(:,:) :: dF_drho_red
      real(kind=DP),allocatable,dimension(:,:,:) :: dfdrho_vdw,dfddrho_vdw

      idp = fft_box_size_CD_nonpara(1,0)
      mmp = fft_box_size_CD_nonpara(2,0)

      nfftx = fft_box_size_CD(1,1)
      nffty = fft_box_size_CD(2,1)
      nfftz = fft_box_size_CD(3,1)
      nfftcd = nfftx*nffty*nfftz

      allocate(dfdrho_vdw (nfftx,nffty,nfftz));dfdrho_vdw=0.d0
      allocate(dfddrho_vdw(nfftx,nffty,nfftz));dfddrho_vdw=0.d0
      allocate(chgrhr_red(nfftcd));chgrhr_red=0.0d0
      allocate(grad_rho_red(nfftcd));grad_rho_red=0.d0

      if(ispin==1)then
         do i=ista_fftph,iend_fftph
            chgrhr_red(i)   = chgrhr_l(i,1)
            grad_rho_red(i) = grad_rho(i,1)
         enddo
      else
         do i=ista_fftph,iend_fftph
            chgrhr_red(i)   = chgrhr_l(i,1) +chgrhr_l(i,2)
            grad_rho_red(i) = grad_rho(i,1) +grad_rho(i,2)
         enddo
      endif
      call mpi_allreduce( MPI_IN_PLACE, chgrhr_red, nfftcd, mpi_double_precision, &
           &              mpi_sum, mpi_cdfft_world(myrank_ggacmp), ierr )
      call mpi_allreduce( MPI_IN_PLACE, grad_rho_red, nfftcd, mpi_double_precision,&
           &              mpi_sum, mpi_cdfft_world(myrank_ggacmp), ierr )
      call vdW_scf( nspin, ispin, nfftx, nffty, nfftz, chgrhr_red, grad_rho_red, &
           &        ecnl, dfdrho_vdw, dfddrho_vdw, version_vdwdf, vflag )

      if(kimg == 1) then
         nlphf = idp/2
      else
         nlphf = idp
      end if

      if ( kimg == 1 ) then
        do ix=1,nfftx
            do iy=1,nffty
               do iz=1,nfftz
                  if ( ix > nlphf ) then
                     ix2 = idp -ix
                     iy2 = nffty +2 -iy
                     iz2 = nfftz +2 -iz
                     if ( iy2 > nffty ) iy2 = iy2 -nffty
                     if ( iz2 > nfftz ) iz2 = iz2 -nfftz
                     cycle
                  else
                     ix2 = ix;  iy2 = iy;   iz2 = iz
                  endif
                  ixyz = (iz2-1)*mmp*nlphf +(iy2-1)*nlphf +ix2

                  if(ixyz<ista_fftph.or.ixyz>iend_fftph) cycle

                  dF_drho(ixyz,1:ispin) = dF_drho(ixyz,1:ispin) &
                       &                  +dfdrho_vdw(ix,iy,iz)
                  dF_dgradrho(ixyz,1:ispin) = dF_dgradrho(ixyz,1:ispin) &
                       &                     +dfddrho_vdw(ix,iy,iz)
               enddo
            enddo
         enddo
      else
         do ix=1, nfftx
            do iy=1,nffty
               do iz=1,nfftz
                  ixyz = (iz-1)*nffty*nfftx +(iy-1)*nfftx +ix
                  if ( ixyz<ista_fftph .or. ixyz>iend_fftph ) cycle

                  dF_drho(ixyz,1:ispin) = dF_drho(ixyz,1:ispin) &
                       &                  +dfdrho_vdw(ix,iy,iz)
                  dF_dgradrho(ixyz,1:ispin) = dF_dgradrho(ixyz,1:ispin) &
                       &                     +dfddrho_vdw(ix,iy,iz)
               enddo
            enddo
         enddo
      endif

      deallocate(dfdrho_vdw); deallocate(dfddrho_vdw)
      deallocate(chgrhr_red); deallocate(grad_rho_red)
    end subroutine add_vdwdf_nonlocal_energy
#endif

    subroutine scale_exchange(alpha)
      implicit none
      real(kind=DP), intent(in) :: alpha

      real(kind=DP) :: fexx

      fexx = 1.d0-alpha
      exc = fexx * exc
      dF_drho = fexx * dF_drho
      dF_dgradrho = fexx * dF_dgradrho
    end subroutine scale_exchange

    subroutine sum_s_gga12
      real(kind=DP), allocatable, dimension(:,:)    :: s_gga_mpi        ! MPI d(3,3)

      if(npes>=2 .and. (nrank_ggacmp>1.or.npes_cdfft>1)) allocate(s_gga_mpi(3,3))
      if(nrank_ggacmp > 1) then
         call mpi_allreduce(s_gga1,s_gga_mpi,9,mpi_double_precision &
              & , mpi_sum,mpi_ggacmp_cross_world(myrank_cdfft),ierr)
         s_gga1 = s_gga_mpi

         call mpi_allreduce(s_gga2,s_gga_mpi,9,mpi_double_precision &
              & , mpi_sum,mpi_ggacmp_cross_world(myrank_cdfft),ierr)
         s_gga2 = s_gga_mpi
      end if

      if(npes_cdfft >= 2)  then
         if(myrank_ggacmp < nrank_ggacmp) then
            call mpi_allreduce(s_gga1,s_gga_mpi,9,mpi_double_precision &
                 & ,mpi_sum,mpi_cdfft_world(myrank_ggacmp),ierr)
            s_gga1 = s_gga_mpi
            call mpi_allreduce(s_gga2,s_gga_mpi,9,mpi_double_precision &
                 & ,mpi_sum,mpi_cdfft_world(myrank_ggacmp),ierr)
            s_gga2 = s_gga_mpi
         end if
      end if
      if(nrest_cdfft >= 1) then
         if(mype > npes_cdfft*nrank_ggacmp - 1) then
            call mpi_recv(s_gga1,9,mpi_double_precision &
                 & ,mype - npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,istatus,ierr)
         end if
         if(mype < nrest_cdfft) then
            call mpi_send(s_gga1,9,mpi_double_precision &
                 & ,mype + npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,ierr)
         end if
         if(mype > npes_cdfft*nrank_ggacmp - 1) then
            call mpi_recv(s_gga2,9,mpi_double_precision &
                 & ,mype - npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,istatus,ierr)
         end if
         if(mype < nrest_cdfft) then
            call mpi_send(s_gga2,9,mpi_double_precision &
                 & ,mype + npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,ierr)
         end if
      end if

!!$         call mpi_allreduce(s_gga1,s_gga_mpi,9 &
!!$              &, mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)   ! MPI
!!$         s_gga1 = s_gga_mpi                  ! MPI
!!$
!!$         call mpi_allreduce(s_gga2,s_gga_mpi,9 &
!!$              &, mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)   ! MPI
!!$         s_gga2 = s_gga_mpi                  ! MPI

      if(npes >= 2 .and. (nrank_ggacmp>1.or.npes_cdfft>1)) deallocate(s_gga_mpi)

    end subroutine sum_s_gga12

#ifndef _XC_SAVE_MEMORY_

    subroutine stress_exchange_part(i1,i2,in,is)
      integer, intent(in) :: i1,i2,in,is
      integer :: m

      do m = ista_fftph, iend_fftph       ! MPI
         s_gga1(i1,i2)=s_gga1(i1,i2) &
              &  - univol * rinplw * dF_dgradrho(m,is) &
              &    * cgrad_rho(m,in,is) * dF_drho(m,is) * f2or1(m)
      end do

      ! correlation part -->
      if(nspin.eq.1) then
         do m = ista_fftph, iend_fftph    ! MPI
            s_gga2(i1,i2)=s_gga2(i1,i2) &
                 &  - univol * rinplw * grad_trho(m) &
                 &    * cgrad_rho(m,in,is) * dF_drho(m,1) * f2or1(m)
         end do
      endif
      ! <--
    end subroutine stress_exchange_part

    subroutine stress_correlation_part(i1,i2,in)
      integer, intent(in) :: i1,i2,in
      integer :: m
      real(kind=DP) :: strs
      strs = 0.d0
      do m = ista_fftph, iend_fftph         ! MPI
         strs=strs &
              &  - univol * rinplw * grad_trho(m) &
              &    * (cgrad_rho(m,in,1)+cgrad_rho(m,in,2)) * dF_drho(m,1) * f2or1(m)
      end do
      s_gga2(i1,i2) = s_gga2(i1,i2) + strs
    end subroutine stress_correlation_part
#else
    subroutine stress_exchange_part(i1,i2,is)
      integer, intent(in) :: i1,i2,is
      integer :: m

      do m = ista_fftph, iend_fftph       ! MPI
         s_gga1(i1,i2)=s_gga1(i1,i2) &
              &  + univol * rinplw * dF_dgradrho(m,is) &
              &    * cggawk13(m) * dF_drho(m,is) * f2or1(m)
      end do

      if(nspin.eq.1) then
         do m = ista_fftph, iend_fftph    ! MPI
            s_gga2(i1,i2)=s_gga2(i1,i2) &
                 &  + univol * rinplw * grad_trho(m) &
                 &    * cggawk13(m) * dF_drho(m,1) * f2or1(m)
         end do
      endif
    end subroutine stress_exchange_part

    subroutine stress_correlation_part(i1,i2)
      integer, intent(in) :: i1,i2
      integer :: m
      real(kind=DP) :: strs
      strs=0.d0
      do m = ista_fftph, iend_fftph         ! MPI
         strs=strs &
              &  + univol * rinplw * grad_trho(m) &
              &    * cggawk13(m) * dF_drho(m,1) * f2or1(m)
      end do

      s_gga2(i1,i2) = s_gga2(i1,i2) + strs
    end subroutine stress_correlation_part
#endif



! ================== added by K. Tagami ================ 11.0
    subroutine get_chgsoft_on_Rspace
      integer :: iloop

      do iloop = 1, ndim_magmom
         afft = 0.0d0
#ifdef _OLD_MAP_CHARGE_
         call map_charge_onto_a_fft_box( chgsoft, ndim_magmom, iloop, &
        &                                Valence_Charge_Only )
                                             !-(m_XC_Pot.) -> afft(*) (xcchg2)
#else
         call map_charge_onto_a_fft_box2( nfout, chgsoft, ndim_magmom, iloop, &
        &                                 Valence_Charge_Only )
                                              !-(m_XC_Pot.) -> afft(*) (xcchg2)
#endif
!         if(iprixc >= 2) write(nfout,*) " just after map_charge_onto_a_fft_box p1"
!         if(iprixc >= 2) call wd_small_part_of_afft(7,'G space',120)

         call m_FFT_CD_inverse_c(nfout,afft)  ! afft(G_sp.) -> afft(R_sp.)
! --
         do i = ista_fftph, iend_fftph    ! MPI
            chgsoft_on_Rspace(i,iloop) = afft(i*2-1)
         end do
      end do

    end subroutine get_chgsoft_on_Rspace

    subroutine alloc_chgsoft_project_noncol
      allocate( chgsoft_projected_on_Gspace(ista_kngp:iend_kngp,kimg,ndim_magmom) )
      chgsoft_projected_on_Gspace = 0.0d0
    end subroutine alloc_chgsoft_project_noncol

    subroutine dealloc_chgsoft_project_noncol
      deallocate( chgsoft_projected_on_Gspace )
    end subroutine dealloc_chgsoft_project_noncol

    subroutine set_chgsoft_projected_noncol
      integer :: iloop

      allocate( chgsoft_on_Rspace(ista_fftph:iend_fftph,ndim_magmom) );
      allocate( chgsoft_along_QuantzAxis(ista_fftph:iend_fftph,nspin) );
      chgsoft_on_Rspace  = 0.0d0
      chgsoft_along_QuantzAxis  = 0.0d0

      call get_chgsoft_on_Rspace
      call m_ES_Chgsoft_Along_QuantzAxis( chgsoft_on_Rspace, &
           &                              chgsoft_along_QuantzAxis, &
           &                              quantz_axis_inversion_flg_mesh )
      !
      deallocate( chgsoft_on_Rspace )
      !
      Do iloop=1, nspin
         afft = 0.0d0
         do i = ista_fftph, iend_fftph    ! MPI
            afft(2*i-1) = chgsoft_along_QuantzAxis(i,iloop)
         end do

         call m_FFT_CD_direct_c( nfout, afft )
!!         afft = afft * rinplw
         call cpafft_to_chgsoft_proj( iloop )

      End do
      deallocate( chgsoft_along_QuantzAxis )

    end subroutine set_chgsoft_projected_noncol

    subroutine cpafft_to_chgsoft_proj( iloop )
      integer, intent(in) :: iloop

      integer :: i, ik, ip
      real(kind=DP),allocatable,dimension(:)    :: afft_mpi1   ! MPI d(nfftp)
      integer :: id_sname = -1

      call tstatc0_begin('cpafft_to_chgsoft_proj ',id_sname)

      allocate(afft_mpi1(nfftp))

      call afft_allgatherv(afft,afft_mpi1)

      do ik = 1, kimg
         do i = ista_kngp, min(kgp_reduced,iend_kngp)  ! mpi
#ifdef _MPIFFT_
            ip = (igfp_l_c(i)-1)*kimg + ik
#else
            ip = (igfp_l(i)-1)*kimg + ik
#endif
            chgsoft_projected_on_Gspace(i,ik,iloop) = afft_mpi1(ip)*rinplw
         end do
      end do

      deallocate(afft_mpi1)
      call tstatc0_end(id_sname)
    end subroutine cpafft_to_chgsoft_proj

    subroutine alloc_chghard_project_noncol
      allocate( hsr_along_QuantzAxis(natm,nlmt,nlmt,nspin) )
      allocate( hsrd_along_QuantzAxis(natm,nlmt,nlmt,nspin,3,3) )
      hsr_along_QuantzAxis = 0.0d0
      hsrd_along_QuantzAxis = 0.0d0
    end subroutine alloc_chghard_project_noncol

    subroutine dealloc_chghard_project_noncol
      deallocate( hsr_along_QuantzAxis )
      deallocate( hsrd_along_QuantzAxis )
    end subroutine dealloc_chghard_project_noncol

    subroutine set_chghard_projected_noncol
      integer :: iloop
      integer :: ia
!
! -----
!      if ( .not. flg_paw ) then
!         call m_CD_estim_magmom_local(nfout)
!      endif
!
!      allocate( quantz_axis_inversion_flg_atm(natm) )
!
!      call m_ES_QuantzAxis_inv_flg_atm( magmom_local_now, &
!       &                                quantz_axis_inversion_flg_atm )
!      call m_ES_Chghard_Along_QuantzAxis( hsr, hsrd,  hsr_along_QuantzAxis, &
!       &                                  hsrd_along_QuantzAxis, &
!       &                                  quantz_axis_inversion_flg_atm )
! --------
      call m_ES_Chghard_Along_QuantzAxis( hsr, hsrd,  hsr_along_QuantzAxis, &
       &                                  hsrd_along_QuantzAxis )

!!      deallocate( quantz_axis_inversion_flg_atm )

    end subroutine set_chghard_projected_noncol

! ====================================================== 11.0

    subroutine rhos_diff(i1,i2)
      integer, intent(in) :: i1,i2
      integer :: is,ig
      integer :: iend  !mpi
      drhodh = 0.d0
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
    end subroutine rhos_diff


! ======================================== added by K. Tagami ============= 11.0
    subroutine rhos_diff_noncl(i1,i2)
      integer, intent(in) :: i1,i2
      integer :: is,ig
      integer :: iend  !mpi
      drhodh = 0.d0

      if ( sw_constrain_on_grad_correction == ON ) then
         do is = 1, nspin
            iend = iend_kngp

            if( iend_kngp > kg ) iend = kg
            if( ista_kngp <= iend ) then
               do ig = ista_kngp, iend  !for mpi
                  drhodh(ig,1,is) = -chgsoft_projected_on_Gspace(ig,1,is) &
                       &             *alinvt(i1,i2)
                  drhodh(ig,kimg,is) = -chgsoft_projected_on_Gspace(ig,kimg,is) &
                       &               *alinvt(i1,i2)
               enddo
            endif
         end do
      else
         do is = 1, ndim_magmom
            iend = iend_kngp

            if( iend_kngp > kg ) iend = kg
            if( ista_kngp <= iend ) then
               do ig = ista_kngp, iend  !for mpi
                  drhodh(ig,1,is) = - chgsoft(ig,1,is) * alinvt(i1,i2)
                  drhodh(ig,kimg,is) = - chgsoft(ig,kimg,is) * alinvt(i1,i2)
               enddo
            endif
         end do
      endif
    end subroutine rhos_diff_noncl
! ====================================================================== 11.0

    subroutine rhopc_diff(i1,i2)
      integer, intent(in) :: i1,i2
      integer :: ia,it,ipc,ig
      real(kind=DP) :: pc
      do ia=1,natm
      it=ityp(ia)
      if(itpcc(it)==0) cycle
      call calc_phase2(natm,pos,ia,kgp,ngabc,ista_kngp,iend_kngp,zfcos,zfsin)
                       ! -(b_Elec.)  -> zfcos, zfsin
      ipc=itpcc(it)
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
      enddo
    end subroutine rhopc_diff

! ========================================= added by K. Tagami ============ 11.0
    subroutine rhopc_diff_noncl(i1,i2)
      integer, intent(in) :: i1,i2
      integer :: ia,it,ipc,ig
      real(kind=DP) :: pc, ctmp

      if ( sw_constrain_on_grad_correction == ON ) then
         do ia=1, natm
            it = ityp(ia);      ipc = itpcc(it)
            if ( itpcc(it)==0 ) cycle

            call calc_phase2( natm, pos, ia, kgp, ngabc, ista_kngp, iend_kngp, &
                 &            zfcos, zfsin )
                                      ! -(b_Elec.)  -> zfcos, zfsin

            do ig = ista_kngp, iend_kngp  !for mpi
               ctmp = g(ig,1)*alinvt(1,i2) +g(ig,2)*alinvt(2,i2) +g(ig,3)*alinvt(3,i2)
               pc = rhpcg_l(ig,ipc)* alinvt(i1,i2) &
                    &   +rhpcg_diff_l(ig,ipc) *grinv(ig) *g(ig,i1) *ctmp
               pc = pc *dble(iwei(ia))
               drhodh(ig,1,1)        = drhodh(ig,1,1)        -zfcos(ig) *pc /nspin
               drhodh(ig,1,nspin)    = drhodh(ig,1,nspin)    -zfcos(ig) *pc /nspin
               drhodh(ig,kimg,1)     = drhodh(ig,kimg,1)     +zfsin(ig) *pc /nspin
               drhodh(ig,kimg,nspin) = drhodh(ig,kimg,nspin) +zfsin(ig) *pc /nspin
            end do
         end do

      else
         do ia=1, natm
            it = ityp(ia);   ipc = itpcc(it)
            if ( itpcc(it)==0 ) cycle

            call calc_phase2( natm, pos, ia, kgp, ngabc, ista_kngp, iend_kngp, &
                 &            zfcos, zfsin )
                                      ! -(b_Elec.)  -> zfcos, zfsin

            do ig = ista_kngp, iend_kngp  !for mpi
               ctmp = g(ig,1)*alinvt(1,i2) +g(ig,2)*alinvt(2,i2) +g(ig,3)*alinvt(3,i2)
               pc = rhpcg_l(ig,ipc)* alinvt(i1,i2) &
                    &   +rhpcg_diff_l(ig,ipc) *grinv(ig) *g(ig,i1) *ctmp
               pc = pc *dble(iwei(ia))
               drhodh(ig,1,1)    = drhodh(ig,1,1)    -zfcos(ig) *pc
               drhodh(ig,kimg,1) = drhodh(ig,kimg,1) +zfsin(ig) *pc
            enddo
         end do
      endif

    end subroutine rhopc_diff_noncl
! ======================================================================== 11.0

    subroutine rhoh_diff(i1,i2)
      integer, intent(in) :: i1,i2
      integer    :: n,ilm3,i,it,mdvdb,ia
      integer    :: lmt1,lmt2,il1,il2,tau1,tau2,l3,iiqitg
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

      if(iprixc >= 2) write(6,'(" << rhoh_diff >>")')
      call m_PP_find_maximum_l(n)   ! n-1: maximum l
      n = (n-1) + (n-1) + 1
      allocate(il3(n**2)); call substitute_il3(n**2,il3) ! -(b_Elec..)

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
!!$!mpi         call m_pwBS_sphrp2(ilm3,rltv,ylm2(1,ilm3))
!!$         call m_pwBS_sphrp2(ilm3,rltv,ylm2_mpi)
!!$         do i = ista_kngp, iend_kngp  !for mpi
!!$            ylm2(i,ilm3) = ylm2_mpi(i)
!!$         enddo
!mpi         call m_pwBS_sphrp2_diff(ilm3,rltv,ylmd(1,1,ilm3))
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
      allocate(mzfsin_blk(ibl1:ibl2))
      do ia = 1, natm
         it = ityp(ia)
         mdvdb = m_PP_include_vanderbilt_pot(it)
         if(mdvdb == SKIP) cycle
         call calc_phase2(natm,pos,ia,kgp,ngabc,ibl1,ibl2,zfcos_blk,zfsin_blk)
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
                  if ( noncol ) then
                     !if(mod(il1+il2,2) == 0) then
                     if(mod(l3,2) == 0) then
                        call even_case_noncl(ibl1,ibl2,i1,i2,ilm3,iq,l3,ia,lmt1,lmt2,fac,dga)
                     else
                        call odd_case_noncl(ibl1,ibl2,i1,i2,ilm3,iq,l3,ia,lmt1,lmt2,fac,dga)
                     endif
                  else
                     !if(mod(il1+il2,2) == 0) then
                     if(mod(l3,2) == 0) then
                        call even_case(ibl1,ibl2,i1,i2,ilm3,iq,l3,ia,lmt1,lmt2,fac,dga)
                     else
                        call odd_case(ibl1,ibl2,i1,i2,ilm3,iq,l3,ia,lmt1,lmt2,fac,dga)
                     endif
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
         call calc_phase2(natm,pos,ia,kgp,ngabc,ista_kngp,iend_kngp,zfcos,zfsin)
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
                   ylm => ylm_l(ista_kngp:iend_kngp,ilm3)
                else
                   call m_pwBS_sphrp2_3D(ilm3,rltv,ista_kngp,iend_kngp,ylm_t)
                   ylm => ylm_t(ista_kngp:iend_kngp)
                end if
                dga = dl2p(lmt1,lmt2,n,it)

! ===============================- modified by K. Tagami ================ 11.0
!                if(mod(il1+il2,2) == 0) then
!                   call even_case(i1,i2,ilm3,iiqitg,l3,ia,lmt1,lmt2,fac,dga)
!                else
!                   call odd_case(i1,i2,ilm3,iiqitg,l3,ia,lmt1,lmt2,fac,dga)
!                endif
                if ( noncol ) then
                   if(mod(il1+il2,2) == 0) then
                      call even_case_noncl(i1,i2,ilm3,iiqitg,l3,ia,lmt1,lmt2,fac,dga)
                   else
                      call odd_case_noncl(i1,i2,ilm3,iiqitg,l3,ia,lmt1,lmt2,fac,dga)
                   endif
                else
                   if(mod(il1+il2,2) == 0) then
                      call even_case(ibl1,ibl2,i1,i2,ilm3,iiqitg,l3,ia,lmt1,lmt2,fac,dga)
                   else
                      call odd_case(ibl1,ibl2,i1,i2,ilm3,iiqitg,l3,ia,lmt1,lmt2,fac,dga)
                   endif
                endif
! ========================================================================== 11.0

              enddo
            enddo
         enddo
      enddo
#endif
      enddo
#ifndef _m_XC_no_loop_exchange_
      deallocate(nc2n,nc2lmt2,nc2lmt1,nc,iq2l3,nqitg_sp,nqitg_sp0)
#endif
!      deallocate(il3);deallocate(ylm_t);deallocate(ylmd)
      deallocate(il3);deallocate(ylmd)
      if(allocated(ylm_ext)) deallocate(ylm_ext)

    end subroutine rhoh_diff

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

    subroutine even_case(ibl1,ibl2,i1,i2,ilm3,iiqitg,l3,ia,lmt1,lmt2,fac,dga)
      integer, intent(in) :: ibl1,ibl2,i1,i2,ilm3,iiqitg,l3,ia,lmt1,lmt2
      real(kind=DP), intent(in) :: fac,dga
      integer :: is

      do is=1,nspin,af+1
        flchgq(is) = fac*real(zi**(-l3))*dga &
                   & *hsr(ia,lmt1,lmt2,is)
        flchgqd(i1,i2,is) = fac*real(zi**(-l3))*dga &
                          & *hsrd(ia,lmt1,lmt2,is,i1,i2)
      enddo
      if(kimg == 1) then
         call real_case(ibl1,ibl2,i1,i2,ilm3,iiqitg,ia,zfcos_blk)
      else
         call complex_case(ibl1,ibl2,i1,i2,ilm3,iiqitg,zfcos_blk,mzfsin_blk)
      endif
    end subroutine even_case

! =================================== added by K. Tagami ================= 11.0
    subroutine even_case_noncl( ibl1, ibl2, i1, i2, ilm3, iiqitg, l3, ia, lmt1, lmt2, fac, dga )
      integer, intent(in) :: ibl1,ibl2,i1,i2,ilm3,iiqitg,l3,ia,lmt1,lmt2
      real(kind=DP), intent(in) :: fac,dga
      integer :: is

      if ( sw_constrain_on_grad_correction == ON ) then
         do is=1,nspin,af+1
            flchgq(is) = fac *real(zi**(-l3)) *dga &
                 &           *hsr_along_QuantzAxis(ia,lmt1,lmt2,is)
            flchgqd(i1,i2,is) = fac *real(zi**(-l3)) *dga &
                 &                  *hsrd_along_QuantzAxis(ia,lmt1,lmt2,is,i1,i2)
         enddo
         call complex_case(ibl1,ibl2,i1,i2,ilm3,iiqitg,zfcos_blk,mzfsin_blk)
      else
         do is=1,ndim_magmom
            flchgq(is) = fac *real(zi**(-l3)) *dga &
                 &           *hsr(ia,lmt1,lmt2,is)
            flchgqd(i1,i2,is) = fac *real(zi**(-l3)) *dga &
                 &                  *hsrd(ia,lmt1,lmt2,is,i1,i2)
         enddo
         call complex_case_noncl(ibl1,ibl2,i1,i2,ilm3,iiqitg,zfcos_blk,mzfsin_blk)
      endif

    end subroutine even_case_noncl
! ================================================================= 11.0

    subroutine odd_case(ibl1,ibl2,i1,i2,ilm3,iiqitg,l3,ia,lmt1,lmt2,fac,dga)
      integer, intent(in) :: ibl1,ibl2,i1,i2,ilm3,iiqitg,l3,ia,lmt1,lmt2
      real(kind=DP), intent(in) :: fac,dga
      integer :: is

      do is=1,nspin,af+1
        flchgq(is) = fac*aimag(zi**(-l3))*dga &
                   & *hsr(ia,lmt1,lmt2,is)
        flchgqd(i1,i2,is) = fac*aimag(zi**(-l3))*dga &
                          & *hsrd(ia,lmt1,lmt2,is,i1,i2)
      enddo
      if(kimg == 1) then
         call real_case(ibl1,ibl2,i1,i2,ilm3,iiqitg,ia,zfsin_blk)
      else
         call complex_case(ibl1,ibl2,i1,i2,ilm3,iiqitg,zfsin_blk,zfcos_blk)
      endif
    end subroutine odd_case

! =================================== added by K. Tagami ================= 11.0
    subroutine odd_case_noncl(ibl1, ibl2, i1, i2, ilm3, iiqitg, l3, ia, lmt1, lmt2, fac, dga )
      integer, intent(in) :: ibl1, ibl2, i1,i2,ilm3,iiqitg,l3,ia,lmt1,lmt2
      real(kind=DP), intent(in) :: fac,dga
      integer :: is

      if ( sw_constrain_on_grad_correction == ON ) then
         do is=1, nspin, af+1
            flchgq(is) = fac *aimag(zi**(-l3)) *dga &
                 &           *hsr_along_QuantzAxis(ia,lmt1,lmt2,is)
            flchgqd(i1,i2,is) = fac *aimag(zi**(-l3)) *dga &
                 &                  *hsrd_along_QuantzAxis(ia,lmt1,lmt2,is,i1,i2)
         enddo
         call complex_case(ibl1,ibl2,i1,i2,ilm3,iiqitg,zfsin_blk,zfcos_blk)
      else
         do is=1,ndim_magmom
            flchgq(is) = fac *aimag(zi**(-l3)) *dga &
                 &           *hsr(ia,lmt1,lmt2,is)
            flchgqd(i1,i2,is) = fac *aimag(zi**(-l3)) *dga &
                 &                  *hsrd(ia,lmt1,lmt2,is,i1,i2)
         enddo
         call complex_case_noncl(ibl1,ibl2,i1,i2,ilm3,iiqitg,zfsin_blk,zfcos_blk)
      endif

    end subroutine odd_case_noncl
! ================================================================= 11.0

    subroutine real_case(ibl1,ibl2,i1,i2,ilm3,iiqitg,ia,zf)
      integer,       intent(in) :: ibl1,ibl2,i1,i2,ilm3,iiqitg,ia
      real(kind=DP), intent(in), dimension(ibl1:ibl2) :: zf
      integer :: is, i, iy
      real(kind=DP) :: c1,c2,c3up,c3dn

      do is=1,nspin,af+1
        flchgq(is) = flchgq(is) * dble(iwei(ia))
        flchgqd(i1,i2,is) = flchgqd(i1,i2,is) * dble(iwei(ia))
      enddo
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
    end subroutine real_case

    subroutine complex_case(is,ie,i1,i2,ilm3,iiqitg,zf1,zf2)
      integer,       intent(in) :: is,ie,i1,i2,ilm3,iiqitg
!      real(kind=DP), intent(in), dimension(ista_kngp:iend_kngp) &
!           &        :: zf1,zf2      ! MPI
      real(kind=DP), intent(in), dimension(is:ie) :: zf1,zf2
      real(kind=DP) :: c1,c2,c3up,c3dn
      integer       :: i, iy
      integer :: id_sname = -1

      if(nspin==1 .or. af==1) then
!         do i = ista_kngp, iend_kngp  !for mpi
         do i = is, ie !for mpi
!            iy = i - ista_kngp+1
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
!         do i = ista_kngp, iend_kngp  !for mpi
         do i = is, ie  !for mpi
!            iy = i - ista_kngp+1
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
    end subroutine complex_case

! ======================================= added by K. Tagami ================= 11.0
    subroutine complex_case_noncl(ibl1,ibl2,i1,i2,ilm3,iiqitg,zf1,zf2)
      integer,       intent(in) :: ibl1,ibl2,i1,i2,ilm3,iiqitg
      real(kind=DP), intent(in), dimension(ibl1:ibl2) &
           &        :: zf1,zf2      ! MPI
      real(kind=DP) :: c1,c2,c3
      integer       :: i, is, iy

      do i = ibl1, ibl2  !for mpi
         iy = i - ibl1+1
         c1=(  g(i,1)*alinvt(1,i2) &
              &  + g(i,2)*alinvt(2,i2) &
              &  + g(i,3)*alinvt(3,i2) )  *g(i,i1) * grinv(i)
         c2=(  ylmd(i,1,ilm3)*alinvt(1,i2) &
              &  + ylmd(i,2,ilm3)*alinvt(2,i2) &
              &  + ylmd(i,3,ilm3)*alinvt(3,i2) )  *g(i,i1)

         Do is=1, ndim_magmom
            c3 =  -flchgq(is) &
                 &  *(ylm(iy)*(qitg_l(i,iiqitg)*alinvt(i1,i2) &
                 &  +qitg_diff_l(i,iiqitg)*c1) &
                 &  +c2*qitg_l(i,iiqitg)) &
                 &  + flchgqd(i1,i2,is)*ylm(iy)*qitg_l(i,iiqitg)
            drhodh(i,1,is)    = drhodh(i,1,is)    +c3 *zf1(i)
            drhodh(i,kimg,is) = drhodh(i,kimg,is) +c3 *zf2(i)
         End do
      enddo

    end subroutine complex_case_noncl
! ==================================================================== 11.0

    subroutine map_drhodh1(is)
      integer, intent(in) :: is
      integer :: i,ip,ri,j
      real(kind=DP), allocatable, dimension(:) :: chgfft_mpi4
      real(kind=DP),allocatable,dimension(:)    :: afft_mpi2,afft_mpi3,afft_mpi4   ! MPI d(mp_fftp)

      chgfft_mpi = 0.d0
      do ri = 1, kimg
         do i = ista_kngp, min(kgp_reduced, iend_kngp)  !for mpi
#ifdef _MPIFFT_
            ip = igfp_l_c(i)*kimg - kimg + ri
#else
            ip = igfp_l(i)*kimg - kimg + ri
#endif
            chgfft_mpi(ip)  = chgfft_mpi(ip) + drhodh(i,ri,is)
         end do
      end do

      if(npes_cdfft > 1) then
         call mpi_barrier(MPI_CommGroup,ierr)
         allocate(afft_mpi2(mp_fftp))
         allocate(afft_mpi3(mp_fftp))
         do j = 0, npes_cdfft-1
#ifdef _MPIFFT_
            call set_cdata(j,chgfft_mpi,afft_mpi2)
#else
            do i = nis_fftp(j), nie_fftp(j)
               afft_mpi2(i-nis_fftp(j)+1) = chgfft_mpi(i)
            end do
#endif
            call mpi_allreduce(afft_mpi2,afft_mpi3,mp_fftp &
                 &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)

            if(j == myrank_cdfft) then
               do i = ista_fftp, iend_fftp
                  chgfft(i) = afft_mpi3(i - ista_fftp + 1)
               end do
            end if
         end do
         deallocate(afft_mpi3,afft_mpi2)
      else
         if(nrank_ggacmp > 1) then
            allocate(afft_mpi4(nfftp))
            call mpi_allreduce(chgfft_mpi, afft_mpi4, nfftp &
                 & , mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
            chgfft = afft_mpi4
            deallocate(afft_mpi4)
         else
            chgfft = chgfft_mpi
         end if
      end if
!!$      call subst_A_into_A_reduced  ! -(contained in subrt. m_XC_cal_potential) MPI
!!$      !   A=chgfft_mpi,A_reduced=chgfft, using work arrays; chgfft_mpi -> chgfft
    end subroutine map_drhodh1

!!$    subroutine subst_A_into_A_reduced
!!$      real(kind=DP), pointer, dimension(:)      :: chgfft_mpi_reduced &
!!$           &          ,                            chgfft_mpi_reduced2   ! MPI
!!$      integer :: i, ip
!!$      if(npes_cdfft >= 2) then
!!$         allocate(chgfft_mpi_reduced(mp_fftp))   ! MPI
!!$         allocate(chgfft_mpi_reduced2(mp_fftp))  ! MPI
!!$         do i = 0, npes_cdfft-1                                            ! MPI
!!$            do ip = 1, nel_fftp(i)                                        ! MPI
!!$               chgfft_mpi_reduced(ip) = chgfft_mpi(ip + nis_fftp(i) - 1)  ! MPI
!!$            end do! MPI
!!$            if(myrank_ggacmp < nrank_ggacmp) then
!!$               call mpi_allreduce(chgfft_mpi_reduced,chgfft_mpi_reduced2, nel_fftp(i) &
!!$                    &, mpi_double_precision, mpi_sum, mpi_cdfft_world(myrank_ggacmp), ierr)  ! MPI
!!$            end if
!!$
!!$            if(myrank_cdfft == i) then                                            ! MPI
!!$               do ip = 1, nel_fftp(i)                                     ! MPI
!!$                  chgfft(ip + ista_fftp - 1) = chgfft_mpi_reduced2(ip)    ! MPI
!!$               end do                                                     ! MPI
!!$            end if                                                        ! MPI
!!$         end do                                                           ! MPI
!!$         deallocate(chgfft_mpi_reduced)
!!$         deallocate(chgfft_mpi_reduced2)
!!$      else if(npes_cdfft == 1) then
!!$         chgfft = chgfft_mpi
!!$      end if
!!$    end subroutine subst_A_into_A_reduced

    subroutine map_drhodh2
      integer :: i,ip,ri,j
      real(kind=DP), allocatable, dimension(:) :: chgfft_mpi4
      real(kind=DP),allocatable,dimension(:)    :: afft_mpi2,afft_mpi3,afft_mpi4   ! MPI d(mp_fftp)

      chgfft_mpi = 0.d0
      do ri = 1, kimg
         do i = ista_kngp, min(kgp_reduced, iend_kngp)  ! for mpi
#ifdef _MPIFFT_
            ip = igfp_l_c(i)*kimg - kimg + ri
#else
            ip = igfp_l(i)*kimg - kimg + ri
#endif
            chgfft_mpi(ip) = chgfft_mpi(ip) + drhodh(i,ri,1) + drhodh(i,ri,nspin)
         end do
      end do

      if(npes_cdfft > 1) then
         call mpi_barrier(MPI_CommGroup,ierr)
         allocate(afft_mpi2(mp_fftp))
         allocate(afft_mpi3(mp_fftp))
         do j = 0, npes_cdfft-1
#ifdef _MPIFFT_
            call set_cdata(j,chgfft_mpi,afft_mpi2)
#else
            do i = nis_fftp(j), nie_fftp(j)
               afft_mpi2(i-nis_fftp(j)+1) = chgfft_mpi(i)
            end do
#endif
            call mpi_allreduce(afft_mpi2,afft_mpi3,mp_fftp &
                 &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)

            if(j == myrank_cdfft) then
               do i = ista_fftp, iend_fftp
                  chgfft(i) = afft_mpi3(i - ista_fftp + 1)
               end do
            end if
         end do
         deallocate(afft_mpi3,afft_mpi2)
      else
         if(nrank_ggacmp > 1) then
            allocate(afft_mpi4(nfftp))
            call mpi_allreduce(chgfft_mpi, afft_mpi4, nfftp &
                 & , mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
            chgfft = afft_mpi4
            deallocate(afft_mpi4)
         else
            chgfft = chgfft_mpi
         end if
      end if

!!$      call subst_A_into_A_reduced   ! -(contained in subrt. m_XC_cal_potential) MPI
!!$      !   A=chgfft_mpi,A_reduced=chgfft, using work arrays; chgfft_mpi -> chgfft
    end subroutine map_drhodh2

    subroutine dgrhodh1(i1,i2,in,is)
      integer, intent(in) :: i1,i2,in,is
      integer       :: n
      real(kind=DP) :: g1,g2

      do n = ista_fftp, iend_fftp                     ! MPI
         g1 = rltv(in,1)*inx(n)+rltv(in,2)*jnx(n)+rltv(in,3)*knx(n)
         g2 = rltv(i1,1)*inx(n)+rltv(i1,2)*jnx(n)+rltv(i1,3)*knx(n)
         afft(n) = g1*chgfft(n) - chden_l(n,is)*g2*alinvt(in,i2)
      enddo
      call boundary_zero_into_afft(in)  ! -(contained in subr. m_XC_cal_potential)

      call m_FFT_CD_inverse_c(nfout,afft)

      do n = ista_fftph, iend_fftph     ! MPI
         dF_drho(n,is) = afft(2*n)
      enddo
    end subroutine dgrhodh1

! ======================================== added by K. Tagami ================ 11.0
    subroutine dgrhodh1_noncl(i1,i2,in,is)
      integer, intent(in) :: i1,i2,in,is
      integer       :: n
      real(kind=DP) :: g1,g2

      do n = ista_fftp, iend_fftp                     ! MPI
         g1 = rltv(in,1)*inx(n)+rltv(in,2)*jnx(n)+rltv(in,3)*knx(n)
         g2 = rltv(i1,1)*inx(n)+rltv(i1,2)*jnx(n)+rltv(i1,3)*knx(n)
         afft(n) = g1*chgfft(n) - chden_l(n,is)*g2*alinvt(in,i2)
      enddo
      call boundary_zero_into_afft(in)  ! -(contained in subr. m_XC_cal_potential)

      call m_FFT_CD_inverse_c(nfout,afft)

    end subroutine dgrhodh1_noncl
! ===================================================================== 11.0

    subroutine dgrhodh2(i1,i2,in)
      integer, intent(in) :: i1,i2,in
      integer       :: n
      real(kind=DP) :: g1,g2

      do n = ista_fftp, iend_fftp
         g1 = rltv(in,1)*inx(n)+rltv(in,2)*jnx(n)+rltv(in,3)*knx(n)
         g2 = rltv(i1,1)*inx(n)+rltv(i1,2)*jnx(n)+rltv(i1,3)*knx(n)
         afft(n) = g1*chgfft(n) &
              &   - (chden_l(n,1)+chden_l(n,nspin))*g2*alinvt(in,i2)
      enddo
      call boundary_zero_into_afft(in)   ! -(contained in subr. m_XC_cal_potential)

      call m_FFT_CD_inverse_c(nfout,afft)

      do n = ista_fftph, iend_fftph      ! MPI
         dF_drho(n,1) = afft(2*n)
      enddo
    end subroutine dgrhodh2

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

    subroutine allocate_for_stress
      integer :: i
      real(kind=DP) :: ga,gb,gc
      allocate(alinvt(3,3)); alinvt = rltv / PAI2
      allocate(zfcos(ista_kngp:iend_kngp)); zfcos = 0.d0
      allocate(zfsin(ista_kngp:iend_kngp)); zfsin = 0.d0

! ========== KT_mod ============ 2013/10/29
!      allocate(drhodh(ista_kngp:iend_kngp,kimg,nspin)); drhodh = 0.d0
!      allocate(flchgq(nspin)); flchgq = 0.d0
!      allocate(flchgqd(3,3,nspin)); flchgqd = 0.d0
!
      if ( noncol .and. sw_constrain_on_grad_correction == OFF ) then
         allocate(drhodh(ista_kngp:iend_kngp,kimg,ndim_magmom)); drhodh = 0.d0
         allocate(flchgq(ndim_magmom)); flchgq = 0.d0
         allocate(flchgqd(3,3,ndim_magmom)); flchgqd = 0.d0
      else
         allocate(drhodh(ista_kngp:iend_kngp,kimg,nspin)); drhodh = 0.d0
         allocate(flchgq(nspin)); flchgq = 0.d0
         allocate(flchgqd(3,3,nspin)); flchgqd = 0.d0
      endif
! ============================= 2013/10/29

      allocate(chgfft(ista_fftp:iend_fftp)); chgfft = 0.d0  ! MPI
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
        ga = ngabc(i,1)
        gb = ngabc(i,2)
        gc = ngabc(i,3)
        g(i,1) = rltv(1,1)*ga+rltv(1,2)*gb+rltv(1,3)*gc
        g(i,2) = rltv(2,1)*ga+rltv(2,2)*gb+rltv(2,3)*gc
        g(i,3) = rltv(3,1)*ga+rltv(3,2)*gb+rltv(3,3)*gc
      enddo
      allocate(chgfft_mpi(nfftp))             ! MPI

    end subroutine allocate_for_stress

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
    end subroutine deallocate_for_stress

    subroutine abs_grad_rho_up_down_total()
      real(kind=DP),allocatable,dimension(:)    :: grad_rho_c2  ! MPI d(ista_fftph:iend_fftph)

      integer  :: is, in, i
      real(kind=DP) :: x,y,z

#ifdef _DETAIL_GGA_TIMING_
      integer :: id_sname = -1
#endif

      if(nrank_ggacmp > 1) allocate(grad_rho_c2(ista_fftph:iend_fftph))    ! MPI

      do is = 1, ispin
         grad_rho(:,is) = 0.d0
         do in = 1, 3
            if(map_ggacmp(in) /= myrank_ggacmp) cycle
            call g_xyz_chden_l(in,is)        ! G_xyz * rho_{up|down}(G)-->afft
            if(iprixc >= 2) call wd_small_part_of_afft(34,'G space before<m_FFT_CD_inverse_c)',120)
            call m_FFT_CD_inverse_c(nfout,afft)! (-i)*d(rho_{up|down}(r))/d(x|y|z)
            if(iprixc >= 2) call wd_small_part_of_afft(34,'R space after <m_FFT_CD_inverse_c)',120)
#ifndef _XC_SAVE_MEMORY_
            call cp_afft_to_cgrad_rho(is,in) ! -> cgrad_rho(i,in,is) = -afft(i*2)
#endif
            call add_sq_afft_to_grad_rho(is) ! grad_rho <--  + afft**2
            if(iprixc >= 2) then
               write(nfout,'(" !XC after add_sq_afft_to_grad_rho: is, in = ",2i8)') is, in
               write(nfout,'(" !XC -- grad_rho -- <<abs_grad_rho_up_down_total>>")')
               write(nfout,'(" !XC ",6f12.6)') (grad_rho(i,1),i=ista_fftph,ista_fftph+17)
               write(nfout,'(" !XC ",6f12.6)') (grad_rho(i,1),i=iend_fftph-17,iend_fftph)
#ifndef _XC_SAVE_MEMORY_
               write(nfout,'(" !XC -- cgrad_rho -- <<abs_grad_rho_up_down_total>>")')
               write(nfout,'(" !XC ",6f12.6)') (cgrad_rho(i,in,1),i=ista_fftph,ista_fftph+17)
               write(nfout,'(" !XC ",6f12.6)') (cgrad_rho(i,in,1),i=iend_fftph-17,iend_fftph)
#endif
            end if

         end do

         if(iprixc >= 2) write(nfout,'(" !XC after add_sq_afft_to_grad_rho")')

         ! --> grad_rho_c
         if(nrank_ggacmp > 1) then
            grad_rho_c2(:) = grad_rho(:,is)
#ifdef _DETAIL_GGA_TIMING_
            call tstatc0_begin('mpi_allreduce(ggaxcp) ',id_sname)
#endif
            call mpi_allreduce(grad_rho_c2,grad_trho,iend_fftph-ista_fftph+1 &
                 & ,mpi_double_precision, mpi_sum, mpi_ggacmp_cross_world(myrank_cdfft),ierr)
#ifdef _DETAIL_GGA_TIMING_
            call tstatc0_end(id_sname)
#endif
            grad_rho(:,is) = dsqrt(grad_trho(:))
         else
            grad_rho(:,is) = dsqrt(grad_rho(:,is)) ! grad_rho(s) <-- |grad(rho_s(r))|
         end if
      enddo
      if(ispin == 2) then
#ifndef _XC_SAVE_MEMORY_
         if(nrank_ggacmp > 1) then
            grad_rho_c2 = 0.d0
!!$            if(myrank_ggacmp >= 0 .and. myrank_ggacmp <=2) then
            do i = ista_fftph, iend_fftph
               x = cgrad_rho(i,myrank_ggacmp+1,1) + cgrad_rho(i,myrank_ggacmp+1,2)
               grad_rho_c2(i) = grad_rho_c2(i) + x*x
            end do
#ifdef _DETAIL_GGA_TIMING_
            call tstatc0_begin('mpi_allreduce(ggaxcp) ',id_sname)
#endif
            call mpi_allreduce(grad_rho_c2,grad_trho,iend_fftph-ista_fftph+1 &
                 & , mpi_double_precision, mpi_sum, mpi_ggacmp_cross_world(myrank_cdfft),ierr)
#ifdef _DETAIL_GGA_TIMING_
            call tstatc0_end(id_sname)
#endif
            grad_trho = dsqrt(grad_trho)
         else
            do in = ista_fftph, iend_fftph
               x = cgrad_rho(in,1,1) + cgrad_rho(in,1,2)
               y = cgrad_rho(in,2,1) + cgrad_rho(in,2,2)
               z = cgrad_rho(in,3,1) + cgrad_rho(in,3,2)
               grad_trho(in) = dsqrt(x*x+y*y+z*z)
            end do
         end if
#else
         do in = 1, 3
            call g_xyz_total_chden_l(in)     ! G_xyz*(rho(G)up+rho(G)down) -> afft
            call m_FFT_CD_inverse_c(nfout,afft)!(-i)*d(rho_total(r))/d(x|y|z)
            call add_sq_afft_to_grad_trho
         end do
         grad_trho = dsqrt(grad_trho)     ! grad_trho <-- |grad(rho_total(r))|
#endif
      else
         grad_trho(:) = grad_rho(:,1)
      end if

! ========================== added by K. Tagami ============ experimental === 11.0
#ifdef APPROX_GGA_NONCL
      if ( noncol ) then
         grad_trho(:) = grad_rho(:,1) + grad_rho(:,2)
      endif
#endif
! ============================================================================ 11.0

      if(nrank_ggacmp > 1) deallocate(grad_rho_c2)
    end subroutine abs_grad_rho_up_down_total

! ======================================== added by K. Tagami =============== 11.0
    subroutine abs_grad_rho_updown_total_noncl()
      real(kind=DP),allocatable,dimension(:)    :: grad_rho_c2
                                                    ! MPI d(ista_fftph:iend_fftph)

      integer  :: is, in, i
      real(kind=DP) :: x,y,z

#ifdef _DETAIL_GGA_TIMING_
      integer :: id_sname = -1
#endif

! ==== KT_add === 2014/08/04
      if ( sw_monitor_magnetic_vorticity == ON ) then
         allocate( MagVorticity(ista_fftph:iend_fftph, 3 ) )
         MagVorticity = 0.0d0
      endif
! =============== 2014/08/04

      allocate( bfft_kt( ista_fftp:iend_fftp,ndim_magmom )); bfft_kt = 0.0d0
      if (nrank_ggacmp > 1) allocate( grad_rho_c2(ista_fftph:iend_fftph) )    ! MPI

      grad_rho = 0.0d0;  grad_trho = 0.0d0

      do in = 1, 3
         if(map_ggacmp(in) /= myrank_ggacmp) cycle

         bfft_kt = 0.0d0

         Do is=1, ndim_magmom
            call g_xyz_chden_l(in,is)        ! G_xyz * rho_{up|down}(G)-->afft
            if ( iprixc >= 2 ) then
               call wd_small_part_of_afft(34,'G space before<m_FFT_CD_inverse_c)',120)
            endif
            call m_FFT_CD_inverse_c(nfout,afft)! (-i)*d(rho_{up|down}(r))/d(x|y|z)

            if (iprixc >= 2) then
               call wd_small_part_of_afft(34,'R space after <m_FFT_CD_inverse_c)',120)
            endif
            bfft_kt(:,is) = afft(:)
         End Do

! ====== KT_add ==== 2014/08/04
         if ( sw_monitor_magnetic_vorticity == ON ) then
            call m_ES_add_contrib_to_Vorticity( in, bfft_kt, MagVorticity )
         endif
! ================== 2014/08/04

         call m_ES_GradCorr_Along_QuantzAxis2( RhoMag_R, bfft_kt, &
              &                                quantz_axis_inversion_flg_mesh )

         afft = 0.0d0

         do is = 1, ispin
            afft(:) = bfft_kt(:,is)
#ifndef _XC_SAVE_MEMORY_
            call cp_afft_to_cgrad_rho(is,in) ! -> cgrad_rho(i,in,is) = -afft(i*2)
#endif
            call add_sq_afft_to_grad_rho(is) ! grad_rho <--  + afft**2

            if(iprixc >= 2) then
               write(nfout,'(" !XC after add_sq_afft_to_grad_rho: is, in = ",2i8)') is,in
               write(nfout,'(" !XC -- grad_rho -- <<abs_grad_rho_up_down_total>>")')
               write(nfout,'(" !XC ",6f12.6)') (grad_rho(i,1),i=ista_fftph,ista_fftph+17)
               write(nfout,'(" !XC ",6f12.6)') (grad_rho(i,1),i=iend_fftph-17,iend_fftph)

#ifndef _XC_SAVE_MEMORY_
               write(nfout,'(" !XC -- cgrad_rho -- <<abs_grad_rho_up_down_total>>")')
               write(nfout,'(" !XC ",6f12.6)') &
                    &               (cgrad_rho(i,in,1),i=ista_fftph,ista_fftph+17)
               write(nfout,'(" !XC ",6f12.6)') &
                    &               (cgrad_rho(i,in,1),i=iend_fftph-17,iend_fftph)
#endif
            end if

         end do

      end do

      if(iprixc >= 2) write(nfout,'(" !XC after add_sq_afft_to_grad_rho")')

      Do is=1, ispin
! --> grad_rho_c
         if(nrank_ggacmp > 1) then
            grad_rho_c2(:) = grad_rho(:,is)
#ifdef _DETAIL_GGA_TIMING_
            call tstatc0_begin('mpi_allreduce(ggaxcp) ',id_sname)
#endif
            call mpi_allreduce( grad_rho_c2, grad_trho, iend_fftph -ista_fftph +1, &
                 &              mpi_double_precision, mpi_sum, &
                 &              mpi_ggacmp_cross_world(myrank_cdfft), ierr )
#ifdef _DETAIL_GGA_TIMING_
            call tstatc0_end(id_sname)
#endif
            grad_rho(:,is) = dsqrt(grad_trho(:))
         else
            grad_rho(:,is) = dsqrt(grad_rho(:,is)) ! grad_rho(s) <-- |grad(rho_s(r))|
         end if
      enddo

      if (ispin == 2) then
#ifndef _XC_SAVE_MEMORY_
         if (nrank_ggacmp > 1) then
            grad_rho_c2 = 0.d0
!!$            if(myrank_ggacmp >= 0 .and. myrank_ggacmp <=2) then
            do i = ista_fftph, iend_fftph
               x = cgrad_rho(i,myrank_ggacmp+1,1) + cgrad_rho(i,myrank_ggacmp+1,2)
               grad_rho_c2(i) = grad_rho_c2(i) + x*x
            end do
#ifdef _DETAIL_GGA_TIMING_
            call tstatc0_begin('mpi_allreduce(ggaxcp) ',id_sname)
#endif
            call mpi_allreduce( grad_rho_c2, grad_trho, iend_fftph -ista_fftph +1, &
                 &              mpi_double_precision, mpi_sum, &
                 &              mpi_ggacmp_cross_world(myrank_cdfft), ierr )
#ifdef _DETAIL_GGA_TIMING_
            call tstatc0_end(id_sname)
#endif
            grad_trho = dsqrt(grad_trho)

         else
            do in = ista_fftph, iend_fftph
               x = cgrad_rho(in,1,1) + cgrad_rho(in,1,2)
               y = cgrad_rho(in,2,1) + cgrad_rho(in,2,2)
               z = cgrad_rho(in,3,1) + cgrad_rho(in,3,2)
               grad_trho(in) = dsqrt(x*x+y*y+z*z)
            end do
         end if
#else
         do in = 1, 3
!            call g_xyz_total_chden_l(in)     ! G_xyz*(rho(G)up+rho(G)down) -> afft
            call g_xyz_chden_l(in,1)     ! G_xyz*(rho(G)up+rho(G)down) -> afft
            call m_FFT_CD_inverse_c(nfout,afft)!(-i)*d(rho_total(r))/d(x|y|z)
            call add_sq_afft_to_grad_trho
         end do
         grad_trho = dsqrt(grad_trho)     ! grad_trho <-- |grad(rho_total(r))|
#endif
      else
         grad_trho(:) = grad_rho(:,1)
      end if

#ifdef APPROX_GGA_NONCL
      if ( noncol ) then
         grad_trho(:) = grad_rho(:,1) + grad_rho(:,2)
      endif
#endif

      if (nrank_ggacmp > 1) deallocate(grad_rho_c2)
      deallocate( bfft_kt )

    end subroutine abs_grad_rho_updown_total_noncl
! ============================================================================= 11.0

! ==== KT_add ================================= 13.0XX
    subroutine grad2_rho_up_down()
      integer  :: is, i

      real(kind=DP) :: x,y,z
      real(kind=DP), allocatable :: grad_rho_c2(:)   ! d(ista_fftph:iend_fftph)

      if (nrank_ggacmp > 1) allocate(grad_rho_c2(ista_fftph:iend_fftph))    ! MPI

      do is = 1, ispin
         grad2_rho(:,is) = 0.d0

         call g2_xyz_chden_l(is)        ! G2 * rho_{up|down}(G)-->afft
         call m_FFT_CD_inverse_c(nfout,afft)   !  Laplacian of (rho_up/dn (r) )
         call cp_afft_to_grad2_rho(is)

         ! --> grad_rho_c
         if (nrank_ggacmp > 1) then
            grad_rho_c2(:) = grad2_rho(:,is)
            call mpi_allreduce( grad_rho_c2, grad2_rho(:,is), iend_fftph-ista_fftph+1, &
                 &              mpi_double_precision, mpi_sum, &
                 &              mpi_ggacmp_cross_world(myrank_cdfft), ierr )
         end if

      enddo
      if (nrank_ggacmp > 1) deallocate(grad_rho_c2)

    end subroutine grad2_rho_up_down
! =============================================== 13.0XX

!!$    subroutine c2g2r_sfftph(grad_rho_c2,grad_rho)
!!$      real(kind=DP),intent(in), dimension(ista_fftph:iend_fftph) :: grad_rho_c2
!!$      real(kind=DP),intent(out),dimension(ista_sfftph: iend_sfftph)  :: grad_rho
!!$      real(kind=DP), allocatable, dimension(:) :: grad_rho_g ! d(nfftp)
!!$
!!$      if(npes_cdfft >= 2) then
!!$         allocate(grad_rho_g(nfftp))
!!$         if(myrank_ggacmp < nrank_ggacmp ) then
!!$            call mpi_allgatherv(grad_rho_c2,nel_fftph(myrank_cdfft),mpi_double_precision &
!!$                 & , grad_rho_g, nel_fftph, idisp_fftph &
!!$                 & , mpi_double_precision,mpi_cdfft_world(myrank_ggacmp),ierr)
!!$         end if
!!$         if(nrest_cdfft >= 1) then
!!$            if(mype > npes_cdfft*nrank_ggacmp - 1) then
!!$               call mpi_recv(grad_rho_g,nfftp,mpi_double_precision &
!!$                    & ,mype - npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,istatus,ierr)
!!$            end if
!!$            if(mype < nrest_cdfft) then
!!$               call mpi_send(grad_rho_g,nfftp,mpi_double_precision &
!!$                    & ,mype + npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,ierr)
!!$            end if
!!$         end if
!!$               ! -- scatter --
!!$         grad_rho(ista_sfftph:iend_sfftph) = grad_rho_g(ista_sfftph:iend_sfftph)
!!$         deallocate(grad_rho_g)
!!$      else
!!$         grad_rho(ista_sfftph:iend_sfftph) = grad_rho_c2(ista_sfftph:iend_sfftph)
!!$      end if
!!$    end subroutine c2g2r_sfftph

    subroutine dFxc_over_ddgradrho()
      real(kind=DP),allocatable,dimension(:)    :: grad_rho_c2  ! MPI d(ista_fftph:iend_fftph)
      real(kind=DP),allocatable,dimension(:)    :: grad_rho_c4  ! MPI d(ista_fftph:iend_fftph)
      integer  :: is, in
      real(kind=DP),pointer,dimension(:)  :: x
#ifdef _DETAIL_GGA_TIMING_
      integer :: id_sname = -1
#endif

!!$      if(myrank_ggacmp >= nrank_ggacmp) goto 1003

      if(nrank_ggacmp > 1)  allocate(grad_rho_c2(ista_fftph:iend_fftph),grad_rho_c4(ista_fftph:iend_fftph))   ! MPI

#ifndef _XC_SAVE_MEMORY_
      do is = 1, ispin
         grad_rho(:,is) = 0.d0
         do in = 1, 3
            if(map_ggacmp(in) /= myrank_ggacmp) cycle
            call dFxc_dgradrho_dot_gradrho2(rinplw,in,ispin,is)
            ! grad_{in}(rho_{is})*dFx/d|grad(rho_{is})| + grad_{in}(rho)*dFc/d|grad(rho)|  -> afft
            if(iprixc >= 2) write(nfout,'(" is, in = ",2i8)') is, in
            if(iprixc >= 2) call wd_small_part_of_afft(33,'R space before<m_FFT_CD_direct_c)',120)
            call m_FFT_CD_direct_c(nfout,afft)     ! afft (= q_{in}(G))
            if(iprixc >= 2) call wd_small_part_of_afft(33,'G space after <m_FFT_CD_direct_c)',120)
            call G_xyz_afft(in)                    ! G_{in}q_{in}(G)  -> afft
            call m_FFT_CD_inverse_c(nfout,afft)    ! afft -> afft
            if(iprixc >= 2) call wd_small_part_of_afft(33,'R space after <m_FFT_CD_direct_c)',120)
            call add_negative_afft_to_grad_rho(is) ! afft -> grad_rho
         end do
         if(nrank_ggacmp > 1) then
            grad_rho_c2(:) = grad_rho(:,is)
            call mpi_allreduce(grad_rho_c2,grad_rho_c4,iend_fftph-ista_fftph+1 &
                 & ,mpi_double_precision, mpi_sum, mpi_ggacmp_cross_world(myrank_cdfft),ierr)
            grad_rho(:,is) = grad_rho_c4(:)
         end if
      end do
#else
      if(ispin == 1) then
         is = 1
         grad_rho(:,is) = 0.d0
         do in = 1, 3
            if(map_ggacmp(in) /= myrank_ggacmp) cycle
            call g_xyz_chden_l(in,is)          ! G_xyz * rho(G) --> afft
            call m_FFT_CD_inverse_c(nfout,afft)  ! (-i)*d(rho(r))/d(x|y|z)
            call cp_afft_to_cggawk13           ! -afft --> cggawk13
            call dFxc_dgradrho_dot_gradrho(rinplw) ! grad(rho)*dFxc/d|grad(rho)| -> afft
            ! Calculation i-th component df/d(d|grad n|) in recip. space
            call m_FFT_CD_direct_c(nfout,afft)
            call G_xyz_afft(in)
            call m_FFT_CD_inverse_c(nfout,afft)
            call add_negative_afft_to_grad_rho(1)
            !  --> grad_rho
         end do
         if(nrank_ggacmp > 1) then
            grad_rho_c2(:) = grad_rho(:,is)
#ifdef _DETAIL_GGA_TIMING_
            call tstatc0_begin('mpi_allreduce(ggaxcp) ',id_sname)
#endif
            call mpi_allreduce(grad_rho_c2,grad_rho_c4,iend_fftph-ista_fftph+1 &
                 & ,mpi_double_precision, mpi_sum, mpi_ggacmp_cross_world(myrank_cdfft),ierr)
#ifdef _DETAIL_GGA_TIMING_
            call tstatc0_end(id_sname)
#endif
            grad_rho(:,is) = grad_rho_c4(:)
         end if
      else if(ispin == 2) then
         do is = 1, 2
            grad_rho(:,is) = 0.d0
            do in = 1, 3
               if(map_ggacmp(in) /= myrank_ggacmp) cycle
               call g_xyz_chden_l(in,is)          ! G_xyz * rho_{is}(G) --> afft
               call m_FFT_CD_inverse_c(nfout,afft)  ! (-i)*d(rho_{is}(r))/d(x|y|z)
               call cp_afft_to_cggawk13           ! -afft --> cggawk13
               x => dF_dgradrho(ista_fftph:iend_fftph,is)
               call x_dot_cggawk13_into_afft(rinplw,x)
               !  grad(rho) * dFx/d|grad(rho)|

               ! Calculation i-th component df/d(d|grad n|) in recip. space
               call m_FFT_CD_direct_c(nfout,afft)
               call G_xyz_afft(in)
               call m_FFT_CD_inverse_c(nfout,afft)
               call add_negative_afft_to_grad_rho(is)
               !  --> grad_rho
            end do
         end do
         !--> correlation part
         do in = 1, 3
            if(map_ggacmp(in) /= myrank_ggacmp) cycle
            call g_xyz_total_chden_l(in)
            call m_FFT_CD_inverse_c(nfout,afft) ! (-i)*d(rho(r))/d(x|y|z)
            call cp_afft_to_cggawk13        ! -afft --> cggawk13
            x => grad_trho(ista_fftph:iend_fftph)
            call x_dot_cggawk13_into_afft(rinplw,x)
            call m_FFT_CD_direct_c(nfout,afft)
            call G_xyz_afft(in)             !
            call m_FFT_CD_inverse_c(nfout,afft)
            do is = 1, 2
               call add_negative_afft_to_grad_rho(is)
            end do
         end do
         if(nrank_ggacmp > 1) then
            Do is=1, 2
               grad_rho_c2(:) = grad_rho(:,is)
#ifdef _DETAIL_GGA_TIMING_
               call tstatc0_begin('mpi_allreduce(ggaxcp) ',id_sname)
#endif
               call mpi_allreduce(grad_rho_c2,grad_rho_c4,iend_fftph-ista_fftph+1 &
                    & ,mpi_double_precision, mpi_sum, mpi_ggacmp_cross_world(myrank_cdfft),ierr)
#ifdef _DETAIL_GGA_TIMING_
               call tstatc0_end(id_sname)
#endif
               grad_rho(:,is) = grad_rho_c4(:)
            End Do
         end if
      end if
#endif
!!$1003 continue
      if(nrank_ggacmp > 1)  deallocate(grad_rho_c4,grad_rho_c2)    ! MPI

    end subroutine dFxc_over_ddgradrho

! ============================= added by K. Tagami ================= 11.0
    subroutine dFxc_over_ddgradrho_noncl()
      real(kind=DP),allocatable,dimension(:)    :: grad_rho_c2
                                                   ! MPI d(ista_fftph:iend_fftph)
      real(kind=DP),allocatable,dimension(:)    :: grad_rho_c4
                                                   ! MPI d(ista_fftph:iend_fftph)
      integer  :: is, in
      real(kind=DP),pointer,dimension(:)  :: x

#ifdef _DETAIL_GGA_TIMING_
      integer :: id_sname = -1
#endif

      if ( nrank_ggacmp >1 ) then
         allocate( grad_rho_c2(ista_fftph:iend_fftph) )
         allocate( grad_rho_c4(ista_fftph:iend_fftph) )       ! MPI
      endif

#ifndef _XC_SAVE_MEMORY_
      grad_rho = 0.0d0

      do is = 1, ispin
         do in = 1, 3
            if (map_ggacmp(in) /= myrank_ggacmp) cycle
            call dFxc_dgradrho_dot_gradrho2(rinplw,in,ispin,is)
                                     ! grad_{in}(rho_{is})*dFx/d|grad(rho_{is})|
                                     ! + grad_{in}(rho)*dFc/d|grad(rho)|  -> afft
            if (iprixc >= 2) write(nfout,'(" is, in = ",2i8)') is, in
            if (iprixc >= 2) then
               call wd_small_part_of_afft(33,'R space before<m_FFT_CD_direct_c)',120)
            endif

            call m_FFT_CD_direct_c(nfout,afft)     ! afft (= q_{in}(G))
            if (iprixc >= 2) then
               call wd_small_part_of_afft(33,'G space after <m_FFT_CD_direct_c)',120)
            endif

            call G_xyz_afft(in)                    ! G_{in}q_{in}(G)  -> afft
            call m_FFT_CD_inverse_c(nfout,afft)    ! afft -> afft
            if (iprixc >= 2) then
               call wd_small_part_of_afft(33,'R space after <m_FFT_CD_direct_c)',120)
            endif

            call add_negative_afft_to_grad_rho(is) ! afft -> grad_rho
         end do

         if (nrank_ggacmp > 1) then
            grad_rho_c2(:) = grad_rho(:,is)
            call mpi_allreduce( grad_rho_c2, grad_rho_c4, iend_fftph -ista_fftph +1, &
                 &              mpi_double_precision, mpi_sum, &
                 &              mpi_ggacmp_cross_world(myrank_cdfft), ierr )
            grad_rho(:,is) = grad_rho_c4(:)
         end if
      end do
#else
      allocate( bfft_kt( ista_fftp:iend_fftp,ndim_magmom )); bfft_kt = 0.0d0

      grad_rho = 0.0d0

      do in = 1, 3
         if(map_ggacmp(in) /= myrank_ggacmp) cycle

         bfft_kt = 0.0d0
         Do is=1, ndim_magmom
            call g_xyz_chden_l(in,is)          ! G_xyz * rho_{is}(G) --> afft
            call m_FFT_CD_inverse_c(nfout,afft)  ! (-i)*d(rho_{is}(r))/d(x|y|z)
            bfft_kt(:,is) = afft(:)
         End do

         call m_ES_GradCorr_Along_QuantzAxis2( RhoMag_R, bfft_kt, &
              &                                quantz_axis_inversion_flg_mesh )

         do is = 1, nspin
            afft(:) = bfft_kt(:,is)
            call cp_afft_to_cggawk13           ! -afft --> cggawk13
            x => dF_dgradrho(ista_fftph:iend_fftph,is)
            call x_dot_cggawk13_into_afft(rinplw,x)
                                               !  grad(rho) * dFx/d|grad(rho)|

            ! Calculation i-th component df/d(d|grad n|) in recip. space
            call m_FFT_CD_direct_c(nfout,afft)
            call G_xyz_afft(in)
            call m_FFT_CD_inverse_c(nfout,afft)
            call add_negative_afft_to_grad_rho(is)
            !  --> grad_rho
         end do
      end do

!--> correlation part
      do in = 1, 3
         if(map_ggacmp(in) /= myrank_ggacmp) cycle

         is = 1
         call g_xyz_chden_l(in,is)          ! G_xyz * rho_{is}(G) --> afft
         call m_FFT_CD_inverse_c(nfout,afft)  ! (-i)*d(rho_{is}(r))/d(x|y|z)

         call cp_afft_to_cggawk13        ! -afft --> cggawk13
         x => grad_trho(ista_fftph:iend_fftph)
         call x_dot_cggawk13_into_afft(rinplw,x)
         call m_FFT_CD_direct_c(nfout,afft)
         call G_xyz_afft(in)             !
         call m_FFT_CD_inverse_c(nfout,afft)
         do is = 1, nspin
            call add_negative_afft_to_grad_rho(is)
         end do
      end do

      if (nrank_ggacmp > 1) then
         Do is=1, nspin
            grad_rho_c2(:) = grad_rho(:,is)
#ifdef _DETAIL_GGA_TIMING_
            call tstatc0_begin('mpi_allreduce(ggaxcp) ',id_sname)
#endif
            call mpi_allreduce( grad_rho_c2, grad_rho_c4, iend_fftph -ista_fftph +1, &
                 &              mpi_double_precision, mpi_sum, &
                 &              mpi_ggacmp_cross_world(myrank_cdfft), ierr )

#ifdef _DETAIL_GGA_TIMING_
            call tstatc0_end(id_sname)
#endif
            grad_rho(:,is) = grad_rho_c4(:)
         End do
      end if

      deallocate( bfft_kt )
#endif

      if(nrank_ggacmp > 1)  deallocate(grad_rho_c4,grad_rho_c2)    ! MPI

    end subroutine dFxc_over_ddgradrho_noncl
! ======================================================================= 11.0

    subroutine finally_gga_xc_pot(is)
      integer, intent(in) :: is
      integer       :: m

! ==== KT_add ==== 2014/08/04
      if ( noncol ) then
         if ( sw_neglect_low_helicity == ON ) then
            call m_ES_neglect_low_Helicity( RhoMag_R, MagVorticity, grad_rho(:,is) )
         endif
         if ( sw_neglect_low_vorticity == ON ) then
            call m_ES_neglect_low_Vorticity( RhoMag_R, MagVorticity, grad_rho(:,is) )
         endif
      endif
! ================ 2014/08/04

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

    subroutine add_negative_afft_to_grad_rho(is)
      integer, intent(in) :: is
      integer             :: m
      do m = ista_fftph, iend_fftph      ! MPI
         grad_rho(m,is) = grad_rho(m,is) - afft(2*m-1)
      end do
      if(iprixc >= 2) then
         write(6,'(" !XC -- afft -- ,is = ",i8," <<add_negative_afft_to_grad_rho>>")') is
         write(6,'(" !XC ",6f12.6)') (afft(m*2-1),m=ista_fftph,ista_fftph+17)

         write(6,'(" !XC -- grad_rho -- ,is = ",i8," <<add_negative_afft_to_grad_rho>>")') is
         write(6,'(" !XC ",6f12.6)') (grad_rho(m,is),m=ista_fftph,ista_fftph+17)
      end if

    end subroutine add_negative_afft_to_grad_rho

    subroutine G_xyz_afft(in)
      integer, intent(in) :: in
      integer       :: n, i,j,k,ip
      real(kind=DP) :: gxyz

      do n = ista_fftp, iend_fftp   ! MPI
         gxyz = rltv(in,1)*inx(n)+rltv(in,2)*jnx(n)+rltv(in,3)*knx(n)
         afft(n) = gxyz*afft(n)
      enddo
      call boundary_zero_into_afft(in)  ! -(contained in subr. m_XC_cal_potential)

!!$      do k = 1, lz
!!$         do j = min(ly-ny_d*myrank_cdfft,ny_d)+1, ly_d
!!$            do i = 1, lx*2
!!$               ip = 2*((j-1)*lx+(k-1)*lx*ly_d)+ista_fftp-1 + i
!!$               afft(ip) = 0.d0
!!$            end do
!!$         end do
!!$      end do
!!$      do k = nz_d*npes_cdfft+1, lz
!!$         do j = 1, ly_d
!!$            do i = 1, lx*2
!!$               ip = 2*((j-1)*lx+(k-1)*lx*ly_d)+ista_fftp-1 + i
!!$               afft(ip) = 0.d0
!!$            end do
!!$         end do
!!$      end do
!!$      do k = 1, lz
!!$         do j = 1, ly_d
!!$            do i = fft_box_size_CD(1,1)+1, fft_box_size_CD_c(1,0)
!!$               ip = 2*((j-1)*lx+(k-1)*lx*ly_d)+ista_fftp-1 + 2*i-1
!!$               afft(ip) = 0.d0
!!$               afft(ip+1) = 0.d0
!!$            end do
!!$         end do
!!$      end do
    end subroutine G_xyz_afft

#ifndef _XC_SAVE_MEMORY_
    subroutine dFxc_dgradrho_dot_gradrho2(rinplw,in,ispin,is)
      real(kind=DP), intent(in) :: rinplw
      integer, intent(in)       :: ispin,in,is
      integer                   :: m
      afft = 0.d0
      if(ispin == 2) then
         do m = ista_fftph, iend_fftph
            afft(2*m) = rinplw*(dF_dgradrho(m,is)*cgrad_rho(m,in,is) &
                 &            +grad_trho(m)*(cgrad_rho(m,in,1)+cgrad_rho(m,in,2)) )
         end do
      else if(ispin == 1) then
         do m = ista_fftph, iend_fftph
            afft(2*m) = rinplw*(dF_dgradrho(m,is)+grad_trho(m))*cgrad_rho(m,in,is)
         end do
      end if
    end subroutine dFxc_dgradrho_dot_gradrho2

#else
    subroutine dFxc_dgradrho_dot_gradrho(rinplw)
      real(kind=DP), intent(in) :: rinplw
      integer                   :: m
      afft = 0.d0
      do m = ista_fftph, iend_fftph
         afft(2*m) = rinplw*(dF_dgradrho(m,1)+grad_trho(m))*cggawk13(m)
      end do
    end subroutine dFxc_dgradrho_dot_gradrho

    subroutine x_dot_cggawk13_into_afft(rinplw,x)
      real(kind=DP), intent(in) :: rinplw
      real(kind=DP), intent(in) :: x(ista_fftph:iend_fftph)
      integer                   :: m
      afft = 0.d0
      do m = ista_fftph, iend_fftph   ! MPI
         afft(2*m) = rinplw*x(m) *cggawk13(m)
      end do
    end subroutine x_dot_cggawk13_into_afft

    subroutine cp_afft_to_cggawk13
      integer  :: m

      do m = ista_fftph, iend_fftph  ! MPI
         cggawk13(m) = - afft(2*m)
      end do
    end subroutine cp_afft_to_cggawk13
#endif


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

! ==== KT_add ======= 13.0XX
    subroutine g2_xyz_chden_l(is)
      integer, intent(in) :: is
      integer       :: n
      real(kind=DP) :: gx, gy, gz, g2

      afft = 0.d0
      do n = ista_fftp, iend_fftp  ! MPI
         gx = rltv(1,1)*inx(n) +rltv(1,2)*jnx(n) +rltv(1,3)*knx(n)
         gy = rltv(2,1)*inx(n) +rltv(2,2)*jnx(n) +rltv(2,3)*knx(n)
         gz = rltv(3,1)*inx(n) +rltv(3,2)*jnx(n) +rltv(3,3)*knx(n)
         g2 = gx**2 +gy**2 +gz**2
         afft(n) = g2 *chden_l(n,is)
      enddo

    end subroutine g2_xyz_chden_l
! =================== 13.0XX

    subroutine check_lmn_even
      integer             :: nlmn, i

      do i = 1, 3
         nlmn = fft_box_size_CD(i,1)/2
         if(2*nlmn == fft_box_size_CD(i,1)) then
            lmn_even(i) = .true.
         else
            lmn_even(i) = .false.
         end if
      end do
    end subroutine check_lmn_even

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

    subroutine initialize_cggawk_arrays
      grad_trho = 0.d0
#ifndef _XC_SAVE_MEMORY_
      cgrad_rho = 0.d0
#else
      cggawk13 = 0.d0
#endif
      grad_rho = 0.d0
      dF_drho = 0.d0
      dF_dgradrho = 0.d0
    end subroutine initialize_cggawk_arrays

!!$    subroutine cp_chgrhr_onto_afft(is)
!!$      integer, intent(in) :: is
!!$      integer   :: i
!!$
!!$      afft = 0.d0
!!$      do i = ista_fftph, iend_fftph   ! MPI
!!$         afft(2*i-1) = chgrhr_l(i,is)*rinplw
!!$      end do
!!$    end subroutine cp_chgrhr_onto_afft

#ifndef _XC_SAVE_MEMORY_
    subroutine cp_afft_to_cgrad_rho(is,in)
      integer, intent(in) :: is, in
      integer :: i
      do i = ista_fftph, iend_fftph
         cgrad_rho(i,in,is) = -afft(i*2)
      end do
    end subroutine cp_afft_to_cgrad_rho

#endif

! ===== KT_add ====== 13.0XX
    subroutine cp_afft_to_grad2_rho(is)
      integer, intent(in) :: is
      integer :: i

      do i = ista_fftph, iend_fftph
         grad2_rho(i,is) = -afft(2*i-1)
      end do
    end subroutine cp_afft_to_grad2_rho
! =================== 13.0XX

! ================================== modified by K. Tagami =============== 11.0
!    subroutine cp_afft_to_chgrhr
!
    subroutine cp_afft_to_chgrhr(iloop)
      integer, intent(in) :: iloop
! ======================================================================== 11.0

      integer  :: i
      if(myrank_ggacmp < nrank_ggacmp) then
         do i = ista_fftph, iend_fftph    ! MPI
            chgrhr_l(i,iloop) = afft(i*2-1)
         end do
      else
         do i = ista_fftph, iend_fftph    ! MPI
            chgrhr_l(i,iloop) = 0.d0
         end do
      end if

      if(iprixc >= 2) then
! =============================== added by K. Tagami ============ 11.0
         write(nfout,*) "!XC "
         write(nfout,*) "!XC         iloop = ", iloop
! =============================================================== 11.0
         write(nfout,'(" !XC -- chgrhr_l -- <<cp_afft_to_chgrhr>>")')

! ============================== modified by K. Tagami ============= 11.0
!         write(nfout,'(" !XC ",6f12.6)') (chgrhr_l(i,1),i=ista_fftph,ista_fftph+17)
         write(nfout,'(" !XC ",6f12.6)') (chgrhr_l(i,iloop),i=ista_fftph,ista_fftph+17)
! ================================================================== 11.0

         write(nfout,'(" !XC -- afft --     <<cp_afft_to_chgrhr>>")')
         write(nfout,'(" !XC ",6f12.6)') (afft(i*2-1),i=ista_fftph,ista_fftph+17)
      end if
    end subroutine cp_afft_to_chgrhr

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

! =============================== added by K. Tagami ==================== 11.0
    subroutine wd_small_part_of_vxc_ssrep

      integer :: i

!!$      if(mype /= 0) return
      if(kimg==1) then
         write(nfout,*) 'vxc'
         if(input_charge == Partial_Core_Charge) then
            do i = ista_kngp, min(ista_kngp+14, iend_kngp)
               write(nfout,200) vxcpc_l(i,1), i
            end do
         else
            do i = ista_kngp, min(ista_kngp+14, iend_kngp)
               write(nfout,201) vxc_l(i,1,1),vxc_l(i,1,ndim_chgpot), i
            end do
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
            do i = ista_kngp, min(ista_kngp+29,iend_kngp)
               write(nfout,202) vxc_l(i,1   ,1),vxc_l(i,1,   ndim_chgpot), i,1
               write(nfout,202) vxc_l(i,kimg,1),vxc_l(i,kimg,ndim_chgpot), i,2
            end do
         end if
203      format(' ','(',d12.4,',',d12.4,')',i8)
202      format(' ','(',d12.4,',',d12.4,')',i8, ' kimg= ',i3)
         write(nfout,*) 'exc'
         write(nfout,200) exc
      end if
    end subroutine wd_small_part_of_vxc_ssrep
! ================================================================= 11.0

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

! ==================================== modified by K. Tagami =============== 11.0
!!!    subroutine map_charge_onto_a_fft_box(chgq_l)
!!      real(DP),intent(in) :: chgq_l(ista_kngp:iend_kngp,kimg,nspin)

    subroutine map_charge_onto_a_fft_box( chgq_l, ndim_dens, iloop, operation_mode )
      integer, intent(in) :: ndim_dens, iloop
      integer, intent(in) :: operation_mode
      real(DP),intent(in) :: chgq_l(ista_kngp:iend_kngp,kimg,ndim_dens)
! ========================================================================== 11.0

      real(kind=DP)      :: fac
      integer            :: mm, it, i, j, ip
      real(kind=DP),allocatable,dimension(:)    :: afft_mpi2,afft_mpi3   ! MPI d(mp_fftp)
      real(kind=DP),allocatable,dimension(:)    :: afft_mpi1,afft_mpi4   ! MPI d(nfftp)
      integer :: id_sname = -1
      call tstatc0_begin('map_charge_onto_a_fft_box ',id_sname)

      allocate(afft_mpi1(nfftp))

! =========================== modified by K. Tagami ============== 11.0
!      if(input_charge == Valence_plus_PC_Charge) then
!         fac = 1.d0/nspin
!      else if(input_charge == Partial_Core_Charge) then
!         fac = 1.d0
!      else
!        stop ' Error of ixc in <<<map_charge_onto_a_fft_box>>>'
!      endif

      if ( operation_mode == Valence_plus_PC_Charge ) then
         fac = 1.d0/nspin
      else if ( operation_mode == Valence_plus_Core_Charge ) then
         fac = 1.d0/nspin
      else if ( operation_mode == Valence_Charge_Only ) then
         fac = 1.d0
      else if ( operation_mode == Partial_Core_Charge ) then
         fac = 1.d0
      else
        call phase_error_with_msg(nfout,' Error of ixc in <<<map_charge_onto_a_fft_box>>>',__LINE__,__FILE__)
      endif
! ================================================================= 11.0

      mm = 0
      afft_mpi1 = 0.d0

! ================================= added by K. Tagami =================== 11.0
      if ( operation_mode == Valence_Charge_Only ) goto 1200
! ======================================================================== 11.0

      do it = 1, ntyp
         if(itpcc(it) == 0) cycle
         mm = mm + 1
         do j = 1, kimg
#ifdef NEC_TUNE_MXCP
!CDIR NODEP
!CDIR INNER
#endif
            if ( operation_mode /= Valence_plus_Core_Charge ) then
               do i = ista_kngp, min(kgp_reduced, iend_kngp)  !for mpi
#ifdef _MPIFFT_
                  ip = (igfp_l_c(i)-1)*kimg + j
#else
                  ip = (igfp_l(i)-1)*kimg + j
#endif
                  afft_mpi1(ip) = afft_mpi1(ip) &   !mpi
                       & + fac*zfm3_l(i,it,j)*rhpcg_l(i,mm)   !mpi
               end do
            else
               do i = ista_kngp, min(kgp_reduced, iend_kngp)  !for mpi
#ifdef _MPIFFT_
                  ip = (igfp_l_c(i)-1)*kimg + j
#else
                  ip = (igfp_l(i)-1)*kimg + j
#endif
                  afft_mpi1(ip) = afft_mpi1(ip) &   !mpi
                       & + fac*zfm3_l(i,it,j)*rhpcg_l(i,mm)   !mpi
               end do
            end if
         end do
      end do

! ==================================== added by K. Tagami ============== 11.0
1200  continue
! ====================================================================== 11.0

! ================================== modified by K. Tagami ============== 11.0
!!      if(input_charge == Valence_plus_PC_Charge) then
      if ( operation_mode /= Partial_Core_Charge ) then
! ======================================================================== 11.0

         do j = 1, kimg
#ifdef NEC_TUNE_MXCP
!CDIR NODEP
!CDIR INNER
#endif
            do i = ista_kngp, min(kgp_reduced, iend_kngp)  !for mpi
#ifdef _MPIFFT_
               ip = (igfp_l_c(i)-1)*kimg + j
#else
               ip = (igfp_l(i)-1)*kimg + j
#endif
               afft_mpi1(ip) = afft_mpi1(ip) + chgq_l(i,j,iloop) !mpi
            end do
         end do
      end if

      if(npes_cdfft > 1) then
         call mpi_barrier(MPI_CommGroup,ierr)
         allocate(afft_mpi2(mp_fftp))
         allocate(afft_mpi3(mp_fftp))
         do j = 0, npes_cdfft-1
#ifdef _MPIFFT_
            call set_cdata(j,afft_mpi1,afft_mpi2)
#else
            do i = nis_fftp(j), nie_fftp(j)
               afft_mpi2(i-nis_fftp(j)+1) = afft_mpi1(i)
            end do
#endif
            call mpi_allreduce(afft_mpi2,afft_mpi3,mp_fftp &
                 &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)

            if(j == myrank_cdfft) then
               do i = ista_fftp, iend_fftp
                  afft(i) = afft_mpi3(i - ista_fftp + 1)
               end do
            end if
         end do
         deallocate(afft_mpi3,afft_mpi2)
      else
         if(nrank_ggacmp > 1) then
            allocate(afft_mpi4(nfftp))
            call mpi_allreduce(afft_mpi1, afft_mpi4, nfftp &
                 & , mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
            afft = afft_mpi4
            deallocate(afft_mpi4)
         else
            afft = afft_mpi1
         end if
      end if

      deallocate(afft_mpi1)
      call tstatc0_end(id_sname)
    end subroutine map_charge_onto_a_fft_box

! ============================= modified by K. Tagami ============== 11.0
!    subroutine map_charge_onto_a_fft_box2(nfout,chgq_l)
!      integer, intent(in) :: nfout
!      real(DP),intent(in) :: chgq_l(ista_kngp:iend_kngp,kimg,nspin)
!
    subroutine map_charge_onto_a_fft_box2( nfout, chgq_l, ndim_dens, iloop, &
	&                                  operation_mode )
      integer, intent(in) :: nfout, ndim_dens, iloop
      integer, intent(in) :: operation_mode
      real(DP),intent(in) :: chgq_l(ista_kngp:iend_kngp,kimg,ndim_dens)
! ================================================================== 11.0

      real(kind=DP)      :: fac
      integer            :: mm, it, i, j
      real(kind=DP),allocatable,dimension(:)    :: tmp_s, tmp_r ! MPI d(nfftp/npes)
      real(kind=DP),allocatable,dimension(:,:)  :: chgqplus_l ! MPI d(ista_kngp:iend_kngp)
!      integer, save, allocatable, dimension(:)  :: ip ! d(ista_kngp:iend_kngp,3)
      integer, allocatable, dimension(:)        :: igfp_s, igfp_r ! d(nfftp/npes)
!      integer, save, allocatable, dimension(:,:):: igfp_ijk ! d(ista_kngp:iend_kngp+(0,1),3)
      integer, allocatable, dimension(:)        :: ipout_s, ipout_r ! d(nfftp/npes)
      integer, allocatable, dimension(:)        :: np ! d(0:npes_cdfft-1)
!      integer, save, allocatable, dimension(:,:):: np_send ! d(0:npes-1,0:npes-1)
      integer, allocatable, dimension(:,:)      :: np_tmp ! d(0:npes-1,0:npes-1)
      integer :: pe_s, pe_r, datasize_s, datasize_r, nfp, igp, nb, npc, maxdatasize
      integer :: nstart, nend, lx, lxy, ipout, j0, k, l, np0, idp, icdfft
      integer, allocatable, dimension(:) :: req_r,req_s

      integer :: ia
      real(kind=DP) :: weight, grt

      integer :: id_sname = -1
#ifdef _DETAIL_MAP_CHARGE_TIMING_
      integer :: id_sname2 = -1, id_sname3=-1, id_sname4=-1, id_sname5=-1
#endif
      call tstatc0_begin('map_charge_onto_a_fft_box2 ',id_sname)

#ifdef _DETAIL_MAP_CHARGE_TIMING_
      call tstatc0_begin('map_charge2fftbox2_part1 ',id_sname2)
#endif

      allocate(chgqplus_l(ista_kngp:iend_kngp,kimg))

! ============================== modified by K. Tagami ================ 11.0
!      if(input_charge == Valence_plus_PC_Charge) then
!         fac = 1.d0/nspin
!      else if(input_charge == Partial_Core_Charge) then
!         fac = 1.d0
!      else
!        stop ' Error of ixc in <<<map_charge_onto_a_fft_box2>>>'
!      endif

      if ( operation_mode == Valence_plus_PC_Charge) then
         fac = 1.d0/nspin
      else if ( operation_mode == Valence_plus_Core_Charge) then
         fac = 1.d0/nspin
      else if ( operation_mode == Valence_Charge_Only ) then
         fac = 1.d0
      else if( operation_mode == Partial_Core_Charge ) then
         fac = 1.d0
      else
        call phase_error_with_msg(nfout,' Error of ixc in <<<map_charge_onto_a_fft_box2>>>',__LINE__,__FILE__)
      endif
! ====================================================================== 11.0

      chgqplus_l = 0.d0
      mm = 0

! ================================= added by K. Tagami =================== 11.0
      if ( operation_mode == Valence_Charge_Only ) goto 1200
! ======================================================================== 11.0

      do it = 1, ntyp
         if(itpcc(it) == 0) cycle
         mm = mm + 1
         do j = 1, kimg
#ifdef NEC_TUNE_MXCP
!CDIR NODEP
!CDIR INNER
#endif
            if ( operation_mode /= Valence_plus_Core_Charge ) then
               do i = ista_kngp, iend_kngp  !for mpi
                  chgqplus_l(i,j) = chgqplus_l(i,j) &
                       & + fac*zfm3_l(i,it,j)*rhpcg_l(i,mm)   !mpi
               end do
            else
               do i = ista_kngp, iend_kngp  !for mpi
                  chgqplus_l(i,j) = chgqplus_l(i,j) &
                       & + fac*zfm3_l(i,it,j)*rhcg_l(i,mm)   !mpi
               end do
            endif
         end do
      end do

1200  continue

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

! ================================== modified by K. Tagami ============== 11.0
!!      if(input_charge == Valence_plus_PC_Charge) then
      if ( operation_mode /= Partial_Core_Charge ) then
! ======================================================================== 11.0

         do j = 1, kimg
#ifdef NEC_TUNE_MXCP
!CDIR NODEP
!CDIR INNER
#endif
            do i = ista_kngp, iend_kngp  !for mpi
               chgqplus_l(i,j) = chgqplus_l(i,j) + chgq_l(i,j,iloop)
            end do
         end do
      end if

      afft = 0.d0
      if(npes > 1) then
         if(.not.np_send_is_set) then
            allocate(ip(ista_kngp:iend_kngp)); ip = -1
            allocate(np_send(0:npes-1,0:npes-1))
            allocate(np(0:npes-1)); np = 0
         end if
!!$         allocate(ip(ista_kngp:iend_kngp)); ip = -1
!!$         allocate(np(0:npes-1)); np = 0
!!$         allocate(np_send(0:npes-1,0:npes-1))
      end if

#ifdef _MPIFFT_
      if(.not.np_send_is_set) then
         if(mod(iend_kngp-ista_kngp+1,2) == 1) then
            allocate(igfp_ijk(ista_kngp:min(kgp_reduced,iend_kngp),3))
         else
            allocate(igfp_ijk(ista_kngp:min(kgp_reduced,iend_kngp)+1,3))
         end if
         idp = fft_box_size_CD_c(1,0)
         do l = ista_kngp, min(kgp_reduced,iend_kngp)
            np0 = igfp_l_c(l)
            i   = mod(np0-1,idp)+1
            j   = mod((np0-i)/idp, ly) + 1
            k   = (np0-i-(j-1)*idp)/(idp*ly)+1
            igfp_ijk(l,1)  = i
            igfp_ijk(l,2)  = j
            igfp_ijk(l,3)  = k
         end do
      end if
      lx     = fft_box_size_CD_c(1,0)
      lxy    = lx*ly_d
#endif

      if(npes > 1) then
         if(.not.np_send_is_set) then
            np_send = 0
            np = 0
            do k = 0, npes-1
               if(map_pe2ggacmp(k) >= nrank_ggacmp) cycle
#ifdef _MPIFFT_
               icdfft = map_pe2cdfft(k)
               nstart = icdfft*ny_d+1
               nend   = min(fft_box_size_CD(2,1),(icdfft+1)*ny_d)
               do l = ista_kngp, min(kgp_reduced, iend_kngp)
                  j   = igfp_ijk(l,2)
                  if(j >= nstart .and. j <= nend) then
                     if(map_pe2ggacmp(k) == 0) ip(l) = k
                     np(k) = np(k) + 1
                  end if
               end do
#else
               icdfft = mod(k,npes_cdfft)
               do l = ista_kngp, min(kgp_reduced, iend_kngp)
                  igp = (igfp_l(l)-1)*kimg+1
                  if(igp >= nis_fftp(icdfft) .and. igp <= nie_fftp(icdfft)) then
                     ip(l) = map_pe2cdfft(k)
                     np(k) = np(k) + 1
                  end if
               end do
#endif
               if(iprixc >=2)  write(nfout,'(" np(",i5,") = ",i8)') k, np(k)
            end do

            np_send(:,mype) = np(:)  ! number element sent from mype to j

            allocate(np_tmp(0:npes-1,0:npes-1))
            call mpi_allreduce(np_send,np_tmp,npes*npes,mpi_integer, mpi_sum, MPI_CommGroup,ierr)
            np_send = np_tmp
            deallocate(np_tmp)

         end if
         maxdatasize = 0
         do j = 0, npes-1
            do i = 0, npes-1
               if(maxdatasize < np_send(i,j)) maxdatasize = np_send(i,j)
            end do
         end do
         if(iprixc>=2) write(nfout,'(" maxdatasize = ",i8)') maxdatasize
      end if
      np_send_is_set = .true.

#ifdef _DETAIL_MAP_CHARGE_TIMING_
      call tstatc0_end(id_sname2)
      call tstatc0_begin('map_charge2fftbox2_part2 ',id_sname3)
#endif

#ifdef _MPIFFT_
      allocate(ipout_s(ista_kngp:min(kgp_reduced,iend_kngp)))
      do l = ista_kngp, min(kgp_reduced,iend_kngp)
         i = igfp_ijk(l,1)
         j = igfp_ijk(l,2)
         k = igfp_ijk(l,3)
         icdfft = (j-1)/ny_d
         j0 = j - ny_d*icdfft
         ipout_s(l) = kimg*(i + (j0-1)*lx + (k-1)*lxy - 1)  + nis_fftp(icdfft)
      end do
#endif

      ! mype -> mype
      if(myrank_ggacmp < nrank_ggacmp) then
#ifdef _MPIFFT_
         nstart = myrank_cdfft*ny_d+1
         nend   = min(fft_box_size_CD(2,1),(myrank_cdfft+1)*ny_d)
         if(kimg == 1) then
            if(npes > 1) then
               do l = ista_kngp, min(kgp_reduced,iend_kngp)
                  i   = igfp_ijk(l,1)
                  if(i < 1 .or. i > fft_box_size_CD(1,1)) cycle
                  j   = igfp_ijk(l,2)
                  if(j >= nstart .and. j <= nend) then
                     ipout = ipout_s(l)
                     afft(ipout) = chgqplus_l(l,1)
                  end if
               end do
            else
               do l = ista_kngp, min(kgp_reduced,iend_kngp)
                  ipout = ipout_s(l)
                  afft(ipout) = chgqplus_l(l,1)
               end do
            end if
         else
            if(npes > 1) then
               do l = ista_kngp, min(kgp_reduced,iend_kngp)
                  i   = igfp_ijk(l,1)
                  if(i < 1 .or. i > fft_box_size_CD(1,1)) cycle
                  j   = igfp_ijk(l,2)
                  if(j >= nstart .and. j <= nend) then
                     ipout = ipout_s(l)
                     afft(ipout)   = chgqplus_l(l,1)
                     afft(ipout+1) = chgqplus_l(l,2)
                  end if
               end do
            else
               do l = ista_kngp, min(kgp_reduced,iend_kngp)
                  ipout = ipout_s(l)
                  afft(ipout)   = chgqplus_l(l,1)
                  afft(ipout+1) = chgqplus_l(l,2)
               end do
            end if
         end if
#else
         if(kimg == 1) then
            do l = ista_kngp, min(kgp_reduced,iend_kngp)
               igp = igfp_l(l)
               if(igp >= ista_fftp .and. igp <= iend_fftp) then
                  afft(igp) = chgqplus_l(l,1)
               end if
            end do
         else
#ifdef NEC_TUNE_MXCP
!CDIR NODEP
#endif
            do l = ista_kngp, min(kgp_reduced,iend_kngp)
               igp = (igfp_l(l)-1)*kimg+1
               if(igp >= ista_fftp .and. igp <= iend_fftp) then
                  afft(igp)   = chgqplus_l(l,1)
                  afft(igp+1) = chgqplus_l(l,2)
               end if
            end do
         end if
#endif
      end if

#ifdef _DETAIL_MAP_CHARGE_TIMING_
      call tstatc0_end(id_sname3)
      call tstatc0_begin('map_charge2fftbox2_part3 ',id_sname4)
#endif

      if(npes > 1) then
         allocate(tmp_s(maxdatasize*kimg))
         allocate(tmp_r(maxdatasize*kimg))
         allocate(igfp_s(maxdatasize),igfp_r(maxdatasize))
         allocate(req_r(npes-1))
         allocate(req_s(npes-1))
         do nb = 1, npes-1
            pe_s = mod(mype+nb,      npes)
            pe_r = mod(mype-nb+npes, npes)
#ifdef _MPIFFT_
            icdfft = map_pe2cdfft(pe_s)
#endif

            if(map_pe2ggacmp(pe_s) < nrank_ggacmp .and. np_send(pe_s,mype) > 0) then
               tmp_s = 0.d0
               igfp_s = 0
               nfp = 0

#ifdef _DETAIL_MAP_CHARGE_TIMING_
          call tstatc0_begin('map_charge2fftbox2_part4 ',id_sname5)
#endif
          if(kimg==1) then
               do l = ista_kngp, min(kgp_reduced,iend_kngp)
                  if(ip(l) == map_pe2cdfft(pe_s)) then
                     nfp = nfp+1
#ifdef _MPIFFT_
                     igfp_s(nfp) = ipout_s(l)
#else
                     igfp_s(nfp) = kimg*(igfp_l(l)-1)+1
#endif
                        tmp_s(nfp) = chgqplus_l(l,1)
                  end if
               end do
          else
               do l = ista_kngp, min(kgp_reduced,iend_kngp)
                  if(ip(l) == map_pe2cdfft(pe_s)) then
                     nfp = nfp+1
#ifdef _MPIFFT_
                     igfp_s(nfp) = ipout_s(l)
#else
                     igfp_s(nfp) = kimg*(igfp_l(l)-1)+1
#endif
                        tmp_s(2*nfp-1) = chgqplus_l(l,1)
                        tmp_s(2*nfp)   = chgqplus_l(l,2)
                  end if
               end do
          endif

#ifdef _DETAIL_MAP_CHARGE_TIMING_
          call tstatc0_end(id_sname5)
#endif

               datasize_s = np_send(pe_s,mype)*kimg
               if(iprixc>=2) then
                  write(nfout,'(" nb, pe_s,pe_r = ",3i8)') nb, pe_s,pe_r
                  write(nfout,'(" nfp*kimg   = ", i8)') nfp*kimg
                  write(nfout,'(" datasize_s = ",i8)') datasize_s
               end if
               if(datasize_s /= nfp*kimg) call phase_error_with_msg(nfout,' datasize_s /= nfp*kimg',__LINE__,__FILE__)
               call mpi_isend(tmp_s,datasize_s,mpi_double_precision, &
                    & pe_s, mype, MPI_CommGroup,req_s(nb),ierr)
            end if

            if(myrank_ggacmp < nrank_ggacmp .and. np_send(mype,pe_r) > 0) then
               datasize_r = np_send(mype,pe_r)*kimg
               if(iprixc>=2) write(nfout,'(" datasize_r = ",i8)') datasize_r
               call mpi_irecv(tmp_r,datasize_r,mpi_double_precision, &
                    & pe_r, pe_r,  MPI_CommGroup,req_r(nb),ierr)
            end if
            if(map_pe2ggacmp(pe_s) < nrank_ggacmp .and. np_send(pe_s,mype) > 0) &
                 & call mpi_wait(req_s(nb),istatus,ierr)
            if(myrank_ggacmp < nrank_ggacmp .and. np_send(mype,pe_r) > 0) &
                 & call mpi_wait(req_r(nb),istatus,ierr)
            if(iprixc>=2) write(nfout,'( " tmp_s , tmp_r  completed")')

            datasize_s = np_send(pe_s,mype)
            datasize_r = np_send(mype,pe_r)
            if(map_pe2ggacmp(pe_s) < nrank_ggacmp .and. np_send(pe_s,mype) > 0) then
               call mpi_isend(igfp_s,datasize_s,mpi_integer, &
                    & pe_s, mype, MPI_CommGroup, req_s(nb),ierr)
            end if
            if(myrank_ggacmp < nrank_ggacmp .and. np_send(mype,pe_r) > 0) then
               call mpi_irecv(igfp_r,datasize_r,mpi_integer, &
                    & pe_r, pe_r,  MPI_CommGroup, req_r(nb),ierr)
            end if
            if(map_pe2ggacmp(pe_s) <nrank_ggacmp .and. np_send(pe_s,mype) > 0) &
                 & call mpi_wait(req_s(nb),istatus,ierr)
            if(myrank_ggacmp < nrank_ggacmp .and. np_send(mype,pe_r) > 0) &
                 & call mpi_wait(req_r(nb),istatus,ierr)
            if(iprixc>=2) write(nfout,'( " igfp_s , igfp_r  completed")')

            if(myrank_ggacmp < nrank_ggacmp) then
               if(kimg == 1) then
                  do l = 1, np_send(mype,pe_r)
                     ipout = igfp_r(l)
                     afft(ipout) = tmp_r(l)
                  end do
               else
#ifdef NEC_TUNE_MXCP
!CDIR NODEP
#endif
                  do l = 1, np_send(mype,pe_r)
                     ipout = igfp_r(l)
                     afft(ipout)   = tmp_r(2*l-1)
                     afft(ipout+1) = tmp_r(2*l)
                  end do
               end if
            end if
         end do
         deallocate(req_s, req_r)
         deallocate(igfp_r, igfp_s)
         deallocate(tmp_r,tmp_s)
      end if
#ifdef _MPIFFT_
      deallocate(ipout_s)
!!$      deallocate(igfp_ijk)
#endif
      if(npes > 1) then
         if(allocated(np)) deallocate(np)
!!$         deallocate(np_send,np, ip)
      end if
      deallocate(chgqplus_l)
#ifdef _DETAIL_MAP_CHARGE_TIMING_
      call tstatc0_end(id_sname4)
#endif

      call tstatc0_end(id_sname)

    end subroutine map_charge_onto_a_fft_box2

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

    subroutine cpafft_to_vxc_or_vxcpc( is )
      integer, intent(in) :: is
      integer :: i, ik, ip
      real(kind=DP),allocatable,dimension(:)    :: afft_mpi1   ! MPI d(nfftp)
      integer :: id_sname = -1
      call tstatc0_begin('cpafft_to_vxc_or_vxcpc ',id_sname)

      allocate(afft_mpi1(nfftp))

      call afft_allgatherv(afft,afft_mpi1)

      do ik = 1, kimg
         if(input_charge == Partial_Core_Charge) then
            do i = ista_kngp, min(kgp_reduced,iend_kngp)  ! mpi
#ifdef _MPIFFT_
               ip = (igfp_l_c(i)-1)*kimg + ik
#else
               ip = (igfp_l(i)-1)*kimg + ik
#endif
               vxcpc_l(i,ik) = afft_mpi1(ip)*rinplw
            end do
         else if(input_charge == Valence_plus_PC_Charge) then
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

    subroutine cpafft_to_vtau(is,afft)
      integer, intent(in) :: is
      real(kind=DP), intent(in) :: afft(ista_fftp,iend_fftp)

      integer :: i, ik, ip
      real(kind=DP),allocatable,dimension(:)    :: afft_mpi1   ! MPI d(nfftp)
      integer :: id_sname = -1
      call tstatc0_begin('cpafft_to_vtau ',id_sname)

      allocate(afft_mpi1(nfftp))
      call afft_allgatherv(afft,afft_mpi1)

      do ik = 1, kimg
         if(input_charge == Partial_Core_Charge) then

         else if(input_charge == Valence_plus_PC_Charge) then
            do i = ista_kngp, min(kgp_reduced,iend_kngp)  ! mpi
#ifdef _MPIFFT_
               ip = (igfp_l_c(i)-1)*kimg + ik
#else
               ip = (igfp_l(i)-1)*kimg + ik
#endif
               vtau_l(i,ik,is)   = afft_mpi1(ip)*rinplw
            end do
         end if
      end do
      deallocate(afft_mpi1)
      call tstatc0_end(id_sname)
    end subroutine cpafft_to_vtau

    subroutine cpafft(is)
      integer, intent(in) :: is
      integer             :: i

      afft = 0.d0
      do i = ista_fftph, iend_fftph                  ! MPI
         afft(2*i-1) = chgrhr_l(i,is)
      end do

    end subroutine cpafft

    subroutine cp_vals_to_afft_rspace1( is, vals, afft )
      integer, intent(in) :: is
      real(kind=DP), intent(in) :: vals(ista_fftph:iend_fftph,nspin)
      real(kind=DP), intent(out) :: afft(ista_fftp:iend_fftp)

      integer             :: i

      afft = 0.d0
      do i = ista_fftph, iend_fftph                  ! MPI
         afft(2*i-1) = vals(i,is)
      end do
    end subroutine cp_vals_to_afft_rspace1

#ifdef __EDA__
! -----  ascat stars modifying  -----
    subroutine cp_exc_for_EDA
      use m_FFT, only : nfftp_nonpara
      integer             :: i, ik, ip, ierr
      integer :: ix,iy,iz,ip1,ip2,idp,nd2p,nd3p
      integer, dimension(3) :: n1,n2
      real(kind=DP) :: fac
      n1(1:3) = fft_box_size_CD_nonpara(1:3,0)
      n2(1:3) = fft_box_size_CD_c(1:3,0)
      exc_on_a_grid = 0.d0

      do iz=1,n1(3)
        do iy=1,n1(2)
          do ix=1,n1(1)
            ip1 = (iz-1)*n1(1)*n1(2)+(iy-1)*n1(1)+ix
            ip2 = (iz-1)*n2(1)*n2(2)+(iy-1)*n2(1)+ix
            if (ip2>=ista_fftph .and. ip2<=iend_fftph) then
              exc_on_a_grid(2*ip1-1) = exc_on_a_grid_wk(ip2)*univol/f2or1(ip2)
            endif
          enddo
        enddo
      enddo

      call mpi_allreduce(mpi_in_place, exc_on_a_grid, nfftp_nonpara, mpi_double_precision,mpi_sum, &
      &                  MPI_CommGroup, ierr)

    end subroutine cp_exc_for_EDA
! -----  ascat ceases modifying  -----
#endif

  end subroutine m_XC_cal_potential

  subroutine xcpotf(ispin,input_charge)
    integer, intent(in) :: ispin,input_charge

! #1) 1994/11/08 by T.Yamasaki
!    Coding for the case of xctype='PERZUN ' and  'XALPHA ' are done.
! #2) Spin-polarization is introduced by T. Yamasaki at 15th Dec. 1994
! #3) f77 -> f90     4th April 1999  by T. Yamasaki

    real(kind=DP) :: rinplw
    real(kind=DP) :: exc_mpi     ! MPI

    real(kind=DP) :: DELTA
    data DELTA/1.d-40/

    rinplw = 1.d0/product(fft_box_size_CD(1:3,1))

    exc = 0.d0
    if(xctype == 'wign   '.or. xctype == 'wigner ')then
       call xcpotf_wigner(ispin,ista_fftph,iend_fftph,input_charge,DELTA,chgrhr_l,f2or1,exc)
    else if(  xctype ==  'pzold  ') then
       call xcpotf_pzold(ispin,ista_fftph,iend_fftph,input_charge,DELTA,chgrhr_l,f2or1,exc)
    else if(  xctype ==  'xalfa  ') then
       call xcpotf_xalfa(ispin,ista_fftph,iend_fftph,input_charge,DELTA,chgrhr_l,f2or1,exc)
    else if(  xctype == 'perzun '.or. xctype == 'pz     ') then
       call xcpotf_pz(ispin,ista_fftph,iend_fftph,input_charge,DELTA,chgrhr_l,f2or1,exc)
    else if(  xctype == 'vwn    ') then
       call xcpotf_vwn(ispin,ista_fftph,iend_fftph,input_charge,DELTA,chgrhr_l,f2or1,exc)
    else if(  xctype=='mjw    '.or. xctype=='bh     ' .or.xctype=='gl     ') then
       call xcpotf_mjw_bh_gl(len_xctype,xctype,ispin,ista_fftph,iend_fftph,input_charge,DELTA,chgrhr_l,f2or1,exc)
    endif
                        ! all xcpotf_* subroutines are in -(b_XC_Potential) ->chgrhr_l,exc
    exc = exc*univol*rinplw

    if(npes >= 2) then               ! MPI
       call mpi_allreduce(exc,exc_mpi,1,mpi_double_precision,mpi_sum &
            & ,MPI_CommGroup,ierr)  ! MPI
       exc = exc_mpi                 ! MPI
    end if

  end subroutine xcpotf

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


  subroutine m_XC_rst_npsend()
     np_send_is_set = .false.
     if(allocated(ip)) deallocate(ip)
     if(allocated(igfp_ijk)) deallocate(igfp_ijk)
     if(allocated(np_send)) deallocate(np_send)
  end subroutine m_XC_rst_npsend

end module m_XC_Potential_2D
