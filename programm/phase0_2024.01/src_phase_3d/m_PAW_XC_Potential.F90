!#define VDWDF_PAW_BEFORE_2017

#ifdef NEC_TIMER
#  define START_TIMER(a) call start_timer(a)
#  define STOP_TIMER(a)  call stop_timer(a)
#else
#  define START_TIMER(a)
#  define STOP_TIMER(a)
#endif
!=======================================================================
!
!  SOFTWARE NAME : PHASE ($Revision: 633 $)
!
!  MODULE: m_PAW_XC_Potential
!
!  AUTHOR(S): T. Yamasaki, T. Yamamoto and T. Ohno November/2009
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
!#========================================================================
!
! Bug fix:  2016/10/13
!    *Antiferromagnetic calculation
!        The variable flg_done is neglected for the moment
!        because the skipping mechanism does not work properly in some cases.
!
!#========================================================================
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
module m_PAW_XC_Potential
    use m_Const_Parameters,     only : DP,GGA,Valence_plus_PC_Charge &
         &                            ,VXC_AND_EXC,PAI4,PAI2, NEWTON, LAGRANGE, EXC_ONLY, DELTA10 &
         &                            ,SphericalHarmonicsExpansion, GaussLegendre
    use m_Control_Parameters,   only : nspin,xctype,len_xctype,printable,ipripp,iprixc &
         &                           , paw_density_gradient, af, exchange_pot_type
    use m_Timing,               only : tstatc0_begin, tstatc0_end
    use m_PseudoPotential,      only : mmesh,wf_mnrc,radr_paw &
         &                            ,m_PP_find_maximum_l &
         &                            ,dion_vxc,ipaw &
         &                            ,ilmt,ltp,mtp,taup &
         &                            ,nltpw,iltpw,lppw,tppw,index_lmt2lt &
         &                            ,il2p,isph,dl2p,iqitg &
         &                            ,wf_mnrc,psirpw,phirpw,qrspspw,xh &
         &                            ,ia2ia_symmtry_op
    use m_XC_Potential,         only : check_of_xctype
    use m_PAW_ChargeDensity,    only : ntheta,nphi,omg_wght &
                                        ,m_PAWCD_set_ae_cd &
                                        ,m_PAWCD_set_ps_cd &
                                        ,m_PAWCD_set_ae_der_cd &
                                        ,m_PAWCD_set_ps_der_cd &
                                        ,m_PAWCD_set_ae_cd_sym &
                                        ,m_PAWCD_set_ps_cd_sym &
                                        ,m_PAWCD_set_ae_der_cd_sym &
                                        ,m_PAWCD_set_ps_der_cd_sym &
                                        ,ylm,dylm_dth,dylm_dph &
                                        ,cos_theta &
                                        ,m_PAWCD_set_ae_cd_sphex2 &
                                        ,m_PAWCD_set_ps_cd_sphex2 &
                                        ,m_PAWCD_set_sq_der_cd_sdphex2 &
                                        ,m_PAWCD_set_cr2_isph2_mmt2 &
                                        ,paw_cr2,paw_isph2,paw_mmt2,paw_dnr &
                                        ,surface_integral_method &
                                        ,calcGaussLegendreIntegration &
                                        ,calcSphericalHarmonicsExpansion
    use m_Ionic_System,         only : ntyp,natm,ityp,iwei
    use m_Crystal_Structure,    only : nopr
    use m_Parallelization,      only : mype,MPI_CommGroup, npes &
                                      , ista_atm_ke, iend_atm_ke, nrank_g, myrank_g


! ========================= added by K. Tagami ================ 11.0
#ifdef __EDA__
  use m_Control_Parameters,   only : noncol, ndim_magmom, ON, sw_eda
#else
  use m_Control_Parameters,   only : noncol, ndim_magmom, ON
#endif
  use m_PAW_ChargeDensity,    only : m_PAWCD_ae_cd_sphex2_nonclA, &
       &                             m_PAWCD_ps_cd_sphex2_nonclA, &
       &                             max_sph_expansion
  use m_Parallelization,      only : ierr
! ============================================================= 11.0

! ============================= added by K. Tagami =========== 11.0
  use m_Ionic_System,        only : magmom_local_now
  use m_PseudoPotential,     only : nlmt
  use m_ES_NonCollinear,     only :  Global_Quantz_Axis_now
  use m_Crystal_Structure,   only :  level_of_projection_paw_charge
! ============================================================ 11.0

#ifdef LIBXC
  use m_Control_Parameters,  only : xc_family_exch, xc_family_corr
  use xc_f03_lib_m
#endif

#ifdef NEC_TIMER
    use nec_timer
#endif
!$$#if defined(PARA3D) && defined(PAW3D)
  use m_Parallelization,        only : nrank_natm, nrank_nrc, myrank_natm, myrank_nrc &
                                     , mpi_natm_world, mpi_nrc_world &
                                     , is_natm, ie_natm, nel_natm, is_nrc, ie_nrc, nel_nrc &
                                     , ista_natm, iend_natm, ne_natm, ista_nrc, iend_nrc, ne_nrc
! === DEBUG by tkato 2011/10/01 ================================================
! use z_interface_3D,           only : decomp_vxc_ae_k_r_3D
! ==============================================================================
! === For nrc decomposion. by takto 2012/12/05 =================================
    use m_PAW_ChargeDensity,    only : m_PAWCD_set_ae_cd_sphex2 &
                                      ,m_PAWCD_set_ps_cd_sphex2 &
                                      ,m_PAWCD_set_sq_der_cd_sdphex2
! ==============================================================================
! === For nrc decomposion. by takto 2012/12/07 =================================
    use m_Parallelization,      only : paw_no_work_to_do, paw_last_rank_on_nrc
! ==============================================================================
    use mpi
    implicit none
!    include 'mpif.h'
    private

    real(DP),allocatable,dimension(:,:,:,:):: vxc_ae_k
    real(DP),allocatable,dimension(:,:,:,:):: vxc_ps_k
!$$#if defined(PARA3D) && defined(PAW3D)
    real(DP),allocatable,dimension(:,:,:,:):: vxc_ae_k_3D
    real(DP),allocatable,dimension(:,:,:,:):: vxc_ps_k_3D
    real(DP),allocatable,dimension(:,:)  :: vxc_ae_m_ps
    real(DP)                                :: exc_ae,exc_ps,texc
    integer                                :: msph

    real(DP),allocatable,dimension(:,:)       :: grad_nae,grad_nps
    real(DP),allocatable,dimension(:)         :: grad_tnae,grad_tnps
    real(DP),allocatable,dimension(:)         :: dgrad_tnae_dr,dgrad_tnps_dr
    real(DP),allocatable,dimension(:,:)       :: dF_dnae,dF_dgradnae,dF_dgradnae_dr
    real(DP),allocatable,dimension(:,:)       :: dF_dnps,dF_dgradnps,dF_dgradnps_dr
    real(DP),allocatable,dimension(:)         :: wos

    real(DP),allocatable,dimension(:,:)       :: nae,nps
    real(DP),allocatable,dimension(:,:)       :: dnae_dr,dnae_dth,dnae_dph
    real(DP),allocatable,dimension(:,:)       :: dnps_dr,dnps_dth,dnps_dph
    real(DP),allocatable,dimension(:,:)       :: ddnae_ddr
    real(DP),allocatable,dimension(:,:)       :: ddnps_ddr

! ***** New member for spherical harmonics expansion *****

    real(DP),allocatable,dimension(:,:,:)     :: nae_sph,nps_sph
    real(DP),allocatable,dimension(:,:,:)     :: dnae_dr_sph,dnps_dr_sph
    real(DP),allocatable,dimension(:,:,:)     :: ddnae_ddr_sph,ddnps_ddr_sph
    real(DP),allocatable,dimension(:,:,:)     :: grad_nae2_sph,grad_nps2_sph
    real(DP),allocatable,dimension(:,:)       :: grad_tnae2_sph,grad_tnps2_sph
    real(DP),allocatable,dimension(:)         :: exc_ae_field,exc_ps_field
    real(DP),allocatable,dimension(:,:)       :: dFx_dnnae,dFx_dnnps
    real(DP),allocatable,dimension(:,:)       :: dFx_dngae,dFx_dngps
    real(DP),allocatable,dimension(:,:)       :: dFx_dggae,dFx_dggps
    real(DP),allocatable,dimension(:,:)       :: dFx_dnnnae,dFx_dnnnps
    real(DP),allocatable,dimension(:,:)       :: dFx_dnngae,dFx_dnngps
    real(DP),allocatable,dimension(:,:)       :: dFx_dnggae,dFx_dnggps
    real(DP),allocatable,dimension(:,:)       :: dFx_dgggae,dFx_dgggps
    real(DP),allocatable,dimension(:)         :: dF_dgradtnae,dF_dgradtnps
    real(DP),allocatable,dimension(:)         :: dFc_daa_ae,dFc_daa_ps
    real(DP),allocatable,dimension(:)         :: dFc_dbb_ae,dFc_dbb_ps
    real(DP),allocatable,dimension(:)         :: dFc_dgg_ae,dFc_dgg_ps
    real(DP),allocatable,dimension(:)         :: dFc_dab_ae,dFc_dab_ps
    real(DP),allocatable,dimension(:)         :: dFc_dag_ae,dFc_dag_ps
    real(DP),allocatable,dimension(:)         :: dFc_dbg_ae,dFc_dbg_ps
    real(DP),allocatable,dimension(:)         :: dFc_daaa_ae,dFc_daaa_ps
    real(DP),allocatable,dimension(:)         :: dFc_dbbb_ae,dFc_dbbb_ps
    real(DP),allocatable,dimension(:)         :: dFc_dggg_ae,dFc_dggg_ps
    real(DP),allocatable,dimension(:)         :: dFc_daab_ae,dFc_daab_ps
    real(DP),allocatable,dimension(:)         :: dFc_daag_ae,dFc_daag_ps
    real(DP),allocatable,dimension(:)         :: dFc_dabb_ae,dFc_dabb_ps
    real(DP),allocatable,dimension(:)         :: dFc_dbbg_ae,dFc_dbbg_ps
    real(DP),allocatable,dimension(:)         :: dFc_dagg_ae,dFc_dagg_ps
    real(DP),allocatable,dimension(:)         :: dFc_dbgg_ae,dFc_dbgg_ps
    real(DP),allocatable,dimension(:)         :: dFc_dabg_ae,dFc_dabg_ps
    real(DP),allocatable,dimension(:,:)       :: nana_ae_sph,nana_ps_sph
    real(DP),allocatable,dimension(:,:)       :: nbnb_ae_sph,nbnb_ps_sph
    real(DP),allocatable,dimension(:,:)       :: gaga_ae_sph,gaga_ps_sph
    real(DP),allocatable,dimension(:,:)       :: gbgb_ae_sph,gbgb_ps_sph
    real(DP),allocatable,dimension(:,:)       :: gg_ae_sph,gg_ps_sph
    real(DP),allocatable,dimension(:,:)       :: nanb_ae_sph,nanb_ps_sph
    real(DP),allocatable,dimension(:,:)       :: naga_ae_sph,naga_ps_sph
    real(DP),allocatable,dimension(:,:)       :: nbgb_ae_sph,nbgb_ps_sph
    real(DP),allocatable,dimension(:,:)       :: nag_ae_sph,nag_ps_sph
    real(DP),allocatable,dimension(:,:)       :: nbg_ae_sph,nbg_ps_sph

    real(DP),allocatable,dimension(:,:)       :: dFxcdna_ae_sph,dFxcdna_ps_sph
    real(DP),allocatable,dimension(:,:)       :: dFxcdnb_ae_sph,dFxcdnb_ps_sph
    real(DP),allocatable,dimension(:,:)       :: dFxdgaovrga_ae_sph,dFxdgaovrga_ps_sph
    real(DP),allocatable,dimension(:,:)       :: dFxdgbovrgb_ae_sph,dFxdgbovrgb_ps_sph
    real(DP),allocatable,dimension(:,:)       :: dFcdgovrg_ae_sph,dFcdgovrg_ps_sph

    real(DP),allocatable,dimension(:)         :: dFadga_ae,   dFadga_ps
    real(DP),allocatable,dimension(:)         :: dFadg_ae,    dFadg_ps
    real(DP),allocatable,dimension(:)         :: dFadgaga_ae, dFadgaga_ps
    real(DP),allocatable,dimension(:)         :: dFadgg_ae,   dFadgg_ps
    real(DP),allocatable,dimension(:)         :: dFadnaga_ae, dFadnaga_ps
    real(DP),allocatable,dimension(:)         :: dFadnag_ae,  dFadnag_ps
!    real(DP),pointer,dimension(:)         :: dFadbg_ae,   dFadbg_ps

    real(DP),allocatable,dimension(:)         :: dFbdgb_ae,   dFbdgb_ps
    real(DP),allocatable,dimension(:)         :: dFbdg_ae,    dFbdg_ps
    real(DP),allocatable,dimension(:)         :: dFbdgbgb_ae, dFbdgbgb_ps
    real(DP),allocatable,dimension(:)         :: dFbdgg_ae,   dFbdgg_ps
    real(DP),allocatable,dimension(:)         :: dFbdnbgb_ae, dFbdnbgb_ps
    real(DP),allocatable,dimension(:)         :: dFbdnbg_ae,  dFbdnbg_ps
    real(DP),allocatable,dimension(:)         :: dFbdag_ae,   dFbdag_ps

    real(DP),allocatable,dimension(:)         :: dGadna_ae,   dGadna_ps
    real(DP),allocatable,dimension(:)         :: dGadga_ae,   dGadga_ps
    real(DP),allocatable,dimension(:)         :: dGadnana_ae, dGadnana_ps
    real(DP),allocatable,dimension(:)         :: dGadnaga_ae, dGadnaga_ps
    real(DP),allocatable,dimension(:)         :: dGadgaga_ae, dGadgaga_ps

    real(DP),allocatable,dimension(:)         :: dGbdnb_ae,   dGbdnb_ps
    real(DP),allocatable,dimension(:)         :: dGbdgb_ae,   dGbdgb_ps
    real(DP),allocatable,dimension(:)         :: dGbdnbnb_ae, dGbdnbnb_ps
    real(DP),allocatable,dimension(:)         :: dGbdnbgb_ae, dGbdnbgb_ps
    real(DP),allocatable,dimension(:)         :: dGbdgbgb_ae, dGbdgbgb_ps

    real(DP),allocatable,dimension(:)         :: dGdna_ae,    dGdna_ps
    real(DP),allocatable,dimension(:)         :: dGdnb_ae,    dGdnb_ps
    real(DP),allocatable,dimension(:)         :: dGdg_ae,     dGdg_ps
    real(DP),allocatable,dimension(:)         :: dGdnana_ae,  dGdnana_ps
    real(DP),allocatable,dimension(:)         :: dGdnbnb_ae,  dGdnbnb_ps
    real(DP),allocatable,dimension(:)         :: dGdgg_ae,    dGdgg_ps
    real(DP),allocatable,dimension(:)         :: dGdnanb_ae,  dGdnanb_ps
    real(DP),allocatable,dimension(:)         :: dGdnag_ae,   dGdnag_ps
    real(DP),allocatable,dimension(:)         :: dGdnbg_ae,   dGdnbg_ps

    logical,allocatable,dimension(:)         :: flg_done

    integer,allocatable,dimension(:)         :: irs

! ========================================= added by K. Tagami ================= 11.0
    real(DP), allocatable,dimension(:,:,:) :: magmom_local_ae
    real(DP), allocatable,dimension(:,:,:) :: magmom_local_ps
    real(DP), allocatable,dimension(:,:) :: magmom_local_wk

    real(DP), allocatable,dimension(:,:,:,:) :: rho_rad_ae
    real(DP), allocatable,dimension(:,:,:,:) :: rho_rad_ps
    real(DP), allocatable,dimension(:,:,:) :: rho_rad_wk
! ============================================================================== 11.0

    public:: m_PAW_XC_cal_potential
!!$    public:: m_PAW_XC_cal_potential_sym
    public:: m_PAW_XC_alloc_vxc
    public:: m_PAW_XC_get_dion_vxc
    public:: exc_ae,exc_ps

    public:: m_PAW_XC_cal_potential_sphex2
    public:: m_PAW_XC_dealloc_vxc
  public:: m_PAW_XC_allocation
  public:: m_PAW_XC_deallocation
  public:: vxc_ae_k_3D, vxc_ps_k_3D

! ======================- added by K. Tagami =============== 11.0
  public :: vxc_ae_k, vxc_ps_k
! ============================================================ 11.0

#ifdef __EDA__
    real(DP),allocatable,dimension(:):: exc_ae_on_atom, exc_ps_on_atom
    public :: exc_ae_on_atom, exc_ps_on_atom
#endif

interface
#ifdef __EDA__
  subroutine ex_ggapw91(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc &
     &                    ,dFx_drho,dFx_dgradrho, exc_on_a_grid_wk, ist, ien)
#else
  subroutine ex_ggapw91(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc &
     &                    ,dFx_drho,dFx_dgradrho, ist, ien)
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
  logical, intent(in),optional :: revPBE
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif
  end subroutine ex_ggapbe
#else
  subroutine ex_ggapbe(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dFx_drho,dFx_dgradrho,ien)
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
  real(kind=DP),intent(out),optional :: ecor
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif
  end subroutine cr_ggapbe
#else
  subroutine cr_ggapbe(nspin,ispin,chgrhr_l,grad_trho,f2or1,exc,dF_drho)
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
  end subroutine cr_ggapbe
#endif
  subroutine xclda(nspin,ispin,ista_r,iend_r,chgrhr_l,wos,exc,dF_drho,ien)
  use m_Const_Parameters,  only : DP, PAI
  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in),optional :: ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)
  real(kind=DP),intent(out) :: exc
!!$  real(kind=DP),intent(inout) :: dF_drho(ista_r:iend_r,nspin)
  real(kind=DP),intent(out) :: dF_drho(ista_r:iend_r,nspin)
  end subroutine xclda

  subroutine ggabek(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dF_drho,dF_dgradrho,ien)
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
  end subroutine ggabek

  subroutine ggaprd(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dF_drho,dF_dgradrho,ien)
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
  end subroutine ggaprd

  subroutine cr_lda(nspin,ispin,ista_r,iend_r,chgrhr_l,exc,dF_drho,ien)
  use m_Const_Parameters,  only : PAI,DP
  integer, intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in),optional :: ien
  real(kind=DP), intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP), intent(inout) :: exc
  real(kind=DP), intent(inout) :: dF_drho(ista_r:iend_r,ispin)
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
  real(kind=DP),intent(out),optional :: ecor
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif
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
end interface

contains

  subroutine m_PAW_XC_alloc_vxc
      integer:: n
      call m_PP_find_maximum_l(n)
      n=(n-1)+(n-1)+1
!!!      msph=n**2                 ! 2020/05/12
      msph=min( n**2, 25 )
      if(.not.allocated(vxc_ae_m_ps)) allocate(vxc_ae_m_ps(nspin,natm))
      vxc_ae_m_ps=0.d0

#ifdef __EDA__
      if ( sw_eda == ON ) then
         if ( .not. allocated(exc_ae_on_atom) ) then
            allocate( exc_ae_on_atom(natm) )
         endif
         if ( .not. allocated(exc_ps_on_atom) ) then
            allocate( exc_ps_on_atom(natm) )
         endif
      endif
#endif

        return
    end subroutine m_PAW_XC_alloc_vxc

  subroutine m_PAW_XC_dealloc_vxc
      if (allocated(vxc_ae_m_ps)) deallocate(vxc_ae_m_ps)

! ===================================== added by K. Tagami ================ 11.0
      if ( noncol ) then
         if ( allocated(magmom_local_ae) ) deallocate( magmom_local_ae )
         if ( allocated(magmom_local_ps) ) deallocate( magmom_local_ps )
         if ( allocated(rho_rad_ae) )      deallocate( rho_rad_ae )
         if ( allocated(rho_rad_ps) )      deallocate( rho_rad_ps )
!
         if ( allocated(magmom_local_wk) ) deallocate( magmom_local_wk )
         if ( allocated(rho_rad_wk) ) deallocate( rho_rad_wk )
      endif
! ========================================================================= 11.0
#ifdef __EDA__
      if ( allocated(exc_ae_on_atom) ) deallocate( exc_ae_on_atom )
      if ( allocated(exc_ps_on_atom) ) deallocate( exc_ps_on_atom )
#endif

  end subroutine m_PAW_XC_dealloc_vxc

  subroutine paw_xc_allocate
    allocate(nae(ista_nrc:iend_nrc,nspin));nae=0
    allocate(nps(ista_nrc:iend_nrc,nspin));nps=0
    allocate(wos(ista_nrc:iend_nrc));wos=0.d0
    if(check_of_xctype()==GGA) then
       allocate(grad_nae(ista_nrc:iend_nrc,nspin));grad_nae=0
       allocate(grad_tnae(ista_nrc:iend_nrc));grad_tnae=0
       allocate(dgrad_tnae_dr(ista_nrc:iend_nrc));dgrad_tnae_dr=0
       allocate(dF_dnae(ista_nrc:iend_nrc,nspin));dF_dnae=0
       allocate(dF_dgradnae(ista_nrc:iend_nrc,nspin));dF_dgradnae=0
       allocate(dF_dgradnae_dr(ista_nrc:iend_nrc,nspin));dF_dgradnae_dr=0
       allocate(grad_nps(ista_nrc:iend_nrc,nspin));grad_nps=0
       allocate(grad_tnps(ista_nrc:iend_nrc));grad_tnps=0
       allocate(dgrad_tnps_dr(ista_nrc:iend_nrc));dgrad_tnps_dr=0
       allocate(dF_dnps(ista_nrc:iend_nrc,nspin));dF_dnps=0
       allocate(dF_dgradnps(ista_nrc:iend_nrc,nspin));dF_dgradnps=0
       allocate(dF_dgradnps_dr(ista_nrc:iend_nrc,nspin));dF_dgradnps_dr=0
       allocate(dnae_dr(ista_nrc:iend_nrc,nspin));dnae_dr=0
       allocate(dnae_dth(ista_nrc:iend_nrc,nspin));dnae_dth=0
       allocate(dnae_dph(ista_nrc:iend_nrc,nspin));dnae_dph=0
       allocate(ddnae_ddr(ista_nrc:iend_nrc,nspin));ddnae_ddr=0
       allocate(dnps_dr(ista_nrc:iend_nrc,nspin));dnps_dr=0
       allocate(dnps_dth(ista_nrc:iend_nrc,nspin));dnps_dth=0
       allocate(dnps_dph(ista_nrc:iend_nrc,nspin));dnps_dph=0
       allocate(ddnps_ddr(ista_nrc:iend_nrc,nspin));ddnps_ddr=0
    end if
    return
  end subroutine paw_xc_allocate

  subroutine paw_xc_deallocate
      deallocate(nae)
      deallocate(nps)
      deallocate(wos)
      if(check_of_xctype()==GGA) then
          deallocate(grad_nae)
          deallocate(grad_tnae)
          deallocate(dF_dnae)
          deallocate(dF_dgradnae)
          deallocate(grad_nps)
          deallocate(grad_tnps)
          deallocate(dF_dnps)
          deallocate(dF_dgradnps)
          deallocate(dnae_dr)
          deallocate(dnae_dth)
          deallocate(dnae_dph)
          deallocate(ddnae_ddr)
          deallocate(dnps_dr)
          deallocate(dnps_dth)
          deallocate(dnps_dph)
          deallocate(ddnps_ddr)
          deallocate(dgrad_tnae_dr)
          deallocate(dgrad_tnps_dr)
          deallocate(dF_dgradnae_dr)
          deallocate(dF_dgradnps_dr)
      end if
      return
  end subroutine paw_xc_deallocate

  subroutine m_PAW_XC_allocation(nfout)
    integer :: nfout
    allocate(vxc_ae_k_3D(ista_nrc:iend_nrc,nspin,msph,ne_natm))
    allocate(vxc_ps_k_3D(ista_nrc:iend_nrc,nspin,msph,ne_natm))
    vxc_ae_k_3D=0.d0
    vxc_ps_k_3D=0.d0
  end subroutine m_PAW_XC_allocation

  subroutine m_PAW_XC_deallocation(nfout)
    integer :: nfout
    deallocate(vxc_ae_k_3D)
    deallocate(vxc_ps_k_3D)
  end subroutine m_PAW_XC_deallocation

  subroutine m_PAW_XC_cal_potential(nfout,vflag,flg_symmetry)
    integer,intent(in):: nfout,vflag
    logical,intent(in):: flg_symmetry

    integer,parameter :: PRINTLEVEL = 2
#ifdef _TIME_DETAIL_
    integer :: id_sname2 = -1, id_sname3=-1, id_sname4=-1,id_sname5=-1, ierr
#endif
    integer:: ia,it,ith,iph,iord, ia_add
    integer:: ir,is,nsph, ier
    integer:: nrc, tmpind
    real(DP):: costh,sinth,phi,dvec(3), radr, dylmt, dylmp
    real(DP):: dtnae_dr,dtnae_dth,dtnae_dph,ddtnae_ddr
    real(DP):: dtnps_dr,dtnps_dth,dtnps_dph,ddtnps_ddr ! h
    real(DP),allocatable, dimension(:,:) :: wost ! wost(mmesh,ntyp)
    real(kind=DP), allocatable, dimension(:,:) :: tae,tps ! === For nrc decomposion. by T. Yamasaki 2020/02/24 ===
    integer :: dnr
    integer :: id_sname = -1
!!$#if 0
!!$    real(kind=DP), allocatable, dimension(:,:) :: nae_g, dnae_dr_g, ddnae_ddr_g
!!$#endif

    call tstatc0_begin('m_PAW_XC_cal_potential',id_sname,level=1)
START_TIMER('m_PAW_XC_cal_potential')

    call paw_xc_allocate
    iord = paw_density_gradient%order

    exc_ae=0.d0
    exc_ps=0.d0

#ifdef __EDA__
    if ( sw_eda == ON ) then
       exc_ae_on_atom = 0.0d0;       exc_ps_on_atom = 0.0d0
    endif
#endif

    vxc_ae_k_3D=0.d0 ! == For nrc decomposion. by takto 2012/12/05 ==
    vxc_ps_k_3D=0.d0 ! == For nrc decomposion. by takto 2012/12/05 ==

    allocate(wost(mmesh,ntyp));wost=0.d0
    do it = 1, ntyp
       nrc = wf_mnrc(it)
       call set_weight_exp(ier,1,nrc,radr_paw(:,it),wost(:,it))
       do ir=1, nrc
          wost(ir,it) = wost(ir,it)*radr_paw(ir,it)**2
       end do
    end do

    if(iprixc >= PRINTLEVEL) then
       write(nfout,'(" ista_natm = ",i4," iend_natm = ",i4)') ista_natm, iend_natm
       write(nfout,'(" ntheta, nphi = ",2i8)') ntheta, nphi
    end if
    do ia = ista_natm, iend_natm   !     do ia=ista_atm_ke,iend_atm_ke
       ia_add = ia - ista_natm + 1

       it=ityp(ia)
       if(ipaw(it)/=1) cycle

       nrc = wf_mnrc(it)
       wos = 0.d0
       wos(ista_nrc:min(nrc,iend_nrc)) = wost(ista_nrc:min(nrc,iend_nrc),it)*iwei(ia)

       do ith=1,ntheta
          if(iprixc >= PRINTLEVEL) write(nfout,'(" ith = ",i4)') ith
          if(flg_symmetry) then
             costh=cos_theta(ith)
             sinth=sqrt(1.d0-costh**2)
          else
             sinth=sqrt(1.d0-cos_theta(ith)**2)
          end if
          do iph=1,nphi
!                    if(it .ne. ityp(ia)) cycle
#ifdef _TIME_DETAIL_
             call tstatc0_begin('m_PAW_XC_cal_pot(former)',id_sname4)
#endif
!!$             tmpind = (ith-1)*nphi+iph
!!$             if (mod(tmpind,nrank_g) /= myrank_g) cycle
             if(flg_symmetry) then
                phi=PAI2/dble(nphi)*dble(iph-1)
                dvec(1)=sinth*cos(phi)
                dvec(2)=sinth*sin(phi)
                dvec(3)=costh
                call m_PAWCD_set_ae_cd_sym( dvec,ia,nspin,nrc,nae,ista_nrc,iend_nrc)
                call m_PAWCD_set_ps_cd_sym( dvec,ia,nspin,nrc,nps,ista_nrc,iend_nrc)
             else
                call m_PAWCD_set_ae_cd(ith,iph,ia,nspin,nrc,nae,ista_nrc,iend_nrc)
                call m_PAWCD_set_ps_cd(ith,iph,ia,nspin,nrc,nps,ista_nrc,iend_nrc)
             end if
#ifdef _TIME_DETAIL_
             call tstatc0_end(id_sname4)
             call tstatc0_begin('m_PAW_XC_cal_pot(ggaxcp_paw)',id_sname3)
#endif
             if(check_of_xctype()==GGA) then
                call ggaxcp_paw()          ! --> texc,dF_dnae,dF_dgradnae,dF_dnps,dF_dgradnps
             else
                call xcpotf_paw(nrc,nspin,ith,Valence_plus_PC_Charge)
             end if
#ifdef _TIME_DETAIL_
             call tstatc0_end(id_sname3)
#endif
             if(iprixc >= PRINTLEVEL)  write(nfout,'(" iph/nphi = ",i4,"/",i4 &
                  &              ," , exc_ae, exc_ps = ",2d20.8 )') iph,nphi,exc_ae, exc_ps

             if(vflag /= VXC_AND_EXC) cycle
!!$             if(vflag == VXC_AND_EXC) then
             if(check_of_xctype()==GGA) then
#ifdef _TIME_DETAIL_
                call tstatc0_begin('m_PAW_XC_cal_pot(latter)',id_sname2)
#endif

                dnr = 1
                allocate(tae(ista_nrc-6*dnr:iend_nrc+6*dnr,nspin))   ! == For nrc decomposion. by T. Yamasaki 2020/02/24 ==
                call boundary_exchange_dim2(tae, DF_dgradnae, dnr, nspin)
                allocate(tps(ista_nrc-6*dnr:iend_nrc+6*dnr,nspin))
                call boundary_exchange_dim2(tps, DF_dgradnps, dnr, nspin)

                do is=1,nspin
                   call calc_diff_exp2_3D(1,iord,nrc,dnr,xh(it),radr_paw(:,it),tae(:,is),dF_dgradnae_dr(:,is) &
                        &    ,ista_nrc,iend_nrc)
                   call calc_diff_exp2_3D(1,iord,nrc,dnr,xh(it),radr_paw(:,it),tps(:,is),dF_dgradnps_dr(:,is) &
                        &    ,ista_nrc,iend_nrc)
                end do

                deallocate(tae)
                allocate(tae(ista_nrc-6*dnr:iend_nrc+6*dnr,1))   ! == For nrc decomposion. by T. Yamasaki 2020/02/24 ==
                call boundary_exchange_dim2(tae, grad_tnae, dnr, 1)
                deallocate(tps)
                allocate(tps(ista_nrc-6*dnr:iend_nrc+6*dnr,1))
                call boundary_exchange_dim2(tps, grad_tnps, dnr, 1)
                call calc_diff_exp2_3D(1,iord,nrc,dnr,xh(it),radr_paw(:,it),tae(:,1),dgrad_tnae_dr,ista_nrc,iend_nrc)
                call calc_diff_exp2_3D(1,iord,nrc,dnr,xh(it),radr_paw(:,it),tps(:,1),dgrad_tnps_dr,ista_nrc,iend_nrc)
                deallocate(tae,tps)
#ifdef _TIME_DETAIL_
                call tstatc0_end(id_sname2)
                call tstatc0_begin('latter2 ',id_sname5)
#endif
                do nsph=1,msph
                   dylmt = dylm_dth(ith,iph,nsph)
                   dylmp = dylm_dph(ith,iph,nsph)
                   do is=1,nspin
                      do ir=ista_nrc,iend_nrc  !    do ir=1,nrc

                         if(nspin==2) then
                            dtnae_dr = dnae_dr(ir,1) + dnae_dr(ir,2)
                            ddtnae_ddr = ddnae_ddr(ir,1) + ddnae_ddr(ir,2)
                            dtnae_dth = dnae_dth(ir,1) + dnae_dth(ir,2)
                            dtnae_dph = dnae_dph(ir,1) + dnae_dph(ir,2)
                            dtnps_dr = dnps_dr(ir,1) + dnps_dr(ir,2)
                            ddtnps_ddr = ddnps_ddr(ir,1) + ddnps_ddr(ir,2)
                            dtnps_dth = dnps_dth(ir,1) + dnps_dth(ir,2)
                            dtnps_dph = dnps_dph(ir,1) + dnps_dph(ir,2)
                         else
                            dtnae_dr = dnae_dr(ir,1)
                            ddtnae_ddr = ddnae_ddr(ir,1)
                            dtnae_dth = dnae_dth(ir,1)
                            dtnae_dph = dnae_dph(ir,1)
                            dtnps_dr = dnps_dr(ir,1)
                            ddtnps_ddr = ddnps_ddr(ir,1)
                            dtnps_dth = dnps_dth(ir,1)
                            dtnps_dph = dnps_dph(ir,1)

                         end if
                         radr = radr_paw(ir,it)
                         vxc_ae_k_3D(ir,is,nsph,ia_add) =  vxc_ae_k_3D(ir,is,nsph,ia_add) &
                              &  + ( ylm(ith,iph,nsph) &
                              &     * (dF_dnae(ir,is) - dF_dgradnae_dr(ir,is)*dnae_dr(ir,is) &
                              &        - dF_dgradnae(ir,is) * (ddnae_ddr(ir,is) + 2.d0/radr*dnae_dr(ir,is)) &
                              &       - dgrad_tnae_dr(ir) * dtnae_dr &
                              &        - grad_tnae(ir) * (ddtnae_ddr + 2.d0/radr*dtnae_dr) ) &
                              &    + dF_dgradnae(ir,is) &
                              &     * (dylmt*dnae_dth(ir,is) + dylmp*dnae_dph(ir,is)/sinth) /radr &
                              &    + grad_tnae(ir) * (dylmt*dtnae_dth + dylmp*dtnae_dph/sinth) /radr &
                              &    ) * omg_wght(ith)
                         vxc_ps_k_3D(ir,is,nsph,ia_add) =  vxc_ps_k_3D(ir,is,nsph,ia_add) &
                              &  + ( ylm(ith,iph,nsph) &
                              &     * (dF_dnps(ir,is) - dF_dgradnps_dr(ir,is)*dnps_dr(ir,is) &
                              &        - dF_dgradnps(ir,is) * (ddnps_ddr(ir,is) + 2.d0/radr*dnps_dr(ir,is)) &
                              &       - dgrad_tnps_dr(ir) * dtnps_dr &
                              &        - grad_tnps(ir) * (ddtnps_ddr + 2.d0/radr*dtnps_dr) ) &
                              &    + dF_dgradnps(ir,is) &
                              &        * (dylmt*dnps_dth(ir,is) + dylmp*dnps_dph(ir,is)/sinth)/ radr &
                              &    + grad_tnps(ir)* (  dylmt*dtnps_dth + dylmp*dtnps_dph/sinth)/ radr &
                              &    ) * omg_wght(ith)
                      end do
                   end do
                end do
#ifdef _TIME_DETAIL_
                call tstatc0_end(id_sname5)
#endif
             else

                do nsph=1,msph
                   do is=1,nspin
                      do ir=1,nrc
                         vxc_ae_k_3D(ir,is,nsph,ia_add) = vxc_ae_k_3D(ir,is,nsph,ia_add) &
                              & + nae(ir,is) * ylm(ith,iph,nsph) * omg_wght(ith)
                         vxc_ps_k_3D(ir,is,nsph,ia_add) = vxc_ps_k_3D(ir,is,nsph,ia_add) &
                              & + nps(ir,is) * ylm(ith,iph,nsph) * omg_wght(ith)
                      end do
                   end do
                end do
             end if
          end do
       end do
    end do

    exc_ae=exc_ae*PAI4
    exc_ps=exc_ps*PAI4
    vxc_ae_k_3D=vxc_ae_k_3D*PAI4
    vxc_ps_k_3D=vxc_ps_k_3D*PAI4

    if(iprixc >= PRINTLEVEL) then
       write(nfout,'(" exc_ae, exc_ps = ",2d20.8 )') exc_ae, exc_ps
    end if
    if(nrank_g>1) then
!!$       call mpi_allreduce(MPI_IN_PLACE,exc_ae,1,mpi_double_precision,mpi_sum,mpi_nrc_world,ier)
!!$       call mpi_allreduce(MPI_IN_PLACE,exc_ps,1,mpi_double_precision,mpi_sum,mpi_nrc_world,ier)
!!$#if 0
!!$       call mpi_allreduce(MPI_IN_PLACE,exc_ae,1,mpi_double_precision,mpi_sum,mpi_nrc_world,ier)
!!$#else
       call mpi_allreduce(MPI_IN_PLACE,exc_ae,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ier)
!!$#endif
       call mpi_allreduce(MPI_IN_PLACE,exc_ps,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ier)
    end if

    if(iprixc >= PRINTLEVEL) then
       write(nfout,'(" exc_ae, exc_ps = ",2d20.8 )') exc_ae, exc_ps
       write(nfout,'(" vxc_ae_k_3D(",i4,":",i4,",1,1,1) = ",6f12.4)') &
            & ista_nrc, ista_nrc+5,(vxc_ae_k_3D(ir,1,1,1),ir=ista_nrc,ista_nrc+5)
       write(nfout,'(" vxc_ps_k_3D(",i4,":",i4,",1,1,1) = ",6f12.4)') &
            & ista_nrc, ista_nrc+5,(vxc_ps_k_3D(ir,1,1,1),ir=ista_nrc,ista_nrc+5)
    end if

    deallocate(wost)
    call paw_xc_deallocate

STOP_TIMER('m_PAW_XC_cal_potential')
    call tstatc0_end(id_sname)
  contains

    subroutine ggaxcp0_paw
      integer :: id_sname = -1
      integer :: pot_type, ist, ien
      real(kind=DP) :: ecor, texc_ae, texc_ps
      real(kind=DP), allocatable :: tmp_work1(:), tmp_work2(:)
      real(kind=DP), allocatable :: dummy1(:)
!!$#if 0
!!$      real(kind=DP), allocatable, dimension(:,:) :: dF_dnae_g, dF_dgradnae_g,grad_nae_g
!!$      real(kind=DP), allocatable, dimension(:) :: wos_g, grad_tnae_g
!!$      integer :: ier
!!$#endif
      integer :: i,j

      call tstatc0_begin('ggaxcp0_paw ',id_sname)

      grad_tnae = 0.d0
      grad_nae = 0.d0
      dF_dnae = 0.d0
      dF_dgradnae = 0.d0
      grad_tnps = 0.d0
      grad_nps = 0.d0
      dF_dnps = 0.d0
      dF_dgradnps = 0.d0
      dnae_dr = 0.d0
      dnae_dth = 0.d0
      dnae_dph = 0.d0
      ddnae_ddr = 0.d0
      dnps_dr = 0.d0
      dnps_dth = 0.d0
      dnps_dph = 0.d0
      ddnps_ddr = 0.d0

!      if(xctype == 'ldapw91' .or. xctype == 'ldapbe '.or. xctype == 'vdwdf' ) then
      if(xctype == 'ldapw91' .or. xctype == 'ldapbe' ) then
         grad_nae=0.d0;grad_tnae=0.d0
         grad_nps=0.d0;grad_tnps=0.d0
      else
         call abs_grad_rho_up_down_total_paw_3D
      end if

      dnr = 1
      ist = (ista_nrc+dnr-2)/dnr
      ist = dnr*ist+1
      ien = min(iend_nrc, nrc)
      if ( ist <= 0 ) ist = 1

      if(xctype == 'ggapw91' .or. xctype == 'ldapw91') then
#ifdef __EDA__
         allocate( dummy1(ista_nrc:iend_nrc) )
         call ex_ggapw91( nspin, nspin, ista_nrc, iend_nrc, nae, grad_nae, &
              &           wos, texc_ae, dF_dnae, dF_dgradnae, dummy1, ist, ien )
         call cr_ggapw91( nspin, nspin, ista_nrc, iend_nrc, nae, grad_tnae, &
              &           wos, texc_ae, dF_dnae, dummy1, ist, ien )
         call ex_ggapw91( nspin, nspin, ista_nrc, iend_nrc, nps, grad_nps, &
              &           wos, texc_ps, dF_dnps, dF_dgradnps, dummy1, ist, ien )
         call cr_ggapw91( nspin, nspin, ista_nrc, iend_nrc, nps, grad_tnps, &
              &           wos, texc_ps, dF_dnps, dummy1, ist, ien )
         exc_ae = exc_ae +texc_ae*omg_wght(ith)
         exc_ps = exc_ps +texc_ps*omg_wght(ith)
         deallocate( dummy1 )
#else
         call ex_ggapw91( nspin, nspin, ista_nrc, iend_nrc, nae, grad_nae, &
              &           wos, texc_ae, dF_dnae, dF_dgradnae, ist, ien )
         call cr_ggapw91( nspin, nspin, ista_nrc, iend_nrc, nae, grad_tnae, &
              &           wos, texc_ae, dF_dnae, ist, ien )
         call ex_ggapw91( nspin, nspin, ista_nrc, iend_nrc, nps, grad_nps, &
              &           wos, texc_ps, dF_dnps, dF_dgradnps, ist, ien )
         call cr_ggapw91( nspin, nspin, ista_nrc, iend_nrc, nps, grad_tnps, &
              &           wos, texc_ps, dF_dnps, ist, ien )
         exc_ae = exc_ae +texc_ae*omg_wght(ith)
         exc_ps = exc_ps +texc_ps*omg_wght(ith)
#endif
!!      else if(xctype == 'ggapbe ' .or. xctype == 'ldapbe '.or. xctype == 'vdwdf') then
      else if(xctype == 'ggapbe ' .or. xctype == 'ldapbe ' ) then
!!$#if 0
!!$         allocate(wos_g(nrc));wos_g = 0.d0
!!$         allocate(grad_nae_g(nrc,nspin)); grad_nae_g = 0.d0
!!$         allocate(grad_tnae_g(nrc)); grad_tnae_g = 0.d0
!!$         do i = ista_nrc, min(nrc,iend_nrc)
!!$            wos_g(i) = wos(i)
!!$            grad_tnae_g(i) = grad_tnae(i)
!!$         end do
!!$         do j = 1,nspin
!!$            do i = ista_nrc, min(nrc,iend_nrc)
!!$               grad_nae_g(i,j) = grad_nae(i,j)
!!$            end do
!!$         end do
!!$         call mpi_allreduce(MPI_IN_PLACE,wos_g,nrc,mpi_double_precision,mpi_sum,mpi_natm_world,ier)
!!$         call mpi_allreduce(MPI_IN_PLACE,grad_tnae_g,nrc,mpi_double_precision,mpi_sum,mpi_natm_world,ier)
!!$         call mpi_allreduce(MPI_IN_PLACE,grad_nae_g,nrc*nspin,mpi_double_precision,mpi_sum,mpi_natm_world,ier)
!!$         allocate(dF_dnae_g(1:nrc,nspin))
!!$         allocate(dF_dgradnae_g(1:nrc,nspin))
!!$         call ex_ggapbe (nspin,nspin,1,nrc,nae_g(1:nrc,1:nspin),grad_nae_g, wos_g,texc,dF_dnae_g,dF_dgradnae_g,.false.)
!!$         call cr_ggapbe (nspin,nspin,1,nrc,nae_g(1:nrc,1:nspin),grad_tnae_g,wos_g,texc,dF_dnae_g,ecor)
!!$         do j = 1, nspin
!!$            do i = ista_nrc, min(nrc,iend_nrc)
!!$               dF_dnae(i,j) = dF_dnae_g(i,j)
!!$               dF_dgradnae(i,j) = dF_dgradnae_g(i,j)
!!$            end do
!!$         end do
!!$         deallocate(dF_dgradnae_g,dF_dnae_g)
!!$         deallocate(grad_tnae_g,grad_nae_g,wos_g)
!!$#else
!!$         grad_nae = 0.d0
!!$#endif

#ifdef __EDA__
         allocate( dummy1(ista_nrc:iend_nrc) )
         call ex_ggapbe( nspin, nspin, ista_nrc, iend_nrc, nae, grad_nae, &
              &          wos, texc_ae, dF_dnae, dF_dgradnae, &
              &          dummy1, .false., ist, ien )
         call cr_ggapbe( nspin, nspin, ista_nrc, iend_nrc, nae, grad_tnae, &
              &          wos, texc_ae, dF_dnae, dummy1, ecor, ist, ien )
         call ex_ggapbe( nspin, nspin, ista_nrc, iend_nrc, nps, grad_nps, &
              &          wos, texc_ps, dF_dnps, dF_dgradnps, &
              &          dummy1, .false., ist, ien )
         call cr_ggapbe( nspin, nspin, ista_nrc, iend_nrc, nps, grad_tnps, &
              &          wos, texc_ps, dF_dnps, dummy1, ecor, ist, ien )
         exc_ae = exc_ae +texc_ae*omg_wght(ith)
         exc_ps = exc_ps +texc_ps*omg_wght(ith)
         deallocate( dummy1 )
#else
         call ex_ggapbe( nspin, nspin, ista_nrc, iend_nrc, nae, grad_nae, &
              &          wos, texc_ae, dF_dnae, dF_dgradnae, &
              &          .false., ist, ien )
         call cr_ggapbe( nspin, nspin, ista_nrc, iend_nrc, nae, grad_tnae, &
              &          wos, texc_ae, dF_dnae, ecor, ist, ien )
         call ex_ggapbe( nspin, nspin, ista_nrc, iend_nrc, nps, grad_nps, &
              &          wos, texc_ps, dF_dnps, dF_dgradnps, &
              &          .false., ist, ien )
         call cr_ggapbe( nspin, nspin, ista_nrc, iend_nrc, nps, grad_tnps, &
              &          wos, texc_ps, dF_dnps, ecor, ist, ien )
         exc_ae = exc_ae +texc_ae*omg_wght(ith)
         exc_ps = exc_ps +texc_ps*omg_wght(ith)
#endif

      else if(xctype == 'katopbe' .or. xctype == 'ggapbek') then
#ifdef __EDA__
         allocate( dummy1(ista_nrc:iend_nrc) )
         call ex_ggapbe( nspin, nspin, ista_nrc, iend_nrc, nae, grad_nae, &
              &          wos, texc_ae, dF_dnae, dF_dgradnae, &
              &          dummy1, .false., ist, ien )
         call cr_ggapbe( nspin, nspin, ista_nrc, iend_nrc, nae, grad_tnae, &
              &          wos, texc_ae, dF_dnae, dummy1, ecor, ist, ien )
         call ex_ggapbe( nspin, nspin, ista_nrc, iend_nrc, nps, grad_nps, &
              &          wos, texc_ps, dF_dnps, dF_dgradnps, &
              &          dummy1, .false., ist, ien )
         call cr_ggapbe( nspin, nspin, ista_nrc, iend_nrc, nps, grad_tnps, &
              &          wos, texc_ps, dF_dnps, dummy1, ecor, ist, ien )
         exc_ae = exc_ae +texc_ae*omg_wght(ith)
         exc_ps = exc_ps +texc_ps*omg_wght(ith)
         deallocate( dummy1 )
#else
         call ex_ggapbe( nspin, nspin, ista_nrc, iend_nrc, nae, grad_nae, &
              &          wos, texc_ae, dF_dnae, dF_dgradnae, &
              &          .false., ist, ien )
         call cr_ggapbe( nspin, nspin, ista_nrc, iend_nrc, nae, grad_tnae, &
              &          wos, texc_ae, dF_dnae, ecor, ist, ien )
         call ex_ggapbe( nspin, nspin, ista_nrc, iend_nrc, nps, grad_nps, &
              &          wos, texc_ps, dF_dnps, dF_dgradnps, &
              &          .false., ist, ien )
         call cr_ggapbe( nspin, nspin, ista_nrc, iend_nrc, nps, grad_tnps, &
              &          wos, texc_ps, dF_dnps, ecor, ist, ien )
         exc_ae = exc_ae +texc_ae*omg_wght(ith)
         exc_ps = exc_ps +texc_ps*omg_wght(ith)
#endif

      else if(xctype == 'ggabp  ') then
         call xclda(nspin,nspin,ista_nrc,iend_nrc,nae,wos,texc,dF_dnae,nrc)
         call ggabek(nspin,nspin,ista_nrc,iend_nrc,nae,grad_nae,wos,texc,dF_dnae,dF_dgradnae,nrc)
         call ggaprd(nspin,nspin,ista_nrc,iend_nrc,nae,grad_nae,wos,texc,dF_dnae,dF_dgradnae,nrc)
         exc_ae=exc_ae+texc*omg_wght(ith)
         call xclda(nspin,nspin,ista_nrc,iend_nrc,nps,wos,texc,dF_dnps,nrc)
         call ggabek(nspin,nspin,ista_nrc,iend_nrc,nps,grad_nps,wos,texc,dF_dnps,dF_dgradnps,nrc)
         call ggaprd(nspin,nspin,ista_nrc,iend_nrc,nps,grad_nps,wos,texc,dF_dnps,dF_dgradnps,nrc)
         exc_ps=exc_ps+texc*omg_wght(ith)

      else if ( xctype == 'revpbe' .or. xctype == 'rpbe' .or. &
           &    xctype == 'wc06'   .or. xctype == 'htbs' .or. &
           &    xctype == 'pbesol' .or. xctype == 'pbeint' .or. &
           &    xctype == 'ev93' .or. xctype == 'evpw91' .or. &
           &    xctype == 'lb94' ) then

         if ( xctype == 'ggapbe' )  pot_type = 1
         if ( xctype == 'revpbe' )  pot_type = 2
         if ( xctype == 'rpbe' )    pot_type = 3
         if ( xctype == 'wc06' )    pot_type = 4
         if ( xctype == 'htbs' )    pot_type = 5
         if ( xctype == 'pbesol' )    pot_type = 6
         if ( xctype == 'pbeint' )    pot_type = 7
         if ( xctype == 'ev93' )    pot_type = 20
         if ( xctype == 'evpw91' )    pot_type = 20
         if ( xctype == 'lb94' )  pot_type = 25

#ifdef __EDA__
         allocate( dummy1(ista_nrc:iend_nrc) )
         call ex_gga_library( nspin, nspin, ista_nrc, iend_nrc, nae, grad_nae, &
              &               wos, texc_ae, dF_dnae, dF_dgradnae, dummy1, &
              &               pot_type, ist, ien )
         if ( xctype == 'ev93' .or. xctype == 'lb94' ) then
            grad_tnae = 0.d0
            call cr_gga_library( nspin, nspin, ista_nrc, iend_nrc, nae, grad_tnae, &
                 &               wos, texc_ae, dF_dnae, dummy1, ecor, &
                 &               pot_type, ist, ien )
         else if ( xctype == 'evpw91' ) then
            call cr_ggapw91( nspin, nspin, ista_nrc, iend_nrc, nae, grad_tnae, &
                 &           wos, texc_ae, dF_dnae, dummy1, ist, ien )
         else
            call cr_gga_library( nspin, nspin, ista_nrc, iend_nrc, nae, grad_tnae, &
                 &               wos, texc_ae, dF_dnae, dummy1, ecor, &
                 &               pot_type, ist, ien )
         endif

         call ex_gga_library( nspin, nspin, ista_nrc, iend_nrc, nps, grad_nps, &
              &               wos, texc_ps, dF_dnps, dF_dgradnps, dummy1, &
              &               pot_type, ist, ien )
         if ( xctype == 'ev93' .or. xctype == 'lb94' ) then
            grad_tnps = 0.d0
            call cr_gga_library( nspin, nspin, ista_nrc, iend_nrc, nps, grad_tnps, &
                 &               wos, texc_ps, dF_dnps, dummy1, ecor, &
                 &               pot_type, ist, ien )
         else if ( xctype == 'evpw91' ) then
            call cr_ggapw91( nspin, nspin, ista_nrc, iend_nrc, nps, grad_tnps, &
                 &           wos, texc_ps, dF_dnps, dummy1, ist, ien )
         else
            call cr_gga_library( nspin, nspin, ista_nrc, iend_nrc, nps, grad_tnps, &
                 &               wos, texc_ps, dF_dnps, dummy1, ecor, &
                 &               pot_type, ist, ien )
         endif

         exc_ae = exc_ae +texc_ae*omg_wght(ith)
         exc_ps = exc_ps +texc_ps*omg_wght(ith)
         deallocate( dummy1 )
#else
         call ex_gga_library( nspin, nspin, ista_nrc, iend_nrc, nae, grad_nae, &
              &               wos, texc_ae, dF_dnae, dF_dgradnae, &
              &               pot_type, ist, ien )
         if ( xctype == 'ev93' .or. xctype == 'lb94' ) then
            grad_tnae = 0.d0
            call cr_gga_library( nspin, nspin, ista_nrc, iend_nrc, nae, grad_tnae, &
                 &               wos, texc_ae, dF_dnae, ecor, &
                 &               pot_type, ist, ien )
         else if ( xctype == 'evpw91' ) then
            call cr_ggapw91( nspin, nspin, ista_nrc, iend_nrc, nae, grad_tnae, &
                 &           wos, texc_ae, dF_dnae, ist, ien )
         else
            call cr_gga_library( nspin, nspin, ista_nrc, iend_nrc, nae, grad_tnae, &
                 &               wos, texc_ae, dF_dnae, ecor, &
                 &               pot_type, ist, ien )
         endif

         call ex_gga_library( nspin, nspin, ista_nrc, iend_nrc, nps, grad_nps, &
              &               wos, texc_ps, dF_dnps, dF_dgradnps, &
              &               pot_type, ist, ien )
         if ( xctype == 'ev93' .or. xctype == 'lb94' ) then
            grad_tnps = 0.d0
            call cr_gga_library( nspin, nspin, ista_nrc, iend_nrc, nps, grad_tnps, &
                 &               wos, texc_ps, dF_dnps, ecor, &
                 &               pot_type, ist, ien )
         else if ( xctype == 'evpw91' ) then
            call cr_ggapw91( nspin, nspin, ista_nrc, iend_nrc, nps, grad_tnps, &
                 &           wos, texc_ps, dF_dnps, ist, ien )
         else
            call cr_gga_library( nspin, nspin, ista_nrc, iend_nrc, nps, grad_tnps, &
                 &               wos, texc_ps, dF_dnps, ecor, &
                 &               pot_type, ist, ien )
         endif

         exc_ae = exc_ae +texc_ae*omg_wght(ith)
         exc_ps = exc_ps +texc_ps*omg_wght(ith)
#endif

      else if ( xctype == 'vdwdf' ) then
         if ( exchange_pot_type == 'pbe' )  pot_type = 1
         if ( exchange_pot_type == 'revpbe' )  pot_type = 2
         if ( exchange_pot_type == 'b86r' )  pot_type = 11
         if ( exchange_pot_type == 'optpbe' )  pot_type = 12
         if ( exchange_pot_type == 'optb86b' )  pot_type = 13
         if ( exchange_pot_type == 'pw86r' )  pot_type = 14
         if ( exchange_pot_type == 'c09x' )  pot_type = 15
         if ( exchange_pot_type == 'lvpw86r' )  pot_type = 16

         allocate( tmp_work1(ista_nrc:iend_nrc) ); tmp_work1 = 0.0d0
#ifdef __EDA__
         allocate( dummy1(ista_nrc:iend_nrc) )
         call ex_gga_library( nspin, nspin, ista_nrc, iend_nrc, nae, grad_nae, &
              &               wos, texc_ae, dF_dnae, dF_dgradnae, &
              &               dummy1, pot_type, ist, ien )
         call cr_gga_library( nspin, nspin, ista_nrc, iend_nrc, nae, tmp_work1, &
              &               wos, texc_ae, dF_dnae, &
              &               dummy1, ecor, pot_type, ist, ien )
         call ex_gga_library( nspin, nspin, ista_nrc, iend_nrc, nps, grad_nps, &
              &               wos, texc_ps, dF_dnps, dF_dgradnps, &
              &               dummy1, pot_type, ist, ien )
         call cr_gga_library( nspin, nspin, ista_nrc, iend_nrc, nps, tmp_work1, &
              &               wos, texc_ps, dF_dnps, &
              &               dummy1, ecor, pot_type, ist, ien )
         exc_ae = exc_ae +texc_ae*omg_wght(ith)
         exc_ps = exc_ps +texc_ps*omg_wght(ith)
         deallocate( dummy1 )
#else
         call ex_gga_library( nspin, nspin, ista_nrc, iend_nrc, nae, grad_nae, &
              &               wos, texc_ae, dF_dnae, dF_dgradnae, &
              &               pot_type, ist, ien )
         call cr_gga_library( nspin, nspin, ista_nrc, iend_nrc, nae, tmp_work1, &
              &               wos, texc_ae, dF_dnae, ecor, &
              &               pot_type, ist, ien )
         call ex_gga_library( nspin, nspin, ista_nrc, iend_nrc, nps, grad_nps, &
              &               wos, texc_ps, dF_dnps, dF_dgradnps, &
              &               pot_type, ist, ien )
         call cr_gga_library( nspin, nspin, ista_nrc, iend_nrc, nps, tmp_work1, &
              &               wos, texc_ps, dF_dnps, ecor, &
              &               pot_type, ist, ien )
         exc_ae = exc_ae +texc_ae*omg_wght(ith)
         exc_ps = exc_ps +texc_ps*omg_wght(ith)
#endif
         deallocate( tmp_work1 )

#ifdef LIBXC
      else if ( xctype == 'libxc' ) then
         select case (xc_family_exch)
         case (XC_FAMILY_LDA)
#ifdef __EDA__
            allocate( dummy1(ista_nrc:iend_nrc) )
            call ex_lda_libxc( nspin, nspin, ista_nrc, iend_nrc, &
                 &             nae, wos, texc_ae, dF_dnae, dummy1, ist, ien )
            call ex_lda_libxc( nspin, nspin, ista_nrc, iend_nrc, &
                 &             nps, wos, texc_ps, dF_dnps, dummy1, ist, ien )
            deallocate( dummy1 )
#else
            call ex_lda_libxc( nspin, nspin, ista_nrc, iend_nrc, &
                 &             nae, wos, texc_ae, dF_dnae, ist, ien )
            call ex_lda_libxc( nspin, nspin, ista_nrc, iend_nrc, &
                 &             nps, wos, texc_ps, dF_dnps, ist, ien )
#endif
         case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
#ifdef __EDA__
            allocate( dummy1(ista_nrc:iend_nrc) )
            call ex_gga_libxc( nspin, nspin, ista_nrc, iend_nrc, &
                 &             nae, grad_nae, grad_tnae, wos, &
                 &             texc_ae, dF_dnae, dF_dgradnae, dummy1, ist, ien )
            call ex_gga_libxc( nspin, nspin, ista_nrc, iend_nrc, &
                 &             nps, grad_nps, grad_tnps, wos, &
                 &             texc_ps, dF_dnps, dF_dgradnps, dummy1, ist, ien )

            deallocate( dummy1 )
#else
            call ex_gga_libxc( nspin, nspin, ista_nrc, iend_nrc, &
                 &             nae, grad_nae, grad_tnae, wos, &
                 &             texc_ae, dF_dnae, dF_dgradnae, ist, ien )
            call ex_gga_libxc( nspin, nspin, ista_nrc, iend_nrc, &
                 &             nps, grad_nps, grad_tnps, wos, &
                 &             texc_ps, dF_dnps, dF_dgradnps, ist, ien )
#endif
         end select

         select case (xc_family_corr)
         case (XC_FAMILY_LDA)
#ifdef __EDA__
            allocate( dummy1(ista_nrc:iend_nrc) )
            call cr_lda_libxc( nspin, nspin, ista_nrc, iend_nrc, &
                 &             nae, wos, texc_ae, dF_dnae, dummy1, ist, ien )
            call cr_lda_libxc( nspin, nspin, ista_nrc, iend_nrc, &
                 &             nps, wos, texc_ps, dF_dnps, dummy1, ist, ien )
            deallocate( dummy1 )
#else
            call cr_lda_libxc( nspin, nspin, ista_nrc, iend_nrc, &
                 &             nae, wos, texc_ae, dF_dnae, ist, ien )
            call cr_lda_libxc( nspin, nspin, ista_nrc, iend_nrc, &
                 &             nps, wos, texc_ps, dF_dnps, ist, ien )
#endif
         case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
            allocate( tmp_work1(ista_nrc:iend_nrc) )
            allocate( tmp_work2(ista_nrc:iend_nrc) )
#ifdef __EDA__
            allocate( dummy1(ista_nrc:iend_nrc) )
            call cr_gga_libxc( nspin, nspin, ista_nrc, iend_nrc, &
                 &             nae, grad_nae, grad_tnae, wos, &
                 &             texc_ae, dF_dnae, tmp_work1, dummy1, ist, ien )
            call cr_gga_libxc( nspin, nspin, ista_nrc, iend_nrc, &
                 &             nps, grad_nps, grad_tnps, wos, &
                 &             texc_ps, dF_dnps, tmp_work2, dummy1, ist, ien )
            deallocate( dummy1 )
#else
            call cr_gga_libxc( nspin, nspin, ista_nrc, iend_nrc, &
                 &             nae, grad_nae, grad_tnae, wos, &
                 &             texc_ae, dF_dnae, tmp_work1, ist, ien )
            call cr_gga_libxc( nspin, nspin, ista_nrc, iend_nrc, &
                 &             nps, grad_nps, grad_tnps, wos, &
                 &             texc_ps, dF_dnps, tmp_work2, ist, ien )
#endif
            grad_tnae = tmp_work1;    grad_tnps = tmp_work2
            deallocate( tmp_work1 );    deallocate( tmp_work2 )
         end select

         exc_ae = exc_ae +texc_ae *omg_wght(ith)
         exc_ps = exc_ps +texc_ps *omg_wght(ith)
#endif
      else
         write(nfout,'(" xctype = ",a7)') xctype
         call phase_error_with_msg(nfout,' xctype is not set properly (ggaxcp0_paw)',__LINE__,__FILE__)
      end if
             !   dF/d|rho(r)| (=vxc) --> dF_drho
             !   dFx/d|grad(rho(r))| --> dF_dgradrho
             !   dFc/d|grad(rho(r))| --> grad_trho
1000  continue

      call tstatc0_end(id_sname)
    end subroutine ggaxcp0_paw

    subroutine ggaxcp_paw
      call ggaxcp0_paw

      return
    end subroutine ggaxcp_paw

    subroutine abs_grad_rho_up_down_total_paw_3D
      integer:: ir,is,iord,dnr
      real(kind=DP), allocatable, dimension(:,:) :: tmp ! === For nrc decomposion. by T. Yamasaki 2020/02/24 ===
      integer       :: id_sname = -1,ier, i, isdiff
!!$#if 0
!!$      real(kind=DP), allocatable, dimension(:,:) :: t_d_dr, t_dd_ddr
!!$#endif

      call tstatc0_begin('abs_grad_rho_up_down_total_paw_3D ',id_sname)
!!$#if 0
!!$      write(nfout,'(" -- abs_grad_rho_up_down_total_paw_3D, nrc = ",i4," --")') nrc
!!$      call flush(nfout)
!!$
!!$      if(.not.allocated(nae_g)) allocate(nae_g(mmesh,nspin)); nae_g = 0.d0
!!$      if(.not.allocated(dnae_dr_g)) allocate(dnae_dr_g(mmesh,nspin)); dnae_dr_g = 0.d0
!!$      if(.not.allocated(ddnae_ddr_g)) allocate(ddnae_ddr_g(mmesh,nspin)); ddnae_ddr_g = 0.d0
!!$      do is=1,nspin
!!$         do i = max(1,ista_nrc),min(mmesh,iend_nrc)
!!$            nae_g(i,is) = nae(i,is)
!!$         end do
!!$      end do
!!$      call mpi_allreduce(MPI_IN_PLACE,nae_g,mmesh*nspin,mpi_double_precision,mpi_sum,mpi_natm_world,ier)
!!$      iord = paw_density_gradient%order
!!$! ------ for AE CD -------
!!$      if(vflag == VXC_AND_EXC) then
!!$         do is=1,nspin
!!$            call calc_diff_exp2(2,iord,nrc,xh(it),radr_paw(:,it), &
!!$                 & nae_g(:,is), dnae_dr_g(:,is), ddnae_ddr_g(:,is))  ! nae(:,is), dnae_dr(:,is), ddnae_ddr(:,is))
!!$         end do
!!$      else
!!$         do is=1,nspin
!!$            call calc_diff_exp2(1,iord,nrc,xh(it),radr_paw(:,it), &
!!$                 & nae_g(:,is), dnae_dr_g(:,is), ddnae_ddr_g(:,is))  ! nae(:,is), dnae_dr(:,is), ddnae_ddr(:,is))
!!$         end do
!!$      end if
!!$      dnae_dr(ista_nrc:iend_nrc,1:nspin) = dnae_dr_g(ista_nrc:iend_nrc,1:nspin)
!!$      ddnae_ddr(ista_nrc:iend_nrc,1:nspin) = ddnae_ddr_g(ista_nrc:iend_nrc,1:nspin)
!!$
!!!!$      deallocate(tmp)
!!!!$      deallocate(t_d_dr)
!!!!$      deallocate(t_dd_ddr)
!!$#else
      dnr = 1
      allocate(tmp(ista_nrc-6*dnr:iend_nrc+6*dnr,nspin))   ! == For nrc decomposion. by T. Yamasaki 2020/02/24 ==
      call boundary_exchange_dim2(tmp, nae, dnr, nspin) ! nae -> tmp


      iord = paw_density_gradient%order
! ------ for AE CD -------
      isdiff = 1
      if(vflag == VXC_AND_EXC) isdiff = 2
      do is=1,nspin
         call calc_diff_exp2_3D(isdiff,iord,nrc,dnr,xh(it),radr_paw(:,it),tmp(:,is) &
              &             ,dnae_dr(:,is),ista_nrc,iend_nrc, ddnae_ddr(:,is))
      end do
!!$#endif

      if(flg_symmetry) then
         call m_PAWCD_set_ae_der_cd_sym(dvec,ith,iph,ia,nspin,nrc,dnae_dth,dnae_dph,ista_nrc,iend_nrc)
      else
         call m_PAWCD_set_ae_der_cd(ith,iph,ia,nspin,nrc, dnae_dth,dnae_dph,ista_nrc,iend_nrc)
      end if

!!$      dnae_dth = 0d0; dnae_dph = 0.d0
      if(iprixc >= PRINTLEVEL) then
         write(nfout,'(" size of dnae_dth = [",i8,":",i8,",",i8,":",i8,"]")') &
              &                              lbound(dnae_dth,1),ubound(dnae_dth,1) &
              &                             ,lbound(dnae_dth,2),ubound(dnae_dth,2)
         write(nfout,'(" dvec = ",3f8.4)') dvec(1),dvec(2),dvec(3)
         write(nfout,'(" --- dnae_dth --- ir=[",i4,":",i4,"]")') ista_nrc,ista_nrc+11
         write(nfout,'(6f16.8)') (dnae_dth(ir,1), ir=ista_nrc,ista_nrc+11)
         write(nfout,'(6f16.8)') (dnae_dth(ir,2), ir=ista_nrc,ista_nrc+11)
         write(nfout,'(" --- dnae_dph ---")')
         write(nfout,'(6f16.8)') (dnae_dph(ir,1), ir=ista_nrc,ista_nrc+11)
         write(nfout,'(6f16.8)') (dnae_dph(ir,2), ir=ista_nrc,ista_nrc+11)
      end if
      do is=1,nspin
         do ir=ista_nrc,min(nrc,iend_nrc)
            grad_nae(ir,is) = dsqrt(  dnae_dr(ir,is)*dnae_dr(ir,is) &
                 &                  + dnae_dth(ir,is)*dnae_dth(ir,is) &
                 &                  + dnae_dph(ir,is)*dnae_dph(ir,is))
         end do
      end do
      if(nspin==2) then
         do ir=ista_nrc,min(nrc,iend_nrc)
            grad_tnae(ir) = dsqrt(  (dnae_dr(ir,1)+dnae_dr(ir,2))**2 &
                 &                + (dnae_dth(ir,1)+dnae_dth(ir,2))**2 &
                 &                + (dnae_dph(ir,1)+dnae_dph(ir,2))**2)
         end do
      else
         grad_tnae(:)=grad_nae(:,1)
      end if

!!$#if 0
!!$      allocate(tmp(mmesh,nspin)); tmp = 0.d0
!!$      allocate(t_d_dr(mmesh,nspin)); t_d_dr = 0.d0
!!$      allocate(t_dd_ddr(mmesh,nspin)); t_dd_ddr = 0.d0
!!$      do is=1,nspin
!!$         tmp(ista_nrc:iend_nrc,is) = nps(ista_nrc:iend_nrc,is)
!!$      end do
!!$      call mpi_allreduce(MPI_IN_PLACE,tmp,mmesh*nspin,mpi_double_precision,mpi_sum,mpi_natm_world,ier)
!!$      iord = paw_density_gradient%order
!!$! ------ for PS CD -------
!!$      if(vflag == VXC_AND_EXC) then
!!$         do is=1,nspin
!!$            call calc_diff_exp2(2,iord,nrc,xh(it),radr_paw(:,it), &
!!$                 & tmp(:,is), t_d_dr(:,is), t_dd_ddr(:,is))  ! nae(:,is), dnae_dr(:,is), ddnae_ddr(:,is))
!!$         end do
!!$      else
!!$         do is=1,nspin
!!$            call calc_diff_exp2(1,iord,nrc,xh(it),radr_paw(:,it), &
!!$                 & tmp(:,is), t_d_dr(:,is), t_dd_ddr(:,is))  ! nae(:,is), dnae_dr(:,is), ddnae_ddr(:,is))
!!$         end do
!!$      end if
!!$      dnps_dr(ista_nrc:iend_nrc,1:nspin) = t_d_dr(ista_nrc:iend_nrc,1:nspin)
!!$      ddnps_ddr(ista_nrc:iend_nrc,1:nspin) = t_dd_ddr(ista_nrc:iend_nrc,1:nspin)
!!$
!!$      deallocate(tmp)
!!$      deallocate(t_d_dr)
!!$      deallocate(t_dd_ddr)
!!$#else
! ------ for PS CD -------
      dnr = 1
      if(.not.allocated(tmp)) &
           & allocate(tmp(ista_nrc-6*dnr:iend_nrc+6*dnr,nspin))   ! == For nrc decomposion. by T. Yamasaki 2020/02/24 ==
      call boundary_exchange_dim2(tmp, nps, dnr, nspin) ! nps -> tmp

      isdiff = 1
      if(vflag == VXC_AND_EXC) isdiff = 2
      do is=1,nspin
         call calc_diff_exp2_3D(isdiff,iord,nrc,dnr,xh(it),radr_paw(:,it), tmp(:,is) &
              &             ,dnps_dr(:,is),ista_nrc,iend_nrc,ddnps_ddr(:,is))
      end do
      deallocate(tmp)
!!$#endif
      if(flg_symmetry) then
         call m_PAWCD_set_ps_der_cd_sym(dvec,ith,iph,ia,nspin,nrc,dnps_dth,dnps_dph,ista_nrc,iend_nrc)
      else
         call m_PAWCD_set_ps_der_cd(ith,iph,ia,nspin,nrc,dnps_dth,dnps_dph,ista_nrc,iend_nrc)
      end if
      do is=1,nspin
         do ir=ista_nrc,min(nrc,iend_nrc)   !!$         do ir=1,nrc
            grad_nps(ir,is)=dsqrt( dnps_dr(ir,is)*dnps_dr(ir,is)  &
                 &               + dnps_dth(ir,is)*dnps_dth(ir,is) &
                 &               + dnps_dph(ir,is)*dnps_dph(ir,is) )
         end do
      end do
      if(nspin==2) then
         do ir=ista_nrc,min(nrc,iend_nrc)   !!$         do ir=1,nrc
            grad_tnps(ir)=dsqrt(  (dnps_dr(ir,1)+dnps_dr(ir,2))**2 &
                 &              + (dnps_dth(ir,1)+dnps_dth(ir,2))**2 &
                 &              + (dnps_dph(ir,1)+dnps_dph(ir,2))**2 )
         end do
      else
         grad_tnps(:)=grad_nps(:,1)  !!$       grad_tnps(1:nrc)=grad_nps(1:nrc,1)
      end if

      call tstatc0_end(id_sname)

    end subroutine abs_grad_rho_up_down_total_paw_3D

  end subroutine m_PAW_XC_cal_potential
! === For nrc decomposion. by takto 2012/12/05 =================================
    subroutine paw_xc_sphex_allocate_3D
        allocate(nae_sph(ista_nrc:iend_nrc,nspin,25))  ;nae_sph=0.d0
        allocate(nps_sph(ista_nrc:iend_nrc,nspin,25))  ;nps_sph=0.d0
        allocate(wos(ista_nrc:iend_nrc));wos=0.d0
        allocate(irs(ista_nrc:iend_nrc));irs=0
        allocate(exc_ae_field(ista_nrc:iend_nrc))      ;exc_ae_field=0.d0
        allocate(exc_ps_field(ista_nrc:iend_nrc))      ;exc_ps_field=0.d0

        allocate(nana_ae_sph(ista_nrc:iend_nrc,25))    ;nana_ae_sph=0.d0
        allocate(nana_ps_sph(ista_nrc:iend_nrc,25))    ;nana_ps_sph=0.d0

        allocate(dFxcdna_ae_sph(ista_nrc:iend_nrc,25)) ;dFxcdna_ae_sph=0.d0
        allocate(dFxcdna_ps_sph(ista_nrc:iend_nrc,25)) ;dFxcdna_ps_sph=0.d0

            allocate(nanb_ae_sph(ista_nrc:iend_nrc,25))    ;nanb_ae_sph=0.d0
            allocate(nanb_ps_sph(ista_nrc:iend_nrc,25))    ;nanb_ps_sph=0.d0

        if(nspin.eq.2) then
            allocate(nbnb_ae_sph(ista_nrc:iend_nrc,25))    ;nbnb_ae_sph=0.d0
            allocate(nbnb_ps_sph(ista_nrc:iend_nrc,25))    ;nbnb_ps_sph=0.d0

            allocate(dFxcdnb_ae_sph(ista_nrc:iend_nrc,25)) ;dFxcdnb_ae_sph=0.d0
            allocate(dFxcdnb_ps_sph(ista_nrc:iend_nrc,25)) ;dFxcdnb_ps_sph=0.d0
        end if

        if(check_of_xctype()==GGA) then
            allocate(grad_nae2_sph(ista_nrc:iend_nrc,nspin,25))        ;grad_nae2_sph=0.d0
            allocate(grad_tnae2_sph(ista_nrc:iend_nrc,25))             ;grad_tnae2_sph=0.d0
            allocate(grad_nps2_sph(ista_nrc:iend_nrc,nspin,25))        ;grad_nps2_sph=0.d0
            allocate(grad_tnps2_sph(ista_nrc:iend_nrc,25))             ;grad_tnps2_sph=0.d0
            allocate(dnae_dr_sph(ista_nrc:iend_nrc,nspin,25))          ;dnae_dr_sph=0.d0
            allocate(ddnae_ddr_sph(ista_nrc:iend_nrc,nspin,25))        ;ddnae_ddr_sph=0.d0
            allocate(dnps_dr_sph(ista_nrc:iend_nrc,nspin,25))          ;dnps_dr_sph=0.d0
            allocate(ddnps_ddr_sph(ista_nrc:iend_nrc,nspin,25))        ;ddnps_ddr_sph=0.d0
            allocate(grad_nae(ista_nrc:iend_nrc,nspin))                ;grad_nae=0.d0
            allocate(grad_tnae(ista_nrc:iend_nrc))                     ;grad_tnae=0.d0
            allocate(grad_nps(ista_nrc:iend_nrc,nspin))                ;grad_nps=0.d0
            allocate(grad_tnps(ista_nrc:iend_nrc))                     ;grad_tnps=0.d0
            allocate(dF_dnae(ista_nrc:iend_nrc,nspin))                 ;dF_dnae=0.d0
            allocate(dF_dgradnae(ista_nrc:iend_nrc,nspin))             ;dF_dgradnae=0.d0
            allocate(dF_dnps(ista_nrc:iend_nrc,nspin))                 ;dF_dnps=0.d0
            allocate(dF_dgradnps(ista_nrc:iend_nrc,nspin))             ;dF_dgradnps=0.d0
            allocate(dFx_dnnae(ista_nrc:iend_nrc,nspin))               ;dFx_dnnae=0.d0
            allocate(dFx_dngae(ista_nrc:iend_nrc,nspin))               ;dFx_dngae=0.d0
            allocate(dFx_dggae(ista_nrc:iend_nrc,nspin))               ;dFx_dggae=0.d0
            allocate(dFx_dnnnae(ista_nrc:iend_nrc,nspin))              ;dFx_dnnnae=0.d0
            allocate(dFx_dnngae(ista_nrc:iend_nrc,nspin))              ;dFx_dnngae=0.d0
            allocate(dFx_dnggae(ista_nrc:iend_nrc,nspin))              ;dFx_dnggae=0.d0
            allocate(dFx_dgggae(ista_nrc:iend_nrc,nspin))              ;dFx_dgggae=0.d0
            allocate(dFx_dnnps(ista_nrc:iend_nrc,nspin))               ;dFx_dnnps=0.d0
            allocate(dFx_dngps(ista_nrc:iend_nrc,nspin))               ;dFx_dngps=0.d0
            allocate(dFx_dggps(ista_nrc:iend_nrc,nspin))               ;dFx_dggps=0.d0
            allocate(dFx_dnnnps(ista_nrc:iend_nrc,nspin))              ;dFx_dnnnps=0.d0
            allocate(dFx_dnngps(ista_nrc:iend_nrc,nspin))              ;dFx_dnngps=0.d0
            allocate(dFx_dnggps(ista_nrc:iend_nrc,nspin))              ;dFx_dnggps=0.d0
            allocate(dFx_dgggps(ista_nrc:iend_nrc,nspin))              ;dFx_dgggps=0.d0

            allocate(dF_dgradtnae(ista_nrc:iend_nrc))                  ;dF_dgradtnae=0.d0
            allocate(dFc_daa_ae(ista_nrc:iend_nrc))                    ;dFc_daa_ae=0.d0
            allocate(dFc_dbb_ae(ista_nrc:iend_nrc))                    ;dFc_dbb_ae=0.d0
            allocate(dFc_dgg_ae(ista_nrc:iend_nrc))                    ;dFc_dgg_ae=0.d0
            allocate(dFc_dab_ae(ista_nrc:iend_nrc))                    ;dFc_dab_ae=0.d0
            allocate(dFc_dag_ae(ista_nrc:iend_nrc))                    ;dFc_dag_ae=0.d0
            allocate(dFc_dbg_ae(ista_nrc:iend_nrc))                    ;dFc_dbg_ae=0.d0
            allocate(dFc_daaa_ae(ista_nrc:iend_nrc))                   ;dFc_daaa_ae=0.d0
            allocate(dFc_dbbb_ae(ista_nrc:iend_nrc))                   ;dFc_dbbb_ae=0.d0
            allocate(dFc_dggg_ae(ista_nrc:iend_nrc))                   ;dFc_dggg_ae=0.d0
            allocate(dFc_daab_ae(ista_nrc:iend_nrc))                   ;dFc_daab_ae=0.d0
            allocate(dFc_daag_ae(ista_nrc:iend_nrc))                   ;dFc_daag_ae=0.d0
            allocate(dFc_dabb_ae(ista_nrc:iend_nrc))                   ;dFc_dabb_ae=0.d0
            allocate(dFc_dbbg_ae(ista_nrc:iend_nrc))                   ;dFc_dbbg_ae=0.d0
            allocate(dFc_dagg_ae(ista_nrc:iend_nrc))                   ;dFc_dagg_ae=0.d0
            allocate(dFc_dbgg_ae(ista_nrc:iend_nrc))                   ;dFc_dbgg_ae=0.d0
            allocate(dFc_dabg_ae(ista_nrc:iend_nrc))                   ;dFc_dabg_ae=0.d0

            allocate(dF_dgradtnps(ista_nrc:iend_nrc))                  ;dF_dgradtnps=0.d0
            allocate(dFc_daa_ps(ista_nrc:iend_nrc))                    ;dFc_daa_ps=0.d0
            allocate(dFc_dbb_ps(ista_nrc:iend_nrc))                    ;dFc_dbb_ps=0.d0
            allocate(dFc_dgg_ps(ista_nrc:iend_nrc))                    ;dFc_dgg_ps=0.d0
            allocate(dFc_dab_ps(ista_nrc:iend_nrc))                    ;dFc_dab_ps=0.d0
            allocate(dFc_dag_ps(ista_nrc:iend_nrc))                    ;dFc_dag_ps=0.d0
            allocate(dFc_dbg_ps(ista_nrc:iend_nrc))                    ;dFc_dbg_ps=0.d0
            allocate(dFc_daaa_ps(ista_nrc:iend_nrc))                   ;dFc_daaa_ps=0.d0
            allocate(dFc_dbbb_ps(ista_nrc:iend_nrc))                   ;dFc_dbbb_ps=0.d0
            allocate(dFc_dggg_ps(ista_nrc:iend_nrc))                   ;dFc_dggg_ps=0.d0
            allocate(dFc_daab_ps(ista_nrc:iend_nrc))                   ;dFc_daab_ps=0.d0
            allocate(dFc_daag_ps(ista_nrc:iend_nrc))                   ;dFc_daag_ps=0.d0
            allocate(dFc_dabb_ps(ista_nrc:iend_nrc))                   ;dFc_dabb_ps=0.d0
            allocate(dFc_dbbg_ps(ista_nrc:iend_nrc))                   ;dFc_dbbg_ps=0.d0
            allocate(dFc_dagg_ps(ista_nrc:iend_nrc))                   ;dFc_dagg_ps=0.d0
            allocate(dFc_dbgg_ps(ista_nrc:iend_nrc))                   ;dFc_dbgg_ps=0.d0
            allocate(dFc_dabg_ps(ista_nrc:iend_nrc))                   ;dFc_dabg_ps=0.d0

            allocate(gaga_ae_sph(ista_nrc:iend_nrc,25))                ;gaga_ae_sph=0.d0
            allocate(naga_ae_sph(ista_nrc:iend_nrc,25))                ;naga_ae_sph=0.d0
            allocate(gg_ae_sph(ista_nrc:iend_nrc,25))                  ;gg_ae_sph=0.d0
            allocate(nag_ae_sph(ista_nrc:iend_nrc,25))                 ;nag_ae_sph=0.d0

            allocate(gaga_ps_sph(ista_nrc:iend_nrc,25))                ;gaga_ps_sph=0.d0
            allocate(naga_ps_sph(ista_nrc:iend_nrc,25))                ;naga_ps_sph=0.d0
            allocate(gg_ps_sph(ista_nrc:iend_nrc,25))                  ;gg_ps_sph=0.d0
            allocate(nag_ps_sph(ista_nrc:iend_nrc,25))                 ;nag_ps_sph=0.d0

            allocate(dFxdgaovrga_ae_sph(ista_nrc:iend_nrc,25))         ;dFxdgaovrga_ae_sph=0.d0
            allocate(dFxdgaovrga_ps_sph(ista_nrc:iend_nrc,25))         ;dFxdgaovrga_ps_sph=0.d0
            allocate(dFcdgovrg_ae_sph(ista_nrc:iend_nrc,25))           ;dFcdgovrg_ae_sph=0.d0
            allocate(dFcdgovrg_ps_sph(ista_nrc:iend_nrc,25))           ;dFcdgovrg_ps_sph=0.d0

            allocate(dFadga_ae(ista_nrc:iend_nrc))                     ;dFadga_ae=0.d0
            allocate(dFadga_ps(ista_nrc:iend_nrc))                     ;dFadga_ps=0.d0
            allocate(dFadg_ae(ista_nrc:iend_nrc))                      ;dFadg_ae=0.d0
            allocate(dFadg_ps(ista_nrc:iend_nrc))                      ;dFadg_ps=0.d0
            allocate(dFadgaga_ae(ista_nrc:iend_nrc))                   ;dFadgaga_ae=0.d0
            allocate(dFadgaga_ps(ista_nrc:iend_nrc))                   ;dFadgaga_ps=0.d0
            allocate(dFadgg_ae(ista_nrc:iend_nrc))                     ;dFadgg_ae=0.d0
            allocate(dFadgg_ps(ista_nrc:iend_nrc))                     ;dFadgg_ps=0.d0
            allocate(dFadnaga_ae(ista_nrc:iend_nrc))                   ;dFadnaga_ae=0.d0
            allocate(dFadnaga_ps(ista_nrc:iend_nrc))                   ;dFadnaga_ps=0.d0
            allocate(dFadnag_ae(ista_nrc:iend_nrc))                    ;dFadnag_ae=0.d0
            allocate(dFadnag_ps(ista_nrc:iend_nrc))                    ;dFadnag_ps=0.d0

            allocate(dGadna_ae(ista_nrc:iend_nrc))                     ;dGadna_ae=0.d0
            allocate(dGadna_ps(ista_nrc:iend_nrc))                     ;dGadna_ps=0.d0
            allocate(dGadga_ae(ista_nrc:iend_nrc))                     ;dGadga_ae=0.d0
            allocate(dGadga_ps(ista_nrc:iend_nrc))                     ;dGadga_ps=0.d0
            allocate(dGadnana_ae(ista_nrc:iend_nrc))                   ;dGadnana_ae=0.d0
            allocate(dGadnana_ps(ista_nrc:iend_nrc))                   ;dGadnana_ps=0.d0
            allocate(dGadnaga_ae(ista_nrc:iend_nrc))                   ;dGadnaga_ae=0.d0
            allocate(dGadnaga_ps(ista_nrc:iend_nrc))                   ;dGadnaga_ps=0.d0
            allocate(dGadgaga_ae(ista_nrc:iend_nrc))                   ;dGadgaga_ae=0.d0
            allocate(dGadgaga_ps(ista_nrc:iend_nrc))                   ;dGadgaga_ps=0.d0

            allocate(dGdna_ae(ista_nrc:iend_nrc))                      ;dGdna_ae=0.d0
            allocate(dGdna_ps(ista_nrc:iend_nrc))                      ;dGdna_ps=0.d0
            allocate(dGdg_ae(ista_nrc:iend_nrc))                       ;dGdg_ae=0.d0
            allocate(dGdg_ps(ista_nrc:iend_nrc))                       ;dGdg_ps=0.d0
            allocate(dGdnana_ae(ista_nrc:iend_nrc))                    ;dGdnana_ae=0.d0
            allocate(dGdnana_ps(ista_nrc:iend_nrc))                    ;dGdnana_ps=0.d0
            allocate(dGdgg_ae(ista_nrc:iend_nrc))                      ;dGdgg_ae=0.d0
            allocate(dGdgg_ps(ista_nrc:iend_nrc))                      ;dGdgg_ps=0.d0
            allocate(dGdnag_ae(ista_nrc:iend_nrc))                     ;dGdnag_ae=0.d0
            allocate(dGdnag_ps(ista_nrc:iend_nrc))                     ;dGdnag_ps=0.d0

            if(nspin.eq.2) then
                allocate(gbgb_ae_sph(ista_nrc:iend_nrc,25))            ;gbgb_ae_sph=0.d0
                allocate(nbgb_ae_sph(ista_nrc:iend_nrc,25))            ;nbgb_ae_sph=0.d0
                allocate(nbg_ae_sph(ista_nrc:iend_nrc,25))             ;nbg_ae_sph=0.d0

                allocate(gbgb_ps_sph(ista_nrc:iend_nrc,25))            ;gbgb_ps_sph=0.d0
                allocate(nbgb_ps_sph(ista_nrc:iend_nrc,25))            ;nbgb_ps_sph=0.d0
                allocate(nbg_ps_sph(ista_nrc:iend_nrc,25))             ;nbg_ps_sph=0.d0

                allocate(dFxdgbovrgb_ae_sph(ista_nrc:iend_nrc,25))     ;dFxdgbovrgb_ae_sph=0.d0
                allocate(dFxdgbovrgb_ps_sph(ista_nrc:iend_nrc,25))     ;dFxdgbovrgb_ps_sph=0.d0

                allocate(dFbdgb_ae(ista_nrc:iend_nrc))                 ;dFbdgb_ae=0.d0
                allocate(dFbdgb_ps(ista_nrc:iend_nrc))                 ;dFbdgb_ps=0.d0
                allocate(dFbdg_ae(ista_nrc:iend_nrc))                  ;dFbdg_ae=0.d0
                allocate(dFbdg_ps(ista_nrc:iend_nrc))                  ;dFbdg_ps=0.d0
                allocate(dFbdgbgb_ae(ista_nrc:iend_nrc))               ;dFbdgbgb_ae=0.d0
                allocate(dFbdgbgb_ps(ista_nrc:iend_nrc))               ;dFbdgbgb_ps=0.d0
                allocate(dFbdgg_ae(ista_nrc:iend_nrc))                 ;dFbdgg_ae=0.d0
                allocate(dFbdgg_ps(ista_nrc:iend_nrc))                 ;dFbdgg_ps=0.d0
                allocate(dFbdnbgb_ae(ista_nrc:iend_nrc))               ;dFbdnbgb_ae=0.d0
                allocate(dFbdnbgb_ps(ista_nrc:iend_nrc))               ;dFbdnbgb_ps=0.d0
                allocate(dFbdnbg_ae(ista_nrc:iend_nrc))                ;dFbdnbg_ae=0.d0
                allocate(dFbdnbg_ps(ista_nrc:iend_nrc))                ;dFbdnbg_ps=0.d0
                allocate(dFbdag_ae(ista_nrc:iend_nrc))                 ;dFbdag_ae=0.d0
                allocate(dFbdag_ps(ista_nrc:iend_nrc))                 ;dFbdag_ps=0.d0

                allocate(dGbdnb_ae(ista_nrc:iend_nrc))                 ;dGbdnb_ae=0.d0
                allocate(dGbdnb_ps(ista_nrc:iend_nrc))                 ;dGbdnb_ps=0.d0
                allocate(dGbdgb_ae(ista_nrc:iend_nrc))                 ;dGbdgb_ae=0.d0
                allocate(dGbdgb_ps(ista_nrc:iend_nrc))                 ;dGbdgb_ps=0.d0
                allocate(dGbdnbnb_ae(ista_nrc:iend_nrc))               ;dGbdnbnb_ae=0.d0
                allocate(dGbdnbnb_ps(ista_nrc:iend_nrc))               ;dGbdnbnb_ps=0.d0
                allocate(dGbdnbgb_ae(ista_nrc:iend_nrc))               ;dGbdnbgb_ae=0.d0
                allocate(dGbdnbgb_ps(ista_nrc:iend_nrc))               ;dGbdnbgb_ps=0.d0
                allocate(dGbdgbgb_ae(ista_nrc:iend_nrc))               ;dGbdgbgb_ae=0.d0
                allocate(dGbdgbgb_ps(ista_nrc:iend_nrc))               ;dGbdgbgb_ps=0.d0

                allocate(dGdnb_ae(ista_nrc:iend_nrc))                  ;dGdnb_ae=0.d0
                allocate(dGdnb_ps(ista_nrc:iend_nrc))                  ;dGdnb_ps=0.d0
                allocate(dGdnbnb_ae(ista_nrc:iend_nrc))                ;dGdnbnb_ae=0.d0
                allocate(dGdnbnb_ps(ista_nrc:iend_nrc))                ;dGdnbnb_ps=0.d0
                allocate(dGdnanb_ae(ista_nrc:iend_nrc))                ;dGdnanb_ae=0.d0
                allocate(dGdnanb_ps(ista_nrc:iend_nrc))                ;dGdnanb_ps=0.d0
                allocate(dGdnbg_ae(ista_nrc:iend_nrc))                 ;dGdnbg_ae=0.d0
                allocate(dGdnbg_ps(ista_nrc:iend_nrc))                 ;dGdnbg_ps=0.d0

            end if


!            allocate(dF_dnae(mmesh,nspin))
!            allocate(dF_dgradnae(mmesh,nspin))
!            allocate(dF_dgradnae_dr(mmesh,nspin))
!            allocate(dgrad_tnps_dr(mmesh))
!            allocate(dF_dnps(mmesh,nspin))
!            allocate(dF_dgradnps(mmesh,nspin))
!            allocate(dF_dgradnps_dr(mmesh,nspin))
!            allocate(dnae_dr(mmesh,nspin))
!            allocate(dnae_dth(mmesh,nspin))
!            allocate(dnae_dph(mmesh,nspin))
!            allocate(ddnae_ddr(mmesh,nspin))
!            allocate(dnps_dr(mmesh,nspin))
!            allocate(dnps_dth(mmesh,nspin))
!            allocate(dnps_dph(mmesh,nspin))
!            allocate(ddnps_ddr(mmesh,nspin))
!            allocate(dgrad_tnae_dr(mmesh))
        end if

        if(af /= 0) then
            allocate(flg_done(natm))
            flg_done=.false.
        end if
        return
    end subroutine paw_xc_sphex_allocate_3D
! ==============================================================================

    subroutine paw_xc_sphex_deallocate
        deallocate(nae_sph)
        deallocate(nps_sph)
        deallocate(wos)
        deallocate(irs)
        deallocate(exc_ae_field)
        deallocate(exc_ps_field)

        deallocate(nana_ae_sph)
        deallocate(nana_ps_sph)

        deallocate(dFxcdna_ae_sph)
        deallocate(dFxcdna_ps_sph)

            deallocate(nanb_ae_sph)
            deallocate(nanb_ps_sph)

        if(nspin.eq.2) then
            deallocate(nbnb_ae_sph)
            deallocate(nbnb_ps_sph)

            deallocate(dFxcdnb_ae_sph)
            deallocate(dFxcdnb_ps_sph)
        end if

        if(check_of_xctype()==GGA) then
            deallocate(grad_nae2_sph)
            deallocate(grad_tnae2_sph)
            deallocate(grad_nps2_sph)
            deallocate(grad_tnps2_sph)
            deallocate(dnae_dr_sph)
            deallocate(ddnae_ddr_sph)
            deallocate(dnps_dr_sph)
            deallocate(ddnps_ddr_sph)
            deallocate(grad_nae)
            deallocate(grad_tnae)
            deallocate(grad_nps)
            deallocate(grad_tnps)
            deallocate(dF_dnae)
            deallocate(dF_dgradnae)
            deallocate(dF_dnps)
            deallocate(dF_dgradnps)
            deallocate(dFx_dnnae)
            deallocate(dFx_dngae)
            deallocate(dFx_dggae)
            deallocate(dFx_dnnnae)
            deallocate(dFx_dnngae)
            deallocate(dFx_dnggae)
            deallocate(dFx_dgggae)
            deallocate(dFx_dnnps)
            deallocate(dFx_dngps)
            deallocate(dFx_dggps)
            deallocate(dFx_dnnnps)
            deallocate(dFx_dnngps)
            deallocate(dFx_dnggps)
            deallocate(dFx_dgggps)

            deallocate(dF_dgradtnae)
            deallocate(dFc_daa_ae)
            deallocate(dFc_dbb_ae)
            deallocate(dFc_dgg_ae)
            deallocate(dFc_dab_ae)
            deallocate(dFc_dag_ae)
            deallocate(dFc_dbg_ae)
            deallocate(dFc_daaa_ae)
            deallocate(dFc_dbbb_ae)
            deallocate(dFc_dggg_ae)
            deallocate(dFc_daab_ae)
            deallocate(dFc_daag_ae)
            deallocate(dFc_dabb_ae)
            deallocate(dFc_dbbg_ae)
            deallocate(dFc_dagg_ae)
            deallocate(dFc_dbgg_ae)
            deallocate(dFc_dabg_ae)

            deallocate(dF_dgradtnps)
            deallocate(dFc_daa_ps)
            deallocate(dFc_dbb_ps)
            deallocate(dFc_dgg_ps)
            deallocate(dFc_dab_ps)
            deallocate(dFc_dag_ps)
            deallocate(dFc_dbg_ps)
            deallocate(dFc_daaa_ps)
            deallocate(dFc_dbbb_ps)
            deallocate(dFc_dggg_ps)
            deallocate(dFc_daab_ps)
            deallocate(dFc_daag_ps)
            deallocate(dFc_dabb_ps)
            deallocate(dFc_dbbg_ps)
            deallocate(dFc_dagg_ps)
            deallocate(dFc_dbgg_ps)
            deallocate(dFc_dabg_ps)

            deallocate(gaga_ae_sph)
            deallocate(naga_ae_sph)
            deallocate(gg_ae_sph)
            deallocate(nag_ae_sph)

            deallocate(gaga_ps_sph)
            deallocate(naga_ps_sph)
            deallocate(gg_ps_sph)
            deallocate(nag_ps_sph)

            deallocate(dFxdgaovrga_ae_sph)
            deallocate(dFxdgaovrga_ps_sph)
            deallocate(dFcdgovrg_ae_sph)
            deallocate(dFcdgovrg_ps_sph)

            deallocate(dFadga_ae)
            deallocate(dFadga_ps)
            deallocate(dFadg_ae)
            deallocate(dFadg_ps)
            deallocate(dFadgaga_ae)
            deallocate(dFadgaga_ps)
            deallocate(dFadgg_ae)
            deallocate(dFadgg_ps)
            deallocate(dFadnaga_ae)
            deallocate(dFadnaga_ps)
            deallocate(dFadnag_ae)
            deallocate(dFadnag_ps)

            deallocate(dGadna_ae)
            deallocate(dGadna_ps)
            deallocate(dGadga_ae)
            deallocate(dGadga_ps)
            deallocate(dGadnana_ae)
            deallocate(dGadnana_ps)
            deallocate(dGadnaga_ae)
            deallocate(dGadnaga_ps)
            deallocate(dGadgaga_ae)
            deallocate(dGadgaga_ps)

            deallocate(dGdna_ae)
            deallocate(dGdna_ps)
            deallocate(dGdg_ae)
            deallocate(dGdg_ps)
            deallocate(dGdnana_ae)
            deallocate(dGdnana_ps)
            deallocate(dGdgg_ae)
            deallocate(dGdgg_ps)
            deallocate(dGdnag_ae)
            deallocate(dGdnag_ps)

            if(nspin.eq.2) then
                deallocate(gbgb_ae_sph)
                deallocate(nbgb_ae_sph)
                deallocate(nbg_ae_sph)

                deallocate(gbgb_ps_sph)
                deallocate(nbgb_ps_sph)
                deallocate(nbg_ps_sph)

                deallocate(dFxdgbovrgb_ae_sph)
                deallocate(dFxdgbovrgb_ps_sph)

                deallocate(dFbdgb_ae)
                deallocate(dFbdgb_ps)
                deallocate(dFbdg_ae)
                deallocate(dFbdg_ps)
                deallocate(dFbdgbgb_ae)
                deallocate(dFbdgbgb_ps)
                deallocate(dFbdgg_ae)
                deallocate(dFbdgg_ps)
                deallocate(dFbdnbgb_ae)
                deallocate(dFbdnbgb_ps)
                deallocate(dFbdnbg_ae)
                deallocate(dFbdnbg_ps)
                deallocate(dFbdag_ae)
                deallocate(dFbdag_ps)

                deallocate(dGbdnb_ae)
                deallocate(dGbdnb_ps)
                deallocate(dGbdgb_ae)
                deallocate(dGbdgb_ps)
                deallocate(dGbdnbnb_ae)
                deallocate(dGbdnbnb_ps)
                deallocate(dGbdnbgb_ae)
                deallocate(dGbdnbgb_ps)
                deallocate(dGbdgbgb_ae)
                deallocate(dGbdgbgb_ps)

                deallocate(dGdnb_ae)
                deallocate(dGdnb_ps)
                deallocate(dGdnbnb_ae)
                deallocate(dGdnbnb_ps)
                deallocate(dGdnanb_ae)
                deallocate(dGdnanb_ps)
                deallocate(dGdnbg_ae)
                deallocate(dGdnbg_ps)

            end if

!            deallocate(grad_nae)
!            deallocate(grad_tnae)
!            deallocate(dF_dnae)
!            deallocate(dF_dgradnae)
!            deallocate(grad_nps)
!            deallocate(grad_tnps)
!            deallocate(dF_dnps)
!            deallocate(dF_dgradnps)
!            deallocate(dnae_dr)
!            deallocate(dnae_dth)
!            deallocate(dnae_dph)
!            deallocate(ddnae_ddr)
!            deallocate(dnps_dr)
!            deallocate(dnps_dth)
!            deallocate(dnps_dph)
!            deallocate(ddnps_ddr)
!            deallocate(dgrad_tnae_dr)
!            deallocate(dgrad_tnps_dr)
!            deallocate(dF_dgradnae_dr)
!            deallocate(dF_dgradnps_dr)
        end if
        if(af /= 0) deallocate(flg_done)
        return
    end subroutine paw_xc_sphex_deallocate

! === For nrc decomposion. by takto 2012/12/05 =================================
    subroutine mult_sphex_element3_3D(nrc,dnr,mode  &
                                            ,n1,msphmx1,num_isph1,isph1 &
                                            ,n2,msphmx2,num_isph2,isph2 &
                                            ,n3,msphmx3,num_isph3,isph3 &
                                            ,ist, ien)
        integer,intent(in):: nrc,dnr,mode
        integer,intent(in):: msphmx1,num_isph1,isph1(25)
        integer,intent(in):: msphmx2,num_isph2,isph2(25)
        real(DP),intent(in):: n1(ista_nrc:iend_nrc,25),n2(ista_nrc:iend_nrc,25)
        real(DP),intent(out):: n3(ista_nrc:iend_nrc,25)
        integer,intent(out):: msphmx3,num_isph3,isph3(25)
        integer, intent(in) :: ist, ien

        integer:: ir,isp,isp1,isp2,n
        integer:: nsp,nsp1,nsp2
        real(DP):: fac,cijk
        logical:: flg_isp(25)

        n3=0.d0
        flg_isp=.false.
!        do isp=2,msphmx
        do nsp=2,num_isph1
            isp=isph1(nsp)
            do ir=ist,ien,dnr
                n3(ir,1)=n3(ir,1) + n1(ir,isp)*n2(ir,isp)
            end do
        end do

        do ir=ist,ien,dnr
            n3(ir,1)=n3(ir,1)/PAI4
        end do

        if(mode.eq.0) return

        flg_isp(1)=.true.

        msphmx3=0
!        do isp1=2,min(16,msphmx)                                    ! sphset2
!        do isp1=2,msphmx                                            ! sphset3
        do nsp1=2,num_isph1
            isp1=isph1(nsp1)
!            do isp2=isp1,min(16,msphmx)                             ! sphset2
!            do isp2=isp1,msphmx                                     ! sphset3
            do nsp2=nsp1,num_isph2
                isp2=isph2(nsp2)
                fac=1.d0;if(isp1.eq.isp2) fac=0.5d0
                do n=1,paw_mmt2(isp1,isp2)
                    isp=paw_isph2(isp1,isp2,n)
                    if(isp.eq.1) cycle
                    if(isp.gt.msphmx3) msphmx3=isp
                    flg_isp(isp)=.true.
                    cijk=paw_cr2(isp1,isp2,n)
! print *,isp,isp1,isp2
                    do ir=ist,ien,dnr
                        n3(ir,isp)=n3(ir,isp) + &
                            fac*(n1(ir,isp1)*n2(ir,isp2)+n1(ir,isp2)*n2(ir,isp1))*cijk
                    end do
                end do
            end do
        end do
        num_isph3=0
        isph3=0
        do isp=1,25
            if(flg_isp(isp)) then
                num_isph3=num_isph3+1
                isph3(num_isph3)=isp
            end if
        end do
        return

    end subroutine mult_sphex_element3_3D
! ==============================================================================

    ! mode = 0 : 0-th only 1 : all element
    subroutine mult_sphex_element_full(nrc,msphmx,mode,n1,n2,n3,msphmx2)
        integer,intent(in):: nrc,msphmx,mode
        real(DP),intent(in):: n1(nrc,25),n2(nrc,25)
        real(DP),intent(out):: n3(nrc,25)
        integer,intent(out):: msphmx2

        integer:: ir,isp,isp2,isp3,n
        real(DP):: fac,cijk

        n3=0.d0
        do isp=2,msphmx
            do ir=1,nrc
                n3(ir,1)=n3(ir,1) + n1(ir,isp)*n2(ir,isp)
            end do
        end do
        do ir=1,nrc
            n3(ir,1)=n3(ir,1)/PAI4+n1(ir,1)*n2(ir,1)
        end do

        if(mode.eq.0) return

        msphmx2=0
!        do isp2=2,min(16,msphmx)                                            ! sphset2
        do isp2=2,msphmx                                                    ! sphset3
!            do isp3=isp2,min(16,msphmx)                                     ! sphset2
            do isp3=isp2,msphmx                                             ! sphset3
                fac=1.d0;if(isp2.eq.isp3) fac=0.5d0
                do n=1,paw_mmt2(isp2,isp3)
                    isp=paw_isph2(isp2,isp3,n)
                    if(isp.eq.1) cycle
                    if(isp.gt.msphmx2) msphmx2=isp
                    cijk=paw_cr2(isp2,isp3,n)
! print *,isp,isp2,isp3
                    do ir=1,nrc
                        n3(ir,isp)=n3(ir,isp) + &
                            fac*(n1(ir,isp2)*n2(ir,isp3)+n1(ir,isp3)*n2(ir,isp2))*cijk
                    end do
                end do
            end do
        end do
        do isp=2,msphmx
            do ir=1,nrc
                n3(ir,isp)=n3(ir,isp) + n1(ir,1)*n2(ir,isp) + n1(ir,isp)*n2(ir,1)
            end do
        end do

        return

    end subroutine mult_sphex_element_full

   subroutine get_paw_sphex_integral(nrc,msphmx,rd,n1,n2,n3,msphmx2)
        integer,intent(in):: nrc,msphmx
        real(DP),intent(in):: n1(nrc,25),n2(nrc,25),rd(nrc)
        real(DP),intent(out):: n3(nrc,25)
        integer,intent(out):: msphmx2

        integer:: isp,isp2,isp3,n,ir
        real(DP):: dl,dl2,dl3,fac,cijk
        integer,allocatable,dimension(:):: il3

        allocate(il3(25));call substitute_il3(25,il3)

        n3=0.d0

!        do isp2=2,min(16,msphmx)                            ! sphset2
        do isp2=2,msphmx                                    ! sphset3
            dl2=dble(il3(isp2))
            dl2=dl2*(dl2+1.d0)
!            do isp3=isp2,min(16,msphmx)                     ! sphset2
            do isp3=isp2,msphmx                             ! sphset3
                dl3=dble(il3(isp3))
                dl3=dl3*(dl3+1.d0)
                fac=0.5d0;if(isp2.eq.isp3) fac=0.25d0
!                fac=0.5d0
                do n=1,paw_mmt2(isp2,isp3)
                    isp=paw_isph2(isp2,isp3,n)
                    if(isp.eq.1) cycle
                    if(isp.gt.msphmx2) msphmx2=isp
                    cijk=paw_cr2(isp2,isp3,n)
                    dl=dble(il3(isp))
                    dl=dl*(dl+1.d0)

                    do ir=1,nrc
                        n3(ir,isp)=n3(ir,isp) + &
                            fac*( &
                            n1(ir,isp2)*n2(ir,isp3)*(dl3+dl-dl2) + &
                            n1(ir,isp3)*n2(ir,isp2)*(dl2+dl-dl3))*cijk
                    end do

                end do
            end do
        end do

        do isp=2,msphmx
            dl=dble(il3(isp))
            dl=dl*(dl+1.d0)
            do ir=1,nrc
                n3(ir,isp)=(dl*n1(ir,1)*n2(ir,isp)+n3(ir,isp))/rd(ir)/rd(ir)
            end do
        end do
        return

    end subroutine get_paw_sphex_integral


    subroutine m_PAW_XC_cal_potential_sphex2(nfout,vflag)
        integer,intent(in):: nfout,vflag

        integer:: ia,it
        integer:: ir,is,nsph,ksph
!integer:: i
        integer:: ier,mode,itmp,ja
        integer:: nrc,msphmx_chg,msphmx_grd
        real(DP):: iga_ae,iga_ps,igb_ae,igb_ps,ig_ae,ig_ps
        real(DP):: dtnae_dr,dtnae_dth,dtnae_dph,ddtnae_ddr
        real(DP):: dtnps_dr,dtnps_dth,dtnps_dph,ddtnps_ddr
        real(DP):: sum1,sum2,sum3,sum4,sum5,sum6,sq4pi,zz
        integer:: nrc0,dnr
        integer:: num_isph_chg,num_isph_grd
        integer:: isph_chg(25),isph_grd(25)
        integer:: num_isph_n_n,isph_n_n(25)
        integer:: num_isph_n_g,isph_n_g(25)
        integer:: num_isph_g_g,isph_g_g(25)
        integer:: num_isph_nnn,isph_nnn(25)
        integer:: num_isph_ngg,isph_ngg(25)
        integer:: num_isph_ggg,isph_ggg(25)
        integer:: num_isph_all,isph_all(25)
        integer:: num_isph_2tm,isph_2tm(25)

        real(kind=DP) :: gnae,gnae2, gtnae, gtnae2, gnps, gnps2, gtnps,gtnps2

! =========================== added by K. Tagami ======================== 11.0
        real(kind=DP), allocatable :: magmom_tmp(:,:,:)
        real(kind=DP), allocatable :: rho_rad_tmp(:,:,:,:)
! ======================================================================= 11.0
        real(kind=DP), allocatable :: tmp_work(:)
        real(kind=DP), allocatable :: tmp_work2(:,:)

        integer, parameter :: PRINTLEVEL = 2
        integer, parameter :: DEBUGPRINTLEVEL = 3
        integer        :: id_sname = -1
        real(DP), allocatable, dimension(:,:,:,:) :: vxc_mpi
        real(DP) :: exc_mpi
        integer :: ista,iend
!$$#if defined(PARA3D) && defined(PAW3D)
        integer :: ia_add, ierr
! === For nrc decomposion. by takto 2012/12/07 =================================
        real(DP), allocatable, dimension(:,:) :: tmp0, tmp1, tmp2, tmp3, tmp4, tmp5
        integer :: ist, ien
! ==============================================================================
!        call tstatc0_begin('m_PAW_XC_cal_potential_sphex2',id_sname)
START_TIMER('m_PAW_XC_cal_potential_sphex2')
        call tstatc0_begin('m_PAW_XC_cal_potential_sphex2 ',id_sname,level=1)

! === For nrc decomposion. by takto 2012/12/05 =================================
        call paw_xc_sphex_allocate_3D
! ==============================================================================

! =========================== added by K. Tagami ======================== 11.0
! ======================================================================= 11.0

        exc_ae=0.d0;   exc_ps=0.d0
#ifdef __EDA__
        if ( sw_eda == ON ) then
           exc_ae_on_atom = 0.0d0;       exc_ps_on_atom = 0.0d0
        endif
#endif

!$$#if defined(PARA3D) && defined(PAW3D)
! === For nrc decomposion. by takto 2012/12/05 =================================
        vxc_ae_k_3D=0.d0
        vxc_ps_k_3D=0.d0

! === For nrc decomposion. by takto 2012/12/07 =================================
        if(paw_no_work_to_do) goto 999 ! This proc. has no nrc elements!
! ==============================================================================
        do ia = ista_natm, iend_natm
           ia_add = ia - ista_natm + 1
            if(af /= 0) then
                if(flg_done(ia)) cycle
            end if

            it=ityp(ia)
            if(ipaw(it)/=1 .or. &
                 & surface_integral_method(it).ne.SphericalHarmonicsExpansion) cycle
!!$                 & surface_integral_method(it).ne."SphericalHarmonicsExpansion") cycle
            dnr=paw_dnr(it)
            if(dnr.gt.1) then
                nrc0=wf_mnrc(it)
                nrc=1+int((nrc0-1)/dnr)*dnr
!                zz = (radr_paw(nrc0,it)-radr_paw(nrc,it))/ &
!                        (radr_paw(nrc+dnr,it)-radr_paw(nrc,it))
!                zz = log(radr_paw(nrc0,it)/radr_paw(nrc,it))/ &
!                        log(radr_paw(nrc+dnr,it)/radr_paw(nrc,it))
                zz = dble(nrc0-nrc)/dble(dnr)
                nrc=nrc+2*dnr                                                 ! 3rd
!               nrc=nrc+dnr                                                   !  1st
! === For nrc decomposion. by takto 2012/12/05 =================================
                ist = (ista_nrc+dnr-2)/dnr
                ist = dnr*ist+1
                ien = min(iend_nrc, nrc)
                call set_weight_exp3_3D(ier,1,nrc,dnr,radr_paw(1,it),zz,wos, &
                                        ista_nrc,iend_nrc,ist,ien)
! ==============================================================================
!               call set_weight_exp4(ier,1,nrc,dnr,radr_paw(:,it),zz,wos)      ! 1st
            else
               nrc=wf_mnrc(it)
! === For nrc decomposion. by takto 2012/12/05 =================================
               ist = (ista_nrc+dnr-2)/dnr
               ist = dnr*ist+1
               ien = min(iend_nrc, nrc)
               call set_weight_exp_3D(ier,1,nrc,radr_paw(1,it),wos, &
                                      ista_nrc,iend_nrc,ist,ien)
! ==============================================================================
            end if

!            nrc=1+int((nrc0-1)/dnr)*dnr
!            zz = (radr_paw(nrc0,it)-radr_paw(nrc,it))/(radr_paw(nrc+dnr,it)-radr_paw(nrc,it))
!            nrc=nrc+2*dnr
!!            nrc=nrc+dnr
!            call set_weight_exp3(ier,1,nrc,dnr,radr_paw(:,it),zz,wos)
!!            call set_weight_exp4(ier,1,nrc,dnr,radr_paw(:,it),zz,wos)

! === For nrc decomposion. by takto 2012/12/05 =================================
            do ir=ist,ien,dnr
! ==============================================================================
                wos(ir)=wos(ir)*radr_paw(ir,it)**2
            end do
            wos=wos*iwei(ia)
!!$            msphmx_chg=0

! ============================= modified by K. Tagami ================ 11.0
!       call m_PAWCD_set_ae_cd_sphex2 &
!            (ia,nspin,nrc,dnr,msph,nae_sph(1:nrc,1:nspin,1:25) &
!            ,msphmx_chg,num_isph_chg,isph_chg)

! === For nrc decomposion. by takto 2012/12/05 =================================
          call m_PAWCD_set_ae_cd_sphex2( ia, nspin, nrc, dnr, msph, nae_sph, msphmx_chg &
               &                       , num_isph_chg, isph_chg, ista_nrc,iend_nrc,ist )
! ==============================================================================
! ====================================================================== 11.0

            if(iprixc >= PRINTLEVEL) then
               write(nfout,'(" num_isph_chg = ",i8)') num_isph_chg
            end if
!do ir=1,nrc
!print '(25e19.6)',nae_sph(ir,1,1:25)
!end do
!stop

! ============================= modified by K. Tagami ================ 11.0
!       call m_PAWCD_set_ps_cd_sphex2 &
!            (ia,nspin,nrc,dnr,msph,nps_sph(1:nrc,1:nspin,1:25) &
!            ,msphmx_chg,num_isph_chg,isph_chg)

! === For nrc decomposion. by takto 2012/12/05 =================================
          call m_PAWCD_set_ps_cd_sphex2( ia, nspin, nrc, dnr, msph, nps_sph, msphmx_chg &
               &     , num_isph_chg, isph_chg, ista_nrc, iend_nrc, ist)
! ==============================================================================
! ====================================================================== 11.0

            if(iprixc >= PRINTLEVEL) then
               write(nfout,'(" num_isph_chg = ",i8)') num_isph_chg
            end if
!print *,num_isph_chg,isph_chg
!stop
!            call m_PAWCD_set_ae_cd_sphex &
!                            (ia,nspin,nrc,msph,nae_sph(1:nrc,1:nspin,1:25),msphmx_chg)
!            call m_PAWCD_set_ps_cd_sphex &
!                            (ia,nspin,nrc,msph,nps_sph(1:nrc,1:nspin,1:25),msphmx_chg)
!do i=1,nrc
!print '(7e19.7)' &
!,radr_paw(i,it) &
!,nae_sph(i,1,5)*radr_paw(i,it)**2 &
!,nps_sph(i,1,5)*radr_paw(i,it)**2 &
!,nae_sph(i,1,6)*radr_paw(i,it)**2 &
!,nps_sph(i,1,6)*radr_paw(i,it)**2 &
!,nae_sph(i,1,7)*radr_paw(i,it)**2 &
!,nps_sph(i,1,7)*radr_paw(i,it)**2
!end do
!stop

            if(check_of_xctype()==GGA) then
                call ggaxcp_paw_sphex2_3D()
! ==============================================================================
            else
!                call xcpotf_paw(nrc,nspin,ith,Valence_plus_PC_Charge)
            end if

            if(vflag == EXC_ONLY) then
                mode=0
            else if(vflag == VXC_AND_EXC) then
                mode=1
            end if

!            call set_sphex_elements2(nrc,dnr,msphmx_chg,msphmx_grd,mode)
! === For nrc decomposion. by takto 2012/12/05 =================================
            call set_sphex_elements2_3D(nrc,dnr,mode &
                                    ,msphmx_chg,num_isph_chg,isph_chg &
                                    ,msphmx_grd,num_isph_grd,isph_grd &
                                    ,num_isph_n_n,isph_n_n &
                                    ,num_isph_n_g,isph_n_g &
                                    ,num_isph_g_g,isph_g_g,ist,ien)
! ==============================================================================
!do ir=1,nrc
!print '(4e19.6)',radr_paw(ir,it) &
!                    ,exc_ae_field(ir) &
!                    ,gaga_ae_sph(ir,1) &
!                    ,gg_ae_sph(ir,1)
!end do
!do ir=1,nrc
!print '(9e19.6)',radr_paw(ir,it) &
!                    ,dFx_dnnae(ir,1) &
!                    ,dFc_daa_ae(ir) &
!                    ,dFx_dggae(ir,1) &
!                    ,dF_dgradnae(ir,1)/grad_nae(ir,1) &
!                    ,dFx_dngae(ir,1) &
!                    ,dFc_dgg_ae(ir) &
!                    ,dF_dgradtnae(ir)/grad_tnae(ir) &
!                    ,dFc_dag_ae(ir)
!end do
!stop
!            exc_ae=0.d0
!            exc_ps=0.d0

! === For nrc decomposion. by takto 2012/12/05 =================================
            do ir=ist,ien,dnr
! ==============================================================================
                sum1 = exc_ae_field(ir) + &
                        0.5d0*(dFx_dnnae(ir,1)+dFc_daa_ae(ir))*nana_ae_sph(ir,1)
                sum2 = exc_ps_field(ir) + &
                        0.5d0*(dFx_dnnps(ir,1)+dFc_daa_ps(ir))*nana_ps_sph(ir,1)
                if(nspin.eq.2) then
                    sum1 = sum1 + &
                        0.5d0*(dFx_dnnae(ir,2)+dFc_dbb_ae(ir))*nbnb_ae_sph(ir,1) + &
                        dFc_dab_ae(ir)*nanb_ae_sph(ir,1)
                    sum2 = sum2 + &
                        0.5d0*(dFx_dnnps(ir,2)+dFc_dbb_ps(ir))*nbnb_ps_sph(ir,1) + &
                        dFc_dab_ps(ir)*nanb_ps_sph(ir,1)
                end if

! ==== ASMS === 2017/11/08
!                if(check_of_xctype()==GGA .and. xctype /= 'ldapw91' &
!                     &                    .and. xctype /= 'ldapbe ' &
!                     &                    .and. xctype /= 'vdwdf') then
                if(check_of_xctype()==GGA .and. xctype /= 'ldapw91' &
                     &                    .and. xctype /= 'ldapbe ' ) then
! ==== ASMS === 2017/11/08

                    if(dabs(grad_nae(ir,1)) < DELTA10) cycle
                    if(dabs(grad_tnae(ir)) < DELTA10) cycle
                    if(dabs(grad_nps(ir,1)) < DELTA10) cycle
                    if(dabs(grad_tnps(ir)) < DELTA10) cycle
                    sum1 = sum1 + &
                        0.125d0*(dFx_dggae(ir,1)-dF_dgradnae(ir,1)/grad_nae(ir,1))/ &
                            grad_nae(ir,1)/grad_nae(ir,1)*gaga_ae_sph(ir,1) + &
                        0.5d0*dFx_dngae(ir,1)/grad_nae(ir,1)*naga_ae_sph(ir,1) + &
                        0.125d0*(dFc_dgg_ae(ir)-dF_dgradtnae(ir)/grad_tnae(ir))/ &
                            grad_tnae(ir)/grad_tnae(ir)*gg_ae_sph(ir,1) + &
                        0.5d0*dFc_dag_ae(ir)/grad_tnae(ir)*nag_ae_sph(ir,1)

                    sum2 = sum2 + &
                        0.125d0*(dFx_dggps(ir,1)-dF_dgradnps(ir,1)/grad_nps(ir,1))/ &
                            grad_nps(ir,1)/grad_nps(ir,1)*gaga_ps_sph(ir,1) + &
                        0.5d0*dFx_dngps(ir,1)/grad_nps(ir,1)*naga_ps_sph(ir,1) + &
                        0.125d0*(dFc_dgg_ps(ir)-dF_dgradtnps(ir)/grad_tnps(ir))/ &
                            grad_tnps(ir)/grad_tnps(ir)*gg_ps_sph(ir,1) + &
                        0.5d0*dFc_dag_ps(ir)/grad_tnps(ir)*nag_ps_sph(ir,1)
                    if(nspin.eq.2) then
                        if(dabs(grad_nae(ir,2)) < DELTA10) cycle
                        if(dabs(grad_nps(ir,2)) < DELTA10) cycle
                        sum1 = sum1 + &
                            0.125d0*(dFx_dggae(ir,2)-dF_dgradnae(ir,2)/grad_nae(ir,2))/ &
                            grad_nae(ir,2)/grad_nae(ir,2)*gbgb_ae_sph(ir,1) + &
                            0.5d0*dFx_dngae(ir,2)/grad_nae(ir,2)*nbgb_ae_sph(ir,1) + &
                            0.5d0*dFc_dbg_ae(ir)/grad_tnae(ir)*nbg_ae_sph(ir,1)
                        sum2 = sum2 + &
                            0.125d0*(dFx_dggps(ir,2)-dF_dgradnps(ir,2)/grad_nps(ir,2))/ &
                            grad_nps(ir,2)/grad_nps(ir,2)*gbgb_ps_sph(ir,1) + &
                            0.5d0*dFx_dngps(ir,2)/grad_nps(ir,2)*nbgb_ps_sph(ir,1) + &
                            0.5d0*dFc_dbg_ps(ir)/grad_tnps(ir)*nbg_ps_sph(ir,1)
                    end if
                end if
                exc_ae = exc_ae + sum1*wos(ir)
                exc_ps = exc_ps + sum2*wos(ir)
#ifdef __EDA__
                if ( sw_eda == ON ) then
                   exc_ae_on_atom(ia) = exc_ae_on_atom(ia) +sum1 *wos(ir)
                   exc_ps_on_atom(ia) = exc_ps_on_atom(ia) +sum2 *wos(ir)
                endif
#endif
             end do

!            exc_ae=exc_ae*PAI4*dble(af+1)
!            exc_ps=exc_ps*PAI4*dble(af+1)

!print *,'exc_ae = ',exc_ae
!print *,'exc_ps = ',exc_ps
!stop

            if(vflag == VXC_AND_EXC) then

! ==== ASMS === 2017/11/08
!                if(check_of_xctype()==GGA   .and. xctype /= 'ldapw91' &
!                                            .and. xctype /= 'ldapbe ' &
!                                            .and. xctype /= 'vdwdf' ) then
                if(check_of_xctype()==GGA   .and. xctype /= 'ldapw91' &
                                            .and. xctype /= 'ldapbe ' ) then
! ==== ASMS === 2017/11/08
! === For nrc decomposion. by takto 2012/12/05 =================================
                    do ir=ist,ien,dnr
! ==============================================================================
                        if(dabs(grad_nae(ir,1)) > 1.d-9) then
                            iga_ae = 1.d0/grad_nae(ir,1)
                        else
                            iga_ae = 0.d0
                        end if
                        if(dabs(grad_tnae(ir)) > 1.d-9) then
                            ig_ae = 1.d0/grad_tnae(ir)
                        else
                            ig_ae = 0.d0
                        end if
                        if(dabs(grad_nps(ir,1)) > 1.d-9) then
                            iga_ps = 1.d0/grad_nps(ir,1)
                        else
                            iga_ps = 0.d0
                        end if
                        if(dabs(grad_tnps(ir)) > 1.d-9) then
                            ig_ps = 1.d0/grad_tnps(ir)
                        else
                            ig_ps = 0.d0
                        end if
                        dFadga_ae(ir)   = 0.5d0*dFx_dngae(ir,1)*iga_ae
                        dFadga_ps(ir)   = 0.5d0*dFx_dngps(ir,1)*iga_ps
                        dFadg_ae(ir)    = 0.5d0*dFc_dag_ae(ir)*ig_ae
                        dFadg_ps(ir)    = 0.5d0*dFc_dag_ps(ir)*ig_ps
                        dFadgaga_ae(ir) = 0.25d0*(dFx_dnggae(ir,1) - dFx_dngae(ir,1)*iga_ae)*iga_ae**2
                        dFadgaga_ps(ir) = 0.25d0*(dFx_dnggps(ir,1) - dFx_dngps(ir,1)*iga_ps)*iga_ps**2
                        dFadgg_ae(ir)   = 0.25d0*(dFc_dagg_ae(ir) - dFc_dag_ae(ir)*ig_ae)*ig_ae**2
                        dFadgg_ps(ir)   = 0.25d0*(dFc_dagg_ps(ir) - dFc_dag_ps(ir)*ig_ps)*ig_ps**2
                        dFadnaga_ae(ir) = 0.5d0*dFx_dnngae(ir,1)*iga_ae
                        dFadnaga_ps(ir) = 0.5d0*dFx_dnngps(ir,1)*iga_ps
                        dFadnag_ae(ir)  = 0.5d0*dFc_daag_ae(ir)*ig_ae
                        dFadnag_ps(ir)  = 0.5d0*dFc_daag_ps(ir)*ig_ps

                        dGadna_ae(ir)   = dFx_dngae(ir,1)*iga_ae
                        dGadna_ps(ir)   = dFx_dngps(ir,1)*iga_ps
                        dGadga_ae(ir)   = 0.5d0*(dFx_dggae(ir,1) - dF_dgradnae(ir,1)*iga_ae)*iga_ae**2
                        dGadga_ps(ir)   = 0.5d0*(dFx_dggps(ir,1) - dF_dgradnps(ir,1)*iga_ps)*iga_ps**2
                        dGadnana_ae(ir) = dFx_dnngae(ir,1)*iga_ae
                        dGadnana_ps(ir) = dFx_dnngps(ir,1)*iga_ps
                        dGadnaga_ae(ir) = 0.5d0*(dFx_dnggae(ir,1) - dFx_dngae(ir,1)*iga_ae)*iga_ae**2
                        dGadnaga_ps(ir) = 0.5d0*(dFx_dnggps(ir,1) - dFx_dngps(ir,1)*iga_ps)*iga_ps**2
                        dGadgaga_ae(ir) = 0.25d0*(dFx_dgggae(ir,1) &
                                                    - 3.d0*dFx_dggae(ir,1)*iga_ae &
                                                    + 3.d0*dF_dgradnae(ir,1)*iga_ae**2)*iga_ae**3
                        dGadgaga_ps(ir) = 0.25d0*(dFx_dgggps(ir,1) &
                                                    - 3.d0*dFx_dggps(ir,1)*iga_ps &
                                                    + 3.d0*dF_dgradnps(ir,1)*iga_ps**2)*iga_ps**3

                        dGdna_ae(ir)    = dFc_dag_ae(ir)*ig_ae
                        dGdna_ps(ir)    = dFc_dag_ps(ir)*ig_ps
                        dGdg_ae(ir)     = 0.5d0*(dFc_dgg_ae(ir) - dF_dgradtnae(ir)*ig_ae)*ig_ae**2
                        dGdg_ps(ir)     = 0.5d0*(dFc_dgg_ps(ir) - dF_dgradtnps(ir)*ig_ps)*ig_ps**2
                        dGdnana_ae(ir)  = dFc_daag_ae(ir)*ig_ae
                        dGdnana_ps(ir)  = dFc_daag_ps(ir)*ig_ps
                        dGdgg_ae(ir)    = 0.25d0*(dFc_dggg_ae(ir) &
                                                    - 3.d0*dFc_dgg_ae(ir)*ig_ae &
                                                    + 3.d0*dF_dgradtnae(ir)*ig_ae**2)*ig_ae**3
                        dGdgg_ps(ir)    = 0.25d0*(dFc_dggg_ps(ir) &
                                                    - 3.d0*dFc_dgg_ps(ir)*ig_ps &
                                                    + 3.d0*dF_dgradtnps(ir)*ig_ps**2)*ig_ps**3
                        dGdnag_ae(ir)   = 0.5d0*(dFc_dagg_ae(ir) - dFc_dag_ae(ir)*ig_ae)*ig_ae**2
                        dGdnag_ps(ir)   = 0.5d0*(dFc_dagg_ps(ir) - dFc_dag_ps(ir)*ig_ps)*ig_ps**2

                        if(nspin.eq.2) then
                            if(dabs(grad_nae(ir,2)) > 1.d-9) then
                                igb_ae=1.d0/grad_nae(ir,2)
                            else
                                igb_ae=0.d0
                            end if
                            if(dabs(grad_nps(ir,2)) > 1.d-9) then
                                igb_ps=1.d0/grad_nps(ir,2)
                            else
                                igb_ps=0.d0
                            end if
                            dFbdgb_ae(ir)   = 0.5d0*dFx_dngae(ir,2)*igb_ae
                            dFbdgb_ps(ir)   = 0.5d0*dFx_dngps(ir,2)*igb_ps
                            dFbdg_ae(ir)    = 0.5d0*dFc_dbg_ae(ir)*ig_ae
                            dFbdg_ps(ir)    = 0.5d0*dFc_dbg_ps(ir)*ig_ps
                            dFbdgbgb_ae(ir) = 0.25d0*(dFx_dnggae(ir,2) - dFx_dngae(ir,2)*igb_ae)*igb_ae**2
                            dFbdgbgb_ps(ir) = 0.25d0*(dFx_dnggps(ir,2) - dFx_dngps(ir,2)*igb_ps)*igb_ps**2
                            dFbdgg_ae(ir)   = 0.25d0*(dFc_dbgg_ae(ir) - dFc_dbg_ae(ir)*ig_ae)*ig_ae**2
                            dFbdgg_ps(ir)   = 0.25d0*(dFc_dbgg_ps(ir) - dFc_dbg_ps(ir)*ig_ps)*ig_ps**2
                            dFbdnbgb_ae(ir) = 0.5d0*dFx_dnngae(ir,2)*igb_ae
                            dFbdnbgb_ps(ir) = 0.5d0*dFx_dnngps(ir,2)*igb_ps
                            dFbdnbg_ae(ir)  = 0.5d0*dFc_dbbg_ae(ir)*ig_ae
                            dFbdnbg_ps(ir)  = 0.5d0*dFc_dbbg_ps(ir)*ig_ps
                            dFbdag_ae(ir)   = 0.5d0*dFc_dabg_ae(ir)*ig_ae
                            dFbdag_ps(ir)   = 0.5d0*dFc_dabg_ps(ir)*ig_ps

                            dGbdnb_ae(ir)   = dFx_dngae(ir,2)*igb_ae
                            dGbdnb_ps(ir)   = dFx_dngps(ir,2)*igb_ps
                            dGbdgb_ae(ir)   = 0.5d0*(dFx_dggae(ir,2) - dF_dgradnae(ir,2)*igb_ae)*igb_ae**2
                            dGbdgb_ps(ir)   = 0.5d0*(dFx_dggps(ir,2) - dF_dgradnps(ir,2)*igb_ps)*igb_ps**2
                            dGbdnbnb_ae(ir) = dFx_dnngae(ir,2)*igb_ae
                            dGbdnbnb_ps(ir) = dFx_dnngps(ir,2)*igb_ps
                            dGbdnbgb_ae(ir) = 0.5d0*(dFx_dnggae(ir,2) - dFx_dngae(ir,2)*igb_ae)*igb_ae**2
                            dGbdnbgb_ps(ir) = 0.5d0*(dFx_dnggps(ir,2) - dFx_dngps(ir,2)*igb_ps)*igb_ps**2
                            dGbdgbgb_ae(ir) = 0.25d0*(dFx_dgggae(ir,2) &
                                                    - 3.d0*dFx_dggae(ir,2)*igb_ae &
                                                    + 3.d0*dF_dgradnae(ir,2)*igb_ae**2)*igb_ae**3
                            dGbdgbgb_ps(ir) = 0.25d0*(dFx_dgggps(ir,2) &
                                                    - 3.d0*dFx_dggps(ir,2)*igb_ps &
                                                    + 3.d0*dF_dgradnps(ir,2)*igb_ps**2)*igb_ps**3

                            dGdnb_ae(ir)    = dFc_dbg_ae(ir)*ig_ae
                            dGdnb_ps(ir)    = dFc_dbg_ps(ir)*ig_ps
                            dGdnbnb_ae(ir)  = dFc_dbbg_ae(ir)*ig_ae
                            dGdnbnb_ps(ir)  = dFc_dbbg_ps(ir)*ig_ps
                            dGdnanb_ae(ir)  = dFc_dabg_ae(ir)*ig_ae
                            dGdnanb_ps(ir)  = dFc_dabg_ps(ir)*ig_ps
                            dGdnbg_ae(ir)   = 0.5d0*(dFc_dbgg_ae(ir) - dFc_dbg_ae(ir)*ig_ae)*ig_ae**2
                            dGdnbg_ps(ir)   = 0.5d0*(dFc_dbgg_ps(ir) - dFc_dbg_ps(ir)*ig_ps)*ig_ps**2

                        end if

                    end do
                 end if

!                dFxcdna_ae_sph=0.d0
!                dFxcdna_ps_sph=0.d0
!                dFxcdnb_ae_sph=0.d0
!                dFxcdnb_ps_sph=0.d0
! ==== ASMS === 2017/11/08
!                if(check_of_xctype()==GGA   .and. xctype /= 'ldapw91' &
!                                            .and. xctype /= 'ldapbe ' &
!                                            .and. xctype /= 'vdwdf' ) then
                if(check_of_xctype()==GGA   .and. xctype /= 'ldapw91' &
                                            .and. xctype /= 'ldapbe ' ) then
! ==== ASMS === 2017/11/08
                    num_isph_2tm=num_isph_g_g
                    isph_2tm=isph_g_g
                else
                    num_isph_2tm=num_isph_n_n
                    isph_2tm=isph_n_n
                end if

! === For nrc decomposion. by takto 2012/12/05 =================================
                do ir=ist,ien,dnr
! ==============================================================================
!                    do nsph=1,msphmx_grd
!                    do ksph=1,num_isph_g_g
!                        nsph=isph_g_g(ksph)
                    do ksph=1,num_isph_2tm
                        nsph=isph_2tm(ksph)
                        sum1 = 0.5d0*(dFx_dnnnae(ir,1) + dFc_daaa_ae(ir))*nana_ae_sph(ir,nsph)
                        sum2 = 0.5d0*(dFx_dnnnps(ir,1) + dFc_daaa_ps(ir))*nana_ps_sph(ir,nsph)
                        if(nspin.eq.2) then
                            sum1 = sum1 + &
                                0.5d0*dFc_dabb_ae(ir)*nbnb_ae_sph(ir,nsph) + &
                                dFc_daab_ae(ir)*nanb_ae_sph(ir,nsph)
                            sum2 = sum2 + &
                                0.5d0*dFc_dabb_ps(ir)*nbnb_ps_sph(ir,nsph) + &
                                dFc_daab_ps(ir)*nanb_ps_sph(ir,nsph)

                            sum3 = 0.5d0*(dFx_dnnnae(ir,2) + dFc_dbbb_ae(ir))*nbnb_ae_sph(ir,nsph) + &
                                0.5d0*dFc_daab_ae(ir)*nana_ae_sph(ir,nsph) + &
                                dFc_dabb_ae(ir)*nanb_ae_sph(ir,nsph)
                            sum4 = 0.5d0*(dFx_dnnnps(ir,2) + dFc_dbbb_ps(ir))*nbnb_ps_sph(ir,nsph) + &
                                0.5d0*dFc_daab_ps(ir)*nana_ps_sph(ir,nsph) + &
                                dFc_dabb_ps(ir)*nanb_ps_sph(ir,nsph)
                        end if

! ==== ASMS === 2017/11/08
!                        if(check_of_xctype()==GGA   .and. xctype /= 'ldapw91' &
!                                                    .and. xctype /= 'ldapbe ' &
!                                                    .and. xctype /= 'vdwdf' ) then
                        if(check_of_xctype()==GGA   .and. xctype /= 'ldapw91' &
                                                    .and. xctype /= 'ldapbe ') then
! ==== ASMS === 2017/11/08
                            sum1 = sum1 + &
                                0.5d0*dFadgaga_ae(ir)*gaga_ae_sph(ir,nsph) + &
                                0.5d0*dFadgg_ae(ir)*gg_ae_sph(ir,nsph) + &
                                dFadnaga_ae(ir)*naga_ae_sph(ir,nsph) + &
                                dFadnag_ae(ir)*nag_ae_sph(ir,nsph)
                            sum2 = sum2 + &
                                0.5d0*dFadgaga_ps(ir)*gaga_ps_sph(ir,nsph) + &
                                0.5d0*dFadgg_ps(ir)*gg_ps_sph(ir,nsph) + &
                                dFadnaga_ps(ir)*naga_ps_sph(ir,nsph) + &
                                dFadnag_ps(ir)*nag_ps_sph(ir,nsph)

                            if(nspin.eq.2) then
                                sum1 = sum1 + dFbdag_ae(ir)*nbg_ae_sph(ir,nsph)
                                sum2 = sum2 + dFbdag_ps(ir)*nbg_ps_sph(ir,nsph)

                                sum3 = sum3 + &
                                    0.5d0*dFbdgbgb_ae(ir)*gbgb_ae_sph(ir,nsph) + &
                                    0.5d0*dFbdgg_ae(ir)*gg_ae_sph(ir,nsph) + &
                                    dFbdnbgb_ae(ir)*nbgb_ae_sph(ir,nsph) + &
                                    dFbdnbg_ae(ir)*nbg_ae_sph(ir,nsph) + &
                                    dFbdag_ae(ir)*nag_ae_sph(ir,nsph)
                                sum4 = sum4 + &
                                    0.5d0*dFbdgbgb_ps(ir)*gbgb_ps_sph(ir,nsph) + &
                                    0.5d0*dFbdgg_ps(ir)*gg_ps_sph(ir,nsph) + &
                                    dFbdnbgb_ps(ir)*nbgb_ps_sph(ir,nsph) + &
                                    dFbdnbg_ps(ir)*nbg_ps_sph(ir,nsph) + &
                                    dFbdag_ps(ir)*nag_ps_sph(ir,nsph)
                            end if
                        end if

                        if(nsph.eq.1) then
                            dFxcdna_ae_sph(ir,nsph) = dF_dnae(ir,1) + sum1
                            dFxcdna_ps_sph(ir,nsph) = dF_dnps(ir,1) + sum2
                        else
                            dFxcdna_ae_sph(ir,nsph) = (dFx_dnnae(ir,1)+dFc_daa_ae(ir))* &
                                                    nae_sph(ir,1,nsph) + &
                                                    dFadga_ae(ir)*grad_nae2_sph(ir,1,nsph) + &
                                                    dFadg_ae(ir)*grad_tnae2_sph(ir,nsph) + &
                                                    sum1
                            dFxcdna_ps_sph(ir,nsph) = (dFx_dnnps(ir,1)+dFc_daa_ps(ir))* &
                                                    nps_sph(ir,1,nsph) + &
                                                    dFadga_ps(ir)*grad_nps2_sph(ir,1,nsph) + &
                                                    dFadg_ps(ir)*grad_tnps2_sph(ir,nsph) + &
                                                    sum2
                        end if

                        if(nspin.eq.2) then
                            if(nsph.eq.1) then
                                dFxcdnb_ae_sph(ir,nsph) = dF_dnae(ir,2) + sum3
                                dFxcdnb_ps_sph(ir,nsph) = dF_dnps(ir,2) + sum4
                            else
                                dFxcdna_ae_sph(ir,nsph) = dFxcdna_ae_sph(ir,nsph) + &
                                                    dFc_dab_ae(ir)*nae_sph(ir,2,nsph)
                                dFxcdna_ps_sph(ir,nsph) = dFxcdna_ps_sph(ir,nsph) + &
                                                    dFc_dab_ps(ir)*nps_sph(ir,2,nsph)

                                dFxcdnb_ae_sph(ir,nsph) = (dFx_dnnae(ir,2)+dFc_dbb_ae(ir))* &
                                                    nae_sph(ir,2,nsph) + &
                                                    dFbdgb_ae(ir)*grad_nae2_sph(ir,2,nsph) + &
                                                    dFbdg_ae(ir)*grad_tnae2_sph(ir,nsph) + &
                                                    dFc_dab_ae(ir)*nae_sph(ir,1,nsph) + &
                                                    sum3
                                dFxcdnb_ps_sph(ir,nsph) = (dFx_dnnps(ir,2)+dFc_dbb_ps(ir))* &
                                                    nps_sph(ir,2,nsph) + &
                                                    dFbdgb_ps(ir)*grad_nps2_sph(ir,2,nsph) + &
                                                    dFbdg_ps(ir)*grad_tnps2_sph(ir,nsph) + &
                                                    dFc_dab_ps(ir)*nps_sph(ir,1,nsph) + &
                                                    sum4
                            end if
                        end if

                    end do
                end do

!                dFxdgaovrga_ae_sph=0.d0
!                dFxdgaovrga_ps_sph=0.d0
!                dFcdgovrg_ae_sph=0.d0
!                dFcdgovrg_ps_sph=0.d0

                do ir=ist,ien,dnr
! ==============================================================================
!                    do nsph=1,msphmx_grd
!                    do ksph=1,num_isph_g_g
!                        nsph=isph_g_g(ksph)
                    do ksph=1,num_isph_2tm
                        nsph=isph_2tm(ksph)
                        sum1 = 0.5d0*dGadnana_ae(ir)*nana_ae_sph(ir,nsph)
                        sum2 = 0.5d0*dGadnana_ps(ir)*nana_ps_sph(ir,nsph)
                        sum3 = 0.5d0*dGdnana_ae(ir)*nana_ae_sph(ir,nsph)
                        sum4 = 0.5d0*dGdnana_ps(ir)*nana_ps_sph(ir,nsph)

                        if(nspin.eq.2) then
                            sum3 = sum3 + &
                                0.5d0*dGdnbnb_ae(ir)*nbnb_ae_sph(ir,nsph) + &
                                dGdnanb_ae(ir)*nanb_ae_sph(ir,nsph)
                            sum4 = sum4 + &
                                0.5d0*dGdnbnb_ps(ir)*nbnb_ps_sph(ir,nsph) + &
                                dGdnanb_ps(ir)*nanb_ps_sph(ir,nsph)
                            sum5 = 0.5d0*dGbdnbnb_ae(ir)*nbnb_ae_sph(ir,nsph)
                            sum6 = 0.5d0*dGbdnbnb_ps(ir)*nbnb_ps_sph(ir,nsph)
                        end if

! ==== ASMS === 2017/11/08
!                        if(check_of_xctype()==GGA   .and. xctype /= 'ldapw91' &
!                                                    .and. xctype /= 'ldapbe ' &
!                                                    .and. xctype /= 'vdwdf' ) then
                        if(check_of_xctype()==GGA   .and. xctype /= 'ldapw91' &
                                                    .and. xctype /= 'ldapbe ') then
! ==== ASMS === 2017/11/08
                            sum1 = sum1 + &
                                0.5d0*dGadgaga_ae(ir)*gaga_ae_sph(ir,nsph) + &
                                dGadnaga_ae(ir)*naga_ae_sph(ir,nsph)
                            sum2 = sum2 + &
                                0.5d0*dGadgaga_ps(ir)*gaga_ps_sph(ir,nsph) + &
                                dGadnaga_ps(ir)*naga_ps_sph(ir,nsph)
                            sum3 = sum3 + &
                                0.5d0*dGdgg_ae(ir)*gg_ae_sph(ir,nsph) + &
                                dGdnag_ae(ir)*nag_ae_sph(ir,nsph)
                            sum4 = sum4 + &
                                0.5d0*dGdgg_ps(ir)*gg_ps_sph(ir,nsph) + &
                                dGdnag_ps(ir)*nag_ps_sph(ir,nsph)
                            if(nspin.eq.2) then
                                sum3 = sum3 + dGdnbg_ae(ir)*nbg_ae_sph(ir,nsph)
                                sum4 = sum4 + dGdnbg_ps(ir)*nbg_ps_sph(ir,nsph)
                                sum5 = sum5 + &
                                    0.5d0*dGbdgbgb_ae(ir)*gbgb_ae_sph(ir,nsph) + &
                                    dGbdnbgb_ae(ir)*nbgb_ae_sph(ir,nsph)
                                sum6 = sum6 + &
                                    0.5d0*dGbdgbgb_ps(ir)*gbgb_ps_sph(ir,nsph) + &
                                    dGbdnbgb_ps(ir)*nbgb_ps_sph(ir,nsph)
                            end if

                        end if

                        if(nsph.eq.1) then
                            if(dabs(grad_nae(ir,1)) > 1.d-9) then
                                dFxdgaovrga_ae_sph(ir,nsph) = dF_dgradnae(ir,1)/grad_nae(ir,1) + sum1
                            else
                                dFxdgaovrga_ae_sph(ir,nsph) = sum1
                            end if
                            if(dabs(grad_nps(ir,1)) > 1.d-9) then
                                dFxdgaovrga_ps_sph(ir,nsph) = dF_dgradnps(ir,1)/grad_nps(ir,1) + sum2
                            else
                                dFxdgaovrga_ps_sph(ir,nsph) = sum2
                            end if
                            if(dabs(grad_tnae(ir)) > 1.d-9) then
                                dFcdgovrg_ae_sph(ir,nsph) = dF_dgradtnae(ir)/grad_tnae(ir) + sum3
                            else
                                dFcdgovrg_ae_sph(ir,nsph) = sum3
                            end if
                            if(dabs(grad_tnps(ir)) > 1.d-9) then
                                dFcdgovrg_ps_sph(ir,nsph) = dF_dgradtnps(ir)/grad_tnps(ir) + sum4
                            else
                                dFcdgovrg_ps_sph(ir,nsph) = sum4
                            end if
                        else
                            dFxdgaovrga_ae_sph(ir,nsph) = dGadna_ae(ir)*nae_sph(ir,1,nsph) + &
                                                            dGadga_ae(ir)*grad_nae2_sph(ir,1,nsph) + &
                                                            sum1
                            dFxdgaovrga_ps_sph(ir,nsph) = dGadna_ps(ir)*nps_sph(ir,1,nsph) + &
                                                            dGadga_ps(ir)*grad_nps2_sph(ir,1,nsph) + &
                                                            sum2
                            dFcdgovrg_ae_sph(ir,nsph) = dGdna_ae(ir)*nae_sph(ir,1,nsph) + &
                                                        dGdg_ae(ir)*grad_tnae2_sph(ir,nsph) + &
                                                        sum3
                            dFcdgovrg_ps_sph(ir,nsph) = dGdna_ps(ir)*nps_sph(ir,1,nsph) + &
                                                        dGdg_ps(ir)*grad_tnps2_sph(ir,nsph) + &
                                                        sum4
                        end if

                        if(nspin.eq.2) then
                            if(nsph.eq.1) then
                                if(dabs(grad_nae(ir,2)) > 1.d-9) then
                                    dFxdgbovrgb_ae_sph(ir,nsph) = dF_dgradnae(ir,2)/grad_nae(ir,2) + sum5
                                else
                                    dFxdgbovrgb_ae_sph(ir,nsph) = sum5
                                end if
                                if(dabs(grad_nps(ir,2)) > 1.d-9) then
                                    dFxdgbovrgb_ps_sph(ir,nsph) = dF_dgradnps(ir,2)/grad_nps(ir,2) + sum6
                                else
                                    dFxdgbovrgb_ps_sph(ir,nsph) = sum6
                                end if
                            else
                                dFxdgbovrgb_ae_sph(ir,nsph) = dGbdnb_ae(ir)*nae_sph(ir,2,nsph) + &
                                                            dGbdgb_ae(ir)*grad_nae2_sph(ir,2,nsph) + &
                                                            sum5
                                dFxdgbovrgb_ps_sph(ir,nsph) = dGbdnb_ps(ir)*nps_sph(ir,2,nsph) + &
                                                            dGbdgb_ps(ir)*grad_nps2_sph(ir,2,nsph) + &
                                                            sum6
                                dFcdgovrg_ae_sph(ir,nsph) = dFcdgovrg_ae_sph(ir,nsph) + &
                                                            dGdnb_ae(ir)*nae_sph(ir,2,nsph)
                                dFcdgovrg_ps_sph(ir,nsph) = dFcdgovrg_ps_sph(ir,nsph) + &
                                                            dGdnb_ps(ir)*nps_sph(ir,2,nsph)
                            end if
                        end if

                    end do
                end do

                nana_ae_sph=0.d0
                nana_ps_sph=0.d0
                nag_ae_sph=0.d0
                nag_ps_sph=0.d0
                naga_ae_sph=0.d0
                naga_ps_sph=0.d0
                gg_ae_sph=0.d0
                gg_ps_sph=0.d0
                nanb_ae_sph=0.d0
                nanb_ps_sph=0.d0
                gaga_ae_sph=0.d0
                gaga_ps_sph=0.d0
                if(nspin.eq.2) then
                    nbnb_ae_sph=0.d0
                    nbnb_ps_sph=0.d0
                    nbg_ae_sph=0.d0
                    nbg_ps_sph=0.d0
                    nbgb_ae_sph=0.d0
                    nbgb_ps_sph=0.d0
                end if

!                do nsph=1,msphmx_grd
!                do ksph=1,num_isph_g_g
!                    nsph=isph_g_g(ksph)
! === For nrc decomposion. by takto 2012/12/07 =================================
                allocate(tmp0(ista_nrc-6*dnr:iend_nrc+6*dnr,25))
                allocate(tmp1(ista_nrc-6*dnr:iend_nrc+6*dnr,25))
                allocate(tmp2(ista_nrc-6*dnr:iend_nrc+6*dnr,25))
                allocate(tmp3(ista_nrc-6*dnr:iend_nrc+6*dnr,25))
                call boundary_exchange_dim2(tmp0, dFxdgaovrga_ae_sph, dnr, 25)
                call boundary_exchange_dim2(tmp1, dFxdgaovrga_ps_sph, dnr, 25)
                call boundary_exchange_dim2(tmp2, dFcdgovrg_ae_sph,   dnr, 25)
                call boundary_exchange_dim2(tmp3, dFcdgovrg_ps_sph,   dnr, 25)
                if(nspin.eq.2) then
                   allocate(tmp4(ista_nrc-6*dnr:iend_nrc+6*dnr,25))
                   allocate(tmp5(ista_nrc-6*dnr:iend_nrc+6*dnr,25))
                   call boundary_exchange_dim2(tmp4, dFxdgbovrgb_ae_sph, dnr, 25)
                   call boundary_exchange_dim2(tmp5, dFxdgbovrgb_ps_sph, dnr, 25)
                end if
! ==============================================================================
                do ksph=1,num_isph_2tm
                   nsph=isph_2tm(ksph)
                   call calc_diff_exp3_3D(ier,3,nrc,dnr,radr_paw(1,it), &
                                          tmp0(ista_nrc-6*dnr,nsph),    &
                                          nana_ae_sph(ista_nrc,nsph),   &
                                          ista_nrc, iend_nrc, ist, ien)
                   call calc_diff_exp3_3D(ier,3,nrc,dnr,radr_paw(1,it), &
                                          tmp1(ista_nrc-6*dnr,nsph),    &
                                          nana_ps_sph(ista_nrc,nsph),   &
                                          ista_nrc, iend_nrc, ist, ien)
                   call calc_diff_exp3_3D(ier,3,nrc,dnr,radr_paw(1,it), &
                                          tmp2(ista_nrc-6*dnr,nsph),    &
                                          gg_ae_sph(ista_nrc,nsph),     &
                                          ista_nrc, iend_nrc, ist, ien)
                   call calc_diff_exp3_3D(ier,3,nrc,dnr,radr_paw(1,it), &
                                          tmp3(ista_nrc-6*dnr,nsph),    &
                                          gg_ps_sph(ista_nrc,nsph),     &
                                          ista_nrc, iend_nrc, ist, ien)
                end do

                if(nspin.eq.2) then
                    do ksph=1,num_isph_2tm
                        nsph=isph_2tm(ksph)
                        call calc_diff_exp3_3D(ier,3,nrc,dnr,radr_paw(1,it), &
                                               tmp4(ista_nrc-6*dnr,nsph),    &
                                               nbnb_ae_sph(ista_nrc,nsph),   &
                                               ista_nrc, iend_nrc, ist, ien)
                        call calc_diff_exp3_3D(ier,3,nrc,dnr,radr_paw(1,it), &
                                               tmp5(ista_nrc-6*dnr,nsph),    &
                                               nbnb_ps_sph(ista_nrc,nsph),   &
                                               ista_nrc, iend_nrc, ist, ien)
                    end do
                end if
! === For nrc decomposion. by takto 2012/12/07 =================================
                deallocate(tmp0)
                deallocate(tmp1)
                deallocate(tmp2)
                deallocate(tmp3)
                if(nspin.eq.2) then
                   deallocate(tmp4)
                   deallocate(tmp5)
                end if
! ==============================================================================

! === For nrc decomposion. by takto 2012/12/05 =================================
                call mult_sphex_element_full3_3D(nrc,dnr,1                        &
                                                ,dnae_dr_sph(:,1,1:25)            &
                                                ,itmp,num_isph_chg,isph_chg(1:25) &
                                                ,nana_ae_sph                      &
                                                ,itmp,num_isph_n_n,isph_n_n(1:25) &
                                                ,nag_ae_sph                       &
                                                ,itmp,num_isph_nnn,isph_nnn(1:25),ist,ien)
                call mult_sphex_element_full3_3D(nrc,dnr,1                        &
                                                ,dnps_dr_sph(:,1,1:25)            &
                                                ,itmp,num_isph_chg,isph_chg(1:25) &
                                                ,nana_ps_sph                      &
                                                ,itmp,num_isph_n_n,isph_n_n(1:25) &
                                                ,nag_ps_sph                       &
                                                ,itmp,num_isph_nnn,isph_nnn(1:25),ist,ien)
                if(nspin.eq.1) then
                    call mult_sphex_element_full3_3D(nrc,dnr,1                        &
                                                    ,dnae_dr_sph(:,1,1:25)            &
                                                    ,itmp,num_isph_chg,isph_chg(1:25) &
                                                    ,gg_ae_sph                        &
                                                    ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                                    ,nanb_ae_sph                      &
                                                    ,itmp,num_isph_ngg,isph_ngg(1:25),ist,ien)
                    call mult_sphex_element_full3_3D(nrc,dnr,1                        &
                                                    ,dnps_dr_sph(:,1,1:25)            &
                                                    ,itmp,num_isph_chg,isph_chg(1:25) &
                                                    ,gg_ps_sph                        &
                                                    ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                                    ,nanb_ps_sph                      &
                                                    ,itmp,num_isph_ngg,isph_ngg(1:25),ist,ien)
                else
                    call mult_sphex_element_full3_3D(nrc,dnr,1                        &
                                                    ,dnae_dr_sph(:,2,1:25)            &
                                                    ,itmp,num_isph_chg,isph_chg(1:25) &
                                                    ,nbnb_ae_sph                      &
                                                    ,itmp,num_isph_n_n,isph_n_n(1:25) &
                                                    ,nbg_ae_sph                       &
                                                    ,itmp,num_isph_nnn,isph_nnn(1:25),ist,ien)
                    call mult_sphex_element_full3_3D(nrc,dnr,1                        &
                                                    ,dnps_dr_sph(:,2,1:25)            &
                                                    ,itmp,num_isph_chg,isph_chg(1:25) &
                                                    ,nbnb_ps_sph                      &
                                                    ,itmp,num_isph_n_n,isph_n_n(1:25) &
                                                    ,nbg_ps_sph                       &
                                                    ,itmp,num_isph_nnn,isph_nnn(1:25),ist,ien)
                    call mult_sphex_element_full3_3D(nrc,dnr,1                        &
                                                    ,dnae_dr_sph(:,1,1:25) +          &
                                                     dnae_dr_sph(:,2,1:25)            &
                                                    ,itmp,num_isph_chg,isph_chg(1:25) &
                                                    ,gg_ae_sph                        &
                                                    ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                                    ,nanb_ae_sph                      &
                                                    ,itmp,num_isph_ngg,isph_ngg(1:25),ist,ien)
                    call mult_sphex_element_full3_3D(nrc,dnr,1                        &
                                                    ,dnps_dr_sph(:,1,1:25) +          &
                                                     dnps_dr_sph(:,2,1:25)            &
                                                    ,itmp,num_isph_chg,isph_chg(1:25) &
                                                    ,gg_ps_sph                        &
                                                    ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                                    ,nanb_ps_sph                      &
                                                    ,itmp,num_isph_ngg,isph_ngg(1:25),ist,ien)
                end if
! ==============================================================================

                nana_ae_sph=0.d0
                nana_ps_sph=0.d0
!                do nsph=1,msphmx_grd
                do ksph=1,num_isph_chg
                    nsph=isph_chg(ksph)
! === For nrc decomposion. by takto 2012/12/05 =================================
                do ir=ist,ien,dnr
! ==============================================================================
                    nana_ae_sph(ir,nsph) = ddnae_ddr_sph(ir,1,nsph) + &
                                            2.d0*dnae_dr_sph(ir,1,nsph)/radr_paw(ir,it)
                    nana_ps_sph(ir,nsph) = ddnps_ddr_sph(ir,1,nsph) + &
                                            2.d0*dnps_dr_sph(ir,1,nsph)/radr_paw(ir,it)
                end do
                end do

                gg_ae_sph=0.d0
                gg_ps_sph=0.d0
! === Debug by Intel "-check all" option! by T.Kato 2011/03/01 =================
                if(nspin.eq.2) then
! ==============================================================================
                nbnb_ae_sph=0.d0
                nbnb_ps_sph=0.d0
! === Debug by Intel "-check all" option! by T.Kato 2011/03/01 =================
                endif
! ==============================================================================
                if(nspin.eq.1) then
!                    do nsph=1,msphmx_grd
                    do ksph=1,num_isph_chg
                        nsph=isph_chg(ksph)
! === For nrc decomposion. by takto 2012/12/05 =================================
                    do ir=ist,ien,dnr
! ==============================================================================
                        gg_ae_sph(ir,nsph) = ddnae_ddr_sph(ir,1,nsph) + &
                                                2.d0*dnae_dr_sph(ir,1,nsph)/radr_paw(ir,it)
                        gg_ps_sph(ir,nsph) = ddnps_ddr_sph(ir,1,nsph) + &
                                                2.d0*dnps_dr_sph(ir,1,nsph)/radr_paw(ir,it)
                    end do
                    end do
                else
!                    do nsph=1,msphmx_grd
                    do ksph=1,num_isph_chg
                        nsph=isph_chg(ksph)
! === For nrc decomposion. by takto 2012/12/05 =================================
                    do ir=ist,ien,dnr
! ==============================================================================
                        nbnb_ae_sph(ir,nsph) = ddnae_ddr_sph(ir,2,nsph) + &
                                                2.d0*dnae_dr_sph(ir,2,nsph)/radr_paw(ir,it)
                        nbnb_ps_sph(ir,nsph) = ddnps_ddr_sph(ir,2,nsph) + &
                                                2.d0*dnps_dr_sph(ir,2,nsph)/radr_paw(ir,it)
                    end do
                    end do

!                    do nsph=1,msphmx_grd
                    do ksph=1,num_isph_chg
                        nsph=isph_chg(ksph)
! === For nrc decomposion. by takto 2012/12/05 =================================
                    do ir=ist,ien,dnr
! ==============================================================================
                        gg_ae_sph(ir,nsph) = ddnae_ddr_sph(ir,1,nsph) + ddnae_ddr_sph(ir,2,nsph) + &
                                                2.d0*(dnae_dr_sph(ir,1,nsph)+dnae_dr_sph(ir,2,nsph)) &
                                                /radr_paw(ir,it)
                        gg_ps_sph(ir,nsph) = ddnps_ddr_sph(ir,1,nsph) + ddnps_ddr_sph(ir,2,nsph) + &
                                                2.d0*(dnps_dr_sph(ir,1,nsph)+dnps_dr_sph(ir,2,nsph)) &
                                                /radr_paw(ir,it)
                    end do
                    end do
                end if

! === For nrc decomposion. by takto 2012/12/05 =================================
                call mult_sphex_element_full3_3D(nrc,dnr,1                        &
                                                ,nana_ae_sph                      &
                                                ,itmp,num_isph_chg,isph_chg(1:25) &
                                                ,dFxdgaovrga_ae_sph               &
                                                ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                                ,naga_ae_sph                      &
                                                ,itmp,num_isph_ngg,isph_ngg(1:25),ist,ien)
                call mult_sphex_element_full3_3D(nrc,dnr,1                        &
                                                ,nana_ps_sph                      &
                                                ,itmp,num_isph_chg,isph_chg(1:25) &
                                                ,dFxdgaovrga_ps_sph               &
                                                ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                                ,naga_ps_sph                      &
                                                ,itmp,num_isph_ngg,isph_ngg(1:25),ist,ien)
                call mult_sphex_element_full3_3D(nrc,dnr,1                        &
                                                ,gg_ae_sph                        &
                                                ,itmp,num_isph_chg,isph_chg(1:25) &
                                                ,dFcdgovrg_ae_sph                 &
                                                ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                                ,gaga_ae_sph                      &
                                                ,itmp,num_isph_ngg,isph_ngg(1:25),ist,ien)
                call mult_sphex_element_full3_3D(nrc,dnr,1                        &
                                                ,gg_ps_sph                        &
                                                ,itmp,num_isph_chg,isph_chg(1:25) &
                                                ,dFcdgovrg_ps_sph                 &
                                                ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                                ,gaga_ps_sph                      &
                                                ,itmp,num_isph_ngg,isph_ngg(1:25),ist,ien)
                if(nspin.eq.2) then
                    call mult_sphex_element_full3_3D(nrc,dnr,1                        &
                                                    ,nbnb_ae_sph                      &
                                                    ,itmp,num_isph_chg,isph_chg(1:25) &
                                                    ,dFxdgbovrgb_ae_sph               &
                                                    ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                                    ,nbgb_ae_sph                      &
                                                    ,itmp,num_isph_ngg,isph_ngg(1:25),ist,ien)
                    call mult_sphex_element_full3_3D(nrc,dnr,1                        &
                                                    ,nbnb_ps_sph                      &
                                                    ,itmp,num_isph_chg,isph_chg(1:25) &
                                                    ,dFxdgbovrgb_ps_sph               &
                                                    ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                                    ,nbgb_ps_sph                      &
                                                    ,itmp,num_isph_ngg,isph_ngg(1:25),ist,ien)
                end if
! ==============================================================================

! === For nrc decomposion. by takto 2012/12/05 =================================
                call get_paw_sphex_integral3_3D(nrc,dnr,radr_paw(1,it)           &
                                               ,dFxdgaovrga_ae_sph               &
                                               ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                               ,nae_sph(:,1,1:25)                &
                                               ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                               ,nana_ae_sph                      &
                                               ,itmp,num_isph_ggg,isph_ggg(1:25),ist,ien)
                call get_paw_sphex_integral3_3D(nrc,dnr,radr_paw(1,it)           &
                                               ,dFxdgaovrga_ps_sph               &
                                               ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                               ,nps_sph(:,1,1:25)                &
                                               ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                               ,nana_ps_sph                      &
                                               ,itmp,num_isph_ggg,isph_ggg(1:25),ist,ien)
                if(nspin.eq.1) then
                    call get_paw_sphex_integral3_3D(nrc,dnr,radr_paw(1,it)           &
                                                   ,dFcdgovrg_ae_sph                 &
                                                   ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                                   ,nae_sph(:,1,1:25)                &
                                                   ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                                   ,gg_ae_sph                        &
                                                   ,itmp,num_isph_ggg,isph_ggg(1:25),ist,ien)
                    call get_paw_sphex_integral3_3D(nrc,dnr,radr_paw(1,it)           &
                                                   ,dFcdgovrg_ps_sph                 &
                                                   ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                                   ,nps_sph(:,1,1:25)                &
                                                   ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                                   ,gg_ps_sph                        &
                                                   ,itmp,num_isph_ggg,isph_ggg(1:25),ist,ien)
                else
                    call get_paw_sphex_integral3_3D(nrc,dnr,radr_paw(1,it)           &
                                                   ,dFxdgbovrgb_ae_sph               &
                                                   ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                                   ,nae_sph(:,2,1:25)                &
                                                   ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                                   ,nbnb_ae_sph                      &
                                                   ,itmp,num_isph_ggg,isph_ggg(1:25),ist,ien)
                    call get_paw_sphex_integral3_3D(nrc,dnr,radr_paw(1,it)           &
                                                   ,dFxdgbovrgb_ps_sph               &
                                                   ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                                   ,nps_sph(:,2,1:25)                &
                                                   ,itmp,num_isph_g_g,isph_g_g(1:25) &
                                                   ,nbnb_ps_sph                      &
                                                   ,itmp,num_isph_ggg,isph_ggg(1:25),ist,ien)
                    call get_paw_sphex_integral3_3D(nrc,dnr,radr_paw(1,it)                  &
                                                   ,dFcdgovrg_ae_sph                        &
                                                   ,itmp,num_isph_g_g,isph_g_g(1:25)        &
                                                   ,(nae_sph(:,1,1:25)+nae_sph(:,2,1:25))   &
                                                   ,itmp,num_isph_g_g,isph_g_g(1:25)        &
                                                   ,gg_ae_sph                               &
                                                   ,itmp,num_isph_ggg,isph_ggg(1:25),ist,ien)
                    call get_paw_sphex_integral3_3D(nrc,dnr,radr_paw(1,it)                  &
                                                   ,dFcdgovrg_ps_sph                        &
                                                   ,itmp,num_isph_g_g,isph_g_g(1:25)        &
                                                   ,(nps_sph(:,1,1:25)+nps_sph(:,2,1:25))   &
                                                   ,itmp,num_isph_g_g,isph_g_g(1:25)        &
                                                   ,gg_ps_sph                               &
                                                   ,itmp,num_isph_ggg,isph_ggg(1:25),ist,ien)
                end if
! ==============================================================================

                call merge_isph_flgs(num_isph_n_n,isph_n_n &
                                    ,num_isph_n_g,isph_n_g &
                                    ,num_isph_g_g,isph_g_g &
                                    ,num_isph_nnn,isph_nnn &
                                    ,num_isph_ngg,isph_ngg &
                                    ,num_isph_ggg,isph_ggg &
                                    ,num_isph_all,isph_all)

                sq4pi=sqrt(PAI4)

! === For nrc decomposion. by takto 2012/12/05 =================================
                do ksph=1,num_isph_all
                    nsph=isph_all(ksph)
                    if(nsph.gt.msph) cycle
                    do ir=ist,ien,dnr
                    vxc_ae_k_3D(ir,1,nsph,ia_add)  = dFxcdna_ae_sph(ir,nsph) - &
                                              nag_ae_sph(ir,nsph) - &
                                              naga_ae_sph(ir,nsph) + &
                                              nana_ae_sph(ir,nsph) - &
                                              nanb_ae_sph(ir,nsph) - &
                                              gaga_ae_sph(ir,nsph) + &
                                              gg_ae_sph(ir,nsph)
                    vxc_ps_k_3D(ir,1,nsph,ia_add)  = dFxcdna_ps_sph(ir,nsph)  - &
                                              nag_ps_sph(ir,nsph) - &
                                              naga_ps_sph(ir,nsph) + &
                                              nana_ps_sph(ir,nsph) - &
                                              nanb_ps_sph(ir,nsph) - &
                                              gaga_ps_sph(ir,nsph) + &
                                              gg_ps_sph(ir,nsph)
                    end do
                end do
                if(nspin.eq.2) then
                    do ksph=1,num_isph_all
                        nsph=isph_all(ksph)
                        if(nsph.gt.msph) cycle
                        do ir=ist,ien,dnr
                        vxc_ae_k_3D(ir,2,nsph,ia_add)  = dFxcdnb_ae_sph(ir,nsph) - &
                                                  nbg_ae_sph(ir,nsph) - &
                                                  nbgb_ae_sph(ir,nsph) + &
                                                  nbnb_ae_sph(ir,nsph) - &
                                                  nanb_ae_sph(ir,nsph) - &
                                                  gaga_ae_sph(ir,nsph) + &
                                                  gg_ae_sph(ir,nsph)
                        vxc_ps_k_3D(ir,2,nsph,ia_add)  = dFxcdnb_ps_sph(ir,nsph) - &
                                                  nbg_ps_sph(ir,nsph) - &
                                                  nbgb_ps_sph(ir,nsph) + &
                                                  nbnb_ps_sph(ir,nsph) - &
                                                  nanb_ps_sph(ir,nsph) - &
                                                  gaga_ps_sph(ir,nsph) + &
                                                  gg_ps_sph(ir,nsph)
                        end do
                    end do
                end if

                do ir=ist,ien,dnr
                    vxc_ae_k_3D(ir,1,1,ia_add)= vxc_ae_k_3D(ir,1,1,ia_add)*sq4pi
                    vxc_ps_k_3D(ir,1,1,ia_add)= vxc_ps_k_3D(ir,1,1,ia_add)*sq4pi
                end do

                if(nspin.eq.2) then
                    do ir=ist,ien,dnr
                        vxc_ae_k_3D(ir,2,1,ia_add)= vxc_ae_k_3D(ir,2,1,ia_add)*sq4pi
                        vxc_ps_k_3D(ir,2,1,ia_add)= vxc_ps_k_3D(ir,2,1,ia_add)*sq4pi
                    end do
                end if
! ==============================================================================

             end if

        end do
! === For nrc decomposion. by takto 2012/12/07 =================================
999     continue ! This proc. has no nrc elements!
! ==============================================================================
!$$#if defined(PARA3D) && defined(PAW3D)
! === For nrc decomposion. by takto 2012/12/05 =================================
!       call mpi_allreduce(MPI_IN_PLACE,exc_ae,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_nrc_world,ierr)
!       call mpi_allreduce(MPI_IN_PLACE,exc_ps,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_nrc_world,ierr)
        call mpi_allreduce(MPI_IN_PLACE,exc_ae,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_CommGroup,ierr)
        call mpi_allreduce(MPI_IN_PLACE,exc_ps,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_CommGroup,ierr)
! ==============================================================================
!!$        call decomp_vxc_ae_k_r_3D(vxc_ae_k,vxc_ae_k_3D,mmesh,msph,natm)
!!$        call decomp_vxc_ae_k_r_3D(vxc_ps_k,vxc_ps_k_3D,mmesh,msph,natm)

!ASMS        exc_ae=exc_ae*PAI4*dble(af+1)
!ASMS        exc_ps=exc_ps*PAI4*dble(af+1)
        exc_ae=exc_ae*PAI4 !ASMS
        exc_ps=exc_ps*PAI4 !ASMS

#ifdef __EDA__
        if ( sw_eda == ON ) then
           exc_ae_on_atom = exc_ae_on_atom *PAI4
           exc_ps_on_atom = exc_ps_on_atom *PAI4

           if ( npes>1 ) then !ASMS
              call mpi_allreduce( MPI_IN_PLACE, exc_ae_on_atom, natm, &
                   &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
              call mpi_allreduce( MPI_IN_PLACE, exc_ps_on_atom, natm, &
                   &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
           endif
        endif
#endif

        if(iprixc >= PRINTLEVEL) then
           write(nfout,'(" exc_ae, exc_ps = ",2d20.8)') exc_ae, exc_ps
        end if
        call paw_xc_sphex_deallocate

! ============================== added by K. Tagami ======================== 11.0
! === DEBUG by tkato 2012/11/30 ================================================
! ==============================================================================
! ========================================================================= 11.0

STOP_TIMER('m_PAW_XC_cal_potential_sphex2')
        call tstatc0_end(id_sname)

contains

! === For nrc decomposion. by takto 2012/12/05 =================================
        subroutine ggaxcp_paw_sphex2_3D
            integer :: pot_type
    !integer:: i
            grad_tnae2_sph = 0.d0
            grad_nae2_sph = 0.d0
            grad_tnps2_sph = 0.d0
            grad_nps2_sph = 0.d0
            dF_dnae = 0.d0
            dF_dgradnae = 0.d0
            dF_dnps = 0.d0
            dF_dgradnps = 0.d0

!            if(xctype == 'ldapw91' .or. xctype == 'ldapbe '.or. xctype == 'vdwdf' ) then
            if(xctype == 'ldapw91' .or. xctype == 'ldapbe ' ) then
                grad_nae2_sph=0.d0;grad_tnae2_sph=0.d0
                grad_nps2_sph=0.d0;grad_tnps2_sph=0.d0
                grad_nae=0.d0;grad_tnae=0.d0
                grad_nps=0.d0;grad_tnps=0.d0
                call m_PAWCD_set_cr2_isph2_mmt2()
                msphmx_grd=0
                num_isph_grd=0
            else
                call abs_grad_rho_ud_paw_sphex2_3D
            end if

          if(xctype == 'ggapw91' .or. xctype == 'ldapw91') then
              call ex_ggapw91_paw_drv2_3D(nrc,dnr,nspin    &
                                    ,nae_sph(ista_nrc,1,1) &
                                    ,grad_nae              &
                                    ,exc_ae_field          &
                                    ,dF_dnae               &
                                    ,dF_dgradnae           &
                                    ,dFx_dnnae             &
                                    ,dFx_dngae             &
                                    ,dFx_dggae             &
                                    ,dFx_dnnnae            &
                                    ,dFx_dnngae            &
                                    ,dFx_dnggae            &
                                    ,dFx_dgggae            &
                                    ,ista_nrc, iend_nrc, ist, ien)
              call cr_ggapw91_paw_drv2_3D(nrc,dnr,nspin &
                                    ,nae_sph(ista_nrc,1,1) &
                                    ,grad_tnae             &
                                    ,exc_ae_field          &
                                    ,dF_dnae               &
                                    ,dF_dgradtnae          &
                                    ,dFc_daa_ae            &
                                    ,dFc_dbb_ae            &
                                    ,dFc_dgg_ae            &
                                    ,dFc_dab_ae            &
                                    ,dFc_dag_ae            &
                                    ,dFc_dbg_ae            &
                                    ,dFc_daaa_ae           &
                                    ,dFc_dbbb_ae           &
                                    ,dFc_dggg_ae           &
                                    ,dFc_daab_ae           &
                                    ,dFc_daag_ae           &
                                    ,dFc_dabb_ae           &
                                    ,dFc_dbbg_ae           &
                                    ,dFc_dagg_ae           &
                                    ,dFc_dbgg_ae           &
                                    ,dFc_dabg_ae           &
                                    ,ista_nrc, iend_nrc, ist, ien)

              call ex_ggapw91_paw_drv2_3D(nrc,dnr,nspin    &
                                    ,nps_sph(ista_nrc,1,1) &
                                    ,grad_nps              &
                                    ,exc_ps_field          &
                                    ,dF_dnps               &
                                    ,dF_dgradnps           &
                                    ,dFx_dnnps             &
                                    ,dFx_dngps             &
                                    ,dFx_dggps             &
                                    ,dFx_dnnnps            &
                                    ,dFx_dnngps            &
                                    ,dFx_dnggps            &
                                    ,dFx_dgggps            &
                                    ,ista_nrc, iend_nrc, ist, ien)
              call cr_ggapw91_paw_drv2_3D(nrc,dnr,nspin    &
                                    ,nps_sph(ista_nrc,1,1) &
                                    ,grad_tnps             &
                                    ,exc_ps_field          &
                                    ,dF_dnps               &
                                    ,dF_dgradtnps          &
                                    ,dFc_daa_ps            &
                                    ,dFc_dbb_ps            &
                                    ,dFc_dgg_ps            &
                                    ,dFc_dab_ps            &
                                    ,dFc_dag_ps            &
                                    ,dFc_dbg_ps            &
                                    ,dFc_daaa_ps           &
                                    ,dFc_dbbb_ps           &
                                    ,dFc_dggg_ps           &
                                    ,dFc_daab_ps           &
                                    ,dFc_daag_ps           &
                                    ,dFc_dabb_ps           &
                                    ,dFc_dbbg_ps           &
                                    ,dFc_dagg_ps           &
                                    ,dFc_dbgg_ps           &
                                    ,dFc_dabg_ps           &
                                    ,ista_nrc, iend_nrc, ist, ien)
!          else if(xctype == 'ggapbe ' .or. xctype == 'ldapbe ' .or. xctype == 'vdwdf') then
          else if(xctype == 'ggapbe ' .or. xctype == 'ldapbe ' ) then
              call ex_ggapbe_paw_drv2_3D(nrc,dnr,nspin     &
                                    ,nae_sph(ista_nrc,1,1) &
                                    ,grad_nae              &
                                    ,exc_ae_field          &
                                    ,dF_dnae               &
                                    ,dF_dgradnae           &
                                    ,dFx_dnnae             &
                                    ,dFx_dngae             &
                                    ,dFx_dggae             &
                                    ,dFx_dnnnae            &
                                    ,dFx_dnngae            &
                                    ,dFx_dnggae            &
                                    ,dFx_dgggae            &
                                    ,ista_nrc, iend_nrc, ist, ien)

                 if(iprixc >= DEBUGPRINTLEVEL) then
                    write(nfout,'(" ia = ", i8)') ia
                    write(nfout,'(" -- exc_ae_field(1:20) -- ")')
                    write(nfout,'(8d16.4)') exc_ae_field(1:20)
                    write(nfout,'(" -- exc_ae_field(nrc-19:nrc) --")')
                    write(nfout,'(8d16.4)') exc_ae_field(nrc-19:nrc)
                    write(nfout,'(" -- nae_sph(1:20,1,1) -- ")')
                    write(nfout,'(8d16.4)') nae_sph(1:20,1,1)
                    write(nfout,'(" -- nae_sph(nrc-19:nrc,1,1) --")')
                    write(nfout,'(8d16.4)') nae_sph(nrc-19:nrc,1,1)
                    if(nspin == 2) then
                       write(nfout,'(" -- nae_sph(1:20,2,1) --")')
                       write(nfout,'(8d16.4)') nae_sph(1:20,2,1)
                       write(nfout,'(" -- nae_sph(nrc-19:nrc,2,1) --")')
                       write(nfout,'(8d16.4)') nae_sph(nrc-19:nrc,2,1)
                    end if
                    write(nfout,'(" -- dFx_dggae(1:20,1) --")')
                    write(nfout,'(8d16.4)') dFx_dggae(1:20,1)
                    write(nfout,'(" -- dFx_dggae(nrc-19:nrc,1) --")')
                    write(nfout,'(8d16.4)') dFx_dggae(nrc-19:nrc,1)
                    if(nspin == 2) then
                       write(nfout,'(" -- dFx_dggae(1:20,2) --")')
                       write(nfout,'(8d16.4)') dFx_dggae(1:20,2)
                       write(nfout,'(" -- dFx_dggae(nrc-19:nrc,2) --")')
                       write(nfout,'(8d16.4)') dFx_dggae(nrc-19:nrc,2)
                    end if
                 end if

              call cr_ggapbe_paw_drv2_3D(nrc,dnr,nspin     &
                                    ,nae_sph(ista_nrc,1,1) &
                                    ,grad_tnae             &
                                    ,exc_ae_field          &
                                    ,dF_dnae               &
                                    ,dF_dgradtnae          &
                                    ,dFc_daa_ae            &
                                    ,dFc_dbb_ae            &
                                    ,dFc_dgg_ae            &
                                    ,dFc_dab_ae            &
                                    ,dFc_dag_ae            &
                                    ,dFc_dbg_ae            &
                                    ,dFc_daaa_ae           &
                                    ,dFc_dbbb_ae           &
                                    ,dFc_dggg_ae           &
                                    ,dFc_daab_ae           &
                                    ,dFc_daag_ae           &
                                    ,dFc_dabb_ae           &
                                    ,dFc_dbbg_ae           &
                                    ,dFc_dagg_ae           &
                                    ,dFc_dbgg_ae           &
                                    ,dFc_dabg_ae           &
                                    ,ista_nrc, iend_nrc, ist, ien)
                 if(iprixc >= DEBUGPRINTLEVEL) then
                    write(nfout,'(" -- grad_tnae(1:20) -- ")')
                    write(nfout,'(8d16.4)') grad_tnae(1:20)
                    write(nfout,'(" -- grad_tnae(nrc-19:nrc) --")')
                    write(nfout,'(8d16.4)') grad_tnae(nrc-19:nrc)
                    write(nfout,'(" -- dFc_dbg_ae(1:20) --")')
                    write(nfout,'(8d16.4)') dFc_dbg_ae(1:20)
                    write(nfout,'(" -- dFc_dbg_ase(nrc-19:nrc) --")')
                    write(nfout,'(8d16.4)') dFc_dbg_ae(nrc-19:nrc)
                 end if

              call ex_ggapbe_paw_drv2_3D(nrc,dnr,nspin     &
                                    ,nps_sph(ista_nrc,1,1) &
                                    ,grad_nps              &
                                    ,exc_ps_field          &
                                    ,dF_dnps               &
                                    ,dF_dgradnps           &
                                    ,dFx_dnnps             &
                                    ,dFx_dngps             &
                                    ,dFx_dggps             &
                                    ,dFx_dnnnps            &
                                    ,dFx_dnngps            &
                                    ,dFx_dnggps            &
                                    ,dFx_dgggps            &
                                    ,ista_nrc, iend_nrc, ist, ien)

                 if(iprixc >= DEBUGPRINTLEVEL) then
                    write(nfout,'(" -- exc_ps_field(1:20) -- ")')
                    write(nfout,'(8d16.4)') exc_ps_field(1:20)
                    write(nfout,'(" -- exc_ps_field(nrc-19:nrc) --")')
                    write(nfout,'(8d16.4)') exc_ps_field(nrc-19:nrc)
                    write(nfout,'(" -- dFx_dggps(1:20,1) --")')
                    write(nfout,'(8d16.4)') dFx_dggps(1:20,1)
                    write(nfout,'(" -- dFx_dggps(nrc-19:nrc,1) --")')
                    write(nfout,'(8d16.4)') dFx_dggps(nrc-19:nrc,1)
                    if(nspin == 2) then
                       write(nfout,'(" -- dFx_dggps(1:20,2) --")')
                       write(nfout,'(8d16.4)') dFx_dggps(1:20,2)
                       write(nfout,'(" -- dFx_dggps(nrc-19:nrc,2) --")')
                       write(nfout,'(8d16.4)') dFx_dggps(nrc-19:nrc,2)
                    end if
                 end if

              call cr_ggapbe_paw_drv2_3D(nrc,dnr,nspin     &
                                    ,nps_sph(ista_nrc,1,1) &
                                    ,grad_tnps             &
                                    ,exc_ps_field          &
                                    ,dF_dnps               &
                                    ,dF_dgradtnps          &
                                    ,dFc_daa_ps            &
                                    ,dFc_dbb_ps            &
                                    ,dFc_dgg_ps            &
                                    ,dFc_dab_ps            &
                                    ,dFc_dag_ps            &
                                    ,dFc_dbg_ps            &
                                    ,dFc_daaa_ps           &
                                    ,dFc_dbbb_ps           &
                                    ,dFc_dggg_ps           &
                                    ,dFc_daab_ps           &
                                    ,dFc_daag_ps           &
                                    ,dFc_dabb_ps           &
                                    ,dFc_dbbg_ps           &
                                    ,dFc_dagg_ps           &
                                    ,dFc_dbgg_ps           &
                                    ,dFc_dabg_ps           &
                                    ,ista_nrc, iend_nrc, ist, ien)

                 if(iprixc >= DEBUGPRINTLEVEL) then
                    write(nfout,'(" -- exc_ps_field(1:20) -- ")')
                    write(nfout,'(8d16.4)') exc_ps_field(1:20)
                    write(nfout,'(" -- exc_ps_field(nrc-19:nrc) --")')
                    write(nfout,'(8d16.4)') exc_ps_field(nrc-19:nrc)
                    write(nfout,'(" -- grad_tnps(1:20) -- ")')
                    write(nfout,'(8d16.4)') grad_tnps(1:20)
                    write(nfout,'(" -- grad_tnps(nrc-19:nrc) --")')
                    write(nfout,'(8d16.4)') grad_tnps(nrc-19:nrc)
                    write(nfout,'(" -- dFc_dbg_ps(1:20) --")')
                    write(nfout,'(8d16.4)') dFc_dbg_ps(1:20)
                    write(nfout,'(" -- dFc_dbg_ps(nrc-19:nrc) --")')
                    write(nfout,'(8d16.4)') dFc_dbg_ps(nrc-19:nrc)
                 end if


!             call ex_ggapbe_paw (nrc,nspin,nae(1:nrc,1:nspin) &
!                                        ,grad_nae(1:nrc,1:nspin) &
!                                        ,wos(1:nrc),texc &
!                                        ,dF_dnae(1:nrc,1:nspin) &
!                                        ,dF_dgradnae(1:nrc,1:nspin))
!             call cr_ggapbe_paw (nrc,nspin,nae(1:nrc,1:nspin) &
!                                        ,grad_tnae(1:nrc) &
!                                        ,wos(1:nrc),texc &
!                                        ,dF_dnae(1:nrc,1:nspin))
!             exc_ae=exc_ae+texc*omg_wght(ith)
!             call ex_ggapbe_paw (nrc,nspin,nps(1:nrc,1:nspin) &
!                                        ,grad_nps(1:nrc,1:nspin) &
!                                        ,wos(1:nrc),texc &
!                                        ,dF_dnps(1:nrc,1:nspin) &
!                                        ,dF_dgradnps(1:nrc,1:nspin))
!             call cr_ggapbe_paw (nrc,nspin,nps(1:nrc,1:nspin) &
!                                        ,grad_tnps(1:nrc) &
!                                        ,wos(1:nrc),texc &
!                                        ,dF_dnps(1:nrc,1:nspin))
!             exc_ps=exc_ps+texc*omg_wght(ith)
          else if(xctype == 'katopbe' .or. xctype == 'ggapbek') then
!             call ex_ggapbe_paw (nrc,nspin,nae(1:nrc,1:nspin) &
!                                        ,grad_nae(1:nrc,1:nspin) &
!                                        ,wos(1:nrc),texc &
!                                        ,dF_dnae(1:nrc,1:nspin) &
!                                        ,dF_dgradnae(1:nrc,1:nspin))
!             call cr_ggapbe_paw (nrc,nspin,nae(1:nrc,1:nspin) &
!                                        ,grad_tnae(1:nrc) &
!                                        ,wos(1:nrc),texc &
!                                        ,dF_dnae(1:nrc,1:nspin))
!             exc_ae=exc_ae+texc*omg_wght(ith)
!             call ex_ggapbe_paw (nrc,nspin,nps(1:nrc,1:nspin) &
!                                        ,grad_nps(1:nrc,1:nspin) &
!                                        ,wos(1:nrc),texc &
!                                        ,dF_dnps(1:nrc,1:nspin) &
!                                        ,dF_dgradnps(1:nrc,1:nspin))
!             call cr_ggapbe_paw (nrc,nspin,nps(1:nrc,1:nspin) &
!                                        ,grad_tnps(1:nrc) &
!                                        ,wos(1:nrc),texc &
!                                        ,dF_dnps(1:nrc,1:nspin))
!             exc_ps=exc_ps+texc*omg_wght(ith)
          else if(xctype == 'ggabp  ') then
!             call xclda_paw(nrc,nspin,nae(1:nrc,1:nspin) &
!                                        ,wos(1:nrc),texc &
!                                        ,dF_dnae(1:nrc,1:nspin))
!             call ggabek_paw(nrc,nspin,nae(1:nrc,1:nspin) &
!                                        ,grad_nae(1:nrc,1:nspin) &
!                                        ,wos(1:nrc),texc &
!                                        ,dF_dnae(1:nrc,1:nspin) &
!                                        ,dF_dgradnae(1:nrc,1:nspin))
!             call ggaprd_paw(nrc,nspin,nae(1:nrc,1:nspin) &
!                                        ,grad_nae(1:nrc,1:nspin) &
!                                        ,wos(1:nrc),texc &
!                                        ,dF_dnae(1:nrc,1:nspin) &
!                                        ,dF_dgradnae(1:nrc,1:nspin))
!             exc_ae=exc_ae+texc*omg_wght(ith)
!             call xclda_paw(nrc,nspin,nps(1:nrc,1:nspin) &
!                                        ,wos(1:nrc),texc &
!                                        ,dF_dnps(1:nrc,1:nspin))
!             call ggabek_paw(nrc,nspin,nps(1:nrc,1:nspin) &
!                                        ,grad_nps(1:nrc,1:nspin) &
!                                        ,wos(1:nrc),texc &
!                                        ,dF_dnps(1:nrc,1:nspin) &
!                                        ,dF_dgradnps(1:nrc,1:nspin))
!             call ggaprd_paw(nrc,nspin,nps(1:nrc,1:nspin) &
!                                        ,grad_nps(1:nrc,1:nspin) &
!                                        ,wos(1:nrc),texc &
!                                        ,dF_dnps(1:nrc,1:nspin) &
!                                        ,dF_dgradnps(1:nrc,1:nspin))
!             exc_ps=exc_ps+texc*omg_wght(ith)
          else if ( xctype == 'revpbe' .or. xctype == 'rpbe' .or. &
               &    xctype == 'wc06'   .or. xctype == 'htbs' .or. &
               &    xctype == 'pbesol' .or. xctype == 'pbeint' .or. &
               &    xctype == 'ev93' .or. xctype == 'evpw91' .or. &
               &    xctype == 'lb94' ) then

             if ( xctype == 'ggapbe' )  pot_type = 1
             if ( xctype == 'revpbe' )  pot_type = 2
             if ( xctype == 'rpbe' )    pot_type = 3
             if ( xctype == 'wc06' )    pot_type = 4
             if ( xctype == 'htbs' )    pot_type = 5
             if ( xctype == 'pbesol' )    pot_type = 6
             if ( xctype == 'pbeint' )    pot_type = 7
             if ( xctype == 'ev93' )    pot_type = 20
             if ( xctype == 'evpw91' )    pot_type = 20
             if ( xctype == 'lb94' )    pot_type = 25

             call ex_gga_paw_library_3D( nrc, dnr, nspin &
                  ,nae_sph(ista_nrc,1,1), grad_nae &
                  ,exc_ae_field &
                  ,dF_dnae,    dF_dgradnae &
                  ,dFx_dnnae,  dFx_dngae   &
                  ,dFx_dggae,  dFx_dnnnae  &
                  ,dFx_dnngae, dFx_dnggae  &
                  ,dFx_dgggae, ista_nrc, iend_nrc, ist, ien, pot_type )
             call ex_gga_paw_library_3D(nrc,dnr,nspin &
                  ,nps_sph(ista_nrc,1,1), grad_nps &
                  ,exc_ps_field &
                  ,dF_dnps,    dF_dgradnps &
                  ,dFx_dnnps,  dFx_dngps   &
                  ,dFx_dggps,  dFx_dnnnps  &
                  ,dFx_dnngps, dFx_dnggps  &
                  ,dFx_dgggps, ista_nrc, iend_nrc, ist, ien, pot_type )

             if ( xctype == 'ev93' .or. xctype == 'lb94' ) then
                allocate( tmp_work(ista_nrc:iend_nrc) ); tmp_work = 0.0d0
                call cr_gga_paw_library_3D(nrc,dnr,nspin &
                     ,nae_sph(ista_nrc,1,1), tmp_work, exc_ae_field &
                     ,dF_dnae, dF_dgradtnae &
                     ,dFc_daa_ae, dFc_dbb_ae, dFc_dgg_ae &
                     ,dFc_dab_ae, dFc_dag_ae, dFc_dbg_ae &
                     ,dFc_daaa_ae, dFc_dbbb_ae &
                     ,dFc_dggg_ae, dFc_daab_ae &
                     ,dFc_daag_ae, dFc_dabb_ae &
                     ,dFc_dbbg_ae, dFc_dagg_ae &
                     ,dFc_dbgg_ae, dFc_dabg_ae, ista_nrc, iend_nrc, ist, ien, pot_type)
                call cr_gga_paw_library_3D(nrc,dnr,nspin &
                     ,nps_sph(ista_nrc,1,1), tmp_work, exc_ps_field &
                     ,dF_dnps, dF_dgradtnps &
                     ,dFc_daa_ps, dFc_dbb_ps, dFc_dgg_ps &
                     ,dFc_dab_ps, dFc_dag_ps, dFc_dbg_ps &
                     ,dFc_daaa_ps, dFc_dbbb_ps &
                     ,dFc_dggg_ps, dFc_daab_ps &
                     ,dFc_daag_ps, dFc_dabb_ps &
                     ,dFc_dbbg_ps, dFc_dagg_ps &
                     ,dFc_dbgg_ps, dFc_dabg_ps, ista_nrc, iend_nrc, ist, ien, pot_type )
             else if ( xctype == 'evpw91' ) then
                call cr_ggapw91_paw_drv2_3D(nrc,dnr,nspin &
                     ,nae_sph(ista_nrc,1,1), grad_tnae, exc_ae_field  &
                     ,dF_dnae, dF_dgradtnae &
                     ,dFc_daa_ae, dFc_dbb_ae, dFc_dgg_ae &
                     ,dFc_dab_ae, dFc_dag_ae, dFc_dbg_ae &
                     ,dFc_daaa_ae, dFc_dbbb_ae &
                     ,dFc_dggg_ae, dFc_daab_ae &
                     ,dFc_daag_ae, dFc_dabb_ae &
                     ,dFc_dbbg_ae, dFc_dagg_ae &
                     ,dFc_dbgg_ae, dFc_dabg_ae, ista_nrc, iend_nrc, ist, ien )
                call cr_ggapw91_paw_drv2_3D(nrc,dnr,nspin &
                     ,nps_sph(ista_nrc,1,1), grad_tnps, exc_ps_field  &
                     ,dF_dnps, dF_dgradtnps &
                     ,dFc_daa_ps, dFc_dbb_ps, dFc_dgg_ps &
                     ,dFc_dab_ps, dFc_dag_ps, dFc_dbg_ps &
                     ,dFc_daaa_ps, dFc_dbbb_ps &
                     ,dFc_dggg_ps, dFc_daab_ps &
                     ,dFc_daag_ps, dFc_dabb_ps &
                     ,dFc_dbbg_ps, dFc_dagg_ps &
                     ,dFc_dbgg_ps, dFc_dabg_ps, ista_nrc, iend_nrc, ist, ien )
             else
                call cr_gga_paw_library_3D(nrc,dnr,nspin &
                     ,nae_sph(ista_nrc,1,1), grad_tnae, exc_ae_field &
                     ,dF_dnae, dF_dgradtnae &
                     ,dFc_daa_ae, dFc_dbb_ae, dFc_dgg_ae &
                     ,dFc_dab_ae, dFc_dag_ae, dFc_dbg_ae &
                     ,dFc_daaa_ae, dFc_dbbb_ae &
                     ,dFc_dggg_ae, dFc_daab_ae &
                     ,dFc_daag_ae, dFc_dabb_ae &
                     ,dFc_dbbg_ae, dFc_dagg_ae &
                     ,dFc_dbgg_ae, dFc_dabg_ae, ista_nrc, iend_nrc, ist, ien, pot_type)
                call cr_gga_paw_library_3D(nrc,dnr,nspin &
                     ,nps_sph(ista_nrc,1,1), grad_tnps, exc_ps_field &
                     ,dF_dnps, dF_dgradtnps &
                     ,dFc_daa_ps, dFc_dbb_ps, dFc_dgg_ps &
                     ,dFc_dab_ps, dFc_dag_ps, dFc_dbg_ps &
                     ,dFc_daaa_ps, dFc_dbbb_ps &
                     ,dFc_dggg_ps, dFc_daab_ps &
                     ,dFc_daag_ps, dFc_dabb_ps &
                     ,dFc_dbbg_ps, dFc_dagg_ps &
                     ,dFc_dbgg_ps, dFc_dabg_ps, ista_nrc, iend_nrc, ist, ien, pot_type )
             endif

          else if ( xctype == 'vdwdf' ) then
             if ( exchange_pot_type == 'pbe' )  pot_type = 1
             if ( exchange_pot_type == 'revpbe' )  pot_type = 2
             if ( exchange_pot_type == 'b86r' )  pot_type = 11
             if ( exchange_pot_type == 'optpbe' )  pot_type = 12
             if ( exchange_pot_type == 'optb86b' )  pot_type = 13
             if ( exchange_pot_type == 'pw86r' )  pot_type = 14
             if ( exchange_pot_type == 'c09x' )  pot_type = 15
             if ( exchange_pot_type == 'lvpw86r' )  pot_type = 16

             allocate( tmp_work(nrc) ); tmp_work = 0.0d0

             call ex_gga_paw_library_3D( nrc, dnr, nspin &
                  ,nae_sph(ista_nrc,1,1), grad_nae &
                  ,exc_ae_field &
                  ,dF_dnae,    dF_dgradnae &
                  ,dFx_dnnae,  dFx_dngae   &
                  ,dFx_dggae,  dFx_dnnnae  &
                  ,dFx_dnngae, dFx_dnggae  &
                  ,dFx_dgggae, ista_nrc, iend_nrc, ist, ien, pot_type )
             call cr_gga_paw_library_3D(nrc,dnr,nspin &
                  ,nae_sph(ista_nrc,1,1), tmp_work, exc_ae_field &
                  ,dF_dnae, dF_dgradtnae &
                  ,dFc_daa_ae, dFc_dbb_ae, dFc_dgg_ae &
                  ,dFc_dab_ae, dFc_dag_ae, dFc_dbg_ae &
                  ,dFc_daaa_ae, dFc_dbbb_ae &
                  ,dFc_dggg_ae, dFc_daab_ae &
                  ,dFc_daag_ae, dFc_dabb_ae &
                  ,dFc_dbbg_ae, dFc_dagg_ae &
                  ,dFc_dbgg_ae, dFc_dabg_ae, ista_nrc, iend_nrc, ist, ien, pot_type)

             call ex_gga_paw_library_3D(nrc,dnr,nspin &
                  ,nps_sph(ista_nrc,1,1), grad_nps &
                  ,exc_ps_field &
                  ,dF_dnps,    dF_dgradnps &
                  ,dFx_dnnps,  dFx_dngps   &
                  ,dFx_dggps,  dFx_dnnnps  &
                  ,dFx_dnngps, dFx_dnggps  &
                  ,dFx_dgggps, ista_nrc, iend_nrc, ist, ien, pot_type )
             call cr_gga_paw_library_3D(nrc,dnr,nspin &
                  ,nps_sph(ista_nrc,1,1), tmp_work, exc_ps_field &
                  ,dF_dnps, dF_dgradtnps &
                  ,dFc_daa_ps, dFc_dbb_ps, dFc_dgg_ps &
                  ,dFc_dab_ps, dFc_dag_ps, dFc_dbg_ps &
                  ,dFc_daaa_ps, dFc_dbbb_ps &
                  ,dFc_dggg_ps, dFc_daab_ps &
                  ,dFc_daag_ps, dFc_dabb_ps &
                  ,dFc_dbbg_ps, dFc_dagg_ps &
                  ,dFc_dbgg_ps, dFc_dabg_ps, ista_nrc, iend_nrc, ist, ien, pot_type)
             deallocate( tmp_work )

#ifdef LIBXC
          else if ( xctype == 'libxc' ) then
             select case (xc_family_exch)
             case (XC_FAMILY_LDA)
                call ex_lda_paw_libxc( nrc, dnr, nspin, &
                     &    nae_sph(ista_nrc,1,1), &
                     &    exc_ae_field, &
                     &    dF_dnae, dFx_dnnae, dFx_dnnnae, &
                     &    ista_nrc, iend_nrc, ist, ien )
                call ex_lda_paw_libxc( nrc, dnr, nspin, &
                     &    nps_sph(ista_nrc,1,1), &
                     &    exc_ps_field, &
                     &    dF_dnps, dFx_dnnps, dFx_dnnnps, &
                     &    ista_nrc, iend_nrc, ist, ien )
             case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
                call ex_gga_paw_libxc( nrc, dnr, nspin, &
                     &    nae_sph(ista_nrc,1,1), grad_nae, grad_tnae, &
                     &    exc_ae_field, &
                     &    dF_dnae,    dF_dgradnae, &
                     &    dFx_dnnae,  dFx_dngae,   dFx_dggae, &
                     &    dFx_dnnnae, dFx_dnngae,  dFx_dnggae,  dFx_dgggae, &
                     &    ista_nrc, iend_nrc, ist, ien )
                call ex_gga_paw_libxc( nrc, dnr, nspin, &
                     &    nps_sph(ista_nrc,1,1), grad_nps, grad_tnps, &
                     &    exc_ps_field, &
                     &    dF_dnps,    dF_dgradnps, &
                     &    dFx_dnnps,  dFx_dngps,   dFx_dggps, &
                     &    dFx_dnnnps, dFx_dnngps,  dFx_dnggps,  dFx_dgggps, &
                     &    ista_nrc, iend_nrc, ist, ien )
             end select

             select case (xc_family_corr)
             case (XC_FAMILY_LDA)
                call cr_lda_paw_libxc( nrc, dnr, nspin, &
                     &    nae_sph(ista_nrc,1,1), &
                     &    exc_ae_field, &
                     &    dF_dnae,    &
                     &    dFc_daa_ae, dFc_dbb_ae, &
                     &    dFc_dab_ae, &
                     &    dFc_daaa_ae, dFc_dbbb_ae, dFc_daab_ae, &
                     &    dFc_dabb_ae, &
                     &    ista_nrc, iend_nrc, ist, ien )
                call cr_lda_paw_libxc( nrc, dnr, nspin, &
                     &    nps_sph(ista_nrc,1,1), &
                     &    exc_ps_field, &
                     &    dF_dnps,    &
                     &    dFc_daa_ps, dFc_dbb_ps, &
                     &    dFc_dab_ps, &
                     &    dFc_daaa_ps, dFc_dbbb_ps, dFc_daab_ps, &
                     &    dFc_dabb_ps, &
                     &    ista_nrc, iend_nrc, ist, ien )
             case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
                call cr_gga_paw_libxc( nrc, dnr, nspin, &
                     &    nae_sph(ista_nrc,1,1), grad_nae, grad_tnae, &
                     &    exc_ae_field, &
                     &    dF_dnae,    dF_dgradtnae, &
                     &    dFc_daa_ae, dFc_dbb_ae, dFc_dgg_ae, &
                     &    dFc_dab_ae, dFc_dag_ae, dFc_dbg_ae, &
                     &    dFc_daaa_ae, dFc_dbbb_ae, dFc_dggg_ae, dFc_daab_ae, &
                     &    dFc_daag_ae, dFc_dabb_ae, dFc_dbbg_ae, dFc_dagg_ae, &
                     &    dFc_dbgg_ae, dFc_dabg_ae, &
                     &    ista_nrc, iend_nrc, ist, ien )
                call cr_gga_paw_libxc( nrc, dnr, nspin, &
                     &    nps_sph(ista_nrc,1,1), grad_nps, grad_tnps, &
                     &    exc_ps_field, &
                     &    dF_dnps,    dF_dgradtnps, &
                     &    dFc_daa_ps, dFc_dbb_ps, dFc_dgg_ps, &
                     &    dFc_dab_ps, dFc_dag_ps, dFc_dbg_ps, &
                     &    dFc_daaa_ps, dFc_dbbb_ps, dFc_dggg_ps, dFc_daab_ps, &
                     &    dFc_daag_ps, dFc_dabb_ps, dFc_dbbg_ps, dFc_dagg_ps, &
                     &    dFc_dbgg_ps, dFc_dabg_ps, &
                     &    ista_nrc, iend_nrc, ist, ien )
             end select
#endif
          else
             write(nfout,'(" xctype = ",a7)') xctype
             call phase_error_with_msg(nfout,' xctype is not set properly (ggaxcp_paw_sphex)',__LINE__,__FILE__)
         end if
                 !   dF/d|rho(r)| (=vxc) --> dF_drho
                 !   dFx/d|grad(rho(r))| --> dF_dgradrho
                 !   dFc/d|grad(rho(r))| --> grad_trho
        end subroutine ggaxcp_paw_sphex2_3D
! ==============================================================================

    subroutine wd_grad_etc(nelm)
!   ---- 2020.01.23 T. Yamasaki --
      integer, intent(in) :: nelm
      integer :: ij
      ij = nelm
      if(iend_nrc-ista_nrc+1<2*ij) ij = (iend_nrc-ista_nrc+1)/2
      if(ij == 0) ij = 1

      write(nfout,'(" -- grad_tnae(",i6,":",i6,") --")') ista_nrc,ista_nrc+ij-1
      write(nfout,'(8d16.4)') grad_tnae(ista_nrc:ista_nrc+ij-1)
      write(nfout,'(" -- grad_tnae(",i6,":",i6,") --")') iend_nrc-ij+1,iend_nrc
      write(nfout,'(8d16.4)') grad_tnae(iend_nrc-ij+1:iend_nrc)
      write(nfout,'(" -- dFc_dbg_ae(",i6,":",i6,") --")') ista_nrc,ista_nrc+ij-1
      write(nfout,'(8d16.4)') dFc_dbg_ae(ista_nrc:ista_nrc+ij-1)
      write(nfout,'(" -- dFc_dbg_ae(",i6,":",i6,") --")') iend_nrc-ij+1,iend_nrc
      write(nfout,'(8d16.4)') dFc_dbg_ae(iend_nrc-ij+1:iend_nrc)
    end subroutine wd_grad_etc

    subroutine wd_dnae_dr_sph(nelm)
!   ---- 2020.01.23 T. Yamasaki --
      integer, intent(in) :: nelm
      integer :: ij, is
      ij = nelm
      if(iend_nrc-ista_nrc+1<2*ij) ij = (iend_nrc-ista_nrc+1)/2
      if(ij == 0) ij = 1

      do is = 1, nspin
         write(nfout,'(" -- dnae_dr_sph(",i6,":",i6,",",i2,",1) --")') ista_nrc,ista_nrc+ij-1,is
         write(nfout,'(8d16.4)') dnae_dr_sph(ista_nrc:ista_nrc+ij-1,is,1)
         write(nfout,'(" -- dnae_dr_sph(",i6,":",i6,",",i2,",1) --")') iend_nrc-ij+1,iend_nrc,is
         write(nfout,'(8d16.4)') dnae_dr_sph(iend_nrc-ij+1:iend_nrc,is,1)
      end do
    end subroutine wd_dnae_dr_sph

    subroutine wd_dnps_dr_sph_etc(nelm)
!   ---- 2020.01.23 T. Yamasaki --
      integer, intent(in) :: nelm
      integer :: ij, is
      ij = nelm
      if(iend_nrc-ista_nrc+1<2*ij) ij = (iend_nrc-ista_nrc+1)/2
      if(ij == 0) ij = 1

      do is = 1, nspin
         write(nfout,'(" -- dnps_dr_sph(",i6,":",i6,",",i2,",1) --")') ista_nrc,ista_nrc+ij-1,is
         write(nfout,'(8d16.4)') dnps_dr_sph(ista_nrc:ista_nrc+ij-1,is,1)
         write(nfout,'(" -- dnps_dr_sph(",i6,":",i6,",",i2,",1) --")') iend_nrc-ij+1,iend_nrc,is
         write(nfout,'(8d16.4)') dnps_dr_sph(iend_nrc-ij+1:iend_nrc,is,1)
      end do
      write(nfout,'(" -- grad_tnps2_sph(",i6,":",i6,",1) --")') ista_nrc,ista_nrc+ij-1
      write(nfout,'(8d16.4)') grad_tnps2_sph(ista_nrc:ista_nrc+ij-1,1)
      write(nfout,'(" -- grad_tnps2_sph(",i6,":",i6,",1) --")') iend_nrc-ij+1,iend_nrc
      write(nfout,'(8d16.4)') grad_tnps2_sph(iend_nrc-ij+1:iend_nrc,1)
    end subroutine wd_dnps_dr_sph_etc

    subroutine wd_exc_etc(nelm)
!   ---- 2020.01.23 T. Yamasaki --
      integer, intent(in) :: nelm
      integer :: ij, is
      ij = nelm
      if(iend_nrc-ista_nrc+1<2*ij) ij = (iend_nrc-ista_nrc+1)/2
      if(ij == 0) ij = 1

      write(nfout,'(" ia = ", i8)') ia
      write(nfout,'(" -- exc_ae_field(",i6,":",i6,") -- ")') ista_nrc,ista_nrc+ij-1
      write(nfout,'(8d16.4)') exc_ae_field(ista_nrc:ista_nrc+ij-1)
      write(nfout,'(" -- exc_ae_field(",i6,":",i6,") -- ")') iend_nrc-ij+1,iend_nrc
      write(nfout,'(8d16.4)') exc_ae_field(iend_nrc-ij+1:iend_nrc)

      do is = 1, nspin
         write(nfout,'(" -- nae_sph(",i6,":",i6,",",i2,",1) -- ")') ista_nrc,ista_nrc+ij-1,is
         write(nfout,'(8d16.4)') nae_sph(ista_nrc:ista_nrc+ij-1, is, 1)
         write(nfout,'(" -- nae_sph(",i6,":",i6,",",i2,",1) -- ")') iend_nrc-ij+1,iend_nrc,is
         write(nfout,'(8d16.4)') nae_sph(iend_nrc-ij+1:iend_nrc, is, 1)
      end do

      do is = 1, nspin
         write(nfout,'(" -- dFx_dggae(",i6,":",i6,",",i2,") --")') ista_nrc,ista_nrc+ij-1,is
         write(nfout,'(8d16.4)') dFx_dggae(ista_nrc:ista_nrc+ij-1,is)
         write(nfout,'(" -- dFx_dggae(",i6,":",i6,",",i2,") --")') iend_nrc-ij+1,iend_nrc,is
         write(nfout,'(8d16.4)') dFx_dggae(iend_nrc-ij+1:iend_nrc,is)
      end do
    end subroutine wd_exc_etc

    subroutine wd_exc_ps_etc(nelm)
!   ---- 2020.01.23 T. Yamasaki --
      integer, intent(in) :: nelm
      integer :: ij, is
      ij = nelm
      if(iend_nrc-ista_nrc+1<2*ij) ij = (iend_nrc-ista_nrc+1)/2
      if(ij == 0) ij = 1

      write(nfout,'(" -- exc_ps_field(",i6,":",i6,") -- ")') ista_nrc,ista_nrc+ij-1
      write(nfout,'(8d16.4)') exc_ps_field(ista_nrc:ista_nrc+ij-1)
      write(nfout,'(" -- exc_ps_field(",i6,":",i6,") -- ")') iend_nrc-ij+1,iend_nrc
      write(nfout,'(8d16.4)') exc_ps_field(iend_nrc-ij+1:iend_nrc)

      do is = 1, nspin
         write(nfout,'(" -- dFx_dggps(",i6,":",i6,",",i2,") --")') ista_nrc,ista_nrc+ij-1,is
         write(nfout,'(8d16.4)') dFx_dggps(ista_nrc:ista_nrc+ij-1,is)
         write(nfout,'(" -- dFx_dggps(",i6,":",i6,",",i2,") --")') iend_nrc-ij+1,iend_nrc,is
         write(nfout,'(8d16.4)') dFx_dggps(iend_nrc-ij+1:iend_nrc,is)
      end do
    end subroutine wd_exc_ps_etc

    subroutine wd_exc_ps_etc2(nelm)
!   ---- 2020.01.23 T. Yamasaki --
      integer, intent(in) :: nelm
      integer :: ij, is
      ij = nelm
      if(iend_nrc-ista_nrc+1<2*ij) ij = (iend_nrc-ista_nrc+1)/2
      if(ij == 0) ij = 1

      write(nfout,'(" -- exc_ps_field(",i6,":",i6,") -- ")') ista_nrc,ista_nrc+ij-1
      write(nfout,'(8d16.4)') exc_ps_field(ista_nrc:ista_nrc+ij-1)
      write(nfout,'(" -- exc_ps_field(",i6,":",i6,") -- ")') iend_nrc-ij+1,iend_nrc
      write(nfout,'(8d16.4)') exc_ps_field(iend_nrc-ij+1:iend_nrc)

      write(nfout,'(" -- grad_tnps(",i6,":",i6,") --")') ista_nrc,ista_nrc+ij-1
      write(nfout,'(8d16.4)') grad_tnps(ista_nrc:ista_nrc+ij-1)
      write(nfout,'(" -- grad_tnps(",i6,":",i6,") --")') iend_nrc-ij+1,iend_nrc
      write(nfout,'(8d16.4)') grad_tnps(iend_nrc-ij+1:iend_nrc)

      write(nfout,'(" -- dFc_dbg_ps(",i6,":",i6,") --")') ista_nrc,ista_nrc+ij-1
      write(nfout,'(8d16.4)') dFc_dbg_ps(ista_nrc:ista_nrc+ij-1)
      write(nfout,'(" -- dFc_dbg_ps(",i6,":",i6,") --")') iend_nrc-ij+1,iend_nrc
      write(nfout,'(8d16.4)') dFc_dbg_ps(iend_nrc-ij+1:iend_nrc)
    end subroutine wd_exc_ps_etc2

! === For nrc decomposion. by takto 2012/12/05 =================================
        subroutine abs_grad_rho_ud_paw_sphex2_3D
      integer:: ir,is,ier,isp,nsp

      real(DP), allocatable, dimension(:,:,:) :: tmp ! === For nrc decomposion. by takto 2012/12/07 ===

      dnae_dr_sph = 0.0d0

      allocate(tmp(ista_nrc-6*dnr:iend_nrc+6*dnr,nspin,25))     ! == For nrc decomposion. by takto 2012/12/07 ==
      call boundary_exchange_dim3(tmp, nae_sph, dnr, nspin, 25) ! == For nrc decomposion ==

    ! ------ for AE CD -------
      if(vflag == VXC_AND_EXC) then
         do nsp=1,num_isph_chg
            isp=isph_chg(nsp)
            do is=1,nspin
               call calc_ddiff_exp3_3D(ier,3,nrc,dnr,radr_paw(:,it)   &
                    &  , tmp(ista_nrc-6*dnr,is,isp), dnae_dr_sph(ista_nrc,is,isp)  &
                    &  , ddnae_ddr_sph(ista_nrc,is,isp), ista_nrc, iend_nrc, ist, ien)
            end do
         end do
      else
         do nsp=1,num_isph_chg
            isp=isph_chg(nsp)
            do is=1,nspin
               call calc_diff_exp3_3D(ier,3,nrc,dnr,radr_paw(1,it) &
                    &  , tmp(ista_nrc-6*dnr,is,isp),  dnae_dr_sph(ista_nrc,is,isp) &
                    &  , ista_nrc, iend_nrc, ist, ien)
            end do
         end do
      end if

      deallocate(tmp) ! === For nrc decomposion. by takto 2012/12/07 ===

      call m_PAWCD_set_sq_der_cd_sdphex2(ia,nspin,nrc,dnr &
           &   , num_isph_chg,isph_chg, nae_sph, dnae_dr_sph &
           &   , grad_nae2_sph, grad_tnae2_sph      &
           &   , msphmx_grd,num_isph_grd,isph_grd,ista_nrc,iend_nrc,ist)
      if(iprixc >= DEBUGPRINTLEVEL) call wd_dnae_dr_sph(20)

      do is=1,nspin
         do ir=ist,ien,dnr
            grad_nae(ir,is)=sqrt(grad_nae2_sph(ir,is,1))
         end do
      end do
      do ir=ist,ien,dnr
         grad_tnae(ir)=sqrt(grad_tnae2_sph(ir,1))
      end do

      dnps_dr_sph = 0.0d0

      allocate(tmp(ista_nrc-6*dnr:iend_nrc+6*dnr,nspin,25))     ! == For nrc decomposion. by takto 2012/12/07 ==
      call boundary_exchange_dim3(tmp, nps_sph, dnr, nspin, 25) ! == For nrc decomposion. ==

    !    ! ------ for PS CD -------
      if(vflag == VXC_AND_EXC) then
         do nsp=1,num_isph_chg
            isp=isph_chg(nsp)
            do is=1,nspin
               call calc_ddiff_exp3_3D(ier,3,nrc,dnr,radr_paw(:,it)   &
                    &      , tmp(ista_nrc-6*dnr,is,isp)    &
                    &      , dnps_dr_sph(ista_nrc,is,isp), ddnps_ddr_sph(ista_nrc,is,isp) & ! -> dnps_dr_sph, ddnps_ddr_sph
                    &      , ista_nrc, iend_nrc, ist, ien)
            end do
         end do
      else
         do nsp=1,num_isph_chg
            isp=isph_chg(nsp)
            do is=1,nspin
               call calc_diff_exp3_3D(ier,3,nrc,dnr,radr_paw(:,it) &
                    &      , tmp(ista_nrc-6*dnr,is,isp)   &
                    &      , dnps_dr_sph(ista_nrc,is,isp) &
                    &      , ista_nrc, iend_nrc, ist, ien)
            end do
         end do
      end if

      deallocate(tmp) ! == For nrc decomposion ==

      call m_PAWCD_set_sq_der_cd_sdphex2(ia,nspin,nrc,dnr &
           &     , num_isph_chg,isph_chg, nps_sph             &
           &     , dnps_dr_sph, grad_nps2_sph &
           &     , grad_tnps2_sph, msphmx_grd &
           &     , num_isph_grd,isph_grd,ista_nrc,iend_nrc,ist)

      if(iprixc >= DEBUGPRINTLEVEL) call wd_dnps_dr_sph_etc(20)

      do is=1,nspin
         do ir=ist,ien,dnr
            grad_nps(ir,is)=sqrt(grad_nps2_sph(ir,is,1))
         end do
      end do
      do ir=ist,ien,dnr
         grad_tnps(ir)=sqrt(grad_tnps2_sph(ir,1))
      end do

      return
        end subroutine abs_grad_rho_ud_paw_sphex2_3D
! ==============================================================================


    end subroutine m_PAW_XC_cal_potential_sphex2

! === For nrc decomposion. by takto 2012/12/05 =================================
    subroutine set_sphex_elements2_3D(nrc,dnr,mode &
                                    ,msphmx_chg,num_isph_chg,isph_chg &
                                    ,msphmx_grd,num_isph_grd,isph_grd &
                                    ,num_isph_n_n,isph_n_n &
                                    ,num_isph_n_g,isph_n_g &
                                    ,num_isph_g_g,isph_g_g,ist,ien)
        integer,intent(in):: nrc,dnr,mode
        integer,intent(in):: msphmx_chg,num_isph_chg,isph_chg(25)
        integer,intent(in):: msphmx_grd,num_isph_grd,isph_grd(25)
        integer,intent(out):: num_isph_n_n,isph_n_n(25)
        integer,intent(out):: num_isph_n_g,isph_n_g(25)
        integer,intent(out):: num_isph_g_g,isph_g_g(25)
        integer,intent(in):: ist, ien
        integer:: itmp

            call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                       ,nae_sph(:,1,1:25)                      &
                                       ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                       ,nae_sph(:,1,1:25)                      &
                                       ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                       ,nana_ae_sph                            &
                                       ,itmp,num_isph_n_n,isph_n_n(1:25),ist,ien)
            call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                       ,nps_sph(:,1,1:25)                      &
                                       ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                       ,nps_sph(:,1,1:25)                      &
                                       ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                       ,nana_ps_sph                            &
                                       ,itmp,num_isph_n_n,isph_n_n(1:25),ist,ien)
            if(nspin.eq.2) then
                call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                           ,nae_sph(:,2,1:25)                      &
                                           ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                           ,nae_sph(:,2,1:25)                      &
                                           ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                           ,nbnb_ae_sph                            &
                                           ,itmp,num_isph_n_n,isph_n_n(1:25),ist,ien)
                call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                           ,nae_sph(:,1,1:25)                      &
                                           ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                           ,nae_sph(:,2,1:25)                      &
                                           ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                           ,nanb_ae_sph                            &
                                           ,itmp,num_isph_n_n,isph_n_n(1:25),ist,ien)
                call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                           ,nps_sph(:,2,1:25)                      &
                                           ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                           ,nps_sph(:,2,1:25)                      &
                                           ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                           ,nbnb_ps_sph                            &
                                           ,itmp,num_isph_n_n,isph_n_n(1:25),ist,ien)
                call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                           ,nps_sph(:,1,1:25)                      &
                                           ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                           ,nps_sph(:,2,1:25)                      &
                                           ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                           ,nanb_ps_sph                            &
                                           ,itmp,num_isph_n_n,isph_n_n(1:25),ist,ien)
            end if

            if(check_of_xctype()==GGA) then
                call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                           ,grad_nae2_sph(:,1,1:25)                &
                                           ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                           ,grad_nae2_sph(:,1,1:25)                &
                                           ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                           ,gaga_ae_sph                            &
                                           ,itmp,num_isph_g_g,isph_g_g(1:25),ist,ien)
                call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                           ,nae_sph(:,1,1:25)                      &
                                           ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                           ,grad_nae2_sph(:,1,1:25)                &
                                           ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                           ,naga_ae_sph                            &
                                           ,itmp,num_isph_n_g,isph_n_g(1:25),ist,ien)
                call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                           ,grad_tnae2_sph                         &
                                           ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                           ,grad_tnae2_sph                         &
                                           ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                           ,gg_ae_sph                              &
                                           ,itmp,num_isph_g_g,isph_g_g(1:25),ist,ien)
                call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                           ,nae_sph(:,1,1:25)                      &
                                           ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                           ,grad_tnae2_sph                         &
                                           ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                           ,nag_ae_sph                             &
                                           ,itmp,num_isph_n_g,isph_n_g(1:25),ist,ien)
                call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                           ,grad_nps2_sph(:,1,1:25)                &
                                           ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                           ,grad_nps2_sph(:,1,1:25)                &
                                           ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                           ,gaga_ps_sph                            &
                                           ,itmp,num_isph_g_g,isph_g_g(1:25),ist,ien)
                call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                           ,nps_sph(:,1,1:25)                      &
                                           ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                           ,grad_nps2_sph(:,1,1:25)                &
                                           ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                           ,naga_ps_sph                            &
                                           ,itmp,num_isph_n_g,isph_n_g(1:25),ist,ien)
                call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                           ,grad_tnps2_sph                         &
                                           ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                           ,grad_tnps2_sph                         &
                                           ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                           ,gg_ps_sph                              &
                                           ,itmp,num_isph_g_g,isph_g_g(1:25),ist,ien)
                call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                           ,nps_sph(:,1,1:25)                      &
                                           ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                           ,grad_tnps2_sph                         &
                                           ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                           ,nag_ps_sph                             &
                                           ,itmp,num_isph_n_g,isph_n_g(1:25),ist,ien)
                if(nspin.eq.2) then
                    call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                               ,grad_nae2_sph(:,2,1:25)                &
                                               ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                               ,grad_nae2_sph(:,2,1:25)                &
                                               ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                               ,gbgb_ae_sph                            &
                                               ,itmp,num_isph_g_g,isph_g_g(1:25),ist,ien)
                    call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                               ,nae_sph(:,2,1:25)                      &
                                               ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                               ,grad_nae2_sph(:,2,1:25)                &
                                               ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                               ,nbgb_ae_sph                            &
                                               ,itmp,num_isph_n_g,isph_n_g(1:25),ist,ien)
                    call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                               ,nae_sph(:,2,1:25)                      &
                                               ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                               ,grad_tnae2_sph                         &
                                               ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                               ,nbg_ae_sph                             &
                                               ,itmp,num_isph_n_g,isph_n_g(1:25),ist,ien)
                    call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                               ,grad_nps2_sph(:,2,1:25)                &
                                               ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                               ,grad_nps2_sph(:,2,1:25)                &
                                               ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                               ,gbgb_ps_sph                            &
                                               ,itmp,num_isph_g_g,isph_g_g(1:25),ist,ien)
                    call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                               ,nps_sph(:,2,1:25)                      &
                                               ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                               ,grad_nps2_sph(:,2,1:25)                &
                                               ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                               ,nbgb_ps_sph                            &
                                               ,itmp,num_isph_n_g,isph_n_g(1:25),ist,ien)
                    call mult_sphex_element3_3D(nrc,dnr,mode                           &
                                               ,nps_sph(:,2,1:25)                      &
                                               ,msphmx_chg,num_isph_chg,isph_chg(1:25) &
                                               ,grad_tnps2_sph                         &
                                               ,msphmx_grd,num_isph_grd,isph_grd(1:25) &
                                               ,nbg_ps_sph                             &
                                               ,itmp,num_isph_n_g,isph_n_g(1:25),ist,ien)

                end if
            end if
            return

    end subroutine set_sphex_elements2_3D
! ==============================================================================

! === For nrc decomposion. by takto 2012/12/05 =================================
    subroutine mult_sphex_element_full3_3D(nrc,dnr,mode  &
                                            ,n1,msphmx1,num_isph1,isph1 &
                                            ,n2,msphmx2,num_isph2,isph2 &
                                            ,n3,msphmx3,num_isph3,isph3 &
                                            ,ist,ien)
        integer,intent(in):: nrc,dnr,mode
        integer,intent(in):: msphmx1,num_isph1,isph1(25)
        integer,intent(in):: msphmx2,num_isph2,isph2(25)
        real(DP),intent(in):: n1(ista_nrc:iend_nrc,25),n2(ista_nrc:iend_nrc,25)
        real(DP),intent(out):: n3(ista_nrc:iend_nrc,25)
        integer,intent(out):: msphmx3,num_isph3,isph3(25)
        integer,intent(in) :: ist,ien

        integer:: ir,isp,isp1,isp2,n
        integer:: nsp,nsp1,nsp2
        real(DP):: fac,cijk
        logical:: flg_isp(25)

        n3=0.d0
!        do isp=2,msphmx
        do nsp=2,num_isph1
            isp=isph1(nsp)
            do ir=ist,ien,dnr
                n3(ir,1)=n3(ir,1) + n1(ir,isp)*n2(ir,isp)
            end do
        end do
        do ir=ist,ien,dnr
            n3(ir,1)=n3(ir,1)/PAI4+n1(ir,1)*n2(ir,1)
        end do

        if(mode.eq.0) return

        flg_isp=.false.
        flg_isp(1)=.true.

        msphmx3=0
!        do isp1=2,min(16,msphmx)                                    ! sphset2
!        do isp1=2,msphmx                                            ! sphset3
        do nsp1=2,num_isph1
            isp1=isph1(nsp1)
!            do isp2=isp1,min(16,msphmx)                             ! sphset2
!            do isp2=isp1,msphmx                                     ! sphset3
            do nsp2=nsp1,num_isph2
                isp2=isph2(nsp2)
                fac=1.d0;if(isp1.eq.isp2) fac=0.5d0
                do n=1,paw_mmt2(isp1,isp2)
                    isp=paw_isph2(isp1,isp2,n)
                    if(isp.eq.1) cycle
                    if(isp.gt.msphmx3) msphmx3=isp
                    flg_isp(isp)=.true.
                    cijk=paw_cr2(isp1,isp2,n)
! print *,isp,isp1,isp2
                    do ir=ist,ien,dnr
                        n3(ir,isp)=n3(ir,isp) + &
                            fac*(n1(ir,isp1)*n2(ir,isp2)+n1(ir,isp2)*n2(ir,isp1))*cijk
                    end do
                end do
            end do
        end do
!        do isp=2,msphmx
        do nsp=2,num_isph1
            isp=isph1(nsp)
            do ir=ist,ien,dnr
                n3(ir,isp)=n3(ir,isp) + n1(ir,1)*n2(ir,isp) + n1(ir,isp)*n2(ir,1)
            end do
        end do

        num_isph3=0
        isph3=0
        do isp=1,25
            if(flg_isp(isp)) then
                num_isph3=num_isph3+1
                isph3(num_isph3)=isp
            end if
        end do

        return

    end subroutine mult_sphex_element_full3_3D
! ==============================================================================

! === For nrc decomposion. by takto 2012/12/05 =================================
    subroutine get_paw_sphex_integral3_3D(nrc,dnr,rd  &
                                            ,n1,msphmx1,num_isph1,isph1 &
                                            ,n2,msphmx2,num_isph2,isph2 &
                                            ,n3,msphmx3,num_isph3,isph3,ist,ien)
        integer,intent(in):: nrc,dnr
        integer,intent(in):: msphmx1,num_isph1,isph1(25)
        integer,intent(in):: msphmx2,num_isph2,isph2(25)
        real(DP),intent(in):: n1(ista_nrc:iend_nrc,25),n2(ista_nrc:iend_nrc,25),rd(nrc)
        real(DP),intent(out):: n3(ista_nrc:iend_nrc,25)
        integer,intent(out):: msphmx3,num_isph3,isph3(25)
        integer,intent(in) :: ist,ien

        integer:: isp,isp1,isp2,n,ir
        integer:: nsp,nsp1,nsp2
        real(DP):: dl,dl1,dl2,fac,cijk
        integer,allocatable,dimension(:):: il3
        logical:: flg_isp(25)

        allocate(il3(25));call substitute_il3(25,il3)

        n3=0.d0
        flg_isp=.false.
!        do isp1=2,min(16,msphmx)                                        ! sphset2
!        do isp1=2,msphmx                                                ! sphset3
        do nsp1=2,num_isph1
            isp1=isph1(nsp1)
            dl1=dble(il3(isp1))
            dl1=dl1*(dl1+1.d0)
!            do isp2=isp1,min(16,msphmx)                                 ! sphset2
!            do isp2=isp1,msphmx                                         ! sphset3
            do nsp2=nsp1,num_isph2
                isp2=isph2(nsp2)
                dl2=dble(il3(isp2))
                dl2=dl2*(dl2+1.d0)
                fac=0.5d0;if(isp1.eq.isp2) fac=0.25d0
!                fac=0.5d0
                do n=1,paw_mmt2(isp1,isp2)
                    isp=paw_isph2(isp1,isp2,n)
                    if(isp.eq.1) cycle
                    if(isp.gt.msphmx3) msphmx3=isp
                    flg_isp(isp)=.true.
                    cijk=paw_cr2(isp1,isp2,n)
                    dl=dble(il3(isp))
                    dl=dl*(dl+1.d0)

                    do ir=ist,ien,dnr
                        n3(ir,isp)=n3(ir,isp) + &
                            fac*( &
                            n1(ir,isp1)*n2(ir,isp2)*(dl2+dl-dl1) + &
                            n1(ir,isp2)*n2(ir,isp1)*(dl1+dl-dl2))*cijk
                    end do

                end do
            end do
        end do

        num_isph3=0
        isph3=0
        do isp=1,25
            if(flg_isp(isp)) then
                num_isph3=num_isph3+1
                isph3(num_isph3)=isp
            end if
        end do

        do isp=2,25
            dl=dble(il3(isp))
            dl=dl*(dl+1.d0)
            do ir=ist,ien,dnr
                n3(ir,isp)=(dl*n1(ir,1)*n2(ir,isp)+n3(ir,isp))/rd(ir)/rd(ir)
            end do
        end do
        return

    end subroutine get_paw_sphex_integral3_3D
! ==============================================================================

    subroutine merge_isph_flgs(n1,isph1,n2,isph2,n3,isph3 &
                                ,n4,isph4,n5,isph5,n6,isph6,n7,isph7)
        integer,intent(in):: n1,isph1(25)
        integer,intent(in):: n2,isph2(25)
        integer,intent(in):: n3,isph3(25)
        integer,intent(in):: n4,isph4(25)
        integer,intent(in):: n5,isph5(25)
        integer,intent(in):: n6,isph6(25)
        integer,intent(out):: n7,isph7(25)

        integer:: nsp
        logical:: flg(25)

        flg=.false.
        do nsp=1,n1
            flg(isph1(nsp))=.true.
        end do
        do nsp=1,n2
            flg(isph2(nsp))=.true.
        end do
        do nsp=1,n3
            flg(isph3(nsp))=.true.
        end do
        do nsp=1,n4
            flg(isph4(nsp))=.true.
        end do
        do nsp=1,n5
            flg(isph5(nsp))=.true.
        end do
        do nsp=1,n6
            flg(isph6(nsp))=.true.
        end do
        n7=0
        isph7=0
        do nsp=1,25
            if(flg(nsp)) then
                n7=n7+1
                isph7(n7)=nsp
            end if
        end do
        return
    end subroutine merge_isph_flgs

    subroutine xcpotf_paw(nrc,ispin,ith,input_charge)
        integer, intent(in) :: nrc,ispin,ith,input_charge

! #1) 1994/11/08 by T.Yamasaki
!    Coding for the case of xctype='PERZUN ' and  'XALPHA ' are done.
! #2) Spin-polarization is introduced by T. Yamasaki at 15th Dec. 1994
! #3) f77 -> f90     4th April 1999  by T. Yamasaki

    real(kind=DP) :: DELTA
    data DELTA/1.d-40/

    if(xctype == 'wign   '.or. xctype == 'wigner ')then
        call xcpotf_wigner(ispin,1,nrc,input_charge,DELTA &
                                ,nae(1:nrc,1:ispin),wos(1:nrc),texc)
        exc_ae=exc_ae+texc*omg_wght(ith)
        call xcpotf_wigner(ispin,1,nrc,input_charge,DELTA &
                                ,nps(1:nrc,1:ispin),wos(1:nrc),texc)
        exc_ps=exc_ps+texc*omg_wght(ith)
    else if(  xctype ==  'pzold  ') then
        call xcpotf_pzold(ispin,1,nrc,input_charge,DELTA &
                                ,nae(1:nrc,1:nspin),wos(1:nrc),texc)
        exc_ae=exc_ae+texc*omg_wght(ith)
        call xcpotf_pzold(ispin,1,nrc,input_charge,DELTA &
                                ,nps(1:nrc,1:nspin),wos(1:nrc),texc)
        exc_ps=exc_ps+texc*omg_wght(ith)
    else if(  xctype ==  'xalfa  ') then
       call xcpotf_xalfa(ispin,1,nrc,input_charge,DELTA &
                                ,nae(1:nrc,1:nspin),wos(1:nrc),texc)
        exc_ae=exc_ae+texc*omg_wght(ith)
       call xcpotf_xalfa(ispin,1,nrc,input_charge,DELTA &
                                ,nps(1:nrc,1:nspin),wos(1:nrc),texc)
        exc_ps=exc_ps+texc*omg_wght(ith)
    else if(  xctype == 'perzun '.or. xctype == 'pz     ') then
       call xcpotf_pz(ispin,1,nrc,input_charge,DELTA &
                                ,nae(1:nrc,1:nspin),wos(1:nrc),texc)
        exc_ae=exc_ae+texc*omg_wght(ith)
       call xcpotf_pz(ispin,1,nrc,input_charge,DELTA &
                                ,nps(1:nrc,1:nspin),wos(1:nrc),texc)
        exc_ps=exc_ps+texc*omg_wght(ith)
    else if(  xctype == 'vwn    ') then
       call xcpotf_vwn(ispin,1,nrc,input_charge,DELTA &
                                ,nae(1:nrc,1:nspin),wos(1:nrc),texc)
        exc_ae=exc_ae+texc*omg_wght(ith)
       call xcpotf_vwn(ispin,1,nrc,input_charge,DELTA &
                                ,nps(1:nrc,1:nspin),wos(1:nrc),texc)
        exc_ps=exc_ps+texc*omg_wght(ith)
    else if(  xctype=='mjw    '.or. xctype=='bh     ' .or.xctype=='gl     ') then
       call xcpotf_mjw_bh_gl(len_xctype,xctype,ispin,1,nrc,input_charge,DELTA &
                                ,nae(1:nrc,1:nspin),wos(1:nrc),texc)
        exc_ae=exc_ae+texc*omg_wght(ith)
       call xcpotf_mjw_bh_gl(len_xctype,xctype,ispin,1,nrc,input_charge,DELTA &
                                ,nps(1:nrc,1:nspin),wos(1:nrc),texc)
        exc_ps=exc_ps+texc*omg_wght(ith)
    else
        write(*,'(" xctype = ",a7)') xctype
        call phase_error_with_msg(6,' xctype is not set properly (xcpotf_paw)',__LINE__,__FILE__)
    endif
                        ! all xcpotf_* subroutines are in -(b_XC_Potential) ->chgrhr_l,exc

    end subroutine xcpotf_paw

#if 0
    subroutine m_PAW_XC_get_dion_vxc(nfout)
        integer,intent(in):: nfout
        integer:: ia,it,lmt1,lmt2,is
        integer:: n,ilm3,l3,iiqitg
        integer:: ilt1,ilt2,il1,il2,ilk1
        real(DP):: fac,sum
        real(DP),allocatable,dimension(:,:,:,:,:):: pipjvxc_k

        allocate(pipjvxc_k(nltpw,nltpw,msph,nspin,natm))
        pipjvxc_k=0.d0
        call cnstrct_of_PiPjVxc_k(nltpw,msph,natm,nspin,pipjvxc_k)

        do ia=1,natm
            it=ityp(ia)
            if(ipaw(it)/=1) then
                dion_vxc(:,:,:,ia)=0.d0
                cycle
            end if

            do lmt1=1,ilmt(it)
                ilt1=index_lmt2lt(lmt1,it)
                do lmt2=lmt1,ilmt(it)
                    ilt2=index_lmt2lt(lmt2,it)
                    do is=1,nspin
                        sum=0.d0
                        do n=1,il2p(lmt1,lmt2,it)
                            ilm3=isph(lmt1,lmt2,n,it)
                            sum=sum+dl2p(lmt1,lmt2,n,it)* &
                                            pipjvxc_k(ilt1,ilt2,ilm3,is,ia)
                        end do
                        dion_vxc(lmt1,lmt2,is,ia)=sum
                    end do
                end do
            end do
        end do

        if(ipripp>=2.and.printable)then
        write(nfout,*)
        write(nfout,*) ' -- dion_vxc ---'
        do ia=1,natm
            it=ityp(ia)
            if(ipaw(it)/=1) cycle
            do is=1,nspin
                write(nfout,'(a,i2,a,i2,a)') '(ia,is)=(',ia,',',is,')'
                do lmt1 = 1, ilmt(it)
                    write(nfout,'(i3,15f12.9/15f12.9)') lmt1 &
                      &               ,(dion_vxc(lmt1,lmt2,is,ia),lmt2 = 1, ilmt(it))
                enddo
            end do
        end do
        endif

        deallocate(pipjvxc_k)
        return
    end subroutine m_PAW_XC_get_dion_vxc
#else
    subroutine m_PAW_XC_get_dion_vxc(nfout)
        integer,intent(in):: nfout
        integer:: ia,it,lmt1,lmt2,is
        integer:: n,ilm3,l3,iiqitg
        integer:: ilt1,ilt2,il1,il2,ilk1, ja
        real(DP):: fac,sum
        real(DP),allocatable,dimension(:,:,:,:,:):: pipjvxc_k

        allocate(pipjvxc_k(nltpw,nltpw,msph,nspin,natm))
        pipjvxc_k=0.d0

        if(af /= 0) allocate(flg_done(natm))

!!$        call cnstrct_of_PiPjVxc_k(nltpw,msph,natm,nspin,pipjvxc_k)
        call cnstrct_of_PiPjVxc_k2(nltpw,msph,natm,nspin,pipjvxc_k)

        if(af /= 0) flg_done=.false.

        do ia=1,natm
! === DEBUG by tkato 2011/10/01 ================================================
!           if(af /= 0 .and. flg_done(ia)) cycle
            if(af /= 0) then
               if(flg_done(ia)) then
                  cycle
               endif
            endif
! ==============================================================================
            it=ityp(ia)
            if(ipaw(it)/=1) then
                dion_vxc(:,:,:,ia)=0.d0
                cycle
            end if

            do lmt1=1,ilmt(it)
                ilt1=index_lmt2lt(lmt1,it)
                do lmt2=lmt1,ilmt(it)
                    ilt2=index_lmt2lt(lmt2,it)
                    do is=1,nspin
                        sum=0.d0
                        do n=1,il2p(lmt1,lmt2,it)
                            ilm3=isph(lmt1,lmt2,n,it)
                            sum=sum+dl2p(lmt1,lmt2,n,it)* &
                                            pipjvxc_k(ilt1,ilt2,ilm3,is,ia)
                        end do
                        dion_vxc(lmt1,lmt2,is,ia)=sum
                    end do
                end do
            end do

!ASMS            if(af /= 0) then
!ASMS                ja=ia2ia_symmtry_op(ia,nopr+af)
!ASMS                dion_vxc(:,:,2,ja)=dion_vxc(:,:,1,ia)
!ASMS                dion_vxc(:,:,1,ja)=dion_vxc(:,:,2,ia)
!ASMS                flg_done(ja)=.true.
!ASMS            end if

        end do

        if(ipripp>=2.and.printable)then
        write(nfout,*)
        write(nfout,*) ' -- dion_vxc ---'
        do ia=1,natm
            it=ityp(ia)
            if(ipaw(it)/=1) cycle
            do is=1,nspin
                write(nfout,'(a,i2,a,i2,a)') '(ia,is)=(',ia,',',is,')'
                do lmt1 = 1, ilmt(it)
                    write(nfout,'(i3,15f12.9/15f12.9)') lmt1 &
                      &               ,(dion_vxc(lmt1,lmt2,is,ia),lmt2 = 1, ilmt(it))
                enddo
            end do
        end do
        endif

        deallocate(pipjvxc_k)
        if(af /= 0) deallocate(flg_done)

        return
    end subroutine m_PAW_XC_get_dion_vxc
#endif


  subroutine cnstrct_of_PiPjVxc_k2(n,m,na,ns,mat)
    integer,intent(in):: n,m,na,ns
    real(DP),intent(out):: mat(n,n,m,ns,na)

    integer:: ia,it,is,ksph,nrc,ier,ir,dnr,nrc0
    integer:: ilt1,ilt2,il1,il2,it1,it2,ilk
    integer:: iiqitg
    integer:: ist,ien, ia_add, ip, ierr
    real(DP):: wos(mmesh),sum,zz
    real(DP),allocatable,dimension(:,:,:,:) :: matb

    mat = 0.d0
    allocate(matb(n,n,m,ns))
    do ia = ista_natm, iend_natm  !      do ia=1,natm
       ia_add = ia - ista_natm + 1

       do it=1,ntyp
          if(ityp(ia)/=it) cycle
          if(ipaw(it)/=1) then
             mat(:,:,:,:,ia)=0.d0
             cycle
          end if

! **** set dnr from input ****
          dnr=paw_dnr(it)
          if(dnr.gt.1) then
             nrc0=wf_mnrc(it)
             nrc=1+int((nrc0-1)/dnr)*dnr
!                    zz = (radr_paw(nrc0,it)-radr_paw(nrc,it))/ &
!                            (radr_paw(nrc+dnr,it)-radr_paw(nrc,it))
!                    zz = log(radr_paw(nrc0,it)/radr_paw(nrc,it))/ &
!                        log(radr_paw(nrc+dnr,it)/radr_paw(nrc,it))
             zz = dble(nrc0-nrc)/dble(dnr)
             nrc=nrc+2*dnr                                                  ! 3rd interpolation
             call set_weight_exp3(ier,1,nrc,dnr,radr_paw(:,it),zz,wos)      !  3rd
          else
             nrc=wf_mnrc(it)
             call set_weight_exp(ier,1,nrc,radr_paw(:,it),wos)
          end if

          ist = (ista_nrc+dnr-2)/dnr
          ist = dnr*ist+1
          ien = min(iend_nrc, nrc)

          if ( ist <= 0 ) ist = 1
#ifdef DEBUG_PAW3D
          write(6,'(" <<cnstrct_of_PiPjVxc_k2>> ia = ",i8)') ia
#endif
          do ilt1=1,iltpw(it)
             il1=lppw(ilt1,it)
             it1=tppw(ilt1,it)
             do ilt2=ilt1,iltpw(it)
                il2=lppw(ilt2,it)
                it2=tppw(ilt2,it)
                do ilk=abs(il1-il2),il1+il2-2,2
                   if(ilk > 4) cycle
                   iiqitg=iqitg(il1,it1,il2,it2,ilk+1,it)
                   if(iiqitg==0) then
                   end if
                   do ksph=ilk**2+1,ilk**2+2*ilk+1
                      do is=1,nspin
                         sum=0.d0
                         do ir=ist,ien,dnr   !!$   do ir=ista_nrc,iend_nrc
                            sum=sum &
                                 & + ( psirpw(ir,il1,it1,it)*psirpw(ir,il2,it2,it)* vxc_ae_k_3D(ir,is,ksph,ia_add) &
                                 &    - (phirpw(ir,il1,it1,it)*phirpw(ir,il2,it2,it)+ &
                                 &                             qrspspw(ir,iiqitg))*vxc_ps_k_3D(ir,is,ksph,ia_add) &
                                 &    ) * wos(ir)
                         end do
!!$                             call mpi_allreduce(MPI_IN_PLACE,sum,1,mpi_double_precision,mpi_sum,mpi_natm_world,ierr)
!!$                             mat(ilt1,ilt2,ksph,is,ia)=sum
                         matb(ilt1,ilt2,ksph,is)=sum
                      end do
                   end do
                end do
             end do
          end do
       end do
       call mpi_allreduce(MPI_IN_PLACE,matb,n*n*m*ns,mpi_double_precision,mpi_sum,mpi_natm_world,ierr)
       mat(:,:,:,:,ia) = matb(:,:,:,:)
    end do

    ip = 0
    do ia = 1, natm
       do ir = 0,nrank_natm-1
          if(is_natm(ir)<=ia .and. ia <=ie_natm(ir)) then
             ip = ir
             if(ip == myrank_natm) matb(:,:,:,:) = mat(:,:,:,:,ia)
             cycle
          end if
       end do
       call mpi_bcast(matb,n*n*m*ns,mpi_double_precision,ip,mpi_nrc_world,ierr)
       mat(:,:,:,:,ia) = matb(:,:,:,:)
    end do

    deallocate(matb)

    return
  end subroutine cnstrct_of_PiPjVxc_k2

!
!    subroutine m_PAW_XC_get_vxc_of_old_CD
!        integer:: lmt1,lmt2,is,ia,it
!        real(DP):: fac,sum
!        do ia=1,natm
!            do is=1,nspin
!                vxc_ae_m_ps(is,ia)=0.d0
!                do it=1,ntyp
!                    if(ityp(ia)/=it .or. .not.ipaw(it)) cycle
!
!                    do lmt1=1,ilmt(it)
!                        do lmt2=lmt1,ilmt(it)
!    !                        fac=2.d0*iwei(ia);if(lmt1.eq.lmt2) fac=iwei(ia)
!    !                            eho_paw(ia)=eho_paw(ia)+ &
!    !                                        fac*hsr(ia,lmt1,lmt2,is)* &
!    !                                        dion_hartree(lmt1,lmt2,ia)
!                        end do
!                    end do
!    !                eho_paw(ia)=eho_paw(ia)/2.d0
!                end do
!            end do
!        end do
!
!
!        return
!    end subroutine m_PAW_XC_get_vxc_of_old_CD

   subroutine calc_diff_exp2(isdiff,iord,n,xh,rn,fn,dfn,ddfn)
      implicit none
   integer,intent(in)  :: isdiff,iord, n
   real(8),intent(in)  :: xh,rn(*), fn(*)
!!$   integer,intent(out) :: ier
   real(8),intent(out) :: dfn(*)
   real(8),optional,intent(inout):: ddfn(*)
   real(8) :: r, f, df,ddf,h
   real(8),allocatable :: p(:,:), dp(:,:), ddp(:,:)
   integer :: n1, n2, i, pdim, indx,ier
!!$   ier = 0
!!$   isdiff = 1
   integer       :: id_sname = -1
   call tstatc0_begin('calc_diff_exp2 ',id_sname)
   pdim = 1+2*iord
   h = 1.d0/xh
   allocate(p(pdim,pdim),dp(pdim,pdim))
   if(isdiff >=2)  allocate(ddp(pdim,pdim))

   do i = 1,n
      r = rn(i)
      if (i < iord+1) then
         n1 = 1 ; n2 = 1 + 2*iord
      else if (i > n-iord) then
         n1 = n - 2*iord ; n2 = n
      else
         n1 = i - iord ; n2 = i + iord
      end if
      indx = i
      ! indx-n1 = i-(i-iord) = iord
      call diff_exp2(isdiff,n1,n2,rn,fn,r,df,ddf)
!!$      call diff_exp(ier,isdiff,n1,n2,rn,fn,r,f,df,ddf)
      dfn(i) = df
      if(isdiff==2) ddfn(i) = ddf
   end do
   if(isdiff >=2) deallocate(ddp)
   deallocate(p,dp)
   call tstatc0_end(id_sname)
 contains
   subroutine diff_exp2(isdiff,n1,n2,rn,fn,r,df,ddf)
!---------------------------------------------------------------------
!
!   Program written by Masakuni Okamoto
!  Revised for logarithmic grid by Takahiro Yamasaki, 2010/06/10
!
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: isdiff, n1, n2
   real(8),intent(in)  :: rn(*), fn(*), r
   real(8),intent(out) :: df,ddf
   real(8) :: denom_inv, a, b, c, x, x1, x2, dn2minusn1
   integer :: n, i, j, n1_delta
   n = n2-n1+1
   p(1:n,1) = fn(n1:n2) ;  dp(:,:) = 0.d0
   if(isdiff>=2) ddp(:,:) = 0.d0
   c = h*(n2-n1)
   n1_delta =indx-n1
   dn2minusn1 = 1.d0/dble(n2-n1)

   if ((isdiff >= 0).and.(n >= 2)) then
      do j = 2,n
         denom_inv = -dble(n-1)/dble(j-1)
         do i = 1,n+1-j
            p(i,j) = ((n1_delta-i-j+2)*p(i,j-1) + (i-1-n1_delta)*p(i+1,j-1))*dn2minusn1 * denom_inv
         end do
      end do
   end if
   if (isdiff >= 1) then
      do j = 2,n
         denom_inv = -dble(n-1)/dble(j-1)
         do i = 1,n+1-j
            dp(i,j) = ((n1_delta-i-j+2)*dp(i,j-1) + (i-1-n1_delta)*dp(i+1,j-1))*dn2minusn1 &
                 &     + (p(i,j-1)-p(i+1,j-1))
            dp(i,j) = dp(i,j) * denom_inv
         end do
      end do
   end if
   if (isdiff >= 2) then
      do j = 2,n
         denom_inv = -dble(n-1)/dble(j-1)
         do i = 1,n+1-j
            ddp(i,j) = ((n1_delta-i-j+2)*ddp(i,j-1) + (i-1-n1_delta)*ddp(i+1,j-1))*dn2minusn1 &
                     + 2.d0*(dp(i,j-1)-dp(i+1,j-1))
            ddp(i,j) = ddp(i,j) * denom_inv
         end do
      end do
   end if
   df  = dp (1,n) * (1.d0/r/c)
   if(isdiff >=2) ddf = ddp(1,n) * (1.d0/r/c)**2 - dp(1,n) * (1.d0/r/r/c)
 end subroutine diff_exp2

 end subroutine calc_diff_exp2
! === For nrc decomposion. by takto 2012/12/07 =================================
  subroutine boundary_exchange_dim2(array_out, array_in, dn, dim2)
  implicit none
!!  include 'mpif.h'
  integer, intent(in) :: dn, dim2
  real(kind=DP), intent(in)  :: array_in (ista_nrc     :iend_nrc,     dim2)
  real(kind=DP), intent(out) :: array_out(ista_nrc-6*dn:iend_nrc+6*dn,dim2)
  integer :: pn, pp, stat(MPI_STATUS_SIZE)

  integer :: i, srcount
  real(kind=DP), allocatable, dimension(:,:) :: sendbuf, recvbuf

  pn = myrank_nrc - 1
  pp = myrank_nrc + 1
  if(pn == -1)        pn = MPI_PROC_NULL
! if(pp == nrank_nrc) pp = MPI_PROC_NULL
  if(paw_last_rank_on_nrc) pp = MPI_PROC_NULL

  srcount = 6*dn*dim2

  allocate(sendbuf(6*dn,dim2))
  allocate(recvbuf(6*dn,dim2))

  array_out(ista_nrc:iend_nrc,:) = array_in(ista_nrc:iend_nrc,:)

  do i = 1, dim2
     sendbuf(1:6*dn,i) = array_in(iend_nrc-6*dn+1:iend_nrc,i)
  end do
  call MPI_Sendrecv(sendbuf, srcount, MPI_DOUBLE_PRECISION, pp, 1, &
                    recvbuf, srcount, MPI_DOUBLE_PRECISION, pn, 1, &
                    mpi_natm_world, stat, ierr)
  do i = 1, dim2
     array_out(ista_nrc-6*dn:ista_nrc-1,i) = recvbuf(1:6*dn,i)
  end do

  do i = 1, dim2
     sendbuf(1:6*dn,i) = array_in(ista_nrc:ista_nrc+6*dn-1,i)
  end do
  call MPI_Sendrecv(sendbuf, srcount, MPI_DOUBLE_PRECISION, pn, 1, &
                    recvbuf, srcount, MPI_DOUBLE_PRECISION, pp, 1, &
                    mpi_natm_world, stat, ierr)
  do i = 1, dim2
     array_out(iend_nrc+1:iend_nrc+6*dn,i) = recvbuf(1:6*dn,i)
  end do

  deallocate(sendbuf)
  deallocate(recvbuf)

  return
  end subroutine boundary_exchange_dim2

  subroutine boundary_exchange_dim3(array_out, array_in, dn, dim2, dim3)
  implicit none
!!  include 'mpif.h'
  integer, intent(in) :: dn, dim2, dim3
  real(kind=DP), intent(in)  :: array_in (ista_nrc     :iend_nrc,     dim2,dim3)
  real(kind=DP), intent(out) :: array_out(ista_nrc-6*dn:iend_nrc+6*dn,dim2,dim3)
  integer :: pn, pp, stat(MPI_STATUS_SIZE)

  integer :: i, j, srcount
  real(kind=DP), allocatable, dimension(:,:,:) :: sendbuf, recvbuf

  pn = myrank_nrc - 1
  pp = myrank_nrc + 1
  if(pn == -1)        pn = MPI_PROC_NULL
! if(pp == nrank_nrc) pp = MPI_PROC_NULL
  if(paw_last_rank_on_nrc) pp = MPI_PROC_NULL

  srcount = 6*dn*dim2*dim3

  allocate(sendbuf(6*dn,dim2,dim3))
  allocate(recvbuf(6*dn,dim2,dim3))

  array_out(ista_nrc:iend_nrc,:,:) = array_in(ista_nrc:iend_nrc,:,:)

  do j = 1, dim3
     do i = 1, dim2
        sendbuf(1:6*dn,i,j) = array_in(iend_nrc-6*dn+1:iend_nrc,i,j)
     end do
  end do
  call MPI_Sendrecv(sendbuf, srcount, MPI_DOUBLE_PRECISION, pp, 1, &
                    recvbuf, srcount, MPI_DOUBLE_PRECISION, pn, 1, &
                    mpi_natm_world, stat, ierr)
  do j = 1, dim3
     do i = 1, dim2
        array_out(ista_nrc-6*dn:ista_nrc-1,i,j) = recvbuf(1:6*dn,i,j)
     end do
  end do

  do j = 1, dim3
     do i = 1, dim2
        sendbuf(1:6*dn,i,j) = array_in(ista_nrc:ista_nrc+6*dn-1,i,j)
     end do
  end do
  call MPI_Sendrecv(sendbuf, srcount, MPI_DOUBLE_PRECISION, pn, 1, &
                    recvbuf, srcount, MPI_DOUBLE_PRECISION, pp, 1, &
                    mpi_natm_world, stat, ierr)
  do j = 1, dim3
     do i = 1, dim2
        array_out(iend_nrc+1:iend_nrc+6*dn,i,j) = recvbuf(1:6*dn,i,j)
     end do
  end do

  deallocate(sendbuf)
  deallocate(recvbuf)

  return
  end subroutine boundary_exchange_dim3

  subroutine decomp_vxc_ae_k_r_3D(vxc_ae_k, vxc_ae_k_3D, mmesh, msph, natm)
  implicit none
!!  include 'mpif.h'
  real(kind=DP), intent(out) :: vxc_ae_k   (mmesh,            nspin,msph,natm)
  real(kind=DP), intent(in)  :: vxc_ae_k_3D(ista_nrc:iend_nrc,nspin,msph,ne_natm)
  integer, intent(in) :: mmesh, msph, natm
  integer :: sendcount, recvcounts(0:npes-1), displs(0:npes-1)
  integer(kind=4) :: ierr, p, rank_natm, rank_nrc, ia, isp, is, ir, ind
  real(kind=DP), allocatable, dimension(:) :: tmp(:)

  allocate(tmp(mmesh*nspin*msph*natm))
  sendcount = ne_nrc*nspin*msph*ne_natm

  do p = 0, npes - 1
     rank_natm = mod(p,nrank_natm)
     rank_nrc  = p/nrank_natm
     recvcounts(p) = nel_nrc(rank_nrc)*nspin*msph*nel_natm(rank_natm)
     if(p == 0) then
        displs(p) = 0
     else
        displs(p) = displs(p-1) + recvcounts(p-1)
     end if
  end do
  call mpi_allgatherv(vxc_ae_k_3D,sendcount,        MPI_DOUBLE_PRECISION, &
                      tmp,        recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_CommGroup,ierr)
  do p = 0, npes - 1
     rank_natm = mod(p,nrank_natm)
     rank_nrc  = p/nrank_natm
     ind = displs(p) + 1
     if(is_nrc(rank_nrc) /= 0) then
        do ia = is_natm(rank_natm), ie_natm(rank_natm)
           do isp = 1, msph
              do is = 1, nspin
                 do ir = is_nrc(rank_nrc), ie_nrc(rank_nrc)
                    vxc_ae_k(ir,is,isp,ia) = tmp(ind)
                    ind = ind + 1
                 end do
              end do
           end do
        end do
     end if
  end do

  deallocate(tmp)

  end subroutine decomp_vxc_ae_k_r_3D
! ==============================================================================

end module m_PAW_XC_Potential
