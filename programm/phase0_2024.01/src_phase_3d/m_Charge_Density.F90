!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
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
#ifdef __TIMER_IODO__
#   define __TIMER_IODO_START(a)   call timer_sta(a)
#   define __TIMER_IODO_STOP(a)    call timer_end(a)
#else
#   define __TIMER_IODO_START(a)
#   define __TIMER_IODO_STOP(a)
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
#ifdef __TIMER_IOCOMM__
#   define __TIMER_IOCOMM_START_w_BARRIER(str,a)   call timer_barrier(str) ;   call timer_sta(a)
#   define __TIMER_IOCOMM_START(a)       call timer_sta(a)
#   define __TIMER_IOCOMM_STOP(a)        call timer_end(a)
#else
#   define __TIMER_IOCOMM_START_w_BARRIER(str,a)
#   define __TIMER_IOCOMM_START(a)
#   define __TIMER_IOCOMM_STOP(a)
#endif

#ifdef VPP
#define _VECTOR_TUNING_
#endif
#ifdef SX
#define _VECTOR_TUNING_
#endif
#ifdef HIUX
#define _VECTOR_TUNING_
#endif

#define MOMENT_AS_PSEUDO_VECTOR

module m_Charge_Density
! $Id: m_Charge_Density.F90 633 2020-12-01 05:11:03Z jkoga $

#ifdef MPI_FFTW
  use, intrinsic :: iso_c_binding
#endif

  use m_Const_Parameters,    only : BUCS, DP, PAI2, DIRECT,OFF,zi,SKIP &
       &                          , EXECUT,SIMPLE_CUBIC,NO,ANTIFERRO, ON &
       &                          , OLD, NEXT, PAI, VTK, CUBE, DENSITY_ONLY, Gauss_distrib_func &
       &                          , VERY_NARROW, from_PSEUDOPOTENTIAL_FILE &
!!$       &                          , unit_conv_byname &
       &                          , GAMMA, DELTA, DELTA10, ELECTRON, INVERSE, YES &
       &                          , SPLINE_INTERPOLATION
  use m_IterationNumbers,    only : iteration,iteration_electronic,iteration_unit_cell
  use m_Parallelization,     only : MPI_CommGroup,ista_e,iend_e,istep_e,map_z,np_e &
#ifdef NEC_TUNE_SOFT
       &                          , itask,ista_e_smp,iend_e_smp &
#endif
       &                          , map_k,map_ek,myrank_k &
       &                          , ista_kngp,iend_kngp,is_kngp,ie_kngp,np_kngp,mp_kngp &
       &                          , ista_kngp_prev, iend_kngp_prev, np_kngp_prev &
       &                          , np_kngp_gw     &
       &                          , npes,mype,ierr &
!!$       &                          , is_kgpm,ie_kgpm,ista_kgpm,iend_kgpm,mp_kgpm &
       &                          , is_atm, ie_atm, nel_atm, np_atm, mp_atm &
       &                          , ista_atm, iend_atm &
       &                          , ista_kngp_gw, iend_kngp_gw &
       &                          , ista_spin, iend_spin, nrank_s, mpi_keg_world &
       &                          , mpi_skg_world, mpi_sge_world
  use m_Timing,              only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters,  only : nspin,ipri,ipriwf,iprichargedensity,cdel_critical, nel_Ylm &
       &                          , istress,sw_fine_STM_simulation,kimg,af,neg &
       &                          , charge_filetype, charge_title, initial_chg &
       &                          , sw_add_corecharge_rspace, eval_corecharge_on_Gspace &
       &                          , iprichargemixing, ipritotalcharge &
       &                          , initial_charge_filetype &
       &                          , m_CtrlP_cachesize &
       &                          , sw_fft_xzy &
       &                          , sw_subset_only, minxyz, maxxyz &
       &                          , sw_extrapolate_charge, rms_threshold, sw_communicator_for_chg &
       &                          , sw_interpolate_charge, interpolation_method_chg, read_charge_hardpart
  use m_Crystal_Structure,   only : nopr, tau, univol, rltv, nbztyp, altv, total_spin, sw_fix_total_spin, additional_charge &
       &                          , univol_prev
  use m_Ionic_System,        only : ntyp, natm, natm2, ityp, iwei, iatomn, iatom, alfa &
       &                          , zeta1, zfm3_l, pos, cps, qex &
       &                          , nopr_supercell, iop_supercell &
#ifdef CHARGE_AVERAGE_WITH_SUPERCELLOPERATIONS
       &                          , tau_supercell &
#endif
       &                          , m_IS_pack_all_ions_in_uc, speciesname
  use m_PseudoPotential,     only : ival,nlmt,ilmt,ltp,mtp,qitg_l,dl2p&
       &                          , taup,il2p,isph,iqitg,lmta,modnrm &
       &                          , nlmta,nqitg, flg_paw,ipaw &
       &                          , nylm_paw,iylm_paw,crotylm_paw &
       &                          , ia2ia_symmtry_op_inv &
       &                          , rhvg_l, rhpcg_l, itpcc, rhcg_l &
       &                          , m_PP_find_maximum_l &
       &                          , m_PP_include_vanderbilt_pot &
       &                          , m_PP_set_index_arrays1 &
       &                          , m_PP_set_index_arrays2 &
       &                          , nmesh, rmax, xh, rhcorpw, radr_paw, rhcg_l &
       &                          , m_PP_rd_PAW_parameters
  use m_PlaneWaveBasisSet,   only : kg,kgp,kgp_reduced, ngpt_l,ngabc,gr_l,igf,kgpm,ylm_l &
       &                          , m_pwBS_sphrp2, igfp_l, igfp_nonpara, igfp_l_prev, kgp_reduced_prev
#ifdef _MPIFFT_
  use m_PlaneWaveBasisSet,   only : igfp_l_c
#endif
  use m_PlaneWaveBasisSet,   only : ngabc_kngp_l, ngabc_kngp_B_l,fp_l
  use m_FFT,                 only : m_FFT_alloc_WF_work &
       &                          , m_FFT_dealloc_WF_work &
       &                          , m_FFT_alloc_CD_box, m_FFT_dealloc_CD_box &
!fj$$       &                          , m_FFT_WF,fft_box_size_WF, nfft,nfftp,nfftps &
       &                          , fft_box_size_WF, nfft,nfftp,nfftps &
       &                          , fft_box_size_CD_nonpara, nfftp_nonpara &
       &                          , m_FFT_CD0, m_FFT_CD_inverse0, fft_box_size_CD &
       &                          , fft_box_size_CD_prev, fft_box_size_CD_nonpara_prev
#ifdef MPI_FFTW
  use m_FFT,                 only : afft_mpifftw, bfft_mpifftw, m_FFT_Direct_MPI_FFTW &
       &                          , afft_mpifftw_kimg1, bfft_mpifftw_kimg1
  use m_Electronic_Structure,only : m_ES_WF_in_Rspace_mpifftw, m_ES_fftbox_map
#endif
  use m_Kpoints,             only : k_symmetry
  use m_Electronic_Structure,only : totch,neordr,occup_l,fsr_l,fsi_l
  use m_Control_Parameters,  only : nblocksize_fftw, nblocksize_fftw_is_given     &
 &                                , nblocksize_gather_f_is_given                  &
 &                                , nblocksize_gather_f
#ifdef MPI_FFTW
  use m_Control_Parameters,  only : sw_mpi_fftw
#endif
  use m_Electronic_Structure,only : nrvf_ordr &
 &                                , m_ES_WF_in_Rspace_3D                          &
 &                                , m_ES_gather_f_3d_to_2d_blk                    &
 &                                , m_ES_alloc_fsr_l_2d, m_ES_dealloc_fsr_l_2d    &
 &                                , m_ES_alloc_fsi_l_2d, m_ES_dealloc_fsi_l_2d    &
 &                                , fsr_l_2D, fsi_l_2D, zaj_l &
 &                                , m_ES_dealloc_chgsoft, chg_softpart, chg_has_been_calculated
  use m_Parallelization,     only : nel_fft_z , nel_fft_y, nel_fft_x &
       &                          , np_fft_x, np_fft_y, np_fft_z  &
       &                          , fft_X_x_nel, fft_X_y_nel, fft_X_z_nel  &
       &                          , xyz_fft_x, mp_fft_x &
       &                          , nrank_e, myrank_g, nrank_g, ista_k, iend_k   &
       &                          , nrank_chg, myrank_chg   &
       &                          , mpi_ke_world, mpi_chg_world, mpi_kg_world, mpi_k_world, mpi_ge_world, mpi_ske_world
  use m_FFT,                 only : m_FFT_Direct_3D, m_FFT_Direct_XYZ_3D
#ifdef FFT_3D_DIVISION
  use m_FFT,                 only : m_FFT_Direct_3DIV_3D
#endif
  use m_PlaneWaveBasisSet,   only : m_pwBS_sphrp2_3D

! ===================== added by  K. Tagami ================== 5.0
  use m_Control_Parameters,   only : sw_mix_charge_hardpart
  use m_Control_Parameters,   only : sw_update_charge_hsr
! =========================================================== 5.0


! ====================================== added by K. Tagami ================== 11.0
  use m_Control_Parameters,  only : noncol, ndim_spinor, ndim_magmom, &
       &                            import_collinear_spindensity
  use m_Ionic_System,         only : mag_direction0_atomtyp, mag_moment0_atomtyp, &
       &                             ionic_charge_atomtyp, &
       &                             mag_moment0_atomtyp_is_defined, &
       &                             ionic_charge_atoms, mag_moment0_atoms, &
       &                             mag_moment0_atoms_is_defined
  use m_Crystal_Structure,      only :  op, sw_fix_global_quantz_axis, &
       &                                Global_Quantz_Axis_Fixed
! ============================================================================ 11.0
  use m_Parallelization, only : ista_fs, iend_fs, myrank_e, map_e
  use m_ErrorMessages,        only : EOF_REACHED

  use m_PlaneWaveBasisSet,  only : kgp_prev
! ===== ASMS ==== 2018/02/26
  use m_Control_Parameters,   only : charge_symm_mode
  use m_Kpoints,             only : nopr_from_fbz_to_ibz, flg_opr_from_fbz_to_ibz, &
       &                            num_star_of_k, star_of_k,  iopr_k_fbz_to_ibz
  use m_Const_Parameters,    only : chg_symm_level1, chg_symm_level2
! ===== ASMS ==== 2018/02/26
  use mpi


  implicit none

  real(kind=DP),public,target,allocatable,dimension(:,:,:) :: chgq_l, chgqo_l  ! d(ista_kngp:iend_kngp,kimg,nspin)
  real(kind=DP),public,target,allocatable,dimension(:,:,:) :: chgq_l_prev
  real(kind=DP),private, allocatable, dimension(:,:,:) ::  chgq_tmp ! d(ista_kngp:iend_kngp,kimg,nspin)
  real(kind=DP),public,allocatable,dimension(:,:,:)     :: chgsoft
  real(kind=DP),public,allocatable,target, dimension(:,:,:,:)   :: hsr ! d(natm,nlmt,nlmt,nspin)
  real(kind=DP),public,allocatable,target, dimension(:,:,:,:)   :: hsro ! d(natm,nlmt,nlmt,nspin)
  real(kind=DP),public,allocatable,target, dimension(:,:,:,:)   :: hsi, hsio ! == added by K. Tagami == 11.0
! ==== KT_add === 2014/08/25
  real(kind=DP),public,target, allocatable :: hsr_add(:,:,:,:) ! d(natm,nlmt,nlmt,nspin)
  real(kind=DP),public,target, allocatable :: hsi_add(:,:,:,:) ! d(natm,nlmt,nlmt,nspin)
! =============== 2014/08/25

  real(kind=DP),public ::                                  total_charge

  real(kind=DP),private,allocatable, dimension(:,:,:) ::       chgq_mpi1, chgq_mpi2  ! MPI

  real(kind=DP),private                                 :: cdelt

  real(kind=DP),private,allocatable, dimension(:)       :: afft, bfft
  real(kind=DP),private,allocatable, dimension(:)       :: afft_mpi1

  real(kind=DP), public, pointer, dimension(:,:)        :: work
  real(kind=DP),private, allocatable,target,dimension(:):: zfcos, zfsin
  integer,private, pointer, dimension(:)                :: il3

  ! -- For Broyden and DFP mixing method --
!!$  real(DP),private,allocatable,dimension(:,:,:)     :: f      !d(nbxmix,nbxmix,nspin)
!!$  real(DP),private,allocatable,dimension(:)         :: g      !d(nbxmix)

  real(kind=DP),allocatable,dimension(:,:) ::          chgq_enl     ! d(kgp,kimg)
  real(kind=DP),allocatable,dimension(:,:) ::          qitg_enl     ! d(kgp,nqitg)
  real(kind=DP),allocatable,target,dimension(:,:) ::   ylm_enl      ! d(kgp,:)
  integer, allocatable,dimension(:) ::                 igfp_enl     ! d(kgp)
  integer, allocatable,dimension(:,:) ::               ngpt_enl     ! d(kgp,nopr+af)

  real(kind=DP), allocatable, dimension(:,:,:,:) :: extpl_target_history

  real(kind=DP), allocatable, dimension(:,:,:)   :: chgq_l_ek
  real(kind=DP), allocatable, dimension(:,:,:,:) :: chgq_l_pc
  integer                                        :: nEwindow_ek

#ifdef MPI_FFTW
  include 'fftw3-mpi.f03'
#endif

!  include 'mpif.h'
  integer istatus(mpi_status_size)

! --- contained subroutines ---
!   1. m_CD_set_ylm_enl_etc      <-(m_Ldos)
!   2. m_CD_dealloc_ylm_enl_etc  <-(m_Ldos)
!   3. m_CD_alloc_chgq           <-(Preparation_for_mpi)
!   4. m_CD_alloc_hsr            <-(Preparation_for_mpi)
!   5. m_CD_dealloc_chgq         <-(Finalization_of_mpi)
!   6. m_CD_check                <-(ChargeDensity_Mixing)
!   8. m_CD_is_converged_core    <-(Convergence_Check)
!   9. m_CD_convergence_check     <-(ChargeDensity_Construction)
!  11. m_CD_hardpart             <-(ChargeDensity_Construction)
!  12. summation_of_ff           <-(@11)
!  13. add_hardpart_to_chgq_l    <-(@11)
!     - add_hardpart_to_chgq_l_core
!  14. m_CD_softpart             <-(ChargeDensity_Construction)
!     - substitute_CD_for_chgq  - add_occupied_densities
!  15. m_CD_cp_chgq_to_chgqo          <-(Renewal_of_WaveFunctions)
!  16. m_CD_wd_chgq_l_small_portion   <-(@11),(@13),(@14),(Initial_Electronic_Structure)
!  17. m_CD_initial_CD_by_Gauss_func  <-(Initial_Electronic_Structure)
!  18. m_CD_rd_chgq              <-(Initial_Electronic_Structure)
!  19. m_CD_wd_chgq              <-(WriteDownData_onto_Files)
!  20. charge_average            <-(@10), (@13), (@53), (@54), (@55), (@60)
!  21. cp_chg_to_work            <-(@20)
!  61. m_CD_alloc_rspace_charge       <-(Postprocessing)
!  62. m_CD_dealloc_rspace_charge     <-(Postprocessing)
!  63. m_CD_rspace_charge             <-(Postprocessing)
!    - wdchgr - map_valence_charge_to_fft_box
!  64. m_CD_hardpart_sub              <-(m_Ldos)
!  65. summation_of_ff_sub            <-(@64)
!  66. m_CD_map_chgq_to_fft_box       <-(m_Ldos)
!  67. m_CD_map_chgqenl_to_fft_box    <-(m_Ldos)
!  68. m_CD_hardpart_sub2             <-(m_Ldos)
!  69. summation_of_ff_sub_enl        <-(@68)
!  70. add_hardpart_to_chgq_enl       <-(@68)
!    - add_hardpart_to_chgq_e_core
! *: Primary subroutines

contains

  subroutine m_CD_set_ylm_enl_etc()
    integer :: n, ilm3 ,i, j
    real(kind=DP),allocatable,dimension(:) :: ylm_t, ylm_mpi, qitg_t, qitg_mpi
    integer, allocatable, dimension(:) ::     igfp_t, igfp_mpi, ngpt_t, ngpt_mpi

    call m_PP_find_maximum_l(n)
    n = (n-1) + (n-1) + 1
    allocate(ylm_enl(kgp,n**2))
    allocate(igfp_enl(kgp))
    allocate(chgq_enl(kgp,kimg))
    allocate(qitg_enl(kgp,nqitg)); qitg_enl = 0.d0

! ============================== Added by K. Tagami ==========
        ylm_enl = 0.0d0; igfp_enl = 0.0d0;
        chgq_enl = 0.0d0
! ==========================================================
    if(ipri >= 2) then
       write(6,'(" !** n*n = ",i7," nel_Ylm = ",i7," <<m_CD_set_ylm_enl_etc>>")') n**2, nel_Ylm
       write(6,'(" !** nqitg = ",i7," <<m_CD_set_ylm_enl_etc>>")') nqitg
    end if

    ! -- ylm_enl --
    ylm_enl = 0.d0
    if(npes > 1) then
       allocate(ylm_t(kgp))
       allocate(ylm_mpi(kgp))
! ================================= Added by K. Tagami ========
        ylm_t = 0.0d0; ylm_mpi = 0.0d0
! =============================================================
    else
       if(nel_Ylm < n**2) then
          allocate(ylm_t(kgp))
! ================================ Added by K. Tagami ========
          ylm_t = 0.0d0
! ============================================================
       end if
    end if
    do ilm3 = 1, nel_Ylm
       if(npes > 1) then
          ylm_t = 0.d0
          do i = ista_kngp, iend_kngp
             ylm_t(i) = ylm_l(i,ilm3)
          end do
          call mpi_allreduce(ylm_t, ylm_mpi, kgp,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
          ylm_enl(:,ilm3) = ylm_mpi(:)
       else
          ylm_enl(:,ilm3) = ylm_l(:,ilm3)
       end if
    end do

    do ilm3 = nel_Ylm+1, n**2
!!$       call m_pwBS_sphrp2(ilm3,rltv,1,kgp,ylm_t)
       ylm_t = 0.d0
       call m_pwBS_sphrp2(ilm3,rltv,ista_kngp,iend_kngp,ylm_t)
       if(npes > 1 ) then
!!$          call mpi_allreduce(ylm_t,ylm_mpi,kgp,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
          call mpi_allreduce(ylm_t,ylm_mpi,kgp,mpi_double_precision,mpi_sum,mpi_chg_world,ierr)
          ylm_enl(:,ilm3) = ylm_mpi(:)
       else
          ylm_enl(:,ilm3) = ylm_t
       end if
    end do
    if(npes > 1) then
       deallocate(ylm_t)
       deallocate(ylm_mpi)
    else
       if(nel_Ylm < n**2) then
          deallocate(ylm_t)
       end if
    end if

    ! -- igfp_enl --
    if(npes > 1) then
       allocate(igfp_t(kgp))
       allocate(igfp_mpi(kgp))
! ==============================Added by K. Tagami========
	igfp_mpi = 0
! ======================================================
       igfp_t = 0
       do i = ista_kngp, iend_kngp
#ifdef _MPIFFT_
          igfp_t(i) = igfp_nonpara(i)
#else
          igfp_t(i) = igfp_l(i)
#endif
       end do
!!$       call mpi_allreduce(igfp_t,igfp_mpi,kgp,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
       call mpi_allreduce(igfp_t,igfp_mpi,kgp,mpi_integer,mpi_sum,mpi_chg_world,ierr)
       igfp_enl(:) = igfp_mpi(:)
       deallocate(igfp_t)
       deallocate(igfp_mpi)
    else
       do i = ista_kngp, iend_kngp
#ifdef _MPIFFT_
          igfp_enl(i) = igfp_nonpara(i)
#else
          igfp_enl(i) = igfp_l(i)
#endif
       end do
    end if

    ! -- qitg_enl --
    if(npes > 1) then
       allocate(qitg_t(kgp))
       allocate(qitg_mpi(kgp))
       qitg_enl = 0.d0
! ============================ Added by K. Tagami ==================
       qitg_mpi = 0.0d0
! ==================================================================
       do j = 1, nqitg
          qitg_t = 0.d0
          do i = ista_kngp, iend_kngp
             qitg_t(i) = qitg_l(i,j)
          end do
!!$          call mpi_allreduce(qitg_t,qitg_mpi,kgp,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
          call mpi_allreduce(qitg_t,qitg_mpi,kgp,mpi_double_precision,mpi_sum,mpi_chg_world,ierr)
          qitg_enl(:,j) = qitg_mpi(:)
       end do
       deallocate(qitg_t)
       deallocate(qitg_mpi)
    else
       qitg_enl(:,:) = qitg_l(:,:)
    end if

    ! -- ngpt_enl --
    allocate(ngpt_enl(kgp,nopr+af))
! ======================================== Added by K. Tagami =====
    ngpt_enl = 0
! =================================================================
    if(npes > 1) then
       allocate(ngpt_t(kgp))
       allocate(ngpt_mpi(kgp))
       do j = 1, nopr+af
          ngpt_t = 0
! ================================= Aded by K. Tagami =======
          ngpt_mpi = 0
! =============================================================
          do i = ista_kngp, iend_kngp
            ngpt_t(i) = ngpt_l(i,j)
          end do
!!$          call mpi_allreduce(ngpt_t,ngpt_mpi,kgp,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
          call mpi_allreduce(ngpt_t,ngpt_mpi,kgp,mpi_integer,mpi_sum,mpi_chg_world,ierr)
          ngpt_enl(:,j) = ngpt_mpi(:)
       end do
       deallocate(ngpt_t, ngpt_mpi)
     else
       ngpt_enl(:,:) = ngpt_l(:,:)
     end if

  end subroutine m_CD_set_ylm_enl_etc

  subroutine m_CD_dealloc_ylm_enl_etc()
    if(allocated(ylm_enl)) then
       deallocate(ylm_enl)
    else
       if(ipri >= 2) write(6,'(" ylm_enl is not allocated <<m_CD_dealloc_ylm_enl_etc>>")')
    end if
    if(allocated(igfp_enl)) then
       deallocate(igfp_enl)
    else
       if(ipri >= 2) write(6,'(" igfp_enl is not allocated <<m_CD_dealloc_ylm_enl_etc>>")')
    end if

    if(allocated(chgq_enl)) then
       deallocate(chgq_enl)
    else
       if(ipri >= 2) write(6,'(" chgq_enl is not allocated <<m_CD_dealloc_ylm_enl_etc>>")')
    end if
    if(allocated(qitg_enl)) then
       deallocate(qitg_enl)
    else
       if(ipri >= 2) write(6,'(" qitg_enl is not allocated <<m_CD_dealloc_ylm_enl_etc>>")')
    end if

    if(allocated(ngpt_enl)) then
       deallocate(ngpt_enl)
    else
       if(ipri >= 2) write(6,'(" ngpt_enl is not allocated <<m_CD_dealloc_ylm_enl_etc>>")')
    end if

  end subroutine m_CD_dealloc_ylm_enl_etc


  subroutine m_CD_alloc_chgq
    allocate(chgq_l(ista_kngp:iend_kngp,kimg,nspin)); chgq_l = 0.d0
    allocate(chgqo_l(ista_kngp:iend_kngp,kimg,nspin)); chgqo_l = 0.d0
    if(istress == ON .or. sw_fine_STM_simulation == ON) then
       allocate(chgsoft(ista_kngp:iend_kngp,kimg,nspin))
       chgsoft = 0.d0
    end if
    if(sw_extrapolate_charge==ON)then
       if(.not.allocated(extpl_target_history)) allocate(extpl_target_history(ista_kngp:iend_kngp,kimg,nspin,3))
       extpl_target_history = 0.d0
    endif
  end subroutine m_CD_alloc_chgq

  subroutine m_CD_alloc_chgsoft
    if(.not.allocated(chgsoft)) then
       allocate(chgsoft(ista_kngp:iend_kngp,kimg,nspin))
       chgsoft = 0.d0
    end if
  end subroutine m_CD_alloc_chgsoft


  subroutine m_CD_alloc_hsr()
!!$    print '(" natm, nlmt, nspin = ",3i5)',natm,nlmt,nspin
    allocate(hsr(natm,nlmt,nlmt,nspin)); hsr = 0.d0

! ======================== modified by K. Tagami ======================= 5.0
!    if(flg_paw) then
!       allocate(hsro(natm,nlmt,nlmt,nspin)); hsro = 0.d0
!    endif
!
!    if ( flg_paw .or. sw_mix_charge_hardpart == ON ) then
       allocate(hsro(natm,nlmt,nlmt,nspin)); hsro = 0.d0
!    endif
  end subroutine m_CD_alloc_hsr

  subroutine m_CD_dealloc_chgq(store_prev)
    logical, intent(in), optional :: store_prev
    if(present(store_prev)) then
      if(store_prev) then
        if(allocated(chgq_l_prev)) deallocate(chgq_l_prev)
        allocate(chgq_l_prev(ista_kngp:iend_kngp,kimg,ndim_magmom))
        chgq_l_prev = chgq_l
        univol_prev = univol
      endif
    endif
    if(allocated(chgsoft)) deallocate(chgsoft)
    if(allocated(chgqo_l)) deallocate(chgqo_l)
    if(allocated(chgq_l)) deallocate(chgq_l)
    if( allocated(extpl_target_history) ) deallocate(extpl_target_history)
  end subroutine m_CD_dealloc_chgq

  subroutine m_CD_check(nfout)
    integer, intent(in) :: nfout
    integer :: i, ispin, ip
    integer :: ipritotalcharge_0
    real(kind=DP) :: totch_old, totch_new
    real(kind=DP),allocatable,dimension(:,:) :: chg_t ! d(2,OLD:NEXT)

    if(myrank_chg == 0) ipritotalcharge_0 = ipritotalcharge
    if(nrank_chg > 1) call mpi_bcast(ipritotalcharge_0,1,mpi_integer,0,mpi_chg_world,ierr)

    if(.not.(ipritotalcharge_0 >= 2 .or. (ipritotalcharge_0 >= 1 .and. (nspin == 2 .and. af == 0)))) return

                                                  __TIMER_SUB_START(1145)
    if(nspin == 2 .and. af == 0) then
        allocate(chg_t(nspin,OLD:NEXT))

! ====================================== Added by K. Tagami =======
! === DEBUG by tkato 2013/08/28 ================================================
!    chg_t = 0.0d0
    if(nspin == 2 .and. af == 0) chg_t = 0.0d0
! ==============================================================================
! ==============================================================
    endif

    totch_old = 0.d0
    totch_new = 0.d0
    ip = 0
                                                  __TIMER_DO_START(1193)
    if(nrank_chg > 1) then
       do i = 0, nrank_chg-1
          if( is_kngp(i) <= 1 .and. 1 <= ie_kngp(i)) then
             ip = i
             exit
          end if
       end do
    end if
                                                  __TIMER_DO_STOP(1193)
    i = 1
    if(ista_kngp <= i .and. i <= ista_kngp) then
                                                  __TIMER_DO_START(1194)
       do ispin = 1, nspin, af+1
          totch_old = totch_old + chgqo_l(i,1,ispin)
          totch_new = totch_new + chgq_l(i,1,ispin)
       end do
                                                  __TIMER_DO_STOP(1194)
       if(nspin == 2 .and. af == 0) then
                                                  __TIMER_DO_START(1195)
          do ispin = 1, nspin
             chg_t(ispin,NEXT) = chgq_l(i,1,ispin)
             chg_t(ispin,OLD) = chgqo_l(i,1,ispin)
          end do
                                                  __TIMER_DO_STOP(1195)
       end if

       if(af == 1) then
          totch_old = totch_old*2.d0
          totch_new = totch_new*2.d0
       end if
    end if
    if(nrank_chg > 1) then
       call mpi_barrier(mpi_chg_world,ierr)
                                                 __TIMER_COMM_START_w_BARRIER(mpi_chg_world,1196)
       call mpi_bcast(totch_old,1,mpi_double_precision,ip,mpi_chg_world,ierr)
       call mpi_bcast(totch_new,1,mpi_double_precision,ip,mpi_chg_world,ierr)
                                                 __TIMER_COMM_STOP(1196)
    end if
    if(nspin == 2 .and. af == 0) then
       if(nrank_chg > 1) then
                                                 __TIMER_COMM_START_w_BARRIER(mpi_chg_world,1197)
          call mpi_barrier(mpi_chg_world,ierr)
          call mpi_bcast(chg_t,4,mpi_double_precision,ip,mpi_chg_world,ierr)
                                                 __TIMER_COMM_STOP(1197)
       end if
       if((ipritotalcharge==1 .and. iteration_electronic==1).or.ipritotalcharge>=2) &
            & write(nfout,96) "OLD",chg_t(1,OLD)*univol,chg_t(2,OLD)*univol,totch_old*univol
       if(ipritotalcharge>=1) &
            & write(nfout,96) "NEW",chg_t(1,NEXT)*univol,chg_t(2,NEXT)*univol,totch_new*univol
    else
       if((ipritotalcharge >= 2 .and. iteration_electronic==1).or.ipritotalcharge>=3) &
            &                 write(nfout,98) "OLD",totch_old*univol
       if(ipritotalcharge>=2) write(nfout,98) "NEW",totch_new*univol
    end if

96  format(' !',a3,' total charge (UP, DOWN, SUM) = ' &
         & ,f14.8,' (+)',f14.8,' (=)',f14.8)
98  format(' !',a3,' total charge = ', f14.8)

    if(nspin == 2 .and. af == 0) deallocate(chg_t)
                                                 __TIMER_SUB_STOP(1145)
  end subroutine m_CD_check


  logical function m_CD_is_converged_core(nfout)
    integer, intent(in) :: nfout
    m_CD_is_converged_core = .false.
!!$       if(cdelt < cdel_critical) then
    if(ipri >= 2) write(nfout,'(" cdelt = ",d20.8, " , iter = ",i7)') cdelt, iteration
          m_CD_is_converged_core = .true.
!!$       end if
  end function m_CD_is_converged_core

  subroutine m_CD_convergence_check(nfout)
    integer, intent(in) :: nfout
    real(kind=DP) :: total_charge_t
    real(kind=DP), parameter :: totch_critical = 1.d-15
    integer :: is, ik, i

! =================================== added by K. Tagami =============== 11.0
    integer :: ni
! ====================================================================== 11.0

    real(kind=DP) :: cdelt_mpi
    integer :: idp, nlp, nmp, nnp, nlphf, ip
    real(kind=DP) :: s1, s2
    integer :: id_sname = -1
    call tstatc0_begin('m_CD_convergence_check ',id_sname,1)
                                                 __TIMER_SUB_START(733)

! --> T. Yamasaki, 18 July 2008
!!$    call m_CD_alloc_rspace_charge()
!!$    cdelt = 0.d0
!!$    do is = 1, nspin, af+1
!!$       call map_valence_charge_to_fftbox(is)
!!$       call m_FFT_CD_inverse0(nfout,afft)
!!$
!!$       idp = fft_box_size_CD(1,0)
!!$       nlp = fft_box_size_CD(1,1)
!!$       nmp = fft_box_size_CD(2,1)
!!$       nnp = fft_box_size_CD(3,1)
!!$
!!$       if(kimg == 1) then
!!$          nlphf = idp/2
!!$       else
!!$          nlphf = idp
!!$       end if
!!$
!!$       s1 = 0.d0; s2 = 0.d0
!!$       do ip = 1, nlphf*nmp*nnp,2
!!$          s1 = s1 + afft(ip)
!!$          s2 = s2 + afft(ip+1)
!!$       end do
!!$       s1 = s1*univol/product(fft_box_size(1:3,1))
!!$       s2 = s1*univol/product(fft_box_size(1:3,1))
!!$       if(iprichargemixing >= 1) &
!!$            & write(nfout,'(" !cdelt realspace-summation          s1, s2 = ",2d20.8)') s1, s2
!!$
!!$       call map_valence_charge_dif_to_fftbox(is)
!!$       call m_FFT_CD_inverse0(nfout,afft)
!!$
!!$       s1 = 0.d0; s2 = 0.d0
!!$       do ip = 1, nlphf*nmp*nnp,2
!!$          s1 = s1 + dabs(afft(ip))
!!$          s2 = s2 + dabs(afft(ip+1))
!!$       end do
!!$       s1 = s1*univol/product(fft_box_size(1:3,1))
!!$       s2 = s2*univol/product(fft_box_size(1:3,1))
!!$       if(iprichargemixing >= 1) &
!!$            &  write(nfout,'(" !cdelt realspace-summation of diff, s1, s2 = ",2d20.8)') s1, s2
!!$       cdelt = cdelt + s1/totch
!!$    end do
!!$    if(iprichargemixing >= 1) &
!!$         & write(nfout,'(" !cdelt delta_charge(",i7,") = ",d20.8)') iteration, cdelt
!!$
    cdelt_mpi = 0.d0
                                                  __TIMER_DO_START(857)
! ==================================== modified by K. Tagami ================ 11.0
!    do is = 1, nspin, af+1
!       do ik = 1, kimg
!          do i = ista_kngp, iend_kngp  !for mpi
!             cdelt_mpi = cdelt_mpi + (chgq_l(i,ik,is)-chgqo_l(i,ik,is))**2
!          end do
!       end do
!    end do
!
    if ( noncol ) then
      do is = 1, ndim_magmom
         do ik = 1, kimg
            do i = ista_kngp, iend_kngp  !for mpi
               cdelt_mpi = cdelt_mpi + (chgq_l(i,ik,is)-chgqo_l(i,ik,is))**2
            end do
         end do
      end do
    else
      do is = 1, nspin, af+1
         do ik = 1, kimg
            do i = ista_kngp, iend_kngp  !for mpi
               cdelt_mpi = cdelt_mpi + (chgq_l(i,ik,is)-chgqo_l(i,ik,is))**2
            end do
         end do
      end do
    endif
! ===================================================================== 11.0
                                                  __TIMER_DO_STOP(857)
                                                 __TIMER_COMM_START_w_BARRIER(mpi_chg_world,858)
    call mpi_allreduce(cdelt_mpi,cdelt,1 &
                   &  ,mpi_double_precision,mpi_sum,mpi_chg_world,ierr)
                                                 __TIMER_COMM_STOP(858)

!!$    call m_CD_dealloc_rspace_charge()
!!$
    total_charge_t = 0.d0
                                                 __TIMER_DO_START(859)
    if(myrank_chg==0) then

! ===================================== modified by K. Tagami ======== 11.0
!       do is = 1, nspin, af+1
!          total_charge_t = total_charge_t + chgq_l(1,1,is)
!       end do
!
       if ( noncol ) then
         ni = 1
         total_charge_t = total_charge_t + chgq_l(1,1,ni)
       else
         do is = 1, nspin, af+1
            total_charge_t = total_charge_t + chgq_l(1,1,is)
         end do
       endif
! ==================================================================== 11.0

    endif
                                                 __TIMER_DO_STOP(859)
                                                 __TIMER_COMM_START_w_BARRIER(mpi_chg_world,860)
    call mpi_bcast(total_charge_t,1,mpi_double_precision,0,mpi_chg_world,ierr)
                                                 __TIMER_COMM_STOP(860)
    if(total_charge_t < totch_critical) then
       cdelt = 0.d0
    else
       cdelt = dsqrt(cdelt)/totch
    end if
    if(iprichargemixing >= 2) &
         & write(nfout,'(" !cdelt delta_charge(",i7,") = ",d20.8)') iteration, cdelt
    call tstatc0_end(id_sname)
!!$  contains
!!$    subroutine map_valence_charge_to_fftbox(iloop)
!!$      integer, intent(in) :: iloop
!!$      integer :: j, i, ip
!!$
!!$      afft_mpi1 = 0.d0
!!$      do j = 1, kimg
!!$         do i = ista_kngp, iend_kngp
!!$#ifdef _MPIFFT_
!!$            ip = (igfp_nonpara(i)-1)*kimg + j
!!$#else
!!$            ip = (igfp_l(i)-1)*kimg + j
!!$#endif
!!$            afft_mpi1(ip) = afft_mpi1(ip) + chgq_l(i,j,iloop)
!!$         end do
!!$      end do
!!$      if(npes >= 2) then
!!$         call mpi_allreduce(afft_mpi1,afft,nfftp_nonpara,mpi_double_precision &
!!$              &  ,mpi_sum,MPI_CommGroup,ierr)
!!$      else
!!$         afft = afft_mpi1
!!$      end if
!!$    end subroutine map_valence_charge_to_fftbox
!!$
!!$    subroutine map_valence_charge_dif_to_fftbox(iloop)
!!$      integer, intent(in) :: iloop
!!$      integer :: j, i, ip
!!$      afft_mpi1 = 0.d0
!!$      do j = 1, kimg
!!$         do i = ista_kngp, iend_kngp
!!$#ifdef _MPIFFT_
!!$            ip = (igfp_nonpara(i)-1)*kimg + j
!!$#else
!!$            ip = (igfp_l(i)-1)*kimg + j
!!$#endif
!!$            afft_mpi1(ip) = afft_mpi1(ip) + (chgq_l(i,j,iloop)-chgqo_l(i,j,iloop))
!!$         end do
!!$      end do
!!$      if(npes >= 2) then
!!$         call mpi_allreduce(afft_mpi1,afft,nfftp_nonpara,mpi_double_precision &
!!$              &  ,mpi_sum,MPI_CommGroup,ierr)
!!$      else
!!$         afft = afft_mpi1
!!$      end if
!!$    end subroutine map_valence_charge_dif_to_fftbox
!!$! <--
                                                  __TIMER_SUB_STOP(733)
  end subroutine m_CD_convergence_check

  subroutine m_CD_hardpart(nfout,kv3)
    integer, intent(in) :: nfout, kv3

    integer, parameter :: DEBUGPRINTLEVEL = 3
    integer :: id_sname = -1
                                                  __TIMER_SUB_START(720)
#ifdef __FAPP__
    call fapp_start('cd_hardpart',1,1)
#endif
    call tstatc0_begin('m_CD_hardpart ',id_sname,1)

! ====================== added by K. Tagami =================== 5.0
    if ( sw_update_charge_hsr == OFF ) goto 100
! ============================================================= 5.0

    call summation_of_ff_3D(kv3,1) ! -(m_C.D.) (vnlsum) fsr_l, fsi_l, occup_l --> hsr
    if(flg_paw .and. iprichargedensity >= DEBUGPRINTLEVEL) then
       if(flg_paw) write(nfout,'(" -- hsr before symmtrz_of_ff --")')
       call wd_hsr(nfout)
    end if
!!$    if(flg_paw) call symmtrz_of_ff
!    if ( flg_paw .or. sw_mix_charge_hardpart == ON ) then
!       call symmtrz_of_ff
!    endif
! ============================================================== 5.0
100 continue

    call add_hardpart_to_chgq_l_3D(nfout,nspin,hsr, NO) ! (lclchg) add hardpart and make average

    if ( sw_update_charge_hsr == OFF ) goto 200

    if ( flg_paw .or. sw_mix_charge_hardpart == ON ) then
       call symmtrz_of_ff
    endif

200 continue

    if(iprichargedensity >= DEBUGPRINTLEVEL) then
       write(nfout,'(" -- hardpart summed --")')
       call m_CD_wd_chgq_l_small_portion(nfout)
    end if

    if(iprichargedensity >= DEBUGPRINTLEVEL) then
       write(nfout,'(" -- hsr --")')
       if(flg_paw) write(nfout,'(" -- hsr after symmtrz_of_ff --")')
       call wd_hsr(nfout)
    end if

#ifdef __FAPP__
    call fapp_stop('cd_hardpart',1,1)
#endif
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(720)
  end subroutine m_CD_hardpart


  subroutine m_CD_hardpart_hsr(nfout,kv3)
    integer, intent(in) :: nfout, kv3
    integer, parameter :: DEBUGPRINTLEVEL = 3

    integer :: id_sname = -1
                                                  __TIMER_SUB_START(734)
    call tstatc0_begin('m_CD_hardpart_hsr ',id_sname,1)

    call summation_of_ff_3D(kv3,0) ! -(m_C.D.) (vnlsum) fsr_l, fsi_l, occup_l --> hsr

    if(iprichargedensity >= DEBUGPRINTLEVEL) then
       write(nfout,'(" -- hsr --")')
       call wd_hsr(nfout)
    end if
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(734)
  end subroutine m_CD_hardpart_hsr


#define _PARA_ATOM_
  subroutine summation_of_ff_3D(kv3,iflag)
    integer, intent(in) :: kv3,iflag

    real(kind=DP)   :: w_n, d_factor
    integer         :: ia, is, k, i, lmt1, lmt2, p, q, it
    integer         :: lnblck, lnb1, lnb2, lnsize
    real(kind=DP),pointer,dimension(:,:,:,:) :: hsr_mpi
                                                  __TIMER_SUB_START(721)
#if 1
    if(nblocksize_gather_f_is_given) then
       lnblck = nblocksize_gather_f
       if(np_e < lnblck) then
          lnblck = np_e
       end if
       if(1 > lnblck) then
          lnblck = np_e
       end if
    else
       lnblck = np_e
    end if

    call m_ES_alloc_fsr_l_2d(lnblck, nlmta)

    d_factor = 2.d0/kv3
    hsr = 0.d0

                                                  __TIMER_DO_START(831)
!    do is = 1, nspin, af+1
    do is = ista_spin, iend_spin, af+1
       do k = is, kv3+is-nspin, nspin
          if(map_k(k) /= myrank_k) cycle            ! MPI
          if(k_symmetry(k) == GAMMA) then

             do lnb1 = 1, np_e , lnblck
                lnb2=min( lnb1+lnblck-1,np_e )
                lnsize = lnb2-lnb1 + 1
                call m_ES_gather_f_3d_to_2d_blk(fsr_l, fsr_l_2D, k, lnblck, lnb1, lnsize)
#ifdef _PARA_ATOM_
                do ia = ista_atm, iend_atm
#else
                do ia = 1, natm
#endif
                   it = ityp(ia)
                   do i = lnb1, lnb2                            ! MPI
                      w_n = occup_l(i,k)*d_factor
                      do lmt1 = 1, ilmt(it)
                         p = lmta(lmt1,ia)
                         do lmt2 = lmt1, ilmt(it)
                            q = lmta(lmt2,ia)
                            hsr(ia,lmt1,lmt2,is) = hsr(ia,lmt1,lmt2,is) + w_n *   &
                           &                       (fsr_l_2D(i-lnb1+1,p)*fsr_l_2D(i-lnb1+1,q))
                         end do! lmt2
                      end do! lmt1
                   end do! i
                end do! ia
             end do! lnb1
          else
             if (.not.allocated(fsi_l_2D)) call m_ES_alloc_fsi_l_2d(lnblck, nlmta)
             do lnb1 = 1, np_e , lnblck
                lnb2=min( lnb1+lnblck-1,np_e )
                lnsize = lnb2-lnb1 + 1
                call m_ES_gather_f_3d_to_2d_blk(fsr_l, fsr_l_2D, k, lnblck, lnb1, lnsize)
                call m_ES_gather_f_3d_to_2d_blk(fsi_l, fsi_l_2D, k, lnblck, lnb1, lnsize)
#ifdef _PARA_ATOM_
                do ia = ista_atm, iend_atm
#else
                do ia = 1, natm
#endif
                   it = ityp(ia)
                   do i = lnb1, lnb2                         ! MPI
                      w_n = occup_l(i,k)*d_factor
                      do lmt1 = 1, ilmt(it)
                         p = lmta(lmt1,ia)
                         do lmt2 = lmt1, ilmt(it)
                            q = lmta(lmt2,ia)
                            hsr(ia,lmt1,lmt2,is) = hsr(ia,lmt1,lmt2,is) + w_n *                 &
                           &                      (fsr_l_2D(i-lnb1+1,p)*fsr_l_2D(i-lnb1+1,q) +  &
                           &                       fsi_l_2D(i-lnb1+1,p)*fsi_l_2D(i-lnb1+1,q))
                         end do! lmt2
                      end do! lmt1
                   end do! i
                end do! ia
             end do! lnb1
          end if

       end do! ik
    end do! is
                                                  __TIMER_DO_STOP(831)

    call m_ES_dealloc_fsr_l_2d()
    call m_ES_dealloc_fsi_l_2d()

#else

    d_factor = 2.d0/kv3
    hsr = 0.d0
                                                  __TIMER_DO_START(831)
!    do is = 1, nspin, af+1
    do is = ista_spin, iend_spin, af+1
       do k = is, kv3+is-nspin, nspin
          if(map_k(k) /= myrank_k) cycle            ! MPI
          if(k_symmetry(k) == GAMMA) then
#ifdef _PARA_ATOM_
             do ia = ista_atm, iend_atm
#else
             do ia = 1, natm
#endif
                it = ityp(ia)
                do i = 1, np_e                            ! MPI
                   w_n = occup_l(i,k)*d_factor
                   do lmt1 = 1, ilmt(it)
                      p = lmta(lmt1,ia)
                         do lmt2 = lmt1, ilmt(it)
                            q = lmta(lmt2,ia)
                            hsr(ia,lmt1,lmt2,is) = hsr(ia,lmt1,lmt2,is) + w_n * &
                           &                       (fsr_gall(i,p,k)*fsr_gall(i,q,k))
                         end do! lmt2
                   end do! lmt1
                end do! i
             end do! ia
          else
#ifdef _PARA_ATOM_
             do ia = ista_atm, iend_atm
#else
             do ia = 1, natm
#endif
                it = ityp(ia)
                do i = 1, np_e                            ! MPI
                   w_n = occup_l(i,k)*d_factor
                   do lmt1 = 1, ilmt(it)
                      p = lmta(lmt1,ia)
                         do lmt2 = lmt1, ilmt(it)
                            q = lmta(lmt2,ia)
                            hsr(ia,lmt1,lmt2,is) = hsr(ia,lmt1,lmt2,is) + w_n * &
                           &   (fsr_gall(i,p,k)*fsr_gall(i,q,k) + fsi_gall(i,p,k)*fsi_gall(i,q,k))
                         end do! lmt2
                   end do! lmt1
                end do! i
             end do! ia
          end if

       end do! ik
    end do! is
                                                  __TIMER_DO_STOP(831)

#endif

#ifdef _PARA_ATOM_
    if(npes > 1) then
       allocate(hsr_mpi(natm,nlmt,nlmt,nspin))
! === DEBUG by tkato 2011/12/07 ================================================
!       hsr_mpi = 0
       hsr_mpi = 0.0d0
! ==============================================================================
                                                 __TIMER_COMM_START_w_BARRIER(MPI_CommGroup,832)
       call mpi_allreduce(hsr, hsr_mpi, natm*nlmt*nlmt*nspin, &
      &                   mpi_double_precision, mpi_sum, MPI_CommGroup, ierr)
                                                 __TIMER_COMM_STOP(832)
       hsr = hsr_mpi
       deallocate(hsr_mpi)
    end if
#else
    if(nrank_e > 1) then
       allocate(hsr_mpi(natm,nlmt,nlmt,nspin))
! === DEBUG by tkato 2011/12/07 ================================================
!	hsr_mpi = 0
       hsr_mpi = 0.0d0
! ==============================================================================
                                                 __TIMER_COMM_START_w_BARRIER(mpi_kg_world,833)
       call mpi_allreduce(hsr, hsr_mpi, natm*nlmt*nlmt*nspin, &
      &                   mpi_double_precision, mpi_sum, mpi_kg_world, ierr)
                                                 __TIMER_COMM_STOP(833)
       hsr = hsr_mpi
       deallocate(hsr_mpi)
    end if
#endif
                                                 __TIMER_SUB_STOP(721)
  end subroutine summation_of_ff_3D



  subroutine symmtrz_of_ff

    integer         :: ia, is, iopr, i, lmt1, lmt2, lmt3, lmt4, it
    real(kind=DP),pointer,dimension(:,:,:,:) :: hsr_mpi
    real(kind=DP),pointer,dimension(:,:,:,:) :: hsr_tmp
    integer :: il1,im1,it1,il2,im2,it2,il3,im3,it3,il4,im4,it4
    integer :: ii,jj,kk,ll,n,m,iii,jjj
    integer :: ja
    real(kind=DP) :: fi

!    allocate(hsr_mpi(natm,nlmt,nlmt,nspin))   ! MPI
    allocate(hsr_tmp(natm,nlmt,nlmt,nspin))

    hsr_tmp = hsr

    do ia=1,natm
        it=ityp(ia)
        do is =1,nspin,af+1
            do lmt2=1,ilmt(it)
                do lmt1=lmt2+1,ilmt(it)
                    hsr_tmp(ia,lmt1,lmt2,is)=hsr_tmp(ia,lmt2,lmt1,is)
                end do
            end do
        end do
    end do

    fi = 1.d0/nopr
    if ( charge_symm_mode >= chg_symm_level1 ) then
       fi = 1.0d0 /dble(nopr_from_fbz_to_ibz)
    endif

    hsr = 0.d0

    do iopr=1,nopr
       if ( charge_symm_mode >= chg_symm_level1 )then
         if(flg_opr_from_fbz_to_ibz(iopr) == 0 ) cycle
       endif

       do ia = 1, natm
          it = ityp(ia)
          ja=abs(ia2ia_symmtry_op_inv(ia,iopr))
          do is = 1, nspin, af+1
             do lmt1 = 1, ilmt(it)
                il1=ltp(lmt1,it)
                im1=mtp(lmt1,it)
                it1=taup(lmt1,it)
                ii=(il1-1)**2+im1
                do lmt2 = lmt1, ilmt(it)
                   il2=ltp(lmt2,it)
                   im2=mtp(lmt2,it)
                   it2=taup(lmt2,it)
                   jj=(il2-1)**2+im2

                   do n=1,nylm_paw(ii,iopr,ia)
                      iii=iylm_paw(n,ii,iopr,ia)
                      do m=1,nylm_paw(jj,iopr,ia)
                         jjj=iylm_paw(m,jj,iopr,ia)

                         do lmt3=1,ilmt(it)
                            il3=ltp(lmt3,it)
                            im3=mtp(lmt3,it)
                            it3=taup(lmt3,it)
                            kk=(il3-1)**2+im3
                            if(kk.ne.iii .or. it1.ne.it3) cycle
                            do lmt4=1,ilmt(it)
                               il4=ltp(lmt4,it)
                               im4=mtp(lmt4,it)
                               it4=taup(lmt4,it)
                               ll=(il4-1)**2+im4
                               if(ll.ne.jjj .or. it2.ne.it4) cycle

                               hsr(ia,lmt1,lmt2,is) = &
                                    hsr(ia,lmt1,lmt2,is) + &
                                    hsr_tmp(ja,lmt3,lmt4,is)* &
                                    crotylm_paw(n,ii,iopr,ia)* &
                                    crotylm_paw(m,jj,iopr,ia)

                            end do! lmt4
                         end do! lmt3

                      end do! jjj
                   end do! iii

                end do! lmt2
             end do! lmt1
          end do! is
       end do! ia
    end do! iopr

!    hsr = hsr /nopr
    hsr = hsr *fi

!ASMS modified from here 2016/09/09
    Do ia = 1, natm
       it = ityp(ia)
       Do is = 1, nspin
          Do lmt1=1, ilmt(it)
            Do lmt2=lmt1, ilmt(it)
               hsr(ia,lmt2,lmt1,is) = hsr(ia,lmt1,lmt2,is)
            End do
         End do
       End do
    End do

    if(af /= 0 .and. flg_paw) then
       iopr = nopr +af

       do ia = 1, natm
          it = ityp(ia)
          ja=abs(ia2ia_symmtry_op_inv(ia,iopr))

          if ( ja <= 0 ) cycle

          do is = 1, nspin, af+1
             do lmt1 = 1, ilmt(it)
                il1=ltp(lmt1,it)
                im1=mtp(lmt1,it)
                it1=taup(lmt1,it)
                ii=(il1-1)**2+im1

                do lmt2 = lmt1, ilmt(it)
                   il2=ltp(lmt2,it)
                   im2=mtp(lmt2,it)
                   it2=taup(lmt2,it)
                   jj=(il2-1)**2+im2

                   do n=1,nylm_paw(ii,iopr,ia)
                     iii=iylm_paw(n,ii,iopr,ia)

                      do m=1,nylm_paw(jj,iopr,ia)
                         jjj=iylm_paw(m,jj,iopr,ia)

                         do lmt3=1,ilmt(it)
                            il3=ltp(lmt3,it)
                            im3=mtp(lmt3,it)
                            it3=taup(lmt3,it)
                            kk=(il3-1)**2+im3
                            if(kk.ne.iii .or. it1.ne.it3) cycle
                            do lmt4=1,ilmt(it)
                               il4=ltp(lmt4,it)
                               im4=mtp(lmt4,it)
                               it4=taup(lmt4,it)
                               ll=(il4-1)**2+im4
                               if(ll.ne.jjj .or. it2.ne.it4) cycle

                               hsr(ia,lmt1,lmt2,nspin) = &
                                    hsr(ia,lmt1,lmt2,nspin) + &
                                    hsr(ja,lmt3,lmt4,is)* &
                                    crotylm_paw(n,ii,iopr,ia)* &
                                    crotylm_paw(m,jj,iopr,ia)

                            end do! lmt4
                         end do! lmt3

                      end do! jjj
                   end do! iii

                end do! lmt2
             end do! lmt1
          end do! is
       end do! ia

       do ia=1,natm
          it=ityp(ia)
          do is =1,nspin,af+1
             do lmt2=1,ilmt(it)
                do lmt1=lmt2+1,ilmt(it)
                   hsr(ia,lmt1,lmt2,is+1)=hsr(ia,lmt2,lmt1,is+1)
                end do
             end do
          end do
       end do
    end if
!ASMS modified to   here 2016/09/09

!    deallocate(hsr_mpi,hsr_tmp) ! MPI
    deallocate(hsr_tmp)

  end subroutine symmtrz_of_ff


#ifdef _HARDPART_3D_
!===============================================================================
  subroutine add_hardpart_to_chgq_l_3D(nfout,kspin,hsr,singlemode)
    !  The total operation number has been reduced not only for the gamma-point
    ! but also for other k-points by T. Yamasaki in April 2006.
    !  ----
    ! (Rev) T. Yamaskai, 31, Aug, 2007
    !     1. 'call set_index_arrays1' that included a bug is replaced
    !       by 'call m_PP_set_index_arrays1', whose bug is fixed.
    !     2. 'call set_index_arrays2' is also replaced by 'call
    !       m_PP_set_index_arrays2' that can be referred from other modules.
    !     3. contained subroutines, set_index_arrays1 and set_index_arrays2 were
    !       deleted.
    use m_Parallelization,     only : ista_kngp_B ,iend_kngp_B      &
  &                                 , np_kngp_B, mp_kngp_B, is_kngp_B, nel_kngp_B &
  &                                 , ista_atm_B, iend_atm_B, mp_atm_B, is_atm_B, ie_atm_B &
  &                                 , nel_atm_B, mem_atm_B, mpi_ke_world, mpi_chg_world

    integer, intent(in)      :: nfout, kspin, singlemode
    real(kind=DP), intent(in), dimension(natm,nlmt,nlmt,nspin):: hsr

    real(kind=DP), pointer, dimension(:)               :: ylm
    real(kind=DP), allocatable, target, dimension(:)   :: ylm_t
    real(kind=DP), allocatable, target, dimension(:,:) :: ylm_ext

    integer :: is,it,n,ia,mdvdb,ilm3,l3

    integer :: kngp_adj
    real(kind=DP), allocatable, target, dimension(:) :: zfcos_x, zfsin_x

    integer, allocatable, dimension(:) :: ia_list
    real(kind=DP), allocatable, dimension(:,:,:) :: shdg_x    ! d(n_ialist0,maxm,nqitg)
    real(kind=DP), allocatable, dimension(:,:,:) :: shdg_l ! d(natm,maxm,nqitg)
    real(kind=DP), allocatable, dimension(:,:,:,:) :: wk_gather
    integer :: m, maxm, iq
    integer, parameter :: mcritical = 4*2+1
    integer, allocatable, dimension(:) :: nqitg_sp, nqitg_sp0 !d(ntyp)
    integer, allocatable, dimension(:) :: iq2l3 ! d(nqitg)
    integer, allocatable, dimension(:,:) :: nc  ! d(maxm,nqitg)
    integer :: mc ! maxval(nc)
    integer, allocatable, dimension(:,:,:) :: nc2lmt1, nc2lmt2, nc2n ! d(mc,maxm,nqitg)
    integer :: ibl1,ibl2,ibsize,ncache,iwidth
    integer :: j,k,l,irank
  real(kind=DP) :: fx, fy, fz, ph
  integer       ::  iend  ,i ,iab, ip !mpi

    real(kind=DP), allocatable, dimension(:,:,:) :: chgq_add, chgq_mpi
    real(kind=DP), allocatable, dimension(:,:,:) :: chgq_sum
    real(kind=DP), allocatable, dimension(:,:) :: recvbuf
    integer :: nel, ispi_start

    if(singlemode==YES) then
       ispi_start = kspin
    else
       ispi_start = 1
    end if

    if(modnrm == EXECUT) then
       call m_PP_find_maximum_l(n)   !  n-1: maximum l
       n = (n-1) + (n-1) + 1
       allocate(il3(n**2)); il3=0; call substitute_il3(n**2,il3) ! -(b_Elec..)

       allocate(nqitg_sp(ntyp)); allocate(nqitg_sp0(ntyp))
       allocate(iq2l3(nqitg))
       allocate(nc(mcritical,nqitg));nc=0

	nqitg_sp = 0; nqitg_sp0 = 0; iq2l3 = 0

       call m_PP_set_index_arrays1(nfout,ntyp,nqitg,mcritical,n**2,il3 &
            & ,maxm,mc,nqitg_sp,nqitg_sp0,iq2l3,nc)
       allocate(nc2lmt1(mc,maxm,nqitg))
       allocate(nc2lmt2(mc,maxm,nqitg))
       allocate(nc2n(mc,maxm,nqitg))

       nc2lmt1 = 0;  nc2lmt2 = 0; nc2n = 0

       call m_PP_set_index_arrays2(nfout,mc,maxm,nqitg,mcritical,n**2,il3,iq2l3 &
            & ,nc2lmt1,nc2lmt2,nc2n,nc) ! -> nc2lmt1, nc2lmt2, nc2n, nc

!      allocate(zfcos_x(ista_kngp:iend_kngp)); zfcos_x = 0.d0
!      allocate(zfsin_x(ista_kngp:iend_kngp)); zfsin_x = 0.d0
       allocate(shdg_l(mp_atm,maxm,nqitg))
       shdg_l = 0.0d0

       if(n**2 > nel_Ylm) then
          allocate(ylm_ext(ista_kngp:iend_kngp,nel_Ylm+1:n**2)); ylm_ext = 0.d0
       end if
       allocate(ylm_t(ista_kngp:iend_kngp)); ylm_t = 0.d0
       do ilm3 = nel_Ylm+1, n**2
          call m_pwBS_sphrp2_3D(ilm3,rltv,ista_kngp,iend_kngp,ylm_t)
          ylm_ext(:,ilm3) = ylm_t(:)
       end do

       do it = 1, ntyp
          mdvdb = m_PP_include_vanderbilt_pot(it)
          if(mdvdb == SKIP) cycle

          do ia = ista_atm, iend_atm
             if(ityp(ia) == it) then
!!$                do is = 1, kspin, af+1
                do is = ispi_start, kspin, af+1
                   do iq = nqitg_sp0(it), nqitg_sp(it)
                      l3 = iq2l3(iq)
                      do m = 1, 2*l3+1
                         call sum_hsr_dot_gauntc_3D(ia,is,it,iq,m)
                         !!    hsr, dl2p -> shdg_l(ista_atm:iend_atm)
                      end do
                   end do
                end do
             end if
          end do! ia_g
       end do! it
       allocate(shdg_x(1:natm,maxm,nqitg))
!!     shdg_x = 0.0d0
       if(nrank_g > 1) then
          allocate(wk_gather(mp_atm,maxm,nqitg,0:nrank_g-1))
          call mpi_allgather(shdg_l, mp_atm*maxm*nqitg, mpi_double_precision,  &
         &                   wk_gather, mp_atm*maxm*nqitg, mpi_double_precision, mpi_ke_world, ierr)
          do irank = 0, nrank_g - 1
             do l = 1, nqitg
                do k = 1, maxm
                   do j = 1, nel_atm(irank)
                      shdg_x(is_atm(irank)+j-1,k,l) = wk_gather(j,k,l,irank)
                   end do
                end do
             end do
          end do
          deallocate(wk_gather)
       else
          shdg_x = shdg_l
       end if

#ifdef _PARA_KNGP_B_
       allocate(zfcos_x(ista_kngp_B:iend_kngp_B))
       allocate(zfsin_x(ista_kngp_B:iend_kngp_B))
!fj --------------------
!xx    allocate(chgq_add(ista_kngp_B:iend_kngp_B,kimg,nspin))
       allocate(chgq_add(ista_kngp_B:ista_kngp_B+mp_kngp_B-1,kimg,nspin))
!fj --------------------
#else
       allocate(zfcos_x(ista_kngp:iend_kngp))
       allocate(zfsin_x(ista_kngp:iend_kngp))
       allocate(chgq_add(ista_kngp:iend_kngp,kimg,nspin))
#endif
!!     zfcos_x  = 0.d0
!!     zfsin_x  = 0.d0
       chgq_add = 0.0d0

       do it = 1, ntyp
          mdvdb = m_PP_include_vanderbilt_pot(it)
          if(mdvdb == SKIP) cycle

          ncache = (m_CtrlP_cachesize()*1024)*3/4
          if(ncache == 0) then
#ifdef _PARA_KNGP_B_
             ibsize = iend_kngp_B - ista_kngp_B + 1
#else
             ibsize = iend_kngp - ista_kngp + 1
#endif
          else
             iwidth = nqitg_sp(it) - nqitg_sp0(it)
             if(kimg == 1) then ! qitg_l(i,iq)*ylm(iy)*zfsc1_1(iy)
               ibsize=ncache/(8*(2+iwidth))
             else ! qitg_l(i,iq),ylm(iy),zfsc1_1(iy),zfsc2_1(iy)
               ibsize=ncache/(8*(3+iwidth))
             end if
          end if
          do ia = 1, natm
             if(ityp(ia) == it) then

                fx = pos(ia,1)*PAI2
                fy = pos(ia,2)*PAI2
                fz = pos(ia,3)*PAI2
#ifdef _PARA_KNGP_B_
                iend = iend_kngp_B
                if( iend_kngp_B > kgp ) iend = kgp
                if( ista_kngp_B <= iend ) then
                   do i = ista_kngp_B, iend  !for mp
!                      ph = ngabc(i,1)*fx+ngabc(i,2)*fy+ngabc(i,3)*fz
                      ph = ngabc_kngp_B_l(i,1)*fx+ngabc_kngp_B_l(i,2)*fy+ngabc_kngp_B_l(i,3)*fz
                      zfcos_x(i) = dcos(ph)
                      zfsin_x(i) = dsin(ph)
                   end do
                end if
#else
                iend = iend_kngp
                if( iend_kngp > kgp ) iend = kgp
                if( ista_kngp <= iend ) then
                   do i = ista_kngp, iend  !for mp
!                      ph = ngabc(i,1)*fx+ngabc(i,2)*fy+ngabc(i,3)*fz
                      ph = ngabc_kngp_l(i,1)*fx+ngabc_kngp_l(i,2)*fy+ngabc_kngp_l(i,3)*fz
                      zfcos_x(i) = dcos(ph)
                      zfsin_x(i) = dsin(ph)
                   end do
                end if
#endif
#ifdef _PARA_KNGP_B_
                if (ista_kngp_B <= iend) then
                do ibl1=ista_kngp_B,iend,ibsize
#else
                if (ista_kngp <= iend) then
                do ibl1=ista_kngp,iend,ibsize
#endif
                   ibl2=ibl1+ibsize-1
                   if(ibl2.gt.iend) ibl2=iend

!oo                if(iprichargedensity >= 2) write(nfout,'(" !mCD:    it,  iq,  l3,   m,ilm3")')
!!$                   do is = 1, kspin, af+1
                   do is = ispi_start, kspin, af+1
                      do iq = nqitg_sp0(it), nqitg_sp(it)
                         l3 = iq2l3(iq)
                         do m = 1, 2*l3+1
                            ilm3 = l3*l3+m
!oo                         if(iprichargedensity >= 2) then
!oo                            write(nfout,'(" !mCD: ",9i5)') it, iq, l3, m, ilm3
!oo                         end if
                            if(ilm3 <= nel_Ylm) then
                               ylm => ylm_l(:,ilm3)
                            else
                               ylm => ylm_ext(:,ilm3)
                            end if
                            !! iq, shdg_x, exp(-iGR), qitg_l, ylm -> chgq_l
                            call add_hardpart_to_chgq_l_core4_3D(iq,chgq_add)
                         end do
                      end do
                   end do
                end do ! ibl1
                end if
             end if
          end do! ia
       end do! it

#ifdef _PARA_KNGP_B_
       if (nrank_e > 1) then
!         allocate(chgq_mpi(ista_kngp:iend_kngp,kimg,nspin))
!         call mpi_allreduce(chgq_add, chgq_mpi, np_kngp*kimg*nspin,                &
!        &                   mpi_double_precision, mpi_sum, mpi_kg_world, ierr)
!         chgq_l = chgq_l + chgq_mpi
!         deallocate(chgq_mpi)
          allocate(chgq_mpi(ista_kngp:iend_kngp,kimg,kspin))
          allocate(recvbuf(mp_kngp_B*kimg*kspin,0:nrank_e-1))
!fj --------------------
!xx       call mpi_allgather(chgq_add, np_kngp_B*kimg*kspin, mpi_double_precision,   &
!xx      &                  recvbuf, mp_kngp_B*kimg*kspin, mpi_double_precision, mpi_kg_world, ierr)
          call mpi_allgather(chgq_add, mp_kngp_B*kimg*kspin, mpi_double_precision,   &
         &                  recvbuf, mp_kngp_B*kimg*kspin, mpi_double_precision, mpi_kg_world, ierr)
!fj --------------------
          do irank = 0, nrank_e - 1
             nel = nel_kngp_B(irank)
!!$             do is = 1, kspin
             do is = ispi_start, kspin
                do j = 1, kimg
                   do i = 1, nel
!fj --------------------
!xx                   chgq_mpi(is_kngp_B(irank)-1+i,j,is) = &
!xx                  &    recvbuf(nel*kimg*(is-1)+nel*(j-1)+i,irank)
                      chgq_mpi(is_kngp_B(irank)-1+i,j,is) = &
                     &    recvbuf(mp_kngp_B*kimg*(is-1)+mp_kngp_B*(j-1)+i,irank)
!fj --------------------
                   end do
                end do
             end do
          end do
          chgq_l = chgq_l + chgq_mpi
          deallocate(recvbuf)
          deallocate(chgq_mpi)
       else
          chgq_l = chgq_l + chgq_add
       end if
#else
       chgq_l = chgq_l + chgq_add
#endif
       deallocate(chgq_add)

       deallocate(ylm_t)
       if(allocated(ylm_ext)) deallocate(ylm_ext)
       deallocate(shdg_x)
       deallocate(shdg_l)
       deallocate(zfsin_x,zfcos_x,il3)

       deallocate(nc2n,nc2lmt2,nc2lmt1,nc,iq2l3,nqitg_sp,nqitg_sp0)
    end if

    if(iprichargedensity >= 2) then
       write(nfout,*) ' -- before average --'
       call m_CD_wd_chgq_l_small_portion(nfout)
    endif

    if(nbztyp >= SIMPLE_CUBIC .or. af /=0 ) then
!xx    allocate(work(kgp,kimg))
! =========================================== Added by K. Tagami ====
!xx	work = 0
! ===================================================================
       if(nbztyp >= SIMPLE_CUBIC)  call charge_average_3D(NO,chgq_l)
       if(af /= 0) then
          call charge_average_3D(ANTIFERRO,chgq_l)
          if(istress == 1 .or. sw_fine_STM_simulation == ON) then
             call charge_average_3D(NO,chgsoft)         ! average of chgsoft
             call charge_average_3D(ANTIFERRO,chgsoft)
          end if
       else if(sw_fine_STM_simulation == ON) then
          call charge_average_3D(NO,chgsoft)         ! average of chgsoft
       endif

       if(af==0 .and. sw_fine_STM_simulation /= ON .and. flg_paw) &
                            call charge_average_3D(NO,chgsoft)

       if(singlemode==YES) then
          if(kspin==2) chgq_l(:,:,1) = chgq_l(:,:,kspin)
       end if

!xx    deallocate(work)
    end if
  contains

    subroutine sum_hsr_dot_gauntc_3D(ia,is,it,iq,m)
      integer, intent(in) :: ia,is,it,iq,m
      integer :: ip, lmt1, lmt2, np
      real(kind=DP) :: fac
      do ip = 1, nc(m,iq)
         lmt1 = nc2lmt1(ip,m,iq)
         lmt2 = nc2lmt2(ip,m,iq)
         np = nc2n(ip,m,iq)
         fac = 2.d0
         if(lmt1 == lmt2) fac = 1.d0
         shdg_l(ia-ista_atm+1,m,iq) = shdg_l(ia-ista_atm+1,m,iq) + &
                  & fac*iwei(ia)*hsr(ia,lmt1,lmt2,is)*dl2p(lmt1,lmt2,np,it)
      end do
    end subroutine sum_hsr_dot_gauntc_3D

    subroutine add_hardpart_to_chgq_l_core4_3D(iq,chgq)
      integer, intent(in) :: iq
#ifdef _PARA_KNGP_B_
!fj --------------------
!xx   real(kind=DP),intent(inout),dimension(ista_kngp_B:iend_kngp_B,kimg,nspin) :: chgq
     real(kind=DP),intent(inout),dimension(ista_kngp_B:ista_kngp_B+mp_kngp_B-1,kimg,nspin) :: chgq
!fj --------------------
#else
      real(kind=DP),intent(inout),dimension(ista_kngp:iend_kngp,kimg,nspin) :: chgq
#endif
      integer       :: i,iy, iyy
      real(kind=DP) :: f,f2, flchgq, zdga
      real(kind=DP), pointer, dimension(:) :: zfsc1,zfsc2

      if(mod(l3,2) == 0) then
         zdga = real(zi**(-l3))
         if(kimg == 1) then
            zfsc1 => zfcos_x(:)
         else
            f2 = -1
            zfsc1 => zfcos_x(:)
            zfsc2 => zfsin_x(:)
         end if
      else
          zdga = aimag(zi**(-l3))
         if(kimg == 1) then
            zfsc1 => zfsin_x(:)
         else
            f2 = 1
            zfsc1 => zfsin_x(:)
            zfsc2 => zfcos_x(:)
         end if
      end if
      flchgq = zdga*shdg_x(ia,m,iq)

      if(kimg == 1) then
         do i = ibl1, ibl2
            iy = i - ista_kngp+1
#ifdef _PARA_KNGP_B_
            iyy = i - ista_kngp_B+1
            chgq(i,1,is) = chgq(i,1,is)+flchgq*qitg_l(i,iq)*ylm(iy)*zfsc1(iyy)
#else
            chgq(i,1,is) = chgq(i,1,is)+flchgq*qitg_l(i,iq)*ylm(iy)*zfsc1(iy)
#endif
         end do
      else
         do i = ibl1, ibl2
            iy = i - ista_kngp+1
            f = flchgq*qitg_l(i,iq)*ylm(iy)
#ifdef _PARA_KNGP_B_
            iyy = i - ista_kngp_B+1
            chgq(i,1,is) = chgq(i,1,is) + f * zfsc1(iyy)
            chgq(i,2,is) = chgq(i,2,is) + f2 * f * zfsc2(iyy)
#else
            chgq(i,1,is) = chgq(i,1,is) + f * zfsc1(iy)
            chgq(i,2,is) = chgq(i,2,is) + f2 * f * zfsc2(iy)
#endif
         end do
      end if

    end subroutine add_hardpart_to_chgq_l_core4_3D

  end subroutine add_hardpart_to_chgq_l_3D

#else
!! else _HARDPART_3D_
!===============================================================================
  subroutine add_hardpart_to_chgq_l_3D(nfout,kspin,hsr,singlemode)
    !  The total operation number has been reduced not only for the gamma-point
    ! but also for other k-points by T. Yamasaki in April 2006.
    !  ----
    ! (Rev) T. Yamaskai, 31, Aug, 2007
    !     1. 'call set_index_arrays1' that included a bug is replaced
    !       by 'call m_PP_set_index_arrays1', whose bug is fixed.
    !     2. 'call set_index_arrays2' is also replaced by 'call
    !       m_PP_set_index_arrays2' that can be referred from other modules.
    !     3. contained subroutines, set_index_arrays1 and set_index_arrays2 were
    !       deleted.
    use m_Parallelization,     only : mpi_ke_world, mpi_chg_world &
#ifndef _NO_PARA_KNGP_B_
  &                                 , ista_kngp_B ,iend_kngp_B      &
  &                                 , np_kngp_B, mp_kngp_B, is_kngp_B, nel_kngp_B
#else
  &                                 , ista_atm_B, iend_atm_B, mp_atm_B  &
  &                                 , nel_atm_B, is_atm_B, ie_atm_B, mem_atm_B
#endif

    integer, intent(in)      :: nfout, kspin, singlemode
    real(kind=DP), intent(in), dimension(natm,nlmt,nlmt,nspin):: hsr

    real(kind=DP), pointer, dimension(:)               :: ylm
    real(kind=DP), allocatable, target, dimension(:)   :: ylm_t
    real(kind=DP), allocatable, target, dimension(:,:) :: ylm_ext

    integer :: is,it,lmt1,lmt2,n,ia,mdvdb,il1,tau1,il2,tau2,ilm3,l3,iiqitg
    real(kind=DP) :: fac !, tpos(3)

    integer :: kngp_adj, n_ialist, n_ialist0, ia_start, ia_end, n_iagroup, n_ia, ia_g
    real(kind=DP), allocatable, target, dimension(:,:) :: zfcos_x, zfsin_x

    integer, allocatable, dimension(:) :: ia_list
    real(kind=DP), allocatable, dimension(:,:) :: ylm_red, qitg_red
    real(kind=DP), allocatable, dimension(:) :: ylm_sum
!   real(kind=DP), allocatable, dimension(:,:,:) :: chgq_red
    real(kind=DP), allocatable, dimension(:,:) :: shdg  ! d(max,nqitg,nspin)
    real(kind=DP)                              :: shdg_s(2)
    integer ::          iqm, iqmmax
    real(kind=DP) :: zdga
    integer :: m, maxm, ip, np, iq, sw_spin
    integer, parameter :: mcritical = 4*2+1
    integer, allocatable, dimension(:) :: nqitg_sp, nqitg_sp0 !d(ntyp)
    integer, allocatable, dimension(:) :: iq2l3 ! d(nqitg)
    integer, allocatable, dimension(:,:) :: nc  ! d(maxm,nqitg)
    integer :: mc ! maxval(nc)
    integer, allocatable, dimension(:,:,:) :: nc2lmt1, nc2lmt2, nc2n ! d(mc,maxm,nqitg)
! NEC tune
    integer :: ibl1,ibl2,ibsize,ncache,iwidth
    real(kind=DP), allocatable, dimension(:,:,:) :: chgq_add, chgq_mpi
!fj
    real(kind=DP) :: fx, fy, fz, ph
    integer :: iy, ri, ilm, ispi_start, ispi_end
    integer :: ista, iend, mp
    real(kind=DP) :: f
!$  real(kind=DP), allocatable, dimension(:,:,:) :: chgq_th
!$  real(kind=DP), allocatable, dimension(:)     :: zfcos_th, zfsin_th, ylm_th
!$  real(kind=DP), allocatable, dimension(:,:)   :: qitg_red_th
!$  integer  OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
!$  external OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
!!!!$  include 'omp_lib.h'
!fj

#ifndef _NO_PARA_KNGP_B_
    real(kind=DP), allocatable, dimension(:,:) :: recvbuf
    integer :: nel, irank, i, j
#endif
    if(singlemode==YES) then
       if(af==0) then
          ispi_start = kspin
          ispi_end   = kspin
       else if(af/=0) then
          ispi_start = 1
          ispi_end   = 1
       end if
    else
       ispi_start = 1
       ispi_end   = kspin
    end if

                                                 __TIMER_SUB_START(722)

    if(modnrm == EXECUT) then
       if(sw_communicator_for_chg == OFF)then
       ista = ista_kngp_B
       iend = iend_kngp_B
       mp = mp_kngp_B
       else
       ista = ista_kngp
       iend = iend_kngp
       mp = mp_kngp
       endif
       call m_PP_find_maximum_l(n)   !  n-1: maximum l
       n = (n-1) + (n-1) + 1
! =============================== Modified by K. Tagami =========
!       allocate(il3(n**2)); call substitute_il3(n**2,il3) ! -(b_Elec..)
       allocate(il3(n**2)); il3=0; call substitute_il3(n**2,il3) ! -(b_Elec..)
! ==============================================================

       allocate(nqitg_sp(ntyp)); allocate(nqitg_sp0(ntyp))
       allocate(iq2l3(nqitg))
       allocate(nc(mcritical,nqitg));nc=0
! =================================== Added by K. Tagami =======
        nqitg_sp = 0; nqitg_sp0 = 0; iq2l3 = 0
! ==============================================================
       call m_PP_set_index_arrays1(nfout,ntyp,nqitg,mcritical,n**2,il3 &
            & ,maxm,mc,nqitg_sp,nqitg_sp0,iq2l3,nc)
       allocate(nc2lmt1(mc,maxm,nqitg))
       allocate(nc2lmt2(mc,maxm,nqitg))
       allocate(nc2n(mc,maxm,nqitg))
! ==================================== Added by K. Tagami ======
       nc2lmt1 = 0;  nc2lmt2 = 0; nc2n = 0
! =============================================================
       call m_PP_set_index_arrays2(nfout,mc,maxm,nqitg,mcritical,n**2,il3,iq2l3 &
            & ,nc2lmt1,nc2lmt2,nc2n,nc) ! -> nc2lmt1, nc2lmt2, nc2n, nc

       ncache = (m_CtrlP_cachesize()*1024)*3/4
       if(ncache == 0) then
#ifdef _OPENMP
          ibsize = iend - ista + 1
#else
#ifndef _NO_PARA_KNGP_B_
          ibsize = iend - ista + 1
#else
          ibsize = iend_kngp - ista_kngp + 1
#endif
#endif
       else
          iwidth = 0.d0
          do it = 1, ntyp
             iwidth = max(iwidth,nqitg_sp(it)-nqitg_sp0(it)+1)
          end do
          if(kimg == 1) then ! qitg_l(i,iq)*ylm(iy)*zfsc1_1(iy)
             ibsize=ncache/(8*(iwidth + 1))
          else ! qitg_l(i,iq),ylm(iy),zfsc1_1(iy),zfsc2_1(iy)
             ibsize=ncache/(8*(iwidth + 2))
          endif
       end if
       if(iprichargedensity >= 2) write(nfout,'(" !mCD: ibsize, iwidth = ",2i8)') ibsize, iwidth

       allocate(zfcos(ibsize))
       allocate(zfsin(ibsize))
!      zfcos = 0.d0
!      zfsin = 0.d0
       allocate(qitg_red(ibsize,nqitg))
       allocate(ylm_red(ibsize,n**2))
       allocate(ylm_sum(ibsize))
       if(kspin == 2 .and. af == 0) then
          sw_spin = ON
       else
          sw_spin = OFF
       end if
!      iqmmax = 0
!      do it = 1, ntyp
!         iqm = 0
!         do iq = nqitg_sp0(it), nqitg_sp(it)
!            l3 = iq2l3(iq)
!            do m = 1, 2*l3+1
!               iqm = iqm+1
!            end do
!         end do
!         if(iqmmax < iqm) iqmmax = iqm
!      end do

!      if(sw_spin == ON) then
!ad       allocate(chgq_red(ibsize,kimg,2))
!d        allocate(shdg(iqmmax,2))
!      else if(sw_spin == OFF) then
!ad       allocate(chgq_red(ibsize,kimg,1))
!d        allocate(shdg(iqmmax,1))
!      end if

#ifdef _OPENMP
       allocate(chgq_add(ista:ista+mp-1,kimg,nspin))
#else
#ifndef _NO_PARA_KNGP_B_
!fj --------------------
!xx    allocate(chgq_add(ista_kngp_B:iend_kngp_B,kimg,nspin))
       allocate(chgq_add(ista:ista+mp-1,kimg,nspin))
!fj --------------------
#else
       allocate(chgq_add(ista_kngp:iend_kngp,kimg,nspin))
#endif
#endif
       chgq_add = 0.0d0

                                                 __TIMER_SUB_START(727)
                                                 __TIMER_SUB_STOP(727)
                                                 __TIMER_SUB_START(728)
!                                                 __TIMER_DO_START(845)
!                                                 __TIMER_DO_STOP(845)
                                                 __TIMER_SUB_STOP(728)
                                                 __TIMER_SUB_START(729)
                                                 __TIMER_DO_START(846)
                                                 __TIMER_DO_STOP(846)
                                                 __TIMER_SUB_STOP(729)
                                                 __TIMER_SUB_START(730)
                                                 __TIMER_DO_START(847)
                                                 __TIMER_DO_STOP(847)
                                                 __TIMER_SUB_STOP(730)
                                                 __TIMER_DO_START(834)
#ifdef _OPENMP

!$OMP PARALLEL DEFAULT(NONE)        &
!$OMP          SHARED(natm,ityp,pos,          ngabc_kngp_l,kspin,af,nqitg_sp0,nqitg_sp,iq2l3,     &
!$OMP                 sw_spin,nc,nc2lmt1,nc2lmt2,nc2n,iwei,hsr,dl2p,qitg_red,ylm_red,chgq_add,    &
!$OMP                 ista,mp,nspin,kimg,ibsize,                                    &
!$OMP                 iend,kgp,nqitg,qitg_l,nel_Ylm,n,ylm_l,ylm_t,rltv,ispi_start,ispi_end)       &
!$OMP          PRIVATE(ia,it,mdvdb,fx,fy,fz,i,iy,ph,zfcos_th,zfsin_th,is,iqm,iq,l3,ylm_th,ilm3,   &
!$OMP                  shdg_s,ip,lmt1,lmt2,np,fac,f,zdga,chgq_th,ri,m,                            &
!$OMP                  ibl1,ibl2,ilm)

#ifdef __TIMER_DO__
  !$OMP BARRIER
#endif
                                                 __TIMER_DO_START(844)
       allocate(chgq_th(ista:ista+mp-1,kimg,nspin))
       allocate(zfsin_th(ibsize),zfcos_th(ibsize),ylm_th(ibsize))
#ifdef __TIMER_DO__
  !$OMP BARRIER
#endif
                                                 __TIMER_DO_STOP(844)

       do ibl1=ista,iend,ibsize
          ibl2=ibl1+ibsize-1
          if(ibl2.gt.iend) ibl2=iend
          if(ibl2.gt.kgp) ibl2 = kgp

                                                 __TIMER_SUB_START(725)
                                                 __TIMER_DO_START(840)
!$OMP DO
          do iq = 1, nqitg
             do i = 1, ibl2-ibl1+1
                qitg_red(i, iq) = qitg_l(i+ibl1-1,iq)
             end do
          end do
!$OMP END DO NOWAIT
                                                 __TIMER_DO_STOP(840)
                                                 __TIMER_SUB_STOP(725)
                                                 __TIMER_SUB_START(726)
                                                 __TIMER_DO_START(841)
!$OMP DO
          do ilm = 1, nel_Ylm
             do i = 1, ibl2-ibl1+1
                ylm_red(i,ilm) = ylm_l(i+ibl1-1,ilm)
             end do
          end do
!$OMP END DO
                                                 __TIMER_DO_STOP(841)
!$OMP MASTER
          if(n**2 > nel_Ylm) then
             allocate(ylm_t(ibl1:ibl2)); ylm_t = 0.d0
             do ilm = nel_ylm+1, n**2
                call m_pwBS_sphrp2_3D(ilm,rltv,ibl1,ibl2,ylm_t)
                                                 __TIMER_DO_START(842)
                do i = 1, ibl2-ibl1+1
                   ylm_red(i,ilm) = ylm_t(i+ibl1-1)
                end do
                                                 __TIMER_DO_STOP(842)
             end do
             deallocate(ylm_t)
          end if
!$OMP END MASTER
!$OMP BARRIER
                                                 __TIMER_SUB_STOP(726)
          chgq_th = 0.d0

!$OMP DO SCHEDULE(RUNTIME)
          do ia = 1, natm
             it = ityp(ia)
             mdvdb = m_PP_include_vanderbilt_pot(it)
             if(mdvdb == SKIP) cycle

             fx = pos(ia,1)*PAI2
             fy = pos(ia,2)*PAI2
             fz = pos(ia,3)*PAI2
             do i = 1, ibl2-ibl1+1
                iy = i + ibl1 - 1
                ph = ngabc_kngp_l(iy,1)*fx+ngabc_kngp_l(iy,2)*fy+ngabc_kngp_l(iy,3)*fz
                zfcos_th(i) = dcos(ph)
                zfsin_th(i) = dsin(ph)
             end do


!!$             do is = 1, kspin, af+1
             do is = ispi_start, ispi_end, af+1
                iqm = 0
                do iq = nqitg_sp0(it), nqitg_sp(it)
                   l3 = iq2l3(iq)
                   ylm_th = 0.d0
                   do m = 1, 2*l3+1
                      ilm3 = l3*l3+m
                      iqm = iqm+1
                      if(sw_spin == OFF) then
                         shdg_s(1) = 0.d0
                         do ip = 1, nc(m,iq)
                            lmt1 = nc2lmt1(ip,m,iq)
                            lmt2 = nc2lmt2(ip,m,iq)
                            np = nc2n(ip,m,iq)
                            fac = 2.d0; if(lmt1 == lmt2) fac = 1.d0
                            shdg_s(1) = shdg_s(1) + fac*iwei(ia)*hsr(ia,lmt1,lmt2,1)*dl2p(lmt1,lmt2,np,it)
                         end do
                      else if(sw_spin == ON) then
                         shdg_s(1) = 0.d0; shdg_s(2) = 0.d0
                         do ip = 1, nc(m,iq)
                            lmt1 = nc2lmt1(ip,m,iq)
                            lmt2 = nc2lmt2(ip,m,iq)
                            np = nc2n(ip,m,iq)
                            fac = 2.d0; if(lmt1 == lmt2) fac = 1.d0
                            shdg_s(1) = shdg_s(1) + fac*iwei(ia)*hsr(ia,lmt1,lmt2,1)*dl2p(lmt1,lmt2,np,it)
                            shdg_s(2) = shdg_s(2) + fac*iwei(ia)*hsr(ia,lmt1,lmt2,2)*dl2p(lmt1,lmt2,np,it)
                         end do
                      end if
                      do i = 1, ibl2-ibl1+1
                         ylm_th(i) = ylm_th(i) + shdg_s(is)*ylm_red(i,ilm3)
                      end do
                   end do

                   if(mod(l3,2) == 0) then
                      zdga = real(zi**(-l3))
                           !! iq, shdg_x, exp(-iGR), qitg_l, ylm -> chgq_l
                      if(kimg == 1) then
                         do i = 1, ibl2-ibl1+1
                            chgq_th(ibl1-1+i,1,is) = chgq_th(ibl1-1+i,1,is) +zdga*qitg_red(i,iq)*ylm_th(i)*zfcos_th(i)
                         end do
                      else
                         do i = 1, ibl2-ibl1+1
                            f = zdga*qitg_red(i,iq)*ylm_th(i)
                            chgq_th(ibl1-1+i,1,is) = chgq_th(ibl1-1+i,1,is) + f * zfcos_th(i)
                            chgq_th(ibl1-1+i,2,is) = chgq_th(ibl1-1+i,2,is) - f * zfsin_th(i)
                         end do
                      end if
                   else
                      zdga = aimag(zi**(-l3))
                           !! iq, shdg_x, exp(-iGR), qitg_l, ylm -> chgq_l
                      if(kimg == 1) then
                         do i = 1, ibl2-ibl1+1
                            chgq_th(ibl1-1+i,1,is) = chgq_th(ibl1-1+i,1,is) +zdga*qitg_red(i,iq)*ylm_th(i)*zfsin_th(i)
                         end do
                      else
                         do i = 1, ibl2-ibl1+1
                            f = zdga*qitg_red(i,iq)*ylm_th(i)
                            chgq_th(ibl1-1+i,1,is) = chgq_th(ibl1-1+i,1,is) +  f * zfsin_th(i)
                            chgq_th(ibl1-1+i,2,is) = chgq_th(ibl1-1+i,2,is) +  f * zfcos_th(i)
                         end do
                      end if

                   end if
                end do
             end do
          end do
!!!$OMP END DO NOWAIT
!$OMP END DO

                                                 __TIMER_DO_START(843)
!$OMP CRITICAL
!!$         do is = 1, kspin
         do is = ispi_start, ispi_end, af+1
            do ri = 1, kimg
               do i = 1, ibl2-ibl1+1
                  chgq_add(ibl1-1+i,ri,is) = chgq_add(ibl1-1+i,ri,is) + chgq_th(ibl1-1+i,ri,is)
               end do
            end do
         end do
!$OMP END CRITICAL
                                                 __TIMER_DO_STOP(843)

       end do

                                                 __TIMER_DO_START(844)
       deallocate(chgq_th,zfsin_th,zfcos_th,ylm_th)
                                                 __TIMER_DO_STOP(844)

!$OMP END PARALLEL

#else
!else ifdef _OPENMP

                                                 __TIMER_DO_START(844)
                                                 __TIMER_DO_STOP(844)
#ifndef _NO_PARA_KNGP_B_
       do ibl1=ista,iend,ibsize
          ibl2=ibl1+ibsize-1
          if(ibl2.gt.iend) ibl2=iend
#else
       do ibl1=ista_kngp,iend_kngp,ibsize
          ibl2=ibl1+ibsize-1
          if(ibl2.gt.iend_kngp) ibl2=iend_kngp
#endif
          if(ibl2.gt.kgp) ibl2 = kgp

!ad       chgq_red = 0.d0
          call substitute_qitgred()  ! qitg_l -> qitg_red
          call substitute_ylmred() ! ylm_l, ylm_ext -> ylm_red
#ifndef _NO_PARA_KNGP_B_
          do ia = 1, natm
#else
          do ia = ista_atm_B, iend_atm_B
#endif
             it = ityp(ia)
             mdvdb = m_PP_include_vanderbilt_pot(it)
             if(mdvdb == SKIP) cycle

!fj          call calc_phase_div(ia) ! -> zfsin, zfcos
             fx = pos(ia,1)*PAI2
             fy = pos(ia,2)*PAI2
             fz = pos(ia,3)*PAI2
             do i = 1, ibl2-ibl1+1
                iy = i + ibl1 - 1
                ph = ngabc_kngp_l(iy,1)*fx+ngabc_kngp_l(iy,2)*fy+ngabc_kngp_l(iy,3)*fz
                zfcos(i) = dcos(ph)
                zfsin(i) = dsin(ph)
             end do


!!$             do is = 1, kspin, af+1
!!$             do is = ispi_start, kspin, af+1
             do is = ispi_start, ispi_end, af+1
                iqm = 0
                do iq = nqitg_sp0(it), nqitg_sp(it)
                   l3 = iq2l3(iq)
                   ylm_sum = 0.d0
                   do m = 1, 2*l3+1
                      ilm3 = l3*l3+m
                      iqm = iqm+1
!fj                   call sum_hsr_dot_gauntc0(it,ia,iq,m,iqm) ! hsr, dl2p -> shdg
                      if(sw_spin == OFF) then
                         shdg_s(1) = 0.d0
                         do ip = 1, nc(m,iq)
                            lmt1 = nc2lmt1(ip,m,iq)
                            lmt2 = nc2lmt2(ip,m,iq)
                            np = nc2n(ip,m,iq)
                            fac = 2.d0; if(lmt1 == lmt2) fac = 1.d0
                            shdg_s(1) = shdg_s(1) + fac*iwei(ia)*hsr(ia,lmt1,lmt2,1)*dl2p(lmt1,lmt2,np,it)
                         end do
                      else if(sw_spin == ON) then
                         shdg_s(1) = 0.d0; shdg_s(2) = 0.d0
                         do ip = 1, nc(m,iq)
                            lmt1 = nc2lmt1(ip,m,iq)
                            lmt2 = nc2lmt2(ip,m,iq)
                            np = nc2n(ip,m,iq)
                            fac = 2.d0; if(lmt1 == lmt2) fac = 1.d0
                            shdg_s(1) = shdg_s(1) + fac*iwei(ia)*hsr(ia,lmt1,lmt2,1)*dl2p(lmt1,lmt2,np,it)
                            shdg_s(2) = shdg_s(2) + fac*iwei(ia)*hsr(ia,lmt1,lmt2,2)*dl2p(lmt1,lmt2,np,it)
                         end do
                      end if
!                     ylm_sum(:) = ylm_sum(:) + shdg_s(is)*ylm_red(:,ilm3)
                      do i = 1, ibl2-ibl1+1
                         ylm_sum(i) = ylm_sum(i) + shdg_s(is)*ylm_red(i,ilm3)
                      end do
                   end do

                   if(mod(l3,2) == 0) then
                      zdga = real(zi**(-l3))
!fj                   call add_hardpart_to_chgq_l_div0(zdga,iq)
                           !! iq, shdg_x, exp(-iGR), qitg_l, ylm -> chgq_l
                      if(kimg == 1) then
                         do i = 1, ibl2-ibl1+1
                            chgq_add(ibl1-1+i,1,is) = chgq_add(ibl1-1+i,1,is) +zdga*qitg_red(i,iq)*ylm_sum(i)*zfcos(i)
                         end do
                      else
                         do i = 1, ibl2-ibl1+1
                            f = zdga*qitg_red(i,iq)*ylm_sum(i)
                            chgq_add(ibl1-1+i,1,is) = chgq_add(ibl1-1+i,1,is) + f * zfcos(i)
                            chgq_add(ibl1-1+i,2,is) = chgq_add(ibl1-1+i,2,is) - f * zfsin(i)
                         end do
                      end if
                   else
                      zdga = aimag(zi**(-l3))
!fj                   call add_hardpart_to_chgq_l_div1(zdga,iq)
                           !! iq, shdg_x, exp(-iGR), qitg_l, ylm -> chgq_l
                      if(kimg == 1) then
                         do i = 1, ibl2-ibl1+1
                            chgq_add(ibl1-1+i,1,is) = chgq_add(ibl1-1+i,1,is) +zdga*qitg_red(i,iq)*ylm_sum(i)*zfsin(i)
                         end do
                      else
                         do i = 1, ibl2-ibl1+1
                            f = zdga*qitg_red(i,iq)*ylm_sum(i)
                            chgq_add(ibl1-1+i,1,is) = chgq_add(ibl1-1+i,1,is) +  f * zfsin(i)
                            chgq_add(ibl1-1+i,2,is) = chgq_add(ibl1-1+i,2,is) +  f * zfcos(i)
                         end do
                      end if

                   end if
                end do
             end do
          end do
!ad       call cp_chgqred2chgq()
       end do

                                                 __TIMER_DO_START(843)
                                                 __TIMER_DO_STOP(843)
#endif
!endif ifdef _OPENMP

!#ifdef USE_PROF
!  call timer_end(834)
!#endif
                                                 __TIMER_DO_STOP(834)

#ifndef _NO_PARA_KNGP_B_
       if (nrank_e > 1 .and. sw_communicator_for_chg == OFF) then
          allocate(chgq_mpi(ista_kngp:iend_kngp,kimg,kspin))
          allocate(recvbuf(mp_kngp_B*kimg*kspin,0:nrank_e-1))
                                                 __TIMER_COMM_START_w_BARRIER(mpi_kg_world,835)
!fj --------------------
!xx       call mpi_allgather(chgq_add, np_kngp_B*kimg*kspin, mpi_double_precision,   &
!xx      &                  recvbuf, mp_kngp_B*kimg*kspin, mpi_double_precision, mpi_kg_world, ierr)
          call mpi_allgather(chgq_add, mp_kngp_B*kimg*kspin, mpi_double_precision,   &
         &                  recvbuf, mp_kngp_B*kimg*kspin, mpi_double_precision, mpi_kg_world, ierr)
!fj --------------------
                                                 __TIMER_COMM_STOP(835)
                                                 __TIMER_DO_START(836)
          do irank = 0, nrank_e - 1
             nel = nel_kngp_B(irank)
             do is = 1, kspin
                do j = 1, kimg
                   do i = 1, nel
!fj --------------------
!xx                   chgq_mpi(is_kngp_B(irank)-1+i,j,is) = &
!xx                  &    recvbuf(nel*kimg*(is-1)+nel*(j-1)+i,irank)
                      chgq_mpi(is_kngp_B(irank)-1+i,j,is) = &
                     &    recvbuf(mp_kngp_B*kimg*(is-1)+mp_kngp_B*(j-1)+i,irank)
!fj --------------------
                   end do
                end do
             end do
          end do
                                                 __TIMER_DO_STOP(836)
          if(singlemode==YES.and.af==0) then
             chgq_l(:,:,kspin) = chgq_mpi(:,:,kspin)
          else
             chgq_l = chgq_l + chgq_mpi
          end if
          deallocate(recvbuf)
          deallocate(chgq_mpi)
       else
          if(singlemode==YES.and.af==0) then
             chgq_l(:,:,kspin) = chgq_add(:,:,kspin)
          else
             chgq_l = chgq_l + chgq_add
          end if
       end if
#else
       if (nrank_e > 1 .and. sw_communicator_for_chg == OFF) then
          allocate(chgq_mpi(ista_kngp:iend_kngp,kimg,nspin))
                                                 __TIMER_COMM_START_w_BARRIER(mpi_kg_world,837)
          call mpi_allreduce(chgq_add, chgq_mpi, np_kngp*kimg*nspin,                &
         &                   mpi_double_precision, mpi_sum, mpi_kg_world, ierr)
                                                 __TIMER_COMM_STOP(837)
          if(singlemode==YES.and.af==0) then
             chgq_l(:,:,kspin) = chgq_mpi(:,:,kspin)
          else
             chgq_l = chgq_l + chgq_mpi
          end if
          deallocate(chgq_mpi)
       else
          if(singlemode==YES.and.af==0) then
             chgq_l(:,:,kspin) = chgq_add(:,:,kspin)
          else
             chgq_l = chgq_l + chgq_add
          end if
       end if
#endif
       deallocate(chgq_add)

       deallocate(ylm_sum,ylm_red,qitg_red)
!d     deallocate(shdg)
!ad    deallocate(chgq_red)
       deallocate(zfsin,zfcos)

       deallocate(il3,nc2n,nc2lmt2,nc2lmt1,nc,iq2l3,nqitg_sp,nqitg_sp0)
    end if

    if(iprichargedensity >= 2) then
       write(nfout,*) ' -- before average --'
       call m_CD_wd_chgq_l_small_portion(nfout)
    endif

    if(nbztyp >= SIMPLE_CUBIC .or. af /=0 ) then
!xx    allocate(work(kgp,kimg))
! =========================================== Added by K. Tagami ====
!xx     work = 0
! ===================================================================
       if(nbztyp >= SIMPLE_CUBIC)  call charge_average_3D(NO,chgq_l)
       if(af /= 0) then
          call charge_average_3D(ANTIFERRO,chgq_l)
          if(singlemode==NO) then
             if(istress == 1 .or. sw_fine_STM_simulation == ON) then
                call charge_average_3D(NO,chgsoft)         ! average of chgsoft
                call charge_average_3D(ANTIFERRO,chgsoft)
             end if
          end if
       else if(sw_fine_STM_simulation == ON) then
          call charge_average_3D(NO,chgsoft)         ! average of chgsoft
       endif

!!$       if(af==0 .and. sw_fine_STM_simulation /= ON .and. flg_paw) &
       if(singlemode==NO .and. af==0 .and.&
            & sw_fine_STM_simulation /= ON .and. flg_paw) &
                            call charge_average_3D(NO,chgsoft)

!xx    deallocate(work)
    end if

    if(singlemode==YES) then
       if(kspin==2)  chgq_l(:,:,1) = chgq_l(:,:,2)
    end if
                                                 __TIMER_SUB_STOP(722)
  contains

    subroutine substitute_qitgred()
      integer :: iq, i
                                                 __TIMER_SUB_START(725)
                                                 __TIMER_DO_START(840)
      do iq = 1, nqitg
         do i = 1, ibl2-ibl1+1
            qitg_red(i, iq) = qitg_l(i+ibl1-1,iq)
         end do
      end do
                                                 __TIMER_DO_STOP(840)
                                                 __TIMER_SUB_STOP(725)
    end subroutine substitute_qitgred

    subroutine substitute_ylmred
      integer :: ilm, i
                                                 __TIMER_SUB_START(726)
                                                 __TIMER_DO_START(841)
      do ilm = 1, nel_Ylm
         do i = 1, ibl2-ibl1+1
            ylm_red(i,ilm) = ylm_l(i+ibl1-1,ilm)
         end do
      end do
                                                 __TIMER_DO_STOP(841)
      if(n**2 > nel_Ylm) then
         allocate(ylm_t(ibl1:ibl2)); ylm_t = 0.d0
         do ilm = nel_ylm+1, n**2
            call m_pwBS_sphrp2_3D(ilm,rltv,ibl1,ibl2,ylm_t)
                                                 __TIMER_DO_START(842)
            do i = 1, ibl2-ibl1+1
               ylm_red(i,ilm) = ylm_t(i+ibl1-1)
            end do
                                                 __TIMER_DO_STOP(842)
         end do
         deallocate(ylm_t)
      end if
                                                 __TIMER_SUB_STOP(726)
    end subroutine substitute_ylmred

    subroutine sum_hsr_dot_gauntc0(it,ia,iq,m,iqm)
      integer, intent(in) :: it,ia,iq,m,iqm
      integer :: ip, lmt1, lmt2, np
      real(kind=DP) :: fac
                                                 __TIMER_SUB_START(728)
      if(sw_spin == OFF) then
!d       shdg(iqm,1) = 0.d0
         shdg_s(1) = 0.d0
                                                 __TIMER_DO_START(844)
         do ip = 1, nc(m,iq)
            lmt1 = nc2lmt1(ip,m,iq)
            lmt2 = nc2lmt2(ip,m,iq)
            np = nc2n(ip,m,iq)
            fac = 2.d0; if(lmt1 == lmt2) fac = 1.d0
!d          shdg(iqm,1) = shdg(iqm,1) + fac*iwei(ia)*hsr(ia,lmt1,lmt2,1)*dl2p(lmt1,lmt2,np,it)
            shdg_s(1) = shdg_s(1) + fac*iwei(ia)*hsr(ia,lmt1,lmt2,1)*dl2p(lmt1,lmt2,np,it)
         end do
                                                 __TIMER_DO_STOP(844)
      else if(sw_spin == ON) then
!d       shdg(iqm,1) = 0.d0; shdg(iqm,2) = 0.d0
         shdg_s(1) = 0.d0; shdg_s(2) = 0.d0
                                                 __TIMER_DO_START(845)
         do ip = 1, nc(m,iq)
            lmt1 = nc2lmt1(ip,m,iq)
            lmt2 = nc2lmt2(ip,m,iq)
            np = nc2n(ip,m,iq)
            fac = 2.d0; if(lmt1 == lmt2) fac = 1.d0
!d          shdg(iqm,1) = shdg(iqm,1) + fac*iwei(ia)*hsr(ia,lmt1,lmt2,1)*dl2p(lmt1,lmt2,np,it)
!d          shdg(iqm,2) = shdg(iqm,2) + fac*iwei(ia)*hsr(ia,lmt1,lmt2,2)*dl2p(lmt1,lmt2,np,it)
            shdg_s(1) = shdg_s(1) + fac*iwei(ia)*hsr(ia,lmt1,lmt2,1)*dl2p(lmt1,lmt2,np,it)
            shdg_s(2) = shdg_s(2) + fac*iwei(ia)*hsr(ia,lmt1,lmt2,2)*dl2p(lmt1,lmt2,np,it)
         end do
                                                 __TIMER_DO_STOP(845)
      end if
                                                 __TIMER_SUB_STOP(728)
    end subroutine sum_hsr_dot_gauntc0

    subroutine calc_phase_div(ia)
      integer, intent(in) :: ia
      real(kind=DP) :: fx, fy, fz, ph
      integer :: i, iy
                                                 __TIMER_SUB_START(727)
      fx = pos(ia,1)*PAI2
      fy = pos(ia,2)*PAI2
      fz = pos(ia,3)*PAI2
                                                 __TIMER_DO_START(843)
      do i = 1, ibl2-ibl1+1
         iy = i + ibl1 - 1
!         ph = ngabc(iy,1)*fx+ngabc(iy,2)*fy+ngabc(iy,3)*fz
         ph = ngabc_kngp_l(iy,1)*fx+ngabc_kngp_l(iy,2)*fy+ngabc_kngp_l(iy,3)*fz
         zfcos(i) = dcos(ph)
         zfsin(i) = dsin(ph)
      end do
                                                 __TIMER_DO_STOP(843)
                                                 __TIMER_SUB_STOP(727)
    end subroutine calc_phase_div

    subroutine add_hardpart_to_chgq_l_div0(zdga,iq)
      real(kind=DP), intent(in) :: zdga
      integer, intent(in) :: iq
      integer       :: i
      real(kind=DP) :: f
                                                 __TIMER_SUB_START(729)
                                                 __TIMER_DO_START(846)
      if(kimg == 1) then
         do i = 1, ibl2-ibl1+1
            chgq_add(ibl1-1+i,1,is) = chgq_add(ibl1-1+i,1,is) +zdga*qitg_red(i,iq)*ylm_sum(i)*zfcos(i)
         end do
      else
         do i = 1, ibl2-ibl1+1
            f = zdga*qitg_red(i,iq)*ylm_sum(i)
            chgq_add(ibl1-1+i,1,is) = chgq_add(ibl1-1+i,1,is) + f * zfcos(i)
            chgq_add(ibl1-1+i,2,is) = chgq_add(ibl1-1+i,2,is) - f * zfsin(i)
         end do
      end if
                                                 __TIMER_DO_STOP(846)
                                                 __TIMER_SUB_STOP(729)
    end subroutine add_hardpart_to_chgq_l_div0

    subroutine add_hardpart_to_chgq_l_div1(zdga,iq)
      real(kind=DP), intent(in) :: zdga
      integer, intent(in) :: iq
      integer       :: i
      real(kind=DP) :: f
                                                 __TIMER_SUB_START(730)
                                                 __TIMER_DO_START(847)
      if(kimg == 1) then
         do i = 1, ibl2-ibl1+1
            chgq_add(ibl1-1+i,1,is) = chgq_add(ibl1-1+i,1,is) +zdga*qitg_red(i,iq)*ylm_sum(i)*zfsin(i)
         end do
      else
         do i = 1, ibl2-ibl1+1
            f = zdga*qitg_red(i,iq)*ylm_sum(i)
            chgq_add(ibl1-1+i,1,is) = chgq_add(ibl1-1+i,1,is) +  f * zfsin(i)
            chgq_add(ibl1-1+i,2,is) = chgq_add(ibl1-1+i,2,is) +  f * zfcos(i)
         end do
      end if
                                                 __TIMER_DO_STOP(847)
                                                 __TIMER_SUB_STOP(730)
    end subroutine add_hardpart_to_chgq_l_div1

    subroutine cp_chgqred2chgq()
!ad   integer :: i, ir

!ad   if(sw_spin == ON) then
!ad      do ir = 1, kimg
!ad         do i = ibl1, ibl2
!ad            chgq_l(i,ir,1) = chgq_l(i,ir,1) + chgq_red(i-ibl1+1,ir,1)
!ad            chgq_l(i,ir,2) = chgq_l(i,ir,2) + chgq_red(i-ibl1+1,ir,2)
!ad         end do
!ad      end do
!ad   else if(sw_spin == OFF) then
!ad      do ir = 1, kimg
!ad         do i = ibl1, ibl2
!ad            chgq_l(i,ir,1) = chgq_l(i,ir,1) + chgq_red(i-ibl1+1,ir,1)
!ad         end do
!ad      end do
!ad   end if
    end subroutine cp_chgqred2chgq

  end subroutine add_hardpart_to_chgq_l_3D
#endif

!===============================================================================
  subroutine add_hardpart_to_chgq_l_in_keworld(nfout,ispin,hsr)
    !  The total operation number has been reduced not only for the gamma-point
    ! but also for other k-points by T. Yamasaki in April 2006.
    !  ----
    ! (Rev) T. Yamaskai, 31, Aug, 2007
    !     1. 'call set_index_arrays1' that included a bug is replaced
    !       by 'call m_PP_set_index_arrays1', whose bug is fixed.
    !     2. 'call set_index_arrays2' is also replaced by 'call
    !       m_PP_set_index_arrays2' that can be referred from other modules.
    !     3. contained subroutines, set_index_arrays1 and set_index_arrays2 were
    !       deleted.
    integer, intent(in)      :: nfout, ispin
    real(kind=DP), intent(in), dimension(natm,nlmt,nlmt,nspin):: hsr

    real(kind=DP), pointer, dimension(:)               :: ylm

    integer :: is,it,lmt1,lmt2,n,ia,mdvdb,ilm3,l3
    real(kind=DP) :: fac !, tpos(3)

    real(kind=DP), allocatable, dimension(:,:) :: ylm_red, qitg_red
    real(kind=DP), allocatable, dimension(:) :: ylm_sum
    real(kind=DP), allocatable, dimension(:,:) :: shdg  ! d(max,nqitg,nspin)
    real(kind=DP)                              :: shdg_s, zdga
    integer :: m, maxm, ip, np, iq, sw_spin
    integer :: mc ! maxval(nc)
    integer, parameter :: mcritical = 4*2+1
    integer, allocatable, dimension(:) :: nqitg_sp, nqitg_sp0 !d(ntyp)
    integer, allocatable, dimension(:) :: iq2l3 ! d(nqitg)
    integer, allocatable, dimension(:,:) :: nc  ! d(maxm,nqitg)
    integer, allocatable, dimension(:,:,:) :: nc2lmt1, nc2lmt2, nc2n ! d(mc,maxm,nqitg)
! NEC tune
    integer :: ibl1,ibl2,ibsize,ncache,iwidth
    real(kind=DP), allocatable, dimension(:,:,:) :: chgq_add
    real(kind=DP), allocatable, dimension(:,:) :: recvbuf
!fj
    real(kind=DP) :: fx, fy, fz, ph, f
    integer :: iy, ri, ispi_start, ispi_end
    integer :: ista, iend, mp

    integer ::  i, j

    if(modnrm == EXECUT) then
       ista = ista_kngp
       iend = iend_kngp
!       mp = mp_kngp
       mp = np_kngp

       is = ispin
       if(af/=0) is = 1

       call m_PP_find_maximum_l(n)   !  n-1: maximum l
       n = (n-1) + (n-1) + 1
       allocate(il3(n**2)); il3=0; call substitute_il3(n**2,il3) ! -(b_Elec..)

       allocate(nqitg_sp(ntyp)); allocate(nqitg_sp0(ntyp))
       allocate(iq2l3(nqitg))
       allocate(nc(mcritical,nqitg));nc=0

       nqitg_sp = 0; nqitg_sp0 = 0; iq2l3 = 0
       call m_PP_set_index_arrays1(nfout,ntyp,nqitg,mcritical,n**2,il3 &
            & ,maxm,mc,nqitg_sp,nqitg_sp0,iq2l3,nc)
       allocate(nc2lmt1(mc,maxm,nqitg)); nc2lmt1 = 0
       allocate(nc2lmt2(mc,maxm,nqitg)); nc2lmt2 = 0
       allocate(nc2n(mc,maxm,nqitg));    nc2n = 0

! =============================================================
       call m_PP_set_index_arrays2(nfout,mc,maxm,nqitg,mcritical,n**2,il3,iq2l3 &
            & ,nc2lmt1,nc2lmt2,nc2n,nc) ! -> nc2lmt1, nc2lmt2, nc2n, nc

       ncache = (m_CtrlP_cachesize()*1024)*3/4
       if(ncache == 0) then
          ibsize = iend - ista + 1
       else
          iwidth = 0.d0
          do it = 1, ntyp
             iwidth = max(iwidth,nqitg_sp(it)-nqitg_sp0(it)+1)
          end do
          if(kimg == 1) then ! qitg_l(i,iq)*ylm(iy)*zfsc1_1(iy)
             ibsize=ncache/(8*(iwidth + 1))
          else ! qitg_l(i,iq),ylm(iy),zfsc1_1(iy),zfsc2_1(iy)
             ibsize=ncache/(8*(iwidth + 2))
          endif
       end if
       if(iprichargedensity >= 2) write(nfout,'(" !mCD: ibsize, iwidth = ",2i8)') ibsize, iwidth

       allocate(zfcos(ibsize))
       allocate(zfsin(ibsize))
!      zfcos = 0.d0
!      zfsin = 0.d0
       allocate(qitg_red(ibsize,nqitg))
       allocate(ylm_red(ibsize,n**2))
       allocate(ylm_sum(ibsize))

       allocate(chgq_add(ista:ista+mp-1,kimg,1));   chgq_add = 0.0d0

#ifdef DEBUG_LDOS
       write(6,'(" size of ylm_red = ",2i8)') ubound(ylm_red,1)-lbound(ylm_red,1)+1 &
            &                           ,ubound(ylm_red,2)-lbound(ylm_red,2)+1
       write(6,'(" size of ylm     = ",2i8)') ubound(ylm_sum,1)-lbound(ylm_sum,1)+1
       write(6,'(" ibsize = ",i8)') ibsize

       write(6,'(" size of chgq_add = ",3i8)') ubound(chgq_add,1)-lbound(chgq_add,1)+1 &
            &                           ,ubound(chgq_add,2)-lbound(chgq_add,2)+1  &
            &                           ,ubound(chgq_add,3)-lbound(chgq_add,3)+1
       write(6,'(" size of chgq_l   = ",3i8)') ubound(chgq_l,1)-lbound(chgq_l,1)+1 &
            &                           ,ubound(chgq_l,2)-lbound(chgq_l,2)+1  &
            &                           ,ubound(chgq_l,3)-lbound(chgq_l,3)+1
#endif
!!$#ifdef _OPENMP
!!$#else
!!$!else ifdef _OPENMP

       do ibl1=ista,iend,ibsize ! G-space
          ibl2=ibl1+ibsize-1
          if(ibl2.gt.iend) ibl2=iend
          if(ibl2.gt.kgp) ibl2 = kgp

#ifdef DEBUG_LDOS
          write(6,'(" ibl1, ibl2 = ",2i8)') ibl1, ibl2
          call flush(6)
#endif

          call substitute_qitgred()  ! qitg_l -> qitg_red
          call substitute_ylmred() ! ylm_l -> ylm_red

          do ia = 1, natm
             it = ityp(ia)
             mdvdb = m_PP_include_vanderbilt_pot(it)
             if(mdvdb == SKIP) cycle

             fx = pos(ia,1)*PAI2
             fy = pos(ia,2)*PAI2
             fz = pos(ia,3)*PAI2
             do i = 1, ibl2-ibl1+1
                iy = i + ibl1 - 1
                ph = ngabc_kngp_l(iy,1)*fx+ngabc_kngp_l(iy,2)*fy+ngabc_kngp_l(iy,3)*fz
                zfcos(i) = dcos(ph)
                zfsin(i) = dsin(ph)
             end do

             do iq = nqitg_sp0(it), nqitg_sp(it)
                l3 = iq2l3(iq)
                ylm_sum = 0.d0
                do m = 1, 2*l3+1
                   ilm3 = l3*l3+m
!fj                   call sum_hsr_dot_gauntc0(it,ia,iq,m,iqm) ! hsr, dl2p -> shdg
                   shdg_s = 0.d0
                   do ip = 1, nc(m,iq)
                      lmt1 = nc2lmt1(ip,m,iq)
                      lmt2 = nc2lmt2(ip,m,iq)
                      np = nc2n(ip,m,iq)
                      fac = 2.d0; if(lmt1 == lmt2) fac = 1.d0
                      shdg_s = shdg_s + fac*iwei(ia)*hsr(ia,lmt1,lmt2,is)*dl2p(lmt1,lmt2,np,it)
                   end do
!                     ylm_sum(:) = ylm_sum(:) + shdg_s(is)*ylm_red(:,ilm3)
                   do i = 1, ibl2-ibl1+1
                      ylm_sum(i) = ylm_sum(i) + shdg_s*ylm_red(i,ilm3)
                   end do
                end do

                if(mod(l3,2) == 0) then
                   zdga = real(zi**(-l3))
!fj                   call add_hardpart_to_chgq_l_div0(zdga,iq)
                           !! iq, shdg_x, exp(-iGR), qitg_l, ylm -> chgq_l
                   if(kimg == 1) then
                      do i = 1, ibl2-ibl1+1  ! G-space
                         chgq_add(ibl1-1+i,1,1) = chgq_add(ibl1-1+i,1,1) &
                              & + zdga*qitg_red(i,iq)*ylm_sum(i)*zfcos(i)
                      end do
                   else
                      do i = 1, ibl2-ibl1+1  ! G-space
                         f = zdga*qitg_red(i,iq)*ylm_sum(i)
                         chgq_add(ibl1-1+i,1,1) = chgq_add(ibl1-1+i,1,1) + f * zfcos(i)
                         chgq_add(ibl1-1+i,2,1) = chgq_add(ibl1-1+i,2,1) - f * zfsin(i)
                      end do
                   end if
                else
                   zdga = aimag(zi**(-l3))
!fj                   call add_hardpart_to_chgq_l_div1(zdga,iq)
                           !! iq, shdg_x, exp(-iGR), qitg_l, ylm -> chgq_l
                   if(kimg == 1) then
                      do i = 1, ibl2-ibl1+1  ! G-space
                         chgq_add(ibl1-1+i,1,1) = chgq_add(ibl1-1+i,1,1) +zdga*qitg_red(i,iq)*ylm_sum(i)*zfsin(i)
                      end do
                   else
                      do i = 1, ibl2-ibl1+1  ! G-space
                         f = zdga*qitg_red(i,iq)*ylm_sum(i)
                         chgq_add(ibl1-1+i,1,1) = chgq_add(ibl1-1+i,1,1) +  f * zfsin(i)
                         chgq_add(ibl1-1+i,2,1) = chgq_add(ibl1-1+i,2,1) +  f * zfcos(i)
                      end do
                   end if
                end if
             end do
          end do
       end do

!!$       write(6,'(" ! almost end of <<add_hardpart_to_chgq_l_in_keworld>>")')
!!$       call flush(6)
!!$
       chgq_l(:,:,ispin) = chgq_add(:,:,1)

       deallocate(chgq_add)

       deallocate(ylm_sum,ylm_red,qitg_red)
       deallocate(zfsin,zfcos)

       deallocate(nc2n,nc2lmt2,nc2lmt1,nc,iq2l3,nqitg_sp,nqitg_sp0)
    end if

    if(nbztyp >= SIMPLE_CUBIC .or. af /=0 ) then
!!$       write(6,'(" ! symmetrising in <<add_hardpart_to_chgq_l_in_keworld>>")')
!!$       write(6,'(" ! nbztyp = ",i8)') nbztyp
!!$       call flush(6)
!xx    allocate(work(kgp,kimg))
! =========================================== Added by K. Tagami ====
!xx     work = 0
! ===================================================================
!!$       if(nbztyp >= SIMPLE_CUBIC)  call charge_average_3D(NO,chgq_l)
!!$       if(af /= 0) then
!!$          call charge_average_3D(ANTIFERRO,chgq_l)
!!$          if(istress == 1 .or. sw_fine_STM_simulation == ON) then
!!$             call charge_average_3D(NO,chgsoft)         ! average of chgsoft
!!$             call charge_average_3D(ANTIFERRO,chgsoft)
!!$          end if
!!$       else if(sw_fine_STM_simulation == ON) then
!!$          call charge_average_3D(NO,chgsoft)         ! average of chgsoft
!!$       endif
!!$
!!$       if(af==0 .and. sw_fine_STM_simulation /= ON .and. flg_paw) &
!!$            & call charge_average_3D(NO,chgsoft)

!xx    deallocate(work)
    end if

!!$    if(singlemode==YES) then
!!$       if(kspin==2)  chgq_l(:,:,1) = chgq_l(:,:,2)
!!$    end if
  contains

    subroutine substitute_qitgred()
      integer :: iq, i
      do iq = 1, nqitg
         do i = 1, ibl2-ibl1+1
            qitg_red(i, iq) = qitg_l(i+ibl1-1,iq)
         end do
      end do
    end subroutine substitute_qitgred

    subroutine substitute_ylmred
      integer :: ilm, i
      real(kind=DP), allocatable, target, dimension(:)   :: ylm_t
      do ilm = 1, nel_Ylm
         do i = 1, ibl2-ibl1+1
            ylm_red(i,ilm) = ylm_l(i+ibl1-1,ilm)
         end do
      end do
      if(n**2 > nel_Ylm) then
         allocate(ylm_t(ibl1:ibl2)); ylm_t = 0.d0
         do ilm = nel_ylm+1, n**2
            call m_pwBS_sphrp2_3D(ilm,rltv,ibl1,ibl2,ylm_t)
            do i = 1, ibl2-ibl1+1
               ylm_red(i,ilm) = ylm_t(i+ibl1-1)
            end do
         end do
         deallocate(ylm_t)
      end if
    end subroutine substitute_ylmred
  end subroutine add_hardpart_to_chgq_l_in_keworld
!! endif _HARDPART_3D_

  subroutine m_CD_softpart_3D(nfout,kv3)
    integer, intent(in) :: nfout, kv3
    real(kind=DP), allocatable, dimension(:)  :: wf_phase
    real(kind=DP), allocatable,dimension(:,:) ::  afft_l_mpi
    real(kind=DP), allocatable,dimension(:,:,:) :: afft_l
    real(kind=DP), allocatable,dimension(:,:) :: bfft_l, bfft_l_in
    real(kind=DP), allocatable,dimension(:,:) :: map_afft_l
    integer ispin, ib1, ik, i, ip, max_elements, icolumn, istart, iend, icycle, ic, is
    integer :: isrsize, fft_l_size, mm, i1
    integer :: ib2, ibsize, ibesize, lsize
    integer :: id_sname = -1
#ifdef MPI_FFTW
    integer :: ix,iy,iz,lxh
    integer(C_INTPTR_T)  :: local_n, local_n_offset, alloc_local, lx, ly, lz, mx,my,mz
#endif
    call tstatc0_begin('m_CD_softpart_3D ',id_sname,1)
                                                 __TIMER_SUB_START(716)

#ifdef FFT_3D_DIVISION
    lsize = fft_X_x_nel*fft_X_y_nel*fft_X_z_nel
!!$#ifdef FFT_USE_SSL2_PAD
!!$    lsize = max(nel_fft_x(myrank_g),nel_fft_y(myrank_g),nel_fft_z(myrank_g))
#else
    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
#ifdef MPI_FFTW
    if(sw_mpi_fftw==ON) then
      lx = fft_box_size_WF(1,0)
      ly = fft_box_size_WF(2,0)
      lz = fft_box_size_WF(3,0)
      if(kimg==2) then
        alloc_local = fftw_mpi_local_size_3d(ly,lz,lx,mpi_ke_world,local_n,local_n_offset)
      else
        alloc_local = fftw_mpi_local_size_3d(ly,lz,lx/2,mpi_ke_world,local_n,local_n_offset)
      endif
      lsize = local_n*lx*lz
    endif
#endif
#endif

    ibsize = 1
    if (nblocksize_fftw_is_given) then
       ibsize = nblocksize_fftw
       if (ibsize < 1) ibsize = 1
    endif

    isrsize = min(lsize,mp_kngp)
    fft_l_size  = nel_fft_x(myrank_g)

#ifdef FFT_3D_DIVISION
    allocate(afft_l(lsize*2,1,nrank_s) ,stat=ierr)
    allocate(bfft_l(lsize*2,ibsize) ,stat=ierr)
    afft_l = 0.0d0
    bfft_l = 0.0d0
#else
!    allocate(afft_l(lsize*kimg,1) ,stat=ierr)
!    allocate(bfft_l(lsize*kimg,ibsize) ,stat=ierr)
    allocate(afft_l(lsize*kimg,1,nrank_s));afft_l=0.d0
    allocate(bfft_l(lsize*kimg,ibsize))
#endif

    if(ipriwf >= 3) then
       max_elements = 200
       icolumn = 10
       allocate(wf_phase(nfft/2))
! ==================================== Added by K. Tagami ========
       wf_phase = 0
! ================================================================
       icycle = ceiling(dble(min(max_elements,nfft/2))/icolumn)
    end if

#ifdef FFT_3D_DIVISION
    allocate(map_afft_l(np_kngp_gw*2,1) ,stat=ierr)
#else
    allocate(map_afft_l(np_kngp_gw*kimg,1) ,stat=ierr)
#endif
    chgq_l = 0.d0
!    do ispin = 1, nspin, af + 1
    is = 1
    do ispin = ista_spin, iend_spin, af + 1
       if(nrank_s==2) is = ispin
       if(.not.chg_has_been_calculated) then
       afft_l(:,:,is) = 0.d0
       do ik = ispin, kv3+ispin-nspin, nspin
          if(map_k(ik) /= myrank_k) cycle! MPI
          allocate(bfft_l_in(lsize*kimg,ibsize) ,stat=ierr);bfft_l_in=0.d0
#ifdef MPI_FFTW
          if(sw_mpi_fftw==ON) call m_ES_fftbox_map(ik)
#endif
          do ib1 = 1, np_e, ibsize     ! MPI
             ib2 = min(ib1+ibsize-1,np_e)
             ibesize = ib2 - ib1 + 1
#ifdef __TIMER_COMM__
             call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,bfft_l, 7) ! (swffft)
#else
#ifdef MPI_FFTW
             if(sw_mpi_fftw==ON) then
               call m_ES_WF_in_Rspace_mpifftw(ista_k,iend_k,ik,ib1,zaj_l)
             else
               call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,bfft_l) ! (swffft)
             endif
#else
             call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,bfft_l) ! (swffft)
!             call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,bfft_l_in, bfft_l) ! (swffft)
#endif
#endif
             if(ipri >= 2 .or. ipriwf >= 3) call wd_wf_phase()
#ifdef MPI_FFTW
             if(sw_mpi_fftw==ON) then
               call add_occupied_densities_mpifftw()
             else
               call add_occupied_densities()  ! -(this module) occup_l, bfft -> afft
             endif
#else
             call add_occupied_densities()  ! -(this module) occup_l, bfft -> afft
#endif
          end do
          deallocate(bfft_l_in)
       end do
       else
          afft_l(:,1,is) = chg_softpart(:,is)
       endif
#ifdef MPI_FFTW
       if(sw_mpi_fftw==ON) then
       if(kimg==2) then
         call mpi_allreduce(mpi_in_place,afft_l(:,:,is),lsize*2,mpi_double_precision,   &
         &                   mpi_sum,mpi_skg_world,ierr)
         call mpi_allreduce(mpi_in_place,afft_l(:,:,is),lsize*2,mpi_double_precision,   &
         &                   mpi_sum,mpi_sge_world,ierr)
         if(nrank_s==1) then
           afft_mpifftw=(0.d0,0.d0)
           do iy=1,local_n
             do iz=1,lz
               do ix=1,lx
                 i1 = (iy-1)*lx*lz+(iz-1)*lx+ix
                 afft_mpifftw(ix,iz,iy) = afft_l(i1*2-1,1,is)
               enddo
             enddo
           enddo
         endif
       else
         call mpi_allreduce(mpi_in_place,afft_l(:,:,is),lsize,mpi_double_precision,   &
         &                   mpi_sum,mpi_skg_world,ierr)
         call mpi_allreduce(mpi_in_place,afft_l(:,:,is),lsize,mpi_double_precision,   &
         &                   mpi_sum,mpi_sge_world,ierr)
         if(nrank_s==1) then
           afft_mpifftw_kimg1=(0.d0,0.d0)
           lxh = lx/2
           do iy=1,local_n
             do iz=1,lz
               do ix=1,lxh
                 i1 = (iy-1)*lxh*lz+(iz-1)*lxh+ix
                 afft_mpifftw_kimg1(ix,iz,iy) = afft_l(i1*2-1,1,is)
               enddo
             enddo
           enddo
         endif
       endif
       else
#ifdef FFT_3D_DIVISION
       allocate(afft_l_mpi(lsize*2,1) ,stat=ierr)
#else
       allocate(afft_l_mpi(lsize*kimg,1) ,stat=ierr)
#endif
!      if(nrank_e >= 2) then
                                                 __TIMER_COMM_START_w_BARRIER(mpi_kg_world,825)
!!$#ifdef FFT_USE_SSL2_3D
!!$          call mpi_allreduce(afft_l,afft_l_mpi,np_fft_x*kimg,mpi_double_precision,   &
#ifdef FFT_3D_DIVISION
          call mpi_allreduce(afft_l(:,:,is),afft_l_mpi,lsize*2,mpi_double_precision,   &
         &                   mpi_sum,mpi_skg_world,ierr)
! === DEBUG by tkato 2012/04/04 ================================================
          afft_l(:,:,is) = afft_l_mpi(:,:)
!!$          call mpi_allreduce(afft_l,afft_l_mpi,np_fft_x*kimg,mpi_double_precision,   &
          call mpi_allreduce(afft_l(:,:,is),afft_l_mpi,lsize*2,mpi_double_precision,   &
         &                   mpi_sum,mpi_sge_world,ierr)
! ==============================================================================
#else
!!$          call mpi_allreduce(afft_l,afft_l_mpi,np_fft_y*kimg,mpi_double_precision,   &
!!$         &                   mpi_sum,mpi_kg_world,ierr)
          if (sw_fft_xzy > 0) then
             call mpi_allreduce(afft_l(:,:,is),afft_l_mpi,np_fft_y*kimg,mpi_double_precision,   &
            &                   mpi_sum,mpi_skg_world,ierr)
! === DEBUG by tkato 2012/04/04 ================================================
             afft_l(:,:,is) = afft_l_mpi(:,:)
             call mpi_allreduce(afft_l(:,:,is),afft_l_mpi,np_fft_y*kimg,mpi_double_precision,   &
            &                   mpi_sum,mpi_sge_world,ierr)
! ==============================================================================
          else
             call mpi_allreduce(afft_l(:,:,is),afft_l_mpi,lsize*kimg,mpi_double_precision,   &
            &                   mpi_sum,mpi_skg_world,ierr)
! === DEBUG by tkato 2012/04/04 ================================================
             afft_l(:,:,is) = afft_l_mpi(:,:)
             call mpi_allreduce(afft_l(:,:,is),afft_l_mpi(:,:),lsize*kimg,mpi_double_precision,   &
            &                   mpi_sum,mpi_sge_world,ierr)
! ==============================================================================
          end if
#endif
       endif
#else
#ifdef FFT_3D_DIVISION
       allocate(afft_l_mpi(lsize*2,1) ,stat=ierr)
#else
       allocate(afft_l_mpi(lsize*kimg,1) ,stat=ierr)
#endif
!      if(nrank_e >= 2) then
                                                 __TIMER_COMM_START_w_BARRIER(mpi_kg_world,825)
!!$#ifdef FFT_USE_SSL2_3D
!!$          call mpi_allreduce(afft_l,afft_l_mpi,np_fft_x*kimg,mpi_double_precision,   &
#ifdef FFT_3D_DIVISION
          call mpi_allreduce(afft_l(:,:,is),afft_l_mpi,lsize*2,mpi_double_precision,   &
         &                   mpi_sum,mpi_skg_world,ierr)
! === DEBUG by tkato 2012/04/04 ================================================
          afft_l(:,:,is) = afft_l_mpi(:,:)
!!$          call mpi_allreduce(afft_l,afft_l_mpi,np_fft_x*kimg,mpi_double_precision,   &
          call mpi_allreduce(afft_l(:,:,is),afft_l_mpi,lsize*2,mpi_double_precision,   &
         &                   mpi_sum,mpi_sge_world,ierr)
! ==============================================================================
#else
!!$          call mpi_allreduce(afft_l,afft_l_mpi,np_fft_y*kimg,mpi_double_precision,   &
!!$         &                   mpi_sum,mpi_kg_world,ierr)
          if (sw_fft_xzy > 0) then
             call mpi_allreduce(afft_l(:,:,is),afft_l_mpi,np_fft_y*kimg,mpi_double_precision,   &
            &                   mpi_sum,mpi_skg_world,ierr)
! === DEBUG by tkato 2012/04/04 ================================================
             afft_l(:,:,is) = afft_l_mpi(:,:)
             call mpi_allreduce(afft_l(:,:,is),afft_l_mpi,np_fft_y*kimg,mpi_double_precision,   &
            &                   mpi_sum,mpi_sge_world,ierr)
! ==============================================================================
          else
             call mpi_allreduce(afft_l(:,:,is),afft_l_mpi,lsize*kimg,mpi_double_precision,   &
            &                   mpi_sum,mpi_skg_world,ierr)
! === DEBUG by tkato 2012/04/04 ================================================
             afft_l(:,:,is) = afft_l_mpi(:,:)
             call mpi_allreduce(afft_l(:,:,is),afft_l_mpi(:,:),lsize*kimg,mpi_double_precision,   &
            &                   mpi_sum,mpi_sge_world,ierr)
! ==============================================================================
          end if
#endif
#endif
                                                 __TIMER_COMM_STOP(825)
#ifdef MPI_FFTW
          if(sw_mpi_fftw/=ON) then
            afft_l(:,:,is) = afft_l_mpi(:,:)
            deallocate(afft_l_mpi)
          endif
#else
          afft_l(:,:,is) = afft_l_mpi(:,:)
          deallocate(afft_l_mpi)
#endif
!      end if
                                                 __TIMER_COMM_START(692)
#ifdef FFT_3D_DIVISION
       call m_FFT_Direct_3DIV_3D (nfout, afft_l(:,:,is), lsize, 1)
#else
#ifdef MPI_FFTW
       if(sw_mpi_fftw==ON.and.nrank_s==1) then
          call m_FFT_Direct_MPI_FFTW(nfout)
       else if (sw_mpi_fftw==OFF) then
          call m_FFT_Direct_XYZ_3D (nfout, afft_l(:,:,is), lsize, 1)
       endif
#else
       if (sw_fft_xzy > 0) then
          call m_FFT_Direct_3D (nfout, afft_l(:,:,is), lsize, 1)
       else
          call m_FFT_Direct_XYZ_3D (nfout, afft_l(:,:,is), lsize, 1)
       end if
#endif
#endif
!!$       call m_FFT_Direct_3D (nfout, afft_l, lsize, 1)
                                                 __TIMER_COMM_STOP(692)
       if(ipri >= 2 .or. ipriwf >= 3) call wd_wf_phase()
       if(nrank_s==1) then
#ifdef MPI_FFTW
         if(sw_mpi_fftw==ON) then
           call map_fft_to_chgq_mpifftw(lsize, 1, map_afft_l, nfout)
         else
           call map_fft_to_chgq_3D(lsize, 1, afft_l(:,:,is), map_afft_l, nfout)
         endif
#else
         call map_fft_to_chgq_3D(lsize, 1, afft_l(:,:,is), map_afft_l, nfout)
#endif
         call substitute_CD_for_chgq()
         !deallocate(map_afft_l)
       endif
    end do
    if(nrank_s==2) then
      call mpi_allreduce(mpi_in_place, afft_l, lsize*kimg*nrank_s, mpi_double_precision, &
                         mpi_sum, mpi_keg_world, ierr)
#ifdef MPI_FFTW
      if(sw_mpi_fftw==ON) then
        do ispin=1, nspin
          if(kimg==2) then
            afft_mpifftw=(0.d0,0.d0)
            do iy=1,local_n
              do iz=1,lz
                do ix=1,lx
                  i1 = (iy-1)*lx*lz+(iz-1)*lx+ix
                  afft_mpifftw(ix,iz,iy) = afft_l(i1*2-1,1,ispin)
                enddo
              enddo
            enddo
          else
            afft_mpifftw_kimg1=(0.d0,0.d0)
            lxh = lx/2
            do iy=1,local_n
              do iz=1,lz
                do ix=1,lxh
                  i1 = (iy-1)*lxh*lz+(iz-1)*lxh+ix
                  afft_mpifftw_kimg1(ix,iz,iy) = afft_l(i1*2-1,1,ispin)
                enddo
              enddo
            enddo
          endif
          call m_FFT_Direct_MPI_FFTW(nfout)
          call map_fft_to_chgq_mpifftw(lsize, 1, map_afft_l, nfout)
          call substitute_CD_for_chgq()
        enddo
      else
#endif
        do ispin=1, nspin
          call map_fft_to_chgq_3D(lsize, 1, afft_l(:,:,ispin), map_afft_l, nfout)
          call substitute_CD_for_chgq()
        end do
#ifdef MPI_FFTW
      endif
#endif
    endif
    deallocate(map_afft_l)
    if(ipriwf >= 3) deallocate(wf_phase)
    if(istress == 1 .or. sw_fine_STM_simulation == ON .or. flg_paw) chgsoft = chgq_l
    if(iprichargedensity >= 2) call m_CD_wd_chgq_l_small_portion(nfout)
    deallocate(afft_l)
    deallocate(bfft_l)
    if(chg_has_been_calculated) call m_ES_dealloc_chgsoft()
    call tstatc0_end(id_sname)
                                                 __TIMER_SUB_STOP(716)
  contains
    subroutine wd_wf_phase()
      ! -----------------
      if(ipri >= 2 .or. ipriwf >= 3 ) write(6,'(" !! ik = ",i8," ib1 = ",i8)') ik,ib1
      if(ipriwf >= 3) then
         ip = 0
         do i = 1, lsize-1, 2
            ip = ip + 1
            if(dabs(bfft_l(i,1)) < 1.d-12) then
               wf_phase(ip) = 0.d0
            else if( dabs(bfft_l(i+1,1)) < 1.d-12) then
               wf_phase(ip) = 0.d0
            else
               wf_phase(ip) = bfft_l(i+1,1)/bfft_l(i,1)
            end if
         end do
         istart = 1
         iend = max_elements
         do ic = 1, icycle
            iend = min(istart+icolumn-1,max_elements,lsize/2)
            write(nfout,'(" !bfft (R)   ",10d12.4)') (bfft_l(2*i-1,1),i=istart,iend)
            write(nfout,'(" !bfft (I)   ",10d12.4)') (bfft_l(2*i,1)  ,i=istart,iend)
            write(nfout,'(" !phase (nz) ",10f12.8)') (wf_phase(i),i=istart,iend)
            istart = iend+1
         end do
      end if
      ! ----------------
    end subroutine wd_wf_phase

    subroutine substitute_CD_for_chgq()
      real(kind=DP) :: fac
      integer       :: i, iend !mpi
      integer       :: ist,ien
      real(kind=DP), pointer, dimension(:,:,:) :: chgq_p
      real(kind=DP), allocatable, target, dimension(:,:,:) :: chgqtmp

      if(sw_communicator_for_chg == ON)then
          allocate(chgqtmp(ista_kngp_gw:iend_kngp_gw,kimg,nspin));chgqtmp = 0.d0
          chgq_p => chgqtmp
          ist = ista_kngp_gw
          ien = iend_kngp_gw
      else
          chgq_p => chgq_l
          ist = ista_kngp
          ien = iend_kngp
      endif
                                                 __TIMER_SUB_START(719)
      fac = 2.d0/(univol*kv3*product(fft_box_size_WF(1:3,1)))
      iend = ien
      if( iend > kg ) iend = kg
                                                 __TIMER_DO_START(830)
      if( ist <= iend ) then
         if (kimg == 1) then
            do i = ist, iend  !for mpi
               chgq_p(i,1,ispin) = map_afft_l(i-ist+1,1)*fac
            end do
         else
            do i = ist, iend  !for mpi
               chgq_p(i,1,ispin) = map_afft_l((i-ist+1)*2-1,1)*fac
               chgq_p(i,2,ispin) = map_afft_l((i-ist+1)*2  ,1)*fac
            end do
         end if
      endif
                                                 __TIMER_DO_STOP(830)
                                                 __TIMER_SUB_STOP(719)
      if(sw_communicator_for_chg == ON)then
        do i=ista_kngp,iend_kngp
           chgq_l(i,:,:) = chgq_p(i,:,:)
        enddo
        deallocate(chgqtmp)
      endif
    end subroutine substitute_CD_for_chgq

    subroutine add_occupied_densities
      integer  :: i, ib
      real(kind=DP) :: occupation
                                                 __TIMER_SUB_START(717)
                                                 __TIMER_DO_START(826)
      do ib = ib1, ib2
         occupation = occup_l(ib,ik)
         if(abs(occupation) < DELTA) cycle
#ifdef FFT_3D_DIVISION
         do i = 1, lsize*2-1, 2
            afft_l(i,1,is) = afft_l(i,1,is) + occupation*(bfft_l(i,ib-ib1+1)**2+bfft_l(i+1,ib-ib1+1)**2) ! MPI
         end do
#else
         if (sw_fft_xzy > 0) then
            do i = 1, np_fft_y*kimg-1, 2
               afft_l(i,1,is) = afft_l(i,1,is) + occupation*(bfft_l(i,ib-ib1+1)**2+bfft_l(i+1,ib-ib1+1)**2) ! MPI
            end do
         else
            do i = 1, np_fft_z*kimg-1, 2
               afft_l(i,1,is) = afft_l(i,1,is) + occupation*(bfft_l(i,ib-ib1+1)**2+bfft_l(i+1,ib-ib1+1)**2) ! MPI
            end do
         end if
#endif
         if(ipri >= 2) then
            write(nfout,'(" !cdsoft    ik, ib1 = ",2i8)') ik, ib
            write(nfout,'(" !cdsoft     occupation = ",d16.8)') occupation
            write(nfout,'(" !cdsoft     afft : ",5d16.8)') (afft_l(i,1,is),i=1,15)
         end if
      end do
                                                 __TIMER_DO_STOP(826)
                                                 __TIMER_SUB_STOP(717)
    end subroutine add_occupied_densities

#ifdef MPI_FFTW
    subroutine add_occupied_densities_mpifftw()
      integer  :: i1, ib, mm
      integer  :: ix, iy, iz
      real(kind=DP) :: occupation
      integer(C_INTPTR_T)  :: local_n, local_n_offset, alloc_local, lx, ly, lz, mx,my,mz
      integer :: lxh
      lx = fft_box_size_WF(1,0)
      ly = fft_box_size_WF(2,0)
      lz = fft_box_size_WF(3,0)
      if(kimg==2) then
        alloc_local = fftw_mpi_local_size_3d(ly,lz,lx,mpi_ke_world,local_n,local_n_offset)
      else
        alloc_local = fftw_mpi_local_size_3d(ly,lz,lx/2,mpi_ke_world,local_n,local_n_offset)
      endif
      do ib = ib1, ib2
         occupation = occup_l(ib,ik)
         if(abs(occupation) < DELTA) cycle
         if(kimg==2) then
           do iy=1,local_n
           do iz=1,lz
           do ix=1,lx
             i1=(iy-1)*lx*lz+(iz-1)*lx+ix
             afft_l(i1*2-1,1,is) = afft_l(i1*2-1,1,is)+occupation*(real(afft_mpifftw(ix,iz,iy))**2 &
                                 +                                aimag(afft_mpifftw(ix,iz,iy))**2)
           enddo
           enddo
           enddo
         else
           lxh = lx/2
           do iy=1,local_n
           do iz=1,lz
           do ix=1,lxh
             i1=(iy-1)*lxh*lz+(iz-1)*lxh+ix
             afft_l(i1*2-1,1,is) = afft_l(i1*2-1,1,is)+occupation*(real(afft_mpifftw_kimg1(ix,iz,iy))**2 &
                                 +                                aimag(afft_mpifftw_kimg1(ix,iz,iy))**2)
           enddo
           enddo
           enddo
         endif
         if(ipri >= 2) then
            write(nfout,'(" !cdsoft    ik, ib1 = ",2i8)') ik, ib
            write(nfout,'(" !cdsoft     occupation = ",d16.8)') occupation
            write(nfout,'(" !cdsoft     afft : ",5d16.8)') (afft_l(i1,1,is),i1=1,15)
         end if
      end do
    end subroutine add_occupied_densities_mpifftw
#endif

  end subroutine m_CD_softpart_3D


  subroutine m_CD_cp_chgq_to_chgqo
                                                 __TIMER_SUB_START(1225)
    chgqo_l = chgq_l
! ============================= modified by K. Tagami ======================= 5.0&11.0
!    if(flg_paw) hsro=hsr
!    if ( flg_paw .or. sw_mix_charge_hardpart == ON) then
       hsro = hsr
       if ( noncol ) hsio = hsi
!    endif
! =========================================================================== 5.0&11.0
                                                 __TIMER_SUB_STOP(1225)
  end subroutine m_CD_cp_chgq_to_chgqo

!=============================== modified by K. Tagami ====================== 5.0&11.0
  subroutine m_CD_cp_hsr_to_hsro
!    if(flg_paw) hsro=hsr
    if ( flg_paw .or. sw_mix_charge_hardpart == ON) then
       hsro = hsr
       if ( noncol ) hsio = hsi
    endif
  end subroutine m_CD_cp_hsr_to_hsro
! =========================================================================== 5.0&11.0

  subroutine m_CD_wd_chgq_l_small_portion(nfout)
    integer, intent(in) :: nfout
    integer :: ispin,i,is,ie,ri, nnspin, aaf
    real(kind=DP) :: total_charge
    if(ipri >= 1) then
       do ispin = 1, nspin, af+1
          if(nspin == 2) write(nfout,'(" !D ispin = ",i5)') ispin
          write(nfout,*) ' zchg in chgfft '
          write(nfout,'(" ! first 16 elements")')
          is = ista_kngp               ! MPI
          ie = min(is + 15, iend_kngp) ! MPI
          do ri = 1, kimg
             if(kimg == 1 .and. ri == 1) write(nfout,*) '       real part'
             if(ri == 2)                 write(nfout,*) '       imaginary part'
             write(nfout,'(" ",4d20.12)') (chgq_l(i,ri,ispin),i=is,ie) ! MPI
          end do
          write(nfout,'(" ! last 16 elements")')
          ie = iend_kngp
          is = max(ista_kngp, ie-15)
          do ri = 1, kimg
             if(kimg == 1 .and. ri == 1) write(nfout,*) '       real part'
             if(ri == 2)                 write(nfout,*) '       imaginary part'
             write(nfout,'(" ",4d20.12)') (chgq_l(i,ri,ispin),i=is,ie) ! MPI
          end do
       end do
       total_charge = 0.d0
       if(mype == 0) then
          do ispin = 1, nspin, af+1
             total_charge = total_charge + chgq_l(1,1,ispin)*univol
          end do
          write(nfout,'(" !! total_charge = ",f10.6," (m_CD_wd_chgq_l_small_portion)")') total_charge
       end if
    end if
  end subroutine m_CD_wd_chgq_l_small_portion

! =============================== added by K. Tagami ===================== 11.0
  subroutine m_CD_wd_chgq_l_portion_noncl(nfout)
    integer, intent(in) :: nfout
    integer :: ni, i,is,ie,ri
    real(kind=DP) :: total_charge

    if(ipri >= 1) then
       do ni=1, ndim_magmom

          write(nfout,*) ' zchg in chgfft '
          write(nfout,'(" ! first 16 elements")')
          is = ista_kngp               ! MPI
          ie = min(is + 15, iend_kngp) ! MPI
          do ri = 1, kimg
             if(kimg == 1 .and. ri == 1) write(nfout,*) '       real part'
             if(ri == 2)                 write(nfout,*) '       imaginary part'
             write(nfout,'(" ",4d20.12)') (chgq_l(i,ri,ni),i=is,ie) ! MPI
          end do

          write(nfout,'(" ! last 16 elements")')
          ie = iend_kngp
          is = max(ista_kngp, ie-15)
          do ri = 1, kimg
             if(kimg == 1 .and. ri == 1) write(nfout,*) '       real part'
             if(ri == 2)                 write(nfout,*) '       imaginary part'
             write(nfout,'(" ",4d20.12)') (chgq_l(i,ri,ni),i=is,ie) ! MPI
          end do

       end do
       total_charge = 0.d0
       if(mype == 0) then
          ni = 1
          total_charge = total_charge + chgq_l(1,1,ni)*univol
          write(nfout,'(" !! total_charge = ",f10.6," (m_CD_wd_chgq_l_small_portion)")') total_charge
       end if
    end if
  end subroutine m_CD_wd_chgq_l_portion_noncl
! ======================================================================== 11.0

  subroutine wd_hsr(nfout)
    integer, intent(in) :: nfout
    integer :: ispin,ia,nl1,nl2,ri
    if(ipri >= 1) then
       do ispin = 1, nspin, af+1
          if(nspin == 2) write(nfout,'(" !D ispin = ",i5)') ispin
          write(nfout,*) ' hsr '
          do ri = 1, kimg
             if(kimg == 1 .and. ri == 1) write(nfout,*) '       real part'
             if(ri == 2)                 write(nfout,*) '       imaginary part'
             do nl2 = 1, nlmt
                do nl1 = 1, nlmt
                   write(nfout,'(" ",4d20.12)') &
                        & (hsr(ia,nl1,nl2,ispin),ia=1,natm)
                end do
             end do
          end do
       end do
    end if
  end subroutine wd_hsr

! ========================================= added by K. Tagami ============ 11.0
  subroutine wd_hsr_noncl(nfout)
    integer, intent(in) :: nfout

    integer :: is, ia, nl1, nl2, ni

    if(ipri >= 1) then
      ia = 1
      write(nfout,'(" !D ia = ", i8)') ia
      do ni=1, ndim_magmom
         write(nfout,*) ' ni = ', ni
         do nl2 = 1, nlmt
           write(nfout,'(i8,14f12.4)') nl2,(hsr(ia,nl1,nl2,ni),nl1=1,min(14,nl2))
         end do
      end do
   endif
  end subroutine wd_hsr_noncl
! ==================================================================== 11.0

  subroutine m_CD_wd_hsr(nfcntn_bin_paw)
    integer, intent(in) :: nfcntn_bin_paw

!    if(mype==0) rewind nfchgt
    if(mype==0) write(nfcntn_bin_paw) hsr,hsro
! ==================== added by K. Tagami ================== 11.0
    if ( noncol ) then
       if(mype==0) write(nfcntn_bin_paw) hsi, hsio
    endif
! ========================================================== 11.0
  end subroutine m_CD_wd_hsr

  subroutine m_CD_rd_hsr(nfcntn_bin_paw)
    integer, intent(in) :: nfcntn_bin_paw

!    if(mype==0) rewind nfchgt
    if(mype==0) read(nfcntn_bin_paw) hsr,hsro
! ==================== added by K. Tagami ================== 11.0
    if ( noncol ) then
       if (mype==0) read(nfcntn_bin_paw) hsi, hsio
    endif
! ========================================================== 11.0
    call bcast_nfcntn_bin_paw

    contains

    subroutine bcast_nfcntn_bin_paw
! ============================== modified by K. Tagami ============== 11.0
!!        call mpi_bcast(hsr,natm*nlmt*nlmt*nspin,mpi_double_precision,0,MPI_CommGroup,ierr)
!!        call mpi_bcast(hsro,natm*nlmt*nlmt*nspin,mpi_double_precision,0,MPI_CommGroup,ierr)
!
      if ( noncol ) then
        call mpi_bcast( hsr, natm*nlmt*nlmt*ndim_magmom, mpi_double_precision,&
                   &    0, MPI_CommGroup, ierr )
        call mpi_bcast(hsro, natm*nlmt*nlmt*ndim_magmom, mpi_double_precision,&
                   &    0, MPI_CommGroup, ierr )
        call mpi_bcast( hsi, natm*nlmt*nlmt*ndim_magmom, mpi_double_precision,&
                   &    0, MPI_CommGroup, ierr )
        call mpi_bcast(hsio, natm*nlmt*nlmt*ndim_magmom, mpi_double_precision,&
                   &    0, MPI_CommGroup, ierr )
      else
        call mpi_bcast( hsr, natm*nlmt*nlmt*nspin, mpi_double_precision,&
                   &    0, MPI_CommGroup, ierr )
        call mpi_bcast(hsro, natm*nlmt*nlmt*nspin, mpi_double_precision,&
                   &    0, MPI_CommGroup, ierr )
      endif
! ================================================================== 11.0
    end subroutine bcast_nfcntn_bin_paw

  end subroutine m_CD_rd_hsr

  subroutine m_CD_adjust_spindensity(nfout)
!      Coded by T. Yamasaki,  17th July 2009
    integer, intent(in) :: nfout
    integer             :: i, ip

    real(kind=DP) :: total_spin_in, total_spin_out, total_charge_in, f1, f2

    if(nspin == 2 .and. sw_fix_total_spin == ON .and. af==0) then
       total_spin_in = 0.d0
       total_charge_in = 0.d0

       if(npes > 1) then
          do i = 0, npes-1
             if( is_kngp(i) <= 1 .and. 1 <= ie_kngp(i)) then
                ip = i
                exit
             end if
          end do
       end if
       i = 1
       if(ista_kngp <= i .and. i <= ista_kngp) then
          total_charge_in = chgq_l(i,1,1) + chgq_l(i,1,2)
          total_spin_in   = chgq_l(i,1,1) - chgq_l(i,1,2)
       end if
       if(npes > 1) then
          call mpi_barrier(MPI_CommGroup,ierr)
          call mpi_bcast(total_charge_in,1,mpi_double_precision,ip,MPI_CommGroup,ierr)
          call mpi_bcast(total_spin_in,  1,mpi_double_precision,ip,MPI_CommGroup,ierr)
       end if

       if(ipri >= 1) then
          write(nfout,'(" total_spin_in, total_charge_in = ",2f12.6)') &
               & total_spin_in*univol, total_charge_in*univol
          write(nfout,'(" total_spin,    total_charge,  totch = ",3f12.6)') &
               & total_spin, total_charge, totch
       end if
       if(dabs(total_spin_in - total_spin)>DELTA10) then
          if(mype == 0) then
             if(total_spin_in*total_spin > 0) then
                chgq_l(1,1,1) = (totch + total_spin)*0.5d0/univol
                chgq_l(1,1,2) = (totch - total_spin)*0.5d0/univol
             else
                chgq_l(1,1,1) = (totch - total_spin)*0.5d0/univol
                chgq_l(1,1,2) = (totch + total_spin)*0.5d0/univol
             end if
          end if
          if(ipri >= 1) then
             if(mype == 0) then
                do i = 1, nspin
                   write(nfout,'(" chgq_l(1,1,",i2,")*univol = ",f12.6)') i,chgq_l(1,1,i)*univol
                end do
             end if
          end if
       end if
    end if
  end subroutine m_CD_adjust_spindensity

! ===================================== added by K. Tagami ============== 11.0
  subroutine m_CD_adjust_spindensity_noncl(nfout)
!      Coded by T. Yamasaki,  17th July 2009
    integer, intent(in) :: nfout
    integer             :: i, ip, ni, is

    real(kind=DP) :: total_spin_in, total_spin_out, total_charge_in, f1, f2

    if ( sw_fix_total_spin == ON ) then
       total_spin_in = 0.d0
       total_charge_in = 0.d0

       if(npes > 1) then
          do i = 0, npes-1
             if( is_kngp(i) <= 1 .and. 1 <= ie_kngp(i)) then
                ip = i
                exit
             end if
          end do
       end if
       i = 1
       if(ista_kngp <= i .and. i <= ista_kngp) then
          total_charge_in = chgq_l(i,1,1)
          total_spin_in   = chgq_l(i,1,ndim_magmom)     !mz
       end if
       if(npes > 1) then
          call mpi_barrier(MPI_CommGroup,ierr)
          call mpi_bcast(total_charge_in,1,mpi_double_precision,ip,MPI_CommGroup,ierr)
          call mpi_bcast(total_spin_in,  1,mpi_double_precision,ip,MPI_CommGroup,ierr)
       end if

       if(ipri >= 1) then
          write(nfout,'(" total_spin_in, total_charge_in = ",2f12.6)') &
               & total_spin_in*univol, total_charge_in*univol
          write(nfout,'(" total_spin,    total_charge,  totch = ",3f12.6)') &
               & total_spin, total_charge, totch
       end if
       if(dabs(total_spin_in - total_spin)>DELTA10) then
          if(mype == 0) then
             if(total_spin_in*total_spin > 0) then
                chgq_l(1,1,1) = totch/univol
                chgq_l(1,1,ndim_magmom) = total_spin/univol
             else
                chgq_l(1,1,1) = totch/univol
                chgq_l(1,1,ndim_magmom) = total_spin/univol
             end if
          end if
          if(ipri >= 1) then
             if(mype == 0) then
                do is = 1, ndim_spinor
                   ni = ( is-1 )*ndim_spinor + is
                   write(nfout,'(" chgq_l(1,1,",i2,")*univol = ",f12.6)') ni,chgq_l(1,1,ni)*univol
                end do
             end if
          end if
       end if
    end if
  end subroutine m_CD_adjust_spindensity_noncl
! ===================================================================== 11.0

  subroutine m_CD_initial_CD_by_Gauss_func(nfout)
!      Revised by T. Yamasaki  Mar 9, 2007 (zeta2)
    integer, intent(in) :: nfout

    real(kind=DP) :: ratio, f1, f2, f3
    integer ::       ispin, it, i, ri, is

    integer ::       ista_kngp0
    real(kind=DP) :: sqrtpi, eta, etainv2, f4, derf
    real(kind=DP),allocatable,dimension(:) :: zeta2 !d(ntyp)

    integer ::       id_sname = -1

                                                 __TIMER_SUB_START(1226)
    call tstatc0_begin('m_CD_initial_CD_by_Gauss_func ',id_sname,1)

    chgq_l = 0.d0

    if(ipri >= 2 .and. nrank_chg == 0) then
       write(nfout,*) ' <<< m_CD_initial_CD_by_Gauss_func >>>'
       do it = 1, ntyp
          write(nfout,'(" -- zfm3_l (it = ",i5," ) --")')
          write(nfout,'(4f20.10)') (zfm3_l(i,it,1),i=1,30)
       end do
    end if

    allocate(zeta2(ntyp))
    if(sw_fix_total_spin == YES .and. nspin == 2) then
       if(dabs(total_spin) < DELTA) then
          zeta2(:) = 0.d0
       else
          f1 = 0.d0
          do it = 1, ntyp
             f1 = f1 + iatom(it)*ival(it)*zeta1(it)
          end do
          if(dabs(f1) < DELTA) then
             do it = 1, ntyp
                zeta2(it) = total_spin/(iatom(it)*ival(it)*ntyp)
             end do
          else
             f2 = total_spin/f1
             do it = 1, ntyp
                zeta2(it) = zeta1(it)*f2
             end do
          end if
       end if
    else
       zeta2 = zeta1
    end if

    if(initial_chg == Gauss_distrib_func) then
       if(ipri>=1) write(nfout,*) 'Charge density initialization: Gauss distrib func'
       do ispin = 1, nspin
          do it = 1, ntyp
             ratio = (nspin-1)*(1+(3-2*ispin)*zeta2(it))*0.5d0 + (2-nspin) ! =1 (when nspin == 1)
             f1 = 0.25d0/alfa(it)
             f2 = ratio * (ival(it)+ (qex(it)/iatom(it)) )/univol

             do ri = 1, kimg
                do i = ista_kngp, iend_kngp      !for mpi
                   f3 = f2*dexp(-gr_l(i)**2*f1)
                   chgq_l(i,ri,ispin) = chgq_l(i,ri,ispin) + zfm3_l(i,it,ri)*f3
                end do
             end do
          end do
       end do
    else if(initial_chg == VERY_NARROW) then ! initial_chg == very_broad
       if(ipri>=1) write(nfout,*) 'Charge density initialization: Very narrow func'
       sqrtpi = dsqrt(PAI)
       ista_kngp0 = ista_kngp
       if(ista_kngp0 == 1) ista_kngp0 = 2
       do ispin = 1, nspin
          do it = 1, ntyp
             eta = dsqrt(alfa(it))
             etainv2 = 0.5/eta
             ratio = (nspin-1)*(1+(3-2*ispin)*zeta2(it))*0.5d0 + (2-nspin) ! =1 (when nspin == 1)
!!$             f1 = 0.25d0/alfa(it)
             f2 = ratio * (ival(it)+ (qex(it)/iatom(it)) )/univol
             f4 = f2 * sqrtpi*eta
             do ri = 1, kimg
                do i = ista_kngp0, iend_kngp      !for mpi
!!$                   f3 = f2*dexp(-gr_l(i)**2*f1)
                   f3 = f4*derf(gr_l(i)*etainv2)/gr_l(i)
                   chgq_l(i,ri,ispin) = chgq_l(i,ri,ispin) + zfm3_l(i,it,ri)*f3
                end do
             end do
!!$          if(ista_kngp == 1) chgq_l(1,1,ispin) = chgq_l(1,1,ispin) + f2
             if(ista_kngp == 1) chgq_l(1,1,ispin) = chgq_l(1,1,ispin) + zfm3_l(1,it,1)*f2
          end do
       end do
    else if(initial_chg == from_PseudoPotential_FILE) then
       if(ipri>=1) write(nfout,*) 'Charge density initialization: from_PseudoPotential_FILE(atomic charge density)'
       ista_kngp0 = ista_kngp
       if(ista_kngp0 == 1) ista_kngp0 = 2
       do ispin = 1, nspin
          do it = 1, ntyp
             ratio = (nspin-1)*(1+(3-2*ispin)*zeta2(it))*0.5d0 + (2-nspin) ! =1 (when nspin == 1)
             f2 = ratio * ( 1.d0 + qex(it)/(iatom(it)*ival(it)) )
             do ri = 1, kimg
                do i = ista_kngp0, iend_kngp      !for mpi
                   chgq_l(i,ri,ispin) = chgq_l(i,ri,ispin) + zfm3_l(i,it,ri)*f2*rhvg_l(i,it)
                end do
             end do
             if(ista_kngp == 1) then
                f2 = ratio * (ival(it)+ (qex(it)/iatom(it)) )/univol
                chgq_l(1,1,ispin) = chgq_l(1,1,ispin) + zfm3_l(1,it,1)*f2
             endif
          end do
       end do
    end if
    deallocate(zeta2)

    if(ipri >= 2 .and. myrank_chg == 0) then
       do ispin = 1, nspin
! ============================= modified by K. Tagami ========
!!          if(kimg == 1) then
!             write(nfout,'(" -- chgq_l --")')  !mpi
!          else
!             if(ri == 1) write(nfout,'(" -- chgq_l -- (real part)")')
!             if(ri == 2) write(nfout,'(" -- chgq_l -- (imag part)")')
!             write(nfout,'(30f20.10)') (chgq_l(i,ri,ispin),i=1,30)
!          end if
          do ri = 1, kimg
             if(kimg == 1) then
                write(nfout,'(" -- chgq_l --")')  !mpi
             else
                if(ri == 1) write(nfout,'(" -- chgq_l -- (real part)")')
                if(ri == 2) write(nfout,'(" -- chgq_l -- (imag part)")')
             end if
             write(nfout,'(4f20.10)') (chgq_l(i,ri,ispin),i=1,30)
          end do
! ===========================================================
       end do
    end if

    total_charge = 0.d0
    if(myrank_chg == 0) then
       do is = 1, nspin, af+1
          total_charge = total_charge + chgq_l(1,1,is)*univol
       end do
    end if
    if(nrank_chg > 1) call mpi_bcast(total_charge,1,mpi_double_precision,0,mpi_chg_world,ierr)
    if(ipritotalcharge >= 2) &
         & write(nfout,'(" !! total_charge = ",f15.6," <<m_CD_initial_CD_by_Gauss_func>>")') total_charge

    if(af /= 0) then
       call charge_average_3D(ANTIFERRO,chgq_l)
    endif

    call tstatc0_end(id_sname)
                                                 __TIMER_SUB_STOP(1226)
  end subroutine m_CD_initial_CD_by_Gauss_func

! ================================= added by K. Tagami ================== 11.0&13.0U
  subroutine m_CD_initial_CD_by_Gauss_kt(nfout)
    integer, intent(in) :: nfout
    integer ::       ista_kngp0
    real(kind=DP),allocatable,dimension(:) :: zeta2, ival2 !d(ntyp)
    logical :: ionize_strict = .false.

    integer ::       id_sname = -1
! --------------------------------- start ---------------
    call tstatc0_begin('m_CD_initial_CD_by_Gauss_func_kt ',id_sname,1)
    chgq_l = 0.d0

    if ( mag_moment0_atoms_is_defined ) then
       select case ( initial_chg )
       case (Gauss_distrib_func)
          call goto_Gauss_distrib_func2
       case (from_PseudoPotential_FILE)
          call goto_from_PseudoPot_file2
       end select

    else
       call print_header

       allocate(ival2(ntyp)); ival2 = 0.0d0
       allocate(zeta2(ntyp)); zeta2 = 0.0d0

       call set_val_ival2
       call reset_val_zeta1
       call set_val_zeta2

       call print_ival2_zeta2

       select case ( initial_chg )
       case (Gauss_distrib_func)
          call goto_Gauss_distrib_func
       case (VERY_NARROW)
          call goto_VERY_NARROW
       case (from_PseudoPotential_FILE)
          call goto_from_PseudoPot_file
       end select

       deallocate(ival2);  deallocate(zeta2)
    endif

    call print_chgql
    call print_total_charge
    if ( noncol ) call print_magnetic_moment
    if(af/=0) then
       call charge_average_3D(ANTIFERRO,chgq_l)
    endif
    call tstatc0_end(id_sname)

  contains

    subroutine print_header
      integer :: it, i

      if(ipri >= 2 .and. mype == 0) then
         write(nfout,*) ' <<< m_CD_initial_CD_by_Gauss_kt >>>'
         do it = 1, ntyp
            write(nfout,'(" -- zfm3_l (it = ",i5," ) --")')
            write(nfout,'(4f20.10)') (zfm3_l(i,it,1),i=1,30)
         end do
      end if
    end subroutine print_header

    subroutine print_chgql
      integer :: is, ni, ri, i

      if(ipri >= 2 .and. mype == 0) then
         if ( noncol ) then
           do is = 1, ndim_magmom
              write(nfout,*) '! --- is = ', is
             do ri = 1, kimg
                if(kimg == 1) then
                   write(nfout,'(" -- chgq_l --")')  !mpi
                else
                   if(ri == 1) write(nfout,'(" -- chgq_l -- (real part)")')
                   if(ri == 2) write(nfout,'(" -- chgq_l -- (imag part)")')
                end if
                write(nfout,'(4f20.10)') (chgq_l(i,ri,is),i=1,30)
             end do
           end do
         else
           do is = 1, nspin
             do ri = 1, kimg
                if(kimg == 1) then
                   write(nfout,'(" -- chgq_l --")')  !mpi
                else
                   if(ri == 1) write(nfout,'(" -- chgq_l -- (real part)")')
                   if(ri == 2) write(nfout,'(" -- chgq_l -- (imag part)")')
                end if
                write(nfout,'(4f20.10)') (chgq_l(i,ri,is),i=1,30)
             end do
           end do
         endif
      endif
    end subroutine print_chgql

! ================================ KT_add ===================== 13.0U
    subroutine print_ival2_zeta2
      integer :: it

      if ( ipri < 2 ) return

      write(nfout,*) '------------ info ival2 zeta2 --------'
      write(nfout,*) '     id     ival      ival2     zeta2'
      Do it=1, ntyp
         write(nfout,'(I8,2F10.4,F12.6)') it, ival(it), ival2(it), zeta2(it)
      End Do
      write(nfout,*) '--------------------------------------'

    end subroutine print_ival2_zeta2

    subroutine set_val_ival2
      integer :: it, ia
      real(kind=DP) :: csum1, csum2, csum3, cfactor, c1

      csum1 = 0.0d0
      Do it=1, ntyp
         csum1 = csum1 + iatom(it) *ionic_charge_atomtyp(it)
      End do
      write(nfout,*) '!! total ionic_charge is ', csum1, 'e'

      csum2 = 0.0d0
      Do it=1, ntyp
         csum2 = csum2 + qex(it)
      End do
      csum2 = csum2 -additional_charge

      write(nfout,*) '!! extra charge is ', -csum2, 'e'

      csum2 = csum2 / dble(natm)

      if ( ionize_strict ) then
         cfactor = 1.0D0
      else
         cfactor = 0.9D0
      endif

      Do it=1, ntyp
         ival2(it) = ival(it) - ionic_charge_atomtyp(it)*cfactor +csum2
      End do

! -- redistribute remaining charge
      if ( .not. ionize_strict ) then
         csum3 = 0.0d0
         Do ia=1, natm
            it = ityp(ia)
            csum3 = csum3 +ival2(it)*iwei(ia)
         End do
         csum3 = ( totch -csum3 )/dble(natm2)
         Do it=1, ntyp
            ival2(it) = ival2(it) +csum3
         End do
      endif

    end subroutine set_val_ival2

    subroutine reset_val_zeta1
      integer :: it
      real(kind=DP) :: cnorm, f1, nelec_up, nelec_down
      real(kind=DP), parameter :: criterion = 1.0D-2

      if ( mag_moment0_atomtyp_is_defined ) then
         if ( noncol ) then
            Do it=1, ntyp
               cnorm = mag_moment0_atomtyp(it,1)**2 &
                    & +mag_moment0_atomtyp(it,2)**2 &
                    & +mag_moment0_atomtyp(it,3)**2
               cnorm = sqrt( cnorm )
!
               if ( ival2(it) > criterion .and.  cnorm > criterion ) then
                  if ( zeta1(it) >= 0.0d0 ) then
                     zeta1(it) = cnorm /ival2(it)
                  else
                     zeta1(it) = -cnorm /ival2(it)
                  endif
               else
                  zeta1(it) = 0.0d0
               endif
            End do

         else
            Do it=1, ntyp
               cnorm = mag_moment0_atomtyp(it,1)
!
               if ( ival2(it) > criterion .and.  abs(cnorm) > criterion ) then
                  zeta1(it) = cnorm /ival2(it)
               else
                  zeta1(it) = 0.0d0
               endif
            End do
         endif
      endif
! ----------------------------------
      f1 = total_spin /dble(natm) /2.0d0
!
      Do it=1, ntyp
         nelec_up   = ival2(it)/2.0d0 *( 1.0d0 +zeta1(it) )
         nelec_down = ival2(it)/2.0d0 *( 1.0d0 -zeta1(it) )

         nelec_up   = nelec_up    + f1
         nelec_down = nelec_down - f1

         if ( ival2(it) > criterion ) then
            zeta1(it) = ( nelec_up -nelec_down )/ dble( ival2(it) )
         endif
      End do

    end subroutine reset_val_zeta1
! =================================================== 13.0U

    subroutine set_val_zeta2
      integer :: it
      real(kind=DP) :: f1, f2
      real(kind=DP), parameter :: criterion = 1.0D-2
      logical :: flag = .false.

      if ( noncol ) then
        if ( sw_fix_total_spin == YES ) then
          flag = .true.
        endif
      else
        if ( sw_fix_total_spin == YES .and. nspin == 2) then
          flag = .true.
        endif
      endif
      if ( flag ) then
         if(dabs(total_spin) < DELTA) then
            zeta2(:) = 0.d0
         else
            f1 = 0.d0
            do it = 1, ntyp
               f1 = f1 + iatom(it)*ival2(it)*zeta1(it)
            end do

            if(dabs(f1) < DELTA) then
               do it = 1, ntyp
                  if ( ival2(it) > criterion ) then
                     zeta2(it) = total_spin/(iatom(it)*ival2(it)*ntyp)
                  endif
               end do
            else
               f2 = total_spin/f1
               do it = 1, ntyp
                  zeta2(it) = zeta1(it)*f2
               end do
            end if
         end if
      else
         zeta2 = zeta1
      end if
    end subroutine set_val_zeta2

    subroutine goto_Gauss_distrib_func
      integer :: is, it, ri, i
      real(kind=DP) :: f1, f2, f3, ratio

      if ( ipri>=1 ) write(nfout,*) 'Charge density initialization: Gauss distrib func'

      if ( noncol ) then
         if ( ndim_spinor /=2 ) then
            write(*,*) 'Not supported : ndim_spinor /=2 '
            call phase_error_with_msg(nfout,'Not supported : ndim_spinor /=2 ',__LINE__,__FILE__)
         endif
         do is = 1, ndim_magmom
            do it = 1, ntyp
               if ( is == 1 ) then
                 ratio = 1.0
               else
                 ratio = zeta2(it) *mag_direction0_atomtyp(it,is-1)
               endif

               f1 = 0.25d0/alfa(it)
!               f2 = ratio * (ival(it)+ (qex(it)/iatom(it)) )/univol
               f2 = ratio * ival2(it) /univol

               do ri = 1, kimg
                  do i = ista_kngp, iend_kngp      !for mpi
                     f3 = f2*dexp(-gr_l(i)**2*f1)
                     chgq_l(i,ri,is) = chgq_l(i,ri,is) +zfm3_l(i,it,ri)*f3
                  end do
               end do
            end do
         end do
      else
         do is = 1, nspin
            do it = 1, ntyp
               ratio = (nspin-1)*(1+(3-2*is)*zeta2(it))*0.5d0 + (2-nspin)
                                                       ! =1 (when nspin == 1)
               f1 = 0.25d0/alfa(it)
!!               f2 = ratio * (ival(it)+ (qex(it)/iatom(it)) )/univol
               f2 = ratio * ival2(it) /univol

               do ri = 1, kimg
                  do i = ista_kngp, iend_kngp      !for mpi
                     f3 = f2*dexp(-gr_l(i)**2*f1)
                     chgq_l(i,ri,is) = chgq_l(i,ri,is) + zfm3_l(i,it,ri)*f3
                  end do
               end do
            end do
         end do
       endif
    end subroutine goto_Gauss_distrib_func

    subroutine goto_Gauss_distrib_func2
      integer :: is, it, ri, i, ia, t
      real(kind=DP) :: f1, f2, f3, ratio, c1, s1
      real(kind=DP) :: grt, nchg(4), ncharge

      if ( ipri>=1 ) write(nfout,*) 'Charge density initialization: Gauss distrib func2'
! --
      if ( noncol ) then
         if ( ndim_spinor /=2 ) then
            write(*,*) 'Not supported: ndim_spinor /=2'
            call phase_error_with_msg(nfout,'Not supported : ndim_spinor /=2 ',__LINE__,__FILE__)
         endif
         do ia=1, natm
            it = ityp(ia)

            ncharge = ival(it) +ionic_charge_atoms(ia)
            nchg(1) = ncharge
            nchg(2) = mag_moment0_atoms(ia,1)
            nchg(3) = mag_moment0_atoms(ia,2)
            nchg(4) = mag_moment0_atoms(ia,3)

            f1 = 0.25d0/alfa(it)
            Do is=1, ndim_magmom
               f2 = nchg(is) /univol
               do i = ista_kngp, iend_kngp  !for mpi
                  f3 = f2*dexp(-gr_l(i)**2*f1)
                  grt = pos(ia,1)*ngabc(i,1) + pos(ia,2)*ngabc_kngp_l(i,2) &
                       &                     + pos(ia,3)*ngabc_kngp_l(i,3)
                  grt = grt *PAI2
                  if ( kimg == 1 ) then
                     c1 = iwei(ia) *cos(grt);
                     chgq_l(i,1,is) = chgq_l(i,1,is) + c1 *f3
                  else
                     c1 = cos(grt);  s1 = -sin(grt)
                     chgq_l(i,1,is) = chgq_l(i,1,is) + c1 *f3
                     chgq_l(i,2,is) = chgq_l(i,2,is) + s1 *f3
                  endif
               end do
            end do
         end do

      else
         do ia = 1, natm
            it = ityp(ia)

            ncharge = ival(it) +ionic_charge_atoms(ia)
            if ( nspin == 1 ) then
               nchg(1) = ncharge
            else
               nchg(1) = ( ncharge +mag_moment0_atoms(ia,1) ) /2.0d0
               nchg(2) = ( ncharge -mag_moment0_atoms(ia,1) ) /2.0d0
            endif

            f1 = 0.25d0/alfa(it)
            Do is=1, nspin
               f2 = nchg(is) /univol
               do i = ista_kngp, iend_kngp  !for mpi
                  f3 = f2*dexp(-gr_l(i)**2*f1)
                  grt = pos(ia,1)*ngabc(i,1) + pos(ia,2)*ngabc_kngp_l(i,2) &
                       &                     + pos(ia,3)*ngabc_kngp_l(i,3)
                  grt = grt *PAI2
                  if ( kimg == 1 ) then
                     c1 = iwei(ia) *cos(grt);
                     chgq_l(i,1,is) = chgq_l(i,1,is) + c1 *f3
                  else
                     c1 = cos(grt);  s1 = -sin(grt)
                     chgq_l(i,1,is) = chgq_l(i,1,is) + c1 *f3
                     chgq_l(i,2,is) = chgq_l(i,2,is) + s1 *f3
                  endif
               end do
            end do
         end do
      endif
    end subroutine goto_Gauss_distrib_func2

    subroutine goto_very_narrow
      integer :: is, it, ri, i
      real(kind=DP) :: f1, f2, f3, ratio
      real(kind=DP) :: sqrtpi, eta, etainv2, f4
      real(kind=DP) :: derf

      if(ipri>=1) write(nfout,*) 'Charge density initialization: Very narrow func'
      sqrtpi = dsqrt(PAI)
      ista_kngp0 = ista_kngp
      if(ista_kngp0 == 1) ista_kngp0 = 2

      if ( noncol ) then
         if ( ndim_spinor /=2 ) then
            write(*,*) 'Not supported : ndim_spinor /=2 '
            call phase_error_with_msg(nfout,'Not supported : ndim_spinor /=2 ',__LINE__,__FILE__)
         endif
         do is = 1, ndim_magmom
            do it = 1, ntyp
               eta = dsqrt(alfa(it));    etainv2 = 0.5/eta

               if ( is == 1 ) then
                 ratio = 1.0
               else
                 ratio = zeta2(it) *mag_direction0_atomtyp(it,is-1)
               endif

!!!!!               f2 = ratio * (ival(it)+ (qex(it)/iatom(it)) )/univol
               f2 = ratio * ival2(it)/univol

               f4 = f2 * sqrtpi*eta
               do ri = 1, kimg
                  do i = ista_kngp0, iend_kngp      !for mpi
!!$                   f3 = f2*dexp(-gr_l(i)**2*f1)
                     f3 = f4*derf(gr_l(i)*etainv2)/gr_l(i)
                     chgq_l(i,ri,is) = chgq_l(i,ri,is) +zfm3_l(i,it,ri)*f3
                  end do
               end do
               if ( ista_kngp == 1) then
                  chgq_l(1,1,is) = chgq_l(1,1,is) + zfm3_l(1,it,1)*f2
               endif
            end do
         end do
      else
         do is = 1, nspin
            do it = 1, ntyp
               eta = dsqrt(alfa(it));  etainv2 = 0.5/eta
               ratio = (nspin-1)*(1+(3-2*is)*zeta2(it))*0.5d0 + (2-nspin)
                                                 ! =1 (when nspin == 1)
!!$             f1 = 0.25d0/alfa(it)
!!!               f2 = ratio * (ival(it)+ (qex(it)/iatom(it)) )/univol
               f2 = ratio * ival2(it)/univol

               f4 = f2 * sqrtpi*eta
               do ri = 1, kimg
                  do i = ista_kngp0, iend_kngp      !for mpi
!!$                   f3 = f2*dexp(-gr_l(i)**2*f1)
                     f3 = f4*derf(gr_l(i)*etainv2)/gr_l(i)
                     chgq_l(i,ri,is) = chgq_l(i,ri,is) + zfm3_l(i,it,ri)*f3
                  end do
               end do
!!$          if(ista_kngp == 1) chgq_l(1,1,ispin) = chgq_l(1,1,ispin) + f2
               if ( ista_kngp == 1) chgq_l(1,1,is) = chgq_l(1,1,is) &
                    &                                 + zfm3_l(1,it,1)*f2
            end do
         end do
       endif
    end subroutine goto_very_narrow

    subroutine goto_from_PseudoPot_file
      integer :: is, it, ri, i
      real(kind=DP) :: f1, f2, f3, ratio

      if (ipri>=1) then
         write(nfout,*) 'Charge density initialization: from_PseudoPotential_FILE(atomic charge density)'
      endif

      ista_kngp0 = ista_kngp
      if(ista_kngp0 == 1) ista_kngp0 = 2

      if ( noncol ) then
         if ( ndim_spinor /=2 ) then
            write(*,*) 'Not supported: ndim_spinor /=2'
            call phase_error_with_msg(nfout,'Not supported : ndim_spinor /=2 ',__LINE__,__FILE__)
         endif
         do is = 1, ndim_magmom
            do it = 1, ntyp
               if ( is == 1 ) then
                 ratio = 1.0
               else
                 ratio = zeta2(it) *mag_direction0_atomtyp(it,is-1)
               endif

               f1 = 0.25d0/alfa(it)
!!!               f2 = ratio * ( 1.d0 + qex(it)/(iatom(it)*ival(it)) )
               f2 = ratio * ival2(it) /ival(it)

               do ri = 1, kimg
                  do i = ista_kngp0, iend_kngp      !for mpi
                     chgq_l(i,ri,is) = chgq_l(i,ri,is) +zfm3_l(i,it,ri)*f2*rhvg_l(i,it)

! -- used for debug --
!                     f3 = f2 *dexp( -gr_l(i)**2 *f1 )
!                     chgq_l(i,ri,is) = chgq_l(i,ri,is) +zfm3_l(i,it,ri)*f3*rhvg_l(i,it)
! ---
                  end do
               end do
               if ( ista_kngp == 1) then
                  f2 = ratio * ival2(it)/univol
                  chgq_l(1,1,is) = chgq_l(1,1,is) + zfm3_l(1,it,1)*f2
               endif
            end do
         end do

! -- for debug--
!         chgq_l(:,:,2:4) = 0.0d0
! --
      else
         do is = 1, nspin
            do it = 1, ntyp
               ratio = (nspin-1)*(1+(3-2*is)*zeta2(it))*0.5d0 + (2-nspin)
                                                   ! =1 (when nspin == 1)
!               f2 = ratio * ( 1.d0 + qex(it)/(iatom(it)*ival(it)) )
               f2 = ratio * ival2(it) /ival(it)

               do ri = 1, kimg
                  do i = ista_kngp0, iend_kngp      !for mpi
                     chgq_l(i,ri,is) = chgq_l(i,ri,is) + zfm3_l(i,it,ri)*f2*rhvg_l(i,it)
                  end do
               end do
               if ( ista_kngp == 1) then
                  f2 = ratio * ival2(it)/univol
                  chgq_l(1,1,is) = chgq_l(1,1,is) + zfm3_l(1,it,1)*f2
               endif
            end do
         end do
      endif
    end subroutine goto_from_PseudoPot_file

    subroutine goto_from_PseudoPot_file2
      integer :: is, it, ri, i, ia, t
      real(kind=DP) :: f1, f2, f3, ratio, c1, s1
      real(kind=DP) :: grt, nchg(4), ncharge

      if (ipri>=1) then
         write(nfout,*) 'Charge density initialization: from_PseudoPotential_FILE2 (atomic charge density)'
      endif
! --
      ista_kngp0 = ista_kngp
      if(ista_kngp0 == 1) ista_kngp0 = 2

      if ( noncol ) then
         if ( ndim_spinor /=2 ) then
            write(*,*) 'Not supported: ndim_spinor /=2'
            call phase_error_with_msg(nfout,'Not supported : ndim_spinor /=2 ',__LINE__,__FILE__)
         endif
         do ia=1, natm
            it = ityp(ia)

            ncharge = ival(it) +ionic_charge_atoms(ia)
            nchg(1) = ncharge
            nchg(2) = mag_moment0_atoms(ia,1)
            nchg(3) = mag_moment0_atoms(ia,2)
            nchg(4) = mag_moment0_atoms(ia,3)

            Do is=1, ndim_magmom
               f2 = nchg(is) /ival(it)
               do i = ista_kngp0, iend_kngp  !for mpi
                  grt = pos(ia,1)*ngabc(i,1) + pos(ia,2)*ngabc_kngp_l(i,2) &
                       &                     + pos(ia,3)*ngabc_kngp_l(i,3)
                  grt = grt *PAI2
                  if ( kimg == 1 ) then
                     c1 = iwei(ia) *cos(grt);
                     chgq_l(i,1,is) = chgq_l(i,1,is) + c1 *f2 *rhvg_l(i,it)
                  else
                     c1 = cos(grt);  s1 = -sin(grt)
                     chgq_l(i,1,is) = chgq_l(i,1,is) + c1 *f2 *rhvg_l(i,it)
                     chgq_l(i,2,is) = chgq_l(i,2,is) + s1 *f2 *rhvg_l(i,it)
                  endif
               end do

               f3 = nchg(is) /univol
               if ( ista_kngp == 1) then
                  c1 = iwei(ia);
                  chgq_l(1,1,is) = chgq_l(1,1,is) + c1 *f3
               endif
            end do
         end do

      else
         do ia = 1, natm
            it = ityp(ia)

            ncharge = ival(it) +ionic_charge_atoms(ia)
            if ( nspin == 1 ) then
               nchg(1) = ncharge
            else
               nchg(1) = ( ncharge +mag_moment0_atoms(ia,1) ) /2.0d0
               nchg(2) = ( ncharge -mag_moment0_atoms(ia,1) ) /2.0d0
            endif

            Do is=1, nspin
               f2 = nchg(is) /ival(it)
               do i = ista_kngp0, iend_kngp  !for mpi
                  grt = pos(ia,1)*ngabc_kngp_l(i,1) + pos(ia,2)*ngabc_kngp_l(i,2) &
                       &                            + pos(ia,3)*ngabc_kngp_l(i,3)
                  grt = grt *PAI2
                  if ( kimg == 1 ) then
                     c1 = iwei(ia) *cos(grt);
                     chgq_l(i,1,is) = chgq_l(i,1,is) + c1 *f2 *rhvg_l(i,it)
                  else
                     c1 = cos(grt);  s1 = -sin(grt)
                     chgq_l(i,1,is) = chgq_l(i,1,is) + c1 *f2 *rhvg_l(i,it)
                     chgq_l(i,2,is) = chgq_l(i,2,is) + s1 *f2 *rhvg_l(i,it)
                  endif
               end do

               f3 = nchg(is) /univol
               if ( ista_kngp == 1) then
                  c1 = iwei(ia);
                  chgq_l(1,1,is) = chgq_l(1,1,is) + c1 *f3
               endif
            end do
         end do
      endif
    end subroutine goto_from_PseudoPot_file2

    subroutine print_total_charge
      integer :: is, ni

      total_charge = 0.d0
      if(mype == 0) then
         if ( noncol ) then
            ni = 1
            total_charge = total_charge + chgq_l(1,1,ni)*univol
         else
            do is = 1, nspin, af+1
               total_charge = total_charge + chgq_l(1,1,is)*univol
            end do
         endif
      end if
      if ( npes > 1 ) call mpi_bcast( total_charge, 1, mpi_double_precision, 0, &
           &                          MPI_CommGroup, ierr )
      if ( ipritotalcharge >= 2 ) &
           & write(nfout,'(" !! total_charge = ",f15.6," <<m_CD_initial_CD_by_Gauss_kt>>")') total_charge
    end subroutine print_total_charge

    subroutine print_magnetic_moment
      integer :: is
      real(kind=DP) ::  magmom0(3)

      magmom0 = 0.d0
      if(mype == 0) then
         if ( noncol ) then
            Do is=2, ndim_magmom
               magmom0(is-1) = chgq_l(1,1,is)*univol
            End do
         else
            magmom0(3) = ( chgq_l(1,1,1)-chgq_l(1,1,2) ) *univol
         endif
      end if
!
      if ( npes > 1 ) call mpi_bcast( magmom0, 3, mpi_double_precision, 0, &
           &                          MPI_CommGroup, ierr )
      if ( ipritotalcharge >= 2 ) then
         write(nfout,'(" !! mag mom x    = ",f15.6," <<m_CD_initial_CD_by_Gauss_kt>>")') &
              &     magmom0(1)
         write(nfout,'(" !! mag mom y    = ",f15.6," <<m_CD_initial_CD_by_Gauss_kt>>")') &
              &     magmom0(2)
         write(nfout,'(" !! mag mom z    = ",f15.6," <<m_CD_initial_CD_by_Gauss_kt>>")') &
              &     magmom0(3)
      endif
    end subroutine print_magnetic_moment

  end subroutine m_CD_initial_CD_by_Gauss_kt
! ===============================================================================11.0

  subroutine m_CD_rd_chgq(nfout,nfchgt, F_CHGT_partitioned,prev)
    integer, intent(in) :: nfout, nfchgt
    logical, intent(in) :: F_CHGT_partitioned
    logical, intent(in), optional :: prev
    logical :: prv
    integer  :: i,j,k,is, ip
    real(kind=DP), allocatable, dimension(:,:,:) :: chgq_mpi

    integer :: id_sname = -1
    integer :: ierror
    integer :: npk,kgpt,istat,iendt
                                                 __TIMER_SUB_START(1374)
    call tstatc0_begin('m_CD_rd_chgq ',id_sname,1)
    prv = .false.
    if(present(prev)) prv = prev
    npk = np_kngp
    kgpt = kgp
    istat = ista_kngp
    iendt = iend_kngp
    if(prv) then
      npk = np_kngp_prev
      kgpt = kgp_prev
      istat = ista_kngp_prev
      iendt = iend_kngp_prev
      if(allocated(chgq_l_prev)) deallocate(chgq_l_prev)
      allocate(chgq_l_prev(ista_kngp_prev:iend_kngp_prev,kimg,ndim_magmom))
    endif

    if(F_CHGT_partitioned) then
       rewind nfchgt
       if(npes > 1) then
          allocate(chgq_mpi(npk,kimg,nspin))
! ============================================== by K. Tagami ============
        chgq_mpi = 0.0d0
! =======================================================================
                                                 __TIMER_IODO_START(1445)
          read(nfchgt) chgq_mpi
                                                 __TIMER_IODO_STOP(1445)
                                                 __TIMER_IODO_START(1446)
          if(prv) then
            do k = 1, nspin
               do j = 1, kimg
                  do i = 1, npk
                     ip = i + ista_kngp-1
                     chgq_l_prev(ip,j,k) = chgq_mpi(i,j,k)
                  end do
               end do
            end do
          else
            do k = 1, nspin
               do j = 1, kimg
                  do i = 1, npk
                     ip = i + ista_kngp-1
                     chgq_l(ip,j,k) = chgq_mpi(i,j,k)
                  end do
               end do
            end do
          endif
                                                 __TIMER_IODO_STOP(1446)
          deallocate(chgq_mpi)
       else
                                                 __TIMER_IODO_START(1447)
          if(prv)then
            read(nfchgt) chgq_l_prev
          else
            read(nfchgt) chgq_l
          endif
                                                 __TIMER_IODO_STOP(1447)
       end if
    else
       if(mype==0) rewind nfchgt
       if(npes > 1) then
          allocate(chgq_mpi1(kgpt,kimg,nspin)); chgq_mpi1 = 0.d0
          if(ipri >= 2) write(nfout,*) ' !D Reading chgq_g.'
                                                 __TIMER_IODO_START(1448)
          if(mype==0) read(nfchgt, end = 9999, err = 9999) chgq_mpi1
                                                 __TIMER_IODO_STOP(1448)
                                                 __TIMER_IOCOMM_START_w_BARRIER(MPI_CommGroup,1449)
          call mpi_bcast(chgq_mpi1,kgpt*kimg*nspin &
               & ,mpi_double_precision,0,MPI_CommGroup,ierr)
                                                 __TIMER_IOCOMM_STOP(1449)
                                                 __TIMER_IODO_START(1450)
          if(prv)then
             do k = 1, nspin
                do j = 1, kimg
                   do i = istat, iendt  !for mpi
                      chgq_l_prev(i,j,k) = chgq_mpi1(i,j,k)
                   enddo
                enddo
             enddo
          else
             do k = 1, nspin
                do j = 1, kimg
                   do i = istat, iendt  !for mpi
                      chgq_l(i,j,k) = chgq_mpi1(i,j,k)
                   enddo
                enddo
             enddo
          endif
                                                 __TIMER_IODO_STOP(1450)
          deallocate(chgq_mpi1)
       else
                                                 __TIMER_IODO_START(1451)
          if(prv)then
            read(nfchgt, end = 9999, err = 9999) chgq_l_prev
          else
            read(nfchgt, end = 9999, err = 9999) chgq_l
          endif
                                                 __TIMER_IODO_STOP(1451)
       end if
    end if

#ifdef FFTW3
    if(prv)then
      call m_CD_gen_chgq_from_prev_chgq(nfout)
    endif
#endif

    if(ipri >= 3) write(nfout,*) ' !D Reading chgq_g finished'
    total_charge = 0.d0
    if(myrank_chg == 0) then
       do is = 1, nspin, af+1
          total_charge = total_charge + chgq_l(1,1,is)*univol
       end do
    end if
    if(nrank_chg > 1) call mpi_bcast(total_charge,1,mpi_double_precision,0,mpi_chg_world,ierr)

    if(iprichargedensity >= 1) write(nfout,'("* total_charge = ",d20.8," <<m_CD_rd_chgq>>")') total_charge
    if(iprichargedensity >= 2) then
       write(nfout,'(" !D chgq_l <<m_CD_rd_chgq>>")')
       write(nfout,'(5f16.8)') (chgq_l(ista_kngp:min(ista_kngp+20,iend_kngp),1,1))
    end if
    call tstatc0_end(id_sname)
                                                 __TIMER_SUB_STOP(1374)
    return
9999 continue
    ierror = EOF_REACHED
    call phase_error_wo_filename(ierror, nfout, nfchgt, __LINE__, __FILE__)
  end subroutine m_CD_rd_chgq



  subroutine m_CD_wd_chgq(nfchgt,F_CHGT_partitioned)
    integer, intent(in) :: nfchgt
    logical, intent(in) :: F_CHGT_partitioned
    real(kind=DP), allocatable, dimension(:,:,:) :: chgq_mpi
    integer :: i, ip, ri, is
    integer :: id_sname = -1
                                                 __TIMER_SUB_START(1375)
    call tstatc0_begin('m_CD_wd_chgq ',id_sname,1)

    if(F_CHGT_partitioned) then
       rewind nfchgt
       if(npes > 1) then
          allocate(chgq_mpi(np_kngp,kimg,nspin))
! ============================================ by K. Tagami =======
	chgq_mpi = 0
! =================================================================
                                                 __TIMER_IODO_START(1452)
          do is = 1, nspin
             do ri = 1, kimg
                do i = ista_kngp, iend_kngp
                   ip = i - ista_kngp+1
                   chgq_mpi(ip,ri,is) = chgq_l(i,ri,is)
                end do
             end do
          end do
                                                 __TIMER_IODO_STOP(1452)
                                                 __TIMER_IODO_START(1453)
          write(nfchgt) chgq_mpi
                                                 __TIMER_IODO_STOP(1453)
          deallocate(chgq_mpi)
       else
                                                 __TIMER_IODO_START(1454)
          write(nfchgt) chgq_l
                                                 __TIMER_IODO_STOP(1454)
       end if
    else
       if(mype==0) rewind nfchgt
       if(npes > 1) then
          allocate(chgq_mpi1(kgp,kimg,nspin)); chgq_mpi1 = 0.d0
          allocate(chgq_mpi2(kgp,kimg,nspin)); chgq_mpi2 = 0.d0

          chgq_mpi1(ista_kngp:iend_kngp,:,:) = chgq_l(ista_kngp:iend_kngp,:,:)
                                                 __TIMER_IOCOMM_START_w_BARRIER(mpi_chg_world,1455)
          call mpi_allreduce(chgq_mpi1,chgq_mpi2,kgp*kimg*nspin &
               &  ,mpi_double_precision,mpi_sum,mpi_chg_world,ierr)
                                                 __TIMER_IOCOMM_STOP(1455)
                                                 __TIMER_IODO_START(1456)
          if(mype==0) write(nfchgt) chgq_mpi2
                                                 __TIMER_IODO_STOP(1456)
          deallocate(chgq_mpi1); deallocate(chgq_mpi2)
       else
                                                 __TIMER_IODO_START(1457)
          write(nfchgt) chgq_l
                                                 __TIMER_IODO_STOP(1457)
       end if
    end if
    call tstatc0_end(id_sname)
                                                 __TIMER_SUB_STOP(1375)
  end subroutine m_CD_wd_chgq


  subroutine charge_average_3D(mode,chg)

    integer,intent(in)           :: mode
    real(kind=DP), intent(inout) :: chg(ista_kngp:iend_kngp,kimg,nspin)
    integer ::       ispin, ng, no, ngp, no1, no2
    real(kind=DP) :: fi, tx,ty,tz, fp, fc, fs, zcr, zci
    real(kind=DP), pointer, dimension(:,:) :: work2
    logical, save                           :: firstcall = .true.
    integer, save                           :: sendmax, recvmax, sendranks, recvranks
    integer                                 :: mymin, mymax, lrank
    integer,                    dimension(2)   :: my_from_to
    integer,       allocatable, dimension(:,:) :: all_from_to
    integer,       allocatable, dimension(:)   :: is_ngpt, ie_ngpt
    integer, save, allocatable, dimension(:,:) :: sendinfo, recvinfo
                                                 __TIMER_SUB_START(731)
    allocate(work2(ista_kngp:iend_kngp,kimg)); work2 = 0.d0

    if(mode == ANTIFERRO) then
       fi = 1.d0/af
       no1 = nopr + 1; no2 = nopr + af
    else
       fi = 1.d0/nopr
       no1 = 1; no2 = nopr
    end if

!XX!if (firstcall) then
    allocate(all_from_to(2,0:nrank_chg-1))
    allocate(is_ngpt(0:nrank_chg-1))
    allocate(ie_ngpt(0:nrank_chg-1))
    allocate(sendinfo(2,0:nrank_chg-1))
    allocate(recvinfo(2,0:nrank_chg-1))
    sendinfo = 0
    recvinfo = 0
    sendmax = 0
    recvmax = 0
    sendranks = 0
    recvranks = 0

    mymin = kgp
    mymax = 1
                                                 __TIMER_DO_START(848)
    do no = no1, no2
       do ng = ista_kngp, iend_kngp !for mp
          ngp = ngpt_l(ng,no)
          if (mymin > ngp) mymin = ngp
          if (mymax < ngp) mymax = ngp
       end do
    end do
                                                 __TIMER_DO_STOP(848)

    my_from_to(1) = mymin
    my_from_to(2) = mymax

                                                 __TIMER_COMM_START_w_BARRIER(mpi_chg_world,849)
    call mpi_allgather(my_from_to , 2, MPI_INTEGER, &
   &                   all_from_to, 2, MPI_INTEGER, mpi_chg_world, ierr)
                                                 __TIMER_COMM_STOP(849)
    is_ngpt(0:nrank_chg-1) = all_from_to(1,0:nrank_chg-1)
    ie_ngpt(0:nrank_chg-1) = all_from_to(2,0:nrank_chg-1)

                                                 __TIMER_DO_START(850)
    do lrank = 0, nrank_chg-1
      if(lrank == myrank_chg) cycle
      if(ie_ngpt(lrank) < ista_kngp) cycle
      if(iend_kngp < is_ngpt(lrank)) cycle
      if(ista_kngp <= ie_ngpt(lrank)) then
         sendinfo(1,lrank) = max(ista_kngp,is_ngpt(lrank))
         sendinfo(2,lrank) = min(iend_kngp,ie_ngpt(lrank))
         sendranks = sendranks + 1
      else if(is_ngpt(lrank) <= iend_kngp) then
         sendinfo(1,lrank) = max(ista_kngp,is_ngpt(lrank))
         sendinfo(2,lrank) = min(iend_kngp,ie_ngpt(lrank))
         sendranks = sendranks + 1
      endif
      sendmax = max(sendmax,sendinfo(2,lrank)-sendinfo(1,lrank)+1)
    end do
                                                 __TIMER_DO_STOP(850)
                                                 __TIMER_DO_START(851)
    do lrank = 0, nrank_chg-1
      if(lrank == myrank_chg) cycle
      if(ie_kngp(lrank) < mymin) cycle
      if(mymax < is_kngp(lrank)) cycle
      if(mymin <= ie_kngp(lrank)) then
         recvinfo(1,lrank) = max(mymin,is_kngp(lrank))
         recvinfo(2,lrank) = min(mymax,ie_kngp(lrank))
         recvranks = recvranks + 1
      else if(is_kngp(lrank) <= mymax) then
         recvinfo(1,lrank) = max(mymin,is_kngp(lrank))
         recvinfo(2,lrank) = min(mymax,ie_kngp(lrank))
         recvranks = recvranks + 1
      endif
      recvmax = max(recvmax,recvinfo(2,lrank)-recvinfo(1,lrank)+1)
    end do
                                                 __TIMER_DO_STOP(851)
    deallocate(all_from_to)
    deallocate(is_ngpt)
    deallocate(ie_ngpt)

!XX!firstcall = .false.
!XX!end if

    allocate(work(mymin:mymax,kimg))

    do ispin = 1, nspin, af+1
       call cp_chgq_by_ngpt() ! chg -> work
       work2 = 0.d0                   ! initialization
                                                 __TIMER_DO_START(852)
       do no = no1, no2
!!$          tx = tau(1,no,BUCS)*PAI2
!!$          ty = tau(2,no,BUCS)*PAI2
!!$          tz = tau(3,no,BUCS)*PAI2
          if(kimg == 1) then
             do ng = ista_kngp, iend_kngp !for mpi
                ngp = ngpt_l(ng,no)
!                fp = ngabc(ngp,1)*tx + ngabc(ngp,2)*ty + ngabc(ngp,3)*tz
!                fp = ngabc_kngp_l(ngp,1)*tx + ngabc_kngp_l(ngp,2)*ty + ngabc_kngp_l(ngp,3)*tz
                fp = fp_l(ng,no)
                work2(ng,1)        = work2(ng,1) + dcos(fp)*work(ngp,1)
             end do
          else if(kimg == 2) then
             do ng = ista_kngp, iend_kngp !for mpi
                ngp= ngpt_l(ng,no)
!                fp = ngabc(ngp,1)*tx + ngabc(ngp,2)*ty + ngabc(ngp,3)*tz
!                fp = ngabc_kngp_l(ngp,1)*tx + ngabc_kngp_l(ngp,2)*ty + ngabc_kngp_l(ngp,3)*tz
                fp = fp_l(ng,no)
                fc = dcos(fp);     fs = dsin(fp)
                zcr= work(ngp,1);  zci= work(ngp,kimg)
                work2(ng,1)        = work2(ng,1) + fc*zcr - fs*zci
                work2(ng,2)        = work2(ng,2) + fc*zci + fs*zcr
             end do
          end if
       end do
                                                 __TIMER_DO_STOP(852)
       if(mode /= ANTIFERRO) chg(:,:,ispin) = work2(:,:)*fi
    end do

    if(mode == ANTIFERRO) chg(:,:,nspin) = work2(:,:)*fi

    deallocate(work2)
    deallocate(sendinfo)
    deallocate(recvinfo)
    deallocate(work)
                                                 __TIMER_SUB_STOP(731)
   contains

    subroutine cp_chgq_by_ngpt()
      real(DP), allocatable, dimension(:,:) :: sendbuf
      real(DP), allocatable, dimension(:,:) :: recvbuf
      integer, dimension(sendranks) :: req_s
      integer, dimension(recvranks) :: req_r, src
      integer, dimension(MPI_STATUS_SIZE,sendranks) :: sta_s
      integer, dimension(MPI_STATUS_SIZE,recvranks) :: sta_r

      integer :: i, j, ri, ista, iend, nel
      integer :: icnt_s, icnt_r, ierr, itag=30
! === DEBUG by tkato 2012/06/04 ================================================
#ifndef USE_NONBLK_COMM
#ifndef USE_ALLTOALLV
      real(kind=DP), allocatable, dimension(:,:) :: sbuf, rbuf
      integer :: maxbuf
#else
      integer, allocatable, dimension(:) :: sdsp, rdsp
      integer, allocatable, dimension(:) :: scnt, rcnt
      real(kind=DP), allocatable, dimension(:,:) :: sbuf, rbuf
#endif
#endif
! ==============================================================================

                                                 __TIMER_SUB_START(732)
      allocate(recvbuf(recvmax*kimg,recvranks))
      allocate(sendbuf(sendmax*kimg,sendranks))

      icnt_r = 0
#ifdef USE_NONBLK_COMM
                                                 __TIMER_COMM_START_w_BARRIER(mpi_chg_world,853)
#endif
      do lrank = 0, nrank_chg - 1
         if(recvinfo(1,lrank) == 0) cycle
         ista = recvinfo(1,lrank)
         iend = recvinfo(2,lrank)
         nel  = iend-ista+1
         icnt_r = icnt_r + 1
         src(icnt_r) = lrank
#ifdef USE_NONBLK_COMM
         call mpi_irecv(recvbuf(1,icnt_r), nel*kimg, mpi_double_precision, &
        &               lrank, itag, mpi_chg_world, req_r(icnt_r), ierr)
          if (ierr /= 0) then
             call mpi_abort(mpi_comm_world, 171, ierr)
          endif
#endif
      end do

      icnt_s = 0
      do lrank = 0, nrank_chg-1
         if(sendinfo(1,lrank) == 0) cycle
         ista = sendinfo(1,lrank)
         iend = sendinfo(2,lrank)
         nel  = iend-ista+1
         icnt_s = icnt_s + 1
                                                 __TIMER_DO_START(854)
         do ri = 1, kimg
            do i = ista, iend
               sendbuf((i-ista+1+nel*(ri-1)),icnt_s) = chg(i,ri,ispin)
            end do
         end do
                                                 __TIMER_DO_STOP(854)
#ifdef USE_NONBLK_COMM
         call mpi_isend(sendbuf(1,icnt_s), nel*kimg, mpi_double_precision, &
        &               lrank, itag, mpi_chg_world, req_s(icnt_s), ierr)
          if (ierr /= 0) then
             call mpi_abort(mpi_comm_world, 172, ierr)
          endif
#endif
      end do
                                                 __TIMER_DO_START(855)
      do ri = 1, kimg
         do i = ista_kngp, iend_kngp
            if (i < mymin) cycle
            if (mymax < i) cycle
            work(i,ri) = chg(i,ri,ispin)
         end do
      end do
                                                 __TIMER_DO_STOP(855)
#ifdef USE_NONBLK_COMM
      call mpi_waitall(icnt_s, req_s, sta_s, ierr)
       if (ierr /= 0) then
          call mpi_abort(mpi_comm_world, 173, ierr)
       endif
      call mpi_waitall(icnt_r, req_r, sta_r, ierr)
       if (ierr /= 0) then
          call mpi_abort(mpi_comm_world, 174, ierr)
       endif
                                                 __TIMER_COMM_STOP(853)
#else
                                                 __TIMER_COMM_START_w_BARRIER(mpi_chg_world,694)
#ifndef USE_ALLTOALLV
! === DEBUG by tkato 2012/06/05 ================================================
!      if(sendmax/=0)then
! ==============================================================================
! === DEBUG by tkato 2012/06/04 ================================================
!      real(kind=DP), allocatable, dimension(:,:) :: sbuf, rbuf
! ==============================================================================
! === DEBUG by tkato 2012/06/04 ================================================
!      allocate(sbuf(sendmax*kimg,0:sendranks-1), stat=ierr)
!      allocate(rbuf(recvmax*kimg,0:recvranks-1), stat=ierr)
!      do i = 0, nrank_g - 1
!         sbuf(:,i)=sendbuf(:,i+1)
!      enddo
!      call MPI_ALLTOALL( sbuf, nel*kimg, mpi_double_precision, &
!     &                   rbuf, nel*kimg, mpi_double_precision, &
!     &                                    mpi_ke_world, ierr )
!      if (ierr /= 0) then
!         call mpi_abort(mpi_comm_world, 175, ierr)
!      endif
!      do i = 0, nrank_g - 1
!         recvbuf(:,i+1)=rbuf(:,i)
!      enddo
       maxbuf = max(sendmax,recvmax)
       call mpi_allreduce(MPI_IN_PLACE,maxbuf,1,mpi_integer,mpi_max,mpi_chg_world,ierr)
       allocate(sbuf(maxbuf*kimg,0:nrank_chg-1), stat=ierr)
       allocate(rbuf(maxbuf*kimg,0:nrank_chg-1), stat=ierr)
       icnt_s = 0
       do i = 0, nrank_chg - 1
          if(sendinfo(1,i) == 0) cycle
          icnt_s = icnt_s + 1
          sbuf(1:sendmax*kimg,i)=sendbuf(1:sendmax*kimg,icnt_s)
       enddo
       call MPI_ALLTOALL( sbuf, maxbuf*kimg, mpi_double_precision, &
      &                   rbuf, maxbuf*kimg, mpi_double_precision, &
      &                                    mpi_chg_world, ierr )
       if (ierr /= 0) then
          call mpi_abort(mpi_comm_world, 175, ierr)
       endif
       icnt_r = 0
       do i = 0, nrank_chg - 1
          if(recvinfo(1,i) == 0) cycle
          icnt_r = icnt_r + 1
          recvbuf(1:recvmax*kimg,icnt_r)=rbuf(1:recvmax*kimg,i)
       enddo
! ==============================================================================
       deallocate(sbuf)
       deallocate(rbuf)
! === DEBUG by tkato 2012/06/05 ================================================
!      endif
! ==============================================================================
#else
! === DEBUG by tkato 2012/06/05 ================================================
!      if(sendmax/=0)then
! ==============================================================================
! === DEBUG by tkato 2012/06/05 ================================================
!      integer, allocatable, dimension(:) :: sdsp, rdsp
!      integer, allocatable, dimension(:) :: scnt, rcnt
!      real(kind=DP), allocatable, dimension(:,:) :: sbuf, rbuf
!      allocate(sdsp(0:sendranks-1), stat=ierr)
!      allocate(rdsp(0:recvranks-1), stat=ierr)
!      allocate(scnt(0:sendranks-1), stat=ierr)
!      allocate(rcnt(0:recvranks-1), stat=ierr)
!      allocate(sbuf(sendmax*kimg,0:sendranks-1), stat=ierr)
!      allocate(rbuf(recvmax*kimg,0:recvranks-1), stat=ierr)
!      write(6,*) sendmax, recvmax, kimg
!      do i = 0, nrank_g - 1
!         sdsp(i)=sendmax*kimg*i
!         rdsp(i)=recvmax*kimg*i
!         scnt(i)=(sendinfo(2,i)-sendinfo(1,i)+1)*kimg
!         rcnt(i)=(recvinfo(2,i)-recvinfo(1,i)+1)*kimg
!         sbuf(:,i)=sendbuf(:,i+1)
!      enddo
!      call MPI_ALLTOALLV(      sendb, scnt, sdsp, &
!     &   mpi_double_precision, recvb, rcnt, rdsp, &
!     &   mpi_double_precision, mpi_ke_world, ierr )
!      if (ierr /= 0) then
!         call mpi_abort(mpi_comm_world, 175, ierr)
!      endif
!      do i = 0, nrank_g - 1
!         recvbuf(:,i+1)=rbuf(:,i)
!      enddo
       allocate(sdsp(0:nrank_chg-1), stat=ierr)
       allocate(rdsp(0:nrank_chg-1), stat=ierr)
       allocate(scnt(0:nrank_chg-1), stat=ierr)
       allocate(rcnt(0:nrank_chg-1), stat=ierr)
! === DEBUG by tkato 2012/06/05 ================================================
!      allocate(sbuf(sendmax*kimg,0:nrank_g-1), stat=ierr)
!      allocate(rbuf(recvmax*kimg,0:nrank_g-1), stat=ierr)
       if(sendmax /= 0) then
          allocate(sbuf(sendmax*kimg,0:nrank_chg-1), stat=ierr)
       else
          allocate(sbuf(1,0:nrank_chg-1), stat=ierr)
       endif
       if(recvmax /= 0) then
          allocate(rbuf(recvmax*kimg,0:nrank_chg-1), stat=ierr)
       else
          allocate(rbuf(1,0:nrank_chg-1), stat=ierr)
       endif
! ==============================================================================
       write(6,*) sendmax, recvmax, kimg
       icnt_s = 0
       do i = 0, nrank_chg - 1
          sdsp(i)=sendmax*kimg*i
          rdsp(i)=recvmax*kimg*i
          scnt(i)=(sendinfo(2,i)-sendinfo(1,i)+1)*kimg
          rcnt(i)=(recvinfo(2,i)-recvinfo(1,i)+1)*kimg
          if(sendinfo(1,i) == 0) cycle
          icnt_s = icnt_s + 1
          sbuf(:,i)=sendbuf(:,icnt_s)
       enddo
       call MPI_ALLTOALLV(      sbuf, scnt, sdsp, &
      &   mpi_double_precision, rbuf, rcnt, rdsp, &
      &   mpi_double_precision, mpi_chg_world, ierr )
       if (ierr /= 0) then
          call mpi_abort(mpi_comm_world, 175, ierr)
       endif
       icnt_r = 0
       do i = 0, nrank_chg - 1
          if(recvinfo(1,i) == 0) cycle
          icnt_r = icnt_r + 1
          recvbuf(:,icnt_r)=rbuf(:,i)
       enddo
! ==============================================================================
       deallocate(sdsp)
       deallocate(rdsp)
       deallocate(scnt)
       deallocate(rcnt)
       deallocate(sbuf)
       deallocate(rbuf)
! === DEBUG by tkato 2012/06/05 ================================================
!      endif
! ==============================================================================
#endif
                                                 __TIMER_COMM_STOP(694)
#endif
                                                 __TIMER_DO_START(856)
      do j = 1, icnt_r
         lrank = src(j)
         ista = recvinfo(1,lrank)
         iend = recvinfo(2,lrank)
         nel  = iend-ista+1
         do ri = 1, kimg
            do i = ista, iend
               work(i,ri) = recvbuf((i-ista+1+nel*(ri-1)),j)
            end do
         end do
      end do
                                                 __TIMER_DO_STOP(856)
      deallocate(sendbuf)
      deallocate(recvbuf)
                                                 __TIMER_SUB_STOP(732)
    end subroutine cp_chgq_by_ngpt

  end subroutine charge_average_3D


  subroutine cp_chg_to_work(ispin,chg)
    integer, intent(in) :: ispin
    real(DP),intent(in),dimension(ista_kngp:iend_kngp,kimg,nspin) :: chg
    integer :: ng,ri
    real(kind=DP), allocatable, dimension(:,:) :: work_mpi
    work = 0.d0
    do ri = 1, kimg
       do ng = ista_kngp, iend_kngp  !for mpi
!!$          work_mpi(ng,ri) = chg(ng,ri,ispin)
          work(ng,ri) = chg(ng,ri,ispin)
       end do
    end do
    if(npes >= 2) then
       allocate(work_mpi(kgp,kimg)); work_mpi = 0.d0
       call mpi_allreduce(work,work_mpi,kgp*kimg &
            &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       work = work_mpi
       deallocate(work_mpi)
    end if
  end subroutine cp_chg_to_work

! ==================================== added by K. Tagami ============== 11.0
  subroutine cp_chg_to_work_noncl(ispin,chg)
    integer, intent(in) :: ispin
    real(DP),intent(in),dimension(ista_kngp:iend_kngp,kimg,ndim_magmom) :: chg
    integer :: ng,ri
    real(kind=DP), allocatable, dimension(:,:) :: work_mpi
    work = 0.d0
    do ri = 1, kimg
       do ng = ista_kngp, iend_kngp  !for mpi
!!$          work_mpi(ng,ri) = chg(ng,ri,ispin)
          work(ng,ri) = chg(ng,ri,ispin)
       end do
    end do
    if(npes >= 2) then
       allocate(work_mpi(kgp,kimg)); work_mpi = 0.d0
       call mpi_allreduce(work,work_mpi,kgp*kimg &
            &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       work = work_mpi
       deallocate(work_mpi)
    end if
  end subroutine cp_chg_to_work_noncl
! ======================================================================== 11.0
!!$  subroutine initialize_chgq_l(ispin)
!!$    integer, intent(in) :: ispin
!!$    chgq_l(:,:,ispin) = 0.d0
!!$  end subroutine initialize_chgq_l

!!$  subroutine devide_v_with_vdF(v_l,v_dF)
!!$    real(DP),intent(inout),dimension(ista_kgpm:iend_kgpm,kimg,nspin_m):: v_l
!!$    real(DP),intent(in),   dimension(nspin)                           :: v_dF
!!$    integer :: is
!!$    do is = 1, nspin, af+1
!!$       v_l(:,:,is) = v_l(:,:,is)/v_dF(is)
!!$    end do
!!$  end subroutine devide_v_with_vdF

  subroutine m_CD_alloc_rspace_charge()
    allocate(afft(nfftp_nonpara)); afft = 0.d0
    allocate(afft_mpi1(nfftp_nonpara))
! ======================================== Added by K. Tagami =======
    afft = 0.0d0; afft_mpi1 = 0.0d0
! ==================================================================
    call m_FFT_alloc_CD_box()
  end subroutine m_CD_alloc_rspace_charge

  subroutine m_CD_dealloc_rspace_charge()
    deallocate(afft); deallocate(afft_mpi1)
    call m_FFT_dealloc_CD_box()
  end subroutine m_CD_dealloc_rspace_charge

  subroutine m_CD_initial_CD_by_file_rspace(nspin, iloop,nfout,nfchr)
! Coded by T. Yamasaki, 28 July 2008
    integer, intent(in) :: nspin, iloop, nfout, nfchr
    integer :: is

    call m_CD_alloc_rspace_charge()
    if(mype == 0) call rdchgr(nfout,nfchr) ! -> afft
    if(npes >= 2) then
       call mpi_allreduce(afft,afft_mpi1,nfftp_nonpara,mpi_double_precision &
            & , mpi_sum, MPI_CommGroup,ierr)
       afft = afft_mpi1
       if(mype == 0) write(nfout,'(" after <<mpi_allreduce>>")')
    end if
    call m_FFT_CD0(nfout,afft,DIRECT)
    if(iprichargedensity >= 1) write(nfout,'(" after m_FFT_CD0")')
    call cpafft_CD_to_valencecharge(iloop) ! afft -> chgq_l
    if(iprichargedensity >= 1) write(nfout,'(" after cpafft_CD_to_valencecharge")')

    call m_CD_dealloc_rspace_charge()
  contains

    subroutine rdchgr(nfout,nfchr)
      integer, intent(in) :: nfout, nfchr
      real(kind=DP),allocatable,dimension(:,:,:) :: wkchr
      integer :: idp,mmp,nlp,nmp,nnp,nlp_t,nmp_t,nnp_t, nfftp_t,inew,jnew,knew,nlphf
      integer :: i,j,k,ip,natm2_t
      real(kind=DP) :: totch_from_fft

      idp = fft_box_size_CD_nonpara(1,0)
      mmp = fft_box_size_CD_nonpara(2,0)
      nlp = fft_box_size_CD(1,1)
      nmp = fft_box_size_CD(2,1)
      nnp = fft_box_size_CD(3,1)

      nlphf = idp/(3-kimg)  ! = idp/2 (if kimg == 1), or idp (if kimg == 2)

      if(initial_charge_filetype == DENSITY_ONLY) then
         allocate(wkchr(nlp,nmp,nnp)); wkchr = 0.d0
      else if(charge_filetype == CUBE) then
         allocate(wkchr(nnp,nmp,nlp)); wkchr = 0.d0
      end if

      if(ipri >= 2) write(nfout,'(" !D Charge density ne = ",i8 &
           & ,"(",3i5," ) <<m_CD_initial_CD_by_file_rspace>>")') nlp*nmp*nnp, nlp, nmp, nnp
      if(ipri >= 2) write(nfout,'(" !D charge density is being read")')
      if(initial_charge_filetype == DENSITY_ONLY) then
         read(nfchr,'(21a,i8,1a,3i5)') nfftp_t, nlp_t, nmp_t, nnp_t
         call check_fftsize(nlp_t,nmp_t,nnp_t,nlp,nmp,nnp)
         read(nfchr,'(6e13.5)') wkchr
         if(iprichargedensity >= 2) &
              & write(nfout,'(" !CD initial_charge_filetype = DENSITY_ONLY")')
      else if(initial_charge_filetype == CUBE) then
         read(nfchr,*) ! comment
         read(nfchr,*) ! a comment of "SCF Total Density"
         read(nfchr,'(i6)') natm2_t ! natm2,x,y,z
         read(nfchr,'(i6)') nlp_t
         read(nfchr,'(i6)') nmp_t
         read(nfchr,'(i6)') nnp_t
         nfftp_t =  nlp_t*nmp_t*nnp_t
         call check_fftsize(nlp_t,nmp_t,nnp_t,nlp,nmp,nnp)
         do i = 1, natm2
            read(nfchr,*)
         end do
! --> T. Yamasaki, 1 Aug. 2008
         read(nfchr,*) wkchr
!!$         read(nfchr,'(6e13.5)') wkchr
!!$         do i = 1, nlp
!!$            do j = 1, nmp
!!$               read(nfchr,*) (wkchr(k,j,i),k=1,nnp)
!!$            end do
!!$         end do
         if(iprichargedensity>=2) then
            write(nfout,'(" -- wkchr --")')
            write(nfout,'(" wkchr(1:6,nmp,nlp)= ",6f13.5)') wkchr(1:min(6,nnp),nmp,nlp)
         end if
! <--
         if(iprichargedensity >= 1) &
              & write(nfout,'(" !CD initial_charge_filetype = CUBE")')
      else
         if(ipri >= 1) write(nfout,'(" initila_charge_filetype is invalid")')
         call phase_error_with_msg(nfout,' initial_charge_filetype is invalid << m_CD_initial_CD_by_file_rspace>>'&
                                  ,__LINE__,__FILE__)
      end if
!!$      if(wkchr(1,1,1) > DELTA10) then
!!$         if(ipri >= 1) write(nfout,'(" The charge density is analyzed as symmetric")')
!!$         up_down = UP
!!$      else


      if(iprichargedensity >= 2) then
         write(nfout,'(" !CD nfftp, nlp, nmp, nnp = ",4i8)') nfftp_t,nlp_t,nmp_t,nnp_t
         write(nfout,'(" !CD wkchr <<m_CD_initial_CD_file_rspace>>")')
         if(initial_charge_filetype == DENSITY_ONLY) then
            write(nfout,'(5f16.8)') wkchr(1:min(nlp,20),1,1)
         else if(initial_charge_filetype == CUBE) then
            write(nfout,'(5f16.8)') wkchr(1:min(nnp,20),1,1)
         end if
         write(nfout,'(" nfftp, nlp, nmp, nnp, nlphf = ",i10,4i8)') nfftp,nlp,nmp,nnp,nlphf
      end if

      afft = 0.d0
      do i = 1, nmp
         do j = 1, nnp
            do k = 1, nlp
               if(kimg == 1 .and. k > nlphf) then
                  knew = idp - k
                  jnew = nnp+2 - j
                  inew = nmp+2 - i
                  if(jnew > nnp) then
                     jnew = jnew - nnp
                  end if
                  if(inew > nmp) then
                     inew = inew - nmp
                  end if
               else
                  knew = k; jnew = j; inew = i
               end if
               ip = nlphf*mmp*(jnew-1) + nlphf*(inew-1) + knew
!!$               if(ip*2 > nfftp) then
!!$                  write(nfout,'(" !CD <<rdchgr>> ip*2=",i10," > nfftp")') ip*2
!!$                  stop ' ip*2 > nfftp'
!!$               end if
               if(charge_filetype == DENSITY_ONLY) then
                  afft(ip*2-1) = wkchr(k,i,j)
               else if(charge_filetype == CUBE) then
                  afft(ip*2-1) = wkchr(j,i,k)
               end if
            end do
         end do
      end do

      if(iprichargedensity >= 1 .and. (charge_filetype == DENSITY_ONLY .or. charge_filetype == CUBE)) then
         totch_from_fft = 0.d0
         if(charge_filetype == DENSITY_ONLY) then
            do j = 1, nnp
               do i = 1, nmp
                  do k = 1, nlp
                     totch_from_fft = totch_from_fft + wkchr(k,i,j)
                  end do
               end do
            end do
         else if(charge_filetype == CUBE) then
            do k = 1, nlp
               do i = 1, nmp
                  do j = 1, nnp
                     totch_from_fft = totch_from_fft + wkchr(j,i,k)
                  end do
               end do
            end do
         end if
         if(nnp*nmp*nlp /= 0) then
            totch_from_fft = totch_from_fft*univol/(nnp*nmp*nlp)
            write(nfout,'(" totch_from_fft = ",d20.8,"<<rdchgr>>")') totch_from_fft
         else
            write(nfout,'(" nnp*nmp*nlp = ",i8,"<<rdchgr>>")') nnp*nmp*nlp
         end if
      end if

      if(allocated(wkchr)) deallocate(wkchr)
      if(iprichargedensity >= 2) write(nfout,'(" end of <<rdchgr>>")')
    end subroutine rdchgr

    subroutine check_fftsize(nlp_t,nmp_t,nnp_t,nlp,nmp,nnp)
      integer, intent(in) :: nlp_t, nmp_t, nnp_t, nlp, nmp, nnp
      character(34) :: chsubroutinename = "<<m_CD_initial_CD_by_file_rspace>>"
      if(ipri >= 1) &
           & write(nfout,'(" nlp_t, nmp_t, nnp_t, nlp, nmp, nnp = ",6i6)') nlp_t,nmp_t,nnp_t,nlp,nmp,nnp
      if(nlp_t /= nlp) then
         if(ipri >=1 ) write(nfout,'(" nlp_t /= nlp ", a34)') chsubroutinename
      end if
      if(nmp_t /= nmp) then
         if(ipri >=1 ) write(nfout,'(" nmp_t /= nmp ", a34)') chsubroutinename
      end if
      if(nnp_t /= nnp) then
         if(ipri >=1 ) write(nfout,'(" nnp_t /= nnp ", a34)') chsubroutinename
      end if
      if(nlp_t /= nlp .or. nmp_t /= nmp .or. nnp_t /= nnp) then
         call phase_error_with_msg(nfout,&
       ' a set of (nlp,nmp,nnp) in the nfchr-file is invalid <<m_CD_initial_CD_by_file_rspace>>'&
       ,__LINE__,__FILE__)
      end if
    end subroutine check_fftsize

  end subroutine m_CD_initial_CD_by_file_rspace

    subroutine cpafft_CD_to_valencecharge(iloop)
      use m_Files, only : nfout
      integer, intent(in) :: iloop
      integer :: i,j, ip, is
      real(kind=DP) :: rinplw

#ifdef _MPIFFTTEST_
      rinplw = 1.d0/product(fft_box_size_CD_c(1:3,1))
#else
      rinplw = 1.d0/product(fft_box_size_CD(1:3,1))
#endif

      chgq_l(:,:,iloop) = 0.d0
      do  j = 1, kimg
         do i = ista_kngp, iend_kngp
            if( kgp_reduced < i ) cycle
#ifdef _MPIFFT_
            ip = (igfp_nonpara(i)-1)*kimg+j
#elif _MPIFFTTEST_
            ip = (igfp_l_c(i)-1)*kimg + j
#else
            ip = (igfp_l(i)-1)*kimg + j
#endif
            chgq_l(i,j,iloop) = afft(ip)*rinplw
         end do
      end do
      total_charge = 0.d0
      if(mype == 0) total_charge = chgq_l(1,1,iloop)*univol
      if(npes > 1) call mpi_bcast(total_charge,1,mpi_double_precision,0,MPI_CommGroup,ierr)

!!$      if(dabs(total_charge-totch)>1.d-5) then
!!$         chgq_l(:,:,iloop) = chgq_l(:,:,iloop)*totch/total_charge
!!$      end if

      if(iprichargedensity >= 1) then
         write(nfout,'("* total_charge(",i3,") = ",d20.8," <<m_CD_initial_CD_by_file_rspace>>")') iloop,total_charge
      end if

      if(iprichargedensity >= 1) write(nfout,'("* total_charge = ",d20.8," <<m_CD_initial_CD_by_file_rspace>>")') total_charge
      if(iprichargedensity >= 2) then
         write(nfout,'(" !D chgq_l <<m_CD_initial_CD_by_file_rspace>>")')
         write(nfout,'(5f16.8)') (chgq_l(ista_kngp:min(ista_kngp+20,iend_kngp),1,1))
      end if

    end subroutine cpafft_CD_to_valencecharge
! ============================== added by K. Tagami ======================== 11.0
  subroutine m_CD_initCD_by_file_rsp_noncl( iloop,nfout,nfchr )
! Coded by T. Yamasaki, 28 July 2008
    integer, intent(in) :: iloop, nfout, nfchr
    integer :: is

    call m_CD_alloc_rspace_charge()         ! prepare afft
    if(mype == 0) call rdchgr(nfout,nfchr) ! -> afft
    if(npes >= 2) then
       call mpi_allreduce(afft,afft_mpi1,nfftp_nonpara,mpi_double_precision &
            & , mpi_sum, MPI_CommGroup,ierr)
       afft = afft_mpi1
       if(mype == 0) write(nfout,'(" after <<mpi_allreduce>>")')
    end if
    call m_FFT_CD0(nfout,afft,DIRECT)
    if(iprichargedensity >= 1) write(nfout,'(" after m_FFT_CD0")')
    call cpafft_CD_to_valencecharge(iloop) ! afft -> chgq_l
    if(iprichargedensity >= 1) write(nfout,'(" after cpafft_CD_to_valencecharge")')

    if ( import_collinear_spindensity == ON ) then
       if ( iloop == ndim_magmom ) then
          call redistrib_charge_density
       endif
    endif

  contains

    subroutine redistrib_charge_density
      integer :: i,j
      real(kind=DP) :: c1, c2

      if ( sw_fix_global_quantz_axis == ON ) then
         do j = 1, kimg
            do i = ista_kngp, iend_kngp
               c1 = chgq_l( i,j,1 ) + chgq_l( i,j,ndim_magmom )
               c2 = chgq_l( i,j,1 ) - chgq_l( i,j,ndim_magmom )
               chgq_l( i,j,1 ) = c1
               chgq_l( i,j,2:4 ) = c2 *Global_Quantz_Axis_Fixed(1:3)
            end do
         end do
      else
         do j = 1, kimg
            do i = ista_kngp, iend_kngp
               c1 = chgq_l( i,j,1 ) + chgq_l( i,j,ndim_magmom )
               c2 = chgq_l( i,j,1 ) - chgq_l( i,j,ndim_magmom )
               chgq_l( i,j,1 ) = c1
               chgq_l( i,j,ndim_magmom ) = c2
            end do
         end do
      endif
    end subroutine redistrib_charge_density

    subroutine cpafft_CD_to_valencecharge(iloop)
      integer, intent(in) :: iloop
      integer :: i,j, ip, is
      real(kind=DP) :: rinplw

#ifdef _MPIFFTTEST_
      rinplw = 1.d0/product(fft_box_size_CD_c(1:3,1))
#else
      rinplw = 1.d0/product(fft_box_size_CD(1:3,1))
#endif

      chgq_l(:,:,iloop) = 0.d0
      do  j = 1, kimg
         do i = ista_kngp, iend_kngp
            if( kgp_reduced < i ) cycle
#ifdef _MPIFFT_
            ip = (igfp_nonpara(i)-1)*kimg+j
#elif _MPIFFTTEST_
            ip = (igfp_l_c(i)-1)*kimg + j
#else
            ip = (igfp_l(i)-1)*kimg + j
#endif
            chgq_l(i,j,iloop) = afft(ip)*rinplw
         end do
      end do

      total_charge = 0.d0
      if(mype == 0) total_charge = chgq_l(1,1,iloop)*univol

      if(npes > 1) call mpi_bcast(total_charge,1,mpi_double_precision,0,MPI_CommGroup,ierr)

!!$      if(dabs(total_charge-totch)>1.d-5) then
!!$         chgq_l(:,:,iloop) = chgq_l(:,:,iloop)*totch/total_charge
!!$      end if

      if(iprichargedensity >= 1) then
         write(nfout,'("* total_charge(",i3,") = ",d20.8," <<m_CD_initial_CD_by_file_rspace>>")') iloop,total_charge
      end if

      if(iprichargedensity >= 1) write(nfout,'("* total_charge = ",d20.8," <<m_CD_initial_CD_by_file_rspace>>")') total_charge
      if(iprichargedensity >= 2) then
         write(nfout,'(" !D chgq_l <<m_CD_initial_CD_by_file_rspace>>")')
         write(nfout,'(5f16.8)') (chgq_l(ista_kngp:min(ista_kngp+20,iend_kngp),1,1))
      end if

    end subroutine cpafft_CD_to_valencecharge

    subroutine rdchgr( nfout,nfchr )
      integer, intent(in) :: nfout, nfchr
      real(kind=DP),allocatable,dimension(:,:,:) :: wkchr
      integer :: idp,mmp,nlp,nmp,nnp,nlp_t,nmp_t,nnp_t, nfftp_t,inew,jnew,knew,nlphf
      integer :: i,j,k,ip,natm2_t
      real(kind=DP) :: totch_from_fft

      idp = fft_box_size_CD_nonpara(1,0)
      mmp = fft_box_size_CD_nonpara(2,0)
      nlp = fft_box_size_CD(1,1)
      nmp = fft_box_size_CD(2,1)
      nnp = fft_box_size_CD(3,1)

      nlphf = idp/(3-kimg)  ! = idp/2 (if kimg == 1), or idp (if kimg == 2)

      if(initial_charge_filetype == DENSITY_ONLY) then
         allocate(wkchr(nlp,nmp,nnp)); wkchr = 0.d0
      else if(charge_filetype == CUBE) then
         allocate(wkchr(nnp,nmp,nlp)); wkchr = 0.d0
      end if

      if(ipri >= 2) write(nfout,'(" !D Charge density ne = ",i8 &
           & ,"(",3i5," ) <<m_CD_initial_CD_by_file_rspace>>")') nlp*nmp*nnp, nlp, nmp, nnp
      if(ipri >= 2) write(nfout,'(" !D charge density is being read")')
      if(initial_charge_filetype == DENSITY_ONLY) then
         read(nfchr,'(21a,i8,1a,3i5)') nfftp_t, nlp_t, nmp_t, nnp_t
         call check_fftsize(nlp_t,nmp_t,nnp_t,nlp,nmp,nnp)
         read(nfchr,'(6e13.5)') wkchr
         if(iprichargedensity >= 2) &
              & write(nfout,'(" !CD initial_charge_filetype = DENSITY_ONLY")')
      else if(initial_charge_filetype == CUBE) then
         read(nfchr,*) ! comment
         read(nfchr,*) ! a comment of "SCF Total Density"
         read(nfchr,'(i6)') natm2_t ! natm2,x,y,z
         read(nfchr,'(i6)') nlp_t
         read(nfchr,'(i6)') nmp_t
         read(nfchr,'(i6)') nnp_t
         nfftp_t =  nlp_t*nmp_t*nnp_t
         call check_fftsize(nlp_t,nmp_t,nnp_t,nlp,nmp,nnp)
         do i = 1, natm2
            read(nfchr,*)
         end do
! --> T. Yamasaki, 1 Aug. 2008
         read(nfchr,*) wkchr
!!$         read(nfchr,'(6e13.5)') wkchr
!!$         do i = 1, nlp
!!$            do j = 1, nmp
!!$               read(nfchr,*) (wkchr(k,j,i),k=1,nnp)
!!$            end do
!!$         end do
         if(iprichargedensity>=2) then
            write(nfout,'(" -- wkchr --")')
            write(nfout,'(" wkchr(1:6,nmp,nlp)= ",6f13.5)') wkchr(1:min(6,nnp),nmp,nlp)
         end if
! <--
         if(iprichargedensity >= 1) &
              & write(nfout,'(" !CD initial_charge_filetype = CUBE")')
      else
         if(ipri >= 1) write(nfout,'(" initila_charge_filetype is invalid")')
         call phase_error_with_msg(nfout,' initial_charge_filetype is invalid << m_CD_initial_CD_by_file_rspace>>'&
                                  ,__LINE__,__FILE__)
      end if
!!$      if(wkchr(1,1,1) > DELTA10) then
!!$         if(ipri >= 1) write(nfout,'(" The charge density is analyzed as symmetric")')
!!$         up_down = UP
!!$      else


      if(iprichargedensity >= 2) then
         write(nfout,'(" !CD nfftp, nlp, nmp, nnp = ",4i8)') nfftp_t,nlp_t,nmp_t,nnp_t
         write(nfout,'(" !CD wkchr <<m_CD_initial_CD_file_rspace>>")')
         if(initial_charge_filetype == DENSITY_ONLY) then
            write(nfout,'(5f16.8)') wkchr(1:min(nlp,20),1,1)
         else if(initial_charge_filetype == CUBE) then
            write(nfout,'(5f16.8)') wkchr(1:min(nnp,20),1,1)
         end if
         write(nfout,'(" nfftp, nlp, nmp, nnp, nlphf = ",i10,4i8)') nfftp,nlp,nmp,nnp,nlphf
      end if

      afft = 0.d0
      do i = 1, nmp
         do j = 1, nnp
            do k = 1, nlp
               if(kimg == 1 .and. k > nlphf) then
                  knew = idp - k
                  jnew = nnp+2 - j
                  inew = nmp+2 - i
                  if(jnew > nnp) then
                     jnew = jnew - nnp
                  end if
                  if(inew > nmp) then
                     inew = inew - nmp
                  end if
               else
                  knew = k; jnew = j; inew = i
               end if
               ip = nlphf*mmp*(jnew-1) + nlphf*(inew-1) + knew
!!$               if(ip*2 > nfftp) then
!!$                  write(nfout,'(" !CD <<rdchgr>> ip*2=",i10," > nfftp")') ip*2
!!$                  stop ' ip*2 > nfftp'
!!$               end if
               if(charge_filetype == DENSITY_ONLY) then
                  afft(ip*2-1) = wkchr(k,i,j)
               else if(charge_filetype == CUBE) then
                  afft(ip*2-1) = wkchr(j,i,k)
               end if
            end do
         end do
      end do

      if(iprichargedensity >= 1 .and. (charge_filetype == DENSITY_ONLY .or. charge_filetype == CUBE)) then
         totch_from_fft = 0.d0
         if(charge_filetype == DENSITY_ONLY) then
            do j = 1, nnp
               do i = 1, nmp
                  do k = 1, nlp
                     totch_from_fft = totch_from_fft + wkchr(k,i,j)
                  end do
               end do
            end do
         else if(charge_filetype == CUBE) then
            do k = 1, nlp
               do i = 1, nmp
                  do j = 1, nnp
                     totch_from_fft = totch_from_fft + wkchr(j,i,k)
                  end do
               end do
            end do
         end if
         if(nnp*nmp*nlp /= 0) then
            totch_from_fft = totch_from_fft*univol/(nnp*nmp*nlp)
            write(nfout,'(" totch_from_fft = ",d20.8,"<<rdchgr>>")') totch_from_fft
         else
            write(nfout,'(" nnp*nmp*nlp = ",i8,"<<rdchgr>>")') nnp*nmp*nlp
         end if
      end if

      if(allocated(wkchr)) deallocate(wkchr)
      if(iprichargedensity >= 2) write(nfout,'(" end of <<rdchgr>>")')

    end subroutine rdchgr

    subroutine check_fftsize(nlp_t,nmp_t,nnp_t,nlp,nmp,nnp)
      integer, intent(in) :: nlp_t, nmp_t, nnp_t, nlp, nmp, nnp
      character(34) :: chsubroutinename = "<<m_CD_initial_CD_by_file_rspace>>"
      if(ipri >= 1) &
           & write(nfout,'(" nlp_t, nmp_t, nnp_t, nlp, nmp, nnp = ",6i6)') nlp_t,nmp_t,nnp_t,nlp,nmp,nnp
      if(nlp_t /= nlp) then
         if(ipri >=1 ) write(nfout,'(" nlp_t /= nlp ", a34)') chsubroutinename
      end if
      if(nmp_t /= nmp) then
         if(ipri >=1 ) write(nfout,'(" nmp_t /= nmp ", a34)') chsubroutinename
      end if
      if(nnp_t /= nnp) then
         if(ipri >=1 ) write(nfout,'(" nnp_t /= nnp ", a34)') chsubroutinename
      end if
      if(nlp_t /= nlp .or. nmp_t /= nmp .or. nnp_t /= nnp) then
         call phase_error_with_msg(nfout,&
         ' a set of (nlp,nmp,nnp) in the nfchr-file is invalid <<m_CD_initial_CD_by_file_rspace>>',__LINE__,__FILE__)
      end if
    end subroutine check_fftsize

  end subroutine m_CD_initCD_by_file_rsp_noncl
! ========================================================================= 11.0

  subroutine m_CD_rspace_charge(nspin,iloop,nfchr,nfout,add_core)
    integer,intent(in) :: nspin,iloop,nfchr,nfout
    logical,intent(in), optional :: add_core

    integer,dimension(3) :: boxsize
    real(kind=DP),dimension(3,3) :: cellsize
    integer, dimension(3,2) :: nind
    logical :: flg_add_core

    if ( present(add_core) ) then
       flg_add_core = add_core
    else
       flg_add_core = .false.
    endif

    call map_valence_charge_to_fft_box(iloop,flg_add_core)
    call m_FFT_CD_inverse0(nfout,afft)
!!$    call m_FFT_CD0(nfout,afft,INVERSE)
    if(mype == 0) call wdchgr(nfout,nfchr)

  contains
    subroutine wdchgr(nfout,nfchr)
      ! Revised according to an indication by Momita-san,
      !     T. Yamasaki (FUJITSU Laboratories ltd.), 2003/07/28
      !  The writing order of wkchr is reversed in the case of
      !  'charge_filetype == CUBE'.
      !
      integer, intent(in) :: nfout,nfchr
      integer :: i,j,k, idp, nlp, nmp, nnp, nlphf,inew,jnew,knew,ip,mmp
      real(kind=DP),allocatable,dimension(:,:,:) :: wkchr
      real(kind=DP) ::      s1, s2, sratio,x,y,z
      integer, parameter :: UP = 1 , DOWN = 2
      integer ::            up_down
      real(kind=DP),allocatable,dimension(:,:) :: cps_full,pos_full
      real(kind=DP), allocatable, dimension(:,:) :: rltv_t
      logical, allocatable, dimension(:) :: ignore
      integer :: natm3
      integer, allocatable,dimension(:) :: ityp_full
      real(kind=DP), dimension(3) :: r_wk    !!! K.Mae 040315
      integer :: ucret, m   !!! K.Mae 040315
      integer :: n1,n2,n3
      real(kind=DP) :: dn1,dn2,dn3
      real(kind=DP),dimension(3) :: trans

      idp = fft_box_size_CD_nonpara(1,0)
      mmp = fft_box_size_CD_nonpara(2,0)
      nlp = fft_box_size_CD(1,1)
      nmp = fft_box_size_CD(2,1)
      nnp = fft_box_size_CD(3,1)

      if(kimg == 1) then
         nlphf = idp/2
      else
         nlphf = idp
      end if

      if(charge_filetype == DENSITY_ONLY .or. charge_filetype == VTK) then
         allocate(wkchr(nlp,nmp,nnp)); wkchr = 0.d0
      else if(charge_filetype == CUBE) then
         allocate(wkchr(nnp,nmp,nlp)); wkchr = 0.d0
      end if

! -- checking of symmetric or anti-symmetric about the charge densities
      s1 = 0.d0; s2 = 0.d0
      do ip = 1, nlphf*mmp*nnp,2
         s1 = s1 + dabs(afft(ip))
         s2 = s2 + dabs(afft(ip+1))
      end do
      if(ipri >= 2) write(nfout,'(" s1 = ",d14.6," s2 = ",d14.6)') s1, s2
      if(s1 > s2) then
         if(ipri >= 2) &
              & write(nfout,*) ' The function of charge density is analyzed as symmetric.'
         up_down = UP
         sratio = s2/s1
      else
         if(ipri >= 2) &
              & write(nfout,*) ' The function of charge density is analyzed as anti-symmetric.'
         up_down = DOWN
         sratio = s1/s2
      endif
      if(ipri >= 2) write(nfout,*) ' !ratio = ', sratio

      if(ipri >= 2) write(nfout,9001) nlp*nmp*nnp, nlp, nmp, nnp
9001  format(' CHARGE DENSITY NE = ',i8,'(',3i5,')')

      if(ipri >= 2) write(nfout,*) ' !D FFT cube mapping start'
      do i = 1, nmp
         do j = 1, nnp
            do k = 1, nlp
               if(kimg == 1 .and. k > nlphf) then
                  knew = idp - k
                  jnew = nnp+2 - j
                  inew = nmp+2 - i
                  if(jnew > nnp) then
                     jnew = jnew - nnp
                  end if
                  if(inew > nmp) then
                     inew = inew - nmp
                  end if
               else
                  knew = k; jnew = j; inew = i
               end if
               ip = nlphf*mmp*(jnew-1) + nlphf*(inew-1) + knew
               if(charge_filetype == DENSITY_ONLY .or. charge_filetype == VTK) then
                  wkchr(k,i,j) = afft(ip*2-2+up_down)
               else if(charge_filetype == CUBE) then
                  wkchr(j,i,k) = afft(ip*2-2+up_down)
               end if
            end do
         end do
      end do

      if ( flg_add_core ) then
         if ( eval_corecharge_on_Gspace == OFF .and. charge_filetype == CUBE ) then
            call add_corecharge_rspace( nnp, nmp, nlp, wkchr, nspin, iloop )
         endif
      endif

      if(charge_filetype == DENSITY_ONLY) then
         write(nfchr,9001) nlp*nmp*nnp, nlp, nmp, nnp
         write(nfchr,'(6e13.5)') wkchr
      else if(charge_filetype == VTK) then
         write(nfchr,'("# vtk DataFile Version 2.0")')
         write(nfchr,'("Electronic-charge density data created by PHASE")')
         write(nfchr,'("ASCII")')
         write(nfchr,'("DATASET STRUCTURED_GRID")')
         write(nfchr,'("DIMENSIONS",3(1x,i5))') nlp+1,nmp+1,nnp+1
         write(nfchr,'("POINTS",1x,i7,1x,"float")') (nlp+1)*(nmp+1)*(nnp+1)
         do n1=0,nlp
            do n2=0,nmp
               do n3=0,nnp
                  dn1 = n1/dble(nlp)
                  dn2 = n2/dble(nmp)
                  dn3 = n3/dble(nnp)
                  x = altv(1,1)*dn1 + altv(1,2)*dn2 + altv(1,3)*dn3
                  y = altv(2,1)*dn1 + altv(2,2)*dn2 + altv(2,3)*dn3
                  z = altv(3,1)*dn1 + altv(3,2)*dn2 + altv(3,3)*dn3
                  write(nfchr,'(3(1x,e13.5))') x,y,z
               end do
            end do
         end do
         write(nfchr,'("")')
         write(nfchr,'("POINT_DATA",1x,i7)') (nlp+1)*(nmp+1)*(nnp+1)
         write(nfchr,'("SCALARS scalars float")')
         write(nfchr,'("LOOKUP_TABLE default")')
         do n1=0,nlp
            i=n1+1
             if(n1==nlp) i=1
            do n2=0,nmp
               j=n2+1
               if(n2==nmp) j=1
               do n3=0,nnp
                  k=n3+1
                  if(n3==nnp) k=1
                  write(nfchr,'(e13.5)') wkchr(i,j,k)
               end do
            end do
         end do
         !!$write(nfchr,'("")')
         !!$write(nfchr,'("# Unit cell vectors")')
         !!$do i=1,3
         !!$   write(nfchr,'(3(1x,f10.6))') altv(1:3,i)
         !!$end do
         !!$write(nfchr,'("")')
         !!$write(nfchr,'("# Atomic structure")')
         !!$allocate(cps_full(natm2,3))
         !!$allocate(ityp_full(natm2))
         !!$call m_IS_pack_all_ions_in_uc(ityp_full,cps_full)
         !!$do i = 1, natm2
         !!$   m = ityp_full(i)
         !!$   write(nfchr,'(a4,3(1x,f10.6))') speciesname(m), cps_full(i,1:3)
         !!$end do
         !!$deallocate(ityp_full,cps_full)
      else if(charge_filetype == CUBE) then
         if(sw_subset_only==ON) call build_adjusted_chr_index(nlp,nmp,nnp)
         if(len_trim(charge_title) >= 1) then
            write(nfchr,*) trim(charge_title)
         else
            write(nfchr,'(" Calculated by phase")')
         end if
         if(nspin == 2) then
            if(iloop == 1) then
               write(nfchr,'(" SCF Total Density  UP")')
            else if ( iloop == 2 ) then
               write(nfchr,'(" SCF Total Density  DOWN")')
            end if
! ====== KT_add ===== 2014/06/07
            if(iloop == -1) then
               write(nfchr,'(" SCF Total Density")')
            else if ( iloop == -2 ) then
               write(nfchr,'(" SCF Spin Manetic Moment Density")')
            end if
! =================== 2014/06/07
         else
            write(nfchr,'(" SCF Total Density")')
         end if
         allocate(cps_full(natm2,3))
         allocate(ityp_full(natm2))
         call m_IS_pack_all_ions_in_uc(ityp_full,cps_full)
         natm3 = natm2
         if(sw_subset_only==ON)then
            allocate(pos_full(natm2,3))
            allocate(rltv_t(3,3))
            rltv_t = transpose(rltv)/PAI2
            call change_of_coordinate_system(rltv_t,cps_full,natm2,natm2,pos_full)
            allocate(ignore(natm2))
            natm3 = 0
            do i=1,natm2
               ignore(i) = (pos_full(i,1).lt.minxyz(1).or.(pos_full(i,1).gt.maxxyz(1).and.maxxyz(1).gt.0)) &
      &               .or. (pos_full(i,2).lt.minxyz(2).or.(pos_full(i,2).gt.maxxyz(2).and.maxxyz(2).gt.0)) &
      &               .or. (pos_full(i,3).lt.minxyz(3).or.(pos_full(i,3).gt.maxxyz(3).and.maxxyz(3).gt.0))
               if (.not.ignore(i)) natm3 = natm3+1
            enddo
         endif
         x = 0.d0; y = 0.d0; z = 0.d0
         write(nfchr,'(i6,3f15.4)') natm3, x,y,z
         if(sw_subset_only==ON)then
            do i=1,3
               boxsize(i) = nind(i,2)-nind(i,1)+1
            enddo
         else
            do i=1,3
               boxsize(i) = fft_box_size_CD(i,1)
            enddo
         endif
         do i = 1, 3
            !!!! K.Mae 040315
            !!do m = 1, 3
            !!   ucret = unit_conv_byname( altv(m,i), r_wk(m), 'bohr', 'angstrom' )
            !!end do
            !!write(nfchr,'(i6,3f10.6)') fft_box_size_CD(i,1), r_wk(1:3)/dble(fft_box_size_CD(i,1))
            !!write(nfchr,'(i6,3f10.6)') fft_box_size_CD(i,1), altv(1:3,i)/dble(fft_box_size_CD(i,1))
            write(nfchr,'(i6,3f25.15)') boxsize(i), altv(1:3,i)/dble(fft_box_size_CD(i,1))
            !!!! end K.Mae 040315
         end do

         trans = 0.d0
         if(sw_subset_only==ON)then
            do i=1,3
               do j=1,3
                  if(minxyz(j).ge.0) trans(i) = trans(i) - altv(i,j)*minxyz(j)
               enddo
            enddo
         endif
         do i = 1, natm2
            if(sw_subset_only==ON)then
               if(ignore(i)) cycle
            endif
            m = ityp_full(i)
            cps_full(i,:) = cps_full(i,:)+trans(:)
!!$            write(nfchr,'(i6,4f10.6)') nint(iatomn(m)), ival(m), cps_full(i,1:3)
            if ( flg_add_core ) then
               write(nfchr,'(i6,1x,f10.6,3(1x,f18.6))') &
                    &        nint(iatomn(m)), iatomn(m), cps_full(i,1:3)
            else
               write(nfchr,'(i6,1x,f10.6,3(1x,f18.6))') &
                    &        nint(iatomn(m)), ival(m), cps_full(i,1:3)
            endif
         end do
         deallocate(ityp_full,cps_full)
         if(sw_subset_only==ON)then
            deallocate(ignore)
            deallocate(pos_full)
            deallocate(rltv_t)
         endif

! --> T. Yamasaki, 1 Aug. 2008
         if(sw_subset_only==OFF)then
            do i = 1, nlp
               do j = 1, nmp
                  write(nfchr,'(6e25.15)') (wkchr(k,j,i),k=1,nnp)
               end do
            end do
         else
            do i = nind(1,1),nind(1,2)
               do j = nind(2,1),nind(2,2)
                  write(nfchr,'(6e25.15)') (wkchr(k,j,i),k=nind(3,1),nind(3,2))
               end do
            end do
         endif
!!$         write(nfchr,'(6e13.5)') wkchr
! <--
      end if
      if(allocated(wkchr)) deallocate(wkchr)
    end subroutine wdchgr

    subroutine map_valence_charge_to_fft_box(iloop,flg_add_core)
      integer, intent(in) :: iloop
      logical, intent(in) :: flg_add_core
      integer :: j, i, ip, it

      afft_mpi1 = 0.d0
      do j = 1, kimg
         do i = ista_kngp, iend_kngp
            ip = (igfp_l(i)-1)*kimg + j
! === KT_mod === 2014/06/07
!            afft_mpi1(ip) = afft_mpi1(ip) + chgq_l(i,j,iloop)
!
            if ( iloop > 0 ) then
               afft_mpi1(ip) = afft_mpi1(ip) + chgq_l(i,j,iloop)
            else
               if ( iloop == -1 ) then
                  afft_mpi1(ip) = afft_mpi1(ip) + chgq_l(i,j,1) +chgq_l(i,j,2)
               else if ( iloop == -2 ) then
                  afft_mpi1(ip) = afft_mpi1(ip) + chgq_l(i,j,1) -chgq_l(i,j,2)
               endif
            endif
! ============== 2014/06/07
         end do
      end do

      if ( flg_add_core .and. eval_corecharge_on_Gspace == ON ) then
         Do it=1, ntyp
            do j = 1, kimg
               do i = ista_kngp, iend_kngp
                  if(kgp_reduced < i) cycle
#ifdef _MPIFFT_
                  ip = (igfp_nonpara(i)-1)*kimg + j
#else
                  ip = (igfp_l(i)-1)*kimg + j
#endif
                  if ( iloop > 0 ) then
                     afft_mpi1(ip) = afft_mpi1(ip) + zfm3_l(i,it,j)*rhcg_l(i,it)/nspin
                  else
                     if ( iloop == -1 ) then
                        afft_mpi1(ip) = afft_mpi1(ip) + zfm3_l(i,it,j)*rhcg_l(i,it)
                     endif
                  endif

               end do
            end do
         end Do
      end if

      if(.not.allocated(afft)) allocate(afft(nfftp_nonpara))
      if(npes >= 2) then
         call mpi_allreduce(afft_mpi1,afft,nfftp_nonpara,mpi_double_precision &
              &  ,mpi_sum,mpi_chg_world,ierr)
      else
         afft = afft_mpi1
      end if
    end subroutine map_valence_charge_to_fft_box

    subroutine build_adjusted_chr_index(nlp,nmp,nnp)
      integer, intent(in) :: nlp,nmp,nnp
      real(kind=DP) :: x,y,z
      integer :: n1,n2,n3,i
      real(kind=DP) :: dn1,dn2,dn3
      do i=1,3
         nind(i,1) = 1
      enddo
      nind(1,2) = nlp
      nind(2,2) = nmp
      nind(3,2) = nnp
      if(minxyz(1).gt.0) nind(1,1) = floor(minxyz(1)*nlp)+1
      if(minxyz(2).gt.0) nind(2,1) = floor(minxyz(2)*nmp)+1
      if(minxyz(3).gt.0) nind(3,1) = floor(minxyz(3)*nnp)+1
      if(maxxyz(1).gt.0) nind(1,2) = floor(maxxyz(1)*nlp)+1
      if(maxxyz(2).gt.0) nind(2,2) = floor(maxxyz(2)*nmp)+1
      if(maxxyz(3).gt.0) nind(3,2) = floor(maxxyz(3)*nnp)+1
      if(mype==0 .and. iprichargedensity>=1) then
         write(nfout,'(a)') ' -- cube subset indices --'
         write(nfout,'(a,i8,a,i8)') '  a-axis : ',nind(1,1),' to ',nind(1,2)
         write(nfout,'(a,i8,a,i8)') '  b-axis : ',nind(2,1),' to ',nind(2,2)
         write(nfout,'(a,i8,a,i8)') '  c-axis : ',nind(3,1),' to ',nind(3,2)
      endif
    end subroutine build_adjusted_chr_index

  end subroutine m_CD_rspace_charge

  subroutine m_CD_get_rspace_charge(nfout,na,nb,nc,rho,is)
    integer,intent(in) :: na,nb,nc,nfout,is
    real(kind=8), dimension(na,nb,nc), intent(out) :: rho
    integer,dimension(3) :: boxsize
    real(kind=DP),dimension(3,3) :: cellsize
    integer, dimension(3,2) :: nind
    integer :: iloop

    call m_CD_alloc_rspace_charge()
    !do iloop=1,nspin
      call map_valence_charge_to_fft_box(is)
      call m_FFT_CD_inverse0(nfout,afft)
!!$    call m_FFT_CD0(nfout,afft,INVERSE)
      call chgr()
    !enddo
    call m_CD_dealloc_rspace_charge()

    contains

    subroutine chgr()
      ! Revised according to an indication by Momita-san,
      !     T. Yamasaki (FUJITSU Laboratories ltd.), 2003/07/28
      !  The writing order of wkchr is reversed in the case of
      !  'charge_filetype == CUBE'.
      !
      integer :: i,j,k, idp, nlp, nmp, nnp, nlphf,inew,jnew,knew,ip,mmp
      real(kind=DP),allocatable,dimension(:,:,:) :: wkchr
      real(kind=DP) ::      s1, s2, sratio,x,y,z
      integer, parameter :: UP = 1 , DOWN = 2
      integer ::            up_down
      real(kind=DP),allocatable,dimension(:,:) :: cps_full,pos_full
      real(kind=DP), allocatable, dimension(:,:) :: rltv_t
      logical, allocatable, dimension(:) :: ignore
      integer :: natm3
      integer, allocatable,dimension(:) :: ityp_full
      real(kind=DP), dimension(3) :: r_wk    !!! K.Mae 040315
      integer :: ucret, m   !!! K.Mae 040315
      integer :: n1,n2,n3
      real(kind=DP) :: dn1,dn2,dn3
      real(kind=DP),dimension(3) :: trans

      idp = fft_box_size_CD_nonpara(1,0)
      mmp = fft_box_size_CD_nonpara(2,0)
      nlp = fft_box_size_CD(1,1)
      nmp = fft_box_size_CD(2,1)
      nnp = fft_box_size_CD(3,1)

      if(kimg == 1) then
         nlphf = idp/2
      else
         nlphf = idp
      end if

! -- checking of symmetric or anti-symmetric about the charge densities
      s1 = 0.d0; s2 = 0.d0
      do ip = 1, nlphf*mmp*nnp,2
         s1 = s1 + dabs(afft(ip))
         s2 = s2 + dabs(afft(ip+1))
      end do
      if(ipri >= 2) write(nfout,'(" s1 = ",d14.6," s2 = ",d14.6)') s1, s2
      if(s1 > s2) then
         if(ipri >= 2) &
              & write(nfout,*) ' The function of charge density is analyzed as symmetric.'
         up_down = UP
         sratio = s2/s1
      else
         if(ipri >= 2) &
              & write(nfout,*) ' The function of charge density is analyzed as anti-symmetric.'
         up_down = DOWN
         sratio = s1/s2
      endif
      if(ipri >= 2) write(nfout,*) ' !ratio = ', sratio

      if(ipri >= 2) write(nfout,9001) nlp*nmp*nnp, nlp, nmp, nnp
9001  format(' CHARGE DENSITY NE = ',i8,'(',3i5,')')

      if(ipri >= 2) write(nfout,*) ' !D FFT cube mapping start'
      do i = 1, nmp
         do j = 1, nnp
            do k = 1, nlp
               if(kimg == 1 .and. k > nlphf) then
                  knew = idp - k
                  jnew = nnp+2 - j
                  inew = nmp+2 - i
                  if(jnew > nnp) then
                     jnew = jnew - nnp
                  end if
                  if(inew > nmp) then
                     inew = inew - nmp
                  end if
               else
                  knew = k; jnew = j; inew = i
               end if
               ip = nlphf*mmp*(jnew-1) + nlphf*(inew-1) + knew
               rho(k,i,j) = afft(ip*2-2+up_down)
            end do
         end do
      end do
    end subroutine chgr

    subroutine map_valence_charge_to_fft_box(iloop)
      integer, intent(in) :: iloop
      integer :: j, i, ip, it, mm

      afft_mpi1 = 0.d0
      do j = 1, kimg
         do i = ista_kngp, iend_kngp
            ip = (igfp_l(i)-1)*kimg + j
! === KT_mod === 2014/06/07
!            afft_mpi1(ip) = afft_mpi1(ip) + chgq_l(i,j,iloop)
!
            if ( iloop > 0 ) then
               afft_mpi1(ip) = afft_mpi1(ip) + chgq_l(i,j,iloop)
            else
               if ( iloop == -1 ) then
                  afft_mpi1(ip) = afft_mpi1(ip) + chgq_l(i,j,1) +chgq_l(i,j,2)
               else if ( iloop == -2 ) then
                  afft_mpi1(ip) = afft_mpi1(ip) + chgq_l(i,j,1) -chgq_l(i,j,2)
               endif
            endif
! ============== 2014/06/07
         end do
      end do
      mm = 0.d0
      do it = 1, ntyp
         if(itpcc(it) == 0) cycle
         mm = mm + 1
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
               afft_mpi1(ip) = afft_mpi1(ip) &   !mpi
                    & + zfm3_l(i,it,j)*rhpcg_l(i,mm)   !mpi
            end do
         end do
      end do
      if(npes >= 2) then
         call mpi_allreduce(afft_mpi1,afft,nfftp_nonpara,mpi_double_precision &
              &  ,mpi_sum,mpi_chg_world,ierr)
      else
         afft = afft_mpi1
      end if
    end subroutine map_valence_charge_to_fft_box

  end subroutine m_CD_get_rspace_charge

! =================================== added by K. Tagami ============= 11.0
  subroutine m_CD_rspace_charge_noncl( iloop, nfchr, nfout, charge_filetype, add_core )
    integer,intent(in) :: nfchr, iloop, nfout, charge_filetype
    logical,intent(in), optional :: add_core
    logical :: flg_add_core

    if ( present(add_core) ) then
       flg_add_core = add_core
    else
       flg_add_core = .false.
    endif

    call map_valcharge_to_fft_box(iloop)
    call m_FFT_CD_inverse0(nfout,afft)          ! G-space -> R-space
    if ( mype == 0 ) call wdchgr_noncl( iloop, nfout, nfchr )

  contains

    subroutine wdchgr_noncl( iloop, nfout, nfchr )
      ! Revised according to an indication by Momita-san,
      !     T. Yamasaki (FUJITSU Laboratories ltd.), 2003/07/28
      !  The writing order of wkchr is reversed in the case of
      !  'charge_filetype == CUBE'.
      !
      integer, intent(in) :: nfout,nfchr
      integer, intent(in) :: iloop

      integer :: i,j,k, idp, nlp, nmp, nnp, nlphf,inew,jnew,knew,ip,mmp
      real(kind=DP),allocatable,dimension(:,:,:) :: wkchr
      real(kind=DP) ::      s1, s2, sratio,x,y,z
      integer, parameter :: UP = 1 , DOWN = 2
      integer ::            up_down
      real(kind=DP),allocatable,dimension(:,:) :: cps_full
      integer, allocatable,dimension(:) :: ityp_full
      real(kind=DP), dimension(3) :: r_wk    !!! K.Mae 040315
      integer :: ucret, m   !!! K.Mae 040315
      integer :: n1,n2,n3
      real(kind=DP) :: dn1,dn2,dn3

      idp = fft_box_size_CD_nonpara(1,0)
      mmp = fft_box_size_CD_nonpara(2,0)
      nlp = fft_box_size_CD(1,1)
      nmp = fft_box_size_CD(2,1)
      nnp = fft_box_size_CD(3,1)

      if(kimg == 1) then
         nlphf = idp/2
      else
         nlphf = idp
      end if

      if(charge_filetype == DENSITY_ONLY .or. charge_filetype == VTK) then
         allocate(wkchr(nlp,nmp,nnp)); wkchr = 0.d0
      else if(charge_filetype == CUBE) then
         allocate(wkchr(nnp,nmp,nlp)); wkchr = 0.d0
      end if

! -- checking of symmetric or anti-symmetric about the charge densities
      s1 = 0.d0; s2 = 0.d0
      do ip = 1, nlphf*mmp*nnp,2
         s1 = s1 + dabs(afft(ip))
         s2 = s2 + dabs(afft(ip+1))
      end do
      if(ipri >= 2) write(nfout,'(" s1 = ",d14.6," s2 = ",d14.6)') s1, s2
      if(s1 > s2) then
         if(ipri >= 2) &
              & write(nfout,*) ' The function of charge density is analyzed as symmetric.'
         up_down = UP
         sratio = s2/s1
      else
         if(ipri >= 2) &
              & write(nfout,*) ' The function of charge density is analyzed as anti-symmetric.'
         up_down = DOWN
         sratio = s1/s2
      endif
      if(ipri >= 2) write(nfout,*) ' !ratio = ', sratio

      if(ipri >= 2) write(nfout,9001) nlp*nmp*nnp, nlp, nmp, nnp
9001  format(' CHARGE DENSITY NE = ',i8,'(',3i5,')')

      if(ipri >= 2) write(nfout,*) ' !D FFT cube mapping start'
      do i = 1, nmp
         do j = 1, nnp
            do k = 1, nlp
               if(kimg == 1 .and. k > nlphf) then
                  knew = idp - k
                  jnew = nnp+2 - j
                  inew = nmp+2 - i
                  if(jnew > nnp) then
                     jnew = jnew - nnp
                  end if
                  if(inew > nmp) then
                     inew = inew - nmp
                  end if
               else
                  knew = k; jnew = j; inew = i
               end if
               ip = nlphf*mmp*(jnew-1) + nlphf*(inew-1) + knew
               if(charge_filetype == DENSITY_ONLY .or. charge_filetype == VTK) then
                  wkchr(k,i,j) = afft(ip*2-2+up_down)
               else if(charge_filetype == CUBE) then
                  wkchr(j,i,k) = afft(ip*2-2+up_down)
               end if
            end do
         end do
      end do

      if(charge_filetype == DENSITY_ONLY) then
         write(nfchr,9001) nlp*nmp*nnp, nlp, nmp, nnp
         write(nfchr,'(6e13.5)') wkchr
      else if(charge_filetype == VTK) then
         write(nfchr,'("# vtk DataFile Version 2.0")')
         write(nfchr,'("Electronic-charge density data created by PHASE")')
         write(nfchr,'("ASCII")')
         write(nfchr,'("DATASET STRUCTURED_GRID")')
         write(nfchr,'("DIMENSIONS",3(1x,i5))') nlp+1,nmp+1,nnp+1
         write(nfchr,'("POINTS",1x,i7,1x,"float")') (nlp+1)*(nmp+1)*(nnp+1)
         do n1=0,nlp
            do n2=0,nmp
               do n3=0,nnp
                  dn1 = n1/dble(nlp)
                  dn2 = n2/dble(nmp)
                  dn3 = n3/dble(nnp)
                  x = altv(1,1)*dn1 + altv(1,2)*dn2 + altv(1,3)*dn3
                  y = altv(2,1)*dn1 + altv(2,2)*dn2 + altv(2,3)*dn3
                  z = altv(3,1)*dn1 + altv(3,2)*dn2 + altv(3,3)*dn3
                  write(nfchr,'(3(1x,e13.5))') x,y,z
               end do
            end do
         end do
         write(nfchr,'("")')
         write(nfchr,'("POINT_DATA",1x,i7)') (nlp+1)*(nmp+1)*(nnp+1)
         write(nfchr,'("SCALARS scalars float")')
         write(nfchr,'("LOOKUP_TABLE default")')
         do n1=0,nlp
            i=n1+1
             if(n1==nlp) i=1
            do n2=0,nmp
               j=n2+1
               if(n2==nmp) j=1
               do n3=0,nnp
                  k=n3+1
                  if(n3==nnp) k=1
                  write(nfchr,'(e13.5)') wkchr(i,j,k)
               end do
            end do
         end do
         !!$write(nfchr,'("")')
         !!$write(nfchr,'("# Unit cell vectors")')
         !!$do i=1,3
         !!$   write(nfchr,'(3(1x,f10.6))') altv(1:3,i)
         !!$end do
         !!$write(nfchr,'("")')
         !!$write(nfchr,'("# Atomic structure")')
         !!$allocate(cps_full(natm2,3))
         !!$allocate(ityp_full(natm2))
         !!$call m_IS_pack_all_ions_in_uc(ityp_full,cps_full)
         !!$do i = 1, natm2
         !!$   m = ityp_full(i)
         !!$   write(nfchr,'(a4,3(1x,f10.6))') speciesname(m), cps_full(i,1:3)
         !!$end do
         !!$deallocate(ityp_full,cps_full)
      else if(charge_filetype == CUBE) then
         if(len_trim(charge_title) >= 1) then
            write(nfchr,*) trim(charge_title)
         else
            write(nfchr,'(" Calculated by phase")')
         end if
! -----
         select case (iloop)
         case(1)
               write(nfchr,'(" SCF Total Density")')
         case(2)
               write(nfchr,'(" SCF Magnetic Density : X" )')
         case(3)
               write(nfchr,'(" SCF Magnetic Density : Y" )')
         case(4)
               write(nfchr,'(" SCF Magnetic Density : Z" )')
         end select
! --
         x = 0.d0; y = 0.d0; z = 0.d0
         write(nfchr,'(i6,3f10.4)') natm2, x,y,z
         do i = 1, 3
            !!!! K.Mae 040315
            !!do m = 1, 3
            !!   ucret = unit_conv_byname( altv(m,i), r_wk(m), 'bohr', 'angstrom' )
            !!end do
            !!write(nfchr,'(i6,3f10.6)') fft_box_size_CD(i,1), r_wk(1:3)/dble(fft_box_size_CD(i,1))
            write(nfchr,'(i6,3f10.6)') fft_box_size_CD(i,1), altv(1:3,i)/dble(fft_box_size_CD(i,1))
            !!!! end K.Mae 040315
         end do

         allocate(cps_full(natm2,3))
         allocate(ityp_full(natm2))
         call m_IS_pack_all_ions_in_uc(ityp_full,cps_full)
         do i = 1, natm2
            m = ityp_full(i)
!!$            write(nfchr,'(i6,4f10.6)') nint(iatomn(m)), ival(m), cps_full(i,1:3)
            if ( flg_add_core ) then
               write(nfchr,'(i6,1x,f10.6,3(1x,f18.6))') &
                    &        nint(iatomn(m)), iatomn(m), cps_full(i,1:3)
            else
               write(nfchr,'(i6,1x,f10.6,3(1x,f18.6))') &
                    &        nint(iatomn(m)), ival(m), cps_full(i,1:3)
            endif
         end do
         deallocate(ityp_full,cps_full)

         do i = 1, nlp
            do j = 1, nmp
               write(nfchr,'(6e13.5)') (wkchr(k,j,i),k=1,nnp)
            end do
         end do
!!$         write(nfchr,'(6e13.5)') wkchr

      end if
      if(allocated(wkchr)) deallocate(wkchr)
    end subroutine wdchgr_noncl

    subroutine map_valcharge_to_fft_box(iloop)
      integer, intent(in) :: iloop
      integer :: j, i, ip, it

      afft_mpi1 = 0.d0
      do j = 1, kimg
         do i = ista_kngp, iend_kngp
            if(kgp_reduced < i) cycle
#ifdef _MPIFFT_
            ip = (igfp_nonpara(i)-1)*kimg + j
#else
            ip = (igfp_l(i)-1)*kimg + j
#endif
            afft_mpi1(ip) = afft_mpi1(ip) + chgq_l(i,j,iloop)
         end do
      end do

      if ( flg_add_core .and. eval_corecharge_on_Gspace == ON ) then
         Do it=1, ntyp
            do j = 1, kimg
               do i = ista_kngp, iend_kngp
                  if(kgp_reduced < i) cycle
#ifdef _MPIFFT_
                  ip = (igfp_nonpara(i)-1)*kimg + j
#else
                  ip = (igfp_l(i)-1)*kimg + j
#endif
                  if ( iloop == 1 ) then
                     afft_mpi1(ip) = afft_mpi1(ip) + zfm3_l(i,it,j)*rhcg_l(i,it)
                  endif
               end do
            end do
         end Do
      end if

      if(npes >= 2) then
         call mpi_allreduce(afft_mpi1,afft,nfftp_nonpara,mpi_double_precision &
              &  ,mpi_sum,MPI_CommGroup,ierr)
      else
         afft = afft_mpi1
      end if
    end subroutine map_valcharge_to_fft_box

  end subroutine m_CD_rspace_charge_noncl
! ======================================================================= 11.0

  subroutine add_corecharge_rspace( nc, nb, na, rho, nspin, iloop )  ! cube only  ??
    integer, intent(in) :: na, nb, nc, nspin, iloop
    real(kind=DP), intent(inout) :: rho(nc,nb,na)

    integer :: i, ind, cix, ciy, ciz
    integer :: ia, it, nx, ny, nz
    integer :: ir, ier
    real(kind=DP) :: aa(3,3), vec1(3), vec2(3), vec3(3), dist
    real(kind=DP) :: f1, f2, c1, ctmp, fac
    real(kind=DP) :: dist_min = 1.0D-16
    real(kind=DP) :: dist_max = 3.0d0
    real(kind=DP) :: csum

    fac = 1.0d0
    if ( nspin == 2 ) then
       if ( iloop > 0 ) fac = 0.5d0
       if ( iloop == -2 ) return
    endif
    fac = fac /PAI /4.0d0

    do i=1,3
       aa(i,1:3) = altv(1:3,i)/dble(fft_box_size_CD(i,1))
    enddo

    Do ia=1, natm
       it = ityp(ia)

       Do nx=-1, 1
          Do ny=-1, 1
             Do nz=-1, 1
                ! ----
                vec1(1) = pos(ia,1) -floor( pos(ia,1) ) +nx
                vec1(2) = pos(ia,2) -floor( pos(ia,2) ) +ny
                vec1(3) = pos(ia,3) -floor( pos(ia,3) ) +nz

                Do cix=1, na
                   Do ciy=1, nb
                      Do ciz=1, nc
                         vec2(1) = dble( cix -1 )/dble(fft_box_size_CD(1,1) )
                         vec2(2) = dble( ciy -1 )/dble(fft_box_size_CD(2,1) )
                         vec2(3) = dble( ciz -1 )/dble(fft_box_size_CD(3,1) )

                         vec3(1:3) = ( vec2(1) -vec1(1) )*altv(1:3,1) &
                              &  + ( vec2(2) -vec1(2) )*altv(1:3,2)  &
                              &  + ( vec2(3) -vec1(3) )*altv(1:3,3)
                         dist = sqrt( vec3(1)**2 +vec3(2)**2 +vec3(3)**2 )
                         if ( dist > dist_max ) cycle

                         if ( dist < dist_min ) then
                            ind = 1
                            c1 = rhcorpw(ind,it) /radr_paw(ind,it)**2
                         else
                            ctmp = log( dist /rmax(it) ) *xh(it) +nmesh(it)
                            ind = int(ctmp)
                            if ( ind < 1 ) then
                               ind = 1
                               c1 = rhcorpw(ind,it) /radr_paw(ind,it)**2
                            else
                               f1 = ctmp -int(ctmp)
                               f2 = 1.0d0 -f1
                               c1 = ( rhcorpw(ind,it)/radr_paw(ind,it)**2 )**f2 &
                                    & *( rhcorpw(ind+1,it)/radr_paw(ind+1,it)**2) **f1
                            endif
                         endif
!                         c1 = c1 *fac /dist**2
                         c1 = c1 *fac
                         rho(ciz,ciy,cix) = rho(ciz,ciy,cix) +c1
                      End do
                   End Do
                End Do
                ! ----
             ENd Do
          End do
       End do
    End Do
  end subroutine add_corecharge_rspace

  subroutine m_CD_map_valence_charge_to_fft_box(iloop,nfft,aff)
    integer, intent(in) :: iloop
    integer, intent(in) :: nfft
    real(kind=DP), dimension(nfft), intent(out) :: aff
    real(kind=DP), allocatable, dimension(:) :: affttmp
    integer :: j, i, ip

    allocate(affttmp(nfft))
    affttmp = 0.d0
    aff = 0.d0
    do j = 1, kimg
       do i = ista_kngp, iend_kngp
          ip = (igfp_l(i)-1)*kimg + j
          affttmp(ip) = affttmp(ip) + chgq_l(i,j,iloop)
       end do
    end do

    if(npes >= 2) then
!!$       call mpi_allreduce(affttmp,aff,nfftp,mpi_double_precision &
       call mpi_allreduce(affttmp,aff,nfft,mpi_double_precision &
            &  ,mpi_sum,mpi_chg_world,ierr)
    else
       aff = affttmp
    end if
    deallocate(affttmp)
  end subroutine m_CD_map_valence_charge_to_fft_box

  subroutine m_CD_hardpart_sub(nfout,is,ik,ib,chgq0)
    integer, intent(in) :: nfout, is,ik, ib
    real(kind=DP), intent(out) :: chgq0

    real(kind=DP) :: z
    integer :: id_sname = -1
    call tstatc0_begin('m_CD_hardpart_sub ',id_sname,1)

    call summation_of_ff_sub(is,ik,ib) ! -(m_C.D.) (vnlsum) fsr_l, fsi_l, occup_l --> hsr

    chgq_l = 0.d0

! ================================ modidied by K. Tagami =========== 11.0
!    call add_hardpart_to_chgq_l(nfout,is,hsr) ! (lclchg) add hardpart and make average
#ifdef DEBUG_LDOS
    write(6,'(" size of hsr  = ",4i8)') &
         &              ubound(hsr,1)-lbound(hsr,1)+1 &, ubound(hsr,2)-lbound(hsr,2)+1 &
         &             ,ubound(hsr,3)-lbound(hsr,3)+1, ubound(hsr,4)-lbound(hsr,4)+1
    write(6,'(" size of dl2p = ",4i8)') &
         &              ubound(dl2p,1)-lbound(dl2p,1)+1, ubound(dl2p,2)-lbound(dl2p,2)+1 &
         &             ,ubound(dl2p,3)-lbound(dl2p,3)+1, ubound(dl2p,4)-lbound(dl2p,4)+1
#endif
    call add_hardpart_to_chgq_l_in_keworld( nfout, is, hsr) ! -> chgq_l
! =================================================================== 11.0

!!$    if(iprichargedensity >= 2) then
!!$       write(nfout,'(" -- hardpart summed --")')
!!$       call m_CD_wd_chgq_l_small_portion(nfout)
!!$    end if

    if(npes > 1) then
       chgq0 = 0.d0
       if(ista_kngp <= 1 .and. iend_kngp >= 1) then
!============================================ modified by K. Tagami ===== 11.0
!          chgq0 = chgq_l(1,1,1)
          chgq0 = chgq_l(1,1,is)
! ======================================================================= 11.0
       end if
!       call mpi_allreduce(chgq0,z,1,mpi_double_precision,mpi_sum,mpi_chg_world,ierr)
       call mpi_allreduce(chgq0,z,1,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       chgq0 = z
    else
!============================================ modified by K. Tagami ===== 11.0
!       chgq0 = chgq_l(1,1,1)
       chgq0 = chgq_l(1,1,is)
! ======================================================================= 11.0
    end if

    call tstatc0_end(id_sname)
  end subroutine m_CD_hardpart_sub


  subroutine m_CD_restore_chgq()
    chgq_l = chgqo_l
  end subroutine m_CD_restore_chgq

  subroutine summation_of_ff_sub(is,ik,ib)
    integer, intent(in) :: is, ik, ib

    real(kind=DP)   :: w_n
    integer         :: ia, lmt1, lmt2, p, q, it
    real(kind=DP),pointer, dimension(:,:) :: fs
    integer :: i, iadd

    if(k_symmetry(ik) == GAMMA) then
       allocate(fs(nlmta,1))
       fs = 0.d0
       do i = ista_fs, iend_fs
          iadd = i - ista_fs + 1
!!$          fs(i,1) = fsr_l(map_z(ib),iadd,ik)
          fs(i,1) = fsr_l(ib,iadd,ik)
       end do
       w_n = 2.d0
       hsr(:,:,:,is) = 0.d0
       do ia = 1, natm
          it = ityp(ia)
          do lmt1 = 1, ilmt(it)
             p = lmta(lmt1,ia)
             if(ista_fs<=p .and. p<=iend_fs) then
                do lmt2 = lmt1, ilmt(it)
                   q = lmta(lmt2,ia)
                   if(ista_fs<=q .and. p<=iend_fs) then
                      hsr(ia,lmt1,lmt2,is) = hsr(ia,lmt1,lmt2,is) + w_n * ( fs(p,1)*fs(q,1))
                   end if
                end do! lmt2
             end if
          end do! lmt1
       end do! ia
    else
       allocate(fs(nlmta,2))
       fs = 0.d0
       do i = ista_fs, iend_fs
!!$       if(map_ek(ib,ik) == myrank_e+myrank_k*nrank_e) then
!!$          do i = ista_fs, iend_fs
          iadd = i - ista_fs + 1
!!$          fs(i,1) = fsr_l(map_z(ib),iadd,ik)
!!$          fs(i,2) = fsi_l(map_z(ib),iadd,ik)
          fs(i,1) = fsr_l(ib,iadd,ik)
          fs(i,2) = fsi_l(ib,iadd,ik)
       end do
!!$          call mpi_allreduce(MPI_IN_PLACE,fs,nlmta*2,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
!!$       end if
       w_n = 2.d0
       hsr(:,:,:,is) = 0.d0
       do ia = 1, natm
          it = ityp(ia)
          do lmt1 = 1, ilmt(it)
             p = lmta(lmt1,ia)
             if(ista_fs<=p .and. p<=iend_fs) then
                do lmt2 = lmt1, ilmt(it)
                   q = lmta(lmt2,ia)
                   if(ista_fs<=q .and. p<=iend_fs) then
                      hsr(ia,lmt1,lmt2,is) = hsr(ia,lmt1,lmt2,is) + w_n * ( fs(p,1)*fs(q,1)+ fs(p,2)*fs(q,2))
                   end if
                end do! lmt2
             end if
          end do! lmt1
       end do! ia
    end if
    call mpi_allreduce(MPI_IN_PLACE,hsr(:,:,:,is),natm*nlmt*nlmt,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)

!!$    if(npes >= 1) then
!!$       allocate(hsr_mpi(natm,nlmt,nlmt))
!!$       call mpi_allreduce(hsr(1,1,1,is),hsr_mpi,natm*nlmt*nlmt,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr) ! MPI
!!$       hsr(:,:,:,is) = hsr_mpi
!!$       deallocate(hsr_mpi)
!!$    end if

    deallocate(fs)
  end subroutine summation_of_ff_sub


  subroutine m_CD_map_chgq_to_fft_box(is,nfftp,bfft)
    integer, intent(in) :: is,nfftp
    real(kind=DP), intent(out) ,dimension(nfftp) :: bfft

    real(kind=DP), allocatable, dimension(:) :: bfft_mpi

    integer :: i, j, ip


    bfft = 0.d0
    do j = 1, kimg
#ifdef VPP
*vocl loop, vector, novrec(bfft)
#endif
       do i = ista_kngp, iend_kngp
          if(kgp_reduced < i) cycle
#ifdef _MPIFFT_
          ip = (igfp_nonpara(i)-1)*kimg + j
#else
          ip = (igfp_l(i)-1)*kimg + j
#endif
          bfft(ip) = bfft(ip) + chgq_l(i,j,is)
       end do
    end do
    if(npes >= 2) then
       allocate(bfft_mpi(nfftp)); bfft_mpi = 0.d0
       call mpi_allreduce(bfft,bfft_mpi,nfftp,mpi_double_precision,mpi_sum,mpi_chg_world,ierr)
       bfft = bfft_mpi
       deallocate(bfft_mpi)
    end if
  end subroutine m_CD_map_chgq_to_fft_box

  subroutine m_CD_map_chgqenl_to_fft_box(is,nfftp,bfft)
    integer, intent(in) :: is,nfftp
    real(kind=DP), intent(out) ,dimension(nfftp) :: bfft

    integer :: i, j, ip, ip_min, ip_max

    ! igfp_enl check

    ip_min = igfp_enl(1)
    ip_max = igfp_enl(1)
    do i = 2, kgp
       if(ip_min > igfp_enl(i)) ip_min = igfp_enl(i)
       if(ip_max < igfp_enl(i)) ip_max = igfp_enl(i)
    end do
    if(ip_min < 0) then
       write(6,'(" !! ip_min = ",i20," <<m_CD_map_chgqnl_to_fft_box>>")') ip_min
       call phase_error_with_msg(6,'<<m_CD_map_chgqenl_to_fft_box>>',__LINE__,__FILE__)
    end if
    if(ip_max*kimg > nfftp) then
       write(6,'(" !! ip_max * kimg = ",i20," > nfftp (= ", i20," ) <<m_CD_map_chgqnl_to_fft_box>>")') ip_max*kimg, nfftp
       call phase_error_with_msg(6,'<<m_CD_map_chgqenl_to_fft_box>>',__LINE__,__FILE__)
    end if


    bfft = 0.d0
    do j = 1, kimg
#ifdef VPP
*vocl loop, vector, novrec(bfft)
#endif
       do i = 1, kgp
          ip = (igfp_enl(i)-1)*kimg + j
          bfft(ip) = bfft(ip) + chgq_enl(i,j)
       end do
    end do
  end subroutine m_CD_map_chgqenl_to_fft_box

  subroutine m_CD_map_fft_box_to_chgqenl(nfftp,bfft)
    integer, intent(in) :: nfftp
    real(kind=DP), intent(in) ,dimension(nfftp) :: bfft
    integer :: i,j,ip
    real(kind=DP) :: rinplw
    rinplw = 1.d0/product(fft_box_size_CD(1:3,1))
    do j = 1, kimg
       do i = 1, kgp
          ip = (igfp_enl(i)-1)*kimg + j
          chgq_enl(i,j) = bfft(ip)*rinplw
       end do
    end do
  end subroutine m_CD_map_fft_box_to_chgqenl

! ================================ added by K. Tagami ================ 11.0
  subroutine m_CD_map_chgqenl_to_fft_box_kt( is, nfftp, bfft, chgq_enl_kt )
    integer, intent(in) :: is,nfftp
    real(kind=DP), intent(in) :: chgq_enl_kt( kgp, kimg, ndim_magmom )

    real(kind=DP), intent(out) ,dimension(nfftp) :: bfft

    integer :: i, j, ip, ip_min, ip_max

    ! igfp_enl check

    ip_min = igfp_enl(1)
    ip_max = igfp_enl(1)
    do i = 2, kgp
       if(ip_min > igfp_enl(i)) ip_min = igfp_enl(i)
       if(ip_max < igfp_enl(i)) ip_max = igfp_enl(i)
    end do
    if(ip_min < 0) then
       write(6,'(" !! ip_min = ",i20," <<m_CD_map_chgqnl_to_fft_box>>")') ip_min
       call phase_error_with_msg(6,'<<m_CD_map_chgqenl_to_fft_box>>',__LINE__,__FILE__)
    end if
    if(ip_max*kimg > nfftp) then
       write(6,'(" !! ip_max * kimg = ",i20," > nfftp (= ", i20," ) <<m_CD_map_chgqnl_to_fft_box>>")') ip_max*kimg, nfftp
       call phase_error_with_msg(6,'<<m_CD_map_chgqenl_to_fft_box>>',__LINE__,__FILE__)
    end if


    bfft = 0.d0
    do j = 1, kimg
#ifdef VPP
*vocl loop, vector, novrec(bfft)
#endif
       do i = 1, kgp
          ip = (igfp_enl(i)-1)*kimg + j
          bfft(ip) = bfft(ip) + chgq_enl_kt(i,j,is)
       end do
    end do
  end subroutine m_CD_map_chgqenl_to_fft_box_kt
! =================================================================== 11.0

!
  subroutine m_CD_keep_retrieve_hsr(keep)
    logical, intent(in) :: keep
    real(kind=DP),allocatable,dimension(:,:,:,:),save   :: hsr_tmp
    real(kind=DP),allocatable,dimension(:,:,:,:),save   :: hsi_tmp
    if(keep)then
      if ( noncol ) then
        if (.not.allocated(hsr_tmp)) allocate(hsr_tmp(natm,nlmt,nlmt,ndim_magmom)); hsr_tmp = hsr
        if (.not.allocated(hsi_tmp)) allocate(hsi_tmp(natm,nlmt,nlmt,ndim_magmom))
! === ASMS 2017/09/04
        hsi_tmp = hsi
! === ASMS 2017/09/04
      else
        if(.not.allocated(hsr_tmp)) allocate(hsr_tmp(natm,nlmt,nlmt,nspin)); hsr_tmp = hsr
      endif
    else
      hsr = hsr_tmp
      deallocate(hsr_tmp)
      if ( noncol ) then
        hsi = hsi_tmp
        deallocate(hsi_tmp)
      endif
    endif
  end subroutine m_CD_keep_retrieve_hsr

  subroutine m_CD_keep_chgq_l()
! =============================== modified y K. Tagami =============== 11.0
!!    allocate(chgq_tmp(ista_kngp:iend_kngp,kimg,nspin)); chgq_tmp = 0.d0
!
    if ( noncol ) then
       allocate(chgq_tmp(ista_kngp:iend_kngp,kimg,ndim_magmom));
    else
       allocate(chgq_tmp(ista_kngp:iend_kngp,kimg,nspin));
    endif
    chgq_tmp = 0.d0
! =================================================================== 11.0
    chgq_tmp = chgq_l
  end subroutine m_CD_keep_chgq_l

  subroutine m_CD_retrieve_chgq()
    chgq_l = chgq_tmp
    deallocate(chgq_tmp)
  end subroutine m_CD_retrieve_chgq

  subroutine m_CD_rspace_put_headermark(nfchr,nspin,iloop,i)
    integer, intent(in) :: nfchr, nspin, iloop, i

    character(len=80) :: line
    logical::            tf

    if(mype==0) then
       backspace nfchr
       read(nfchr,'(a80)',end=1001,err=1001) line
       call strncmp0(trim(line),'END',tf)
       if(tf) goto 1002
1001    write(nfchr,'(a3)') 'END'
1002    continue

! ========================================= modified by K. Tagami ====== 11.0
!!    if(nspin == 1) then
!!!!$       write(nfchr,'(" ---- partial charge ---- ", i5," -th energy-windows")') i
!!       write(nfchr,'("PARTIALCHARGE  energy_window = ",i5)') i
!!    else if(nspin == 2) then
!!!!$       write(nfchr,'(" ---- partial charge ---- ", i5 &
!!!!$            &       ," -th energy-windows, spin-state = ",i3)') i, iloop
!!       write(nfchr,'("PARTIALCHARGE  energy_windows = ",i5," spin-state = ",i3)') i, iloop
!!    end if

       if ( noncol ) then
          write(nfchr,'("PARTIALCHARGE  energy_windows = ",i5," mag-direction = ",i3)') &
                                    & i, iloop
       else
         if(nspin == 1) then
            write(nfchr,'("PARTIALCHARGE  energy_window = ",i5)') i
         else if(nspin == 2) then
! ============================= KT_mod === 2014/06/07
!         write(nfchr,'("PARTIALCHARGE  energy_windows = ",i5," spin-state = ",i3)') &
!                   &          i, iloop
!
            if ( iloop > 0 ) then
               write(nfchr,'("PARTIALCHARGE  energy_windows = ",i5," spin-state = ",i3)') &
                    &          i, iloop
            else
               if ( iloop == -1 ) then
                  write(nfchr,'("PARTIALCHARGE  energy_windows = ",i5," TOTAL")') i
               else if ( iloop == -2 ) then
                  write(nfchr,'("PARTIALCHARGE  energy_windows = ",i5," MAGMOM")') i
               endif
            endif
! ======================================= 2014/06/07
         endif
       end if
    endif
! ======================================================================= 11.0
  end subroutine m_CD_rspace_put_headermark

  subroutine m_CD_rspace_put_endmark(nfchr)
    integer, intent(in) :: nfchr
    if(mype==0) then
       write(nfchr,'(a3)') 'END'
    endif
  end subroutine m_CD_rspace_put_endmark


  subroutine m_CD_den_mat(nfout,kv3)
    integer, intent(in) :: nfout,kv3
    !! debug
                                                 __TIMER_SUB_START(736)
    write(*,*) 'calling summation_of_ff_3D'
    call summation_of_ff_3D(kv3,0)
                                                 __TIMER_SUB_STOP(736)
  end subroutine m_CD_den_mat

  subroutine m_CD_cp_chgq_prev_to_chgq(nfout)
    integer, intent(in) :: nfout
    if(mype == 0) then
      write(nfout,'(a)') ' !** copied previous charge density to current charge density'
    endif
    chgq_l = chgq_l_prev*(univol_prev/univol)
  end subroutine m_CD_cp_chgq_prev_to_chgq

#ifdef FFTW3
  subroutine m_CD_gen_chgq_from_prev_chgq(nfout)
    use m_Files, only : m_Files_nfcntn_bin_paw_exists,m_Files_open_nfcntn_bin_paw &
                      , m_Files_close_nfcntn_bin_paw, nfcntn_bin_paw
    integer, intent(in) :: nfout
    integer :: i,j,k, idp, nlp, nmp, nnp, nlphf,inew,jnew,knew,ip,mmp,is,nlp2,nmp2,nnp2
    integer :: ii,jj,kk
    integer, parameter :: UP = 1 , DOWN = 2
    integer ::            up_down
    integer :: nfftp_nonpara_prev
    real(kind=DP) ::      s1, s2, sratio,x,y,z
    real(kind=DP), allocatable, dimension(:) :: xmat,ymat,zmat
    real(kind=DP), allocatable, dimension(:,:,:,:) :: cmat, cmat2
    real(kind=DP) :: trilinear_interpolation,tricubic_spline
    real(kind=DP) :: total_charge,dy1,dy2,dy3
    integer(kind=DP) :: plan
    integer :: FFTW_MEASURE = 0
    integer :: id_sname = -1
    call tstatc0_begin('m_CD_gen_chgq_from_prev_chgq ',id_sname,1)
    idp = fft_box_size_CD_nonpara_prev(1,0)
    mmp = fft_box_size_CD_nonpara_prev(2,0)
    nfftp_nonpara_prev  = product(fft_box_size_CD_nonpara_prev(1:3,0)) * kimg
    nlp = fft_box_size_CD_prev(1,1)
    nmp = fft_box_size_CD_prev(2,1)
    nnp = fft_box_size_CD_prev(3,1)
    nlp2 = nlp+2;nmp2 = nmp+2;nnp2 = nnp+2
    allocate(xmat(nlp2));allocate(ymat(nmp2));allocate(zmat(nnp2))
    allocate(cmat(nlp2,nmp2,nnp2,nspin))
    if(interpolation_method_chg == SPLINE_INTERPOLATION) allocate(cmat2(nlp2,nmp2,nnp2,nspin))
    if(kimg == 1) then
       nlphf = idp/2
    else
       nlphf = idp
    end if
    call alloc_rspace_charge_prev()
    if(kimg == 1)then
      call dfftw_plan_dft_r2c_3d(plan,nlp,nmp,nnp,afft(1),afft(1),FFTW_MEASURE)
    else
      call dfftw_plan_dft_3d(plan,nlp,nmp,nnp,afft(1),afft(1),+1,FFTW_MEASURE)
    endif
    call dealloc_rspace_charge_prev()

    do is=1,nspin
      call alloc_rspace_charge_prev()
      call map_valence_charge_to_fft_box_prev(is,.false.)
      if(kimg==1) then
         call dfftw_execute_dft_r2c(plan,afft,afft)
      else
         call dfftw_execute_dft(plan,afft,afft)
      endif
      s1 = 0.d0; s2 = 0.d0
      do ip = 1, nlphf*mmp*nnp,2
         s1 = s1 + dabs(afft(ip))
         s2 = s2 + dabs(afft(ip+1))
      end do
      if(s1 > s2) then
         up_down = UP
         sratio = s2/s1
      else
         up_down = DOWN
         sratio = s1/s2
      endif
      do i = 0, nmp+1
         do j = 0, nnp+1
            do k = 0, nlp+1
               ii=i;jj=j;kk=k
               if(i.eq.0) ii = nmp
               if(j.eq.0) jj = nnp
               if(k.eq.0) kk = nlp
               if(i.eq.nmp+1) ii = 1
               if(j.eq.nnp+1) jj = 1
               if(k.eq.nlp+1) kk = 1
               if(kimg == 1 .and. kk > nlphf) then
                  knew = idp - kk
                  jnew = nnp+2 - jj
                  inew = nmp+2 - ii
                  if(jnew > nnp) then
                     jnew = jnew - nnp
                  end if
                  if(inew > nmp) then
                     inew = inew - nmp
                  end if
               else
                  knew = kk; jnew = jj; inew = ii
               end if
               ip = nlphf*mmp*(jnew-1) + nlphf*(inew-1) + knew
               xmat(k+1) = dble(k)/dble(nlp)
               ymat(i+1) = dble(i)/dble(nmp)
               zmat(j+1) = dble(j)/dble(nnp)
               cmat(k+1,i+1,j+1,is) = afft(ip*2-2+up_down)
            enddo
         enddo
      enddo
      call dealloc_rspace_charge_prev()
      if(interpolation_method_chg== SPLINE_INTERPOLATION) &
   &   call init_tricubic_spline(xmat,ymat,zmat,cmat(:,:,:,is),nlp2,nmp2,nnp2,cmat2(:,:,:,is))
    enddo

    idp = fft_box_size_CD_nonpara(1,0)
    mmp = fft_box_size_CD_nonpara(2,0)
    nlp = fft_box_size_CD(1,1)
    nmp = fft_box_size_CD(2,1)
    nnp = fft_box_size_CD(3,1)
    if(kimg == 1) then
       nlphf = idp/2
    else
       nlphf = idp
    end if
    do is=1,nspin
      call m_CD_alloc_rspace_charge()
      do i = 1, nmp
         do j = 1, nnp
            do k = 1, nlp
               if(kimg == 1 .and. k > nlphf) then
                  knew = idp - k
                  jnew = nnp+2 - j
                  inew = nmp+2 - i
                  if(jnew > nnp) then
                     jnew = jnew - nnp
                  end if
                  if(inew > nmp) then
                     inew = inew - nmp
                  end if
               else
                  knew = k; jnew = j; inew = i
               end if
               ip = nlphf*mmp*(jnew-1) + nlphf*(inew-1) + knew
               if(interpolation_method_chg == SPLINE_INTERPOLATION) then
                 afft(ip*2-2+up_down) = tricubic_spline( &
               & xmat,ymat,zmat,cmat(:,:,:,is),cmat2(:,:,:,is),nlp2,nmp2,nnp2, &
               & dble(k)/dble(nlp),dble(i)/dble(nmp),dble(j)/dble(nnp), &
               & dy1,dy2,dy3)
               else
                 afft(ip*2-2+up_down) = trilinear_interpolation( &
               & nlp2,nmp2,nnp2,cmat(:,:,:,is),xmat,ymat,zmat, &
               & dble(k)/dble(nlp),dble(i)/dble(nmp),dble(j)/dble(nnp))
               endif
            enddo
         enddo
      enddo
      afft = afft * (univol_prev/univol)
      call m_FFT_CD0(nfout,afft,DIRECT)
      call cpafft_CD_to_valencecharge(is) ! afft -> chgq_l
      call m_CD_dealloc_rspace_charge()
    enddo

    total_charge = 0.d0
    if(mype == 0) then
       total_charge = chgq_l(1,1,1)*univol
       if(nspin>1) total_charge = total_charge + chgq_l(1,1,2)*univol
    endif
    if(npes > 1) call mpi_bcast(total_charge,1,mpi_double_precision,0,MPI_CommGroup,ierr)

    chgq_l(:,:,:) = chgq_l(:,:,:)*totch/total_charge

    if ( flg_paw .and. read_charge_hardpart == YES .and.  m_Files_nfcntn_bin_paw_exists()) then
       call m_Files_open_nfcntn_bin_paw()
       call m_PP_rd_PAW_parameters(nfout,nfcntn_bin_paw)
       call m_CD_rd_hsr(nfcntn_bin_paw)
       call m_Files_close_nfcntn_bin_paw()
    endif

    deallocate(xmat)
    deallocate(ymat)
    deallocate(zmat)
    deallocate(cmat)
    if(interpolation_method_chg==SPLINE_INTERPOLATION) deallocate(cmat2)
    call dfftw_destroy_plan(plan)
    call tstatc0_end(id_sname)

    contains

    subroutine map_valence_charge_to_fft_box_prev(iloop,flg_add_core)
      integer, intent(in) :: iloop
      logical, intent(in) :: flg_add_core
      integer :: j, i, ip, it

      afft_mpi1 = 0.d0
      do j = 1, kimg
         do i = ista_kngp_prev, iend_kngp_prev
            if (kgp_reduced_prev < i) cycle
            ip = (igfp_l_prev(i)-1)*kimg + j
! === KT_mod === 2014/06/07
!            afft_mpi1(ip) = afft_mpi1(ip) + chgq_l(i,j,iloop)
!
            if ( iloop > 0 ) then
               afft_mpi1(ip) = afft_mpi1(ip) + chgq_l_prev(i,j,iloop)
            else
               if ( iloop == -1 ) then
                  afft_mpi1(ip) = afft_mpi1(ip) + chgq_l_prev(i,j,1) +chgq_l_prev(i,j,2)
               else if ( iloop == -2 ) then
                  afft_mpi1(ip) = afft_mpi1(ip) + chgq_l_prev(i,j,1) -chgq_l_prev(i,j,2)
               endif
            endif
! ============== 2014/06/07
         end do
      end do


      if(npes >= 2) then
         call mpi_allreduce(afft_mpi1,afft,nfftp_nonpara_prev,mpi_double_precision &
              &  ,mpi_sum,mpi_chg_world,ierr)
      else
         afft = afft_mpi1
      end if
    end subroutine map_valence_charge_to_fft_box_prev

    subroutine alloc_rspace_charge_prev()
      allocate(afft(nfftp_nonpara_prev)); afft = 0.d0
      allocate(afft_mpi1(nfftp_nonpara_prev))
! ======================================== Added by K. Tagami =======
      afft = 0.0d0; afft_mpi1 = 0.0d0
    end subroutine alloc_rspace_charge_prev

    subroutine dealloc_rspace_charge_prev()
      deallocate(afft); deallocate(afft_mpi1)
    end subroutine dealloc_rspace_charge_prev
  end subroutine m_CD_gen_chgq_from_prev_chgq
#endif

  subroutine m_CD_dealloc
    use m_Files, only : nfout
    logical :: unitcell_can_change
    integer :: isp

    if(unitcell_can_change() .and. sw_interpolate_charge == ON)then
      if(allocated(chgq_l_prev)) deallocate(chgq_l_prev)
      allocate(chgq_l_prev(ista_kngp:iend_kngp,kimg,ndim_magmom))
      chgq_l_prev = chgq_l
      univol_prev = univol
    endif
    if(allocated(chgq_l)) deallocate(chgq_l)
    if(allocated(chgqo_l)) deallocate(chgqo_l)
    if (allocated(chgsoft)) deallocate(chgsoft)

    if(allocated(hsr)) deallocate(hsr)

! ==================== added by K. Tagami ======================== 5.0&11.0
    if (allocated(hsi)) deallocate(hsi)
    if (allocated(hsio)) deallocate(hsio)

    if (allocated(hsro)) deallocate(hsro)
! ================================================================ 5.0&11.0

! === KT_add === 2014/06/10
    if ( allocated(extpl_target_history) ) deallocate(extpl_target_history)
! ============== 2014/06/10

  end subroutine m_CD_dealloc

  subroutine map_fft_to_chgq_3D(lsize, ibesize, inn_fft_l, out_fft_l, nfout)
    use m_Parallelization,     only : fft_chgq_scnt, fft_chgq_rcnt &
   &                                , fft_chgq_send, fft_chgq_recv &
   &                                , fft_chgq_index, fft_chgq_dist &
   &                                , fft_chgq_maxrecv, fft_chgq_maxsend
    integer, intent(in)  :: lsize, ibesize, nfout
#ifdef FFT_3D_DIVISION
    real(kind=DP), dimension(lsize*2   ,ibesize), intent(in)    :: inn_fft_l
#else
    real(kind=DP), dimension(lsize*kimg,ibesize), intent(in)    :: inn_fft_l
#endif
    real(kind=DP), dimension(np_kngp_gw*kimg,ibesize), intent(out) :: out_fft_l
    integer, dimension(0:nrank_g-1)                       ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)       ::sta_r, sta_s
    integer, parameter :: itag = 11

    real(kind=DP), allocatable, dimension(:,:) :: sendbuf, recvbuf
    integer :: icnt_send, icnt_recv, ierr, lrank, i, j, k, iadd
! === DEBUG by tkato 2012/06/04 ================================================
#ifndef USE_NONBLK_COMM
    integer, allocatable, dimension(:) :: sdsp, rdsp
#endif
! ==============================================================================
                                                 __TIMER_SUB_START(718)
    ierr = 0
    if (fft_chgq_maxsend /= 0) then
       allocate(sendbuf(fft_chgq_maxsend*kimg*ibesize,0:nrank_g-1), stat=ierr)
       sendbuf = 0.0d0
! === DEBUG by tkato 2012/06/05 ================================================
    else
       allocate(sendbuf(1,1), stat=ierr)
! ==============================================================================
    endif
    if (fft_chgq_maxrecv /= 0) then
       allocate(recvbuf(fft_chgq_maxrecv*kimg*ibesize,0:nrank_g-1), stat=ierr)
       recvbuf = 0.0d0
! === DEBUG by tkato 2012/06/05 ================================================
    else
       allocate(recvbuf(1,1), stat=ierr)
! ==============================================================================
    endif
     if (ierr /= 0) then
        write(nfout,*)' map_info_fft_to_WF_3D :  Not allocate '
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 261, ierr)
     endif

#ifdef USE_NONBLK_COMM
                                                 __TIMER_COMM_START_w_BARRIER(mpi_kg_world,827)
    if (fft_chgq_maxrecv /= 0) then
       icnt_recv = 0
       lrank = mod(myrank_g,nrank_g)
       do i = 0, nrank_g - 1
          lrank = lrank + 1
          if (lrank > (nrank_g-1)) lrank = 0
          if (fft_chgq_rcnt(lrank) /= 0) then
             call mpi_irecv(recvbuf(1,lrank), fft_chgq_rcnt(lrank)*kimg*ibesize, mpi_double_precision,&
            &               lrank, itag, mpi_ske_world, req_r(icnt_recv), ierr)
              if (ierr /= 0) then
                 write(nfout,*)' map_info_fft_to_WF_3D :  mpi_irecv error'
                 call flush(nfout)
                 call mpi_abort(mpi_comm_world, 262, ierr)
              endif
             icnt_recv = icnt_recv + 1
          endif
       enddo
    endif
#endif

                                                 __TIMER_DO_START(828)
    if (fft_chgq_maxsend /= 0) then
       if (kimg == 1) then
#ifdef FFT_3D_DIVISION
       integer :: i1, lx, ly, lz, mx, my, mz, mm, kx1p, kx2p, kx3p, jadd
          lx = fft_box_size_WF(1,0)
          ly = fft_box_size_WF(2,0)
          lz = fft_box_size_WF(3,0)
          kx1p = fft_X_x_nel
          kx2p = fft_X_y_nel
          kx3p = fft_X_z_nel
          do k = 1, nel_fft_x(myrank_g)
             if(fft_chgq_index(k) == 0) cycle

             i1 = mp_fft_x(k)
             mz = (i1-1)/(lx*ly)+1
             mm = mod(i1,(lx*ly))
             if (mm==0) mm=lx*ly
             my = (mm-1)/lx+1
             mx = mod(mm,lx)
             if (mx==0) mx = lx
             jadd = mx-xyz_fft_x(1,1)+1+kx1p*(my-xyz_fft_x(1,2))+kx1p*kx2p*(mz-xyz_fft_x(1,3))

             do i = 1, ibesize
                sendbuf(ibesize*(fft_chgq_index(k)-1)+i,fft_chgq_dist(k)) = inn_fft_l(jadd*2-1,i)
             enddo
          end do
#else
          do k = 1, nel_fft_x(myrank_g)
             if(fft_chgq_index(k) == 0) cycle
             do i = 1, ibesize
                sendbuf(ibesize*(fft_chgq_index(k)-1)+i,fft_chgq_dist(k)) = inn_fft_l(k,i)
             enddo
          end do
#endif
       else
#ifdef FFT_3D_DIVISION
!      integer :: i1, lx, ly, lz, mx, my, mz, mm, kx1p, kx2p, kx3p, jadd
          lx = fft_box_size_WF(1,0)
          ly = fft_box_size_WF(2,0)
          lz = fft_box_size_WF(3,0)
          kx1p = fft_X_x_nel
          kx2p = fft_X_y_nel
          kx3p = fft_X_z_nel
          do k = 1, nel_fft_x(myrank_g)
             if(fft_chgq_index(k) == 0) cycle

             i1 = mp_fft_x(k)
             mz = (i1-1)/(lx*ly)+1
             mm = mod(i1,(lx*ly))
             if (mm==0) mm=lx*ly
             my = (mm-1)/lx+1
             mx = mod(mm,lx)
             if (mx==0) mx = lx
             jadd = mx-xyz_fft_x(1,1)+1+kx1p*(my-xyz_fft_x(1,2))+kx1p*kx2p*(mz-xyz_fft_x(1,3))

             do i = 1, ibesize
                iadd = ibesize*2*(fft_chgq_index(k)-1)+i*2
                sendbuf(iadd-1,fft_chgq_dist(k)) = inn_fft_l(jadd*2-1,i)
                sendbuf(iadd,  fft_chgq_dist(k)) = inn_fft_l(jadd*2  ,i)
             enddo
          end do
#else
          do k = 1, nel_fft_x(myrank_g)
             if(fft_chgq_index(k) == 0) cycle
             do i = 1, ibesize
                iadd = ibesize*2*(fft_chgq_index(k)-1)+i*2
                sendbuf(iadd-1,fft_chgq_dist(k)) = inn_fft_l(k*2-1,i)
                sendbuf(iadd,  fft_chgq_dist(k)) = inn_fft_l(k*2  ,i)
             enddo
          end do
#endif
       end if
    end if
                                                 __TIMER_DO_STOP(828)
#ifdef USE_NONBLK_COMM
    if (fft_chgq_maxsend /= 0) then
       icnt_send = 0
       lrank = mod((myrank_g+1),nrank_g)
       do i = 0, nrank_g - 1
          lrank = lrank + 1
          if (lrank > (nrank_g - 1)) lrank = 0
          if (fft_chgq_scnt(lrank) /= 0) then
             call mpi_isend(sendbuf(1,lrank), fft_chgq_scnt(lrank)*kimg*ibesize, mpi_double_precision,&
            &               lrank, itag, mpi_ske_world, req_s(icnt_send), ierr)
              if (ierr /= 0) then
                 write(nfout,*)' map_info_fft_to_WF_3D :  mpi_isend error'
                 call flush(nfout)
                 call mpi_abort(mpi_comm_world, 263, ierr)
              endif
             icnt_send = icnt_send + 1
          endif
       enddo
    endif

    if (fft_chgq_maxrecv /= 0) then
       call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
        if (ierr /= 0) then
           write(nfout,*)' map_info_fft_to_WF_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 264, ierr)
        endif
    endif

    if (fft_chgq_maxsend /= 0) then
       call mpi_waitall(icnt_send, req_s, sta_s, ierr)
        if (ierr /= 0) then
           write(nfout,*)' map_info_fft_to_WF_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 265, ierr)
        endif
    endif
                                                 __TIMER_COMM_STOP(827)
#else
                                                 __TIMER_COMM_START_w_BARRIER(mpi_kg_world,693)
! === DEBUG by tkato 2012/06/04 ================================================
!      integer, allocatable, dimension(:) :: sdsp, rdsp
! ==============================================================================
       allocate(sdsp(0:nrank_g-1), stat=ierr)
       allocate(rdsp(0:nrank_g-1), stat=ierr)
       do i = 0, nrank_g - 1
          sdsp(i)=fft_chgq_maxsend*kimg*ibesize*i
          rdsp(i)=fft_chgq_maxrecv*kimg*ibesize*i
       enddo
       call MPI_ALLTOALLV(      sendbuf, fft_chgq_scnt*kimg*ibesize, sdsp, &
      &   mpi_double_precision, recvbuf, fft_chgq_rcnt*kimg*ibesize, rdsp, &
      &   mpi_double_precision, mpi_ske_world, ierr )
       if (ierr /= 0) then
        write(nfout,*)' map_info_fft_to_WF_3D : mpi_alltoallv error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 266, ierr)
       endif
       deallocate(sdsp)
       deallocate(rdsp)
                                                 __TIMER_COMM_STOP(693)
#endif

    out_fft_l = 0.0d0
                                                 __TIMER_DO_START(829)
    if (kimg == 1) then
       do i = 0, nrank_g - 1
          if (fft_chgq_rcnt(i) /= 0) then
             do k = 1, fft_chgq_rcnt(i)
                do j = 1, ibesize
                   out_fft_l(fft_chgq_recv(k,i),j) = recvbuf(ibesize*(k-1)+j,i)
                enddo
             end do
          end if
       end do
    else
       do i = 0, nrank_g - 1
          if (fft_chgq_rcnt(i) /= 0) then
             do k = 1, fft_chgq_rcnt(i)
                do j = 1, ibesize
                   iadd = ibesize*2*(k-1)+j*2
                   out_fft_l(fft_chgq_recv(k,i)*2-1,j) = recvbuf(iadd-1,i)
                   out_fft_l(fft_chgq_recv(k,i)*2  ,j) = recvbuf(iadd,  i)
                enddo
             end do
          end if
       end do
    end if
                                                 __TIMER_DO_STOP(829)
    if (allocated(sendbuf)) deallocate(sendbuf)
    if (allocated(recvbuf)) deallocate(recvbuf)
                                                 __TIMER_SUB_STOP(718)
  end subroutine map_fft_to_chgq_3D

#ifdef MPI_FFTW
  subroutine map_fft_to_chgq_mpifftw(lsize, ibesize, out_fft_l, nfout)
    use m_Parallelization,     only : fft_chgq_scnt_mpifftw, fft_chgq_rcnt_mpifftw &
   &                                , fft_chgq_send_mpifftw, fft_chgq_recv_mpifftw &
   &                                , fft_chgq_index_mpifftw, fft_chgq_dist_mpifftw &
   &                                , fft_chgq_maxrecv_mpifftw, fft_chgq_maxsend_mpifftw
    integer, intent(in)  :: lsize, ibesize, nfout
    real(kind=DP), dimension(np_kngp_gw*kimg,ibesize), intent(out) :: out_fft_l
    integer, dimension(0:nrank_g-1)                       ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)       ::sta_r, sta_s
    integer, parameter :: itag = 11

    real(kind=DP), allocatable, dimension(:,:) :: sendbuf, recvbuf
    integer :: icnt_send, icnt_recv, ierr, lrank, i, j, k, iadd, mm, i1
    integer(C_INTPTR_T)  :: local_n, local_n_offset, alloc_local, lx, ly, lz, mx,my,mz
    lx = fft_box_size_WF(1,0)
    ly = fft_box_size_WF(2,0)
    lz = fft_box_size_WF(3,0)
    alloc_local = fftw_mpi_local_size_3d(lz,ly,lx,mpi_ke_world,local_n,local_n_offset)

    ierr = 0
    if (fft_chgq_maxsend_mpifftw /= 0) then
       allocate(sendbuf(fft_chgq_maxsend_mpifftw*kimg*ibesize,0:nrank_g-1), stat=ierr)
       sendbuf = 0.0d0
! === DEBUG by tkato 2012/06/05 ================================================
    else
       allocate(sendbuf(1,1), stat=ierr)
! ==============================================================================
    endif
    if (fft_chgq_maxrecv_mpifftw /= 0) then
       allocate(recvbuf(fft_chgq_maxrecv_mpifftw*kimg*ibesize,0:nrank_g-1), stat=ierr)
       recvbuf = 0.0d0
! === DEBUG by tkato 2012/06/05 ================================================
    else
       allocate(recvbuf(1,1), stat=ierr)
! ==============================================================================
    endif
     if (ierr /= 0) then
        write(nfout,*)' map_info_fft_to_WF_3D :  Not allocate '
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 261, ierr)
     endif

    if (fft_chgq_maxrecv_mpifftw /= 0) then
       icnt_recv = 0
       lrank = mod(myrank_g,nrank_g)
       do i = 0, nrank_g - 1
          lrank = lrank + 1
          if (lrank > (nrank_g-1)) lrank = 0
          if (fft_chgq_rcnt_mpifftw(lrank) /= 0) then
             call mpi_irecv(recvbuf(1,lrank), fft_chgq_rcnt_mpifftw(lrank)*kimg*ibesize, mpi_double_precision,&
            &               lrank, itag, mpi_ke_world, req_r(icnt_recv), ierr)
              if (ierr /= 0) then
                 write(nfout,*)' map_info_fft_to_WF_3D :  mpi_irecv error'
                 call flush(nfout)
                 call mpi_abort(mpi_comm_world, 262, ierr)
              endif
             icnt_recv = icnt_recv + 1
          endif
       enddo
    endif

    if (fft_chgq_maxsend_mpifftw /= 0) then
       if (kimg == 1) then
          do k = 1, local_n*lx*ly
             if(fft_chgq_index_mpifftw(k) == 0) cycle
             i1 = k
             mz = (i1-1)/(lx*ly)+1
             mm = mod(i1,lx*ly)
             if (mm==0) then
                my = ly
             else
                my = (mm-1)/lx+1
             end if
             mx = mod(mm,lx)
             if(mx==0) mx = lx
             do i = 1, ibesize
                sendbuf(ibesize*(fft_chgq_index_mpifftw(k)-1)+i,fft_chgq_dist_mpifftw(k)) = bfft_mpifftw_kimg1(mx,my,mz)
             enddo
          end do
       else
          do k = 1, local_n*lx*ly
             if(fft_chgq_index_mpifftw(k) == 0) cycle
             i1 = k
             mz = (i1-1)/(lx*ly)+1
             mm = mod(i1,(lx*ly))
             if (mm==0) mm=lx*ly
             my = (mm-1)/lx+1
             mx = mod(mm,lx)
             if (mx==0) mx = lx
             do i = 1, ibesize
                iadd = ibesize*2*(fft_chgq_index_mpifftw(k)-1)+i*2
                sendbuf(iadd-1,fft_chgq_dist_mpifftw(k)) = real (bfft_mpifftw(mx,my,mz))
                sendbuf(iadd,  fft_chgq_dist_mpifftw(k)) = aimag(bfft_mpifftw(mx,my,mz))
             enddo
          end do
       end if
    end if

    if (fft_chgq_maxsend_mpifftw /= 0) then
       icnt_send = 0
       lrank = mod((myrank_g+1),nrank_g)
       do i = 0, nrank_g - 1
          lrank = lrank + 1
          if (lrank > (nrank_g - 1)) lrank = 0
          if (fft_chgq_scnt_mpifftw(lrank) /= 0) then
             call mpi_isend(sendbuf(1,lrank), fft_chgq_scnt_mpifftw(lrank)*kimg*ibesize, mpi_double_precision,&
            &               lrank, itag, mpi_ke_world, req_s(icnt_send), ierr)
              if (ierr /= 0) then
                 write(nfout,*)' map_info_fft_to_WF_3D :  mpi_isend error'
                 call flush(nfout)
                 call mpi_abort(mpi_comm_world, 263, ierr)
              endif
             icnt_send = icnt_send + 1
          endif
       enddo
    endif

    if (fft_chgq_maxrecv_mpifftw /= 0) then
       call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
        if (ierr /= 0) then
           write(nfout,*)' map_info_fft_to_WF_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 264, ierr)
        endif
    endif

    if (fft_chgq_maxsend_mpifftw /= 0) then
       call mpi_waitall(icnt_send, req_s, sta_s, ierr)
        if (ierr /= 0) then
           write(nfout,*)' map_info_fft_to_WF_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 265, ierr)
        endif
    endif

    out_fft_l = 0.0d0
    if (kimg == 1) then
       do i = 0, nrank_g - 1
          if (fft_chgq_rcnt_mpifftw(i) /= 0) then
             do k = 1, fft_chgq_rcnt_mpifftw(i)
                do j = 1, ibesize
                   out_fft_l(fft_chgq_recv_mpifftw(k,i),j) = recvbuf(ibesize*(k-1)+j,i)
                enddo
             end do
          end if
       end do
    else
       do i = 0, nrank_g - 1
          if (fft_chgq_rcnt_mpifftw(i) /= 0) then
             do k = 1, fft_chgq_rcnt_mpifftw(i)
                do j = 1, ibesize
                   iadd = ibesize*2*(k-1)+j*2
                   out_fft_l(fft_chgq_recv_mpifftw(k,i)*2-1,j) = recvbuf(iadd-1,i)
                   out_fft_l(fft_chgq_recv_mpifftw(k,i)*2  ,j) = recvbuf(iadd,  i)
                enddo
             end do
          end if
       end do
    end if

    if (allocated(sendbuf)) deallocate(sendbuf)
    if (allocated(recvbuf)) deallocate(recvbuf)

  end subroutine map_fft_to_chgq_mpifftw
#endif

  subroutine m_CD_predictor_pre(nfout,printable)
    use m_Ionic_System,  only : cpo_l

    integer, intent(in) :: nfout
    logical, intent(in) :: printable

    if ( mag_moment0_atoms_is_defined ) then
       call case_pre_by_atoms
    else
       call case_pre_by_atomtype
    endif

    if(sw_extrapolate_charge==ON)then
       extpl_target_history(:,:,:,3) = extpl_target_history(:,:,:,2)
       extpl_target_history(:,:,:,2) = extpl_target_history(:,:,:,1)
       extpl_target_history(:,:,:,1) = chgq_l(:,:,:)
    endif

  contains

    subroutine case_pre_by_atomtype
      integer :: ik,it,is,ig,i
      real(kind=DP) :: fac,pm

      if ( noncol ) then
         Do is=1, ndim_magmom
            do it=1, ntyp
               if ( is == 1 ) then
                  fac = 1.0
               else
                  fac = zeta1(it) *mag_direction0_atomtyp(it,is-1)
               endif
               
               do ik=1,kimg
                  do ig=ista_kngp,iend_kngp
                     chgq_l(ig,ik,is) = chgq_l(ig,ik,is) &
                          &           -fac*rhvg_l(ig,it)*zfm3_l(ig,it,ik)
                  end do
               end do
               if(ista_kngp==1) then
                  chgq_l(1,1,is) = chgq_l(1,1,is) - zfm3_l(1,it,1)*fac*ival(it)/univol
               endif
            end do
         end do
      else
         
         fac=1.d0;   pm=1.d0
         do is=1,nspin
            if(is==2) pm=-1.d0
            do it=1,ntyp
               if(nspin>1) fac = (1.d0+pm*zeta1(it))*0.5d0
               do ik=1,kimg
                  do ig=ista_kngp,iend_kngp
                     chgq_l(ig,ik,is) = chgq_l(ig,ik,is) &
                          &           -fac*rhvg_l(ig,it)*zfm3_l(ig,it,ik)
                  enddo
               enddo
               if(ista_kngp==1) then
                  chgq_l(1,1,is) = chgq_l(1,1,is) - zfm3_l(1,it,1)*fac*ival(it)/univol
               endif
            enddo
         enddo
         
      endif
    end subroutine case_pre_by_atomtype

    subroutine case_pre_by_atoms
      integer :: ia, it, is, i
      real(kind=DP) :: f2, f3, c1, s1
      real(kind=DP) :: grt, nchg(4), ncharge, rltv_t(3,3)
      real(kind=DP), allocatable :: pos_wk(:,:)

      allocate( pos_wk(natm,3) )
      rltv_t = transpose(rltv)/PAI2
      call change_of_coordinate_system( rltv_t, cpo_l(:,:,1), &
           &                            natm, natm, pos_wk )

      if ( noncol ) then
         do ia=1, natm
            it = ityp(ia)
            ncharge = ival(it) +ionic_charge_atoms(ia)
            nchg(1) = ncharge
            nchg(2) = mag_moment0_atoms(ia,1)
            nchg(3) = mag_moment0_atoms(ia,2)
            nchg(4) = mag_moment0_atoms(ia,3)

            Do is=1, ndim_magmom
               f2 = nchg(is) /ival(it)
               do i=ista_kngp,iend_kngp
                  grt = pos_wk(ia,1) *ngabc_kngp_l(i,1) &
                       &   +pos_wk(ia,2) *ngabc_kngp_l(i,2) &
                       &   +pos_wk(ia,3 )*ngabc_kngp_l(i,3)
                  grt = grt *PAI2

                  if ( kimg == 1 ) then
                     c1 = iwei(ia) *cos(grt);
                     chgq_l(i,1,is) = chgq_l(i,1,is) -c1 *f2 *rhvg_l(i,it)
                  else
                     c1 = cos(grt);  s1 = -sin(grt)
                     chgq_l(i,1,is) = chgq_l(i,1,is) -c1 *f2 *rhvg_l(i,it)
                     chgq_l(i,2,is) = chgq_l(i,2,is) -s1 *f2 *rhvg_l(i,it)
                  endif
               end do

               f3 = nchg(is) /univol
               if ( ista_kngp == 1) then
                  c1 = iwei(ia);
                  chgq_l(1,1,is) = chgq_l(1,1,is) -c1 *f3
               endif
            end do
         end do

      else
         do ia = 1, natm
            it = ityp(ia)

            ncharge = ival(it) +ionic_charge_atoms(ia)
            if ( nspin == 1 ) then
               nchg(1) = ncharge
            else
               nchg(1) = ( ncharge +mag_moment0_atoms(ia,1) ) /2.0d0
               nchg(2) = ( ncharge -mag_moment0_atoms(ia,1) ) /2.0d0
            endif

            Do is=1, nspin
               f2 = nchg(is) /ival(it)
               do i = ista_kngp, iend_kngp
                  grt = pos_wk(ia,1) *ngabc_kngp_l(i,1) &
                       &   +pos_wk(ia,2) *ngabc_kngp_l(i,2) &
                       &   +pos_wk(ia,3 )*ngabc_kngp_l(i,3)
                  grt = grt *PAI2
                  if ( kimg == 1 ) then
                     c1 = iwei(ia) *cos(grt);
                     chgq_l(i,1,is) = chgq_l(i,1,is) -c1 *f2 *rhvg_l(i,it)
                  else
                     c1 = cos(grt);  s1 = -sin(grt)
                     chgq_l(i,1,is) = chgq_l(i,1,is) -c1 *f2 *rhvg_l(i,it)
                     chgq_l(i,2,is) = chgq_l(i,2,is) -s1 *f2 *rhvg_l(i,it)
                  endif
               end do

               f3 = nchg(is) /univol
               if ( ista_kngp == 1) then
                  c1 = iwei(ia);
                  chgq_l(1,1,is) = chgq_l(1,1,is) -c1 *f3
               endif
            end do
         end do
      end if
      deallocate( pos_wk )

    end subroutine case_pre_by_atoms
    
  end subroutine m_CD_predictor_pre

  subroutine m_CD_predictor_post(alpha,beta,rms,nextpl,nfout,printable)
    real(kind=DP), intent(in) :: alpha,beta,rms
    integer, intent(in) :: nextpl
    integer, intent(in) :: nfout
    logical, intent(in) :: printable

    if(sw_extrapolate_charge==ON)then
      call extrapolate_charge()
    endif

    if ( mag_moment0_atoms_is_defined ) then
       call case_post_by_atoms
    else
       call case_post_by_atomtype
    endif

    call m_CD_cp_chgq_to_chgqo()

  contains

    subroutine case_post_by_atomtype
      integer :: ik,it,is,ig
      real(kind=DP) :: fac,pm

      if ( noncol ) then
         Do is=1, ndim_magmom
            do it=1, ntyp
               if ( is == 1 ) then
                  fac = 1.0
               else
                  fac = zeta1(it) *mag_direction0_atomtyp(it,is-1)
               endif
               
               do ik=1,kimg
                  do ig=ista_kngp,iend_kngp
                     chgq_l(ig,ik,is) = chgq_l(ig,ik,is) &
                          &         +fac*rhvg_l(ig,it)*zfm3_l(ig,it,ik)
                  end do
               end do
               if(ista_kngp==1) then
                  chgq_l(1,1,is) = chgq_l(1,1,is) + zfm3_l(1,it,1)*fac*(ival(it)/univol)
               endif
            end do
         end do

      else
         fac=1.d0; pm=1.d0
         do is=1,nspin
            if(is==2) pm=-1.d0
            do it=1,ntyp
               if(nspin>1) fac = (1.d0+pm*zeta1(it))*0.5d0
               do ik=1,kimg
                  do ig=ista_kngp,iend_kngp
                     chgq_l(ig,ik,is) = chgq_l(ig,ik,is) &
                          &            +fac*rhvg_l(ig,it)*zfm3_l(ig,it,ik)
                  enddo
               enddo
               if(ista_kngp==1) then
                  chgq_l(1,1,is) = chgq_l(1,1,is) + zfm3_l(1,it,1)*fac*(ival(it)/univol)
               endif
            enddo
         enddo
      endif
    end subroutine case_post_by_atomtype

    subroutine case_post_by_atoms
      integer :: ia, it, is, i
      real(kind=DP) :: f2, f3, c1, s1
      real(kind=DP) :: grt, nchg(4), ncharge
      real(kind=DP) :: fac,pm

      if ( noncol ) then
         do ia=1, natm
            it = ityp(ia)
            ncharge = ival(it) +ionic_charge_atoms(ia)
            nchg(1) = ncharge
            nchg(2) = mag_moment0_atoms(ia,1)
            nchg(3) = mag_moment0_atoms(ia,2)
            nchg(4) = mag_moment0_atoms(ia,3)

            Do is=1, ndim_magmom
               f2 = nchg(is) /ival(it)
               do i=ista_kngp,iend_kngp
                  grt = pos(ia,1) *ngabc_kngp_l(i,1) &
                       &   +pos(ia,2) *ngabc_kngp_l(i,2) &
                       &   +pos(ia,3 )*ngabc_kngp_l(i,3)
                  grt = grt *PAI2

                  if ( kimg == 1 ) then
                     c1 = iwei(ia) *cos(grt);
                     chgq_l(i,1,is) = chgq_l(i,1,is) +c1 *f2 *rhvg_l(i,it)
                  else
                     c1 = cos(grt);  s1 = -sin(grt)
                     chgq_l(i,1,is) = chgq_l(i,1,is) +c1 *f2 *rhvg_l(i,it)
                     chgq_l(i,2,is) = chgq_l(i,2,is) +s1 *f2 *rhvg_l(i,it)
                  endif
               end do

               f3 = nchg(is) /univol
               if ( ista_kngp == 1) then
                  c1 = iwei(ia);
                  chgq_l(1,1,is) = chgq_l(1,1,is) +c1 *f3
               endif
            end do
         end do

      else
         do ia = 1, natm
            it = ityp(ia)

            ncharge = ival(it) +ionic_charge_atoms(ia)
            if ( nspin == 1 ) then
               nchg(1) = ncharge
            else
               nchg(1) = ( ncharge +mag_moment0_atoms(ia,1) ) /2.0d0
               nchg(2) = ( ncharge -mag_moment0_atoms(ia,1) ) /2.0d0
            endif

            Do is=1, nspin
               f2 = nchg(is) /ival(it)
               do i = ista_kngp, iend_kngp
                  grt = pos(ia,1) *ngabc_kngp_l(i,1) &
                       &   +pos(ia,2) *ngabc_kngp_l(i,2) &
                       &   +pos(ia,3 )*ngabc_kngp_l(i,3)
                  grt = grt *PAI2
                  if ( kimg == 1 ) then
                     c1 = iwei(ia) *cos(grt);
                     chgq_l(i,1,is) = chgq_l(i,1,is) +c1 *f2 *rhvg_l(i,it)
                  else
                     c1 = cos(grt);  s1 = -sin(grt)
                     chgq_l(i,1,is) = chgq_l(i,1,is) +c1 *f2 *rhvg_l(i,it)
                     chgq_l(i,2,is) = chgq_l(i,2,is) +s1 *f2 *rhvg_l(i,it)
                  endif
               end do

               f3 = nchg(is) /univol
               if ( ista_kngp == 1) then
                  c1 = iwei(ia);
                  chgq_l(1,1,is) = chgq_l(1,1,is) +c1 *f3
               endif
            end do
         end do
      end if
    end subroutine case_post_by_atoms

    subroutine extrapolate_charge()
       if(nextpl<2) return ! not ready
       if(rms<rms_threshold)then
          chgq_l(:,:,:) = chgq_l(:,:,:)  &
      & + alpha*(extpl_target_history(:,:,:,1)-extpl_target_history(:,:,:,2)) &
      & + beta* (extpl_target_history(:,:,:,2)-extpl_target_history(:,:,:,3))
       endif
    end subroutine extrapolate_charge

  end subroutine m_CD_predictor_post

  subroutine m_CD_map_chgq_l_to_fft_l_box(nfftp_l,afft_l,is,chgqopt_l)
    use m_Parallelization, only :  chgq_fftcd_maxrecv, chgq_fftcd_maxsend &
         &                       , chgq_fftcd_index, chgq_fftcd_dist &
         &                       , chgq_fftcd_send, chgq_fftcd_recv &
         &                       , chgq_fftcd_scnt, chgq_fftcd_rcnt
    use m_Files,           only : nfout
    integer, intent(in) :: nfftp_l, is
    real(kind=DP), dimension(nfftp_l),                  intent(out) :: afft_l
    real(kind=DP), optional, dimension(ista_kngp:iend_kngp,kimg), intent(in)  :: chgqopt_l

    integer :: i, j
    real(kind=DP), allocatable, dimension(:,:) :: sendbuf,recvbuf
    integer :: icnt_send, icnt_recv, lrank
    integer :: iend, iadd, mpi_comm_all
    integer, parameter :: itag = 19
#ifdef USE_ALLTOALLV
    integer, allocatable, dimension(:) :: sdsp, rdsp
#endif
    integer, dimension(0:nrank_g-1) ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)::sta_r, sta_s

    mpi_comm_all = mpi_ke_world
    if (chgq_fftcd_maxsend /= 0) then
       allocate(sendbuf(chgq_fftcd_maxsend*kimg,0:nrank_g-1))
       sendbuf = 0.0d0
       if (kimg == 1) then
          iend = iend_kngp
          if (iend > kgp_reduced) iend = kgp_reduced
          if(present(chgqopt_l)) then
             do i = ista_kngp, iend
                iadd = i-ista_kngp+1
#ifdef CD_FFT_ALL
                if (chgq_fftcd_dist(iadd) < 0) cycle
#endif
                sendbuf(chgq_fftcd_index(iadd),chgq_fftcd_dist(iadd)) = chgqopt_l(i,1)
             end do
          else
             do i = ista_kngp, iend
                iadd = i-ista_kngp+1
#ifdef CD_FFT_ALL
                if (chgq_fftcd_dist(iadd) < 0) cycle
#endif
                sendbuf(chgq_fftcd_index(iadd),chgq_fftcd_dist(iadd)) = chgq_l(i,1,is)
             end do
          end if
       else
          iend = min(iend_kngp,kgp_reduced)
          if(present(chgqopt_l)) then
             do i = ista_kngp, iend
                iadd = i-ista_kngp+1
#ifdef CD_FFT_ALL
                if (chgq_fftcd_dist(iadd) < 0) cycle
#endif
                sendbuf(chgq_fftcd_index(iadd)*2-1,chgq_fftcd_dist(iadd)) = chgqopt_l(i,1)
                sendbuf(chgq_fftcd_index(iadd)*2,  chgq_fftcd_dist(iadd)) = chgqopt_l(i,2)
             end do
          else
             do i = ista_kngp, iend
                iadd = i-ista_kngp+1
#ifdef CD_FFT_ALL
                if (chgq_fftcd_dist(iadd) < 0) cycle
#endif
                sendbuf(chgq_fftcd_index(iadd)*2-1,chgq_fftcd_dist(iadd)) = chgq_l(i,1,is)
                sendbuf(chgq_fftcd_index(iadd)*2,  chgq_fftcd_dist(iadd)) = chgq_l(i,2,is)
             end do
          end if
       end if
    end if
    if (chgq_fftcd_maxrecv /= 0) then
       allocate(recvbuf(chgq_fftcd_maxrecv*kimg,0:nrank_g-1))
       recvbuf = 0.0d0
    endif

#ifndef USE_ALLTOALLV
    icnt_recv = 0

    do lrank = 0, nrank_g - 1
       if (chgq_fftcd_rcnt(lrank) /= 0) then
          call mpi_irecv(recvbuf(1,lrank),chgq_fftcd_rcnt(lrank)*kimg,mpi_double_precision, &
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

#else

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
#endif

    afft_l = 0.0d0
    if (kimg == 1) then
#ifdef FFT_3D_DIVISION_CD
       do i = 0, nrank_g - 1
          do j = 1, chgq_fftcd_rcnt(i)
             afft_l(chgq_fftcd_recv(j,i)*2-1) = recvbuf(j,i)
             afft_l(chgq_fftcd_recv(j,i)*2  ) = 0.0d0
          enddo
       enddo
#else
       do i = 0, nrank_g - 1
          do j = 1, chgq_fftcd_rcnt(i)
             afft_l(chgq_fftcd_recv(j,i)) = recvbuf(j,i)
          enddo
       enddo
#endif
    else
       do i = 0, nrank_g - 1
          do j = 1, chgq_fftcd_rcnt(i)
             afft_l(chgq_fftcd_recv(j,i)*2-1) = recvbuf(j*2-1,i)
             afft_l(chgq_fftcd_recv(j,i)*2  ) = recvbuf(j*2  ,i)
          enddo
       enddo
    endif
    if (allocated(sendbuf)) deallocate(sendbuf)
    if (allocated(recvbuf)) deallocate(recvbuf)
  end subroutine m_CD_map_chgq_l_to_fft_l_box

  subroutine m_CD_partial_charge_ek(nfout,kv3)
    use m_ES_occup,           only : m_ESoc_set_nEwindows_pc, m_ESoc_keep_occup, m_ESoc_if_elec_state &
                  , m_ESoc_substitute_occup, m_ESoc_retrieve_occup, m_ESoc_free_nEwindows &
                  , m_ESoc_check_if_metalic
    use m_Control_Parameters, only : sw_spin_magmom_rspace
    use m_Kpoints,            only : kv3_ek
    use m_Electronic_Structure, only : efermi
    integer, intent(in) :: nfout,kv3
    integer :: i, iloop, iloop2, nf, ig,ir
    real(kind=DP) :: emin, emax, fac
    nf = 0 ! dummy
    call m_ESoc_check_if_metalic(nfout)
    if(.not.allocated(chgq_l_pc)) then
      nEwindow_ek = 0
      call m_ESoc_set_nEwindows_pc(nfout,nEwindow_ek)
      if ( nEwindow_ek <= 0 )then
         write(nfout,'(a)') '!** nEwindow could not be resolved. Partial charge will NOT be evaluated.'
         return
      else
         allocate(chgq_l_pc(ista_kngp:iend_kngp,kimg,nspin,nEwindow_ek));chgq_l_pc = 0.d0
      endif
    endif
    call m_ESoc_keep_occup()
    call m_CD_keep_chgq_l()
    call m_CD_keep_retrieve_hsr(.true.)
    fac = dble(kv3)/dble(kv3_ek)
    do i = 1, nEwindow_ek
       if(m_ESoc_if_elec_state(nfout,i,emin,emax) /= YES) cycle
       call m_ESoc_substitute_occup(nfout,i)
       call m_CD_softpart_3D(nfout,kv3)
       call m_CD_hardpart(nfout,kv3)
       chgq_l_pc(ista_kngp:iend_kngp,1:kimg,1:nspin,i) = &
       chgq_l_pc(ista_kngp:iend_kngp,1:kimg,1:nspin,i) + &
       chgq_l   (ista_kngp:iend_kngp,1:kimg,1:nspin)*fac
    enddo
    call m_ESoc_retrieve_occup()
    call m_CD_retrieve_chgq()
    call m_CD_keep_retrieve_hsr(.false.)
  end subroutine m_CD_partial_charge_ek

  subroutine m_CD_partial_charge_ek_finalize()
    use m_ES_occup, only : m_ESoc_free_nEwindows
    if(allocated(chgq_l_pc)) deallocate(chgq_l_pc)
    call m_ESoc_free_nEwindows()
  end subroutine m_CD_partial_charge_ek_finalize

end module m_Charge_Density
