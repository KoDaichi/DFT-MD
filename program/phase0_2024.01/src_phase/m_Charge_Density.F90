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

!!$#define _VECTOR_TUNING_

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
!!$#define DEBUG_LDOS

module m_Charge_Density
! $Id: m_Charge_Density.F90 633 2020-12-01 05:11:03Z jkoga $
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
       &                          , ista_fftp, iend_fftp, mp_fftp &
       &                          , nis_fftp, nie_fftp &
       &                          , is_atm, ie_atm, nel_atm, np_atm, mp_atm &
       &                          , ista_atm, iend_atm &
       &                          , ista_kngp_gw, iend_kngp_gw &
       &                          , mpi_spin_group, ista_spin, iend_spin, nrank_s
  use m_Timing,              only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters,  only : nspin,ipri,ipriwf,iprichargedensity,cdel_critical, nel_Ylm &
       &                          , istress,sw_fine_STM_simulation,kimg,af,neg &
       &                          , charge_filetype, charge_title, initial_chg &
       &                          , sw_add_corecharge_rspace, eval_corecharge_on_Gspace &
       &                          , iprichargemixing, ipritotalcharge &
       &                          , initial_charge_filetype &
       &                          , m_CtrlP_cachesize &
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
  use m_FFT,                 only : m_FFT_alloc_WF_work &
       &                          , m_FFT_dealloc_WF_work &
       &                          , m_FFT_alloc_CD_box, m_FFT_dealloc_CD_box &
!fj$$       &                          , m_FFT_WF,fft_box_size_WF, nfft,nfftp,nfftps &
       &                          , fft_box_size_WF, nfft,nfftp,nfftps &
       &                          , fft_box_size_CD_nonpara, nfftp_nonpara &
       &                          , m_FFT_CD0, m_FFT_CD_inverse0, fft_box_size_CD &
       &                          , fft_box_size_CD_prev, fft_box_size_CD_nonpara_prev
  use m_FFT,                 only : m_FFT_WF
  use m_Kpoints,             only : k_symmetry
  use m_Electronic_Structure,only : totch,neordr,occup_l,fsr_l,fsi_l, zaj_l
  use m_Electronic_Structure,only : m_ES_WF_in_Rspace, m_ES_dealloc_chgsoft, chg_has_been_calculated, chg_softpart

! ===================== added by  K. Tagami ================== 5.0
  use m_Control_Parameters,   only : sw_mix_charge_hardpart
  use m_Control_Parameters,   only : sw_update_charge_hsr
! =========================================================== 5.0

! ====================================== added by K. Tagami ============ 11.0&13.0U
  use m_Control_Parameters,  only : noncol, ndim_spinor, ndim_chgpot, ndim_magmom, &
       &                            import_collinear_spindensity, &
       &                            previous_nspin_collinear, &
       &                            SpinOrbit_mode
  use m_ES_NonCollinear,      only : m_ES_DensMat_to_MagMom_Gspace, &
       &                             m_ES_DensMat_to_MagMom_hsr, &
       &                             m_ES_set_Mat_hsr_with_soc
  use m_Ionic_System,         only : mag_direction0_atomtyp, mag_moment0_atomtyp, &
       &                             ionic_charge_atomtyp, &
       &                             mag_moment0_atomtyp_is_defined, &
       &                             ionic_charge_atoms, mag_moment0_atoms, &
       &                             mag_moment0_atoms_is_defined
  use m_Const_Parameters,    only : BuiltIn
  use m_Crystal_Structure,      only :  op, sw_fix_global_quantz_axis, &
       &                                Global_Quantz_Axis_Fixed
  use m_CS_Magnetic,      only  : invop, magmom_dir_inversion_opr_flag, determinant_op
! ===================================================================== 11.0&13.0U

  use m_ErrorMessages,        only : EOF_REACHED

  use m_Realspace, only : nmesh_rs_aug,nmesh_rs_aug_max,meshxyz_rs_aug,qr_clm_ylm,nlmtpair,plmt1,plmt2

! ==== KT_add ======= 2014/08/25
  use m_PseudoPotential,     only : nlmt_add, ilmt_add, ltp_add, mtp_add, lmta_add
  use m_Electronic_Structure,only : fsr_add_l, fsi_add_l
! ==================== 2014/08/25

! ==== EXP_CELLOPT === 2015/09/24
  use m_Parallelization, only : ista_fftph, iend_fftph, idisp_fftp, nel_fftp, npes_cdfft
  use m_FFT,                only : m_FFT_CD_direct, m_FFT_CD_inverse_c
! ==================== 2015/09/24
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
  real(kind=DP),private,allocatable, dimension(:,:)     :: afft_spin
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
!  71. add_hardpart_to_chgq_enl_IA    <-(@68)
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
          call mpi_allreduce(ylm_t,ylm_mpi,kgp,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
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
       call mpi_allreduce(igfp_t,igfp_mpi,kgp,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
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
          call mpi_allreduce(qitg_t,qitg_mpi,kgp,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
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
          call mpi_allreduce(ngpt_t,ngpt_mpi,kgp,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
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

! =============================== modified by K. Tagami ============== 11.0
!  subroutine m_CD_alloc_chgq
!    allocate(chgq_l(ista_kngp:iend_kngp,kimg,nspin)); chgq_l = 0.d0
!    allocate(chgqo_l(ista_kngp:iend_kngp,kimg,nspin)); chgqo_l = 0.d0
!    if(istress == ON .or. sw_fine_STM_simulation == ON) then
!       allocate(chgsoft(ista_kngp:iend_kngp,kimg,nspin))
!       chgsoft = 0.d0
!    end if
!  end subroutine m_CD_alloc_chgq

!  subroutine m_CD_alloc_chgsoft
!    if(.not.allocated(chgsoft)) then
!       allocate(chgsoft(ista_kngp:iend_kngp,kimg,nspin))
!       chgsoft = 0.d0
!    end if
!  end subroutine m_CD_alloc_chgsoft

  subroutine m_CD_alloc_chgq
    if(.not.allocated(chgq_l)) then
      allocate(chgq_l(ista_kngp:iend_kngp,kimg,ndim_magmom)); chgq_l = 0.d0
    endif
    if(.not.allocated(chgqo_l)) then
      allocate(chgqo_l(ista_kngp:iend_kngp,kimg,ndim_magmom)); chgqo_l = 0.d0
    endif

    if(istress == ON .or. sw_fine_STM_simulation == ON) then
       if(.not.allocated(chgsoft)) then
         allocate(chgsoft(ista_kngp:iend_kngp,kimg,ndim_magmom))
         chgsoft = 0.d0
       endif
    end if
    if(sw_extrapolate_charge==ON)then
       if(.not.allocated(extpl_target_history)) allocate(extpl_target_history(ista_kngp:iend_kngp,kimg,ndim_magmom,3))
       extpl_target_history = 0.d0
    endif
  end subroutine m_CD_alloc_chgq

  subroutine m_CD_alloc_chgsoft
    if(.not.allocated(chgsoft)) then
       allocate(chgsoft(ista_kngp:iend_kngp,kimg,ndim_magmom))
       chgsoft = 0.d0
    end if
  end subroutine m_CD_alloc_chgsoft
! ================================================================ 11.0


! ============================== modified by K. Tagami ========== 11.0
!  subroutine m_CD_alloc_hsr()
!!!$    print '(" natm, nlmt, nspin = ",3i5)',natm,nlmt,nspin
!    if(flg_paw) then
!       allocate(hsr(natm,nlmt,nlmt,nspin)); hsr = 0.d0
!    endif
!    if ( flg_paw .or. sw_mix_charge_hardpart == ON ) then
!       allocate(hsro(natm,nlmt,nlmt,nspin)); hsro = 0.d0
!    endif
!  end subroutine m_CD_alloc_hsr

  subroutine m_CD_alloc_hsr()
    if ( noncol ) then
      if (.not.allocated(hsr)) allocate(hsr(natm,nlmt,nlmt,ndim_magmom)); hsr = 0.d0
      if (.not.allocated(hsi)) allocate(hsi(natm,nlmt,nlmt,ndim_magmom)); hsi = 0.d0
!      if ( flg_paw .or. sw_mix_charge_hardpart == ON ) then
         if(.not.allocated(hsro)) allocate(hsro(natm,nlmt,nlmt,ndim_magmom)); hsro = 0.d0
         if(.not.allocated(hsio)) allocate(hsio(natm,nlmt,nlmt,ndim_magmom)); hsio = 0.d0
!      endif
    else
      if(.not.allocated(hsr)) allocate(hsr(natm,nlmt,nlmt,nspin)); hsr = 0.d0
!      if ( flg_paw .or. sw_mix_charge_hardpart == ON ) then
         if(.not.allocated(hsro)) allocate(hsro(natm,nlmt,nlmt,nspin)); hsro = 0.d0
!      endif
    endif
  end subroutine m_CD_alloc_hsr
! ================================================================= 11.0

! ====== KT_add ===== 2014/08/25
  subroutine m_CD_alloc_hsr_add()
    if (.not.allocated(hsr_add)) then
       allocate(hsr_add(natm,nlmt_add,nlmt_add,ndim_magmom)); hsr_add = 0.d0
    endif
    if (.not.allocated(hsi_add)) then
       allocate(hsi_add(natm,nlmt_add,nlmt_add,ndim_magmom)); hsr_add = 0.d0
    endif
  end subroutine m_CD_alloc_hsr_add
! =================== 2014/08/25


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

    if(mype == 0) ipritotalcharge_0 = ipritotalcharge
    if(npes > 1) call mpi_bcast(ipritotalcharge_0,1,mpi_integer,0,MPI_CommGroup,ierr)

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
    if(npes > 1) then
       do i = 0, npes-1
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
    if(npes > 1) then
       call mpi_barrier(MPI_CommGroup,ierr)
       call mpi_bcast(totch_old,1,mpi_double_precision,ip,MPI_CommGroup,ierr)
       call mpi_bcast(totch_new,1,mpi_double_precision,ip,MPI_CommGroup,ierr)
    end if
    if(nspin == 2 .and. af == 0) then
       if(npes > 1) then
          call mpi_barrier(MPI_CommGroup,ierr)
          call mpi_bcast(chg_t,4,mpi_double_precision,ip,MPI_CommGroup,ierr)
       end if
       if((ipritotalcharge==1 .and. iteration_electronic==1).or.ipritotalcharge>=2) &
            & write(nfout,96) "OLD",chg_t(1,OLD)*univol,chg_t(2,OLD)*univol,totch_old*univol
       if(ipritotalcharge>=1) &
            & write(nfout,96) "NEW",chg_t(1,NEXT)*univol,chg_t(2,NEXT)*univol,totch_new*univol
       if ( ipritotalcharge_0 >=2 ) call m_CD_calc_abs_magetization( nfout )
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

! ================================ added by K. Tagami ================ 11.0
  subroutine m_CD_check_noncl(nfout)
    integer, intent(in) :: nfout
    integer :: i, ispin, ip
    integer :: is, ni

    integer :: ipritotalcharge_0
    real(kind=DP) :: totch_old(ndim_magmom), totch_new(ndim_magmom)

    real(kind=DP),allocatable,dimension(:,:) :: chg_t ! d(ndim_magmom,OLD:NEXT)

    if(mype == 0) ipritotalcharge_0 = ipritotalcharge
    if(npes > 1) call mpi_bcast(ipritotalcharge_0,1,mpi_integer,0,MPI_CommGroup,ierr)

!    if ( .not.(ipritotalcharge_0 >= 2 .or. &
!       &     (ipritotalcharge_0 >= 1 .and. (ndim_spinor==2)) )) return
!!!!!!    if ( .not.(ipritotalcharge_0 >= 2) ) return

!!!!!!    allocate( chg_t(ndim_chgpot,OLD:NEXT) );   chg_t = 0.0d0
    totch_old = 0.d0;     totch_new = 0.d0

    ip = 0
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
       ni = 1
       totch_old(:) = totch_old(:) + chgqo_l(i,1,:)
       totch_new(:) = totch_new(:) + chgq_l(i,1,:)
    end if
    if(npes > 1) then
       call mpi_barrier(MPI_CommGroup,ierr)
       call mpi_bcast( totch_old, ndim_magmom, mpi_double_precision, &
            &          ip, MPI_CommGroup, ierr )
       call mpi_bcast( totch_new, ndim_magmom, mpi_double_precision,&
            &          ip, MPI_CommGroup, ierr )
    end if

    if(ipritotalcharge >= 1) then
       write(nfout,'(A,F14.8,A,F14.8,A,F14.8,A,F14.8)') &
            &        '!OLD Chg **  Tot:', totch_old(1)*univol, &
            &                     ' Mx:', totch_old(2)*univol, &
            &                     ' My:', totch_old(3)*univol, &
            &                     ' Mz:', totch_old(4)*univol
       write(nfout,'(A,F14.8,A,F14.8,A,F14.8,A,F14.8)') &
            &        '!NEW Chg **  Tot:', totch_new(1)*univol, &
            &                     ' Mx:', totch_new(2)*univol, &
            &                     ' My:', totch_new(3)*univol, &
            &                     ' Mz:', totch_new(4)*univol
    endif
    if ( ipritotalcharge_0 >=2 ) call m_CD_calc_abs_magetization( nfout )

  end subroutine m_CD_check_noncl
! ===================================================================== 11.0

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
    call mpi_allreduce(cdelt_mpi,cdelt,1 &
                   &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)

!!$    call m_CD_dealloc_rspace_charge()
!!$
    total_charge_t = 0.d0
                                                 __TIMER_DO_START(859)
    if(mype==0) then

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
    call mpi_bcast(total_charge_t,1,mpi_double_precision,0,MPI_CommGroup,ierr)
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
    call tstatc0_begin('m_CD_hardpart ',id_sname,1)

! ====================== added by K. Tagami =================== 5.0
    if ( sw_update_charge_hsr == OFF ) goto 100
! ============================================================= 5.0

    call summation_of_ff(kv3) ! -(m_C.D.) (vnlsum) fsr_l, fsi_l, occup_l --> hsr
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

! ========================= modified by K. Tagami ==================== 11.0
!    call add_hardpart_to_chgq_l(nfout,nspin,hsr)
    call add_hardpart_to_chgq_l( nfout, nspin, hsr, NO )
                        ! (lclchg) add hardpart and make average
! ==================================================================== 11.0


    if ( sw_update_charge_hsr == OFF ) goto 200

    if ( flg_paw .or. sw_mix_charge_hardpart == ON ) then
       if ( charge_symm_mode <= chg_symm_level1 ) then
          call symmtrz_of_ff
       else if ( charge_symm_mode == chg_symm_level2 ) then
          call summation_of_ff_with_symmtrz(kv3)
       endif
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

    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(720)
  end subroutine m_CD_hardpart

! ========================== added by K. Tagami ================= 11.0
  subroutine m_CD_hardpart_noncl(nfout,kv3)
    integer, intent(in) :: nfout, kv3

    integer, parameter :: DEBUGPRINTLEVEL = 2
    integer :: id_sname = -1
    call tstatc0_begin('m_CD_hardpart_noncl ',id_sname,1)

    call summation_of_ff_noncl(kv3)
                   ! -(m_C.D.) (vnlsum) fsr_l, fsi_l, occup_l --> hsr

    if (flg_paw .and. iprichargedensity >= DEBUGPRINTLEVEL) then
       if(flg_paw) write(nfout,'(" -- hsr before symmtrz_of_ff --")')
       call wd_hsr_noncl(nfout)
    end if

!    if ( flg_paw .or. sw_mix_charge_hardpart == ON ) then
!      call m_CD_symmtrz_of_ff_noncl_C( hsr, hsi )
!    endif

    if(iprichargedensity >= DEBUGPRINTLEVEL) then
       write(nfout,'(" -- hsr --")')
       if(flg_paw) write(nfout,'(" -- hsr after symmtrz_of_ff --")')
       call wd_hsr_noncl(nfout)
    end if

    call add_hardpart_to_chgq_l( nfout, ndim_magmom, hsr, NO )

    if ( flg_paw .or. sw_mix_charge_hardpart == ON ) then
      call m_CD_symmtrz_of_ff_noncl_C( hsr, hsi )
    endif

    if(iprichargedensity >= DEBUGPRINTLEVEL) then
       write(nfout,'(" -- hardpart summed --")')
       call m_CD_wd_chgq_l_portion_noncl(nfout)
    endif

    call tstatc0_end(id_sname)
  end subroutine m_CD_hardpart_noncl
! ================================================================ 11.0

  subroutine m_CD_hardpart_hsr(nfout,kv3)
    integer, intent(in) :: nfout, kv3
    integer, parameter :: DEBUGPRINTLEVEL = 3

    integer :: id_sname = -1
                                                  __TIMER_SUB_START(734)
    call tstatc0_begin('m_CD_hardpart_hsr ',id_sname,1)

    call summation_of_ff(kv3) ! -(m_C.D.) (vnlsum) fsr_l, fsi_l, occup_l --> hsr

    if(iprichargedensity >= DEBUGPRINTLEVEL) then
       write(nfout,'(" -- hsr --")')
       call wd_hsr(nfout)
    end if
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(734)
  end subroutine m_CD_hardpart_hsr

! ==================================== added by K. Tagami =============== 11.0
  subroutine m_CD_hardpart_hsr_noncl(nfout,kv3)
    integer, intent(in) :: nfout, kv3
    integer, parameter :: DEBUGPRINTLEVEL = 2

    integer :: id_sname = -1
                                                  __TIMER_SUB_START(734)
    call tstatc0_begin('m_CD_hardpart_hsr_noncl ',id_sname,1)

    call summation_of_ff_noncl(kv3) ! -(m_C.D.) (vnlsum) fsr_l, fsi_l, occup_l --> hsr

    if(iprichargedensity >= DEBUGPRINTLEVEL) then
       write(nfout,'(" -- hsr --")')
       call wd_hsr_noncl(nfout)
    end if
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(734)
  end subroutine m_CD_hardpart_hsr_noncl
! ========================================================================= 11.0

  subroutine summation_of_ff(kv3)
    integer, intent(in) :: kv3

    real(kind=DP)   :: w_n, d_factor
    integer         :: ia, is, k, i, lmt1, lmt2, p, q, it
    real(kind=DP),pointer,dimension(:,:,:,:) :: hsr_mpi

    d_factor = 2.d0/kv3
    hsr = 0.d0
    do ia = 1, natm
       it = ityp(ia)
!       do is = 1, nspin, af+1
       do is = ista_spin, iend_spin, af+1
          do i = 1, np_e                                  ! MPI
             do k = is, kv3+is-nspin, nspin
                if(map_k(k) /= myrank_k) cycle            ! MPI
                w_n = occup_l(i,k)*d_factor
                do lmt1 = 1, ilmt(it)
                   p = lmta(lmt1,ia)
                   if(k_symmetry(k) == GAMMA) then
                      do lmt2 = lmt1, ilmt(it)
                         q = lmta(lmt2,ia)
                         hsr(ia,lmt1,lmt2,is) = hsr(ia,lmt1,lmt2,is) + w_n * &
                              &  (fsr_l(i,p,k)*fsr_l(i,q,k))
                      end do! lmt2
                   else
                      do lmt2 = lmt1, ilmt(it)
                         q = lmta(lmt2,ia)
                         hsr(ia,lmt1,lmt2,is) = hsr(ia,lmt1,lmt2,is) + w_n * &
                              &  (fsr_l(i,p,k)*fsr_l(i,q,k) + fsi_l(i,p,k)*fsi_l(i,q,k))
                      end do! lmt2
                   end if
                end do! lmt1
             end do! ik
          end do! i
       end do! is
    end do! ia
! ============================= added by K. Tagami ===================== 5.0
! ------------------------------------------------- just in case
    Do ia = 1, natm
       it = ityp(ia)
       Do is = 1, nspin, af+1
          Do lmt1=1, ilmt(it)
            Do lmt2=lmt1, ilmt(it)
               hsr(ia,lmt2,lmt1,is) = hsr(ia,lmt1,lmt2,is)
            End do
         End do
       End do
    End do
! ======================================================================= 5.0

    if(npes >= 2) then
       allocate(hsr_mpi(natm,nlmt,nlmt,nspin))
! ======================================= Adde by K. Tagami ===========
! === DEBUG by tkato 2011/12/07 ================================================
!	hsr_mpi = 0
       hsr_mpi = 0.0d0
! ==============================================================================
! ======================================================================
       call mpi_allreduce(hsr,hsr_mpi,natm*nlmt*nlmt*nspin,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
!       call mpi_allreduce(hsr,hsr_mpi,natm*nlmt*nlmt*nspin,mpi_double_precision,mpi_sum,mpi_spin_group,ierr)
       hsr = hsr_mpi
       deallocate(hsr_mpi)
    end if
  end subroutine summation_of_ff

  subroutine summation_of_ff_with_symmtrz(kv3)
    integer, intent(in) :: kv3

    real(kind=DP)   :: w_n, d_factor
    integer         :: ia, is, ik, i, lmt1, lmt2, lmt3, lmt4, p, q, it, j
    integer :: il1,im1,it1,il2,im2,it2,il3,im3,it3,il4,im4,it4
    integer :: ii,jj,kk,ll,n,m,iii,jjj
    integer :: ja
    integer :: ikbz, iopr
    real(kind=DP),pointer,dimension(:,:,:,:) :: hsr_tmp1, hsr_tmp2

    d_factor = 2.d0 /dble(kv3)
    hsr = 0.d0
!
    allocate(hsr_tmp1(natm,nlmt,nlmt,nspin))
    allocate(hsr_tmp2(natm,nlmt,nlmt,nspin))

!    Do is=1, nspin
    Do is=ista_spin, iend_spin
       Do ik=is, kv3, nspin
          if(map_k(ik) /= myrank_k) cycle            ! MPI

          hsr_tmp1 = 0.0d0

! ----- part 1 -------
          do ia = 1, natm
             it = ityp(ia)
             do i = 1, np_e                                  ! MPI
                w_n = occup_l(i,ik)*d_factor
                do lmt1 = 1, ilmt(it)
                   p = lmta(lmt1,ia)

                   if(k_symmetry(ik) == GAMMA) then
                      do lmt2 = lmt1, ilmt(it)
                         q = lmta(lmt2,ia)
                         hsr_tmp1(ia,lmt1,lmt2,is) = hsr_tmp1(ia,lmt1,lmt2,is) &
                              &  + w_n *(fsr_l(i,p,ik)*fsr_l(i,q,ik))
                      end do! lmt2
                   else
                      do lmt2 = lmt1, ilmt(it)
                         q = lmta(lmt2,ia)
                         hsr_tmp1(ia,lmt1,lmt2,is) = hsr_tmp1(ia,lmt1,lmt2,is) &
                              & + w_n *( fsr_l(i,p,ik)*fsr_l(i,q,ik) &
                              &          +fsi_l(i,p,ik)*fsi_l(i,q,ik) )
                      end do! lmt2
                   end if
                end do! lmt1
             end do ! i
          end do! ia

          Do ia = 1, natm
             it = ityp(ia)
             Do lmt1=1, ilmt(it)
                Do lmt2=lmt1, ilmt(it)
                   hsr_tmp1(ia,lmt2,lmt1,is) = hsr_tmp1(ia,lmt1,lmt2,is)
                End do
             End do
          End do
! ----- part 2 -------
          hsr_tmp2 = 0.0d0
          Loop_j: Do j=1, num_star_of_k(ik)
             ikbz = star_of_k(ik,j)
             iopr = iopr_k_fbz_to_ibz(ikbz)

             do ia = 1, natm
                it = ityp(ia)
                ja=abs(ia2ia_symmtry_op_inv(ia,iopr))
                do lmt1 = 1, ilmt(it)
                   il1=ltp(lmt1,it);  im1=mtp(lmt1,it);   it1=taup(lmt1,it)
                   ii=(il1-1)**2+im1
                   do lmt2 = lmt1, ilmt(it)
                      il2=ltp(lmt2,it);  im2=mtp(lmt2,it); it2=taup(lmt2,it)
                      jj=(il2-1)**2+im2

                      do n=1,nylm_paw(ii,iopr,ia)
                         iii=iylm_paw(n,ii,iopr,ia)
                         do m=1,nylm_paw(jj,iopr,ia)
                            jjj=iylm_paw(m,jj,iopr,ia)

                            do lmt3=1,ilmt(it)
                               il3=ltp(lmt3,it); im3=mtp(lmt3,it); it3=taup(lmt3,it)
                               kk=(il3-1)**2+im3
                               if(kk.ne.iii .or. it1.ne.it3) cycle
                               do lmt4=1,ilmt(it)
                                  il4=ltp(lmt4,it); im4=mtp(lmt4,it); it4=taup(lmt4,it)
                                  ll=(il4-1)**2+im4
                                  if(ll.ne.jjj .or. it2.ne.it4) cycle

                                  hsr_tmp2(ia,lmt1,lmt2,is) = &
                                       hsr_tmp2(ia,lmt1,lmt2,is) &
                                       & + hsr_tmp1(ja,lmt3,lmt4,is) &
                                       &  *crotylm_paw(n,ii,iopr,ia) &
                                       &  *crotylm_paw(m,jj,iopr,ia)
                               end do! lmt4
                            end do! lmt3

                         end do! jjj
                      end do! iii

                   end do! lmt2
                end do! lmt1
             end do! ia
          End do Loop_j
          hsr_tmp2 = hsr_tmp2 /dble( num_star_of_k(ik) )
          hsr = hsr +hsr_tmp2
       End Do
    End Do
    deallocate( hsr_tmp1 )

    if(npes >= 2) then
       hsr_tmp2 = 0.0d0
!       call mpi_allreduce( hsr, hsr_tmp2, natm*nlmt*nlmt*nspin, &
!            &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
       call mpi_allreduce( hsr, hsr_tmp2, natm*nlmt*nlmt*nspin, &
            &              mpi_double_precision, mpi_sum, mpi_spin_group, ierr )
       hsr = hsr_tmp2
    end if
    deallocate(hsr_tmp2)

  end subroutine summation_of_ff_with_symmtrz

! ======================================= added by K. Tagami ============ 11.0
  subroutine summation_of_ff_noncl(kv3)
    integer, intent(in) :: kv3

    real(kind=DP)   :: w_n, d_factor
    integer         :: ia, is, k, i, lmt1, lmt2, p, q, it

    integer :: is1, is2, is_tmp, k1, k2

    real(kind=DP), allocatable, dimension(:,:,:,:) :: hsr_or_hsi_mpi
    real(kind=DP), allocatable, dimension(:,:,:,:) :: hsr_ssrep
    real(kind=DP), allocatable, dimension(:,:,:,:) :: hsi_ssrep
    real(kind=DP), allocatable, dimension(:,:,:,:) :: hsr_with_soc
    real(kind=DP), allocatable, dimension(:,:,:,:) :: hsi_with_soc
!
! ---------------------- kt : uncertain ---------
!!!!!    d_factor = 2.d0/ ( kv3 /ndim_spinor )
!
    d_factor = 1.d0/ ( kv3 /ndim_spinor )
! -------------------------------------------

    allocate( hsr_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsr_ssrep = 0.0d0
    allocate( hsi_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsi_ssrep = 0.0d0

    do ia = 1, natm
      it = ityp(ia)

      do is1 = 1, ndim_spinor
         do is2 = 1, ndim_spinor
            is_tmp = ( is1 -1 )*ndim_spinor + is2
            do i = 1, np_e                                  ! MPI

               do k = 1, kv3, ndim_spinor
                  if ( map_k(k) /= myrank_k ) cycle            ! MPI
                  w_n = occup_l(i,k) *d_factor

                  k1 = k + is1 -1;  k2 = k + is2 -1
                  do lmt1 = 1, ilmt(it)
                     p = lmta(lmt1,ia)
!
! -------------- kt caution : ---- Unfamiliar with treatment of GAMMA point ----
!
!!!!!!!!                     if ( k_symmetry(k) == GAMMA ) then
!!!!!!                        do lmt2 = lmt1, ilmt(it)
!!                        do lmt2 = 1, ilmt(it)
!!                           q = lmta(lmt2,ia)
!!!                          hsr(ia,lmt1,lmt2,is_tmp) = hsr(ia,lmt1,lmt2,is_tmp) &
!!                                & + w_n *( fsr_l(i,p,k1)*fsr_l(i,q,k2) )
!!                        end do! lmt2
!!                     else
!!!!                        do lmt2 = lmt1, ilmt(it)
                        do lmt2 = 1, ilmt(it)
                           q = lmta(lmt2,ia)
                           hsr_ssrep(ia,lmt1,lmt2,is_tmp) &
                                &  = hsr_ssrep(ia,lmt1,lmt2,is_tmp) &
                                &    + w_n * ( fsr_l(i,p,k1)*fsr_l(i,q,k2) &
                                &            + fsi_l(i,p,k1)*fsi_l(i,q,k2) )
                           hsi_ssrep(ia,lmt1,lmt2,is_tmp) &
                                &  = hsi_ssrep(ia,lmt1,lmt2,is_tmp) &
                                &    + w_n * ( -fsr_l(i,p,k1)*fsi_l(i,q,k2) &
                                &              +fsi_l(i,p,k1)*fsr_l(i,q,k2) )
                        end do! lmt2

!!!!                     end if

                  end do! lmt1
               end do! k
            end do! i
         end do! is2
      end do! is1
    end do! ia
!
    if (npes >= 2) then
       allocate(hsr_or_hsi_mpi(natm,nlmt,nlmt,ndim_chgpot))
       hsr_or_hsi_mpi = 0.0d0

!       call mpi_allreduce( hsr_ssrep, hsr_or_hsi_mpi, natm*nlmt*nlmt*ndim_chgpot, &
!        &                  mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
       call mpi_allreduce( hsr_ssrep, hsr_or_hsi_mpi, natm*nlmt*nlmt*ndim_chgpot, &
        &                  mpi_double_precision, mpi_sum, mpi_spin_group, ierr )
       hsr_ssrep = hsr_or_hsi_mpi

!       call mpi_allreduce( hsi_ssrep, hsr_or_hsi_mpi, natm*nlmt*nlmt*ndim_chgpot, &
!        &                  mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
       call mpi_allreduce( hsi_ssrep, hsr_or_hsi_mpi, natm*nlmt*nlmt*ndim_chgpot, &
        &                  mpi_double_precision, mpi_sum, mpi_spin_group, ierr )
       hsi_ssrep = hsr_or_hsi_mpi

       deallocate(hsr_or_hsi_mpi)
    end if
! ---------------
    if ( SpinOrbit_mode == BuiltIn ) then
       allocate( hsr_with_soc( natm, nlmt, nlmt, ndim_chgpot ) ); hsr_with_soc = 0.0d0
       allocate( hsi_with_soc( natm, nlmt, nlmt, ndim_chgpot ) ); hsi_with_soc = 0.0d0

       call m_ES_set_Mat_hsr_with_soc( hsr_ssrep, hsi_ssrep, hsr_with_soc, hsi_with_soc )

       hsr_ssrep = hsr_with_soc;    hsi_ssrep = hsi_with_soc
       deallocate( hsr_with_soc, hsi_with_soc )
    endif
! ---------------

    call m_ES_DensMat_To_MagMom_hsr( natm, nlmt, hsr_ssrep, hsi_ssrep, hsr, hsi )
!-
    deallocate( hsr_ssrep, hsi_ssrep )
  end subroutine summation_of_ff_noncl
! ========================================================================== 11.0

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

! ================================= added by K. Tagami ================== 11.0
  subroutine m_CD_symmtrz_of_ff_noncl_C( hsr_wk, hsi_wk )
    real(kind=DP), intent(inout) :: hsr_wk(natm,nlmt,nlmt,ndim_magmom)
    real(kind=DP), intent(inout) :: hsi_wk(natm,nlmt,nlmt,ndim_magmom)

    integer         :: ia, is, iopr, i, lmt1, lmt2, lmt3, lmt4, it
    real(kind=DP), allocatable, dimension(:,:,:,:) :: hsr_tmp
    real(kind=DP), allocatable, dimension(:,:,:,:) :: hsi_tmp

    integer :: il1,im1,it1,il2,im2,it2,il3,im3,it3,il4,im4,it4
    integer :: ii,jj,kk,ll,n,m,iii,jjj
    integer :: ja

    integer :: ixyz1, ixyz2, is_tmp
    real(kind=DP) :: ctmp1, weight, ctmp2, weight2, fi

    fi = 1.d0 /dble(nopr)
    if ( charge_symm_mode >= chg_symm_level1 ) then
       fi = 1.0d0 /dble(nopr_from_fbz_to_ibz)
    endif

    allocate(hsr_tmp(natm,nlmt,nlmt,ndim_magmom)); hsr_tmp = 0.0d0
    allocate(hsi_tmp(natm,nlmt,nlmt,ndim_magmom)); hsi_tmp = 0.0d0
!
    hsr_tmp = hsr_wk;  hsi_tmp = hsi_wk

    do ia=1,natm
       it=ityp(ia)
       do is =1, ndim_magmom
          do lmt2=1,ilmt(it)
             do lmt1=lmt2+1,ilmt(it)
                hsr_tmp(ia,lmt1,lmt2,is) =  hsr_tmp(ia,lmt2,lmt1,is)
                hsi_tmp(ia,lmt1,lmt2,is) = -hsi_tmp(ia,lmt2,lmt1,is)
             end do
          end do
       end do
    end do

    hsr_wk = 0.d0; hsi_wk = 0.0d0

    do iopr=1,nopr
       if ( charge_symm_mode >= chg_symm_level1 )then
         if(flg_opr_from_fbz_to_ibz(iopr) == 0 ) cycle
       endif

! === KT_add ==== 2014/08/14
       if ( magmom_dir_inversion_opr_flag(iopr) == -1 ) then
          weight = -1.0d0
       else
          weight = 1.0d0
       endif
! ============== 2014/08/14

! === KT_add ==== 2015/07/18  !!! test
       if ( determinant_op(iopr) > 0 ) then
          weight2 = 1.0d0
       else
          weight2 = -1.0d0
       endif
! =============== 2015/07/18

       do ia = 1, natm
          it = ityp(ia)
          ja=abs(ia2ia_symmtry_op_inv(ia,iopr))

          do lmt1 = 1, ilmt(it)
             il1=ltp(lmt1,it)
             im1=mtp(lmt1,it)
             it1=taup(lmt1,it)
             ii=(il1-1)**2+im1

             do lmt2 = lmt1, ilmt(it)
!             do lmt2 = 1, ilmt(it)
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

! ---- modified 2023/05/29
!                            hsr_wk(ia,lmt1,lmt2,1) = &
!                                 hsr_wk(ia,lmt1,lmt2,1) + &
!                                 hsr_tmp(ja,lmt3,lmt4,1)* &
!                                 crotylm_paw(n,ii,iopr,ia)* &
!                                 crotylm_paw(m,jj,iopr,ia)
!                            hsi_wk(ia,lmt1,lmt2,1) = &
!                                 hsi_wk(ia,lmt1,lmt2,1) + &
!                                 weight *weight2 * &
!                                 hsi_tmp(ja,lmt3,lmt4,1)* &
!                                 crotylm_paw(n,ii,iopr,ia)* &
!                                 crotylm_paw(m,jj,iopr,ia)
!
                            ctmp1 = 1.0d0;   ctmp2 = weight

                            hsr_wk(ia,lmt1,lmt2,1) = &
                                 hsr_wk(ia,lmt1,lmt2,1) + &
                                 ctmp1 * &
                                 hsr_tmp(ja,lmt3,lmt4,1)* &
                                 crotylm_paw(n,ii,iopr,ia)* &
                                 crotylm_paw(m,jj,iopr,ia)
                            hsi_wk(ia,lmt1,lmt2,1) = &
                                 hsi_wk(ia,lmt1,lmt2,1) + &
                                 ctmp2 * &
                                 hsi_tmp(ja,lmt3,lmt4,1)* &
                                 crotylm_paw(n,ii,iopr,ia)* &
                                 crotylm_paw(m,jj,iopr,ia)
! ---  modified 2023/05/29

                            Do ixyz1=1, 3
                               Do ixyz2=1, 3
                                  ctmp1 = op(ixyz2, ixyz1, iopr) *weight
#ifdef MOMENT_AS_PSEUDO_VECTOR
                                  ctmp1 = ctmp1 *weight2
#endif
                                  ctmp2 = op(ixyz2, ixyz1, iopr) *weight2
!
                                  hsr_wk(ia,lmt1,lmt2,ixyz2+1) &
                                       & = hsr_wk(ia,lmt1,lmt2,ixyz2+1)  &
                                       &  + ctmp1 &
                                       &    *hsr_tmp(ja,lmt3,lmt4,ixyz1+1) &
                                       &    *crotylm_paw(n,ii,iopr,ia)  &
                                       &    *crotylm_paw(m,jj,iopr,ia)
                                  hsi_wk(ia,lmt1,lmt2,ixyz2+1) &
                                       & = hsi_wk(ia,lmt1,lmt2,ixyz2+1)  &
                                       &  + ctmp2 &
                                       &    *hsi_tmp(ja,lmt3,lmt4,ixyz1+1) &
                                       &    *crotylm_paw(n,ii,iopr,ia)  &
                                       &    *crotylm_paw(m,jj,iopr,ia)
                               End do
                            End do

                         end do! lmt4
                      end do! lmt3

                   end do! jjj
                end do! iii

             end do! lmt2
          end do! lmt1
       end do! ia
    end do! iopr

!    hsr_wk = hsr_wk/nopr;  hsi_wk = hsi_wk /nopr
    hsr_wk = hsr_wk *fi;   hsi_wk = hsi_wk *fi

    do ia=1,natm
       it=ityp(ia)
       do is =1, ndim_magmom
          do lmt2=1,ilmt(it)
             do lmt1=lmt2+1,ilmt(it)
                hsr_wk(ia,lmt1,lmt2,is) =  hsr_wk(ia,lmt2,lmt1,is)
                hsi_wk(ia,lmt1,lmt2,is) = -hsi_wk(ia,lmt2,lmt1,is)
             end do
          end do
       end do
    end do

    deallocate(hsr_tmp); deallocate(hsi_tmp)

  end subroutine m_CD_symmtrz_of_ff_noncl_C
! ========================================================================= 11.0

! =================================== modified by K. Tagami ============ 11.0
!!  subroutine add_hardpart_to_chgq_l(nfout,kspin,hsr)
  subroutine add_hardpart_to_chgq_l( nfout, kspin, hsr, singlemode )
! ===================================================================== 11.0
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

    integer, intent(in)      :: nfout, kspin

! =============================== added by K. Tagami ================ 11.0
    integer, intent(in) ::  singlemode
! =================================================================== 11.0

! ================================ modified by K. Tagami ================= 11.0
!    real(kind=DP), intent(in), dimension(natm,nlmt,nlmt,nspin):: hsr
    real(kind=DP), intent(in), dimension(natm,nlmt,nlmt,ndim_magmom):: hsr
! ======================================================================== 11.0

    real(kind=DP), pointer, dimension(:)               :: ylm
    real(kind=DP), allocatable, target, dimension(:)   :: ylm_t
    real(kind=DP), allocatable, target, dimension(:,:) :: ylm_ext

    integer :: is,it,lmt1,lmt2,n,ia,mdvdb,il1,tau1,il2,tau2,ilm3,l3,iiqitg
    real(kind=DP) :: fac !, tpos(3)

    integer :: kngp_adj, n_ialist, n_ialist0, ia_start, ia_end, n_iagroup, n_ia, ia_g
    real(kind=DP), allocatable, target, dimension(:,:) :: zfcos_x, zfsin_x

    integer, allocatable, dimension(:) :: ia_list
#ifdef _VECTOR_TUNING_
    real(kind=DP), allocatable, dimension(:,:,:) :: shdg_x ! d(n_ialist0,maxm,nqitg)
#else
    real(kind=DP), allocatable, dimension(:,:) :: ylm_red, qitg_red
    real(kind=DP), allocatable, dimension(:) :: ylm_sum
    real(kind=DP), allocatable, dimension(:,:,:) :: chgq_red
    real(kind=DP), allocatable, dimension(:,:) :: shdg  ! d(max,nqitg,nspin)
    integer ::          iqm, iqmmax
    real(kind=DP) :: zdga
#endif

    integer :: m, maxm, ip, np, iq, sw_spin
    integer, parameter :: mcritical = 4*2+1
    integer, allocatable, dimension(:) :: nqitg_sp, nqitg_sp0 !d(ntyp)
    integer, allocatable, dimension(:) :: iq2l3 ! d(nqitg)
    integer, allocatable, dimension(:,:) :: nc  ! d(maxm,nqitg)
    integer :: mc ! maxval(nc)
    integer, allocatable, dimension(:,:,:) :: nc2lmt1, nc2lmt2, nc2n ! d(mc,maxm,nqitg)
! NEC tune
    integer :: ibl1,ibl2,ibsize,ncache,iwidth

! =============================== added by K. Tagami ================ 11.0
    integer :: nspin_kt, ispi_start
! =================================================================== 11.0

    integer :: id_sname = -1
    call tstatc0_begin('add_hardpart_to_chgq_l ',id_sname,1)

! =============================== added by K. Tagami =============== 11.0
    if ( singlemode == YES ) then
       ispi_start = kspin
    else
       ispi_start = 1
    endif
! ================================================================== 11.0

    if(modnrm == EXECUT) then
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

#ifdef _VECTOR_TUNING_

       n_ialist = 1

#ifdef HIUX
       n_ialist = 4
#endif
#ifdef VPP
       n_ialist = 8
#endif
#ifdef SX
       n_ialist = 8
#endif

#ifdef DEBUG_LDOS
       if(ipri>=1)  write(nfout,'(" n_ialist = ",i8)') n_ialist
#endif
       kngp_adj = iend_kngp - ista_kngp + 1
       if(mod(kngp_adj,2) == 0) kngp_adj = kngp_adj + 1
       allocate(zfcos_x(kngp_adj,n_ialist)); zfcos_x = 0.d0
       allocate(zfsin_x(kngp_adj,n_ialist)); zfsin_x = 0.d0
       allocate(ia_list(n_ialist)); ia_list = 0
! NEC tune
!       allocate(shdg_x(n_ialist))
       allocate(shdg_x(n_ialist,maxm,nqitg))
! ================================================ by K. Tagami =======
        shdg_x = 0.0d0
! ===================================================================

       if(n**2 > nel_Ylm) then
          allocate(ylm_ext(ista_kngp:iend_kngp,nel_Ylm+1:n**2)); ylm_ext = 0.d0
       end if
       allocate(ylm_t(ista_kngp:iend_kngp)); ylm_t = 0.d0
       do ilm3 = nel_Ylm+1, n**2
          call m_pwBS_sphrp2(ilm3,rltv,ista_kngp,iend_kngp,ylm_t)
          ylm_ext(:,ilm3) = ylm_t(:)
       end do

!!$       do is = 1, kspin, af+1
          do it = 1, ntyp
             mdvdb = m_PP_include_vanderbilt_pot(it)
             if(mdvdb == SKIP) cycle

             n_ia = 0
             do ia = 1, natm
                if(ityp(ia) == it) n_ia = n_ia + 1
             end do

             if(n_ialist <=0) call phase_error_with_msg(nfout, &
             'n_ialist is illegal <<m_Charge_Density.add_hardpart_to_chgq_l>>', &
             __LINE__,__FILE__)
             n_iagroup = n_ia/n_ialist + 1
             ia_start = 1
             if(ipri >= 2) write(nfout,'(" !m_CD.add_hardpart_to_chgq_l: n_iagroup = ",i8, " ityp = ",i8)') n_iagroup,it
             do ia_g = 1, n_iagroup
                n_ialist0 = 0
                ia_list = 0
                AtomcountLoop: do ia = ia_start, natm
                   if(ityp(ia) == it) then
                      n_ialist0 = n_ialist0 + 1
                      ia_list(n_ialist0) = ia
                   end if
                   if(n_ialist0 >= n_ialist) exit AtomcountLoop
                end do AtomcountLoop
                ia_start = ia+1
                if(n_ialist0 >= 1 )then
                   if(ipri >= 2) write(nfout,'(" !m_CD.add_hardpart_to_chgq_l: ia_list = ",8i8)') (ia_list(ia),ia=1,n_ialist0)


! NEC tune ------------------------------------------------------------------->
!!$      ncache = (cachesize(3)*1024)*3/4
      ncache = (m_CtrlP_cachesize()*1024)*3/4
!!$      ncache = 1024*1024*3/4
      if(ncache == 0) then
         ibsize = iend_kngp - ista_kngp + 1
      else

      iwidth = nqitg_sp(it) - nqitg_sp0(it)
      if(n_ialist0 == 1) then
         if(kimg == 1) then ! qitg_l(i,iq)*ylm(iy)*zfsc1_1(iy)
           ibsize=ncache/(8*(2+iwidth))
         else ! qitg_l(i,iq),ylm(iy),zfsc1_1(iy),zfsc2_1(iy)
           ibsize=ncache/(8*(3+iwidth))
         endif
      else if(n_ialist0 == 2) then
         if(kimg == 1) then ! qitg_l(i,iq),ylm(iy),zfsc1_1(iy),zfsc1_2(iy)
           ibsize=ncache/(8*(3+iwidth))
         else ! qitg_l(i,iq)
              ! ylm(iy),zfsc1_1(iy),zfsc1_2(iy),zfsc2_1(iy),zfsc2_2(iy)
           ibsize=ncache/(8*(5+iwidth))
         endif
      else if(n_ialist0 == 3) then
         if(kimg == 1) then ! qitg_l(i,iq)
                            ! ylm(iy),zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy)
           ibsize=ncache/(8*(4+iwidth))
         else ! qitg_l(i,iq),ylm(iy),zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy)
              !                      zfsc2_1(iy),zfsc2_2(iy),zfsc2_3(iy)
           ibsize=ncache/(8*(7+iwidth))
         endif
      else if(n_ialist0 == 4) then
         if(kimg == 1) then ! qitg_l(i,iq)
                   ! ylm(iy),zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy)
           ibsize=ncache/(8*(5+iwidth))
         else ! qitg_l(i,iq),ylm(iy)
              ! zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy)
              ! zfsc2_1(iy),zfsc2_2(iy),zfsc2_3(iy),zfsc2_4(iy)
           ibsize=ncache/(8*(9+iwidth))
         endif
      else if(n_ialist0 == 5) then
         if(kimg == 1) then ! qitg_l(i,iq),ylm(iy)
               ! zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy),zfsc1_5(iy))
           ibsize=ncache/(8*(6+iwidth))
         else ! qitg_l(i,iq),ylm(iy)
              ! zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy),zfsc1_5(iy)
              ! zfsc2_1(iy),zfsc2_2(iy),zfsc2_3(iy),zfsc2_4(iy),zfsc2_5(iy)
           ibsize=ncache/(8*(11+iwidth))
         endif
      else if(n_ialist0 == 6) then
         if(kimg == 1) then ! qitg_l(i,iq),ylm(iy)
                 ! zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy)
                 ! zfsc1_5(iy),zfsc1_6(iy)
           ibsize=ncache/(8*(7+iwidth))
         else ! qitg_l(i,iq),ylm(iy)
              ! zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy)
              ! zfsc1_5(iy),zfsc1_6(iy))
              ! zfsc2_1(iy),zfsc2_2(iy),zfsc2_3(iy),zfsc2_4(iy)
              ! zfsc2_5(iy),zfsc2_6(iy))
           ibsize=ncache/(8*(13+iwidth))
         endif
      else if(n_ialist0 == 7) then
         if(kimg == 1) then
             !  qitg_l(i,iq),ylm(iy)
             !  zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy)
             !  zfsc1_5(iy),zfsc1_6(iy),zfsc1_7(iy)
           ibsize=ncache/(8*(8+iwidth))
         else
             !  qitg_l(i,iq),ylm(iy)
             !  zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy)
             !  zfsc1_5(iy),zfsc1_6(iy),zfsc1_7(iy)
             !  zfsc2_1(iy),zfsc2_2(iy),zfsc2_3(iy),zfsc2_4(iy)
             !  zfsc2_5(iy),zfsc2_6(iy),zfsc2_7(iy)
           ibsize=ncache/(8*(15+iwidth))
         endif
      else if(n_ialist0 >= 8) then
         if(kimg == 1) then
             !  qitg_l(i,iq),ylm(iy)
             !  zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy)
             !  zfsc1_5(iy),zfsc1_6(iy),zfsc1_7(iy),zfsc1_8(iy)
           ibsize=ncache/(8*(9+iwidth))
         else
            !   qitg_l(i,iq),ylm(iy)
            !   zfsc1_1(iy),zfsc1_2(iy),zfsc1_3(iy),zfsc1_4(iy)
            !   zfsc1_5(iy),zfsc1_6(iy),zfsc1_7(iy),zfsc1_8(iy)
            !   zfsc2_1(iy),zfsc2_2(iy),zfsc2_3(iy),zfsc2_4(iy)
            !   zfsc2_5(iy),zfsc2_6(iy),zfsc2_7(iy),zfsc2_8(iy)
           ibsize=ncache/(8*(17+iwidth))
         endif
      end if
      endif
! debug
!write(6,990) 'n_ialist0,kimg,ibsize,ista_kngp,iend_kngp,iwidth=',&
!n_ialist0,kimg,ibsize,ista_kngp,iend_kngp,iwidth
!990 format(a,i2,i2,4i8)
      call calc_phase_b(natm,pos,ia_list,n_ialist0,kgp,ngabc,ista_kngp,iend_kngp,1,kngp_adj,zfcos_x,zfsin_x)
      do ibl1=ista_kngp,iend_kngp,ibsize
        ibl2=ibl1+ibsize-1
        if(ibl2.gt.iend_kngp) ibl2=iend_kngp
! NEC tune <-------------------------------------------------------------------

!!$          do ia = 1, natm
!!$             it = ityp(ia)
!!$             mdvdb = m_PP_include_vanderbilt_pot(it)
!!$             if(mdvdb == SKIP) cycle
! NEC tune (move 1 line to up)
!                   call calc_phase_b(natm,pos,ia_list,n_ialist0,kgp,ngabc,ista_kngp,iend_kngp,1,kngp_adj,zfcos_x,zfsin_x)
!!$             call calc_phase2(natm,pos,ia,kgp,ngabc,ista_kngp,iend_kngp,zfcos,zfsin)
                   ! -(b_Elec.)  -> zfcos, zfsin

                   if(iprichargedensity >= 2) write(nfout,'(" !mCD:    it,  iq,  l3,   m,ilm3")')

! =========================== modified by K. Tagami ============= 11.0
!                do is = 1, kspin, af+1
                do is = ispi_start, kspin, af+1
! =============================================================== 11.0

! NEC tune --------------------------------------------------->
                   do iq = nqitg_sp0(it), nqitg_sp(it)
                      l3 = iq2l3(iq)
                      do m = 1, 2*l3+1
                         call sum_hsr_dot_gauntc(is,it,iq,m) ! hsr, dl2p -> shdg_x(n_ialist0)
                      end do
                   end do
! NEC tune <---------------------------------------------------

                   do iq = nqitg_sp0(it), nqitg_sp(it)
                      l3 = iq2l3(iq)
                      do m = 1, 2*l3+1
                         ilm3 = l3*l3+m
                         if(iprichargedensity >= 2) write(nfout,'(" !mCD: ",9i5)') it, iq, l3, m, ilm3

! NEC tune
!                         call sum_hsr_dot_gauntc(is,it,iq,m) ! hsr, dl2p -> shdg_x(n_ialist0)
                         if(ilm3 <= nel_Ylm) then
                            ylm => ylm_l(ista_kngp:iend_kngp,ilm3)
                         else
                            ylm => ylm_ext(ista_kngp:iend_kngp,ilm3)
                         end if
                         call add_hardpart_to_chgq_l_core4(iq) ! iq, shdg_x, exp(-iGR), qitg_l, ylm -> chgq_l
                         !    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      end do
                   end do
! NEC tune
                end do
             end do ! ibl1
                end if
             end do! ia_g
          end do! it
!!$       end do! is
       deallocate(ylm_t)
       if(allocated(ylm_ext)) deallocate(ylm_ext)
       deallocate(shdg_x)
       deallocate(ia_list,zfsin_x,zfcos_x,il3)
#else

       ncache = (m_CtrlP_cachesize()*1024)*3/4
       if(ncache == 0) then
          ibsize = iend_kngp - ista_kngp + 1
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

       allocate(zfcos(ibsize)); zfcos = 0.d0
       allocate(zfsin(ibsize)); zfsin = 0.d0
       allocate(qitg_red(ibsize,nqitg))
       allocate(ylm_red(ibsize,n**2))
       allocate(ylm_sum(ibsize))

! =========================== modified by K. Tagami =================== 11.0
!       if(kspin == 2 .and. af == 0) then
!          sw_spin = ON
!       else
!          sw_spin = OFF
!       end if
! ===================================================================== 11.0

! ============================== added by K. Tagami ================== 11.0
       if ( noncol ) then
          nspin_kt = kspin
       else
          nspin_kt = nspin /(af+1)
       endif
! ================================================================ 11.0

       iqmmax = 0
       do it = 1, ntyp
          iqm = 0
          do iq = nqitg_sp0(it), nqitg_sp(it)
             l3 = iq2l3(iq)
             do m = 1, 2*l3+1
                iqm = iqm+1
             end do
          end do
          if(iqmmax < iqm) iqmmax = iqm
       end do

! ===================== modified by K. Tagami ============== 11.0
!       if(sw_spin == ON) then
!          allocate(chgq_red(ibsize,kimg,2))
!          allocate(shdg(iqmmax,2))
!       else if(sw_spin == OFF) then
!          allocate(chgq_red(ibsize,kimg,1))
!          allocate(shdg(iqmmax,1))
!       end if
! ========================================================== 11.0
! ========================= added by K. Tagami =========== 11.0
       allocate( chgq_red( ibsize, kimg, nspin_kt ))
       allocate( shdg( iqmmax,nspin_kt ))
! ========================================================= 11.0

       do ibl1=ista_kngp,iend_kngp,ibsize
          ibl2=ibl1+ibsize-1
          if(ibl2.gt.iend_kngp) ibl2=iend_kngp
          if(ibl2.gt.kgp) ibl2 = kgp

          chgq_red = 0.d0
          call substitute_qitgred()  ! qitg_l -> qitg_red
          call substitute_ylmred() ! ylm_l, ylm_ext -> ylm_red
          do ia = 1, natm
             it = ityp(ia)
             mdvdb = m_PP_include_vanderbilt_pot(it)
             if(mdvdb == SKIP) cycle

             call calc_phase_div(ia) ! -> zfsin, zfcos

!!$             iqm = 0
!!$             do iq = nqitg_sp0(it), nqitg_sp(it)
!!$                l3 = iq2l3(iq)
!!$                do m = 1, 2*l3+1
!!$                   iqm = iqm+1
!!$                   call sum_hsr_dot_gauntc0(it,ia,iq,m,iqm) ! hsr, dl2p -> shdg
!!$                end do
!!$             end do
!!$             iqmmax = iqm

! ================================= modified by K. Tagami ========== 11.0
!             do is = 1, kspin, af+1
             do is = ispi_start, kspin, af+1
! ================================================================== 11.0

                iqm = 0
                do iq = nqitg_sp0(it), nqitg_sp(it)
                   l3 = iq2l3(iq)
                   ylm_sum = 0.d0
                   do m = 1, 2*l3+1
                      ilm3 = l3*l3+m
                      iqm = iqm+1
                      call sum_hsr_dot_gauntc0(it,ia,iq,m,iqm) ! hsr, dl2p -> shdg
                      ylm_sum(:) = ylm_sum(:) + shdg(iqm,is)*ylm_red(:,ilm3)
                   end do
                   if(mod(l3,2) == 0) then
                      zdga = real(zi**(-l3))
                      call add_hardpart_to_chgq_l_div0(zdga,iq) ! iq, shdg_x, exp(-iGR), qitg_l, ylm -> chgq_l
                   else
                      zdga = aimag(zi**(-l3))
                      call add_hardpart_to_chgq_l_div1(zdga,iq) ! iq, shdg_x, exp(-iGR), qitg_l, ylm -> chgq_l
                   end if
!!$                   call add_hardpart_to_chgq_l_div(iq) ! iq, shdg_x, exp(-iGR), qitg_l, ylm -> chgq_l
                   !    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                end do
             end do
          end do
          call cp_chgqred2chgq()
       end do
       deallocate(ylm_sum,ylm_red,qitg_red,shdg)
       deallocate(chgq_red)
       deallocate(zfsin,zfcos,il3)
#endif

       deallocate(nc2n,nc2lmt2,nc2lmt1,nc,iq2l3,nqitg_sp,nqitg_sp0)
    end if

    if(iprichargedensity >= 2) then
       write(nfout,*) ' -- before average --'
! =============================================== modified by K. Tagami ======== 11.0
!       call m_CD_wd_chgq_l_small_portion(nfout)
       if ( noncol ) then
          call m_CD_wd_chgq_l_portion_noncl(nfout)
       else
         call m_CD_wd_chgq_l_small_portion(nfout)
      endif
! ============================================================================= 11.0
    endif

    if(nbztyp >= SIMPLE_CUBIC .or. af /=0 ) then
       allocate(work(kgp,kimg))
! =========================================== Added by K. Tagami ====
        work = 0.0d0
! ===================================================================

! ==================================== modified by K. Tagami ============== 11.0
!!       if(nbztyp >= SIMPLE_CUBIC)  call charge_average(NO,chgq_l)
!
      if(nbztyp >= SIMPLE_CUBIC) then
         if ( noncol ) then
!!            call charge_average_noncl( NO,chgq_l )
#ifndef USE_CHGAVG_NONCL3
            call charge_average_noncl2( NO,chgq_l )
#else
            call charge_average_noncl3( NO,chgq_l )
#endif
         else
            call charge_average( NO,chgq_l )
         endif
      endif
! ========================================================================== 11.0

       if(af /= 0) then
          call charge_average(ANTIFERRO,chgq_l)
          if(istress == 1 .or. sw_fine_STM_simulation == ON) then
             call charge_average(NO,chgsoft)         ! average of chgsoft
             call charge_average(ANTIFERRO,chgsoft)
          end if
       else if(sw_fine_STM_simulation == ON) then
! ==================================== modified by K. Tagami ============== 11.0
!!          call charge_average(NO,chgsoft)         ! average of chgsoft
!
         if ( noncol ) then
!!            call charge_average_noncl(NO,chgsoft)         ! average of chgsoft
#ifndef USE_CHGAVG_NONCL3
            call charge_average_noncl2(NO,chgsoft)
#else
            call charge_average_noncl3(NO,chgsoft)
#endif
         else
            call charge_average(NO,chgsoft)         ! average of chgsoft
         endif
! ========================================================================= 11.0
       endif

! ==================================== modified by K. Tagami ============== 11.0
!       if(af==0 .and. sw_fine_STM_simulation /= ON .and. flg_paw) &
!                            call charge_average(NO,chgsoft)
!
      if (af==0 .and. sw_fine_STM_simulation /= ON .and. flg_paw) then
         if ( noncol ) then
!!            call charge_average_noncl(NO,chgsoft)
#ifndef USE_CHGAVG_NONCL3
            call charge_average_noncl2(NO,chgsoft)
#else
            call charge_average_noncl3(NO,chgsoft)
#endif
         else
            call charge_average(NO,chgsoft)
         endif
      endif
! ========================================================================= 11.0

       deallocate(work)
    end if

    call tstatc0_end(id_sname)
  contains

#ifdef _VECTOR_TUNING_
    subroutine sum_hsr_dot_gauntc(is,it,iq,m)
      integer, intent(in) :: is,it,iq,m
      integer :: ip, lmt1, lmt2, np
      real(kind=DP) :: fac
! NEC tune
!      shdg_x(1:n_ialist0) = 0.d0
      shdg_x(1:n_ialist0,m,iq) = 0.d0
      do ip = 1, nc(m,iq)
         lmt1 = nc2lmt1(ip,m,iq)
         lmt2 = nc2lmt2(ip,m,iq)
         np = nc2n(ip,m,iq)
         fac = 2.d0; if(lmt1 == lmt2) fac = 1.d0
         do ia = 1, n_ialist0
! NEC tune
!            shdg_x(ia) = shdg_x(ia) + &
            shdg_x(ia,m,iq) = shdg_x(ia,m,iq) + &
                 & fac*iwei(ia_list(ia))*hsr(ia_list(ia),lmt1,lmt2,is)*dl2p(lmt1,lmt2,np,it)
         end do
      end do
    end subroutine sum_hsr_dot_gauntc
    subroutine add_hardpart_to_chgq_l_core4(iq)
      integer, intent(in) :: iq
      integer       :: i,iy, ia
!!$      real(kind=DP) :: dga, flchgq, f,f2
      real(kind=DP) :: flchgq, f,f2, flchgq_1, flchgq_2, flchgq_3, flchgq_4 &
           &                            , flchgq_5, flchgq_6, flchgq_7, flchgq_8, zdga
      real(kind=DP) :: qf1,qf2, qy
!!$      real(kind=DP), pointer, dimension(:) :: zfsc1,zfsc2
      real(kind=DP), pointer, dimension(:) :: zfsc1_1,zfsc2_1, zfsc1_2, zfsc2_2 &
           &         , zfsc1_3,zfsc2_3,zfsc1_4,zfsc2_4, zfsc1_5,zfsc2_5,zfsc1_6,zfsc2_6 &
           &         , zfsc1_7,zfsc2_7,zfsc1_8,zfsc2_8
!!$      real(kind=DP), allocatable, dimension(:) :: w_f ! d(ista_kngp:iend_kngp)
!!$      if(n_ialist == 1) allocate(w_f(ista_kngp:iend_kngp))
      if(mod(l3,2) == 0) then
         if(n_ialist0 >= 1) zdga = real(zi**(-l3))
         if(kimg == 1) then
!!$            zfsc1 => zfcos
            if(n_ialist0 >= 1) zfsc1_1 => zfcos_x(:,1)
            if(n_ialist0 >= 2) zfsc1_2 => zfcos_x(:,2)
            if(n_ialist0 >= 3) zfsc1_3 => zfcos_x(:,3)
            if(n_ialist0 >= 4) zfsc1_4 => zfcos_x(:,4)
            if(n_ialist0 >= 5) zfsc1_5 => zfcos_x(:,5)
            if(n_ialist0 >= 6) zfsc1_6 => zfcos_x(:,6)
            if(n_ialist0 >= 7) zfsc1_7 => zfcos_x(:,7)
            if(n_ialist0 >= 8) zfsc1_8 => zfcos_x(:,8)
         else
!!$            f2 = -1; zfsc1 => zfcos; zfsc2 => zfsin
            f2 = -1
            if(n_ialist0 >= 1) then
               zfsc1_1 => zfcos_x(:,1); zfsc2_1 => zfsin_x(:,1)
            end if
            if(n_ialist0 >= 2) then
               zfsc1_2 => zfcos_x(:,2); zfsc2_2 => zfsin_x(:,2)
            end if
            if(n_ialist0 >= 3) then
               zfsc1_3 => zfcos_x(:,3); zfsc2_3 => zfsin_x(:,3)
            end if
            if(n_ialist0 >= 4) then
               zfsc1_4 => zfcos_x(:,4); zfsc2_4 => zfsin_x(:,4)
            end if
            if(n_ialist0 >= 5) then
               zfsc1_5 => zfcos_x(:,5); zfsc2_5 => zfsin_x(:,5)
            end if
            if(n_ialist0 >= 6) then
               zfsc1_6 => zfcos_x(:,6); zfsc2_6 => zfsin_x(:,6)
            end if
            if(n_ialist0 >= 7) then
               zfsc1_7 => zfcos_x(:,7); zfsc2_7 => zfsin_x(:,7)
            end if
            if(n_ialist0 >= 8) then
               zfsc1_8 => zfcos_x(:,8); zfsc2_8 => zfsin_x(:,8)
            end if
         end if
      else
!!$         flchgq = fac*aimag(zi**(-l3))*dga*hsr(ia,lmt1,lmt2,is)
         if(n_ialist0 >= 1) zdga = aimag(zi**(-l3))
         if(kimg == 1) then
!!$            zfsc1 => zfsin
            if(n_ialist0 >= 1) zfsc1_1 => zfsin_x(:,1)
            if(n_ialist0 >= 2) zfsc1_2 => zfsin_x(:,2)
            if(n_ialist0 >= 3) zfsc1_3 => zfsin_x(:,3)
            if(n_ialist0 >= 4) zfsc1_4 => zfsin_x(:,4)
            if(n_ialist0 >= 5) zfsc1_5 => zfsin_x(:,5)
            if(n_ialist0 >= 6) zfsc1_6 => zfsin_x(:,6)
            if(n_ialist0 >= 7) zfsc1_7 => zfsin_x(:,7)
            if(n_ialist0 >= 8) zfsc1_8 => zfsin_x(:,8)
         else
!!$            f2 = 1; zfsc1 => zfsin; zfsc2 => zfcos
            f2 = 1
            if(n_ialist0 >= 1) then
               zfsc1_1 => zfsin_x(:,1); zfsc2_1 => zfcos_x(:,1)
            end if
            if(n_ialist0 >= 2) then
               zfsc1_2 => zfsin_x(:,2); zfsc2_2 => zfcos_x(:,2)
            end if
            if(n_ialist0 >= 3) then
               zfsc1_3 => zfsin_x(:,3); zfsc2_3 => zfcos_x(:,3)
            end if
            if(n_ialist0 >= 4) then
               zfsc1_4 => zfsin_x(:,4); zfsc2_4 => zfcos_x(:,4)
            end if
            if(n_ialist0 >= 5) then
               zfsc1_5 => zfsin_x(:,5); zfsc2_5 => zfcos_x(:,5)
            end if
            if(n_ialist0 >= 6) then
               zfsc1_6 => zfsin_x(:,6); zfsc2_6 => zfcos_x(:,6)
            end if
            if(n_ialist0 >= 7) then
               zfsc1_7 => zfsin_x(:,7); zfsc2_7 => zfcos_x(:,7)
            end if
            if(n_ialist0 >= 8) then
               zfsc1_8 => zfsin_x(:,8); zfsc2_8 => zfcos_x(:,8)
            end if
         end if
      end if
! NEC tune ---------------------------------------->
!      if(n_ialist0 >= 1) flchgq_1 = zdga*shdg_x(1)
!      if(n_ialist0 >= 2) flchgq_2 = zdga*shdg_x(2)
!      if(n_ialist0 >= 3) flchgq_3 = zdga*shdg_x(3)
!      if(n_ialist0 >= 4) flchgq_4 = zdga*shdg_x(4)
!      if(n_ialist0 >= 5) flchgq_5 = zdga*shdg_x(5)
!      if(n_ialist0 >= 6) flchgq_6 = zdga*shdg_x(6)
!      if(n_ialist0 >= 7) flchgq_7 = zdga*shdg_x(7)
!      if(n_ialist0 >= 8) flchgq_8 = zdga*shdg_x(8)
! NEC tune -----------------------------------------
      if(n_ialist0 >= 1) flchgq_1 = zdga*shdg_x(1,m,iq)
      if(n_ialist0 >= 2) flchgq_2 = zdga*shdg_x(2,m,iq)
      if(n_ialist0 >= 3) flchgq_3 = zdga*shdg_x(3,m,iq)
      if(n_ialist0 >= 4) flchgq_4 = zdga*shdg_x(4,m,iq)
      if(n_ialist0 >= 5) flchgq_5 = zdga*shdg_x(5,m,iq)
      if(n_ialist0 >= 6) flchgq_6 = zdga*shdg_x(6,m,iq)
      if(n_ialist0 >= 7) flchgq_7 = zdga*shdg_x(7,m,iq)
      if(n_ialist0 >= 8) flchgq_8 = zdga*shdg_x(8,m,iq)
! NEC tune <----------------------------------------

      if(n_ialist0 == 1) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp     !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               chgq_l(i,1,is) = chgq_l(i,1,is) &
                    & +flchgq_1*qitg_l(i,iq)*ylm(iy)*zfsc1_1(iy)
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp   !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               f = flchgq_1*qitg_l(i,iq)*ylm(iy)
               chgq_l(i,1,is) = chgq_l(i,1,is) + f * zfsc1_1(iy)
               chgq_l(i,2,is) = chgq_l(i,2,is) + f2 * f * zfsc2_1(iy)
            end do
         end if
      else if(n_ialist0 == 2) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp     !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               chgq_l(i,1,is) = chgq_l(i,1,is) + qitg_l(i,iq)*ylm(iy) &
                    &         *( flchgq_1*zfsc1_1(iy) + flchgq_2*zfsc1_2(iy) )
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp   !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               chgq_l(i,1,is) = chgq_l(i,1,is) + qitg_l(i,iq)*ylm(iy) &
                    &         *( flchgq_1*zfsc1_1(iy) + flchgq_2*zfsc1_2(iy) )
               chgq_l(i,2,is) = chgq_l(i,2,is) + f2 * qitg_l(i,iq)*ylm(iy) &
                    &         *( flchgq_1*zfsc2_1(iy) + flchgq_2*zfsc2_2(iy) )
            end do
         end if
      else if(n_ialist0 == 3) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp     !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               chgq_l(i,1,is) = chgq_l(i,1,is) + qitg_l(i,iq)*ylm(iy) &
                    &     *( flchgq_1*zfsc1_1(iy) + flchgq_2*zfsc1_2(iy) + flchgq_3*zfsc1_3(iy) )
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp   !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               chgq_l(i,1,is) = chgq_l(i,1,is) + qitg_l(i,iq)*ylm(iy) &
                    &     *( flchgq_1*zfsc1_1(iy) + flchgq_2*zfsc1_2(iy) + flchgq_3*zfsc1_3(iy))
               chgq_l(i,2,is) = chgq_l(i,2,is) + f2 * qitg_l(i,iq)*ylm(iy) &
                    &     *( flchgq_1*zfsc2_1(iy) + flchgq_2*zfsc2_2(iy) + flchgq_3*zfsc2_3(iy))
            end do
         end if
      else if(n_ialist0 == 4) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp     !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               chgq_l(i,1,is) = chgq_l(i,1,is) + qitg_l(i,iq)*ylm(iy) &
                    &   *(flchgq_1*zfsc1_1(iy) + flchgq_2*zfsc1_2(iy) + flchgq_3*zfsc1_3(iy)+ flchgq_4*zfsc1_4(iy))
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp   !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               chgq_l(i,1,is) = chgq_l(i,1,is) + qitg_l(i,iq)*ylm(iy) &
                    &   *( flchgq_1*zfsc1_1(iy) + flchgq_2*zfsc1_2(iy) + flchgq_3*zfsc1_3(iy) + flchgq_4*zfsc1_4(iy))
               chgq_l(i,2,is) = chgq_l(i,2,is) + f2 * qitg_l(i,iq)*ylm(iy) &
                    &   *( flchgq_1*zfsc2_1(iy) + flchgq_2*zfsc2_2(iy) + flchgq_3*zfsc2_3(iy) + flchgq_4*zfsc2_4(iy))
            end do
         end if
      else if(n_ialist0 == 5) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp     !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               chgq_l(i,1,is) = chgq_l(i,1,is) + qitg_l(i,iq)*ylm(iy) &
                    &   *(flchgq_1*zfsc1_1(iy) + flchgq_2*zfsc1_2(iy) + flchgq_3*zfsc1_3(iy)+ flchgq_4*zfsc1_4(iy) &
                    &    +flchgq_5*zfsc1_5(iy))
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp   !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               chgq_l(i,1,is) = chgq_l(i,1,is) + qitg_l(i,iq)*ylm(iy) &
                    &   *( flchgq_1*zfsc1_1(iy) + flchgq_2*zfsc1_2(iy) + flchgq_3*zfsc1_3(iy) + flchgq_4*zfsc1_4(iy)&
                    &     +flchgq_5*zfsc1_5(iy))
               chgq_l(i,2,is) = chgq_l(i,2,is) + f2 * qitg_l(i,iq)*ylm(iy) &
                    &   *( flchgq_1*zfsc2_1(iy) + flchgq_2*zfsc2_2(iy) + flchgq_3*zfsc2_3(iy) + flchgq_4*zfsc2_4(iy)&
                    &     +flchgq_5*zfsc2_5(iy))
            end do
         end if
      else if(n_ialist0 == 6) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp     !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               chgq_l(i,1,is) = chgq_l(i,1,is) + qitg_l(i,iq)*ylm(iy) &
                    &   *(flchgq_1*zfsc1_1(iy) + flchgq_2*zfsc1_2(iy) + flchgq_3*zfsc1_3(iy)+ flchgq_4*zfsc1_4(iy) &
                    &    +flchgq_5*zfsc1_5(iy) + flchgq_6*zfsc1_6(iy))
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp   !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               chgq_l(i,1,is) = chgq_l(i,1,is) + qitg_l(i,iq)*ylm(iy) &
                    &   *( flchgq_1*zfsc1_1(iy) + flchgq_2*zfsc1_2(iy) + flchgq_3*zfsc1_3(iy) + flchgq_4*zfsc1_4(iy)&
                    &     +flchgq_5*zfsc1_5(iy) + flchgq_6*zfsc1_6(iy))
               chgq_l(i,2,is) = chgq_l(i,2,is) + f2 * qitg_l(i,iq)*ylm(iy) &
                    &   *( flchgq_1*zfsc2_1(iy) + flchgq_2*zfsc2_2(iy) + flchgq_3*zfsc2_3(iy) + flchgq_4*zfsc2_4(iy)&
                    &     +flchgq_5*zfsc2_5(iy) + flchgq_6*zfsc2_6(iy))
            end do
         end if
      else if(n_ialist0 == 7) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp     !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               chgq_l(i,1,is) = chgq_l(i,1,is) + qitg_l(i,iq)*ylm(iy) &
                    &   *(flchgq_1*zfsc1_1(iy) + flchgq_2*zfsc1_2(iy) + flchgq_3*zfsc1_3(iy) + flchgq_4*zfsc1_4(iy) &
                    &    +flchgq_5*zfsc1_5(iy) + flchgq_6*zfsc1_6(iy) + flchgq_7*zfsc1_7(iy))
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp   !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               chgq_l(i,1,is) = chgq_l(i,1,is) + qitg_l(i,iq)*ylm(iy) &
                    &   *( flchgq_1*zfsc1_1(iy) + flchgq_2*zfsc1_2(iy) + flchgq_3*zfsc1_3(iy) + flchgq_4*zfsc1_4(iy)&
                    &     +flchgq_5*zfsc1_5(iy) + flchgq_6*zfsc1_6(iy) + flchgq_7*zfsc1_7(iy))
               chgq_l(i,2,is) = chgq_l(i,2,is) + f2 * qitg_l(i,iq)*ylm(iy) &
                    &   *( flchgq_1*zfsc2_1(iy) + flchgq_2*zfsc2_2(iy) + flchgq_3*zfsc2_3(iy) + flchgq_4*zfsc2_4(iy)&
                    &     +flchgq_5*zfsc2_5(iy) + flchgq_6*zfsc2_6(iy) + flchgq_7*zfsc2_7(iy))
            end do
         end if
      else if(n_ialist0 >= 8) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp     !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               chgq_l(i,1,is) = chgq_l(i,1,is) + qitg_l(i,iq)*ylm(iy) &
                    &   *(flchgq_1*zfsc1_1(iy) + flchgq_2*zfsc1_2(iy) + flchgq_3*zfsc1_3(iy) + flchgq_4*zfsc1_4(iy) &
                    &    +flchgq_5*zfsc1_5(iy) + flchgq_6*zfsc1_6(iy) + flchgq_7*zfsc1_7(iy) + flchgq_8*zfsc1_8(iy))
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
#ifdef NEC_TUNE_HARD
!CDIR PARALLEL DO PRIVATE(i,iy)
#endif
! NEC tune
!            do i = ista_kngp, iend_kngp   !for mpi
            do i = ibl1, ibl2
               iy = i - ista_kngp+1
               chgq_l(i,1,is) = chgq_l(i,1,is) + qitg_l(i,iq)*ylm(iy) &
                    &   *( flchgq_1*zfsc1_1(iy) + flchgq_2*zfsc1_2(iy) + flchgq_3*zfsc1_3(iy) + flchgq_4*zfsc1_4(iy)&
                    &     +flchgq_5*zfsc1_5(iy) + flchgq_6*zfsc1_6(iy) + flchgq_7*zfsc1_7(iy) + flchgq_8*zfsc1_8(iy))
               chgq_l(i,2,is) = chgq_l(i,2,is) + f2 * qitg_l(i,iq)*ylm(iy) &
                    &   *( flchgq_1*zfsc2_1(iy) + flchgq_2*zfsc2_2(iy) + flchgq_3*zfsc2_3(iy) + flchgq_4*zfsc2_4(iy)&
                    &     +flchgq_5*zfsc2_5(iy) + flchgq_6*zfsc2_6(iy) + flchgq_7*zfsc2_7(iy) + flchgq_8*zfsc2_8(iy))
            end do
         end if
      end if
!!$      if(n_ialist == 1) deallocate(w_f)
    end subroutine add_hardpart_to_chgq_l_core4
#else
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
      do ilm = 1, nel_Ylm
         do i = 1, ibl2-ibl1+1
            ylm_red(i,ilm) = ylm_l(i+ibl1-1,ilm)
         end do
      end do
      if(n**2 > nel_Ylm) then
         allocate(ylm_t(ibl1:ibl2)); ylm_t = 0.d0
         do ilm = nel_ylm+1, n**2
            call m_pwBS_sphrp2(ilm,rltv,ibl1,ibl2,ylm_t)
            do i = 1, ibl2-ibl1+1
               ylm_red(i,ilm) = ylm_t(i+ibl1-1)
            end do
         end do
         deallocate(ylm_t)
      end if
    end subroutine substitute_ylmred

    subroutine sum_hsr_dot_gauntc0(it,ia,iq,m,iqm)
      integer, intent(in) :: it,ia,iq,m,iqm
      integer :: ip, lmt1, lmt2, np
      real(kind=DP) :: fac

! ========================= modified by K. Tagami ====================== 11.0
!      if(sw_spin == OFF) then
!         shdg(iqm,1) = 0.d0
!         do ip = 1, nc(m,iq)
!            lmt1 = nc2lmt1(ip,m,iq)
!            lmt2 = nc2lmt2(ip,m,iq)
!            np = nc2n(ip,m,iq)
!            fac = 2.d0; if(lmt1 == lmt2) fac = 1.d0
!            shdg(iqm,1) = shdg(iqm,1) + fac*iwei(ia)*hsr(ia,lmt1,lmt2,1)*dl2p(lmt1,lmt2,np,it)
!         end do
!      else if(sw_spin == ON) then
!         shdg(iqm,1) = 0.d0; shdg(iqm,2) = 0.d0
!         do ip = 1, nc(m,iq)
!            lmt1 = nc2lmt1(ip,m,iq)
!            lmt2 = nc2lmt2(ip,m,iq)
!            np = nc2n(ip,m,iq)
!            fac = 2.d0; if(lmt1 == lmt2) fac = 1.d0
!            shdg(iqm,1) = shdg(iqm,1) + fac*iwei(ia)*hsr(ia,lmt1,lmt2,1) &
!                           &           *dl2p(lmt1,lmt2,np,it)
!            shdg(iqm,2) = shdg(iqm,2) + fac*iwei(ia)*hsr(ia,lmt1,lmt2,2) &
!                           &          *dl2p(lmt1,lmt2,np,it)
!         end do
!      end if

!      shdg(iqm,:) = 0.d0
      shdg(iqm,is) = 0.d0

      do ip = 1, nc(m,iq)
         lmt1 = nc2lmt1(ip,m,iq);    lmt2 = nc2lmt2(ip,m,iq)
         np = nc2n(ip,m,iq)
         fac = 2.d0; if(lmt1 == lmt2) fac = 1.d0
!         shdg(iqm,:) = shdg(iqm,:) &
!              &  + fac *iwei(ia) *hsr(ia,lmt1,lmt2,:) *dl2p(lmt1,lmt2,np,it)
         shdg(iqm,is) = shdg(iqm,is) &
              &  + fac *iwei(ia) *hsr(ia,lmt1,lmt2,is) *dl2p(lmt1,lmt2,np,it)
      end do
! ============================================================== 11.0

    end subroutine sum_hsr_dot_gauntc0

    subroutine calc_phase_div(ia)
      integer, intent(in) :: ia
      real(kind=DP) :: fx, fy, fz, ph
      integer :: i, iy
      fx = pos(ia,1)*PAI2
      fy = pos(ia,2)*PAI2
      fz = pos(ia,3)*PAI2
      do i = 1, ibl2-ibl1+1
         iy = i + ibl1 - 1
         ph = ngabc(iy,1)*fx+ngabc(iy,2)*fy+ngabc(iy,3)*fz
         zfcos(i) = dcos(ph)
         zfsin(i) = dsin(ph)
      end do
    end subroutine calc_phase_div

    subroutine add_hardpart_to_chgq_l_div(iq)
      integer, intent(in) :: iq
      integer       :: i,iy, iy2
      real(kind=DP) :: flchgq, f,f2,  zdga, flchgq_up, flchgq_dw, f_up, f_dw
      real(kind=DP), pointer, dimension(:) :: zfsc1,zfsc2

      if(mod(l3,2) == 0) then
         zdga = real(zi**(-l3))
         if(kimg == 1) then
            zfsc1 => zfcos
         else
            f2 = -1
            zfsc1 => zfcos; zfsc2 => zfsin
         end if
      else
         zdga = aimag(zi**(-l3))
         if(kimg == 1) then
            zfsc1 => zfsin
         else
            f2 = 1
            zfsc1 => zfsin; zfsc2 => zfcos
         end if
      end if

!!$      if(sw_spin == ON) then
!!$         flchgq_up = zdga*shdg(iqm,1)
!!$         flchgq_dw = zdga*shdg(iqm,2)
!!$      else if(sw_spin == OFF) then
!!$         flchgq = zdga*shdg(iqm,1)
!!$      end if

      if(kimg == 1) then
!!$         if(sw_spin == ON) then
!!$            do i = 1, ibl2-ibl1+1
!!$               f = qitg_red(i,iq)*ylm_sum(i)*zfsc1(i)
!!$               chgq_red(i,1,1) = chgq_red(i,1,1) +zdga*f
!!$               chgq_red(i,1,2) = chgq_red(i,1,2) +zdga*f
!!$            end do
!!$         else if(sw_spin == OFF) then
         do i = 1, ibl2-ibl1+1
            chgq_red(i,1,is) = chgq_red(i,1,is) +zdga*qitg_red(i,iq)*ylm_sum(i)*zfsc1(i)
         end do
!!$         end if
      else
!!$         if(sw_spin == ON) then
!!$            do i = 1, ibl2-ibl1+1
!!$               f = qitg_red(i,iq)*ylm_red(i,ilm3)
!!$               chgq_red(i,1,1) = chgq_red(i,1,1) +      flchgq_up*f * zfsc1(i)
!!$               chgq_red(i,2,1) = chgq_red(i,2,1) + f2 * flchgq_up*f * zfsc2(i)
!!$               chgq_red(i,1,2) = chgq_red(i,1,2) +      flchgq_dw*f * zfsc1(i)
!!$               chgq_red(i,2,2) = chgq_red(i,2,2) + f2 * flchgq_dw*f * zfsc2(i)
!!$            end do
!!$         else if(sw_spin == OFF) then
         do i = 1, ibl2-ibl1+1
            f = zdga*qitg_red(i,iq)*ylm_sum(i)
            chgq_red(i,1,is) = chgq_red(i,1,is) +      f * zfsc1(i)
            chgq_red(i,2,is) = chgq_red(i,2,is) + f2 * f * zfsc2(i)
         end do
!!$         end if
      end if
    end subroutine add_hardpart_to_chgq_l_div

    subroutine add_hardpart_to_chgq_l_div0(zdga,iq)
      real(kind=DP), intent(in) :: zdga
      integer, intent(in) :: iq
      integer       :: i
      real(kind=DP) :: f

      if(kimg == 1) then
         do i = 1, ibl2-ibl1+1
            chgq_red(i,1,is) = chgq_red(i,1,is) +zdga*qitg_red(i,iq)*ylm_sum(i)*zfcos(i)
         end do
      else
         do i = 1, ibl2-ibl1+1
            f = zdga*qitg_red(i,iq)*ylm_sum(i)
            chgq_red(i,1,is) = chgq_red(i,1,is) + f * zfcos(i)
            chgq_red(i,2,is) = chgq_red(i,2,is) - f * zfsin(i)
         end do
      end if
    end subroutine add_hardpart_to_chgq_l_div0

    subroutine add_hardpart_to_chgq_l_div1(zdga,iq)
      real(kind=DP), intent(in) :: zdga
      integer, intent(in) :: iq
      integer       :: i
      real(kind=DP) :: f

      if(kimg == 1) then
         do i = 1, ibl2-ibl1+1
            chgq_red(i,1,is) = chgq_red(i,1,is) +zdga*qitg_red(i,iq)*ylm_sum(i)*zfsin(i)
         end do
      else
         do i = 1, ibl2-ibl1+1
            f = zdga*qitg_red(i,iq)*ylm_sum(i)
            chgq_red(i,1,is) = chgq_red(i,1,is) +  f * zfsin(i)
            chgq_red(i,2,is) = chgq_red(i,2,is) +  f * zfcos(i)
         end do
      end if
    end subroutine add_hardpart_to_chgq_l_div1

    subroutine cp_chgqred2chgq()
      integer :: i, ir

! =============================== modified by K. Tagami ============= 11.0
!      if(sw_spin == ON) then
!         do ir = 1, kimg
!            do i = ibl1, ibl2
!               chgq_l(i,ir,1) = chgq_l(i,ir,1) + chgq_red(i-ibl1+1,ir,1)
!               chgq_l(i,ir,2) = chgq_l(i,ir,2) + chgq_red(i-ibl1+1,ir,2)
!            end do
!         end do
!      else if(sw_spin == OFF) then
!         do ir = 1, kimg
!            do i = ibl1, ibl2
!               chgq_l(i,ir,1) = chgq_l(i,ir,1) + chgq_red(i-ibl1+1,ir,1)
!            end do
!         end do
!      end if

      do ir = 1, kimg
         do i = ibl1, ibl2
!            chgq_l(i,ir,:) = chgq_l(i,ir,:) + chgq_red(i-ibl1+1,ir,:)
            chgq_l(i,ir,1:nspin_kt) = chgq_l(i,ir,1:nspin_kt) &
                 &                  + chgq_red(i-ibl1+1,ir,1:nspin_kt)   !ASMS
         end do
      end do
! ==================================================================== 11.0

    end subroutine cp_chgqred2chgq
#endif
  end subroutine add_hardpart_to_chgq_l

  subroutine m_CD_softpart(nfout,kv3)
    integer, intent(in) :: nfout, kv3
    integer ispin, ib1, ik, i, ip, max_elements, icolumn, istart, iend, icycle, ic
    integer :: id_sname = -1
    real(kind=DP), allocatable, dimension(:) :: wf_phase
    real(kind=DP) :: occupation

!!$    real(kind=DP), allocatable, dimension(:,:) :: occup_keep
#ifdef NEC_TUNE_SOFT
    real(kind=DP), dimension(nfft) :: bfft
    real(kind=DP), allocatable, dimension(:,:) :: tmp_afft,tmp_bfft
!!$    real(kind=DP) :: occupation
    integer i,ismp
#endif
    call tstatc0_begin('m_CD_softpart ',id_sname,1)
    if(allocated(afft)) deallocate(afft)
#ifdef NEC_TUNE_SOFT
! ==================================== Modified by K. Tagami ==========
!    allocate(afft(nfft)); call m_FFT_alloc_WF_work()
    allocate(afft(nfft)); afft = 0.0d0; call m_FFT_alloc_WF_work()
! ======================================================================
#else
! ==================================== Modified by K. Tagami ==========
!    allocate(afft(nfft)); allocate(bfft(nfft)); call m_FFT_alloc_WF_work()
!    allocate(afft(nfft)); allocate(bfft(nfft)); afft =0.0d0; bfft= 0.0d0
    allocate(afft_spin(nfft,nrank_s)); allocate(bfft(nfft)); afft_spin =0.0d0; bfft= 0.0d0
    call m_FFT_alloc_WF_work()
! ==========================================================================
#endif

    if(ipriwf >= 3) then
       max_elements = 200
       icolumn = 10
       allocate(wf_phase(nfft/2));        wf_phase = 0.0d0  !  ==== Modified by K. Tagami ==========
       icycle = ceiling(dble(min(max_elements,nfft/2))/icolumn)
    end if

#ifdef NEC_TUNE_SOFT
    if(nfft*itask < 67108865) then  ! less than 1GB
    allocate(tmp_afft(nfft,itask),tmp_bfft(nfft,itask));  tmp_afft = 0.0d0; tmp_bfft = 0.0d0 ! ===== Added by K. Tagami ====
    chgq_l = 0.d0
!    do ispin = 1, nspin, af + 1
    do ispin = ista_spin, iend_spin, (af+1)
       afft = 0.d0
       tmp_afft = 0.d0
!CDIR PARALLEL DO PRIVATE(i,ik,ib1,occupation)
       do ismp=1,itask
       do ib1 = ista_e_smp(ismp), iend_e_smp(ismp), istep_e     ! MPI
          do ik = ispin, kv3+ispin-nspin, nspin
             if(map_k(ik) /= myrank_k) cycle! MPI
             occupation = occup_l(map_z(ib1),ik)
             if(abs(occupation) < DELTA) cycle
             call m_ES_WF_in_Rspace(ik,ib1,tmp_bfft(1,ismp)) ! (swffft)
             if(ipri >= 2 .or. ipriwf >= 3) call wd_wf_phase()
!            call add_occupied_densities()  ! -(this module) occup_l, bfft -> afft
!### inline add_occupied_densities ###
             do i = 1, nfft-1, 2
                tmp_afft(i,ismp) = tmp_afft(i,ismp) + occupation*(tmp_bfft(i,ismp)**2+tmp_bfft(i+1,ismp)**2)
             end do
!#####################################
          end do
       end do
       end do
!CDIR NOCONCUR
       do ismp = 1, itask
!CDIR INNER
          do i = 1, nfft-1, 2
             afft(i)=afft(i)+tmp_afft(i,ismp)
          end do
       end do
       if(npes >= 2) then
       !   call mpi_allreduce(afft,bfft,nfft,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr) ! MPI
          call mpi_allreduce(afft,bfft,nfft,mpi_double_precision,mpi_sum,mpi_spin_group,ierr) ! MPI
          afft = bfft                          ! MPI
       end if
       call m_FFT_WF(ELECTRON,nfout,afft,DIRECT,OFF)
       call substitute_CD_for_chgq()
    end do
    if(istress == 1 .or. sw_fine_STM_simulation == ON) chgsoft = chgq_l
    if(iprichargedensity >= 2) call m_CD_wd_chgq_l_small_portion(nfout)
    deallocate(tmp_afft,tmp_bfft)

    else

    allocate(tmp_afft(nfft,itask))
! ============================================= by K. Tagami =========
        tmp_afft = 0.0d0
! ===================================================================
#endif
!!!!!!!!!!!!!!!!!!!!!!!
!   allocate(occup_keep(np_e,ista_k:iend_k))
!   occup_keep(:,:) = occup_l(:,:)
!   occup_l(:,:) = 1.0d0
!!!!!!!!!!!!!!!!!!!!!!!
    chgq_l = 0.d0
!    do ispin = 1, nspin, af + 1
    do ispin = ista_spin, iend_spin, (af+1)
    if(.not.chg_has_been_calculated) then
!       afft = 0.d0
       if(nrank_s==1) then
         afft_spin(:,1) = 0.d0
       endif
#ifdef NEC_TUNE_SOFT
       tmp_afft = 0.d0
!CDIR PARALLEL DO PRIVATE(i,ik,ib1,occupation,bfft)
       do ismp=1,itask
       do ib1 = ista_e_smp(ismp), iend_e_smp(ismp), istep_e     ! MPI
#else
       do ib1 = ista_e, iend_e, istep_e     ! MPI
#endif
          do ik = ispin, kv3+ispin-nspin, nspin
             if(map_k(ik) /= myrank_k) cycle! MPI
             occupation = occup_l(map_z(ib1),ik)
             if(abs(occupation) < DELTA) cycle
             call m_ES_WF_in_Rspace(ik,ib1,bfft) ! (swffft)
             if(ipri >= 2 .or. ipriwf >= 3) call wd_wf_phase()
#ifdef NEC_TUNE_SOFT
!### inline add_occupied_densities ###
      do i = 1, nfft-1, 2
         tmp_afft(i,ismp) = tmp_afft(i,ismp) + occupation*(bfft(i)**2+bfft(i+1)**2)
      end do
!#####################################
#else
 !            call add_occupied_densities()  ! -(this module) occup_l, bfft -> afft
             call add_occupied_densities_afft(afft_spin, ispin)  ! -(this module) occup_l, bfft -> afft
#endif
          end do
       end do
#ifdef NEC_TUNE_SOFT
       end do
#endif

!!!!!!!!!!!!!!!!!!
!      occup_keep = occup_keep
!      deallocate(occup_keep)
!!!!!!!!!!!!!!!!!!

#ifdef NEC_TUNE_SOFT
!CDIR NOCONCUR
      do ismp = 1, itask
!CDIR INNER
        do i = 1, nfft-1, 2
         afft(i)=afft(i)+tmp_afft(i,ismp)
        end do
      end do
#endif

       if(npes >= 2 .and. nrank_s==1) then
          call mpi_allreduce(afft_spin(:,1),bfft,nfft,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr) ! MPI
          afft_spin(:,1) = bfft(:)                          ! MPI
       end if
       if(nrank_s==1) then
         call m_FFT_WF(ELECTRON,nfout,afft_spin(:,1),DIRECT,OFF)
         if(ipri >= 2 .or. ipriwf >= 3) call wd_wf_phase()
         call substitute_CD_for_chgq_afft(afft_spin,ispin)
       endif
    else
       if(nrank_s==1) then
         call mpi_allreduce(chg_softpart(:,ispin),bfft,nfft,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr) ! MPI
         afft = bfft                          ! MPI
         call m_FFT_WF(ELECTRON,nfout,afft,DIRECT,OFF)
         if(ipri >= 2 .or. ipriwf >= 3) call wd_wf_phase()
         call substitute_CD_for_chgq()
       endif
    endif
    end do
    if(nrank_s==2) then
       if(chg_has_been_calculated) afft_spin(:,ista_spin) = chg_softpart(:,ista_spin)
       call mpi_allreduce(mpi_in_place,afft_spin,nfft*nrank_s,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr) ! MPI
       do ispin=1,nrank_s
         call m_FFT_WF(ELECTRON,nfout,afft_spin(:,ispin),DIRECT,OFF)
         call substitute_CD_for_chgq_afft(afft_spin,ispin)
       enddo
    endif
    if(ipriwf >= 3) deallocate(wf_phase)
    if(istress == 1 .or. sw_fine_STM_simulation == ON .or. flg_paw) chgsoft = chgq_l
    if(iprichargedensity >= 2) call m_CD_wd_chgq_l_small_portion(nfout)
#ifdef NEC_TUNE_SOFT
    deallocate(tmp_afft)
    endif
#endif
#ifdef NEC_TUNE_SOFT
    deallocate(afft); call m_FFT_dealloc_WF_work()
#else
!    deallocate(afft); deallocate(bfft); call m_FFT_dealloc_WF_work()
    deallocate(afft_spin); deallocate(bfft); call m_FFT_dealloc_WF_work()
#endif
    if(chg_has_been_calculated) call m_ES_dealloc_chgsoft()
    call tstatc0_end(id_sname)
  contains
    subroutine wd_wf_phase()
      ! -----------------
      if(ipri >= 2 .or. ipriwf >= 3 ) write(6,'(" !! ik = ",i8," ib1 = ",i8)') ik,ib1
      if(ipriwf >= 3) then
         ip = 0
         do i = 1, nfft-1, 2
            ip = ip + 1
            if(dabs(bfft(i)) < 1.d-12) then
               wf_phase(ip) = 0.d0
            else if( dabs(bfft(i+1)) < 1.d-12) then
               wf_phase(ip) = 0.d0
            else
               wf_phase(ip) = bfft(i+1)/bfft(i)
            end if
         end do
         istart = 1
         iend = max_elements
         do ic = 1, icycle
            iend = min(istart+icolumn-1,max_elements,nfft/2)
            write(nfout,'(" !bfft (R)   ",10d12.4)') (bfft(2*i-1),i=istart,iend)
            write(nfout,'(" !bfft (I)   ",10d12.4)') (bfft(2*i)  ,i=istart,iend)
            write(nfout,'(" !phase (nz) ",10f12.8)') (wf_phase(i),i=istart,iend)
            istart = iend+1
         end do
      end if
      ! ----------------
    end subroutine wd_wf_phase

    subroutine substitute_CD_for_chgq_afft(afft,ispin)
      real(kind=DP), intent(in),dimension(nfft,nrank_s) :: afft
      integer, intent(in) :: ispin
      integer       :: i, ri, i1, is
      real(kind=DP) :: fac
      integer       :: iend !mpi
      fac = 2.d0/(univol*kv3*product(fft_box_size_WF(1:3,1)))
      is = ispin
      if(nrank_s==1) is=1
      do ri = 1, kimg
         iend = iend_kngp
         if( iend_kngp > kg ) iend = kg
         if( ista_kngp <= iend ) then
            do i = ista_kngp, iend  !for mpi
               i1 = kimg*igf(i) + (ri - kimg)
!               chgq_l(i,ri,ispin) = afft(i1)*fac
               chgq_l(i,ri,ispin) = afft(i1,is)*fac
            end do
         endif
      end do
    end subroutine substitute_CD_for_chgq_afft

    subroutine substitute_CD_for_chgq
      integer       :: i, ri, i1
      real(kind=DP) :: fac
      integer       :: iend !mpi
      fac = 2.d0/(univol*kv3*product(fft_box_size_WF(1:3,1)))
      do ri = 1, kimg
         iend = iend_kngp
         if( iend_kngp > kg ) iend = kg
         if( ista_kngp <= iend ) then
            do i = ista_kngp, iend  !for mpi
               i1 = kimg*igf(i) + (ri - kimg)
               chgq_l(i,ri,ispin) = afft(i1)*fac
            end do
         endif
      end do
    end subroutine substitute_CD_for_chgq

#ifndef NEC_TUNE_SOFT
    subroutine add_occupied_densities_afft(afft,is)
      integer, intent(in) :: is
      real(kind=DP), intent(inout), dimension(nfft,nrank_s) :: afft
      integer  :: i,iss
      real(kind=DP) :: occupation
      occupation = occup_l(map_z(ib1),ik)
      iss = is
      if(nrank_s==1) iss = 1
      do i = 1, nfft-1, 2
         afft(i,iss) = afft(i,iss) + occupation*(bfft(i)**2+bfft(i+1)**2) ! MPI
      end do
      if(ipri >= 2) then
         write(nfout,'(" !cdsoft    ik, ib1 = ",2i8)') ik, ib1
         write(nfout,'(" !cdsoft     occupation = ",d16.8)') occupation
         write(nfout,'(" !cdsoft     afft : ",5d16.8)') (afft(i,iss),i=1,15)
      end if
    end subroutine add_occupied_densities_afft

    subroutine add_occupied_densities
      integer  :: i
      real(kind=DP) :: occupation
      occupation = occup_l(map_z(ib1),ik)
      do i = 1, nfft-1, 2
         afft(i) = afft(i) + occupation*(bfft(i)**2+bfft(i+1)**2) ! MPI
      end do
      if(ipri >= 2) then
         write(nfout,'(" !cdsoft    ik, ib1 = ",2i8)') ik, ib1
         write(nfout,'(" !cdsoft     occupation = ",d16.8)') occupation
         write(nfout,'(" !cdsoft     afft : ",5d16.8)') (afft(i),i=1,15)
      end if
    end subroutine add_occupied_densities
#endif

  end subroutine m_CD_softpart

! =================================== added by K. Tagami ================ 11.0
  subroutine m_CD_softpart_noncl( nfout,kv3 )
    integer, intent(in) :: nfout, kv3
    integer ispin, ib1, ik, i, ip, max_elements, icolumn, istart, iend, icycle, ic
    integer :: id_sname = -1
    real(kind=DP), allocatable, dimension(:) :: wf_phase
    real(kind=DP) :: occupation
!
    integer :: is, is1, is2, is_tmp
! ----------------
    real(kind=DP), allocatable :: afft_kt(:,:)
    real(kind=DP), allocatable :: bfft_kt(:,:)
    real(kind=DP), allocatable :: chgq_magmom( :,:,: )
! ----------------------- start -------------------

    call tstatc0_begin('m_CD_softpart_noncl ',id_sname,1)

    allocate(afft_kt(nfft,ndim_chgpot)); afft_kt = 0.0d0;
    allocate(bfft_kt(nfft,ndim_spinor)); bfft_kt = 0.0d0;
    call m_FFT_alloc_WF_work()

    if ( ipriwf >= 3 ) then
       max_elements = 200
       icolumn = 10
       allocate(wf_phase(nfft/2)); wf_phase = 0.0d0
       icycle = ceiling(dble(min(max_elements,nfft/2))/icolumn)
    end if

    chgq_l = 0.d0

    Do ib1 = ista_e, iend_e, istep_e     ! MPI
       Do ik = 1, kv3, ndim_spinor
          if ( map_k(ik) /= myrank_k ) cycle! MPI

          occupation = occup_l( map_z(ib1),ik )
          if ( abs(occupation) < DELTA ) cycle

          Do is=1, ndim_spinor
            call m_ES_WF_in_Rspace( ik +is-1, ib1, bfft_kt(:,is) )
          End do

          Do is1=1, ndim_spinor
             Do is2=1, ndim_spinor
                is_tmp = ( is1 -1 )*ndim_spinor + is2
                call add_occupied_density_matrix()
             End do
          End do
       End do
    End do

    bfft_kt = 0.0d0
    Do is1= 1, ndim_spinor
       Do is2=1, ndim_spinor
          is_tmp = ( is1 -1 )*ndim_spinor + is2
          if ( npes >= 2 ) then
!             call mpi_allreduce( afft_kt(:,is_tmp), bfft_kt(:,1), nfft, &
!                  &              mpi_double_precision, mpi_sum, &
!                  &              MPI_CommGroup, ierr )
             call mpi_allreduce( afft_kt(:,is_tmp), bfft_kt(:,1), nfft, &
                  &              mpi_double_precision, mpi_sum, &
                  &              mpi_spin_group, ierr )
             afft_kt(:,is_tmp) = bfft_kt(:,1)
          endif
          call m_FFT_WF( ELECTRON,nfout,afft_kt(:,is_tmp),DIRECT,OFF )
          call substitute_CD_for_chgq()
       end do
    end do
! --
    allocate( chgq_magmom( ista_kngp:iend_kngp,kimg,ndim_magmom ) )
    chgq_magmom = 0.0d0
!
    call m_ES_DensMat_To_MagMom_Gspace( chgq_l, chgq_magmom )
    chgq_l = chgq_magmom
!
!    write(*,*) 'CHG 1 ', chgq_l( 1,1,1 )
!    write(*,*) 'CHG 1 ', chgq_l( 1,1,2 )
!    write(*,*) 'CHG 1 ', chgq_l( 1,1,3 )
!    write(*,*) 'CHG 1 ', chgq_l( 1,1,4 )
!    stop

    deallocate( chgq_magmom )

! --
    if(istress == 1 .or. sw_fine_STM_simulation == ON .or. flg_paw) chgsoft = chgq_l
    if(iprichargedensity >= 2)  call m_CD_wd_chgq_l_portion_noncl(nfout)

    deallocate(afft_kt); deallocate(bfft_kt);
    call m_FFT_dealloc_WF_work()

    call tstatc0_end(id_sname)
  contains

    subroutine substitute_CD_for_chgq
      integer       :: i, ri, i1
      real(kind=DP) :: fac
      integer       :: iend !mpi

!!!      fac = 2.d0/( univol *(kv3/ndim_spinor) *product(fft_box_size_WF(1:3,1)))
      fac = 1.d0/( univol *(kv3/ndim_spinor) *product(fft_box_size_WF(1:3,1)))

      do ri = 1, kimg
         iend = iend_kngp
         if( iend_kngp > kg ) iend = kg
         if( ista_kngp <= iend ) then
            do i = ista_kngp, iend  !for mpi
               i1 = kimg*igf(i) + (ri - kimg)
               chgq_l(i,ri,is_tmp) = afft_kt(i1,is_tmp)*fac
            end do
         endif
      end do
    end subroutine substitute_CD_for_chgq

    subroutine add_occupied_density_matrix
      integer  :: i
      real(kind=DP) :: cr, ci

      do i = 1, nfft-1, 2
         cr =  bfft_kt(i,  is1) *bfft_kt(i,  is2) &
         &    +bfft_kt(i+1,is1) *bfft_kt(i+1,is2)
!
! ---------------------------------------------
         ci = -bfft_kt(i,  is1) *bfft_kt(i+1,is2) &
         &    +bfft_kt(i+1,is1) *bfft_kt(i,  is2)
!
! ----
!         ci = bfft_kt(i,  is1) *bfft_kt(i+1,is2) &
!           &  -bfft_kt(i+1,is1) *bfft_kt(i,  is2)
!
! --
         afft_kt(i,  is_tmp) = afft_kt(i,  is_tmp) + occupation *cr
         afft_kt(i+1,is_tmp) = afft_kt(i+1,is_tmp) + occupation *ci
      end do
    end subroutine add_occupied_density_matrix

  end subroutine m_CD_softpart_noncl
! =========================================================== 11.0

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

    if(ipri >= 2 .and. mype == 0) then
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

    if(ipri >= 2 .and. mype == 0) then
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
    if(mype == 0) then
       do is = 1, nspin, af+1
          total_charge = total_charge + chgq_l(1,1,is)*univol
       end do
    end if
    if(npes > 1) call mpi_bcast(total_charge,1,mpi_double_precision,0,MPI_CommGroup,ierr)
    if(ipritotalcharge >= 2) &
         & write(nfout,'(" !! total_charge = ",f15.6," <<m_CD_initial_CD_by_Gauss_func>>")') total_charge

    if(af /= 0) then
       call charge_average(ANTIFERRO,chgq_l)
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
       allocate(work(kgp,kimg))
       work = 0.0d0
       call charge_average(ANTIFERRO,chgq_l)
       deallocate(work)
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
            !write(*,*) 'Not supported : ndim_spinor /=2 '
            !stop
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
            !write(*,*) 'Not supported: ndim_spinor /=2'
            !stop
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
                  grt = pos(ia,1)*ngabc(i,1) + pos(ia,2)*ngabc(i,2) &
                       &                     + pos(ia,3)*ngabc(i,3)
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
                  grt = pos(ia,1)*ngabc(i,1) + pos(ia,2)*ngabc(i,2) &
                       &                     + pos(ia,3)*ngabc(i,3)
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
            !write(*,*) 'Not supported : ndim_spinor /=2 '
            !stop
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
            !write(*,*) 'Not supported: ndim_spinor /=2'
            !stop
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
            !write(*,*) 'Not supported: ndim_spinor /=2'
            !stop
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
                  grt = pos(ia,1)*ngabc(i,1) + pos(ia,2)*ngabc(i,2) &
                       &                     + pos(ia,3)*ngabc(i,3)
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
                  grt = pos(ia,1)*ngabc(i,1) + pos(ia,2)*ngabc(i,2) &
                       &                     + pos(ia,3)*ngabc(i,3)
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
    if(mype == 0) then
       do is = 1, nspin, af+1
          total_charge = total_charge + chgq_l(1,1,is)*univol
       end do
    end if
    if(npes > 1) call mpi_bcast(total_charge,1,mpi_double_precision,0,MPI_CommGroup,ierr)

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

! ===== EXP_CELLOPT ==== 2015/09/24
  subroutine m_CD_import_chgq_prev_cell(nfout,nfchgt, F_CHGT_partitioned)
    integer, intent(in) :: nfout, nfchgt
    logical, intent(in) :: F_CHGT_partitioned
    integer  :: i,j,k,is, ip
    real(kind=DP), allocatable, dimension(:,:,:) :: chgq_mpi

    real(kind=DP) :: totch_here

    integer :: id_sname = -1
    integer :: ierror

    call tstatc0_begin('m_CD_import_chgq_prev_cell ',id_sname,1)
    chgq_l = 0.0d0

    if(F_CHGT_partitioned) then
       !stop "Not supported"
       call phase_error_with_msg(nfout,'Not supported',__LINE__,__FILE__)
    else
       allocate(chgq_mpi1(kgp_prev,kimg,ndim_magmom)); chgq_mpi1 = 0.d0
       if (mype==0) then
          rewind nfchgt
          read(nfchgt, end = 9999, err = 9999) chgq_mpi1
       endif

       if (npes > 1) then
          call mpi_bcast( chgq_mpi1, kgp_prev*kimg*ndim_magmom, &
               &          mpi_double_precision, 0, MPI_CommGroup, ierr )
       endif

       do k = 1, ndim_magmom
          do j = 1, kimg
             do i = ista_kngp, iend_kngp  !for mpi
                if ( i > kgp_prev ) cycle
                chgq_l(i,j,k) = chgq_mpi1(i,j,k)
             enddo
          enddo
       enddo
       deallocate(chgq_mpi1)
    end if
    call remove_imaginary_charge

    total_charge = 0.d0
    if(mype == 0) then
       if ( noncol ) then
          do is = 1, 1
             total_charge = total_charge + chgq_l(1,1,is)*univol
          end do
       else
          do is = 1, nspin, af+1
             total_charge = total_charge + chgq_l(1,1,is)*univol
          end do
       endif
    end if
    if ( mype == 0 ) then
       chgq_l(1,:,:) = chgq_l(1,:,:) *totch /total_charge
    endif
!
    total_charge = totch

    call tstatc0_end(id_sname)

    return
9999 continue
    ierror = EOF_REACHED
    call phase_error_wo_filename(ierror, nfout, nfchgt, __LINE__, __FILE__)

  contains

    subroutine remove_imaginary_charge
      integer :: iloop, i, j, ip
      real(kind=DP) :: rinplw
      real(kind=DP), allocatable :: afft(:), afft_mpi1(:), afft_mpi2(:), afft_mpi3(:)
      real(kind=DP),parameter   :: D_min  = 1.d-40

      call m_FFT_alloc_CD_box()

      allocate(afft(ista_fftp:iend_fftp)); afft =0.0d0
      allocate(afft_mpi1(nfftp)); afft_mpi1 = 0.0d0
      if(npes >= 2) then
         allocate(afft_mpi2(mp_fftp)); afft_mpi2 = 0.0d0
         allocate(afft_mpi3(mp_fftp)); afft_mpi3 = 0.0d0
      end if

      rinplw = 1.d0 /product(fft_box_size_CD(1:3,1))

      do iloop = 1, ndim_magmom
         afft = 0.0d0;   afft_mpi1 = 0.d0

         do j = 1, kimg
            do i = ista_kngp, iend_kngp  !for mpi
               ip = (igfp_l(i)-1)*kimg + j
               !               afft_mpi1(ip) = afft_mpi1(ip) + chgq_l(i,j,iloop) !mpi
               afft_mpi1(ip) = chgq_l(i,j,iloop) !mpi
            end do
         end do

         if (npes >= 2) then
            call mpi_barrier(MPI_CommGroup,ierr)
            do j = 0, npes-1

               do i = nis_fftp(j),nie_fftp(j)
                  afft_mpi2(i-nis_fftp(j)+1) = afft_mpi1(i)
               end do
               call mpi_allreduce(afft_mpi2,afft_mpi3,mp_fftp &
                    &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)

               if(j == mype) then
                  do i = ista_fftp, iend_fftp
                     afft(i) = afft_mpi3(i - ista_fftp + 1)
                  end do
               end if
            end do
         else
            afft = afft_mpi1
         end if

         call m_FFT_CD_inverse_c(nfout,afft)        ! G-->R space
         Do i=ista_fftp, iend_fftp, 2
            afft(i+1) = 0.0d0
         End Do
!
         if ( .not. noncol ) then
            Do i=ista_fftp, iend_fftp, 2
               afft(i) = max( afft(i), D_min )
            End Do
         else if ( iloop == 1 ) then
            Do i=ista_fftp, iend_fftp, 2
               afft(i) = max( afft(i), D_min )
            End Do
         endif
         call m_FFT_CD_direct( nfout, afft )      ! R-- >G space
         !
         if(npes >= 2) then
            call mpi_allgatherv( afft, nel_fftp(mype), mpi_double_precision, &
                 &               afft_mpi1, nel_fftp, idisp_fftp, &
                 &               mpi_double_precision, MPI_CommGroup, ierr )
         else
            afft_mpi1 = afft
         end if

         do j = 1, kimg
            do i = ista_kngp, iend_kngp  !for mpi
               ip = (igfp_l(i)-1)*kimg + j
               chgq_l(i,j,iloop) = afft_mpi1(ip)
            end do
         end do
      End do
      chgq_l = chgq_l *rinplw
!
      deallocate( afft );    deallocate( afft_mpi1 )
      if ( npes >=2 ) then
         deallocate( afft_mpi2 ); deallocate( afft_mpi3 )
      endif
      call m_FFT_dealloc_CD_box()

    end subroutine remove_imaginary_charge

  end subroutine m_CD_import_chgq_prev_cell

  subroutine m_CD_calc_abs_magetization( nfout )
    use m_FFT, only : m_FFT_coef_CD_integration_kt

    integer, intent(in) :: nfout
    real(kind=DP), allocatable :: f2or1(:)

    call m_FFT_alloc_CD_box()

    allocate( f2or1(ista_fftph:iend_fftph) ); f2or1 = 0.0d0
    call m_FFT_coef_CD_integration_kt( ista_fftph, iend_fftph, f2or1 )
!
    if ( noncol ) then
       call case_noncollinear
    else
       call case_collinear
    endif

    call m_FFT_dealloc_CD_box()
    deallocate( f2or1 )

  contains

    subroutine case_noncollinear
      integer :: iloop, i, j, ip
      real(kind=DP) :: rinplw, csum, csum_mpi
      real(kind=DP), allocatable :: afft(:), afft_mpi1(:), afft_mpi2(:), afft_mpi3(:)

      allocate(afft(ista_fftp:iend_fftp)); afft =0.0d0
      allocate(afft_mpi1(nfftp)); afft_mpi1 = 0.0d0
      if(npes >= 2) then
         allocate(afft_mpi2(mp_fftp)); afft_mpi2 = 0.0d0
         allocate(afft_mpi3(mp_fftp)); afft_mpi3 = 0.0d0
      end if

      rinplw = 1.d0 /product(fft_box_size_CD(1:3,1))

      do iloop = 2, ndim_magmom
         afft = 0.0d0;   afft_mpi1 = 0.d0

         do j = 1, kimg
            do i = ista_kngp, iend_kngp  !for mpi
               ip = (igfp_l(i)-1)*kimg + j
               afft_mpi1(ip) = abs( chgq_l(i,j,iloop)  )
            end do
         end do

         if (npes >= 2) then
            call mpi_barrier(MPI_CommGroup,ierr)
            do j = 0, npes-1
               do i = nis_fftp(j),nie_fftp(j)
                  afft_mpi2(i-nis_fftp(j)+1) = afft_mpi1(i)
               end do
               call mpi_allreduce(afft_mpi2,afft_mpi3,mp_fftp &
                    &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
               if(j == mype) then
                  do i = ista_fftp, iend_fftp
                     afft(i) = afft_mpi3(i - ista_fftp + 1)
                  end do
               end if
            end do
         else
            afft = afft_mpi1
         end if

         call m_FFT_CD_inverse_c(nfout,afft)        ! G-->R space
!
         csum = 0.0d0
         Do i=ista_fftph, iend_fftph
            csum = csum + f2or1(i) *abs( afft(2*i-1) )
         End Do
         if ( npes > 1 ) then
            call mpi_allreduce( csum, csum_mpi, 1, mpi_double_precision, mpi_sum, &
                 &              MPI_CommGroup,ierr )
            csum = csum_mpi *rinplw *univol
         else
            csum = csum *rinplw *univol
         endif
!
         if ( iloop == 2 ) then
            write(nfout,'(A,F14.8)') ' !      absolute magnetization (x) = ', csum
         else if ( iloop == 3 ) then
            write(nfout,'(A,F14.8)') ' !      absolute magnetization (y) = ', csum
         else if ( iloop == 4 ) then
            write(nfout,'(A,F14.8)') ' !      absolute magnetization (z) = ', csum
         endif
      End do

      write(nfout,*)

      deallocate( afft );    deallocate( afft_mpi1 )
      if ( npes >=2 ) then
         deallocate( afft_mpi2 ); deallocate( afft_mpi3 )
      endif

    end subroutine case_noncollinear

    subroutine case_collinear
      integer :: iloop, i, j, ip
      real(kind=DP) :: rinplw, csum, csum_mpi
      real(kind=DP), allocatable :: afft(:), afft_mpi1(:), afft_mpi2(:), afft_mpi3(:)
      real(kind=DP), allocatable :: bfft(:)

      allocate(afft(ista_fftp:iend_fftp)); afft =0.0d0
      allocate(bfft(ista_fftp:iend_fftp)); bfft =0.0d0

      allocate(afft_mpi1(nfftp)); afft_mpi1 = 0.0d0
      if(npes >= 2) then
         allocate(afft_mpi2(mp_fftp)); afft_mpi2 = 0.0d0
         allocate(afft_mpi3(mp_fftp)); afft_mpi3 = 0.0d0
      end if

      rinplw = 1.d0 /product(fft_box_size_CD(1:3,1))

      Do iloop=1, nspin
         afft = 0.0d0;   afft_mpi1 = 0.d0
         do j = 1, kimg
            do i = ista_kngp, iend_kngp  !for mpi
               ip = (igfp_l(i)-1)*kimg + j
               afft_mpi1(ip) = chgq_l(i,j,iloop)
            end do
         end do

         if (npes >= 2) then
            call mpi_barrier(MPI_CommGroup,ierr)
            do j = 0, npes-1
               do i = nis_fftp(j),nie_fftp(j)
                  afft_mpi2(i-nis_fftp(j)+1) = afft_mpi1(i)
               end do
               call mpi_allreduce(afft_mpi2,afft_mpi3,mp_fftp &
                    &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
               if(j == mype) then
                  do i = ista_fftp, iend_fftp
                     afft(i) = afft_mpi3(i - ista_fftp + 1)
                  end do
               end if
            end do
         else
            afft = afft_mpi1
         end if

         call m_FFT_CD_inverse_c(nfout,afft)        ! G-->R space
         if ( iloop == 1 ) then
            bfft = bfft +afft
         else
            bfft = bfft -afft
!            bfft = bfft +afft
         endif
      End Do

      csum = 0.0d0
      if ( kimg == 1 ) then
         Do i=ista_fftph, iend_fftph
            csum = csum + f2or1(i) *abs( bfft(2*i-1) )
         End Do
      else
         Do i=ista_fftph, iend_fftph
            csum = csum + f2or1(i) *abs( bfft(2*i-1) )
         End Do
!         Do i=ista_fftp, iend_fftp, kimg
!            csum = csum +abs( bfft(i) )
!         End Do
      endif

      if ( npes > 1 ) then
         call mpi_allreduce( csum, csum_mpi, 1, mpi_double_precision, mpi_sum, &
              &              MPI_CommGroup,ierr )
         csum = csum_mpi *rinplw *univol
      else
         csum = csum *rinplw *univol
      endif
!
      write(nfout,'(A,F14.8)') ' !      absolute magnetization (z) = ', csum
      write(nfout,*)

      deallocate( afft );    deallocate( afft_mpi1 );   deallocate( bfft )
      if ( npes >=2 ) then
         deallocate( afft_mpi2 ); deallocate( afft_mpi3 )
      endif

    end subroutine case_collinear

  end subroutine m_CD_calc_abs_magetization
! ====================== 2015/09/24

! ============================ added by K. Tagami ======================= 11.0
  subroutine m_CD_rd_chgq_import_frm_collin(nfout,nfchgt, F_CHGT_partitioned)
    integer, intent(in) :: nfout, nfchgt
    logical, intent(in) :: F_CHGT_partitioned
    integer  :: i,j,k,is, ip
    real(kind=DP) :: c1

    integer :: id_sname = -1
    call tstatc0_begin('m_CD_rd_chgq_import_frm_collin ',id_sname,1)

    if(F_CHGT_partitioned) then
       write(*,*) &
            & 'Not supported : importing collinear Charge when F_CHGT_partitioned = true'
    else

       write(nfout,*) '******************************** '
       write(nfout,*) '!! Collinear spin (charge) density are used. '
       write(nfout,*) '******************************** '

       if(mype==0) rewind nfchgt
       if(npes > 1) then

          allocate(chgq_mpi1(kgp,kimg,previous_nspin_collinear)); chgq_mpi1 = 0.d0
          if(ipri >= 2) write(nfout,*) ' !D Reading chgq_g.'
          if(mype==0) read(nfchgt) chgq_mpi1

          call mpi_bcast( chgq_mpi1, kgp*kimg*previous_nspin_collinear &
               & ,mpi_double_precision,0,MPI_CommGroup,ierr)

          if ( previous_nspin_collinear == 2 ) then
             do j = 1, kimg
                do i = ista_kngp, iend_kngp  !for mpi
                   chgq_l(i,j,1) = chgq_mpi1(i,j,1) + chgq_mpi1(i,j,2 )
                   chgq_l(i,j,4) = chgq_mpi1(i,j,1) - chgq_mpi1(i,j,2 )
                end do
             enddo
          else
             do j = 1, kimg
                do i = ista_kngp, iend_kngp  !for mpi
                   chgq_l(i,j,1) = chgq_mpi1(i,j,1)
                end do
             enddo
          endif
          deallocate(chgq_mpi1)
       else

          allocate(chgq_mpi1(kgp,kimg,previous_nspin_collinear)); chgq_mpi1 = 0.d0
          read(nfchgt) chgq_mpi1
! --
          if ( previous_nspin_collinear == 2 ) then
             chgq_l(:,:,1) = chgq_mpi1(:,:,1) + chgq_mpi1(:,:,2)
             chgq_l(:,:,4) = chgq_mpi1(:,:,1) - chgq_mpi1(:,:,2)
          else
             chgq_l(:,:,1) = chgq_mpi1(:,:,1)
          endif

          deallocate(chgq_mpi1)
       end if
    end if

    if ( sw_fix_global_quantz_axis == ON ) then
       do j = 1, kimg
          do i = ista_kngp, iend_kngp  !for mpi
             c1 = chgq_l( i,j,4 )
             chgq_l( i,j,2:4 ) = c1 *Global_Quantz_Axis_Fixed(1:3)
          end do
       end do
    endif

    if(ipri >= 3) write(nfout,*) ' !D Reading chgq_g finished'
    total_charge = 0.d0

    if(mype == 0) then
       total_charge = total_charge + chgq_l(1,1,1)*univol
    end if
    if(npes > 1) call mpi_bcast(total_charge,1,mpi_double_precision,0,MPI_CommGroup,ierr)

    if(iprichargedensity >= 1) write(nfout,'("* total_charge = ",d20.8," <<m_CD_rd_chgq>>")') total_charge
    if(iprichargedensity >= 2) then
       write(nfout,'(" !D chgq_l <<m_CD_rd_chgq_import_frm_collin>>")')
       write(nfout,'(5f16.8)') (chgq_l(ista_kngp:min(ista_kngp+20,iend_kngp),1,1))
    end if
    call tstatc0_end(id_sname)

  end subroutine m_CD_rd_chgq_import_frm_collin

  subroutine m_CD_rd_chgq_noncl(nfout,nfchgt, F_CHGT_partitioned)
    integer, intent(in) :: nfout, nfchgt
    logical, intent(in) :: F_CHGT_partitioned
    integer  :: i,j,k,is, ip, ni
    real(kind=DP), allocatable, dimension(:,:,:) :: chgq_mpi

    integer :: id_sname = -1
    call tstatc0_begin('m_CD_rd_chgq_noncl ',id_sname,1)

    if(F_CHGT_partitioned) then
       rewind nfchgt
       if(npes > 1) then
          allocate(chgq_mpi(np_kngp,kimg,ndim_magmom))
          chgq_mpi = 0.0d0
          read(nfchgt) chgq_mpi

          do k = 1, ndim_magmom
             do j = 1, kimg
                do i = 1, np_kngp
                   ip = i + ista_kngp-1
                   chgq_l(ip,j,k) = chgq_mpi(i,j,k)
                end do
             end do
          end do
          deallocate(chgq_mpi)
       else
          read(nfchgt) chgq_l
       end if
    else
       if(mype==0) rewind nfchgt
       if(npes > 1) then

          allocate(chgq_mpi1(kgp,kimg,ndim_magmom)); chgq_mpi1 = 0.d0

          if(ipri >= 2) write(nfout,*) ' !D Reading chgq_g.'
          if(mype==0) read(nfchgt) chgq_mpi1

          call mpi_bcast(chgq_mpi1,kgp*kimg*ndim_magmom &
               & ,mpi_double_precision,0,MPI_CommGroup,ierr)

          do k = 1, ndim_magmom
             do j = 1, kimg
                do i = ista_kngp, iend_kngp  !for mpi
                   chgq_l(i,j,k) = chgq_mpi1(i,j,k)
                enddo
             enddo
          enddo
          deallocate(chgq_mpi1)
       else
          read(nfchgt) chgq_l
       end if
    end if

    if(ipri >= 3) write(nfout,*) ' !D Reading chgq_g finished'
    total_charge = 0.d0

    if(mype == 0) then
       ni =  1
       total_charge = total_charge + chgq_l(1,1,ni)*univol
    end if

    if(npes > 1) call mpi_bcast(total_charge,1,mpi_double_precision,0,MPI_CommGroup,ierr)

    if(iprichargedensity >= 1) write(nfout,'("* total_charge = ",d20.8," <<m_CD_rd_chgq_noncl>>")') total_charge
    if(iprichargedensity >= 2) then
       write(nfout,'(" !D chgq_l <<m_CD_rd_chgq_noncl>>")')
       write(nfout,'(5f16.8)') (chgq_l(ista_kngp:min(ista_kngp+20,iend_kngp),1,1))
    end if
    call tstatc0_end(id_sname)

  end subroutine m_CD_rd_chgq_noncl
! ======================================================================= 11.0

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
          call mpi_allreduce(chgq_mpi1,chgq_mpi2,kgp*kimg*nspin &
               &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
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

! ====================================== added by K. Tagami ============== 11.0
  subroutine m_CD_wd_chgq_noncl(nfchgt,F_CHGT_partitioned)
    integer, intent(in) :: nfchgt
    logical, intent(in) :: F_CHGT_partitioned
    real(kind=DP), allocatable, dimension(:,:,:) :: chgq_mpi
    integer :: i, ip, ri, is
    integer :: id_sname = -1
    call tstatc0_begin('m_CD_wd_chgq ',id_sname,1)

    if(F_CHGT_partitioned) then
       rewind nfchgt
       if(npes > 1) then
          allocate(chgq_mpi(np_kngp,kimg,ndim_magmom)); chgq_mpi = 0.0d0

          do is = 1, ndim_magmom
             do ri = 1, kimg
                do i = ista_kngp, iend_kngp
                   ip = i - ista_kngp+1
                   chgq_mpi(ip,ri,is) = chgq_l(i,ri,is)
                end do
             end do
          end do
          write(nfchgt) chgq_mpi
          deallocate(chgq_mpi)
       else
          write(nfchgt) chgq_l
       end if
    else
       if(mype==0) rewind nfchgt
       if(npes > 1) then

          allocate(chgq_mpi1(kgp,kimg,ndim_magmom)); chgq_mpi1 = 0.d0
          allocate(chgq_mpi2(kgp,kimg,ndim_magmom)); chgq_mpi2 = 0.d0

          chgq_mpi1(ista_kngp:iend_kngp,:,:) = chgq_l(ista_kngp:iend_kngp,:,:)

          call mpi_allreduce(chgq_mpi1,chgq_mpi2,kgp*kimg*ndim_magmom &
               &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)

          if(mype==0) write(nfchgt) chgq_mpi2
          deallocate(chgq_mpi1); deallocate(chgq_mpi2)
       else
          write(nfchgt) chgq_l
       end if
    end if
    call tstatc0_end(id_sname)
  end subroutine m_CD_wd_chgq_noncl
! ===================================================================== 11.0

  subroutine charge_average(mode,chg)
    integer,intent(in)           :: mode
    real(kind=DP), intent(inout) :: chg(ista_kngp:iend_kngp,kimg,nspin)
    integer ::       ispin, ng, no, ngp, no1, no2
    real(kind=DP) :: fi, tx,ty,tz, fp, fc, fs, zcr, zci
    real(kind=DP), pointer, dimension(:,:) :: work2
    real(kind=DP), dimension(3) :: txyz
    integer :: nopr_t, not

    allocate(work2(ista_kngp:iend_kngp,kimg)); work2 = 0.d0

    if(mode == ANTIFERRO) then
       fi = 1.d0/af
       no1 = nopr + 1; no2 = nopr + af
    else
#ifdef CHARGE_AVERAGE_WITH_SUPERCELLOPERATIONS
       if(nopr_supercell > nopr) then
          fi = 1.d0/nopr_supercell
          no1 = 1; no2 = nopr_supercell
          if(ipri>=1) write(6,'(" no1,no2 = ",2i8)') no1, no2
       else
#endif
          fi = 1.d0/nopr
          if ( charge_symm_mode >= chg_symm_level1 ) then
             fi = 1.0d0 /dble(nopr_from_fbz_to_ibz)
          endif
          no1 = 1; no2 = nopr
#ifdef CHARGE_AVERAGE_WITH_SUPERCELLOPERATIONS
       end if
#endif
    end if

    do ispin = 1, nspin, af+1
       call cp_chg_to_work(ispin,chg) ! chg -> work
       work2 = 0.d0                   ! initialization

       do not = no1, no2

#ifdef CHARGE_AVERAGE_WITH_SUPERCELLOPERATIONS
          if(mode /= ANTIFERRO .and. nopr_supercell > nopr) then
             txyz(1:3) = tau_supercell(1:3,not)*PAI2
             no = iop_supercell(not)
          else
#endif
             no = not
             txyz(1:3) = tau(1:3,no,BUCS)*PAI2
#ifdef CHARGE_AVERAGE_WITH_SUPERCELLOPERATIONS
          end if
#endif
          if(kimg == 1) then
             do ng = ista_kngp, iend_kngp !for mpi
                ngp = ngpt_l(ng,no)
                fp = ngabc(ngp,1)*txyz(1) + ngabc(ngp,2)*txyz(2) + ngabc(ngp,3)*txyz(3)
                work2(ng,1)        = work2(ng,1) + dcos(fp)*work(ngp,1)
             end do
          else if(kimg == 2) then
             do ng = ista_kngp, iend_kngp !for mpi
                ngp= ngpt_l(ng,no)
                fp = ngabc(ngp,1)*txyz(1) + ngabc(ngp,2)*txyz(2) + ngabc(ngp,3)*txyz(3)
                fc = dcos(fp);     fs = dsin(fp)
                zcr= work(ngp,1);  zci= work(ngp,kimg)
                work2(ng,1)        = work2(ng,1) + fc*zcr - fs*zci
                work2(ng,2)        = work2(ng,2) + fc*zci + fs*zcr
             end do
          end if
       end do
       if(mode /= ANTIFERRO) chg(:,:,ispin) = work2(:,:)*fi
    end do

    if(mode == ANTIFERRO) chg(:,:,nspin) = work2(:,:)*fi

    deallocate(work2)
  end subroutine charge_average

! ================================== added by K. Tagami =================== 11.0
  subroutine charge_average_noncl(mode,chg)
    integer,intent(in)           :: mode
    real(kind=DP), intent(inout) :: chg(ista_kngp:iend_kngp,kimg,ndim_magmom)
    integer ::       ng, no, ngp, no1, no2
    integer ::       is, ni
    real(kind=DP) :: fi, tx,ty,tz, fp, fc, fs, zcr, zci
    real(kind=DP), pointer, dimension(:,:) :: work2

    allocate(work2(ista_kngp:iend_kngp,kimg)); work2 = 0.d0

!    if(mode == ANTIFERRO) then
!       fi = 1.d0/af
!       no1 = nopr + 1; no2 = nopr + af
!    else
       fi = 1.d0/nopr
       no1 = 1; no2 = nopr
!    end if

     if ( charge_symm_mode >= chg_symm_level1 ) then
        fi = 1.d0 /dble(nopr_from_fbz_to_ibz)
     endif

      do ni = 1, ndim_magmom

         call cp_chg_to_work_noncl( ni,chg ) ! chg -> work
         work2 = 0.d0                   ! initialization

         do no = no1, no2
          if ( charge_symm_mode >= chg_symm_level1 ) then
            if( flg_opr_from_fbz_to_ibz(no) == 0 ) cycle
          endif

            tx = tau(1,no,BUCS)*PAI2
            ty = tau(2,no,BUCS)*PAI2
            tz = tau(3,no,BUCS)*PAI2
            if(kimg == 1) then
               do ng = ista_kngp, iend_kngp !for mpi
                  ngp = ngpt_l(ng,no)
                  fp = ngabc(ngp,1)*tx + ngabc(ngp,2)*ty + ngabc(ngp,3)*tz
                  work2(ng,1)        = work2(ng,1) + dcos(fp)*work(ngp,1)
               end do
            else if(kimg == 2) then
               do ng = ista_kngp, iend_kngp !for mpi
                  ngp= ngpt_l(ng,no)
                  fp = ngabc(ngp,1)*tx + ngabc(ngp,2)*ty + ngabc(ngp,3)*tz
                  fc = dcos(fp);     fs = dsin(fp)
                  zcr= work(ngp,1);  zci= work(ngp,kimg)
                  work2(ng,1)        = work2(ng,1) + fc*zcr - fs*zci
                  work2(ng,2)        = work2(ng,2) + fc*zci + fs*zcr
               end do
            end if
         end do
         chg(:,:,ni) = work2(:,:)*fi
      end do

    deallocate(work2)
  end subroutine charge_average_noncl

  subroutine charge_average_noncl2( mode, chg )
    integer,intent(in)           :: mode
    real(kind=DP), intent(inout) :: chg(ista_kngp:iend_kngp,kimg,ndim_magmom)
    integer ::       ng, no, ngp, no1, no2
    integer ::       is, ii
    real(kind=DP) :: tx,ty,tz, fp, fc, fs, factor
!
    real(kind=DP), allocatable :: work2(:,:)
    real(kind=DP), allocatable :: mag_work(:,:,:)
    real(kind=DP), allocatable :: mag_work2(:,:,:)
!
    real(kind=DP) :: mag_tmp(kimg,3)
!!!    real(kind=DP) :: op_in_gsp( 3,3,nopr )
! --
    allocate(work2(ista_kngp:iend_kngp,kimg)); work2 = 0.d0
    allocate(mag_work(kgp,kimg,3)); mag_work = 0.d0
    allocate(mag_work2(ista_kngp:iend_kngp,kimg,3)); mag_work2 = 0.d0

    do is = 2, ndim_magmom
       call cp_chg_to_work_noncl( is,chg ) ! chg -> work   ( mag_mom )
       mag_work(:,1:kimg,is-1) = work(:,1:kimg)
    End do

    call cp_chg_to_work_noncl( 1,chg )       ! chg -> work ( toal_chg )

! ------------------------------
    factor = 1.d0/nopr;  no1 = 1;     no2 = nopr
    if ( charge_symm_mode >= chg_symm_level1 ) then
       factor = 1.0d0 /dble(nopr_from_fbz_to_ibz)
    endif
!
! ------------------- begin ------
    do ng = ista_kngp, iend_kngp                    !for mpi
       do no = no1, no2
          if ( charge_symm_mode >= chg_symm_level1 )then
               if(flg_opr_from_fbz_to_ibz(no) == 0 ) cycle
          endif

          tx = tau(1,no,BUCS)*PAI2
          ty = tau(2,no,BUCS)*PAI2
          tz = tau(3,no,BUCS)*PAI2

          ngp= ngpt_l(ng,no)

          fp = ngabc(ngp,1)*tx + ngabc(ngp,2)*ty + ngabc(ngp,3)*tz
          fc = dcos(fp);     fs = dsin(fp)

!          Do ii=1, 3
!             mag_tmp(1:kimg,ii) = op_in_gsp(ii,1,no) *mag_work(ngp,1:kimg,1) &
!                  &            + op_in_gsp(ii,2,no) *mag_work(ngp,1:kimg,2) &
!                  &            + op_in_gsp(ii,3,no) *mag_work(ngp,1:kimg,3)
!          End Do
          Do ii=1, 3
!               ! Following three lines are revised according to a report from ASMS Co.ltd, 10 March 2016.
             mag_tmp(1:kimg,ii) = op(ii,1,int(invop(no))) *mag_work(ngp,1:kimg,1) &
                  &            + op(ii,2,int(invop(no))) *mag_work(ngp,1:kimg,2) &
                  &            + op(ii,3,int(invop(no))) *mag_work(ngp,1:kimg,3)
!!$             mag_tmp(1:kimg,ii) = op(ii,1,invop(no)) *mag_work(ngp,1:kimg,1) &
!!$                  &            + op(ii,2,invop(no)) *mag_work(ngp,1:kimg,2) &
!!$                  &            + op(ii,3,invop(no)) *mag_work(ngp,1:kimg,3)
          End Do
#ifdef MOMENT_AS_PSEUDO_VECTOR
          if ( determinant_op(invop(no)) < 0 ) mag_tmp = -mag_tmp
#endif
! == KT_add === 2014/12/29
          if ( magmom_dir_inversion_opr_flag(no) == -1 ) mag_tmp = -mag_tmp
! ============= 2014/12/29
          work2(ng,1) = work2(ng,1) +fc *work(ngp,1) -fs *work(ngp,2)
          work2(ng,2) = work2(ng,2) +fc *work(ngp,2) +fs *work(ngp,1)

          mag_work2(ng,1,1:3) = mag_work2(ng,1,1:3) &
               &               +fc *mag_tmp(1,1:3) -fs *mag_tmp(2,1:3)
          mag_work2(ng,2,1:3) = mag_work2(ng,2,1:3) &
               &               +fc *mag_tmp(2,1:3) +fs *mag_tmp(1,1:3)
       end do
    end do
!
    chg(:,:,1) = work2(:,:) *factor
    Do is=2, ndim_magmom
       chg(:,:,is) = mag_work2(:,:,is-1) *factor
    End do
!
    deallocate( mag_work ); deallocate( mag_work2 )
    deallocate( work2 )

  end subroutine charge_average_noncl2

  subroutine charge_average_noncl3( mode, chg )
    integer,intent(in)           :: mode
    real(kind=DP), intent(inout) :: chg(ista_kngp:iend_kngp,kimg,ndim_magmom)
    integer ::       ng, no, ngp, no1, no2
    integer ::       is, ii
    real(kind=DP) :: tx,ty,tz, fp, fc, fs, factor
!
    real(kind=DP), allocatable :: work2(:,:)
    real(kind=DP), allocatable :: mag_work(:,:,:)
    real(kind=DP), allocatable :: mag_work2(:,:,:)
!
    real(kind=DP) :: mag_tmp(3,kimg), mag_sum(3,kimg), chg_sum(kimg)
    integer :: flag(ista_kngp:iend_kngp)
! --
    allocate(work2(ista_kngp:iend_kngp,kimg)); work2 = 0.d0
    allocate(mag_work(3,kgp,kimg)); mag_work = 0.d0
    allocate(mag_work2(3,ista_kngp:iend_kngp,kimg)); mag_work2 = 0.d0

    do is = 2, ndim_magmom
       call cp_chg_to_work_noncl( is,chg ) ! chg -> work   ( mag_mom )
       mag_work(is-1,:,1:kimg) = work(:,1:kimg)
    End do

    call cp_chg_to_work_noncl( 1,chg )       ! chg -> work ( toal_chg )
! ------------------------------
    factor = 1.d0/nopr;  no1 = 1;     no2 = nopr
    if ( charge_symm_mode >= chg_symm_level1 ) then
       factor = 1.0d0 /dble(nopr_from_fbz_to_ibz)
    endif
!
! ------------------- begin ------
    flag = 0

    do ng = ista_kngp, iend_kngp                    !for mpi
       mag_sum = 0.0d0
       chg_sum = 0.0d0

       do no = no1, no2
          if ( charge_symm_mode >= chg_symm_level1 )then
            if(flg_opr_from_fbz_to_ibz(no) == 0 ) cycle
          endif

!          tx = tau(1,no,BUCS)*PAI2
!          ty = tau(2,no,BUCS)*PAI2
!          tz = tau(3,no,BUCS)*PAI2
          tx = tau(1,int(invop(no)),BUCS)*PAI2
          ty = tau(2,int(invop(no)),BUCS)*PAI2
          tz = tau(3,int(invop(no)),BUCS)*PAI2

          ngp= ngpt_l( ng, int(invop(no)) )
          fp = ngabc(ngp,1)*tx + ngabc(ngp,2)*ty + ngabc(ngp,3)*tz
          fc = dcos(fp);     fs = -dsin(fp)

          Do ii=1, 3
             mag_tmp(ii,1:kimg) = op(ii,1,int(invop(no))) *mag_work(1,ngp,1:kimg) &
                  &             + op(ii,2,int(invop(no))) *mag_work(2,ngp,1:kimg) &
                  &             + op(ii,3,int(invop(no))) *mag_work(3,ngp,1:kimg)
          End Do
#ifdef MOMENT_AS_PSEUDO_VECTOR
          if ( determinant_op(invop(no)) < 0 ) mag_tmp = -mag_tmp
#endif
! == KT_add === 2014/08/14
          if ( magmom_dir_inversion_opr_flag(no) == -1 ) mag_tmp = -mag_tmp
! ============= 2014/08/14
          mag_sum(:,1) = mag_sum(:,1) + fc *mag_tmp(:,1) -fs *mag_tmp(:,2)
          mag_sum(:,2) = mag_sum(:,2) + fc *mag_tmp(:,2) +fs *mag_tmp(:,1)
          chg_sum(1) = chg_sum(1) + fc *work(ngp,1) - fs *work(ngp,2)
          chg_sum(2) = chg_sum(2) + fc *work(ngp,2) + fs *work(ngp,1)
       end do

!       chg_sum = chg_sum /dble(nopr)
!       mag_sum = mag_sum /dble(nopr)
       chg_sum = chg_sum *factor
       mag_sum = mag_sum *factor

       do no = no1, no2
          tx = tau(1,no,BUCS)*PAI2
          ty = tau(2,no,BUCS)*PAI2
          tz = tau(3,no,BUCS)*PAI2

          ngp= ngpt_l( ng, no )
          if ( ngp < ista_kngp ) cycle
          if ( ngp > iend_kngp ) cycle
          if ( flag(ngp) == 1 ) cycle

          fp = ngabc(ngp,1)*tx + ngabc(ngp,2)*ty + ngabc(ngp,3)*tz
          fc = dcos(fp);     fs = dsin(fp)

          work2(ngp,1) = fc *chg_sum(1) -fs *chg_sum(2)
          work2(ngp,2) = fc *chg_sum(2) +fs *chg_sum(1)

          Do ii=1, 3
             mag_tmp(ii,1:kimg) = op(ii,1,no) *mag_sum(1,1:kimg) &
                  &             + op(ii,2,no) *mag_sum(2,1:kimg) &
                  &             + op(ii,3,no) *mag_sum(3,1:kimg)
          End Do
#ifdef MOMENT_AS_PSEUDO_VECTOR
          if ( determinant_op(no) < 0 ) mag_tmp = -mag_tmp
#endif
          mag_work2(:,ngp,1) = fc *mag_tmp(:,1) - fs *mag_tmp(:,2)
          mag_work2(:,ngp,2) = fc *mag_tmp(:,2) + fs *mag_tmp(:,1)

          if ( magmom_dir_inversion_opr_flag(no) == -1 ) then
             mag_work2(:,ngp,:) = -mag_work2(:,ngp,:)
          endif

          flag(ngp) = 1
       end do
    end do
! --
    do ng = ista_kngp, iend_kngp
       if ( flag(ng) == 0 ) then
          !write(*,*) 'error : charge_average_noncl3 ', ng
          !stop
          call phase_error_with_msg(6,'error : charge_average_noncl3 ',__LINE__,__FILE__)
       endif
    end do
!
    chg(:,:,1) = work2(:,:)
    Do is=2, ndim_magmom
       chg(:,:,is) = mag_work2(is-1,:,:)
    End do
!
    deallocate( mag_work ); deallocate( mag_work2 )
    deallocate( work2 )

  end subroutine charge_average_noncl3
! ====================================================================== 11.0

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
         !stop ' initial_charge_filetype is invalid << m_CD_initial_CD_by_file_rspace>>'
         call phase_error_with_msg(nfout,'initial_charge_filetype is invalid << m_CD_initial_CD_by_file_rspace>>' &
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
         call phase_error_with_msg(nfout, ' a set of (nlp,nmp,nnp) in the nfchr-file is invalid '//&
        '<<m_CD_initial_CD_by_file_rspace>>',__LINE__,__FILE__)

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
         call phase_error_with_msg(nfout, ' initial_charge_filetype is invalid <<'// &
         & 'm_CD_initial_CD_by_file_rspace>>',__LINE__,__FILE__)
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
         call phase_error_with_msg(nfout, ' a set of (nlp,nmp,nnp) in the nfchr-file is invalid'//&
         ' <<m_CD_initial_CD_by_file_rspace>>',__LINE__,__FILE__)
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
            if(kgp_reduced < i) cycle
#ifdef _MPIFFT_
            ip = (igfp_nonpara(i)-1)*kimg + j
#else
            ip = (igfp_l(i)-1)*kimg + j
#endif
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
              &  ,mpi_sum,MPI_CommGroup,ierr)
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
            if(kgp_reduced < i) cycle
#ifdef _MPIFFT_
            ip = (igfp_nonpara(i)-1)*kimg + j
#else
            ip = (igfp_l(i)-1)*kimg + j
#endif
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
              &  ,mpi_sum,MPI_CommGroup,ierr)
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
          if(kgp_reduced < i) cycle
#ifdef _MPIFFT_
          ip = (igfp_nonpara(i)-1)*kimg + j
#else
          ip = (igfp_l(i)-1)*kimg + j
#endif
          affttmp(ip) = affttmp(ip) + chgq_l(i,j,iloop)
       end do
    end do

    if(npes >= 2) then
!!$       call mpi_allreduce(affttmp,aff,nfftp,mpi_double_precision &
       call mpi_allreduce(affttmp,aff,nfft,mpi_double_precision &
            &  ,mpi_sum,MPI_CommGroup,ierr)
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
    call add_hardpart_to_chgq_l( nfout, is, hsr, YES )   ! (lclchg) add hardpart and make average
!              ==== modidied by K. Tagami === 11.0  ! (previous)   call add_hardpart_to_chgq_l(nfout,is,hsr)
    if(npes > 1) then
       chgq0 = 0.d0
       if(ista_kngp <= 1 .and. iend_kngp >= 1) then
          chgq0 = chgq_l(1,1,is)   !=== modified by K. Tagami === 11.0 ! (previous)  chgq0 = chgq_l(1,1,1)
       end if
       call mpi_allreduce(chgq0,z,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       chgq0 = z
    else
       chgq0 = chgq_l(1,1,is)   !=== modified by K. Tagami === 11.0 ! (previous)  chgq0 = chgq_l(1,1,1)
    end if

    call tstatc0_end(id_sname)
  end subroutine m_CD_hardpart_sub

  subroutine m_CD_hardpart_sub3(nfout,is,ik,ib,chgq0)
    integer, intent(in) :: nfout, is,ik, ib
    real(kind=DP), intent(out) :: chgq0

    real(kind=DP) :: z
    integer :: id_sname = -1
    call tstatc0_begin('m_CD_hardpart_sub ',id_sname,1)

    call summation_of_ff_sub_enl(nfout,is,ik,ib) ! -(m_C.D.) (vnlsum) fsr_l, fsi_l, occup_l --> hsr

    chgq_l = 0.d0
    call add_hardpart_to_chgq_l( nfout, is, hsr, YES )

    if(npes > 1) then
       chgq0 = 0.d0
       if(ista_kngp <= 1 .and. iend_kngp >= 1) then
          chgq0 = chgq_l(1,1,is)  ! ===== modified by K. Tagami ===== 11.0  !          chgq0 = chgq_l(1,1,1)
       end if
       call mpi_allreduce(chgq0,z,1,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       chgq0 = z
    else
       chgq0 = chgq_l(1,1,is)     ! ==== modified by K. Tagami ===== 11.0   !       chgq0 = chgq_l(1,1,1)
    end if

    call tstatc0_end(id_sname)
  end subroutine m_CD_hardpart_sub3

!  ================================= added by K. Tagami =============== 11.0
  subroutine m_CD_hardpart_sub_noncl( nfout, ik, ib, chgq0 )
    integer, intent(in) :: nfout, ik, ib
    real(kind=DP), intent(out) :: chgq0( ndim_magmom )

    real(kind=DP) :: ztmp( ndim_magmom )

    integer :: id_sname = -1

    call tstatc0_begin('m_CD_hardpart_sub_noncl ',id_sname,1)

    call summation_of_ff_sub_noncl(ik,ib)
                           ! -(m_C.D.) (vnlsum) fsr_l, fsi_l, occup_l --> hsr

    chgq_l = 0.d0
    call add_hardpart_to_chgq_l( nfout, ndim_magmom, hsr, NO )
                   ! (lclchg) add hardpart and make average
!!$    if(iprichargedensity >= 2) then
!!$       write(nfout,'(" -- hardpart summed --")')
!!$       call m_CD_wd_chgq_l_small_portion(nfout)
!!$    end if

    if (npes > 1) then
       chgq0 = 0.d0
       if(ista_kngp <= 1 .and. iend_kngp >= 1) then
          chgq0(:) = chgq_l(1,1,:)
       end if
       call mpi_allreduce( chgq0, ztmp, ndim_magmom, mpi_double_precision, &
            &              mpi_sum, MPI_CommGroup, ierr )
       chgq0 = ztmp
    else
       chgq0(:) = chgq_l(1,1,:)
    end if

    call tstatc0_end(id_sname)
  end subroutine m_CD_hardpart_sub_noncl
! ============================================================ 11.0

  subroutine m_CD_restore_chgq()
    chgq_l = chgqo_l
  end subroutine m_CD_restore_chgq

  subroutine summation_of_ff_sub(is,ik,ib)
    integer, intent(in) :: is, ik, ib

    real(kind=DP)   :: w_n
    integer         :: ia, lmt1, lmt2, p, q, it
    real(kind=DP),pointer, dimension(:,:) :: fs

    if(k_symmetry(ik) == GAMMA) then
       allocate(fs(nlmta,1))
       fs = 0.d0
       if(map_ek(ib,ik) == mype) then
          fs(1:nlmta,1) = fsr_l(map_z(ib),1:nlmta,ik)
       end if
       if(npes > 1) then
          call mpi_bcast(fs,nlmta,mpi_double_precision,map_ek(ib,ik),MPI_CommGroup,ierr)
       end if
       w_n = 2.d0
       hsr(:,:,:,is) = 0.d0
       do ia = 1, natm
          it = ityp(ia)
          do lmt1 = 1, ilmt(it)
             p = lmta(lmt1,ia)
             do lmt2 = lmt1, ilmt(it)
                q = lmta(lmt2,ia)
                hsr(ia,lmt1,lmt2,is) = hsr(ia,lmt1,lmt2,is) &
                     &  + w_n * (  fs(p,1)*fs(q,1))
             end do! lmt2
          end do! lmt1
       end do! ia
    else
       allocate(fs(nlmta,2))
       fs = 0.d0
       if(map_ek(ib,ik) == mype) then
          fs(1:nlmta,1) = fsr_l(map_z(ib),1:nlmta,ik)
          fs(1:nlmta,2) = fsi_l(map_z(ib),1:nlmta,ik)
       end if
       if(npes > 1) then
          call mpi_bcast(fs,nlmta*2,mpi_double_precision,map_ek(ib,ik),MPI_CommGroup,ierr)
       end if
       w_n = 2.d0
       hsr(:,:,:,is) = 0.d0
       do ia = 1, natm
          it = ityp(ia)
          do lmt1 = 1, ilmt(it)
             p = lmta(lmt1,ia)
             do lmt2 = lmt1, ilmt(it)
                q = lmta(lmt2,ia)
                hsr(ia,lmt1,lmt2,is) = hsr(ia,lmt1,lmt2,is) &
                     &  + w_n * (  fs(p,1)*fs(q,1)+ fs(p,2)*fs(q,2))
             end do! lmt2
          end do! lmt1
       end do! ia
    end if

!!$    if(npes >= 1) then
!!$       allocate(hsr_mpi(natm,nlmt,nlmt))
!!$       call mpi_allreduce(hsr(1,1,1,is),hsr_mpi,natm*nlmt*nlmt,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr) ! MPI
!!$       hsr(:,:,:,is) = hsr_mpi
!!$       deallocate(hsr_mpi)
!!$    end if

    deallocate(fs)
  end subroutine summation_of_ff_sub

!=================================== added by K. Tagami ================= 11.0
  subroutine summation_of_ff_sub_noncl(ik,ib)
    integer, intent(in) :: ik, ib

    real(kind=DP)   :: w_n
    integer         :: ia, is1, is2, istmp, lmt1, lmt2, p, q, it
    real(kind=DP), allocatable, dimension(:,:,:) :: fs
    real(kind=DP), allocatable, dimension(:,:,:,:) :: hsr_magmom
    real(kind=DP), allocatable, dimension(:,:,:,:) :: hsi_magmom

    allocate(fs(nlmta,2,ndim_spinor));   fs = 0.d0

    if(map_ek(ib,ik) == mype) then
       Do is1=1, ndim_spinor
          fs(1:nlmta,1,is1) = fsr_l(map_z(ib),1:nlmta,ik+is1-1)
          fs(1:nlmta,2,is1) = fsi_l(map_z(ib),1:nlmta,ik+is1-1)
       End do
    end if

    if(npes > 1) then
       call mpi_bcast( fs, nlmta*2*ndim_spinor, mpi_double_precision, &
            &          map_ek(ib,ik), MPI_CommGroup, ierr )
    end if

    w_n = 1.d0
    hsr = 0.d0; hsi = 0.0d0

    do ia = 1, natm
       it = ityp(ia)
       Do is1=1, ndim_spinor
          Do is2=1, ndim_spinor
             istmp = ( is1 -1 )*ndim_spinor + is2

             do lmt1 = 1, ilmt(it)
                p = lmta(lmt1,ia)
                do lmt2 = lmt1, ilmt(it)
                   q = lmta(lmt2,ia)

                   hsr(ia,lmt1,lmt2,istmp) = hsr(ia,lmt1,lmt2,istmp) &
                        &  + w_n * ( fs(p,1,is1)*fs(q,1,is2) &
                        &          + fs(p,2,is1)*fs(q,2,is2) )
                   hsi(ia,lmt1,lmt2,istmp) = hsi(ia,lmt1,lmt2,istmp) &
                        &  + w_n * (-fs(p,1,is1)*fs(q,2,is2) &
                        &          + fs(p,2,is1)*fs(q,1,is2) )
                end do! lmt2
             end do! lmt1
          End do
       End do
    end do! ia

    allocate( hsr_magmom( natm, nlmt, nlmt, ndim_magmom ) ); hsr_magmom = 0.0d0
    allocate( hsi_magmom( natm, nlmt, nlmt, ndim_magmom ) ); hsi_magmom = 0.0d0

    call m_ES_DensMat_To_MagMom_hsr( natm, nlmt, hsr, hsi, hsr_magmom, hsi_magmom )
    hsr = hsr_magmom; hsi = hsi_magmom

    deallocate( hsr_magmom, hsi_magmom )
    deallocate(fs)

  end subroutine summation_of_ff_sub_noncl
! ===================================================================== 11.0

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
       call mpi_allreduce(bfft,bfft_mpi,nfftp,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
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
       call phase_error_with_msg(6, '<<m_CD_map_chgqenl_to_fft_box>>',__LINE__,__FILE__)
    end if
    if(ip_max*kimg > nfftp) then
       write(6,'(" !! ip_max * kimg = ",i20," > nfftp (= ", i20," ) <<m_CD_map_chgqnl_to_fft_box>>")') ip_max*kimg, nfftp
       call phase_error_with_msg(6, '<<m_CD_map_chgqenl_to_fft_box>>',__LINE__,__FILE__)
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
       call phase_error_with_msg(6, '<<m_CD_map_chgqenl_to_fft_box>>',__LINE__,__FILE__)
    end if
    if(ip_max*kimg > nfftp) then
       write(6,'(" !! ip_max * kimg = ",i20," > nfftp (= ", i20," ) <<m_CD_map_chgqnl_to_fft_box>>")') ip_max*kimg, nfftp
       call phase_error_with_msg(6, '<<m_CD_map_chgqenl_to_fft_box>>',__LINE__,__FILE__)
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

  subroutine m_CD_hardpart_sub2(nfout,is,ik,ib,chgq0)
    integer, intent(in) :: nfout, is,ik, ib
    real(kind=DP), intent(out) :: chgq0

    integer :: id_sname = -1
    call tstatc0_begin('m_CD_hardpart_sub2 ',id_sname,1)

    call summation_of_ff_sub_enl(nfout,is,ik,ib) ! -(m_C.D.) (vnlsum) fsr_l, fsi_l, occup_l --> hsr
    chgq_enl = 0.d0

    if(.not.noncol) then
       call add_hardpart_to_chgq_enl_IA(nfout,is,hsr,ndim_magmom)    ! (lclchg) add hardpart and make average
    else
       call add_hardpart_to_chgq_enl   (nfout,is,hsr,ndim_magmom)    ! (lclchg) add hardpart and make average
    end if
    chgq0 = chgq_enl(1,1)

    call tstatc0_end(id_sname)
  end subroutine m_CD_hardpart_sub2

  subroutine m_CD_hardpart_sub2_rs(nfout,is,ik,ib,bff,chgq0)
    integer, intent(in)      :: nfout, is, ik, ib
    real(kind=DP), intent(inout) :: bff( nfftp_nonpara)
    real(kind=DP), intent(out), optional :: chgq0
    real(kind=DP), allocatable, dimension(:) :: bfftmp
    integer :: ia,it,lmt1,lmt2,ilmta1,ilmta2,i,ind,lmtp
    integer :: nma
    real(kind=DP) :: fac
    real(kind=DP), allocatable, dimension(:) :: chgqtmp
    real(kind=DP) :: fac0,prod
    integer :: id_sname = -1
    call tstatc0_begin('m_CD_hardpart_sub2_rs ',id_sname,1)
    call summation_of_ff_sub_enl(nfout,is,ik,ib)
    bff=0.d0
    allocate(chgqtmp(nmesh_rs_aug_max));chgqtmp=0.d0
    do ia=1,natm
       it = ityp(ia)
       if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
       nma = nmesh_rs_aug(ia)
       chgqtmp = 0.d0
       do lmtp=1,nlmtpair(ia)
          lmt1 = plmt1(lmtp,ia)
          lmt2 = plmt2(lmtp,ia)
          fac0=2.0d0
          if(lmt1.eq.lmt2) fac0=1.d0
          fac = hsr(ia,lmt1,lmt2,is)*fac0
          do i=1,nma
             chgqtmp(i) = chgqtmp(i)+fac*qr_clm_ylm(i,ia,lmtp)
          enddo
       enddo
       do i=1,nma
          ind = kimg*(meshxyz_rs_aug(i,ia)-1)+1
          bff(ind) = bff(ind) + chgqtmp(i)/univol
          bff(ind+1) = 0.d0
       enddo
    enddo
    if(present(chgq0))then
       allocate(bfftmp(nfftp_nonpara));bfftmp(1:nfftp_nonpara) = bff(1:nfftp_nonpara)
       call m_FFT_CD0(nfout,bfftmp,DIRECT)
       chgq0 = bfftmp(kimg*(igfp_enl(1)-1)+1)/product(fft_box_size_CD(1:3,1))
       deallocate(bfftmp)
    endif
    deallocate(chgqtmp)
    call tstatc0_end(id_sname)
  end subroutine m_CD_hardpart_sub2_rs

! ============================= added by K. Tagami ===================== 11.0
  subroutine m_CD_hardpart_sub2_noncl(nfout,ik,ib,chgq0,chgq_enl_kt)
    integer, intent(in) :: nfout, ik, ib
    real(kind=DP), intent(out) :: chgq0( ndim_magmom )
    real(kind=DP), intent(out) :: chgq_enl_kt( kgp,kimg,ndim_magmom )
    real(kind=DP) ::  wct_start

    integer :: istmp
    integer :: id_sname = -1

    call tstatc0_begin('m_CD_hardpart_sub2_noncl ',id_sname,1)

    call summation_of_ff_sub_enl_noncl(nfout,ik,ib)
                  ! -(m_C.D.) (vnlsum) fsr_l, fsi_l, occup_l --> hsr

    chgq_enl_kt = 0.d0
    Do istmp=1, ndim_magmom
       call add_hardpart_to_chgq_enl( nfout, istmp, hsr, ndim_magmom, chgq_enl_kt )
    End do
    chgq0(1:ndim_magmom) = chgq_enl_kt(1,1,1:ndim_magmom)

    call tstatc0_end(id_sname)

  end subroutine m_CD_hardpart_sub2_noncl
! ================================================================== 11.0

  subroutine summation_of_ff_sub_enl(nfout,is,ik,ib)
    integer, intent(in) :: nfout,is, ik, ib

    real(kind=DP)   :: w_n
    integer         :: ia, lmt1, lmt2, p, q, it
!!$    real(kind=DP),pointer,dimension(:,:,:) :: hsr_mpi
    real(kind=DP),pointer, dimension(:,:) :: fs
!!$    integer ::         lmt_max
!!$    real(kind=DP),pointer, dimension(:,:,:) :: fx

!!$    lmt_max = maxval(ilmt(1:ntyp))
!!$    write(nfout,'(" lmt_max = ",i5)') lmt_max

!!$    allocate(fx(natm,lmt_max,2)); fx = 0.d0
    integer :: id_sname = -1
    call tstatc0_begin('summation_of_ff_sub_enl ',id_sname,1)

    if(k_symmetry(ik) == GAMMA) then
       allocate(fs(nlmta,1));    fs = 0.d0

       fs(1:nlmta,1) = fsr_l(map_z(ib),1:nlmta,ik)

       w_n = 2.d0
       hsr = 0.d0
       do ia = 1, natm
          it = ityp(ia)
          do lmt1 = 1, ilmt(it)
             p = lmta(lmt1,ia)
             do lmt2 = lmt1, ilmt(it)
                q = lmta(lmt2,ia)
                hsr(ia,lmt1,lmt2,is) = hsr(ia,lmt1,lmt2,is) &
                     &  + w_n * (  fs(p,1)*fs(q,1))
             end do! lmt2
          end do! lmt1
       end do! ia
    else
       allocate(fs(nlmta,2));    fs = 0.d0

       fs(1:nlmta,1) = fsr_l(map_z(ib),1:nlmta,ik)
       fs(1:nlmta,2) = fsi_l(map_z(ib),1:nlmta,ik)

       w_n = 2.d0
       hsr = 0.d0
       do ia = 1, natm
          it = ityp(ia)
          do lmt1 = 1, ilmt(it)
             p = lmta(lmt1,ia)
             do lmt2 = lmt1, ilmt(it)
                q = lmta(lmt2,ia)
                hsr(ia,lmt1,lmt2,is) = hsr(ia,lmt1,lmt2,is) &
                     &  + w_n * (  fs(p,1)*fs(q,1)+ fs(p,2)*fs(q,2))
             end do! lmt2
          end do! lmt1
       end do! ia
    end if

    deallocate(fs)
    call tstatc0_end(id_sname)
  end subroutine summation_of_ff_sub_enl

! ========================== added by K. Tagami ===================== 11.0
  subroutine summation_of_ff_sub_enl_noncl(nfout,ik,ib)
    integer, intent(in) :: nfout, ik, ib

    real(kind=DP)   :: w_n
    integer         :: ia, lmt1, lmt2, p, q, it
    integer :: is1, is2, istmp

    real(kind=DP), allocatable, dimension(:,:,:) :: fs
    real(kind=DP), allocatable, dimension(:,:,:,:) :: hsr_magmom
    real(kind=DP), allocatable, dimension(:,:,:,:) :: hsi_magmom

    allocate(fs(nlmta,2,ndim_spinor));   fs = 0.d0

    Do is1=1, ndim_spinor
       fs(1:nlmta,1,is1) = fsr_l(map_z(ib),1:nlmta,ik+is1-1)
       fs(1:nlmta,2,is1) = fsi_l(map_z(ib),1:nlmta,ik+is1-1)
    End do

    w_n = 1.d0
    hsr = 0.d0; hsi = 0.0d0

    do ia = 1, natm
       it = ityp(ia)
       Do is1=1, ndim_spinor
          Do is2=1, ndim_spinor
             istmp = ( is1 -1 )*ndim_spinor + is2

             do lmt1 = 1, ilmt(it)
                p = lmta(lmt1,ia)
                do lmt2 = lmt1, ilmt(it)
                   q = lmta(lmt2,ia)

                   hsr(ia,lmt1,lmt2,istmp) = hsr(ia,lmt1,lmt2,istmp) &
                        &  + w_n * ( fs(p,1,is1)*fs(q,1,is2) &
                        &          + fs(p,2,is1)*fs(q,2,is2) )
                   hsi(ia,lmt1,lmt2,istmp) = hsi(ia,lmt1,lmt2,istmp) &
                        &  + w_n * (-fs(p,1,is1)*fs(q,2,is2) &
                        &          + fs(p,2,is1)*fs(q,1,is2) )
                end do! lmt2
             end do! lmt1
          End do
       End do
    end do! ia

    allocate( hsr_magmom( natm, nlmt, nlmt, ndim_magmom ) ); hsr_magmom = 0.0d0
    allocate( hsi_magmom( natm, nlmt, nlmt, ndim_magmom ) ); hsi_magmom = 0.0d0

    call m_ES_DensMat_To_MagMom_hsr( natm, nlmt, hsr, hsi, hsr_magmom, hsi_magmom )
    hsr = hsr_magmom; hsi = hsi_magmom

    deallocate( hsr_magmom, hsi_magmom )
    deallocate(fs)

  end subroutine summation_of_ff_sub_enl_noncl
! ======================================================================= 11.0

! ===================================== modified by K. Tagami ================ 11.0
!  subroutine add_hardpart_to_chgq_enl(nfout,is,hsr)
!    integer, intent(in)      :: nfout, is
!    real(kind=DP), intent(in), dimension(natm,nlmt,nlmt,nspin):: hsr
  subroutine add_hardpart_to_chgq_enl( nfout, is, hsr, kspin, chgq_enl_kt )
    integer, intent(in)      :: nfout, is, kspin
    real(kind=DP), intent(in) :: hsr(natm,nlmt,nlmt,kspin)
    real(kind=DP), intent(inout), optional :: chgq_enl_kt( kgp,kimg,ndim_magmom )
! ========================================================================== 11.0

    real(kind=DP), pointer, dimension(:)               :: ylm

    integer :: it,lmt1,lmt2,n,ia,mdvdb,il1,tau1,il2,tau2,ilm3,l3,iiqitg
    real(kind=DP) :: fac !, tpos(3)

    integer ::                            kgp_adj, n_ialist, n_ialist0, ia_start, ia_end, n_iagroup, n_ia, ia_g
    integer, allocatable, dimension(:) :: ia_list
    real(kind=DP), allocatable, target, dimension(:,:) :: zfcos_x, zfsin_x
    real(kind=DP), allocatable, dimension(:) :: fac_x

    integer :: id_sname = -1
    call tstatc0_begin('add_hardpart_to_chgq_enl ',id_sname,1)

    if(modnrm == EXECUT) then

       call m_PP_find_maximum_l(n)   !  n-1: maximum l
       n = (n-1) + (n-1) + 1
       allocate(il3(n**2)); il3 = 0; call substitute_il3(n**2,il3) ! -(b_Elec..)! === Modified by K. Tagami === , "il3 = 0" is added
!       n_ialist = 1
       n_ialist = 4
       if(ipri>=1) write(nfout,'(" n_ialist = ",i8)') n_ialist
#ifdef HIUX
       n_ialist = 4
#endif
#ifdef VPP
       n_ialist = 8
#endif
#ifdef SX
       n_ialist = 8
#endif
       kgp_adj = kgp
       if(mod(kgp_adj,2) == 0) kgp_adj = kgp_adj + 1
       allocate(zfcos_x(kgp_adj,n_ialist)); zfcos_x = 0.d0
       allocate(zfsin_x(kgp_adj,n_ialist)); zfsin_x = 0.d0
       allocate(ia_list(n_ialist)); ia_list = 0
       allocate(fac_x(n_ialist))
!!$       allocate(zfcos(kgp))
!!$       allocate(zfsin(kgp))
!!$       zfcos = 0.d0; zfsin = 0.d0

! ====================================== Added by K. Tagami ==========
	fac_x = 0
! =====================================================================

       do it = 1, ntyp
          mdvdb = m_PP_include_vanderbilt_pot(it)
          if(mdvdb == SKIP) cycle

          n_ia = 0
          do ia = 1, natm
             if(ityp(ia) == it) n_ia = n_ia + 1
          end do

          if(n_ialist <= 0) call phase_error_with_msg(nfout, ' n_ialist is illegal',__LINE__,__FILE__)
          n_iagroup = n_ia/n_ialist + 1
          ia_start = 1
          do ia_g = 1, n_iagroup
             n_ialist0 = 0
             ia_list = 0
             AtomcountLoop: do ia = ia_start, natm
                if(ityp(ia) == it) then
                   n_ialist0 = n_ialist0 + 1
                   ia_list(n_ialist0) = ia
                end if
                if(n_ialist0 >= n_ialist) exit AtomcountLoop
             end do AtomcountLoop
             ia_start = ia+1
             if(n_ialist0 >= 1) then
                if(ipri >= 2) write(nfout,'(" !m_CD ia_list = ",8i8)') (ia_list(ia),ia=1,n_ialist0)

!!$          do ia = 1, natm
!!$             it = ityp(ia)
                call calc_phase_b(natm,pos,ia_list,n_ialist0,kgp,ngabc,1,kgp,1,kgp_adj,zfcos_x,zfsin_x)
!!$                   do ia = 1, n_ialist0
!!$                      call calc_phase2(natm,pos,ia_list(ia),kgp,ngabc,1,kgp,zfcos_x(1,ia),zfsin_x(1,ia))
!!$                call calc_phase2(natm,pos,ia,kgp,ngabc,1,kgp,zfcos,zfsin)
!!$                   end do
                ! -(b_Elec.)  -> zfcos, zfsin
                do lmt1 = 1,ilmt(it)
                   il1 = ltp(lmt1,it); tau1 = taup(lmt1,it)
                   do lmt2 = lmt1, ilmt(it)
                      il2 = ltp(lmt2,it); tau2 = taup(lmt2,it)
!!$                   fac = 2*iwei(ia); if(lmt1 == lmt2) fac = iwei(ia)
                      fac = 2.d0; if(lmt1 == lmt2 ) fac = 1.d0
                      do ia = 1, n_ialist0
                         fac_x(ia) = fac*iwei(ia_list(ia))
                      end do
                      do n = 1, il2p(lmt1,lmt2,it)
                         ilm3 = isph(lmt1,lmt2,n,it);    l3   =  il3(ilm3)
                         iiqitg = iqitg(il1,tau1,il2,tau2,l3+1,it)
                         if(iiqitg == 0) cycle
                         ylm => ylm_enl(:,ilm3)
                         call add_hardpart_to_chgq_e_core
                         !    ~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      end do! n
                   end do! lmt2
                end do! lmt1
             end if
          end do! ia_g
       end do! it

!!$       deallocate(zfsin); deallocate(zfcos); deallocate(il3)
       deallocate(fac_x,ia_list,zfsin_x,zfcos_x,il3)
    end if

    if(iprichargedensity >= 2) then
       write(nfout,*) ' -- before average --'
! ========================================= modified by K. Tagami ============== 11.0
!       call m_CD_wd_chgq_l_small_portion(nfout)
       if ( noncol ) then
         call m_CD_wd_chgq_l_portion_noncl(nfout)
       else
         call m_CD_wd_chgq_l_small_portion(nfout)
       endif
! ============================================================================= 11.0
    endif

! ========================== added by K. Tagami ============= 11.0
    if ( noncol ) then
       chgq_enl_kt(:,:,is) = chgq_enl(:,:)
    endif
! ============================================================ 11.0

! ========================== modified by K. Tagami ============== 11.0
!    if(nbztyp >= SIMPLE_CUBIC)  call charge_average_enl()
!
    if (nbztyp >= SIMPLE_CUBIC) then
       if ( noncol ) then
          if ( is == ndim_magmom ) then
             call charge_average_enl_noncl2( chgq_enl_kt )
          endif
       else
          call charge_average_enl()
       endif
    endif
! ================================================================ 11.0

    call tstatc0_end(id_sname)
  contains
    subroutine add_hardpart_to_chgq_e_core
      integer       :: i
      real(kind=DP) :: dga, flchgq, f,f2, flchgq_1, flchgq_2, flchgq_3, flchgq_4 &
           &                            , flchgq_5, flchgq_6, flchgq_7, flchgq_8, zdga
!!$      real(kind=DP), pointer, dimension(:) :: zfsc1,zfsc2
      real(kind=DP), pointer, dimension(:) :: zfsc1_1,zfsc2_1, zfsc1_2, zfsc2_2 &
           &         , zfsc1_3,zfsc2_3,zfsc1_4,zfsc2_4, zfsc1_5,zfsc2_5,zfsc1_6,zfsc2_6 &
           &         , zfsc1_7,zfsc2_7,zfsc1_8,zfsc2_8
      dga = dl2p(lmt1,lmt2,n,it)
      if(mod(l3,2) == 0) then
         if(n_ialist0 >= 1) zdga = real(zi**(-l3))*dga
         if(n_ialist0 >= 1) flchgq_1 = fac_x(1)*zdga*hsr(ia_list(1),lmt1,lmt2,is)
         if(n_ialist0 >= 2) flchgq_2 = fac_x(2)*zdga*hsr(ia_list(2),lmt1,lmt2,is)
         if(n_ialist0 >= 3) flchgq_3 = fac_x(3)*zdga*hsr(ia_list(3),lmt1,lmt2,is)
         if(n_ialist0 >= 4) flchgq_4 = fac_x(4)*zdga*hsr(ia_list(4),lmt1,lmt2,is)
         if(n_ialist0 >= 5) flchgq_5 = fac_x(5)*zdga*hsr(ia_list(5),lmt1,lmt2,is)
         if(n_ialist0 >= 6) flchgq_6 = fac_x(6)*zdga*hsr(ia_list(6),lmt1,lmt2,is)
         if(n_ialist0 >= 7) flchgq_7 = fac_x(7)*zdga*hsr(ia_list(7),lmt1,lmt2,is)
         if(n_ialist0 >= 8) flchgq_8 = fac_x(8)*zdga*hsr(ia_list(8),lmt1,lmt2,is)
         if(kimg == 1) then
            if(n_ialist0 >= 1) zfsc1_1 => zfcos_x(:,1)
            if(n_ialist0 >= 2) zfsc1_2 => zfcos_x(:,2)
            if(n_ialist0 >= 3) zfsc1_3 => zfcos_x(:,3)
            if(n_ialist0 >= 4) zfsc1_4 => zfcos_x(:,4)
            if(n_ialist0 >= 5) zfsc1_5 => zfcos_x(:,5)
            if(n_ialist0 >= 6) zfsc1_6 => zfcos_x(:,6)
            if(n_ialist0 >= 7) zfsc1_7 => zfcos_x(:,7)
            if(n_ialist0 >= 8) zfsc1_8 => zfcos_x(:,8)
         else
            f2 = -1
            if(n_ialist0 >= 1) then
               zfsc1_1 => zfcos_x(:,1); zfsc2_1 => zfsin_x(:,1)
            end if
            if(n_ialist0 >= 2) then
               zfsc1_2 => zfcos_x(:,2); zfsc2_2 => zfsin_x(:,2)
            end if
            if(n_ialist0 >= 3) then
               zfsc1_3 => zfcos_x(:,3); zfsc2_3 => zfsin_x(:,3)
            end if
            if(n_ialist0 >= 4) then
               zfsc1_4 => zfcos_x(:,4); zfsc2_4 => zfsin_x(:,4)
            end if
            if(n_ialist0 >= 5) then
               zfsc1_5 => zfcos_x(:,5); zfsc2_5 => zfsin_x(:,5)
            end if
            if(n_ialist0 >= 6) then
               zfsc1_6 => zfcos_x(:,6); zfsc2_6 => zfsin_x(:,6)
            end if
            if(n_ialist0 >= 7) then
               zfsc1_7 => zfcos_x(:,7); zfsc2_7 => zfsin_x(:,7)
            end if
            if(n_ialist0 >= 8) then
               zfsc1_8 => zfcos_x(:,8); zfsc2_8 => zfsin_x(:,8)
            end if
         end if
      else
!!$         flchgq = fac*aimag(zi**(-l3))*dga*hsr(ia,lmt1,lmt2,is)
         if(n_ialist0 >= 1) zdga = aimag(zi**(-l3))*dga
         if(n_ialist0 >= 1) flchgq_1 = fac_x(1)*zdga*hsr(ia_list(1),lmt1,lmt2,is)
         if(n_ialist0 >= 2) flchgq_2 = fac_x(2)*zdga*hsr(ia_list(2),lmt1,lmt2,is)
         if(n_ialist0 >= 3) flchgq_3 = fac_x(3)*zdga*hsr(ia_list(3),lmt1,lmt2,is)
         if(n_ialist0 >= 4) flchgq_4 = fac_x(4)*zdga*hsr(ia_list(4),lmt1,lmt2,is)
         if(n_ialist0 >= 5) flchgq_5 = fac_x(5)*zdga*hsr(ia_list(5),lmt1,lmt2,is)
         if(n_ialist0 >= 6) flchgq_6 = fac_x(6)*zdga*hsr(ia_list(6),lmt1,lmt2,is)
         if(n_ialist0 >= 7) flchgq_7 = fac_x(7)*zdga*hsr(ia_list(7),lmt1,lmt2,is)
         if(n_ialist0 >= 8) flchgq_8 = fac_x(8)*zdga*hsr(ia_list(8),lmt1,lmt2,is)
         if(kimg == 1) then
!!$            zfsc1 => zfsin
            if(n_ialist0 >= 1) zfsc1_1 => zfsin_x(:,1)
            if(n_ialist0 >= 2) zfsc1_2 => zfsin_x(:,2)
            if(n_ialist0 >= 3) zfsc1_3 => zfsin_x(:,3)
            if(n_ialist0 >= 4) zfsc1_4 => zfsin_x(:,4)
            if(n_ialist0 >= 5) zfsc1_5 => zfsin_x(:,5)
            if(n_ialist0 >= 6) zfsc1_6 => zfsin_x(:,6)
            if(n_ialist0 >= 7) zfsc1_7 => zfsin_x(:,7)
            if(n_ialist0 >= 8) zfsc1_8 => zfsin_x(:,8)
         else
!!$            f2 = 1; zfsc1 => zfsin; zfsc2 => zfcos
            f2 = 1
            if(n_ialist0 >= 1) then
               zfsc1_1 => zfsin_x(:,1); zfsc2_1 => zfcos_x(:,1)
            end if
            if(n_ialist0 >= 2) then
               zfsc1_2 => zfsin_x(:,2); zfsc2_2 => zfcos_x(:,2)
            end if
            if(n_ialist0 >= 3) then
               zfsc1_3 => zfsin_x(:,3); zfsc2_3 => zfcos_x(:,3)
            end if
            if(n_ialist0 >= 4) then
               zfsc1_4 => zfsin_x(:,4); zfsc2_4 => zfcos_x(:,4)
            end if
            if(n_ialist0 >= 5) then
               zfsc1_5 => zfsin_x(:,5); zfsc2_5 => zfcos_x(:,5)
            end if
            if(n_ialist0 >= 6) then
               zfsc1_6 => zfsin_x(:,6); zfsc2_6 => zfcos_x(:,6)
            end if
            if(n_ialist0 >= 7) then
               zfsc1_7 => zfsin_x(:,7); zfsc2_7 => zfcos_x(:,7)
            end if
            if(n_ialist0 >= 8) then
               zfsc1_8 => zfsin_x(:,8); zfsc2_8 => zfcos_x(:,8)
            end if
         end if
      end if


      if(n_ialist0 == 1) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
            do i = 1, kgp     !for mpi
               chgq_enl(i,1) = chgq_enl(i,1) &
                    & +flchgq_1*qitg_enl(i,iiqitg)*ylm(i)*zfsc1_1(i)
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
            do i = 1, kgp   !for mpi
               f = flchgq_1*qitg_enl(i,iiqitg)*ylm(i)
               chgq_enl(i,1) = chgq_enl(i,1) + f * zfsc1_1(i)
               chgq_enl(i,2) = chgq_enl(i,2) + f2 * f * zfsc2_1(i)
            end do
         end if
      else if(n_ialist0 == 2) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
            do i = 1, kgp     !for mpi
               chgq_enl(i,1) = chgq_enl(i,1) &
                    & + qitg_enl(i,iiqitg)*ylm(i)*( &
                    &        flchgq_1*zfsc1_1(i) &
                    &       +flchgq_2*zfsc1_2(i) )
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
            do i = 1, kgp   !for mpi
               chgq_enl(i,1) = chgq_enl(i,1) &
                    &  + qitg_enl(i,iiqitg)*ylm(i)*( &
                    &         flchgq_1 * zfsc1_1(i) &
                    &        +flchgq_2 * zfsc1_2(i) )
               chgq_enl(i,2) = chgq_enl(i,2) &
                    &  + f2 * qitg_enl(i,iiqitg)*ylm(i)*( &
                    &         flchgq_1 * zfsc2_1(i) &
                    &        +flchgq_2 * zfsc2_2(i))
            end do
         end if
      else if(n_ialist0 == 3) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
            do i = 1, kgp     !for mpi
               chgq_enl(i,1) = chgq_enl(i,1) &
                    & + qitg_enl(i,iiqitg)*ylm(i)*( &
                    &        flchgq_1*zfsc1_1(i) + flchgq_2*zfsc1_2(i) &
                    &       +flchgq_3*zfsc1_3(i) )
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
            do i = 1, kgp   !for mpi
               chgq_enl(i,1) = chgq_enl(i,1) &
                    &  + qitg_enl(i,iiqitg)*ylm(i)*( &
                    &         flchgq_1*zfsc1_1(i) + flchgq_2*zfsc1_2(i) &
                    &        +flchgq_3*zfsc1_3(i) )
               chgq_enl(i,2) = chgq_enl(i,2) &
                    &  + f2 * qitg_enl(i,iiqitg)*ylm(i)*( &
                    &         flchgq_1*zfsc2_1(i) + flchgq_2*zfsc2_2(i) &
                    &        +flchgq_3*zfsc2_3(i) )
            end do
         end if
      else if(n_ialist0 == 4) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
            do i = 1, kgp     !for mpi
               chgq_enl(i,1) = chgq_enl(i,1) &
                    & + qitg_enl(i,iiqitg)*ylm(i)*( &
                    &        flchgq_1*zfsc1_1(i) + flchgq_2*zfsc1_2(i) &
                    &       +flchgq_3*zfsc1_3(i) + flchgq_4*zfsc1_4(i))
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
            do i = 1, kgp   !for mpi
               chgq_enl(i,1) = chgq_enl(i,1) &
                    &  + qitg_enl(i,iiqitg)*ylm(i)*( &
                    &         flchgq_1*zfsc1_1(i) + flchgq_2*zfsc1_2(i) &
                    &        +flchgq_3*zfsc1_3(i) + flchgq_4*zfsc1_4(i))
               chgq_enl(i,2) = chgq_enl(i,2) &
                    &  + f2 * qitg_enl(i,iiqitg)*ylm(i)*( &
                    &         flchgq_1*zfsc2_1(i) + flchgq_2*zfsc2_2(i) &
                    &        +flchgq_3*zfsc2_3(i) + flchgq_4*zfsc2_4(i) )
            end do
         end if
      else if(n_ialist0 == 5) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
            do i = 1, kgp     !for mpi
               chgq_enl(i,1) = chgq_enl(i,1) &
                    & + qitg_enl(i,iiqitg)*ylm(i)*( &
                    &        flchgq_1*zfsc1_1(i) + flchgq_2*zfsc1_2(i) &
                    &       +flchgq_3*zfsc1_3(i) + flchgq_4*zfsc1_4(i) + flchgq_5*zfsc1_5(i))
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
            do i = 1, kgp   !for mpi
               chgq_enl(i,1) = chgq_enl(i,1) &
                    &  + qitg_enl(i,iiqitg)*ylm(i)*( &
                    &         flchgq_1*zfsc1_1(i) + flchgq_2*zfsc1_2(i) &
                    &        +flchgq_3*zfsc1_3(i) + flchgq_4*zfsc1_4(i) &
                    &        +flchgq_5*zfsc1_5(i))
               chgq_enl(i,2) = chgq_enl(i,2) &
                    &  + f2 * qitg_enl(i,iiqitg)*ylm(i)*( &
                    &         flchgq_1*zfsc2_1(i) + flchgq_2*zfsc2_2(i) &
                    &        +flchgq_3*zfsc2_3(i) + flchgq_4*zfsc2_4(i) &
                    &        +flchgq_5*zfsc2_5(i))
            end do
         end if
      else if(n_ialist0 == 6) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
            do i = 1, kgp     !for mpi
               chgq_enl(i,1) = chgq_enl(i,1) &
                    & + qitg_enl(i,iiqitg)*ylm(i)*( &
                    &        flchgq_1*zfsc1_1(i) + flchgq_2*zfsc1_2(i) &
                    &       +flchgq_3*zfsc1_3(i) + flchgq_4*zfsc1_4(i) &
                    &       +flchgq_5*zfsc1_5(i) + flchgq_6*zfsc1_6(i))
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
            do i = 1, kgp   !for mpi
               chgq_enl(i,1) = chgq_enl(i,1) &
                    &  + qitg_enl(i,iiqitg)*ylm(i)*( &
                    &         flchgq_1*zfsc1_1(i) + flchgq_2*zfsc1_2(i) &
                    &        +flchgq_3*zfsc1_3(i) + flchgq_4*zfsc1_4(i) &
                    &        +flchgq_5*zfsc1_5(i) + flchgq_6*zfsc1_6(i))
               chgq_enl(i,2) = chgq_enl(i,2) &
                    &  + f2 * qitg_enl(i,iiqitg)*ylm(i)*( &
                    &         flchgq_1*zfsc2_1(i) + flchgq_2*zfsc2_2(i) &
                    &        +flchgq_3*zfsc2_3(i) + flchgq_4*zfsc2_4(i) &
                    &        +flchgq_5*zfsc2_5(i) + flchgq_6*zfsc2_6(i) )
            end do
         end if
      else if(n_ialist0 == 7) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
            do i = 1, kgp     !for mpi
               chgq_enl(i,1) = chgq_enl(i,1) &
                    & + qitg_enl(i,iiqitg)*ylm(i)*( &
                    &        flchgq_1*zfsc1_1(i) + flchgq_2*zfsc1_2(i) &
                    &       +flchgq_3*zfsc1_3(i) + flchgq_4*zfsc1_4(i) &
                    &       +flchgq_5*zfsc1_5(i) + flchgq_6*zfsc1_6(i) &
                    &       +flchgq_7*zfsc1_7(i) )
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
            do i = 1, kgp   !for mpi
               chgq_enl(i,1) = chgq_enl(i,1) &
                    &  + qitg_enl(i,iiqitg)*ylm(i)*( &
                    &         flchgq_1*zfsc1_1(i) + flchgq_2*zfsc1_2(i) &
                    &        +flchgq_3*zfsc1_3(i) + flchgq_4*zfsc1_4(i) &
                    &        +flchgq_5*zfsc1_5(i) + flchgq_6*zfsc1_6(i) &
                    &        +flchgq_7*zfsc1_7(i))
               chgq_enl(i,2) = chgq_enl(i,2) &
                    &  + f2 * qitg_enl(i,iiqitg)*ylm(i)*( &
                    &         flchgq_1*zfsc2_1(i) + flchgq_2*zfsc2_2(i) &
                    &        +flchgq_3*zfsc2_3(i) + flchgq_4*zfsc2_4(i) &
                    &        +flchgq_5*zfsc2_5(i) + flchgq_6*zfsc2_6(i) &
                    &        +flchgq_7*zfsc2_7(i))
            end do
         end if
      else if(n_ialist0 == 8) then
         if(kimg == 1) then
#ifdef HIUX
*POPTION PARALLEL
#endif
            do i = 1, kgp     !for mpi
               chgq_enl(i,1) = chgq_enl(i,1) &
                    & + qitg_enl(i,iiqitg)*ylm(i)*( &
                    &        flchgq_1*zfsc1_1(i) + flchgq_2*zfsc1_2(i) &
                    &       +flchgq_3*zfsc1_3(i) + flchgq_4*zfsc1_4(i) &
                    &       +flchgq_5*zfsc1_5(i) + flchgq_6*zfsc1_6(i) &
                    &       +flchgq_7*zfsc1_7(i) + flchgq_8*zfsc1_8(i) )
            end do
         else
#ifdef HIUX
*POPTION PARALLEL
#endif
            do i = 1, kgp   !for mpi
               chgq_enl(i,1) = chgq_enl(i,1) &
                    &  + qitg_enl(i,iiqitg)*ylm(i)*( &
                    &         flchgq_1*zfsc1_1(i) + flchgq_2*zfsc1_2(i) &
                    &        +flchgq_3*zfsc1_3(i) + flchgq_4*zfsc1_4(i) &
                    &        +flchgq_5*zfsc1_5(i) + flchgq_6*zfsc1_6(i) &
                    &        +flchgq_7*zfsc1_7(i) + flchgq_8*zfsc1_8(i))
               chgq_enl(i,2) = chgq_enl(i,2) &
                    &  + f2 * qitg_enl(i,iiqitg)*ylm(i)*( &
                    &         flchgq_1*zfsc2_1(i) + flchgq_2*zfsc2_2(i) &
                    &        +flchgq_3*zfsc2_3(i) + flchgq_4*zfsc2_4(i) &
                    &        +flchgq_5*zfsc2_5(i) + flchgq_6*zfsc2_6(i) &
                    &        +flchgq_7*zfsc2_7(i) + flchgq_8*zfsc2_8(i))
            end do
         end if
      end if

    end subroutine add_hardpart_to_chgq_e_core
  end subroutine add_hardpart_to_chgq_enl

  subroutine add_hardpart_to_chgq_enl_IA( nfout, ispin, hsr, kspin)
!      Translated from "add_hardpart_to_chgq_l_in_keworld" in the 3d-version m_Charge_Density.F90, by T. Yamasaki,  19 Jul. 2021
    integer, intent(in)      :: nfout, ispin, kspin
    real(kind=DP), intent(in), dimension(natm,nlmt,nlmt,kspin) :: hsr

    integer :: is, it,lmt1,lmt2,n,ia,mdvdb,il1,il2,ilm3,l3
    real(kind=DP) :: fac

    real(kind=DP), allocatable, dimension(:) :: ylm_sum
    real(kind=DP)                              :: shdg_s, zdga
    integer :: m, maxm, ip, np, iq
    integer :: mc ! maxval(nc)
    integer, parameter :: mcritical = 4*2+1
    integer, allocatable, dimension(:) :: nqitg_sp, nqitg_sp0 !d(ntyp)
    integer, allocatable, dimension(:) :: iq2l3 ! d(nqitg)
    integer, allocatable, dimension(:,:) :: nc  ! d(maxm,nqitg)
    integer, allocatable, dimension(:,:,:) :: nc2lmt1, nc2lmt2, nc2n ! d(mc,maxm,nqitg)
    integer :: ibl1,ibl2,ibsize,ncache,iwidth
    real(kind=DP) :: ph, f
    integer :: iy, i
    integer :: ista, iend
    integer,save :: printcount = 0
#ifdef DEBUG_LDOS
    integer,parameter :: critical_printoutlevel = 1
#else
    integer,parameter :: critical_printoutlevel = 2
#endif

    integer :: id_sname = -1
    call tstatc0_begin('add_hardpart_to_chgq_IA ',id_sname,1)

    if(modnrm == EXECUT) then
       ista = 1
       iend = kgp
       is = ispin

       call m_PP_find_maximum_l(n)   !  n-1: maximum l
       n = (n-1) + (n-1) + 1
       allocate(il3(n**2)); il3 = 0; call substitute_il3(n**2,il3) ! -(b_Elec..)

       allocate(nqitg_sp(ntyp)); allocate(nqitg_sp0(ntyp)); nqitg_sp = 0; nqitg_sp0 = 0
       allocate(iq2l3(nqitg))  ; iq2l3 = 0
       allocate(nc(mcritical,nqitg));nc=0

       call m_PP_set_index_arrays1(nfout,ntyp,nqitg,mcritical,n**2,il3,maxm,mc,nqitg_sp,nqitg_sp0,iq2l3,nc)
       allocate(nc2lmt1(mc,maxm,nqitg)); nc2lmt1 = 0
       allocate(nc2lmt2(mc,maxm,nqitg)); nc2lmt2 = 0
       allocate(nc2n(mc,maxm,nqitg));    nc2n = 0

       call m_PP_set_index_arrays2(nfout,mc,maxm,nqitg,mcritical,n**2,il3,iq2l3 &
            & ,nc2lmt1,nc2lmt2,nc2n,nc) ! -> nc2lmt1, nc2lmt2, nc2n, nc

!!$       ncache = 0
!!$       ncache = (1024*1024)*3/4
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
       if(iprichargedensity >= critical_printoutlevel-1 .and. printcount == 0) then
           write(nfout,'(" !mCD:(add_hardpart_to_chgq_enl_with_index_arrays) ibsize, iwidth = ",2i8 &
                &                  ," ncache = ",i8)') ibsize, iwidth, ncache
           printcount = 1
        end if

       allocate(zfcos(ibsize),zfsin(ibsize))
       allocate(ylm_sum(ibsize))

       do ibl1=ista,iend,ibsize         ! G-space
          ibl2=min(ibl1+ibsize-1,kgp)

          do ia = 1, natm
             it = ityp(ia)
             mdvdb = m_PP_include_vanderbilt_pot(it)
             if(mdvdb == SKIP) cycle

             do i = 1, ibl2-ibl1+1
                iy = ibl1-1+i
                ph = (ngabc(iy,1)*pos(ia,1)+ngabc(iy,2)*pos(ia,2)+ngabc(iy,3)*pos(ia,3))*PAI2
                zfcos(i) = dcos(ph)
                zfsin(i) = dsin(ph)
             end do

             do iq = nqitg_sp0(it), nqitg_sp(it)
                l3 = iq2l3(iq)
                ylm_sum = 0.d0
                do m = 1, 2*l3+1
                   ilm3 = l3*l3+m
!fj                   call sum_hsr_dot_gauntc0(it,ia,iq,m,iqm) ! hsr, dl2p -> shdg_s
                   shdg_s = 0.d0
                   do ip = 1, nc(m,iq)
                      lmt1 = nc2lmt1(ip,m,iq)
                      lmt2 = nc2lmt2(ip,m,iq)
                      np = nc2n(ip,m,iq)
                      fac = 2.d0; if(lmt1 == lmt2) fac = 1.d0
                      shdg_s = shdg_s + fac*iwei(ia)*hsr(ia,lmt1,lmt2,is)*dl2p(lmt1,lmt2,np,it)
                   end do
                   do i = 1, ibl2-ibl1+1
                      ylm_sum(i) = ylm_sum(i) + shdg_s*ylm_enl(ibl1-1+i,ilm3)  ! ylm_enl, shdg_s -> ylm_sum
                   end do
                end do

                if(mod(l3,2) == 0) then
                   zdga = real(zi**(-l3))
!fj                   call add_hardpart_to_chgq_l_div0(zdga,iq)
                           !! iq, ylm_sum, exp(-iGR), qitg_enl -> chgq_enl
                   if(kimg == 1) then
                      do i = 1, ibl2-ibl1+1  ! G-space
                         chgq_enl(ibl1-1+i,1) = chgq_enl(ibl1-1+i,1) + zdga*qitg_enl(ibl1-1+i,iq)*ylm_sum(i)*zfcos(i)
                      end do
                   else
                      do i = 1, ibl2-ibl1+1  ! G-space
                         f = zdga*qitg_enl(ibl1-1+i,iq)*ylm_sum(i)
                         chgq_enl(ibl1-1+i,1) = chgq_enl(ibl1-1+i,1) + f * zfcos(i)
                         chgq_enl(ibl1-1+i,2) = chgq_enl(ibl1-1+i,2) - f * zfsin(i)
                      end do
                   end if
                else
                   zdga = aimag(zi**(-l3))
!fj                   call add_hardpart_to_chgq_l_div1(zdga,iq)
                           !! iq, ylm_sum, exp(-iGR), qitg_enl -> chgq_enl
                   if(kimg == 1) then
                      do i = 1, ibl2-ibl1+1  ! G-space
                         chgq_enl(ibl1-1+i,1) = chgq_enl(ibl1-1+i,1) + zdga*qitg_enl(ibl1-1+i,iq)*ylm_sum(i)*zfsin(i)
                      end do
                   else
                      do i = 1, ibl2-ibl1+1  ! G-space
                         f = zdga*qitg_enl(ibl1-1+i,iq)*ylm_sum(i)
                         chgq_enl(ibl1-1+i,1) = chgq_enl(ibl1-1+i,1) +  f * zfsin(i)
                         chgq_enl(ibl1-1+i,2) = chgq_enl(ibl1-1+i,2) +  f * zfcos(i)
                      end do
                   end if
                end if
             end do
          end do
       end do


       deallocate(ylm_sum)
       deallocate(zfsin,zfcos)

       deallocate(nc2n,nc2lmt2,nc2lmt1,nc,iq2l3,nqitg_sp,nqitg_sp0)
    end if

    if(nbztyp >= SIMPLE_CUBIC .or. af /=0 ) then
       call charge_average_enl()
    end if
    call tstatc0_end(id_sname)
  end subroutine add_hardpart_to_chgq_enl_IA

  subroutine charge_average_enl()
    integer ::       ng, no, ngp, no1, no2
    real(kind=DP) :: fi, tx,ty,tz, fp, fc, fs, zcr, zci
    real(kind=DP), pointer, dimension(:,:) :: work2

    allocate(work2(1:kgp,kimg)); work2 = 0.d0

    fi = 1.d0/nopr
    no1 = 1; no2 = nopr

    work2 = 0.d0                   ! initialization
    do no = no1, no2

       tx = tau(1,no,BUCS)*PAI2
       ty = tau(2,no,BUCS)*PAI2
       tz = tau(3,no,BUCS)*PAI2

       if(kimg == 1) then
          do ng = 1, kgp
             ngp = ngpt_enl(ng,no)
             fp = ngabc(ngp,1)*tx + ngabc(ngp,2)*ty + ngabc(ngp,3)*tz
             work2(ng,1)        = work2(ng,1) + dcos(fp)*chgq_enl(ngp,1)
          end do
       else if(kimg == 2) then
          do ng = 1, kgp
             ngp= ngpt_enl(ng,no)
             fp = ngabc(ngp,1)*tx + ngabc(ngp,2)*ty + ngabc(ngp,3)*tz
             fc = dcos(fp);     fs = dsin(fp)
             zcr= chgq_enl(ngp,1);  zci= chgq_enl(ngp,kimg)
             work2(ng,1)        = work2(ng,1) + fc*zcr - fs*zci
             work2(ng,2)        = work2(ng,2) + fc*zci + fs*zcr
          end do
       end if
    end do

    chgq_enl(:,:) = work2(:,:)*fi

    deallocate(work2)
  end subroutine charge_average_enl

! =============================== added by K. Tagami ================= 11.0
  subroutine charge_average_enl_noncl2( chgq_enl_kt )
    real(kind=DP), intent(inout) :: chgq_enl_kt( kgp,kimg,ndim_magmom )

    integer ::       ng, no, ngp, no1, no2
    real(kind=DP) :: tx,ty,tz, fp, fc, fs, factor
    integer ::       is, ii

    real(kind=DP), allocatable :: work2(:,:)
    real(kind=DP), allocatable :: mag_work2(:,:,:)
!
    real(kind=DP) :: mag_tmp(kimg,3)

    allocate(work2(1:kgp,kimg)); work2 = 0.d0
    allocate(mag_work2(1:kgp,kimg,3)); mag_work2 = 0.d0

! ----------------------------
    factor = 1.d0 /dble(nopr);     no1 = 1;       no2 = nopr
! --------------------------

    do ng = 1, kgp
       do no = no1, no2

          tx = tau(1,no,BUCS)*PAI2
          ty = tau(2,no,BUCS)*PAI2
          tz = tau(3,no,BUCS)*PAI2

          ngp= ngpt_enl(ng,no)
          fp = ngabc(ngp,1)*tx + ngabc(ngp,2)*ty + ngabc(ngp,3)*tz
          fc = dcos(fp);     fs = dsin(fp)

          Do ii=1, 3
             mag_tmp(1:kimg,ii) = op(ii,1,no) *chgq_enl_kt(ngp,1:kimg,2) &
                  &             + op(ii,2,no) *chgq_enl_kt(ngp,1:kimg,3) &
                  &             + op(ii,3,no) *chgq_enl_kt(ngp,1:kimg,4)
          End Do

          work2(ng,1) = work2(ng,1) +fc *chgq_enl_kt(ngp,1,1) -fs*chgq_enl_kt(ngp,2,1)
          work2(ng,2) = work2(ng,2) +fc *chgq_enl_kt(ngp,2,1) +fs*chgq_enl_kt(ngp,1,1)

          mag_work2(ng,1,1:3) = mag_work2(ng,1,1:3) &
               &               +fc *mag_tmp(1,1:3) -fs *mag_tmp(2,1:3)
          mag_work2(ng,2,1:3) = mag_work2(ng,2,1:3) &
               &               +fc *mag_tmp(2,1:3) +fs *mag_tmp(1,1:3)
       end do
    end do

    chgq_enl_kt(:,:,1) = work2(:,:) *factor
    Do is=2, ndim_magmom
       chgq_enl_kt(:,:,is) = mag_work2(:,:,is-1) *factor
    End do

    deallocate(work2); deallocate( mag_work2 )
  end subroutine charge_average_enl_noncl2
! =============================================================== 11.0

! ============== KT_add ======= 2014/08/25
  subroutine m_CD_hardpart_hsr_add( nfout, kv3 )
    integer, intent(in) :: nfout, kv3

    call summation_of_ff_add(kv3)
    call symmtrz_of_ff_add

  contains

    subroutine summation_of_ff_add(kv3)
      integer, intent(in) :: kv3

      real(kind=DP)   :: w_n, d_factor
      integer         :: ia, is, k, i, lmt1, lmt2, p, q, it
      real(kind=DP),pointer,dimension(:,:,:,:) :: hsr_mpi

      d_factor = 2.d0/kv3
      hsr_add = 0.d0

      do ia = 1, natm
         it = ityp(ia)
         do is = 1, nspin, af+1
            do i = 1, np_e                                  ! MPI
               do k = is, kv3+is-nspin, nspin
                  if(map_k(k) /= myrank_k) cycle            ! MPI
                  w_n = occup_l(i,k)*d_factor

                  do lmt1 = 1, ilmt_add(it)
                     p = lmta_add(lmt1,ia)
                     if(k_symmetry(k) == GAMMA) then
                        do lmt2 = lmt1, ilmt_add(it)
                           q = lmta_add(lmt2,ia)
                           hsr_add(ia,lmt1,lmt2,is)  &
                                &  = hsr_add(ia,lmt1,lmt2,is) &
                                &   + w_n *fsr_add_l(i,p,k) *fsr_add_l(i,q,k)
                        end do! lmt2
                     else
                        do lmt2 = lmt1, ilmt_add(it)
                           q = lmta_add(lmt2,ia)
                           hsr_add(ia,lmt1,lmt2,is) &
                                &  = hsr_add(ia,lmt1,lmt2,is) &
                                &  + w_n * ( fsr_add_l(i,p,k)*fsr_add_l(i,q,k) &
                                &          + fsi_add_l(i,p,k)*fsi_add_l(i,q,k) )
                        end do! lmt2
                     end if
                  end do! lmt1
               end do! ik
            end do! i
         end do! is
      end do! ia

      Do ia = 1, natm
         it = ityp(ia)
         Do is = 1, nspin, af+1
            Do lmt1=1, ilmt_add(it)
               Do lmt2=lmt1, ilmt_add(it)
                  hsr_add(ia,lmt2,lmt1,is) = hsr_add(ia,lmt1,lmt2,is)
               End do
            End do
         End do
      End do

      if(npes >= 2) then
         allocate(hsr_mpi(natm,nlmt_add,nlmt_add,nspin)); hsr_mpi = 0.0d0
!         call mpi_allreduce( hsr_add, hsr_mpi,natm*nlmt_add*nlmt_add*nspin, &
!              &              mpi_double_precision, mpi_sum, MPI_CommGroup,ierr )
         call mpi_allreduce( hsr_add, hsr_mpi,natm*nlmt_add*nlmt_add*nspin, &
              &              mpi_double_precision, mpi_sum, mpi_spin_group,ierr )
         hsr_add = hsr_mpi
         deallocate(hsr_mpi)
      end if
    end subroutine summation_of_ff_add

    subroutine symmtrz_of_ff_add
      integer         :: ia, is, iopr, i, lmt1, lmt2, lmt3, lmt4, it
      real(kind=DP),pointer,dimension(:,:,:,:) :: hsr_mpi
      real(kind=DP),pointer,dimension(:,:,:,:) :: hsr_tmp
      integer :: il1,im1,it1,il2,im2,it2,il3,im3,it3,il4,im4,it4
      integer :: ii,jj,kk,ll,n,m,iii,jjj
      integer :: ja

      allocate(hsr_mpi(natm,nlmt_add,nlmt_add,nspin))
      allocate(hsr_tmp(natm,nlmt_add,nlmt_add,nspin))

      hsr_tmp = hsr_add

      do ia=1,natm
         it=ityp(ia)
         do is =1,nspin,af+1
            do lmt2=1,ilmt_add(it)
               do lmt1=lmt2+1,ilmt_add(it)
                  hsr_tmp(ia,lmt1,lmt2,is)=hsr_tmp(ia,lmt2,lmt1,is)
               end do
            end do
         end do
      end do

      hsr_add = 0.d0

      do iopr=1,nopr
         do ia = 1, natm
            it = ityp(ia)
            ja=abs(ia2ia_symmtry_op_inv(ia,iopr))

            do is = 1, nspin, af+1
               do lmt1 = 1, ilmt_add(it)
                  il1=ltp_add(lmt1,it); im1=mtp_add(lmt1,it); it1=1
                  ii=(il1-1)**2+im1

                  do lmt2 = lmt1, ilmt_add(it)
                     il2=ltp_add(lmt2,it); im2=mtp_add(lmt2,it); it2=1
                     jj=(il2-1)**2+im2

                     do n=1,nylm_paw(ii,iopr,ia)
                        iii=iylm_paw(n,ii,iopr,ia)

                        do m=1,nylm_paw(jj,iopr,ia)
                           jjj=iylm_paw(m,jj,iopr,ia)

                           do lmt3=1,ilmt_add(it)
                              il3=ltp_add(lmt3,it); im3=mtp_add(lmt3,it); it3=1
                              kk=(il3-1)**2+im3
                              if(kk.ne.iii .or. it1.ne.it3) cycle

                              do lmt4=1,ilmt_add(it)
                                 il4=ltp_add(lmt4,it); im4=mtp_add(lmt4,it); it4=1
                                 ll=(il4-1)**2+im4
                                 if(ll.ne.jjj .or. it2.ne.it4) cycle

                                 hsr_add(ia,lmt1,lmt2,is) = hsr_add(ia,lmt1,lmt2,is) &
                                      &                   + hsr_tmp(ja,lmt3,lmt4,is) &
                                      &                     *crotylm_paw(n,ii,iopr,ia) &
                                      &                     *crotylm_paw(m,jj,iopr,ia)
                              end do! lmt4
                           end do! lmt3

                        end do! jjj
                     end do! iii

                  end do! lmt2
               end do! lmt1
            end do! is
         end do! ia
      end do! iopr

      hsr_add = hsr_add /nopr

      if(af /= 0 .and. flg_paw) then
         do ia = 1, natm
            ja=abs(ia2ia_symmtry_op_inv(ia,nopr+af))
            if(ja <= 0 .or. natm < ia) cycle
            hsr_add(ia,:,:,nspin) = hsr_add(ja,:,:,1)
         end do
      end if

      deallocate(hsr_mpi,hsr_tmp) ! MPI

    end subroutine symmtrz_of_ff_add

  end subroutine m_CD_hardpart_hsr_add

  subroutine m_CD_hardpart_hsr_add_noncl( nfout, kv3 )
    integer, intent(in) :: nfout, kv3

    call summation_of_ff_add_noncl(kv3)
    call symmtrz_of_ff_add_noncl

  contains

    subroutine summation_of_ff_add_noncl(kv3)
      integer, intent(in) :: kv3

      real(kind=DP)   :: w_n, d_factor
      integer         :: ia, is, k, i, lmt1, lmt2, p, q, it

      integer :: is1, is2, is_tmp, k1, k2

      real(kind=DP), allocatable, dimension(:,:,:,:) :: hsr_or_hsi_mpi
      real(kind=DP), allocatable, dimension(:,:,:,:) :: hsr_ssrep
      real(kind=DP), allocatable, dimension(:,:,:,:) :: hsi_ssrep
      real(kind=DP), allocatable, dimension(:,:,:,:) :: hsr_with_soc
      real(kind=DP), allocatable, dimension(:,:,:,:) :: hsi_with_soc
!
      d_factor = 1.d0/ ( kv3 /ndim_spinor )

      allocate( hsr_ssrep( natm, nlmt_add, nlmt_add, ndim_chgpot ) ); hsr_ssrep = 0.0d0
      allocate( hsi_ssrep( natm, nlmt_add, nlmt_add, ndim_chgpot ) ); hsi_ssrep = 0.0d0

      do ia = 1, natm
         it = ityp(ia)

         do is1 = 1, ndim_spinor
            do is2 = 1, ndim_spinor
               is_tmp = ( is1 -1 )*ndim_spinor + is2
               do i = 1, np_e                                  ! MPI

                  do k = 1, kv3, ndim_spinor
                     if ( map_k(k) /= myrank_k ) cycle            ! MPI
                     w_n = occup_l(i,k) *d_factor

                     k1 = k + is1 -1;  k2 = k + is2 -1
                     do lmt1 = 1, ilmt_add(it)
                        p = lmta_add(lmt1,ia)

                        do lmt2 = 1, ilmt_add(it)
                           q = lmta_add(lmt2,ia)
                           hsr_ssrep(ia,lmt1,lmt2,is_tmp) &
                                &  = hsr_ssrep(ia,lmt1,lmt2,is_tmp) &
                                &    + w_n * ( fsr_add_l(i,p,k1)*fsr_add_l(i,q,k2) &
                                &            + fsi_add_l(i,p,k1)*fsi_add_l(i,q,k2) )
                           hsi_ssrep(ia,lmt1,lmt2,is_tmp) &
                                &  = hsi_ssrep(ia,lmt1,lmt2,is_tmp) &
                                &    + w_n * ( -fsr_add_l(i,p,k1)*fsi_add_l(i,q,k2) &
                                &              +fsi_add_l(i,p,k1)*fsr_add_l(i,q,k2) )
                        end do! lmt2

                     end do! lmt1
                  end do! k
               end do! i
            end do! is2
         end do! is1
      end do! ia
      !
      if (npes >= 2) then
         allocate(hsr_or_hsi_mpi(natm,nlmt_add,nlmt_add,ndim_chgpot))
         hsr_or_hsi_mpi = 0.0d0

!         call mpi_allreduce( hsr_ssrep, hsr_or_hsi_mpi, &
!              &              natm*nlmt_add*nlmt_add*ndim_chgpot, &
!              &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
         call mpi_allreduce( hsr_ssrep, hsr_or_hsi_mpi, &
              &              natm*nlmt_add*nlmt_add*ndim_chgpot, &
              &              mpi_double_precision, mpi_sum, mpi_spin_group, ierr )
         hsr_ssrep = hsr_or_hsi_mpi

!         call mpi_allreduce( hsi_ssrep, hsr_or_hsi_mpi, &
!              &              natm*nlmt_add*nlmt_add*ndim_chgpot, &
!              &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
         call mpi_allreduce( hsi_ssrep, hsr_or_hsi_mpi, &
              &              natm*nlmt_add*nlmt_add*ndim_chgpot, &
              &              mpi_double_precision, mpi_sum, mpi_spin_group, ierr )
         hsi_ssrep = hsr_or_hsi_mpi

         deallocate(hsr_or_hsi_mpi)
      end if
! ---------------

      call m_ES_DensMat_To_MagMom_hsr( natm, nlmt_add, hsr_ssrep, hsi_ssrep, &
           &                           hsr_add, hsi_add )
!-
      deallocate( hsr_ssrep, hsi_ssrep )
    end subroutine summation_of_ff_add_noncl

    subroutine symmtrz_of_ff_add_noncl
      integer         :: ia, is, iopr, i, lmt1, lmt2, lmt3, lmt4, it
      real(kind=DP), allocatable, dimension(:,:,:,:) :: hsr_tmp
      real(kind=DP), allocatable, dimension(:,:,:,:) :: hsi_tmp

      integer :: il1,im1,it1,il2,im2,it2,il3,im3,it3,il4,im4,it4
      integer :: ii,jj,kk,ll,n,m,iii,jjj
      integer :: ja

      integer :: ixyz1, ixyz2, is_tmp
      real(kind=DP) :: ctmp1, weight

      allocate(hsr_tmp(natm,nlmt_add,nlmt_add,ndim_magmom)); hsr_tmp = 0.0d0
      allocate(hsi_tmp(natm,nlmt_add,nlmt_add,ndim_magmom)); hsi_tmp = 0.0d0
!
      hsr_tmp = hsr_add;  hsi_tmp = hsi_add

      do ia=1,natm
         it=ityp(ia)
         do is =1, ndim_magmom
            do lmt2=1,ilmt_add(it)
               do lmt1=lmt2+1,ilmt_add(it)
                  hsr_tmp(ia,lmt1,lmt2,is) =  hsr_tmp(ia,lmt2,lmt1,is)
                  hsi_tmp(ia,lmt1,lmt2,is) = -hsi_tmp(ia,lmt2,lmt1,is)
               end do
            end do
         end do
      end do

      hsr_add = 0.d0; hsi_add = 0.0d0

      do iopr=1,nopr
         if ( magmom_dir_inversion_opr_flag(iopr) == -1 ) then
            weight = -1.0d0
         else
            weight = 1.0d0
         endif

         do ia = 1, natm
            it = ityp(ia)
            ja=abs(ia2ia_symmtry_op_inv(ia,iopr))

            do lmt1 = 1, ilmt_add(it)
               il1=ltp_add(lmt1,it); im1=mtp_add(lmt1,it);  it1=1
               ii=(il1-1)**2+im1

               do lmt2 = lmt1, ilmt_add(it)
                  !             do lmt2 = 1, ilmt_add(it)
                  il2=ltp_add(lmt2,it); im2=mtp_add(lmt2,it); it2=1
                  jj=(il2-1)**2+im2

                  do n=1,nylm_paw(ii,iopr,ia)
                     iii=iylm_paw(n,ii,iopr,ia)
                     do m=1,nylm_paw(jj,iopr,ia)
                        jjj=iylm_paw(m,jj,iopr,ia)

                        do lmt3=1,ilmt_add(it)
                           il3=ltp_add(lmt3,it); im3=mtp_add(lmt3,it); it3=1
                           kk=(il3-1)**2+im3

                           if(kk.ne.iii .or. it1.ne.it3) cycle

                           do lmt4=1,ilmt_add(it)
                              il4=ltp_add(lmt4,it); im4=mtp_add(lmt4,it); it4=1
                              ll=(il4-1)**2+im4
                              if(ll.ne.jjj .or. it2.ne.it4) cycle

                              hsr_add(ia,lmt1,lmt2,1) = &
                                   hsr_add(ia,lmt1,lmt2,1) + &
                                   hsr_tmp(ja,lmt3,lmt4,1)* &
                                   crotylm_paw(n,ii,iopr,ia)* &
                                   crotylm_paw(m,jj,iopr,ia)
                              hsi_add(ia,lmt1,lmt2,1) = &
                                   hsi_add(ia,lmt1,lmt2,1) + &
                                   weight * &
                                   hsi_tmp(ja,lmt3,lmt4,1)* &
                                   crotylm_paw(n,ii,iopr,ia)* &
                                   crotylm_paw(m,jj,iopr,ia)

                              Do ixyz1=1, 3
                                 Do ixyz2=1, 3
                                    ctmp1 = op(ixyz2, ixyz1, iopr) *weight

                                    hsr_add(ia,lmt1,lmt2,ixyz2+1) &
                                         & = hsr_add(ia,lmt1,lmt2,ixyz2+1)  &
                                         &  + ctmp1 &
                                         &    *hsr_tmp(ja,lmt3,lmt4,ixyz1+1) &
                                         &    *crotylm_paw(n,ii,iopr,ia)  &
                                         &    *crotylm_paw(m,jj,iopr,ia)
                                    hsi_add(ia,lmt1,lmt2,ixyz2+1) &
                                         & = hsi_add(ia,lmt1,lmt2,ixyz2+1)  &
                                         &  + op(ixyz2, ixyz1, iopr) &
                                         &    *hsi_tmp(ja,lmt3,lmt4,ixyz1+1) &
                                         &    *crotylm_paw(n,ii,iopr,ia)  &
                                         &    *crotylm_paw(m,jj,iopr,ia)
                                 End do
                              End do

                           end do! lmt4
                        end do! lmt3

                     end do! jjj
                  end do! iii

               end do! lmt2
            end do! lmt1
         end do! ia
      end do! iopr

      hsr_add = hsr_add/nopr;  hsi_add = hsi_add /nopr

      do ia=1,natm
         it=ityp(ia)
         do is =1, ndim_magmom
            do lmt2=1,ilmt_add(it)
               do lmt1=lmt2+1,ilmt_add(it)
                  hsr_add(ia,lmt1,lmt2,is) =  hsr_add(ia,lmt2,lmt1,is)
                  hsi_add(ia,lmt1,lmt2,is) = -hsi_add(ia,lmt2,lmt1,is)
               end do
            end do
         end do
      end do

      deallocate(hsr_tmp);  deallocate(hsi_tmp)

    end subroutine symmtrz_of_ff_add_noncl

  end subroutine m_CD_hardpart_hsr_add_noncl
! ======================== 2014/08/25
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
! ======================================================================= 11.0
    endif
  end subroutine m_CD_rspace_put_headermark

  subroutine m_CD_rspace_put_endmark(nfchr)
    integer, intent(in) :: nfchr
    if(mype==0) then
       write(nfchr,'(a3)') 'END'
    endif
  end subroutine m_CD_rspace_put_endmark

  subroutine m_CD_softpart_rwf2(nfout,nfrwf2,kv3,ik,ib,center,nr,radius)
    integer, intent(in) :: nfout, nfrwf2, kv3, ik, ib
    real(kind=DP), intent(in) :: center(3)
    integer, intent(in) :: nr
    real(kind=DP), intent(in) :: radius(1:nr)

    real(kind=DP) :: wf2q_l(ista_kngp:iend_kngp,kimg)
    real(kind=DP) :: zfcos(ista_kngp:iend_kngp)
    real(kind=DP) :: zfsin(ista_kngp:iend_kngp)
    real(kind=DP) :: rwf2(nr)
    integer ib1
    integer :: id_sname = -1
    call tstatc0_begin('m_CD_softpart_wfn ',id_sname,1)

    call calc_phase2(1,center(1),1,kgp,ngabc,ista_kngp,iend_kngp,zfcos,zfsin)

! ================================== Modified by K. Tagami =============
!    allocate(afft(nfft)); allocate(bfft(nfft)); call m_FFT_alloc_WF_work()
    allocate(afft(nfft)); allocate(bfft(nfft));
    afft = 0; bfft = 0 ; call m_FFT_alloc_WF_work()
! =========================================================================
    wf2q_l = 0.d0
    afft = 0.d0
    ib1 = neordr(ib,ik)
    if(map_ek(ib1,ik) == mype) then
       call m_ES_WF_in_Rspace(ik,ib1,bfft) ! (swffft)
       if(ipri >= 2) write(6,'(" !! ik = ,",i8," ib1 = ",i8)') ik,ib
       call add_density()  ! -(this module) bfft -> afft
    end if
    if(npes >= 2) then
!       call mpi_allreduce(afft,bfft,nfft,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr) ! MPI
       call mpi_allreduce(afft,bfft,nfft,mpi_double_precision,mpi_sum,mpi_spin_group,ierr) ! MPI
       afft = bfft                          ! MPI
    end if
    call m_FFT_WF(ELECTRON,nfout,afft,DIRECT,OFF)
    call substitute_CD_for_wf2q()
    call radial_wf2(nr,radius,rwf2)
    if(mype == 0) call wd_radial_wf2(nr,radius,rwf2,nfrwf2)
    deallocate(afft); deallocate(bfft); call m_FFT_dealloc_WF_work()
    call tstatc0_end(id_sname)
  contains
    subroutine substitute_CD_for_wf2q
      integer       :: i, ri, i1
      real(kind=DP) :: fac
      integer       :: iend !mpi
      real(kind=DP) :: chgr, chgi

      fac = 4.d0*PAI/(univol*product(fft_box_size_WF(1:3,1)))
      iend = iend_kngp
      if( iend_kngp > kg ) iend = kg
      if( ista_kngp <= iend ) then
         if(kimg==1) then
            do i = ista_kngp, iend  !for mpi
               i1 = igf(i)
               wf2q_l(i,kimg) = afft(i1)*fac*zfcos(i)
            end do
         else
           do i = ista_kngp, iend  !for mpi
              i1 = 2*igf(i)
               chgr = afft(i1-1)*fac
               chgi = afft(i1)*fac
               wf2q_l(i,1)    = chgr*zfcos(i) - chgi*zfsin(i)
               wf2q_l(i,kimg) =  chgr*zfsin(i) + chgi*zfcos(i)
            end do
         end if
      endif
    end subroutine substitute_CD_for_wf2q

    subroutine add_density
      integer  :: i
      do i = 1, nfft-1, 2
         afft(i) = afft(i) + (bfft(i)**2+bfft(i+1)**2) ! MPI
      end do
    end subroutine add_density

    subroutine radial_wf2(nr,radius,rwf2)
      integer, intent(in) :: nr
      real(kind=DP), intent(in) :: radius(nr)
      real(kind=DP), intent(out) :: rwf2(nr)

      real(kind=DP) :: jgr(ista_kngp:iend_kngp)
      real(kind=DP) :: gr(ista_kngp:iend_kngp)
      real(kind=DP) :: rwf2_mpi(nr)
      integer :: ri, ir, i, iend, nsize

      rwf2 = 0.d0
      do ir=1,nr
         do i = ista_kngp, iend  !for mpi
            gr(i) = gr_l(i)*radius(ir)
         end do
         nsize = iend_kngp-ista_kngp+1
         call dsjnv(0,nsize,gr,jgr)
         iend = iend_kngp
         if( iend_kngp > kg ) iend = kg
         if( ista_kngp <= iend ) then
            do i = ista_kngp, iend  !for mpi
               rwf2(ir) = rwf2(ir) + wf2q_l(i,1)*jgr(i)
            end do
         endif
      end do
      if(npes >= 2) then
         call mpi_allreduce(rwf2,rwf2_mpi,nr,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr) ! MPI
         rwf2 = rwf2_mpi                          ! MPI
      end if
    end subroutine radial_wf2

    subroutine wd_radial_wf2(nr,radius,rwf2,nfrwf2)
      integer, intent(in) :: nr
      real(kind=DP), intent(in) :: radius(1:nr),rwf2(1:nr)
      integer, intent(in) :: nfrwf2

      integer :: ir

      do ir=1,nr
         write(nfrwf2,'(i7,2(1x,e20.12))') ir, radius(ir), rwf2(ir)
      end do
    end subroutine wd_radial_wf2

  end subroutine m_CD_softpart_rwf2

! ============================= added by K. Tagami ============== 11.0
  subroutine m_CD_softpart_rwf2_noncl(nfout,nfrwf2,kv3,ik,ib,center,nr,radius)
    integer, intent(in) :: nfout, nfrwf2, kv3, ik, ib
    real(kind=DP), intent(in) :: center(3)
    integer, intent(in) :: nr
    real(kind=DP), intent(in) :: radius(1:nr)

  end subroutine m_CD_softpart_rwf2_noncl
!=============================================================== 11.0

  subroutine m_CD_den_mat(nfout,kv3)
    integer, intent(in) :: nfout,kv3
    !! debug
                                                 __TIMER_SUB_START(736)

! ========================== modiifed by K. Tagami ============= 11.0
!    write(*,*) 'calling summation_of_ff'
!    call summation_of_ff(kv3)
!
    if ( noncol ) then
!       write(*,*) 'calling summation_of_ff_noncl'
       call summation_of_ff_noncl(kv3)
    else
!       write(*,*) 'calling summation_of_ff'
       call summation_of_ff(kv3)
    endif
! ============================================================= 11.0

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
    allocate(cmat(nlp2,nmp2,nnp2,ndim_magmom))
    if(interpolation_method_chg == SPLINE_INTERPOLATION) allocate(cmat2(nlp2,nmp2,nnp2,ndim_magmom))
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

    do is=1,ndim_magmom
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
    do is=1,ndim_magmom
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
       if ( .not. noncol ) then
          if (nspin>1) total_charge = total_charge + chgq_l(1,1,2)*univol
       endif
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
              &  ,mpi_sum,MPI_CommGroup,ierr)
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
                  grt = pos_wk(ia,1)*ngabc(i,1) + pos_wk(ia,2)*ngabc(i,2) &
                       &                        + pos_wk(ia,3)*ngabc(i,3)
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
                  grt = pos_wk(ia,1)*ngabc(i,1) + pos_wk(ia,2)*ngabc(i,2) &
                       &                        + pos_wk(ia,3)*ngabc(i,3)
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
                  grt = pos(ia,1)*ngabc(i,1) + pos(ia,2)*ngabc(i,2) &
                       &                     + pos(ia,3)*ngabc(i,3)
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
                  grt = pos(ia,1)*ngabc(i,1) + pos(ia,2)*ngabc(i,2) &
                       &                     + pos(ia,3)*ngabc(i,3)
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

! ==== KT_add ======= 2015/05/16
!
!  m_CD_softpart_ktsub_noncl and m_CD_hardpart_ktsub_noncl
!  are used for writing squared wavefuctions at ik_spec, ib_spec
!
! -------------
 subroutine m_CD_softpart_ktsub_noncl( nfout, kv3, ik_spec, ib_spec )
    integer, intent(in) :: ik_spec, ib_spec, nfout, kv3

    integer ispin, ib1, ik, i, ip, max_elements, icolumn, istart, iend, icycle, ic
    integer :: id_sname = -1

    real(kind=DP), allocatable, dimension(:) :: wf_phase
    real(kind=DP) :: occupation
!
    integer :: is, is1, is2, is_tmp
! ----------------
    real(kind=DP), allocatable :: afft_kt(:,:)
    real(kind=DP), allocatable :: bfft_kt(:,:)
    real(kind=DP), allocatable :: chgq_magmom( :,:,: )

! ----------------------- start -------------------

    allocate(afft_kt(nfft,ndim_chgpot)); afft_kt = 0.0d0;
    allocate(bfft_kt(nfft,ndim_spinor)); bfft_kt = 0.0d0;
    call m_FFT_alloc_WF_work()

    chgq_l = 0.d0

    Do ib1 = ista_e, iend_e, istep_e     ! MPI
       if ( ib1 /= ib_spec ) cycle

       Do ik = 1, kv3, ndim_spinor
          if ( map_k(ik) /= myrank_k ) cycle! MPI

          if ( ik /= ik_spec ) cycle

          occupation = 1.0d0

          Do is=1, ndim_spinor
            call m_ES_WF_in_Rspace( ik +is-1, ib1, bfft_kt(:,is) )
          End do

          Do is1=1, ndim_spinor
             Do is2=1, ndim_spinor
                is_tmp = ( is1 -1 )*ndim_spinor + is2
                call add_occupied_density_matrix()
             End do
          End do
       End do
    End do

    bfft_kt = 0.0d0
    Do is1= 1, ndim_spinor
       Do is2=1, ndim_spinor
          is_tmp = ( is1 -1 )*ndim_spinor + is2
          if ( npes >= 2 ) then
!             call mpi_allreduce( afft_kt(:,is_tmp), bfft_kt(:,1), nfft, &
!                  &              mpi_double_precision, mpi_sum, &
!                  &              MPI_CommGroup, ierr )
             call mpi_allreduce( afft_kt(:,is_tmp), bfft_kt(:,1), nfft, &
                  &              mpi_double_precision, mpi_sum, &
                  &              mpi_spin_group, ierr )
             afft_kt(:,is_tmp) = bfft_kt(:,1)
          endif
          call m_FFT_WF( ELECTRON,nfout,afft_kt(:,is_tmp),DIRECT,OFF )
          call substitute_CD_for_chgq()
       end do
    end do
! --
    allocate( chgq_magmom( ista_kngp:iend_kngp,kimg,ndim_magmom ) )
    chgq_magmom = 0.0d0
!
    call m_ES_DensMat_To_MagMom_Gspace( chgq_l, chgq_magmom )
    chgq_l = chgq_magmom

    deallocate( chgq_magmom )
    deallocate(afft_kt); deallocate(bfft_kt);
    call m_FFT_dealloc_WF_work()

  contains

    subroutine substitute_CD_for_chgq
      integer       :: i, ri, i1
      real(kind=DP) :: fac
      integer       :: iend !mpi

      fac = 1.d0 /( univol *product(fft_box_size_WF(1:3,1)) )

      do ri = 1, kimg
         iend = iend_kngp
         if( iend_kngp > kg ) iend = kg
         if( ista_kngp <= iend ) then
            do i = ista_kngp, iend  !for mpi
               i1 = kimg*igf(i) + (ri - kimg)
               chgq_l(i,ri,is_tmp) = afft_kt(i1,is_tmp)*fac
            end do
         endif
      end do
    end subroutine substitute_CD_for_chgq

    subroutine add_occupied_density_matrix
      integer  :: i
      real(kind=DP) :: cr, ci

      do i = 1, nfft-1, 2
         cr =  bfft_kt(i,  is1) *bfft_kt(i,  is2) &
         &    +bfft_kt(i+1,is1) *bfft_kt(i+1,is2)
!
! ---------------------------------------------
         ci = -bfft_kt(i,  is1) *bfft_kt(i+1,is2) &
         &    +bfft_kt(i+1,is1) *bfft_kt(i,  is2)
!
! ----
!         ci = bfft_kt(i,  is1) *bfft_kt(i+1,is2) &
!           &  -bfft_kt(i+1,is1) *bfft_kt(i,  is2)
!
! --
         afft_kt(i,  is_tmp) = afft_kt(i,  is_tmp) + occupation *cr
         afft_kt(i+1,is_tmp) = afft_kt(i+1,is_tmp) + occupation *ci
      end do
    end subroutine add_occupied_density_matrix

  end subroutine m_CD_softpart_ktsub_noncl

  subroutine m_CD_hardpart_ktsub_noncl( nfout, ik_spec, ib_spec )
    integer, intent(in) :: ik_spec, ib_spec, nfout

    call summation_of_ff_sub_noncl( ik_spec, ib_spec )
                               ! -(m_C.D.) (vnlsum) fsr_l, fsi_l, occup_l --> hsr
    call add_hardpart_to_chgq_l( nfout, ndim_magmom, hsr, NO )

  end subroutine m_CD_hardpart_ktsub_noncl
! ================ 2015/05/16


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
       call m_CD_softpart(nfout,kv3)
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
