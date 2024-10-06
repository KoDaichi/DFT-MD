!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 606 $)
!
!  MODULE:  m_PlaneWaveBasisSet
!
!  AUTHOR(S): T. Yamasaki, M. Okamoto,   August/20/2003
!
!  FURTHER MODIFICATION: T. Yamasaki, April/15/2006
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
#ifdef FJ_TIMER
#   define __TIMER_FJ_START_w_BARRIER(str,a)   call mpi_barrier(str,ierr) ;   call timer_sta(a)
#   define __TIMER_FJ_START(a)   call timer_sta(a)
#   define __TIMER_FJ_STOP(a)    call timer_end(a)
#else
#   define __TIMER_FJ_START_w_BARRIER(str,a)
#   define __TIMER_FJ_START(a)
#   define __TIMER_FJ_STOP(a)
#endif
#ifdef __TIMER_INIDO__
#   define __TIMER_INIDO_START(a)   call timer_sta(a)
#   define __TIMER_INIDO_STOP(a)    call timer_end(a)
#else
#   define __TIMER_INIDO_START(a)
#   define __TIMER_INIDO_STOP(a)
#endif

module m_PlaneWaveBasisSet
!    ( m_pwBS )
! $Id: m_PlaneWaveBasisSet.F90 606 2020-04-15 06:45:49Z ktagami $
  use m_Crystal_Structure,  only : nopr, altv,rltv, m_CS_op_in_PUCV,univol, univol_prev
  use m_Kpoints,            only : kv3,kv3_ek,vkxyz,qwgt, k_symmetry
  use m_FFT,                only : fft_box_size_WF, fft_box_size_CD, fft_box_size_pWF &
       &                        , fft_box_size_CD_c, fft_box_size_CD_nonpara, fft_box_size_WF_prev &
       &                        , fft_box_size_CD_prev, fft_box_size_CD_nonpara_prev

  use m_Files,              only : nfout
  use m_Timing,             only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters, only : ipri,iprifftmap,printable,gmax,gmaxp,gmaxp_reduced, gmax_positron &
       &                         , intzaj,cutoff_mix,af,neg,nel_Ylm, sw_positron &
       &                         , fixed_charge_k_parallel, sw_interpolate_wfs, neg_previous, nspin &
       &                         , sw_modified_kpoint_increment
  use m_Const_Parameters,   only : CARTS, BUCS, CRDTYP, OFF &
       &                         , SURFACE, WHOLE_BZ, SIMPLE_CUBIC &
       &                         , BCC,  FCC, DIAMOND, HEXAGONAL, ORTHORHOMBIC, RUTILE &
       &                         , C2v_SURFACE, Csy_SURFACE, GENERAL, GENERAL_LARGER, DP &
       &                         , by_matrix_diagon, by_random_numbers &
       &                         , PAI, PAI4,SMALL,MEDIUM,LARGE, DELTA, DELTA10, DELTA07 &
       &                         , SmallestPositiveNumber, GAMMA, GAMMA_base_symmetrization &
       &                         , ONE_BY_ONE, ON
  use m_Parallelization,    only : is_kngp,ie_kngp,ista_k,iend_k,ista_kngp,iend_kngp &
       &                         , npes,mype,ierr, MPI_CommGroup, myrank_g, nrank_g &
       &                         , m_Parallel_init_mpi_kngp_prev, ista_kngp_prev, iend_kngp_prev, nrank_k &
       &                         , nis_kv3_ek

! =================== KT_add ========================= 13.0F&13.0U2
  use m_Control_Parameters,  only : sw_hybrid_functional, gmax_exx, gmaxp_exx, use_fft_exx, &
       &                            sw_modified_TFW_functional
  use m_FFT, only : fft_box_size_exx, fft_box_size_CD_exx, fft_box_size_CD_nonpara
! ==================================================== 13.0F&13.0U2

  use m_IterationNumbers, only : nkgroup
  use mpi

  implicit none
!  include 'mpif.h'

  integer, parameter                    :: how_give_initial_zaj = 0
  integer, public,allocatable,target,dimension(:,:)  :: ngabc  !d(kgp,3)
!!$  integer, public,allocatable,target,dimension(:,:)  :: ngabc_kngp_l
!!$  integer, public,allocatable,target,dimension(:,:)  :: ngabc_kngp_B_l
  real(DP),allocatable, dimension(:,:)    :: fp_l
  integer, allocatable, dimension(:,:)  :: ngpt_l !d(ista_kngp:iend_kngp,nopr+af)
  real(DP),allocatable, dimension(:)    :: gr_l   !d(ista_kngp:iend_kngp)
  integer, allocatable, dimension(:)    :: igfp_l !d(ista_kngp:iend_kngp)
  integer, allocatable, dimension(:)    :: igfp_l_prev !d(ista_kngp:iend_kngp)
  integer, allocatable, dimension(:)    :: igfp_nonpara !d(ista_kngp:iend_kngp)
#ifdef _MPIFFT_
  integer, allocatable, dimension(:)    :: igfp_l_c !d(ista_kngp:iend_kngp)
#endif
  integer, allocatable, dimension(:)    :: igf    !d(kg)
  integer, allocatable, dimension(:)    :: igf_prev
  integer, target                       :: kg, kgp, kgp_reduced, kg1_ext=0, kg1=0, kgpm &
	&                                , kg0, kg_gamma = 0, kg2_gamma = 0, kgp_reduced_prev
  integer                               :: ngshell ! depending on mype
  integer, allocatable, dimension(:,:)  :: ngshell_range ! d(ngshell,2)
  integer                               :: ngshell_GAMMA ! not depend on mype
  ! ------- Positron start
  integer ::                               kg_pwf, kg1_pwf=0
  integer, dimension(3) ::                 n_rGv_pstrn
  integer, allocatable, dimension(:) ::    igf_pstrn !d(kg_pwf)
  ! ------- Positron end

! ======= KT_add ==================================== 13.0F
  integer ::                               kg_exx, kgp_exx, kg1_exx=0,kg1p_exx=0
  integer, dimension(3) ::                 n_rGv_exx,n_rGpv_exx
  integer, allocatable, dimension(:) ::    igf_exx !d(kg_exx)
  integer, allocatable, dimension(:) ::    igfp_exx !d(kgp_exx)
! =================================================== 13.0F

! ====== KT_add ============ 13.0U2
  integer :: kg_tfw = 0
! ========================== 13.0U2

! ===== EXP_CELLOPT ==== 2015/09/24
  integer :: kg1_prev = 0
  integer :: kgp_prev = 0
! ====================== 2015/09/24

! === Band_Unfolding
  integer, allocatable :: GVec_on_refcell(:,:)
! ====

  integer, dimension(3)                 :: n_rGv
  integer, dimension(3)                 :: n_rGpv
  integer, dimension(3)                 :: n_rGpv_reduced
  integer, dimension(3)                 :: n_rGv_prev
  integer, dimension(3)                 :: n_rGpv_prev
  ! n_rG(p)v is the size of a rhombohedron that contains
  ! a sphere whose radius is 2*Gmax (Gmaxp) .

  integer, allocatable, target, dimension(:,:)  :: nbase !d(kg1_ext,kv3) PW Basis Set for each k-points
  integer, allocatable, target, dimension(:,:)  :: nbase_prev !d(kg1_ext,kv3) PW Basis Set for each k-points
  integer, allocatable, dimension(:,:)  :: nbmat ! sub PWBS for each k-points
  !                                                  d(nmatsz,ista_k:iend_k)
  integer, allocatable, dimension(:,:)  :: nbmat2! sub PWBS for each k-points
  !                                                  d(nmatsz,ista_k:iend_k)
  integer, allocatable, target, dimension(:) :: iba    !d(kv3)
  integer, allocatable, target, dimension(:) :: iba_prev    !d(kv3)
  integer, allocatable, dimension(:)    :: iba_ek !d(kv3_ek)
  integer, allocatable, dimension(:)    :: iba2   !d(kv3*how_give_initial_zaj)
  integer, allocatable, target, dimension(:,:)  :: nbase_gamma !d(kg_gamma,2)
  integer, allocatable, target, dimension(:,:)  :: nbase_gamma_prev !d(kg_gamma,2)
  integer                               :: nbmx   ! maximum # of elements in PWBSs
  real(kind=DP),private                 :: gmaxs
  integer                               :: nmatsz ! subspace Hermitian matrix size
  integer                               :: nmatsz2
  integer, dimension(3)                 :: n_rGsv
  ! n_rGhv(1:3) is the size of a rhombohedron that contains a sphare whose
  ! redius is 2*Gmaxs
  integer, allocatable, dimension(:,:,:):: igpo
  !             d(-n_rGsv(1):n_rGsv(1),-n_rGsv(2):n_rGsv(2),-n_rGsv(3):n_rGsv(3))

  real(kind=DP), dimension(6) :: ttr
  real(kind=DP), allocatable,         dimension(:,:,:) :: op_br
  real(kind=DP), allocatable, target, dimension(:,:)   :: ylm_l !d(ista_kngp:iend_kngp, nel_Ylm)

!  1-1.  m_pwBS_alloc_ylm_l            <- Initilal_Electronic_Structure
!  1-2.  m_pwBS_dealloc_ylm_l          <- Initilal_Electronic_Structure
!  1-3.  m_pwBS_alloc_igpo             <- m_ES_WF_by_MatDiagon
!  1-4.  m_pwBS_dealloc_igpo           <- m_ES_WF_by_MatDiagon
!  1-5.  alloc_nbase                   <- (3)
!  1-6.  reduce_nbase                  <- (3)
!  1-7.  alloc_iba                     <- (3)
!  1-8.  G_vectors_alloc               <- (7)
!  1-9.  m_pwBS_alloc_nbmat_and_iba2   <- m_ES_WF_by_MatDiagon
!  1-10. m_pwBS_dealloc_nbmat_and_iba2 <- m_ES_WF_by_MatDiagon
!  1-11. m_pwBS_alloc_ngpt_igfp_gr     <- Preparation
!  1-12. m_pwBS_dealloc_ngpt_igfp_gr   <- Finalization
!  2.    m_pwBS_mat_for_each_WF        <- m_ES_WF_by_MatDiagon
!  3.    m_pwBS_for_each_WF            <- Preparation
!        - wd_ngshell_range
!  4.    m_pwBS_positronWF             <- Preparation
!  5.    m_pwBS_set_gmaxs              <- m_ES_WF_by_MatDiagon
!        - check_gmax
!  6.    m_pwBS_assume_G_rhombohedron  <- Preparation
!        - size_of_Gvector_rhombohedron, - check_the_size_of_rhombohedron
!        - wd_n_rGv_and_n_rGpv
!  7.    m_pwBS_generate_G_vectors     <- Preparation
!        - get_Gvectors_in_a_gmaxp_sphere, - count_Gvectors_in_spheres
!        - adjust_n_rGv_to_2l3m5n
!  8.    m_pwBS_calc_length_of_G       <- Preparation
!  9.    m_pwBS_setup_FFTmapfunctions  <- Preparation
!  10.   m_pwBS_GminusGmapfunction     <- m_ES_WF_by_MatDiagon
!  11.   decide_g_list_size            <- (12)
!  12.   m_pwBS_G_trans_functions      <- Preparation
!        - substitute_g_list, - symmetric_G_points_using_g_list
!        - counting_error_points, - symmetric_G_points
!  13.   m_pwBS_sphrp_l                <- Initial_Electronic_Structure
!  14.   m_pwBS_sphrp2                 <- m_Chargek_Density, m_ES_Intgr_VlhxcQlm,
!                                         m_Force, m_Stress, m_XC_Potential
!  15.   m_pwBS_sphrp2_diff            <- m_Stress, m_XC_Potential, m_epc_potential
!  16.   m_pwBS_kinetic_energies       <- m_ES_WF_by_Davidson, m_ES_WF_by_RMM,
!                                         m_ES_WF_by_SDorCG, m_ES_WF_by_submat,
!                                         m_Electronic_Structure
!  17.   m_pwBS_pstrn_kinetic_energies <- m_Positron_Wave_Functions
!  18.   m_pwBS_find_min_max_G         <- m_ES_WF_by_RMM
!  19.   m_pwBS_decide_cutoff_mix      <- Preparation
!  20.   m_pwBS_wd_ngabc_etc           <- Postprocessing
!  21.   m_pwBS_increase_kg1           <- Preparation

contains
  subroutine m_pwBS_alloc_ylm_l()
    if(nel_Ylm >= 1) then
       if(.not.allocated(ylm_l)) allocate(ylm_l(ista_kngp:iend_kngp,nel_Ylm))
    end if
  end subroutine m_pwBS_alloc_ylm_l

  subroutine m_pwBS_dealloc_ylm_l()
    if(allocated(ylm_l)) deallocate(ylm_l)
  end subroutine m_pwBS_dealloc_ylm_l

  subroutine m_pwBS_alloc_igpo()
    integer :: i,j,k
    if(.not.allocated(igpo)) then
       i = n_rGsv(1); j = n_rGsv(2); k = n_rGsv(3)
       allocate(igpo(-i:i,-j:j,-k:k))
    end if
  end subroutine m_pwBS_alloc_igpo

  subroutine m_pwBS_dealloc_igpo()
    if(allocated(igpo)) deallocate(igpo)
  end subroutine m_pwBS_dealloc_igpo

  function get_kg1_ext() result (res)
    integer :: res
    if(kg1 < kg1_ext) then
      res = kg1
    else
      res = kg1_ext
    endif
  end function get_kg1_ext

  subroutine alloc_nbase
    integer :: ik
    logical :: tf
    if(.not.allocated(nbase)) then
       allocate(nbase(kg1_ext,kv3))
!!$       write(6,'(" nbase is allocated, kg1_ext = ",i8,"<<alloc_nbase>>")') kg1_ext
    else
!!$       write(6,'(" nbase has already been allocated, kg1_ext = ",i8,"<<alloc_nbase>>")') kg1_ext
    end if
    tf = .false.
    do ik = 1, kv3
       if(k_symmetry(ik) == GAMMA .or. k_symmetry(ik) == GAMMA_base_symmetrization) then
          tf = .true.
       end if
    end do
    if(tf) then
       if(.not.allocated(nbase_gamma)) then
          allocate(nbase_gamma(kg_gamma,2)); nbase_gamma = -1
       end if
! === Debug by Intel "-check all" option! by T.Kato 2011/03/28 =================
!   end if
    else
       if(.not.allocated(nbase_gamma)) then
          allocate(nbase_gamma(1,1)); nbase_gamma = -1
       end if
    end if
! ==============================================================================
  end subroutine alloc_nbase

  subroutine reduce_nbase()
    real(kind=DP), allocatable, dimension(:,:) :: nbase_tmp
    integer :: i, k
    allocate(nbase_tmp(kg1,kv3)); nbase_tmp = 0
    do k = 1, kv3
       do i = 1, iba(k)
          nbase_tmp(i,k) = nbase(i,k)
       end do
    end do
    deallocate(nbase)
    allocate(nbase(kg1,kv3))
    nbase(:,:) = nbase_tmp(:,:)

    deallocate(nbase_tmp)
  end subroutine reduce_nbase

  subroutine alloc_iba
    if(.not.allocated(iba)) then
       allocate(iba(kv3))
       if(ipri >= 3) write(nfout,'(" iba is allocated")')
    else
       if(ipri >= 3) write(nfout,'(" iba is not allocated")')
    end if
  end subroutine alloc_iba

  subroutine G_vectors_alloc
    allocate(ngabc(kgp,3))
  end subroutine G_vectors_alloc

  subroutine m_pwBS_alloc_nbmat_and_iba2()
    if(nmatsz <= 0) call phase_error_with_msg(nfout,' nmatsz <= 0 (alloc_nbmat_and_iba2)',__LINE__,__FILE__)
    if(ipri >= 2 ) write(6,'(" nmatsz = ",i5," kv3 = ",i5," (alloc_nbmat_and_iba2)")') nmatsz,kv3
!!$    allocate(nbmat(nmatsz,ista_k:iend_k))
!!$    allocate(iba2(ista_k:iend_k))
    allocate(nbmat(nmatsz,kv3)); nbmat = 0.d0
    allocate(nbmat2(nmatsz,kv3)); nbmat2 = 0.d0
    allocate(iba2(kv3))
  end subroutine m_pwBS_alloc_nbmat_and_iba2

  subroutine m_pwBS_dealloc_nbmat_and_iba2()
    if(allocated(nbmat)) deallocate(nbmat)
    if(allocated(nbmat2)) deallocate(nbmat2)
    if(allocated(iba2))  deallocate(iba2)
  end subroutine m_pwBS_dealloc_nbmat_and_iba2

  subroutine m_pwBS_alloc_ngpt_igfp_gr
    allocate(ngpt_l(ista_kngp:iend_kngp,nopr+af)); ngpt_l = 0.d0
    allocate(igfp_l(ista_kngp:iend_kngp));      igfp_l = 0
#ifdef _MPIFFT_
    allocate(igfp_nonpara(ista_kngp:iend_kngp)); igfp_nonpara = 0
    allocate(igfp_l_c(ista_kngp:iend_kngp));    igfp_l_c = 0
#endif
    allocate(igf(kg));                          igf    = 0
    allocate(gr_l(ista_kngp:iend_kngp));        gr_l   = 0.d0
    ! ------- Positron start
    if(sw_positron /= OFF) then
       allocate(igf_pstrn(kg_pwf)) ; igf_pstrn = 0
    end if
    ! ------- Positron end

! ========== KT_add ========================== 13.0F
    if (sw_hybrid_functional /= OFF .and. use_fft_exx ) then
       allocate(igf_exx(kg_exx)) ; igf_exx = 0
    end if
    if (sw_hybrid_functional /= OFF) then
       allocate(igfp_exx(kgp_exx)) ; igfp_exx = 0
    endif
! ============================================ 13.0F

  end subroutine m_pwBS_alloc_ngpt_igfp_gr

  subroutine m_pwBS_dealloc_ngpt_igfp_gr
    if(allocated(gr_l)) deallocate(gr_l)
    if(allocated(ngpt_l)) deallocate(ngpt_l)
    if(allocated(igfp_l)) deallocate(igfp_l)
#ifdef _MPIFFT_
    if(allocated(igfp_nonpara)) deallocate(igfp_nonpara)
    if(allocated(igfp_l_c)) deallocate(igfp_l_c)
#endif
    if(allocated(igf)) deallocate(igf)
    ! ------- Positron start
    if(sw_positron /= OFF) deallocate(igf_pstrn)
    ! ------- Positron end

! ========== KT_add ========================== 13.0F
    if (sw_hybrid_functional /= OFF .and. use_fft_exx ) deallocate(igf_exx)
    if (sw_hybrid_functional /= OFF ) deallocate(igfp_exx)
! ============================================ 13.0F

  end subroutine m_pwBS_dealloc_ngpt_igfp_gr

  subroutine m_pwBS_mat_for_each_WF()
    integer          :: i,j,k,kk, icycle, istart, iend, ic, ig, max_elements, icolumn &
         &            , kk_gamma, n1,n2,n3
    real(kind=DP)    :: ga,gb,gc,grvv
    integer          :: id_sname = -1
    call tstatc0_begin('m_pwBS_mat_for_each_WF ',id_sname)

    do i = 1, kv3
       kk = 0
       if(k_symmetry(i) == GAMMA .or. k_symmetry(i) == GAMMA_base_symmetrization) then
          kk_gamma = 0
          do j = 1, kg
             n1 = ngabc(j,1); n2 = ngabc(j,2); n3 = ngabc(j,3)
             grvv = dsqrt(ttr(1)*n1*n1 + ttr(2)*n2*n2 + ttr(3)*n3*n3 &
                  &   +   ttr(4)*n1*n2 + ttr(5)*n2*n3 + ttr(6)*n3*n1)
             if(grvv <= gmaxs) then
                kk = kk + 1
                nbmat(kk,i) = j
                if(n3 >  SmallestPositiveNumber) then
                   kk_gamma = kk_gamma + 1
                else if(n3 == 0) then
                   if((n1 >= 0 .and. n2 >= 0 ) .or. (n1 < 0 .and. n2 > 0)) then
                      kk_gamma = kk_gamma + 1
                   end if
                end if
             end if
          end do
          kg2_gamma = kk_gamma
       else
          do j = 1, kg
             ga = vkxyz(i,1,BUCS) + ngabc(j,1)
             gb = vkxyz(i,2,BUCS) + ngabc(j,2)
             gc = vkxyz(i,3,BUCS) + ngabc(j,3)
             grvv = dsqrt(ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
                  &   +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga)
             if(grvv <= gmaxs) then
                kk = kk + 1
                nbmat(kk,i) = j
             end if
          end do
       end if
       iba2(i) = kk
       if(iba2(i) < neg) then
          write(nfout,'(" iba2(",i4,") < neg = ",i4)') i,neg
          call phase_error_with_msg(nfout,' m_pwBS_mat_for_each_WF',__LINE__,__FILE__)
       end if
       if(k_symmetry(i) == GAMMA .or. k_symmetry(i) == GAMMA) then
          do j=1,iba2(i)
             nbmat2(j,i) = nbmat(j,i)
!!$             if(nbmat(j,i) == nbase_gamma(1,1)) nbmat2(j,i) = 1
!!$             if(nbmat2(j,i) == 0) then
!!$                do k=2, kg_gamma
!!$                   if(nbmat(j,i) == nbase_gamma(k,1)) then
!!$                      nbmat2(j,i) = nbase_gamma(k,1)
!!$                      exit
!!$                   else if(nbmat(j,i) == nbase_gamma(k,2)) then
!!$                      nbmat2(j,i) = nbase_gamma(k,2)
!!$                      exit
!!$                   end if
!!$                end do
!!$             end if
          end do
       else
          do j = 1, iba2(i)
             do k=1,iba(i)
                if(nbmat(j,i) == nbase(k,i)) then
                   nbmat2(j,i) = k
                   exit
                end if
             end do
          end do
       end if
       if(ipri >= 2) then
          ! --- nbmat, nbmat2 ---
          max_elements = iba2(i)
          icolumn = 10
          write(nfout,'(" !pwbs -- nbmat, nbase --")')
          write(nfout,'(" !pwbs k = ",i4," iba2= ",i6)') i,iba2(i)
          icycle = ceiling(dble(min(max_elements,iba2(i)))/icolumn)
          istart = 1
          do ic = 1, icycle
             iend = min(istart+icolumn-1,max_elements,iba2(i))
             if(ipri >= 3 .or. (ipri>=2 .and. (ic <= 3 .or. ic >= icycle-2))) then
                write(nfout,'(" ! nbmat  ",10i10)') (nbmat(ig,i),ig=istart,iend)
                write(nfout,'(" ! nbmat2 ",10i10)') (nbmat2(ig,i),ig=istart,iend)
             end if
             istart = iend+1
          end do

          ! --- nbase ---
!!$          max_elements = kg1_ext
          max_elements = iba(i)
          icolumn = 10
          write(nfout,'(" !pwbs -- nbase --")')
!!$          icycle = ceiling(dble(min(max_elements,kg1_ext))/icolumn)
          icycle = ceiling(dble(min(max_elements,iba(i)))/icolumn)
          istart = 1
          do ic = 1, icycle
!!$             iend = min(istart+icolumn-1,max_elements,kg1_ext)
             iend = min(istart+icolumn-1,max_elements,iba(i))
             if(ipri >= 2 .or. (ipri>=1 .and. (ic <= 3 .or. ic >= icycle-2))) then
                write(nfout,'(" ! nbase ",10i10)') (nbase(ig,i),ig=istart,iend)
             end if
             istart = iend+1
          end do

          ! --- nbase_gamma ---
          if(kg_gamma > 0 .and. (k_symmetry(i) == GAMMA .or. k_symmetry(i) == GAMMA_base_symmetrization)) then
             max_elements = kg_gamma
             icolumn = 6
             write(nfout,'(" !pwbs -- nbase_gamma --")')
             write(nfout,'(" !pwbs k = ",i4," kg_gamma = ",i6)') i,kg_gamma
             if(kg2_gamma > 0) write(nfout,'(" !pwbs kg2_gamma = ",i6)') kg2_gamma
             icycle = ceiling(dble(min(max_elements,kg_gamma))/icolumn)
             istart = 1
             do ic = 1, icycle
                iend = min(istart+icolumn-1,max_elements,kg_gamma)
                if(ipri >= 2 .or. (ipri>=1 .and. (ic <= 3 .or. ic >= icycle-2))) then
                   write(nfout,'(" !nbase_gamma1 (",i9,"-",i9,") ",6i10)') &
                        &  istart, iend, (nbase_gamma(ig,1),ig=istart,iend)
                   write(nfout,'(" !nbase_gamma2 (",i9,"-",i9,") ",6i10)') &
                        &  istart, iend, (nbase_gamma(ig,2),ig=istart,iend)
                end if
                istart = iend+1
             end do
          end if
       end if
    end do

    call tstatc0_end(id_sname)
  end subroutine m_pwBS_mat_for_each_WF

  recursive subroutine m_pwBS_for_each_WF(preallocation)
!   Revised on 19th September 2006 by T. Yamasaki
!      A finding algorism of nbase_gamma(*,2) is revised for the case that
!     ngshell_range_G's are not properly prepared.
    logical, intent(in)       :: preallocation

    integer :: jjt, i, jj, jj_gamma, j, n1, n2, n3, j2, jg, k, iloop &
         &   , j2s, j2e, jg1, jg2, kg_gamma_part, max_elements, icolumn, icycle &
         &   , istart, ic, iend, ig, idelta
    real(kind=DP)    :: ga,gb,gc, grvv, grvt, grvmx, length_start, length
    logical :: found
    integer, allocatable, dimension(:)   :: ngshell_range_G !d(ngshell_GAMMA)
    integer, allocatable, dimension(:)   :: ngbshell !d(kg0)
    integer, parameter :: CRITICAL_VECTOR_LENGTH = 10000

    integer          :: id_sname = -1
                                                     __TIMER_SUB_START(1221)
    call tstatc0_begin('m_pwBS_for_each_WF ',id_sname)

    if(ipri >= 2) then
       if(     preallocation) write(nfout,*) ' ! preallocation = .true.'
       if(.not.preallocation) write(nfout,*) ' ! preallocation = .false.'
    end if

    if(allocated(iba))then
       if(ipri >= 3) write(nfout,'(" iba is allocated <<m_pwBS_for_each_WF>>")')
    else
       if(ipri >= 3) write(nfout,'(" iba is not allocated <<m_pwBS_for_each_WF>>")')
    end if

    call alloc_iba()
    if(.not.preallocation) then
       !!$call alloc_iba()
       if(kg1_ext == 0) call m_pwBS_for_each_WF(.true.)
       call alloc_nbase()
!!$       write(6,'(" alloc_nbase")')
    else
!!$       write(6,'(" alloc_nbase is not called")')
    end if

    if(preallocation .and. ipri >= 1) write(nfout,'(/" <<< m_pwBS_for_each_WF >>>")')

    call getttr(rltv,ttr)

    if(preallocation) kg1_ext = 0; if(preallocation) nbmx = 0
    if(preallocation) kg1 = 0;

    if(preallocation) then
       jjt = 0
       grvt = 0.d0
    else if(.not.preallocation) then
       nbase = 0
    end if

    do i = 1, kv3
       if(k_symmetry(i) == GAMMA .or. k_symmetry(i) == GAMMA_base_symmetrization) then
          if(preallocation) then
             grvmx = 0.d0
             ngshell_GAMMA = 1
             length_start = 0.d0
             jj_gamma = 0
!CDIR NOVECTOR
             g_search_loop:do j = 1, kg0
                n1 = ngabc(j,1)
                n2 = ngabc(j,2)
                n3 = ngabc(j,3)
                grvv = dsqrt(ttr(1)*n1*n1 + ttr(2)*n2*n2 + ttr(3)*n3*n3 &
                     &   +   ttr(4)*n1*n2 + ttr(5)*n2*n3 + ttr(6)*n3*n1)
!!$                grvv = ( ttr(1)*n1*n1 + ttr(2)*n2*n2 + ttr(3)*n3*n3 &
!!$                     & + ttr(4)*n1*n2 + ttr(5)*n2*n3 + ttr(6)*n3*n1)
                if(grvv > length_start + DELTA10) then
!!$                if(grvv > length_start + DELTA07) then
                   ngshell_GAMMA = ngshell_GAMMA+1
                   length_start = grvv
                end if
                if(grvv > grvmx) grvmx=grvv
                if(n3 >  SmallestPositiveNumber) then
                   jj_gamma = jj_gamma + 1
                else if(n3 == 0) then
                   if((n1 >= 0 .and. n2 >= 0 ) .or. (n1 < 0 .and. n2 > 0)) then
                      jj_gamma = jj_gamma + 1
                   end if
                end if
             end do g_search_loop
!!$             grvmx = dsqrt(grvmx)

             if(ipri >= 1 ) write(nfout,'(" !pwbs ngshell_GAMMA = ",i8)') ngshell_GAMMA
             jj = kg0
!!$             kg_gamma = jj_gamma
             kg_gamma = (kg0+1)/2

             if(nbmx < kg0) nbmx = kg0
             iba(i) = kg0
             if(ipri >= 1) write(nfout,'(" !# iba(",i4,") = ",i8," (preallocation)<<m_pwBS_for_each_WF>>")')&
                  & i, iba(i)
          else if(.not.preallocation) then
             allocate(ngshell_range_G(ngshell_GAMMA))
             allocate(ngbshell(kg0))
             jj = 0
             jg = 1
             ngshell_range_G(jg) = 1
             length_start = 0.d0
             g_search_loop2:do j = 1, kg0
                n1 = ngabc(j,1)
                n2 = ngabc(j,2)
                n3 = ngabc(j,3)
                grvv = dsqrt(ttr(1)*n1*n1 + ttr(2)*n2*n2 + ttr(3)*n3*n3 &
                     &   +   ttr(4)*n1*n2 + ttr(5)*n2*n3 + ttr(6)*n3*n1)
                nbase(j,i)=j
                if(grvv > length_start + DELTA10) then
!!$                if(grvv > length_start + DELTA07) then
                   jg = jg+1
                   ngshell_range_G(jg) = j
                   length_start = grvv
                end if
                ngbshell(j) = jg
                if(n3 >  SmallestPositiveNumber) then
                   jj = jj + 1
                   nbase_gamma(jj,1) = j
                else if(n3 == 0) then
                   if((n1 >= 0 .and. n2 >= 0 ) .or. (n1 < 0 .and. n2 > 0)) then
                      jj = jj + 1
                      nbase_gamma(jj,1) = j
                   end if
                end if
             end do g_search_loop2

             do j = 1, ngshell_GAMMA-1
                if(ngshell_range_G(j+1) < ngshell_range_G(j)) then
                   write(nfout,'("!! ngshell_range order illegal <<m_pwBS_for_each_WF>>")')
                   write(nfout,'("!! j = ",i5," ngshell_range_G(j+1), ngshell_range_G(j) = ",2i8)') &
                        & ngshell_range_G(j+1), ngshell_range_G(j)
                   call phase_error_with_msg(nfout,' ngshell_range_G order is illegl << m_pwBS_for_each_WF>>', &
                   __LINE__,__FILE__)
                end if
             end do
             if(ipri >= 1) call wd_ngshell_range()

             kg_gamma_part = ceiling(dble(kg_gamma)/npes)
             if(kg_gamma_part > CRITICAL_VECTOR_LENGTH) then
                if(kg_gamma_part == 0) kg_gamma_part = 1
                jg1 = kg_gamma_part*mype + 1
                jg2 = kg_gamma_part*(mype+1)
                if(jg2 > kg_gamma) jg2 = kg_gamma
                if(ipri >= 2) write(nfout,'(" !pwBS jg1, jg2 = ", 2i5)') jg1, jg2
             else
                jg1 = 1; jg2 = kg_gamma
             end if

             nbase_gamma(1:kg_gamma,2) = 0
             if(jg1 == 1) then
                nbase_gamma(1,2) = -1
                jg1 = 2
             end if
             do jg = jg1, jg2
                j = nbase_gamma(jg,1)
                k = ngbshell(j)
                n1 = ngabc(j,1); n2 = ngabc(j,2); n3 = ngabc(j,3)
                j2s = ngshell_range_G(k)
                if(k<=ngshell_GAMMA-1) then; j2e = ngshell_range_G(k+1)-1; else; j2e = iba(i); end if ! j2e
                   idelta = 0
1003               found = .false.
                G_search: do j2 = j2s, j2e
                   if(ngabc(j2,1) == -n1 .and. ngabc(j2,2) == -n2 .and. ngabc(j2,3) == -n3) then
                      nbase_gamma(jg,2) = j2
                      if(j2 == j) nbase_gamma(jg,2) = -1
                      found = .true.
                      exit G_search
                   end if
                end do G_search
                if(.not.found) then
                   idelta = idelta + 1
                   if(k-idelta >= 1) then
                      j2s = ngshell_range_G(k-idelta)
                   else
                      j2s = 1
                   end if
                   if(k+idelta <= ngshell_GAMMA-1) then
                      j2e = ngshell_range_G(k+idelta+1)-1
                   else
                      j2e = iba(i)
                   end if
                   goto 1003
                end if
             end do
             if(npes > 1) then
                if(kg_gamma_part > CRITICAL_VECTOR_LENGTH) then
                   call mpi_allreduce(nbase_gamma(1:kg_gamma,2),ngbshell,kg_gamma,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
                   nbase_gamma(1:kg_gamma,2) = ngbshell(1:kg_gamma)
                end if
             end if

             if(k_symmetry(i) == GAMMA) then
                do jg = 1, kg_gamma
                   nbase(jg,i) = nbase_gamma(jg,1)
                end do
                do jg = kg_gamma+1, kg1_ext
                   nbase(jg,i) = 0
                end do
                iba(i) = kg_gamma
             end if
             deallocate(ngbshell,ngshell_range_G)
          end if
          if(k_symmetry(i) == GAMMA) then
             if(preallocation .and. kg_gamma > kg1 ) kg1 = kg_gamma
          else
!!$             if(preallocation .and. jj > kg1) kg1 = jj
             if(preallocation .and. kg0 > kg1) kg1 = kg0
          end if
       else
          if(preallocation) grvmx = 0.d0
          jj = 0
          do j = 1, kg
             ga = vkxyz(i,1,BUCS) + ngabc(j,1)
             gb = vkxyz(i,2,BUCS) + ngabc(j,2)
             gc = vkxyz(i,3,BUCS) + ngabc(j,3)
             grvv = dsqrt(ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
                  &   +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga)
             if(grvv <= gmax) then
                jj = jj + 1
                if(preallocation .and. nbmx < j ) nbmx = j
                if(preallocation) then
                   if(grvv > grvmx) grvmx=grvv
                end if

! ===== KT_mod ===================== 13.0S [OK??]
!
!  The values of vkxyz when prealloc=true are different
!  from those values when prealloc=false
!                     ( see Preparation.F90, ekmode=on )
!
!  Thus, "jj" (when prealloc=false) may be over "kg1_ext" set when prealloc=true.
!
! ------------------------------
!                if(.not.preallocation) nbase(jj,i)=j
                if(.not.preallocation) then
                   if ( jj <= kg1_ext ) nbase(jj,i)=j
                endif
! ================================== 13.0S

             end if
          enddo
          if(preallocation) then
             if(jj > kg1) kg1 = jj
             iba(i) = jj
             jjt    = jjt + jj
             grvt   = grvt + grvmx*qwgt(i)
          end if
       end if

       if(ipri >= 3) then
          if(preallocation) then
             write(nfout,'(" !# K= ",i4," iba= ",i6," GRV= ",F12.8)') i,iba(i),grvmx
          else
             write(nfout,'(" !# K= ",i4," iba= ",i6," nbase(",i6,",",i4,")= ",i7)') &
                  & i,iba(i),iba(i),i,nbase(iba(i),i)
             if(iba(i) > kg1_ext) then
                write(nfout,'(" ! iba(",i3,") > kg1_ext (= ",i9)') i,kg1_ext
                call phase_error_with_msg(nfout,' iba(i) > kg1_ext <<m_pwBS_for_each_WF>>',__LINE__,__FILE__)
             else if(iba(i) == 0) then
                write(nfout,'(" ! iba(",i3,") = ",i6)') i,iba(i)
                call phase_error_with_msg(nfout,' iba(i) == 0 <<m_pwBS_for_each_WF>>',__LINE__,__FILE__)
             end if
             if(ipri >= 2) then
                if(kg_gamma > 0 .and.(k_symmetry(i) == GAMMA &
                     &              .or. k_symmetry(i) == GAMMA_base_symmetrization)) then
                   max_elements = kg_gamma
                   icolumn = 10
                   write(nfout,'(" !pwbs -- nbase_gamma --")')
                   write(nfout,'(" !pwbs k = ",i4," kg_gamma = ",i6)') i,kg_gamma
                   if(kg2_gamma > 0) write(nfout,'(" !pwbs kg2_gamma = ",i6)') kg2_gamma
                   icycle = ceiling(dble(min(max_elements,kg_gamma))/icolumn)
                   istart = 1
                   do ic = 1, icycle
                      iend = min(istart+icolumn-1,max_elements,kg_gamma)
                      if(ipri >= 2 .or. (ipri>=1 .and. (ic <= 3 .or. ic >= icycle-2))) then
                         write(nfout,'(" ! nbase_gamma 1 ",10i10)') (nbase_gamma(ig,1),ig=istart,iend)
                         write(nfout,'(" ! nbase_g1 (nx) ",10i10)') (ngabc(nbase_gamma(ig,1),1),ig=istart,iend)
                         write(nfout,'(" ! nbase_g1 (ny) ",10i10)') (ngabc(nbase_gamma(ig,1),2),ig=istart,iend)
                         write(nfout,'(" ! nbase_g1 (nz) ",10i10)') (ngabc(nbase_gamma(ig,1),3),ig=istart,iend)
                         write(nfout,'(" ! nbase_gamma 2 ",10i10)') (nbase_gamma(ig,2),ig=istart,iend)
                         write(nfout,'(" ! nbase_g2 (nx) ",10i10)') (ngabc(nbase_gamma(ig,2),1),ig=istart,iend)
                         write(nfout,'(" ! nbase_g2 (ny) ",10i10)') (ngabc(nbase_gamma(ig,2),2),ig=istart,iend)
                         write(nfout,'(" ! nbase_g2 (nz) ",10i10)') (ngabc(nbase_gamma(ig,2),3),ig=istart,iend)
                      end if
                      istart = iend+1
                   end do
                end if
             end if
          end if
       end if
       if(preallocation .and. jj > kg1_ext ) kg1_ext = jj
    end do ! k-point loop

    if(.not.preallocation) then
       if(kg1 < kg1_ext) then
          call reduce_nbase()
       end if
    end if

    if(ipri >= 1) then
       if(preallocation) then
          write(nfout,130) kg1_ext,kg1,nbmx
          write(nfout,140) jjt,grvt
          write(nfout,'(" !# pwbs kg_gamma = ",i8)') kg_gamma
130       format(' ','!# ** kg1_ext, kg1, nbmx (=matrix size) = ',3i8)
140       format(' ','!#  JJT(=sum(iba)) = ',I8,' MEAN GRV = ',f12.8)
          write(nfout,*)
       else
          if(ipri >= 2 .or. (ipri>=1 .and. kv3 <=10)) then
             do i = 1, kv3
                write(nfout,'(" ! iba(",i6,") = ",i6,",  nbase(",i6,",",i6,") = ",i7)') i, iba(i),iba(i),i,nbase(iba(i),i)
             end do
          else
             if(kv3 > 10) then
                do iloop = 1, 2
                   do i = (kv3-4)*(iloop-1)+1, (kv3-4)*(iloop-1)+4
                      write(nfout,'(" ! iba(",i6,") = ",i6,",  nbase(",i6,",",i6,") = ",i7)') i, iba(i),iba(i),i,nbase(iba(i),i)
                   end do
                   if(iloop == 1) write(nfout,'(" ! ......")')
                end do
             end if
          end if
       end if
    end if

    call tstatc0_end(id_sname)
                                                     __TIMER_SUB_STOP(1221)
  contains
    subroutine wd_ngshell_range
      integer :: jg, j,n1,n2,n3
      real(kind=DP) :: gr1, gr2

      write(nfout,'(" !pwbs ngshell_range_G")')
      write(nfout,'(" !pwbs  ngshell_G, range1 - range2, nelement, gr1, gr2")')
      do jg = 1, ngshell_GAMMA
         if(ipri >= 2 .or.&
              & (ipri==1 .and. (jg <= 5 .or. jg >= ngshell_GAMMA - 4))) then
            j = ngshell_range_G(jg)
            n1 = ngabc(j,1);  n2 = ngabc(j,2);  n3 = ngabc(j,3)
            gr1 = dsqrt(ttr(1)*n1*n1 + ttr(2)*n2*n2 + ttr(3)*n3*n3 &
                 &  +   ttr(4)*n1*n2 + ttr(5)*n2*n3 + ttr(6)*n3*n1)
            if(jg == ngshell_GAMMA) then
               j = iba(i)
            else
               j = ngshell_range_G(jg+1)-1
            end if
            n1 = ngabc(j,1);  n2 = ngabc(j,2);  n3 = ngabc(j,3)
            gr2 = dsqrt(ttr(1)*n1*n1 + ttr(2)*n2*n2 + ttr(3)*n3*n3 &
                 &  +   ttr(4)*n1*n2 + ttr(5)*n2*n3 + ttr(6)*n3*n1)
            write(nfout,'(" !pwbs ",i8, " [",i8," : ",i8," ] ", i8, 2f20.8)') &
                 & jg, ngshell_range_G(jg),j &
                 & ,j-ngshell_range_G(jg)+1,gr1, gr2
         else if(ipri == 1 .and. jg == 6) then
            write(nfout,'(" !pwbs  ..........")')
         end if
      end do
      if(ipri >= 3) then
         write(nfout,'(" !pwBS -- ngbshell --")')
         write(nfout,'(5("(",i8,i6,")"))') (j,ngbshell(j),j=1,iba(i))
      end if
    end subroutine wd_ngshell_range
  end subroutine m_pwBS_for_each_WF

  ! ------- Positron start
  recursive subroutine m_pwBS_positronWF()

    integer          :: j
    real(kind=DP)    :: ga,gb,gc, grvv
    integer          :: id_sname = -1
    call tstatc0_begin('m_pwBS_positronWF ',id_sname)


    call getttr(rltv,ttr)

    kg1_pwf = 0

    do j = 1, kg_pwf
       ga = ngabc(j,1)
       gb = ngabc(j,2)
       gc = ngabc(j,3)
       grvv = dsqrt(ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
            &     + ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga )
       if(grvv <= gmax_positron) then
          kg1_pwf = kg1_pwf + 1
       end if
    enddo

    if(ipri>=1) write(nfout,'(" !pwbs kg1_pwf = ",i7)') kg1_pwf
    call tstatc0_end(id_sname)
  end subroutine m_pwBS_positronWF
  ! ------- Positron end

! =================== KT_add ====================== 13.0F
  recursive subroutine m_pwBS_exxWF()

    integer          :: j
    real(kind=DP)    :: ga,gb,gc, grvv
    integer          :: id_sname = -1
    call tstatc0_begin('m_pwBS_exxWF ',id_sname)

    call getttr(rltv,ttr)
    kg1_exx = 0

    do j = 1, kg_exx
       ga = ngabc(j,1);    gb = ngabc(j,2);    gc = ngabc(j,3)
       grvv = dsqrt(ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
            &     + ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga )
       if(grvv <= gmax_exx) then
          kg1_exx = kg1_exx + 1
       end if
    enddo

    if(ipri>=1) write(nfout,'(" !pwbs kg1_exx = ",i7)') kg1_exx
    call tstatc0_end(id_sname)

  end subroutine m_pwBS_exxWF

  recursive subroutine m_pwBS_exxCD()

    integer          :: j
    real(kind=DP)    :: ga,gb,gc, grvv
    integer          :: id_sname = -1
    call tstatc0_begin('m_pwBS_exxCD ',id_sname)

    call getttr(rltv,ttr)
    kg1p_exx = 0

    do j = 1, kgp_exx
       ga = ngabc(j,1);    gb = ngabc(j,2);    gc = ngabc(j,3)
       grvv = dsqrt(ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
            &     + ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga )
       if(grvv <= gmaxp_exx) then
          kg1p_exx = kg1p_exx + 1
       end if
    enddo

    if(ipri>=1) write(nfout,'(" !pwbs kg1p_exx = ",i7)') kg1p_exx
    call tstatc0_end(id_sname)

  end subroutine m_pwBS_exxCD

! ================================================ 13.0F

  subroutine m_pwBS_set_gmaxs(n_matrix_size,gmaxs_given)
    integer, intent(in) ::       n_matrix_size
    real(kind=DP), intent(in) :: gmaxs_given

    real(kind=DP)    :: grvv,ga,gb,gc
    integer          :: negg,i,j,kk,nmatsz_min,counter

!!$    if(intzaj == by_matrix_diagon) then
       call check_gmax(gmaxs_given)
       gmaxs = gmaxs_given
       if(gmaxs < SmallestPositiveNumber*1.d5 .and. n_matrix_size == 0) then
          if(kg < 15000) then
             negg = kg/15
          else if(neg < 1000) then
             negg = 1000
          else
             negg = neg
          end if
          if(negg < neg) negg = neg+20
          gmaxs = 2.05*gmax*(dble(neg)/dble(kg))**(1.d0/3.d0)
          if(printable) write(nfout,*) ' ! Both gmaxs and n_matrix_size were not given'
       else if(n_matrix_size > 0) then
          if(n_matrix_size < neg) then
             negg = neg + 20
          else
             negg = n_matrix_size
          end if
          gmaxs = 2.05*gmax*(dble(neg)/dble(kg))**(1.d0/3.d0)
          if(printable) write(nfout,*) ' ! n_matrix_size is given'
       else
          negg = kg * (gmaxs/(2.05*gmax))**3
          if(negg < neg) then
             negg = neg + 20
          end if
!!$          write(6,*) ' ! gmaxs is given'
       end if
!!$       write(nfout,'(" ! negg  = ",i8)') negg
!!$       write(nfout,'(" ! gmaxs = ",f8.4)') gmaxs

       counter = 0
1      nmatsz = 0
       nmatsz_min = kg
       n_rGsv = 0
       counter = counter+1
       do i = 1, kv3
          kk = 0
          do j = 1, kg
             ga = vkxyz(i,1,BUCS) + ngabc(j,1)
             gb = vkxyz(i,2,BUCS) + ngabc(j,2)
             gc = vkxyz(i,3,BUCS) + ngabc(j,3)
             grvv = dsqrt(ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
                  &   +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga)
             if(grvv <= gmaxs) then
                kk = kk + 1
                if(n_rGsv(1) < abs(ngabc(j,1))) n_rGsv(1) = abs(ngabc(j,1))
                if(n_rGsv(2) < abs(ngabc(j,2))) n_rGsv(2) = abs(ngabc(j,2))
                if(n_rGsv(3) < abs(ngabc(j,3))) n_rGsv(3) = abs(ngabc(j,3))
             end if
          end do
          if(kk > nmatsz) nmatsz = kk
          if(kk < nmatsz_min) nmatsz_min = kk
!!$          if(nmatsz_min < neg) exit
       end do
!!$       if(nmatsz < negg) then

!!$       write(nfout,'(" ! negg = ",i8," kg = ",i8)') negg,kg
       if(ipri >= 1) then
          write(nfout,'(" ! nmatsz = ",i8," nmatsz_min = ",i8)') nmatsz, nmatsz_min
          write(nfout,'(" ! gmaxs = ",f20.8)') gmaxs
       end if
       if(nmatsz_min < negg) then
          gmaxs = gmaxs*1.05
          if(ipri >= 1) then
             write(nfout,'(" ---- trial ",i4," ----")') counter
             write(nfout,'(" ! gmaxs is enlarged")')
          end if
          goto 1
       end if

       if(gmaxs > gmax) then
          if(ipri >= 1) write(nfout,'(" gmaxs(=",f10.4,") > gmax(=",f10.4,")")') gmaxs,gmax
          gmaxs = gmax
          if(ipri >= 1) write(nfout,'(" gmaxs is set gmax")')
       end if
       n_rGsv = 2*n_rGsv
       if(ipri >= 1) then
          write(nfout,'(" n_rGsv(1:3) = ",3i8)') n_rGsv(1:3)
          write(nfout,'(" nmatsz (m_pwBS_set_gmaxs) = ",i8)') nmatsz
       end if

  contains
    subroutine check_gmax(givenvalue)
      real(kind=DP),intent(in) :: givenvalue
      if(givenvalue < SmallestPositiveNumber*1.d5) then
         if(printable) write(nfout,'(" gmax_given is not given")')
!!$         stop ' illegal value of gmax (m_pwBS_set_gmaxs)'
      end if
    end subroutine check_gmax
  end subroutine m_pwBS_set_gmaxs

  subroutine m_pwBS_assume_G_rhombohedron
    real(kind=DP)   :: gmax2, gmax2_pstrn

! ============== KT_add ============================= 13.0F
    real(kind=DP)   :: gmax2_exx
! =================================================== 13.0F

    gmax2 = gmax*2
    call size_of_Gvector_rhombohedron  ! ->n_rGv, n_rGpv
    if(ipri >=1) call check_the_size_of_rhombohedron
  contains
    subroutine size_of_Gvector_rhombohedron
      integer i
      real(kind=DP)   :: alen
      do i = 1, 3
         alen = ( rltv(1,i)*altv(1,i) + rltv(2,i)*altv(2,i) + rltv(3,i)*altv(3,i) ) &
              & /dsqrt(altv(1,i)*altv(1,i) + altv(2,i)*altv(2,i) + altv(3,i)*altv(3,i))
         n_rGv(i)  = abs(int(gmax2/alen)) + 1
         n_rGpv(i) = abs(int(gmaxp/alen)) + 1
         n_rGpv_reduced(i) = abs(int(gmaxp_reduced/alen)) + 1
  ! ------- Positron start
         if(sw_positron /= OFF) then
            gmax2_pstrn = gmax_positron * 2
            n_rGv_pstrn(i) = abs(int(gmax2_pstrn/alen)) + 1
         end if
  ! ------- Positron end

! ============== KT_add ============================= 13.0F
         if ( sw_hybrid_functional /= OFF .and. use_fft_exx ) then
            gmax2_exx = gmax_exx * 2
            n_rGv_exx(i) = abs(int(gmax2_exx/alen)) + 1
         end if
         if ( sw_hybrid_functional /= OFF ) then
            n_rGpv_exx(i) = abs(int(gmaxp_exx/alen)) + 1
         endif
! =================================================== 13.0F

      enddo
    end subroutine size_of_Gvector_rhombohedron

    subroutine check_the_size_of_rhombohedron

      integer, dimension(3) :: kmax = (/0,0,0/), kmaxp = (/0,0,0/)
      real(kind=DP)         :: length_of_G
      integer i, j, k

      if(ipri >= 1) call wd_n_rGv_and_n_rGpv
      call getttr(rltv,ttr)
      Kloop: do k = - n_rGpv(3), n_rGpv(3)
         Jloop: do j = - n_rGpv(2), n_rGpv(2)
            Iloop: do i = - n_rGpv(1), n_rGpv(1)
               length_of_G = dsqrt(ttr(1)*i*i + ttr(2)*j*j + ttr(3)*k*k &
                    &            +      ttr(4)*i*j + ttr(5)*j*k + ttr(6)*k*i )
               if(length_of_G <= gmax2) then
                  if(kmax(1) < abs(i)) kmax(1) = abs(i)
                  if(kmax(2) < abs(j)) kmax(2) = abs(j)
                  if(kmax(3) < abs(k)) kmax(3) = abs(k)
               else if(length_of_G <= gmaxp) then
                  if(kmaxp(1) < abs(i)) kmaxp(1) = abs(i)
                  if(kmaxp(2) < abs(j)) kmaxp(2) = abs(j)
                  if(kmaxp(3) < abs(k)) kmaxp(3) = abs(k)
               endif
            enddo Iloop
         enddo Jloop
      enddo Kloop

      do i = 1, 3
         kmaxp(i) = max(kmax(i),kmaxp(i))
         kmax(i)  = kmax(i) + 1
         kmaxp(i) = kmaxp(i) + 1
      enddo

    end subroutine check_the_size_of_rhombohedron

    subroutine wd_n_rGv_and_n_rGpv
      if(ipri >= 1) then
         write(nfout,*) ' ! Size of rhombohedrons which contain spheres of G_vectors being smaller than gmax and gmaxp'
         write(nfout,110) n_rGv(1), n_rGv(2), n_rGv(3)
         write(nfout,111) n_rGpv(1), n_rGpv(2), n_rGpv(3)
         write(nfout,112) n_rGpv_reduced(1:3)
      end if
110   FORMAT(' ',' KNX ,KNY ,KNZ  = ',3I6)
111   FORMAT(' ',' KNXP,KNYP,KNZP = ',3I6)
112   format(' ',' knxp_reduced, knyp_reduced, knzp_reduced = ',3i6)
      if(sw_positron /= OFF .and. ipri >= 1) &
           & write(nfout,'("  knx, kny, knz (positron) = ",3i6)') n_rGv_pstrn(1:3)
    end subroutine wd_n_rGv_and_n_rGpv
  end subroutine m_pwBS_assume_G_rhombohedron

  subroutine m_pwBS_generate_G_vectors
    real(kind=DP)         :: gmax2, gmax2_pstrn

! =========== KT_add ================================ 13.0F
    real(kind=DP)         :: gmax2_exx
! =================================================== 13.0F

    real(kind=DP), allocatable, dimension(:,:,:) :: G_length_in_cube
    real(kind=DP), allocatable, dimension(:)     :: gr_t !d(kgp)
    integer               :: id_sname = -1

    call tstatc0_begin('m_pwBS_generate_G_vectors ',id_sname)

    gmax2 = gmax*2
    gmax2_pstrn = gmax_positron*2

! =========== KT_add ================================ 13.0F
    gmax2_exx = gmax_exx*2
! =================================================== 13.0F

    call adjust_n_rGv_to_2l3m5n  ! -(contained here) fft_box_size_WF,CD -> n_rGV,n_rGpv

    allocate(G_length_in_cube(2*n_rGpv(1)+3,2*n_rGpv(2)+3, 2*n_rGpv(3)+3))
    call count_Gvectors_in_spheres

    call G_vectors_alloc           ! ngabc
    allocate(gr_t(kgp))

    call get_Gvectors_in_a_gmaxp_sphere

    deallocate(gr_t)
    deallocate(G_length_in_cube)

    if(ipri>=1) write(nfout,*) ' !kg  = ', kg
    if(ipri>=1) write(nfout,*) ' !kgp = ', kgp
    if(ipri>=1) write(nfout,*) ' !kgp_reduced = ',kgp_reduced

    if(sw_positron /= OFF .and. ipri>=1) write(nfout,*) ' !kg_pwf = ', kg_pwf

! =========== KT_add ================================ 13.0F
    if (sw_hybrid_functional /= OFF .and. use_fft_exx ) then
       if (ipri>=1) write(nfout,*) ' !kg_exx  = ', kg_exx
    endif
    if (sw_hybrid_functional /= OFF ) then
       if (ipri>=1) write(nfout,*) ' !kgp_exx = ', kgp_exx
    endif
! =================================================== 13.0F

    call tstatc0_end(id_sname)
  contains
    subroutine get_Gvectors_in_a_gmaxp_sphere
      integer             :: i, j, k
      real(kind=DP)       :: length_of_Gvector

      call getttr(rltv,ttr)

      kg = 0
      kgp = 0
!!$      if(sw_positron /= OFF) kg_pwf = 0
      do k = -n_rGpv(3)-1, n_rGpv(3)+1
         do j = -n_rGpv(2)-1, n_rGpv(2)+1
            do i = -n_rGpv(1)-1, n_rGpv(1)+1
               length_of_Gvector = G_length_in_cube(i+n_rGpv(1)+2,j+n_rGpv(2)+2, k+n_rGpv(3)+2)
!!$               length_of_Gvector = dsqrt(ttr(1)*i*i+ttr(2)*j*j+ttr(3)*k*k &
!!$                    &                  + ttr(4)*i*j+ttr(5)*j*k+ttr(6)*k*i)
               if(length_of_Gvector <= gmax2) kg  = kg + 1
!!$               if(sw_positron/=OFF .and. length_of_Gvector<=gmax2_pstrn) kg_pwf=kg_pwf+1
               if(length_of_Gvector <= gmaxp) then
                  kgp = kgp + 1
                  ngabc(kgp,1) = i
                  ngabc(kgp,2) = j
                  ngabc(kgp,3) = k
                  gr_t(kgp) = length_of_Gvector
               endif
            enddo
         enddo
      enddo
      if(ipri >= 2) then
         write(nfout,'(" -- gr_t -- <<m_pwBS_generate_G_vectors>>")')
         do i = 1, min(20, kgp)
            write(nfout,'(" gr(",3i4,") = ",d12.4)') ngabc(i,1),ngabc(i,2),ngabc(i,3),gr_t(i)
         end do
      end if

#ifdef _HEAP_SORT_
     ! Sorting G using heap sorting (M.Okamoto)
!!$      call sort_gvec_heap(ttr,kgp,kgp,ngabc)
!!$      call shellsort(nfout,ipri,ttr,kgp,kgp,ngabc)
      !!$call sort_gvec_heap2(ttr,kgp,kgp,ngabc,gr_t)
      call sort_gvec_heap3(ttr,kgp,kgp,ngabc,gr_t)
      if(ipri >= 2) then
         write(nfout,'(" -- gr after (sort_gvec_heap2) -- <<m_pwBS_generate_G_vectors>>")')
         do i = 1, min(20, kgp)
            write(nfout,'(" gr(",3i4,") = ",d12.4)') ngabc(i,1),ngabc(i,2),ngabc(i,3),gr_t(i)
         end do
      end if
      call shellsort2(nfout,ipri,ttr,kgp,kgp,ngabc,gr_t)
#elif _SIMPLE_SORT_
     ! Sorting G simply (M.Okamoto)
      call sort_gvec_simple(ttr,kgp,kgp,ngabc)
#else
      call sort_gvec_heap3(ttr,kgp,kgp,ngabc,gr_t)
      if(ipri >= 2) then
         write(nfout,'(" -- gr after (sort_gvec_heap2) -- <<m_pwBS_generate_G_vectors>>")')
         do i = 1, min(20, kgp)
            write(nfout,'(" gr(",3i4,") = ",d12.4)') ngabc(i,1),ngabc(i,2),ngabc(i,3),gr_t(i)
         end do
      end if
      call shellsort2(nfout,ipri,ttr,kgp,kgp,ngabc,gr_t)
#endif
#ifndef _SIMPLE_SORT_      
      if(ipri >= 2) then
         write(nfout,'(" -- gr after sorted -- <<m_pwBS_generate_G_vectors>>")')
         do i = 1, min(20, kgp)
            write(nfout,'(" gr(",3i4,") = ",d12.4)') ngabc(i,1),ngabc(i,2),ngabc(i,3),gr_t(i)
         end do
      end if
#endif      

    end subroutine get_Gvectors_in_a_gmaxp_sphere

    subroutine count_Gvectors_in_spheres()
      integer i, j, k
      integer, parameter  :: DP = kind(1.d0)
      real(kind=DP)       :: length_of_Gvector
      call getttr(rltv,ttr)
      kg = 0
      kgp = 0
      kgp_reduced = 0
      kg0 = 0
      if(sw_positron /= OFF) kg_pwf = 0

! =========== KT_add ================================ 13.0F
      if (sw_hybrid_functional /= OFF .and. use_fft_exx ) then
         kg_exx = 0
      endif
      if (sw_hybrid_functional /= OFF ) kgp_exx = 0
! =================================================== 13.0F

! ==== KT_add ============= 13.0U2
      if ( sw_modified_TFW_functional /= OFF ) kg_tfw = 0
! ========================= 13.0U2

      do k = -n_rGpv(3)-1, n_rGpv(3)+1
         do j = -n_rGpv(2)-1, n_rGpv(2)+1
            do i = -n_rGpv(1)-1, n_rGpv(1)+1
               length_of_Gvector = dsqrt(ttr(1)*i*i + ttr(2)*j*j &
                    &                  + ttr(3)*k*k + ttr(4)*i*j &
                    &                  + ttr(5)*j*k + ttr(6)*k*i)
               if(length_of_Gvector <= gmax2) kg = kg + 1
               if(length_of_Gvector <= gmax)  kg0 = kg0 + 1
               if(sw_positron/=OFF .and. length_of_Gvector<=gmax2_pstrn) kg_pwf=kg_pwf+1

! =========== KT_add ================================ 13.0F
               if(sw_hybrid_functional/=OFF .and. use_fft_exx ) then
#ifdef FFTW3
                  if ( length_of_Gvector<=gmax2_exx) then
#else
                  if ( length_of_Gvector<=gmax2) then
#endif
                     kg_exx = kg_exx +1
                  endif
               endif
               if(sw_hybrid_functional /= OFF) then
                  if( length_of_Gvector <= gmaxp_exx) then
                     kgp_exx = kgp_exx+1
                  endif
               endif
! =================================================== 13.0F

! ============ KT_add ========== 13.0U2
               if ( sw_modified_TFW_functional /= OFF ) then
                  if ( length_of_Gvector<= gmaxp/2 ) then
                     kg_tfw = kg_tfw +1
                  endif
               endif
! ============================== 13.0U2

               if(length_of_Gvector <= gmaxp) kgp = kgp + 1
               if(length_of_Gvector <= gmaxp_reduced) kgp_reduced = kgp_reduced + 1
               G_length_in_cube(i+n_rGpv(1)+2,j+n_rGpv(2)+2, k+n_rGpv(3)+2) = length_of_Gvector
            enddo
         enddo
      enddo
      if(ipri >= 1) then
         write(nfout,'(" !pwBS kg0, kg, kgp = ",3i10)') kg0, kg, kgp
         write(nfout,'(" !pwBS kgp_reduced  = ",i10)') kgp_reduced
         write(nfout,'(" !pwBS   kg0 = (#G(<=Gmax)), kg = (#G(<=2Gmax)), kgp = (#G(<=Gmaxp))")')
      end if
    end subroutine count_Gvectors_in_spheres

    subroutine adjust_n_rGv_to_2l3m5n
      integer :: i
      do i = 1, 3
         if(n_rGv(i) > fft_box_size_WF(i,1)/2) n_rGv(i)   = fft_box_size_WF(i,1)/2
!!$         if(n_rGpv(i) > fft_box_size_CD(i,1)/2) n_rGpv(i) = fft_box_size_CD(i,1)/2
         if(n_rGpv_reduced(i) > fft_box_size_CD(i,1)/2) n_rGpv_reduced(i) = fft_box_size_CD(i,1)/2
      enddo
    end subroutine adjust_n_rGv_to_2l3m5n
  end subroutine m_pwBS_generate_G_vectors

  subroutine m_pwBS_rebuild_gr_l()
    integer :: i
    call getttr(rltv,ttr)
    gr_l = 0.d0
    do i = ista_kngp, iend_kngp       !for mpi
       gr_l(i) =      dsqrt(ttr(1)*ngabc(i,1)*ngabc(i,1) &
            &             + ttr(2)*ngabc(i,2)*ngabc(i,2) &
            &             + ttr(3)*ngabc(i,3)*ngabc(i,3) &
            &             + ttr(4)*ngabc(i,1)*ngabc(i,2) &
            &             + ttr(5)*ngabc(i,2)*ngabc(i,3) &
            &             + ttr(6)*ngabc(i,3)*ngabc(i,1))
    enddo
    if(ipri >= 2 ) then
       do i = ista_kngp, iend_kngp  !for mpi
          write(nfout,'(" ", i9," = ",f10.5, 3i5)') i, gr_l(i) &
               &, ngabc(i,1),ngabc(i,2),ngabc(i,3)
       enddo
    end if
  end subroutine m_pwBS_rebuild_gr_l

  subroutine m_pwBS_calc_length_of_G
    integer :: id_sname = -1
    integer :: i, k, ngs
    real(kind=DP) :: length_start, length

    call tstatc0_begin('m_pwBS_calc_length_of_G ',id_sname)

    do i = ista_kngp, iend_kngp       !for mpi
       gr_l(i) =      dsqrt(ttr(1)*ngabc(i,1)*ngabc(i,1) &
            &             + ttr(2)*ngabc(i,2)*ngabc(i,2) &
            &             + ttr(3)*ngabc(i,3)*ngabc(i,3) &
            &             + ttr(4)*ngabc(i,1)*ngabc(i,2) &
            &             + ttr(5)*ngabc(i,2)*ngabc(i,3) &
            &             + ttr(6)*ngabc(i,3)*ngabc(i,1))
    enddo

    if(ipri >= 2 ) then
       do i = ista_kngp, iend_kngp  !for mpi
          write(nfout,'(" ", i9," = ",f10.5, 3i5)') i, gr_l(i) &
               &, ngabc(i,1),ngabc(i,2),ngabc(i,3)
       enddo
    end if

    ngshell = 1
    k = ista_kngp+1
    if(iend_kngp < k) goto 11
    length_start = gr_l(ista_kngp)
1   continue
    shellloop: do i = k, iend_kngp
       length = gr_l(i)
       if(length > length_start + DELTA) then
          exit shellloop
       else if(i == iend_kngp) then
          go to 11
       end if
    end do shellloop
    length_start = length
    k = i
    if(k <= iend_kngp) then
       ngshell = ngshell + 1
       goto 1
    end if
11 continue
    if(ipri >= 1 ) write(nfout,'(" !pwbs  ngshell = ",i8)') ngshell

    allocate(ngshell_range(ngshell,2))
    ngs = 1
    k = ista_kngp+1
    length_start = gr_l(ista_kngp)
    ngshell_range(1,1) = ista_kngp
    if(iend_kngp < k) goto 22

2   continue
    shellloop2: do i = k, iend_kngp
       length = gr_l(i)
       if(length > length_start + DELTA) then
          ngshell_range(ngs,2)   = i - 1
          if(ngs+1 <= ngshell) ngshell_range(ngs+1,1) = i
          exit shellloop2
       else if( i == iend_kngp) then
          go to 22
       end if
    end do shellloop2
    length_start = length
    k = i
    if( k <= iend_kngp) then
       ngs = ngs+1
       goto 2
    end if
22 continue
    ngshell_range(ngshell,2) = iend_kngp
    if(ipri >= 2) then
       write(nfout,'(" !pwbs ngshell_range")')
       write(nfout,'(" !pwbs  ngshell, range1 - range2, nelement, gr")')
       do ngs = 1, ngshell
          write(nfout,'(" !pwbs ",i8, i8," - ",i8, i8, f20.8)') ngs, ngshell_range(ngs,1)&
               & ,ngshell_range(ngs,2),ngshell_range(ngs,2)-ngshell_range(ngs,1)+1,gr_l(ngshell_range(ngs,1))
       end do
    end if

    call tstatc0_end(id_sname)
  end subroutine m_pwBS_calc_length_of_G

  subroutine m_pwBS_setup_FFTmapfunctions
    integer :: id, i
    integer :: igf1, igf2, igf3
    integer :: id_sname = -1
    integer, parameter :: NPRINTLINE = 20
    call tstatc0_begin('m_pwBS_setup_FFTmapfunctions ',id_sname)

    id = fft_box_size_WF(1,0)
    do i = 1, kg
       igf1 = ngabc(i,1) + 1
       igf2 = ngabc(i,2) + 1
       igf3 = ngabc(i,3) + 1
       if(ngabc(i,1) <= -1) igf1 = igf1 + fft_box_size_WF(1,1)
       if(ngabc(i,2) <= -1) igf2 = igf2 + fft_box_size_WF(2,1)
       if(ngabc(i,3) <= -1) igf3 = igf3 + fft_box_size_WF(3,1)
       igf(i) = igf1 + (igf2-1)*id + (igf3-1)*id*fft_box_size_WF(2,0)
    enddo

    ! ------- Positron start
    if(sw_positron /= OFF) then
       if(    fft_box_size_WF(1,1) == fft_box_size_pWF(1,1) .and. &
            & fft_box_size_WF(2,1) == fft_box_size_pWF(2,1) .and. &
            & fft_box_size_WF(3,1) == fft_box_size_pWF(3,1) .and. &
            & kg == kg_pwf ) then
          igf_pstrn(1:kg_pwf) = igf(1:kg_pwf)
       else
          id = fft_box_size_pWF(1,0)
          do i = 1, kg_pwf
             igf1 = ngabc(i,1) + 1
             igf2 = ngabc(i,2) + 1
             igf3 = ngabc(i,3) + 1
             if(ngabc(i,1) <= -1) igf1 = igf1 + fft_box_size_pWF(1,1)
             if(ngabc(i,2) <= -1) igf2 = igf2 + fft_box_size_pWF(2,1)
             if(ngabc(i,3) <= -1) igf3 = igf3 + fft_box_size_pWF(3,1)
             igf_pstrn(i) = igf1 + (igf2-1)*id + (igf3-1)*id*fft_box_size_pWF(2,0)
          end do
       end if
    end if
    ! ------- Positron end

! ==============  KT_add ============================= 13.0F
    if ( sw_hybrid_functional /= OFF .and. use_fft_exx ) then
#ifndef FFTW3
       fft_box_size_exx = fft_box_size_WF
#endif
       if(    fft_box_size_WF(1,1) == fft_box_size_exx(1,1) .and. &
            & fft_box_size_WF(2,1) == fft_box_size_exx(2,1) .and. &
            & fft_box_size_WF(3,1) == fft_box_size_exx(3,1) .and. &
            & kg == kg_exx ) then
          igf_exx(1:kg_exx) = igf(1:kg_exx)
       else
          id = fft_box_size_exx(1,0)
          do i = 1, kg_exx
             igf1 = ngabc(i,1) + 1
             igf2 = ngabc(i,2) + 1
             igf3 = ngabc(i,3) + 1
             if(ngabc(i,1) <= -1) igf1 = igf1 + fft_box_size_exx(1,1)
             if(ngabc(i,2) <= -1) igf2 = igf2 + fft_box_size_exx(2,1)
             if(ngabc(i,3) <= -1) igf3 = igf3 + fft_box_size_exx(3,1)
             igf_exx(i) = igf1 + (igf2-1)*id + (igf3-1)*id*fft_box_size_exx(2,0)
          end do
       end if
    endif
! ==================================================== 13.0F

    id = fft_box_size_CD(1,0)
    do i = ista_kngp, iend_kngp  !for mpi
       if(i <= kgp_reduced) then
          igf1 = ngabc(i,1) + 1
          igf2 = ngabc(i,2) + 1
          igf3 = ngabc(i,3) + 1
          if(ngabc(i,1) <= -1) igf1 = igf1 + fft_box_size_CD(1,1)
          if(ngabc(i,2) <= -1) igf2 = igf2 + fft_box_size_CD(2,1)
          if(ngabc(i,3) <= -1) igf3 = igf3 + fft_box_size_CD(3,1)
          igfp_l(i) = igf1 + (igf2-1)*id + (igf3-1)*id*fft_box_size_CD(2,0)
#ifdef _MPIFFT_
          igfp_l_c(i) = igf1 + fft_box_size_CD_c(1,0)*( &
               &                    (igf2-1)+ (igf3-1)*fft_box_size_CD_c(2,0))
          igfp_nonpara(i) = igf1 + fft_box_size_CD_nonpara(1,0)*( &
               &                    (igf2-1)+ (igf3-1)*fft_box_size_CD_nonpara(2,0))
#endif
       else
          igfp_l(i) = -1
#ifdef _MPIFFT_
          igfp_l_c(i) = -1
          igfp_nonpara(i) = -1
#endif
       end if
    enddo

    if(iprifftmap >= 2) then
       if(mype == 0) then
          if(sw_positron == OFF) then
             write(nfout,*) ' -- ngabc(1-3),igf,igfp_l --'
             id = NPRINTLINE
             if(id > iend_kngp) id = iend_kngp
             do i = 1, id
                write(nfout,'(" ",8i8)') &
                     & i,ngabc(i,1),ngabc(i,2),ngabc(i,3),igf(i),igfp_l(i)
             end do
          else
             write(nfout,*) ' -- ngabc(1-3),igf,igfp_l,igf_pstrn --'
             id = NPRINTLINE
             if(id > iend_kngp) id = iend_kngp
             do i = 1, id
                write(nfout,'(" ",8i8)') &
                     & i,ngabc(i,1),ngabc(i,2),ngabc(i,3),igf(i),igfp_l(i),igf_pstrn(i)
             end do
          end if
#ifdef _MPIFFT_
          write(nfout,*) ' -- ngabc(1-3),igf,igfp_nonpara --'
          id = NPRINTLINE
          if(id > iend_kngp) id = iend_kngp
          do i = 1, id
             write(nfout,'(" ",8i8)') &
                  & i,ngabc(i,1),ngabc(i,2),ngabc(i,3),igf(i),igfp_nonpara(i)
          end do
#endif
       endif
       if(mype == npes-1) then
          id = kg-(NPRINTLINE-1)
          if(id < ista_kngp) id = ista_kngp
          if(sw_positron == OFF) then
             do i = id,kg
                write(nfout,'(" ",8i8)') &
                     & i,ngabc(i,1),ngabc(i,2),ngabc(i,3),igf(i),igfp_l(i)
             end do
          else
             if(kg > kg_pwf) kg = kg_pwf
             do i = id,kg
                write(nfout,'(" ",8i8)') &
                     & i,ngabc(i,1),ngabc(i,2),ngabc(i,3),igf(i),igfp_l(i),igf_pstrn(i)
             end do
          end if
       endif
    end if

    if ( sw_hybrid_functional /= OFF ) then
       if(    fft_box_size_CD(1,1) == fft_box_size_CD_exx(1,1) .and. &
            & fft_box_size_CD(2,1) == fft_box_size_CD_exx(2,1) .and. &
            & fft_box_size_CD(3,1) == fft_box_size_CD_exx(3,1) .and. &
            & kgp == kgp_exx ) then
          igfp_exx(ista_kngp:iend_kngp) = igfp_l(ista_kngp:iend_kngp)
          call mpi_allreduce(MPI_IN_PLACE,igfp_exx,kgp_exx,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
       else
          id = fft_box_size_CD_exx(1,0)
          do i = 1, kgp_exx
             igf1 = ngabc(i,1) + 1
             igf2 = ngabc(i,2) + 1
             igf3 = ngabc(i,3) + 1
             if(ngabc(i,1) <= -1) igf1 = igf1 + fft_box_size_CD_exx(1,1)
             if(ngabc(i,2) <= -1) igf2 = igf2 + fft_box_size_CD_exx(2,1)
             if(ngabc(i,3) <= -1) igf3 = igf3 + fft_box_size_CD_exx(3,1)
             igfp_exx(i) = igf1 + (igf2-1)*id + (igf3-1)*id*fft_box_size_CD_exx(2,0)
          end do
       end if
    endif
    call tstatc0_end(id_sname)
  end subroutine m_pwBS_setup_FFTmapfunctions

  subroutine m_pwBS_GminusGmapfunction()
    integer :: i,igf1,igf2,igf3
    igpo = 0
    nmatsz2 = 0
    do i = 1, kg
       igf1 = ngabc(i,1); igf2 = ngabc(i,2); igf3 = ngabc(i,3)
       if(          igf1 >= -n_rGsv(1) .and. igf1 <= n_rGsv(1) &
            & .and. igf2 >= -n_rGsv(2) .and. igf2 <= n_rGsv(2) &
            & .and. igf3 >= -n_rGsv(3) .and. igf3 <= n_rGsv(3)) then
          igpo(igf1,igf2,igf3) = i
          nmatsz2 = i
       end if
    end do
    if(ipri >= 1) write(nfout,'(" !pwbs nmatsz2 = ",i8 &
         &                ," <<m_pwBS_GminusGmapfunction>>")') nmatsz2
  end subroutine m_pwBS_GminusGmapfunction

  subroutine decide_g_list_size(nfout, mmdim1, mmdim2, mmdim3)
    integer, intent(in)   :: nfout
    integer, intent(out)  :: mmdim1, mmdim2, mmdim3
    mmdim1 = 2*n_rGpv(1)
    mmdim2 = 2*n_rGpv(2)
    mmdim3 = 2*n_rGpv(3)
    if(ipri>=1) then
       write(nfout,*) ' n_rGpv = ', n_rGpv
       write(nfout,*) ' mmdim  = ', mmdim1, mmdim2, mmdim3
    end if
  end subroutine decide_g_list_size

  subroutine m_pwBS_G_trans_functions
!              -> ngpt_l
    integer, pointer ,dimension(:,:,:) ::g_list
    integer :: mmdim1, mmdim2, mmdim3
    integer error_count

    if(ipri >= 2) then
       write(nfout,'(" !pwBS: nopr+af = ",i6," <<m_pwBS_G_trans_functions>>")') nopr+af
       write(nfout,'(" !pwBS: nopr, af = ",2i6,"<<m_pwBS_G_trans_functions>>")') nopr,af
    end if
    if(nopr+af > 100) call phase_error_with_msg(nfout,' !! illegal number of nopr+af',__LINE__,__FILE__)
    allocate(op_br(3,3,nopr+af))
    call m_CS_op_in_PUCV(nfout,op_br,nopr+af)
    call decide_g_list_size(nfout, mmdim1, mmdim2,mmdim3 )
    if(ipri >= 1) write(nfout,'(" !pwBS: g_list size = ",3i9)') mmdim1, mmdim2, mmdim3
    allocate(g_list(0:mmdim1,0:mmdim2,0:mmdim3))

    call substitute_g_list
!TS** conjecture of the symmetry point by using the list vector g_list...
    call symmetric_G_points_using_g_list
    call counting_error_points(error_count)
    if(error_count > 0) then
       if(ipri>=1) write(nfout &
            & ,'(" !pwBS: (G_vector_transformation_functions) error_count = ",i9)') error_count
       call symmetric_G_points
    else
       if(ipri>=1) write(nfout,'(" !pwBS: (ngpt_l)s are all decided by using g_list")')
    endif

    deallocate(g_list)
    deallocate(op_br)
  contains

    subroutine substitute_g_list
      integer i, int_gx, int_gy, int_gz
      real(kind=DP) :: gxL, gyL, gzL

      g_list = 0

      if(2*n_rGpv(1)-mmdim1 /=0 .and. 2*n_rGpv(2)-mmdim2 /= 0 .and. 2*n_rGpv(3)-mmdim3 /=0) then
         gxL = 2.d0*n_rGpv(1); gyL = 2.d0*n_rGpv(2); gzL = 2.d0*n_rGpv(3)
!** (This do-loop can be vectorizable.)
!xocl spread do/ind_kngp
         do  i = 1, kgp
            int_gx = (dble(ngabc(i,1))+n_rGpv(1))/gxL*dble(mmdim1)
            int_gy = (dble(ngabc(i,2))+n_rGpv(2))/gyL*dble(mmdim2)
            int_gz = (dble(ngabc(i,3))+n_rGpv(3))/gzL*dble(mmdim3)
            int_gx = max(int_gx,0)
            int_gy = max(int_gy,0)
            int_gz = max(int_gz,0)
            g_list(int_gx,int_gy,int_gz) = i
         enddo
!xocl end spread max(g_list)
      else
!xocl spread do/ind_kngp
         do  i = 1, kgp
            int_gx = ngabc(i,1)+n_rGpv(1)
            int_gy = ngabc(i,2)+n_rGpv(2)
            int_gz = ngabc(i,3)+n_rGpv(3)
            int_gx = max(int_gx,0)
            int_gy = max(int_gy,0)
            int_gz = max(int_gz,0)
            g_list(int_gx,int_gy,int_gz) = i
         enddo
!xocl end spread max(g_list)
      endif
    end subroutine substitute_g_list

    subroutine symmetric_G_points_using_g_list
      integer i
      integer no
      real(kind=DP) :: fx,fy,fz, px,py,pz, distance
      integer       :: int_fx, int_fy, int_fz, nng
      real(kind=DP),  parameter :: epsilon = 1.0d-5
      real(kind=DP) :: gxL, gyL, gzL

      if(2*n_rGpv(1)-mmdim1 /=0 .and. 2*n_rGpv(2)-mmdim2 /= 0 .and. 2*n_rGpv(3)-mmdim3 /=0) then
         gxL = 2.d0*n_rGpv(1); gyL = 2.d0*n_rGpv(2); gzL = 2.d0*n_rGpv(3)
         do i = ista_kngp, iend_kngp  !for mpi
            do no = 1, nopr + af

               px = ngabc(i,1)
               py = ngabc(i,2)
               pz = ngabc(i,3)

               fx=op_br(1,1,no)*px +op_br(1,2,no)*py +op_br(1,3,no)*pz
               fy=op_br(2,1,no)*px +op_br(2,2,no)*py +op_br(2,3,no)*pz
               fz=op_br(3,1,no)*px +op_br(3,2,no)*py +op_br(3,3,no)*pz

               int_fx = (fx + n_rGpv(1)+epsilon)/gxL*dble(mmdim1)
               int_fy = (fy + n_rGpv(2)+epsilon)/gyL*dble(mmdim2)
               int_fz = (fz + n_rGpv(3)+epsilon)/gzL*dble(mmdim3)

               int_fx = max(int_fx,0)
               int_fy = max(int_fy,0)
               int_fz = max(int_fz,0)

               int_fx = min(int_fx,mmdim1)
               int_fy = min(int_fy,mmdim2)
               int_fz = min(int_fz,mmdim3)

               nng = g_list(int_fx,int_fy,int_fz)

               if(nng /= 0) then
                  distance =  abs(ngabc(nng,1)-fx) + abs(ngabc(nng,2)-fy) &
                       &    + abs(ngabc(nng,3)-fz)
                  if(distance < epsilon) then
                     ngpt_l(i,no) = nng
                  else
                     ngpt_l(i,no) = -1
                  end if
               else
                  ngpt_l(i,no) = -1
               end if
            enddo
         enddo
      else

         do i = ista_kngp, iend_kngp  !for mpi
            do no = 1, nopr + af

               px = ngabc(i,1)
               py = ngabc(i,2)
               pz = ngabc(i,3)

               fx=op_br(1,1,no)*px +op_br(1,2,no)*py +op_br(1,3,no)*pz
               fy=op_br(2,1,no)*px +op_br(2,2,no)*py +op_br(2,3,no)*pz
               fz=op_br(3,1,no)*px +op_br(3,2,no)*py +op_br(3,3,no)*pz

               int_fx = fx + n_rGpv(1) + epsilon
               int_fy = fy + n_rGpv(2) + epsilon
               int_fz = fz + n_rGpv(3) + epsilon

               int_fx = max(int_fx,0)
               int_fy = max(int_fy,0)
               int_fz = max(int_fz,0)

               int_fx = min(int_fx,mmdim1)
               int_fy = min(int_fy,mmdim2)
               int_fz = min(int_fz,mmdim3)

               nng = g_list(int_fx,int_fy,int_fz)

               if(nng /= 0) then
                  distance =  abs(ngabc(nng,1)-fx) + abs(ngabc(nng,2)-fy) &
                       &    + abs(ngabc(nng,3)-fz)
                  if(distance < epsilon) then
                     ngpt_l(i,no) = nng
                  else
                     ngpt_l(i,no) = -1
                  end if
               else
                  ngpt_l(i,no) = -1
               end if
            enddo
         enddo
      endif
    end subroutine symmetric_G_points_using_g_list

    subroutine counting_error_points(error_count)
      integer, intent(out) :: error_count

      integer :: i, no
      integer :: e_mpi  ! mpi
      error_count = 0

      do i = ista_kngp, iend_kngp  !for mpi
         do no=1,nopr + af
            if(ngpt_l(i,no) == -1) error_count = error_count + 1
         enddo
      enddo

      if(npes > 1 ) then
         call mpi_allreduce(error_count,e_mpi,1,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
         error_count = e_mpi
      end if

    end subroutine counting_error_points

    subroutine symmetric_G_points
      integer :: ierr, i, no, loop, jstart, jend, jstep, j
      integer, parameter :: Nequiv = 50
      integer, parameter :: MAXERR = 100
      real(kind=DP) :: px,py,pz, fx,fy,fz, ff1
      ierr = 0

      do i = ista_kngp, iend_kngp  !for mpi
         px = ngabc(i,1)
         py = ngabc(i,2)
         pz = ngabc(i,3)
         do no=1,nopr + af

            if(ngpt_l(i,no).eq.-1) then

               fx=op_br(1,1,no)*px +op_br(1,2,no)*py +op_br(1,3,no)*pz
               fy=op_br(2,1,no)*px +op_br(2,2,no)*py +op_br(2,3,no)*pz
               fz=op_br(3,1,no)*px +op_br(3,2,no)*py +op_br(3,3,no)*pz
               do loop = 1, 2
                  if(loop.eq.1) then
                     jstart = i - Nequiv + 1
                     if(jstart.le.0) jstart = 1
                     jend = kgp
                     jstep = 1
                  else if(loop.eq.2) then
                     jstart = i - Nequiv
                     jend = 1
                     jstep = -1
                  endif
!*vocl loop, vi(j)
                  do j = jstart, jend, jstep
                     ff1 = abs(ngabc(j,1)-fx) +abs(ngabc(j,2)-fy) +abs(ngabc(j,3)-fz)
                     if (ff1.le.1.d-5) then
                        ngpt_l(i,no) = j
                        goto 111
                     end if
                  enddo
               enddo
               if(ipri >= 1) write(nfout,130) i,no
               ierr = ierr + 1
               if(ierr.ge.MAXERR) goto 1001
111            continue

            endif

         enddo
      enddo

130   FORMAT(' ','THERE IS NO PAIR FOR NG,NOP=',2I8)
1001  if(ierr.ge.1) then
         if(ipri >= 1 ) write(nfout,*) ' No pair was discovered for at least ',ierr &
              &       ,'s combinations'
         call phase_error_with_msg(nfout,' No pair was discovered',__LINE__,__FILE__)
      endif

    end subroutine symmetric_G_points

  end subroutine m_pwBS_G_trans_functions

  subroutine m_pwBS_sphrp_l()
    integer :: is
    real(kind=DP), allocatable, dimension(:) :: ylm
    allocate(ylm(ista_kngp:iend_kngp))
    if(ipri >= 1) write(nfout,'(" !pwBS: nel_Ylm = ",i5," (m_pwBS_sphrp)")') nel_Ylm
    if(nel_Ylm >= 1 .and. .not.allocated(ylm_l)) call phase_error_with_msg(nfout,' ylm_l is not allocated (m_pwBS_sphrp_l)'&
                                                                          ,__LINE__,__FILE__)
    do is = 1, nel_Ylm
       call m_pwBS_sphrp2(is,rltv,ista_kngp,iend_kngp,ylm)
       ylm_l(:,is) = ylm(:)
    end do
    deallocate(ylm)
  end subroutine m_pwBS_sphrp_l

  subroutine m_pwBS_sphrp2(is,rltv,ista_ylm,iend_ylm,ylm)
    integer, intent(in) :: is
    real(kind=DP), intent(in), dimension(3,3)  :: rltv
    integer, intent(in) :: ista_ylm, iend_ylm
    real(kind=DP), intent(inout), dimension(ista_ylm:iend_ylm) :: ylm
!!$C spherical harmonics.
!!$c     Some lines are rewritten by T.Yamasaki according to Y.Morikawa's
!!$c   e-mail. at 4th May 1994

    integer :: i,ni
    real(kind=DP) a,b,c,d,e,f,gx,gy,gz
    real(kind=DP) b1x,b1y,b1z,b2x,b2y,b2z,b3x,b3y,b3z

    b1x = rltv(1,1);  b1y = rltv(2,1);  b1z = rltv(3,1)
    b2x = rltv(1,2);  b2y = rltv(2,2);  b2z = rltv(3,2)
    b3x = rltv(1,3);  b3y = rltv(2,3);  b3z = rltv(3,3)

    ni = ista_ylm
    if(ngabc(1,1) == 0 .and. ngabc(1,2) == 0 .and. ngabc(1,3) == 0) then
       if(ni == 1) then
          ni = 2
          ylm(1) = dsqrt(1.d0/PAI4)
       endif
    end if
    if(is == 1) then
       a = dsqrt(1.d0/PAI4)
       do i = ni, iend_ylm  !for mpi
          ylm(i) = a
       end do
    else if(is == 2) then
       a = dsqrt(3.d0/PAI4)
       do i = ni, iend_ylm  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          ylm(i) = a*gx/gr_l(i)
       end do
    else if(is == 3) then
       a = dsqrt(3.d0/PAI4)
       do i = ni, iend_ylm  !for mpi
          gy =   b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          ylm(i) = a*gy/gr_l(i)
       end do
    else if(is == 4) then
       a = dsqrt(3.d0/PAI4)
       do i = ni, iend_ylm  !for mpi
          gz =   b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          ylm(i) = a*gz/gr_l(i)
       end do
    else if(is == 5) then
       a = dsqrt(5.d0/(16*PAI))
       do i = ni, iend_ylm  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2)  + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2)  + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2)  + b3z*ngabc(i,3)
          d = gz**2
          e = gx**2+gy**2+d
          ylm(i) = a*(3*d-e)/e
       end do
    else if(is == 6) then
       a = dsqrt(15.d0/(16*PAI))
       do i = ni, iend_ylm  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2)  + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2)  + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2)  + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          e = b+c+gz**2
          ylm(i) = a*(b-c)/e
       end do
    else if(is == 7) then
       a = dsqrt(15.d0/PAI4)
       do i = ni, iend_ylm  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          e = gr_l(i)**2
          ylm(i) = a*gx*gy/e
       end do
    else if(is == 8) then
       a = dsqrt(15.d0/PAI4)
       do i = ni, iend_ylm  !for mpi
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2)  + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2)  + b3z*ngabc(i,3)
          e = gr_l(i)**2
          ylm(i) = a*gy*gz/e
       end do
    else if(is == 9) then
       a = dsqrt(15.d0/PAI4)
       do i = ni, iend_ylm  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          e = gr_l(i)**2
          ylm(i) = a*gz*gx/e
       end do
    else if(is == 10) then
       a = dsqrt(7.d0/(16*PAI))
       do i = ni, iend_ylm  !for mpi
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          d = gz**2
          e = gr_l(i)**2
          f = e * gr_l(i)
          ylm(i) = a*gz*(5*d-3*e)/f
       end do
    else if(is == 11) then
       a = dsqrt(21.d0/(32*PAI))
       do i = ni, iend_ylm  !for mpi
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2)  + b3z*ngabc(i,3)
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2)  + b3x*ngabc(i,3)
          d = gz**2
          e = gr_l(i)**2
          f = e * gr_l(i)
          ylm(i) = a*gx*(5*d-e)/f
       end do
    else if(is == 12) then
       a = dsqrt(21.d0/(32*PAI))
       do i = ni, iend_ylm  !for mpi
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          d = gz**2
          e = gr_l(i)**2
          f = e * gr_l(i)
          ylm(i) = a*gy*(5*d-e)/f
       end do
    else if(is == 13) then
       a = dsqrt(105.d0/(16*PAI))
       do i = ni, iend_ylm  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          e = gr_l(i)**2
          f = e * gr_l(i)
          ylm(i) = a*gz*(b-c)/f
       end do
    else if(is == 14) then
       a = dsqrt(105.d0/PAI4)
       do i = ni, iend_ylm  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2)  + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2)  + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2)  + b3z*ngabc(i,3)
          f = gr_l(i)**3
          ylm(i) = a*gx*gy*gz/f
       end do
    else if(is == 15) then
       a = dsqrt(35.d0/(32*PAI))
       do i = ni, iend_ylm  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          f = gr_l(i)**3
          ylm(i) = a*gx*(b-3*c)/f
       end do
    else if(is == 16) then
       a = dsqrt(35.d0/(32*pai))
       do i = ni, iend_ylm !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          f = gr_l(i)**3
          ylm(i) = a*gy*(3*b-c)/f
       end do
    else if(is == 17) then
       a = 3.d0/8.d0/dsqrt(PAI4)
       do i = ni, iend_ylm  !for mpi
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          d = gz**2
          e = gr_l(i)**2
          f = e**2
          ylm(i) = a*(5*d*(7*d-6*e)/f+3.d0)
       end do
    else if(is == 18) then
       a = 15.d0/4.d0/dsqrt(10*PAI)
       do i = ni, iend_ylm  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          d = gz**2
          e = gr_l(i)**2
          f = e**2
          ylm(i) = a*gz*gx*(7*d-3*e)/f
       end do
    else if(is == 19) then
       a = 15.d0/4.d0/dsqrt(10.d0*PAI)
       do i = ni, iend_ylm  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2)  + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2)  + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2)  + b3z*ngabc(i,3)
          d = gz**2
          e = gr_l(i)**2
          f = e**2
          ylm(i) = a*gy*gz*(7*d-3*e)/f
       end do
    else if(is == 20) then
       a = 15.d0/8.d0/dsqrt(5*PAI)
       do i = ni, iend_ylm  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2)  + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2)  + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2)  + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b+c+d
          f = e**2
          ylm(i) = a*(7*d-e)*(b-c)/f
       end do
    else if(is == 21) then
       a = 15.d0/4.d0/dsqrt(5*PAI)
       do i = ni, iend_ylm  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          d = gz**2
          e = gr_l(i)**2
          f = e**2
          ylm(i) = a*(7*d-e)*gx*gy/f
       end do
    else if(is == 22) then
       a = 105.d0/4.d0/dsqrt(70*pai)
       do i = ni, iend_ylm  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          f = (b+c+gz**2)**2
          ylm(i) = a*(b-3*c)*gz*gx/f
       end do
    else if(is == 23) then
       a = 105.d0/4.d0/dsqrt(70*PAI)
       do i = ni, iend_ylm  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          f = (b+c+gz**2)**2
          ylm(i) = a*(3.d0*b-c)*gy*gz/f
       end do
    else if(is == 24) then
       a = 105.d0/16.d0/dsqrt(35*PAI)
       do i = ni, iend_ylm  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2)  + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2)  + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2)  + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          f = (b+c+gz**2)**2
          ylm(i) = a*((b-c)**2-4*b*c)/f
       end do
    else if(is == 25) then
       a = 105.d0/4.d0/dsqrt(35*PAI)
       do i = ni, iend_ylm  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2)  + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2)  + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2)  + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          f = (b+c+gz**2)**2
          ylm(i) = a*(b-c)*gx*gy/f
       end do
    end if
  end subroutine m_pwBS_sphrp2

  subroutine m_pwBS_sphrp2_diff(is,rltv,ylmd)
!!$C derivative of spherical harmonics coded by H. Sawada 8th May 1997

    integer, intent(in) :: is
    real(kind=DP), intent(in),  dimension(3,3)                   :: rltv
    real(kind=DP), intent(out), dimension(ista_kngp:iend_kngp,3) :: ylmd

    integer       :: i,ni
    real(kind=DP) :: a,b,c,d,e,e2,e3,gr,gr3,gr5,gx,gy,gz,t
    real(kind=DP) :: b1x,b1y,b1z,b2x,b2y,b2z,b3x,b3y,b3z

    b1x = rltv(1,1);  b1y = rltv(2,1);  b1z = rltv(3,1)
    b2x = rltv(1,2);  b2y = rltv(2,2);  b2z = rltv(3,2)
    b3x = rltv(1,3);  b3y = rltv(2,3);  b3z = rltv(3,3)

    ni = ista_kngp
    if(ngabc(1,1) == 0 .and. ngabc(1,2) == 0 .and. ngabc(1,3) == 0) then
       if(ista_kngp ==1) then
          ni = 2
          ylmd(1,1) = 0.d0
          ylmd(1,2) = 0.d0
          ylmd(1,3) = 0.d0
       endif
    end if
    if(is == 1) then
       a = 0.d0
       do i = ni, iend_kngp !for mpi
          ylmd(i,1) = a
          ylmd(i,2) = a
          ylmd(i,3) = a
       end do
    else if(is == 2) then
       a = dsqrt(3.d0/PAI4)
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          gr = dsqrt( e )
          gr3 = gr * e
          ylmd(i,1) = a* ( -b/gr3 + 1.d0/gr )
          ylmd(i,2) = a* ( -gx*gy/gr3 )
          ylmd(i,3) = a* ( -gx*gz/gr3 )
       enddo
    else if(is == 3) then
       a = dsqrt(3.d0/PAI4)
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          gr = dsqrt( e )
          gr3 = gr * e
          ylmd(i,1) = a* ( -gx*gy/gr3 )
          ylmd(i,2) = a* ( -c/gr3 + 1.d0/gr )
          ylmd(i,3) = a* ( -gy*gz/gr3 )
       enddo
    else if(is == 4) then
       a = dsqrt(3.d0/PAI4)
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          gr = sqrt( e )
          gr3 = gr * e
          ylmd(i,1) = a* ( -gx*gz/gr3 )
          ylmd(i,2) = a* ( -gy*gz/gr3 )
          ylmd(i,3) = a* ( -d/gr3 + 1.d0/gr )
       enddo
    else if(is == 5) then
       a = dsqrt(5.d0/(16*PAI))
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          e2 = e**2
          t  = 2*(3*d - e)/e2
          ylmd(i,1) = a* ( -2*gx/e - t*gx)
          ylmd(i,2) = a* ( -2*gy/e - t*gy)
          ylmd(i,3) = a* (  4*gz/e - t*gz)
       enddo
    else if(is == 6) then
       a = dsqrt(15.d0/(16*PAI))
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          e2 = e**2
          t = -2*(b-c)/e2
          ylmd(i,1) = a* ( t*gx + 2*gx/e )
          ylmd(i,2) = a* ( t*gy - 2*gy/e )
          ylmd(i,3) = a* ( t*gz )
       enddo
    else if(is == 7) then
       a = dsqrt(15.d0/PAI4)
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          e2 = e**2
          ylmd(i,1) = a* ( -2*b*gy/e2 + gy/e )
          ylmd(i,2) = a* ( -2*c*gx/e2 + gx/e )
          ylmd(i,3) = a* ( -2*gx*gy*gz/e2 )
       enddo
    else if(is == 8) then
       a = dsqrt(15.d0/PAI4)
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          e2 = e**2
          ylmd(i,1) = a* ( -2*gx*gy*gz/e2 )
          ylmd(i,2) = a* ( -2*c*gz/e2 + gz/e )
          ylmd(i,3) = a* ( -2*d*gy/e2 + gy/e )
       enddo
    else if(is == 9) then
       a = dsqrt(15.d0/PAI4)
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          e2 = e**2
          ylmd(i,1) = a* ( -2*b*gz/e2 + gz/e )
          ylmd(i,2) = a* ( -2*gx*gy*gz/e2 )
          ylmd(i,3) = a* ( -2*d*gx/e2 + gx/e )
       enddo
    else if(is == 10) then
       a = dsqrt(7.d0/(16*PAI))
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          gr = dsqrt( e )
          gr3 = gr * e
          gr5 = gr3 * e
          t   = 3*(5*d-3*e)/gr5
          ylmd(i,1) = a* ( -6*gx*gz/gr3 - gx*gz*t )
          ylmd(i,2) = a* ( -6*gy*gz/gr3 - gy*gz*t )
          ylmd(i,3) = a* ( 4*d/gr3 - d*t + (5*d-3*e)/gr3 )
       enddo
    else if(is == 11) then
       a = dsqrt(21.d0/(32*PAI))
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          gr = dsqrt( e )
          gr3 = gr * e
          gr5 = gr3 * e
          t   = 3*(5*d-e)/gr5
          ylmd(i,1) = a* ( -2*b/gr3      - b*t     + (5*d-e)/gr3 )
          ylmd(i,2) = a* ( -2*gx*gy/gr3  - gx*gy*t )
          ylmd(i,3) = a* ( 8*gx*gz/gr3   - gx*gz*t )
       enddo
    else if(is == 12) then
       a = dsqrt(21.d0/(32*PAI))
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          gr = dsqrt( e )
          gr3 = gr * e
          gr5 = gr3 * e
          t   = 3*(5*d-e)/gr5
          ylmd(i,1) = a* ( -2*gx*gy/gr3 - gx*gy*t )
          ylmd(i,2) = a* ( -2*c/gr3     - c*t     + (5*d-e)/gr3 )
          ylmd(i,3) = a* ( 8*gy*gz/gr3  - gy*gz*t )
       enddo
    else if(is == 13) then
       a = dsqrt(105.d0/(16*PAI))
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          gr = dsqrt( e )
          gr3 = gr * e
          gr5 = gr3 * e
          t = 3*(b-c)/gr5
          ylmd(i,1) = a* ( -gx*gz*t + 2*gx*gz/gr3 )
          ylmd(i,2) = a* ( -gy*gz*t - 2*gy*gz/gr3 )
          ylmd(i,3) = a* ( -d*t     + (b-c)/gr3 )
       enddo
    else if(is == 14) then
       a = dsqrt(105.d0/PAI4)
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          gr = dsqrt( e )
          gr3 = gr * e
          gr5 = gr3 * e
          ylmd(i,1) = a* ( -3*b*gy*gz/gr5   + gy*gz/gr3 )
          ylmd(i,2) = a* ( -3*gx*c*gz/gr5   + gx*gz/gr3 )
          ylmd(i,3) = a* ( -3*gx*gy*d/gr5   + gx*gy/gr3 )
       enddo
    else if(is == 15) then
       a = dsqrt(35.d0/(32*PAI))
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          gr = dsqrt( e )
          gr3 = gr * e
          gr5 = gr3 * e
          t   = 3*(b-3*c)/gr5
          ylmd(i,1) = a* ( -b*t     + 3*(b-c)/gr3 )
          ylmd(i,2) = a* ( -gx*gy*t - 6*gx*gy/gr3 )
          ylmd(i,3) = a* ( -gx*gz*t )
       enddo
    else if(is == 16) then
       a = dsqrt(35.d0/(32*PAI))
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          gr = dsqrt( e )
          gr3 = gr * e
          gr5 = gr3 * e
          t   = 3*(3*b-c)/gr5
          ylmd(i,1) = a* ( -t*gx*gy  + 6*gx*gy/gr3 )
          ylmd(i,2) = a* ( -t*c      + 3*(b-c)/gr3 )
          ylmd(i,3) = a* ( -t*gy*gz )
       enddo
    else if(is == 17) then
       a = 3.d0/8.d0/dsqrt(PAI4)
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          e2 = e**2
          e3 = e * e2
          t  = 20*d*(7*d - 6*e)/e3
          ylmd(i,1) = a* ( -60*gx*d/e2 - t*gx )
          ylmd(i,2) = a* ( -60*gy*d/e2 - t*gy )
          ylmd(i,3) = a* ( 10*gz*d/e2  - t*gz + 10*gz*(7*d-6*e)/e2 )
       enddo
    else if(is == 18) then
       a = 15.d0/4.d0/dsqrt(10*PAI)
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          e2 = e**2
          e3 = e * e2
          t  = 4*(7*d - 3*e)/e3
          ylmd(i,1) = a* ( -6*b*gz/e2     - t*b*gz + gz*(7*d-3*e)/e2 )
          ylmd(i,2) = a* ( -6*gx*gy*gz/e2 - t*gx*gy*gz )
          ylmd(i,3) = a* ( 8*gx*d/e2      - t*gx*d + gx*(7*d-3*e)/e2 )
       enddo
    else if(is == 19) then
       a = 15.d0/4.d0/dsqrt(10*PAI)
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          e2 = e**2
          e3 = e * e2
          t  = 4*(7*d-3*e)/e3
          ylmd(i,1) = a* ( -6*gx*gy*gz/e2 - t*gx*gy*gz )
          ylmd(i,2) = a* ( -6*c*gz/e2     - t*c*gz     + gz*(7*d-3*e)/e2 )
          ylmd(i,3) = a* ( 8*gy*d/e2      - t*gy*d     + gy*(7*d-3*e)/e2 )
       enddo
    else if(is == 20) then
       a = 15.d0/8.d0/dsqrt(5*PAI)
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          e2 = e**2
          e3 = e * e2
          t  = 4*(b-c)*(7*d-e)/e3
          ylmd(i,1) = a* ( -2*gx*(b-c)/e2 - t*gx + 2*gx*(7*d-e)/e2 )
          ylmd(i,2) = a* ( -2*gy*(b-c)/e2 - t*gy - 2*gy*(7*d-e)/e2 )
          ylmd(i,3) = a* ( 12*gz*(b-c)/e2 - t*gz )
       enddo
    else if(is == 21) then
       a = 15.d0/4.d0/dsqrt(5*PAI)
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          e2 = e**2
          e3 = e * e2
          t  = 4*(7*d-e)/e3
          ylmd(i,1) = a* ( -2*b*gy/e2 - t*b*gy + gy*(7*d-e)/e2 )
          ylmd(i,2) = a* ( -2*gx*c/e2 - t*gx*c + gx*(7*d-e)/e2 )
          ylmd(i,3) = a* ( 12*gx*gy*gz/e2 - t*gx*gy*gz )
       enddo
    else if(is == 22) then
       a = 105.d0/4.d0/dsqrt(70*PAI)
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          e2 = e**2
          e3 = e * e2
          t  = 4*(b-3*c)/e3
          ylmd(i,1) = a* ( -t*b*gz     + 3*gz*(b-c)/e2 )
          ylmd(i,2) = a* ( -t*gx*gy*gz - 6*gx*gy*gz/e2 )
          ylmd(i,3) = a* ( -t*gx*d     + gx*(b-3*c)/e2 )
       enddo
    else if(is == 23) then
       a = 105.d0/4.d0/dsqrt(70*PAI)
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          e2 = e**2
          e3 = e * e2
          t  = 4*(3*b-c)/e3
          ylmd(i,1) = a* ( -t*gx*gy*gz + 6*gx*gy*gz/e2 )
          ylmd(i,2) = a* ( -t*c*gz     + 3*gz*(b-c)/e2 )
          ylmd(i,3) = a* ( -t*gy*d     + gy*(3*b-c)/e2 )
       enddo
    else if(is == 24) then
       a = 105.d0/16.d0/dsqrt(35*PAI)
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          e2 = e**2
          e3 = e * e2
          t = 4*(-4*b*c+(b-c)**2)/e3
          ylmd(i,1) = a* ( -t*gx + (-8*gx*c+4*gx*(b-c))/e2 )
          ylmd(i,2) = a* ( -t*gy + (-8*b*gy-4*gy*(b-c))/e2 )
          ylmd(i,3) = a* ( -t*gz )
       enddo
    else if(is == 25) then
       a = 105.d0/4.d0/dsqrt(35*PAI)
       do i = ni, iend_kngp  !for mpi
          gx = b1x*ngabc(i,1) + b2x*ngabc(i,2) + b3x*ngabc(i,3)
          gy = b1y*ngabc(i,1) + b2y*ngabc(i,2) + b3y*ngabc(i,3)
          gz = b1z*ngabc(i,1) + b2z*ngabc(i,2) + b3z*ngabc(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b + c + d
          e2 = e**2
          e3 = e * e2
          t = 4*(b-c)/e3
          ylmd(i,1) = a* ( -t*b*gy + gy*(3*b-c)/e2 )
          ylmd(i,2) = a* ( -t*gx*c + gx*(b-3*c)/e2 )
          ylmd(i,3) = a* ( -t*gx*gy*gz )
       enddo
    end if
  end subroutine m_pwBS_sphrp2_diff

! === necessary to make 3D_Parallel, too!!! by tkato ===========================
!!BRANCH_P ORG_Parallel
! ==============================================================================
  subroutine m_pwBS_kinetic_energies(ik,vkxyz,ekin)
    integer, intent(in)       :: ik
    real(kind=DP), intent(in) :: vkxyz(kv3,3,CRDTYP)
    real(kind=DP), intent(out):: ekin(kg1)

    integer i, nb
    real(kind=DP) :: ga,gb,gc
    integer       :: id_sname = -1
    call tstatc0_begin('m_pwBS_kinetic_energies ',id_sname)
    call getttr(rltv,ttr)

    do i = 1, iba(ik)
       nb = nbase(i,ik)
       ga = vkxyz(ik,1,BUCS) + ngabc(nb,1)
       gb = vkxyz(ik,2,BUCS) + ngabc(nb,2)
       gc = vkxyz(ik,3,BUCS) + ngabc(nb,3)
       ekin(i) = ( ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
            &  +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga)*0.5d0
    end do
    call tstatc0_end(id_sname)
  end subroutine m_pwBS_kinetic_energies
! === necessary to make 3D_Parallel, too!!! by tkato ===========================
!!BRANCH_P_END ORG_Parallel
! ==============================================================================

  ! ------- Positron start
  subroutine m_pwBS_pstrn_kinetic_energies(ekin)
    real(kind=DP), intent(out):: ekin(kg1_pwf)

    integer i, nb
    real(kind=DP) :: ga,gb,gc
    integer       :: id_sname = -1
    call tstatc0_begin('m_pwBS_pstrn_kinetic_energies ',id_sname)
    call getttr(rltv,ttr)

    do i = 1, kg1_pwf
       nb = i
       ga = ngabc(nb,1)
       gb = ngabc(nb,2)
       gc = ngabc(nb,3)
       ekin(i) = ( ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
            &  +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga)*0.5d0
    end do
    call tstatc0_end(id_sname)
  end subroutine m_pwBS_pstrn_kinetic_energies
  ! ------- Positron end

  subroutine m_pwBS_find_min_max_G(nfout,nspin,ng_max,ng_min,n_max,n_min)
    integer, intent(in)    :: nfout,nspin
    integer, intent(inout) :: ng_max, ng_min
    integer, intent(out)   :: n_max(3),n_min(3)

    integer             :: nb, ik, i

    if(ng_max == 0 .and. ng_min == 0) then
       nb = nbase(1,1)
       ng_max=max(ngabc(nb,1),ngabc(nb,2),ngabc(nb,3))
       ng_min=min(ngabc(nb,1),ngabc(nb,2),ngabc(nb,3))

       n_max(1)=ngabc(nb,1); n_max(2)=ngabc(nb,2); n_max(3)=ngabc(nb,3)
       n_min(1)=ngabc(nb,1); n_min(2)=ngabc(nb,2); n_min(3)=ngabc(nb,3)
       do ik = 1, kv3, nspin
          do i=1,iba(ik)
             nb = nbase(i,ik)
             ng_max=max(ng_max,ngabc(nb,1),ngabc(nb,2),ngabc(nb,3))
             ng_min=min(ng_min,ngabc(nb,1),ngabc(nb,2),ngabc(nb,3))
             n_max(1) = max(n_max(1), ngabc(nb,1))
             n_max(2) = max(n_max(2), ngabc(nb,2))
             n_max(3) = max(n_max(3), ngabc(nb,3))
             n_min(1) = min(n_min(1), ngabc(nb,1))
             n_min(2) = min(n_min(2), ngabc(nb,2))
             n_min(3) = min(n_min(3), ngabc(nb,3))
          end do
       end do
       if(ipri >= 1) then
          write(nfout,*) ' n_min1,n_max1 = ',n_min(1),n_max(1)
          write(nfout,*) ' n_min2,n_max2 = ',n_min(2),n_max(2)
          write(nfout,*) ' n_min3,n_max3 = ',n_min(3),n_max(3)
          write(nfout,*) ' ng_min,ng_max = ',ng_min,  ng_max
       end if
    end if
  end subroutine m_pwBS_find_min_max_G

  subroutine m_pwBS_decide_cutoff_mix
    if(cutoff_mix == SMALL) then
       kgpm = kg1_ext
    else if(cutoff_mix == MEDIUM) then
       kgpm = kg
    else if(cutoff_mix == LARGE) then
       kgpm  = kgp
    end if
    if(mype == 0) then
      !!$ print '(" !! cutoff_mix = ",i4)', cutoff_mix
      !!$ print '(" !! kgpm       = ",i8)', kgpm
    end if
  end subroutine m_pwBS_decide_cutoff_mix

  subroutine m_pwBS_wd_ngabc_etc(nfcntn_bin)
    integer,intent(in) :: nfcntn_bin

    integer,pointer,dimension(:) :: igfp_mpi,igfp_mpi2
    integer :: i

    if(mype == 0) write(nfcntn_bin) ngabc
    if(mype == 0) write(nfcntn_bin) igf
    if(npes > 1) then
       allocate(igfp_mpi(kgp));allocate(igfp_mpi2(kgp))
       igfp_mpi = 0
#ifdef _MPIFFT_
       do i = ista_kngp, iend_kngp
          igfp_mpi(i) = igfp_nonpara(i)
       end do
#else
       do i = ista_kngp, iend_kngp
          igfp_mpi(i) = igfp_l(i)
       end do
#endif
       call mpi_allreduce(igfp_mpi,igfp_mpi2,kgp,mpi_integer,mpi_sum &
            &            ,MPI_CommGroup,ierr)
       if(mype == 0) write(nfcntn_bin) igfp_mpi2
    else
#ifdef _MPIFFT_
       write(nfcntn_bin) igfp_nonpara
#else
       write(nfcntn_bin) igfp_l
#endif
    end if
    if(mype == 0) write(nfcntn_bin) kg,kgp
    if(mype == 0) write(nfcntn_bin) nbase
    if(mype == 0) write(nfcntn_bin) iba
    if(mype == 0) write(nfcntn_bin) ttr
  end subroutine m_pwBS_wd_ngabc_etc

  subroutine m_pwBS_increase_kg1(iadd)
    integer, intent(in) :: iadd
    kg1 = kg1 + iadd
    kg1_ext = kg1_ext + iadd
    if(ipri >= 1) write(nfout,'(" !pwBS: kg1 is increased to be ",i9)') kg1
    if(ipri >= 1) write(nfout,'(" !pwBS: kg1_ext is increased to be ",i9)') kg1_ext
  end subroutine m_pwBS_increase_kg1

  subroutine m_pwBS_cp_iba_to_iba_ek()
     integer :: i
     if(.not.allocated(iba_ek)) allocate(iba_ek(kv3_ek))
     iba_ek = iba

     if(ipri >= 3) then
        do i = 1, kv3_ek
           write(nfout,'(" iba_ek(",i6,") = ",i6 &
                & ," <<m_pwBS_cp_iba_to_iba_ek>>")') i, iba_ek(i)
        end do
     end if

  end subroutine m_pwBS_cp_iba_to_iba_ek

  subroutine m_pwBS_cp_iba_ek_to_iba(nk)
     integer, intent(in) :: nk
     integer :: kvt, i, nks, ikt, is, kv3t

     call mpi_barrier(MPI_CommGroup,ierr)
     if(allocated(iba)) then
        deallocate(iba)
        allocate(iba(kv3))
        if(ipri >= 3) write(nfout,'(" iba is allocated <<m_pwBS_cp_iba_ek_to_iba>>")')
     end if

     if(nk > kv3_ek ) then
        if(ipri >= 1) write(nfout,'(" nk = ",i6," > kv3_ek = ",i6)') nk, kv3_ek
        !stop ' nk > kv3_ek <<m_pwBS_cp_iba_ek_to_iba>>'
        write(nfout,'(a,2i10)') '!** nk > kv3_ek', nk, kv3_ek
     end if

     if(sw_modified_kpoint_increment == ON) then
       do i=0,nrank_k-1
         do is=1,nspin
           ikt = nspin*(nis_kv3_ek(i)-1)+(nkgroup-1)*nspin+is
           if(ikt>kv3_ek) ikt = kv3_ek
!           iba((i-1)*nspin+is) = iba_ek(ikt)
           iba(i*nspin+is) = iba_ek(ikt)
         enddo
       enddo
     else
       if(fixed_charge_k_parallel == ONE_BY_ONE .and. nk+kv3-1 > kv3_ek) then
          kvt = kv3_ek - nk + 1
          iba(1:kvt) = iba_ek(nk:kv3_ek)
          do i = kvt+1, kv3
             iba(i) = iba_ek(kv3_ek)
          end do
       else
          iba(1:kv3) = iba_ek(nk:nk+kv3-1)
       end if
     endif
     if(ipri >= 2) then
        do i = 1, kv3
           write(nfout,'(" iba(",i6,") = ",i6 &
                & ," <<m_pwBS_cp_iba_ek_to_iba>>")') i, iba(i)
        end do
     end if
   end subroutine m_pwBS_cp_iba_ek_to_iba

   subroutine m_pwBS_get_igfp(ista,iend,igfp,idp,idp2,nboxsize)
     integer, intent(in) :: ista,iend,idp,idp2,nboxsize(3)
     integer, intent(out),dimension(ista:iend) :: igfp

     integer :: igf1, igf2, igf3, i, nlp,nnm,nnn
     nlp = nboxsize(1)
     nnm = nboxsize(2)
     nnn = nboxsize(3)

     do i = ista, iend  !for mpi
        igf1 = ngabc(i,1) + 1
        igf2 = ngabc(i,2) + 1
        igf3 = ngabc(i,3) + 1
        if(ngabc(i,1) <= -1) igf1 = igf1 + nlp
        if(ngabc(i,2) <= -1) igf2 = igf2 + nnm
        if(ngabc(i,3) <= -1) igf3 = igf3 + nnn
        igfp(i) = igf1 + (igf2-1)*idp + (igf3-1)*idp*idp2
     end do
   end subroutine m_pwBS_get_igfp

   subroutine m_pwBS_dealloc(dealloc_all)
     logical, intent(in), optional :: dealloc_all
     logical :: dealloc, unitcell_can_change
     dealloc = unitcell_can_change()
     if(present(dealloc_all)) then
       dealloc = dealloc_all
     endif
     if(dealloc)then
        if(allocated(igf_prev)) deallocate(igf_prev)
        allocate(igf_prev(kg));igf_prev = igf

        if(allocated(nbase_prev)) deallocate(nbase_prev)
        allocate(nbase_prev(size(nbase,1),size(nbase,2)))
        nbase_prev = nbase

        if(allocated(nbase_gamma_prev)) deallocate(nbase_gamma_prev)
        if(allocated(nbase_gamma))then
          allocate(nbase_gamma_prev(size(nbase_gamma,1),size(nbase_gamma,2)))
          nbase_gamma_prev = nbase_gamma
        endif
        if(allocated(iba_prev)) deallocate(iba_prev)
        allocate(iba_prev(kv3));iba_prev = iba

        if(allocated(igfp_l_prev)) deallocate(igfp_l_prev)
        allocate(igfp_l_prev(ista_kngp:iend_kngp));igfp_l_prev = igfp_l
        kgp_reduced_prev = kgp_reduced
     endif
     if(allocated(ngabc)) deallocate(ngabc)

     if(allocated(ngpt_l)) deallocate(ngpt_l)
     if(allocated(igfp_l)) deallocate(igfp_l)
     if(allocated(igf)) deallocate(igf)
     if(allocated(gr_l)) deallocate(gr_l)

     if(allocated(ngshell_range)) deallocate(ngshell_range)

     if(allocated(nbase)) deallocate(nbase)

     if(allocated(ylm_l)) deallocate(ylm_l)

     if(allocated(nbase_gamma)) deallocate(nbase_gamma)

#ifdef _MPIFFT_
     if(allocated(igfp_l_c)) deallocate(igfp_l_c)
     if(allocated(igfp_nonpara)) deallocate(igfp_nonpara)
#endif

     if(allocated(iba)) deallocate(iba)

     kg1_ext = 0

   end subroutine m_pwBS_dealloc

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


  subroutine m_pwBS_sphrp_exx(is,rltv,ista_ylm,iend_ylm,gqmk,gqmkr,ylm)
    integer, intent(in) :: is
    real(kind=DP), intent(in), dimension(3,3)  :: rltv
    integer, intent(in) :: ista_ylm, iend_ylm
    real(kind=DP), intent(in), dimension(ista_ylm:iend_ylm,3) :: gqmk
    real(kind=DP), intent(inout), dimension(ista_ylm:iend_ylm) :: gqmkr
    real(kind=DP), intent(out), dimension(ista_ylm:iend_ylm) :: ylm
!!$C spherical harmonics.
!!$c     Some lines are rewritten by T.Yamasaki according to Y.Morikawa's
!!$c   e-mail. at 4th May 1994

    integer :: i,ni
    real(kind=DP) a,b,c,d,e,f,gx,gy,gz
    real(kind=DP) b1x,b1y,b1z,b2x,b2y,b2z,b3x,b3y,b3z
    real(kind=DP) :: gqmkr0

    b1x = rltv(1,1);  b1y = rltv(2,1);  b1z = rltv(3,1)
    b2x = rltv(1,2);  b2y = rltv(2,2);  b2z = rltv(3,2)
    b3x = rltv(1,3);  b3y = rltv(2,3);  b3z = rltv(3,3)

    ni = 0
    do i=ista_ylm,iend_ylm
       if(gqmkr(i) < DELTA) then
          ni = i
          gqmkr0 = gqmkr(i)
          gqmkr(i) = 1.d0
       end if
    end do

    if(is == 1) then
       a = dsqrt(1.d0/PAI4)
       do i = ista_ylm, iend_ylm  !for mpi
          ylm(i) = a
       end do
    else if(is == 2) then
       a = dsqrt(3.d0/PAI4)
       do i = ista_ylm, iend_ylm  !for mpi
          gx = b1x*gqmk(i,1) + b2x*gqmk(i,2) + b3x*gqmk(i,3)
          ylm(i) = a*gx/gqmkr(i)
       end do
    else if(is == 3) then
       a = dsqrt(3.d0/PAI4)
       do i = ista_ylm, iend_ylm  !for mpi
          gy =   b1y*gqmk(i,1) + b2y*gqmk(i,2) + b3y*gqmk(i,3)
          ylm(i) = a*gy/gqmkr(i)
       end do
    else if(is == 4) then
       a = dsqrt(3.d0/PAI4)
       do i = ista_ylm, iend_ylm  !for mpi
          gz =   b1z*gqmk(i,1) + b2z*gqmk(i,2) + b3z*gqmk(i,3)
          ylm(i) = a*gz/gqmkr(i)
       end do
    else if(is == 5) then
       a = dsqrt(5.d0/(16*PAI))
       do i = ista_ylm, iend_ylm  !for mpi
          gx = b1x*gqmk(i,1) + b2x*gqmk(i,2)  + b3x*gqmk(i,3)
          gy = b1y*gqmk(i,1) + b2y*gqmk(i,2)  + b3y*gqmk(i,3)
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2)  + b3z*gqmk(i,3)
          d = gz**2
          e = gx**2+gy**2+d
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*(3*d-e)/e
          if(abs(e) > tiny(e)) then
             ylm(i) = a*(3*d-e)/e
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 6) then
       a = dsqrt(15.d0/(16*PAI))
       do i = ista_ylm, iend_ylm  !for mpi
          gx = b1x*gqmk(i,1) + b2x*gqmk(i,2)  + b3x*gqmk(i,3)
          gy = b1y*gqmk(i,1) + b2y*gqmk(i,2)  + b3y*gqmk(i,3)
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2)  + b3z*gqmk(i,3)
          b = gx**2
          c = gy**2
          e = b+c+gz**2
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*(b-c)/e
          if(abs(e) > tiny(e)) then
             ylm(i) = a*(b-c)/e
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 7) then
       a = dsqrt(15.d0/PAI4)
       do i = ista_ylm, iend_ylm  !for mpi
          gx = b1x*gqmk(i,1) + b2x*gqmk(i,2) + b3x*gqmk(i,3)
          gy = b1y*gqmk(i,1) + b2y*gqmk(i,2) + b3y*gqmk(i,3)
          e = gqmkr(i)**2
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*gx*gy/e
          if(abs(e) > tiny(e)) then
             ylm(i) = a*gx*gy/e
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 8) then
       a = dsqrt(15.d0/PAI4)
       do i = ista_ylm, iend_ylm  !for mpi
          gy = b1y*gqmk(i,1) + b2y*gqmk(i,2)  + b3y*gqmk(i,3)
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2)  + b3z*gqmk(i,3)
          e = gqmkr(i)**2
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*gy*gz/e
          if(abs(e) > tiny(e)) then
             ylm(i) = a*gy*gz/e
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 9) then
       a = dsqrt(15.d0/PAI4)
       do i = ista_ylm, iend_ylm  !for mpi
          gx = b1x*gqmk(i,1) + b2x*gqmk(i,2) + b3x*gqmk(i,3)
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2) + b3z*gqmk(i,3)
          e = gqmkr(i)**2
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*gz*gx/e
          if(abs(e) > tiny(e)) then
             ylm(i) = a*gz*gx/e
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 10) then
       a = dsqrt(7.d0/(16*PAI))
       do i = ista_ylm, iend_ylm  !for mpi
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2) + b3z*gqmk(i,3)
          d = gz**2
          e = gqmkr(i)**2
          f = e * gqmkr(i)
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*gz*(5*d-3*e)/f
          if(abs(f) > tiny(f)) then
             ylm(i) = a*gz*(5*d-3*e)/f
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 11) then
       a = dsqrt(21.d0/(32*PAI))
       do i = ista_ylm, iend_ylm  !for mpi
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2)  + b3z*gqmk(i,3)
          gx = b1x*gqmk(i,1) + b2x*gqmk(i,2)  + b3x*gqmk(i,3)
          d = gz**2
          e = gqmkr(i)**2
          f = e * gqmkr(i)
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*gx*(5*d-e)/f
          if(abs(f) > tiny(f)) then
             ylm(i) = a*gx*(5*d-e)/f
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 12) then
       a = dsqrt(21.d0/(32*PAI))
       do i = ista_ylm, iend_ylm  !for mpi
          gy = b1y*gqmk(i,1) + b2y*gqmk(i,2) + b3y*gqmk(i,3)
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2) + b3z*gqmk(i,3)
          d = gz**2
          e = gqmkr(i)**2
          f = e * gqmkr(i)
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*gy*(5*d-e)/f
          if(abs(f) > tiny(f)) then
             ylm(i) = a*gy*(5*d-e)/f
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 13) then
       a = dsqrt(105.d0/(16*PAI))
       do i = ista_ylm, iend_ylm  !for mpi
          gx = b1x*gqmk(i,1) + b2x*gqmk(i,2) + b3x*gqmk(i,3)
          gy = b1y*gqmk(i,1) + b2y*gqmk(i,2) + b3y*gqmk(i,3)
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2) + b3z*gqmk(i,3)
          b = gx**2
          c = gy**2
          e = gqmkr(i)**2
          f = e * gqmkr(i)
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*gz*(b-c)/f
          if(abs(f) > tiny(f)) then
             ylm(i) = a*gz*(b-c)/f
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 14) then
       a = dsqrt(105.d0/PAI4)
       do i = ista_ylm, iend_ylm  !for mpi
          gx = b1x*gqmk(i,1) + b2x*gqmk(i,2)  + b3x*gqmk(i,3)
          gy = b1y*gqmk(i,1) + b2y*gqmk(i,2)  + b3y*gqmk(i,3)
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2)  + b3z*gqmk(i,3)
          f = gqmkr(i)**3
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*gx*gy*gz/f
          if(abs(f) > tiny(f)) then
             ylm(i) = a*gx*gy*gz/f
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 15) then
       a = dsqrt(35.d0/(32*PAI))
       do i = ista_ylm, iend_ylm  !for mpi
          gx = b1x*gqmk(i,1) + b2x*gqmk(i,2) + b3x*gqmk(i,3)
          gy = b1y*gqmk(i,1) + b2y*gqmk(i,2) + b3y*gqmk(i,3)
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2) + b3z*gqmk(i,3)
          b = gx**2
          c = gy**2
          f = gqmkr(i)**3
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*gx*(b-3*c)/f
          if(abs(f) > tiny(f)) then
             ylm(i) = a*gx*(b-3*c)/f
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 16) then
       a = dsqrt(35.d0/(32*pai))
       do i = ista_ylm, iend_ylm !for mpi
          gx = b1x*gqmk(i,1) + b2x*gqmk(i,2) + b3x*gqmk(i,3)
          gy = b1y*gqmk(i,1) + b2y*gqmk(i,2) + b3y*gqmk(i,3)
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2) + b3z*gqmk(i,3)
          b = gx**2
          c = gy**2
          f = gqmkr(i)**3
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*gy*(3*b-c)/f
          if(abs(f) > tiny(f)) then
             ylm(i) = a*gy*(3*b-c)/f
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 17) then
       a = 3.d0/8.d0/dsqrt(PAI4)
       do i = ista_ylm, iend_ylm  !for mpi
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2) + b3z*gqmk(i,3)
          d = gz**2
          e = gqmkr(i)**2
          f = e**2
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*(5*d*(7*d-6*e)/f+3.d0)
          if(abs(f) > tiny(f)) then
             ylm(i) = a*(5*d*(7*d-6*e)/f+3.d0)
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 18) then
       a = 15.d0/4.d0/dsqrt(10*PAI)
       do i = ista_ylm, iend_ylm  !for mpi
          gx = b1x*gqmk(i,1) + b2x*gqmk(i,2) + b3x*gqmk(i,3)
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2) + b3z*gqmk(i,3)
          d = gz**2
          e = gqmkr(i)**2
          f = e**2
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*gz*gx*(7*d-3*e)/f
          if(abs(f) > tiny(f)) then
             ylm(i) = a*gz*gx*(7*d-3*e)/f
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 19) then
       a = 15.d0/4.d0/dsqrt(10.d0*PAI)
       do i = ista_ylm, iend_ylm  !for mpi
          gx = b1x*gqmk(i,1) + b2x*gqmk(i,2)  + b3x*gqmk(i,3)
          gy = b1y*gqmk(i,1) + b2y*gqmk(i,2)  + b3y*gqmk(i,3)
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2)  + b3z*gqmk(i,3)
          d = gz**2
          e = gqmkr(i)**2
          f = e**2
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*gy*gz*(7*d-3*e)/f
          if(abs(f) > tiny(f)) then
             ylm(i) = a*gy*gz*(7*d-3*e)/f
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 20) then
       a = 15.d0/8.d0/dsqrt(5*PAI)
       do i = ista_ylm, iend_ylm  !for mpi
          gx = b1x*gqmk(i,1) + b2x*gqmk(i,2)  + b3x*gqmk(i,3)
          gy = b1y*gqmk(i,1) + b2y*gqmk(i,2)  + b3y*gqmk(i,3)
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2)  + b3z*gqmk(i,3)
          b = gx**2
          c = gy**2
          d = gz**2
          e = b+c+d
          f = e**2
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*(7*d-e)*(b-c)/f
          if(abs(f) > tiny(f)) then
             ylm(i) = a*(7*d-e)*(b-c)/f
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 21) then
       a = 15.d0/4.d0/dsqrt(5*PAI)
       do i = ista_ylm, iend_ylm  !for mpi
          gx = b1x*gqmk(i,1) + b2x*gqmk(i,2) + b3x*gqmk(i,3)
          gy = b1y*gqmk(i,1) + b2y*gqmk(i,2) + b3y*gqmk(i,3)
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2) + b3z*gqmk(i,3)
          d = gz**2
          e = gqmkr(i)**2
          f = e**2
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*(7*d-e)*gx*gy/f
          if(abs(f) > tiny(f)) then
             ylm(i) = a*(7*d-e)*gx*gy/f
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 22) then
       a = 105.d0/4.d0/dsqrt(70*pai)
       do i = ista_ylm, iend_ylm  !for mpi
          gx = b1x*gqmk(i,1) + b2x*gqmk(i,2) + b3x*gqmk(i,3)
          gy = b1y*gqmk(i,1) + b2y*gqmk(i,2) + b3y*gqmk(i,3)
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2) + b3z*gqmk(i,3)
          b = gx**2
          c = gy**2
          f = (b+c+gz**2)**2
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*(b-3*c)*gz*gx/f
          if(abs(f) > tiny(f)) then
             ylm(i) = a*(b-3*c)*gz*gx/f
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 23) then
       a = 105.d0/4.d0/dsqrt(70*PAI)
       do i = ista_ylm, iend_ylm  !for mpi
          gx = b1x*gqmk(i,1) + b2x*gqmk(i,2) + b3x*gqmk(i,3)
          gy = b1y*gqmk(i,1) + b2y*gqmk(i,2) + b3y*gqmk(i,3)
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2) + b3z*gqmk(i,3)
          b = gx**2
          c = gy**2
          f = (b+c+gz**2)**2
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*(3.d0*b-c)*gy*gz/f
          if(abs(f) > tiny(f)) then
             ylm(i) = a*(3.d0*b-c)*gy*gz/f
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 24) then
       a = 105.d0/16.d0/dsqrt(35*PAI)
       do i = ista_ylm, iend_ylm  !for mpi
          gx = b1x*gqmk(i,1) + b2x*gqmk(i,2)  + b3x*gqmk(i,3)
          gy = b1y*gqmk(i,1) + b2y*gqmk(i,2)  + b3y*gqmk(i,3)
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2)  + b3z*gqmk(i,3)
          b = gx**2
          c = gy**2
          f = (b+c+gz**2)**2
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*((b-c)**2-4*b*c)/f
          if(abs(f) > tiny(f)) then
             ylm(i) = a*((b-c)**2-4*b*c)/f
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    else if(is == 25) then
       a = 105.d0/4.d0/dsqrt(35*PAI)
       do i = ista_ylm, iend_ylm  !for mpi
          gx = b1x*gqmk(i,1) + b2x*gqmk(i,2)  + b3x*gqmk(i,3)
          gy = b1y*gqmk(i,1) + b2y*gqmk(i,2)  + b3y*gqmk(i,3)
          gz = b1z*gqmk(i,1) + b2z*gqmk(i,2)  + b3z*gqmk(i,3)
          b = gx**2
          c = gy**2
          f = (b+c+gz**2)**2
! === DEBUG by tkato 2015/03/19 ================================================
!         ylm(i) = a*(b-c)*gx*gy/f
          if(abs(f) > tiny(f)) then
             ylm(i) = a*(b-c)*gx*gy/f
          else
             ylm(i) = 0.0d0
          end if
! ==============================================================================
       end do
    end if

    if(ni > 0) then
       ylm(ni) = dsqrt(1.d0/PAI4)
       gqmkr(ni) = gqmkr0
    end if
  end subroutine m_pwBS_sphrp_exx


! ==== EXP_CELLOPT === 2015/09/24
  subroutine m_pwBS_store_prev_kg1_kgp
    kg1_prev = kg1;   kgp_prev = kgp
    if ( mype == 0 ) then
       write(nfout,*) '** kg1_prev is ', kg1_prev
       write(nfout,*) '** kgp_prev is ', kgp_prev
    endif
  end subroutine m_pwBS_store_prev_kg1_kgp
! =================== 2015/09/24

  subroutine m_pwBS_wd_curr_pws()
    use m_Files, only : m_Files_open_pwbs, m_Files_close_pwbs, nfpwbs
    integer, allocatable, dimension(:) :: igfp_l_t
    integer :: i,ierrr
    allocate(igfp_l_t(kgp));igfp_l_t=0
    do i=ista_kngp, iend_kngp
       igfp_l_t(i) = igfp_l(i)
    enddo
    call mpi_allreduce(mpi_in_place,igfp_l_t,kgp,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
    call m_Files_open_pwbs()
    if(mype==0)then
    write(nfpwbs) kg1,kg,size(nbase,1),kgp_reduced,neg
    write(nfpwbs) fft_box_size_WF
    write(nfpwbs) nbase
    if(kg_gamma>0)then
    write(nfpwbs) kg_gamma
    write(nfpwbs) nbase_gamma
    endif
    write(nfpwbs) igf
    write(nfpwbs) iba
    write(nfpwbs) fft_box_size_CD
    write(nfpwbs) fft_box_size_CD_nonpara
    write(nfpwbs) kgp
    write(nfpwbs) igfp_l_t
    write(nfpwbs) univol
    endif
    call m_Files_close_pwbs()
    deallocate(igfp_l_t)
  end subroutine m_pwBS_wd_curr_pws

  subroutine m_pwBS_rd_prev_pws(rd)
    use m_Files, only : m_Files_open_pwbs, m_Files_close_pwbs, nfpwbs, F_PWBS
    logical, intent(out) :: rd
    integer :: i,ierr,kgtmp,kggammatmp,kgtmp2
    integer, allocatable, dimension(:) :: igfp_l_t
    logical :: logi
    rd = .false.
    inquire(file=F_PWBS,exist=logi)
    if(.not. logi) then
      if(printable) write(nfout,'(a)') '!** pwbs file '//F_PWBS//' does not exist'
      return
    endif
    call m_Files_open_pwbs()
    if(mype==0)then
    read(nfpwbs,err=100,end=100) kg1_prev,kgtmp,kgtmp2,kgp_reduced_prev,neg_previous
    read(nfpwbs,err=100,end=100) fft_box_size_WF_prev
    endif
    call mpi_bcast(kg1_prev,1,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(kgtmp,1,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(kgtmp2,1,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(kgp_reduced_prev,1,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(neg_previous,1,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(fft_box_size_WF_prev,6,mpi_integer,0,MPI_CommGroup,ierr)

    if(allocated(nbase_prev)) deallocate(nbase_prev)
    allocate(nbase_prev(kgtmp2,kv3))
    if(allocated(igf_prev)) deallocate(igf_prev)
    allocate(igf_prev(kgtmp))
    if(allocated(iba_prev)) deallocate(iba_prev)
    allocate(iba_prev(kv3))
    if(mype==0) then
      read(nfpwbs,err=100,end=100) nbase_prev
    endif
    call mpi_bcast(nbase_prev,kgtmp2*kv3,mpi_integer,0,MPI_CommGroup,ierr)
    if(kg_gamma>0)then
      if(mype==0) then
        read(nfpwbs,err=100,end=100) kggammatmp
      endif
      call mpi_bcast(kggammatmp,1,mpi_integer,0,MPI_CommGroup,ierr)
      if(allocated(nbase_gamma_prev)) deallocate(nbase_gamma_prev)
      allocate(nbase_gamma_prev(kggammatmp,2))
      if(mype==0)then
        read(nfpwbs,err=100,end=100) nbase_gamma_prev
      endif
      call mpi_bcast(nbase_gamma_prev,kggammatmp*2,mpi_integer,0,MPI_CommGroup,ierr)
    endif
    if(mype==0) then
      read(nfpwbs,err=100,end=100) igf_prev
      read(nfpwbs,err=100,end=100) iba_prev
    endif
    call mpi_bcast(igf_prev,kgtmp,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(iba_prev,kv3,mpi_integer,0,MPI_CommGroup,ierr)

    if(mype==0)then
      read(nfpwbs,err=100,end=100) fft_box_size_CD_prev
      read(nfpwbs,err=100,end=100) fft_box_size_CD_nonpara_prev
      read(nfpwbs,err=100,end=100) kgp_prev
    endif
    call mpi_bcast(fft_box_size_CD_prev,6,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(fft_box_size_CD_nonpara_prev,3,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(kgp_prev,1,mpi_integer,0,MPI_CommGroup,ierr)
    call m_Parallel_init_mpi_kngp_prev(kgp_prev)
    allocate(igfp_l_t(kgp_prev));igfp_l_t=0
    if(allocated(igfp_l_prev)) deallocate(igfp_l_prev)
    allocate(igfp_l_prev(ista_kngp_prev:iend_kngp_prev));igfp_l_prev=0
    if(mype==0) read(nfpwbs,err=100,end=100) igfp_l_t
    call mpi_bcast(igfp_l_t,kgp_prev,mpi_integer,0,MPI_CommGroup,ierr)
    do i=ista_kngp_prev,iend_kngp_prev
       igfp_l_prev(i) = igfp_l_t(i)
    enddo
    deallocate(igfp_l_t)
    if(mype==0) read(nfpwbs,err=100,end=100) univol_prev
    call mpi_bcast(univol_prev,1,mpi_double_precision,0,MPI_CommGroup,ierr)

    call m_Files_close_pwbs()
    if(printable) write(nfout,'(a)') ' !** read pwbs from file'
    rd = .true.
    return
100 continue
    if(mype==0) write(nfout,'(a)') ' !** encountered error while reading data from '//F_PWBS
  end subroutine m_pwBS_rd_prev_pws

end module m_PlaneWaveBasisSet
