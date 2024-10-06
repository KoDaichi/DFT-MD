!!$#define DEBUG_WRITE
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 597 $)
!
!  MODULE:  m_ES_WF_by_submat
!
!  AUTHOR(S): Tsuyoshi Miyazaki, T. Yamasaki   August/20/2003
!             Takenori Yamamoto June/22/2005
!  FURTHER MODIFICATION: T. Yamasaki, September 2009
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
!
#ifdef __TIMER__
#   define __TIMER_START_w_BARRIER(str,a)  call mpi_barrier(str,ierr); call timer_sta(a)
#   define __TIMER_STOP_w_BARRIER(str,a)   call timer_sta(a); call mpi_barrier(str,ierr)
#   define __TIMER_START(a)                call timer_sta(a)
#   define __TIMER_STOP(a)                 call timer_end(a)
#else
#   define __TIMER_START_w_BARRIER(str,a)
#   define __TIMER_STOP_w_BARRIER(str,a)
#   define __TIMER_START(a)
#   define __TIMER_STOP(a)
#endif
#ifdef __TIMER_SUB__
#   define __TIMER_SUB_START(a)  call timer_sta(a)
#   define __TIMER_SUB_STOP(a)   call timer_end(a)
#else
#   define __TIMER_SUB_START(a)
#   define __TIMER_SUB_STOP(a)
#endif
#ifdef __TIMER_DO__
#   define __TIMER_DO_START(a)  call timer_sta(a)
#   define __TIMER_DO_STOP(a)   call timer_end(a)
#else
#   define __TIMER_DO_START(a)
#   define __TIMER_DO_STOP(a)
#endif
#ifdef __TIMER_DGEMM__
#   define __TIMER_DGEMM_START(a)  call timer_sta(a)
#   define __TIMER_DGEMM_STOP(a)   call timer_end(a)
#else
#   define __TIMER_DGEMM_START(a)
#   define __TIMER_DGEMM_STOP(a)
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

module m_ES_WF_by_submat
! $Id: m_ES_WF_by_submat.F90 597 2019-07-26 03:08:56Z jkoga $

!!$#ifdef TRANSPOSE
#ifdef VPP
#define _ODD_BOUNDARY_
#endif
#ifdef SX
#define _ODD_BOUNDARY_
#endif
#ifdef NEC_TUNE1
#define _ODD_BOUNDARY_
#endif

#ifndef NO_SUBMAT_DGEMM
#define SUBMAT_DGEMM
#endif

  use m_Const_Parameters,    only : DP,SP,DIRECT,ON,OFF,SCF, SmallestPositiveNumber,GAMMA &
       &                          , HOUSEHOLDER, DIVIDEandCONQUER, ELECTRON, CMPLDP, OLD, INVERSE, EXECUT
  use m_Parallelization,     only : MPI_CommGroup &
       &                          , myrank_e,myrank_k,map_e,map_k,ista_e,iend_e,istep_e &
       &                          , ista_k,iend_k,np_g1k,mpi_k_world,map_z &
       &                          , np_e,npes, nrank_k, nrank_e,mype &
       &                          , np_fs, mpi_spin_group, ista_spin, iend_spin, nrank_s
  use m_Control_Parameters,  only : nspin,iprisubmat,ldiag,kimg,neg,af &
       &                          , submat_critical_ratio,printable &
       &                          , m_CtrlP_set_submat, sw_scalapack &
       &                          , method_scalapack, block_size, nprow, npcol &
       &                          , msize_submat, sw_hybrid_functional, sw_fef &
       &                          , nblocksize_submat_is_given, nblocksize_submat, nb_submat_default &
       &                          , sw_gep &
       &                          , submat_uncalled &
       &                          , sw_serial_fft                    &
       &                          , sw_keep_hloc_phi &
#ifdef SAVE_FFT_TIMES
       &                          , nb_mgs_default, sw_save_fft
#else
       &                          , nb_mgs_default
#endif
  use m_Files,               only : nfout
  use m_Timing,              only : tstatc0_begin, tstatc0_end
  use m_FFT,                 only : nfft,fft_box_size_WF
  use m_FFT,                 only : m_FFT_Vlocal_W, m_FFT_WF
  use m_Kpoints,             only : kv3,vkxyz, k_symmetry
  use m_PlaneWaveBasisSet,   only : kg1, iba, igf, nbase, m_pwBS_kinetic_energies
  use m_Electronic_Structure,only : zaj_l, zaj_l_buf,neordr,           nrvf_ordr &
#ifdef SAVE_FFT_TIMES
       &                          , status_saved_phifftr &
#endif
       &                          , eko_l,vnlph_l
  use m_Electronic_Structure,only : m_ES_Vlocal_in_Rspace &
       &                          , m_ES_WF_in_Rspace &
       &                          , m_ES_wd_zaj_small_portion &
       &                          , m_ES_wd_eko &
       &                          , hlocphi_l, vtau_l
  use m_ES_nonlocal,         only : sc ,ss ,qc                                        &
 &                                , m_ES_Vnonlocal_W                                  &
 &                                , m_ES_betar_dot_WFs_4_each_k                       &
 &                                , m_ES_alloc_scss_etc                               &
 &                                , m_ES_dealloc_scss_etc
  use m_ES_ortho,            only : np_g1k_x
  use m_ES_ortho,            only : m_ES_W_transpose_r                                &
 &                                , m_ES_W_transpose_back_r                           &
 &                                , m_ES_W_transpose_back2                            &
 &                                , m_ES_W_transpose_fft_cmplx                        &
 &                                , m_ES_W_transpose_back_fft_cmplx

  use m_ES_ExactExchange,    only : m_ES_Vexx_W
  use m_FiniteElectricField, only : m_FEF_add_grad_to_vnlph

! ================================= added by K. Tagami ============== 11.0
  use m_Control_Parameters,    only : ndim_spinor, ndim_chgpot
  use m_Electronic_Structure,  only : m_ES_Vlocal_in_Rspace_noncl
  use m_FFT,                   only : m_FFT_Vlocal_W_noncl
  use m_Electronic_Structure,      only : m_ES_sort_eigen_vals_noncl
! =================================================================== 11.0

  use m_Control_Parameters,  only : use_metagga, vtau_exists
  use m_ES_WF_by_SDorCG, only : m_ES_contrib_kindens_to_vnlph

  use m_Control_Parameters, only : sw_rsb
  use mpi

  implicit none
  logical, allocatable, dimension(:) :: non_diagonal_part_is_small

  integer, save :: npart,npart2
  integer, allocatable, save :: isp(:),iep(:)
  integer, allocatable, save :: isp2(:),iep2(:)

  real(kind=DP), allocatable, target, dimension(:,:) :: unimat,unimat_h

  integer, private, save :: lda,occ
#ifdef SINGLE_CONTEXT
   integer, private, allocatable, save :: usermap(:,:)
#else
   integer, private, allocatable, save :: usermap(:,:,:) ! d(nprow,npcol,0:nrank_k-1)
#endif

   logical :: disable_utransform = .false.
   integer :: ierr

#ifdef DEBUG_WRITE
  integer, parameter :: DEBUGPRINTLEVEL = 1
#else
  integer, parameter :: DEBUGPRINTLEVEL = 2
#endif
  
!  include 'mpif.h'
contains
  subroutine m_ESsubmat_alloc()
    if(.not.allocated(non_diagonal_part_is_small)) then
       allocate(non_diagonal_part_is_small(kv3))
!!$       non_diagonal_part_is_small = .false.
    end if
    non_diagonal_part_is_small = .false.
  end subroutine m_ESsubmat_alloc

  subroutine m_ESsubmat_dealloc()
    if(allocated(non_diagonal_part_is_small)) deallocate(non_diagonal_part_is_small)
  end subroutine m_ESsubmat_dealloc

  subroutine m_ESsubmat_reset_non_diagon()
    if(allocated(non_diagonal_part_is_small)) non_diagonal_part_is_small = .false.
  end subroutine m_ESsubmat_reset_non_diagon

  subroutine m_ESsubmat_Renew_WF(nfout,meg,damp)
    integer, intent(in) :: nfout,meg
    real(kind=DP)       :: damp

    integer             :: ispin, iksnl, ik, ie
    real(kind=DP), dimension(:), pointer ::  afft
    real(kind=DP), dimension(:), allocatable :: bfft, cfft
    real(kind=DP), allocatable, dimension(:,:), target :: afft_spin
    real(kind=DP), pointer, dimension(:) :: ekin, p, vnldi
    real(kind=DP), allocatable, dimension(:,:) :: ekotmp

    integer :: n_all_kpoints, n_submat, ipri0
    logical :: submat_is_done

    submat_is_done = .false.
    call m_ES_alloc_scss_etc()

    if(nrank_s==1) then
      allocate(afft(nfft))
    endif
    allocate(bfft(nfft))
    if(nrank_s==2) allocate(afft_spin(nfft,2))
    if ( use_metagga .and. vtau_exists ) allocate( cfft(nfft) )

    ekin => sc(1:kg1); p => ss(1:kg1); vnldi => qc(1:kg1)

    n_all_kpoints = 0
    n_submat = 0
!!$    if(m_CtrlP_ntcnvg_clear()) non_diagonal_part_is_small = .false.
    if(iprisubmat >=DEBUGPRINTLEVEL) write(nfout,'(" !! <<m_ESsubmat_Renew_WF>>")')
    call get_ipri0(iprisubmat,ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko(nfout,mode=SCF)

!    do ispin = 1, nspin, (af+1)
    if(nrank_s==2) then
      do ispin = 1, nspin, (af+1)
        call m_ES_Vlocal_in_Rspace(ispin,afft_spin(:,ispin))     ! (ptfft1) vlhxc_l->afft
      enddo
    endif
    do ispin = ista_spin, iend_spin, (af+1)
       if(nrank_s==1) then
         call m_ES_Vlocal_in_Rspace(ispin,afft)      ! (ptfft1) vlhxc_l->afft
         if ( use_metagga .and. vtau_exists ) then
            call m_ES_Vlocal_in_Rspace( ispin, cfft, vtau_l )    ! r space
         endif
       else
         afft => afft_spin(:,ispin)
       endif

       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) /= myrank_k) cycle          ! MPI
          n_all_kpoints = n_all_kpoints + 1
          if(non_diagonal_part_is_small(ik)) then
             if(iprisubmat >=DEBUGPRINTLEVEL) write(nfout,'(" !!! submat is not done for ik=",i4)') ik
             cycle
          else
             if(iprisubmat >=DEBUGPRINTLEVEL) write(nfout,'(" !!! submat is done for ik=",i4)') ik
             n_submat = n_submat + 1
             submat_is_done = .true.
          end if
          iksnl = (ik-1)/nspin + 1
          call m_ES_Vnonlocal_W(ik,iksnl,ispin,switch_of_eko_part=OFF)  ! (nonloc) ->(vnlph_l)
          if(sw_hybrid_functional==ON) call m_ES_Vexx_W(ik) ! (exx) ->(vnlph_l)
          if(sw_fef==ON) call m_FEF_add_grad_to_vnlph(ik) !->vnlph_l

          if ( use_metagga .and. vtau_exists ) then
             call m_ES_contrib_kindens_to_vnlph( ispin, ik, cfft ) ! -> vnlph_l
          endif
          call m_pwBS_kinetic_energies(ik,vkxyz,ekin) ! (diakin) ->ekin
          call evolve_WFs_in_subspace&     !-(m_ES_WF_by_submat)
                  &(ik,iksnl,ispin,meg,damp,ekin,afft,bfft)

          if(ik==1.and.iprisubmat>= 2) &
               & call m_ES_wd_zaj_small_portion(nfout,1," -- after subspace-rotation --",30)
          call m_ES_betar_dot_WFs_4_each_k(nfout,ik)   ! -> fsr_l,fsi_l
       enddo      ! k-point loop
    enddo      ! spin loop

    if ( allocated(cfft) ) deallocate( cfft )

!!$    if(npes > 1) then
!!$       call mpi_allreduce(n_submat,n_submat_mpi,1,mpi_integer,mpi_sum,mpi_k_world(myrank_k),ierr)
!!$       n_submat = n_submat_mpi
!!$    end if

!!  ( in case of af=1 )
    if(af /= 0) then
       call cp_eigen_values_for_af       !-(contained here)
       call expand_neordr_and_nrvf_ordr  !-(contained here)
    end if

    if(ipri0 >= 2) call m_ES_wd_eko(nfout,mode=SCF)
    deallocate(bfft)
    if(nrank_s==1) then
      deallocate(afft)
    else
      deallocate(afft_spin)
    endif
    call m_ES_dealloc_scss_etc()
    if(iprisubmat >= 2) then
       if(submat_is_done) write(nfout,'(" !solver -- subspace_rotation is done for ",i4," of ",i4," kpoints")') &
            & n_submat, n_all_kpoints
       if(.not.submat_is_done) write(nfout,'(" !solver --  no subspace_rotation")')
    end if

    if(nrank_k > 1) then
       call mpi_allreduce(MPI_IN_PLACE, n_submat,1,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
       call mpi_allreduce(MPI_IN_PLACE, n_all_kpoints,1,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
!       call mpi_allreduce(MPI_IN_PLACE, n_submat,1,mpi_integer,mpi_sum,mpi_spin_group,ierr)
!       call mpi_allreduce(MPI_IN_PLACE, n_all_kpoints,1,mpi_integer,mpi_sum,mpi_spin_group,ierr)
    end if

    if(submat_is_done .and. n_submat < n_all_kpoints) submat_is_done = .false.
    call m_CtrlP_set_submat(submat_is_done)

  contains
    subroutine get_ipri0(ipri_in, ipri_out)
      integer, intent(in)  :: ipri_in
      integer, intent(out) :: ipri_out
      if(npes > 1) then
         if(mype == 0) ipri_out = ipri_in
         call mpi_bcast(ipri_out,1,mpi_integer,0,MPI_CommGroup,ierr)
      else
         ipri_out = ipri_in
      end if
    end subroutine get_ipri0

    subroutine cp_eigen_values_for_af
      integer :: ik,ib
      do ik = 1, kv3, af+1
         if(map_k(ik) /= myrank_k) cycle    ! MPI
         do ib = 1, np_e                    ! MPI
            eko_l(ib,ik+af) = eko_l(ib,ik)
         enddo
      enddo
    end subroutine cp_eigen_values_for_af

    subroutine expand_neordr_and_nrvf_ordr
      integer :: ik
      do ik = 1, kv3, af+1
         if(map_k(ik) /= myrank_k) cycle     ! MPI
         neordr(1:neg,ik+af) = neordr(1:neg,ik)
         nrvf_ordr(1:neg,ik+af) = nrvf_ordr(1:neg,ik)
      end do
    end subroutine expand_neordr_and_nrvf_ordr

  end subroutine m_ESsubmat_Renew_WF


! ===================================== added by K. Tagami ================ 11.0
  subroutine m_ESsubmat_Renew_WF_noncl(nfout,meg,damp)
    integer, intent(in) :: nfout,meg
    real(kind=DP)       :: damp

    integer             :: is, iksnl, ik

    real(kind=DP), allocatable ::  afft_kt(:,:)
    real(kind=DP), allocatable ::  bfft_kt(:,:)

    real(kind=DP), pointer, dimension(:) :: ekin, p, vnldi

    real(kind=DP), allocatable :: vnlph_noncl(:,:,:,:)

    integer :: n_all_kpoints, n_submat, ipri0
    logical :: submat_is_done

    integer :: is1, is2, k2, istmp

    submat_is_done = .false.

    call m_ES_alloc_scss_etc()
    allocate(afft_kt(nfft,ndim_chgpot))
    allocate(bfft_kt(nfft,ndim_spinor))

    ekin => sc(1:kg1); p => ss(1:kg1); vnldi => qc(1:kg1)

    allocate(vnlph_noncl(kg1,np_e,kimg,ndim_spinor)); vnlph_noncl = 0.0d0

    n_all_kpoints = 0;     n_submat = 0

    if (iprisubmat >= 2) write(nfout,'(" !! <<m_ESsubmat_Renew_WF_noncl>>")')
    call get_ipri0(iprisubmat,ipri0)

    if(ipri0 >= 2) call m_ES_wd_eko(nfout,mode=SCF)

    call m_ES_Vlocal_in_Rspace_noncl( afft_kt )    ! (ptfft1) vlhxc_l->afft

    Do ik = 1, kv3, ndim_spinor
       if(map_k(ik) /= myrank_k) cycle          ! MPI
       iksnl = (ik-1) /ndim_spinor + 1

       n_all_kpoints = n_all_kpoints + 1

       if ( non_diagonal_part_is_small(ik) ) then
          if(iprisubmat >= 2) write(nfout,'(" !!! submat is not done for ik=",i4)') ik
          cycle
       else
          if(iprisubmat >= 2) write(nfout,'(" !!! submat is done for ik=",i4)') ik
          n_submat = n_submat + 1
          submat_is_done = .true.
       end if

       vnlph_noncl = 0.0d0

       Do is1=1, ndim_spinor
          Do is2=1, ndim_spinor
            istmp = ( is1-1 )*ndim_spinor +is2
            k2 = ik +is2 -1

            call m_ES_Vnonlocal_W( k2, iksnl, istmp, switch_of_eko_part=OFF )
!!!!           call m_ES_Vnonlocal_W( k2, iksnl, istmp, switch_of_eko_part=ON )
                                                   ! (nonloc) ->(vnlph_l)
            vnlph_noncl(:,:,:,is1) = vnlph_noncl(:,:,:,is1) &
             &                     + vnlph_l(:,:,:)
          End do
       End do

       if ( sw_hybrid_functional==ON ) call m_ES_Vexx_W(ik) ! (exx) ->(vnlph_l)

       call m_pwBS_kinetic_energies( ik, vkxyz, ekin )       ! (diakin) ->ekin

       call evolve_WFs_in_subspace_noncl( ik, meg, damp, ekin, afft_kt, bfft_kt, &
	&                                 vnlph_noncl )

       if ( ik==1 .and. iprisubmat>= 2 ) then
          Do is1=1, ndim_spinor
            call m_ES_wd_zaj_small_portion( nfout, ik+is1-1, &
                   &                        " -- after subspace-rotation --",30 )
          End do
       endif
       Do is1=1, ndim_spinor
          call m_ES_betar_dot_WFs_4_each_k( nfout,ik+is1-1 )   ! -> fsr_l,fsi_l
       End do
    End do

! ****************************************** This is important ******
    call m_ES_sort_eigen_vals_noncl()
! *******************************************************************

    if(ipri0 >= 2) call m_ES_wd_eko(nfout,mode=SCF)
    deallocate(bfft_kt);   deallocate(afft_kt)

    call m_ES_dealloc_scss_etc()

    if (iprisubmat >= 2) then
       if(submat_is_done) write(nfout,'(" !solver -- subspace_rotation is done for ",i4," of ",i4," kpoints")') &
            & n_submat, n_all_kpoints
       if(.not.submat_is_done) write(nfout,'(" !solver --  no subspace_rotation")')
    end if

    if (nrank_k > 1) then
!       call mpi_allreduce(MPI_IN_PLACE, n_submat,     1,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
!       call mpi_allreduce(MPI_IN_PLACE, n_all_kpoints,1,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
       call mpi_allreduce(MPI_IN_PLACE, n_submat,     1,mpi_integer,mpi_sum,mpi_spin_group,ierr)
       call mpi_allreduce(MPI_IN_PLACE, n_all_kpoints,1,mpi_integer,mpi_sum,mpi_spin_group,ierr)
    end if

    if (submat_is_done .and. n_submat < n_all_kpoints) submat_is_done = .false.
    call m_CtrlP_set_submat(submat_is_done)

  contains
    subroutine get_ipri0(ipri_in, ipri_out)
      integer, intent(in)  :: ipri_in
      integer, intent(out) :: ipri_out
      if(npes > 1) then
         if(mype == 0) ipri_out = ipri_in
         call mpi_bcast(ipri_out,1,mpi_integer,0,MPI_CommGroup,ierr)
      else
         ipri_out = ipri_in
      end if
    end subroutine get_ipri0

  end subroutine m_ESsubmat_Renew_WF_noncl
!======================================================================= 11.0

 subroutine evolve_WFs_in_subspace&
       &(ik,iksnl,ispin,meg,damp,ekin,afft,bfft)
!
! Revised by T. Yamasaki, September 2009 : SUBMAT_DGEMM
!
    integer, intent(in) :: ik,iksnl,ispin,meg
    real(kind=DP), intent(in)  :: damp,ekin(kg1)
    real(kind=DP), intent(in)  :: afft(nfft)
    real(kind=DP), intent(out) :: bfft(nfft)
! (allocatable variables)
    real(kind=DP), allocatable,dimension(:) ::     eig
    real(kind=DP), allocatable,dimension(:) ::     eko_d
    real(kind=DP), allocatable,dimension(:) ::   w1hw2
    real(kind=DP), allocatable,dimension(:,:) ::   zz
    real(kind=DP), allocatable,dimension(:,:,:) :: zat_l
    real(kind=DP), allocatable,dimension(:,:,:) :: zah_l
    real(kind=DP), allocatable,dimension(:,:,:) :: zahloc_l
    real(kind=DP), allocatable,dimension(:,:,:) :: zat_wk
    real(kind=DP), allocatable,dimension(:,:,:) :: zahloc_wk

    real(kind=DP), allocatable :: amat(:,:), zmat(:,:)
!    integer, save :: lwork,lrwork,liwork,lda,occ
    integer, save :: lwork,lrwork,liwork
#ifdef _USE_SCALAPACK_
    integer, save :: lwork1,lrwork1,liwork1
    integer, save :: lwork2,lrwork2,liwork2
#endif
#ifdef SINGLE_CONTEXT
    integer, save :: ictxt, myrow, mycol
#else
    integer, save :: myrow, mycol, tmpctxt
    integer, allocatable, save :: ictxt(:) ! d(0:nrank_k-1)
#endif
    integer, dimension(9), save :: desca,descz
!#ifdef SINGLE_CONTEXT
!    integer, allocatable, save :: usermap(:,:)
!#else
!    integer, allocatable, save :: usermap(:,:,:) ! d(nprow,npcol,0:nrank_k-1)
!#endif
!!$    integer, save :: npart,npart2
!!$    integer, allocatable, save :: isp(:),iep(:)
!!$    integer, allocatable, save :: isp2(:),iep2(:)
!    logical, save :: first = .true.
    integer :: is,ie,m,jp
    integer :: msize

    integer       :: ib1,ib2,ib1to,i1,ii,ri,ib, i, kimg_t, ig1
    real(kind=DP) :: denom, eko1, eko2, ekod
    real(kind=DP) :: dr1,dr2,di1,di2, dd
    real(kind=DP) :: sum_sq_diagonal, sum_sq_non_diagonal &
         & , sum_abs_diagonal, sum_abs_non_diagonal
    integer :: id_sname = -1, id_sname2 = -1, ipri0
    integer :: max_block_size
#ifdef SUBMAT_DGEMM
    real(kind=DP),allocatable,dimension(:,:) :: w1hw2r,w1hw2i
    real(kind=DP) :: alpha, beta
#endif
! <-- T.Kokubo & D.Fukata, Feb. 2010
    integer :: ibsta, ibend, nblk1, nblk2, ibsize, ibsize2
    if(nblocksize_submat_is_given) then
      ibsize = nblocksize_submat
    else
      ibsize = nb_submat_default
    endif
! -->
    call tstatc0_begin('evolve_WFs_in_subspace ', id_sname,1)

    denom = 1.d0/product(fft_box_size_WF(1:3,1))

! (allocate)  zat_l will be divided by G-index
    if(k_symmetry(ik) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = kimg
    end if

    allocate(eko_d(neg));     eko_d = 0.d0
    allocate(eig(meg)); eig=0.d0
#ifdef _ODD_BOUNDARY_
    if(mod(np_g1k(ik),2) == 0) then
       np_g1k_x = np_g1k(ik) + 1
    else
       np_g1k_x = np_g1k(ik)
    end if
#else
    np_g1k_x = np_g1k(ik)
#endif
    allocate(zat_l(np_g1k_x,neg,kimg))
    !!$allocate(zah_l(np_g1k_x,neg,kimg))

    do ib1 = 1, neg
       if(map_e(ib1) == myrank_e) eko_d(ib1) = eko_l(map_z(ib1),ik)
    end do
    call mpi_allreduce(MPI_IN_PLACE,eko_d,neg,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
    eko1 = sum(eko_d(1:neg))

!( tenchi ) (zat_l <- zaj_l)
! --> T. Yamasaki, 28th Aug. 2009
!!$    call m_ES_W_transpose(ista_k,iend_k,ik,zaj_l,zat_l)
    call m_ES_W_transpose_r(.true.,ista_k,iend_k,ik,zaj_l,zat_l)
! <--

! (zaj_l <- H |phi> )
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(bfft,ib1to,ib,ii,i1,dr1,di1,dr2)
!CDIR CONCUR(BY=1)
#endif
    do ib1 = ista_e, iend_e, istep_e     ! MPI
       ib1to=nrvf_ordr(ib1,ik)
       if(ib1to>meg) cycle
       call m_ES_WF_in_Rspace(ik,ib1,bfft)!(swffft)
       call m_FFT_Vlocal_W(afft,bfft)    ! (afft, bfft)-> (bfft)
       call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
       ib = map_z(ib1)                 ! MPI
       if(kimg == 1) then
          do ii=1,iba(ik)
             i1  = igf(nbase(ii,ik))
             dr1 = zaj_l(ii,ib,ik,1)
             dr2 = bfft(i1)*denom
             zaj_l(ii,ib,ik,1)= ekin(ii)*dr1+vnlph_l(ii,ib,1)+dr2
             if(sw_keep_hloc_phi==ON) hlocphi_l(ii,ib,ik,1) = ekin(ii)*dr1+dr2
          enddo
       else
          do ii=1,iba(ik)
             i1  = igf(nbase(ii,ik))
             dr1  = zaj_l(ii,ib,ik,1)
             di1  = zaj_l(ii,ib,ik,2)
             zaj_l(ii,ib,ik,1)= ekin(ii)*dr1+bfft(2*i1-1)*denom+vnlph_l(ii,ib,1)
             zaj_l(ii,ib,ik,2)= ekin(ii)*di1+bfft(2*i1  )*denom+vnlph_l(ii,ib,2)
             if(sw_keep_hloc_phi==ON) then
               hlocphi_l(ii,ib,ik,1)= ekin(ii)*dr1+bfft(2*i1-1)*denom
               hlocphi_l(ii,ib,ik,2)= ekin(ii)*di1+bfft(2*i1  )*denom
             endif
          enddo
       end if
#ifdef SAVE_FFT_TIMES
       if(sw_save_fft == ON) status_saved_phifftr(ib,ik) = OLD
#endif
    end do

!( tenchi ) (zah_l <- zaj_l)
    allocate(zah_l(np_g1k_x,neg,kimg))                                 ! MPI
    allocate(zahloc_l(np_g1k_x,neg,kimg))                                 ! MPI
! --> T. Yamasaki 28th Aug. 009
!!$    call m_ES_W_transpose(ista_k,iend_k,ik,zaj_l,zah_l)
    call m_ES_W_transpose_r(.true.,ista_k,iend_k,ik,zaj_l,zah_l)
    if(sw_keep_hloc_phi==ON) &
    call m_ES_W_transpose_r(.true.,ista_k,iend_k,ik,hlocphi_l,zahloc_l)
! <--

! (make matrix elements )
#ifdef _USE_SCALAPACK_
    if(iprisubmat>=2) write(nfout,'(" sw_scalapack = ",i3," <<evolve_WFs_in_subspace>>")') sw_scalapack
    if(sw_scalapack == ON) then
       if(submat_uncalled) then
! === DEBUG by tkato 2012/01/23 ================================================
          if(block_size == 0) then
             max_block_size = int(real(meg)/real(nrank_e))
             block_size =  nb_mgs_default
             if(block_size.ge.max_block_size) block_size = max_block_size
          end if
! ==============================================================================
          call set_nprow_npcol(nprow,npcol)
#ifdef SINGLE_CONTEXT
          if(allocated(usermap)) deallocate(usermap)
          allocate(usermap(nprow,npcol))
#else
          if(allocated(ictxt)) deallocate(ictxt)
          if(allocated(usermap)) deallocate(usermap)
          allocate(ictxt(0:nrank_k-1))
          allocate(usermap(nprow,npcol,0:nrank_k-1))
#endif
!         call scalapack_setup(meg,lwork,lrwork,liwork,lda,occ,ictxt,myrow,mycol,desca,descz,usermap)
#ifdef SINGLE_CONTEXT
          call scalapack_setup(meg,lwork1,lrwork1,liwork1,lwork2,lrwork2,liwork2,lda,occ,ictxt,myrow,mycol,desca,descz,usermap)
#else
          call scalapack_setup(meg,lwork1,lrwork1,liwork1,lwork2,lrwork2,liwork2,lda,occ,ictxt(myrank_k),&
           & myrow,mycol,desca,descz,usermap)
#endif
! ===== We need more work area!!! ==============================================
          lwork1  = 10*lwork1
          lrwork1 = 10*lrwork1
          liwork1 = 10*liwork1
          lwork2  = 10*lwork2
          lrwork2 = 10*lrwork2
          liwork2 = 10*liwork2
! ==============================================================================
       end if
       allocate(amat(lda*kimg_t,occ))
       allocate(zmat(lda*kimg_t,occ))
       if(sw_rsb==ON) then
          if(.not.allocated(unimat)) then
             allocate(unimat(lda*kimg_t,occ))
             allocate(unimat_h(lda*kimg_t,occ))
          endif
       endif
    else
       allocate(zmat(meg*kimg_t,meg))
       if(sw_rsb==ON) then
          if(.not.allocated(unimat)) then
             allocate(unimat(meg*kimg_t,meg))
             allocate(unimat_h(meg*kimg_t,meg))
          endif
       endif
    end if

#else
    allocate(zmat(meg*kimg_t,meg))
    if(sw_rsb==ON) then
       if(.not.allocated(unimat)) then
          allocate(unimat(meg*kimg_t,meg))
          allocate(unimat_h(meg*kimg_t,meg))
       endif
    endif
    zmat=0.0d0
#endif
    if(submat_uncalled) then
!!$       call set_col_partition ! -> npart,isp,iep, npart2, isp2, iep2
       call set_col_partition_G(kimg_t,meg,ibsize) ! -> npart,isp,iep, npart2, isp2, iep2
       submat_uncalled = .false.
    end if

    call tstatc0_begin('evolve_WFs_in_subspace(PART1) ', id_sname2)
!!$    if(iprisubmat>=2) then
!!$       write(nfout,'(" --- neordr ---")')
!!$       write(nfout,'(" ik = ",i8)') ik
!!$       write(nfout,'(10i8)') (neordr(ib1,ik),ib1=1, neg)
!!$       write(nfout,'(" --- nrvf_ordr ---")')
!!$       write(nfout,'(" ik = ",i8)') ik
!!$       write(nfout,'(10i8)') (nrvf_ordr(ib1,ik),ib1=1, neg)
!!$    end if

    ekod = 0.d0
    sum_abs_diagonal = 0.d0
    sum_sq_diagonal = 0.d0
    sum_sq_non_diagonal = 0.d0
    sum_abs_non_diagonal = 0.d0
    PART1: do jp = 1,npart
       is = isp(jp)
       ie = iep(jp)
       msize = is*(ie-is+1) + (ie-is+1)*(ie-is)/2
       allocate(w1hw2(msize*kimg_t))
       w1hw2 = 0.d0
       if(iprisubmat>=2) then
          write(nfout,'(" jp, is, ie = ",3i8)') jp, is, ie
       end if
#ifdef SUBMAT_DGEMM
! <-- T.Kokubo & D.Fukata Feb. 2010
       nblk1 = ie-is+1
#ifdef SX
       ibsize2 = min(ie,meg)
#else
       ibsize2 = ibsize
#endif

       allocate(w1hw2r(ibsize2,nblk1)); w1hw2r = 0.d0
       if(kimg==2 .and. k_symmetry(ik) /= GAMMA) then
          allocate(w1hw2i(ibsize2,nblk1)); w1hw2i = 0.d0
       end if

       do ibsta=1,min(ie,meg), ibsize2
          ibend=min( ibsta+ibsize2-1, min(ie,meg) )
          nblk2=ibend-ibsta+1
          w1hw2r = 0.d0
          if(kimg==2 .and. k_symmetry(ik) /= GAMMA) w1hw2i = 0.d0

          if(kimg==1) then
             call dgemm('T','N',nblk2,nblk1,np_g1k(ik),1.d0,zat_l(1,ibsta,1),np_g1k_x, &
                  &     zah_l(1,is,1),np_g1k_x,1.d0,w1hw2r,ibsize2)
             do ib1 = max(is,ibsta), min(ie,meg)
                do ib2 = ibsta,min(ibend,ib1)
                   m = ib1*(ib1-1)/2-is*(is-1)/2+ib2
                   w1hw2(m)=w1hw2r(ib2-ibsta+1, ib1-max(is,ibsta)+1)
                end do
             end do
          else
             if(k_symmetry(ik) == GAMMA) then
                ig1=1;if(myrank_e==0) ig1=2
                call dgemm('T','N',nblk2,nblk1,np_g1k(ik),1.d0,zat_l(1,ibsta,1),np_g1k_x, &
                     &     zah_l(1,is,1),np_g1k_x,1.d0,w1hw2r,ibsize2)
                call dgemm('T','N',nblk2,nblk1,np_g1k(ik),1.d0,zat_l(1,ibsta,2),np_g1k_x, &
                     &     zah_l(1,is,2),np_g1k_x,1.d0,w1hw2r,ibsize2)
                w1hw2r = w1hw2r*2.d0
                if(myrank_e == 0) then
                   do ib1 = max(is,ibsta), min(ie,meg)
                      do ib2 = ibsta,min(ibend,ib1)
                         w1hw2r(ib2-ibsta+1, ib1-max(is,ibsta)+1) =  w1hw2r(ib2-ibsta+1, ib1-max(is,ibsta)+1) &
                        &                                          - 2.d0*zah_l(1,ib1,2)*zat_l(1,ib2,2)       &
                        &                                          -      zah_l(1,ib1,1)*zat_l(1,ib2,1)
                      end do
                   end do
                end if
                do ib1 = max(is,ibsta), min(ie,meg)
                   do ib2 = ibsta,min(ibend,ib1)
                      m = ib1*(ib1-1)/2-is*(is-1)/2+ib2
                      w1hw2(m)=w1hw2r(ib2-ibsta+1, ib1-max(is,ibsta)+1)
                   end do
                end do
             else
                call dgemm('T','N',nblk2,nblk1,np_g1k(ik), 1.d0,zat_l(1,ibsta,1),np_g1k_x, &
                     &     zah_l(1,is,1),np_g1k_x,1.d0,w1hw2r,ibsize2)
                call dgemm('T','N',nblk2,nblk1,np_g1k(ik), 1.d0,zat_l(1,ibsta,2),np_g1k_x, &
                     &     zah_l(1,is,2),np_g1k_x,1.d0,w1hw2r,ibsize2)
                call dgemm('T','N',nblk2,nblk1,np_g1k(ik), 1.d0,zat_l(1,ibsta,2),np_g1k_x, &
                     &     zah_l(1,is,1),np_g1k_x,1.d0,w1hw2i,ibsize2)
                call dgemm('T','N',nblk2,nblk1,np_g1k(ik),-1.d0,zat_l(1,ibsta,1),np_g1k_x, &
                     &     zah_l(1,is,2),np_g1k_x,1.d0,w1hw2i,ibsize2)
                do ib1 = max(is,ibsta), min(ie,meg)
                   do ib2 = ibsta,min(ibend,ib1)
                      m = ib1*(ib1-1)/2-is*(is-1)/2+ib2
                      w1hw2(2*m-1) = w1hw2r(ib2-ibsta+1, ib1-max(is,ibsta)+1)
                      w1hw2(2*m  ) = w1hw2i(ib2-ibsta+1, ib1-max(is,ibsta)+1)
                   end do
                end do
             end if
          end if
       enddo
! -->
       if(kimg==2 .and. k_symmetry(ik) /= GAMMA) deallocate(w1hw2i)
       deallocate(w1hw2r)

#else
! parallel loop
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(ib2,ii,dr1,dr2,di1,di2,m)
!CDIR CONCUR(BY=1)
#endif
! --> T. Yamasaki 28th Aug. 2009
!   ib1to --> ib1,  instead of ib1to = neordr(ib1,ik)
!   ib2to --> ib2
       do ib1 = is, ie
          if(ib1>meg) cycle
          do ib2 = 1, ib1
             m = ib1*(ib1-1)/2-is*(is-1)/2+ib2
             if(kimg == 1) then
                do ii = 1, np_g1k(ik)              ! MPI
                   dr1 = zah_l(ii,ib1,1)
                   dr2 = zat_l(ii,ib2,1)
                   w1hw2(m) = w1hw2(m) + dr1*dr2
                end do
             else
                if(k_symmetry(ik) == GAMMA) then
                   dd = 0.d0
!!$                   ig1 = 1;  if(mype == 0) ig1 = 2
                   ig1 = 1;  if(myrank_e == 0) ig1 = 2
                   do ii = ig1, np_g1k(ik)            ! MPI
                      dr1 = zah_l(ii,ib1,1)
                      di1 = zah_l(ii,ib1,2)
                      dr2 = zat_l(ii,ib2,1)
                      di2 = zat_l(ii,ib2,2)
                      dd = dd + (dr1*dr2+di1*di2)*2.d0
                   end do
!!$                   if(mype == 0) dd = dd + zah_l(1,ib1,1)*zat_l(1,ib2,1)
                   if(myrank_e == 0) dd = dd + zah_l(1,ib1,1)*zat_l(1,ib2,1)
                   w1hw2(m)=dd
                else
                   do ii = 1, np_g1k(ik)             ! MPI
                      dr1 = zah_l(ii,ib1,1)
                      di1 = zah_l(ii,ib1,2)
                      dr2 = zat_l(ii,ib2,1)
                      di2 = zat_l(ii,ib2,2)
                      w1hw2(2*m-1)=w1hw2(2*m-1)+dr1*dr2+di1*di2
                      w1hw2(2*m  )=w1hw2(2*m  )+dr1*di2-di1*dr2
                   end do
                end if
             end if
          end do
       end do
#endif

       if(iprisubmat >= 2) call wd_w1hw2(" -- just after making w1hw2 --")

!! (off-diagonal matrix elements will be damped.)
       if(kimg_t == 1) then
          do ib1=max(2,is),ie
             do ib2=1,ib1-1
                m = ib1*(ib1-1)/2-is*(is-1)/2+ib2
                sum_sq_non_diagonal = sum_sq_non_diagonal + w1hw2(m)*w1hw2(m)
                sum_abs_non_diagonal = sum_abs_non_diagonal + abs(w1hw2(m))
                w1hw2(m) = damp*w1hw2(m)
             enddo
          enddo
          do ib = is, ie
             m = ib*(ib-1)/2-is*(is-1)/2+ib
             sum_sq_diagonal = sum_sq_diagonal + w1hw2(m)**2
             sum_abs_diagonal = sum_abs_diagonal + abs(w1hw2(m))
             ekod = ekod + w1hw2(m)
          end do
       else
          do ib1=max(2,is),ie
             do ib2=1,ib1-1
                m = ib1*(ib1-1)/2-is*(is-1)/2+ib2
                dd = w1hw2(2*m-1)*w1hw2(2*m-1) + w1hw2(2*m)*w1hw2(2*m)
                sum_sq_non_diagonal = sum_sq_non_diagonal + dd
                if(dd >= SmallestPositiveNumber) sum_abs_non_diagonal = sum_abs_non_diagonal + sqrt(dd)
                w1hw2(2*m-1) = damp*w1hw2(2*m-1)
                w1hw2(2*m  ) = damp*w1hw2(2*m  )

             enddo
          enddo
          do ib = is, ie
             m = ib*(ib-1)/2-is*(is-1)/2+ib
             sum_sq_diagonal = sum_sq_diagonal + w1hw2(2*m-1)**2
             sum_abs_diagonal = sum_abs_diagonal + abs(w1hw2(2*m-1))
             ekod = ekod + w1hw2(2*m-1)
          end do
       endif
       if(iprisubmat >= 2) call wd_w1hw2(" -- after damping --")

#ifdef _USE_SCALAPACK_
       if(sw_scalapack==ON) then
          call set_amat(lda,kimg_t,occ,amat,meg,nprow,npcol,block_size,usermap,is,ie,w1hw2)
       else
#endif
          call set_hmat(kimg_t,meg,is,ie,w1hw2,zmat)
#ifdef _USE_SCALAPACK_
       end if
#endif

       deallocate(w1hw2)
    end do PART1
    call tstatc0_end(id_sname2)
    deallocate(zah_l)

    sum_abs_diagonal = sum_abs_diagonal/meg
    sum_abs_non_diagonal = sum_abs_non_diagonal/(meg*(meg-1)/2)

    if(iprisubmat >= 2) then
       write(nfout,*) 'neordr for ik = ',ik
       write(nfout,9002) (neordr(ib1,ik),ib1=1,neg)
       write(nfout,*) 'nrvf_ordr for ik = ',ik
       write(nfout,9002) (nrvf_ordr(ib1,ik),ib1=1,neg)
       write(nfout,*) 'eig'
       write(nfout,'(5x,10f8.4)') (eko_l(ib1,ik),ib1=1,neg)
    endif
9002 format(5x,10i8)

!! (Diagonalization )  !!
#ifdef _USE_SCALAPACK_
     if(sw_scalapack == ON) then
        if(kimg_t == 1) then
!          call pdsyev_driver(meg,eig,desca,descz,lda,occ,amat,zmat,lwork,liwork)
           call pdsyev_driver(meg,eig,desca,descz,lda,occ,amat,zmat,lwork1,liwork1)
        else
!          call pzheev_driver(meg,eig,desca,descz,lda,occ,amat,zmat,lwork,lrwork,liwork)
           call pzheev_driver(meg,eig,desca,descz,lda,occ,amat,zmat,lwork2,lrwork2,liwork2)
        endif
     else
#endif
        if(kimg_t == 1) then
           call dsyev_driver(meg,eig,zmat)
        else
           call zheev_driver(meg,eig,zmat)
        endif
#ifdef _USE_SCALAPACK_
     endif
#endif

     if(iprisubmat >= 2) then
        write(nfout,*) 'eko_d for ik = ',ik
        write(nfout,9001) (eko_d(ib),ib=1,meg)
        write(nfout,*) 'eig for ik = ',ik
        write(nfout,9001) (eig(ib),ib=1,meg)
        !!$call wd_w1hw2(" -- after diagonalization --")
        !!$call wd_zmat(" -- after diagonalization --")
!sum eko
     end if

     call get_ipri0(iprisubmat,ipri0)
     if(ipri0>=2) then
        dr1=0.d0;dr2=0.d0
        do ib1=1,meg
           ib1to = neordr(ib1,ik)
           if(map_e(ib1to) == myrank_e) dr1=dr1+eko_l(map_z(ib1to),ik)
           dr2=dr2+eig(ib1)
        enddo
        call mpi_allreduce(MPI_IN_PLACE, dr1,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
     end if
     if(iprisubmat>=2) then
        write(nfout,*) ' sum of eko_l & eig =',dr1,dr2
        write(nfout,'(" sum_sq_diagonal , sum_sq_non_diagonal  = ",2d20.12)') sum_sq_diagonal, sum_sq_non_diagonal
        write(nfout,'(" sum_abs_diagonal ,sum_abs_non_diagonal = ",2d20.12)') sum_abs_diagonal, sum_abs_non_diagonal
     end if

!! (subspace rotation) !!
     allocate(zat_wk(np_g1k_x,meg,kimg))
     allocate(zahloc_wk(np_g1k_x,meg,kimg))
     zat_wk=0.d0
     zahloc_wk=0.d0
! <-- T.Kokubo & D.Fukata, Feb. 2010
        PART2: do jp = 1,npart2
           is = isp2(jp)
           ie = iep2(jp)
           nblk1 = ie-is+1
#ifdef SX
           ibsize2 = meg
#else
           ibsize2 = ibsize
#endif
#ifdef _USE_SCALAPACK_
     if(sw_scalapack == ON) then
           allocate(zz(meg*kimg_t,nblk1))
           call put_zmat(kimg_t,lda,occ,zmat,meg,nprow,npcol,block_size,usermap,is,ie,zz)
           if(kimg == 1) then
              call subspace_rotation_real(is,ie,zat_wk,zz,nblk1,ibsize2,ik,meg,zat_l,zahloc_wk, zahloc_l)
           else
              call subspace_rotation_imag(is,ie,kimg_t,zat_wk,zz,nblk1,ibsize2,ik,meg,zat_l,zahloc_wk, zahloc_l)
           endif
           deallocate(zz)
     else
#endif
           if(kimg == 1) then
              call subspace_rotation_real(is,ie,zat_wk,zmat(1,is),nblk1,ibsize2,ik,meg,zat_l,zahloc_wk, zahloc_l)
           else
              call subspace_rotation_imag(is,ie,kimg_t,zat_wk,zmat(1,is),nblk1,ibsize2,ik,meg,zat_l,zahloc_wk, zahloc_l)
           endif
#ifdef _USE_SCALAPACK_
     end if
#endif
        end do PART2

        if(sw_rsb==ON) then
           unimat(:,:) = zmat(:,:)
           do ib1=1,neg
              do ib2=1,neg
                 unimat_h(kimg_t*ib2-1,ib1) =  unimat(kimg_t*ib1-1,ib2)
                 unimat_h(kimg_t*ib2,  ib1) = -unimat(kimg_t*ib1,  ib2)
              enddo
           enddo
        endif
!
!#ifdef _USE_SCALAPACK_
!     if(sw_scalapack == ON) then
!        PART2: do jp = 1,npart2
!           is = isp2(jp)
!           ie = iep2(jp)
!           allocate(zz(meg*kimg_t,is:ie))
!           call put_zmat(lda,occ,zmat,meg,nprow,npcol,block_size,usermap,is,ie,zz)
!           if(kimg == 1) then
!              call subspace_rotation_real(is,ie,zat_wk,zz)
!           else
!              call subspace_rotation_imag(is,ie,zat_wk,zz)
!           endif
!           deallocate(zz)
!        end do PART2
!     else
!#endif
!        is=1; ie=meg
!        if(kimg == 1) then
!           call subspace_rotation_real(is,ie,zat_wk,zmat)
!        else
!           call subspace_rotation_imag(is,ie,zat_wk,zmat)
!        endif
!#ifdef _USE_SCALAPACK_
!     end if
!#endif
! -->
     zat_l(1:np_g1k(ik),1:meg,1:kimg) = zat_wk(1:np_g1k(ik),1:meg,1:kimg)
     if(sw_keep_hloc_phi==ON) zahloc_l(1:np_g1k(ik),1:meg,1:kimg) = zahloc_wk(1:np_g1k(ik),1:meg,1:kimg)
     deallocate(zat_wk)
     deallocate(zahloc_wk)

!! (zaj_l) by tenchi
     call m_ES_W_transpose_back_r(.false.,ista_k,iend_k,ik,zaj_l,zat_l)            ! zat_l -> zaj_l
     if(sw_keep_hloc_phi==ON) &
     call m_ES_W_transpose_back_r(.false.,ista_k,iend_k,ik,hlocphi_l,zahloc_l)            ! zat_l -> zaj_l
!!$     call m_ES_W_transpose_back2(ista_k,iend_k,ik,0,neordr,zaj_l,zat_l)
     if(meg < neg) then
        call m_ES_W_transpose_back2(ista_k,iend_k,ik,meg,neordr,zaj_l,zat_l)  ! zat_l -> zaj_l
        if(sw_keep_hloc_phi==ON) &
        call m_ES_W_transpose_back2(ista_k,iend_k,ik,meg,neordr,hlocphi_l,zahloc_l)  ! zat_l -> zaj_l
     end if
#ifdef SAVE_FFT_TIMES
     if(sw_save_fft == ON) then
        do ib = 1, np_e
           status_saved_phifftr(ib,ik) = OLD
        end do
     end if
#endif

!! (eko_l)
     if(meg < neg) then
        do ib1=meg+1,neg
           if(map_e(ib1) == myrank_e) then        ! MPI
              eko_l(map_z(ib1),ik)=eko_d(neordr(ib1,ik))
           end if
        enddo
     endif
     do ib1 = 1, meg
        if(map_e(ib1) == myrank_e) then         ! MPI
           eko_l(map_z(ib1),ik)=eig(ib1)
        end if
     end do
     if(iprisubmat >= 2) then
        eko2 = sum(eig(1:meg))
        write(nfout,1201) ik,eko1,ekod,eko2

        write(nfout,*) 'eko_l'
        write(nfout,9001) (eko_l(ib1,ik),ib1=1,meg)
     endif

1201 format(' %% for ik = ',i4,4x,' eko1&ekod&eko2 = ',3f14.7)
9001 format(5x,6f12.5)
!! (neordr & nrvf_ordr)

     neordr(1:neg,ik) = (/(ib1,ib1=1,neg)/)
     nrvf_ordr(1:neg,ik) = (/(ib1,ib1=1,neg)/)

! (deallocate)
    deallocate(eko_d)
    deallocate(eig)
    deallocate(zat_l)
    !!$deallocate(zah_l)
    if(sw_scalapack==ON) deallocate(amat)
    deallocate(zmat)

! === DEBUG by tkato 2012/12/11 ================================================
    call mpi_allreduce(MPI_IN_PLACE,sum_abs_non_diagonal,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
    call mpi_allreduce(MPI_IN_PLACE,sum_abs_diagonal,    1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
! ==============================================================================
    if(.not.non_diagonal_part_is_small(ik)  &
         & .and. sum_abs_non_diagonal/sum_abs_diagonal <= submat_critical_ratio) &
         & non_diagonal_part_is_small(ik) = .true.

    if(iprisubmat >= 2) write(nfout,'(" !submat  non_diagonal/diagonal = ",d20.8)') &
         & sum_abs_non_diagonal/sum_abs_diagonal

    call tstatc0_end(id_sname)

  contains
    subroutine get_ipri0(ipri_in, ipri_out)
      integer, intent(in)  :: ipri_in
      integer, intent(out) :: ipri_out
      if(nrank_e>1) then
         if(myrank_e == 0) ipri_out = ipri_in
         call mpi_bcast(ipri_out,1,mpi_integer,0,mpi_k_world(myrank_k),ierr)
      else
         ipri_out = ipri_in
      end if
    end subroutine get_ipri0

    subroutine wd_w1hw2(somecomment)
      character(len=*), intent(in) :: somecomment
      integer :: ib1, ib2, neg_wd1, neg_wd2
      write(nfout,'(a35)') somecomment
      write(nfout,*) 'w1hw2 for ik = ',ik
      neg_wd1 = 8
      neg_wd2 = 8
      if(neg_wd1 > ie) neg_wd1 = ie
      if(neg_wd2 > meg) neg_wd2 = meg
      write(nfout,'(" real part (",i8,":",i8,",1:",i8,")")') is,neg_wd1,neg_wd2
      if(kimg_t == 1) then
         do ib1=is,neg_wd1
            write(nfout,9001) (w1hw2(ib1*(ib1-1)/2-is*(is-1)/2+ib2),ib2=1,neg_wd2)
         enddo
         write(nfout,'(" real part (",i8,":",i8,",1:",i8,")")') ie-neg_wd1+1,ie, neg_wd2
         do ib1=ie-neg_wd1+1, ie
            write(nfout,9001) (w1hw2(ib1*(ib1-1)/2-is*(is-1)/2+ib2),ib2=1,neg_wd2)
         enddo
      else
         do ib1=is,neg_wd1
            write(nfout,9001) (w1hw2(2*(ib1*(ib1-1)/2-is*(is-1)/2+ib2)-1),ib2=1,neg_wd2)
         enddo
         write(nfout,'(" real part (",i8,":",i8,",1:",i8,")")') ie-neg_wd1+1,ie, neg_wd2
         do ib1=ie-neg_wd1+1, ie
            write(nfout,9001) (w1hw2(2*(ib1*(ib1-1)/2-is*(is-1)/2+ib2)-1),ib2=1,neg_wd2)
         enddo
         write(nfout,'(" imag part")')
         do ib1=is,neg_wd1
            write(nfout,9001) (w1hw2(2*(ib1*(ib1-1)/2-is*(is-1)/2+ib2)),ib2=1,neg_wd2)
         enddo
      end if
9001  format(1x,8e14.6)
      write(nfout,*) 'eko_l for ik = ',ik
      write(nfout,9001) (eko_d(neordr(ib1,ik)),ib1=1,neg_wd2)
    end subroutine wd_w1hw2

#ifdef _USE_SCALAPACK_
#endif
! _USE_SCALAPACK_

  end subroutine evolve_WFs_in_subspace
!fj$$#endif

    subroutine subspace_rotation_real(is,ie,zat_wk,zz,ibsize1,ibsize2,ik,meg,zat_l,zahloc_wk,zahloc_l)
!    Revised by T. Yamasaki 1st Sep. 2009
!    modified by T.Kokubo & D.Fukata, Feb. 2010
      integer, intent(in) :: is,ie, ibsize1,ibsize2,ik,meg
      real(kind=DP), intent(inout) :: zat_wk(np_g1k_x,meg,kimg),zahloc_wk(np_g1k_x,meg,kimg)
      real(kind=DP), intent(in) :: zat_l(np_g1k_x,meg,kimg), zahloc_l(np_g1k_x,meg,kimg)
      real(kind=DP), intent(in) :: zz(meg,ibsize1)
      integer :: ibsta, ibend
      integer :: ib1,ib2
#ifdef SUBMAT_DGEMM
      integer :: nblk2
      real(kind=DP) :: alpha, beta
#endif
      integer :: id_sname = -1

      call tstatc0_begin('subspace_roation_real ', id_sname)
#ifdef SUBMAT_DGEMM
      alpha=1.d0; beta=1.d0
      do ibsta=1,meg, ibsize2
         ibend=min(ibsta+ibsize2-1,meg)
         nblk2=ibend-ibsta+1
            call dgemm('N','N',np_g1k(ik),ibsize1,nblk2,alpha,zat_l(1,ibsta,1),np_g1k_x,zz(ibsta,1),meg, &
           &           beta, zat_wk(1,is,1),np_g1k_x)
            if(sw_keep_hloc_phi==ON) &
            call dgemm('N','N',np_g1k(ik),ibsize1,nblk2,alpha,zahloc_l(1,ibsta,1),np_g1k_x,zz(ibsta,1),meg, &
           &           beta, zahloc_wk(1,is,1),np_g1k_x)
      enddo
#else
      do ibsta=1,meg, ibsize2
         ibend=min(ibsta+ibsize2-1,meg)
         do ib2=is,ie
            do ib1=ibsta,ibend
               do ii=1,np_g1k(ik)
                  zat_wk(ii,ib2,1)=zat_wk(ii,ib2,1)+zz(ib1,ib2-is+1)*zat_l(ii,ib1,1)
               enddo
            enddo
         enddo
      enddo
#endif
      call tstatc0_end(id_sname)
    end subroutine subspace_rotation_real

    subroutine subspace_rotation_imag(is,ie,kimg_t,zat_wk,zz,ibsize1,ibsize2,ik,meg,zat_l,zahloc_wk,zahloc_l)
!    Revised by T. Yamasaki 1st Sep. 2009
!    modified by T.Kokubo & D.Fukata, Feb. 2010
      integer, intent(in) :: is,ie,kimg_t, ibsize1,ibsize2,ik,meg
      real(kind=DP), intent(inout) :: zat_wk(np_g1k_x,meg,kimg),zahloc_wk(np_g1k_x,meg,kimg)
      real(kind=DP), intent(in) :: zat_l(np_g1k_x,meg,kimg),zahloc_l(np_g1k_x,meg,kimg)
      real(kind=DP), intent(in) :: zz(meg*kimg_t,ibsize1)
      integer :: id_sname = -1
      integer :: ibsta, ibend
      integer :: ib1,ib2,ii
#ifdef SUBMAT_DGEMM
      real(kind=DP) :: alpha, beta
      integer :: nblk2
      real(kind=DP), allocatable, dimension(:,:) :: zzr,zzi
#endif

      call tstatc0_begin('subspace_roation_imag ', id_sname)

      if(k_symmetry(ik) == GAMMA) then
#ifdef SUBMAT_DGEMM
         alpha = 1.d0; beta = 1.d0
         do ibsta=1,meg,ibsize2
            ibend=min(ibsta+ibsize2-1,meg)
            nblk2=ibend-ibsta+1
            call dgemm('N','N',np_g1k(ik),ibsize1,nblk2,alpha,zat_l(1,ibsta,1),np_g1k_x,zz(ibsta,1),meg, beta &
                 &                         , zat_wk(1,is,1),np_g1k_x)
            call dgemm('N','N',np_g1k(ik),ibsize1,nblk2,alpha,zat_l(1,ibsta,2),np_g1k_x,zz(ibsta,1),meg, beta &
                 &                         , zat_wk(1,is,2),np_g1k_x)
            if(sw_keep_hloc_phi==ON) then
              call dgemm('N','N',np_g1k(ik),ibsize1,nblk2,alpha,zahloc_l(1,ibsta,1),np_g1k_x,zz(ibsta,1),meg, beta &
                   &                         , zahloc_wk(1,is,1),np_g1k_x)
              call dgemm('N','N',np_g1k(ik),ibsize1,nblk2,alpha,zahloc_l(1,ibsta,2),np_g1k_x,zz(ibsta,1),meg, beta &
                   &                         , zahloc_wk(1,is,2),np_g1k_x)
            endif
         enddo
         !!if(iprisubmat>=2) write(nfout,'(" !!<<subspace_rotation_imag>> ibsize = ",i8)') ibsize
#else
         do ibsta=1,meg,ibsize2
            ibend=min(ibsta+ibsize2-1,meg)
            do ib2=is,ie
!CDIR OUTERUNROLL=4
               do ib1=ibsta,ibend
                  dr1=zz(ib1,ib2-is+1)
                  do ii=1,np_g1k(ik)
                     zat_wk(ii,ib2,1)=zat_wk(ii,ib2,1)+dr1*zat_l(ii,ib1,1)
                     zat_wk(ii,ib2,2)=zat_wk(ii,ib2,2)+dr1*zat_l(ii,ib1,2)
                  enddo
               enddo! ib1 loop
            enddo! ib2 loop
         enddo
#endif
         if(iprisubmat >= 2) then
            do ib2 = is, ie
               write(nfout,'(" --- zat_l (input, output) ---, ib =",i5)') ib2
               write(nfout,'(" zat_l(input ) Re : ",5e14.6)') (zat_l(ii,ib2,1),ii=1,5)
               write(nfout,'(" zat_l(input ) Im : ",5e14.6)') (zat_l(ii,ib2,2),ii=1,5)
               write(nfout,'(" zat_l(output) Re : ",5e14.6)') (zat_wk(ii,ib2,1),ii=1,5)
               write(nfout,'(" zat_l(output) Im : ",5e14.6)') (zat_wk(ii,ib2,2),ii=1,5)
            end do
         end if
      else
#ifdef SUBMAT_DGEMM
         allocate(zzr(ibsize2,ibsize1)); allocate(zzi(ibsize2,ibsize1))
         do ibsta=1,meg,ibsize2
            ibend=min(ibsta+ibsize2-1,meg)
            nblk2=ibend-ibsta+1
            do ib2=is,ie
               do ib1=ibsta, ibend
                  zzr(ib1-ibsta+1,ib2-is+1)=zz(2*ib1-1,ib2-is+1)
                  zzi(ib1-ibsta+1,ib2-is+1)=zz(2*ib1  ,ib2-is+1)
               end do
            end do

            alpha=1.d0; beta=1.d0
            call dgemm('N','N',np_g1k(ik),ibsize1,nblk2,alpha,zat_l(1,ibsta,1),np_g1k_x,zzr,ibsize2 &
                 &     ,beta, zat_wk(1,is,1),np_g1k_x )
            if(sw_keep_hloc_phi==ON) &
            call dgemm('N','N',np_g1k(ik),ibsize1,nblk2,alpha,zahloc_l(1,ibsta,1),np_g1k_x,zzr,ibsize2 &
                 &     ,beta, zahloc_wk(1,is,1),np_g1k_x )
            alpha=-1.d0; beta=1.d0
            call dgemm('N','N',np_g1k(ik),ibsize1,nblk2,alpha,zat_l(1,ibsta,2),np_g1k_x,zzi,ibsize2 &
                 &     ,beta, zat_wk(1,is,1),np_g1k_x )
            if(sw_keep_hloc_phi==ON) &
            call dgemm('N','N',np_g1k(ik),ibsize1,nblk2,alpha,zahloc_l(1,ibsta,2),np_g1k_x,zzi,ibsize2 &
                 &     ,beta, zahloc_wk(1,is,1),np_g1k_x )
            alpha=1.d0; beta=1.d0
            call dgemm('N','N',np_g1k(ik),ibsize1,nblk2,alpha,zat_l(1,ibsta,2),np_g1k_x,zzr,ibsize2 &
                 &     ,beta, zat_wk(1,is,2),np_g1k_x )
            if(sw_keep_hloc_phi==ON) &
            call dgemm('N','N',np_g1k(ik),ibsize1,nblk2,alpha,zahloc_l(1,ibsta,2),np_g1k_x,zzr,ibsize2 &
                 &     ,beta, zahloc_wk(1,is,2),np_g1k_x )
            alpha=1.d0; beta=1.d0
            call dgemm('N','N',np_g1k(ik),ibsize1,nblk2,alpha,zat_l(1,ibsta,1),np_g1k_x,zzi,ibsize2 &
                 &     ,beta, zat_wk(1,is,2),np_g1k_x )
            if(sw_keep_hloc_phi==ON) &
            call dgemm('N','N',np_g1k(ik),ibsize1,nblk2,alpha,zahloc_l(1,ibsta,1),np_g1k_x,zzi,ibsize2 &
                 &     ,beta, zahloc_wk(1,is,2),np_g1k_x )
         enddo
         deallocate(zzi,zzr)
#else
         do ibsta=1,meg,ibsize2
            ibend=min(ibsta+ibsize2-1,meg)
            do ib2=is,ie
!CDIR OUTERUNROLL=4
               do ib1=ibsta,ibend
                  dr1=zz(2*ib1-1,ib2-is+1)
                  di1=zz(2*ib1  ,ib2-is+1)
                  do ii=1,np_g1k(ik)
                     dr2=zat_l(ii,ib1,1)
                     di2=zat_l(ii,ib1,2)
                     zat_wk(ii,ib2,1)=zat_wk(ii,ib2,1)+dr1*dr2-di1*di2
                     zat_wk(ii,ib2,2)=zat_wk(ii,ib2,2)+dr1*di2+di1*dr2
                  enddo
              enddo! ib1 loop
            enddo! ib2 loop
         enddo
#endif
       end if
       call tstatc0_end(id_sname)
    end subroutine subspace_rotation_imag

  subroutine set_col_partition_G(kimg_t,meg,ibsize)
    integer, intent(in) :: kimg_t,meg
    integer, intent(inout) :: ibsize
!!$integer,parameter :: maxmatsize = 67108864 ! 512MB
    integer :: maxmatsize
    integer :: trisize, size
    integer :: rowsize, block
    integer :: i,is,ie
    integer :: maxsize
! <-- T.Kokubo & D.Fukata, Feb. 2010
#ifdef SUBMAT_DGEMM
    integer :: ibsize_tmp, modnum
    integer,allocatable  :: tmpmem(:)
#endif
! -->

    integer :: id_sname = -1
    call tstatc0_begin('set_col_partition_G ', id_sname)

    maxmatsize = msize_submat * 1024 * 1024 / 8 ! bytes

! <-- T.Kokubo & D.Fukata, Feb. 2010
#ifdef SUBMAT_DGEMM
    ibsize_tmp = ibsize
    npart = meg/ibsize_tmp
    modnum=mod(meg,ibsize_tmp)
    if( modnum.ne.0 ) npart = npart + 1
    if(iprisubmat >= 2) write(nfout,*)' !! SUBMAT :: ibsize=', ibsize, npart

    if(allocated(isp)) deallocate(isp)
    if(allocated(iep)) deallocate(iep)
    if(allocated(tmpmem)) deallocate(tmpmem)
    do
       allocate( isp(npart)   ) ; isp=1
       allocate( iep(npart)   ) ; iep=1
       allocate( tmpmem(npart)) ; tmpmem=0
       do i=1,npart
          is = ibsize_tmp*(i-1)+1           ; isp(i)=is
          ie = min(isp(i)+ibsize_tmp-1,meg) ; iep(i)=ie
#ifdef SX
          tmpmem(i) = kimg*(is*(ie-is+1)+(ie-is+1)*(ie-is)/2) +  2*ie*(ie-is+1)
          !           <----------      w1hw2      ----------> <- w1hw2r,w1hw2i ->
#else
          tmpmem(i) = kimg*(is*(ie-is+1)+(ie-is+1)*(ie-is)/2) + 2*ibsize_tmp*(ie-is+1)
          !           <----------      w1hw2      ----------> <--- w1hw2r, w1hw2i --->
#endif
       enddo
       maxsize=maxval(tmpmem)
       if( maxsize <= maxmatsize ) exit

       ibsize_tmp = ibsize_tmp - 1
       npart = meg/ibsize_tmp
       modnum=mod(meg,ibsize_tmp)
       if( modnum.ne.0 ) npart = npart + 1
       deallocate(isp,iep,tmpmem)
    enddo
    deallocate(tmpmem)
    if( ibsize .ne. ibsize_tmp ) then
       write(nfout,*)' !! SUBMAT_DGEMM : ibsize is changed', ibsize, ' --> ', ibsize_tmp
       ibsize = ibsize_tmp
    endif
#else
    trisize = meg*(meg+1)/2*kimg_t
    if(trisize <= maxmatsize.or.sw_scalapack==OFF) then
       npart = 1
       allocate(isp(npart))
       allocate(iep(npart))
       isp(npart) = 1
       iep(npart) = meg
       maxsize  = trisize
    else
       npart = trisize/maxmatsize
       allocate(isp(npart))
       allocate(iep(npart))
       maxsize  = trisize/npart
       ie = 0
       do i=1,npart
          is = ie+1
          ie = is
          size = 0
          do while(size*kimg_t <= maxsize)
             ie = ie+1
             size = is*(ie-is+1)+(ie-is+2)*(ie-is+1)/2
          end do
          isp(i) = is
          iep(i) = ie
          if(ie>meg) then
             npart = i
             exit
          end if
       end do
       iep(npart) = meg
    end if
!!$#endif
#endif

#ifdef SUBMAT_DGEMM
    if(allocated(isp2)) deallocate(isp2)
    if(allocated(iep2)) deallocate(iep2)
    npart2 = npart
    allocate(isp2(npart2))
    allocate(iep2(npart2))
    isp2 = isp
    iep2 = iep
#else
    if(allocated(isp2)) deallocate(isp2)
    if(allocated(iep2)) deallocate(iep2)
    rowsize = meg*kimg_t
    if(rowsize*meg<=maxmatsize.or.sw_scalapack==OFF) then
       npart2 = 1
       allocate(isp2(npart2))
       allocate(iep2(npart2))
       isp2(npart2) = 1
       iep2(npart2) = meg
       block = meg
    else
       npart2 = (rowsize*meg)/maxmatsize
       if(npart2 > meg) npart2 = meg
       block = meg/npart2
       allocate(isp2(npart2))
       allocate(iep2(npart2))
       is = 0
       ie = 0
       do i=1,npart2
          is = ie+1
          ie = is+block-1
          if(i<=mod(meg,npart2)) ie=ie+1
          isp2(i) = is
          iep2(i) = ie
       end do
       iep2(npart2) = meg
    end if
#endif
! -->

    if(printable.and.(sw_scalapack==ON.or.iprisubmat>=2)) then
       write(nfout,*) '== col partition == '
       write(nfout,*) 'target memory(MB)=',msize_submat
       write(nfout,*) 'npart=',npart
       write(nfout,*) 'maxsize=',maxsize
       write(nfout,*) 'memory(MB)=',maxsize*8/1024/1024
       write(nfout,*) 'n, isp, iep, size'
       do i=1,npart
          is = isp(i)
          ie = iep(i)
          size = is*(ie-is+1)+(ie-is+2)*(ie-is+1)/2
          write(nfout,'(4i7)') i, is, ie, size
       end do
       write(nfout,*) 'npart2=',npart2
#ifndef SUBMAT_DGEMM
       write(nfout,*) 'block=',block
       write(nfout,*) 'memory(MB)=',block*rowsize*8/1024/1024
#endif
       write(nfout,*) 'n, isp2, iep2, block'
       do i=1,npart2
          write(nfout,'(4i7)') i, isp2(i), iep2(i), iep2(i)-isp2(i)+1
       end do
    else if(iprisubmat >=2) then
       write(nfout,'(" !target memory(MB) =",i10)') msize_submat
       write(nfout,'(" npart, npart2 = ",2i8)') npart, npart2
       write(nfout,'(" -- isp, iep --")')
       do i = 1, npart
          write(nfout,'(1x,i3," isp(",i3,"), iep(",i3,") = ",2i8)') i,i,i,isp(i),iep(i)
       end do
       write(nfout,'(" -- isp2, iep2 --")')
       do i = 1, npart
          write(nfout,'(1x,i3," isp2(",i3,"), iep2(",i3,") = ",2i8)') i,i,i,isp2(i),iep2(i)
       end do
    end if
    call tstatc0_end(id_sname)
  end subroutine set_col_partition_G


#ifdef _USE_SCALAPACK_
#endif


! ======================================= added by K. Tagami ================= 11.0
  subroutine evolve_WFs_in_subspace_noncl( ik, meg, damp, ekin, afft_kt, bfft_kt, &
	&                                  vnlph_noncl )
!
! Revised by T. Yamasaki, September 2009 : SUBMAT_DGEMM
!
    integer, intent(in) :: ik, meg
    real(kind=DP), intent(in)  :: damp,ekin(kg1)
    real(kind=DP), intent(inout)  :: afft_kt(nfft,ndim_chgpot)
    real(kind=DP), intent(out) :: bfft_kt(nfft,ndim_spinor)
    real(kind=DP), intent(in)  :: vnlph_noncl(kg1,np_e,kimg,ndim_spinor)

! (allocatable variables)
    real(kind=DP), allocatable,dimension(:) ::     eig
    real(kind=DP), allocatable,dimension(:) ::     eko_d
    real(kind=DP), allocatable,dimension(:) ::   w1hw2
    real(kind=DP), allocatable,dimension(:,:) ::   zz

    real(kind=DP), allocatable,dimension(:,:,:,:) :: zat_l_noncl
    real(kind=DP), allocatable,dimension(:,:,:,:) :: zah_l_noncl

    real(kind=DP), allocatable,dimension(:,:,:,:) :: zat_wk_noncl

    real(kind=DP), allocatable :: amat(:,:), zmat(:,:)
    integer, save :: lwork,lrwork,liwork,lda,occ

#ifdef _USE_SCALAPACK_
    integer, save :: lwork1,lrwork1,liwork1
    integer, save :: lwork2,lrwork2,liwork2
#endif

#ifdef SINGLE_CONTEXT
    integer, save :: ictxt, myrow, mycol
#else
    integer, save :: myrow, mycol, tmpctxt
    integer, allocatable, save :: ictxt(:) ! d(0:nrank_k-1)
#endif
    integer, dimension(9), save :: desca,descz
#ifdef SINGLE_CONTEXT
    integer, allocatable, save :: usermap(:,:)
#else
    integer, allocatable, save :: usermap(:,:,:) ! d(nprow,npcol,0:nrank_k-1)
#endif
!!$    integer, save :: npart,npart2
!!$    integer, allocatable, save :: isp(:),iep(:)
!!$    integer, allocatable, save :: isp2(:),iep2(:)
!    logical, save :: first = .true.
    integer :: is,ie,m,jp
    integer :: msize

    integer       :: ib1,ib2,ib1to,i1,ii,ri,ib, i, kimg_t, ig1
    real(kind=DP) :: denom, eko1, eko2, ekod
    real(kind=DP) :: dr1,dr2,di1,di2, dd
    real(kind=DP) :: sum_sq_diagonal, sum_sq_non_diagonal &
         & , sum_abs_diagonal, sum_abs_non_diagonal

    integer :: k1, ispinor
    real(kind=DP) :: dd_max

    integer :: id_sname = -1, id_sname2 = -1, ipri0
    integer :: max_block_size
#ifdef SUBMAT_DGEMM
    real(kind=DP),allocatable,dimension(:,:) :: w1hw2r,w1hw2i
    real(kind=DP) :: alpha, beta
#endif
! <-- T.Kokubo & D.Fukata, Feb. 2010
    integer :: ibsta, ibend, nblk1, nblk2, ibsize, ibsize2

    if(nblocksize_submat_is_given) then
      ibsize = nblocksize_submat
    else
      ibsize = nb_submat_default
    endif
! -->
    call tstatc0_begin('evolve_WFs_in_subspace_noncl ', id_sname,1)

    denom = 1.d0/product(fft_box_size_WF(1:3,1))

! (allocate)  zat_l will be divided by G-index
    if(k_symmetry(ik) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = kimg
    end if

    allocate(eko_d(neg));     eko_d = 0.d0
    allocate(eig(meg)); eig=0.d0
#ifdef _ODD_BOUNDARY_
    if(mod(np_g1k(ik),2) == 0) then
       np_g1k_x = np_g1k(ik) + 1
    else
       np_g1k_x = np_g1k(ik)
    end if
#else
    np_g1k_x = np_g1k(ik)
#endif
    allocate( zat_l_noncl(np_g1k_x,neg,kimg,ndim_spinor))
    zat_l_noncl = 0.0d0


    do ib1 = 1, neg
       if ( map_e(ib1) == myrank_e) eko_d(ib1) = eko_l(map_z(ib1),ik)  ! MPI
    end do

    call mpi_allreduce(MPI_IN_PLACE, eko_d,neg,mpi_double_precision,mpi_sum, mpi_k_world(myrank_k),ierr)

    eko1 = sum(eko_d(1:neg))

!( tenchi ) (zat_l <- zaj_l)
! --> T. Yamasaki, 28th Aug. 2009
!!$    call m_ES_W_transpose(ista_k,iend_k,ik,zaj_l,zat_l)


    Do is=1, ndim_spinor
      call m_ES_W_transpose_r(.true., ista_k, iend_k, ik+is-1, zaj_l, zat_l_noncl(:,:,:,is) )
    End do

! (zaj_l <- H |phi> )
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(bfft,ib1to,ib,ii,i1,dr1,di1,dr2)
!CDIR CONCUR(BY=1)
#endif
    do ib1 = ista_e, iend_e, istep_e     ! MPI
       ib1to=nrvf_ordr(ib1,ik)
       if(ib1to>meg) cycle

       Do is=1, ndim_spinor
         call m_ES_WF_in_Rspace( ik+is-1, ib1, bfft_kt(:,is) )   !(swffft)
       End do
       call m_FFT_Vlocal_W_noncl( afft_kt, bfft_kt, ndim_chgpot, ndim_spinor )
                                                 ! (afft, bfft)-> (bfft)
       Do is=1, ndim_spinor
         call m_FFT_WF( ELECTRON, nfout, bfft_kt(:,is), DIRECT, ON )
       End do

       ib = map_z(ib1)                 ! MPI
       if ( kimg == 1 ) then
          Do is=1, ndim_spinor
            k1 = ik + is -1
            do ii=1,iba(ik)
               i1  = igf(nbase(ii,ik))
               dr1 = zaj_l( ii, ib, k1, 1 )
               dr2 = bfft_kt(i1,is) *denom
               zaj_l( ii, ib, k1, 1 ) = ekin(ii)*dr1 &
	&                              +vnlph_noncl( ii, ib, 1, is ) +dr2
            enddo
          End do
       else
          Do is=1, ndim_spinor
            k1 = ik + is -1
            do ii=1,iba(ik)
               i1  = igf(nbase(ii,ik))
               dr1  = zaj_l( ii, ib, k1, 1 )
               di1  = zaj_l( ii, ib, k1, 2 )
               zaj_l( ii, ib, k1, 1 ) = ekin(ii)*dr1 +bfft_kt( 2*i1-1,is )*denom &
	&                              +vnlph_noncl( ii, ib, 1, is )
               zaj_l( ii, ib, k1, 2 ) = ekin(ii)*di1 +bfft_kt( 2*i1,  is )*denom &
	&                              +vnlph_noncl( ii, ib, 2, is )
            enddo
          End do
       end if
#ifdef SAVE_FFT_TIMES
       if(sw_save_fft == ON) then
          do is =1,ndim_spinor
             k1=ik+is-1
             status_saved_phifftr(ib,k1) = OLD
          end do
       end if
#endif
    end do
! ------------------ok  -


!( tenchi ) (zah_l <- zaj_l)
    allocate( zah_l_noncl( np_g1k_x, neg, kimg, ndim_spinor ) )
    zah_l_noncl = 0.0d0

! --> T. Yamasaki 28th Aug. 009
!!$    call m_ES_W_transpose(ista_k,iend_k,ik,zaj_l,zah_l)

    Do is=1, ndim_spinor
      call m_ES_W_transpose_r(.true., ista_k, iend_k, ik+is-1, zaj_l, zah_l_noncl(:,:,:,is) )
    End do

! <--

! (make matrix elements )
#ifdef _USE_SCALAPACK_
    if(sw_scalapack == ON) then
       if(submat_uncalled) then
! ==== ASMS 2021/03/19
          if(block_size == 0) then
             max_block_size = int(real(meg)/real(nrank_e))
             block_size =  nb_mgs_default
             if(block_size.ge.max_block_size) block_size = max_block_size
          end if
! ==== ASMS 2021/03/19
          call set_nprow_npcol(nprow,npcol)
#ifdef SINGLE_CONTEXT
          if(allocated(usermap)) deallocate(usermap)
          allocate(usermap(nprow,npcol))
#else
          if(allocated(ictxt)) deallocate(ictxt)
          if(allocated(usermap)) deallocate(usermap)
          allocate(ictxt(0:nrank_k-1))
          allocate(usermap(nprow,npcol,0:nrank_k-1))
#endif
!         call scalapack_setup(meg,lwork,lrwork,liwork,lda,occ,ictxt,myrow,mycol,desca,descz,usermap)
#ifdef SINGLE_CONTEXT
          call scalapack_setup(meg,lwork1,lrwork1,liwork1,lwork2,lrwork2,liwork2,lda,occ,ictxt,myrow,mycol,desca,descz,usermap)
#else
          call scalapack_setup(meg,lwork1,lrwork1,liwork1,lwork2,lrwork2,liwork2,lda,occ,ictxt(myrank_k),&
          & myrow,mycol,desca,descz,usermap)
#endif
! ===== We need more work area!!! ==============================================
          lwork1  = 10*lwork1
          lrwork1 = 10*lrwork1
          liwork1 = 10*liwork1
          lwork2  = 10*lwork2
          lrwork2 = 10*lrwork2
          liwork2 = 10*liwork2
! ==============================================================================
       end if
       allocate(amat(lda*kimg_t,occ))
       allocate(zmat(lda*kimg_t,occ))
    else
       allocate(zmat(meg*kimg_t,meg))
    end if
#else
    allocate(zmat(meg*kimg_t,meg))
    zmat=0.0d0
#endif
    if(submat_uncalled) then
!!$       call set_col_partition ! -> npart,isp,iep, npart2, isp2, iep2
       call set_col_partition_G(kimg_t,meg,ibsize) ! -> npart,isp,iep, npart2, isp2, iep2
       submat_uncalled = .false.
    end if

    call tstatc0_begin('evolve_WFs_in_subspace(PART1) ', id_sname2)
!!$    if(iprisubmat>=2) then
!!$       write(nfout,'(" --- neordr ---")')
!!$       write(nfout,'(" ik = ",i8)') ik
!!$       write(nfout,'(10i8)') (neordr(ib1,ik),ib1=1, neg)
!!$       write(nfout,'(" --- nrvf_ordr ---")')
!!$       write(nfout,'(" ik = ",i8)') ik
!!$       write(nfout,'(10i8)') (nrvf_ordr(ib1,ik),ib1=1, neg)
!!$    end if

    ekod = 0.d0
    sum_abs_diagonal = 0.d0
    sum_sq_diagonal = 0.d0
    sum_sq_non_diagonal = 0.d0
    sum_abs_non_diagonal = 0.d0

! ------------------------------------ part 1 -----------
    PART1: do jp = 1,npart
       is = isp(jp)
       ie = iep(jp)
       msize = is*(ie-is+1) + (ie-is+1)*(ie-is)/2
       allocate(w1hw2(msize*kimg_t))
       w1hw2 = 0.d0
       if(iprisubmat>=2) then
          write(nfout,'(" jp, is, ie = ",3i8)') jp, is, ie
       end if

#ifdef SUBMAT_DGEMM
! <-- T.Kokubo & D.Fukata Feb. 2010
       nblk1 = ie-is+1
#ifdef SX
       ibsize2 = min(ie,meg)
#else
       ibsize2 = ibsize
#endif

       allocate(w1hw2r(ibsize2,nblk1)); w1hw2r = 0.d0
       if(kimg==2 .and. k_symmetry(ik) /= GAMMA) then
          allocate(w1hw2i(ibsize2,nblk1)); w1hw2i = 0.d0
       end if

       do ibsta=1,min(ie,meg), ibsize2
          ibend=min( ibsta+ibsize2-1, min(ie,meg) )
          nblk2=ibend-ibsta+1

          w1hw2r = 0.d0
          if ( kimg==2 .and. k_symmetry(ik) /= GAMMA ) w1hw2i = 0.d0

          if ( kimg==1 ) then
             Do ispinor=1, ndim_spinor
               call dgemm( 'T', 'N', nblk2, nblk1, np_g1k(ik), 1.d0,&
	&                  zat_l_noncl(1,ibsta,1,ispinor), np_g1k_x,  &
                  &        zah_l_noncl(1,is,1,ispinor),    np_g1k_x, &
                  &        1.d0, w1hw2r, ibsize2 )
             End do

             do ib1 = max(is,ibsta), min(ie,meg)
                do ib2 = ibsta,min(ibend,ib1)
                   m = ib1*(ib1-1)/2-is*(is-1)/2+ib2
                   w1hw2(m) = w1hw2r(ib2-ibsta+1, ib1-max(is,ibsta)+1)
                end do
             end do
          else
             if ( k_symmetry(ik) == GAMMA ) then
                ig1=1
                if ( myrank_e==0 ) ig1=2

                Do ispinor=1, ndim_spinor
                  call dgemm( 'T', 'N', nblk2, nblk1, np_g1k(ik), 1.d0, &
	&                      zat_l_noncl(1,ibsta,1,ispinor), np_g1k_x, &
        &                      zah_l_noncl(1,is,1,ispinor),    np_g1k_x, &
	&                      1.d0, w1hw2r, ibsize2 )
                  call dgemm( 'T', 'N', nblk2, nblk1, np_g1k(ik), 1.d0, &
	&                      zat_l_noncl(1,ibsta,2,ispinor), np_g1k_x, &
        &                      zah_l_noncl(1,is,2,ispinor),    np_g1k_x, &
	&                      1.d0, w1hw2r, ibsize2 )
                End do

                w1hw2r = w1hw2r*2.d0
                if ( myrank_e == 0 ) then
                   do ib1 = max(is,ibsta), min(ie,meg)
                      do ib2 = ibsta,min(ibend,ib1)
                         Do ispinor=1, ndim_spinor
                            w1hw2r(ib2-ibsta+1, ib1-max(is,ibsta)+1) &
	&                   =  w1hw2r(ib2-ibsta+1, ib1-max(is,ibsta)+1) &
        &                     - 2.d0 *zah_l_noncl(1,ib1,2,ispinor) &
	&                            *zat_l_noncl(1,ib2,2,ispinor) &
        &                     -       zah_l_noncl(1,ib1,1,ispinor) &
	&                            *zat_l_noncl(1,ib2,1,ispinor)
                         End do
                      end do
                   end do
                end if

                do ib1 = max(is,ibsta), min(ie,meg)
                   do ib2 = ibsta,min(ibend,ib1)
                      m = ib1*(ib1-1)/2-is*(is-1)/2+ib2
                      w1hw2(m)=w1hw2r(ib2-ibsta+1, ib1-max(is,ibsta)+1)
                   end do
                end do
             else
                Do ispinor=1, ndim_spinor
                  call dgemm( 'T', 'N', nblk2, nblk1, np_g1k(ik), 1.d0, &
	&                     zat_l_noncl(1,ibsta,1,ispinor), np_g1k_x, &
        &                     zah_l_noncl(1,is,1,ispinor),    np_g1k_x, &
	&                     1.d0, w1hw2r, ibsize2 )
                  call dgemm( 'T', 'N', nblk2, nblk1, np_g1k(ik), 1.d0, &
	&                     zat_l_noncl(1,ibsta,2,ispinor), np_g1k_x, &
        &                     zah_l_noncl(1,is,2,ispinor),    np_g1k_x, &
	&                     1.d0, w1hw2r, ibsize2 )
                  call dgemm( 'T', 'N', nblk2, nblk1, np_g1k(ik), 1.d0, &
	&                     zat_l_noncl(1,ibsta,2,ispinor), np_g1k_x, &
        &                     zah_l_noncl(1,is,1,ispinor),    np_g1k_x, &
	&                     1.d0, w1hw2i, ibsize2 )
                  call dgemm( 'T', 'N', nblk2, nblk1, np_g1k(ik), -1.d0, &
	&                     zat_l_noncl(1,ibsta,1,ispinor), np_g1k_x, &
        &                     zah_l_noncl(1,is,2,ispinor),    np_g1k_x, &
	&                     1.d0, w1hw2i, ibsize2 )
                End do

                do ib1 = max(is,ibsta), min(ie,meg)
                   do ib2 = ibsta,min(ibend,ib1)
                      m = ib1*(ib1-1)/2-is*(is-1)/2+ib2
                      w1hw2(2*m-1) = w1hw2r(ib2-ibsta+1, ib1-max(is,ibsta)+1)
                      w1hw2(2*m  ) = w1hw2i(ib2-ibsta+1, ib1-max(is,ibsta)+1)
                   end do
                end do
             end if
          end if
       enddo

       if (kimg==2 .and. k_symmetry(ik) /= GAMMA) deallocate(w1hw2i)
       deallocate(w1hw2r)

#else
! parallel loop
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(ib2,ii,dr1,dr2,di1,di2,m)
!CDIR CONCUR(BY=1)
#endif
! --> T. Yamasaki 28th Aug. 2009
!   ib1to --> ib1,  instead of ib1to = neordr(ib1,ik)
!   ib2to --> ib2
       do ib1 = is, ie
          if(ib1>meg) cycle
          do ib2 = 1, ib1
             m = ib1*(ib1-1)/2-is*(is-1)/2+ib2

             if (kimg == 1) then
	        Do ispinor=1, ndim_spinor
                  do ii = 1, np_g1k(ik)              ! MPI
                     dr1 = zah_l_noncl( ii, ib1, 1, ispinor )
                     dr2 = zat_l_noncl( ii, ib2, 1, ispinor )
                     w1hw2(m) = w1hw2(m) + dr1*dr2
                  end do
                End do
             else
                if (k_symmetry(ik) == GAMMA) then
                   dd = 0.d0
                   ig1 = 1
	           if(myrank_e == 0) ig1 = 2

                   Do ispinor=1, ndim_spinor
                     do ii = ig1, np_g1k(ik)            ! MPI
                        dr1 = zah_l_noncl( ii, ib1, 1, ispinor )
                        di1 = zah_l_noncl( ii, ib1, 2, ispinor )
                        dr2 = zat_l_noncl( ii, ib2, 1, ispinor )
                        di2 = zat_l_noncl( ii, ib2, 2, ispinor )
                        dd = dd + ( dr1*dr2 +di1*di2 )*2.d0
                     end do
                   End do

                   if ( myrank_e == 0 ) then
	              Do ispinor=1, ndim_spinor
                         dd = dd + zah_l_noncl(1,ib1,1,ispinor) &
	&                         *zat_l_noncl(1,ib2,1,ispinor)
                      End do
                   endif
                   w1hw2(m)=dd
                else
	           Do ispinor=1, ndim_spinor
                     do ii = 1, np_g1k(ik)             ! MPI
                        dr1 = zah_l_noncl( ii, ib1, 1, ispinor )
                        di1 = zah_l_noncl( ii, ib1, 2, ispinor )
                        dr2 = zat_l_noncl( ii, ib2, 1, ispinor )
                        di2 = zat_l_noncl( ii, ib2, 2, ispinor )
                        w1hw2(2*m-1)=w1hw2(2*m-1) +dr1*dr2 +di1*di2
                        w1hw2(2*m  )=w1hw2(2*m  ) +dr1*di2 -di1*dr2
                     end do
                   End do
                end if
             end if
          end do
       end do
#endif

       if (iprisubmat >= 2) call wd_w1hw2(" -- just after making w1hw2 --")

!! (off-diagonal matrix elements will be damped.)
       if (kimg_t == 1) then
          do ib1=max(2,is),ie
             do ib2=1,ib1-1
                m = ib1*(ib1-1)/2-is*(is-1)/2+ib2
                sum_sq_non_diagonal = sum_sq_non_diagonal + w1hw2(m)*w1hw2(m)
                sum_abs_non_diagonal = sum_abs_non_diagonal + abs(w1hw2(m))
                w1hw2(m) = damp*w1hw2(m)
             enddo
          enddo
          do ib = is, ie
             m = ib*(ib-1)/2-is*(is-1)/2+ib
             sum_sq_diagonal = sum_sq_diagonal + w1hw2(m)**2
             sum_abs_diagonal = sum_abs_diagonal + abs(w1hw2(m))
             ekod = ekod + w1hw2(m)
          end do
       else
          do ib1=max(2,is),ie
             do ib2=1,ib1-1
                m = ib1*(ib1-1)/2-is*(is-1)/2+ib2
                dd = w1hw2(2*m-1)*w1hw2(2*m-1) + w1hw2(2*m)*w1hw2(2*m)
                sum_sq_non_diagonal = sum_sq_non_diagonal + dd
                if(dd >= SmallestPositiveNumber) sum_abs_non_diagonal = sum_abs_non_diagonal + sqrt(dd)
                w1hw2(2*m-1) = damp*w1hw2(2*m-1)
                w1hw2(2*m  ) = damp*w1hw2(2*m  )

             enddo
          enddo
          do ib = is, ie
             m = ib*(ib-1)/2-is*(is-1)/2+ib
             sum_sq_diagonal = sum_sq_diagonal + w1hw2(2*m-1)**2
             sum_abs_diagonal = sum_abs_diagonal + abs(w1hw2(2*m-1))
             ekod = ekod + w1hw2(2*m-1)
          end do
       endif
       if(iprisubmat >= 2) call wd_w1hw2(" -- after damping --")

#ifdef _USE_SCALAPACK_
       if (sw_scalapack==ON) then
          call set_amat(lda,kimg_t,occ,amat,meg,nprow,npcol,block_size,usermap,is,ie,w1hw2)
       else
#endif
          call set_hmat(kimg_t,meg,is,ie,w1hw2,zmat)
#ifdef _USE_SCALAPACK_
       end if
#endif
       deallocate(w1hw2)
    end do PART1

    call tstatc0_end(id_sname2)

    deallocate(zah_l_noncl)

    sum_abs_diagonal = sum_abs_diagonal/meg
    sum_abs_non_diagonal = sum_abs_non_diagonal/(meg*(meg-1)/2)

    if(iprisubmat >= 2) then
       write(nfout,*) 'neordr for ik = ',ik
       write(nfout,9002) (neordr(ib1,ik),ib1=1,neg)
       write(nfout,*) 'nrvf_ordr for ik = ',ik
       write(nfout,9002) (nrvf_ordr(ib1,ik),ib1=1,neg)
       write(nfout,*) 'eig'
! ===================== ktDEBUG
!       write(nfout,'(5x,10f8.4)') (eko_l(ib1,ik),ib1=1,neg)
       write(nfout,'(5x,10f8.4)') (eko_l(ib1,ik),ib1=1,np_e)
! =====================
    endif
9002 format(5x,10i8)


! ==============================================================

!! (Diagonalization )  !!
#ifdef _USE_SCALAPACK_
     if(sw_scalapack == ON) then
        if(kimg_t == 1) then
!          call pdsyev_driver(meg,eig,desca,descz,lda,occ,amat,zmat,lwork,liwork)
           call pdsyev_driver(meg,eig,desca,descz,lda,occ,amat,zmat,lwork1,liwork1)
        else
!          call pzheev_driver(meg,eig,desca,descz,lda,occ,amat,zmat,lwork,lrwork,liwork)
           call pzheev_driver(meg,eig,desca,descz,lda,occ,amat,zmat,lwork2,lrwork2,liwork2)
        endif
     else
#endif
        if(kimg_t == 1) then
           call dsyev_driver(meg,eig,zmat)
        else
           call zheev_driver(meg,eig,zmat)
        endif
#ifdef _USE_SCALAPACK_
     endif
#endif
!!$     if(kimg_t == 1) then
!!$       call dsyev_driver(meg,eig,zmat)
!!$     else
!!$       call zheev_driver(meg,eig,zmat)
!!$     endif
!!$#endif

     if(iprisubmat >= 2) then
        write(nfout,*) 'eko_d for ik = ',ik
        write(nfout,9001) (eko_d(ib),ib=1,meg)
        write(nfout,*) 'eig for ik = ',ik
        write(nfout,9001) (eig(ib),ib=1,meg)
        !!$call wd_w1hw2(" -- after diagonalization --")
        !!$call wd_zmat(" -- after diagonalization --")
!sum eko
     end if

     call get_ipri0(iprisubmat,ipri0)
     if(ipri0>=2) then
        dr1=0.d0;dr2=0.d0
        do ib1=1,meg
           ib1to = neordr(ib1,ik)
           if(map_e(ib1to) == myrank_e) dr1=dr1+eko_l(map_z(ib1to),ik)
           dr2=dr2+eig(ib1)
        enddo
        call mpi_allreduce(MPI_IN_PLACE, dr1,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
     end if

     if(iprisubmat>=2) then
        write(nfout,*) ' sum of eko_l & eig =',dr1,dr2
        write(nfout,'(" sum_sq_diagonal , sum_sq_non_diagonal  = ",2d20.12)') sum_sq_diagonal, sum_sq_non_diagonal
        write(nfout,'(" sum_abs_diagonal ,sum_abs_non_diagonal = ",2d20.12)') sum_abs_diagonal, sum_abs_non_diagonal
     end if

! ------

!! (subspace rotation) !!
     allocate( zat_wk_noncl( np_g1k_x, meg, kimg, ndim_spinor ))
     zat_wk_noncl = 0.d0
! <-- T.Kokubo & D.Fukata, Feb. 2010
        PART2: do jp = 1,npart2
           is = isp2(jp)
           ie = iep2(jp)
           nblk1 = ie-is+1
#ifdef SX
           ibsize2 = meg
#else
           ibsize2 = ibsize
#endif
#ifdef _USE_SCALAPACK_
     if (sw_scalapack == ON) then
           allocate(zz(meg*kimg_t,nblk1))
           call put_zmat(kimg_t,lda,occ,zmat,meg,nprow,npcol,block_size,usermap,is,ie,zz)

           if ( kimg == 1) then
              call subspace_rotation_real_noncl( is, ie, zat_wk_noncl, zz, &
	&                                        nblk1, ibsize2 )
           else
              call subspace_rotation_imag_noncl( is, ie, zat_wk_noncl, zz, &
	&                                        nblk1, ibsize2 )
           endif
           deallocate(zz)
     else
#endif
           if (kimg == 1) then
              call subspace_rotation_real_noncl( is, ie, zat_wk_noncl, zmat(1,is),&
	&                                        nblk1, ibsize2 )
           else
              call subspace_rotation_imag_noncl( is, ie, zat_wk_noncl, zmat(1,is),&
	&                                        nblk1, ibsize2 )
           endif
#ifdef _USE_SCALAPACK_
     end if
#endif
        end do PART2
!
!#ifdef _USE_SCALAPACK_
!     if(sw_scalapack == ON) then
!        PART2: do jp = 1,npart2
!           is = isp2(jp)
!           ie = iep2(jp)
!           allocate(zz(meg*kimg_t,is:ie))
!           call put_zmat(lda,occ,zmat,meg,nprow,npcol,block_size,usermap,is,ie,zz)
!           if(kimg == 1) then
!              call subspace_rotation_real(is,ie,zat_wk,zz)
!           else
!              call subspace_rotation_imag(is,ie,zat_wk,zz)
!           endif
!           deallocate(zz)
!        end do PART2
!     else
!#endif
!        is=1; ie=meg
!        if(kimg == 1) then
!           call subspace_rotation_real(is,ie,zat_wk,zmat)
!        else
!           call subspace_rotation_imag(is,ie,zat_wk,zmat)
!        endif
!#ifdef _USE_SCALAPACK_
!     end if
!#endif
! -->

     Do ispinor=1, ndim_spinor
       zat_l_noncl(1:np_g1k(ik),1:meg,1:kimg,ispinor) &
	&  = zat_wk_noncl(1:np_g1k(ik),1:meg,1:kimg,ispinor)
     End do
     deallocate(zat_wk_noncl)

!! (zaj_l) by tenchi

     Do ispinor=1, ndim_spinor
       call m_ES_W_transpose_back_r(.false., ista_k, iend_k, ik+ispinor-1, &
	&                          zaj_l, zat_l_noncl(:,:,:,ispinor) )
                                                     ! zat_l -> zaj_l
     End do


     if (meg < neg) then
        Do ispinor=1, ndim_spinor
           call m_ES_W_transpose_back2( ista_k, iend_k, ik +ispinor -1, &
	&                               meg, neordr, zaj_l, &
	&                               zat_l_noncl(:,:,:,ispinor)  )  ! zat_l -> zaj_l
        End do
     end if

#ifdef SAVE_FFT_TIMES
     if(sw_save_fft == ON) then
        Do ispinor=1, ndim_spinor
           do ib = 1, np_e
              status_saved_phifftr(ib,ik+ispinor-1) = OLD
           end do
        end Do
     end if
#endif

!! (eko_l)
     if(meg < neg) then
        do ib1=meg+1,neg
           if(map_e(ib1) == myrank_e) then        ! MPI
              eko_l(map_z(ib1),ik)=eko_d(neordr(ib1,ik))
           end if
        enddo
     endif
     do ib1 = 1, meg
        if(map_e(ib1) == myrank_e) then         ! MPI
           eko_l(map_z(ib1),ik)=eig(ib1)
        end if
     end do
     if(iprisubmat >= 2) then
        eko2 = sum(eig(1:meg))
        write(nfout,1201) ik,eko1,ekod,eko2

        write(nfout,*) 'eko_l'
        write(nfout,9001) (eko_l(ib1,ik),ib1=1,meg)
     endif
1201 format(' %% for ik = ',i4,4x,' eko1&ekod&eko2 = ',3f14.7)
9001 format(5x,6f12.5)
!! (neordr & nrvf_ordr)

     neordr(1:neg,ik) = (/(ib1,ib1=1,neg)/)
     nrvf_ordr(1:neg,ik) = (/(ib1,ib1=1,neg)/)

! (deallocate)
    deallocate(eko_d)
    deallocate(eig)
    deallocate(zat_l_noncl)

    !!$deallocate(zah_l)
    if (sw_scalapack==ON) deallocate(amat)
    deallocate(zmat)

! === DEBUG by tkato 2012/12/11 ================================================
    call mpi_allreduce(MPI_IN_PLACE,sum_abs_non_diagonal,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
    call mpi_allreduce(MPI_IN_PLACE,sum_abs_diagonal,    1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
! ==============================================================================

    if(.not.non_diagonal_part_is_small(ik)  &
         & .and. sum_abs_non_diagonal/sum_abs_diagonal <= submat_critical_ratio) &
         & non_diagonal_part_is_small(ik) = .true.

    if(iprisubmat >= 2) write(nfout,'(" !submat  non_diagonal/diagonal = ",d20.8)') &
         & sum_abs_non_diagonal/sum_abs_diagonal

    call tstatc0_end(id_sname)

  contains
    subroutine get_ipri0(ipri_in, ipri_out)
      integer, intent(in)  :: ipri_in
      integer, intent(out) :: ipri_out
      if(nrank_e>1) then
         if(myrank_e == 0) ipri_out = ipri_in
         call mpi_bcast(ipri_out,1,mpi_integer,0,mpi_k_world(myrank_k),ierr)
      else
         ipri_out = ipri_in
      end if
    end subroutine get_ipri0

    subroutine wd_w1hw2(somecomment)
      character(len=*), intent(in) :: somecomment
      integer :: ib1, ib2, neg_wd1, neg_wd2
      write(nfout,'(a35)') somecomment
      write(nfout,*) 'w1hw2 for ik = ',ik
      neg_wd1 = 8
      neg_wd2 = 8
      if(neg_wd1 > ie) neg_wd1 = ie
      if(neg_wd2 > meg) neg_wd2 = meg
      write(nfout,'(" real part (",i8,":",i8,",1:",i8,")")') is,neg_wd1,neg_wd2
      if(kimg_t == 1) then
         do ib1=is,neg_wd1
            write(nfout,9001) (w1hw2(ib1*(ib1-1)/2-is*(is-1)/2+ib2),ib2=1,neg_wd2)
         enddo
         write(nfout,'(" real part (",i8,":",i8,",1:",i8,")")') ie-neg_wd1+1,ie, neg_wd2
         do ib1=ie-neg_wd1+1, ie
            write(nfout,9001) (w1hw2(ib1*(ib1-1)/2-is*(is-1)/2+ib2),ib2=1,neg_wd2)
         enddo
      else
         do ib1=is,neg_wd1
            write(nfout,9001) (w1hw2(2*(ib1*(ib1-1)/2-is*(is-1)/2+ib2)-1),ib2=1,neg_wd2)
         enddo
         write(nfout,'(" real part (",i8,":",i8,",1:",i8,")")') ie-neg_wd1+1,ie, neg_wd2
         do ib1=ie-neg_wd1+1, ie
            write(nfout,9001) (w1hw2(2*(ib1*(ib1-1)/2-is*(is-1)/2+ib2)-1),ib2=1,neg_wd2)
         enddo
         write(nfout,'(" imag part")')
         do ib1=is,neg_wd1
            write(nfout,9001) (w1hw2(2*(ib1*(ib1-1)/2-is*(is-1)/2+ib2)),ib2=1,neg_wd2)
         enddo
      end if
9001  format(1x,8e14.6)
      write(nfout,*) 'eko_l for ik = ',ik
      write(nfout,9001) (eko_d(neordr(ib1,ik)),ib1=1,neg_wd2)
    end subroutine wd_w1hw2


    subroutine subspace_rotation_real_noncl( is, ie, zat_wk_noncl, zz, ibsize1, ibsize2 )

      integer, intent(in) :: is,ie, ibsize1,ibsize2
      real(kind=DP), intent(inout) :: zat_wk_noncl( np_g1k_x, meg, kimg, ndim_spinor )
      real(kind=DP), intent(in) :: zz(meg,ibsize1)
      integer :: ibsta, ibend
#ifdef SUBMAT_DGEMM
      integer :: nblk2
      real(kind=DP) :: alpha, beta
#endif
      integer :: id_sname = -1

      call tstatc0_begin('subspace_roation_real_noncl ', id_sname)

#ifdef SUBMAT_DGEMM
      alpha=1.d0; beta=1.d0
      do ibsta=1,meg, ibsize2
         ibend=min(ibsta+ibsize2-1,meg)
         nblk2=ibend-ibsta+1

         Do ispinor=1, ndim_spinor
            call dgemm( 'N', 'N', np_g1k(ik), ibsize1, nblk2, alpha, &
	&               zat_l_noncl(1,ibsta,1,ispinor), np_g1k_x, &
	&               zz(ibsta,1), meg, beta, &
	&               zat_wk_noncl(1,is,1,ispinor), np_g1k_x )
         End do
      enddo
#else
      do ibsta=1,meg, ibsize2
         ibend=min(ibsta+ibsize2-1,meg)
         do ib2=is,ie
            do ib1=ibsta,ibend
               Do ispinor=1, ndim_spinor
                 do ii=1,np_g1k(ik)
                    zat_wk_noncl( ii, ib2, 1, ispinor ) &
	&              = zat_wk_noncl( ii, ib2, 1, ispinor ) &
	&               + zz( ib1, ib2-is+1 ) *zat_l_noncl( ii, ib1, 1, ispinor )
                 enddo
               End do
            enddo
         enddo
      enddo
#endif
      call tstatc0_end(id_sname)
    end subroutine subspace_rotation_real_noncl

    subroutine subspace_rotation_imag_noncl( is, ie, zat_wk_noncl, zz, ibsize1, ibsize2 )
      integer, intent(in) :: is,ie, ibsize1,ibsize2
      real(kind=DP), intent(inout) :: zat_wk_noncl( np_g1k_x, meg, kimg, ndim_spinor )
      real(kind=DP), intent(in) :: zz(meg*kimg_t,ibsize1)
      integer :: id_sname = -1
      integer :: ibsta, ibend

      integer :: ispinor

#ifdef SUBMAT_DGEMM
      real(kind=DP) :: alpha, beta
      integer :: nblk2
      real(kind=DP), allocatable, dimension(:,:) :: zzr,zzi
#endif

      call tstatc0_begin('subspace_roation_imag_noncl ', id_sname)

!      stop

      if (k_symmetry(ik) == GAMMA) then
#ifdef SUBMAT_DGEMM
         alpha = 1.d0; beta = 1.d0
         do ibsta=1,meg,ibsize2
            ibend=min(ibsta+ibsize2-1,meg)
            nblk2=ibend-ibsta+1

            Do ispinor=1, ndim_spinor
              call dgemm( 'N', 'N', np_g1k(ik), ibsize1, nblk2, alpha, &
	&                  zat_l_noncl(1,ibsta,1,ispinor), np_g1k_x, &
	&                  zz(ibsta,1), meg, beta, &
        &                  zat_wk_noncl(1,is,1,ispinor), np_g1k_x )
              call dgemm( 'N', 'N', np_g1k(ik), ibsize1, nblk2, alpha, &
	&                  zat_l_noncl(1,ibsta,2,ispinor), np_g1k_x, &
	&                  zz(ibsta,1), meg, beta, &
        &                  zat_wk_noncl(1,is,2,ispinor), np_g1k_x )
            End do
         enddo
         !!if(iprisubmat>=2) write(nfout,'(" !!<<subspace_rotation_imag>> ibsize = ",i8)') ibsize
#else
         do ibsta=1,meg,ibsize2
            ibend=min(ibsta+ibsize2-1,meg)
            do ib2=is,ie
!CDIR OUTERUNROLL=4
               do ib1=ibsta,ibend
                  dr1=zz(ib1,ib2-is+1)
	          Do ispinor=1, ndim_spinor
                     do ii=1,np_g1k(ik)
                       zat_wk_noncl( ii, ib2, 1, ispinor ) &
	&               = zat_wk_noncl( ii, ib2, 1, ispinor ) &
	&                + dr1 *zat_l_noncl(ii, ib1, 1, ispinor )
                       zat_wk_noncl( ii, ib2, 2, ispinor ) &
	&               = zat_wk_noncl( ii, ib2, 2, ispinor )&
	&                + dr1 *zat_l_noncl( ii, ib1, 2, ispinor )
                     enddo
                   ENd do
               enddo! ib1 loop
            enddo! ib2 loop
         enddo
#endif
         if ( iprisubmat >= 2) then
            do ib2 = is, ie
	       Do ispinor=1, ndim_spinor
	          write(nfout,*) '*** ispinor = ', ispinor
                  write(nfout,'(" --- zat_l_noncl (input, output) ---, ib =",i5)') ib2
                  write(nfout,'(" zat_l_noncl(input ) Re : ",5e14.6)') &
	&                  ( zat_l_noncl(ii,ib2,1,ispinor),ii=1,5 )
                  write(nfout,'(" zat_l_noncl(input ) Im : ",5e14.6)') &
	&                  ( zat_l_noncl(ii,ib2,2,ispinor),ii=1,5 )
                 write(nfout,'(" zat_l_noncl(output) Re : ",5e14.6)') &
	&                  ( zat_wk_noncl(ii,ib2,1,ispinor),ii=1,5)
                 write(nfout,'(" zat_l_noncl(output) Im : ",5e14.6)') &
	&                  ( zat_wk_noncl(ii,ib2,2,ispinor),ii=1,5)
	       End do
            end do
         end if
      else
#ifdef SUBMAT_DGEMM

         allocate(zzr(ibsize2,ibsize1)); allocate(zzi(ibsize2,ibsize1))
         do ibsta=1,meg,ibsize2
            ibend=min(ibsta+ibsize2-1,meg)
            nblk2=ibend-ibsta+1
            do ib2=is,ie
               do ib1=ibsta, ibend
                  zzr(ib1-ibsta+1,ib2-is+1)=zz(2*ib1-1,ib2-is+1)
                  zzi(ib1-ibsta+1,ib2-is+1)=zz(2*ib1  ,ib2-is+1)
               end do
            end do

	    Do ispinor=1, ndim_spinor
              alpha=1.d0; beta=1.d0
              call dgemm( 'N', 'N', np_g1k(ik), ibsize1, nblk2, alpha, &
	&                 zat_l_noncl(1,ibsta,1,ispinor), np_g1k_x, &
	&                 zzr, ibsize2, beta, &
	&                 zat_wk_noncl(1,is,1,ispinor), np_g1k_x )

              alpha=-1.d0; beta=1.d0
              call dgemm( 'N', 'N', np_g1k(ik), ibsize1, nblk2, alpha, &
	&                  zat_l_noncl(1,ibsta,2,ispinor), np_g1k_x, &
	&                  zzi, ibsize2, beta, &
	&                  zat_wk_noncl(1,is,1,ispinor), np_g1k_x )

              alpha=1.d0; beta=1.d0
              call dgemm( 'N', 'N', np_g1k(ik), ibsize1, nblk2, alpha, &
	&                 zat_l_noncl(1,ibsta,2,ispinor), np_g1k_x, &
	&                 zzr, ibsize2, beta, &
	&                 zat_wk_noncl(1,is,2,ispinor), np_g1k_x )

              alpha=1.d0; beta=1.d0
              call dgemm( 'N', 'N', np_g1k(ik), ibsize1, nblk2, alpha, &
	&                 zat_l_noncl(1,ibsta,1,ispinor), np_g1k_x, &
	&                 zzi, ibsize2, beta, &
	&                 zat_wk_noncl(1,is,2,ispinor), np_g1k_x )
            End do
         enddo
         deallocate(zzi,zzr)
#else
         do ibsta=1,meg,ibsize2
            ibend=min(ibsta+ibsize2-1,meg)
            do ib2=is,ie
!CDIR OUTERUNROLL=4
               do ib1=ibsta,ibend
                  dr1=zz(2*ib1-1,ib2-is+1)
                  di1=zz(2*ib1  ,ib2-is+1)

	          Do ispinor=1, ndim_spinor
                    do ii=1,np_g1k(ik)
                       dr2 = zat_l_noncl( ii, ib1, 1, ispinor )
                       di2 = zat_l_noncl( ii, ib1, 2, ispinor )
                       zat_wk_noncl( ii, ib2, 1, ispinor ) &
	&                = zat_wk_noncl( ii, ib2, 1, ispinor ) +dr1*dr2 -di1*di2
                       zat_wk_noncl( ii, ib2, 2, ispinor ) &
	&                = zat_wk_noncl( ii, ib2, 2, ispinor ) +dr1*di2 +di1*dr2
                    enddo
                  End do
               enddo! ib1 loop
            enddo! ib2 loop
         enddo
#endif
       end if
       call tstatc0_end(id_sname)
    end subroutine subspace_rotation_imag_noncl

#ifdef _USE_SCALAPACK_
#endif
! _USE_SCALAPACK_

  end subroutine evolve_WFs_in_subspace_noncl
! ==================================================================== 11.0



!!$#endif

!! ++++++++++++++++++++ global(=not contained) subroutines in this module +++++++++++++++
#ifndef _USE_DSYEVD_
  subroutine dsyev_driver(ndim,eig,w1hw2)
    integer, intent(in):: ndim
    real(kind=DP), intent(out) ,dimension(ndim) :: eig
    real(kind=DP), intent(inout) ,dimension(ndim,ndim) :: w1hw2

    integer, parameter :: PRINTLEVEL = 2
    character(len=1) :: JOBZ,UPLO
!   integer, parameter :: lda=ndim
    integer :: lda
!!$    integer :: lwork,lwork_min
    integer, save :: lwork=0
    integer       :: lwork_min
    integer,parameter :: nb = 64
    integer, parameter :: lwork_huge=10**8
    real(kind=DP),allocatable,dimension(:) :: work_lapack
    integer :: info

    integer :: id_sname = -1
    call tstatc0_begin('dsyev_driver ', id_sname)

    lda=ndim
!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
    JOBZ = 'V'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
    UPLO = 'U'
! (define lwork)
!!$    lwork = 0
!!$    lwork_min=max(1,3*ndim-1)
    lwork_min=max(1,(nb+2)*ndim)
    if(lwork == 0) then
       lwork = -1
       allocate(work_lapack(1))
       call dsyev(JOBZ,UPLO,ndim,w1hw2,lda,eig,work_lapack,lwork,info)
       lwork = int(work_lapack(1))+1
       if(iprisubmat >= PRINTLEVEL) then
          write(nfout,'(" lwork = ",i8)') lwork
       end if
       deallocate(work_lapack)
    end if
    if(lwork < lwork_min) lwork=lwork_min
    if(lwork > lwork_huge) lwork=lwork_min
    allocate(work_lapack(lwork))
    if(iprisubmat >= PRINTLEVEL .and. lwork == lwork_min) then
       write(nfout,'(" lwork = ",i8)') lwork
       write(nfout,'(" lda, ndim = ",2i8)') lda, ndim
    end if

    call dsyev(JOBZ,UPLO,ndim,w1hw2,lda,eig,work_lapack,lwork,info)
!!$    write(*,*) 'lwork& info & suggested lwork = ',lwork,info,work_lapack(1)
    if(iprisubmat >= 2) then
       write(nfout,'(" info = ",i8)') info
    end if
!!$    lwork= int(work_lapack(1))+1
    deallocate(work_lapack)

    call tstatc0_end(id_sname)
  end subroutine dsyev_driver
#else
  subroutine dsyev_driver(ndim,eig,w1hw2)
    integer, intent(in):: ndim
    real(kind=DP), intent(out) ,dimension(ndim) :: eig
    real(kind=DP), intent(inout) ,dimension(ndim,ndim) :: w1hw2

    character(len=1) :: JOBZ,UPLO
!   integer, parameter :: lda=ndim
    integer :: lda
    integer :: lwork,liwork
    complex(kind=CMPLDP),allocatable,dimension(:) :: work
    real(kind=DP),allocatable,dimension(:) :: iwork
    integer :: info

    integer :: id_sname = -1
    call tstatc0_begin('dsyev_driver ', id_sname)

    lda=ndim
!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
    JOBZ = 'V'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
    UPLO = 'U'
! (define lwork)

    lwork  = 1 + ndim*(6+2*ndim)
    liwork = 3 + 5*ndim

    allocate(work(lwork)); work = 0.d0
    allocate(iwork(liwork)); iwork = 0.d0

#ifdef NEC_ITER_REG
    call FTRACE_REGION_BEGIN("DSYEVD")
#endif
    call dsyevd(JOBZ,UPLO,ndim,w1hw2,lda,eig,work,lwork,iwork,liwork,info)
#ifdef NEC_ITER_REG
    call FTRACE_REGION_END("DSYEVD")
#endif

    deallocate(work)
    deallocate(iwork)

    call tstatc0_end(id_sname)
  end subroutine dsyev_driver
#endif

#ifndef _USE_ZHEEVD_
  subroutine zheev_driver(ndim,eig,w1hw2)
   integer, intent(in):: ndim
   real(kind=DP), intent(out) ,dimension(ndim) :: eig
   real(kind=DP), intent(inout) ,dimension(ndim*2,ndim) :: w1hw2
   character(len=1) :: JOBZ,UPLO
!   integer, parameter :: lda=ndim
   integer :: lda
   integer :: i,j
! --> T. Yamasaki 2009/07/24 revised according to a Dr. Katagiri's report
!!$   integer :: lwork,lwork_min,rwork_size
   integer :: lwork_min,rwork_size
   integer, save :: lwork = 0
! <--
   integer, parameter :: lwork_huge=10**8
   real(kind=DP),allocatable,dimension(:) :: work_lapack
   real(kind=DP),allocatable,dimension(:) :: rwork_lapack
   integer :: info

   integer :: id_sname = -1
   call tstatc0_begin('zheev_driver ', id_sname)

      lda=ndim
!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
     JOBZ = 'V'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
     UPLO = 'U'
      lwork_min=max(1,2*ndim-1)
      rwork_size=max(3*ndim-2,1)
     if(lwork < lwork_min)  lwork=lwork_min
     if(lwork > lwork_huge) lwork=lwork_min

     allocate(work_lapack(lwork*2))
     allocate(rwork_lapack(rwork_size))

     call zheev(JOBZ,UPLO,ndim,w1hw2,lda,eig,work_lapack,lwork,rwork_lapack,info)
!!$     write(*,*) 'lwork& info & suggested lwork = ',lwork,info,work_lapack(1)
     lwork= int(work_lapack(1))+1

     deallocate(work_lapack)
     deallocate(rwork_lapack)

     call tstatc0_end(id_sname)
   end subroutine zheev_driver
#else
  subroutine zheev_driver(ndim,eig,w1hw2)
   integer, intent(in):: ndim
   real(kind=DP), intent(out) ,dimension(ndim) :: eig
   real(kind=DP), intent(inout) ,dimension(ndim*2,ndim) :: w1hw2
   character(len=1) :: JOBZ,UPLO
!   integer, parameter :: lda=ndim
   integer :: lda
   integer :: lwork,lrwork,liwork
   complex(kind=CMPLDP),allocatable,dimension(:) :: work
   real(kind=DP),allocatable,dimension(:) :: rwork
   real(kind=DP),allocatable,dimension(:) :: iwork
   integer :: info

   integer :: id_sname = -1
   call tstatc0_begin('zheev_driver ', id_sname)

      lda=ndim
!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
     JOBZ = 'V'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
     UPLO = 'U'
      lwork=ndim*(2+ndim)
      lrwork=1+ndim*(5+2*ndim)
      liwork=3+5*ndim

     allocate(work(lwork))
     allocate(rwork(lrwork))
     allocate(iwork(liwork))

#ifdef NEC_ITER_REG
    call FTRACE_REGION_BEGIN("ZHEEVD")
#endif
     call zheevd(JOBZ,UPLO,ndim,w1hw2,lda,eig,work,lwork,rwork,lrwork,iwork,liwork,info)
#ifdef NEC_ITER_REG
    call FTRACE_REGION_END("ZHEEVD")
#endif

     deallocate(work)
     deallocate(rwork)
     deallocate(iwork)

     call tstatc0_end(id_sname)
   end subroutine zheev_driver
#endif
! _NO_USE_ZHEEVD_


#ifndef _USE_ZHEGVXD_
  subroutine zhegvx_driver(ndim,eig,w1hw2,w1sw2)
   integer, intent(in):: ndim
   real(kind=DP), intent(out) ,dimension(ndim) :: eig
   real(kind=DP), intent(inout) ,dimension(ndim*2,ndim) :: w1hw2,w1sw2
   character(len=1) :: JOBZ,UPLO
!   integer, parameter :: lda=ndim
   integer :: lda
! --> T. Yamasaki 2009/07/24 revised according to a Dr. Katagiri's report
!!$   integer :: lwork,lwork_min,rwork_size
   integer :: lwork_min,rwork_size
   integer, save :: lwork = 0
! <--
   integer, parameter :: lwork_huge=10**8
   real(kind=DP),allocatable,dimension(:) :: work_lapack
   real(kind=DP),allocatable,dimension(:) :: rwork_lapack
   integer :: info

   real(kind=DP),allocatable,dimension(:,:) :: z
   integer,allocatable,dimension(:) :: ifail
   integer,allocatable,dimension(:) :: iwork
   integer :: ITYPE,ldb,IL,IU,ldz
   integer :: i,j
   integer :: meig
   real(kind=DP) :: VL,VU,ABSTOL
   character(len=1) :: RANGE

   integer :: id_sname = -1
   call tstatc0_begin('zhegvx_driver ', id_sname)

!write(90,*) 'fail ::    zhegvx_driver'

      ITYPE=1
      RANGE='A'
      ABSTOL=-1.d0

      lda=ndim
      ldb=ndim
      ldz=ndim
!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
     JOBZ = 'V'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
     UPLO = 'U'

      lwork=ndim*(2+ndim)
      lwork_min=max(1,2*ndim-1)
     if(lwork < lwork_min)  lwork=lwork_min
     if(lwork > lwork_huge) lwork=lwork_min

      rwork_size=7*ndim !max(3*ndim-2,1)

     allocate(work_lapack(lwork*2))
     allocate(rwork_lapack(rwork_size))

     allocate(z(ldz*2,ndim))
     allocate(ifail(ndim))
     allocate(iwork(5*ndim))

     call zhegvx(ITYPE,JOBZ,RANGE,UPLO,ndim,w1hw2,lda,w1sw2,ldb,VL,VU,IL,IU,ABSTOL &
                ,meig,eig,z,ldz,work_lapack,lwork,rwork_lapack,iwork,ifail,info)
     w1hw2=z
     lwork= int(work_lapack(1))+1

     deallocate(work_lapack)
     deallocate(rwork_lapack)

     deallocate(z)
     deallocate(ifail)
     deallocate(iwork)

     call tstatc0_end(id_sname)
   end subroutine zhegvx_driver
#else
  subroutine zhegvx_driver(ndim,eig,w1hw2,w1sw2)
   integer, intent(in):: ndim
   real(kind=DP), intent(out) ,dimension(ndim) :: eig
   real(kind=DP), intent(inout) ,dimension(ndim*2,ndim) :: w1hw2,w1sw2
   character(len=1) :: JOBZ,UPLO
!   integer, parameter :: lda=ndim
   integer :: lda
   integer :: lwork,lrwork,liwork
   complex(kind=CMPLDP),allocatable,dimension(:) :: work
   real(kind=DP),allocatable,dimension(:) :: rwork
   real(kind=DP),allocatable,dimension(:) :: iwork
   integer :: info
   integer :: meig

   complex(kind=CMPLDP),allocatable,dimension(:,:) :: z
   integer,allocatable,dimension(:) :: ifail
   integer :: ITYPE,ldb,IL,IU,ldz
   real(kind=DP) :: VL,VU,ABSTOL
   character(len=1) :: RANGE

   integer :: id_sname = -1
   call tstatc0_begin('zhegvx_driver ', id_sname)

!write(90,*) 'zhegvx_driver'

      ITYPE=1
      RANGE='A'
      ABSTOL=-1.d0

      lda=ndim
      ldb=ndim
      ldz=ndim
!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
     JOBZ = 'V'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
     UPLO = 'U'
      lwork=ndim*(2+ndim)
      lrwork=1+ndim*(5+2*ndim)
      liwork=5*ndim !3+5*ndim

     allocate(work(lwork))
     allocate(rwork(lrwork))
     allocate(iwork(liwork))

     allocate(z(ldz,ndim))
     allocate(ifail(ndim))

#ifdef NEC_ITER_REG
    call FTRACE_REGION_BEGIN("ZHEGVXD")
#endif
!     call zheevd(JOBZ,UPLO,ndim,w1hw2,lda,eig,work,lwork,rwork,lrwork,iwork,liwork,info)
     call zhegvx(ITYPE,JOBZ,RANGE,UPLO,ndim,w1hw2,lda,w1sw2,ldb,VL,VU,IL,IU,ABSTOL &
                ,meig,eig,z,ldz,work,lwork,rwork,iwork,ifail,info)
     w1hw2=z
#ifdef NEC_ITER_REG
    call FTRACE_REGION_END("ZHEGVXD")
#endif

     deallocate(work)
     deallocate(rwork)
     deallocate(iwork)

     deallocate(z)
     deallocate(ifail)

     call tstatc0_end(id_sname)
   end subroutine zhegvx_driver
#endif


! <-- T.Kokubo & D.Fukata, Feb. 2010
  subroutine set_hmat(kimg_t,meg,is,ie,w1hw2,zmat)
    integer, intent(in) :: kimg_t,meg,is,ie
    real(kind=DP), intent(inout) :: w1hw2(*)
    real(kind=DP), intent(out) :: zmat(meg*kimg_t,meg)
    integer, parameter         :: msize_allreduce_max = 1000000

    real(kind=DP), allocatable :: w1hw2_mpi(:), w1hw2_mpi2(:)
    integer :: i,j, msize_partitioned
    integer :: m,msize
    integer :: id_sname = -1
    call tstatc0_begin('set_hmat ', id_sname)

    msize = is*(ie-is+1) + (ie-is+1)*(ie-is)/2

#ifdef NEC_ITER_REG
    call FTRACE_REGION_BEGIN("SET_HMAT")
#endif
!! (spread sum of w1hw2)
    if(npes > 1) then
       allocate(w1hw2_mpi(msize*kimg_t))
       if(msize*kimg_t <= msize_allreduce_max) then
          if(iprisubmat >= 2) write(nfout,'(" msize*kimg_t = ",i12,",  msize_allreduce_max = ",i12," <<set_hmat>>")') &
         & msize*kimg_t, msize_allreduce_max
          call mpi_allreduce(w1hw2,w1hw2_mpi,msize*kimg_t,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI
          w1hw2(1:msize*kimg_t) = w1hw2_mpi         ! MPI
       else
          if(iprisubmat >= 2)  &
               & write(nfout,'(" msize*kimg_t = ",i12," > msize_allreduce_max = ",i12," <<set_hmat>>")') &
               & msize*kimg_t, msize_allreduce_max
          allocate(w1hw2_mpi2(msize_allreduce_max))
          do j = 1, ceiling(dble(msize*kimg_t)/msize_allreduce_max)
             msize_partitioned = min(msize_allreduce_max, msize*kimg_t-(j-1)*msize_allreduce_max)
             do m = 1, msize_partitioned
                i = m + (j-1)*msize_allreduce_max
                w1hw2_mpi2(m) = w1hw2(i)
             end do
             call mpi_allreduce(MPI_IN_PLACE,w1hw2_mpi2,msize_partitioned,mpi_double_precision &
                  &                             ,mpi_sum,mpi_k_world(myrank_k),ierr)
             do m = 1, msize_partitioned
                i = m + (j-1)*msize_allreduce_max
                w1hw2_mpi(i) = w1hw2_mpi2(m)
             end do
          end do
          w1hw2(1:msize*kimg_t) = w1hw2_mpi
          deallocate(w1hw2_mpi2)
       end if
       deallocate(w1hw2_mpi)
    end if
    if(kimg_t==1) then
!CDIR PARALLEL DO PRIVATE (i,m)
       do j=is,ie
          do i=1,j
             m = j*(j-1)/2-is*(is-1)/2+i
             zmat(i,j) =  w1hw2(m)
          end do
       end do
    else
!CDIR PARALLEL DO PRIVATE (i,m)
       do j=is,ie
          do i=1,j
             m = j*(j-1)/2-is*(is-1)/2+i
             zmat(2*i-1,j) =  w1hw2(2*m-1)
             zmat(2*i  ,j) = -w1hw2(2*m)
          end do
       end do
    end if
#ifdef NEC_ITER_REG
    call FTRACE_REGION_END("SET_HMAT")
#endif
    call tstatc0_end(id_sname)
  end subroutine set_hmat
#ifdef _USE_SCALAPACK_
  subroutine set_nprow_npcol(nprow,npcol)
    integer, intent(inout) :: nprow,npcol
    logical :: cflag
    integer :: j,n

    cflag=.true.

    if(nprow*npcol /= nrank_e) then
       nprow=1
       npcol=1
       j = nrank_e

 120   continue
       if(mod(j,2) == 0) then
         j=j/2
         n=2
         call prod(cflag,n,nprow,npcol)
         goto 120
       else if(mod(j,3) == 0) then
         j=j/3
         n=3
         call prod(cflag,n,nprow,npcol)
         goto 120
       else if(mod(j,5) == 0) then
         j=j/5
         n=5
         call prod(cflag,n,nprow,npcol)
         goto 120
       else if(j == 1) then
         n=1
       else
         n=j
       end if
       call prod(cflag,n,nprow,npcol)
       write(*,'("set nprow,npcol=",2i5)') nprow,npcol
       if(printable) write(nfout,'("set nprow,npcol=",2i5)') nprow,npcol
    else
       if(iprisubmat >= 2) write(nfout,'("nprow,npcol=",2i5)') nprow,npcol
    end if
  end subroutine set_nprow_npcol
  subroutine prod(cflag,n,nr,nc)
    logical,intent(inout) :: cflag
    integer,intent(in)    :: n
    integer,intent(inout) :: nr,nc
    if(cflag) then
       nr = nr * n
    else
       nc = nc * n
    endif
    cflag= .not.cflag
  end subroutine prod
#endif


#ifdef _USE_SCALAPACK_
  subroutine set_amat(lda,kimg_t,occ,amat,ndim,nprow,npcol,nb,usermap,is,ie,w1hw2)
    integer, intent(in) :: lda,kimg_t,occ
    real(kind=DP), intent(out) :: amat(lda*kimg_t,occ)
    integer, intent(in) :: ndim,nprow,npcol,nb
#ifdef SINGLE_CONTEXT
    integer, intent(in) :: usermap(nprow,npcol)
#else
    integer, intent(in) :: usermap(nprow,npcol,0:nrank_k-1)
#endif
    integer, intent(in) :: is,ie
    real(kind=DP), intent(inout) :: w1hw2(*)

    integer :: i,j
    integer :: jj,ipcol,ilcol,ixcol
    integer :: ii,iprow,ilrow,ixrow
    integer :: m,msize
    real(kind=DP), allocatable :: w1hw2_mpi(:)
    integer :: id_sname = -1
    call tstatc0_begin('set_amat ', id_sname)

    msize = is*(ie-is+1) + (ie-is+1)*(ie-is)/2
#ifdef NEC_ITER_REG
    call FTRACE_REGION_BEGIN("SET_AMAT")
#endif

!! (spread sum of w1hw2)
    if(npes > 1) then
       if(iprisubmat >= 2) write(nfout,'(" msize*kimg_t = ",i12," <<set_amat>>")') msize*kimg_t
       allocate(w1hw2_mpi(msize*kimg_t))
       call mpi_allreduce(w1hw2,w1hw2_mpi,msize*kimg_t,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
       w1hw2(1:msize*kimg_t) = w1hw2_mpi
       deallocate(w1hw2_mpi)
    end if

    if(kimg_t==1) then
       !!$!CDIR PARALLEL DO PRIVATE ( ipcol,ilcol,ixcol,jj,i,iprow,ilrow,ii,m)
       do j=is,ie
          ipcol = mod((j-1)/nb,npcol)
          ilcol = (j-1)/(npcol*nb)
          ixcol = mod(j-1,nb)+1
          jj = ilcol*nb+ixcol
          do i=1,j
             iprow = mod((i-1)/nb,nprow)
             ilrow = (i-1)/(nprow*nb)
             ixrow = mod(i-1,nb)+1
             ii = ilrow*nb+ixrow
#ifdef SINGLE_CONTEXT
             if(usermap(iprow+1,ipcol+1) /= mype) cycle
#else
             if(usermap(iprow+1,ipcol+1,myrank_k) /= mype) cycle
#endif
             m = j*(j-1)/2-is*(is-1)/2+i
             amat(ii,jj) =  w1hw2(m)
          end do
       end do
    else
       !!$!CDIR PARALLEL DO PRIVATE ( ipcol,ilcol,ixcol,jj,i,iprow,ilrow,ii,m)
       do j=is,ie
          ipcol = mod((j-1)/nb,npcol)
          ilcol = (j-1)/(npcol*nb)
          ixcol = mod(j-1,nb)+1
          jj = ilcol*nb+ixcol
          do i=1,j
             iprow = mod((i-1)/nb,nprow)
             ilrow = (i-1)/(nprow*nb)
             ixrow = mod(i-1,nb)+1
             ii = ilrow*nb+ixrow
#ifdef SINGLE_CONTEXT
             if(usermap(iprow+1,ipcol+1) /= mype) cycle
#else
             if(usermap(iprow+1,ipcol+1,myrank_k) /= mype) cycle
#endif
             m = j*(j-1)/2-is*(is-1)/2+i
             amat(2*ii-1,jj) =  w1hw2(2*m-1)
             amat(2*ii  ,jj) = -w1hw2(2*m)
          end do
       end do
    end if
#ifdef NEC_ITER_REG
    call FTRACE_REGION_END("SET_AMAT")
#endif
    call tstatc0_end(id_sname)
  end subroutine set_amat
#endif

#ifdef _USE_SCALAPACK_
  subroutine put_zmat(kimg_t,lda,occ,zmat,ndim,nprow,npcol,nb,usermap,is,ie,zz)
    integer, intent(in) :: kimg_t,lda,occ
    real(kind=DP), intent(in) :: zmat(lda*kimg_t,occ)
    integer, intent(in) :: ndim,nprow,npcol,nb
#ifdef SINGLE_CONTEXT
    integer, intent(in) :: usermap(nprow,npcol)
#else
    integer, intent(in) :: usermap(nprow,npcol,0:nrank_k-1)
#endif
    integer, intent(in) :: is,ie
    real(kind=DP), intent(out) :: zz(ndim*kimg_t,is:ie)

    integer :: i,j
    integer :: jj,ipcol,ilcol,ixcol
    integer :: ii,iprow,ilrow,ixrow
    integer :: jsize
    real(kind=DP), allocatable :: zz_mpi(:,:)
    integer :: id_sname = -1
    call tstatc0_begin('put_zmat ', id_sname)

    jsize = ie-is+1

#ifdef NEC_ITER_REG
    call FTRACE_REGION_BEGIN("PUT_ZMAT")
#endif
    zz = 0.d0
    if(kimg_t==1) then
       do j=is,ie
          ipcol = mod((j-1)/nb,npcol)
          ilcol = (j-1)/(npcol*nb)
          ixcol = mod(j-1,nb)+1
          jj = ilcol*nb+ixcol
          do i=1,ndim
             iprow = mod((i-1)/nb,nprow)
             ilrow = (i-1)/(nprow*nb)
             ixrow = mod(i-1,nb)+1
             ii = ilrow*nb+ixrow
#ifdef SINGLE_CONTEXT
             if(usermap(iprow+1,ipcol+1) /= mype) cycle
#else
             if(usermap(iprow+1,ipcol+1,myrank_k) /= mype) cycle
#endif
             zz(i,j) = zmat(ii,jj)
          end do
       end do
    else
       do j=is,ie
          ipcol = mod((j-1)/nb,npcol)
          ilcol = (j-1)/(npcol*nb)
          ixcol = mod(j-1,nb)+1
          jj = ilcol*nb+ixcol
          do i=1,ndim
             iprow = mod((i-1)/nb,nprow)
             ilrow = (i-1)/(nprow*nb)
             ixrow = mod(i-1,nb)+1
             ii = ilrow*nb+ixrow
#ifdef SINGLE_CONTEXT
             if(usermap(iprow+1,ipcol+1) /= mype) cycle
#else
             if(usermap(iprow+1,ipcol+1,myrank_k) /= mype) cycle
#endif
             zz(2*i-1,j) = zmat(2*ii-1,jj)
             zz(2*i  ,j) = zmat(2*ii  ,jj)
          end do
       end do
    end if
    if(npes > 1) then
       allocate(zz_mpi(ndim*kimg_t,is:ie))
       call mpi_allreduce(zz,zz_mpi,ndim*kimg_t*jsize,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI
       zz = zz_mpi         ! MPI
       deallocate(zz_mpi)
    end if
#ifdef NEC_ITER_REG
    call FTRACE_REGION_END("PUT_ZMAT")
#endif
    call tstatc0_end(id_sname)
  end subroutine put_zmat
#endif

! ------------------------------
#ifdef _USE_SCALAPACK_
  subroutine scalapack_setup(ndim,lwork1,lrwork1,liwork1,lwork2,lrwork2,liwork2,lda,occ,ictxt,myrow,mycol,desca,descz,usermap)
    integer, intent(in) :: ndim
#ifdef SINGLE_CONTEXT
    integer, intent(out) :: lwork1,lrwork1,liwork1,lda,occ,ictxt,myrow,mycol
#else
    integer, intent(out) :: lwork1,lrwork1,liwork1,lda,occ,ictxt(0:nrank_k-1),myrow,mycol
#endif
    integer, intent(out) :: lwork2,lrwork2,liwork2
    integer, dimension(9), intent(out) :: desca,descz
#ifdef SINGLE_CONTEXT
    integer, intent(inout) :: usermap(nprow,npcol)
#else
    integer, intent(inout) :: usermap(nprow,npcol,0:nrank_k-1)
#endif

    integer :: i,j
#ifndef SINGLE_CONTEXT
    integer :: ik, pid, blacs_pnum
    integer, save ::  tmpctxt
#endif
    integer :: nb
    integer :: np0,nq0,np,nq,iarow,iacol
    integer :: iroffa, icoffa, nrc, ldc
    integer :: sizemqrleft, qrmem, trilwmin
    integer :: iam,nprocs,info
    integer, external :: indxg2l, indxg2p, numroc
    integer, parameter :: BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1, &
     &                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
     &                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9
    integer :: id_sname = -1
    call tstatc0_begin('scalapack_setup ', id_sname)

    nb = block_size

    if(mod(ndim,nb)>0) then
       lda = ndim/nb +1
    else
       lda = ndim/nb
    end if
    occ = lda
    if(mod(lda,nprow)>0) then
       lda = lda/nprow+1
    else
       lda = lda/nprow
    end if
    if(mod(occ,npcol)>0) then
       occ = occ/npcol+1
    else
       occ = occ/npcol
    end if
    lda = lda*nb
    occ = occ*nb
    if(iprisubmat>=2) write(nfout,'("lda,occ=",2i5)') lda,occ

#ifdef SINGLE_CONTEXT
    do j=1,npcol
       do i=1,nprow
          usermap(i,j) = myrank_k*nrank_e + (j-1)*nprow + i - 1
       end do
    end do
    if(iprisubmat>=1) then
       write(nfout,'("USERMPAP")')
       do i=1,nprow
          write(nfout,'(1x,i3)') (usermap(i,j),j=1,npcol)
       end do
    end if
    call blacs_setup(iam,nprocs)
    call blacs_get(-1,0,ictxt)
    call blacs_gridmap(ictxt,usermap,nprow,nprow,npcol)
    call blacs_gridinfo(ictxt,nprow,npcol,myrow,mycol)
#else
    call blacs_setup(iam,nprocs)
    call blacs_get(0,0,tmpctxt)
    call blacs_gridinit(tmpctxt,'R',nprocs,1)
    do ik=0,nrank_k-1
       pid = ik*nrank_e
       do j=1,npcol
          do i=1,nprow
             usermap(i,j,ik) = blacs_pnum(tmpctxt,pid,0)
             pid = pid + 1
          end do
       end do
       call blacs_get(0,0,ictxt(ik))
       call blacs_gridmap(ictxt(ik),usermap(1,1,ik),nprow,nprow,npcol)
    end do
    call blacs_gridinfo(ictxt(myrank_k),nprow,npcol,myrow,mycol)
    if(iprisubmat>=1) then
       do ik=0,nrank_k-1
          write(nfout,'("USERMPAP: ik=",i5)') ik
          do i=1,nprow
             write(nfout,'(128(1x,i3))') (usermap(i,j,ik),j=1,npcol)
          end do
       end do
    end if
#endif
    if(myrow == -1) then
       if(printable) write(nfout,*) 'BLACS init failed.'
       stop 'BLACS init failed.'
    else
       if(iprisubmat>=2) write(nfout,'("nprow,npcol,myrow,mycol=",4i5)') nprow,npcol,myrow,mycol
    end if
#ifdef SINGLE_CONTEXT
    call descinit(desca,ndim,ndim,nb,nb,0,0,ictxt,lda,info)
    if(iprisubmat>=2) write(nfout,*) 'descinit for amat: info=',info
    call descinit(descz,ndim,ndim,nb,nb,0,0,ictxt,lda,info)
    if(iprisubmat>=2) write(nfout,*) 'descinit for zmat: info=',info
#else
    call descinit(desca,ndim,ndim,nb,nb,0,0,ictxt(myrank_k),lda,info)
    if(iprisubmat>=2) write(nfout,*) 'descinit for amat: info=',info
    call descinit(descz,ndim,ndim,nb,nb,0,0,ictxt(myrank_k),lda,info)
    if(iprisubmat>=2) write(nfout,*) 'descinit for zmat: info=',info
#endif

!   if(kimg_t==1) then
       ! setup for pdsyev
       if(method_scalapack==HOUSEHOLDER) then
          IROFFA = MOD( 0, nb )
          ICOFFA = MOD( 0, nb )
          iarow = INDXG2P( 1, nb, myrow, desca( RSRC_ ), nprow )
          iacol = INDXG2P( 1, nb, mycol, desca( RSRC_ ), npcol )
          np = NUMROC( ndim+iroffa, nb, myrow, iarow, nprow )
          nq = NUMROC( ndim+icoffa, nb, mycol, iacol, npcol )
          sizemqrleft = max( (nb*(nb-1))/2, (np+nq)*nb ) + nb*nb
          if(printable) write(nfout,*) 'np,nq=',np,nq
          if(printable) write(nfout,*) 'sizemqrleft=',sizemqrleft
          np0 = NUMROC( max(ndim,nb,2), nb, 0, 0, nprow )
          nq0 = NUMROC( max(ndim,nb,2), nb, 0, 0, npcol )
          nrc = NUMROC( ndim, nb, myrow, 0, nprocs)
          ldc = max( 1, nrc )
          if(printable) write(nfout,*) 'np0,nq0=',np0,nq0
          if(printable) write(nfout,*) 'nrc,ldc=',nrc,ldc
          qrmem = 2*ndim-2
          if(printable) write(nfout,'("qrmem,=",i9)') qrmem
          lwork1 = 5*ndim + ndim*ldc + max( sizemqrleft, qrmem ) + 1
          if(printable) write(nfout,'("lwork1,=",i9)') lwork1
          lrwork1 = 0
          liwork1 = 0
       else
          iarow = INDXG2P( 1, nb, myrow, desca( RSRC_ ), nprow )
          iacol = INDXG2P( 1, nb, mycol, desca( CSRC_ ), npcol )
          np0 = NUMROC( ndim, nb, myrow, iarow, nprow )
          nq0 = NUMROC( ndim, nb, mycol, iacol, npcol )
          if(printable) then
             write(nfout,*) 'iarow,iacol=',iarow,iacol
             write(nfout,*) 'np0,nq0=',np0,nq0
             write(nfout,*) 'nb_a=', desca( NB_ )
          end if
          trilwmin = 3*ndim + max( nb*( np0+1 ), 3*nb )
          lwork1 = max( 1+6*ndim+2*np0*nq0, trilwmin ) + 2*ndim
          liwork1 = 7*ndim + 8*npcol + 2
          if(printable) write(nfout,'("lwork1,liwork1,=",2i9)') lwork1,liwork1
          lrwork1 = 0
       end if
!   else
       ! setup for pzheev
       if(method_scalapack==HOUSEHOLDER) then
          np0 = NUMROC( max(ndim,nb,2), nb, 0, 0, nprow )
          nq0 = NUMROC( max(ndim,nb,2), nb, 0, 0, npcol )
          if(printable) write(nfout,*) 'np0,nq0=',np0,nq0
          lwork2 = ndim*(ndim+3) + ( np0 + nq0 + nb ) * nb
          lrwork2 = 2*ndim + 2*ndim-2
          if(printable) write(nfout,'("lwork2,lrwork2,=",3i9)') lwork2,lrwork2
          liwork2 = 0
       else
          iarow = INDXG2P( 1, nb, myrow, desca( RSRC_ ), nprow )
          iacol = INDXG2P( 1, nb, mycol, desca( CSRC_ ), npcol )
          np0 = NUMROC( ndim, nb, myrow, iarow, nprow )
          nq0 = NUMROC( ndim, nb, mycol, iacol, npcol )
          nq = NUMROC( MAX( ndim, nb, 2 ), nb, 0, 0, npcol )
          if(printable) then
             write(nfout,*) 'iarow,iacol=',iarow,iacol
             write(nfout,*) 'np0,nq0,nq=',np0,nq0,nq
             write(nfout,*) 'nb_a=', desca( NB_ )
          end if
          lwork2 = ndim + ( np0 + nq + nb ) * nb
          lrwork2 = 1 + 9*ndim + 3*np0*nq0
          liwork2 = 7*ndim + 8*npcol + 2
          if(printable) write(nfout,'("lwork2,lrwork2,liwork2,=",3i9)') lwork2,lrwork2,liwork2
       end if
!   end if

    call tstatc0_end(id_sname)
  end subroutine scalapack_setup
#endif

#ifdef _USE_SCALAPACK_
  subroutine pdsyev_driver(ndim,eig,desca,descz,lda,occ,amat,zmat,lwork,liwork)
    integer, intent(in):: ndim
    real(kind=DP), intent(out) ,dimension(ndim) :: eig
    integer, dimension(9), intent(in) :: desca,descz
    integer, intent(in) :: lda,occ
    real(kind=DP), intent(inout) ,dimension(lda,occ) :: amat
    real(kind=DP), intent(out) ,dimension(lda,occ) :: zmat
    integer, intent(in):: lwork,liwork

    character(len=1) :: JOBZ,UPLO
    real(kind=DP),allocatable,dimension(:) :: work
    integer,allocatable,dimension(:) :: iwork
    integer :: info
    integer :: id_sname = -1
    call tstatc0_begin('pdsyev_driver ', id_sname)

!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
    JOBZ = 'V'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
    UPLO = 'U'

    allocate(work(lwork))
    if(method_scalapack==DIVIDEandCONQUER) allocate(iwork(liwork))

    if(method_scalapack==HOUSEHOLDER) then
#ifdef NEC_ITER_REG
       call FTRACE_REGION_BEGIN("PDSYEV")
#endif
       call pdsyev(JOBZ,UPLO,ndim,amat,1,1,desca,eig,zmat,1,1 &
                & ,descz,work,lwork,info)
#ifdef NEC_ITER_REG
       call FTRACE_REGION_END("PDSYEV")
#endif
       if(iprisubmat>=2) write(nfout,*) 'pdsyev: info=',info
    else
#ifdef NEC_ITER_REG
       call FTRACE_REGION_BEGIN("PDSYEVD")
#endif
       call pdsyevd(JOBZ,UPLO,ndim,amat,1,1,desca,eig,zmat,1,1 &
                 & ,descz,work,lwork,iwork,liwork,info)
#ifdef NEC_ITER_REG
       call FTRACE_REGION_END("PDSYEVD")
#endif
       if(iprisubmat>=2) write(nfout,*) 'pdsyevd: info=',info
    end if


    deallocate(work)
    if(method_scalapack==DIVIDEandCONQUER) deallocate(iwork)

    call tstatc0_end(id_sname)
  end subroutine pdsyev_driver

  subroutine pzheev_driver(ndim,eig,desca,descz,lda,occ,amat,zmat,lwork,lrwork,liwork)
    integer, intent(in):: ndim
    real(kind=DP), intent(out) ,dimension(ndim) :: eig
    integer, dimension(9), intent(in) :: desca,descz
    integer, intent(in) :: lda,occ
    real(kind=DP), intent(inout) ,dimension(lda*2,occ) :: amat
    real(kind=DP), intent(out) ,dimension(lda*2,occ) :: zmat
    integer, intent(in):: lwork,lrwork,liwork

    character(len=1) :: JOBZ,UPLO
    complex(kind=CMPLDP),allocatable,dimension(:) :: work
    real(kind=DP),allocatable,dimension(:) :: rwork
    integer,allocatable,dimension(:) :: iwork
    integer :: info

    integer :: id_sname = -1
    call tstatc0_begin('pzheev_driver ', id_sname)

!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
    JOBZ = 'V'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
    UPLO = 'U'

    allocate(work(lwork))
    allocate(rwork(lrwork))
    if(method_scalapack==DIVIDEandCONQUER) allocate(iwork(liwork))

    if(method_scalapack==HOUSEHOLDER) then
#ifdef NEC_ITER_REG
       call FTRACE_REGION_BEGIN("PZHEEV")
#endif
       call pzheev(JOBZ,UPLO,ndim,amat,1,1,desca,eig,zmat,1,1 &
                & ,descz,work,lwork,rwork,lrwork,info)
#ifdef NEC_ITER_REG
       call FTRACE_REGION_END("PZHEEV")
#endif
       if(iprisubmat>=2) write(nfout,*) 'pzheev: info=',info
    else
#ifdef NEC_ITER_REG
       call FTRACE_REGION_BEGIN("PZHEEVD")
#endif
       call pzheevd(JOBZ,UPLO,ndim,amat,1,1,desca,eig,zmat,1,1 &
                 & ,descz,work,lwork,rwork,lrwork,iwork,liwork,info)
#ifdef NEC_ITER_REG
       call FTRACE_REGION_END("PZHEEVD")
#endif
       if(iprisubmat>=2) write(nfout,*) 'pzheevd: info=',info
    end if

    deallocate(work)
    deallocate(rwork)
    if(method_scalapack==DIVIDEandCONQUER) deallocate(iwork)

    call tstatc0_end(id_sname)
  end subroutine pzheev_driver
#endif
! _USE_SCALAPACK_



   subroutine m_ESsubmat_utransform_wf(imode,ik,zaj_l,update_fsr,wf_mode)
     integer, intent(in) :: imode,ik
     real(kind=DP), dimension(kg1,np_e,ista_k:iend_k,kimg),intent(inout) :: zaj_l
     logical, intent(in) :: update_fsr
     logical, intent(in) :: wf_mode
     real(kind=DP), allocatable, dimension(:,:,:) :: zat_l,zat_wk
     real(kind=DP), allocatable, dimension(:,:,:) :: zahloc_l,zahloc_wk
     real(kind=DP), allocatable, dimension(:,:) :: zz
     integer :: jp,is,ie,nblk1,ibsize2,ibsize,kimg_t
     real(kind=DP),allocatable,dimension(:,:) :: mat
     if(.not.allocated(unimat)) return
     if(disable_utransform) return
#ifdef _ODD_BOUNDARY_
     if(mod(np_g1k(ik),2) == 0) then
        np_g1k_x = np_g1k(ik) + 1
     else
        np_g1k_x = np_g1k(ik)
     end if
#else
     np_g1k_x = np_g1k(ik)
#endif
     np_g1k_x = maxval(np_g1k)

     if(k_symmetry(ik) == GAMMA) then
        kimg_t = 1
     else
        kimg_t = kimg
     end if

     if(nblocksize_submat_is_given) then
       ibsize = nblocksize_submat
     else
       ibsize = nb_submat_default
     endif

!     if(imode==DIRECT .and. wf_mode)then
!        zaj_l(1:kg1,1:np_e,ik,1:kimg) = zaj_l_buf(1:kg1,1:np_e,ik,1:kimg)
!        if(update_fsr) call m_ES_betar_dot_WFs_4_each_k(nfout,ik)
!        return
!     endif
!     if(imode==INVERSE .and. wf_mode) then
!        zaj_l_buf(1:kg1,1:np_e,ik,1:kimg) = zaj_l(1:kg1,1:np_e,ik,1:kimg)
!     endif

     allocate(zat_l(np_g1k_x,neg,kimg));zat_l=0.d0
     allocate(zahloc_l(np_g1k_x,neg,kimg))
     call m_ES_W_transpose_r(.false.,ista_k,iend_k,ik,zaj_l,zat_l)
     allocate(zat_wk(np_g1k_x,neg,kimg));zat_wk=0.d0
     allocate(zahloc_wk(np_g1k_x,neg,kimg))
     allocate(mat(neg*kimg_t,neg))
     if(imode==DIRECT)then
        mat(:,:)=unimat(:,:)
     else
        mat(:,:)=unimat_h(:,:)
     endif
! <-- T.Kokubo & D.Fukata, Feb. 2010
     PART2: do jp = 1,npart2
        is = isp2(jp)
        ie = iep2(jp)
        nblk1 = ie-is+1
#ifdef SX
        ibsize2 = neg
#else
        ibsize2 = ibsize
#endif
#ifdef _USE_SCALAPACK_
        if(sw_scalapack == ON) then
           allocate(zz(neg*kimg_t,nblk1))
           call put_zmat(kimg_t,lda,occ,mat,neg,nprow,npcol,block_size,usermap,is,ie,zz)
           if(kimg == 1) then
              call subspace_rotation_real(is,ie,zat_wk,zz,nblk1,ibsize2,ik,neg,zat_l,zahloc_wk,zahloc_l)
           else
              call subspace_rotation_imag(is,ie,kimg_t,zat_wk,zz,nblk1,ibsize2,ik,neg,zat_l,zahloc_wk,zahloc_l)
           endif
           deallocate(zz)
        else
#endif
           if(kimg == 1) then
              call subspace_rotation_real(is,ie,zat_wk,mat(1,is),nblk1,ibsize2,ik,neg,zat_l,zahloc_wk,zahloc_l)
           else
              call subspace_rotation_imag(is,ie,kimg_t,zat_wk,mat(1,is),nblk1,ibsize2,ik,neg,zat_l,zahloc_wk,zahloc_l)
           endif
#ifdef _USE_SCALAPACK_
        end if
#endif
     end do PART2
     call m_ES_W_transpose_back_r(.false.,ista_k,iend_k,ik,zaj_l,zat_wk)            ! zat_l -> zaj_l
     deallocate(zat_l)
     deallocate(zat_wk)
     if(update_fsr) call m_ES_betar_dot_WFs_4_each_k(nfout,ik)
   end subroutine m_ESsubmat_utransform_wf


end module m_ES_WF_by_submat
