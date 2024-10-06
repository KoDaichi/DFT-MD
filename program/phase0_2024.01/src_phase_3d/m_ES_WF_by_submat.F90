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
#define SINGLE_CONTEXT

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

#ifdef MPI_FFTW
  use, intrinsic :: iso_c_binding
#endif

  use m_Const_Parameters,    only : DP,SP,DIRECT,ON,OFF,SCF, SmallestPositiveNumber,GAMMA &
       &                          , HOUSEHOLDER, DIVIDEandCONQUER, ELECTRON, CMPLDP, OLD, INVERSE, EXECUT
  use m_Parallelization,     only : MPI_CommGroup &
       &                          , myrank_e,myrank_k,map_e,map_k,ista_e,iend_e,istep_e &
       &                          , ista_k,iend_k,np_g1k,mpi_k_world,map_z &
       &                          , np_e,npes, nrank_k, nrank_e,mype &
       &                          , np_fs, mpi_ge_world, mype_conf &
       &                          , nrank_conf, mpi_kr_world &
       &                          , ista_spin, iend_spin, nrank_s, myrank_spin
  use m_Control_Parameters,  only : nspin,iprisubmat,ldiag,kimg,neg,af &
       &                          , submat_critical_ratio,printable &
#ifdef USE_EIGENLIB
       &                          , m_CtrlP_set_submat, sw_scalapack, sw_eigen_exa &
#else
       &                          , m_CtrlP_set_submat, sw_scalapack &
#endif
       &                          , method_scalapack, block_size, nprow, npcol &
       &                          , msize_submat, sw_hybrid_functional, sw_fef &
       &                          , nblocksize_submat_is_given, nblocksize_submat, nb_submat_default &
       &                          , sw_gep &
       &                          , submat_uncalled &
       &                          , sw_serial_fft                    &
       &                          , nb_mgs_default &
       &                          , nblocksize_submat_latter_is_given, nblocksize_submat_latter &
       &                          , nblocksize_mgs_is_given                    &
       &                          , sw_keep_hloc_phi &
       &                          , m_CtrlP_get_conftag &
#ifdef MPI_FFTW
       &                          , sw_mpi_fftw &
#endif
#ifdef SAVE_FFT_TIMES
       &                          , nblocksize_mgs, divide_square, sw_save_fft
#else
       &                          , nblocksize_mgs, divide_square
#endif
  use m_Files,               only : nfout
  use m_Timing,              only : tstatc0_begin, tstatc0_end
  use m_FFT,                 only : nfft,fft_box_size_WF
#ifdef MPI_FFTW
  use m_FFT,                 only : bfft_mpifftw, m_FFT_Vlocal_W_mpifftw3d, m_FFT_Direct_MPI_FFTW, afft_mpifftw &
       &                          , bfft_mpifftw_kimg1, afft_mpifftw_kimg1, afft_mpifftw_vlocal
#endif
  use m_Kpoints,             only : kv3,vkxyz, k_symmetry
  use m_PlaneWaveBasisSet,   only : kg1, iba, igf, nbase, m_pwBS_kinetic_energies
  use m_Electronic_Structure,only : zaj_l, zaj_l_buf,neordr,           nrvf_ordr &
#ifdef SAVE_FFT_TIMES
       &                          , status_saved_phifftr &
#endif
       &                          , eko_l, vnlph_l, vtau_l
  use m_ES_nonlocal,         only : sc ,ss ,qc                                        &
 &                                , m_ES_dealloc_scss_etc
  use m_ES_ortho,            only : np_g1k_x, WSW_for_overlap_t_g
  use m_Electronic_Structure,only : zaj_ball, zah_ball         &
       &                          , fsr_l_2d,fsi_l_2d,fsr_l,fsi_l  &
       &                          , fsr_ball, fsi_ball &
       &                          , vlhxc_l                                &
       &                          , m_ES_Vlocal_in_Rspace_3D                          &
#ifdef MPI_FFTW
       &                          , m_ES_Vlocal_in_Rspace_mpifftw3d                   &
       &                          , m_ES_WF_in_Rspace_mpifftw                         &
#endif
       &                          , m_ES_WF_in_Rspace_3D                              &
       &                          , m_ES_wd_zaj_small_portion_3D                      &
       &                          , m_ES_wd_eko_3D                                    &
       &                          , nblocksize_mgs_default                            &
       &                          , m_ES_WF_2D                                        &
       &                          , hlocphi_l
  use m_ES_nonlocal,         only : sc_l                                              &
 &                                , m_ES_Vnonlocal_W_3D                               &
 &                                , m_ES_betar_dot_WFs_4_each_k_3D                    &
 &                                , m_ES_alloc_scss_etc_3D
  use m_PseudoPotential,     only : nlmtt, modnrm
  use m_FFT,                 only : m_FFT_Direct_3D, m_FFT_Direct_XYZ_3D              &
#ifdef FFT_3D_DIVISION
                                  , m_FFT_Direct_3DIV_3D, m_FFT_Vlocal_W_3DIV_3D      &
#endif
       &                          , m_FFT_Vlocal_W_3D
  use m_Parallelization,   only : nel_fft_z, nel_fft_y, nel_fft_x &
       &                        , fft_X_x_nel, fft_X_y_nel, fft_X_z_nel   &
       &                        , mp_fft_x, xyz_fft_x &
       &                        , myrank_g, nrank_g   &
       &                        , nis_e, nie_e, nel_e &
       &                        , ista_g1k, iend_g1k, mp_g1k, nel_g1k  &
       &                        , mpi_ke_world,  mpi_kg_world, ista_g1k, iend_g1k &
       &                        , nis_g1k, nie_g1k, neg_g, neg_g_all &
       &                        , nel_eg, nis_eg, nie_eg, neg_gg, neg_gg_all, mp_e
!#ifdef _USE_SCALAPACK_
#ifdef EIGEN_TRANS
  use m_Parallelization,     only : make_index_band_for_scalapack                &
       &                          , neg_col, nel_col, nis_col, nie_col, scl_col  &
       &                          , neg_row, nel_row, nis_row, nie_row, scl_row  &
       &                          , scl_comm_max, scl_comm_max_r                 &
       &                          , scl_comm_rank, scl_comm_rank_r               &
       &                          , scl_comm_rno,  scl_comm_rno_r
#endif
!#endif
  use m_Control_Parameters,  only : nblocksize_fftw, nblocksize_fftw_is_given &
       &                          , nblocksize_submat, nblocksize_subspace, sw_fft_xzy
  use m_IterationNumbers,     only : iteration, iteration_electronic, iteration_ionic &
       &                           , m_Iter_set_rmm_start, nk_in_the_process
  use m_ES_ExactExchange,    only : m_ES_Vexx_W
  use m_FiniteElectricField, only : m_FEF_add_grad_to_vnlph


! ================================= added by K. Tagami ============== 11.0
  use m_Control_Parameters,    only : ndim_spinor, ndim_chgpot
!  use m_FFT,                   only : m_FFT_Vlocal_W_noncl
!  use m_Electronic_Structure,      only : m_ES_sort_eigen_vals_noncl
! =================================================================== 11.0

  use m_Control_Parameters,  only : use_metagga, vtau_exists
  use m_ES_WF_by_SDorCG, only : m_ES_contrib_kindens_to_vnlph

  use m_Control_Parameters, only : sw_rsb
#ifdef MPI_FFTW
  use m_ES_WF_by_SDorCG,    only : gen_fft_to_WF_map, mmx, mmy, mmz, iiadd, wfdist, maxwfdist
  use m_Electronic_Structure, only : m_ES_fftbox_map
  use m_ES_WF_by_SDorCG,    only : m_ES_con_kindens_to_vnlph_mpfw
#endif
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

#ifdef MPI_FFTW
  include 'fftw3-mpi.f03'
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

    integer             :: ispin, iksnl, ik
    real(kind=DP), allocatable, dimension(:,:) ::  afft_l, cfft_l
    real(kind=DP), allocatable, dimension(:,:,:) ::  cfft_mpifftw_vlocal
    real(kind=DP), pointer, dimension(:) :: ekin_l
    integer             :: lsize

    integer :: n_all_kpoints, n_submat, ipri0
    logical :: submat_is_done
#ifdef MPI_FFTW
    integer(C_INTPTR_T)  :: local_n, local_n_offset, alloc_local, lx, ly, lz, mx,my,mz
#endif

    submat_is_done = .false.
    call m_ES_alloc_scss_etc_3D()
#ifdef FFT_3D_DIVISION
    lsize = fft_X_x_nel*fft_X_y_nel*fft_X_z_nel
    allocate(afft_l(lsize*2,1), stat=ierr)
    if ( use_metagga .and. vtau_exists ) allocate(cfft_l(lsize*2,1), stat=ierr)
#else
#ifdef MPI_FFTW
    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
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
      allocate(afft_mpifftw_vlocal(lx,lz,local_n));afft_mpifftw_vlocal=0.d0
      if ( use_metagga .and. vtau_exists ) then
         allocate(cfft_mpifftw_vlocal(lx,lz,local_n));cfft_mpifftw_vlocal=0.d0
      endif
    else
      allocate(afft_l(lsize*kimg,1), stat=ierr)
      if ( use_metagga .and. vtau_exists ) then
         allocate(cfft_l(lsize*kimg,1), stat=ierr)
      endif
    endif
#else
    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
    allocate(afft_l(lsize*kimg,1), stat=ierr)
    if ( use_metagga .and. vtau_exists ) allocate(cfft_l(lsize*kimg,1), stat=ierr)
#endif

#endif
     if(ierr /= 0) then
        write(nfout,*)' m_ESsubmat_Renew_WF : Not allocated afft_l array'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 101, ierr)
     endif
    ekin_l => sc_l

    n_all_kpoints = 0
    n_submat = 0
!!$    if(m_CtrlP_ntcnvg_clear()) non_diagonal_part_is_small = .false.
    if(iprisubmat >= 2) write(nfout,'(" !! <<m_ESsubmat_Renew_WF>>")')
    call get_ipri0(iprisubmat,ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko_3D(nfout,mode=SCF)


!    do ispin = 1, nspin, (af+1)
    do ispin = ista_spin, iend_spin, (af+1)
#ifdef MPI_FFTW
       if(sw_mpi_fftw==ON) then
         call m_ES_Vlocal_in_Rspace_mpifftw3d(ispin,lx,local_n,lz,afft_mpifftw_vlocal)
         if ( use_metagga .and. vtau_exists ) then
            call m_ES_Vlocal_in_Rspace_mpifftw3d( ispin, lx, local_n, lz, &
                 &                            cfft_mpifftw_vlocal, vtau_l )
         endif
       else
         call m_ES_Vlocal_in_Rspace_3D(ispin,afft_l,lsize,1,OFF)      ! (ptfft1) vlhxc_l->afft
         if ( use_metagga .and. vtau_exists ) then
            call m_ES_Vlocal_in_Rspace_3D( ispin, cfft_l, lsize, 1, OFF, vtau_l ) ! r space
         endif
       endif
#else
                                                  __TIMER_START_w_BARRIER(MPI_CommGroup,1)
       call m_ES_Vlocal_in_Rspace_3D(ispin,afft_l,lsize,1,OFF)   ! (ptfft1) vlhxc_l->afft
                                                  __TIMER_STOP(1)
       if ( use_metagga .and. vtau_exists ) then
          call m_ES_Vlocal_in_Rspace_3D( ispin, cfft_l, lsize, 1, OFF, vtau_l ) ! r space
       endif
#endif

       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) /= myrank_k) cycle          ! MPI
          n_all_kpoints = n_all_kpoints + 1
          if(non_diagonal_part_is_small(ik)) then
             if(iprisubmat >= 2) write(nfout,'(" !!! submat is not done for ik=",i4)') ik
             cycle
          else
             if(iprisubmat >= 2) write(nfout,'(" !!! submat is done for ik=",i4)') ik
             n_submat = n_submat + 1
             submat_is_done = .true.
          end if
          iksnl = (ik-1)/nspin + 1
                                                  __TIMER_START_w_BARRIER(mpi_k_world(myrank_k),2)
          call m_ES_Vnonlocal_W_3D(ik,iksnl,ispin,switch_of_eko_part=OFF)  ! (nonloc) ->(vnlph_l)
                                                  __TIMER_STOP_w_BARRIER(mpi_k_world(myrank_k),2)
! === Support Hybrid on 3D_Parallel by tkato 2013/02/10 ========================
          if(sw_hybrid_functional==ON) call m_ES_Vexx_W(ik) ! (exx) ->(vnlph_l)
! ==============================================================================

          if ( use_metagga .and. vtau_exists ) then
#ifdef MPI_FFTW
             if (sw_mpi_fftw==ON) then
                call m_ES_con_kindens_to_vnlph_mpfw( ispin, ik, lsize, &
                     &                    lx, ly, lz, local_n, cfft_mpifftw_vlocal )
             else
                call m_ES_contrib_kindens_to_vnlph( ispin, ik, lsize, cfft_l ) ! -> vnlph_l
             endif
#else
             call m_ES_contrib_kindens_to_vnlph( ispin, ik, lsize, cfft_l ) ! -> vnlph_l
#endif
          endif

          call m_pwBS_kinetic_energies(ik,vkxyz,ekin_l) ! (diakin) ->ekin
          call evolve_WFs_in_subspace_3D&     !-(m_ES_WF_by_submat)
                  &(ik,iksnl,ispin,meg,damp,ekin_l,afft_l(1,1),lsize)
          if(ik==1.and.iprisubmat>= 2) &
               & call m_ES_wd_zaj_small_portion_3D(nfout,1," -- after subspace-rotation --",30)
                                                  __TIMER_START_w_BARRIER(mpi_k_world(myrank_k),10)
!          call m_ES_betar_dot_WFs_4_each_k_3D(nfout,ik)   ! -> fsr_l,fsi_l
                                                  __TIMER_STOP_w_BARRIER(mpi_k_world(myrank_k),10)
       enddo      ! k-point loop
    enddo      ! spin loop

    if ( allocated(cfft_l) ) deallocate( cfft_l )
    if ( allocated(cfft_mpifftw_vlocal) ) deallocate( cfft_mpifftw_vlocal )

!!$    if(npes > 1) then
!!$       call mpi_allreduce(n_submat,n_submat_mpi,1,mpi_integer,mpi_sum,mpi_k_world(myrank_k),ierr)
!!$       n_submat = n_submat_mpi
!!$    end if

!!  ( in case of af=1 )
    if(af /= 0) then
       call cp_eigen_values_for_af       !-(contained here)
       call expand_neordr_and_nrvf_ordr  !-(contained here)
    end if

    if(ipri0 >= 2) call m_ES_wd_eko_3D(nfout,mode=SCF)
#ifdef MPI_FFTW
    if(sw_mpi_fftw==ON) then
      deallocate(afft_mpifftw_vlocal)
    else
      deallocate(afft_l)
    endif
#else
    deallocate(afft_l)
#endif
    call m_ES_dealloc_scss_etc()
    if(iprisubmat >= 2) then
       if(submat_is_done) write(nfout,'(" !solver -- subspace_rotation is done for ",i4," of ",i4," kpoints")') &
            & n_submat, n_all_kpoints
       if(.not.submat_is_done) write(nfout,'(" !solver --  no subspace_rotation")')
    end if

    if(nrank_k > 1) then
       call mpi_allreduce(MPI_IN_PLACE, n_submat,1,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
       call mpi_allreduce(MPI_IN_PLACE, n_all_kpoints,1,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
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




!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine evolve_WFs_in_subspace_3D&
       &(ik,iksnl,ispin,meg,damp,ekin_l,afft_l,lsize,dryrun)
!
! Revised by FUJITSU , March 2010 : MULTI DIMENSION
!
    use m_Parallelization,     only : fft_wf_scnt, fft_wf_rcnt &
   &                                , fft_wf_recv &
   &                                , fft_wf_index, fft_wf_dist &
   &                                , fft_wf_maxrecv, fft_wf_maxsend

#ifdef MPI_FFTW
    use m_Parallelization,     only : fft_wf_scnt_mfftw, fft_wf_rcnt_mfftw &
   &                                , fft_wf_recv_mfftw &
   &                                , fft_wf_index_mfftw, fft_wf_dist_mfftw &
   &                                , fft_wf_maxrecv_mfftw, fft_wf_maxsend_mfftw
#endif

    integer, intent(in) :: ik,iksnl,ispin,meg,lsize
    real(kind=DP), intent(in)  :: damp,ekin_l(np_g1k(ik))
#ifdef FFT_3D_DIVISION
    real(kind=DP), intent(inout) :: afft_l(lsize*2)
#else
    real(kind=DP), intent(inout) :: afft_l(lsize*kimg)
#endif
    logical, intent(in), optional :: dryrun
    logical :: dry
! (allocatable variables)
! === Debug by Intel "-check all" option! by T.Kato 2011/01/24 =================
!   real(kind=DP), allocatable,dimension(:,:) ::   bfft
! ==============================================================================
    real(kind=DP), allocatable,dimension(:) ::     eig
    real(kind=DP), allocatable,dimension(:) ::     eko_d
    real(kind=DP), allocatable,dimension(:,:) ::   w1hw2
!******************** modified by RIST_16
    real(kind=DP), allocatable,dimension(:,:) ::   w1sw2
!******************** modified by RIST_16
!   real(kind=DP), allocatable,dimension(:,:,:) :: zat_wk
    real(kind=DP), allocatable,dimension(:,:) :: wk_bfft_l
    real(kind=DP), allocatable,dimension(:,:) :: wk_bfft_l_in
    real(kind=DP), allocatable,dimension(:,:) :: bfft_l

!fj --------------------
!fj real(kind=DP), allocatable :: amat_l_3D(:,:), zmat_l_3D(:,:)
!fj real(kind=DP), allocatable :: wk_amat(:,:)
!fj --------------------
    real(kind=DP), allocatable :: amat_l(:,:)
    real(kind=DP), allocatable :: amat1_l(:,:)
    real(kind=DP), allocatable :: zmat(:,:), zmat_l(:,:), zmat_l_2D(:,:)
    real(kind=DP), allocatable :: zmat1(:,:), zmat1_l(:,:)
    real(kind=DP), allocatable :: zat_l(:,:,:,:)
    real(kind=DP), allocatable :: wk_mpi(:,:), wk_gather(:,:,:), wk_ball(:,:,:), wk_meg(:,:,:)
!fj integer, save :: lwork,lrwork,liwork
    integer, save :: lda,occ,lrwork
!#ifdef _USE_SCALAPACK_
    integer, save :: lwork1,lrwork1,liwork1
    integer, save :: lwork2,lrwork2,liwork2
    integer, save :: nb_mgs
    logical, save :: nsame
    integer, save :: nsclrow, nsclcol
    integer, save :: icolor
!#endif
#ifdef SINGLE_CONTEXT
    integer, save :: ictxt, myrow, mycol, ictxtt
#else
    integer, save :: myrow, mycol, tmpctxt
    integer, allocatable, save :: ictxt(:) ! d(0:nrank_k-1)
#endif
    integer, dimension(9), save :: desca,descz,descb
#ifdef SINGLE_CONTEXT
    integer, allocatable, save :: usermap(:,:)
#else
    integer, allocatable, save :: usermap(:,:,:) ! d(nprow,npcol,0:nrank_k-1)
#endif
    integer, allocatable, save :: icol(:)
#ifdef EIGEN_TRANS
    integer, allocatable, save :: irank_c(:)
    integer, allocatable, save :: irank_r(:)
#endif
!!$    integer, save :: npart,npart2
!!$    integer, allocatable, save :: isp(:),iep(:)
!!$    integer, allocatable, save :: isp2(:),iep2(:)
!    logical, save :: first = .true.
    integer :: is,ie,m,jp,ig1
    integer :: maxe, maxeg

!fj integer, allocatable, save, dimension(:) :: G_nel_e,G_nis_e,G_nie_e
!fj integer, dimension(meg) :: G_map_e

    integer       :: ib1,ib2,ib1to,i1,ii, i, kimg_t, iadd, j, k, l, ib1ad, ib2ad,ibb2,ibb1
    integer       :: ibsize, iesize, jj
    real(kind=DP) :: denom, eko1, ekod
    real(kind=DP) :: dr1,dr2,di1,di2, dd
    real(kind=DP) :: sum_sq_diagonal, sum_sq_non_diagonal &
         & , sum_abs_diagonal, sum_abs_non_diagonal
    integer :: id_sname = -1, id_sname2 = -1, ipri0
#ifdef SUBMAT_DGEMM
    real(kind=DP),allocatable,dimension(:,:) :: w1hw2r,w1hw2i
!fj real(kind=DP) :: alpha, beta
    integer :: imax
!******************** modified by RIST_16
    real(kind=DP),allocatable,dimension(:,:) :: w1sw2r,w1sw2i
    real(kind=DP) ::fr,fi
!******************** modified by RIST_16

    real(kind=DP),allocatable,dimension(:,:) :: fsr_in, fsi_in
#endif

    integer :: isrsize, num
!**** TYamasaki 2019/06/10
    integer :: nb,kb,jb,ib
!**** <===
    integer :: ilmta

    integer, dimension(0:nrank_g-1) ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)::sta_r, sta_s
     real(kind=DP), allocatable, dimension(:,:) :: sendbuf, recvbuf
     integer :: icnt_send, icnt_recv, lrank
     integer :: itag = 10
! === DEBUG by tkato 2011/10/24 ================================================
     ! Revised by T. Yamasaki 2013/05/14
     !  real(kind=DP) :: wk1(maxval(np_g1k),neg,kimg)
     real(kind=DP),allocatable,dimension(:,:,:) :: wk1, wk2 ! wk1(maxval(np_g1k),neg,kimg)
! ==============================================================================
! === DEBUG by tkato 2012/06/05 ================================================
#ifdef USE_ALLTOALLV
     integer, allocatable, dimension(:) :: sdsp, rdsp
#endif
! ==============================================================================
#ifdef EIGEN_TRANS
#else
! === DEBUG by tkato 2013/09/18 ================================================
     integer, save, dimension(9) :: desc0, desc1
#ifdef SINGLE_CONTEXT
     integer, save               :: ictxt0 ! , ictxt1
#else
     integer, allocatable, save  :: ictxt0(:) ! , ictxt1
#endif
! ==============================================================================
#endif
#ifdef MPI_FFTW
    integer :: iwf
    integer :: mm,nsize
    integer(C_INTPTR_T)  :: local_n, local_n_offset, alloc_local, lx, ly, lz, mx,my,mz
#endif
#ifdef _USE_SCALAPACK_
    integer :: sys2blacs_handle
#endif
                                                  __TIMER_START(8)
                                                  __TIMER_SUB_START(901)

#ifdef __FAPP__
    call fapp_start('submat_pre',1,1)
#endif
    call tstatc0_begin('evolve_WFs_in_subspace_3D ', id_sname,1)

    dry = .false.
    if(present(dryrun)) dry = dryrun

    denom = 1.d0/product(fft_box_size_WF(1:3,1))

    if(k_symmetry(ik) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = kimg
    end if

    allocate(eko_d(neg))
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

!! Not use !!
!   do ib1 = 1, neg
!      if(map_e(ib1) == myrank_e) eko_d(ib1) = eko_l(map_z(ib1),ik)  ! MPI
!   end do
!   call mpi_allreduce(eko_d,eko_d_mpi,neg,mpi_double_precision,mpi_sum &
!        & ,mpi_k_world(myrank_k),ierr)       ! MPI
!   eko_d = eko_d_mpi                         ! MPI
    eko_d = 0.0d0
    do ib1 = 1, np_e
       eko_d(neordr(neg_g(ib1),ik)) = eko_l(ib1,ik)  ! MPI
    end do
                                                  __TIMER_COMM_START_w_BARRIER(mpi_kg_world,914)
    call mpi_allreduce(MPI_IN_PLACE,eko_d,neg,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
                                                  __TIMER_COMM_STOP(914)
    eko1 = sum(eko_d(1:neg))
! == revised by T.Yamasaki 2019/06/10 ==============================================================================
!!$    allocate(wk1(maxval(np_g1k),neg,kimg))
    if(nrank_e == 1) then
       zaj_ball(:,:,ik,:) = zaj_l(:,:,ik,:)
    else
!!$       write(nfout,'(" *** mpi_allreduce size = ",i12," *** ")') maxval(np_g1k)*neg*kimg*8
!!$       write(nfout,'(" ***       reduced size = ",i12," *** ")') maxval(np_g1k)*mp_e*kimg
!!$       call flush(nfout)
       allocate(wk1(maxval(np_g1k),mp_e,kimg),stat=ierr)
       do nb = 0, nrank_e-1
          if(nb == myrank_e) then
             wk1(:,1:np_e,1:kimg) = zaj_l(:,1:np_e,ik,1:kimg)
             if(np_e+1<mp_e) wk1(:,np_e+1:mp_e,:) = 0.d0
          end if
          call mpi_bcast(wk1,maxval(np_g1k)*mp_e*kimg,mpi_double_precision,nb,mpi_kg_world,ierr)
          do kb = 1, kimg
             do jb = nis_e(nb),nie_e(nb)
                do ib = 1, np_g1k(ik)
                   zaj_ball(ib,jb,ik,kb) = wk1(ib,jb-nis_e(nb)+1,kb)
                end do
             end do
          end do
       end do
       deallocate(wk1)
    end if

    allocate(fsr_in(neg,np_fs));fsr_in=0.d0
    if(k_symmetry(ik) /= GAMMA) then
      allocate(fsi_in(neg,np_fs))
      fsi_in = 0.d0
    endif
    if(nrank_e == 1)then
      do ib=1,neg
        do ilmta=1,np_fs
          fsr_in(ib,ilmta) = fsr_l(ib,ilmta,ik)
        enddo
      enddo
      if(k_symmetry(ik) /= GAMMA) then
        do ib=1,neg
          do ilmta=1,np_fs
            fsi_in(ib,ilmta) = fsi_l(ib,ilmta,ik)
          enddo
        enddo
      endif
    else
      do ib = ista_e, iend_e
        do ilmta=1,np_fs
          fsr_in(ib,ilmta) = fsr_l(ib-ista_e+1,ilmta,ik)
          if(k_symmetry(ik) /= GAMMA) fsi_in(ib,ilmta) = fsi_l(ib-ista_e+1,ilmta,ik)
        enddo
      enddo
      call mpi_allreduce(mpi_in_place,fsr_in,neg*np_fs,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
      if(k_symmetry(ik) /= GAMMA) &
      call mpi_allreduce(mpi_in_place,fsi_in,neg*np_fs,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
    endif

!!$! === DEBUG by tkato 2011/10/24 ================================================
!!$    wk1 = 0.0d0
!!$    write(nfout,'(" *** mpi_allreduce size = ",i12," *** ")') maxval(np_g1k)*neg*kimg*8
!!$    do i = ista_e, iend_e
!!$       wk1(:,i,:) = zaj_l(:,i-ista_e+1,ik,:)
!!$    enddo
!!$    call mpi_allreduce(MPI_IN_PLACE,wk1,maxval(np_g1k)*neg*kimg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
!!$    do i = 1, neg
!!$       zaj_ball(:,i,ik,:) = wk1(:,i,:)
!!$    enddo
!!$! ==============================================================================
!!$    deallocate(wk1)
! ==============================================================================


    if(sw_gep == ON)then
       allocate(zat_l(maxval(np_g1k),np_e,ista_k:iend_k,kimg)); zat_l = zaj_l
    endif

    ibsize = 1
    if (nblocksize_fftw_is_given) then
       ibsize = nblocksize_fftw
       if (ibsize < 1) ibsize = 1
    endif
    isrsize = min(lsize,mp_g1k(ik))
! === DEBUG by tkato 2011/07/12 ================================================
!   bfft = 0.0d0
! ==============================================================================

    maxe = mp_g1k(ik)
#ifdef FFT_3D_DIVISION
    allocate(wk_bfft_l(lsize*2,ibsize) ,stat=ierr)
    allocate(bfft_l(lsize*2,ibsize) ,stat=ierr)
#else
!    allocate(wk_bfft_l(lsize*kimg,ibsize) ,stat=ierr)
!    allocate(wk_bfft_l_in(lsize*kimg,ibsize) ,stat=ierr);wk_bfft_l_in = 0.d0
!    allocate(bfft_l(lsize*kimg,ibsize) ,stat=ierr)
    allocate(wk_bfft_l(lsize*kimg,ibsize))
!    allocate(wk_bfft_l_in(lsize*kimg,ibsize));wk_bfft_l_in = 0.d0
    allocate(bfft_l(lsize*kimg,ibsize))
#endif
!   allocate(bfft(nfft,ibsize), stat=ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_ES_WF_in_Rspace_3D :  Not allocate '
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 109, ierr)
     endif

! (zaj_l <- H |phi> )
    if(.not.dry) then
                                                  __TIMER_DO_START(960)
#ifdef MPI_FFTW
    if(sw_mpi_fftw==ON) then
      call gen_fft_to_WF_map(ik)
      call m_ES_fftbox_map(ik)
    endif
#endif
    do ib1 = 1, np_e, ibsize
       ib2 = min(ib1+ibsize-1,np_e)
       iesize = ib2 - ib1 + 1

#ifdef __TIMER_COMM__
       call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,wk_bfft_l, 8)
#else
#ifdef MPI_FFTW
       if(sw_mpi_fftw==ON) then
         call m_ES_WF_in_Rspace_mpifftw(ista_k,iend_k,ik,ib1,zaj_l)
       else
         call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,wk_bfft_l)
       endif
#else
       call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,wk_bfft_l)
!       call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,wk_bfft_l_in,wk_bfft_l)
#endif
#endif
#ifdef FFT_3D_DIVISION
       call m_FFT_Vlocal_W_3DIV_3D(afft_l,wk_bfft_l,lsize,iesize,nel_fft_x(myrank_g))
#else
#ifdef MPI_FFTW
       if(sw_mpi_fftw==ON) then
         lx = fft_box_size_WF(1,0)
         ly = fft_box_size_WF(2,0)
         lz = fft_box_size_WF(3,0)
         alloc_local = fftw_mpi_local_size_3d(lz,ly,lx,mpi_ke_world,local_n,local_n_offset)
         nsize = local_n*lx*ly
         call m_FFT_Vlocal_W_mpifftw3d(afft_mpifftw_vlocal,lx,local_n,lz)
       else
         call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_z(myrank_g))
       endif
#else
       if (sw_fft_xzy > 0) then
          call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,iesize,nel_fft_y(myrank_g))
       else
          call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,iesize,nel_fft_z(myrank_g))
       end if
#endif
#endif
                                                  __TIMER_COMM_START(962)
#ifdef FFT_3D_DIVISION
       call m_FFT_Direct_3DIV_3D (nfout, wk_bfft_l, lsize, iesize)
#else
#ifdef MPI_FFTW
       if(sw_mpi_fftw==ON) then
         call m_FFT_Direct_MPI_FFTW(nfout)
       else
         call m_FFT_Direct_XYZ_3D (nfout, wk_bfft_l, lsize, iesize)
       endif
#else
       if (sw_fft_xzy > 0) then
          call m_FFT_Direct_3D (nfout, wk_bfft_l, lsize, iesize)
       else
          if(sw_serial_fft == ON) then
             call m_ES_WF_2D(ik,wk_bfft_l,ib2,ib1,ibsize,lsize,DIRECT)
          else
             call m_FFT_Direct_XYZ_3D (nfout, wk_bfft_l, lsize, iesize)
          endif
       endif
#endif
#endif
                                                  __TIMER_COMM_STOP(962)
#ifdef MPI_FFTW
       if(sw_mpi_fftw==ON) then
       if (fft_wf_maxsend_mfftw(ik) /= 0) then
          allocate(sendbuf(fft_wf_maxsend_mfftw(ik)*kimg*ibsize,0:nrank_g-1), stat=ierr)
!          sendbuf = 0.0d0
! === DEBUG by tkato 2012/06/05 ================================================
       else
          allocate(sendbuf(1,1), stat=ierr)
! ==============================================================================
       endif
       if (fft_wf_maxrecv_mfftw(ik) /= 0) then
          allocate(recvbuf(fft_wf_maxrecv_mfftw(ik)*kimg*ibsize,0:nrank_g-1), stat=ierr)
!          recvbuf = 0.0d0
! === DEBUG by tkato 2012/06/05 ================================================
       else
          allocate(recvbuf(1,1), stat=ierr)
! ==============================================================================
       endif
        if (ierr /= 0) then
           write(nfout,*)' evolve_WFs_in_subspace_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 110, ierr)
        endif

        req_r = 0
        req_s = 0
        sta_r = 0
        sta_s = 0
       if (fft_wf_maxrecv_mfftw(ik) /= 0) then
          icnt_recv = 0
          lrank = mod(myrank_g,nrank_g)
          do i = 0, nrank_g - 1
             lrank = lrank + 1
             if (lrank > nrank_g -1) lrank = 0
             if (fft_wf_rcnt_mfftw(lrank,ik) /= 0) then
                call mpi_irecv(recvbuf(1,lrank), fft_wf_rcnt_mfftw(lrank,ik)*kimg*iesize, &
               &     mpi_double_precision, lrank, itag, mpi_ke_world, req_r(icnt_recv), ierr)
                 if (ierr /= 0) then
                    write(nfout,*)' evolve_WFs_in_subspace_3D :  mpi_irecv error'
                    call flush(nfout)
                    call mpi_abort(mpi_comm_world, 111, ierr)
                 endif
                icnt_recv = icnt_recv + 1
             endif
          enddo
       endif
                                                  __TIMER_DO_START(916)
       alloc_local = fftw_mpi_local_size_3d(lz,ly,lx,mpi_ke_world,local_n,local_n_offset)
       if (fft_wf_maxsend_mfftw(ik) /= 0) then
          if (kimg == 1) then
!OCL NORECURRENCE
             !do k = 1, nel_fft_x(myrank_g)
             do k = 1, local_n*lx*ly
                if(fft_wf_index_mfftw(k,ik) == 0) cycle
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
                sendbuf(fft_wf_index_mfftw(k,ik),fft_wf_dist_mfftw(k,ik)) = bfft_mpifftw_kimg1(mx,my,mz)
             end do
          else
!OCL NORECURRENCE
!             do k = 1, local_n*lx*ly
!                if(fft_wf_index_mfftw(k,ik) == 0) cycle
!                iadd = 2*(fft_wf_index_mfftw(k,ik))
!                i1 = k
!                mz = (i1-1)/(lx*ly)+1
!                mm = mod(i1,(lx*ly))
!                if (mm==0) mm=lx*ly
!                my = (mm-1)/lx+1
!                mx = mod(mm,lx)
!                if (mx==0) mx = lx
!                sendbuf(iadd-1,fft_wf_dist_mfftw(k,ik)) = real (bfft_mpifftw(mx,my,mz))
!                sendbuf(iadd,  fft_wf_dist_mfftw(k,ik)) = aimag(bfft_mpifftw(mx,my,mz))
!             end do
             do k=1, maxwfdist
               iadd = iiadd(k)
               iwf  = wfdist(k)
               mx   = mmx(k)
               my   = mmy(k)
               mz   = mmz(k)
               sendbuf(iadd-1,iwf) = real (bfft_mpifftw(mx,my,mz))
               sendbuf(iadd,  iwf) = aimag(bfft_mpifftw(mx,my,mz))
             enddo
          end if
       end if
                                                  __TIMER_DO_STOP(916)

       if (fft_wf_maxsend_mfftw(ik) /= 0) then
          icnt_send = 0
          lrank = mod((myrank_g+1),nrank_g)
          do i = 0, nrank_g - 1
             lrank = lrank + 1
             if (lrank > (nrank_g - 1)) lrank = 0
             if (fft_wf_scnt_mfftw(lrank,ik) /= 0) then
                call mpi_isend(sendbuf(1,lrank), fft_wf_scnt_mfftw(lrank,ik)*kimg*iesize, &
               &     mpi_double_precision, lrank, itag, mpi_ke_world, req_s(icnt_send), ierr)
                 if (ierr /= 0) then
                    write(nfout,*)' evolve_WFs_in_subspace_3D :  mpi_isend error'
                    call flush(nfout)
                    call mpi_abort(mpi_comm_world, 112, ierr)
                 endif
                icnt_send = icnt_send + 1
             endif
          enddo
       endif

       if (fft_wf_maxrecv_mfftw(ik) /= 0) then
          call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
           if (ierr /= 0) then
              write(nfout,*)' evolve_WFs_in_subspace_3D :  mpi_waitall error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 113, ierr)
           endif
       endif

       if (fft_wf_maxsend_mfftw(ik) /= 0) then
          call mpi_waitall(icnt_send, req_s, sta_s, ierr)
           if (ierr /= 0) then
              write(nfout,*)' evolve_WFs_in_subspace_3D :  mpi_waitall error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 114, ierr)
           endif
       endif
                                                  __TIMER_COMM_STOP(915)
!       bfft_l = 0.0d0
                                                  __TIMER_DO_START(916)
       if (kimg == 1) then
!x!OCL NORECURRENCE
          do i = 0, nrank_g - 1
             if (fft_wf_rcnt_mfftw(i,ik) /= 0) then
                do k = 1, fft_wf_rcnt_mfftw(i,ik)
                   do j = 1, iesize
                      bfft_l(fft_wf_recv_mfftw(k,ik,i),j) = recvbuf(iesize*(k-1)+j,i)
                   enddo
                end do
             end if
          end do
       else
!x!OCL NORECURRENCE
          do i = 0, nrank_g - 1
             if (fft_wf_rcnt_mfftw(i,ik) /= 0) then
                do k = 1, fft_wf_rcnt_mfftw(i,ik)
                   do j = 1, iesize
                      iadd = iesize*2*(k-1)+j*2
                      bfft_l(fft_wf_recv_mfftw(k,ik,i)*2-1,j) = recvbuf(iadd-1,i)
                      bfft_l(fft_wf_recv_mfftw(k,ik,i)*2  ,j) = recvbuf(iadd,  i)
                   enddo
                end do
             end if
          end do
       end if
                                                  __TIMER_DO_STOP(916)

       if (allocated(sendbuf)) deallocate(sendbuf)
       if (allocated(recvbuf)) deallocate(recvbuf)

                                                  __TIMER_DO_START(918)
       if(kimg == 1) then
!x!OCL NORECURRENCE
          do ii=ista_g1k(ik), iend_g1k(ik)
             i1  = igf(nbase(ii,ik))
             iadd = ii-ista_g1k(ik)+1
             do jj = ib1, ib2
!fj --------------------
                if(neg_g(jj)>meg) cycle
!fj --------------------
                dr1 = zaj_l(iadd,jj,ik,1)
                dr2 = bfft_l(iadd,jj-ib1+1)*denom
                zaj_l(iadd,jj,ik,1)= ekin_l(iadd)*dr1+vnlph_l(iadd,jj,1)+dr2
                if(sw_keep_hloc_phi==ON) hlocphi_l(iadd,jj,ik,1)= ekin_l(iadd)*dr1+dr2
             enddo
          enddo
       else
!x!OCL NORECURRENCE
           do ii=ista_g1k(ik), iend_g1k(ik)
              i1  = igf(nbase(ii,ik))
              iadd = ii-ista_g1k(ik)+1
              do jj = ib1, ib2
!!fj --------------------
                 if(neg_g(jj)>meg) cycle
!!fj --------------------
                 dr1  = zaj_l(iadd,jj,ik,1)
                 di1  = zaj_l(iadd,jj,ik,2)
                 zaj_l(iadd,jj,ik,1)= ekin_l(iadd)*dr1+bfft_l(2*iadd-1,jj-ib1+1)*denom+vnlph_l(iadd,jj,1)
                 zaj_l(iadd,jj,ik,2)= ekin_l(iadd)*di1+bfft_l(2*iadd  ,jj-ib1+1)*denom+vnlph_l(iadd,jj,2)
                 if(sw_keep_hloc_phi==ON) then
                   hlocphi_l(iadd,jj,ik,1)= ekin_l(iadd)*dr1+bfft_l(2*iadd-1,jj-ib1+1)*denom
                   hlocphi_l(iadd,jj,ik,2)= ekin_l(iadd)*di1+bfft_l(2*iadd  ,jj-ib1+1)*denom
                 endif
              enddo
           enddo
       end if
       else
#endif
       if (fft_wf_maxsend(ik) /= 0) then
          allocate(sendbuf(fft_wf_maxsend(ik)*kimg*ibsize,0:nrank_g-1), stat=ierr)
!          sendbuf = 0.0d0
! === DEBUG by tkato 2012/06/05 ================================================
       else
          allocate(sendbuf(1,1), stat=ierr)
! ==============================================================================
       endif
       if (fft_wf_maxrecv(ik) /= 0) then
          allocate(recvbuf(fft_wf_maxrecv(ik)*kimg*ibsize,0:nrank_g-1), stat=ierr)
!          recvbuf = 0.0d0
! === DEBUG by tkato 2012/06/05 ================================================
       else
          allocate(recvbuf(1,1), stat=ierr)
! ==============================================================================
       endif
        if (ierr /= 0) then
           write(nfout,*)' evolve_WFs_in_subspace_3D :  Not allocate '
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 110, ierr)
        endif

    req_r = 0
    req_s = 0
    sta_r = 0
    sta_s = 0
       if (fft_wf_maxrecv(ik) /= 0) then
          icnt_recv = 0
          lrank = mod(myrank_g,nrank_g)
#ifndef USE_ALLTOALLV
                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,915)
#endif
          do i = 0, nrank_g - 1
             lrank = lrank + 1
             if (lrank > nrank_g -1) lrank = 0
             if (fft_wf_rcnt(lrank,ik) /= 0) then
#ifndef USE_ALLTOALLV
                call mpi_irecv(recvbuf(1,lrank), fft_wf_rcnt(lrank,ik)*kimg*iesize, &
               &     mpi_double_precision, lrank, itag, mpi_ke_world, req_r(icnt_recv), ierr)
                 if (ierr /= 0) then
                    write(nfout,*)' evolve_WFs_in_subspace_3D :  mpi_irecv error'
                    call flush(nfout)
                    call mpi_abort(mpi_comm_world, 111, ierr)
                 endif
#endif
                icnt_recv = icnt_recv + 1
             endif
          enddo
       endif
                                                  __TIMER_DO_START(916)
       if (fft_wf_maxsend(ik) /= 0) then
          if (kimg == 1) then
#ifdef FFT_3D_DIVISION
          integer :: i1, lx, ly, lz, mx, my, mz, mm, kx1p, kx2p, kx3p, jadd
             lx = fft_box_size_WF(1,0)
             ly = fft_box_size_WF(2,0)
             lz = fft_box_size_WF(3,0)
             kx1p = fft_X_x_nel
             kx2p = fft_X_y_nel
             kx3p = fft_X_z_nel
!OCL NORECURRENCE
             do k = 1, nel_fft_x(myrank_g)
                if(fft_wf_index(k,ik) == 0) cycle
                i1 = mp_fft_x(k)
                mz = (i1-1)/(lx*ly)+1
                mm = mod(i1,(lx*ly))
                if (mm==0) mm=lx*ly
                my = (mm-1)/lx+1
                mx = mod(mm,lx)
                if (mx==0) mx = lx
                jadd = mx-xyz_fft_x(1,1)+1+kx1p*(my-xyz_fft_x(1,2))+kx1p*kx2p*(mz-xyz_fft_x(1,3))
                do i = 1, iesize
                   sendbuf(iesize*(fft_wf_index(k,ik)-1)+i,fft_wf_dist(k,ik)) = wk_bfft_l(jadd*2-1,i)
                enddo
             end do

#else
!OCL NORECURRENCE
             do k = 1, nel_fft_x(myrank_g)
                if(fft_wf_index(k,ik) == 0) cycle
                do i = 1, iesize
                   sendbuf(iesize*(fft_wf_index(k,ik)-1)+i,fft_wf_dist(k,ik)) = wk_bfft_l(k,i)
                enddo
             end do
#endif
          else
#ifdef FFT_3D_DIVISION
             lx = fft_box_size_WF(1,0)
             ly = fft_box_size_WF(2,0)
             lz = fft_box_size_WF(3,0)
             kx1p = fft_X_x_nel
             kx2p = fft_X_y_nel
             kx3p = fft_X_z_nel
!OCL NORECURRENCE
             do k = 1, nel_fft_x(myrank_g)
                if(fft_wf_index(k,ik) == 0) cycle
                i1 = mp_fft_x(k)
                mz = (i1-1)/(lx*ly)+1
                mm = mod(i1,(lx*ly))
                if (mm==0) mm=lx*ly
                my = (mm-1)/lx+1
                mx = mod(mm,lx)
                if (mx==0) mx = lx
                jadd = mx-xyz_fft_x(1,1)+1+kx1p*(my-xyz_fft_x(1,2))+kx1p*kx2p*(mz-xyz_fft_x(1,3))
                do i = 1, iesize
                   iadd = iesize*2*(fft_wf_index(k,ik)-1)+i*2
                   sendbuf(iadd-1,fft_wf_dist(k,ik)) = wk_bfft_l(jadd*2-1,i)
                   sendbuf(iadd,  fft_wf_dist(k,ik)) = wk_bfft_l(jadd*2  ,i)
                enddo
             end do
#else
!OCL NORECURRENCE
             do k = 1, nel_fft_x(myrank_g)
                if(fft_wf_index(k,ik) == 0) cycle
                do i = 1, iesize
                   iadd = iesize*2*(fft_wf_index(k,ik)-1)+i*2
                   sendbuf(iadd-1,fft_wf_dist(k,ik)) = wk_bfft_l(k*2-1,i)
                   sendbuf(iadd,  fft_wf_dist(k,ik)) = wk_bfft_l(k*2  ,i)
                enddo
             end do
#endif
          end if
       end if
                                                  __TIMER_DO_STOP(916)

       if (fft_wf_maxsend(ik) /= 0) then
          icnt_send = 0
          lrank = mod((myrank_g+1),nrank_g)
          do i = 0, nrank_g - 1
             lrank = lrank + 1
             if (lrank > (nrank_g - 1)) lrank = 0
             if (fft_wf_scnt(lrank,ik) /= 0) then
#ifndef USE_ALLTOALLV
                call mpi_isend(sendbuf(1,lrank), fft_wf_scnt(lrank,ik)*kimg*iesize, &
               &     mpi_double_precision, lrank, itag, mpi_ke_world, req_s(icnt_send), ierr)
                 if (ierr /= 0) then
                    write(nfout,*)' evolve_WFs_in_subspace_3D :  mpi_isend error'
                    call flush(nfout)
                    call mpi_abort(mpi_comm_world, 112, ierr)
                 endif
#endif
                icnt_send = icnt_send + 1
             endif
          enddo
       endif

#ifndef USE_ALLTOALLV
       if (fft_wf_maxrecv(ik) /= 0) then
          call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
           if (ierr /= 0) then
              write(nfout,*)' evolve_WFs_in_subspace_3D :  mpi_waitall error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 113, ierr)
           endif
       endif

       if (fft_wf_maxsend(ik) /= 0) then
          call mpi_waitall(icnt_send, req_s, sta_s, ierr)
           if (ierr /= 0) then
              write(nfout,*)' evolve_WFs_in_subspace_3D :  mpi_waitall error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 114, ierr)
           endif
       endif
                                                  __TIMER_COMM_STOP(915)
#else
                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,994)
! === DEBUG by tkato 2012/06/05 ================================================
!      integer, allocatable, dimension(:) :: sdsp, rdsp
! ==============================================================================
       allocate(sdsp(0:nrank_g-1), stat=ierr)
       allocate(rdsp(0:nrank_g-1), stat=ierr)
       do i = 0, nrank_g - 1
          sdsp(i)=fft_wf_maxsend(ik)*kimg*ibsize*i
          rdsp(i)=fft_wf_maxrecv(ik)*kimg*ibsize*i
       enddo
       call MPI_ALLTOALLV(      sendbuf, fft_wf_scnt(:,ik)*kimg*iesize, sdsp, &
      &   mpi_double_precision, recvbuf, fft_wf_rcnt(:,ik)*kimg*iesize, rdsp, &
      &   mpi_double_precision, mpi_ke_world, ierr )
       if (ierr /= 0) then
          write(nfout,*)' evolve_WFs_in_subspace_3D :  mpi_alltoallv error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 115, ierr)
       endif
       deallocate(sdsp)
       deallocate(rdsp)
                                                  __TIMER_COMM_STOP(994)
#endif

!       bfft_l = 0.0d0
                                                  __TIMER_DO_START(916)
       if (kimg == 1) then
!x!OCL NORECURRENCE
          do i = 0, nrank_g - 1
             if (fft_wf_rcnt(i,ik) /= 0) then
                do k = 1, fft_wf_rcnt(i,ik)
                   do j = 1, iesize
                      bfft_l(fft_wf_recv(k,ik,i),j) = recvbuf(iesize*(k-1)+j,i)
                   enddo
                end do
             end if
          end do
       else
!x!OCL NORECURRENCE
          do i = 0, nrank_g - 1
             if (fft_wf_rcnt(i,ik) /= 0) then
                do k = 1, fft_wf_rcnt(i,ik)
                   do j = 1, iesize
                      iadd = iesize*2*(k-1)+j*2
                      bfft_l(fft_wf_recv(k,ik,i)*2-1,j) = recvbuf(iadd-1,i)
                      bfft_l(fft_wf_recv(k,ik,i)*2  ,j) = recvbuf(iadd,  i)
                   enddo
                end do
             end if
          end do
       end if
                                                  __TIMER_DO_STOP(916)

       if (allocated(sendbuf)) deallocate(sendbuf)
       if (allocated(recvbuf)) deallocate(recvbuf)

                                                  __TIMER_DO_START(918)
       if(kimg == 1) then
!x!OCL NORECURRENCE
          do ii=ista_g1k(ik), iend_g1k(ik)
             i1  = igf(nbase(ii,ik))
             iadd = ii-ista_g1k(ik)+1
             do jj = ib1, ib2
!fj --------------------
                if(neg_g(jj)>meg) cycle
!fj --------------------
                dr1 = zaj_l(iadd,jj,ik,1)
                dr2 = bfft_l(iadd,jj-ib1+1)*denom
                zaj_l(iadd,jj,ik,1)= ekin_l(iadd)*dr1+vnlph_l(iadd,jj,1)+dr2
                if(sw_keep_hloc_phi==ON) hlocphi_l(iadd,jj,ik,1)= ekin_l(iadd)*dr1+dr2
             enddo
          enddo
       else
!x!OCL NORECURRENCE
           do ii=ista_g1k(ik), iend_g1k(ik)
              i1  = igf(nbase(ii,ik))
              iadd = ii-ista_g1k(ik)+1
              do jj = ib1, ib2
!!fj --------------------
                 if(neg_g(jj)>meg) cycle
!!fj --------------------
                 dr1  = zaj_l(iadd,jj,ik,1)
                 di1  = zaj_l(iadd,jj,ik,2)
                 zaj_l(iadd,jj,ik,1)= ekin_l(iadd)*dr1+bfft_l(2*iadd-1,jj-ib1+1)*denom+vnlph_l(iadd,jj,1)
                 zaj_l(iadd,jj,ik,2)= ekin_l(iadd)*di1+bfft_l(2*iadd  ,jj-ib1+1)*denom+vnlph_l(iadd,jj,2)
                 if(sw_keep_hloc_phi==ON) then
                   hlocphi_l(iadd,jj,ik,1)= ekin_l(iadd)*dr1+bfft_l(2*iadd-1,jj-ib1+1)*denom
                   hlocphi_l(iadd,jj,ik,2)= ekin_l(iadd)*di1+bfft_l(2*iadd  ,jj-ib1+1)*denom
                 endif
              enddo
           enddo
       end if
#ifdef SAVE_FFT_TIMES
       if(sw_save_fft == ON) status_saved_phifftr(ib1:ib2,ik) = OLD
#endif
                                                  __TIMER_DO_STOP(918)
#ifdef MPI_FFTW
    endif
#endif
    end do
                                                  __TIMER_DO_STOP(960)
    end if

    deallocate(wk_bfft_l)
!    deallocate(wk_bfft_l_in)

    if(sw_keep_hloc_phi==ON) then
      if(nrank_e == 1) then
         zah_ball(:,:,:) = hlocphi_l(:,:,ik,:)
      else
         allocate(wk1(maxval(np_g1k),mp_e,kimg),stat=ierr)
         do nb = 0, nrank_e-1
            if(nb == myrank_e) then
               wk1(:,1:np_e,1:kimg) = hlocphi_l(:,1:np_e,ik,1:kimg)
               if(np_e+1<mp_e) wk1(:,np_e+1:mp_e,:) = 0.d0
            end if
            call mpi_bcast(wk1,maxval(np_g1k)*mp_e*kimg,mpi_double_precision,nb,mpi_kg_world,ierr)
            do kb = 1, kimg
               do jb = nis_e(nb),nie_e(nb)
                  do ib = 1, np_g1k(ik)
                     zah_ball(ib,jb,kb) = wk1(ib,jb-nis_e(nb)+1,kb)
                  end do
               end do
            end do
         end do
         deallocate(wk1)
      end if
    endif

    if(submat_uncalled) then
       call set_col_partition_3D(neg) ! -> npart,isp,iep, npart2, isp2, iep2
!      call set_block_range(meg,nrank_g,G_nel_e,G_nis_e,G_nie_e,G_map_e)
    end if
    maxe = maxval(nel_e(:))
    maxeg = maxval(nel_eg(:))

! (make matrix elements )
#ifdef _USE_SCALAPACK_
    if(iprisubmat>=2) write(nfout,'(" sw_scalapack = ",i3," <<evolve_WFs_in_subspace>>")') sw_scalapack
    if(sw_scalapack == ON) then
       if(submat_uncalled) then
          if(iprisubmat >= 2) write(nfout,'("first nprow,npcol=",2i5)') nprow,npcol
!fj --------------------
          if(nblocksize_mgs_is_given) then
             nb_mgs = nblocksize_mgs
          else
             nb_mgs = nblocksize_mgs_default
          end if
!fj    -----------------
          if (nprow==0 .or. npcol==0 .or. nprow*npcol > nrank_g*nrank_e) then

             nprow = int(sqrt(real(nrank_e*nrank_g)))

             if (divide_square > 0) then
                npcol = nprow
             else
                do
                   if ( nprow <= 1 ) exit
                   if ( mod(nrank_g*nrank_e,nprow) == 0 ) exit
                   nprow = nprow-1
                end do
                npcol = nrank_g*nrank_e/nprow
             endif
          endif
!xx       if(nprow*npcol /= nrank_g*nrank_e) then
!xx          method_scalapack=HOUSEHOLDER
!xx          write(nfout,'("ScaLAPACK METHOD ---> HOUSEHOLDER")')
!xx       end if
!fj    -----------------
          call set_nprow_npcol(nprow,npcol)
#ifdef SINGLE_CONTEXT
          if(allocated(usermap)) deallocate(usermap)
          allocate(usermap(nprow,npcol))
#else
          if(allocated(ictxt)) deallocate(ictxt)
          if(allocated(ictxt0)) deallocate(ictxt0)
          if(allocated(usermap)) deallocate(usermap)
          allocate(ictxt(0:nrank_k-1))
          allocate(ictxt0(0:nrank_k-1))
          allocate(usermap(nprow,npcol,0:nrank_k-1))
#endif
          if(allocated(icol)) deallocate(icol)
          allocate(icol(0:nrank_s*nrank_e*nrank_g*nrank_k-1))
#ifdef EIGEN_TRANS
          if(allocated(irank_c)) deallocate(irank_c)
          allocate(irank_c(0:nrank_s*nrank_e*nrank_g*nrank_k-1))
          if(allocated(irank_r)) deallocate(irank_r)
          allocate(irank_r(0:nrank_s*nrank_e*nrank_g*nrank_k-1))
#endif
          nsclrow = nprow
          nsclcol = npcol
#ifdef USE_EIGENLIB
          if(sw_eigen_exa==ON) then
            nsame = .false.
            block_size = 1

#ifdef EIGEN_TRANS
            call eigen_setup(meg,block_size,nprow,npcol,lda,occ,icolor,ictxt,usermap,icol,irank_c,irank_r)
#else
            call eigen_setup(meg,block_size,nprow,npcol,lda,occ,icolor,ictxt,usermap,icol)
#endif
            if(iprisubmat >= 2) then
               write(nfout,'("eigenlib_setup")')
               write(nfout,'("nprow=",i4,", npcol=",i4,", neg=",i6,", meg=",i6)') nprow,npcol,neg,meg
               write(nfout,'("lda=",i4,", occ=",i4,", block_size=",i4)') lda,occ,block_size
               write(nfout,'("maxeg=",i4,", mp_e=",i4,", nsame=",l4)')   maxeg, mp_e ,nsame
               write(nfout,'("nsclrow=",i4,", nsclcol=",i4)') nsclrow, nsclcol
               call flush(nfout)
            endif
          else
#endif
          if(block_size == 0) then
             block_size = nb_mgs
          end if
#if defined(ASSIGN_G_PREVIOUS) || defined(_ASSIGN_ROW_)
          nsame = .false.
#else
          if(nprow==nrank_g .and. npcol==nrank_e .and. block_size==nb_mgs) then
             nsame = .true.
          else
             nsame = .false.
          endif
#endif

!s        call scalapack_setup_3D(meg,lwork,lrwork,liwork,lda,occ,ictxt,myrow,mycol,desca,descz,usermap)
#ifdef EIGEN_TRANS
          call scalapack_setup_3D(meg,lwork1,lrwork1,liwork1,lwork2,lrwork2,liwork2,lda,occ,&
         & ictxt,myrow,mycol,desca,descb,descz,usermap,icolor,icol,irank_c,irank_r)
#else
          call scalapack_setup_3D(meg,lwork1,lrwork1,liwork1,lwork2,lrwork2,liwork2,lda,occ,&
         & ictxt,myrow,mycol,desca,descb,descz,usermap,icolor,icol)
#endif
          if(iprisubmat >= 2) then
             write(nfout,'("scalapack_setup")')
             write(nfout,'("nprow=",i4,", npcol=",i4,", neg=",i6,", meg=",i6)') nprow,npcol,neg,meg
             write(nfout,'("lda=",i4,", occ=",i4,", block_size=",i4)') lda,occ,block_size
             write(nfout,'("maxeg=",i4,", mp_e=",i4,", nsame=",l4)')   maxeg, mp_e ,nsame
             write(nfout,'("nsclrow=",i4,", nsclcol=",i4)') nsclrow, nsclcol
             write(nfout,'("lwork1=",i9,", lrwork1=",i9,", liwork1=",i9)') lwork1,lrwork1,liwork1
             write(nfout,'("lwork2=",i9,", lrwork2=",i9,", liwork2=",i9)') lwork2,lrwork2,liwork2
             write(nfout,'("scalapack nprow=",i4,", npcol=",i4,", block_size=",i4)') nprow,npcol,block_size
             call flush(nfout)
          endif
#ifdef USE_EIGENLIB
          endif
#endif
!endif ifdef USE_EIGENLIB
!fj --------------------
          if (.not. nsame) then
                                                  __TIMER_SUB_START(1470)
#ifdef EIGEN_TRANS
#ifdef EIGEN_6D
          call make_index_band_for_scalapack(neg, meg, nb_mgs, block_size, nsclrow, nsclcol, usermap, irank_c, irank_r)
#else
          call make_index_band_for_scalapack(neg, meg, nb_mgs, block_size, nsclrow, nsclcol)
#endif
#else
! === DEBUG by tkato 2013/09/18 ================================================
!            integer, save, dimension(9) :: desc0, desc1
!            integer, save               :: ictxt0 ! , ictxt1
! ==============================================================================

#ifdef SINGLE_CONTEXT
!!             call blacs_get(-1, 0, ictxt0)
             ictxt0 = sys2blacs_handle(mpi_k_world(myrank_k))

             call blacs_gridinit(ictxt0, "R", nrank_g, nrank_e)
#else
!!$             call blacs_get(-1, 0, ictxt0(myrank_k))
             ictxt0(myrank_k) = sys2blacs_handle(mpi_k_world(myrank_k))
             call blacs_gridinit(ictxt0(myrank_k), "R", nrank_g, nrank_e)
#endif

             ierr=0
#ifdef SINGLE_CONTEXT
             call descinit(desc0, neg, neg, nb_mgs,     nb_mgs,     0, 0, ictxt0, maxeg, ierr)
#else
             call descinit(desc0, neg, neg, nb_mgs,     nb_mgs,     0, 0, ictxt0(myrank_k), maxeg, ierr)
#endif
             desc1(1) = 1
#ifdef SINGLE_CONTEXT
             desc1(2) = ictxt
#else
             desc1(2) = ictxt(myrank_k)
#endif
             desc1(3) = meg
             desc1(4) = meg
             desc1(5) = block_size
             desc1(6) = block_size
             desc1(7) = 0
             desc1(8) = 0
             desc1(9) = lda
             if(iprisubmat >= 2 ) then
#ifdef SINGLE_CONTEXT
               write(nfout,*) 'usermap=',usermap
#else
               write(nfout,*) 'usermap=',usermap(:,:,myrank_k)
#endif
               write(nfout,'("neg=",i6,", meg=",i6)') neg,meg
               write(nfout,'("nprow  =",i6,", npcol  =",i6)') nprow,npcol
               write(nfout,'("nsclrow=",i6,", nsclcol=",i6)') nsclrow,nsclcol
               write(nfout,'("nb_mgs    =",i6,", maxeg=",i6,", maxe=",i6)') nb_mgs,maxeg,maxe
               write(nfout,'("block_size=",i6,", lda  =",i6,", occ =",i6)') block_size,lda,occ
               write(nfout,'("desc0=",9(1x,i6))') desc0
               write(nfout,'("desc1=",9(1x,i6))') desc1
             end if
#endif
                                                  __TIMER_SUB_STOP(1470)
          end if
       end if
!fj --------------------
!fj    allocate(zmat_l(lda*kimg_t,occ))
       allocate(zmat_l(maxeg*kimg_t,maxe))
!fj --------------------
       zmat_l = 0.0d0
!fj --------------------
!s     allocate(amat_l(maxeg*kimg_t,maxe))
       allocate(amat_l(lda*kimg_t,occ))
!fj --------------------
       amat_l = 0.0d0
       if(sw_gep == ON)then
         allocate(zmat1_l(maxeg*kimg_t,maxe));zmat1_l = 0.d0
         allocate(amat1_l(lda*kimg_t,occ));amat1_l = 0.d0
       endif
    else
       allocate(zmat(meg*kimg_t,meg))
       zmat = 0.0d0
       allocate(zmat_l(maxeg*kimg_t,maxe))
! === DEBUG by Tkato 2011/06/28 ================================================
       lda = maxeg
       occ = maxe
! ==============================================================================
       zmat_l = 0.0d0
       if(sw_gep == ON)then
         allocate(zmat1(meg*kimg_t,meg));zmat1 = 0.d0
         allocate(zmat1_l(maxeg*kimg_t,maxe));zmat1_l = 0.d0
       endif
    end if
#else
    sw_scalapack = OFF
    allocate(zmat(meg*kimg_t,meg))
    zmat = 0.0d0
    allocate(zmat_l(nel_eg(myrank_g)*kimg_t,np_e))
    zmat_l = 0.0d0
    if(sw_gep == ON)then
       allocate(zmat1(meg*kimg_t,meg));zmat1=0.d0
       allocate(zmat1_l(nel_eg(myrank_g)*kimg_t,np_e));zmat1_l=0.d0
    endif
#endif
    if(submat_uncalled) then
       submat_uncalled = .false.
    end if
!    call tstatc0_begin('evolve_WFs_in_subspace_3D(PART1) ', id_sname2)

    allocate(wk_ball(maxval(np_g1k),neg,kimg), stat=ierr)
!!
!!  order of eigenvalue B-block cyclic  ->  order of eigenvalue
!!
                                                  __TIMER_DO_START(919)
    do j = 1, neg
      do i = 1, np_g1k(ik)
         wk_ball(i,neg_g_all(j),1) = zaj_ball(i,j,ik,1)
      enddo
    enddo
    if (kimg == 2) then
      do j = 1, neg
        do i = 1, np_g1k(ik)
          wk_ball(i,neg_g_all(j),2) = zaj_ball(i,j,ik,2)
        enddo
      enddo
    end if
                                                  __TIMER_DO_STOP(919)
!!
!!  order of eigenvalue  ->  order of eigenvalue G-block cyclic
!!
                                                  __TIMER_DO_START(920)
    do j = 1, neg
      do i = 1, np_g1k(ik)
        zaj_ball(i,j,ik,1) = wk_ball(i,neg_gg_all(j),1)
      enddo
    enddo
    if (kimg == 2) then
      do j = 1, neg
        do i = 1, np_g1k(ik)
          zaj_ball(i,j,ik,2) = wk_ball(i,neg_gg_all(j),2)
        enddo
      enddo
    end if
                                                  __TIMER_DO_STOP(920)
    deallocate(wk_ball)
#ifdef __FAPP__
    call fapp_stop('submat_pre',1,1)
#endif

#ifdef __FAPP__
    call fapp_start('submat_part1',1,1)
#endif
    ekod = 0.d0
    sum_abs_diagonal = 0.d0
    sum_sq_diagonal = 0.d0
    sum_sq_non_diagonal = 0.d0
    sum_abs_non_diagonal = 0.d0
    PART1: do jp = 1, npart
       is = isp(jp)
       ie = iep(jp)
       allocate(w1hw2((ie-is+1)*kimg_t,np_e))
       w1hw2 = 0.d0
       if(sw_gep == ON)then
         allocate(w1sw2((ie-is+1)*kimg_t,np_e))
         w1sw2 = 0.d0
       endif
#ifdef SUBMAT_DGEMM
       ibsize = min(ie,neg)-is+1
       imax = maxval(np_g1k(:))
       if(kimg==1) then
                                                  __TIMER_DGEMM_START(921)
          call dgemm('T','N',ibsize,np_e,np_g1k(ik),1.d0,zaj_ball(1,is,ik,1),imax &
               &            ,zaj_l(1,1,ik,1),imax,1.d0,w1hw2(1,1),ibsize)
                                                  __TIMER_DGEMM_STOP(921)
          if(sw_gep == ON)then
               call dgemm('T','N',ibsize,np_e,np_g1k(ik),1.d0,zaj_ball(1,is,ik,1),imax &
               &            ,zat_l(1,1,ik,1),imax,1.d0,w1sw2(1,1),ibsize)
          endif
       else
          if(k_symmetry(ik) == GAMMA) then
                                                  __TIMER_DGEMM_START(922)
             call dgemm('T','N',ibsize,np_e,np_g1k(ik),1.d0,zaj_ball(1,is,ik,1),imax &
                  &         ,zaj_l(1,1,ik,1),imax,1.d0,w1hw2(1,1),ibsize)
             call dgemm('T','N',ibsize,np_e,np_g1k(ik),1.d0,zaj_ball(1,is,ik,2),imax &
                  &         ,zaj_l(1,1,ik,2),imax,1.d0,w1hw2(1,1),ibsize)
                                                  __TIMER_DGEMM_STOP(922)
             w1hw2 = w1hw2*2.d0
             if(myrank_g == 0) then
                                                  __TIMER_DO_START(923)
                do ib1 = is, ie
                   do ib2 = 1, np_e
                      w1hw2(ib1-is+1,ib2) = w1hw2(ib1-is+1,ib2) &
                     &                     - 2.d0*zaj_l(1,ib2,ik,2)*zaj_ball(1,ib1,ik,2) &
                     &                          - zaj_l(1,ib2,ik,1)*zaj_ball(1,ib1,ik,1)
                   end do
                end do
                                                  __TIMER_DO_STOP(923)
             end if
             if(sw_gep == ON)then
               call dgemm('T','N',ibsize,np_e,np_g1k(ik),1.d0,zaj_ball(1,is,ik,1),imax &
                    &         ,zat_l(1,1,ik,1),imax,1.d0,w1sw2(1,1),ibsize)
               call dgemm('T','N',ibsize,np_e,np_g1k(ik),1.d0,zaj_ball(1,is,ik,2),imax &
                    &         ,zat_l(1,1,ik,2),imax,1.d0,w1sw2(1,1),ibsize)
               w1sw2 = w1sw2*2.d0
               if(myrank_g == 0) then
                  do ib1 = is, ie
                     do ib2 = 1, np_e
                        w1sw2(ib1-is+1,ib2) = w1sw2(ib1-is+1,ib2) &
                       &                     - 2.d0*zat_l(1,ib2,ik,2)*zaj_ball(1,ib1,ik,2) &
                       &                          - zat_l(1,ib2,ik,1)*zaj_ball(1,ib1,ik,1)
                     end do
                  end do
               end if
             endif
          else
             allocate(w1hw2r((ie-is+1),np_e))
              if (ierr /= 0) then
                 call mpi_abort(mpi_comm_world, 115, ierr)
              end if
             w1hw2r = 0.d0
             allocate(w1hw2i((ie-is+1),np_e))
              if (ierr /= 0) then
                 call mpi_abort(mpi_comm_world, 116, ierr)
              end if
             w1hw2i = 0.d0
                                                  __TIMER_DGEMM_START(924)
             call dgemm('T','N',ibsize,np_e,np_g1k(ik),1.d0,zaj_ball(1,is,ik,1),imax &
                  &                 , zaj_l(1,1,ik,1),imax,1.d0,w1hw2r(1,1),ibsize)
             call dgemm('T','N',ibsize,np_e,np_g1k(ik),1.d0,zaj_ball(1,is,ik,2),imax &
                  &                 , zaj_l(1,1,ik,2),imax,1.d0,w1hw2r(1,1),ibsize)
             call dgemm('T','N',ibsize,np_e,np_g1k(ik),1.d0,zaj_ball(1,is,ik,2),imax &
                  &                 , zaj_l(1,1,ik,1),imax,1.d0,w1hw2i(1,1),ibsize)
             call dgemm('T','N',ibsize,np_e,np_g1k(ik),-1.d0,zaj_ball(1,is,ik,1),imax &
                  &                 , zaj_l(1,1,ik,2),imax,1.d0,w1hw2i(1,1),ibsize)
                                                  __TIMER_DGEMM_STOP(924)
                                                  __TIMER_DO_START(925)
             do ib1 = is, min(ie,neg)
                do ib2 = 1, np_e
                   w1hw2((ib1-is+1)*2-1,ib2) = w1hw2r(ib1-is+1,ib2)
                   w1hw2((ib1-is+1)*2  ,ib2) = w1hw2i(ib1-is+1,ib2)
                end do
             end do
                                                  __TIMER_DO_STOP(925)
             deallocate(w1hw2i)
             deallocate(w1hw2r)
             if(sw_gep == ON)then
               allocate(w1sw2r((ie-is+1),np_e))
               if (ierr /= 0) then
                  call mpi_abort(mpi_comm_world, 115, ierr)
               end if
               w1sw2r = 0.d0
               allocate(w1sw2i((ie-is+1),np_e))
               if (ierr /= 0) then
                 call mpi_abort(mpi_comm_world, 116, ierr)
               end if
               w1sw2i = 0.d0
               call dgemm('T','N',ibsize,np_e,np_g1k(ik),1.d0,zaj_ball(1,is,ik,1),imax &
                    &                 , zat_l(1,1,ik,1),imax,1.d0,w1sw2r(1,1),ibsize)
               call dgemm('T','N',ibsize,np_e,np_g1k(ik),1.d0,zaj_ball(1,is,ik,2),imax &
                    &                 , zat_l(1,1,ik,2),imax,1.d0,w1sw2r(1,1),ibsize)
               call dgemm('T','N',ibsize,np_e,np_g1k(ik),1.d0,zaj_ball(1,is,ik,2),imax &
                    &                 , zat_l(1,1,ik,1),imax,1.d0,w1sw2i(1,1),ibsize)
               call dgemm('T','N',ibsize,np_e,np_g1k(ik),-1.d0,zaj_ball(1,is,ik,1),imax &
                    &                 , zat_l(1,1,ik,2),imax,1.d0,w1sw2i(1,1),ibsize)
               do ib1 = is, min(ie,neg)
                  do ib2 = 1, np_e
                     w1sw2((ib1-is+1)*2-1,ib2) = w1sw2r(ib1-is+1,ib2)
                     w1sw2((ib1-is+1)*2  ,ib2) = w1sw2i(ib1-is+1,ib2)
                  end do
               end do
               deallocate(w1sw2i)
               deallocate(w1sw2r)
             endif
          end if
       end if
#else
                                                  __TIMER_DO_START(926)
       do ib1 = is, ie
!fj --------------------
          if (neg_g_all(ib1)>meg) cycle
!fj --------------------
          ib1ad = ib1-is+1
          do ib2 = 1, np_e
!fj --------------------
             if (neg_g(ib2)>meg) cycle
!fj --------------------
             if(kimg == 1) then
                do ii = 1, np_g1k(ik)              ! MPI
                   dr1 = zaj_l(ii,ib2,ik,1)
                   dr2 = zaj_ball(ii,ib1,ik,1)
                   w1hw2(ib1ad,ib2) = w1hw2(ib1ad,ib2) + dr1*dr2
                end do
             else
                if(k_symmetry(ik) == GAMMA) then
                   dd = 0.d0
                   ig1 = 1
                   if(myrank_g == 0) ig1 = 2
                   do ii = ig1, np_g1k(ik)            ! MPI
                      dr1 = zaj_l(ii,ib2,ik,1)
                      di1 = zaj_l(ii,ib2,ik,2)
                      dr2 = zaj_ball(ii,ib1,ik,1)
                      di2 = zaj_ball(ii,ib1,ik,2)
                      dd = dd + (dr1*dr2+di1*di2)*2.d0
                   end do
                   if(myrank_g == 0) dd = dd + zaj_l(1,ib2,ik,1)*zaj_ball(1,ib1,ik,1)
                   w1hw2(ib1ad,ib2)=dd
                else
                   do ii = 1, np_g1k(ik)             ! MPI
                      dr1 = zaj_l(ii,ib2,ik,1)
                      di1 = zaj_l(ii,ib2,ik,2)
                      dr2 = zaj_ball(ii,ib1,ik,1)
                      di2 = zaj_ball(ii,ib1,ik,2)
                      w1hw2(ib1ad*2-1,ib2)=w1hw2(ib1ad*2-1,ib2)+dr1*dr2+di1*di2
                      w1hw2(ib1ad*2  ,ib2)=w1hw2(ib1ad*2  ,ib2)+dr1*di2-di1*dr2
                   end do
                end if
             end if
          end do
       end do
                                                  __TIMER_DO_STOP(926)
#endif

!******************** modified by RIST_16
!       call m_ES_S_4_each_k(nfout,ik,is,ie,kimg,meg,w1sw2)

       if(modnrm == EXECUT .and. sw_gep == ON) then
          do ib1 = is, min(ie,neg)
             if (neg_g_all(ib1)>meg) cycle
             ib1ad = ib1-is+1
             do ib2 = 1,np_e
              ibb2 = neg_g(ib2)
              if (ibb2>meg) cycle
                if((k_symmetry(ik) == GAMMA .and. kimg == 2)) then
                   call WSW_for_overlap_t_g(ik,ibb2,neg_gg_all(ib1),fr,fi,fsr_ball(:,:,ik))
                   w1sw2(ib1ad,ib2) = w1sw2(ib1ad,ib2) + fr
                else
                   call WSW_for_overlap_t_g(ik,ibb2,neg_gg_all(ib1),fr,fi,fsr_ball(:,:,ik),fsi_ball(:,:,ik))
                   w1sw2(ib1ad*2-1,ib2) = w1sw2(ib1ad*2-1,ib2) + fr
                   w1sw2(ib1ad*2,  ib2) = w1sw2(ib1ad*2,  ib2) - fi
                end if
             end do
          end do
       end if
!******************** modified by RIST_16
!! (off-diagonal matrix elements will be damped.)
                                                  __TIMER_DO_START(927)
       if(kimg_t == 1) then
!OCL PARALLEL
          do ib1=is,ie
             ib1ad = ib1-is+1
             do ib2=ista_e,iend_e
                ib2ad = ib2-ista_e+1
                if (neg_g(ib2ad) /= neg_gg_all(ib1)) then
                   sum_sq_non_diagonal = sum_sq_non_diagonal + w1hw2(ib1ad,ib2ad)*w1hw2(ib1ad,ib2ad)
                   sum_abs_non_diagonal = sum_abs_non_diagonal + abs(w1hw2(ib1ad,ib2ad))
                   w1hw2(ib1ad,ib2ad) = damp*w1hw2(ib1ad,ib2ad)
                else
                   sum_sq_diagonal = sum_sq_diagonal + w1hw2(ib1ad,ib2ad)**2
                   sum_abs_diagonal = sum_abs_diagonal + abs(w1hw2(ib1ad,ib2ad))
                   ekod = ekod + w1hw2(ib1ad,ib2ad)
                endif
             enddo
          enddo
       else
!OCL PARALLEL
          do ib1=is,ie
             ib1ad = ib1-is+1
             do ib2=ista_e,iend_e
                ib2ad = ib2-ista_e+1
                if (neg_g(ib2ad) /= neg_gg_all(ib1)) then
                   dd = w1hw2(2*ib1ad-1,ib2ad)*w1hw2(2*ib1ad-1,ib2ad) &
                  &   + w1hw2(2*ib1ad,  ib2ad)*w1hw2(2*ib1ad,  ib2ad)
                   sum_sq_non_diagonal = sum_sq_non_diagonal + dd
                   if(dd >= SmallestPositiveNumber) then
                      sum_abs_non_diagonal = sum_abs_non_diagonal + sqrt(dd)
                   endif
                   w1hw2(2*ib1ad-1,ib2ad) = damp*w1hw2(2*ib1ad-1,ib2ad)
                   w1hw2(2*ib1ad  ,ib2ad) = damp*w1hw2(2*ib1ad  ,ib2ad)
                else
                   sum_sq_diagonal = sum_sq_diagonal + w1hw2(2*ib1ad-1,ib2ad)**2
                   sum_abs_diagonal = sum_abs_diagonal + abs(w1hw2(2*ib1ad-1,ib2ad))
                   ekod = ekod + w1hw2(2*ib1ad-1,ib2ad)
                endif
             enddo
          enddo
       endif
                                                  __TIMER_DO_STOP(927)

#ifdef _USE_SCALAPACK_
       if(sw_scalapack==ON) then
          call set_hmat_3D_scl(neg,w1hw2,zmat_l)
          if(sw_gep == ON)then
             call set_hmat_3D_scl(neg,w1sw2,zmat1_l)
          endif
       else
          call set_hmat_3D(neg,w1hw2,zmat_l)
          if(sw_gep == ON)then
             call set_hmat_3D(neg,w1sw2,zmat1_l)
          endif
       end if
#else
       call set_hmat_3D(neg,w1hw2,zmat_l)
       if(sw_gep == ON)then
          call set_hmat_3D(neg,w1sw2,zmat1_l)
       endif
#endif
       deallocate(w1hw2)
       if(sw_gep==ON)then
          deallocate(w1sw2)
       endif
    end do PART1

#ifdef __FAPP__
    call fapp_stop('submat_part1',1,1)
#endif

#ifdef __FAPP__
    call fapp_start('submat_scalapack',1,1)
#endif
                                                  __TIMER_STOP_w_BARRIER(mpi_k_world(myrank_k),8)
                                                  __TIMER_START(9)
                                                  __TIMER_START(1500)
    if(.not.dry) then
    if(sw_scalapack == ON) then
       if(nsame) then
                                                  __TIMER_DO_START(981)
          if (kimg_t == 1) then
             do j = 1, np_e
                if (neg_g(j) > meg) cycle
                do i = 1, nel_eg(myrank_g)
                   if (neg_gg(i) > meg) cycle
                   amat_l(i,j) = zmat_l(i,j)
                enddo
                if(sw_gep == ON)then
                   do i = 1, nel_eg(myrank_g)
                      if (neg_gg(i) > meg) cycle
                      amat1_l(i,j) = zmat1_l(i,j)
                   enddo
                endif
             enddo
          else
             do j = 1, np_e
                if (neg_g(j) > meg) cycle
                do i = 1, nel_eg(myrank_g)
                   if (neg_gg(i) > meg) cycle
                   amat_l(i*2-1,j) = zmat_l(i*2-1,j)
                   amat_l(i*2  ,j) = zmat_l(i*2  ,j)
                enddo
                if(sw_gep == ON)then
                   do i = 1, nel_eg(myrank_g)
                      if (neg_gg(i) > meg) cycle
                      amat1_l(i*2-1,j) = zmat1_l(i*2-1,j)
                      amat1_l(i*2  ,j) = zmat1_l(i*2  ,j)
                   enddo
                endif
             enddo
          end if
                                                  __TIMER_DO_STOP(981)
       else
!fj --------------------s
                                                  __TIMER_SUB_START(964)
#ifdef EIGEN_TRANS
#ifdef EIGEN_6D
          call trans_scalapack(neg, meg, nb_mgs, block_size, zmat_l, amat_l, &
         &                     maxeg, maxe, lda, occ, kimg_t, usermap, irank_c, irank_r)
          if(sw_gep == ON)then
             call trans_scalapack(neg, meg, nb_mgs, block_size, zmat1_l, amat1_l, &
            &                     maxeg, maxe, lda, occ, kimg_t, usermap, irank_c, irank_r)
          endif
#else
          call trans_scalapack(neg, meg, nb_mgs, block_size, zmat_l, amat_l, &
         &                     maxeg, maxe, lda, occ, kimg_t)
          if(sw_gep == ON)then
             call trans_scalapack(neg, meg, nb_mgs, block_size, zmat1_l, amat1_l, &
            &                     maxeg, maxe, lda, occ, kimg_t)
          endif
#endif
#else
#ifdef SINGLE_CONTEXT
          if(kimg_t == 1) then
            call pdgemr2d(meg, meg, zmat_l, 1, 1, desc0, amat_l, 1, 1, desc1, ictxt0)
            if(sw_gep == ON)then
               call pdgemr2d(meg, meg, zmat1_l, 1, 1, desc0, amat1_l, 1, 1, desc1, ictxt0)
            endif
          else
            call pzgemr2d(meg, meg, zmat_l, 1, 1, desc0, amat_l, 1, 1, desc1, ictxt0)
            if(sw_gep == ON)then
               call pzgemr2d(meg, meg, zmat1_l, 1, 1, desc0, amat1_l, 1, 1, desc1, ictxt0)
            endif
          endif
#else
          if(kimg_t == 1) then
            call pdgemr2d(meg, meg, zmat_l, 1, 1, desc0, amat_l, 1, 1, desc1, ictxt0(myrank_k))
            if(sw_gep == ON)then
               call pdgemr2d(meg, meg, zmat1_l, 1, 1, desc0, amat1_l, 1, 1, desc1, ictxt0(myrank_k))
            endif
          else
            call pzgemr2d(meg, meg, zmat_l, 1, 1, desc0, amat_l, 1, 1, desc1, ictxt0(myrank_k))
            if(sw_gep == ON)then
               call pzgemr2d(meg, meg, zmat1_l, 1, 1, desc0, amat1_l, 1, 1, desc1, ictxt0(myrank_k))
            endif
          endif
#endif
#endif
                                                  __TIMER_SUB_STOP(964)
!!$#if 0
!!$          allocate(zmat(neg*kimg_t,neg),stat=j)
!!$          zmat = 0.0d0
!!$          call gather_zmat_all(zmat_l,zmat)
!!$          call get_zmat(lda,occ,amat_l,neg,nsclrow,nsclcol,block_size,usermap,1,neg,zmat)
!!$          deallocate(zmat)
!!$#else
!!$          call trans_scalapack(neg, meg, nb_mgs, block_size, zmat_l, amat_l, &
!!$         &                     maxeg, maxe, lda, occ, kimg_t)
!!$#endif
!fj --------------------t
       end if
       deallocate(zmat_l)
       allocate(zmat_l(lda*kimg_t,occ))
       zmat_l = 0.0d0
       if(sw_gep == ON)then
          deallocate(zmat1_l)
          allocate(zmat1_l(lda*kimg_t,occ))
          zmat1_l = 0.0d0
       endif
    else
       call gather_zmat_all(zmat_l,zmat)
       if(sw_gep == ON)then
          call gather_zmat_all(zmat1_l,zmat1)
       endif
    endif
    endif

    sum_abs_diagonal = sum_abs_diagonal/meg
    sum_abs_non_diagonal = sum_abs_non_diagonal/(meg*(meg-1)/2)

!! (Diagonalization )  !!

     if(.not.dry)then
#ifdef _USE_SCALAPACK_
     if(sw_scalapack == ON) then
                                                  __TIMER_STOP_w_BARRIER(mpi_k_world(myrank_k),1500)
#ifdef USE_EIGENLIB
     if(sw_eigen_exa==ON) then
       if (icolor == 1) then
          if(kimg_t == 1) then
             call eigen_solver(lda,occ)
          else
!!!!#ifdef EIGEN_EXA_2_9
!!!!             call phase_error_with_msg(nfout,'EigenExa supports only real symmetric matrices'&
!!!!             ,__LINE__,__FILE__)
!!!!#else
             call eigen_solver_h(lda,occ)
!!!!#endif
          endif
       endif
     else
#endif
       if (myrow /= -1) then
          if(sw_gep == ON)then
             if(kimg_t == 1) then
                call pdsygvx_driver_3D(meg,eig,desca,descb,descz,lda,occ,amat_l,amat1_l,zmat_l,lwork1,liwork1)
             else
                call pzhegvx_driver_3D(meg,eig,desca,descb,descz,lda,occ,amat_l,amat1_l,zmat_l,lwork2,lrwork2,liwork2)
             end if
          else
             if(kimg_t == 1) then
                call pdsyev_driver_3D(meg,eig,desca,descz,lda,occ,amat_l,zmat_l,lwork1,liwork1)
             else
                call pzheev_driver_3D(meg,eig,desca,descz,lda,occ,amat_l,zmat_l,lwork2,lrwork2,liwork2)
             end if
          endif
       end if
#ifdef USE_EIGENLIB
     endif
#endif
                                                  __TIMER_START_w_BARRIER(mpi_k_world(myrank_k),1500)
        if (nrank_e*nrank_g /= nsclcol*nsclrow) then
           call eigsend()
        end if
! ** dummy ** !
                                                  __TIMER_SUB_START(998)
                                                  __TIMER_SUB_STOP(998)
                                                  __TIMER_COMM_START(989)
                                                  __TIMER_COMM_START(990)
                                                  __TIMER_COMM_STOP(989)
                                                  __TIMER_COMM_STOP(990)
     else
        if(sw_gep == ON)then
           if(kimg_t == 1) then
              call dsygvx_driver(meg,eig,zmat,zmat1)
           else
              call zhegvx_driver(meg,eig,zmat,zmat1)
           endif
        else
           if(kimg_t == 1) then
              call dsyev_driver(meg,eig,zmat)
           else
              call zheev_driver(meg,eig,zmat)
           endif
        endif
     endif
#else
     if(kimg_t == 1) then
       call dsyev_driver(meg,eig,zmat)
     else
       call zheev_driver(meg,eig,zmat)
     endif
#endif
     endif

     if (sw_scalapack == ON) then
        if (.not.nsame) then
!fj --------------------s
#if 0
           allocate(zmat(neg*kimg_t,neg),stat=j)
           zmat = 0.0d0
           call put_zmat(kimg_t,lda,maxe,zmat_l,neg,nsclrow,nsclcol,block_size,usermap,1,neg,zmat) ! zmat_l->zmat
           allocate(zmat_l_2D(neg*kimg_t,np_e), stat=ierr)
           zmat_l_2D = 0.0d0
           if(kimg_t == 1) then
              do i = 1, np_e
                 if (neg_g(i) > meg) cycle
                 do j = 1, neg
                    if (neg_g_all(j) > meg) cycle
                    zmat_l_2D(j,i) = zmat(neg_g_all(j),neg_g(i))
                 enddo
              end do
           else
              do i = 1, np_e
                 if (neg_g(i) > meg) cycle
                 do j = 1, neg
                    if (neg_g_all(j) > meg) cycle
                    zmat_l_2D(j*2-1,i) = zmat(neg_g_all(j)*2-1,neg_g(i))
                    zmat_l_2D(j*2  ,i) = zmat(neg_g_all(j)*2  ,neg_g(i))
                 enddo
              end do
           endif
           deallocate(zmat)
#else
           deallocate(amat_l)
           allocate(amat_l(maxeg*kimg_t,maxe))
           amat_l = 0.0d0
           if(sw_gep == ON)then
              deallocate(amat1_l)
              allocate(amat1_l(maxeg*kimg_t,maxe))
              amat1_l = 0.0d0
           endif
                                                  __TIMER_SUB_START(965)
#ifdef EIGEN_TRANS
#ifdef EIGEN_6D
           call trans_scalapack_r(neg, meg, nb_mgs, block_size, zmat_l, amat_l, maxeg, maxe, lda, occ, kimg_t, usermap, irank_c, irank_r)
#else
           call trans_scalapack_r(neg, meg, nb_mgs, block_size, zmat_l, amat_l, maxeg, maxe, lda, occ, kimg_t)
#endif
#else
#ifdef SINGLE_CONTEXT
           if(kimg_t == 1) then
              call pdgemr2d(meg, meg, zmat_l, 1, 1, desc1, amat_l, 1, 1, desc0, ictxt0)
              if(sw_gep == ON)then
                 call pdgemr2d(meg, meg, zmat1_l, 1, 1, desc1, amat1_l, 1, 1, desc0, ictxt0)
              endif
           else
              call pzgemr2d(meg, meg, zmat_l, 1, 1, desc1, amat_l, 1, 1, desc0, ictxt0)
              if(sw_gep == ON)then
                 call pzgemr2d(meg, meg, zmat1_l, 1, 1, desc1, amat1_l, 1, 1, desc0, ictxt0)
              endif
           endif
#endif
           if(kimg_t == 1) then
#ifdef SINGLE_CONTEXT
              call pdgemr2d(meg, meg, zmat_l, 1, 1, desc1, amat_l, 1, 1, desc0, ictxt0)
#else
              call pdgemr2d(meg, meg, zmat_l, 1, 1, desc1, amat_l, 1, 1, desc0, ictxt0(myrank_k))
#endif
              if(sw_gep == ON)then
#ifdef SINGLE_CONTEXT
                 call pdgemr2d(meg, meg, zmat1_l, 1, 1, desc1, amat1_l, 1, 1, desc0, ictxt0)
#else
                 call pdgemr2d(meg, meg, zmat1_l, 1, 1, desc1, amat1_l, 1, 1, desc0, ictxt0(myrank_k))
#endif
              endif
           else
#ifdef SINGLE_CONTEXT
              call pzgemr2d(meg, meg, zmat_l, 1, 1, desc1, amat_l, 1, 1, desc0, ictxt0)
#else
              call pzgemr2d(meg, meg, zmat_l, 1, 1, desc1, amat_l, 1, 1, desc0, ictxt0(myrank_k))
#endif
              if(sw_gep == ON)then
#ifdef SINGLE_CONTEXT
                 call pzgemr2d(meg, meg, zmat1_l, 1, 1, desc1, amat1_l, 1, 1, desc0, ictxt0)
#else
                 call pzgemr2d(meg, meg, zmat1_l, 1, 1, desc1, amat1_l, 1, 1, desc0, ictxt0(myrank_k))
#endif
              endif
           endif
#endif
                                                  __TIMER_SUB_STOP(965)
           allocate(wk_gather(maxeg*kimg_t, maxe, 0:nrank_g-1), stat=ierr)
           allocate(wk_mpi(neg*kimg_t,maxe), stat=ierr)
                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,982)
           call mpi_allgather(amat_l,    maxeg*kimg_t*maxe, mpi_double_precision,&
          &                   wk_gather, maxeg*kimg_t*maxe, mpi_double_precision,&
          &                   mpi_ke_world, ierr)
                                                  __TIMER_COMM_STOP(982)
                                                  __TIMER_DO_START(983)
           if (kimg_t == 1) then
              do k = 0, nrank_g-1
                 do j = 1, np_e
                    do i = nis_eg(k), nie_eg(k)
                       wk_mpi(neg_gg_all(i),j) = wk_gather(i-nis_eg(k)+1,j,k)
                    enddo
                 enddo
              enddo
           else
              do k = 0, nrank_g-1
                 do j = 1, np_e
                    do i = nis_eg(k), nie_eg(k)
                       wk_mpi(neg_gg_all(i)*2-1,j) = wk_gather((i-nis_eg(k)+1)*2-1,j,k)
                       wk_mpi(neg_gg_all(i)*2  ,j) = wk_gather((i-nis_eg(k)+1)*2  ,j,k)
                    enddo
                 enddo
              enddo
           endif
                                                  __TIMER_DO_STOP(983)
           deallocate(wk_gather)
           allocate(zmat_l_2D(neg*kimg_t,np_e), stat=ierr)
                                                  __TIMER_DO_START(984)
           if (kimg_t == 1) then
              do j = 1, np_e
                 do i = 1, neg
                    zmat_l_2D(i,j) = wk_mpi(neg_g_all(i),j)
                 enddo
              enddo
           else
              do j = 1, np_e
                 do i = 1, neg
                    zmat_l_2D(i*2-1,j) = wk_mpi(neg_g_all(i)*2-1,j)
                    zmat_l_2D(i*2  ,j) = wk_mpi(neg_g_all(i)*2  ,j)
                 enddo
              enddo
           endif
                                                  __TIMER_DO_STOP(984)
           deallocate(wk_mpi)
#endif
!fj --------------------e
        else
! ( zmat_l  G b-cyclic -> B b-cyclic )
           allocate(wk_gather(maxeg*kimg_t, maxe, 0:nrank_g-1), stat=ierr)
           allocate(wk_mpi(maxeg*kimg_t, maxe), stat=ierr)
                                                  __TIMER_DO_START(984)
           wk_mpi = 0.0d0
           if (kimg_t == 1) then
              do j = 1,np_e
                 if (neg_g(j) > meg) cycle
                 do i = 1,nel_eg(myrank_g)
                   if (neg_gg(i) > meg) cycle
                    wk_mpi(i,j) = zmat_l(i,j)
                 enddo
              enddo
           else
              do j = 1,np_e
                 if (neg_g(j) > meg) cycle
                 do i = 1,nel_eg(myrank_g)
                    if (neg_gg(i) > meg) cycle
                    wk_mpi(i*2-1,j) = zmat_l(i*2-1,j)
                    wk_mpi(i*2  ,j) = zmat_l(i*2  ,j)
                 enddo
              enddo
           endif
                                                  __TIMER_DO_STOP(984)
                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,987)
           call mpi_allgather(wk_mpi,   maxeg*kimg_t*maxe, mpi_double_precision,&
          &                   wk_gather, maxeg*kimg_t*maxe, mpi_double_precision,&
          &                   mpi_ke_world, ierr)
                                                  __TIMER_COMM_STOP(987)

           deallocate(wk_mpi)
           allocate(wk_mpi(neg*kimg_t,maxe), stat=ierr)
                                                  __TIMER_DO_START(985)
           wk_mpi = 0.0d0
           if (kimg_t == 1) then
              do k = 0, nrank_g-1
                 do j = 1, np_e
                    do i = nis_eg(k), nie_eg(k)
                       wk_mpi(neg_gg_all(i),j) = wk_gather(i-nis_eg(k)+1,j,k)
                    enddo
                 enddo
              enddo
           else
              do k = 0, nrank_g-1
                 do j = 1, np_e
                    do i = nis_eg(k), nie_eg(k)
                       wk_mpi(neg_gg_all(i)*2-1,j) = wk_gather((i-nis_eg(k)+1)*2-1,j,k)
                       wk_mpi(neg_gg_all(i)*2  ,j) = wk_gather((i-nis_eg(k)+1)*2  ,j,k)
                    enddo
                 enddo
              enddo
           endif
                                                  __TIMER_DO_STOP(985)
           deallocate(wk_gather)
           allocate(zmat_l_2D(neg*kimg_t,np_e), stat=ierr)
                                                  __TIMER_DO_START(986)
           if (kimg_t == 1) then
              do j = 1, np_e
                 do i = 1, neg
                    zmat_l_2D(i,j) = wk_mpi(neg_g_all(i),j)
                 enddo
              enddo
           else
              do j = 1, np_e
                 do i = 1, neg
                    zmat_l_2D(i*2-1,j) = wk_mpi(neg_g_all(i)*2-1,j)
                    zmat_l_2D(i*2  ,j) = wk_mpi(neg_g_all(i)*2  ,j)
                 enddo
              enddo
           endif
                                                  __TIMER_DO_STOP(986)
           deallocate(wk_mpi)
        end if
     else
       allocate(zmat_l_2D(neg*kimg_t,np_e), stat=ierr)
       zmat_l_2D = 0.0d0
                                                  __TIMER_DO_START(928)
       if(kimg_t == 1) then
          do i = 1, np_e
             if (neg_g(i) > meg) cycle
             do j = 1, neg
                if (neg_g_all(j) > meg) cycle
                zmat_l_2D(j,i) = zmat(neg_g_all(j),neg_g(i))
             enddo
          end do
       else
          do i = 1, np_e
             if (neg_g(i) > meg) cycle
             do j = 1, neg
                if (neg_g_all(j) > meg) cycle
                zmat_l_2D(j*2-1,i) = zmat(neg_g_all(j)*2-1,neg_g(i))
                zmat_l_2D(j*2  ,i) = zmat(neg_g_all(j)*2  ,neg_g(i))
             enddo
          end do
       endif
                                                  __TIMER_DO_STOP(928)
    end if

! ( zaj_ball  Gb-cyclic -> Bb-cyclic )
    allocate(wk_ball(np_g1k(ik),neg,kimg), stat=ierr)
    wk_ball = 0.0d0
                                                  __TIMER_DO_START(929)
!OCL SERIAL
    do k = 1, kimg
!OCL PARALLEL
       do j = 1, neg
          do i = 1, np_g1k(ik)
             wk_ball(i,neg_gg_all(j),k) = zaj_ball(i,j,ik,k)
          enddo
       enddo
    enddo
                                                  __TIMER_DO_STOP(929)
                                                  __TIMER_DO_START(930)
!OCL SERIAL
    do k = 1, kimg
!OCL PARALLEL
       do j = 1, neg
          do i = 1, np_g1k(ik)
             zaj_ball(i,j,ik,k) = wk_ball(i,neg_g_all(j),k)
          enddo
       enddo
    enddo
                                                  __TIMER_DO_STOP(930)
#ifdef __FAPP__
    call fapp_stop('submat_scalapack',1,1)
#endif
#ifdef __FAPP__
    call fapp_start('submat_rotation',1,1)
#endif
    zaj_l(:,:,ik,:) = 0.0d0
    if(sw_keep_hloc_phi==ON) hlocphi_l(:,:,ik,:) = 0.0d0
    fsr_l(:,:,ik) = 0.d0
    if(k_symmetry(ik) /= GAMMA) fsi_l(:,:,ik) = 0.d0
    PART2: do jp = 1, npart2
       is = isp2(jp)
       ie = iep2(jp)
       if(kimg == 1) then
          call subspace_rotation_real_3D(is,ie,neg)
       else
          call subspace_rotation_imag_3D(is,ie,neg)
       endif
    end do PART2

    deallocate(zmat_l_2D)
    deallocate(fsr_in)
    if(k_symmetry(ik) /= GAMMA) deallocate(fsi_in)
    if (meg < neg) then
       allocate(wk_meg(1:np_g1k(ik),meg+1:neg,kimg), stat=ierr)
                                                  __TIMER_DO_START(931)
!OCL SERIAL
       do k = 1, kimg
          do j = meg+1, neg
!OCL PARALLEL
             do i = 1, np_g1k(ik)
                wk_meg(i,j,k) = wk_ball(i,j,k)
             end do
          end do
       end do
                                                  __TIMER_DO_STOP(931)
       wk_ball = 0.0d0
                                                  __TIMER_DO_START(932)
!OCL SERIAL
       do k = 1, kimg
!OCL PARALLEL
          do j = 1, np_e
             do i = 1, np_g1k(ik)
                wk_ball(i,neg_g(j),k) = zaj_l(i,j,ik,k)
             end do
          end do
       end do
                                                  __TIMER_DO_STOP(932)

       allocate(wk_gather(np_g1k(ik),neg,kimg), stat=ierr)
                                                  __TIMER_COMM_START_w_BARRIER(mpi_kg_world,933)
       call mpi_allreduce(wk_ball,wk_gather,np_g1k(ik)*neg*kimg,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
                                                  __TIMER_COMM_STOP(933)
                                                  __TIMER_DO_START(934)
!OCL SERIAL
       do k = 1, kimg
!OCL PARALLEL
          do j = meg+1, neg
             do i = 1, np_g1k(ik)
                wk_gather(i,j,k) = wk_meg(i,j,k)
             end do
          end do
       end do
                                                  __TIMER_DO_STOP(934)
                                                  __TIMER_DO_START(935)
!OCL SERIAL
       do k = 1, kimg
!OCL PARALLEL
          do j = 1, np_e
             if (neg_g(j) > meg) then
                do i = 1, np_g1k(ik)
                   zaj_l(i,j,ik,k) = wk_gather(i,neordr(neg_g(j),ik),k)
                end do
             end if
          end do
       end do
#ifdef SAVE_FFT_TIMES
       if(sw_save_fft == ON) then
          do j = 1, np_e
             status_saved_phifftr(j,ik) = OLD
          end do
       end if
#endif

                                                  __TIMER_DO_STOP(935)
       deallocate(wk_meg)
       deallocate(wk_gather)
    end if

    deallocate(wk_ball)
#ifdef __FAPP__
    call fapp_stop('submat_rotation',1,1)
#endif
!!  maxe = maxval(nel_e(:))
!!  maxeg = maxval(nel_g1k(:,ik))
!!  allocate(wk_mpi(maxeg*kimg,maxe), stat=ierr)
!!  allocate(wk_gather(maxeg*kimg,maxe,0:nrank_e-1), stat=ierr)
!!  wk_mpi = 0.0d0
!!  if (kimg == 1) then
!!     do j = 1, np_e
!!        do i = 1, np_g1k(ik)
!!           wk_mpi(i,j) = zaj_l(i,j,ik,1)
!!        enddo
!!     enddo
!!  else
!!     do j = 1, np_e
!!        do i = 1, np_g1k(ik)
!!           wk_mpi(i*2-1,j) = zaj_l(i,j,ik,1)
!!           wk_mpi(i*2  ,j) = zaj_l(i,j,ik,2)
!!        enddo
!!     enddo
!!  endif
!!
!!  call mpi_allgather(wk_mpi,   maxeg*kimg_t*maxe, mpi_double_precision,&
!! &                   wk_gather, maxeg*kimg_t*maxe, mpi_double_precision,&
!! &                   mpi_kg_world, ierr)
!!
!!  if (kimg_t == 1) then
!!     do k = 0, nrank_e-1
!!        do j = nis_e(k), nie_e(k)
!!           do i = 1, np_g1k(ik)
!!              zaj_ball(i,j,ik,1) = wk_gather(i,j-nis_e(k)+1,k)
!!           enddo
!!        enddo
!!     enddo
!!  else
!!     do k = 0, nrank_e-1
!!        do j = nis_e(k), nie_e(k)
!!           do i = 1, np_g1k(ik)
!!              zaj_ball(i,j,ik,1) = wk_gather(i*2-1,j-nis_e(k)+1,k)
!!              zaj_ball(i,j,ik,2) = wk_gather(i*2  ,j-nis_e(k)+1,k)
!!           enddo
!!        enddo
!!     enddo
!!  endif
!!  deallocate(wk_mpi)
!!  deallocate(wk_gather)

!! (eko_l)
                                                  __TIMER_DO_START(936)
     if(meg < neg) then
        do ib1=meg+1,neg
           if(map_e(ib1) == myrank_e) then        ! MPI
!             eko_l(map_z(ib1),ik)=eko_d_mpi(neordr(ib1,ik))
              eko_l(map_z(ib1),ik)=eko_d(neordr(ib1,ik))
           end if
        enddo
     endif
                                                  __TIMER_DO_STOP(936)
                                                  __TIMER_DO_START(937)
     do ib1 = 1, meg
        if(map_e(ib1) == myrank_e) then         ! MPI
           eko_l(map_z(ib1),ik)=eig(ib1)
        end if
     end do
                                                  __TIMER_DO_STOP(937)
#ifdef _USE_SCALAPACK_
!!   if(sw_scalapack==ON) then
!!   else
!!   end if
#else
                                                  __TIMER_DO_START(938)
     if(meg < neg) then
        do ib1 = 1, np_e
           if(neg_g(ib1) > meg) then
!             eko_l(ib1,ik)=eko_d_mpi(neordr(neg_g(ib1),ik))
              eko_l(ib1,ik)=eko_d(neordr(neg_g(ib1),ik))
           end if
        enddo
     endif
                                                  __TIMER_DO_STOP(938)
                                                  __TIMER_DO_START(939)
     do ib1 = 1, np_e
        if(neg_g(ib1) > meg) cycle
        eko_l(ib1,ik) = eig(neg_g(ib1))
     end do
                                                  __TIMER_DO_STOP(939)
#endif

!! (neordr & nrvf_ordr)
     neordr(1:neg,ik) = (/(ib1,ib1=1,neg)/)
     nrvf_ordr(1:neg,ik) = (/(ib1,ib1=1,neg)/)

! (deallocate)
    deallocate(eko_d)
    deallocate(eig)
    if(sw_scalapack==ON) deallocate(amat_l)
!fj --------------------
!s  deallocate(zmat)
    if (allocated(zmat)) deallocate(zmat)
!fj --------------------
    deallocate(zmat_l)
    if(sw_gep == ON)then
      deallocate(zmat1_l)
      if(allocated(zmat1)) deallocate(zmat1)
      if(sw_scalapack == ON)then
         deallocate(amat1_l)
      endif
      deallocate(zat_l)
    endif

! === DEBUG by tkato 2012/12/11 ================================================
    call mpi_allreduce(MPI_IN_PLACE,sum_abs_non_diagonal,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
    call mpi_allreduce(MPI_IN_PLACE,sum_abs_diagonal,    1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
! ==============================================================================
    if(.not.non_diagonal_part_is_small(ik)  &
         & .and. sum_abs_non_diagonal/sum_abs_diagonal <= submat_critical_ratio) &
         & non_diagonal_part_is_small(ik) = .true.

    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(901)
                                                  __TIMER_STOP(1500)
                                                  __TIMER_STOP(9)
  contains

#ifdef USE_EIGENLIB

#ifdef EIGEN_TRANS
  subroutine eigen_setup(ndim,nb,nprow,npcol,lda,occ,icolor,ictxt,usermap,icol,irank_c,irank_r)
#else
  subroutine eigen_setup(ndim,nb,nprow,npcol,lda,occ,icolor,ictxt,usermap,icol)
#endif

#ifdef EIGEN_EXA
    use eigen_libs,  only: eigen_init
    use eigen_blacs, only: eigen_set_blacs_context
#endif
#ifdef EIGEN_EXA_2_9
    use eigen_libs_mod,  only: eigen_init
    use eigen_blacs_mod, only: eigen_get_blacs_context
#endif
    implicit none
    integer, intent(in)  :: ndim, nb, nprow, npcol
    integer, intent(out) :: lda, occ, icolor
    integer, intent(out) :: ictxt, usermap(nprow,npcol)
    integer, intent(out) :: icol(0:nrank_e*nrank_g*nrank_k-1)
    integer              :: ikey(0:nrank_e*nrank_g*nrank_k-1)
#ifdef EIGEN_TRANS
    integer, intent(out) :: irank_c(0:nrank_e*nrank_g*nrank_k-1)
    integer, intent(out) :: irank_r(0:nrank_e*nrank_g*nrank_k-1)
#endif
    integer :: myrow, mycol
    integer :: i, k
    integer :: id_sname = -1
    integer :: iam, nprocs, info

    include 'commtxt.h'

                                                  __TIMER_SUB_START(903)
    call tstatc0_begin('eigen_setup ', id_sname)

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

#ifdef EIGEN_TRANS
    call make_usermap(usermap,icol,ikey,irank_c,irank_r)
#else
    call make_usermap(usermap,icol,ikey)
#endif

    call MPI_COMM_RANK(mpi_k_world(myrank_k), iam, ierr)
    call blacs_setup(iam,nprocs)

#ifdef EIGEN_6D
    icolor=icol(iam)
    call MPI_comm_split(MPI_COMM_WORLD,icolor,ikey(iam),MPI_COMM_EIGEN,ierr)
#else
    if ( iam < nprow*npcol ) then
      icolor=1
    else
      icolor=0
    endif
    call MPI_comm_split(mpi_k_world(myrank_k),icolor,iam,MPI_COMM_EIGEN,ierr)
#endif

#ifdef EIGEN_EXA
    call eigen_init( mpi_comm_eigen, order='C' )
    call blacs_get(-1,0,ictxt)
    call blacs_gridmap(ictxt,usermap,nprow,nprow,npcol)
    call eigen_set_blacs_context(ictxt)
#elif EIGEN_EXA_2_9
    call eigen_init(mpi_comm_eigen, order='C')
    ictxt = eigen_get_blacs_context()
#else
    call blacs_get(-1,0,ictxt)
    call blacs_gridmap(ictxt,usermap,nprow,nprow,npcol)
    call blacs_gridinfo(ictxt,nprow,npcol,myrow,mycol)

    ICTXTtoDC=ictxt
    MYROWtoDC=myrow
    MYCOLtoDC=mycol
! === Specify NPROWtoDC and NPCOLtoDC in commtxt.h! ============================
    NPROWtoDC = nprow
    NPCOLtoDC = npcol
! ==============================================================================
#endif

    if(iprisubmat>=2) then
      write(nfout,'("nb=",i6,", lda=",i6,", occ=",i6,", ndim=",i6,", nprow=",i6,", npcol=",i6)') nb, lda, occ, ndim, nprow, npcol
      if(iprisubmat>=2) write(nfout,'("myrow,mycol=",2i8)') myrow,mycol
      call flush(nfout)
    end if

    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(903)
  end subroutine eigen_setup

  subroutine eigen_solver(nn1,nn2)
#ifdef EIGEN_EXA
       use eigen_libs, only: eigen_sx
#endif
#ifdef EIGEN_EXA_2_9
       use eigen_libs_mod, only: eigen_sx, eigen_get_matdims
#endif
    implicit none
    integer :: ierr
    integer,intent(in) :: nn1, nn2
    integer :: nm_cache1,nm_cache2, nm_block
    real(kind=DP), allocatable :: utmp(:,:), vtmp(:,:)
    integer :: NB, nmz, nmw
    integer :: LLD_RX,LLD_CX,MXLLDX,LLD_R,LLD_C,MXLLD
    integer :: iam, lddz, lddw
    integer :: nx, ny
    integer :: id_sname=-1

    call tstatc0_begin('eigen_solver', id_sname)
    if(iprisubmat>=2) write(nfout,*) "Submat: eigen_solver"
    call eigen_get_matdims(meg, nx, ny)
    allocate (utmp(nx,ny), vtmp(nx,ny) )
    utmp(1:nn1,1:nn2) = amat_l(1:nn1,1:nn2)
    vtmp = 0.d0
    call flush(nfout)
#if defined(EIGEN_EXA) || defined(EIGEN_EXA_2_9)
    call eigen_sx(meg, meg, utmp, nx, eig, vtmp, nx)
#else
    call eigen_sx(meg, utmp(1), nn1, eig(1), vtmp(1), nn1, 32)
#endif
    zmat_l(1:nn1,1:nn2) = vtmp(1:nn1,1:nn2)

    deallocate (utmp, vtmp)
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(906)
    return
  end subroutine eigen_solver

#ifdef EIGEN_EXA_2_9
  subroutine eigen_solver_h(nn1,nn2)
    use eigen_libs_mod, only: eigen_h, eigen_get_matdims
    integer,intent(in) :: nn1, nn2
    integer :: ierr, nx, ny
    complex(kind=DP), allocatable, dimension(:,:) :: utmp,vtmp
    integer :: id_sname=-1
    call tstatc0_begin('eigen_solver_h', id_sname)
    call eigen_get_matdims(meg, nx, ny)
    allocate (utmp(nx,ny))
    allocate (vtmp(nx,ny))
    utmp = dcmplx(0.d0,0.d0)
    call mset_h(amat_l,utmp,nn1,nn2,nx,ny)
    call eigen_h(meg, meg, utmp, nx, eig, vtmp, nx)
    call mset2_h(zmat_l, vtmp, nn1, nn2, nx, ny)
    deallocate(utmp,vtmp)
    call tstatc0_end(id_sname)

  end subroutine eigen_solver_h

  subroutine mset_h(a,b,n1,n2,nx,ny)
    implicit none
    integer, intent(in) :: n1,n2,nx,ny
    real(kind=DP),intent(in)    :: a(n1*kimg_t,n2)
    complex(kind=DP),intent(inout) :: b(nx,ny)
    integer i,j

    l=0
    do j=1,n2
       do i=1,n1
          b(i,j)= dcmplx(a(i*2-1,j), a(i*2,j))
       enddo
    enddo

!    if (l>n3) write(*,*) "***** error in mset"

    return
  end subroutine mset_h

  subroutine mset2_h(a,b,n1,n2,nx,ny)
    implicit none
    integer,intent(in) :: n1,n2,nx,ny
    real(kind=DP),intent(inout) :: a(n1*kimg_t,n2)
    complex(kind=DP),intent(in) :: b(nx,ny)
    integer i,j

    l=0
    do j=1,n2
       do i=1,n1
          a(i*2-1,j)=dreal(b(i,j))
          a(i*2  ,j)=dimag(b(i,j))
       enddo
    enddo

!    if (l>n3) write(*,*) "***** error in mset"

    return
  end subroutine mset2_h

#else
  subroutine eigen_solver_h(nn1,nn2)

!.M  use communication_h, only : eigen_init, eigen_free
    implicit none

    integer :: nx, ierr, lda, ldz, larray
    integer,intent(in) :: nn1, nn2
    integer :: nm_cache, nm_block
    real(8), allocatable :: utmp_r(:), utmp_i(:)
    real(8), allocatable :: vtmp_r(:), vtmp_i(:)
    integer :: npos, nb, nmz, nmw, nm, m0, nmz1, nmw1
    integer :: trilwmin, np, nq, lddz, lddw, numroc
    integer :: iam

    include 'commtxt.h'
                                                  __TIMER_SUB_START(907)
    if(iprisubmat>=2) write(nfout,*) "Submat: eigen_h_solver"

    npos = npcol*nprow

    nb  = 48
    nmz = ((meg-1)/nprow+1)
    nmz = ((nmz-1)/nb+1)*nb+1
    nmw = ((meg-1)/npcol+1)
    nmw = ((nmw-1)/nb+1)*nb+1

    nm = (meg/2)*2+1
    nmz1 = nm

    nmw1 = (meg-1)/npos+1
    m0 = 32

    lda = (meg-1)/npcol+1
    larray = MAX(nmz*nmw+nm+1,nmz1*nmw1)

    np = numroc( meg, nb, myrow, 0, nprow )
    nq = numroc( meg, nb, mycol, 0, npcol )
    lddz = (meg-1)/nprow+1
    lddz = ((lddz-1)/nb+1)*nb+1
    lddw = (meg-1)/npcol+1
    lddw = ((lddw-1)/nb+1)*nb+1

    trilwmin = 3*meg + max( nb*( np+1 ), 3*nb )
    ldz = max( max( 1+6*meg+2*np*nq, trilwmin ) + 2*meg, lddz*lddw)

    allocate(utmp_r(ldz), utmp_i(ldz))
    allocate(vtmp_r(ldz), vtmp_i(ldz))

    utmp_r = 0.d0
    utmp_i = 0.d0

    call mset_h(amat_l,utmp_r,utmp_i,nn1,nn2,larray)
    call eigen_h(meg,utmp_r(1),utmp_i(1),lda,vtmp_r(1),vtmp_i(1),larray,eig(1),m0,0)
    call mset2_h(zmat_l,utmp_r,utmp_i,nn1,nn2,larray)

    deallocate(utmp_r, utmp_i)
!   call MPI_comm_free(MPI_COMM_EIGEN,ierr)
                                                  __TIMER_SUB_STOP(907)
    return
  end subroutine eigen_solver_h

  subroutine mset_h(a,b,c,n1,n2,n3)
    implicit none
    real(8),intent(in)    :: a(n1*kimg_t,n2)
    real(8),intent(inout) :: b(n3), c(n3)
    integer n1,n2,n3
    integer i,j,l

    l=0
    do j=1,n2
       do i=1,n1
          l=l+1
          b(l)=a(i*2-1,j)
          c(l)=a(i*2  ,j)
       enddo
    enddo

    if (l>n3) write(*,*) "***** error in mset"

    return
  end subroutine mset_h

  subroutine mset2_h(a,b,c,n1,n2,n3)
    implicit none
    integer,intent(in) :: n1,n2,n3
    real(8),intent(inout) :: a(n1*kimg_t,n2)
    real(8),intent(in)    :: b(n3),c(n3)
    integer i,j,l

    l=0
    do j=1,n2
       do i=1,n1
          l=l+1
          a(i*2-1,j)=b(l)
          a(i*2  ,j)=c(l)
       enddo
    enddo

    if (l>n3) write(*,*) "***** error in mset"

    return
  end subroutine mset2_h
#endif
#endif
!endif if USE_EIGENLIB

  subroutine eigsend()
    integer :: i, ierr, itag=100
    integer, dimension(MPI_STATUS_SIZE) :: stat
                                                  __TIMER_SUB_START(988)
                                                  __TIMER_COMM_START(989)
    if ((myrank_e==0) .and. (myrank_g==0)) then
#ifdef EIGEN_6D
       do i=0,nrank_e*nrank_g*nrank_k-1
       if ( icol(i) == 1 ) cycle
#else
       do i = nrank_e*nrank_g-1, nsclcol*nsclrow, -1
#endif
          call mpi_send(eig, meg, mpi_double_precision, i, &
         &              itag, mpi_k_world(myrank_k), ierr)
           if (ierr /= 0) then
              write(nfout,*)' eigsend :  mpi_send error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 10010, ierr)
           end if
       end do
    end if
                                                  __TIMER_COMM_STOP(989)
                                                  __TIMER_COMM_START(990)
#ifdef EIGEN_6D
    if (icolor==0) then
#else
#ifdef ASSIGN_G_PREVIOUS
    if ((nsclcol*nsclrow-1) < (nrank_g*myrank_e+myrank_g)) then
#else
    if ((nsclcol*nsclrow-1) < (nrank_e*myrank_g+myrank_e)) then
#endif
#endif
       call mpi_recv(eig, meg, mpi_double_precision, 0, &
      &              itag, mpi_k_world(myrank_k), stat, ierr)
        if (ierr /= 0) then
           write(nfout,*)' eigsend :  mpi_irecv error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 10012, ierr)
        end if
    end if
                                                  __TIMER_COMM_STOP(990)
                                                  __TIMER_SUB_STOP(988)
    end subroutine eigsend

    subroutine gather_zmat_all(zmat_l, zmat)
    integer(kind=4) :: i, j

    !!real(kind=DP), allocatable :: zmat_l(:,:), zmat(:,:)
#ifdef _USE_SCALAPACK_
! === DEBUG by Tkato 2011/06/28 ================================================
!   real(kind=DP), intent(out) :: zmat_l(maxeg*kimg_t,maxe)
    real(kind=DP), intent(in) :: zmat_l(lda*kimg_t,occ)
! ==============================================================================
#else
    real(kind=DP), intent(in) :: zmat_l(nel_eg(myrank_g)*kimg_t,np_e)
#endif
    real(kind=DP), intent(out) :: zmat(meg*kimg_t,meg)
    real(kind=DP), allocatable :: zmat_wk(:,:),zmat_all(:,:)
                                                  __TIMER_SUB_START(905)
    allocate(wk_mpi(maxeg*kimg_t,maxe), stat=ierr)
     if (ierr /= 0) then
        call mpi_abort(mpi_comm_world, 117, ierr)
     end if
    allocate(wk_gather(maxeg*kimg_t,maxe,0:nrank_e*nrank_g-1), stat=ierr)
     if (ierr /= 0) then
        call mpi_abort(mpi_comm_world, 118, ierr)
     end if
!   wk_mpi = 0.0d0
!   wk_gather = 0.0d0
                                                  __TIMER_DO_START(948)
    if(kimg_t==1) then
       do j = 1, np_e
          do i = 1, nel_eg(myrank_g)
             wk_mpi(i,j) = zmat_l(i,j)
          enddo
       enddo
    else
       do j = 1, np_e
          do i = 1, nel_eg(myrank_g)
             wk_mpi(i*2-1,j) = zmat_l(i*2-1,j)
             wk_mpi(i*2  ,j) = zmat_l(i*2  ,j)
          enddo
       enddo
    endif
                                                  __TIMER_DO_STOP(948)
                                                  __TIMER_COMM_START_w_BARRIER(mpi_k_world(myrank_k),949)
    call mpi_allgather(wk_mpi,    maxeg*kimg_t*maxe, mpi_double_precision,&
   &                   wk_gather, maxeg*kimg_t*maxe, mpi_double_precision,&
   &                   mpi_k_world(myrank_k), ierr)
                                                  __TIMER_COMM_STOP(949)
     if (ierr /= 0) then
        call mpi_abort(mpi_comm_world, 119, ierr)
     end if

    allocate(zmat_all(neg*kimg_t,neg), stat=ierr)
     if (ierr /= 0) then
        call mpi_abort(mpi_comm_world, 120, ierr)
     end if

                                                  __TIMER_DO_START(950)
    if(kimg_t==1) then
       do l = 0, nrank_g-1
          do k = 0, nrank_e-1
             do j = 1, nel_e(k)
                do i = 1, nel_eg(l)
! === DEBUG by tkato 2013/08/28 ================================================
!                  zmat_all(i+nis_eg(l)-1,j+nis_e(k)-1) = wk_gather(i,j,l*nrank_e+k)
#ifndef ASSIGN_G_PREVIOUS
                   zmat_all(i+nis_eg(l)-1,j+nis_e(k)-1) = wk_gather(i,j,l*nrank_e+k)
#else
                   zmat_all(i+nis_eg(l)-1,j+nis_e(k)-1) = wk_gather(i,j,k*nrank_g+l)
#endif
! ==============================================================================
                enddo
             enddo
          enddo
       enddo
    else
       do l = 0, nrank_g-1
          do k = 0, nrank_e-1
             do j = 1, nel_e(k)
                do i = 1, nel_eg(l)
! === DEBUG by tkato 2013/08/28 ================================================
!                  zmat_all((i+nis_eg(l)-1)*2-1,j+nis_e(k)-1) = wk_gather(i*2-1,j,l*nrank_e+k)
!                  zmat_all((i+nis_eg(l)-1)*2  ,j+nis_e(k)-1) = wk_gather(i*2  ,j,l*nrank_e+k)
#ifndef ASSIGN_G_PREVIOUS
                   zmat_all((i+nis_eg(l)-1)*2-1,j+nis_e(k)-1) = wk_gather(i*2-1,j,l*nrank_e+k)
                   zmat_all((i+nis_eg(l)-1)*2  ,j+nis_e(k)-1) = wk_gather(i*2  ,j,l*nrank_e+k)
#else
                   zmat_all((i+nis_eg(l)-1)*2-1,j+nis_e(k)-1) = wk_gather(i*2-1,j,k*nrank_g+l)
                   zmat_all((i+nis_eg(l)-1)*2  ,j+nis_e(k)-1) = wk_gather(i*2  ,j,k*nrank_g+l)
#endif
! ==============================================================================
                enddo
             enddo
          enddo
       enddo
    endif
                                                  __TIMER_DO_STOP(950)
    deallocate(wk_gather)
    deallocate(wk_mpi)

    allocate(zmat_wk(neg*kimg_t,neg))
     if (ierr /= 0) then
        call mpi_abort(mpi_comm_world, 121, ierr)
     end if
    zmat_wk = 0.0d0
                                                  __TIMER_DO_START(951)
    if(kimg_t==1) then
!OCL NOFLTLD
       do j = 1,neg
          if (neg_g_all(j) > meg) cycle
          do i = 1,neg
             if (neg_gg_all(i) > meg) cycle
             zmat_wk(neg_gg_all(i),neg_g_all(j)) = zmat_all(i,j)
          enddo
       enddo
    else
!OCL NOFLTLD
       do j = 1,neg
          if (neg_g_all(j) > meg) cycle
          do i = 1,neg
             if (neg_gg_all(i) > meg) cycle
             zmat_wk(neg_gg_all(i)*2-1,neg_g_all(j)) = zmat_all(i*2-1,j)
             zmat_wk(neg_gg_all(i)*2  ,neg_g_all(j)) = zmat_all(i*2  ,j)
          enddo
       enddo
    endif
                                                  __TIMER_DO_STOP(951)
    deallocate(zmat_all)
                                                  __TIMER_DO_START(952)
    zmat=0.d0
    if(kimg_t==1) then
       do j= 1, meg
          do i = 1, j
             zmat(i,j) = zmat_wk(i,j)
          enddo
       enddo
    else
       do j= 1, meg
          do i = 1, j
             zmat(i*2-1,j) = zmat_wk(i*2-1,j)
             zmat(i*2  ,j) = zmat_wk(i*2  ,j)
          enddo
       enddo
    endif
                                                  __TIMER_DO_STOP(952)
    deallocate(zmat_wk)
                                                  __TIMER_SUB_STOP(905)
    end subroutine gather_zmat_all

#ifdef EIGEN_TRANS
  subroutine make_usermap(usermap,icol,ikey,irank_c,irank_r)
#else
  subroutine make_usermap(usermap,icol,ikey)
#endif

    implicit none
    integer, intent(out) :: usermap(nprow,npcol)
    integer, intent(out) :: icol(0:nrank_s*nrank_e*nrank_g*nrank_k-1)
    integer, intent(out) :: ikey(0:nrank_s*nrank_e*nrank_g*nrank_k-1)
#ifdef EIGEN_TRANS
    integer, intent(out) :: irank_c(0:nrank_s*nrank_e*nrank_g*nrank_k-1)
    integer, intent(out) :: irank_r(0:nrank_s*nrank_e*nrank_g*nrank_k-1)
#endif

    integer :: ix, iy, iz, ia, ib, ic, irank, imode, icont
    integer :: nx, ny, nz, mx, my, mz, nb
    logical :: EIGEN_6D_
    integer :: id_sname = -1
                                                  __TIMER_SUB_START(942)
    call tstatc0_begin('make_usermap', id_sname)

    imode = 2
    icont = 0
    icol(:) = 0
    do i=0,nrank_e*nrank_g*nrank_k-1
      ikey(i)= nrank_e*nrank_g*nrank_k+i
    enddo
#ifdef EIGEN_TRANS
    irank_c(:) = -1
    irank_r(:) = -1
#endif

#ifdef EIGEN_6D
    EIGEN_6D_=.true.
! === DEBUG by tkato 2013/09/18 ================================================
!#else
!    EIGEN_6D_=.false.
!#endif
! ==============================================================================
    if ( EIGEN_6D_ ) then
      imode=0
      if ( nprow == 32 .and. npcol == 32 ) then
        nb=1; nx=7; ny=3; nz=3
        mx=7; my=7; mz=15
      endif
      if ( nprow == 16 .and. npcol == 16 ) then
        nb=1; nx=3; ny=1; nz=3
        mx=7; my=3; mz=7
      endif
      if ( nprow == 8 .and. npcol == 8 ) then
        nb=1; nx=1; ny=1; nz=1
        mx=3; my=3; mz=3
      endif
      if ( nprow == 4 .and. npcol == 4 ) then
        nb=0; nx=1; ny=0; nz=1
        mx=3; my=1; mz=1
      endif
      if ( nprow == 2 .and. npcol == 2 ) then
        nb=0; nx=0; ny=0; nz=0
        mx=1; my=0; mz=1
      endif

      do ib=0,nb
        do iy=0,ny
          do iz=0,nz
            do ic=0,1
              do ix=0,nx
                do ia=0,1
                  j=icont/npcol
                  i=mod(icont,npcol)
                  call get6d2rank(ix,iy,iz,ia,ib,ic,irank)
                  if (irank < 0) imode=1
                  usermap(i+1,j+1)=myrank_k*nrank_e*nrank_g+irank
#ifdef EIGEN_TRANS
                  irank_c(irank)=j
                  irank_r(irank)=i
#endif
                  icol(irank)=1
                  ikey(irank)=icont
                  icont=icont+1
                end do
              end do
            end do
          end do
        end do
      end do
      if ( imode == 1 ) then
        icont = 0
        do iz=0,mz
          do iy=0,my
            do ix=0,mx
              j=icont/npcol
              i=mod(icont,npcol)
              call get3d2rank(ix,iy,iz,irank)
              if (irank < 0) imode=2
              usermap(i+1,j+1)=myrank_k*nrank_e*nrank_g+irank
#ifdef EIGEN_TRANS
              irank_c(irank)=j
              irank_r(irank)=i
#endif
              icol(irank)=1
              ikey(irank)=icont
              icont=icont+1
            end do
          end do
        end do
      endif
    endif
! === DEBUG by tkato 2013/09/18 ================================================
#else
    EIGEN_6D_=.false.
#endif
! ==============================================================================

    if ( imode == 2 ) then
    do j=1,npcol
       do i=1,nprow
!fj --------------------
!fj       usermap(i,j) = myrank_k*nrank_e + (j-1)*nprow + i - 1
          if(m_CtrlP_get_conftag()>=0) then
#if defined(_ASSIGN_ROW_) || defined(USE_EIGENLIB)
            usermap(i,j) = m_CtrlP_get_conftag()*nrank_k*nrank_e*nrank_g + myrank_k*nrank_e*nrank_g &
                         + myrank_spin*nrank_k*nrank_e*nrank_g + (j-1)*nprow + i - 1
#else
            usermap(i,j) = m_CtrlP_get_conftag()*nrank_k*nrank_e*nrank_g + myrank_k*nrank_e*nrank_g &
                         + myrank_spin*nrank_k*nrank_e*nrank_g + (i-1)*npcol + j - 1
#endif
          else if(nrank_conf>1) then
#if defined(_ASSIGN_ROW_) || defined(USE_EIGENLIB)
            usermap(i,j) = mype_conf*nrank_k*nrank_e*nrank_g*nrank_s + myrank_k*nrank_e*nrank_g &
                         + myrank_spin*nrank_k*nrank_e*nrank_g + (j-1)*nprow + i - 1
#else
            usermap(i,j) = mype_conf*nrank_k*nrank_e*nrank_g*nrank_s + myrank_k*nrank_e*nrank_g &
                         + myrank_spin*nrank_k*nrank_e*nrank_g + (i-1)*npcol + j - 1
#endif
          else
#if defined(_ASSIGN_ROW_) || defined(USE_EIGENLIB)
            usermap(i,j) = myrank_k*nrank_e*nrank_g + myrank_spin*nrank_k*nrank_e*nrank_g + (j-1)*nprow + i - 1
#else
            usermap(i,j) = myrank_k*nrank_e*nrank_g + myrank_spin*nrank_k*nrank_e*nrank_g + (i-1)*npcol + j - 1
#endif
          endif
!fj --------------------
       end do
    end do
    endif

    if(iprisubmat>=2) then
       write(nfout,'("USERMAP mode=",i1)') imode
!      do i=1,nprow
!         write(nfout,'(16(1x,i5))') (usermap(i,j),j=1,npcol)
!      end do
       do j=1,npcol
          write(nfout,'(16(1x,i5))') (usermap(i,j),i=1,nprow)
       end do
       call flush(nfout)
    end if

    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(942)
  end subroutine make_usermap

!!#ifdef _USE_SCALAPACK_
#ifdef EIGEN_TRANS
#ifdef EIGEN_6D
    subroutine trans_scalapack(neg, meg, mgs_nb,scl_nb,zmat_l,amat_l,maxeg,maxe,lda,occ,kimg_t, usermap, irank_c, irank_r)
#else
    subroutine trans_scalapack(neg, meg, mgs_nb,scl_nb,zmat_l,amat_l,maxeg,maxe,lda,occ,kimg_t)
#endif

     integer, intent(in) :: neg, meg, mgs_nb, scl_nb, maxeg, maxe, lda, occ, kimg_t
     real(kind=DP), intent(in) :: zmat_l(maxeg*kimg_t,maxe)
     real(kind=DP), intent(inout) :: amat_l(lda*kimg_t,occ)

     integer, allocatable, dimension(:) :: wk_rk
     integer :: i, k, j,  nn, mm
     integer :: lrk, lad
     integer :: lno_r, ng_r, lrk_r, rank_r
     integer :: lno_c, ng_c, lrk_c, rank_c

     real(kind=DP), allocatable, dimension(:,:) :: send_buf, recv_buf

     integer, allocatable, dimension(:)   :: req_r, req_s
     integer, allocatable, dimension(:,:) :: sta_r, sta_s
     integer :: lrank, itag, icnt_send, icnt_recv, mpi_comm, ierr
#ifdef EIGEN_6D
     integer, intent(in) :: irank_c(0:nrank_e*nrank_g*nrank_k-1)
     integer, intent(in) :: irank_r(0:nrank_e*nrank_g*nrank_k-1)
     integer, intent(in) :: usermap(scl_row,scl_col)
#endif
! === DEBUG by tkato 2012/06/04 ================================================
#ifndef USE_NONBLK_COMM
#ifndef USE_ALLTOALLV
       real(kind=DP), allocatable, dimension(:,:) :: sendb, recvb
#else
       real(kind=DP), allocatable, dimension(:,:) :: sbuf, rbuf
       integer, allocatable, dimension(:) :: scnt, rcnt
       integer, allocatable, dimension(:) :: sdsp, rdsp
#endif
#endif
! ==============================================================================

     if (iprisubmat>=2) then
        write(nfout,'("lda,occ=",2i5)') lda,occ
        write(nfout,'("==== trans_scalapack ====")')
        write(nfout,'(" neg=",i4,", meg=",i4,", mgs_nb=",i4,", scl_nb=",i4,", maxeg=",i4,", maxe=",i4,", &
           & lda=",i4,", occ=",i4,", kimg_t=",i4)') neg, meg, mgs_nb, scl_nb, maxeg, maxe, lda, occ, kimg_t
        call flush(nfout)
      end if

      mpi_comm = mpi_k_world(myrank_k)
      itag = 10001
      allocate(req_r(scl_comm_rank_r), stat=ierr)
      allocate(req_s(scl_comm_rank  ), stat=ierr)
      allocate(sta_r(MPI_STATUS_SIZE,scl_comm_rank_r), stat=ierr)
      allocate(sta_s(MPI_STATUS_SIZE,scl_comm_rank  ), stat=ierr)
      allocate(recv_buf(scl_comm_max_r*kimg_t,scl_comm_rank_r), stat=ierr)
      allocate(send_buf(scl_comm_max  *kimg_t,scl_comm_rank  ), stat=ierr)
      icnt_recv = 0
      icnt_send = 0

#ifdef USE_NONBLK_COMM
                                                  __TIMER_DO_START(971)
#endif
      do i = 1, scl_comm_rank_r
         lrank = scl_comm_rno_r(i)
         if (lrank /= mype) then
            icnt_recv = icnt_recv + 1
#ifdef USE_NONBLK_COMM
            call mpi_irecv(recv_buf(1,i), scl_comm_max_r*kimg_t, &
           &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
             if (ierr /= 0) then
                write(nfout,*)' trans_scalapack :  mpi_irecv error'
                call flush(nfout)
                call mpi_abort(mpi_comm_world, 10001, ierr)
             endif
#endif
         endif
      enddo
#ifdef USE_NONBLK_COMM
                                                  __TIMER_DO_STOP(971)
#endif

      allocate(wk_rk(max(scl_comm_rank,scl_comm_rank_r)))
      wk_rk = 0
      rank_c = myrank_e
      rank_r = myrank_g
                                                  __TIMER_DO_START(972)
      do k = nis_e(rank_c), nie_e(rank_c)
         ng_c = neg_g_all(k)
         if (ng_c > meg) cycle
         if (mod(ng_c,scl_nb) > 0) then
            nn = ng_c / scl_nb + 1
         else
            nn = ng_c / scl_nb
         end if
         lrk_c = mod((nn-1),scl_col)
         do j = nis_eg(rank_r), nie_eg(rank_r)
            ng_r = neg_gg_all(j)
            if (ng_r > meg) cycle
            if (mod(ng_r,scl_nb) > 0) then
               nn = ng_r / scl_nb + 1
            else
               nn = ng_r / scl_nb
            end if
            lrk_r = mod((nn-1),scl_row)
#ifdef EIGEN_6D
            lrk = usermap(lrk_r+1,lrk_c+1)
#else
#if defined(_ASSIGN_ROW_) || defined(USE_EIGENLIB)
            lrk = scl_row*lrk_c+lrk_r
#else
            lrk = scl_col*lrk_r+lrk_c
#endif
#endif
!!          if (lrk /= mype) then
               do lrank = 1, scl_comm_rank
                  if (scl_comm_rno(lrank) == lrk) then
                     wk_rk(lrank) = wk_rk(lrank) + 1
                     if (kimg_t == 1) then
                        send_buf(wk_rk(lrank),lrank) = zmat_l(j-nis_eg(rank_r)+1,k-nis_e(rank_c)+1)
                     else
                        send_buf(wk_rk(lrank)*2-1,lrank) = zmat_l((j-nis_eg(rank_r)+1)*2-1,k-nis_e(rank_c)+1)
                        send_buf(wk_rk(lrank)*2  ,lrank) = zmat_l((j-nis_eg(rank_r)+1)*2  ,k-nis_e(rank_c)+1)
                     end if
                     exit
                  end if
               end do
!!          end if
         end do
      end do
                                                  __TIMER_DO_STOP(972)
                                                  __TIMER_COMM_START_w_BARRIER(mpi_com,975)
#ifdef USE_NONBLK_COMM
                                                  __TIMER_DO_START(971)
      do i = 1, scl_comm_rank
         lrank = scl_comm_rno(i)
         if (lrank /= mype) then
            icnt_send = icnt_send + 1
            call mpi_isend(send_buf(1,i), scl_comm_max*kimg_t, &
           &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
             if (ierr /= 0) then
                write(nfout,*)' trans_scalapack :  mpi_isend error'
                call flush(nfout)
                call mpi_abort(mpi_comm_world, 10002, ierr)
             endif
         endif
      enddo
                                                  __TIMER_DO_STOP(971)
#else
      do i = 1, scl_comm_rank
         lrank = scl_comm_rno(i)
         if (lrank /= mype) icnt_send = icnt_send + 1
      enddo
#endif

#ifdef EIGEN_6D
      rank_c = irank_c(mype-nrank_e*nrank_g*myrank_k)
      rank_r = irank_r(mype-nrank_e*nrank_g*myrank_k)
   if (rank_c > -1 .and. rank_r > -1 ) then
#else
   if ((mype-nrank_e*nrank_g*myrank_k) < (scl_col*scl_row)) then
#if defined(_ASSIGN_ROW_) || defined(USE_EIGENLIB)
      rank_c = (mype-nrank_e*nrank_g*myrank_k)/scl_row
      rank_r = mod((mype-nrank_e*nrank_g*myrank_k),scl_row)
#else
      rank_c = mod((mype-nrank_e*nrank_g*myrank_k),scl_col)
      rank_r = (mype-nrank_e*nrank_g*myrank_k)/scl_col
#endif
#endif
                                                  __TIMER_DO_START(973)
      do k = nis_col(rank_c), nie_col(rank_c)
         ng_c = neg_col(k)
         if( ng_c > meg ) cycle
         if (mod(ng_c,mgs_nb) > 0) then
            nn = ng_c / mgs_nb + 1
         else
            nn = ng_c / mgs_nb
         end if
         lrk_c = mod((nn-1),nrank_e)
         lad = mod(ng_c,mgs_nb)
         if (lad == 0) lad = mgs_nb
         mm  = ((ng_c-1)/mgs_nb)/nrank_e
         lno_c = mm * mgs_nb + lad

         do j = nis_row(rank_r), nie_row(rank_r)
            ng_r = neg_row(j)
            if( ng_r > meg ) cycle
            if (mod(ng_r,mgs_nb) > 0) then
               nn = ng_r / mgs_nb + 1
            else
               nn = ng_r / mgs_nb
            end if
            lrk_r = mod((nn-1),nrank_g)
            lad = mod(ng_r,mgs_nb)
            if (lad == 0) lad = mgs_nb
            mm  = ((ng_r-1)/mgs_nb)/nrank_g
            lno_r = mm * mgs_nb + lad
#ifdef ASSIGN_G_PREVIOUS
            lrk = nrank_g*lrk_c+lrk_r
#else
            lrk = nrank_e*lrk_r+lrk_c
#endif
            if (lrk == mype) then

               if(ng_r> ng_c) cycle

               if (kimg_t == 1) then
                  amat_l(j-nis_row(rank_r)+1,k-nis_col(rank_c)+1) = zmat_l(lno_r,lno_c)
               else
                  amat_l((j-nis_row(rank_r)+1)*2-1,k-nis_col(rank_c)+1) = zmat_l(lno_r*2-1,lno_c)
                  amat_l((j-nis_row(rank_r)+1)*2  ,k-nis_col(rank_c)+1) = zmat_l(lno_r*2  ,lno_c)
               end if
            end if
         end do
      end do
                                                  __TIMER_DO_STOP(973)
   end if

#ifdef USE_NONBLK_COMM
      call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
       if (ierr /= 0) then
          write(nfout,*)' trans_scalapack :  mpi_waitall error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 10003, ierr)
       endif

      call mpi_waitall(icnt_send, req_s, sta_s, ierr)
       if (ierr /= 0) then
          write(nfout,*)' trans_scalapack :  mpi_waitall error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 10004, ierr)
       endif
                                                  __TIMER_COMM_STOP(975)
#else
                                                  __TIMER_COMM_START_w_BARRIER(mpi_k_world,992)
#ifndef USE_ALLTOALLV
! === DEBUG by tkato 2012/06/04 ================================================
!      real(kind=DP), allocatable, dimension(:,:) :: sendb, recvb
! ==============================================================================
       allocate(sendb(scl_comm_max  *kimg_t,0:nrank_e*nrank_g-1), stat=ierr)
       allocate(recvb(scl_comm_max_r*kimg_t,0:nrank_e*nrank_g-1), stat=ierr)
       sendb=0
       recvb=0
       do i = 1, scl_comm_rank
          sendb(:,scl_comm_rno(i))=send_buf(:,i)
       enddo
       call mpi_alltoall(sendb, scl_comm_max  *kimg_t, mpi_double_precision, &
      &                  recvb, scl_comm_max_r*kimg_t, mpi_double_precision, &
      &                                                MPI_COMM_WORLD, ierr)
       if (ierr /= 0) then
          write(nfout,*)' trans_scalapack :  mpi_alltoall(List) error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 10005, ierr)
       endif
       do i = 1, scl_comm_rank_r
          recv_buf(:,i)=recvb(:,scl_comm_rno_r(i))
       enddo
       deallocate(sendb)
       deallocate(recvb)
#else
! === DEBUG by tkato 2012/06/05 ================================================
!      real(kind=DP), allocatable, dimension(:,:) :: sbuf, rbuf
!      integer, allocatable, dimension(:) :: scnt, rcnt
!      integer, allocatable, dimension(:) :: sdsp, rdsp
! ==============================================================================
       allocate(sbuf(scl_comm_max  *kimg_t,0:nrank_e*nrank_g-1), stat=ierr)
       allocate(rbuf(scl_comm_max_r*kimg_t,0:nrank_e*nrank_g-1), stat=ierr)
       allocate(scnt(0:nrank_e*nrank_g-1), stat=ierr)
       allocate(rcnt(0:nrank_e*nrank_g-1), stat=ierr)
       allocate(sdsp(0:nrank_e*nrank_g-1), stat=ierr)
       allocate(rdsp(0:nrank_e*nrank_g-1), stat=ierr)
       sbuf=0
       rbuf=0
       scnt=0
       rcnt=0
       do i = 1, scl_comm_rank
          lrank = scl_comm_rno(i)
          sbuf(:,lrank)=send_buf(:,i)
          scnt(lrank)=scl_comm_max  *kimg_t
       enddo
       do i = 1, scl_comm_rank_r
          rcnt(scl_comm_rno_r(i))=scl_comm_max_r*kimg_t
       enddo
       do i = 0, nrank_e*nrank_g - 1
          sdsp(i)=scl_comm_max  *kimg_t*i
          rdsp(i)=scl_comm_max_r*kimg_t*i
       enddo
       call mpi_alltoallv(sbuf, scnt, sdsp, mpi_double_precision, &
      &                   rbuf, rcnt, rdsp, mpi_double_precision, &
      &                                            MPI_COMM_WORLD, ierr)
       if (ierr /= 0) then
          write(nfout,*)' trans_scalapack :  mpi_alltoallv error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 10005, ierr)
       endif
       do i = 1, scl_comm_rank_r
          recv_buf(:,i)=rbuf(:,scl_comm_rno_r(i))
       enddo
       deallocate(sbuf)
       deallocate(rbuf)
       deallocate(scnt)
       deallocate(rcnt)
       deallocate(sdsp)
       deallocate(rdsp)
#endif
                                                  __TIMER_COMM_STOP(992)
#endif

#ifdef EIGEN_6D
      rank_c = irank_c(mype-nrank_e*nrank_g*myrank_k)
      rank_r = irank_r(mype-nrank_e*nrank_g*myrank_k)
      if (rank_c > -1 .and. rank_r > -1 ) then
         wk_rk = 0
#else
      if ((mype-nrank_e*nrank_g*myrank_k) < (scl_col*scl_row)) then

         wk_rk = 0
#if defined(_ASSIGN_ROW_) || defined(USE_EIGENLIB)
         rank_c = (mype-nrank_e*nrank_g*myrank_k)/scl_row
         rank_r = mod((mype-nrank_e*nrank_g*myrank_k),scl_row)
#else
         rank_c = mod((mype-nrank_e*nrank_g*myrank_k),scl_col)
         rank_r = (mype-nrank_e*nrank_g*myrank_k)/scl_col
#endif
#endif
                                                  __TIMER_DO_START(974)
         do k = nis_col(rank_c), nie_col(rank_c)
            ng_c = neg_col(k)
            if( ng_c > meg ) cycle
            if (mod(ng_c,mgs_nb) > 0) then
               nn = ng_c / mgs_nb + 1
            else
               nn = ng_c / mgs_nb
            end if
            lrk_c = mod((nn-1),nrank_e)
            do j = nis_row(rank_r), nie_row(rank_r)
               ng_r = neg_row(j)
               if( ng_r > meg ) cycle
               if (mod(ng_r,mgs_nb) > 0) then
                  nn = ng_r / mgs_nb + 1
               else
                  nn = ng_r / mgs_nb
               end if
               lrk_r = mod((nn-1),nrank_g)
#ifdef ASSIGN_G_PREVIOUS
               lrk = nrank_g*lrk_c+lrk_r
#else
               lrk = nrank_e*lrk_r+lrk_c
#endif
!!             if (lrk /= mype) then
                  do lrank = 1, scl_comm_rank_r
                     if (scl_comm_rno_r(lrank) == lrk) then
                        wk_rk(lrank) = wk_rk(lrank) + 1

                        if(ng_r> ng_c) cycle

                        if (kimg_t == 1) then
                           amat_l(j-nis_row(rank_r)+1,k-nis_col(rank_c)+1) = recv_buf(wk_rk(lrank),lrank)
                        else
                           amat_l((j-nis_row(rank_r)+1)*2-1,k-nis_col(rank_c)+1) = recv_buf(wk_rk(lrank)*2-1,lrank)
                           amat_l((j-nis_row(rank_r)+1)*2  ,k-nis_col(rank_c)+1) = recv_buf(wk_rk(lrank)*2  ,lrank)
                        end if
                        exit
                     end if
                  end do
!!             end if
            end do
         end do
                                                  __TIMER_DO_STOP(974)
      end if

      deallocate(req_r, stat=ierr)
      deallocate(req_s, stat=ierr)
      deallocate(sta_r, stat=ierr)
      deallocate(sta_s, stat=ierr)
      deallocate(send_buf, stat=ierr)
      deallocate(recv_buf, stat=ierr)
      deallocate(wk_rk, stat=ierr)

    end subroutine trans_scalapack

#ifdef EIGEN_6D
    subroutine trans_scalapack_r(neg, meg, mgs_nb,scl_nb,zmat_l,amat_l,maxeg,maxe,lda,occ,kimg_t, usermap, irank_c, irank_r)
#else
    subroutine trans_scalapack_r(neg, meg, mgs_nb,scl_nb,zmat_l,amat_l,maxeg,maxe,lda,occ,kimg_t)
#endif

     integer, intent(in) :: neg, meg, mgs_nb, scl_nb, maxeg, maxe, lda, occ, kimg_t
     real(kind=DP), intent(in) :: zmat_l(lda*kimg_t,occ)
     real(kind=DP), intent(inout) :: amat_l(maxeg*kimg_t,maxe)

     integer , allocatable, dimension(:) :: wk_rk
     integer :: i, k, j, nn, mm
     integer :: lrk, lad
     integer :: lno_r, ng_r, lrk_r, rank_r
     integer :: lno_c, ng_c, lrk_c, rank_c

     real(kind=DP), allocatable, dimension(:,:) :: send_buf, recv_buf

     integer, allocatable, dimension(:)   :: req_r, req_s
     integer, allocatable, dimension(:,:) :: sta_r, sta_s
     integer :: lrank, itag, icnt_send, icnt_recv, mpi_comm, ierr
#ifdef EIGEN_6D
     integer, intent(in) :: usermap(scl_row,scl_col)
     integer, intent(in) :: irank_c(0:nrank_e*nrank_g*nrank_k-1)
     integer, intent(in) :: irank_r(0:nrank_e*nrank_g*nrank_k-1)
#endif

! === DEBUG by tkato 2012/06/04 ================================================
#ifndef USE_NONBLK_COMM
#ifndef USE_ALLTOALLV
     real(kind=DP), allocatable, dimension(:,:) :: sbuf, rbuf
#else
     real(kind=DP), allocatable, dimension(:,:) :: sbuf, rbuf
     integer, allocatable, dimension(:) :: scnt, rcnt
     integer, allocatable, dimension(:) :: sdsp, rdsp
#endif
#endif
! ==============================================================================

      mpi_comm = mpi_k_world(myrank_k)
      itag = 10002
      if (iprisubmat>=2) then
         write(nfout,'("==== trans_scalapack_r ====")')
         write(nfout,'(" neg=",i4,", meg=",i4,", mgs_nb=",i4,", scl_nb=",i4,", maxeg=",i4,", maxe=",i4,", &
           & lda=",i4,", occ=",i4,", kimg_t=",i4)') neg, meg, mgs_nb, scl_nb, maxeg, maxe, lda, occ, kimg_t
         write(nfout,'(" scl_comm_rank=",i4,", scl_comm_rank_r=",i4)') scl_comm_rank, scl_comm_rank_r
         call flush(nfout)
      end if

      allocate(req_r(scl_comm_rank  ), stat=ierr)
      allocate(req_s(scl_comm_rank_r), stat=ierr)
      allocate(sta_r(MPI_STATUS_SIZE,scl_comm_rank  ), stat=ierr)
      allocate(sta_s(MPI_STATUS_SIZE,scl_comm_rank_r), stat=ierr)
      allocate(recv_buf(scl_comm_max  *kimg_t,scl_comm_rank  ), stat=ierr)
      allocate(send_buf(scl_comm_max_r*kimg_t,scl_comm_rank_r), stat=ierr)
      icnt_recv = 0
      icnt_send = 0
!     send_buf=0.0d0
!     recv_buf=0.0d0

#ifdef USE_NONBLK_COMM
                                                  __TIMER_DO_START(976)
      do i = 1, scl_comm_rank
         lrank = scl_comm_rno(i)
         if (lrank /= mype) then
            icnt_recv = icnt_recv + 1
            call mpi_irecv(recv_buf(1,i), scl_comm_max*kimg_t, &
           &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
             if (ierr /= 0) then
                write(nfout,*)' trans_scalapack_r :  mpi_irecv error'
                call flush(nfout)
                call mpi_abort(mpi_comm_world, 10001, ierr)
             endif
         endif
      enddo
                                                  __TIMER_DO_STOP(976)
#else
      do i = 1, scl_comm_rank
         lrank = scl_comm_rno(i)
         if (lrank /= mype) icnt_recv = icnt_recv + 1
      enddo
#endif

      allocate(wk_rk(max(scl_comm_rank,scl_comm_rank_r)))

#ifdef EIGEN_6D
      rank_c = irank_c(mype-nrank_e*nrank_g*myrank_k)
      rank_r = irank_r(mype-nrank_e*nrank_g*myrank_k)
      if (rank_c > -1 .and. rank_r > -1 ) then
         wk_rk = 0
#else
      if ((mype-nrank_e*nrank_g*myrank_k) < (scl_col*scl_row)) then

         wk_rk = 0

#if defined(_ASSIGN_ROW_) || defined(USE_EIGENLIB)
         rank_c = (mype-nrank_e*nrank_g*myrank_k)/scl_row
         rank_r = mod((mype-nrank_e*nrank_g*myrank_k),scl_row)
#else
         rank_c = mod((mype-nrank_e*nrank_g*myrank_k),scl_col)
         rank_r = (mype-nrank_e*nrank_g*myrank_k)/scl_col
#endif
#endif
                                                  __TIMER_DO_START(977)
         do k = nis_col(rank_c), nie_col(rank_c)
            ng_c = neg_col(k)
            if (ng_c > meg) cycle
            if (mod(ng_c,mgs_nb) > 0) then
               nn = ng_c / mgs_nb + 1
            else
               nn = ng_c / mgs_nb
            end if
            lrk_c = mod((nn-1),nrank_e)
            do j = nis_row(rank_r), nie_row(rank_r)
               ng_r = neg_row(j)
               if (ng_r > meg) cycle
               if (mod(ng_r,mgs_nb) > 0) then
                  nn = ng_r / mgs_nb + 1
               else
                  nn = ng_r / mgs_nb
               end if
               lrk_r = mod((nn-1),nrank_g)
#ifdef ASSIGN_G_PREVIOUS
               lrk = nrank_g*lrk_c+lrk_r
#else
               lrk = nrank_e*lrk_r+lrk_c
#endif
!!             if (lrk /= mype) then
                  do lrank = 1, scl_comm_rank_r
                     if(scl_comm_rno_r(lrank) == lrk) then
                        wk_rk(lrank) = wk_rk(lrank) + 1
                        if (kimg_t == 1) then
                           send_buf(wk_rk(lrank),lrank) = zmat_l(j-nis_row(rank_r)+1,k-nis_col(rank_c)+1)
                        else
                           send_buf(wk_rk(lrank)*2-1,lrank) = zmat_l((j-nis_row(rank_r)+1)*2-1,k-nis_col(rank_c)+1)
                           send_buf(wk_rk(lrank)*2  ,lrank) = zmat_l((j-nis_row(rank_r)+1)*2  ,k-nis_col(rank_c)+1)
                        end if
                        exit
                     end if
                  end do
!!             end if
            end do
         end do
                                                  __TIMER_DO_STOP(977)
      end if

#ifdef USE_NONBLK_COMM
                                                  __TIMER_COMM_START_w_BARRIER(mpi_comm,980)
                                                  __TIMER_DO_START(976)
      do i = 1, scl_comm_rank_r
         lrank = scl_comm_rno_r(i)
         if (lrank /= mype) then
            icnt_send = icnt_send + 1
            call mpi_isend(send_buf(1,i), scl_comm_max*kimg_t, &
           &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
             if (ierr /= 0) then
                write(nfout,*)' trans_scalapack_r :  mpi_isend error'
                call flush(nfout)
                call mpi_abort(mpi_comm_world, 10002, ierr)
             endif
         endif
      enddo
                                                  __TIMER_DO_STOP(976)
#else
      do i = 1, scl_comm_rank_r
         lrank = scl_comm_rno_r(i)
         if (lrank /= mype)  icnt_send = icnt_send + 1
      enddo
#endif

      rank_c = myrank_e
      rank_r = myrank_g
                                                  __TIMER_DO_START(978)
      do k = nis_e(rank_c), nie_e(rank_c)
         ng_c = neg_g_all(k)
         if( ng_c > meg ) cycle
         if (mod(ng_c,scl_nb) > 0) then
            nn = ng_c / scl_nb + 1
         else
            nn = ng_c / scl_nb
         end if
         lrk_c = mod((nn-1),scl_col)
         lad = mod(ng_c,scl_nb)
         if (lad == 0) lad = scl_nb
         mm  = ((ng_c-1)/scl_nb)/scl_col
         lno_c = mm * scl_nb + lad

         do j = nis_eg(rank_r), nie_eg(rank_r)
            ng_r = neg_gg_all(j)
            if( ng_r > meg ) cycle
            if (mod(ng_r,scl_nb) > 0) then
               nn = ng_r / scl_nb + 1
            else
               nn = ng_r / scl_nb
            end if
            lrk_r = mod((nn-1),scl_row)
            lad = mod(ng_r,scl_nb)
            if (lad == 0) lad = scl_nb
            mm  = ((ng_r-1)/scl_nb)/scl_row
            lno_r = mm * scl_nb + lad

#ifdef EIGEN_6D
            lrk = usermap(lrk_r+1,lrk_c+1)
#else
#if defined(_ASSIGN_ROW_) || defined(USE_EIGENLIB)
            lrk = scl_row*lrk_c+lrk_r
#else
            lrk = scl_col*lrk_r+lrk_c
#endif
#endif
            if (lrk == mype) then
!fj --------------------
!fj            if(ng_r> ng_c) cycle
!fj --------------------
               if (kimg_t == 1) then
                  amat_l(j-nis_eg(rank_r)+1,k-nis_e(rank_c)+1) = zmat_l(lno_r,lno_c)
               else
                  amat_l((j-nis_eg(rank_r)+1)*2-1,k-nis_e(rank_c)+1) = zmat_l(lno_r*2-1,lno_c)
                  amat_l((j-nis_eg(rank_r)+1)*2  ,k-nis_e(rank_c)+1) = zmat_l(lno_r*2  ,lno_c)
               end if
            end if
         end do
      end do
                                                  __TIMER_DO_STOP(978)

#ifdef USE_NONBLK_COMM
      call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
       if (ierr /= 0) then
          write(nfout,*)' trans_scalapack_r :  mpi_waitall error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 10003, ierr)
       endif

      call mpi_waitall(icnt_send, req_s, sta_s, ierr)
       if (ierr /= 0) then
          write(nfout,*)' trans_scalapack_r :  mpi_waitall error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 10004, ierr)
       endif
                                                  __TIMER_COMM_STOP(980)
#else
                                                  __TIMER_COMM_START_w_BARRIER(mpi_comm,993)
#ifndef USE_ALLTOALLV
! === DEBUG by tkato 2012/06/04 ================================================
!      real(kind=DP), allocatable, dimension(:,:) :: sbuf, rbuf
! ==============================================================================
       allocate(sbuf(scl_comm_max_r*kimg_t,0:nrank_e*nrank_g-1), stat=ierr)
       allocate(rbuf(scl_comm_max  *kimg_t,0:nrank_e*nrank_g-1), stat=ierr)
       sbuf=0
       rbuf=0
       do i = 1, scl_comm_rank_r
          sbuf(:,scl_comm_rno_r(i))=send_buf(:,i)
       enddo
       call mpi_alltoall(sbuf, scl_comm_max_r*kimg_t, mpi_double_precision, &
      &                  rbuf, scl_comm_max  *kimg_t, mpi_double_precision, &
      &                                                MPI_COMM_WORLD, ierr)
       if (ierr /= 0) then
          write(nfout,*)' trans_scalapack :  mpi_alltoall(List) error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 10005, ierr)
       endif
       do i = 1, scl_comm_rank
          recv_buf(:,i)=rbuf(:,scl_comm_rno(i))
       enddo
       deallocate(sbuf)
       deallocate(rbuf)
#else
! === DEBUG by tkato 2012/06/05 ================================================
!      real(kind=DP), allocatable, dimension(:,:) :: sbuf, rbuf
!      integer, allocatable, dimension(:) :: scnt, rcnt
!      integer, allocatable, dimension(:) :: sdsp, rdsp
! ==============================================================================
       allocate(sbuf(scl_comm_max_r*kimg_t,0:nrank_e*nrank_g-1), stat=ierr)
       allocate(rbuf(scl_comm_max  *kimg_t,0:nrank_e*nrank_g-1), stat=ierr)
       allocate(scnt(0:nrank_e*nrank_g-1), stat=ierr)
       allocate(rcnt(0:nrank_e*nrank_g-1), stat=ierr)
       allocate(sdsp(0:nrank_e*nrank_g-1), stat=ierr)
       allocate(rdsp(0:nrank_e*nrank_g-1), stat=ierr)
       sbuf=0
       rbuf=0
       scnt=0
       rcnt=0
       do i = 1, scl_comm_rank_r
          lrank = scl_comm_rno_r(i)
          sbuf(:,lrank)=send_buf(:,i)
          scnt(lrank)=scl_comm_max_r*kimg_t
       enddo
       do i = 1, scl_comm_rank
          rcnt(scl_comm_rno(i))=scl_comm_max  *kimg_t
       enddo
       do i = 0, nrank_e*nrank_g - 1
          sdsp(i)=scl_comm_max_r*kimg_t*i
          rdsp(i)=scl_comm_max  *kimg_t*i
       enddo
       call mpi_alltoallv(sbuf, scnt, sdsp, mpi_double_precision, &
      &                   rbuf, rcnt, rdsp, mpi_double_precision, &
      &                                     MPI_COMM_WORLD, ierr)
       if (ierr /= 0) then
          write(nfout,*)' trans_scalapack :  mpi_alltoallv error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 10005, ierr)
       endif
       do i = 1, scl_comm_rank
          recv_buf(:,i)=rbuf(:,scl_comm_rno(i))
       enddo
       deallocate(sbuf)
       deallocate(rbuf)
       deallocate(scnt)
       deallocate(rcnt)
       deallocate(sdsp)
       deallocate(rdsp)
#endif
                                                  __TIMER_COMM_STOP(993)
#endif

      wk_rk = 0
      rank_c = myrank_e
      rank_r = myrank_g
                                                  __TIMER_DO_START(979)
      do k = nis_e(rank_c), nie_e(rank_c)
         ng_c = neg_g_all(k)
         if( ng_c > meg ) cycle
         if (mod(ng_c,scl_nb) > 0) then
            nn = ng_c / scl_nb + 1
         else
            nn = ng_c / scl_nb
         end if
         lrk_c = mod((nn-1),scl_col)
         do j = nis_eg(rank_r), nie_eg(rank_r)
            ng_r = neg_gg_all(j)
            if( ng_r > meg ) cycle
            if (mod(ng_r,scl_nb) > 0) then
               nn = ng_r / scl_nb + 1
            else
               nn = ng_r / scl_nb
            end if
            lrk_r = mod((nn-1),scl_row)
#ifdef EIGEN_6D
            lrk = usermap(lrk_r+1,lrk_c+1)
#else
#if defined(_ASSIGN_ROW_) || defined(USE_EIGENLIB)
            lrk = scl_row*lrk_c+lrk_r
#else
            lrk = scl_col*lrk_r+lrk_c
#endif
#endif
!!          if (lrk /= mype) then
               do lrank = 1, scl_comm_rank
                  if (scl_comm_rno(lrank) == lrk) then
                     wk_rk(lrank) = wk_rk(lrank) + 1
!fj --------------------
!fj                  if(ng_r> ng_c) cycle
!fj --------------------
                     if (kimg_t == 1) then
                        amat_l(j-nis_eg(rank_r)+1,k-nis_e(rank_c)+1) = recv_buf(wk_rk(lrank),lrank)
                     else
                        amat_l((j-nis_eg(rank_r)+1)*2-1,k-nis_e(rank_c)+1) = recv_buf(wk_rk(lrank)*2-1,lrank)
                        amat_l((j-nis_eg(rank_r)+1)*2  ,k-nis_e(rank_c)+1) = recv_buf(wk_rk(lrank)*2  ,lrank)
                     end if
                     exit
                  end if
               end do
!!          end if
         end do
      end do
                                                  __TIMER_DO_STOP(979)

      deallocate(req_r, stat=ierr)
      deallocate(req_s, stat=ierr)
      deallocate(sta_r, stat=ierr)
      deallocate(sta_s, stat=ierr)
      deallocate(send_buf, stat=ierr)
      deallocate(recv_buf, stat=ierr)
      deallocate(wk_rk, stat=ierr)

    end subroutine trans_scalapack_r
#endif
!#endif

!   subroutine set_block_range(ne,np,nel_p,nis_p,nie_p,map_p)
!     integer, intent(in)                     :: ne ! number of total elements
!     integer, intent(in)                     :: np ! number of ranks (or processors)
!     integer, allocatable, dimension(:) :: nel_p,nis_p,nie_p
!     integer, intent(out), dimension(ne)  :: map_p
!     integer :: j,i
!
!     allocate(nel_p(0:np-1))
!     allocate(nis_p(0:np-1))
!     allocate(nie_p(0:np-1))

!     nel_p = ne/np
!     j = mod(ne,np)
!     do i = 0, j-1
!        nel_p(i) = nel_p(i) + 1
!     end do

!     nis_p(0) = 1
!     do i = 1, np-1
!        nis_p(i)   = nis_p(i-1) + nel_p(i-1)
!        nie_p(i-1) = nis_p(i) - 1
!     end do
!     nie_p(np-1) = ne

!     j = 0
!     do i = 1, ne
!        if(nie_p(j) < i) j = j + 1
!        map_p(i) = j
!     end do

!   end subroutine set_block_range

    subroutine subspace_rotation_real_3D(is,ie,meg)
      integer, intent(in) :: is,ie,meg
#ifdef SUBMAT_DGEMM
      real(kind=DP) :: alpha, beta
      integer       :: ib,ilmta
      real(kind=DP), allocatable, dimension(:,:) :: fsr_wk, fsi_wk
#endif
      integer :: id_sname = -1
                                                  __TIMER_SUB_START(912)
      call tstatc0_begin('subspace_roation_real_3D ', id_sname)
#ifdef SUBMAT_DGEMM
      alpha=1.d0
      beta=1.d0
                                                  __TIMER_DGEMM_START(953)
!xx   call dgemm('N','N',np_g1k(ik),np_e,meg,alpha,zaj_ball(1,1,ik,1),maxval(np_g1k(:)), &
!xx  &            zmat_l_2D(1,1),meg,beta,zaj_l(1,1,ik,1),maxval(np_g1k(:)))
      call dgemm('N','N',np_g1k(ik),np_e,ie-is+1,alpha,zaj_ball(1,is,ik,1),maxval(np_g1k(:)), &
     &            zmat_l_2D(is,1),meg,beta,zaj_l(1,1,ik,1),maxval(np_g1k(:)))
      if(sw_keep_hloc_phi==ON) then
        call dgemm('N','N',np_g1k(ik),np_e,ie-is+1,alpha,zah_ball(1,is,1),maxval(np_g1k(:)), &
     &              zmat_l_2D(is,1),meg,beta,hlocphi_l(1,1,ik,1),maxval(np_g1k(:)))
      endif

      call dgemm('T','N',np_e,np_fs,ie-is+1,alpha,zmat_l_2D(is,1),neg,&
                  fsr_in(is,1),neg,beta,fsr_l(1,1,ik),np_e)
      call dgemm('T','N',np_e,np_fs,ie-is+1,alpha,zmat_l_2D(is,1),neg,&
                  fsi_in(is,1),neg,beta,fsi_l(1,1,ik),np_e)

                                                  __TIMER_DGEMM_STOP(953)
#else
                                                  __TIMER_DO_START(954)
      do ib2=1, np_e
         do ib1=1,np_g1k(ik)
            do ii=is, ie
               zaj_l(ib1,ib2,ik,1) = zaj_l(ib1,ib2,ik,1) &
              &                       + zmat_l_2D(ii,ib2) * zaj_ball(ib1,ii,ik,1)
            enddo
         enddo
      enddo
#ifdef SAVE_FFT_TIMES
      if(sw_save_fft == ON) then
         do ib2 = 1, np_e
            status_saved_phifftr(ib2,ik) = OLD
         end do
      end if
#endif
                                                  __TIMER_DO_STOP(954)
#endif
      call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(912)
    end subroutine subspace_rotation_real_3D

    subroutine subspace_rotation_imag_3D(is,ie,meg)
      integer, intent(in) :: is,ie,meg
      integer :: id_sname = -1
#ifdef SUBMAT_DGEMM
      real(kind=DP) :: alpha, beta
      integer :: ibsize
      real(kind=DP), allocatable, dimension(:,:) :: zzr,zzi
#endif
                                                  __TIMER_SUB_START(913)
      call tstatc0_begin('subspace_roation_imag_3D ', id_sname)

      if(k_symmetry(ik) == GAMMA) then
#ifdef SUBMAT_DGEMM
         ibsize = maxval(np_g1k(:))
         alpha=1.d0
         beta=1.d0
                                                  __TIMER_DGEMM_START(955)
!x       call dgemm('N','N',np_g1k(ik),np_e,meg,alpha,zaj_ball(1,1,ik,1),ibsize, &
!x      &           zmat_l_2D(1,1),meg,beta,zaj_l(1,1,ik,1),ibsize)
         call dgemm('N','N',np_g1k(ik),np_e,ie-is+1,alpha,zaj_ball(1,is,ik,1),ibsize, &
        &           zmat_l_2D(is,1),meg,beta,zaj_l(1,1,ik,1),ibsize)
!x       call dgemm('N','N',np_g1k(ik),np_e,meg,alpha,zaj_ball(1,1,ik,2),ibsize, &
!x      &           zmat_l_2D(1,1),meg,beta,zaj_l(1,1,ik,2),ibsize)
         call dgemm('N','N',np_g1k(ik),np_e,ie-is+1,alpha,zaj_ball(1,is,ik,2),ibsize, &
        &           zmat_l_2D(is,1),meg,beta,zaj_l(1,1,ik,2),ibsize)

         if(sw_keep_hloc_phi==ON) then
           call dgemm('N','N',np_g1k(ik),np_e,ie-is+1,alpha,zah_ball(1,is,1),ibsize, &
        &             zmat_l_2D(is,1),meg,beta,hlocphi_l(1,1,ik,1),ibsize)
           call dgemm('N','N',np_g1k(ik),np_e,ie-is+1,alpha,zah_ball(1,is,2),ibsize, &
        &             zmat_l_2D(is,1),meg,beta,hlocphi_l(1,1,ik,2),ibsize)
         endif


         call dgemm('T','N',np_e,np_fs,ie-is+1,alpha,zmat_l_2D(is,1),neg,&
                     fsr_in(is,1),neg,beta,fsr_l(1,1,ik),np_e)

                                                  __TIMER_DGEMM_STOP(955)
#else
                                                  __TIMER_DO_START(956)
         do ib2=1, np_e
            do ib1=1,np_g1k(ik)
               do ii=is, ie
                  zaj_l(ib1,ib2,ik,1) = zaj_l(ib1,ib2,ik,1) &
                 &                       + zmat_l_2D(ii,ib2) * zaj_ball(ib1,ii,ik,1)
                  zaj_l(ib1,ib2,ik,2) = zaj_l(ib1,ib2,ik,2) &
                 &                       + zmat_l_2D(ii,ib2) * zaj_ball(ib1,ii,ik,2)
               enddo
            enddo
         enddo
#ifdef SAVE_FFT_TIMES
         if(sw_save_fft == ON) then
            do ib2 = 1, np_e
               status_saved_phifftr(ib2,ik) = OLD
            end do
         end if
#endif
                                                  __TIMER_DO_STOP(956)
#endif
         if(iprisubmat >= 2) then
            do ib2 = is, ie
            end do
         end if
      else
#ifdef SUBMAT_DGEMM
         allocate(zzr(ie-is+1,np_e)); allocate(zzi(ie-is+1,np_e))
                                                  __TIMER_DO_START(957)
         do ib2=1,np_e
            do ib1=is,ie
               zzr(ib1-is+1,ib2)=zmat_l_2D(2*ib1-1,ib2)
               zzi(ib1-is+1,ib2)=zmat_l_2D(2*ib1  ,ib2)
            end do
         end do
                                                  __TIMER_DO_STOP(957)
         ibsize = maxval(np_g1k(:))
         alpha=1.d0
         beta=1.d0
                                                  __TIMER_DGEMM_START(958)
!x       call dgemm('N','N',np_g1k(ik),np_e,meg,alpha,zaj_ball(1,1,ik,1),ibsize, &
!x      &           zzr(1,1),meg,beta,zaj_l(1,1,ik,1),ibsize)
         call dgemm('N','N',np_g1k(ik),np_e,ie-is+1,alpha,zaj_ball(1,is,ik,1),ibsize, &
        &           zzr(1,1),ie-is+1,beta,zaj_l(1,1,ik,1),ibsize)
         alpha=-1.d0
         beta=1.d0
!x       call dgemm('N','N',np_g1k(ik),np_e,meg,alpha,zaj_ball(1,1,ik,2),ibsize, &
!x      &           zzi(1,1),meg,beta,zaj_l(1,1,ik,1),ibsize)
         call dgemm('N','N',np_g1k(ik),np_e,ie-is+1,alpha,zaj_ball(1,is,ik,2),ibsize, &
        &           zzi(1,1),ie-is+1,beta,zaj_l(1,1,ik,1),ibsize)
         alpha=1.d0
         beta=1.d0
!x       call dgemm('N','N',np_g1k(ik),np_e,meg,alpha,zaj_ball(1,1,ik,2),ibsize, &
!x      &           zzr(1,1),meg,beta,zaj_l(1,1,ik,2),ibsize)
         call dgemm('N','N',np_g1k(ik),np_e,ie-is+1,alpha,zaj_ball(1,is,ik,2),ibsize, &
        &           zzr(1,1),ie-is+1,beta,zaj_l(1,1,ik,2),ibsize)
         alpha=1.d0
         beta=1.d0
!x       call dgemm('N','N',np_g1k(ik),np_e,meg,alpha,zaj_ball(1,1,ik,1),ibsize, &
!x      &           zzi(1,1),meg,beta,zaj_l(1,1,ik,2),ibsize)
         call dgemm('N','N',np_g1k(ik),np_e,ie-is+1,alpha,zaj_ball(1,is,ik,1),ibsize, &
        &           zzi(1,1),ie-is+1,beta,zaj_l(1,1,ik,2),ibsize)


         if(sw_keep_hloc_phi==ON) then
           alpha=1.d0
           beta=1.d0
           call dgemm('N','N',np_g1k(ik),np_e,ie-is+1,alpha,zah_ball(1,is,1),ibsize, &
        &             zzr(1,1),ie-is+1,beta,hlocphi_l(1,1,ik,1),ibsize)
           alpha=-1.d0
           beta=1.d0
           call dgemm('N','N',np_g1k(ik),np_e,ie-is+1,alpha,zah_ball(1,is,2),ibsize, &
        &             zzi(1,1),ie-is+1,beta,hlocphi_l(1,1,ik,1),ibsize)
           alpha=1.d0
           beta=1.d0
           call dgemm('N','N',np_g1k(ik),np_e,ie-is+1,alpha,zah_ball(1,is,2),ibsize, &
        &             zzr(1,1),ie-is+1,beta,hlocphi_l(1,1,ik,2),ibsize)
           alpha=1.d0
           beta=1.d0
           call dgemm('N','N',np_g1k(ik),np_e,ie-is+1,alpha,zah_ball(1,is,1),ibsize, &
        &             zzi(1,1),ie-is+1,beta,hlocphi_l(1,1,ik,2),ibsize)
         endif
                                                  __TIMER_DGEMM_STOP(958)

         alpha=1.d0
         beta=1.d0
         call dgemm('T','N',np_e,np_fs,ie-is+1,alpha,zzr(1,1),neg,&
                     fsr_in(is,1),neg,beta,fsr_l(1,1,ik),np_e)

         alpha=-1.d0
         beta=1.d0
         call dgemm('T','N',np_e,np_fs,ie-is+1,alpha,zzi(1,1),neg,&
                     fsi_in(is,1),neg,beta,fsr_l(1,1,ik),np_e)

         alpha=1.d0
         beta=1.d0
         call dgemm('T','N',np_e,np_fs,ie-is+1,alpha,zzr(1,1),neg,&
                     fsi_in(is,1),neg,beta,fsi_l(1,1,ik),np_e)
         alpha=1.d0
         beta=1.d0
         call dgemm('T','N',np_e,np_fs,ie-is+1,alpha,zzi(1,1),neg,&
                     fsr_in(is,1),neg,beta,fsi_l(1,1,ik),np_e)

         deallocate(zzi,zzr)
#else
                                                  __TIMER_DO_START(959)
         do ib2=1, np_e
            do ib1=1,np_g1k(ik)
               do ii=is,ie
                  zaj_l(ib1,ib2,ik,1) = zaj_l(ib1,ib2,ik,1) &
                 &           + zmat_l_2D(ii*2-1,ib2) * zaj_ball(ib1,ii,ik,1) &
                 &           - zmat_l_2D(ii*2,  ib2) * zaj_ball(ib1,ii,ik,2)
                  zaj_l(ib1,ib2,ik,2) = zaj_l(ib1,ib2,ik,2) &
                 &           + zmat_l_2D(ii*2-1,ib2) * zaj_ball(ib1,ii,ik,2) &
                 &           + zmat_l_2D(ii*2,  ib2) * zaj_ball(ib1,ii,ik,1)
               enddo
            enddo
         enddo
#ifdef SAVE_FFT_TIMES
         if(sw_save_fft == ON) then
            do ib2 = 1, np_e
               status_saved_phifftr(ib2,ik) = OLD
            end do
         end if
#endif
                                                  __TIMER_DO_STOP(959)
#endif
       end if
       call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(913)
    end subroutine subspace_rotation_imag_3D

    subroutine set_col_partition_3D(meg)
      !!$integer,parameter :: maxmatsize = 67108864 ! 512MB
      integer, intent(in) :: meg
      integer :: maxmatsize
      integer :: trisize, size
      integer :: rowsize, block
      integer :: i,is,ie
      integer :: maxsize, nblocksize

      integer :: id_sname = -1
                                                  __TIMER_SUB_START(902)
      call tstatc0_begin('set_col_partition_3D ', id_sname)

      maxmatsize = msize_submat * 1024 * 1024 / 8 ! bytes

      trisize = meg*kimg_t*np_e
      maxsize = trisize
!fj --------------------
    if (nblocksize_submat > 0) then
      if  (mod(meg,nblocksize_submat) > 0) then
        npart = meg / nblocksize_submat + 1
      else
        npart = meg / nblocksize_submat
      end if
      if(allocated(isp)) deallocate(isp)
      if(allocated(iep)) deallocate(iep)
      allocate(isp(npart))
      allocate(iep(npart))
      isp(1) = 1
      iep(1) = isp(1) + nblocksize_submat - 1
      do i = 2, npart
        isp(i) = iep(i-1) + 1
        iep(i) = isp(i) + nblocksize_submat - 1
      end do
      iep(npart) = min(iep(npart),meg)
    else
!fj --------------------
      if(trisize <= maxmatsize.or.sw_scalapack==OFF) then
         npart = 1
         if(.not.allocated(isp)) allocate(isp(npart))
         if(.not.allocated(iep)) allocate(iep(npart))
         isp(npart) = 1
         iep(npart) = meg
         maxsize  = trisize
      else
         npart = (trisize - 1)/maxmatsize + 1
         if(.not.allocated(isp)) allocate(isp(npart))
         if(.not.allocated(iep)) allocate(iep(npart))
         maxsize  = trisize/npart
         ie = 0
                                                  __TIMER_DO_START(940)
         do i=1,npart
            is = ie+1
            ie = is
            size = 0
            do while(size*kimg_t <= maxsize)
               ie = ie+1
               size = (ie-is+1)*np_e
            end do
            isp(i) = is
            iep(i) = ie
            if(ie>meg) then
               npart = i
               exit
            end if
         end do
                                                  __TIMER_DO_STOP(940)
         iep(npart) = meg
      end if
!fj --------------------
    end if
!fj --------------------

!fj --------------------
    if(nblocksize_submat_latter_is_given) then
       nblocksize = nblocksize_submat_latter
    else if(nblocksize_submat_is_given) then
       nblocksize = nblocksize_submat
    else
       nblocksize = 0
    end if

!    block = meg
!    rowsize = meg*kimg
!    if (nblocksize > 0) then
!      if  (mod(meg,nblocksize) > 0) then
!        npart2 = meg / nblocksize + 1
!      else
!        npart2 = meg / nblocksize
!      end if
!      if(allocated(isp2)) deallocate(isp2)
!      if(allocated(iep2)) deallocate(iep2)
!      allocate(isp2(npart2))
!      allocate(iep2(npart2))
!      isp2(1) = 1
!      iep2(1) = isp2(1) + nblocksize_submat_latter - 1
!      do i = 2, npart2
!        isp2(i) = iep2(i-1) + 1
!        iep2(i) = isp2(i) + nblocksize_submat_latter - 1
!      end do
!      iep2(npart2) = meg
!    else
!fj --------------------
      rowsize = meg*kimg_t
      if(allocated(isp2)) deallocate(isp2)
      if(allocated(iep2)) deallocate(iep2)
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
                                                  __TIMER_DO_START(941)
         do i=1,npart2
            is = ie+1
            ie = is+block-1
            if(i<=mod(meg,npart2)) ie=ie+1
            isp2(i) = is
            iep2(i) = ie
         end do
                                                  __TIMER_DO_STOP(941)
         iep2(npart2) = meg
      end if
!fj --------------------
!    end if
!fj --------------------

      if (iprisubmat >= 2) then
         write(nfout,'("npart = ",i4,", nblocksize_submat= ",i4)') npart, nblocksize_submat
         write(nfout,'("isp   = ",10(i4,1x))') isp(:)
         write(nfout,'("iep   = ",10(i4,1x))') iep(:)
         write(nfout,'("npart2= ",i4,", nblocksize_submat_latter = ",i4)') npart2, nblocksize_submat_latter
         write(nfout,'("isp2  = ",10(i4,1x))') isp2(:)
         write(nfout,'("iep2  = ",10(i4,1x))') iep2(:)
      endif

      if(printable.and. sw_scalapack==ON.and.iprisubmat>=2) then
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
            write(nfout,'(4i10)') i, is, ie, size
         end do
         write(nfout,*) 'npart2=',npart2
         write(nfout,*) 'block=',block
         write(nfout,*) 'memory(MB)=',block*rowsize*8/1024/1024
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
                                                  __TIMER_SUB_STOP(902)
    end subroutine set_col_partition_3D


  subroutine set_hmat_3D(meg,w1hw2,zmat_l)
    integer, intent(in) :: meg
    real(kind=DP), intent(in) :: w1hw2((ie-is+1)*kimg_t,np_e)
#ifdef _USE_SCALAPACK_
! === DEBUG by Tkato 2011/06/28 ================================================
!   real(kind=DP), intent(out) :: zmat_l(maxeg*kimg_t,maxe)
    real(kind=DP), intent(out) :: zmat_l(lda*kimg_t,occ)
! ==============================================================================
#else
    real(kind=DP), intent(out) :: zmat_l(nel_eg(myrank_g)*kimg_t,np_e)
#endif

    real(kind=DP), allocatable :: w1hw2_mpi(:,:), zmat_mpi(:,:)
    integer(kind=4) :: irecv_num(0:nrank_g-1)
    integer(kind=4) :: i, j, G_ista_e, G_iend_e
    integer :: id_sname = -1
                                                  __TIMER_SUB_START(904)
    call tstatc0_begin('set_hmat_3D ', id_sname)

    if (npart == 1) then
       allocate(w1hw2_mpi(np_e,meg*kimg_t), stat=ierr)
        if (ierr /= 0) then
           call mpi_abort(mpi_comm_world, 122, ierr)
        end if
       w1hw2_mpi=0.d0
       allocate(zmat_mpi(np_e,meg*kimg_t), stat=ierr)
        if (ierr /= 0) then
           call mpi_abort(mpi_comm_world, 123, ierr)
        end if
        zmat_mpi=0.d0
                                                  __TIMER_DO_START(943)
       do i = 1, np_e
          do j = 1, meg*kimg_t
             w1hw2_mpi(i,j) = w1hw2(j,i)
          enddo
       enddo
                                                  __TIMER_DO_STOP(943)
       irecv_num(:) = nel_eg(:)*kimg_t*np_e
                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,944)
       call mpi_reduce_scatter(w1hw2_mpi(1,1), zmat_mpi(1,1), irecv_num, mpi_double_precision, &
      &                        mpi_sum, mpi_ke_world, ierr)
                                                  __TIMER_COMM_STOP(944)
                                                  __TIMER_DO_START(945)
       if(kimg_t==1) then
          do i = 1, np_e
             if (neg_g(i) > meg) cycle
             do j = 1, nel_eg(myrank_g)
                if (neg_gg(j) > meg) cycle
                zmat_l(j,i) = zmat_mpi(i,j)
             enddo
          enddo
       else
          do i = 1, np_e
             if (neg_g(i) > meg) cycle
             do j = 1, nel_eg(myrank_g)
                if (neg_gg(j) > meg) cycle
                zmat_l(j*2-1,i) = zmat_mpi(i,j*2-1)
                zmat_l(j*2  ,i) = -zmat_mpi(i,j*2  )
             enddo
          enddo
       endif
                                                  __TIMER_DO_STOP(945)
       deallocate(w1hw2_mpi)
       deallocate(zmat_mpi)
    else
       allocate(w1hw2_mpi((ie-is+1)*kimg_t,np_e))
                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,946)
       call mpi_allreduce(w1hw2, w1hw2_mpi,np_e*(ie-is+1)*kimg_t,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
                                                  __TIMER_COMM_STOP(946)
       G_ista_e = nis_eg(myrank_g)
       G_iend_e = nie_eg(myrank_g)
                                                  __TIMER_DO_START(947)
       if(kimg_t==1) then
          do j = 1, np_e
             if (neg_g(j) > meg) cycle
             do i = is, ie
                if (neg_gg_all(i) > meg) cycle
                if ((G_ista_e <= i) .and. (i <= G_iend_e)) then
                   zmat_l(i-G_ista_e+1,j) = w1hw2_mpi(i-is+1,j)
                endif
             enddo
          enddo
       else
          do j = 1, np_e
             if (neg_g(j) > meg) cycle
             do i = is, ie
                if (neg_gg_all(i) > meg) cycle
                if ((G_ista_e <= i) .and. (i <= G_iend_e)) then
                   zmat_l((i-G_ista_e+1)*2-1,j) =  w1hw2_mpi((i-is+1)*2-1,j)
                   zmat_l((i-G_ista_e+1)*2  ,j) = -w1hw2_mpi((i-is+1)*2  ,j)
                endif
             enddo
          enddo
       end if
                                                  __TIMER_DO_STOP(947)
       deallocate(w1hw2_mpi)
    endif
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(904)
  end subroutine set_hmat_3D

#ifdef _USE_SCALAPACK_
! ScaLapack

  subroutine set_hmat_3D_scl(neg,w1hw2,zmat_l)
    integer, intent(in) :: neg
    real(kind=DP), intent(inout) :: w1hw2((ie-is+1)*kimg_t,np_e)
!fj real(kind=DP), intent(out) :: zmat_l(lda*kimg_t,occ)
    real(kind=DP), intent(inout) ::zmat_l(maxeg*kimg_t,maxe)

    real(kind=DP), allocatable :: w1hw2_mpi(:,:), zmat_mpi(:,:)
    integer(kind=4) :: irecv_num(0:nrank_g-1)
    integer(kind=4) :: i, j, G_ista_e, G_iend_e
    integer :: id_sname = -1
                                                  __TIMER_SUB_START(963)
    call tstatc0_begin('set_hmat_3D ', id_sname)

    if (npart == 1) then
       allocate(w1hw2_mpi(np_e,neg*kimg_t), stat=ierr)
        if (ierr /= 0) then
           call mpi_abort(mpi_comm_world, 122, ierr)
        end if
       allocate(zmat_mpi(np_e,neg*kimg_t), stat=ierr)
        if (ierr /= 0) then
           call mpi_abort(mpi_comm_world, 123, ierr)
        end if
                                                  __TIMER_DO_START(966)
       do i = 1, np_e
          do j = 1, neg*kimg_t
             w1hw2_mpi(i,j) = w1hw2(j,i)
          enddo
       enddo
                                                  __TIMER_DO_STOP(966)
       irecv_num(:) = nel_eg(:)*kimg_t*np_e
                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,967)
       call mpi_reduce_scatter(w1hw2_mpi(1,1), zmat_mpi(1,1), irecv_num, mpi_double_precision, &
      &                        mpi_sum, mpi_ke_world, ierr)
                                                  __TIMER_COMM_STOP(967)
                                                  __TIMER_DO_START(968)
       if(kimg_t==1) then
          do i = 1, np_e
!            if (neg_g(i) > meg) cycle
             do j = 1, nel_eg(myrank_g)
!               if (neg_gg(j) > meg) cycle
                zmat_l(j,i) = zmat_mpi(i,j)
             enddo
          enddo
       else
          do i = 1, np_e
!            if (neg_g(i) > meg) cycle
             do j = 1, nel_eg(myrank_g)
!               if (neg_gg(j) > meg) cycle
                zmat_l(j*2-1,i) = zmat_mpi(i,j*2-1)
                zmat_l(j*2  ,i) = -zmat_mpi(i,j*2  )
             enddo
          enddo
       endif
                                                  __TIMER_DO_STOP(968)
       deallocate(w1hw2_mpi)
       deallocate(zmat_mpi)
    else
       allocate(w1hw2_mpi((ie-is+1)*kimg_t,np_e))
                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,969)
       call mpi_allreduce(w1hw2,w1hw2_mpi, np_e*(ie-is+1)*kimg_t, mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
                                                  __TIMER_COMM_STOP(969)
       G_ista_e = nis_eg(myrank_g)
       G_iend_e = nie_eg(myrank_g)
                                                  __TIMER_DO_START(970)
       if(kimg_t==1) then
          do j = 1, np_e
!            if (neg_g(j) > meg) cycle
             do i = is, ie
!               if (neg_gg_all(i) > meg) cycle
                if ((G_ista_e <= i) .and. (i <= G_iend_e)) then
                   zmat_l(i-G_ista_e+1,j) = w1hw2_mpi(i-is+1,j)
                endif
             enddo
          enddo
       else
          do j = 1, np_e
!            if (neg_g(j) > meg) cycle
             do i = is, ie
!               if (neg_gg_all(i) > meg) cycle
                if ((G_ista_e <= i) .and. (i <= G_iend_e)) then
                   zmat_l((i-G_ista_e+1)*2-1,j) =  w1hw2_mpi((i-is+1)*2-1,j)
                   zmat_l((i-G_ista_e+1)*2  ,j) = -w1hw2_mpi((i-is+1)*2  ,j)
                endif
             enddo
          enddo
       end if
                                                  __TIMER_DO_STOP(970)
       deallocate(w1hw2_mpi)
    endif
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(963)
  end subroutine set_hmat_3D_scl

! ScaLapack

! === DEBUG by Tkato 2011/06/28 ================================================
  subroutine get_zmat(lda,occ,zmat,ndim,nprow,npcol,nb,usermap,is,ie,zz)
    integer, intent(in) :: lda,occ
    real(kind=DP), intent(out) :: zmat(lda*kimg_t,occ)
    integer, intent(in) :: ndim,nprow,npcol
    integer             ::                  nb
    integer, intent(in) :: usermap(nprow,npcol)
    integer, intent(in) :: is,ie
    real(kind=DP), intent(in) :: zz(ndim*kimg_t,is:ie)

    integer :: i,j
    integer :: jj,ipcol,ilcol,ixcol
    integer :: ii,iprow,ilrow,ixrow
    integer :: jsize
!   real(kind=DP), allocatable :: zz_mpi(:,:)
    integer :: id_sname = -1
    call tstatc0_begin('get_zmat ', id_sname)

    nb = block_size

    jsize = ie-is+1

#ifdef NEC_ITER_REG
    call FTRACE_REGION_BEGIN("GET_ZMAT")
#endif
!   zz = 0.d0
    zmat = 0.d0
    if(kimg_t==1) then
       do j=is,ie
!fj --------------------
          if (j > meg) cycle
!fj --------------------
          ipcol = mod((j-1)/nb,npcol)
          ilcol = (j-1)/(npcol*nb)
          ixcol = mod(j-1,nb)+1
          jj = ilcol*nb+ixcol
          do i=1,ndim
!fj --------------------
             if (i > meg) cycle
!fj --------------------
             iprow = mod((i-1)/nb,nprow)
             ilrow = (i-1)/(nprow*nb)
             ixrow = mod(i-1,nb)+1
             ii = ilrow*nb+ixrow
             if(usermap(iprow+1,ipcol+1) /= mype) cycle
!            zz(i,j) = zmat(ii,jj)
             zmat(ii,jj) = zz(i,j)
          end do
       end do
    else
       do j=is,ie
!fj --------------------
          if (j > meg) cycle
!fj --------------------
          ipcol = mod((j-1)/nb,npcol)
          ilcol = (j-1)/(npcol*nb)
          ixcol = mod(j-1,nb)+1
          jj = ilcol*nb+ixcol
          do i=1,ndim
!fj --------------------
             if (i > meg) cycle
!fj --------------------
             iprow = mod((i-1)/nb,nprow)
             ilrow = (i-1)/(nprow*nb)
             ixrow = mod(i-1,nb)+1
             ii = ilrow*nb+ixrow
             if(usermap(iprow+1,ipcol+1) /= mype) cycle
!            zz(2*i-1,j) = zmat(2*ii-1,jj)
!            zz(2*i  ,j) = zmat(2*ii  ,jj)
             zmat(2*ii-1,jj) = zz(2*i-1,j)
             zmat(2*ii  ,jj) = zz(2*i  ,j)
          end do
       end do
    end if
!   if(npes > 1) then
!      allocate(zz_mpi(ndim*kimg_t,is:ie))
!      call mpi_allreduce(zz,zz_mpi,ndim*kimg_t*jsize,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI
!      zz = zz_mpi         ! MPI
!      deallocate(zz_mpi)
!   end if
#ifdef NEC_ITER_REG
    call FTRACE_REGION_END("GET_ZMAT")
#endif
    call tstatc0_end(id_sname)
  end subroutine get_zmat
! ==============================================================================

#ifdef EIGEN_TRANS
  subroutine scalapack_setup_3D(ndim,lwork1,lrwork1,liwork1,lwork2,lrwork2,liwork2,lda,occ,&
      & ictxt,myrow,mycol,desca,descz,usermap,icolor,icol,irank_c,irank_r)
#else
  subroutine scalapack_setup_3D(ndim,lwork1,lrwork1,liwork1,lwork2,lrwork2,liwork2,lda,occ,&
      & ictxt,myrow,mycol,desca,descb,descz,usermap,icolor,icol)
#endif

    integer, intent(in) :: ndim
#ifdef SINGLE_CONTEXT
    integer, intent(out) :: lwork1,lrwork1,liwork1,lwork2,lrwork2,liwork2,lda,occ,ictxt,myrow,mycol
#else
    integer, intent(out) :: lwork1,lrwork1,liwork1,lwork2,lrwork2,liwork2,lda,occ,ictxt(0:nrank_k-1),myrow,mycol
#endif
    integer, dimension(9), intent(out) :: desca,descb,descz
#ifdef SINGLE_CONTEXT
    integer, intent(out) :: usermap(nprow,npcol)
#else
    integer, intent(out) :: usermap(nprow,npcol,0:nrank_k-1)
#endif
    integer, intent(out) :: icolor
    integer, intent(out) :: icol(0:nrank_s*nrank_e*nrank_g*nrank_k-1)
    integer              :: ikey(0:nrank_s*nrank_e*nrank_g*nrank_k-1)
#ifdef EIGEN_TRANS
    integer, intent(out) :: irank_c(0:nrank_s*nrank_e*nrank_g*nrank_k-1)
    integer, intent(out) :: irank_r(0:nrank_s*nrank_e*nrank_g*nrank_k-1)
#endif

    integer :: i,j
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
!fj+
    integer      :: lwork0,lrwork0,liwork0
    real(DP)    :: eig
    integer     :: itmp(1)
    real(DP)    :: rtmp(1)
    complex(DP) :: ctmp(1)
    real(DP)    :: utmp2(1),vtmp2(1)
!fj-
! === Change block_size if it's too large. by tkato 2014/=======================
    integer            :: max_block_size
! ==============================================================================
                                                  __TIMER_SUB_START(903)
    call tstatc0_begin('scalapack_setup ', id_sname)

    lwork1  = 0
    lrwork1 = 0
    liwork1 = 0
    lwork2  = 0
    lrwork2 = 0
    liwork2 = 0
    desca = 0
    descb = 0
    descz = 0

! === Change block_size if it's too large. by tkato 2014/=======================
    if(max(nprow,npcol) /= 1) then
       max_block_size = int(real(ndim)/real(max(nprow,npcol) - 1))
       if(block_size > max_block_size) then
          if(mype == 0) then
             write(0,'(a)') '=== WARNING!!! =============================================='
             write(0,'(a)') 'Block size for block-cyclic division on ScaLAPACK is too large!'
             write(0,'(a,i5,a)') 'Block size should be not greater than ', max_block_size, '!'
             write(0,'(a,i5,a)') 'So, block_size is changed into ', max_block_size, '!'
             write(0,'(a)') 'FYI: '
             write(0,'(a,i8)') '   ndim: ', ndim
             write(0,'(a,i5)') '   Specified block Size: ', block_size
             write(0,'(a,2i5)') '   nprow, npcol: ', nprow, npcol
             write(0,'(a)') '=== WARNING!!! =============================================='
          end if
          block_size = max_block_size
       end if
    end if
! ==============================================================================
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

#ifdef EIGEN_TRANS
    call make_usermap(usermap,icol,ikey,irank_c,irank_r)
#else
#ifdef SINGLE_CONTEXT
    call make_usermap(usermap(:,:),icol,ikey)
#else
    call make_usermap(usermap(:,:,myrank_k),icol,ikey)
#endif
#endif
#ifdef SINGLE_CONTEXT
    call MPI_COMM_RANK(MPI_CommGroup, iam, ierr)
#else
    call MPI_COMM_RANK(mpi_ge_world, iam, ierr)
#endif
    icolor=icol(iam)

    if(iprisubmat>=2) then
       write(nfout,'("nb=",i4,", lda=",i4,", occ=",i4,", nb=",i4,", ndim=",i4,", nprow=",i4,", npcol=",i4)') &
           & nb, lda, occ, nb, ndim, nprow, npcol
    end if
!    call blacs_setup(iam,nprocs)
    call blacs_pinfo(iam,nprocs)
#ifdef SINGLE_CONTEXT
    call blacs_get(-1,0,ictxt)
    call blacs_gridmap(ictxt,usermap,nprow,nprow,npcol)
    call blacs_gridinfo(ictxt,nprow,npcol,myrow,mycol)
#else
    call blacs_get(-1,0,ictxt(myrank_k))
    call blacs_gridmap(ictxt(myrank_k),usermap(:,:,myrank_k),nprow,nprow,npcol)
    call blacs_gridinfo(ictxt(myrank_k),nprow,npcol,myrow,mycol)
#endif
    if(iprisubmat>=2) write(nfout,'("nprow,npcol,myrow,mycol=",4i5)') nprow,npcol,myrow,mycol
    if(myrow == -1) then
!xx    if(printable) write(nfout,*) 'BLACS init failed.'
!xx    stop 'BLACS init failed.'
         return
    end if
#ifdef SINGLE_CONTEXT
    call descinit(desca,ndim,ndim,nb,nb,0,0,ictxt,lda,info)
    if(iprisubmat>=2) write(nfout,*) 'descinit for amat: info=',info
    call descinit(descb,ndim,ndim,nb,nb,0,0,ictxt,lda,info)
    if(iprisubmat>=2) write(nfout,*) 'descinit for bmat: info=',info
    call descinit(descz,ndim,ndim,nb,nb,0,0,ictxt,lda,info)
    if(iprisubmat>=2) write(nfout,*) 'descinit for zmat: info=',info
#else
    call descinit(desca,ndim,ndim,nb,nb,0,0,ictxt(myrank_k),lda,info)
    if(iprisubmat>=2) write(nfout,*) 'descinit for amat: info=',info
    call descinit(descb,ndim,ndim,nb,nb,0,0,ictxt(myrank_k),lda,info)
    if(iprisubmat>=2) write(nfout,*) 'descinit for bmat: info=',info
    call descinit(descz,ndim,ndim,nb,nb,0,0,ictxt(myrank_k),lda,info)
    if(iprisubmat>=2) write(nfout,*) 'descinit for zmat: info=',info
#endif

!   if(kimg_t==1) then
       ! setup for pdsyev
       if(method_scalapack==HOUSEHOLDER.and.sw_gep == OFF) then
          IROFFA = MOD( 0, nb )
          ICOFFA = MOD( 0, nb )
          iarow = INDXG2P( 1, nb, myrow, desca( RSRC_ ), nprow )
          iacol = INDXG2P( 1, nb, mycol, desca( RSRC_ ), npcol )
          np = NUMROC( ndim+iroffa, nb, myrow, iarow, nprow )
          nq = NUMROC( ndim+icoffa, nb, mycol, iacol, npcol )
          sizemqrleft = max( (nb*(nb-1))/2, (np+nq)*nb ) + nb*nb
          if(iprisubmat>=2) write(nfout,*) 'np,nq=',np,nq
          if(iprisubmat>=2) write(nfout,*) 'sizemqrleft=',sizemqrleft
          np0 = NUMROC( max(ndim,nb,2), nb, 0, 0, nprow )
          nq0 = NUMROC( max(ndim,nb,2), nb, 0, 0, npcol )
          nrc = NUMROC( ndim, nb, myrow, 0, nprocs)
          ldc = max( 1, nrc )
          if(iprisubmat>=2) write(nfout,*) 'np0,nq0=',np0,nq0
          if(iprisubmat>=2) write(nfout,*) 'nrc,ldc=',nrc,ldc
          qrmem = 2*ndim-2
          if(iprisubmat>=2) write(nfout,'("qrmem,=",i9)') qrmem
          lwork1 = 5*ndim + ndim*ldc + max( sizemqrleft, qrmem ) + 1
! === DEBUG by tkato 2012/04/05 ================================================
          call pdsyev('V','U',ndim,utmp2,1,1,desca,eig,vtmp2,1,1 &
                   & ,descz,rtmp,-1,info)
          lwork1 = max(lwork1, nint(rtmp(1)))
! ==============================================================================
          if(iprisubmat>=2) write(nfout,'("lwork1,=",i9)') lwork1
          lrwork1 = 0
          liwork1 = 0
       else
          iarow = INDXG2P( 1, nb, myrow, desca( RSRC_ ), nprow )
          iacol = INDXG2P( 1, nb, mycol, desca( CSRC_ ), npcol )
          np0 = NUMROC( ndim, nb, myrow, iarow, nprow )
! === DEBUG by tkato 2012/04/05 ================================================
          nq0 = NUMROC( ndim, nb, mycol, iacol, npcol )
! ==============================================================================
          if(iprisubmat>=2) then
             write(nfout,*) 'iarow,iacol=',iarow,iacol
             write(nfout,*) 'np0,nq0=',np0,nq0
             write(nfout,*) 'nb_a=', desca( NB_ )
          end if
          trilwmin = 3*ndim + max( nb*( np0+1 ), 3*nb )
          lwork0  = max( 1+6*ndim+2*np0*nq0, trilwmin ) + 2*ndim
          liwork0 = 7*ndim + 8*npcol + 2
          call pdsyevd('V','U',ndim,utmp2,1,1,desca,eig,vtmp2,1,1 &
                    & ,descz,rtmp,-1,itmp,liwork1,info)
          lwork1  = max( lwork0*10, nint(rtmp(1))*10 )
          liwork1 = max( liwork0  , liwork1 )
          if(iprisubmat>=2) write(nfout,'("lwork1,liwork1,=",2i16)') lwork1,liwork1
          lrwork1 = 0
       end if
!   else
       ! setup for pzheev
       if(method_scalapack==HOUSEHOLDER) then
          np0 = NUMROC( max(ndim,nb,2), nb, 0, 0, nprow )
          nq0 = NUMROC( max(ndim,nb,2), nb, 0, 0, npcol )
          if(iprisubmat>=2) write(nfout,*) 'np0,nq0=',np0,nq0
          lwork2 = ndim*(ndim+3) + ( np0 + nq0 + nb ) * nb
          lrwork2 = 2*ndim + 2*ndim-2
          if(iprisubmat>=2) write(nfout,'("lwork2,lrwork2,=",3i9)') lwork2,lrwork2
          liwork2 = 0
       else
          iarow = INDXG2P( 1, nb, myrow, desca( RSRC_ ), nprow )
          iacol = INDXG2P( 1, nb, mycol, desca( CSRC_ ), npcol )
          np0 = NUMROC( ndim, nb, myrow, iarow, nprow )
          nq0 = NUMROC( ndim, nb, mycol, iacol, npcol )
          nq  = NUMROC( MAX( ndim, nb, 2 ), nb, 0, 0, npcol )
          if(iprisubmat>=2) then
             write(nfout,*) 'iarow,iacol=',iarow,iacol
             write(nfout,*) 'np0,nq0,nq=',np0,nq0,nq
             write(nfout,*) 'nb_a=', desca( NB_ )
          end if
          lwork0  = ndim + ( np0 + nq + nb ) * nb
!fj       lrwork0 = 1 + 9*ndim + 3*np0*nq0
          lrwork0 = (1 + 8*ndim + 2*np0*nq0)*2
          liwork0 = 7*ndim + 8*npcol + 2

          call pzheevd('V','U',ndim,utmp2,1,1,desca,eig,vtmp2,1,1 &
                    & ,descz,ctmp,-1,rtmp,-1,itmp,-1,info)
          lwork2  = max( lwork0 , nint(real(ctmp(1))) )
          lrwork2 = max( lrwork0, nint(rtmp(1)) )
          liwork2 = max( liwork0, itmp(1) )

          lrwork2 =lrwork2*10
          if(iprisubmat>=2) write(nfout,'("lwork2,lrwork2,liwork2,=",3i9)') lwork2,lrwork2,liwork2
       end if
!   end if
    if(sw_gep == ON)then
    endif
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(903)
  end subroutine scalapack_setup_3D

  subroutine pdsyev_driver_3D(ndim,eig,desca,descz,lda,occ,amat,zmat,lwork,liwork)
    integer, intent(in):: ndim
    real(kind=DP), intent(out) ,dimension(ndim) :: eig
    integer, dimension(9), intent(in) :: desca,descz
    integer, intent(in) :: lda,occ
!   real(kind=DP), intent(inout) ,dimension(maxeg,maxe) :: amat
!   real(kind=DP), intent(out) ,dimension(maxeg,maxe) :: zmat
! === DEBUG by Tkato 2011/06/28 ================================================
!   real(kind=DP), intent(inout) ,dimension(nel_eg(myrank_g),np_e) :: amat
!   real(kind=DP), intent(out) ,dimension(nel_eg(myrank_g),np_e) :: zmat
    real(kind=DP), intent(inout) ,dimension(lda,occ) :: amat
    real(kind=DP), intent(out) ,dimension(lda,occ) :: zmat
! ==============================================================================
    integer, intent(in):: lwork,liwork

    character(len=1) :: JOBZ,UPLO
    real(kind=DP),allocatable,dimension(:) :: work
    integer,allocatable,dimension(:) :: iwork
    integer :: info
    integer :: id_sname = -1
    call tstatc0_begin('pdsyev_driver ', id_sname)
                                                  __TIMER_SUB_START(906)

    if (iprisubmat >= 2) then
       write(nfout,'("-- pdsyev_driver_3D --")')
       write(nfout,'("ndim=",i4,", lda=",i4,", occ=",i4)') ndim,lda,occ
       write(nfout,'("desca : ",9(i4," ,"))') desca
       write(nfout,'("descz : ",9(i4," ,"))') descz
       call flush(nfout)
    end if

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
                                                  __TIMER_SUB_STOP(906)
  end subroutine pdsyev_driver_3D

  subroutine pzheev_driver_3D(ndim,eig,desca,descz,lda,occ,amat,zmat,lwork,lrwork,liwork)
    integer, intent(in):: ndim
    real(kind=DP), intent(out) ,dimension(ndim) :: eig
    integer, dimension(9), intent(in) :: desca,descz
    integer, intent(in) :: lda,occ
!   real(kind=DP), intent(inout) ,dimension(maxeg*2,maxe) :: amat
!   real(kind=DP), intent(out) ,dimension(maxeg*2,maxe) :: zmat
! === DEBUG by Tkato 2011/06/28 ================================================
!   real(kind=DP), intent(inout) ,dimension(nel_eg(myrank_g)*2,np_e) :: amat
!   real(kind=DP), intent(out) ,dimension(nel_eg(myrank_g)*2,np_e) :: zmat
    real(kind=DP), intent(inout) ,dimension(lda*2,occ) :: amat
    real(kind=DP), intent(out) ,dimension(lda*2,occ) :: zmat
! ==============================================================================
    integer, intent(in):: lwork,lrwork,liwork

    character(len=1) :: JOBZ,UPLO
    complex(kind=CMPLDP),allocatable,dimension(:) :: work
    real(kind=DP),allocatable,dimension(:) :: rwork
    integer,allocatable,dimension(:) :: iwork
    integer :: info

    integer :: id_sname = -1
                                                  __TIMER_SUB_START(907)
    call tstatc0_begin('pzheev_driver ', id_sname)

    if(iprisubmat>=2) write(nfout,*) "Submat: pzheev_driver_3D"
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
                                                  __TIMER_SUB_STOP(907)
  end subroutine pzheev_driver_3D

  subroutine pdsygvx_driver_3D(ndim,eig,desca,descb,descz,lda,occ,amat,smat,zmat,lwork,liwork)
    integer, intent(in):: ndim
    real(kind=DP), intent(out) ,dimension(ndim) :: eig
    integer, dimension(9), intent(in) :: desca,descb,descz
    integer, intent(in) :: lda,occ
!   real(kind=DP), intent(inout) ,dimension(maxeg,maxe) :: amat
!   real(kind=DP), intent(out) ,dimension(maxeg,maxe) :: zmat
! === DEBUG by Tkato 2011/06/28 ================================================
!   real(kind=DP), intent(inout) ,dimension(nel_eg(myrank_g),np_e) :: amat
!   real(kind=DP), intent(out) ,dimension(nel_eg(myrank_g),np_e) :: zmat
    real(kind=DP), intent(inout) ,dimension(lda,occ) :: amat
    real(kind=DP), intent(inout) ,dimension(lda,occ) :: smat
    real(kind=DP), intent(out) ,dimension(lda,occ) :: zmat
! ==============================================================================
    integer, intent(in):: lwork,liwork

    character(len=1) :: JOBZ,UPLO
    real(kind=DP),allocatable,dimension(:) :: work
    integer,allocatable,dimension(:) :: iwork
    real(kind=DP) :: vl,vu
    integer :: il,iu
    integer :: info
    integer :: id_sname = -1
    integer :: m,nz
    integer,allocatable,dimension(:) :: ifail
    integer,allocatable,dimension(:) :: iclustr
    real,allocatable,dimension(:) :: gap
    integer,save :: lwork0=-1
    integer,save :: liwork0=-1
    real(kind=DP) :: ABSTOL
    ABSTOL=-1.d0
    call tstatc0_begin('pdsyev_driver ', id_sname)
                                                  __TIMER_SUB_START(906)

    if (iprisubmat >= 2) then
       write(nfout,'("-- pdsygvx_driver_3D --")')
       write(nfout,'("ndim=",i4,", lda=",i4,", occ=",i4)') ndim,lda,occ
       write(nfout,'("desca : ",9(i4," ,"))') desca
       write(nfout,'("descb : ",9(i4," ,"))') descb
       write(nfout,'("descz : ",9(i4," ,"))') descz
       call flush(nfout)
    end if

!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
    JOBZ = 'V'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
    UPLO = 'U'
    if(liwork0<0)then
       allocate(ifail(ndim)); ifail = 0
       allocate(iclustr(npcol*nprow*2));iclustr=0
       allocate(gap(npcol*nprow*2));iclustr=0
       allocate(work(lwork))
       allocate(iwork(liwork))
       call pdsygvx(1,JOBZ,'A',UPLO,ndim,amat,1,1,desca,smat,1,1,descb  &
                & ,vl,vu,il,iu,abstol,m,nz,eig,-1,zmat,1,1,descz,work,lwork,iwork,liwork &
                & ,ifail,iclustr,gap,info)
       liwork0 = iwork(1)
       deallocate(ifail,iclustr,gap,work,iwork)
    endif
    if(lwork0<0)then
       allocate(ifail(ndim)); ifail = 0
       allocate(iclustr(npcol*nprow*2));iclustr=0
       allocate(gap(npcol*nprow*2));iclustr=0
       allocate(work(lwork))
       allocate(iwork(liwork))
       call pdsygvx(1,JOBZ,'A',UPLO,ndim,amat,1,1,desca,smat,1,1,descb  &
                & ,vl,vu,il,iu,abstol,m,nz,eig,-1,zmat,1,1,descz,work,lwork,iwork,liwork &
                & ,ifail,iclustr,gap,info)
       lwork0 = work(1)
       deallocate(ifail,iclustr,gap,work,iwork)
    endif
    allocate(ifail(ndim)); ifail = 0
    allocate(iclustr(npcol*nprow*2));iclustr=0
    allocate(gap(npcol*nprow*2));iclustr=0
    allocate(work(lwork0))
    allocate(iwork(liwork0))

#ifdef NEC_ITER_REG
    call FTRACE_REGION_BEGIN("PDSYGVX")
#endif
    call pdsygvx(1,JOBZ,'A',UPLO,ndim,amat,1,1,desca,smat,1,1,descb  &
             & ,vl,vu,il,iu,abstol,m,nz,eig,-1,zmat,1,1,descz,work,lwork0,iwork,liwork0 &
             & ,ifail,iclustr,gap,info)
#ifdef NEC_ITER_REG
    call FTRACE_REGION_END("PDSYGVX")
#endif
       if(iprisubmat>=2) write(nfout,*) 'pdsygvx: info=',info

    deallocate(ifail)
    deallocate(iclustr)
    deallocate(gap)
    deallocate(work)
    deallocate(iwork)

    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(906)
  end subroutine pdsygvx_driver_3D

  subroutine pzhegvx_driver_3D(ndim,eig,desca,descb,descz,lda,occ,amat,smat,zmat,lwork,lrwork,liwork)
    integer, intent(in):: ndim
    real(kind=DP), intent(out) ,dimension(ndim) :: eig
    integer, dimension(9), intent(in) :: desca,descb,descz
    integer, intent(in) :: lda,occ
!   real(kind=DP), intent(inout) ,dimension(maxeg,maxe) :: amat
!   real(kind=DP), intent(out) ,dimension(maxeg,maxe) :: zmat
! === DEBUG by Tkato 2011/06/28 ================================================
!   real(kind=DP), intent(inout) ,dimension(nel_eg(myrank_g),np_e) :: amat
!   real(kind=DP), intent(out) ,dimension(nel_eg(myrank_g),np_e) :: zmat
    real(kind=DP), intent(inout) ,dimension(lda,occ) :: amat
    real(kind=DP), intent(inout) ,dimension(lda,occ) :: smat
    real(kind=DP), intent(out) ,dimension(lda,occ) :: zmat
! ==============================================================================
    integer, intent(in):: lwork,lrwork,liwork

    character(len=1) :: JOBZ,UPLO
    complex(kind=DP),allocatable,dimension(:) :: work
    integer,allocatable,dimension(:) :: iwork
    real(kind=DP) :: vl,vu
    integer :: il,iu
    integer :: info
    integer :: i
    integer :: id_sname = -1
    integer,allocatable,dimension(:) :: ifail
    integer,allocatable,dimension(:) :: iclustr
    real(kind=DP), allocatable, dimension(:) :: gap
    real(kind=DP),allocatable,dimension(:) :: rwork
    integer :: lrwrok
    integer :: m,nz
    integer,save :: lwork0=-1
    integer,save :: lrwork0=-1
    integer,save :: liwork0=-1
    real(kind=DP) :: ABSTOL
    ABSTOL=-1.d0
    call tstatc0_begin('pzhegvx_driver ', id_sname)
                                                  __TIMER_SUB_START(906)

!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
    JOBZ = 'V'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
    UPLO = 'U'
    if(lwork0<0)then
        allocate(ifail(ndim)); ifail = 0
        allocate(iclustr(npcol*nprow*2));iclustr=0
        allocate(gap(npcol*nprow));gap=0
        allocate(work(lwork))
        allocate(iwork(1))
        allocate(rwork(lrwork))
        call pzhegvx(1,JOBZ,'A',UPLO,ndim,amat,1,1,desca,smat,1,1,descb  &
             & ,vl,vu,il,iu,abstol,m,nz,eig,-1,zmat,1,1,descz,work,lwork0,rwork &
             & , lrwork, iwork, liwork, ifail,iclustr,gap,info)
        lwork0 = work(1)
        deallocate(ifail,iclustr,gap,work,iwork,rwork)
    endif
    if(liwork0<0)then
        allocate(ifail(ndim)); ifail = 0
        allocate(iclustr(npcol*nprow*2));iclustr=0
        allocate(gap(npcol*nprow));gap=0
        allocate(work(lwork))
        allocate(iwork(1))
        allocate(rwork(lrwork))
        call pzhegvx(1,JOBZ,'A',UPLO,ndim,amat,1,1,desca,smat,1,1,descb  &
             & ,vl,vu,il,iu,abstol,m,nz,eig,-1.d0,zmat,1,1,descz,work,lwork0,rwork &
             & , lrwork, iwork, liwork0, ifail,iclustr,gap,info)
        liwork0 = iwork(1)
        deallocate(ifail,iclustr,gap,work,iwork,rwork)
    endif
    if(lrwork0<0)then
        allocate(ifail(ndim)); ifail = 0
        allocate(iclustr(npcol*nprow*2));iclustr=0
        allocate(gap(npcol*nprow));gap=0
        allocate(work(lwork))
        allocate(iwork(1))
        allocate(rwork(lrwork))
        call pzhegvx(1,JOBZ,'A',UPLO,ndim,amat,1,1,desca,smat,1,1,descb  &
             & ,vl,vu,il,iu,abstol,m,nz,eig,-1,zmat,1,1,descz,work,lwork0,rwork &
             & , lrwork, iwork, lrwork0, ifail,iclustr,gap,info)
        lrwork0 = rwork(1)
        deallocate(ifail,iclustr,gap,work,iwork,rwork)
    endif
    if (iprisubmat >= 2) then
       write(nfout,'("-- pdsygvx_driver_3D --")')
       write(nfout,'("ndim=",i4,", lda=",i4,", occ=",i4)') ndim,lda,occ
       write(nfout,'("desca : ",9(i4," ,"))') desca
       write(nfout,'("descb : ",9(i4," ,"))') descb
       write(nfout,'("descz : ",9(i4," ,"))') descz
       call flush(nfout)
    end if

    allocate(ifail(ndim)); ifail = 0
    allocate(iclustr(npcol*nprow*2));iclustr=0
    allocate(gap(npcol*nprow));gap=0
    allocate(work(lwork0))
    allocate(iwork(liwork0))
    allocate(rwork(lrwork0))
    !if(mype == 0)then
    !do i=1,occ
    !   do j=1,lda
    !      write(0,'(2i5,4f15.10)') j,i,amat(2*j-1,i),amat(2*j,i),smat(2*j-1,i),smat(2*j,i)
    !   enddo
    !enddo
    !endif
    call pzhegvx(1,JOBZ,'A',UPLO,ndim,amat,1,1,desca,smat,1,1,descb  &
             & ,vl,vu,il,iu,abstol,m,nz,eig,-1,zmat,1,1,descz,work,lwork0,rwork &
             & , lrwork0, iwork, liwork0, ifail,iclustr,gap,info)
    deallocate(ifail)
    deallocate(iclustr)
    deallocate(gap)
    deallocate(work)
    deallocate(iwork)
    deallocate(rwork)
    !if(mype == 0)then
    !   do i=1,ndim
    !   write(0,*) i,eig(i)
    !   enddo
    !   write(0,*) 'info',info
    !endif
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(906)
  end subroutine pzhegvx_driver_3D
#endif
! _USE_SCALAPACK_

  end subroutine evolve_WFs_in_subspace_3D

!!$#endif




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

!******************** modified by RIST_16
#ifndef _USE_DSYGVXD_
  subroutine dsygvx_driver(ndim,eig,w1hw2,w1sw2)
    integer, intent(in):: ndim
    real(kind=DP), intent(out) ,dimension(ndim) :: eig
    real(kind=DP), intent(inout) ,dimension(ndim,ndim) :: w1hw2,w1sw2

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

    integer :: ITYPE,ldb,IL,IU,ldz !,liwork
    character(len=1) :: RANGE
    real(kind=DP) ::  VL,VU,ABSTOL
    real(kind=DP),allocatable,dimension(:,:) :: z
    integer,allocatable,dimension(:) :: ifail
    integer,allocatable,dimension(:) :: iwork

    integer :: it1,it2
    integer :: meig
    integer :: id_sname = -1
    call tstatc0_begin('dsygvx_driver ', id_sname)

!write(90,*) 'fail ::    dsygvx_driver'

    ITYPE=1
    RANGE='A'
    ABSTOL=-1.d0

    lda=ndim
    ldb=ndim
    ldz=ndim

    allocate(iwork(5*ndim)); iwork = 0.d0
    allocate(z(ldz,ndim)); z = 0.d0
    allocate(ifail(ndim)); ifail = 0

!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
    JOBZ = 'V'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
    UPLO = 'U'

! (define lwork)
!!$    lwork = 0
!!$    lwork_min=max(1,3*ndim-1)

!**************************************testtesttest
lwork  = 1 + ndim*(6+2*ndim)
!**************************************testtesttest
    if(lwork == 0) then
       lwork = -1
       allocate(work_lapack(1))
!       call dsyev(JOBZ,UPLO,ndim,w1hw2,lda,eig,work_lapack,lwork,info)
       call dsygvx(ITYPE,JOBZ,RANGE,UPLO,ndim,w1hw2,lda,w1sw2,ldb,VL,VU,IL,IU,ABSTOL &
                  ,meig,eig,z,ldz,work_lapack,LWORK,IWORK,ifail,info)
       lwork = int(work_lapack(1))+1
!write(90,*) "lwork=",lwork,work_lapack(1)
       if(iprisubmat >= PRINTLEVEL) then
          write(nfout,'(" lwork = ",i8)') lwork
       end if
       deallocate(work_lapack)
    end if

    lwork_min=max(1,(nb+2)*ndim)
!lwork_min  = 1 + ndim*(6+2*ndim)
    if(lwork < lwork_min) lwork=lwork_min
    if(lwork > lwork_huge) lwork=lwork_min

!write(90,*) "lwork=",lwork, lwork_min,lwork_huge, 1 + ndim*(6+2*ndim)

    allocate(work_lapack(lwork))
    if(iprisubmat >= PRINTLEVEL .and. lwork == lwork_min) then
       write(nfout,'(" lwork = ",i8)') lwork
       write(nfout,'(" lda, ndim = ",2i8)') lda, ndim
    end if

!    call dsyev(JOBZ,UPLO,ndim,w1hw2,lda,eig,work_lapack,lwork,info)
    call dsygvx(ITYPE,JOBZ,RANGE,UPLO,ndim,w1hw2,lda,w1sw2,ldb,VL,VU,IL,IU,ABSTOL &
               ,ndim,eig,z,ldz,work_lapack,LWORK,IWORK,ifail,info)

!do it1=1,ndim
!  write(80,*) it1,eig(it1),info
!end do

    w1hw2=z

!!$    write(*,*) 'lwork& info & suggested lwork = ',lwork,info,work_lapack(1)
    if(iprisubmat >= 2) then
       write(nfout,'(" info = ",i8)') info
    end if
!!$    lwork= int(work_lapack(1))+1
    deallocate(work_lapack)

    deallocate(z)
    deallocate(ifail)
    deallocate(iwork)

    call tstatc0_end(id_sname)
  end subroutine dsygvx_driver
#else
  subroutine dsygvx_driver(ndim,eig,w1hw2,w1sw2)
    integer, intent(in):: ndim
    real(kind=DP), intent(out) ,dimension(ndim) :: eig
    real(kind=DP), intent(inout) ,dimension(ndim,ndim) :: w1hw2,w1sw2

    character(len=1) :: JOBZ,UPLO
!   integer, parameter :: lda=ndim
    integer :: lda
    integer :: lwork,liwork
    real(kind=DP),allocatable,dimension(:) :: work
    integer,allocatable,dimension(:) :: iwork
    integer :: info

    integer :: id_sname = -1
    integer :: meig
    integer :: ITYPE,ldb,IL,IU,ldz
    character(len=1) :: RANGE
    real(kind=DP) :: VL,VU,ABSTOL
    real(kind=DP),allocatable,dimension(:,:) :: z
    real(kind=DP),allocatable,dimension(:,:) :: ifail

    call tstatc0_begin('dsygvx_driver ', id_sname)

!write(90,*) 'dsygvx_driver'

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
! (define lwork)

    lwork  = 1 + ndim*(6+2*ndim)
    liwork = 5*ndim !3 + 5*ndim

    allocate(work(lwork)); work = 0.d0
    allocate(iwork(liwork)); iwork = 0.d0

    allocate(z(ldz,ndim)); z = 0.d0
    allocate(ifail(ndim)); ifail = 0

#ifdef NEC_ITER_REG
    call FTRACE_REGION_BEGIN("DSYGVXD")
#endif
!    call dsyevd(JOBZ,UPLO,ndim,w1hw2,lda,eig,work,lwork,iwork,liwork,info)
    call dsygvx(ITYPE,JOBZ,RANGE,UPLO, ndim,w1hw2,lda,w1sw2,ldb,VL,VU,IL,IU,ABSTOL &
               ,meig,eig,z,ldz,WORK,LWORK,IWORK,ifail,info)
    w1hw2=z
#ifdef NEC_ITER_REG
    call FTRACE_REGION_END("DSYGVXD")
#endif

    deallocate(work)
    deallocate(iwork)

    deallocate(z)
    deallocate(ifail)

    call tstatc0_end(id_sname)
  end subroutine dsygvx_driver
#endif

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



#ifdef _USE_SCALAPACK_
  subroutine set_nprow_npcol(nprow,npcol)
    integer, intent(inout) :: nprow,npcol
!
!   if(nprow /= nrank_e) then
!      nprow = 1
!      npcol = nrank_e
!      if(printable) write(nfout,'("set nprow,npcol=",2i5)') nprow,npcol
!   else
       if(iprisubmat >= 2) write(nfout,'("nprow,npcol=",2i5)') nprow,npcol
!   end if
  end subroutine set_nprow_npcol
#endif


#ifndef _USE_NUMPAC_
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
#endif

! ------------------------------

end module m_ES_WF_by_submat
