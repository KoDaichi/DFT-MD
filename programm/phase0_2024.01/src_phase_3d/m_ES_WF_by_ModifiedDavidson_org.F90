#define NEC_TUNE
!#undef NEC_TUNE
! AAS for Modified Davidson and Modified Kosugi
module m_ES_WF_by_ModifiedDavidson
#ifdef TRANSPOSE
#ifdef VPP
#define _ODD_BOUNDARY_
#endif
#ifdef SX
#define _ODD_BOUNDARY_
#endif
#ifdef NEC_TUNE1
#define _ODD_BOUNDARY_
#endif
  use m_Const_Parameters,    only : DP,SP,DIRECT,ON,OFF,SCF,GAMMA,OTHER_BANDS,ALL_BANDS &
        &                         , SmallestPositiveNumber,SKIP,EXECUT,zi, ELECTRON
  use m_Parallelization,     only : MPI_CommGroup &
       &                          , myrank_e,myrank_k,map_e,map_k,ista_e,iend_e,istep_e &
       &                          , ista_k,iend_k,np_g1k,ista_g1,mpi_k_world,ierr,map_z &
       &                          , np_e,npes,nrank_e,mype
  use m_Control_Parameters,  only : nspin,iprimddavidson,kimg,neg,af,npartition_mddavid &
       &                          , max_iter_mddavid, delta_eig_occup_md, delta_eig_empty_md &
       &                          , eps_mddavid, eps_residual_mddavid &
       &                          , iprimdkosugi, npartition_mdkosugi, max_iter_mdkosugi &
       &                          , delta_eig_occup_mdkosugi, delta_eig_empty_mdkosugi &
       &                          , eps_mdkosugi, eps_residual_mdkosugi
  use m_Files,               only : nfout
  use m_Timing,              only : tstatc0_begin, tstatc0_end
  use m_FFT,                 only : nfft,fft_box_size_WF, m_FFT_Vlocal_W_3D, m_FFT_Direct_3D
  use m_ES_WF_by_SDorCG,     only : map_fft_to_WF_3D
  use m_Kpoints,             only : kv3,vkxyz, k_symmetry
  use m_PlaneWaveBasisSet,   only : kg1, iba, igf, nbase, m_pwBS_kinetic_energies_3D
!!$  use m_Electronic_Structure,only : zaj_l, neordr, sc, nrvf_ordr, eko_l, vlhxcQ &
  use m_Electronic_Structure,only : zaj_l, neordr, nrvf_ordr, eko_l, vlhxcQ &
       &                          , occup_l &
       &                          , fsr_l,fsi_l, vnlph_l, vlhxc_l &
       &                          , m_ES_Vlocal_in_Rspace_3D &
       &                          , m_ES_WF_in_Rspace_3D &
       &                          , m_ES_wd_zaj_small_portion_3D &
       &                          , m_ES_wd_eko_3D &
       &                          , m_ES_sort_eigen_values_3D
  use m_ES_ortho,           only : np_g1k_x                                          &
       &                         , m_ES_orthogonalize_SD_to_WFs_3D
  use m_ES_nonlocal,        only : m_ES_Vnonlocal_W_3D                               &
       &                         , m_ES_betar_dot_WFs_4_each_k_3D                    &
       &                         , m_ES_alloc_scss_etc_3D                            &
       &                         , m_ES_dealloc_scss_etc                             &
       &                         , m_ES_betar_dot_Psi_4_each_k_3D
  use m_Ionic_System,       only : natm, iwei, ityp, ntyp
  use m_PseudoPotential,    only : ilmt,nlmta,lmta,q,dion &
       &                         , lmtt,ltp,mtp &
       &                         , m_PP_include_vanderbilt_pot &
       &                         , ipaw,dion_paw,modnrm
  use m_NonLocal_Potential, only : snl
  use m_Parallelization,    only : nel_fft_x , nel_fft_y, nel_fft_z &
                                 , mp_g1k, myrank_g
  use m_Parallelization,    only : ista_g1k, iend_g1k, neg_g, mpi_ke_world &
       &                         , ista_k, iend_k, mpi_kg_world            &
       &                         , np_fs, myrank_g, nis_fs
#ifdef NEC_TUNE
  use m_Parallelization,    only : ista_atm, iend_atm, ista_fs, iend_fs, np_fs
#endif


! ============================== added by K. Tagami ================== 11.0
  use m_Const_Parameters,   only : CMPLDP
  use m_Control_Parameters,  only : ndim_spinor, noncol, ndim_chgpot
  use m_PseudoPotential,     only : q_noncl, dion_scr_noncl
  use m_Electronic_Structure,  only : m_ES_Vlocal_in_Rspace_noncl

!  use m_FFT,                 only : m_FFT_Vlocal_W_noncl
!  use m_ES_ortho,              only : m_ES_orthogonl_SD_to_WFs_noncl
!  use m_Electronic_Structure,      only : m_ES_sort_eigen_vals_noncl
! ==================================================================== 11.0

  implicit none
  private

  integer, allocatable, dimension(:) :: nsize_subspace, nsize_matrix
  integer, allocatable, dimension(:) :: ista_e_l,iend_e_l,ielm_e_l
  integer :: nsize_sb_now, nsize_mt_now, nsize_mt_old,msize_matrix
  integer :: nblock,msize_subspace
  real(kind=DP), allocatable, dimension(:) :: w1hw2,w1hw2_mpi
  real(kind=DP), allocatable, dimension(:) :: w1sw2,w1sw2_mpi
  real(kind=DP), allocatable, dimension(:,:,:,:) :: zat_l
  real(kind=DP), allocatable, dimension(:,:,:,:) :: zah_l
  real(kind=DP), allocatable, dimension(:,:,:) :: fsr,fsi
  real(kind=DP), allocatable, dimension(:,:,:) :: zajold_l
  real(kind=DP), allocatable, dimension(:,:,:) :: zaj_l_backup
  real(kind=DP), allocatable, dimension(:,:) :: fsrold_l,fsiold_l
  real(kind=DP), allocatable, dimension(:,:,:,:):: wfsd_l  !d(kg1,np_e,ik:ik,kimg)
  real(kind=DP), allocatable, dimension(:,:,:)  :: bsdr_l, bsdi_l !d(np_e,nlmta,1)
  real(kind=DP) :: eps_residual
  logical, allocatable, dimension(:)   :: feigconv
  integer, allocatable, dimension(:,:) :: ibover


  public :: m_ESmddavid_Renew_WF
  public :: m_ESmddavid_Subspace_Rotation
  public :: m_ESmdkosugi_Renew_WF


! ================================= added by K. Tagami ================= 11.0
  real(kind=DP), allocatable, dimension(:,:,:,:,:) :: zat_l_noncl
  real(kind=DP), allocatable, dimension(:,:,:,:,:) :: zah_l_noncl
  real(kind=DP), allocatable, dimension(:,:,:,:) :: fsr_noncl, fsi_noncl
  real(kind=DP), allocatable, dimension(:,:,:,:) :: zaj_l_backup_noncl
!
  real(kind=DP), allocatable, dimension(:,:,:,:) :: zajold_l_noncl
  real(kind=DP), allocatable, dimension(:,:,:) :: fsrold_l_noncl,fsiold_l_noncl
! ===================================================================== 11.0

  include 'mpif.h'

contains

  subroutine m_ESmddavid_Renew_WF(nfout,precon)
    integer, intent(in) :: nfout,precon
    integer             :: ispin, ik, iksnl, switch_of_eko_part
    integer :: iblock,itot
    real(kind=DP), allocatable, dimension(:) ::  afft, bfft
!!$    real(kind=DP), dimension(kg1) :: ekin,vnldi,hdiag,sdiag
    real(kind=DP), allocatable, dimension(:) :: ekin,p
    real(kind=DP), allocatable, dimension(:) :: afft_l
    real(kind=DP), allocatable, dimension(:,:) :: wk_bfft_l
    real(kind=DP), allocatable, dimension(:,:) :: bfft_l
    integer :: lsize, ibsize, isrsize, fft_l_size
    real(kind=DP), allocatable, dimension(:) :: ekin_l,p_l
    logical :: frestart
    integer :: max_itot
    integer :: iblock_now, itot_now, ipri0
    integer :: n_unconv
#ifdef FFT_USE_SSL2_PAD
    lsize = max(nel_fft_x(myrank_g),nel_fft_y(myrank_g),nel_fft_z(myrank_g))
#else
    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
#endif
    allocate(afft_l(lsize*kimg), stat=ierr)
    if(ierr /= 0) then
       write(nfout,*)' m_ESmddavid_Renew_WF: Not allocated afft_l array'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 201, ierr)
    endif
    ibsize = 1

    max_itot = 1

    allocate(ekin_l(maxval(np_g1k)),p_l(maxval(np_g1k)))
    call m_ES_alloc_scss_etc_3D()
    allocate(afft(nfft)); allocate(bfft(nfft))

    call allocate_matrix
    do ispin = 1, nspin, (af+1)
       call m_ES_Vlocal_in_Rspace_3D(ispin,afft_l,lsize,1,OFF)      ! (ptfft1) vlhxc_l->afft
       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) /= myrank_k) cycle          ! MPI
          isrsize = min(lsize,mp_g1k(ik))
          fft_l_size  = nel_fft_x(myrank_g)
          allocate(wk_bfft_l(lsize*kimg,ibsize) ,stat=ierr)
          allocate(bfft_l(lsize*kimg,ibsize) ,stat=ierr)
          if (ierr /= 0) then
             write(nfout,*)' m_ESmddavid_Renew_WF:  Not allocate '
             call flush(nfout)
             call mpi_abort(mpi_comm_world, 205, ierr)
          endif
          iksnl = (ik-1)/nspin + 1
          call allocate_t_matrix_3D(ik)
          call m_pwBS_kinetic_energies_3D(ik,vkxyz,ekin_l) ! (diakin) ->ekin
          max_itot=0
          feigconv = .false.
          Loop: do itot=1,max_iter_mddavid
             itot_now = itot
             if(itot>max_itot) max_itot=itot
             call m_ES_Vnonlocal_W_3D(ik,iksnl,ispin,switch_of_eko_part=ON) ! -> vnlph_l
             call decide_correction_vector_3D(precon,ik,ekin_l,afft_l,bfft_l, &
                                              wk_bfft_l,lsize,ibsize,isrsize,fft_l_size,p_l)   ! -> wfsd_l
             call prepare_Hloc_phi_3D(ik,ekin_l,afft_l,bfft_l, &
                                      wk_bfft_l,lsize,ibsize,isrsize,fft_l_size,itot)
!print *,ik,itot
write(nfout,*) 'B',feigconv
             Block_Loop: do iblock=1,nblock
                iblock_now = iblock
                call evolve_WFs_in_subspace_3D &
                          (ik,ispin,ekin_l,afft_l,bfft_l, &
                           wk_bfft_l,lsize,ibsize,isrsize,fft_l_size,iblock,itot,frestart)
!                !!$write(nfout,*) 'debug loop itot,idavid=',itot,idavid
!                if(ik==1.and.iprimddavidson>= 2) &
!                  & call m_ES_wd_zaj_small_portion(nfout,ik," -- after md davidson --",21)
!                if(frestart) exit David_Loop
!                if(eigenvalues_are_converged(n_unconv)) exit Loop
             end do Block_Loop
write(nfout,*) 'A',feigconv
!             call m_ES_betar_dot_WFs_4_each_k(nfout,ik)   ! -> fsr_l,fsi_l
          end do Loop
          if(iprimddavidson>=2) then
             write(nfout,'("MdDavidson: ik=",i5," itot=",i5," subspace=",i5)') &
                                                                 ik, itot_now, nsize_sb_now
          end if
          call deallocate_t_matrix
! === DEBUG by tkato 2012/06/14 ================================================
          deallocate(wk_bfft_l)
          deallocate(bfft_l)
! ==============================================================================
       enddo      ! k-point loop
    enddo      ! spin loop
    call deallocate_matrix
    if(iprimddavidson>=2) then
       write(nfout,'("Modified Davidson: max_itot=",i5)') max_itot
    end if
!
! ========================================================================================
! === NOTE: m_ES_sort_eigen_values_3D causes difference with ORG_Parallel!!! =============
! ========================================================================================
!   call m_ES_sort_eigen_values_3D()
! ========================================================================================
! ========================================================================================
! ========================================================================================
!!!  ( in case of af=1 )
!    if(af /= 0) then
!       call cp_eigen_values_for_af       !-(contained here)
!       call expand_neordr_and_nrvf_ordr  !-(contained here)
!    end if

    call get_ipri0(iprimddavidson,ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko_3D(nfout,mode=SCF)
!
    deallocate(bfft);   deallocate(afft)
    call m_ES_dealloc_scss_etc()
    deallocate(ekin_l,p_l)
!deallocate(vnldi,hdiag,sdiag)

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

    logical function eigenvalues_are_converged(n_unconv)
       integer, intent(out) :: n_unconv
       integer :: ib

       n_unconv = 0
       eigenvalues_are_converged = .true.
       do ib=1,np_e
          if(.not.feigconv(ib)) then
             eigenvalues_are_converged = .false.
             n_unconv = n_unconv + 1
          end if
       end do
    end function eigenvalues_are_converged

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

  end subroutine m_ESmddavid_Renew_WF


  subroutine allocate_matrix
    integer:: i,j

    nblock = npartition_mddavid
    if(np_e .lt. npartition_mddavid) nblock = np_e

    allocate(nsize_subspace(nblock))
    allocate(nsize_matrix(nblock))
    allocate(ista_e_l(nblock))
    allocate(iend_e_l(nblock))
    allocate(ielm_e_l(nblock))
    ielm_e_l=np_e/nblock
    j = mod(np_e,nblock)
    do i = 1, j
       ielm_e_l(i) = ielm_e_l(i) + 1
    end do
    ista_e_l(1) = 1
    do i = 2, nblock
       ista_e_l(i) = ista_e_l(i-1) + ielm_e_l(i-1)
       iend_e_l(i-1) = ista_e_l(i) - 1
    end do
    iend_e_l(nblock) = np_e
    do i=1,nblock
      nsize_subspace(i)=ielm_e_l(i)*2
      nsize_matrix(i)=nsize_subspace(i)*(nsize_subspace(i)+1)/2
    end do
    msize_subspace=maxval(nsize_subspace)
    msize_matrix=maxval(nsize_matrix)
    allocate(feigconv(np_e))
    allocate(ibover(msize_subspace,nblock))
    allocate(fsr(np_e,nlmta,2))
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       allocate(fsi(np_e,nlmta,2))
    end if
    eps_residual = eps_residual_mddavid
!!    allocate(zajold_l(kg1,np_e,kimg,ndavid)) ! MPI
!!    allocate(zaj_l_backup(kg1,np_e,kimg)) ! MPI
  end subroutine allocate_matrix

! =========================== added by K. Tagami ============== 11.0
  subroutine allocate_matrix_noncl
    integer:: i,j

    nblock = npartition_mddavid
    if(np_e .lt. npartition_mddavid) nblock = np_e

    allocate(nsize_subspace(nblock))
    allocate(nsize_matrix(nblock))
    allocate(ista_e_l(nblock))
    allocate(iend_e_l(nblock))
    allocate(ielm_e_l(nblock))
    ielm_e_l=np_e/nblock
    j = mod(np_e,nblock)
    do i = 1, j
       ielm_e_l(i) = ielm_e_l(i) + 1
    end do
    ista_e_l(1) = 1
    do i = 2, nblock
       ista_e_l(i) = ista_e_l(i-1) + ielm_e_l(i-1)
       iend_e_l(i-1) = ista_e_l(i) - 1
    end do
    iend_e_l(nblock) = np_e
    do i=1,nblock
      nsize_subspace(i)=ielm_e_l(i)*2
      nsize_matrix(i)=nsize_subspace(i)*(nsize_subspace(i)+1)/2
    end do
    msize_subspace=maxval(nsize_subspace)
    msize_matrix=maxval(nsize_matrix)
    allocate(feigconv(np_e))
    allocate(ibover(msize_subspace,nblock))

    allocate( fsr_noncl(np_e,nlmta,2,ndim_spinor ) )
    if(.not.(kv3/ndim_spinor == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       allocate(fsi_noncl(np_e,nlmta,2,ndim_spinor))
    end if
    eps_residual = eps_residual_mddavid
!!    allocate(zajold_l(kg1,np_e,kimg,ndavid)) ! MPI
!!    allocate(zaj_l_backup(kg1,np_e,kimg)) ! MPI
  end subroutine allocate_matrix_noncl
! ======================================================== 11.0

  subroutine deallocate_matrix
    deallocate(nsize_subspace)
    deallocate(nsize_matrix)
    deallocate(ista_e_l)
    deallocate(iend_e_l)
    deallocate(ielm_e_l)
    deallocate(feigconv)
    deallocate(ibover)
    deallocate(fsr)
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       deallocate(fsi)
    end if
!!    deallocate(zajold_l)
!!    deallocate(zaj_l_backup)
  end subroutine deallocate_matrix

! =========================== addd by K. Tagami ============= 11.0
  subroutine deallocate_matrix_noncl
    deallocate(nsize_subspace)
    deallocate(nsize_matrix)
    deallocate(ista_e_l)
    deallocate(iend_e_l)
    deallocate(ielm_e_l)
    deallocate(feigconv)
    deallocate(ibover)
    deallocate(fsr_noncl)
    if ( allocated(fsi_noncl) )  deallocate( fsi_noncl )
!!    deallocate(zajold_l)
!!    deallocate(zaj_l_backup)
  end subroutine deallocate_matrix_noncl
! ============================================================ 11.0


  subroutine allocate_t_matrix_3D(ik)
    integer, intent(in) :: ik
    integer :: kimg_t
    if(k_symmetry(ik) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = kimg
    end if
    allocate(zat_l(maxval(np_g1k),np_e,kimg,2)) ! MPI
    allocate(zah_l(maxval(np_g1k),np_e,kimg,2)) ! MPI
    allocate(w1hw2(msize_matrix*kimg_t))
    allocate(w1sw2(msize_matrix*kimg_t))
    if(npes>1) then
       allocate(w1hw2_mpi(msize_matrix*kimg_t))
       allocate(w1sw2_mpi(msize_matrix*kimg_t))
    end if
    allocate(wfsd_l(maxval(np_g1k),np_e,ik:ik,kimg)); wfsd_l = 0.d0
    allocate(bsdr_l(np_e,nlmta,ik:ik)); bsdr_l = 0.d0
    allocate(bsdi_l(np_e,nlmta,ik:ik)); bsdi_l = 0.d0
  end subroutine allocate_t_matrix_3D

  subroutine deallocate_t_matrix
! ======================================== modified by K. Tgagami ======= 11.0
!    deallocate(zat_l) ! MPI
!    deallocate(zah_l) ! MPI

    if ( noncol ) then
       deallocate(zat_l_noncl) ! MPI
       deallocate(zah_l_noncl) ! MPI
    else
       deallocate(zat_l) ! MPI
       deallocate(zah_l) ! MPI
    endif
! ======================================================================== 11.0
    deallocate(w1hw2)
    deallocate(w1sw2)
    if(npes>1) then
       deallocate(w1hw2_mpi)
       deallocate(w1sw2_mpi)
    end if
    deallocate(wfsd_l)
    deallocate(bsdr_l)
    deallocate(bsdi_l)
  end subroutine deallocate_t_matrix

  subroutine decide_correction_vector_3D(precon,ik,ekin_l,afft_l,bfft_l, &
                                         wk_bfft_l,lsize,ibsize,isrsize,fft_l_size,p_l)
    integer, intent(in)       :: precon, ik, lsize, ibsize, isrsize, fft_l_size
    real(kind=DP), intent(in)  :: ekin_l(maxval(np_g1k))
    real(kind=DP), intent(in)  :: afft_l(lsize*kimg)
    real(kind=DP), intent(out) :: bfft_l(lsize*kimg,1)
    real(kind=DP), intent(inout) :: wk_bfft_l(lsize*kimg,ibsize)
    real(kind=DP)              :: p_l(maxval(np_g1k))

    integer :: ib
    integer :: id_sname = -1
    call tstatc0_begin('decide_correction_vector ', id_sname,1)

    do ib = 1, np_e ! MPI
#ifdef __TIMER_COMM__
       call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l,0)
#else
       call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l)
#endif
       call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_y(myrank_g))
       call m_FFT_Direct_3D(nfout,wk_bfft_l,lsize,ibsize)
       call map_fft_to_WF_3D(ik,lsize,ibsize,wk_bfft_l,bfft_l,isrsize,fft_l_size)
       call SD_direction_3D(precon,ik,ib,ekin_l,bfft_l,p_l,lsize) !-here
    end do

    call orthogonalize_SD_drctns_3D(ik,to=OTHER_BANDS)  ! -(m_ES_WF_by_SDorCG) ->(wfsd_l, bsd(ri)_l)
!    call orthogonalize_SD_drctns(ik,to=ALL_BANDS)  ! -(m_ES_WF_by_SDorCG) ->(wfsd_l, bsd(ri)_l)

    call normalize_wfsd_3D(ik)

    call tstatc0_end(id_sname)
  end subroutine decide_correction_vector_3D

  subroutine SD_direction_3D(precon,ik,ibo,ekin_l,VlocalW,p_l,lsize)
    integer     , intent(in)                   :: precon,ik,ibo,lsize
    real(kind=DP), intent(in), dimension(maxval(np_g1k))  :: ekin_l
    real(kind=DP), intent(in), dimension(lsize*kimg,1) :: VlocalW
    real(kind=DP)             , dimension(maxval(np_g1k))  :: p_l

    integer       :: i, ib, iadd
    real(kind=DP) :: devr,denom, e1, devi, norm

!   ib = map_z(ibo)                                  ! MPI
    ib = ibo
    denom = 1.d0/product(fft_box_size_WF(1:3,1))
    norm = 0.d0

    if(kimg == 1) then
       do i = ista_g1k(ik), iend_g1k(ik)
          iadd = i - ista_g1k(ik) + 1
          devr  = (ekin_l(iadd)-eko_l(ib,ik))*zaj_l(iadd,ib,ik,1)&
               & + VlocalW(iadd,1)*denom + vnlph_l(iadd,ib,1)
          wfsd_l(iadd,ib,ik,1) = - devr
          norm = norm + devr*devr
       end do
    else if(kimg == 2) then
       do i = ista_g1k(ik), iend_g1k(ik)
          iadd = i - ista_g1k(ik) + 1
          e1    = ekin_l(iadd) - eko_l(ib,ik)
          devr  = e1*zaj_l(iadd,ib,ik,1) + VlocalW(2*iadd-1,1)*denom+vnlph_l(iadd,ib,1)
          devi  = e1*zaj_l(iadd,ib,ik,2) + VlocalW(2*iadd,  1)*denom+vnlph_l(iadd,ib,2)
          wfsd_l(iadd,ib,ik,1) = - devr
          wfsd_l(iadd,ib,ik,2) = - devi
          norm = norm + devr*devr + devi*devi
       end do
       if(k_symmetry(ik) == GAMMA) then
          if(ista_g1k(ik) == 1) then
             devr=wfsd_l(1,ib,ik,1)
             devi=wfsd_l(1,ib,ik,2)
             norm = norm*2.d0 - devr*devr - devi*devi
          end if
       end if
    end if
    call mpi_allreduce(MPI_IN_PLACE,norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)

    feigconv(ib)=.false.
    if(sqrt(norm) .lt. eps_residual) feigconv(ib)=.true.

    if(precon==ON) then
      call decide_precon_factor_wfsd_3D(ik,ibo,ekin_l,p_l)
!call decide_precon_factor_david(ik,hdiag,sdiag,eko_l(ib,ik),p)
      if(kimg == 1) then
         do i = ista_g1k(ik), iend_g1k(ik)
            iadd = i - ista_g1k(ik) + 1
            wfsd_l(iadd,ib,ik,1) = p_l(iadd)*wfsd_l(iadd,ib,ik,1)
         end do
      else if(kimg == 2) then
         do i = ista_g1k(ik), iend_g1k(ik)
            iadd = i - ista_g1k(ik) + 1
            wfsd_l(iadd,ib,ik,1) = p_l(iadd)*wfsd_l(iadd,ib,ik,1)
            wfsd_l(iadd,ib,ik,2) = p_l(iadd)*wfsd_l(iadd,ib,ik,2)
         end do
      end if
    end if
  end subroutine SD_direction_3D

  subroutine decide_precon_factor_wfsd_3D(ik,ibo,ekin_l,p_l)
    integer, intent(in)                         :: ik,ibo
    real(kind=DP), intent(in),  dimension(maxval(np_g1k))  :: ekin_l
    real(kind=DP), intent(out), dimension(maxval(np_g1k))  :: p_l

    integer       :: i, iadd
    real(kind=DP) :: ektot, x, x1, x2, d_ektot

!    call kinetic_energy_wfsd(ik,ibo,ekin,ektot)   ! -here
    call kinetic_energy_3D(ik,ibo,ekin_l,ektot)   ! -here
    d_ektot = 4.d0/ektot/3.d0
    p_l = 0.d0
    do i = ista_g1k(ik), iend_g1k(ik)
       iadd = i - ista_g1k(ik) + 1
       x = ekin_l(iadd)*d_ektot
       x1 = (x*x+9.d0)*(x+3.d0)
       x2 = (x*x)*(x*x)
       p_l(iadd)  = x1/(x1 + x2 )
    end do
!    p=p*d_ektot
  end subroutine decide_precon_factor_wfsd_3D

  subroutine kinetic_energy_wfsd(ik,ibo,dekin,ektot)
    integer, intent(in) :: ik, ibo
    real(kind=DP), intent(in), dimension(kg1)  :: dekin
    real(kind=DP), intent(out)                 :: ektot
    integer  :: i, ri, ib
    ektot = 0.d0
    ib=map_z(ibo)
    do ri = 1, kimg
       do i = 1, iba(ik)
          ektot = ektot + dekin(i)*wfsd_l(i,ib,ik,ri)**2   ! MPI
       end do
    end do

    if(k_symmetry(ik) == GAMMA) ektot = ektot*2.d0

  end subroutine kinetic_energy_wfsd

  subroutine kinetic_energy_3D(ik,ibo,dekin,ektot)
    integer, intent(in) :: ik, ibo
    real(kind=DP), intent(in), dimension(maxval(np_g1k)) :: dekin
    real(kind=DP), intent(out)                 :: ektot
    integer  :: i, ib, iadd
    ektot = 0.d0
!   ib = map_z(ibo)
    ib = ibo
    if(kimg == 1) then
       do i = ista_g1k(ik), iend_g1k(ik)
          iadd = i - ista_g1k(ik) + 1
          ektot = ektot + dekin(iadd)*zaj_l(iadd,ib,ik,1)**2
       end do
    else
       do i = ista_g1k(ik), iend_g1k(ik)
          iadd = i - ista_g1k(ik) + 1
          ektot = ektot + dekin(iadd)*( zaj_l(iadd,ib,ik,1)**2 &
               &                      + zaj_l(iadd,ib,ik,2)**2)
       end do
    end if
    call mpi_allreduce(MPI_IN_PLACE,ektot,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
    if(k_symmetry(ik) == GAMMA)  ektot = ektot*2.d0
  end subroutine kinetic_energy_3D

  subroutine orthogonalize_SD_drctns_3D(ik,to)
    integer, intent(in) :: ik,to

    integer :: itmp
    integer :: id_sname = -1
    call tstatc0_begin('orthogonalize_SD_drctns in Modified Davidson ', id_sname)

!    if(modnrm == EXECUT) call m_ES_betar_dot_Psi_4_each_k(wfsd_l,ik,ik,ik,bsdr_l,bsdi_l)
    !                                                          ->bsd(ri)_l
    itmp=modnrm
    modnrm=EXECUT
    call m_ES_betar_dot_Psi_4_each_k_3D(wfsd_l,ik,ik,ik,bsdr_l,bsdi_l)
    call m_ES_orthogonalize_SD_to_WFs_3D(ik,to,wfsd_l,bsdr_l,bsdi_l) ! ->(wfsd_l,bsd(ri)_l)
    modnrm=itmp
!    call m_ES_betar_dot_Psi_4_each_k(wfsd_l,ik,ik,ik,bsdr_l,bsdi_l)
    call tstatc0_end(id_sname)
  end subroutine orthogonalize_SD_drctns_3D


  subroutine normalize_wfsd_3D(ik)
    integer,intent(in) :: ik
    real(kind=DP) :: norm, wfsdr, wfsdi
    integer :: ib1,ii,ib,iadd

    do ib1 = 1, np_e
!      ib=map_z(ib1)
       ib=ib1
       norm = 0.d0
       if(kimg==1) then
          do ii=ista_g1k(ik),iend_g1k(ik)
             iadd = ii - ista_g1k(ik) + 1
             wfsdr = wfsd_l(iadd,ib,ik,kimg)
             norm = norm + wfsdr*wfsdr
          end do
          call mpi_allreduce(MPI_IN_PLACE,norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
          norm = 1.d0/sqrt(norm)
          do ii=ista_g1k(ik),iend_g1k(ik)
             iadd = ii - ista_g1k(ik) + 1
             wfsd_l(iadd,ib,ik,1) = wfsd_l(iadd,ib,ik,1)*norm
          end do
          bsdr_l(ib,1:nlmta,ik) = bsdr_l(ib,1:nlmta,ik)*norm
          bsdi_l(ib,1:nlmta,ik) = bsdi_l(ib,1:nlmta,ik)*norm
       else
          do ii=ista_g1k(ik),iend_g1k(ik)
             iadd = ii - ista_g1k(ik) + 1
             wfsdr = wfsd_l(iadd,ib,ik,1   )
             wfsdi = wfsd_l(iadd,ib,ik,kimg)
             norm = norm + wfsdr*wfsdr+wfsdi*wfsdi
          end do
          call mpi_allreduce(MPI_IN_PLACE,norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
          norm = 1.d0/sqrt(norm)
          do ii=ista_g1k(ik),iend_g1k(ik)
             iadd = ii - ista_g1k(ik) + 1
             wfsd_l(iadd,ib,ik,1) = wfsd_l(iadd,ib,ik,1)*norm
             wfsd_l(iadd,ib,ik,2) = wfsd_l(iadd,ib,ik,2)*norm
          enddo
          if(k_symmetry(ik) == GAMMA) then
            bsdr_l(ib,1:nlmta,ik) = bsdr_l(ib,1:nlmta,ik)*norm
          else
            bsdr_l(ib,1:nlmta,ik) = bsdr_l(ib,1:nlmta,ik)*norm
            bsdi_l(ib,1:nlmta,ik) = bsdi_l(ib,1:nlmta,ik)*norm
          end if
       end if
    end do

  end subroutine normalize_wfsd_3D

  subroutine prepare_Hloc_phi_3D(ik,ekin_l,afft_l,bfft_l, &
                                 wk_bfft_l,lsize,ibsize,isrsize,fft_l_size,itot)
    integer, intent(in) :: ik,lsize,ibsize,isrsize,fft_l_size
    integer, intent(in) :: itot
    real(kind=DP), intent(in)  :: ekin_l(maxval(np_g1k))
    real(kind=DP), intent(in)  :: afft_l(lsize*kimg)
    real(kind=DP), intent(out) :: bfft_l(lsize*kimg,1)
    real(kind=DP), intent(inout) :: wk_bfft_l(lsize*kimg,ibsize)

    integer       :: ib1,i1,ii,ib,iadd
    real(kind=DP), allocatable, dimension(:,:) :: fs_mpi,fs_mpi2
    integer :: is, is1
    real(kind=DP) :: denom
    real(kind=DP) :: dr1,dr2,di1,di2,dd
    integer :: id_sname = -1, ipri0
    call tstatc0_begin('prepare_Hloc_phi (mddavidson) ', id_sname,1)

    call get_ipri0(iprimddavidson,ipri0)
    allocate(fs_mpi(np_e,nlmta))
    allocate(fs_mpi2(np_e,nlmta))

    denom = 1.d0/product(fft_box_size_WF(1:3,1))

! (zaj_l <- (T+Vloc)|phi> )
!!    zaj_l(:,:,ik,:) = zajold_l(:,:,:,idavid)
!( tenchi ) (zat_l <- zaj_l)

    zat_l(:,:,:,1) = zaj_l(:,:,ik,:)
    zat_l(:,:,:,2) = wfsd_l(:,:,ik,:)
    fs_mpi=0.d0
    do is = 1, np_fs ! MPI
       is1=nis_fs(myrank_g)+is-1
       do ib = 1, np_e ! MPI
          fs_mpi(ib,is1) = fsr_l(ib,is,ik)
       end do
    end do
    fs_mpi2=0.d0
    call mpi_allreduce(fs_mpi,fs_mpi2,np_e*nlmta &
      & ,mpi_double_precision,mpi_sum &
      & ,mpi_ke_world,ierr)       ! MPI
    fsr(1:np_e,1:nlmta,1) = fs_mpi2(1:np_e,1:nlmta)
    fsr(:,:,2)=bsdr_l(:,:,ik)
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       fs_mpi=0.d0
       do is = 1, np_fs ! MPI
          is1=nis_fs(myrank_g)+is-1
          do ib = 1, np_e ! MPI
             fs_mpi(ib,is1) = fsi_l(ib,is,ik)
          end do
       end do
       fs_mpi2=0.d0
       call mpi_allreduce(fs_mpi,fs_mpi2,np_e*nlmta &
         & ,mpi_double_precision,mpi_sum &
         & ,mpi_ke_world,ierr)       ! MPI
       fsi(1:np_e,1:nlmta,1) = fs_mpi2(1:np_e,1:nlmta)
       fsi(:,:,2)=bsdi_l(:,:,ik)
    end if
    deallocate(fs_mpi)
    deallocate(fs_mpi2)

    if(itot == 1) then
      do ib1 = 1, np_e ! MPI
         ib = ib1
#ifdef __TIMER_COMM__
         call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l,0)
#else
         call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l)
#endif
         call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_y(myrank_g))
         call m_FFT_Direct_3D(nfout,wk_bfft_l,lsize,ibsize)
         call map_fft_to_WF_3D(ik,lsize,ibsize,wk_bfft_l,bfft_l,isrsize,fft_l_size)
         if(kimg == 1) then
            do ii=ista_g1k(ik),iend_g1k(ik)
               iadd = ii - ista_g1k(ik) + 1
               dr1 = zaj_l(iadd,ib,ik,1)
               dr2 = bfft_l(iadd,1)*denom
               zah_l(iadd,ib,1,1) = ekin_l(iadd)*dr1+dr2
            enddo
         else
            do ii=ista_g1k(ik),iend_g1k(ik)
               iadd = ii - ista_g1k(ik) + 1
               dr1  = zaj_l(iadd,ib,ik,1)
               di1  = zaj_l(iadd,ib,ik,kimg)
               zah_l(iadd,ib,1,   1) = ekin_l(iadd)*dr1+bfft_l(2*iadd-1,1)*denom
               zah_l(iadd,ib,kimg,1) = ekin_l(iadd)*di1+bfft_l(2*iadd,  1)*denom
            enddo
         endif
      enddo
    end if

!!( tenchi ) (zah_l <- zaj_l)
    zaj_l(:,:,ik,:) = wfsd_l(:,:,ik,:)

    do ib1 = 1, np_e ! MPI
       ib = ib1
#ifdef __TIMER_COMM__
       call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l,0)
#else
       call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l)
#endif
       call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_y(myrank_g))
       call m_FFT_Direct_3D(nfout,wk_bfft_l,lsize,ibsize)
       call map_fft_to_WF_3D(ik,lsize,ibsize,wk_bfft_l,bfft_l,isrsize,fft_l_size)
       if(kimg == 1) then
          do ii=ista_g1k(ik),iend_g1k(ik)
             iadd = ii - ista_g1k(ik) + 1
             dr1 = zaj_l(iadd,ib,ik,1)
             dr2 = bfft_l(iadd,1)*denom
             zah_l(iadd,ib,1,2) = ekin_l(iadd)*dr1+dr2
          enddo
       else
          do ii=ista_g1k(ik),iend_g1k(ik)
             iadd = ii - ista_g1k(ik) + 1
             dr1  = zaj_l(iadd,ib,ik,1)
             di1  = zaj_l(iadd,ib,ik,kimg)
             zah_l(iadd,ib,1,   2) = ekin_l(iadd)*dr1+bfft_l(2*iadd-1,1)*denom
             zah_l(iadd,ib,kimg,2) = ekin_l(iadd)*di1+bfft_l(2*iadd,  1)*denom
          enddo
       endif
    enddo

    call tstatc0_end(id_sname)

  contains

    subroutine get_ipri0(ipri_in, ipri_out)
      integer, intent(in)  :: ipri_in
      integer, intent(out) :: ipri_out
      if(npes > 1) then
         if(mype == 0) ipri_out = ipri_in
         call mpi_bcast(ipri_out,1,mpi_integer,0,mpi_k_world(myrank_k),ierr)
      else
         ipri_out = ipri_in
      end if
    end subroutine get_ipri0

  end subroutine prepare_Hloc_phi_3D


  subroutine evolve_WFs_in_subspace_3D &
                                (ik,ispin,ekin_l,afft_l,bfft_l, &
                                 wk_bfft_l,lsize,ibsize,isrsize,fft_l_size,iblock,itot,frestart)
    integer, intent(in) :: ik,ispin,lsize,ibsize,isrsize,fft_l_size
    integer, intent(in) :: iblock,itot
    logical, intent(out) :: frestart
    real(kind=DP), intent(in)  :: ekin_l(maxval(np_g1k))
    real(kind=DP), intent(in)  :: afft_l(lsize*kimg)
    real(kind=DP), intent(out) :: bfft_l(lsize*kimg,1)
    real(kind=DP), intent(inout) :: wk_bfft_l(lsize*kimg,ibsize)
    integer :: iadd
! (allocatable variables)
    real(kind=DP), allocatable,dimension(:) ::     eig
    real(kind=DP), allocatable,dimension(:,:) ::   vec
    real(kind=DP), allocatable,dimension(:) ::     eko_d
    real(kind=DP), allocatable,dimension(:) ::     eko_d_mpi
    integer, allocatable,dimension(:) ::     occup
    integer, allocatable,dimension(:) ::     occup_mpi

    integer       :: ib1,ib2,ib1to,ib2to,i1,ii,ri,ib
    integer       :: ibb1,ibb2
    integer       :: ii1,ii2,iter,iter1,iter2
    real(kind=DP) :: eko1, eko2, ekod
    real(kind=DP) :: hr2,hi2,dr1,dr2,di1,di2,dd
    integer :: ip0,ip0b,ip1,ip1b,ib1n,ib2n,ndata,nshift,kimg_t,ig1
    integer :: noffset
    integer :: nsize_max_sb_now
    integer :: ierr_diag
    integer :: nel,nsta,nend
    integer :: id_sname = -1, ipri0
    call tstatc0_begin('evolve_WFs_in_subspace (modified davidson) ', id_sname,1)

    call get_ipri0(iprimddavidson,ipri0)

! ==================== added by K. Tagami ========== 11.0
#ifdef forsafe
    ekod = 0.0d0
#endif
! ================================================== 11.0

    if(k_symmetry(ik) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = kimg
    end if

    nel =ielm_e_l(iblock)
    nsta=ista_e_l(iblock)
    nend=iend_e_l(iblock)

!    ip0=nsize_subspace(iblock)
!    do ib=1,ip0
!      ibover(ib,iblock) = ib
!    end do

    ip0=nel
    do ib=1,nel
      ibover(ib,iblock)=ib
      if(.not.feigconv(nsta+ib-1)) then
        ip0=ip0+1
        ibover(nel+ib,iblock) = ip0
      else
        ibover(nel+ib,iblock) = -1
      end if
    end do

    nsize_sb_now = ip0
    nsize_mt_now = nsize_sb_now*(nsize_sb_now+1)/2
!     nsize_max_sb_now = neg*idavid
!     nsize_sb_now = nsize_subspace(iblock)          ! For the moment
!     nsize_mt_now = nsize_matrix(iblock)
    nsize_max_sb_now = nsize_subspace(iblock)
    if(iprimddavidson >=2) then
       write(nfout,*) 'ibover=',ibover(1:nsize_max_sb_now,iblock)
    end if
!
    allocate(eig(nsize_sb_now)); eig=0.d0
    allocate(vec(nsize_sb_now*kimg_t,nsize_sb_now))
!    allocate(eko_d(neg));     eko_d = 0.d0
!    allocate(eko_d_mpi(neg))
!    allocate(occup(neg)); occup=0
!    allocate(occup_mpi(neg))

    if(iprimddavidson >=2) then
       write(nfout,*) 'Modified Davidson:ik,iblock,nsize_sb_now=', ik,iblock,nsize_sb_now
    end if

!    do ib1 = 1, neg
!       if(map_e(ib1) == myrank_e) eko_d(ib1) = eko_l(map_z(ib1),ik)  ! MPI
!    end do
!    call mpi_allreduce(eko_d,eko_d_mpi,neg,mpi_double_precision,mpi_sum &
!         & ,mpi_k_world(myrank_k),ierr)       ! MPI
!    eko_d = eko_d_mpi                         ! MPI
    eko1 = sum(eko_l(nsta:nend,ik))
!
!    do ib1=1,neg
!       if(map_e(ib1) == myrank_e) then !MPI
!          if( occup_l(map_z(ib1),ik) > 0.d0 ) occup(ib1) = 1  ! MPI
!       end if
!    end do
!    call mpi_allreduce(occup,occup_mpi,neg,mpi_integer,mpi_sum &
!         & ,mpi_k_world(myrank_k),ierr)       ! MPI
!    occup = occup_mpi                         ! MPI
!
!! (zaj_l <- (T+Vloc)|phi> )
!!!    zaj_l(:,:,ik,:) = zajold_l(:,:,:,idavid)
!!( tenchi ) (zat_l <- zaj_l)
!    call m_ES_W_transpose(ista_k,iend_k,ik,zaj_l,zat_l(1,1,1,idavid))
!    do ib1 = ista_e, iend_e, istep_e     ! MPI
!       call m_ES_WF_in_Rspace(ik,ib1,bfft)!(swffft)
!       call m_FFT_Vlocal_W(afft,bfft)    ! (afft, bfft)-> (bfft)
!       call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
!       ib = map_z(ib1)                 ! MPI
!       if(kimg == 1) then
!          do ii=1,iba(ik)
!             i1  = igf(nbase(ii,ik))
!             dr1 = zaj_l(ii,ib,ik,1)
!             dr2 = bfft(i1)*denom
!             zaj_l(ii,ib,ik,1)= ekin(ii)*dr1+dr2
!          enddo
!       else
!          do ii=1,iba(ik)
!             i1  = igf(nbase(ii,ik))
!             dr1  = zaj_l(ii,ib,ik,1)
!             di1  = zaj_l(ii,ib,ik,kimg)
!             zaj_l(ii,ib,ik,1)= ekin(ii)*dr1+bfft(2*i1-1)*denom
!             zaj_l(ii,ib,ik,kimg)= ekin(ii)*di1+bfft(2*i1)*denom
!          enddo
!       endif
!    enddo
!!( tenchi ) (zah_l <- zaj_l)
!    call m_ES_W_transpose(ista_k,iend_k,ik,zaj_l,zah_l(1,1,1))

!! (make matrix elements )
!    ! <n|T+Vloc|m> !
    do ibb2 = 1,nsize_max_sb_now
       if(ibover(ibb2,iblock)<0) cycle
       ib2 = ibover(ibb2,iblock)
       iter2 = (ibb2-1)/nel+1
       ii2  = ibb2-nel*(iter2-1)
       ii2  = nsta+ii2-1
       ip0b = ib2*(ib2-1)/2
       do ibb1 = 1,ibb2
          if(ibover(ibb1,iblock)<0) cycle
          ib1 = ibover(ibb1,iblock)
          iter1 = (ibb1-1)/nel+1
          ii1 = ibb1-nel*(iter1-1)
          ii1 = nsta+ii1-1
          ip0 = ip0b + ib1
          if(kimg == 1) then
             w1hw2(ip0) = 0.d0
             w1sw2(ip0) = 0.d0
             do ii = ista_g1k(ik),iend_g1k(ik) ! MPI
                iadd = ii - ista_g1k(ik) + 1
                hr2 = zah_l(iadd,ii2,1,iter2)
                dr2 = zat_l(iadd,ii2,1,iter2)
                dr1 = zat_l(iadd,ii1,1,iter1)
                w1hw2(ip0) = w1hw2(ip0) + dr1*hr2
                w1sw2(ip0) = w1sw2(ip0) + dr1*dr2
             end do
             call mpi_allreduce(MPI_IN_PLACE,w1hw2(ip0),1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
             call mpi_allreduce(MPI_IN_PLACE,w1sw2(ip0),1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
          else
             if(k_symmetry(ik) == GAMMA) then
                w1hw2(ip0) = 0.d0
                w1sw2(ip0) = 0.d0
                do ii = max(2,ista_g1k(ik)), iend_g1k(ik) ! MPI
                   iadd = ii - ista_g1k(ik) + 1
                   hr2 = zah_l(iadd,ii2,1,iter2) ! MPI
                   hi2 = zah_l(iadd,ii2,2,iter2) ! MPI
                   dr2 = zat_l(iadd,ii2,1,iter2) ! MPI
                   di2 = zat_l(iadd,ii2,2,iter2) ! MPI
                   dr1 = zat_l(iadd,ii1,1,iter1) ! MPI
                   di1 = zat_l(iadd,ii1,2,iter1) ! MPI
                   w1hw2(ip0) =w1hw2(ip0)+(dr1*hr2+di1*hi2)*2.d0
                   w1sw2(ip0) =w1sw2(ip0)+(dr1*dr2+di1*di2)*2.d0
                end do
                if(ista_g1k(ik) == 1) then
                   hr2 = zah_l(1,ii2,1,iter2) ! MPI
                   hi2 = zah_l(1,ii2,2,iter2) ! MPI
                   dr2 = zat_l(1,ii2,1,iter2) ! MPI
                   di2 = zat_l(1,ii2,2,iter2) ! MPI
                   dr1 = zat_l(1,ii1,1,iter1) ! MPI
                   di1 = zat_l(1,ii1,2,iter1) ! MPI
                   w1hw2(ip0) =w1hw2(ip0)+dr1*hr2+di1*hi2
                   w1sw2(ip0) =w1sw2(ip0)+dr1*dr2+di1*di2
                end if
                call mpi_allreduce(MPI_IN_PLACE,w1hw2(ip0),1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
                call mpi_allreduce(MPI_IN_PLACE,w1sw2(ip0),1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
             else
                w1hw2(2*ip0-1:2*ip0) = 0.d0
                w1sw2(2*ip0-1:2*ip0) = 0.d0
                do ii = ista_g1k(ik), iend_g1k(ik) ! MPI
                   iadd = ii - ista_g1k(ik) + 1
                   hr2 = zah_l(iadd,ii2,1,iter2) ! MPI
                   hi2 = zah_l(iadd,ii2,2,iter2) ! MPI
                   dr2 = zat_l(iadd,ii2,1,iter2) ! MPI
                   di2 = zat_l(iadd,ii2,2,iter2) ! MPI
                   dr1 = zat_l(iadd,ii1,1,iter1) ! MPI
                   di1 = zat_l(iadd,ii1,2,iter1) ! MPI
                   w1hw2(2*ip0-1) =w1hw2(2*ip0-1)+dr1*hr2+di1*hi2
                   w1hw2(2*ip0  ) =w1hw2(2*ip0  )+dr1*hi2-di1*hr2
                   w1sw2(2*ip0-1) =w1sw2(2*ip0-1)+dr1*dr2+di1*di2
                   w1sw2(2*ip0  ) =w1sw2(2*ip0  )+dr1*di2-di1*dr2
                end do
                call mpi_allreduce(MPI_IN_PLACE,w1hw2(2*ip0-1),2,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
                call mpi_allreduce(MPI_IN_PLACE,w1sw2(2*ip0-1),2,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
             end if
          end if
       end do
    end do
    if(iprimddavidson >= 3) call wd_w1hw2(" -- w1hw2 without nl part--",iblock)
    ! <n|Vnl|m>
    call add_nonlocal_part ! w1hw2 = w1hw2 + w1Vnlw2
                           ! w1sw2 = w1sw2 + w1qw2
    if(iprimddavidson >= 3) call wd_w1hw2(" -- w1hw2 with nl part--",iblock)
!
    if(iprimddavidson >= 2) call wd_w1hw2(" -- just after making w1hw2 --",iblock)
!    if(ipridavidson >= 2) then
!       write(nfout,*) 'neordr for ik = ',ik
!       write(nfout,9002) (neordr(ib1,ik),ib1=1,neg)
!       write(nfout,*) 'nrvf_ordr for ik = ',ik
!       write(nfout,9002) (nrvf_ordr(ib1,ik),ib1=1,neg)
!       write(nfout,*) 'eig'
!       write(nfout,'(5x,10f8.4)') (eko_l(ib1,ik),ib1=1,np_e)
!    endif
!9002 format(5x,10i8)

!! (Diagonalization )  !!

    if(kimg_t == 1) then
       call dspgvx_driver_loc(eig,vec,w1hw2,w1sw2,ierr_diag,nel)
    else
       call zhpgvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag,nel)
    endif

    frestart = .false.
    if(ierr_diag /= 0) then
!       zaj_l(:,:,ik,:) = zaj_l_backup(:,:,:)
!       do ib1 = 1, neg
!          if(map_e(ib1) == myrank_e) then         ! MPI
!             eko_l(map_z(ib1),ik)=eko_d(ib1)
!          end if
!       end do
      frestart = .true.
      zaj_l(:,nsta:nend,ik,:)=zat_l(:,nsta:nend,:,1)
      if(iprimddavidson >= 2) then
        write(nfout,*) '** restart Modified Davidson iteration **'
        write(nfout,*) 'ik=',ik,' iblock=',iblock,' itot=',itot
      end if
!!$print *,'Restart'
      goto 9000
    else

      feigconv(nsta:nend) = .false.
      do ib=1,nel
         if(occup_l(nsta+ib-1,ik) > 0.d0) then
            if(abs(eko_l(nsta+ib-1,ik)-eig(ib)) < delta_eig_occup_md) &
                                                   feigconv(nsta+ib-1) = .true.
         else
            if(abs(eko_l(nsta+ib-1,ik)-eig(ib)) < delta_eig_empty_md) &
                                                   feigconv(nsta+ib-1) = .true.
         end if
      end do
       if(ipri0 >= 2) then
          write(nfout,*) 'eko_l for ik = ',ik
          write(nfout,*) 'iblock       = ',iblock
          write(nfout,9001) (eko_l(nsta+ib-1,ik),ib=1,nel)
          write(nfout,*) 'eig for ik = ',ik
          write(nfout,9001) (eig(ib),ib=1,nel)
          call wd_w1hw2(" -- after diagonalization --",iblock)
!sum eko
          dr1=0.d0;dr2=0.d0
          do ib1=1,nel
             dr1=dr1+eko_l(nsta+ib1-1,ik) ! MPI
             dr2=dr2+eig(ib1)
          enddo
          write(nfout,'(" sum of eko_l, eig, abs diff =",3e25.10)') dr1,dr2,abs(dr2-dr1)
       endif
!!! (subspace rotation) !!
       call subspace_rotation ! vec,zat_l -> zat_l
!!( tenchi ) (zaj_l <- zat_l)
!       iter = min(idavid+1,ndavid)
!       call m_ES_W_transpose_back(ista_k,iend_k,ik,zaj_l,zat_l(1,1,1,iter))
!!       zaj_l_backup(:,nsta:nend,:) = zaj_l(:,nsta:nend,ik,:)
!! (eko_l)
       do ib1 = 1, nel
             eko_l(nsta+ib1-1,ik)=eig(ib1)
       end do
       if(iprimddavidson >= 2) then
          eko2 = sum(eig(1:nel))
          write(nfout,1201) ik,eko1,ekod,eko2

          write(nfout,*) 'eko_l'
          write(nfout,9001) (eko_l(nsta+ib1-1,ik),ib1=1,nel)
       endif
1201 format(' %% for ik = ',i4,4x,' eko1&ekod&eko2 = ',3f14.7)
9001 format(5x,6f12.5)
!!! (neordr & nrvf_ordr)
!
    end if
9000 continue
!    neordr(1:neg,ik) = (/(ib1,ib1=1,neg)/)
!    nrvf_ordr(1:neg,ik) = (/(ib1,ib1=1,neg)/)
!
!! (deallocate)
!    deallocate(eko_d)
!    deallocate(eko_d_mpi)
    deallocate(eig)
    deallocate(vec)
!    deallocate(occup)
!    deallocate(occup_mpi)
!
    call tstatc0_end(id_sname)
!
  contains
    subroutine get_ipri0(ipri_in, ipri_out)
      integer, intent(in)  :: ipri_in
      integer, intent(out) :: ipri_out
      if(npes > 1) then
         if(mype == 0) ipri_out = ipri_in
         call mpi_bcast(ipri_out,1,mpi_integer,0,mpi_k_world(myrank_k),ierr)
      else
         ipri_out = ipri_in
      end if
    end subroutine get_ipri0
!
    subroutine wd_w1hw2(somecomment,iblock)
      character(len=*), intent(in) :: somecomment
      integer,intent(in) :: iblock
      integer :: ib1, ib2, nel_wd, nsb_wd
      write(nfout,'(a35)') somecomment
      write(nfout,*) 'w1hw2 for ik = ',ik
      nel_wd = 8
      nsb_wd = 8
      if(nel_wd > nel) nel_wd = nel
      if(nsb_wd > nsize_sb_now) nsb_wd = nsize_sb_now
      if(kimg_t==1) then
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1hw2(ip0),ip0=ip0b+1,ip0b+ib2)
         end do
      else
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1hw2(2*ip0-1),ip0=ip0b+1,ip0b+ib2)
         end do
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1hw2(2*ip0),ip0=ip0b+1,ip0b+ib2)
         end do
      end if
      write(nfout,*) 'w1sw2 for ik = ',ik
      if(kimg_t==1) then
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1sw2(ip0),ip0=ip0b+1,ip0b+ib2)
         end do
      else
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1sw2(2*ip0-1),ip0=ip0b+1,ip0b+ib2)
         end do
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1sw2(2*ip0),ip0=ip0b+1,ip0b+ib2)
         end do
      end if
9001  format(5x,9f12.5)
      write(nfout,*) 'eko_l for ik = ',ik
      write(nfout,9001) (eko_l(nsta+ib1-1,ik),ib1=1,nel_wd)
    end subroutine wd_w1hw2

    subroutine add_nonlocal_part
      integer :: ip,ib1,ib2,ibb1,ibb2
      integer       :: ia, lmt1, lmt2, it, p, s, ib
      real(kind=DP) :: facv,facq,vr,vi,qr,qi
      real(kind=DP) :: tmpr,tmpi
! ========================== added by K. Tagami ========== 11.0
#ifdef forsafe
      integer :: ipaw_tmp
#endif
! ======================================================== 11.0
      do ibb2 = 1,nsize_max_sb_now
         if(ibover(ibb2,iblock)<0) cycle
         ib2 = ibover(ibb2,iblock)
         iter2= (ibb2-1)/nel+1
         ii2  = ibb2-nel*(iter2-1)
         ii2  = nsta+ii2-1
         ip0b = ib2*(ib2-1)/2
         do ibb1 = 1,ibb2
            if(ibover(ibb1,iblock)<0) cycle
            ib1 = ibover(ibb1,iblock)
            iter1=(ibb1-1)/nel+1
            ii1  = ibb1-nel*(iter1-1)
            ii1  = nsta+ii1-1
            ip0 = ip0b + ib1
!            if(mod(ip0-1,nrank_e)/=myrank_e) cycle
            if(kimg_t==1) then
               vr=0.d0
               qr=0.d0
            else
               vr=0.d0
               vi=0.d0
               qr=0.d0
               qi=0.d0
            end if
            do ia = 1, natm
               it = ityp(ia)
! ========================== added by K. Tagami =================== 11.0
#ifdef forsafe
               ipaw_tmp = ipaw(it)
#endif
! ================================================================= 11.0
               do lmt1 = 1, ilmt(it)
                  p = lmta(lmt1,ia)
                  do lmt2 = 1, ilmt(it)
                     s = lmta(lmt2,ia)
! ========================== modified by K. Tagami =================== 11.0
!                     if(ipaw(it)==0)then
#ifdef forsafe
                     if ( ipaw_tmp == 0 ) then
#else
                     if (ipaw(it)==0 ) then
#endif
! ================================================================= 11.0
                        facv   = iwei(ia)*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin))
                     else
                        facv   = iwei(ia)*(dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin))
                     endif
                     facq   = iwei(ia)*q(lmt1,lmt2,it)
                     if(kimg==1) then
                        tmpr = fsr(ii1,p,iter1)*fsr(ii2,s,iter2)&
                    &        + fsi(ii1,p,iter1)*fsi(ii2,s,iter2)
                        vr = vr + facv*tmpr
                        qr = qr + facq*tmpr
                     else
                        if(k_symmetry(ik) == GAMMA) then
                           tmpr = fsr(ii1,p,iter1)*fsr(ii2,s,iter2)
                           vr = vr + facv*tmpr
                           qr = qr + facq*tmpr
                        else
                           tmpr = fsr(ii1,p,iter1)*fsr(ii2,s,iter2)&
                    &        + fsi(ii1,p,iter1)*fsi(ii2,s,iter2)
                           tmpi = fsr(ii1,p,iter1)*fsi(ii2,s,iter2)&
                    &        - fsi(ii1,p,iter1)*fsr(ii2,s,iter2)
                           vr = vr + facv*tmpr
                           vi = vi + facv*tmpi
                           qr = qr + facq*tmpr
                           qi = qi + facq*tmpi
                        end if
                     end if
                  end do
               end do
            end do
            if(kimg_t==1) then
               w1hw2(ip0) = w1hw2(ip0) + vr
               w1sw2(ip0) = w1sw2(ip0) + qr
            else
               w1hw2(2*ip0-1) = w1hw2(2*ip0-1) + vr
               w1hw2(2*ip0  ) = w1hw2(2*ip0  ) + vi
               w1sw2(2*ip0-1) = w1sw2(2*ip0-1) + qr
               w1sw2(2*ip0  ) = w1sw2(2*ip0  ) + qi
            end if
         end do
      end do
    end subroutine add_nonlocal_part

    subroutine subspace_rotation
      integer :: ib1,ib2,ibb2,iadd,is,is1
!!$      real(kind=DP), dimension(np_g1k_x,neg,kimg) :: zaj_wk
      real(kind=DP), allocatable, dimension(:,:,:) :: zaj_wk
      real(kind=DP), allocatable, dimension(:,:,:) :: zah_wk
      real(kind=DP), allocatable, dimension(:,:)   :: fsr_wk
      real(kind=DP), allocatable, dimension(:,:)   :: fsi_wk
      allocate(zaj_wk(maxval(np_g1k),nel,kimg))
      if(itot /= max_iter_mddavid) then
        allocate(zah_wk(maxval(np_g1k),nel,kimg))
        zah_wk(:,:,:) = 0.d0
      end if
      allocate(fsr_wk(nel,nlmta))
      if(k_symmetry(ik) /= GAMMA) then
        allocate(fsi_wk(nel,nlmta))
        fsi_wk(:,:)=0.d0
      end if

      zaj_wk(:,:,:) = 0.d0
      fsr_wk(:,:)=0.d0
      if(kimg==1) then
         do ib1=1,nel
            do ibb2=1,nsize_max_sb_now
               if(ibover(ibb2,iblock)<0) cycle
               ib2 = ibover(ibb2,iblock)
               iter2=(ibb2-1)/nel+1
               ii2=ibb2-nel*(iter2-1)
               ii2=nsta+ii2-1
               hr2=vec(ib2,ib1)
               do ii=ista_g1k(ik),iend_g1k(ik)
                  iadd = ii - ista_g1k(ik) + 1
                  zaj_wk(iadd,ib1,kimg) = zaj_wk(iadd,ib1,kimg) + zat_l(iadd,ii2,kimg,iter2)*hr2
               end do
               do ii=1,nlmta
                 fsr_wk(ib1,ii) = fsr_wk(ib1,ii) + fsr(ii2,ii,iter2)*hr2
                 fsi_wk(ib1,ii) = fsi_wk(ib1,ii) + fsi(ii2,ii,iter2)*hr2
               end do
               if(itot /= max_iter_mddavid) then
                 do ii=ista_g1k(ik),iend_g1k(ik)
                  iadd = ii - ista_g1k(ik) + 1
                    zah_wk(iadd,ib1,kimg) = zah_wk(iadd,ib1,kimg) + zah_l(iadd,ii2,kimg,iter2)*hr2
                 end do
               end if
            end do
         end do
      else
         if(k_symmetry(ik) == GAMMA) then
            do ib1=1,nel
               do ibb2=1,nsize_max_sb_now
                  if(ibover(ibb2,iblock)<0) cycle
                  ib2 = ibover(ibb2,iblock)
                  iter2=(ibb2-1)/nel+1
                  ii2=ibb2-nel*(iter2-1)
                  ii2=nsta+ii2-1
                  hr2=vec(ib2,ib1)
                  do ii=ista_g1k(ik),iend_g1k(ik)
                     iadd = ii - ista_g1k(ik) + 1
                     dr1=zat_l(iadd,ii2,1   ,iter2)
                     di1=zat_l(iadd,ii2,kimg,iter2)
                     zaj_wk(iadd,ib1,1   ) = zaj_wk(iadd,ib1,1   ) + dr1*hr2
                     zaj_wk(iadd,ib1,kimg) = zaj_wk(iadd,ib1,kimg) + di1*hr2
                  end do
                  do ii=1,nlmta
                    fsr_wk(ib1,ii) = fsr_wk(ib1,ii) + fsr(ii2,ii,iter2)*hr2
                  end do
                  if(itot /= max_iter_mddavid) then
                    do ii=ista_g1k(ik),iend_g1k(ik)
                       iadd = ii - ista_g1k(ik) + 1
                       dr1=zah_l(iadd,ii2,1   ,iter2)
                       di1=zah_l(iadd,ii2,kimg,iter2)
                       zah_wk(iadd,ib1,1   ) = zah_wk(iadd,ib1,1   ) + dr1*hr2
                       zah_wk(iadd,ib1,kimg) = zah_wk(iadd,ib1,kimg) + di1*hr2
                    end do
                  end if
               end do
            end do
         else
            do ib1=1,nel
               do ibb2=1,nsize_max_sb_now
                  if(ibover(ibb2,iblock)<0) cycle
                  ib2 = ibover(ibb2,iblock)
                  iter2=(ibb2-1)/nel+1
                  ii2=ibb2-nel*(iter2-1)
                  ii2=nsta+ii2-1
                  hr2=vec(2*ib2-1,ib1)
                  hi2=vec(2*ib2  ,ib1)
                  do ii=ista_g1k(ik),iend_g1k(ik)
                     iadd = ii - ista_g1k(ik) + 1
                     dr1=zat_l(iadd,ii2,1   ,iter2)
                     di1=zat_l(iadd,ii2,kimg,iter2)
                     zaj_wk(iadd,ib1,1   ) = zaj_wk(iadd,ib1,1   ) + dr1*hr2 - di1*hi2
                     zaj_wk(iadd,ib1,kimg) = zaj_wk(iadd,ib1,kimg) + dr1*hi2 + di1*hr2
                  end do
                  do ii=1,nlmta
                     dr1=fsr(ii2,ii,iter2)
                     di1=fsi(ii2,ii,iter2)
                     fsr_wk(ib1,ii) = fsr_wk(ib1,ii) + dr1*hr2 - di1*hi2
                     fsi_wk(ib1,ii) = fsi_wk(ib1,ii) + dr1*hi2 + di1*hr2
                  end do
                  if(itot /= max_iter_mddavid) then
                    do ii=ista_g1k(ik),iend_g1k(ik)
                       iadd = ii - ista_g1k(ik) + 1
                       dr1=zah_l(iadd,ii2,1   ,iter2)
                       di1=zah_l(iadd,ii2,kimg,iter2)
                       zah_wk(iadd,ib1,1   ) = zah_wk(iadd,ib1,1   ) + dr1*hr2 - di1*hi2
                       zah_wk(iadd,ib1,kimg) = zah_wk(iadd,ib1,kimg) + dr1*hi2 + di1*hr2
                    end do
                  end if
               end do
            end do
         end if
      end if
!print *,itot,  itot /= max_iter_mddavid
      zaj_l(:,nsta:nend,ik,:) = zaj_wk(:,:,:)
      if(itot /= max_iter_mddavid) zah_l(:,nsta:nend,:,1) = zah_wk(:,:,:)
!     fsr_l(nsta:nend,:,ik) = fsr_wk(:,:)
      do is = 1, np_fs ! MPI
         is1=nis_fs(myrank_g)+is-1
         do ib = nsta, nend ! MPI
             fsr_l(ib,is,ik) = fsr_wk(ib-nsta+1,is1)
         end do
      end do
!     if(k_symmetry(ik) /= GAMMA) fsi_l(nsta:nend,:,ik) = fsi_wk(:,:)
      if(k_symmetry(ik) /= GAMMA) then
         do is = 1, np_fs ! MPI
            is1=nis_fs(myrank_g)+is-1
            do ib = nsta, nend ! MPI
               fsi_l(ib,is,ik) = fsi_wk(ib-nsta+1,is1)
            end do
         end do
      end if

      deallocate(zaj_wk,fsr_wk)
      if(itot /= max_iter_mddavid) deallocate(zah_wk)
      if(k_symmetry(ik) /= GAMMA) deallocate(fsi_wk)
    end subroutine subspace_rotation

  end subroutine evolve_WFs_in_subspace_3D





  subroutine m_ESmddavid_Subspace_Rotation(nfout)
    integer, intent(in) :: nfout

    integer             :: ispin, ik, iksnl, switch_of_eko_part
    real(kind=DP), allocatable, dimension(:) ::  afft, bfft
    real(kind=DP), allocatable, dimension(:) :: ekin
    real(kind=DP), allocatable, dimension(:) :: afft_l
    real(kind=DP), allocatable, dimension(:,:) :: wk_bfft_l
    real(kind=DP), allocatable, dimension(:,:) :: bfft_l
    integer :: lsize, ibsize, isrsize, fft_l_size
    real(kind=DP), allocatable, dimension(:) :: ekin_l
    integer :: ipri0
#ifdef FFT_USE_SSL2_PAD
    lsize = max(nel_fft_x(myrank_g),nel_fft_y(myrank_g),nel_fft_z(myrank_g))
#else
    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
#endif
    allocate(afft_l(lsize*kimg), stat=ierr)
    if(ierr /= 0) then
       write(nfout,*)' m_ESmddavid_Subspace_Rotation: Not allocated afft_l array'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 201, ierr)
    endif
    ibsize = 1

    allocate(ekin_l(maxval(np_g1k)))
    call m_ES_alloc_scss_etc_3D()
    allocate(afft(nfft)); allocate(bfft(nfft))
    call allocate_fsri_3D

    do ispin = 1, nspin, (af+1)
       call m_ES_Vlocal_in_Rspace_3D(ispin,afft_l,lsize,1,OFF)      ! (ptfft1) vlhxc_l->afft
       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) /= myrank_k) cycle          ! MPI
          isrsize = min(lsize,mp_g1k(ik))
          fft_l_size  = nel_fft_x(myrank_g)
          allocate(wk_bfft_l(lsize*kimg,ibsize) ,stat=ierr)
          allocate(bfft_l(lsize*kimg,ibsize) ,stat=ierr)
          if (ierr /= 0) then
             write(nfout,*)' m_ESmddavid_Subspace_Rotation:  Not allocate '
             call flush(nfout)
             call mpi_abort(mpi_comm_world, 205, ierr)
          endif
          iksnl = (ik-1)/nspin + 1

          call allocate_t_matrix_sr_3D(ik) ! -> np_g1k_x
          call m_pwBS_kinetic_energies_3D(ik,vkxyz,ekin_l) ! (diakin) ->ekin
          call m_ES_Vnonlocal_W_3D(ik,iksnl,ispin,switch_of_eko_part=OFF) ! -> vnlph_l
          call allreduce_fs_sr_3D(ik) ! -> fsr,fsi
          call evolve_WFs_in_subspace_sr_3D& !-(m_ES_WF_by_ModifiedDavidson)
                                      &(ik,ispin,ekin_l,afft_l,bfft_l, &
                                        wk_bfft_l,lsize,ibsize,isrsize,fft_l_size) !-> zaj_l
          if(ik==1.and.iprimddavidson>= 2) &
            & call m_ES_wd_zaj_small_portion_3D(nfout,ik," -- after davidson subspace rotation --",21)
          call m_ES_betar_dot_WFs_4_each_k_3D(nfout,ik)   ! -> fsr_l,fsi_l
          if(iprimddavidson>=2) then
             write(nfout,'("Davidson Subspace Rotation: ik=",i5," subspace=",i5)') ik, nsize_sb_now
          end if
          call deallocate_t_matrix_sr
          deallocate(wk_bfft_l)
          deallocate(bfft_l)

       enddo      ! k-point loop
    enddo      ! spin loop

    call deallocate_fsri

!!  ( in case of af=1 )
    if(af /= 0) then
       call cp_eigen_values_for_af       !-(contained here)
       call expand_neordr_and_nrvf_ordr  !-(contained here)
    end if

    call get_ipri0(iprimddavidson,ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko_3D(nfout,mode=SCF)

    deallocate(bfft);   deallocate(afft)
    call m_ES_dealloc_scss_etc()
    deallocate(ekin_l)

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

  end subroutine m_ESmddavid_Subspace_Rotation


  subroutine allocate_fsri_3D
    nsize_sb_now = neg
    nsize_mt_now =  nsize_sb_now*(nsize_sb_now+1)/2
    allocate(fsr(neg,nlmta,1))
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       allocate(fsi(neg,nlmta,1))
    end if
    allocate(zaj_l_backup(maxval(np_g1k),np_e,kimg)) ! MPI
  end subroutine allocate_fsri_3D

  subroutine deallocate_fsri
    deallocate(fsr)
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       deallocate(fsi)
    end if
    deallocate(zaj_l_backup)
  end subroutine deallocate_fsri

! =============================== added by K. Tagami ============= 11.0
  subroutine deallocate_fsri_noncl
    deallocate(fsr_noncl)
    if(.not.(kv3/ndim_spinor == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       deallocate(fsi_noncl)
    end if
    deallocate(zaj_l_backup_noncl)
  end subroutine deallocate_fsri_noncl
! ================================================================= 11.0

  subroutine allocate_t_matrix_sr_3D(ik)
    integer, intent(in) :: ik
    integer :: kimg_t
    if(k_symmetry(ik) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = kimg
    end if
#ifdef _ODD_BOUNDARY_
    if(mod(np_g1k(ik),2) == 0) then
       np_g1k_x = np_g1k(ik) + 1
    else
       np_g1k_x = np_g1k(ik)
    end if
#else
    np_g1k_x = np_g1k(ik)
#endif
    allocate(zat_l(maxval(np_g1k),neg,kimg,1)) ! MPI
    allocate(zah_l(maxval(np_g1k),np_e,kimg,1)) ! MPI
    allocate(w1hw2(nsize_mt_now*kimg_t))
    allocate(w1sw2(nsize_mt_now*kimg_t))
    if(npes>1) then
       allocate(w1hw2_mpi(nsize_mt_now*kimg_t))
       allocate(w1sw2_mpi(nsize_mt_now*kimg_t))
    end if
    zaj_l_backup(:,:,:) = zaj_l(:,:,ik,:)
  end subroutine allocate_t_matrix_sr_3D

  subroutine deallocate_t_matrix_sr
! ============================ modified by K. Tagami ========== 11.0
!    deallocate(zat_l) ! MPI
!    deallocate(zah_l) ! MPI

    if ( noncol ) then
       deallocate(zat_l_noncl) ! MPI
       deallocate(zah_l_noncl) ! MPI
    else
       deallocate(zat_l) ! MPI
       deallocate(zah_l) ! MPI
    endif
! ============================================================== 11.0
    deallocate(w1hw2)
    deallocate(w1sw2)
    if(npes>1) then
       deallocate(w1hw2_mpi)
       deallocate(w1sw2_mpi)
    end if
  end subroutine deallocate_t_matrix_sr

  subroutine allreduce_fs_sr_3D(ik)
    integer, intent(in) :: ik
    integer :: ib,ib1,kimg_t,is,is1
    real(kind=DP), allocatable, dimension(:,:) :: fs_mpi,fs_mpi2

    allocate(fs_mpi(neg,nlmta))
    allocate(fs_mpi2(neg,nlmta))

    fs_mpi=0.d0
    do is = 1, np_fs ! MPI
       is1=nis_fs(myrank_g)+is-1
       do ib = 1, np_e ! MPI
          ib1=neg_g(ib)
          fs_mpi(ib1,is1) = fsr_l(ib,is,ik)
       end do
    end do
    fs_mpi2=0.d0
    call mpi_allreduce(fs_mpi,fs_mpi2,neg*nlmta &
      & ,mpi_double_precision,mpi_sum &
      & ,mpi_k_world(myrank_k),ierr)       ! MPI
    fsr(1:neg,1:nlmta,1) = fs_mpi2(1:neg,1:nlmta)
    if(.not. k_symmetry(ik) == GAMMA) then
       fs_mpi=0.d0
       do is = 1, np_fs ! MPI
          is1=nis_fs(myrank_g)+is-1
          do ib = 1, np_e ! MPI
             ib1=neg_g(ib)
             fs_mpi(ib1,is1) = fsi_l(ib,is,ik)
          end do
       end do
       fs_mpi2=0.d0
       call mpi_allreduce(fs_mpi,fs_mpi2,neg*nlmta &
         & ,mpi_double_precision,mpi_sum &
         & ,mpi_k_world(myrank_k),ierr)       ! MPI
       fsi(1:neg,1:nlmta,1) = fs_mpi2(1:neg,1:nlmta)
    end if
    deallocate(fs_mpi)
    deallocate(fs_mpi2)
  end subroutine allreduce_fs_sr_3D

  subroutine evolve_WFs_in_subspace_sr_3D(ik,ispin,ekin_l,afft_l,bfft_l, &
                                          wk_bfft_l,lsize,ibsize,isrsize,fft_l_size)
    integer, intent(in) :: ik,ispin,lsize,ibsize,isrsize,fft_l_size
    real(kind=DP), intent(in)  :: ekin_l(maxval(np_g1k))
    real(kind=DP), intent(in)  :: afft_l(lsize*kimg)
    real(kind=DP), intent(out) :: bfft_l(lsize*kimg,1)
    real(kind=DP), intent(inout) :: wk_bfft_l(lsize*kimg,ibsize)
    integer :: iadd, ib2_
! (allocatable variables)
    real(kind=DP), allocatable,dimension(:) ::     eig
    real(kind=DP), allocatable,dimension(:,:) ::   vec
    real(kind=DP), allocatable,dimension(:) ::     eko_d
    real(kind=DP), allocatable,dimension(:) ::     eko_d_mpi
    integer, allocatable,dimension(:) ::     occup
    integer, allocatable,dimension(:) ::     occup_mpi

    integer       :: ib1,ib2,ib1to,ib2to,i1,ii,ri,ib
!    integer       :: ibb1,ibb2
!    integer       :: ii1,ii2,iter,iter1,iter2
    real(kind=DP) :: denom, eko1, eko2, ekod
    real(kind=DP) :: hr2,hi2,dr1,dr2,di1,di2,dd
    integer :: ip0,ip0b,ip1,ip1b,ib1n,ib2n,ndata,nshift,kimg_t,ig1
    integer :: noffset
!    integer :: nsize_max_sb_now
    integer :: ierr_diag
    integer :: id_sname = -1, ipri0
    call tstatc0_begin('evolve_WFs_in_subspace_sr_3D(davidson) ', id_sname,1)

    call get_ipri0(iprimddavidson,ipri0)

    denom = 1.d0/product(fft_box_size_WF(1:3,1))
    if(k_symmetry(ik) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = kimg
    end if

!    nsize_sb_now = nsize_subspace(1)
!    nsize_mt_now = nsize_matrix(1)

    allocate(eig(nsize_sb_now)); eig=0.d0
    allocate(vec(nsize_sb_now*kimg_t,nsize_sb_now))
    allocate(eko_d(neg));     eko_d = 0.d0
    allocate(eko_d_mpi(neg))
    allocate(occup(neg)); occup=0
    allocate(occup_mpi(neg))

    if(iprimddavidson >=2) then
       write(nfout,*) 'MdDavidson Subspace Rotation:ik,nsize_sb_now=', ik,nsize_sb_now
    end if


    do ib1 = 1, np_e
       eko_d(neg_g(ib1)) = eko_l(ib1,ik)  ! MPI
    end do
    call mpi_allreduce(eko_d,eko_d_mpi,neg,mpi_double_precision,mpi_sum &
         & ,mpi_kg_world,ierr) ! MPI
    eko_d = eko_d_mpi                         ! MPI
    eko1 = sum(eko_d(1:neg))

    do ib1=1,neg
       if(map_e(ib1) == myrank_e) then !MPI
          if( occup_l(map_z(ib1),ik) > 0.d0 ) occup(ib1) = 1  ! MPI
       end if
    end do
    call mpi_allreduce(occup,occup_mpi,neg,mpi_integer,mpi_sum &
         & ,mpi_k_world(myrank_k),ierr)       ! MPI
    occup = occup_mpi                         ! MPI

! (zaj_l <- (T+Vloc)|phi> )
!!    zaj_l(:,:,ik,:) = zajold_l(:,:,:,idavid)
!( tenchi ) (zat_l <- zaj_l)
    zat_l(:,:,:,1) = 0.0d0
    do ib1 = 1, np_e
       zat_l(:,neg_g(ib1),:,1) = zaj_l(:,ib1,ik,:)
    enddo
    call mpi_allreduce(MPI_IN_PLACE,zat_l(:,:,:,1),maxval(np_g1k)*neg*kimg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
    do ib1 = 1, np_e ! MPI
       ib = ib1 ! MPI
#ifdef __TIMER_COMM__
       call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l,0)
#else
       call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l)
#endif
       call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_y(myrank_g))
       call m_FFT_Direct_3D(nfout,wk_bfft_l,lsize,ibsize)
       call map_fft_to_WF_3D(ik,lsize,ibsize,wk_bfft_l,bfft_l,isrsize,fft_l_size)
       if(kimg == 1) then
          do ii=ista_g1k(ik),iend_g1k(ik)
             iadd = ii - ista_g1k(ik) + 1
             dr1 = zaj_l(iadd,ib,ik,1)
             dr2 = bfft_l(iadd,1)*denom
             zaj_l(iadd,ib,ik,1)= ekin_l(iadd)*dr1+dr2
          enddo
       else
          do ii=ista_g1k(ik),iend_g1k(ik)
             iadd = ii - ista_g1k(ik) + 1
             dr1  = zaj_l(iadd,ib,ik,1)
             di1  = zaj_l(iadd,ib,ik,kimg)
             zaj_l(iadd,ib,ik,1)    = ekin_l(iadd)*dr1+bfft_l(2*iadd-1,1)*denom
             zaj_l(iadd,ib,ik,kimg) = ekin_l(iadd)*di1+bfft_l(2*iadd,  1)*denom
          enddo
       endif
    enddo
!( tenchi ) (zah_l <- zaj_l)
    zah_l(:,:,:,1) = zaj_l(:,:,ik,:)

! (make matrix elements )
! parallel loop
    ! <n|T+Vloc|m> G-wise parallel
    w1hw2 = 0.d0
    w1sw2 = 0.d0
!   do ib2 = 1,nsize_sb_now
    do ib2_ = 1, np_e
       ib2 = neg_g(ib2_)
       ip0b = ib2*(ib2-1)/2
       do ib1 = 1,ib2
          ip0 = ip0b + ib1
          if(kimg == 1) then
             do ii = ista_g1k(ik), iend_g1k(ik) ! MPI
                iadd = ii - ista_g1k(ik) + 1
                hr2 = zah_l(iadd,ib2_,1,1)
                dr2 = zat_l(iadd,ib2, 1,1)
                dr1 = zat_l(iadd,ib1, 1,1)
                w1hw2(ip0) = w1hw2(ip0) + dr1*hr2
                w1sw2(ip0) = w1sw2(ip0) + dr1*dr2
             end do
          else
             if(k_symmetry(ik) == GAMMA) then
                do ii = max(ista_g1k(ik),2), iend_g1k(ik) ! MPI
                   iadd = ii - ista_g1k(ik) + 1
                   hr2 = zah_l(iadd,ib2_,1,1) ! MPI
                   hi2 = zah_l(iadd,ib2_,2,1) ! MPI
                   dr2 = zat_l(iadd,ib2, 1,1) ! MPI
                   di2 = zat_l(iadd,ib2, 2,1) ! MPI
                   dr1 = zat_l(iadd,ib1, 1,1) ! MPI
                   di1 = zat_l(iadd,ib1, 2,1) ! MPI
                   w1hw2(ip0) =w1hw2(ip0)+(dr1*hr2+di1*hi2)*2.d0
                   w1sw2(ip0) =w1sw2(ip0)+(dr1*dr2+di1*di2)*2.d0
                end do
                if(ista_g1k(ik) == 1) then
                   hr2 = zah_l(1,ib2_,1,1) ! MPI
                   hi2 = zah_l(1,ib2_,2,1) ! MPI
                   dr2 = zat_l(1,ib2, 1,1) ! MPI
                   di2 = zat_l(1,ib2, 2,1) ! MPI
                   dr1 = zat_l(1,ib1, 1,1) ! MPI
                   di1 = zat_l(1,ib1, 2,1) ! MPI
                   w1hw2(ip0) =w1hw2(ip0)+dr1*hr2+di1*hi2
                   w1sw2(ip0) =w1sw2(ip0)+dr1*dr2+di1*di2
                end if
             else
                do ii = ista_g1k(ik), iend_g1k(ik) ! MPI
                   iadd = ii - ista_g1k(ik) + 1
                   hr2 = zah_l(iadd,ib2_,1,1) ! MPI
                   hi2 = zah_l(iadd,ib2_,2,1) ! MPI
                   dr2 = zat_l(iadd,ib2, 1,1) ! MPI
                   di2 = zat_l(iadd,ib2, 2,1) ! MPI
                   dr1 = zat_l(iadd,ib1, 1,1) ! MPI
                   di1 = zat_l(iadd,ib1, 2,1) ! MPI
                   w1hw2(2*ip0-1) =w1hw2(2*ip0-1)+dr1*hr2+di1*hi2
                   w1hw2(2*ip0  ) =w1hw2(2*ip0  )+dr1*hi2-di1*hr2
                   w1sw2(2*ip0-1) =w1sw2(2*ip0-1)+dr1*dr2+di1*di2
                   w1sw2(2*ip0  ) =w1sw2(2*ip0  )+dr1*di2-di1*dr2
                end do
             end if
          end if
       end do
    end do
    if(iprimddavidson >= 3) call wd_w1hw2(" -- w1hw2 without nl part--")
    ! <n|Vnl|m> G-wise parallel
    if(myrank_g == 0) then
    call add_nonlocal_part ! w1hw2 = w1hw2 + w1Vnlw2
                           ! w1sw2 = w1sw2 + w1qw2
    endif
    if(iprimddavidson >= 3) call wd_w1hw2(" -- w1hw2 with nl part--")

!! (spread sum of w1hw2 and w1sw2)
    if(npes > 1) then
       w1hw2_mpi = 0.d0
       w1sw2_mpi = 0.d0
       nshift = 0
       ndata = nsize_mt_now*kimg_t
       call mpi_allreduce(w1hw2(nshift+1),w1hw2_mpi,ndata,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI
       w1hw2(nshift+1:nshift+ndata) = w1hw2_mpi(1:ndata) ! MPI
       call mpi_allreduce(w1sw2(nshift+1),w1sw2_mpi,ndata,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI
       w1sw2(nshift+1:nshift+ndata) = w1sw2_mpi(1:ndata) ! MPI
    end if

    if(iprimddavidson >= 2) call wd_w1hw2(" -- just after making w1hw2 --")
    if(iprimddavidson >= 2) then
       write(nfout,*) 'neordr for ik = ',ik
       write(nfout,9002) (neordr(ib1,ik),ib1=1,neg)
       write(nfout,*) 'nrvf_ordr for ik = ',ik
       write(nfout,9002) (nrvf_ordr(ib1,ik),ib1=1,neg)
       write(nfout,*) 'eig'
       write(nfout,'(5x,10f8.4)') (eko_l(ib1,ik),ib1=1,np_e)
    endif
9002 format(5x,10i8)

!! (Diagonalization )  !!

    if(kimg_t == 1) then
       call dspgvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
    else
       call zhpgvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag)
    endif

    if(ierr_diag /= 0) then
       zaj_l(:,:,ik,:) = zaj_l_backup(:,:,:)
       do ib1 = 1, np_e
          eko_l(ib1,ik)=eko_d(neg_g(ib1))
       end do
       if(iprimddavidson >= 2) then
          write(nfout,*) '** Mod Davidson SR error **'
          stop
       end if
    else

!!$       if(ipridavidson >= 2) then
       if(ipri0 >= 2) then
          write(nfout,*) 'eko_d for ik = ',ik
          write(nfout,9001) (eko_d(ib),ib=1,neg)
          write(nfout,*) 'eig for ik = ',ik
          write(nfout,9001) (eig(ib),ib=1,neg)
          call wd_w1hw2(" -- after diagonalization --")
!sum eko
          dr1=0.d0;dr2=0.d0
          do ib1=1,np_e
             ib1to = neordr(neg_g(ib1),ik)
             if(map_e(ib1to) == myrank_e) dr1=dr1+eko_l(ib1,ik) ! MPI
             dr2=dr2+eig(neg_g(ib1))
          enddo
          call mpi_allreduce(dr1,di1,1,mpi_double_precision,mpi_sum,mpi_kg_world,ierr) ! MPI
          dr1 = di1  ! MPI
          write(nfout,'(" sum of eko_l, eig, abs diff =",3e25.10)') dr1,dr2,abs(dr2-dr1)
       endif
!! (subspace rotation) !!
       call subspace_rotation ! vec,zat_l -> zat_l
!( tenchi ) (zaj_l <- zat_l)
       do ib1 = 1, np_e
          zaj_l(:,ib1,ik,:) = zat_l(:,neg_g(ib1),:,1)
       enddo
       zaj_l_backup(:,:,:) = zaj_l(:,:,ik,:)
!! (eko_l)
       do ib1 = 1, np_e
          eko_l(ib1,ik)=eig(neg_g(ib1))
       end do
       if(iprimddavidson >= 2) then
          eko2 = sum(eig(1:neg))
          write(nfout,1201) ik,eko1,ekod,eko2

          write(nfout,*) 'eko_l'
          write(nfout,9001) (eko_l(ib1,ik),ib1=1,np_e)
       endif
1201 format(' %% for ik = ',i4,4x,' eko1&ekod&eko2 = ',3f14.7)
9001 format(5x,6f12.5)
!! (neordr & nrvf_ordr)

    end if

    neordr(1:neg,ik) = (/(ib1,ib1=1,neg)/)
    nrvf_ordr(1:neg,ik) = (/(ib1,ib1=1,neg)/)

! (deallocate)
    deallocate(eko_d)
    deallocate(eko_d_mpi)
    deallocate(eig)
    deallocate(vec)
    deallocate(occup)
    deallocate(occup_mpi)

    call tstatc0_end(id_sname)

  contains
    subroutine get_ipri0(ipri_in, ipri_out)
      integer, intent(in)  :: ipri_in
      integer, intent(out) :: ipri_out
      if(npes > 1) then
         if(mype == 0) ipri_out = ipri_in
         call mpi_bcast(ipri_out,1,mpi_integer,0,mpi_k_world(myrank_k),ierr)
      else
         ipri_out = ipri_in
      end if
    end subroutine get_ipri0

    subroutine wd_w1hw2(somecomment)
      character(len=*), intent(in) :: somecomment
      integer :: ib1, ib2, neg_wd, nsb_wd
      write(nfout,'(a35)') somecomment
      write(nfout,*) 'w1hw2 for ik = ',ik
      neg_wd = 8
      nsb_wd = 8
      if(neg_wd > neg) neg_wd = neg
      if(nsb_wd > nsize_sb_now) nsb_wd = nsize_sb_now
      if(kimg_t==1) then
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1hw2(ip0),ip0=ip0b+1,ip0b+ib2)
         end do
      else
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1hw2(2*ip0-1),ip0=ip0b+1,ip0b+ib2)
         end do
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1hw2(2*ip0),ip0=ip0b+1,ip0b+ib2)
         end do
      end if
      write(nfout,*) 'w1sw2 for ik = ',ik
      if(kimg_t==1) then
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1sw2(ip0),ip0=ip0b+1,ip0b+ib2)
         end do
      else
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1sw2(2*ip0-1),ip0=ip0b+1,ip0b+ib2)
         end do
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1sw2(2*ip0),ip0=ip0b+1,ip0b+ib2)
         end do
      end if
9001  format(5x,9f12.5)
      write(nfout,*) 'eko_l for ik = ',ik
      write(nfout,9001) (eko_d(neordr(ib1,ik)),ib1=1,neg_wd)
    end subroutine wd_w1hw2

    subroutine add_nonlocal_part
      integer :: ip,ib1,ib2
      integer       :: ia, lmt1, lmt2, it, p, s, ib
      real(kind=DP) :: facv,facq,vr,vi,qr,qi
      real(kind=DP) :: tmpr,tmpi
! ========================== added by K. Tagami ========== 11.0
#ifdef forsafe
      integer :: ipaw_tmp
#endif
! ======================================================== 11.0
      do ib2 = 1,nsize_sb_now
         ip0b = ib2*(ib2-1)/2
         do ib1 = 1,ib2
            ip0 = ip0b + ib1
            if(mod(ip0-1,nrank_e)/=myrank_e) cycle
            if(kimg_t==1) then
               vr=0.d0
               qr=0.d0
            else
               vr=0.d0
               vi=0.d0
               qr=0.d0
               qi=0.d0
            end if
            do ia = 1, natm
               it = ityp(ia)
! ========================== added by K. Tagami =================== 11.0
#ifdef forsafe
               ipaw_tmp = ipaw(it)
#endif
! ================================================================= 11.0
               do lmt1 = 1, ilmt(it)
                  p = lmta(lmt1,ia)
                  do lmt2 = 1, ilmt(it)
                     s = lmta(lmt2,ia)
! ========================== modified by K. Tagami =================== 11.0
!                     if(ipaw(it)==0)then
#ifdef forsafe
                     if ( ipaw_tmp == 0 ) then
#else
                     if (ipaw(it)==0 ) then
#endif
! ================================================================= 11.0
                        facv   = iwei(ia)*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin))
                     else
                        facv   = iwei(ia)*(dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin))
                     endif
                     facq   = iwei(ia)*q(lmt1,lmt2,it)
                     if(kimg==1) then
                        tmpr = fsr(ib1,p,1)*fsr(ib2,s,1)&
                    &        + fsi(ib1,p,1)*fsi(ib2,s,1)
                        vr = vr + facv*tmpr
                        qr = qr + facq*tmpr
                     else
                        if(k_symmetry(ik) == GAMMA) then
                           tmpr = fsr(ib1,p,1)*fsr(ib2,s,1)
                           vr = vr + facv*tmpr
                           qr = qr + facq*tmpr
                        else
                           tmpr = fsr(ib1,p,1)*fsr(ib2,s,1)&
                    &        + fsi(ib1,p,1)*fsi(ib2,s,1)
                           tmpi = fsr(ib1,p,1)*fsi(ib2,s,1)&
                    &        - fsi(ib1,p,1)*fsr(ib2,s,1)
                           vr = vr + facv*tmpr
                           vi = vi + facv*tmpi
                           qr = qr + facq*tmpr
                           qi = qi + facq*tmpi
                        end if
                     end if
                  end do
               end do
            end do
            if(kimg_t==1) then
               w1hw2(ip0) = w1hw2(ip0) + vr
               w1sw2(ip0) = w1sw2(ip0) + qr
            else
               w1hw2(2*ip0-1) = w1hw2(2*ip0-1) + vr
               w1hw2(2*ip0  ) = w1hw2(2*ip0  ) + vi
               w1sw2(2*ip0-1) = w1sw2(2*ip0-1) + qr
               w1sw2(2*ip0  ) = w1sw2(2*ip0  ) + qi
            end if
         end do
      end do
    end subroutine add_nonlocal_part

    subroutine subspace_rotation
      integer :: ib1,ib2,ibb2,iadd
!!$      real(kind=DP), dimension(np_g1k_x,neg,kimg) :: zaj_wk
      real(kind=DP), allocatable, dimension(:,:,:) :: zaj_wk
!     allocate(zaj_wk(np_g1k_x,neg,kimg))
      allocate(zaj_wk(maxval(np_g1k),np_e,kimg))

      zaj_wk(:,:,:) = 0.d0
      if(kimg==1) then
         do ib1=1,np_e
            do ib2=1,neg
               do ii=ista_g1k(ik),iend_g1k(ik)
                  iadd = ii - ista_g1k(ik) + 1
                  zaj_wk(iadd,ib1,kimg) = zaj_wk(iadd,ib1,kimg) + zat_l(iadd,ib2,kimg,1)*vec(ib2,neg_g(ib1))
               end do
            end do
         end do
      else
         if(k_symmetry(ik) == GAMMA) then
            do ib1=1,np_e
               do ib2=1,neg
                  hr2=vec(ib2,neg_g(ib1))
                  do ii=ista_g1k(ik),iend_g1k(ik)
                     iadd = ii - ista_g1k(ik) + 1
                     dr1=zat_l(iadd,ib2,1   ,1)
                     di1=zat_l(iadd,ib2,kimg,1)
                     zaj_wk(iadd,ib1,1   ) = zaj_wk(iadd,ib1,1   ) + dr1*hr2
                     zaj_wk(iadd,ib1,kimg) = zaj_wk(iadd,ib1,kimg) + di1*hr2
                  end do
               end do
            end do
         else
            do ib1=1,np_e
               do ib2=1,neg
                  hr2=vec(2*ib2-1,neg_g(ib1))
                  hi2=vec(2*ib2  ,neg_g(ib1))
                  do ii=ista_g1k(ik),iend_g1k(ik)
                     iadd = ii - ista_g1k(ik) + 1
                     dr1=zat_l(iadd,ib2,1   ,1)
                     di1=zat_l(iadd,ib2,kimg,1)
                     zaj_wk(iadd,ib1,1   ) = zaj_wk(iadd,ib1,1   ) + dr1*hr2 - di1*hi2
                     zaj_wk(iadd,ib1,kimg) = zaj_wk(iadd,ib1,kimg) + dr1*hi2 + di1*hr2
                  end do
               end do
            end do
         end if
      end if
      zat_l(:,:,:,1) = 0.0d0
      do ib1 = 1, np_e
         zat_l(:,neg_g(ib1),:,1) = zaj_wk(:,ib1,:)
      enddo
      call mpi_allreduce(MPI_IN_PLACE,zat_l(:,:,:,1),maxval(np_g1k)*neg*kimg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
      deallocate(zaj_wk)
    end subroutine subspace_rotation

  end subroutine evolve_WFs_in_subspace_sr_3D


  subroutine dspgvx_driver(eig,vec,w1hw2,w1sw2,ierr,nel)
    real(kind=DP), intent(out) ,dimension(nsize_sb_now) :: eig
    real(kind=DP), intent(out) ,dimension(nsize_sb_now*neg) :: vec
    real(kind=DP), intent(inout) ,dimension(nsize_mt_now) :: w1hw2,w1sw2
    integer, intent(out) :: ierr
    integer, intent(in), optional :: nel

    integer :: ITYPE
    character(len=1) :: JOBZ,RANGE,UPLO
    integer :: il,iu
    real(kind=DP),allocatable,dimension(:) :: work_lapack
    integer, allocatable, dimension(:) :: iwork_lapack, ifail_lapack
    real(kind=DP) :: vl,vu,abstol
    integer :: info,m
    real(kind=DP), external :: dlamch
    !!$real(kind=DP), dimension(nsize_mt_now) :: ap,bp
    real(kind=DP), allocatable, dimension(:) :: ap,bp

    allocate(ap(nsize_mt_now),bp(nsize_mt_now))
    abstol = 2*dlamch('S')

    il=1; iu=neg
    if(present(nel)) iu=nel
!(LAPACK)  ITYPE = 1:  A*x = (lambda)*B*x, 2:  A*B*x = (lambda)*x, 3:  B*A*x = (lambda)*x
    ITYPE = 1
!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
    JOBZ = 'V'
!(LAPACK)  RANGE = A : all eigenvalues, V: all eigenvalues in (VL, VU], I: the IL-th through IU-th eigenvalues
    RANGE = 'I'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
    UPLO = 'U'
    allocate(work_lapack(8*nsize_sb_now))
    allocate(iwork_lapack(5*nsize_sb_now))
    allocate(ifail_lapack(nsize_sb_now))

    ap = w1hw2
    bp = w1sw2

    call dspgvx(ITYPE,JOBZ,RANGE,UPLO,nsize_sb_now,ap,bp &
    &          ,vl,vu,il,iu,abstol,m,eig,vec,nsize_sb_now &
    &          ,work_lapack,iwork_lapack,ifail_lapack,info)

    if(iprimddavidson >=2 .and. info/=0) then
       write(nfout,*) "debug(dspgvx) info=",info
       write(nfout,*) "debug(dspgvx) ifail"
       write(nfout,'(8(1x,i3))') ifail_lapack
       write(nfout,*) "debug(dspgvx) eig"
       write(nfout,'(8(1x,f10.5))') eig
    end if

    deallocate(work_lapack)
    deallocate(iwork_lapack)
    deallocate(ifail_lapack)

    deallocate(ap,bp)

    if(info/=0) then
       !!write(nfout,*) "dspgvx: info=",info
       !!stop 'error in dspgvx_driver'
       ierr = 1
    else
       ierr = 0
    end if

  end subroutine dspgvx_driver

  subroutine dspgvx_driver_loc(eig,vec,w1hw2,w1sw2,ierr,nel)
    integer, intent(in) :: nel
    real(kind=DP), intent(out) ,dimension(nsize_sb_now) :: eig
    real(kind=DP), intent(out) ,dimension(nsize_sb_now*nel) :: vec
    real(kind=DP), intent(inout) ,dimension(nsize_mt_now) :: w1hw2,w1sw2
    integer, intent(out) :: ierr

    integer :: ITYPE
    character(len=1) :: JOBZ,RANGE,UPLO
    integer :: il,iu
    real(kind=DP),allocatable,dimension(:) :: work_lapack
    integer, allocatable, dimension(:) :: iwork_lapack, ifail_lapack
    real(kind=DP) :: vl,vu,abstol
    integer :: info,m
    real(kind=DP), external :: dlamch
    !!$real(kind=DP), dimension(nsize_mt_now) :: ap,bp
    real(kind=DP), allocatable, dimension(:) :: ap,bp

    allocate(ap(nsize_mt_now),bp(nsize_mt_now))
    abstol = 2*dlamch('S')

    il=1; iu=nel
!(LAPACK)  ITYPE = 1:  A*x = (lambda)*B*x, 2:  A*B*x = (lambda)*x, 3:  B*A*x = (lambda)*x
    ITYPE = 1
!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
    JOBZ = 'V'
!(LAPACK)  RANGE = A : all eigenvalues, V: all eigenvalues in (VL, VU], I: the IL-th through IU-th eigenvalues
    RANGE = 'I'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
    UPLO = 'U'
    allocate(work_lapack(8*nsize_sb_now))
    allocate(iwork_lapack(5*nsize_sb_now))
    allocate(ifail_lapack(nsize_sb_now))

    ap = w1hw2
    bp = w1sw2

    call dspgvx(ITYPE,JOBZ,RANGE,UPLO,nsize_sb_now,ap,bp &
    &          ,vl,vu,il,iu,abstol,m,eig,vec,nsize_sb_now &
    &          ,work_lapack,iwork_lapack,ifail_lapack,info)

!    if(iprimddavidson >=2 .and. info/=0) then
    if(info/=0) then
       write(nfout,*) "debug(dspgvx) info=",info
       write(nfout,*) "debug(dspgvx) ifail"
       write(nfout,'(8(1x,i3))') ifail_lapack
       write(nfout,*) "debug(dspgvx) eig"
       write(nfout,'(8(1x,f10.5))') eig
    end if

    deallocate(work_lapack)
    deallocate(iwork_lapack)
    deallocate(ifail_lapack)

    deallocate(ap,bp)

    if(info/=0) then
       !!write(nfout,*) "dspgvx: info=",info
       !!stop 'error in dspgvx_driver'
       ierr = 1
    else
       ierr = 0
    end if

  end subroutine dspgvx_driver_loc
#ifdef NEC_TUNE
  subroutine dsygvx_driver(eig, vec, w1hw2, w1sw2, ierr, nel)
    integer, intent(in) :: nel
    real(kind=DP), intent(out),   dimension(nsize_sb_now)               :: eig
    real(kind=DP), intent(out),   dimension(nsize_sb_now*nel)           :: vec
    real(kind=DP), intent(inout), dimension(nsize_sb_now, nsize_sb_now) :: w1hw2, w1sw2
    integer, intent(out) :: ierr
    integer :: ITYPE, il, iu, lwork, info, m
    character(len=1) :: JOBZ, RANGE, UPLO
    real(kind=DP) :: vl, vu, abstol, work_tmp
    real(kind=DP), allocatable, dimension(:) :: work
    integer, allocatable, dimension(:) :: iwork, ifail
    real(kind=DP), external :: dlamch

    abstol = 2*dlamch('S')
    il = 1; iu = nel
    ITYPE = 1; JOBZ = 'V'; RANGE = 'I'; UPLO = 'U'

    allocate(iwork(5*nsize_sb_now))
    allocate(ifail(nsize_sb_now))

    call dsygvx(ITYPE, JOBZ, RANGE, UPLO, nsize_sb_now, &
   &            w1hw2, nsize_sb_now, w1sw2, nsize_sb_now, &
   &            vl, vu, il, iu, abstol, m, eig, vec, nsize_sb_now, &
   &            work_tmp, -1, iwork, ifail, info)

    lwork = int(work_tmp)
    allocate(work(lwork))

    call dsygvx(ITYPE, JOBZ, RANGE, UPLO, nsize_sb_now, &
   &            w1hw2, nsize_sb_now, w1sw2, nsize_sb_now, &
   &            vl, vu, il, iu, abstol, m, eig, vec, nsize_sb_now, &
   &            work, lwork, iwork, ifail, info)

    if(info /= 0) then
       write(nfout,*) "debug(dsygvx) info=",info
       write(nfout,*) "debug(dsygvx) ifail"
       write(nfout,'(8(1x,i3))') ifail
       write(nfout,*) "debug(dsygvx) eig"
       write(nfout,'(8(1x,f10.5))') eig
    end if

    deallocate(work)
    deallocate(iwork)
    deallocate(ifail)

    if(info /= 0) then
       ierr = 1
    else
       ierr = 0
    end if

  end subroutine dsygvx_driver
#endif

  subroutine zhpgvx_driver(eig,vec,w1hw2,w1sw2,ierr,nel)
    real(kind=DP), intent(out) ,dimension(nsize_sb_now) :: eig
    real(kind=DP), intent(out) ,dimension(nsize_sb_now*kimg,nsize_sb_now) :: vec
    real(kind=DP), intent(in) ,dimension(nsize_mt_now*kimg) :: w1hw2,w1sw2
    integer, intent(out) :: ierr
    integer, intent(in), optional :: nel
    integer :: ITYPE
    character(len=1) :: JOBZ,RANGE,UPLO
    integer :: il,iu
    real(kind=DP),allocatable,dimension(:) :: work_lapack
    real(kind=DP),allocatable,dimension(:) :: rwork_lapack
    integer, allocatable, dimension(:) :: iwork_lapack, ifail_lapack
    real(kind=DP) :: vl,vu,abstol
    integer :: info,m
    real(kind=DP), external :: dlamch
!!$    real(kind=DP), dimension(nsize_mt_now*kimg) :: ap,bp
    real(kind=DP), allocatable, dimension(:) :: ap,bp
    integer :: ib,i

    abstol = 2*dlamch('S')

    il=1; iu=neg
    if(present(nel)) iu=nel
!(LAPACK)  ITYPE = 1:  A*x = (lambda)*B*x, 2:  A*B*x = (lambda)*x, 3:  B*A*x = (lambda)*x
    ITYPE = 1
!(LAPACK)  JOBZ = N : eigenvalue, V : eigenvalue + eigenvector
    JOBZ = 'V'
!(LAPACK)  RANGE = A : all eigenvalues, V: all eigenvalues in (VL, VU], I: the IL-th through IU-th eigenvalues
    RANGE = 'I'
!(LAPACK)  UPLO = U : upper triangle matrix,  L : lower triangle matrix
    UPLO = 'U'
    allocate(work_lapack(2*nsize_sb_now*kimg)); work_lapack = 0.d0
    allocate(rwork_lapack(7*nsize_sb_now)); rwork_lapack = 0.d0
    allocate(iwork_lapack(5*nsize_sb_now)); iwork_lapack = 0
    allocate(ifail_lapack(nsize_sb_now)); ifail_lapack=0
    allocate(ap(nsize_mt_now*kimg));ap=0.d0
    allocate(bp(nsize_mt_now*kimg));bp=0.d0

    if(iprimddavidson >=3 ) then
       write(nfout,*) "debug(zhpgvx) i,w1hw2,w1sw2"
       do ib=1,nsize_sb_now
          i = ib*(ib-1)/2 + ib
          write(nfout,*) ib,w1hw2(2*i-1),w1sw2(2*i-1)
       end do
       ap = w1sw2
       write(nfout,*) "debug(zhpgvx) i,ap,w1sw2"
       do i=1,nsize_mt_now
          write(nfout,*) i,ap(2*i-1),w1sw2(2*i-1)
       end do
       call zhpevx('N','A','U',nsize_sb_now,ap,vl,vu,il,iu,abstol,m &
    &          ,eig,vec,nsize_sb_now &
    &          ,work_lapack,rwork_lapack,iwork_lapack,ifail_lapack,info)
       write(nfout,*) "debug(zhpgvx) eig of w1sw2"
       do ib=1,nsize_sb_now
          write(nfout,*) ib,eig(ib)
       end do
    end if

    ap = w1hw2
    bp = w1sw2

    call zhpgvx(ITYPE,JOBZ,RANGE,UPLO,nsize_sb_now,ap,bp &
    &          ,vl,vu,il,iu,abstol,m,eig,vec,nsize_sb_now &
    &          ,work_lapack,rwork_lapack,iwork_lapack,ifail_lapack,info)

!    if(iprimddavidson >=2 .and. info/=0) then
    if(info/=0) then
       write(nfout,*) "debug(zhpgvx) info=",info
       write(nfout,*) "debug(zhpgvx) ifail"
       write(nfout,'(8(1x,i3))') ifail_lapack
       write(nfout,*) "debug(zhpgvx) eig"
       write(nfout,'(8(1x,f10.5))') eig
    end if

    deallocate(work_lapack)
    deallocate(rwork_lapack)
    deallocate(iwork_lapack)
    deallocate(ifail_lapack)
    deallocate(ap,bp)

    if(info/=0) then
       !!write(nfout,*) "zhpgvx: info=",info
       !!stop 'error in zhpgvx_driver'
       ierr = 1
    else
       ierr = 0
    end if

  end subroutine zhpgvx_driver
#ifdef NEC_TUNE
  subroutine zhegvx_driver(eig, vec, w1hw2, w1sw2, ierr, nel)
    integer, intent(in):: nel
    real(kind=DP), intent(out),   dimension(nsize_sb_now) :: eig
    real(kind=DP), intent(out),   dimension(nsize_sb_now*kimg, nel) :: vec
    real(kind=DP), intent(inout), dimension(nsize_sb_now*kimg, nsize_sb_now) :: w1hw2, w1sw2
    integer, intent(out) :: ierr
    integer :: ITYPE, il, iu, lwork, info, m
    character(len=1) :: JOBZ, RANGE, UPLO
    complex(kind=DP) :: work_tmp
    complex(kind=DP), allocatable, dimension(:) :: work
    real(kind=DP), allocatable, dimension(:) :: rwork
    integer, allocatable, dimension(:) :: iwork, ifail
    real(kind=DP) :: vl, vu, abstol
    real(kind=DP), external :: dlamch

    abstol = 2*dlamch('S')
    il = 1; iu = nel
    ITYPE = 1; JOBZ = 'V'; RANGE = 'I'; UPLO = 'U'

    allocate(rwork(7*nsize_sb_now))
    allocate(iwork(5*nsize_sb_now))
    allocate(ifail(nsize_sb_now))

    call zhegvx(ITYPE, JOBZ, RANGE, UPLO, nsize_sb_now, &
   &            w1hw2, nsize_sb_now, w1sw2, nsize_sb_now, &
   &            vl, vu, il, iu, abstol, m, eig, vec, nsize_sb_now, &
   &            work_tmp, -1, rwork, iwork, ifail, info)

    lwork = int(real(work_tmp))
    allocate(work(lwork))

    call zhegvx(ITYPE, JOBZ, RANGE, UPLO, nsize_sb_now, &
   &            w1hw2, nsize_sb_now, w1sw2, nsize_sb_now, &
   &            vl, vu, il, iu, abstol, m, eig, vec, nsize_sb_now, &
   &            work, lwork, rwork, iwork, ifail, info)

    if(info /= 0) then
       write(nfout,*) "debug(zhegvx) info=",info
       write(nfout,*) "debug(zhegvx) ifail"
       write(nfout,'(8(1x,i3))') ifail
       write(nfout,*) "debug(zhegvx) eig"
       write(nfout,'(8(1x,f10.5))') eig
    end if

    deallocate(work)
    deallocate(rwork)
    deallocate(iwork)
    deallocate(ifail)

    if(info /= 0) then
       ierr = 1
    else
       ierr = 0
    end if

  end subroutine zhegvx_driver
#endif

!  subroutine Vnonlocal_Diagonal_part(ispin,ik,iksnl,vnldi)
!    integer, intent(in)                        :: ispin, ik, iksnl
!    real(kind=DP), intent(out), dimension(kg1) :: vnldi
!
!    integer :: it,mdvdb
!
!    vnldi = 0.d0
!    do it = 1, ntyp
!       mdvdb = m_PP_include_vanderbilt_pot(it)
!       if(mdvdb == SKIP) then
!          call Vnonlocal_D_norm_conserve_case
!       else if(mdvdb == EXECUT) then
!          call Vnonlocal_D_vanderbilt_case
!       end if
!    end do
!  contains
!    subroutine Vnonlocal_D_vanderbilt_case
!      integer       :: ia, p1,p2,lmtt1,il1,lmtt2,il2,i
!      real(kind=DP) :: ph,fac
!
!      do p1 = 1, ilmt(it)
!         lmtt1 = lmtt(p1,it); il1 = ltp(p1,it)
!         do p2 = p1, ilmt(it)
!            lmtt2 = lmtt(p2,it); il2 = ltp(p2,it)
!            if( p1 /= p2) then
!               ph = 2.d0*real(zi**(il2-il1))
!            else
!               ph = 1.d0
!            endif
!            if(mod(il1+il2,2) == 1) cycle
!            do ia = 1, natm
!               if(ityp(ia) /= it) cycle
!               if(ipaw(it)==0) then
!                  fac = ph*iwei(ia) * (dion(p1,p2,it)+vlhxcQ(p1,p2,ia,ispin))
!               else
!                  fac = ph*iwei(ia) * (dion_paw(p1,p2,ispin,ia)+vlhxcQ(p1,p2,ia,ispin))
!               endif
!               do i = 1, iba(ik)
!                  vnldi(i) = vnldi(i)+fac*snl(i,lmtt1,iksnl)*snl(i,lmtt2,iksnl)
!               end do
!            end do
!         end do
!      end do
!    end subroutine Vnonlocal_D_vanderbilt_case
!
!    subroutine Vnonlocal_D_norm_conserve_case
!      integer       :: ia, lmt1,lmt2,lmtt1,il1,im1,il2,im2,i
!      real(kind=DP) :: ph,fac
!
!      ph = 0.d0
!      do ia = 1, natm
!         if(ityp(ia) /= it) cycle
!         ph = ph + iwei(ia)
!      end do
!      do lmt1 = 1, ilmt(it)
!         lmtt1 = lmtt(lmt1,it); il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
!         do lmt2 = lmt1, ilmt(it)
!            il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
!            if(il1 /= il2 .or. im1 /= im2) cycle
!            if(mod(il1+il2,2) == 1) cycle
!            if(ipaw(it)==0)then
!               fac = ph * dion(lmt1,lmt2,it)
!            else
!               fac = ph * dion_paw(lmt1,lmt2,ispin,ia)
!            endif
!            do i = 1, iba(ik)
!               vnldi(i)  = vnldi(i) + fac * snl(i,lmtt1,iksnl)*snl(i,lmtt1,iksnl)
!            end do
!         end do
!      end do
!    end subroutine Vnonlocal_D_norm_conserve_case
!  end subroutine Vnonlocal_Diagonal_part
!
!  subroutine S_Diagonal_part(ik,iksnl,sdiag)
!    integer, intent(in)                        :: ik, iksnl
!    real(kind=DP), intent(out), dimension(kg1) :: sdiag
!
!    integer :: it,mdvdb
!
!    sdiag = 1.d0
!    do it = 1, ntyp
!       mdvdb = m_PP_include_vanderbilt_pot(it)
!       if(mdvdb /= SKIP) call Vanderbilt_case
!    end do
!  contains
!    subroutine Vanderbilt_case
!      integer       :: ia, p1,p2,lmtt1,il1,lmtt2,il2,i
!      real(kind=DP) :: ph,fac
!
!      do p1 = 1, ilmt(it)
!         lmtt1 = lmtt(p1,it); il1 = ltp(p1,it)
!         do p2 = p1, ilmt(it)
!            lmtt2 = lmtt(p2,it); il2 = ltp(p2,it)
!            if( p1 /= p2) then
!               ph = 2.d0*real(zi**(il2-il1))
!            else
!               ph = 1.d0
!            endif
!            if(mod(il1+il2,2) == 1) cycle
!            do ia = 1, natm
!               if(ityp(ia) /= it) cycle
!               fac = ph*iwei(ia) * q(p1,p2,it)
!               do i = 1, iba(ik)
!                  sdiag(i) = sdiag(i)+fac*snl(i,lmtt1,iksnl)*snl(i,lmtt2,iksnl)
!               end do
!            end do
!         end do
!      end do
!    end subroutine Vanderbilt_case
!  end subroutine S_Diagonal_part
!
!  subroutine vlhxc_l_zero_term(vlhxc0,ispin)
!    real(kind=DP), intent(out) :: vlhxc0
!    integer, intent(in)        :: ispin
!
!    if(mype == 0) vlhxc0 = vlhxc_l(1,1,ispin)
!    call mpi_bcast(vlhxc0,1,mpi_double_precision,0,MPI_CommGroup,ierr)
!  end subroutine vlhxc_l_zero_term
!
!  subroutine decide_precon_factor_david(ik,hdiag,sdiag,eig,p)
!    integer, intent(in) :: ik
!    real(kind=DP), intent(in) :: hdiag(kg1),sdiag(kg1),eig
!    real(kind=DP), intent(out) :: p(kg1)
!
!    integer :: ii
!    real(kind=DP) :: denom
!
!    do ii=1,iba(ik)
!       denom = hdiag(ii)-eig*sdiag(ii)
!       if(abs(denom) < eps_mddavid) then
!          denom = sign(eps_mddavid,denom)
!       end if
!       p(ii) = 1.d0/denom
!    end do
!
!  end subroutine decide_precon_factor_david

  subroutine m_ESmdkosugi_Renew_WF(nfout,precon)
    integer, intent(in) :: nfout,precon
    integer             :: ispin, ik, iksnl, switch_of_eko_part
    integer :: iblock,itot
    real(kind=DP), allocatable, dimension(:) ::  afft, bfft
    real(kind=DP), allocatable, dimension(:) :: ekin,p
    real(kind=DP), allocatable, dimension(:) :: afft_l
    real(kind=DP), allocatable, dimension(:,:) :: wk_bfft_l
    real(kind=DP), allocatable, dimension(:,:) :: bfft_l
    integer :: lsize, ibsize, isrsize, fft_l_size
    real(kind=DP), allocatable, dimension(:) :: ekin_l,p_l
    logical :: frestart
    integer :: iblock_now, itot_now, ipri0
    integer :: n_unconv
#ifdef FFT_USE_SSL2_PAD
    lsize = max(nel_fft_x(myrank_g),nel_fft_y(myrank_g),nel_fft_z(myrank_g))
#else
    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
#endif
    allocate(afft_l(lsize*kimg), stat=ierr)
    if(ierr /= 0) then
       write(nfout,*)' m_ESmdkosugi_Renew_WF: Not allocated afft_l array'
       call flush(nfout)
       call mpi_abort(mpi_comm_world, 201, ierr)
    endif
    ibsize = 1

    allocate(ekin_l(maxval(np_g1k)),p_l(maxval(np_g1k)))
    call m_ES_alloc_scss_etc_3D()
    allocate(afft(nfft)); allocate(bfft(nfft))

    call allocate_matrix_ksg_3D
    do ispin = 1, nspin, (af+1)
       call m_ES_Vlocal_in_Rspace_3D(ispin,afft_l,lsize,1,OFF)      ! (ptfft1) vlhxc_l->afft
       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) /= myrank_k) cycle          ! MPI
          isrsize = min(lsize,mp_g1k(ik))
          fft_l_size  = nel_fft_x(myrank_g)
          allocate(wk_bfft_l(lsize*kimg,ibsize) ,stat=ierr)
          allocate(bfft_l(lsize*kimg,ibsize) ,stat=ierr)
          if (ierr /= 0) then
             write(nfout,*)' m_ESmdkosugi_Renew_WF:  Not allocate '
             call flush(nfout)
             call mpi_abort(mpi_comm_world, 205, ierr)
          endif
          iksnl = (ik-1)/nspin + 1
          call allocate_t_matrix_ksg_3D(ik)
          call m_pwBS_kinetic_energies_3D(ik,vkxyz,ekin_l) ! (diakin) ->ekin
          feigconv = .false.
          Loop: do itot=1,max_iter_mdkosugi
             itot_now = itot
             call m_ES_Vnonlocal_W_3D(ik,iksnl,ispin,switch_of_eko_part=ON) ! -> vnlph_l
             call decide_correction_vector_ksg_3D(precon,ik,ekin_l,afft_l,bfft_l, &
                                                  wk_bfft_l,lsize,ibsize,isrsize,fft_l_size,p_l)
             call prepare_Hloc_phi_ksg_3D(ik,ekin_l,afft_l,bfft_l, &
                                          wk_bfft_l,lsize,ibsize,isrsize,fft_l_size,itot)
!print *,ik,itot
write(nfout,*) 'B',feigconv
             Block_Loop: do iblock=1,nblock
                iblock_now = iblock
                call evolve_WFs_in_subspace_ksg_3D &
                          (ik,ispin,ekin_l,afft_l,bfft_l, &
                           wk_bfft_l,lsize,ibsize,isrsize,fft_l_size,iblock,itot,frestart)
             end do Block_Loop
write(nfout,*) 'A',feigconv
!             call m_ES_betar_dot_WFs_4_each_k(nfout,ik)   ! -> fsr_l,fsi_l
          end do Loop
          if(iprimdkosugi>=2) then
             write(nfout,'("MdKosugi: ik=",i5," itot=",i5," subspace=",i5)') &
                                                                 ik, itot_now, nsize_sb_now
          end if
          call deallocate_t_matrix_ksg
! === DEBUG by tkato 2012/06/14 ================================================
          deallocate(wk_bfft_l)
          deallocate(bfft_l)
! ==============================================================================
       enddo      ! k-point loop
    enddo      ! spin loop
    call deallocate_matrix_ksg
!    if(iprimddavidson>=2) then
!       write(nfout,'("Modified Davidson: max_itot=",i5)') max_itot
!    end if
!
! ========================================================================================
! === NOTE: m_ES_sort_eigen_values_3D causes difference with ORG_Parallel!!! =============
! ========================================================================================
!   call m_ES_sort_eigen_values_3D()
! ========================================================================================
! ========================================================================================
! ========================================================================================
!!!  ( in case of af=1 )
!    if(af /= 0) then
!       call cp_eigen_values_for_af       !-(contained here)
!       call expand_neordr_and_nrvf_ordr  !-(contained here)
!    end if

    call get_ipri0(iprimdkosugi,ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko_3D(nfout,mode=SCF)
!
    deallocate(bfft);   deallocate(afft)
    call m_ES_dealloc_scss_etc()
    deallocate(ekin_l,p_l)
!deallocate(vnldi,hdiag,sdiag)

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

    logical function eigenvalues_are_converged(n_unconv)
       integer, intent(out) :: n_unconv
       integer :: ib

       n_unconv = 0
       eigenvalues_are_converged = .true.
       do ib=1,np_e
          if(.not.feigconv(ib)) then
             eigenvalues_are_converged = .false.
             n_unconv = n_unconv + 1
          end if
       end do
    end function eigenvalues_are_converged

  end subroutine m_ESmdkosugi_Renew_WF


  subroutine allocate_matrix_ksg_3D
    integer:: i,j

    nblock = npartition_mdkosugi
    if(np_e .lt. npartition_mdkosugi) nblock = np_e

    allocate(nsize_subspace(nblock))
    allocate(nsize_matrix(nblock))
    allocate(ista_e_l(nblock))
    allocate(iend_e_l(nblock))
    allocate(ielm_e_l(nblock))
    ielm_e_l=np_e/nblock
    j = mod(np_e,nblock)
    do i = 1, j
       ielm_e_l(i) = ielm_e_l(i) + 1
    end do
    ista_e_l(1) = 1
    do i = 2, nblock
       ista_e_l(i) = ista_e_l(i-1) + ielm_e_l(i-1)
       iend_e_l(i-1) = ista_e_l(i) - 1
    end do
    iend_e_l(nblock) = np_e
    do i=1,nblock
      nsize_subspace(i)=ielm_e_l(i)*(max_iter_mdkosugi+1)
      nsize_matrix(i)=nsize_subspace(i)*(nsize_subspace(i)+1)/2
    end do
    msize_subspace=maxval(nsize_subspace)
    msize_matrix=maxval(nsize_matrix)
    allocate(feigconv(np_e))
    allocate(ibover(msize_subspace,nblock))
    allocate(fsr(np_e,nlmta,max_iter_mdkosugi+1))
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       allocate(fsi(np_e,nlmta,max_iter_mdkosugi+1))
    end if
    eps_residual = eps_residual_mdkosugi
allocate(zajold_l(maxval(np_g1k),np_e,kimg))
allocate(fsrold_l(np_e,np_fs))
if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
  allocate(fsiold_l(np_e,np_fs))
end if
  end subroutine allocate_matrix_ksg_3D

  subroutine deallocate_matrix_ksg
    deallocate(nsize_subspace)
    deallocate(nsize_matrix)
    deallocate(ista_e_l)
    deallocate(iend_e_l)
    deallocate(ielm_e_l)
    deallocate(feigconv)
    deallocate(ibover)
    deallocate(fsr)
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       deallocate(fsi)
    end if
deallocate(zajold_l)
deallocate(fsrold_l)
if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
  deallocate(fsiold_l)
end if
  end subroutine deallocate_matrix_ksg

! ================================== added by K. Tagami =============== 11.0
  subroutine deallocate_matrix_ksg_noncl
    deallocate(nsize_subspace)
    deallocate(nsize_matrix)
    deallocate(ista_e_l)
    deallocate(iend_e_l)
    deallocate(ielm_e_l)
    deallocate(feigconv)
    deallocate(ibover)
    deallocate(fsr_noncl)
    if(.not.(kv3/ndim_spinor == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       deallocate(fsi_noncl)
    end if
    deallocate(zajold_l_noncl)
    deallocate(fsrold_l_noncl)
    if(.not.(kv3/ndim_spinor == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       deallocate(fsiold_l_noncl)
    end if
  end subroutine deallocate_matrix_ksg_noncl
! ====================================================================== 11.0

  subroutine allocate_t_matrix_ksg_3D(ik)
    integer, intent(in) :: ik
    integer :: kimg_t
    if(k_symmetry(ik) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = kimg
    end if
    allocate(zat_l(maxval(np_g1k),np_e,kimg,max_iter_mdkosugi+1)) ! MPI
    allocate(zah_l(maxval(np_g1k),np_e,kimg,max_iter_mdkosugi+1)) ! MPI
    allocate(w1hw2(msize_matrix*kimg_t))
    allocate(w1sw2(msize_matrix*kimg_t))
    allocate(wfsd_l(maxval(np_g1k),np_e,ik:ik,kimg)); wfsd_l = 0.d0
    allocate(bsdr_l(np_e,nlmta,ik:ik)); bsdr_l = 0.d0
    allocate(bsdi_l(np_e,nlmta,ik:ik)); bsdi_l = 0.d0
zajold_l(:,:,:)=zaj_l(:,:,ik,:)
fsrold_l(:,:)=fsr_l(:,:,ik)
if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
  fsiold_l(:,:)=fsi_l(:,:,ik)
end if

  end subroutine allocate_t_matrix_ksg_3D

  subroutine deallocate_t_matrix_ksg
! ================================== modified by K. Tagami =========== 11.0
!    deallocate(zat_l) ! MPI
!    deallocate(zah_l) ! MPI

    if ( noncol ) then
       deallocate(zat_l_noncl) ! MPI
       deallocate(zah_l_noncl) ! MPI
    else
       deallocate(zat_l) ! MPI
       deallocate(zah_l) ! MPI
    endif
! ==================================================================== 11.0
    deallocate(w1hw2)
    deallocate(w1sw2)
    deallocate(wfsd_l)
    deallocate(bsdr_l)
    deallocate(bsdi_l)
  end subroutine deallocate_t_matrix_ksg

  subroutine prepare_Hloc_phi_ksg_3D(ik,ekin_l,afft_l,bfft_l, &
                                     wk_bfft_l,lsize,ibsize,isrsize,fft_l_size,itot)
    integer, intent(in) :: ik,lsize,ibsize,isrsize,fft_l_size
    integer, intent(in) :: itot
    real(kind=DP), intent(in)  :: ekin_l(maxval(np_g1k))
    real(kind=DP), intent(in)  :: afft_l(lsize*kimg)
    real(kind=DP), intent(out) :: bfft_l(lsize*kimg,1)
    real(kind=DP), intent(inout) :: wk_bfft_l(lsize*kimg,ibsize)

    integer       :: ib1,i1,ii,ib,iadd
    real(kind=DP) :: denom
    real(kind=DP) :: dr1,dr2,di1,di2,dd
    integer :: id_sname = -1, ipri0
    call tstatc0_begin('prepare_Hloc_phi (mdkosugi) ', id_sname,1)

    call get_ipri0(iprimdkosugi,ipri0)

    denom = 1.d0/product(fft_box_size_WF(1:3,1))

! (zaj_l <- (T+Vloc)|phi> )
!!    zaj_l(:,:,ik,:) = zajold_l(:,:,:,idavid)
!( tenchi ) (zat_l <- zaj_l)

!    zat_l(:,:,:,1) = zaj_l(:,:,ik,:)
    zat_l(:,:,:,itot+1) = wfsd_l(:,:,ik,:)
!    fsr(:,:,1)=fsr_l(:,:,ik)
    fsr(:,:,itot+1)=bsdr_l(:,:,ik)
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
!       fsi(:,:,1)=fsi_l(:,:,ik)
       fsi(:,:,itot+1)=bsdi_l(:,:,ik)
    end if

    if(itot == 1) then
      do ib1 = 1, np_e ! MPI
         ib = ib1
#ifdef __TIMER_COMM__
         call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l,0)
#else
         call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l)
#endif
         call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_y(myrank_g))
         call m_FFT_Direct_3D(nfout,wk_bfft_l,lsize,ibsize)
         call map_fft_to_WF_3D(ik,lsize,ibsize,wk_bfft_l,bfft_l,isrsize,fft_l_size)
         if(kimg == 1) then
            do ii=ista_g1k(ik),iend_g1k(ik)
               iadd = ii - ista_g1k(ik) + 1
               dr1 = zaj_l(iadd,ib,ik,1)
               dr2 = bfft_l(iadd,1)*denom
               zah_l(iadd,ib,1,1) = ekin_l(iadd)*dr1+dr2
            enddo
         else
            do ii=ista_g1k(ik),iend_g1k(ik)
               iadd = ii - ista_g1k(ik) + 1
               dr1  = zaj_l(iadd,ib,ik,1)
               di1  = zaj_l(iadd,ib,ik,kimg)
               zah_l(iadd,ib,1,   1)= ekin_l(iadd)*dr1+bfft_l(2*iadd-1,1)*denom
               zah_l(iadd,ib,kimg,1)= ekin_l(iadd)*di1+bfft_l(2*iadd,  1)*denom
            enddo
         endif
      enddo
    end if

!!( tenchi ) (zah_l <- zaj_l)
    zaj_l(:,:,ik,:) = wfsd_l(:,:,ik,:)

    do ib1 = 1, np_e ! MPI
       ib = ib1
#ifdef __TIMER_COMM__
       call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l,0)
#else
       call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l)
#endif
       call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_y(myrank_g))
       call m_FFT_Direct_3D(nfout,wk_bfft_l,lsize,ibsize)
       call map_fft_to_WF_3D(ik,lsize,ibsize,wk_bfft_l,bfft_l,isrsize,fft_l_size)
       if(kimg == 1) then
          do ii=ista_g1k(ik),iend_g1k(ik)
             iadd = ii - ista_g1k(ik) + 1
             dr1 = zaj_l(iadd,ib,ik,1)
             dr2 = bfft_l(iadd,1)*denom
             zah_l(iadd,ib,1,itot+1) = ekin_l(iadd)*dr1+dr2
          enddo
       else
          do ii=ista_g1k(ik),iend_g1k(ik)
             iadd = ii - ista_g1k(ik) + 1
             dr1  = zaj_l(iadd,ib,ik,1)
             di1  = zaj_l(iadd,ib,ik,kimg)
             zah_l(iadd,ib,1,   itot+1) = ekin_l(iadd)*dr1+bfft_l(2*iadd-1,1)*denom
             zah_l(iadd,ib,kimg,itot+1) = ekin_l(iadd)*di1+bfft_l(2*iadd,  1)*denom
          enddo
       endif
    enddo

    call tstatc0_end(id_sname)

  contains

    subroutine get_ipri0(ipri_in, ipri_out)
      integer, intent(in)  :: ipri_in
      integer, intent(out) :: ipri_out
      if(npes > 1) then
         if(mype == 0) ipri_out = ipri_in
         call mpi_bcast(ipri_out,1,mpi_integer,0,mpi_k_world(myrank_k),ierr)
      else
         ipri_out = ipri_in
      end if
    end subroutine get_ipri0

  end subroutine prepare_Hloc_phi_ksg_3D


  subroutine evolve_WFs_in_subspace_ksg_3D &
                                (ik,ispin,ekin_l,afft_l,bfft_l, &
                                 wk_bfft_l,lsize,ibsize,isrsize,fft_l_size,iblock,itot,frestart)
    integer, intent(in) :: ik,ispin,lsize,ibsize,isrsize,fft_l_size
    integer, intent(in) :: iblock,itot
    logical, intent(out) :: frestart
    real(kind=DP), intent(in)  :: ekin_l(maxval(np_g1k))
    real(kind=DP), intent(in)  :: afft_l(lsize*kimg)
    real(kind=DP), intent(out) :: bfft_l(lsize*kimg,1)
    real(kind=DP), intent(inout) :: wk_bfft_l(lsize*kimg,ibsize)
! (allocatable variables)
    real(kind=DP), allocatable,dimension(:) ::     eig
    real(kind=DP), allocatable,dimension(:,:) ::   vec
    real(kind=DP), allocatable,dimension(:) ::     eko_d
    real(kind=DP), allocatable,dimension(:) ::     eko_d_mpi
    integer, allocatable,dimension(:) ::     occup
    integer, allocatable,dimension(:) ::     occup_mpi

    integer       :: ib1,ib2,ib1to,ib2to,i1,ii,ri,ib,iadd
    integer       :: ibb1,ibb2
    integer       :: ii1,ii2,iter,iter1,iter2
    real(kind=DP) :: eko1, eko2, ekod
    real(kind=DP) :: hr2,hi2,dr1,dr2,di1,di2,dd
    integer :: ip0,ip0b,ip1,ip1b,ib1n,ib2n,ndata,nshift,kimg_t,ig1
    integer :: noffset
    integer :: nsize_max_sb_now
    integer :: ierr_diag
    integer :: nel,nsta,nend
    integer :: id_sname = -1, ipri0
#ifdef NEC_TUNE
    integer :: num, k, lda, ldb, ldc
    real(kind=DP), allocatable, dimension(:,:) :: w1hw2_, w1sw2_
    real(kind=DP), allocatable, dimension(:,:,:) :: zat_t, zah_t
    real(kind=DP), allocatable, dimension(:,:) :: zat_t1, zah_t1
    complex(kind=DP), parameter :: c0 = (0.0d0, 0.0d0), c1 = (1.0d0, 0.0d0)
#endif
    call tstatc0_begin('evolve_WFs_in_subspace (modified kosugi) ', id_sname,1)

    call get_ipri0(iprimdkosugi,ipri0)

    if(k_symmetry(ik) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = kimg
    end if

    nel =ielm_e_l(iblock)
    nsta=ista_e_l(iblock)
    nend=iend_e_l(iblock)

    if(itot==1) then
      do ib=1,nel
        ibover(ib,iblock)=ib
      end do
      nsize_subspace(iblock)=nel
    end if

!    ip0=nsize_subspace(iblock)
!    do ib=1,ip0
!      ibover(ib,iblock) = ib
!    end do

    noffset=nel*itot
    ip0=nsize_subspace(iblock)
    do ib=1,nel
      if(.not.feigconv(nsta+ib-1)) then
        ip0=ip0+1
        ibover(noffset+ib,iblock) = ip0
      else
        ibover(noffset+ib,iblock) = -1
      end if
    end do

    nsize_sb_now = ip0
    nsize_subspace(iblock) = ip0
    nsize_mt_now = nsize_sb_now*(nsize_sb_now+1)/2
    nsize_max_sb_now = nel*(itot+1)                   !nsize_subspace(iblock)

    if(iprimdkosugi >=2) then
       write(nfout,*) 'ibover=',ibover(1:nsize_max_sb_now,iblock)
    end if
!
    allocate(eig(nsize_sb_now)); eig=0.d0
    allocate(vec(nsize_sb_now*kimg_t,nsize_sb_now))

    if(iprimdkosugi >=2) then
       write(nfout,*) 'Modified Kosugi:ik,iblock,nsize_sb_now=', ik,iblock,nsize_sb_now
    end if

    eko1 = sum(eko_l(nsta:nend,ik))

!! (make matrix elements )
!    ! <n|T+Vloc|m> !
#ifdef NEC_TUNE
    allocate(w1hw2_(nsize_sb_now*kimg_t,nsize_sb_now))
    allocate(w1sw2_(nsize_sb_now*kimg_t,nsize_sb_now))

    if(kimg == 1) then
       allocate(zah_t(maxval(np_g1k),nsize_sb_now,1))
       allocate(zat_t(maxval(np_g1k),nsize_sb_now,1))
       num = 0
       do ibb1 = 1, nsize_max_sb_now
          if(ibover(ibb1,iblock) < 0) cycle
          iter1 = (ibb1-1)/nel+1
          ii1  = ibb1-nel*(iter1-1)
          ii1  = nsta+ii1-1
          num = num + 1
          zah_t(:,num,1) = zah_l(:,ii1,1,iter1)
          zat_t(:,num,1) = zat_l(:,ii1,1,iter1)
       end do
       k = np_g1k(ik); lda = maxval(np_g1k); ldb = maxval(np_g1k); ldc = nsize_sb_now
       call dgemm('T','N',num,num,k,1.0d0,zat_t,lda,zah_t,ldb,0.0d0,w1hw2_,ldc)
       call dgemm('T','N',num,num,k,1.0d0,zat_t,lda,zat_t,ldb,0.0d0,w1sw2_,ldc)
    else
      if(k_symmetry(ik) == GAMMA) then
         allocate(zah_t(maxval(np_g1k),nsize_sb_now,kimg))
         allocate(zat_t(maxval(np_g1k),nsize_sb_now,kimg))
         allocate(zah_t1(nsize_sb_now,kimg))
         allocate(zat_t1(nsize_sb_now,kimg))
         num = 0
         do ibb1 = 1, nsize_max_sb_now
            if(ibover(ibb1,iblock) < 0) cycle
            iter1 = (ibb1-1)/nel+1
            ii1  = ibb1-nel*(iter1-1)
            ii1  = nsta+ii1-1
            num = num + 1
            zah_t(:,num,1) = zah_l(:,ii1,1,iter1)
            zah_t(:,num,2) = zah_l(:,ii1,2,iter1)
            zat_t(:,num,1) = zat_l(:,ii1,1,iter1)
            zat_t(:,num,2) = zat_l(:,ii1,2,iter1)
         end do
         if(ista_g1k(ik) == 1) then
            zat_t1(:,:) = zat_t(1,:,:)
            zah_t1(:,:) = zah_t(1,:,:)
            zat_t(1,:,:) = zat_t(1,:,:)/sqrt(2.0d0)
            zah_t(1,:,:) = zah_t(1,:,:)/sqrt(2.0d0)
         end if
         k = np_g1k(ik); lda = maxval(np_g1k); ldb = maxval(np_g1k); ldc = nsize_sb_now
         call dgemm('T','N',num,num,k,2.0d0,zat_t(1,1,1),lda,zah_t(1,1,1),ldb,0.0d0,w1hw2_,ldc)
         call dgemm('T','N',num,num,k,2.0d0,zat_t(1,1,2),lda,zah_t(1,1,2),ldb,1.0d0,w1hw2_,ldc)
         call dgemm('T','N',num,num,k,2.0d0,zat_t(1,1,1),lda,zat_t(1,1,1),ldb,0.0d0,w1sw2_,ldc)
         call dgemm('T','N',num,num,k,2.0d0,zat_t(1,1,2),lda,zat_t(1,1,2),ldb,1.0d0,w1sw2_,ldc)
         if(ista_g1k(ik) == 1) then
            zat_t(1,:,:) = zat_t1(:,:)
            zah_t(1,:,:) = zah_t1(:,:)
         end if
         deallocate(zah_t1)
         deallocate(zat_t1)
      else
         allocate(zah_t(maxval(np_g1k)*kimg,nsize_sb_now,1))
         allocate(zat_t(maxval(np_g1k)*kimg,nsize_sb_now,1))
         num = 0
         do ibb1 = 1, nsize_max_sb_now
            if(ibover(ibb1,iblock) < 0) cycle
            iter1 = (ibb1-1)/nel+1
            ii1  = ibb1-nel*(iter1-1)
            ii1  = nsta+ii1-1
            num = num + 1
            do ii = 1, np_g1k(ik)
               zah_t(2*ii-1,num,1) = zah_l(ii,ii1,1,iter1)
               zah_t(2*ii,  num,1) = zah_l(ii,ii1,2,iter1)
               zat_t(2*ii-1,num,1) = zat_l(ii,ii1,1,iter1)
               zat_t(2*ii,  num,1) = zat_l(ii,ii1,2,iter1)
            end do
         end do
         k = np_g1k(ik); lda = maxval(np_g1k); ldb = maxval(np_g1k); ldc = nsize_sb_now
         call zgemm('C','N',num,num,k,c1,zat_t,lda,zah_t,ldb,c0,w1hw2_,ldc)
         call zgemm('C','N',num,num,k,c1,zat_t,lda,zat_t,ldb,c0,w1sw2_,ldc)
      end if
    end if
#else
    do ibb2 = 1,nsize_max_sb_now
       if(ibover(ibb2,iblock)<0) cycle
       ib2 = ibover(ibb2,iblock)
       iter2 = (ibb2-1)/nel+1
       ii2  = ibb2-nel*(iter2-1)
       ii2  = nsta+ii2-1
       ip0b = ib2*(ib2-1)/2
       do ibb1 = 1,ibb2
          if(ibover(ibb1,iblock)<0) cycle
          ib1 = ibover(ibb1,iblock)
          iter1 = (ibb1-1)/nel+1
          ii1 = ibb1-nel*(iter1-1)
          ii1 = nsta+ii1-1
          ip0 = ip0b + ib1
          if(kimg == 1) then
             w1hw2(ip0) = 0.d0
             w1sw2(ip0) = 0.d0
             do ii = ista_g1k(ik),iend_g1k(ik) ! MPI
                iadd = ii - ista_g1k(ik) + 1
                hr2 = zah_l(iadd,ii2,1,iter2)
                dr2 = zat_l(iadd,ii2,1,iter2)
                dr1 = zat_l(iadd,ii1,1,iter1)
                w1hw2(ip0) = w1hw2(ip0) + dr1*hr2
                w1sw2(ip0) = w1sw2(ip0) + dr1*dr2
             end do
             call mpi_allreduce(MPI_IN_PLACE,w1hw2(ip0),1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
             call mpi_allreduce(MPI_IN_PLACE,w1sw2(ip0),1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
          else
             if(k_symmetry(ik) == GAMMA) then
                w1hw2(ip0) = 0.d0
                w1sw2(ip0) = 0.d0
                do ii = max(2,ista_g1k(ik)), iend_g1k(ik) ! MPI
                   iadd = ii - ista_g1k(ik) + 1
                   hr2 = zah_l(iadd,ii2,1,iter2) ! MPI
                   hi2 = zah_l(iadd,ii2,2,iter2) ! MPI
                   dr2 = zat_l(iadd,ii2,1,iter2) ! MPI
                   di2 = zat_l(iadd,ii2,2,iter2) ! MPI
                   dr1 = zat_l(iadd,ii1,1,iter1) ! MPI
                   di1 = zat_l(iadd,ii1,2,iter1) ! MPI
                   w1hw2(ip0) =w1hw2(ip0)+(dr1*hr2+di1*hi2)*2.d0
                   w1sw2(ip0) =w1sw2(ip0)+(dr1*dr2+di1*di2)*2.d0
                end do
                if(ista_g1k(ik) == 1) then
                   hr2 = zah_l(1,ii2,1,iter2) ! MPI
                   hi2 = zah_l(1,ii2,2,iter2) ! MPI
                   dr2 = zat_l(1,ii2,1,iter2) ! MPI
                   di2 = zat_l(1,ii2,2,iter2) ! MPI
                   dr1 = zat_l(1,ii1,1,iter1) ! MPI
                   di1 = zat_l(1,ii1,2,iter1) ! MPI
                   w1hw2(ip0) =w1hw2(ip0)+dr1*hr2+di1*hi2
                   w1sw2(ip0) =w1sw2(ip0)+dr1*dr2+di1*di2
                end if
                call mpi_allreduce(MPI_IN_PLACE,w1hw2(ip0),1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
                call mpi_allreduce(MPI_IN_PLACE,w1sw2(ip0),1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
             else
                w1hw2(2*ip0-1:2*ip0) = 0.d0
                w1sw2(2*ip0-1:2*ip0) = 0.d0
                do ii = ista_g1k(ik), iend_g1k(ik)
                   iadd = ii - ista_g1k(ik) + 1
                   hr2 = zah_l(iadd,ii2,1,iter2) ! MPI
                   hi2 = zah_l(iadd,ii2,2,iter2) ! MPI
                   dr2 = zat_l(iadd,ii2,1,iter2) ! MPI
                   di2 = zat_l(iadd,ii2,2,iter2) ! MPI
                   dr1 = zat_l(iadd,ii1,1,iter1) ! MPI
                   di1 = zat_l(iadd,ii1,2,iter1) ! MPI
                   w1hw2(2*ip0-1) =w1hw2(2*ip0-1)+dr1*hr2+di1*hi2
                   w1hw2(2*ip0  ) =w1hw2(2*ip0  )+dr1*hi2-di1*hr2
                   w1sw2(2*ip0-1) =w1sw2(2*ip0-1)+dr1*dr2+di1*di2
                   w1sw2(2*ip0  ) =w1sw2(2*ip0  )+dr1*di2-di1*dr2
                end do
                call mpi_allreduce(MPI_IN_PLACE,w1hw2(2*ip0-1),2,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
                call mpi_allreduce(MPI_IN_PLACE,w1sw2(2*ip0-1),2,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
             end if
          end if
       end do
    end do
#endif
    if(iprimdkosugi >= 3) call wd_w1hw2(" -- w1hw2 without nl part--",iblock)
    ! <n|Vnl|m>
!   if(myrank_g == 0) then
    call add_nonlocal_part ! w1hw2 = w1hw2 + w1Vnlw2
                           ! w1sw2 = w1sw2 + w1qw2
!   endif
#ifdef NEC_TUNE
    call mpi_allreduce(MPI_IN_PLACE,w1hw2_,nsize_sb_now*nsize_sb_now*kimg_t,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
    call mpi_allreduce(MPI_IN_PLACE,w1sw2_,nsize_sb_now*nsize_sb_now*kimg_t,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
#endif
    if(iprimdkosugi >= 3) call wd_w1hw2(" -- w1hw2 with nl part--",iblock)
!
    if(iprimdkosugi >= 2) call wd_w1hw2(" -- just after making w1hw2 --",iblock)
!    if(ipridavidson >= 2) then
!       write(nfout,*) 'neordr for ik = ',ik
!       write(nfout,9002) (neordr(ib1,ik),ib1=1,neg)
!       write(nfout,*) 'nrvf_ordr for ik = ',ik
!       write(nfout,9002) (nrvf_ordr(ib1,ik),ib1=1,neg)
!       write(nfout,*) 'eig'
!       write(nfout,'(5x,10f8.4)') (eko_l(ib1,ik),ib1=1,np_e)
!    endif
!9002 format(5x,10i8)

!! (Diagonalization )  !!

#ifdef NEC_TUNE
    if(kimg_t == 1) then
       call dsygvx_driver(eig, vec, w1hw2_, w1sw2_, ierr_diag, nel)
    else
       call zhegvx_driver(eig, vec, w1hw2_, w1sw2_, ierr_diag, nel)
    endif
    deallocate(w1hw2_)
    deallocate(w1sw2_)
#else
    if(kimg_t == 1) then
       call dspgvx_driver_loc(eig,vec,w1hw2,w1sw2,ierr_diag,nel)
    else
       call zhpgvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag,nel)
    endif
#endif

    frestart = .false.
    if(ierr_diag /= 0) then
!       zaj_l(:,:,ik,:) = zaj_l_backup(:,:,:)
!       do ib1 = 1, neg
!          if(map_e(ib1) == myrank_e) then         ! MPI
!             eko_l(map_z(ib1),ik)=eko_d(ib1)
!          end if
!       end do
      frestart = .true.
      zaj_l(:,nsta:nend,ik,:)=zat_l(:,nsta:nend,:,1)
      if(iprimdkosugi >= 2) then
        write(nfout,*) '** restart Modified Kosugi iteration **'
        write(nfout,*) 'ik=',ik,' iblock=',iblock,' itot=',itot
      end if
      write(nfout,*) '** restart Modified Kosugi **'
      write(nfout,*) 'ik=',ik,' iblock=',iblock,' itot=',itot
      call wd_w1hw2(" -- restart Modified Kosugi iteration --",iblock)
!!$ print *,'Restart'
      goto 9000
    else

      feigconv(nsta:nend) = .false.
      do ib=1,nel
         if(occup_l(nsta+ib-1,ik) > 0.d0) then
            if(abs(eko_l(nsta+ib-1,ik)-eig(ib)) < delta_eig_occup_mdkosugi) &
                                                   feigconv(nsta+ib-1) = .true.
         else
            if(abs(eko_l(nsta+ib-1,ik)-eig(ib)) < delta_eig_empty_mdkosugi) &
                                                   feigconv(nsta+ib-1) = .true.
         end if
      end do
       if(ipri0 >= 2) then
          write(nfout,*) 'eko_l for ik = ',ik
          write(nfout,*) 'iblock       = ',iblock
          write(nfout,9001) (eko_l(nsta+ib-1,ik),ib=1,nel)
          write(nfout,*) 'eig for ik = ',ik
          write(nfout,9001) (eig(ib),ib=1,nel)
          call wd_w1hw2(" -- after diagonalization --",iblock)
!sum eko
          dr1=0.d0;dr2=0.d0
          do ib1=1,nel
             dr1=dr1+eko_l(nsta+ib1-1,ik) ! MPI
             dr2=dr2+eig(ib1)
          enddo
          write(nfout,'(" sum of eko_l, eig, abs diff =",3e25.10)') dr1,dr2,abs(dr2-dr1)
       endif
!!! (subspace rotation) !!
       call subspace_rotation ! vec,zat_l -> zat_l
!!( tenchi ) (zaj_l <- zat_l)
!       iter = min(idavid+1,ndavid)
!       call m_ES_W_transpose_back(ista_k,iend_k,ik,zaj_l,zat_l(1,1,1,iter))
!!       zaj_l_backup(:,nsta:nend,:) = zaj_l(:,nsta:nend,ik,:)
!! (eko_l)
       do ib1 = 1, nel
             eko_l(nsta+ib1-1,ik)=eig(ib1)
       end do
       if(iprimdkosugi >= 2) then
          eko2 = sum(eig(1:nel))
          write(nfout,1201) ik,eko1,ekod,eko2

          write(nfout,*) 'eko_l'
          write(nfout,9001) (eko_l(nsta+ib1-1,ik),ib1=1,nel)
       endif
1201 format(' %% for ik = ',i4,4x,' eko1&ekod&eko2 = ',3f14.7)
9001 format(5x,6f12.5)
!!! (neordr & nrvf_ordr)
!
    end if
9000 continue
!    neordr(1:neg,ik) = (/(ib1,ib1=1,neg)/)
!    nrvf_ordr(1:neg,ik) = (/(ib1,ib1=1,neg)/)
!
!! (deallocate)
#ifdef NEC_TUNE
    deallocate(zah_t)
    deallocate(zat_t)
#endif
    deallocate(eig)
    deallocate(vec)
!
    call tstatc0_end(id_sname)
!
  contains
    subroutine get_ipri0(ipri_in, ipri_out)
      integer, intent(in)  :: ipri_in
      integer, intent(out) :: ipri_out
      if(npes > 1) then
         if(mype == 0) ipri_out = ipri_in
         call mpi_bcast(ipri_out,1,mpi_integer,0,mpi_k_world(myrank_k),ierr)
      else
         ipri_out = ipri_in
      end if
    end subroutine get_ipri0
!
    subroutine wd_w1hw2(somecomment,iblock)
      character(len=*), intent(in) :: somecomment
      integer,intent(in) :: iblock
      integer :: ib1, ib2, nel_wd, nsb_wd
      write(nfout,'(a35)') somecomment
      write(nfout,*) 'w1hw2 for ik = ',ik
      nel_wd = 8
      nsb_wd = 8
      if(nel_wd > nel) nel_wd = nel
      if(nsb_wd > nsize_sb_now) nsb_wd = nsize_sb_now
      if(kimg_t==1) then
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1hw2(ip0),ip0=ip0b+1,ip0b+ib2)
         end do
      else
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1hw2(2*ip0-1),ip0=ip0b+1,ip0b+ib2)
         end do
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1hw2(2*ip0),ip0=ip0b+1,ip0b+ib2)
         end do
      end if
      write(nfout,*) 'w1sw2 for ik = ',ik
      if(kimg_t==1) then
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1sw2(ip0),ip0=ip0b+1,ip0b+ib2)
         end do
      else
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1sw2(2*ip0-1),ip0=ip0b+1,ip0b+ib2)
         end do
         do ib2=1,nsb_wd
            ip0b = ib2*(ib2-1)/2
            write(nfout,9001) (w1sw2(2*ip0),ip0=ip0b+1,ip0b+ib2)
         end do
      end if
9001  format(5x,9f12.5)
      write(nfout,*) 'eko_l for ik = ',ik
      write(nfout,9001) (eko_l(nsta+ib1-1,ik),ib1=1,nel_wd)
    end subroutine wd_w1hw2

    subroutine add_nonlocal_part
      integer :: ip,ib1,ib2,ibb1,ibb2
      integer       :: ia, lmt1, lmt2, it, p, s, ib
      real(kind=DP) :: facv,facq,vr,vi,qr,qi
      real(kind=DP) :: tmpr,tmpi
#ifdef NEC_TUNE
      real(kind=DP), dimension(nlmta,nel,itot+1) :: fsr_t, fsr_v, fsr_q
      real(kind=DP), dimension(nlmta,nel,itot+1) :: fsi_t, fsi_v, fsi_q
      real(kind=DP), allocatable, dimension(:,:) :: fsr_tt, fsr_vt, fsr_qt
      real(kind=DP), allocatable, dimension(:,:) :: fsi_tt, fsi_vt, fsi_qt
      integer :: iter, ncount, i

      if(kimg == 1) then
         fsr_v = 0.0d0
         fsr_q = 0.0d0
         fsi_v = 0.0d0
         fsi_q = 0.0d0
         do iter = 1, itot+1
            ncount = 0
            do ia = ista_atm, iend_atm
               it = ityp(ia)
               do lmt1 = 1, ilmt(it)
                  ncount = ncount + 1
                  p = lmta(lmt1,ia)
                  do ibb1 = 1, nel
                     ibb2 = nsta+ibb1-1
                     fsr_t(ncount,ibb1,iter) = fsr(ibb2,p,iter)
                     fsi_t(ncount,ibb1,iter) = fsi(ibb2,p,iter)
                  end do ! ibb1
                  do lmt2 = 1, ilmt(it)
                     s = lmta(lmt2,ia)
                     if(ipaw(it) == 0)then
                        facv   = iwei(ia)*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin))
                     else
                        facv   = iwei(ia)*(dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin))
                     endif
                     facq   = iwei(ia)*q(lmt1,lmt2,it)
                     do ibb1 = 1, nel
                        ibb2 = nsta+ibb1-1
                        fsr_v(ncount,ibb1,iter) = fsr_v(ncount,ibb1,iter) + facv*fsr(ibb2,s,iter)
                        fsr_q(ncount,ibb1,iter) = fsr_q(ncount,ibb1,iter) + facq*fsr(ibb2,s,iter)
                        fsi_v(ncount,ibb1,iter) = fsi_v(ncount,ibb1,iter) + facv*fsi(ibb2,s,iter)
                        fsi_q(ncount,ibb1,iter) = fsi_q(ncount,ibb1,iter) + facq*fsi(ibb2,s,iter)
                     end do ! ibb1
                  end do ! lmt2
               end do ! lmt1
            end do ! ia
         end do ! iter
         allocate(fsr_tt(ncount,nsize_sb_now))
         allocate(fsr_vt(ncount,nsize_sb_now))
         allocate(fsr_qt(ncount,nsize_sb_now))
         allocate(fsi_tt(ncount,nsize_sb_now))
         allocate(fsi_vt(ncount,nsize_sb_now))
         allocate(fsi_qt(ncount,nsize_sb_now))
         num = 0
         do ibb1 = 1, nsize_max_sb_now
            if(ibover(ibb1,iblock)<0) cycle
            iter1=(ibb1-1)/nel+1
            ii1  = ibb1-nel*(iter1-1)
            num = num + 1
            fsr_tt(1:ncount,num) = fsr_t(1:ncount,ii1,iter1)
            fsr_vt(1:ncount,num) = fsr_v(1:ncount,ii1,iter1)
            fsr_qt(1:ncount,num) = fsr_q(1:ncount,ii1,iter1)
            fsi_tt(1:ncount,num) = fsi_t(1:ncount,ii1,iter1)
            fsi_vt(1:ncount,num) = fsi_v(1:ncount,ii1,iter1)
            fsi_qt(1:ncount,num) = fsi_q(1:ncount,ii1,iter1)
         end do
         ldc = nsize_sb_now
         call dgemm('T','N',num,num,ncount,1.0d0,fsr_tt,ncount,fsr_vt,ncount,1.0d0,w1hw2_,ldc)
         call dgemm('T','N',num,num,ncount,1.0d0,fsi_tt,ncount,fsi_vt,ncount,1.0d0,w1hw2_,ldc)
         call dgemm('T','N',num,num,ncount,1.0d0,fsr_tt,ncount,fsr_qt,ncount,1.0d0,w1sw2_,ldc)
         call dgemm('T','N',num,num,ncount,1.0d0,fsi_tt,ncount,fsi_qt,ncount,1.0d0,w1sw2_,ldc)
         deallocate(fsr_tt)
         deallocate(fsr_vt)
         deallocate(fsr_qt)
         deallocate(fsi_tt)
         deallocate(fsi_vt)
         deallocate(fsi_qt)
      else ! if(kimg == 1)
         if(k_symmetry(ik) == GAMMA) then
            fsr_v = 0.0d0
            fsr_q = 0.0d0
            do iter = 1, itot+1
               ncount = 0
               do ia = ista_atm, iend_atm
                  it = ityp(ia)
                  do lmt1 = 1, ilmt(it)
                     ncount = ncount + 1
                     p = lmta(lmt1,ia)
                     do ibb1 = 1, nel
                        ibb2 = nsta+ibb1-1
                        fsr_t(ncount,ibb1,iter) = fsr(ibb2,p,iter)
                     end do ! ibb1
                     do lmt2 = 1, ilmt(it)
                        s = lmta(lmt2,ia)
                        if(ipaw(it) == 0)then
                           facv   = iwei(ia)*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin))
                        else
                           facv   = iwei(ia)*(dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin))
                        endif
                        facq   = iwei(ia)*q(lmt1,lmt2,it)
                        do ibb1 = 1, nel
                           ibb2 = nsta+ibb1-1
                           fsr_v(ncount,ibb1,iter) = fsr_v(ncount,ibb1,iter) + facv*fsr(ibb2,s,iter)
                           fsr_q(ncount,ibb1,iter) = fsr_q(ncount,ibb1,iter) + facq*fsr(ibb2,s,iter)
                        end do ! ibb1
                     end do ! lmt2
                  end do ! lmt1
               end do ! ia
            end do ! iter
            allocate(fsr_tt(ncount,nsize_sb_now))
            allocate(fsr_vt(ncount,nsize_sb_now))
            allocate(fsr_qt(ncount,nsize_sb_now))
            num = 0
            do ibb1 = 1, nsize_max_sb_now
               if(ibover(ibb1,iblock)<0) cycle
               iter1=(ibb1-1)/nel+1
               ii1  = ibb1-nel*(iter1-1)
               num = num + 1
               fsr_tt(1:ncount,num) = fsr_t(1:ncount,ii1,iter1)
               fsr_vt(1:ncount,num) = fsr_v(1:ncount,ii1,iter1)
               fsr_qt(1:ncount,num) = fsr_q(1:ncount,ii1,iter1)
            end do
            ldc = nsize_sb_now
            call dgemm('T','N',num,num,ncount,1.0d0,fsr_tt,ncount,fsr_vt,ncount,1.0d0,w1hw2_,ldc)
            call dgemm('T','N',num,num,ncount,1.0d0,fsr_tt,ncount,fsr_qt,ncount,1.0d0,w1sw2_,ldc)
            deallocate(fsr_tt)
            deallocate(fsr_vt)
            deallocate(fsr_qt)
         else ! if(k_symmetry(ik) == GAMMA)
            fsr_v = 0.0d0
            fsr_q = 0.0d0
            fsi_v = 0.0d0
            fsi_q = 0.0d0
            do iter = 1, itot+1
               ncount = 0
               do ia = ista_atm, iend_atm
                  it = ityp(ia)
                  do lmt1 = 1, ilmt(it)
                     ncount = ncount + 1
                     p = lmta(lmt1,ia)
                     do ibb1 = 1, nel
                        ibb2 = nsta+ibb1-1
                        fsr_t(ncount,ibb1,iter) = fsr(ibb2,p,iter)
                        fsi_t(ncount,ibb1,iter) = fsi(ibb2,p,iter)
                     end do ! ibb1
                     do lmt2 = 1, ilmt(it)
                        s = lmta(lmt2,ia)
                        if(ipaw(it) == 0)then
                           facv   = iwei(ia)*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin))
                        else
                           facv   = iwei(ia)*(dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin))
                        endif
                        facq   = iwei(ia)*q(lmt1,lmt2,it)
                        do ibb1 = 1, nel
                           ibb2 = nsta+ibb1-1
                           fsr_v(ncount,ibb1,iter) = fsr_v(ncount,ibb1,iter) + facv*fsr(ibb2,s,iter)
                           fsr_q(ncount,ibb1,iter) = fsr_q(ncount,ibb1,iter) + facq*fsr(ibb2,s,iter)
                           fsi_v(ncount,ibb1,iter) = fsi_v(ncount,ibb1,iter) + facv*fsi(ibb2,s,iter)
                           fsi_q(ncount,ibb1,iter) = fsi_q(ncount,ibb1,iter) + facq*fsi(ibb2,s,iter)
                        end do ! ibb1
                     end do ! lmt2
                  end do ! lmt1
               end do ! ia
            end do ! iter
            allocate(fsr_tt(ncount*kimg,nsize_sb_now))
            allocate(fsr_vt(ncount*kimg,nsize_sb_now))
            allocate(fsr_qt(ncount*kimg,nsize_sb_now))
            num = 0
            do ibb1 = 1, nsize_max_sb_now
               if(ibover(ibb1,iblock)<0) cycle
               iter1=(ibb1-1)/nel+1
               ii1  = ibb1-nel*(iter1-1)
               num = num + 1
               do i = 1, ncount
                  fsr_tt(i*2-1,num) = fsr_t(i,ii1,iter1)
                  fsr_tt(i*2,  num) = fsi_t(i,ii1,iter1)
                  fsr_vt(i*2-1,num) = fsr_v(i,ii1,iter1)
                  fsr_vt(i*2,  num) = fsi_v(i,ii1,iter1)
                  fsr_qt(i*2-1,num) = fsr_q(i,ii1,iter1)
                  fsr_qt(i*2,  num) = fsi_q(i,ii1,iter1)
               end do
            end do
            ldc = nsize_sb_now
            call zgemm('C','N',num,num,ncount,c1,fsr_tt,ncount,fsr_vt,ncount,c1,w1hw2_,ldc)
            call zgemm('C','N',num,num,ncount,c1,fsr_tt,ncount,fsr_qt,ncount,c1,w1sw2_,ldc)
            deallocate(fsr_tt)
            deallocate(fsr_vt)
            deallocate(fsr_qt)
         end if ! if(k_symmetry(ik) == GAMMA)
      end if ! if(kimg == 1)
#else
      do ibb2 = 1,nsize_max_sb_now
         if(ibover(ibb2,iblock)<0) cycle
         ib2 = ibover(ibb2,iblock)
         iter2= (ibb2-1)/nel+1
         ii2  = ibb2-nel*(iter2-1)
         ii2  = nsta+ii2-1
         ip0b = ib2*(ib2-1)/2
         do ibb1 = 1,ibb2
            if(ibover(ibb1,iblock)<0) cycle
            ib1 = ibover(ibb1,iblock)
            iter1=(ibb1-1)/nel+1
            ii1  = ibb1-nel*(iter1-1)
            ii1  = nsta+ii1-1
            ip0 = ip0b + ib1
!            if(mod(ip0-1,nrank_e)/=myrank_e) cycle
            if(kimg_t==1) then
               vr=0.d0
               qr=0.d0
            else
               vr=0.d0
               vi=0.d0
               qr=0.d0
               qi=0.d0
            end if
            do ia = 1, natm
               it = ityp(ia)
               do lmt1 = 1, ilmt(it)
                  p = lmta(lmt1,ia)
                  do lmt2 = 1, ilmt(it)
                     s = lmta(lmt2,ia)
                     if(ipaw(it)==0)then
                        facv   = iwei(ia)*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin))
                     else
                        facv   = iwei(ia)*(dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin))
                     endif
                     facq   = iwei(ia)*q(lmt1,lmt2,it)
                     if(kimg==1) then
                        tmpr = fsr(ii1,p,iter1)*fsr(ii2,s,iter2)&
                    &        + fsi(ii1,p,iter1)*fsi(ii2,s,iter2)
                        vr = vr + facv*tmpr
                        qr = qr + facq*tmpr
                     else
                        if(k_symmetry(ik) == GAMMA) then
                           tmpr = fsr(ii1,p,iter1)*fsr(ii2,s,iter2)
                           vr = vr + facv*tmpr
                           qr = qr + facq*tmpr
                        else
                           tmpr = fsr(ii1,p,iter1)*fsr(ii2,s,iter2)&
                    &        + fsi(ii1,p,iter1)*fsi(ii2,s,iter2)
                           tmpi = fsr(ii1,p,iter1)*fsi(ii2,s,iter2)&
                    &        - fsi(ii1,p,iter1)*fsr(ii2,s,iter2)
                           vr = vr + facv*tmpr
                           vi = vi + facv*tmpi
                           qr = qr + facq*tmpr
                           qi = qi + facq*tmpi
                        end if
                     end if
                  end do
               end do
            end do
            if(kimg_t==1) then
               w1hw2(ip0) = w1hw2(ip0) + vr
               w1sw2(ip0) = w1sw2(ip0) + qr
            else
               w1hw2(2*ip0-1) = w1hw2(2*ip0-1) + vr
               w1hw2(2*ip0  ) = w1hw2(2*ip0  ) + vi
               w1sw2(2*ip0-1) = w1sw2(2*ip0-1) + qr
               w1sw2(2*ip0  ) = w1sw2(2*ip0  ) + qi
            end if
         end do
      end do
#endif
    end subroutine add_nonlocal_part

    subroutine subspace_rotation
      integer :: ib1,ib2,ibb2,iadd,is,is1
!!$      real(kind=DP), dimension(np_g1k_x,neg,kimg) :: zaj_wk
      real(kind=DP), allocatable, dimension(:,:,:) :: zaj_wk
      real(kind=DP), allocatable, dimension(:,:,:) :: zah_wk
      real(kind=DP), allocatable, dimension(:,:)   :: fsr_wk
      real(kind=DP), allocatable, dimension(:,:)   :: fsi_wk
#ifdef NEC_TUNE
      integer :: i
      real(kind=DP), allocatable, dimension(:,:)   :: fsr_t, fsi_t

      if(kimg==1) then
         allocate(zaj_wk(maxval(np_g1k),nel,kimg))
         if(itot /= max_iter_mdkosugi) allocate(zah_wk(maxval(np_g1k),nel,kimg))
         allocate(fsr_wk(nel,np_fs),fsi_wk(nel,np_fs))
         allocate(fsr_t(nsize_sb_now,np_fs),fsi_t(nsize_sb_now,np_fs))
         num = 0
         do ibb2 = 1, nsize_max_sb_now
            if(ibover(ibb2,iblock) < 0) cycle
            iter2=(ibb2-1)/nel+1
            ii2=ibb2-nel*(iter2-1)
            ii2=nsta+ii2-1
            num = num + 1
            do ii = ista_fs, iend_fs
               i = ii - ista_fs + 1
               fsr_t(num,i) = fsr(ii2,ii,iter2)
               fsi_t(num,i) = fsi(ii2,ii,iter2)
            end do
         end do
         call dgemm('N','N',np_g1k(ik),nel,nsize_sb_now,1.0d0,zat_t,maxval(np_g1k),vec,nsize_sb_now,0.0d0,zaj_wk,maxval(np_g1k))
         call dgemm('T','N',nel,np_fs,nsize_sb_now,1.0d0,vec,nsize_sb_now,fsr_t,nsize_sb_now,0.0d0,fsr_wk,nel)
         call dgemm('T','N',nel,np_fs,nsize_sb_now,1.0d0,vec,nsize_sb_now,fsi_t,nsize_sb_now,0.0d0,fsi_wk,nel)
         if(itot /= max_iter_mdkosugi) then
            call dgemm('N','N',np_g1k(ik),nel,nsize_sb_now,1.0d0,zah_t,maxval(np_g1k),vec,nsize_sb_now,0.0d0,zah_wk,maxval(np_g1k))
         end if
         zaj_l(:,nsta:nend,ik,:) = zaj_wk(:,:,:)
         if(itot /= max_iter_mdkosugi) zah_l(:,nsta:nend,:,1) = zah_wk(:,:,:)
         do is = 1, np_fs ! MPI
            do ib = nsta, nend ! MPI
               fsr_l(ib,is,ik) = fsr_wk(ib-nsta+1,is)
               fsi_l(ib,is,ik) = fsi_wk(ib-nsta+1,is)
            end do
         end do
         deallocate(fsr_t,fsi_t)
         deallocate(zaj_wk,fsr_wk,fsi_wk)
         if(itot /= max_iter_mdkosugi) deallocate(zah_wk)
      else
         if(k_symmetry(ik) == GAMMA) then
            allocate(zaj_wk(maxval(np_g1k),nel,kimg))
            if(itot /= max_iter_mdkosugi) allocate(zah_wk(maxval(np_g1k),nel,kimg))
            allocate(fsr_wk(nel,np_fs))
            allocate(fsr_t(nsize_sb_now,np_fs))
            num = 0
            do ibb2 = 1, nsize_max_sb_now
               if(ibover(ibb2,iblock) < 0) cycle
               iter2=(ibb2-1)/nel+1
               ii2=ibb2-nel*(iter2-1)
               ii2=nsta+ii2-1
               num = num + 1
               do ii = ista_fs, iend_fs
                  i = ii - ista_fs + 1
                  fsr_t(num,i) = fsr(ii2,ii,iter2)
               end do
            end do
            call dgemm('N','N',np_g1k(ik),nel,nsize_sb_now,1.0d0,zat_t(1,1,1),   maxval(np_g1k),vec,nsize_sb_now,0.0d0,zaj_wk(1,1,1),   maxval(np_g1k))
            call dgemm('N','N',np_g1k(ik),nel,nsize_sb_now,1.0d0,zat_t(1,1,kimg),maxval(np_g1k),vec,nsize_sb_now,0.0d0,zaj_wk(1,1,kimg),maxval(np_g1k))
            call dgemm('T','N',nel,np_fs,nsize_sb_now,1.0d0,vec,nsize_sb_now,fsr_t,nsize_sb_now,0.0d0,fsr_wk,nel)
            if(itot /= max_iter_mdkosugi) then
               call dgemm('N','N',np_g1k(ik),nel,nsize_sb_now,1.0d0,zah_t(1,1,1),   maxval(np_g1k),vec,nsize_sb_now,0.0d0,zah_wk(1,1,1),   maxval(np_g1k))
               call dgemm('N','N',np_g1k(ik),nel,nsize_sb_now,1.0d0,zah_t(1,1,kimg),maxval(np_g1k),vec,nsize_sb_now,0.0d0,zah_wk(1,1,kimg),maxval(np_g1k))
            end if
            zaj_l(:,nsta:nend,ik,:) = zaj_wk(:,:,:)
            if(itot /= max_iter_mdkosugi) zah_l(:,nsta:nend,:,1) = zah_wk(:,:,:)
            do is = 1, np_fs ! MPI
               do ib = nsta, nend ! MPI
                  fsr_l(ib,is,ik) = fsr_wk(ib-nsta+1,is)
               end do
            end do
            deallocate(fsr_t)
            deallocate(zaj_wk,fsr_wk)
            if(itot /= max_iter_mdkosugi) deallocate(zah_wk)
         else
            allocate(zaj_wk(maxval(np_g1k)*kimg,nel,1))
            if(itot /= max_iter_mdkosugi) allocate(zah_wk(maxval(np_g1k)*kimg,nel,1))
            allocate(fsr_wk(nel*kimg,np_fs))
            allocate(fsr_t(nsize_sb_now*kimg,np_fs))
            num = 0
            do ibb2 = 1, nsize_max_sb_now
               if(ibover(ibb2,iblock) < 0) cycle
               iter2=(ibb2-1)/nel+1
               ii2=ibb2-nel*(iter2-1)
               ii2=nsta+ii2-1
               num = num + 1
               do ii = ista_fs, iend_fs
                  i = ii - ista_fs + 1
                  fsr_t(2*num-1,i) = fsr(ii2,ii,iter2)
                  fsr_t(2*num,  i) = fsi(ii2,ii,iter2)
               end do
            end do
            call zgemm('N','N',np_g1k(ik),nel,nsize_sb_now,c1,zat_t,maxval(np_g1k),vec,nsize_sb_now,c0,zaj_wk,maxval(np_g1k))
            call zgemm('T','N',nel,np_fs,nsize_sb_now,c1,vec,nsize_sb_now,fsr_t,nsize_sb_now,c0,fsr_wk,nel)
            if(itot /= max_iter_mdkosugi) then
               call zgemm('N','N',np_g1k(ik),nel,nsize_sb_now,c1,zah_t,maxval(np_g1k),vec,nsize_sb_now,c0,zah_wk,maxval(np_g1k))
            end if
            do ib = nsta, nend
               i = ib - nsta + 1
               do ii = 1, np_g1k(ik)
                  zaj_l(ii,ib,ik,1)    = zaj_wk(2*ii-1,i,1)
                  zaj_l(ii,ib,ik,kimg) = zaj_wk(2*ii,  i,1)
               end do
            end do
            if(itot /= max_iter_mdkosugi) then
               do ib = nsta, nend
                  i = ib - nsta + 1
                  do ii = 1, np_g1k(ik)
                     zah_l(ii,ib,1,   1) = zah_wk(2*ii-1,i,1)
                     zah_l(ii,ib,kimg,1) = zah_wk(2*ii,  i,1)
                  end do
               end do
            end if
            do is = 1, np_fs ! MPI
               do ib = nsta, nend ! MPI
                  i = ib - nsta + 1
                  fsr_l(ib,is,ik) = fsr_wk(2*i-1,is)
                  fsi_l(ib,is,ik) = fsr_wk(2*i,  is)
               end do
            end do
            deallocate(fsr_t)
            deallocate(zaj_wk,fsr_wk)
            if(itot /= max_iter_mdkosugi) deallocate(zah_wk)
         end if
      end if
#else
      allocate(zaj_wk(maxval(np_g1k),nel,kimg))
      if(itot /= max_iter_mdkosugi) then
        allocate(zah_wk(maxval(np_g1k),nel,kimg))
        zah_wk(:,:,:) = 0.d0
      end if
      allocate(fsr_wk(nel,nlmta))
      if(k_symmetry(ik) /= GAMMA) then
        allocate(fsi_wk(nel,nlmta))
        fsi_wk(:,:)=0.d0
      end if

      zaj_wk(:,:,:) = 0.d0
      fsr_wk(:,:)=0.d0
      if(kimg==1) then
         do ib1=1,nel
            do ibb2=1,nsize_max_sb_now
               if(ibover(ibb2,iblock)<0) cycle
               ib2 = ibover(ibb2,iblock)
               iter2=(ibb2-1)/nel+1
               ii2=ibb2-nel*(iter2-1)
               ii2=nsta+ii2-1
               hr2=vec(ib2,ib1)
               do ii=ista_g1k(ik),iend_g1k(ik)
                  iadd = ii - ista_g1k(ik) + 1
                  zaj_wk(iadd,ib1,kimg) = zaj_wk(iadd,ib1,kimg) + zat_l(iadd,ii2,kimg,iter2)*hr2
               end do
               do ii=1,nlmta
                 fsr_wk(ib1,ii) = fsr_wk(ib1,ii) + fsr(ii2,ii,iter2)*hr2
                 fsi_wk(ib1,ii) = fsi_wk(ib1,ii) + fsi(ii2,ii,iter2)*hr2
               end do
               if(itot /= max_iter_mdkosugi) then
                 do ii=ista_g1k(ik),iend_g1k(ik)
                    iadd = ii - ista_g1k(ik) + 1
                    zah_wk(iadd,ib1,kimg) = zah_wk(iadd,ib1,kimg) + zah_l(iadd,ii2,kimg,iter2)*hr2
                 end do
               end if
            end do
         end do
      else
         if(k_symmetry(ik) == GAMMA) then
            do ib1=1,nel
               do ibb2=1,nsize_max_sb_now
                  if(ibover(ibb2,iblock)<0) cycle
                  ib2 = ibover(ibb2,iblock)
                  iter2=(ibb2-1)/nel+1
                  ii2=ibb2-nel*(iter2-1)
                  ii2=nsta+ii2-1
                  hr2=vec(ib2,ib1)
                  do ii=ista_g1k(ik),iend_g1k(ik)
                     iadd = ii - ista_g1k(ik) + 1
                     dr1=zat_l(iadd,ii2,1   ,iter2)
                     di1=zat_l(iadd,ii2,kimg,iter2)
                     zaj_wk(iadd,ib1,1   ) = zaj_wk(iadd,ib1,1   ) + dr1*hr2
                     zaj_wk(iadd,ib1,kimg) = zaj_wk(iadd,ib1,kimg) + di1*hr2
                  end do
                  do ii=1,nlmta
                    fsr_wk(ib1,ii) = fsr_wk(ib1,ii) + fsr(ii2,ii,iter2)*hr2
                  end do
                  if(itot /= max_iter_mdkosugi) then
                    do ii=ista_g1k(ik),iend_g1k(ik)
                       iadd = ii - ista_g1k(ik) + 1
                       dr1=zah_l(iadd,ii2,1   ,iter2)
                       di1=zah_l(iadd,ii2,kimg,iter2)
                       zah_wk(iadd,ib1,1   ) = zah_wk(iadd,ib1,1   ) + dr1*hr2
                       zah_wk(iadd,ib1,kimg) = zah_wk(iadd,ib1,kimg) + di1*hr2
                    end do
                  end if
               end do
            end do
         else
            do ib1=1,nel
               do ibb2=1,nsize_max_sb_now
                  if(ibover(ibb2,iblock)<0) cycle
                  ib2 = ibover(ibb2,iblock)
                  iter2=(ibb2-1)/nel+1
                  ii2=ibb2-nel*(iter2-1)
                  ii2=nsta+ii2-1
                  hr2=vec(2*ib2-1,ib1)
                  hi2=vec(2*ib2  ,ib1)
                  do ii=ista_g1k(ik),iend_g1k(ik)
                     iadd = ii - ista_g1k(ik) + 1
                     dr1=zat_l(iadd,ii2,1   ,iter2)
                     di1=zat_l(iadd,ii2,kimg,iter2)
                     zaj_wk(iadd,ib1,1   ) = zaj_wk(iadd,ib1,1   ) + dr1*hr2 - di1*hi2
                     zaj_wk(iadd,ib1,kimg) = zaj_wk(iadd,ib1,kimg) + dr1*hi2 + di1*hr2
                  end do
                  do ii=1,nlmta
                     dr1=fsr(ii2,ii,iter2)
                     di1=fsi(ii2,ii,iter2)
                     fsr_wk(ib1,ii) = fsr_wk(ib1,ii) + dr1*hr2 - di1*hi2
                     fsi_wk(ib1,ii) = fsi_wk(ib1,ii) + dr1*hi2 + di1*hr2
                  end do
                  if(itot /= max_iter_mdkosugi) then
                    do ii=ista_g1k(ik),iend_g1k(ik)
                       iadd = ii - ista_g1k(ik) + 1
                       dr1=zah_l(iadd,ii2,1   ,iter2)
                       di1=zah_l(iadd,ii2,kimg,iter2)
                       zah_wk(iadd,ib1,1   ) = zah_wk(iadd,ib1,1   ) + dr1*hr2 - di1*hi2
                       zah_wk(iadd,ib1,kimg) = zah_wk(iadd,ib1,kimg) + dr1*hi2 + di1*hr2
                    end do
                  end if
               end do
            end do
         end if
      end if
!print *,itot,  itot /= max_iter_mdkosugi
      zaj_l(:,nsta:nend,ik,:) = zaj_wk(:,:,:)
      if(itot /= max_iter_mdkosugi) zah_l(:,nsta:nend,:,1) = zah_wk(:,:,:)
!     fsr_l(nsta:nend,:,ik) = fsr_wk(:,:)
      do is = 1, np_fs ! MPI
         is1=nis_fs(myrank_g)+is-1
         do ib = nsta, nend ! MPI
             fsr_l(ib,is,ik) = fsr_wk(ib-nsta+1,is1)
         end do
      end do
!     if(k_symmetry(ik) /= GAMMA) fsi_l(nsta:nend,:,ik) = fsi_wk(:,:)
      if(k_symmetry(ik) /= GAMMA) then
         do is = 1, np_fs ! MPI
            is1=nis_fs(myrank_g)+is-1
            do ib = nsta, nend ! MPI
               fsi_l(ib,is,ik) = fsi_wk(ib-nsta+1,is1)
            end do
         end do
      end if

      deallocate(zaj_wk,fsr_wk)
      if(itot /= max_iter_mdkosugi) deallocate(zah_wk)
      if(k_symmetry(ik) /= GAMMA) deallocate(fsi_wk)
#endif
    end subroutine subspace_rotation

  end subroutine evolve_WFs_in_subspace_ksg_3D


  subroutine decide_correction_vector_ksg_3D(precon,ik,ekin_l,afft_l,bfft_l, &
                                             wk_bfft_l,lsize,ibsize,isrsize,fft_l_size,p_l)
    integer, intent(in)       :: precon, ik, lsize, ibsize, isrsize, fft_l_size
    real(kind=DP), intent(in)  :: ekin_l(maxval(np_g1k))
    real(kind=DP), intent(in)  :: afft_l(lsize*kimg)
    real(kind=DP), intent(out) :: bfft_l(lsize*kimg,1)
    real(kind=DP), intent(inout) :: wk_bfft_l(lsize*kimg,ibsize)
    real(kind=DP)              :: p_l(maxval(np_g1k))
    integer :: ib
    integer :: id_sname = -1
    call tstatc0_begin('decide_correction_vector ', id_sname,1)

    do ib = 1, np_e ! MPI
#ifdef __TIMER_COMM__
       call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l,0)
#else
       call m_ES_WF_in_Rspace_3D(ik,ib,ib,ibsize,lsize,wk_bfft_l)
#endif
       call m_FFT_Vlocal_W_3D(afft_l,wk_bfft_l,lsize,ibsize,nel_fft_y(myrank_g))
       call m_FFT_Direct_3D(nfout,wk_bfft_l,lsize,ibsize)
       call map_fft_to_WF_3D(ik,lsize,ibsize,wk_bfft_l,bfft_l,isrsize,fft_l_size)
       call SD_direction_3D(precon,ik,ib,ekin_l,bfft_l,p_l,lsize) !-here
    end do

!    call orthogonalize_SD_drctns_ksg(ik,to=OTHER_BANDS)  ! -(m_ES_WF_by_SDorCG) ->(wfsd_l, bsd(ri)_l)
    call orthogonalize_SD_drctns_ksg(ik,to=OTHER_BANDS)
!    call orthogonalize_SD_drctns(ik,to=ALL_BANDS)  ! -(m_ES_WF_by_SDorCG) ->(wfsd_l, bsd(ri)_l)

    call normalize_wfsd_3D(ik)

    call tstatc0_end(id_sname)
  end subroutine decide_correction_vector_ksg_3D

  subroutine orthogonalize_SD_drctns_ksg(ik,to)
    integer, intent(in) :: ik,to

    integer :: itmp
    integer :: id_sname = -1
    real(kind=DP), allocatable, dimension(:,:) :: fs_mpi,fs_mpi2
    integer :: is, is1, ib
    call tstatc0_begin('orthogonalize_SD_drctns in Modified Davidson ', id_sname)

!    if(modnrm == EXECUT) call m_ES_betar_dot_Psi_4_each_k(wfsd_l,ik,ik,ik,bsdr_l,bsdi_l)
    !                                        ->bsd(ri)_l

    zat_l(:,:,:,1) = zaj_l(:,:,ik,:)
    allocate(fs_mpi(np_e,nlmta))
    allocate(fs_mpi2(np_e,nlmta))

    fs_mpi=0.d0
    do is = 1, np_fs ! MPI
       is1=nis_fs(myrank_g)+is-1
       do ib = 1, np_e ! MPI
          fs_mpi(ib,is1) = fsr_l(ib,is,ik)
       end do
    end do
    fs_mpi2=0.d0
    call mpi_allreduce(fs_mpi,fs_mpi2,np_e*nlmta &
      & ,mpi_double_precision,mpi_sum &
      & ,mpi_ke_world,ierr)       ! MPI
    fsr(1:np_e,1:nlmta,1) = fs_mpi2(1:np_e,1:nlmta)
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       fs_mpi=0.d0
       do is = 1, np_fs ! MPI
          is1=nis_fs(myrank_g)+is-1
          do ib = 1, np_e ! MPI
             fs_mpi(ib,is1) = fsi_l(ib,is,ik)
          end do
       end do
       fs_mpi2=0.d0
       call mpi_allreduce(fs_mpi,fs_mpi2,np_e*nlmta &
         & ,mpi_double_precision,mpi_sum &
         & ,mpi_ke_world,ierr)       ! MPI
       fsi(1:np_e,1:nlmta,1) = fs_mpi2(1:np_e,1:nlmta)
    end if
    deallocate(fs_mpi)
    deallocate(fs_mpi2)

    zaj_l(:,:,ik,:) = zajold_l(:,:,:)
    fsr_l(:,:,ik) = fsrold_l(:,:)
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       fsi_l(:,:,ik) = fsiold_l(:,:)
    end if

    itmp=modnrm
    modnrm=EXECUT
    call m_ES_betar_dot_Psi_4_each_k_3D(wfsd_l,ik,ik,ik,bsdr_l,bsdi_l)
    call m_ES_orthogonalize_SD_to_WFs_3D(ik,to,wfsd_l,bsdr_l,bsdi_l) ! ->(wfsd_l,bsd(ri)_l)
    modnrm=itmp
!    call m_ES_betar_dot_Psi_4_each_k(wfsd_l,ik,ik,ik,bsdr_l,bsdi_l)
    call tstatc0_end(id_sname)
  end subroutine orthogonalize_SD_drctns_ksg



#endif
end module m_ES_WF_by_ModifiedDavidson

