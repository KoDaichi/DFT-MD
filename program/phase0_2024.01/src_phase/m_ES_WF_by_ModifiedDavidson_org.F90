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
  use m_FFT,                 only : nfft,fft_box_size_WF, m_FFT_Vlocal_W, m_FFT_WF
  use m_Kpoints,             only : kv3,vkxyz, k_symmetry
  use m_PlaneWaveBasisSet,   only : kg1, iba, igf, nbase, m_pwBS_kinetic_energies
!!$  use m_Electronic_Structure,only : zaj_l, neordr, sc, nrvf_ordr, eko_l, vlhxcQ &
  use m_Electronic_Structure,only : zaj_l, neordr, nrvf_ordr, eko_l, vlhxcQ &
       &                          , occup_l &
       &                          , fsr_l,fsi_l, vnlph_l, vlhxc_l &
       &                          , m_ES_Vlocal_in_Rspace &
       &                          , m_ES_WF_in_Rspace &
       &                          , m_ES_wd_zaj_small_portion &
       &                          , m_ES_wd_eko &
       &                          , m_ES_sort_eigen_values
  use m_ES_ortho,           only : np_g1k_x                                          &
       &                         , m_ES_W_transpose                                  &
       &                         , m_ES_W_transpose_back                             &
       &                         , m_ES_orthogonalize_SD_to_WFs
  use m_ES_nonlocal,        only : m_ES_Vnonlocal_W                                  &
       &                         , m_ES_betar_dot_WFs_4_each_k                       &
       &                         , m_ES_alloc_scss_etc                               &
       &                         , m_ES_dealloc_scss_etc                             &
       &                         , m_ES_betar_dot_Psi_4_each_k
  use m_Ionic_System,       only : natm, iwei, ityp, ntyp
  use m_PseudoPotential,    only : ilmt,nlmta,lmta,q,dion &
       &                         , lmtt,ltp,mtp &
       &                         , m_PP_include_vanderbilt_pot &
       &                         , ipaw,dion_paw,modnrm
  use m_NonLocal_Potential, only : snl

! ============================== added by K. Tagami ================== 11.0
  use m_Const_Parameters,   only : CMPLDP
  use m_Control_Parameters,  only : ndim_spinor, noncol, ndim_chgpot
  use m_PseudoPotential,     only : q_noncl, dion_scr_noncl
  use m_Electronic_Structure,  only : m_ES_Vlocal_in_Rspace_noncl

  use m_FFT,                 only : m_FFT_Vlocal_W_noncl
  use m_ES_ortho,              only : m_ES_orthogonl_SD_to_WFs_noncl
  use m_Electronic_Structure,      only : m_ES_sort_eigen_vals_noncl
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

! =========================== added by K. Tagami ============== 11.0
  public :: m_ESmddavid_Renew_WF_noncl
  public :: m_ESmddavid_Subspace_Rot_noncl
  public :: m_ESmdkosugi_Renew_WF_noncl
! ============================================================= 11.0

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
    logical :: frestart
    integer :: max_itot
    integer :: iblock_now, itot_now, ipri0
    integer :: n_unconv

    max_itot = 1

    allocate(ekin(kg1),p(kg1))
    call m_ES_alloc_scss_etc()
    allocate(afft(nfft)); allocate(bfft(nfft))

    call allocate_matrix
    do ispin = 1, nspin, (af+1)
       call m_ES_Vlocal_in_Rspace(ispin,afft)      ! (ptfft1) vlhxc_l->afft
       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) /= myrank_k) cycle          ! MPI
          iksnl = (ik-1)/nspin + 1
          call allocate_t_matrix(ik)
          call m_pwBS_kinetic_energies(ik,vkxyz,ekin) ! (diakin) ->ekin
          max_itot=0
          feigconv = .false.
          Loop: do itot=1,max_iter_mddavid
             itot_now = itot
             if(itot>max_itot) max_itot=itot
             call m_ES_Vnonlocal_W(ik,iksnl,ispin,switch_of_eko_part=ON) ! -> vnlph_l
             call decide_correction_vector(precon,ik,ekin,afft,bfft,p)   ! -> wfsd_l
             call prepare_Hloc_phi(ik,ekin,afft,bfft,itot)
!print *,ik,itot
if(iprimddavidson>=2) write(nfout,*) 'B',feigconv
             Block_Loop: do iblock=1,nblock
                iblock_now = iblock
                call evolve_WFs_in_subspace &
                          (ik,ispin,ekin,afft,bfft,iblock,itot,frestart)
!                !!$write(nfout,*) 'debug loop itot,idavid=',itot,idavid
!                if(ik==1.and.iprimddavidson>= 2) &
!                  & call m_ES_wd_zaj_small_portion(nfout,ik," -- after md davidson --",21)
!                if(frestart) exit David_Loop
!                if(eigenvalues_are_converged(n_unconv)) exit Loop
             end do Block_Loop
if(iprimddavidson>=2) write(nfout,*) 'A',feigconv
!             call m_ES_betar_dot_WFs_4_each_k(nfout,ik)   ! -> fsr_l,fsi_l
          end do Loop
          if(iprimddavidson>=2) then
             write(nfout,'("MdDavidson: ik=",i5," itot=",i5," subspace=",i5)') &
                                                                 ik, itot_now, nsize_sb_now
          end if
          call deallocate_t_matrix
       enddo      ! k-point loop
    enddo      ! spin loop
    call deallocate_matrix
    if(iprimddavidson>=2) then
       write(nfout,'("Modified Davidson: max_itot=",i5)') max_itot
    end if
!
    call m_ES_sort_eigen_values()
!!!  ( in case of af=1 )
!    if(af /= 0) then
!       call cp_eigen_values_for_af       !-(contained here)
!       call expand_neordr_and_nrvf_ordr  !-(contained here)
!    end if

    call get_ipri0(iprimddavidson,ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko(nfout,mode=SCF)
!
    deallocate(bfft);   deallocate(afft)
    call m_ES_dealloc_scss_etc()
    deallocate(ekin,p)
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

! ================================= added by K. Tagami ================ 11.0
  subroutine m_ESmddavid_Renew_WF_noncl(nfout,precon)
    integer, intent(in) :: nfout,precon
    integer             :: ispin, ik, iksnl, switch_of_eko_part
    integer :: iblock,itot

    real(kind=DP), allocatable ::  afft_kt(:,:)
    real(kind=DP), allocatable ::  bfft_kt(:,:)
    real(kind=DP), allocatable :: vnlph_noncl(:,:,:,:)

    real(kind=DP), allocatable, dimension(:) :: ekin,p

    logical :: frestart
    integer :: max_itot
    integer :: iblock_now, itot_now, ipri0
    integer :: n_unconv

    integer :: is1, is2, istmp, k2

    max_itot = 1

    allocate(ekin(kg1),p(kg1))

    call m_ES_alloc_scss_etc()
    allocate(afft_kt(nfft,ndim_chgpot)); afft_kt = 0.0d0
    allocate(bfft_kt(nfft,ndim_spinor)); bfft_kt = 0.0d0
    allocate(vnlph_noncl(kg1,np_e,kimg,ndim_spinor)); vnlph_noncl = 0.0d0

    call allocate_matrix_noncl
    call m_ES_Vlocal_in_Rspace_noncl( afft_kt )    ! (ptfft1) vlhxc_l->afft

    Do ik=1, kv3, ndim_spinor
       if ( map_k(ik) /= myrank_k )  cycle          ! MPI
       iksnl = (ik-1)/ndim_spinor + 1

       call allocate_t_matrix_noncl( ik, ik +ndim_spinor -1 ) ! -> np_g1k_x
       call m_pwBS_kinetic_energies( ik, vkxyz, ekin ) ! (diakin) ->ekin


       max_itot=0
       feigconv = .false.

       Loop: do itot=1,max_iter_mddavid
          itot_now = itot
          if (itot>max_itot) max_itot=itot

          vnlph_noncl = 0.0d0
          Do is1=1, ndim_spinor
             Do is2=1, ndim_spinor
                istmp = ( is1-1 )*ndim_spinor +is2
                k2 = ik +is2 -1

                if ( precon==ON ) then
                   call m_ES_Vnonlocal_W( k2, iksnl, istmp, switch_of_eko_part=OFF )
                   ! -> vnlph_l
                else
                   call m_ES_Vnonlocal_W( k2, iksnl, istmp, switch_of_eko_part=ON )
                   ! -> vnlph_l
                endif
                vnlph_noncl(:,:,:,is1) = vnlph_noncl(:,:,:,is1) &
                     &                     + vnlph_l(:,:,:)
             End do
          End do

          call decide_correction_vector_noncl( precon, ik, ekin, &
               &                               afft_kt, bfft_kt, vnlph_noncl, p )
                                                              ! -> wfsd_l

          call prepare_Hloc_phi_noncl( ik, ekin, afft_kt, bfft_kt, itot )

          Block_Loop: do iblock=1,nblock
             iblock_now = iblock
             call evolve_WFs_in_subspace_noncl( ik, ekin, afft_kt, bfft_kt, &
                  &                             iblock, itot, frestart )

             if (ik==1.and.iprimddavidson>= 2) then
                Do is1=1, ndim_spinor
                   call m_ES_wd_zaj_small_portion( nfout, ik+is1-1, &
                        &                          " -- after md davidson --",21 )
                End do
             endif

!             if(frestart) exit David_Loop
!            if(eigenvalues_are_converged(n_unconv)) exit Loop

          end do Block_Loop

! -
          Do is1=1, ndim_spinor
             call m_ES_betar_dot_WFs_4_each_k( nfout, ik+is1-1 )   ! -> fsr_l,fsi_l
          End do

       end do Loop

       if(iprimddavidson>=2) then
          write(nfout,'("MdDavidson: ik=",i5," itot=",i5," subspace=",i5)') &
               &               ik, itot_now, nsize_sb_now
       end if

       call deallocate_t_matrix
    enddo      ! k-point loop

    call deallocate_matrix_noncl

    if (iprimddavidson>=2) then
       write(nfout,'("Modified Davidson: max_itot=",i5)') max_itot
    end if
!
    call m_ES_sort_eigen_vals_noncl()


    call get_ipri0(iprimddavidson,ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko(nfout,mode=SCF)
!
    deallocate(bfft_kt);   deallocate(afft_kt)
    call m_ES_dealloc_scss_etc()
    deallocate(ekin,p)

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

  end subroutine m_ESmddavid_Renew_WF_noncl
! ================================================================= 11.0

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

  subroutine allocate_t_matrix(ik)
    integer, intent(in) :: ik
    integer :: kimg_t
    if(k_symmetry(ik) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = kimg
    end if
    allocate(zat_l(kg1,np_e,kimg,2)) ! MPI
    allocate(zah_l(kg1,np_e,kimg,2)) ! MPI
    allocate(w1hw2(msize_matrix*kimg_t))
    allocate(w1sw2(msize_matrix*kimg_t))
    if(npes>1) then
       allocate(w1hw2_mpi(msize_matrix*kimg_t))
       allocate(w1sw2_mpi(msize_matrix*kimg_t))
    end if
    allocate(wfsd_l(kg1,np_e,ik:ik,kimg)); wfsd_l = 0.d0
    allocate(bsdr_l(np_e,nlmta,ik:ik)); bsdr_l = 0.d0
    allocate(bsdi_l(np_e,nlmta,ik:ik)); bsdi_l = 0.d0
  end subroutine allocate_t_matrix

! ====================== added by K. Tagami ================== 11.0
  subroutine allocate_t_matrix_noncl( ik1, ik2 )
    integer, intent(in) :: ik1, ik2
    integer :: kimg_t

    if(k_symmetry(ik1) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = kimg
    end if

    allocate(zat_l_noncl(kg1,np_e,kimg,2,ndim_spinor)) ! MPI
    allocate(zah_l_noncl(kg1,np_e,kimg,2,ndim_spinor)) ! MPI

    allocate(w1hw2(msize_matrix*kimg_t))
    allocate(w1sw2(msize_matrix*kimg_t))
    if(npes>1) then
       allocate(w1hw2_mpi(msize_matrix*kimg_t))
       allocate(w1sw2_mpi(msize_matrix*kimg_t))
    end if

    allocate(wfsd_l(kg1,np_e,ik1:ik2,kimg)); wfsd_l = 0.d0
    allocate(bsdr_l(np_e,nlmta,ik1:ik2)); bsdr_l = 0.d0
    allocate(bsdi_l(np_e,nlmta,ik1:ik2)); bsdi_l = 0.d0
  end subroutine allocate_t_matrix_noncl
! ============================================================ 11.0


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

  subroutine decide_correction_vector(precon,ik,ekin,afft,bfft,p)
    integer, intent(in)       :: precon, ik
    real(kind=DP), intent(in)  :: ekin(kg1)
    real(kind=DP), intent(in)  :: afft(nfft)
    real(kind=DP), intent(out) :: bfft(nfft)
    real(kind=DP)              :: p(kg1)

    integer :: ib
    integer :: id_sname = -1
    call tstatc0_begin('decide_correction_vector ', id_sname,1)

    do ib = ista_e, iend_e, istep_e      ! MPI
       call m_ES_WF_in_Rspace(ik,ib,bfft)! (swffft)
       call m_FFT_Vlocal_W(afft,bfft)    ! (afft, bfft) -> (bfft)
       call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON) ! bfft: G-space repres.
       call SD_direction(precon,ik,ib,ekin,bfft,p) !-here
    end do

    call orthogonalize_SD_drctns(ik,to=OTHER_BANDS)  ! -(m_ES_WF_by_SDorCG) ->(wfsd_l, bsd(ri)_l)
!    call orthogonalize_SD_drctns(ik,to=ALL_BANDS)  ! -(m_ES_WF_by_SDorCG) ->(wfsd_l, bsd(ri)_l)

    call normalize_wfsd(ik)

    call tstatc0_end(id_sname)
  end subroutine decide_correction_vector

! =========================== added by K. Tagami ============== 11.0
  subroutine decide_correction_vector_noncl( precon, ik, ekin, &
       &                                     afft_kt, bfft_kt, vnlph_noncl, p )
    integer, intent(in)       :: precon, ik
    real(kind=DP), intent(in)  :: ekin(kg1)
    real(kind=DP), intent(inout)  :: afft_kt(nfft,ndim_chgpot)
    real(kind=DP), intent(out) :: bfft_kt(nfft,ndim_spinor)
    real(kind=DP), intent(in), dimension(kg1,np_e,kimg,ndim_spinor) :: vnlph_noncl

    real(kind=DP)              :: p(kg1)

    integer :: ib, is
    integer :: id_sname = -1

    call tstatc0_begin('decide_correction_vector_noncl ', id_sname,1)

    do ib = ista_e, iend_e, istep_e      ! MPI

       Do is=1, ndim_spinor
          call m_ES_WF_in_Rspace( ik+is-1, ib, bfft_kt(:,is) ) ! (swffft)
       End do
       call m_FFT_Vlocal_W_noncl( afft_kt, bfft_kt, ndim_chgpot, ndim_spinor )
                                           ! (afft, bfft)-> (bfft)
       Do is=1, ndim_spinor
         call m_FFT_WF( ELECTRON, nfout, bfft_kt(:,is), DIRECT, ON )
       End do

       call SD_direction_noncl( precon, ik, ib, ekin, bfft_kt, &
            &                   vnlph_noncl, p ) !-here
    end do

    call orthogonalize_SD_drctns_noncl(ik,to=OTHER_BANDS)  ! -(m_ES_WF_by_SDorCG) ->(wfsd_l, bsd(ri)_l)
!    call orthogonalize_SD_drctns(ik,to=ALL_BANDS)  ! -(m_ES_WF_by_SDorCG) ->(wfsd_l, bsd(ri)_l)

    call normalize_wfsd_noncl(ik)

    call tstatc0_end(id_sname)
  end subroutine decide_correction_vector_noncl
! ============================================================== 11.0

  subroutine SD_direction(precon,ik,ibo,ekin,VlocalW,p)
    integer     , intent(in)                   :: precon,ik,ibo
    real(kind=DP), intent(in), dimension(kg1)  :: ekin
    real(kind=DP), intent(in), dimension(nfft) :: VlocalW
    real(kind=DP)             , dimension(kg1)  :: p

    integer       :: i, i1, ib
    real(kind=DP) :: devr,denom, e1, devi, norm

    ib = map_z(ibo)                                  ! MPI
    denom = 1.d0/product(fft_box_size_WF(1:3,1))
    norm = 0.d0

    if(kimg == 1) then
       do i = 1, iba(ik)
          i1    = igf(nbase(i,ik))
          devr  = (ekin(i)-eko_l(ib,ik))*zaj_l(i,ib,ik,1)&
               & + VlocalW(i1)*denom + vnlph_l(i,ib,1)
          wfsd_l(i,ib,ik,1) = - devr
          norm = norm + devr*devr
       end do
    else if(kimg == 2) then
       do i = 1, iba(ik)
          i1    = igf(nbase(i,ik))
          e1    = ekin(i) - eko_l(ib,ik)
          devr  = e1*zaj_l(i,ib,ik,1) + VlocalW(2*i1-1)*denom+vnlph_l(i,ib,1)
          devi  = e1*zaj_l(i,ib,ik,2) + VlocalW(2*i1  )*denom+vnlph_l(i,ib,2)
          wfsd_l(i,ib,ik,1) = - devr
          wfsd_l(i,ib,ik,2) = - devi
          norm = norm + devr*devr + devi*devi
       end do
       if(k_symmetry(ik) == GAMMA) then
          devr=wfsd_l(1,ib,ik,1)
          devi=wfsd_l(1,ib,ik,2)
          norm = norm*2.d0 - devr*devr - devi*devi
       end if
    end if

    feigconv(ib)=.false.
    if(sqrt(norm) .lt. eps_residual) feigconv(ib)=.true.

    if(precon==ON) then
      call decide_precon_factor_wfsd(ik,ibo,ekin,p)
!call decide_precon_factor_david(ik,hdiag,sdiag,eko_l(ib,ik),p)
      if(kimg == 1) then
         do i = 1, iba(ik)
            wfsd_l(i,ib,ik,1) = p(i)*wfsd_l(i,ib,ik,1)
         end do
      else if(kimg == 2) then
         do i = 1, iba(ik)
            wfsd_l(i,ib,ik,1) = p(i)*wfsd_l(i,ib,ik,1)
            wfsd_l(i,ib,ik,2) = p(i)*wfsd_l(i,ib,ik,2)
         end do
      end if
    end if
  end subroutine SD_direction

! ============================= added by K. Tagami ================ 11.0
  subroutine SD_direction_noncl( precon, ik, ibo, ekin, VlocalW_noncl, &
       &                         vnlph_noncl, p )
    integer     , intent(in)                   :: precon,ik,ibo
    real(kind=DP), intent(in), dimension(kg1)  :: ekin
    real(kind=DP), intent(in), dimension(nfft,ndim_spinor) :: VlocalW_noncl
    real(kind=DP), intent(in), dimension(kg1,np_e,kimg,ndim_spinor) :: vnlph_noncl

    real(kind=DP)             , dimension(kg1)  :: p

    integer       :: i, i1, ib, is, k1
    real(kind=DP) :: devr,denom, e1, devi, norm, ctmp

    ib = map_z(ibo)                                  ! MPI
    denom = 1.d0/product(fft_box_size_WF(1:3,1))
    norm = 0.d0

    if(kimg == 1) then
       Do is=1, ndim_spinor
          k1 = ik + is -1
          do i = 1, iba(ik)
             i1    = igf(nbase(i,ik))
             devr  = (ekin(i)-eko_l(ib,ik))*zaj_l(i,ib,k1,1)&
                  & + VlocalW_noncl(i1,is)*denom + vnlph_noncl(i,ib,1,is)
             wfsd_l(i,ib,k1,1) = - devr
             norm = norm + devr*devr
          end do
       End do
    else if(kimg == 2) then
       Do is=1, ndim_spinor
          k1 = ik + is -1
          do i = 1, iba(ik)
             i1    = igf(nbase(i,ik))
             e1    = ekin(i) - eko_l(ib,ik)
             devr  = e1*zaj_l(i,ib,k1,1) + VlocalW_noncl(2*i1-1,is)*denom &
                  &  +vnlph_noncl(i,ib,1,is)
             devi  = e1*zaj_l(i,ib,k1,2) + VlocalW_noncl(2*i1  ,is)*denom &
                  &  +vnlph_noncl(i,ib,2,is)
             wfsd_l(i,ib,k1,1) = - devr
             wfsd_l(i,ib,k1,2) = - devi
             norm = norm + devr*devr + devi*devi
          end do
       End do

       if(k_symmetry(ik) == GAMMA) then
          ctmp = 0.0d0
          Do is=1, ndim_spinor
             k1 = ik + is -1
             devr = wfsd_l(1,ib,k1,1)
             devi = wfsd_l(1,ib,k1,2)
             ctmp = ctmp + devr*devr + devi*devi
          End do
          norm = norm*2.d0 - ctmp
       end if
    end if

    feigconv(ib)=.false.
    if(sqrt(norm) .lt. eps_residual) feigconv(ib)=.true.

    if(precon==ON) then
      call decide_precon_factor_wfsd(ik,ibo,ekin,p)
!call decide_precon_factor_david(ik,hdiag,sdiag,eko_l(ib,ik),p)
      if(kimg == 1) then
         Do is=1, ndim_spinor
            k1 = ik + is -1
            do i = 1, iba(ik)
               wfsd_l(i,ib,k1,1) = p(i)*wfsd_l(i,ib,k1,1)
            end do
         End do
      else if(kimg == 2) then
         Do is=1, ndim_spinor
            k1 = ik + is -1
            do i = 1, iba(ik)
               wfsd_l(i,ib,k1,1) = p(i)*wfsd_l(i,ib,k1,1)
               wfsd_l(i,ib,k1,2) = p(i)*wfsd_l(i,ib,k1,2)
            end do
         End do
      end if
    end if
  end subroutine SD_direction_noncl
! ============================================================= 11.0

  subroutine decide_precon_factor_wfsd(ik,ibo,ekin,p)
    integer, intent(in)                         :: ik,ibo
    real(kind=DP), intent(in),  dimension(kg1)  :: ekin
    real(kind=DP), intent(out), dimension(kg1)  :: p

    integer       :: i
    real(kind=DP) :: ektot, x, x1, x2, d_ektot

!    call kinetic_energy_wfsd(ik,ibo,ekin,ektot)   ! -here

! ====================== modified by K. Tagami =============== 11.0
!    call kinetic_energy(ik,ibo,ekin,ektot)   ! -here

    if ( noncol ) then
       call kinetic_energy_noncl(ik,ibo,ekin,ektot)   ! -here
    else
       call kinetic_energy(ik,ibo,ekin,ektot)   ! -here
    end if
! ============================================================ 11.0

    d_ektot = 4.d0/ektot/3.d0
    p = 0.d0
    do i = 1, iba(ik)
       x = ekin(i)*d_ektot
       x1 = (x*x+9.d0)*(x+3.d0)
       x2 = (x*x)*(x*x)
       p(i)  = x1/(x1 + x2 )
    end do
!    p=p*d_ektot
  end subroutine decide_precon_factor_wfsd

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

  subroutine kinetic_energy(ik,ibo,dekin,ektot)
    integer, intent(in) :: ik, ibo
    real(kind=DP), intent(in), dimension(kg1) :: dekin
    real(kind=DP), intent(out)                 :: ektot
    integer  :: i, ib
    ektot = 0.d0
    ib = map_z(ibo)
    if(kimg == 1) then
       do i = 1, iba(ik)
          ektot = ektot + dekin(i)*zaj_l(i,ib,ik,1)**2
       end do
    else
       do i = 1, iba(ik)
          ektot = ektot + dekin(i)*( zaj_l(i,ib,ik,1)**2 &
               &                   + zaj_l(i,ib,ik,2)**2)
       end do
    end if
    if(k_symmetry(ik) == GAMMA)  ektot = ektot*2.d0
  end subroutine kinetic_energy

! ======================= added by K. Tagami =============== 11.0
  subroutine kinetic_energy_noncl( ik, ibo, dekin, ektot )
    integer, intent(in) :: ik, ibo
    real(kind=DP), intent(in), dimension(kg1) :: dekin
    real(kind=DP), intent(out)                 :: ektot
    integer  :: i, ib, is, k1

    ektot = 0.d0
    ib = map_z(ibo)
    if(kimg == 1) then
       Do is=1, ndim_spinor
          k1 = ik + is -1
          do i = 1, iba(ik)
             ektot = ektot + dekin(i)*zaj_l(i,ib,k1,1)**2
          end do
       End do
    else
       Do is=1, ndim_spinor
          k1 = ik + is -1
          do i = 1, iba(ik)
             ektot = ektot + dekin(i)*( zaj_l(i,ib,k1,1)**2 &
                  &                   + zaj_l(i,ib,k1,2)**2)
          end do
       End do
    end if
    if(k_symmetry(ik) == GAMMA)  ektot = ektot*2.d0
  end subroutine kinetic_energy_noncl
! ======================================================== 11.0


  subroutine orthogonalize_SD_drctns(ik,to)
    integer, intent(in) :: ik,to

    integer :: itmp
    integer :: id_sname = -1
    call tstatc0_begin('orthogonalize_SD_drctns in Modified Davidson ', id_sname)

!    if(modnrm == EXECUT) call m_ES_betar_dot_Psi_4_each_k(wfsd_l,ik,ik,ik,bsdr_l,bsdi_l)
    !                                                          ->bsd(ri)_l
    itmp=modnrm
    modnrm=EXECUT
    call m_ES_betar_dot_Psi_4_each_k(wfsd_l,ik,ik,ik,bsdr_l,bsdi_l)
    call m_ES_orthogonalize_SD_to_WFs(ik,to,wfsd_l,bsdr_l,bsdi_l) ! ->(wfsd_l,bsd(ri)_l)
    modnrm=itmp
!    call m_ES_betar_dot_Psi_4_each_k(wfsd_l,ik,ik,ik,bsdr_l,bsdi_l)
    call tstatc0_end(id_sname)
  end subroutine orthogonalize_SD_drctns

! ===================== added by K. Tagami =================== 11.0
  subroutine orthogonalize_SD_drctns_noncl(ik,to)
    integer, intent(in) :: ik,to

    integer :: is, itmp
    integer :: id_sname = -1

    call tstatc0_begin('orthogonalize_SD_drctns_noncl in Modified Davidson ', id_sname)

    itmp = modnrm;  modnrm = EXECUT
    Do is=1, ndim_spinor
       call m_ES_betar_dot_Psi_4_each_k( wfsd_l, ik, ik+ndim_spinor-1, &
            &                            ik +is -1, bsdr_l, bsdi_l )
                                        !           ->bsd(ri)_l
    End do

    call m_ES_orthogonl_SD_to_WFs_noncl( ik, to, ik, ik+ndim_spinor-1, &
        &                                wfsd_l, bsdr_l, bsdi_l)
                                             ! ->(wfsd_l,bsd(ri)_l)
    modnrm=itmp

    call tstatc0_end(id_sname)

  end subroutine orthogonalize_SD_drctns_noncl
! ================================================================== 11.0

  subroutine normalize_wfsd(ik)
    integer,intent(in) :: ik
    real(kind=DP) :: norm, wfsdr, wfsdi
    integer :: ib1,ii,ib

    do ib1 = ista_e, iend_e, istep_e     ! MPI
       ib=map_z(ib1)
       norm = 0.d0
       if(kimg==1) then
          do ii=1,iba(ik)
             wfsdr = wfsd_l(ii,ib,ik,kimg)
             norm = norm + wfsdr*wfsdr
          end do
          norm = 1.d0/sqrt(norm)
          do ii=1,iba(ik)
             wfsd_l(ii,ib,ik,1) = wfsd_l(ii,ib,ik,1)*norm
          end do
          bsdr_l(ib,1:nlmta,ik) = bsdr_l(ib,1:nlmta,ik)*norm
          bsdi_l(ib,1:nlmta,ik) = bsdi_l(ib,1:nlmta,ik)*norm
       else
          do ii=1,iba(ik)
             wfsdr = wfsd_l(ii,ib,ik,1   )
             wfsdi = wfsd_l(ii,ib,ik,kimg)
             norm = norm + wfsdr*wfsdr+wfsdi*wfsdi
          end do
          norm = 1.d0/sqrt(norm)
          do ii=1,iba(ik)
             wfsd_l(ii,ib,ik,1) = wfsd_l(ii,ib,ik,1)*norm
             wfsd_l(ii,ib,ik,2) = wfsd_l(ii,ib,ik,2)*norm
          enddo
          if(k_symmetry(ik) == GAMMA) then
            bsdr_l(ib,1:nlmta,ik) = bsdr_l(ib,1:nlmta,ik)*norm
          else
            bsdr_l(ib,1:nlmta,ik) = bsdr_l(ib,1:nlmta,ik)*norm
            bsdi_l(ib,1:nlmta,ik) = bsdi_l(ib,1:nlmta,ik)*norm
          end if
       end if
    end do

  end subroutine normalize_wfsd

! ===================== added by K. Tagami  ================== 11.0
  subroutine normalize_wfsd_noncl(ik)
    integer,intent(in) :: ik
    real(kind=DP) :: norm, wfsdr, wfsdi
    integer :: ib1,ii,ib, is, k1

    do ib1 = ista_e, iend_e, istep_e     ! MPI
       ib=map_z(ib1)
       norm = 0.d0
       if (kimg==1) then
          Do is=1, ndim_spinor
             k1 = ik + is -1
             do ii=1,iba(ik)
                wfsdr = wfsd_l(ii,ib,k1,kimg)
                norm = norm + wfsdr*wfsdr
             end do
          End do
          norm = 1.d0 /sqrt(norm)
          Do is=1, ndim_spinor
             k1 = ik + is -1
             do ii=1,iba(ik)
                wfsd_l(ii,ib,k1,1) = wfsd_l(ii,ib,k1,1)*norm
             end do
             bsdr_l(ib,1:nlmta,k1) = bsdr_l(ib,1:nlmta,k1)*norm
             bsdi_l(ib,1:nlmta,k1) = bsdi_l(ib,1:nlmta,k1)*norm
          End do
       else
          Do is=1, ndim_spinor
             k1 = ik + is -1
             do ii=1,iba(ik)
                wfsdr = wfsd_l(ii,ib,k1,1   )
                wfsdi = wfsd_l(ii,ib,k1,kimg)
                norm = norm + wfsdr*wfsdr+wfsdi*wfsdi
             end do
          End do
          norm = 1.d0/sqrt(norm)
          Do is=1, ndim_spinor
             k1 = ik + is -1
             do ii=1,iba(ik)
                wfsd_l(ii,ib,k1,1) = wfsd_l(ii,ib,k1,1)*norm
                wfsd_l(ii,ib,k1,2) = wfsd_l(ii,ib,k1,2)*norm
             enddo
             if(k_symmetry(ik) == GAMMA) then
                bsdr_l(ib,1:nlmta,k1) = bsdr_l(ib,1:nlmta,k1)*norm
             else
                bsdr_l(ib,1:nlmta,k1) = bsdr_l(ib,1:nlmta,k1)*norm
                bsdi_l(ib,1:nlmta,k1) = bsdi_l(ib,1:nlmta,k1)*norm
             end if
          End do
       end if
    end do

  end subroutine normalize_wfsd_noncl
! =========================================================== 11.0


  subroutine prepare_Hloc_phi(ik,ekin,afft,bfft,itot)
    integer, intent(in) :: ik
    integer, intent(in) :: itot
    real(kind=DP), intent(in)  :: ekin(kg1)
    real(kind=DP), intent(in)  :: afft(nfft)
    real(kind=DP), intent(out) :: bfft(nfft)

    integer       :: ib1,i1,ii,ib
    real(kind=DP) :: denom
    real(kind=DP) :: dr1,dr2,di1,di2,dd
    integer :: id_sname = -1, ipri0
    call tstatc0_begin('prepare_Hloc_phi (mddavidson) ', id_sname,1)

    call get_ipri0(iprimddavidson,ipri0)

    denom = 1.d0/product(fft_box_size_WF(1:3,1))

! (zaj_l <- (T+Vloc)|phi> )
!!    zaj_l(:,:,ik,:) = zajold_l(:,:,:,idavid)
!( tenchi ) (zat_l <- zaj_l)

    zat_l(:,:,:,1) = zaj_l(:,:,ik,:)
    zat_l(:,:,:,2) = wfsd_l(:,:,ik,:)
    fsr(:,:,1)=fsr_l(:,:,ik)
    fsr(:,:,2)=bsdr_l(:,:,ik)
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       fsi(:,:,1)=fsi_l(:,:,ik)
       fsi(:,:,2)=bsdi_l(:,:,ik)
    end if

    if(itot == 1) then
      do ib1 = ista_e, iend_e, istep_e     ! MPI
         call m_ES_WF_in_Rspace(ik,ib1,bfft)!(swffft)
         call m_FFT_Vlocal_W(afft,bfft)    ! (afft, bfft)-> (bfft)
         call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
         ib = map_z(ib1)                 ! MPI
         if(kimg == 1) then
            do ii=1,iba(ik)
               i1  = igf(nbase(ii,ik))
               dr1 = zaj_l(ii,ib,ik,1)
               dr2 = bfft(i1)*denom
  !             zaj_l(ii,ib,ik,1)= ekin(ii)*dr1+dr2
               zah_l(ii,ib,1,1) = ekin(ii)*dr1+dr2
            enddo
         else
            do ii=1,iba(ik)
               i1  = igf(nbase(ii,ik))
               dr1  = zaj_l(ii,ib,ik,1)
               di1  = zaj_l(ii,ib,ik,kimg)
               zah_l(ii,ib,1,1)= ekin(ii)*dr1+bfft(2*i1-1)*denom
               zah_l(ii,ib,kimg,1)= ekin(ii)*di1+bfft(2*i1)*denom
            enddo
         endif
      enddo
    end if

!!( tenchi ) (zah_l <- zaj_l)
    zaj_l(:,:,ik,:) = wfsd_l(:,:,ik,:)

    do ib1 = ista_e, iend_e, istep_e     ! MPI
       call m_ES_WF_in_Rspace(ik,ib1,bfft)!(swffft)
       call m_FFT_Vlocal_W(afft,bfft)    ! (afft, bfft)-> (bfft)
       call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
       ib = map_z(ib1)                 ! MPI
       if(kimg == 1) then
          do ii=1,iba(ik)
             i1  = igf(nbase(ii,ik))
             dr1 = zaj_l(ii,ib,ik,1)
             dr2 = bfft(i1)*denom
!             zaj_l(ii,ib,ik,1)= ekin(ii)*dr1+dr2
             zah_l(ii,ib,1,2) = ekin(ii)*dr1+dr2
          enddo
       else
          do ii=1,iba(ik)
             i1  = igf(nbase(ii,ik))
             dr1  = zaj_l(ii,ib,ik,1)
             di1  = zaj_l(ii,ib,ik,kimg)
!             zaj_l(ii,ib,ik,1)= ekin(ii)*dr1+bfft(2*i1-1)*denom
!             zaj_l(ii,ib,ik,kimg)= ekin(ii)*di1+bfft(2*i1)*denom
             zah_l(ii,ib,1,2)= ekin(ii)*dr1+bfft(2*i1-1)*denom
             zah_l(ii,ib,kimg,2)= ekin(ii)*di1+bfft(2*i1)*denom
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

  end subroutine prepare_Hloc_phi

! =========================== added by K. Tagami =============== 11.0
  subroutine prepare_Hloc_phi_noncl( ik, ekin, afft_kt, bfft_kt, itot )
    integer, intent(in) :: ik
    integer, intent(in) :: itot
    real(kind=DP), intent(in)  :: ekin(kg1)
    real(kind=DP), intent(inout)  :: afft_kt(nfft,ndim_chgpot)
    real(kind=DP), intent(out) :: bfft_kt(nfft,ndim_spinor)

    integer       :: ib1,i1,ii,ib, is, k1
    real(kind=DP) :: denom
    real(kind=DP) :: dr1,dr2,di1,di2,dd
    integer :: id_sname = -1, ipri0

    call tstatc0_begin('prepare_Hloc_phi_noncl (mddavidson) ', id_sname,1)

    call get_ipri0(iprimddavidson,ipri0)

    denom = 1.d0/product(fft_box_size_WF(1:3,1))


! (zaj_l <- (T+Vloc)|phi> )
!!    zaj_l(:,:,ik,:) = zajold_l(:,:,:,idavid)
!( tenchi ) (zat_l <- zaj_l)

    Do is=1, ndim_spinor
       zat_l_noncl(:,:,:,1,is) =  zaj_l(:,:,ik+is-1,:)
       zat_l_noncl(:,:,:,2,is) = wfsd_l(:,:,ik+is-1,:)
    End do

    Do is=1, ndim_spinor
       fsr_noncl(:,:,1,is) = fsr_l(:,:,ik+is-1)
       fsr_noncl(:,:,2,is) = bsdr_l(:,:,ik+is-1)
    End do
    if(.not.(kv3/ndim_spinor == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       Do is=1, ndim_spinor
          fsi_noncl(:,:,1,is) =  fsi_l(:,:,ik+is-1)
          fsi_noncl(:,:,2,is) = bsdi_l(:,:,ik+is-1)
       End do
    end if

    if(itot == 1) then
       do ib1 = ista_e, iend_e, istep_e     ! MPI
          Do is=1, ndim_spinor
             k1 = ik +is- 1
             call m_ES_WF_in_Rspace( k1, ib1, bfft_kt(:,is) )   !(swffft)
          End do
          call m_FFT_Vlocal_W_noncl( afft_kt, bfft_kt, ndim_chgpot, ndim_spinor )
                                       ! (afft, bfft)-> (bfft)
          Do is=1, ndim_spinor
             call m_FFT_WF( ELECTRON, nfout, bfft_kt(:,is), DIRECT, ON )
          End do

          ib = map_z(ib1)                 ! MPI
          if(kimg == 1) then
             Do is=1, ndim_spinor
                k1 = ik + is -1
                do ii=1,iba(ik)
                   i1  = igf(nbase(ii,ik))
                   dr1 = zaj_l( ii, ib, k1,1 )
                   dr2 = bfft_kt( i1,is )*denom
                   zah_l_noncl(ii,ib,1,1,is) = ekin(ii)*dr1+dr2
                enddo
             End do
          else
             Do is=1, ndim_spinor
                k1 = ik + is -1
                do ii=1,iba(ik)
                   i1  = igf(nbase(ii,ik))
                   dr1  = zaj_l( ii, ib, k1, 1 )
                   di1  = zaj_l( ii, ib, k1, kimg )
                   zah_l_noncl(ii,ib,   1,1,is) = ekin(ii)*dr1 &
                        &                        +bfft_kt(2*i1-1,is)*denom
                   zah_l_noncl(ii,ib,kimg,1,is) = ekin(ii)*di1 &
                        &                        +bfft_kt(2*i1,is)*denom
                enddo
             End do
          endif
       enddo
    end if

!!( tenchi ) (zah_l <- zaj_l)
    Do is=1, ndim_spinor
       k1 = ik + is -1
       zaj_l(:,:,k1,:) = wfsd_l(:,:,k1,:)
    End do

    do ib1 = ista_e, iend_e, istep_e     ! MPI
       Do is=1, ndim_spinor
          k1 = ik +is- 1
          call m_ES_WF_in_Rspace( k1, ib1, bfft_kt(:,is) )!(swffft)
       End do
       call m_FFT_Vlocal_W_noncl( afft_kt, bfft_kt, ndim_chgpot, ndim_spinor )
                                       ! (afft, bfft)-> (bfft)
       Do is=1, ndim_spinor
          call m_FFT_WF( ELECTRON, nfout, bfft_kt(:,is), DIRECT, ON )
       End do

       ib = map_z(ib1)                 ! MPI
       if (kimg == 1) then
          Do is=1, ndim_spinor
             k1 = ik + is -1
             do ii=1,iba(ik)
                i1  = igf(nbase(ii,ik))
                dr1 = zaj_l(ii,ib,k1,1)
                dr2 = bfft_kt(i1,is)*denom
                zah_l_noncl(ii,ib,1,2,is) = ekin(ii)*dr1+dr2
             enddo
          End Do
       else
          Do is=1, ndim_spinor
             k1 = ik + is -1
             do ii=1,iba(ik)
                i1  = igf(nbase(ii,ik))
                dr1  = zaj_l(ii,ib,k1,1)
                di1  = zaj_l(ii,ib,k1,kimg)
                zah_l_noncl(ii,ib,1,   2,is) = ekin(ii)*dr1 &
                     &                        +bfft_kt(2*i1-1,is)*denom
                zah_l_noncl(ii,ib,kimg,2,is) = ekin(ii)*di1 &
                     &                        +bfft_kt(2*i1,is)*denom
             enddo
          End do
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

  end subroutine prepare_Hloc_phi_noncl
! ================================================================= 11.0

  subroutine evolve_WFs_in_subspace &
                                (ik,ispin,ekin,afft,bfft,iblock,itot,frestart)
    integer, intent(in) :: ik,ispin
    integer, intent(in) :: iblock,itot
    logical, intent(out) :: frestart
    real(kind=DP), intent(in)  :: ekin(kg1)
    real(kind=DP), intent(in)  :: afft(nfft)
    real(kind=DP), intent(out) :: bfft(nfft)
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
             do ii = 1, iba(ik)            ! MPI
                hr2 = zah_l(ii,ii2,1,iter2)
                dr2 = zat_l(ii,ii2,1,iter2)
                dr1 = zat_l(ii,ii1,1,iter1)
                w1hw2(ip0) = w1hw2(ip0) + dr1*hr2
                w1sw2(ip0) = w1sw2(ip0) + dr1*dr2
             end do
          else
             if(k_symmetry(ik) == GAMMA) then
                w1hw2(ip0) = 0.d0
                w1sw2(ip0) = 0.d0
!!$                ig1 = 1;  if(mype == 0) ig1 = 2
!!                ig1 = 1;  if(myrank_e == 0) ig1 = 2
                ig1 = 2
                do ii = ig1, iba(ik)            ! MPI
                   hr2 = zah_l(ii,ii2,1,iter2) ! MPI
                   hi2 = zah_l(ii,ii2,2,iter2) ! MPI
                   dr2 = zat_l(ii,ii2,1,iter2) ! MPI
                   di2 = zat_l(ii,ii2,2,iter2) ! MPI
                   dr1 = zat_l(ii,ii1,1,iter1) ! MPI
                   di1 = zat_l(ii,ii1,2,iter1) ! MPI
                   w1hw2(ip0) =w1hw2(ip0)+(dr1*hr2+di1*hi2)*2.d0
                   w1sw2(ip0) =w1sw2(ip0)+(dr1*dr2+di1*di2)*2.d0
                end do
!!$                if(mype == 0) then
!!                if(myrank_e == 0) then
                hr2 = zah_l(1,ii2,1,iter2) ! MPI
                hi2 = zah_l(1,ii2,2,iter2) ! MPI
                dr2 = zat_l(1,ii2,1,iter2) ! MPI
                di2 = zat_l(1,ii2,2,iter2) ! MPI
                dr1 = zat_l(1,ii1,1,iter1) ! MPI
                di1 = zat_l(1,ii1,2,iter1) ! MPI
                w1hw2(ip0) =w1hw2(ip0)+dr1*hr2+di1*hi2
                w1sw2(ip0) =w1sw2(ip0)+dr1*dr2+di1*di2
!!                end if
             else
                w1hw2(2*ip0-1:2*ip0) = 0.d0
                w1sw2(2*ip0-1:2*ip0) = 0.d0
                do ii = 1, iba(ik)           ! MPI
                   hr2 = zah_l(ii,ii2,1,iter2) ! MPI
                   hi2 = zah_l(ii,ii2,2,iter2) ! MPI
                   dr2 = zat_l(ii,ii2,1,iter2) ! MPI
                   di2 = zat_l(ii,ii2,2,iter2) ! MPI
                   dr1 = zat_l(ii,ii1,1,iter1) ! MPI
                   di1 = zat_l(ii,ii1,2,iter1) ! MPI
                   w1hw2(2*ip0-1) =w1hw2(2*ip0-1)+dr1*hr2+di1*hi2
                   w1hw2(2*ip0  ) =w1hw2(2*ip0  )+dr1*hi2-di1*hr2
                   w1sw2(2*ip0-1) =w1sw2(2*ip0-1)+dr1*dr2+di1*di2
                   w1sw2(2*ip0  ) =w1sw2(2*ip0  )+dr1*di2-di1*dr2
                end do
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
      integer :: ib1,ib2,ibb2
!!$      real(kind=DP), dimension(np_g1k_x,neg,kimg) :: zaj_wk
      real(kind=DP), allocatable, dimension(:,:,:) :: zaj_wk
      real(kind=DP), allocatable, dimension(:,:,:) :: zah_wk
      real(kind=DP), allocatable, dimension(:,:)   :: fsr_wk
      real(kind=DP), allocatable, dimension(:,:)   :: fsi_wk
      allocate(zaj_wk(kg1,nel,kimg))
      if(itot /= max_iter_mddavid) then
        allocate(zah_wk(kg1,nel,kimg))
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
               do ii=1,iba(ik)
!                  zaj_wk(ii,ib1,kimg) = zaj_wk(ii,ib1,kimg) + zat_l(ii,ii2,kimg,iter2)*vec(ib2,ib1)
                  zaj_wk(ii,ib1,kimg) = zaj_wk(ii,ib1,kimg) + zat_l(ii,ii2,kimg,iter2)*hr2
               end do
               do ii=1,nlmta
                 fsr_wk(ib1,ii) = fsr_wk(ib1,ii) + fsr(ii2,ii,iter2)*hr2
                 fsi_wk(ib1,ii) = fsi_wk(ib1,ii) + fsi(ii2,ii,iter2)*hr2
               end do
               if(itot /= max_iter_mddavid) then
                 do ii=1,iba(ik)
                    zah_wk(ii,ib1,kimg) = zah_wk(ii,ib1,kimg) + zah_l(ii,ii2,kimg,iter2)*hr2
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
                  do ii=1,iba(ik)
                     dr1=zat_l(ii,ii2,1   ,iter2)
                     di1=zat_l(ii,ii2,kimg,iter2)
                     zaj_wk(ii,ib1,1   ) = zaj_wk(ii,ib1,1   ) + dr1*hr2
                     zaj_wk(ii,ib1,kimg) = zaj_wk(ii,ib1,kimg) + di1*hr2
                  end do
                  do ii=1,nlmta
                    fsr_wk(ib1,ii) = fsr_wk(ib1,ii) + fsr(ii2,ii,iter2)*hr2
                  end do
                  if(itot /= max_iter_mddavid) then
                    do ii=1,iba(ik)
                       dr1=zah_l(ii,ii2,1   ,iter2)
                       di1=zah_l(ii,ii2,kimg,iter2)
                       zah_wk(ii,ib1,1   ) = zah_wk(ii,ib1,1   ) + dr1*hr2
                       zah_wk(ii,ib1,kimg) = zah_wk(ii,ib1,kimg) + di1*hr2
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
                  do ii=1,iba(ik)
                     dr1=zat_l(ii,ii2,1   ,iter2)
                     di1=zat_l(ii,ii2,kimg,iter2)
                     zaj_wk(ii,ib1,1   ) = zaj_wk(ii,ib1,1   ) + dr1*hr2 - di1*hi2
                     zaj_wk(ii,ib1,kimg) = zaj_wk(ii,ib1,kimg) + dr1*hi2 + di1*hr2
                  end do
                  do ii=1,nlmta
                     dr1=fsr(ii2,ii,iter2)
                     di1=fsi(ii2,ii,iter2)
                     fsr_wk(ib1,ii) = fsr_wk(ib1,ii) + dr1*hr2 - di1*hi2
                     fsi_wk(ib1,ii) = fsi_wk(ib1,ii) + dr1*hi2 + di1*hr2
                  end do
                  if(itot /= max_iter_mddavid) then
                    do ii=1,iba(ik)
                       dr1=zah_l(ii,ii2,1   ,iter2)
                       di1=zah_l(ii,ii2,kimg,iter2)
                       zah_wk(ii,ib1,1   ) = zah_wk(ii,ib1,1   ) + dr1*hr2 - di1*hi2
                       zah_wk(ii,ib1,kimg) = zah_wk(ii,ib1,kimg) + dr1*hi2 + di1*hr2
                    end do
                  end if
               end do
            end do
         end if
      end if
!print *,itot,  itot /= max_iter_mddavid
      zaj_l(:,nsta:nend,ik,:) = zaj_wk(:,:,:)
      if(itot /= max_iter_mddavid) zah_l(:,nsta:nend,:,1) = zah_wk(:,:,:)
      fsr_l(nsta:nend,:,ik) = fsr_wk(:,:)
      if(k_symmetry(ik) /= GAMMA) fsi_l(nsta:nend,:,ik) = fsi_wk(:,:)

      deallocate(zaj_wk,fsr_wk)
      if(itot /= max_iter_mddavid) deallocate(zah_wk)
      if(k_symmetry(ik) /= GAMMA) deallocate(fsi_wk)
    end subroutine subspace_rotation

  end subroutine evolve_WFs_in_subspace


! ======================================= added by K. Tagami ============= 11.0
  subroutine evolve_WFs_in_subspace_noncl( ik, ekin, afft_kt, bfft_kt, &
       &                                   iblock,itot,frestart )
    integer, intent(in) :: ik
    integer, intent(in) :: iblock,itot
    logical, intent(out) :: frestart
    real(kind=DP), intent(in)  :: ekin(kg1)
    real(kind=DP), intent(inout)  :: afft_kt(nfft,ndim_chgpot)
    real(kind=DP), intent(out) :: bfft_kt(nfft,ndim_spinor)

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

    integer :: is

    call tstatc0_begin('evolve_WFs_in_subspace_noncl (modified davidson) ', id_sname,1)

    call get_ipri0(iprimddavidson,ipri0)

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

          if (kimg == 1) then
             w1hw2(ip0) = 0.d0
             w1sw2(ip0) = 0.d0
             Do is=1, ndim_spinor
                do ii = 1, iba(ik)            ! MPI
                   hr2 = zah_l_noncl(ii,ii2,1,iter2,is)
                   dr2 = zat_l_noncl(ii,ii2,1,iter2,is)
                   dr1 = zat_l_noncl(ii,ii1,1,iter1,is)
                   w1hw2(ip0) = w1hw2(ip0) + dr1*hr2
                   w1sw2(ip0) = w1sw2(ip0) + dr1*dr2
                end do
             End do
          else
             if(k_symmetry(ik) == GAMMA) then
                w1hw2(ip0) = 0.d0
                w1sw2(ip0) = 0.d0
!!$                ig1 = 1;  if(mype == 0) ig1 = 2
!!                ig1 = 1;  if(myrank_e == 0) ig1 = 2
                ig1 = 2
                Do is=1, ndim_spinor
                   do ii = ig1, iba(ik)            ! MPI
                      hr2 = zah_l_noncl(ii,ii2,1,iter2,is) ! MPI
                      hi2 = zah_l_noncl(ii,ii2,2,iter2,is) ! MPI
                      dr2 = zat_l_noncl(ii,ii2,1,iter2,is) ! MPI
                      di2 = zat_l_noncl(ii,ii2,2,iter2,is) ! MPI
                      dr1 = zat_l_noncl(ii,ii1,1,iter1,is) ! MPI
                      di1 = zat_l_noncl(ii,ii1,2,iter1,is) ! MPI
                      w1hw2(ip0) =w1hw2(ip0)+(dr1*hr2+di1*hi2)*2.d0
                      w1sw2(ip0) =w1sw2(ip0)+(dr1*dr2+di1*di2)*2.d0
                   end do
                End do
!!$                if(mype == 0) then
!!                if(myrank_e == 0) then

                Do is=1, ndim_spinor
                   hr2 = zah_l_noncl(1,ii2,1,iter2,is) ! MPI
                   hi2 = zah_l_noncl(1,ii2,2,iter2,is) ! MPI
                   dr2 = zat_l_noncl(1,ii2,1,iter2,is) ! MPI
                   di2 = zat_l_noncl(1,ii2,2,iter2,is) ! MPI
                   dr1 = zat_l_noncl(1,ii1,1,iter1,is) ! MPI
                   di1 = zat_l_noncl(1,ii1,2,iter1,is) ! MPI
                   w1hw2(ip0) =w1hw2(ip0)+dr1*hr2+di1*hi2
                   w1sw2(ip0) =w1sw2(ip0)+dr1*dr2+di1*di2
                End do
!!                end if

             else
                w1hw2(2*ip0-1:2*ip0) = 0.d0
                w1sw2(2*ip0-1:2*ip0) = 0.d0
                Do is=1, ndim_spinor
                   do ii = 1, iba(ik)           ! MPI
                      hr2 = zah_l_noncl(ii,ii2,1,iter2,is) ! MPI
                      hi2 = zah_l_noncl(ii,ii2,2,iter2,is) ! MPI
                      dr2 = zat_l_noncl(ii,ii2,1,iter2,is) ! MPI
                      di2 = zat_l_noncl(ii,ii2,2,iter2,is) ! MPI
                      dr1 = zat_l_noncl(ii,ii1,1,iter1,is) ! MPI
                      di1 = zat_l_noncl(ii,ii1,2,iter1,is) ! MPI
                      w1hw2(2*ip0-1) =w1hw2(2*ip0-1)+dr1*hr2+di1*hi2
                      w1hw2(2*ip0  ) =w1hw2(2*ip0  )+dr1*hi2-di1*hr2
                      w1sw2(2*ip0-1) =w1sw2(2*ip0-1)+dr1*dr2+di1*di2
                      w1sw2(2*ip0  ) =w1sw2(2*ip0  )+dr1*di2-di1*dr2
                   end do
                End do
             end if
          end if
       end do
    end do
    if(iprimddavidson >= 3) call wd_w1hw2(" -- w1hw2 without nl part--",iblock)
    ! <n|Vnl|m>
    call add_nonlocal_part_noncl ! w1hw2 = w1hw2 + w1Vnlw2
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

       Do is=1, ndim_spinor
          zaj_l(:,nsta:nend,ik+is-1,:) = zat_l_noncl(:,nsta:nend,:,1,is)
       End do

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
       call subspace_rotation_noncl ! vec,zat_l -> zat_l

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
1201   format(' %% for ik = ',i4,4x,' eko1&ekod&eko2 = ',3f14.7)
9001   format(5x,6f12.5)
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

    subroutine add_nonlocal_part_noncl
      integer :: ip,ib1,ib2,ibb1,ibb2
      integer       :: ia, lmt1, lmt2, it, p, s, ib
      real(kind=DP) :: tmpr,tmpi

      integer :: is1, is2, is_tmp

      complex(kind=CMPLDP) :: facv,facq
      real(kind=DP) :: vr,vi,qr,qi, cv1, cv2, cq1, cq2

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

                     Do is1=1, ndim_spinor
                        Do is2=1, ndim_spinor
                           is_tmp = 2*( is1 -1 ) +is2
                           facv = iwei(ia) *dion_scr_noncl( lmt1, lmt2, is_tmp, ia )
                           facq = iwei(ia) *q_noncl( lmt1, lmt2, is_tmp, it )

                           if ( kimg==1 ) then
                              tmpr =  fsr_noncl( ii1, p, iter1, is1 ) &
                                   & *fsr_noncl( ii2, s, iter2, is2 ) &
                                   & +fsi_noncl( ii1, p, iter1, is1 ) &
                                   & *fsi_noncl( ii2, s, iter2, is2 )
                              tmpi = fsr_noncl( ii1, p, iter1, is1 )&
        &                           *fsi_noncl( ii2, s, iter2, is2 )&
        &                          - fsi_noncl( ii1, p, iter1, is1 )&
        &                           *fsr_noncl( ii2, s, iter2, is2 )

                              vr = vr + real( facv )*tmpr -aimag( facv )*tmpi
                              qr = qr + real( facq )*tmpr -aimag( facq )*tmpi

                           else
                              if ( k_symmetry(ik) == GAMMA ) then
                                 tmpr = fsr_noncl( ii1, p, iter1, is1 )&
        &                              *fsr_noncl( ii2, s, iter2, is2 )
                                 vr = vr + real( facv )*tmpr
                                 qr = qr + real( facq )*tmpr
                              else
                                 tmpr = fsr_noncl( ii1, p, iter1, is1 )&
        &                              *fsr_noncl( ii2, s, iter2, is2 )&
        &                             + fsi_noncl( ii1, p, iter1, is1 )&
        &                              *fsi_noncl( ii2, s, iter2, is2 )
                                 tmpi = fsr_noncl( ii1, p, iter1, is1 )&
        &                              *fsi_noncl( ii2, s, iter2, is2 )&
        &                             - fsi_noncl( ii1, p, iter1, is1 )&
        &                              *fsr_noncl( ii2, s, iter2, is2 )

                                 cv1 = real(facv);  cv2 = aimag(facv)
                                 cq1 = real(facq);  cq2 = aimag(facq)

                                 vr = vr + cv1 *tmpr - cv2 *tmpi
                                 vi = vi + cv1 *tmpi + cv2 *tmpr
                                 qr = qr + cq1 *tmpr - cq2 *tmpi
                                 qi = qi + cq1 *tmpi + cq2 *tmpr
                              endif
                           endif
                        End do
                     End do

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
    end subroutine add_nonlocal_part_noncl

    subroutine subspace_rotation_noncl
      integer :: ib1,ib2,ibb2
!!$      real(kind=DP), dimension(np_g1k_x,neg,kimg) :: zaj_wk
      real(kind=DP), allocatable, dimension(:,:,:,:) :: zaj_wk
      real(kind=DP), allocatable, dimension(:,:,:,:) :: zah_wk
      real(kind=DP), allocatable, dimension(:,:,:)   :: fsr_wk
      real(kind=DP), allocatable, dimension(:,:,:)   :: fsi_wk

      integer :: is, k1

      allocate(zaj_wk(kg1,nel,kimg,ndim_spinor))
      if(itot /= max_iter_mddavid) then
        allocate(zah_wk(kg1,nel,kimg,ndim_spinor))
        zah_wk(:,:,:,:) = 0.d0
      end if

      allocate(fsr_wk(nel,nlmta,ndim_spinor))
      if(k_symmetry(ik) /= GAMMA) then
        allocate(fsi_wk(nel,nlmta,ndim_spinor))
        fsi_wk(:,:,:)=0.d0
      end if

      zaj_wk(:,:,:,:) = 0.d0
      fsr_wk(:,:,:)=0.d0

      if (kimg==1) then
         do ib1=1,nel
            do ibb2=1,nsize_max_sb_now
               if(ibover(ibb2,iblock)<0) cycle
               ib2 = ibover(ibb2,iblock)
               iter2=(ibb2-1)/nel+1
               ii2=ibb2-nel*(iter2-1)
               ii2=nsta+ii2-1
               hr2=vec(ib2,ib1)

               Do is=1, ndim_spinor
                  do ii=1,iba(ik)
                     zaj_wk(ii,ib1,kimg,is) = zaj_wk(ii,ib1,kimg,is) &
                          &                 + zat_l_noncl(ii,ii2,kimg,iter2,is)*hr2
                  end do
               End do

               Do is=1, ndim_spinor
                  do ii=1,nlmta
                     fsr_wk(ib1,ii,is) = fsr_wk(ib1,ii,is) &
                          &            + fsr_noncl(ii2,ii,iter2,is) *hr2
                     fsi_wk(ib1,ii,is) = fsi_wk(ib1,ii,is) &
                          &            + fsi_noncl(ii2,ii,iter2,is)*hr2
                  end do
               End do
               if (itot /= max_iter_mddavid) then
                  Do is=1, ndim_spinor
                     do ii=1,iba(ik)
                        zah_wk(ii,ib1,kimg,is) = zah_wk(ii,ib1,kimg,is)  &
                             &                 + zah_l_noncl(ii,ii2,kimg,iter2,is)*hr2
                     end do
                  End do
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

                  Do is=1, ndim_spinor
                     do ii=1,iba(ik)
                        dr1 = zat_l_noncl(ii,ii2,1   ,iter2,is)
                        di1 = zat_l_noncl(ii,ii2,kimg,iter2,is)
                        zaj_wk(ii,ib1,1,   is) = zaj_wk(ii,ib1,1,   is) + dr1*hr2
                        zaj_wk(ii,ib1,kimg,is) = zaj_wk(ii,ib1,kimg,is) + di1*hr2
                     end do
                  End do
                  Do is=1, ndim_spinor
                     do ii=1,nlmta
                        fsr_wk(ib1,ii,is) = fsr_wk(ib1,ii,is) &
                             &            + fsr_noncl(ii2,ii,iter2,is)*hr2
                     end do
                  End Do
                  if(itot /= max_iter_mddavid) then
                     Do is=1, ndim_spinor
                        do ii=1,iba(ik)
                           dr1 = zah_l_noncl(ii,ii2,1   ,iter2,is)
                           di1 = zah_l_noncl(ii,ii2,kimg,iter2,is)
                           zah_wk(ii,ib1,1,   is) = zah_wk(ii,ib1,1,   is) +dr1*hr2
                           zah_wk(ii,ib1,kimg,is) = zah_wk(ii,ib1,kimg,is) +di1*hr2
                        end do
                     End do
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

                  Do is=1, ndim_spinor
                     do ii=1,iba(ik)
                        dr1 = zat_l_noncl(ii,ii2,1   ,iter2,is)
                        di1 = zat_l_noncl(ii,ii2,kimg,iter2,is)
                        zaj_wk(ii,ib1,1   ,is) = zaj_wk(ii,ib1,1,   is) +dr1*hr2 -di1*hi2
                        zaj_wk(ii,ib1,kimg,is) = zaj_wk(ii,ib1,kimg,is) +dr1*hi2 +di1*hr2
                     end do
                  End Do
                  Do is=1, ndim_spinor
                     do ii=1,nlmta
                        dr1 = fsr_noncl(ii2,ii,iter2,is)
                        di1 = fsi_noncl(ii2,ii,iter2,is)
                        fsr_wk(ib1,ii,is) = fsr_wk(ib1,ii,is) + dr1*hr2 - di1*hi2
                        fsi_wk(ib1,ii,is) = fsi_wk(ib1,ii,is) + dr1*hi2 + di1*hr2
                     end do
                  End do
                  if (itot /= max_iter_mddavid) then
                     Do is=1, ndim_spinor
                        do ii=1,iba(ik)
                           dr1 = zah_l_noncl(ii,ii2,1   ,iter2,is)
                           di1 = zah_l_noncl(ii,ii2,kimg,iter2,is)
                           zah_wk(ii,ib1,1,   is) = zah_wk(ii,ib1,1,   is) &
                                &                  +dr1*hr2 -di1*hi2
                           zah_wk(ii,ib1,kimg,is) = zah_wk(ii,ib1,kimg,is) &
                                &                  +dr1*hi2 +di1*hr2
                        end do
                     End do
                  end if
               end do
            end do
         end if
      end if
!print *,itot,  itot /= max_iter_mddavid

      Do is=1, ndim_spinor
         k1 = ik + is -1
         zaj_l(:,nsta:nend,k1,:) = zaj_wk(:,:,:,is)
      End do

      if(itot /= max_iter_mddavid) then
         Do is=1, ndim_spinor
            zah_l_noncl(:,nsta:nend,:,1,is) = zah_wk(:,:,:,is)
         End do
      endif
      Do is=1, ndim_spinor
         k1 = ik + is -1
         fsr_l(nsta:nend,:,k1) = fsr_wk(:,:,is)
         if (k_symmetry(ik) /= GAMMA) then
            fsi_l(nsta:nend,:,k1) = fsi_wk(:,:,is)
         endif
      End do

      deallocate(zaj_wk,fsr_wk)
      if(itot /= max_iter_mddavid) deallocate(zah_wk)
      if(k_symmetry(ik) /= GAMMA) deallocate(fsi_wk)

    end subroutine subspace_rotation_noncl

  end subroutine evolve_WFs_in_subspace_noncl
! ============================================================== 11.0



  subroutine m_ESmddavid_Subspace_Rotation(nfout)
    integer, intent(in) :: nfout

    integer             :: ispin, ik, iksnl, switch_of_eko_part
    real(kind=DP), allocatable, dimension(:) ::  afft, bfft
    real(kind=DP), allocatable, dimension(:) :: ekin
    integer :: ipri0

    allocate(ekin(kg1))
    call m_ES_alloc_scss_etc()
    allocate(afft(nfft)); allocate(bfft(nfft))
    call allocate_fsri

    do ispin = 1, nspin, (af+1)
       call m_ES_Vlocal_in_Rspace(ispin,afft)      ! (ptfft1) vlhxc_l->afft
       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) /= myrank_k) cycle          ! MPI
          iksnl = (ik-1)/nspin + 1

          call allocate_t_matrix_sr(ik) ! -> np_g1k_x
          call m_pwBS_kinetic_energies(ik,vkxyz,ekin) ! (diakin) ->ekin
          call m_ES_Vnonlocal_W(ik,iksnl,ispin,switch_of_eko_part=OFF) ! -> vnlph_l
          call allreduce_fs_sr(ik) ! -> fsr,fsi
          call evolve_WFs_in_subspace_sr&     !-(m_ES_WF_by_ModifiedDavidson)
                                      &(ik,ispin,ekin,afft,bfft) !-> zaj_l
          if(ik==1.and.iprimddavidson>= 2) &
            & call m_ES_wd_zaj_small_portion(nfout,ik," -- after davidson subspace rotation --",21)
          call m_ES_betar_dot_WFs_4_each_k(nfout,ik)   ! -> fsr_l,fsi_l
          if(iprimddavidson>=2) then
             write(nfout,'("Davidson Subspace Rotation: ik=",i5," subspace=",i5)') ik, nsize_sb_now
          end if
          call deallocate_t_matrix_sr

       enddo      ! k-point loop
    enddo      ! spin loop

    call deallocate_fsri

!!  ( in case of af=1 )
    if(af /= 0) then
       call cp_eigen_values_for_af       !-(contained here)
       call expand_neordr_and_nrvf_ordr  !-(contained here)
    end if

    call get_ipri0(iprimddavidson,ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko(nfout,mode=SCF)

    deallocate(bfft);   deallocate(afft)
    call m_ES_dealloc_scss_etc()
    deallocate(ekin)

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

! ================================ added by K. Tagami ================ 11.0
  subroutine m_ESmddavid_Subspace_Rot_noncl(nfout)
    integer, intent(in) :: nfout

    integer             :: ik, iksnl, switch_of_eko_part
    real(kind=DP), allocatable ::  afft_kt(:,:)
    real(kind=DP), allocatable ::  bfft_kt(:,:)
    real(kind=DP), allocatable :: vnlph_noncl(:,:,:,:)
    real(kind=DP), allocatable, dimension(:) :: ekin
    integer :: ipri0

    integer :: precon, is1, is2, istmp, k2

    precon = ON

    allocate(ekin(kg1))
    call m_ES_alloc_scss_etc()
    allocate(afft_kt(nfft,ndim_chgpot)); afft_kt = 0.0d0
    allocate(bfft_kt(nfft,ndim_spinor)); bfft_kt = 0.0d0
    allocate(vnlph_noncl(kg1,np_e,kimg,ndim_spinor)); vnlph_noncl = 0.0d0

    call allocate_fsri_noncl

    call m_ES_Vlocal_in_Rspace_noncl( afft_kt )    ! (ptfft1) vlhxc_l->afft

    Do ik=1, kv3, ndim_spinor
       if ( map_k(ik) /= myrank_k )  cycle          ! MPI
       iksnl = (ik-1)/ndim_spinor + 1

       call allocate_t_matrix_sr_noncl(ik)         ! -> np_g1k_x
       call m_pwBS_kinetic_energies( ik, vkxyz, ekin ) ! (diakin) ->ekin

       vnlph_noncl = 0.0d0
       Do is1=1, ndim_spinor
          Do is2=1, ndim_spinor
             istmp = ( is1-1 )*ndim_spinor +is2
             k2 = ik +is2 -1

             if ( precon==ON ) then
                call m_ES_Vnonlocal_W( k2, iksnl, istmp, switch_of_eko_part=OFF )
                ! -> vnlph_l
             else
                call m_ES_Vnonlocal_W( k2, iksnl, istmp, switch_of_eko_part=ON )
                ! -> vnlph_l
             endif
             vnlph_noncl(:,:,:,is1) = vnlph_noncl(:,:,:,is1) &
                  &                     + vnlph_l(:,:,:)
          End do
       End do

       call allreduce_fs_sr_noncl(ik) ! -> fsr,fsi
       call evolve_WFs_in_subspace_sr_noncl( ik, ekin, afft_kt, bfft_kt )
                                             !-> zaj_l
                                             !-(m_ES_WF_by_ModifiedDavidson)

       if (ik==1.and.iprimddavidson>= 2) then
          Do is1=1, ndim_spinor
             call m_ES_wd_zaj_small_portion( nfout, ik+is1-1, &
                  &                         " -- after davidson subspace rotation --",21 )
          End do
       endif

       Do is1=1, ndim_spinor
          call m_ES_betar_dot_WFs_4_each_k( nfout,ik+is1-1 )   ! -> fsr_l,fsi_l
       End do

       if(iprimddavidson>=2) then
          write(nfout,'("Davidson Subspace Rotation: ik=",i5," subspace=",i5)') &
               &        ik, nsize_sb_now
       end if
       call deallocate_t_matrix_sr

    enddo      ! k-point loop

    call deallocate_fsri_noncl

    call m_ES_sort_eigen_vals_noncl()

    call get_ipri0(iprimddavidson,ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko(nfout,mode=SCF)

    deallocate(bfft_kt);   deallocate(afft_kt)
    call m_ES_dealloc_scss_etc()
    deallocate(ekin)

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

  end subroutine m_ESmddavid_Subspace_Rot_noncl
! ================================================================= 11.0

  subroutine allocate_fsri
    nsize_sb_now = neg
    nsize_mt_now =  nsize_sb_now*(nsize_sb_now+1)/2
    allocate(fsr(neg,nlmta,1))
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       allocate(fsi(neg,nlmta,1))
    end if
    allocate(zaj_l_backup(kg1,np_e,kimg)) ! MPI
  end subroutine allocate_fsri

! ===================== added by K. Tagami ==================== 11.0
  subroutine allocate_fsri_noncl
    nsize_sb_now = neg
    nsize_mt_now =  nsize_sb_now*(nsize_sb_now+1)/2
    allocate(fsr_noncl(neg,nlmta,1,ndim_spinor))
    if(.not.(kv3/ndim_spinor == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       allocate(fsi_noncl(neg,nlmta,1,ndim_spinor))
    end if
    allocate(zaj_l_backup_noncl(kg1,np_e,kimg,ndim_spinor)) ! MPI
  end subroutine allocate_fsri_noncl
! ============================================================== 11.0

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

  subroutine allocate_t_matrix_sr(ik)
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
    allocate(zat_l(np_g1k_x,neg,kimg,1)) ! MPI
    allocate(zah_l(np_g1k_x,neg,kimg,1)) ! MPI
    allocate(w1hw2(nsize_mt_now*kimg_t))
    allocate(w1sw2(nsize_mt_now*kimg_t))
    if(npes>1) then
       allocate(w1hw2_mpi(nsize_mt_now*kimg_t))
       allocate(w1sw2_mpi(nsize_mt_now*kimg_t))
    end if
    zaj_l_backup(:,:,:) = zaj_l(:,:,ik,:)
  end subroutine allocate_t_matrix_sr

! =============================== added by K. Tagami ============= 11.0
  subroutine allocate_t_matrix_sr_noncl(ik)
    integer, intent(in) :: ik
    integer :: kimg_t, is

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

    allocate(zat_l_noncl(np_g1k_x,neg,kimg,1,ndim_spinor)) ! MPI
    allocate(zah_l_noncl(np_g1k_x,neg,kimg,1,ndim_spinor)) ! MPI
    allocate(w1hw2(nsize_mt_now*kimg_t))
    allocate(w1sw2(nsize_mt_now*kimg_t))
    if(npes>1) then
       allocate(w1hw2_mpi(nsize_mt_now*kimg_t))
       allocate(w1sw2_mpi(nsize_mt_now*kimg_t))
    end if
    Do is=1, ndim_spinor
       zaj_l_backup_noncl(:,:,:,is) = zaj_l(:,:,ik+is-1,:)
    End Do
  end subroutine allocate_t_matrix_sr_noncl
! ============================================================ 11.0

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

  subroutine allreduce_fs_sr(ik)
    integer, intent(in) :: ik
    integer :: ib,ib1,kimg_t
    real(kind=DP), allocatable, dimension(:,:) :: fs_mpi,fs_mpi2

    allocate(fs_mpi(neg,nlmta))
    allocate(fs_mpi2(neg,nlmta))

    fs_mpi=0.d0
    do ib1 = ista_e, iend_e, istep_e     ! MPI
       ib=map_z(ib1)
       fs_mpi(ib1,1:nlmta) = fsr_l(ib,1:nlmta,ik)
    end do
    fs_mpi2=0.d0
    call mpi_allreduce(fs_mpi,fs_mpi2,neg*nlmta &
      & ,mpi_double_precision,mpi_sum &
      & ,mpi_k_world(myrank_k),ierr)       ! MPI
    fsr(1:neg,1:nlmta,1) = fs_mpi2(1:neg,1:nlmta)
    if(.not. k_symmetry(ik) == GAMMA) then
       fs_mpi=0.d0
       do ib1 = ista_e, iend_e, istep_e     ! MPI
          ib=map_z(ib1)
          fs_mpi(ib1,1:nlmta) = fsi_l(ib,1:nlmta,ik)
       end do
       fs_mpi2=0.d0
       call mpi_allreduce(fs_mpi,fs_mpi2,neg*nlmta &
         & ,mpi_double_precision,mpi_sum &
         & ,mpi_k_world(myrank_k),ierr)       ! MPI
       fsi(1:neg,1:nlmta,1) = fs_mpi2(1:neg,1:nlmta)
    end if
    deallocate(fs_mpi)
    deallocate(fs_mpi2)
  end subroutine allreduce_fs_sr

! ================================= added by K. Tagami ============= 11.0
  subroutine allreduce_fs_sr_noncl(ik)
    integer, intent(in) :: ik

    integer :: ib, ib1, kimg_t, is
    real(kind=DP), allocatable, dimension(:,:,:) :: fs_mpi,fs_mpi2

    allocate(fs_mpi(neg,nlmta,ndim_spinor))
    allocate(fs_mpi2(neg,nlmta,ndim_spinor))

    fs_mpi=0.d0
    do ib1 = ista_e, iend_e, istep_e     ! MPI
       ib=map_z(ib1)
       Do is=1, ndim_spinor
          fs_mpi(ib1,1:nlmta,is) = fsr_l(ib,1:nlmta,ik+is-1)
       End Do
    end do

    fs_mpi2=0.d0
    call mpi_allreduce( fs_mpi, fs_mpi2, neg*nlmta*ndim_spinor, &
         &              mpi_double_precision, mpi_sum, &
         &              mpi_k_world(myrank_k),ierr )       ! MPI

    fsr_noncl(1:neg,1:nlmta,1,:) = fs_mpi2(1:neg,1:nlmta,:)

    if (.not. k_symmetry(ik) == GAMMA) then
       fs_mpi=0.d0
       do ib1 = ista_e, iend_e, istep_e     ! MPI
          ib=map_z(ib1)
          Do is=1, ndim_spinor
             fs_mpi(ib1,1:nlmta,is) = fsi_l(ib,1:nlmta,ik+is-1)
          End do
       end do

       fs_mpi2=0.d0
       call mpi_allreduce( fs_mpi, fs_mpi2, neg*nlmta*ndim_spinor, &
            &              mpi_double_precision, mpi_sum, &
            &              mpi_k_world(myrank_k),ierr )       ! MPI

       fsi_noncl(1:neg,1:nlmta,1,:) = fs_mpi2(1:neg,1:nlmta,:)
    end if

    deallocate(fs_mpi);    deallocate(fs_mpi2)

  end subroutine allreduce_fs_sr_noncl
! ============================================================= 11.0

  subroutine evolve_WFs_in_subspace_sr(ik,ispin,ekin,afft,bfft)
    integer, intent(in) :: ik,ispin
    real(kind=DP), intent(in)  :: ekin(kg1)
    real(kind=DP), intent(in)  :: afft(nfft)
    real(kind=DP), intent(out) :: bfft(nfft)
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
    call tstatc0_begin('evolve_WFs_in_subspace_sr(davidson) ', id_sname,1)

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


    do ib1 = 1, neg
       if(map_e(ib1) == myrank_e) eko_d(ib1) = eko_l(map_z(ib1),ik)  ! MPI
    end do
    call mpi_allreduce(eko_d,eko_d_mpi,neg,mpi_double_precision,mpi_sum &
         & ,mpi_k_world(myrank_k),ierr)       ! MPI
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
    call m_ES_W_transpose(ista_k,iend_k,ik,zaj_l,zat_l(1,1,1,1))
    do ib1 = ista_e, iend_e, istep_e     ! MPI
       call m_ES_WF_in_Rspace(ik,ib1,bfft)!(swffft)
       call m_FFT_Vlocal_W(afft,bfft)    ! (afft, bfft)-> (bfft)
       call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
       ib = map_z(ib1)                 ! MPI
       if(kimg == 1) then
          do ii=1,iba(ik)
             i1  = igf(nbase(ii,ik))
             dr1 = zaj_l(ii,ib,ik,1)
             dr2 = bfft(i1)*denom
             zaj_l(ii,ib,ik,1)= ekin(ii)*dr1+dr2
          enddo
       else
          do ii=1,iba(ik)
             i1  = igf(nbase(ii,ik))
             dr1  = zaj_l(ii,ib,ik,1)
             di1  = zaj_l(ii,ib,ik,kimg)
             zaj_l(ii,ib,ik,1)= ekin(ii)*dr1+bfft(2*i1-1)*denom
             zaj_l(ii,ib,ik,kimg)= ekin(ii)*di1+bfft(2*i1)*denom
          enddo
       endif
    enddo
!( tenchi ) (zah_l <- zaj_l)
    call m_ES_W_transpose(ista_k,iend_k,ik,zaj_l,zah_l(1,1,1,1))

! (make matrix elements )
! parallel loop
    ! <n|T+Vloc|m> G-wise parallel
    do ib2 = 1,nsize_sb_now
       ip0b = ib2*(ib2-1)/2
       do ib1 = 1,ib2
          ip0 = ip0b + ib1
          if(kimg == 1) then
             w1hw2(ip0) = 0.d0
             w1sw2(ip0) = 0.d0
             do ii = 1, np_g1k(ik)            ! MPI
                hr2 = zah_l(ii,ib2,1,1)
                dr2 = zat_l(ii,ib2,1,1)
                dr1 = zat_l(ii,ib1,1,1)
                w1hw2(ip0) = w1hw2(ip0) + dr1*hr2
                w1sw2(ip0) = w1sw2(ip0) + dr1*dr2
             end do
          else
             if(k_symmetry(ik) == GAMMA) then
                w1hw2(ip0) = 0.d0
                w1sw2(ip0) = 0.d0
!!$                ig1 = 1;  if(mype == 0) ig1 = 2
                ig1 = 1;  if(myrank_e == 0) ig1 = 2
                do ii = ig1, np_g1k(ik)            ! MPI
                   hr2 = zah_l(ii,ib2,1,1) ! MPI
                   hi2 = zah_l(ii,ib2,2,1) ! MPI
                   dr2 = zat_l(ii,ib2,1,1) ! MPI
                   di2 = zat_l(ii,ib2,2,1) ! MPI
                   dr1 = zat_l(ii,ib1,1,1) ! MPI
                   di1 = zat_l(ii,ib1,2,1) ! MPI
                   w1hw2(ip0) =w1hw2(ip0)+(dr1*hr2+di1*hi2)*2.d0
                   w1sw2(ip0) =w1sw2(ip0)+(dr1*dr2+di1*di2)*2.d0
                end do
!!$                if(mype == 0) then
                if(myrank_e == 0) then
                   hr2 = zah_l(1,ib2,1,1) ! MPI
                   hi2 = zah_l(1,ib2,2,1) ! MPI
                   dr2 = zat_l(1,ib2,1,1) ! MPI
                   di2 = zat_l(1,ib2,2,1) ! MPI
                   dr1 = zat_l(1,ib1,1,1) ! MPI
                   di1 = zat_l(1,ib1,2,1) ! MPI
                   w1hw2(ip0) =w1hw2(ip0)+dr1*hr2+di1*hi2
                   w1sw2(ip0) =w1sw2(ip0)+dr1*dr2+di1*di2
                end if
             else
                w1hw2(2*ip0-1:2*ip0) = 0.d0
                w1sw2(2*ip0-1:2*ip0) = 0.d0
                do ii = 1, np_g1k(ik)           ! MPI
                   hr2 = zah_l(ii,ib2,1,1) ! MPI
                   hi2 = zah_l(ii,ib2,2,1) ! MPI
                   dr2 = zat_l(ii,ib2,1,1) ! MPI
                   di2 = zat_l(ii,ib2,2,1) ! MPI
                   dr1 = zat_l(ii,ib1,1,1) ! MPI
                   di1 = zat_l(ii,ib1,2,1) ! MPI
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
    call add_nonlocal_part ! w1hw2 = w1hw2 + w1Vnlw2
                           ! w1sw2 = w1sw2 + w1qw2
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
       do ib1 = 1, neg
          if(map_e(ib1) == myrank_e) then         ! MPI
             eko_l(map_z(ib1),ik)=eko_d(ib1)
          end if
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
          do ib1=1,neg
             ib1to = neordr(ib1,ik)
             if(map_e(ib1to) == myrank_e) dr1=dr1+eko_l(map_z(ib1to),ik) ! MPI
             dr2=dr2+eig(ib1)
          enddo
          call mpi_allreduce(dr1,di1,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI
          dr1 = di1  ! MPI
          write(nfout,'(" sum of eko_l, eig, abs diff =",3e25.10)') dr1,dr2,abs(dr2-dr1)
       endif
!! (subspace rotation) !!
       call subspace_rotation ! vec,zat_l -> zat_l
!( tenchi ) (zaj_l <- zat_l)
       call m_ES_W_transpose_back(ista_k,iend_k,ik,zaj_l,zat_l(1,1,1,1))
       zaj_l_backup(:,:,:) = zaj_l(:,:,ik,:)
!! (eko_l)
       do ib1 = 1, neg
          if(map_e(ib1) == myrank_e) then         ! MPI
             eko_l(map_z(ib1),ik)=eig(ib1)
          end if
       end do
       if(iprimddavidson >= 2) then
          eko2 = sum(eig(1:neg))
          write(nfout,1201) ik,eko1,ekod,eko2

          write(nfout,*) 'eko_l'
          write(nfout,9001) (eko_l(ib1,ik),ib1=1,neg)
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
      integer :: ib1,ib2,ibb2
!!$      real(kind=DP), dimension(np_g1k_x,neg,kimg) :: zaj_wk
      real(kind=DP), allocatable, dimension(:,:,:) :: zaj_wk
      allocate(zaj_wk(np_g1k_x,neg,kimg))

      zaj_wk(:,:,:) = 0.d0
      if(kimg==1) then
         do ib1=1,neg
            do ib2=1,neg
               do ii=1,np_g1k(ik)
                  zaj_wk(ii,ib1,kimg) = zaj_wk(ii,ib1,kimg) + zat_l(ii,ib2,kimg,1)*vec(ib2,ib1)
               end do
            end do
         end do
      else
         if(k_symmetry(ik) == GAMMA) then
            do ib1=1,neg
               do ib2=1,neg
                  hr2=vec(ib2,ib1)
                  do ii=1,np_g1k(ik)
                     dr1=zat_l(ii,ib2,1   ,1)
                     di1=zat_l(ii,ib2,kimg,1)
                     zaj_wk(ii,ib1,1   ) = zaj_wk(ii,ib1,1   ) + dr1*hr2
                     zaj_wk(ii,ib1,kimg) = zaj_wk(ii,ib1,kimg) + di1*hr2
                  end do
               end do
            end do
         else
            do ib1=1,neg
               do ib2=1,neg
                  hr2=vec(2*ib2-1,ib1)
                  hi2=vec(2*ib2  ,ib1)
                  do ii=1,np_g1k(ik)
                     dr1=zat_l(ii,ib2,1   ,1)
                     di1=zat_l(ii,ib2,kimg,1)
                     zaj_wk(ii,ib1,1   ) = zaj_wk(ii,ib1,1   ) + dr1*hr2 - di1*hi2
                     zaj_wk(ii,ib1,kimg) = zaj_wk(ii,ib1,kimg) + dr1*hi2 + di1*hr2
                  end do
               end do
            end do
         end if
      end if
      zat_l(:,:,:,1) = zaj_wk(:,:,:)
      deallocate(zaj_wk)
    end subroutine subspace_rotation

  end subroutine evolve_WFs_in_subspace_sr

! =========================== added by K. Tagami ================= 11.0
  subroutine evolve_WFs_in_subspace_sr_noncl( ik, ekin, afft_kt, bfft_kt )
    integer, intent(in) :: ik
    real(kind=DP), intent(in)  :: ekin(kg1)
    real(kind=DP), intent(inout)  :: afft_kt(nfft,ndim_chgpot)
    real(kind=DP), intent(out) :: bfft_kt(nfft,ndim_spinor)

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

    integer :: is, k1

    call tstatc0_begin('evolve_WFs_in_subspace_sr_noncl(davidson) ', id_sname,1)

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


    do ib1 = 1, neg
       if(map_e(ib1) == myrank_e) eko_d(ib1) = eko_l(map_z(ib1),ik)  ! MPI
    end do
    call mpi_allreduce(eko_d,eko_d_mpi,neg,mpi_double_precision,mpi_sum &
         & ,mpi_k_world(myrank_k),ierr)       ! MPI
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

    Do is=1, ndim_spinor
       k1 = ik + is -1
       call m_ES_W_transpose( ista_k, iend_k, k1, zaj_l, &
            &                 zat_l_noncl(1,1,1,1,is))
    End do

    do ib1 = ista_e, iend_e, istep_e     ! MPI
       Do is=1, ndim_spinor
          k1 = ik +is -1
         call m_ES_WF_in_Rspace( k1, ib1, bfft_kt(:,is) )!(swffft)
       End do
       call m_FFT_Vlocal_W_noncl( afft_kt, bfft_kt, ndim_chgpot, ndim_spinor )
                                       ! (afft, bfft)-> (bfft)
       Do is=1, ndim_spinor
         call m_FFT_WF( ELECTRON, nfout, bfft_kt(:,is), DIRECT, ON )
       End do

       ib = map_z(ib1)                 ! MPI
       if (kimg == 1) then
          Do is=1, ndim_spinor
             k1 = ik + is -1
             do ii=1,iba(ik)
                i1  = igf(nbase(ii,ik))
                dr1 = zaj_l(ii,ib,k1,1)
                dr2 = bfft_kt(i1,is)*denom
                zaj_l(ii,ib,k1,1) = ekin(ii)*dr1 +dr2
             enddo
          End do
       else
          Do is=1, ndim_spinor
             k1 = ik + is -1
             do ii=1,iba(ik)
                i1  = igf(nbase(ii,ik))
                dr1  = zaj_l(ii,ib,k1,1)
                di1  = zaj_l(ii,ib,k1,kimg)
                zaj_l(ii,ib,k1,1)   = ekin(ii)*dr1 +bfft_kt(2*i1-1,is)*denom
                zaj_l(ii,ib,k1,kimg)= ekin(ii)*di1 +bfft_kt(2*i1,  is)*denom
             enddo
          End Do
       endif
    enddo
!( tenchi ) (zah_l <- zaj_l)
    Do is=1, ndim_spinor
       k1 = ik + is -1
       call m_ES_W_transpose( ista_k, iend_k, k1, zaj_l, &
            &                 zah_l_noncl(1,1,1,1,is) )
    End do

! (make matrix elements )
! parallel loop
    ! <n|T+Vloc|m> G-wise parallel
    do ib2 = 1,nsize_sb_now
       ip0b = ib2*(ib2-1)/2
       do ib1 = 1,ib2
          ip0 = ip0b + ib1
          if (kimg == 1) then
             w1hw2(ip0) = 0.d0
             w1sw2(ip0) = 0.d0

             Do is=1, ndim_spinor
                do ii = 1, np_g1k(ik)            ! MPI
                   hr2 = zah_l_noncl(ii,ib2,1,1,is)
                   dr2 = zat_l_noncl(ii,ib2,1,1,is)
                   dr1 = zat_l_noncl(ii,ib1,1,1,is)
                   w1hw2(ip0) = w1hw2(ip0) + dr1*hr2
                   w1sw2(ip0) = w1sw2(ip0) + dr1*dr2
                end do
             End Do
          else
             if(k_symmetry(ik) == GAMMA) then
                w1hw2(ip0) = 0.d0
                w1sw2(ip0) = 0.d0
!!$                ig1 = 1;  if(mype == 0) ig1 = 2
                ig1 = 1;  if(myrank_e == 0) ig1 = 2

                Do is=1, ndim_spinor
                   do ii = ig1, np_g1k(ik)            ! MPI
                      hr2 = zah_l_noncl(ii,ib2,1,1,is) ! MPI
                      hi2 = zah_l_noncl(ii,ib2,2,1,is) ! MPI
                      dr2 = zat_l_noncl(ii,ib2,1,1,is) ! MPI
                      di2 = zat_l_noncl(ii,ib2,2,1,is) ! MPI
                      dr1 = zat_l_noncl(ii,ib1,1,1,is) ! MPI
                      di1 = zat_l_noncl(ii,ib1,2,1,is) ! MPI
                      w1hw2(ip0) =w1hw2(ip0)+(dr1*hr2+di1*hi2)*2.d0
                      w1sw2(ip0) =w1sw2(ip0)+(dr1*dr2+di1*di2)*2.d0
                   end do
                End do
!!$                if(mype == 0) then

                if(myrank_e == 0) then
                   Do is=1, ndim_spinor
                      hr2 = zah_l_noncl(1,ib2,1,1,is) ! MPI
                      hi2 = zah_l_noncl(1,ib2,2,1,is) ! MPI
                      dr2 = zat_l_noncl(1,ib2,1,1,is) ! MPI
                      di2 = zat_l_noncl(1,ib2,2,1,is) ! MPI
                      dr1 = zat_l_noncl(1,ib1,1,1,is) ! MPI
                      di1 = zat_l_noncl(1,ib1,2,1,is) ! MPI
                      w1hw2(ip0) =w1hw2(ip0)+dr1*hr2+di1*hi2
                      w1sw2(ip0) =w1sw2(ip0)+dr1*dr2+di1*di2
                   End do
                end if
             else
                w1hw2(2*ip0-1:2*ip0) = 0.d0
                w1sw2(2*ip0-1:2*ip0) = 0.d0
                Do is=1, ndim_spinor
                   do ii = 1, np_g1k(ik)           ! MPI
                      hr2 = zah_l_noncl(ii,ib2,1,1,is) ! MPI
                      hi2 = zah_l_noncl(ii,ib2,2,1,is) ! MPI
                      dr2 = zat_l_noncl(ii,ib2,1,1,is) ! MPI
                      di2 = zat_l_noncl(ii,ib2,2,1,is) ! MPI
                      dr1 = zat_l_noncl(ii,ib1,1,1,is) ! MPI
                      di1 = zat_l_noncl(ii,ib1,2,1,is) ! MPI
                      w1hw2(2*ip0-1) =w1hw2(2*ip0-1)+dr1*hr2+di1*hi2
                      w1hw2(2*ip0  ) =w1hw2(2*ip0  )+dr1*hi2-di1*hr2
                      w1sw2(2*ip0-1) =w1sw2(2*ip0-1)+dr1*dr2+di1*di2
                      w1sw2(2*ip0  ) =w1sw2(2*ip0  )+dr1*di2-di1*dr2
                   end do
                End do
             end if
          end if
       end do
    end do

    if(iprimddavidson >= 3) call wd_w1hw2(" -- w1hw2 without nl part--")
    ! <n|Vnl|m> G-wise parallel
    call add_nonlocal_part_noncl  ! w1hw2 = w1hw2 + w1Vnlw2
                                  ! w1sw2 = w1sw2 + w1qw2
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
       Do is=1, ndim_spinor
          zaj_l(:,:,ik+is-1,:) = zaj_l_backup_noncl(:,:,:,is)
       End do
       do ib1 = 1, neg
          if(map_e(ib1) == myrank_e) then         ! MPI
             eko_l(map_z(ib1),ik)=eko_d(ib1)
          end if
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
          do ib1=1,neg
             ib1to = neordr(ib1,ik)
             if(map_e(ib1to) == myrank_e) dr1=dr1+eko_l(map_z(ib1to),ik) ! MPI
             dr2=dr2+eig(ib1)
          enddo
          call mpi_allreduce(dr1,di1,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI
          dr1 = di1  ! MPI
          write(nfout,'(" sum of eko_l, eig, abs diff =",3e25.10)') dr1,dr2,abs(dr2-dr1)
       endif


!! (subspace rotation) !!
       call subspace_rotation_noncl          ! vec,zat_l -> zat_l

!( tenchi ) (zaj_l <- zat_l)
       Do is=1, ndim_spinor
          k1 = ik + is -1
          call m_ES_W_transpose_back( ista_k, iend_k, k1, zaj_l, &
               &                      zat_l_noncl(1,1,1,1,is) )
          zaj_l_backup_noncl(:,:,:,is) = zaj_l(:,:,k1,:)
       End do

!! (eko_l)
       do ib1 = 1, neg
          if(map_e(ib1) == myrank_e) then         ! MPI
             eko_l(map_z(ib1),ik)=eig(ib1)
          end if
       end do
       if(iprimddavidson >= 2) then
          eko2 = sum(eig(1:neg))
          write(nfout,1201) ik,eko1,ekod,eko2

          write(nfout,*) 'eko_l'
          write(nfout,9001) (eko_l(ib1,ik),ib1=1,neg)
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

    subroutine add_nonlocal_part_noncl
      integer :: ip,ib1,ib2
      integer       :: ia, lmt1, lmt2, it, p, s, ib
      integer :: is1, is2, is_tmp
      complex(kind=CMPLDP) :: facv,facq
      real(kind=DP) :: vr,vi,qr,qi, cv1, cv2, cq1, cq2
      real(kind=DP) :: tmpr,tmpi

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
               do lmt1 = 1, ilmt(it)
                  p = lmta(lmt1,ia)
                  do lmt2 = 1, ilmt(it)
                     s = lmta(lmt2,ia)

                     Do is1=1, ndim_spinor
                        Do is2=1, ndim_spinor
                           is_tmp = 2*( is1 -1 ) +is2
                           facv = iwei(ia) *dion_scr_noncl( lmt1, lmt2, is_tmp, ia )
                           facq = iwei(ia) *q_noncl( lmt1, lmt2, is_tmp, it )

                           if ( kimg==1 ) then
                              tmpr = fsr_noncl( ib1, p, 1, is1 ) &
        &                           *fsr_noncl( ib2, s, 1, is2 ) &
        &                          + fsi_noncl( ib1, p, 1, is1 ) &
        &                           *fsi_noncl( ib2, s, 1, is2 )
                              tmpi = fsr_noncl( ib1, p, 1, is1 )&
        &                           *fsi_noncl( ib2, s, 1, is2 )&
        &                          - fsi_noncl( ib1, p, 1, is1 )&
        &                           *fsr_noncl( ib2, s, 1, is2 )
                              vr = vr + real( facv )*tmpr -aimag( facv )*tmpi
                              qr = qr + real( facq )*tmpr -aimag( facq )*tmpi
                           else
                              if ( k_symmetry(ik) == GAMMA ) then
                                 tmpr = fsr_noncl( ib1, p, 1, is1 )&
        &                              *fsr_noncl( ib2, s, 1, is2 )
                                 vr = vr + real( facv )*tmpr
                                 qr = qr + real( facq )*tmpr
                              else
                                 tmpr = fsr_noncl( ib1, p, 1, is1 )&
        &                              *fsr_noncl( ib2, s, 1, is2 )&
                    &                 + fsi_noncl( ib1, p, 1, is1 )&
        &                              *fsi_noncl( ib2, s, 1, is2 )
                                 tmpi = fsr_noncl( ib1, p, 1, is1 )&
        &                              *fsi_noncl( ib2, s, 1, is2 )&
        &                             - fsi_noncl( ib1, p, 1, is1 )&
        &                              *fsr_noncl( ib2, s, 1, is2 )

                                 cv1 = real(facv);  cv2 = aimag(facv)
                                 cq1 = real(facq);  cq2 = aimag(facq)

                                 vr = vr + cv1 *tmpr - cv2 *tmpi
                                 vi = vi + cv1 *tmpi + cv2 *tmpr
                                 qr = qr + cq1 *tmpr - cq2 *tmpi
                                 qi = qi + cq1 *tmpi + cq2 *tmpr

                              endif
                           end if
                        End do
                     End do
                  end do
               end do
            end do

            if ( kimg_t==1 ) then
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
    end subroutine add_nonlocal_part_noncl

    subroutine subspace_rotation_noncl
      integer :: ib1,ib2,ibb2, is
!!$      real(kind=DP), dimension(np_g1k_x,neg,kimg) :: zaj_wk
      real(kind=DP), allocatable, dimension(:,:,:,:) :: zaj_wk
      allocate(zaj_wk(np_g1k_x,neg,kimg,ndim_spinor))

      zaj_wk = 0.d0
      if (kimg==1) then
         do ib1=1,neg
            do ib2=1,neg
               Do is=1, ndim_spinor
                  do ii=1,np_g1k(ik)
                     zaj_wk(ii,ib1,kimg,is) = zaj_wk(ii,ib1,kimg,is) &
                          &                 + zat_l_noncl(ii,ib2,kimg,1,is)*vec(ib2,ib1)
                  end do
               End Do
            end do
         end do
      else
         if (k_symmetry(ik) == GAMMA) then
            do ib1=1,neg
               do ib2=1,neg
                  hr2=vec(ib2,ib1)
                  Do is=1, ndim_spinor
                     do ii=1,np_g1k(ik)
                        dr1 = zat_l_noncl(ii,ib2,1   ,1,is)
                        di1 = zat_l_noncl(ii,ib2,kimg,1,is)
                        zaj_wk(ii,ib1,1   ,is) = zaj_wk(ii,ib1,1   ,is) +dr1*hr2
                        zaj_wk(ii,ib1,kimg,is) = zaj_wk(ii,ib1,kimg,is) +di1*hr2
                     end do
                  End do
               end do
            end do
         else
            do ib1=1,neg
               do ib2=1,neg
                  hr2=vec(2*ib2-1,ib1)
                  hi2=vec(2*ib2  ,ib1)
                  Do is=1, ndim_spinor
                     do ii=1,np_g1k(ik)
                        dr1 = zat_l_noncl(ii,ib2,1   ,1,is)
                        di1 = zat_l_noncl(ii,ib2,kimg,1,is)
                        zaj_wk(ii,ib1,1   ,is) = zaj_wk(ii,ib1,1   ,is) +dr1*hr2 -di1*hi2
                        zaj_wk(ii,ib1,kimg,is) = zaj_wk(ii,ib1,kimg,is) +dr1*hi2 +di1*hr2
                     end do
                  End do

               end do
            end do
         end if
      end if

      Do is=1, ndim_spinor
         zat_l_noncl(:,:,:,1,is) = zaj_wk(:,:,:,is)
      End Do
      deallocate(zaj_wk)

    end subroutine subspace_rotation_noncl

  end subroutine evolve_WFs_in_subspace_sr_noncl
! =================================================================== 11.0

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
    logical :: frestart
    integer :: iblock_now, itot_now, ipri0
    integer :: n_unconv

    allocate(ekin(kg1),p(kg1))
    call m_ES_alloc_scss_etc()
    allocate(afft(nfft)); allocate(bfft(nfft))

    call allocate_matrix_ksg
    do ispin = 1, nspin, (af+1)
       call m_ES_Vlocal_in_Rspace(ispin,afft)      ! (ptfft1) vlhxc_l->afft
       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) /= myrank_k) cycle          ! MPI
          iksnl = (ik-1)/nspin + 1
          call allocate_t_matrix_ksg(ik)
          call m_pwBS_kinetic_energies(ik,vkxyz,ekin) ! (diakin) ->ekin
          feigconv = .false.
          Loop: do itot=1,max_iter_mdkosugi
             itot_now = itot
             call m_ES_Vnonlocal_W(ik,iksnl,ispin,switch_of_eko_part=ON) ! -> vnlph_l
             call decide_correction_vector_ksg(precon,ik,ekin,afft,bfft,p)   ! -> wfsd_l
             call prepare_Hloc_phi_ksg(ik,ekin,afft,bfft,itot)
!print *,ik,itot
if(iprimdkosugi>=2) write(nfout,*) 'B',feigconv
             Block_Loop: do iblock=1,nblock
                iblock_now = iblock
                call evolve_WFs_in_subspace_ksg &
                          (ik,ispin,ekin,afft,bfft,iblock,itot,frestart)
             end do Block_Loop
if(iprimdkosugi>=2) write(nfout,*) 'A',feigconv
!             call m_ES_betar_dot_WFs_4_each_k(nfout,ik)   ! -> fsr_l,fsi_l
          end do Loop
          if(iprimdkosugi>=2) then
             write(nfout,'("MdKosugi: ik=",i5," itot=",i5," subspace=",i5)') &
                                                                 ik, itot_now, nsize_sb_now
          end if
          call deallocate_t_matrix_ksg
       enddo      ! k-point loop
    enddo      ! spin loop
    call deallocate_matrix_ksg
!    if(iprimddavidson>=2) then
!       write(nfout,'("Modified Davidson: max_itot=",i5)') max_itot
!    end if
!
    call m_ES_sort_eigen_values()
!!!  ( in case of af=1 )
!    if(af /= 0) then
!       call cp_eigen_values_for_af       !-(contained here)
!       call expand_neordr_and_nrvf_ordr  !-(contained here)
!    end if

    call get_ipri0(iprimdkosugi,ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko(nfout,mode=SCF)
!
    deallocate(bfft);   deallocate(afft)
    call m_ES_dealloc_scss_etc()
    deallocate(ekin,p)
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

! ============================== added by K. Tagami ================ 11.0
  subroutine m_ESmdkosugi_Renew_WF_noncl(nfout,precon)
    integer, intent(in) :: nfout,precon
    integer             :: ispin, ik, iksnl, switch_of_eko_part
    integer :: iblock,itot

    real(kind=DP), allocatable ::  afft_kt(:,:)
    real(kind=DP), allocatable ::  bfft_kt(:,:)
    real(kind=DP), allocatable :: vnlph_noncl(:,:,:,:)
    real(kind=DP), allocatable, dimension(:) :: ekin,p

    logical :: frestart
    integer :: iblock_now, itot_now, ipri0
    integer :: n_unconv

    integer :: is1, is2, k2, istmp

    allocate(ekin(kg1),p(kg1))
    call m_ES_alloc_scss_etc()
    allocate(afft_kt(nfft,ndim_chgpot)); afft_kt = 0.0d0
    allocate(bfft_kt(nfft,ndim_spinor)); bfft_kt = 0.0d0
    allocate(vnlph_noncl(kg1,np_e,kimg,ndim_spinor)); vnlph_noncl = 0.0d0

    call allocate_matrix_ksg_noncl

    call m_ES_Vlocal_in_Rspace_noncl( afft_kt )    ! (ptfft1) vlhxc_l->afft

    Do ik=1, kv3, ndim_spinor
       if ( map_k(ik) /= myrank_k )  cycle          ! MPI
       iksnl = (ik-1)/ndim_spinor + 1

       call allocate_t_matrix_ksg_noncl( ik,ik+ndim_spinor-1 ) ! -> np_g1k_x
       call m_pwBS_kinetic_energies( ik, vkxyz, ekin ) ! (diakin) ->ekin

       feigconv = .false.
       Loop: do itot=1,max_iter_mdkosugi
          itot_now = itot

          vnlph_noncl = 0.0d0
          Do is1=1, ndim_spinor
             Do is2=1, ndim_spinor
                istmp = ( is1-1 )*ndim_spinor +is2
                k2 = ik +is2 -1

                if ( precon==ON ) then
                   call m_ES_Vnonlocal_W( k2, iksnl, istmp, switch_of_eko_part=OFF )
                   ! -> vnlph_l
                else
                   call m_ES_Vnonlocal_W( k2, iksnl, istmp, switch_of_eko_part=ON )
                   ! -> vnlph_l
                endif
                vnlph_noncl(:,:,:,is1) = vnlph_noncl(:,:,:,is1) &
                     &                     + vnlph_l(:,:,:)
             End do
          End do

          call decide_correction_vec_ksg_noncl( precon, ik, ekin, &
               &                                afft_kt, bfft_kt, vnlph_noncl, p )
                                                                    ! -> wfsd_l
          call prepare_Hloc_phi_ksg_noncl( ik, ekin, afft_kt, bfft_kt, itot )

          Block_Loop: do iblock=1,nblock
             iblock_now = iblock
             call evolve_WFs_in_subsp_ksg_noncl( ik, ekin, afft_kt, bfft_kt, &
                  &                              iblock, itot, frestart )
          end do Block_Loop

!             call m_ES_betar_dot_WFs_4_each_k(nfout,ik)   ! -> fsr_l,fsi_l
       end do Loop

       if (iprimdkosugi>=2) then
          write(nfout,'("MdKosugi: ik=",i5," itot=",i5," subspace=",i5)') &
               ik, itot_now, nsize_sb_now
       end if
       call deallocate_t_matrix_ksg

    enddo      ! k-point loop

    call deallocate_matrix_ksg_noncl

!    if(iprimddavidson>=2) then
!       write(nfout,'("Modified Davidson: max_itot=",i5)') max_itot
!    end if
!
    call m_ES_sort_eigen_vals_noncl()

    call get_ipri0(iprimdkosugi,ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko(nfout,mode=SCF)
!
    deallocate(bfft_kt);   deallocate(afft_kt)
    call m_ES_dealloc_scss_etc()
    deallocate(ekin,p)

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

  end subroutine m_ESmdkosugi_Renew_WF_noncl
! ======================================================================= 11.0

  subroutine allocate_matrix_ksg
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
allocate(zajold_l(kg1,np_e,kimg))
allocate(fsrold_l(np_e,nlmta))
if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
  allocate(fsiold_l(np_e,nlmta))
end if
  end subroutine allocate_matrix_ksg

! =============================== added by K. Tagami ====================== 11.0
  subroutine allocate_matrix_ksg_noncl
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

    allocate(fsr_noncl(np_e,nlmta,max_iter_mdkosugi+1,ndim_spinor))
    if(.not.(kv3/ndim_spinor == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       allocate(fsi_noncl(np_e,nlmta,max_iter_mdkosugi+1,ndim_spinor))
    end if

    eps_residual = eps_residual_mdkosugi

    allocate(zajold_l_noncl(kg1,np_e,kimg,ndim_spinor))
    allocate(fsrold_l_noncl(np_e,nlmta,ndim_spinor))
    if(.not.(kv3/ndim_spinor == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       allocate(fsiold_l_noncl(np_e,nlmta,ndim_spinor))
    end if
  end subroutine allocate_matrix_ksg_noncl
! ==============================================================  11.0


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

  subroutine allocate_t_matrix_ksg(ik)
    integer, intent(in) :: ik
    integer :: kimg_t
    if(k_symmetry(ik) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = kimg
    end if
    allocate(zat_l(kg1,np_e,kimg,max_iter_mdkosugi+1)) ! MPI
    allocate(zah_l(kg1,np_e,kimg,max_iter_mdkosugi+1)) ! MPI
    allocate(w1hw2(msize_matrix*kimg_t))
    allocate(w1sw2(msize_matrix*kimg_t))
    allocate(wfsd_l(kg1,np_e,ik:ik,kimg)); wfsd_l = 0.d0
    allocate(bsdr_l(np_e,nlmta,ik:ik)); bsdr_l = 0.d0
    allocate(bsdi_l(np_e,nlmta,ik:ik)); bsdi_l = 0.d0
zajold_l(:,:,:)=zaj_l(:,:,ik,:)
fsrold_l(:,:)=fsr_l(:,:,ik)
if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
  fsiold_l(:,:)=fsi_l(:,:,ik)
end if

  end subroutine allocate_t_matrix_ksg

! ============================= added by K. Tagami ===================== 11.0
  subroutine allocate_t_matrix_ksg_noncl( ik1, ik2 )
    integer, intent(in) :: ik1, ik2
    integer :: kimg_t
    integer :: is, k1

    if(k_symmetry(ik1) == GAMMA) then
       kimg_t = 1
    else
       kimg_t = kimg
    end if
    allocate(zat_l_noncl(kg1,np_e,kimg,max_iter_mdkosugi+1,ndim_spinor)) ! MPI
    allocate(zah_l_noncl(kg1,np_e,kimg,max_iter_mdkosugi+1,ndim_spinor)) ! MPI

    allocate(w1hw2(msize_matrix*kimg_t))
    allocate(w1sw2(msize_matrix*kimg_t))

    allocate(wfsd_l(kg1,np_e,ik1:ik2,kimg)); wfsd_l = 0.d0
    allocate(bsdr_l(np_e,nlmta,ik1:ik2)); bsdr_l = 0.d0
    allocate(bsdi_l(np_e,nlmta,ik1:ik2)); bsdi_l = 0.d0

    Do is=1, ndim_spinor
       k1 = ik1 + is -1
       zajold_l_noncl(:,:,:,is)= zaj_l(:,:,k1,:)
       fsrold_l_noncl(:,:,is)  = fsr_l(:,:,k1)

       if (.not.(kv3/ndim_spinor == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
          fsiold_l_noncl(:,:,is) = fsi_l(:,:,k1)
       end if
    End do

  end subroutine allocate_t_matrix_ksg_noncl
! ====================================================================== 11.0


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

  subroutine prepare_Hloc_phi_ksg(ik,ekin,afft,bfft,itot)
    integer, intent(in) :: ik
    integer, intent(in) :: itot
    real(kind=DP), intent(in)  :: ekin(kg1)
    real(kind=DP), intent(in)  :: afft(nfft)
    real(kind=DP), intent(out) :: bfft(nfft)

    integer       :: ib1,i1,ii,ib
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
      do ib1 = ista_e, iend_e, istep_e     ! MPI
         call m_ES_WF_in_Rspace(ik,ib1,bfft)!(swffft)
         call m_FFT_Vlocal_W(afft,bfft)    ! (afft, bfft)-> (bfft)
         call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
         ib = map_z(ib1)                 ! MPI
         if(kimg == 1) then
            do ii=1,iba(ik)
               i1  = igf(nbase(ii,ik))
               dr1 = zaj_l(ii,ib,ik,1)
               dr2 = bfft(i1)*denom
  !             zaj_l(ii,ib,ik,1)= ekin(ii)*dr1+dr2
               zah_l(ii,ib,1,1) = ekin(ii)*dr1+dr2
            enddo
         else
            do ii=1,iba(ik)
               i1  = igf(nbase(ii,ik))
               dr1  = zaj_l(ii,ib,ik,1)
               di1  = zaj_l(ii,ib,ik,kimg)
               zah_l(ii,ib,1,1)= ekin(ii)*dr1+bfft(2*i1-1)*denom
               zah_l(ii,ib,kimg,1)= ekin(ii)*di1+bfft(2*i1)*denom
            enddo
         endif
      enddo
    end if

!!( tenchi ) (zah_l <- zaj_l)
    zaj_l(:,:,ik,:) = wfsd_l(:,:,ik,:)

    do ib1 = ista_e, iend_e, istep_e     ! MPI
       call m_ES_WF_in_Rspace(ik,ib1,bfft)!(swffft)
       call m_FFT_Vlocal_W(afft,bfft)    ! (afft, bfft)-> (bfft)
       call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
       ib = map_z(ib1)                 ! MPI
       if(kimg == 1) then
          do ii=1,iba(ik)
             i1  = igf(nbase(ii,ik))
             dr1 = zaj_l(ii,ib,ik,1)
             dr2 = bfft(i1)*denom
!             zaj_l(ii,ib,ik,1)= ekin(ii)*dr1+dr2
             zah_l(ii,ib,1,itot+1) = ekin(ii)*dr1+dr2
          enddo
       else
          do ii=1,iba(ik)
             i1  = igf(nbase(ii,ik))
             dr1  = zaj_l(ii,ib,ik,1)
             di1  = zaj_l(ii,ib,ik,kimg)
!             zaj_l(ii,ib,ik,1)= ekin(ii)*dr1+bfft(2*i1-1)*denom
!             zaj_l(ii,ib,ik,kimg)= ekin(ii)*di1+bfft(2*i1)*denom
             zah_l(ii,ib,1,itot+1)= ekin(ii)*dr1+bfft(2*i1-1)*denom
             zah_l(ii,ib,kimg,itot+1)= ekin(ii)*di1+bfft(2*i1)*denom
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

  end subroutine prepare_Hloc_phi_ksg

! ============================= aded by K. Tagami ==================== 11.0
  subroutine prepare_Hloc_phi_ksg_noncl( ik, ekin, afft_kt, bfft_kt, itot )
    integer, intent(in) :: ik
    integer, intent(in) :: itot
    real(kind=DP), intent(in)  :: ekin(kg1)
    real(kind=DP), intent(inout)  :: afft_kt(nfft,ndim_chgpot)
    real(kind=DP), intent(out) :: bfft_kt(nfft,ndim_spinor)

    integer       :: ib1,i1,ii,ib
    real(kind=DP) :: denom
    real(kind=DP) :: dr1,dr2,di1,di2,dd

    integer :: is, k1
    integer :: id_sname = -1, ipri0

    call tstatc0_begin('prepare_Hloc_phi_noncl (mdkosugi) ', id_sname,1)

    call get_ipri0(iprimdkosugi,ipri0)

    denom = 1.d0/product(fft_box_size_WF(1:3,1))

! (zaj_l <- (T+Vloc)|phi> )
!!    zaj_l(:,:,ik,:) = zajold_l(:,:,:,idavid)
!( tenchi ) (zat_l <- zaj_l)

    Do is=1, ndim_spinor
!!       zat_l_noncl(:,:,:,1,is) = zaj_l(:,:,ik+is-1,:)
       zat_l_noncl(:,:,:,itot+1,is) = wfsd_l(:,:,ik+is-1,:)
    End do

    Do is=1, ndim_spinor
!       fsr_noncl(:,:,1,is) = fsr_l(:,:,ik+is-1)
       fsr_noncl(:,:,itot+1,is)=bsdr_l(:,:,ik+is-1)
    End do

    if(.not.(kv3/ndim_spinor == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       Do is=1, ndim_spinor
!          fsi_noncl(:,:,1,is) = fsi_l(:,:,ik+is-1)
          fsi_noncl(:,:,itot+1,is) = bsdi_l(:,:,ik+is-1)
       End do
    end if

    if(itot == 1) then
       do ib1 = ista_e, iend_e, istep_e     ! MPI
          Do is=1, ndim_spinor
             k1 = ik +is- 1
             call m_ES_WF_in_Rspace( k1, ib1, bfft_kt(:,is) )   !(swffft)
          End do
          call m_FFT_Vlocal_W_noncl( afft_kt, bfft_kt, ndim_chgpot, ndim_spinor )
                                       ! (afft, bfft)-> (bfft)
          Do is=1, ndim_spinor
             call m_FFT_WF( ELECTRON, nfout, bfft_kt(:,is), DIRECT, ON )
          End do

         ib = map_z(ib1)                 ! MPI
         if (kimg == 1) then
            Do is=1, ndim_spinor
               k1 = ik + is -1
               do ii=1,iba(ik)
                  i1  = igf(nbase(ii,ik))
                  dr1 = zaj_l(ii,ib,k1,1)
                  dr2 = bfft_kt(i1,is)*denom
                  zah_l_noncl(ii,ib,1,1,is) = ekin(ii)*dr1 +dr2
               enddo
            End do
         else
            Do is=1, ndim_spinor
               k1 = ik + is -1
               do ii=1,iba(ik)
                  i1  = igf(nbase(ii,ik))
                  dr1  = zaj_l(ii,ib,k1,1)
                  di1  = zaj_l(ii,ib,k1,kimg)
                  zah_l_noncl(ii,ib,1,   1,is) = ekin(ii)*dr1 &
                       &                        +bfft_kt(2*i1-1,is)*denom
                  zah_l_noncl(ii,ib,kimg,1,is) = ekin(ii)*di1 &
                       &                        +bfft_kt(2*i1,  is)*denom
               enddo
            End do
         endif
      enddo
    end if

!!( tenchi ) (zah_l <- zaj_l)
    Do is=1, ndim_spinor
       k1 = ik + is -1
       zaj_l(:,:,k1,:) = wfsd_l(:,:,k1,:)
    End do

    do ib1 = ista_e, iend_e, istep_e     ! MPI
       Do is=1, ndim_spinor
          k1 = ik +is- 1
         call m_ES_WF_in_Rspace( k1, ib1, bfft_kt(:,is) )!(swffft)
       End do
       call m_FFT_Vlocal_W_noncl( afft_kt, bfft_kt, ndim_chgpot, ndim_spinor )
                                       ! (afft, bfft)-> (bfft)
       Do is=1, ndim_spinor
         call m_FFT_WF( ELECTRON, nfout, bfft_kt(:,is), DIRECT, ON )
       End do

       ib = map_z(ib1)                 ! MPI
       if (kimg == 1) then
          Do is=1, ndim_spinor
             k1 = ik + is -1
             do ii=1,iba(ik)
                i1  = igf(nbase(ii,ik))
                dr1 = zaj_l(ii,ib,k1,1)
                dr2 = bfft_kt(i1,is) *denom
                zah_l_noncl(ii,ib,1,itot+1,is) = ekin(ii)*dr1 +dr2
             enddo
          End do
       else
          Do is=1, ndim_spinor
             k1 = ik + is -1
             do ii=1,iba(ik)
                i1  = igf(nbase(ii,ik))
                dr1  = zaj_l(ii,ib,k1,1)
                di1  = zaj_l(ii,ib,k1,kimg)
                zah_l_noncl(ii,ib,1,   itot+1,is) = ekin(ii)*dr1 &
                     &                             +bfft_kt(2*i1-1,is) *denom
                zah_l_noncl(ii,ib,kimg,itot+1,is) = ekin(ii)*di1 &
                     &                             +bfft_kt(2*i1,  is) *denom
             enddo
          End do
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

  end subroutine prepare_Hloc_phi_ksg_noncl
! =========================================================== 11.0

  subroutine evolve_WFs_in_subspace_ksg &
                                (ik,ispin,ekin,afft,bfft,iblock,itot,frestart)
    integer, intent(in) :: ik,ispin
    integer, intent(in) :: iblock,itot
    logical, intent(out) :: frestart
    real(kind=DP), intent(in)  :: ekin(kg1)
    real(kind=DP), intent(in)  :: afft(nfft)
    real(kind=DP), intent(out) :: bfft(nfft)
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
             do ii = 1, iba(ik)            ! MPI
                hr2 = zah_l(ii,ii2,1,iter2)
                dr2 = zat_l(ii,ii2,1,iter2)
                dr1 = zat_l(ii,ii1,1,iter1)
                w1hw2(ip0) = w1hw2(ip0) + dr1*hr2
                w1sw2(ip0) = w1sw2(ip0) + dr1*dr2
             end do
          else
             if(k_symmetry(ik) == GAMMA) then
                w1hw2(ip0) = 0.d0
                w1sw2(ip0) = 0.d0
!!$                ig1 = 1;  if(mype == 0) ig1 = 2
!!                ig1 = 1;  if(myrank_e == 0) ig1 = 2
                ig1 = 2
                do ii = ig1, iba(ik)            ! MPI
                   hr2 = zah_l(ii,ii2,1,iter2) ! MPI
                   hi2 = zah_l(ii,ii2,2,iter2) ! MPI
                   dr2 = zat_l(ii,ii2,1,iter2) ! MPI
                   di2 = zat_l(ii,ii2,2,iter2) ! MPI
                   dr1 = zat_l(ii,ii1,1,iter1) ! MPI
                   di1 = zat_l(ii,ii1,2,iter1) ! MPI
                   w1hw2(ip0) =w1hw2(ip0)+(dr1*hr2+di1*hi2)*2.d0
                   w1sw2(ip0) =w1sw2(ip0)+(dr1*dr2+di1*di2)*2.d0
                end do
!!$                if(mype == 0) then
!!                if(myrank_e == 0) then
                hr2 = zah_l(1,ii2,1,iter2) ! MPI
                hi2 = zah_l(1,ii2,2,iter2) ! MPI
                dr2 = zat_l(1,ii2,1,iter2) ! MPI
                di2 = zat_l(1,ii2,2,iter2) ! MPI
                dr1 = zat_l(1,ii1,1,iter1) ! MPI
                di1 = zat_l(1,ii1,2,iter1) ! MPI
                w1hw2(ip0) =w1hw2(ip0)+dr1*hr2+di1*hi2
                w1sw2(ip0) =w1sw2(ip0)+dr1*dr2+di1*di2
!!                end if
             else
                w1hw2(2*ip0-1:2*ip0) = 0.d0
                w1sw2(2*ip0-1:2*ip0) = 0.d0
                do ii = 1, iba(ik)           ! MPI
                   hr2 = zah_l(ii,ii2,1,iter2) ! MPI
                   hi2 = zah_l(ii,ii2,2,iter2) ! MPI
                   dr2 = zat_l(ii,ii2,1,iter2) ! MPI
                   di2 = zat_l(ii,ii2,2,iter2) ! MPI
                   dr1 = zat_l(ii,ii1,1,iter1) ! MPI
                   di1 = zat_l(ii,ii1,2,iter1) ! MPI
                   w1hw2(2*ip0-1) =w1hw2(2*ip0-1)+dr1*hr2+di1*hi2
                   w1hw2(2*ip0  ) =w1hw2(2*ip0  )+dr1*hi2-di1*hr2
                   w1sw2(2*ip0-1) =w1sw2(2*ip0-1)+dr1*dr2+di1*di2
                   w1sw2(2*ip0  ) =w1sw2(2*ip0  )+dr1*di2-di1*dr2
                end do
             end if
          end if
       end do
    end do
    if(iprimdkosugi >= 3) call wd_w1hw2(" -- w1hw2 without nl part--",iblock)
    ! <n|Vnl|m>
    call add_nonlocal_part ! w1hw2 = w1hw2 + w1Vnlw2
                           ! w1sw2 = w1sw2 + w1qw2
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
      integer :: ib1,ib2,ibb2
!!$      real(kind=DP), dimension(np_g1k_x,neg,kimg) :: zaj_wk
      real(kind=DP), allocatable, dimension(:,:,:) :: zaj_wk
      real(kind=DP), allocatable, dimension(:,:,:) :: zah_wk
      real(kind=DP), allocatable, dimension(:,:)   :: fsr_wk
      real(kind=DP), allocatable, dimension(:,:)   :: fsi_wk
      allocate(zaj_wk(kg1,nel,kimg))
      if(itot /= max_iter_mdkosugi) then
        allocate(zah_wk(kg1,nel,kimg))
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
               do ii=1,iba(ik)
!                  zaj_wk(ii,ib1,kimg) = zaj_wk(ii,ib1,kimg) + zat_l(ii,ii2,kimg,iter2)*vec(ib2,ib1)
                  zaj_wk(ii,ib1,kimg) = zaj_wk(ii,ib1,kimg) + zat_l(ii,ii2,kimg,iter2)*hr2
               end do
               do ii=1,nlmta
                 fsr_wk(ib1,ii) = fsr_wk(ib1,ii) + fsr(ii2,ii,iter2)*hr2
                 fsi_wk(ib1,ii) = fsi_wk(ib1,ii) + fsi(ii2,ii,iter2)*hr2
               end do
               if(itot /= max_iter_mdkosugi) then
                 do ii=1,iba(ik)
                    zah_wk(ii,ib1,kimg) = zah_wk(ii,ib1,kimg) + zah_l(ii,ii2,kimg,iter2)*hr2
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
                  do ii=1,iba(ik)
                     dr1=zat_l(ii,ii2,1   ,iter2)
                     di1=zat_l(ii,ii2,kimg,iter2)
                     zaj_wk(ii,ib1,1   ) = zaj_wk(ii,ib1,1   ) + dr1*hr2
                     zaj_wk(ii,ib1,kimg) = zaj_wk(ii,ib1,kimg) + di1*hr2
                  end do
                  do ii=1,nlmta
                    fsr_wk(ib1,ii) = fsr_wk(ib1,ii) + fsr(ii2,ii,iter2)*hr2
                  end do
                  if(itot /= max_iter_mdkosugi) then
                    do ii=1,iba(ik)
                       dr1=zah_l(ii,ii2,1   ,iter2)
                       di1=zah_l(ii,ii2,kimg,iter2)
                       zah_wk(ii,ib1,1   ) = zah_wk(ii,ib1,1   ) + dr1*hr2
                       zah_wk(ii,ib1,kimg) = zah_wk(ii,ib1,kimg) + di1*hr2
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
                  do ii=1,iba(ik)
                     dr1=zat_l(ii,ii2,1   ,iter2)
                     di1=zat_l(ii,ii2,kimg,iter2)
                     zaj_wk(ii,ib1,1   ) = zaj_wk(ii,ib1,1   ) + dr1*hr2 - di1*hi2
                     zaj_wk(ii,ib1,kimg) = zaj_wk(ii,ib1,kimg) + dr1*hi2 + di1*hr2
                  end do
                  do ii=1,nlmta
                     dr1=fsr(ii2,ii,iter2)
                     di1=fsi(ii2,ii,iter2)
                     fsr_wk(ib1,ii) = fsr_wk(ib1,ii) + dr1*hr2 - di1*hi2
                     fsi_wk(ib1,ii) = fsi_wk(ib1,ii) + dr1*hi2 + di1*hr2
                  end do
                  if(itot /= max_iter_mdkosugi) then
                    do ii=1,iba(ik)
                       dr1=zah_l(ii,ii2,1   ,iter2)
                       di1=zah_l(ii,ii2,kimg,iter2)
                       zah_wk(ii,ib1,1   ) = zah_wk(ii,ib1,1   ) + dr1*hr2 - di1*hi2
                       zah_wk(ii,ib1,kimg) = zah_wk(ii,ib1,kimg) + dr1*hi2 + di1*hr2
                    end do
                  end if
               end do
            end do
         end if
      end if
!print *,itot,  itot /= max_iter_mdkosugi
      zaj_l(:,nsta:nend,ik,:) = zaj_wk(:,:,:)
      if(itot /= max_iter_mdkosugi) zah_l(:,nsta:nend,:,1) = zah_wk(:,:,:)
      fsr_l(nsta:nend,:,ik) = fsr_wk(:,:)
      if(k_symmetry(ik) /= GAMMA) fsi_l(nsta:nend,:,ik) = fsi_wk(:,:)

      deallocate(zaj_wk,fsr_wk)
      if(itot /= max_iter_mdkosugi) deallocate(zah_wk)
      if(k_symmetry(ik) /= GAMMA) deallocate(fsi_wk)
    end subroutine subspace_rotation

  end subroutine evolve_WFs_in_subspace_ksg

! ================================ added by K. Tagami ============ 11.0
  subroutine evolve_WFs_in_subsp_ksg_noncl( ik, ekin, afft_kt, bfft_kt, &
       &                                    iblock, itot,frestart )
    integer, intent(in) :: ik
    integer, intent(in) :: iblock,itot
    logical, intent(out) :: frestart
    real(kind=DP), intent(in)  :: ekin(kg1)
    real(kind=DP), intent(inout)  :: afft_kt(nfft,ndim_chgpot)
    real(kind=DP), intent(out) :: bfft_kt(nfft,ndim_spinor)

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

    integer :: is, k1

    integer :: id_sname = -1, ipri0

    call tstatc0_begin('evolve_WFs_in_subspace_noncl (modified kosugi) ', id_sname,1)

! ==================== added by K. Tagami ========== 11.0
#ifdef forsafe
    ekod = 0.0d0
#endif
! ================================================== 11.0

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

          if (kimg == 1) then
             w1hw2(ip0) = 0.d0
             w1sw2(ip0) = 0.d0

             Do is=1, ndim_spinor
                do ii = 1, iba(ik)            ! MPI
                   hr2 = zah_l_noncl(ii,ii2,1,iter2,is)
                   dr2 = zat_l_noncl(ii,ii2,1,iter2,is)
                   dr1 = zat_l_noncl(ii,ii1,1,iter1,is)
                   w1hw2(ip0) = w1hw2(ip0) + dr1*hr2
                   w1sw2(ip0) = w1sw2(ip0) + dr1*dr2
                end do
             End do
          else
             if(k_symmetry(ik) == GAMMA) then
                w1hw2(ip0) = 0.d0
                w1sw2(ip0) = 0.d0
!!$                ig1 = 1;  if(mype == 0) ig1 = 2
!!                ig1 = 1;  if(myrank_e == 0) ig1 = 2
                ig1 = 2

                Do is=1, ndim_spinor
                   do ii = ig1, iba(ik)            ! MPI
                      hr2 = zah_l_noncl(ii,ii2,1,iter2,is) ! MPI
                      hi2 = zah_l_noncl(ii,ii2,2,iter2,is) ! MPI
                      dr2 = zat_l_noncl(ii,ii2,1,iter2,is) ! MPI
                      di2 = zat_l_noncl(ii,ii2,2,iter2,is) ! MPI
                      dr1 = zat_l_noncl(ii,ii1,1,iter1,is) ! MPI
                      di1 = zat_l_noncl(ii,ii1,2,iter1,is) ! MPI
                      w1hw2(ip0) =w1hw2(ip0)+(dr1*hr2+di1*hi2)*2.d0
                      w1sw2(ip0) =w1sw2(ip0)+(dr1*dr2+di1*di2)*2.d0
                   end do
                End do

!!$                if(mype == 0) then
!!                if(myrank_e == 0) then

                Do is=1, ndim_spinor
                   hr2 = zah_l_noncl(1,ii2,1,iter2,is) ! MPI
                   hi2 = zah_l_noncl(1,ii2,2,iter2,is) ! MPI
                   dr2 = zat_l_noncl(1,ii2,1,iter2,is) ! MPI
                   di2 = zat_l_noncl(1,ii2,2,iter2,is) ! MPI
                   dr1 = zat_l_noncl(1,ii1,1,iter1,is) ! MPI
                   di1 = zat_l_noncl(1,ii1,2,iter1,is) ! MPI
                   w1hw2(ip0) =w1hw2(ip0)+dr1*hr2+di1*hi2
                   w1sw2(ip0) =w1sw2(ip0)+dr1*dr2+di1*di2
                End Do
                !!                end if
             else
                w1hw2(2*ip0-1:2*ip0) = 0.d0
                w1sw2(2*ip0-1:2*ip0) = 0.d0
                Do is=1, ndim_spinor
                   do ii = 1, iba(ik)           ! MPI
                      hr2 = zah_l_noncl(ii,ii2,1,iter2,is) ! MPI
                      hi2 = zah_l_noncl(ii,ii2,2,iter2,is) ! MPI
                      dr2 = zat_l_noncl(ii,ii2,1,iter2,is) ! MPI
                      di2 = zat_l_noncl(ii,ii2,2,iter2,is) ! MPI
                      dr1 = zat_l_noncl(ii,ii1,1,iter1,is) ! MPI
                      di1 = zat_l_noncl(ii,ii1,2,iter1,is) ! MPI
                      w1hw2(2*ip0-1) =w1hw2(2*ip0-1)+dr1*hr2+di1*hi2
                      w1hw2(2*ip0  ) =w1hw2(2*ip0  )+dr1*hi2-di1*hr2
                      w1sw2(2*ip0-1) =w1sw2(2*ip0-1)+dr1*dr2+di1*di2
                      w1sw2(2*ip0  ) =w1sw2(2*ip0  )+dr1*di2-di1*dr2
                   end do
                End do

             end if
          end if
       end do
    end do

    if(iprimdkosugi >= 3) call wd_w1hw2(" -- w1hw2 without nl part--",iblock)
    ! <n|Vnl|m>
    call add_nonlocal_part_noncl ! w1hw2 = w1hw2 + w1Vnlw2
                                ! w1sw2 = w1sw2 + w1qw2
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

    if(kimg_t == 1) then
       call dspgvx_driver_loc(eig,vec,w1hw2,w1sw2,ierr_diag,nel)
    else
       call zhpgvx_driver(eig,vec,w1hw2,w1sw2,ierr_diag,nel)
    endif

    frestart = .false.
    if (ierr_diag /= 0) then
!       zaj_l(:,:,ik,:) = zaj_l_backup(:,:,:)
!       do ib1 = 1, neg
!          if(map_e(ib1) == myrank_e) then         ! MPI
!             eko_l(map_z(ib1),ik)=eko_d(ib1)
!          end if
!       end do
       frestart = .true.

       Do is=1, ndim_spinor
          k1 = ik + is -1
          zaj_l(:,nsta:nend,k1,:)=zat_l_noncl(:,nsta:nend,:,1,is)
       End do
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

          write(nfout,'(" sum of eko_l, eig, abs diff =",3e15.10)') dr1,dr2,abs(dr2-dr1)
       endif


!!! (subspace rotation) !!
       call subspace_rotation_noncl            ! vec,zat_l -> zat_l


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
1201   format(' %% for ik = ',i4,4x,' eko1&ekod&eko2 = ',3f14.7)
9001   format(5x,6f12.5)
!!! (neordr & nrvf_ordr)
!
    end if
9000 continue
!    neordr(1:neg,ik) = (/(ib1,ib1=1,neg)/)
!    nrvf_ordr(1:neg,ik) = (/(ib1,ib1=1,neg)/)
!
!! (deallocate)
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

    subroutine add_nonlocal_part_noncl
      integer :: ip,ib1,ib2,ibb1,ibb2
      integer       :: ia, lmt1, lmt2, it, p, s, ib

      integer :: is1, is2, is_tmp

      complex(kind=CMPLDP) :: facv,facq
      real(kind=DP) :: vr,vi,qr,qi, cv1, cv2, cq1, cq2
      real(kind=DP) :: tmpr,tmpi

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

            if (kimg_t==1) then
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


                     Do is1=1, ndim_spinor
                        Do is2=1, ndim_spinor
                           is_tmp = 2*( is1 -1 ) +is2
                           facv = iwei(ia) *dion_scr_noncl( lmt1, lmt2, is_tmp, ia )
                           facq = iwei(ia) *q_noncl( lmt1, lmt2, is_tmp, it )

                           if ( kimg==1 ) then
                              tmpr = fsr_noncl( ii1, p, iter1, is1 ) &
        &                           *fsr_noncl( ii2, s, iter2, is2 ) &
        &                          + fsi_noncl( ii1, p, iter1, is1 ) &
        &                           *fsi_noncl( ii2, s, iter2, is2 )
                              tmpi = fsr_noncl( ii1, p, iter1, is1 )&
        &                           *fsi_noncl( ii2, s, iter2, is2 )&
        &                          - fsi_noncl( ii1, p, iter1, is1 )&
        &                           *fsr_noncl( ii2, s, iter2, is2 )
                              vr = vr + real( facv )*tmpr -aimag( facv )*tmpi
                              qr = qr + real( facq )*tmpr -aimag( facq )*tmpi
                           else
                              if ( k_symmetry(ik) == GAMMA ) then
                                 tmpr = fsr_noncl( ii1, p, iter1, is1 )&
        &                              *fsr_noncl( ii2, s, iter2, is2 )
                                 vr = vr + real( facv )*tmpr
                                 qr = qr + real( facq )*tmpr
                              else
                                 tmpr = fsr_noncl( ii1, p, iter1, is1 )&
        &                              *fsr_noncl( ii2, s, iter2, is2 )&
                    &                 + fsi_noncl( ii1, p, iter1, is1 )&
        &                              *fsi_noncl( ii2, s, iter2, is2 )
                                 tmpi = fsr_noncl( ii1, p, iter1, is1 )&
        &                              *fsi_noncl( ii2, s, iter2, is2 )&
        &                             - fsi_noncl( ii1, p, iter1, is1 )&
        &                              *fsr_noncl( ii2, s, iter2, is2 )

                                 cv1 = real(facv);  cv2 = aimag(facv)
                                 cq1 = real(facq);  cq2 = aimag(facq)

                                 vr = vr + cv1 *tmpr - cv2 *tmpi
                                 vi = vi + cv1 *tmpi + cv2 *tmpr
                                 qr = qr + cq1 *tmpr - cq2 *tmpi
                                 qi = qi + cq1 *tmpi + cq2 *tmpr

                              endif
                           end if
                        End do
                     End do

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
    end subroutine add_nonlocal_part_noncl

    subroutine subspace_rotation_noncl
      integer :: ib1,ib2,ibb2
!!$      real(kind=DP), dimension(np_g1k_x,neg,kimg) :: zaj_wk
      real(kind=DP), allocatable, dimension(:,:,:,:) :: zaj_wk
      real(kind=DP), allocatable, dimension(:,:,:,:) :: zah_wk
      real(kind=DP), allocatable, dimension(:,:,:)   :: fsr_wk
      real(kind=DP), allocatable, dimension(:,:,:)   :: fsi_wk

      integer :: is

      allocate(zaj_wk(kg1,nel,kimg,ndim_spinor))
      if(itot /= max_iter_mdkosugi) then
        allocate(zah_wk(kg1,nel,kimg,ndim_spinor))
        zah_wk = 0.d0
      end if
      allocate(fsr_wk(nel,nlmta,ndim_spinor))
      if (k_symmetry(ik) /= GAMMA) then
        allocate(fsi_wk(nel,nlmta,ndim_spinor))
        fsi_wk = 0.d0
      end if

      zaj_wk = 0.d0
      fsr_wk = 0.d0
      if (kimg==1) then
         do ib1=1,nel
            do ibb2=1,nsize_max_sb_now
               if(ibover(ibb2,iblock)<0) cycle
               ib2 = ibover(ibb2,iblock)
               iter2=(ibb2-1)/nel+1
               ii2=ibb2-nel*(iter2-1)
               ii2=nsta+ii2-1
               hr2=vec(ib2,ib1)

               Do is=1, ndim_spinor
                  do ii=1,iba(ik)
                     zaj_wk(ii,ib1,kimg,is) = zaj_wk(ii,ib1,kimg,is) &
                          &                  + zat_l_noncl(ii,ii2,kimg,iter2,is) *hr2
                  end do
                  do ii=1,nlmta
                     fsr_wk(ib1,ii,is) = fsr_wk(ib1,ii,is) &
                          &             +fsr_noncl(ii2,ii,iter2,is)*hr2
                     fsi_wk(ib1,ii,is) = fsi_wk(ib1,ii,is) &
                          &             +fsi_noncl(ii2,ii,iter2,is)*hr2
                  end do
               End Do

               if ( itot /= max_iter_mdkosugi ) then
                  Do is=1, ndim_spinor
                     do ii=1,iba(ik)
                        zah_wk(ii,ib1,kimg,is) = zah_wk(ii,ib1,kimg,is) &
                             &                 + zah_l_noncl(ii,ii2,kimg,iter2,is) *hr2
                     end do
                  End do
               end if
            end do
         end do
      else
         if (k_symmetry(ik) == GAMMA) then
            do ib1=1,nel
               do ibb2=1,nsize_max_sb_now
                  if(ibover(ibb2,iblock)<0) cycle
                  ib2 = ibover(ibb2,iblock)
                  iter2=(ibb2-1)/nel+1
                  ii2=ibb2-nel*(iter2-1)
                  ii2=nsta+ii2-1
                  hr2=vec(ib2,ib1)

                  Do is=1, ndim_spinor
                     do ii=1,iba(ik)
                        dr1 = zat_l_noncl(ii,ii2,1   ,iter2,is)
                        di1 = zat_l_noncl(ii,ii2,kimg,iter2,is)
                        zaj_wk(ii,ib1,1,   is) = zaj_wk(ii,ib1,1,   is) +dr1*hr2
                        zaj_wk(ii,ib1,kimg,is) = zaj_wk(ii,ib1,kimg,is) +di1*hr2
                     end do
                     do ii=1,nlmta
                        fsr_wk(ib1,ii,is) = fsr_wk(ib1,ii,is) &
                             &             +fsr_noncl(ii2,ii,iter2,is)*hr2
                     end do
                  End do
                  if (itot /= max_iter_mdkosugi) then
                     Do is=1, ndim_spinor
                        do ii=1,iba(ik)
                           dr1 = zah_l_noncl(ii,ii2,1   ,iter2,is)
                           di1 = zah_l_noncl(ii,ii2,kimg,iter2,is)
                           zah_wk(ii,ib1,1,   is) = zah_wk(ii,ib1,1,   is) +dr1*hr2
                           zah_wk(ii,ib1,kimg,is) = zah_wk(ii,ib1,kimg,is) +di1*hr2
                        end do
                     End do
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

                  Do is=1, ndim_spinor
                     do ii=1,iba(ik)
                        dr1 = zat_l_noncl(ii,ii2,1   ,iter2,is)
                        di1 = zat_l_noncl(ii,ii2,kimg,iter2,is)
                        zaj_wk(ii,ib1,1,   is) = zaj_wk(ii,ib1,1,   is) &
                             &                  +dr1*hr2 -di1*hi2
                        zaj_wk(ii,ib1,kimg,is) = zaj_wk(ii,ib1,kimg,is) &
                             &                  +dr1*hi2 +di1*hr2
                     end do
                  End do
                  Do is=1, ndim_spinor
                     do ii=1,nlmta
                        dr1 = fsr_noncl(ii2,ii,iter2,is)
                        di1 = fsi_noncl(ii2,ii,iter2,is)
                        fsr_wk(ib1,ii,is) = fsr_wk(ib1,ii,is) +dr1*hr2 -di1*hi2
                        fsi_wk(ib1,ii,is) = fsi_wk(ib1,ii,is) +dr1*hi2 +di1*hr2
                     end do
                  End do

                  if (itot /= max_iter_mdkosugi) then
                     Do is=1, ndim_spinor
                        do ii=1,iba(ik)
                           dr1 = zah_l_noncl(ii,ii2,1   ,iter2,is)
                           di1 = zah_l_noncl(ii,ii2,kimg,iter2,is)
                           zah_wk(ii,ib1,1,   is) = zah_wk(ii,ib1,1,   is) &
                                &                  +dr1*hr2 -di1*hi2
                           zah_wk(ii,ib1,kimg,is) = zah_wk(ii,ib1,kimg,is) &
                                &                  +dr1*hi2 +di1*hr2
                        end do
                     End Do
                  end if

               end do
            end do
         end if
      end if

      Do is=1, ndim_spinor
         k1 = ik + is -1
         zaj_l(:,nsta:nend,k1,:) = zaj_wk(:,:,:,is)
      End do

      if (itot /= max_iter_mdkosugi) then
         Do is=1, ndim_spinor
            zah_l_noncl(:,nsta:nend,:,1,is) = zah_wk(:,:,:,is)
         End do
      endif

      Do is=1, ndim_spinor
         k1 = ik + is -1
         fsr_l(nsta:nend,:,k1) = fsr_wk(:,:,is)
         if (k_symmetry(ik) /= GAMMA) then
            fsi_l(nsta:nend,:,k1) = fsi_wk(:,:,is)
         endif
      End do

      deallocate(zaj_wk,fsr_wk)
      if(itot /= max_iter_mdkosugi) deallocate(zah_wk)
      if(k_symmetry(ik) /= GAMMA) deallocate(fsi_wk)

    end subroutine subspace_rotation_noncl

  end subroutine evolve_WFs_in_subsp_ksg_noncl
! ================================================================ 11.0

  subroutine decide_correction_vector_ksg(precon,ik,ekin,afft,bfft,p)
    integer, intent(in)       :: precon, ik
    real(kind=DP), intent(in)  :: ekin(kg1)
    real(kind=DP), intent(in)  :: afft(nfft)
    real(kind=DP), intent(out) :: bfft(nfft)
    real(kind=DP)              :: p(kg1)

    integer :: ib
    integer :: id_sname = -1
    call tstatc0_begin('decide_correction_vector ', id_sname,1)

    do ib = ista_e, iend_e, istep_e      ! MPI
       call m_ES_WF_in_Rspace(ik,ib,bfft)! (swffft)
       call m_FFT_Vlocal_W(afft,bfft)    ! (afft, bfft) -> (bfft)
       call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON) ! bfft: G-space repres.
       call SD_direction(precon,ik,ib,ekin,bfft,p) !-here
    end do

!    call orthogonalize_SD_drctns_ksg(ik,to=OTHER_BANDS)  ! -(m_ES_WF_by_SDorCG) ->(wfsd_l, bsd(ri)_l)
    call orthogonalize_SD_drctns_ksg(ik,to=OTHER_BANDS)
!    call orthogonalize_SD_drctns(ik,to=ALL_BANDS)  ! -(m_ES_WF_by_SDorCG) ->(wfsd_l, bsd(ri)_l)

    call normalize_wfsd(ik)

    call tstatc0_end(id_sname)
  end subroutine decide_correction_vector_ksg

! ========================= added by K. Tagami ======================== 11.0
  subroutine decide_correction_vec_ksg_noncl( precon, ik, ekin, &
       &                                   afft_kt, bfft_kt, vnlph_noncl, p )

    integer, intent(in)       :: precon, ik
    real(kind=DP), intent(in)  :: ekin(kg1)
    real(kind=DP), intent(inout)  :: afft_kt(nfft,ndim_chgpot)
    real(kind=DP), intent(out) :: bfft_kt(nfft,ndim_spinor)
    real(kind=DP), intent(in), dimension(kg1,np_e,kimg,ndim_spinor) :: vnlph_noncl

    real(kind=DP)              :: p(kg1)

    integer :: ib, is
    integer :: id_sname = -1

    call tstatc0_begin('decide_correction_vec_ksg_noncl ', id_sname,1)

    do ib = ista_e, iend_e, istep_e      ! MPI

       Do is=1, ndim_spinor
          call m_ES_WF_in_Rspace( ik+is-1, ib, bfft_kt(:,is) ) ! (swffft)
       End do
       call m_FFT_Vlocal_W_noncl( afft_kt, bfft_kt, ndim_chgpot, ndim_spinor )
                                           ! (afft, bfft)-> (bfft)
       Do is=1, ndim_spinor
         call m_FFT_WF( ELECTRON, nfout, bfft_kt(:,is), DIRECT, ON )
       End do

       call SD_direction_noncl( precon,ik, ib, ekin, bfft_kt, &
            &                   vnlph_noncl, p ) !-here
    end do

    call orthogonl_SD_drctns_ksg_noncl(ik,to=OTHER_BANDS)

!    call orthogonalize_SD_drctns(ik,to=ALL_BANDS)  ! -(m_ES_WF_by_SDorCG) ->(wfsd_l, bsd(ri)_l)

    call normalize_wfsd_noncl(ik)

    call tstatc0_end(id_sname)
  end subroutine decide_correction_vec_ksg_noncl
! ================================================================ 11.0

  subroutine orthogonalize_SD_drctns_ksg(ik,to)
    integer, intent(in) :: ik,to

    integer :: itmp
    integer :: id_sname = -1
    call tstatc0_begin('orthogonalize_SD_drctns in Modified Davidson ', id_sname)

!    if(modnrm == EXECUT) call m_ES_betar_dot_Psi_4_each_k(wfsd_l,ik,ik,ik,bsdr_l,bsdi_l)
    !                                        ->bsd(ri)_l

    zat_l(:,:,:,1) = zaj_l(:,:,ik,:)
    fsr(:,:,1)=fsr_l(:,:,ik)
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       fsi(:,:,1)=fsi_l(:,:,ik)
    end if

    zaj_l(:,:,ik,:) = zajold_l(:,:,:)
    fsr_l(:,:,ik) = fsrold_l(:,:)
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       fsi_l(:,:,ik) = fsiold_l(:,:)
    end if

    itmp=modnrm
    modnrm=EXECUT
    call m_ES_betar_dot_Psi_4_each_k(wfsd_l,ik,ik,ik,bsdr_l,bsdi_l)
    call m_ES_orthogonalize_SD_to_WFs(ik,to,wfsd_l,bsdr_l,bsdi_l) ! ->(wfsd_l,bsd(ri)_l)
    modnrm=itmp
!    call m_ES_betar_dot_Psi_4_each_k(wfsd_l,ik,ik,ik,bsdr_l,bsdi_l)
    call tstatc0_end(id_sname)
  end subroutine orthogonalize_SD_drctns_ksg

! ================================= added by K. Tagami ================= 11.0
  subroutine orthogonl_SD_drctns_ksg_noncl(ik,to)
    integer, intent(in) :: ik,to

    integer :: is, itmp, k1
    integer :: id_sname = -1
    call tstatc0_begin('orthogonalize_SD_drctns_noncl in Modified Davidson ', id_sname)

!    if(modnrm == EXECUT) call m_ES_betar_dot_Psi_4_each_k(wfsd_l,ik,ik,ik,bsdr_l,bsdi_l)
    !                                        ->bsd(ri)_l

    Do is=1, ndim_spinor
       k1 = ik + is -1
       zat_l_noncl(:,:,:,1,is) = zaj_l(:,:,k1,:)
       fsr_noncl(:,:,1,is)     = fsr_l(:,:,k1)
    End do

    if (.not.(kv3/ndim_spinor == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       Do is=1, ndim_spinor
          k1 = ik + is -1
          fsi_noncl(:,:,1,is) = fsi_l(:,:,k1)
       End do
    end if

    Do is=1, ndim_spinor
       k1 = ik + is -1
       zaj_l(:,:,k1,:) = zajold_l_noncl(:,:,:,is)
       fsr_l(:,:,k1) = fsrold_l_noncl(:,:,is)
       if(.not.(kv3/ndim_spinor == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
          fsi_l(:,:,k1) = fsiold_l_noncl(:,:,is)
       end if
    End do

    itmp=modnrm
    modnrm=EXECUT
    Do is=1, ndim_spinor
       call m_ES_betar_dot_Psi_4_each_k( wfsd_l, ik, ik +ndim_spinor -1, &
            &                            ik+is-1, bsdr_l, bsdi_l)
    End do

    call m_ES_orthogonl_SD_to_WFs_noncl( ik, to, ik, ik+ndim_spinor-1, &
         &                               wfsd_l, bsdr_l, bsdi_l )
                                                ! ->(wfsd_l,bsd(ri)_l)

    modnrm=itmp

    call tstatc0_end(id_sname)

  end subroutine orthogonl_SD_drctns_ksg_noncl
! =================================================================== 11.0


#endif
end module m_ES_WF_by_ModifiedDavidson

