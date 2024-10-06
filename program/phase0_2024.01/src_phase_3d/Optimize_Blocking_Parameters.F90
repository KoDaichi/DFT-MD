
subroutine Optimize_Blocking_Parameters()
  use m_Const_Parameters,     only : DP, ORTHOGONALIZATION, EXECUT, VDB, NORMCONSERVATION, ON, OFF
  use m_Files,                only : nfout
  use m_PseudoPotential,      only : modnrm
  use m_Control_Parameters ,  only : printable, kimg, ipriblsize &
  &                                , nblsizecand_betar, blsizecand_betar, nblsizecand_mgs, blsizecand_mgs &
  &                                , nblocksize_mgs, nblocksize_mgs_is_given, nblocksize_betar_dot_wfs    &
  &                                , nblocksize_betar_is_given, ipriparallel,neg,nb_mgs_default &
  &                                , nblsizecand_vnonlocal_w, blsizecand_vnonolocal_w &
  &                                , nblocksize_vnonlocal_is_given, nblocksize_vnonlocal_w &
  &                                , nblsizecand_submat, blsizecand_submat &
  &                                , nblocksize_submat_is_given, nblocksize_submat &
  &                                , m_CntrlP_rst_submat_call_stat,nspin &
  &                                , sw_betar_dot_wfs_exp, sw_precalculate_phase_vnonlocal, damp, meg &
  &                                , sw_scalapack
  use m_Kpoints,              only : kv3
  use m_Parallelization,      only : np_g1k, np_e, np_fs
#ifdef FFT_3D_DIVISION
  use m_Parallelization,      only : fft_X_x_nel, fft_X_y_nel, fft_X_z_nel
#else
  use m_Parallelization,      only  :nel_fft_x, nel_fft_y, nel_fft_z
#endif
  use m_Parallelization,      only : mype, myrank_e, m_Parallel_dealloc_mpi_elec, m_Parallel_init_mpi_elec_3D &
  &                                , make_ball_buff, MPI_CommGroup, map_k, myrank_k
  use m_ES_ortho,             only : mgs_4_each_k_G_3D, m_ESortho_mgs_alloc, m_ESortho_mgs_dealloc &
  &                                , m_ES_modified_gram_schmidt
  use m_ES_nonlocal,          only : m_ES_betar_dot_Psi_4_each_k_3D, m_ES_Vnonlocal_W_3D
  use m_Parallelization,      only : make_index_band_3D,make_index_band_for_Gdiv_3D, nrank_e, nrank_k
  use m_ES_WF_by_submat,      only : evolve_WFs_in_subspace_3D, m_ESsubmat_alloc, m_ESsubmat_dealloc, m_ESsubmat_renew_WF
  use m_Timing,               only : tstatc0_begin, tstatc0_end
  use m_PlaneWaveBasisSet,    only : kg1, iba
  use m_Electronic_Structure, only : m_ES_dealloc, m_ES_alloc_zaj_etc, m_ES_alloc_vlhxc, m_ES_alloc_vlhxcQ
  use mpi


  implicit none

!  include 'mpif.h'             ! MPI

  integer :: i, t1, nbuf, lsize
  real(kind=DP) :: elapsed_time
  real(kind=DP), allocatable, dimension(:,:,:,:) :: phi_l
  real(kind=DP), allocatable, dimension(:,:,:) :: bsdr_l, bsdi_l
  real(kind=DP), allocatable, dimension(:,:) :: afft_l
  real(kind=DP), allocatable, dimension(:) :: ekin_l
  real(kind=DP) :: best
  logical :: logi
  integer :: max_blsize, ierr
  integer :: id_sname = -1

  call tstatc0_begin('Optimize_Blocking_Parameters ',id_sname,-1)
  if (nrank_k>1 .and. sw_scalapack == on ) then
    call phase_error_with_msg(nfout, &
   & 'optimization of the blocking parameters, kpoint parallelization and ScaLAPACK cannot be used simultaneously', &
   & __LINE__, __FILE__)
  endif
  nblocksize_mgs_is_given = .true.
  best = HUGE(0)
  nbuf = blsizecand_mgs(1)
  max_blsize = int(real(neg)/real(nrank_e))
  blsizecand_mgs(1) = max_blsize
  do i=1,nblsizecand_mgs
    blsizecand_mgs(i) = max_blsize/(2**(i-1))
    if (blsizecand_mgs(i) < 4) blsizecand_mgs(i) = 4
    if (blsizecand_betar(i)>max_blsize) cycle
    call system_clock(t1)
    nblocksize_mgs = blsizecand_mgs(i)
    call m_CntrlP_rst_submat_call_stat()
    call make_index_band_3D(nfout,ipriparallel,printable,kv3,neg,nblocksize_mgs,nblocksize_mgs_is_given,nb_mgs_default)
    call make_index_band_for_Gdiv_3D(neg, nblocksize_mgs,nblocksize_mgs_is_given,nb_mgs_default)
    allocate(phi_l(maxval(np_g1k),np_e,1,kimg))
    if(modnrm == EXECUT) then
      allocate(bsdr_l(np_e, np_fs, 1))
      allocate(bsdi_l(np_e, np_fs, 1))
    endif

    if(modnrm == EXECUT) then
      call mgs_4_each_k_G_3D(1,1,1,phi_l,ORTHOGONALIZATION,bsdr_l,bsdi_l,mod_pot=VDB,dryrun=.true.)
    else
      call mgs_4_each_k_G_3D(1,1,1,phi_l,ORTHOGONALIZATION,mod_pot=NORMCONSERVATION,dryrun=.true.)
    endif
    elapsed_time = get_time_from(t1)
    if(elapsed_time<best) then
      best = elapsed_time
      nbuf = nblocksize_mgs
    endif
    deallocate(phi_l)
    if(modnrm == EXECUT) then
      deallocate(bsdr_l)
      deallocate(bsdi_l)
    endif
    if(ipriblsize>=2) write(nfout,'(a,i5,1f12.4)') '!** MGS blocksize and elapsed time ',blsizecand_mgs(i),elapsed_time
  enddo
  call mpi_bcast(nbuf,1,mpi_integer,0,MPI_CommGroup,ierr)
  nblocksize_mgs = nbuf
  call m_CntrlP_rst_submat_call_stat()
  call make_index_band_3D(nfout,ipriparallel,printable,kv3,neg,nblocksize_mgs,nblocksize_mgs_is_given,nb_mgs_default)
  call make_index_band_for_Gdiv_3D(neg, nblocksize_mgs,nblocksize_mgs_is_given,nb_mgs_default)
  call make_ball_buff()
  call m_ES_dealloc()
  call m_ES_alloc_zaj_etc()
  call m_ES_alloc_vlhxc()
  call m_ES_alloc_vlhxcQ()

  allocate(phi_l(maxval(np_g1k),np_e,1,kimg))
  allocate(bsdr_l(np_e, np_fs, 1))
  allocate(bsdi_l(np_e, np_fs, 1))
  if(sw_betar_dot_wfs_exp==OFF) then
  best = HUGE(0)
  nbuf = blsizecand_betar(1)
  nblocksize_betar_is_given = .true.
  do i=1,nblsizecand_betar
    call system_clock(t1)
    nblocksize_betar_dot_wfs = blsizecand_betar(i)
    call m_ES_betar_dot_Psi_4_each_k_3D(nfout,phi_l,1,1,1,bsdr_l,bsdi_l,mod_ball=OFF)
    elapsed_time = get_time_from(t1)
    if(elapsed_time<best) then
      best = elapsed_time
      nbuf = nblocksize_betar_dot_wfs
    endif
    if(ipriblsize>=2) write(nfout,'(a,i5,1f12.4)') '!** betar_dot_WFs blocksize and elapsed time ' &
    & ,blsizecand_betar(i),elapsed_time
  enddo
  call mpi_bcast(nbuf,1,mpi_integer,0,MPI_CommGroup,ierr)
  nblocksize_betar_dot_wfs = nbuf
  endif

  if(sw_precalculate_phase_vnonlocal==OFF) then
  best = HUGE(0)
  nbuf = blsizecand_vnonolocal_w(1)
  nblocksize_vnonlocal_is_given = .true.
  do i=1,nblsizecand_vnonlocal_w
    call system_clock(t1)
    nblocksize_vnonlocal_w = blsizecand_vnonolocal_w(i)
    call m_ES_Vnonlocal_W_3D(1,1,1,switch_of_eko_part=OFF)
    elapsed_time = get_time_from(t1)
    if(elapsed_time<best) then
      best = elapsed_time
      nbuf = nblocksize_vnonlocal_w
    endif
    if(ipriblsize>=2) write(nfout,'(a,i5,1f12.4)') '!** vnonlocal_w blocksize and elapsed time ' &
    & ,blsizecand_vnonolocal_w(i),elapsed_time
  enddo
  call mpi_bcast(nbuf,1,mpi_integer,0,MPI_CommGroup,ierr)
  nblocksize_vnonlocal_w = nbuf
  endif

  if(map_k(1) == myrank_k) then
    best = HUGE(0)
    nbuf = blsizecand_submat(1)
    nblocksize_submat_is_given = .true.
#ifdef FFT_3D_DIVISION
    lsize = fft_X_x_nel*fft_X_y_nel*fft_X_z_nel
    allocate(afft_l(lsize*2,1))
#else
    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
    allocate(afft_l(lsize*kimg,1))
#endif
    allocate(ekin_l(np_g1k(1)))
    call m_ESsubmat_alloc()
    do i=1,nblsizecand_submat
      blsizecand_submat(i) = neg/i
      if(blsizecand_submat(i)==0) blsizecand_submat(i) = 1
      if (blsizecand_submat(i)>neg) exit
      call system_clock(t1)
      nblocksize_submat = blsizecand_submat(i)
      call m_CntrlP_rst_submat_call_stat()
      call evolve_WFs_in_subspace_3D(1,1,1,meg,damp,ekin_l,afft_l(1,1),lsize,dryrun=.true.)
      elapsed_time = get_time_from(t1)
      if(elapsed_time<best) then
        best = elapsed_time
        nbuf = nblocksize_submat
      endif
      if(ipriblsize>=2) write(nfout,'(a,i5,1f12.4)') '!** submat blocksize and elapsed time ' &
      & ,blsizecand_submat(i),elapsed_time
    enddo
    call m_CntrlP_rst_submat_call_stat()
    call m_ESsubmat_dealloc()
    deallocate(afft_l)
    deallocate(ekin_l)
  endif
  call mpi_bcast(nbuf,1,mpi_integer,0,MPI_CommGroup,ierr)
  nblocksize_submat = nbuf

  if(printable) then
    write(nfout,'(a)')    '!** estimated optimal blocking parameters'
    write(nfout,'(a,i5)') '    nblocksize_mgs           = ',nblocksize_mgs
    if(sw_betar_dot_wfs_exp==OFF) write(nfout,'(a,i5)') '    nblocksize_betar_dot_wfs = ',nblocksize_betar_dot_wfs
    if(sw_precalculate_phase_vnonlocal==OFF)write(nfout,'(a,i5)') '    nblocksize_vnonlocal_w   = ',nblocksize_vnonlocal_w
    write(nfout,'(a,i5)') '    nblocksize_submat        = ',nblocksize_submat
  endif

  deallocate(phi_l)
  deallocate(bsdr_l)
  deallocate(bsdi_l)

  call tstatc0_end(id_sname)

  contains

  real(kind=DP) function get_time_from(t1)
    integer, intent(in) :: t1
    integer :: t2, t_rate, t_max, diff
    call system_clock(t2, t_rate, t_max)   ! 終了時を記録
    if ( t2 < t1 ) then
      diff = (t_max - t1) + t2 + 1
    else
      diff = t2 - t1
    endif
    get_time_from = diff/dble(t_rate)
  end function get_time_from

end subroutine Optimize_Blocking_Parameters
