#ifdef HIUX
*option MP(P(0))
#endif
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 635 $)
!
!  MODULE: m_Parallelization
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!  FURTHER MODIFICATION: T. Yamasaki, August/23/2006
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
! ***************************************************************
!
! This module "m_Parallelization" is for mpi-parallelization and was
! coded by T. Yamasaki (JRCAT-ATP, Fujitsu Laboratories Ltd.) in 1999.
!
#ifdef __TIMER__
#   define __TIMER_START(a) call timer_init() ; call timer_sta(a)
#   define __TIMER_STOP     call timer_fin()
#else
#   define __TIMER_START(a)
#   define __TIMER_STOP
#endif
#ifdef __TIMER_SUB__
#   define __TIMER_SUB_START(a)  call timer_sta(a)
#   define __TIMER_SUB_STOP(a)   call timer_end(a)
#else
#   define __TIMER_SUB_START(a)
#   define __TIMER_SUB_STOP(a)
#endif
#ifdef __TIMER_INICOMM__
#   define __TIMER_INICOMM_START_w_BARRIER(str,a)   call timer_barrier(str) ; call timer_sta(a)
#   define __TIMER_INICOMM_STOP(a)         call timer_end(a)
#else
#   define __TIMER_INICOMM_START_w_BARRIER(str,a)
#   define __TIMER_INICOMM_STOP(a)
#endif

module m_Parallelization
!   (m_Parallel)
! $Id: m_Parallelization.F90 635 2021-02-26 07:16:10Z jkoga $
  use m_Const_Parameters, only       : ON, OFF, tag_npes_etc
  use m_ErrorMessages
  use mpi
#ifdef _CUDA_
  use cudafor
#endif

  implicit none
  integer :: npes, mype, nrank_e, nrank_k, nrank_g, nrank_s

  integer :: nrank_g1
  integer, allocatable, dimension(:) :: map_e, map_z              ! d(neg)  --> expl.1 and 5
  integer                            :: myrank_e,ista_e,iend_e,istep_e,np_e,mp_e ! --> expl.2 to 4

#ifdef NEC_TUNE_SOFT
  integer, allocatable, dimension(:) :: ista_e_smp, iend_e_smp
  integer                            :: itask
#elif NEC_TUNE_FFT
  integer                            :: itask
#endif
  integer, allocatable, dimension(:) :: nis_e, nie_e, nel_e, idisp_e   ! d(nrank_e)
  integer, allocatable, dimension(:) :: map_k                     ! d(kv3)
  integer, allocatable, dimension(:,:):: map_ek                   ! d(neg,kv3)
  integer, allocatable, dimension(:,:,:):: map_eks                  ! d(neg,kv3/nspin,nspin)
  integer, allocatable, dimension(:) :: map_s                   ! d(neg,kv3)
  integer                            :: myrank_k,ista_k,iend_k
  integer, allocatable, dimension(:,:,:) :: map_rank_gek  ! d(nrank_g,nrank_e,nrank_k)
  integer, allocatable, dimension(:) :: nis_k,  nie_k,  nel_k


!!$#ifdef TRANSPOSE
  integer                            :: ista_g1, iend_g1, np_g1, mp_g1
  integer, allocatable, dimension(:) :: nis_g1, nie_g1, nel_g1      !d(0:nrank_g1-1)
  integer, allocatable, dimension(:) :: ista_g1k, iend_g1k, np_g1k, mp_g1k, np_g1k_prev !d(kv3)
  integer, allocatable, dimension(:) :: ista_g1k_prev, iend_g1k_prev
  integer, allocatable, dimension(:,:):: nis_g1k, nie_g1k, nel_g1k  !d(0:nrank_g1-1,kv3)

  integer                            :: ista_fs, iend_fs, np_fs, mp_fs
  integer, allocatable, dimension(:) :: nis_fs, nie_fs, nel_fs      !d(0:nrank_g-1)
  integer                            :: ista_fs_atm, iend_fs_atm, np_fs_atm, mp_fs_atm
  integer, allocatable, dimension(:) :: nis_fs_atm, nie_fs_atm, nel_fs_atm !d(0:nrank_e-1)
  integer                            :: myrank_g
  integer                            :: ista_lmta,iend_lmta
  integer, allocatable, dimension(:) :: lmta_atm

!!$#endif
  integer, allocatable, dimension(:) :: mpi_k_world  ! kd(0:nrank_k-1)
  integer, allocatable, dimension(:) :: mpi_e_world  ! kd(0:nrank_k-1)
! === DEBUG by tkato 2012/04/04 ================================================
  integer                            :: mpi_ge_world
! ==============================================================================
  integer                            :: ista_kngp, iend_kngp, np_kngp, mp_kngp
  integer                            :: ista_kngp_prev, iend_kngp_prev, np_kngp_prev
  integer, allocatable, dimension(:) :: is_kngp, ie_kngp, nel_kngp

  integer                            :: ista_kngp_gw, iend_kngp_gw, np_kngp_gw, mp_kngp_gw
  integer, allocatable, dimension(:) :: is_kngp_gw, ie_kngp_gw, nel_kngp_gw
  integer                            :: ista_kngp_exx, iend_kngp_exx, np_kngp_exx, mp_kngp_exx

  integer, allocatable, dimension(:) :: is_kngp_exx, ie_kngp_exx, nel_kngp_exx
! natm for fxyzew_l
  integer                            :: ista_atm, iend_atm, np_atm, mp_atm
  integer, allocatable, dimension(:) :: is_atm, ie_atm, nel_atm

  integer                            :: ista_atm_ke, iend_atm_ke, np_atm_ke, mp_atm_ke
  integer, allocatable, dimension(:) :: is_atm_ke, ie_atm_ke, nel_atm_ke,map_atm_ke

  integer                            :: ista_atm_f, iend_atm_f, np_atm_f, mp_atm_f
  integer, allocatable, dimension(:) :: is_atm_f, ie_atm_f, nel_atm_f

  integer                            :: ista_rspace_aug, iend_rspace_aug, np_rspace_aug, mp_rspace_aug
  integer, allocatable, dimension(:) :: ista_rspace_aug_atm, iend_rspace_aug_atm
  integer, allocatable, dimension(:) :: is_rspace_aug, ie_rspace_aug, nel_rspace_aug

  integer                            :: ista_nn,iend_nn,np_nn,mp_nn
  integer, allocatable, dimension(:) :: is_nn, ie_nn, nel_nn

  integer                            :: ista_nq, iend_nq, np_nq, mp_nq
  integer, allocatable, dimension(:) :: is_nq, ie_nq, nel_nq, map_nq, map_z_nq

!  integer                            :: ista_nval,iend_nval,np_nval,mp_nval
!  integer, allocatable, dimension(:) :: is_nval, ie_nval, nel_nval

! natm2 for s_ew in m_Ionic_System
  integer                            :: ista_atm2, iend_atm2, np_atm2, mp_atm2
  integer, allocatable, dimension(:) :: is_atm2, ie_atm2, nel_atm2

! BROYDEN or DFP MIXING METHOD
  integer                            :: ista_kgpm,iend_kgpm, np_kgpm, mp_kgpm
  integer, allocatable, dimension(:) :: is_kgpm,  ie_kgpm, nel_kgpm
  integer                            :: ista_urec_hsr, iend_urec_hsr
  logical, save ::                      ista_and_iend_urec_hsr_set = .false.

! SPIN
  integer                            :: ista_spin, iend_spin, np_spin, mp_spin
  integer, allocatable, dimension(:) :: nis_spin, nie_spin, nel_spin
  integer                            :: myrank_spin
  integer                            :: mpi_spin_group

! FFT BOX
  integer                            :: npes_cdfft, nrank_ggacmp, nrest_cdfft
!!$  integer                            :: myrank_cdfft, myrank_ggacmp, max_ggacmp
  integer                            :: myrank_cdfft, myrank_ggacmp
  integer, allocatable, dimension(:) :: map_ggacmp ! d(3)
  integer, allocatable, dimension(:) :: mpi_ggacmp_cross_world ! d(0:npes_cdfft-1)
  integer, allocatable, dimension(:) :: mpi_cdfft_world        ! d(0:nrank_ggacmp-1)
  integer, allocatable, dimension(:) :: map_pe2ggacmp, map_pe2cdfft    ! d(0:npes-1)

  integer                            :: ista_sfftp, iend_sfftp, np_sfftp, mp_sfftp
  integer                            :: ista_sfftph,iend_sfftph
  integer, allocatable, dimension(:) :: nis_sfftp,  nie_sfftp,  nel_sfftp, idisp_sfftp
  integer, allocatable, dimension(:) :: nis_sfftph, nie_sfftph, nel_sfftph
  integer                            :: ista_fftp, iend_fftp, np_fftp, mp_fftp
  integer                            :: ista_fftph,iend_fftph
  integer, allocatable, dimension(:) :: nis_fftp,  nie_fftp,  nel_fftp,  idisp_fftp
  integer, allocatable, dimension(:) :: nis_fftph, nie_fftph, nel_fftph, idisp_fftph


! NonLocalPotential snl
  integer                            :: ista_snl, iend_snl

! nbmx
  integer  :: ng_nbmx, myrank_nbmx, nrank_nbmx, np_nbmx, mp_nbmx,ista_nbmx, iend_nbmx
  integer  :: nbmx_ext
  integer, allocatable, dimension(:) :: nis_nbmx,nie_nbmx,nel_nbmx, idisp_nbmx
  integer, allocatable, dimension(:) :: mpi_nbmx_world  ! kd(0:nrank_nbmx-1)
  integer  :: ng_nbmx_k, myrank_nbmx_k, nrank_nbmx_k, np_nbmx_k, mp_nbmx_k &
       & ,ista_nbmx_k, iend_nbmx_k
  integer, allocatable, dimension(:) :: nis_nbmx_k, nie_nbmx_k, nel_nbmx_k, idisp_nbmx_k
! kg1
  integer, allocatable, dimension(:) :: mpi_nbmx_world_k  ! kd(0:nrank_nbmx_k-1)
  integer  :: np_kg1_k, mp_kg1_k, ista_kg1_k, iend_kg1_k
  integer  :: kg1_ext
  integer, allocatable, dimension(:) :: nis_kg1_k, nie_kg1_k, nel_kg1_k, idisp_kg1_k

  integer                            :: ierr
  character(len=50)                  :: workdir
  integer                            :: sw_wdir = 0

  integer                            :: MPI_CommGroup
  integer                            :: color,key

  logical,save                       :: conf_para = .false.
  integer,save                       :: mype_conf = 0
  integer                            :: nrank_conf

  integer :: nrank_nvale,myrank_nvale
  integer :: ista_nvale,iend_nvale,np_nvale,mp_nvale
  integer, allocatable, dimension(:) :: nis_nvale,nie_nvale,nel_nvale
  integer, allocatable, dimension(:) :: map_nvale,map_z_nvale
  logical,private :: mpi_nvale_enabled = .false.

  integer :: nrank_kv3_ek,myrank_kv3_ek
  integer :: ista_kv3_ek,iend_kv3_ek,np_kv3_ek,mp_kv3_ek
  integer, allocatable, dimension(:) :: nis_kv3_ek,nie_kv3_ek,nel_kv3_ek
  integer, allocatable, dimension(:) :: map_kv3_ek,map_z_kv3_ek


! for nfft
  integer                            :: ista_ffth, iend_ffth, np_ffth, mp_ffth
  integer, allocatable, dimension(:) :: is_ffth, ie_ffth, nel_ffth

! whether parallelization was configured from the command line or not
  logical :: read_from_args = .false.

! 1. map_e(neg)
! rank# = map_e(ieg) : rank# is the rank number which operates ieg-th data
!  -- examples --
!  (case 1) cyclic partitioning case (assuming nrank_e == 3 and neg == 8)
!   ieg :   1   2   3   4   5   6   7   8
! rank# :   0   1   2   0   1   2   0   1
!
!  (case 2) block  partitioning case (also assuming nrank_e == 3 and neg == 8)
!   ieg :   1   2   3   4   5   6   7   8
! rank# :   0   0   0   1   1   1   2   2
!
! 2. myrank_e
!  This lets each node know its rank number according to an axis-e.
!  -- an example --
!  Assuming nrank_e == 4 and nrank_k == 2
!   #rank (or #pe)  0  1  2  3  4  5  6  7
!   myrank_e        0  1  2  3  0  1  2  3
!
! 3. ista_e, iend_e, istep_e, nis_e, nie_e
!  ista_e, iend_e : Bounds of a window which is a range 'e' runs in each nrank_e.
!  nis_e,  nie_e  : Bounds of windows. This shows all windows ranges.
!  istep_e        : = nrank_e (the cyclic partitioning case)
!                   = 1       (the block partitioning case)
!  -- examples --
!  Assuming nrank_e == 3, nrank_k == 2 and neg ==  8
!  (case 1) cyclic partitioning case
!   #rank (or #pe)  0  1  2  3  4  5
!   ista_e          1  2  3  1  2  3
!   iend_e          7  8  6  7  8  6
!   istep_e         3  3  3  3  3  3
!   nis_e(#rank)    1  2  3  1  2  3
!   nie_e(#rank)    7  8  6  7  8  6
!
!  (case 2) block partitioning case
!   #rank (or #pe)  0  1  2  3  4  5
!   ista_e          1  4  7  1  4  7
!   iend_e          3  6  8  3  6  8
!   istep_e         1  1  1  1  1  1
!   nis_e(#rank)    1  4  7  1  4  7
!   nie_e(#rank)    3  6  8  3  6  8
!
! 4. np_e, mp_e, nel_e
!  np_e is the number of elements of e which each rank_e operates
!  mp_e is the maximum of np_e's over all nodes
!  -- examples --
!  Assuming nrank_e == 3, nrank_k == 2 and neg ==  8
!   #rank (or #pe)  0  1  2  3  4  5
!   np_e            3  3  2  3  3  2
!   mp_e            3  3  3  3  3  3
!   nel_e(#rank)    3  3  2  3  3  2
!
! 5. map_z(neg)
! j = map_z(ieg) : ieg-th data is mapped on j-th position on the
!                 map_e(i)-th array. Common on all ranks
!  -- examples --
! Assuming nrank_e = 3 and neg = 8
!  (case 1) cyclic partitioning case
!   ieg    :   1   2   3   4   5   6   7   8
!  rank_e# :   0   1   2   0   1   2   0   1
!    j     :   1   1   1   2   2   2   3   3
!  (case 2) block partitioning case
!   ieg    :   1   2   3   4   5   6   7   8
!  rank_e# :   0   0   0   1   1   1   2   2
!    j     :   1   2   3   1   2   3   1   2
!
! 6. map_k(kv3)
! rank_k# = map_k(ik) : rank_k# is the rank number which operates ik-th data
!  -- an examples --
!  Assuming nrank_k == 2 and kv3 == 9)
!    ik   :   1   2   3   4   5   6   7   8   9
! rank_k# :   0   0   0   0   0   1   1   1   1
!
! 7. map_ek(neg,kv3)
! rank# = map_ek(ie,ik) : rank# is,  in the whole using nodes, the rank number which operates
!                         ie-th and -k-th data
!  -- an example --
!  Assuming nrank_e == 2, nrank_k = 2, npes == 4, neg == 8 and kv3 = 5
!  and taking block partitioning for an axis-'e'.
!   ik   : 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5
!   ie   : 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8
! rank_k#: 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
! rank_e#: 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1
! rank#  : 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3 2 2 2 2 3 3 3 3
!
! 8. myrank_k
!  This lets each node know its rank number according to an axis-k.
!  -- an example --
!  Assuming nrank_e == 4 and nrank_k == 2
!   #rank (or #pe)  0  1  2  3  4  5  6  7
!   myrank_k        0  0  0  0  1  1  1  1
!   myrank_e        0  1  2  3  0  1  2  3  (for reference, assuming block partitioning)
!
! 9. ista_k, iend_k
!  ista_k, iend_k : bounds of a window which is a range each nrank_k node operates
!
!10. ista_g1, iend_g1, nis_g1, nie_g1, nel_g1
!  ista_g1,iend_g1: bounds of a window which is a range 'g' runs in each nrank_e after transpose.
!  -- An example --
!  Assuming nrank_e == 3, nrank_k == 2 and kg1 == 101
!   #rank (or #pe)  0   1   2   3   4   5
!   ista_e          1  35  69   1  35  69
!   iend_e         34  68 101  34  68 101
!
!11. np_g1, mp_g1
!  np_g1 is the number of elements of 'g' which each rank_e operates
!  mp_g1 is the maximum of np_g1's over all nodes
!  -- examples --
!  Assuming nrank_e == 3, nrank_k == 2 and kg1 ==  101
!   #rank (or #pe)  0  1  2  3  4  5
!   np_g1          34 34 33 34 34 33
!   mp_g1          34 34 34 34 34 34
!
!12. mpi_k_world(0:nrank_k-1)
! Sub-group of rank_k
!
!13. ng_nbmx, myrank_nbmx, nrank_nbmx,np_nbmx, mp_nbmx, nel_nbmx(#rank_nbmx),
!    ista_nbmx, and iend_nbmx
!  np_nbmx is the number of elements of G which each rank_nbmx operates
!  mp_nbmx is the maximum of np_nbmx's over all nodes in a mpi_nbmx_world
!  -- example --
!  Assuming nbmx = 152, nrank_nbmx = 3 or 2, npe = 5, and ngnode_nbmx = 2
!  #rank (or #pe)               0   1   2   3   4
!  ng_nbmx                      0   1   0   1   0
!  myrank_nbmx                  0   0   1   1   2
!  nrank_nbmx                   3   2   3   2   3
!  np_nbmx                     51  76  51  76  50
!  mp_nbmx                     51  76  51  76  51
!  nel_nbmx(#rank_nbmx)        51  76  51  76  50
!  ista_nbmx                    1   1  52  77 103
!  iend_nbmx                   51  76 102 152 152
!
!                          ||
!
!  ng_nbmx            (      0       ) (    1    )
!  #rank ( or #node)    0    2    4       1   3
!  nrank_nbmx           3    3    3       2   2
!  myrank_nbmx          0    1    2       0   1
!    np_nbmx           51   51   50      76  76
!    mp_nbmx           51   51   51      76  76
! nel_nbmx(#rank_nbmx) 51   51   50      76  76
!    ista_nbmx          1   52  103       1  77
!    iend_nbmx         51  102  152      76 152
!
!14. mpi_nbmx_world(0:nrank_nbmx-1)
! Sub-group of rank_nbmx
!
!  include 'mpif.h'             ! MPI
  integer                            :: istatus(MPI_STATUS_SIZE)
contains

  subroutine m_Parallel_init_comm_world(init_mpi)
!!$    character :: F_OUT*13, name_mype*3
!!$    character :: name_jobstep*3
!!$    integer ::   i, js
!!$    logical ::   existence
    integer, optional, intent(in) :: init_mpi
    integer :: inimpi
    integer :: i,lenc,lencom,ntmp,numdir,np,hispe
    logical :: exist_dirlist
    character(len=50) :: buf
    character(len=50), allocatable :: workdir_array(:)

    integer, allocatable :: color(:),key(:),nproc(:)
    integer :: my_color, my_key
    integer :: maxproc, count, idir, npes1, numdir1
    logical, save :: firstcall = .true.
#ifdef _CUDA_
    integer :: nDevices, istat, iDevice, mype_loc, gpuID, shmcomm
#endif

    workdir = ''
    inimpi=ON
    if(present(init_mpi)) inimpi = init_mpi

    if(inimpi==ON) then

    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world, npes, ierr)
    call mpi_comm_rank(mpi_comm_world, mype, ierr)
                                    __TIMER_START(17)
#ifdef FJ_TIMER
                 call timer_sta(21)
#endif

!! Group creation
    inquire(file='dirlist',exist=exist_dirlist)
    if(.not.exist_dirlist) then

       !if(firstcall) call mpi_comm_dup(mpi_comm_world,MPI_CommGroup,ierr)
      if(.not.firstcall) call mpi_comm_free(MPI_CommGroup, ierr)
      call mpi_comm_dup(mpi_comm_world,MPI_CommGroup,ierr)

       !!write(*,*) 'mpi_comm_world=',mpi_comm_world
       !!write(*,*) 'MPI_CommGroup=',MPI_CommGroup

    else
       !!write(*,*) 'dirlist exists.'
       !! switching working directory
       sw_wdir = 1

       if(mype == 0) then
          open(unit=100,file="dirlist",form='formatted')
          read(100,*) numdir
          !! debug
          !!write(*,*) 'numdir=',numdir
          allocate(workdir_array(numdir))
          ntmp=0
          do i=1,numdir
             buf(1:50) = ' '
             read(100,*,err=500,end=500) buf
             ntmp=ntmp+1
             lenc = len_trim(buf)
             buf(lenc+1:lenc+1) = '/'
             workdir_array(i) = buf
          end do
500       continue
          numdir = ntmp

          allocate(color(npes))
          allocate(key(npes))
          allocate(nproc(npes))

          npes1 = npes
          numdir1 = numdir
          count = 0
          my_color = 0
          do while(npes1 > 0)
             maxproc = npes1/numdir1
             if(npes1-maxproc*numdir1>0) maxproc=maxproc+1
             my_color = my_color + 1
             do i=1,maxproc
                nproc(i+count) = maxproc
                color(i+count) = my_color
                key(i+count) = i-1
             end do
             count = count + maxproc
             npes1 = npes1 - maxproc
             numdir1 = numdir1 -1
          end do

          !!do i=1,npes
          !!   key(i) = mod(i-1,maxproc)
          !!   color(i) = (i-key(i))/maxproc
          !!end do
          !!write(*,'("nproc:",50(1x,i2))') nproc
          !!write(*,'("color:",50(1x,i2))') color
          !!write(*,'("key:",50(1x,i2))') key

          do i=1,npes
             idir = color(i)
             buf = workdir_array(idir)
             !!$write(*,*) buf(1:lenc)
             hispe = i-1
             if(mype /= hispe) then
                call MPI_Send(workdir_array(idir),50,MPI_CHARACTER &
                           & ,hispe,1,MPI_COMM_WORLD,ierr)
                call MPI_Send(color(i),1,MPI_INTEGER &
                           & ,hispe,1,MPI_COMM_WORLD,ierr)
                call MPI_Send(key(i),1,MPI_INTEGER &
                           & ,hispe,1,MPI_COMM_WORLD,ierr)
             else
                workdir = workdir_array(idir)
                my_color = color(i)
                my_key = key(i)
             end if
          end do

          deallocate(workdir_array)
       else
          call MPI_Recv(workdir,50,MPI_CHARACTER, &
                   & MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,istatus,ierr)
          call MPI_Recv(my_color,1,MPI_INTEGER, &
                   & MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,istatus,ierr)
          call MPI_Recv(my_key,1,MPI_INTEGER, &
                   & MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,istatus,ierr)
       end if
       call MPI_barrier(MPI_COMM_WORLD,ierr)
       !!write(*,*) 'workdir=',trim(workdir)
       !!my_color = mype
       !!my_key = 0
       !if(firstcall) call mpi_comm_split(MPI_COMM_WORLD,my_color,my_key,MPI_CommGroup,ierr)
       if(.not.firstcall) call mpi_comm_free(MPI_CommGroup, ierr)
       call mpi_comm_split(MPI_COMM_WORLD,my_color,my_key,MPI_CommGroup,ierr)
       call mpi_comm_size(MPI_CommGroup, npes, ierr)
       call mpi_comm_rank(MPI_CommGroup, mype, ierr)
       !!write(*,*) 'MPI_CommGroup=',MPI_CommGroup
       !!write(*,*) 'npes=',npes
       !!write(*,*) 'myes=',mype
       !!write(*,*) 'workdir=',trim(workdir)

       !! debug
       !!call mpi_Finalize(ierr)
       !!stop 'group_creation end'

    end if
!! Group creation

    else

    call mpi_comm_size(MPI_CommGroup, npes, ierr)
    call mpi_comm_rank(MPI_CommGroup, mype, ierr)

    endif
#ifdef _CUDA_
    call mpi_comm_split_type(mpi_comm_world, mpi_comm_type_shared, 0, mpi_info_null, shmcomm, ierr)
    call mpi_comm_rank(shmcomm, mype_loc, ierr)
    istat = cudaGetDeviceCount(nDevices)
    gpuID = mod(mype_loc, nDevices)
    istat = cudaSetDevice(gpuID)
#endif

    firstcall = .false.

!!$    if(npes > 1) then
!!$       write(name_mype,'(i3)') mype
!!$       do i = 1, 3
!!$          if(name_mype(i:i) == ' ') name_mype(i:i)='0'
!!$       enddo
!!$
!!$       if(mype==0) then
!!$          do js = 0, 999
!!$             write(name_jobstep,'(i3)') js
!!$             do i = 1, 3
!!$                if(name_jobstep(i:i) == ' ') name_jobstep(i:i)='0'
!!$             enddo
!!$             F_OUT = 'output'//name_jobstep
!!$             inquire(file=F_OUT, exist=existence)
!!$             if(.not.existence) goto 1001
!!$          end do
!!$1001      if(js >=1000) js = 999
!!$       end if
!!$       if(npes > 1) call mpi_bcast(js,1,mpi_integer,0,mpi_comm_world,ierr)
!!$       write(name_jobstep,'(i3)') js
!!$       do i = 1, 3
!!$          if(name_jobstep(i:i) == ' ') name_jobstep(i:i)='0'
!!$       enddo
!!$
!!$       if(mype==0) then
!!$          F_OUT = 'output'//name_jobstep
!!$       else
!!$          F_OUT = 'output'//name_jobstep//'_'//name_mype
!!$       end if
!!$    else if(npes == 1) then
!!$       do js = 0, 999
!!$          write(name_jobstep,'(i3)') js
!!$          do i = 1, 3
!!$             if(name_jobstep(i:i) == ' ') name_jobstep(i:i)='0'
!!$          enddo
!!$          F_OUT = 'output'//name_jobstep
!!$          inquire(file=F_OUT, exist=existence)
!!$          if(.not.existence) goto 1000
!!$       end do
!!$1000   continue
!!$    end if
!!$    open(6, file=F_OUT, status='unknown', form='formatted')

  end subroutine m_Parallel_init_comm_world

  subroutine m_Parallel_init_comm_group(mpi_comm_root)
!!$    character :: F_OUT*13, name_mype*3
!!$    character :: name_jobstep*3
!!$    integer ::   i, js
!!$    logical ::   existence
    integer,  intent(in) :: mpi_comm_root
    integer :: inimpi

    logical, save :: firstcall = .true.

    workdir = ''

    inimpi=ON

    if(inimpi==ON) then

      if(.not.firstcall) call mpi_comm_free(MPI_CommGroup, ierr)
      call mpi_comm_dup(mpi_comm_root,MPI_CommGroup,ierr)
      call MPI_barrier(MPI_CommGroup,ierr)
      call mpi_comm_size(MPI_CommGroup, npes, ierr)
      call mpi_comm_rank(MPI_CommGroup, mype, ierr)
                                    __TIMER_START(17)
#ifdef FJ_TIMER
                 call timer_sta(21)
#endif


    else

      call mpi_comm_size(MPI_CommGroup, npes, ierr)
      call mpi_comm_rank(MPI_CommGroup, mype, ierr)

    endif

    firstcall = .false.

  end subroutine m_Parallel_init_comm_group

  subroutine m_Parallel_init_mpi(printable)
    logical, intent(in) :: printable
    integer   :: i
    integer, dimension(8) :: ipresent_time

!!$    call mpi_init(ierr)
!!$    call mpi_comm_size(mpi_comm_world, npes, ierr)
!!$    call mpi_comm_rank(mpi_comm_world, mype, ierr)

    call date_and_time(values=ipresent_time)

    if(printable) then
       write(6,'(" !|| program start    ",i2,":",i2,":",i2,"  ",i2,"/",i2,"/",i4)') &
            & ipresent_time(5:7),ipresent_time(3),ipresent_time(2),ipresent_time(1)

       write(6,'(" npes = ",i6)') npes
       write(6,'(" mype = ",i3," mpi_comm_world = ",i3)') mype,MPI_CommGroup
       write(6,'(" mype = ",i3," mpi_comm_world = ",i3)') mype,MPI_CommGroup
    end if

    call mpi_barrier(MPI_CommGroup,ierr)
    if(printable) write(6,'(" mpi_comm_world = ",i5," mype = ",i5)') MPI_CommGroup, mype
    call mpi_barrier(MPI_CommGroup,ierr)

  end subroutine m_Parallel_init_mpi

  subroutine m_Parallel_get_nproc_from_arg(printable)
    logical, intent(in) :: printable
#ifndef IRIX64
    character(100) :: q1
    integer     :: narg, iargc
#endif
    integer     :: n1, n2, n3
!!!!!!!!!!!!!!!!! modified by mizouchi@adv 2003.02.27 !!!!!!!!!!!!!!
    integer     :: nn , neflag, nkflag, nsflag
!!!!!!!!!!!!!!!!! modified by mizouchi@adv 2003.02.27 !!!!!!!!!!!!!!
    logical     :: wrong_ne, wrong_nk

    if(printable) then
      if(conf_para) then
        write(6,'(a)') ''
        write(6,'(a)') '-- configuration-parallelization scheme enabled --'
        write(6,'(a,i8)')  '   nrank_config : ',nrank_conf
        write(6,'(a,i8)')  '   mype_config  : ',mype_conf
        write(6,'(a,i20)') '   communicator : ',MPI_CommGroup
        write(6,'(a)') ''
      endif
      if(sw_wdir==1)then
        write(6,'(a)') ''
        write(6,'(a)') '-- directory-parallelization scheme enabled --'
        write(6,'(a)') '   current working directory : '//trim(workdir)
        write(6,'(a)') ''
      endif
!!$       write(6,'(" << m_Parallel_get_nproc_from_arg>>")')
       write(6,'(" npes = ",i6," << m_Parallel_get_nproc_from_arg>>")') npes
!!$       write(6,'(" mype = ",i3," mpi_comm_world = ",i3)') mype,mpi_comm_world
!!$       write(6,'(" mype = ",i3," mpi_comm_world = ",i3)') mype,mpi_comm_world
    end if

#ifndef IRIX64

!!!!!!!!!!!!!!!!! modified by mizouchi@adv 2003.02.27 !!!!!!!!!!!!!!
    read_from_args = .false.
    if(mype==0) then
       wrong_ne = .false.
       wrong_nk = .false.
       neflag=0
       nkflag=0
       nsflag=0
       narg = iargc()
!!$       if(printable) write(6,'(" narg = ", i3)') narg

       do nn = 1,narg
          call getarg(nn,q1)

!!!!!!!!!!!!!!!!! modified by mizouchi@adv 2003.03.20 !!!!!!!!!!!!!!
          q1 = trim(adjustl(q1))
!!!!!!!!!!!!!!!!! modified by mizouchi@adv 2003.03.20 !!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!! modified by mizouchi@adv 2003.03.20 !!!!!!!!!!!!!!
!!          if(index(trim(adjustl(q1)),"ne").ne.0) then

          if(q1(1:1) == "n".and. q1(2:2) == "e") then
!!!!!!!!!!!!!!!!! modified by mizouchi@adv 2003.03.20 !!!!!!!!!!!!!!
             neflag = neflag+1
             if((q1(3:3) == "=".or.q1(3:3) == ":") .and. len_trim(q1(4:)).ne.0) then
                 q1 = q1(4:)
                 if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                    read(q1,*) n1
                 else
                    if(printable) write(6,*) 'wrong ne'
                    wrong_ne = .true.
                 end if
             else if(q1(3:3) == "=".or. q1(3:3) == ":") then
                 call getarg(nn+1,q1)
                 if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                    read(q1,*) n1
                 else
                    if(printable) write(6,*) 'wrong ne'
                    wrong_ne = .true.
                 end if
             else
                  call getarg(nn+1,q1)
                  if(q1(1:1) == ":".and. len_trim(q1(2:)).ne.0) then
                      q1 = q1(2:)
                      if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                         read(q1,*) n1
                      else
                         if(printable) write(6,*) 'wrong ne'
                         wrong_ne = .true.
                      end if
                  else if((q1(1:1) == "=".or.q1(1:1) == ":").and. len_trim(q1(2:)).eq.0) then
                      call getarg(nn+2,q1)
                      if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                         read(q1,*) n1
                      else
                         if(printable) write(6,*) 'wrong ne'
                         wrong_ne = .true.
                      end if
                   else
                      if(printable) write(6,*) 'wrong ne'
                      wrong_ne = .true.
                   end if                   
             end if
          end if

!!!!!!!!!!!!!!!!! modified by mizouchi@adv 2003.03.20 !!!!!!!!!!!!!!
!!          if(index(trim(adjustl(q1)),"nk").ne.0) then

          if(q1(1:1) == "n".and. q1(2:2) == "k") then
!!!!!!!!!!!!!!!!! modified by mizouchi@adv 2003.03.20 !!!!!!!!!!!!!!

             nkflag = nkflag+1
             if((q1(3:3) == "=".or.q1(3:3) == ":") .and. len_trim(q1(4:)).ne.0) then
                 q1 = q1(4:)
                 if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                    read(q1,*) n2
                 else
                    if(printable) write(6,*) 'wrong nk'
                    wrong_nk = .true.
                 end if
             else if(q1(3:3) == "=".or.q1(3:3) == ":") then
                 call getarg(nn+1,q1)
                 if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                    read(q1,*) n2
                 else
                    if(printable) write(6,*) 'wrong nk'
                    wrong_nk = .true.
                 end if
             else
                  call getarg(nn+1,q1)
                  if(q1(1:1) == ":".and. len_trim(q1(2:)).ne.0) then
                      q1 = q1(2:)
                      if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                         read(q1,*) n2
                      else
                         if(printable) write(6,*) 'wrong nk'
                         wrong_nk = .true.
                      end if
                  else if((q1(1:1) == "=".or.q1(1:1) == ":").and. len_trim(q1(2:)).eq.0) then
                      call getarg(nn+2,q1)
                      if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                         read(q1,*) n2
                      else
                         if(printable) write(6,*) 'wrong nk'
                         wrong_nk = .true.
                      end if
                   else
                      if(printable) write(6,*) 'wrong nk'
                      wrong_nk = .true.
                   end if                   
             end if
          end if

          if(q1(1:1) == "n".and. q1(2:2) == "s") then

             nsflag = nsflag+1
             if((q1(3:3) == "=".or.q1(3:3) == ":") .and. len_trim(q1(4:)).ne.0) then
                 q1 = q1(4:)
                 if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                    read(q1,*) n3
                 else
                    if(printable) write(6,*) 'wrong ns'
                    !stop "wrong nk "
                    call phase_error_with_msg(6,'wrong ns',__LINE__,__FILE__)
                 end if
             else if(q1(3:3) == "=".or.q1(3:3) == ":") then
                 call getarg(nn+1,q1)
                 if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                    read(q1,*) n3
                 else
                    if(printable) write(6,*) 'wrong ns'
                    !stop "wrong nk "
                    call phase_error_with_msg(6,'wrong ns',__LINE__,__FILE__)
                 end if
             else
                  call getarg(nn+1,q1)
                  if(q1(1:1) == ":".and. len_trim(q1(2:)).ne.0) then
                      q1 = q1(2:)
                      if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                         read(q1,*) n3
                      else
                         if(printable) write(6,*) 'wrong ns'
                         !stop "wrong nk "
                         call phase_error_with_msg(6,'wrong ns',__LINE__,__FILE__)
                      end if
                  else if((q1(1:1) == "=".or.q1(1:1) == ":").and. len_trim(q1(2:)).eq.0) then
                      call getarg(nn+2,q1)
                      if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                         read(q1,*) n3
                      else
                         if(printable) write(6,*) 'wrong ns'
                         !stop "wrong nk "
                         call phase_error_with_msg(6,'wrong ns',__LINE__,__FILE__)
                      end if
                   else
                      if(printable) write(6,*) 'wrong ns'
                      !stop "wrong nk "
                      call phase_error_with_msg(6,'wrong ns',__LINE__,__FILE__)
                   end if
             end if
          end if

       end do

       if(nsflag == 0) n3 = 1
       if(neflag == 0 .or. nkflag == 0) then
!          if(printable) write(6,*) 'set default ne and nk'
          n1 = npes
          n2 = npes/n1
          n3 = 1
       else
          read_from_args = .true.
       end if

       if(neflag >= 2 .or. nkflag >= 2 .or. nsflag >= 2) then
          !if(printable) write(6,*) 'wrong ne and nk'
          !stop "wrong ne and nk "
          call phase_error_with_msg(6,'wrong ne nk and ns',__LINE__,__FILE__)
       end if

    end if

    if(npes>1) call mpi_bcast(neflag,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(npes>1) call mpi_bcast(nkflag,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(npes>1) call mpi_bcast(wrong_nk,1,mpi_logical,0,MPI_CommGroup,ierr)
    if(npes>1) call mpi_bcast(wrong_ne,1,mpi_logical,0,MPI_CommGroup,ierr)

    if(neflag == 0 .or. nkflag == 0) then
!       if(printable) write(6,*) 'set default ne and nk'
       n1 = npes
       n2 = npes/n1
    else
       read_from_args = .true.
    end if

    if(neflag >= 2 .or. nkflag >= 2 .or. (wrong_nk .and. wrong_ne)) then
       !if(printable) write(6,*) 'wrong ne and nk'
       call phase_error_with_msg(6,'wrong ne and nk',__LINE__,__FILE__)
    end if

    if(wrong_nk) then
       call phase_error_with_msg(6,'wrong nk',__LINE__,__FILE__)
    endif

    if(wrong_ne) then
       call phase_error_with_msg(6,'wrong ne',__LINE__,__FILE__)
    endif

!  <---- ! modified by mizouchi@adv 2003.02.27 !!!!!!!!!!!!!!

    if(npes>1) call mpi_bcast(n1,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(npes>1) call mpi_bcast(n2,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(npes>1) call mpi_bcast(n3,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(npes>1) call mpi_bcast(read_from_args,1,mpi_logical,0,MPI_CommGroup,ierr)
!!$    if(printable) write(6,*) 'ne, nk = ',n1,n2
#else
    n1 = 1       ! in-chiki
    n1 = npes/2  ! in-chiki
    if(n1 == 0 .or. n1 > npes) n1 = npes
    n2 = npes/n1
#endif
    if(n1*n2*n3 /= npes) then
       call phase_execution_error(PARALLELIZATION_INVALID_2D)
    else if (n3>2) then
        call phase_error_with_msg(6, 'wrong ns',__LINE__,__FILE__)
    else
       nrank_e  = n1
       nrank_k  = n2
       nrank_s  = n3
       if(printable) then
          if(read_from_args)then
            !!$ if(mype == 0) print '(" nrank_e = ", i3)', nrank_e
            write(6,'(" nrank_e = ",i3)') nrank_e
            !!$ if(mype == 0) print '(" nrank_k = ", i3)', nrank_k
            write(6,'(" nrank_k = ",i3)') nrank_k
            write(6,'(" nrank_s = ",i3)') nrank_s
          endif
       end if
    end if
!!$#ifdef TRANSPOSE
    nrank_g1 = nrank_e
!!$#endif
  end subroutine m_Parallel_get_nproc_from_arg

  subroutine m_Parallel_wd_npes_etc(nfcntn)
    integer, intent(in) :: nfcntn
    if(mype==0) then
       write(nfcntn,*) tag_npes_etc
       write(nfcntn,'(3i10)') npes,nrank_e,nrank_k
    end if
  end subroutine m_Parallel_wd_npes_etc


  subroutine m_Parallel_init_mpi_nbmx(nfout,ipri,printable,nbmx,kg1,ngnode_nbmx,flag_mpi,flag_mpi_k)
    integer, intent(in) :: nfout,ipri,nbmx, kg1, ngnode_nbmx
    logical, intent(in) :: printable,flag_mpi,flag_mpi_k
    integer             :: i, j, ip, icolor, key, np
    integer             :: newpes, newmype
    logical, save       :: firstcall=.true.
                                                  __TIMER_SUB_START(1235)

    if(.not.flag_mpi) then
! (( ng_nbmx, myrank_nbmx, nrank_nbmx, np_nbmx, mp_nbmx, ista_nbmx, iend_nbmx
       ng_nbmx = mype
       myrank_nbmx = 0
       nrank_nbmx  = 1
       np_nbmx = nbmx
       mp_nbmx = nbmx
       ista_nbmx = 1
       iend_nbmx = nbmx
       nbmx_ext = nbmx
       if(ipri >= 2 .and. printable) then
          write(nfout,'(" !|| -- parallelization parameters for nbmx --")')
          write(nfout,'(" !||  - ng_nbmx     = ",i12)') ng_nbmx
          write(nfout,'(" !||  - myrank_nbmx = ",i12)') myrank_nbmx
          write(nfout,'(" !||  - nrank_nbmx  = ",i12)') nrank_nbmx
          write(nfout,'(" !||  - nbmx, nbmx_ext = ",2i12)') nbmx, nbmx_ext
          write(nfout,'(" !||  - np_nbmx     = ",i12)') np_nbmx
          write(nfout,'(" !||  - mp_nbmx     = ",i12)') mp_nbmx
          write(nfout,'(" !||  - ista_nbmx   = ",i12)') ista_nbmx
          write(nfout,'(" !||  - iend_nbmx   = ",i12)') iend_nbmx
       end if
    else
! (( ng_nbmx, myrank_nbmx, nrank_nbmx, np_nbmx, mp_nbmx, ista_nbmx, iend_nbmx
!      , nis_nbmx, nie_nbmx, nel_nbmx, idisp_nbmx))
       ng_nbmx     = floor(mod(mype+1, ngnode_nbmx) - 0.5)
       if(ng_nbmx < 0) ng_nbmx = ng_nbmx + ngnode_nbmx
       myrank_nbmx = floor((mype+0.5)/ngnode_nbmx)
       nrank_nbmx  = (npes - ng_nbmx-1)/ngnode_nbmx + 1

       if(.not.allocated(mpi_nbmx_world)) allocate(mpi_nbmx_world(0:ngnode_nbmx-1))

       allocate(nis_nbmx(0:nrank_nbmx-1))
       allocate(nie_nbmx(0:nrank_nbmx-1))
       allocate(nel_nbmx(0:nrank_nbmx-1))
       allocate(idisp_nbmx(0:nrank_nbmx-1))

       call set_block_range4allgather(nbmx,nrank_nbmx,nel_nbmx,nis_nbmx,nie_nbmx,idisp_nbmx)
       ista_nbmx = nis_nbmx(myrank_nbmx)
       iend_nbmx = nie_nbmx(myrank_nbmx)
       np_nbmx = nel_nbmx(myrank_nbmx)
       mp_nbmx = maxval(nel_nbmx)
       nbmx_ext = mp_nbmx*nrank_nbmx

       if(ipri >= 1 .and. printable) then
          write(nfout,'(" !|| -- parallelization parameters for nbmx --")')
          write(nfout,'(" !||  - ng_nbmx     = ",i12)') ng_nbmx
          write(nfout,'(" !||  - myrank_nbmx = ",i12)') myrank_nbmx
          write(nfout,'(" !||  - nrank_nbmx  = ",i12)') nrank_nbmx
          write(nfout,'(" !||  - nbmx, nbmx_ext = ",2i12)') nbmx, nbmx_ext
          write(nfout,'(" !||  - np_nbmx     = ",i12)') np_nbmx
          write(nfout,'(" !||  - mp_nbmx     = ",i12)') mp_nbmx
          write(nfout,'(" !||  - ista_nbmx   = ",i12)') ista_nbmx
          write(nfout,'(" !||  - iend_nbmx   = ",i12)') iend_nbmx
          write(nfout,'(" !||   ( nis_nbmx )",8i10)') (nis_nbmx(i),i=0,nrank_nbmx-1)
          write(nfout,'(" !||   ( nie_nbmx )",8i10)') (nie_nbmx(i),i=0,nrank_nbmx-1)
          write(nfout,'(" !||   ( nel_nbmx )",8i10)') (nel_nbmx(i),i=0,nrank_nbmx-1)
          write(nfout,'(" !||   (idisp_nbmx)",8i10)') (idisp_nbmx(i),i=0,nrank_nbmx-1)
       end if

! (( mpi_nbmx_world ))
       do j = 0, ngnode_nbmx-1
          icolor = 0
          key = 0
          if(ng_nbmx == j) then
             icolor = 1
             key = myrank_nbmx
          end if
!          if(firstcall) call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_nbmx_world(j),ierr)
          if(.not.firstcall) call mpi_comm_free(mpi_nbmx_world(j), ierr)
          call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_nbmx_world(j),ierr)
          call mpi_comm_size(mpi_nbmx_world(j), newpes, ierr)
          call mpi_comm_rank(mpi_nbmx_world(j), newmype,ierr)
       end do
    end if

    if(.not.flag_mpi_k) then
! (( ng_nbmx_k, myrank_nbmx_k, nrank_nbmx_k, np_nbmx_k, mp_nbmx_k, ista_nbmx_k, iend_nbmx_k,
!    ista_kg1_k, iend_kg1_k ))
       ng_nbmx_k = myrank_k
       myrank_nbmx_k = 0
       nrank_nbmx_k  = 1
       np_nbmx_k = nbmx
       mp_nbmx_k = nbmx
       ista_nbmx_k = 1
       iend_nbmx_k = nbmx

       ista_kg1_k = 1
       iend_kg1_k = kg1
       kg1_ext = kg1
    else
! (( ng_nbmx_k, myrank_nbmx_k, nrank_nbmx_k, np_nbmx_k, mp_nbmx_k, nel_nbmx_k,
!  , ista_nbmx_k, iend_nbmx_k,  nis_kg1_k, nie_kg1_k, nel_kg1_k, ista_kg1_k, iend_kg1_k ))
       ng_nbmx_k     = floor(mod(myrank_e+1, ngnode_nbmx) - 0.5)
       if(ng_nbmx_k < 0) ng_nbmx_k = ng_nbmx_k + ngnode_nbmx
       myrank_nbmx_k = floor((myrank_e+0.5)/ngnode_nbmx)
       nrank_nbmx_k  = (nrank_e - ng_nbmx-1)/ngnode_nbmx + 1

       allocate(nel_nbmx_k(0:nrank_nbmx_k-1))
       allocate(nis_nbmx_k(0:nrank_nbmx_k-1))
       allocate(nie_nbmx_k(0:nrank_nbmx_k-1))
       allocate(idisp_nbmx_k(0:nrank_nbmx_k-1))
       if(.not.allocated(mpi_nbmx_world_k)) then
       allocate(mpi_nbmx_world_k(0:ngnode_nbmx-1))
       endif

       call set_block_range4allgather(nbmx,nrank_nbmx_k,nel_nbmx_k,nis_nbmx_k,nie_nbmx_k,idisp_nbmx_k)
       ista_nbmx_k = nis_nbmx_k(myrank_nbmx_k)
       iend_nbmx_k = nie_nbmx_k(myrank_nbmx_k)
       np_nbmx_k = nel_nbmx_k(myrank_nbmx_k)
       mp_nbmx_k = maxval(nel_nbmx_k)
       np = mp_nbmx_k*nrank_nbmx_k
       if(np .gt. nbmx_ext) nbmx_ext = np

! (( np_kg1_k, mp_kg1_k, nel_kg1_k, ista_kg1_k, iend_kg1_k, idisp_kg1_k ))
       allocate(nel_kg1_k(0:nrank_nbmx_k-1))
       allocate(nis_kg1_k(0:nrank_nbmx_k-1))
       allocate(nie_kg1_k(0:nrank_nbmx_k-1))
       allocate(idisp_kg1_k(0:nrank_nbmx_k-1))

       call set_block_range4allgather(kg1,nrank_nbmx_k,nel_kg1_k,nis_kg1_k,nie_kg1_k,idisp_kg1_k)
       np_kg1_k = nel_kg1_k(myrank_nbmx_k)
       mp_kg1_k = maxval(nel_kg1_k)
       ista_kg1_k = nis_kg1_k(myrank_nbmx_k)
       iend_kg1_k = nie_kg1_k(myrank_nbmx_k)
       kg1_ext = mp_kg1_k*nrank_nbmx_k

       if(ipri >= 1 .and. printable) then
          write(nfout,'(" !|| -- parallelization parameters for nbmx in mpi_k_world--")')
          write(nfout,'(" !||  - ng_nbmx_k     = ",i12)') ng_nbmx_k
          write(nfout,'(" !||  - myrank_nbmx_k = ",i12)') myrank_nbmx_k
          write(nfout,'(" !||  - nrank_nbmx_k  = ",i12)') nrank_nbmx_k
          write(nfout,'(" !||  - nbmx, nbmx_ext = ",2i12)') nbmx, nbmx_ext
          write(nfout,'(" !||  - np_nbmx_k     = ",i12)') np_nbmx_k
          write(nfout,'(" !||  - mp_nbmx_k     = ",i12)') mp_nbmx_k
          write(nfout,'(" !||  - ista_nbmx_k   = ",i12)') ista_nbmx_k
          write(nfout,'(" !||  - iend_nbmx_k   = ",i12)') iend_nbmx_k
          write(nfout,'(" !||   ( nis_nbmx_k )",8i10)') (nis_nbmx_k(i),i=0,nrank_nbmx_k-1)
          write(nfout,'(" !||   ( nie_nbmx_k )",8i10)') (nie_nbmx_k(i),i=0,nrank_nbmx_k-1)
          write(nfout,'(" !||   ( nel_nbmx_k )",8i10)') (nel_nbmx_k(i),i=0,nrank_nbmx_k-1)
          write(nfout,'(" !||   (idisp_nbmx_k)",8i10)') (idisp_nbmx_k(i),i=0,nrank_nbmx_k-1)
          write(nfout,'(" !||  - kg1, kg1_ext  = ",2i12)') kg1, kg1_ext
          write(nfout,'(" !||  - np_kg1_k      = ",i12)') np_kg1_k
          write(nfout,'(" !||  - mp_kg1_k      = ",i12)') mp_kg1_k
          write(nfout,'(" !||  - ista_kg1_k    = ",i12)') ista_kg1_k
          write(nfout,'(" !||  - iend_kg1_k    = ",i12)') iend_kg1_k
          write(nfout,'(" !||   ( nis_kg1_k )",8i10)') (nis_kg1_k(i),i=0,nrank_nbmx_k-1)
          write(nfout,'(" !||   ( nie_kg1_k )",8i10)') (nie_kg1_k(i),i=0,nrank_nbmx_k-1)
          write(nfout,'(" !||   ( nel_kg1_k )",8i10)') (nel_kg1_k(i),i=0,nrank_nbmx_k-1)
          write(nfout,'(" !||   (idisp_kg1_k)",8i10)') (idisp_kg1_k(i),i=0,nrank_nbmx_k-1)
       end if

! (( mpi_nbmx_world ))
       do i = 0, nrank_k-1
          if(i == myrank_k) then
             do j = 0, ngnode_nbmx-1
                icolor = 0
                key = 0
                if(ng_nbmx_k == j) then
                   icolor = 1
                   key = myrank_nbmx_k
                end if
!                if(firstcall) call mpi_comm_split(mpi_k_world(myrank_k), icolor, key, mpi_nbmx_world_k(j),ierr)
                if(.not.firstcall) call mpi_comm_free(mpi_nbmx_world_k(j), ierr)
                call mpi_comm_split(mpi_k_world(myrank_k), icolor, key, mpi_nbmx_world_k(j),ierr)
                call mpi_comm_size(mpi_nbmx_world_k(j), newpes, ierr)
                call mpi_comm_rank(mpi_nbmx_world_k(j), newmype,ierr)
             end do
          end if
       end do
    end if
    firstcall = .false.
                                                  __TIMER_SUB_STOP(1235)
  end subroutine m_Parallel_init_mpi_nbmx

  subroutine m_Parallel_init_mpi_urec_hsr(nfout,nsize_rho_hsr)
    ! Coded by T. Ymasaki, 2023/07/07
    integer, intent(in) :: nfout, nsize_rho_hsr
    integer, allocatable, dimension(:) :: is_hsr, ie_hsr
    integer :: iwork, k
    allocate(is_hsr(0:npes-1),ie_hsr(0:npes-1))
    iwork = ( nsize_rho_hsr - 1)/npes + 1
    do k = 0, npes-1
       is_hsr(k) = min(k*iwork+1,nsize_rho_hsr+1)
       ie_hsr(k) = min(is_hsr(k)+iwork-1,nsize_rho_hsr)
    end do
    ista_urec_hsr = is_hsr(mype)
    iend_urec_hsr = ie_hsr(mype)
    write(nfout,'(" ista_urec_hsr = ",i8, " iend_urec_hsr = ",i8)') ista_urec_hsr, iend_urec_hsr
    call flush(nfout)
    deallocate(is_hsr, ie_hsr)
    ista_and_iend_urec_hsr_set = .true.
  end subroutine m_Parallel_init_mpi_urec_hsr
  
  subroutine set_block_range4allgather(ne,np,nel_p,nis_p,nie_p,idisp_p)
    integer, intent(in)                     :: ne ! number of total elements
    integer, intent(in)                     :: np ! number of ranks (or processors)
    integer, intent(out), dimension(0:np-1) :: nel_p,nis_p,nie_p,idisp_p

    integer :: i, j, npf

    nel_p = 0
    i = ne/np
    j = i*np
    if(j == ne) then
       nel_p = i
    else
       npf = ne/(i+1)
       nel_p(0:npf-1) = i+1
       nel_p(npf)     = ne - (i+1)*npf
    end if
    nis_p(0) = 1
    do i = 1, np-1
       nis_p(i)   = nis_p(i-1) + nel_p(i-1)
       nie_p(i-1) = nis_p(i)   - 1
    end do
    nie_p(np-1) = ne
    do i = 0, np-1
       idisp_p(i) = nis_p(i)-1
    end do
  end subroutine set_block_range4allgather

! === necessary to make 3D_Parallel, too!!! by tkato ===========================
! ==============================================================================
  subroutine m_Parallel_init_mpi_elec(nfout,ipri,printable,neg,kv3,nspin,kg1)
    integer, intent(in) :: nfout,ipri,neg, kv3, nspin, kg1
    logical, intent(in) :: printable
    logical, save       :: firstcall=.true.

#ifdef NEC_TUNE_SOFT
    character*4 F_RSVTASK
    integer :: tmp_a1,tmp_a2,tmp_b1,tmp_b2
#elif NEC_TUNE_FFT
    character*4 F_RSVTASK
#endif

    integer             :: i, j, k, ip, icolor, key, kv3_half, ispin

!!$! =========================================== added by K. Tagami ======== 11.0
!!$    integer ::           kv3_div
!!$! ======================================================================= 11.0

    integer             :: newpes, newmype
    logical             :: set_mapping_func
    integer, pointer, dimension(:) :: nproc2rank_k, nproc2rank_e, nproc2rank_s

    integer, parameter  :: nc = 24
    integer, parameter  :: nc1 = 16
!!$    integer, parameter  :: nc2 = 18
    integer             :: i0, i1, sw_title
    character*7         :: strmap

    sw_title = 1
    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| << m_Parallel_init_mpi_elec >>")')
       write(nfout,'(" !|| neg, kv3 = ",i5, i8)') neg, kv3
    end if
    if(neg <= 0 ) call phase_error_with_msg(nfout," neg is not positive value ",__LINE__,__FILE__)
    if(kv3 <= 0 ) call phase_error_with_msg(nfout," kv3 is not positive value ",__LINE__,__FILE__)
    if(nrank_s>nspin)  &
    call phase_error_with_msg(nfout,"ns>nspin",__LINE__,__FILE__)

!!!!!!!!!!!!!!!!! added by mizouchi@adv 2003.02.26 !!!!!!!!!!!!!!
#ifdef IRIX64

! ====================================== modified by K. Tagami =============== 11.0
!!    if((kv3 == 1.and.nspin ==1).or.(kv3 == 2.and.nspin ==2)) then
    if ( kv3/nspin == 1 ) then
! ===========================================================================  11.0

       nrank_e  = npes
       nrank_k  = 1
       if(ipri >= 1 .and. printable) then
          write(nfout,'(" !| modified nrank_e and nrank_k due to 1 kpoint calculation  ")')
          write(nfout,'(" !| modified nrank_e = ", i4)') nrank_e
          write(nfout,'(" !| modified nrank_k = ", i4)') nrank_k
       end if
#endif
!!!!!!!!!!!!!!!!! added by mizouchi@adv 2003.02.26 !!!!!!!!!!!!!!

!    if(nrank_k > kv3) stop ' nrank_k > kv3 (m_Parallel_init_mpi_elec)'
    if(nrank_k > kv3) call phase_execution_error(PARALLELIZATION_INVALID_NK)
!!$    if(nrank_e > neg) call phase_execution_error(PARALLELIZATION_INVALID_NE) !! not here !!

    allocate(map_e(neg))
    allocate(map_z(neg))
    allocate(map_k(kv3))
    allocate(map_ek(neg,kv3))
    allocate(map_eks(neg,kv3,nspin))
    allocate(map_s(nspin))
    if(.not.allocated(mpi_k_world)) allocate(mpi_k_world(0:nrank_k-1))
    if(.not.allocated(mpi_e_world)) allocate(mpi_e_world(0:nrank_e-1))
    allocate(nis_e(0:nrank_e-1)); nis_e = neg+1
    allocate(nie_e(0:nrank_e-1)); nie_e = 0
    allocate(nel_e(0:nrank_e-1)); nel_e = 0
    allocate(idisp_e(0:nrank_e-1)); idisp_e = 0
!!$#ifdef TRANSPOSE
    allocate(nis_g1(0:nrank_g1-1)); nis_g1 = kg1 + 1
    allocate(nie_g1(0:nrank_g1-1)); nie_g1 = 0
    allocate(nel_g1(0:nrank_g1-1)); nel_g1 = 0
!!$#endif
    allocate(nis_k(0:nrank_k-1)); nis_k = 0
    allocate(nie_k(0:nrank_k-1))
    allocate(nel_k(0:nrank_k-1))

    allocate(nis_spin(0:nrank_s-1)); nis_spin = 0
    allocate(nie_spin(0:nrank_s-1))
    allocate(nel_spin(0:nrank_s-1))

! (( myrank_e, myrank_k ))
    allocate(nproc2rank_e(0:npes-1))
    allocate(nproc2rank_k(0:npes-1))
    allocate(nproc2rank_s(0:npes-1))
    ip = 0
    do k = 0, nrank_s-1
       do i = 0, nrank_k-1
          do j = 0, nrank_e-1
             nproc2rank_k(ip) = i
             nproc2rank_e(ip) = j
             nproc2rank_s(ip) = k
             ip = ip + 1
          end do
       end do
    enddo
    myrank_k = nproc2rank_k(mype)
    myrank_e = nproc2rank_e(mype)
    myrank_spin = nproc2rank_s(mype)
    deallocate(nproc2rank_k)
    deallocate(nproc2rank_e)
    deallocate(nproc2rank_s)
! (( map_e, nis_e, nie_e, nel_e, ista_e, iend_e, istep_e, np_e, mp_e  ))
!!$#ifdef TRANSPOSE
    set_mapping_func = .true.
    call set_block_range(neg,nrank_e,nel_e,nis_e,nie_e, set_mapping_func, map_e)
    do i = 0, nrank_e-1
       idisp_e(i) = nis_e(i)-1
    end do
    if(ipri >= 2) call wd_e_range
    istep_e = 1
!!$#else
!!$    do i = 1, neg
!!$       map_e(i) = mod(i-1,nrank_e)
!!$    end do
!!$    do i = 1, neg
!!$       ip = map_e(i)
!!$       nel_e(ip) = nel_e(ip) + 1
!!$       if(nis_e(ip) > i) nis_e(ip) = i
!!$       if(nie_e(ip) < i) nie_e(ip) = i
!!$    end do
!!$    istep_e = nrank_e
!!$
!!$#endif
    if(ipri >= 2 .and. printable) then
       strmap = "   i   "
       call wd_maparray(neg,map_e,strmap," map_e ",sw_title)
    end if

    ista_e = nis_e(myrank_e)
    iend_e = nie_e(myrank_e)
    np_e = nel_e(myrank_e)
    mp_e = maxval(nel_e)
    if(ipri==1 .and. printable) then
       if(mp_e < 100000) then
          write(nfout,'(" !|| ista_e, iend_e = ",i5,",",i5,", mp_e = ",i5,", neg = ",i12,", kv3 = ",i5 &
               & ," << m_Parallel_init_mpi_elec >>")') ista_e,iend_e,mp_e,neg,kv3
       else if(mp_e <100000000) then
          write(nfout,'(" !|| ista_e, iend_e = ",i8,",",i8,", mp_e = ",i8,", neg = ",i13,", kv3 = ",i5 &
               & ," << m_Parallel_init_mpi_elec >>")') ista_e,iend_e,mp_e,neg,kv3
       else
          write(nfout,'(" !|| ista_e, iend_e = ",i0,",",i0,", mp_e = ",i0,", neg = ",i0,", kv3 = ",i5 &
               & ," << m_Parallel_init_mpi_elec >>")') ista_e,iend_e,mp_e,neg,kv3
       end if
    end if


#ifdef NEC_TUNE_SOFT
    call getenv('F_RSVTASK',F_RSVTASK)
    read (F_RSVTASK,'(i4)') itask
    allocate(ista_e_smp(itask))
    allocate(iend_e_smp(itask))

    tmp_a1 = (iend_e-ista_e+1)/istep_e
    tmp_a2 = mod(iend_e-ista_e+1,istep_e)
    if (tmp_a2 .ne. 0) tmp_a1=tmp_a1+1
    tmp_b1 = tmp_a1/itask
    tmp_b2 = mod(tmp_a1,itask)
    if (tmp_b2 .ne. 0) tmp_b1=tmp_b1+1

    do i=1,itask
      ista_e_smp(i) = ista_e + tmp_b1 * istep_e * (i -1)
      iend_e_smp(i) = min((ista_e + tmp_b1 * istep_e * i -1),iend_e)
    end do
#elif NEC_TUNE_FFT
    call getenv('F_RSVTASK',F_RSVTASK)
    read (F_RSVTASK,'(i4)') itask
#endif

    if(ipri >= 2 .and. printable) &
         & write(nfout,'(" !|| -- ista_e, iend_e, istep_e, np_e, mp_e = ",5i6)') &
         &  ista_e,iend_e,istep_e, np_e, mp_e

! (( map_z ))
!!$#ifdef TRANSPOSE
    j = 0
    do ip = 1, nrank_e
       do i = 1, nel_e(ip-1)
          j = j + 1
          map_z(j) = i
       end do
    end do
!!$#else
!!$    ip = 1
!!$    do i = 1, neg
!!$       map_z(i) = ip
!!$       if(mod(i,nrank_e) == 0) ip = ip + 1
!!$    end do
!!$#endif
    if(ipri >= 2 .and. printable) then
       strmap = "   i   "
       call wd_maparray(neg,map_z,strmap," map_z ",sw_title)
    end if
! (( map_k, nis_k, nie_k, nel_k, ista_k, iend_k ))
    if(kv3/nspin < nrank_k) then
       if(ipri >= 1 .and. printable) then
          write(nfout,'(" ******")')
          write(nfout,'(" ** The nrank_k that you have specified is smaller than the number of k-points (kv3/nspin).")')
          write(nfout,'(" ** kv3/nspin = ",i8," nrank_k = ",i8)') kv3/nspin, nrank_k
          write(nfout,'(" ** Reduce the size of the rank for k-point division (nrank_k) &
  &                     which is set in nk= option when you run this program, please")')
          call flush(nfout)
       end if
       !stop ' ** The nrank_k that is specified is smaller than the number of generated or read k-points (=kv3/nspin).'
       call phase_error_with_msg(nfout, &
       'The nrank_k that is specified is smaller than the number of generated or read k-points (=kv3/nspin).', &
       __LINE__, __FILE__)
    end if
    if(nspin == 1) then
       set_mapping_func = .true.
       call set_block_range(kv3,nrank_k,nel_k,nis_k,nie_k, set_mapping_func,map_k)
       ista_k = nis_k(myrank_k)
       iend_k = nie_k(myrank_k)
    else if(nspin == 2) then
       set_mapping_func = .false.
       kv3_half = kv3/nspin
       call set_block_range(kv3_half,nrank_k,nel_k,nis_k,nie_k, set_mapping_func)
       nel_k = nel_k*nspin
       nie_k = nie_k*nspin
       nis_k = nis_k*nspin-1
       j = 0
       do i = 1, kv3
          if(nie_k(j) < i) j = j + 1
          map_k(i) = j
       end do
       ista_k = nis_k(myrank_k)
       iend_k = nie_k(myrank_k)
    end if

    set_mapping_func = .true.
    call set_block_range(nspin,nrank_s,nel_spin,nis_spin,nie_spin, set_mapping_func, map_s)
    ista_spin = nis_spin(myrank_spin)
    iend_spin = nie_spin(myrank_spin)
    np_spin   = nel_spin(myrank_spin)
    mp_spin   = np_spin
    do i = 0, nrank_e-1
       idisp_e(i) = nis_e(i)-1
    end do
    if(ipri >= 2) call wd_e_range
    istep_e = 1
!!$#else
!!$    do i = 1, neg
!!$       map_e(i) = mod(i-1,nrank_e)
!!$    end do
!!$    do i = 1, neg
!!$       ip = map_e(i)
!!$       nel_e(ip) = nel_e(ip) + 1
!!$       if(nis_e(ip) > i) nis_e(ip) = i
!!$       if(nie_e(ip) < i) nie_e(ip) = i
!!$    end do
!!$    istep_e = nrank_e
!!$
!!$#endif
    if(ipri >= 2 .and. printable) then
       strmap = "   i   "
       call wd_maparray(neg,map_e,strmap," map_e ",sw_title)
    end if

    ista_e = nis_e(myrank_e)
    iend_e = nie_e(myrank_e)
    np_e = nel_e(myrank_e)
    mp_e = maxval(nel_e)
    if(ipri==1 .and. printable) then
       if(mp_e < 100000) then
          write(nfout,'(" !|| ista_e, iend_e = ",i5,",",i5,", mp_e = ",i5,", neg = ",i12,", kv3 = ",i5 &
               & ," << m_Parallel_init_mpi_elec >>")') ista_e,iend_e,mp_e,neg,kv3
       else if(mp_e <100000000) then
          write(nfout,'(" !|| ista_e, iend_e = ",i8,",",i8,", mp_e = ",i8,", neg = ",i13,", kv3 = ",i5 &
               & ," << m_Parallel_init_mpi_elec >>")') ista_e,iend_e,mp_e,neg,kv3
       else
          write(nfout,'(" !|| ista_e, iend_e = ",i0,",",i0,", mp_e = ",i0,", neg = ",i0,", kv3 = ",i5 &
               & ," << m_Parallel_init_mpi_elec >>")') ista_e,iend_e,mp_e,neg,kv3
       end if
    end if

    if(ipri >= 2 .and. printable) then
       strmap = "   i   "
       call wd_maparray(kv3,map_k,strmap," map_k ",sw_title)
       write(nfout,'(" !|| -- myrank_k = ",i8)') myrank_k
       write(nfout,'(" !|| -- ista_k, iend_k = ",2i8)') ista_k,iend_k
    end if
   if(ista_k > iend_k) call phase_error_with_msg(nfout, ' !! illegal combination of ista_k, and iend_k', &
   __LINE__,__FILE__)

! (( map_ek ))
    do j = 1, kv3
       do i = 1, neg
          map_ek(i,j) = map_e(i) + map_k(j)*(nrank_e)
       end do
    end do
! (( map_eks ))
    do ispin = 1, nspin
    !do j = 1, kv3
    do j = ispin, kv3-nspin+ispin, nspin
       do i = 1, neg
          map_eks(i,j,ispin) = map_e(i) + map_k(j)*nrank_e + map_s(ispin)*nrank_e*nrank_k
       end do
    end do
    end do

    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| --- map_ek ---")')
       sw_title = 0
       do j = 1, kv3
          write(strmap,'("ik=",i4)') j
          call wd_maparray(neg,map_ek(1,j),strmap,"map_ek ",sw_title)
       end do
    end if

!!$#ifdef TRANSPOSE
! (( ista_g1, iend_g1 ))
    set_mapping_func = .false.
    if(ipri >= 2 .and. printable) write(nfout,'(" !|| kg1 = ",i10)') kg1
    call set_block_range(kg1,nrank_e,nel_g1,nis_g1,nie_g1, set_mapping_func)
    ista_g1 = nis_g1(myrank_e)
    iend_g1 = nie_g1(myrank_e)
    np_g1   = nel_g1(myrank_e)
    mp_g1   = maxval(nel_g1)
    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| --- mis_g1,nie_g1,nel_g1 ---")')
       write(nfout,'(" !|| ( rank_e  )",15i7)')(i,i=1,nrank_e)
       write(nfout,'(" !|| ( nis_g1  )",15i7)')(nis_g1(i),i=0,nrank_e-1)
       write(nfout,'(" !|| ( nie_g1  )",15i7)')(nie_g1(i),i=0,nrank_e-1)
       write(nfout,'(" !|| ( nel_g1  )",15i7)')(nel_g1(i),i=0,nrank_e-1)
       write(nfout,'(" !|| (myrank_e = ",i3," kg1,ista_g1,iend_g1,np_g1,mp_g1 = ",5i7)') &
            & myrank_e,kg1,ista_g1,iend_g1,np_g1,mp_g1
    end if
!!$#endif

! (( mpi_k_world ))
    icolor = myrank_spin
    key    = myrank_e + myrank_k*nrank_e
    if(.not.firstcall) call mpi_comm_free(mpi_spin_group, ierr)
    call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_spin_group,ierr)
    call mpi_comm_size(mpi_spin_group, newpes,ierr)
    call mpi_comm_rank(mpi_spin_group, newmype,ierr)
    do j = 0, nrank_k-1
       icolor = 0
       key = 0
       do i = 0, nrank_e-1
          if(mype == i + nrank_e*j + myrank_spin*nrank_e*nrank_k) then
             icolor = 1
             key = i
          end if
       end do
!       if(firstcall) call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_k_world(j),ierr)
       if(.not.firstcall) call mpi_comm_free(mpi_k_world(j), ierr)
!       call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_k_world(j),ierr)
       call mpi_comm_split(mpi_spin_group, icolor, key, mpi_k_world(j),ierr)
       call mpi_comm_size(mpi_k_world(j), newpes,ierr)
       call mpi_comm_rank(mpi_k_world(j), newmype,ierr)
!!$       print '(" mype = ",i3," newmype = ", i3," newpes = ",i3," mpi_k_world = ", i3 &
!!$            &, " icolor = ",i3," key = ",i3)' &
!!$            & ,mype,newmype, newpes,mpi_k_world(j),icolor, key
    end do

    do j=0,nrank_e-1
       icolor = myrank_e
       key = mype
!       if(firstcall) call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_e_world(j),ierr)
       if(.not.firstcall) call mpi_comm_free(mpi_e_world(j), ierr)
       call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_e_world(j),ierr)
       call mpi_comm_size(mpi_e_world(j), newpes,ierr)
       call mpi_comm_rank(mpi_e_world(j), newmype,ierr)
    enddo
    icolor = myrank_e
    key    = mype
    call mpi_comm_split(MPI_CommGroup, icolor, key, mpi_ge_world, ierr)
    call mpi_comm_size(mpi_ge_world, newpes, ierr)
    call mpi_comm_rank(mpi_ge_world, newmype, ierr)
    firstcall = .false.
  contains
    subroutine wd_maparray(nelement,map_a,map_title0,map_title,sw_title)
      integer, intent(in) ::                 nelement
      integer, dimension(nelement), intent(in) :: map_a
      character*7, intent(in) ::             map_title0,map_title
      integer, intent(in) ::                 sw_title
      integer :: nca, j, i0, i1, i, iloop

      if(sw_title==1)  write(nfout,'(" !|| --- ",a7," ---")') map_title
      if(nelement < 100) then
         nca = nc      ! = 24
      else if(nelement < 1000) then
         nca = nc*3/4  ! = 18
      else if(nelement < 10000) then
         nca = nc*3/5  ! = 14
      else if(nelement < 100000) then
         nca = nc*3/6  ! = 12
      else
         nca = nc*3/9  ! = 8
      end if
      iloop = ceiling(dble(nelement)/nca)
      do j = 1, iloop
         i0 = (j-1)*nca + 1
         i1 = i0 + nca -1
         if(i1 > nelement) i1 = nelement
         if(nelement < 100) then
            write(nfout,'(" !|| (",a7,")",24i3)') map_title0,(i,i=i0,i1)
            write(nfout,'(" !|| (",a7,")",24i3)') map_title, (map_a(i),i=i0,min(nelement,i1))
         else if(nelement < 1000) then
            write(nfout,'(" !|| (",a7,")",18i4)') map_title0,(i,i=i0,i1)
            write(nfout,'(" !|| (",a7,")",18i4)') map_title, (map_a(i),i=i0,i1)
         else if(nelement < 10000) then
            write(nfout,'(" !|| (",a7,")",14i5)') map_title0,(i,i=i0,i1)
            write(nfout,'(" !|| (",a7,")",14i5)') map_title, (map_a(i),i=i0,i1)
         else if(nelement < 100000) then
            write(nfout,'(" !|| (",a7,")",12i5)') map_title0,(i,i=i0,i1)
            write(nfout,'(" !|| (",a7,")",12i5)') map_title, (map_a(i),i=i0,i1)
         else
            write(nfout,'(" !|| (",a7,")",8i9)') map_title0,(i,i=i0,i1)
            write(nfout,'(" !|| (",a7,")",8i9)') map_title, (map_a(i),i=i0,i1)
         end if
      end do
    end subroutine wd_maparray

    subroutine wd_e_range
      if(printable) then
         write(nfout,'(" !|| nrank_e")')
         write(nfout,'(" !||    i    : ",20i4)')(i,i=0,nrank_e-1)
         write(nfout,'(" !|| nis_e   : ",20i4)')(nis_e(i),i=0,nrank_e-1)
         write(nfout,'(" !|| nie_e   : ",20i4)')(nie_e(i),i=0,nrank_e-1)
         write(nfout,'(" !|| nel_e   : ",20i4)')(nel_e(i),i=0,nrank_e-1)
!!$         write(nfout,'(" !||    i : ",20i4)')(i,i=0,nrank_e-1)
!!$         write(nfout,'(" !|| nis_e : ",20i4)')(nis_e(i),i=0,nrank_e-1)
!!$         write(nfout,'(" !|| nie_e : ",20i4)')(nie_e(i),i=0,nrank_e-1)
!!$         write(nfout,'(" !|| nel_e : ",20i4)')(nel_e(i),i=0,nrank_e-1)
!!$#ifdef TRANSPOSE
         write(nfout,'(" !|| idisp_e : ",20i4)')(idisp_e(i),i=0,nrank_e-1)
!!$#endif
      end if
    end subroutine wd_e_range

  end subroutine m_Parallel_init_mpi_elec
! === necessary to make 3D_Parallel, too!!! by tkato ===========================

! ==============================================================================
  subroutine m_Parallel_dealloc_mpi_elec()

    if(allocated(map_e)) deallocate(map_e)
    if(allocated(map_z)) deallocate(map_z)
    if(allocated(map_k)) deallocate(map_k)
    if(allocated(map_ek)) deallocate(map_ek)
    if(allocated(map_eks)) deallocate(map_eks)
!fj    deallocate(mpi_k_world)
    if(allocated(nis_e)) deallocate(nis_e)
    if(allocated(nie_e)) deallocate(nie_e)
    if(allocated(nel_e)) deallocate(nel_e)
    if(allocated(idisp_e)) deallocate(idisp_e)
    if(allocated(nis_g1)) deallocate(nis_g1)
    if(allocated(nie_g1)) deallocate(nie_g1)
    if(allocated(nel_g1)) deallocate(nel_g1)
    if(allocated(nis_k)) deallocate(nis_k)
    if(allocated(nie_k)) deallocate(nie_k)
    if(allocated(nel_k)) deallocate(nel_k)
    if(allocated(nis_spin)) deallocate(nis_spin)
    if(allocated(nie_spin)) deallocate(nie_spin)
    if(allocated(nel_spin)) deallocate(nel_spin)
    if(allocated(map_s)) deallocate(map_s)

#ifdef NEC_TUNE_SOFT
    if(allocated(ista_e_smp)) deallocate(ista_e_smp)
    if(allocated(iend_e_smp)) deallocate(iend_e_smp)
#endif


  end subroutine m_Parallel_dealloc_mpi_elec

  subroutine m_Parallel_init_mpi_iba(nfout,ipri,printable,kv3,iba)
    integer, intent(in) :: nfout,ipri, kv3, iba(kv3)
    logical, intent(in) :: printable

    logical             :: set_mapping_func
    integer             :: i, ik
    integer, parameter  :: nc = 24
    integer, parameter  :: nc1 = 16

!!$#ifdef TRANSPOSE
    if(nrank_k > kv3) call phase_error_with_msg(nfout, ' nrank_k > kv3 (m_Parallel_init_mpi_elec)', __LINE__, __FILE__)

! (( ista_g1k, iend_g1k ))
    if(allocated(np_g1k)) deallocate(np_g1k);  allocate(np_g1k(kv3))
    if(allocated(mp_g1k)) deallocate(mp_g1k);  allocate(mp_g1k(kv3))
    if(allocated(ista_g1k)) deallocate(ista_g1k); allocate(ista_g1k(kv3))
    if(allocated(iend_g1k)) deallocate(iend_g1k); allocate(iend_g1k(kv3))

    if(allocated(nis_g1k)) deallocate(nis_g1k);  allocate(nis_g1k(0:nrank_g1-1,kv3))
    do ik = 1, kv3
       nis_g1k(:,ik) = iba(ik)
    end do

    if(allocated(nie_g1k)) deallocate(nie_g1k); allocate(nie_g1k(0:nrank_g1-1,kv3))
    nie_g1k = 0

    if(allocated(nel_g1k)) deallocate(nel_g1k); allocate(nel_g1k(0:nrank_g1-1,kv3))
    nel_g1k = 0

    do ik = 1, kv3
       set_mapping_func = .false.
       if((ipri >= 2 .and. printable).or. &
            & (ipri >= 1 .and. printable .and. (ik<=3 .or. ik>=kv3-2))) &
            & write(nfout,'(" !|| ik = ",i10 &
            & ," iba(",i6,") = ",i10," <<m_Parallel_init_mpi_iba>>")') ik, ik, iba(ik)

       call set_block_range(iba(ik),nrank_e,nel_g1k(0,ik),nis_g1k(0,ik) &
            & ,nie_g1k(0,ik), set_mapping_func)
       ista_g1k(ik) = nis_g1k(myrank_e,ik)
       iend_g1k(ik) = nie_g1k(myrank_e,ik)
       np_g1k(ik)   = nel_g1k(myrank_e,ik)
       mp_g1k(ik)   = maxval(nel_g1k(0:nrank_g1-1,ik))
       if((ipri >= 3 .and. printable) .or. &
            & (ipri >= 2 .and. printable .and. (ik<=3 .or. ik>=kv3-2))) then
          write(nfout,'(" !|| --- mis_g1k,nie_g1k,nel_g1k ---")')
          write(nfout,'(" !|| ( rank_e  )",15i7)')(i,i=1,nrank_e)
          write(nfout,'(" !|| ( nis_g1k )",15i7)')(nis_g1k(i,ik),i=0,nrank_e-1)
          write(nfout,'(" !|| ( nie_g1k )",15i7)')(nie_g1k(i,ik),i=0,nrank_e-1)
          write(nfout,'(" !|| ( nel_g1k )",15i7)')(nel_g1k(i,ik),i=0,nrank_e-1)
          write(nfout,'(" !|| (myrank_e = ",i3," iba(ik),ista_g1k,iend_g1k,np_g1k,mp_g1k = ",5i7)') &
               & myrank_e,iba(ik),ista_g1k(ik),iend_g1k(ik),np_g1k(ik),mp_g1k(ik)
       end if
       if(ipri>=1 .and. printable .and. kv3 >=7 .and. ik == 3) then
          write(nfout,'(" !||  ......")')
       end if
    end do
!!$#endif
  end subroutine m_Parallel_init_mpi_iba

  subroutine m_Parallel_end_mpi
                                                  __TIMER_STOP
    call mpi_finalize(ierr)
  end subroutine m_Parallel_end_mpi

  subroutine set_block_range(ne,np,nel_p,nis_p,nie_p,set_mapping_func,map_p)
    integer, intent(in)                     :: ne ! number of total elements
    integer, intent(in)                     :: np ! number of ranks (or processors)
    integer, intent(out), dimension(0:np-1) :: nel_p,nis_p,nie_p
    logical, intent(in)                     :: set_mapping_func
    integer, intent(out), optional, dimension(ne)  :: map_p

    integer :: j,i

    if(np == 0) call phase_error_with_msg(6,' np == 0',__LINE__,__FILE__)
    nel_p = ne/np
    j = mod(ne,np)
    do i = 0, j-1
       nel_p(i) = nel_p(i) + 1
    end do

    nis_p(0) = 1
    do i = 1, np-1
       nis_p(i)   = nis_p(i-1) + nel_p(i-1)
       nie_p(i-1) = nis_p(i) - 1
    end do
    nie_p(np-1) = ne

    if(set_mapping_func) then
       j = 0
       do i = 1, ne
          if(nie_p(j) < i) j = j + 1
          map_p(i) = j
       end do
    end if

  end subroutine set_block_range

  subroutine m_Parallel_init_mpi_kngp(nfout,ipri,kngp)
    integer, intent(in) :: nfout,ipri,kngp
    integer :: iwork, i
    if(allocated(is_kngp)) deallocate(is_kngp)
    if(allocated(ie_kngp)) deallocate(ie_kngp)
    if(allocated(nel_kngp)) deallocate(nel_kngp)
    allocate(is_kngp(0:npes-1))
    allocate(ie_kngp(0:npes-1))
    allocate(nel_kngp(0:npes-1))
    iwork = ( kngp - 1 ) / npes + 1
    do i = 0, npes-1
       is_kngp(i) = min(i*iwork+1, kngp+1)
       ie_kngp(i) = min(is_kngp(i)+iwork-1, kngp)
       nel_kngp(i) = ie_kngp(i) - is_kngp(i) + 1
    enddo
    ista_kngp = is_kngp(mype)
    iend_kngp = ie_kngp(mype)
    np_kngp   = nel_kngp(mype)
    mp_kngp   = maxval(nel_kngp)
    if(ipri==1) then
       if(mp_kngp < 100000) then
          write(nfout,'(" !|| ista_kngp, iend_kngp = ",i5,",",i5,", mp_kngp = ",i5,", kngp = ",i12&
               & ," << m_Parallel_init_mpi_kngp >>")') ista_kngp,iend_kngp,mp_kngp,kngp
       else
          write(nfout,'(" !|| ista_kngp, iend_kngp = ",i10,",",i10,", mp_kngp = ",i10,", kngp = ",i12&
               & ," << m_Parallel_init_mpi_kngp >>")') ista_kngp,iend_kngp,mp_kngp,kngp
       end if
    else if(ipri>=2) then
       write(nfout,'(" << m_Parallel_init_mpi_kngp >>")')
       write(nfout,'(" !|| kngp = ",i12)') kngp
       write(nfout,'(" !|| -- npe, is_kngp(ista_kngp),  ie_kngp(iend_kngp) --")')
       do i = 0, npes-1
          write(nfout,'(" !|| ",i0,2i10)') i,is_kngp(i),ie_kngp(i)
       end do
    end if
  end subroutine m_Parallel_init_mpi_kngp

  subroutine m_Parallel_init_mpi_kngp_prev(kngp)
    integer, intent(in) :: kngp
    integer :: iwork, i
    integer, allocatable, dimension(:) :: is_kngp,ie_kngp,nel_kngp
    allocate(is_kngp(0:npes-1))
    allocate(ie_kngp(0:npes-1))
    allocate(nel_kngp(0:npes-1))
    iwork = ( kngp - 1 ) / npes + 1
    do i = 0, npes-1
       is_kngp(i) = min(i*iwork+1, kngp+1)
       ie_kngp(i) = min(is_kngp(i)+iwork-1, kngp)
       nel_kngp(i) = ie_kngp(i) - is_kngp(i) + 1
    enddo
    ista_kngp_prev = is_kngp(mype)
    iend_kngp_prev = ie_kngp(mype)
    np_kngp_prev   = nel_kngp(mype)
    deallocate(is_kngp)
    deallocate(ie_kngp)
    deallocate(nel_kngp)
  end subroutine m_Parallel_init_mpi_kngp_prev

  subroutine m_Parallel_init_mpi_cdfft(nfout,ipri,ggacmp_parallel)
    integer, intent(in)  :: nfout,ipri,ggacmp_parallel
    integer :: i, iprilevel

    if(ggacmp_parallel == ON) then
       call split_into_ggablock_cdfft(npes,nrank_ggacmp,npes_cdfft,nrest_cdfft)
    else
       nrank_ggacmp = 1
       npes_cdfft   = npes
       nrest_cdfft  = 0
    end if
    if(ipri >= 1) write(nfout,'(" !|| ggacmp_parallel = ",i3)') ggacmp_parallel
    if(ggacmp_parallel==ON) then
       iprilevel=1
    else
       iprilevel=2
    end if
    if(ipri >= iprilevel) write(nfout,'(" !|| nrank_ggacmp, npes_cdfft, nrest_cdfft",/ &
         & ," !||   = ",3i10," <<m_Parallel_init_mpi_cdfft>>")') &
         & nrank_ggacmp, npes_cdfft, nrest_cdfft
    myrank_ggacmp = mype/npes_cdfft
    myrank_cdfft  = mype - myrank_ggacmp*npes_cdfft
    if(ipri >= iprilevel) write(nfout,'(" !|| myrank_ggacmp, myrank_cdfft",/ &
         & ," !|| = ",2i10," <<m_Parallel_init_mpi_cdfft>>")') &
         & myrank_ggacmp, myrank_cdfft
    if(allocated(map_pe2ggacmp)) deallocate(map_pe2ggacmp)
    if(allocated(map_pe2cdfft)) deallocate(map_pe2cdfft)
    allocate(map_pe2ggacmp(0:npes-1))
    allocate(map_pe2cdfft(0:npes-1))
    do i = 0, npes-1
       map_pe2ggacmp(i) = i/npes_cdfft
       map_pe2cdfft(i) = i - map_pe2ggacmp(i)*npes_cdfft
    end do
    if(ipri >= iprilevel) write(nfout,'(" !|| map_pe2ggacmp",/,99( " !|| ",10i5,/))') &
         & (map_pe2ggacmp(i),i=0,npes-1)
    if(ipri >= iprilevel) write(nfout,'(" !|| map_pe2cdfft",/,99( " !|| ",10i5,/))') &
         & (map_pe2cdfft(i),i=0,npes-1)
  end subroutine m_Parallel_init_mpi_cdfft

  subroutine m_Parallel_dealloc_mpi_gga()
    if(allocated(map_ggacmp))  deallocate(map_ggacmp)
    if(allocated(nis_fftp))    deallocate(nis_fftp)
    if(allocated(nie_fftp))    deallocate(nie_fftp)
    if(allocated(nel_fftp))    deallocate(nel_fftp)
    if(allocated(idisp_fftp))  deallocate(idisp_fftp)
    if(allocated(nis_fftph))   deallocate(nis_fftph)
    if(allocated(nie_fftph))   deallocate(nie_fftph)
    if(allocated(nel_fftph))   deallocate(nel_fftph)
    if(allocated(idisp_fftph)) deallocate(idisp_fftph)
    if(allocated(nis_sfftp))   deallocate(nis_sfftp)
    if(allocated(nie_sfftp))   deallocate(nie_sfftp)
    if(allocated(nel_sfftp))   deallocate(nel_sfftp)
    if(allocated(idisp_sfftp)) deallocate(idisp_sfftp)
    if(allocated(nis_sfftph))  deallocate(nis_sfftph)
    if(allocated(nie_sfftph))  deallocate(nie_sfftph)
    if(allocated(nel_sfftph))  deallocate(nel_sfftph)
  end subroutine m_Parallel_dealloc_mpi_gga

  subroutine m_Parallel_init_mpi_gga(nfout,ipri,printable,nfftp,nfftps)
    integer, intent(in)  :: nfout,ipri,nfftp,nfftps
    logical, intent(in)  :: printable
    integer              :: nfftph, i
    logical              :: set_mapping_func
    integer              :: icolor, key, j, newpes, newmype
    logical, save        :: firstcall = .true.
                                                  __TIMER_SUB_START(1236)
                                                  __TIMER_SUB_STOP(1236)

    allocate(map_ggacmp(3)); map_ggacmp = 0
    allocate(nis_fftp(0:npes_cdfft-1)); nis_fftp = 0
    allocate(nie_fftp(0:npes_cdfft-1)); nie_fftp = 0
    allocate(nel_fftp(0:npes_cdfft-1)); nel_fftp = 0
    allocate(idisp_fftp(0:npes_cdfft-1)); idisp_fftp = 0
    allocate(nis_fftph(0:npes_cdfft-1));nis_fftph = 0
    allocate(nie_fftph(0:npes_cdfft-1)); nie_fftph = 0
    allocate(nel_fftph(0:npes_cdfft-1)); nel_fftph = 0
    allocate(idisp_fftph(0:npes_cdfft-1)); idisp_fftph = 0
    if(.not.allocated(mpi_cdfft_world))then
    allocate(mpi_cdfft_world(0:nrank_ggacmp-1)); mpi_cdfft_world = 0
    endif
    if(.not.allocated(mpi_ggacmp_cross_world)) then
    allocate(mpi_ggacmp_cross_world(0:npes_cdfft-1)); mpi_ggacmp_cross_world = 0
    endif

    if(nrank_ggacmp > 1) then
       do i = 1, 3
          map_ggacmp(i) = i-1
       end do
    else
       map_ggacmp(1:3) = 0
    end if

    set_mapping_func = .false.
    nfftph = nfftp/2
    if(ipri >= 1) write(nfout,'(" !|| nfftp, nfftph = ",2i10," <<m_Parallel_init_mpi_gga>>")') nfftp, nfftph
    call set_block_range(nfftph,npes_cdfft,nel_fftph,nis_fftph,nie_fftph,set_mapping_func)

    ista_fftph = nis_fftph(myrank_cdfft)
    iend_fftph = nie_fftph(myrank_cdfft)

    if(ipri >= 1) write(nfout,'(" !|| ista_fftph, iend_fftph = ",2i10 &
         & ," <<m_Parallel_init_mpi_gga>>")') ista_fftph, iend_fftph

    do i = 0, npes_cdfft -1
       nis_fftp(i) = nis_fftph(i)*2 - 1
       nie_fftp(i) = nie_fftph(i)*2
       nel_fftp(i) = nel_fftph(i)*2
       idisp_fftp(i) = nis_fftp(i) - 1
       idisp_fftph(i) = nis_fftph(i) - 1
    end do
    ista_fftp = nis_fftp(myrank_cdfft)
    iend_fftp = nie_fftp(myrank_cdfft)
    mp_fftp = maxval(nel_fftp)
    np_fftp = nel_fftp(myrank_cdfft)

    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| --- nis_fftp,nie_fftp,nel_fftp,idisp_fftp ---")')
       write(nfout,'(" !||  npes_cdfft = ",i12)') npes_cdfft
       write(nfout,'(" !|| ( npe_cdfft )",10i12)') (i,i=1,npes_cdfft)
       write(nfout,'(" !|| ( nis_fftp )",10i12)') (nis_fftp(i),i=0,npes_cdfft-1)
       write(nfout,'(" !|| ( nie_fftp )",10i12)') (nie_fftp(i),i=0,npes_cdfft-1)
       write(nfout,'(" !|| ( nel_fftp )",10i12)') (nel_fftp(i),i=0,npes_cdfft-1)
       write(nfout,'(" !|| (idisp_fftp)",10i12)') (idisp_fftp(i),i=0,npes_cdfft-1)
       write(nfout,'(" !|| ( ista_fftp,  iend_fftp )",2i12)') ista_fftp, iend_fftp
       write(nfout,'(" !|| --- nis_fftph,nie_fftph,nel_fftph ---")')
       write(nfout,'(" !|| ( npe_cdfft )",10i12)') (i,i=1,npes_cdfft)
       write(nfout,'(" !|| ( nis_fftph)",10i12)') (nis_fftph(i),i=0,npes_cdfft-1)
       write(nfout,'(" !|| ( nie_fftph)",10i12)') (nie_fftph(i),i=0,npes_cdfft-1)
       write(nfout,'(" !|| ( nel_fftph)",10i12)') (nel_fftph(i),i=0,npes_cdfft-1)
       write(nfout,'(" !|| ( ista_fftph, iend_fftph)",2i12)') ista_fftph, iend_fftph
    end if

    allocate(nis_sfftp(0:npes-1)); nis_sfftp = 0
    allocate(nie_sfftp(0:npes-1)); nie_sfftp = 0
    allocate(nel_sfftp(0:npes-1)); nel_sfftp = 0
    allocate(idisp_sfftp(0:npes-1)); idisp_sfftp = 0
    allocate(nis_sfftph(0:npes-1)); nis_sfftph = 0
    allocate(nie_sfftph(0:npes-1)); nie_sfftph = 0
    allocate(nel_sfftph(0:npes-1)); nel_sfftph = 0
    set_mapping_func = .false.
    if(ipri >= 1 .and. printable) write(nfout,'(" !|| nfftp = ",i10," <<m_Parallel_init_mpi_gga>>")') nfftp

!!$    call set_block_range(nfftp,npes,nel_sfftp,nis_sfftp,nie_sfftp,set_mapping_func)
!!$    ista_sfftp = nis_sfftp(mype)
!!$    iend_sfftp = nie_sfftp(mype)
!!$    idisp_sfftp(0) = 0
!!$    do i = 1, npes-1
!!$       idisp_sfftp(i) = idisp_sfftp(i-1) + nel_sfftp(i-1)
!!$    end do

    nfftph = nfftps/2
    call set_block_range(nfftph,npes,nel_sfftph,nis_sfftph,nie_sfftph,set_mapping_func)
    ista_sfftph = nis_sfftph(mype)
    iend_sfftph = nie_sfftph(mype)

    do i = 0, npes -1
       nis_sfftp(i) = nis_sfftph(i)*2 - 1
       nie_sfftp(i) = nie_sfftph(i)*2
       nel_sfftp(i) = nel_sfftph(i)*2
       idisp_sfftp(i) = nis_sfftp(i) - 1
    end do
    ista_sfftp = nis_sfftp(mype)
    iend_sfftp = nie_sfftp(mype)
    mp_sfftp = maxval(nel_sfftp)
    np_sfftp = nel_sfftp(mype)

    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| --- nis_sfftp,nie_sfftp,nel_sfftp,idisp_sfftp ---")')
       write(nfout,'(" !|| ( npe      )",10i12)') (i,i=1,npes)
       write(nfout,'(" !|| ( nis_sfftp )",10i12)') (nis_sfftp(i),i=0,npes-1)
       write(nfout,'(" !|| ( nie_sfftp )",10i12)') (nie_sfftp(i),i=0,npes-1)
       write(nfout,'(" !|| ( nel_sfftp )",10i12)') (nel_sfftp(i),i=0,npes-1)
       write(nfout,'(" !|| (idisp_sfftp)",10i12)') (idisp_sfftp(i),i=0,npes-1)
       write(nfout,'(" !|| ( ista_sfftp,  iend_sfftp )",2i12)') ista_sfftp, iend_sfftp
       write(nfout,'(" !|| --- nis_sfftph,nie_sfftph,nel_sfftph ---")')
       write(nfout,'(" !|| ( npe       )",10i12)') (i,i=1,npes)
       write(nfout,'(" !|| ( nis_sfftph )",10i12)') (nis_sfftph(i),i=0,npes-1)
       write(nfout,'(" !|| ( nie_sfftph )",10i12)') (nie_sfftph(i),i=0,npes-1)
       write(nfout,'(" !|| ( nel_sfftph )",10i12)') (nel_sfftph(i),i=0,npes-1)
       write(nfout,'(" !|| ( ista_sfftph, iend_sfftph)",2i12)') ista_sfftph, iend_sfftph
    end if

    do j = 0, nrank_ggacmp-1
       icolor = 0
       key = 0
       if(myrank_ggacmp == j) then
          icolor = 1
          key = myrank_cdfft
       end if
!!$       do i = 0, npes_cdfft-1
!!$          if(mype == i + npes_cdfft*j) then
!!$             icolor = 1
!!$             key = i
!!$          end if
!!$       end do
!       if(firstcall) call mpi_comm_split(MPI_CommGroup,icolor,key,mpi_cdfft_world(j),ierr)
       if(.not.firstcall) call mpi_comm_free(mpi_cdfft_world(j), ierr)
       call mpi_comm_split(MPI_CommGroup,icolor,key,mpi_cdfft_world(j),ierr)
       call mpi_comm_size(mpi_cdfft_world(j), newpes, ierr)
       call mpi_comm_rank(mpi_cdfft_world(j), newmype, ierr)
    end do

    do j = 0, npes_cdfft-1
       icolor = 0
       key = 0
       if(myrank_cdfft == j) then
          icolor = 1
          key = myrank_ggacmp
       end if
!!$       do i = 0, nrank_ggacmp-1
!!$          if(mype == i + nrank_ggacmp*j) then
!!$             icolor = 1
!!$             key = i
!!$          end if
!!$       end do
!       if(firstcall) call mpi_comm_split(MPI_CommGroup,icolor,key,mpi_ggacmp_cross_world(j),ierr)
       if(.not.firstcall) call mpi_comm_free(mpi_ggacmp_cross_world(j), ierr)
       call mpi_comm_split(MPI_CommGroup,icolor,key,mpi_ggacmp_cross_world(j),ierr)
       call mpi_comm_size(mpi_ggacmp_cross_world(j),newpes,ierr)
       call mpi_comm_rank(mpi_ggacmp_cross_world(j),newmype,ierr)
       if(ipri >= 2) then
! === DEBUG by tkato 2011/07/12 ================================================
!         write(nfout,'(" !|| mype, newmype, newpes, mpi_ggacmp_cross_world = ",4i8)') &
!              & mype, newmype,newpes,mpi_ggacmp_cross_world(j)
          write(nfout,'(" !|| mype, newmype, newpes, mpi_ggacmp_cross_world = ",4i12)') &
               & mype, newmype,newpes,mpi_ggacmp_cross_world(j)
! ==============================================================================
       end if
    end do
    firstcall = .false.
  end subroutine m_Parallel_init_mpi_gga


  subroutine m_Parallel_init_mpi_mix(nfout,ipri,printable,kgpm)
    integer, intent(in) :: nfout,ipri,kgpm
    logical, intent(in) :: printable
    integer :: iwork, i
    allocate(is_kgpm(0:npes-1))
    allocate(ie_kgpm(0:npes-1))
    allocate(nel_kgpm(0:npes-1))
    iwork = ( kgpm - 1 ) / npes + 1
    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| << init_mpi_mix >>")')
       write(nfout,'(" !|| kgpm = ",i12)') kgpm
       write(nfout,'(" !|| -- is_kgpm, ie_kgpm --")')
    end if
    do i = 0, npes-1
       is_kgpm(i) = min(i*iwork+1, kgpm+1)
       ie_kgpm(i) = min(is_kgpm(i)+iwork-1, kgpm)
       nel_kgpm(i) = ie_kgpm(i) - is_kgpm(i) + 1
       if(ipri >= 2 .and. printable) write(nfout,'(" !|| ",2i12)') is_kgpm(i),ie_kgpm(i)
    enddo
    ista_kgpm = is_kgpm(mype)
    iend_kgpm = ie_kgpm(mype)
    np_kgpm   = nel_kgpm(mype)
    mp_kgpm   = maxval(nel_kgpm)
    if(ipri == 1 .and. printable) then
       if(iend_kgpm <1000000) then
          write(nfout,'(" !|| ista_kgpm, iend_kgpm = ",i6,",",i6,", mp_kgpm = ",i6,", kgpm = ",i12 &
               & ,"  << m_Parallel_init_mpi_mix>>")') ista_kgpm,iend_kgpm,mp_kgpm,kgpm
       elseif(iend_kgpm <100000000)then
          write(nfout,'(" !|| ista_kgpm, iend_kgpm = ",i8,",",i8,", mp_kgpm = ",i8,", kgpm = ",i12 &
               & ,"  << m_Parallel_init_mpi_mix>>")') ista_kgpm,iend_kgpm,mp_kgpm,kgpm
       else
          write(nfout,'(" !|| ista_kgpm, iend_kgpm = ",i12,",",i12,", mp_kgpm = ",i12,", kgpm = ",i0 &
               & ,"  << m_Parallel_init_mpi_mix>>")') ista_kgpm,iend_kgpm,mp_kgpm,kgpm
       end if
    end if
  end subroutine m_Parallel_init_mpi_mix

  subroutine m_Parallel_init_mpi_snl(nfout,ipri,printable,nspin)
    integer, intent(in) :: nfout,ipri,nspin
    logical, intent(in) :: printable

    ista_snl = (ista_k + nspin - 1)/nspin
    iend_snl = iend_k/nspin
    if(ipri >= 1 .and. printable) then
       write(nfout,'(" !|| ista_snl, iend_snl = ",i5,",",i5," << m_Parallel_init_mpi_snl >>")') &
            & ista_snl,iend_snl
    end if
  end subroutine m_Parallel_init_mpi_snl

  subroutine m_Parallel_init_mpi_ffth(nfout,ipri,printable,nfft)
    integer, intent(in) :: nfout,ipri,nfft
    logical, intent(in) :: printable
    integer :: iwork, i, nffth
    nffth = nfft/2

    if ( .not. allocated(is_ffth) ) allocate(is_ffth(0:npes-1))
    if ( .not. allocated(ie_ffth) ) allocate(ie_ffth(0:npes-1))
    if ( .not. allocated(nel_ffth) ) allocate(nel_ffth(0:npes-1))

    iwork = ( nffth - 1 ) / npes + 1
    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| << init_mpi_ffth >>")')
       write(nfout,'(" !|| nffth = ",i12)') nffth
       write(nfout,'(" !|| -- is_ffth, ie_ffth --")')
    end if
    do i = 0, npes-1
       is_ffth(i) = min(i*iwork+1, nffth+1)
       ie_ffth(i) = min(is_ffth(i)+iwork-1, nffth)
       nel_ffth(i) = ie_ffth(i) - is_ffth(i) + 1
       if(ipri >= 2 .and. printable) write(nfout,'(" !|| ",2i12)') is_ffth(i),ie_ffth(i)
    enddo
    ista_ffth = is_ffth(mype)
    iend_ffth = ie_ffth(mype)
    np_ffth   = nel_ffth(mype)
    mp_ffth   = maxval(nel_ffth)
    if(ipri == 1 .and. printable) then
       if(iend_ffth <1000) then
          write(nfout,'(" !|| ista_ffth, iend_ffth = ",i3,",",i3,", mp_ffth = ",i3,", nffth = ",i5 &
               & ,"  << m_Parallel_init_mpi_ffth >>")') ista_ffth,iend_ffth,mp_ffth,nffth
       elseif(iend_ffth <100000)then
          write(nfout,'(" !|| ista_ffth, iend_ffth = ",i5,",",i5,", mp_ffth = ",i5,", ffth = ",i6 &
               & ,"  << m_Parallel_init_mpi_ffth >>")') ista_ffth,iend_ffth,mp_ffth,nffth
       else
          write(nfout,'(" !|| ista_ffth, iend_ffth = ",i0,",",i0,", mp_atm = ",i0,", natm = ",i0 &
               & ,"  << m_Parallel_init_mpi_ffth >>")') ista_ffth,iend_ffth,mp_ffth,nffth
       end if
    end if
  end subroutine m_Parallel_init_mpi_ffth


  subroutine m_Parallel_init_mpi_nn(nfout,ipri,printable,nn)
    integer, intent(in) :: nfout,ipri,nn
    logical, intent(in) :: printable
    integer :: iwork, i,ia

    allocate(is_nn(0:nrank_e-1))
    allocate(ie_nn(0:nrank_e-1))
    allocate(nel_nn(0:nrank_e-1))
    iwork = ( nn - 1 ) / nrank_e + 1
    do i = 0, nrank_e-1
       is_nn(i) = min(i*iwork+1, nn+1)
       ie_nn(i) = min(is_nn(i)+iwork-1, nn)
       nel_nn(i) = ie_nn(i) - is_nn(i) + 1
       if(ipri >= 2 .and. printable) write(nfout,'(" !|| ",2i12)') is_nn(i),ie_nn(i)
    enddo
    ista_nn = is_nn(myrank_e)
    iend_nn = ie_nn(myrank_e)
  end subroutine m_Parallel_init_mpi_nn

  subroutine m_Parallel_init_mpi_rspace_aug(nfout,ipri,printable,natm,nmesh_rs_aug)
    integer, intent(in) :: nfout,ipri,natm
    logical, intent(in) :: printable
    integer, dimension(natm), intent(in) :: nmesh_rs_aug
    integer :: iwork, i,ia

    if(allocated(is_rspace_aug))deallocate(is_rspace_aug)
    if(allocated(ie_rspace_aug))deallocate(ie_rspace_aug)
    if(allocated(nel_rspace_aug))deallocate(nel_rspace_aug)
    if(allocated(ista_rspace_aug_atm))deallocate(ista_rspace_aug_atm)
    if(allocated(iend_rspace_aug_atm))deallocate(iend_rspace_aug_atm)
    allocate(is_rspace_aug(0:npes-1))
    allocate(ie_rspace_aug(0:npes-1))
    allocate(nel_rspace_aug(0:npes-1))
    allocate(ista_rspace_aug_atm(natm));ista_rspace_aug_atm=0
    allocate(iend_rspace_aug_atm(natm));iend_rspace_aug_atm=0
    do ia=1,natm
    iwork = ( nmesh_rs_aug(ia) - 1 ) / npes + 1
    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| << init_mpi_rspace_aug >>")')
       write(nfout,'(" !|| nmesh_rs_aug = ",i12)') nmesh_rs_aug(ia)
       write(nfout,'(" !|| -- is_rspace_aug, ie_rspace_aug --")')
    end if
    do i = 0, npes-1
       is_rspace_aug(i) = min(i*iwork+1, nmesh_rs_aug(ia)+1)
       ie_rspace_aug(i) = min(is_rspace_aug(i)+iwork-1, nmesh_rs_aug(ia))
       nel_rspace_aug(i) = ie_rspace_aug(i) - is_rspace_aug(i) + 1
       if(ipri >= 2 .and. printable) write(nfout,'(" !|| ",2i12)') is_rspace_aug(i),ie_rspace_aug(i)
    enddo
    ista_rspace_aug_atm(ia) = is_rspace_aug(mype)
    iend_rspace_aug_atm(ia) = ie_rspace_aug(mype)
    if(ipri == 1 .and. printable) then
       if(iend_rspace_aug_atm(ia) <1000) then
          write(nfout,'(" !|| ista_rspace_aug, iend_rspace_aug = ",i3,",",i3,", mp_rspace_aug = ",i3,", nmesh_rs_aug_max = ",i5 &
          & ,"  << m_Parallel_init_mpi_rspace_aug >>")') &
          & ista_rspace_aug_atm(ia),iend_rspace_aug_atm(ia),mp_rspace_aug,nmesh_rs_aug(ia)
       else
          write(nfout,'(" !|| ista_rspace_aug, iend_rspace_aug = ",i0,",",i0,", mp_rspace_aug = ",i0,", nmesh_rs_aug_max = ",i0 &
          & ,"  << m_Parallel_init_mpi_atm >>")') ista_rspace_aug_atm(ia),iend_rspace_aug_atm(ia),mp_rspace_aug,nmesh_rs_aug(ia)
       end if
    end if
    enddo
  end subroutine m_Parallel_init_mpi_rspace_aug


  subroutine m_Parallel_init_mpi_nq(nfout,ipri,printable,nq)
    integer, intent(in) :: nfout, ipri,nq
    logical, intent(in) :: printable
    integer :: i,j,ip

    allocate(is_nq(0:npes-1));ista_nq=nq+1
    allocate(ie_nq(0:npes-1));iend_nq=0
    allocate(nel_nq(0:npes-1));nel_nq=0
    allocate(map_nq(nq));map_nq=0
    allocate(map_z_nq(nq));map_z_nq=0

    call set_block_range(nq,npes,nel_nq,is_nq,ie_nq,.true.,map_nq)
    ista_nq = is_nq(mype)
    iend_nq = ie_nq(mype)
    np_nq = nel_nq(mype)
    mp_nq = maxval(nel_nq)
    j = 0
    do ip = 1, npes
       do i = 1, nel_nq(ip-1)
          j = j + 1
          map_z_nq(j) = i
       end do
    end do
  end subroutine m_Parallel_init_mpi_nq


  subroutine m_Parallel_init_mpi_atm(nfout,ipri,printable,natm)
    integer, intent(in) :: nfout,ipri,natm
    logical, intent(in) :: printable
    integer :: iwork, i

    allocate(is_atm(0:npes-1))
    allocate(ie_atm(0:npes-1))
    allocate(nel_atm(0:npes-1))
    iwork = ( natm - 1 ) / npes + 1
    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| << init_mpi_atm >>")')
       write(nfout,'(" !|| natm = ",i12)') natm
       write(nfout,'(" !|| -- is_natm, ie_natm --")')
    end if
    do i = 0, npes-1
       is_atm(i) = min(i*iwork+1, natm+1)
       ie_atm(i) = min(is_atm(i)+iwork-1, natm)
       nel_atm(i) = ie_atm(i) - is_atm(i) + 1
       if(ipri >= 2 .and. printable) write(nfout,'(" !|| ",2i12)') is_atm(i),ie_atm(i)
    enddo
    ista_atm = is_atm(mype)
    iend_atm = ie_atm(mype)
    np_atm   = nel_atm(mype)
    mp_atm   = maxval(nel_atm)
! ====================================================================
    write(nfout,*) '===== m_Parallel_init_mpi_atm ====='
    write(nfout,*) 'npes: ', npes, ' = nrank_g : ',nrank_g
    write(nfout,*) 'mype: ', mype
    write(nfout,*) 'ista_atm, iend_atm, np_natm, mp_atm: ', ista_atm, iend_atm, np_atm, mp_atm
    write(nfout,*) '----- natm Decomposed as... -----'
    write(nfout,*) 'natm:  ', natm
    write(nfout,*) ' i, is_atm, ie_atm, nel_atm'
    do i = 0, npes-1
      write(nfout,*) i, is_atm(i), ie_atm(i), nel_atm(i)
    enddo
! =====================================================================
    if(ipri == 1 .and. printable) then
       if(iend_atm <1000) then
          write(nfout,'(" !|| ista_atm, iend_atm = ",i3,",",i3,", mp_atm = ",i3,", natm = ",i5 &
               & ,"  << m_Parallel_init_mpi_atm >>")') ista_atm,iend_atm,mp_atm,natm
       elseif(iend_kgpm <100000)then
          write(nfout,'(" !|| ista_atm, iend_atm = ",i5,",",i5,", mp_atm = ",i5,", natm = ",i6 &
               & ,"  << m_Parallel_init_mpi_atm >>")') ista_atm,iend_atm,mp_atm,natm
       else
          write(nfout,'(" !|| ista_atm, iend_atm = ",i0,",",i0,", mp_atm = ",i0,", natm = ",i0 &
               & ,"  << m_Parallel_init_mpi_atm >>")') ista_atm,iend_atm,mp_atm,natm
       end if
    end if
  end subroutine m_Parallel_init_mpi_atm

  subroutine m_Parallel_init_mpi_atm2(nfout,ipri,printable,natm2)
    integer, intent(in) :: nfout,ipri,natm2
    logical, intent(in) :: printable
    integer :: iwork, i
    if(allocated(is_atm2)) deallocate(is_atm2)
    if(allocated(ie_atm2)) deallocate(ie_atm2)
    if(allocated(nel_atm2)) deallocate(nel_atm2)
    allocate(is_atm2(0:npes-1))
    allocate(ie_atm2(0:npes-1))
    allocate(nel_atm2(0:npes-1))
    iwork = ( natm2 - 1 ) / npes + 1
    if(ipri >= 2 .and. printable) then
       write(nfout,'(" !|| << init_mpi_atm2 >>")')
       write(nfout,'(" !|| natm2 = ",i12)') natm2
       write(nfout,'(" !|| -- is_natm2, ie_natm2 --")')
    end if
    do i = 0, npes-1
       is_atm2(i) = min(i*iwork+1, natm2+1)
       ie_atm2(i) = min(is_atm2(i)+iwork-1, natm2)
       nel_atm2(i) = ie_atm2(i) - is_atm2(i) + 1
       if(ipri >= 2 .and. printable) write(nfout,'(" !|| ",2i12)') is_atm2(i),ie_atm2(i)
    enddo
    ista_atm2 = is_atm2(mype)
    iend_atm2 = ie_atm2(mype)
    np_atm2   = nel_atm2(mype)
    mp_atm2   = maxval(nel_atm2)

    if(ipri == 1 .and. printable) then
       if(iend_atm2 <1000) then
          write(nfout,'(" !|| ista_atm2, iend_atm2 = ",i3,",",i3,", mp_atm2 = ",i3,", natm2 = ",i5 &
               & ,"  << m_Parallel_init_mpi_atm2 >>")') ista_atm2,iend_atm2,mp_atm2,natm2
       elseif(iend_kgpm <100000)then
          write(nfout,'(" !|| ista_atm2, iend_atm2 = ",i5,",",i5,", mp_atm2 = ",i5,", natm2 = ",i6 &
               & ,"  << m_Parallel_init_mpi_atm2 >>")') ista_atm2,iend_atm2,mp_atm2,natm2
       else
          write(nfout,'(" !|| ista_atm2, iend_atm2 = ",i0,",",i0,", mp_atm2 = ",i0,", natm2 = ",i0 &
               & ,"  << m_Parallel_init_mpi_atm2 >>")') ista_atm2,iend_atm2,mp_atm2,natm2
       end if
    end if
  end subroutine m_Parallel_init_mpi_atm2

!!$#ifdef TRANSPOSE
  subroutine m_Parallel_init_mpi_nlmta(nfout,ipri,nlmta,nlmt,natm,lmta,ntyp,ilmt,ityp)
    integer, intent(in) :: nfout, ipri, nlmta, nlmt,natm,ntyp
!!$ASASASASAS
!!$    integer, intent(in), dimension(nlmt,ntyp) :: lmta
    integer, intent(in), dimension(nlmt,natm) :: lmta
!!$ASASASASAS
    integer, intent(in), dimension(ntyp)      :: ilmt
    integer, intent(in), dimension(natm)      :: ityp

    integer :: ip, ia, it, it0, n1, i, nsum
!!$ASASASASAS
    integer :: itmp
!!$ASASASASAS
    real(kind=8), allocatable, dimension(:) :: np
    real(kind=8) :: n_l, n_r
    allocate(np(0:nrank_e-1)); np = 0.d0


    if(allocated(nis_fs)) deallocate(nis_fs)
    if(allocated(nie_fs)) deallocate(nie_fs)
    if(allocated(nel_fs)) deallocate(nel_fs)
    if(allocated(nis_fs_atm)) deallocate(nis_fs_atm)
    if(allocated(nie_fs_atm)) deallocate(nie_fs_atm)
    if(allocated(nel_fs_atm)) deallocate(nel_fs_atm)
    allocate(nis_fs(0:nrank_e-1))
    allocate(nie_fs(0:nrank_e-1))
    allocate(nel_fs(0:nrank_e-1))
    allocate(nis_fs_atm(0:nrank_e-1))
    allocate(nie_fs_atm(0:nrank_e-1))
    allocate(nel_fs_atm(0:nrank_e-1))

    do ip = 1, nrank_e
       np(ip-1) = dble(nlmta)/nrank_e *ip
       if(ipri >= 2) write(nfout,'(" !|| ip = ",i5," np = ",f12.4)') ip,np(ip-1)
    end do

    ip = 0
    do ia = 1, natm
       it = ityp(ia)
!!$ASASASASAS
       itmp=ilmt(it)
       if(itmp==0)cycle
!!$ASASASASAS
       n1 = lmta(ilmt(it),ia)
       if(abs(n1-np(ip)) < 1.d-10) then
          ip = ip + 1
          nie_fs(ip-1) = n1
          nie_fs_atm(ip-1) = ia
       else if(n1 > np(ip)) then
          ip = ip+1
          if(ia == 1) then
             n_l = np(ip-1)+1
          else if(ia >= 2) then
             it0 = ityp(ia-1)
             itmp=ilmt(it0)
!!$ASASASASAS
             if(itmp/=0) n_l = np(ip-1)-lmta(ilmt(it0),ia-1)
!!$             n_l = np(ip-1)-lmta(ilmt(it0),ia-1)
!!$ASASASASAS
          end if
          n_r = n1 - np(ip-1)
          if(n_l >= n_r) then
             nie_fs(ip-1) = n1
             nie_fs_atm(ip-1) = ia
          else
             nie_fs(ip-1) = np(ip-1) - n_l
!!$ASASASASAS
!!$             if(nie_fs(ip-1) < 0) nie_fs(ip-1) = n1
             if(nie_fs(ip-1) < 0) nie_fs(ip-1) = 0
!!$ASASASASAS
             nie_fs_atm(ip-1) = ia-1
          end if
       end if
    end do
    do i = ip+1, nrank_e
       nie_fs(i-1) = nlmta
       nie_fs_atm(i-1) = natm
    end do
    nis_fs_atm(0) = 1
    nis_fs(0) = 1
    nel_fs(0) = nie_fs(0)
    nel_fs_atm(0) = nie_fs_atm(0)
    do ip = 1, nrank_e-1
       nis_fs_atm(ip) = nie_fs_atm(ip-1)+1
       nis_fs(ip)  = nie_fs(ip-1)+1
       nel_fs(ip) = nie_fs(ip) - nie_fs(ip-1)
       nel_fs_atm(ip) = nie_fs_atm(ip) - nie_fs_atm(ip-1)
    end do

    ista_fs = nis_fs(myrank_e)
    iend_fs = nie_fs(myrank_e)
    np_fs   = nel_fs(myrank_e)
    mp_fs   = maxval(nel_fs)
    ista_fs_atm = nis_fs_atm(myrank_e)
    iend_fs_atm = nie_fs_atm(myrank_e)
    np_fs_atm = nel_fs_atm(myrank_e)
    mp_fs_atm = maxval(nel_fs_atm)

    if(ipri >= 2) then
       nsum = 0
       write(nfout,'(" !||",9x,"ip",4x,"nis_fs",4x,"nie_fs",7x,"nel",11x &
         & ,"atm(nis,nie)",5x,"lmta")')
       do ip = 1, nrank_e
          it = ityp(nie_fs_atm(ip-1))
          if(nie_fs_atm(ip-1) > 0 .and. nie_fs_atm(ip-1) <= natm ) then
             write(nfout,'(" !|| ",4i10,2x,3i10)') &
                  &   ip,nis_fs(ip-1),nie_fs(ip-1),nel_fs(ip-1), nis_fs_atm(ip-1), nie_fs_atm(ip-1) &
                  & , lmta(ilmt(it),nie_fs_atm(ip-1))
          else
             write(nfout,'(" !|| ",4i10,2x,2i10)') &
                  &   ip,nis_fs(ip-1),nie_fs(ip-1),nel_fs(ip-1), nis_fs_atm(ip-1), nie_fs_atm(ip-1)
          end if

          nsum = nsum + nel_fs(ip-1)
       end do
       write(nfout,'(" !|| nsum = ", i5)') nsum

       write(nfout,'(" !|| --- nis_fs,nie_fs,nel_fs ---")')
       write(nfout,'(" !|| ( rank_e  )",15i7)')(i,i=1,nrank_e)
       write(nfout,'(" !|| ( nis_fs  )",15i7)')(nis_fs(i),i=0,nrank_e-1)
       write(nfout,'(" !|| ( nie_fs  )",15i7)')(nie_fs(i),i=0,nrank_e-1)
       write(nfout,'(" !|| ( nel_fs  )",15i7)')(nel_fs(i),i=0,nrank_e-1)
       write(nfout,'(" !|| myrank_e = ",i3 &
            & ," nlmta,ista_fs,iend_fs,np_fs,mp_fs = ",/,10x,5i7)') &
            & myrank_e,nlmta,ista_fs,iend_fs,np_fs,mp_fs

       write(nfout,'(" !|| --- nis_fs_atm,nie_fs_atm,nel_fs_atm ---")')
       write(nfout,'(" !|| ( rank_e  )",15i7)')(i,i=1,nrank_e)
       write(nfout,'(" !|| ( nis_fs_atm  )",15i7)')(nis_fs_atm(i),i=0,nrank_e-1)
       write(nfout,'(" !|| ( nie_fs_atm  )",15i7)')(nie_fs_atm(i),i=0,nrank_e-1)
       write(nfout,'(" !|| ( nel_fs_atm  )",15i7)')(nel_fs_atm(i),i=0,nrank_e-1)
       write(nfout,'(" !|| myrank_e = ",i3 &
            & ," nlmta,ista_fs_atm,iend_fs_atm,np_fs_atm,mp_fs_atm = ",/,10x,5i7)') &
            & myrank_e,nlmta,ista_fs_atm,iend_fs_atm,np_fs_atm,mp_fs_atm

    end if

    deallocate(np)
  end subroutine m_Parallel_init_mpi_nlmta
!!$#endif

  subroutine split_into_ggablock_cdfft(npes,nrank_ggacmp,npes_cdfft,nrest_cdfft)
    integer, intent(in)  :: npes
    integer, intent(out) :: nrank_ggacmp, npes_cdfft, nrest_cdfft

#ifdef _NO_GGA_SPLIT_
    nrank_ggacmp = 1
    npes_cdfft  = npes
    nrest_cdfft = 0
#else
    if(npes <= 2) then
       nrank_ggacmp = 1
       npes_cdfft  = npes
       nrest_cdfft = 0
    else
       nrank_ggacmp = 3
       npes_cdfft  = npes/nrank_ggacmp
       nrest_cdfft = npes - npes_cdfft*nrank_ggacmp
    end if
#endif
  end subroutine split_into_ggablock_cdfft

  subroutine m_Parallel_dealloc_mpi_nlmta()
    if(allocated(nis_fs)) deallocate(nis_fs)
    if(allocated(nie_fs)) deallocate(nie_fs)
    if(allocated(nel_fs)) deallocate(nel_fs)
    if(allocated(nis_fs_atm)) deallocate(nis_fs_atm)
    if(allocated(nie_fs_atm)) deallocate(nie_fs_atm)
    if(allocated(nel_fs_atm)) deallocate(nel_fs_atm)
  end subroutine m_Parallel_dealloc_mpi_nlmta

  subroutine m_Parallel_dealloc_mpi_exx()
    if(allocated(is_kngp_exx)) deallocate(is_kngp_exx)
    if(allocated(ie_kngp_exx)) deallocate(ie_kngp_exx)
    if(allocated(nel_kngp_exx)) deallocate(nel_kngp_exx)
  end subroutine m_Parallel_dealloc_mpi_exx

  subroutine m_Parallel_cp_g1k(kv3)
    integer, intent(in) :: kv3
    ista_kngp_prev = ista_kngp
    iend_kngp_prev = iend_kngp
  end subroutine m_Parallel_cp_g1k

  subroutine m_Parallel_dealloc(neb_mode)
    logical, intent(in), optional :: neb_mode
    logical :: neb
    neb = .false.
    if(present(neb_mode)) neb = neb_mode

    if(allocated(is_kngp)) deallocate(is_kngp)
    if(allocated(ie_kngp)) deallocate(ie_kngp)
    if(allocated(nel_kngp)) deallocate(nel_kngp)

!    if(.not.neb)then
    if(allocated(map_e)) deallocate(map_e)
    if(allocated(map_z)) deallocate(map_z)
    if(allocated(map_k)) deallocate(map_k)
    if(allocated(map_ek)) deallocate(map_ek)
    if(allocated(map_eks)) deallocate(map_eks)

! <---- ASMS 2023/08/30 ---
    if(allocated(map_s)) deallocate(map_s)
    if(allocated(nis_spin)) deallocate(nis_spin)
    if(allocated(nie_spin)) deallocate(nie_spin)
    if(allocated(nel_spin)) deallocate(nel_spin)
! ---- ASMS 2023/08/30 --->

!    if(allocated(mpi_k_world)) deallocate(mpi_k_world)
    if(allocated(nis_e)) deallocate(nis_e)
    if(allocated(nie_e)) deallocate(nie_e)
    if(allocated(nel_e)) deallocate(nel_e)
    if(allocated(idisp_e)) deallocate(idisp_e)
#ifdef NEC_TUNE_SOFT
    if(allocated(ista_e_smp)) deallocate(ista_e_smp)
    if(allocated(iend_e_smp)) deallocate(iend_e_smp)
#endif
    if(allocated(nis_k)) deallocate(nis_k)
    if(allocated(nie_k)) deallocate(nie_k)
    if(allocated(nel_k)) deallocate(nel_k)

    if(allocated(ista_g1k)) deallocate(ista_g1k)
    if(allocated(iend_g1k)) deallocate(iend_g1k)
    if(allocated(np_g1k)) deallocate(np_g1k)
    if(allocated(np_g1k_prev)) deallocate(np_g1k_prev)
    if(allocated(mp_g1k)) deallocate(mp_g1k)

    if(allocated(nis_fftp)) deallocate(nis_fftp)
    if(allocated(nie_fftp)) deallocate(nie_fftp)
    if(allocated(nel_fftp)) deallocate(nel_fftp)
    if(allocated(idisp_fftp)) deallocate(idisp_fftp)
    if(allocated(nis_fftph)) deallocate(nis_fftph)
    if(allocated(nie_fftph)) deallocate(nie_fftph)
    if(allocated(nel_fftph)) deallocate(nel_fftph)

    if(allocated(is_kgpm)) deallocate(is_kgpm)
    if(allocated(ie_kgpm)) deallocate(ie_kgpm)
    if(allocated(nel_kgpm)) deallocate(nel_kgpm)

    if(allocated(is_atm)) deallocate(is_atm)
    if(allocated(ie_atm)) deallocate(ie_atm)
    if(allocated(nel_atm)) deallocate(nel_atm)

    if(allocated(is_atm2)) deallocate(is_atm2)
    if(allocated(ie_atm2)) deallocate(ie_atm2)
    if(allocated(nel_atm2)) deallocate(nel_atm2)
!!$#ifndef TRANSPOSE
!!$    endif
!!$#endif
!!$#ifdef TRANSPOSE
    if(allocated(nis_g1)) deallocate(nis_g1)
    if(allocated(nie_g1)) deallocate(nie_g1)
    if(allocated(nel_g1)) deallocate(nel_g1)
    if(.not.neb)then
    if(allocated(nis_fs)) deallocate(nis_fs)
    if(allocated(nie_fs)) deallocate(nie_fs)
    if(allocated(nel_fs)) deallocate(nel_fs)
    if(allocated(nis_fs_atm)) deallocate(nis_fs_atm)
    if(allocated(nie_fs_atm)) deallocate(nie_fs_atm)
    if(allocated(nel_fs_atm)) deallocate(nel_fs_atm)
    endif

    if(allocated(nis_g1k)) deallocate(nis_g1k)
    if(allocated(nie_g1k)) deallocate(nie_g1k)
    if(allocated(nel_g1k)) deallocate(nel_g1k)
!!$#endif

    if(allocated(map_pe2ggacmp)) deallocate(map_pe2ggacmp)
    if(allocated(map_pe2cdfft)) deallocate(map_pe2cdfft)

!    if(.not.neb)then
    if(allocated(map_ggacmp)) deallocate(map_ggacmp)
    if(allocated(idisp_fftph)) deallocate(idisp_fftph)
!    if(allocated(mpi_cdfft_world)) deallocate(mpi_cdfft_world)
    !if(allocated(mpi_ggacmp_cross_world)) deallocate(mpi_ggacmp_cross_world)

    if(allocated(nis_sfftp)) deallocate(nis_sfftp)
    if(allocated(nie_sfftp)) deallocate(nie_sfftp)
    if(allocated(nel_sfftp)) deallocate(nel_sfftp)
    if(allocated(idisp_sfftp)) deallocate(idisp_sfftp)
    if(allocated(nis_sfftph)) deallocate(nis_sfftph)
    if(allocated(nie_sfftph)) deallocate(nie_sfftph)
    if(allocated(nel_sfftph)) deallocate(nel_sfftph)

!    if(allocated(mpi_nbmx_world)) deallocate(mpi_nbmx_world)
    if(allocated(nel_nbmx)) deallocate(nel_nbmx)
    if(allocated(nis_nbmx)) deallocate(nis_nbmx)
    if(allocated(nie_nbmx)) deallocate(nie_nbmx)
    if(allocated(idisp_nbmx)) deallocate(idisp_nbmx)

    if(allocated(nel_nbmx_k)) deallocate(nel_nbmx_k)
    if(allocated(nis_nbmx_k)) deallocate(nis_nbmx_k)
    if(allocated(nie_nbmx_k)) deallocate(nie_nbmx_k)
    if(allocated(idisp_nbmx_k)) deallocate(idisp_nbmx_k)
!    if(allocated(mpi_nbmx_world_k)) deallocate(mpi_nbmx_world_k)
    if(allocated(nel_kg1_k)) deallocate(nel_kg1_k)
    if(allocated(nis_kg1_k)) deallocate(nis_kg1_k)
    if(allocated(nie_kg1_k)) deallocate(nie_kg1_k)
    if(allocated(idisp_kg1_k)) deallocate(idisp_kg1_k)
!    endif
  end subroutine m_Parallel_dealloc

 logical function m_Parallel_resolve_conf_para()
   integer :: i,nn
   integer :: n1
   integer iargc,narg
   character(100) :: q1
   character(100) arg
   nrank_conf=-1
   m_Parallel_resolve_conf_para = .false.
   narg = iargc()
   n1 = -1
   do nn=1,narg
      call getarg(nn,q1)

!!!!!!!!!!!!!!!!! modified by mizouchi@adv 2003.03.20 !!!!!!!!!!!!!!
      q1 = trim(adjustl(q1))
!!!!!!!!!!!!!!!!! modified by mizouchi@adv 2003.03.20 !!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!! modified by mizouchi@adv 2003.03.20 !!!!!!!!!!!!!!
!!          if(index(trim(adjustl(q1)),"ne").ne.0) then

      if(q1(1:1) == "n".and. (q1(2:2) == "r".or.q1(2:2) == "c" .or. q1(2:2) == "i")) then
!!!!!!!!!!!!!!!!! modified by mizouchi@adv 2003.03.20 !!!!!!!!!!!!!!
         if((q1(3:3) == "=".or.q1(3:3) == ":") .and. len_trim(q1(4:)).ne.0) then
             q1 = q1(4:)
             if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                read(q1,*) n1
             else
                call phase_error_with_msg(6,"wrong nr",__LINE__,__FILE__)
             end if
         else if(q1(3:3) == "=".or. q1(3:3) == ":") then
             call getarg(nn+1,q1)
             if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                read(q1,*) n1
             else
                call phase_error_with_msg(6,"wrong nr",__LINE__,__FILE__)
             end if
         else
              call getarg(nn+1,q1)
              if(q1(1:1) == ":".and. len_trim(q1(2:)).ne.0) then
                  q1 = q1(2:)
                  if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                     read(q1,*) n1
                  else
                     call phase_error_with_msg(6,"wrong nr",__LINE__,__FILE__)
                  end if
              else if((q1(1:1) == "=".or.q1(1:1) == ":").and. len_trim(q1(2:)).eq.0) then
                  call getarg(nn+2,q1)
                  if(iachar(q1(1:1)).ge.iachar("0").and.iachar(q1(1:1)).le.iachar("9")) then
                     read(q1,*) n1
                  else
                     call phase_error_with_msg(6,"wrong nr",__LINE__,__FILE__)
                  end if
               else
                  call phase_error_with_msg(6,"wrong nr",__LINE__,__FILE__)
               end if
         end if
      end if
   enddo

   nrank_conf = n1

   if(nrank_conf>=1) m_Parallel_resolve_conf_para = .true.
   conf_para = m_Parallel_resolve_conf_para
   if(conf_para)then
     call config_parallel_init()
   endif

   contains

   subroutine config_parallel_init()
     integer :: i,j
     integer npestmp,mypetmp,icolor,ikey,mpi_err
     integer mype_e, npes_e
     integer, allocatable ::  new_comm_world(:)
     logical, save        ::  firstcall=.true.
     call mpi_init(mpi_err)
     call mpi_comm_size(mpi_comm_world,npestmp,mpi_err)
     call mpi_comm_rank(mpi_comm_world,mypetmp,mpi_err)
     allocate(new_comm_world(0:nrank_conf-1))
     do i=0,nrank_conf-1
       icolor = 0
       ikey = 0
       do j=0, npestmp/nrank_conf-1
         if(mypetmp == j + (npestmp/nrank_conf)*i ) then
           icolor = 1
           ikey = i
         end if
       end do
!       if(firstcall) call mpi_comm_split( mpi_comm_world, icolor, ikey, new_comm_world(i), mpi_err )
       if(.not.firstcall) call mpi_comm_free(new_comm_world(i), mpi_err)
       call mpi_comm_split( mpi_comm_world, icolor, ikey, new_comm_world(i), mpi_err )
     end do

     mype_conf = mypetmp/(npestmp/nrank_conf)
     MPI_CommGroup = new_comm_world(mype_conf)
     firstcall = .false.
   end subroutine config_parallel_init

 end function m_Parallel_resolve_conf_para

  subroutine m_Parallel_init_mpi_nvale(nfout,ipri,printable,nvale)
    integer, intent(in) :: nfout, ipri,nvale
    logical, intent(in) :: printable
    integer :: i,j,ip
    nrank_nvale = nrank_e
    myrank_nvale = myrank_e
    if(mpi_nvale_enabled) call m_Parallel_dealloc_mpi_nvale()
    call m_Parallel_alloc_mpi_nvale(nvale)
    call set_block_range(nvale,nrank_nvale,nel_nvale,nis_nvale,nie_nvale,.true.,map_nvale)
    mpi_nvale_enabled = .true.
    if(ipri>=2) call wd_val_range_3D()
    ista_nvale = nis_nvale(myrank_nvale)
    iend_nvale = nie_nvale(myrank_nvale)
    np_nvale = nel_nvale(myrank_nvale)
    mp_nvale = maxval(nel_nvale)

    j = 0
    do ip = 1, nrank_nvale
       do i = 1, nel_nvale(ip-1)
          j = j + 1
          map_z_nvale(j) = i
       end do
    end do

    contains

    subroutine wd_val_range_3D()
      if(printable) then
         write(nfout,'(" !|| nrank_nvale")')
         write(nfout,'(" !||    i    : ",20i4)')(i,i=0,nrank_nvale-1)
         write(nfout,'(" !|| ista_nvale   : ",20i4)')(nis_nvale(i),i=0,nrank_nvale-1)
         write(nfout,'(" !|| iend_nvale   : ",20i4)')(nie_nvale(i),i=0,nrank_nvale-1)
         write(nfout,'(" !|| nel_nvale    : ",20i4)')(nel_nvale(i),i=0,nrank_nvale-1)
      end if
    end subroutine wd_val_range_3D
  end subroutine m_Parallel_init_mpi_nvale

  subroutine m_Parallel_alloc_mpi_nvale(nvale)
    integer, intent(in) :: nvale
    allocate(nis_nvale(0:nrank_nvale-1));ista_nvale=nvale+1
    allocate(nie_nvale(0:nrank_nvale-1));iend_nvale=0
    allocate(nel_nvale(0:nrank_nvale-1));nel_nvale=0
    allocate(map_nvale(nvale));map_nvale=0
    allocate(map_z_nvale(nvale));map_z_nvale=0
  end subroutine m_Parallel_alloc_mpi_nvale

  subroutine m_Parallel_dealloc_mpi_nvale()
    if(allocated(nis_nvale)) deallocate(nis_nvale)
    if(allocated(nie_nvale)) deallocate(nie_nvale)
    if(allocated(nel_nvale)) deallocate(nel_nvale)
    if(allocated(map_nvale)) deallocate(map_nvale)
    if(allocated(map_z_nvale)) deallocate(map_z_nvale)
    mpi_nvale_enabled = .false.
  end subroutine m_Parallel_dealloc_mpi_nvale

  subroutine m_Parallel_init_mpi_kv3_ek(nfout,ipri,printable,kv3_ek,nspin)
    integer, intent(in) :: nfout, ipri,kv3_ek,nspin
    logical, intent(in) :: printable
    integer :: i,j,ip
    integer :: siz
    siz = kv3_ek/nspin
    nrank_kv3_ek = nrank_k
    myrank_kv3_ek = myrank_k
    call m_Parallel_alloc_mpi_kv3_ek(siz)
    call set_block_range(siz,nrank_kv3_ek,nel_kv3_ek,nis_kv3_ek,nie_kv3_ek,.true.,map_kv3_ek)
    ista_kv3_ek = nis_kv3_ek(myrank_kv3_ek)
    iend_kv3_ek = nie_kv3_ek(myrank_kv3_ek)
    np_kv3_ek = nel_kv3_ek(myrank_kv3_ek)
    mp_kv3_ek = maxval(nel_kv3_ek)

    j = 0
    do ip = 1, nrank_kv3_ek
       do i = 1, nel_kv3_ek(ip-1)
          j = j + 1
          map_z_kv3_ek(j) = i
       end do
    end do
  end subroutine m_Parallel_init_mpi_kv3_ek

  subroutine m_Parallel_alloc_mpi_kv3_ek(kv3_ek)
    integer, intent(in) :: kv3_ek
    allocate(nis_kv3_ek(0:nrank_kv3_ek-1));ista_kv3_ek=kv3_ek+1
    allocate(nie_kv3_ek(0:nrank_kv3_ek-1));iend_kv3_ek=0
    allocate(nel_kv3_ek(0:nrank_kv3_ek-1));nel_kv3_ek=0
    allocate(map_kv3_ek(kv3_ek));map_kv3_ek=0
    allocate(map_z_kv3_ek(kv3_ek));map_z_kv3_ek=0
  end subroutine m_Parallel_alloc_mpi_kv3_ek


  subroutine m_Parallel_resolve_decomp(nfout,kv3,done_something)
    integer, intent(in) :: nfout,kv3
    logical, intent(out) :: done_something
    integer :: i,neng,icount,itarget,nk,ne
    integer, allocatable, dimension(:) :: cands1,cands2
    real(kind=8) :: closest,fac
    done_something = .false.
    if(read_from_args) return
    done_something = .true.
    neng = 1;nk=npes
    if (npes>kv3) then
       do i=kv3,1,-1
          if(mod(npes,i)==0 .and. mod(kv3,i)==0) then
            neng = npes/i
            nk = npes/neng
            exit
          endif
       enddo
    else if(npes<kv3) then
       do i=npes,1,-1
          if(mod(npes,i)==0 .and. mod(kv3,i)==0) then
            neng = npes/i
            nk = npes/neng
            exit
          endif
       enddo
    endif
    ne = neng

    nrank_k = nk
    nrank_e = ne
    if(mype==0) then
       !!$ if(mype == 0) print '(" nrank_e = ", i3)', nrank_e
       write(nfout,'(" nrank_e = ",i3)') nrank_e
       !!$ if(mype == 0) print '(" nrank_k = ", i3)', nrank_k
       write(nfout,'(" nrank_k = ",i3)') nrank_k
    endif
  end subroutine m_Parallel_resolve_decomp

end module m_Parallelization
