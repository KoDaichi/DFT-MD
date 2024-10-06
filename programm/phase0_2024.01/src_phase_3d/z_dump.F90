#if 0

#define ORIGINAL

!===============================================================================
  subroutine dump_zaj(zinn)

  use m_PlaneWaveBasisSet,  only : kg1, iba
  use m_Control_Parameters, only : nspin, kimg, neg, meg
  use m_Parallelization,    only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype &
       &                         , nrank_e,myrank_e,map_e,ista_e,iend_e,istep_e,idisp_e &
       &                         , map_z,np_e,mpi_k_world,myrank_k,map_k,ista_k,iend_k &
       &                         , nis_e, nie_e
  use m_IterationNumbers,   only : iteration_electronic, nk_in_the_process, iteration, nkgroup

  integer :: ierr,err

  real(kind=8),intent(in) :: zinn(kg1,np_e,ista_k:iend_k,kimg)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, is, js, nn, ifrom
  real(kind=8), allocatable, dimension(:,:,:,:) :: t_wk1, t_wk2
  character(len=4) :: sort
  logical :: existence
  character(len=3) :: c_num
  character(len=2) :: c_img
  character(len=80) :: fname
  character(len=8) :: psuffix = 'dump_zaj'
#ifdef ORIGINAL
  character(len=5)  :: suffix = '.dat3'
#else
  character(len=5)  :: suffix = '.dat1'
#endif

  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

   allocate( t_wk1(kg1,neg,ista_k:iend_k,kimg) )
   allocate( t_wk2(kg1,neg,ista_k:iend_k,kimg) )

   t_wk1   = 0.0d0
   t_wk2   = 0.0d0

   do l = 1, kimg
   do k = ista_k, iend_k
   do j = 1,np_e
   do i = 1,iba(k)
     js = nis_e(myrank_e) + j-1
     t_wk2(i,js,k,l) = zinn(i,j,k,l)
   enddo
   enddo
   enddo
   enddo

   nn = kg1*neg*(iend_k-ista_k+1)*kimg
   call mpi_allreduce(t_wk2,t_wk1,nn,mpi_real8,mpi_sum, mpi_comm_world,err)

   if (myrank_e .eq. 0) then
      write(*,*)'called dump_zaj'
      c_img='_1'
      if (kimg == 2) c_img='_2'
      do i = 0,999
         write(c_num,'(i3.3)') i
         fname=psuffix//'_'//C_NUM//c_img//suffix
         inquire(file=trim(fname), exist=existence)
         if(existence) cycle
         if(.not.existence) exit
      enddo
      if (.not.existence) then
         open(unit=2100,file=fname,form='UNFORMATTED',status='NEW')
         write(2100) 4
         write(2100) kg1
         write(2100) neg
         write(2100) (iend_k-ista_k+1)
         write(2100) kimg
         write(2100) t_wk1
         close(2100)
      endif
   endif
   deallocate( t_wk1, t_wk2)
  end subroutine dump_zaj

!===============================================================================
  subroutine dump_eko(eko_l)

  use m_PlaneWaveBasisSet,  only : kg1, iba
  use m_Control_Parameters, only : nspin, kimg, neg, meg
  use m_Parallelization,    only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype &
       &                         , nrank_e,myrank_e,map_e,ista_e,iend_e,istep_e,idisp_e &
       &                         , map_z,np_e,mpi_k_world,myrank_k,map_k,ista_k,iend_k &
       &                         , nis_e, nie_e
  use m_IterationNumbers,   only : iteration_electronic, nk_in_the_process, iteration, nkgroup

  integer :: ierr,err

  real(kind=8),intent(in) :: eko_l(np_e,ista_k:iend_k)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, is, js, nn, ifrom
  real(kind=8), allocatable, dimension(:,:) :: t_wk1, t_wk2
  character(len=4) :: sort
  logical :: existence
  character(len=3) :: c_num
  character(len=80) :: fname
  character(len=8) :: psuffix = 'dump_eko'
#ifdef ORIGINAL
  character(len=5)  :: suffix = '.dat3'
#else
  character(len=5)  :: suffix = '.data'
#endif

  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

   allocate( t_wk1(neg,ista_k:iend_k) )
   allocate( t_wk2(neg,ista_k:iend_k) )

   t_wk1   = 0.0d0
   t_wk2   = 0.0d0

   do j = ista_k, iend_k
   do i = 1,np_e
     t_wk2(ista_e+i-1,j) = eko_l(i,j)
   enddo
   enddo

   nn = neg*(iend_k-ista_k+1)
   call mpi_allreduce(t_wk2,t_wk1,nn,mpi_real8,mpi_sum, mpi_k_world(myrank_k),err)

   if (myrank_e .eq. 0) then
      write(*,*)'called dump_eko'
      do i = 0,999
         write(c_num,'(i3.3)') i
         fname=psuffix//'_'//C_NUM//suffix
         inquire(file=trim(fname), exist=existence)
         if(existence) cycle
         if(.not.existence) exit
      enddo
      if (.not.existence) then
         open(unit=2100,file=fname,form='UNFORMATTED',status='NEW')
         write(2100) 2
         write(2100) neg
         write(2100) (iend_k-ista_k+1)
         write(2100) t_wk1
         close(2100)
      endif
   endif
   deallocate( t_wk1, t_wk2)
  end subroutine dump_eko

!===============================================================================
  subroutine dump_occup(eko_l)

  use m_PlaneWaveBasisSet,  only : kg1, iba
  use m_Control_Parameters, only : nspin, kimg, neg, meg
  use m_Parallelization,    only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype &
       &                         , nrank_e,myrank_e,map_e,ista_e,iend_e,istep_e,idisp_e &
       &                         , map_z,np_e,mpi_k_world,myrank_k,map_k,ista_k,iend_k &
       &                         , nis_e, nie_e
  use m_IterationNumbers,   only : iteration_electronic, nk_in_the_process, iteration, nkgroup

  integer :: ierr,err

  real(kind=8),intent(in) :: eko_l(np_e,ista_k:iend_k)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, is, js, nn, ifrom
  real(kind=8), allocatable, dimension(:,:) :: t_wk1, t_wk2
  character(len=4) :: sort
  logical :: existence
  character(len=3) :: c_num
  character(len=80) :: fname
  character(len=8) :: psuffix = 'dump_occup'
#ifdef ORIGINAL
  character(len=5)  :: suffix = '.dat3'
#else
  character(len=5)  :: suffix = '.data'
#endif

  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

   allocate( t_wk1(neg,ista_k:iend_k) )
   allocate( t_wk2(neg,ista_k:iend_k) )

   t_wk1   = 0.0d0
   t_wk2   = 0.0d0

   do j = ista_k, iend_k
   do i = 1,np_e
     t_wk2(ista_e+i-1,j) = eko_l(i,j)
   enddo
   enddo

   nn = neg*(iend_k-ista_k+1)
   call mpi_allreduce(t_wk2,t_wk1,nn,mpi_real8,mpi_sum, mpi_k_world(myrank_k),err)

   if (myrank_e .eq. 0) then
      write(*,*)'called dump_eko'
      do i = 0,999
         write(c_num,'(i3.3)') i
         fname=psuffix//'_'//C_NUM//suffix
         inquire(file=trim(fname), exist=existence)
         if(existence) cycle
         if(.not.existence) exit
      enddo
      if (.not.existence) then
         open(unit=2100,file=fname,form='UNFORMATTED',status='NEW')
         write(2100) 2
         write(2100) neg
         write(2100) (iend_k-ista_k+1)
         write(2100) t_wk1
         close(2100)
      endif
   endif
   deallocate( t_wk1, t_wk2)
  end subroutine dump_occup

!===============================================================================
  subroutine dump_neordr(neordr)

  use m_PlaneWaveBasisSet,  only : kg1, iba
  use m_Control_Parameters, only : nspin, kimg, neg, meg
  use m_Parallelization,    only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype &
       &                         , nrank_e,myrank_e,map_e,ista_e,iend_e,istep_e,idisp_e &
       &                         , map_z,np_e,mpi_k_world,myrank_k,map_k,ista_k,iend_k &
       &                         , nis_e, nie_e
  use m_IterationNumbers,   only : iteration_electronic, nk_in_the_process, iteration, nkgroup

  integer :: ierr,err

  real(kind=8),intent(in) :: neordr(neg,ista_k:iend_k)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, is, js, nn, ifrom
  real(kind=8), allocatable, dimension(:,:) :: t_wk1, t_wk2
  character(len=4) :: sort
  logical :: existence
  character(len=3) :: c_num
  character(len=80) :: fname
  character(len=8) :: psuffix = 'dump_neordr'
#ifdef ORIGINAL
  character(len=5)  :: suffix = '.dat3'
#else
  character(len=5)  :: suffix = '.data'
#endif

  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI


   if (myrank_e .eq. 0) then
      write(*,*)'called dump_neordr'
      do i = 0,999
         write(c_num,'(i3.3)') i
         fname=psuffix//'_'//C_NUM//suffix
         inquire(file=trim(fname), exist=existence)
         if(existence) cycle
         if(.not.existence) exit
      enddo
      if (.not.existence) then
         open(unit=2100,file=fname,form='UNFORMATTED',status='NEW')
         write(2100) 2
         write(2100) neg
         write(2100) (iend_k-ista_k+1)
         write(2100) neordr
         close(2100)
      endif
   endif
  end subroutine dump_neordr

!===============================================================================
  subroutine dump_nrvf_ordr(nrvf_ordr)

  use m_PlaneWaveBasisSet,  only : kg1, iba
  use m_Control_Parameters, only : nspin, kimg, neg, meg
  use m_Parallelization,    only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype &
       &                         , nrank_e,myrank_e,map_e,ista_e,iend_e,istep_e,idisp_e &
       &                         , map_z,np_e,mpi_k_world,myrank_k,map_k,ista_k,iend_k &
       &                         , nis_e, nie_e
  use m_IterationNumbers,   only : iteration_electronic, nk_in_the_process, iteration, nkgroup

  integer :: ierr,err

  real(kind=8),intent(in) :: nrvf_ordr(neg,ista_k:iend_k)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, is, js, nn, ifrom
  real(kind=8), allocatable, dimension(:,:) :: t_wk1, t_wk2
  character(len=4) :: sort
  logical :: existence
  character(len=3) :: c_num
  character(len=80) :: fname
  character(len=8) :: psuffix = 'dump_nrvf_ordr'
#ifdef ORIGINAL
  character(len=5)  :: suffix = '.dat3'
#else
  character(len=5)  :: suffix = '.data'
#endif

  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI


   if (myrank_e .eq. 0) then
      write(*,*)'called dump_nrvf_ordr'
      do i = 0,999
         write(c_num,'(i3.3)') i
         fname=psuffix//'_'//C_NUM//suffix
         inquire(file=trim(fname), exist=existence)
         if(existence) cycle
         if(.not.existence) exit
      enddo
      if (.not.existence) then
         open(unit=2100,file=fname,form='UNFORMATTED',status='NEW')
         write(2100) 2
         write(2100) neg
         write(2100) (iend_k-ista_k+1)
         write(2100) nrvf_ordr
         close(2100)
      endif
   endif
  end subroutine dump_nrvf_ordr


!===============================================================================
  subroutine dump_chgq(chgq_l)

  use m_PlaneWaveBasisSet,  only : kg1, iba, kgp
  use m_Control_Parameters, only : nspin, kimg, neg, meg
  use m_Parallelization,    only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype &
       &                         , nrank_e,myrank_e,map_e,ista_e,iend_e,istep_e,idisp_e &
       &                         , map_z,np_e,mpi_k_world,myrank_k,map_k,ista_k,iend_k &
       &                         , nis_e, nie_e, ista_kngp, iend_kngp
  use m_IterationNumbers,   only : iteration_electronic, nk_in_the_process, iteration, nkgroup

  integer :: ierr,err

  real(kind=8),intent(in) :: chgq_l(ista_kngp:iend_kngp,kimg,nspin)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, is, js, nn, ifrom
  real(kind=8), allocatable, dimension(:,:,:) :: t_wk1, t_wk2
  character(len=4) :: sort
  logical :: existence
  character(len=3) :: c_num, c_iter
  character(len=80) :: fname
  character(len=8) :: psuffix = 'dump_chgq'
#ifdef ORIGINAL
  character(len=5)  :: suffix = '.dat3'
#else
  character(len=5)  :: suffix = '.data'
#endif

  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

   allocate( t_wk1(kgp,kimg,nspin) )
   allocate( t_wk2(kgp,kimg,nspin) )

   t_wk1   = 0.0d0
   t_wk2   = 0.0d0

   do k = 1, nspin
   do j = 1, kimg
   do i = ista_kngp, iend_kngp
     t_wk2(i,j,k) = chgq_l(i,j,k)
   enddo
   enddo
   enddo

   call mpi_allreduce(t_wk2,t_wk1,kgp*kimg*nspin,mpi_real8,mpi_sum, mpi_k_world(myrank_k),err)

   if (myrank_e .eq. 0) then
      write(*,*)'called dump_chgq'
      do i = 0,999
         write(c_num,'(i3.3)') i
         write(c_iter,'(i3.3)') iteration
         fname=psuffix//'_'//c_iter//'_'//C_NUM//suffix
         inquire(file=trim(fname), exist=existence)
         if(existence) cycle
         if(.not.existence) exit
      enddo
      if (.not.existence) then
         open(unit=2100,file=fname,form='UNFORMATTED',status='NEW')
         write(2100) 3
         write(2100) kgp
         write(2100) kimg
         write(2100) nspin
         write(2100) t_wk1
         close(2100)
      endif
   endif
   deallocate( t_wk1, t_wk2)
  end subroutine dump_chgq

!===============================================================================
  subroutine dump_vxc(vxc_l)

  use m_PlaneWaveBasisSet,  only : kg1, iba, kgp
  use m_Control_Parameters, only : nspin, kimg, neg, meg
  use m_Parallelization,    only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype &
       &                         , nrank_e,myrank_e,map_e,ista_e,iend_e,istep_e,idisp_e &
       &                         , map_z,np_e,mpi_k_world,myrank_k,map_k,ista_k,iend_k &
       &                         , nis_e, nie_e, ista_kngp, iend_kngp
  use m_IterationNumbers,   only : iteration_electronic, nk_in_the_process, iteration, nkgroup

  integer :: ierr,err

  real(kind=8),intent(in) :: vxc_l(ista_kngp:iend_kngp,kimg,nspin)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, is, js, nn, ifrom
  real(kind=8), allocatable, dimension(:,:,:) :: t_wk1, t_wk2
  character(len=4) :: sort
  logical :: existence
  character(len=3) :: c_num, c_iter
  character(len=80) :: fname
  character(len=8) :: psuffix = 'dump_vxc'
#ifdef ORIGINAL
  character(len=5)  :: suffix = '.dat3'
#else
  character(len=5)  :: suffix = '.data'
#endif

  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

   allocate( t_wk1(kgp,kimg,nspin) )
   allocate( t_wk2(kgp,kimg,nspin) )

   t_wk1   = 0.0d0
   t_wk2   = 0.0d0

   do k = 1, nspin
   do j = 1, kimg
   do i = ista_kngp, iend_kngp
     t_wk2(i,j,k) = vxc_l(i,j,k)
   enddo
   enddo
   enddo

   call mpi_allreduce(t_wk2,t_wk1,kgp*kimg*nspin,mpi_real8,mpi_sum, mpi_k_world(myrank_k),err)

   if (myrank_e .eq. 0) then
      write(*,*)'called dump_vxc'
      do i = 0,999
         write(c_num,'(i3.3)') i
         write(c_iter,'(i3.3)') iteration
         fname=psuffix//'_'//c_iter//'_'//C_NUM//suffix
         inquire(file=trim(fname), exist=existence)
         if(existence) cycle
         if(.not.existence) exit
      enddo
      if (.not.existence) then
         open(unit=2100,file=fname,form='UNFORMATTED',status='NEW')
         write(2100) 3
         write(2100) kgp
         write(2100) kimg
         write(2100) nspin
         write(2100) t_wk1
         close(2100)
      endif
   endif
   deallocate( t_wk1, t_wk2)
  end subroutine dump_vxc

!===============================================================================
  subroutine dump_vlhxcQ(vlhxcQ)

  use m_PseudoPotential,    only : nlmt
  use m_Ionic_System,       only : natm
  use m_Control_Parameters, only : nspin, kimg
  use m_IterationNumbers,   only : iteration
  use m_Parallelization,    only : myrank_e

  integer :: ierr,err

  real(kind=8),intent(in) :: vlhxcQ(nlmt,nlmt,natm,nspin)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, is, js, nn, ifrom
  character(len=4) :: sort
  logical :: existence
  character(len=3) :: c_num, c_iter
  character(len=80) :: fname
  character(len=11) :: psuffix = 'dump_vlhxcQ'
#ifdef ORIGINAL
  character(len=5)  :: suffix = '.dat3'
#else
  character(len=5)  :: suffix = '.data'
#endif

   if (myrank_e .eq. 0) then
      write(*,*)'called dump_vlhxcQ'
      do i = 0,999
         write(c_num,'(i3.3)') i
         write(c_iter,'(i3.3)') iteration
         fname=psuffix//'_'//c_iter//'_'//C_NUM//suffix
         inquire(file=trim(fname), exist=existence)
         if(existence) cycle
         if(.not.existence) exit
      enddo
      if (.not.existence) then
         open(unit=2100,file=fname,form='UNFORMATTED',status='NEW')
         write(2100) 4
         write(2100) nlmt
         write(2100) nlmt
         write(2100) natm
         write(2100) nspin
         write(2100) vlhxcQ
         close(2100)
      endif
   endif
  end subroutine dump_vlhxcQ

!===============================================================================
  subroutine dump_snl(zinn)

  use m_PlaneWaveBasisSet,  only : kg1, iba
  use m_Control_Parameters, only : nspin, kimg, neg, meg
  use m_Parallelization,    only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype &
       &                         , nrank_e,myrank_e,map_e,ista_e,iend_e,istep_e,idisp_e &
       &                         , map_z,np_e,mpi_k_world,myrank_k,map_k,ista_k,iend_k &
       &                         , nis_e, nie_e,mype
  use m_IterationNumbers,   only : iteration_electronic, nk_in_the_process, iteration, nkgroup
  use m_Parallelization,    only : ista_snl, iend_snl
  use m_PseudoPotential,    only : nlmtt

  integer :: ierr,err
  real(kind=8) :: zinn(kg1,nlmtt,ista_snl:iend_snl)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, ik, is, js, nn, ifrom
  real(kind=8), allocatable, dimension(:,:,:,:) :: t_wk1, t_wk2
  character(len=4) :: sort
  logical :: existence
  character(len=3) :: c_num
  character(len=2) :: c_img
  character(len=80) :: fname



#if 1
!   if (myrank_e .eq. 0) then
   if (mype .eq. 0) then
      write(*,*)'called dump snl'
      c_img='_1'
      if (kimg == 2) c_img='_2'
      do i = 0,999
         write(c_num,'(i3.3)') i
!         fname='dump_snl'//'.'//C_NUM//c_img//'.dat'
         fname='dump_snl'//'_'//C_NUM//'.data'
         inquire(file=trim(fname), exist=existence)
         if(existence) cycle
         if(.not.existence) exit
      enddo
      if (.not.existence) then
         open(unit=2100,file=fname,form='UNFORMATTED',status='NEW')
         write(2100) 3
         write(2100) kg1
         write(2100) nlmtt
         write(2100) (iend_snl-ista_snl+1)
         write(2100) zinn
         close(2100)
      endif
   endif
#endif
  end subroutine dump_snl

!===============================================================================
  subroutine dump_qitg(qitg_l,n)

  use m_PlaneWaveBasisSet,  only : kg1, iba, kgp
  use m_Control_Parameters, only : nspin, kimg, neg, meg
  use m_Parallelization,    only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype &
       &                         , nrank_e,myrank_e,map_e,ista_e,iend_e,istep_e,idisp_e &
       &                         , map_z,np_e,mpi_k_world,myrank_k,map_k,ista_k,iend_k &
       &                         , nis_e, nie_e, ista_kngp, iend_kngp
  use m_IterationNumbers,   only : iteration_electronic, nk_in_the_process, iteration, nkgroup

  integer :: ierr,err

  real(kind=8),intent(in) :: qitg_l(ista_kngp:iend_kngp,n)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, is, js, nn, ifrom
  real(kind=8), allocatable, dimension(:,:) :: t_wk1, t_wk2
  character(len=4) :: sort
  logical :: existence
  character(len=3) :: c_num, c_iter
  character(len=80) :: fname
  character(len=8) :: psuffix = 'dump_qitg'
#ifdef ORIGINAL
  character(len=5)  :: suffix = '.dat3'
#else
  character(len=5)  :: suffix = '.data'
#endif

  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

   allocate( t_wk1(kgp,n) )
   allocate( t_wk2(kgp,n) )

   t_wk1   = 0.0d0
   t_wk2   = 0.0d0

   do j = 1, n
   do i = ista_kngp, iend_kngp
     t_wk2(i,j) = qitg_l(i,j)
   enddo
   enddo

   call mpi_allreduce(t_wk2,t_wk1,kgp*n,mpi_real8,mpi_sum, mpi_k_world(myrank_k),err)

   if (myrank_e .eq. 0) then
      write(*,*)'called dump_qitg'
      do i = 0,999
         write(c_num,'(i3.3)') i
         write(c_iter,'(i3.3)') iteration
         fname=psuffix//'_'//c_iter//'_'//C_NUM//suffix
         inquire(file=trim(fname), exist=existence)
         if(existence) cycle
         if(.not.existence) exit
      enddo
      if (.not.existence) then
         open(unit=2100,file=fname,form='UNFORMATTED',status='NEW')
         write(2100) 2
         write(2100) kgp
         write(2100) n
         write(2100) t_wk1
         close(2100)
      endif
   endif
   deallocate( t_wk1, t_wk2)
  end subroutine dump_qitg

!===============================================================================
  subroutine dump_amat(amat,ix,iy)

  use m_PseudoPotential,    only : nlmt
  use m_Ionic_System,       only : natm
  use m_Control_Parameters, only : nspin, kimg
  use m_IterationNumbers,   only : iteration
  use m_Parallelization,    only : myrank_e,mype

  integer :: ierr,err

  real(kind=8),intent(in) :: amat(ix,iy)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, is, js, nn, ifrom
  character(len=4) :: sort
  logical :: existence
  character(len=3) :: c_num, c_iter
  character(len=80) :: fname
  character(len=9) :: psuffix = 'dump_amat'

  character(len=5)  :: suffix = '.data'

   if (mype .eq. 0) then
      write(*,*)'called dump_amat'
      do i = 0,999
         write(c_num,'(i3.3)') i
         write(c_iter,'(i3.3)') iteration
         fname=psuffix//'_'//c_iter//'_'//C_NUM//suffix
         inquire(file=trim(fname), exist=existence)
         if(existence) cycle
         if(.not.existence) exit
      enddo
      if (.not.existence) then
         open(unit=2100,file=fname,form='UNFORMATTED',status='NEW')
         write(2100) 2
         write(2100) ix
         write(2100) iy
         write(2100) amat
         close(2100)
      endif
   endif
  end subroutine dump_amat

  subroutine dump_zajmat3D(zinn)

  use m_PlaneWaveBasisSet,  only : kg1, iba,iba2,iba,nmatsz,nmatsz2
  use m_Control_Parameters, only : nspin, kimg, neg, meg
  use m_Parallelization,    only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype &
       &                         , nrank_e,myrank_e,map_e,ista_e,iend_e,istep_e,idisp_e &
       &                         , map_z,np_e_3D,mpi_k_world,myrank_k,map_k,ista_k,iend_k &
       &                         , nis_e, nie_e,mype,mpi_kg_world,neg_g
  use m_IterationNumbers,   only : iteration_electronic, nk_in_the_process, iteration, nkgroup

  integer :: ierr,err

  real(kind=8),intent(in) :: zinn(nmatsz,np_e_3D,kimg)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, is, js, nn, ifrom
  real(kind=8), allocatable, dimension(:,:,:) :: t_wk1, t_wk2
  character(len=4) :: sort
  logical :: existence
  character(len=3) :: c_num
  character(len=2) :: c_img
  character(len=80) :: fname
  character(len=9) :: psuffix = 'dump_zajm'
#ifdef ORIGINAL
  character(len=5)  :: suffix = '.data'
#else
  character(len=5)  :: suffix = '.dat1'
#endif

  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

   allocate( t_wk1(nmatsz,neg,kimg) )
   allocate( t_wk2(nmatsz,neg,kimg) )

   t_wk1   = 0.0d0
   t_wk2   = 0.0d0

   do l = 1, kimg
   do j = 1,np_e_3D
   do i = 1,nmatsz
     js = neg_g(j)
     t_wk2(i,js,l) = zinn(i,j,l)
   enddo
   enddo
   enddo

   nn = nmatsz*neg*kimg
   call mpi_allreduce(t_wk2,t_wk1,nn,mpi_real8,mpi_sum, mpi_kg_world,err)

   if (mype .eq. 0) then
      write(*,*)'called dump_zaj'
      c_img='_1'
      if (kimg == 2) c_img='_2'
      do i = 0,999
         write(c_num,'(i3.3)') i
         fname=psuffix//'_'//C_NUM//c_img//suffix
         inquire(file=trim(fname), exist=existence)
         if(existence) cycle
         if(.not.existence) exit
      enddo
      if (.not.existence) then
         open(unit=2100,file=fname,form='UNFORMATTED',status='NEW')
         write(2100) 3
         write(2100) nmatsz
         write(2100) neg
         write(2100) kimg
         write(2100) t_wk1
         close(2100)
      endif
   endif
   deallocate( t_wk1, t_wk2)
  end subroutine dump_zajmat3D

  subroutine dump_eko3D(eko_l)

  use m_PlaneWaveBasisSet,  only : kg1, iba
  use m_Control_Parameters, only : nspin, kimg, neg, meg
  use m_Parallelization,    only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype &
       &                         , nrank_e,myrank_e,map_e,ista_e,iend_e,istep_e,idisp_e &
       &                         , map_z,np_e_3D,mpi_k_world,myrank_k,map_k,ista_k_3D,iend_k_3D &
       &                         , nis_e, nie_e,mype,mpi_kg_world,neg_g
  use m_IterationNumbers,   only : iteration_electronic, nk_in_the_process, iteration, nkgroup

  integer :: ierr,err

  real(kind=8),intent(in) :: eko_l(np_e_3D,ista_k_3D:iend_k_3D)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, is, js, nn, ifrom
  real(kind=8), allocatable, dimension(:,:) :: t_wk1, t_wk2
  character(len=4) :: sort
  logical :: existence
  character(len=3) :: c_num
  character(len=80) :: fname
  character(len=8) :: psuffix = 'dump_eko3'
#ifdef ORIGINAL
  character(len=5)  :: suffix = '.data'
#else
  character(len=5)  :: suffix = '.data'
#endif

  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

   allocate( t_wk1(neg,ista_k_3D:iend_k_3D) )
   allocate( t_wk2(neg,ista_k_3D:iend_k_3D) )

   t_wk1   = 0.0d0
   t_wk2   = 0.0d0

   do j = ista_k_3D, iend_k_3D
   do i = 1,np_e_3D
     t_wk2(neg_g(i),j) = eko_l(i,j)
   enddo
   enddo

   nn = neg*(iend_k_3D-ista_k_3D+1)
   call mpi_allreduce(t_wk2,t_wk1,nn,mpi_real8,mpi_sum, mpi_kg_world,err)

   if (mype .eq. 0) then
      write(*,*)'called dump_eko'
      do i = 0,999
         write(c_num,'(i3.3)') i
         fname=psuffix//'_'//C_NUM//suffix
         inquire(file=trim(fname), exist=existence)
         if(existence) cycle
         if(.not.existence) exit
      enddo
      if (.not.existence) then
         open(unit=2100,file=fname,form='UNFORMATTED',status='NEW')
         write(2100) 2
         write(2100) neg
         write(2100) (iend_k_3D-ista_k_3D+1)
         write(2100) t_wk1
         close(2100)
      endif
   endif
   deallocate( t_wk1, t_wk2)
  end subroutine dump_eko3D


  subroutine dump_dbg(amat,ix,iy)

  use m_PseudoPotential,    only : nlmt
  use m_Ionic_System,       only : natm
  use m_Control_Parameters, only : nspin, kimg
  use m_IterationNumbers,   only : iteration
  use m_Parallelization,    only : myrank_e,mype

  integer :: ierr,err

  real(kind=8),intent(in) :: amat(ix,iy)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, is, js, nn, ifrom
  character(len=4) :: sort
  logical :: existence
  character(len=3) :: c_num, c_iter
  character(len=80) :: fname
  character(len=8) :: psuffix = 'dump_dbg'

  character(len=5)  :: suffix = '.data'

   if (mype .eq. 0) then
      write(*,*)'called dump_amat'
      do i = 0,999
         write(c_num,'(i3.3)') i
         write(c_iter,'(i3.3)') iteration
         fname=psuffix//'_'//c_iter//'_'//C_NUM//suffix
         inquire(file=trim(fname), exist=existence)
         if(existence) cycle
         if(.not.existence) exit
      enddo
      if (.not.existence) then
         open(unit=2100,file=fname,form='UNFORMATTED',status='NEW')
         write(2100) 2
         write(2100) ix
         write(2100) iy
         write(2100) amat
         close(2100)
      endif
   endif
  end subroutine dump_dbg

  subroutine dump_dbg_int(amat,ix,iy)

  use m_PseudoPotential,    only : nlmt
  use m_Ionic_System,       only : natm
  use m_Control_Parameters, only : nspin, kimg
  use m_IterationNumbers,   only : iteration
  use m_Parallelization,    only : myrank_e,mype

  integer :: ierr,err

  integer,intent(in) :: amat(ix,iy)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, is, js, nn, ifrom
  character(len=4) :: sort
  logical :: existence
  character(len=3) :: c_num, c_iter
  character(len=80) :: fname
  character(len=8) :: psuffix = 'dump_int'

  character(len=5)  :: suffix = '.data'

   if (mype .eq. 0) then
      write(*,*)'called dump_amat'
      do i = 0,999
         write(c_num,'(i3.3)') i
         write(c_iter,'(i3.3)') iteration
         fname=psuffix//'_'//c_iter//'_'//C_NUM//suffix
         inquire(file=trim(fname), exist=existence)
         if(existence) cycle
         if(.not.existence) exit
      enddo
      if (.not.existence) then
         open(unit=2100,file=fname,form='UNFORMATTED',status='NEW')
         write(2100) 2
         write(2100) ix
         write(2100) iy
         write(2100) amat
         close(2100)
      endif
   endif
  end subroutine dump_dbg_int

  subroutine dump_amat_c(amat,ix,iy)

  use m_PseudoPotential,    only : nlmt
  use m_Ionic_System,       only : natm
  use m_Control_Parameters, only : nspin, kimg
  use m_IterationNumbers,   only : iteration
  use m_Parallelization,    only : myrank_e,mype
  use m_Const_Parameters,    only : CMPLDP

  integer :: ierr,err

  complex(kind=CMPLDP),intent(in) :: amat(ix,iy)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, is, js, nn, ifrom
  character(len=4) :: sort
  logical :: existence
  character(len=3) :: c_num, c_iter
  character(len=80) :: fname
  character(len=9) :: psuffix = 'dump_amat'

  character(len=5)  :: suffix = '.data'

   if (mype .eq. 0) then
      write(*,*)'called dump_amat'
      do i = 0,999
         write(c_num,'(i3.3)') i
         write(c_iter,'(i3.3)') iteration
         fname=psuffix//'_'//c_iter//'_'//C_NUM//suffix
         inquire(file=trim(fname), exist=existence)
         if(existence) cycle
         if(.not.existence) exit
      enddo
      if (.not.existence) then
         open(unit=2100,file=fname,form='UNFORMATTED',status='NEW')
         write(2100) 2
         write(2100) ix
         write(2100) iy
         write(2100) amat
         close(2100)
      endif
   endif
  end subroutine dump_amat_c

!===============================================================================

!===============================================================================
#ifndef ORIGINAL
  subroutine dump_many(string)
  use m_Electronic_Structure,only : zaj_l_3D,zaj_ball &
 &                                , eko_l_3D,vnlph_l_3D,nrvf_ordr,neordr &
 &                                , fsr_l_3D,fsi_l_3D ,fsr_l_2d,fsi_l_2d,fsr_l,fsi_l  &
 &                                , vlhxc_l, vlhxc_l_3D, occup_l_3D
  use m_NonLocal_Potential,  only : snl, snl_l_3D
  use m_Parallelization,     only : np_e_3D, ista_k_3D, iend_k_3D, ista_snl, iend_snl, np_g1k_3D &
 &                                , ista_kngp_3D, iend_kngp_3D
  use m_Control_Parameters,  only : neg
  use m_NonLocal_Potential,  only : snl, snl_l_3D
  use m_ES_WF_by_SDorCG,     only : zaj_old, zaj_old_3D

  character(len=32),intent(in) :: string
  integer :: i, j, k, l, m, n, is, js, nn, ifrom, IJK=40
#if 0
  write(1002,'("***** ",32a)') string

  write(1002,'("zaj_l_3D(1,:,1,1)")')
  write(1002,'(10(1x,f9.6))')  (zaj_l_3D(1,i,1,1),i=1,np_e_3D)
  write(1002,'("zaj_l_3D(2,:,1,1)")')
  write(1002,'(10(1x,f9.6))')  (zaj_l_3D(2,i,1,1),i=1,np_e_3D)

  if (allocated(zaj_old_3D)) then
  write(1002,'("zaj_old_3D(1,:,1,1)")')
  write(1002,'(10(1x,f9.6))')  (zaj_old_3D(1,i,1,1),i=1,np_e_3D)
  write(1002,'("zaj_old_3D(2,:,1,1)")')
  write(1002,'(10(1x,f9.6))')  (zaj_old_3D(2,i,1,1),i=1,np_e_3D)
  endif

  write(1002,'("zaj_ball(1,:,1,1)")')
  write(1002,'(10(1x,f9.6))')  (zaj_ball(1,i,1,1),i=1,min(neg,IJK))
  write(1002,'("zaj_ball(2,:,1,1)")')
  write(1002,'(10(1x,f9.6))')  (zaj_ball(2,i,1,1),i=1,min(neg,IJK))

  write(1002,'("fsr_l_3D(:,1,1)")')
  write(1002,'(10(1x,f9.6))')  (fsr_l_3D(i,1,1),i=1,min(np_e_3D,IJK))
  write(1002,'("fsr_l_3D(:,2,1)")')
  write(1002,'(10(1x,f9.6))')  (fsr_l_3D(i,2,1),i=1,min(np_e_3D,IJK))

  write(1002,'("eko_l_3D(:,1)")')
  write(1002,'(10(1x,f9.6))')  (eko_l_3D(i,ista_k_3D),i=1,min(np_e_3D,IJK))

  write(1002,'("occup_l_3D(:,1)")')
  write(1002,'(10(1x,f9.6))')  (occup_l_3D(i,ista_k_3D),i=1,min(np_e_3D,IJK))

  write(1002,'("vnlph_l_3D(1,:,1)")')
  write(1002,'(10(1x,f9.6))')  (vnlph_l_3D(1,i,1),i=1,min(np_e_3D,IJK))
  write(1002,'("vnlph_l_3D(2,:,1)")')
  write(1002,'(10(1x,f9.6))')  (vnlph_l_3D(2,i,1),i=1,min(np_e_3D,IJK))

! write(1002,'("snl_l_3D(:,1,1)")')
! write(1002,'(10(1x,f9.6))')  (snl_l_3D(i,1,ista_snl),i=1,min(np_g1k_3D(1),IJK))
! write(1002,'("snl_l_3D(:,2,1)")')
! write(1002,'(10(1x,f9.6))')  (snl_l_3D(i,2,ista_snl),i=1,min(np_g1k_3D(1),IJK))
! write(1002,'("snl_l_3D(:,3,1)")')
! write(1002,'(10(1x,f9.6))')  (snl_l_3D(i,3,ista_snl),i=1,min(np_g1k_3D(1),IJK))

  write(1002,'("vlhxc_l_3D(:,1,1)")')
  write(1002,'(10(1x,f9.6))')  (vlhxc_l_3D(i,1,1),i=ista_kngp_3D,min(ista_kngp_3D+IJK,iend_kngp_3D))
#endif
  end subroutine dump_many

!===============================================================================
  subroutine dump_many_ik(ik,string)
  use m_Electronic_Structure,only : zaj_l_3D,zaj_ball &
 &                                , eko_l_3D,vnlph_l_3D,nrvf_ordr,neordr &
 &                                , fsr_l_3D,fsi_l_3D ,fsr_l_2d,fsi_l_2d,fsr_l,fsi_l  &
 &                                , vlhxc_l, vlhxc_l_3D, occup_l_3D
  use m_NonLocal_Potential,  only : snl, snl_l_3D
  use m_Parallelization,     only : np_e_3D, ista_k_3D, iend_k_3D, ista_snl, iend_snl, np_g1k_3D &
 &                                , ista_kngp_3D, iend_kngp_3D
  use m_Control_Parameters,  only : neg
  use m_NonLocal_Potential,  only : snl, snl_l_3D
  use m_ES_WF_by_SDorCG,     only : zaj_old, zaj_old_3D

  character(len=32),intent(in) :: string
  integer          ,intent(in) :: ik
  integer :: i, j, k, l, m, n, is, js, nn, ifrom, IJK=40
#if 0
  write(1002,'("***** ",32a)') string

  write(1002,'("IK:",i4)') ik
  write(1002,'("zaj_l_3D(1,:,ik,1)")')
  write(1002,'(10(1x,f9.6))')  (zaj_l_3D(1,i,ik,1),i=1,np_e_3D)
  write(1002,'("zaj_l_3D(2,:,ik,1)")')
  write(1002,'(10(1x,f9.6))')  (zaj_l_3D(2,i,ik,1),i=1,np_e_3D)

  if (allocated(zaj_old_3D)) then
  write(1002,'("zaj_old_3D(1,:,ik,1)")')
  write(1002,'(10(1x,f9.6))')  (zaj_old_3D(1,i,ik,1),i=1,np_e_3D)
  write(1002,'("zaj_old_3D(2,:,ik,1)")')
  write(1002,'(10(1x,f9.6))')  (zaj_old_3D(2,i,ik,1),i=1,np_e_3D)
  endif

  write(1002,'("zaj_ball(1,:,ik,1)")')
  write(1002,'(10(1x,f9.6))')  (zaj_ball(1,i,ik,1),i=1,min(neg,IJK))
  write(1002,'("zaj_ball(2,:,ik,1)")')
  write(1002,'(10(1x,f9.6))')  (zaj_ball(2,i,ik,1),i=1,min(neg,IJK))

  write(1002,'("fsr_l_3D(:,1,ik)")')
  write(1002,'(10(1x,f9.6))')  (fsr_l_3D(i,1,ik),i=1,min(np_e_3D,IJK))
  write(1002,'("fsr_l_3D(:,2,ik)")')
  write(1002,'(10(1x,f9.6))')  (fsr_l_3D(i,2,ik),i=1,min(np_e_3D,IJK))

  write(1002,'("eko_l_3D(:,ik)")')
  write(1002,'(10(1x,f9.6))')  (eko_l_3D(i,ik),i=1,min(np_e_3D,IJK))

  write(1002,'("occup_l_3D(:,ik)")')
  write(1002,'(10(1x,f9.6))')  (occup_l_3D(i,ik),i=1,min(np_e_3D,IJK))

  write(1002,'("vnlph_l_3D(1,:,1)")')
  write(1002,'(10(1x,f9.6))')  (vnlph_l_3D(1,i,1),i=1,min(np_e_3D,IJK))
  write(1002,'("vnlph_l_3D(2,:,1)")')
  write(1002,'(10(1x,f9.6))')  (vnlph_l_3D(2,i,1),i=1,min(np_e_3D,IJK))

! write(1002,'("snl_l_3D(:,1,1)")')
! write(1002,'(10(1x,f9.6))')  (snl_l_3D(i,1,ik),i=1,min(np_g1k_3D(ik),IJK))
! write(1002,'("snl_l_3D(:,2,1)")')
! write(1002,'(10(1x,f9.6))')  (snl_l_3D(i,2,ik),i=1,min(np_g1k_3D(ik),IJK))
! write(1002,'("snl_l_3D(:,3,1)")')
! write(1002,'(10(1x,f9.6))')  (snl_l_3D(i,3,ik),i=1,min(np_g1k_3D(ik),IJK))

  write(1002,'("vlhxc_l_3D(:,1,1)")')
  write(1002,'(10(1x,f9.6))')  (vlhxc_l_3D(i,1,1),i=ista_kngp_3D,min(ista_kngp_3D+IJK,iend_kngp_3D))
#endif
  end subroutine dump_many_ik

#endif


#endif  ! if0
!===============================================================================
  subroutine dump_fsr(zinn)

  use m_Control_Parameters, only : neg
  use m_PseudoPotential,    only : nlmta
  use m_Parallelization,    only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype &
       &                         , nrank_e,myrank_e,map_e,ista_e,iend_e,istep_e,idisp_e &
       &                         , map_z,np_e,mpi_k_world,myrank_k,map_k,ista_k,iend_k &
       &                         , nis_e, nie_e,ista_fs,iend_fs,np_fs
  use m_IterationNumbers,   only : iteration_electronic, nk_in_the_process, iteration, nkgroup

  integer :: ierr,err

  real(kind=8),intent(in) :: zinn(np_e,np_fs,ista_k:iend_k)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, is, js, nn, ifrom
  real(kind=8), allocatable, dimension(:,:,:) :: t_wk1, t_wk2
  character(len=4) :: sort
  logical :: existence
  character(len=3) :: c_num
  character(len=2) :: c_img
  character(len=80) :: fname
  character(len=8) :: psuffix = 'dump_fsr'
#ifdef ORIGINAL
  character(len=5)  :: suffix = '.data'
#else
  character(len=5)  :: suffix = '.data'
#endif

  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

   allocate( t_wk1(neg,nlmta,ista_k:iend_k) )
   allocate( t_wk2(neg,nlmta,ista_k:iend_k) )

   t_wk1   = 0.0d0
   t_wk2   = 0.0d0

   do k = ista_k, iend_k
   do j = 1,np_fs
   do i = 1,np_e
     js = ista_fs + j - 1
     is = ista_e  + i - 1
     t_wk2(is,js,k) = zinn(i,j,k)
   enddo
   enddo
   enddo

   nn = neg*nlmta*(iend_k-ista_k+1)
   call mpi_allreduce(t_wk2,t_wk1,nn,mpi_real8,mpi_sum, mpi_comm_world,err)

!   if (myrank_e .eq. 0) then
   if (mype .eq. 0) then
      write(*,*)'called dump_fsr'
      c_img='_1'
      if (kimg == 2) c_img='_2'
      do i = 0,999
         write(c_num,'(i3.3)') i
         fname=psuffix//'_'//C_NUM//c_img//suffix
         inquire(file=trim(fname), exist=existence)
         if(existence) cycle
         if(.not.existence) exit
      enddo
      if (.not.existence) then
         open(unit=2100,file=fname,form='UNFORMATTED',status='NEW')
         write(2100) 3
         write(2100) neg
         write(2100) nlmta
         write(2100) (iend_k-ista_k+1)
         write(2100) t_wk1
         close(2100)
      endif
   endif
   deallocate( t_wk1, t_wk2)
  end subroutine dump_fsr

!===============================================================================

