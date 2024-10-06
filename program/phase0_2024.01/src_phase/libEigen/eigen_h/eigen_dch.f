      subroutine eigen_dch(n,d,e,z,ldz,info)
!
      use communication_h, only : eigen_init
     &               , eigen_free
     &               , bcast_dbl
     &               , cyc2d_cyc1d
!
      integer, intent(in)           :: n, ldz
      real(8), intent(inout)        :: d(1:n), e(1:n)
      real(8), intent(out)          :: z(1:ldz,*)

! parameters blacs array descritor(the position of entry tags), etc
      integer, parameter            :: block_cyclic_2d = 1
      integer, parameter            :: dlen_  = 9
      integer, parameter            :: dtype_ = 1
      integer, parameter            :: ctxt_  = 2
      integer, parameter            :: m_     = 3
      integer, parameter            :: n_     = 4
      integer, parameter            :: mb_    = 5
      integer, parameter            :: nb_    = 6
      integer, parameter            :: rsrc_  = 7
      integer, parameter            :: csrc_  = 8
      integer, parameter            :: lld_   = 9

      integer                       :: descz( dlen_ )
      integer                       :: descw( dlen_ )
 
      integer                       :: world_size, my_rank
      integer                       :: trilwmin, lwork, liwork
 
      real(8), pointer              :: work(:)
      integer, pointer              :: iwork(:)

      real(8)                       :: s0, s1, d1, d2, d3

      include 'mpif.h'
      include 'trd.h'
 
      s0=mpi_wtime()
!! blacs/pblas/scalapack initialization
!      call blacs_pinfo( iam, nprocs )
!      if ( nprocs < 1 ) then
!!  mpi group setup
!         call mpi_comm_size( mpi_comm_eigen, world_size, ierr )
!         call mpi_comm_rank( mpi_comm_eigen, my_rank, ierr )
!         nprocs = world_size
!         iam = my_rank
!         call blacs_setup( iam, nprocs )
!      end if
!      call blacs_get( -1, 0, ictxt )
!      call eigen_init(2)
      nprow = size_of_col
      npcol = size_of_row
!      call blacs_gridinit( ictxt, 'column-major', nprow, npcol )
!      call blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )
      ICTXT = ICTXTtoDC
      MYROW = MYROWtoDC
      MYCOL = MYCOLtoDC


! blacs array registration
      nb = 48
      np = numroc( n, nb, myrow, 0, nprow )
      nq = numroc( n, nb, mycol, 0, npcol )
      lddz = (n-1)/nprow+1
      lddz = ((lddz-1)/nb+1)*nb+1
      lddw = (n-1)/npcol+1
      lddw = ((lddw-1)/nb+1)*nb+1
      call descinit( descz, n, n, nb, nb, 0, 0, ictxt, lddz, info )
 
! preparing working arrays
      trilwmin = 3*n + max( nb*( np+1 ), 3*nb )
      lwork  = max( max( 1+6*n+2*np*nq, trilwmin ) + 2*n, lddz*lddw)
      liwork = 2+7*n+8*npcol
      allocate(work(lwork), iwork(liwork), stat=istat)
      if(istat /= 0) then
         print*,"memory exhausted"
         call flush(6)
         call mpi_abort( mpi_comm_eigen, 1, ierr )
      end if


! somehow, z must be nullified (originally next loop is not required.)
      z(1:lddz*lddw,1) = 0.0d+00

      call pdstedc( 'i', n, d(1), e(1), z(1,1), 1, 1, descz,
     &              work(1), lwork, iwork(1), liwork, info )

      call descinit( descw, n, n, 1, 1, 0, 0, ictxt, lddz, info )
      call pdgemr2d( n, n, z, 1, 1, descz, work, 1, 1, descw, ictxt )
      d1=mpi_wtime()
      call cyc2d_cyc1d(n,z,ldz,work,lddz, nprow, npcol, ierr)
      d2=mpi_wtime()
      call bcast_dbl( d(1), n, 1, mpi_comm_eigen )
      d3=mpi_wtime()
!      call eigen_free(2)

! freeing working arrays
      deallocate(work, iwork)

 
! blacs/pblas/scalapack finalize
!      call blacs_gridexit( ictxt )

#ifdef DETAIL
      if(myrank==1) then
         print*,"   bcast_dbl   :: ",d2-d1,"(sec)"
         print*,"   cyc2d_cyc1d :: ",d3-d2,"(sec)"
      endif
#endif
#ifdef TIMER
      s1=mpi_wtime()
      if(myrank==1)then
         print*,"Exectime of \"eigen_dch\" routine =",s1-s0,"(sec)"
      endif
#endif

      return
      end subroutine

