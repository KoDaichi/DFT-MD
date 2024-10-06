      subroutine dcx(n, d, e, lde, z, ldz, INFO)

      integer, intent(in)           :: n, ldz
      real(8), intent(inout)        :: d(1:n), e(1:2*lde)
      real(8), intent(out)          :: z(1:ldz,*)

! Parameters BLACS array descritor(the position of entry tags), etc
      INTEGER, PARAMETER            :: BLOCK_CYCLIC_2D = 1
      INTEGER, PARAMETER            :: DLEN_  = 9
      INTEGER, PARAMETER            :: DTYPE_ = 1
      INTEGER, PARAMETER            :: CTXT_  = 2
      INTEGER, PARAMETER            :: M_     = 3
      INTEGER, PARAMETER            :: N_     = 4
      INTEGER, PARAMETER            :: MB_    = 5
      INTEGER, PARAMETER            :: NB_    = 6
      INTEGER, PARAMETER            :: RSRC_  = 7
      INTEGER, PARAMETER            :: CSRC_  = 8
      INTEGER, PARAMETER            :: LLD_   = 9

      integer                       :: DESCZ( DLEN_ )
      integer                       :: DESCW( DLEN_ )
 
      integer                       :: world_size, my_rank
      integer                       :: TRILWMIN, LWORK, LIWORK
 
      real(8), pointer              :: work(:)
      integer, pointer              :: iwork(:)

      real(8)  :: d1,d2

      include 'mpif.h'
      include 'trd.h'
 
 
      d1 = MPI_Wtime()

#ifdef PHASE
      CALL BLACS_PINFO( IAM, NPROCS )
      call datacast_init(2)
      NPROW = x_nnod
      NPCOL = y_nnod
      ICTXT = ICTXTtoDC
      MYROW = MYROWtoDC
      MYCOL = MYCOLtoDC
!     call BLACS_GET( -1, 0, ICTXT )
!     call BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
#else
! BLACS/PBLAS/SCALAPACK initialization
      CALL BLACS_PINFO( IAM, NPROCS )
      if ( NPROCS < 1 ) then
!  MPI group setup
         call MPI_COMM_SIZE( MPI_COMM_EIGEN, NPROCS, ierr )
         call MPI_COMM_RANK( MPI_COMM_EIGEN, IAM, ierr )
         CALL BLACS_SETUP( IAM, NPROCS )
      end if
      call BLACS_GET( -1, 0, ICTXT )
!     call BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      call datacast_init(2)
      NPROW = x_nnod
      NPCOL = y_nnod
!     if (myrow<0.or.myrow>nprow.or.mycol<0.or.mycol>npcol)then
         call BLACS_GRIDINIT( ICTXT, 'Column-major', nprow, npcol)
         call BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
!     endif
#endif

      d2 = MPI_Wtime()


! BLACS array registration
!      NB = 48
      NB = 64+32
!      NB = 1
!      NB = 4
      NP = NUMROC( n, NB, MYROW, 0, NPROW )
      NQ = NUMROC( n, NB, MYCOL, 0, NPCOL )
      lddz = (n-1)/NPROW+1
      lddz = ((lddz-1)/NB+1)*NB
      lddw = (n-1)/NPCOL+1
      lddw = ((lddw-1)/NB+1)*NB
      call DESCINIT( DESCZ, n, n, NB, NB, 0, 0, ICTXT, lddz, INFO )

!
! preparing working arrays
      nx     = (N-1)/NPCOL+1
      LWORK  = MAX(1+6*N+2*NP*(NQ+MAX(NQ,NB)), lddz*lddw, ldz*nx)
      LIWORK = 2+7*n+8*NPCOL
      allocate(work(lwork), iwork(liwork), stat=istat)
      if(istat.ne.0) then
           print*,"Memory exhausted"
           call flush(6)
           call MPI_Abort( MPI_COMM_EIGEN, 1, ierr )
      end if


! Somehow, Z must be nullified (Originally next loop is not required.)
      z(1:lddz*lddw,1) = 0.0D+00


!     MKL_MODE = MKL_GET_DYNAMIC()
!     call MKL_SET_DYNAMIC(0)

      d2 = MPI_Wtime()
      IF(IAM==0)print*,"before PDSTEDC",d2-d1



      d1 = MPI_Wtime ()
      CALL MY_PDSXEDC('I',2, n, d(1), e(1), lde, z(1,1), 1, 1, DESCZ,
     $               WORK(1), LWORK, IWORK(1), LIWORK, INFO)

      d2 = MPI_Wtime()
!     call MKL_SET_DYNAMIC(MKL_MODE)

      IF(IAM==0)print*,"PDSTEDC",d2-d1

      IF(NB==1)THEN
         do i=nx,1,-1
            work(1:lddz)=z(1+(i-1)*lddz:lddz+(i-1)*lddz,1)
            z(1:lddz,i)=work(1:lddz)
         enddo
         do i=1,nx
            z(lddz+1:ldz,i)=0.0D0
         enddo
      ELSE
         call DESCINIT( DESCW, n, n, 1, 1, 0, 0, ICTXT, ldz, INFO )
         call PDGEMR2D( n, n, z, 1, 1, DESCZ, work, 1, 1, DESCW, ICTXT )
         z(1:ldz*nx,1)=work(1:ldz*nx)
      ENDIF

      call bcast_dbl( d(1), n, 1, TRD_COMM_WORLD )
      call datacast_free(1)

! freeing working arrays
      deallocate(work)
      deallocate(iwork)

 
! BLACS/PBLAS/SCALAPACK finalize
!      call BLACS_GRIDEXIT( ICTXT )

      return
      end subroutine

