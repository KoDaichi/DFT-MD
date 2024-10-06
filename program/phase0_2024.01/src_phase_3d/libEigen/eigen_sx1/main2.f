      implicit double precision (a-h,o-z)

      real(8), pointer ::
     &           a(:),z(:),w(:)
      integer, parameter :: NDIM = 2

!      real(8), parameter :: PAI = 3.141592653589794D0
      real(8), parameter :: ONE = 1D0
      real(8) :: PAI
!
      include 'mpif.h'
      include 'trd.h'
!
      PAI = 4*ATAN(ONE)
!Kro  call errset(320, 256, 256, 0, 0, 320)
!
      call mpi_init(ierr)
      MPI_COMM_EIGEN=MPI_COMM_WORLD
      call mpi_comm_rank(MPI_COMM_EIGEN,i$inod,ierr)
      call mpi_comm_size(MPI_COMM_EIGEN,i$nnod,ierr)
*
      if(i$inod.eq.0)then
         open(10,file='IN')
      endif
*
      DO

         if(i$inod.eq.0)then
            read(10,*) n,m, nall, mtype
         endif

         call mpi_bcast(n,1,mpi_integer,0,MPI_COMM_EIGEN,ierr)
         if ( n < 1 .or. n > 400000 ) EXIT
         call mpi_bcast(m,1,mpi_integer,0,MPI_COMM_EIGEN,ierr)
         call mpi_bcast(nall,1,mpi_integer,0,MPI_COMM_EIGEN,ierr)
         call mpi_bcast(mtype,1,mpi_integer,0,MPI_COMM_EIGEN,ierr)


         call datacast_init(NDIM)
         NPROW = x_nnod
         NPCOL = y_nnod
! BLACS/PBLAS/SCALAPACK initialization
         CALL BLACS_PINFO( IAM, NPROCS )
         if ( NPROCS < 1 ) then
!  MPI group setup
            NPROCS = i$nnod
            IAM = i$inod
            CALL BLACS_SETUP( IAM, NPROCS )
         end if

         nx = ((n-1)/NPROW+1)
         call CSTAB_get_optdim(nx, 6, 16*4, 16*4*2, nm)
         if(i$inod.eq.0)then
            PRINT*,"N=",n, "NM=",nm
         endif
         call datacast_free(0)


!         NB  = 48
         NB  = 64+32
         nmz = ((n-1)/NPROW+1)
         nmz = ((nmz-1)/NB+1)*NB+1
         nmw = ((n-1)/NPCOL+1)
         nmw = ((nmw-1)/NB+1)*NB+1

         larray = MAX(nmz,nm)*nmw
         allocate(
     &               a(larray), z(larray),
     &               w(n),
     &               stat=istat)
         if(istat.ne.0) then
            print*,"Memory exhausted"
            call flush(6)
            EXIT
         endif
*
         call mat_set(n, a(1), nm, mtype, NDIM)
*-
         call eigen_sx(n, a(1), nm, w(1), z(1), nm, m)
*-
         if(i$inod.eq.0)then
            ax=0D0
            do i=1,n
               j=n-i
               theta=PAI*(2*j+1)/(2*n+1)
               x=5D-1/(1D0-COS(theta))
               ax=MAX(ax,ABS(w(i)-x)/ABS(x))
            enddo
            print*,"max|w(i)-w(i).true|/|w|=",ax
         endif
         call flush(6)
*-
         if(IAND(nall,4)/=0)then
            call mat_set(n, a(1), nm, mtype, NDIM)
            call ev_test_2D(n, a(1), nm, w(1), z(1), nm, t)
         endif
*
         deallocate(a)
         deallocate(z)
         deallocate(w)

      ENDDO

      if(i$inod.eq.0)then
         close(10)
      endif
*
*
      call MPI_Finalize(ierr)
      end

      subroutine flush(i)
      end

