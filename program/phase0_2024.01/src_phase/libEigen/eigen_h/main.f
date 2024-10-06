      implicit double precision (a-h,o-z)

      real(8), allocatable :: a(:,:),ai(:,:)
      real(8), pointer :: b(:),bi(:),c(:),ci(:),w(:)
!
      include 'mpif.h'
!
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_eigen,i$inod,ierr)
      call mpi_comm_size(mpi_comm_eigen,i$nnod,ierr)
*
      n=1000
      m=32
      nm=(n/2)*2+1
      allocate( a(n,n), ai(n,n) )
      allocate(
     &         b((nm)*((n-1)/i$nnod+1+nm+1)),
     &         bi((nm)*((n-1)/i$nnod+1+nm+1)),
     &         c((nm)*((n-1)/i$nnod+1+nm+1)),
     &         ci((nm)*((n-1)/i$nnod+1+nm+1)),
     &         w(nm),
     &         stat=istat )
      if(istat.ne.0) then
         print*,"Memory exhausted"
         call flush(6)
         stop
      endif
*
      call matrix_set(a,n)
      call matrix_seti(ai,n)
      call matrix_adjust_h(a,n,b,nm)
      call matrix_adjust_hi(ai,n,bi,nm)
      call eigen_h(n,b,bi,w,nm,m,0)
*
#ifdef CHECK
      call matrix_set(a,n)
      call matrix_seti(ai,n)
      call matrix_adjust_h(a,n,c,nm)
      call matrix_adjust_hi(ai,n,ci,nm)
      call ev_test2(c,ci,w,b,bi,n,nm)
#endif
*
      deallocate(a,ai,b,bi,c,ci,w)

      call MPI_Finalize(ierr)
      end

