         subroutine eigen_s(n, a, lda, w, z, ldz, m)
         implicit NONE
*-
         integer, intent(in)    :: n, lda, ldz, m
         real(8), intent(inout) :: a(lda,*), w(*), z(ldz,*)
         real(8), pointer       :: d(:), e(:), e2(:)
*-
         real(8)                :: hs0, hs1, s0, s1, FPS
         integer                :: m0, my_rank, world_size, INFO, ierr,i
         include 'mpif.h'
         include 'commtxt.h'
*-
         allocate(d(1:n), e(1:n), e2(1:n))
*-
         hs0 = MPI_Wtime()
         call MPI_Comm_size(MPI_COMM_EIGEN, world_size, ierr)
         call MPI_Comm_rank(MPI_COMM_EIGEN, my_rank,    ierr)
*-
            call papi_timer_start()

         s0 = MPI_Wtime()

            call tred1(n, a(1,1), lda, d(1), e(1), m)

         s1 = MPI_Wtime()

            call papi_timer_stop()

         if(my_rank==0)then
            print*,"TRD-BLK",n," ",s1-s0," ",
     &               1*dble(n)**3*4/3/(s1-s0)*1D-9,"GFLOPS"
            print*,"TRD-BLK-INFO",n," ",m
         endif

            call papi_timer_report()

         call flush(6)
*-
            call papi_timer_start()

         s0 = MPI_Wtime()

            e2(1:n) = e(1:n)
            w(1:n)  = d(1:n)

            call dc(n, w(1), e2(2), z(1,1), ldz, INFO, FPS)

         s1 = MPI_Wtime()

            call papi_timer_stop()

         if(my_rank==0)then
            print*,"D&C",s1-s0,"ERRCODE=",INFO
         endif
            call papi_timer_report()

         call flush(6)
*-
            call papi_timer_start()

         s0 = MPI_Wtime()

            m0 = 64
            m0 = m
            m0 = 128
            call trbakwy(n, a(1,1), lda, z(1,1), ldz, e(1), m0)

         s1 = MPI_Wtime()

            call papi_timer_stop()

         if(my_rank==0)then
             print*,"TRDBAK",n,s1-s0,2e-9*dble(n)**3/(s1-s0),"GFLOPS"
         endif
            call papi_timer_report()

         call flush(6)
*-
         hs1 = MPI_Wtime()
*
         if(my_rank==0)then
            print*,"Total",hs1-hs0,10e-9*dble(n)**3/3/(hs1-hs0),"GFLOPS"
         endif
         call flush(6)
*-
         deallocate(d, e, e2)
*-
         end subroutine

