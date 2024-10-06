         subroutine eigen_sx(n, a, lda, w, z, ldz, m)
         implicit NONE
*-
         integer, intent(in)    :: n, lda, ldz, m
         real(8), intent(inout) :: a(lda,*), w(*), z(ldz,*)
         real(8), pointer       :: d(:), e(:), e2(:)
*-
         real(8)                :: hs0, hs1, s0, s1
         integer                :: m0, nb, nme
         integer                :: my_rank, world_size, INFO, ierr

         include 'mpif.h'
         include 'trd.h'
*-
!Kro
         call errset(320, 256, 256, 0, 0, 320)
         call errset(323, 256, 256, 0, 0, 323)
         hs0 = MPI_Wtime()
*-
         nme = ((n-1)/2+1)*2
         allocate(d(1:n), e(1:nme*2), e2(1:nme*2))
*-
#IFDEF PHASE
!        call datacast_init(2)
#ENDIF
         world_size = TRD_nnod
         my_rank    = TRD_inod-1
*-
         s0 = MPI_Wtime()

            call tred1(n, a(1,1), lda, d(1), e(1), nme, m)

         s1 = MPI_Wtime()

         if(my_rank==0)then
            print*,"TRD-BLK",n," ",s1-s0," ",
     &               1*dble(n)**3*4/3/(s1-s0)*1D-9,"GFLOPS"
            print*,"TRD-BLK-INFO",n," ",m
         endif
         call flush(6)
*-
         s0 = MPI_Wtime()

            e2(0*nme+1:0*nme+N-1) = e(0*nme+2:0*nme+N)
            e2(0*nme+N) = 0
            e2(1*nme+1:1*nme+N-2) = e(1*nme+3:1*nme+N)
            e2(1*nme+N-1) = 0
            e2(1*nme+N) = 0
            w(1:n)=d(1:n)

            call dcx(n, w(1), e2(1), nme, z(1,1), ldz, INFO)

         s1 = MPI_Wtime()

         if(my_rank==0)then
            print*,"D&C",s1-s0,"ERRCODE=",INFO
         endif
         call flush(6)
*-
         s0 = MPI_Wtime()

            m0 = 128
            nb = 1
            call trbakwy(n, a(1,1), lda, z(1,1), ldz,
     $                   e(1+nme), m0)

         s1 = MPI_Wtime()

         if(my_rank==0)then
             print*,"TRDBAK",n,s1-s0,2e-9*dble(n)**3/(s1-s0),"GFLOPS"
         endif
         call flush(6)
*-
         hs1 = MPI_Wtime()
*
         if(my_rank==0)then
            print*,"Total",hs1-hs0,10e-9*dble(n)**3/3/(hs1-hs0),"GFLOPS"
         endif
         call flush(6)
*-
#IFDEF PHASE
!        call datacast_free(0)
#ENDIF
         deallocate(d, e, e2)
*-
         end subroutine

