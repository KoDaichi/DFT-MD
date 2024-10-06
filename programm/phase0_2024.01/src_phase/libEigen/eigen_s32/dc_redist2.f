      subroutine dc_redist2(n, NB, a, lda, b, ldb, wk, lwk)
      implicit NONE

      integer, intent(in)    :: n, NB, lda, ldb, lwk
      real(8), intent(in)    :: a(lda,*)
      real(8), intent(out)   :: b(ldb,*)
      real(8), intent(inout) :: wk(lwk)

      real(8), pointer       :: wk1(:), wk2(:)
      integer, pointer       :: ir_sz(:,:)
      integer, pointer       :: l_(:)

      integer :: i,i0,j,k,k0,l,lx,m,m0
      integer :: iblk, jblk
      integer :: idist, ir_size, is_size
      integer :: his_rank, her_rank

      integer            :: PACK, MPACK

      include 'trd.h'
      include 'mpif.h'

      real(8) :: d1, d2, aa(10)


      aa = 0

      iblk = (n-1)/y_nnod+1
      jblk = (n-1)/x_nnod+1
      jblk = (jblk-1)/NB+1


      d1 = MPI_Wtime()
      allocate ( ir_sz(x_nnod,jblk) )
      d2 = MPI_Wtime()
      aa(1) = aa(1) + (d2-d1)

      d1=MPI_Wtime()

      do m = 4,1,-1
         lx = 0
         PACK = m
         do i0=1,jblk,PACK
         do idist = 1, x_nnod-1

            her_rank = MOD(x_inod-1+idist+x_nnod,x_nnod)+1
            his_rank = MOD(x_inod-1-idist+x_nnod,x_nnod)+1

            l = 0
            do i=i0,MIN(jblk,i0+PACK-1)
            do j=1,NB
              K=j+((x_inod-1)+(i-1)*x_nnod)*NB
              if ( K <= n .AND. MOD(K-1,x_nnod)+1 == her_rank ) then
                 l = l + 1
              end if
            end do
            end do

            lx = MAX(l, lx)

            l = 0
            do i=i0,MIN(jblk,i0+PACK-1)
            do j=1,NB
              K=j+((his_rank-1)+(i-1)*x_nnod)*NB
              if ( K <= n .AND. MOD(K-1,x_nnod)+1 == x_inod ) then
                 l = l + 1
              end if
            end do
            end do
            ir_sz(idist,i0) = l

            lx = MAX(l, lx)

         end do
         end do
         lx = MAX(1,lx)
         MPACK = MAX(1,lwk/(2*lx))
         MPACK = MIN(MPACK, iblk)
         if ( 2*lx*MPACK <= lwk ) EXIT
      enddo

      d2 = MPI_Wtime()
      aa(3) = aa(3) + (d2-d1)


      MPACK = MAX(1,lwk/(2*lx))
      MPACK = MIN(MPACK, iblk)
      if ( 2*lx*MPACK <= lwk ) then
         call dc_redist2_sub(n, NB, a, lda, b, ldb,
     &                       wk(1), wk(1+lx*MPACK),
     &                       ir_sz, x_nnod, PACK, MPACK, aa)
      else
      MPACK = MAX(1,lwk/(1*lx))
      MPACK = MIN(MPACK, iblk)
      if ( lx*MPACK <= lwk ) then
         d1 = MPI_Wtime()
         allocate ( wk2(1:lx*MPACK) )
         d2 = MPI_Wtime()
         aa(1) = aa(1) + (d2-d1)
         call dc_redist2_sub(n, NB, a, lda, b, ldb,
     &                       wk(1), wk2(1),
     &                       ir_sz, x_nnod, PACK, MPACK, aa)
         d1 = MPI_Wtime()
         deallocate ( wk1 )
         deallocate ( wk2 )
         d2 = MPI_Wtime()
         aa(1) = aa(1) + (d2-d1)
      else
      MPACK = MAX(1,lwk/(2*lx))
      MPACK = MIN(MPACK, iblk)
         d1 = MPI_Wtime()
         allocate ( wk1(1:lx*MPACK) )
         allocate ( wk2(1:lx*MPACK) )
         d2 = MPI_Wtime()
         aa(1) = aa(1) + (d2-d1)
         call dc_redist2_sub(n, NB, a, lda, b, ldb,
     &                       wk1(1), wk2(1),
     &                       ir_sz, x_nnod, PACK, MPACK, aa)
         d1 = MPI_Wtime()
         deallocate ( wk1 )
         deallocate ( wk2 )
         d2 = MPI_Wtime()
         aa(1) = aa(1) + (d2-d1)
      endif
      endif
      
      d1 = MPI_Wtime()
      deallocate ( ir_sz )
      d2 = MPI_Wtime()
      aa(1) = aa(1) + (d2-d1)

      if ( TRD_inod == 1 ) print *, PACK, MPACK, aa(1:6)


      end subroutine

      subroutine dc_redist2_sub(n, NB, a, lda, b, ldb,
     &                      wk1, wk2, ir_sz, ldr, PACK, MPACK, aa)
      implicit NONE

      integer, intent(in)    :: n, NB, lda, ldb, ldr, PACK, MPACK
      real(8), intent(in)    :: a(lda,*)
      real(8), intent(out)   :: b(ldb,*)
      real(8), intent(inout) :: wk1(*), wk2(*)
      integer, intent(in)    :: ir_sz(ldr,*)
      real(8), intent(inout) :: aa(*)

      integer, pointer       :: l_(:)

      integer :: i,i0,j,k,k0,l,lx,m,m0
      integer :: iblk, jblk
      integer :: idist, ir_size, is_size
      integer :: iq(2), his_rank, her_rank

      include 'trd.h'
      include 'mpif.h'

      real(8) :: d1, d2


      iblk = (n-1)/y_nnod+1
      jblk = (n-1)/x_nnod+1
      jblk = (jblk-1)/NB+1


      allocate ( l_(1:iblk+1) )

      do i0=1,jblk,PACK

            d1 = MPI_Wtime()
!$OMP PARALLEL DO private(i,j,K,K0,m)
      do m=1,iblk
         do i=i0,MIN(jblk,i0+PACK-1)
         do j=1,NB
            K=j+((x_inod-1)+(i-1)*x_nnod)*NB
            if ( K <= n .AND. MOD(K-1,x_nnod)+1 == x_inod ) then
               K0=(K-1)/x_nnod+1
               b(K0,m) = a(j+NB*(i-1),m)
            end if
         end do
         end do
      end do
!$OMP END PARALLEL DO
            d2 = MPI_Wtime()
            aa(2) = aa(2) + (d2-d1)

      do idist = 1, x_nnod-1

         her_rank = MOD(x_inod-1+idist+x_nnod,x_nnod)+1
         his_rank = MOD(x_inod-1-idist+x_nnod,x_nnod)+1


      do m0=1,iblk,MPACK

         ir_size=ir_sz(idist,i0)*MIN(iblk-m0+1,MPACK)
         if ( ir_size > 0 ) then

            d1 = MPI_Wtime()
            call irecv_dbl(wk2, ir_size, his_rank, iq(2), x_COMM_WORLD)
            d2 = MPI_Wtime()
            aa(4) = aa(4) + (d2-d1)

         end if


         l_(1)=0
!$OMP PARALLEL DO private(i,j,K,l,m)
         do m=m0,MIN(iblk,m0+MPACK-1)
            l = 0
            do i=i0,MIN(jblk,i0+PACK-1)
            do j=1,NB
               K=j+((x_inod-1)+(i-1)*x_nnod)*NB
               if ( K <= n .AND. MOD(K-1,x_nnod)+1 == her_rank ) then
                  l = l + 1
               end if
            end do
            end do
            l_(m-m0+1+1)=l
         end do
!$OMP END PARALLEL DO
         do m=m0+1,MIN(iblk,m0+MPACK-1)+1
            l=l_(m-m0+1)+l_(m-m0+1-1)
            l_(m-m0+1)=l
         end do


         is_size=l
         if ( is_size > 0 ) then

            d1 = MPI_Wtime()
!$OMP PARALLEL DO private(i,j,K,l,m)
            do m=m0,MIN(iblk,m0+MPACK-1)
               l=l_(m-m0+1)
               do i=i0,MIN(jblk,i0+PACK-1)
               do j=1,NB
                  K=j+((x_inod-1)+(i-1)*x_nnod)*NB
                  if ( K <= n .AND. MOD(K-1,x_nnod)+1 == her_rank ) then
                     l = l + 1
                     wk1(l) = a(j+NB*(i-1),m)
                  end if
               end do
               end do
            end do
!$OMP END PARALLEL DO
            d2 = MPI_Wtime()
            aa(5) = aa(5) + (d2-d1)

            d1 = MPI_Wtime()
            call isend_dbl(wk1, is_size, her_rank, iq(1), x_COMM_WORLD)
            d2 = MPI_Wtime()
            aa(4) = aa(4) + (d2-d1)

         end if


         ir_size=ir_sz(idist,i0)*MIN(iblk-m0+1,MPACK)
         if ( ir_size > 0 ) then

            l_(1)=0
!$OMP PARALLEL DO private(i,j,K,K0,l,m)
            do m=m0,MIN(iblk,m0+MPACK-1)
               l = 0
               do i=i0,MIN(jblk,i0+PACK-1)
               do j=1,NB
                  K=j+((his_rank-1)+(i-1)*x_nnod)*NB
                  if ( K <= n .AND. MOD(K-1,x_nnod)+1 == x_inod ) then
                     l = l + 1
                  end if
               end do
               end do
               l_(m-m0+1+1)=l
            end do
!$OMP END PARALLEL DO
            do m=m0+1,MIN(iblk,m0+MPACK-1)+1
               l=l_(m-m0+1)+l_(m-m0+1-1)
               l_(m-m0+1)=l
            end do

            d1 = MPI_Wtime()
            call wait_dbl(iq(2))
            d2 = MPI_Wtime()
            aa(4) = aa(4) + (d2-d1)

            d1 = MPI_Wtime()
!$OMP PARALLEL DO private(i,j,K,K0,l,m)
            do m=m0,MIN(iblk,m0+MPACK-1)
               l=l_(m-m0+1)
               do i=i0,MIN(jblk,i0+PACK-1)
               do j=1,NB
                  K=j+((his_rank-1)+(i-1)*x_nnod)*NB
                  if ( K <= n .AND. MOD(K-1,x_nnod)+1 == x_inod ) then
                     l = l + 1
                     K0=(K-1)/x_nnod+1
                     b(K0,m) = wk2(l)
                  end if
               end do
               end do
            end do
!$OMP END PARALLEL DO
            d2 = MPI_Wtime()
            aa(6) = aa(6) + (d2-d1)

         end if

         if ( is_size > 0 ) then

            d1 = MPI_Wtime()
            call wait_dbl(iq(1))
            d2 = MPI_Wtime()
            aa(4) = aa(4) + (d2-d1)

         end if

      end do

      end do

      end do

      d1 = MPI_Wtime()
      deallocate ( l_ )
      d2 = MPI_Wtime()
      aa(1) = aa(1) + (d2-d1)

      return
      end subroutine

