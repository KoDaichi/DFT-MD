      subroutine dc_redist1(n, NB, a, b, ldm, wk, lwk)
      implicit NONE

      integer, intent(in)  :: n, NB, ldm, lwk
      real(8), intent(in)  :: a(ldm, *)
      real(8), intent(out) :: b(ldm, *)
      real(8), intent(inout) :: wk(*)

      real(8), pointer     :: wk1(:, :), wk2(:, :)

      integer :: i,j,j0,k,k0,l,lx
      integer :: iNQ
      integer :: iblk, jblk
      integer :: idist, ir_size, is_size
      integer :: his_rank, her_rank
      integer, pointer :: ir_sz(:, :, :)

      integer :: NQ, NBQ

      include 'trd.h'


      iblk = (n-1)/y_nnod+1
      iblk = (iblk-1)/NB+1
      jblk = (n-1)/x_nnod+1
      jblk = (iblk-1)/NB+1


      do iNQ=1,NB

         NQ = iNQ
         NBQ = (NB-1)/NQ+1

         lx = 0
         do i=1,iblk
         do j0=1,NB,NBQ
         do idist = 1, y_nnod-1

            his_rank = MOD(y_inod-1-idist+y_nnod,y_nnod)+1
            her_rank = MOD(y_inod-1+idist+y_nnod,y_nnod)+1

            l = 0
            do j=j0,MIN(NB,j0+NBQ-1)
              K=j+((y_inod-1)+(i-1)*y_nnod)*NB
              if ( K <= n .AND. MOD(K-1,y_nnod)+1 == her_rank ) then
                 l = l + 1
              end if
            end do
            lx = MAX(l, lx)

            l = 0
            do j=j0,MIN(NB,j0+NBQ-1)
              K=j+((his_rank-1)+(i-1)*y_nnod)*NB
              if ( K <= n .AND. MOD(K-1,y_nnod)+1 == y_inod ) then
                 l = l + 1
              end if
            end do
            lx = MAX(l, lx)

            if ( lwk >= y_nnod*iblk*NQ ) then
               K = idist + ((i-1) + ((j0-1)/NBQ) * iblk) * y_nnod
               wk(K) = l
            end if

         end do
         end do
         end do

         if ( ldm*lx*2 <= lwk ) then
            exit
         end if

      end do

      allocate ( ir_sz(y_nnod, iblk, NQ) )

      if ( lwk >= y_nnod*iblk*NQ ) then

         do i=1,iblk
         do j0=1,NB,NBQ
         do idist = 1, y_nnod-1
            K = idist + ((i-1) + ((j0-1)/NBQ) * iblk) * y_nnod
            ir_sz(idist, i, (j0-1)/NBQ+1) = wk(K)
         enddo
         enddo
         enddo

      else

         do i=1,iblk
         do j0=1,NB,NBQ
         do idist = 1, y_nnod-1

            his_rank = MOD(y_inod-1-idist+y_nnod,y_nnod)+1
            l = 0
            do j=j0,MIN(NB,j0+NBQ-1)
              K=j+((his_rank-1)+(i-1)*y_nnod)*NB
              if ( K <= n .AND. MOD(K-1,y_nnod)+1 == y_inod ) then
                 l = l + 1
              end if
            end do
            ir_sz(idist,i,(j0-1)/NBQ+1) = l

         end do
         end do
         end do

      end if

      if ( ldm*lx*2 <= lwk ) then
         call dc_redist1_sub ( n, NB, a, b, ldm,
     &                       wk(1), wk(1+ldm*lx),
     &                       ir_sz, y_nnod, iblk, NBQ )
      else if ( ldm*lx <= lwk ) then
         allocate ( wk2(1:ldm, 1:lx) )
         call dc_redist1_sub ( n, NB, a, b, ldm,
     &                       wk, wk2(1,1),
     &                       ir_sz, y_nnod, iblk, NBQ )
         deallocate ( wk2 )
      else
         allocate ( wk1(1:ldm, 1:lx) )
         allocate ( wk2(1:ldm, 1:lx) )
         call dc_redist1_sub ( n, NB, a, b, ldm,
     &                       wk1(1,1), wk2(1,1),
     &                       ir_sz, y_nnod, iblk, NBQ )
         deallocate ( wk1 )
         deallocate ( wk2 )
      end if

      deallocate ( ir_sz )

      return
      end subroutine


      subroutine dc_redist1_sub(n, NB, a, b, ldm,
     $               wk1, wk2, ir_sz, l1, l2, NBQ)
      implicit NONE

      integer, intent(in)  :: n, NB, ldm, l1, l2, NBQ
      real(8), intent(in)  :: a(ldm,*)
      real(8), intent(out) :: b(ldm,*)
      real(8), intent(inout) :: wk1(ldm,*), wk2(ldm,*)
      integer, intent(in) :: ir_sz(l1,l2,*)

      integer :: i,j,j0,k,k0,l,lx
      integer :: iblk, jblk
      integer :: idist, ir_size, is_size
      integer :: iq_r, iq_s, his_rank, her_rank

      include 'trd.h'
      integer, pointer :: l_(:)

      iblk = (n-1)/y_nnod+1
      iblk = (iblk-1)/NB+1
      jblk = (n-1)/x_nnod+1
      jblk = (iblk-1)/NB+1


      allocate ( l_(0:NB) )


      do i=1,iblk
      do j0=1,NB,NBQ

!$OMP PARALLEL DO PRIVATE(j,K,K0,l)
         do j=j0,MIN(NB,j0+NBQ-1)
           K=j+((y_inod-1)+(i-1)*y_nnod)*NB
           if ( K <= n .AND. MOD(K-1,y_nnod)+1 == y_inod ) then
              K0=(K-1)/y_nnod+1
              b(1:ldm, K0) = a(1:ldm, j+NB*(i-1))
           end if
         end do
!$OMP END PARALLEL DO

      do idist = 1, y_nnod-1

         her_rank = MOD(y_inod-1+idist+y_nnod,y_nnod)+1
         his_rank = MOD(y_inod-1-idist+y_nnod,y_nnod)+1

         ir_size=ldm*ir_sz(idist,i,(j0-1)/NBQ+1)
         if ( ir_size > 0 ) then
         call irecv_dbl ( wk2, ir_size, his_rank, iq_r, y_COMM_WORLD )
         end if

         l = 0; l_(0)=0
         do j=j0,MIN(NB,j0+NBQ-1)
           K=j+((y_inod-1)+(i-1)*y_nnod)*NB
           if ( K <= n .AND. MOD(K-1,y_nnod)+1 == her_rank ) then
              l = l + 1
           end if
           l_(j-j0+1) = l
         end do
         lx = l
!$OMP PARALLEL DO PRIVATE(j,K,K0,l)
         do j=j0,MIN(NB,j0+NBQ-1)
           l = l_(j-j0)
           K=j+((y_inod-1)+(i-1)*y_nnod)*NB
           if ( K <= n .AND. MOD(K-1,y_nnod)+1 == her_rank ) then
              l = l + 1
              wk1(1:ldm, l) = a(1:ldm, j+NB*(i-1))
           end if
         end do
!$OMP END PARALLEL DO

         is_size=ldm*lx
         if ( is_size > 0 ) then
         call isend_dbl ( wk1, is_size, her_rank, iq_s, y_COMM_WORLD )
         end if

         if ( ir_size > 0 ) then
         call wait_dbl ( iq_r )

         l = 0; l_(0)=0
         do j=j0,MIN(NB,j0+NBQ-1)
           K=j+((his_rank-1)+(i-1)*y_nnod)*NB
           if ( K <= n .AND. MOD(K-1,y_nnod)+1 == y_inod ) then
              l = l + 1
           end if
           l_(j-j0+1) = l
         end do
         lx = l
!$OMP PARALLEL DO PRIVATE(j,K,K0,l)
         do j=j0,MIN(NB,j0+NBQ-1)
           l = l_(j-j0)
           K=j+((his_rank-1)+(i-1)*y_nnod)*NB
           if ( K <= n .AND. MOD(K-1,y_nnod)+1 == y_inod ) then
              l = l + 1
              K0=(K-1)/y_nnod+1
              b(1:ldm, K0) = wk2(1:ldm, l)
           end if
         end do
!$OMP END PARALLEL DO
         end if

         if ( is_size > 0 ) then
         call wait_dbl ( iq_s )
         end if

      enddo

      end do
      end do


      deallocate ( l_ )


      return
      end subroutine

