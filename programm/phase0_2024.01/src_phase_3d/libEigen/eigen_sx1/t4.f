      subroutine tred1_u(
     &              a, nm,
     &              u_x, u_y, nv,
     &              u_t, v_t, i, c, e, ne)

      implicit NONE

      integer, intent(in)    :: nm, nv, ne, i
      real(8), intent(inout) :: a(1:nm,*)
      real(8), intent(out)   :: u_x(1:nv,*), u_y(1:nv,*)
      real(8), intent(out)   :: u_t(*), v_t(*)

      include 'trd.h'

      real(8), intent(out)   :: c(MBAND,MBAND), e(1:ne,*)

      real(8)                :: anorm2, g_n(MBAND), a_n
      real(8)                :: t, tt(2*MBAND**2), ss(2*MBAND**2)
      real(8)                :: alpha(MBAND), beta(MBAND), scale(MBAND)

      real(8), parameter     :: ZERO = 0.0D+00, ONE = 1.0D+00

      integer                :: x_owner_nod, x_pos
      integer                :: y_owner_nod, y_pos
      integer                :: i_1, i_2, i_3, i_4, i_5
      integer                :: j_1, j_2, j_3, j_4, j_5
      integer                :: k_1, k_2, k_3, k_4, k_5
      integer                :: L, LL

      real(8) :: t11, t12, t22, tol
      real(8) :: a1, a2
      real(8) :: aL(MBAND,MBAND), aX(MBAND, MBAND)
      real(8), parameter :: EPS = 2D-16


      L  = i-MBAND
      LL = L-MBAND

      alpha(1:MBAND) = ONE
      beta (1:MBAND) = ONE
      g_n  (1:MBAND) = ZERO
      scale(1:MBAND) = ZERO

      c(1:MBAND, 1:MBAND)  = ZERO


      do i_1=MBAND,1,-1
      if ( LL + i_1 >= 1  ) then

         y_owner_nod = loop_c_node     (L+i_1, y_nnod, y_inod)
         x_pos       = loop_c_g2l_depth(i-1,   x_nnod, x_inod)
         call bcast_dbl(a(1,i_1), x_pos, y_owner_nod, y_COMM_WORLD)

      end if
      end do


      aL(1:MBAND, 1:MBAND) = ZERO

      do i_1=MBAND,1,-1
      if ( LL + i_1 >= 1  ) then

         x_owner_nod = loop_c_node     (LL+i_1, x_nnod, x_inod)
         x_pos       = loop_c_g2l_depth(LL+i_1, x_nnod, x_inod)

         do j_1=MBAND,1,-1
            if ( x_inod == x_owner_nod ) then
               aL(i_1,j_1) = a(x_pos, j_1)
            endif
         end do

      end if
      end do

      j_2 = loop_c_start(1, x_nnod, x_inod)
      j_3 = loop_c_end  (L, x_nnod, x_inod)

      t11 = ZERO; t12 = ZERO; t22 = ZERO
      do j_1=j_2,j_3
         a2  =a(j_1, 2)
         a1  =a(j_1, 1)
         t22 = t22 + a2 * a2
         t12 = t12 + a1 * a2
         t11 = t11 + a1 * a1
      end do

      tt(1) =  t22
      tt(2) =  t12
      tt(3) =  t11
      tt(4) =  aL(2,2)
      tt(5) =  aL(1,2)
      tt(6) =  aL(2,1)
      tt(7) =  aL(1,1)
      call reduce_dbl(tt, ss, 7, 1, x_COMM_WORLD)
      aX(2,2) = tt(1)
      aX(1,2) = tt(2)
      aX(2,1) = aX(1,2)
      aX(1,1) = tt(3)
      aL(2,2) = tt(4)
      aL(1,2) = tt(5)
      aL(2,1) = tt(6)
      aL(1,1) = tt(7)

      do i_1=MBAND,1,-1
      if ( LL + i_1 >= 1  ) then

         j_2 = loop_c_start(1,      x_nnod, x_inod)
         j_3 = loop_c_end  (LL+i_1, x_nnod, x_inod)

         x_owner_nod = loop_c_node     (LL+i_1, x_nnod, x_inod)
         x_pos       = loop_c_g2l_depth(LL+i_1, x_nnod, x_inod)

         if ( i_1 == 2 ) then

            alpha(i_1) = aX(i_1, i_1)
            a_n        = aL(i_1, i_1)

            do j_1=j_2,j_3
               u_x(j_1, i_1) = a(j_1, i_1)
            end do
         else
            tol = 16 * aX(1, 1) * EPS
            aX(1, 1) = MAX(ZERO, aX(1, 1) - aL(2, 1)**2)
            aX(1, 2) =           aX(1, 2) - aL(2, 1)*aL(2, 2)

            alpha(i_1) = aX(i_1, i_1)
            a_n        = aL(i_1, i_1)

! Singularity found, so re-calculate norm and inner products
            if ( alpha(i_1) < tol ) then
               do j_1=j_2,j_3
                  u_x(j_1, 1) = a(j_1, 1)
               end do
               do k_1=1,MBAND
                  t = ZERO
                  do j_1=j_2,j_3
                     t = t + u_x(j_1, i_1)*u_x(j_1, k_1)
                  enddo
                  tt(1+k_1)=t
               end do
               if ( x_inod == x_owner_nod ) then
                  a_n = a(x_pos, i_1)
               else
                  a_n = ZERO
               end if
               tt(1) =  a_n
               call reduce_dbl(tt, ss, 1+MBAND, 1, x_COMM_WORLD)
               a_n        = tt(1)
               do k_1=1,MBAND
                  aX(i_1,k_1) = tt(1+k_1)
               end do
               do k_1=1,MBAND
                  if ( k_1 /= i_1 ) then
                     aX(k_1,i_1) = aX(i_1,k_1)
                  end if
               end do
               alpha(i_1) = aX(i_1,i_1)
            else
               do j_1=j_2,j_3
                  u_x(j_1, i_1) = a(j_1, i_1)
               end do
            end if

         endif

         if ( alpha(i_1) > ZERO ) then

            scale(i_1) = SQRT(alpha(i_1))
!
! scaling for evading numerical instability
!
            t = ONE / scale(i_1)
            do j_1=j_2,j_3
               u_x(j_1, i_1) = u_x(j_1, i_1) * t
            end do
            do k_1 = 1, MBAND
               aL (k_1, i_1) = aL (k_1, i_1) * t
               aX (k_1, i_1) = aX (k_1, i_1) * t
               aX (i_1, k_1) = aX (i_1, k_1) * t
            end do

            a_n          =  a_n * t
            g_n(i_1)     = -SIGN(ONE, a_n)
            beta(i_1)    =  ONE - a_n * g_n(i_1)
            aL(i_1, i_1) =  a_n - g_n(i_1)

            if ( x_inod == x_owner_nod ) then
               u_x(x_pos, i_1) =  aL(i_1, i_1)
            end if

            do k_1=1,i_1-1
               aX(k_1, i_1) = aX(k_1, i_1) - g_n(i_1) * aL(i_1, k_1)
               aX(i_1, k_1) = aX(k_1, i_1)
            end do
               aX(i_1, i_1) = 2 * beta(i_1)
            do k_1=i_1+1, MBAND
               aX(i_1, k_1) = aX(i_1, k_1) - g_n(i_1) * aL(i_1, k_1)
               aX(k_1, i_1) = aX(i_1, k_1)
            end do

            do k_1=1,i_1-1
               t            = aX(k_1, i_1) / beta(i_1)
               j_2 = loop_c_start(1,      x_nnod, x_inod)
               j_3 = loop_c_end  (LL+i_1, x_nnod, x_inod)
               do j_1=j_2,j_3
                  a (j_1, k_1) = a (j_1, k_1) - t * u_x(j_1, i_1)
               end do! j_1
               do j_1=1,i_1
                  aL(j_1, k_1) = aL(j_1, k_1) - t * aL (j_1, i_1)
               end do
               aX(k_1, i_1) = - aX(k_1, i_1)
               aX(i_1, k_1) =   aX(k_1, i_1)
            end do

         end if

      end if
      end do
!
      do i_1=1,MBAND
      if ( LL + i_1 >= 1  ) then
!
! scaling back
!
         if ( alpha(i_1) > ZERO ) then
            g_n(i_1)   = g_n(i_1)  * scale(i_1)
            beta(i_1)  = beta(i_1) * alpha(i_1)
         else
            g_n(i_1)   = ZERO
            beta(i_1)  = ONE
         endif
         j_2 = loop_c_start(1,      x_nnod, x_inod)
         j_3 = loop_c_end  (LL+i_1, x_nnod, x_inod)
         t = scale(i_1)
         do j_1=j_2,j_3
            u_x(j_1, i_1) = u_x(j_1, i_1) * t
         end do
      end if
      end do
!
      y_owner_nod = loop_c_node(i, y_nnod, y_inod)
      if ( y_owner_nod == y_inod ) then
         if ( i > 2 ) then
            x_owner_nod = loop_c_node(i-2, x_nnod, x_inod)
            if ( x_owner_nod == x_inod ) then
               x_pos = loop_c_g2l_depth(i-2, x_nnod, x_inod)
               e(i-1, 1) = a(x_pos, 1)
            end if
         end if
         if ( i > 1 ) then
            x_owner_nod = loop_c_node(i-1, x_nnod, x_inod)
            if ( x_owner_nod == x_inod ) then
               x_pos = loop_c_g2l_depth(i-1, x_nnod, x_inod)
               e(i, 1) = a(x_pos, 2)
            end if
         end if
         x_owner_nod = loop_c_node(i,   x_nnod, x_inod)
         if ( x_owner_nod == x_inod ) then
            e(i-1, 2) = g_n(1)
            e(i  , 2) = g_n(2)
         end if
      end if

      do i_1=1,MBAND
      if ( LL + i_1 >= 1  ) then
         j_2 = loop_c_start(LL+i_1+1, x_nnod, x_inod)
         j_3 = loop_c_end  (i-1,      x_nnod, x_inod)
         u_x(j_2:j_3,i_1) = ZERO

         j_3 = loop_c_end  (LL+i_1+MBAND-1, x_nnod, x_inod)
         a(1:j_3,i_1) = u_x(1:j_3,i_1)

!         x_pos = loop_c_g2l_depth(i-1, x_nnod, x_inod)
!         call datacast_dbl(u_y(1, i_1), u_x(1, i_1), u_t, v_t, x_pos)
      end if
      end do
         x_pos = loop_c_g2l_depth(i-1, x_nnod, x_inod)
      if ( MOD(MBAND,2)==1 ) then; i_1 = 1
         call datacast_dbl(u_y(1, i_1), u_x(1, i_1),   u_t, v_t, x_pos)
      end if
      do i_1=MOD(MBAND,2)+1,MBAND,2
         call datacast_dbl2(u_y(1, i_1),u_y(1, i_1+1),
     &                      u_x(1, i_1),u_x(1, i_1+1), u_t, v_t, x_pos)
      end do

      c(1:MBAND, 1:MBAND) = ZERO
      do i_1=1,MBAND
         c(i_1, i_1) = ONE / beta(i_1)
      end do
      do i_1=1,MBAND
      do j_1=1,i_1-1
         c(j_1,i_1) = aX(j_1,i_1)*scale(j_1)*scale(i_1)
      end do
      end do

      return
      end subroutine tred1_u

