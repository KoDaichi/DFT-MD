      subroutine eigen_prd_vo(
     &              ux, vx, vy, nv,
     &              u_t, v_t, c,
     &              i, i_base, m0)
!$    use OMP_LIB
      use eigen_devel
      use eigen_libs
      implicit NONE
!
      integer, intent(in)    ::  nv, i, i_base, m0
      real(8), intent(inout) ::  ux(1:nv, *)
      real(8), intent(inout) ::  vx(1:nv, *), vy(1:nv, *)
      real(8), intent(inout) ::  u_t(4,*)
      real(8), intent(inout) ::  v_t(*)
!
      real(8), intent(inout) ::  c(1:2,1:2)
!
      real(8)                ::  g_g(1:2,1:2)
      real(8)                ::  cu_cv(1:2,1:2)
!
      integer                ::  j, k, L, n, LL
      integer                ::  k_1, k_2
!
      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  l_1, l_2, l_3, l_4
      integer                ::  jj_1, jj_2, jj_3, jj_4
      integer                ::  kk_1, kk_2, kk_3, kk_4
      integer                ::  LX
!
!-       integer, parameter     ::  VTOL = 2048
       integer, parameter     ::  VTOL = 1024
!===       integer, parameter     ::  VTOL = 204800
!
      include 'param.h'
      include 'CSTAB.h'
!
      real(8)                ::  w0, w1, w2, w3
      real(8)                ::  u0_0, v0_0
      real(8)                ::  u1_0, v1_0
      real(8)                ::  u0_1, v0_1
      real(8)                ::  u1_1, v1_1
      real(8)                ::  uu_x
      real(8)                ::  ux0, vx0
      real(8)                ::  ux1, vx1

      real(8)                ::  tt(4), ss(4)
      real(8)                ::  sx(2,2), tx(2,2)
      real(8)                ::  s11, s21, s12, s22
      real(8)                ::  c11, c21, c12, c22
      real(8)                ::  u0, u1, v0, v1, u12


      k_1 = i - i_base
      k_2 = m0

      L = i - MBAND
      n  = eigen_translate_g2l(L, x_nnod,x_inod)


! FOR attention to unexpected overflow or NAN
          j_3 = eigen_loop_end(L, x_nnod,x_inod)
          if ( j_3 < n ) then
          do i_1=1,MBAND
             vx(j_3+1:n,k_1-i_1+1) = ZERO ! in case
          end do
          end if


      l_4 = MOD(k_2-k_1,2)+k_1+1
      LX  = L1_LSIZE*L1_WAY/16

      if ( k_2 <= k_1 ) then

         u_t(1:MBAND*MBAND, 1:2)=ZERO

         jj_2 = 1
         jj_3 = n

         do jj_1=jj_2,jj_3,LX
            j_2 = jj_1; j_3 = MIN(jj_1+LX-1,jj_3)

            j = 1
            u0_0 = u_t(1, j+0)
            u0_1 = u_t(2, j+0)
            u1_0 = u_t(3, j+0)
            u1_1 = u_t(4, j+0)
            uu_x = u_t(1, j+1)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
            do j_1=j_2,j_3
               u0 = ux(j_1, k_1-0)
               u1 = ux(j_1, k_1-1)
               v0 = vx(j_1, k_1-0)
               v1 = vx(j_1, k_1-1)
               u0_0 = u0_0 + u1 * v1
               u0_1 = u0_1 + u1 * v0
               u1_0 = u1_0 + u0 * v1
               u1_1 = u1_1 + u0 * v0
               uu_x = uu_x + u0 * u1
            end do
            u_t(1, j+0) = u0_0
            u_t(2, j+0) = u0_1
            u_t(3, j+0) = u1_0
            u_t(4, j+0) = u1_1
            u_t(1, j+1) = uu_x

         end do

!$OMP BARRIER

         call reduce_dbl(u_t, v_t, 5, 1, x_COMM_WORLD)

         g_g(1,1) =  u_t(1,1)
         g_g(1,2) = (u_t(2,1) + u_t(3,1))/2
         g_g(2,2) =  u_t(4,1)
         g_g(2,1) =  u_t(1,2)

      else

!
! v=v-(UV+VU)u
!
         l_2 = k_2-k_1+2
         u_t(1:MBAND*MBAND, 1:l_2) = ZERO

         l_4 = MOD(k_2-k_1, 2)+k_1+1
         LX  = L1_LSIZE*L1_WAY/16

         LL = (n-1)/y_nnod+1
         LL = ((LL-1)/2+1)*2

         if ( n > VTOL ) then
            jj_2 = 1+LL*(y_inod-1)
            jj_3 = MIN(n, LL*y_inod)
         else
            jj_2 = 1
            jj_3 = n
         endif

         do jj_1=jj_2,jj_3,LX
            j_2 = jj_1; j_3 = MIN(jj_1+LX-1, jj_3)

            j = 1
            u0_0 = u_t(1, j+0)
            u0_1 = u_t(2, j+0)
            u1_0 = u_t(3, j+0)
            u1_1 = u_t(4, j+0)
            uu_x = u_t(1, j+1)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
            do j_1=j_2,j_3
               u0 = ux(j_1, k_1-0)
               u1 = ux(j_1, k_1-1)
               v0 = vx(j_1, k_1-0)
               v1 = vx(j_1, k_1-1)
               u0_0 = u0_0 + u1 * v1
               u0_1 = u0_1 + u1 * v0
               u1_0 = u1_0 + u0 * v1
               u1_1 = u1_1 + u0 * v0
               uu_x = uu_x + u0 * u1
            end do
            u_t(1, j+0) = u0_0
            u_t(2, j+0) = u0_1
            u_t(3, j+0) = u1_0
            u_t(4, j+0) = u1_1
            u_t(1, j+1) = uu_x

         end do

         do jj_1=jj_2,jj_3,LX
            j_2 = jj_1; j_3 = MIN(jj_1+LX-1, jj_3)

            do l_1=k_1+1,l_4-1

               j = l_1-k_1 +2
               u0_0 = u_t(1, j+0)
               v0_0 = u_t(2, j+0)
               u1_0 = u_t(3, j+0)
               v1_0 = u_t(4, j+0)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
               do j_1=j_2,j_3
                  w0 = ux(j_1, k_1-0)
                  w1 = ux(j_1, k_1-1)
                  u0_0 = u0_0 + vx(j_1, l_1+0) * w0
                  v0_0 = v0_0 + ux(j_1, l_1+0) * w0
                  u1_0 = u1_0 + vx(j_1, l_1+0) * w1
                  v1_0 = v1_0 + ux(j_1, l_1+0) * w1
               end do! j_1
               u_t(1, j+0) = u0_0
               u_t(2, j+0) = v0_0
               u_t(3, j+0) = u1_0
               u_t(4, j+0) = v1_0

            end do! l_1
            do l_1=l_4,k_2,2                      ! 1

               j = l_1-k_1 +2
               u0_0 = u_t(1, j+0)
               v0_0 = u_t(2, j+0)
               u1_0 = u_t(3, j+0)
               v1_0 = u_t(4, j+0)
               u0_1 = u_t(1, j+1)
               v0_1 = u_t(2, j+1)
               u1_1 = u_t(3, j+1)
               v1_1 = u_t(4, j+1)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
               do j_1=j_2,j_3
                  w0 = ux(j_1, k_1-0)
                  u0_0 = u0_0 + vx(j_1, l_1+0) * w0
                  v0_0 = v0_0 + ux(j_1, l_1+0) * w0
                  u0_1 = u0_1 + vx(j_1, l_1+1) * w0
                  v0_1 = v0_1 + ux(j_1, l_1+1) * w0
                  w1 = ux(j_1, k_1-1)
                  u1_0 = u1_0 + vx(j_1, l_1+0) * w1
                  v1_0 = v1_0 + ux(j_1, l_1+0) * w1
                  u1_1 = u1_1 + vx(j_1, l_1+1) * w1
                  v1_1 = v1_1 + ux(j_1, l_1+1) * w1
               end do! j_1
               u_t(1, j+0) = u0_0
               u_t(2, j+0) = v0_0
               u_t(3, j+0) = u1_0
               u_t(4, j+0) = v1_0
               u_t(1, j+1) = u0_1
               u_t(2, j+1) = v0_1
               u_t(3, j+1) = u1_1
               u_t(4, j+1) = v1_1

            end do! l_1

         end do! jj_1

!$OMP BARRIER

         l_2 = (k_2-k_1+2)*2*MBAND
         if ( n > VTOL ) then
            call reduce_dbl(u_t, v_t, l_2, 1, y_COMM_WORLD)
         end if
         call reduce_dbl(u_t, v_t, l_2, 1, x_COMM_WORLD)

         g_g(1,1) = u_t(1,1)
         g_g(1,2) = (u_t(2,1) + u_t(3,1))/2
         g_g(2,2) = u_t(4,1)
         g_g(2,1) = u_t(1,2)

         cu_cv(1:2,1:2)=ZERO

         do l_1=k_1+1,k_2
            j = l_1-k_1 +2

            u0_0 = u_t(1, j+0)
            v0_0 = u_t(2, j+0)
            u1_0 = u_t(3, j+0)
            v1_0 = u_t(4, j+0)

            cu_cv(1,1) = cu_cv(1,1) + v1_0*u1_0
            cu_cv(2,1) = cu_cv(2,1) + v1_0*u0_0
            cu_cv(1,2) = cu_cv(1,2) + v0_0*u1_0
            cu_cv(2,2) = cu_cv(2,2) + v0_0*u0_0
         end do

         g_g(1,1) = g_g(1,1) - 2*cu_cv(1,1)
         g_g(1,2) = g_g(1,2) -  (cu_cv(1,2) + cu_cv(2,1))
         g_g(2,2) = g_g(2,2) - 2*cu_cv(2,2)

         if ( n > VTOL ) then
            jj_2 = 1+LL*(y_inod-1)
            jj_3 = MIN(n, LL*y_inod)
         else
            jj_2 = 1
            jj_3 = n
         endif

         do jj_1=jj_2,jj_3,LX
            j_2 = jj_1; j_3 = MIN(jj_1+LX-1, jj_3)

            do l_1=k_1+1,l_4-1

               j = l_1-k_1 +2

               u0_0 = u_t(1, j+0)
               v0_0 = u_t(2, j+0)
               u1_0 = u_t(3, j+0)
               v1_0 = u_t(4, j+0)

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
               do j_1=j_2,j_3
                  w0  = ZERO ! vx(j_1, k_1-0)
                  w1  = ZERO ! vx(j_1, k_1-1)
                  ux0 = ux(j_1, l_1+0)
                  vx0 = vx(j_1, l_1+0)
                  w0 = w0
     &               - ux0 * u0_0
     &               - vx0 * v0_0
                  w1 = w1
     &               - ux0 * u1_0
     &               - vx0 * v1_0
                  vx(j_1, k_1-0) = vx(j_1, k_1-0) + w0
                  vx(j_1, k_1-1) = vx(j_1, k_1-1) + w1
               end do! j_1
            end do! l_1
            do l_1=l_4,k_2,2

               j = l_1-k_1 +2

               u0_0 = u_t(1, j+0)
               v0_0 = u_t(2, j+0)
               u1_0 = u_t(3, j+0)
               v1_0 = u_t(4, j+0)
               u0_1 = u_t(1, j+1)
               v0_1 = u_t(2, j+1)
               u1_1 = u_t(3, j+1)
               v1_1 = u_t(4, j+1)

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
               do j_1=j_2,j_3
                  w0  = ZERO ! vx(j_1, k_1-0)
                  w1  = ZERO ! vx(j_1, k_1-1)
                  ux0 = ux(j_1, l_1+0)
                  vx0 = vx(j_1, l_1+0)
                  w0 = w0
     &               - ux0 * u0_0
     &               - vx0 * v0_0
                  w1 = w1
     &               - ux0 * u1_0
     &               - vx0 * v1_0
                  ux1 = ux(j_1, l_1+1)
                  vx1 = vx(j_1, l_1+1)
                  w0 = w0
     &               - ux1 * u0_1
     &               - vx1 * v0_1
                  w1 = w1
     &               - ux1 * u1_1
     &               - vx1 * v1_1
                  vx(j_1, k_1-0) = vx(j_1, k_1-0) + w0
                  vx(j_1, k_1-1) = vx(j_1, k_1-1) + w1
               end do! j_1
            end do! l_1

         end do! jj_1

      end if

      ss(1)   = g_g(1,1)
      ss(2)   = g_g(1,2)
      ss(3)   = g_g(2,2)

      c11     = c(1,1)
      c21     = -c(2,2)*(c(1,1)*g_g(2,1))
      c12     = ZERO
      c22     = c(2,2)

! SX:=SS
      sx(1,1) = ss(1)
      sx(2,1) = ss(2)
      sx(1,2) = ss(2)
      sx(2,2) = ss(3)

! TX:=SX*C
      tx(1,1) = sx(1,1)*c11 + sx(1,2)*c21
      tx(2,1) = sx(2,1)*c11 + sx(2,2)*c21
      tx(1,2) =               sx(1,2)*c22
      tx(2,2) =               sx(2,2)*c22
            
! S:=C^T*TX
      s11     = c11*tx(1,1) + c21*tx(2,1)
      s21     =               c22*tx(2,1)
      s12     = c11*tx(1,2) + c21*tx(2,2)
      s22     =               c22*tx(2,2)

! S:=diag(S)+upper(S)/2
      s11     = (s11      )/2
      s21     = (s21 + s12)/2
      s12     = ZERO
      s22     = (s22      )/2

! V:=VC^T
!DIR$ UNROLL(2)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do j_1=jj_2,jj_3
         v0 = vx(j_1, k_1-1)
         v1 = vx(j_1, k_1-0)
         v0            =      v0 * c11 + v1 * c21
         v1            =                 v1 * c22
!
! V:=V-US
!
         u0 = ux(j_1, k_1-1)
         u1 = ux(j_1, k_1-0)
         vx(j_1, k_1-1) = v0 - s11 * u0 - s21 * u1
         vx(j_1, k_1-0) = v1            - s22 * u1
      end do! j_1

      if ( k_2 > k_1 ) then
      if ( n > VTOL ) then
         j_3 = eigen_loop_end(L, x_nnod, x_inod)
         call allgather_dbl(vx(jj_2, k_1-0), v_t, LL, y_COMM_WORLD)
         vx(1:j_3, k_1-0) = v_t(1:j_3)
         call allgather_dbl(vx(jj_2, k_1-1), v_t, LL, y_COMM_WORLD)
         vx(1:j_3, k_1-1) = v_t(1:j_3)
      end if
      end if

      call datacast_dbl2(vy(1, k_1-0), vy(1, k_1-1),
     &                   vx(1, k_1-0), vx(1, k_1-1), u_t, v_t, n)

! FOR attention to unexpected overflow or NAN
      j_3 = eigen_loop_end     (L, x_nnod, x_inod)
      n   = eigen_translate_g2l(L, x_nnod, x_inod)
      if ( j_3 < n ) then
         vx(j_3+1:n, k_1-0) = ZERO ! in case
         vx(j_3+1:n, k_1-1) = ZERO ! in case
      end if

      j_3 = eigen_loop_end     (L, y_nnod, y_inod)
      n   = eigen_translate_g2l(L, y_nnod, y_inod)
      if ( j_3 < n ) then
         vy(j_3+1:n, k_1-0) = ZERO ! in case
         vy(j_3+1:n, k_1-1) = ZERO ! in case
      end if


      return
      end subroutine  eigen_prd_vo

