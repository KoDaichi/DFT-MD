      subroutine tred1_vorhto(
     &              ux, vx, nv,
     &              u_t, v_t,
     &              i, i_base, m0)

!
      implicit none
!
      integer, intent(in)    ::  nv, i, i_base, m0
      real(8), intent(inout) ::  ux(1:nv,*)
      real(8), intent(inout) ::  vx(1:nv,*)
      real(8), intent(inout) ::  u_t(4,*)
      real(8), intent(inout) ::  v_t(*)
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
      include 'trd.h'
      include 'param.h'
      include 'CSTAB.h'
!
      real(8)                ::  w0, w1
      real(8)                ::  u0_0, v0_0
      real(8)                ::  u1_0, v1_0
      real(8)                ::  u0_1, v0_1
      real(8)                ::  u1_1, v1_1
      real(8)                ::  ux0, vx0
      real(8)                ::  ux1, vx1


      k_1 = i - i_base
      k_2 = m0

      if ( k_2 <= k_1 ) return

      L = i - MBAND
      n  = loop_c_g2l_depth(L, x_nnod,x_inod)

! FOR attention to unexpected overflow or NAN
          j_3 = loop_c_end(L, x_nnod,x_inod)
          if ( j_3 < n ) then
             vx(j_3+1:n,k_1-0) = ZERO ! in case
             vx(j_3+1:n,k_1-1) = ZERO ! in case
          end if
!
! v=v-(UV+VU)u
!
      l_2 = k_2-k_1
      do j=1,l_2*(2*MBAND)
         u_t(j,1) = ZERO
      end do

      l_4 = MOD(k_2-k_1,2)+k_1+1
      LX  = L1_LSIZE*L1_WAY/16

      LL = (n-1)/y_nnod+1
      LL = ((LL-1)/2+1)*2

      jj_2 = 1+LL*(y_inod-1)
      jj_3 = MIN(n,LL*y_inod)
      do jj_1=jj_2,jj_3,LX
         j_2 = jj_1; j_3 = MIN(jj_1+LX-1,jj_3)

         do l_1=k_1+1,l_4-1                        ! 0

            j = l_1-k_1
            u0_0 = u_t(1,j+0)
            v0_0 = u_t(2,j+0)
            u1_0 = u_t(3,j+0)
            v1_0 = u_t(4,j+0)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
            do j_1=j_2,j_3
               w0 = ux(j_1,k_1-0)
               w1 = ux(j_1,k_1-1)
               u0_0 = u0_0 + vx(j_1,l_1+0) * w0
               v0_0 = v0_0 + ux(j_1,l_1+0) * w0
               u1_0 = u1_0 + vx(j_1,l_1+0) * w1
               v1_0 = v1_0 + ux(j_1,l_1+0) * w1
            end do! j_1
            u_t(1,j+0) = u0_0
            u_t(2,j+0) = v0_0
            u_t(3,j+0) = u1_0
            u_t(4,j+0) = v1_0
         end do! l_1
         do l_1=l_4,k_2,2                      ! 1

            j = l_1-k_1
            u0_0 = u_t(1,j+0)
            v0_0 = u_t(2,j+0)
            u1_0 = u_t(3,j+0)
            v1_0 = u_t(4,j+0)
            u0_1 = u_t(1,j+1)
            v0_1 = u_t(2,j+1)
            u1_1 = u_t(3,j+1)
            v1_1 = u_t(4,j+1)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
            do j_1=j_2,j_3
               w0 = ux(j_1,k_1-0)
               w1 = ux(j_1,k_1-1)
               u0_0 = u0_0 + vx(j_1,l_1+0) * w0
               v0_0 = v0_0 + ux(j_1,l_1+0) * w0
               u1_0 = u1_0 + vx(j_1,l_1+0) * w1
               v1_0 = v1_0 + ux(j_1,l_1+0) * w1
               u0_1 = u0_1 + vx(j_1,l_1+1) * w0
               v0_1 = v0_1 + ux(j_1,l_1+1) * w0
               u1_1 = u1_1 + vx(j_1,l_1+1) * w1
               v1_1 = v1_1 + ux(j_1,l_1+1) * w1
            end do! j_1
            u_t(1,j+0) = u0_0
            u_t(2,j+0) = v0_0
            u_t(3,j+0) = u1_0
            u_t(4,j+0) = v1_0
            u_t(1,j+1) = u0_1
            u_t(2,j+1) = v0_1
            u_t(3,j+1) = u1_1
            u_t(4,j+1) = v1_1
         end do! l_1

      end do! jj_1

      call reduce_dbl(u_t, v_t, (k_2-k_1)*(2*MBAND), 1, y_COMM_WORLD)
      call reduce_dbl(u_t, v_t, (k_2-k_1)*(2*MBAND), 1, x_COMM_WORLD)

      if ( n > 100000 ) then
         jj_2 = 1+LL*(y_inod-1)
         jj_3 = MIN(n,LL*y_inod)
      else
         jj_2 = 1
         jj_3 = n
      end if
      do jj_1=jj_2,jj_3,LX
         j_2 = jj_1; j_3 = MIN(jj_1+LX-1,jj_3)

         do l_1=k_1+1,l_4-1                        ! 0

            j = l_1-k_1

            u0_0 = u_t(1,j+0)
            v0_0 = u_t(2,j+0)
            u1_0 = u_t(3,j+0)
            v1_0 = u_t(4,j+0)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
            do j_1=j_2,j_3
               w0 = vx(j_1,k_1-0)
               w1 = vx(j_1,k_1-1)
               ux0 = ux(j_1,l_1+0)
               vx0 = vx(j_1,l_1+0)
               w0 = w0
     &            - ux0 * u0_0
     &            - vx0 * v0_0
               w1 = w1
     &            - ux0 * u1_0
     &            - vx0 * v1_0
               vx(j_1,k_1-0) = w0
               vx(j_1,k_1-1) = w1
            end do! j_1
         end do! l_1
         do l_1=l_4,k_2,2                      ! 1

            j = l_1-k_1

            u0_0 = u_t(1,j+0)
            v0_0 = u_t(2,j+0)
            u1_0 = u_t(3,j+0)
            v1_0 = u_t(4,j+0)
            u0_1 = u_t(1,j+1)
            v0_1 = u_t(2,j+1)
            u1_1 = u_t(3,j+1)
            v1_1 = u_t(4,j+1)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
            do j_1=j_2,j_3
               w0 = vx(j_1,k_1-0)
               w1 = vx(j_1,k_1-1)
               ux0 = ux(j_1,l_1+0)
               vx0 = vx(j_1,l_1+0)
               w0 = w0
     &            - ux0 * u0_0
     &            - vx0 * v0_0
               w1 = w1
     &            - ux0 * u1_0
     &            - vx0 * v1_0
               ux1 = ux(j_1,l_1+1)
               vx1 = vx(j_1,l_1+1)
               w0 = w0
     &            - ux1 * u0_1
     &            - vx1 * v0_1
               w1 = w1
     &            - ux1 * u1_1
     &            - vx1 * v1_1
               vx(j_1,k_1-0) = w0
               vx(j_1,k_1-1) = w1
            end do! j_1
         end do! l_1

      end do! jj_1

      if ( n > 100000 ) then
         do j_1=1,jj_2-1
            vx(j_1,k_1-0) = ZERO
            vx(j_1,k_1-1) = ZERO
         end do! j_1
         do j_1=jj_3+1,n
            vx(j_1,k_1-0) = ZERO
            vx(j_1,k_1-1) = ZERO
         end do! j_1

         do j_1=1,n
            v_t(j_1)   = vx(1,k_1-0)
            v_t(j_1+n) = vx(1,k_1-1)
         end do! j_1
         call reduce_dbl(v_t, u_t, MBAND*n, 1, y_COMM_WORLD)
         do j_1=1,n
            vx(1,k_1-0) = v_t(j_1)
            vx(1,k_1-1) = v_t(j_1+n)
         end do! j_1
      end if

!  FOR attention to unexpected overflow or NAN
          j_3 = loop_c_end(L, x_nnod,x_inod)
          if ( j_3 < n ) then
             vx(j_3+1:n,k_1-0) = ZERO ! in case
             vx(j_3+1:n,k_1-1) = ZERO ! in case
          end if


      return
      end subroutine ! tred1_vorhto

