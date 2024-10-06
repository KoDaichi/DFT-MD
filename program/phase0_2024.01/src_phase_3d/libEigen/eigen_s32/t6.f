       subroutine tred1_vortho(
     &               u_x, v_x,
     &               ux, vx, nv,
     &               u_t, v_t, prod_uv,
     &               i, i_base, m0)
!
       implicit none
!
       integer, intent(in)    ::  nv, i, i_base, m0
       real(8), intent(in)    ::  u_x(*)
       real(8), intent(out)   ::  v_x(*)
       real(8), intent(in)    ::  ux(1:nv,*)
       real(8), intent(in)    ::  vx(1:nv,*)
       real(8), intent(inout) ::  u_t(*)
       real(8), intent(inout) ::  v_t(*)
       real(8), intent(out)   ::  prod_uv
!
       integer                ::  j, k, L, n, LL
       integer                ::  k_1, k_2, k_3
!
       integer                ::  i_1, i_2, i_3, i_4
       integer                ::  j_1, j_2, j_3, j_4
       integer                ::  l_1, l_2, l_3, l_4
       integer                ::  jj_1, jj_2, jj_3, jj_4
       integer                ::  LX
!
       integer, parameter     ::  VTOL = 2048
!
       include 'trd.h'
       include 'param.h'
       include 'CSTAB.h'
!
       real(8)                ::  w0
       real(8)                ::  u0, v0
       real(8)                ::  u1, v1
       real(8)                ::  u2, v2
       real(8)                ::  ux0, vx0
       real(8)                ::  ux1, vx1
       real(8)                ::  ux2, vx2


          k_1 = i - i_base
          k_2 = m0

          L = i-1
          n  = loop_c_g2l_depth(L, x_nnod,x_inod)

! FOR attention to unexpected overflow or NAN
          j_3 = loop_c_end(L, x_nnod,x_inod)
          if ( j_3 < n ) then
             v_x(j_3+1:n) = ZERO ! in case
          end if

          prod_uv=0.0D+00
          if ( k_2 <= k_1 ) then
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
             do j_1=1,n
                prod_uv = prod_uv+v_x(j_1)*u_x(j_1)
             end do! j_1
             call reduce_dbl(prod_uv, v_t, 1, 1, x_COMM_WORLD)
             return
          end if
!
! v=v-(UV+VU)u
!
          l_2 = k_2-k_1
          do j=1,l_2*2+1
             u_t(j) = ZERO
          end do

          l_4 = MOD(k_2-k_1, 3)+k_1+1
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

          k_3 = loop_c_end(L, x_nnod,x_inod)

          do jj_1=jj_2,jj_3,LX
             j_2 = jj_1; j_3 = MIN(jj_1+LX-1, jj_3)

             if(l_4-1==k_1+1)then
             l_1 = k_1+1                           ! 0

                j = l_1-k_1
                u0 = u_t(2*(j+0)-1+1)
                v0 = u_t(2*(j+0)-0+1)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
                do j_1=j_2,j_3
                   w0 = u_x(j_1)
                   u0 = u0+vx(j_1,l_1+0)*w0
                   v0 = v0+ux(j_1,l_1+0)*w0
                end do! j_1
                u_t(2*(j+0)-1+1) = u0
                u_t(2*(j+0)-0+1) = v0
             end if
             if(l_4-2==k_1+1)then
             l_1 = k_1+1                           ! 1

                j = l_1-k_1
                u0 = u_t(2*(j+0)-1+1)
                v0 = u_t(2*(j+0)-0+1)
                u1 = u_t(2*(j+1)-1+1)
                v1 = u_t(2*(j+1)-0+1)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
                do j_1=j_2,j_3
                   w0 = u_x(j_1)
                   u0 = u0+vx(j_1,l_1+0)*w0
                   v0 = v0+ux(j_1,l_1+0)*w0
                   u1 = u1+vx(j_1,l_1+1)*w0
                   v1 = v1+ux(j_1,l_1+1)*w0
                end do! j_1
                u_t(2*(j+0)-1+1) = u0
                u_t(2*(j+0)-0+1) = v0
                u_t(2*(j+1)-1+1) = u1
                u_t(2*(j+1)-0+1) = v1
             end if
             do l_1=l_4,k_2,3                  ! 2

                j = l_1-k_1
                u0 = u_t(2*(j+0)-1+1)
                v0 = u_t(2*(j+0)-0+1)
                u1 = u_t(2*(j+1)-1+1)
                v1 = u_t(2*(j+1)-0+1)
                u2 = u_t(2*(j+2)-1+1)
                v2 = u_t(2*(j+2)-0+1)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
                do j_1=j_2,j_3
                   w0 = u_x(j_1)
                   u0 = u0+vx(j_1,l_1+0)*w0
                   v0 = v0+ux(j_1,l_1+0)*w0
                   u1 = u1+vx(j_1,l_1+1)*w0
                   v1 = v1+ux(j_1,l_1+1)*w0
                   u2 = u2+vx(j_1,l_1+2)*w0
                   v2 = v2+ux(j_1,l_1+2)*w0
                end do! j_1
                u_t(2*(j+0)-1+1) = u0
                u_t(2*(j+0)-0+1) = v0
                u_t(2*(j+1)-1+1) = u1
                u_t(2*(j+1)-0+1) = v1
                u_t(2*(j+2)-1+1) = u2
                u_t(2*(j+2)-0+1) = v2
             end do! l_1

             u0 = u_t(1)
             do j_1=j_2,j_3
                u0 = u0+v_x(j_1)*u_x(j_1)
             end do! j_1
             u_t(1)=u0

          end do! jj_1

       if ( n > VTOL ) then
          call reduce_dbl(u_t, v_t, (k_2-k_1)*2+1, 1, y_COMM_WORLD)
       end if
          call reduce_dbl(u_t, v_t, (k_2-k_1)*2+1, 1, x_COMM_WORLD)

       prod_uv = u_t(1)
       do l_1=k_1+1,k_2
          j = l_1-k_1
          prod_uv =prod_uv - 2*u_t(2*(j+0)-1+1)*u_t(2*(j+0)-0+1)
       enddo

       if ( n > VTOL ) then
          jj_2 = 1+LL*(y_inod-1)
          jj_3 = MIN(n, LL*y_inod)
       else
          jj_2 = 1
          jj_3 = n
       endif
          do jj_1=jj_2,jj_3,LX
             j_2 = jj_1; j_3 = MIN(jj_1+LX-1, jj_3)

             if(l_4-1==k_1+1)then
             l_1 = k_1+1                           ! 0

                j = l_1-k_1

                u0 = u_t(2*(j+0)-1+1)
                v0 = u_t(2*(j+0)-0+1)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
                do j_1=j_2,j_3
                   w0 = v_x(j_1)
                   ux0 = ux(j_1,l_1+0)
                   vx0 = vx(j_1,l_1+0)
                   w0 = w0
     &                -ux0*u0
     &                -vx0*v0
                   v_x(j_1) = w0
                end do! j_1
             end if
             if(l_4-2==k_1+1)then
             l_1 = k_1+1                           ! 1

                j = l_1-k_1

                u0 = u_t(2*(j+0)-1+1)
                v0 = u_t(2*(j+0)-0+1)
                u1 = u_t(2*(j+1)-1+1)
                v1 = u_t(2*(j+1)-0+1)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
                do j_1=j_2,j_3
                   w0 = v_x(j_1)
                   ux0 = ux(j_1,l_1+0)
                   vx0 = vx(j_1,l_1+0)
                   w0 = w0
     &                -ux0*u0
     &                -vx0*v0
                   ux1 = ux(j_1,l_1+1)
                   vx1 = vx(j_1,l_1+1)
                   w0 = w0
     &                -ux1*u1
     &                -vx1*v1
                   v_x(j_1) = w0
                end do! j_1
             end if
             do l_1=l_4,k_2,3                  ! 2

                j = l_1-k_1

                u0 = u_t(2*(j+0)-1+1)
                v0 = u_t(2*(j+0)-0+1)
                u1 = u_t(2*(j+1)-1+1)
                v1 = u_t(2*(j+1)-0+1)
                u2 = u_t(2*(j+2)-1+1)
                v2 = u_t(2*(j+2)-0+1)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
                do j_1=j_2,j_3
                   w0 = v_x(j_1)
                   ux0 = ux(j_1,l_1+0)
                   vx0 = vx(j_1,l_1+0)
                   w0 = w0
     &                -ux0*u0
     &                -vx0*v0
                   ux1 = ux(j_1,l_1+1)
                   vx1 = vx(j_1,l_1+1)
                   w0 = w0
     &                -ux1*u1
     &                -vx1*v1
                   ux2 = ux(j_1,l_1+2)
                   vx2 = vx(j_1,l_1+2)
                   w0 = w0
     &                -ux2*u2
     &                -vx2*v2
                   v_x(j_1) = w0
                end do! j_1
             end do! l_1

          end do! jj_1

       if ( n > VTOL ) then
!          do j_1=1,jj_2-1
!             vx(j_1,k_1) = ZERO
!          enddo
!          do j_1=jj_3+1,n
!             vx(j_1,k_1) = ZERO
!          enddo
!          call reduce_dbl(vx(1,k_1), v_t, n, 1, y_COMM_WORLD)
          call allgather_dbl(v_x(jj_2), v_t, LL, y_COMM_WORLD)
          j_3 = loop_c_end(L, x_nnod,x_inod)
          v_x(1:j_3)=v_t(1:j_3)
       endif

! FOR attention to unexpected overflow or NAN
          j_3 = loop_c_end(L, x_nnod,x_inod)
          if ( j_3 < n ) then
             v_x(j_3+1:n) = ZERO ! in case
          end if

      return
      end subroutine ! tred1_vorhto

