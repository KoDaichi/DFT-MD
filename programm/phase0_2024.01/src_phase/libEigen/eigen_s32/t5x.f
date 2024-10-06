       subroutine tred1_local_2update0(
     &               w, nm,
     &               ux, uy, vx, vy, nv, uxx, anorm2,
     &               i_base, i)
!
       implicit none
!
       integer, intent(in)    ::  nm, nv
       real(8), intent(inout) ::  w(1:nm)
       real(8), intent(in)    ::  ux(1:nv),uy(1:nv)
       real(8), intent(in)    ::  vx(1:nv),vy(1:nv)
       real(8), intent(out)   ::  uxx(1:nv),anorm2
       integer, intent(in)    ::  i_base, i
!
       integer                ::  k_1
       integer                ::  j, k, l
!
       integer                ::  i_1, i_2, i_3, i_4
       integer                ::  j_1, j_2, j_3, j_4
       integer                ::  l_1, l_2, l_3, l_4
       integer                ::  jj_1, jj_2, jj_3, jj_4
       integer                ::  LX
!
       include 'trd.h'
       include 'param.h'
       include 'CSTAB.h'
!
       real(8)                :: u_x, v_x
       real(8)                :: uy0, vy0
       real(8)                :: uy1, vy1
       real(8)                :: uy2, vy2
       real(8)                :: uy3, vy3
       real(8)                :: w0
       real(8)                :: w1
       real(8)                :: w2
       real(8)                :: w3
       real(8)                :: s, temp
       integer                :: ierr


          k_1 = i - i_base
          if ( k_1 <= 1 ) return
!
          L = i - 1
          i_2 = loop_c_start(L, y_nnod,y_inod)
          i_3 = loop_c_end  (L, y_nnod,y_inod)
          if ( i_2 > i_3 ) return

          j_2 = loop_c_start(1, x_nnod,x_inod)
          j_3 = loop_c_end  (L-1, x_nnod,x_inod)

          i_1=i_2

          uy0 = uy(i_1+0)
          vy0 = vy(i_1+0)

          anorm2 = ZERO
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do j_1=j_2,j_3
             u_x = ux(j_1)
             v_x = vx(j_1)
             w0 = w(j_1)
     &              -(u_x*vy0)
     &              -(v_x*uy0)
!             anorm2 = anorm2 + w0**2
             w(j_1) = w0
             uxx(j_1) = w0
          end do! j_1

          j_2 = loop_c_start(L, x_nnod,x_inod)
          j_3 = loop_c_end  (L, x_nnod,x_inod)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do j_1=j_2,j_3
             u_x = ux(j_1)
             v_x = vx(j_1)
             w0 = w(j_1)
     &              -(u_x*vy0)
     &              -(v_x*uy0)
             w(j_1) = w0
          end do! j_1

       return
       end subroutine ! tred1_local_2update

