       subroutine tred1_v(u_x, v_x, v_y, nv, u_t, v_t, beta, prod_uv, i)

       implicit NONE

       integer, intent(in)    ::  nv, i
       real(8), intent(in)    ::  u_x(1:nv)
       real(8), intent(inout) ::  v_x(1:nv), v_y(1:nv)
       real(8), intent(out)   ::  u_t(*), v_t(*)
       real(8), intent(in)    ::  beta, prod_uv

       real(8)                ::  alpha, temp
!       save                       alpha

       integer                ::  x_pos, x_owner_nod
       integer                ::  y_pos, y_owner_nod
       integer                ::  j_1, j_2, j_3
       integer                ::  L
       integer                ::  n, LL, jj_1, jj_2, jj_3

!       integer, parameter     ::  VTOL = 512
       integer, parameter     ::  VTOL = 51200

       include 'trd.h'
       include 'param.h'


       L = i-1

       x_owner_nod = loop_c_node (L, x_nnod,x_inod)
       y_owner_nod = loop_c_node (i, y_nnod,y_inod)

       x_pos       = loop_c_g2l_depth(L, x_nnod,x_inod)
       y_pos       = loop_c_g2l_depth(i, y_nnod,y_inod)

       j_2         = loop_c_start(1, x_nnod,x_inod)
       j_3         = loop_c_end  (L, x_nnod,x_inod)

       n = j_3-j_2+1
       LL = (n-1)/y_nnod+1
       LL = ((LL-1)/2+1)*2

       if ( n > VTOL ) then
          jj_2 = j_2+(1+LL*(y_inod-1))-1
          jj_3 = j_2+(MIN(n,LL*y_inod))-1
       else
          jj_2 = j_2
          jj_3 = j_3
       endif
!
! v':= v-((u,v)/2|u|^2)u
!
       if ( beta /= ZERO ) then

          alpha = prod_uv/(2*beta)

!DIR$ IVDEP
          do j_1=jj_2,jj_3
             v_x(j_1) = (v_x(j_1)-alpha*u_x(j_1))/beta
          end do! j_1

          if ( n > VTOL ) then
             call allgather_dbl(v_x(jj_2), v_t(j_2), LL, y_COMM_WORLD)
             v_x(j_2:j_3) = v_t(j_2:j_3)
          end if

          call datacast_dbl(v_y(1), v_x(1), u_t(1), v_t(1), x_pos)

       else

          x_pos       = loop_c_g2l_depth(i, x_nnod,x_inod)
          y_pos       = loop_c_g2l_depth(i, y_nnod,y_inod)

          v_x(1:x_pos) = ZERO
          v_y(1:y_pos) = ZERO

       end if

! FOR attention to unexpected overflow or NAN
       j_3 = loop_c_end      (L, x_nnod, x_inod)
       n   = loop_c_g2l_depth(L, x_nnod, x_inod)
       if ( j_3 < n ) then
          v_x(j_3+1:n) = ZERO ! in case
       end if

       j_3 = loop_c_end      (L, y_nnod, y_inod)
       n   = loop_c_g2l_depth(L, y_nnod, y_inod)
       if ( j_3 < n ) then
          v_y(j_3+1:n) = ZERO ! in case
       end if


       return
       end subroutine tred1_v

