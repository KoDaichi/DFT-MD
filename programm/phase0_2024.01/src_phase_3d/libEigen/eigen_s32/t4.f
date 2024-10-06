       subroutine  tred1_u(
     &               a, nm,
     &               u_x, u_y, nv,
     &               u_t, v_t, i, beta)
       implicit NONE

       integer, intent(in)    ::  nm, nv, i
       real(8), intent(inout) ::  a(1:nm)
       real(8), intent(inout) ::  u_x(1:nv), u_y(1:nv)
       real(8), intent(inout) ::  u_t(*), v_t(*)
       real(8), intent(inout) ::  beta

       include 'trd.h'

       real(8)                ::  anorm2, a_n, g_n
       real(8)                ::  tt(4), ss(4), t, s

       integer                ::  x_owner_nod, x_pos
       integer                ::  y_owner_nod, y_pos
       integer                ::  j_1, j_2, j_3
       integer                ::  jj_1, jj_2, jj_3
       integer                ::  L

       include 'param.h'


       L = i-1

       x_owner_nod = loop_c_node (L, x_nnod, x_inod)
       y_owner_nod = loop_c_node (i, y_nnod, y_inod)

       x_pos       = loop_c_g2l_depth(L, x_nnod, x_inod)
       y_pos       = loop_c_g2l_depth(i, y_nnod, y_inod)

       j_2         = loop_c_start(1, x_nnod, x_inod)
       j_3         = loop_c_end  (L, x_nnod, x_inod)

!
! u=...
!
       if ( y_owner_nod == y_inod ) then
       if ( beta < ZERO ) then
          u_x(j_2:j_3) = a(j_2:j_3)
       end if
       end if

       call bcast_dbl(u_x(1), x_pos, y_owner_nod, y_COMM_WORLD)

       if ( x_owner_nod == x_inod ) then
          u_x(j_3+1:x_pos) = ZERO
       end if

       call datacast_dbl(u_y(1), u_x(1), u_t(1), v_t(1), x_pos)


       return
       end subroutine  tred1_u

