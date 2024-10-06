      subroutine  tred1_init(a, nm, n,
     &              d_out, e_out,
     &              u_t, v_t, nv)
!$    use omp_lib
      implicit NONE

      integer, intent(in)    ::  nm, n, nv
      real(8), intent(inout) ::  a(1:nm, *)
      real(8), intent(out)   ::  d_out(*)
      real(8), intent(out)   ::  e_out(*)
      real(8), intent(in)    ::  u_t(*), v_t(*)

      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  k_1, k_2, k_3
      integer                ::  i, j

      integer                ::  local_size, local_rank

      include 'trd.h'
      include 'trd_au.h'
      include 'CSTAB.h'
      include 'param.h'


      d_out(1:n) = ZERO
      e_out(1:n) = ZERO

      if ( diag_0 > 0 ) then
         i_2 = diag_0
         i_3 = loop_c_end  (n, y_nnod, y_inod)
         if ( i_2 <= i_3 ) then
         i_4 = x_nnod/n_common
         j_2 = diag_1
         j_3 = loop_c_end  (n, x_nnod, x_inod)
         j_4 = y_nnod/n_common
         k_2 = 0
         k_3 = (i_3-i_2)/i_4
         do  k_1 = k_2, k_3
            i_1 = i_2 + k_1 * i_4
            j_1 = j_2 + k_1 * j_4
            j   = (i_1-1)*y_nnod+y_inod
            d_out(j) = a(j_1, i_1)
         end do ! k_1
         end if
      end if

      i_2 = loop_c_start(1, y_nnod, y_inod)
      i_3 = loop_c_end  (n, y_nnod, y_inod)
      do  i_1 = i_2, i_3
         j   = loop_c_l2g_depth(i_1, y_nnod, y_inod)
         j_2 = loop_c_start    (j+1, x_nnod, x_inod)
         if ( j <= n ) then
            a(j_2:nm, i_1) = ZERO
         else
            a(1:nm, i_1) = ZERO
         end if
      end do ! i_1

      local_rank = 0
      local_size = 1
!$    local_rank = omp_get_thread_num()
!$    local_size = omp_get_num_threads()

      if ( local_size > 1 ) then

         allocate(u0_z(nv*local_size+n_columns),
     &            v0_z(nv*local_size+n_columns))
         u0_z = 0.0D0
         v0_z = 0.0D0
         call CSTAB_adjust_base(u0_z(1), u_t(1), offset1)
         call CSTAB_adjust_base(v0_z(1), v_t(1), offset2)

      end if


      return
      end subroutine  tred1_init


      subroutine  tred1_final(a, nm, n, d_out, e_out, u_t)

!$    use omp_lib
      implicit NONE

      integer, intent(in)    ::  nm, n
      real(8), intent(inout) ::  a(1:nm, *)
      real(8), intent(out)   ::  d_out(*)
      real(8), intent(out)   ::  e_out(*)
      real(8), intent(out)   ::  u_t(*)

      integer                ::  x_pos, x_owner_nod
      integer                ::  y_pos, y_owner_nod
      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  k_1, k_2, k_3, k_4
      integer                ::  i, j, k, L
      real(8)                ::  t

      integer                ::  local_size, local_rank

      include 'trd.h'
      include 'trd_au.h'
      include 'CSTAB.h'


      if ( n >= 2 ) then
         i = 2; L = i-1
         y_owner_nod = loop_c_node (i, y_nnod, y_inod)
         x_owner_nod = loop_c_node (L, x_nnod, x_inod)
         i_1 = loop_c_g2l_depth(i, y_nnod, y_owner_nod)
         j_1 = loop_c_g2l_depth(L, x_nnod, x_owner_nod)
         if ( y_owner_nod == y_inod
     &         .AND. x_owner_nod == x_inod ) then
            e_out(i)    = -a(j_1, i_1)
            a(j_1, i_1) =2*a(j_1, i_1)
         end if
         j = x_nnod * (y_owner_nod-1) + (x_owner_nod-1) + 1
         call bcast_dbl( e_out(2), 1, j, TRD_COMM_WORLD )
      end if

      if ( diag_0 > 0 ) then
         i_2 = diag_0
         i_3 = loop_c_end  (n, y_nnod, y_inod)
         if ( i_2 <= i_3 ) then
         i_4 = x_nnod/n_common
         j_2 = diag_1
         j_3 = loop_c_end  (n, x_nnod, x_inod)
         j_4 = y_nnod/n_common
         k_2 = 0
         k_3 = (i_3-i_2)/i_4
         do  k_1 = k_2, k_3
            i_1 = i_2 + k_1 * i_4
            j_1 = j_2 + k_1 * j_4
            j   = (i_1-1)*y_nnod+y_inod
            t           = d_out(j)
            d_out(j)    = a(j_1, i_1)
            a(j_1, i_1) = t
         end do ! k_1
         end if
      end if

      call reduce_dbl(d_out(1),   u_t(1), n, 1, TRD_COMM_WORLD)

      local_rank = 0
      local_size = 1
!$    local_rank = omp_get_thread_num()
!$    local_size = omp_get_num_threads()

      if ( local_size > 1 ) then

         deallocate(u0_z,v0_z)

      end if

      return
      end subroutine tred1_final

