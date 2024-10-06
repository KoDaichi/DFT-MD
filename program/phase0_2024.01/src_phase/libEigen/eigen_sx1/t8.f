      subroutine tred1_init(a, nm, n,
     &              d_out, e_out, ne,
     &              u_t, v_t, nv)
!$    use OMP_LIB
      implicit NONE

      integer, intent(in)    :: nm, n, ne, nv
      real(8), intent(inout) :: a(1:nm, *)
      real(8), intent(out)   :: d_out(*)
      real(8), intent(out)   :: e_out(1:ne, *)
      real(8), intent(in)    :: u_t(*), v_t(*)

      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  i, j

      integer                :: kx
      integer                :: local_size, local_rank

      real(8), pointer       :: u0_z(:),v0_z(:)
      real(8), pointer       :: u1_z(:),v1_z(:)
      integer                :: offset1, offset2
      integer                :: offset3, offset4
      common/tred_au_common/    u0_z, v0_z, u1_z, v1_z,
     &                          offset1, offset2, offset3, offset4

      include 'trd.h'
      include 'param.h'
      include 'CSTAB.h'


      d_out(1:n)    = ZERO
      e_out(1:n,1)  = ZERO
      e_out(1:n,2)  = ZERO

      if ( diag_0 > 0 ) then
         i_2 = diag_0
         i_3 = loop_c_end  (n, y_nnod,y_inod)
         i_4 = x_nnod/n_common
         do i_1=i_2,i_3,i_4
!           j   = loop_c_l2g_depth(i_1, y_nnod,y_inod)
!           j_1 = loop_c_g2l_depth(j, x_nnod,x_inod)
            j   = (i_1-1)*y_nnod+y_inod
            j_1 = (j-1)/x_nnod+1
            d_out(j) = a(j_1,i_1)
         end do! i_1
      end if

      i_2 = loop_c_start(1, y_nnod,y_inod)
      i_3 = loop_c_end  (n, y_nnod,y_inod)
      do i_1=i_2,i_3
         j   = loop_c_l2g_depth(i_1, y_nnod,y_inod)
         j_2 = loop_c_start(j+1, x_nnod,x_inod)
         do j_1=j_2,nm
            a(j_1,i_1) = ZERO
         end do! j_1
         if ( j > n ) then
         do j_1=j_2,nm
            a(1:nm,i_1) = ZERO
         end do! j_1
         end if
      end do! i_1

      local_rank = 0
      local_size = 1
!$    local_rank = omp_get_thread_num()
!$    local_size = omp_get_num_threads()

      allocate(u0_z(nv*local_size+n_columns),
     &         v0_z(nv*local_size+n_columns))
      allocate(u1_z(nv*local_size+n_columns),
     &         v1_z(nv*local_size+n_columns))
      call CSTAB_adjust_base(u0_z(1), u_t(1), offset1)
      call CSTAB_adjust_base(v0_z(1), v_t(1), offset2)
      call CSTAB_adjust_base(u1_z(1), u_t(1), offset3)
      call CSTAB_adjust_base(v1_z(1), v_t(1), offset4)
      kx =  (L1_WINDOW/8)
!     &           +(L1_WINDOW)
!     &           +(L1_LSIZE/8)
     &           +(L1_LSIZE)
     &           +(L2_LSIZE/8)
      offset1 = offset1 + kx * 1
      offset2 = offset2 + kx * 2
      offset3 = offset3 + kx * 3
      offset4 = offset4 + kx * 4
      call CSTAB_round_offset(offset1)
      call CSTAB_round_offset(offset2)
      call CSTAB_round_offset(offset3)
      call CSTAB_round_offset(offset4)

      return
      end subroutine tred1_init


      subroutine tred1_final(a, nm, n, d_out, e_out, ne, u_t)

      implicit NONE

      integer, intent(in)    :: nm, ne, n
      real(8), intent(inout) :: a(1:nm, *)
      real(8), intent(out)   :: d_out(*)
      real(8), intent(out)   :: e_out(1:ne, *)
      real(8), intent(out)   :: u_t(*)

      integer                ::  x_pos, x_owner_nod
      integer                ::  y_pos, y_owner_nod
      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  i, j, k, L
      real(8)                ::  t

      real(8), pointer       :: u0_z(:),v0_z(:)
      real(8), pointer       :: u1_z(:),v1_z(:)
      integer                :: offset1, offset2
      integer                :: offset3, offset4
      common/tred_au_common/    u0_z, v0_z, u1_z, v1_z,
     &                          offset1, offset2, offset3, offset4

      include 'trd.h'
      include 'param.h'


      if ( n >= 2 ) then
         do i = 2+MOD(n,2), 2, -1
            L = i-1
            y_owner_nod = loop_c_node (i, y_nnod,y_inod)
            x_owner_nod = loop_c_node (L, x_nnod,x_inod)
            if ( y_owner_nod == y_inod
     &         .AND. x_owner_nod == x_inod ) then
               i_1 = loop_c_g2l_depth(i, y_nnod,y_inod)
               j_1 = loop_c_g2l_depth(L, x_nnod,x_inod)
               e_out(i,1)  =  a(j_1,i_1)
               a(j_1,i_1)  =  ZERO
            end if
            L = i-2; if ( L < 1 ) cycle
            y_owner_nod = loop_c_node (i, y_nnod,y_inod)
            x_owner_nod = loop_c_node (L, x_nnod,x_inod)
            if ( y_owner_nod == y_inod
     &         .AND. x_owner_nod == x_inod ) then
               i_1 = loop_c_g2l_depth(i, y_nnod,y_inod)
               j_1 = loop_c_g2l_depth(L, x_nnod,x_inod)
               e_out(i,2)  =  a(j_1,i_1)
               a(j_1,i_1)  =  ZERO
            end if
         end do
      end if

      if ( diag_0 > 0 ) then
         i_2 = diag_0
         i_3 = loop_c_end  (n, y_nnod,y_inod)
         i_4 = x_nnod/n_common
         do i_1=i_2,i_3,i_4
            j   = loop_c_l2g_depth(i_1, y_nnod,y_inod)
            j_1 = loop_c_g2l_depth(j, x_nnod,x_inod)
            t          = d_out(j)
            d_out(j)   = a(j_1,i_1)
            a(j_1,i_1) = t
         end do! i_1
      end if

      call reduce_dbl(d_out(1),   u_t(1), n, 1, TRD_COMM_WORLD)
      call reduce_dbl(e_out(1,1), u_t(1), n, 1, TRD_COMM_WORLD)
      call reduce_dbl(e_out(1,2), u_t(1), n, 1, TRD_COMM_WORLD)

      deallocate(u0_z, v0_z)
      deallocate(u1_z, v1_z)


      return
      end subroutine tred1_final

