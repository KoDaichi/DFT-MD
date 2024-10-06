       subroutine  tred1_au(
     &                a, nm,
     &                u_x, u_y, v_x, nv,
     &                u_t, v_t, d_t,
     &                i, i_base, m,
     &                w, e, beta)
!$     use omp_lib
       implicit NONE

       integer, intent(in)    ::  nm, nv, i, i_base, m
       real(8), intent(inout) ::  a(1:nm,*)
       real(8), intent(inout) ::  u_x(1:nv), u_y(1:nv)
       real(8), intent(inout) ::  v_x(1:nv)
       real(8), intent(inout) ::  u_t(*), v_t(*)
       real(8), intent(in)    ::  d_t(*)
       real(8), intent(inout) ::  w(*), e(*), beta

       integer                ::  n, blk_0
       integer                ::  x_pos, y_pos, x_root, y_root
       integer                ::  n1, n2
       integer                ::  i_1, i_2, i_3
       integer                ::  j_1, j_2, j_3
       integer                ::  k_1, k_2, k_3
       integer                ::  L

       real(8) :: anorm2, a_n, g_n

       real(8) :: d1, d2

       integer                ::  local_size, local_rank

       include 'trd_au.h'

       real(8)      t2_reduce1, t2_reduce2
       common /t2/  t2_reduce1, t2_reduce2


       include 'mpif.h'
       include 'trd.h'
       include 'CSTAB.h'
       include 'param.h'


          local_rank = 0
          local_size = 1
!$        local_rank = omp_get_thread_num()
!$        local_size = omp_get_num_threads()

             L = i-1

             i_2 = loop_c_start(1, y_nnod,y_inod)
             i_3 = loop_c_end  (L, y_nnod,y_inod)
             j_2 = loop_c_start(1, x_nnod,x_inod)
             j_3 = loop_c_end  (L, x_nnod,x_inod)

             y_root = loop_c_node (L, y_nnod,y_inod)
             y_pos  = loop_c_g2l_depth(i, y_nnod,y_inod)

             x_root = loop_c_node (L, x_nnod,x_inod)
             x_pos  = loop_c_g2l_depth(L, x_nnod,x_inod)

             k_1   = i-i_base
             k_2   = m
             n     = i_3
             blk_0 = k_1-k_2


       if ( local_size > 1 ) then

          n1 = offset1+nv*local_rank
          n2 = offset2+nv*local_rank

          do j_1=1,x_pos
             u0_z(j_1+n1) = ZERO
          enddo
          do j_1=1,y_pos
             v0_z(j_1+n2) = ZERO
          enddo

          call  tred1_au_body1(
     &             a, nm,
     &             u_x, u_y, u0_z(1+n1), v0_z(1+n2),
     &             1, n, x_pos, y_pos, nv, blk_0
     &             ,local_rank, local_size
     &             )

!$OMP BARRIER

          call  tred1_au_body2(
     &             v_x, u0_z(1+offset1),
     &             v_t, v0_z(1+offset2),
     &             nv, x_pos, y_pos
     &             ,local_rank, local_size
     &             )

!$OMP BARRIER

       else

          do j_1=1,x_pos
             v_x(j_1) = ZERO
          enddo
          do j_1=1,y_pos
             v_t(j_1) = ZERO
          enddo

          call  tred1_au_body1(
     &             a, nm,
     &             u_x, u_y, v_x, v_t,
     &             1, n, x_pos, y_pos, nv, blk_0
     &             ,local_rank, local_size
     &             )

       endif

!          if ( local_rank == 0 ) then
!$OMP MASTER
             anorm2 = ZERO
             do  j_1=1,x_pos
                anorm2   =  anorm2 + u_x(j_1)**2
             end do ! j_1

             if ( x_root == x_inod ) then
                a_n = u_x(x_pos)
             else
                a_n = ZERO
             endif

             v_t(y_pos+1) = anorm2
             v_t(y_pos+2) = a_n

             if ( TRD_nnod > 1 ) then
            d1=MPI_Wtime()
                call reduce_dbl(v_t, u_t, y_pos+2, 1, x_COMM_WORLD)
            d2=MPI_Wtime()
            t2_reduce1=t2_reduce1+(d2-d1)
             end if

             anorm2 =  v_t(y_pos+1)
             a_n    =  v_t(y_pos+2)
             g_n    = -SIGN(SQRT(anorm2), a_n)
             beta   =  anorm2 - a_n * g_n
             e (i)  =  g_n

             if ( x_root == x_inod ) then
                u_x(x_pos) =  a_n - g_n
                i_1 = loop_c_node (i, y_nnod,y_inod)
                if ( i_1 == y_inod ) then
                   w(x_pos) =  u_x(x_pos)
                end if
             end if

             if ( y_root == y_inod ) then
                i_1  = loop_c_g2l_depth(L, y_nnod,y_inod)
                u_y(i_1) =  a_n - g_n
                do  j_1=1,x_pos
                   v_x(j_1) = v_x(j_1) - g_n * a(j_1,i_1)
                end do
             end if

             call  tred1_au_body3(
     &                u_x, u_y, v_x,
     &                u_t, v_t, d_t,
     &                1, n, x_pos, y_pos, nv
     &                )

             if ( TRD_nnod > 1 ) then
            d1=MPI_Wtime()
                call reduce_dbl(v_x, v_t, x_pos, x_nnod, y_COMM_WORLD)
            d2=MPI_Wtime()
            t2_reduce2=t2_reduce2+(d2-d1)
             end if

!$OMP END MASTER
!          end if


       return
       end subroutine ! tred1_au

       subroutine  tred1_au_body1(
     &                a, nm,
     &                u_x, u_y, u_t, v_t,
     &                n1, n2, x_pos, y_pos, nv, blk_0
     &                ,local_rank, local_size
     &                )
!$     use omp_lib
       implicit NONE

       integer, intent(in)    :: nm, nv, n1, n2
       real(8), intent(inout) :: a(1:nm,*)
       real(8), intent(inout) :: u_x(1:nv), u_y(1:nv)
       real(8), intent(inout) :: u_t(1:nv), v_t(1:nv)
       integer, intent(in)    :: blk_0
       integer                :: x_pos, y_pos, x_root, y_root
       integer                :: local_size, local_rank

       integer                :: i_0
       integer                :: i_1, i_2, i_3, i_4
       integer                :: j_1, j_2, j_3, j_4
       integer                :: k_1, k_2, k_3, k_4
       integer                :: l_1, l_2, l_3, l_4
       integer                :: i, j, k

       real(8)                :: v0, u0
       real(8)                :: a0_0
       real(8)                :: a0_1
       real(8)                :: a0_2
       real(8)                :: a0_3
       real(8)                :: a0_4
       real(8)                :: a0_5
       real(8)                :: a0_6
       real(8)                :: a0_7
       real(8)                :: v_t0, u_y0
       real(8)                :: v_t1, u_y1
       real(8)                :: v_t2, u_y2
       real(8)                :: v_t3, u_y3
       real(8)                :: v_t4, u_y4
       real(8)                :: v_t5, u_y5
       real(8)                :: v_t6, u_y6
       real(8)                :: v_t7, u_y7

       integer                ::   LX, LY
       integer                ::   ii_1, ii_2, ii_3, ii_4, ii_5
       integer                ::   jj_1, jj_2, jj_3, jj_4, jj_5
       integer                ::   kk_1, kk_2, kk_3, kk_4, kk_5

       include 'trd.h'
       include 'param.h'


          i_2 = n1
          i_3 = n2
!
! v:= Au
!
          if ( blk_0 == 0 ) then
             do i_1=i_2+local_rank,i_3,local_size
                j   = loop_c_l2g_depth(i_1, y_nnod,y_inod)
                j_3 = loop_c_end  (j, x_nnod,x_inod)
                j   = j+y_nnod*(6*2-1)
                j_4 = loop_c_end  (j, x_nnod,x_inod)
                do j_1=j_3+1,MIN(j_4,nm)
                   a(j_1,i_1) = ZERO
                end do! j_1
             end do! i_1
!$OMP BARRIER
          end if
!
!
          LX = 32*1000; LY = 32*1
!
          l_2 = i_2; l_3 = i_3
          l_1 = l_2; l_4 = l_3

             k_2 = 1
             k   = loop_c_l2g_depth(l_4, y_nnod,y_inod)
             k_3 = loop_c_end  (k-1, x_nnod,x_inod)
             do k_1=k_2,k_3,LX; k_4 = MIN(k_3,k_1+LX-1)

                j    = loop_c_l2g_depth(k_1, x_nnod,x_inod)
                ii_2 = loop_c_start(j,y_nnod,y_inod)
                ii_2 = MAX(l_1,ii_2)
                ii_3 = l_4
                if ( ii_2 > ii_3 ) cycle
                ii_4 = MOD(ii_3-ii_2+1,6*2)+ii_2


                do i_1=ii_2+local_rank,ii_4-1,local_size

                   j    = loop_c_l2g_depth(i_1, y_nnod,y_inod)
                   j    = j+(1-1)*y_nnod
                   jj_2 = k_1
                   jj_3 = loop_c_end  (j-1, x_nnod,x_inod)
                   jj_3 = MIN(k_4,jj_3)
                   if ( jj_2 > jj_3 ) cycle

             do kk_1=jj_2,jj_3,64*64; kk_4=MIN(kk_1+64*64-1,jj_3)

                   u_y0 = u_y(i_1+0)

!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
                   do j_1=kk_1,kk_4
!                   do j_1=jj_2,jj_3

                      u0 = u_x(j_1+0)
                      a0_0 = a(j_1+0,i_1+0)
                      v0 = u_t(j_1+0)


                      v0 = v0
     &                            + (a0_0*u_y0)
                      v_t(i_1+0) = v_t(i_1+0)
     &                            + (a0_0*u0)

                      u_t(j_1+0) = v0

                   end do! j_1

            end do! kk_1

                end do! i_1

                do i_0=ii_4+local_rank*6*2,ii_3, local_size*6*2

                   j    = loop_c_l2g_depth(i_0, y_nnod,y_inod)
                   j    = j+(6*2-1)*y_nnod
                   jj_2 = k_1
                   jj_3 = loop_c_end  (j-1, x_nnod,x_inod)
                   jj_3 = MIN(k_4,jj_3)
                   if ( jj_2 > jj_3 ) cycle

             do kk_1=jj_2,jj_3,16*31; kk_4=MIN(kk_1+16*31-1,jj_3)
             do i_1 = i_0,i_0+6*2-1,6

                   u_y0 = u_y(i_1+0)
                   u_y1 = u_y(i_1+1)
                   u_y2 = u_y(i_1+2)
                   u_y3 = u_y(i_1+3)
                   u_y4 = u_y(i_1+4)
                   u_y5 = u_y(i_1+5)

!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
!DIR$ UNROLL(8)
                   do j_1=kk_1,kk_4

                      u0 = u_x(j_1+0)

                      a0_0 = a(j_1+0,i_1+0)
                      a0_1 = a(j_1+0,i_1+1)

                      v_t(i_1+0) = v_t(i_1+0)
     &                            + (a0_0*u0)
                      v_t(i_1+1) = v_t(i_1+1)
     &                            + (a0_1*u0)

                   enddo

!              call MM_PREFETCH(a(kk_1,i_1),1)

!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
!DIR$ UNROLL(4)
                   do j_1=kk_1,kk_4

                      u0 = u_x(j_1+0)
                      v0 = ZERO

                      a0_0 = a(j_1+0,i_1+0)
                      a0_1 = a(j_1+0,i_1+1)

                      v0 = v0
     &                            + (a0_0*u_y0)
     &                            + (a0_1*u_y1)

                      a0_2 = a(j_1+0,i_1+2)
                      a0_3 = a(j_1+0,i_1+3)

                      v_t(i_1+2) = v_t(i_1+2)
     &                            + (a0_2*u0)
                      v_t(i_1+3) = v_t(i_1+3)
     &                            + (a0_3*u0)

                      v0 = v0
     &                            + (a0_2*u_y2)
     &                            + (a0_3*u_y3)

                      a0_4 = a(j_1+0,i_1+4)
                      a0_5 = a(j_1+0,i_1+5)

                      v_t(i_1+4) = v_t(i_1+4)
     &                            + (a0_4*u0)
                      v_t(i_1+5) = v_t(i_1+5)
     &                            + (a0_5*u0)

                      v0 = v0
     &                            + (a0_4*u_y4)
     &                            + (a0_5*u_y5)

                      u_t(j_1+0) = v0 + u_t(j_1+0)

                   end do! j_1

            enddo! i_1
            enddo! kk_1

                end do! i_1
             end do! k_1


       return
       end subroutine ! tred1_au_body1

       subroutine  tred1_au_body2(
     &                u_t, u_z, v_t, v_z, nv,
     &                x_pos, y_pos
     &                ,local_rank, local_size
     &                )
       implicit NONE

       integer, intent(in)    :: nv, x_pos, y_pos
       real(8), intent(out)   :: u_t(1:nv),   v_t(1:nv)
       real(8), intent(in)    :: u_z(1:nv,*), v_z(1:nv,*)

       integer                :: local_rank, local_size

       integer                :: i_1, i_2, i_3, i_4
       integer                :: j_1, j_2, j_3, j_4
       integer                :: jj_1, jj_2, jj_3, jj_4
       integer                :: i, j, k

       integer, parameter     :: LX = 1024

       include 'param.h'


          jj_1 = x_pos
          jj_2 = MAX(512,(jj_1-1)/local_size+1)
          jj_3 =    (jj_2*(local_rank+0)     )+1
          jj_4 = MIN(jj_2*(local_rank+1),jj_1)

          do jj_1=jj_3,jj_4,LX
             j_3=jj_1; j_4=MIN(jj_1+LX-1,jj_4)
           if ( local_size == 4 ) then
                do j_1=j_3,j_4
                  u_t(j_1) = u_z(j_1,1)+u_z(j_1,2)+u_z(j_1,3)+u_z(j_1,4)
                end do
           else
             if ( MOD(local_size,2) == 1 ) then
                do j_1=j_3,j_4
                   u_t(j_1) = u_z(j_1,1)
                end do
             else
                do j_1=j_3,j_4
                   u_t(j_1) = ZERO
                end do
             end if
             do j=MOD(local_size,2)+1,local_size,2
                do j_1=j_3,j_4
                   u_t(j_1) = u_t(j_1)+u_z(j_1,j+0)+u_z(j_1,j+1)
                end do
             end do
            end if
          end do

          jj_1 = y_pos
          jj_2 = MAX((jj_1-1)/local_size+1,512)
          jj_3 =    (jj_2*(local_rank+0)     )+1
          jj_4 = MIN(jj_2*(local_rank+1),jj_1)

          do jj_1=jj_3,jj_4,LX
             j_3=jj_1; j_4=MIN(jj_1+LX-1,jj_4)
           if ( local_size == 4 ) then
                do j_1=j_3,j_4
                  v_t(j_1) = v_z(j_1,1)+v_z(j_1,2)+v_z(j_1,3)+v_z(j_1,4)
                end do
           else
             if ( MOD(local_size,2) == 1 ) then
                do j_1=j_3,j_4
                   v_t(j_1) = v_z(j_1,1)
                end do
             else
                do j_1=j_3,j_4
                   v_t(j_1) = ZERO
                end do
             end if
             do j=MOD(local_size,2)+1,local_size,2
                do j_1=j_3,j_4
                   v_t(j_1) = v_t(j_1)+v_z(j_1,j+0)+v_z(j_1,j+1)
                end do
             end do
           end if
          end do


       return
       end subroutine ! tred1_au_body2

       subroutine  tred1_au_body3(
     &                u_x, u_y, v_x,
     &                u_t,v_t,d_t,
     &                n1, n2, x_pos, y_pos, nv
     &                )
       implicit NONE

       integer, intent(in)    :: nv, n1, n2
       real(8), intent(in)    :: u_x(1:nv),u_y(1:nv)
       real(8), intent(inout) :: v_x(1:nv)
       real(8), intent(in)    :: u_t(*),v_t(*)
       real(8), intent(in)    :: d_t(1:nv)
       integer                :: x_pos, y_pos, x_root, y_root

       integer                :: i_1, i_2, i_3, i_4
       integer                :: j_1, j_2, j_3, j_4
       integer                :: i, j, k
       integer                :: nm1, nm2

       include 'trd.h'


          i_2 = n1
          i_3 = n2

          if ( diag_0 > 0 ) then

             j = loop_c_l2g_depth(diag_0, y_nnod,y_inod)
             j = loop_c_g2l_depth(j, x_nnod,x_inod)
             if ( j > nv ) return

             nm1 = y_nnod/n_common
             nm2 = x_nnod/n_common

             if ( nm2 == 1 ) then
             if ( nm1 == 1 ) then
                call tred1_au_body3_sub11(
     &               v_x(j), v_t(diag_0), d_t(diag_0), u_y(diag_0),
     &               (i_3-diag_0)/(x_nnod/n_common)+1, nm1, nm2)
             else
                call tred1_au_body3_subX1(
     &               v_x(j), v_t(diag_0), d_t(diag_0), u_y(diag_0),
     &               (i_3-diag_0)/(x_nnod/n_common)+1, nm1, nm2)
             end if
             else
                call tred1_au_body3_subXX(
     &               v_x(j), v_t(diag_0), d_t(diag_0), u_y(diag_0),
     &               (i_3-diag_0)/(x_nnod/n_common)+1, nm1, nm2)

             end if

          end if


       return
       end subroutine ! tred1_au_body3

       subroutine tred1_au_body3_subXX(v_x,v_t,d_t,u_y, n,nm1,nm2)
       implicit NONE

       integer, intent(in)    :: n, nm1, nm2
       real(8), intent(inout) :: v_x(nm1,*)
       real(8), intent(in)    :: v_t(nm2,*)
       real(8), intent(in)    :: d_t(nm2,*)
       real(8), intent(in)    :: u_y(nm2,*)

       integer                :: i


!DIR$ VECTOR ALWAYS
          do i=1,n
             v_x(1,i) = v_x(1,i)+v_t(1,i)+d_t(1,i)*u_y(1,i)
          end do! i


       return
       end subroutine ! tred1_au_body3_subXX

       subroutine tred1_au_body3_subX1(v_x,v_t,d_t,u_y, n,nm1,nm2)
       implicit NONE

       integer, intent(in)    :: n, nm1, nm2
       real(8), intent(inout) :: v_x(nm1,*)
       real(8), intent(in)    :: v_t(*)
       real(8), intent(in)    :: d_t(*)
       real(8), intent(in)    :: u_y(*)

       integer                :: i


!DIR$ VECTOR ALWAYS
          do i=1,n
             v_x(1,i) = v_x(1,i)+v_t(i)+d_t(i)*u_y(i)
          end do! i


       return
       end subroutine ! tred1_au_body3_subX1

       subroutine tred1_au_body3_sub11(v_x,v_t,d_t,u_y, n,nm1,nm2)
       implicit NONE

       integer, intent(in)    :: n, nm1, nm2
       real(8), intent(inout) :: v_x(1:n)
       real(8), intent(in)    :: v_t(1:n)
       real(8), intent(in)    :: d_t(1:n)
       real(8), intent(in)    :: u_y(1:n)

       integer                :: i


!DIR$ VECTOR ALWAYS
          do i=1,n
             v_x(i) = v_x(i)+v_t(i)+d_t(i)*u_y(i)
          end do! i


       return
       end subroutine ! tred1_au_body3_sub11

