      subroutine  eigen_hrd_Au(
     &              ar, ai, nm,
     &              ur_x, ui_x,
     &              ur_y, ui_y,
     &              vr_x, vi_x,
     &              ur_t, ui_t, vr_t, vi_t,
     &              d_t,
     &              n, x_pos,y_pos, nv,blk_0,
     &              ar_,ai_,n_,nm_
     &            )
!----
      use communication_h, only : send_dbl
     &               , recv_dbl
     &               , reduce_dbl
!----
      implicit none

      integer, intent(in)    :: nm, nv, n
      real(8), intent(inout) :: ar(1:nm,*), ai(1:nm,*)
      real(8), intent(inout) :: ur_x(1:nv), ui_x(1:nv)
      real(8), intent(inout) :: ur_y(1:nv), ui_y(1:nv)
      real(8), intent(inout) :: vr_x(1:nv), vi_x(1:nv)
      real(8), intent(inout) :: vr_t(1:nv), vi_t(1:nv)
      real(8), intent(inout) :: ur_t(1:nv), ui_t(1:nv)
      real(8), intent(in)    :: d_t(1:nv)
      integer, intent(in)    :: blk_0
      integer, intent(in)    :: n_,nm_
      real(8), intent(inout) :: ar_(1:nm_,*), ai_(1:nm_,*)
      integer                :: x_pos, y_pos
      integer                :: local_size, local_rank
      real(8), pointer       :: ur_z(:),ui_z(:),vr_z(:),vi_z(:)
      save                   :: ur_z,   ui_z,   vr_z   ,vi_z

      include 'mpif.h'
!$    include 'omp_lib.h'
      include 'trd.h'

      local_rank = 0
      local_size = 1
!$    local_rank = omp_get_thread_num()
!$    local_size = omp_get_num_threads()

      if ( n == -1 ) then
         if ( local_rank == 0 ) then
            allocate(ur_z(nv*local_size), ui_z(nv*local_size),
     &               vr_z(nv*local_size), vi_z(nv*local_size))
         end if
         return
      end if
      if ( n <  -1 ) then
         if ( local_rank == 0 ) then
            deallocate(ur_z, ui_z, vr_z, vi_z)
         end if
         return
      end if

      call eigen_hrd_Au_step1(
     &       ar,ai, nm,
     &       ur_x,ui_x, ur_y,ui_y,
     &       ur_z(1+nv*local_rank),
     &       ui_z(1+nv*local_rank),
     &       vr_z(1+nv*local_rank),
     &       vi_z(1+nv*local_rank),
     &       1,n,x_pos,y_pos,nv,blk_0,
     &       local_rank, local_size
     &     )

!$omp barrier

      call eigen_hrd_Au_step2(
     &       vr_x,vi_x, ur_z,ui_z,
     &       vr_t,vi_t, vr_z,vi_z,
     &       nv, x_pos, y_pos,
     &       local_rank, local_size
     &     )

!$omp barrier

      if ( local_rank == 0 ) then

         if(n+1<=y_pos)then
            vr_t(n+1:y_pos) = 0.0d0
            vi_t(n+1:y_pos) = 0.0d0
         endif
         if ( nprocs > 1 ) then
            call reduce_dbl(vr_t, ur_t, y_pos, mpi_comm_col)
            call reduce_dbl(vi_t, ui_t, y_pos, mpi_comm_col)
         end if

         call eigen_hrd_Au_step3(
     &          ur_x,ui_x,ur_y,ui_y,vr_x,vi_x,
     &          ur_t,ui_t,vr_t,vi_t,d_t,
     &          1,n,x_pos,y_pos,nv
     &        )

         if ( nprocs > 1 ) then
            if ( nprocs >= es_cluster**2 .and.
     &                       size_of_col > 1 .and. es_icnod > 1 ) then
               call recv_dbl(vr_x(1), 0, es_icnod-1, es_ccomm_world)
            endif
            call reduce_dbl(vr_x, vr_t, x_pos, mpi_comm_row)
            call reduce_dbl(vi_x, vi_t, x_pos, mpi_comm_row)
            if ( nprocs >= es_cluster**2 .and.
     &           size_of_col > 1 .and. 
     &           es_icnod < es_ncnod ) then
               call send_dbl(vr_x(1), 0, es_icnod+1, es_ccomm_world)
            endif
         end if

      endif

      return
      end subroutine ! eigen_hrd_Au

      subroutine  eigen_hrd_Au_step1(
     &              ar,ai,nm,
     &              ur_x,ui_x,ur_y,ui_y,
     &              ur_t,ui_t,vr_t,vi_t,
     &              n1,n2,x_pos,y_pos,nv,blk_0,
     &              local_rank, local_size
     &            )
!----
      use communication_h, only : get_loop_end
     &               , translate_l2g
!----
      implicit none

      integer, intent(in)    :: nm, nv, n1, n2
      real(8), intent(inout) :: ar(1:nm,*), ai(1:nm,*)
      real(8), intent(inout) :: ur_x(1:nv), ui_x(1:nv)
      real(8), intent(inout) :: ur_y(1:nv), ui_y(1:nv)
      real(8), intent(inout) :: ur_t(1:nv), ui_t(1:nv)
      real(8), intent(inout) :: vr_t(1:nv), vi_t(1:nv)
      integer, intent(in)    :: blk_0
      integer                :: x_pos, y_pos
      integer                :: local_size, local_rank
      integer                :: i_1, i_2, i_3, i_4
      integer                :: j_1, j_2, j_3, j_4
      integer                :: k_1, k_2, k_3, k_4
      integer                :: i, j, k

      real(8)                :: vr0, vi0, ur0, ui0
      real(8)                :: ar0_0, ai0_0
      real(8)                :: ar0_1, ai0_1
      real(8)                :: ar0_2, ai0_2
      real(8)                :: ar0_3, ai0_3
      real(8)                :: ar0_4, ai0_4
      real(8)                :: ar0_5, ai0_5
      real(8)                :: ar0_6, ai0_6
      real(8)                :: ar0_7, ai0_7
      real(8)                :: ar0_8, ai0_8
      real(8)                :: ar0_9, ai0_9
      real(8)                :: vr_t0,vi_t0, ur_y0,ui_y0
      real(8)                :: vr_t1,vi_t1, ur_y1,ui_y1
      real(8)                :: vr_t2,vi_t2, ur_y2,ui_y2
      real(8)                :: vr_t3,vi_t3, ur_y3,ui_y3
      real(8)                :: vr_t4,vi_t4, ur_y4,ui_y4
      real(8)                :: vr_t5,vi_t5, ur_y5,ui_y5
      real(8)                :: vr_t6,vi_t6, ur_y6,ui_y6
      real(8)                :: vr_t7,vi_t7, ur_y7,ui_y7
      real(8)                :: vr_t8,vi_t8, ur_y8,ui_y8
      real(8)                :: vr_t9,vi_t9, ur_y9,ui_y9

      integer                ::   ii_1, ii_2, ii_3, ii_4
      integer                ::   jj_1, jj_2, jj_3, jj_4

      include 'mpif.h'
!$    include 'omp_lib.h'
      include 'trd.h'

      i_2 = n1
      i_3 = n2
!----
! v:= au
!----
      if ( n1 == 1 ) then
         ur_t(1:x_pos) = 0.0d+00
         ui_t(1:x_pos) = 0.0d+00
         vr_t(1:y_pos) = 0.0d+00
         vi_t(1:y_pos) = 0.0d+00
      end if

      if ( blk_0 == 0 ) then
         do i_1=i_2+local_rank,i_3,local_size
            j   = translate_l2g(i_1, size_of_row,my_row)
            j_3 = get_loop_end  (j, size_of_col,my_col)
            j   = j+size_of_row*(10-1)
            j_4 = get_loop_end  (j, size_of_col,my_col)
            do j_1=j_3+1,min(j_4,nm)
               ar(j_1,i_1) = 0.0d+00
               ai(j_1,i_1) = 0.0d+00
            end do! j_1
         end do! i_1
!$omp barrier
      end if
!----
      k_1 = 1
      ii_2 = i_2; ii_3 = i_3
      ii_4 = mod(ii_3-ii_2+1,10)+ii_2

      do i_1=ii_2+local_rank,ii_4-1,local_size

         j    = translate_l2g(i_1, size_of_row,my_row)
         j    = j+(1-1)*size_of_row
         jj_2 = k_1
         jj_3 = get_loop_end  (j-1, size_of_col,my_col)
         if ( jj_2 > jj_3 ) cycle

         do j_1=jj_2,jj_3

            ur0 = ur_x(j_1+0)
            ui0 = ui_x(j_1+0)
            vr0 = 0
            vi0 = 0
!----       vr0 = ur_t(j_1+0)
!----       vi0 = ui_t(j_1+0)

            ar0_0 = ar(j_1+0,i_1+0)
            ai0_0 = ai(j_1+0,i_1+0)
            vr_t(i_1+0) = vr_t(i_1+0)
     &                  + (ar0_0*ur0)
     &                  + (ai0_0*ui0)

            vi_t(i_1+0) = vi_t(i_1+0)
     &                  + (ar0_0*ui0)
     &                  - (ai0_0*ur0)
            vr0 = vr0
     &          + (ar0_0*ur_y(i_1+0))
     &          - (ai0_0*ui_y(i_1+0))
            vi0 = vi0
     &          + (ar0_0*ui_y(i_1+0))
     &          + (ai0_0*ur_y(i_1+0))

            ur_t(j_1+0) = vr0 + ur_t(j_1+0)
            ui_t(j_1+0) = vi0 + ui_t(j_1+0)

         end do! j_1

      end do! i_1

      do i_1=ii_4+local_rank*10,ii_3, local_size*10

         j    = translate_l2g(i_1, size_of_row,my_row)
         j    = j+(10-1)*size_of_row
         jj_2 = k_1
         jj_3 = get_loop_end  (j-1, size_of_col,my_col)
         if ( jj_2 > jj_3 ) cycle

         do j_1=jj_2,jj_3

            ur0 = ur_x(j_1+0)
            ui0 = ui_x(j_1+0)
            vr0 = 0
            vi0 = 0
!----       vr0 = ur_t(j_1+0)
!----       vi0 = ui_t(j_1+0)

            ar0_0 = ar(j_1+0,i_1+0)
            ai0_0 = ai(j_1+0,i_1+0)
            ar0_1 = ar(j_1+0,i_1+1)
            ai0_1 = ai(j_1+0,i_1+1)
            ar0_2 = ar(j_1+0,i_1+2)
            ai0_2 = ai(j_1+0,i_1+2)
            ar0_3 = ar(j_1+0,i_1+3)
            ai0_3 = ai(j_1+0,i_1+3)
            vr_t(i_1+0) = vr_t(i_1+0)
     &                  + (ar0_0*ur0)
     &                  + (ai0_0*ui0)

            vi_t(i_1+0) = vi_t(i_1+0)
     &                  + (ar0_0*ui0)
     &                  - (ai0_0*ur0)
            vr_t(i_1+1) = vr_t(i_1+1)
     &                  + (ar0_1*ur0)
     &                  + (ai0_1*ui0)

            vi_t(i_1+1) = vi_t(i_1+1)
     &                  + (ar0_1*ui0)
     &                  - (ai0_1*ur0)
            vr_t(i_1+2) = vr_t(i_1+2)
     &                  + (ar0_2*ur0)
     &                  + (ai0_2*ui0)

            vi_t(i_1+2) = vi_t(i_1+2)
     &                  + (ar0_2*ui0)
     &                  - (ai0_2*ur0)
            vr_t(i_1+3) = vr_t(i_1+3)
     &                  + (ar0_3*ur0)
     &                  + (ai0_3*ui0)

            vi_t(i_1+3) = vi_t(i_1+3)
     &                  + (ar0_3*ui0)
     &                  - (ai0_3*ur0)
            vr0 = vr0
     &                  + (ar0_0*ur_y(i_1+0))
     &                  - (ai0_0*ui_y(i_1+0))
     &                  + (ar0_1*ur_y(i_1+1))
     &                  - (ai0_1*ui_y(i_1+1))
     &                  + (ar0_2*ur_y(i_1+2))
     &                  - (ai0_2*ui_y(i_1+2))
     &                  + (ar0_3*ur_y(i_1+3))
     &                  - (ai0_3*ui_y(i_1+3))
            vi0 = vi0
     &                  + (ar0_0*ui_y(i_1+0))
     &                  + (ai0_0*ur_y(i_1+0))
     &                  + (ar0_1*ui_y(i_1+1))
     &                  + (ai0_1*ur_y(i_1+1))
     &                  + (ar0_2*ui_y(i_1+2))
     &                  + (ai0_2*ur_y(i_1+2))
     &                  + (ar0_3*ui_y(i_1+3))
     &                  + (ai0_3*ur_y(i_1+3))
            ar0_4 = ar(j_1+0,i_1+4)
            ai0_4 = ai(j_1+0,i_1+4)
            ar0_5 = ar(j_1+0,i_1+5)
            ai0_5 = ai(j_1+0,i_1+5)
            ar0_6 = ar(j_1+0,i_1+6)
            ai0_6 = ai(j_1+0,i_1+6)
            ar0_7 = ar(j_1+0,i_1+7)
            ai0_7 = ai(j_1+0,i_1+7)
            vr_t(i_1+4) = vr_t(i_1+4)
     &                  + (ar0_4*ur0)
     &                  + (ai0_4*ui0)

            vi_t(i_1+4) = vi_t(i_1+4)
     &                  + (ar0_4*ui0)
     &                  - (ai0_4*ur0)
            vr_t(i_1+5) = vr_t(i_1+5)
     &                  + (ar0_5*ur0)
     &                  + (ai0_5*ui0)

            vi_t(i_1+5) = vi_t(i_1+5)
     &                  + (ar0_5*ui0)
     &                  - (ai0_5*ur0)
            vr_t(i_1+6) = vr_t(i_1+6)
     &                  + (ar0_6*ur0)
     &                  + (ai0_6*ui0)

            vi_t(i_1+6) = vi_t(i_1+6)
     &                  + (ar0_6*ui0)
     &                  - (ai0_6*ur0)
            vr_t(i_1+7) = vr_t(i_1+7)
     &                  + (ar0_7*ur0)
     &                  + (ai0_7*ui0)

            vi_t(i_1+7) = vi_t(i_1+7)
     &                  + (ar0_7*ui0)
     &                  - (ai0_7*ur0)
            vr0 = vr0
     &                  + (ar0_4*ur_y(i_1+4))
     &                  - (ai0_4*ui_y(i_1+4))
     &                  + (ar0_5*ur_y(i_1+5))
     &                  - (ai0_5*ui_y(i_1+5))
     &                  + (ar0_6*ur_y(i_1+6))
     &                  - (ai0_6*ui_y(i_1+6))
     &                  + (ar0_7*ur_y(i_1+7))
     &                  - (ai0_7*ui_y(i_1+7))
            vi0 = vi0
     &                  + (ar0_4*ui_y(i_1+4))
     &                  + (ai0_4*ur_y(i_1+4))
     &                  + (ar0_5*ui_y(i_1+5))
     &                  + (ai0_5*ur_y(i_1+5))
     &                  + (ar0_6*ui_y(i_1+6))
     &                  + (ai0_6*ur_y(i_1+6))
     &                  + (ar0_7*ui_y(i_1+7))
     &                  + (ai0_7*ur_y(i_1+7))
            ar0_8 = ar(j_1+0,i_1+8)
            ai0_8 = ai(j_1+0,i_1+8)
            ar0_9 = ar(j_1+0,i_1+9)
            ai0_9 = ai(j_1+0,i_1+9)
            vr_t(i_1+8) = vr_t(i_1+8)
     &                  + (ar0_8*ur0)
     &                  + (ai0_8*ui0)

            vi_t(i_1+8) = vi_t(i_1+8)
     &                  + (ar0_8*ui0)
     &                  - (ai0_8*ur0)
            vr_t(i_1+9) = vr_t(i_1+9)
     &                  + (ar0_9*ur0)
     &                  + (ai0_9*ui0)

            vi_t(i_1+9) = vi_t(i_1+9)
     &                  + (ar0_9*ui0)
     &                  - (ai0_9*ur0)
            vr0 = vr0
     &                  + (ar0_8*ur_y(i_1+8))
     &                  - (ai0_8*ui_y(i_1+8))
     &                  + (ar0_9*ur_y(i_1+9))
     &                  - (ai0_9*ui_y(i_1+9))
            vi0 = vi0
     &                  + (ar0_8*ui_y(i_1+8))
     &                  + (ai0_8*ur_y(i_1+8))
     &                  + (ar0_9*ui_y(i_1+9))
     &                  + (ai0_9*ur_y(i_1+9))

            ur_t(j_1+0) = vr0 + ur_t(j_1+0)
            ui_t(j_1+0) = vi0 + ui_t(j_1+0)

         end do! j_1

      end do! i_1

      return
      end subroutine ! eigen_hrd_Au_step1

      subroutine  eigen_hrd_Au_step2(
     &            ur_t,ui_t, ur_z,ui_z,
     &            vr_t,vi_t, vr_z,vi_z, nv, x_pos, y_pos
     &           ,local_rank, local_size
     &            )
      implicit none

      integer, intent(in)    :: nv, x_pos, y_pos
      real(8), intent(out)   :: ur_t(1:nv),   ui_t(1:nv)
      real(8), intent(out)   :: vr_t(1:nv),   vi_t(1:nv)
      real(8), intent(in)    :: ur_z(1:nv,*), ui_z(1:nv,*)
      real(8), intent(in)    :: vr_z(1:nv,*), vi_z(1:nv,*)

      integer                :: local_rank, local_size
      integer                :: j_1, j_2, j_3, j_4
      integer                :: jj_1, jj_2, jj_3, jj_4
      integer                :: i, j, k
      integer, parameter     :: lx = 1024

      jj_1 = x_pos
      jj_2 = (jj_1-1)/local_size+1
      jj_3 =    (jj_2*(local_rank+0)     )+1
      jj_4 = min(jj_2*(local_rank+1),jj_1)

      if ( jj_2 < 1024 ) then
         if(local_rank==0)then
            jj_3=1; jj_4=x_pos
         else
            jj_3=1; jj_4=0
         endif
      endif

      j_3=jj_3; j_4=jj_4
      if ( mod(local_size,2) == 1 ) then
         do j_1=j_3,j_4
            ur_t(j_1) = ur_z(j_1,1)
            ui_t(j_1) = ui_z(j_1,1)
         end do
      else
         do j_1=j_3,j_4
            ur_t(j_1) = ur_z(j_1,1)+ur_z(j_1,2)
            ui_t(j_1) = ui_z(j_1,1)+ui_z(j_1,2)
         end do
      end if
      do j=mod(local_size-1,2)+2,local_size,2
         do j_1=j_3,j_4
            ur_t(j_1) = ur_t(j_1)+ur_z(j_1,j+0)+ur_z(j_1,j+1)
            ui_t(j_1) = ui_t(j_1)+ui_z(j_1,j+0)+ui_z(j_1,j+1)
         end do
      end do

      jj_1 = y_pos
      jj_2 = (jj_1-1)/local_size+1
      jj_3 =    (jj_2*(local_rank+0)     )+1
      jj_4 = min(jj_2*(local_rank+1),jj_1)

      if ( jj_2 < 1024 ) then
         if(local_rank==0)then
            jj_3=1; jj_4=y_pos
         else
            jj_3=1; jj_4=0
         endif
      endif

      j_3=jj_3; j_4=jj_4
      if ( mod(local_size,2) == 1 ) then
         do j_1=j_3,j_4
            vr_t(j_1) = vr_z(j_1,1)
            vi_t(j_1) = vi_z(j_1,1)
         end do
      else
         do j_1=j_3,j_4
            vr_t(j_1) = vr_z(j_1,1)+vr_z(j_1,2)
            vi_t(j_1) = vi_z(j_1,1)+vi_z(j_1,2)
         end do
      end if
      do j=mod(local_size-1,2)+2,local_size,2
         do j_1=j_3,j_4
            vr_t(j_1) = vr_t(j_1)+vr_z(j_1,j+0)+vr_z(j_1,j+1)
            vi_t(j_1) = vi_t(j_1)+vi_z(j_1,j+0)+vi_z(j_1,j+1)
         end do
      end do

      return
      end subroutine ! eigen_hrd_Au_step2

      subroutine  eigen_hrd_Au_step3(
     &              ur_x,ui_x,ur_y,ui_y,
     &              vr_x,vi_x,
     &              ur_t,ui_t,vr_t,vi_t,d_t,
     &              n1,n2,x_pos,y_pos,nv
     &            )
!----
      use communication_h, only : translate_l2g
     &               , translate_g2l
!----
      implicit none

      integer, intent(in)    :: nv, n1, n2
      real(8), intent(in)    :: ur_x(1:nv),ui_x(1:nv)
      real(8), intent(in)    :: ur_y(1:nv),ui_y(1:nv)
      real(8), intent(inout) :: vr_x(1:nv),vi_x(1:nv)
      real(8), intent(in)    :: ur_t(*),ui_t(*),vr_t(*),vi_t(*)
      real(8), intent(in)    :: d_t(1:nv)
      integer                :: x_pos, y_pos

      integer                :: i_2, i_3
      integer                :: i, j, k
      integer                :: nm1, nm2

      include 'mpif.h'
      include 'trd.h'

      i_2 = n1
      i_3 = n2

      if ( diag_0 > 0 ) then

         j = translate_l2g(diag_0, size_of_row,my_row)
         j = translate_g2l(j, size_of_col,my_col)
         if ( j > nv ) return

         nm1 = size_of_row/n_common
         nm2 = size_of_col/n_common

         if ( nm2 == 1 ) then
            if ( nm1 == 1 ) then
               call eigen_hrd_Au_step3_sub3(
     &                vr_x(j), vi_x(j),
     &                vr_t(diag_0), vi_t(diag_0),
     &                d_t(diag_0),
     &                ur_y(diag_0), ui_y(diag_0),
     &                (i_3-diag_0)/(size_of_col/n_common)+1, nm1, nm2
     &              )
            else
               call eigen_hrd_Au_step3_sub2(
     &                vr_x(j), vi_x(j),
     &                vr_t(diag_0), vi_t(diag_0),
     &                d_t(diag_0),
     &                ur_y(diag_0), ui_y(diag_0),
     &                (i_3-diag_0)/(size_of_col/n_common)+1, nm1, nm2
     &              )
            end if
         else
            call eigen_hrd_Au_step3_sub1(
     &             vr_x(j), vi_x(j),
     &             vr_t(diag_0), vi_t(diag_0),
     &             d_t(diag_0),
     &             ur_y(diag_0), ui_y(diag_0),
     &             (i_3-diag_0)/(size_of_col/n_common)+1, nm1, nm2
     &           )

         end if

      end if

      return
      end subroutine ! eigen_hrd_Au_step3

      subroutine eigen_hrd_Au_step3_sub1(
     &             vr_x,vi_x,vr_t,vi_t,d_t,ur_y,ui_y,
     &             n,nm1,nm2
     &           )
      implicit none

      integer, intent(in)    :: n, nm1, nm2
      real(8), intent(inout) :: vr_x(nm1,*), vi_x(nm1,*)
      real(8), intent(in)    :: vr_t(nm2,*), vi_t(nm2,*)
      real(8), intent(in)    :: d_t(nm2,*)
      real(8), intent(in)    :: ur_y(nm2,*), ui_y(nm2,*)

      integer                :: i

      do i=1,n
         vr_x(1,i) = vr_x(1,i)+vr_t(1,i)+d_t(1,i)*ur_y(1,i)
         vi_x(1,i) = vi_x(1,i)+vi_t(1,i)+d_t(1,i)*ui_y(1,i)
      end do! i

      return
      end subroutine ! eigen_hrd_Au_step3_sub1

      subroutine eigen_hrd_Au_step3_sub2(
     &             vr_x,vi_x,vr_t,vi_t,d_t,ur_y,ui_y,
     &             n,nm1,nm2
     &           )
      implicit none

      integer, intent(in)    :: n, nm1, nm2
      real(8), intent(inout) :: vr_x(nm1,*), vi_x(nm1,*)
      real(8), intent(in)    :: vr_t(*), vi_t(*)
      real(8), intent(in)    :: d_t(*)
      real(8), intent(in)    :: ur_y(*), ui_y(*)

      integer                :: i

      do i=1,n
         vr_x(1,i) = vr_x(1,i)+vr_t(i)+d_t(i)*ur_y(i)
         vi_x(1,i) = vi_x(1,i)+vi_t(i)+d_t(i)*ui_y(i)
      end do! i

      return
      end subroutine ! eigen_hrd_Au_step3_sub2

      subroutine eigen_hrd_Au_step3_sub3(
     &             vr_x,vi_x,vr_t,vi_t,d_t,ur_y,ui_y,
     &             n,nm1,nm2
     &           )
      implicit none

      integer, intent(in)    :: n, nm1, nm2
      real(8), intent(inout) :: vr_x(1:n), vi_x(1:n)
      real(8), intent(in)    :: vr_t(1:n), vi_t(1:n)
      real(8), intent(in)    :: d_t(1:n)
      real(8), intent(in)    :: ur_y(1:n), ui_y(1:n)

      integer                :: i

      do i=1,n
         vr_x(i) = vr_x(i)+vr_t(i)+d_t(i)*ur_y(i)
         vi_x(i) = vi_x(i)+vi_t(i)+d_t(i)*ui_y(i)
      end do! i

      return
      end subroutine ! eigen_hrd_Au_step3_sub3
