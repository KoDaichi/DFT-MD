!--
       subroutine send_dbl( buf, n, idest, icom )
       use MPI
       implicit NONE

       integer, intent(in)    :: n, idest, icom
       real(8), intent(inout) :: buf(1:n)

       integer                :: ierr


          call MPI_Send( buf, n, MPI_DOUBLE_PRECISION,
     &                  idest-1, 1, icom, ierr )


       return
       end subroutine  send_dbl
!--
       subroutine send_dblt( buf, n, idest, itag, icom )
       use MPI
       implicit NONE

       integer, intent(in)    :: n, idest, itag, icom
       real(8), intent(inout) :: buf(1:n)

       integer                :: ierr


          call MPI_Send( buf, n, MPI_DOUBLE_PRECISION,
     &                  idest-1, itag, icom, ierr )


       return
       end subroutine  send_dblt
!--
       subroutine isend_dbl( buf, n, idest, ireq, icom )
       use MPI
       implicit NONE

       integer, intent(in)    :: n, idest, icom
       integer, intent(inout) :: ireq
       real(8), intent(inout) :: buf(1:n)

       integer                :: ierr


          call MPI_Isend( buf, n, MPI_DOUBLE_PRECISION,
     &                   idest-1, 1, icom, ireq, ierr )


       return
       end subroutine  isend_dbl
!--
       subroutine isend_dblt( buf, n, idest, itag, ireq, icom )
       use MPI
       implicit NONE

       integer, intent(in)    :: n, idest, itag, icom
       integer, intent(inout) :: ireq
       real(8), intent(inout) :: buf(1:n)

       integer                :: ierr


          call MPI_Isend( buf, n, MPI_DOUBLE_PRECISION,
     &                   idest-1, itag, icom, ireq, ierr )


       return
       end subroutine  isend_dblt
!--
       subroutine recv_dbl( buf, n, isrc, icom )
       use MPI
       implicit NONE

       integer, intent(in)    :: n, isrc, icom
       real(8), intent(inout) :: buf(1:n)

       integer                :: ierr


          call MPI_Recv( buf, n, MPI_DOUBLE_PRECISION,
     &                  isrc-1, 1, icom, MPI_STATUS_IGNORE, ierr )


       return
       end subroutine  recv_dbl
!--
       subroutine recv_dblt( buf, n, isrc, itag, icom )
       use MPI
       implicit NONE

       integer, intent(in)    :: n, isrc, itag, icom
       real(8), intent(inout) :: buf(1:n)

       integer                :: ierr


          call MPI_Recv( buf, n, MPI_DOUBLE_PRECISION,
     &                  isrc-1, itag, icom, MPI_STATUS_IGNORE, ierr )


       return
       end subroutine  recv_dblt
!--
       subroutine irecv_dbl( buf, n, isrc, ireq, icom )
       use MPI
       implicit NONE

       integer, intent(in)    :: n, isrc, icom
       integer, intent(inout) :: ireq
       real(8), intent(inout) :: buf(1:n)

       integer                :: ierr


          call MPI_Irecv( buf, n, MPI_DOUBLE_PRECISION,
     &                   isrc-1, 1, icom, ireq, ierr )


       return
       end subroutine  irecv_dbl
!--
       subroutine irecv_dblt( buf, n, isrc, itag, ireq, icom )
       use MPI
       implicit NONE

       integer, intent(in)    :: n, isrc, itag, icom
       integer, intent(inout) :: ireq
       real(8), intent(inout) :: buf(1:n)

       integer                :: ierr


          call MPI_Irecv( buf, n, MPI_DOUBLE_PRECISION,
     &                   isrc-1, itag, icom, ireq, ierr )


       return
       end subroutine  irecv_dblt
!--
       subroutine wait_dbl( ireq )
       use MPI
       implicit NONE

       integer, intent(inout) :: ireq

       integer                :: ierr


          call MPI_Wait( ireq, MPI_STATUS_IGNORE, ierr )


       return
       end subroutine  wait_dbl
!--
       subroutine waitall_dbl( n, ireq )
       use MPI
       implicit NONE

       integer, intent(in   ) :: n
       integer, intent(inout) :: ireq(n)

       integer                :: ierr


          call MPI_Waitall( n, ireq, MPI_STATUSES_IGNORE, ierr )


       return
       end subroutine  waitall_dbl
!--
       subroutine barrier( icom )
       use MPI
       implicit NONE

       integer, intent(in)    :: icom

       integer                :: ierr


          call MPI_Barrier( icom, ierr )


       return
       end subroutine  barrier
!--
       subroutine bcast_dbl( buf, n, iroot, icom )
       use MPI
       use eigen_devel
       implicit NONE

       integer, intent(in)    :: n, iroot, icom
       real(8), intent(inout) :: buf(1:n)

       integer                :: ierr, i, j, k, tag
       integer, parameter     :: BCAST_TAG = 100000

       real(8) :: d1,d2

       integer :: my_rank, world_size, local_rank
       integer, pointer :: ireq(:)


          d1 = eigen_get_wtime( )

! On Fujitsu K and FX10 computers, the following message dividing
! on broadcasting reduces the unexpectd communication overhead
! when the message size exeeds 4096, somehow.

#if defined(__INTEL_COMPILER)
          call MPI_Bcast( buf, n, MPI_DOUBLE_PRECISION,
     &                   iroot-1, icom, ierr )
#else
       if ( icom == TRD_COMM_WORLD .OR. n <= 4096 ) then

          do i=1,n,4096
             j = MIN(n-i+1, 4096)
             call MPI_Bcast( buf(i), j, MPI_DOUBLE_PRECISION,
     &                      iroot-1, icom, ierr )
          end do

       else

          call MPI_Comm_rank( icom, my_rank, ierr )
          my_rank = my_rank + 1
          call MPI_Comm_size( icom, world_size, ierr )

#if 0
!
! Sequential algorithm
!
          if ( iroot == my_rank ) then
             allocate( ireq(1:world_size) )
             j = 0
             do i=1, world_size
                if ( i /= iroot ) then
                   j = j + 1
                   call isend_dblt( buf, n, i, i, ireq(j), icom )
                end if
             end do
             if ( j > 0 ) then
                call waitall_dbl( j, ireq )
             end if
             deallocate( ireq )
          else
             call recv_dblt( buf, n, iroot, my_rank, icom )
          end if
#else
!
! Binomial algorithm
!
          local_rank = MOD(world_size+my_rank-iroot, world_size) ! relative rank
          allocate( ireq(1:world_size) )
          i=1; j=0
          do
             if ( i >= world_size ) exit
             if ( local_rank < i ) then
                if ( local_rank+i < world_size ) then
                   k = MOD(my_rank+i-1, world_size)+1
                   j = j + 1
                   tag = BCAST_TAG + k
                   call isend_dblt( buf, n, k, tag, ireq(j), icom )
                end if
             else if ( local_rank < 2*i ) then
                   k = MOD(world_size+my_rank-i-1, world_size)+1
                   tag = BCAST_TAG + my_rank
                   call recv_dblt( buf, n, k, tag, icom )
             end if
             i = i * 2
          end do
          if ( j > 0 ) then
             call waitall_dbl( j, ireq )
          end if
          deallocate( ireq )
#endif

       end if
#endif

          d2 = eigen_get_wtime( )
          time_bcast = time_bcast + (d2-d1)


       return
       end subroutine  bcast_dbl
!--
       subroutine bcastw_dbl( buf, n, iroot, lda, px, buffer, icom )
       use MPI
       use eigen_devel
       implicit NONE

       integer, intent(in)    :: n, iroot, lda, px, icom
       real(8), intent(inout) :: buf(1:lda, 1:px)
       real(8), intent(out)   :: buffer(1:n, *)

       integer                :: ierr, i, j, k, tag
       integer, parameter     :: BCASTW_TAG = 200000

       real(8) :: d1,d2
       real(8) :: t1,t2

       integer :: my_rank, world_size, local_rank, new_type
       integer, pointer :: ireq(:)


          d1 = eigen_get_wtime( )
          if ( px == 1 ) then
             call bcast_dbl( buf, n, iroot, icom )
          else

             call MPI_Comm_size( icom, world_size, ierr )

             if ( n <= 4096 .OR. world_size < 64/px ) then
                do i=0,px-1
                   k = MOD((iroot)-1+(i), world_size)+1
                   call bcast_dbl( buf(1,i+1), n, k, icom )
                end do
                return
             end if

             call MPI_Comm_rank( icom, my_rank, ierr )
             my_rank = my_rank + 1

             allocate( ireq(1:world_size) )

             local_rank = MOD(world_size+my_rank-iroot, world_size)

             if ( local_rank ==0 ) then
                do i = 1, px-1
                   k = MOD(iroot+(i)-1, world_size)+1
                   tag = BCASTW_TAG + k
                   call irecv_dblt(buffer(1,i+1), n, k, tag,
     &                             ireq(i), icom)
                end do
                buffer(1:n, 1) = buf(1:n, 1)
                call waitall_dbl( px-1, ireq )
             else
                if ( local_rank < px ) then
                   i = local_rank + 1
                   tag = BCASTW_TAG + my_rank
                   call isend_dblt(buf(1,i), n, iroot, tag,
     &                             ireq(1), icom)
                end if
             end if

          t1 = eigen_get_wtime( )
             call bcast_dbl( buffer, n*px, iroot, icom )
          t2 = eigen_get_wtime( )

             do i=0, px-1
                if ( i /= local_rank ) then
                   buf(1:n, i+1) = buffer(1:n, i+1)
                end if
             end do

             if ( local_rank /=0 ) then
                if ( local_rank < px ) then
                   call wait_dbl( ireq(1) )
                end if
             end if

             deallocate( ireq )

          end if
          d2 = eigen_get_wtime( )
          time_bcast = time_bcast+(d2-d1)-(t2-t1)


       return
       end subroutine  bcastw_dbl
!--
       subroutine reduce_dbl( buf, wrk, n, dist, icom )
       use MPI
       use eigen_devel
       implicit NONE

       integer, intent(in)    :: n, dist, icom
       real(8), intent(inout) :: buf(1:n), wrk(1:n)

       integer                :: ierr, i, j

       real(8) :: d1,d2


          d1 = eigen_get_wtime( )

          call MPI_Allreduce( buf, wrk, n, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, icom, ierr )
          buf(1:n) = wrk(1:n)

          d2 = eigen_get_wtime( )
          time_reduce = time_reduce + (d2-d1)


       return
       end subroutine  reduce_dbl
!--
       subroutine allgather_dbl( buf, wrk, n, icom )
       use MPI
       use eigen_devel
       implicit NONE

       integer, intent(in)    :: n, icom
       real(8), intent(inout) :: buf(1:n), wrk(1:n)

       integer                :: ierr, i, j

       real(8) :: d1,d2


          d1 = eigen_get_wtime( )

          call MPI_Allgather(buf, n, MPI_DOUBLE_PRECISION,
     &                       wrk, n, MPI_DOUBLE_PRECISION,
     &                       icom, ierr)

          d2 = eigen_get_wtime( )
          time_gather = time_gather + (d2-d1)


       return
       end subroutine  allgather_dbl
!--
       subroutine datacast_dbl( u_y, u_x, u_t, u_s, n )
       use eigen_devel
       implicit NONE

       integer, intent(in)    :: n
       real(8), intent(inout) :: u_y(1:n), u_x(1:n), u_t(1:n), u_s(1:n)

       integer :: nx, ny, ic, i, j, k
       integer, pointer :: req(:)
       integer :: reqr(2), reqs(2), x_snod, y_snod
       integer :: his_rank, her_rank

       real(8) :: d1,d2
       real(8) :: e1,e2
       real(8) :: e3,e4


          d1 = eigen_get_wtime( )

          if ( x_nnod == 1 ) then
             if ( y_nnod == 1 ) then
                u_y(1:n) = u_x(1:n)
             else
                ny = (n-1)/y_nnod+1
                do i=1,ny
                   j = y_inod+y_nnod*(i-1)
                   u_y(i) = u_x(j)
                end do! i
             end if
             return
          end if

          if ( x_nnod == y_nnod ) then
             if ( x_inod == y_inod ) then
                u_y(1:n) = u_x(1:n)
             end if
             call bcast_dbl(u_y, n, y_inod, x_COMM_WORLD)
             return
          end if

          x_snod = x_nnod/n_common
          y_snod = y_nnod/n_common

          if ( p0_(x_inod) > 0 ) then

             nx = (n-1)/y_snod+1
             do i=1,nx
                j = p0_(x_inod)+y_snod*(i-1)
                k = q0_(x_inod)+x_snod*(i-1)
                u_t(i) = u_x(j)
                u_y(k) = u_x(j)
             end do

             e3 = eigen_get_wtime( )
             call allgather_dbl( u_t(1), u_s(1), nx, w_COMM_WORLD )
             e4 = eigen_get_wtime( )

             do ic=1,x_snod-1
                his_rank = MOD(x_inod-1 +x_nnod +ic*n_common,x_nnod)+1
                ny = (n-p0_(his_rank))/y_snod+1
                do i=1,ny
                   k = q0_(his_rank)+x_snod*(i-1)
                   j = (his_rank-1)*nx+i
                   u_y(k) = u_s(j)
                end do! i
             end do! ic

             e1 = eigen_get_wtime( )
             call bcast_dbl( u_y, n, z_inod, z_COMM_WORLD )
             e2 = eigen_get_wtime( )

          else

             i = MOD(y_inod-1,n_common)
             j = MOD(x_inod-1,n_common)
             ic = MOD(j-i+n_common,n_common)
             his_rank = ((x_inod-1)/n_common)*n_common + i + 1

             e1 = eigen_get_wtime( )
             call bcast_dbl( u_y, n, his_rank, z_COMM_WORLD )
             e2 = eigen_get_wtime( )

             e3 = 0; e4 = 0

          end if

9999      continue

          d2 = eigen_get_wtime( )
          time_redist = time_redist + (d2-d1) - (e2-e1) - (e4-e3)


       return
       end subroutine  datacast_dbl
!--
       subroutine datacast_dbl2(ur_y,ui_y, ur_x,ui_x, u_t,u_s, n)
       use eigen_devel
       implicit NONE

       integer, intent(in)    :: n
       real(8), intent(inout) :: ur_y(1:n), ui_y(1:n)
       real(8), intent(inout) :: ur_x(1:n), ui_x(1:n)
       real(8), intent(inout) :: u_t(1:2*n), u_s(1:2*n)

       integer :: nx, ny, ic, i, j, k
       integer, pointer :: req(:)
       integer :: reqr(2), reqs(2), x_snod, y_snod
       integer :: his_rank, her_rank

       real(8) :: d1,d2
       real(8) :: e1,e2
       real(8) :: e3,e4


          d1=eigen_get_wtime( )

          if ( x_nnod == 1 ) then
             if ( y_nnod == 1 ) then
                ur_y(1:n)=ur_x(1:n)
                ui_y(1:n)=ui_x(1:n)
             else
                ny=(n-1)/y_nnod+1
                do i=1,ny
                   j=y_inod+y_nnod*(i-1)
                   ur_y(i)=ur_x(j)
                   ui_y(i)=ui_x(j)
                end do! i
             end if
             return
          end if

          if ( x_nnod == y_nnod ) then
             if ( x_inod == y_inod ) then
                u_t(  1:  n) = ur_x(1:n)
                u_t(n+1:n+n) = ui_x(1:n)
             end if
             call bcast_dbl(u_t, 2*n, y_inod, x_COMM_WORLD)
             ur_y(1:n) = u_t(  1:  n)
             ui_y(1:n) = u_t(n+1:n+n)
             return
          end if

          x_snod = x_nnod/n_common
          y_snod = y_nnod/n_common

          if ( p0_(x_inod) > 0 ) then

             nx = (n-1)/y_snod+1
             do i=1,nx
                j = p0_(x_inod)+y_snod*(i-1)
                k = q0_(x_inod)+x_snod*(i-1)
                u_t(   i) = ur_x(j)
                u_t(nx+i) = ui_x(j)
                ur_y(k) = ur_x(j)
                ui_y(k) = ui_x(j)
             end do

             e3 = eigen_get_wtime( )
             call allgather_dbl( u_t(1), u_s(1), 2*nx, w_COMM_WORLD )
             e4 = eigen_get_wtime( )

             do ic=1,x_snod-1
                his_rank = MOD(x_inod-1 +x_nnod +ic*n_common,x_nnod)+1
                ny = (n-p0_(his_rank))/y_snod+1
                do i=1,ny
                   k = q0_(his_rank)+x_snod*(i-1)
                   j = 2*(his_rank-1)*nx+i
                   ur_y(k) = u_s(   j)
                   ui_y(k) = u_s(nx+j)
                end do! i
             end do! ic

             u_t(  1:  n) = ur_y(1:n)
             u_t(n+1:n+n) = ui_y(1:n)

             e1=eigen_get_wtime( )
             call bcast_dbl(u_t, 2*n, z_inod, z_COMM_WORLD)
             e2=eigen_get_wtime( )

          else

             i = MOD(y_inod-1,n_common)
             j = MOD(x_inod-1,n_common)
             ic = MOD(j-i+n_common,n_common)
             his_rank = ((x_inod-1)/n_common)*n_common + i + 1

             e1=eigen_get_wtime( )
             call bcast_dbl(u_s, 2*n, his_rank, z_COMM_WORLD)
             e2=eigen_get_wtime( )

             ur_y(1:n) = u_s(  1:  n)
             ui_y(1:n) = u_s(n+1:n+n)

          end if

9999      continue

          d2=eigen_get_wtime( )
          time_redist=time_redist+(d2-d1)-(e2-e1)


       return
       end subroutine  datacast_dbl2

