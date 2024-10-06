!--
       subroutine send_dbl(buf, n, idest, icom)
       implicit NONE

       integer, intent(in)    :: n, idest, icom
       real(8), intent(inout) :: buf(1:n)

       include 'mpif.h'
       integer                :: ierr


          call MPI_Send(buf, n, MPI_DOUBLE_PRECISION,
     &                  idest-1, 1, icom, ierr)


       return
       end subroutine ! send_dbl
!--
       subroutine send_dblt(buf, n, idest, itag, icom)
       implicit NONE

       integer, intent(in)    :: n, idest, itag, icom
       real(8), intent(inout) :: buf(1:n)

       include 'mpif.h'
       integer                :: ierr


          call MPI_Send(buf, n, MPI_DOUBLE_PRECISION,
     &                  idest-1, itag, icom, ierr)


       return
       end subroutine ! send_dblt
!--
       subroutine isend_dbl(buf, n, idest, ireq, icom)
       implicit NONE

       integer, intent(in)    :: n, idest, icom
       integer, intent(inout) :: ireq
       real(8), intent(inout) :: buf(1:n)

       include 'mpif.h'
       integer                :: ierr


          call MPI_Isend(buf,n,MPI_DOUBLE_PRECISION,
     &                   idest-1, 1, icom, ireq, ierr)


       return
       end subroutine ! isend_dbl
!--
       subroutine recv_dbl(buf, n, isrc, icom)
       implicit NONE

       integer, intent(in)    :: n, isrc, icom
       real(8), intent(inout) :: buf(1:n)

       include 'mpif.h'
       integer                :: ierr
       integer                :: status(MPI_STATUS_SIZE)


          call MPI_Recv(buf, n, MPI_DOUBLE_PRECISION,
     &                  isrc-1, 1, icom, status, ierr)


       return
       end subroutine ! recv_dbl
!--
       subroutine recv_dblt(buf, n, isrc, itag, icom)
       implicit NONE

       integer, intent(in)    :: n, isrc, itag, icom
       real(8), intent(inout) :: buf(1:n)

       include 'mpif.h'
       integer                :: ierr
       integer                :: status(MPI_STATUS_SIZE)


          call MPI_Recv(buf, n, MPI_DOUBLE_PRECISION,
     &                  isrc-1, itag, icom, status, ierr)


       return
       end subroutine ! recv_dblt
!--
       subroutine irecv_dbl(buf, n, isrc, ireq, icom)
       implicit NONE

       integer, intent(in)    :: n, isrc, icom
       integer, intent(inout) :: ireq
       real(8), intent(inout) :: buf(1:n)

       include 'mpif.h'
       integer                :: ierr


          call MPI_Irecv(buf, n, MPI_DOUBLE_PRECISION,
     &                   isrc-1, 1, icom, ireq, ierr)


       return
       end subroutine ! irecv_dbl
!--
       subroutine wait_dbl(ireq)
       implicit NONE

       integer, intent(inout) :: ireq

       include 'mpif.h'
       integer                :: status(MPI_STATUS_SIZE)
       integer                :: ierr


          call MPI_Wait(ireq, status, ierr)


       return
       end subroutine ! wait_dbl
!--
       subroutine waitall_dbl(n, ireq)
       implicit NONE

       integer, intent(in   ) :: n
       integer, intent(inout) :: ireq(n)

       include 'mpif.h'
       integer, pointer       :: status(:)
       integer                :: ierr


          allocate(status(MPI_STATUS_SIZE*n))
          call MPI_Waitall(n, ireq, status(1), ierr)
          deallocate(status)


       return
       end subroutine ! waitall_dbl
!--
       subroutine bcast_dbl(buf, n, iroot, icom)
       implicit NONE

       integer, intent(in)    :: n, iroot, icom
       real(8), intent(inout) :: buf(1:n)

       include 'mpif.h'
       integer                :: ierr

       integer :: world_size, my_rank, i, j

       real(8) :: d1,d2
       real(8) ::    time_bcast, time_reduce, time_redist, time_gather
       common /stat/ time_bcast, time_reduce, time_redist, time_gather


          d1=MPI_Wtime()

          call MPI_Bcast(buf, n, MPI_DOUBLE_PRECISION,
     &                   iroot-1, icom, ierr)

!          call d_bBcast(buf, n, iroot-1, icom, ierr)

          d2=MPI_Wtime()
          time_bcast=time_bcast+(d2-d1)


       return
       end subroutine ! bcast_dbl
!--
       subroutine reduce_dbl(buf, wrk, n, dist, icom)
       implicit NONE

       integer, intent(in)    :: n, dist, icom
       real(8), intent(inout) :: buf(1:n), wrk(1:n)

       include 'mpif.h'
       integer                :: ierr

       real(8) :: d1,d2
       real(8) ::    time_bcast, time_reduce, time_redist, time_gather
       common /stat/ time_bcast, time_reduce, time_redist, time_gather

          d1=MPI_Wtime()


          call MPI_Allreduce(buf, wrk, n, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, icom, ierr)
          buf(1:n) = wrk(1:n)
!          call MPI_Allreduce(MPI_IN_PLACE, buf, n, MPI_DOUBLE_PRECISION,
!     &                       MPI_SUM, icom, ierr)

          d2=MPI_Wtime()
          time_reduce=time_reduce+(d2-d1)

       return
       end subroutine ! reduce_dbl
!--
       subroutine allgather_dbl(buf, wrk, n, icom)
       implicit NONE

       integer, intent(in)    :: n, icom
       real(8), intent(inout) :: buf(1:n), wrk(1:n)

       include 'mpif.h'
       integer                :: ierr

       real(8) :: d1,d2
       real(8) ::    time_bcast, time_reduce, time_redist, time_gather
       common /stat/ time_bcast, time_reduce, time_redist, time_gather

          d1=MPI_Wtime()

          call MPI_Allgather(buf, n, MPI_DOUBLE_PRECISION,
     &                       wrk, n, MPI_DOUBLE_PRECISION,
     &                       icom, ierr)

          d2=MPI_Wtime()
          time_gather=time_gather+(d2-d1)

       return
       end subroutine ! allgather_dbl
!--
       integer function loop_c_start(istart, nnod, inod)
       implicit NONE

       integer, intent(in)    :: istart, nnod, inod


          loop_c_start = (istart+nnod-1-inod)/nnod+1


       return
       end function ! loop_c_start
!--
       integer function loop_c_end(iend, nnod, inod)
       implicit NONE

       integer, intent(in)    :: iend, nnod, inod


          loop_c_end   = (iend  +nnod  -inod)/nnod+0


       return
       end function ! loop_c_end
!--
       integer function loop_c_l2g_depth(ictr, nnod, inod)
       implicit NONE

       integer, intent(in)    :: ictr, nnod, inod


          loop_c_l2g_depth   = (ictr-1)*nnod+inod


       return
       end function ! loop_c_end
!--
       integer function loop_c_node(ictr, nnod, inod)
       implicit NONE

       integer, intent(in)    :: ictr, nnod, inod


          loop_c_node   = MOD(ictr-1,nnod)+1


       return
       end function ! loop_c_end
!--
       integer function loop_c_g2l_depth(ictr, nnod, inod)
       implicit NONE

       integer, intent(in)    :: ictr, nnod, inod


          loop_c_g2l_depth   = (ictr-1)/nnod+1


       return
       end function ! loop_c_end
!--
       subroutine loop_c(istart, iend, idoctr, isize, nnod, inod)
       implicit NONE

       integer, intent(inout) :: istart, iend, idoctr
       integer, intent(in)    :: isize, nnod, inod

       integer                :: istart0, iend0, idoctr0
       integer                :: istart1, iend1

       integer, external      :: loop_c_start, loop_c_end


          if ( istart > iend ) then
             istart=1; iend=0; idoctr=0
             return
          end if
          if ( nnod == 1 ) then
             idoctr=istart
             return
          end if

          istart1 = MAX(istart,1    )
          iend1   = MIN(iend  ,isize)

          istart = loop_c_start(istart1,nnod,inod)
          iend   = loop_c_end  (iend1,  nnod,inod)
          idoctr = (istart-1)*nnod+inod


       return
       end subroutine ! loop_c
!--

       subroutine datacast_init(ndim)
       implicit NONE

       include 'mpif.h'
       include 'trd.h'

       integer                ::  ndim
       integer                ::  n1, n2, n3, i, j, k, ierr
       integer :: old_grp, new_grp, ranks(0:4095)

       real(8) :: d1,d2
       real(8) ::    time_bcast, time_reduce, time_redist, time_gather
       common /stat/ time_bcast, time_reduce, time_redist, time_gather


       time_bcast = 0d0
       time_reduce= 0d0
       time_redist= 0d0
       time_gather= 0d0


          TRD_COMM_WORLD = MPI_COMM_EIGEN

          call MPI_Comm_size(TRD_COMM_WORLD, TRD_nnod, ierr)
          call MPI_Comm_rank(TRD_COMM_WORLD, TRD_inod, ierr)
          TRD_inod = TRD_inod+1

          if ( ndim == 1 ) then
             x_nnod = 1 !! fixed
          else
             x_nnod = INT(SQRT(DBLE(TRD_nnod)))
             i = 1
             if ( MOD(TRD_nnod,i) == 0 ) then
                k = i
             else
                k = 1
             end if
             do
                if ( x_nnod <= k ) exit
                if ( MOD(x_nnod,k) == 0 .AND.
     &               MOD(TRD_nnod,x_nnod) == 0 ) exit
                x_nnod = x_nnod-1
             end do!!
          endif
          y_nnod = TRD_nnod/x_nnod
!
!          print*,TRD_inod,"PE partition = ",x_nnod,y_nnod
!
!          x_inod =    (TRD_inod-1)/y_nnod +1
!          y_inod = MOD(TRD_inod-1, y_nnod)+1
          x_inod = MOD(TRD_inod-1, x_nnod)+1
          y_inod =    (TRD_inod-1)/x_nnod +1

          call MPI_Comm_split(TRD_COMM_WORLD,y_inod,x_inod,
     1         x_COMM_WORLD,ierr)
          call MPI_Comm_split(TRD_COMM_WORLD,x_inod,y_inod,
     1         y_COMM_WORLD,ierr)


          n1 = MAX(x_nnod,y_nnod)
          n2 = MIN(x_nnod,y_nnod)
          do
             if ( n1 == n2 ) then
                n_common = n1
                exit
             end if
             n3 = n1-n2
             n1 = MAX(n2,n3)
             n2 = MIN(n2,n3)
          end do!!

          p0_ = -1
          q0_ = -1
          do i=1,x_nnod
             if ( MOD(i-1,n_common) == MOD(y_inod-1,n_common) ) then
                n1 = y_inod-i
                if ( n1 >= 0 ) then
                   do j=1,x_nnod
                      k = +n1+(j-1)*y_nnod
                      if ( MOD(k,x_nnod) == 0 ) then
                         p0_(i) = k/x_nnod
                         q0_(i) = (j-1)
                         exit
                      end if
                   end do! j
                else
                   do j=1,y_nnod
                      k = -n1+(j-1)*x_nnod
                      if ( MOD(k,y_nnod) == 0 ) then
                         q0_(i) = k/y_nnod
                         p0_(i) = (j-1)
                         exit
                      end if
                   end do! j
                end if
             end if
          end do! i
          p0_ = p0_+1
          q0_ = q0_+1


          diag_0 = 0
          diag_1 = 0
          do i=1,y_nnod/n_common
             j = (i-1)*y_nnod+y_inod
             k = MOD(j-1,x_nnod)+1
             if ( k == x_inod ) then
                diag_0 = i
                diag_1 = (j-1)/x_nnod+1
                exit
             end if
          end do! i_1

!          call MPI_Barrier(TRD_COMM_WORLD, ierr)

       return
       end subroutine ! datacast_init
!--
       subroutine datacast_free(flag)
       implicit NONE

       include 'mpif.h'
       include 'trd.h'

       integer                :: flag, ierr

       real(8) :: d1,d2
       real(8) ::    time_bcast, time_reduce, time_redist, time_gather
       common /stat/ time_bcast, time_reduce, time_redist, time_gather


          call MPI_Comm_free(x_COMM_WORLD,ierr)
          call MPI_Comm_free(y_COMM_WORLD,ierr)

          TRD_COMM_WORLD=MPI_COMM_EIGEN

       if ( flag == 1 ) then
       if ( TRD_inod == 1 ) then
          print*, "COMM_STAT"
          print*, "   BCAST  :: ", time_bcast
          print*, "   REDUCE :: ", time_reduce
          print*, "   REDIST :: ", time_redist
          print*, "   GATHER :: ", time_gather
       endif
       endif


       return
       end subroutine ! datacast_free
!--
       subroutine datacast_dbl(u_y, u_x, u_t, u_s, n)
       implicit NONE

       integer, intent(in)    :: n
       real(8), intent(inout) :: u_y(1:n), u_x(1:n), u_t(1:n), u_s(1:n)

       include 'trd.h'
       include 'mpif.h'

       integer :: nx, ny, ic, i, j, k
       integer :: req(1024), reqr(2), reqs(2), x_snod, y_snod
       integer :: his_rank, her_rank

       real(8) :: d1,d2

       real(8) ::    time_bcast, time_reduce, time_redist, time_gather
       common /stat/ time_bcast, time_reduce, time_redist, time_gather


          d1=MPI_Wtime()

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

             nx = (n-p0_(x_inod))/y_snod+1
             do i=1,nx
                j = p0_(x_inod)+y_snod*(i-1)
                k = q0_(x_inod)+x_snod*(i-1)
                u_t(i) = u_x(j)
                u_y(k) = u_x(j)
!             print*,"I am",TRD_inod,"SdRv ",TRD_inod
             end do! i

             do ic=1,x_snod-1

!
! receiving message : length    = ifloor(n/y_snod)
!                   : #sender   = x_nnod-1
!
                his_rank = MOD(x_inod-1 +x_nnod +ic*n_common,x_nnod)+1
                ny = (n-p0_(his_rank))/y_snod+1
!             print*,"I am",TRD_inod,"Recv ",(his_rank-1)*y_nnod+y_inod
                call irecv_dbl(u_s, ny, his_rank, reqr(1), x_COMM_WORLD)

!
! sending message   : length    = ifloor(n/y_snod)
!                   : #receiver = x_nnod-1
!
                her_rank = MOD(x_inod-1 +x_nnod -ic*n_common,x_nnod)+1
                call isend_dbl(u_t, nx, her_rank, req(ic), x_COMM_WORLD)
!             print*,"I am",TRD_inod,"Send ",(her_rank-1)*y_nnod+y_inod

                call wait_dbl(reqr(1))

                do i=1,ny
                   k = q0_(his_rank)+x_snod*(i-1)
                   u_y(k) = u_s(i)
                end do! i

             end do! ic
                call waitall_dbl(x_snod-1, req)

             do ic=1,n_common-1
                her_rank = MOD(x_inod-1 +x_nnod +ic,x_nnod)+1
!             print*,"I am",TRD_inod,"Send ",(her_rank-1)*y_nnod+y_inod
                call isend_dbl(u_y, n, her_rank, req(ic), x_COMM_WORLD)
!                call wait_dbl(req(ic))
             end do! ic
                call waitall_dbl(n_common-1, req)

          else

                i = MOD(y_inod-1,n_common)
                j = MOD(x_inod-1,n_common)
                ic = MOD(j-i+n_common,n_common)
                his_rank = MOD(x_inod-1 +x_nnod -ic,x_nnod)+1
!             print*,"I am",TRD_inod,"Recv ",(his_rank-1)*y_nnod+y_inod
                call recv_dbl(u_y, n, his_rank, x_COMM_WORLD)

          end if

9999      continue

          d2=MPI_Wtime()
          time_redist=time_redist+(d2-d1)


       return
       end subroutine ! datacast_dbl
!--
       subroutine datacast_dbl2(ur_y,ui_y, ur_x,ui_x, u_t,u_s, n)
       implicit NONE

       integer, intent(in)    :: n
       real(8), intent(inout) :: ur_y(1:n), ui_y(1:n)
       real(8), intent(inout) :: ur_x(1:n), ui_x(1:n)
       real(8), intent(inout) :: u_t(1:2*n), u_s(1:2*n)

       include 'trd.h'
       include 'mpif.h'

       integer :: nx, ny, ic, i, j, k
       integer :: req(1024), reqr(2), reqs(2), x_snod, y_snod
       integer :: his_rank, her_rank

       real(8) :: d1,d2

       real(8) ::    time_bcast, time_reduce, time_redist, time_gather
       common /stat/ time_bcast, time_reduce, time_redist, time_gather


          d1=MPI_Wtime()

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

             nx = (n-p0_(x_inod))/y_snod+1
             do i=1,nx
                j = p0_(x_inod)+y_snod*(i-1)
                k = q0_(x_inod)+x_snod*(i-1)
                u_t(   i) = ur_x(j)
                u_t(nx+i) = ui_x(j)
                ur_y(k) = ur_x(j)
                ui_y(k) = ui_x(j)
!             print*,"I am",TRD_inod,"SdRv ",TRD_inod
             end do! i

             do ic=1,x_snod-1

!
! receiving message : length    = ifloor(n/y_snod)
!                   : #sender   = x_nnod-1
!
                his_rank = MOD(x_inod-1 +x_nnod +ic*n_common,x_nnod)+1
                ny = (n-p0_(his_rank))/y_snod+1
!             print*,"I am",TRD_inod,"Recv ",(his_rank-1)*y_nnod+y_inod
                call irecv_dbl(u_s, 2*ny,
     &                        his_rank, reqr(1), x_COMM_WORLD)

!
! sending message   : length    = ifloor(n/y_snod)
!                   : #receiver = x_nnod-1
!
                her_rank = MOD(x_inod-1 +x_nnod -ic*n_common,x_nnod)+1
                call isend_dbl(u_t, 2*nx,
     &                         her_rank, req(ic), x_COMM_WORLD)
!             print*,"I am",TRD_inod,"Send ",(her_rank-1)*y_nnod+y_inod

                call wait_dbl(reqr(1))

                do i=1,ny
                   k = q0_(his_rank)+x_snod*(i-1)
                   ur_y(k) = u_s(   i)
                   ui_y(k) = u_s(ny+i)
                end do! i

             end do! ic
                call waitall_dbl(x_snod-1, req)

             do ic=1,n_common-1
                her_rank = MOD(x_inod-1 +x_nnod +ic,x_nnod)+1
!             print*,"I am",TRD_inod,"Send ",(her_rank-1)*y_nnod+y_inod
                u_t(  1:  n) = ur_y(1:n)
                u_t(n+1:n+n) = ui_y(1:n)
!                call send_dbl(u_t, 2*n, her_rank, x_COMM_WORLD)
                call isend_dbl(u_t, 2*n,
     &                         her_rank, req(ic), x_COMM_WORLD)
             end do! ic
                call waitall_dbl(n_common-1, req)

          else

                i = MOD(y_inod-1,n_common)
                j = MOD(x_inod-1,n_common)
                ic = MOD(j-i+n_common,n_common)
                his_rank = MOD(x_inod-1 +x_nnod -ic,x_nnod)+1
!             print*,"I am",TRD_inod,"Recv ",(his_rank-1)*y_nnod+y_inod
                call recv_dbl(u_s, 2*n,
     &                        his_rank, x_COMM_WORLD)
                ur_y(1:n) = u_s(  1:  n)
                ui_y(1:n) = u_s(n+1:n+n)

          end if

9999      continue

          d2=MPI_Wtime()
          time_redist=time_redist+(d2-d1)


       return
       end subroutine ! datacast_dbl2

      subroutine CYC1D_CYC2D(n,a,lda,b,ldb, x_nnod, y_nnod, iret)
      implicit NONE
*
      include 'mpif.h'
      include 'commtxt.h'
*
      integer                 :: n,lda,ldb,iret
      real(8)                 :: a(lda,1), b(ldb,1)
      integer                 :: x_nnod, y_nnod
*
      real(8), pointer        :: work(:,:)
      integer                 :: world_size, my_rank, sq_size
      integer                 :: his_rank, her_rank
      integer                 :: his_posx, her_posx
      integer                 :: his_posy, her_posy
      integer                 :: my_posx, my_posy, my_posu, my_posv
      integer                 :: size_s,size_r
      integer                 :: nws, nqs
      integer                 :: i,j,k, ic,jc,kc,lc,id
*
      integer                 :: status(MPI_STATUS_SIZE)
      integer                 :: type
      integer                 :: comm
      integer                 :: tag
      integer                 :: req
      integer                 :: error
*
      iret = 0
*
      if ( n > lda ) then
         iret = 1
         return
      end if
*
*    [*,world_size] --> [y_nnod,x_nnod]
*     A(lda,nws) --> B(ldb,nqs)
*
      comm = MPI_COMM_EIGEN
      type = MPI_DOUBLE_PRECISION
      tag  = 10000
*
      call MPI_Comm_size(comm,world_size,error)
      call MPI_Comm_rank(comm,my_rank,error)
      my_rank = my_rank+1
!      sq_size=int(sqrt(dble(world_size)))
*
!      y_nnod = sq_size
!      x_nnod = world_size / y_nnod
      if ( x_nnod*y_nnod /= world_size ) then
         iret = 3
         return
      end if
*
      nws = (n-1)/world_size+1
      nqs = (n-1)/x_nnod+1
      if ( nqs > ldb ) then
         iret = 2
         return
      end if
*
!      if ( x_nnod > 1 ) then
         allocate(work(ldb,x_nnod), stat=error)
         if ( error /= 0 ) then
            iret = -1
            return
         end if
!      end if
*
!      my_posy = MOD(my_rank-1, y_nnod) + 1
!      my_posx =    (my_rank-1)/y_nnod  + 1
      my_posx = MOD(my_rank-1, x_nnod) + 1
      my_posy =    (my_rank-1)/x_nnod  + 1
      my_posu = MOD(my_rank-1, y_nnod) + 1
      my_posv =    (my_rank-1)/y_nnod  + 1
*
*
      call MPI_Barrier(comm,error)
*
      do ic=1,nws
*
         kc = my_rank+(ic-1)*world_size
         if ( kc <= n ) then
*
            do id=1,x_nnod
               his_posx = MOD(my_posv+id-2,x_nnod)+1
               size_s   =    (n-his_posx)/x_nnod+1
               do i=1,size_s
                  j = his_posx+(i-1)*x_nnod
                  work(i,id) = A(j,ic)
               end do
            end do
*
         end if
*
* BUCKET Algorithm
*
         do id=1,x_nnod
*-
            his_posx = MOD(my_posv+id-2,x_nnod)+1
            his_rank =    (his_posx-1)*1+(my_posu-1)*x_nnod+1

            if ( kc <= n ) then
               size_s = (n-his_posx)/x_nnod+1
            else
               size_s = 0
            end if

            if ( size_s > 0 ) then
               call MPI_Isend(work(1,id),size_s,type,
     &                     his_rank-1,tag,comm,req,error)
!             print*,"I am",my_rank,"Send ",his_rank,work(1:size_s,id)
            end if
*-
            her_posx = MOD(my_posx-1+x_nnod-id+1,x_nnod)+1
            her_rank =    (her_posx-1)*y_nnod+(my_posy-1)*1+1

            lc = her_rank+(ic-1)*world_size

            if ( lc <= n ) then
               size_r = (n-my_posx)/x_nnod+1
            else
               size_r = 0
            end if

            if ( size_r > 0 ) then
               jc = her_posx+(ic-1)*x_nnod
               call MPI_Recv(B(1,jc),size_r,type,
     &                     her_rank-1,tag,comm,status,error)
!             print*,"I am",my_rank,"Recv ",her_rank,B(1:size_r,jc),jc
            end if
*-
            if ( size_s > 0 ) then
               call MPI_Wait(req,status,error)
            end if
*-
         end do

      end do
*
      call MPI_Barrier(comm,error)
*
!      if ( x_nnod > 1 ) then
         deallocate(work)
!      end if
*
      return
      end
      subroutine CYC2D_CYC1D(n,a,lda,b,ldb, x_nnod, y_nnod, iret)
      implicit NONE
*
      include 'mpif.h'
      include 'commtxt.h'
*
      integer                 :: n,lda,ldb,iret
      real(8)                 :: a(lda,1), b(ldb,1)
      integer                 :: x_nnod, y_nnod
*
      real(8), pointer        :: work(:,:)
      integer                 :: world_size, my_rank, sq_size
      integer                 :: his_rank, her_rank
      integer                 :: his_posx, her_posx
      integer                 :: his_posy, her_posy
      integer                 :: my_posx, my_posy, my_posu, my_posv
      integer                 :: size_s,size_r
      integer                 :: nws, nqs
      integer                 :: i,j,k, ic,jc,kc,lc,id
*
      integer                 :: status(MPI_STATUS_SIZE)
      integer                 :: type
      integer                 :: comm
      integer                 :: tag
      integer                 :: req
      integer                 :: error
*
      iret = 0
*
      if ( n > lda ) then
         iret = 1
         return
      end if
*
*    [*,world_size] <-- [y_nnod,x_nnod]
*     A(lda,nws) <-- B(ldb,nqs)
*
      comm = MPI_COMM_EIGEN
      type = MPI_DOUBLE_PRECISION
      tag  = 10000
*
      call MPI_Comm_size(comm,world_size,error)
      call MPI_Comm_rank(comm,my_rank,error)
      my_rank = my_rank+1
!      sq_size = int(sqrt(dble(world_size)))
*
!      y_nnod = sq_size
!      x_nnod = world_size / y_nnod
      if ( x_nnod*y_nnod /= world_size ) then
         iret = 3
         return
      end if
*
      nws = (n-1)/world_size+1
      nqs = (n-1)/x_nnod+1
      if ( nqs > ldb ) then
         iret = 2
         return
      end if
*
!      if ( x_nnod > 1 ) then
         allocate(work(ldb,x_nnod), stat=error)
         if ( error /= 0 ) then
            iret = -1
            return
         end if
!      end if
*
!      my_posy = MOD(my_rank-1, y_nnod) + 1
!      my_posx =    (my_rank-1)/y_nnod  + 1
      my_posx = MOD(my_rank-1, x_nnod) + 1
      my_posy =    (my_rank-1)/x_nnod  + 1
      my_posu = MOD(my_rank-1, y_nnod) + 1
      my_posv =    (my_rank-1)/y_nnod  + 1
*
*
      call MPI_Barrier(comm,error)
*
      do ic=1,nws
*
         kc = my_rank+(ic-1)*world_size
*
* BUCKET Algorithm
*
         do id=1,x_nnod
*-
            her_posx = MOD(my_posx-1+x_nnod-id+1,x_nnod)+1
            her_rank =    (her_posx-1)*y_nnod+(my_posy-1)*1+1

            lc = her_rank+(ic-1)*world_size
            if ( lc <= n ) then
               size_r = (n-my_posx)/x_nnod+1
            else
               size_r = 0
            end if

            if ( size_r > 0 ) then
               jc = her_posx+(ic-1)*x_nnod
               call MPI_Isend(B(1,jc),size_r,type,
     &                     her_rank-1,tag,comm,req,error)
!              print*,"I am",my_rank,"Recv ",her_rank,B(1:size_r,jc)
            end if
*-
            his_posx = MOD(my_posv+id-2,x_nnod)+1
            his_rank =    (his_posx-1)*1+(my_posu-1)*x_nnod+1

            if ( kc <= n ) then
               size_s = (n-his_posx)/x_nnod+1
            else
               size_s = 0
            end if

            if ( size_s > 0 ) then
               call MPI_Recv(work(1,id),size_s,type,
     &                     his_rank-1,tag,comm,status,error)
!              print*,"I am",my_rank,"Send ",his_rank,work(1:size_s,id)
            end if
*-
            if ( size_r>0)then
               call MPI_Wait(req,status,error)
            end if
*-
         end do

         if ( kc <= n ) then
            do id=1,x_nnod
               his_posx = MOD(my_posv+id-2,x_nnod)+1
               size_s   =    (n-his_posx)/x_nnod+1
               do i=1,size_s
                  j = his_posx+(i-1)*x_nnod
                  A(j,ic) = work(i,id)
               end do
            end do
         end if
*
      end do
*
      call MPI_Barrier(comm,error)
*
!      if ( x_nnod > 1 ) then
         deallocate(work)
!      end if
*
      return
      end
      subroutine pcast_cyclic_real(u, a, nma, i_st, i_ed, comm)

      implicit NONE
      integer, intent(in)     :: nma, i_st, i_ed, comm
      real(8), intent(inout)  :: u(*),a(nma,*)

      include'mpif.h'

      real(8), pointer        :: buff(:)

      integer                 :: i_inod, i_nnod
      integer                 :: i_status(MPI_STATUS_SIZE), i_req
      integer, parameter      :: i_dtype = MPI_DOUBLE_PRECISION
      integer                 :: i_stat, i_err

      integer                 :: i_1, j_1
      integer                 :: i1, i2
      integer                 :: j1, j2, j3, j4, j5, j6
      integer                 :: im0, im1, im2, im3
      integer                 :: imid0, imid2
      integer                 :: iwidth, idisp, ilength
      integer                 :: idispr, idisps, icross
      integer                 :: ilhs, irhs


      if ( i_st > i_ed ) return

      call MPI_Comm_size(comm, i_nnod, i_err)
      call MPI_Comm_rank(comm, i_inod, i_err)
      i_inod  = i_inod+1

      if ( i_nnod == 1 ) then
         do i_1=i_st,i_ed
            u(i_1) = a(1,i_1)
         end do! i_1
         return
      end if

      i1 =    ((i_st-1)/i_nnod)+1
      j1 = MOD((i_st-1),i_nnod)+1
      i2 =    ((i_ed-1)/i_nnod)+1
      j2 = MOD((i_ed-1),i_nnod)+1
      iwidth = i2-i1+1

      imid2 = 2**(INT(LOG(DBLE(2*i_nnod-1))/LOG(2.00D+00)))
      imid0 = imid2/2

      allocate(buff(iwidth*imid2), stat=i_stat)
      if ( i_stat /= 0 ) then
         print*,"memory allocation error pcast"
         call mpi_abort(comm,1,i_err)
      end if

      idisp = (i_inod-1)*2*iwidth
      if ( i_inod > imid0 ) idisp = 0
      do i_1=1,iwidth
         buff(i_1+idisp) = a(1,i_1+i1-1)
      end do! i_1

      if ( i_inod <= imid0 ) then
         if ( i_inod+imid0 <= i_nnod ) then
            call mpi_recv(
     1          buff(1+idisp+iwidth),
     2          iwidth, i_dtype, (i_inod+imid0)-1,
     3          90, comm, i_status, i_err)
         end if

         im0 = 1
         do
            if ( im0 >= imid0 ) exit

            im2 = im0
            im0 = im0*2
            im1 = MOD((i_inod-1),im0)+1
            im3 =    ((i_inod-1)/im0)+1
            ilength = (iwidth*im0)

            if ( im1 <= im2 ) then
               idisps = (im3-1)*ilength*2
               idispr = idisps +ilength
               icross = i_inod+im2-1
            else
               idispr = (im3-1)*ilength*2
               idisps = idispr +ilength
               icross = i_inod-im2-1
            end if

            call mpi_isend(
     1          buff(1+idisps),
     2          ilength, i_dtype, icross,
     3          100+im0, comm, i_req, i_err)
            call mpi_recv(
     1          buff(1+idispr),
     2          ilength, i_dtype, icross,
     3          100+im0, comm, i_status, i_err)
            call mpi_wait(i_req, i_status, i_err)

         end do!!

         if ( i_inod+imid0 <= i_nnod ) then
            ilength = iwidth*imid2
            call mpi_send(
     1          buff(1),
     2          ilength, i_dtype, (i_inod+imid0)-1,
     3          95, comm, i_err)
         endif
      else
            call mpi_send(
     1          buff(1),
     2          iwidth, i_dtype, (i_inod-imid0)-1,
     3          90, comm, i_err)
            ilength = iwidth*imid2
            call mpi_recv(
     1          buff(1),
     2          ilength, i_dtype, (i_inod-imid0)-1,
     3          95, comm, i_status, i_err)
      endif

      if ( iwidth /= 1 ) then

         do j_1=1,imid0
         do i_1=2,iwidth-1
            ilhs = (i_1-1)*i_nnod+j_1
            irhs = i_1+(j_1-1)*(2*iwidth)
            u(ilhs) = buff(irhs)
         end do! i_1
         end do! j_1
         i_1=1
         do j_1=j1,imid0
            ilhs = (i_1-1)*i_nnod+j_1
            irhs = i_1+(j_1-1)*(2*iwidth)
            u(ilhs) = buff(irhs)
         end do! j_1
         i_1=iwidth
         do j_1=1,MIN(imid0,j2)
            ilhs = (i_1-1)*i_nnod+j_1
            irhs = i_1+(j_1-1)*(2*iwidth)
            u(ilhs) = buff(irhs)
         end do! j_1

         j4 = i_nnod-imid0
         j5 = MAX(1,j1-imid0)
         j6 = MIN(j4,j2-imid0)

         do j_1=1,j4
         do i_1=2,iwidth-1
            ilhs = (i_1-1)*i_nnod+j_1+imid0
            irhs = i_1+(j_1-1)*(2*iwidth)+iwidth
            u(ilhs) = buff(irhs)
         end do! i_1
         end do! j_1
         i_1=1
         do j_1=j5,j4
            ilhs = (i_1-1)*i_nnod+j_1+imid0
            irhs = i_1+(j_1-1)*(2*iwidth)+iwidth
            u(ilhs) = buff(irhs)
         end do! j_1
         i_1=iwidth
         do j_1=1,j6
            ilhs = (i_1-1)*i_nnod+j_1+imid0
            irhs = i_1+(j_1-1)*(2*iwidth)+iwidth
            u(ilhs) = buff(irhs)
         end do! j_1

      else

         do i_1=j1,MIN(j2,imid0)
            u(i_1      ) = buff(2*i_1-1)
         end do! i_1
         do i_1=1,j2-imid0
            u(i_1+imid0) = buff(2*i_1  )
         end do! i_1

      end if

      deallocate(buff)

      return
      end subroutine ! pcast_cyclic_real

      subroutine pcast_cyclic_int(u, a, nma, i_st, i_ed, comm)

      implicit NONE
      integer, intent(in)     :: nma, i_st, i_ed, comm
      integer, intent(inout)  :: u(*),a(nma,*)

      include'mpif.h'

      integer, pointer        :: buff(:)

      integer                 :: i_inod, i_nnod
      integer                 :: i_status(MPI_STATUS_SIZE), i_req
      integer, parameter      :: i_dtype = MPI_INTEGER
      integer                 :: i_stat, i_err

      integer                 :: i_1, j_1
      integer                 :: i1, i2
      integer                 :: j1, j2, j3, j4, j5, j6
      integer                 :: im0, im1, im2, im3
      integer                 :: imid0, imid2
      integer                 :: iwidth, idisp, ilength
      integer                 :: idispr, idisps, icross
      integer                 :: ilhs, irhs


      if ( i_st > i_ed ) return

      call MPI_Comm_size(comm, i_nnod, i_err)
      call MPI_Comm_rank(comm, i_inod, i_err)
      i_inod  = i_inod+1

      if ( i_nnod == 1 ) then
         do i_1=i_st,i_ed
            u(i_1) = a(1,i_1)
         end do! i_1
         return
      end if

      i1 =    ((i_st-1)/i_nnod)+1
      j1 = MOD((i_st-1),i_nnod)+1
      i2 =    ((i_ed-1)/i_nnod)+1
      j2 = MOD((i_ed-1),i_nnod)+1
      iwidth = i2-i1+1

      imid2 = 2**(INT(LOG(DBLE(2*i_nnod-1))/LOG(2.00D+00)))
      imid0 = imid2/2

      allocate(buff(iwidth*imid2), stat=i_stat)
      if ( i_stat /= 0 ) then
         print*,"memory allocation error pcast"
         call mpi_abort(comm,1,i_err)
      end if

      idisp = (i_inod-1)*2*iwidth
      if ( i_inod > imid0 ) idisp = 0
      do i_1=1,iwidth
         buff(i_1+idisp) = a(1,i_1+i1-1)
      end do! i_1

      if ( i_inod <= imid0 ) then
         if ( i_inod+imid0 <= i_nnod ) then
            call mpi_recv(
     1          buff(1+idisp+iwidth),
     2          iwidth, i_dtype, (i_inod+imid0)-1,
     3          90, comm, i_status, i_err)
         end if

         im0 = 1
         do
            if ( im0 >= imid0 ) exit

            im2 = im0
            im0 = im0*2
            im1 = MOD((i_inod-1),im0)+1
            im3 =    ((i_inod-1)/im0)+1
            ilength = (iwidth*im0)

            if ( im1 <= im2 ) then
               idisps = (im3-1)*ilength*2
               idispr = idisps +ilength
               icross = i_inod+im2-1
            else
               idispr = (im3-1)*ilength*2
               idisps = idispr +ilength
               icross = i_inod-im2-1
            end if

            call mpi_isend(
     1          buff(1+idisps),
     2          ilength, i_dtype, icross,
     3          100+im0, comm, i_req, i_err)
            call mpi_recv(
     1          buff(1+idispr),
     2          ilength, i_dtype, icross,
     3          100+im0, comm, i_status, i_err)
            call mpi_wait(i_req, i_status, i_err)

         end do!!

         if ( i_inod+imid0 <= i_nnod ) then
            ilength = iwidth*imid2
            call mpi_send(
     1          buff(1),
     2          ilength, i_dtype, (i_inod+imid0)-1,
     3          95, comm, i_err)
         endif
      else
            call mpi_send(
     1          buff(1),
     2          iwidth, i_dtype, (i_inod-imid0)-1,
     3          90, comm, i_err)
            ilength = iwidth*imid2
            call mpi_recv(
     1          buff(1),
     2          ilength, i_dtype, (i_inod-imid0)-1,
     3          95, comm, i_status, i_err)
      endif

      if ( iwidth /= 1 ) then

         do j_1=1,imid0
         do i_1=2,iwidth-1
            ilhs = (i_1-1)*i_nnod+j_1
            irhs = i_1+(j_1-1)*(2*iwidth)
            u(ilhs) = buff(irhs)
         end do! i_1
         end do! j_1
         i_1=1
         do j_1=j1,imid0
            ilhs = (i_1-1)*i_nnod+j_1
            irhs = i_1+(j_1-1)*(2*iwidth)
            u(ilhs) = buff(irhs)
         end do! j_1
         i_1=iwidth
         do j_1=1,MIN(imid0,j2)
            ilhs = (i_1-1)*i_nnod+j_1
            irhs = i_1+(j_1-1)*(2*iwidth)
            u(ilhs) = buff(irhs)
         end do! j_1

         j4 = i_nnod-imid0
         j5 = MAX(1,j1-imid0)
         j6 = MIN(j4,j2-imid0)

         do j_1=1,j4
         do i_1=2,iwidth-1
            ilhs = (i_1-1)*i_nnod+j_1+imid0
            irhs = i_1+(j_1-1)*(2*iwidth)+iwidth
            u(ilhs) = buff(irhs)
         end do! i_1
         end do! j_1
         i_1=1
         do j_1=j5,j4
            ilhs = (i_1-1)*i_nnod+j_1+imid0
            irhs = i_1+(j_1-1)*(2*iwidth)+iwidth
            u(ilhs) = buff(irhs)
         end do! j_1
         i_1=iwidth
         do j_1=1,j6
            ilhs = (i_1-1)*i_nnod+j_1+imid0
            irhs = i_1+(j_1-1)*(2*iwidth)+iwidth
            u(ilhs) = buff(irhs)
         end do! j_1

      else

         do i_1=j1,MIN(j2,imid0)
            u(i_1      ) = buff(2*i_1-1)
         end do! i_1
         do i_1=1,j2-imid0
            u(i_1+imid0) = buff(2*i_1  )
         end do! i_1

      end if

      deallocate(buff)

      return
      end subroutine ! pcast_cyclic_int

