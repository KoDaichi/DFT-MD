      module communication_h
      implicit none

      include 'mpif.h'

      public  :: send_dbl
      private :: isend_dbl
      public  :: recv_dbl
      private :: irecv_dbl
      private :: wait_dbl
      public  :: bcast_dbl
      public  :: reduce_dbl
      public  :: get_loop_start
      public  :: get_loop_end
      public  :: translate_l2g
      public  :: translate_g2l
      public  :: get_owner_node
      public  :: loop_c
      public  :: eigen_init
      public  :: eigen_free
      public  :: datacast_dbl
      public  :: cyc1d_cyc2d
      public  :: cyc2d_cyc1d
      private :: ierr

      integer :: ierr

      contains
!----
      subroutine send_dbl(buf, n, idest, icom)
      implicit none

      integer, intent(in)    :: n, idest, icom
      real(8), intent(inout) :: buf(1:n)

      call mpi_send(buf, n, mpi_double_precision,
     &              idest-1, 1, icom, ierr)

      return
      end subroutine ! send_dbl
!----
      subroutine isend_dbl(buf, n, idest, ireq, icom)
      implicit none

      integer, intent(in)    :: n, idest, icom
      integer, intent(inout) :: ireq
      real(8), intent(inout) :: buf(1:n)

      call mpi_isend(buf,n,mpi_double_precision,
     &               idest-1, 1, icom, ireq, ierr)

      return
      end subroutine ! isend_dbl
!----
      subroutine recv_dbl(buf, n, isrc, icom)
      implicit none

      integer, intent(in)    :: n, isrc, icom
      real(8), intent(inout) :: buf(1:n)
      integer                :: status(mpi_status_size)

      call mpi_recv(buf, n, mpi_double_precision,
     &              isrc-1, 1, icom, status, ierr)

      return
      end subroutine ! recv_dbl
!----
      subroutine irecv_dbl(buf, n, isrc, ireq, icom)
      implicit none

      integer, intent(in)    :: n, isrc, icom
      integer, intent(inout) :: ireq
      real(8), intent(inout) :: buf(1:n)

      call mpi_irecv(buf, n, mpi_double_precision,
     &               isrc-1, 1, icom, ireq, ierr)

      return
      end subroutine ! irecv_dbl
!----
      subroutine wait_dbl(ireq)
      implicit none

      integer, intent(inout) :: ireq

      integer                :: status(mpi_status_size)

      call mpi_wait(ireq, status, ierr)

      return
      end subroutine ! wait_dbl
!----
      subroutine bcast_dbl(buf, n, iroot, icom)
      implicit none

      integer, intent(in)    :: n, iroot, icom
      real(8), intent(inout) :: buf(1:n)

      real(8) ::    d1, d2
      real(8) ::    time_bcast, time_reduce, time_redist
      common /stat/ time_bcast, time_reduce, time_redist

#ifdef TIMER
      call mpi_barrier(icom,ierr)
      d1=mpi_wtime()
#endif

      call mpi_bcast(buf, n, mpi_double_precision,
     &               iroot-1, icom, ierr)

#ifdef TIMER
      d2=mpi_wtime()
      time_bcast=time_bcast+(d2-d1)
#endif

      return
      end subroutine ! bcast_dbl
!----
      subroutine reduce_dbl(buf, wrk, n, icom)
      implicit none

      integer, intent(in)    :: n, icom
      real(8), intent(inout) :: buf(1:n), wrk(1:n)

      real(8) ::    d1, d2
      real(8) ::    time_bcast, time_reduce, time_redist
      common /stat/ time_bcast, time_reduce, time_redist

#ifdef TIMER
      call mpi_barrier(icom,ierr)
      d1=mpi_wtime()
#endif

      call mpi_allreduce(buf, wrk, n, mpi_double_precision,
     &                   mpi_sum, icom, ierr)
      buf(1:n) = wrk(1:n)

#ifdef TIMER
      d2=mpi_wtime()
      time_reduce=time_reduce+(d2-d1)
#endif

      return
      end subroutine ! reduce_dbl
!----
      integer function get_loop_start(istart, nnod, inod)
      implicit none

      integer, intent(in)    :: istart, nnod, inod

      get_loop_start = (istart+nnod-1-inod)/nnod+1

      return
      end function ! get_loop_start
!----
      integer function get_loop_end(iend, nnod, inod)
      implicit none

      integer, intent(in)    :: iend, nnod, inod

      get_loop_end   = (iend  +nnod  -inod)/nnod+0

      return
      end function ! get_loop_end
!----
      integer function translate_l2g(ictr, nnod, inod)
      implicit none

      integer, intent(in)    :: ictr, nnod, inod

      translate_l2g   = (ictr-1)*nnod+inod

      return
      end function ! translate_l2g
!----
      integer function get_owner_node(ictr, nnod, inod)
      implicit none

      integer, intent(in)    :: ictr, nnod, inod

      get_owner_node   = mod(ictr-1,nnod)+1

      return
      end function ! get_owner_node
!----
      integer function translate_g2l(ictr, nnod, inod)
      implicit none

      integer, intent(in)    :: ictr, nnod, inod

      translate_g2l   = (ictr-1)/nnod+1

      return
      end function ! translate_g2l
!----
      subroutine loop_c(istart, iend, idoctr, isize, nnod, inod)
      implicit none

      integer, intent(inout) :: istart, iend, idoctr
      integer, intent(in)    :: isize, nnod, inod

      integer                :: istart0, iend0, idoctr0
      integer                :: istart1, iend1

!---- integer, external      :: get_loop_start, get_loop_end

      if ( istart > iend ) then
         istart=1; iend=0; idoctr=0
         return
      end if
      if ( nnod == 1 ) then
         idoctr=istart
         return
      end if

      istart1 = max(istart,1    )
      iend1   = min(iend  ,isize)

      istart = get_loop_start(istart1,nnod,inod)
      iend   = get_loop_end  (iend1,  nnod,inod)
      idoctr = (istart-1)*nnod+inod

      return
      end subroutine ! loop_c
!----

      subroutine eigen_init(ndim)
      implicit none

      include 'trd.h'

      integer                ::  ndim
      integer                ::  n1, n2, n3, i, j, k, ierr
      integer                ::  base_comm, base_inod, base_nnod

!---- integer, parameter     ::  cluster_const = 2 ! test case
      integer, parameter     ::  cluster_const = 8 ! es

      call mpi_comm_size(mpi_comm_eigen, nprocs, ierr)
      call mpi_comm_rank(mpi_comm_eigen, myrank, ierr)
      myrank = myrank+1

      if ( ndim == 1 ) then
         size_of_col = 1 !! fixed
      else if ( nprocs < cluster_const**2 ) then
         size_of_col = int(sqrt(dble(nprocs)))
         do
            if ( size_of_col == 1 ) exit
            if ( mod(nprocs,size_of_col) == 0 ) exit
            size_of_col = size_of_col-1
         end do!!
         else
!----
            size_of_col = int(sqrt(dble(nprocs)))
         if ( mod(nprocs,cluster_const) == 0 ) then
            do
               if ( size_of_col <= cluster_const ) exit
               if ( mod(size_of_col,cluster_const) == 0 .and.
     &              mod(nprocs,size_of_col) == 0 ) exit
                  size_of_col = size_of_col-1
            end do!!
         else
            size_of_col = 1
         end if
!----
      endif
      size_of_row = nprocs/size_of_col
!----
!---- if(myrank==1)print*,"pe partition = ",size_of_col,size_of_row
!----
      my_col = mod(myrank-1, size_of_col)+1
      my_row =    (myrank-1)/size_of_col +1

      call mpi_comm_split(mpi_comm_eigen,my_row,my_col,
     &                    mpi_comm_col,ierr)
      call mpi_comm_split(mpi_comm_eigen,my_col,my_row,
     &                    mpi_comm_row,ierr)

!----
      if ( mod(size_of_col,cluster_const) == 0 ) then
         base_nnod = size_of_col
         base_inod = my_col
         base_comm = mpi_comm_col
      else
         base_nnod = nprocs
         base_inod = myrank
         base_comm = mpi_comm_eigen
      endif

      es_cluster = cluster_const
      es_icnod = mod(base_inod-1, es_cluster)+1
      es_ignod =    (base_inod-1)/es_cluster +1
      call mpi_comm_split(base_comm,es_ignod,es_icnod,
     &                    es_ccomm_world,ierr)
      call mpi_comm_split(base_comm,es_icnod,es_ignod,
     &                    es_gcomm_world,ierr)
      call mpi_comm_size(es_ccomm_world, es_ncnod, ierr)
      call mpi_comm_size(es_gcomm_world, es_ngnod, ierr)
!----

      n1 = max(size_of_col,size_of_row)
      n2 = min(size_of_col,size_of_row)
      do
         if ( n1 == n2 ) then
            n_common = n1
            exit
         end if
         n3 = n1-n2
         n1 = max(n2,n3)
         n2 = min(n2,n3)
      end do!!


      p0_ = -1
      q0_ = -1
      do i=1,size_of_col
         if ( mod(i-1,n_common) == mod(my_row-1,n_common) ) then
            n1 = my_row-i
            if ( n1 >= 0 ) then
               do j=1,size_of_col
                  k = +n1+(j-1)*size_of_row
                  if ( mod(k,size_of_col) == 0 ) then
                     p0_(i) = k/size_of_col
                     q0_(i) = (j-1)
                     exit
                  end if
               end do! j
            else
               do j=1,size_of_row
                  k = -n1+(j-1)*size_of_col
                  if ( mod(k,size_of_row) == 0 ) then
                     q0_(i) = k/size_of_row
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
      do i=1,size_of_row/n_common
         j = (i-1)*size_of_row+my_row
         k = mod(j-1,size_of_col)+1
         if ( k == my_col ) then
            diag_0 = i
            exit
         end if
      end do! i_1


      return
      end subroutine ! eigen_init
!----
      subroutine eigen_free(flag)
      implicit none

      include 'trd.h'

      integer :: flag

      real(8) ::    time_bcast, time_reduce, time_redist
      common /stat/ time_bcast, time_reduce, time_redist

      call mpi_comm_free(mpi_comm_col,ierr)
      call mpi_comm_free(mpi_comm_row,ierr)
!----
      call mpi_comm_free(es_ccomm_world,ierr)
      call mpi_comm_free(es_gcomm_world,ierr)
!----
#ifdef DETAIL
      if ( flag == 2 .and. myrank == 1 ) then
         print*, " "
         print*, "detail of exectime in \"eigen_dch\""
         print*, "  communication time in \"eigen_dch\""
      endif
      if ( flag == 3 .and. myrank == 1 ) then
         print*, "  communication time in \"eigen_hbk\""
      endif
      if ( flag >= 1 .and. flag <=3 .and. myrank == 1 ) then
         print*, "   bcast       :: ", time_bcast,"(sec)"
         print*, "   reduce      :: ", time_reduce,"(sec)"
         print*, "   redist      :: ", time_redist,"(sec)"
      endif
#endif

      return
      end subroutine ! eigen_free
!----
      subroutine datacast_dbl(u_y, u_x, u_t, u_s, n)
      implicit none

      integer, intent(in)    :: n
      real(8), intent(inout) :: u_y(1:n), u_x(1:n), u_t(1:n), u_s(1:n)

      include 'trd.h'

      integer                :: nx, ny, ic, i, j, k
      integer                :: req(1024), x_snod, y_snod
      integer                :: his_rank, her_rank

      real(8) ::    d1, d2
      real(8) ::    time_bcast, time_reduce, time_redist
      common /stat/ time_bcast, time_reduce, time_redist

#ifdef TIMER
      call mpi_barrier(mpi_comm_col,ierr)
      d1=mpi_wtime()
#endif

      if ( size_of_col == 1 ) then
         if ( size_of_row == 1 ) then
            u_y(1:n) = u_x(1:n)
         else
            ny = (n-1)/size_of_row+1
            do i=1,ny
               j = my_row+size_of_row*(i-1)
               u_y(i) = u_x(j)
            end do! i
         end if
         return
      end if

      if ( size_of_col == size_of_row ) then
         if ( my_col == my_row ) then
            u_y(1:n) = u_x(1:n)
         end if
         call bcast_dbl(u_y, n, my_row, mpi_comm_col)
         return
      end if

      x_snod = size_of_col/n_common
      y_snod = size_of_row/n_common

      if ( p0_(my_col) > 0 ) then

         nx = (n-p0_(my_col))/y_snod+1
         do i=1,nx
            j = p0_(my_col)+y_snod*(i-1)
            k = q0_(my_col)+x_snod*(i-1)
            u_t(i) = u_x(j)
            u_y(k) = u_x(j)
!----       print*,"i am",myrank,"sdrv ",myrank
         end do! i

         do ic=1,x_snod-1

!----
! receiving message : length    = ifloor(n/y_snod)
!                   : #sender   = size_of_col-1
!----
            his_rank = mod(my_col-1+size_of_col+ic*n_common,
     &                     size_of_col)+1
            ny = (n-p0_(his_rank))/y_snod+1
!----       print*,"i am",myrank,"recv ",(his_rank-1)*size_of_row+my_row
            call irecv_dbl(u_s, ny, his_rank, req(ic), mpi_comm_col)

!----
! sending message   : length    = ifloor(n/y_snod)
!                   : #receiver = size_of_col-1
!----
            her_rank = mod(my_col-1+size_of_col-ic*n_common,
     &                     size_of_col)+1
            call send_dbl(u_t, nx, her_rank, mpi_comm_col)
!----       print*,"i am",myrank,"send ",(her_rank-1)*size_of_row+my_row

            call wait_dbl(req(ic))

            do i=1,ny
               k = q0_(his_rank)+x_snod*(i-1)
               u_y(k) = u_s(i)
            end do! i

         end do! ic

         do ic=1,n_common-1
            her_rank = mod(my_col-1 +size_of_col +ic,size_of_col)+1
!----       print*,"i am",myrank,"send ",(her_rank-1)*size_of_row+my_row
            call isend_dbl(u_y, n, her_rank, req(ic), mpi_comm_col)
            call wait_dbl(req(ic))
         end do! ic

      else

         i = mod(my_row-1,n_common)
         j = mod(my_col-1,n_common)
         ic = mod(j-i+n_common,n_common)
         his_rank = mod(my_col-1 +size_of_col -ic,size_of_col)+1
!----    print*,"i am",myrank,"recv ",(his_rank-1)*size_of_row+my_row
         call recv_dbl(u_y, n, his_rank, mpi_comm_col)

      end if

#ifdef TIMER
      d2=mpi_wtime()
      time_redist=time_redist+(d2-d1)
#endif

      return
      end subroutine ! datacast_dbl
!----
      subroutine cyc1d_cyc2d(n,a,lda,b,ldb,size_of_col,size_of_row,iret)
      implicit none
!----
      integer                 :: n,lda,ldb,iret
      real(8)                 :: a(lda,1), b(ldb,1)
      integer                 :: size_of_col, size_of_row
!----
      real(8), pointer        :: work(:,:)
      integer                 :: world_size, my_rank, sq_size
      integer                 :: his_rank, her_rank
      integer                 :: his_posx, her_posx
      integer                 :: his_posy, her_posy
      integer                 :: my_posx, my_posy, my_posu, my_posv
      integer                 :: size_s,size_r
      integer                 :: nws, nqs
      integer                 :: i,j,k, ic,jc,kc,lc,id
!----
      integer                 :: status(mpi_status_size)
      integer                 :: type
      integer                 :: comm
      integer                 :: tag
      integer                 :: req
      integer                 :: error
      include 'commtxt.h'
!----
      iret = 0
!----
      if ( n > lda ) then
         iret = 1
         return
      end if
!----
!---- [*,world_size] --> [size_of_row,size_of_col]
!----     a(lda,nws) --> b(ldb,nqs)
!----
      comm = mpi_comm_eigen
      type = mpi_double_precision
      tag  = 10000
!----
      call mpi_comm_size(comm,world_size,error)
      call mpi_comm_rank(comm,my_rank,error)
      my_rank = my_rank+1
!---- sq_size=int(sqrt(dble(world_size)))
!----
!---- size_of_row = sq_size
!---- size_of_col = world_size / size_of_row
      if ( size_of_col*size_of_row /= world_size ) then
         iret = 3
         return
      end if
!----
      nws = (n-1)/world_size+1
      nqs = (n-1)/size_of_col+1
      if ( nqs > ldb ) then
         iret = 2
         return
      end if
!----
!---- if ( size_of_col > 1 ) then
      allocate(work(ldb,size_of_col), stat=error)
      if ( error /= 0 ) then
         iret = -1
         return
      end if
!---- end if
!----
      my_posx = mod(my_rank-1, size_of_col) + 1
      my_posy =    (my_rank-1)/size_of_col  + 1
      my_posu = mod(my_rank-1, size_of_row) + 1
      my_posv =    (my_rank-1)/size_of_row  + 1
!----
!----
      call mpi_barrier(comm,error)
!----
      do ic=1,nws
!----
         kc = my_rank+(ic-1)*world_size
         if ( kc <= n ) then
!----
            do id=1,size_of_col
               his_posx = mod(my_posv+id-2,size_of_col)+1
               size_s   =    (n-his_posx)/size_of_col+1
               do i=1,size_s
                  j = his_posx+(i-1)*size_of_col
                  work(i,id) = a(j,ic)
               end do
            end do
!----
         end if
!----
!---- bucket algorithm
!----
         do id=1,size_of_col
!----
            his_posx = mod(my_posv+id-2,size_of_col)+1
            his_rank =    (his_posx-1)*1+(my_posu-1)*size_of_col+1

            if ( kc <= n ) then
               size_s = (n-his_posx)/size_of_col+1
            else
               size_s = 0
            end if

            if ( size_s > 0 ) then
               call mpi_isend(work(1,id),size_s,type,
     &                        his_rank-1,tag,comm,req,error)
!----          print*,"i am",my_rank,"send ",his_rank,work(1:size_s,id)
            end if
!----
            her_posx = mod(my_posx-1+size_of_col-id+1,size_of_col)+1
            her_rank =    (her_posx-1)*size_of_row+(my_posy-1)*1+1

            lc = her_rank+(ic-1)*world_size

            if ( lc <= n ) then
               size_r = (n-my_posx)/size_of_col+1
            else
               size_r = 0
            end if

            if ( size_r > 0 ) then
               jc = her_posx+(ic-1)*size_of_col
               call mpi_recv(b(1,jc),size_r,type,
     &                       her_rank-1,tag,comm,status,error)
!----          print*,"i am",my_rank,"recv ",her_rank,b(1:size_r,jc),jc
            end if
!----
            if ( size_s > 0 ) then
               call mpi_wait(req,status,error)
            end if
!----
         end do

      end do
!----
!---- call mpi_barrier(comm,error)
!----
!---- if ( size_of_col > 1 ) then
      deallocate(work)
!---- end if
!----
      return
      end subroutine
!----
      subroutine cyc2d_cyc1d(n,a,lda,b,ldb,size_of_col,size_of_row,iret)
      implicit none
!----
      integer                 :: n,lda,ldb,iret
      real(8)                 :: a(lda,1), b(ldb,1)
      integer                 :: size_of_col, size_of_row
!----
      real(8), pointer        :: work(:,:)
      integer                 :: world_size, my_rank, sq_size
      integer                 :: his_rank, her_rank
      integer                 :: his_posx, her_posx
      integer                 :: his_posy, her_posy
      integer                 :: my_posx, my_posy, my_posu, my_posv
      integer                 :: size_s,size_r
      integer                 :: nws, nqs
      integer                 :: i,j,k, ic,jc,kc,lc,id
!----
      integer                 :: status(mpi_status_size)
      integer                 :: type
      integer                 :: comm
      integer                 :: tag
      integer                 :: req
      integer                 :: error

      include 'commtxt.h'
!----
      iret = 0
!----
      if ( n > lda ) then
         iret = 1
         return
      end if
!----
!---- [*,world_size] <-- [size_of_row,size_of_col]
!---- a(lda,nws) <-- b(ldb,nqs)
!----
      comm = mpi_comm_eigen
      type = mpi_double_precision
      tag  = 10000
!----
      call mpi_comm_size(comm,world_size,error)
      call mpi_comm_rank(comm,my_rank,error)
      my_rank = my_rank+1
!---- sq_size = int(sqrt(dble(world_size)))
!----
!---- size_of_row = sq_size
!---- size_of_col = world_size / size_of_row
      if ( size_of_col*size_of_row /= world_size ) then
         iret = 3
         return
      end if
!----
      nws = (n-1)/world_size+1
      nqs = (n-1)/size_of_col+1
      if ( nqs > ldb ) then
         iret = 2
         return
      end if
!----
!---- if ( size_of_col > 1 ) then
      allocate(work(ldb,size_of_col), stat=error)
      if ( error /= 0 ) then
         iret = -1
         return
      end if
!---- end if
!----
      my_posx = mod(my_rank-1, size_of_col) + 1
      my_posy =    (my_rank-1)/size_of_col  + 1
      my_posu = mod(my_rank-1, size_of_row) + 1
      my_posv =    (my_rank-1)/size_of_row  + 1
!----
!----
      call mpi_barrier(comm,error)
!----
      do ic=1,nws
!----
         kc = my_rank+(ic-1)*world_size
!----
!---- bucket algorithm
!----
         do id=1,size_of_col
!----
            her_posx = mod(my_posx-1+size_of_col-id+1,size_of_col)+1
            her_rank =    (her_posx-1)*size_of_row+(my_posy-1)*1+1

            lc = her_rank+(ic-1)*world_size
            if ( lc <= n ) then
               size_r = (n-my_posx)/size_of_col+1
            else
               size_r = 0
            end if

            if ( size_r > 0 ) then
               jc = her_posx+(ic-1)*size_of_col
               call mpi_isend(b(1,jc),size_r,type,
     &                        her_rank-1,tag,comm,req,error)
!----          print*,"i am",my_rank,"recv ",her_rank,b(1:size_r,jc)
            end if
!----
            his_posx = mod(my_posv+id-2,size_of_col)+1
            his_rank =    (his_posx-1)*1+(my_posu-1)*size_of_col+1

            if ( kc <= n ) then
               size_s = (n-his_posx)/size_of_col+1
            else
               size_s = 0
            end if

            if ( size_s > 0 ) then
               call mpi_recv(work(1,id),size_s,type,
     &                       his_rank-1,tag,comm,status,error)
!----          print*,"i am",my_rank,"send ",his_rank,work(1:size_s,id)
            end if
!----
            if ( size_r>0)then
               call mpi_wait(req,status,error)
            end if
!----
         end do

         if ( kc <= n ) then
            do id=1,size_of_col
               his_posx = mod(my_posv+id-2,size_of_col)+1
               size_s   =    (n-his_posx)/size_of_col+1
               do i=1,size_s
                  j = his_posx+(i-1)*size_of_col
                  a(j,ic) = work(i,id)
               end do
            end do
         end if
!----
      end do
!----
      call mpi_barrier(comm,error)
!----
!---- if ( size_of_col > 1 ) then
      deallocate(work)
!---- end if
!----
      return
      end subroutine

      end module communication_h
