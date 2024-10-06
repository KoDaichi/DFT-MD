      subroutine eigen_hrd(
     &             n,ar,ai,lda,
     &             d_out,e_out,e2_out,
     &             tau,ldtau, m,zr,zi
     &           )
!----
      use communication_h, only : eigen_init
     &               , eigen_free
     &               , cyc1d_cyc2d
     &               , cyc2d_cyc1d
!----
      implicit none
!----
      include 'mpif.h'
!----
      integer  ::  n, lda, ldtau, m
      real(8)  ::  ar(*),ai(*)
      real(8)  ::  zr(*),zi(*)
      real(8)  ::  d_out(*),e_out(*),e2_out(*)
      real(8)  ::  tau(*)
!----
      include 'trd.h'
!----
      integer  ::  ierr
      integer  ::  n_, nm_

      real(8), pointer :: ar_(:), ai_(:)

      real(8)  :: d0, d1, d2, d3, d4

      integer  :: ldz, nz
!----
      call eigen_init(2)
!----
      d0=mpi_wtime()
      ldz = (n-1)/size_of_col+1
      nz = (n-1)/size_of_row+1
      d1=mpi_wtime()
      call cyc1d_cyc2d(n,ar,lda,zr,ldz, size_of_col, size_of_row, ierr)
      call cyc1d_cyc2d(n,ai,lda,zi,ldz, size_of_col, size_of_row, ierr)
      d2=mpi_wtime()
!----
      n_ = 0; nm_ = 1

      allocate(ar_(1),ai_(1))
      call eigen_hrd_main1(
     &       n,zr(1),zi(1),ldz,
     &       d_out(1),e_out(1),e2_out(1),tau(1),ldtau,m,
     &       ar_(1),ai_(1),n_,nm_
     &     )
      deallocate(ar_,ai_)

      d3=mpi_wtime()
      call cyc2d_cyc1d(n,ar,lda,zr,ldz, size_of_col, size_of_row, ierr)
      call cyc2d_cyc1d(n,ai,lda,zi,ldz, size_of_col, size_of_row, ierr)
      d4=mpi_wtime()
#ifdef DETAIL
      if(myrank==1) then
         print*,"  communication time in \"eigen_hrd\""
         print*,"   cyc1d_cyc2d :: ",d2-d1,"(sec)"
         print*,"   cyc2d_cyc1d :: ",d4-d3,"(sec)"
      endif
#endif
!----
      call eigen_free(1)
!----
#ifdef TIMER
      if(myrank==1)then
         print*,"Exectime of \"eigen_hrd\" routine =",d4-d0,"(sec)"
!----&            4*dble(n)**3*4/3/(s1-s0)*1d-6,"MFLOPS"
      endif
#endif

      return
      end subroutine
!----
      subroutine eigen_hrd_main1(
     &             n,ar,ai,lda,d_out,e_out,e2_out,tau,ldtau,m,
     &             ar_,ai_,n_,nm_
     &           )
!----
      implicit none
!----
      integer  ::  n,  lda, nv, ldtau, m
      real(8)  ::  ar(lda,*),ai(lda,*)
      real(8)  ::  d_out(n),e_out(n),e2_out(n)
      real(8)  ::  tau(ldtau,2)
      integer  ::  n_, nm_
      real(8)  ::  ar_(*), ai_(*)
!----
      include 'trd.h'
!----
      real(8), pointer ::  uv_x(:), uv_y(:)
      real(8), pointer ::  w(:)
      real(8), pointer ::  u_t(:), v_t(:)
      real(8), pointer ::  d_t(:)

      integer :: i,j,k
!----
      nv=lda
!----
      allocate(
     &         uv_x(4*nv*m),
     &         uv_y(4*nv*m),
     &         w(2*nv*m),
     &         u_t(1+8*2*max(nv+4,m*4)),
     &         v_t(1+8*2*max(nv+4,m*4)),
     &         d_t(nv)
     &        )
!----
      j = 0
!----
      call eigen_hrd_main2(
     &       n,ar,ai,lda,d_out,e_out,e2_out,
     &       tau,ldtau,nv,m,
     &       uv_x(1+0*nv*m), uv_x(1+1*nv*m),
     &       uv_x(1+2*nv*m), uv_x(1+3*nv*m),
     &       uv_y(1+0*nv*m), uv_y(1+1*nv*m),
     &       uv_y(1+2*nv*m), uv_y(1+3*nv*m),
     &       w(1), w(1+nv*m),
     &       u_t(1), u_t(1+8*max(nv+4,m*4)),
     &       v_t(1), v_t(1+8*max(nv+4,m*4)),
     &       d_t(1),
     &       ar_,ai_,n_,nm_
     &     )
!----
      deallocate(uv_x, uv_y, w, u_t, v_t, d_t )

      return
      end subroutine
!----
      subroutine eigen_hrd_main2(
     &             n,ar,ai,lda,d_out,e_out,e2_out,
     &             tau,ldtau,nv,m,
     &             ur_x, ui_x, vr_x, vi_x,
     &             ur_y, ui_y, vr_y, vi_y,
     &             wr, wi,
     &             ur_t, ui_t,
     &             vr_t, vi_t,
     &             d_t,
     &             ar_,ai_,n_,nm_
     &           )
!----
      use communication_h, only : reduce_dbl
     &               , bcast_dbl
     &               , get_loop_start
     &               , get_loop_end
     &               , get_owner_node
     &               , translate_l2g
     &               , translate_g2l
!----
      implicit integer         (i-n)
      implicit real(8)(a-h,o-z)
!----
      integer  ::  i,j,l,n,lda,ldtau,nv
      real(8)  ::  ar(lda,n),ai(lda,n)
      real(8)  ::  d_out(n),e_out(n),e2_out(n)
      real(8)  ::  tau(ldtau,2)
      real(8)  ::  g,h,scale
      real(8)  ::  sr,si,tr,ti
!----
      integer  ::  m_orig
      integer  ::  m,mm,m0
      integer  ::  i_block,i_base
!----
      integer  ::  y_root, y_pos
      integer  ::  x_root, x_pos
!----
      real(8)  ::  ur_x(nv,m),ui_x(nv,m)
      real(8)  ::  ur_y(nv,m),ui_y(nv,m)
      real(8)  ::  vr_x(nv,m),vi_x(nv,m)
      real(8)  ::  vr_y(nv,m),vi_y(nv,m)
      real(8)  ::  wr(nv,m),wi(nv,m)
      real(8)  ::  ur_t(1:8*max(nv+4,4*m))
      real(8)  ::  ui_t(1:8*max(nv+4,4*m))
      real(8)  ::  vr_t(1:8*max(nv+4,4*m))
      real(8)  ::  vi_t(1:8*max(nv+4,4*m))
      real(8)  ::  d_t(nv)
!----
      integer  ::  n_, nm_
      real(8)  ::  ar_(nm_,n_), ai_(nm_,n_)

      integer  :: local_rank, local_size
!----
      include 'mpif.h'
!$    include 'omp_lib.h'
      include 'trd.h'
!----
      real(8)  :: dd(10),d1,d2
!----
      dd=0d0
!----
!$omp parallel private(local_rank,local_size)
#ifdef DETAIL
!$omp master
      if ( myrank == 1 ) then
         local_rank = 0
         local_size = 1
!$       local_rank = omp_get_thread_num()
!$       local_size = omp_get_num_threads()
!$       if ( local_rank == 0 ) then
            print*,"num.of.process=",nprocs,
     &             "(",size_of_col,size_of_row,")"
!$          print*,"num.of.threads=",local_size
!$       endif
      endif
!$omp end master
#endif
      call eigen_hrd_Au(
     &       ar,ai,lda,
     &       ur_x,ui_x,ur_y,ui_y,vr_x,vi_x,
     &       ur_t,ui_t, vr_t,vi_t, d_t,
     &       -1, x_pos,y_pos, nv,k_1-k_2,
     &       ar_,ai_,n_,nm_
     &     )
      call eigen_hrd_2update(
     &       ar,ai,lda,
     &       ur_x,ui_x,ur_y,ui_y,
     &       vr_x,vi_x,vr_y,vi_y,   
     &       nv,m0,-1,
     &       ar_,ai_,n_,nm_
     &     )
!$omp end parallel

!----
      d_out(1:n)  = 0.0d0
      e_out(1:n)  = 0.0d0
      e2_out(1:n) = 0.0d0
      tau(1:n,1) = 1.0d0
      tau(1:n,2) = 0.0d0
!----
      if(diag_0>0)then
         i_2=diag_0
         i_3=((n)+size_of_row-my_row)/size_of_row
         i_4=size_of_col/n_common
         if(i_2<=n_)then
            do i_1=i_2,min(i_3,n_),i_4
               j=(i_1-1)*size_of_row+my_row
               j_1=(j-1)/size_of_col+1
               d_out(j)     = ar_(j_1,i_1)
               ai_(j_1,i_1) = 0.0d0
            enddo! i_1
            i_2=((n_-i_2)/i_4+1)*i_4+i_2
         endif
*poption indep(a,d_out)
         do i_1=i_2,i_3,i_4
            j=(i_1-1)*size_of_row+my_row
            j_1=(j-1)/size_of_col+1
            d_out(j)     = ar(j_1,i_1)
            ai(j_1,i_1)  = 0.0d0
         enddo! i_1
      endif
!----
      i_2=1
      i_3=((n)+size_of_row-my_row)/size_of_row
      do i_1=i_2,min(n_,i_3)
         j=(i_1-1)*size_of_row+my_row
         j_2=((j+1)+size_of_col-1-my_col)/size_of_col+1
         j_3=((n)+size_of_col-my_col)/size_of_col
         do j_1=j_2,nm_
            ar_(j_1,i_1) = 0.0d+00
            ai_(j_1,i_1) = 0.0d+00
         enddo! j_1
      enddo! i_1
      do i_1=max(i_2,n_+1),i_3
         j=(i_1-1)*size_of_row+my_row
         j_2=((j+1)+size_of_col-1-my_col)/size_of_col+1
         j_3=((n)+size_of_col-my_col)/size_of_col
         do j_1=j_2,lda
            ar(j_1,i_1)  = 0.0d+00
            ai(j_1,i_1)  = 0.0d+00
         enddo! j_1
      enddo! i_1
!----
      m_orig = m
      if ( m > n ) then
         m = n
      endif
      mm = (n-1)/m+1

      do  i_block = mm,max(1,3*(2-m)),-1

         i_base = (i_block-1)*m
         m0 = min(m,n-i_base)

         i = i_base + m0
         l = i_base + m0 - 1

         i_2=((i_base+1)+size_of_row-1-my_row)/size_of_row+1
         i_3=((i_base+m0)+size_of_row-my_row)/size_of_row
         j_2=1
         j_3=((i_base+m0)+size_of_col-my_col)/size_of_col

         do i_1=i_2,min(n_,i_3)
            j=(i_1-1)*size_of_row+my_row
            do j_1=j_2,min(nm_,j_3)
               wr(j_1,j-i_base) = ar_(j_1,i_1)
               wi(j_1,j-i_base) = ai_(j_1,i_1)
            enddo
            do j_1=min(nm_,j_3)+1,j_3
               wr(j_1,j-i_base) = 0.0d0+00
               wi(j_1,j-i_base) = 0.0d0+00
            enddo
         enddo
         do i_1=max(i_2,n_+1),i_3
            j=(i_1-1)*size_of_row+my_row
            do j_1=j_2,j_3
               wr(j_1,j-i_base) = ar(j_1,i_1)
               wi(j_1,j-i_base) = ai(j_1,i_1)
            enddo
         enddo
         if(diag_0>0)then
            i_2=diag_0
            i_3=((i_base+m0-1)+size_of_row-my_row)/size_of_row
            i_4=size_of_col/n_common
            if(i_2<=n_)then
               if(size_of_row==size_of_col)then
*poption indep(a,d_out)
                  do i_1=i_2,min(i_3,n_),i_4
                     d_t(i_1)     = ar_(i_1,i_1)
                     ar_(i_1,i_1) = 0.0d+00
                     ai_(i_1,i_1) = 0.0d+00
                  enddo! i_1
               else
*poption indep(a,d_out)
                  do i_1=i_2,min(i_3,n_),i_4
                     j=(i_1-1)*size_of_row+my_row
                     j_1 = (j-1)/size_of_col +1
                     d_t(i_1)     = ar_(j_1,i_1)
                     ar_(j_1,i_1) = 0.0d+00
                     ai_(j_1,i_1) = 0.0d+00
                  enddo! i_1
               endif
               i_2=((n_-i_2)/i_4+1)*i_4+i_2
            endif
            if(size_of_row==size_of_col)then
*poption indep(a,d_out)
               do i_1=i_2,i_3,i_4
                  d_t(i_1)    = ar(i_1,i_1)
                  ar(i_1,i_1) = 0.0d+00
                  ai(i_1,i_1) = 0.0d+00
               enddo! i_1
            else
*poption indep(a,d_out)
               do i_1=i_2,i_3,i_4
                  j=(i_1-1)*size_of_row+my_row
                  j_1 = (j-1)/size_of_col +1
                  d_t(i_1)    = ar(j_1,i_1)
                  ar(j_1,i_1) = 0.0d+00
                  ai(j_1,i_1) = 0.0d+00
               enddo! i_1
            endif
         endif

         k_2=m0
         k_3=max(1,3*(2-i_block))

         i=i_base+k_2
         l=i-1
         y_pos =   (i-1)/size_of_row +1
         x_pos =   (l-1)/size_of_col +1
         j_2=max(1,x_pos-m0)
         j_3=max(1,y_pos-m0)
         do  k_1 = k_3,k_2,1
            if(j_2<=x_pos)then
               ur_x(j_2:x_pos,k_1)=0.0d+00
               ui_x(j_2:x_pos,k_1)=0.0d+00
               vr_x(j_2:x_pos,k_1)=0.0d+00
               vi_x(j_2:x_pos,k_1)=0.0d+00
            endif
            if(j_3<=y_pos)then
               ur_y(j_3:y_pos,k_1)=0.0d+00
               ui_y(j_3:y_pos,k_1)=0.0d+00
               vr_y(j_3:y_pos,k_1)=0.0d+00
               vi_y(j_3:y_pos,k_1)=0.0d+00
            endif
         enddo! k_1

         do  k_1 = k_2,k_3,-1

            i = i_base+k_1
            l = i-1

            i_2 = get_loop_start(1, size_of_row,my_row)
            i_3 = get_loop_end  (l, size_of_row,my_row)
            j_2 = get_loop_start(1, size_of_col,my_col)
            j_3 = get_loop_end  (l, size_of_col,my_col)
!----
            y_root = get_owner_node (i, size_of_row,my_row)
            y_pos  = translate_g2l(i, size_of_row,my_row)

            x_root = get_owner_node (l, size_of_col,my_col)
            x_pos  = translate_g2l(l, size_of_col,my_col)

            d1=mpi_wtime()
            call eigen_hrd_u(
     &             wr(1,k_1),wi(1,k_1),
     &             ur_x(1,k_1),ui_x(1,k_1),
     &             ur_y(1,k_1),ui_y(1,k_1),
     &             ur_t(1),vr_t(1),
     &             x_root,y_root,x_pos, j_3,lda,nv
     &           )

            scale = ur_t(x_pos+1)
            if ( scale == 0.0d+00 ) then
               do j=1,x_pos
                  ur_x(j,k_1)=0.0d+00
                  ui_x(j,k_1)=0.0d+00
                  vr_x(j,k_1)=0.0d+00
                  vi_x(j,k_1)=0.0d+00
               enddo
               do j=1,y_pos
                  ur_y(j,k_1)=0.0d+00
                  ui_y(j,k_1)=0.0d+00
                  vr_y(j,k_1)=0.0d+00
                  vi_y(j,k_1)=0.0d+00
               enddo
               goto 1000
            endif
            d2=mpi_wtime()
            dd(1)=dd(1)+(d2-d1)

            if ( vr_t(x_pos+3) == 0.0d+00 ) then
               tau(l,1) = -tau(i,1)
               tau(l,2) = -tau(i,2)
               wr(x_pos,k_1) = ur_t(x_pos+3)
            else
               sr       =  vr_t(x_pos+1)
               si       =  vr_t(x_pos+2)
               tr       =  si*tau(i,2) - sr*tau(i,1)
               ti       = -sr*tau(i,2) - si*tau(i,1)
               tau(l,1) =  tr / vr_t(x_pos+3)
               tau(l,2) =  ti / vr_t(x_pos+3)
            endif

            beta  = ur_t(x_pos+4)
            if ( x_root == my_col .and. y_root == my_row ) then
               e_out(i)      = scale*ur_t(x_pos+3)
               e2_out(i)     = scale*scale*ur_t(x_pos+2)
               wr(x_pos,k_1) = scale*ur_t(x_pos)
               wi(x_pos,k_1) = scale*vr_t(x_pos)
            endif
            i_1 = get_owner_node(i, size_of_col,my_col)
            if ( i_1 == my_col .and. y_root == my_row ) then
               j_1 = translate_g2l(i, size_of_col,my_col)
               wi(j_1,k_1) = scale*sqrt(ur_t(x_pos+4))
            endif
!----
!---- v = au
!----
            d1=mpi_wtime() 
!$omp parallel
            call eigen_hrd_Au(
     &             ar(1,1),ai(1,1),lda,
     &             ur_x(1,k_1),ui_x(1,k_1),
     &             ur_y(1,k_1),ui_y(1,k_1),
     &             vr_x(1,k_1),vi_x(1,k_1),
     &             ur_t(1),ui_t(1), vr_t(1),vi_t(1),
     &             d_t(1),
     &             i_3, x_pos,y_pos, nv,k_1-k_2,
     &             ar_(1,1),ai_(1,1),n_,nm_
     &           )
!$omp end parallel

            d2=mpi_wtime() 
            dd(2)=dd(2)+d2-d1
!----
!---- v = (a+uv+vu)u = au+(vu)u+(uu)v
!----
            d1=mpi_wtime() 
            if(k_2>k_1)then
*soption unroll(2)
               do l_1=k_2,k_1+1,-1
                  j=k_2-l_1
                  sr = 0.0d+00
                  si = 0.0d+00
                  tr = 0.0d+00
                  ti = 0.0d+00
*soption unroll(4)
!dir$ ivdep
!dir$ vector always
                  do  j_1 = 1, x_pos
                     wr_x0 = ur_x(j_1,k_1)
                     wi_x0 = ui_x(j_1,k_1)
                     ur_x0 = ur_x(j_1,l_1)
                     ui_x0 = ui_x(j_1,l_1)
                     vr_x0 = vr_x(j_1,l_1)
                     vi_x0 = vi_x(j_1,l_1)

                     sr = sr
     &                      + vr_x0*wr_x0
                     sr = sr
     &                      + vi_x0*wi_x0
                     si = si
     &                      + vi_x0*wr_x0
                     si = si
     &                      - vr_x0*wi_x0
                     tr = tr
     &                      + ur_x0*wr_x0
                     tr = tr
     &                      + ui_x0*wi_x0
                     ti = ti
     &                      + ui_x0*wr_x0
                     ti = ti
     &                      - ur_x0*wi_x0
                  enddo! k
                  ur_t(4*j+1) = sr
                  ur_t(4*j+2) = si
                  ur_t(4*j+3) = tr
                  ur_t(4*j+4) = ti
               enddo! j0
               j=(k_2-k_1)*4
               call reduce_dbl(ur_t(1), vr_t(1), j, mpi_comm_col)
*soption unroll(2)
               do l_1=k_2,k_1+1,-1
                  j=k_2-l_1
                  sr = ur_t(4*j+1)
                  si = ur_t(4*j+2)
                  tr = ur_t(4*j+3)
                  ti = ur_t(4*j+4)
*soption unroll(4)
!dir$ ivdep
!dir$ vector always
                  do  j_1 = 1, x_pos
                     wr_x0 = 0
                     wi_x0 = 0
                     ur_x0 = ur_x(j_1,l_1)
                     ui_x0 = ui_x(j_1,l_1)
                     vr_x0 = vr_x(j_1,l_1)
                     vi_x0 = vi_x(j_1,l_1)

                     wr_x0 = wr_x0
     &                            - tr*vr_x0
                     wi_x0 = wi_x0
     &                            + ti*vr_x0
                     wr_x0 = wr_x0
     &                            - sr*ur_x0
                     wi_x0 = wi_x0
     &                            + si*ur_x0
                     wr_x0 = wr_x0
     &                            - ti*vi_x0
                     wi_x0 = wi_x0
     &                            - tr*vi_x0
                     wr_x0 = wr_x0
     &                            - si*ui_x0
                     wi_x0 = wi_x0
     &                            - sr*ui_x0

                     vr_x(j_1,k_1) = wr_x0+vr_x(j_1,k_1)
                     vi_x(j_1,k_1) = wi_x0+vi_x(j_1,k_1)
                  enddo! k
               enddo! j0
            endif
            d2=mpi_wtime()
            dd(5)=dd(5)+d2-d1

            d1=mpi_wtime()
            call  eigen_hrd_v(
     &              ur_x(1,k_1),ui_x(1,k_1),
     &              vr_x(1,k_1),vi_x(1,k_1),
     &              vr_y(1,k_1),vi_y(1,k_1),
     &              ur_t(1),vr_t(1),
     &              beta, j_3,x_pos,nv
     &            )
            d2=mpi_wtime()
            dd(4)=dd(4)+d2-d1

            d1=mpi_wtime()
            if(k_1>1)then

               i_2 = get_loop_start(i_base+1,     size_of_row,my_row)
               i_3 = get_loop_end  (i_base+k_1-1, size_of_row,my_row)

               do i_1=i_2,i_3

                  j   = translate_l2g(i_1, size_of_row,my_row)
                  l_1 = j-i_base

                  ur_y0=ur_y(i_1,k_1)
                  vr_y0=vr_y(i_1,k_1)
                  ui_y0=ui_y(i_1,k_1)
                  vi_y0=vi_y(i_1,k_1)

*soption unroll(4)
!dir$ ivdep
!dir$ vector always
                  do j_1=j_2,j_3

                     wr_x0 = 0
                     wi_x0 = 0
                     ur_x0 = ur_x(j_1,k_1)
                     ui_x0 = ui_x(j_1,k_1)
                     vr_x0 = vr_x(j_1,k_1)
                     vi_x0 = vi_x(j_1,k_1)

                     wr_x0 = wr_x0
     &                            - vr_x0*ur_y0
                     wi_x0 = wi_x0
     &                            - vi_x0*ur_y0
                     wr_x0 = wr_x0
     &                            - ur_x0*vr_y0
                     wi_x0 = wi_x0
     &                            - ui_x0*vr_y0
                     wr_x0 = wr_x0
     &                            - vi_x0*ui_y0
                     wi_x0 = wi_x0
     &                            + vr_x0*ui_y0
                     wr_x0 = wr_x0
     &                            - ui_x0*vi_y0
                     wi_x0 = wi_x0
     &                            + ur_x0*vi_y0

                     wr(j_1,l_1) = wr_x0+wr(j_1,l_1)
                     wi(j_1,l_1) = wi_x0+wi(j_1,l_1)

                  enddo! j_1
               enddo! i_1
            endif
            d2=mpi_wtime()
            dd(6)=dd(6)+d2-d1

1000        continue
         enddo! k_1

         i_2 = get_loop_start(i_base+1,  size_of_row,my_row)
         i_3 = get_loop_end  (i_base+m0, size_of_row,my_row)
         do i_1=i_2,min(n_,i_3)
            j   = translate_l2g(i_1, size_of_row,my_row)
            j_2 = 1
            j_3 = get_loop_end  (j, size_of_col,my_col)
            do j_1=j_2,j_3
               ar_(j_1,i_1) = wr(j_1,j-i_base)
               ai_(j_1,i_1) = wi(j_1,j-i_base)
            enddo! j_1
         enddo! i_1
         do i_1=max(i_2,n_+1),i_3
            j   = translate_l2g(i_1, size_of_row,my_row)
            j_2 = 1
            j_3 = get_loop_end  (j, size_of_col,my_col)
            do j_1=j_2,j_3
               ar(j_1,i_1) = wr(j_1,j-i_base)
               ai(j_1,i_1) = wi(j_1,j-i_base)
            enddo! j_1
         enddo! i_1

         if ( diag_0 > 0 ) then
            i_2 = diag_0
            i_3 = get_loop_end  (i_base, size_of_row,my_row)
            i_4 = size_of_col/n_common
            if ( i_2 <= n_ ) then
               if ( size_of_row == size_of_col ) then
                  do i_1=i_2,min(i_3,n_),i_4
                     j_1 = i_1
                     ar_(j_1,i_1) = d_t(i_1)
                  enddo! i_1
               else
                  do i_1=i_2,min(i_3,n_),i_4
                     j   = (i_1-1)*size_of_row+my_row
                     j_1 = (j-1)/size_of_col +1
                     ar_(j_1,i_1) = d_t(i_1)
                  enddo! i_1
               endif
               i_2=((n_-i_2)/i_4+1)*i_4+i_2
            endif
            if ( size_of_row == size_of_col ) then
               do i_1=i_2,i_3,i_4
                  j_1 = i_1
                  ar(j_1,i_1)  = d_t(i_1)
               enddo! i_1
            else
               do i_1=i_2,i_3,i_4
                  j   = (i_1-1)*size_of_row+my_row
                  j_1 = (j-1)/size_of_col +1
                  ar(j_1,i_1)  = d_t(i_1)
               enddo! i_1
            endif
         endif

         if ( i_block > 1 ) then
            i_2 = 1
!----       i_3 = get_loop_end  (i_base, size_of_row,my_row)
            i_3 = i_base
            d1=mpi_wtime() 
!$omp parallel
            call eigen_hrd_2update(
     &             ar(1,1),ai(1,1),lda,
     &             ur_x(1,1),ui_x(1,1),
     &             ur_y(1,1),ui_y(1,1),
     &             vr_x(1,1),vi_x(1,1),
     &             vr_y(1,1),vi_y(1,1),   
     &             nv,m0,i_3,
     &             ar_(1,1),ai_(1,1),n_,nm_
     &           )
!$omp end parallel
            d2=mpi_wtime() 
            dd(3)=dd(3)+d2-d1

         endif

      enddo! i_block

      if ( n >= 2 ) then
         i = 2; l = i-1
         y_root = get_owner_node (i, size_of_row,my_row)
         x_root = get_owner_node (l, size_of_col,my_col)
         if ( y_root == my_row .and. x_root == my_col ) then
            i_1 = translate_g2l(i, size_of_row,my_row)
            j_1 = translate_g2l(l, size_of_col,my_col)
            if(i_1<=n_)then
               ur_t(1) = ar_(j_1,i_1)
               ur_t(2) = ai_(j_1,i_1)
            else
               ur_t(1) = ar(j_1,i_1)
               ur_t(2) = ai(j_1,i_1)
            endif
         endif
         call bcast_dbl(ur_t(1), 2, x_root, mpi_comm_col)
         call bcast_dbl(ur_t(1), 2, y_root, mpi_comm_row)
         scale = abs(ur_t(1)) + abs(ur_t(2))

         if ( scale > 0.0d0 ) then
            sr = ur_t(1) / scale
            si = ur_t(2) / scale
            h = sr**2 + si**2
            g = sqrt(h)
            tau(l,1) =  (si * tau(i,2) - sr * tau(i,1)) / g
            tau(l,2) = -(sr * tau(i,2) + si * tau(i,1)) / g
            y_root = get_owner_node (i, size_of_row,my_row)
            x_root = get_owner_node (l, size_of_col,my_col)
            if ( y_root == my_row .and. x_root == my_col ) then
               i_1 = translate_g2l(i, size_of_row,my_row)
               j_1 = translate_g2l(l, size_of_col,my_col)
               e_out(i)  = scale * g
               e2_out(i) = scale * scale * h
               if ( i_1 <= n_ ) then
                  ar_(j_1,i_1) = scale * 2 * sr
                  ai_(j_1,i_1) = scale * 2 * si
               else
                  ar(j_1,i_1)  = scale * 2 * sr
                  ai(j_1,i_1)  = scale * 2 * si
               endif
            endif
            y_root = get_owner_node (i, size_of_row,my_row)
            x_root = get_owner_node (i, size_of_col,my_col)
            if ( y_root == my_row .and. x_root == my_col ) then
               i_1 = translate_g2l(i, size_of_row,my_row)
               j_1 = translate_g2l(i, size_of_col,my_col)
               if ( i_1 <= n_ ) then
                  ai_(j_1,i_1) = scale * sqrt(2*h)
               else
                  ai(j_1,i_1)  = scale * sqrt(2*h)
               endif
            endif
            y_root = get_owner_node (l, size_of_row,my_row)
            x_root = get_owner_node (l, size_of_col,my_col)
            if ( y_root == my_row .and. x_root == my_col ) then
               i_1 = translate_g2l(l, size_of_row,my_row)
               j_1 = translate_g2l(l, size_of_col,my_col)
               if ( i_1 <= n_ ) then
                  ai_(j_1,i_1) = 0.0d0
               else
                  ai(j_1,i_1)  = 0.0d0
               endif
            endif
         endif
      endif

      if ( diag_0 > 0 ) then
         i_2 = diag_0
         i_3 = get_loop_end  (n, size_of_row,my_row)
         i_4 = size_of_col/n_common
         if ( i_2 <= n_ ) then
         do i_1=i_2,min(i_3,n_),i_4
            j   = (i_1-1)*size_of_row+my_row
            j_1 = (j-1)/size_of_col +1
            t1           = d_out(j)
            d_out(j)     = ar_(j_1,i_1)
            ar_(j_1,i_1) = t1
         enddo! i_1
         i_2=((n_-i_2)/i_4+1)*i_4+i_2
         endif
         do i_1=i_2,i_3,i_4
            j   = (i_1-1)*size_of_row+my_row
            j_1 = (j-1)/size_of_col +1
            t1           = d_out(j)
            d_out(j)     = ar(j_1,i_1)
            ar(j_1,i_1)  = t1
         enddo! i_1
      endif

      call reduce_dbl(d_out(1), ur_t(1), n, mpi_comm_eigen)
      call reduce_dbl(e_out(1), ur_t(1), n, mpi_comm_eigen)
      call reduce_dbl(e2_out(1),ur_t(1), n, mpi_comm_eigen)

!$omp parallel
      call eigen_hrd_Au(
     &       ar,ai,lda,
     &       ur_x,ui_x,ur_y,ui_y,vr_x,vi_x,
     &       ur_t,ui_t, vr_t,vi_t, d_t,
     &       -2, x_pos,y_pos, nv,k_1-k_2,
     &       ar_,ai_,n_,nm_
     &     )
      call eigen_hrd_2update(
     &       ar,ai,lda,
     &       ur_x,ui_x,ur_y,ui_y,
     &       vr_x,vi_x,vr_y,vi_y,   
     &       nv,m0,-2,
     &       ar_,ai_,n_,nm_
     &     )
!$omp end parallel

!----
#ifdef DETAIL
      if(myrank==1)then
         print*," "
         print*,"detail of exectime in \"eigen_hrd\""
         print*,"   calc (u,beta)    ",dd(1),"(sec)"
         print*,"   mat-vec (au)     ",dd(2),"(sec)"
         print*,"   2update (a-uv-vu)",dd(3),"(sec)"
!----    print*,"   mat-vec (au)     ",dd(2),(dble(n)**3*8d-9/3)/dd(2)
!----    print*,"   2update (a-uv-vu)",dd(3),(dble(n)**3*8d-9/3)/dd(3)
         print*,"   calc v           ",dd(4),"(sec)"
         print*,"   v=v-(uv+vu)u     ",dd(5),"(sec)"
         print*,"   uv post reduction",dd(6),"(sec)"
      endif
#endif
!----
      m = m_orig

      return
      end

