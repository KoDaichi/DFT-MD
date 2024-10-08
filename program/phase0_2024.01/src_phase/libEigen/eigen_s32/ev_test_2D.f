       subroutine ev_test_2D(n, a, nma, w, z, nmz, ret)

       integer                :: n, nma, nmz
       real(8)                :: a(*), w(*), z(*), ret

       include 'mpif.h'
       include 'trd.h'

       call datacast_init(2)
       call MPI_Barrier(MPI_COMM_EIGEN,ierr)
       call ev_check_2D_body(n, a(1), nma, w(1), z(1), nmz, ret)
       call MPI_Barrier(MPI_COMM_EIGEN,ierr)
       call datacast_free(0)

       return
       end
*
*
*
       subroutine ev_check_2D_body(n, a, nma, w, z, nmz, ret)

       implicit double precision(a-h,o-z),integer(i-n)

       integer                :: n, nma, nmz
       real(8)                :: a(1:nma,*)
       real(8)                :: w(1:n)
       real(8)                :: z(1:nmz,*)
       real(8)                :: ret

       real(8), pointer       :: b(:)
       real(8), pointer       :: c1(:), c2(:)
       real(8), pointer       :: d(:)
       real(8), pointer       :: e(:)

       real(8)                :: m_epsilon
       real(8), external      :: machine_epsilon

       include 'mpif.h'
       include 'trd.h'
*
       m_epsilon = machine_epsilon()
*
       ra=0d0
       rz=0d0
       rr=0d0
       r =0d0

       j_2 = loop_c_start(1, x_nnod,x_inod)
       j_3 = loop_c_end  (n, x_nnod,x_inod)
       j_4 = loop_c_g2l_depth(n, x_nnod,x_inod)
       i_2 = loop_c_start(1, y_nnod,y_inod)
       i_3 = loop_c_end  (n, y_nnod,y_inod)
       i_4 = loop_c_g2l_depth(n, y_nnod,y_inod)
       k_4 = MAX(i_4,j_4)+8

       larray = n
       allocate(b(1:larray))

       larray = j_4
       allocate(c1(1:larray), c2(1:larray))

       larray = MAX(n+1,3*k_4)
       allocate(d(1:larray), e(1:larray))

       b(1:n)=0d0

       i_m = MOD(i_3-i_2+1,4)+i_2
       do i=1,n,1

          jj_0 = loop_c_node(i, x_nnod,x_inod)
          ii_0 = loop_c_node(i, y_nnod,y_inod)
          jj_1 = loop_c_g2l_depth(i, x_nnod,x_inod)
          ii_1 = loop_c_g2l_depth(i, y_nnod,y_inod)

          if ( ii_0 == y_inod ) then
             do j_1 = j_2, j_3
                c1(j_1) = z(j_1,ii_1)
             end do
          end if
          call MPI_Bcast(c1(1),j_4,MPI_DOUBLE_PRECISION,
     &                   ii_0-1,y_COMM_WORLD,ierr)

          do i_1 = i_2, i_m-1
             t0 = 0.0D+00
             s0 = 0.0D+00
             u0 = 0.0D+00
             do j_1 = j_2, j_3
                a0 = a(j_1,i_1)
                b0 = c1(j_1)
                t0 = t0 + b0 * a0          ! t = Az(i)
                s0 = s0 + b0 * z(j_1,i_1)  ! s = Zz(i)
                u0 = u0 + ABS(a0)          ! u = |A|_1
             end do! k
             d(1+(i_1-1)*3) = t0
             d(2+(i_1-1)*3) = s0
             d(3+(i_1-1)*3) = u0
          end do
          do i_1 = i_m, i_3, 4
             t0 = 0.0D+00
             t1 = 0.0D+00
             t2 = 0.0D+00
             t3 = 0.0D+00
             s0 = 0.0D+00
             s1 = 0.0D+00
             s2 = 0.0D+00
             s3 = 0.0D+00
             u0 = 0.0D+00
             u1 = 0.0D+00
             u2 = 0.0D+00
             u3 = 0.0D+00
             do j_1 = j_2, j_3
                a0 = a(j_1,i_1+0)
                a1 = a(j_1,i_1+1)
                a2 = a(j_1,i_1+2)
                a3 = a(j_1,i_1+3)
                b0 = c1(j_1+0)
                t0 = t0 + b0 * a0          ! t = Az(i)
                t1 = t1 + b0 * a1          ! t = Az(i)
                t2 = t2 + b0 * a2          ! t = Az(i)
                t3 = t3 + b0 * a3          ! t = Az(i)
                s0 = s0 + b0 * z(j_1,i_1+0)! s = Zz(i)
                s1 = s1 + b0 * z(j_1,i_1+1)! s = Zz(i)
                s2 = s2 + b0 * z(j_1,i_1+2)! s = Zz(i)
                s3 = s3 + b0 * z(j_1,i_1+3)! s = Zz(i)
                u0 = u0 + ABS(a0)          ! u = |A|_1
                u1 = u1 + ABS(a1)          ! u = |A|_1
                u2 = u2 + ABS(a2)          ! u = |A|_1
                u3 = u3 + ABS(a3)          ! u = |A|_1
             end do! k
             d(1+(i_1-1)*3) = t0
             d(2+(i_1-1)*3) = s0
             d(3+(i_1-1)*3) = u0
             d(1+(i_1  )*3) = t1
             d(2+(i_1  )*3) = s1
             d(3+(i_1  )*3) = u1
             d(1+(i_1+1)*3) = t2
             d(2+(i_1+1)*3) = s2
             d(3+(i_1+1)*3) = u2
             d(1+(i_1+2)*3) = t3
             d(2+(i_1+2)*3) = s3
             d(3+(i_1+2)*3) = u3
          end do

          call MPI_Allreduce(d(1),e(1),3*(i_3-i_2+1),
     &                    MPI_DOUBLE_PRECISION,
     &                    MPI_SUM,x_COMM_WORLD,ierr)

          call datacast_dbl(d(1), c1(1), d(1+k_4), d(1+2*k_4), j_4)

          if ( jj_0 == x_inod ) then
             do i_1 = i_2, i_m-1
                t0 = e(1+(i_1-1)*3)
                s0 = e(2+(i_1-1)*3)
                u0 = e(3+(i_1-1)*3)
                t0 = t0 - w(i) * d(i_1)     ! t = Az(i) - wz(i)
             if ( ii_0 == y_inod .AND. i_1 ==ii_1 ) then
                s0 = s0 - 1.0D+00           ! s = Zz(i) - 1(i)
             end if
                b(i) = b(i) + ABS(t0)       ! b(i)=|Az(i)-wz(i)|_1
                rz   = rz + s0**2           ! rz  =|ZZ-I|_F
                ra   = MAX(ra,u0)           ! ra  =|A|_1
             end do
             do i_1 = i_m, i_3, 4
                t0 = e(1+(i_1-1)*3)
                s0 = e(2+(i_1-1)*3)
                u0 = e(3+(i_1-1)*3)
                t1 = e(1+(i_1  )*3)
                s1 = e(2+(i_1  )*3)
                u1 = e(3+(i_1  )*3)
                t2 = e(1+(i_1+1)*3)
                s2 = e(2+(i_1+1)*3)
                u2 = e(3+(i_1+1)*3)
                t3 = e(1+(i_1+2)*3)
                s3 = e(2+(i_1+2)*3)
                u3 = e(3+(i_1+2)*3)
                t0 = t0 - w(i) * d(i_1+0)   ! t = Az(i) - wz(i)
                t1 = t1 - w(i) * d(i_1+1)   ! t = Az(i) - wz(i)
                t2 = t2 - w(i) * d(i_1+2)   ! t = Az(i) - wz(i)
                t3 = t3 - w(i) * d(i_1+3)   ! t = Az(i) - wz(i)
             if ( ii_0 == y_inod ) then
             if ( i_1+0 ==ii_1 ) then
                s0 = s0 - 1.0D+00           ! s = Zz(i) - 1(i)
             else if ( i_1+1 ==ii_1 ) then
                s1 = s1 - 1.0D+00           ! s = Zz(i) - 1(i)
             else if ( i_1+2 ==ii_1 ) then
                s2 = s2 - 1.0D+00           ! s = Zz(i) - 1(i)
             else if ( i_1+3 ==ii_1 ) then
                s3 = s3 - 1.0D+00           ! s = Zz(i) - 1(i)
             end if
             end if
                b(i) = b(i) + ABS(t0)       ! b(i)=|Az(i)-wz(i)|_1
                b(i) = b(i) + ABS(t1)       ! b(i)=|Az(i)-wz(i)|_1
                b(i) = b(i) + ABS(t2)       ! b(i)=|Az(i)-wz(i)|_1
                b(i) = b(i) + ABS(t3)       ! b(i)=|Az(i)-wz(i)|_1
                rz   = rz + s0**2           ! rz  =|ZZ-I|_F
                rz   = rz + s1**2           ! rz  =|ZZ-I|_F
                rz   = rz + s2**2           ! rz  =|ZZ-I|_F
                rz   = rz + s3**2           ! rz  =|ZZ-I|_F
                ra   = MAX(ra,u0)           ! ra  =|A|_1
                ra   = MAX(ra,u1)           ! ra  =|A|_1
                ra   = MAX(ra,u2)           ! ra  =|A|_1
                ra   = MAX(ra,u3)           ! ra  =|A|_1
             end do
          end if

       end do

*>
       d(1:n)=b(1:n)
       d(n+1)=rz
       call MPI_Allreduce(d(1),e(1),n+1,MPI_DOUBLE_PRECISION,
     &                    MPI_SUM,MPI_COMM_EIGEN,ierr)
       b(1:n)=e(1:n)
       rz    =e(n+1)

       r=ra
       call MPI_Allreduce(ra,r,1,MPI_DOUBLE_PRECISION,
     &                    MPI_MAX,MPI_COMM_EIGEN,ierr)
       ra=r

       rr=0d0
       do i=1,n
          if(b(i)>rr)then
             rr=b(i); k=i
          endif
       enddo! i

       if(TRD_inod==1)then
          print*,"max|Ax-wx|=",rr,k
          print*,"|A|=",ra
          print*,"epsilon=",m_epsilon
          print*,"max|Ax-wx|/Ne|A|=",rr/(n*m_epsilon*ra)
          print*,"|ZZ-I|=",SQRT(rz),rz
       endif
*
       deallocate(e)
       deallocate(b)
       deallocate(c1,c2)
       deallocate(d)
*
       return
       end
