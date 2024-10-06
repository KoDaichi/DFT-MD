       subroutine mat_set(n, a, nm, mtype, ndim)

       implicit double precision(a-h,o-z),integer(i-n)
       real(8)                :: a(1:nm,*)

       include 'mpif.h'
       include 'trd.h'


          call datacast_init(ndim)

          j_2=loop_c_start(1, x_nnod, x_inod)
          j_3=loop_c_end  (n, x_nnod, x_inod)
          i_2=loop_c_start(1, y_nnod, y_inod)
          i_3=loop_c_end  (n, y_nnod, y_inod)

          do i_1=i_2,i_3
             i = loop_c_l2g_depth(i_1, y_nnod, y_inod)
             if(mtype==0)then
                do j_1=j_2,j_3
                   j = loop_c_l2g_depth(j_1, x_nnod, x_inod)
                   a(j_1,i_1)=(n+1-Max(n+1-i,n+1-j))*1.0D+00
                enddo
             endif
!             if(mtype==1)then
!                a(1:n,i_1)=0D0
!                a(i,i_1)=2D0
!                if(i-1>=1)a(i-1,i_1)=-1D0
!                if(i+1<=n)a(i+1,i_1)=-1D0
!             endif
             if(mtype==1)then
                do j_1=j_2,j_3
                   j = loop_c_l2g_depth(j_1, x_nnod, x_inod)
                   if(i==j)then
                      a(j_1,i_1)=-7.2D0
                   else
                      a(j_1,i_1)=-3.0D0/(i-j)**2
                   endif
                enddo
             endif
             if(mtype==2)then
                do j_1=j_2,j_3
                   CALL RANDOM_NUMBER(t)
                   a(j_1,i_1)=1.D0-2*t
                enddo
             endif
          enddo

          call datacast_free(0)


       return
       end

