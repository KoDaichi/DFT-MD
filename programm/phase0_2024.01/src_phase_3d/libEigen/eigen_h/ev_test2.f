       subroutine ev_test2(ar,ai,w,zr,zi,n,nm)
!
      use communication_h, only : eigen_init
     &               , eigen_free
!
       integer                :: n, nm
       real(8)                :: ar(*), ai(*), w(*), zr(*), zi(*)

       include 'mpif.h'
       include 'trd.h'

       call eigen_init(1)
       call MPI_Barrier(MPI_COMM_EIGEN,ierr)
       call ev_check2_body(ar(1),ai(1),w(1),zr(1),zi(1),n,nm)
       call MPI_Barrier(MPI_COMM_EIGEN,ierr)
       call eigen_free(0)

       return
       end
*
*
*
       subroutine ev_check2_body(ar,ai,w,zr,zi,n,nm)
!
      use communication_h, only : loop_c
!
       implicit double precision(a-h,o-z),integer(i-n)

       integer                :: n, nm
       real(8)                :: ar(1:nm,*), ai(1:nm,*)
       real(8)                :: w(1:n)
       real(8)                :: zr(1:nm,*), zi(1:nm,*)

       real(8)                :: b(1:n)
       real(8)                :: b1r(1:n), b1i(1:n)

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
       b(1:n)=0d0

       do i=1,n,1

             i_3=MOD(i-1,nprocs)+1
             if(i_3==myrank)then
                i_1=(i-1)/nprocs+1
                do k=1,n
                   b1r(k)=zr(k,i_1)
                   b1i(k)=zi(k,i_1)
                enddo
             endif
             call MPI_Bcast(b1r(1),n,MPI_DOUBLE_PRECISION,i_3-1
     1       ,MPI_COMM_EIGEN,ierr)
             call MPI_Bcast(b1i(1),n,MPI_DOUBLE_PRECISION,i_3-1
     1       ,MPI_COMM_EIGEN,ierr)

             i_2=1
             i_3=n
             call loop_c(i_2,i_3,i_4,n, nprocs,myrank)
             do i_1=i_2,i_3
                j=(i_1-1)*nprocs+myrank
                t1=0.0D+00
                t2=0.0D+00
                s1=0.0D+00
                s2=0.0D+00
                u =0.0D+00
                do k=1,n,1
                   t1=t1+b1r(k)*ar(k,i_1) ! t = Az(i)
     &                  +b1i(k)*ai(k,i_1) ! t = Az(i)
                   t2=t2-b1r(k)*ai(k,i_1) ! t = Az(i)
     &                  +b1i(k)*ar(k,i_1) ! t = Az(i)
                   s1=s1+b1r(k)*zr(k,i_1) ! s = Zz(i)
     &                  +b1i(k)*zi(k,i_1) ! s = Zz(i)
                   s2=s2+b1r(k)*zi(k,i_1) ! s = Zz(i)
     &                  -b1i(k)*zr(k,i_1) ! s = Zz(i)
                   u =u +SQRT(ar(k,i_1)**2   ! u = |A|_1
     &                       +ai(k,i_1)**2)  ! u = |A|_1
                enddo! k
                t1=t1-w(i)*b1r(j) ! t = Az(i) - wz(i)
                t2=t2-w(i)*b1i(j) ! t = Az(i) - wz(i)
                if(i==j) s1=s1-1.0D+00 ! s = Zz(i)-1(i)
                b(i)=b(i)+SQRT(t1**2+t2**2) ! b(i)=|Az(i)-wz(i)|_1
                rz=rz+s1**2+s2**2
                ra=MAX(ra,u)
             enddo
       enddo
*>

       b1r(1:n)=b(1:n)
       call MPI_Allreduce(b,b1r,n,MPI_DOUBLE_PRECISION,
     &                    MPI_SUM,MPI_COMM_EIGEN,ierr)
       b(1:n)=b1r(1:n)
       r=ra
       call MPI_Allreduce(ra,r,1,MPI_DOUBLE_PRECISION,
     &                    MPI_MAX,MPI_COMM_EIGEN,ierr)
       ra=r
       r=rz
       call MPI_Allreduce(rz,r,1,MPI_DOUBLE_PRECISION,
     1                    MPI_SUM,MPI_COMM_EIGEN,ierr)
       rz=r

       rr=0d0
       do i=1,n
          if(b(i)>rr)then
             rr=b(i); k=i
          endif
       enddo! i

!       print*,b(1:n)

       if(myrank==1)then
          print*," "
          print*,"----- Check a calculation result. -----"
          print*,"max|Ax-wx|=",rr,k
          print*,"|A|=",ra
          print*,"epsilon=",m_epsilon
          print*,"max|Ax-wx|/Ne|A|=",rr/(n*m_epsilon*ra)
          print*,"|ZZ-I|=",SQRT(rz),rz
          print*,"---------------------------------------"
       endif
*
       return
       end
