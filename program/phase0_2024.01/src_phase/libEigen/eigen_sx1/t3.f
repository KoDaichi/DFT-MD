      subroutine tred1_v(u_x, v_x, v_y, nv, u_t, v_t, c, i)

      implicit NONE

      integer, intent(in)    :: nv, i
      real(8), intent(in)    :: u_x(1:nv,*)
      real(8), intent(inout) :: v_x(1:nv,*), v_y(1:nv,*)
      real(8), intent(out)   :: u_t(*), v_t(*)
      real(8), intent(inOut) :: c(2,2)

      real(8)                ::  tt(4), ss(4)
      real(8)                ::  s11,s12,s21,s22
      real(8)                ::  c11,c12,c21,c22
      real(8)                ::  s,t,u, u12
      real(8)                ::  sx(2,2), tx(2,2)
      real(8)                ::  u1, u2
      real(8)                ::  v1, v2

      integer                ::  x_pos, x_owner_nod
      integer                ::  y_pos, y_owner_nod
      integer                ::  j_1, j_2, j_3
      integer                ::  L
      integer                ::  n, LL, jj_1, jj_2, jj_3

!      integer, parameter     ::  VTOL = 512
      integer, parameter     ::  VTOL = 51200

      include 'trd.h'
      include 'param.h'

      L = i-2

      x_owner_nod = loop_c_node (L, x_nnod,x_inod)
      y_owner_nod = loop_c_node (i, y_nnod,y_inod)

      x_pos       = loop_c_g2l_depth(L, x_nnod,x_inod)
      y_pos       = loop_c_g2l_depth(i, y_nnod,y_inod)

      j_2         = loop_c_start(1, x_nnod,x_inod)
      j_3         = loop_c_end  (L, x_nnod,x_inod)

      n = j_3-j_2+1
!Kro  n = nv*2
      LL = (n-1)/y_nnod+1
      LL = ((LL-1)/2+1)*2

      if ( n > VTOL ) then
         jj_2 = j_2+(1+LL*(y_inod-1))-1
         jj_3 = j_2+(MIN(n,LL*y_inod))-1
      else
         jj_2 = j_2
         jj_3 = j_3
      end if

! S:=V^TU=U^TAU
      s11 = ZERO
      s12 = ZERO
      s22 = ZERO
      u12 = ZERO
!DIR$ UNROLL(2)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do j_1=jj_2,jj_3
         v1 = v_x(j_1,1)
         v2 = v_x(j_1,2)
         u1 = u_x(j_1,1)
         u2 = u_x(j_1,2)
         s11  = s11 + v1 * u1
         s12  = s12 + v1 * u2
         s22  = s22 + v2 * u2
!         u12  = u12 + u1 * u2
      end do! j_1
      ss(1) = s11
      ss(2) = s12
      ss(3) = s22
!      ss(4) = u12


      if ( n > VTOL ) then
         call reduce_dbl(ss(1), tt(1), 4, 1, y_COMM_WORLD)
      end if
      call reduce_dbl(ss(1), tt(1), 4, 1, x_COMM_WORLD)

!      u12 = ss(4)
      u12 = c(1,2)
      c(2,1) = -c(2,2)*(c(1,1)*u12)

!!       print*,"V2'=",v_x(1:L,2)

! SX:=S
      sx(1,1) = ss(1)
      sx(2,1) = ss(2)
      sx(1,2) = ss(2)
      sx(2,2) = ss(3)

! TX:=SX*C
      tx(1,1) = sx(1,1)*c(1,1) + sx(1,2)*c(2,1)
      tx(2,1) = sx(2,1)*c(1,1) + sx(2,2)*c(2,1)
      tx(1,2) =                  sx(1,2)*c(2,2)
      tx(2,2) =                  sx(2,2)*c(2,2)
            
! SX:=C^T*TX
      sx(1,1) = c(1,1)*tx(1,1) + c(2,1)*tx(2,1)
      sx(2,1) =                  c(2,2)*tx(2,1)
      sx(1,2) = c(1,1)*tx(1,2) + c(2,1)*tx(2,2)
      sx(2,2) =                  c(2,2)*tx(2,2)

      tx(1,1) = sx(1,1)
      tx(2,2) = sx(2,2)
      tx(2,1) = sx(1,2)+sx(2,1)
      tx(1,2) = ZERO

      s11 = -tx(1,1)/2
      s21 = -tx(2,1)/2
      s12 = -tx(1,2)/2
      s22 = -tx(2,2)/2

      c11 = c(1,1)
      c21 = c(2,1)
      c12 = ZERO
      c22 = c(2,2)

!!       print*,"S=[",s11,s12
!!       print*,"  [",s21,s22

! V:=VC^T
!DIR$ UNROLL(2)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      do j_1=jj_2,jj_3
         v1 = v_x(j_1,1)
         v2 = v_x(j_1,2)
         v1         =      v1 * c(1,1) + v2 * c(2,1)
         v2         =                    v2 * c(2,2)
!      end do! j_1
!
!!       print*,"V2''=",v_x(1:L,2)
!
! V:=V-US
!      do j_1=jj_2,jj_3
         u1 = u_x(j_1,1)
         u2 = u_x(j_1,2)
         v_x(j_1,1) = v1 + s11 * u1 + s21 * u2
         v_x(j_1,2) = v2            + s22 * u2
      end do! j_1

!       print*,"U2'''=",u_x(1:L,2)
!       print*,"U1'''=",u_x(1:L,1)
!       print*,"V2'''=",v_x(1:L,2)
!       print*,"V1'''=",v_x(1:L,1)

      if ( n > VTOL ) then
         do j_1=j_2,jj_2-1
            v_x(j_1,1)=ZERO
            v_x(j_1,2)=ZERO
         end do! j_1
         do j_1=jj_3+1,j_3
            v_x(j_1,1)=ZERO
            v_x(j_1,2)=ZERO
         end do! j_1

         if ( y_nnod > 1 ) then
         do j_1=j_2,j_3
            v_t(j_1)  =v_x(j_1,1)
            v_t(j_1+n)=v_x(j_1,2)
         end do
         call reduce_dbl(v_t(1), u_t(1), 2*n, 1, y_COMM_WORLD)
         do j_1=j_2,j_3
            v_x(j_1,1) = v_t(j_1)
            v_x(j_1,2) = v_t(j_1+n)
         end do
         end if
      end if

      call datacast_dbl(v_y(1,1), v_x(1,1), u_t(1), v_t(1), x_pos)
      call datacast_dbl(v_y(1,2), v_x(1,2), u_t(1), v_t(1), x_pos)

! FOR attention to unexpected overflow or NAN
       j_3 = loop_c_end      (L, x_nnod, x_inod)
       n   = loop_c_g2l_depth(L, x_nnod, x_inod)
       if ( j_3 < n ) then
          v_x(j_3+1:n,1) = ZERO ! in case
          v_x(j_3+1:n,2) = ZERO ! in case
       end if

       j_3 = loop_c_end      (L, y_nnod, y_inod)
       n   = loop_c_g2l_depth(L, y_nnod, y_inod)
       if ( j_3 < n ) then
          v_y(j_3+1:n,1) = ZERO ! in case
          v_y(j_3+1:n,2) = ZERO ! in case
       end if


      return
      end subroutine tred1_v

