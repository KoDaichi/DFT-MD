       subroutine eigen_hrd_u(
     &              ar,ai,
     &              ur_x,ui_x,
     &              ur_y,ui_y,
     &              u_t,v_t,
     &              x_root,y_root,x_pos,
     &              n,nm,nv
     &            )
!
      use communication_h, only : send_dbl
     &               , recv_dbl
     &               , reduce_dbl
     &               , bcast_dbl
     &               , datacast_dbl
!
      implicit none

      integer, intent(in)    :: nm, nv, n
      real(8), intent(in)    :: ar(1:nm)
      real(8), intent(in)    :: ai(1:nm)
      real(8), intent(inout) :: ur_x(1:nv), ur_y(1:nv)
      real(8), intent(inout) :: ui_x(1:nv), ui_y(1:nv)
      real(8), intent(inout) :: u_t(1:nv+4), v_t(1:nv+4)

      integer ::  j_1, j_2, j_3
      integer ::  x_root, x_pos
      integer ::  y_root

!...  real(8) ::  scale, theta
      real(8) ::  scale(1), theta(1)
      real(8) ::  alpha, beta, gamma(2), kappa(2)
      real(8) ::  wk1(3), wk2(3)

      include 'trd.h'


      j_2 = 1
      j_3 = n
!
! u=...
!
!     print*,"eigen_hrd_u start"
      if(y_root==my_row)then
         scale(1) = 0.0d+00
         do j_1=j_2,j_3
            scale(1) = scale(1)+abs(ar(j_1))+abs(ai(j_1))
         enddo! j_1
         call reduce_dbl(scale, theta, 1, mpi_comm_col)
         if(scale(1)/=0.0d+00)then
            alpha    = 0.0d+00
            gamma(1) = 0.0d+00
            gamma(2) = 0.0d+00
            do j_1=j_2,j_3
               u_t(j_1) = ar(j_1)/scale(1)
               v_t(j_1) = ai(j_1)/scale(1)
               alpha     = alpha+u_t(j_1)**2+v_t(j_1)**2
            enddo! j_1
            if(x_root==my_col)then
               gamma(1) = u_t(x_pos)
               gamma(2) = v_t(x_pos)
            endif
            wk1(1) = gamma(1)
            wk1(2) = gamma(2)
            wk1(3) = alpha
            call reduce_dbl(wk1, wk2, 3, mpi_comm_col)
            gamma(1) = wk1(1)
            gamma(2) = wk1(2)
            alpha    = wk1(3)
         endif
         u_t(x_pos+1) = scale(1)
         u_t(x_pos+2) = alpha
         u_t(x_pos+3) = gamma(1)
         v_t(x_pos+3) = gamma(2)
      endif
!
      if ( nprocs >= es_cluster**2 .and.
     &       size_of_col > 1 .and. es_icnod > 1 ) then
         call recv_dbl(u_t(1), 0, es_icnod-1, es_ccomm_world)
      endif
      call bcast_dbl(u_t(1), x_pos+3, y_root, mpi_comm_row)
      call bcast_dbl(v_t(1), x_pos+3, y_root, mpi_comm_row)
      if ( nprocs >= es_cluster**2 .and.
     &       size_of_col > 1 .and. es_icnod < es_ncnod ) then
         call send_dbl(u_t(1), 0, es_icnod+1, es_ccomm_world)
      endif
!
      scale(1)=u_t(x_pos+1)
      if(scale(1)/=0.0d+00)then
         alpha = u_t(x_pos+2)
         gamma(1) = u_t(x_pos+3)
         gamma(2) = v_t(x_pos+3)
         theta(1) = sqrt(alpha)
         kappa(1) = sqrt(u_t(x_pos+3)**2+v_t(x_pos+3)**2)
         ur_x(j_2:j_3) = u_t(j_2:j_3)
         ui_x(j_2:j_3) = v_t(j_2:j_3)
         if(x_root==my_col)then
            if(kappa(1)/=0.0d+00)then
               kappa(2) = 1.0d+00+theta(1)/kappa(1)
            else
               kappa(2) = 1.0d+00
            endif
            ur_x(x_pos) = kappa(2)*ur_x(x_pos)
            ui_x(x_pos) = kappa(2)*ui_x(x_pos)
         endif
         call datacast_dbl(ur_y,ur_x, u_t,v_t, x_pos)
         call datacast_dbl(ui_y,ui_x, u_t,v_t, x_pos)
         beta = alpha+kappa(1)*theta(1)
         if(x_root==my_col)then
            u_t(x_pos  ) = ur_x(x_pos  )
            v_t(x_pos  ) = ui_x(x_pos  )
         endif
         v_t(x_pos+1) = gamma(1)
         v_t(x_pos+2) = gamma(2)
         v_t(x_pos+3) = kappa(1)
         v_t(x_pos+4) = kappa(2)
         u_t(x_pos+1) = scale(1)
         u_t(x_pos+2) = alpha
         u_t(x_pos+3) = theta(1)
         u_t(x_pos+4) = beta
      endif
!     print*,"eigen_hrd_u passed"

      return
      end subroutine

