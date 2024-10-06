      subroutine eigen_hrd_v(
     &             ur_x,ui_x, vr_x,vi_x, vr_y,vi_y,
     &             u_t,v_t,
     &             beta, n,x_pos,nv
     &           )
!----
      use communication_h, only : reduce_dbl
     &               , datacast_dbl
!----
      implicit none

      integer, intent(in)    :: n, nv
      real(8), intent(in)    :: ur_x(1:nv)
      real(8), intent(inout) :: vr_x(1:nv), vr_y(1:nv)
      real(8), intent(in)    :: ui_x(1:nv)
      real(8), intent(inout) :: vi_x(1:nv), vi_y(1:nv)
      real(8), intent(inout) :: u_t(1:nv),  v_t(1:nv)

!---- real(8) ::  alpha, beta, gamma, theta
      real(8) ::  alpha(1), gamma(1), theta(1)
      real(8) ::  beta
      integer ::  j_1, j_2, j_3, x_pos

      real(8) ::  vr_x0, vi_x0
      real(8) ::  ur_x0, ui_x0

      include 'trd.h'

!----
      j_2 = 1
      j_3 = n
!----
! v':= v-((u,v)/2|u|^2)u
!----
      alpha(1) = 0.0d+00
      gamma(1) = 0.0d+00
      do j_1=j_2,j_3
         gamma(1) = gamma(1)+abs(vr_x(j_1))
     &                      +abs(vi_x(j_1))
      enddo! j_1
      call reduce_dbl(gamma, theta, 1, mpi_comm_col)
      if(gamma(1)/=0.0d+00)then
         do j_1=j_2,j_3
            vr_x0 = vr_x(j_1)
            vi_x0 = vi_x(j_1)
            vr_x0 = vr_x0/(gamma(1)*beta)
            vi_x0 = vi_x0/(gamma(1)*beta)
            ur_x0 = ur_x(j_1)
            ui_x0 = ui_x(j_1)
            alpha(1)     = alpha(1)+vr_x0*ur_x0
     &                             +vi_x0*ui_x0
            vr_x(j_1) = vr_x0
            vi_x(j_1) = vi_x0
         enddo! j_1
         alpha(1) = alpha(1)*gamma(1)
      endif
      call reduce_dbl(alpha, theta, 1, mpi_comm_col)
      alpha(1) = alpha(1)/(beta+beta)
      do j_1=j_2,j_3
         vr_x0 = vr_x(j_1)
         ur_x0 = ur_x(j_1)
         vi_x0 = vi_x(j_1)
         ui_x0 = ui_x(j_1)
         vr_x0 = gamma(1)*vr_x0-alpha(1)*ur_x0
         vi_x0 = gamma(1)*vi_x0-alpha(1)*ui_x0
         vr_x(j_1) = vr_x0
         vi_x(j_1) = vi_x0
      enddo! j_1
      call datacast_dbl(vr_y(1), vr_x(1), u_t(1), v_t(1), x_pos)
      call datacast_dbl(vi_y(1), vi_x(1), u_t(1), v_t(1), x_pos)

      return
      end subroutine

