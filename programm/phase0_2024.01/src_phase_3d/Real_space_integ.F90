!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: Gauss_Legendre
!
!  AUTHOR: Takenori Yamamotoi   June/01/2004
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
!=======================================================================
!
! coded by Takenori YAMAMOTO (Univ. Tokyo) at Jun. 1, 2004
!
subroutine Gauss_Legendre(x,w,n)
! $Id: Real_space_integ.f90 570 2017-04-21 20:34:50Z yamasaki $
implicit none
!
! int root(1-x**2) * f(x) dex = sum_{i} w_i * f(x_i)
!
integer, intent(in) :: n
real(kind=8), intent(out) :: x(n),w(n)

integer, parameter :: itermax = 100
real(kind=8), parameter :: eps = 3.d-14

integer :: i,j,iter,m
real(kind=8) :: z,z1,p1,p2,p3,pp
logical :: conv

m=(n+1)/2

!!! write(*,'("m=",i3)') m

do i=1,m
   z = cos(4.d0*atan(1.d0)*(i-0.25d0)/(n+0.5d0))
   do iter=1,itermax
      p1=1.d0
      p2=0.d0
      do j=1,n
         p3=p2
         p2=p1
         p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
      end do
      pp=n*(z*p1-p2)/(z*z-1.d0)
      z1 = z
      z = z1 - p1/pp
      ! debug
      !  write(*,*) 'i=',i,'iter=',iter,' z=',z
      ! end debug
      if(abs(z-z1) <= eps) then
         conv = .true.
         exit
      end if
   end do
   if(.not.conv) write(6,*) 'too many iterations in Gauss_Legendre'
   x(i) = -z
   x(n+1-i) = z
   w(i) = 2.d0/((1.d0-z*z)*pp*pp)
   w(n+1-i) = w(i)
end do
      
end subroutine Gauss_Legendre

!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: int_spherical_surf
!
!  AUTHOR: Takenori Yamamotoi   June/01/2004
!
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)
!
!=======================================================================
!
! coded by Takenori YAMAMOTO (Univ. Tokyo) at Jun. 1, 2004
!
subroutine int_spherical_surf(lmax,nmax,n,x,y,z,w)
! $Id: Real_space_integ.f90 570 2017-04-21 20:34:50Z yamasaki $
!
!  Int F(x,y,z) dOmega = sum_i w(i) * F(x(i),y(i),z(i))
!  ( x^2 + y^2 + z^2 = 1 )
!
!  nt = 1 + int(3*lmax/2)
!  np = 1 + 3*lmax
!  n = np*nt
!
implicit none
integer, intent(in) :: lmax,nmax
integer, intent(out) :: n
real(kind=8), intent(out) :: x(nmax),y(nmax),z(nmax),w(nmax)

integer :: i,j,m,nt,np
real(kind=8), allocatable :: xt(:),wt(:)
real(kind=8) :: pi,wp,dt,dp,p

nt = 1+int(3*lmax/2)
np = 1+3*lmax
n  = nt*np

if(nmax < n) call phase_error_with_msg(6,'nmax is too samll.',__LINE__,__FILE__)

allocate(xt(nt),wt(nt))

call Gauss_Legendre(xt,wt,nt)
! debug
!write(*,'("debug: nt=",i3))') nt
!do i=1,nt
!  write(*,'("debug:",2(1x,f10.5))') xt(i),wt(i)
!end do
! end debug

pi = 4.d0*atan(1.d0)
wp = 2.d0*pi/dble(np)
dp = 2.d0*pi/dble(np)

do i=1,nt
   do j=1,np
      m = (j-1)*nt + (i-1) + 1
      p = dp*(j-1)
      x(m)=sin(acos(xt(i)))*cos(p)
      y(m)=sin(acos(xt(i)))*sin(p)
      z(m)=xt(i)
      w(m)=wt(i)*wp
   end do
end do

deallocate(xt)
deallocate(wt)

end subroutine int_spherical_surf
