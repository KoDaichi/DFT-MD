!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 612 $)
!
!  SUBROUINE: get_crotylm
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
subroutine get_crotylm(l1max,mmax,nsph,nopr,crotylm,iylm,nylm,op)
! $Id: crotylm.f90 612 2020-04-25 10:36:25Z ktagami $
  implicit none
  integer, intent(in) :: l1max,mmax,nsph,nopr
  integer, intent(out) :: iylm(mmax,nsph,nopr),nylm(nsph,nopr)
  real(kind=8), intent(out) :: crotylm(mmax,nsph,nopr)
  real(kind=8), intent(in) :: op(3,3,nopr)

  integer :: llmax,lmax,nmax,n,np,nt
  integer :: m,m2,is,is1,is2,i,l,iopr
  integer :: mm
  real(kind=8) :: dint, cwk
  real(kind=8),allocatable :: x(:),y(:),z(:),w(:)
  real(kind=8),allocatable :: x2(:),y2(:),z2(:)
  real(kind=8),allocatable :: ylm(:,:),ylmr(:)
  real(kind=8), parameter  :: eps=1.d-6

  lmax = l1max-1
  nt = 1 + int(3*lmax/2)
  np = 1 + 3*lmax
  nmax = np*nt
  allocate(x(nmax),y(nmax),z(nmax),w(nmax))
  call int_spherical_surf(lmax,nmax,n,x,y,z,w)

  lmax = l1max-1
  llmax=(lmax+1)**2
  allocate(ylm(n,llmax),ylmr(n))
  do is=1,llmax
     call real_spherical_harmonics(is,n,x,y,z,ylm(1,is),cwk)
! debug
!     write(*,'(i3,10(1x,f25.10))') is,ylm(1:10,is)
! end debug
  end do

  allocate(x2(nmax),y2(nmax),z2(nmax))
  do iopr=1,nopr
     do i=1,n
        x2(i) = op(1,1,iopr)*x(i)+op(1,2,iopr)*y(i)+op(1,3,iopr)*z(i)
        y2(i) = op(2,1,iopr)*x(i)+op(2,2,iopr)*y(i)+op(2,3,iopr)*z(i)
        z2(i) = op(3,1,iopr)*x(i)+op(3,2,iopr)*y(i)+op(3,3,iopr)*z(i)
     end do
     do l=0,lmax
        do m2=1,2*l+1
           is2 = l**2+m2
           call real_spherical_harmonics(is2,n,x2,y2,z2,ylmr,cwk)
! debug
!           write(*,'(4(1x,i3))') iopr,l,m2,is2
! end debug
           mm=0
           do m=1,2*l+1
              is1 = l**2+m
              dint = 0.d0
              do i=1,n
                 dint = dint + ylm(i,is1)*ylmr(i)*w(i)
              end do
              if(dabs(dint)>eps) then
                 mm = mm+1
                 nylm(is2,iopr) = mm
                 iylm(mm,is2,iopr) = is1
                 crotylm(mm,is2,iopr) = dint
! debug
!                 write(*,'(1x,i3,"=",f25.10)') is1,dint
! end debug
              end if
           end do
        end do
     end do
  end do

  deallocate(x,y,z,w)
  deallocate(x2,y2,z2)
  deallocate(ylm,ylmr)

end subroutine get_crotylm
