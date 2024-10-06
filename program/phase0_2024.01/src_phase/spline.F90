
subroutine init_cubic_spline(n,x,y,y2)
  implicit none
  integer,intent(in) :: n
  real(kind=8), dimension(n), intent(in) :: x,y
  real(kind=8), dimension(n), intent(out) :: y2
  integer i,j
  real(kind=8) a,b
  real(kind=8), allocatable, dimension(:) :: tmp

  allocate(tmp(n));tmp=0.d0
  y2 = 0.0d0

  do i=2,n-1
    a=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    b=a*y2(i-1)+2
    y2(i)=(a-1)/b
    tmp(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
    & /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-a*tmp(i-1))/b
  enddo

  do j=n-1,1,-1
    y2(j) = y2(j)*y2(j+1)+tmp(j)
  enddo
  deallocate(tmp)

end subroutine init_cubic_spline

subroutine cubic_spline(n,x,y,y2,xin,yout,y1out)
  implicit none
  integer, intent(in) :: n
  real(kind=8), dimension(n), intent(in) :: x,y,y2
  real(kind=8), intent(in) :: xin
  real(kind=8), intent(out) :: yout,y1out
  real(kind=8) :: a,b,h
  integer :: i,i0,i1

  i0=1
  i1=n

  do while ( (i1-i0).gt.1 ) 
    i=(i1+i0)/2
    if (x(i).gt.xin) then
        i1=i
    else
        i0=i
    endif
  enddo

  h=x(i1)-x(i0)
  a=(x(i1)-xin)/h
  b=(xin-x(i0))/h
  yout  = a*y(i0)+b*y(i1)+((a*a*a-a)*y2(i0)+(b*b*b-b)*y2(i1))*(h**2)/6.0d0
  y1out = (y(i1)-y(i0))/h-(3.0d0*a*a-1)/6.0d0*h*y2(i0)+(3.0d0*b*b-1)/6.0d0*h*y2(i1)
end subroutine cubic_spline

function cubic_splinef(n,x,y,y2,xin,y1out) result (res)
  implicit none
  integer, intent(in) :: n
  real(kind=8), dimension(n), intent(in) :: x,y,y2
  real(kind=8), intent(in) :: xin
  real(kind=8), intent(out) :: y1out
  real(kind=8) :: a,b,h
  integer :: i,i0,i1
  real(kind=8) :: res

  i0=1
  i1=n

  do while ( (i1-i0).gt.1 ) 
    i=(i1+i0)/2
    if (x(i).gt.xin) then
        i1=i
    else
        i0=i
    endif
  enddo

  h=x(i1)-x(i0)
  a=(x(i1)-xin)/h
  b=(xin-x(i0))/h
  res  = a*y(i0)+b*y(i1)+((a*a*a-a)*y2(i0)+(b*b*b-b)*y2(i1))*(h**2)/6.0d0
  y1out = (y(i1)-y(i0))/h-(3.0d0*a*a-1)/6.0d0*h*y2(i0)+(3.0d0*b*b-1)/6.0d0*h*y2(i1)
end function cubic_splinef

subroutine init_bicubic_spline(x1a,x2a,ya,m,n,y2a)
  implicit none
  integer, intent(in) :: m,n
  real(kind(1.d0)), dimension(m), intent(in) :: x1a
  real(kind(1.d0)), dimension(n), intent(in) :: x2a
  real(kind(1.d0)), dimension(m,n), intent(in) :: ya
  real(kind(1.d0)), dimension(m,n), intent(out) :: y2a
  integer :: i,j,k
  real(kind(1.d0)), allocatable, dimension(:) :: y2tmp,ytmp

  allocate(y2tmp(n))
  allocate(ytmp(n))
  y2tmp=0.0d0;ytmp=0.0d0;y2a=0.0d0

  do j=1,m
    do k=1,n
      ytmp(k) = ya(j,k)
    enddo
    call init_cubic_spline(n,x2a,ytmp,y2tmp)
    do k=1,n
      y2a(j,k)=y2tmp(k)
    enddo
  enddo
  deallocate(y2tmp)
  deallocate(ytmp)
end subroutine init_bicubic_spline

subroutine init_tricubic_spline(x1a,x2a,x3a,ya,m,n,o,y2a)
  implicit none
  integer, intent(in) :: m,n,o
  real(kind(1.d0)), dimension(m), intent(in) :: x1a
  real(kind(1.d0)), dimension(n), intent(in) :: x2a
  real(kind(1.d0)), dimension(o), intent(in) :: x3a
  real(kind(1.d0)), dimension(m,n,o), intent(in) :: ya
  real(kind(1.d0)), dimension(m,n,o), intent(out) :: y2a

  integer :: i
  real(kind(1.d0)), allocatable, dimension(:,:) :: y2atmp,yatmp

  allocate(y2atmp(m,n))
  allocate(yatmp(m,n))
  y2atmp=0.0d0;yatmp=0.0d0;y2a=0.0d0

  do i=1,o
    yatmp(:,:) = ya(:,:,i)
    call init_bicubic_spline(x1a,x2a,yatmp,m,n,y2atmp)
    y2a(:,:,i) = y2atmp(:,:)
  enddo
  deallocate(y2atmp)
  deallocate(yatmp)
end subroutine init_tricubic_spline

real(kind(1.d0)) function bicubic_spline(x1a,x2a,ya,y2a,m,n,x1,x2,dy1,dy2)
  implicit none
  integer, intent(in) :: m,n
  real(kind(1.d0)), dimension(m), intent(in) :: x1a
  real(kind(1.d0)), dimension(n), intent(in) :: x2a
  real(kind(1.d0)), dimension(m,n), intent(in) :: ya,y2a
  real(kind(1.d0)), intent(in) :: x1,x2
  real(kind(1.d0)), intent(out) :: dy1,dy2

  integer j,k
  real(kind(1.d0)) :: dytmp,y1
  real(kind(1.d0)), allocatable, dimension(:) :: y2tmp, ytmp,yytmp,dytmpmat
  real(kind(1.d0)) :: cubic_splinef

  allocate(y2tmp(n))
  allocate(ytmp(n))
  allocate(yytmp(n))
  allocate(dytmpmat(n))
  y2tmp=0.0d0;ytmp=0.0d0;yytmp=0.0d0;dytmpmat=0.0d0

  do j=1,m
    ytmp(:) = ya(j,:)
    y2tmp(:) = y2a(j,:)
    yytmp(j) = cubic_splinef(n,x2a,ytmp,y2tmp,x2,dytmp)
    dytmpmat(j) = dytmp
  enddo

  dy1=0.0d0;dy2=0.0d0

!    y2tmp=0.0d0
  call init_cubic_spline(m,x1a,dytmpmat,y2tmp)
  dy2 = cubic_splinef(m,x1a,dytmpmat,y2tmp,x1,dytmp)

!    y2tmp=0.0d0
  call init_cubic_spline(m,x1a,yytmp,y2tmp)
  bicubic_spline = cubic_splinef(m,x1a,yytmp,y2tmp,x1,dy1)
  deallocate(y2tmp)
  deallocate(ytmp)
  deallocate(yytmp)
  deallocate(dytmpmat)
end function bicubic_spline

real(kind(1.d0)) function tricubic_spline(x1a,x2a,x3a,ya,y2a,m,n,o,x1,x2,x3,&
& dy1,dy2,dy3)
  implicit none
  integer, intent(in) :: m,n,o
  real(kind(1.d0)), dimension(m), intent(in) :: x1a
  real(kind(1.d0)), dimension(n), intent(in) :: x2a
  real(kind(1.d0)), dimension(o), intent(in) :: x3a
  real(kind(1.d0)), dimension(m,n,o), intent (in) :: ya
  real(kind(1.d0)), dimension(m,n,o), intent(in) :: y2a
  real(kind(1.d0)), intent(in) :: x1,x2,x3
  real(kind(1.d0)), intent(out) :: dy1,dy2,dy3

  integer j,k
  real(kind(1.d0)), allocatable, dimension(:) :: y2tmp
  real(kind(1.d0)), allocatable, dimension(:) :: yytmp
  real(kind(1.d0)), allocatable, dimension(:) :: dy1tmpmat
  real(kind(1.d0)), allocatable, dimension(:) :: dy2tmpmat
  real(kind(1.d0)) :: dy1tmp,dy2tmp,y3
  real(kind(1.d0)), allocatable, dimension(:,:) :: ytmp,y2atmp

  real(kind(1.d0)) :: bicubic_spline,cubic_splinef

  allocate(y2tmp(o))
  allocate(yytmp(o))
  allocate(dy1tmpmat(o))
  allocate(dy2tmpmat(o))
  allocate(ytmp(m,n))
  allocate(y2atmp(m,n))
  y2tmp=0.0d0;yytmp=0.0d0;dy1tmpmat=0.0d0
  dy2tmpmat=0.0d0;ytmp=0.0d0;y2atmp=0.0d0

  dy1=0.0d0;dy2=0.0d0;dy3=0.0d0
  do j=1,o
    ytmp(:,:) = ya(:,:,j)
    y2atmp(:,:) = y2a(:,:,j)
    yytmp(j) = bicubic_spline(x1a,x2a,ytmp,y2atmp,m,n,x1,x2,&
    & dy1tmp,dy2tmp)
    dy1tmpmat(j) = dy1tmp
    dy2tmpmat(j) = dy2tmp
  enddo

  call init_cubic_spline(o,x3a,dy1tmpmat,y2tmp)
  dy1 = cubic_splinef(o,x3a,dy1tmpmat,y2tmp,x3,dy1tmp)

  call init_cubic_spline(o,x3a,dy2tmpmat,y2tmp)
  dy2 = cubic_splinef(o,x3a,dy2tmpmat,y2tmp,x3,dy2tmp)

  call init_cubic_spline(o,x3a,yytmp,y2tmp)
  tricubic_spline = cubic_splinef(o,x3a,yytmp,y2tmp,x3,dy3)
  deallocate(y2tmp)
  deallocate(yytmp)
  deallocate(dy1tmpmat)
  deallocate(dy2tmpmat)
  deallocate(ytmp)
  deallocate(y2atmp)
end function tricubic_spline

