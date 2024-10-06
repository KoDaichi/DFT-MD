#ifdef ACMLFFT
  subroutine dzfft3d(mode,l,m,n,x,ldl,ldm,ldn,comm,info)
    !! ldl >= l+2 (l+3) when l is even (odd), and ldl is even.
    !! ldm >= m
    !! ldn >= n
    !! comm : real*8 dimension(3*l+6*m+6*n+2*m*n+500)
    implicit none
    integer, intent(in) :: mode !! 0: initialization, 1: transform
    integer, intent(in) :: l,m,n
    integer, intent(in) :: ldl,ldm,ldn
    real(kind=8), intent(inout) :: x(ldl*ldm*ldn)
    real(kind=8), intent(inout), target :: comm(*)
    integer, intent(out) :: info

    integer :: ldh,incm,incn,i
    integer :: lcomm1, lcomm2
    real(kind=8), pointer :: comm1(:)
    real(kind=8), pointer :: comm2(:)

    lcomm1 = 3*l+100
    lcomm2 = 2*(m*n+3*m+3*n+200)

    comm1 => comm(1:lcomm1)
    comm2 => comm(lcomm1+1:lcomm1+lcomm2)

    ldh = ldl/2
    incm = ldh
    incn = ldh*ldm

    if(mode==0) then
       call dzfft(0,l,x,comm1,info)
       if(info/=0) return
       call zfft2dx(0,1.d0,.true.,.true.,m,n, &
                  & x,incm,incn,x,incm,incn, &
                  & comm2,info)
       return
    else if(mode/=1) then
       write(*,*) 'irregular mode = ',mode
       return
    end if

    call dzfftm2(m*n,l,x,ldl,comm1,info)
    if(info/=0) return
    do i=1,l/2+1
       call zfft2dx(-1,1.d0,.true.,.true.,m,n, &
                  & x(2*i-1),incm,incn,x,incm,incn, &
                  & comm2,info)
       if(info/=0) return
    end do
  end subroutine dzfft3d
 
  subroutine zdfft3d(mode,l,m,n,x,ldl,ldm,ldn,comm,info)
    !! ldl >= l+2 (l+3) when l is even (odd), and ldl is even.
    !! ldm >= m
    !! ldn >= n
    !! comm : real*8 dimension(3*l+6*m+6*n+2*m*n+500)
    implicit none
    integer, intent(in) :: mode !! 0: initialization, 1: transform
    integer, intent(in) :: l,m,n
    integer, intent(in) :: ldl,ldm,ldn
    real(kind=8), intent(inout) :: x(ldl*ldm*ldn)
    real(kind=8), intent(inout), target :: comm(*)
    integer, intent(out) :: info

    integer :: ldh,incm,incn,i
    integer :: lcomm1, lcomm2
    real(kind=8), pointer :: comm1(:)
    real(kind=8), pointer :: comm2(:)

    lcomm1 = 3*l+100
    lcomm2 = 2*(m*n+3*m+3*n+200)

    comm1 => comm(1:lcomm1)
    comm2 => comm(lcomm1+1:lcomm1+lcomm2)

    ldh = ldl/2
    incm = ldh
    incn = ldh*ldm

    if(mode==0) then
       call dzfft(0,l,x,comm1,info)
       if(info/=0) return
       call zfft2dx(0,1.d0,.true.,.true.,m,n, &
                  & x,incm,incn,x,incm,incn, &
                  & comm2,info)
       return
    else if(mode/=1) then
       write(*,*) 'irregular mode = ',mode
       return
    end if

    do i=1,l/2+1
       call zfft2dx(1,1.d0,.true.,.true.,m,n, &
                  & x(2*i-1),incm,incn,x,incm,incn, &
                  & comm2,info)
       if(info/=0) return
    end do
    call zdfftm2(m*n,l,x,ldl,comm1,info)
  end subroutine zdfft3d

  subroutine dzfftm2(m,n,x,ld,comm,info)
    !! ld >= n+2 when n is even.
    !! ld >= n+1 when n is odd.
    implicit none
    integer, intent(in) :: m,n,ld
    real(kind=8), intent(inout) :: x(ld*m)
    real(kind=8), intent(inout) :: comm(3*n+100)
    integer, intent(out) :: info

    real(kind=8) :: y(n)
    integer :: i,j
    real(kind=8) :: s

    s = sqrt(dble(n))

    do j=0,m-1
       do i=1,n
          y(i) = x(i+ld*j)
       end do
       call dzfftm(1,n,y,comm,info)
       do i=0,n/2
          x(1+2*i+ld*j) = s*y(1+i)
       end do
       do i=1,(n-1)/2
          x(2+2*i+ld*j) = s*y(1+n-i)
       end do
       x(2+ld*j) = 0.d0
       x(n+2+ld*j) = 0.d0
    end do
  end subroutine dzfftm2

  subroutine zdfftm2(m,n,x,ld,comm,info)
    !! ld >= n+2 when n is even.
    !! ld >= n+1 when n is odd.
    implicit none
    integer, intent(in) :: m,n,ld
    real(kind=8), intent(inout) :: x(ld*m)
    real(kind=8), intent(inout) :: comm(3*n+100)
    integer, intent(out) :: info

    real(kind=8) :: y(n)
    integer :: i,j
    real(kind=8) :: s

    s = sqrt(dble(n))

    do j=0,m-1
       do i=0,n/2
          y(1+i) = x(1+2*i+ld*j)
       end do
       do i=1,(n-1)/2
          y(1+n-i) = -x(2+2*i+ld*j)
       end do
       call zdfftm(1,n,y,comm,info)
       do i=1,n
          x(i+ld*j) = s*y(i)
       end do
       x(n+1+ld*j) = 0.d0
       x(n+2+ld*j) = 0.d0
    end do
  end subroutine zdfftm2
#else
  subroutine dummy_acmlfft
  end subroutine dummy_acmlfft
#endif
