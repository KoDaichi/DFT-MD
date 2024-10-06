!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 612 $)
!
!  SUBROUINE: real_spherical_harmonics
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
subroutine real_spherical_harmonics(is,n,x,y,z,ylm,coeff)
!$Id: real_spherical_harmonics.f90 612 2020-04-25 10:36:25Z ktagami $
  integer, intent(in) :: is,n
  real(kind=8), intent(in), dimension(n)  :: x(n),y(n),z(n)
  real(kind=8), intent(out) :: ylm(n)
  real(kind=8), intent(out) :: coeff

  real(kind=8) :: a,b,c,d,e,f,g,r
  real(kind=8) :: PAI4,PAI
 
  PAI  = 4.d0*atan(1.d0)
  PAI4 = 4.d0*PAI

  if(is == 1) then
     do i=1,n        
        a = sqrt(1.d0/PAI4)
        ylm(i) = a
     end do
  else if(is == 2) then
     do i=1,n
        r=sqrt(x(i)**2+y(i)**2+z(i)**2)
        a = sqrt(3.d0/PAI4)
        ylm(i) = a*x(i)/r
     end do 
  else if(is == 3) then
     do i=1,n
        r=sqrt(x(i)**2+y(i)**2+z(i)**2)
        a = sqrt(3.d0/PAI4)
        ylm(i) = a*y(i)/r
     end do
  else if(is == 4) then
     do i=1,n
        r=sqrt(x(i)**2+y(i)**2+z(i)**2)
        a = dsqrt(3.d0/PAI4)
        ylm(i) = a*z(i)/r
     end do
  else if(is == 5) then
     do i=1,n
        a = dsqrt(5.d0/(16*PAI))
        d = z(i)**2
        e = x(i)**2+y(i)**2+d
        ylm(i) = a*(3*d-e)/e
     end do
  else if(is == 6) then
     do i=1,n
        a = dsqrt(15.d0/(16*PAI))
        b = x(i)**2
        c = y(i)**2
        e = b+c+z(i)**2
        ylm(i) = a*(b-c)/e
     end do
  else if(is == 7) then
     do i=1,n
        a = dsqrt(15.d0/PAI4)
        e = x(i)**2+y(i)**2+z(i)**2
        ylm(i) = a*x(i)*y(i)/e
     end do
  else if(is == 8) then
     do i=1,n
        a = dsqrt(15.d0/PAI4)
        e = x(i)**2+y(i)**2+z(i)**2
        ylm(i) = a*y(i)*z(i)/e
     end do
  else if(is == 9) then
     do i=1,n
        a = dsqrt(15.d0/PAI4)
        e = x(i)**2+y(i)**2+z(i)**2
        ylm(i) = a*z(i)*x(i)/e
     end do
  else if(is == 10) then
     do i=1,n
        a = dsqrt(7.d0/(16*PAI))
        d = z(i)**2
        e = x(i)**2+y(i)**2+d
        r=sqrt(e)
        f = e * r
        ylm(i) = a*z(i)*(5*d-3*e)/f
     end do
  else if(is == 11) then
     do i=1,n
        a = dsqrt(21.d0/(32*PAI))
        d = z(i)**2
        e = x(i)**2+y(i)**2+d
        r=sqrt(e)
        f = e * r
        ylm(i) = a*x(i)*(5*d-e)/f
     end do
  else if(is == 12) then
     do i=1,n
        a = dsqrt(21.d0/(32*PAI))
        d = z(i)**2
        e = x(i)**2+y(i)**2+d
        r=sqrt(e)
        f = e * r
        ylm(i) = a*y(i)*(5*d-e)/f
     end do
  else if(is == 13) then
     do i=1,n
        a = dsqrt(105.d0/(16*PAI))
        b = x(i)*x(i)
        c = y(i)*y(i)
        e = b+c+z(i)*z(i)
        r=sqrt(e)
        f = e * r
        ylm(i) = a*z(i)*(b-c)/f
     end do
  else if(is == 14) then
     do i=1,n
        r=sqrt(x(i)**2+y(i)**2+z(i)**2)
        a = dsqrt(105.d0/PAI4)
        f = r**3
        ylm(i) = a*x(i)*y(i)*z(i)/f
     end do
  else if(is == 15) then
     do i=1,n
        a = dsqrt(35.d0/(32*PAI))
        b = x(i)*x(i)
        c = y(i)*y(i)
        r=sqrt(b+c+z(i)*z(i))
        f = r**3
        ylm(i) = a*x(i)*(b-3*c)/f
     end do
  else if(is == 16) then
     do i=1,n
        a = dsqrt(35.d0/(32*pai))
        b = x(i)*x(i)
        c = y(i)*y(i)
        r=sqrt(b+c+z(i)*z(i))
        f = r**3
        ylm(i) = a*y(i)*(3*b-c)/f
     end do
  else if(is == 17) then
     do i=1,n
        a = 3.d0/8.d0/dsqrt(PAI4)
        d = z(i)*z(i)
        e = x(i)*x(i)+y(i)*y(i)+d
        r=sqrt(e)
        f = e*e
        ylm(i) = a*(5*d*(7*d-6*e)/f+3.d0)
     end do
  else if(is == 18) then
     do i=1,n
        a = 15.d0/4.d0/dsqrt(10*PAI)
        d = z(i)*z(i)
        e = x(i)*x(i)+y(i)*y(i)+d
        r=sqrt(e)
        f = e*e
        ylm(i) = a*z(i)*x(i)*(7*d-3*e)/f
     end do
  else if(is == 19) then
     do i=1,n
        a = 15.d0/4.d0/dsqrt(10.d0*PAI)
        d = z(i)*z(i)
        e = x(i)*x(i)+y(i)*y(i)+d
        r=sqrt(e)
        f = e*e
        ylm(i)= a*y(i)*z(i)*(7*d-3*e)/f
     end do
  else if(is == 20) then
     do i=1,n
        a = 15.d0/8.d0/dsqrt(5*PAI)
        b = x(i)*x(i)
        c = y(i)*y(i)
        d = z(i)*z(i)
        e = b+c+d
        f = e*e
        ylm(i) = a*(7*d-e)*(b-c)/f
     end do
  else if(is == 21) then
     do i=1,n
        a = 15.d0/4.d0/dsqrt(5*PAI)
        d = z(i)*z(i)
        e = x(i)*x(i)+y(i)*y(i)+d
        r=sqrt(e)
        f = e*e
        ylm(i) = a*(7*d-e)*x(i)*y(i)/f
     end do
  else if(is == 22) then
     do i=1,n
        a = 105.d0/4.d0/dsqrt(70*pai)
        b = x(i)*x(i)
        c = y(i)*y(i)
        f = (b+c+z(i)*z(i))**2
        ylm(i) = a*(b-3*c)*z(i)*x(i)/f
     end do
  else if(is == 23) then
     do i=1,n 
        a = 105.d0/4.d0/dsqrt(70*PAI)
        b = x(i)*x(i)
        c = y(i)*y(i)
        f = (b+c+z(i)*z(i))**2
        ylm(i) = a*(3.d0*b-c)*y(i)*z(i)/f
     end do
  else if(is == 24) then
     do i=1,n 
        a = 105.d0/16.d0/dsqrt(35*PAI)
        b = x(i)*x(i)
        c = y(i)*y(i)
        f = (b+c+z(i)*z(i))**2
        ylm(i) = a*((b-c)**2-4*b*c)/f
     end do
  else if(is == 25) then
     do i=1,n 
        a = 105.d0/4.d0/dsqrt(35*PAI)
        b = x(i)*x(i)
        c = y(i)*y(i)
        f = (b+c+z(i)*z(i))**2
        ylm(i) = a*(b-c)*x(i)*y(i)/f
     end do
  else if(is == 26) then
     do i=1,n 
        a = dsqrt(11.d0/PAI)/16.d0
        b = z(i)*z(i)
        c = x(i)*x(i)+y(i)*y(i)+b
        r = sqrt(c)
        f = r**5
        ylm(i) = a*z(i)*(63.d0*b*b-70.d0*b*c+15.d0*c*c)/f
     end do
  else if(is == 27) then
     do i=1,n 
        a = dsqrt(165.d0/PAI)/16.d0
        b = z(i)*z(i)
        c = x(i)*x(i)+y(i)*y(i)+b
        r = sqrt(c)
        f = r**5
        ylm(i) = a*x(i)*(21.d0*b*b-14.d0*b*c+c*c)/f
     end do
  else if(is == 28) then
     do i=1,n 
        a = dsqrt(165.d0/PAI)/16.d0
        b = z(i)*z(i)
        c = x(i)*x(i)+y(i)*y(i)+b
        r = sqrt(c)
        f = r**5
        ylm(i) = a*y(i)*(21.d0*b*b-14.d0*b*c+c*c)/f
     end do
  else if(is == 29) then
     do i=1,n 
        a = dsqrt(1155.d0/(64*PAI))
        b = z(i)*z(i)
        c = x(i)*x(i)
        d = y(i)*y(i)
        e = c+d
        r = sqrt(b+e)
        f = r**5
        ylm(i) = a*z(i)*(2.d0*b-e)*(c-d)/f
     end do
  else if(is == 30) then
     do i=1,n 
        a = dsqrt(1155.d0/(16.d0*PAI))
        b = z(i)*z(i)
        c = x(i)*x(i)
        d = y(i)*y(i)
        e = c+d
        r = sqrt(b+e)
        f = r**5
        ylm(i) = a*x(i)*y(i)*z(i)*(2.d0*b-e)/f
     end do
  else if(is == 31) then
     do i=1,n 
        a = dsqrt(385.d0/(512.d0*PAI))
        b = z(i)*z(i)
        c = x(i)*x(i)
        d = y(i)*y(i)
        e = c+d
        r = sqrt(b+e)
        f = r**5
        ylm(i) = a*x(i)*(8.d0*b-e)*(c-3.d0*d)/f
     end do
  else if(is == 32) then
     do i=1,n 
        a = dsqrt(385.d0/(512.d0*PAI))
        b = z(i)*z(i)
        c = x(i)*x(i)
        d = y(i)*y(i)
        e = c+d
        r = sqrt(b+e)
        f = r**5
        ylm(i) = a*y(i)*(8.d0*b-e)*(3.d0*c-d)/f
     end do
  else if(is == 33) then
     do i=1,n 
        a = dsqrt(385.d0/(256.d0*PAI))*3.d0
        b = x(i)*x(i)
        c = y(i)*y(i)
        r = sqrt(b+c+z(i)*z(i))
        f = r**5
        ylm(i) = a*z(i)*(b*b-6.d0*b*c+c*c)/f
     end do
  else if(is == 34) then
     do i=1,n 
        r=sqrt(x(i)**2+y(i)**2+z(i)**2)
        a = dsqrt(385.d0/(16.d0*PAI))*3.d0
        b = x(i)*x(i)
        c = y(i)*y(i)
        d = b-c
        r=sqrt(b+c+z(i)*z(i))
        f = r**5
        ylm(i) = a*x(i)*y(i)*z(i)*d/f
     end do
  else if(is == 35) then
     do i=1,n 
        a = dsqrt(693.d0/(512.d0*PAI))
        b = x(i)*x(i)
        c = y(i)*y(i)
        r=sqrt(b+c+z(i)**2)
        f = r**5
        ylm(i) = a*x(i)*(b*b-10.d0*b*c+5.d0*c*c)/f
     end do
  else if(is == 36) then
     do i=1,n 
        a = dsqrt(693.d0/(512.d0*PAI))
        b = y(i)*y(i)
        c = x(i)*x(i)
        r=sqrt(b+c+z(i)**2)
        f = r**5
        ylm(i) = a*y(i)*(b*b-10.d0*b*c+5.d0*c*c)/f
     end do
  else if(is == 37) then
     do i=1,n 
        a = dsqrt(13.d0/(1024.d0*PAI))
        b = z(i)*z(i)
        c = b*b
        d = c*b
        e = x(i)*x(i)+y(i)*y(i)+b
        c = c*e
        f = e*e
        b = b*f
        f = f*e
        ylm(i) = a*((231.d0*d-315.d0*c+105.d0*b)/f-5.d0)
     end do
  else if(is == 38) then
     do i=1,n 
        a = dsqrt(273.d0/(256.d0*PAI))
        b = z(i)*z(i)
        c = b*b
        e = x(i)*x(i)+y(i)*y(i)+b
        b = b*e
        d = e*e
        f = d*e
        ylm(i) = a*z(i)*x(i)*(33.d0*c-30.d0*b+5.d0*d)/f
     end do
  else if(is == 39) then
     do i=1,n 
        a = dsqrt(273.d0/(256.d0*PAI))
        b = z(i)*z(i)
        c = b*b
        e = x(i)*x(i)+y(i)*y(i)+b
        b = b*e
        d = e*e
        f = d*e
        ylm(i) = a*y(i)*z(i)*(33.d0*c-30.d0*b+5.d0*d)/f
     end do
  else if(is == 40) then
     do i=1,n 
        a = dsqrt(1365.d0/(2048.d0*PAI))
        b = z(i)*z(i)
        c = b*b
        d = x(i)*x(i)
        f = y(i)*y(i)
        e = d+f+b 
        d = d-f
        b = b*e
        g = e*e
        f = g*e
        ylm(i) = a*d*(33.d0*c-18.d0*b+g)/f
     end do
  else if(is == 41) then
     do i=1,n 
        a = dsqrt(1365.d0/(512.d0*PAI))
        b = z(i)*z(i)
        c = b*b
        e = x(i)*x(i)+y(i)*y(i)+b 
        b = b*e
        g = e*e
        f = g*e
        ylm(i) = a*x(i)*y(i)*(33.d0*c-18.d0*b+g)/f
     end do
  else if(is == 42) then
     do i=1,n 
        a = dsqrt(1365.d0/(512.d0*PAI))
        b = z(i)*z(i)
        c = x(i)*x(i)
        d = y(i)*y(i)
        e = c+d+b 
        f = e**3
        ylm(i) = a*z(i)*x(i)*(11.d0*b-3.d0*e)*(c-3.d0*d)/f
     end do
  else if(is == 43) then
     do i=1,n 
        a = dsqrt(1365.d0/(512.d0*PAI))
        b = z(i)*z(i)
        c = y(i)*y(i)
        d = x(i)*x(i)
        e = d+c+b 
        f = e**3
        ylm(i) = a*y(i)*z(i)*(11.d0*b-3.d0*e)*(3.d0*d-c)/f
     end do
  else if(is == 44) then
     do i=1,n 
        a = dsqrt(819.d0/(1024d0*PAI))
        b = x(i)*x(i)
        c = y(i)*y(i)
        d = z(i)*z(i)
        e = b+c+d
        f = e**3
        ylm(i) = a*(11.d0*d-e)*(b*b-6.d0*b*c+c*c)/f
     end do
  else if(is == 45) then
     do i=1,n 
        a = dsqrt(819.d0/(64.d0*PAI))
        b = x(i)*x(i)
        c = y(i)*y(i)
        d = z(i)*z(i)
        e = b+c+d
        f = e**3
        ylm(i) = a*x(i)*y(i)*(11.d0*d-e)*(b-c)/f
     end do
  else if(is == 46) then
     do i=1,n 
        a = dsqrt(9009.d0/(512.d0*PAI))
        b = x(i)*x(i)
        c = y(i)*y(i)
        d = z(i)*z(i)
        e = b+c+d
        f = e**3
        ylm(i) = a*z(i)*x(i)*(b*b-10.d0*b*c+5.d0*c*c)/f
     end do
  else if(is == 47) then
     do i=1,n 
        a = dsqrt(9009.d0/(512.d0*PAI))
        b = y(i)*y(i)
        c = x(i)*x(i)
        e = c+b+z(i)*z(i)
        f = e**3
        ylm(i) = a*y(i)*z(i)*(b*b-10.d0*b*c+5.d0*c*c)/f
     end do
  else if(is == 48) then
     do i=1,n 
        a = dsqrt(3003.d0/(2048.d0*PAI))
        b = x(i)*x(i)
        c = y(i)*y(i)
        e = b+c+z(i)*z(i)
        f = e**3
        ylm(i) = a*(b-c)*(b*b-14.d0*b*c+c*c)/f
     end do
  else if(is == 49) then
     do i=1,n 
        a = dsqrt(3003.d0/(2048.d0*PAI))
        b = x(i)*x(i)
        c = y(i)*y(i)
        d = b+c
        e = b+c+z(i)*z(i)
        f = e**3
        ylm(i) = a*x(i)*y(i)*(6.d0*d*d-32.d0*b*c)/f
     end do
  end if
  coeff = a
end subroutine real_spherical_harmonics
