   subroutine set_weight_unif(ier,n1,n2,x,w)
!=====================================================================
!
!  M. Okamoto
!
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2
   real(8),intent(in)  :: x(*)
   integer,intent(out) :: ier
   real(8),intent(out) :: w(*)
   real(8) :: h
   integer :: n, i
   ier = 0
   n = n2-n1+1
   if (n < 6) then
      ier = 1 ; go to 99
   end if
   h = (x(n2)-x(n1))/dble(n2-n1)
   w(n1  ) =  3.d0/8.d0  * h
   w(n1+1) =  7.d0/6.d0  * h
   w(n1+2) = 23.d0/24.d0 * h
   w(n2-2) = 23.d0/24.d0 * h
   w(n2-1) =  7.d0/6.d0  * h
   w(n2  ) =  3.d0/8.d0  * h
   if (n > 6) then
      do i = n1+3,n2-3
         w(i) = h
      end do
   end if
99 continue
   end subroutine set_weight_unif

!=====================================================================
   subroutine set_weight_exp(ier,n1,n2,r,w)
!=====================================================================
!
!  M. Okamoto
!
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2
   real(8),intent(in)  :: r(*)
   integer,intent(out) :: ier
   real(8),intent(out) :: w(*)
   real(8) :: h
   integer :: n, i
   ier = 0
   n = n2-n1+1
   if (n < 6) then
      ier = 1 ; go to 99
   end if
   h = log(r(n2)/r(n1)) / dble(n2-n1)
   w(n1  ) =  3.d0/8.d0  * r(n1  ) * h
   w(n1+1) =  7.d0/6.d0  * r(n1+1) * h
   w(n1+2) = 23.d0/24.d0 * r(n1+2) * h
   w(n2-2) = 23.d0/24.d0 * r(n2-2) * h
   w(n2-1) =  7.d0/6.d0  * r(n2-1) * h
   w(n2  ) =  3.d0/8.d0  * r(n2  ) * h
   if (n > 6) then
      do i = n1+3,n2-3
         w(i) = r(i) * h
      end do
   end if
99 continue
   end subroutine set_weight_exp
! === For nrc decomposion. by takto 2012/12/05 =================================
   subroutine set_weight_exp_3D(ier,n1,n2,r,w,ista_nrc,iend_nrc,ist,ien)
!=====================================================================
!
!  M. Okamoto
!
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2
   real(8),intent(in)  :: r(n2)
   integer,intent(out) :: ier
   real(8),intent(out) :: w(ista_nrc:iend_nrc)
   integer,intent(in)  :: ista_nrc,iend_nrc,ist,ien
   real(8) :: h
   integer :: n, i
   ier = 0
   n = n2-n1+1
   if (n < 6) then
      ier = 1 ; go to 99
   end if
   h = log(r(n2)/r(n1)) / dble(n2-n1)
   if(n1 >= ista_nrc .and. n1 <= iend_nrc) then
      w(n1  ) =  3.d0/8.d0  * r(n1  ) * h
   end if
   if(n1+1 >= ista_nrc .and. n1+1 <= iend_nrc) then
      w(n1+1) =  7.d0/6.d0  * r(n1+1) * h
   end if
   if(n1+2 >= ista_nrc .and. n1+2 <= iend_nrc) then
      w(n1+2) = 23.d0/24.d0 * r(n1+2) * h
   end if
   if(n2-2 >= ista_nrc .and. n2-2 <= iend_nrc) then
      w(n2-2) = 23.d0/24.d0 * r(n2-2) * h
   end if
   if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
      w(n2-1) =  7.d0/6.d0  * r(n2-1) * h
   end if
   if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
      w(n2  ) =  3.d0/8.d0  * r(n2  ) * h
   end if
   if (n > 6) then
      do i = max(n1+3,ist),min(n2-3,ien)
         w(i) = r(i) * h
      end do
   end if
99 continue
   end subroutine set_weight_exp_3D
! ==============================================================================
!=====================================================================
   subroutine set_open_weight_exp(ier,i0,is,r,w)
!=====================================================================
!
!  M. Okamoto
!
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: i0, is
   real(8),intent(in)  :: r(*)
   integer,intent(out) :: ier
   real(8),intent(out) :: w(*)
   real(8) :: h
   ier = 0
   h = log(r(i0+4*is)/r(i0)) / dble(4)
   if(is.lt.0) h=-h
   w(i0+  is) = 55.d0/24.d0 * r(i0+  is) * h
   w(i0+2*is) =-59.d0/24.d0 * r(i0+2*is) * h
   w(i0+3*is) = 37.d0/24.d0 * r(i0+3*is) * h
   w(i0+4*is) = -9.d0/24.d0 * r(i0+4*is) * h
99 continue
   end subroutine set_open_weight_exp

!=====================================================================
   subroutine set_weight_exp3(ier,n1,n2,dn,r,z,w)
!=====================================================================
!
!  M. Okamoto & AAS
!
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: r(*)
   real(8),intent(in)  :: z
   integer,intent(out) :: ier
   real(8),intent(out) :: w(*)
   real(8) :: h, a1,a2,a3,a4,tw3tw4,svsx,theght
   integer :: n, i
   ier = 0
   n = n2-n1+1
   if (n < 9*dn) then
      ier = 1 ; go to 99
   end if
   h = log(r(n2)/r(n1)) / dble(n2-n1) * dble(dn)
   a1           = -z*(z-1.d0)*(z-2.d0)/6.d0
   a2           = (z+1.d0)*(z-1.d0)*(z-2.d0)/2.d0
   a3           = -(z+1.d0)*z*(z-2.d0)/2.d0 
   a4           = (z+1.d0)*z*(z-1.d0)/6.d0
   tw3tw4       = 23.d0/24.d0
   svsx         = 7.d0/6.d0
   theght       = 3.d0/8.d0
   
   w(n1  )      =  theght     * r(n1        ) * h
   w(n1+1*dn)   =  svsx       * r(n1+1*dn   ) * h
   w(n1+2*dn)   =  tw3tw4     * r(n1+2*dn   ) * h
   
   w(n2-5*dn)   =  (tw3tw4*a1 +         a2 +        a3 +         a4)* r(n2-5*dn   ) * h
   w(n2-4*dn)   =  (  svsx*a1 +  tw3tw4*a2 +        a3 +         a4)* r(n2-4*dn   ) * h
   w(n2-3*dn)   =  (theght*a1 +    svsx*a2 + tw3tw4*a3 +         a4)* r(n2-3*dn   ) * h
   w(n2-2*dn)   =  (             theght*a2 +   svsx*a3 +  tw3tw4*a4)* r(n2-2*dn   ) * h
   w(n2-1*dn)   =  (                         theght*a3 +    svsx*a4)* r(n2-1*dn   ) * h
   w(n2  )      =                                         theght*a4 * r(n2        ) * h
   if (n > 6*dn) then
      do i = n1+3*dn,n2-6*dn,dn
         w(i) = r(i) * h
      end do
   end if
99 continue
   end subroutine set_weight_exp3
! === For nrc decomposion. by takto 2012/12/05 =================================
   subroutine set_weight_exp3_3D(ier,n1,n2,dn,r,z,w,ista_nrc,iend_nrc,ist,ien)
!=====================================================================
!
!  M. Okamoto & AAS
!
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: r(n2)
   real(8),intent(in)  :: z
   integer,intent(out) :: ier
   real(8),intent(out) :: w(ista_nrc:iend_nrc)
   integer, intent(in) :: ista_nrc, iend_nrc, ist, ien
   real(8) :: h, a1,a2,a3,a4,tw3tw4,svsx,theght
   integer :: n, i
   ier = 0
   n = n2-n1+1
   if (n < 9*dn) then
      ier = 1 ; go to 99
   end if
   h = log(r(n2)/r(n1)) / dble(n2-n1) * dble(dn)
   a1           = -z*(z-1.d0)*(z-2.d0)/6.d0
   a2           = (z+1.d0)*(z-1.d0)*(z-2.d0)/2.d0
   a3           = -(z+1.d0)*z*(z-2.d0)/2.d0 
   a4           = (z+1.d0)*z*(z-1.d0)/6.d0
   tw3tw4       = 23.d0/24.d0
   svsx         = 7.d0/6.d0
   theght       = 3.d0/8.d0
   
   if(n1 >= ista_nrc .and. n1 <= iend_nrc) then
      w(n1  )      =  theght     * r(n1        ) * h
   end if
   if(n1+1*dn >= ista_nrc .and. n1+1*dn <= iend_nrc) then
      w(n1+1*dn)   =  svsx       * r(n1+1*dn   ) * h
   end if
   if(n1+2*dn >= ista_nrc .and. n1+2*dn <= iend_nrc) then
      w(n1+2*dn)   =  tw3tw4     * r(n1+2*dn   ) * h
   end if
   
   if(n2-5*dn >= ista_nrc .and. n2-5*dn <= iend_nrc) then
      w(n2-5*dn)   =  (tw3tw4*a1 +         a2 +        a3 +         a4)* r(n2-5*dn   ) * h
   end if
   if(n2-4*dn >= ista_nrc .and. n2-4*dn <= iend_nrc) then
      w(n2-4*dn)   =  (  svsx*a1 +  tw3tw4*a2 +        a3 +         a4)* r(n2-4*dn   ) * h
   end if
   if(n2-3*dn >= ista_nrc .and. n2-3*dn <= iend_nrc) then
      w(n2-3*dn)   =  (theght*a1 +    svsx*a2 + tw3tw4*a3 +         a4)* r(n2-3*dn   ) * h
   end if
   if(n2-2*dn >= ista_nrc .and. n2-2*dn <= iend_nrc) then
      w(n2-2*dn)   =  (             theght*a2 +   svsx*a3 +  tw3tw4*a4)* r(n2-2*dn   ) * h
   end if
   if(n2-1*dn >= ista_nrc .and. n2-1*dn <= iend_nrc) then
      w(n2-1*dn)   =  (                         theght*a3 +    svsx*a4)* r(n2-1*dn   ) * h
   end if
   if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
      w(n2  )      =                                         theght*a4 * r(n2        ) * h
   end if
   if (n > 6*dn) then
      do i = max(n1+3*dn,ist),min(n2-6*dn,ien),dn
         w(i) = r(i) * h
      end do
   end if
99 continue
   end subroutine set_weight_exp3_3D
! ==============================================================================

!=====================================================================
   subroutine set_weight_exp4(ier,n1,n2,dn,r,z,w)
!=====================================================================
!
!  M. Okamoto & AAS
!
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: r(*)
   real(8),intent(in)  :: z
   integer,intent(out) :: ier
   real(8),intent(out) :: w(*)
   real(8) :: h, tw3tw4,svsx,theght,omz
   integer :: n, i
   ier = 0
   n = n2-n1+1
   if (n < 7*dn) then
      ier = 1 ; go to 99
   end if
   h = log(r(n2)/r(n1)) / dble(n2-n1) * dble(dn)
   tw3tw4       = 23.d0/24.d0
   svsx         = 7.d0/6.d0
   theght       = 3.d0/8.d0
   omz          = 1.d0-z
   
   w(n1  )      =  theght     * r(n1        ) * h
   w(n1+1*dn)   =  svsx       * r(n1+1*dn   ) * h
   w(n1+2*dn)   =  tw3tw4     * r(n1+2*dn   ) * h
   
   w(n2-3*dn)   =  (tw3tw4*omz +        z)* r(n2-3*dn   ) * h
   w(n2-2*dn)   =  (  svsx*omz + tw3tw4*z)* r(n2-2*dn   ) * h
   w(n2-1*dn)   =  (theght*omz +   svsx*z)* r(n2-1*dn   ) * h
   w(n2  )      =                theght*z * r(n2        ) * h
   if (n > 6*dn) then
      do i = n1+3*dn,n2-4*dn,dn
         w(i) = r(i) * h
      end do
   end if
99 continue
   end subroutine set_weight_exp4

!=====================================================================
   subroutine diff_exp(ier,isdiff,n1,n2,rn,fn,r,f,df,ddf)
!---------------------------------------------------------------------
!
!   Program written by Masakuni Okamoto
!
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: isdiff, n1, n2
   real(8),intent(in)  :: rn(*), fn(*), r
   integer,intent(out) :: ier
   real(8),intent(out) :: f, df, ddf
   real(8) :: denom, a, b, c, x, x1, x2
   integer :: n, i, j
   real(8),allocatable :: p(:,:), dp(:,:), ddp(:,:)
   ier=0
   if ((isdiff < 0).or.(isdiff > 2)) then
      ier = 1 ; go to 99
   end if
   n = n2-n1+1
   if (n <= isdiff) then
      ier = 2 ; go to 99
   end if
  !++++++++++++++++++++++++++++++++++
   allocate(p(n,n),dp(n,n),ddp(n,n))
  !++++++++++++++++++++++++++++++++++
   p(1:n,1) = fn(n1:n2) ; dp(:,:) = 0.d0 ; ddp(:,:) = 0.d0
   a = rn(n1) ; b = rn(n2) ; c = log(b/a) ; x = log(r/a)/c
   if ((isdiff >= 0).and.(n >= 2)) then
      do j = 2,n
         denom = -dble(j-1) / dble(n-1)
         do i = 1,n+1-j
            x1 = log(rn(i+n1+j-2)/a) / c
            x2 = log(rn(i+n1-1)  /a) / c
            p(i,j) = ((x-x1)*p(i,j-1) + (x2-x)*p(i+1,j-1)) / denom
         end do
      end do
   end if
   if (isdiff >= 1) then
      do j = 2,n
         denom = -dble(j-1)/dble(n-1)
         do i = 1,n+1-j
            x1 = log(rn(i+n1+j-2)/a) / c
            x2 = log(rn(i+n1-1)  /a) / c
            dp(i,j) = (x-x1)*dp(i,j-1) + (x2-x)*dp(i+1,j-1) &
                    + (p(i,j-1)-p(i+1,j-1))
            dp(i,j) = dp(i,j) / denom
         end do
      end do
   end if
   if (isdiff >= 2) then
      do j = 2,n
         denom = -dble(j-1)/dble(n-1)
         do i = 1,n+1-j
            x1 = log(rn(i+n1+j-2)/a)/c
            x2 = log(rn(i+n1-1)  /a)/c
            ddp(i,j) = (x-x1)*ddp(i,j-1) + (x2-x)*ddp(i+1,j-1) &
                     + 2.d0*(dp(i,j-1)-dp(i+1,j-1))
            ddp(i,j) = ddp(i,j) / denom
         end do
      end do
   end if
   f   = p  (1,n)
   df  = dp (1,n) * (1.d0/r/c)
   ddf = ddp(1,n) * (1.d0/r/c)**2 - dp(1,n) * (1.d0/r/r/c)
99 continue
  !+++++++++++++++++++++++
   if (allocated(p)) then
   deallocate(p,dp,ddp)
   end if
  !+++++++++++++++++++++++
   end subroutine diff_exp

!=====================================================================
   subroutine diff3_exp(ier,isdiff,n1,n2,rn,fn,r,f,df,ddf,dddf)
!---------------------------------------------------------------------
!
!   Program written by Masakuni Okamoto
!
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: isdiff, n1, n2
   real(8),intent(in)  :: rn(*), fn(*), r
   integer,intent(out) :: ier
   real(8),intent(out) :: f, df, ddf, dddf
   real(8) :: denom, a, b, c, d, x, x1, x2
   integer :: n, i, j
   real(8),allocatable :: &
      p(:,:), dp(:,:), ddp(:,:), dddp(:,:)
   ier=0
   if ((isdiff < 0).or.(isdiff > 3)) then
      ier = 1 ; go to 99
   end if
   n = n2-n1+1
   if (n <= isdiff) then
      ier = 2 ; go to 99
   end if
  !++++++++++++++++++++++++++++++++++++++++++++
   allocate(p(n,n),dp(n,n),ddp(n,n),dddp(n,n))
  !++++++++++++++++++++++++++++++++++++++++++++
   p(1:n,1) = fn(n1:n2) ; dp(:,:) = 0.d0 ; ddp(:,:) = 0.d0
   dddp(:,:) = 0.d0
   a = rn(n1) ; b = rn(n2) ; c = log(b/a) ; x = log(r/a)/c
   d = 1.d0/c
   if ((isdiff >= 0).and.(n >= 2)) then
      do j = 2,n
         denom = -dble(j-1) / dble(n-1)
         do i = 1,n+1-j
            x1 = log(rn(i+n1+j-2)/a) * d
            x2 = log(rn(i+n1-1)  /a) * d
            p(i,j) = ((x-x1)*p(i,j-1) + (x2-x)*p(i+1,j-1)) / denom
         end do
      end do
   end if
   if (isdiff >= 1) then
      do j = 2,n
         denom = -dble(j-1)/dble(n-1)
         do i = 1,n+1-j
            x1 = log(rn(i+n1+j-2)/a) * d
            x2 = log(rn(i+n1-1)  /a) * d
            dp(i,j) = (x-x1)*dp(i,j-1) + (x2-x)*dp(i+1,j-1) &
                    + (p(i,j-1)-p(i+1,j-1))
            dp(i,j) = dp(i,j) / denom
         end do
      end do
   end if
   if (isdiff >= 2) then
      do j = 2,n
         denom = -dble(j-1)/dble(n-1)
         do i = 1,n+1-j
            x1 = log(rn(i+n1+j-2)/a) * d
            x2 = log(rn(i+n1-1)  /a) * d
            ddp(i,j) = (x-x1)*ddp(i,j-1) + (x2-x)*ddp(i+1,j-1) &
                     + 2.d0*(dp(i,j-1)-dp(i+1,j-1))
            ddp(i,j) = ddp(i,j) / denom
         end do
      end do
   end if
   if (isdiff >= 3) then
      do j = 2,n
         denom = -dble(j-1)/dble(n-1)
         do i = 1,n+1-j
            x1 = log(rn(i+n1+j-2)/a) * d
            x2 = log(rn(i+n1-1)  /a) * d
            dddp(i,j) = (x-x1)*dddp(i,j-1) + (x2-x)*dddp(i+1,j-1) &
                     + 3.d0*(ddp(i,j-1)-ddp(i+1,j-1))
            dddp(i,j) = dddp(i,j) / denom
         end do
      end do
   end if
   f     = p(1,n)
   df    = dp(1,n) * d
   ddf   = ddp(1,n) * d*d - dp(1,n) * d
   dddf  = dddp(1,n) * d*d*d - ddp(1,n) * 3.d0*d*d + dp(1,n) * 2.d0*d
   df    = df / r
   ddf   = ddf / (r*r)
   dddf  = dddf / (r*r*r)
99 continue
  !++++++++++++++++++++++++++
   if (allocated(p)) then
   deallocate(p,dp,ddp,dddp)
   end if
  !++++++++++++++++++++++++++
   end subroutine diff3_exp

!=====================================================================
   subroutine diff4_exp(ier,isdiff,n1,n2,rn,fn,r,f,df,ddf,dddf,ddddf)
!---------------------------------------------------------------------
!
!   Program written by Masakuni Okamoto
!
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: isdiff, n1, n2
   real(8),intent(in)  :: rn(*), fn(*), r
   integer,intent(out) :: ier
   real(8),intent(out) :: f, df, ddf, dddf, ddddf
   real(8) :: denom, a, b, c, d, x, x1, x2
   integer :: n, i, j
   real(8),allocatable :: &
      p(:,:), dp(:,:), ddp(:,:), dddp(:,:), ddddp(:,:)
   ier=0
   if ((isdiff < 0).or.(isdiff > 4)) then
      ier = 1 ; go to 99
   end if
   n = n2-n1+1
   if (n <= isdiff) then
      ier = 2 ; go to 99
   end if
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
   allocate(p(n,n),dp(n,n),ddp(n,n),dddp(n,n),ddddp(n,n))
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
   p(1:n,1) = fn(n1:n2) ; dp(:,:) = 0.d0 ; ddp(:,:) = 0.d0
   dddp(:,:) = 0.d0 ; ddddp(:,:) = 0.d0
   a = rn(n1) ; b = rn(n2) ; c = log(b/a) ; x = log(r/a)/c
   d = 1.d0/c
   if ((isdiff >= 0).and.(n >= 2)) then
      do j = 2,n
         denom = -dble(j-1) / dble(n-1)
         do i = 1,n+1-j
            x1 = log(rn(i+n1+j-2)/a) * d
            x2 = log(rn(i+n1-1)  /a) * d
            p(i,j) = ((x-x1)*p(i,j-1) + (x2-x)*p(i+1,j-1)) / denom
         end do
      end do
   end if
   if (isdiff >= 1) then
      do j = 2,n
         denom = -dble(j-1)/dble(n-1)
         do i = 1,n+1-j
            x1 = log(rn(i+n1+j-2)/a) * d
            x2 = log(rn(i+n1-1)  /a) * d
            dp(i,j) = (x-x1)*dp(i,j-1) + (x2-x)*dp(i+1,j-1) &
                    + (p(i,j-1)-p(i+1,j-1))
            dp(i,j) = dp(i,j) / denom
         end do
      end do
   end if
   if (isdiff >= 2) then
      do j = 2,n
         denom = -dble(j-1)/dble(n-1)
         do i = 1,n+1-j
            x1 = log(rn(i+n1+j-2)/a) * d
            x2 = log(rn(i+n1-1)  /a) * d
            ddp(i,j) = (x-x1)*ddp(i,j-1) + (x2-x)*ddp(i+1,j-1) &
                     + 2.d0*(dp(i,j-1)-dp(i+1,j-1))
            ddp(i,j) = ddp(i,j) / denom
         end do
      end do
   end if
   if (isdiff >= 3) then
      do j = 2,n
         denom = -dble(j-1)/dble(n-1)
         do i = 1,n+1-j
            x1 = log(rn(i+n1+j-2)/a) * d
            x2 = log(rn(i+n1-1)  /a) * d
            dddp(i,j) = (x-x1)*dddp(i,j-1) + (x2-x)*dddp(i+1,j-1) &
                     + 3.d0*(ddp(i,j-1)-ddp(i+1,j-1))
            dddp(i,j) = dddp(i,j) / denom
         end do
      end do
   end if
   if (isdiff >= 4) then
      do j = 2,n
         denom = -dble(j-1)/dble(n-1)
         do i = 1,n+1-j
            x1 = log(rn(i+n1+j-2)/a) * d
            x2 = log(rn(i+n1-1)  /a) * d
            ddddp(i,j) = (x-x1)*ddddp(i,j-1) + (x2-x)*ddddp(i+1,j-1) &
                     + 4.d0*(dddp(i,j-1)-dddp(i+1,j-1))
            ddddp(i,j) = ddddp(i,j) / denom
         end do
      end do
   end if
   f     = p(1,n)
   df    = dp(1,n) * d
   ddf   = ddp(1,n) * d*d - dp(1,n) * d
   dddf  = dddp(1,n) * d*d*d - ddp(1,n) * 3.d0*d*d + dp(1,n) * 2.d0*d
   ddddf = ddddp(1,n) * d*d*d*d - dddp(1,n) * 6.d0*d*d*d &
         + ddp(1,n) * 11.d0*d*d - dp(1,n) * 6.d0*d
   df    = df / r
   ddf   = ddf / (r*r)
   dddf  = dddf / (r*r*r)
   ddddf = ddddf / (r*r*r*r)
99 continue
  !++++++++++++++++++++++++++++++++
   if (allocated(p)) then
   deallocate(p,dp,ddp,dddp,ddddp)
   end if
  !++++++++++++++++++++++++++++++++
   end subroutine diff4_exp

!=====================================================================
   subroutine calc_diff_exp(ier,iord,n,rn,fn,dfn)
!---------------------------------------------------------------------
!
!  M. Okamoto
!
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: iord, n
!   real(8),intent(in)  :: rn(*)
   real(8),intent(in)  :: rn(*), fn(*)
   integer,intent(out) :: ier
!   real(8),intent(out) :: fn(*), dfn(*)
   real(8),intent(out) :: dfn(*)
   real(8) :: r, f, df, ddf
   integer :: isdiff, n1, n2, i
   ier = 0
   isdiff = 1
   do i = 1,n
      r = rn(i)
      if (i < iord+1) then
         n1 = 1 ; n2 = 1 + 2*iord
      else if (i > n-iord) then
         n1 = n - 2*iord ; n2 = n
      else
         n1 = i - iord ; n2 = i + iord
      end if
      call diff_exp(ier,isdiff,n1,n2,rn,fn,r,f,df,ddf)
      dfn(i) = df
   end do
99 continue
   end subroutine calc_diff_exp

!=====================================================================
   subroutine calc_ddiff_exp(ier,iord,n,rn,fn,dfn,ddfn)
!---------------------------------------------------------------------
!
!  M. Okamoto
!
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: iord, n
   real(8),intent(in)  :: rn(*), fn(*)
   integer,intent(out) :: ier
!   real(8),intent(out) :: fn(*), dfn(*), ddfn(*)                   !AAS 2009
   real(8),intent(out) :: dfn(*), ddfn(*)
   real(8) :: r, f, df, ddf
   integer :: isdiff, n1, n2, i
   ier = 0
   isdiff = 2
   do i = 1,n
      r = rn(i)
      if (i < iord+1) then
         n1 = 1 ; n2 = 1 + 2*iord
      else if (i > n-iord) then
         n1 = n - 2*iord ; n2 = n
      else
         n1 = i - iord ; n2 = i + iord
      end if
      call diff_exp(ier,isdiff,n1,n2,rn,fn,r,f,df,ddf)
      dfn(i) = df ; ddfn(i) = ddf
   end do
99 continue
   end subroutine calc_ddiff_exp
   
!=====================================================================
   subroutine calc_ddiff_exp2(ier,iord,n,dn,rn,fn,dfn,ddfn)
!---------------------------------------------------------------------
!
!  M. Okamoto
!
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: iord, n, dn
   real(8),intent(in)  :: rn(*), fn(*)
   integer,intent(out) :: ier
!   real(8),intent(out) :: fn(*), dfn(*), ddfn(*)                   !AAS 2009
   real(8),intent(out) :: dfn(*), ddfn(*)
   real(8) :: r, f, df, ddf
   real(8),pointer,dimension(:):: rnn,fnn
   integer :: isdiff, n1, n2, i, j, k
   ier = 0
   isdiff = 2
   allocate(rnn(1+2*iord))
   allocate(fnn(1+2*iord))
   do i = 1,n,dn
      r = rn(i)
      if (i < 1+iord*dn) then
         n1 = 1 ; n2 = 1 + 2*iord*dn
      else if (i > n-iord*dn) then
         n1 = n - 2*iord*dn ; n2 = n
      else
         n1 = i - iord*dn ; n2 = i + iord*dn
      end if
      k=0
      do j=n1,n2,dn
        k=k+1
        rnn(k)=rn(j)
        fnn(k)=fn(j)
      end do
      call diff_exp(ier,isdiff,1,k,rnn,fnn,r,f,df,ddf)
      dfn(i) = df ; ddfn(i) = ddf
   end do
99 continue
   deallocate(rnn,fnn)
   end subroutine calc_ddiff_exp2
   
!=====================================================================
   subroutine calc_diff_exp2(ier,iord,n,dn,rn,fn,dfn)
!---------------------------------------------------------------------
!
!  M. Okamoto
!
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: iord, n, dn
   real(8),intent(in)  :: rn(*), fn(*)
   integer,intent(out) :: ier
!   real(8),intent(out) :: fn(*), dfn(*), ddfn(*)                   !AAS 2009
   real(8),intent(out) :: dfn(*)
   real(8) :: r, f, df, ddf
   real(8),pointer,dimension(:):: rnn,fnn
   integer :: isdiff, n1, n2, i, j, k
   ier = 0
   isdiff = 1
   allocate(rnn(1+2*iord))
   allocate(fnn(1+2*iord))
   do i = 1,n,dn
      r = rn(i)
      if (i < 1+iord*dn) then
         n1 = 1 ; n2 = 1 + 2*iord*dn
      else if (i > n-iord*dn) then
         n1 = n - 2*iord*dn ; n2 = n
      else
         n1 = i - iord*dn ; n2 = i + iord*dn
      end if
      k=0
      do j=n1,n2,dn
        k=k+1
        rnn(k)=rn(j)
        fnn(k)=fn(j)
      end do
      call diff_exp(ier,isdiff,1,k,rnn,fnn,r,f,df,ddf)
      dfn(i) = df
   end do
99 continue
   deallocate(rnn,fnn)
   end subroutine calc_diff_exp2
   
  subroutine calc_diff_exp2_3D(isdiff,iord,n,dnr,xh,rn,fn,dfn,ista_nrc,iend_nrc,ddfn)
    implicit none
    integer,intent(in)  :: isdiff,iord, n, dnr
    real(8),intent(in)  :: xh,rn(n), fn(ista_nrc-6*dnr:iend_nrc+6*dnr)
    real(8),intent(out) :: dfn(ista_nrc:iend_nrc)
    integer, intent(in) :: ista_nrc, iend_nrc
    real(8),optional,intent(inout):: ddfn(ista_nrc:iend_nrc)
    real(8) :: r, f, df,ddf,h
    real(8),allocatable :: p(:,:), dp(:,:), ddp(:,:)
    integer :: n1, n2, i, pdim, indx
    integer       :: id_sname = -1

    pdim = 1+2*iord
    h = 1.d0/xh
    allocate(p(pdim,pdim),dp(pdim,pdim))
    if(isdiff >=2)  allocate(ddp(pdim,pdim))

!!$    do i = 1,n
    do i = max(1,ista_nrc),min(iend_nrc,n)
       r = rn(i)
       if (i < iord+1) then
          n1 = 1 ; n2 = 1 + 2*iord
       else if (i > n-iord) then
          n1 = n - 2*iord ; n2 = n
       else
          n1 = i - iord ; n2 = i + iord
       end if
       indx = i
      ! indx-n1 = i-(i-iord) = iord
!!$       call diff_exp2(isdiff,n1,n2,rn,fn,r,df,ddf)
       call diff_exp2(isdiff,n1,n2,fn,r,df,ddf)  ! isdiff,n1,n2,fn,r -> df,ddf
!!$      call diff_exp(ier,isdiff,n1,n2,rn,fn,r,f,df,ddf)
       dfn(i) = df
       if(isdiff==2) ddfn(i) = ddf
    end do
    if(isdiff >=2) deallocate(ddp)
    deallocate(p,dp)
  contains
    subroutine diff_exp2(isdiff,n1,n2,fn,r,df,ddf) ! isdiff,n1,n2,fn,r
!---------------------------------------------------------------------
!
!   Program written by Masakuni Okamoto
!  Revised for logarithmic grid by Takahiro Yamasaki, 2010/06/10 
!
!---------------------------------------------------------------------
      integer,intent(in)  :: isdiff,n1,n2
      real(8),intent(in)  :: fn(ista_nrc-6:iend_nrc+6), r
      real(8),intent(out) :: df,ddf
      real(8) :: denom_inv, a, b, c, x, x1, x2, dn2minusn1
      integer :: n, i, j, n1_delta
      n = n2-n1+1
      p(1:n,1) = fn(n1:n2) ;  dp(:,:) = 0.d0 
      if(isdiff>=2) ddp(:,:) = 0.d0
      c = h*(n2-n1)
      n1_delta =indx-n1
      dn2minusn1 = 1.d0/dble(n2-n1)

      if ((isdiff >= 0).and.(n >= 2)) then
         do j = 2,n
            denom_inv = -dble(n-1)/dble(j-1)
            do i = 1,n+1-j
               p(i,j) = ((n1_delta-i-j+2)*p(i,j-1) + (i-1-n1_delta)*p(i+1,j-1))*dn2minusn1 * denom_inv
            end do
         end do
      end if
      if (isdiff >= 1) then
         do j = 2,n
            denom_inv = -dble(n-1)/dble(j-1)
            do i = 1,n+1-j
               dp(i,j) = ((n1_delta-i-j+2)*dp(i,j-1) + (i-1-n1_delta)*dp(i+1,j-1))*dn2minusn1 &
                    &     + (p(i,j-1)-p(i+1,j-1))
               dp(i,j) = dp(i,j) * denom_inv
            end do
         end do
      end if
      if (isdiff >= 2) then
         do j = 2,n
            denom_inv = -dble(n-1)/dble(j-1)
            do i = 1,n+1-j
               ddp(i,j) = ((n1_delta-i-j+2)*ddp(i,j-1) + (i-1-n1_delta)*ddp(i+1,j-1))*dn2minusn1 &
                    &      + 2.d0*(dp(i,j-1)-dp(i+1,j-1))
               ddp(i,j) = ddp(i,j) * denom_inv
            end do
         end do
      end if
      df  = dp (1,n) * (1.d0/r/c)
      if(isdiff >=2) ddf = ddp(1,n) * (1.d0/r/c)**2 - dp(1,n) * (1.d0/r/r/c)
    end subroutine diff_exp2

  end subroutine calc_diff_exp2_3D

!=====================================================================
   subroutine calc_diff_exp3(ier,iord,n,dn,rn,fn,dfn)
!---------------------------------------------------------------------
    implicit none
    integer,intent(in)  :: iord, n, dn
    real(8),intent(in)  :: rn(*), fn(*)
    integer,intent(out) :: ier
    real(8),intent(out) :: dfn(*)
    
    integer:: ir
    real(8):: idh
    
    if(dn.eq.1) then
        select case (iord)
        case (1)
            call diff_3points2(ier,1,n,dn,fn,dfn) 
        case (2)
            call diff_5points2(ier,1,n,dn,fn,dfn) 
        case (3)
            call diff_7points2(ier,1,n,dn,fn,dfn) 
        end select
    else
        select case (iord)
        case (1)
            call diff_3points3(ier,1,n,dn,fn,dfn) 
        case (2)
            call diff_5points3(ier,1,n,dn,fn,dfn) 
        case (3)
            call diff_7points3(ier,1,n,dn,fn,dfn) 
        end select
    end if
        
    idh=dble(n-1)/log(rn(n)/rn(1))/dble(dn)
    
    do ir=1,n,dn
        dfn(ir)=dfn(ir)/rn(ir)*idh
    end do
    return
    end subroutine calc_diff_exp3
! === For nrc decomposion. by takto 2012/12/05 =================================
    subroutine calc_diff_exp3_3D(ier,iord,n,dn,rn,fn,dfn,ista_nrc,iend_nrc,ist,ien)
!---------------------------------------------------------------------
    implicit none
    integer,intent(in)  :: iord, n, dn
!   real(8),intent(in)  :: rn(n), fn(ista_nrc-3*dn:iend_nrc+3*dn)
    real(8),intent(in)  :: rn(n), fn(ista_nrc-6*dn:iend_nrc+6*dn)
    integer,intent(out) :: ier
    real(8),intent(out) :: dfn(ista_nrc:iend_nrc)
    integer,intent(in)  :: ista_nrc, iend_nrc, ist, ien

    integer:: ir
    real(8):: idh

    if(dn.eq.1) then
        select case (iord)
        case (1)
            call diff_3points2_3D(ier,1,n,dn,fn,dfn,ista_nrc,iend_nrc,ist,ien)
        case (2)
            call diff_5points2_3D(ier,1,n,dn,fn,dfn,ista_nrc,iend_nrc,ist,ien)
        case (3)
            call diff_7points2_3D(ier,1,n,dn,fn,dfn,ista_nrc,iend_nrc,ist,ien)
        end select
    else
        select case (iord)
        case (1)
            call diff_3points3_3D(ier,1,n,dn,fn,dfn,ista_nrc,iend_nrc,ist,ien)
        case (2)
            call diff_5points3_3D(ier,1,n,dn,fn,dfn,ista_nrc,iend_nrc,ist,ien)
        case (3)
            call diff_7points3_3D(ier,1,n,dn,fn,dfn,ista_nrc,iend_nrc,ist,ien)
        end select
    end if

    idh=dble(n-1)/log(rn(n)/rn(1))/dble(dn)

    do ir=ist,ien,dn
        dfn(ir)=dfn(ir)/rn(ir)*idh
    end do
    return
    end subroutine calc_diff_exp3_3D
! ==============================================================================
    
!=====================================================================
   subroutine calc_diff_exp4(ier,iord,n1,n2,dn,rn,fn,dfn)
!---------------------------------------------------------------------
    implicit none
    integer,intent(in)  :: iord, n1, n2, dn
    real(8),intent(in)  :: rn(*), fn(*)
    integer,intent(out) :: ier
    real(8),intent(out) :: dfn(*)
    
    integer:: ir
    real(8):: idh
    
    if(dn.eq.1) then
        select case (iord)
        case (1)
            call diff_3points2(ier,n1,n2,dn,fn,dfn) 
        case (2)
            call diff_5points2(ier,n1,n2,dn,fn,dfn) 
        case (3)
            call diff_7points2(ier,n1,n2,dn,fn,dfn) 
        end select
    else
        select case (iord)
        case (1)
            call diff_3points3(ier,n1,n2,dn,fn,dfn) 
        case (2)
            call diff_5points3(ier,n1,n2,dn,fn,dfn) 
        case (3)
            call diff_7points3(ier,n1,n2,dn,fn,dfn) 
        end select
    end if
        
    idh=dble(n2-n1)/log(rn(n2)/rn(n1))/dble(dn)
    
    do ir=n1,n2,dn
        dfn(ir)=dfn(ir)/rn(ir)*idh
    end do
    return
    end subroutine calc_diff_exp4
    
!=====================================================================
   subroutine calc_ddiff_exp3(ier,iord,n,dn,rn,fn,dfn,ddfn)
!---------------------------------------------------------------------
    implicit none
    integer,intent(in)  :: iord, n, dn
    real(8),intent(in)  :: rn(*), fn(*)
    integer,intent(out) :: ier
    real(8),intent(out) :: dfn(*)
    real(8),intent(out) :: ddfn(*)
    
    integer:: ir
    real(8):: idh,irn
   
   if(dn.eq.1) then
        select case (iord)
        case (1)
            call ddiff_3points2(ier,1,n,dn,fn,dfn,ddfn) 
        case (2)
            call ddiff_5points2(ier,1,n,dn,fn,dfn,ddfn) 
        case (3)
            call ddiff_7points2(ier,1,n,dn,fn,dfn,ddfn) 
        end select
    else
        select case (iord)
        case (1)
            call ddiff_3points3(ier,1,n,dn,fn,dfn,ddfn) 
        case (2)
            call ddiff_5points3(ier,1,n,dn,fn,dfn,ddfn) 
        case (3)
            call ddiff_7points3(ier,1,n,dn,fn,dfn,ddfn) 
        end select
    end if
    
!    idh=1.d0/dh
    idh=dble(n-1)/log(rn(n)/rn(1))/dble(dn)
    
    do ir=1,n,dn
        irn=1.d0/rn(ir)
        dfn(ir)=dfn(ir)*irn*idh
        ddfn(ir)=(ddfn(ir)*idh*idh*irn-dfn(ir))*irn
    end do
    return
    end subroutine calc_ddiff_exp3
! === For nrc decomposion. by takto 2012/12/05 =================================
    subroutine calc_ddiff_exp3_3D(ier,iord,n,dn,rn,fn,dfn,ddfn,ista_nrc,iend_nrc,ist,ien)
!---------------------------------------------------------------------
    implicit none
    integer,intent(in)  :: iord, n, dn
    real(8),intent(in)  :: rn(n), fn(ista_nrc-6*dn:iend_nrc+6*dn)
    integer,intent(out) :: ier
    real(8),intent(out) :: dfn(ista_nrc:iend_nrc)
    real(8),intent(out) :: ddfn(ista_nrc:iend_nrc)
    integer,intent(in)  :: ista_nrc, iend_nrc, ist, ien

    integer:: ir
    real(8):: idh,irn

    if(dn.eq.1) then
        select case (iord)
        case (1)
            call ddiff_3points2_3D(ier,1,n,dn,fn,dfn,ddfn,ista_nrc,iend_nrc,ist,ien)
        case (2)
            call ddiff_5points2_3D(ier,1,n,dn,fn,dfn,ddfn,ista_nrc,iend_nrc,ist,ien)
        case (3)
            call ddiff_7points2_3D(ier,1,n,dn,fn,dfn,ddfn,ista_nrc,iend_nrc,ist,ien)
        end select
    else
        select case (iord)
        case (1)
            call ddiff_3points3_3D(ier,1,n,dn,fn,dfn,ddfn,ista_nrc,iend_nrc,ist,ien)
        case (2)
            call ddiff_5points3_3D(ier,1,n,dn,fn,dfn,ddfn,ista_nrc,iend_nrc,ist,ien)
        case (3)
            call ddiff_7points3_3D(ier,1,n,dn,fn,dfn,ddfn,ista_nrc,iend_nrc,ist,ien)
        end select
    end if

    idh=dble(n-1)/log(rn(n)/rn(1))/dble(dn)

    do ir=ist,ien,dn
        irn=1.d0/rn(ir)
        dfn(ir)=dfn(ir)*irn*idh
        ddfn(ir)=(ddfn(ir)*idh*idh*irn-dfn(ir))*irn
    end do
    return
    end subroutine calc_ddiff_exp3_3D
! ==============================================================================
    
!=====================================================================
   subroutine calc_ddiff_exp4(ier,iord,n1,n2,dn,rn,fn,dfn,ddfn)
!---------------------------------------------------------------------
    implicit none
    integer,intent(in)  :: iord, n1, n2, dn
    real(8),intent(in)  :: rn(*), fn(*)
    integer,intent(out) :: ier
    real(8),intent(out) :: dfn(*)
    real(8),intent(out) :: ddfn(*)
    
    integer:: ir
    real(8):: idh,irn
   
   if(dn.eq.1) then
        select case (iord)
        case (1)
            call ddiff_3points2(ier,n1,n2,dn,fn,dfn,ddfn) 
        case (2)
            call ddiff_5points2(ier,n1,n2,dn,fn,dfn,ddfn) 
        case (3)
            call ddiff_7points2(ier,n1,n2,dn,fn,dfn,ddfn) 
        end select
    else
        select case (iord)
        case (1)
            call ddiff_3points3(ier,n1,n2,dn,fn,dfn,ddfn) 
        case (2)
            call ddiff_5points3(ier,n1,n2,dn,fn,dfn,ddfn) 
        case (3)
            call ddiff_7points3(ier,n1,n2,dn,fn,dfn,ddfn) 
        end select
    end if
    
!    idh=1.d0/dh
    idh=dble(n2-n1)/log(rn(n2)/rn(n1))/dble(dn)
    
    do ir=n1,n2,dn
        irn=1.d0/rn(ir)
        dfn(ir)=dfn(ir)*irn*idh
        ddfn(ir)=(ddfn(ir)*idh*idh*irn-dfn(ir))*irn
    end do
    return
    end subroutine calc_ddiff_exp4
    
!=====================================================================
   subroutine diff_3points(ier,n1,n2,fn,df)                       
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2
   real(8),intent(in)  :: fn(n1:n2)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(n1:n2)
   
   integer:: i
   real(8):: f0,f1,f2
   
   ier=0
   f0=fn(n1)
   f1=fn(n1+1)
   f2=fn(n1+2)
   df(n1)=-0.5d0*f2+2.d0*f1-1.5d0*f0
   do i=n1+1,n2-2
        df(i)=0.5d0*(f2-f0)
        f0=f1
        f1=f2
        f2=fn(i+2)
   end do
   df(n2-1)=0.5d0*(f2-f0)
   df(n2)=1.5d0*f2-2.d0*f1+0.5d0*f0
   return
   end subroutine diff_3points
   
!=====================================================================
   subroutine diff_3points2(ier,n1,n2,dn,fn,df)                       
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: fn(*)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(*)
   
   integer:: i,j
   real(8):: f0,f1,f2
   
    ier=0
    f0=fn(n1)
    f1=fn(n1+1)
    f2=fn(n1+2)
    df(n1)=-0.5d0*f2+2.d0*f1-1.5d0*f0
   
    select case (dn)
    case (1)
        do i=n1+1,n2-2
            df(i)=0.5d0*(f2-f0)
            f0=f1
            f1=f2
            f2=fn(i+2)
        end do
        df(n2-1)=0.5d0*(f2-f0)
        df(n2)=1.5d0*f2-2.d0*f1+0.5d0*f0
    case (2)
        f0=f1
        do i=n1+2,n2-1,2
            f2=fn(i+1)
            df(i)=0.5d0*(f2-f0)
            f0=f2
        end do
        if(i.eq.n2) df(n2)=1.5d0*fn(n2)-2.d0*fn(n2-1)+0.5d0*fn(n2-2)
    case default
        do i=n1+dn,n2-1,dn
            df(i)=0.5d0*(fn(i+1)-fn(i-1))
        end do
        if(i.eq.n2) df(n2)=1.5d0*fn(n2)-2.d0*fn(n2-1)+0.5d0*fn(n2-2)
    end select
    return
   end subroutine diff_3points2
! === For nrc decomposion. by takto 2012/12/05 =================================
   subroutine diff_3points2_3D(ier,n1,n2,dn,fn,df,ista_nrc,iend_nrc,ist,ien)
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
!  real(8),intent(in)  :: fn(ista_nrc-3*dn:iend_nrc+3*dn)
   real(8),intent(in)  :: fn(ista_nrc-6*dn:iend_nrc+6*dn)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(ista_nrc:iend_nrc)
   integer,intent(in)  :: ista_nrc, iend_nrc, ist, ien

   integer:: i,j

   ier=0

   if(n1 >= ista_nrc .and. n1 <= iend_nrc) then
      df(n1)=-0.5d0*fn(n1+2)+2.d0*fn(n1+1)-1.5d0*fn(n1)
   end if
   
   select case (dn)
      case (1)
         do i=max(n1+1,ist),min(n2-2,ien)
            df(i)=0.5d0*(fn(i+1)-fn(i-1))
         end do
         i = (n2-2) - (n1+1) + 1
         i = (n1+1) + i
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            df(n2-1)=0.5d0*(fn(i+1)-fn(i-1))
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            df(n2)=1.5d0*fn(i+1)-2.d0*fn(i)+0.5d0*fn(i-1)
         end if
      case (2)
         do i=max(n1+2,ist),min(n2-1,ien),2
            df(i)=0.5d0*(fn(i+1)-fn(i-1))
         end do
         i = ((n2-1) - (n1+2) + 2)/2
         i = (n1+2) + 2*i
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.eq.n2) df(n2)=1.5d0*fn(n2)-2.d0*fn(n2-1)+0.5d0*fn(n2-2)
         end if
      case default
         do i=max(n1+dn,ist),min(n2-1,ien),dn
            df(i)=0.5d0*(fn(i+1)-fn(i-1))
         end do
         i = ((n2-1) - (n1+dn) + dn)/dn
         i = (n1+dn) + dn*i
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.eq.n2) df(n2)=1.5d0*fn(n2)-2.d0*fn(n2-1)+0.5d0*fn(n2-2)
         end if
   end select

   return
   end subroutine diff_3points2_3D
! ==============================================================================
   
!=====================================================================
   subroutine diff_3points3(ier,n1,n2,dn,fn,df)                       
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: fn(*)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(*)
   
   integer:: i
   real(8):: f0,f1,f2
   
   ier=0
   f0=fn(n1)
   f1=fn(n1+dn)
   f2=fn(n1+2*dn)
   df(n1)=-0.5d0*f2+2.d0*f1-1.5d0*f0
   do i=n1+dn,n2-2*dn,dn
        df(i)=0.5d0*(f2-f0)
        f0=f1
        f1=f2
        f2=fn(i+2*dn)
   end do
   df(n2-dn)=0.5d0*(f2-f0)
   df(n2)=1.5d0*f2-2.d0*f1+0.5d0*f0
   return
   end subroutine diff_3points3
! === For nrc decomposion. by takto 2012/12/05 =================================
   subroutine diff_3points3_3D(ier,n1,n2,dn,fn,df,ista_nrc,iend_nrc,ist,ien)                       
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
!  real(8),intent(in)  :: fn(ista_nrc-3*dn:iend_nrc+3*dn)
   real(8),intent(in)  :: fn(ista_nrc-6*dn:iend_nrc+6*dn)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(ista_nrc:iend_nrc)
   integer,intent(in)  :: ista_nrc, iend_nrc, ist, ien

   integer:: i

   ier=0

   if(n1 >= ista_nrc .and. n1 <= iend_nrc) then
      df(n1)=-0.5d0*fn(n1+2*dn)+2.d0*fn(n1+dn)-1.5d0*fn(n1)
   end if
   do i=max(n1+dn,ist),min(n2-2*dn,ien),dn
      df(i)=0.5d0*(fn(i+dn)-fn(i-dn))
   end do
   i = ((n2-2*dn) - (n1+dn) + dn)/dn
   i = (n1+dn) + dn*i
   if(n2-dn >= ista_nrc .and. n2-dn <= iend_nrc) then
      df(n2-dn)=0.5d0*(fn(i+dn)-fn(i-dn))
   end if
   if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
      df(n2)=1.5d0*fn(i+dn)-2.d0*fn(i)+0.5d0*fn(i-dn)
   end if
   return
   end subroutine diff_3points3_3D
! ==============================================================================
   
!=====================================================================
   subroutine ddiff_3points(ier,n1,n2,fn,df,ddf)                       
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2
   real(8),intent(in)  :: fn(n1:n2)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(n1:n2)
   real(8),intent(out) :: ddf(n1:n2)
   
   integer:: i
   real(8):: f0,f1,f2
   
   ier=0
   f0=fn(n1)
   f1=fn(n1+1)
   f2=fn(n1+2)
   df(n1)=-0.5d0*f2+2.d0*f1-1.5d0*f0
   ddf(n1)=-fn(n1+3)+4.d0*f2-5.d0*f1+2.d0*f0
   do i=n1+1,n2-2
        df(i)=0.5d0*(f2-f0)
        ddf(i)=f2-2.d0*f1+f0
        f0=f1
        f1=f2
        f2=fn(i+2)
   end do
   df(n2-1)=0.5d0*(f2-f0)
   ddf(n2-1)=f2-2.d0*f1+f0
   df(n2)=1.5d0*f2-2.d0*f1+0.5d0*f0
   ddf(n2)=2.d0*f2-5.d0*f1+4.d0*f0-fn(n2-3)
   return
   end subroutine ddiff_3points
   
!=====================================================================
   subroutine ddiff_3points2(ier,n1,n2,dn,fn,df,ddf)                       
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: fn(*)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(*)
   real(8),intent(out) :: ddf(*)
   
   integer:: i
   real(8):: f0,f1,f2
   
   ier=0
   f0=fn(n1)
   f1=fn(n1+1)
   f2=fn(n1+2)
   df(n1)=-0.5d0*f2+2.d0*f1-1.5d0*f0
   ddf(n1)=-fn(n1+3)+4.d0*f2-5.d0*f1+2.d0*f0
   
    select case(dn)
    case(1)
        do i=n1+1,n2-2
            df(i)=0.5d0*(f2-f0)
            ddf(i)=f2-2.d0*f1+f0
            f0=f1
            f1=f2
            f2=fn(i+2)
        end do
        df(n2-1)=0.5d0*(f2-f0)
        ddf(n2-1)=f2-2.d0*f1+f0
        df(n2)=1.5d0*f2-2.d0*f1+0.5d0*f0
        ddf(n2)=2.d0*f2-5.d0*f1+4.d0*f0-fn(n2-3)
    case(2)
        f0=f1
        do i=n1+2,n2-1,2
            f1=fn(i)
            f2=fn(i+1)
            df(i)=0.5d0*(f2-f0)
            ddf(i)=f2-2.d0*f1+f0
            f0=f2
        end do
        if(i.eq.n2) then
            df(n2)=1.5d0*fn(n2)-2.d0*fn(n2-1)+0.5d0*fn(n2-2)
            ddf(n2)=2.d0*fn(n2)-5.d0*fn(n2-1)+4.d0*fn(n2-2)-fn(n2-3)
        end if
    case default
        do i=n1+dn,n2-1,dn
            f0=fn(i-1)
            f1=fn(i)
            f2=fn(i+1)
            df(i)=0.5d0*(f2-f0)
            ddf(i)=f2-2.d0*f1+f0
        end do
        if(i.eq.n2) then
            df(n2)=1.5d0*fn(n2)-2.d0*fn(n2-1)+0.5d0*fn(n2-2)
            ddf(n2)=2.d0*fn(n2)-5.d0*fn(n2-1)+4.d0*fn(n2-2)-fn(n2-3)
        end if
    end select
    
   return
   end subroutine ddiff_3points2
! === For nrc decomposion. by takto 2012/12/05 =================================
   subroutine ddiff_3points2_3D(ier,n1,n2,dn,fn,df,ddf,ista_nrc,iend_nrc,ist,ien)
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: fn(ista_nrc-6*dn:iend_nrc+6*dn)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(ista_nrc:iend_nrc)
   real(8),intent(out) :: ddf(ista_nrc:iend_nrc)
   integer,intent(in)  :: ista_nrc, iend_nrc, ist, ien

   integer:: i

   ier=0

   if(n1 >= ista_nrc .and. n1 <= iend_nrc) then
      df(n1)=-0.5d0*fn(n1+2)+2.d0*fn(n1+1)-1.5d0*fn(n1)
      ddf(n1)=-fn(n1+3)+4.d0*fn(n1+2)-5.d0*fn(n1+1)+2.d0*fn(n1)
   end if

   select case(dn)
      case(1)
         do i=max(n1+1,ist),min(n2-2,ien)
            df(i)=0.5d0*(fn(i+1)-fn(i-1))
            ddf(i)=fn(i+1)-2.d0*fn(i)+fn(i-1)
         end do
         i = (n2-2) - (n1+1) + 1
         i = (n1+1) + i
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            df(n2-1)=0.5d0*(fn(i+1)-fn(i-1))
            ddf(n2-1)=fn(i+1)-2.d0*fn(i)+fn(i-1)
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            df(n2)=1.5d0*fn(i+1)-2.d0*fn(i)+0.5d0*fn(i-1)
            ddf(n2)=2.d0*fn(i+1)-5.d0*fn(i)+4.d0*fn(i-1)-fn(n2-3)
         end if
      case(2)
         do i=max(n1+2,ist),min(n2-1,ien),2
            df(i)=0.5d0*(fn(i+1)-fn(i-1))
            ddf(i)=fn(i+1)-2.d0*fn(i)+fn(i-1)
         end do
         i = ((n2-1) - (n1+2) + 2)/2
         i = (n1+2) + 2*i
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.eq.n2) then
               df(n2)=1.5d0*fn(n2)-2.d0*fn(n2-1)+0.5d0*fn(n2-2)
               ddf(n2)=2.d0*fn(n2)-5.d0*fn(n2-1)+4.d0*fn(n2-2)-fn(n2-3)
            end if
         end if
      case default
         do i=max(n1+dn,ist),min(n2-1,ien),dn
            df(i)=0.5d0*(fn(i+1)-fn(i-1))
            ddf(i)=fn(i+1)-2.d0*fn(i)+fn(i-1)
         end do
         i = ((n2-1) - (n1+dn) + dn)/dn
         i = (n1+dn) + dn*i
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.eq.n2) then
               df(n2)=1.5d0*fn(n2)-2.d0*fn(n2-1)+0.5d0*fn(n2-2)
               ddf(n2)=2.d0*fn(n2)-5.d0*fn(n2-1)+4.d0*fn(n2-2)-fn(n2-3)
            end if
         end if
   end select

   return
   end subroutine ddiff_3points2_3D
! ==============================================================================
   
!=====================================================================
   subroutine ddiff_3points3(ier,n1,n2,dn,fn,df,ddf)                       
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: fn(*)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(*)
   real(8),intent(out) :: ddf(*)
   
   integer:: i
   real(8):: f0,f1,f2
   
   ier=0
   f0=fn(n1)
   f1=fn(n1+dn)
   f2=fn(n1+2*dn)
   df(n1)=-0.5d0*f2+2.d0*f1-1.5d0*f0
   ddf(n1)=-fn(n1+3*dn)+4.d0*f2-5.d0*f1+2.d0*f0
   do i=n1+dn,n2-2*dn,dn
        df(i)=0.5d0*(f2-f0)
        ddf(i)=f2-2.d0*f1+f0
        f0=f1
        f1=f2
        f2=fn(i+2*dn)
   end do
   df(n2-dn)=0.5d0*(f2-f0)
   ddf(n2-dn)=f2-2.d0*f1+f0
   df(n2)=1.5d0*f2-2.d0*f1+0.5d0*f0
   ddf(n2)=2.d0*f2-5.d0*f1+4.d0*f0-fn(n2-3*dn)
   return
   end subroutine ddiff_3points3
! === For nrc decomposion. by takto 2012/12/05 =================================
   subroutine ddiff_3points3_3D(ier,n1,n2,dn,fn,df,ddf,ista_nrc,iend_nrc,ist,ien)
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: fn(ista_nrc-6*dn:iend_nrc+6*dn)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(ista_nrc:iend_nrc)
   real(8),intent(out) :: ddf(ista_nrc:iend_nrc)
   integer,intent(in)  :: ista_nrc, iend_nrc, ist, ien

   integer:: i

   ier=0

   if(n1 >= ista_nrc .and. n1 <= iend_nrc) then
      df(n1)=-0.5d0*fn(n1+2*dn)+2.d0*fn(n1+dn)-1.5d0*fn(n1)
      ddf(n1)=-fn(n1+3*dn)+4.d0*fn(n1+2*dn)-5.d0*fn(n1+dn)+2.d0*fn(n1)
   end if
   do i=max(n1+dn,ist),min(n2-2*dn,ien),dn
      df(i)=0.5d0*(fn(i+dn)-fn(i-dn))
      ddf(i)=fn(i+dn)-2.d0*fn(i)+fn(i-dn)
   end do
   i = ((n2-2*dn) - (n1+dn) + dn)/dn
   i = (n1+dn) + dn*i
   if(n2-dn >= ista_nrc .and. n2-dn <= iend_nrc) then
      df(n2-dn)=0.5d0*(fn(i+dn)-fn(i-dn))
      ddf(n2-dn)=fn(i+dn)-2.d0*fn(i)+fn(i-dn)
   end if
   if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
      df(n2)=1.5d0*fn(i+dn)-2.d0*fn(i)+0.5d0*fn(i-dn)
      ddf(n2)=2.d0*fn(i+dn)-5.d0*fn(i)+4.d0*fn(i-dn)-fn(n2-3*dn)
   end if
   return
   end subroutine ddiff_3points3_3D
! ==============================================================================
   
!=====================================================================
   subroutine diff_5points(ier,n1,n2,fn,df)                       
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2
   real(8),intent(in)  :: fn(n1:n2)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(n1:n2)
   
   integer:: i
   real(8):: f0,f1,f2,f3,f4
   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d25_12=2.08333333333333d0
   real(8),parameter:: d1_12=0.08333333333333d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d2_3=0.666666666666d0
   
   ier=0
   f0=fn(n1)
   f1=fn(n1+1)
   f2=fn(n1+2)
   f3=fn(n1+3)
   f4=fn(n1+4)
   df(n1)=-0.25d0*f4+d4_3*f3-3.d0*f2+4.d0*f1-d25_12*f0
   df(n1+1)=d1_12*f4-0.5d0*f3+1.5d0*f2-d5_6*f1-0.25d0*f0
   do i=n1+2,n2-3
        df(i)=d2_3*(f3-f1)+d1_12*(f0-f4)
        f0=f1
        f1=f2
        f2=f3
        f3=f4
        f4=fn(i+3)
   end do
   df(n2-2)=d2_3*(f3-f1)+d1_12*(f0-f4)
   df(n2-1)=0.25d0*f4+d5_6*f3-1.5d0*f2+0.5d0*f1-d1_12*f0
   df(n2)=d25_12*f4-4.d0*f3+3.d0*f2-d4_3*f1+0.25d0*f0
   return
   end subroutine diff_5points
   
!=====================================================================
   subroutine diff_5points2(ier,n1,n2,dn,fn,df)                       
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: fn(*)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(*)
   
   integer:: i
   real(8):: f0,f1,f2,f3,f4
   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d25_12=2.08333333333333d0
   real(8),parameter:: d1_12=0.08333333333333d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d2_3=0.666666666666d0
   
   ier=0
   f0=fn(n1)
   f1=fn(n1+1)
   f2=fn(n1+2)
   f3=fn(n1+3)
   f4=fn(n1+4)
   df(n1)=-0.25d0*f4+d4_3*f3-3.d0*f2+4.d0*f1-d25_12*f0
   
    select case(dn)
    case(1)
        df(n1+1)=d1_12*f4-0.5d0*f3+1.5d0*f2-d5_6*f1-0.25d0*f0
        do i=n1+2,n2-3
            df(i)=d2_3*(f3-f1)+d1_12*(f0-f4)
            f0=f1
            f1=f2
            f2=f3
            f3=f4
            f4=fn(i+3)
        end do
        df(n2-2)=d2_3*(f3-f1)+d1_12*(f0-f4)
        df(n2-1)=0.25d0*f4+d5_6*f3-1.5d0*f2+0.5d0*f1-d1_12*f0
        df(n2)=d25_12*f4-4.d0*f3+3.d0*f2-d4_3*f1+0.25d0*f0
    case(2)
        f4=f2
        f3=f1
        f2=f0
        do i=n1+2,n2-2,2
            f0=f2
            f1=f3
            f2=f4
            f3=fn(i+1)
            f4=fn(i+2)
            df(i)=d2_3*(f3-f1)+d1_12*(f0-f4)
        end do
        if(i.eq.n2-1) then
            df(n2-1)=0.25d0*fn(n2)+d5_6*f4-1.5d0*f3+0.5d0*f2-d1_12*f1
        else
            df(n2)=d25_12*f4-4.d0*f3+3.d0*f2-d4_3*f1+0.25d0*f0
        end if
    case(3)
        f3=f1
        f4=f2
        do i=n1+3,n2-2,3
            f0=f3
            f1=f4
            f2=fn(i)
            f3=fn(i+1)
            f4=fn(i+2)
            df(i)=d2_3*(f3-f1)+d1_12*(f0-f4)
        end do
        if(i.eq.n2-1) then
            df(n2-1)=0.25d0*fn(n2)+d5_6*fn(n2-1)-1.5d0*f4+0.5d0*f3-d1_12*f2
        else
            df(n2)=d25_12*fn(n2)-4.d0*f4+3.d0*f3-d4_3*f2+0.25d0*f1
        end if
    case(4)
        f4=f2
        do i=n1+4,n2-2,4
            f0=f4
            f1=fn(i-1)
            f2=fn(i)
            f3=fn(i+1)
            f4=fn(i+2)
            df(i)=d2_3*(f3-f1)+d1_12*(f0-f4)
        end do
        if(i.eq.n2-1) then
            df(n2-1)=0.25d0*fn(n2)+d5_6*fn(n2-1)-1.5d0*fn(n2-2)+0.5d0*f4-d1_12*f3
        else
            df(n2)=d25_12*fn(n2)-4.d0*fn(n2-1)+3.d0*f4-d4_3*f3+0.25d0*f2
        end if
    case default
        do i=n1+dn,n2-2,dn
            df(i)=d2_3*(fn(i+1)-fn(i-1))+d1_12*(fn(i-2)-fn(i+2))
        end do
        if(i.eq.n2-1) then
            df(n2-1)=0.25d0*fn(n2)+d5_6*fn(n2-1)-1.5d0*fn(n2-2)+0.5d0*fn(n2-3)-d1_12*fn(n2-4)
        else
            df(n2)=d25_12*fn(n2)-4.d0*fn(n2-1)+3.d0*fn(n2-2)-d4_3*fn(n2-3)+0.25d0*fn(n2-4)
        end if
    end select
   return
   end subroutine diff_5points2
! === For nrc decomposion. by takto 2012/12/05 =================================
   subroutine diff_5points2_3D(ier,n1,n2,dn,fn,df,ista_nrc,iend_nrc,ist,ien)
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
!  real(8),intent(in)  :: fn(ista_nrc-3*dn:iend_nrc+3*dn)
   real(8),intent(in)  :: fn(ista_nrc-6*dn:iend_nrc+6*dn)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(ista_nrc:iend_nrc)
   integer,intent(in)  :: ista_nrc, iend_nrc, ist, ien

   integer:: i
   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d25_12=2.08333333333333d0
   real(8),parameter:: d1_12=0.08333333333333d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d2_3=0.666666666666d0

   ier=0

   if(n1 >= ista_nrc .and. n1 <= iend_nrc) then
      df(n1)=-0.25d0*fn(n1+4)+d4_3*fn(n1+3)-3.d0*fn(n1+2)+4.d0*fn(n1+1)-d25_12*fn(n1)
   end if

   select case(dn)
      case(1)
         if(n1+1 >= ista_nrc .and. n1+1 <= iend_nrc) then
            df(n1+1)=d1_12*fn(n1+4)-0.5d0*fn(n1+3)+1.5d0*fn(n1+2)-d5_6*fn(n1+1)-0.25d0*fn(n1)
         end if
         do i=max(n1+2,ist),min(n2-3,ien)
            df(i)=d2_3*(fn(i+1)-fn(i-1))+d1_12*(fn(i-2)-fn(i+2))
         end do
         i = (n2-3) - (n1+2) + 1
         i = (n1+2) + i
         if(n2-2 >= ista_nrc .and. n2-2 <= iend_nrc) then
            df(n2-2)=d2_3*(fn(i+1)-fn(i-1))+d1_12*(fn(i-2)-fn(i+2))
         end if
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            df(n2-1)=0.25d0*fn(i+2)+d5_6*fn(i+1)-1.5d0*fn(i)+0.5d0*fn(i-1)-d1_12*fn(i-2)
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            df(n2)=d25_12*fn(i+2)-4.d0*fn(i+1)+3.d0*fn(i)-d4_3*fn(i-1)+0.25d0*fn(i-2)
         end if
      case(2)
         do i=max(n1+2,ist),min(n2-2,ien),2
            df(i)=d2_3*(fn(i+1)-fn(i-1))+d1_12*(fn(i-2)-fn(i+2))
         end do
         i = ((n2-2) - (n1+2) + 2)/2
         i = (n1+2) + 2*i
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=0.25d0*fn(n2)+d5_6*fn(n2-1)-1.5d0*fn(n2-2)+0.5d0*fn(n2-3)-d1_12*fn(n2-4)
            end if
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.ne.n2-1) then
               df(n2)=d25_12*fn(i)-4.d0*fn(i-1)+3.d0*fn(i-2)-d4_3*fn(i-3)+0.25d0*fn(i-4)
            end if
         end if
      case(3)
         do i=max(n1+3,ist),min(n2-2,ien),3
            df(i)=d2_3*(fn(i+1)-fn(i-1))+d1_12*(fn(i-2)-fn(i+2))
         end do
         i = ((n2-2) - (n1+3) + 3)/3
         i = (n1+3) + 3*i
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=0.25d0*fn(n2)+d5_6*fn(n2-1)-1.5d0*fn(n2-2)+0.5d0*fn(n2-3)-d1_12*fn(n2-4)
            end if
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.ne.n2-1) then
               df(n2)=d25_12*fn(n2)-4.d0*fn(i-1)+3.d0*fn(i-2)-d4_3*fn(i-3)+0.25d0*fn(i-4)
            end if
         end if
      case(4)
         do i=max(n1+4,ist),min(n2-2,ien),4
            df(i)=d2_3*(fn(i+1)-fn(i-1))+d1_12*(fn(i-2)-fn(i+2))
         end do
         i = ((n2-2) - (n1+4) + 4)/4
         i = (n1+4) + 4*i
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=0.25d0*fn(n2)+d5_6*fn(n2-1)-1.5d0*fn(n2-2)+0.5d0*fn(n2-3)-d1_12*fn(n2-4)
            end if
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.ne.n2-1) then
               df(n2)=d25_12*fn(n2)-4.d0*fn(n2-1)+3.d0*fn(i-2)-d4_3*fn(i-3)+0.25d0*fn(i-4)
            end if
         end if
      case default
         do i=max(n1+dn,ist),min(n2-2,ien),dn
            df(i)=d2_3*(fn(i+1)-fn(i-1))+d1_12*(fn(i-2)-fn(i+2))
         end do
         i = ((n2-2) - (n1+dn) + dn)/dn
         i = (n1+dn) + dn*i
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=0.25d0*fn(n2)+d5_6*fn(n2-1)-1.5d0*fn(n2-2)+0.5d0*fn(n2-3)-d1_12*fn(n2-4)
            end if
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.ne.n2-1) then
               df(n2)=d25_12*fn(n2)-4.d0*fn(n2-1)+3.d0*fn(n2-2)-d4_3*fn(n2-3)+0.25d0*fn(n2-4)
            end if
         end if
   end select
   return
   end subroutine diff_5points2_3D
! ==============================================================================
   
!=====================================================================
   subroutine diff_5points3(ier,n1,n2,dn,fn,df)                       
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: fn(*)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(*)
   
   integer:: i
   real(8):: f0,f1,f2,f3,f4
   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d25_12=2.08333333333333d0
   real(8),parameter:: d1_12=0.08333333333333d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d2_3=0.666666666666d0
   
   ier=0
   f0=fn(n1)
   f1=fn(n1+dn)
   f2=fn(n1+2*dn)
   f3=fn(n1+3*dn)
   f4=fn(n1+4*dn)
   df(n1)=-0.25d0*f4+d4_3*f3-3.d0*f2+4.d0*f1-d25_12*f0
   df(n1+dn)=d1_12*f4-0.5d0*f3+1.5d0*f2-d5_6*f1-0.25d0*f0
   do i=n1+2*dn,n2-3*dn,dn
        df(i)=d2_3*(f3-f1)+d1_12*(f0-f4)
        f0=f1
        f1=f2
        f2=f3
        f3=f4
        f4=fn(i+3*dn)
   end do
   df(n2-2*dn)=d2_3*(f3-f1)+d1_12*(f0-f4)
   df(n2-dn)=0.25d0*f4+d5_6*f3-1.5d0*f2+0.5d0*f1-d1_12*f0
   df(n2)=d25_12*f4-4.d0*f3+3.d0*f2-d4_3*f1+0.25d0*f0
   return
   end subroutine diff_5points3
! === For nrc decomposion. by takto 2012/12/05 =================================
   subroutine diff_5points3_3D(ier,n1,n2,dn,fn,df,ista_nrc,iend_nrc,ist,ien)
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
!  real(8),intent(in)  :: fn(ista_nrc-3*dn:iend_nrc+3*dn)
   real(8),intent(in)  :: fn(ista_nrc-6*dn:iend_nrc+6*dn)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(ista_nrc:iend_nrc)
   integer,intent(in)  :: ista_nrc, iend_nrc, ist, ien

   integer:: i
   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d25_12=2.08333333333333d0
   real(8),parameter:: d1_12=0.08333333333333d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d2_3=0.666666666666d0

   ier=0
   if(n1 >= ista_nrc .and. n1 <= iend_nrc) then
      df(n1)=-0.25d0*fn(n1+4*dn)+d4_3*fn(n1+3*dn)-3.d0*fn(n1+2*dn)+4.d0*fn(n1+dn)-d25_12*fn(n1)
   end if
   if(n1+dn >= ista_nrc .and. n1+dn <= iend_nrc) then
      df(n1+dn)=d1_12*fn(n1+4*dn)-0.5d0*fn(n1+3*dn)+1.5d0*fn(n1+2*dn)-d5_6*fn(n1+dn)-0.25d0*fn(n1)
   end if
   do i=max(n1+2*dn,ist),min(n2-3*dn,ien),dn
      df(i)=d2_3*(fn(i+dn)-fn(i-dn))+d1_12*(fn(i-2*dn)-fn(i+2*dn))
   end do
   i = ((n2-3*dn) - (n1+2*dn) + dn)/dn
   i = (n1+2*dn) + dn*i
   if(n2-2*dn >= ista_nrc .and. n2-2*dn <= iend_nrc) then
      df(n2-2*dn)=d2_3*(fn(i+dn)-fn(i-dn))+d1_12*(fn(i-2*dn)-fn(i+2*dn))
   end if
   if(n2-dn >= ista_nrc .and. n2-dn <= iend_nrc) then
      df(n2-dn)=0.25d0*fn(i+2*dn)+d5_6*fn(i+dn)-1.5d0*fn(i)+0.5d0*fn(i-dn)-d1_12*fn(i-2*dn)
   end if
   if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
      df(n2)=d25_12*fn(i+2*dn)-4.d0*fn(i+dn)+3.d0*fn(i)-d4_3*fn(i-dn)+0.25d0*fn(i-2*dn)
   end if
   return
   end subroutine diff_5points3_3D
! ==============================================================================
   
!=====================================================================
   subroutine ddiff_5points(ier,n1,n2,fn,df,ddf)                       
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2
   real(8),intent(in)  :: fn(n1:n2)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(n1:n2)
   real(8),intent(out) :: ddf(n1:n2)
   
   integer:: i
   real(8):: f0,f1,f2,f3,f4
   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d25_12=2.08333333333333d0
   real(8),parameter:: d1_12=0.08333333333333d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d2_3=0.666666666666d0
   real(8),parameter:: d11_12=0.91666666666666d0
   real(8),parameter:: d14_3=4.666666666666d0
   real(8),parameter:: d26_3=8.666666666666d0
   real(8),parameter:: d35_12=2.91666666666666d0
   real(8),parameter:: d5_3=1.666666666666d0
   real(8),parameter:: d1_3=0.333333333333d0
   
   ier=0
   f0=fn(n1)
   f1=fn(n1+1)
   f2=fn(n1+2)
   f3=fn(n1+3)
   f4=fn(n1+4)
   df(n1)=-0.25d0*f4+d4_3*f3-3.d0*f2+4.d0*f1-d25_12*f0
   df(n1+1)=d1_12*f4-0.5d0*f3+1.5d0*f2-d5_6*f1-0.25d0*f0
   ddf(n1)=d11_12*f4-d14_3*f3+9.5d0*f2-d26_3*f1+d35_12*f0
   ddf(n1+1)=-d1_12*f4+d1_3*f3+0.5d0*f2-d5_3*f1+d11_12*f0
   do i=n1+2,n2-3
        df(i)=d2_3*(f3-f1)+d1_12*(f0-f4)
        ddf(i)=-d1_12*(f4+f0)+d4_3*(f3+f1)-2.5d0*f2
        f0=f1
        f1=f2
        f2=f3
        f3=f4
        f4=fn(i+3)
   end do
   df(n2-2)=d2_3*(f3-f1)+d1_12*(f0-f4)
   df(n2-1)=0.25d0*f4+d5_6*f3-1.5d0*f2+0.5d0*f1-d1_12*f0
   df(n2)=d25_12*f4-4.d0*f3+3.d0*f2-d4_3*f1+0.25d0*f0
   ddf(n2-2)=-d1_12*(f4+f0)+d4_3*(f3+f1)-2.5d0*f2
   ddf(n2-1)=d11_12*f4-d5_3*f3+0.5d0*f2+d1_3*f1-d1_12*f0
   ddf(n2)=d35_12*f4-d26_3*f3+9.5d0*f2-d14_3*f1+d11_12*f0
   return
   end subroutine ddiff_5points
   
!=====================================================================
   subroutine ddiff_5points2(ier,n1,n2,dn,fn,df,ddf)                       
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: fn(*)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(*)
   real(8),intent(out) :: ddf(*)
   
   integer:: i
   real(8):: f0,f1,f2,f3,f4
   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d25_12=2.08333333333333d0
   real(8),parameter:: d1_12=0.08333333333333d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d2_3=0.666666666666d0
   real(8),parameter:: d11_12=0.91666666666666d0
   real(8),parameter:: d14_3=4.666666666666d0
   real(8),parameter:: d26_3=8.666666666666d0
   real(8),parameter:: d35_12=2.91666666666666d0
   real(8),parameter:: d5_3=1.666666666666d0
   real(8),parameter:: d1_3=0.333333333333d0
   
   ier=0
   f0=fn(n1)
   f1=fn(n1+1)
   f2=fn(n1+2)
   f3=fn(n1+3)
   f4=fn(n1+4)
   df(n1)=-0.25d0*f4+d4_3*f3-3.d0*f2+4.d0*f1-d25_12*f0
   ddf(n1)=d11_12*f4-d14_3*f3+9.5d0*f2-d26_3*f1+d35_12*f0

    select case (dn)
    case (1)
        df(n1+1)=d1_12*f4-0.5d0*f3+1.5d0*f2-d5_6*f1-0.25d0*f0
        ddf(n1+1)=-d1_12*f4+d1_3*f3+0.5d0*f2-d5_3*f1+d11_12*f0
        do i=n1+2,n2-3
            df(i)=d2_3*(f3-f1)+d1_12*(f0-f4)
            ddf(i)=-d1_12*(f4+f0)+d4_3*(f3+f1)-2.5d0*f2
            f0=f1
            f1=f2
            f2=f3
            f3=f4
            f4=fn(i+3)
        end do
        df(n2-2)=d2_3*(f3-f1)+d1_12*(f0-f4)
        df(n2-1)=0.25d0*f4+d5_6*f3-1.5d0*f2+0.5d0*f1-d1_12*f0
        df(n2)=d25_12*f4-4.d0*f3+3.d0*f2-d4_3*f1+0.25d0*f0
        ddf(n2-2)=-d1_12*(f4+f0)+d4_3*(f3+f1)-2.5d0*f2
        ddf(n2-1)=d11_12*f4-d5_3*f3+0.5d0*f2+d1_3*f1-d1_12*f0
        ddf(n2)=d35_12*f4-d26_3*f3+9.5d0*f2-d14_3*f1+d11_12*f0
    case (2)
        f4=f2
        f3=f1
        f2=f0
        do i=n1+2,n2-2,2
            f0=f2
            f1=f3
            f2=f4
            f3=fn(i+1)
            f4=fn(i+2)
            df(i)=d2_3*(f3-f1)+d1_12*(f0-f4)
            ddf(i)=-d1_12*(f4+f0)+d4_3*(f3+f1)-2.5d0*f2
        end do
        if(i.eq.n2-1) then
            df(n2-1)=0.25d0*fn(n2)+d5_6*f4-1.5d0*f3+0.5d0*f2-d1_12*f1
            ddf(n2-1)=d11_12*fn(n2)-d5_3*f4+0.5d0*f3+d1_3*f2-d1_12*f1
        else
            df(n2)=d25_12*f4-4.d0*f3+3.d0*f2-d4_3*f1+0.25d0*f0
            ddf(n2)=d35_12*f4-d26_3*f3+9.5d0*f2-d14_3*f1+d11_12*f0
        end if
    case (3)
        f3=f1
        f4=f2
        do i=n1+3,n2-2,3
            f0=f3
            f1=f4
            f2=fn(i)
            f3=fn(i+1)
            f4=fn(i+2)
            df(i)=d2_3*(f3-f1)+d1_12*(f0-f4)
            ddf(i)=-d1_12*(f4+f0)+d4_3*(f3+f1)-2.5d0*f2
        end do
        if(i.eq.n2-1) then
            df(n2-1)=0.25d0*fn(n2)+d5_6*fn(n2-1)-1.5d0*f4+0.5d0*f3-d1_12*f2
            ddf(n2-1)=d11_12*fn(n2)-d5_3*fn(n2-1)+0.5d0*f4+d1_3*f3-d1_12*f2
        else
            df(n2)=d25_12*fn(n2)-4.d0*f4+3.d0*f3-d4_3*f2+0.25d0*f1
            ddf(n2)=d35_12*fn(n2)-d26_3*f4+9.5d0*f3-d14_3*f2+d11_12*f1
        end if
    case (4)
        f4=f2
        do i=n1+4,n2-2,4
            f0=f4
            f1=fn(i-1)
            f2=fn(i)
            f3=fn(i+1)
            f4=fn(i+2)
            df(i)=d2_3*(f3-f1)+d1_12*(f0-f4)
            ddf(i)=-d1_12*(f4+f0)+d4_3*(f3+f1)-2.5d0*f2
        end do
        if(i.eq.n2-1) then
            df(n2-1)=0.25d0*fn(n2)+d5_6*fn(n2-1)-1.5d0*fn(n2-2)+0.5d0*f4-d1_12*f3
            ddf(n2-1)=d11_12*fn(n2)-d5_3*fn(n2-1)+0.5d0*fn(n2-2)+d1_3*f4-d1_12*f3
        else
            df(n2)=d25_12*fn(n2)-4.d0*fn(n2-1)+3.d0*f4-d4_3*f3+0.25d0*f2
            ddf(n2)=d35_12*fn(n2)-d26_3*fn(n2-1)+9.5d0*f4-d14_3*f3+d11_12*f2
        end if
    case default
        do i=n1+dn,n2-2,dn
            f0=fn(i-2)
            f1=fn(i-1)
            f2=fn(i)
            f3=fn(i+1)
            f4=fn(i+2)
            df(i)=d2_3*(f3-f1)+d1_12*(f0-f4)
            ddf(i)=-d1_12*(f4+f0)+d4_3*(f3+f1)-2.5d0*f2
        end do
        if(i.eq.n2-1) then
            df(n2-1)=0.25d0*fn(n2)+d5_6*fn(n2-1)-1.5d0*fn(n2-2)+0.5d0*fn(n2-3)-d1_12*fn(n2-4)
            ddf(n2-1)=d11_12*fn(n2)-d5_3*fn(n2-1)+0.5d0*fn(n2-2)+d1_3*fn(n2-3)-d1_12*fn(n2-4)
        else
            df(n2)=d25_12*fn(n2)-4.d0*fn(n2-1)+3.d0*fn(n2-2)-d4_3*fn(n2-3)+0.25d0*fn(n2-4)
            ddf(n2)=d35_12*fn(n2)-d26_3*fn(n2-1)+9.5d0*fn(n2-2)-d14_3*fn(n2-3)+d11_12*fn(n2-4)
        end if
    end select
    
   return
   end subroutine ddiff_5points2
! === For nrc decomposion. by takto 2012/12/05 =================================
   subroutine ddiff_5points2_3D(ier,n1,n2,dn,fn,df,ddf,ista_nrc,iend_nrc,ist,ien)
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: fn(ista_nrc-6*dn:iend_nrc+6*dn)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(ista_nrc:iend_nrc)
   real(8),intent(out) :: ddf(ista_nrc:iend_nrc)
   integer,intent(in)  :: ista_nrc, iend_nrc, ist, ien

   integer:: i
   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d25_12=2.08333333333333d0
   real(8),parameter:: d1_12=0.08333333333333d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d2_3=0.666666666666d0
   real(8),parameter:: d11_12=0.91666666666666d0
   real(8),parameter:: d14_3=4.666666666666d0
   real(8),parameter:: d26_3=8.666666666666d0
   real(8),parameter:: d35_12=2.91666666666666d0
   real(8),parameter:: d5_3=1.666666666666d0
   real(8),parameter:: d1_3=0.333333333333d0
   
   ier=0

   if(n1 >= ista_nrc .and. n1 <= iend_nrc) then
      df(n1)=-0.25d0*fn(n1+4)+d4_3*fn(n1+3)-3.d0*fn(n1+2)+4.d0*fn(n1+1)-d25_12*fn(n1)
      ddf(n1)=d11_12*fn(n1+4)-d14_3*fn(n1+3)+9.5d0*fn(n1+2)-d26_3*fn(n1+1)+d35_12*fn(n1)
   end if

   select case (dn)
      case (1)
         if(n1+1 >= ista_nrc .and. n1+1 <= iend_nrc) then
            df(n1+1)=d1_12*fn(n1+4)-0.5d0*fn(n1+3)+1.5d0*fn(n1+2)-d5_6*fn(n1+1)-0.25d0*fn(n1)
            ddf(n1+1)=-d1_12*fn(n1+4)+d1_3*fn(n1+3)+0.5d0*fn(n1+2)-d5_3*fn(n1+1)+d11_12*fn(n1)
         end if
         do i=max(n1+2,ist),min(n2-3,ien)
            df(i)=d2_3*(fn(i+1)-fn(i-1))+d1_12*(fn(i-2)-fn(i+2))
            ddf(i)=-d1_12*(fn(i+2)+fn(i-2))+d4_3*(fn(i+1)+fn(i-1))-2.5d0*fn(i)
         end do
         i = (n2-3) - (n1+2) + 1
         i = (n1+2) + i
         if(n2-2 >= ista_nrc .and. n2-2 <= iend_nrc) then
            df(n2-2)=d2_3*(fn(i+1)-fn(i-1))+d1_12*(fn(i-2)-fn(i+2))
            ddf(n2-2)=-d1_12*(fn(i+2)+fn(i-2))+d4_3*(fn(i+1)+fn(i-1))-2.5d0*fn(i)
         end if
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            df(n2-1)=0.25d0*fn(i+2)+d5_6*fn(i+1)-1.5d0*fn(i)+0.5d0*fn(i-1)-d1_12*fn(i-2)
            ddf(n2-1)=d11_12*fn(i+2)-d5_3*fn(i+1)+0.5d0*fn(i)+d1_3*fn(i-1)-d1_12*fn(i-2)
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            df(n2)=d25_12*fn(i+2)-4.d0*fn(i+1)+3.d0*fn(i)-d4_3*fn(i-1)+0.25d0*fn(i-2)
            ddf(n2)=d35_12*fn(i+2)-d26_3*fn(i+1)+9.5d0*fn(i)-d14_3*fn(i-1)+d11_12*fn(i-2)
         end if
      case (2)
         do i=max(n1+2,ist),min(n2-2,ien),2
            df(i)=d2_3*(fn(i+1)-fn(i-1))+d1_12*(fn(i-2)-fn(i+2))
            ddf(i)=-d1_12*(fn(i+2)+fn(i-2))+d4_3*(fn(i+1)+fn(i-1))-2.5d0*fn(i)
         end do
         i = ((n2-2) - (n1+2) + 2)/2
         i = (n1+2) + 2*i
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=0.25d0*fn(n2)+d5_6*fn(i)-1.5d0*fn(i-1)+0.5d0*fn(i-2)-d1_12*fn(i-3)
               ddf(n2-1)=d11_12*fn(n2)-d5_3*fn(i)+0.5d0*fn(i-1)+d1_3*fn(i-2)-d1_12*fn(i-3)
            end if
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.ne.n2-1) then
               df(n2)=d25_12*fn(i)-4.d0*fn(i-1)+3.d0*fn(i-2)-d4_3*fn(i-3)+0.25d0*fn(i-4)
               ddf(n2)=d35_12*fn(i)-d26_3*fn(i-1)+9.5d0*fn(i-2)-d14_3*fn(i-3)+d11_12*fn(i-4)
            end if
         end if
      case (3)
         do i=max(n1+3,ist),min(n2-2,ien),3
            df(i)=d2_3*(fn(i+1)-fn(i-1))+d1_12*(fn(i-2)-fn(i+2))
            ddf(i)=-d1_12*(fn(i+2)+fn(i-2))+d4_3*(fn(i+1)+fn(i-1))-2.5d0*fn(i)
         end do
         i = ((n2-2) - (n1+3) + 3)/3
         i = (n1+3) + 3*i
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=0.25d0*fn(n2)+d5_6*fn(n2-1)-1.5d0*fn(i-1)+0.5d0*fn(i-2)-d1_12*fn(i-3)
               ddf(n2-1)=d11_12*fn(n2)-d5_3*fn(n2-1)+0.5d0*fn(i-1)+d1_3*fn(i-2)-d1_12*fn(i-3)
            end if
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.ne.n2-1) then
               df(n2)=d25_12*fn(n2)-4.d0*fn(i-1)+3.d0*fn(i-2)-d4_3*fn(i-3)+0.25d0*fn(i-4)
               ddf(n2)=d35_12*fn(n2)-d26_3*fn(i-1)+9.5d0*fn(i-2)-d14_3*fn(i-3)+d11_12*fn(i-4)
            end if
         end if
      case (4)
         do i=max(n1+4,ist),min(n2-2,ien),4
            df(i)=d2_3*(fn(i+1)-fn(i-1))+d1_12*(fn(i-2)-fn(i+2))
            ddf(i)=-d1_12*(fn(i+2)+fn(i-2))+d4_3*(fn(i+1)+fn(i-1))-2.5d0*fn(i)
         end do
         i = ((n2-2) - (n1+4) + 4)/4
         i = (n1+4) + 4*i
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=0.25d0*fn(n2)+d5_6*fn(n2-1)-1.5d0*fn(n2-2)+0.5d0*fn(i-2)-d1_12*fn(i-3)
               ddf(n2-1)=d11_12*fn(n2)-d5_3*fn(n2-1)+0.5d0*fn(n2-2)+d1_3*fn(i-2)-d1_12*fn(i-3)
            end if
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.ne.n2-1) then
               df(n2)=d25_12*fn(n2)-4.d0*fn(n2-1)+3.d0*fn(i-2)-d4_3*fn(i-3)+0.25d0*fn(i-4)
               ddf(n2)=d35_12*fn(n2)-d26_3*fn(n2-1)+9.5d0*fn(i-2)-d14_3*fn(i-3)+d11_12*fn(i-4)
            end if
         end if
      case default
         do i=max(n1+dn,ist),min(n2-2,ien),dn
            df(i)=d2_3*(fn(i+1)-fn(i-1))+d1_12*(fn(i-2)-fn(i+2))
            ddf(i)=-d1_12*(fn(i+2)+fn(i-2))+d4_3*(fn(i+1)+fn(i-1))-2.5d0*fn(i)
         end do
         i = ((n2-2) - (n1+dn) + dn)/dn
         i = (n1+dn) + dn*i
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=0.25d0*fn(n2)+d5_6*fn(n2-1)-1.5d0*fn(n2-2)+0.5d0*fn(n2-3)-d1_12*fn(n2-4)
               ddf(n2-1)=d11_12*fn(n2)-d5_3*fn(n2-1)+0.5d0*fn(n2-2)+d1_3*fn(n2-3)-d1_12*fn(n2-4)
            end if
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.ne.n2-1) then
               df(n2)=d25_12*fn(n2)-4.d0*fn(n2-1)+3.d0*fn(n2-2)-d4_3*fn(n2-3)+0.25d0*fn(n2-4)
               ddf(n2)=d35_12*fn(n2)-d26_3*fn(n2-1)+9.5d0*fn(n2-2)-d14_3*fn(n2-3)+d11_12*fn(n2-4)
            end if
         end if
   end select

   return
   end subroutine ddiff_5points2_3D
! ==============================================================================
   
!=====================================================================
   subroutine ddiff_5points3(ier,n1,n2,dn,fn,df,ddf)                       
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: fn(*)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(*)
   real(8),intent(out) :: ddf(*)
   
   integer:: i
   real(8):: f0,f1,f2,f3,f4
   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d25_12=2.08333333333333d0
   real(8),parameter:: d1_12=0.08333333333333d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d2_3=0.666666666666d0
   real(8),parameter:: d11_12=0.91666666666666d0
   real(8),parameter:: d14_3=4.666666666666d0
   real(8),parameter:: d26_3=8.666666666666d0
   real(8),parameter:: d35_12=2.91666666666666d0
   real(8),parameter:: d5_3=1.666666666666d0
   real(8),parameter:: d1_3=0.333333333333d0
   
   ier=0
   f0=fn(n1)
   f1=fn(n1+dn)
   f2=fn(n1+2*dn)
   f3=fn(n1+3*dn)
   f4=fn(n1+4*dn)
   df(n1)=-0.25d0*f4+d4_3*f3-3.d0*f2+4.d0*f1-d25_12*f0
   df(n1+dn)=d1_12*f4-0.5d0*f3+1.5d0*f2-d5_6*f1-0.25d0*f0
   ddf(n1)=d11_12*f4-d14_3*f3+9.5d0*f2-d26_3*f1+d35_12*f0
   ddf(n1+dn)=-d1_12*f4+d1_3*f3+0.5d0*f2-d5_3*f1+d11_12*f0
   do i=n1+2*dn,n2-3*dn,dn
        df(i)=d2_3*(f3-f1)+d1_12*(f0-f4)
        ddf(i)=-d1_12*(f4+f0)+d4_3*(f3+f1)-2.5d0*f2
        f0=f1
        f1=f2
        f2=f3
        f3=f4
        f4=fn(i+3*dn)
   end do
   df(n2-2*dn)=d2_3*(f3-f1)+d1_12*(f0-f4)
   df(n2-dn)=0.25d0*f4+d5_6*f3-1.5d0*f2+0.5d0*f1-d1_12*f0
   df(n2)=d25_12*f4-4.d0*f3+3.d0*f2-d4_3*f1+0.25d0*f0
   ddf(n2-2*dn)=-d1_12*(f4+f0)+d4_3*(f3+f1)-2.5d0*f2
   ddf(n2-dn)=d11_12*f4-d5_3*f3+0.5d0*f2+d1_3*f1-d1_12*f0
   ddf(n2)=d35_12*f4-d26_3*f3+9.5d0*f2-d14_3*f1+d11_12*f0
   return
   end subroutine ddiff_5points3
! === For nrc decomposion. by takto 2012/12/05 =================================
   subroutine ddiff_5points3_3D(ier,n1,n2,dn,fn,df,ddf,ista_nrc,iend_nrc,ist,ien)
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: fn(ista_nrc-6*dn:iend_nrc+6*dn)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(ista_nrc:iend_nrc)
   real(8),intent(out) :: ddf(ista_nrc:iend_nrc)
   integer,intent(in)  :: ista_nrc, iend_nrc, ist, ien

   integer:: i
   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d25_12=2.08333333333333d0
   real(8),parameter:: d1_12=0.08333333333333d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d2_3=0.666666666666d0
   real(8),parameter:: d11_12=0.91666666666666d0
   real(8),parameter:: d14_3=4.666666666666d0
   real(8),parameter:: d26_3=8.666666666666d0
   real(8),parameter:: d35_12=2.91666666666666d0
   real(8),parameter:: d5_3=1.666666666666d0
   real(8),parameter:: d1_3=0.333333333333d0
   
   ier=0
   if(n1 >= ista_nrc .and. n1 <= iend_nrc) then
      df(n1)=-0.25d0*fn(n1+4*dn)+d4_3*fn(n1+3*dn)-3.d0*fn(n1+2*dn)+4.d0*fn(n1+dn)-d25_12*fn(n1)
      ddf(n1)=d11_12*fn(n1+4*dn)-d14_3*fn(n1+3*dn)+9.5d0*fn(n1+2*dn)-d26_3*fn(n1+dn)+d35_12*fn(n1)
   end if
   if(n1+dn >= ista_nrc .and. n1+dn <= iend_nrc) then
      df(n1+dn)=d1_12*fn(n1+4*dn)-0.5d0*fn(n1+3*dn)+1.5d0*fn(n1+2*dn)-d5_6*fn(n1+dn)-0.25d0*fn(n1)
      ddf(n1+dn)=-d1_12*fn(n1+4*dn)+d1_3*fn(n1+3*dn)+0.5d0*fn(n1+2*dn)-d5_3*fn(n1+dn)+d11_12*fn(n1)
   end if
   do i=max(n1+2*dn,ist),min(n2-3*dn,ien),dn
      df(i)=d2_3*(fn(i+dn)-fn(i-dn))+d1_12*(fn(i-2*dn)-fn(i+2*dn))
      ddf(i)=-d1_12*(fn(i+2*dn)+fn(i-2*dn))+d4_3*(fn(i+dn)+fn(i-dn))-2.5d0*fn(i)
   end do
   i = ((n2-3*dn) - (n1+2*dn) + dn)/dn
   i = (n1+2*dn) + dn*i
   if(n2-2*dn >= ista_nrc .and. n2-2*dn <= iend_nrc) then
      df(n2-2*dn)=d2_3*(fn(i+dn)-fn(i-dn))+d1_12*(fn(i-2*dn)-fn(i+2*dn))
      ddf(n2-2*dn)=-d1_12*(fn(i+2*dn)+fn(i-2*dn))+d4_3*(fn(i+dn)+fn(i-dn))-2.5d0*fn(i)
   end if
   if(n2-dn >= ista_nrc .and. n2-dn <= iend_nrc) then
      df(n2-dn)=0.25d0*fn(i+2*dn)+d5_6*fn(i+dn)-1.5d0*fn(i)+0.5d0*fn(i-dn)-d1_12*fn(i-2*dn)
      ddf(n2-dn)=d11_12*fn(i+2*dn)-d5_3*fn(i+dn)+0.5d0*fn(i)+d1_3*fn(i-dn)-d1_12*fn(i-2*dn)
   end if
   if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
      df(n2)=d25_12*fn(i+2*dn)-4.d0*fn(i+dn)+3.d0*fn(i)-d4_3*fn(i-dn)+0.25d0*fn(i-2*dn)
      ddf(n2)=d35_12*fn(i+2*dn)-d26_3*fn(i+dn)+9.5d0*fn(i)-d14_3*fn(i-dn)+d11_12*fn(i-2*dn)
   end if
   return
   end subroutine ddiff_5points3_3D
! ==============================================================================
   
!=====================================================================
   subroutine diff_7points(ier,n1,n2,fn,df)                       
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2
   real(8),intent(in)  :: fn(n1:n2)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(n1:n2)
   
   integer:: i
   real(8):: f0,f1,f2,f3,f4,f5,f6
   
   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d5_3=1.666666666666d0
   real(8),parameter:: d20_3=6.666666666666d0
   real(8),parameter:: d1_6=0.1666666666666d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d7_12=0.58333333333333d0
   real(8),parameter:: d2_15=0.1333333333333d0
   real(8),parameter:: d1_30=0.0333333333333d0
   real(8),parameter:: d1_60=0.01666666666666d0
   real(8),parameter:: d77_60=1.28333333333333d0
   real(8),parameter:: d147_60=2.45d0
   
   ier=0
   f0=fn(n1)
   f1=fn(n1+1)
   f2=fn(n1+2)
   f3=fn(n1+3)
   f4=fn(n1+4)
   f5=fn(n1+5)
   f6=fn(n1+6)
   df(n1)=-d1_6*f6+1.2d0*f5-3.75d0*f4+d20_3*f3-7.5d0*f2+6.d0*f1-d147_60*f0
   df(n1+1)=d1_30*f6-0.25d0*f5+d5_6*f4-d5_3*f3+2.5d0*f2-d77_60*f1-d1_6*f0
   df(n1+2)=-d1_60*f6+d2_15*f5-0.5d0*f4+d4_3*f3-d7_12*f2-0.4d0*f1+d1_30*f0
   do i=n1+3,n2-4
        df(i)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
        f0=f1
        f1=f2
        f2=f3
        f3=f4
        f4=f5
        f5=f6
        f6=fn(i+4)
   end do
   df(n2-3)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
   df(n2-2)=-d1_30*f6+0.4d0*f5+d7_12*f4-d4_3*f3+0.5d0*f2-d2_15*f1+d1_60*f0
   df(n2-1)=d1_6*f6+d77_60*f5-2.5d0*f4+d5_3*f3-d5_6*f2+0.25d0*f1-d1_30*f0
   df(n2)=d147_60*f6-6.d0*f5+7.5d0*f4-d20_3*f3+3.75d0*f2-1.2d0*f1+d1_6*f0
   return
   end subroutine diff_7points
   
!=====================================================================
   subroutine diff_7points2(ier,n1,n2,dn,fn,df)                       
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: fn(*)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(*)
   
   integer:: i
   real(8):: f0,f1,f2,f3,f4,f5,f6
   
   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d5_3=1.666666666666d0
   real(8),parameter:: d20_3=6.666666666666d0
   real(8),parameter:: d1_6=0.1666666666666d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d7_12=0.58333333333333d0
   real(8),parameter:: d2_15=0.1333333333333d0
   real(8),parameter:: d1_30=0.0333333333333d0
   real(8),parameter:: d1_60=0.01666666666666d0
   real(8),parameter:: d77_60=1.28333333333333d0
   real(8),parameter:: d147_60=2.45d0
   
   ier=0
   f0=fn(n1)
   f1=fn(n1+1)
   f2=fn(n1+2)
   f3=fn(n1+3)
   f4=fn(n1+4)
   f5=fn(n1+5)
   f6=fn(n1+6)
   df(n1)=-d1_6*f6+1.2d0*f5-3.75d0*f4+d20_3*f3-7.5d0*f2+6.d0*f1-d147_60*f0
   
    select case (dn)
    case (1)
        df(n1+1)=d1_30*f6-0.25d0*f5+d5_6*f4-d5_3*f3+2.5d0*f2-d77_60*f1-d1_6*f0
        df(n1+2)=-d1_60*f6+d2_15*f5-0.5d0*f4+d4_3*f3-d7_12*f2-0.4d0*f1+d1_30*f0
        do i=n1+3,n2-4
            df(i)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
            f0=f1
            f1=f2
            f2=f3
            f3=f4
            f4=f5
            f5=f6
            f6=fn(i+4)
        end do
        df(n2-3)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
        df(n2-2)=-d1_30*f6+0.4d0*f5+d7_12*f4-d4_3*f3+0.5d0*f2-d2_15*f1+d1_60*f0
        df(n2-1)=d1_6*f6+d77_60*f5-2.5d0*f4+d5_3*f3-d5_6*f2+0.25d0*f1-d1_30*f0
        df(n2)=d147_60*f6-6.d0*f5+7.5d0*f4-d20_3*f3+3.75d0*f2-1.2d0*f1+d1_6*f0
    case (2)
        df(n1+2)=-d1_60*f6+d2_15*f5-0.5d0*f4+d4_3*f3-d7_12*f2-0.4d0*f1+d1_30*f0
        f6=f5
        f5=f4
        f4=f3
        f3=f2
        f2=f1
        do i=n1+4,n2-3,2
            f0=f2
            f1=f3
            f2=f4
            f3=f5
            f4=f6
            f5=fn(i+2)
            f6=fn(i+3)
            df(i)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
        end do
        if(i.eq.n2-2) then
            df(n2-2)=-d1_30*fn(n2)+0.4d0*f6+d7_12*f5-d4_3*f4+0.5d0*f3-d2_15*f2+d1_60*f1
            df(n2)=d147_60*fn(n2)-6.d0*f6+7.5d0*f5-d20_3*f4+3.75d0*f3-1.2d0*f2+d1_6*f1
        else if(i.eq.n2-1) then
            df(n2-1)=d1_6*f6+d77_60*f5-2.5d0*f4+d5_3*f3-d5_6*f2+0.25d0*f1-d1_30*f0
        end if
    case (3)
        f6=f3
        f5=f2
        f4=f1
        f3=f0
        do i=n1+3,n2-3,3
            f0=f3
            f1=f4
            f2=f5
            f3=f6
            f4=fn(i+1)
            f5=fn(i+2)
            f6=fn(i+3)
            df(i)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
        end do
        if(i.eq.n2-2) then
            df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*f6-d4_3*f5+0.5d0*f4-d2_15*f3+d1_60*f2
        else if(i.eq.n2-1) then
            df(n2-1)=d1_6*fn(n2)+d77_60*f6-2.5d0*f5+d5_3*f4-d5_6*f3+0.25d0*f2-d1_30*f1
        else if(i.eq.n2) then
            df(n2)=d147_60*f6-6.d0*f5+7.5d0*f4-d20_3*f3+3.75d0*f2-1.2d0*f1+d1_6*f0
        end if
    case (4)
        f6=f3
        f5=f2
        f4=f1
        do i=n1+4,n2-3,4
            f0=f4
            f1=f5
            f2=f6
            f3=fn(i)
            f4=fn(i+1)
            f5=fn(i+2)
            f6=fn(i+3)
            df(i)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
        end do
        if(i.eq.n2-2) then
            df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*fn(n2-2)-d4_3*f6+0.5d0*f5-d2_15*f4+d1_60*f3
        else if(i.eq.n2-1) then
            df(n2-1)=d1_6*fn(n2)+d77_60*fn(n2-1)-2.5d0*f6+d5_3*f5-d5_6*f4+0.25d0*f3-d1_30*f2
        else if(i.eq.n2) then
            df(n2)=d147_60*fn(n2)-6.d0*f6+7.5d0*f5-d20_3*f4+3.75d0*f3-1.2d0*f2+d1_6*f1
        end if
    case (5)
        f6=f3
        f5=f2
        do i=n1+5,n2-3,5
            f0=f5
            f1=f6
            f2=fn(i-1)
            f3=fn(i)
            f4=fn(i+1)
            f5=fn(i+2)
            f6=fn(i+3)
            df(i)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
        end do
        if(i.eq.n2-2) then
            df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*fn(n2-2)-d4_3*fn(n2-3)+0.5d0*f6-d2_15*f5+d1_60*f4
        else if(i.eq.n2-1) then
            df(n2-1)=d1_6*fn(n2)+d77_60*fn(n2-1)-2.5d0*fn(n2-2)+d5_3*f6-d5_6*f5+0.25d0*f4-d1_30*f3
        else if(i.eq.n2) then
            df(n2)=d147_60*fn(n2)-6.d0*fn(n2-1)+7.5d0*f6-d20_3*f5+3.75d0*f4-1.2d0*f3+d1_6*f2
        end if
    case (6)
        f6=f3
        do i=n1+6,n2-3,6
            f0=f6
            f1=fn(i-2)
            f2=fn(i-1)
            f3=fn(i)
            f4=fn(i+1)
            f5=fn(i+2)
            f6=fn(i+3)
            df(i)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
        end do
        if(i.eq.n2-2) then
            df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*fn(n2-2)-d4_3*fn(n2-3)+0.5d0*fn(n2-4)-d2_15*f6+d1_60*f5
        else if(i.eq.n2-1) then
            df(n2-1)=d1_6*fn(n2)+d77_60*fn(n2-1)-2.5d0*fn(n2-2)+d5_3*fn(n2-3)-d5_6*f6+0.25d0*f5-d1_30*f4
        else if(i.eq.n2) then
            df(n2)=d147_60*fn(n2)-6.d0*fn(n2-1)+7.5d0*fn(n2-2)-d20_3*f6+3.75d0*f5-1.2d0*f4+d1_6*f3
        end if
    case default
        do i=n1+dn,n2-3,dn
            df(i)=d1_60*(fn(i+3)-fn(i-3)) &
                            +0.15d0*(fn(i-2)-fn(i+2)) &
                            +0.75d0*(fn(i+1)-fn(i-1))
        end do
        if(i.eq.n2-2) then
            df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*fn(n2-2)-d4_3*fn(n2-3)&
            +0.5d0*fn(n2-4)-d2_15*fn(n2-5)+d1_60*fn(n2-6)
        else if(i.eq.n2-1) then
            df(n2-1)=d1_6*fn(n2)+d77_60*fn(n2-1)-2.5d0*fn(n2-2)+d5_3*fn(n2-3)&
            -d5_6*fn(n2-4)+0.25d0*fn(n2-5)-d1_30*fn(n2-6)
        else if(i.eq.n2) then
            df(n2)=d147_60*fn(n2)-6.d0*fn(n2-1)+7.5d0*fn(n2-2)-d20_3*fn(n2-3)&
            +3.75d0*fn(n2-4)-1.2d0*fn(n2-5)+d1_6*fn(n2-6)
        end if
    end select
   return
   end subroutine diff_7points2
! === For nrc decomposion. by takto 2012/12/05 =================================
   subroutine diff_7points2_3D(ier,n1,n2,dn,fn,df,ista_nrc,iend_nrc,ist,ien)
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
!  real(8),intent(in)  :: fn(ista_nrc-3*dn:iend_nrc+3*dn)
   real(8),intent(in)  :: fn(ista_nrc-6*dn:iend_nrc+6*dn)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(ista_nrc:iend_nrc)
   integer,intent(in)  :: ista_nrc, iend_nrc, ist, ien

   integer:: i

   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d5_3=1.666666666666d0
   real(8),parameter:: d20_3=6.666666666666d0
   real(8),parameter:: d1_6=0.1666666666666d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d7_12=0.58333333333333d0
   real(8),parameter:: d2_15=0.1333333333333d0
   real(8),parameter:: d1_30=0.0333333333333d0
   real(8),parameter:: d1_60=0.01666666666666d0
   real(8),parameter:: d77_60=1.28333333333333d0
   real(8),parameter:: d147_60=2.45d0

   ier=0

   if(n1 >= ista_nrc .and. n1 <= iend_nrc) then
      df(n1)=-d1_6*fn(n1+6)+1.2d0*fn(n1+5)-3.75d0*fn(n1+4)+d20_3*fn(n1+3)-7.5d0*fn(n1+2)+6.d0*fn(n1+1)-d147_60*fn(n1)
   end if

   select case (dn)
      case (1)
         if(n1+1 >= ista_nrc .and. n1+1 <= iend_nrc) then
            df(n1+1)=d1_30*fn(n1+6)-0.25d0*fn(n1+5)+d5_6*fn(n1+4) &
                  & -d5_3*fn(n1+3)+2.5d0*fn(n1+2)-d77_60*fn(n1+1)-d1_6*fn(n1)
         end if
         if(n1+2 >= ista_nrc .and. n1+2 <= iend_nrc) then
            df(n1+2)=-d1_60*fn(n1+6)+d2_15*fn(n1+5)-0.5d0*fn(n1+4)+d4_3*fn(n1+3) &
                  & -d7_12*fn(n1+2)-0.4d0*fn(n1+1)+d1_30*fn(n1)
         end if
         do i=max(n1+3,ist),min(n2-4,ien)
            df(i)=d1_60*(fn(i+3)-fn(i-3))+0.15d0*(fn(i-2)-fn(i+2))+0.75d0*(fn(i+1)-fn(i-1))
         end do
         i = (n2-4) - (n1+3) + 1
         i = (n1+3) + i
         if(n2-3 >= ista_nrc .and. n2-3 <= iend_nrc) then
            df(n2-3)=d1_60*(fn(i+3)-fn(i-3))+0.15d0*(fn(i-2)-fn(i+2))+0.75d0*(fn(i+1)-fn(i-1))
         end if
         if(n2-2 >= ista_nrc .and. n2-2 <= iend_nrc) then
            df(n2-2)=-d1_30*fn(i+3)+0.4d0*fn(i+2)+d7_12*fn(i+1)-d4_3*fn(i)+0.5d0*fn(i-1)-d2_15*fn(i-2)+d1_60*fn(i-3)
         end if
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            df(n2-1)=d1_6*fn(i+3)+d77_60*fn(i+2)-2.5d0*fn(i+1)+d5_3*fn(i)-d5_6*fn(i-1)+0.25d0*fn(i-2)-d1_30*fn(i-3)
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            df(n2)=d147_60*fn(i+3)-6.d0*fn(i+2)+7.5d0*fn(i+1)-d20_3*fn(i)+3.75d0*fn(i-1)-1.2d0*fn(i-2)+d1_6*fn(i-3)
         end if
      case (2)
         if(n1+2 >= ista_nrc .and. n1+2 <= iend_nrc) then
            df(n1+2)=-d1_60*fn(n1+6)+d2_15*fn(n1+5)-0.5d0*fn(n1+4) &
                   & +d4_3*fn(n1+3)-d7_12*fn(n1+2)-0.4d0*fn(n1+1)+d1_30*fn(n1)
         end if
         do i=max(n1+4,ist),min(n2-3,ien),2
            df(i)=d1_60*(fn(i+3)-fn(i-3))+0.15d0*(fn(i-2)-fn(i+2))+0.75d0*(fn(i+1)-fn(i-1))
         end do
         i = ((n2-3) - (n1+4) + 2)/2
         i = (n1+4) + 2*i
         if(n2-2 >= ista_nrc .and. n2-2 <= iend_nrc) then
            if(i.eq.n2-2) then
               df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*fn(n2-2)-d4_3*fn(n2-3)+0.5d0*fn(n2-4)-d2_15*fn(n2-5)+d1_60*fn(n2-6)
               df(n2)=d147_60*fn(n2)-6.d0*fn(n2-1)+7.5d0*fn(n2-2)-d20_3*fn(n2-3)+3.75d0*fn(n2-4)-1.2d0*fn(n2-5)+d1_6*fn(n2-6)
            end if
         end if
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=d1_6*fn(n2)+d77_60*fn(n2-1)-2.5d0*fn(n2-1)+d5_3*fn(n2-2)-d5_6*fn(n2-3)+0.25d0*fn(n2-4)-d1_30*fn(n2-5)
            end if
         end if
      case (3)
         do i=max(n1+3,ist),min(n2-3,ien),3
            df(i)=d1_60*(fn(i+3)-fn(i-3))+0.15d0*(fn(i-2)-fn(i+2))+0.75d0*(fn(i+1)-fn(i-1))
         end do
         i = ((n2-3) - (n1+3) + 3)/3
         i = (n1+3) + 3*i
         if(n2-2 >= ista_nrc .and. n2-2 <= iend_nrc) then
            if(i.eq.n2-2) then
               df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*fn(n2-2)-d4_3*fn(n2-3)+0.5d0*fn(n2-4)-d2_15*fn(n2-5)+d1_60*fn(n2-6)
            end if
         end if
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=d1_6*fn(n2)+d77_60*fn(n2-1)-2.5d0*fn(n2-2)+d5_3*fn(n2-3)-d5_6*fn(n2-4)+0.25d0*fn(n2-5)-d1_30*fn(n2-6)
            end if
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.eq.n2) then
               df(n2)=d147_60*fn(n2)-6.d0*fn(n2-1)+7.5d0*fn(n2-2)-d20_3*fn(n2-3)+3.75d0*fn(n2-4)-1.2d0*fn(n2-5)+d1_6*fn(n2-6)
            end if
         end if
      case (4)
         do i=max(n1+4,ist),min(n2-3,ien),4
            df(i)=d1_60*(fn(i+3)-fn(i-3))+0.15d0*(fn(i-2)-fn(i+2))+0.75d0*(fn(i+1)-fn(i-1))
         end do
         i = ((n2-3) - (n1+4) + 4)/4
         i = (n1+4) + 4*i
         if(n2-2 >= ista_nrc .and. n2-2 <= iend_nrc) then
            if(i.eq.n2-2) then
               df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*fn(n2-2)-d4_3*fn(n2-3)+0.5d0*fn(n2-4)-d2_15*fn(n2-5)+d1_60*fn(n2-6)
            end if
         end if
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=d1_6*fn(n2)+d77_60*fn(n2-1)-2.5d0*fn(n2-2)+d5_3*fn(n2-3)-d5_6*fn(n2-4)+0.25d0*fn(n2-5)-d1_30*fn(n2-6)
            end if
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.eq.n2) then
               df(n2)=d147_60*fn(n2)-6.d0*fn(n2-1)+7.5d0*fn(n2-2)-d20_3*fn(n2-3)+3.75d0*fn(n2-4)-1.2d0*fn(n2-5)+d1_6*fn(n2-6)
            end if
         end if
      case (5)
         do i=max(n1+5,ist),min(n2-3,ien),5
            df(i)=d1_60*(fn(i+3)-fn(i-3))+0.15d0*(fn(i-2)-fn(i+2))+0.75d0*(fn(i+1)-fn(i-1))
         end do
         i = ((n2-3) - (n1+5) + 5)/5
         i = (n1+5) + 5*i
         if(n2-2 >= ista_nrc .and. n2-2 <= iend_nrc) then
            if(i.eq.n2-2) then
               df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*fn(n2-2)-d4_3*fn(n2-3)+0.5d0*fn(n2-4)-d2_15*fn(n2-5)+d1_60*fn(n2-6)
            end if
         end if
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=d1_6*fn(n2)+d77_60*fn(n2-1)-2.5d0*fn(n2-2)+d5_3*fn(n2-3)-d5_6*fn(n2-4)+0.25d0*fn(n2-5)-d1_30*fn(n2-6)
            end if
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.eq.n2) then
               df(n2)=d147_60*fn(n2)-6.d0*fn(n2-1)+7.5d0*fn(n2-2)-d20_3*fn(n2-3)+3.75d0*fn(n2-4)-1.2d0*fn(n2-5)+d1_6*fn(n2-6)
            end if
         end if
      case (6)
         do i=max(n1+6,ist),min(n2-3,ien),6
            df(i)=d1_60*(fn(i+3)-fn(i-3))+0.15d0*(fn(i-2)-fn(i+2))+0.75d0*(fn(i+1)-fn(i-1))
         end do
         i = ((n2-3) - (n1+6) + 6)/6
         i = (n1+6) + 6*i
         if(n2-2 >= ista_nrc .and. n2-2 <= iend_nrc) then
            if(i.eq.n2-2) then
               df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*fn(n2-2)-d4_3*fn(n2-3)+0.5d0*fn(n2-4)-d2_15*fn(n2-5)+d1_60*fn(n2-6)
            end if
         end if
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=d1_6*fn(n2)+d77_60*fn(n2-1)-2.5d0*fn(n2-2)+d5_3*fn(n2-3)-d5_6*fn(n2-4)+0.25d0*fn(n2-5)-d1_30*fn(n2-6)
            end if
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.eq.n2) then
               df(n2)=d147_60*fn(n2)-6.d0*fn(n2-1)+7.5d0*fn(n2-2)-d20_3*fn(n2-3)+3.75d0*fn(n2-4)-1.2d0*fn(n2-5)+d1_6*fn(n2-6)
            end if
         end if
      case default
         do i=max(n1+dn,ist),min(n2-3,ien),dn
            df(i)=d1_60*(fn(i+3)-fn(i-3)) &
                            +0.15d0*(fn(i-2)-fn(i+2)) &
                            +0.75d0*(fn(i+1)-fn(i-1))
         end do
         i = ((n2-3) - (n1+dn) + dn)/dn
         i = (n1+dn) + dn*i
         if(n2-2 >= ista_nrc .and. n2-2 <= iend_nrc) then
            if(i.eq.n2-2) then
               df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*fn(n2-2)-d4_3*fn(n2-3)&
               +0.5d0*fn(n2-4)-d2_15*fn(n2-5)+d1_60*fn(n2-6)
            end if
         end if
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=d1_6*fn(n2)+d77_60*fn(n2-1)-2.5d0*fn(n2-2)+d5_3*fn(n2-3)&
               -d5_6*fn(n2-4)+0.25d0*fn(n2-5)-d1_30*fn(n2-6)
            end if
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.eq.n2) then
               df(n2)=d147_60*fn(n2)-6.d0*fn(n2-1)+7.5d0*fn(n2-2)-d20_3*fn(n2-3)&
               +3.75d0*fn(n2-4)-1.2d0*fn(n2-5)+d1_6*fn(n2-6)
            end if
         end if
   end select
   return
   end subroutine diff_7points2_3D
! ==============================================================================
   
!=====================================================================
   subroutine diff_7points3(ier,n1,n2,dn,fn,df)                       
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: fn(*)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(*)
   
   integer:: i
   real(8):: f0,f1,f2,f3,f4,f5,f6
   
   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d5_3=1.666666666666d0
   real(8),parameter:: d20_3=6.666666666666d0
   real(8),parameter:: d1_6=0.1666666666666d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d7_12=0.58333333333333d0
   real(8),parameter:: d2_15=0.1333333333333d0
   real(8),parameter:: d1_30=0.0333333333333d0
   real(8),parameter:: d1_60=0.01666666666666d0
   real(8),parameter:: d77_60=1.28333333333333d0
   real(8),parameter:: d147_60=2.45d0
   
   ier=0
   f0=fn(n1)
   f1=fn(n1+dn)
   f2=fn(n1+2*dn)
   f3=fn(n1+3*dn)
   f4=fn(n1+4*dn)
   f5=fn(n1+5*dn)
   f6=fn(n1+6*dn)
   df(n1)=-d1_6*f6+1.2d0*f5-3.75d0*f4+d20_3*f3-7.5d0*f2+6.d0*f1-d147_60*f0
   df(n1+dn)=d1_30*f6-0.25d0*f5+d5_6*f4-d5_3*f3+2.5d0*f2-d77_60*f1-d1_6*f0
   df(n1+2*dn)=-d1_60*f6+d2_15*f5-0.5d0*f4+d4_3*f3-d7_12*f2-0.4d0*f1+d1_30*f0
   do i=n1+3*dn,n2-4*dn,dn
        df(i)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
        f0=f1
        f1=f2
        f2=f3
        f3=f4
        f4=f5
        f5=f6
        f6=fn(i+4*dn)
   end do
   df(n2-3*dn)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
   df(n2-2*dn)=-d1_30*f6+0.4d0*f5+d7_12*f4-d4_3*f3+0.5d0*f2-d2_15*f1+d1_60*f0
   df(n2-dn)=d1_6*f6+d77_60*f5-2.5d0*f4+d5_3*f3-d5_6*f2+0.25d0*f1-d1_30*f0
   df(n2)=d147_60*f6-6.d0*f5+7.5d0*f4-d20_3*f3+3.75d0*f2-1.2d0*f1+d1_6*f0
   return
   end subroutine diff_7points3
! === For nrc decomposion. by takto 2012/12/05 =================================
   subroutine diff_7points3_3D(ier,n1,n2,dn,fn,df,ista_nrc,iend_nrc,ist,ien)
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
!  real(8),intent(in)  :: fn(ista_nrc-3*dn:iend_nrc+3*dn)
   real(8),intent(in)  :: fn(ista_nrc-6*dn:iend_nrc+6*dn)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(ista_nrc:iend_nrc)
   integer,intent(in)  :: ista_nrc, iend_nrc, ist, ien

   integer:: i

   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d5_3=1.666666666666d0
   real(8),parameter:: d20_3=6.666666666666d0
   real(8),parameter:: d1_6=0.1666666666666d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d7_12=0.58333333333333d0
   real(8),parameter:: d2_15=0.1333333333333d0
   real(8),parameter:: d1_30=0.0333333333333d0
   real(8),parameter:: d1_60=0.01666666666666d0
   real(8),parameter:: d77_60=1.28333333333333d0
   real(8),parameter:: d147_60=2.45d0
   
   ier=0
   if(n1 >= ista_nrc .and. n1 <= iend_nrc) then
      df(n1)=-d1_6*fn(n1+6*dn)+1.2d0*fn(n1+5*dn)-3.75d0*fn(n1+4*dn)+d20_3*fn(n1+3*dn) &
           & -7.5d0*fn(n1+2*dn)+6.d0*fn(n1+dn)-d147_60*fn(n1)
   end if
   if(n1+dn >= ista_nrc .and. n1+dn <= iend_nrc) then
      df(n1+dn)=d1_30*fn(n1+6*dn)-0.25d0*fn(n1+5*dn)+d5_6*fn(n1+4*dn)-d5_3*fn(n1+3*dn) &
             & +2.5d0*fn(n1+2*dn)-d77_60*fn(n1+dn)-d1_6*fn(n1)
   end if
   if(n1+2*dn >= ista_nrc .and. n1+2*dn <= iend_nrc) then
      df(n1+2*dn)=-d1_60*fn(n1+6*dn)+d2_15*fn(n1+5*dn)-0.5d0*fn(n1+4*dn)+d4_3*fn(n1+3*dn) &
                & -d7_12*fn(n1+2*dn)-0.4d0*fn(n1+dn)+d1_30*fn(n1)
   end if
   do i=max(n1+3*dn,ist),min(n2-4*dn,ien),dn
      df(i)=d1_60*(fn(i+3*dn)-fn(i-3*dn))+0.15d0*(fn(i-2*dn)-fn(i+2*dn))+0.75d0*(fn(i+dn)-fn(i-dn))
   end do
   i = ((n2-4*dn) - (n1+3*dn) + dn)/dn
   i = (n1+3*dn) + dn*i
   if(n2-3*dn >= ista_nrc .and. n2-3*dn <= iend_nrc) then
      df(n2-3*dn)=d1_60*(fn(i+3*dn)-fn(i-3*dn))+0.15d0*(fn(i-2*dn)-fn(i+2*dn))+0.75d0*(fn(i+dn)-fn(i-dn))
   end if
   if(n2-2*dn >= ista_nrc .and. n2-2*dn <= iend_nrc) then
      df(n2-2*dn)=-d1_30*fn(i+3*dn)+0.4d0*fn(i+2*dn)+d7_12*fn(i+dn)-d4_3*fn(i)+0.5d0*fn(i-dn)-d2_15*fn(i-2*dn)+d1_60*fn(i-3*dn)
   end if
   if(n2-dn >= ista_nrc .and. n2-dn <= iend_nrc) then
      df(n2-dn)=d1_6*fn(i+3*dn)+d77_60*fn(i+2*dn)-2.5d0*fn(i+dn)+d5_3*fn(i)-d5_6*fn(i-dn)+0.25d0*fn(i-2*dn)-d1_30*fn(i-3*dn)
   end if
   if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
      df(n2)=d147_60*fn(i+3*dn)-6.d0*fn(i+2*dn)+7.5d0*fn(i+dn)-d20_3*fn(i)+3.75d0*fn(i-dn)-1.2d0*fn(i-2*dn)+d1_6*fn(i-3*dn)
   end if
   return
   end subroutine diff_7points3_3D
! ==============================================================================
   
!=====================================================================
   subroutine ddiff_7points(ier,n1,n2,fn,df,ddf)                       
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2
   real(8),intent(in)  :: fn(n1:n2)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(n1:n2)
   real(8),intent(out) :: ddf(n1:n2)
   
   integer:: i
   real(8):: f0,f1,f2,f3,f4,f5,f6
   
   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d5_3=1.666666666666d0
   real(8),parameter:: d7_3=2.333333333333d0
   real(8),parameter:: d20_3=6.666666666666d0
   real(8),parameter:: d1_6=0.1666666666666d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d10_9=1.111111111111d0
   real(8),parameter:: d254_9=28.222222222222d0
   real(8),parameter:: d1_12=0.08333333333333d0
   real(8),parameter:: d7_12=0.58333333333333d0
   real(8),parameter:: d17_12=1.41666666666666d0
   real(8),parameter:: d19_12=1.58333333333333d0
   real(8),parameter:: d1_15=0.0666666666666d0
   real(8),parameter:: d2_15=0.1333333333333d0
   real(8),parameter:: d19_15=1.2666666666666d0
   real(8),parameter:: d47_18=2.6111111111111d0
   real(8),parameter:: d49_18=2.7222222222222d0
   real(8),parameter:: d1_30=0.0333333333333d0
   real(8),parameter:: d203_45=4.5111111111111d0
   real(8),parameter:: d1_60=0.01666666666666d0
   real(8),parameter:: d31_60=0.51666666666666d0
   real(8),parameter:: d49_60=0.81666666666666d0
   real(8),parameter:: d77_60=1.28333333333333d0
   real(8),parameter:: d147_60=2.45d0
   real(8),parameter:: d1_90=0.0111111111111d0
   real(8),parameter:: d13_180=0.07222222222222d0
   real(8),parameter:: d137_180=0.76111111111111d0
   
   ier=0
   f0=fn(n1)
   f1=fn(n1+1)
   f2=fn(n1+2)
   f3=fn(n1+3)
   f4=fn(n1+4)
   f5=fn(n1+5)
   f6=fn(n1+6)
   df(n1)=-d1_6*f6+1.2d0*f5-3.75d0*f4+d20_3*f3-7.5d0*f2+6.d0*f1-d147_60*f0
   df(n1+1)=d1_30*f6-0.25d0*f5+d5_6*f4-d5_3*f3+2.5d0*f2-d77_60*f1-d1_6*f0
   df(n1+2)=-d1_60*f6+d2_15*f5-0.5d0*f4+d4_3*f3-d7_12*f2-0.4d0*f1+d1_30*f0
   ddf(n1)=d137_180*f6-5.4d0*f5+16.5d0*f4-d254_9*f3+29.25d0*f2-17.4d0*f1+d203_45*f0
   ddf(n1+1)=-d13_180*f6+d31_60*f5-d19_12*f4+d47_18*f3-d17_12*f2-d49_60*f1+d137_180*f0
   ddf(n1+2)=d1_90*f6-d1_15*f5+d1_12*f4+d10_9*f3-d7_3*f2+d19_15*f1-d13_180*f0
   do i=n1+3,n2-4
        df(i)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
        ddf(i)=d1_90*(f6+f0)-0.15d0*(f5+f1)+1.5d0*(f4+f2)-d49_18*f3
        f0=f1
        f1=f2
        f2=f3
        f3=f4
        f4=f5
        f5=f6
        f6=fn(i+4)
   end do
   df(n2-3)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
   df(n2-2)=-d1_30*f6+0.4d0*f5+d7_12*f4-d4_3*f3+0.5d0*f2-d2_15*f1+d1_60*f0
   df(n2-1)=d1_6*f6+d77_60*f5-2.5d0*f4+d5_3*f3-d5_6*f2+0.25d0*f1-d1_30*f0
   df(n2)=d147_60*f6-6.d0*f5+7.5d0*f4-d20_3*f3+3.75d0*f2-1.2d0*f1+d1_6*f0
   ddf(n2-3)=d1_90*(f6+f0)-0.15d0*(f5+f1)+1.5d0*(f4+f2)-d49_18*f3
   ddf(n2-2)=-d13_180*f6+d19_15*f5-d7_3*f4+d10_9*f3+d1_12*f2-d1_15*f1+d1_90*f0
   ddf(n2-1)=d137_180*f6-d49_60*f5-d17_12*f4+d47_18*f3-d19_12*f2+d31_60*f1-d13_180*f0
   ddf(n2)=d203_45*f6-17.4d0*f5+29.25d0*f4-d254_9*f3+16.5d0*f2-5.4d0*f1+d137_180*f0
   return
   end subroutine ddiff_7points
   
!=====================================================================
   subroutine ddiff_7points2(ier,n1,n2,dn,fn,df,ddf)                       
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: fn(*)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(*)
   real(8),intent(out) :: ddf(*)
   
   integer:: i
   real(8):: f0,f1,f2,f3,f4,f5,f6
   
   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d5_3=1.666666666666d0
   real(8),parameter:: d7_3=2.333333333333d0
   real(8),parameter:: d20_3=6.666666666666d0
   real(8),parameter:: d1_6=0.1666666666666d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d10_9=1.111111111111d0
   real(8),parameter:: d254_9=28.222222222222d0
   real(8),parameter:: d1_12=0.08333333333333d0
   real(8),parameter:: d7_12=0.58333333333333d0
   real(8),parameter:: d17_12=1.41666666666666d0
   real(8),parameter:: d19_12=1.58333333333333d0
   real(8),parameter:: d1_15=0.0666666666666d0
   real(8),parameter:: d2_15=0.1333333333333d0
   real(8),parameter:: d19_15=1.2666666666666d0
   real(8),parameter:: d47_18=2.6111111111111d0
   real(8),parameter:: d49_18=2.7222222222222d0
   real(8),parameter:: d1_30=0.0333333333333d0
   real(8),parameter:: d203_45=4.5111111111111d0
   real(8),parameter:: d1_60=0.01666666666666d0
   real(8),parameter:: d31_60=0.51666666666666d0
   real(8),parameter:: d49_60=0.81666666666666d0
   real(8),parameter:: d77_60=1.28333333333333d0
   real(8),parameter:: d147_60=2.45d0
   real(8),parameter:: d1_90=0.0111111111111d0
   real(8),parameter:: d13_180=0.07222222222222d0
   real(8),parameter:: d137_180=0.76111111111111d0
   
   ier=0
   f0=fn(n1)
   f1=fn(n1+1)
   f2=fn(n1+2)
   f3=fn(n1+3)
   f4=fn(n1+4)
   f5=fn(n1+5)
   f6=fn(n1+6)
   df(n1)=-d1_6*f6+1.2d0*f5-3.75d0*f4+d20_3*f3-7.5d0*f2+6.d0*f1-d147_60*f0
   ddf(n1)=d137_180*f6-5.4d0*f5+16.5d0*f4-d254_9*f3+29.25d0*f2-17.4d0*f1+d203_45*f0
   
    select case (dn)
    case (1)
        df(n1+1)=d1_30*f6-0.25d0*f5+d5_6*f4-d5_3*f3+2.5d0*f2-d77_60*f1-d1_6*f0
        df(n1+2)=-d1_60*f6+d2_15*f5-0.5d0*f4+d4_3*f3-d7_12*f2-0.4d0*f1+d1_30*f0
        ddf(n1+1)=-d13_180*f6+d31_60*f5-d19_12*f4+d47_18*f3-d17_12*f2-d49_60*f1+d137_180*f0
        ddf(n1+2)=d1_90*f6-d1_15*f5+d1_12*f4+d10_9*f3-d7_3*f2+d19_15*f1-d13_180*f0
        do i=n1+3,n2-4
            df(i)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
            ddf(i)=d1_90*(f6+f0)-0.15d0*(f5+f1)+1.5d0*(f4+f2)-d49_18*f3
            f0=f1
            f1=f2
            f2=f3
            f3=f4
            f4=f5
            f5=f6
            f6=fn(i+4)
        end do
        df(n2-3)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
        df(n2-2)=-d1_30*f6+0.4d0*f5+d7_12*f4-d4_3*f3+0.5d0*f2-d2_15*f1+d1_60*f0
        df(n2-1)=d1_6*f6+d77_60*f5-2.5d0*f4+d5_3*f3-d5_6*f2+0.25d0*f1-d1_30*f0
        df(n2)=d147_60*f6-6.d0*f5+7.5d0*f4-d20_3*f3+3.75d0*f2-1.2d0*f1+d1_6*f0
        ddf(n2-3)=d1_90*(f6+f0)-0.15d0*(f5+f1)+1.5d0*(f4+f2)-d49_18*f3
        ddf(n2-2)=-d13_180*f6+d19_15*f5-d7_3*f4+d10_9*f3+d1_12*f2-d1_15*f1+d1_90*f0
        ddf(n2-1)=d137_180*f6-d49_60*f5-d17_12*f4+d47_18*f3-d19_12*f2+d31_60*f1-d13_180*f0
        ddf(n2)=d203_45*f6-17.4d0*f5+29.25d0*f4-d254_9*f3+16.5d0*f2-5.4d0*f1+d137_180*f0
    case (2)
        df(n1+2)=-d1_60*f6+d2_15*f5-0.5d0*f4+d4_3*f3-d7_12*f2-0.4d0*f1+d1_30*f0
        ddf(n1+2)=d1_90*f6-d1_15*f5+d1_12*f4+d10_9*f3-d7_3*f2+d19_15*f1-d13_180*f0
        f6=f5
        f5=f4
        f4=f3
        f3=f2
        f2=f1
        do i=n1+4,n2-3,2
            f0=f2
            f1=f3
            f2=f4
            f3=f5
            f4=f6
            f5=fn(i+2)
            f6=fn(i+3)
            df(i)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
            ddf(i)=d1_90*(f6+f0)-0.15d0*(f5+f1)+1.5d0*(f4+f2)-d49_18*f3
        end do
        if(i.eq.n2-2) then
            df(n2-2)=-d1_30*fn(n2)+0.4d0*f6+d7_12*f5-d4_3*f4+0.5d0*f3-d2_15*f2+d1_60*f1
            ddf(n2-2)=-d13_180*fn(n2)+d19_15*f6-d7_3*f5+d10_9*f4+d1_12*f3-d1_15*f2+d1_90*f1
            df(n2)=d147_60*fn(n2)-6.d0*f6+7.5d0*f5-d20_3*f4+3.75d0*f3-1.2d0*f2+d1_6*f1
            ddf(n2)=d203_45*fn(n2)-17.4d0*f6+29.25d0*f5-d254_9*f4+16.5d0*f3-5.4d0*f2+d137_180*f1
        else if(i.eq.n2-1) then
            df(n2-1)=d1_6*f6+d77_60*f5-2.5d0*f4+d5_3*f3-d5_6*f2+0.25d0*f1-d1_30*f0
            ddf(n2-1)=d137_180*f6-d49_60*f5-d17_12*f4+d47_18*f3-d19_12*f2+d31_60*f1-d13_180*f0
        end if
    case (3)
        f6=f3
        f5=f2
        f4=f1
        f3=f0
        do i=n1+3,n2-3,3
            f0=f3
            f1=f4
            f2=f5
            f3=f6
            f4=fn(i+1)
            f5=fn(i+2)
            f6=fn(i+3)
            df(i)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
            ddf(i)=d1_90*(f6+f0)-0.15d0*(f5+f1)+1.5d0*(f4+f2)-d49_18*f3
        end do
        if(i.eq.n2-2) then
            df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*f6-d4_3*f5+0.5d0*f4-d2_15*f3+d1_60*f2
            ddf(n2-2)=-d13_180*fn(n2)+d19_15*fn(n2-1)-d7_3*f6+d10_9*f5+d1_12*f4-d1_15*f3+d1_90*f2
        else if(i.eq.n2-1) then
            df(n2-1)=d1_6*fn(n2)+d77_60*f6-2.5d0*f5+d5_3*f4-d5_6*f3+0.25d0*f2-d1_30*f1
            ddf(n2-1)=d137_180*fn(n2)-d49_60*f6-d17_12*f5+d47_18*f4-d19_12*f3+d31_60*f2-d13_180*f1
        else if(i.eq.n2) then
            df(n2)=d147_60*f6-6.d0*f5+7.5d0*f4-d20_3*f3+3.75d0*f2-1.2d0*f1+d1_6*f0
            ddf(n2)=d203_45*f6-17.4d0*f5+29.25d0*f4-d254_9*f3+16.5d0*f2-5.4d0*f1+d137_180*f0
        end if
    case (4)
        f6=f3
        f5=f2
        f4=f1
        do i=n1+4,n2-3,4
            f0=f4
            f1=f5
            f2=f6
            f3=fn(i)
            f4=fn(i+1)
            f5=fn(i+2)
            f6=fn(i+3)
            df(i)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
            ddf(i)=d1_90*(f6+f0)-0.15d0*(f5+f1)+1.5d0*(f4+f2)-d49_18*f3
        end do
        if(i.eq.n2-2) then
            df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*fn(n2-2)-d4_3*f6+0.5d0*f5-d2_15*f4+d1_60*f3
            ddf(n2-2)=-d13_180*fn(n2)+d19_15*fn(n2-1)-d7_3*fn(n2-2)+d10_9*f6+d1_12*f5-d1_15*f4+d1_90*f3
        else if(i.eq.n2-1) then
            df(n2-1)=d1_6*fn(n2)+d77_60*fn(n2-1)-2.5d0*f6+d5_3*f5-d5_6*f4+0.25d0*f3-d1_30*f2
            ddf(n2-1)=d137_180*fn(n2)-d49_60*fn(n2-1)-d17_12*f6+d47_18*f5-d19_12*f4+d31_60*f3-d13_180*f2
        else if(i.eq.n2) then
            df(n2)=d147_60*fn(n2)-6.d0*f6+7.5d0*f5-d20_3*f4+3.75d0*f3-1.2d0*f2+d1_6*f1
            ddf(n2)=d203_45*fn(n2)-17.4d0*f6+29.25d0*f5-d254_9*f4+16.5d0*f3-5.4d0*f2+d137_180*f1
        end if
    case (5)
        f6=f3
        f5=f2
        do i=n1+5,n2-3,5
            f0=f5
            f1=f6
            f2=fn(i-1)
            f3=fn(i)
            f4=fn(i+1)
            f5=fn(i+2)
            f6=fn(i+3)
            df(i)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
            ddf(i)=d1_90*(f6+f0)-0.15d0*(f5+f1)+1.5d0*(f4+f2)-d49_18*f3
        end do
        if(i.eq.n2-2) then
            df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*fn(n2-2)-d4_3*fn(n2-3)+0.5d0*f6-d2_15*f5+d1_60*f4
            ddf(n2-2)=-d13_180*fn(n2)+d19_15*fn(n2-1)-d7_3*fn(n2-2)+d10_9*fn(n2-3)+d1_12*f6-d1_15*f5+d1_90*f4
        else if(i.eq.n2-1) then
            df(n2-1)=d1_6*fn(n2)+d77_60*fn(n2-1)-2.5d0*fn(n2-2)+d5_3*f6-d5_6*f5+0.25d0*f4-d1_30*f3
            ddf(n2-1)=d137_180*fn(n2)-d49_60*fn(n2-1)-d17_12*fn(n2-2)+d47_18*f6-d19_12*f5+d31_60*f4-d13_180*f3
        else if(i.eq.n2) then
            df(n2)=d147_60*fn(n2)-6.d0*fn(n2-1)+7.5d0*f6-d20_3*f5+3.75d0*f4-1.2d0*f3+d1_6*f2
            ddf(n2)=d203_45*fn(n2)-17.4d0*fn(n2-1)+29.25d0*f6-d254_9*f5+16.5d0*f4-5.4d0*f3+d137_180*f2
        end if
    case (6)
        f6=f3
        do i=n1+6,n2-3,6
            f0=f6
            f1=fn(i-2)
            f2=fn(i-1)
            f3=fn(i)
            f4=fn(i+1)
            f5=fn(i+2)
            f6=fn(i+3)
            df(i)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
            ddf(i)=d1_90*(f6+f0)-0.15d0*(f5+f1)+1.5d0*(f4+f2)-d49_18*f3
        end do
        if(i.eq.n2-2) then
            df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*fn(n2-2)-d4_3*fn(n2-3)+0.5d0*fn(n2-4)-d2_15*f6+d1_60*f5
            ddf(n2-2)=-d13_180*fn(n2)+d19_15*fn(n2-1)-d7_3*fn(n2-2)+d10_9*fn(n2-3)+d1_12*fn(n2-4)-d1_15*f6+d1_90*f5
        else if(i.eq.n2-1) then
            df(n2-1)=d1_6*fn(n2)+d77_60*fn(n2-1)-2.5d0*fn(n2-2)+d5_3*fn(n2-3)-d5_6*f6+0.25d0*f5-d1_30*f4
            ddf(n2-1)=d137_180*fn(n2)-d49_60*fn(n2-1)-d17_12*fn(n2-2)+d47_18*fn(n2-3)-d19_12*f6+d31_60*f5-d13_180*f4
        else if(i.eq.n2) then
            df(n2)=d147_60*fn(n2)-6.d0*fn(n2-1)+7.5d0*fn(n2-2)-d20_3*f6+3.75d0*f5-1.2d0*f4+d1_6*f3
            ddf(n2)=d203_45*fn(n2)-17.4d0*fn(n2-1)+29.25d0*fn(n2-2)-d254_9*f6+16.5d0*f5-5.4d0*f4+d137_180*f3
        end if
    case default
        do i=n1+dn,n2-3,dn
            f0=fn(i-3)
            f1=fn(i-2)
            f2=fn(i-1)
            f3=fn(i)
            f4=fn(i+1)
            f5=fn(i+2)
            f6=fn(i+3)
            df(i)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
            ddf(i)=d1_90*(f6+f0)-0.15d0*(f5+f1)+1.5d0*(f4+f2)-d49_18*f3
        end do
        if(i.eq.n2-2) then
            df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*fn(n2-2)-d4_3*fn(n2-3)&
            +0.5d0*fn(n2-4)-d2_15*fn(n2-5)+d1_60*fn(n2-6)
            ddf(n2-2)=-d13_180*fn(n2)+d19_15*fn(n2-1)-d7_3*fn(n2-2)+d10_9*fn(n2-3)&
            +d1_12*fn(n2-4)-d1_15*fn(n2-5)+d1_90*fn(n2-6)
        else if(i.eq.n2-1) then
            df(n2-1)=d1_6*fn(n2)+d77_60*fn(n2-1)-2.5d0*fn(n2-2)+d5_3*fn(n2-3)&
            -d5_6*fn(n2-4)+0.25d0*fn(n2-5)-d1_30*fn(n2-6)
            ddf(n2-1)=d137_180*fn(n2)-d49_60*fn(n2-1)-d17_12*fn(n2-2)&
            +d47_18*fn(n2-3)-d19_12*fn(n2-4)+d31_60*fn(n2-5)-d13_180*fn(n2-6)
        else if(i.eq.n2) then
            df(n2)=d147_60*fn(n2)-6.d0*fn(n2-1)+7.5d0*fn(n2-2)-d20_3*fn(n2-3)&
            +3.75d0*fn(n2-4)-1.2d0*fn(n2-5)+d1_6*fn(n2-6)
            ddf(n2)=d203_45*fn(n2)-17.4d0*fn(n2-1)+29.25d0*fn(n2-2)&
            -d254_9*fn(n2-3)+16.5d0*fn(n2-4)-5.4d0*fn(n2-5)+d137_180*fn(n2-6)
        end if
    end select
   return
   end subroutine ddiff_7points2
! === For nrc decomposion. by takto 2012/12/05 =================================
   subroutine ddiff_7points2_3D(ier,n1,n2,dn,fn,df,ddf,ista_nrc,iend_nrc,ist,ien)
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: fn(ista_nrc-6*dn:iend_nrc+6*dn)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(ista_nrc:iend_nrc)
   real(8),intent(out) :: ddf(ista_nrc:iend_nrc)
   integer,intent(in)  :: ista_nrc, iend_nrc, ist, ien

   integer:: i

   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d5_3=1.666666666666d0
   real(8),parameter:: d7_3=2.333333333333d0
   real(8),parameter:: d20_3=6.666666666666d0
   real(8),parameter:: d1_6=0.1666666666666d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d10_9=1.111111111111d0
   real(8),parameter:: d254_9=28.222222222222d0
   real(8),parameter:: d1_12=0.08333333333333d0
   real(8),parameter:: d7_12=0.58333333333333d0
   real(8),parameter:: d17_12=1.41666666666666d0
   real(8),parameter:: d19_12=1.58333333333333d0
   real(8),parameter:: d1_15=0.0666666666666d0
   real(8),parameter:: d2_15=0.1333333333333d0
   real(8),parameter:: d19_15=1.2666666666666d0
   real(8),parameter:: d47_18=2.6111111111111d0
   real(8),parameter:: d49_18=2.7222222222222d0
   real(8),parameter:: d1_30=0.0333333333333d0
   real(8),parameter:: d203_45=4.5111111111111d0
   real(8),parameter:: d1_60=0.01666666666666d0
   real(8),parameter:: d31_60=0.51666666666666d0
   real(8),parameter:: d49_60=0.81666666666666d0
   real(8),parameter:: d77_60=1.28333333333333d0
   real(8),parameter:: d147_60=2.45d0
   real(8),parameter:: d1_90=0.0111111111111d0
   real(8),parameter:: d13_180=0.07222222222222d0
   real(8),parameter:: d137_180=0.76111111111111d0

   ier=0

   if(n1 >= ista_nrc .and. n1 <= iend_nrc) then
      df(n1)=-d1_6*fn(n1+6)+1.2d0*fn(n1+5)-3.75d0*fn(n1+4)+d20_3*fn(n1+3)-7.5d0*fn(n1+2)+6.d0*fn(n1+1)-d147_60*fn(n1)
      ddf(n1)=d137_180*fn(n1+6)-5.4d0*fn(n1+5)+16.5d0*fn(n1+4)-d254_9*fn(n1+3)+29.25d0*fn(n1+2)-17.4d0*fn(n1+1)+d203_45*fn(n1)
   end if

   select case (dn)
      case (1)
         if(n1+1 >= ista_nrc .and. n1+1 <= iend_nrc) then
            df(n1+1)=d1_30*fn(n1+6)-0.25d0*fn(n1+5)+d5_6*fn(n1+4) &
                  & -d5_3*fn(n1+3)+2.5d0*fn(n1+2)-d77_60*fn(n1+1)-d1_6*fn(n1)
            ddf(n1+1)=-d13_180*fn(n1+6)+d31_60*fn(n1+5)-d19_12*fn(n1+4) &
                    & +d47_18*fn(n1+3)-d17_12*fn(n1+2)-d49_60*fn(n1+1)+d137_180*fn(n1)
         end if
         if(n1+2 >= ista_nrc .and. n1+2 <= iend_nrc) then
            df(n1+2)=-d1_60*fn(n1+6)+d2_15*fn(n1+5)-0.5d0*fn(n1+4) &
                   & +d4_3*fn(n1+3)-d7_12*fn(n1+2)-0.4d0*fn(n1+1)+d1_30*fn(n1)
            ddf(n1+2)=d1_90*fn(n1+6)-d1_15*fn(n1+5)+d1_12*fn(n1+4) &
                   & +d10_9*fn(n1+3)-d7_3*fn(n1+2)+d19_15*fn(n1+1)-d13_180*fn(n1)
         end if
         do i=max(n1+3,ist),min(n2-4,ien)
            df(i)=d1_60*(fn(i+3)-fn(i-3))+0.15d0*(fn(i-2)-fn(i+2))+0.75d0*(fn(i+1)-fn(i-1))
            ddf(i)=d1_90*(fn(i+3)+fn(i-3))-0.15d0*(fn(i+2)+fn(i-2))+1.5d0*(fn(i+1)+fn(i-1))-d49_18*fn(i)
         end do
         i = (n2-4) - (n1+3) + 1
         i = (n1+3) + i
         if(n2-3 >= ista_nrc .and. n2-3 <= iend_nrc) then
            df(n2-3)=d1_60*(fn(i+3)-fn(i-3))+0.15d0*(fn(i-2)-fn(i+2))+0.75d0*(fn(i+1)-fn(i-1))
            ddf(n2-3)=d1_90*(fn(i+3)+fn(i-3))-0.15d0*(fn(i+2)+fn(i-2))+1.5d0*(fn(i+1)+fn(i-1))-d49_18*fn(i)
         end if
         if(n2-2 >= ista_nrc .and. n2-2 <= iend_nrc) then
            df(n2-2)=-d1_30*fn(i+3)+0.4d0*fn(i+2)+d7_12*fn(i+1)-d4_3*fn(i)+0.5d0*fn(i-1)-d2_15*fn(i-2)+d1_60*fn(i-3)
            ddf(n2-2)=-d13_180*fn(i+3)+d19_15*fn(i+2)-d7_3*fn(i+1)+d10_9*fn(i)+d1_12*fn(i-1)-d1_15*fn(i-2)+d1_90*fn(i-3)
         end if
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            df(n2-1)=d1_6*fn(i+3)+d77_60*fn(i+2)-2.5d0*fn(i+1)+d5_3*fn(i)-d5_6*fn(i-1)+0.25d0*fn(i-2)-d1_30*fn(i-3)
            ddf(n2-1)=d137_180*fn(i+3)-d49_60*fn(i+2)-d17_12*fn(i+1)+d47_18*fn(i)-d19_12*fn(i-1)+d31_60*fn(i-2)-d13_180*fn(i-3)
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            df(n2)=d147_60*fn(i+3)-6.d0*fn(i+2)+7.5d0*fn(i+1)-d20_3*fn(i)+3.75d0*fn(i-1)-1.2d0*fn(i-2)+d1_6*fn(i-3)
            ddf(n2)=d203_45*fn(i+3)-17.4d0*fn(i+2)+29.25d0*fn(i+1)-d254_9*fn(i)+16.5d0*fn(i-1)-5.4d0*fn(i-2)+d137_180*fn(i-3)
         end if
      case (2)
         if(n1+2 >= ista_nrc .and. n1+2 <= iend_nrc) then
            df(n1+2)=-d1_60*fn(n1+6)+d2_15*fn(n1+5)-0.5d0*fn(n1+4)+d4_3*fn(n1+3)-d7_12*fn(n1+2)-0.4d0*fn(n1+1)+d1_30*fn(n1)
            ddf(n1+2)=d1_90*fn(n1+6)-d1_15*fn(n1+5)+d1_12*fn(n1+4)+d10_9*fn(n1+3)-d7_3*fn(n1+2)+d19_15*fn(n1+1)-d13_180*fn(n1)
         end if
         do i=max(n1+4,ist),min(n2-3,ien),2
            df(i)=d1_60*(fn(i+3)-fn(i-3))+0.15d0*(fn(i-2)-fn(i+2))+0.75d0*(fn(i+1)-fn(i-1))
            ddf(i)=d1_90*(fn(i+3)+fn(i-3))-0.15d0*(fn(i+2)+fn(i-2))+1.5d0*(fn(i+1)+fn(i-1))-d49_18*fn(i)
         end do
         i = ((n2-3) - (n1+4) + 2)/2
         i = (n1+4) + 2*i
         if(n2-2 >= ista_nrc .and. n2-2 <= iend_nrc) then
            if(i.eq.n2-2) then
               df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(i+1)+d7_12*fn(i)-d4_3*fn(i-1)+0.5d0*fn(i-2)-d2_15*fn(i-3)+d1_60*fn(i-4)
               ddf(n2-2)=-d13_180*fn(n2)+d19_15*fn(i+1)-d7_3*fn(i)+d10_9*fn(i-1)+d1_12*fn(i-2)-d1_15*fn(i-3)+d1_90*fn(i-4)
               df(n2)=d147_60*fn(n2)-6.d0*fn(i+1)+7.5d0*fn(i)-d20_3*fn(i-1)+3.75d0*fn(i-2)-1.2d0*fn(i-3)+d1_6*fn(i-4)
               ddf(n2)=d203_45*fn(n2)-17.4d0*fn(i+1)+29.25d0*fn(i)-d254_9*fn(i-1)+16.5d0*fn(i-2)-5.4d0*fn(i-3)+d137_180*fn(i-4)
            end if
         end if
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=d1_6*fn(i+1)+d77_60*fn(i)-2.5d0*fn(i-1)+d5_3*fn(i-2)-d5_6*fn(i-3)+0.25d0*fn(i-4)-d1_30*fn(i-5)
               ddf(n2-1)=d137_180*fn(i+1)-d49_60*fn(i)-d17_12*fn(i-1)+d47_18*fn(i-2)-d19_12*fn(i-3)+d31_60*fn(i-4)-d13_180*fn(i-5)
            end if
         end if
      case (3)
         do i=max(n1+3,ist),min(n2-3,ien),3
            df(i)=d1_60*(fn(i+3)-fn(i-3))+0.15d0*(fn(i-2)-fn(i+2))+0.75d0*(fn(i+1)-fn(i-1))
            ddf(i)=d1_90*(fn(i+3)+fn(i-3))-0.15d0*(fn(i+2)+fn(i-2))+1.5d0*(fn(i+1)+fn(i-1))-d49_18*fn(i)
         end do
         i = ((n2-3) - (n1+3) + 3)/3
         i = (n1+3) + 3*i
         if(n2-2 >= ista_nrc .and. n2-2 <= iend_nrc) then
            if(i.eq.n2-2) then
               df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*fn(i)-d4_3*fn(i-1)+0.5d0*fn(i-2)-d2_15*fn(i-3)+d1_60*fn(i-4)
               ddf(n2-2)=-d13_180*fn(n2)+d19_15*fn(n2-1)-d7_3*fn(i)+d10_9*fn(i-1)+d1_12*fn(i-2)-d1_15*fn(i-3)+d1_90*fn(i-4)
            end if
         end if
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=d1_6*fn(n2)+d77_60*fn(i)-2.5d0*fn(i-1)+d5_3*fn(i-2)-d5_6*fn(i-3)+0.25d0*fn(i-4)-d1_30*fn(i-5)
               ddf(n2-1)=d137_180*fn(n2)-d49_60*fn(i)-d17_12*fn(i-1)&
                & +d47_18*fn(i-2)-d19_12*fn(i-3)+d31_60*fn(i-4)-d13_180*fn(i-5)
            end if
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.eq.n2) then
               df(n2)=d147_60*fn(i)-6.d0*fn(i-1)+7.5d0*fn(i-2)-d20_3*fn(i-3)+3.75d0*fn(i-4)-1.2d0*fn(i-5)+d1_6*fn(i-6)
               ddf(n2)=d203_45*fn(i)-17.4d0*fn(i-1)+29.25d0*fn(i-2)-d254_9*fn(i-3)+16.5d0*fn(i-4)-5.4d0*fn(i-5)+d137_180*fn(i-6)
            end if
         end if
      case (4)
         do i=max(n1+4,ist),min(n2-3,ien),4
            df(i)=d1_60*(fn(i+3)-fn(i-3))+0.15d0*(fn(i-2)-fn(i+2))+0.75d0*(fn(i+1)-fn(i-1))
            ddf(i)=d1_90*(fn(i+3)+fn(i-3))-0.15d0*(fn(i+2)+fn(i-2))+1.5d0*(fn(i+1)+fn(i-1))-d49_18*fn(i)
         end do
         i = ((n2-3) - (n1+4) + 4)/4
         i = (n1+4) + 4*i
         if(n2-2 >= ista_nrc .and. n2-2 <= iend_nrc) then
            if(i.eq.n2-2) then
               df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*fn(n2-2)-d4_3*fn(i-1)+0.5d0*fn(i-2)-d2_15*fn(i-3)+d1_60*fn(i-4)
               ddf(n2-2)=-d13_180*fn(n2)+d19_15*fn(n2-1)-d7_3*fn(n2-2)+d10_9*fn(i-1)+d1_12*fn(i-2)-d1_15*fn(i-3)+d1_90*fn(i-4)
            end if
         end if
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=d1_6*fn(n2)+d77_60*fn(n2-1)-2.5d0*fn(i-1)+d5_3*fn(i-2)-d5_6*fn(i-3)+0.25d0*fn(i-4)-d1_30*fn(i-5)
               ddf(n2-1)=d137_180*fn(n2)-d49_60*fn(n2-1)-d17_12*fn(i-1)&
                & +d47_18*fn(i-2)-d19_12*fn(i-3)+d31_60*fn(i-4)-d13_180*fn(i-5)
            end if
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.eq.n2) then
               df(n2)=d147_60*fn(n2)-6.d0*fn(i-1)+7.5d0*fn(i-2)-d20_3*fn(i-3)+3.75d0*fn(i-4)-1.2d0*fn(i-5)+d1_6*fn(i-6)
               ddf(n2)=d203_45*fn(n2)-17.4d0*fn(i-1)+29.25d0*fn(i-2)-d254_9*fn(i-3)+16.5d0*fn(i-4)-5.4d0*fn(i-5)+d137_180*fn(i-6)
            end if
        end if
      case (5)
         do i=max(n1+5,ist),min(n2-3,ien),5
            df(i)=d1_60*(fn(i+3)-fn(i-3))+0.15d0*(fn(i-2)-fn(i+2))+0.75d0*(fn(i+1)-fn(i-1))
            ddf(i)=d1_90*(fn(i+3)+fn(i-3))-0.15d0*(fn(i+2)+fn(i-2))+1.5d0*(fn(i+1)+fn(i-1))-d49_18*fn(i)
         end do
         i = ((n2-3) - (n1+5) + 5)/5
         i = (n1+5) + 5*i
         if(n2-2 >= ista_nrc .and. n2-2 <= iend_nrc) then
            if(i.eq.n2-2) then
               df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*fn(n2-2)-d4_3*fn(n2-3)+0.5d0*fn(i-2)-d2_15*fn(i-3)+d1_60*fn(i-4)
               ddf(n2-2)=-d13_180*fn(n2)+d19_15*fn(n2-1)-d7_3*fn(n2-2)+d10_9*fn(n2-3)+d1_12*fn(i-2)-d1_15*fn(i-3)+d1_90*fn(i-4)
            end if
         end if
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=d1_6*fn(n2)+d77_60*fn(n2-1)-2.5d0*fn(n2-2)+d5_3*fn(i-2)-d5_6*fn(i-3)+0.25d0*fn(i-4)-d1_30*fn(i-5)
               ddf(n2-1)=d137_180*fn(n2)-d49_60*fn(n2-1)-d17_12*fn(n2-2)&
                & +d47_18*fn(i-2)-d19_12*fn(i-3)+d31_60*fn(i-4)-d13_180*fn(i-5)
            end if
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.eq.n2) then
               df(n2)=d147_60*fn(n2)-6.d0*fn(n2-1)+7.5d0*fn(i-2)-d20_3*fn(i-3)+3.75d0*fn(i-4)-1.2d0*fn(i-5)+d1_6*fn(i-6)
               ddf(n2)=d203_45*fn(n2)-17.4d0*fn(n2-1)+29.25d0*fn(i-2)&
                & -d254_9*fn(i-3)+16.5d0*fn(i-4)-5.4d0*fn(i-5)+d137_180*fn(i-6)
            end if
         end if
      case (6)
         do i=max(n1+6,ist),min(n2-3,ien),6
            df(i)=d1_60*(fn(i+3)-fn(i-3))+0.15d0*(fn(i-2)-fn(i+2))+0.75d0*(fn(i+1)-fn(i-1))
            ddf(i)=d1_90*(fn(i+3)+fn(i-3))-0.15d0*(fn(i+2)+fn(i-2))+1.5d0*(fn(i+1)+fn(i-1))-d49_18*fn(i)
         end do
         i = ((n2-3) - (n1+6) + 6)/6
         i = (n1+6) + 6*i
         if(n2-2 >= ista_nrc .and. n2-2 <= iend_nrc) then
            if(i.eq.n2-2) then
               df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*fn(n2-2)-d4_3*fn(n2-3)+0.5d0*fn(n2-4)-d2_15*fn(i-3)+d1_60*fn(i-4)
               ddf(n2-2)=-d13_180*fn(n2)+d19_15*fn(n2-1)-d7_3*fn(n2-2)+d10_9*fn(n2-3)+d1_12*fn(n2-4)-d1_15*fn(i-3)+d1_90*fn(i-4)
            end if
         end if
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=d1_6*fn(n2)+d77_60*fn(n2-1)-2.5d0*fn(n2-2)+d5_3*fn(n2-3)-d5_6*fn(i-3)+0.25d0*fn(i-4)-d1_30*fn(i-5)
               ddf(n2-1)=d137_180*fn(n2)-d49_60*fn(n2-1)-d17_12*fn(n2-2)&
                & +d47_18*fn(n2-3)-d19_12*fn(i-3)+d31_60*fn(i-4)-d13_180*fn(i-5)
            end if
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.eq.n2) then
               df(n2)=d147_60*fn(n2)-6.d0*fn(n2-1)+7.5d0*fn(n2-2)-d20_3*fn(i-3)+3.75d0*fn(i-4)-1.2d0*fn(i-5)+d1_6*fn(i-6)
               ddf(n2)=d203_45*fn(n2)-17.4d0*fn(n2-1)+29.25d0*fn(n2-2)&
                & -d254_9*fn(i-3)+16.5d0*fn(i-4)-5.4d0*fn(i-5)+d137_180*fn(i-6)
            end if
         end if
      case default
         do i=max(n1+dn,ist),min(n2-3,ien),dn
            df(i)=d1_60*(fn(i+3)-fn(i-3))+0.15d0*(fn(i-2)-fn(i+2))+0.75d0*(fn(i+1)-fn(i-1))
            ddf(i)=d1_90*(fn(i+3)+fn(i-3))-0.15d0*(fn(i+2)+fn(i-2))+1.5d0*(fn(i+1)+fn(i-1))-d49_18*fn(i)
         end do
         i = ((n2-3) - (n1+dn) + dn)/dn
         i = (n1+dn) + dn*i
         if(n2-2 >= ista_nrc .and. n2-2 <= iend_nrc) then
            if(i.eq.n2-2) then
               df(n2-2)=-d1_30*fn(n2)+0.4d0*fn(n2-1)+d7_12*fn(n2-2)-d4_3*fn(n2-3)&
               +0.5d0*fn(n2-4)-d2_15*fn(n2-5)+d1_60*fn(n2-6)
               ddf(n2-2)=-d13_180*fn(n2)+d19_15*fn(n2-1)-d7_3*fn(n2-2)+d10_9*fn(n2-3)&
               +d1_12*fn(n2-4)-d1_15*fn(n2-5)+d1_90*fn(n2-6)
            end if
         end if
         if(n2-1 >= ista_nrc .and. n2-1 <= iend_nrc) then
            if(i.eq.n2-1) then
               df(n2-1)=d1_6*fn(n2)+d77_60*fn(n2-1)-2.5d0*fn(n2-2)+d5_3*fn(n2-3)&
               -d5_6*fn(n2-4)+0.25d0*fn(n2-5)-d1_30*fn(n2-6)
               ddf(n2-1)=d137_180*fn(n2)-d49_60*fn(n2-1)-d17_12*fn(n2-2)&
               +d47_18*fn(n2-3)-d19_12*fn(n2-4)+d31_60*fn(n2-5)-d13_180*fn(n2-6)
            end if
         end if
         if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
            if(i.eq.n2) then
               df(n2)=d147_60*fn(n2)-6.d0*fn(n2-1)+7.5d0*fn(n2-2)-d20_3*fn(n2-3)&
               +3.75d0*fn(n2-4)-1.2d0*fn(n2-5)+d1_6*fn(n2-6)
               ddf(n2)=d203_45*fn(n2)-17.4d0*fn(n2-1)+29.25d0*fn(n2-2)&
               -d254_9*fn(n2-3)+16.5d0*fn(n2-4)-5.4d0*fn(n2-5)+d137_180*fn(n2-6)
            end if
         end if
   end select
   return
   end subroutine ddiff_7points2_3D
! ==============================================================================
   
!=====================================================================
   subroutine ddiff_7points3(ier,n1,n2,dn,fn,df,ddf)                       
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: fn(*)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(*)
   real(8),intent(out) :: ddf(*)
   
   integer:: i
   real(8):: f0,f1,f2,f3,f4,f5,f6
   
   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d5_3=1.666666666666d0
   real(8),parameter:: d7_3=2.333333333333d0
   real(8),parameter:: d20_3=6.666666666666d0
   real(8),parameter:: d1_6=0.1666666666666d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d10_9=1.111111111111d0
   real(8),parameter:: d254_9=28.222222222222d0
   real(8),parameter:: d1_12=0.08333333333333d0
   real(8),parameter:: d7_12=0.58333333333333d0
   real(8),parameter:: d17_12=1.41666666666666d0
   real(8),parameter:: d19_12=1.58333333333333d0
   real(8),parameter:: d1_15=0.0666666666666d0
   real(8),parameter:: d2_15=0.1333333333333d0
   real(8),parameter:: d19_15=1.2666666666666d0
   real(8),parameter:: d47_18=2.6111111111111d0
   real(8),parameter:: d49_18=2.7222222222222d0
   real(8),parameter:: d1_30=0.0333333333333d0
   real(8),parameter:: d203_45=4.5111111111111d0
   real(8),parameter:: d1_60=0.01666666666666d0
   real(8),parameter:: d31_60=0.51666666666666d0
   real(8),parameter:: d49_60=0.81666666666666d0
   real(8),parameter:: d77_60=1.28333333333333d0
   real(8),parameter:: d147_60=2.45d0
   real(8),parameter:: d1_90=0.0111111111111d0
   real(8),parameter:: d13_180=0.07222222222222d0
   real(8),parameter:: d137_180=0.76111111111111d0
   
   ier=0
   f0=fn(n1)
   f1=fn(n1+dn)
   f2=fn(n1+2*dn)
   f3=fn(n1+3*dn)
   f4=fn(n1+4*dn)
   f5=fn(n1+5*dn)
   f6=fn(n1+6*dn)
   df(n1)=-d1_6*f6+1.2d0*f5-3.75d0*f4+d20_3*f3-7.5d0*f2+6.d0*f1-d147_60*f0
   df(n1+dn)=d1_30*f6-0.25d0*f5+d5_6*f4-d5_3*f3+2.5d0*f2-d77_60*f1-d1_6*f0
   df(n1+2*dn)=-d1_60*f6+d2_15*f5-0.5d0*f4+d4_3*f3-d7_12*f2-0.4d0*f1+d1_30*f0
   ddf(n1)=d137_180*f6-5.4d0*f5+16.5d0*f4-d254_9*f3+29.25d0*f2-17.4d0*f1+d203_45*f0
   ddf(n1+dn)=-d13_180*f6+d31_60*f5-d19_12*f4+d47_18*f3-d17_12*f2-d49_60*f1+d137_180*f0
   ddf(n1+2*dn)=d1_90*f6-d1_15*f5+d1_12*f4+d10_9*f3-d7_3*f2+d19_15*f1-d13_180*f0
   do i=n1+3*dn,n2-4*dn,dn
        df(i)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
        ddf(i)=d1_90*(f6+f0)-0.15d0*(f5+f1)+1.5d0*(f4+f2)-d49_18*f3
        f0=f1
        f1=f2
        f2=f3
        f3=f4
        f4=f5
        f5=f6
        f6=fn(i+4*dn)
   end do
   df(n2-3*dn)=d1_60*(f6-f0)+0.15d0*(f1-f5)+0.75d0*(f4-f2)
   df(n2-2*dn)=-d1_30*f6+0.4d0*f5+d7_12*f4-d4_3*f3+0.5d0*f2-d2_15*f1+d1_60*f0
   df(n2-dn)=d1_6*f6+d77_60*f5-2.5d0*f4+d5_3*f3-d5_6*f2+0.25d0*f1-d1_30*f0
   df(n2)=d147_60*f6-6.d0*f5+7.5d0*f4-d20_3*f3+3.75d0*f2-1.2d0*f1+d1_6*f0
   ddf(n2-3*dn)=d1_90*(f6+f0)-0.15d0*(f5+f1)+1.5d0*(f4+f2)-d49_18*f3
   ddf(n2-2*dn)=-d13_180*f6+d19_15*f5-d7_3*f4+d10_9*f3+d1_12*f2-d1_15*f1+d1_90*f0
   ddf(n2-dn)=d137_180*f6-d49_60*f5-d17_12*f4+d47_18*f3-d19_12*f2+d31_60*f1-d13_180*f0
   ddf(n2)=d203_45*f6-17.4d0*f5+29.25d0*f4-d254_9*f3+16.5d0*f2-5.4d0*f1+d137_180*f0
   return
   end subroutine ddiff_7points3
! === For nrc decomposion. by takto 2012/12/05 =================================
   subroutine ddiff_7points3_3D(ier,n1,n2,dn,fn,df,ddf,ista_nrc,iend_nrc,ist,ien)
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: fn(ista_nrc-6*dn:iend_nrc+6*dn)
   integer,intent(out) :: ier
   real(8),intent(out) :: df(ista_nrc:iend_nrc)
   real(8),intent(out) :: ddf(ista_nrc:iend_nrc)
   integer,intent(in)  :: ista_nrc, iend_nrc, ist, ien

   integer:: i

   real(8),parameter:: d4_3=1.333333333333d0
   real(8),parameter:: d5_3=1.666666666666d0
   real(8),parameter:: d7_3=2.333333333333d0
   real(8),parameter:: d20_3=6.666666666666d0
   real(8),parameter:: d1_6=0.1666666666666d0
   real(8),parameter:: d5_6=0.8333333333333d0
   real(8),parameter:: d10_9=1.111111111111d0
   real(8),parameter:: d254_9=28.222222222222d0
   real(8),parameter:: d1_12=0.08333333333333d0
   real(8),parameter:: d7_12=0.58333333333333d0
   real(8),parameter:: d17_12=1.41666666666666d0
   real(8),parameter:: d19_12=1.58333333333333d0
   real(8),parameter:: d1_15=0.0666666666666d0
   real(8),parameter:: d2_15=0.1333333333333d0
   real(8),parameter:: d19_15=1.2666666666666d0
   real(8),parameter:: d47_18=2.6111111111111d0
   real(8),parameter:: d49_18=2.7222222222222d0
   real(8),parameter:: d1_30=0.0333333333333d0
   real(8),parameter:: d203_45=4.5111111111111d0
   real(8),parameter:: d1_60=0.01666666666666d0
   real(8),parameter:: d31_60=0.51666666666666d0
   real(8),parameter:: d49_60=0.81666666666666d0
   real(8),parameter:: d77_60=1.28333333333333d0
   real(8),parameter:: d147_60=2.45d0
   real(8),parameter:: d1_90=0.0111111111111d0
   real(8),parameter:: d13_180=0.07222222222222d0
   real(8),parameter:: d137_180=0.76111111111111d0

   ier=0
   if(n1 >= ista_nrc .and. n1 <= iend_nrc) then
      df(n1)=-d1_6*fn(n1+6*dn)+1.2d0*fn(n1+5*dn)-3.75d0*fn(n1+4*dn) &
           & +d20_3*fn(n1+3*dn)-7.5d0*fn(n1+2*dn)+6.d0*fn(n1+dn)-d147_60*fn(n1)
      ddf(n1)=d137_180*fn(n1+6*dn)-5.4d0*fn(n1+5*dn)+16.5d0*fn(n1+4*dn) &
           & -d254_9*fn(n1+3*dn)+29.25d0*fn(n1+2*dn)-17.4d0*fn(n1+dn)+d203_45*fn(n1)
   end if
   if(n1+dn >= ista_nrc .and. n1+dn <= iend_nrc) then
      df(n1+dn)=d1_30*fn(n1+6*dn)-0.25d0*fn(n1+5*dn)+d5_6*fn(n1+4*dn) &
             & -d5_3*fn(n1+3*dn)+2.5d0*fn(n1+2*dn)-d77_60*fn(n1+dn)-d1_6*fn(n1)
      ddf(n1+dn)=-d13_180*fn(n1+6*dn)+d31_60*fn(n1+5*dn)-d19_12*fn(n1+4*dn) &
               & +d47_18*fn(n1+3*dn)-d17_12*fn(n1+2*dn)-d49_60*fn(n1+dn)+d137_180*fn(n1)
   end if
   if(n1+2*dn >= ista_nrc .and. n1+2*dn <= iend_nrc) then
      df(n1+2*dn)=-d1_60*fn(n1+6*dn)+d2_15*fn(n1+5*dn)-0.5d0*fn(n1+4*dn) &
                & +d4_3*fn(n1+3*dn)-d7_12*fn(n1+2*dn)-0.4d0*fn(n1+dn)+d1_30*fn(n1)
      ddf(n1+2*dn)=d1_90*fn(n1+6*dn)-d1_15*fn(n1+5*dn)+d1_12*fn(n1+4*dn) &
                & +d10_9*fn(n1+3*dn)-d7_3*fn(n1+2*dn)+d19_15*fn(n1+dn)-d13_180*fn(n1)
   end if
   do i=max(n1+3*dn,ist),min(n2-4*dn,ien),dn
      df(i)=d1_60*(fn(i+3*dn)-fn(i-3*dn))+0.15d0*(fn(i-2*dn)-fn(i+2*dn))+0.75d0*(fn(i+dn)-fn(i-dn))
      ddf(i)=d1_90*(fn(i+3*dn)+fn(i-3*dn))-0.15d0*(fn(i+2*dn)+fn(i-2*dn))+1.5d0*(fn(i+dn)+fn(i-dn))-d49_18*fn(i)
   end do
   i = ((n2-4*dn) - (n1+3*dn) + dn)/dn
   i = (n1+3*dn) + dn*i
   if(n2-3*dn >= ista_nrc .and. n2-3*dn <= iend_nrc) then
      df(n2-3*dn)=d1_60*(fn(i+3*dn)-fn(i-3*dn))+0.15d0*(fn(i-2*dn)-fn(i+2*dn))+0.75d0*(fn(i+dn)-fn(i-dn))
      ddf(n2-3*dn)=d1_90*(fn(i+3*dn)+fn(i-3*dn))-0.15d0*(fn(i+2*dn)+fn(i-2*dn))+1.5d0*(fn(i+dn)+fn(i-dn))-d49_18*fn(i)
   end if
   if(n2-2*dn >= ista_nrc .and. n2-2*dn <= iend_nrc) then
      df(n2-2*dn)=-d1_30*fn(i+3*dn)+0.4d0*fn(i+2*dn)+d7_12*fn(i+dn)-d4_3*fn(i)+0.5d0*fn(i-dn)-d2_15*fn(i-2*dn)+d1_60*fn(i-3*dn)
      ddf(n2-2*dn)=-d13_180*fn(i+3*dn)+d19_15*fn(i+2*dn)-d7_3*fn(i+dn)&
                & +d10_9*fn(i)+d1_12*fn(i-dn)-d1_15*fn(i-2*dn)+d1_90*fn(i-3*dn)
   end if
   if(n2-dn >= ista_nrc .and. n2-dn <= iend_nrc) then
      df(n2-dn)=d1_6*fn(i+3*dn)+d77_60*fn(i+2*dn)-2.5d0*fn(i+dn)+d5_3*fn(i)-d5_6*fn(i-dn)+0.25d0*fn(i-2*dn)-d1_30*fn(i-3*dn)
      ddf(n2-dn)=d137_180*fn(i+3*dn)-d49_60*fn(i+2*dn)-d17_12*fn(i+dn)&
                & +d47_18*fn(i)-d19_12*fn(i-dn)+d31_60*fn(i-2*dn)-d13_180*fn(i-3*dn)
   end if
   if(n2 >= ista_nrc .and. n2 <= iend_nrc) then
      df(n2)=d147_60*fn(i+3*dn)-6.d0*fn(i+2*dn)+7.5d0*fn(i+dn)-d20_3*fn(i)+3.75d0*fn(i-dn)-1.2d0*fn(i-2*dn)+d1_6*fn(i-3*dn)
      ddf(n2)=d203_45*fn(i+3*dn)-17.4d0*fn(i+2*dn)+29.25d0*fn(i+dn)&
                & -d254_9*fn(i)+16.5d0*fn(i-dn)-5.4d0*fn(i-2*dn)+d137_180*fn(i-3*dn)
   end if
   return
   end subroutine ddiff_7points3_3D
! ==============================================================================
