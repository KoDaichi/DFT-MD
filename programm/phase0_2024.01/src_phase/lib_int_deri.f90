!
!  Software name : CIAO (Code for Investigating Atomic Orbitals)
!  Subroutine(s) : set_weight_unif, set_weight_exp
!                  set_open_weight_exp, diff_exp, diff3_exp, diff4_exp
!                  calc_diff_exp, calc_ddiff_exp
!  Author(s)     : Masakuni Okamoto (August 25, 2003)
!
!  The license of the code and contact address :
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)
!
!=====================================================================
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
   subroutine set_weight_exp2(ier,n1,n2,dn,r,w)
!=====================================================================
!
!  M. Okamoto
!
!---------------------------------------------------------------------
   implicit none
   integer,intent(in)  :: n1, n2, dn
   real(8),intent(in)  :: r(*)
   integer,intent(out) :: ier
   real(8),intent(out) :: w(*)
   real(8) :: h
   integer :: n, i
   ier = 0
   n = n2-n1+1
   if (n < 6*dn) then
      ier = 1 ; go to 99
   end if
   h = log(r(n2)/r(n1)) / dble(n2-n1) * dble(dn)
   w(n1  )      =  3.d0/8.d0  * r(n1        ) * h
   w(n1+1*dn)   =  7.d0/6.d0  * r(n1+1*dn   ) * h
   w(n1+2*dn)   = 23.d0/24.d0 * r(n1+2*dn   ) * h
   w(n2-2*dn)   = 23.d0/24.d0 * r(n2-2*dn   ) * h
   w(n2-1*dn)   =  7.d0/6.d0  * r(n2-1*dn   ) * h
   w(n2  )      =  3.d0/8.d0  * r(n2        ) * h
   if (n > 6*dn) then
      do i = n1+3*dn,n2-3*dn,dn
         w(i) = r(i) * h
      end do
   end if
99 continue
   end subroutine set_weight_exp2

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
      call diff_exp(ier,isdiff,1,5,rnn,fnn,r,f,df,ddf)
      dfn(i) = df ; ddfn(i) = ddf
   end do
99 continue
   deallocate(rnn,fnn)
   end subroutine calc_ddiff_exp2
