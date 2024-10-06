!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 614 $)
!
!  SUBROUINE:  rsreal, rsintg, tsort, cmpstr, skip_lines, getttr, hpsort,
!              dsjnv, dsjnvn, sbess1, k_plus_G_vectors, G_dot_R, error,
!              pucv2cart, sphr, sphr_diff, eqivvl
!
!  AUTHOR(S): T. Yamasaki, M. Okamoto, H. Sawada   
!  August/20/2003
!  
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
!=======================================================================
!
!     The original version of this set of the computer programs "PHASE"
!  was developed by the members of the Theory Group of Joint Research
!  Center for Atom Technology (JRCAT), based in Tsukuba, in the period
!  1993-2001.
!
!     Since 2002, this set has been tuned and new functions have been
!  added to it as a part of the national project "Frontier Simulation 
!  Software for Industrial Science (FSIS)",  which is supported by
!  the IT program of the Ministry of Education, Culture, Sports,
!  Science and Technology (MEXT) of Japan. 
!     Since 2006, this program set has been developed as a part of the
!  national project "Revolutionary Simulation Software (RSS21)", which
!  is supported by the next-generation IT program of MEXT of Japan.
!   Since 2013, this program set has been further developed centering on PHASE System
!  Consortium.
!   The activity of development of this program set has been supervised by Takahisa Ohno.
!
! $Id: bottom_Subroutines.F90 614 2020-05-07 03:24:24Z jkoga $
!
! ***********************************************************
!
subroutine rsreal(n, obj)
  integer, intent(in) :: n
  integer, parameter :: DP = kind(1.d0)
  real(kind=DP), intent(out), dimension(n) :: obj

  integer :: i
  do i = 1, n
     obj(i) = 0.d0
  enddo
end subroutine rsreal

subroutine rsintg(n, obj)
  integer, intent(in) :: n
  integer, intent(out), dimension(n) :: obj

  integer i
  do i = 1, n
     obj(i) = 0
  enddo
end subroutine rsintg

!!$subroutine gettod(x)
!!$  integer, parameter :: DP = kind(1.d0)
!!$  real(kind=DP), intent(out) :: x
!!$  
!!$  real s, secnds
!!$  s = secnds(0.0)
!!$  x = s*1.d6
!!$end subroutine gettod

!---------------------------------------------
subroutine tsort(ecpu,msbrnm,iorder)
  implicit none
  integer, parameter :: DP = kind(1.d0)
  integer      , intent(in)    :: msbrnm
  real(kind=DP), intent(in), dimension(msbrnm) :: ecpu
  integer      , intent(out),dimension(msbrnm) :: iorder

  integer :: i,ib,ia,it
  do i = 1, msbrnm
     iorder(i) = i
  enddo
  do ib = 2, msbrnm
     do ia = 1, ib-1
        if(ecpu(iorder(ia)).lt.ecpu(iorder(ib))) then
           it = iorder(ia)
           iorder(ia) = iorder(ib)
           iorder(ib) = it
        endif
     enddo
  enddo
end subroutine tsort

subroutine cmpstr(name1,name2,nlen,ident)
  implicit none
  integer, intent(in) :: nlen
  character(len=*), intent(in) :: name1, name2
  logical, intent(out) ::  ident

  character :: c
  integer   :: ndif,i

  ndif = ichar('A') - ichar('a')

  ident = .true.
  do i = 1, nlen
     c = name2(i:i)
     if(ichar(c).ge.ichar('a').and.ichar(c).le.ichar('z')) then
        c = char(ichar(c) + ndif)
     endif
     if(name1(i:i) .ne. c) goto 1001
  enddo
  goto 1002
1001 ident = .false.
1002 continue
end subroutine cmpstr

subroutine skip_lines(nfinp, nline)
  integer, intent(in) :: nfinp
  integer, intent(in) :: nline
  integer :: i
  do i = 1, nline
     read(nfinp,*)
  enddo
end subroutine skip_lines

subroutine getttr (rltv,ttr)
  implicit none
  integer, parameter :: DP = kind(1.d0)
  real(kind=DP), intent(in), dimension(3,3) :: rltv
  real(kind=DP), intent(out),dimension(6)   :: ttr

  integer i, j, jj

!!$  call rsreal(6,ttr)
  ttr = 0.d0
  do j = 1, 3
     jj = mod(j+1,3)
     if(jj.eq.0) jj = 3
     do i = 1, 3
        ttr(j) = ttr(j) + rltv(i,j)*rltv(i,j)
        ttr(j+3) = ttr(j+3) + 2.d0*rltv(i,j)*rltv(i,jj)
     enddo
  enddo
end subroutine getttr

#ifdef LIBRARY_BUILD
subroutine hpsort_phase0(k,n,bxyz,br)
#else
subroutine hpsort(k,n,bxyz,br)
#endif
!=======================================================================
!
!  SOFTWARE NAME : PHASE ($Revision: 614 $)
!
!  SUBROUINE:
!
!  AUTHOR(S): M. Okamoto   August/20/2003
!
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)
!
!=======================================================================
  use m_Const_Parameters, only : DP
  implicit none
  integer,intent(in)          :: k, n
  real(kind=DP),intent(inout) :: bxyz(k,3), br(k)

  integer,allocatable       :: idx_gvec(:)
  real(kind=DP),allocatable :: norm_gvec(:), gvec_tmp(:,:)
  integer :: i
 !+++++++++++++++++++++++++++++++++++++++++++++++++
  allocate(norm_gvec(n),idx_gvec(n),gvec_tmp(k,3))
 !+++++++++++++++++++++++++++++++++++++++++++++++++
  do i = 1,n
     norm_gvec(i) = br(i)
     idx_gvec (i) = i
  end do
  do i = n/2,1,-1
     call sift_down(i,n)
  end do
  do i = n,2,-1
     call swap_i(idx_gvec(1),idx_gvec(i))
     call sift_down(1,i-1)
  end do
  gvec_tmp(:,:) = bxyz(:,:)
  do i = 1,n
     bxyz(i,:) = gvec_tmp(idx_gvec(i),:)
     br(i) = norm_gvec(idx_gvec(i))
  end do
 !++++++++++++++++++++++++++++++++++++++++
  deallocate(norm_gvec,idx_gvec,gvec_tmp)
 !++++++++++++++++++++++++++++++++++++++++
!....................................................................
  contains
  subroutine sift_down(l,r)
     integer,intent(in) :: l, r
     integer :: j, j0, l0
     real(kind=DP) :: v0
     l0 = idx_gvec(l)
     v0 = norm_gvec(l0)
     j0 = l
     j  = l+l
     do while (j <= r)
        if (j < r) then
           if (norm_gvec(idx_gvec(j)) < norm_gvec(idx_gvec(j+1))) then
              j = j+1
           end if
        end if
        if (v0 >= norm_gvec(idx_gvec(j))) then
           exit
        end if
        idx_gvec(j0) = idx_gvec(j)
        j0 = j
        j  = j+j
     end do
     idx_gvec(j0) = l0
  end subroutine sift_down
!....................................................................
#ifdef LIBRARY_BUILD
end subroutine hpsort_phase0
#else
end subroutine hpsort
#endif

#ifdef _RPL_BESSEL_

!#####################################################################
! These subroutines, dsjnv and dsjnvn, are replacements of their 
! original ones. (M.Okamoto, December, 2004)

   subroutine dsjnv(n,k,x,y)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 614 $)
!        
!  SUBROUINE: dsjnv
!     
!  AUTHOR(S): M. Okamoto   December/24/2004
!     
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)
!     
!=======================================================================
   implicit none
   integer,intent(in)  :: n,k
   real(8),intent(in)  :: x(k)
   real(8),intent(out) :: y(k)
   integer :: ik
   real(8) :: fn_js
   do ik = 1,k
      y(ik) = fn_js(n,x(ik))
   end do
   end subroutine dsjnv

   subroutine dsjnvn(n,k,idp,x,y)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 614 $)
!        
!  SUBROUINE: dsjnvn
!     
!  AUTHOR(S): M. Okamoto   December/24/2004
!     
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)
!     
!  FURTHER MODIFICATION: T. Yamasaki, March/01/2006
!   idp is introduced
!
!=======================================================================
! This is the same with dsjnv
!
     implicit none
     integer,intent(in)  :: n,k, idp
     real(8),intent(in)  :: x(k)
     real(8),intent(out) :: y(k)
     integer :: ik
     real(8) :: fn_js
     do ik = 1,k
        y(ik) = fn_js(n,x(ik))
     end do
   end subroutine dsjnvn

   function fn_js(n,x) 
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 614 $)
!
!  SUBROUINE: fn_js
!     
!  AUTHOR(S): M. Okamoto   December/25/2004
!     
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)
!
!=======================================================================
   implicit none
   integer,intent(in) :: n
   real(8),intent(in) :: x
   real(8)            :: fn_js
   real(8) :: cc(8,0:8), s, c, x2
   real(8),parameter :: EPS = 1.0d0
   x2 = x*x ; s = sin(x) ; c = cos(x)
   cc(:,:) = 0.d0
   cc(1,0) =  1.0000000000000000000d0  
   cc(2,0) = -1.6666666666666666667d-1 
   cc(3,0) =  8.3333333333333333333d-3 
   cc(4,0) = -1.9841269841269841270d-4 
   cc(5,0) =  2.7557319223985890653d-6 
   cc(6,0) = -2.5052108385441718775d-8 
   cc(7,0) =  1.6059043836821614599d-10
   cc(8,0) = -7.6471637318198164759d-13
   cc(1,1) =  3.3333333333333333333d-1 
   cc(2,1) = -3.3333333333333333333d-2 
   cc(3,1) =  1.1904761904761904762d-3 
   cc(4,1) = -2.2045855379188712522d-5 
   cc(5,1) =  2.5052108385441718775d-7 
   cc(6,1) = -1.9270852604185937519d-9 
   cc(7,1) =  1.0706029224547743066d-11
   cc(8,1) = -4.4983316069528332211d-14
   cc(1,2) =  6.6666666666666666667d-2 
   cc(2,2) = -4.7619047619047619048d-3 
   cc(3,2) =  1.3227513227513227513d-4 
   cc(4,2) = -2.0041686708353375020d-6 
   cc(5,2) =  1.9270852604185937519d-8 
   cc(6,2) = -1.2847235069457291680d-10
   cc(7,2) =  6.2976642497339665096d-13
   cc(8,2) = -2.3675429510278069585d-15
   cc(1,3) =  9.5238095238095238095d-3 
   cc(2,3) = -5.2910052910052910053d-4 
   cc(3,3) =  1.2025012025012025012d-5 
   cc(4,3) = -1.5416682083348750015d-7 
   cc(5,3) =  1.2847235069457291680d-9 
   cc(6,3) = -7.5571970996807598115d-12
   cc(7,3) =  3.3145601314389297419d-14
   cc(8,3) = -1.1274014052513366469d-16
   cc(1,4) =  1.0582010582010582011d-3 
   cc(2,4) = -4.8100048100048100048d-5 
   cc(3,4) =  9.2500092500092500093d-7 
   cc(4,4) = -1.0277788055565833344d-8 
   cc(5,4) =  7.5571970996807598115d-11
   cc(6,4) = -3.9774721577267156903d-13
   cc(7,4) =  1.5783619673518713057d-15
   cc(8,4) = -4.9017452402232028126d-18
   cc(1,5) =  9.6200096200096200096d-5 
   cc(2,5) = -3.7000037000037000037d-6 
   cc(3,5) =  6.1666728333395000062d-8 
   cc(4,5) = -6.0457576797446078492d-10
   cc(5,5) =  3.9774721577267156903d-12
   cc(6,5) = -1.8940343608222455668d-14
   cc(7,5) =  6.8624433363124839376d-17
   cc(8,5) = -1.9606980960892811250d-19
   cc(1,6) =  7.4000074000074000074d-6 
   cc(2,6) = -2.4666691333358000025d-7 
   cc(3,6) =  3.6274546078467647095d-9 
   cc(4,6) = -3.1819777261813725522d-11
   cc(5,6) =  1.8940343608222455668d-13
   cc(6,6) = -8.2349320035749807252d-16
   cc(7,6) =  2.7449773345249935751d-18
   cc(8,6) = -7.2618448003306708335d-21
#ifdef NONE
   cc(1,7) =  4.9333382666716000049d-7 
   cc(2,7) = -1.4509818431387058838d-8 
   cc(3,7) =  1.9091866357088235313d-10
   cc(4,7) = -1.5152274886577964534d-12
   cc(5,7) =  8.2349320035749807252d-15
   cc(6,7) = -3.2939728014299922901d-17
   cc(7,7) =  1.0166582720462939167d-19
   cc(8,7) = -2.5040844139071278736d-22
   cc(1,8) =  2.9019636862774117676d-8 
   cc(2,8) = -7.6367465428352941253d-10
   cc(3,8) =  9.0913649319467787206d-12
   cc(4,8) = -6.5879456028599845801d-14
   cc(5,8) =  3.2939728014299922901d-16
   cc(6,8) = -1.2199899264555527000d-18
   cc(7,8) =  3.5057181794699790231d-21
   cc(8,8) = -8.0776916577649286246d-24
#endif
   if (abs(x) < EPS) then
#ifdef NONE
! spherical bessels with l=7,8 are not performed in double precision.
      if (abs(n) <= 8) then
#else
      if (abs(n) <= 6) then
#endif
         fn_js = (cc(1,n)+x2*(cc(2,n)+x2*(cc(3,n)+x2*(cc(4,n) &
                +x2*(cc(5,n)+x2*(cc(6,n)+x2*(cc(7,n)+x2*cc(8,n) &
                 )))))))*x**n
      else
         fn_js = 0.d0
      end if
   else
      select case (abs(n))
      case (0)
         fn_js = s/x
      case (1)
         fn_js = (s-x*c)/x2
      case (2)
         fn_js = ((3.d0-x2)*s-3.d0*x*c)/x/x2
      case (3)
         fn_js = ((15.d0-6.d0*x2)*s &
                  -(15.d0-x2)*x*c)/x2/x2
      case (4)
         fn_js = ((105.d0+(-45.d0+x2)*x2)*s &
                  -(105.d0-10.d0*x2)*x*c)/x/x2/x2
      case (5)
         fn_js = ((945.d0+(-420.d0+15.d0*x2)*x2)*s &
                  -(945.d0+(-105.d0+x2)*x2)*x*c)/x2/x2/x2
      case (6)
         fn_js = ((10395.d0+(-4725.d0+(210.d0-x2)*x2)*x2)*s &
                  -(10395.d0+(-1260.d0+21.d0*x2)*x2)*x*c &
                 )/x/x2/x2/x2
#ifdef NONE
! spherical bessels with l=7,8 are not performed in double precision.
      case (7)
         fn_js = ((135135.d0+(-62370.d0+(3150.d0-28.d0*x2)*x2)*x2)*s &
                  -(135135.d0+(-17325.d0+(378.d0-x2)*x2)*x2)*x*c &
                 )/x2/x2/x2/x2
      case (8)
         fn_js = ( &
         (2027025.d0+(-945945.d0+(51975.d0+(-630.d0+x2)*x2)*x2)*x2)*s &
         - (2027025.d0+(-270270.d0+(6930.d0-36.d0*x2)*x2)*x2)*x*c &
                 )/x/x2/x2/x2/x2
#endif
      case default
         fn_js = 0.d0
      end select
   end if
   end function fn_js

#elif _RPL_BESSEL_ABCAP_

!#####################################################################
! These subroutines, dsjnv_abcap and dsjnvn_abcap, are replacements of their
! original ones. (M.Okamoto, August, 2003)
         
   subroutine dsjnv(n,k,x,y)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 614 $)
!
!  SUBROUINE: dsjnv
!
!  AUTHOR(S): M. Okamoto   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

   implicit none
   integer,intent(in)  :: n,k
   real(8),intent(in)  :: x(k)
   real(8),intent(out) :: y(k)
   integer :: ik
   real(8) :: ff(0:20),dd(0:20)
   do ik = 1,k
      call sbess1(n,x(ik),ff(0),dd(0))
      y(ik) = ff(n)
   end do
   end subroutine dsjnv

   subroutine dsjnvn(n,k,idp,x,y)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 614 $)
!
!  SUBROUINE: dsjnvn
!
!  AUTHOR(S): M. Okamoto   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================
   implicit none
   integer,intent(in)  :: n,k,idp
   real(8),intent(in)  :: x(k)
   real(8),intent(out) :: y(k)
   integer :: ik
   real(8) :: ff(0:20),dd(0:20)
   do ik = 1,k
      call sbess1(n,x(ik),ff(0),dd(0))
      y(ik) = ff(n)
   end do
   end subroutine dsjnvn
   
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
!
      subroutine sbess1(lmax0,x,ff,d)
!
      implicit none
      integer,intent(in)  :: lmax0
      real(8),intent(in)  :: x
      real(8),intent(out) :: ff(0:lmax0),d(0:lmax0)
      
      integer,parameter :: lmax=20
      real(8) :: f(0:lmax),q(0:lmax)
      real(8) :: fval,fplus,fsum,fmin
      integer :: l
      real(8) :: xn,f1,f2,xl
!
      if(lmax0 > lmax) then
        write(6,*) ' lmax0=',lmax0,'   lmax=',lmax
        call phase_error_with_msg(6, '=== stop in sub.sbess1. (lmax0) ===',__LINE__,__FILE__)
      end if
!                                                 ----- x<0.02 -----
      if(x < 0.02) then
        xn=1.0
        f1=1.0
        f2=3.0
        f(0)=xn/f1-0.5*xn*x*x/f2
        do 10 l=1,lmax
          f1=f2
          f2=f2*(2*l+3)
          q(l)=xn/f1-0.5*xn*x*x/f2
          f(l)=q(l)*x
          xn=xn*x
   10   continue
!
!                                                 ----- x<15.   -----
      else if(x < 15.) then
        fval=1.e-30
        fplus=0.0
        fsum=0.0
        do 20 l=lmax+20,lmax+2,-1
          fmin=fval*(2*l+1)/x-fplus
          fsum=fsum+(2*l-1)*fmin*fmin
          fplus=fval
          fval=fmin
          if(abs(fmin).gt.1.e10) then
            fsum=1./sqrt(fsum)
            fplus=fplus*fsum
            fval=fval*fsum
            fsum=1.0
          end if
   20   continue
        if(x > 1.0) then
          xl=1.0
        else
          xl=x**lmax
        end if
        fsum=xl/sqrt(fsum)
        fplus=fplus*fsum
        fval=fval*fsum
        fsum=xl*xl
        do 30 l=lmax+1,1,-1
          fmin=fval*(2*l+1)/x-fplus
          fsum=fsum+(2*l-1)*fmin*fmin
          fplus=fval
          fval=fmin
          f(l-1)=fval
   30   continue
        fsum=1.0/sqrt(fsum)
        do 40 l=0,lmax
          f(l)=f(l)*fsum
          q(l)=f(l)/x
   40   continue
!
!                                                 ----- x>15.   -----
      else
        f(0)=sin(x)/x
        f(1)=sin(x)/(x*x)-cos(x)/x
        q(1)=f(1)/x
        do 50 l=1,lmax-1
          f(l+1)=(2*l+1)*f(l)/x-f(l-1)
          q(l+1)=f(l+1)/x
   50   continue
      end if
!
      do 55 l=0,lmax0
   55 ff(l)=f(l)
!                                                ----- derivative ----
      d(0)=-f(1)
      do 60 l=1,lmax0
   60 d(l)=f(l-1)-(l+1)*q(l)
!
      end subroutine sbess1
#else


! Vector version
   subroutine dsjnv(n,k,x,y)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 614 $)
!        
!  SUBROUINE: dsjnv
!     
!  AUTHOR(S): M. Okamoto and T. Yamamoto   December/24/2004
!     
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)
!     
!=======================================================================
   implicit none
   integer,intent(in)  :: n,k
   real(8),intent(in)  :: x(k)
   real(8),intent(out) :: y(k)
   integer :: ik

   real(8) :: x2,x4
   real(8),parameter :: EPS = 1.0d0
   real(8),parameter :: &
         & cc10 =  1.0000000000000000000d0, &
         & cc20 = -1.6666666666666666667d-1, &
         & cc30 =  8.3333333333333333333d-3, &
         & cc40 = -1.9841269841269841270d-4, &
         & cc50 =  2.7557319223985890653d-6, & 
         & cc60 = -2.5052108385441718775d-8, &
         & cc70 =  1.6059043836821614599d-10, &
         & cc80 = -7.6471637318198164759d-13
   real(8),parameter :: &
         & cc11 =  3.3333333333333333333d-1, &
         & cc21 = -3.3333333333333333333d-2, &
         & cc31 =  1.1904761904761904762d-3, &
         & cc41 = -2.2045855379188712522d-5, &
         & cc51 =  2.5052108385441718775d-7, &
         & cc61 = -1.9270852604185937519d-9, &
         & cc71 =  1.0706029224547743066d-11, &
         & cc81 = -4.4983316069528332211d-14
   real(8),parameter :: &
         & cc12 =  6.6666666666666666667d-2, &
         & cc22 = -4.7619047619047619048d-3, &
         & cc32 =  1.3227513227513227513d-4, &
         & cc42 = -2.0041686708353375020d-6, &
         & cc52 =  1.9270852604185937519d-8, &
         & cc62 = -1.2847235069457291680d-10, &
         & cc72 =  6.2976642497339665096d-13, &
         & cc82 = -2.3675429510278069585d-15
   real(8),parameter :: &
         & cc13 =  9.5238095238095238095d-3, &
         & cc23 = -5.2910052910052910053d-4, &
         & cc33 =  1.2025012025012025012d-5, &
         & cc43 = -1.5416682083348750015d-7, &
         & cc53 =  1.2847235069457291680d-9, &
         & cc63 = -7.5571970996807598115d-12, &
         & cc73 =  3.3145601314389297419d-14, &
         & cc83 = -1.1274014052513366469d-16
   real(8),parameter :: &
         & cc14 =  1.0582010582010582011d-3, &
         & cc24 = -4.8100048100048100048d-5, &
         & cc34 =  9.2500092500092500093d-7, &
         & cc44 = -1.0277788055565833344d-8, &
         & cc54 =  7.5571970996807598115d-11, &
         & cc64 = -3.9774721577267156903d-13, &
         & cc74 =  1.5783619673518713057d-15, &
         & cc84 = -4.9017452402232028126d-18
   real(8),parameter :: &
         & cc15 =  9.6200096200096200096d-5, &
         & cc25 = -3.7000037000037000037d-6, &
         & cc35 =  6.1666728333395000062d-8, &
         & cc45 = -6.0457576797446078492d-10, &
         & cc55 =  3.9774721577267156903d-12, &
         & cc65 = -1.8940343608222455668d-14, &
         & cc75 =  6.8624433363124839376d-17, &
         & cc85 = -1.9606980960892811250d-19
   real(8),parameter :: &
         & cc16 =  7.4000074000074000074d-6, &
         & cc26 = -2.4666691333358000025d-7, &
         & cc36 =  3.6274546078467647095d-9, &
         & cc46 = -3.1819777261813725522d-11, &
         & cc56 =  1.8940343608222455668d-13, &
         & cc66 = -8.2349320035749807252d-16, &
         & cc76 =  2.7449773345249935751d-18, &
#ifdef NONE
         & cc86 = -7.2618448003306708335d-21
   real(8),parameter :: &
         & cc17 =  4.9333382666716000049d-7, &
         & cc27 = -1.4509818431387058838d-8, &
         & cc37 =  1.9091866357088235313d-10, &
         & cc47 = -1.5152274886577964534d-12, &
         & cc57 =  8.2349320035749807252d-15, &
         & cc67 = -3.2939728014299922901d-17, &
         & cc77 =  1.0166582720462939167d-19, &
         & cc87 = -2.5040844139071278736d-22
   real(8),parameter :: &
         & cc18 =  2.9019636862774117676d-8, &
         & cc28 = -7.6367465428352941253d-10, &
         & cc38 =  9.0913649319467787206d-12, &
         & cc48 = -6.5879456028599845801d-14, &
         & cc58 =  3.2939728014299922901d-16, &
         & cc68 = -1.2199899264555527000d-18, &
         & cc78 =  3.5057181794699790231d-21, &
         & cc88 = -8.0776916577649286246d-24
#else
         & cc86 = -7.2618448003306708335d-21
#endif

   select case (abs(n))
   case (0)
      do ik=1,k
         if (abs(x(ik)) < EPS) then
            x2=x(ik)*x(ik)
            y(ik) = (cc10+x2*(cc20+x2*(cc30+x2*(cc40 &
                +x2*(cc50+x2*(cc60+x2*(cc70+x2*cc80 &
                 )))))))
         else
            y(ik) = sin(x(ik))/x(ik)
         end if
      end do
   case (1)
      do ik=1,k
         x2=x(ik)*x(ik)
         if (abs(x(ik)) < EPS) then
            y(ik) = (cc11+x2*(cc21+x2*(cc31+x2*(cc41 &
                +x2*(cc51+x2*(cc61+x2*(cc71+x2*cc81 &
                 )))))))*x(ik)
         else
            y(ik) = (sin(x(ik))-x(ik)*cos(x(ik)))/x2
         end if
      end do
   case (2)
      do ik=1,k
         x2=x(ik)*x(ik)
         if (abs(x(ik)) < EPS) then
            y(ik) = (cc12+x2*(cc22+x2*(cc32+x2*(cc42 &
                +x2*(cc52+x2*(cc62+x2*(cc72+x2*cc82 &
                 )))))))*x2
         else
            y(ik) = ((3.d0-x2)*sin(x(ik))-3.d0*x(ik)*cos(x(ik)))/(x(ik)*x2)
         end if
      end do
   case (3)
      do ik=1,k
         x2=x(ik)*x(ik)
         if (abs(x(ik)) < EPS) then
            y(ik) = (cc13+x2*(cc23+x2*(cc33+x2*(cc43 &
                +x2*(cc53+x2*(cc63+x2*(cc73+x2*cc83 &
                 )))))))*x2*x(ik)
         else
            y(ik) = ((15.d0-6.d0*x2)*sin(x(ik)) &
                  -(15.d0-x2)*x(ik)*cos(x(ik)))/(x2*x2)
         end if
      end do
   case (4)
      do ik=1,k
         x2=x(ik)*x(ik)
         if (abs(x(ik)) < EPS) then
            y(ik) = (cc14+x2*(cc24+x2*(cc34+x2*(cc44 &
                +x2*(cc54+x2*(cc64+x2*(cc74+x2*cc84 &
                 )))))))*x2*x2
         else
            y(ik) = ((105.d0+(-45.d0+x2)*x2)*sin(x(ik)) &
                  -(105.d0-10.d0*x2)*x(ik)*cos(x(ik)))/(x(ik)*x2*x2)
         end if
      end do
   case (5)
      do ik=1,k
         x2=x(ik)*x(ik)
         if (abs(x(ik)) < EPS) then
            y(ik) = (cc15+x2*(cc25+x2*(cc35+x2*(cc45 &
                +x2*(cc55+x2*(cc65+x2*(cc75+x2*cc85 &
                 )))))))*x2*x2*x(ik)
         else
            y(ik) = ((945.d0+(-420.d0+15.d0*x2)*x2)*sin(x(ik)) &
                  -(945.d0+(-105.d0+x2)*x2)*x(ik)*cos(x(ik)))/(x2*x2*x2)
         end if
      end do
   case (6)
      do ik=1,k
         x2=x(ik)*x(ik)
         if (abs(x(ik)) < EPS) then
            y(ik) = (cc16+x2*(cc26+x2*(cc36+x2*(cc46 &
                +x2*(cc56+x2*(cc66+x2*(cc76+x2*cc86 &
                 )))))))*x2*x2*x2
         else
            y(ik) = ((10395.d0+(-4725.d0+(210.d0-x2)*x2)*x2)*sin(x(ik)) &
                  -(10395.d0+(-1260.d0+21.d0*x2)*x2)*x(ik)*cos(x(ik)) &
                 )/(x(ik)*x2*x2*x2)
         end if
      end do
#ifdef NONE
! spherical bessels with l=7,8 are not performed in double precision.
   case (7)
      do ik=1,k
         x2=x(ik)*x(ik)
         if (abs(x(ik)) < EPS) then
            y(ik) = (cc17+x2*(cc27+x2*(cc37+x2*(cc47 &
                +x2*(cc57+x2*(cc67+x2*(cc77+x2*cc87 &
                 )))))))*x2*x2*x2*x(ik)
         else
            x4=x2*x2
            y(ik) = ((135135.d0+(-62370.d0+(3150.d0-28.d0*x2)*x2)*x2)*sin(x(ik)) &
                  -(135135.d0+(-17325.d0+(378.d0-x2)*x2)*x2)*x(ik)*cos(x(ik)) &
                 )/(x4*x4)
         end if
      end do
   case (8)
      do ik=1,k
         x2=x(ik)*x(ik)
         x4=x2*x2
         if (abs(x(ik)) < EPS) then
            y(ik) = (cc18+x2*(cc28+x2*(cc38+x2*(cc48 &
                +x2*(cc58+x2*(cc68+x2*(cc78+x2*cc88 &
                 )))))))*x4*x4
         else
            y(ik) = ( &
         (2027025.d0+(-945945.d0+(51975.d0+(-630.d0+x2)*x2)*x2)*x2)*sin(x(ik)) &
         - (2027025.d0+(-270270.d0+(6930.d0-36.d0*x2)*x2)*x2)*x(ik)*cos(x(ik)) &
                 )/(x(ik)*x4*x4)
         end if
      end do
#endif
   case default
      y(1:k) = 0.d0
   end select

   end subroutine dsjnv

   subroutine dsjnvn(n,k,idp,x,y)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 614 $)
!        
!  SUBROUINE: dsjnvn
!     
!  AUTHOR(S): M. Okamoto and T. Yamamoto   December/24/2004
!     
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)
!
!  FURTHER MODIFICATION: T. Yamasaki, March/01/2006
!   idp is introduced
!     
!=======================================================================
   implicit none
   integer,intent(in)  :: n,k,idp
   real(8),intent(in)  :: x(k)
   real(8),intent(out) :: y(k)
   integer :: ik, i

   real(8) :: x2,x4
   real(8),parameter :: EPS = 1.0d0
   real(8),parameter :: &
         & cc10 =  1.0000000000000000000d0, &
         & cc20 = -1.6666666666666666667d-1, &
         & cc30 =  8.3333333333333333333d-3, &
         & cc40 = -1.9841269841269841270d-4, &
         & cc50 =  2.7557319223985890653d-6, & 
         & cc60 = -2.5052108385441718775d-8, &
         & cc70 =  1.6059043836821614599d-10, &
         & cc80 = -7.6471637318198164759d-13
   real(8),parameter :: &
         & cc11 =  3.3333333333333333333d-1, &
         & cc21 = -3.3333333333333333333d-2, &
         & cc31 =  1.1904761904761904762d-3, &
         & cc41 = -2.2045855379188712522d-5, &
         & cc51 =  2.5052108385441718775d-7, &
         & cc61 = -1.9270852604185937519d-9, &
         & cc71 =  1.0706029224547743066d-11, &
         & cc81 = -4.4983316069528332211d-14
   real(8),parameter :: &
         & cc12 =  6.6666666666666666667d-2, &
         & cc22 = -4.7619047619047619048d-3, &
         & cc32 =  1.3227513227513227513d-4, &
         & cc42 = -2.0041686708353375020d-6, &
         & cc52 =  1.9270852604185937519d-8, &
         & cc62 = -1.2847235069457291680d-10, &
         & cc72 =  6.2976642497339665096d-13, &
         & cc82 = -2.3675429510278069585d-15
   real(8),parameter :: &
         & cc13 =  9.5238095238095238095d-3, &
         & cc23 = -5.2910052910052910053d-4, &
         & cc33 =  1.2025012025012025012d-5, &
         & cc43 = -1.5416682083348750015d-7, &
         & cc53 =  1.2847235069457291680d-9, &
         & cc63 = -7.5571970996807598115d-12, &
         & cc73 =  3.3145601314389297419d-14, &
         & cc83 = -1.1274014052513366469d-16
   real(8),parameter :: &
         & cc14 =  1.0582010582010582011d-3, &
         & cc24 = -4.8100048100048100048d-5, &
         & cc34 =  9.2500092500092500093d-7, &
         & cc44 = -1.0277788055565833344d-8, &
         & cc54 =  7.5571970996807598115d-11, &
         & cc64 = -3.9774721577267156903d-13, &
         & cc74 =  1.5783619673518713057d-15, &
         & cc84 = -4.9017452402232028126d-18
   real(8),parameter :: &
         & cc15 =  9.6200096200096200096d-5, &
         & cc25 = -3.7000037000037000037d-6, &
         & cc35 =  6.1666728333395000062d-8, &
         & cc45 = -6.0457576797446078492d-10, &
         & cc55 =  3.9774721577267156903d-12, &
         & cc65 = -1.8940343608222455668d-14, &
         & cc75 =  6.8624433363124839376d-17, &
         & cc85 = -1.9606980960892811250d-19
   real(8),parameter :: &
         & cc16 =  7.4000074000074000074d-6, &
         & cc26 = -2.4666691333358000025d-7, &
         & cc36 =  3.6274546078467647095d-9, &
         & cc46 = -3.1819777261813725522d-11, &
         & cc56 =  1.8940343608222455668d-13, &
         & cc66 = -8.2349320035749807252d-16, &
         & cc76 =  2.7449773345249935751d-18, &
#ifdef NONE
         & cc86 = -7.2618448003306708335d-21
   real(8),parameter :: &
         & cc17 =  4.9333382666716000049d-7, &
         & cc27 = -1.4509818431387058838d-8, &
         & cc37 =  1.9091866357088235313d-10, &
         & cc47 = -1.5152274886577964534d-12, &
         & cc57 =  8.2349320035749807252d-15, &
         & cc67 = -3.2939728014299922901d-17, &
         & cc77 =  1.0166582720462939167d-19, &
         & cc87 = -2.5040844139071278736d-22
   real(8),parameter :: &
         & cc18 =  2.9019636862774117676d-8, &
         & cc28 = -7.6367465428352941253d-10, &
         & cc38 =  9.0913649319467787206d-12, &
         & cc48 = -6.5879456028599845801d-14, &
         & cc58 =  3.2939728014299922901d-16, &
         & cc68 = -1.2199899264555527000d-18, &
         & cc78 =  3.5057181794699790231d-21, &
         & cc88 = -8.0776916577649286246d-24
#else
         & cc86 = -7.2618448003306708335d-21
#endif

   select case (abs(n))
   case (0)
      do ik=1,idp-1
         x2=x(ik)*x(ik)
         y(ik) = (cc10+x2*(cc20+x2*(cc30+x2*(cc40 &
              & +x2*(cc50+x2*(cc60+x2*(cc70+x2*cc80 &
              & )))))))
      end do
      do ik = idp, k
         y(ik) = sin(x(ik))/x(ik)
      end do
   case (1)
      do ik=1,idp-1
         x2=x(ik)*x(ik)
         y(ik) = (cc11+x2*(cc21+x2*(cc31+x2*(cc41 &
              & +x2*(cc51+x2*(cc61+x2*(cc71+x2*cc81 &
              & )))))))*x(ik)
      end do
      do ik=idp,k
         x2=x(ik)*x(ik)
         y(ik) = (sin(x(ik))-x(ik)*cos(x(ik)))/x2
      end do
   case (2)
      do ik=1,idp-1
         x2=x(ik)*x(ik)
         y(ik) = (cc12+x2*(cc22+x2*(cc32+x2*(cc42 &
              & +x2*(cc52+x2*(cc62+x2*(cc72+x2*cc82 &
              & )))))))*x2
      end do
      do ik=idp,k
         x2=x(ik)*x(ik)
         y(ik) = ((3.d0-x2)*sin(x(ik))-3.d0*x(ik)*cos(x(ik)))/(x(ik)*x2)
      end do
   case (3)
      do ik=1,idp-1
         x2=x(ik)*x(ik)
         y(ik) = (cc13+x2*(cc23+x2*(cc33+x2*(cc43 &
              & +x2*(cc53+x2*(cc63+x2*(cc73+x2*cc83 &
              & )))))))*x2*x(ik)
      end do
      do ik=idp,k
         x2=x(ik)*x(ik)
         y(ik) = ((15.d0-6.d0*x2)*sin(x(ik)) &
              -(15.d0-x2)*x(ik)*cos(x(ik)))/(x2*x2)
      end do
   case (4)
      do ik=1,idp-1
         x2=x(ik)*x(ik)
         y(ik) = (cc14+x2*(cc24+x2*(cc34+x2*(cc44 &
              & +x2*(cc54+x2*(cc64+x2*(cc74+x2*cc84 &
              & )))))))*x2*x2
      end do
      do ik=idp,k
         x2=x(ik)*x(ik)
         y(ik) = ((105.d0+(-45.d0+x2)*x2)*sin(x(ik)) &
              & -(105.d0-10.d0*x2)*x(ik)*cos(x(ik)))/(x(ik)*x2*x2)
      end do
   case (5)
      do ik=1,idp-1
         x2=x(ik)*x(ik)
         y(ik) = (cc15+x2*(cc25+x2*(cc35+x2*(cc45 &
              & +x2*(cc55+x2*(cc65+x2*(cc75+x2*cc85 &
              & )))))))*x2*x2*x(ik)
      end do
      do ik=idp,k
         x2=x(ik)*x(ik)
         y(ik) = ((945.d0+(-420.d0+15.d0*x2)*x2)*sin(x(ik)) &
              & -(945.d0+(-105.d0+x2)*x2)*x(ik)*cos(x(ik)))/(x2*x2*x2)
      end do
   case (6)
      do ik=1,idp-1
         x2=x(ik)*x(ik)
         y(ik) = (cc16+x2*(cc26+x2*(cc36+x2*(cc46 &
              & +x2*(cc56+x2*(cc66+x2*(cc76+x2*cc86 &
              & )))))))*x2*x2*x2
      end do
      do ik=idp,k
         x2=x(ik)*x(ik)
         y(ik) = ((10395.d0+(-4725.d0+(210.d0-x2)*x2)*x2)*sin(x(ik)) &
              & -(10395.d0+(-1260.d0+21.d0*x2)*x2)*x(ik)*cos(x(ik)) &
              & )/(x(ik)*x2*x2*x2)
      end do
#ifdef NONE
! spherical bessels with l=7,8 are not performed in double precision.
   case (7)
      do ik=1,idp-1
         x2=x(ik)*x(ik)
         y(ik) = (cc17+x2*(cc27+x2*(cc37+x2*(cc47 &
              & +x2*(cc57+x2*(cc67+x2*(cc77+x2*cc87 &
              & )))))))*x2*x2*x2*x(ik)
      end do
      do ik=idp,k
         x2=x(ik)*x(ik)
         x4=x2*x2
         y(ik) = ((135135.d0+(-62370.d0+(3150.d0-28.d0*x2)*x2)*x2)*sin(x(ik)) &
              & -(135135.d0+(-17325.d0+(378.d0-x2)*x2)*x2)*x(ik)*cos(x(ik)) &
              & )/(x4*x4)
      end do
   case (8)
      do ik=1,idp-1
         x2=x(ik)*x(ik)
         x4=x2*x2
         y(ik) = (cc18+x2*(cc28+x2*(cc38+x2*(cc48 &
              & +x2*(cc58+x2*(cc68+x2*(cc78+x2*cc88 &
              & )))))))*x4*x4
      end do
      do ik=idp,k
         x2=x(ik)*x(ik)
         x4=x2*x2
         y(ik) = ( &
              & (2027025.d0+(-945945.d0+(51975.d0+(-630.d0+x2)*x2)*x2)*x2)*sin(x(ik)) &
              & - (2027025.d0+(-270270.d0+(6930.d0-36.d0*x2)*x2)*x2)*x(ik)*cos(x(ik)) &
              & )/(x(ik)*x4*x4)
      end do
#endif
   case default
      y(1:k) = 0.d0
   end select

 end subroutine dsjnvn

#endif

!$$#ifndef PARA3D
subroutine k_plus_G_vectors(ik,kgp,kg1,knv3,iba,nbase,vk,ngabc,rltv&
     &,qx,qy,qz,vlen)
  use m_Const_Parameters, only : DP, CRDTYP, BUCS
  implicit none
  integer, intent(in)        :: ik, kgp,kg1,knv3,iba(knv3),nbase(kg1,knv3)
  real(kind=DP), intent(in)  :: vk(knv3,3,CRDTYP)
  integer, intent(in)        :: ngabc(kgp,3)
  real(kind=DP), intent(in)  :: rltv(3,3)
  real(kind=DP), intent(out) :: qx(kg1),qy(kg1),qz(kg1),vlen(kg1)

  integer i, ip
  real(kind=DP) :: ga, gb, gc
  real(kind=DP), dimension(3) :: kxyz
  kxyz(1:3) = vk(ik,1:3,BUCS)
  !!$do i = 1, 3
  !!$   if(kxyz(i) <= -1.d0) kxyz(i) = kxyz(i) + ceiling(-kxyz(i))
  !!$   if(kxyz(i) >=  1.d0) kxyz(i) = kxyz(i) - floor(kxyz(i))
  !!$end do

  do i = 1, iba(ik)
     ip = nbase(i,ik)
     if (ip<1) cycle
     ga = kxyz(1) + ngabc(ip,1)
     gb = kxyz(2) + ngabc(ip,2)
     gc = kxyz(3) + ngabc(ip,3)
     qx(i)  = rltv(1,1)*ga + rltv(1,2)*gb + rltv(1,3)*gc
     qy(i)  = rltv(2,1)*ga + rltv(2,2)*gb + rltv(2,3)*gc
     qz(i)  = rltv(3,1)*ga + rltv(3,2)*gb + rltv(3,3)*gc
     vlen(i) = dsqrt( qx(i)**2 + qy(i)**2 + qz(i)**2 )
  end do

end subroutine k_plus_G_vectors
!$$#endif

! === KT_add === 2014/09/24
!                                   moved From m_Epsilon_ek.F90
subroutine k_plus_G_vectors_m(ik,kgp,kg1,knv3,iba,nbase,vkxyz,ngabc,rltv, &
     &                        qx,qy,qz)
  !
!  modifed subroutine k_plus_G vectors
!  T. Hamada (Univ. Tokyo) 2003.03.07
!
  use m_Const_Parameters, only : DP, CRDTYP, BUCS
  implicit none
  integer, intent(in)        :: ik, kgp,kg1,knv3,iba(knv3),nbase(kg1,knv3)
  integer, intent(in)        :: ngabc(kgp,3)
  real(kind=DP), intent(in)  :: rltv(3,3)
  real(kind=DP), intent(in)  :: vkxyz(knv3,3,CRDTYP)
  real(kind=DP), intent(out) :: qx(kg1),qy(kg1),qz(kg1)
  integer :: i, ip
  real(kind=DP) :: ga, gb, gc

!!$  write(6,*) ' !! iba(ik) = ', iba(ik)
!!$  write(6,*) ' -- vk --'
!!$  write(6,'(3f20.10)') (vk(ik,i,BUCS),i=1,3)
!!$  write(6,*) ' -- nbase --'
!!$  write(6,'(15i5)') (nbase(i,ik),i=1,100)
  do i = 1, iba(ik)
     ip = nbase(i,ik)
     if (ip<1) cycle
     ga = vkxyz(ik,1,BUCS) + real(ngabc(ip,1),kind=DP)
     gb = vkxyz(ik,2,BUCS) + real(ngabc(ip,2),kind=DP)
     gc = vkxyz(ik,3,BUCS) + real(ngabc(ip,3),kind=DP)
     qx(i)  = rltv(1,1)*ga + rltv(1,2)*gb + rltv(1,3)*gc
     qy(i)  = rltv(2,1)*ga + rltv(2,2)*gb + rltv(2,3)*gc
     qz(i)  = rltv(3,1)*ga + rltv(3,2)*gb + rltv(3,3)*gc
  end do
end subroutine k_plus_G_vectors_m
! ============== 2014/09/24

subroutine G_vectors(kgp,kg,nbase,ngabc,rltv,qx,qy,qz,vlen)
  use m_Const_Parameters, only : DP, CRDTYP, BUCS
  implicit none
  integer, intent(in)        :: kgp,kg,nbase(kg)
  integer, intent(in)        :: ngabc(kgp,3)
  real(kind=DP), intent(in)  :: rltv(3,3)
  real(kind=DP), intent(out) :: qx(kg),qy(kg),qz(kg),vlen(kg)

  integer i, ip
  real(kind=DP) :: ga, gb, gc

  do i = 1, kg
     ip = nbase(i)
     if (ip<1) cycle
     ga = ngabc(ip,1)
     gb = ngabc(ip,2)
     gc = ngabc(ip,3)
     qx(i)  = rltv(1,1)*ga + rltv(1,2)*gb + rltv(1,3)*gc
     qy(i)  = rltv(2,1)*ga + rltv(2,2)*gb + rltv(2,3)*gc
     qz(i)  = rltv(3,1)*ga + rltv(3,2)*gb + rltv(3,3)*gc
     vlen(i) = dsqrt( qx(i)**2 + qy(i)**2 + qz(i)**2 )
  end do
end subroutine G_vectors

subroutine Rotated_k_plus_G_vectors(ik,kgp,kg1,knv3,iba,nbase,vk,ngabc,rltv&
     &,qx,qy,qz,vlen,opr)
  use m_Const_Parameters, only : DP, CRDTYP, BUCS
  implicit none
  integer, intent(in)        :: ik, kgp,kg1,knv3,iba(knv3),nbase(kg1,knv3)
  real(kind=DP), intent(in)  :: vk(knv3,3,CRDTYP)
  integer, intent(in)        :: ngabc(kgp,3)
  real(kind=DP), intent(in)  :: rltv(3,3)
  real(kind=DP), intent(out) :: qx(kg1),qy(kg1),qz(kg1),vlen(kg1)
  real(kind=DP), intent(in)  :: opr(3,3)

  integer :: i,ip
  real(kind=DP) :: ga,gb,gc,gx,gy,gz

  do i = 1, iba(ik)
     ip = nbase(i,ik)
     if (ip<1) cycle
     ga = vk(ik,1,BUCS) + ngabc(ip,1)
     gb = vk(ik,2,BUCS) + ngabc(ip,2)
     gc = vk(ik,3,BUCS) + ngabc(ip,3)
     gx  = rltv(1,1)*ga + rltv(1,2)*gb + rltv(1,3)*gc
     gy  = rltv(2,1)*ga + rltv(2,2)*gb + rltv(2,3)*gc
     gz  = rltv(3,1)*ga + rltv(3,2)*gb + rltv(3,3)*gc
     qx(i) = opr(1,1)*gx+opr(1,2)*gy+opr(1,3)*gz
     qy(i) = opr(2,1)*gx+opr(2,2)*gy+opr(2,3)*gz
     qz(i) = opr(3,1)*gx+opr(3,2)*gy+opr(3,3)*gz
     vlen(i) = dsqrt( qx(i)**2 + qy(i)**2 + qz(i)**2 )
  end do

end subroutine Rotated_k_plus_G_vectors

subroutine G_dot_R(katm,ia,pos,kgp,nbmx,ngabc,zfcos,zfsin)
  use m_Const_Parameters, only : PAI2, DP
  implicit none
  integer, intent(in)          :: katm,ia, kgp, nbmx
  real(kind=DP),intent(in), dimension(katm,3) :: pos
  integer, intent(in),      dimension(kgp,3)  :: ngabc
  real(kind=DP), intent(out), dimension(nbmx) :: zfcos, zfsin

  integer       :: i
  real(kind=DP) :: grt
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
  do i = 1, nbmx
     grt = (pos(ia,1)*ngabc(i,1) + pos(ia,2)*ngabc(i,2)&
          & + pos(ia,3)*ngabc(i,3))*PAI2
     zfcos(i) = dcos(grt)
     zfsin(i) = dsin(grt)
  end do
end subroutine G_dot_R

subroutine Rotated_G_dot_R(katm,ia,cps,kgp,nbmx,ngabc,rltv,opr,zfcos,zfsin)
  use m_Const_Parameters, only : PAI2, DP
  implicit none
  integer, intent(in)          :: katm,ia, kgp, nbmx
  real(kind=DP),intent(in), dimension(katm,3) :: cps
  integer, intent(in),      dimension(kgp,3)  :: ngabc
  real(kind=DP), intent(in), dimension(3,3)  :: rltv
  real(kind=DP), intent(in)  :: opr(3,3)
  real(kind=DP), intent(out), dimension(nbmx) :: zfcos, zfsin

  integer :: i
  real(kind=DP) :: ga,gb,gc,qx,qy,qz,gx,gy,gz,grt

#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
  do i = 1, nbmx
     ga = ngabc(i,1)
     gb = ngabc(i,2)
     gc = ngabc(i,3)
     qx = rltv(1,1)*ga + rltv(1,2)*gb + rltv(1,3)*gc
     qy = rltv(2,1)*ga + rltv(2,2)*gb + rltv(2,3)*gc
     qz = rltv(3,1)*ga + rltv(3,2)*gb + rltv(3,3)*gc
     gx = opr(1,1)*qx+opr(1,2)*qy+opr(1,3)*qz
     gy = opr(2,1)*qx+opr(2,2)*qy+opr(2,3)*qz
     gz = opr(3,1)*qx+opr(3,2)*qy+opr(3,3)*qz
     grt = cps(ia,1)*gx + cps(ia,2)*gy + cps(ia,3)*gz
     zfcos(i) = dcos(grt)
     zfsin(i) = dsin(grt)
  end do

end subroutine Rotated_G_dot_R

subroutine error
  write(6,*) ' ! error '
end subroutine error

subroutine pucv2cart(rltv,x,y,z)
  use m_Const_Parameters, only :DP
  implicit none
  real(kind=DP),intent(in)    :: rltv(3,3)
  real(kind=DP),intent(inout) :: x,y,z
  real(kind=DP) :: tx,ty,tz

  tx = rltv(1,1)*x + rltv(1,2)*y + rltv(1,3)*z
  ty = rltv(2,1)*x + rltv(2,2)*y + rltv(2,3)*z
  tz = rltv(3,1)*x + rltv(3,2)*y + rltv(3,3)*z

  x = tx
  y = ty
  z = tz
end subroutine pucv2cart

!!$C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
!!$C spherical harmonics.
!!$c     Some lines are revised by T.Yamasaki according to Y.Morikawa's
!!$c   e-mail. at 4th May 1994
!!$C*************************************************************
subroutine sphr(kgp, is, gx,gy,gz,ylm)
  use m_Const_Parameters, only  : DP, PAI, PAI4
  implicit none
  integer, intent(in)          :: kgp,is
  real(kind=DP), dimension(kgp), intent(in)  :: gx,gy,gz
  real(kind=DP), dimension(kgp), intent(out) :: ylm

  real(kind=DP) :: g_minimum = 1.d-20
  integer       :: i,ni
  real(kind=DP) :: a,b,c,d,e,f
  if(gx(1)**2+gy(1)**2+gz(1)**2 < g_minimum) then
     ni = 2
     ylm(1) = dsqrt(1.d0/PAI4)
  else
     ni = 1
  end if
  if(is == 1) then
     a = dsqrt(1.d0/PAI4)
     do i = ni, kgp
        ylm(i) = a
     end do
  else if(is == 2) then
     a = dsqrt(3.d0/PAI4)
     do i = ni,kgp
        ylm(i) = a*gx(i)/(dsqrt(gx(i)**2+gy(i)**2+gz(i)**2))
     end do
  else if(is == 3) then
     a = dsqrt(3.d0/PAI4)
     do i = ni,kgp
        ylm(i) = a*gy(i)/(dsqrt(gx(i)**2+gy(i)**2+gz(i)**2))
     end do
  else if(is == 4) then
     a = dsqrt(3.d0/PAI4)
     do i = ni,kgp
        ylm(i) = a*gz(i)/(dsqrt(gx(i)**2+gy(i)**2+gz(i)**2))
     end do
  else if(is == 5) then
     a = dsqrt(5.d0/(16*PAI))
     do i = ni,kgp
        d = gz(i)**2
        e = gx(i)**2+gy(i)**2+d
        ylm(i) = a*(3*d-e)/e
     end do
  else if(is == 6) then
     a = dsqrt(15.d0/(16*PAI))
     do i = ni,kgp
        b = gx(i)**2
        c = gy(i)**2
        e = b+c+gz(i)**2
        ylm(i) = a*(b-c)/e
     end do
  else if(is == 7) then
     a = dsqrt(15.d0/PAI4)
     do i = ni,kgp
        e = gx(i)**2+gy(i)**2+gz(i)**2
        ylm(i) = a*gx(i)*gy(i)/e
     end do
  else if(is == 8) then
     a = dsqrt(15.d0/PAI4)
     do i = ni,kgp
        e = gx(i)**2+gy(i)**2+gz(i)**2
        ylm(i) = a*gy(i)*gz(i)/e
     end do
  else if(is == 9) then
     a = dsqrt(15.d0/PAI4)
     do i = ni,kgp
        e = gx(i)**2+gy(i)**2+gz(i)**2
        ylm(i) = a*gz(i)*gx(i)/e
     end do
  else if(is == 10) then
     a = dsqrt(7.d0/(16*PAI))
     do i = ni,kgp
        d = gz(i)**2
        e = gx(i)**2+gy(i)**2+d
        f = (dsqrt(e))**3
        ylm(i) = a*gz(i)*(5.d0*d-3.d0*e)/f
     end do
  else if(is == 11) then
     a = dsqrt(21.d0/(32*PAI))
     do i = ni,kgp
        d = gz(i)**2
        e = gx(i)**2+gy(i)**2+d
        f = (dsqrt(e)**3)
        ylm(i) = a*gx(i)*(5.d0*d-e)/f
     end do
  else if(is == 12) then
     a = dsqrt(21.d0/(32*PAI))
     do i = ni,kgp
        d = gz(i)**2
        e = gx(i)**2+gy(i)**2+d
        f = (dsqrt(e))**3
        ylm(i) = a*gy(i)*(5.d0*d-e)/f
     end do
  else if(is == 13) then
     a = dsqrt(105.d0/(16*PAI))
     do i = ni,kgp
        b = gx(i)**2
        c = gy(i)**2
        e = b+c+gz(i)**2
        f = dsqrt(e)**3
        ylm(i) = a*gz(i)*(b-c)/f
     end do
  else if(is == 14) then
     a = dsqrt(105.d0/PAI4)
     do i = ni,kgp
        f = (dsqrt(gx(i)**2+gy(i)**2+gz(i)**2))**3
        ylm(i) = a*gx(i)*gy(i)*gz(i)/f
     end do
  else if(is == 15) then
     a = dsqrt(35.d0/(32*PAI))
     do i = ni,kgp
        b = gx(i)**2
        c = gy(i)**2
        f = (dsqrt(b+c+gz(i)**2))**3
        ylm(i) = a*gx(i)*(b-3.d0*c)/f
     end do
  else if(is == 16) then
     a = dsqrt(35.d0/(32*PAI))
     do i = ni,kgp
        b = gx(i)**2
        c = gy(i)**2
        f = (dsqrt(b+c+gz(i)**2))**3
        ylm(i) = a*gy(i)*(3.d0*b-c)/f
     end do
  else if(is == 17) then
     a = 3.d0/8.d0/dsqrt(PAI4)
     do i = ni,kgp
        d = gz(i)**2
        e = gx(i)**2+gy(i)**2+d
        d = d/e
        ylm(i) = a*(5.d0*d*(7.d0*d-6.d0  )  +3.d0)
     end do
  else if(is == 18) then
     a = 15.d0/4.d0/dsqrt(10*PAI)
     do i = ni,kgp
        d = gz(i)**2
        e = gx(i)**2+gy(i)**2+d
        d = d/e
        ylm(i) = a*gz(i)*gx(i)*(7.d0*d-3.d0  )/e
     end do
  else if(is == 19) then
     a = 15.d0/4.d0/dsqrt(10*PAI)
     do i = ni,kgp
        d = gz(i)**2
        e = gx(i)**2+gy(i)**2+d
        d = d/e
        ylm(i) = a*gy(i)*gz(i)*(7.d0*d-3.d0  )/e
     end do
  else if(is == 20) then
     a = 15.d0/8.d0/dsqrt(5*PAI)
     do i = ni,kgp
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b+c+d
        b = b/e
        c = c/e
        d = d/e
        ylm(i) = a*(7.d0*d-1.d0)*(b-c)
     end do
  else if(is == 21) then
     a = 15.d0/4.d0/dsqrt(5*PAI)
     do i = ni,kgp
        d = gz(i)**2
        e = gx(i)**2+gy(i)**2+d
        d = d/e
        ylm(i) = a*(7.d0*d-1.d0)*gx(i)*gy(i)/e
     end do
  else if(is == 22) then
     a = 105.d0/4.d0/dsqrt(70*PAI)
     do i = ni,kgp
        b = gx(i)**2
        c = gy(i)**2
        e = (b+c+gz(i)**2)
        b = b/e
        c = c/e
        ylm(i) = a*(b-3.d0*c)*gz(i)*gx(i)/e
     end do
  else if(is == 23) then
     a = 105.d0/4.d0/dsqrt(70*PAI)
     do i = ni,kgp
        b = gx(i)**2
        c = gy(i)**2
        e = (b+c+gz(i)**2)
        b = b/e
        c = c/e
        ylm(i) = a*(3.d0*b-c)*gy(i)*gz(i)/e
     end do
  else if(is == 24) then
     a = 105.d0/16.d0/dsqrt(35*PAI)
     do i = ni,kgp
        b = gx(i)**2
        c = gy(i)**2
        e = (b+c+gz(i)**2)
        b = b/e
        c = c/e
        ylm(i) = a*((b-c)**2-4.d0*b*c)
     end do
  else if(is == 25) then
     a = 105.d0/4.d0/dsqrt(35*PAI)
     do i = ni,kgp
        b = gx(i)**2
        c = gy(i)**2
        e = (b+c+gz(i)**2)
        b = b/e
        c = c/e
        ylm(i) = a*(b-c)*gx(i)*gy(i)/e
     end do
  end if
end subroutine sphr

! ========================================
subroutine sphr_general( kgp, is, gx, gy, gz, ylm )
  use m_Const_Parameters, only  : DP, PAI, PAI4
  implicit none
  integer, intent(in)          :: kgp, is
  real(kind=DP), dimension(kgp), intent(in)  :: gx, gy, gz
  real(kind=DP), dimension(kgp), intent(out) :: ylm

  real(kind=DP) :: g_minimum = 1.d-20
  integer       :: i,ni
  real(kind=DP) :: a,b,c,d,e,f

  Do i=1, kgp
     if (gx(i)**2+gy(i)**2+gz(i)**2 < g_minimum) then
        ylm(i) = dsqrt(1.d0/PAI4)
        cycle
     endif
     select case (is)
! ****************************** L = 0 ************************
     case (1)
        a = dsqrt(1.d0/PAI4)
        ylm(i) = a
! ****************************** L = 1 ************************
     case (2)
        a = dsqrt(3.d0/PAI4)
        ylm(i) = a*gx(i)/(dsqrt(gx(i)**2+gy(i)**2+gz(i)**2))
     case (3)
        a = dsqrt(3.d0/PAI4)
        ylm(i) = a*gy(i)/(dsqrt(gx(i)**2+gy(i)**2+gz(i)**2))
     case (4)
        a = dsqrt(3.d0/PAI4)
        ylm(i) = a*gz(i)/(dsqrt(gx(i)**2+gy(i)**2+gz(i)**2))
! ****************************** L = 2 ************************
     case (5)
        a = dsqrt(5.d0/(16*PAI))
        d = gz(i)**2
        e = gx(i)**2+gy(i)**2+d
        ylm(i) = a*(3*d-e)/e
     case (6)
        a = dsqrt(15.d0/(16*PAI))
        b = gx(i)**2
        c = gy(i)**2
        e = b+c+gz(i)**2
        ylm(i) = a*(b-c)/e
     case (7)
        a = dsqrt(15.d0/PAI4)
        e = gx(i)**2+gy(i)**2+gz(i)**2
        ylm(i) = a*gx(i)*gy(i)/e
     case (8)
        a = dsqrt(15.d0/PAI4)
        e = gx(i)**2+gy(i)**2+gz(i)**2
        ylm(i) = a*gy(i)*gz(i)/e
     case (9)
        a = dsqrt(15.d0/PAI4)
        e = gx(i)**2+gy(i)**2+gz(i)**2
        ylm(i) = a*gz(i)*gx(i)/e
! ****************************** L = 3 ************************
     case (10)
        a = dsqrt(7.d0/(16*PAI))
        d = gz(i)**2
        e = gx(i)**2+gy(i)**2+d
        f = (dsqrt(e))**3
        ylm(i) = a*gz(i)*(5.d0*d-3.d0*e)/f
     case (11)
        a = dsqrt(21.d0/(32*PAI))
        d = gz(i)**2
        e = gx(i)**2+gy(i)**2+d
        f = (dsqrt(e)**3)
        ylm(i) = a*gx(i)*(5.d0*d-e)/f
     case (12)
        a = dsqrt(21.d0/(32*PAI))
        d = gz(i)**2
        e = gx(i)**2+gy(i)**2+d
        f = (dsqrt(e))**3
        ylm(i) = a*gy(i)*(5.d0*d-e)/f
     case (13)
        a = dsqrt(105.d0/(16*PAI))
        b = gx(i)**2
        c = gy(i)**2
        e = b+c+gz(i)**2
        f = dsqrt(e)**3
        ylm(i) = a*gz(i)*(b-c)/f
     case (14)
        a = dsqrt(105.d0/PAI4)
        f = (dsqrt(gx(i)**2+gy(i)**2+gz(i)**2))**3
        ylm(i) = a*gx(i)*gy(i)*gz(i)/f
     case (15)
        a = dsqrt(35.d0/(32*PAI))
        b = gx(i)**2
        c = gy(i)**2
        f = (dsqrt(b+c+gz(i)**2))**3
        ylm(i) = a*gx(i)*(b-3.d0*c)/f
     case (16)
        a = dsqrt(35.d0/(32*PAI))
        b = gx(i)**2
        c = gy(i)**2
        f = (dsqrt(b+c+gz(i)**2))**3
        ylm(i) = a*gy(i)*(3.d0*b-c)/f
! ****************************** L = 4 ************************
     case (17)
        a = 3.d0/8.d0/dsqrt(PAI4)
        d = gz(i)**2
        e = gx(i)**2+gy(i)**2+d
        d = d/e
        ylm(i) = a*(5.d0*d*(7.d0*d-6.d0  )  +3.d0)
     case (18)
        a = 15.d0/4.d0/dsqrt(10*PAI)
        d = gz(i)**2
        e = gx(i)**2+gy(i)**2+d
        d = d/e
        ylm(i) = a*gz(i)*gx(i)*(7.d0*d-3.d0  )/e
     case (19)
        a = 15.d0/4.d0/dsqrt(10*PAI)
        d = gz(i)**2
        e = gx(i)**2+gy(i)**2+d
        d = d/e
        ylm(i) = a*gy(i)*gz(i)*(7.d0*d-3.d0  )/e
     case (20)
        a = 15.d0/8.d0/dsqrt(5*PAI)
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b+c+d
        b = b/e
        c = c/e
        d = d/e
        ylm(i) = a*(7.d0*d-1.d0)*(b-c)
     case (21)
        a = 15.d0/4.d0/dsqrt(5*PAI)
        d = gz(i)**2
        e = gx(i)**2+gy(i)**2+d
        d = d/e
        ylm(i) = a*(7.d0*d-1.d0)*gx(i)*gy(i)/e
     case (22)
        a = 105.d0/4.d0/dsqrt(70*PAI)
        b = gx(i)**2
        c = gy(i)**2
        e = (b+c+gz(i)**2)
        b = b/e
        c = c/e
        ylm(i) = a*(b-3.d0*c)*gz(i)*gx(i)/e
     case (23)
        a = 105.d0/4.d0/dsqrt(70*PAI)
        b = gx(i)**2
        c = gy(i)**2
        e = (b+c+gz(i)**2)
        b = b/e
        c = c/e
        ylm(i) = a*(3.d0*b-c)*gy(i)*gz(i)/e
     case (24)
        a = 105.d0/16.d0/dsqrt(35*PAI)
        b = gx(i)**2
        c = gy(i)**2
        e = (b+c+gz(i)**2)
        b = b/e
        c = c/e
        ylm(i) = a*((b-c)**2-4.d0*b*c)
     case (25)
        a = 105.d0/4.d0/dsqrt(35*PAI)
        b = gx(i)**2
        c = gy(i)**2
        e = (b+c+gz(i)**2)
        b = b/e
        c = c/e
        ylm(i) = a*(b-c)*gx(i)*gy(i)/e
! ****************************** L = 5 ************************
     case (26)
        a = sqrt( 11.0d0 /256.d0 /PAI )
        b = 63.d0 *gz(i)**4
        c = 70.d0 *gz(i)**2
        d = gx(i)**2 +gy(i)**2 +gz(i)**2
        e = d *d
        f = d **(5./2.)
        ylm(i) = a *gz(i) *( b -c *d +15.d0 *e ) /f
     case (27)
        a = sqrt( 165.d0 /256.d0 /PAI )
        b = 21.d0 *gz(i)**4
        c = 14.d0 *gz(i)**2
        d = gx(i)**2 +gy(i)**2 +gz(i)**2
        e = d *d
        f = d **(5./2.)
        ylm(i) = a *gx(i) *( b -c *d +e ) /f
     case (28)
        a = sqrt( 165.d0 /256.d0 /PAI )
        b = 21.d0 *gz(i)**4
        c = 14.d0 *gz(i)**2
        d = gx(i)**2 +gy(i)**2 +gz(i)**2
        e = d *d
        f = d **(5./2.)
        ylm(i) = a *gy(i) *( b -c *d +e ) /f
     case (29)
        a = sqrt( 1155.d0 /64.d0 /PAI )
        b = 2.d0 *gz(i)**2 -gx(i)**2 -gy(i)**2
        c = gx(i)**2 -gy(i)**2
        d = gx(i)**2 +gy(i)**2 +gz(i)**2
        f = d **(5./2.)
        ylm(i) = a *b *c *gz(i) /f
     case (30)
        a = sqrt( 1155.d0 /64.d0 /PAI ) *2.d0
        b = 2.d0 *gz(i)**2 -gx(i)**2 -gy(i)**2
        d = gx(i)**2 +gy(i)**2 +gz(i)**2
        f = d **(5./2.)
        ylm(i) = a *b *gx(i) *gy(i) *gz(i) /f
     case (31)
        a = sqrt( 385.d0 /512.d0 /PAI )
        b = 8.d0 *gz(i)**2 -gx(i)**2 -gy(i)**2
        c = gx(i)**2 -3.d0 *gy(i)**2
        d = gx(i)**2 +gy(i)**2 +gz(i)**2
        f = d **(5./2.)
        ylm(i) = a *b *c *gx(i) /f
     case (32)
        a = sqrt( 385.d0 /512.d0 /PAI )
        b = 8.d0 *gz(i)**2 -gx(i)**2 -gy(i)**2
        c = 3.d0 *gx(i)**2 -gy(i)**2
        d = gx(i)**2 +gy(i)**2 +gz(i)**2
        f = d **(5./2.)
        ylm(i) = a *b *c *gy(i) /f
     case (33)
        a = sqrt( 3465.d0 /256.d0 /PAI )
        b = gx(i)**4 -6.d0 *gx(i)**2 *gy(i)**2 +gy(i)**4
        d = gx(i)**2 +gy(i)**2 +gz(i)**2
        f = d **(5./2.)
        ylm(i) = a *b *gz(i) /f
     case (34)
        a = sqrt( 3465.d0 /256.d0 /PAI ) *4.d0
        b = gx(i)**2 -gy(i)**2
        d = gx(i)**2 +gy(i)**2 +gz(i)**2
        f = d **(5./2.)
        ylm(i) = a *b *gx(i) *gy(i) *gz(i) /f
     case (35)
        a = sqrt( 693.d0 /512.d0 /PAI )
        b = gx(i)**4 -10.d0 *gx(i)**2 *gy(i)**2 +5.d0 *gy(i)**4
        d = gx(i)**2 +gy(i)**2 +gz(i)**2
        f = d **(5./2.)
        ylm(i) = a *b *gx(i) /f
     case (36)
        a = sqrt( 693.d0 /512.d0 /PAI )
        b = 5.d0 *gx(i)**4 -10.d0 *gx(i)**2 *gy(i)**2 +gy(i)**4
        d = gx(i)**2 +gy(i)**2 +gz(i)**2
        f = d **(5./2.)
        ylm(i) = a *b *gy(i) /f
! ****************************** L = 6 ************************
     case (37)
        a = sqrt( 13.d0 /1024.d0 /PAI )
        b = gx(i)**2 +gy(i)**2 +gz(i)**2
        c = b *b
        d = b *c
        f = 231.d0 *gz(i)**6 -315.d0 *gz(i)**4 *b +105.d0 *gz(i)**2 *c -5.d0 *d
        ylm(i) = a *f /d
     case (38)
        a = sqrt( 273.d0 /256.d0 /PAI )
        b = gx(i)**2 +gy(i)**2 +gz(i)**2
        c = b *b
        d = b *c
        f = 33.d0 *gz(i)**4 -30.d0 *gz(i)**2 *b +5.d0 *c
        ylm(i) = a *f *gx(i) *gz(i) /d
     case (39)
        a = sqrt( 273.d0 /256.d0 /PAI )
        b = gx(i)**2 +gy(i)**2 +gz(i)**2
        c = b *b
        d = b *c
        f = 33.d0 *gz(i)**4 -30.d0 *gz(i)**2 *b +5.d0 *c
        ylm(i) = a *f *gy(i) *gz(i) /d
     case (40)
        a = sqrt( 1365.d0 /2048.d0 /PAI )
        b = gx(i)**2 +gy(i)**2 +gz(i)**2
        c = b *b
        d = b *c
        e = gx(i)**2 -gy(i)**2
        f = 33.d0 *gz(i)**4 -18.d0 *gz(i)**2 *b +c
        ylm(i) = a *e *f /d
     case (41)
        a = sqrt( 1365.d0 /2048.d0 /PAI ) *2.d0
        b = gx(i)**2 +gy(i)**2 +gz(i)**2
        c = b *b
        d = b *c
        f = 33.d0 *gz(i)**4 -18.d0 *gz(i)**2 *b +c
        ylm(i) = a *f *gx(i) *gy(i) /d
     case (42)
        a = sqrt( 1365.d0 /512.d0 /PAI )
        b = gx(i)**2 +gy(i)**2 +gz(i)**2
        c = b *b
        d = b *c
        e = 11.d0 *gz(i)**2 -3.d0 *b
        f = gx(i)**2 -3.d0 *gy(i)**2
        ylm(i) = a *e *f *gx(i) *gz(i) /d
     case (43)
        a = sqrt( 1365.d0 /512.d0 /PAI )
        b = gx(i)**2 +gy(i)**2 +gz(i)**2
        c = b *b
        d = b *c
        e = 11.d0 *gz(i)**2 -3.d0 *b
        f = 3.d0 *gx(i)**2 -gy(i)**2
        ylm(i) = a *e *f *gy(i) *gz(i) /d
     case (44)
        a = sqrt( 819.d0 /1024.d0 /PAI )
        b = gx(i)**2 +gy(i)**2 +gz(i)**2
        c = b *b
        d = b *c
        e = 11.d0 *gz(i)**2 -b
        f = gx(i)**4 -6.d0 *gx(i)**2 *gy(i)**2 +gy(i)**4
        ylm(i) = a *e *f /d
     case (45)
        a = sqrt( 819.d0 /1024.d0 /PAI ) *4.d0
        b = gx(i)**2 +gy(i)**2 +gz(i)**2
        c = b *b
        d = b *c
        e = 11.d0 *gz(i)**2 -b
        f = gx(i)**2 -gy(i)**2
        ylm(i) = a *e *f *gx(i) *gy(i) /d
     case (46)
        a = sqrt( 9009.d0 /512.d0 /PAI )
        b = gx(i)**2 +gy(i)**2 +gz(i)**2
        c = b *b
        d = b *c
        f = gx(i)**4 -10.d0 *gx(i)**2 *gy(i)**2 +5.d0 *gy(i)**4
        ylm(i) = a *f *gx(i) *gz(i) /d
     case (47)
        a = sqrt( 9009.d0 /512.d0 /PAI )
        b = gx(i)**2 +gy(i)**2 +gz(i)**2
        c = b *b
        d = b *c
        f = 5.d0 *gx(i)**4 -10.d0 *gx(i)**2 *gy(i)**2 +gy(i)**4
        ylm(i) = a *f *gy(i) *gz(i) /d
     case (48)
        a = sqrt( 3003.d0 /2048.d0 /PAI )
        b = gx(i)**2 +gy(i)**2 +gz(i)**2
        c = b *b
        d = b *c
        f = gx(i)**6 -15.d0 *gx(i)**4 *gy(i)**2 +15.d0 *gx(i)**2 *gy(i)**4 -gy(i)**6
        ylm(i) = a *f /d
     case (49)
        a = sqrt( 3003.d0 /2048.d0 /PAI )
        b = gx(i)**2 +gy(i)**2 +gz(i)**2
        c = b *b
        d = b *c
        f = 6.d0 *gx(i)**4 -20.d0 *gx(i)**2 *gy(i)**2 +6.d0 *gy(i)**4
        ylm(i) = a * f* gx(i) *gy(i) /d
     end select
  end Do
end subroutine sphr_general
! ========================================

!!$C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
!!$C derivative of spherical harmonics.
!!$c   coded by H. Sawada 8th May 1997
!!$C*************************************************************
subroutine sphr_diff(kg1, kgp, is, gx,gy,gz,ylmd)
  use m_Const_Parameters, only  : DP, PAI, PAI4
  implicit none
  integer, intent(in)          :: kg1, kgp, is
  real(kind=DP), dimension(kg1),   intent(in)  :: gx,gy,gz
  real(kind=DP), dimension(kg1,3), intent(out) :: ylmd

  real(kind=DP) :: g_minimum = 1.d-20
  integer       :: i,ni
  real(kind=DP) :: a,b,c,d,e,e2,e3,gr,gr3,gr5
  if(gx(1)**2+gy(1)**2+gz(1)**2 < g_minimum) then
     ni = 2
     ylmd(1,1:3) = 0.d0
  else
     ni = 1
  end if
  if(is == 1) then
      a = 0.d0
      do i=ni,kgp
      ylmd(i,1) = a
      ylmd(i,2) = a
      ylmd(i,3) = a
      end do
  else if(is == 2) then
      a = dsqrt(3.d0/(4.d0*PAI))
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      gr = sqrt( e )
      gr3 = gr * e
      ylmd(i,1) = a* ( -b/gr3 + 1.d0/gr )
      ylmd(i,2) = a* ( -gx(i)*gy(i)/gr3 )
      ylmd(i,3) = a* ( -gx(i)*gz(i)/gr3 )
      end do
  else if(is == 3) then
      a = dsqrt(3.d0/(4.d0*PAI))
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      gr = sqrt( e )
      gr3 = gr * e
      ylmd(i,1) = a* ( -gx(i)*gy(i)/gr3 )
      ylmd(i,2) = a* ( -c/gr3 + 1.d0/gr )
      ylmd(i,3) = a* ( -gy(i)*gz(i)/gr3 )
      end do
  else if(is == 4) then
      a = dsqrt(3.d0/(4.d0*PAI))
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      gr = sqrt( e )
      gr3 = gr * e
      ylmd(i,1) = a* ( -gx(i)*gz(i)/gr3 )
      ylmd(i,2) = a* ( -gy(i)*gz(i)/gr3 )
      ylmd(i,3) = a* ( -d/gr3 + 1.d0/gr )
      end do
  else if(is == 5) then
      a = dsqrt(5.d0/(16.d0*PAI))
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      e2 = e**2
      ylmd(i,1) = a* ( -2.d0*gx(i)/e - 2.d0*gx(i)*(3.d0*d-e)/e2 )
      ylmd(i,2) = a* ( -2.d0*gy(i)/e - 2.d0*gy(i)*(3.d0*d-e)/e2 )
      ylmd(i,3) = a* (  4.d0*gz(i)/e - 2.d0*gz(i)*(3.d0*d-e)/e2 )
      end do
  else if(is == 6) then
      a = dsqrt(15.d0/(16.d0*PAI))
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      e2 = e**2
      ylmd(i,1) = a* ( -2.d0*gx(i)*(b-c)/e2 + 2.d0*gx(i)/e )
      ylmd(i,2) = a* ( -2.d0*gy(i)*(b-c)/e2 - 2.d0*gy(i)/e )
      ylmd(i,3) = a* ( -2.d0*gz(i)*(b-c)/e2 )
      end do
  else if(is == 7) then
      a = dsqrt(15.d0/(4.d0*PAI))
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      e2 = e**2
      ylmd(i,1) = a* ( -2.d0*b*gy(i)/e2 + gy(i)/e )
      ylmd(i,2) = a* ( -2.d0*c*gx(i)/e2 + gx(i)/e )
      ylmd(i,3) = a* ( -2.d0*gx(i)*gy(i)*gz(i)/e2 )
      end do
  else if(is == 8) then
      a = dsqrt(15.d0/(4.d0*PAI))
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      e2 = e**2
      ylmd(i,1) = a* ( -2.d0*gx(i)*gy(i)*gz(i)/e2 )
      ylmd(i,2) = a* ( -2.d0*c*gz(i)/e2 + gz(i)/e )
      ylmd(i,3) = a* ( -2.d0*d*gy(i)/e2 + gy(i)/e )
      end do
  else if(is == 9) then
      a = dsqrt(15.d0/(4.d0*PAI))
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      e2 = e**2
      ylmd(i,1) = a* ( -2.d0*b*gz(i)/e2 + gz(i)/e )
      ylmd(i,2) = a* ( -2.d0*gx(i)*gy(i)*gz(i)/e2 )
      ylmd(i,3) = a* ( -2.d0*d*gx(i)/e2 + gx(i)/e )
      end do
  else if(is == 10) then
      a = dsqrt(7.d0/(16.d0*PAI))
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      gr = sqrt( e )
      gr3 = gr * e
      gr5 = gr3 * e
      ylmd(i,1) = a* ( -6.d0*gx(i)*gz(i)/gr3 &
     &                 - 3.d0*gx(i)*gz(i)*(5.d0*d-3.d0*e)/gr5 )
      ylmd(i,2) = a* ( -6.d0*gy(i)*gz(i)/gr3 &
     &                 - 3.d0*gy(i)*gz(i)*(5.d0*d-3.d0*e)/gr5 )
      ylmd(i,3) = a* ( 4.d0*d/gr3 &
     &                 - 3.d0*d*(5.d0*d-3.d0*e)/gr5 &
     &                 + (5.d0*d-3.d0*e)/gr3 )
      end do
  else if(is == 11) then
      a = dsqrt(21.d0/(32.d0*PAI))
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      gr = sqrt( e )
      gr3 = gr * e
      gr5 = gr3 * e
      ylmd(i,1) = a* ( -2.d0*b/gr3 &
     &                 - 3.d0*b*(5.d0*d-e)/gr5 &
     &                 + (5.d0*d-e)/gr3 )
      ylmd(i,2) = a* ( -2.d0*gx(i)*gy(i)/gr3 &
     &                 - 3.d0*gx(i)*gy(i)*(5.d0*d-e)/gr5 )
      ylmd(i,3) = a* ( 8.d0*gx(i)*gz(i)/gr3 &
     &                 - 3.d0*gx(i)*gz(i)*(5.d0*d-e)/gr5 )
      end do
  else if(is == 12) then
      a = dsqrt(21.d0/(32.d0*PAI))
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      gr = sqrt( e )
      gr3 = gr * e
      gr5 = gr3 * e
      ylmd(i,1) = a* ( -2.d0*gx(i)*gy(i)/gr3 &
     &                 - 3.d0*gx(i)*gy(i)*(5.d0*d-e)/gr5 )
      ylmd(i,2) = a* ( -2.d0*c/gr3 &
     &                 - 3.d0*c*(5.d0*d-e)/gr5 &
     &                 + (5.d0*d-e)/gr3 )
      ylmd(i,3) = a* ( 8.d0*gy(i)*gz(i)/gr3 &
     &                 - 3.d0*gy(i)*gz(i)*(5.d0*d-e)/gr5 )
      end do
  else if(is == 13) then
      a = dsqrt(105.d0/(16.d0*PAI))
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      gr = sqrt( e )
      gr3 = gr * e
      gr5 = gr3 * e
      ylmd(i,1) = a* ( -3.d0*gx(i)*gz(i)*(b-c)/gr5 &
     &                 + 2.d0*gx(i)*gz(i)/gr3 )
      ylmd(i,2) = a* ( -3.d0*gy(i)*gz(i)*(b-c)/gr5 &
     &                 - 2.d0*gy(i)*gz(i)/gr3 )
      ylmd(i,3) = a* ( -3.d0*d*(b-c)/gr5 &
     &                 + (b-c)/gr3 )
      end do
  else if(is == 14) then
      a = dsqrt(105.d0/(4.d0*PAI))
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      gr = sqrt( e )
      gr3 = gr * e
      gr5 = gr3 * e
      ylmd(i,1) = a* ( -3.d0*b*gy(i)*gz(i)/gr5 &
     &                 + gy(i)*gz(i)/gr3 )
      ylmd(i,2) = a* ( -3.d0*gx(i)*c*gz(i)/gr5 &
     &                 + gx(i)*gz(i)/gr3 )
      ylmd(i,3) = a* ( -3.d0*gx(i)*gy(i)*d/gr5 &
     &                 + gx(i)*gy(i)/gr3 )
      end do
  else if(is == 15) then
      a = dsqrt(35.d0/(32.d0*PAI))
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      gr = sqrt( e )
      gr3 = gr * e
      gr5 = gr3 * e
      ylmd(i,1) = a* ( -3.d0*b*(b-3.d0*c)/gr5 &
     &                 + 3.d0*(b-c)/gr3 )
      ylmd(i,2) = a* ( -3.d0*gx(i)*gy(i)*(b-3.d0*c)/gr5 &
     &                 - 6.d0*gx(i)*gy(i)/gr3 )
      ylmd(i,3) = a* ( -3.d0*gx(i)*gz(i)*(b-3.d0*c)/gr5 )
      end do
  else if(is == 16) then
      a = dsqrt(35.d0/(32.d0*PAI))
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      gr = sqrt( e )
      gr3 = gr * e
      gr5 = gr3 * e
      ylmd(i,1) = a* ( -3.d0*gx(i)*gy(i)*(3.d0*b-c)/gr5 &
     &                 + 6.d0*gx(i)*gy(i)/gr3 )
      ylmd(i,2) = a* ( -3.d0*c*(3.d0*b-c)/gr5 &
     &                 + 3.d0*(b-c)/gr3 )
      ylmd(i,3) = a* ( -3.d0*gy(i)*gz(i)*(3.d0*b-c)/gr5 )
      end do
  else if(is == 17) then
      a = 3.d0/8.d0/dsqrt(4.d0*PAI)
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      e2 = e**2
      e3 = e * e2
      ylmd(i,1) = a* ( -60.d0*gx(i)*d/e2 &
     &                 - 20.d0*gx(i)*d*(7.d0*d-6.d0*e)/e3 )
      ylmd(i,2) = a* ( -60.d0*gy(i)*d/e2 &
     &                 - 20.d0*gy(i)*d*(7.d0*d-6.d0*e)/e3 )
      ylmd(i,3) = a* ( 10.d0*gz(i)*d/e2 &
     &                 - 20.d0*gz(i)*d*(7.d0*d-6.d0*e)/e3 &
     &                 + 10.d0*gz(i)*(7.d0*d-6.d0*e)/e2 )
      end do
  else if(is == 18) then
      a = 15.d0/4.d0/dsqrt(10.d0*PAI)
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      e2 = e**2
      e3 = e * e2
      ylmd(i,1) = a* ( -6.d0*b*gz(i)/e2 &
     &                 - 4.d0*b*gz(i)*(7.d0*d-3.d0*e)/e3 &
     &                 + gz(i)*(7.d0*d-3.d0*e)/e2 )
      ylmd(i,2) = a* ( -6.d0*gx(i)*gy(i)*gz(i)/e2 &
     &                 - 4.d0*gx(i)*gy(i)*gz(i)*(7.d0*d-3.d0*e)/e3 )
      ylmd(i,3) = a* ( 8.d0*gx(i)*d/e2 &
     &                 - 4.d0*gx(i)*d*(7.d0*d-3.d0*e)/e3 &
     &                 + gx(i)*(7.d0*d-3.d0*e)/e2 )
      end do
  else if(is == 19) then
      a = 15.d0/4.d0/dsqrt(10.d0*PAI)
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      e2 = e**2
      e3 = e * e2
      ylmd(i,1) = a* ( -6.d0*gx(i)*gy(i)*gz(i)/e2 &
     &                 - 4.d0*gx(i)*gy(i)*gz(i)*(7.d0*d-3.d0*e)/e3 )
      ylmd(i,2) = a* ( -6.d0*c*gz(i)/e2 &
     &                 - 4.d0*c*gz(i)*(7.d0*d-3.d0*e)/e3 &
     &                 + gz(i)*(7.d0*d-3.d0*e)/e2 )
      ylmd(i,3) = a* ( 8.d0*gy(i)*d/e2 &
     &                 - 4.d0*gy(i)*d*(7.d0*d-3.d0*e)/e3 &
     &                 + gy(i)*(7.d0*d-3.d0*e)/e2 )
      end do
  else if(is == 20) then
      a = 15.d0/8.d0/dsqrt(5.d0*PAI)
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      e2 = e**2
      e3 = e * e2
      ylmd(i,1) = a* ( -2.d0*gx(i)*(b-c)/e2 &
     &                 - 4.d0*gx(i)*(b-c)*(7.d0*d-e)/e3 &
     &                 + 2.d0*gx(i)*(7.d0*d-e)/e2 )
      ylmd(i,2) = a* ( -2.d0*gy(i)*(b-c)/e2 &
     &                 - 4.d0*gy(i)*(b-c)*(7.d0*d-e)/e3 &
     &                 - 2.d0*gy(i)*(7.d0*d-e)/e2 )
      ylmd(i,3) = a* ( 12.d0*gz(i)*(b-c)/e2 &
     &                 - 4.d0*gz(i)*(b-c)*(7.d0*d-e)/e3 )
      end do
  else if(is == 21) then
      a = 15.d0/4.d0/dsqrt(5.d0*PAI)
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      e2 = e**2
      e3 = e * e2
      ylmd(i,1) = a* ( -2.d0*b*gy(i)/e2 &
     &                 - 4.d0*b*gy(i)*(7.d0*d-e)/e3 &
     &                 + gy(i)*(7.d0*d-e)/e2 )
      ylmd(i,2) = a* ( -2.d0*gx(i)*c/e2 &
     &                 - 4.d0*gx(i)*c*(7.d0*d-e)/e3 &
     &                 + gx(i)*(7.d0*d-e)/e2 )
      ylmd(i,3) = a* ( 12.d0*gx(i)*gy(i)*gz(i)/e2 &
     &                 - 4.d0*gx(i)*gy(i)*gz(i)*(7.d0*d-e)/e3 )
      end do
  else if(is == 22) then
      a = 105.d0/4.d0/dsqrt(70.d0*PAI)
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      e2 = e**2
      e3 = e * e2
      ylmd(i,1) = a* ( -4.d0*b*gz(i)*(b-3.d0*c)/e3 &
     &                 + 3.d0*gz(i)*(b-c)/e2 )
      ylmd(i,2) = a* ( -4.d0*gx(i)*gy(i)*gz(i)*(b-3.d0*c)/e3 &
     &                 - 6.d0*gx(i)*gy(i)*gz(i)/e2 )
      ylmd(i,3) = a* ( -4.d0*gx(i)*d*(b-3.d0*c)/e3 &
     &                 + gx(i)*(b-3.d0*c)/e2 )
      end do
  else if(is == 23) then
      a = 105.d0/4.d0/dsqrt(70.d0*PAI)
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      e2 = e**2
      e3 = e * e2
      ylmd(i,1) = a* ( -4.d0*gx(i)*gy(i)*gz(i)*(3.d0*b-c)/e3 &
     &                 + 6.d0*gx(i)*gy(i)*gz(i)/e2 )
      ylmd(i,2) = a* ( -4.d0*c*gz(i)*(3.d0*b-c)/e3 &
     &                 + 3.d0*gz(i)*(b-c)/e2 )
      ylmd(i,3) = a* ( -4.d0*gy(i)*d*(3.d0*b-c)/e3 &
     &                 + gy(i)*(3.d0*b-c)/e2 )
      end do
  else if(is == 24) then
      a = 105.d0/16.d0/dsqrt(35.d0*PAI)
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      e2 = e**2
      e3 = e * e2
      ylmd(i,1) = a* ( -4.d0*gx(i)*(-4.d0*b*c+(b-c)**2)/e3 &
     &                 + (-8.d0*gx(i)*c+4.d0*gx(i)*(b-c))/e2 )
      ylmd(i,2) = a* ( -4.d0*gy(i)*(-4.d0*b*c+(b-c)**2)/e3 &
     &                 + (-8.d0*b*gy(i)-4.d0*gy(i)*(b-c))/e2 )
      ylmd(i,3) = a* ( -4.d0*gz(i)*(-4.d0*b*c+(b-c)**2)/e3 )
      end do
  else if(is == 25) then
      a = 105.d0/4.d0/dsqrt(35.d0*PAI)
      do i=ni,kgp
      b = gx(i)**2
      c = gy(i)**2
      d = gz(i)**2
      e = b + c + d
      e2 = e**2
      e3 = e * e2
      ylmd(i,1) = a* ( -4.d0*b*gy(i)*(b-c)/e3 &
     &                 + gy(i)*(3.d0*b-c)/e2 )
      ylmd(i,2) = a* ( -4.d0*gx(i)*c*(b-c)/e3 &
     &                 + gx(i)*(b-3.d0*c)/e2 )
      ylmd(i,3) = a* ( -4.d0*gx(i)*gy(i)*gz(i)*(b-c)/e3 )
      end do
  end if

end subroutine sphr_diff

subroutine sphr_diff_general(kg1, kgp, is, gx,gy,gz,ylmd)
  use m_Const_Parameters, only  : DP, PAI, PAI4
  implicit none
  integer, intent(in)          :: kg1, kgp, is
  real(kind=DP), dimension(kg1),   intent(in)  :: gx,gy,gz
  real(kind=DP), dimension(kg1,3), intent(out) :: ylmd

  real(kind=DP) :: g_minimum = 1.d-20
  integer       :: i,ni
  real(kind=DP) :: a,b,c,d,e,e2,e3,gr,gr3,gr5

  Do i=1, kgp
     if (gx(i)**2+gy(i)**2+gz(i)**2 < g_minimum) then
        ylmd(i,1:3) = 0.0d0
        cycle
     endif
     select case (is)
     case (1)
        a = 0.d0
        ylmd(i,1:3) = a
     case (2)
        a = dsqrt(3.d0/(4.d0*PAI))
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        gr = sqrt( e )
        gr3 = gr * e
        ylmd(i,1) = a* ( -b/gr3 + 1.d0/gr )
        ylmd(i,2) = a* ( -gx(i)*gy(i)/gr3 )
        ylmd(i,3) = a* ( -gx(i)*gz(i)/gr3 )
     case (3)
        a = dsqrt(3.d0/(4.d0*PAI))
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        gr = sqrt( e )
        gr3 = gr * e
        ylmd(i,1) = a* ( -gx(i)*gy(i)/gr3 )
        ylmd(i,2) = a* ( -c/gr3 + 1.d0/gr )
        ylmd(i,3) = a* ( -gy(i)*gz(i)/gr3 )
     case (4)
        a = dsqrt(3.d0/(4.d0*PAI))
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        gr = sqrt( e )
        gr3 = gr * e
        ylmd(i,1) = a* ( -gx(i)*gz(i)/gr3 )
        ylmd(i,2) = a* ( -gy(i)*gz(i)/gr3 )
        ylmd(i,3) = a* ( -d/gr3 + 1.d0/gr )
     case (5)
        a = dsqrt(5.d0/(16.d0*PAI))
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        e2 = e**2
        ylmd(i,1) = a* ( -2.d0*gx(i)/e - 2.d0*gx(i)*(3.d0*d-e)/e2 )
        ylmd(i,2) = a* ( -2.d0*gy(i)/e - 2.d0*gy(i)*(3.d0*d-e)/e2 )
        ylmd(i,3) = a* (  4.d0*gz(i)/e - 2.d0*gz(i)*(3.d0*d-e)/e2 )
     case (6)
        a = dsqrt(15.d0/(16.d0*PAI))
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        e2 = e**2
        ylmd(i,1) = a* ( -2.d0*gx(i)*(b-c)/e2 + 2.d0*gx(i)/e )
        ylmd(i,2) = a* ( -2.d0*gy(i)*(b-c)/e2 - 2.d0*gy(i)/e )
        ylmd(i,3) = a* ( -2.d0*gz(i)*(b-c)/e2 )

     case (7)
        a = dsqrt(15.d0/(4.d0*PAI))
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        e2 = e**2
        ylmd(i,1) = a* ( -2.d0*b*gy(i)/e2 + gy(i)/e )
        ylmd(i,2) = a* ( -2.d0*c*gx(i)/e2 + gx(i)/e )
        ylmd(i,3) = a* ( -2.d0*gx(i)*gy(i)*gz(i)/e2 )
     case (8)
        a = dsqrt(15.d0/(4.d0*PAI))
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        e2 = e**2
        ylmd(i,1) = a* ( -2.d0*gx(i)*gy(i)*gz(i)/e2 )
        ylmd(i,2) = a* ( -2.d0*c*gz(i)/e2 + gz(i)/e )
        ylmd(i,3) = a* ( -2.d0*d*gy(i)/e2 + gy(i)/e )

     case (9)
        a = dsqrt(15.d0/(4.d0*PAI))
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        e2 = e**2
        ylmd(i,1) = a* ( -2.d0*b*gz(i)/e2 + gz(i)/e )
        ylmd(i,2) = a* ( -2.d0*gx(i)*gy(i)*gz(i)/e2 )
        ylmd(i,3) = a* ( -2.d0*d*gx(i)/e2 + gx(i)/e )
     case (10)
        a = dsqrt(7.d0/(16.d0*PAI))
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        gr = sqrt( e )
        gr3 = gr * e
        gr5 = gr3 * e
        ylmd(i,1) = a* ( -6.d0*gx(i)*gz(i)/gr3 &
             &                 - 3.d0*gx(i)*gz(i)*(5.d0*d-3.d0*e)/gr5 )
        ylmd(i,2) = a* ( -6.d0*gy(i)*gz(i)/gr3 &
             &                 - 3.d0*gy(i)*gz(i)*(5.d0*d-3.d0*e)/gr5 )
        ylmd(i,3) = a* ( 4.d0*d/gr3 &
             &                 - 3.d0*d*(5.d0*d-3.d0*e)/gr5 &
             &                 + (5.d0*d-3.d0*e)/gr3 )

     case (11)
        a = dsqrt(21.d0/(32.d0*PAI))
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        gr = sqrt( e )
        gr3 = gr * e
        gr5 = gr3 * e
        ylmd(i,1) = a* ( -2.d0*b/gr3 &
             &                 - 3.d0*b*(5.d0*d-e)/gr5 &
             &                 + (5.d0*d-e)/gr3 )
        ylmd(i,2) = a* ( -2.d0*gx(i)*gy(i)/gr3 &
             &                 - 3.d0*gx(i)*gy(i)*(5.d0*d-e)/gr5 )
        ylmd(i,3) = a* ( 8.d0*gx(i)*gz(i)/gr3 &
             &                 - 3.d0*gx(i)*gz(i)*(5.d0*d-e)/gr5 )
     case (12)
        a = dsqrt(21.d0/(32.d0*PAI))
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        gr = sqrt( e )
        gr3 = gr * e
        gr5 = gr3 * e
        ylmd(i,1) = a* ( -2.d0*gx(i)*gy(i)/gr3 &
             &                 - 3.d0*gx(i)*gy(i)*(5.d0*d-e)/gr5 )
        ylmd(i,2) = a* ( -2.d0*c/gr3 &
             &                 - 3.d0*c*(5.d0*d-e)/gr5 &
             &                 + (5.d0*d-e)/gr3 )
        ylmd(i,3) = a* ( 8.d0*gy(i)*gz(i)/gr3 &
             &                 - 3.d0*gy(i)*gz(i)*(5.d0*d-e)/gr5 )
     case (13)
        a = dsqrt(105.d0/(16.d0*PAI))
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        gr = sqrt( e )
        gr3 = gr * e
        gr5 = gr3 * e
        ylmd(i,1) = a* ( -3.d0*gx(i)*gz(i)*(b-c)/gr5 &
             &                 + 2.d0*gx(i)*gz(i)/gr3 )
        ylmd(i,2) = a* ( -3.d0*gy(i)*gz(i)*(b-c)/gr5 &
             &                 - 2.d0*gy(i)*gz(i)/gr3 )
        ylmd(i,3) = a* ( -3.d0*d*(b-c)/gr5 &
     &                 + (b-c)/gr3 )
     case (14)
        a = dsqrt(105.d0/(4.d0*PAI))
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        gr = sqrt( e )
        gr3 = gr * e
        gr5 = gr3 * e
        ylmd(i,1) = a* ( -3.d0*b*gy(i)*gz(i)/gr5 &
             &                 + gy(i)*gz(i)/gr3 )
        ylmd(i,2) = a* ( -3.d0*gx(i)*c*gz(i)/gr5 &
             &                 + gx(i)*gz(i)/gr3 )
        ylmd(i,3) = a* ( -3.d0*gx(i)*gy(i)*d/gr5 &
             &                 + gx(i)*gy(i)/gr3 )
     case (15)
        a = dsqrt(35.d0/(32.d0*PAI))
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        gr = sqrt( e )
        gr3 = gr * e
        gr5 = gr3 * e
        ylmd(i,1) = a* ( -3.d0*b*(b-3.d0*c)/gr5 &
             &                 + 3.d0*(b-c)/gr3 )
        ylmd(i,2) = a* ( -3.d0*gx(i)*gy(i)*(b-3.d0*c)/gr5 &
             &                 - 6.d0*gx(i)*gy(i)/gr3 )
        ylmd(i,3) = a* ( -3.d0*gx(i)*gz(i)*(b-3.d0*c)/gr5 )
     case (16)
        a = dsqrt(35.d0/(32.d0*PAI))
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        gr = sqrt( e )
        gr3 = gr * e
        gr5 = gr3 * e
        ylmd(i,1) = a* ( -3.d0*gx(i)*gy(i)*(3.d0*b-c)/gr5 &
             &                 + 6.d0*gx(i)*gy(i)/gr3 )
        ylmd(i,2) = a* ( -3.d0*c*(3.d0*b-c)/gr5 &
             &                 + 3.d0*(b-c)/gr3 )
        ylmd(i,3) = a* ( -3.d0*gy(i)*gz(i)*(3.d0*b-c)/gr5 )
     
     case (17)
        a = 3.d0/8.d0/dsqrt(4.d0*PAI)
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        e2 = e**2
        e3 = e * e2
        ylmd(i,1) = a* ( -60.d0*gx(i)*d/e2 &
             &                 - 20.d0*gx(i)*d*(7.d0*d-6.d0*e)/e3 )
        ylmd(i,2) = a* ( -60.d0*gy(i)*d/e2 &
             &                 - 20.d0*gy(i)*d*(7.d0*d-6.d0*e)/e3 )
        ylmd(i,3) = a* ( 10.d0*gz(i)*d/e2 &
             &                 - 20.d0*gz(i)*d*(7.d0*d-6.d0*e)/e3 &
             &                 + 10.d0*gz(i)*(7.d0*d-6.d0*e)/e2 )
     case (18)
        a = 15.d0/4.d0/dsqrt(10.d0*PAI)
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        e2 = e**2
        e3 = e * e2
        ylmd(i,1) = a* ( -6.d0*b*gz(i)/e2 &
             &                 - 4.d0*b*gz(i)*(7.d0*d-3.d0*e)/e3 &
             &                 + gz(i)*(7.d0*d-3.d0*e)/e2 )
        ylmd(i,2) = a* ( -6.d0*gx(i)*gy(i)*gz(i)/e2 &
             &                 - 4.d0*gx(i)*gy(i)*gz(i)*(7.d0*d-3.d0*e)/e3 )
        ylmd(i,3) = a* ( 8.d0*gx(i)*d/e2 &
             &                 - 4.d0*gx(i)*d*(7.d0*d-3.d0*e)/e3 &
             &                 + gx(i)*(7.d0*d-3.d0*e)/e2 )
     case (19)
        a = 15.d0/4.d0/dsqrt(10.d0*PAI)
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        e2 = e**2
        e3 = e * e2
        ylmd(i,1) = a* ( -6.d0*gx(i)*gy(i)*gz(i)/e2 &
             &                 - 4.d0*gx(i)*gy(i)*gz(i)*(7.d0*d-3.d0*e)/e3 )
        ylmd(i,2) = a* ( -6.d0*c*gz(i)/e2 &
             &                 - 4.d0*c*gz(i)*(7.d0*d-3.d0*e)/e3 &
             &                 + gz(i)*(7.d0*d-3.d0*e)/e2 )
        ylmd(i,3) = a* ( 8.d0*gy(i)*d/e2 &
             &                 - 4.d0*gy(i)*d*(7.d0*d-3.d0*e)/e3 &
             &                 + gy(i)*(7.d0*d-3.d0*e)/e2 )
     case (20)
        a = 15.d0/8.d0/dsqrt(5.d0*PAI)
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        e2 = e**2
        e3 = e * e2
        ylmd(i,1) = a* ( -2.d0*gx(i)*(b-c)/e2 &
             &                 - 4.d0*gx(i)*(b-c)*(7.d0*d-e)/e3 &
             &                 + 2.d0*gx(i)*(7.d0*d-e)/e2 )
        ylmd(i,2) = a* ( -2.d0*gy(i)*(b-c)/e2 &
             &                 - 4.d0*gy(i)*(b-c)*(7.d0*d-e)/e3 &
             &                 - 2.d0*gy(i)*(7.d0*d-e)/e2 )
        ylmd(i,3) = a* ( 12.d0*gz(i)*(b-c)/e2 &
             &                 - 4.d0*gz(i)*(b-c)*(7.d0*d-e)/e3 )
     case (21)
        a = 15.d0/4.d0/dsqrt(5.d0*PAI)
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        e2 = e**2
        e3 = e * e2
        ylmd(i,1) = a* ( -2.d0*b*gy(i)/e2 &
             &                 - 4.d0*b*gy(i)*(7.d0*d-e)/e3 &
             &                 + gy(i)*(7.d0*d-e)/e2 )
        ylmd(i,2) = a* ( -2.d0*gx(i)*c/e2 &
             &                 - 4.d0*gx(i)*c*(7.d0*d-e)/e3 &
             &                 + gx(i)*(7.d0*d-e)/e2 )
        ylmd(i,3) = a* ( 12.d0*gx(i)*gy(i)*gz(i)/e2 &
     &                 - 4.d0*gx(i)*gy(i)*gz(i)*(7.d0*d-e)/e3 )

     case (22)
        a = 105.d0/4.d0/dsqrt(70.d0*PAI)
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        e2 = e**2
        e3 = e * e2
        ylmd(i,1) = a* ( -4.d0*b*gz(i)*(b-3.d0*c)/e3 &
             &                 + 3.d0*gz(i)*(b-c)/e2 )
        ylmd(i,2) = a* ( -4.d0*gx(i)*gy(i)*gz(i)*(b-3.d0*c)/e3 &
             &                 - 6.d0*gx(i)*gy(i)*gz(i)/e2 )
        ylmd(i,3) = a* ( -4.d0*gx(i)*d*(b-3.d0*c)/e3 &
             &                 + gx(i)*(b-3.d0*c)/e2 )
     case (23)
        a = 105.d0/4.d0/dsqrt(70.d0*PAI)
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        e2 = e**2
        e3 = e * e2
        ylmd(i,1) = a* ( -4.d0*gx(i)*gy(i)*gz(i)*(3.d0*b-c)/e3 &
             &                 + 6.d0*gx(i)*gy(i)*gz(i)/e2 )
        ylmd(i,2) = a* ( -4.d0*c*gz(i)*(3.d0*b-c)/e3 &
             &                 + 3.d0*gz(i)*(b-c)/e2 )
        ylmd(i,3) = a* ( -4.d0*gy(i)*d*(3.d0*b-c)/e3 &
     &                 + gy(i)*(3.d0*b-c)/e2 )
        
     case (24)
        a = 105.d0/16.d0/dsqrt(35.d0*PAI)
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        e2 = e**2
        e3 = e * e2
        ylmd(i,1) = a* ( -4.d0*gx(i)*(-4.d0*b*c+(b-c)**2)/e3 &
             &                 + (-8.d0*gx(i)*c+4.d0*gx(i)*(b-c))/e2 )
        ylmd(i,2) = a* ( -4.d0*gy(i)*(-4.d0*b*c+(b-c)**2)/e3 &
             &                 + (-8.d0*b*gy(i)-4.d0*gy(i)*(b-c))/e2 )
        ylmd(i,3) = a* ( -4.d0*gz(i)*(-4.d0*b*c+(b-c)**2)/e3 )

     case (25)
        a = 105.d0/4.d0/dsqrt(35.d0*PAI)
        b = gx(i)**2
        c = gy(i)**2
        d = gz(i)**2
        e = b + c + d
        e2 = e**2
        e3 = e * e2
        ylmd(i,1) = a* ( -4.d0*b*gy(i)*(b-c)/e3 &
             &                 + gy(i)*(3.d0*b-c)/e2 )
        ylmd(i,2) = a* ( -4.d0*gx(i)*c*(b-c)/e3 &
             &                 + gx(i)*(b-3.d0*c)/e2 )
        ylmd(i,3) = a* ( -4.d0*gx(i)*gy(i)*gz(i)*(b-c)/e3 )
     end select
  end Do

end subroutine sphr_diff_general

subroutine eqivvl(value)
!!$                           @(#)bottom_Subroutines.F90 1.3 00/07/11 14:12:19
!!$                by T.Yamasaki (JRCAT,ATP)
  use m_Const_Parameters, only : DP
  implicit none
  real(kind=DP), intent(inout) :: value
  value = value
end subroutine eqivvl

subroutine skip_commentlines(nf,ipri,nfout,ierr)
  use m_ErrorMessages, only : NO_ERROR, EOF_REACHED
  integer, intent(in) :: nf,ipri,nfout
  integer, intent(out) :: ierr
  integer, parameter  ::    len_str = 80
  character(len=len_str) :: str
  logical :: comment_statement
  comment_statement = .true.
  ierr = 0
  do while(comment_statement)
     read(nf,'(a80)',end=1000) str
     if(str(1:1) == '#'.or. str(1:1) == '$' .or. str(1:1) == '!' &
          & .or. str(1:1) == '%' .or. str(1:1) == '*') then
        if(ipri >= 2) write(nfout,'(a80)') str
     else if(len(trim(str)) == 0) then
        if(ipri >= 2) write(nfout,'(a80)') str
     else
        comment_statement = .false.
     endif
  enddo
  backspace(nf)
  ierr = NO_ERROR
  return
1000 ierr = EOF_REACHED
!!$1000 write(nfout,'(" eof is reached")')
!!$  stop ' eof is reached in file reading'
  return
end subroutine skip_commentlines

!=====================================================================
function fn_number_of_words(string)
!=====================================================================
!
!  Calculates number of words in a string
!
!  M. Okamoto
!
!---------------------------------------------------------------------
  implicit none
  integer :: fn_number_of_words
  character(*),intent(in) :: string
  integer        :: ipos, lpos
  character(256) :: tmp_string
  tmp_string = string
  lpos       = len(tmp_string)
  fn_number_of_words = 0
  if (lpos == 0) then
     go to 99
  end if
  tmp_string = adjustl(tmp_string)
  if (index(tmp_string,' ') == 1) then
     go to 99
  end if
  do while (index(tmp_string,' ') > 1)
     fn_number_of_words = fn_number_of_words + 1
     ipos               = index(tmp_string,' ')
     if (ipos == 0) then
        go to 99
     end if
     tmp_string = adjustl(tmp_string(ipos:lpos))
     tmp_string = adjustl(tmp_string)
  end do
99 continue
end function fn_number_of_words

character(len=256) function to_string(i,keta)
  integer, intent(in) :: i
  integer, intent(in) :: keta
  character :: string*256
  character :: ch*256
  integer :: iketa
  integer :: itmp
  iketa = keta
  itmp = -1
  mype_conf = i
  if(mype_conf>0) itmp=int(log10(real(mype_conf)))+1
  if(itmp>keta)iketa=itmp
  write(ch,*) iketa
  write(string,'(i'//trim(adjustl(ch))//'.'//trim(adjustl(ch))//')') i
  to_string = trim(adjustl(string))
end function to_string

!===============================================================================

function get_target_index(l,xmat,xin) result(res)
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in) :: l
  real(kind=DP), intent(in), dimension(l) :: xmat
  real(kind=DP), intent(in) :: xin
  integer :: i,res
  do i=1,l-1
    if((xin .ge. xmat(i)) .and. (xin .le. xmat(i+1)))then
      res = i
      return
    endif
  enddo
  res = 0
  return
end function get_target_index

function trilinear_interpolation(l,m,n,cmat,xmat,ymat,zmat,xin,yin,zin) result(res)
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in) :: l,m,n
  real(kind=DP), intent(in), dimension(l,m,n) :: cmat
  real(kind=DP), intent(in), dimension(l) :: xmat
  real(kind=DP), intent(in), dimension(m) :: ymat
  real(kind=DP), intent(in), dimension(n) :: zmat
  real(kind=DP), intent(in) :: xin, yin, zin
  integer :: i,ltarget,mtarget,ntarget
  integer :: get_target_index
  real(kind=DP) :: xd,x0,x1,yd,y0,y1,zd,z0,z1
  real(kind=DP) :: c000,c100,c001,c101,c010,c110,c011,c111
  real(kind=DP) :: c00,c01,c10,c11
  real(kind=DP) :: c0,c1
  real(kind=DP) :: res
  if(l.le.1 .or. m.le.1 .or. n.le.1)then
    call phase_error_with_msg(6, 'error from trilinear interpolation : l, m, n must be larget than 2',__LINE__,__FILE__)
  endif
  ltarget=get_target_index(l,xmat,xin)
  mtarget=get_target_index(m,ymat,yin)
  ntarget=get_target_index(n,zmat,zin)
  if(ltarget.lt.1 .or. mtarget.lt.1 .or. ntarget.lt.1)then
    call phase_error_with_msg(6, 'error from trilinear interpolation : input coordinate out of range',__LINE__,__FILE__)
  endif
  x0 = xmat(ltarget);x1 = xmat(ltarget+1);xd = (xin-x0)/(x1-x0)
  y0 = ymat(mtarget);y1 = ymat(mtarget+1);yd = (yin-y0)/(y1-y0)
  z0 = zmat(ntarget);z1 = zmat(ntarget+1);zd = (zin-z0)/(z1-z0)
  c000 = cmat(ltarget,mtarget,ntarget)
  c100 = cmat(ltarget+1,mtarget,ntarget)
  c001 = cmat(ltarget,mtarget,ntarget+1)
  c101 = cmat(ltarget+1,mtarget,ntarget+1)
  c010 = cmat(ltarget,mtarget+1,ntarget)
  c110 = cmat(ltarget+1,mtarget+1,ntarget)
  c011 = cmat(ltarget,mtarget+1,ntarget+1)
  c111 = cmat(ltarget+1,mtarget+1,ntarget+1)
  c00 = c000*(1.d0-xd)+c100*xd
  c01 = c001*(1.d0-xd)+c101*xd
  c10 = c010*(1.d0-xd)+c110*xd
  c11 = c011*(1.d0-xd)+c111*xd
  c0 = c00*(1.d0-yd)+c10*yd
  c1 = c01*(1.d0-yd)+c11*yd
  res = c0*(1.d0-zd)+c1*zd
  return
end function trilinear_interpolation

integer function get_nspher_cylm( il, im )       ! complex Ylm
  implicit none
  integer, intent(in) :: il, im
  !
  integer :: num
  num = il**2 +(im+il) +1;
  get_nspher_cylm = num;
  return
end function get_nspher_cylm

! get a pair of normal random number by the Box-Muller method
subroutine normal_random_number(z1, z2)
  use m_Const_Parameters, only : PAI, DP
  implicit none
  real(kind=DP), intent(out) :: z1, z2
  real(kind=DP) :: X, Y, twologx, twopiy
  call random_number(X)
  call random_number(Y)
  twologx = max(-2.d0 * log(X),0.d0)
  twopiy  = 2*PAI*Y
  z1 = sqrt(twologx) * cos(twopiy)
  z2 = sqrt(twologx) * sin(twopiy)

  return

end subroutine normal_random_number

