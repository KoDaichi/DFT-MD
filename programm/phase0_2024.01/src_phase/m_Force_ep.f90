!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  MODULE:  m_Force
!
!  AUTHOR(S): T. Uchiyama, T. Yamasaki   August/20/2003
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
!
module m_Force
! $Id: m_Force_ep.f90 570 2017-04-21 20:34:50Z yamasaki $
  use m_Const_Parameters,     only : DP, RELAX, BONDLENGTH_FIX, HEAT_BATH &
       &                           , WITHOUTTAG, WITHTAG
  use m_Timing,               only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters,   only : forccr,printable
  use m_Crystal_Structure,    only : altv,univol
  use m_Ionic_System,         only : natm,natm2,iwei,imdtyp,pos

  implicit none

!-- For Si atom --
  real(kind=DP), parameter :: a_cut = 1.8d0,sigma = 3.959164919d0, eps = 7.968005097d-2
  real(kind=DP), parameter :: bmax = sigma*a_cut
!--
  real(kind=DP),allocatable,dimension(:,:) :: forc_l !d(natm,3)
  real(kind=DP)                            :: forcmx
  real(kind=DP)                            :: etotal

contains
  subroutine m_Force_alloc()
    allocate(forc_l(natm,3))
  end subroutine m_Force_alloc

  subroutine m_Force_stlweb()
! ** Stillinger-Weber potential **
!                          T. Uchiyama,  June 6, 1998
! ** Transfer from f77 to f90 **
!                      corded by k.iwata July 30,1998
! ** Revised by T. Yamasaki, Oct. 1998
    integer, parameter :: nbd = 24
    real(kind=DP), pointer, dimension(:,:) :: posb, forc   !d(natm2,3)
    integer,       pointer, dimension(:,:) :: ipar         !d(katm*nbd,3)
    real(kind=DP) :: r1(3),r2(3)
    integer       :: npar

    integer       :: id_sname = -1
    call tstatc0_begin('stlweb ',id_sname)

    allocate(posb(natm2,3));  allocate(forc(natm2,3)); forc = 0.d0
    allocate(ipar(natm*nbd,3))

    call full_ionic_coordinates    !-(contained here) pos -> posb
    call two_body_term(npar)       !-(contained here) -> forc
    call three_body_term(nbd,npar) !-(contained here) -> forc
    call copy_back_to_forc_l       !-(contained here) forc -> forc_l
    call get_forc_maximum          !-(contained here) forc_l -> forcmx

    deallocate(ipar)
    deallocate(forc); deallocate(posb)

    call tstatc0_end(id_sname)
  contains
    subroutine full_ionic_coordinates
      integer :: i,ia

      do i = 1, natm
         posb(i,1:3)= pos(i,1:3)
      end do

      ia= natm
      do i= 1, natm
         if(iwei(i) == 2) then
            ia= ia +1
            posb(ia,1:3)= -posb(i,1:3)
         end if
      end do

      do i = 1, natm2
         posb(i,1:3) = dmod( posb(i,1:3), 1.d0 )
      end do
    end subroutine full_ionic_coordinates

! ***** 2-body term ***********************************
    subroutine two_body_term(npar)
      integer, intent(out) :: npar
!     Interaction with self-images is not included,
!   where the cell size is assumed to be larger than
!   the cutoff radius.
!
      real(kind=DP),pointer,dimension(:,:) :: femp2   !d(3,2)
      integer       :: mm, m1, m2,lx,ly,lz, j
      real(kind=DP) :: b, eemp2, cp(3)

      allocate(femp2(3,2))

      etotal = 0.d0
      mm= 0
      do m1 = 1,    natm2-1
         do m2 = m1+1, natm2
            do lx= -1, 1
               do ly= -1, 1
                  do lz= -1, 1

                     cp(1) = posb(m2,1) - posb(m1,1) + lx
                     cp(2) = posb(m2,2) - posb(m1,2) + ly
                     cp(3) = posb(m2,3) - posb(m1,3) + lz

                     r1 = matmul(altv,cp)   ! pucv -> cartesian coordinates

                     b= dsqrt( pvec( r1, r1 ))

                     if(b < bmax) then
                        mm= mm +1
                        ipar(mm,1)= m1
                        ipar(mm,2)= m2
                        ipar(mm,3)= ( lx+1 )*9 +( ly+1 )*3 +lz+1
           !  ------------------------------------------------
                        call emp2(r1,eemp2,femp2)   !-(m_Force_ep)
                        etotal= etotal +eemp2
                        do j= 1, 3
                           forc(m1,j)= forc(m1,j) +femp2(j,1)
                           forc(m2,j)= forc(m2,j) +femp2(j,2)
                        end do
           !  ------------------------------------------------
                     end if

                  end do
               end do
            end do

         end do
      end do
      npar= mm
!!$      print *,' ! npar = ', npar
      deallocate(femp2)
    end subroutine two_body_term

! ***** 3-body term **********************************
    subroutine three_body_term(nbd,npar)
      integer, intent(in) :: nbd,npar
      real(kind=DP)       :: kb(nbd),t(3,nbd)
      real(kind=DP),pointer,dimension(:,:) :: femp3   ! d(3,3)
      integer             :: m, nb, k, m1, m2, lx,ly,lz,j1,j2,j
      real(kind=DP)       :: eemp3
      logical :: over
      allocate(femp3(3,3))

      over = .false.
      do m = 1, natm2

         nb= 0
         do k= 1, npar
            if(ipar(k,1) == m .or. ipar(k,2) == m) then
               nb= nb +1
               if(nb > nbd) then
                  if(printable) write(6,'(" !warning nbd, nb = ",i8,i8)') nbd,nb 
               end if
               m1= ipar(k,1)
               m2= ipar(k,2)

               lx= ipar(k,3)/9 -1
               ly= ( ipar(k,3) -lx*9 -9 )/3 -1
               lz= ipar(k,3) -lx*9 -ly*3 -13

               if(nb <= nbd) then

               t(1,nb)= posb(m2,1) -posb(m1,1) +lx
               t(2,nb)= posb(m2,2) -posb(m1,2) +ly
               t(3,nb)= posb(m2,3) -posb(m1,3) +lz

               if(m1 == m) kb(nb)= m2
               if(m2 == m) then
                  kb(nb)= m1 
                  t(1,nb)= -t(1,nb)
                  t(2,nb)= -t(2,nb)
                  t(3,nb)= -t(3,nb)
               end if

               end if
            end if
         end do

         if(nb > nbd) over = .true.

         do j1= 1,    nb-1
            do j2= j1+1, nb

               r1 = matmul(altv,t(1:3,j1))
               r2 = matmul(altv,t(1:3,j2))

         ! -------------------------------------------------
               call emp3(r1,r2,eemp3,femp3)
               etotal= etotal +eemp3

               m1 = kb(j1)
               m2 = kb(j2)
               do j= 1, 3
                  forc(m, j)= forc(m, j) +femp3(j,1)
                  forc(m1,j)= forc(m1,j) +femp3(j,2)
                  forc(m2,j)= forc(m2,j) +femp3(j,3)
               end do
         ! -------------------------------------------------

            end do
         end do
      end do

      if(over) stop ' stop at three_body_term'

!!$      print *,' ! etotal = ' , etotal
      deallocate(femp3)
    end subroutine three_body_term

    subroutine copy_back_to_forc_l
      integer :: i, j
      do i= 1, natm
         do j= 1,3
            forc_l(i,j)= forc(i,j)
         end do
      end do
    end subroutine copy_back_to_forc_l

    subroutine get_forc_maximum
      integer       :: i
      real(kind=DP) :: fca

      forcmx = 0.d0
      do i = 1, natm
         if(imdtyp(i) == RELAX .or. imdtyp(i) == BONDLENGTH_FIX &
              & .or. imdtyp(i) > HEAT_BATH) then
            fca = dsqrt(forc_l(i,1)**2+forc_l(i,2)**2+forc_l(i,3)**2)
            if(forcmx < fca) forcmx = fca
         end if
      end do
      write(6,'(" !!f forcmx = ", f8.4," <<m_Force_stlweb.get_forc_maximum>>")') forcmx
!!$      print *, ' !D forcmx = ', forcmx
    end subroutine get_forc_maximum

  end subroutine m_Force_stlweb

! ****************************
  function pvec( a, b )
! ****************************

    real(kind=DP) :: pvec
    real(kind=DP) :: a(3), b(3)

    pvec= a(1)*b(1) +a(2)*b(2) +a(3)*b(3)

  end function pvec

! **************************************************************
  subroutine emp2(r,eemp2,femp2)   
    !   Stillinger and Weber potential,    coded by T. Miyake

    real(kind=DP),intent(in), dimension(3)   :: r
    real(kind=DP),intent(out)                :: eemp2
    real(kind=DP),intent(out),dimension(3,2) :: femp2

    real(kind=DP) :: x,y,z, rad, femp2r
    x = r(1)
    y = r(2)
    z = r(3)
    rad = dsqrt(x**2 + y**2 + z**2)

    call vemp2SiSi(rad,eemp2,femp2r)

    femp2(1,1) = - femp2r * x / rad
    femp2(2,1) = - femp2r * y / rad
    femp2(3,1) = - femp2r * z / rad
    femp2(1,2) = - femp2(1,1)
    femp2(2,2) = - femp2(2,1)
    femp2(3,2) = - femp2(3,1)

  end subroutine emp2

! ***********************************************************************
  subroutine emp3(r12,r13,eemp3,femp3)
    real(kind=DP), intent(in), dimension(3)  :: r12,r13
    real(kind=DP), intent(out)               :: eemp3
    real(kind=DP), intent(out),dimension(3,3):: femp3

    real(kind=DP) :: x12,y12,z12,x13,y13,z13

    x12 = r12(1)
    y12 = r12(2)
    z12 = r12(3)
    x13 = r13(1)
    y13 = r13(2)
    z13 = r13(3)

    call vemp3SiSiSi(x12,y12,z12,x13,y13,z13, eemp3,femp3)
  end subroutine emp3

! ***********************************************************************
!  two body contributions to Si-Si interaction 
!  Stillinger and Weber, Phys.Rev.B31,5262(1985)
!  atomic unit
  subroutine vemp2SiSi(r,eemp2,femp2r)
    real(kind=DP),intent(in)  :: r
    real(kind=DP),intent(out) :: eemp2,femp2r

    real(kind=DP),parameter   :: A = 7.049556277d0, B = 0.6022245584d0
    integer, parameter        :: p = 4, q = 0

    real(kind=DP) :: r0, f2, df2

! =======================================================================
!  begin:

    r0 = r / sigma

    if (r0 >=  a_cut) then 

       eemp2 = 0d0
       femp2r = 0d0

    else 

       f2 = A * (B*r0**(-p) - r0**(-q)) * dexp(1d0/(r0-a_cut))
       df2 = A * (-p*B*r0**(-p-1) + q*r0**(-q-1)) * &
            &  dexp(1d0/(r0-a_cut)) - &
            &  A * (B*r0**(-p) - r0**(-q)) * &
            &  dexp(1d0/(r0-a_cut)) * (r0-a_cut)**(-2)

       eemp2 = eps * f2
       femp2r = - eps * df2 / sigma

    endif

  end subroutine vemp2SiSi

! ***********************************************************************
!  three body contributions to Si-Si-Si interaction 
!  Stillinger and Weber, Phys.Rev.B31,5262(1985)
!  atomic unit
  subroutine vemp3SiSiSi(x12,y12,z12,x13,y13,z13,eemp3,femp3)
    real(kind=DP), intent(inout) :: x12,y12,z12,x13,y13,z13
    real(kind=DP), intent(out)   :: eemp3,femp3(3,3)

    real(kind=DP) :: r12,r13
    real(kind=DP) :: h,dh(3,3)

    x12 = x12 / sigma
    y12 = y12 / sigma
    z12 = z12 / sigma
    x13 = x13 / sigma
    y13 = y13 / sigma
    z13 = z13 / sigma

    r12 = dsqrt(x12**2 + y12**2 + z12**2)
    r13 = dsqrt(x13**2 + y13**2 + z13**2)

    if ((r12 < a_cut).and.(r13 < a_cut)) then 

       call sw3body(a_cut,r12,r13,x12,x13,y12,y13,z12,z13,h,dh)
       eemp3 = eps * h
       femp3 = -eps/sigma * dh
    else

        eemp3 = 0.d0
        femp3 = 0.d0

    endif

  end subroutine vemp3SiSiSi

! ***********************************************************************
!  h and its derivatives in Stillinger-Weber potential
  subroutine sw3body(a_cut,r12,r13,x12,x13,y12,y13,z12,z13, &
       &  h,dh)
    implicit none

    real(kind=DP), parameter :: lamda=21.0d0
    real(kind=DP), parameter :: gamma=1.20d0

    real(kind=DP) ::  a_cut
    real(kind=DP) ::  r12,r13,x12,y12,z12,x13,y13,z13
    real(kind=DP) ::  h,dh(3,3)
    real(kind=DP) ::  b,c

! =======================================================================
!  begin:
    b = x12*x13 + y12*y13 + z12*z13
    c = b / (r12*r13) + 1d0 / 3d0
!
    if (c.eq.0d0) then 
       h = 0d0
       dh = 0.d0
    else
       h = lamda * dexp(gamma/(r12-a_cut) + gamma/(r13-a_cut)) & 
            &  * c**2

       dh(1,1) = h * (gamma*x12/r12/(r12-a_cut)**2 + &
            &  gamma*x13/r13/(r13-a_cut)**2) &
            &  + 2d0*h/c * (-(x12+x13)/(r12*r13) - b*(r12*(-x13/r13) &
            &   +r13*(-x12/r12)) / (r12*r13)**2)
       dh(2,1) = h * (gamma*y12/r12/(r12-a_cut)**2 + &
            &  gamma*y13/r13/(r13-a_cut)**2) &
            &  + 2d0*h/c * (-(y12+y13)/(r12*r13) - b*(r12*(-y13/r13) &
            &  +r13*(-y12/r12)) / (r12*r13)**2)
       dh(3,1) = h * (gamma*z12/r12/(r12-a_cut)**2 + &
            &  gamma*z13/r13/(r13-a_cut)**2) &
            &  + 2d0*h/c * (-(z12+z13)/(r12*r13) - b*(r12*(-z13/r13) &
            &  +r13*(-z12/r12)) / (r12*r13)**2)

       dh(1,2) = h * (gamma*(-x12/r12)/(r12-a_cut)**2) & 
            & + 2d0*h/c * (x13/(r12*r13) - b*(x12/r12)/(r12**2*r13))
       dh(2,2) = h * (gamma*(-y12/r12)/(r12-a_cut)**2) &
            & + 2d0*h/c * (y13/(r12*r13) - b*(y12/r12)/(r12**2*r13))
       dh(3,2) = h * (gamma*(-z12/r12)/(r12-a_cut)**2) &
            & + 2d0*h/c * (z13/(r12*r13) - b*(z12/r12)/(r12**2*r13))

       dh(1,3) = h * (gamma*(-x13/r13)/(r13-a_cut)**2) &
            & + 2d0*h/c * (x12/(r12*r13) - b*(x13/r13)/(r12*r13**2))
       dh(2,3) = h * (gamma*(-y13/r13)/(r13-a_cut)**2) &
            & + 2d0*h/c * (y12/(r12*r13) - b*(y13/r13)/(r12*r13**2))
       dh(3,3) = h * (gamma*(-z13/r13)/(r13-a_cut)**2) & 
            & + 2d0*h/c * (z12/(r12*r13) - b*(z13/r13)/(r12*r13**2))

    endif
!
!       if (ideb == 1) write(*,*)'sw3body in'
!
  end subroutine sw3body

  logical function m_Forces_are_Converged_core()
    if(forcmx < forccr) then
       m_Forces_are_Converged_core = .true.
    else
       m_Forces_are_Converged_core = .false.
    end if
!!$    write(6,'(" !!f forcmx, forccr = ",2f12.6)') forcmx, forccr
  end function m_Forces_are_Converged_core

  subroutine m_Force_wd_force_and_cps(nf, withorwithout, cps, ndim)
    integer, intent(in) :: nf, withorwithout,ndim
    real(kind=DP), intent(in) :: cps(ndim,3)
    real(kind=DP) :: fca
    integer :: ia

    if(withorwithout == WITHOUTTAG) then
       do ia = 1, natm
          fca = dsqrt(forc_l(ia,1)**2+forc_l(ia,2)**2+forc_l(ia,3)**2)
          write(nf,'(" ",i4,3f15.9,4f12.6)') ia &
               &, cps(ia,1),cps(ia,2),cps(ia,3) &
               &, forc_l(ia,1), forc_l(ia,2), forc_l(ia,3),fca
       end do
    else if(withorwithout == WITHTAG) then
       do ia = 1, natm
          fca = dsqrt(forc_l(ia,1)**2+forc_l(ia,2)**2+forc_l(ia,3)**2)
          write(nf,'(" !forc ",i4,3f12.6,4f12.6)') ia &
               &, cps(ia,1),cps(ia,2),cps(ia,3) &
               &, forc_l(ia,1), forc_l(ia,2), forc_l(ia,3), fca
       end do
    else
       write(nf,'(" keyword in the argument list of m_Force_wd_force_and_cps is invalid")')
    end if

  end subroutine m_Force_wd_force_and_cps

  subroutine m_Force_wd_force_cps_cpd(nf, withorwithout, cps, cpd,ndim)
    integer, intent(in) :: nf, withorwithout,ndim
    real(kind=DP), intent(in) :: cps(ndim,3), cpd(ndim,3)
    integer :: ia

    if(withorwithout == WITHOUTTAG) then
       do ia = 1, natm
          write(nf,'(" ",i4,3f15.9,3f12.6, 3f12.6)') ia &
               &, cps(ia,1:3), cpd(ia,1:3), forc_l(ia,1:3)
       end do
    else if(withorwithout == WITHTAG) then
       do ia = 1, natm
          write(nf,'(" !forc ",i4,3f12.6,3f12.6,3f12.6)') ia &
               &, cps(ia,1:3),cpd(ia,1:3),forc_l(ia,1:3)
       end do
    else
       write(nf,'(" keyword in the argument list of m_Force_wd_force_cps_cpd is invalid")')
    end if

  end subroutine m_Force_wd_force_cps_cpd

end module m_Force
