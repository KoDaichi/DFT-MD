!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  SUBROUINE:  width2, wd_fermi_error1, wd_fermi_error2
!
!  AUTHOR(S): N. Hamada, T. Yamasaki   August/20/2003
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
subroutine wd_fermi_error1(nfout,emin,emax,tot,totch)
! $Id: b_Fermi.F90 633 2020-12-01 05:11:03Z jkoga $
  use m_Const_Parameters, only : DP
  implicit none
  integer,       intent(in) :: nfout
  real(kind=DP), intent(in) :: emin, emax, tot, totch
  write(nfout,'(" emin = ",f10.4,"  emax = ",f10.4)') emin, emax
  write(nfout,'("  tot = ",f10.4," totch = ",f10.4)') tot, totch
  !stop ' === stop in sub.fermi_parabolic (too few of states) =='
  call phase_error_with_msg(nfout, ' === stop in sub.fermi_parabolic (too few of states) ==',__LINE__,__FILE__)
end subroutine wd_fermi_error1

subroutine wd_fermi_error2(nfout,e1,e2,efermi,emin,emax,tot,totch,neg,MAXITR)
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in) :: nfout, neg, MAXITR
  real(kind=DP), intent(in) :: e1, e2,efermi, emin,emax, tot, totch
  
  write(nfout,'("     e1  = ", f10.5)') e1
  write(nfout,'("     e2  = ", f10.5)') e2
  write(nfout,'(" efermi  = ", f10.5)') efermi
  write(nfout,'("   emin  = ", f10.5, "  emax  = ",f10.5)') emin,emax
  write(nfout,'("    tot  = ", f10.5, " totch  = ",f10.5)') tot,totch
  write(nfout,'("  occup ( i = 1, ",i5,")")') neg
  write(nfout,'(" jcount exceeds MAXITR (= ",i8," (in FERMI)")') MAXITR
  !stop ' stop at <wd_fermi_error2>'
  call phase_error_with_msg(nfout, ' <wd_fermi_error2>',__LINE__,__FILE__)
end subroutine wd_fermi_error2

!!$C SUBROUTINE WIDTH2
!!$c                           @(#)width2.f 9.1 97/05/08 14:49:38 
!!$C
!!$C          1983/05/18 :   NORIAKI HAMADA
!!$C
!!$C----*----1----*----2----*----3----*----4----*----5----*----6----*----7
!!$C
subroutine width2( E,EFF,W, DOS,OCC  )
  use m_Const_Parameters, only : DP
  implicit none
  real(kind=DP), intent(in)  :: E, EFF, W
  real(kind=DP), intent(out) :: DOS, OCC
!!$C          EIGENENERGY IS BROADENED
!!$C
!!$C          E :  EIGENENERGY
!!$C          EF:  FERMI ENERGY
!!$C          W :  4*W IS BOTTOM LENGTH
!!$C          DOS    :  DENSITY OF STATES AT EF
!!$C          OCCUP :  OCUPATION FRACTION OF ELECTRON
  real(kind=DP) :: ee, ww

  if(w <= 0.0d0) STOP ' === STOP AT SUB.WIDTHE. (W<=0.0D0) ==='
#ifdef __TIMER_SUB__
!  call timer_sta(705)
#endif
  ee=(eff-e)/w
  ww=4.0d0*w
  if(ee <= -2.0d0) then
     dos    = 0.0d0
     occ  = 0.0d0
  else if(ee < -1.0d0) then
     dos    = (ee+2.0d0)**2/ww
     occ  = ((ee+2.0)**3)/(12.0d0)
  else if(ee < 1.0d0) then
     dos    = (2.0d0-ee**2)/ww
     occ  = (6.0d0+6.0d0*ee-(ee**3))/12.0d0
  else if(ee < 2.0d0) then
     dos    = ((ee-2.0d0)**2)/ww
     occ  = (12.0d0+(ee-2.0d0)**3)/12.0d0
  else
     dos    = 0.0d0
     occ  = 1.0d0
  end if
#ifdef __TIMER_SUB__
!  call timer_end(705)
#endif
end subroutine width2

! ============= KT_add =========================== 13.0E
subroutine width_fermi_dirac( ene, ef, width, dos, occ )
  use m_Const_Parameters, only : DP
  implicit none

  real(kind=DP), intent(in) :: ene, ef, width
  real(kind=DP), intent(out) :: dos, occ

  real(kind=DP) :: ee, c1, c2

  if ( width < 1.0D-10 ) call phase_error_with_msg(6,"smearing width is too small",__LINE__,__FILE__)

  ee = ( ene -ef )/width
  if(ee.gt.100)then
     occ = 0.d0
     dos = 0.d0
  else if (ee.lt.-100)then
     occ = 1.d0
     dos = 0.d0
  else
  c1 = exp( ee );   c2 = cosh( ee/2.0D0 )
!
  occ = 1.0D0 / ( 1.0D0 +c1 )
  dos = 0.25D0 / c2**2 /width
  endif

end subroutine width_fermi_dirac
! ============================================== 13.0E

subroutine width_methfessel_paxton(n,ene,ef,width,dos,occ,entropy)
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in) :: n
  real(kind=DP), intent(in) :: ene, ef, width
  real(kind=DP), intent(out) :: dos, occ
  real(kind=DP), intent(out) :: entropy
  real(kind=DP) :: ee
  if ( width < 1.0D-10 ) call phase_error_with_msg(6,"smearing width is too small",__LINE__,__FILE__)

  ee = ( ene -ef )/width
  if(ee.gt.7)then
     occ = 0.d0
     dos = 0.d0
     entropy = 0.d0
  else if (ee.lt.-7)then
     occ = 1.d0
     dos = 0.d0
     entropy = 0.0d0
  else
     occ = step_function_mp(n,ee)
     dos = delta_function_mp(n,ee)
     entropy = 0.5d0*coeffA(n)*hermite_polynomial(2*n,ee)*exp(-ee*ee)
  endif

  return

  contains 

  function delta_function_mp(n,x) result(res)
    use m_Const_Parameters, only : DP
    implicit none
    integer, intent(in) :: n
    real(kind=DP), intent(in) :: x
    real(kind=DP) :: res,ex
    integer :: i
    res = 0.d0
    ex = exp(-x*x)
    do i=0,n
       res = res+coeffA(i)*hermite_polynomial(2*i,x)*ex
    enddo
  end function delta_function_mp

  function step_function_mp(n,x) result(res)
    use m_Const_Parameters, only : DP
    implicit none
    integer, intent(in) :: n
    real(kind=DP), intent(in) :: x
    real(kind=DP) :: res,an,h2n,ex
    integer :: i
    res = 0.5d0*erfc(x)
    if (n==0) return
    ex = exp(-x*x)
    do i=1,n
       an  = coeffA(i)
       h2n = hermite_polynomial(2*i-1,x)
       res = res + an*h2n*ex
    enddo
    return
  end function step_function_mp

  recursive function hermite_polynomial(n,x) result(res)
    use m_Const_Parameters, only : DP
    implicit none
    integer, intent(in) :: n
    real(kind=DP), intent(in) :: x
    real(kind=DP) :: res
    real(kind=DP) :: y1,y2
    real(kind=DP) :: x2,x3,x4,x5
    integer :: i
    if (n==0) then
      res = 1.d0
      return
    endif
    if (n==1) then
      res = 2.d0*x
      return
    endif
    if (n==2) then
      res = 4.d0*x*x-2.d0
      return
    endif
    if (n==3) then
      res = 8.d0*x*x*x-12.d0*x
      return
    endif
    if (n==4) then
      x2 = x*x
      res = x2*x2*16.d0-48.d0*x2+12.d0
      return
    endif
    if (n==5) then
      x3=x*x*x
      res = 32.d0*x3*x*x-160.d0*x3+120.d0*x
      return
    endif
    if (n==6) then
      x2 = x*x
      x3 = x2*x
      res = 64.d0*x3*x3-480.d0*x3*x+720.d0*x2-120.d0
      return
    endif
    if (n==7) then
      x3=x*x*x
      x5=x3*x*x
      res = 128.d0*x5*x*x-1344.d0*x5+3360.d0*x3-1680.d0*x
      return
    endif
    if (n==8) then
      x2 = x*x
      x3 = x2*x
      x4 = x3*x
      res = 256.d0*x4*x4-3584.d0*x3*x3+13440.d0*x4-13440.d0*x2+1680.d0
      return
    endif
    if (n==9) then
      x3 = x*x*x
      x4 = x3*x
      x5 = x4*x
      res = 512.d0*x5*x4-9216.d0*x3*x4+48384.d0*x5-80640.d0*x3+30240.d0*x
      return
    endif
    if (n==10) then
      x2 = x*x
      x3 = x2*x
      x4 = x3*x
      x5 = x4*x
      res = 1024.d0*x5*x5-23040.d0*x5*x3+161280*x3*x3-403200.d0*x4+302400.d0*x2-30240.d0
      return
    endif
    y1 = hermite_polynomial(n-1,x)
    y2 = hermite_polynomial(n-2,x)
    res  = 2.d0*x*y1-2.d0*(n-1)*y2
    return
  end function hermite_polynomial

  function coeffA(n) result(res)
    use m_Const_Parameters, only : DP,PAI
    implicit none
    integer, intent(in) :: n
    real(kind=DP) :: res
    integer :: i
    real(kind=DP) :: kaijo
    kaijo = 1
    do i=2,n
       kaijo = kaijo*i
    enddo
    res = (-1)**n/(kaijo*4**n*dsqrt(PAI))
    return
  end function coeffA

end subroutine width_methfessel_paxton

