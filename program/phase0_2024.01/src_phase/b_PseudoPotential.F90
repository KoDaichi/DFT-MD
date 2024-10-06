!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: matrix_inversion, psbhs0, qitgft_diff, qitgft, 
!             read_chgpc_nrc_mord, read_nrc_mord, read_tau_eps_nrc_mord, 
!             read_itau_ivanl, coef_simpson_integration, rmeshs, ggabp,
!             ggacorl_pbe, ggaexch_pbe, ggacorl_pw91, ggaexch_pw91,
!             gga_xc_potential_atomic, get_gradient_of_rho
!
!  AUTHOR(S): T. Yamasaki   June/25/2003
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
! $Id: b_PseudoPotential.F90 570 2017-04-21 20:34:50Z yamasaki $
subroutine get_gradient_of_rho(jsm,mddiff,kmesh,mesh,a &
     &, rad,h,fdiff,coeff1,coeff2,rho,grdnts)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE:  get_gradient_of_rho
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)  :: jsm, mddiff, kmesh, mesh
  real(kind=DP),intent(in),dimension(kmesh,jsm)        :: a
  real(kind=DP),intent(in),dimension(kmesh)            :: rad
  real(kind=DP),intent(in)                             :: h
  real(kind=DP),           dimension(kmesh,0:mddiff)   :: fdiff
  real(kind=DP),           dimension(0:mddiff,mddiff)  :: coeff1
  real(kind=DP),           dimension(0:mddiff,2:mddiff):: coeff2
  real(kind=DP),           dimension(kmesh,jsm)      :: rho
  real(kind=DP),intent(out),dimension(kmesh,jsm*3) :: grdnts

  integer, parameter :: execut = 2, modeex = execut &
       &, abslut = 1, normal = 2
  integer :: nddiff, npcntr, nddriv, cpmod, ngp, i

  nddiff = 8
  if(mddiff < nddiff) nddiff = mddiff
  npcntr = (nddiff+1)/2
#ifdef GGA_CHECK
  write(6,'(" -- get_gradient_of_rho -- ")')
#endif
  do i = 1, jsm
     nddriv = 1
     call getroh(a(1,i),kmesh,mesh,rad,rho(1,i))
#ifdef GGA_CHECK
     do j = 1, 10
        write(6,'(" ! ( ",i5,") rho = ",d16.8)') j, rho(j,i)
     end do
#endif
     cpmod = normal
     call cpval(rho(1,i),kmesh,mesh,cpmod,1,fdiff(1,0))
     call gdiffs(h,rad,fdiff,kmesh,mesh,mddiff,nddiff)
     !               output: fdiff
     ngp = 1 + (i-1)*3
     call gtgrad(h,rad,fdiff(1,1),kmesh,mesh,nddiff,nddriv &
          &       ,npcntr,coeff1 ,mddiff,grdnts(1,ngp),modeex)
     !               output: grdnts(*,1or4) = \nabla \rho
     !                                      = \partial \rho/\partial \r
     nddriv = 2
     call gtgrad(h,rad,fdiff(1,1),kmesh,mesh,nddiff,nddriv &
          &       ,npcntr,coeff2,mddiff,grdnts(1,ngp+1),modeex)
     !               output: grdnts(*,2or5) = \nabla^2 \rho
     cpmod = abslut
     call cpval(grdnts(1,ngp),kmesh,mesh,cpmod,1,fdiff(1,0))
     call gdiffs(h,rad,fdiff,kmesh,mesh,mddiff,nddiff)
     nddriv = 1
     call gtgrad(h,rad,fdiff(1,1),kmesh,mesh,nddiff,nddriv &
          &       ,npcntr,coeff1,mddiff,grdnts(1,ngp+2),modeex)
     !               output: grdnts(*,3or6) = \nabla (abs \nabla \rho)
     call cnggrd(grdnts(1,ngp),rad,kmesh,mesh )
     !               output: grdnts
  end do
end subroutine get_gradient_of_rho

subroutine gga_xc_potential_atomic(jsm,len_str,xctype,kmesh,nmesh,rho,rad&
     &,grdnts,exc,vxc,eps_chg,eps_grad)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: gga_xc_potential_atomic
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)         :: jsm,len_str
  character(len=len_str), intent(in) :: xctype
  integer, intent(in)         :: kmesh,nmesh
  real(kind=DP),intent(in),dimension(kmesh)   :: rho,rad
  real(kind=DP),intent(in),dimension(kmesh,3) :: grdnts
  real(kind=DP),intent(out),dimension(kmesh)  :: exc, vxc
  real(kind=DP),intent(in)                    :: eps_chg, eps_grad

  integer :: i
#ifdef GGA_CHECK
  integer :: j
#endif

  if(jsm == 2) then
     write(6,*) ' -- Warning -- '
     write(6,*) ' Sorry, this program( xc potential calculation) is not'&
          &     , ' corresponding to ferromagnetic case.'
  end if
#ifdef GGA_CHECK
  write(6,'(" -- gga_xc_potential_atomic --")')
  do j = 1, 10
     write(6,'(" ( ",i5,") rho = ",d12.4)') j, rho(j)
  end do
#endif

  exc = 0.d0; vxc = 0.d0
  if(         xctype == 'ggapw91' .or. xctype == 'GGAPW91'&
       & .or. xctype == 'ldapw91' .or. xctype == 'LDAPW91') then
     call ggaexch_pw91(nmesh,rho,grdnts(1,1),grdnts(1,2),grdnts(1,3),exc,vxc)  ! -> exc,vxc
     call ggacorl_pw91(nmesh,rho,grdnts(1,1),grdnts(1,2),grdnts(1,3),exc,vxc)  ! -> exc,vxc
  else if(    xctype == 'ggapbe ' .or. xctype == 'GGAPBE ' &
       & .or. xctype == 'ggapbex' .or. xctype == 'GGAPBEX' &
       & .or. xctype == 'ggapbey' .or. xctype == 'GGAPBEY' &
       & .or. xctype == 'vdwdf  ' .or. xctype == 'VDWDF  ' &
       & .or. xctype == 'ldapbe ' .or. xctype == 'LDAPBE ' &
       & .or. xctype == 'ggapbek' .or. xctype == 'GGAPBEK' &
       & .or. xctype == 'ldapbek' .or. xctype == 'LDAPBEK' &
       & .or. xctype == 'katopbe' .or. xctype == 'KATOPBE' ) then
     call ggaexch_pbe(nmesh,rho,grdnts(1,1),grdnts(1,2),grdnts(1,3),exc,vxc,eps_chg)  ! -> exc,vxc
     call ggacorl_pbe (nmesh,rho,grdnts(1,1),grdnts(1,2),grdnts(1,3),exc,vxc,eps_chg,eps_grad)  ! -> exc,vxc

! ================================= KT_add =================== 13.0A
  else if (   xctype == 'rpbe   ' .or. xctype == 'RPBE   ' ) then
     call ggaexch_rpbe(nmesh,rho,grdnts(1,1),grdnts(1,2),grdnts(1,3),exc,vxc,eps_chg)
     call ggacorl_pbe (nmesh,rho,grdnts(1,1),grdnts(1,2),grdnts(1,3),exc,vxc,eps_chg,&
          &            eps_grad)
! ============================================================= 13.0A

  else if(xctype == 'ggabp  ' .or. xctype == 'GGABP  ') then
     call ggabp(nmesh,rho,grdnts(1,1),grdnts(1,2),grdnts(1,3),exc,vxc)
  else
     write(6,'(" !! xctype does not match with any term")')
     call phase_error_with_msg(6, ' xctype is not properly set',__LINE__,__FILE__)
  end if

#ifdef GGA_CHECK
  write(6,'(" --- after gga[exch,corl]* in ( gga_xc_potential_atomic()) ---")')
#endif
  do i = 1, nmesh
     vxc(i) = vxc(i)*rad(i)
#ifdef GGA_CHECK
     if(i < 10 .or. i > nmesh-10+1 ) &
          & write(6,'(" ( ",i4," ) radr, rho, vxc = ", 3d18.8)') i, rad(i), rho(i), vxc(i)
#endif
  end do
  !  write(6,*) '! after vxc'

end subroutine gga_xc_potential_atomic

subroutine ggaexch_pw91(nmesh,rho,grad1,grad2,grad3,exc,vxc )
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: ggaexch_pw91
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP, PAI
  implicit none
  integer, intent(in)                         :: nmesh
  real(kind=DP),intent(in),   dimension(nmesh):: rho,grad1,grad2,grad3
  real(kind=DP),intent(inout),dimension(nmesh):: exc,vxc
  
!                           @(#)ggaexch_pw91.f 9.1 97/05/08 14:48:26 
!*********************************************************************
!      PW'91 GGA exchange energy and potential.
!      The original code came from Jim White and David Bird.
!           Ref: J.P.Perdew, in "Electronic Structure of Solids '91",
!                                edited by P.Ziesche and E.Eschrig
!                                   (Akademie Verlag, Berlin,1991).
!        modified  by T.Yamasaki  94/11/14
!        modified  by Y.Morikawa '94/11/30 at JRCAT Tsukuba.
!        modified  by T.Yamasaki '94/12/01 
! inputs:
!     nmesh      : number of grid points
!     rho     : charge density
!     grad1   : abs(\nabla(rho))
!     grad2   : \div(\nabla(rho))
!     grad3   : 1/2*(\nabla(rho))\cdot(\nabla(\abs(\nabla(rho))**2))
!                                           /(\abs(\nabla(rho))**2
! outputs:
!     exc     : exchange energy
!     vxc     : exchange potential
!     
!*********************************************************************
  real(kind=DP),parameter :: a1  = 0.19645d0,    a2 = 0.27430d0 &
       &                   , a3  = 0.15084d0,    a4 = 100.0d0   &
       &                   , ax  = -0.7385588d0,  a = 7.7956d0,  b1 = 0.004d0 &
       &                   , thrd = 0.33333333333d0, thrd4 = 1.333333333333333d0
!!$  real(kind=DP),parameter :: thpith = 3.0936677262801d0

  real(kind=DP)  :: d,dd,fk,s,fac,s2,s3,s4,p0,p1,p2,p3,p4,f,ex &
       &          , p5,p6,fs,exc0&
       &          , p5s,p6s,p7,fs2,vx1,vx2,vx3,vx4,vx5, vxall
  integer  ::  i

  do  i = 1, nmesh
     d = rho(i)
     if(d < 1.d-40) cycle
     dd = grad1(i)
     fk = (3*PAI*PAI*d)**thrd
!!$     sk = dsqrt(4*fk/PAI)
     if(d > 1.d-20) then
        s = dd/(d*fk*2.d0)
     else
        s = 0.d0
     endif
!  -------------------------------------
     fac = ax*d**thrd
     s2  = s*s
     s3  = s2*s
     s4  = s2*s2
     p0  = 1.d0/(dsqrt(1.d0+a*a*s2))
     p1  = dlog(a*s+1.d0/p0)
     p2  = dexp(-a4*s2)
     p3  = 1.d0/(1.d0+a1*s*p1+b1*s4)
     p4  = 1.d0 + a1*s*p1 + (a2-a3*p2)*s2
     f   = p3*p4

!  exchange energy

     ex  = fac*f*d

     p5  = 2.d0*(s*(a2-a3*p2)+a3*a4*s3*p2 - 2.d0*b1*s3)
     p6  = (a1*(p1+a*s*p0)+4.d0*b1*s3)*((a2-a3*p2)*s2-b1*s4)
     fs  = (p5*p3-p6*p3*p3)
!!$     exd = thrd4*fac*(f-s*fs)
!!$     exdd = ax*fs*0.5d0/thpith
!!$     p5  = 2*(s*(a2-a3*p2)+a3*a4*s3*p2 - 2*b1*s3)
!!$     p6  = (a1*(p1+a*s*p0)+4*b1*s3)*((a2-a3*p2)*s2-b1*s4)
!!$     fs  = (p5*p3-p6*p3*p3)
!!$     exd = thrd4*fac*(f-s*fs)
!!$     exdd = ax*fs*0.5d0/thpith

     p5s = 2.d0*(a2-a3*p2*(1.d0-2.d0*a4*s2*(1.d0-a4*s2))&
          &              +3.d0*s2*(a3*a4*p2-2.d0*b1))

     p6s = (a1*a*p0*(2.d0-(a*s*p0)**2)+12.d0*b1*s2) &
          &       *((a2-a3*p2)*s2-b1*s4) &
          &       +(a1*(p1+a*s*p0)+4*b1*s3)*2.d0*s&
          &        *(a3*a4*s2*p2+(a2-a3*p2)-2.d0*b1*s2)

!!$     p5s = 2*(a2-a3*p2*(1-2*a4*s2*(1-a4*s2))+3*s2*(a3*a4*p2-2*b1))
!!$
!!$     p6s = (a1*a*p0*(2-(a*s*p0)**2)+12*b1*s2)*((a2-a3*p2)*s2-b1*s4) &
!!$          &   +(a1*(p1+a*s*p0)+4*b1*s3)*2*s*(a3*a4*s2*p2+(a2-a3*p2)-2*b1*s2)

     p7  =  a1*(p1+a*s*p0)+4*b1*s3
     fs2 = p3*(p5s-p3*((p5*p7+p6s)-p3*2*p6))

!  exchange potential

     vx1 =  thrd4*fac*f
     vx2 =  thrd4*fac*s2*fs2
     vx3 = -thrd4*fac*s*fs
     if(dd < 1.d-20) then
        vx4 = 0.d0
        vx5 = 0.d0
     else
        vx4 = -grad2(i)*d/dd**2 * fac*s*fs
        vx5 =  grad3(i)*d/dd**2*fac*(s*fs-s2*fs2)
     endif
     vxall = vx1 + vx2 + vx3 + vx4 + vx5
     vxc(i) = vxc(i) + vxall
!------------------------------------------     
     exc0 = ex
     exc(i)  = exc(i) + exc0
  end do
end subroutine ggaexch_pw91

subroutine ggacorl_pw91(nmesh,rho,grad1,grad2,grad3,exc,vxc)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: ggacorl_pw91
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP, PAI
  implicit none
  integer, intent(in)                         :: nmesh
  real(kind=DP),intent(in),   dimension(nmesh):: rho,grad1,grad2,grad3
  real(kind=DP),intent(inout),dimension(nmesh):: exc,vxc
!                           @(#)ggacorl_pw91.f 9.1 97/05/08 14:48:26 
!*********************************************************************
!      PW'91 GGA correlation energy and potential.
!      The original code came from Jim White and David Bird.
!           Ref: J.P.Perdew, in "Electronic Structure of Solids '91",
!                                edited by P.Ziesche and E.Eschrig
!                                   (Akademie Verlag, Berlin,1991).
!                LDA part is given in
!                J.P.Perdew and Y.Wang, PRB45, 13244 (1992).
!
!        modified  by T.Yamasaki  94/11/15
!        modified  by Y.Morikawa '94/11/30 at JRCAT Tsukuba.
!        modified  by T.Yamasaki '94/12/01 
! inputs:
!     nmesh      : number of grid points
!     rho     : charge density
!     grad1   : abs(\nabla(rho))
!     grad2   : \div(\nabla(rho))
!     grad3   : 1/2*(\nabla(rho))\cdot(\nabla(\abs(\nabla(rho))**2))
!                                           /(\abs(\nabla(rho))**2
! outputs:
!     exc        : correlation energy
!     vxc   : correlation potential
!*********************************************************************
!  parameters for LDA part
  real(kind=DP),parameter :: a   = 0.0310907d0&
       &,                    a1  = 0.21370d0  &
       &,                    b1  = 7.5957d0   &
       &,                    b2n = 3.5876d0   &
       &,                    b3  = 1.6382d0   &
       &,                    b4  = 0.49294d0  &
       &,                    p   = 1.00d0     &
       &,                    p1  = p+1.d0 

!  parameters for GGA part
  real(kind=DP),parameter :: xnu = 15.75592d0,     cc0 = 0.004235d0 &
       &,                    cx  = -0.001667212d0, alf = 0.09d0     &
       &,                    c1  = 0.002568d0,      c2 = 0.023266d0 &
       &,                    c3  = 7.389d-6,        c4 = 8.723d0    &
       &,                    c5  = 0.472d0,         c6 = 7.389d-2   &
       &,                    a4  = 100.d0 &
       &,                  thrd  = 0.333333333333d0, sixth7 = 1.1666666666666666d0

  integer :: i
  real(kind=DP) :: bet,delt,g,g3,g4,facpon &
       &   , d,dd,rs,fk,sk,t,q0,q1,q2,q3,q4,q5,q6,q7,q8 &
       &   , rs12,rs32,rsp,eu,eurs,pon,b,b2,t2,t4,t6,rs2,rs3 &
       &   , h0,h0t,h0b,h0bt,h0rs,h0rst,h0t2 &
       &   , h1,h1t,h1rs,h1rst,h1t2,cc,r0,r1,r2,r3,coeff &
       &   , ccrs,r1rs,h &
       &   , v1,v2,v3,v4,v5

!--------------------------
  bet  = xnu*cc0
  delt = 2*alf/bet
  g    = 1.d0
  g3   = g**3
  g4  =  g3*g
  facpon = -delt/(g3*bet)
!---------------------------
  do i = 1, nmesh
     d = rho(i)
     if(d < 1.d-15) cycle

     dd = grad1(i)
	
     rs = (0.75d0/(PAI*d))**thrd
     fk = (3*PAI*PAI*d)**thrd
     sk = dsqrt(4*fk/PAI)
     if(d > 1.d-15) then
        t = dd/(d*sk*2)
     else
        t = 0.d0
     endif

!   LDA part

     q0   = -2*a*(1+a1*rs)
     rs12 = dsqrt(rs)
     rs32 = rs12**3
     rsp  = rs**p
     q1   = 2*a*(b1*rs12+b2n*rs+b3*rs32+b4*rs*rsp)
     q2   = log(1+1.d0/q1)
     q3   = a*(b1/rs12+2*b2n+3*b3*rs12+2*b4*p1*rsp)
     
     eu   = q0*q2
     eurs = -2*a*a1*q2 - q0*q3/(q1**2 + q1)

!  GGA part

!  H0 part

     pon = facpon*eu
     b   = delt/(dexp(pon)-1)
     b2  = b*b
     t2  = t*t
     t4  = t2*t2
     t6  = t4*t2
     rs2 = rs*rs
     rs3 = rs2*rs
     q4  = 1+b*t2
     q5  = 1+b*t2+b2*t4
     q8  = q5*q5+delt*q4*q5*t2

!!$     q8s = 2*t*(1+2*b*t2)*(2*b*q5+delt*(1+2*b*t2*(1+b*t2)))

     h0  = g3*(bet/delt)*dlog(1+delt*q4*t2/q5)
     h0t = 2*bet*t*(1+2*b*t2)/q8
     h0b = -bet*t6*(2*b+b2*t2)/q8

     h0bt = -bet*b*g3*t4*t/q8**2 &
          &   *(4*q5*(3*(1+b*t2)-5*b2*t4)+2*delt*t2*(4+7*b*t2+4*b2*t4))

     h0rs = h0b*b*eurs*(b+delt)/bet

     h0rst= h0bt*b*eurs*(b+delt)/bet
     h0t2 = 2*bet*g3/q8**2 &
          &    *(q5     *(1+b*t2*(3+b*t2*(-9-b*t2*10))) &
          &     -delt*t2*(1+b*t2*(4+b*t2*(14+b*t2*(19+b*t2*10)))))

!  H1 part

     q6  = c1+c2*rs+c3*rs2
     q7  = 1+c4*rs+c5*rs2+c6*rs3
     cc  = -cx + q6/q7
     r0  = (sk/fk)**2
     r1  = a4*r0*g4
     coeff = cc-cc0-3*cx/7.d0
     r2  = xnu*coeff*g3
     r3  = dexp(-r1*t2)

     h1  = r3*r2*t2
     h1t  = 2*r3*r2*t*(1-r1*t2)
     ccrs = (c2+2*c3*rs)/q7-q6*(c4+2*c5*rs+3*c6*rs2)/q7**2
     r1rs = 100*r0/rs
     h1rs = xnu*t2*r3*(ccrs - coeff*t2*r1rs)

     h1rst= xnu*g3*2*t*r3 &
          &  *((1-r1*t2)*(ccrs-t2*coeff*r1rs)-t2*coeff*r1rs)
     h1t2 = 2*r2*r3*(1-r1*t2*(5-2*r1*t2))

!======================================================
     h   = h0+h1
!!$     ht   = h0t+h1t
!!$     hrs  = h0rs + h1rs

!!$     hrst = h0rst + h1rst
!!$     ht2  = h0t2  + h1t2

!!$     ec1  = d*h
!!$     ec1d = h-thrd*rs*hrs-sixth7*t*ht
!!$     ec1dd = 0.5d0*ht/sk

     v1 = eu+h0+h1-thrd*rs*(eurs+h0rs+h1rs)
     v2 = sixth7*t2*(h0t2+h1t2)
     v3 = -t*((h0t+h1t)-thrd*rs*(h0rst+h1rst))
     if(dd < 1.d-15) then
        v4 = 0.d0
        v5 = 0.d0
     else
        v4 = -grad2(i)/dd**2*t*d*(h0t+h1t)
        v5 =  grad3(i)      /dd**2*d*t*(h0t+h1t-t*(h0t2+h1t2))
     endif

     vxc(i) = vxc(i) + v1+v2+v3+v4+v5
     exc(i) = exc(i) + (eu+h)*d
  end do
end subroutine ggacorl_pw91

! ========================================== KT_add =============== 13.0A
subroutine ggaexch_rpbe(nmesh,rho,grad1,grad2,grad3,exc,vxc,eps_chg )

  use m_Const_Parameters, only : DP, PAI
  implicit none
  integer, intent(in)                         :: nmesh
  real(kind=DP),intent(in),   dimension(nmesh):: rho,grad1,grad2,grad3
  real(kind=DP),intent(inout),dimension(nmesh):: exc,vxc
  real(kind=DP),intent(in)                    :: eps_chg

!***
!       Coded based on the subroutine of ggaexch_pw91 by K. Kato 99/06/25
!       Modified into F90 by T. Yamasaki 99/06/26
! inputs:
!     nmesh      : number of grid points
!     rho     : charge density
!     grad1   : abs(\nabla(rho))
!     grad2   : \div(\nabla(rho))
!     grad3   : 1/2*(\nabla(rho))\cdot(\nabla(\abs(\nabla(rho))**2))
!                                           /(\abs(\nabla(rho))**2
! outputs:
!     exc     : exchange energy
!     vxc     : exchange potential
!
!*********************************************************************
  real(kind=DP),parameter :: ax  = -0.738558766382022405884230032680836d0 &
!!$       &                   , ax  = -0.7385588d0,  a = 7.7956d0,  b1 = 0.004d0 &
       &                   , thrd = 1.d0/3.d0, thrd4 = 4.d0/3.d0 &
!!$  real(kind=DP),parameter :: thpith = 3.0936677262801d0 &
!!$       &                   , cupa = 0.804,      yum = 0.21951
       &                   , cupa = 0.8040d0,      yum = 0.2195149727645171d0

  real(kind=DP)  :: d,dd,fk,s,fac,s2,f,ex,fs,exc0&
       &          , fs2,vx1,vx2,vx3,vx4,vx5, vxall,yums2, ctmp1
  integer  ::  i

#ifdef GGA_CHECK
  write(6,'(" -- ggaexch_rpbe --")')
#endif
  do  i = 1, nmesh
#ifdef GGA_ATOMIC_WITH_NEW_GNCPP
     d  = max(rho(i),eps_chg)
#else
     d = rho(i)
     if(d < 1.d-40) cycle
#endif
     dd = grad1(i)
     fk = (3*PAI*PAI*d)**thrd
!!$     sk = dsqrt(4*fk/PAI)
     if(d > 1.d-20) then
        s = dd/(d*fk*2)
     else
        s = 0.d0
     endif
!  -------------------------------------
     fac = ax*d**thrd
     s2  = s*s
     yums2 = yum * s2

     ctmp1 = exp( -yums2 /cupa )
     f = 1.0D0 + cupa *( 1.0D0 - ctmp1 )
     fs = 2.0D0 *yum * s * ctmp1
     fs2 = 2.0D0 *yum *( 1.0D0 - 2.0D0 *yums2/cupa ) *ctmp1

!  exchange energy

     ex  = fac*f*d
!!$     exd = thrd4*fac*(f-s*fs)
!!$     exdd = ax*fs*0.5/thpith

!  exchange potential

     vx1 =  thrd4*fac*f
     vx2 =  thrd4*fac*s2*fs2
     vx3 = -thrd4*fac*s*fs
     if(dd < 1.d-20) then
        vx4 = 0.d0
        vx5 = 0.d0
     else
        vx4 = -grad2(i)*d/dd**2 * fac*s*fs
        vx5 =  grad3(i)*d/dd**2*fac*(s*fs-s2*fs2)
     endif
     vxall = vx1 + vx2 + vx3 + vx4 + vx5
     vxc(i) = vxc(i) + vxall
#ifdef GGA_CHECK
     write(6,'(" (",i4,") rho, vxc = ",2d18.8," (ggaexch_rpbe)")') i, rho(i), vxc(i)
#endif
!------------------------------------------
     exc0 = ex
     exc(i)  = exc(i) + exc0
  end do
end subroutine ggaexch_rpbe
! ==============================================================================13.0A

subroutine ggaexch_pbe(nmesh,rho,grad1,grad2,grad3,exc,vxc,eps_chg )
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE:  ggaexch_pbe
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP, PAI
  implicit none
  integer, intent(in)                         :: nmesh
  real(kind=DP),intent(in),   dimension(nmesh):: rho,grad1,grad2,grad3
  real(kind=DP),intent(inout),dimension(nmesh):: exc,vxc
  real(kind=DP),intent(in)                    :: eps_chg
  
!***
!       Coded based on the subroutine of ggaexch_pw91 by K. Kato 99/06/25 
!       Modified into F90 by T. Yamasaki 99/06/26
! inputs:
!     nmesh      : number of grid points
!     rho     : charge density
!     grad1   : abs(\nabla(rho))
!     grad2   : \div(\nabla(rho))
!     grad3   : 1/2*(\nabla(rho))\cdot(\nabla(\abs(\nabla(rho))**2))
!                                           /(\abs(\nabla(rho))**2
! outputs:
!     exc     : exchange energy
!     vxc     : exchange potential
!     
!*********************************************************************
  real(kind=DP),parameter :: ax  = -0.738558766382022405884230032680836d0 &
!!$       &                   , ax  = -0.7385588d0,  a = 7.7956d0,  b1 = 0.004d0 &
       &                   , thrd = 1.d0/3.d0, thrd4 = 4.d0/3.d0 &
!!$  real(kind=DP),parameter :: thpith = 3.0936677262801d0 &
!!$       &                   , cupa = 0.804,      yum = 0.21951
       &                   , cupa = 0.8040d0,      yum = 0.2195149727645171d0

  real(kind=DP)  :: d,dd,fk,s,fac,s2,f,ex,fs,exc0&
       &          , fs2,vx1,vx2,vx3,vx4,vx5, vxall,yums2
  integer  ::  i

#ifdef GGA_CHECK
  write(6,'(" -- ggaexch_pbe --")')
#endif
  do  i = 1, nmesh
#ifdef GGA_ATOMIC_WITH_NEW_GNCPP
     d  = max(rho(i),eps_chg)
#else
     d = rho(i)
     if(d < 1.d-40) cycle
#endif
     dd = grad1(i)
     fk = (3*PAI*PAI*d)**thrd
!!$     sk = dsqrt(4*fk/PAI)
     if(d > 1.d-20) then
        s = dd/(d*fk*2)
     else
        s = 0.d0
     endif
!  -------------------------------------
     fac = ax*d**thrd
     s2  = s*s
     yums2 = yum * s2
     f   = 1 + cupa*yums2/(cupa+yums2)
!  exchange energy

     ex  = fac*f*d
     fs = 2 * cupa ** 2 * yum * s / ( cupa + yums2 ) ** 2
!!$     exd = thrd4*fac*(f-s*fs)
!!$     exdd = ax*fs*0.5/thpith

     fs2 = 2 * cupa **2 * (cupa*yum - 3*yum*yums2)/(cupa+yums2)**3

!  exchange potential

     vx1 =  thrd4*fac*f
     vx2 =  thrd4*fac*s2*fs2
     vx3 = -thrd4*fac*s*fs
     if(dd < 1.d-20) then
        vx4 = 0.d0
        vx5 = 0.d0
     else
        vx4 = -grad2(i)*d/dd**2 * fac*s*fs
        vx5 =  grad3(i)*d/dd**2*fac*(s*fs-s2*fs2)
     endif
     vxall = vx1 + vx2 + vx3 + vx4 + vx5
     vxc(i) = vxc(i) + vxall
#ifdef GGA_CHECK
     write(6,'(" (",i4,") rho, vxc = ",2d18.8," (ggaexch_pbe)")') i, rho(i), vxc(i)
#endif
!------------------------------------------     
     exc0 = ex
     exc(i)  = exc(i) + exc0
  end do
end subroutine ggaexch_pbe

subroutine ggacorl_pbe(nmesh,rho,grad1,grad2,grad3,exc,vxc,eps_chg,eps_grad)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE:  ggacorl_pbe
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP, PAI
  implicit none
  integer, intent(in)                         :: nmesh
  real(kind=DP),intent(in),   dimension(nmesh):: rho,grad1,grad2,grad3
  real(kind=DP),intent(inout),dimension(nmesh):: exc,vxc
  real(kind=DP),intent(in)                    :: eps_chg, eps_grad
!***
!       Coded based on the subroutine of ggaexch_pw91 by K. Kato 99/06/25 
!       Modified into F90 by T. Yamasaki 99/06/26
!
! inputs:
!     nmesh      : number of grid points
!     rho     : charge density
!     grad1   : abs(\nabla(rho))
!     grad2   : \div(\nabla(rho))
!     grad3   : 1/2*(\nabla(rho))\cdot(\nabla(\abs(\nabla(rho))**2))
!                                           /(\abs(\nabla(rho))**2
! outputs:
!     exc        : correlation energy
!     vxc   : correlation potential
!*********************************************************************
!  parameters for LDA part
  real(kind=DP),parameter :: a   = 0.0310907d0&
       &,                    a1  = 0.21370d0  &
       &,                    b1  = 7.5957d0   &
       &,                    b2n = 3.5876d0   &
       &,                    b3  = 1.6382d0   &
       &,                    b4  = 0.49294d0  &
       &,                    p   = 1.00d0     &
       &,                    p1  = p+1.d0 

!  parameters for GGA part
!!$  real(kind=DP),parameter :: xnu = 15.75592d0,      cc0 = 0.004235d0 &
  real(kind=DP),parameter :: bet = 0.06672455060314922d0 &
       &                   , thrd = 1.d0/3.d0, sixth7 = 7.d0/6.d0 &
       &                   , gamma = 0.03109069086965489503494086371273d0

  real(kind=DP) :: delt,g,g3,facpon &
       &   , d,dd,rs,fk,sk,t,q0,q1,q2,q3,q4,q5,q8 &
       &   , rs12,rs32,rsp,eu,eurs,pon,b,b2,t2,t4,t6 &
       &   , h0,h0t,h0b,h0bt,h0rs,h0rst,h0t2,h &
       &   , v1,v2,v3,v4,v5
  integer :: i
#ifdef GGA_CHECK
  write(6,'(" -- ggacorl_pbe --")')
#endif

!!$  bet = xnu*cc0
  delt = bet / gamma
  g = 1.d0
  g3  = g**3
  facpon = -delt/(g3*bet)
!---------------------------
  do i = 1, nmesh
#ifdef GGA_ATOMIC_WITH_NEW_GNCPP
     d  = max(rho(i), eps_chg)
     dd = max(grad1(i), eps_grad)
#else
     d  = rho(i)
     if(d < 1.d-15) cycle
     dd = grad1(i)
#endif

	
     rs = (0.75d0/(PAI*d))**thrd
     fk = (3*PAI*PAI*d)**thrd
     sk = dsqrt(4*fk/PAI)
#ifdef GGA_ATOMIC_WITH_NEW_GNCPP
     if(rho(i) > eps_chg) then
#else
     if(d > 1.d-15) then
#endif
        t = dd/(d*sk*2)
     else
        t = 0.d0
     endif

!   LDA part

     q0   = -2*a*(1+a1*rs)
     rs12 = dsqrt(rs)
     rs32 = rs12**3
     rsp  = rs**p
     q1   = 2*a*(b1*rs12+b2n*rs+b3*rs32+b4*rs*rsp)
     q2   = log(1+1.d0/q1)
     q3   = a*(b1/rs12+2*b2n+3*b3*rs12+2*b4*p1*rsp)
     
     eu   = q0*q2
     eurs = -2*a*a1*q2 - q0*q3/(q1**2 + q1)

!  GGA part

!  H0 part

     pon = facpon*eu
     b   = delt/(dexp(pon)-1)
     b2  = b*b
     t2  = t*t
     t4  = t2*t2
     t6  = t4*t2
!!$     rs2 = rs*rs
!!$     rs3 = rs2*rs
     q4  = 1+b*t2
     q5  = 1+b*t2+b2*t4
     q8  = q5*q5+delt*q4*q5*t2

!!$     q8s = 2*t*(1+2*b*t2)*(2*b*q5+delt*(1+2*b*t2*(1+b*t2)))

     h0  = g3*gamma*dlog(1+delt*q4*t2/q5)
#ifdef GGA_ATOMIC_WITH_NEW_GNCPP
     if(abs(rho(i)) < eps_chg) h0 = 0.d0
#endif
     h0t = 2*bet*t*(1+2*b*t2)/q8
     h0b = -bet*t6*(2*b+b2*t2)/q8

     h0bt = -bet*b*g3*t4*t/q8**2 &
          &   *(4*q5*(3*(1+b*t2)-5*b2*t4)+2*delt*t2*(4+7*b*t2+4*b2*t4))

     h0rs = h0b*b*eurs*(b+delt)/bet

     h0rst= h0bt*b*eurs*(b+delt)/bet
     h0t2 = 2*bet*g3/q8**2 &
          &    *(q5     *(1+b*t2*(3+b*t2*(-9-b*t2*10))) &
          &       -delt*t2*(1+b*t2*(4+b*t2*(14+b*t2 &
          &                          *(19+b*t2*10)))))

!  H1 part
     h   = h0
!!$     ht   = h0t
!!$     hrs  = h0rs

     v1 = eu+h0-thrd*rs*(eurs+h0rs)
     v2 = sixth7*t2*h0t2
     v3 = -t*(h0t-thrd*rs*(h0rst))

#ifdef GGA_ATOMIC_WITH_NEW_GNCPP
     if(grad1(i) < eps_grad) then
#else
     if(dd < 1.d-15) then
#endif
        v4 = 0.d0
        v5 = 0.d0
     else
        v4 = -grad2(i)/dd**2*t*d*h0t
        v5 =  grad3(i)/dd**2*d*t*(h0t-t*h0t2)
     endif

     vxc(i) = vxc(i) + v1+v2+v3+v4+v5
     exc(i) = exc(i) + (eu+h)*d
#ifdef GGA_CHECK
     write(6,'(" (",i4,") rho,vxc,v4,v5 = ",4d18.6," ggacorl_pbe")') i, rho(i), vxc(i), v4, v5
#endif
  end do
end subroutine ggacorl_pbe

subroutine ggabp(nmesh,rho,grad1,grad2,grad3,exc,vxc)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: ggabp
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP, PAI
  implicit none
  integer, intent(in)                         :: nmesh
  real(kind=DP),intent(in),   dimension(nmesh):: rho,grad1,grad2,grad3
  real(kind=DP),intent(inout),dimension(nmesh):: exc,vxc
!        by T. Yamasaki on  Dec. 1994
  real(kind=DP),parameter :: a1 = 0.19645d0,   a2 = 0.27430d0 &
       &,                    a3 = 0.15084d0,   a4 = 100.0d0   &
       &,                    ax = -0.7385588d0, a = 7.7956d0,  b1 = 0.004d0 &
       &,                  thrd = 0.33333333333d0, thrd4 = 1.333333333333333d0
!!$       &,                thpith = 3.0936677262801d0

  integer :: i
  real(kind=DP) :: d,dd,fk,s,fac,s2,s3,s4,p0,p1,p2,p3,p4 &
       &         , f,ex,p5,p6,fs &
       &         , p5s,p6s,p7,fs2,vx1,vx2,vx3,vx4,vx5

  do i = 1, nmesh

     d = rho(i)
     dd = grad1(i)
     fk = (3*PAI*PAI*d)**thrd
!!$     sk = dsqrt(4*fk/PAI)
     if(d > 1.d-20) then
        s = dd/(d*fk*2)
     else
        s = 0.d0
     endif
!-------------------------------------
     fac = ax*d**thrd
     s2  = s*s
     s3  = s2*s
     s4  = s2*s2
     p0  = 1.d0/(dsqrt(1+a*a*s2))
     p1  = dlog(a*s+1/p0)
     p2  = dexp(-a4*s2)
     p3  = 1.d0/(1+a1*s*p1+b1*s4)
     p4  = 1+a1*s*p1 + (a2-a3*p2)*s2
     f   = p3*p4

!  exchange energy

     ex  = fac*f*d

     p5  = 2*(s*(a2-a3*p2)+a3*a4*s3*p2 - 2*b1*s3)
     p6  = (a1*(p1+a*s*p0)+4*b1*s3)*((a2-a3*p2)*s2-b1*s4)
     fs  = (p5*p3-p6*p3*p3)
!!$     exd = thrd4*fac*(f-s*fs)
!!$     exdd = ax*fs*0.5d0/thpith

     p5s = 2*(a2-a3*p2*(1-2*a4*s2*(1-a4*s2)) &
          &              +3*s2*(a3*a4*p2-2*b1))
     p6s = (a1*a*p0*(2-(a*s*p0)**2)+12*b1*s2) &
          &       *((a2-a3*p2)*s2-b1*s4) &
          &       +(a1*(p1+a*s*p0)+4*b1*s3)*2*s &
          &        *(a3*a4*s2*p2+(a2-a3*p2)-2*b1*s2)

     p7  =  a1*(p1+a*s*p0)+4*b1*s3
     fs2 = p3*(p5s-p3*((p5*p7+p6s)-p3*2*p6))

!  exchange potential

     vx1 =  thrd4*fac*f
     vx2 =  thrd4*fac*s2*fs2
     vx3 = -thrd4*fac*s*fs

     if(dd.lt.1.d-20) then
        vx4 = 0.d0
        vx5 = 0.d0
     else
        vx4 = -grad2(i)*d/dd**2 * fac*s*fs
        vx5 =        grad3(i)*d/dd**2*fac*(s*fs-s2*fs2)
     end if

     vxc(i) = vx1+vx2+vx3+vx4+vx5
     exc(i) = ex
  end do
end subroutine ggabp

subroutine rmeshs(k,m,xh,rmax,radr,h)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: rmeshs
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)                      :: k, m
  real(kind=DP), intent(in)                :: xh, rmax
  real(kind=DP), intent(out), dimension(k) :: radr
  real(kind=DP), intent(out)               :: h

  integer i

  h = 1.d0/xh
  do i = 1, m
     radr(i) = rmax*exp(h*dble(i-m))
  enddo

end subroutine rmeshs

subroutine coef_simpson_integration(k,m,xh,radr,wos)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: coef_simpson_integration
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)                      :: k, m
  real(kind=DP), intent(in)                :: xh
  real(kind=DP), intent(in),  dimension(k) :: radr
  real(kind=DP), intent(out), dimension(k) :: wos

  integer i

!!$  h = 1.d0/xh
  do i = 2, m-1, 2
     wos(i  ) = radr(i  )*4/(xh*3)
     wos(i+1) = radr(i+1)*2/(xh*3)
  enddo
  wos(1) = radr(1)    /(xh*3)
  wos(m) = radr(m)    /(xh*3)

end subroutine coef_simpson_integration

subroutine read_itau_ivanl(nfp,nfout,iloc,il,itau,ivanl,ipri)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE:  read_itau_ivanl
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  implicit none
  integer, intent(in)  :: nfp, nfout
  integer, intent(in)  :: iloc, il, ipri
  integer, intent(out) :: itau
  integer, intent(out) :: ivanl

  integer il_f

  read(nfp,*)      il_f, itau, ivanl
  if(ipri>=2) write(nfout,201) il_f, itau, ivanl
201 format(' !bPP ',3i4,'  :   IL, ITAU(IL1),IVANL(IL1) ')
  if(il_f /= il) then
     !write(nfout,*) ' IL_F in the psfile is invalid'
     !stop
     call phase_error_with_msg(nfout, 'IL_F in the psfile is invalid',__LINE__,__FILE__)
  endif
  if( (  il_f == iloc                         .and. &
     &       (itau /= 1 .or. ivanl /= 0) )           ) then
     !write(nfout,*) ' iloc,il1,itau,ivanl is not good'
     !stop
     call phase_error_with_msg(nfout, 'iloc,il1,itau,ivanl is not good',__LINE__,__FILE__)
  endif
end subroutine read_itau_ivanl

subroutine read_tau_eps_nrc_mord(nfp,nfout,kord,il,t1,eps,nrc,mord)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: read_tau_eps_nrc_mord
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)  :: nfp,nfout,kord
  integer, intent(in)  :: il
  integer, intent(in)  :: t1
  real(kind=DP), intent(out) :: eps
  integer, intent(out) :: nrc, mord

  integer              :: il0, tau0

  read(nfp,*) il0,tau0,eps,nrc,mord
!!$  write(nfout,202) il0, tau0, eps,nrc,mord
!!$202 format(' ',2i4,f12.8,2i6,' : il, ntau,  eps, nrc, mord')
  if(il0 /= il .or. tau0 /= t1) then
     write(nfout,*) ' il0, il = ', il0, il, ' t1, tau0 = ', t1, tau0
     !stop
     call phase_error_with_msg(nfout,'invalid il0 or ta0',__LINE__,__FILE__)
  endif
  if(mord > kord) then
     write(nfout,*) ' ! MORD is too large : mord = ', mord, ' kord = ', kord
!     stop
     call phase_error_with_msg(nfout,'MORD is too large',__LINE__,__FILE__)
  endif
end subroutine read_tau_eps_nrc_mord

subroutine read_tau_eps_nrc_mord_nrc0(nfp,nfout,kord,il,t1,eps,nrc,mord,nrc0)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: read_tau_eps_nrc_mord_nrc0
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)  :: nfp,nfout,kord
  integer, intent(in)  :: il
  integer, intent(in)  :: t1
  real(kind=DP), intent(out) :: eps
  integer, intent(out) :: nrc, mord, nrc0

  integer              :: il0, tau0, ios

  read(nfp,*,iostat=ios) il0,tau0,eps,nrc,mord,nrc0
!!$  write(nfout,202) il0, tau0, eps,nrc,mord
!!$202 format(' ',2i4,f12.8,2i6,' : il, ntau,  eps, nrc, mord')
  if(il0 /= il .or. tau0 /= t1) then
     write(nfout,*) ' il0, il = ', il0, il, ' t1, tau0 = ', t1, tau0
     !stop
     call phase_error_with_msg(nfout, 'invalid il0 or tau0',__LINE__,__FILE__)
  endif
  if(mord > kord) then
     write(nfout,*) ' ! MORD is too large : mord = ', mord, ' kord = ', kord
     !stop
     call phase_error_with_msg(nfout, 'MORD is too large',__LINE__,__FILE__)
  endif
  if(ios.ne.0) nrc0=0
end subroutine read_tau_eps_nrc_mord_nrc0

subroutine read_nrc_mord(nfp,nfout,kord,mm,il1,tau1,il2,tau2,il3,nrc&
     & ,mord,ipripp)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: read_nrc_mord
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  implicit none
  integer, intent(in)  :: nfp,nfout,kord,mm,il1,tau1,il2,tau2,il3,ipripp
  integer, intent(out) :: nrc, mord

  integer :: l1, t1, l2, t2, l3
  character(len=80) :: str

  if(ipripp >= 2) write(nfout,'(" !bPP nfp,nfout,kord,mm,il1,tau1,il2,tau2,il3 = ",9i5)') &
       &  nfp,nfout,kord,mm,il1,tau1,il2,tau2,il3
  read(nfp,'(a80)') str
  if(ipripp >= 2) write(nfout,'(" !bPP ",a80)') str
  read(str,*)          L1,T1,L2,T2,L3,NRC,MORD
!!$  READ(NFP,*)          L1,T1,L2,T2,L3,NRC,MORD
  if(ipripp >= 2) WRITE(NFOUT,341) MM, L1,T1,L2,T2,L3,NRC,MORD
341 FORMAT(' !bPP ',6I4,I6,I4,' : MM, L1,T1,L2,T2,L3,NRC,MORD')
  IF(IL1.NE.L1.OR.TAU1.NE.T1.OR. &
       &        IL2.NE.L2.OR.TAU2.NE.T2.OR.IL3.NE.L3) THEN
     WRITE(NFOUT,'(" !bPP L1,T1,L2,T2,L3 ARE NOT CORRECT")')
     STOP
  END IF
  if(mord > kord) then
     write(nfout,*) ' ! MORD is too large : mord = ', mord, ' kord = ', kord
     !stop
     call phase_error_with_msg(nfout, 'MORD is too large',__LINE__,__FILE__)
  endif
end subroutine read_nrc_mord

subroutine read_chgpc_nrc_mord(ipri,nfp,nfout,chgpc,nrc,mord)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: read_chgpc_nrc_mord
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)        :: ipri, nfp, nfout
  real(kind=DP), intent(out) :: chgpc
  integer, intent(out)       :: nrc, mord

  read(nfp,*)      chgpc, nrc, mord
  if(ipri>=2) then
     write(nfout,400) chgpc, nrc, mord
400  format(' !bPP',' chgpc = ',f12.8,' nrc, mord = ',2i4)
  end if

end subroutine read_chgpc_nrc_mord

subroutine qitgft &
     &     (ista_kngp,iend_kngp,KMESH,kqitg,NMESH,IPRI,NFOUT,IL3 &
     &     ,RADR,WOS,QRSPS,gr_l,UNIVOL,mmpt,qitg_l &
     &     ,X,Y)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: qitgft
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

! Fouriner transformation of the deficit charge QRSPS.
!***********************************************************
  use m_Const_Parameters, only : DP, PAI4
  implicit none
  integer, intent(in) :: ista_kngp, iend_kngp, kmesh,kqitg,nmesh,ipri,nfout,il3
  real(kind=DP), intent(in), dimension(kmesh) :: radr, wos, qrsps
  real(kind=DP), intent(in)                   :: gr_l(ista_kngp:iend_kngp)
  real(kind=DP), intent(in)                   :: univol
  integer,       intent(in)                   :: mmpt
  real(kind=DP), intent(out)                  :: qitg_l(ista_kngp:iend_kngp,kqitg)
  real(kind=DP), dimension(kmesh), optional   :: x, y

  integer       ::  i, n
  real(kind=DP) ::  time_start, time_end
!!$#ifndef VPP
!!$  integer,parameter :: nblock = 32
!!$  integer           :: nn, kgpcount, nmeshcount
!!$  integer           :: ii,is,ie,nis,nie
!!$  real(kind=DP),pointer, dimension(:) :: fact1
!!$  real(kind=DP),allocatable,dimension(:,:) :: xm, ym
!!$  real(kind=DP) ::  pai4univol
!!$
!!$  allocate(fact1(nblock)); allocate(xm(nblock,nblock)); allocate(ym(nblock,nblock))
!!$#else
  real(kind=DP) :: gabs, f2
!!$#endif

  call gettod(time_start)

!CDIR PARALLEL DO PRIVATE ( gabs, x, y, n, f2 )
  do i = ista_kngp, iend_kngp  !for mpi
     gabs = gr_l(i)
     x(1:nmesh) = gabs*radr(1:nmesh)
     call dsjnv(il3,nmesh,x,y)
     qitg_l(i,mmpt) = 0.d0
     do n = 1, nmesh
#ifdef QITGFT_OLD
        f2 = wos(n)/univol*PAI4*qrsps(n)
#else
        f2 = wos(n)*qrsps(n)
#endif
        qitg_l(i,mmpt) = qitg_l(i,mmpt) + f2*y(n)
     enddo
#ifndef QITGFT_OLD
     qitg_l(i,mmpt) = qitg_l(i,mmpt)/univol*PAI4
#endif
  end do

  call gettod(time_end)
  if(ipri >= 2) write(nfout,'(" !bPP sub qitgft time = ",f20.6,"(sec.)")') &
       & (time_end-time_start)/1.d6
  IF(IPRI.GE.3) THEN
     WRITE(NFOUT,400)
400  FORMAT(' ',' SUB Q I T G F T ')
  end IF
end subroutine qitgft
  
subroutine qitgft_diff &
     &     (ista_kngp,iend_kngp,KMESH,kqitg,NMESH,IPRI,NFOUT,IL3 &
     &     ,RADR,WOS,QRSPS,gr_l,UNIVOL,mmpt &
     &     ,qitg_l,qitg_diff_l &
     &     ,X,Y,Z1,Z2)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: qitgft_diff
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

! Fouriner transformation of the deficit charge QRSPS.
!***********************************************************
  use m_Const_Parameters, only : DP, PAI4
  implicit none
  integer, intent(in) :: ista_kngp, iend_kngp, kmesh,kqitg,nmesh,ipri,nfout,il3
  real(kind=DP), intent(in), dimension(kmesh) :: radr, wos, qrsps
  real(kind=DP), intent(in), dimension(ista_kngp:iend_kngp) :: gr_l
  real(kind=DP), intent(in)                   :: univol
  integer, intent(in)                         :: mmpt
  real(kind=DP), intent(out), dimension(ista_kngp:iend_kngp,kqitg) &
       & :: qitg_l, qitg_diff_l
  real(kind=DP), dimension(kmesh), optional :: x, y, z1, z2

  integer       ::  i, n
  real(kind=DP) ::  gabs, fac
  real(kind=DP) ::  time_start, time_end
!!$#ifndef VPP
!!$  integer,parameter :: nblock = 32
!!$  integer           :: nn, kgpcount, nmeshcount
!!$  integer           :: ii,is,ie,nis,nie
!!$  real(kind=DP),pointer, dimension(:) :: fact1,fact2
!!$  real(kind=DP),allocatable,dimension(:,:) :: xm, ym, zm1, zm2
!!$  real(kind=DP) ::  pai4univol
!!$
!!$  allocate(fact1(nblock)); allocate(xm(nblock,nblock)); allocate(ym(nblock,nblock))
!!$  allocate(fact2(nblock)); allocate(zm1(nblock,nblock)); allocate(zm2(nblock,nblock))
!!$#endif

  call gettod(time_start)

!!$#ifdef VPP
!CDIR PARALLEL DO PRIVATE ( gabs, x, y, z1, z2, fac, n, f2 )
  do i = ista_kngp, iend_kngp  !for mpi
     gabs = gr_l(i)
     qitg_l(i,mmpt) = 0.d0
     qitg_diff_l(i,mmpt) = 0.d0
     x(1:nmesh) = gabs*radr(1:nmesh)
     call dsjnv(il3,  nmesh,x,y)
     if(il3 == 0) then
       z1 = 0.d0
     else
       call dsjnv(il3-1,nmesh,x,z1)
     endif
     call dsjnv(il3+1,nmesh,x,z2)
     do n = 1, nmesh
        fac = wos(n)/univol*PAI4*qrsps(n)
        qitg_l(i,mmpt) = qitg_l(i,mmpt) + fac*y(n)
        qitg_diff_l(i,mmpt) = qitg_diff_l(i,mmpt) + fac * radr(n) &
     &          / (2.d0*il3+1.d0) * (il3*z1(n)-(il3+1)*z2(n))
     enddo
  end do

  call gettod(time_end)
  if(ipri >= 2) write(nfout,'(" sub qitgft_diff time = ",f20.6,"(sec.)")') &
       & (time_end-time_start)/1.d6
  IF(IPRI >= 3 ) THEN
     WRITE(NFOUT,400)
400  FORMAT(' ',' SUB Q I T G F T _ D I F F')
     write(nfout,'('' mmpt ='',i5)') mmpt
     write(nfout,'(8f10.5)') (qitg_l(i,mmpt),i=ista_kngp+1,ista_kngp+20)
     write(nfout,'(8f10.5)') (qitg_diff_l(i,mmpt),i=ista_kngp+1,ista_kngp+20)
  end IF
end subroutine qitgft_diff

subroutine psbhs0(nfout, iatomn,ival,itpcc,alp,cc)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: psbhs0
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in) ::                      nfout, iatomn
  real(kind=DP), intent(out) ::               ival
  integer, intent(out) ::                     itpcc
  real(kind=DP), intent(out), dimension(2) :: alp, cc

  IF (IATOMN.EQ.1) THEN
!--*-- H  --*--
     IVAL   =   1.d0
     ITPCC  =   0
     ALP(1) =  16.2200D0
     ALP(2) =   5.5500D0
     CC(1)  =   1.1924D0
     CC(2)  =  -0.1924D0
  ELSE IF (IATOMN.EQ.2) THEN
!--*-- H  --*--
     ALP(1) =  56.2300D0
     ALP(2) =  19.2400D0
     CC(1)  =   1.1998D0
     CC(2)  =  -0.1998D0
  ELSE IF(IATOMN.EQ.3) THEN
!--*-- Li --*--
     IVAL   =   1.d0
     ITPCC  =   1
     ALP(1) =   1.8400D0
     ALP(2) =   0.7300D0
     CC(1)  =   2.9081D0
     CC(2)  =  -1.9081D0
  ELSE IF(IATOMN.EQ.4) THEN
!--*-- Be --*--
     ALP(1) =   2.6100D0
     ALP(2) =   1.0000D0
     CC(1)  =   1.5280D0
     CC(2)  =  -0.5280D0
  ELSE IF(IATOMN.EQ.5) THEN
!--*-- B  --*--
     IVAL   =   3.d0
     ITPCC  =   0
     ALP(1) =   6.2100D0
     ALP(2) =   2.4700D0
     CC(1)  =   1.6546D0
     CC(2)  =  -0.6546D0
  ELSE IF(IATOMN.EQ.6) THEN
!--*-- C  --*--
     IVAL   =   4.d0
     ITPCC  =   0
     ALP(1) =   9.2800D0
     ALP(2) =   3.6900D0
     CC(1)  =   1.5222D0
     CC(2)  =  -0.5222D0
  ELSE IF(IATOMN.EQ.7) THEN
!--*-- N  --*--
     IVAL   =   5.d0
     ITPCC  =   0
     ALP(1) =  12.8700D0
     ALP(2) =   5.1200D0
     CC(1)  =   1.4504D0
     CC(2)  =  -0.4504D0
  ELSE IF(IATOMN.EQ.8) THEN
!--*-- O  --*--
     IVAL   =   6.d0
     ITPCC  =   0
     ALP(1) =  18.0900D0
     ALP(2) =   7.1900D0
     CC(1)  =   1.4224D0
     CC(2)  =  -0.4224D0
  ELSE IF(IATOMN.EQ.9) THEN
!--*-- F  --*--
     IVAL   =   7.d0
     ITPCC  =   0
     ALP(1) =  23.7800D0
     ALP(2) =   9.4500D0
     CC(1)  =   1.3974D0
     CC(2)  =  -0.3974D0
  ELSE IF(IATOMN.EQ.10) THEN
!--*-- Ne  --*--
     ALP(1) =  29.1300D0
     ALP(2) =  11.5800D0
     CC(1)  =   1.3711D0
     CC(2)  =  -0.3711D0
  ELSE IF(IATOMN.EQ.11) THEN
!--*-- Na --*--
     IVAL   =   1.d0
     ITPCC  =   1
     ALP(1) =   1.7100D0
     ALP(2) =   0.5000D0
     CC(1)  =   5.1815D0
     CC(2)  =  -4.1815D0
  ELSE IF(IATOMN.EQ.12) THEN
!--*-- Mg --*--
     IVAL   =   2.d0
     ITPCC  =   1
     ALP(1) =   2.0400D0
     ALP(2) =   0.8100D0
     CC(1)  =   3.5602D0
     CC(2)  =  -2.5602D0
  ELSE IF(IATOMN.EQ.13) THEN
!--*-- Al --*--
     IVAL   =   3.d0
     ITPCC  =   0
     ALP(1) =   1.7700D0
     ALP(2) =   0.7000D0
     CC(1)  =   1.7905D0
     CC(2)  =  -0.7905D0
  ELSE IF(IATOMN.EQ.14) THEN
!--*-- Si --*--
     IVAL   =   4.d0
     ITPCC  =   0
     ALP(1) =   2.1600D0
     ALP(2) =   0.8600D0
     CC(1)  =   1.6054D0
     CC(2)  =  -0.6054D0
  ELSE IF(IATOMN.EQ.15) THEN
!--*-- P  --*--
     IVAL   =   5.d0
     ITPCC  =   0
     ALP(1) =   2.5900D0
     ALP(2) =   1.0300D0
     CC(1)  =   1.4995D0
     CC(2)  =  -0.4995D0
  ELSE IF(IATOMN.EQ.16) THEN
!--*-- S  --*--
     IVAL   =   6.d0
     ITPCC  =   0
     ALP(1) =   2.9900D0
     ALP(2) =   1.1900D0
     CC(1)  =   1.4261D0
     CC(2)  =  -0.4261D0
  ELSE IF(IATOMN.EQ.17) THEN
!--*-- Cl --*--
     IVAL   =   7.d0
     ITPCC  =   0
     ALP(1) =   3.4800D0
     ALP(2) =   1.3800D0
     CC(1)  =   1.3860D0
     CC(2)  =  -0.3860D0
  ELSE IF(IATOMN.EQ.18) THEN
!--*-- Ar --*--
     ALP(1) =   3.9900D0
     ALP(2) =   1.5900D0
     CC(1)  =   1.3622D0
     CC(2)  =  -0.3622D0
  ELSE IF(IATOMN.EQ.19) THEN
!--*-- K  --*--
     IVAL   =   1.d0
     ITPCC  =   1
     ALP(1) =   1.4200D0
     ALP(2) =   0.2600D0
     CC(1)  =   6.3140D0
     CC(2)  =  -5.3140D0
  ELSE IF(IATOMN.EQ.20) THEN
!--*-- Ca --*--
     IVAL   =   2.d0
     ITPCC  =   1
     ALP(1) =   1.6100D0
     ALP(2) =   0.4500D0
     CC(1)  =   4.8360D0
     CC(2)  =  -3.8360D0
  ELSE IF(IATOMN.EQ.21) THEN
!--*-- Sc --*--
     ALP(1) =   3.9600D0
     ALP(2) =   0.6900D0
     CC(1)  =   3.7703D0
     CC(2)  =  -2.7703D0
  ELSE IF(IATOMN.EQ.22) THEN
!--*-- Ti --*--
     ALP(1) =   4.6800D0
     ALP(2) =   0.9400D0
     CC(1)  =   3.3889D0
     CC(2)  =  -2.3889D0
  ELSE IF(IATOMN.EQ.23) THEN
!--*-- V  --*--
     ALP(1) =   5.1400D0
     ALP(2) =   1.1100D0
     CC(1)  =   2.9680D0
     CC(2)  =  -1.9680D0
  ELSE IF(IATOMN.EQ.24) THEN
!--*-- Cr --*--
     ALP(1) =   5.1900D0
     ALP(2) =   1.3700D0
     CC(1)  =   2.8897D0
     CC(2)  =  -1.8897D0
  ELSE IF(IATOMN.EQ.25) THEN
!--*-- Mn --*--
     ALP(1) =   6.0300D0
     ALP(2) =   1.6300D0
     CC(1)  =   2.7024D0
     CC(2)  =  -1.7024D0
  ELSE IF(IATOMN.EQ.26) THEN
!--*-- Fe --*--
     ALP(1) =   6.5100D0
     ALP(2) =   1.9100D0
     CC(1)  =   2.6179D0
     CC(2)  =  -1.6179D0
  ELSE IF(IATOMN.EQ.27) THEN
!--*-- Co --*--
     ALP(1) =   6.9500D0
     ALP(2) =   2.3800D0
     CC(1)  =   2.7407D0
     CC(2)  =  -1.7407D0
  ELSE IF(IATOMN.EQ.28) THEN
!--*-- Ni --*--
     ALP(1) =   7.6000D0
     ALP(2) =   2.7400D0
     CC(1)  =   2.6949D0
     CC(2)  =  -1.6949D0
  ELSE IF(IATOMN.EQ.29) THEN
!--*-- Cu --*--
     ALP(1) =   7.5900D0
     ALP(2) =   3.0200D0
     CC(1)  =   2.6959D0
     CC(2)  =  -1.6959D0
  ELSE IF(IATOMN.EQ.30) THEN
!--*-- Zn --*--
     ALP(1) =   8.7800D0
     ALP(2) =   3.4900D0
     CC(1)  =   2.6313D0
     CC(2)  =  -1.6313D0
  ELSE IF(IATOMN.EQ.31) THEN
!--*-- Ga --*--
     IVAL   =   3.d0
     ITPCC  =   0
     ALP(1) =   2.0100D0
     ALP(2) =   0.8000D0
     CC(1)  =   4.0433D0
     CC(2)  =  -3.0433D0
  ELSE IF(IATOMN.EQ.32) THEN
!--*-- Ge --*--
     IVAL   =   4.d0
     ITPCC  =   0
     ALP(1) =   2.2800D0
     ALP(2) =   0.9100D0
     CC(1)  =   3.1110D0
     CC(2)  =  -2.1110D0
  ELSE IF(IATOMN.EQ.33) THEN
!--*-- As --*--
     IVAL   =   5.d0
     ITPCC  =   0
     ALP(1) =   2.6000D0
     ALP(2) =   1.0300D0
     CC(1)  =   2.6218D0
     CC(2)  =  -1.6218D0
  ELSE IF(IATOMN.EQ.34) THEN
!--*-- Se --*--
     IVAL   =   6.d0
     ITPCC  =   0
     ALP(1) =   2.8800D0
     ALP(2) =   1.1400D0
     CC(1)  =   2.2934D0
     CC(2)  =  -1.2934D0
  ELSE IF(IATOMN.EQ.35) THEN
!--*-- Br --*--
     ALP(1) =   3.200D0
     ALP(2) =   1.390D0
     CC(1)  =   2.1007D0
     CC(2)  =  -1.1007D0
  ELSE IF(IATOMN.EQ.36) THEN
!--*-- Kr --*--
     ALP(1) =   3.4900D0
     ALP(2) =   1.3900D0
     CC(1)  =   1.9478D0
     CC(2)  =  -0.9478D0
  ELSE IF(IATOMN.EQ.37) THEN
!--*-- Rb --*--
     ALP(1) =   1.3700D0
     ALP(2) =   0.2100D0
     CC(1)  =   6.8301D0
     CC(2)  =  -5.8301D0
  ELSE IF(IATOMN.EQ.38) THEN
!--*-- Sr --*--
     ALP(1) =   1.5200D0
     ALP(2) =   0.3300D0
     CC(1)  =   4.8514D0
     CC(2)  =  -3.8514D0
  ELSE IF(IATOMN.EQ.39) THEN
!--*-- Y  --*--
     ALP(1) =   2.0600D0
     ALP(2) =   0.4900D0
     CC(1)  =   4.1719D0
     CC(2)  =  -3.1719D0
  ELSE IF(IATOMN.EQ.40) THEN
!--*-- Zr --*--
     ALP(1) =   2.2800D0
     ALP(2) =   0.6600D0
     CC(1)  =   3.9162D0
     CC(2)  =  -2.9162D0
  ELSE IF(IATOMN.EQ.41) THEN
!--*-- Nb --*--
     ALP(1) =   2.4100D0
     ALP(2) =   0.8200D0
     CC(1)  =   3.7419D0
     CC(2)  =  -2.7419D0
  ELSE IF(IATOMN.EQ.42) THEN
!--*-- Mo --*--
     ALP(1) =   2.5700D0
     ALP(2) =   1.0200D0
     CC(1)  =   3.8044D0
     CC(2)  =  -2.8044D0
  ELSE IF(IATOMN.EQ.43) THEN
!--*-- Tc --*--
     ALP(1) =   2.8200D0
     ALP(2) =   1.1200D0
     CC(1)  =   3.3669D0
     CC(2)  =  -2.3669D0
  ELSE IF(IATOMN.EQ.44) THEN
!--*-- Ru --*--
     ALP(1) =   3.0000D0
     ALP(2) =   1.1900D0
     CC(1)  =   3.0213D0
     CC(2)  =  -2.0213D0
  ELSE IF(IATOMN.EQ.45) THEN
!--*-- Rh --*--
     ALP(1) =   3.2100D0
     ALP(2) =   1.2800D0
     CC(1)  =   2.7857D0
     CC(2)  =  -1.7857D0
  ELSE IF(IATOMN.EQ.46) THEN
!--*-- Pd --*--
     ALP(1) =   3.3100D0
     ALP(2) =   1.3200D0
     CC(1)  =   2.5256D0
     CC(2)  =  -1.5256D0
  ELSE IF(IATOMN.EQ.47) THEN
!--*-- Ag --*--
     ALP(1) =   3.5300D0
     ALP(2) =   1.4100D0
     CC(1)  =   2.3857D0
     CC(2)  =  -1.3857D0
  ELSE IF(IATOMN.EQ.48) THEN
!--*-- Cd --*--
     ALP(1) =   3.9100D0
     ALP(2) =   1.5600D0
     CC(1)  =   2.3128D0
     CC(2)  =  -1.3128D0
  ELSE IF(IATOMN.EQ.49) THEN
!--*-- In --*--
     ALP(1) =   1.7900D0
     ALP(2) =   0.7100D0
     CC(1)  =   6.7251D0
     CC(2)  =  -5.7251D0
  ELSE IF(IATOMN.EQ.50) THEN
!--*-- Sn --*--
     ALP(1) =   1.9700D0
     ALP(2) =   0.7800D0
     CC(1)  =   5.0086D0
     CC(2)  =  -4.0086D0
  ELSE IF(IATOMN.EQ.51) THEN
!--*-- Sb --*--
     ALP(1) =   2.1200D0
     ALP(2) =   0.8500D0
     CC(1)  =   4.0534D0
     CC(2)  =  -3.0534D0
  ELSE IF(IATOMN.EQ.52) THEN
!--*-- Te --*--
     ALP(1) =   2.3700D0
     ALP(2) =   0.9500D0
     CC(1)  =   3.5696D0
     CC(2)  =  -2.5696D0
  ELSE IF(IATOMN.EQ.53) THEN
!--*-- I  --*--
     ALP(1) =   2.5200D0
     ALP(2) =   1.0100D0
     CC(1)  =   3.0856D0
     CC(2)  =  -2.0856D0
  ELSE IF(IATOMN.EQ.54) THEN
!--*-- Xe --*--
     ALP(1) =   2.6300D0
     ALP(2) =   1.0500D0
     CC(1)  =   2.6837D0
     CC(2)  =  -1.6837D0
  ELSE IF(IATOMN.EQ.55) THEN
!--*-- Cs --*--
     ALP(1) =   1.2900D0
     ALP(2) =   1.0500D0
     CC(1)  =   2.6837D0
     CC(2)  =  -1.6837D0
  ELSE IF(IATOMN.EQ.56) THEN
!--*-- Ba --*--
     ALP(1) =   1.3400D0
     ALP(2) =   0.2500D0
     CC(1)  =   5.3338D0
     CC(2)  =  -4.3338D0
  ELSE IF(IATOMN.EQ.57) THEN
!--*-- La --*--
     ALP(1) =   5.0000D0
     ALP(2) =   2.0000D0
     CC(1)  =   3.2494D0
     CC(2)  =  -2.2494D0
  ELSE IF(IATOMN.EQ.58) THEN
!--*-- Ce --*--
     ALP(1) =   5.3200D0
     ALP(2) =   2.1300D0
     CC(1)  =   3.0538D0
     CC(2)  =  -2.0538D0
  ELSE IF(IATOMN.EQ.59) THEN
!--*-- Pr --*--
     ALP(1) =   5.5300D0
     ALP(2) =   2.2100D0
     CC(1)  =   2.8422D0
     CC(2)  =  -1.8422D0
  ELSE IF(IATOMN.EQ.60) THEN
!--*-- Nd --*--
     ALP(1) =   5.6900D0
     ALP(2) =   2.2800D0
     CC(1)  =   2.6639D0
     CC(2)  =  -1.6639D0
  ELSE IF(IATOMN.EQ.61) THEN
!--*-- Pm --*--
     ALP(1) =   5.9100D0
     ALP(2) =   2.3600D0
     CC(1)  =   2.5202D0
     CC(2)  =  -1.5202D0
  ELSE IF(IATOMN.EQ.62) THEN
!--*-- Sm --*--
     ALP(1) =   6.0700D0
     ALP(2) =   2.4300D0
     CC(1)  =   2.3958D0
     CC(2)  =  -1.3958D0
  ELSE IF(IATOMN.EQ.63) THEN
!--*-- Eu --*--
     ALP(1) =   6.2400D0
     ALP(2) =   2.500D0
     CC(1)  =   2.2886D0
     CC(2)  =  -1.2886D0
  ELSE IF(IATOMN.EQ.64) THEN
!--*-- Gd --*--
     ALP(1) =   6.4700D0
     ALP(2) =   2.5900D0
     CC(1)  =   2.2076D0
     CC(2)  =  -1.2076D0
  ELSE IF(IATOMN.EQ.65) THEN
!--*-- Tb --*--
     ALP(1) =   7.1500D0
     ALP(2) =   2.8600D0
     CC(1)  =   2.2267D0
     CC(2)  =  -1.2267D0
  ELSE IF(IATOMN.EQ.66) THEN
!--*-- Dy --*--
     ALP(1) =   7.3400D0
     ALP(2) =   2.9400D0
     CC(1)  =   2.1510D0
     CC(2)  =  -1.1510D0
  ELSE IF(IATOMN.EQ.67) THEN
!--*-- Ho --*--
     ALP(1) =   7.6100D0
     ALP(2) =   3.0400D0
     CC(1)  =   2.0909D0
     CC(2)  =  -1.0909D0
  ELSE IF(IATOMN.EQ.68) THEN
!--*-- Er --*--
     ALP(1) =   7.8100D0
     ALP(2) =   3.1200D0
     CC(1)  =   2.0304D0
     CC(2)  =  -1.0304D0
  ELSE IF(IATOMN.EQ.69) THEN
!--*-- Tm --*--
     ALP(1) =   8.0000D0
     ALP(2) =   3.2000D0
     CC(1)  =   1.9772D0
     CC(2)  =  -0.9772D0
  ELSE IF(IATOMN.EQ.70) THEN
!--*-- Yb --*--
     ALP(1) =   8.2100D0
     ALP(2) =   3.2800D0
     CC(1)  =   1.9271D0
     CC(2)  =  -0.9271D0
  ELSE IF(IATOMN.EQ.71) THEN
!--*-- Lu --*--
     ALP(1) =   8.5000D0
     ALP(2) =   3.4000D0
     CC(1)  =   1.8919D0
     CC(2)  =  -0.8919D0
  ELSE IF(IATOMN.EQ.72) THEN
!--*-- Hf --*--
     ALP(1) =   2.0200D0
     ALP(2) =   0.8000D0
     CC(1)  =   5.8921D0
     CC(2)  =  -4.8921D0
  ELSE IF(IATOMN.EQ.73) THEN
!--*-- Ta --*--
     ALP(1) =   2.1600D0
     ALP(2) =   0.8600D0
     CC(1)  =   4.7264D0
     CC(2)  =  -3.7264D0
  ELSE IF(IATOMN.EQ.74) THEN
!--*-- W  --*--
     ALP(1) =   2.2600D0
     ALP(2) =   0.9000D0
     CC(1)  =   3.9028D0
     CC(2)  =  -2.9028D0
  ELSE IF(IATOMN.EQ.75) THEN
!--*-- Re --*--
     ALP(1) =   2.4400D0
     ALP(2) =   0.9800D0
     CC(1)  =   3.5075D0
     CC(2)  =  -2.5075D0
  ELSE IF(IATOMN.EQ.76) THEN
!--*-- Os --*--
     ALP(1) =   2.4700D0
     ALP(2) =   0.9900D0
     CC(1)  =   2.9996D0
     CC(2)  =  -1.9996D0
  ELSE IF(IATOMN.EQ.77) THEN
!--*-- Ir --*--
     ALP(1) =   2.6700D0
     ALP(2) =   1.0700D0
     CC(1)  =   2.8089D0
     CC(2)  =  -1.8089D0
  ELSE IF(IATOMN.EQ.78) THEN
!--*-- Pt --*--
     ALP(1) =   2.7100D0
     ALP(2) =   1.0800D0
     CC(1)  =   2.5166D0
     CC(2)  =  -1.5166D0
  ELSE IF(IATOMN.EQ.79) THEN
!--*-- Au --*--
     ALP(1) =   2.8500D0
     ALP(2) =   1.1400D0
     CC(1)  =   2.3778D0
     CC(2)  =  -1.3778D0
  ELSE IF(IATOMN.EQ.80) THEN
!--*-- Hg --*--
     ALP(1) =   2.9900D0
     ALP(2) =   1.2000D0
     CC(1)  =   2.2533D0
     CC(2)  =  -1.2533D0
  ELSE IF(IATOMN.EQ.81) THEN
!--*-- Tl --*--
     ALP(1) =   1.7700D0
     ALP(2) =   0.5900D0
     CC(1)  =   6.7158D0
     CC(2)  =  -5.7158D0
  ELSE IF(IATOMN.EQ.82) THEN
!--*-- Pb --*--
     ALP(1) =   1.9200D0
     ALP(2) =   0.7600D0
     CC(1)  =   6.5444D0
     CC(2)  =  -5.5444D0
  ELSE IF(IATOMN.EQ.83) THEN
!--*-- Bi --*--
     ALP(1) =   2.0300D0
     ALP(2) =   0.8100D0
     CC(1)  =   5.2104D0
     CC(2)  =  -4.2104D0
  ELSE IF(IATOMN.EQ.84) THEN
!--*-- Po --*--
     ALP(1) =   2.1700D0
     ALP(2) =   0.8700D0
     CC(1)  =   4.4190D0
     CC(2)  =  -3.4190D0
  ELSE IF(IATOMN.EQ.85) THEN
!--*-- At --*--
     ALP(1) =   2.2900D0
     ALP(2) =   0.9200D0
     CC(1)  =   3.8179D0
     CC(2)  =  -2.8179D0
  ELSE IF(IATOMN.EQ.86) THEN
!--*-- Rn --*--
     ALP(1) =   2.4500D0
     ALP(2) =   0.9800D0
     CC(1)  =   3.4235D0
     CC(2)  =  -2.4235D0
  ELSE IF(IATOMN.EQ.87) THEN
!--*-- Fr --*--
     ALP(1) =   1.2300D0
     ALP(2) =   0.1600D0
     CC(1)  =   8.4416D0
     CC(2)  =  -7.4416D0
  ELSE IF(IATOMN.EQ.88) THEN
!--*-- Ra --*--
     ALP(1) =   0.9900D0
     ALP(2) =   0.2400D0
     CC(1)  =   6.1114D0
     CC(2)  =  -5.1114D0
  ELSE IF(IATOMN.EQ.89) THEN
!--*-- Ac --*--
     ALP(1) =   1.4900D0
     ALP(2) =   0.3100D0
     CC(1)  =   4.7270D0
     CC(2)  =  -3.7270D0
  ELSE IF(IATOMN.EQ.90) THEN
!--*-- Th --*--
     ALP(1) =   1.5900D0
     ALP(2) =   0.4000D0
     CC(1)  =   4.3902D0
     CC(2)  =  -3.3902D0
  ELSE IF(IATOMN.EQ.91) THEN
!--*-- Pa --*--
     ALP(1) =   2.0600D0
     ALP(2) =   0.4500D0
     CC(1)  =   3.8349D0
     CC(2)  =  -2.8349D0
  ELSE IF(IATOMN.EQ.92) THEN
!--*-- U  --*--
     ALP(1) =   2.0800D0
     ALP(2) =   0.5000D0
     CC(1)  =   3.5183D0
     CC(2)  =  -2.5183D0
  ELSE IF(IATOMN.EQ.93) THEN
!--*-- Np  --*--
     ALP(1) =   2.7200D0
     ALP(2) =   0.6000D0
     CC(1)  =   3.5099D0
     CC(2)  =  -2.5099D0
  ELSE IF(IATOMN.EQ.94) THEN
!--*-- Pu  --*--
     ALP(1) =   2.7600D0
     ALP(2) =   0.6200D0
     CC(1)  =   3.1671D0
     CC(2)  =  -2.1671D0
  ELSE
     WRITE(NFOUT,*) ' IATOMN =',IATOMN,'IS NOT PREPARED IN PSBHS'
  END IF
end subroutine psbhs0

subroutine matrix_inversion(nfout,k,n,a,b)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: matrix_inversion
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)                        :: nfout,k, n
  real(kind=DP), intent(in),  dimension(k,k) :: a
  real(kind=DP), intent(out), dimension(k,k) :: b

  if(n == 1) then
     b(1,1) = 1.d0/a(1,1)
  else if(n == 2) then
     call inver2n(k,a,b)
  else if(n == 3) then
     call inver3n(k,a,b)
  else
     write(nfout,*) ' n >= 4 ', ' n = ',n
     !stop
     call phase_error_with_msg(nfout, 'n>=4',__LINE__,__FILE__)
  endif
end subroutine matrix_inversion

