!=======================================================================
!
!  SOFTWARE NAME : PHASE ($Revision: 618 $)
!
!  SUBROUINE:  ex_ggapw91, cr_ggapw91, ex_ggapbe, cr_ggapbe, xclda, ggabek, 
!            ggaprd, xcpotf_wigner, xcpotf_pzold, xcpotf_xalfa, xcpotf_pz,
!            xcpotf_vwn, xcpotf_mjw_bh_gl, cpden, gdiffs, gtgrad, wdcoef,
!            gcoef1, gcoef2, ggrade, cnggrd1, cpval, cpval_abs, getroh, 
!            getroh2, cnggrd, mkprdc, 
!
!  AUTHOR(S): T. Yamasaki, K. Kato   August/20/2003
!  
!  FURTHER MODIFICATION: M. Saito, T. Yamasaki, December/01/2003
!  FURTHER MODIFICATION: T. Yamasaki, June/03/2010
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
! $Id: b_XC_Potential.F90 618 2020-05-12 08:44:53Z ktagami $
!!====================================================================
!!        modified  by T.Yamasaki  94/11/14
!!
!!        Modified for Spin polrarized GGA by K. Kato 1995/1/18 
!!
! === xclda etc. are necessary to make 3D_Parallel, too!!! by tkato ============
!!BRANCH_P ORG_Parallel
! ==============================================================================
!$$#ifndef PARA3D
#ifdef __EDA__
subroutine ex_ggapw91(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc &
     &                    ,dFx_drho,dFx_dgradrho, exc_on_a_grid_wk, ist, ien )
#else
subroutine ex_ggapw91(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc &
     &                    ,dFx_drho,dFx_dgradrho, ist, ien )
#endif
!!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : PAI,DP
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif
  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in)        :: ist, ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)  :: grad_rho(ista_r:iend_r,nspin)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)  
  real(kind=DP),intent(out) :: exc
  real(kind=DP),intent(out) :: dFx_drho(ista_r:iend_r,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(out) :: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif


  integer :: is, i, iend_true
  real(kind=DP) :: exc1,facw,d,dd,fk,s,fac,s2,s3,s4&
       & ,p0,p1,p2,p3,p4,p5,p6,f,ex,fs,exd,exdd,exc0,excd,excdd

  real(kind=DP), parameter :: a1 = 0.19645d0
  real(kind=DP), parameter :: a2 = 0.27430d0
  real(kind=DP), parameter :: a3 = 0.15084d0
  real(kind=DP), parameter :: a4 = 100.0d0

! === KT_mod === 2014/11/25
!!!  real(kind=DP), parameter :: ax = -0.7385588d0
       real(kind=DP), parameter ::  ax = -0.7385587663820224d0
!=============== 2014/11/25

  real(kind=DP), parameter :: a  =  7.7956d0
  real(kind=DP), parameter :: b1 =  0.004d0

! === KT_mod === 2014/11/25
!  real(kind=DP), parameter :: thrd   = 0.33333333333d0
!  real(kind=DP), parameter :: thrd4  = 1.333333333333333d0
  real(kind=DP), parameter :: thrd   = 0.33333333333333333333333d0
  real(kind=DP), parameter :: thrd4  = 1.33333333333333333333333d0
!=============== 2014/11/25

  real(kind=DP), parameter :: thpith = 3.0936677262801d0

!  iend_true = iend_r
!  if(present(ien)) iend_true = min(iend_r,ien)

!---- Spin dependency

  facw = ispin
  exc = 0.d0
#ifdef __EDA__
  if(sw_eda==ON) then
! -----  ascat starts modifying  -----
  exc_on_a_grid_wk = 0.d0
! -----  ascat ceases modifying  -----
  endif
#endif
#ifdef NEC_TUNE_MXCP
!CDIR NOCONCUR
#endif
  do is = 1, ispin
     exc1 = 0.d0
#ifdef NEC_TUNE_MXCP
!CDIR INNER
#endif
!     do i = ista_r, iend_true       ! MPI

     do i = ist, ien       ! MPI
        d  = facw * chgrhr_l(i, is)
        dd = facw * grad_rho(i, is)
        fk = (3.d0*PAI*PAI*d)**thrd
        if(d > 1.d-05) then
           s = dd/(d*fk*2.d0)
        else
           s = 0.d0
        endif

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
        ex  = fac*f*d
        p5  = 2.d0 *(s*(a2-a3*p2)+a3*a4*s3*p2 - 2.d0*b1*s3)
        p6  = (a1*(p1+a*s*p0)+4.d0*b1*s3)*((a2-a3*p2)*s2-b1*s4)
        fs  = (p5*p3-p6*p3*p3)
        exd = thrd4*fac*(f-s*fs)
        exdd = ax*fs*0.5d0/thpith

        exc0 = ex / facw
        excd = exd
        dFx_drho(i, is) = excd 
        excdd = exdd
!!$c gradient of charge density 95/12/3 H.S.
        if(dabs(grad_rho(i, is)) > 1.d-9) then
           dFx_dgradrho(i, is) = excdd / grad_rho(i, is)
        else
           dFx_dgradrho(i, is) = 0.d0
        endif
        exc1 = exc1 + exc0*wos(i)
#ifdef __EDA__
        if(sw_eda==ON) then
! -----  ascat starts modifying  -----
        exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + exc0!*f2or1(i)
! -----  ascat ceases modifying  -----
        endif
#endif
     end do
     exc = exc + exc1
  end do
end subroutine ex_ggapw91


!!===================================================================
!!        modified  by T.Yamasaki  94/11/15
!!
!!        The original is modified and Spin polrarized GGA programs
!!        are added in this subroutine
!!                                                      by K. Kato 1995/1/18 
!!
#ifdef __EDA__
subroutine cr_ggapw91(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_trho,wos,exc,dF_drho,exc_on_a_grid_wk,ist,ien)
#else
subroutine cr_ggapw91(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_trho,wos,exc,dF_drho,ist,ien)
#endif
!!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : PAI,DP
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif
  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in)        :: ist, ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(inout) :: grad_trho(ista_r:iend_r)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)  
  real(kind=DP),intent(inout) :: exc
  real(kind=DP),intent(inout) :: dF_drho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout) :: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif

  real(kind=DP), parameter :: zeta_minimum    = 1.d-4
  real(kind=DP), parameter :: zeta_minimum2   = 1.d-6
  real(kind=DP), parameter :: density_minimum = 1.d-20
  real(kind=DP), parameter :: density_minimum2= 1.d-5

  real(kind=DP), parameter :: a  = 0.0310907d0,     a1   = 0.21370d0
  real(kind=DP), parameter :: b1 = 7.5957d0,        b2n  = 3.5876d0
  real(kind=DP), parameter :: b3 = 1.6382d0,        b4   = 0.49294d0
  real(kind=DP), parameter :: p  = 1.00d0,          p1   = p+1.d0
  real(kind=DP), parameter :: ap  = 0.015545d0,     a1p  = 0.20548d0
  real(kind=DP), parameter :: b1p = 14.1189d0,      b2np = 6.1977d0
  real(kind=DP), parameter :: b3p = 3.3662d0,       b4p = 0.62517d0
! vwn 95/12/2 Y.M
  real(kind=DP), parameter :: aq  = 0.016887d0,     a1q = 0.11125d0
  real(kind=DP), parameter :: b1q = 10.3570d0,      b2nq = 3.6231d0
  real(kind=DP), parameter :: b3q = 0.88026d0,      b4q = 0.49671d0
! vwn 95/12/2 Y.M

  real(kind=DP), parameter :: xnu = 15.75592d0,     cc0 =  0.004235d0
  real(kind=DP), parameter :: cx  = -0.001667212d0, alf =  0.09d0
  real(kind=DP), parameter :: c1  =  0.002568d0,    c2  =  0.023266d0
  real(kind=DP), parameter :: c3  =  7.389d-6,      c4  =  8.723d0
  real(kind=DP), parameter :: c5  =  0.472d0,       c6  =  7.389d-2
  real(kind=DP), parameter :: a4  = 100.d0

! ==== KT_mod ==== 2014/11/25
!  real(kind=DP), parameter :: thrd   = 0.333333333333d0
!  real(kind=DP), parameter :: sixth7 = 1.1666666666666666d0
  real(kind=DP), parameter :: thrd   = 0.3333333333333333333333d0
  real(kind=DP), parameter :: sixth7 = 1.1666666666666666666666d0
! ================= 2014/11/25

  integer       :: is,i,iend_true
  real(kind=DP) :: facw,bet,delt,exc1,g,g3,g4,facpon,d,dd,rs,fk,sk,t,s&
       & , q0,rs12,rs32,rsp,q1,q2,eu,q3,eurs,ecd,pon,b,b2,t2,s2,t4,t6 &
       & , rs2,rs3,q4,q5,q6,q7,cc,r0,r1,coeff,r2,r3,a4ms,h0,h1,h,q8 &
       & , h0t,h0b,h0rs,h1t,ccrs,r1rs,h1rs,ht,hrs &
       & , ec1,ec1d,ec1dd,exc0,excd,excdd &
       & , thrd2,thrd4,zeta,q0p,q1p,q2p,q3p,q0q,q1q,q2q,q3q,fzd,fzdd,fz&
       & , euzt,zetadxd,bzt,h0zt,h1zt,gzt &
       & , onpzeta, onmzeta, onzt4 

! ==== kT_add ===== 2014/11/25
  real(kind=DP), parameter  :: eta = 1.d-12
  real(kind=DP), parameter  :: sixthm = -0.16666666666666666666d0
! ================= 2014/11/25

  facw = ispin

! === KT_mod === 2014/11/25
!  bet = xnu*cc0
  bet=0.06672455060314922d0
! ============== 2014/11/25

  delt = 2.d0*alf/bet

!  iend_true = iend_r
!  if(present(ien)) iend_true = min(iend_r,ien)

!---- Spin dependency

#ifdef NEC_TUNE_MXCP
!CDIR NOCONCUR
#endif
  do is=1, ispin
     exc1 = 0.d0
     if ( ispin == 1 ) then
        g = 1.d0
        g3  = g**3
        g4  = g3*g
        facpon = -delt/(g3*bet)

#ifdef NEC_TUNE_MXCP
!CDIR INNER
#endif
!        do i = ista_r, iend_true  ! MPI

        do i = ist, ien  ! MPI
           d = facw*chgrhr_l(i, 1)
           if(d < density_minimum) cycle
! gradient of charge density 95/12/3 H.S.
           dd =facw* grad_trho(i)
           rs = (0.75d0/(PAI*d))**thrd
           fk = (3.d0*PAI*PAI*d)**thrd
           sk = dsqrt(4.d0*fk/PAI)
           if(d > density_minimum2 ) then
              t = dd/(d*sk*2.d0)
              s = dd/(d*fk*2.d0)
           else
              t = 0.d0
              s = 0.d0
           endif

           q0 = -2.d0 *a*(1.d0+a1*rs)
           rs12 = dsqrt(rs)
           rs32 = rs12**3
           rsp = rs**p
           q1 = 2.d0*a*(b1*rs12+b2n*rs+b3*rs32+b4*rs*rsp)
           q2 = log(1.d0 + 1.d0/q1)
           eu = q0*q2
           q3 = a*(b1/rs12 +2.d0*b2n +3.d0*b3*rs12 +2.d0*b4*p1*rsp)
           eurs = -2.d0 *a*a1*q2 - q0*q3/(q1**2 + q1)
           ecd = eu-thrd*rs*eurs
           pon = facpon*eu
           b   = delt/(dexp(pon)-1.d0)
           b2  = b*b
           t2  = t*t
           s2  = s*s
           t4  = t2*t2
           t6  = t4*t2
           rs2 = rs*rs
           rs3 = rs2*rs
           q4  = 1.d0+b*t2
           q5  = 1.d0+b*t2+b2*t4
           q6  = c1+c2*rs+c3*rs2
           q7  = 1.d0+c4*rs+c5*rs2+c6*rs3
           cc  = -cx + q6/q7
           r0  = (sk/fk)**2
           r1  = a4*r0*g4
           coeff = cc-cc0-3.d0*cx/7.d0
           r2  = xnu*coeff*g3
           a4ms = -a4*g4*s2
           r3  = dexp(a4ms)
           h0  = g3*(bet/delt)*dlog(1.d0+delt*q4*t2/q5)
           h1  = r3*r2*t2
           h   = h0+h1
!======================================================
           q8  = q5*q5+delt*q4*q5*t2
           h0t = 2.d0*bet*t*( 1.d0 +2.d0*b*t2 )/q8 * g3
           h0b = -bet*t6*( 2.d0*b +b2*t2 )/q8 * g3
           h0rs = h0b*b*eurs*(b+delt)/ bet / g3
           h1t  = 2.d0*r3*r2*t*(1.d0-r1*t2)
           ccrs = ( c2 +2.d0*c3*rs )/q7 -q6*( c4 +2*c5*rs +3.d0*c6*rs2 )/q7**2
           r1rs = 100.d0 *r0/rs
           h1rs = xnu*t2*r3*(ccrs - coeff*t2*r1rs) * g3
           ht   = h0t +h1t
           hrs  = h0rs +h1rs
           ec1  = d*h/facw
           ec1d = h -thrd*rs*hrs -sixth7*t*ht
           ec1dd = 0.5d0*ht/sk

!---------------------------------------------------
           exc0 = eu*d/facw + ec1
           excd = ecd + ec1d
           dF_drho(i, is) = dF_drho(i, is) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(dabs(grad_trho(i)) > 1.0d-9) then
              grad_trho(i) = excdd / grad_trho(i)
           else
              grad_trho(i) = 0.d0
           endif
!!$           exec1 = eu*d/facw
           exc1 = exc1 + exc0*wos(i)
#ifdef __EDA__
           if(sw_eda==ON) then
! -----  ascat starts modifying  -----
           exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + exc0!*f2or1(i)
! -----  ascat ceases modifying  -----
           endif
#endif
        end do

     else if ( ispin  == 2 ) then
        thrd2 = thrd * 2.0d0
        thrd4 = thrd * 4.0d0
        fzdd = 8.d0 /(9.d0*(2.d0**thrd4-2.d0))

        do i = ista_r, iend_r      ! MPI
           d = chgrhr_l(i, 1) + chgrhr_l(i, ispin)
           if(d.lt.density_minimum) cycle
           zeta = ( chgrhr_l(i, 1) - chgrhr_l(i, ispin) ) / d
           onpzeta = 2.d0*chgrhr_l(i,1)     / d
           onmzeta = 2.d0*chgrhr_l(i,ispin) / d
           g  = ( onpzeta**thrd2 + onmzeta**thrd2) * 0.5d0
!!$           g  = ( (1+zeta)**thrd2 + (1-zeta)**thrd2 ) * 0.5

           if(dabs(zeta) < zeta_minimum) then
              fz = 4.d0/9.d0 * zeta*zeta*(1.d0 + 5.d0/54.d0 * zeta*zeta&
                   &                 *(1.d0 +44.d0/133.d0* zeta*zeta ))&
                   &    /(2.d0**thrd4-2.d0)
           else
              fz = ( (1.d0+zeta)**thrd4+(1.d0-zeta)**thrd4-2.d0)/(2.d0**thrd4-2.d0)
           end if

           g3  = g**3
           g4  = g3*g
           facpon = -delt/(g3*bet)

           if(d < density_minimum) cycle
! gradient of charge density 95/12/3 H.S.
           dd = grad_trho(i)
!!$c        dd = grad_rho(i, 1) + grad_rho(i, 2)
! gradient of charge density 95/12/3 H.S.
           rs = (0.75d0/(PAI*d))**thrd
           fk = (3.d0*PAI*PAI*d)**thrd
           sk = dsqrt(4.d0*fk/PAI)

           if(d > density_minimum2) then
! --> check if this is correct or not
!!$              t = dd/(d*sk*2.d0)
              t = dd/(d*sk*2.d0)/g  ! 2011/04/26 revised according to a Katsumata-san's comment.
! <--
              s = dd/(d*fk*2.d0)
           else
              t = 0.d0
              s = 0.d0
           endif

           q0 = -2.d0 *a*(1.d0+a1*rs)
           rs12 = dsqrt(rs)
           rs32 = rs12**3
           rsp = rs**p
           q1 = 2.d0*a*(b1*rs12+b2n*rs+b3*rs32+b4*rs*rsp)
           q2 = log(1.d0+1.d0/q1)
           q3 = a*( b1/rs12 +2.d0*b2n +3.d0*b3*rs12 +2.d0*b4*p1*rsp )

           q0p = -2.d0*ap*(1.d0+a1p*rs)
           q1p =  2.d0*ap*(b1p*rs12+b2np*rs+b3p*rs32+b4p*rs*rsp)
           q2p = log(1.d0+1.d0/q1p)
           q3p = ap*(b1p/rs12+2.d0*b2np+3.d0*b3p*rs12+2.d0*b4p*p1*rsp)
! vwn 95/12/2 Y.M
           q0q = -2.d0 *aq*(1.d0+a1q*rs)
           q1q =  2.d0 *aq*(b1q*rs12+b2nq*rs+b3q*rs32+b4q*rs*rsp)
           q2q = log(1.d0 +1.d0/q1q)
           q3q = aq*(b1q/rs12+2.d0*b2nq+3.d0*b3q*rs12+2.d0*b4q*p1*rsp)

           if(dabs(zeta) < zeta_minimum) then
              fzd = 8.d0/(9.d0*(2.d0**thrd4 - 2.d0 )) &
                   &     *zeta *( 1.d0 + 5.d0/27.d0 *zeta*zeta &
                   &           *( 1.d0 +22.d0/45.d0 *zeta*zeta ))
           else
              fzd = (4.d0/3.d0)/(2.d0**thrd4-2.d0)*(onpzeta**thrd-onmzeta**thrd )
           end if

           onzt4 = onmzeta*onpzeta * (1.d0 + zeta*zeta)
           eu = q0*q2 + ( q0p * q2p - q0 * q2 ) * fz * zeta**4 &
                &        -   q0q * q2q * fz / fzdd * onzt4
           eurs = (-2.d0*a*a1*q2 - q0*q3/(q1**2+q1)) * (1.d0-fz*zeta**4) &
                &+(-2.d0*ap*a1p*q2p-q0p*q3p/(q1p**2+q1p))*fz*zeta**4 &
                &-(-2.d0*aq*a1q*q2q-q0q*q3q/(q1q**2+q1q))*fz/fzdd*onzt4
           euzt = -q0q*q2q/fzdd*(fzd*(1.d0-zeta**4)-4.d0*fz*zeta**3) &
                &+(q0p*q2p-q0*q2)*(fzd*zeta +4.d0*fz )*zeta**3
! vwn 95/12/2 Y.M
           zetadxd = - 2.d0 * chgrhr_l(i, 3-is) * (-1.d0)**is / d
           ecd = eu - thrd*rs*eurs + euzt * zetadxd
!
! ==== KT_mod ==== 2014/11/25
!           if(onpzeta < zeta_minimum2 .or. onmzeta < zeta_minimum2) then
!              ec1 = 0.d0; ec1d = 0.d0; ec1dd= 0.d0
!              goto 10
!           end if
! ================= 2014/11/25

           pon = facpon*eu
           b   = delt/(dexp(pon)-1.d0)
           b2  = b*b
           t2  = t*t
           s2  = s*s
           t4  = t2*t2
           t6  = t4*t2
           rs2 = rs*rs
           rs3 = rs2*rs
           q4  = 1.d0+b*t2
           q5  = 1.d0+b*t2+b2*t4
           q6  = c1+c2*rs+c3*rs2
           q7  = 1.d0+c4*rs+c5*rs2+c6*rs3
           cc  = -cx + q6/q7
           r0  = (sk/fk)**2
           r1  = a4*r0*g4
           coeff = cc-cc0-3.d0*cx/7.d0
           r2  = xnu*coeff*g3
! --> Check if this is correct or not
!!$           a4ms = -a4*g4*s2
           a4ms = -r1*t2   ! 2011/04/26 revised according to a Katsumata-san's comment.
! <--
           r3  = dexp(a4ms)
           h0  = g3*(bet/delt)*dlog(1.d0+delt*q4*t2/q5)
           h1  = r3*r2*t2
           h   = h0+h1

           q8  = q5*q5+delt*q4*q5*t2
           h0t = 2*bet*t*(1.d0+2*b*t2)/q8 * g3
           h0b = -bet*t6*(2*b+b2*t2)/q8 * g3
! --> Check if this is correct or not
!!$           h0rs = h0b * b2 * eurs / bet / g3
           h0rs = h0b * b * eurs * (b+delt) / bet / g3   ! 2011/04/26 revised according to a Katsumata-san's comment.
! <--
           h1t  = 2*r3*r2*t*(1.d0-r1*t2)
           ccrs = (c2+2*c3*rs)/q7-q6*(c4+2*c5*rs+3*c6*rs2)/q7**2
! --> Check if this is correct or not
!!$           r1rs = 100*r0/rs
           r1rs = r1/rs  ! 2011/04/26 revised according to a Katsumata-san's comment.
! <--
           h1rs = xnu*t2*r3*(ccrs - coeff*t2*r1rs) * g3
           ht   = h0t+h1t
           hrs  = h0rs + h1rs

           if(dabs(zeta) < zeta_minimum) then
              gzt = -2.d0/9.d0*zeta *( 1.d0 +14.d0/27.d0 *zeta*zeta &
                   &                *( 1.d0 +13.d0/18.d0 *zeta*zeta ))
           else
! ==== KT_mod ==== 2014/11/25
!!!              gzt = ( 1.d0/onpzeta**thrd -1.d0/onmzeta**thrd )/3.d0
	      gzt = ( (onpzeta**2+eta)**sixthm - (onmzeta**2+eta)**sixthm)/3.d0
! ================ 2014/11/25
           end if

! --> Check if this is correct or not
!!$           bzt = pon*b**2 * (euzt-3*eu*gzt/g) / ( bet * g3 ) 
           bzt = b * (b+delt) * (euzt-3*eu*gzt/g) / ( bet * g3 )  ! 2011/04/26 revised according to a Katsumata-san's comment.
! <--
           h0zt = 3*gzt*h0/g + h0b*bzt
! --> Check if this is correct or not
!!$           h1zt =(3*gzt + 4*gzt*a4ms ) * h1 / g
           h1zt =(3.d0 + 4.d0*a4ms ) * h1 * gzt/ g   ! 2011/04/26 revised according to a Katsumata-san's comment.
! <--
           
           ec1  = d*h * 0.5d0

! --> Check if this is correct or not
!!$           ec1d = h-thrd*rs*hrs-sixth7*t*ht + (h0zt+h1zt) * zetadxd
           ec1d = h-thrd*rs*hrs-sixth7*t*ht + (h0zt+h1zt-gzt*ht/g) * zetadxd    ! 2011/04/26 revised according to a Katsumata-san's comment.
!!$           ec1dd = 0.5*ht/sk
           ec1dd = 0.5d0*ht/sk/g     ! 2011/04/26 revised according to a Katsumata-san's comment.
! <--

10         continue
           exc0 = eu*d * 0.5d0 + ec1
           excd = ecd + ec1d
           dF_drho(i, is) = dF_drho(i, is) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(is == 2) then
              if(dabs(grad_trho(i)) > 1.0d-9) then
                 grad_trho(i) = excdd / grad_trho(i)
              else
                 grad_trho(i) = 0.d0
              endif
           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d * 0.5d0
           exc1 = exc1 + exc0*wos(i)
#ifdef __EDA__
           if(sw_eda==ON) then
! -----  ascat starts modifying  -----
           exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + exc0!*f2or1(i)
! -----  ascat ceases modifying  -----
           endif
#endif
        end do
     end if
     exc = exc + exc1
  end do
end subroutine cr_ggapw91

#ifndef PREV_EX_GGAPBE
!!$====================================================================
!!$        modified  by T.Yamasaki  94/11/14
!!$
!!$        Modified for Spin polrarized GGA by K. Kato 1995/1/18 
!!$
!!$        Modified to Perdew-Burke-Ernzerhof(PBE) GGA by K. Kato 1999/4/1 
!!$
#ifdef __EDA__
subroutine ex_ggapbe(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dFx_drho,dFx_dgradrho,exc_on_a_grid_wk,revPBE,ist,ien)
#else
subroutine ex_ggapbe(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dFx_drho,dFx_dgradrho,revPBE,ist,ien)
#endif
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,PAI
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif
  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in)        :: ist, ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)  :: grad_rho(ista_r:iend_r,nspin)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)  
  real(kind=DP),intent(out) :: exc
  real(kind=DP),intent(out) :: dFx_drho(ista_r:iend_r,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(out) :: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif
  logical, intent(in),optional :: revPBE

!!  real(kind=DP), parameter :: ax = -0.7385588d0
  real(kind=DP), parameter :: ax = -0.7385587663820224d0 !! ax = -(3/4)*(3/PI)^(1/3)
  real(kind=DP), parameter :: thrd = 0.33333333333d0, thrd4 = 1.333333333333333d0
  real(kind=DP), parameter :: thpith = 3.0936677262801d0
!!  real(kind=DP), parameter :: cupa = 0.804, yum = 0.21951
!!  real(kind=DP), parameter :: cupa = 0.8040d0, yum = 0.2195149727645171d0
  real(kind=DP), parameter ::  yum = 0.2195149727645171d0
  real(kind=DP) :: cupa

  real(kind=DP) :: facw,d,dd,fk,s,fac,s2,yums2,f,ex,fs,exd,exdd,exc0,excd,excdd,exc1
  integer       :: is,i,iend_true
  logical       :: repbe

  repbe = .false.
  if(present(revPBE)) repbe = revPBE
  cupa = 0.8040d0
  if(repbe) cupa = 1.245d0

!  iend_true = iend_r
!  if(present(ien)) iend_true = min(iend_r,ien)

!---- Spin dependency

  facw = ispin
  exc  = 0.d0
#ifdef __EDA__
  if(sw_eda==ON) then
! -----  ascat starts modifying  -----
  exc_on_a_grid_wk = 0.d0
! -----  ascat ceases modifying  -----
  endif
#endif
#ifdef NEC_TUNE_MXCP
!CDIR NOCONCUR
#endif
  do is = 1, ispin
     exc1 = 0.d0
#ifdef NEC_TUNE_MXCP
!CDIR INNER
#endif
!     do i = ista_r, iend_true  ! MPI

     do i = ist, ien  ! MPI
        d  = facw * chgrhr_l(i, is)
        dd = facw * grad_rho(i, is)
        fk = (3.d0*PAI*PAI*d)**thrd
!!$        sk = dsqrt(4*fk/PAI)
        if(d > 1.d-05) then
           s = dd/(d*fk*2.d0)
        else
           s = 0.d0
        endif
!-------------------------------------
        fac = ax*d**thrd
        s2  = s*s
        yums2 = yum * s2
        f = 1.d0 + cupa * yums2 / ( cupa + yums2 )
        ex  = fac*f*d
        fs = 2.d0 * cupa ** 2 * yum * s / ( cupa + yums2 ) ** 2
        exd = thrd4*fac*(f-s*fs)
        exdd = ax*fs*0.5d0/thpith
!------------------------------------------     
        exc0 = ex / facw
        excd = exd
        dFx_drho(i, is) = excd 
        excdd = exdd
! gradient of charge density 95/12/3 H.S.
        if(dabs(grad_rho(i, is)) > 1.d-9) then
           dFx_dgradrho(i, is) = excdd / grad_rho(i, is)
        else
           dFx_dgradrho(i, is) = 0.d0
        endif
! gradient of charge density 95/12/3 H.S.
        exc1 = exc1 + exc0*wos(i)
#ifdef __EDA__
        if(sw_eda==ON) then
! -----  ascat starts modifying  -----
        exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + exc0!*f2or1(i)
! -----  ascat ceases modifying  -----
        endif
#endif
     end do
     exc = exc + exc1
  end do
end subroutine ex_ggapbe
#else
!!$====================================================================
!!$        modified  by T.Yamasaki  94/11/14
!!$
!!$        Modified for Spin polrarized GGA by K. Kato 1995/1/18 
!!$
!!$        Modified to Perdew-Burke-Ernzerhof(PBE) GGA by K. Kato 1999/4/1 
!!$
! ==== DEBUG! ==================================================================
#if 1
#ifdef __EDA__
subroutine ex_ggapbe(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dFx_drho,dFx_dgradrho,exc_on_a_grid_wk,ien)
#else
subroutine ex_ggapbe(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dFx_drho,dFx_dgradrho,ien)
#endif
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,PAI
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif
  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in),optional :: ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)  :: grad_rho(ista_r:iend_r,nspin)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)  
  real(kind=DP),intent(out) :: exc
  real(kind=DP),intent(out) :: dFx_drho(ista_r:iend_r,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(out) :: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif

#ifndef PREV_EX_GGAPBE
  real(kind=DP), parameter :: ax = -0.7385587663820224d0 !! ax = -(3/4)*(3/PI)^(1/3)
#else
  real(kind=DP), parameter :: ax = -0.7385588d0
#endif
  real(kind=DP), parameter :: thrd = 0.33333333333d0, thrd4 = 1.333333333333333d0
  real(kind=DP), parameter :: thpith = 3.0936677262801d0
#ifndef PREV_EX_GGAPBE
  real(kind=DP), parameter :: cupa = 0.8040d0, yum = 0.2195149727645171d0
#else
  real(kind=DP), parameter :: cupa = 0.804, yum = 0.21951
#endif

  real(kind=DP) :: facw,d,dd,fk,s,fac,s2,yums2,f,ex,fs,exd,exdd,exc0,excd,excdd,exc1
  integer       :: is,i,iend_true

  iend_true = iend_r
  if(present(ien)) iend_true = min(iend_r,ien)
!---- Spin dependency

  facw = ispin
  exc  = 0.d0
#ifdef __EDA__
  exc_on_a_grid_wk = 0.d0
#endif
#ifdef NEC_TUNE_MXCP
!CDIR NOCONCUR
#endif
  do is = 1, ispin
     exc1 = 0.d0
#ifdef NEC_TUNE_MXCP
!CDIR INNER
#endif
     do i = ista_r, iend_true
        d  = facw * chgrhr_l(i, is)
        dd = facw * grad_rho(i, is)
        fk = (3*PAI*PAI*d)**thrd
!!$        sk = dsqrt(4*fk/PAI)
        if(d > 1.d-05) then
           s = dd/(d*fk*2)
        else
           s = 0.d0
        endif
!-------------------------------------
        fac = ax*d**thrd
        s2  = s*s
        yums2 = yum * s2
        f = 1 + cupa * yums2 / ( cupa + yums2 )
        ex  = fac*f*d
        fs = 2 * cupa ** 2 * yum * s / ( cupa + yums2 ) ** 2
        exd = thrd4*fac*(f-s*fs)
        exdd = ax*fs*0.5d0/thpith
!------------------------------------------     
        exc0 = ex / facw
        excd = exd
        dFx_drho(i, is) = excd 
        excdd = exdd
! gradient of charge density 95/12/3 H.S.
        if(dabs(grad_rho(i, is)) > 1.d-9) then
           dFx_dgradrho(i, is) = excdd / grad_rho(i, is)
        else
           dFx_dgradrho(i, is) = 0.d0
        endif
! gradient of charge density 95/12/3 H.S.
        exc1 = exc1 + exc0*wos(i)
#ifdef __EDA__
        if(sw_eda==ON) then
! -----  ascat starts modifying  -----
        exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + exc0!*f2or1(i)
! -----  ascat ceases modifying  -----
        endif
#endif
     end do
     exc = exc + exc1
  end do
end subroutine ex_ggapbe
#else
#ifdef __EDA__
subroutine ex_ggapbe(nspin,ispin,chgrhr_l,grad_rho,f2or1,exc,dFx_drho,dFx_dgradrho,exc_on_a_grid_wk)
#else
subroutine ex_ggapbe(nspin,ispin,chgrhr_l,grad_rho,f2or1,exc,dFx_drho,dFx_dgradrho)
#endif
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,PAI
  use m_Parallelization,   only : ista_fftph, iend_fftph
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif
  implicit none

  integer,intent(in)        :: nspin,ispin
  real(kind=DP),intent(in)  :: chgrhr_l(ista_fftph:iend_fftph,ispin)
  real(kind=DP),intent(in)  :: grad_rho(ista_fftph:iend_fftph,nspin)
  real(kind=DP),intent(in)  :: f2or1(ista_fftph:iend_fftph)  
  real(kind=DP),intent(out) :: exc
  real(kind=DP),intent(out) :: dFx_drho(ista_fftph:iend_fftph,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(ista_fftph:iend_fftph,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(out) :: exc_on_a_grid_wk(ista_fftph:iend_fftph)
! -----  ascat ceases modifying  -----
#endif

#ifndef PREV_EX_GGAPBE
  real(kind=DP), parameter :: ax = -0.7385587663820224d0 !! ax = -(3/4)*(3/PI)^(1/3)
#else
  real(kind=DP), parameter :: ax = -0.7385588d0
#endif
  real(kind=DP), parameter :: thrd = 0.33333333333d0, thrd4 = 1.333333333333333d0
  real(kind=DP), parameter :: thpith = 3.0936677262801d0
#ifndef PREV_EX_GGAPBE
  real(kind=DP), parameter :: cupa = 0.8040d0, yum = 0.2195149727645171d0
#else
  real(kind=DP), parameter :: cupa = 0.804, yum = 0.21951
#endif

  real(kind=DP) :: facw,d,dd,fk,s,fac,s2,yums2,f,ex,fs,exd,exdd,exc0,excd,excdd,exc1
  integer       :: is,i

!---- Spin dependency

  facw = ispin
  exc  = 0.d0
#ifdef __EDA__
  if(sw_eda==ON) then
! -----  ascat starts modifying  -----
  exc_on_a_grid_wk = 0.d0
! -----  ascat ceases modifying  -----
  endif
#endif
#ifdef NEC_TUNE_MXCP
!CDIR NOCONCUR
#endif
  do is = 1, ispin
     exc1 = 0.d0
#ifdef NEC_TUNE_MXCP
!CDIR INNER
#endif
     do i = ista_fftph, iend_fftph
        d  = facw * chgrhr_l(i, is)
        dd = facw * grad_rho(i, is)
        fk = (3*PAI*PAI*d)**thrd
!!$        sk = dsqrt(4*fk/PAI)
        if(d > 1.d-05) then
           s = dd/(d*fk*2)
        else
           s = 0.d0
        endif
!-------------------------------------
        fac = ax*d**thrd
        s2  = s*s
        yums2 = yum * s2
        f = 1 + cupa * yums2 / ( cupa + yums2 )
        ex  = fac*f*d
        fs = 2 * cupa ** 2 * yum * s / ( cupa + yums2 ) ** 2
        exd = thrd4*fac*(f-s*fs)
        exdd = ax*fs*0.5d0/thpith
!------------------------------------------     
        exc0 = ex / facw
        excd = exd
        dFx_drho(i, is) = excd 
        excdd = exdd
! gradient of charge density 95/12/3 H.S.
        if(dabs(grad_rho(i, is)) > 1.d-9) then
           dFx_dgradrho(i, is) = excdd / grad_rho(i, is)
        else
           dFx_dgradrho(i, is) = 0.d0
        endif
! gradient of charge density 95/12/3 H.S.
        exc1 = exc1 + exc0*f2or1(i)
#ifdef __EDA__
        if(sw_eda==ON) then
! -----  ascat starts modifying  -----
        exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + exc0!*f2or1(i)
! -----  ascat ceases modifying  -----
        endif
#endif
     end do
     exc = exc + exc1
  end do
end subroutine ex_ggapbe
#endif
! ==== DEBUG! ==================================================================
#endif

#ifndef PREV_CR_GGAPBE
!!$===================================================================
!!$        modified  by T.Yamasaki  94/11/15
!!$
!!$        The original is modified and Spin polrarized GGA programs
!!$        are added in this subroutine
!!$                                                      by K. Kato 1995/1/18 
!!$
!!$        Modified to Perdew-Burke-Ernzerhof(PBE) GGA by K. Kato 1999/4/1 
!!$
#ifdef __EDA__
subroutine cr_ggapbe(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_trho,wos,exc,dF_drho,exc_on_a_grid_wk,ecor,ist,ien)
#else
subroutine cr_ggapbe(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_trho,wos,exc,dF_drho,ecor,ist,ien)
#endif
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : PAI,DP
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif
  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in)        :: ist, ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(inout) :: grad_trho(ista_r:iend_r)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)  
  real(kind=DP),intent(inout) :: exc
  real(kind=DP),intent(inout) :: dF_drho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout) :: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif
  real(kind=DP),intent(out),optional :: ecor
  real(kind=DP), parameter :: zeta_minimum    = 1.d-4
  real(kind=DP), parameter :: zeta_minimum2   = 1.d-6
  real(kind=DP), parameter :: density_minimum = 1.d-20
  real(kind=DP), parameter :: density_minimum2= 1.d-5

  real(kind=DP), parameter  :: a  = 0.0310907d0,    a1   = 0.21370d0
  real(kind=DP), parameter  :: b1 = 7.5957d0,       b2n  = 3.5876d0
  real(kind=DP), parameter  :: b3 = 1.6382d0,       b4   = 0.49294d0
  real(kind=DP), parameter  :: p  = 1.00d0,         p1   = p+1.d0
!Fix!  real(kind=DP), parameter  :: ap  = 0.015545d0,    a1p  = 0.20548d0
  real(kind=DP), parameter  :: ap  = 0.01554535d0,    a1p  = 0.20548d0
  real(kind=DP), parameter  :: b1p = 14.1189d0,     b2np = 6.1977d0
  real(kind=DP), parameter  :: b3p = 3.3662d0,      b4p  = 0.62517d0
! vwn 95/12/2 Y.M
!Fix!  real(kind=DP), parameter  :: aq  = 0.016887d0,    a1q  = 0.11125d0
  real(kind=DP), parameter  :: aq  = 0.0168869d0,    a1q  = 0.11125d0
  real(kind=DP), parameter  :: b1q = 10.3570d0,     b2nq = 3.6231d0
  real(kind=DP), parameter  :: b3q = 0.88026d0,     b4q  = 0.49671d0
! vwn 95/12/2 Y.M
  !Fix!real(kind=DP), parameter  :: gamma = 0.031091
  real(kind=DP), parameter  :: gamma = 0.03109069086965489503494086371273d0
!--------------------------
!Fix!  real(kind=DP), parameter  :: xnu = 15.75592d0, cc0 = 0.004235d0
  real(kind=DP), parameter  :: bet=0.06672455060314922d0
  real(kind=DP), parameter  :: delt=bet/gamma
  real(kind=DP), parameter  :: thrd = 0.333333333333d0, sixth7 = 1.1666666666666666d0

!Fix! T.Yamatomto 2011/06/29
  real(kind=DP), parameter  :: eta = 1.d-12
  real(kind=DP), parameter  :: sixthm = -0.16666666666666d0
!Fix! T.Yamatomto 2011/06/29


  integer       :: is,i,iend_true
  real(kind=DP) :: facw,exc1,g,g3,facpon,d,dd,rs,fk,sk,t,q0 &
       &         , rs12,rs32,rsp,q1,q2,eu,q3,eurs,ecd,pon &
       &         , b,b2,t2,t4,t6,q4,q5,h0,h,q8,h0t,h0b,h0rs,ht,hrs &
       &         , ec1,ec1d,ec1dd,exc0,excd,excdd &
       &         , thrd2,thrd4,fzdd,zeta,onpzeta,onmzeta,fz,zetadxd &
       &         , q0p,q1p,q2p,q3p,q0q,q1q,q2q,q3q,fzd,onzt4,gzt,bzt,h0zt,euzt

  facw = ispin

!  iend_true = iend_r
!  if(present(ien)) iend_true = min(iend_r,ien)
!Fix!  bet = xnu*cc0
!Fix!  delt = bet / gamma
!---- Spin dependency
#ifdef NEC_TUNE_MXCP
!CDIR NOCONCUR
#endif
  if(present(ecor)) ecor = 0.d0
  do is = 1, ispin
     exc1 = 0.d0
     if ( ispin == 1 ) then
        g   = 1.d0
        g3  = g**3
!!$        g4  = g3*g
        facpon = -delt/(g3*bet)

#ifdef NEC_TUNE_MXCP
!CDIR INNER
#endif
!        do i = ista_r,iend_true ! MPI

        do i = ist, ien     ! MPI
           d  = facw*chgrhr_l(i, 1)
           if(d < density_minimum) cycle
! gradient of charge density 95/12/3 H.S.
           dd = facw* grad_trho(i)
           !Fix!rs = (0.75/(PAI*d))**thrd
           rs = (0.75d0/(PAI*d))**thrd
           fk = (3.d0*PAI*PAI*d)**thrd
           sk = dsqrt(4.d0*fk/PAI)
           if(d > density_minimum2) then
              t = dd/(d*sk*2.d0)
!!$              s = dd/(d*fk*2)
           else
              t = 0.d0
!!$              s = 0.d0
           endif

           q0 = -2.d0 *a*(1.d0 + a1*rs)
           rs12 = dsqrt(rs)
           rs32 = rs12**3
           rsp = rs**p
           q1 = 2.d0 *a*(b1*rs12+b2n*rs+b3*rs32+b4*rs*rsp)
           q2 = log( 1.d0 + 1.d0/q1)
           eu = q0*q2
           q3 = a*(b1/rs12 +2.d0*b2n +3.d0*b3*rs12 +2.d0*b4*p1*rsp)
           eurs = -2.d0 *a*a1*q2 - q0*q3/(q1**2 + q1)
           ecd = eu-thrd*rs*eurs
           pon = facpon*eu
           b   = delt/(dexp(pon)-1.d0)
           b2  = b*b
           t2  = t*t
!!$           s2  = s*s
           t4  = t2*t2
           t6  = t4*t2
           q4  = 1.d0 + b*t2
           q5  = 1.d0 + b*t2+b2*t4
           h0  = g3*(bet/delt)*dlog(1d0+delt*q4*t2/q5)
           h   = h0
           q8  = q5*q5+delt*q4*q5*t2
           h0t = 2.d0 *bet*t*( 1.d0 +2d0 *b*t2)/q8 * g3
           h0b = -bet*t6*( 2.d0 *b+b2*t2)/q8 * g3
           h0rs = h0b*b*eurs*(b+delt)/ bet / g3
           ht   = h0t
           hrs  = h0rs
           ec1  = d*h/facw
           ec1d = h-thrd*rs*hrs-sixth7*t*ht
           !Fix!ec1dd = 0.5*ht/sk
           ec1dd = 0.5d0*ht/sk

           exc0 = eu*d/facw + ec1
           excd = ecd + ec1d
           dF_drho(i, is) = dF_drho(i, is) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(dabs(grad_trho(i)) > 1.0d-9) then
              grad_trho(i) = excdd / grad_trho(i)
           else
              grad_trho(i) = 0.d0
           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d/facw
           exc1 = exc1 + exc0*wos(i)
#ifdef __EDA__
           if(sw_eda==ON) then
! -----  ascat starts modifying  -----
           exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + exc0!*f2or1(i)
! -----  ascat ceases modifying  -----
           endif
#endif
        end do

     else if ( ispin ==  2 ) then
        thrd2 = thrd * 2.d0
        thrd4 = thrd * 4.d0
        fzdd = 8.d0 /(9.d0*(2.d0**thrd4 - 2.d0 ))

        do i = ista_r,iend_r
           d = chgrhr_l(i, 1) + chgrhr_l(i, ispin)
           if(d < density_minimum ) cycle
           zeta = ( chgrhr_l(i, 1) - chgrhr_l(i, 2) ) / d
           onpzeta = 2.d0 *chgrhr_l(i, 1) / d
           onmzeta = 2.d0 *chgrhr_l(i, 2) / d
           g = ( onpzeta**thrd2 + onmzeta**thrd2 )/2.d0
           if(dabs(zeta) < zeta_minimum ) then
              fz= 4.d0/9.d0 *zeta*zeta*(1.d0 + 5.d0/54.d0 *zeta*zeta &
                   &               *(1.d0 +44.d0/135.d0*zeta*zeta )) &
                   &       /(2.d0**thrd4 -2.d0)
           else
              fz = (onpzeta**thrd4 + onmzeta**thrd4 -2.d0) / (2.d0**thrd4 - 2.d0)
           end if
           g3  = g**3
!!$           g4  = g3*g
           facpon = -delt/(g3*bet)
           if(d < density_minimum) cycle
           dd = grad_trho(i)
           rs = (0.75d0/(PAI*d))**thrd
           fk = (3.d0*PAI*PAI*d)**thrd
           sk = dsqrt(4.d0*fk/PAI)

           if(d > density_minimum2) then
!Fix!              t = dd/(d*sk*2)
              t = dd/(d*sk*2.d0*g)
!!$              s = dd/(d*fk*2)
           else
              t = 0.d0
!!$              s = 0.d0
           endif

           q0 = -2.d0 *a*(1.d0+a1*rs)
           rs12 = dsqrt(rs)
           rs32 = rs12**3
           rsp = rs**p
           q1 = 2.d0 *a*(b1*rs12+b2n*rs+b3*rs32+b4*rs*rsp)
           q2 = log( 1.d0 +1.d0/q1 )
           q3 = a*(b1/rs12 +2.d0*b2n +3.d0*b3*rs12 +2.d0*b4*p1*rsp)

           q0p = -2.d0 *ap*(1.d0+a1p*rs)
           q1p = 2.d0*ap*(b1p*rs12+b2np*rs+b3p*rs32+b4p*rs*rsp)
           q2p = log(1.d0 +1.d0/q1p)
           q3p = ap*(b1p/rs12 +2.d0*b2np +3.d0*b3p*rs12 +2.d0*b4p*p1*rsp)
! vwn 95/12/2 Y.M
           q0q = -2.d0*aq*(1.d0+a1q*rs)
           q1q =  2.d0*aq*(b1q*rs12+b2nq*rs+b3q*rs32+b4q*rs*rsp)
           q2q = log(1.d0 +1.d0/q1q)
           q3q = aq*(b1q/rs12 +2.d0*b2nq +3.d0*b3q*rs12 +2.d0*b4q*p1*rsp)
           if(dabs(zeta) < zeta_minimum) then
              fzd = 8.d0/9.d0 /( 2.d0**thrd4 -2.d0 ) &
                   &     *zeta *( 1.d0 + 5.d0/27.d0 *zeta*zeta &
                   &           *( 1.d0 +22.d0/45.d0 *zeta*zeta ))
           else
              fzd = 4.d0/3.d0 /( 2.d0**thrd4 -2.d0) &
                   &     *( onpzeta**thrd -onmzeta**thrd )
           end if
           onzt4= onmzeta*onpzeta *( 1.d0 + zeta*zeta )
           eu = q0*q2 + ( q0p * q2p - q0 * q2 ) * fz * zeta**4 &
                &             -   q0q * q2q * fz / fzdd * onzt4
           eurs = ( -2.d0*a*a1*q2 -q0*q3/(q1**2 +q1 )) * ( 1.d0 -fz*zeta**4 )&
                &+( -2.d0*ap*a1p*q2p -q0p*q3p/( q1p**2 +q1p )) *fz*zeta**4 &
                &-( -2.d0*aq*a1q*q2q -q0q*q3q/( q1q**2 +q1q )) *fz/fzdd*onzt4
           euzt =  -q0q*q2q /fzdd *( fzd*onzt4 -4.d0*fz*zeta**3 )&
                &        +( q0p*q2p -q0*q2 )*(fzd*zeta +4.d0*fz )*zeta**3 
! vwn 95/12/2 Y.M

           zetadxd = - 2.d0 * chgrhr_l(i, 3-is) * (-1.d0)**is / d
           ecd = eu-thrd*rs*eurs + euzt * zetadxd 

!Fix!           if(onpzeta < zeta_minimum2 .or. onmzeta < zeta_minimum2) then
!Fix!              ec1  = 0.d0;    ec1d = 0.d0;     ec1dd= 0.d0
!Fix!              go to 10
!Fix!           end if
           pon = facpon*eu
           b   = delt/(dexp(pon)-1.d0)
           b2  = b*b
           t2  = t*t
!!$           s2  = s*s
           t4  = t2*t2
           t6  = t4*t2
           q4  = 1.d0 +b*t2
           q5  = 1.d0 +b*t2+b2*t4
           h0  = g3*(bet/delt)*dlog(1.d0+delt*q4*t2/q5)
           h   = h0

           q8  = q5*q5+delt*q4*q5*t2
           h0t = 2.d0 *bet*t*(1.d0 +2.d0*b*t2)/q8 * g3
           h0b = -bet*t6*( 2.d0*b +b2*t2)/q8 * g3
           h0rs = h0b * b * eurs * (b+delt) / bet / g3
           ht   = h0t
           hrs  = h0rs

           if(dabs(zeta) < zeta_minimum) then
              gzt = -2.d0/9.d0*zeta *( 1.d0 +14.d0/27.d0 *zeta*zeta &
                   &             *( 1.d0 +13.d0/18.d0 *zeta*zeta ))
           else
!Fix!              gzt = ( 1.d0/onpzeta**thrd -1.d0/onmzeta**thrd )/3.d0
              gzt = ( (onpzeta**2+eta)**sixthm - (onmzeta**2+eta)**sixthm)/3.d0
           end if
           bzt = b * ( b+delt) * ( euzt - 3.d0*eu*gzt/g ) / ( bet * g3 ) 
           h0zt = 3.d0 * gzt * h0 / g + h0b * bzt

           !Fix!ec1  = d*h * 0.5
           ec1  = d*h * 0.5d0
! --> revised by Kato-san,  1999Jul16
!!$           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + h0zt * zetadxd
           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxd
           !Fix!ec1dd = 0.5*ht/sk / g
           ec1dd = 0.5d0*ht/sk / g

10         continue
           !FiX!EXC0= eu*d * 0.5 + ec1
           exc0 = eu*d * 0.5d0 + ec1
           excd = ecd + ec1d
           dF_drho(i, is) = dF_drho(i, is) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(is ==  2) then
              if(dabs(grad_trho(i)) > 1.0d-9) then
                 grad_trho(i) = excdd / grad_trho(i)
              else
                 grad_trho(i) = 0.d0
              endif
           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d * 0.5
           exc1 = exc1 + exc0*wos(i)
#ifdef __EDA__
           if(sw_eda==ON) then
! -----  ascat starts modifying  -----
           exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + exc0!*f2or1(i)
! -----  ascat ceases modifying  -----
           endif
#endif
        end do
     end if
     exc = exc + exc1
     if(present(ecor)) ecor = ecor + exc1
  end do
end subroutine cr_ggapbe
#else
!!$===================================================================
!!$        modified  by T.Yamasaki  94/11/15
!!$
!!$        The original is modified and Spin polrarized GGA programs
!!$        are added in this subroutine
!!$                                                      by K. Kato 1995/1/18 
!!$
!!$        Modified to Perdew-Burke-Ernzerhof(PBE) GGA by K. Kato 1999/4/1 
!!$
! ==== DEBUG! ==================================================================
#if 1
#ifdef __EDA__
subroutine cr_ggapbe(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_trho,wos,exc,dF_drho,exc_on_a_grid_wk,ien)
#else
subroutine cr_ggapbe(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_trho,wos,exc,dF_drho,ien)
#endif
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : PAI,DP
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif
  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in),optional :: ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(inout) :: grad_trho(ista_r:iend_r)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)  
  real(kind=DP),intent(inout) :: exc
  real(kind=DP),intent(inout) :: dF_drho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout) :: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif

  real(kind=DP), parameter :: zeta_minimum    = 1.d-4
  real(kind=DP), parameter :: zeta_minimum2   = 1.d-6
  real(kind=DP), parameter :: density_minimum = 1.d-20
  real(kind=DP), parameter :: density_minimum2= 1.d-5

  real(kind=DP), parameter  :: a  = 0.0310907d0,    a1   = 0.21370d0
  real(kind=DP), parameter  :: b1 = 7.5957d0,       b2n  = 3.5876d0
  real(kind=DP), parameter  :: b3 = 1.6382d0,       b4   = 0.49294d0
  real(kind=DP), parameter  :: p  = 1.00d0,         p1   = p+1.d0
#ifndef PREV_CR_GGAPBE
  real(kind=DP), parameter  :: ap  = 0.01554535d0,    a1p  = 0.20548d0
#else
  real(kind=DP), parameter  :: ap  = 0.015545d0,    a1p  = 0.20548d0
#endif
  real(kind=DP), parameter  :: b1p = 14.1189d0,     b2np = 6.1977d0
  real(kind=DP), parameter  :: b3p = 3.3662d0,      b4p  = 0.62517d0
! vwn 95/12/2 Y.M
#ifndef PREV_CR_GGAPBE
  real(kind=DP), parameter  :: aq  = 0.0168869d0,    a1q  = 0.11125d0
#else
  real(kind=DP), parameter  :: aq  = 0.016887d0,    a1q  = 0.11125d0
#endif
  real(kind=DP), parameter  :: b1q = 10.3570d0,     b2nq = 3.6231d0
  real(kind=DP), parameter  :: b3q = 0.88026d0,     b4q  = 0.49671d0
! vwn 95/12/2 Y.M
#ifndef PREV_CR_GGAPBE
  real(kind=DP), parameter  :: gamma = 0.03109069086965489503494086371273d0
#else
  real(kind=DP), parameter  :: gamma = 0.031091
#endif
!--------------------------
#ifndef PREV_CR_GGAPBE
  real(kind=DP), parameter  :: bet=0.06672455060314922d0
  real(kind=DP), parameter  :: delt=bet/gamma
#else
  real(kind=DP), parameter  :: xnu = 15.75592d0, cc0 = 0.004235d0
#endif
  real(kind=DP), parameter  :: thrd = 0.333333333333d0, sixth7 = 1.1666666666666666d0

#ifndef PREV_CR_GGAPBE
!Fix! T.Yamatomto 2011/06/29
  real(kind=DP), parameter  :: eta = 1.d-12
  real(kind=DP), parameter  :: sixthm = -0.16666666666666d0
!Fix! T.Yamatomto 2011/06/29
#endif

  integer       :: is,i,iend_true
#ifndef PREV_CR_GGAPBE
  real(kind=DP) :: facw,exc1,g,g3,facpon,d,dd,rs,fk,sk,t,q0 &
#else
  real(kind=DP) :: facw,bet,delt,exc1,g,g3,facpon,d,dd,rs,fk,sk,t,q0 &
#endif
       &         , rs12,rs32,rsp,q1,q2,eu,q3,eurs,ecd,pon &
       &         , b,b2,t2,t4,t6,q4,q5,h0,h,q8,h0t,h0b,h0rs,ht,hrs &
       &         , ec1,ec1d,ec1dd,exc0,excd,excdd &
       &         , thrd2,thrd4,fzdd,zeta,onpzeta,onmzeta,fz,zetadxd &
       &         , q0p,q1p,q2p,q3p,q0q,q1q,q2q,q3q,fzd,onzt4,gzt,bzt,h0zt,euzt

  facw = ispin

  iend_true = iend_r
  if(present(ien)) iend_true = min(iend_r,ien)

#ifdef PREV_CR_GGAPBE
  bet = xnu*cc0
  delt = bet / gamma
#endif
!---- Spin dependency
#ifdef NEC_TUNE_MXCP
!CDIR NOCONCUR
#endif
  do is = 1, ispin
     exc1 = 0.d0
     if ( ispin == 1 ) then
        g   = 1.d0
        g3  = g**3
!!$        g4  = g3*g
        facpon = -delt/(g3*bet)

#ifdef NEC_TUNE_MXCP
!CDIR INNER
#endif
        do i = ista_r,iend_true ! MPI
           d  = facw*chgrhr_l(i, 1)
           if(d < density_minimum) cycle
! gradient of charge density 95/12/3 H.S.
           dd = facw* grad_trho(i)
#ifndef PREV_CR_GGAPBE
           rs = (0.75d0/(PAI*d))**thrd
#else
           rs = (0.75/(PAI*d))**thrd
#endif
           fk = (3*PAI*PAI*d)**thrd
           sk = dsqrt(4*fk/PAI)
           if(d > density_minimum2) then
              t = dd/(d*sk*2)
!!$              s = dd/(d*fk*2)
           else
              t = 0.d0
!!$              s = 0.d0
           endif

           q0 = -2*a*(1 + a1*rs)
           rs12 = dsqrt(rs)
           rs32 = rs12**3
           rsp = rs**p
           q1 = 2*a*(b1*rs12+b2n*rs+b3*rs32+b4*rs*rsp)
           q2 = log(1 + 1.d0/q1)
           eu = q0*q2
           q3 = a*(b1/rs12+2*b2n+3*b3*rs12+2*b4*p1*rsp)
           eurs = -2*a*a1*q2 - q0*q3/(q1**2 + q1)
           ecd = eu-thrd*rs*eurs
           pon = facpon*eu
           b   = delt/(dexp(pon)-1)
           b2  = b*b
           t2  = t*t
!!$           s2  = s*s
           t4  = t2*t2
           t6  = t4*t2
           q4  = 1 + b*t2
           q5  = 1 + b*t2+b2*t4
           h0  = g3*(bet/delt)*dlog(1+delt*q4*t2/q5)
           h   = h0
           q8  = q5*q5+delt*q4*q5*t2
           h0t = 2*bet*t*(1+2*b*t2)/q8 * g3
           h0b = -bet*t6*(2*b+b2*t2)/q8 * g3
           h0rs = h0b*b*eurs*(b+delt)/ bet / g3
           ht   = h0t
           hrs  = h0rs
           ec1  = d*h/facw
           ec1d = h-thrd*rs*hrs-sixth7*t*ht
#ifndef PREV_CR_GGAPBE
           ec1dd = 0.5d0*ht/sk
#else
           ec1dd = 0.5*ht/sk
#endif
           exc0 = eu*d/facw + ec1
           excd = ecd + ec1d
           dF_drho(i, is) = dF_drho(i, is) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(dabs(grad_trho(i)) > 1.0d-9) then
              grad_trho(i) = excdd / grad_trho(i)
           else
              grad_trho(i) = 0.d0
           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d/facw
           exc1 = exc1 + exc0*wos(i)
#ifdef __EDA__
           if(sw_eda==ON) then
! -----  ascat starts modifying  -----
           exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + exc0!*f2or1(i)
! -----  ascat ceases modifying  -----
           endif
#endif
        end do
     else if ( ispin ==  2 ) then
        thrd2 = thrd * 2
        thrd4 = thrd * 4
        fzdd = 8/(9*(2**thrd4 - 2 ))

        do i = ista_r,iend_r
           d = chgrhr_l(i, 1) + chgrhr_l(i, ispin)
           if(d < density_minimum ) cycle
           zeta = ( chgrhr_l(i, 1) - chgrhr_l(i, 2) ) / d
           onpzeta = 2*chgrhr_l(i, 1) / d
           onmzeta = 2*chgrhr_l(i, 2) / d
           g = ( onpzeta**thrd2 + onmzeta**thrd2 )/2
           if(dabs(zeta) < zeta_minimum ) then
              fz= 4.d0/9 *zeta*zeta*(1 + 5.d0/54 *zeta*zeta &
                   &               *(1 +44.d0/135*zeta*zeta )) &
                   &       /(2**thrd4 -2)
           else
              fz = (onpzeta**thrd4 + onmzeta**thrd4 -2) / (2**thrd4 - 2)
           end if
           g3  = g**3
!!$           g4  = g3*g
           facpon = -delt/(g3*bet)
           if(d < density_minimum) cycle
           dd = grad_trho(i)
           rs = (0.75d0/(PAI*d))**thrd
           fk = (3*PAI*PAI*d)**thrd
           sk = dsqrt(4*fk/PAI)

           if(d > density_minimum2) then
#ifndef PREV_CR_GGAPBE
              t = dd/(d*sk*2*g)
#else
              t = dd/(d*sk*2)
#endif
           else
              t = 0.d0
           endif

           q0 = -2*a*(1+a1*rs)
           rs12 = dsqrt(rs)
           rs32 = rs12**3
           rsp = rs**p
           q1 = 2*a*(b1*rs12+b2n*rs+b3*rs32+b4*rs*rsp)
           q2 = log(1+1/q1)
           q3 = a*(b1/rs12+2*b2n+3*b3*rs12+2*b4*p1*rsp)

           q0p = -2*ap*(1+a1p*rs)
           q1p = 2*ap*(b1p*rs12+b2np*rs+b3p*rs32+b4p*rs*rsp)
           q2p = log(1+1/q1p)
           q3p = ap*(b1p/rs12+2*b2np+3*b3p*rs12+2*b4p*p1*rsp)
! vwn 95/12/2 Y.M
           q0q = -2*aq*(1+a1q*rs)
           q1q =  2*aq*(b1q*rs12+b2nq*rs+b3q*rs32+b4q*rs*rsp)
           q2q = log(1+1/q1q)
           q3q = aq*(b1q/rs12+2*b2nq+3*b3q*rs12+2*b4q*p1*rsp)
           if(dabs(zeta) < zeta_minimum) then
              fzd = 8.d0/9.d0 /( 2.d0**thrd4 -2 ) &
                   &     *zeta *( 1 + 5.d0/27 *zeta*zeta &
                   &           *( 1 +22.d0/45 *zeta*zeta ))
           else
              fzd = 4.d0/3.d0 /( 2.d0**thrd4 -2) &
                   &     *( onpzeta**thrd -onmzeta**thrd )
           end if
           onzt4= onmzeta*onpzeta *( 1 + zeta*zeta )
           eu = q0*q2 + ( q0p * q2p - q0 * q2 ) * fz * zeta**4 &
                &             -   q0q * q2q * fz / fzdd * onzt4
           eurs = ( -2*a*a1*q2 -q0*q3/(q1**2 +q1 )) * ( 1 -fz*zeta**4 )&
                &+( -2*ap*a1p*q2p -q0p*q3p/( q1p**2 +q1p )) *fz*zeta**4 &
                &-( -2*aq*a1q*q2q -q0q*q3q/( q1q**2 +q1q )) *fz/fzdd*onzt4
           euzt =  -q0q*q2q /fzdd *( fzd*onzt4 -4*fz*zeta**3 )&
                &        +( q0p*q2p -q0*q2 )*(fzd*zeta +4*fz )*zeta**3 
! vwn 95/12/2 Y.M

           zetadxd = - 2 * chgrhr_l(i, 3-is) * (-1.d0)**is / d
           ecd = eu-thrd*rs*eurs + euzt * zetadxd 

#ifdef NEC_TUNE_MXCP
           if(onpzeta < zeta_minimum2 .or. onmzeta < zeta_minimum2) then
              ec1  = 0.d0;    ec1d = 0.d0;     ec1dd= 0.d0
              go to 10
           end if
#endif
           pon = facpon*eu
           b   = delt/(dexp(pon)-1)
           b2  = b*b
           t2  = t*t
!!$           s2  = s*s
           t4  = t2*t2
           t6  = t4*t2
           q4  = 1+b*t2
           q5  = 1+b*t2+b2*t4
           h0  = g3*(bet/delt)*dlog(1+delt*q4*t2/q5)
           h   = h0

           q8  = q5*q5+delt*q4*q5*t2
           h0t = 2*bet*t*(1+2*b*t2)/q8 * g3
           h0b = -bet*t6*(2*b+b2*t2)/q8 * g3
           h0rs = h0b * b * eurs * (b+delt) / bet / g3
           ht   = h0t
           hrs  = h0rs

           if(dabs(zeta) < zeta_minimum) then
              gzt = -2.d0/9*zeta *( 1 +14.d0/27 *zeta*zeta &
                   &             *( 1 +13.d0/18 *zeta*zeta ))
           else
#ifndef PREV_CR_GGAPBE
              gzt = ( (onpzeta**2+eta)**sixthm - (onmzeta**2+eta)**sixthm)/3.d0
#else              
              gzt = ( 1.d0/onpzeta**thrd -1.d0/onmzeta**thrd )/3.d0
#endif
           end if
           bzt = b * ( b+delt) * ( euzt - 3*eu*gzt/g ) / ( bet * g3 ) 
           h0zt = 3 * gzt * h0 / g + h0b * bzt

#ifndef PREV_CR_GGAPBE
           ec1  = d*h * 0.5d0
#else
           ec1  = d*h * 0.5
#endif
! --> revised by Kato-san,  1999Jul16
!!$           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + h0zt * zetadxd
           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxd
           ec1dd = 0.5*ht/sk / g

10         continue
#ifndef PREV_CR_GGAPBE
           exc0 = eu*d * 0.5d0 + ec1
#else
           exc0 = eu*d * 0.5 + ec1
#endif
           excd = ecd + ec1d
           dF_drho(i, is) = dF_drho(i, is) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(is ==  2) then
              if(dabs(grad_trho(i)) > 1.0d-9) then
                 grad_trho(i) = excdd / grad_trho(i)
              else
                 grad_trho(i) = 0.d0
              endif
           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d * 0.5
           exc1 = exc1 + exc0*wos(i)
#ifdef __EDA__
           if(sw_eda==ON) then
! -----  ascat starts modifying  -----
           exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + exc0!*f2or1(i)
! -----  ascat ceases modifying  -----
           endif
#endif
        end do
     end if
     exc = exc + exc1
  end do
end subroutine cr_ggapbe
#else
#ifdef __EDA__
subroutine cr_ggapbe(nspin,ispin,chgrhr_l,grad_trho,f2or1,exc,dF_drho,exc_on_a_grid_wk)
#else
subroutine cr_ggapbe(nspin,ispin,chgrhr_l,grad_trho,f2or1,exc,dF_drho)
#endif
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : PAI,DP
  use m_Parallelization,   only : ista_fftph, iend_fftph
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif
  implicit none

  integer,intent(in)        :: nspin,ispin
  real(kind=DP),intent(in)  :: chgrhr_l(ista_fftph:iend_fftph,ispin)
  real(kind=DP),intent(inout) :: grad_trho(ista_fftph:iend_fftph)
  real(kind=DP),intent(in)  :: f2or1(ista_fftph:iend_fftph)  
  real(kind=DP),intent(inout) :: exc
  real(kind=DP),intent(inout) :: dF_drho(ista_fftph:iend_fftph,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout) :: exc_on_a_grid_wk(ista_fftph:iend_fftph)
! -----  ascat ceases modifying  -----
#endif

  real(kind=DP), parameter :: zeta_minimum    = 1.d-4
  real(kind=DP), parameter :: zeta_minimum2   = 1.d-6
  real(kind=DP), parameter :: density_minimum = 1.d-20
  real(kind=DP), parameter :: density_minimum2= 1.d-5

  real(kind=DP), parameter  :: a  = 0.0310907d0,    a1   = 0.21370d0
  real(kind=DP), parameter  :: b1 = 7.5957d0,       b2n  = 3.5876d0
  real(kind=DP), parameter  :: b3 = 1.6382d0,       b4   = 0.49294d0
  real(kind=DP), parameter  :: p  = 1.00d0,         p1   = p+1.d0
  real(kind=DP), parameter  :: ap  = 0.015545d0,    a1p  = 0.20548d0
  real(kind=DP), parameter  :: b1p = 14.1189d0,     b2np = 6.1977d0
  real(kind=DP), parameter  :: b3p = 3.3662d0,      b4p  = 0.62517d0
! vwn 95/12/2 Y.M
  real(kind=DP), parameter  :: aq  = 0.016887d0,    a1q  = 0.11125d0
  real(kind=DP), parameter  :: b1q = 10.3570d0,     b2nq = 3.6231d0
  real(kind=DP), parameter  :: b3q = 0.88026d0,     b4q  = 0.49671d0
! vwn 95/12/2 Y.M
  real(kind=DP), parameter  :: gamma = 0.031091
!--------------------------
  real(kind=DP), parameter  :: xnu = 15.75592d0, cc0 = 0.004235d0
  real(kind=DP), parameter  :: thrd = 0.333333333333d0, sixth7 = 1.1666666666666666d0

  integer       :: is,i
  real(kind=DP) :: facw,bet,delt,exc1,g,g3,facpon,d,dd,rs,fk,sk,t,q0 &
       &         , rs12,rs32,rsp,q1,q2,eu,q3,eurs,ecd,pon &
       &         , b,b2,t2,t4,t6,q4,q5,h0,h,q8,h0t,h0b,h0rs,ht,hrs &
       &         , ec1,ec1d,ec1dd,exc0,excd,excdd &
       &         , thrd2,thrd4,fzdd,zeta,onpzeta,onmzeta,fz,zetadxd &
       &         , q0p,q1p,q2p,q3p,q0q,q1q,q2q,q3q,fzd,onzt4,gzt,bzt,h0zt,euzt

  facw = ispin
  bet = xnu*cc0
  delt = bet / gamma
!---- Spin dependency
#ifdef NEC_TUNE_MXCP
!CDIR NOCONCUR
#endif
  do is = 1, ispin
     exc1 = 0.d0
     if ( ispin == 1 ) then
        g   = 1.d0
        g3  = g**3
!!$        g4  = g3*g
        facpon = -delt/(g3*bet)

#ifdef NEC_TUNE_MXCP
!CDIR INNER
#endif
        do i = ista_fftph,iend_fftph ! MPI
           d  = facw*chgrhr_l(i, 1)
           if(d < density_minimum) cycle
! gradient of charge density 95/12/3 H.S.
           dd = facw* grad_trho(i)
           rs = (0.75/(PAI*d))**thrd
           fk = (3*PAI*PAI*d)**thrd
           sk = dsqrt(4*fk/PAI)
           if(d > density_minimum2) then
              t = dd/(d*sk*2)
!!$              s = dd/(d*fk*2)
           else
              t = 0.d0
!!$              s = 0.d0
           endif

           q0 = -2*a*(1 + a1*rs)
           rs12 = dsqrt(rs)
           rs32 = rs12**3
           rsp = rs**p
           q1 = 2*a*(b1*rs12+b2n*rs+b3*rs32+b4*rs*rsp)
           q2 = log(1 + 1.d0/q1)
           eu = q0*q2
           q3 = a*(b1/rs12+2*b2n+3*b3*rs12+2*b4*p1*rsp)
           eurs = -2*a*a1*q2 - q0*q3/(q1**2 + q1)
           ecd = eu-thrd*rs*eurs
           pon = facpon*eu
           b   = delt/(dexp(pon)-1)
           b2  = b*b
           t2  = t*t
!!$           s2  = s*s
           t4  = t2*t2
           t6  = t4*t2
           q4  = 1 + b*t2
           q5  = 1 + b*t2+b2*t4
           h0  = g3*(bet/delt)*dlog(1+delt*q4*t2/q5)
           h   = h0
           q8  = q5*q5+delt*q4*q5*t2
           h0t = 2*bet*t*(1+2*b*t2)/q8 * g3
           h0b = -bet*t6*(2*b+b2*t2)/q8 * g3
           h0rs = h0b*b*eurs*(b+delt)/ bet / g3
           ht   = h0t
           hrs  = h0rs
           ec1  = d*h/facw
           ec1d = h-thrd*rs*hrs-sixth7*t*ht
           ec1dd = 0.5*ht/sk

           exc0 = eu*d/facw + ec1
           excd = ecd + ec1d
           dF_drho(i, is) = dF_drho(i, is) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(dabs(grad_trho(i)) > 1.0d-9) then
              grad_trho(i) = excdd / grad_trho(i)
           else
              grad_trho(i) = 0.d0
           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d/facw
           exc1 = exc1 + exc0*f2or1(i)
#ifdef __EDA__
           if(sw_eda==ON) then
! -----  ascat starts modifying  -----
           exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + exc0!*f2or1(i)
! -----  ascat ceases modifying  -----
           endif
#endif
        end do
     else if ( ispin ==  2 ) then
        thrd2 = thrd * 2
        thrd4 = thrd * 4
        fzdd = 8/(9*(2**thrd4 - 2 ))

        do i = ista_fftph,iend_fftph
           d = chgrhr_l(i, 1) + chgrhr_l(i, ispin)
           if(d < density_minimum ) cycle
           zeta = ( chgrhr_l(i, 1) - chgrhr_l(i, 2) ) / d
           onpzeta = 2*chgrhr_l(i, 1) / d
           onmzeta = 2*chgrhr_l(i, 2) / d
           g = ( onpzeta**thrd2 + onmzeta**thrd2 )/2
           if(dabs(zeta) < zeta_minimum ) then
              fz= 4.d0/9 *zeta*zeta*(1 + 5.d0/54 *zeta*zeta &
                   &               *(1 +44.d0/135*zeta*zeta )) &
                   &       /(2**thrd4 -2)
           else
              fz = (onpzeta**thrd4 + onmzeta**thrd4 -2) / (2**thrd4 - 2)
           end if
           g3  = g**3
!!$           g4  = g3*g
           facpon = -delt/(g3*bet)
           if(d < density_minimum) cycle
           dd = grad_trho(i)
           rs = (0.75d0/(PAI*d))**thrd
           fk = (3*PAI*PAI*d)**thrd
           sk = dsqrt(4*fk/PAI)

           if(d > density_minimum2) then
! --> Check if this is correct or not
              ! 2011/04/26 revised according to a Katsumata-san's comment.
!!$              t = dd/(d*sk*2)
              t = dd/(d*sk*2)
              t = t/g
! <--
!!$              s = dd/(d*fk*2)
           else
              t = 0.d0
!!$              s = 0.d0
           endif

           q0 = -2*a*(1+a1*rs)
           rs12 = dsqrt(rs)
           rs32 = rs12**3
           rsp = rs**p
           q1 = 2*a*(b1*rs12+b2n*rs+b3*rs32+b4*rs*rsp)
           q2 = log(1+1/q1)
           q3 = a*(b1/rs12+2*b2n+3*b3*rs12+2*b4*p1*rsp)

           q0p = -2*ap*(1+a1p*rs)
           q1p = 2*ap*(b1p*rs12+b2np*rs+b3p*rs32+b4p*rs*rsp)
           q2p = log(1+1/q1p)
           q3p = ap*(b1p/rs12+2*b2np+3*b3p*rs12+2*b4p*p1*rsp)
! vwn 95/12/2 Y.M
           q0q = -2*aq*(1+a1q*rs)
           q1q =  2*aq*(b1q*rs12+b2nq*rs+b3q*rs32+b4q*rs*rsp)
           q2q = log(1+1/q1q)
           q3q = aq*(b1q/rs12+2*b2nq+3*b3q*rs12+2*b4q*p1*rsp)
           if(dabs(zeta) < zeta_minimum) then
              fzd = 8.d0/9.d0 /( 2.d0**thrd4 -2 ) &
                   &     *zeta *( 1 + 5.d0/27 *zeta*zeta &
                   &           *( 1 +22.d0/45 *zeta*zeta ))
           else
              fzd = 4.d0/3.d0 /( 2.d0**thrd4 -2) &
                   &     *( onpzeta**thrd -onmzeta**thrd )
           end if
           onzt4= onmzeta*onpzeta *( 1 + zeta*zeta )
           eu = q0*q2 + ( q0p * q2p - q0 * q2 ) * fz * zeta**4 &
                &             -   q0q * q2q * fz / fzdd * onzt4
           eurs = ( -2*a*a1*q2 -q0*q3/(q1**2 +q1 )) * ( 1 -fz*zeta**4 )&
                &+( -2*ap*a1p*q2p -q0p*q3p/( q1p**2 +q1p )) *fz*zeta**4 &
                &-( -2*aq*a1q*q2q -q0q*q3q/( q1q**2 +q1q )) *fz/fzdd*onzt4
           euzt =  -q0q*q2q /fzdd *( fzd*onzt4 -4*fz*zeta**3 )&
                &        +( q0p*q2p -q0*q2 )*(fzd*zeta +4*fz )*zeta**3 
! vwn 95/12/2 Y.M

           zetadxd = - 2 * chgrhr_l(i, 3-is) * (-1.d0)**is / d
           ecd = eu-thrd*rs*eurs + euzt * zetadxd 

           if(onpzeta < zeta_minimum2 .or. onmzeta < zeta_minimum2) then
              ec1  = 0.d0;    ec1d = 0.d0;     ec1dd= 0.d0
              go to 10
           end if
           pon = facpon*eu
           b   = delt/(dexp(pon)-1)
           b2  = b*b
           t2  = t*t
!!$           s2  = s*s
           t4  = t2*t2
           t6  = t4*t2
           q4  = 1+b*t2
           q5  = 1+b*t2+b2*t4
           h0  = g3*(bet/delt)*dlog(1+delt*q4*t2/q5)
           h   = h0

           q8  = q5*q5+delt*q4*q5*t2
           h0t = 2*bet*t*(1+2*b*t2)/q8 * g3
           h0b = -bet*t6*(2*b+b2*t2)/q8 * g3
           h0rs = h0b * b * eurs * (b+delt) / bet / g3
           ht   = h0t
           hrs  = h0rs

           if(dabs(zeta) < zeta_minimum) then
              gzt = -2.d0/9*zeta *( 1 +14.d0/27 *zeta*zeta &
                   &             *( 1 +13.d0/18 *zeta*zeta ))
           else
              gzt = ( 1.d0/onpzeta**thrd -1.d0/onmzeta**thrd )/3.d0
           end if
           bzt = b * ( b+delt) * ( euzt - 3*eu*gzt/g ) / ( bet * g3 ) 
           h0zt = 3 * gzt * h0 / g + h0b * bzt

           ec1  = d*h * 0.5
! --> revised by Kato-san,  1999Jul16
!!$           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + h0zt * zetadxd
           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxd
           ec1dd = 0.5*ht/sk / g

10         continue
           exc0 = eu*d * 0.5 + ec1
           excd = ecd + ec1d
           dF_drho(i, is) = dF_drho(i, is) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(is ==  2) then
              if(dabs(grad_trho(i)) > 1.0d-9) then
                 grad_trho(i) = excdd / grad_trho(i)
              else
                 grad_trho(i) = 0.d0
              endif
           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d * 0.5
           exc1 = exc1 + exc0*f2or1(i)
#ifdef __EDA__
           if(sw_eda==ON) then
! -----  ascat starts modifying  -----
           exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + exc0!*f2or1(i)
! -----  ascat ceases modifying  -----
           endif
#endif
        end do
     end if
     exc = exc + exc1
  end do
end subroutine cr_ggapbe
#endif
! ==== DEBUG! ==================================================================
#endif

!!$c===================================================================
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
!!$     1994/11/24 coded by T.Yamasaki
!!$
#ifdef __EDA__
subroutine xclda(nspin,ispin,ista_r,iend_r,chgrhr_l,wos,exc,dF_drho,exc_on_a_grid_wk,ien)
#else
subroutine xclda(nspin,ispin,ista_r,iend_r,chgrhr_l,wos,exc,dF_drho,ien)
#endif
  use m_Const_Parameters,  only : DP, PAI
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif
  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in),optional :: ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)  
  real(kind=DP),intent(out) :: exc
!!$  real(kind=DP),intent(inout) :: dF_drho(ista_r:iend_r,nspin)
  real(kind=DP),intent(out) :: dF_drho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(out) :: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif

  integer :: is, i, iend_true
  real(kind=DP) :: co2,d,rs,exclds
  real(kind=DP), parameter :: thrd = 0.333333333333333333d0

  co2 = (0.75d0/PAI)**thrd
  exc = 0.d0
#ifdef __EDA__
  if(sw_eda==ON) then
! -----  ascat starts modifying  -----
  exc_on_a_grid_wk = 0.d0
! -----  ascat ceases modifying  -----
  endif
#endif

  iend_true = iend_r
  if(present(ien)) iend_true = min(iend_r,ien)

  is = 1
  if(ispin == 2)  call phase_error_with_msg(6,'!Sorry, ispin == 2 is not supported now',__LINE__,__FILE__)

  do i = ista_r, iend_true
     d = chgrhr_l(i,is)
     if(d < 1.d-40) then
        dF_drho(i,is) = 0.d0
        cycle
     endif
     rs = co2*(1.d0/d)**thrd
     if(rs >= 1.d0) then
        dF_drho(i,is) = &
             &   (-0.61093d0/rs - (0.1423d0+0.1759d0*dsqrt(rs) &
             &    + 0.06366d0*rs)&
             &   /(1 + 1.0529d0*dsqrt(rs) + 0.3334d0*rs)**2 )
        exclds = -0.61093d0*0.75d0/rs &
             &   -0.1423d0/(1 + 1.0529d0*dsqrt(rs) + 0.3334d0*rs)
     else
        dF_drho(i,is) =  &
             &  ( -0.61093d0/rs - 0.05837d0 - 0.0084d0*rs &
             &   + 0.0311d0*log(rs) + 0.00133d0*rs*dlog(rs) )
        exclds = -0.61093d0*0.75d0/rs &
             &  - 0.048d0 -0.0116*rs &
             &  + 0.0311d0*dlog(rs) + 0.0020d0*rs*dlog(rs)
     end if
     exc = exc + exclds*wos(i)*d
#ifdef __EDA__
     if(sw_eda==ON) then
! -----  ascat starts modifying  -----
     exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + exclds*wos(i)*d
! -----  ascat ceases modifying  -----
     endif
#endif
  end do
end subroutine xclda

!!$c ---------------------------------------
!!$c       1994/11/24 coded by T. Yamasaki
!!$c       1994/12/22 modified by T. Yamasaki
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
#ifdef __EDA__
subroutine ggabek(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dF_drho,dF_dgradrho,exc_on_a_grid_wk,ien)
#else
subroutine ggabek(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dF_drho,dF_dgradrho,ien)
#endif
  use m_Const_Parameters,  only : DP
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif
  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in),optional :: ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)  :: grad_rho(ista_r:iend_r,nspin)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)  
  real(kind=DP),intent(out) :: exc
  real(kind=DP),intent(inout) :: dF_drho(ista_r:iend_r,nspin)
  real(kind=DP),intent(out) :: dF_dgradrho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(out) :: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif

  integer        :: i, iend_true
  real(kind=DP)  :: fac,tthrd,fac_t,d,dd,dthrd,dthrd4,x,x2,p0,x2ovp0&
       &, y,f,vxc1,exgga

  real(kind=DP), parameter ::    b = 0.0042d0
  real(kind=DP), parameter :: thrd = 0.333333333333333333d0
        
  iend_true = iend_r
  if(present(ien)) iend_true = min(iend_r,ien)

  fac = 4*thrd * b
  tthrd = 2.d0**thrd
  fac_t  = 1.d0/tthrd
  do i = ista_r, iend_true
     d = chgrhr_l(i,1)
     dd = grad_rho(i,1)
     if(d < 1.d-10 .or. dd > 1.d10) then
        dF_dgradrho(i,1) = 0.d0
        cycle
     endif
     dthrd  = d**thrd
     dthrd4 = dthrd*d
     x = tthrd*dd/dthrd4
     x2 = x*x
     p0 = dsqrt(1+x2)
     x2ovp0 = 6*b*x2/p0
     if(x < 1.d-40) then
        y = 0.d0
     else
        y = dlog(x + p0)
     endif
     f  = 1.d0/(1.d0+6.d0*b*x*y)
     vxc1    = x2*f*f*(1.d0-x2ovp0)
     dF_drho(i,1) = dF_drho(i,1) + fac*dthrd*vxc1*fac_t
     dF_dgradrho(i,1) = fac*dthrd*(vxc1 + x2*f)*fac_t

     exgga  = -b * dthrd4 * x2 * f * fac_t
     exc      = exc      + exgga*wos(i)
#ifdef __EDA__
     if(sw_eda==ON) then
! -----  ascat starts modifying  -----
     exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + exgga*wos(i)
! -----  ascat ceases modifying  -----
     endif
#endif
  end do
end subroutine ggabek

!!$c ---------------------------------------
!!$c     1994/11/24 coded by T.Yamasaki
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
#ifdef __EDA__
subroutine ggaprd(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dF_drho,dF_dgradrho,exc_on_a_grid_wk,ien)
#else
subroutine ggaprd(nspin,ispin,ista_r,iend_r,chgrhr_l,grad_rho,wos,exc,dF_drho,dF_dgradrho,ien)
#endif
  use m_Const_Parameters,  only : PAI,DP
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif
  implicit none

  integer,intent(in)         :: nspin,ispin,ista_r,iend_r
  integer,intent(in),optional :: ien
  real(kind=DP),intent(in)   :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)   :: grad_rho(ista_r:iend_r,nspin)
  real(kind=DP),intent(in)   :: wos(ista_r:iend_r)  
  real(kind=DP),intent(out)  :: exc
  real(kind=DP),intent(inout):: dF_drho(ista_r:iend_r,nspin)
  real(kind=DP),intent(inout)  :: dF_dgradrho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(out)  :: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif

  real(kind=DP) :: cinfin,d,ovd,thrd4,sixth7,co2,dd,rs,dthrd,dthrd4&
       &, x, cdenom,c,phi,expphi,ecgga,dcdrs,dcdrho
  integer       :: i, iend_true

  real(kind=DP), parameter :: alpha = 0.023266d0
  real(kind=DP), parameter ::  beta = 7.389d-6
  real(kind=DP), parameter :: gamma = 8.723d0
  real(kind=DP), parameter :: delta = 0.472d0
  real(kind=DP), parameter ::fchird = 0.11
  real(kind=DP), parameter :: thrd  = 0.333333333333333333d0
  real(kind=DP), parameter :: c0    = 0.001667d0
  real(kind=DP), parameter :: c1    = 0.002568d0

  iend_true = iend_r
  if(present(ien)) iend_true = min(iend_r,ien)

  cinfin = c0 + c1
  d = 1.d0
  ovd = 1.d0/d
  thrd4 = 4*thrd
  sixth7 = 7.d0/6.d0

  co2 = (0.75d0/PAI)**thrd
  do i = ista_r, iend_true
     d = chgrhr_l(i,1)
     dd = grad_rho(i,1)
     if(d < 1.d-40) cycle
     rs = co2*(1.d0/d)**thrd
     dthrd  = d**thrd
     dthrd4 = dthrd*d
     x = dd/dthrd4
     cdenom = 1.d0/(1 + rs*(gamma + rs*(delta + 1.d+4*beta*rs)))
     c = c0 + (c1 + rs*(alpha + beta*rs))*cdenom
     phi = 1.745d0 * fchird * cinfin/c * x *dsqrt(dthrd)
     expphi = dexp(-phi)
     ecgga = ovd*expphi*c*dd*x

     dcdrs = cdenom * ((alpha + 2*beta*rs) &
          &         - cdenom * (c1+rs*(alpha+beta*rs)) &
          &       * (gamma+rs*(2*delta + 3.d+4*beta*rs)))
     dcdrho = -thrd*(rs/d)*dcdrs
     dF_drho(i,1) = dF_drho(i,1) + ovd*(dcdrho*(phi+1) &
          &       +(sixth7*phi - thrd4)*c)*expphi*x*dd
     dF_dgradrho(i,1) = dF_dgradrho(i,1) + ovd*expphi*c*x*(2 - phi)

     exc = exc + ecgga*wos(i)
#ifdef __EDA__
     if(sw_eda==ON) then
! -----  ascat starts modifying  -----
     exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + ecgga*wos(i)
! -----  ascat ceases modifying  -----
     endif
#endif
  end do
end subroutine ggaprd

subroutine cr_lda(nspin,ispin,ista_r,iend_r,chgrhr_l,exc,dF_drho,ien)
  use m_Const_Parameters,  only : PAI,DP
  integer, intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in),optional :: ien
  real(kind=DP), intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP), intent(inout) :: exc
  real(kind=DP), intent(inout) :: dF_drho(ista_r:iend_r,ispin)
  real(kind=DP) :: rhomin = 1.d-9
  integer :: i, iend_true
  real(kind=DP) :: n,rs,x,ec,dedx,dxdr
  iend_true = iend_r
  if(present(ien)) iend_true = min(iend_r,ien)
  ec = 0.d0
  do i=ista_r,iend_true
     n = chgrhr_l(i,ispin)
!     if(n<rhomin) n = rhomin
     if(n<rhomin) cycle
     
     rs = ((3.d0/(4.d0*PAI*n))**(1.d0/3.d0))
     x = rs/11.4d0
     ec = -0.0666d0*0.5d0*((1.d0+x**3)*DLOG(1.d0+1.d0/x)-x**2+x/2.d0-1.d0/3.d0)
     exc = exc+ec*n

     dxdr = -(3.d0/(4.d0*PAI))**(1.d0/3.d0)/11.4d0/3.d0*n**(-4.d0/3.d0)
     dedx = -0.0333d0*(3*x*x*dlog(1.d0+1.d0/x)-(1+x**3)/(x*(x+1))-2*x+0.5d0)
     dF_drho(i,ispin) = dF_drho(i,ispin) + dedx*dxdr*n+ec
  enddo
end subroutine cr_lda

#ifdef __EDA__
subroutine xcpotf_wigner(ispin,ista_r,iend_r,input_charge,DELTA,chgrhr_l,wos,exc,exc_on_a_grid_wk)
#else
subroutine xcpotf_wigner(ispin,ista_r,iend_r,input_charge,DELTA,chgrhr_l,wos,exc)
#endif
!*****<<< WIGNER INTERPOLATION >>>
  use m_Const_Parameters,  only : DP,UP,DOWN,PAI,Valence_plus_PC_Charge
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif
  implicit none

  integer,intent(in)         :: ispin,ista_r,iend_r,input_charge
  real(kind=DP),intent(in)   :: DELTA
  real(kind=DP),intent(inout):: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)   :: wos(ista_r:iend_r)  
  real(kind=DP),intent(inout):: exc
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif

  real(kind=DP) :: rho,rs,zeta,vxc,tmp,WX0,CO2
  integer       :: i, icount
  real(kind=DP), parameter :: C13 = 1.d0/3.d0
  real(kind=DP), parameter :: CO3 = 4.d0/3.d0

  CO2 = (0.75d0/PAI)**C13
  WX0 =  3*(((9*PAI)/4)**C13)/(2*PAI)

  if(ispin == 2 .and. input_charge == Valence_plus_PC_Charge) then
     icount = 0
     do i = ista_r, iend_r
        rho    = chgrhr_l(i,UP) + chgrhr_l(i,DOWN)
        rs     = CO2*( 1.0d0/rho  )**C13
        if(rho < DELTA) cycle
        zeta = (chgrhr_l(i,UP) - chgrhr_l(i,DOWN))/rho
        if(dabs(zeta) > DELTA) icount = icount + 1
        tmp    = - WX0/2.d0*CO3/RS -(0.440d0*CO3*RS + 3.432d0)/(RS+7.8D0)**2
        chgrhr_l(i,UP  ) = tmp
        chgrhr_l(i,DOWN) = tmp
        vxc    = -(WX0/2.d0/RS+0.44d0/(RS+7.8d0))
        exc    = exc + vxc*rho*wos(i)
#ifdef __EDA__
        if(sw_eda==ON) then
! -----  ascat starts modifying  -----
        exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + vxc*rho*wos(i)
! -----  ascat ceases modifying  -----
        endif
#endif
     end do
     if(icount >= 1) then
        write(6,*)' WIGNER XC-ENERGY and POTENTIAL are not' &
             &         , ' appropriate for Spin-Polarized calculation'
        call phase_error_with_msg(6,'WIGNER',__LINE__,__FILE__)
     endif
  else
     do i = ista_r, iend_r
        rho    = chgrhr_l(i,1)
        rs     = CO2*( 1.0d0/rho  )**C13
        chgrhr_l(i,1)   = - WX0/2.d0*CO3/rs &
             &           -(0.440d0*CO3*rs + 3.432d0)/(rs+7.8d0)**2
        vxc    = -(WX0/2.d0/rs+0.44d0/(rs+7.8d0))
        exc    = exc + vxc*rho*wos(i)
#ifdef __EDA__
        if(sw_eda==ON) then
! -----  ascat starts modifying  -----
        exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + vxc*rho*wos(i)
! -----  ascat ceases modifying  -----
        endif
#endif
     end do
  endif
end subroutine xcpotf_wigner

#ifdef __EDA__
subroutine xcpotf_pzold(ispin,ista_r,iend_r,input_charge,DELTA,chgrhr_l,wos,exc,exc_on_a_grid_wk)
#else
subroutine xcpotf_pzold(ispin,ista_r,iend_r,input_charge,DELTA,chgrhr_l,wos,exc)
#endif
  use m_Const_Parameters,  only : DP,UP,DOWN,PAI,Valence_plus_PC_Charge
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif
  implicit none

  integer,intent(in)         :: ispin,ista_r,iend_r,input_charge
  real(kind=DP),intent(in)   :: DELTA
  real(kind=DP),intent(inout):: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)   :: wos(ista_r:iend_r)  
  real(kind=DP),intent(inout):: exc
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif

  real(kind=DP) :: B1P01,B1P02,CP00,CP01,CP02,OVC213,OVC243 &
       &         , B1F01,B1F02,CF00,CF01,CF02
  real(kind=DP) :: rho,ovn,rs,denmp,denmf,wcp,wcf,vcp,vcf,zeta,f00,f2up,f2down &
       &         , t1,wxp,wxf,vxp,x,t2,vxc
  integer       :: i

  real(kind=DP), parameter   :: C13   =  1.d0/3.d0
  real(kind=DP), parameter   :: C43   =  4.d0/3.d0
  real(kind=DP) :: CO2, WX0, C213

  real(kind=DP), parameter :: GAMP   = -0.1423d0, GAMF   = -0.0843d0 &
       &                    , BETA1P =  1.0529d0, BETA2P =  0.3334d0 &
       &                    , BETA1F =  1.3981d0, BETA2F =  0.2611d0 &
       &                    , AP     =  0.0311d0, BP     = -0.0480d0 &
       &                    , AF     =  0.01555d0,BF     = -0.0269d0 &
       &                    , CP     =  0.0020d0, DPx    = -0.0116d0 &
       &                    , CF     =  0.00070d0,DF     = -0.0048d0

  CO2   = (0.75d0/PAI)**C13
  WX0   =  3*(((9*PAI)/4)**C13)/(2*PAI)
  C213  = 2.d0**C13

  B1P01 = 7.d0*BETA1P/6.d0
  B1P02 = 4.d0*BETA2P/3.d0

  CP00 = BP - AP/3.d0
  CP01 = (2*DPx - CP)/3.d0
  CP02 = 2*CP/3.d0

  if(ispin == 2 .and. input_charge == Valence_plus_PC_Charge) then
     OVC213 = 1.D0/(C213 - 1)
     OVC243 = 1.D0/(2.D0**C43 - 2)

     B1F01 = 7.D0*BETA1F/6.D0
     B1F02 = 4.D0*BETA2F/3.D0
     CF00 = BF - AF/3.D0
     CF01 = (2*DF - CF)/3.D0
     CF02 = 2*CF/3.D0

     do i = ista_r, iend_r   ! MPI
        rho    = chgrhr_l(i, UP) + chgrhr_l(i,DOWN)
        if(rho <= DELTA) then
           chgrhr_l(i,UP  ) = 0.D0
           chgrhr_l(i,DOWN) = 0.D0
           cycle
        endif
        ovn = 1.d0/rho
        rs     = CO2*ovn**C13
        if(rs >= 1.d0) then
           denmp = 1 + BETA1P*sqrt(rs) + BETA2P*rs
           denmf = 1 + BETA1F*sqrt(rs) + BETA2F*rs
           wcp = GAMP/denmp
           wcf = GAMF/denmf
           vcp = wcp*(1+B1P01*sqrt(rs)+B1P02*rs)/denmp
           vcf = wcf*(1+B1F01*sqrt(rs)+B1F02*RS)/denmf
        else
           wcp = AP*dlog(rs) + BP + CP*rs*dlog(rs) + DPx*rs
           wcf = AF*dlog(rs) + BF + CF*rs*dlog(rs) + DF*rs
           vcp = (AP+CP02*RS)*dLOG(RS) + CP00 + CP01*RS
           vcf = (AF+CF02*RS)*dLOG(RS) + CF00 + CF01*RS
        end if
        zeta   = (chgrhr_l(i,UP) - chgrhr_l(i,DOWN))*ovn
        f00    = ((1+zeta)**C43 + (1-zeta)**C43 - 2)*OVC243
        f2up   = ((1+zeta)**C13 - 1)*OVC213
        f2down = ((1-zeta)**C13 - 1)*OVC213
        t1 = vcf - vcp - C43*(wcf - wcp)
        wxp = -WX0*0.5d0/rs
        wxf = wxp * C213
        vxp = wxp * C43
        x = vxp + vcp + t1*f00
        t2 = C43 * ((WXF - WXP) + ( WCF - WCP ))
        chgrhr_l(i,UP  ) = x + t2*f2up
        chgrhr_l(i,DOWN) = x + t2*f2down

        vxc = wxp + wcp + ((wxf - wxp) + (wcf - wcp))*f00
        exc = exc + vxc*rho*wos(i)
#ifdef __EDA__
        if(sw_eda==ON) then
! -----  ascat starts modifying  -----
        exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + vxc*rho*wos(i)
! -----  ascat ceases modifying  -----
        endif
#endif
     end do
  else
     do i = ista_r, iend_r
        rho    = chgrhr_l(i,1)
        if(rho < DELTA) then
           chgrhr_l(i,1) = 0.d0
           cycle
        endif
        ovn    = 1.d0/rho
        rs     = CO2*ovn**C13
        wxp    = -WX0*0.5d0/rs
        if(rs >= 1.d0) then
           denmp = 1 + BETA1P*sqrt(rs) + BETA2P*rs
           wcp = GAMP/denmp
           vcp = wcp*(1.D0+B1P01*SQRT(rs)+B1P02*rs)/denmp
        else
           wcp = AP*dlog(rs) + BP + CP*rs*log(rs) + DPx*rs
           vcp = (AP+CP02*rs)*log(rs) + CP00 + CP01*rs
        end if

        wxp = -WX0*0.5D0/rs
        wxf = WXP * C213
        vxp = WXP * C43
        x = vxp + vcp
        chgrhr_l(i,1  ) = x
        vxc = wxp + wcp
        exc = exc + vxc*rho*wos(i)
#ifdef __EDA__
        if(sw_eda==ON) then
! -----  ascat starts modifying  -----
        exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + vxc*rho*wos(i)
! -----  ascat ceases modifying  -----
        endif
#endif
     end do
  endif
end subroutine xcpotf_pzold

#ifdef __EDA__
subroutine xcpotf_xalfa(ispin,ista_r,iend_r,input_charge,DELTA,chgrhr_l,wos,exc,exc_on_a_grid_wk)
#else
subroutine xcpotf_xalfa(ispin,ista_r,iend_r,input_charge,DELTA,chgrhr_l,wos,exc)
#endif
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,UP,DOWN,PAI,Valence_plus_PC_Charge
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif
  implicit none

  integer,intent(in)         :: ispin,ista_r,iend_r,input_charge
  real(kind=DP),intent(in)   :: DELTA
  real(kind=DP),intent(inout):: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)   :: wos(ista_r:iend_r)  
  real(kind=DP),intent(inout):: exc
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif

  real(kind=DP), parameter   :: C13   =  1.d0/3.d0
  real(kind=DP), parameter   :: ALPHA = 0.7d0
  real(kind=DP)              :: CO1, CO2
  real(kind=DP)              :: rho,tmp,vxc
  integer                    :: i

  CO1   = 1.5d0*ALPHA*(3.D0/PAI)**C13
  CO2   = (0.75d0/PAI)**C13

  if(ispin == 2 .and. input_charge == Valence_plus_PC_Charge) then
     do i = ista_r, iend_r
        rho    = chgrhr_l(i,UP) + chgrhr_l(i,DOWN)
!!$        rs     = CO2*( 1.0D0/rho  )**C13
        tmp    = - CO1*(rho**C13)
        chgrhr_l(i,UP)   = tmp
        chgrhr_l(i,DOWN) = tmp
        vxc = 0.25d0*CO1*(rho**C13)
        exc = exc+vxc*rho*wos(i)
#ifdef __EDA__
        if(sw_eda==ON) then
! -----  ascat starts modifying  -----
        exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + vxc*rho*wos(i)
! -----  ascat ceases modifying  -----
        endif
#endif
     end do
  else
     do i = ista_r, iend_r
        rho    = chgrhr_l(i,1)
!!$        rs     = CO2*( 1.0d0/rho  )**C13
        chgrhr_l(i,1)   = - CO1*(rho**C13)
        vxc = 0.25d0*CO1*(rho**C13)
        exc = exc+vxc*rho*wos(i)
#ifdef __EDA__
        if(sw_eda==ON) then
! -----  ascat starts modifying  -----
        exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + vxc*rho*wos(i)
! -----  ascat ceases modifying  -----
        endif
#endif
     end do
  endif
end subroutine xcpotf_xalfa

#ifdef __EDA__
subroutine xcpotf_pz(ispin,ista_r,iend_r,input_charge,DELTA,chgrhr_l,wos,exc,exc_on_a_grid_wk)
#else
subroutine xcpotf_pz(ispin,ista_r,iend_r,input_charge,DELTA,chgrhr_l,wos,exc)
#endif
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,UP,DOWN,PAI,Valence_plus_PC_Charge
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif
  implicit none

  integer,intent(in)         :: ispin,ista_r,iend_r,input_charge
  real(kind=DP),intent(in)   :: DELTA
  real(kind=DP),intent(inout):: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)   :: wos(ista_r:iend_r)  
  real(kind=DP),intent(inout):: exc
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif

  real(kind=DP), parameter :: C13   =  1.d0/3.d0,  C23   =  2.d0/3.d0,  C43   =  4.d0/3.d0
  real(kind=DP), parameter :: GAMP   = -0.1423d0, GAMF   = -0.0843d0 &
       &                    , BETA1P =  1.0529d0, BETA2P =  0.3334d0 &
       &                    , BETA1F =  1.3981d0, BETA2F =  0.2611d0 &
       &                    , AP     =  0.0311d0, BP     = -0.0480d0 &
       &                    , CP     =  0.0020d0, DPx    = -0.0116d0 &
       &                    , AF     =  0.01555d0,BF     = -0.0269d0 &
       &                    , CF     =  0.00070d0,DF     = -0.0048d0

  real(kind=DP)  :: B1P01,B1P02,CP00,CP01,CP02,C243 &
       &          , B1F01,B1F02,CF00,CF01,CF02
  real(kind=DP)  :: CO2, WX0, C213
  real(kind=DP)  :: rho,ovn,rs,zeta,wxp,f00,wx,vx1,vx2,denmp,denmf,wcp,wcf &
       &          , vcp,vcf,wc,t1,f11,t2,vc1,vc2,vxc,vxp
  integer        :: i

  CO2   = (0.75d0/PAI)**C13
  WX0   =  3*(((9*PAI)/4)**C13)/(2*PAI)
  B1P01 = 7*BETA1P/6.d0
  B1P02 = 4*BETA2P/3.d0
  CP00 = BP - AP/3.D0
  CP01 = (2*DPx - CP)/3.D0
  CP02 = 2*CP/3.D0

  if(ispin == 2 .and. input_charge == Valence_plus_PC_Charge) then
     C213 = 2.d0**C13 - 1
!!$     OVC213 = 1.d0/C213
     C243   = 2*C213
!           ovc243 = 1.d0/(2.d0**C43 - 2)

     B1F01 = 7*BETA1F/6.d0
     B1F02 = 4*BETA2F/3.d0
     CF00 = BF - AF/3.d0
     CF01 = (2*DF - CF)/3.d0
     CF02 = 2*CF/3.d0

     do i = ista_r, iend_r
        rho    = chgrhr_l(i, UP) + chgrhr_l(i,DOWN)
        if(rho <= DELTA) then
           chgrhr_l(i,UP  ) = 0.D0
           chgrhr_l(i,DOWN) = 0.D0
           cycle
        end if
        ovn = 1.d0/rho
        rs     = CO2*ovn**C13
        zeta   = (chgrhr_l(i,UP) - chgrhr_l(i,DOWN))*ovn
        wxp    = -WX0*0.5d0/rs
        f00 = ( (1+zeta)**C43 + (1-zeta)**C43 - 2 ) / c243
        wx  = wxp*(1+C213*f00)
        vx1 = c43*wxp*(1+zeta)**c13
        vx2 = c43*wxp*(1-zeta)**c13
        if(rs >= 1.d0) then
           denmp = 1.d0 + BETA1P*sqrt(rs) + BETA2P*rs
           denmf = 1.D0 + BETA1F*sqrt(rs) + BETA2F*rs
           wcp = GAMP/denmp
           wcf = GAMF/denmf
           vcp = wcp*(1.d0+B1P01*sqrt(rs)+B1P02*rs)/denmp
           vcf = wcf*(1.D0+B1F01*sqrt(rs)+B1F02*rs)/denmf
        else
           wcp = AP*dlog(rs) + BP + CP*rs*dlog(rs) + DPx*rs
           wcf = AF*dlog(rs) + BF + CF*rs*dlog(rs) + DF*rs
           vcp = (AP+CP02*rs)*dLOG(rs) + CP00 + CP01*rs
           vcf = (AF+CF02*rs)*dLOG(rs) + CF00 + CF01*rs
        end if
        wc     = wcp + f00*(wcf - wcp)
        t1     = vcp + f00*(vcf - vcp)
        f11    = C23/C213*((1+zeta)**c13-(1-zeta)**c13)
        t2     = f11*(wcf-wcp)
        vc1    = t1+t2*(1-zeta)
        vc2    = t1-t2*(1+zeta)

        chgrhr_l(i,UP  ) = vx1 + vc1
        chgrhr_l(i,DOWN) = vx2 + vc2

        vxc    = wx + wc
        exc = exc + vxc*rho*wos(i)
#ifdef __EDA__
        if(sw_eda==ON) then
! -----  ascat starts modifying  -----
        exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + vxc*rho*wos(i)
! -----  ascat ceases modifying  -----
        endif
#endif
     end do
  else
     do i = ista_r, iend_r
        rho    = chgrhr_l(i, 1)
        if(rho <= DELTA) then
           chgrhr_l(i,UP  ) = 0.D0
           cycle
        end if
        ovn = 1.d0/rho
        rs     = CO2*ovn**C13
        wxp    = -WX0*0.5d0/rs
        vxp = C43*wxp
        if(rs >= 1.d0) then
           denmp = 1 + BETA1P*sqrt(rs) + BETA2P*rs
           wcp = GAMP/denmp
           vcp = wcp*(1+B1P01*sqrt(rs)+B1P02*rs)/denmp
        else
           wcp = AP*dlog(rs) + BP + CP*rs*dlog(rs) + DPx*rs
           vcp = (AP+CP02*rs)*dLOG(rs) + CP00 + CP01*rs
        end if
        wc     = wcp
        chgrhr_l(i,1   ) = vxp + vcp

        vxc    = wxp + wcp
        exc = exc + vxc*rho*wos(i)
#ifdef __EDA__
        if(sw_eda==ON) then
! -----  ascat starts modifying  -----
        exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + vxc*rho*wos(i)
! -----  ascat ceases modifying  -----
        endif
#endif
     end do
  endif
end subroutine xcpotf_pz

#ifdef __EDA__
subroutine xcpotf_vwn(ispin,ista_r,iend_r,input_charge,DELTA,chgrhr_l,wos,exc,exc_on_a_grid_wk)
#else
subroutine xcpotf_vwn(ispin,ista_r,iend_r,input_charge,DELTA,chgrhr_l,wos,exc)
#endif
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,UP,DOWN,PAI,Valence_plus_PC_Charge
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif
  implicit none

  integer,intent(in)         :: ispin,ista_r,iend_r,input_charge
  real(kind=DP),intent(in)   :: DELTA
  real(kind=DP),intent(inout):: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)   :: wos(ista_r:iend_r)  
  real(kind=DP),intent(inout):: exc
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif

  real(kind=DP), parameter :: B1 = 1.13107d0,    C1 = 13.0045d0 &
       &                    , D1 = -0.0047584d0 &
       &                    , A2 =  0.0621814d0, B2 =  3.72744d0 &
       &                    , C2 = 12.9352d0   , D2 = -0.10498d0 &
       &                    , A3 =  0.0310907d0, B3 =  7.06042d0 &
       &                    , C3 = 18.0578d0,    D3 = -0.32500d0
  real(kind=DP), parameter :: C13   =  1.d0/3.d0,  C43   =  4.d0/3.d0
  real(kind=DP) :: CO2, WX0, C213

  real(kind=DP) :: F20,A1,Q1,P01,Q2,P02,Q3,P03,C243
  real(kind=DP) :: rho,ovn,zeta,f00,wx,vx1,vx2,x,p,r,h00,h10,h0p,h1p,h0f,h1f &
       &         , z4,h000,wc,vc,g00,vc0,vc1,vc2,rs,wxp,vxc
  integer       :: i

  CO2   = (0.75d0/PAI)**C13
  WX0   =  3*(((9*PAI)/4)**C13)/(2*PAI)
  C213 = 2**C13 - 1

  F20   = 4/(9*C213)
  A1   = -1/(3*PAI*PAI)
  Q1   = dsqrt(4*C1-B1**2)
  P01  = D1**2 + B1*D1 + C1

  Q2   = dsqrt(4*C2-B2**2)
  P02  = D2**2 + B2*D2 + C2

  Q3   = dsqrt(4*C3-B3**2)
  P03  = D3**2 + B3*D3 + C3

  if(ispin ==  2 .and. input_charge == Valence_plus_PC_Charge) then
!!$     OVC213 = 1.d0/C213
     C243   = 2*C213
!           ovc243 = 1.d0/(2.d0**c43 - 2)

     do i = ista_r, iend_r
        rho    = chgrhr_l(i, UP) + chgrhr_l(i,DOWN)
        if(rho <= DELTA) then
           chgrhr_l(i,UP  ) = 0.D0
           chgrhr_l(i,DOWN) = 0.D0
           cycle
        end if
        ovn = 1.d0/rho
        rs     = CO2*ovn**C13
        if(rs <= DELTA) then
           chgrhr_l(i, UP)   = 0.d0  
           chgrhr_l(i, DOWN) = 0.d0
           cycle
        end if
        zeta   = (chgrhr_l(i,UP) - chgrhr_l(i,DOWN))*ovn
        wxp    = -WX0*0.5d0/rs
        f00 = ( (1+zeta)**C43 + (1-zeta)**C43 - 2 ) / C243
        wx  = wxp*(1+C213*f00)
        vx1 = C43*wxp*(1+zeta)**C13
        vx2 = C43*wxp*(1-zeta)**C13
        x   = dsqrt(rs)
        p   = x**2 + B1*x + C1
        r   = datan(Q1/(2*x+B1))
        h00= (A1/F20) * ( dlog((x**2)/p)+2*B1*r/Q1 &
             &       -B1*D1*(dlog(((x-D1)**2)/p)+2*(B1+2*D1)*r/Q1)/P01 )
        h10= (A1/F20) * ( 2/x - (1-B1*D1/P01)*(2*x+B1)/p - 2*B1*D1/((x-D1)*P01) &
             &                 - 4*B1*(1-(B1+2*D1)*D1/P01)/(Q1**2+(2*x+B1)**2) )

        p   = x**2 + B2*x + C2
        r   = datan(Q2/(2*x+B2))
        h0p= A2 * ( dlog((x**2)/p)+2*B2*r/Q2 &
             &        -B2*D2*(dlog(((x-D2)**2)/p)+2*(B2+2*D2)*r/Q2)/P02 )
        h1p= A2 * ( 2/x - (1-B2*D2/P02)*(2*x+B2)/p - 2*B2*D2/((x-D2)*P02) &
             &        -4*B2*(1-(B2+2*D2)*D2/P02)/(Q2**2+(2*x+B2)**2) )

        p   = x**2 + B3*x + C3
        r   = datan(Q3/(2*x+B3))
        h0f = A3 * ( dlog((x**2)/p)+2*B3*r/Q3 &
             &       -B3*D3*(dlog(((x-D3)**2)/p)+2*(B3+2*D3)*r/Q3)/P03 )
        h1f = A3 * ( 2/x - (1-B3*D3/P03)*(2*x+B3)/p - 2*B3*D3/((x-D3)*P03) &
             &         - 4*B3*(1-(B3+2*D3)*D3/P03)/(Q3**2+(2*x+B3)**2) )

        z4   = zeta**4
        h000 = h0f - h0p - h00
        wc=h0p+h00*f00+h000*z4*f00
        vc=wc-x*((1-z4*f00)*h1p+z4*f00*h1f+(1-z4)*f00*h10)/6
        g00 = 4*((1+zeta)**c13-(1-zeta)**c13) / (6*c213)
        vc0 = 4* h000*(zeta**3)*f00 + (h000*z4+h00)*g00
        vc1=vc+(1-zeta)*vc0
        vc2=vc-(1+zeta)*vc0

        wc  = wc*0.5d0
        vc1 = vc1*0.5d0
        vc2 = vc2*0.5d0

        chgrhr_l(i,UP  ) = vx1 + vc1
        chgrhr_l(i,DOWN) = vx2 + vc2

        vxc    = wx + wc
        exc = exc + vxc*rho*wos(i)
#ifdef __EDA__
        if(sw_eda==ON) then
! -----  ascat starts modifying  -----
        exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + vxc*rho*wos(i)
! -----  ascat ceases modifying  -----
        endif
#endif
     end do
  else
!!$     OVC213 = 1.d0/C213
     C243   = 2*C213

     do i = ista_r, iend_r
        rho    = chgrhr_l(i, 1)
        if(rho <= DELTA) then
           chgrhr_l(i,1  ) = 0.D0
           cycle
        end if
        ovn = 1.d0/rho
        rs     = CO2*ovn**C13
        wxp    = -WX0*0.5d0/rs
        f00 = 0.d0
        wx  = wxp
        vx1 = c43*wxp

        x   = dsqrt(rs)
        p   = x**2 + B2*x + C2
        r   = datan(Q2/(2*x+B2))
        h0p= a2 * ( dlog((x**2)/p)+2*b2*r/q2 &
             &       -b2*d2*(dlog(((x-d2)**2)/p)+2*(b2+2*d2)*r/q2)/p02 )
        h1p= a2 * ( 2/x - (1-b2*d2/p02)*(2*x+b2)/p - 2*b2*d2/((x-d2)*p02) &
             &           - 4*b2*(1-(b2+2*d2)*d2/p02)/(q2**2+(2*x+b2)**2) )

        wc  = h0p
        vc1 = wc-x*h1p/6

        wc  = wc*0.5d0
        vc1 = vc1*0.5d0

        chgrhr_l(i,1  ) = vx1 + vc1
        vxc    = wx + wc
        exc = exc + vxc*rho*wos(i)
#ifdef __EDA__
! -----  ascat starts modifying  -----
        if(sw_eda==ON) then
        exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + vxc*rho*wos(i)
        endif
! -----  ascat ceases modifying  -----
#endif
     end do
  endif
end subroutine xcpotf_vwn

#ifdef __EDA__
subroutine xcpotf_mjw_bh_gl(len_xctype,xctype,ispin,ista_r,iend_r,input_charge,DELTA,chgrhr_l,wos,exc,exc_on_a_grid_wk)
#else
subroutine xcpotf_mjw_bh_gl(len_xctype,xctype,ispin,ista_r,iend_r,input_charge,DELTA,chgrhr_l,wos,exc)
#endif
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,UP,DOWN,PAI,Valence_plus_PC_Charge
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif
  implicit none

  integer,intent(in)         :: len_xctype,ispin,ista_r,iend_r,input_charge
  character(len=len_xctype)  :: xctype
  real(kind=DP),intent(in)   :: DELTA
  real(kind=DP),intent(inout):: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)   :: wos(ista_r:iend_r)  
  real(kind=DP),intent(inout):: exc
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif

  real(kind=DP), parameter ::  CP    =  0.0020d0,  CF     =  0.00070d0
  real(kind=DP), parameter ::  C13   =  1.d0/3.d0, C43   =  4.d0/3.d0

  real(kind=DP) :: CO2, WX0, C213
  real(kind=DP) :: C243, CPX, RP, CFX, RF
  real(kind=DP) :: rho,ovn,rs,zeta,wxp,f00,wx,vx1,vx2,rsp,rsf,aap,aaf &
       &         , wcp,wcf,vcp,vcf,t1,t2,f11,f12,wc,vc1,vc2,vxc
  integer       :: i

  CO2   = (0.75d0/PAI)**C13
  WX0   =  3*(((9*PAI)/4)**C13)/(2*PAI)
  C213 = 2**C13 - 1
  if(xctype ==  'mjw    ') then
! --*   moruzzi-janak-williams; Phys.Rev.12(1975)1257.  
     CPX=0.0450d0
     RP=21.0d0
  else if(xctype == 'bh     ') then                                    
! --*   von barth-hedin; J.Phys.C:Solid State Phys.5(1972)1629.  
     CPX=0.0504d0
     RP=30.0d0
  else if(xctype == 'gl     ') then                                    
! --*   gunnarsson-lundqvist; Phys.Rev.B13(1976)4274. 
     CPX=0.0666d0
     RP=11.4d0
  end if
  if(ispin ==  2 .and. input_charge == Valence_plus_PC_Charge) then
     if(xctype ==  'mjw    ') then
        CFX=cp/2
        RF=RP*(2.d0**C43)
     else if(xctype ==  'bh     ') then                                    
        CFX=0.0254d0
        RF=75.0d0
     else if(xctype ==  'gl     ') then                                    
        CFX=0.0406d0
        RF=15.9d0
     end if

     C243   = 2*C213
     do i = ista_r, iend_r
        rho    = chgrhr_l(i, UP) + chgrhr_l(i,DOWN)
        if(rho <= DELTA) then
           chgrhr_l(i,UP  ) = 0.D0
           chgrhr_l(i,DOWN) = 0.D0
           cycle
        end if
        ovn = 1.d0/rho
        rs     = CO2*ovn**C13
        if(rs <= DELTA) then
           chgrhr_l(i, UP  ) = 0.d0
           chgrhr_l(i, DOWN) = 0.d0
           cycle
        end if
        zeta   = (chgrhr_l(i,UP) - chgrhr_l(i,DOWN))*ovn
        wxp    = -WX0*0.5d0/rs
        f00 = ( (1+zeta)**C43 + (1-zeta)**C43 - 2 ) / C243
        wx  = wxp*(1+C213*f00)
        vx1 = C43*wxp*(1+zeta)**C13
        vx2 = C43*wxp*(1-zeta)**C13
        rsp = rs/RP
        rsf = rs/RF
        aap = dlog(1+1/rsp)                                      
        aaf = dlog(1+1/rsf)                                      
        wcp = - CPX * ((1+rsp**3)*aap - rsp**2 + rsp/2 - C13)
        wcf = - CFX * ((1+rsf**3)*aaf - rsf**2 + rsf/2 - C13)
        vcp = - CP * aap
        vcf = - CF * aaf

        t1  = vcf - vcp - C43*(wcf-wcp) 
        t2  = C43 * (wcf-wcp)

        f11= ( (1+zeta)**C13 -1 ) / C213
        f12= ( (1-zeta)**C13 -1 ) / C213
        wc=wcp + (wcf - wcp)*f00
        vc1=vcp + t1*f00 + t2*f11
        vc2=vcp + t1*f00 + t2*f12

        wc = wc*0.5d0
        vc1 = vc1*0.5d0
        vc2 = vc2*0.5d0

        chgrhr_l(i,UP  ) = vx1 + vc1
        chgrhr_l(i,DOWN) = vx2 + vc2

        vxc    = wx + wc
        exc = exc + vxc*rho*wos(i)
#ifdef __EDA__
        if(sw_eda==ON) then
! -----  ascat starts modifying  -----
        exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + vxc*rho*wos(i)
! -----  ascat ceases modifying  -----
        endif
#endif
     end do
  else
     c243   = 2*C213
     do i = ista_r, iend_r
        rho    = chgrhr_l(i, 1)
        if(rho <= DELTA) then
           chgrhr_l(i,1   ) = 0.D0
           cycle
        end if
        ovn = 1.d0/rho
        rs     = CO2*ovn**C13
        if(rs <= DELTA) then
           chgrhr_l(i, 1   ) = 0.d0
           cycle
        end if
        zeta   = 0.d0
        wxp    = -WX0*0.5d0/rs
        wx  = wxp
        vx1 = c43*wxp
        vx2 = c43*wxp

        rsp = rs/RP
        aap = dlog(1+1/rsp)                                      
        wcp = - CPX * ((1+rsp**3)*aap - rsp**2 + rsp/2 - C13)
        vcp = - CP * aap

        wc=wcp * 0.5d0
        vc1=vcp * 0.5d0

        chgrhr_l(i,1   ) = vx1 + vc1

        vxc    = wx + wc
        exc = exc + vxc*rho*wos(i)
#ifdef __EDA__
        if(sw_eda==ON) then
! -----  ascat starts modifying  -----
        exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + vxc*rho*wos(i)
! -----  ascat ceases modifying  -----
        endif
#endif
     end do
  endif
end subroutine xcpotf_mjw_bh_gl
!$$#endif
! === xclda etc. are necessary to make 3D_Parallel, too!!! by tkato ============
!!BRANCH_P_END ORG_Parallel
! ==============================================================================

!! --- ggasubs ---
!! -------------------------------------------
subroutine cpden(rh,kmesh,mesh,lfp,rad,rho,fdiff)
!      30th Oct. 1994
!             by T.Yamasaki
!
  use m_Const_Parameters, only : DP,PAI4
  implicit none
  integer, intent(in)                          :: kmesh,mesh,lfp
  real(kind=DP),intent(in), dimension(kmesh,2) :: rh
  real(kind=DP),intent(in), dimension(kmesh)   :: rad
  real(kind=DP),intent(out),dimension(kmesh,2) :: rho
  real(kind=DP),intent(out),dimension(kmesh)   :: fdiff

  real(kind=DP)  :: ovpi4

  ovpi4 = 1.d0/PAI4
  if(lfp == 1) then
     rho(:,1) = rh(:,1)*ovpi4/(rad(:)*rad(:))
     rho(:,2) = rh(:,2)*ovpi4/(rad(:)*rad(:))
     fdiff(:) = rho(:,1) + rho(:,2)
  else
     rho(:,1) = rh(:,1)*ovpi4/(rad(:)*rad(:))
     fdiff(:) = rho(:,1)
  end if
end subroutine cpden

! ------------------------------------------------------
subroutine gdiffs(h,rad,fdiff,kmesh,mmesh,mddiff,nddiff)
!             8th Aug. 1992
!                 by T.Yamasaki
!           revised by T.Yamasaki on 30th Oct. 1994
!
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)          :: kmesh,mmesh,mddiff,nddiff
  real(kind=DP), intent(in)    :: h, rad(kmesh)
  real(kind=DP), intent(inout) :: fdiff(kmesh,0:mddiff)

  integer       :: id,j
  real(kind=DP) :: alpha

  if(nddiff > mddiff) write(6,*) ' *** nddiff > mddiff *** '
  do id = 1, nddiff
     alpha = 1.d0/(dexp(h*id)-1.d0)
     do j = 1, mmesh - id
        fdiff(j,id) = (fdiff(j+1,id-1) - fdiff(j,id-1)) * alpha/rad(j)
     end do
  end do
end subroutine gdiffs

! ------------------------------------------------------
!             8th Aug. 1992
!                 by T.Yamasaki
subroutine gtgrad(h,rad,fdiff,kmesh,mmesh,nddiff,nddriv &
     &     ,npcntr,coeff,mddiff,grad,modeex)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters, only : DP
  implicit none
  integer,      intent(in) :: kmesh,mmesh,nddiff,nddriv,npcntr,mddiff,modeex
  real(kind=DP),intent(in) :: h,rad(kmesh),fdiff(kmesh,nddiff)
  real(kind=DP)            :: coeff(0:mddiff,nddriv:mddiff)
  real(kind=DP),intent(out):: grad(kmesh)

  integer, parameter :: EXECUT = 2

  if(npcntr > nddiff) then
     write(6,*) ' *** npcntr should be smaller than nddiff ***'
     write(6,*) ' *** (npcntr, nddiff) = (',npcntr,',',nddiff,') ***'
     call phase_error_with_msg(6,'npcntr should be smaller than nddiff',__LINE__,__FILE__)
  endif
  if(nddriv == 1) then
     call gcoef1(h,coeff,mddiff,nddiff)
  else if(nddriv == 2) then
     call gcoef2(h,coeff,mddiff,nddiff)
  else
     write(6,*) ' *** nddriv = ', nddriv, ' ***'
     call phase_error_with_msg(6,' *** nddriv should be smaller than 3 ***',__LINE__,__FILE__)
  endif
#ifdef DEBUG_WRITE
  call wdcoef(coeff,mddiff,nddiff,nddriv)
#endif
  if(modeex == EXECUT) then
     call ggrade(coeff,mddiff,nddiff,nddriv,rad,fdiff,kmesh,mmesh,npcntr,grad)
  endif
end subroutine gtgrad

#ifdef DEBUG_WRITE
! --------------------------------------------------------------
subroutine wdcoef(coeff,mddiff,nddiff,nddriv)
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)       :: mddiff,nddiff,nddriv
  real(kind=DP), intent(in) :: coeff(0:mddiff,nddriv:mddiff)

  integer :: i

  write(6,*) ' ++++ nddriv = ', nddriv, ' ++++'
  do i = nddriv, nddiff
     write(6,*) ' ++++ i(alpha) = ', i ,' ++++'
     write(6,'(10f8.4)') (coeff(j,i),j = 0, nddiff)
  end do
end subroutine wdcoef
#endif

! --------------------------------------------------------------
subroutine gcoef1(h,coeff,mddiff,nddiff)
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)      :: mddiff,nddiff
  real(kind=DP),intent(in) :: h
  real(kind=DP),intent(out):: coeff(0:mddiff,mddiff)

  integer       :: i,npd,ip,j,isum,imul,ipp
  real(kind=DP) :: x, y
  coeff(:,1) = 1.d0

#ifdef SX
!cdir skip
#endif
  do npd = 2, nddiff
#ifdef SX
!cdir skip
#endif
     do i = 0, npd-1
        ip = -i
        x = 1.d0
#ifdef SX
!cdir skip
#endif
        do j = 1, npd - 1
           if(ip == 0) ip = 1
           x = x*(1.d0 - dexp(h*ip))
           ip = ip + 1
        end do
        coeff(i,npd) = x
     end do
#ifdef SX
!cdir skip
#endif
     do i = npd, nddiff
        ip = -i
        y = 0.d0
#ifdef SX
!cdir skip
#endif
        do isum = 1, npd
           ipp = ip
           x = 1.d0
#ifdef SX
!cdir skip
#endif
           do imul = 1, npd
              if(imul /= isum)  x = x*(1.d0 - dexp(h*ipp))
!                    x = x*(0 - ipp*h)
              ipp = ipp + 1
           end do
           y = y + x
        end do
        coeff(i, npd) = y
     end do
  end do
end subroutine gcoef1

! --------------------------------------------------------------------
subroutine gcoef2(h,coeff,mddiff,nddiff)
  use m_Const_Parameters, only : DP
  implicit none

  integer, intent(in) :: mddiff, nddiff
  real(kind=DP), intent(in)  :: h
  real(kind=DP), intent(out) :: coeff(0:mddiff,2:mddiff)

  integer       :: i, npd,isum,imul,ip,ipp,jsum
  real(kind=DP) :: x, y

  coeff(:,2) = 2.d0

#ifdef SX
!cdir skip
#endif
  do npd = 3, nddiff
#ifdef SX
!cdir skip
#endif
     do i = 0, npd-1
        ip = -i
        y = 0.d0
#ifdef SX
!cdir skip
#endif
        do isum = 1, npd-1
           ipp = ip
           x = 1.d0
#ifdef SX
!cdir skip
#endif
           do imul = 1, npd-1
              if(ipp == 0) ipp = 1
              if(imul /= isum)   x = x*(1 - dexp(h*ipp))
!                    x = x*(0.d0 - ipp*h)
              ipp = ipp + 1
           end do
           y = y + x
        end do
        coeff(i,npd) = 2.d0*y
     end do
#ifdef SX
!cdir skip
#endif
     do i = npd, nddiff
        ip = -i
        y = 0.d0
#ifdef SX
!cdir skip
#endif
        do isum = 1, npd-1
           do jsum = isum+1, npd
              ipp = ip
              x = 1.d0
#ifdef SX
!cdir skip
#endif
              do imul = 1, npd
                 if(isum /= imul .and. jsum /= imul)  x = x*(1 - dexp(h*ipp))
!                       x = x*(0.d0 - ipp*h)
                 ipp = ipp + 1
              end do
              y = y + x
           end do
        end do
        coeff(i, npd) = 2.d0*y
     end do
  end do
end subroutine gcoef2

! --------------------------------------------------------------------
subroutine ggrade(coeff,mddiff,nddiff,nddriv,rad,fdiff,kmesh,mmesh,npcntr,grad)
  use m_Const_Parameters, only : DP
  implicit none

  integer, intent(in) :: mddiff,nddiff,nddriv,kmesh,mmesh,npcntr
  real(DP),intent(in) :: coeff(0:mddiff,nddriv:mddiff),rad(0:kmesh-1),fdiff(0:kmesh-1,mddiff)
  real(DP),intent(out):: grad(0:kmesh-1)

  integer :: j,ipcntr,idf,iddiff,id

  grad = 0.d0

  ipcntr = -1
  idf  = nddiff - (npcntr+1)*2
  do j = 0, npcntr - 1
     ipcntr = ipcntr + 1
     idf = idf + 2
     if(idf < nddriv+1) then
        iddiff = nddriv + 1
     else
        iddiff = idf
     endif
     do id = iddiff, nddriv+1, -1
        grad(j) = (grad(j) + coeff(ipcntr,id)*fdiff(j-ipcntr,id))*rad(j)
     end do
  end do

  do id = nddiff, nddriv+1, -1
     do j = npcntr, mmesh-1-nddiff+npcntr
        grad(j) = (grad(j) + coeff(npcntr,id)*fdiff(j-npcntr,id))*rad(j)
     end do
  end do

  ipcntr = npcntr
  iddiff = nddiff
  do j = mmesh-nddiff+npcntr,mmesh-1
     iddiff = iddiff - 2
     if(iddiff < nddriv+1) then
        iddiff = nddriv + 1
        ipcntr = ipcntr + 1
     else
        ipcntr = ipcntr - 1
     endif
     do id = iddiff, nddriv + 1, -1
        grad(j) = (grad(j) + coeff(ipcntr,id)*fdiff(j-ipcntr,id))*rad(j)
     end do
  end do

  ipcntr = -1
  idf = nddiff - (npcntr + 1)*2
  do j = 0, npcntr-1
     ipcntr = ipcntr + 1
     idf = idf + 2
     if(idf < nddriv+1) then
        iddiff = nddriv+1
     else
        iddiff = idf
     endif
     grad(j) = grad(j) + coeff(j,nddriv)*fdiff(j-ipcntr,nddriv)
  end do

  do j = npcntr, mmesh-nddiff+npcntr-1
     grad(j) = grad(j) + coeff(npcntr,nddriv)*fdiff(j-npcntr,nddriv)
  end do

  ipcntr = npcntr
  iddiff = nddiff
  do j = mmesh-nddiff+npcntr, mmesh-1
     iddiff = iddiff - 2
     if(iddiff < nddriv+1) then
        iddiff = nddriv + 1
        ipcntr = ipcntr + 1
     else
        ipcntr = ipcntr - 1
     endif
     grad(j) = grad(j) + coeff(ipcntr,nddriv)*fdiff(j-ipcntr,nddriv)
  end do

end subroutine ggrade

! ----------------------------------------
!       14th Aug. 1992 by T.Yamasaki
!    Revised by T.Yamasaki on 30th Oct. 1994
subroutine cnggrd1(grdnts,kmesh,mmesh)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)         :: kmesh,mmesh
  real(kind=DP),intent(inout) :: grdnts(kmesh,1)

  integer       :: i
  real(kind=DP) :: fac

  do i = 1, mmesh
     fac = 1.d0
     if(grdnts(i,1) < 0.d0) fac = -1.d0
     grdnts(i,1) = grdnts(i,1)*grdnts(i,1)
     grdnts(i,1) = dsqrt(grdnts(i,1))*fac
  end do
end subroutine cnggrd1

! -------------------------------------------
!      15TH AUG. 1992
!             BY T.YAMASAKI
!
subroutine cpval(a,mel,nel,cpmod,nmode,b)
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)      :: mel,nel,cpmod,nmode
  real(kind=DP),intent(in) :: a(mel,nmode)
  real(kind=DP),intent(out):: b(mel)

  integer, parameter :: ABSLUT = 1, NORMAL = 2

  if(cpmod == NORMAL) then
     if(nmode == 1) then
        b(1:nel) = a(1:nel,1)
     else if(nmode ==  2) then
        b(1:nel) = a(1:nel,1) + a(1:nel,2)
     endif
  else if(cpmod == ABSLUT) then
     if(nmode == 1) then
        b(1:nel) = dsqrt(a(1:nel,1)*a(1:nel,1))
     else if(nmode == 2) then
        b(1:nel) = dsqrt(a(1:nel,1)*a(1:nel,1) + a(1:nel,2)*a(1:nel,2))
     endif
  else
     write(6,*) ' !D invalid value of (cpmod = ',cpmod,')'
  endif

end subroutine cpval

! -------------------------------------------
!      30TH Nov. 1994
!             BY T.YAMASAKI
!
subroutine cpval_abs(a,mel,nel,nmode,b)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)       :: mel, nel, nmode
  real(kind=DP), intent(in) :: a(mel,nmode)
  real(kind=DP), intent(out):: b(mel)

  if(nmode == 1) then
     b(1:nel) = dsqrt(a(1:nel,1)*a(1:nel,1))
  else if(nmode == 2) then
     b(1:nel) = dsqrt(a(1:nel,1)*a(1:nel,1) + a(1:nel,2)*a(1:nel,2))
  end if
end subroutine cpval_abs

! -------------------------------------------
!      15TH AUG. 1992
!             BY T.YAMASAKI
!      1st DEC. 1994
!           Revised by T.Yamasaki 
!
subroutine getroh(rh,kmesh,mesh,rad,rho)
  use m_Const_Parameters, only : DP,PAI4
  implicit none
  integer, intent(in) :: kmesh,mesh
  real(DP),intent(in) :: rh(kmesh),rad(kmesh)
  real(DP),intent(out):: rho(kmesh)

  real(DP) :: ovpi4

  ovpi4 = 1.d0/PAI4
  rho(1:mesh) = rh(1:mesh)*ovpi4/(rad(1:mesh)*rad(1:mesh))
end subroutine getroh

! -------------------------------------------
!      1ST DEC. 1994
!             BY T.YAMASAKI
!
subroutine getroh2(rh,rhpc,kmesh,mesh,rad,rho)
  use m_Const_Parameters, only : DP,PAI4
  implicit none
  integer, intent(in) :: kmesh,mesh
  real(DP),intent(in) :: rh(kmesh),rhpc(kmesh),rad(kmesh)
  real(DP),intent(out):: rho(kmesh)

  real(DP) :: ovpi4

  ovpi4 = 1.d0/PAI4
  rho(1:mesh) = (rh(1:mesh)+rhpc(1:mesh))*ovpi4/(rad(1:mesh)*rad(1:mesh))
end subroutine getroh2

! ----------------------------------------
!       14TH AUG. 1992 BY T.YAMASAKI
!       1ST DEC.  1994  Revised by T.Yamasaki
!       3rd DEC.  1994  Revised by Y.Morikawa
!
subroutine cnggrd(grdnts,rad,kmesh,mesh)
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)   :: kmesh,mesh
  real(DP),intent(in)   :: rad(kmesh)
  real(DP),intent(inout):: grdnts(kmesh,3)

  integer :: i
!!$cccc  changed by Y.M
  do  i = 1, mesh
     grdnts(i,3) = grdnts(i,2)
     grdnts(i,2) = grdnts(i,2) + 2*grdnts(i,1)/rad(i)
  end do
!!$cccc  changed by Y.M

  grdnts(1:mesh,1) = dabs(grdnts(1:mesh,1))

end subroutine cnggrd

! ----------------------------------------
!       15TH AUG. 1992 BY T.YAMASAKI
!
subroutine mkprdc(a,b,c,m,n,ndeg)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in) :: m,n,ndeg
  real(DP),intent(in) :: a(m,ndeg),b(m,ndeg)
  real(DP),intent(out):: c(m)

  integer :: ideg,ip
  c = 0.d0

  do ideg = 1, ndeg
     do ip = 1, n
        c(ip) = c(ip) + a(ip,ideg)*b(ip,ideg)
     end do
  end do
end subroutine mkprdc

subroutine  epcor_0(nspin,DELTA,tchgr_l,f2or1,epc)
  ! Electron-positron correlation, 0-density limit
  ! Coded by Dr. Mineo Saito, 15th Nov. 2003
  use m_Const_Parameters,  only : PAI,DP
  use m_Parallelization,   only : ista_sfftph, iend_sfftph
  implicit none

 integer,intent(in)        :: nspin
 real(kind=DP),intent(in)  :: DELTA
 real(kind=DP),intent(inout)  :: tchgr_l(ista_sfftph:iend_sfftph,nspin)
 real(kind=DP),intent(in)  :: f2or1(ista_sfftph:iend_sfftph) 
 real(kind=DP),intent(out)  :: epc 
 
 real*8 rs,facw,d,alog,e_density,epc_derivative,epc_energy, HARTREE
 integer i_derive,is,i

  HARTREE = 1.d0/2.d0
  epc=0.d0
  facw = nspin
  do is = 1, nspin
   do i = ista_sfftph, iend_sfftph 
     d  = facw * tchgr_l(i, is)
     rs=1.d0/((4.d0/3.d0*PAI*d)**(1.d0/3.d0))  
     if(rs.le.0.d0) call phase_error_with_msg(6,'error epxc',__LINE__,__FILE__)
     i_derive=0
	 if(rs.lt.0.302d0)then
	  alog=log(rs)
	  epc_energy=-1.56d0/dsqrt(rs)+(0.051d0*alog-0.081d0)*alog+1.14d0
	  epc_derivative=1.56d0/2/rs**1.5d0+(0.051d0*alog-0.081d0)/rs+0.051d0*alog/rs
         elseif(rs.lt.0.56d0)then
	  epc_energy=-0.92305d0-0.05459d0/rs/rs
	  epc_derivative=0.05459*2/rs**3
         elseif(rs.lt.8.0)then
	  epc_energy=-13.15111d0/(rs+2.5d0)**2+2.8655d0/(rs+2.5)-0.6298d0
	  epc_derivative=13.15111*2/(rs+2.5d0)**3-2.8655d0/(rs+2.5d0)**2
	 else
          e_density=1.d0/(PAI*4.d0/3.d0*rs*rs*rs) 
	  epc_energy=-179856.2768*e_density*e_density+186.4207*e_density-0.524d0
	  epc_derivative=-179856.2768*3*e_density**2+186.4207d0*2*e_density-0.524d0
	  i_derive=1
	endif 
!       tchgr_l(i,is)=(epc_energy-rs/3.d0*epc_derivative)*HARTREE
       tchgr_l(i,is)=epc_energy*HARTREE
!	    if(i_derive.eq.1) tchgr_l(i,is)=epc_derivative*HARTREE
!		epc=epc+tchgr_l(i,is)*f2or1(i)*HARTREE
		epc=epc+tchgr_l(i,is)*f2or1(i)
    enddo
  enddo
end subroutine epcor_0  
subroutine  epcor_00(DELTA,tchgr_l,nnmesh)
  ! Electron-positron correlation, 0-density limit
  ! Coded by Dr. Mineo Saito, 15th Nov. 2003
  use m_Const_Parameters,  only : PAI,DP
  implicit none

real(kind=DP),intent(in)  :: DELTA
 real(kind=DP),intent(inout)  :: tchgr_l(nnmesh)
 
real(kind=DP)  :: epc 
real*8 rs,d,alog,e_density,epc_derivative,epc_energy,HARTREE
 integer i_derive,i,nnmesh
 
HARTREE=1.d0/2.d0
  epc=0.d0
  do i = 1, nnmesh 
     d  = tchgr_l(i)
     rs=1.d0/((4.d0/3.d0*PAI*d)**(1.d0/3.d0))  
     if(rs.le.0.d0) call phase_error_with_msg(6,'error epxc',__LINE__,__FILE__)
     i_derive=0
	 if(rs.lt.0.302d0)then
	  alog=log(rs)
	  epc_energy=-1.56d0/dsqrt(rs)+(0.051d0*alog-0.081d0)*alog+1.14d0
	  epc_derivative=1.56d0/2/rs**1.5d0+(0.051d0*alog-0.081d0)/rs+0.051d0*alog/rs
         elseif(rs.lt.0.56d0)then
	  epc_energy=-0.92305d0-0.05459d0/rs/rs
	  epc_derivative=0.05459*2/rs**3
         elseif(rs.lt.8.0)then
	  epc_energy=-13.15111d0/(rs+2.5d0)**2+2.8655d0/(rs+2.5)-0.6298d0
	  epc_derivative=13.15111*2/(rs+2.5d0)**3-2.8655d0/(rs+2.5d0)**2
	 else
          e_density=1.d0/(PAI*4.d0/3.d0*rs*rs*rs) 
	  epc_energy=-179856.2768*e_density*e_density+186.4207*e_density-0.524d0
	  epc_derivative=-179856.2768*3*e_density**2+186.4207d0*2*e_density-0.524d0
	  i_derive=1
	endif 
!       tchgr_l(i)=(epc_energy-rs/3.d0*epc_derivative)*HARTREE
       tchgr_l(i)=epc_energy*HARTREE
!	    if(i_derive.eq.1) tchgr_l(i)=epc_derivative*HARTREE
!		epc=epc+tchgr_l(i)*HARTREE
		epc=epc+tchgr_l(i)
    enddo
 end subroutine epcor_00 
subroutine enhance_0(e_density,enhance)
  use m_Const_Parameters,  only : PAI,DP
  real*8 rs,e_density,DEL,enhance
  DEL=1.d-8
  if(e_density.gt.DEL)then
     rs=1.d0/((4.d0/3.d0*PAI*e_density)**(1.d0/3.d0))  
     enhance=1.d0+1.23d0*rs+0.9889d0*rs**1.5d0 &
          -1.4820d0*rs**2+0.3956d0*rs**2.5d0+rs**3/6.d0
  else
     enhance=0.d0
  endif
end subroutine enhance_0

subroutine enhance_01(e_density,enhance,epsilon_ele)
 use m_Const_Parameters,  only : PAI,DP
 real*8 rs,e_density,DEL,enhance,epsilon_ele
 DEL=1.d-8
!  epsilon=6.34466667d0
 if(e_density.gt.DEL)then
    rs=1.d0/((4.d0/3.d0*PAI*e_density)**(1.d0/3.d0))
!!$    enhance=1.d0+1.23d0*rs+0.8295d0*rs**1.5d0 &
!!$         -1.26d0*rs**2+0.3286d0*rs**2.5d0+(rs**3)/6.d0*(1.d0-1.d0/epsilon_ele)
    enhance=1.d0+1.23d0*rs+0.9889d0*rs**1.5d0 &
         -1.4820d0*rs**2+0.3956d0*rs**2.5d0+(rs**3)/6.d0*(1.d0-1.d0/epsilon_ele)
 else
    enhance=0.d0
 endif
end subroutine enhance_01

subroutine enhance_gga_0(e_density,gradient,enhance)
  use m_Const_Parameters,  only : PAI,DP
  real*8 rs,e_density,DEL,enhance,gradient,fermi_momentum,Thomas_length,enhance_lda,e_param
  DEL=1.d-8
  if(e_density.gt.DEL)then
     rs=1.d0/((4.d0/3.d0*PAI*e_density)**(1.d0/3.d0))  
     enhance_lda=1.d0+1.23d0*rs-0.0742d0*rs**2.d0+rs**3.d0/6.d0
	 fermi_momentum=1.92d0/rs
	 Thomas_length=(4.d0/PAI*fermi_momentum)**(0.5d0)
	 e_param=(gradient/(e_density*Thomas_length))**2.d0
	 enhance=1.d0+(enhance_lda-1.d0)*exp(-0.22d0*e_param)
  else
     enhance=0.d0
  endif
end subroutine enhance_gga_0

subroutine epcor_finite(nspin,DELTA,etchgr_l,ptchgr_l,f2or1,epc)
    use m_Const_Parameters,  only : PAI,DP
    use m_Parallelization,   only : ista_sfftph, iend_sfftph
	implicit none

    integer,intent(in)        :: nspin
    real(kind=DP),intent(in)  :: DELTA
    real(kind=DP),intent(inout)  :: etchgr_l(ista_sfftph:iend_sfftph,nspin)
    real(kind=DP),intent(inout)  :: ptchgr_l(ista_sfftph:iend_sfftph,nspin)
    real(kind=DP),intent(in)  :: f2or1(ista_sfftph:iend_sfftph) 
    real(kind=DP),intent(out)  :: epc
    real*8 facw
    integer is,i
	
	real*8 rs_e,rs_p,epcor_inverse,epcor,ed,pd
	real*8 a_rs,b_rs,c_rs,e_energy_ap,p_energy_ap,p_ap_derivative
	real*8 e_ap_derivative,epc_ne_derivative,epc_np_derivative
	real*8 a_derivative,b_derivative,c_derivative
	epc=0.d0
	facw=nspin
	do is=1,nspin
        do i = ista_sfftph, iend_sfftph 
			ed=facw*etchgr_l(i,is)
			pd=facw*ptchgr_l(i,is)
			rs_e=1.d0/((4.d0/3.d0*PAI*ed)**(1.d0/3.d0))
			rs_p=1.d0/((4.d0/3.d0*PAI*pd)**(1.d0/3.d0))
			if(rs_e.le.0.d0) call phase_error_with_msg(6,'error epxc',__LINE__,__FILE__)
			call epcor_abc(rs_e,a_rs,b_rs,c_rs,a_derivative,b_derivative,c_derivative)
			call cal_energy_ap(rs_e,e_energy_ap,e_ap_derivative)
			call cal_energy_ap(rs_p,p_energy_ap,p_ap_derivative)
			epcor_inverse=a_rs+b_rs*rs_p+c_rs*rs_p**2.d0+(4.d0*PAI/3.d0)*(rs_p**3.d0) &
								/e_energy_ap+(4.d0*PAI/3.d0)*(rs_e**3.d0)/p_energy_ap
			epcor=1.d0/epcor_inverse
			epc_ne_derivative=-(a_derivative+b_derivative*rs_p+c_derivative*rs_p**2.d0 &
								-4.d0*PAI/3.d0*e_ap_derivative*rs_p*rs_p*rs_p/(e_energy_ap &
								*e_energy_ap)+(4.d0*PAI*rs_e*rs_e/p_energy_ap)* &
								(-4.d0*PAI/9.d0*rs_e**4.d0))*(epcor*epcor)
			epc_np_derivative=-((b_rs+2.d0*c_rs*rs_p+4.d0*PAI*rs_p*rs_p/e_energy_ap)* &
								(-4.d0*PAI/9.d0*rs_p**4.d0)-4.d0*PAI/3.d0*p_ap_derivative* &
								rs_e*rs_e*rs_e/(p_energy_ap*p_energy_ap))*(epcor*epcor)
			etchgr_l(i,is)=epc_ne_derivative
			ptchgr_l(i,is)=epc_np_derivative
			epc=epc+epcor*f2or1(i)
		end do
	end do
      end subroutine epcor_finite

subroutine cal_energy_ap(rs,energy_ap,ap_derivative)
    use m_Const_Parameters,  only : PAI
	implicit none
	real*8 rs,log,density,energy_ap,ap_derivative
	integer i_derive
	i_derive=0
	if(rs.lt.0.302d0)then
		log=dlog(rs)
		energy_ap=-1.56d0/dsqrt(rs)+(0.051d0*log-0.081d0)*log+1.14d0
		ap_derivative=1.56d0/2.d0/rs**1.5d0+(0.051d0*log-0.081d0)/rs &
						+0.051d0*log/rs
	elseif(rs.lt.0.56d0)then
		energy_ap=-0.92305d0-0.05459d0/rs/rs
		ap_derivative=0.05459*2.d0/rs**3.d0
	elseif(rs.lt.8.0)then
		energy_ap=-13.15111d0/(rs+2.5d0)**2+2.8655d0/(rs+2.5)-0.6298d0
		ap_derivative=13.15111*2.d0/(rs+2.5d0)**3.d0-2.8655d0/(rs+2.5d0)**2.d0
	else
		density=3.d0/(4.d0*PAI*rs**3.d0)
		energy_ap=-179856.2768*density*density+186.4207*density-0.524d0
		ap_derivative=-179856.2768*2.d0*density+186.4207d0
		i_derive=1
	end if
	if(i_derive.eq.0)   ap_derivative=ap_derivative*(-4.d0*PAI/9.d0*rs**4.d0)
end subroutine cal_energy_ap

subroutine epcor_abc(rs_e,a_rs,b_rs,c_rs,a_derivative,b_derivative,c_derivative)
    use m_Const_Parameters,  only : PAI
	implicit none
	real*8 rs_e,a_rs,b_rs,c_rs,a_derivative,b_derivative,c_derivative
	a_rs=69.7029d0-107.4927d0*rs_e+23.7182d0*rs_e**2.d0
	a_derivative=(-107.4927d0+23.7182d0*2.d0*rs_e)*(-4.d0*PAI/9.d0*rs_e**4.d0)
	b_rs=-107.4927d0+141.8458d0*rs_e-33.6472d0*rs_e**2.d0
	b_derivative=(141.8458d0-33.6472d0*2.d0*rs_e)*(-4.d0*PAI/9.d0*rs_e**4.d0)
	c_rs=23.7182d0-33.6472d0*rs_e+5.21152d0*rs_e**2.d0
	c_derivative=(-33.6472d0+5.21152d0*2.d0*rs_e)*(-4.d0*PAI/9.d0*rs_e**4.d0)
end subroutine epcor_abc

!$$#ifdef PARA3D

!===============================================================================

!!====================================================================
!!        modified  by T.Yamasaki  94/11/14
!!
!!        Modified for Spin polrarized GGA by K. Kato 1995/1/18 
!!
subroutine ex_ggapw91_3D(nspin,ispin,chgrhr_l,grad_rho,f2or1,exc,dFx_drho,dFx_dgradrho,nfft_y)
!!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : PAI,DP
  implicit none

  integer,intent(in)        :: nspin,ispin,nfft_y
  real(kind=DP),intent(in)  :: chgrhr_l(1:nfft_y,ispin)
  real(kind=DP),intent(in)  :: grad_rho(1:nfft_y,nspin)
  real(kind=DP),intent(in)  :: f2or1(1:nfft_y)  
  real(kind=DP),intent(out) :: exc
  real(kind=DP),intent(out) :: dFx_drho(1:nfft_y,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(1:nfft_y,nspin)

  integer :: is, i
  real(kind=DP) :: exc1,facw,d,dd,fk,s,fac,s2,s3,s4&
       & ,p0,p1,p2,p3,p4,p5,p6,f,ex,fs,exd,exdd,exc0,excd,excdd

  real(kind=DP), parameter :: a1 = 0.19645d0
  real(kind=DP), parameter :: a2 = 0.27430d0
  real(kind=DP), parameter :: a3 = 0.15084d0
  real(kind=DP), parameter :: a4 = 100.0d0

! === KT_mod === 2014/11/25
!!!  real(kind=DP), parameter :: ax = -0.7385588d0
       real(kind=DP), parameter ::  ax = -0.7385587663820224d0
!=============== 2014/11/25

  real(kind=DP), parameter :: a  =  7.7956d0
  real(kind=DP), parameter :: b1 =  0.004d0

! === KT_mod === 2014/11/25
!  real(kind=DP), parameter :: thrd   = 0.33333333333d0
!  real(kind=DP), parameter :: thrd4  = 1.333333333333333d0
  real(kind=DP), parameter :: thrd   = 0.33333333333333333333333d0
  real(kind=DP), parameter :: thrd4  = 1.33333333333333333333333d0
!=============== 2014/11/25

  real(kind=DP), parameter :: thpith = 3.0936677262801d0

!---- Spin dependency
#ifdef __TIMER_SUB__
    call timer_sta(764)
#endif

  facw = ispin
  exc = 0.d0
#ifdef NEC_TUNE_MXCP
!CDIR NOCONCUR
#endif
#ifdef __TIMER_DO__
  call timer_sta(898)
#endif
  do is = 1, ispin
     exc1 = 0.d0
#ifdef NEC_TUNE_MXCP
!CDIR INNER
#endif
     do i = 1, nfft_y       ! MPI
        d  = facw * chgrhr_l(i, is)
        dd = facw * grad_rho(i, is)
        fk = (3*PAI*PAI*d)**thrd
        if(d > 1.d-05) then
           s = dd/(d*fk*2)
        else
           s = 0.d0
        endif

        fac = ax*d**thrd
        s2  = s*s
        s3  = s2*s
        s4  = s2*s2
        p0  = 1.d0/(dsqrt(1+a*a*s2))
        p1  = dlog(a*s+1/p0)
        p2  = dexp(-a4*s2)
        p3  = 1.d0/(1+a1*s*p1+b1*s4)
        p4  = 1 + a1*s*p1 + (a2-a3*p2)*s2
        f   = p3*p4
        ex  = fac*f*d
        p5  = 2*(s*(a2-a3*p2)+a3*a4*s3*p2 - 2*b1*s3)
        p6  = (a1*(p1+a*s*p0)+4*b1*s3)*((a2-a3*p2)*s2-b1*s4)
        fs  = (p5*p3-p6*p3*p3)
        exd = thrd4*fac*(f-s*fs)
        exdd = ax*fs*0.5/thpith

        exc0 = ex / facw
        excd = exd
        dFx_drho(i, is) = excd 
        excdd = exdd
!!$c gradient of charge density 95/12/3 H.S.
        if(dabs(grad_rho(i, is)) > 1.d-9) then
           dFx_dgradrho(i, is) = excdd / grad_rho(i, is)
        else
           dFx_dgradrho(i, is) = 0.d0
        endif
        exc1 = exc1 + exc0*f2or1(i)
     end do
     exc = exc + exc1
  end do
#ifdef __TIMER_DO__
  call timer_end(898)
#endif
#ifdef __TIMER_SUB__
    call timer_end(764)
#endif
end subroutine ex_ggapw91_3D

!===============================================================================

!!===================================================================
!!        modified  by T.Yamasaki  94/11/15
!!
!!        The original is modified and Spin polrarized GGA programs
!!        are added in this subroutine
!!                                                      by K. Kato 1995/1/18 
!!
subroutine cr_ggapw91_3D(nspin,ispin,chgrhr_l,grad_trho,f2or1,exc,dF_drho,nfft_y)
!!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : PAI,DP
  implicit none

  integer,intent(in)        :: nspin,ispin,nfft_y
  real(kind=DP),intent(in)  :: chgrhr_l(1:nfft_y,ispin)
  real(kind=DP),intent(inout) :: grad_trho(1:nfft_y)
  real(kind=DP),intent(in)  :: f2or1(1:nfft_y)  
  real(kind=DP),intent(inout) :: exc
  real(kind=DP),intent(inout) :: dF_drho(1:nfft_y,nspin)

  real(kind=DP), parameter :: zeta_minimum    = 1.d-4
  real(kind=DP), parameter :: zeta_minimum2   = 1.d-6
  real(kind=DP), parameter :: density_minimum = 1.d-20
  real(kind=DP), parameter :: density_minimum2= 1.d-5

  real(kind=DP), parameter :: a  = 0.0310907d0,     a1   = 0.21370d0
  real(kind=DP), parameter :: b1 = 7.5957d0,        b2n  = 3.5876d0
  real(kind=DP), parameter :: b3 = 1.6382d0,        b4   = 0.49294d0
  real(kind=DP), parameter :: p  = 1.00d0,          p1   = p+1.d0
  real(kind=DP), parameter :: ap  = 0.015545d0,     a1p  = 0.20548d0
  real(kind=DP), parameter :: b1p = 14.1189d0,      b2np = 6.1977d0
  real(kind=DP), parameter :: b3p = 3.3662d0,       b4p = 0.62517d0
! vwn 95/12/2 Y.M
  real(kind=DP), parameter :: aq  = 0.016887d0,     a1q = 0.11125d0
  real(kind=DP), parameter :: b1q = 10.3570d0,      b2nq = 3.6231d0
  real(kind=DP), parameter :: b3q = 0.88026d0,      b4q = 0.49671d0
! vwn 95/12/2 Y.M

  real(kind=DP), parameter :: xnu = 15.75592d0,     cc0 =  0.004235d0
  real(kind=DP), parameter :: cx  = -0.001667212d0, alf =  0.09d0
  real(kind=DP), parameter :: c1  =  0.002568d0,    c2  =  0.023266d0
  real(kind=DP), parameter :: c3  =  7.389d-6,      c4  =  8.723d0
  real(kind=DP), parameter :: c5  =  0.472d0,       c6  =  7.389d-2
  real(kind=DP), parameter :: a4  = 100.d0

! ==== KT_mod ==== 2014/11/25
!  real(kind=DP), parameter :: thrd   = 0.333333333333d0
!  real(kind=DP), parameter :: sixth7 = 1.1666666666666666d0
  real(kind=DP), parameter :: thrd   = 0.3333333333333333333333d0
  real(kind=DP), parameter :: sixth7 = 1.1666666666666666666666d0
! ================= 2014/11/25

  integer       :: is,i
  real(kind=DP) :: facw,bet,delt,exc1,g,g3,g4,facpon,d,dd,rs,fk,sk,t,s&
       & , q0,rs12,rs32,rsp,q1,q2,eu,q3,eurs,ecd,pon,b,b2,t2,s2,t4,t6 &
       & , rs2,rs3,q4,q5,q6,q7,cc,r0,r1,coeff,r2,r3,a4ms,h0,h1,h,q8 &
       & , h0t,h0b,h0rs,h1t,ccrs,r1rs,h1rs,ht,hrs &
       & , ec1,ec1d,ec1dd,exc0,excd,excdd &
       & , thrd2,thrd4,zeta,q0p,q1p,q2p,q3p,q0q,q1q,q2q,q3q,fzd,fzdd,fz&
       & , euzt,zetadxd,bzt,h0zt,h1zt,gzt &
       & , onpzeta, onmzeta, onzt4 

#ifdef __TIMER_SUB__
    call timer_sta(765)
#endif
  facw = ispin

! === KT_mod === 2014/11/25
!  bet = xnu*cc0
  bet=0.06672455060314922d0
! ============== 2014/11/25

  delt = 2.d0*alf/bet

!---- Spin dependency

#ifdef NEC_TUNE_MXCP
!CDIR NOCONCUR
#endif
#ifdef __TIMER_DO__
  call timer_sta(899)
#endif
  do is=1, ispin
     exc1 = 0.d0
     if ( ispin == 1 ) then
        g = 1.d0
        g3  = g**3
        g4  = g3*g
        facpon = -delt/(g3*bet)

#ifdef NEC_TUNE_MXCP
!CDIR INNER
#endif
        do i = 1, nfft_y  ! MPI
           d = facw*chgrhr_l(i, 1)
           if(d < 1.d-20) cycle
! gradient of charge density 95/12/3 H.S.
           dd =facw* grad_trho(i)
           rs = (0.75/(PAI*d))**thrd
           fk = (3*PAI*PAI*d)**thrd
           sk = dsqrt(4*fk/PAI)
           if(d > 1.d-05) then
              t = dd/(d*sk*2)
              s = dd/(d*fk*2)
           else
              t = 0.d0
              s = 0.d0
           endif

           q0 = -2*a*(1+a1*rs)
           rs12 = dsqrt(rs)
           rs32 = rs12**3
           rsp = rs**p
           q1 = 2.d0*a*(b1*rs12+b2n*rs+b3*rs32+b4*rs*rsp)
           q2 = log(1 + 1/q1)
           eu = q0*q2
           q3 = a*(b1/rs12+2*b2n+3*b3*rs12+2*b4*p1*rsp)
           eurs = -2*a*a1*q2 - q0*q3/(q1**2 + q1)
           ecd = eu-thrd*rs*eurs
           pon = facpon*eu
           b   = delt/(dexp(pon)-1.d0)
           b2  = b*b
           t2  = t*t
           s2  = s*s
           t4  = t2*t2
           t6  = t4*t2
           rs2 = rs*rs
           rs3 = rs2*rs
           q4  = 1.d0+b*t2
           q5  = 1.d0+b*t2+b2*t4
           q6  = c1+c2*rs+c3*rs2
           q7  = 1.d0+c4*rs+c5*rs2+c6*rs3
           cc  = -cx + q6/q7
           r0  = (sk/fk)**2
           r1  = a4*r0*g4
           coeff = cc-cc0-3*cx/7
           r2  = xnu*coeff*g3
           a4ms = -a4*g4*s2
           r3  = dexp(a4ms)
           h0  = g3*(bet/delt)*dlog(1+delt*q4*t2/q5)
           h1  = r3*r2*t2
           h   = h0+h1
!======================================================
           q8  = q5*q5+delt*q4*q5*t2
           h0t = 2.d0*bet*t*(1.d0+2.d0*b*t2)/q8 * g3
           h0b = -bet*t6*(2.d0*b+b2*t2)/q8 * g3
           h0rs = h0b*b*eurs*(b+delt)/ bet / g3
           h1t  = 2.d0*r3*r2*t*(1.d0-r1*t2)
           ccrs = (c2+2*c3*rs)/q7-q6*(c4+2*c5*rs+3*c6*rs2)/q7**2
           r1rs = 100*r0/rs
           h1rs = xnu*t2*r3*(ccrs - coeff*t2*r1rs) * g3
           ht   = h0t+h1t
           hrs  = h0rs + h1rs
           ec1  = d*h/facw
           ec1d = h-thrd*rs*hrs-sixth7*t*ht
           ec1dd = 0.5*ht/sk

!---------------------------------------------------
           exc0 = eu*d/facw + ec1
           excd = ecd + ec1d
           dF_drho(i, is) = dF_drho(i, is) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(dabs(grad_trho(i)) > 1.0d-9) then
              grad_trho(i) = excdd / grad_trho(i)
           else
              grad_trho(i) = 0.d0
           endif
!!$           exec1 = eu*d/facw
           exc1 = exc1 + exc0*f2or1(i)
        end do
     else if ( ispin  == 2 ) then
        thrd2 = thrd * 2.0d0
        thrd4 = thrd * 4.0d0
        fzdd = 8/(9*(2**thrd4-2))

        do i = 1, nfft_y      ! MPI
           d = chgrhr_l(i, 1) + chgrhr_l(i, ispin)
           if(d.lt.1.d-20) cycle
           zeta = ( chgrhr_l(i, 1) - chgrhr_l(i, ispin) ) / d
           onpzeta = 2*chgrhr_l(i,1)     / d
           onmzeta = 2*chgrhr_l(i,ispin) / d
           g  = ( onpzeta**thrd2 + onmzeta**thrd2) * 0.5
!!$           g  = ( (1+zeta)**thrd2 + (1-zeta)**thrd2 ) * 0.5
           if(dabs(zeta) < zeta_minimum) then
              fz = 4.d0/9 * zeta*zeta*(1 + 5.d0/54 * zeta*zeta&
                   &                 *(1 +44.d0/133* zeta*zeta ))&
                   &    /(2**thrd4-2)
           else
              fz = ( (1+zeta)**thrd4+(1-zeta)**thrd4-2)/(2**thrd4-2)
           end if
           g3  = g**3
           g4  = g3*g
           facpon = -delt/(g3*bet)
           if(d < density_minimum) cycle
! gradient of charge density 95/12/3 H.S.
           dd = grad_trho(i)
!!$c        dd = grad_rho(i, 1) + grad_rho(i, 2)
! gradient of charge density 95/12/3 H.S.
           rs = (0.75d0/(PAI*d))**thrd
           fk = (3*PAI*PAI*d)**thrd
           sk = dsqrt(4*fk/PAI)

           if(d > density_minimum2) then
              t = dd/(d*sk*2.d0)
              s = dd/(d*fk*2.d0)
           else
              t = 0.d0
              s = 0.d0
           endif

           q0 = -2*a*(1+a1*rs)
           rs12 = dsqrt(rs)
           rs32 = rs12**3
           rsp = rs**p
           q1 = 2*a*(b1*rs12+b2n*rs+b3*rs32+b4*rs*rsp)
           q2 = log(1+1/q1)
           q3 = a*(b1/rs12+2*b2n+3*b3*rs12+2*b4*p1*rsp)

           q0p = -2*ap*(1+a1p*rs)
           q1p =  2*ap*(b1p*rs12+b2np*rs+b3p*rs32+b4p*rs*rsp)
           q2p = log(1+1/q1p)
           q3p = ap*(b1p/rs12+2*b2np+3*b3p*rs12+2*b4p*p1*rsp)
! vwn 95/12/2 Y.M
           q0q = -2*aq*(1+a1q*rs)
           q1q =  2*aq*(b1q*rs12+b2nq*rs+b3q*rs32+b4q*rs*rsp)
           q2q = log(1+1/q1q)
           q3q = aq*(b1q/rs12+2*b2nq+3*b3q*rs12+2*b4q*p1*rsp)
           if(dabs(zeta) < zeta_minimum) then
              fzd = 8/(9*(2**thrd4 - 2 )) &
                   &     *zeta *( 1 + 5.d0/27.d0 *zeta*zeta &
                   &           *( 1 +22.d0/45.d0 *zeta*zeta ))
           else
              fzd = (4.d0/3)/(2**thrd4-2)*(onpzeta**thrd-onmzeta**thrd )
           end if
           onzt4 = onmzeta*onpzeta * (1 + zeta*zeta)
           eu = q0*q2 + ( q0p * q2p - q0 * q2 ) * fz * zeta**4 &
                &        -   q0q * q2q * fz / fzdd * onzt4
           eurs = (-2*a*a1*q2 - q0*q3/(q1**2+q1)) * (1-fz*zeta**4) &
                &+(-2*ap*a1p*q2p-q0p*q3p/(q1p**2+q1p))*fz*zeta**4 &
                &-(-2*aq*a1q*q2q-q0q*q3q/(q1q**2+q1q))*fz/fzdd*onzt4
           euzt = -q0q*q2q/fzdd*(fzd*(1-zeta**4)-4*fz*zeta**3) &
                &+(q0p*q2p-q0*q2)*(fzd*zeta +4*fz )*zeta**3
! vwn 95/12/2 Y.M
           zetadxd = - 2 * chgrhr_l(i, 3-is) * (-1.d0)**is / d
           ecd = eu - thrd*rs*eurs + euzt * zetadxd

!           if(onpzeta < zeta_minimum2 .or. onmzeta < zeta_minimum2) then
!              ec1 = 0.d0; ec1d = 0.d0; ec1dd= 0.d0
!              goto 10
!           end if
           pon = facpon*eu
           b   = delt/(dexp(pon)-1.d0)
           b2  = b*b
           t2  = t*t
           s2  = s*s
           t4  = t2*t2
           t6  = t4*t2
           rs2 = rs*rs
           rs3 = rs2*rs
           q4  = 1.d0+b*t2
           q5  = 1.d0+b*t2+b2*t4
           q6  = c1+c2*rs+c3*rs2
           q7  = 1.d0+c4*rs+c5*rs2+c6*rs3
           cc  = -cx + q6/q7
           r0  = (sk/fk)**2
           r1  = a4*r0*g4
           coeff = cc-cc0-3.d0*cx/7.d0
           r2  = xnu*coeff*g3
           a4ms = -a4*g4*s2
           r3  = dexp(a4ms)
           h0  = g3*(bet/delt)*dlog(1.d0+delt*q4*t2/q5)
           h1  = r3*r2*t2
           h   = h0+h1

           q8  = q5*q5+delt*q4*q5*t2
           h0t = 2*bet*t*(1+2*b*t2)/q8 * g3
           h0b = -bet*t6*(2*b+b2*t2)/q8 * g3
           h0rs = h0b * b2 * eurs / bet / g3
           h1t  = 2*r3*r2*t*(1.d0-r1*t2)
           ccrs = (c2+2*c3*rs)/q7-q6*(c4+2*c5*rs+3*c6*rs2)/q7**2
           r1rs = 100*r0/rs
           h1rs = xnu*t2*r3*(ccrs - coeff*t2*r1rs) * g3
           ht   = h0t+h1t
           hrs  = h0rs + h1rs

           if(dabs(zeta) < zeta_minimum) then
              gzt = -2.d0/9*zeta *( 1+14.d0/27 *zeta*zeta &
                   &             *( 1+13.d0/18 *zeta*zeta ))
           else
              gzt = ( 1.d0/onpzeta**thrd -1.d0/onmzeta**thrd )/3.d0
           end if
           bzt = pon*b**2 * (euzt-3*eu*gzt/g) / ( bet * g3 ) 
           h0zt = 3*gzt*h0/g + h0b*bzt
           h1zt =(3*gzt + 4*gzt*a4ms ) * h1 / g
           
           ec1  = d*h * 0.5
           ec1d = h-thrd*rs*hrs-sixth7*t*ht + (h0zt+h1zt) * zetadxd
           ec1dd = 0.5*ht/sk

10         continue
           exc0 = eu*d * 0.5 + ec1
           excd = ecd + ec1d
           dF_drho(i, is) = dF_drho(i, is) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(is == 2) then
              if(dabs(grad_trho(i)) > 1.0d-9) then
                 grad_trho(i) = excdd / grad_trho(i)
              else
                 grad_trho(i) = 0.d0
              endif
           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d * 0.5d0
           exc1 = exc1 + exc0*f2or1(i)
        end do
     end if
     exc = exc + exc1
  end do
#ifdef __TIMER_DO__
  call timer_end(899)
#endif
#ifdef __TIMER_SUB__
    call timer_end(765)
#endif
end subroutine cr_ggapw91_3D

!===============================================================================

!!$====================================================================
!!$        modified  by T.Yamasaki  94/11/14
!!$
!!$        Modified for Spin polrarized GGA by K. Kato 1995/1/18 
!!$
!!$        Modified to Perdew-Burke-Ernzerhof(PBE) GGA by K. Kato 1999/4/1 
!!$
subroutine ex_ggapbe_3D(nspin,ispin,chgrhr_l,grad_rho,f2or1,exc,dFx_drho,dFx_dgradrho,nfft_y,iteration,revPBE)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,PAI
  implicit none

  integer,intent(in)        :: nspin,ispin,nfft_y,iteration
  real(kind=DP),intent(in)  :: chgrhr_l(1:nfft_y,ispin)
  real(kind=DP),intent(in)  :: grad_rho(1:nfft_y,nspin)
  real(kind=DP),intent(in)  :: f2or1(1:nfft_y)  
  real(kind=DP),intent(out) :: exc
  real(kind=DP),intent(out) :: dFx_drho(1:nfft_y,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(1:nfft_y,nspin)
  logical, intent(in),optional :: revPBE

#ifndef PREV_EX_GGAPBE
  real(kind=DP), parameter :: ax = -0.7385587663820224d0 !! ax = -(3/4)*(3/PI)^(1/3)
#else
  real(kind=DP), parameter :: ax = -0.7385588d0
#endif
  real(kind=DP), parameter :: thrd = 0.33333333333d0, thrd4 = 1.333333333333333d0
  real(kind=DP), parameter :: thpith = 3.0936677262801d0
#ifndef PREV_EX_GGAPBE
!  real(kind=DP), parameter :: cupa = 0.8040d0, yum = 0.2195149727645171d0
  real(kind=DP), parameter :: yum = 0.2195149727645171d0
#else
  real(kind=DP), parameter :: yum = 0.21951
#endif
  real(kind=DP) :: cupa
  real(kind=DP) :: facw,d,dd,fk,s,fac,s2,yums2,f,ex,fs,exd,exdd,exc0,excd,excdd,exc1
  real(kind=DP) :: tmp
  integer       :: is,i
  logical       :: repbe

!---- Spin dependency
  repbe = .false.
  if(present(revPBE)) repbe = revPBE
  cupa = 0.8040d0
  if(repbe) cupa = 1.245d0

#ifdef __TIMER_SUB__
    call timer_sta(766)
#endif
  facw = ispin
  exc  = 0.d0
#ifdef __TIMER_DO__
  call timer_sta(900)
#endif
  do is = 1, ispin
     exc1 = 0.d0
     do i = 1, nfft_y
        d  = facw * chgrhr_l(i, is)
        dd = facw * grad_rho(i, is)
        fk = (3*PAI*PAI*d)**thrd
        if(d > 1.d-05) then
           s = dd/(d*fk*2)
        else
           s = 0.d0
        endif
!-------------------------------------
!       fac = ax*d**thrd
        tmp = d**thrd
        fac = ax*tmp
        s2  = s*s
        yums2 = yum * s2
        f = 1 + cupa * yums2 / ( cupa + yums2 )
        ex  = fac*f*d
        fs = 2 * cupa ** 2 * yum * s / ( cupa + yums2 ) ** 2
        exd = thrd4*fac*(f-s*fs)
        exdd = ax*fs*0.5d0/thpith
!------------------------------------------     
        exc0 = ex / facw
        excd = exd
        dFx_drho(i, is) = excd 
        excdd = exdd
! gradient of charge density 95/12/3 H.S.
        if(dabs(grad_rho(i, is)) > 1.d-9) then
           dFx_dgradrho(i, is) = excdd / grad_rho(i, is)
        else
           dFx_dgradrho(i, is) = 0.d0
        endif
! gradient of charge density 95/12/3 H.S.
        exc1 = exc1 + exc0*f2or1(i)
     end do
     exc = exc + exc1
  end do
#ifdef __TIMER_DO__
  call timer_end(900)
#endif
#ifdef __TIMER_SUB__
    call timer_end(766)
#endif
end subroutine ex_ggapbe_3D

!===============================================================================

!!$===================================================================
!!$        modified  by T.Yamasaki  94/11/15
!!$
!!$        The original is modified and Spin polrarized GGA programs
!!$        are added in this subroutine
!!$                                                      by K. Kato 1995/1/18 
!!$
!!$        Modified to Perdew-Burke-Ernzerhof(PBE) GGA by K. Kato 1999/4/1 
!!$
subroutine cr_ggapbe_3D(nspin,ispin,chgrhr_l,grad_trho,f2or1,exc,dF_drho,nfft_y)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : PAI,DP
  implicit none

  integer,intent(in)        :: nspin,ispin,nfft_y
  real(kind=DP),intent(in)  :: chgrhr_l(1:nfft_y,ispin)
  real(kind=DP),intent(inout) :: grad_trho(1:nfft_y)
  real(kind=DP),intent(in)  :: f2or1(1:nfft_y)  
  real(kind=DP),intent(inout) :: exc
  real(kind=DP),intent(inout) :: dF_drho(1:nfft_y,nspin)

  real(kind=DP), parameter :: zeta_minimum    = 1.d-4
  real(kind=DP), parameter :: zeta_minimum2   = 1.d-6
  real(kind=DP), parameter :: density_minimum = 1.d-20
  real(kind=DP), parameter :: density_minimum2= 1.d-5

  real(kind=DP), parameter  :: a  = 0.0310907d0,    a1   = 0.21370d0
  real(kind=DP), parameter  :: b1 = 7.5957d0,       b2n  = 3.5876d0
  real(kind=DP), parameter  :: b3 = 1.6382d0,       b4   = 0.49294d0
  real(kind=DP), parameter  :: p  = 1.00d0,         p1   = p+1.d0
#ifndef PREV_CR_GGAPBE
  real(kind=DP), parameter  :: ap  = 0.01554535d0,    a1p  = 0.20548d0
#else
  real(kind=DP), parameter  :: ap  = 0.015545d0,    a1p  = 0.20548d0
#endif
  real(kind=DP), parameter  :: b1p = 14.1189d0,     b2np = 6.1977d0
  real(kind=DP), parameter  :: b3p = 3.3662d0,      b4p  = 0.62517d0
! vwn 95/12/2 Y.M
#ifndef PREV_CR_GGAPBE
  real(kind=DP), parameter  :: aq  = 0.0168869d0,    a1q  = 0.11125d0
#else
  real(kind=DP), parameter  :: aq  = 0.016887d0,    a1q  = 0.11125d0
#endif
  real(kind=DP), parameter  :: b1q = 10.3570d0,     b2nq = 3.6231d0
  real(kind=DP), parameter  :: b3q = 0.88026d0,     b4q  = 0.49671d0
! vwn 95/12/2 Y.M
#ifndef PREV_CR_GGAPBE
  real(kind=DP), parameter  :: gamma = 0.03109069086965489503494086371273d0
#else
  real(kind=DP), parameter  :: gamma = 0.031091
#endif
!--------------------------
#ifndef PREV_CR_GGAPBE
  real(kind=DP), parameter  :: bet=0.06672455060314922d0
  real(kind=DP), parameter  :: delt=bet/gamma
#else
  real(kind=DP), parameter  :: xnu = 15.75592d0, cc0 = 0.004235d0
#endif
  real(kind=DP), parameter  :: thrd = 0.333333333333d0, sixth7 = 1.1666666666666666d0

#ifndef PREV_CR_GGAPBE
!Fix! T.Yamatomto 2011/06/29
  real(kind=DP), parameter  :: eta = 1.d-12
  real(kind=DP), parameter  :: sixthm = -0.16666666666666d0
!Fix! T.Yamatomto 2011/06/29
#endif

  integer       :: is,i
#ifndef PREV_CR_GGAPBE
  real(kind=DP) :: facw,exc1,g,g3,facpon,d,dd,rs,fk,sk,t,q0 &
#else
  real(kind=DP) :: facw,bet,delt,exc1,g,g3,facpon,d,dd,rs,fk,sk,t,q0 &
#endif
       &         , rs12,rs32,rsp,q1,q2,eu,q3,eurs,ecd,pon &
       &         , b,b2,t2,t4,t6,q4,q5,h0,h,q8,h0t,h0b,h0rs,ht,hrs &
       &         , ec1,ec1d,ec1dd,exc0,excd,excdd &
       &         , thrd2,thrd4,fzdd,zeta,onpzeta,onmzeta,fz,zetadxd &
       &         , q0p,q1p,q2p,q3p,q0q,q1q,q2q,q3q,fzd,onzt4,gzt,bzt,h0zt,euzt

#ifdef __TIMER_SUB__
    call timer_sta(767)
#endif
  facw = ispin
#ifdef PREV_CR_GGAPBE
  bet = xnu*cc0
  delt = bet / gamma
#endif
!---- Spin dependency
#ifdef NEC_TUNE_MXCP
!CDIR NOCONCUR
#endif
#ifdef __TIMER_DO__
  call timer_sta(1001)
#endif
  do is = 1, ispin
     exc1 = 0.d0
     if ( ispin == 1 ) then
        g   = 1.d0
        g3  = g**3
        facpon = -delt/(g3*bet)

        do i = 1,nfft_y ! MPI
           d  = facw*chgrhr_l(i, 1)
           if(d < density_minimum) cycle
! gradient of charge density 95/12/3 H.S.
           dd = facw* grad_trho(i)
#ifndef PREV_CR_GGAPBE
           rs = (0.75d0/(PAI*d))**thrd
#else
           rs = (0.75/(PAI*d))**thrd
#endif
           fk = (3*PAI*PAI*d)**thrd
           sk = dsqrt(4*fk/PAI)
           if(d > density_minimum2) then
              t = dd/(d*sk*2)
           else
              t = 0.d0
           endif

           q0 = -2*a*(1 + a1*rs)
           rs12 = dsqrt(rs)
           rs32 = rs12**3
           rsp = rs**p
           q1 = 2*a*(b1*rs12+b2n*rs+b3*rs32+b4*rs*rsp)
           q2 = log(1 + 1.d0/q1)
           eu = q0*q2
           q3 = a*(b1/rs12+2*b2n+3*b3*rs12+2*b4*p1*rsp)
           eurs = -2*a*a1*q2 - q0*q3/(q1**2 + q1)
           ecd = eu-thrd*rs*eurs
           pon = facpon*eu
           b   = delt/(dexp(pon)-1)
           b2  = b*b
           t2  = t*t
           t4  = t2*t2
           t6  = t4*t2
           q4  = 1 + b*t2
           q5  = 1 + b*t2+b2*t4
           h0  = g3*(bet/delt)*dlog(1+delt*q4*t2/q5)
           h   = h0
           q8  = q5*q5+delt*q4*q5*t2
           h0t = 2*bet*t*(1+2*b*t2)/q8 * g3
           h0b = -bet*t6*(2*b+b2*t2)/q8 * g3
           h0rs = h0b*b*eurs*(b+delt)/ bet / g3
           ht   = h0t
           hrs  = h0rs
           ec1  = d*h/facw
           ec1d = h-thrd*rs*hrs-sixth7*t*ht
#ifndef PREV_CR_GGAPBE
           ec1dd = 0.5d0*ht/sk
#else
           ec1dd = 0.5*ht/sk
#endif
           exc0 = eu*d/facw + ec1
           excd = ecd + ec1d
           dF_drho(i, is) = dF_drho(i, is) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(dabs(grad_trho(i)) > 1.0d-9) then
              grad_trho(i) = excdd / grad_trho(i)
           else
              grad_trho(i) = 0.d0
           endif
! gradient of charge density 95/12/3 H.S.
           exc1 = exc1 + exc0*f2or1(i)
        end do
     else if ( ispin ==  2 ) then
        thrd2 = thrd * 2
        thrd4 = thrd * 4
        fzdd = 8/(9*(2**thrd4 - 2 ))

        do i = 1,nfft_y
           d = chgrhr_l(i, 1) + chgrhr_l(i, ispin)
           if(d < density_minimum ) cycle
           zeta = ( chgrhr_l(i, 1) - chgrhr_l(i, 2) ) / d
           onpzeta = 2*chgrhr_l(i, 1) / d
           onmzeta = 2*chgrhr_l(i, 2) / d
           g = ( onpzeta**thrd2 + onmzeta**thrd2 )/2
           if(dabs(zeta) < zeta_minimum ) then
              fz= 4.d0/9 *zeta*zeta*(1 + 5.d0/54 *zeta*zeta &
                   &               *(1 +44.d0/135*zeta*zeta )) &
                   &       /(2**thrd4 -2)
           else
              fz = (onpzeta**thrd4 + onmzeta**thrd4 -2) / (2**thrd4 - 2)
           end if
           g3  = g**3
           facpon = -delt/(g3*bet)
           if(d < density_minimum) cycle
           dd = grad_trho(i)
           rs = (0.75d0/(PAI*d))**thrd
           fk = (3*PAI*PAI*d)**thrd
           sk = dsqrt(4*fk/PAI)

           if(d > density_minimum2) then
#ifndef PREV_CR_GGAPBE
              t = dd/(d*sk*2*g)
#else
              t = dd/(d*sk*2)
#endif
           else
              t = 0.d0
           endif

           q0 = -2*a*(1+a1*rs)
           rs12 = dsqrt(rs)
           rs32 = rs12**3
           rsp = rs**p
           q1 = 2*a*(b1*rs12+b2n*rs+b3*rs32+b4*rs*rsp)
           q2 = log(1+1/q1)
           q3 = a*(b1/rs12+2*b2n+3*b3*rs12+2*b4*p1*rsp)

           q0p = -2*ap*(1+a1p*rs)
           q1p = 2*ap*(b1p*rs12+b2np*rs+b3p*rs32+b4p*rs*rsp)
           q2p = log(1+1/q1p)
           q3p = ap*(b1p/rs12+2*b2np+3*b3p*rs12+2*b4p*p1*rsp)
! vwn 95/12/2 Y.M
           q0q = -2*aq*(1+a1q*rs)
           q1q =  2*aq*(b1q*rs12+b2nq*rs+b3q*rs32+b4q*rs*rsp)
           q2q = log(1+1/q1q)
           q3q = aq*(b1q/rs12+2*b2nq+3*b3q*rs12+2*b4q*p1*rsp)
           if(dabs(zeta) < zeta_minimum) then
              fzd = 8.d0/9.d0 /( 2.d0**thrd4 -2 ) &
                   &     *zeta *( 1 + 5.d0/27 *zeta*zeta &
                   &           *( 1 +22.d0/45 *zeta*zeta ))
           else
              fzd = 4.d0/3.d0 /( 2.d0**thrd4 -2) &
                   &     *( onpzeta**thrd -onmzeta**thrd )
           end if
           onzt4= onmzeta*onpzeta *( 1 + zeta*zeta )
           eu = q0*q2 + ( q0p * q2p - q0 * q2 ) * fz * zeta**4 &
                &             -   q0q * q2q * fz / fzdd * onzt4
           eurs = ( -2*a*a1*q2 -q0*q3/(q1**2 +q1 )) * ( 1 -fz*zeta**4 )&
                &+( -2*ap*a1p*q2p -q0p*q3p/( q1p**2 +q1p )) *fz*zeta**4 &
                &-( -2*aq*a1q*q2q -q0q*q3q/( q1q**2 +q1q )) *fz/fzdd*onzt4
           euzt =  -q0q*q2q /fzdd *( fzd*onzt4 -4*fz*zeta**3 )&
                &        +( q0p*q2p -q0*q2 )*(fzd*zeta +4*fz )*zeta**3 
! vwn 95/12/2 Y.M

           zetadxd = - 2 * chgrhr_l(i, 3-is) * (-1.d0)**is / d
           ecd = eu-thrd*rs*eurs + euzt * zetadxd 

!           if(onpzeta < zeta_minimum2 .or. onmzeta < zeta_minimum2) then
!              ec1  = 0.d0;    ec1d = 0.d0;     ec1dd= 0.d0
!              go to 10
!           end if
           pon = facpon*eu
           b   = delt/(dexp(pon)-1)
           b2  = b*b
           t2  = t*t
           t4  = t2*t2
           t6  = t4*t2
           q4  = 1+b*t2
           q5  = 1+b*t2+b2*t4
           h0  = g3*(bet/delt)*dlog(1+delt*q4*t2/q5)
           h   = h0

           q8  = q5*q5+delt*q4*q5*t2
           h0t = 2*bet*t*(1+2*b*t2)/q8 * g3
           h0b = -bet*t6*(2*b+b2*t2)/q8 * g3
           h0rs = h0b * b * eurs * (b+delt) / bet / g3
           ht   = h0t
           hrs  = h0rs

           if(dabs(zeta) < zeta_minimum) then
              gzt = -2.d0/9*zeta *( 1 +14.d0/27 *zeta*zeta &
                   &             *( 1 +13.d0/18 *zeta*zeta ))
           else
#ifndef PREV_CR_GGAPBE
              gzt = ( (onpzeta**2+eta)**sixthm - (onmzeta**2+eta)**sixthm)/3.d0
#else              
              gzt = ( 1.d0/onpzeta**thrd -1.d0/onmzeta**thrd )/3.d0
#endif
           end if
           bzt = b * ( b+delt) * ( euzt - 3*eu*gzt/g ) / ( bet * g3 ) 
           h0zt = 3 * gzt * h0 / g + h0b * bzt

#ifndef PREV_CR_GGAPBE
           ec1  = d*h * 0.5d0
#else
           ec1  = d*h * 0.5
#endif
! --> revised by Kato-san,  1999Jul16
           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxd
           ec1dd = 0.5*ht/sk / g

10         continue
#ifndef PREV_CR_GGAPBE
           exc0 = eu*d * 0.5d0 + ec1
#else
           exc0 = eu*d * 0.5 + ec1
#endif
           excd = ecd + ec1d
           dF_drho(i, is) = dF_drho(i, is) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(is ==  2) then
              if(dabs(grad_trho(i)) > 1.0d-9) then
                 grad_trho(i) = excdd / grad_trho(i)
              else
                 grad_trho(i) = 0.d0
              endif
           endif
! gradient of charge density 95/12/3 H.S.
           exc1 = exc1 + exc0*f2or1(i)
        end do
     end if
     exc = exc + exc1
  end do
#ifdef __TIMER_DO__
  call timer_end(1001)
#endif
#ifdef __TIMER_SUB__
    call timer_end(767)
#endif
end subroutine cr_ggapbe_3D

!===============================================================================

!!$c===================================================================
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
!!$     1994/11/24 coded by T.Yamasaki
!!$
subroutine xclda_3D(nspin,ispin,chgrhr_l,f2or1,exc,dF_drho,nfft_y)
  use m_Const_Parameters,  only : DP, PAI
  implicit none

  integer, intent(in)       :: nspin,ispin,nfft_y
  real(kind=DP),intent(in)  :: chgrhr_l(1:nfft_y,ispin)
  real(kind=DP),intent(in)  :: f2or1(1:nfft_y)  
  real(kind=DP),intent(out) :: exc
  real(kind=DP),intent(inout) :: dF_drho(1:nfft_y,nspin)

  integer :: is, i
  real(kind=DP) :: co2,d,rs,exclds
  real(kind=DP), parameter :: thrd = 0.333333333333333333d0

#ifdef __TIMER_SUB__
    call timer_sta(768)
#endif
  co2 = (0.75d0/PAI)**thrd
  exc = 0.d0

  is = 1
  if(ispin == 2)  call phase_error_with_msg(6,'!Sorry, ispin == 2 is not supported now',__LINE__,__FILE__)

#ifdef __TIMER_DO__
  call timer_sta(1002)
#endif
  do i = 1, nfft_y
     d = chgrhr_l(i,is)
     if(d < 1.d-40) then
        dF_drho(i,is) = 0.d0
        cycle
     endif
     rs = co2*(1.d0/d)**thrd
     if(rs >= 1.d0) then
        dF_drho(i,is) = &
             &   (-0.61093d0/rs - (0.1423d0+0.1759d0*dsqrt(rs) &
             &    + 0.06366d0*rs)&
             &   /(1 + 1.0529d0*dsqrt(rs) + 0.3334d0*rs)**2 )
        exclds = -0.61093d0*0.75d0/rs &
             &   -0.1423d0/(1 + 1.0529d0*dsqrt(rs) + 0.3334d0*rs)
     else
        dF_drho(i,is) =  &
             &  ( -0.61093d0/rs - 0.05837d0 - 0.0084d0*rs &
             &   + 0.0311d0*log(rs) + 0.00133d0*rs*dlog(rs) )
        exclds = -0.61093d0*0.75d0/rs &
             &  - 0.048d0 -0.0116*rs &
             &  + 0.0311d0*dlog(rs) + 0.0020d0*rs*dlog(rs)
     end if
     exc = exc + exclds*f2or1(i)*d
  end do
#ifdef __TIMER_DO__
  call timer_end(1002)
#endif
#ifdef __TIMER_SUB__
    call timer_end(768)
#endif
end subroutine xclda_3D

!===============================================================================

!!$c ---------------------------------------
!!$c       1994/11/24 coded by T. Yamasaki
!!$c       1994/12/22 modified by T. Yamasaki
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
subroutine ggabek_3D(nspin,ispin,chgrhr_l,grad_rho,f2or1,exc,dF_drho,dF_dgradrho,nfft_y)
  use m_Const_Parameters,  only : DP
  implicit none

  integer,intent(in)        :: nspin,ispin,nfft_y
  real(kind=DP),intent(in)  :: chgrhr_l(1:nfft_y,ispin)
  real(kind=DP),intent(in)  :: grad_rho(1:nfft_y,nspin)
  real(kind=DP),intent(in)  :: f2or1(1:nfft_y)  
  real(kind=DP),intent(out) :: exc
  real(kind=DP),intent(inout) :: dF_drho(1:nfft_y,nspin)
  real(kind=DP),intent(out) :: dF_dgradrho(1:nfft_y,nspin)

  integer        :: i
  real(kind=DP)  :: fac,tthrd,fac_t,d,dd,dthrd,dthrd4,x,x2,p0,x2ovp0&
       &, y,f,vxc1,exgga

  real(kind=DP), parameter ::    b = 0.0042d0
  real(kind=DP), parameter :: thrd = 0.333333333333333333d0
        
#ifdef __TIMER_SUB__
    call timer_sta(769)
#endif
  fac = 4*thrd * b
  tthrd = 2.d0**thrd
  fac_t  = 1.d0/tthrd
#ifdef __TIMER_DO__
  call timer_sta(1003)
#endif
  do i = 1, nfft_y
     d = chgrhr_l(i,1)
     dd = grad_rho(i,1)
     if(d < 1.d-10 .or. dd > 1.d10) then
        dF_dgradrho(i,1) = 0.d0
        cycle
     endif
     dthrd  = d**thrd
     dthrd4 = dthrd*d
     x = tthrd*dd/dthrd4
     x2 = x*x
     p0 = dsqrt(1+x2)
     x2ovp0 = 6*b*x2/p0
     if(x < 1.d-40) then
        y = 0.d0
     else
        y = dlog(x + p0)
     endif
     f  = 1.d0/(1.d0+6.d0*b*x*y)
     vxc1    = x2*f*f*(1.d0-x2ovp0)
     dF_drho(i,1) = dF_drho(i,1) + fac*dthrd*vxc1*fac_t
     dF_dgradrho(i,1) = fac*dthrd*(vxc1 + x2*f)*fac_t

     exgga  = -b * dthrd4 * x2 * f * fac_t
     exc      = exc      + exgga*f2or1(i)
  end do
#ifdef __TIMER_DO__
  call timer_end(1003)
#endif
#ifdef __TIMER_SUB__
    call timer_end(769)
#endif
end subroutine ggabek_3D

!===============================================================================

!!$c ---------------------------------------
!!$c     1994/11/24 coded by T.Yamasaki
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
subroutine ggaprd_3D(nspin,ispin,chgrhr_l,grad_rho,f2or1,exc,dF_drho,dF_dgradrho,nfft_y)
  use m_Const_Parameters,  only : PAI,DP
  implicit none

  integer,intent(in)         :: nspin,ispin,nfft_y
  real(kind=DP),intent(in)   :: chgrhr_l(1:nfft_y,ispin)
  real(kind=DP),intent(in)   :: grad_rho(1:nfft_y,nspin)
  real(kind=DP),intent(in)   :: f2or1(1:nfft_y)  
  real(kind=DP),intent(out)  :: exc
  real(kind=DP),intent(inout):: dF_drho(1:nfft_y,nspin)
  real(kind=DP),intent(out)  :: dF_dgradrho(1:nfft_y,nspin)

  real(kind=DP) :: cinfin,d,ovd,thrd4,sixth7,co2,dd,rs,dthrd,dthrd4&
       &, x, cdenom,c,phi,expphi,ecgga,dcdrs,dcdrho
  integer       :: i

  real(kind=DP), parameter :: alpha = 0.023266d0
  real(kind=DP), parameter ::  beta = 7.389d-6
  real(kind=DP), parameter :: gamma = 8.723d0
  real(kind=DP), parameter :: delta = 0.472d0
  real(kind=DP), parameter ::fchird = 0.11
  real(kind=DP), parameter :: thrd  = 0.333333333333333333d0
  real(kind=DP), parameter :: c0    = 0.001667d0
  real(kind=DP), parameter :: c1    = 0.002568d0

#ifdef __TIMER_SUB__
    call timer_sta(770)
#endif
  cinfin = c0 + c1
  d = 1.d0
  ovd = 1.d0/d
  thrd4 = 4*thrd
  sixth7 = 7.d0/6.d0

  co2 = (0.75d0/PAI)**thrd
#ifdef __TIMER_DO__
  call timer_sta(1004)
#endif
  do i = 1, nfft_y          ! MPI
     d = chgrhr_l(i,1)
     dd = grad_rho(i,1)
     if(d < 1.d-40) cycle
     rs = co2*(1.d0/d)**thrd
     dthrd  = d**thrd
     dthrd4 = dthrd*d
     x = dd/dthrd4
     cdenom = 1.d0/(1 + rs*(gamma + rs*(delta + 1.d+4*beta*rs)))
     c = c0 + (c1 + rs*(alpha + beta*rs))*cdenom
     phi = 1.745d0 * fchird * cinfin/c * x *dsqrt(dthrd)
     expphi = dexp(-phi)
     ecgga = ovd*expphi*c*dd*x

     dcdrs = cdenom * ((alpha + 2*beta*rs) &
          &         - cdenom * (c1+rs*(alpha+beta*rs)) &
          &       * (gamma+rs*(2*delta + 3.d+4*beta*rs)))
     dcdrho = -thrd*(rs/d)*dcdrs
     dF_drho(i,1) = dF_drho(i,1) + ovd*(dcdrho*(phi+1) &
          &       +(sixth7*phi - thrd4)*c)*expphi*x*dd
     dF_dgradrho(i,1) = dF_dgradrho(i,1) + ovd*expphi*c*x*(2 - phi)

     exc = exc + ecgga*f2or1(i)
  end do
#ifdef __TIMER_DO__
  call timer_end(1004)
#endif
#ifdef __TIMER_SUB__
    call timer_end(770)
#endif
end subroutine ggaprd_3D

!===============================================================================

subroutine xcpotf_wigner_3D(nspin,ispin,input_charge,DELTA,chgrhr_l,f2or1,exc,nfft_y)
!*****<<< WIGNER INTERPOLATION >>>
  use m_Const_Parameters,  only : DP,UP,DOWN,PAI,Valence_plus_PC_Charge
  implicit none

  integer,intent(in)         :: nspin,ispin,input_charge,nfft_y
  real(kind=DP),intent(in)   :: DELTA
  real(kind=DP),intent(inout):: chgrhr_l(1:nfft_y,ispin)
  real(kind=DP),intent(in)   :: f2or1(1:nfft_y)  
  real(kind=DP),intent(inout):: exc

  real(kind=DP) :: rho,rs,zeta,vxc,tmp,WX0,CO2
  integer       :: i, icount
  real(kind=DP), parameter :: C13 = 1.d0/3.d0
  real(kind=DP), parameter :: CO3 = 4.d0/3.d0

#ifdef __TIMER_SUB__
    call timer_sta(778)
#endif
  CO2 = (0.75d0/PAI)**C13
  WX0 =  3*(((9*PAI)/4)**C13)/(2*PAI)

  if(ispin == 2 .and. input_charge == Valence_plus_PC_Charge) then
     icount = 0
#ifdef __TIMER_DO__
  call timer_sta(1013)
#endif
     do i = 1, nfft_y
        rho    = chgrhr_l(i,UP) + chgrhr_l(i,DOWN)
        rs     = CO2*( 1.0d0/rho  )**C13
        if(rho < DELTA) cycle
        zeta = (chgrhr_l(i,UP) - chgrhr_l(i,DOWN))/rho
        if(dabs(zeta) > DELTA) icount = icount + 1
        tmp    = - WX0/2.d0*CO3/RS -(0.440d0*CO3*RS + 3.432d0)/(RS+7.8D0)**2
        chgrhr_l(i,UP  ) = tmp
        chgrhr_l(i,DOWN) = tmp
        vxc    = -(WX0/2.d0/RS+0.44d0/(RS+7.8d0))
        exc    = exc + vxc*rho*f2or1(i)
     end do
#ifdef __TIMER_DO__
  call timer_end(1013)
#endif
     if(icount >= 1) then
        write(6,*)' WIGNER XC-ENERGY and POTENTIAL are not' &
             &         , ' appropriate for Spin-Polarized calculation'
        call phase_error_with_msg(6,'WIGNER',__LINE__,__FILE__)
     endif
  else
#ifdef __TIMER_DO__
  call timer_sta(1014)
#endif
     do i = 1, nfft_y
        rho    = chgrhr_l(i,1)
        rs     = CO2*( 1.0d0/rho  )**C13
        chgrhr_l(i,1)   = - WX0/2.d0*CO3/rs &
             &           -(0.440d0*CO3*rs + 3.432d0)/(rs+7.8d0)**2
        vxc    = -(WX0/2.d0/rs+0.44d0/(rs+7.8d0))
        exc    = exc + vxc*rho*f2or1(i)
     end do
#ifdef __TIMER_DO__
  call timer_end(1014)
#endif
  endif
#ifdef __TIMER_SUB__
    call timer_end(778)
#endif
end subroutine xcpotf_wigner_3D

!===============================================================================

subroutine xcpotf_pzold_3D(nspin,ispin,input_charge,DELTA,chgrhr_l,f2or1,exc,nfft_y)
  use m_Const_Parameters,  only : DP,UP,DOWN,PAI,Valence_plus_PC_Charge
  implicit none

  integer,intent(in)         :: nspin,ispin,input_charge,nfft_y
  real(kind=DP),intent(in)   :: DELTA
  real(kind=DP),intent(inout):: chgrhr_l(1:nfft_y,ispin)
  real(kind=DP),intent(in)   :: f2or1(1:nfft_y)  
  real(kind=DP),intent(inout):: exc

  real(kind=DP) :: B1P01,B1P02,CP00,CP01,CP02,OVC213,OVC243 &
       &         , B1F01,B1F02,CF00,CF01,CF02
  real(kind=DP) :: rho,ovn,rs,denmp,denmf,wcp,wcf,vcp,vcf,zeta,f00,f2up,f2down &
       &         , t1,wxp,wxf,vxp,x,t2,vxc
  integer       :: i

  real(kind=DP), parameter   :: C13   =  1.d0/3.d0
  real(kind=DP), parameter   :: C43   =  4.d0/3.d0
  real(kind=DP) :: CO2, WX0, C213

  real(kind=DP), parameter :: GAMP   = -0.1423d0, GAMF   = -0.0843d0 &
       &                    , BETA1P =  1.0529d0, BETA2P =  0.3334d0 &
       &                    , BETA1F =  1.3981d0, BETA2F =  0.2611d0 &
       &                    , AP     =  0.0311d0, BP     = -0.0480d0 &
       &                    , AF     =  0.01555d0,BF     = -0.0269d0 &
       &                    , CP     =  0.0020d0, DPx    = -0.0116d0 &
       &                    , CF     =  0.00070d0,DF     = -0.0048d0

#ifdef __TIMER_SUB__
    call timer_sta(779)
#endif
  CO2   = (0.75d0/PAI)**C13
  WX0   =  3*(((9*PAI)/4)**C13)/(2*PAI)
  C213  = 2.d0**C13

  B1P01 = 7.d0*BETA1P/6.d0
  B1P02 = 4.d0*BETA2P/3.d0

  CP00 = BP - AP/3.d0
  CP01 = (2*DPx - CP)/3.d0
  CP02 = 2*CP/3.d0

#ifdef __TIMER_DO__
  call timer_sta(1015)
#endif
  if(ispin == 2 .and. input_charge == Valence_plus_PC_Charge) then
     OVC213 = 1.D0/(C213 - 1)
     OVC243 = 1.D0/(2.D0**C43 - 2)

     B1F01 = 7.D0*BETA1F/6.D0
     B1F02 = 4.D0*BETA2F/3.D0
     CF00 = BF - AF/3.D0
     CF01 = (2*DF - CF)/3.D0
     CF02 = 2*CF/3.D0

     do i = 1, nfft_y   ! MPI
        rho    = chgrhr_l(i, UP) + chgrhr_l(i,DOWN)
        if(rho <= DELTA) then
           chgrhr_l(i,UP  ) = 0.D0
           chgrhr_l(i,DOWN) = 0.D0
           cycle
        endif
        ovn = 1.d0/rho
        rs     = CO2*ovn**C13
        if(rs >= 1.d0) then
           denmp = 1 + BETA1P*sqrt(rs) + BETA2P*rs
           denmf = 1 + BETA1F*sqrt(rs) + BETA2F*rs
           wcp = GAMP/denmp
           wcf = GAMF/denmf
           vcp = wcp*(1+B1P01*sqrt(rs)+B1P02*rs)/denmp
           vcf = wcf*(1+B1F01*sqrt(rs)+B1F02*RS)/denmf
        else
           wcp = AP*dlog(rs) + BP + CP*rs*dlog(rs) + DPx*rs
           wcf = AF*dlog(rs) + BF + CF*rs*dlog(rs) + DF*rs
           vcp = (AP+CP02*RS)*dLOG(RS) + CP00 + CP01*RS
           vcf = (AF+CF02*RS)*dLOG(RS) + CF00 + CF01*RS
        end if
        zeta   = (chgrhr_l(i,UP) - chgrhr_l(i,DOWN))*ovn
        f00    = ((1+zeta)**C43 + (1-zeta)**C43 - 2)*OVC243
        f2up   = ((1+zeta)**C13 - 1)*OVC213
        f2down = ((1-zeta)**C13 - 1)*OVC213
        t1 = vcf - vcp - C43*(wcf - wcp)
        wxp = -WX0*0.5d0/rs
        wxf = wxp * C213
        vxp = wxp * C43
        x = vxp + vcp + t1*f00
        t2 = C43 * ((WXF - WXP) + ( WCF - WCP ))
        chgrhr_l(i,UP  ) = x + t2*f2up
        chgrhr_l(i,DOWN) = x + t2*f2down

        vxc = wxp + wcp + ((wxf - wxp) + (wcf - wcp))*f00
        exc = exc + vxc*rho*f2or1(i)
     end do
  else
     do i = 1, nfft_y
        rho    = chgrhr_l(i,1)
        if(rho < DELTA) then
           chgrhr_l(i,1) = 0.d0
           cycle
        endif
        ovn    = 1.d0/rho
        rs     = CO2*ovn**C13
        wxp    = -WX0*0.5d0/rs
        if(rs >= 1.d0) then
           denmp = 1 + BETA1P*sqrt(rs) + BETA2P*rs
           wcp = GAMP/denmp
           vcp = wcp*(1.D0+B1P01*SQRT(rs)+B1P02*rs)/denmp
        else
           wcp = AP*dlog(rs) + BP + CP*rs*log(rs) + DPx*rs
           vcp = (AP+CP02*rs)*log(rs) + CP00 + CP01*rs
        end if

        wxp = -WX0*0.5D0/rs
        wxf = WXP * C213
        vxp = WXP * C43
        x = vxp + vcp
        chgrhr_l(i,1  ) = x
        vxc = wxp + wcp
        exc = exc + vxc*rho*f2or1(i)
     end do
  endif
#ifdef __TIMER_DO__
  call timer_end(1015)
#endif
#ifdef __TIMER_SUB__
    call timer_end(779)
#endif
end subroutine xcpotf_pzold_3D

!===============================================================================

subroutine xcpotf_xalfa_3D(nspin,ispin,input_charge,DELTA,chgrhr_l,f2or1,exc,nfft_y)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,UP,DOWN,PAI,Valence_plus_PC_Charge
  implicit none

  integer,intent(in)         :: nspin,ispin,input_charge,nfft_y
  real(kind=DP),intent(in)   :: DELTA
  real(kind=DP),intent(inout):: chgrhr_l(1:nfft_y,ispin)
  real(kind=DP),intent(in)   :: f2or1(1:nfft_y)  
  real(kind=DP),intent(inout):: exc

  real(kind=DP), parameter   :: C13   =  1.d0/3.d0
  real(kind=DP), parameter   :: ALPHA = 0.7d0
  real(kind=DP)              :: CO1, CO2
  real(kind=DP)              :: rho,tmp,vxc
  integer                    :: i
#ifdef __TIMER_SUB__
    call timer_sta(780)
#endif

  CO1   = 1.5d0*ALPHA*(3.D0/PAI)**C13
  CO2   = (0.75d0/PAI)**C13

#ifdef __TIMER_DO__
  call timer_sta(1016)
#endif
  if(ispin == 2 .and. input_charge == Valence_plus_PC_Charge) then
     do i = 1, nfft_y
        rho    = chgrhr_l(i,UP) + chgrhr_l(i,DOWN)
!!$        rs     = CO2*( 1.0D0/rho  )**C13
        tmp    = - CO1*(rho**C13)
        chgrhr_l(i,UP)   = tmp
        chgrhr_l(i,DOWN) = tmp
        vxc = 0.25d0*CO1*(rho**C13)
        exc = exc+vxc*rho*f2or1(i)
     end do
  else
     do i = 1, nfft_y
        rho    = chgrhr_l(i,1)
!!$        rs     = CO2*( 1.0d0/rho  )**C13
        chgrhr_l(i,1)   = - CO1*(rho**C13)
        vxc = 0.25d0*CO1*(rho**C13)
        exc = exc+vxc*rho*f2or1(i)
     end do
  endif
#ifdef __TIMER_DO__
  call timer_end(1016)
#endif
#ifdef __TIMER_SUB__
    call timer_end(780)
#endif
end subroutine xcpotf_xalfa_3D

!===============================================================================

subroutine xcpotf_pz_3D(nspin,ispin,input_charge,DELTA,chgrhr_l,f2or1,exc,nfft_y)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,UP,DOWN,PAI,Valence_plus_PC_Charge
  implicit none

  integer,intent(in)         :: nspin,ispin,input_charge,nfft_y
  real(kind=DP),intent(in)   :: DELTA
  real(kind=DP),intent(inout):: chgrhr_l(1:nfft_y,ispin)
  real(kind=DP),intent(in)   :: f2or1(1:nfft_y)  
  real(kind=DP),intent(inout):: exc

  real(kind=DP), parameter :: C13   =  1.d0/3.d0,  C23   =  2.d0/3.d0,  C43   =  4.d0/3.d0
  real(kind=DP), parameter :: GAMP   = -0.1423d0, GAMF   = -0.0843d0 &
       &                    , BETA1P =  1.0529d0, BETA2P =  0.3334d0 &
       &                    , BETA1F =  1.3981d0, BETA2F =  0.2611d0 &
       &                    , AP     =  0.0311d0, BP     = -0.0480d0 &
       &                    , CP     =  0.0020d0, DPx    = -0.0116d0 &
       &                    , AF     =  0.01555d0,BF     = -0.0269d0 &
       &                    , CF     =  0.00070d0,DF     = -0.0048d0

  real(kind=DP)  :: B1P01,B1P02,CP00,CP01,CP02,C243 &
       &          , B1F01,B1F02,CF00,CF01,CF02
  real(kind=DP)  :: CO2, WX0, C213
  real(kind=DP)  :: rho,ovn,rs,zeta,wxp,f00,wx,vx1,vx2,denmp,denmf,wcp,wcf &
       &          , vcp,vcf,wc,t1,f11,t2,vc1,vc2,vxc,vxp
  integer        :: i
#ifdef __TIMER_SUB__
    call timer_sta(781)
#endif

  CO2   = (0.75d0/PAI)**C13
  WX0   =  3*(((9*PAI)/4)**C13)/(2*PAI)
  B1P01 = 7*BETA1P/6.d0
  B1P02 = 4*BETA2P/3.d0
  CP00 = BP - AP/3.D0
  CP01 = (2*DPx - CP)/3.D0
  CP02 = 2*CP/3.D0

#ifdef __TIMER_DO__
  call timer_sta(1017)
#endif
  if(ispin == 2 .and. input_charge == Valence_plus_PC_Charge) then
     C213 = 2.d0**C13 - 1
!!$     OVC213 = 1.d0/C213
     C243   = 2*C213
!           ovc243 = 1.d0/(2.d0**C43 - 2)

     B1F01 = 7*BETA1F/6.d0
     B1F02 = 4*BETA2F/3.d0
     CF00 = BF - AF/3.d0
     CF01 = (2*DF - CF)/3.d0
     CF02 = 2*CF/3.d0

     do i = 1, nfft_y
        rho    = chgrhr_l(i, UP) + chgrhr_l(i,DOWN)
        if(rho <= DELTA) then
           chgrhr_l(i,UP  ) = 0.D0
           chgrhr_l(i,DOWN) = 0.D0
           cycle
        end if
        ovn = 1.d0/rho
        rs     = CO2*ovn**C13
        zeta   = (chgrhr_l(i,UP) - chgrhr_l(i,DOWN))*ovn
        wxp    = -WX0*0.5d0/rs
        f00 = ( (1+zeta)**C43 + (1-zeta)**C43 - 2 ) / c243
        wx  = wxp*(1+C213*f00)
        vx1 = c43*wxp*(1+zeta)**c13
        vx2 = c43*wxp*(1-zeta)**c13
        if(rs >= 1.d0) then
           denmp = 1.d0 + BETA1P*sqrt(rs) + BETA2P*rs
           denmf = 1.D0 + BETA1F*sqrt(rs) + BETA2F*rs
           wcp = GAMP/denmp
           wcf = GAMF/denmf
           vcp = wcp*(1.d0+B1P01*sqrt(rs)+B1P02*rs)/denmp
           vcf = wcf*(1.D0+B1F01*sqrt(rs)+B1F02*rs)/denmf
        else
           wcp = AP*dlog(rs) + BP + CP*rs*dlog(rs) + DPx*rs
           wcf = AF*dlog(rs) + BF + CF*rs*dlog(rs) + DF*rs
           vcp = (AP+CP02*rs)*dLOG(rs) + CP00 + CP01*rs
           vcf = (AF+CF02*rs)*dLOG(rs) + CF00 + CF01*rs
        end if
        wc     = wcp + f00*(wcf - wcp)
        t1     = vcp + f00*(vcf - vcp)
        f11    = C23/C213*((1+zeta)**c13-(1-zeta)**c13)
        t2     = f11*(wcf-wcp)
        vc1    = t1+t2*(1-zeta)
        vc2    = t1-t2*(1+zeta)

        chgrhr_l(i,UP  ) = vx1 + vc1
        chgrhr_l(i,DOWN) = vx2 + vc2

        vxc    = wx + wc
        exc = exc + vxc*rho*f2or1(i)
     end do
  else
     do i = 1, nfft_y
        rho    = chgrhr_l(i, 1)
        if(rho <= DELTA) then
           chgrhr_l(i,UP  ) = 0.D0
           cycle
        end if
        ovn = 1.d0/rho
        rs     = CO2*ovn**C13
        wxp    = -WX0*0.5d0/rs
        vxp = C43*wxp
        if(rs >= 1.d0) then
           denmp = 1 + BETA1P*sqrt(rs) + BETA2P*rs
           wcp = GAMP/denmp
           vcp = wcp*(1+B1P01*sqrt(rs)+B1P02*rs)/denmp
        else
           wcp = AP*dlog(rs) + BP + CP*rs*dlog(rs) + DPx*rs
           vcp = (AP+CP02*rs)*dLOG(rs) + CP00 + CP01*rs
        end if
        wc     = wcp
        chgrhr_l(i,1   ) = vxp + vcp

        vxc    = wxp + wcp
        exc = exc + vxc*rho*f2or1(i)
     end do
  endif
#ifdef __TIMER_DO__
  call timer_end(1017)
#endif
#ifdef __TIMER_SUB__
    call timer_end(781)
#endif
end subroutine xcpotf_pz_3D

!===============================================================================

subroutine xcpotf_vwn_3D(nspin,ispin,input_charge,DELTA,chgrhr_l,f2or1,exc,nfft_y)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,UP,DOWN,PAI,Valence_plus_PC_Charge
  implicit none

  integer,intent(in)         :: nspin,ispin,input_charge,nfft_y
  real(kind=DP),intent(in)   :: DELTA
  real(kind=DP),intent(inout):: chgrhr_l(1:nfft_y,ispin)
  real(kind=DP),intent(in)   :: f2or1(1:nfft_y)  
  real(kind=DP),intent(inout):: exc

  real(kind=DP), parameter :: B1 = 1.13107d0,    C1 = 13.0045d0 &
       &                    , D1 = -0.0047584d0 &
       &                    , A2 =  0.0621814d0, B2 =  3.72744d0 &
       &                    , C2 = 12.9352d0   , D2 = -0.10498d0 &
       &                    , A3 =  0.0310907d0, B3 =  7.06042d0 &
       &                    , C3 = 18.0578d0,    D3 = -0.32500d0
  real(kind=DP), parameter :: C13   =  1.d0/3.d0,  C43   =  4.d0/3.d0
  real(kind=DP) :: CO2, WX0, C213

  real(kind=DP) :: F20,A1,Q1,P01,Q2,P02,Q3,P03,C243
  real(kind=DP) :: rho,ovn,zeta,f00,wx,vx1,vx2,x,p,r,h00,h10,h0p,h1p,h0f,h1f &
       &         , z4,h000,wc,vc,g00,vc0,vc1,vc2,rs,wxp,vxc
  integer       :: i

#ifdef __TIMER_SUB__
    call timer_sta(782)
#endif
  CO2   = (0.75d0/PAI)**C13
  WX0   =  3*(((9*PAI)/4)**C13)/(2*PAI)
  C213 = 2**C13 - 1

  F20   = 4/(9*C213)
  A1   = -1/(3*PAI*PAI)
  Q1   = dsqrt(4*C1-B1**2)
  P01  = D1**2 + B1*D1 + C1

  Q2   = dsqrt(4*C2-B2**2)
  P02  = D2**2 + B2*D2 + C2

  Q3   = dsqrt(4*C3-B3**2)
  P03  = D3**2 + B3*D3 + C3

#ifdef __TIMER_DO__
  call timer_sta(1018)
#endif
  if(ispin ==  2 .and. input_charge == Valence_plus_PC_Charge) then
!!$     OVC213 = 1.d0/C213
     C243   = 2*C213
!           ovc243 = 1.d0/(2.d0**c43 - 2)

     do i = 1, nfft_y
        rho    = chgrhr_l(i, UP) + chgrhr_l(i,DOWN)
        if(rho <= DELTA) then
           chgrhr_l(i,UP  ) = 0.D0
           chgrhr_l(i,DOWN) = 0.D0
           cycle
        end if
        ovn = 1.d0/rho
        rs     = CO2*ovn**C13
        if(rs <= DELTA) then
           chgrhr_l(i, UP)   = 0.d0  
           chgrhr_l(i, DOWN) = 0.d0
           cycle
        end if
        zeta   = (chgrhr_l(i,UP) - chgrhr_l(i,DOWN))*ovn
        wxp    = -WX0*0.5d0/rs
        f00 = ( (1+zeta)**C43 + (1-zeta)**C43 - 2 ) / C243
        wx  = wxp*(1+C213*f00)
        vx1 = C43*wxp*(1+zeta)**C13
        vx2 = C43*wxp*(1-zeta)**C13
        x   = dsqrt(rs)
        p   = x**2 + B1*x + C1
        r   = datan(Q1/(2*x+B1))
        h00= (A1/F20) * ( dlog((x**2)/p)+2*B1*r/Q1 &
             &       -B1*D1*(dlog(((x-D1)**2)/p)+2*(B1+2*D1)*r/Q1)/P01 )
        h10= (A1/F20) * ( 2/x - (1-B1*D1/P01)*(2*x+B1)/p - 2*B1*D1/((x-D1)*P01) &
             &                 - 4*B1*(1-(B1+2*D1)*D1/P01)/(Q1**2+(2*x+B1)**2) )

        p   = x**2 + B2*x + C2
        r   = datan(Q2/(2*x+B2))
        h0p= A2 * ( dlog((x**2)/p)+2*B2*r/Q2 &
             &        -B2*D2*(dlog(((x-D2)**2)/p)+2*(B2+2*D2)*r/Q2)/P02 )
        h1p= A2 * ( 2/x - (1-B2*D2/P02)*(2*x+B2)/p - 2*B2*D2/((x-D2)*P02) &
             &        -4*B2*(1-(B2+2*D2)*D2/P02)/(Q2**2+(2*x+B2)**2) )

        p   = x**2 + B3*x + C3
        r   = datan(Q3/(2*x+B3))
        h0f = A3 * ( dlog((x**2)/p)+2*B3*r/Q3 &
             &       -B3*D3*(dlog(((x-D3)**2)/p)+2*(B3+2*D3)*r/Q3)/P03 )
        h1f = A3 * ( 2/x - (1-B3*D3/P03)*(2*x+B3)/p - 2*B3*D3/((x-D3)*P03) &
             &         - 4*B3*(1-(B3+2*D3)*D3/P03)/(Q3**2+(2*x+B3)**2) )

        z4   = zeta**4
        h000 = h0f - h0p - h00
        wc=h0p+h00*f00+h000*z4*f00
        vc=wc-x*((1-z4*f00)*h1p+z4*f00*h1f+(1-z4)*f00*h10)/6
        g00 = 4*((1+zeta)**c13-(1-zeta)**c13) / (6*c213)
        vc0 = 4* h000*(zeta**3)*f00 + (h000*z4+h00)*g00
        vc1=vc+(1-zeta)*vc0
        vc2=vc-(1+zeta)*vc0

        wc  = wc*0.5d0
        vc1 = vc1*0.5d0
        vc2 = vc2*0.5d0

        chgrhr_l(i,UP  ) = vx1 + vc1
        chgrhr_l(i,DOWN) = vx2 + vc2

        vxc    = wx + wc
        exc = exc + vxc*rho*f2or1(i)
     end do
  else
!!$     OVC213 = 1.d0/C213
     C243   = 2*C213

     do i = 1, nfft_y
        rho    = chgrhr_l(i, 1)
        if(rho <= DELTA) then
           chgrhr_l(i,1  ) = 0.D0
           cycle
        end if
        ovn = 1.d0/rho
        rs     = CO2*ovn**C13
        wxp    = -WX0*0.5d0/rs
        f00 = 0.d0
        wx  = wxp
        vx1 = c43*wxp

        x   = dsqrt(rs)
        p   = x**2 + B2*x + C2
        r   = datan(Q2/(2*x+B2))
        h0p= a2 * ( dlog((x**2)/p)+2*b2*r/q2 &
             &       -b2*d2*(dlog(((x-d2)**2)/p)+2*(b2+2*d2)*r/q2)/p02 )
        h1p= a2 * ( 2/x - (1-b2*d2/p02)*(2*x+b2)/p - 2*b2*d2/((x-d2)*p02) &
             &           - 4*b2*(1-(b2+2*d2)*d2/p02)/(q2**2+(2*x+b2)**2) )

        wc  = h0p
        vc1 = wc-x*h1p/6

        wc  = wc*0.5d0
        vc1 = vc1*0.5d0

        chgrhr_l(i,1  ) = vx1 + vc1
        vxc    = wx + wc
        exc = exc + vxc*rho*f2or1(i)
     end do
  endif
#ifdef __TIMER_DO__
  call timer_end(1018)
#endif
#ifdef __TIMER_SUB__
    call timer_end(782)
#endif
end subroutine xcpotf_vwn_3D

!===============================================================================

subroutine xcpotf_mjw_bh_gl_3D(len_xctype,xctype,nspin,ispin,input_charge,DELTA,chgrhr_l,f2or1,exc,nfft_y)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,UP,DOWN,PAI,Valence_plus_PC_Charge
  implicit none

  integer,intent(in)         :: len_xctype,nspin,ispin,input_charge,nfft_y
  character(len=len_xctype)  :: xctype
  real(kind=DP),intent(in)   :: DELTA
  real(kind=DP),intent(inout):: chgrhr_l(1:nfft_y,ispin)
  real(kind=DP),intent(in)   :: f2or1(1:nfft_y)  
  real(kind=DP),intent(inout):: exc

  real(kind=DP), parameter ::  CP    =  0.0020d0,  CF     =  0.00070d0
  real(kind=DP), parameter ::  C13   =  1.d0/3.d0, C43   =  4.d0/3.d0

  real(kind=DP) :: CO2, WX0, C213
  real(kind=DP) :: C243, CPX, RP, CFX, RF
  real(kind=DP) :: rho,ovn,rs,zeta,wxp,f00,wx,vx1,vx2,rsp,rsf,aap,aaf &
       &         , wcp,wcf,vcp,vcf,t1,t2,f11,f12,wc,vc1,vc2,vxc
  integer       :: i

#ifdef __TIMER_SUB__
    call timer_sta(783)
#endif
  CO2   = (0.75d0/PAI)**C13
  WX0   =  3*(((9*PAI)/4)**C13)/(2*PAI)
  C213 = 2**C13 - 1
  if(xctype ==  'mjw    ') then
! --*   moruzzi-janak-williams; Phys.Rev.12(1975)1257.  
     CPX=0.0450d0
     RP=21.0d0
  else if(xctype == 'bh     ') then                                    
! --*   von barth-hedin; J.Phys.C:Solid State Phys.5(1972)1629.  
     CPX=0.0504d0
     RP=30.0d0
  else if(xctype == 'gl     ') then                                    
! --*   gunnarsson-lundqvist; Phys.Rev.B13(1976)4274. 
     CPX=0.0666d0
     RP=11.4d0
  end if
#ifdef __TIMER_DO__
  call timer_sta(1019)
#endif

  if(ispin ==  2 .and. input_charge == Valence_plus_PC_Charge) then
     if(xctype ==  'mjw    ') then
        CFX=cp/2
        RF=RP*(2.d0**C43)
     else if(xctype ==  'bh     ') then                                    
        CFX=0.0254d0
        RF=75.0d0
     else if(xctype ==  'gl     ') then                                    
        CFX=0.0406d0
        RF=15.9d0
     end if

     C243   = 2*C213
     do i = 1, nfft_y
        rho    = chgrhr_l(i, UP) + chgrhr_l(i,DOWN)
        if(rho <= DELTA) then
           chgrhr_l(i,UP  ) = 0.D0
           chgrhr_l(i,DOWN) = 0.D0
           cycle
        end if
        ovn = 1.d0/rho
        rs     = CO2*ovn**C13
        if(rs <= DELTA) then
           chgrhr_l(i, UP  ) = 0.d0
           chgrhr_l(i, DOWN) = 0.d0
           cycle
        end if
        zeta   = (chgrhr_l(i,UP) - chgrhr_l(i,DOWN))*ovn
        wxp    = -WX0*0.5d0/rs
        f00 = ( (1+zeta)**C43 + (1-zeta)**C43 - 2 ) / C243
        wx  = wxp*(1+C213*f00)
        vx1 = C43*wxp*(1+zeta)**C13
        vx2 = C43*wxp*(1-zeta)**C13
        rsp = rs/RP
        rsf = rs/RF
        aap = dlog(1+1/rsp)                                      
        aaf = dlog(1+1/rsf)                                      
        wcp = - CPX * ((1+rsp**3)*aap - rsp**2 + rsp/2 - C13)
        wcf = - CFX * ((1+rsf**3)*aaf - rsf**2 + rsf/2 - C13)
        vcp = - CP * aap
        vcf = - CF * aaf

        t1  = vcf - vcp - C43*(wcf-wcp) 
        t2  = C43 * (wcf-wcp)

        f11= ( (1+zeta)**C13 -1 ) / C213
        f12= ( (1-zeta)**C13 -1 ) / C213
        wc=wcp + (wcf - wcp)*f00
        vc1=vcp + t1*f00 + t2*f11
        vc2=vcp + t1*f00 + t2*f12

        wc = wc*0.5d0
        vc1 = vc1*0.5d0
        vc2 = vc2*0.5d0

        chgrhr_l(i,UP  ) = vx1 + vc1
        chgrhr_l(i,DOWN) = vx2 + vc2

        vxc    = wx + wc
        exc = exc + vxc*rho*f2or1(i)
     end do
  else
     c243   = 2*C213
     do i = 1, nfft_y
        rho    = chgrhr_l(i, 1)
        if(rho <= DELTA) then
           chgrhr_l(i,1   ) = 0.D0
           cycle
        end if
        ovn = 1.d0/rho
        rs     = CO2*ovn**C13
        if(rs <= DELTA) then
           chgrhr_l(i, 1   ) = 0.d0
           cycle
        end if
        zeta   = 0.d0
        wxp    = -WX0*0.5d0/rs
        wx  = wxp
        vx1 = c43*wxp
        vx2 = c43*wxp

        rsp = rs/RP
        aap = dlog(1+1/rsp)                                      
        wcp = - CPX * ((1+rsp**3)*aap - rsp**2 + rsp/2 - C13)
        vcp = - CP * aap

        wc=wcp * 0.5d0
        vc1=vcp * 0.5d0

        chgrhr_l(i,1   ) = vx1 + vc1

        vxc    = wx + wc
        exc = exc + vxc*rho*f2or1(i)
     end do
  endif
#ifdef __TIMER_DO__
  call timer_end(1019)
#endif
#ifdef __TIMER_SUB__
    call timer_end(783)
#endif
end subroutine xcpotf_mjw_bh_gl_3D

!===============================================================================

!$$#endif

! ========= KT_add === 13.0A
! ==================================== KT_add ========================== 13.0A
#ifdef __EDA__
subroutine ex_gga_library( nspin, ispin, ista_r, iend_r, chgrhr_l, grad_rho, &
     &                     wos, exc, dFx_drho, dFx_dgradrho, exc_on_a_grid_wk, &
     &                     pot_type, ist, ien )
#else
subroutine ex_gga_library( nspin, ispin, ista_r, iend_r, chgrhr_l, grad_rho, &
     &                     wos, exc, dFx_drho, dFx_dgradrho, pot_type, ist, ien )
#endif
  use m_Const_Parameters,  only : DP,PAI
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif
  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  integer,intent(in)        :: ist, ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)  :: grad_rho(ista_r:iend_r,nspin)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)
  real(kind=DP),intent(out) :: exc
  real(kind=DP),intent(out) :: dFx_drho(ista_r:iend_r,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif
  integer, intent(in) :: pot_type

!!  real(kind=DP), parameter :: ax = -0.7385588d0
  real(kind=DP), parameter :: ax = -0.7385587663820224d0 !! ax = -(3/4)*(3/PI)^(1/3)
  real(kind=DP), parameter :: thrd = 0.33333333333d0, thrd4 = 1.333333333333333d0
  real(kind=DP), parameter :: thpith = 3.0936677262801d0

! ----- parameters of enhanced factor --
  real(kind=DP), parameter :: mu_gel    = 0.1234567901234568d0   !  = 10.0 /81.0

  real(kind=DP), parameter :: mu_pbe    = 0.2195149727645171d0
  real(kind=DP), parameter :: mu_revpbe = 0.2195149727645171d0
  real(kind=DP), parameter :: mu_rpbe   = 0.2195149727645171d0
  real(kind=DP), parameter :: mu_wc     = 0.2195149727645171d0
  real(kind=DP), parameter :: mu_pbesol = 0.1234567901234568d0   !  = 10.0 /81.

  real(kind=DP), parameter :: mu_b86r   = 0.1234567901234568d0   !  = 10.0 /81.
  real(kind=DP), parameter :: mu_optpbe = 0.175519d0
  real(kind=DP), parameter :: mu_optb86b = 0.1234567901234568d0   !  = 10.0 /81.0
  real(kind=DP), parameter :: mu_c09    = 0.0617d0
  real(kind=DP), parameter :: mu_lvpw86r = 0.094344444444444444d0  ! = 0.8491 /9.0

  real(kind=DP), parameter :: kappa_pbe    = 0.804d0
  real(kind=DP), parameter :: kappa_revpbe = 1.245d0          ! by Koga-san
  real(kind=DP), parameter :: kappa_rpbe   = 0.804d0
  real(kind=DP), parameter :: kappa_wc     = 0.804d0
  real(kind=DP), parameter :: kappa_pbesol = 0.804d0
  real(kind=DP), parameter :: kappa_pbeint = 0.804d0

  real(kind=DP), parameter :: kappa_b86r   = 0.7114d0
  real(kind=DP), parameter :: kappa_optpbe = 1.04804d0
  real(kind=DP), parameter :: kappa_optb86b = 1.0d0
  real(kind=DP), parameter :: kappa_c09    = 1.245d0

! ---- others
  real(kind=DP), parameter :: coeff_WC = 0.0079325D0      ! ---- Wu and Cohen (2006)
  real(kind=DP), parameter :: alpha_c09 = 0.0483d0

  real(kind=DP), parameter :: s_low_htbs = 0.6d0, s_high_htbs = 2.6d0
!
  real(kind=DP), parameter :: a_pw86r = 0.1234567901234568d0   !  = 10.0 /81.0
  real(kind=DP), parameter :: b_pw86r = 17.33d0
  real(kind=DP), parameter :: c_pw86r = 0.163d0

  real(kind=DP), parameter :: alpha_lvpw86r = 0.02178d0
  real(kind=DP), parameter :: beta_lvpw86r = 1.15d0

  real(kind=DP), parameter :: a1_ev93 = 1.647127D0
  real(kind=DP), parameter :: a2_ev93 = 0.980118D0
  real(kind=DP), parameter :: a3_ev93 = 0.017399D0
  real(kind=DP), parameter :: b1_ev93 = 1.523671D0
  real(kind=DP), parameter :: b2_ev93 = 0.367229D0
  real(kind=DP), parameter :: b3_ev93 = 0.011282D0

  real(kind=DP), parameter :: beta_lb94 = 0.05d0
! ----------------
  integer       :: is, i, pot_mode
  integer       :: istart, iend
  real(kind=DP) :: facw, d, dd, fk, s, fac, s2, mu_s2, f, ex, df_ds, &
       &           exd, exdd, exc0, excd, excdd, exc1, pot_add
  real(kind=DP) :: ctmp1, ctmp2, ctmp3, ctmp4, ctmp5, ctmp6, x, dx_ds, s4, s6, df_dx
  real(kind=DP) :: y, dy_ds
  real(kind=DP) :: alpha_pbeint, mu_this, dmu_ds
!
  logical :: First = .true.
  real(kind=DP), save :: coeff_spline_htbs(0:5)

  istart = ist;  iend = ien

  if ( pot_type == 5 ) then
     if ( First ) then
        call init_spline_htbs( mu_wc, kappa_rpbe, kappa_wc, coeff_wc, &
             &                 s_low_htbs, s_high_htbs, coeff_spline_htbs )
        First = .false.
     endif
  endif

!---- Spin dependency

  facw = ispin
  exc  = 0.d0
#ifdef __EDA__
! -----  ascat starts modifying  -----
  if(sw_eda==ON) then
  exc_on_a_grid_wk = 0.d0
  endif
! -----  ascat ceases modifying  -----
#endif
#ifdef NEC_TUNE_MXCP
!CDIR NOCONCUR
#endif
  do is = 1, ispin
     exc1 = 0.d0
#ifdef NEC_TUNE_MXCP
!CDIR INNER
#endif
     do i = istart, iend
        d  = facw * chgrhr_l(i, is)
        dd = facw * grad_rho(i, is)
        fk = (3.d0*PAI*PAI*d)**thrd
!!$        sk = dsqrt(4*fk/PAI)
        if(d > 1.d-05) then
           s = dd/(d*fk*2.d0)
        else
           s = 0.d0
        endif
!-------------------------------------
        fac = ax*d**thrd
        s2  = s*s
!
        pot_mode = pot_type
        if ( pot_type == 5 ) then
           if ( s <= s_low_htbs ) then
              pot_mode = 4
           else if ( s >= s_high_htbs ) then
              pot_mode = 3
           endif
        endif
!
        select case ( pot_mode )
        case (1)                        ! pbe
           x = mu_pbe *s2;   dx_ds = 2.0d0 *mu_pbe *s
           ctmp1 = kappa_pbe + x

           f = 1.0d0 +kappa_pbe *x /ctmp1
           df_dx = kappa_pbe**2 /ctmp1**2
           df_ds = dx_ds *df_dx

        case (2)                        ! revpbe
           x = mu_revpbe *s2;   dx_ds = 2.0d0 *mu_revpbe *s
           ctmp1 = kappa_revpbe + x

           f = 1.0d0 +kappa_revpbe *x /ctmp1
           df_dx = kappa_revpbe**2 /ctmp1**2
           df_ds = dx_ds *df_dx

        case (3)                        ! rpbe
           x = mu_rpbe *s2;   dx_ds = 2.0d0 *mu_rpbe *s
           ctmp1 = exp( -x/ kappa_rpbe )

           f = 1.0d0 +kappa_rpbe -kappa_rpbe *ctmp1
           df_dx = ctmp1
           df_ds = dx_ds *df_dx

        case (4)                        ! W. Cohen
           s4 = s2 *s2
           ctmp1 = exp( -s2 )
           ctmp2 = coeff_WC *s4

           x = mu_gel *s2 + ( mu_wc -mu_gel )*s2 *ctmp1 + log( 1.0d0 +ctmp2 )
           dx_ds = 2.0d0 *mu_gel *s + ( mu_wc -mu_gel ) *2.0d0 *s *ctmp1 *( 1.0d0 -s2 ) &
                &  +4.0d0 *coeff_WC *s2 *s /( 1.0d0 +ctmp2 )

           ctmp3 = kappa_wc +x

           f = 1.0d0 +kappa_wc *x / ctmp3
           df_dx = kappa_wc**2 /ctmp3**2           
           df_ds = dx_ds *df_dx

        case (5)
           ctmp1 = ( s -s_low_htbs )/( s_high_htbs -s_low_htbs )
           ctmp2 = ctmp1 *ctmp1
           ctmp3 = ctmp1 *ctmp2
           ctmp4 = ctmp1 *ctmp3
           ctmp5 = ctmp1 *ctmp4
           f = coeff_spline_htbs(0) +coeff_spline_htbs(1) *ctmp1 &
                &                   +coeff_spline_htbs(2) *ctmp2 &
                &                   +coeff_spline_htbs(3) *ctmp3 &
                &                   +coeff_spline_htbs(4) *ctmp4 &
                &                   +coeff_spline_htbs(5) *ctmp5
!
           df_ds = coeff_spline_htbs(1) + 2.0d0 *coeff_spline_htbs(2) *ctmp1 &
                &                       + 3.0d0 *coeff_spline_htbs(3) *ctmp2 &
                &                       + 4.0d0 *coeff_spline_htbs(4) *ctmp3 &
                &                       + 5.0d0 *coeff_spline_htbs(5) *ctmp4
           df_ds = df_ds / (s_high_htbs -s_low_htbs )

        case (6)                       ! pbesol
           x = mu_pbesol *s2;   dx_ds = 2.0d0 *mu_pbesol *s
           ctmp1 = kappa_pbesol + x

           f = 1.0d0 +kappa_pbesol *x /ctmp1
           df_dx = kappa_pbesol**2 /ctmp1**2
           df_ds = dx_ds *df_dx

        case (7)                       ! pbeint
           ctmp1 = mu_pbe - mu_gel

           alpha_pbeint = mu_gel**2 /kappa_pbeint /ctmp1
           ctmp2 = alpha_pbeint *s2
           ctmp3 = 1.0d0 +ctmp2

           mu_this = mu_gel +ctmp1 *ctmp2 /ctmp3
           dmu_ds = ctmp1 *alpha_pbeint *s  /ctmp3**2

           x = mu_this *s2;   dx_ds = 2.0d0 *mu_this *s + dmu_ds *s2

           ctmp4 = kappa_pbeint + x

           f = 1.0d0 +kappa_pbeint *x /ctmp4
           df_dx = kappa_pbeint**2 /ctmp4**2
           df_ds = dx_ds *df_dx

        case (11)                      ! B86r  I. Hamada
           x = mu_b86r *s2;   dx_ds = 2.0d0 *mu_b86r *s

           ctmp1 = 1.0d0 + x /kappa_b86r
           ctmp2 = ctmp1**(-4.0d0/5.0d0)
           ctmp3 = -0.8d0 *ctmp2 /ctmp1

           f = 1.0d0 + x *ctmp2
           df_dx = ctmp2 +ctmp3 *x /kappa_b86r
           df_ds = dx_ds *df_dx

        case (12)                      ! optpbe
           x = mu_optpbe *s2;   dx_ds = 2.0d0 *mu_optpbe *s
           ctmp1 = kappa_optpbe + x

           f = 1.0d0 +kappa_optpbe *x /ctmp1
           df_dx = kappa_optpbe**2 /ctmp1**2
           df_ds = dx_ds *df_dx

        case (13)                      ! optb86b
           x = mu_optb86b *s2;   dx_ds = 2.0d0 *mu_optb86b *s
           ctmp1 = 1.0d0 + x /kappa_optb86b
           ctmp2 = ctmp1**(-4.0d0/5.0d0)
           ctmp3 = -0.8d0 *ctmp2 /ctmp1

           f = 1.0d0 + x *ctmp2
           df_dx = ctmp2 +ctmp3 *x /kappa_optb86b
           df_ds = dx_ds *df_dx

        case (14)                      ! pw86r
           s4 = s2 *s2
           ctmp1 = a_pw86r
           ctmp2 = b_pw86r *s2
           ctmp3 = c_pw86r *s4

           ctmp4 = 1.0d0 +( 15.0d0 *ctmp1 +ctmp2 +ctmp3 ) *s2
           ctmp5 = 30.0d0 *ctmp1 +4.0d0 *ctmp2  +6.0d0 *ctmp3

           f = ctmp4 ** (1.0d0/15.0d0)
           df_ds = 1.0d0 /15.0d0 *f /ctmp4 * ctmp5 *s

        case (15)                      ! c09x
           x = mu_c09 *s2;  dx_ds = 2.0d0 *mu_c09 *s
           ctmp1 = alpha_c09 /mu_c09

           ctmp2 = exp( -ctmp1 *x )
           ctmp3 = exp( -ctmp1 *x /2.0d0 )

           ctmp4 = -ctmp1 *ctmp2
           ctmp5 = -ctmp1 *ctmp3 /2.0d0

           f = 1.0d0 + x *ctmp2 + kappa_c09 *( 1.0d0 -ctmp3 )
           df_dx = ctmp2 +x*ctmp4 -kappa_c09* ctmp5
           df_ds = dx_ds *df_dx

        case (16)                      ! lv-pw86r
           s4 = s2 *s2;  s6 = s4 *s2
           ctmp1 = 1.0d0 +15.d0 *a_pw86r *s2 +b_pw86r *s4 + c_pw86r *s6
           ctmp2 = ctmp1**(1.0d0/15.0d0)

           ctmp3 = 1.0d0/15.d0 /ctmp1 *ctmp2 *s &
                &  *( 30.0d0 *a_pw86r +4.0d0 *b_pw86r *s2 + 6.0d0 *c_pw86r *s4 )

           ctmp4 = 1.0d0 +alpha_lvpw86r *s6
           ctmp5 = beta_lvpw86r + alpha_lvpw86r *s6
           ctmp6 = 1.0d0 +mu_lvpw86r *s2

           f = ctmp6 /ctmp4 +alpha_lvpw86r *s6 /ctmp5 *ctmp2
           df_ds = s /ctmp4**2 *( -6.0d0 *alpha_lvpw86r *s4 *ctmp6 &
                &                 +2.0d0 *mu_lvpw86r *ctmp4 ) &
                &  +6.0d0 *alpha_lvpw86r *beta_lvpw86r *s4 *s /ctmp5**2 *ctmp2 &
                &  +alpha_lvpw86r *s6 /ctmp5 *ctmp3

        case (20)                    ! Engel Vosko (EV93)
           s4 = s2 *s2;  s6 = s4 *s2

           x = 1.0D0 +a1_ev93 *s2 +a2_ev93 *s4 +a3_ev93 *s6
           y = 1.0D0 +b1_ev93 *s2 +b2_ev93 *s4 +b3_ev93 *s6
           dx_ds = s *( 2.0D0 *a1_ev93 +4.0D0 *a2_ev93 *s2 +6.0D0 *a3_ev93 *s4 )
           dy_ds = s *( 2.0D0 *b1_ev93 +4.0D0 *b2_ev93 *s2 +6.0D0 *b3_ev93 *s4 )

           ctmp1 = y*dx_ds -x*dy_ds

           f = x /y
           df_ds = ctmp1 /y**2

        case (25)                     ! Leeuwen Baerends (LB94)
           f = 1.0d0;       df_ds = 0.0d0
           ctmp1 = d **(1./3.)
           if (d > 1.d-05) then
              x = dd /d /ctmp1
           else
              x = 0.d0
           endif
           ctmp2 = beta_lb94 *x**2 *ctmp1
           ctmp3 = 1.0d0 +3.d0 *beta_lb94 *x *asinh(x)
           pot_add = -ctmp2 /ctmp3

        case default
           call phase_error_with_msg(6,"xc functinal not supported ",__LINE__,__FILE__)
        end select

        ex  = fac *f *d
        exd = thrd4 *fac* ( f -s *df_ds )
        exdd = ax *df_ds *0.5d0 /thpith
!------------------------------------------
        exc0 = ex / facw
        excd = exd
        dFx_drho(i, is) = excd
        if ( pot_type == 25 ) dFx_drho(i,is) = dFx_drho(i,is) +pot_add

        excdd = exdd
! gradient of charge density 95/12/3 H.S.
        if(dabs(grad_rho(i, is)) > 1.d-9) then
           dFx_dgradrho(i, is) = excdd / grad_rho(i, is)
        else
           dFx_dgradrho(i, is) = 0.d0
        endif
! gradient of charge density 95/12/3 H.S.
        exc1 = exc1 + exc0*wos(i)
#ifdef __EDA__
        if(sw_eda==ON) then
! -----  ascat starts modifying  -----
        exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + exc0!*f2or1(i)
! -----  ascat ceases modifying  -----
        endif
#endif
     end do
     exc = exc + exc1
  end do

end subroutine ex_gga_library

#ifdef __EDA__
subroutine cr_gga_library( nspin, ispin, ista_r, iend_r, chgrhr_l, grad_trho, wos, &
     &                     exc, dF_drho, exc_on_a_grid_wk, ecor, pot_type, ist, ien )
#else
subroutine cr_gga_library( nspin, ispin, ista_r, iend_r, chgrhr_l, grad_trho, wos, &
     &                     exc, dF_drho, ecor, pot_type, ist, ien )
#endif
  use m_Const_Parameters,  only : PAI,DP
#ifdef __EDA__
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_eda
#endif
  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r, pot_type
  integer,intent(in)        :: ist, ien
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(inout) :: grad_trho(ista_r:iend_r)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)  
  real(kind=DP),intent(inout) :: exc
  real(kind=DP),intent(inout) :: dF_drho(ista_r:iend_r,nspin)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),intent(inout):: exc_on_a_grid_wk(ista_r:iend_r)
! -----  ascat ceases modifying  -----
#endif
  real(kind=DP),intent(out),optional :: ecor

  real(kind=DP), parameter :: zeta_minimum    = 1.d-4
  real(kind=DP), parameter :: zeta_minimum2   = 1.d-6
  real(kind=DP), parameter :: density_minimum = 1.d-20
  real(kind=DP), parameter :: density_minimum2= 1.d-5

  real(kind=DP), parameter  :: a  = 0.0310907d0,    a1   = 0.21370d0
  real(kind=DP), parameter  :: b1 = 7.5957d0,       b2n  = 3.5876d0
  real(kind=DP), parameter  :: b3 = 1.6382d0,       b4   = 0.49294d0
  real(kind=DP), parameter  :: p  = 1.00d0,         p1   = p+1.d0

  real(kind=DP), parameter  :: ap  = 0.01554535d0,    a1p  = 0.20548d0
  real(kind=DP), parameter  :: b1p = 14.1189d0,     b2np = 6.1977d0
  real(kind=DP), parameter  :: b3p = 3.3662d0,      b4p  = 0.62517d0

  real(kind=DP), parameter  :: aq  = 0.0168869d0,    a1q  = 0.11125d0
  real(kind=DP), parameter  :: b1q = 10.3570d0,     b2nq = 3.6231d0
  real(kind=DP), parameter  :: b3q = 0.88026d0,     b4q  = 0.49671d0

  real(kind=DP), parameter  :: gamma = 0.03109069086965489503494086371273d0
!--------------------------

  real(kind=DP), parameter  :: thrd = 0.333333333333d0, sixth7 = 1.1666666666666666d0

  real(kind=DP), parameter  :: eta = 1.d-12
  real(kind=DP), parameter  :: sixthm = -0.16666666666666d0

! ----- parameters 
  real(kind=DP), parameter :: beta_pbe   = 0.06672455060314922d0
  real(kind=DP), parameter :: beta_pbesol = 0.046d0
  real(kind=DP), parameter :: beta_pbeint = 0.052d0

  integer       :: is,i
  integer       :: istart, iend
  real(kind=DP) :: facw,exc1,g,g3,facpon,d,dd,rs,fk,sk,t,q0 &
       &         , rs12,rs32,rsp,q1,q2,eu,q3,eurs,ecd,pon &
       &         , b,b2,t2,t4,t6,q4,q5,h0,h,q8,h0t,h0b,h0rs,ht,hrs &
       &         , ec1,ec1d,ec1dd,exc0,excd,excdd &
       &         , thrd2,thrd4,fzdd,zeta,onpzeta,onmzeta,fz,zetadxd &
       &         , q0p,q1p,q2p,q3p,q0q,q1q,q2q,q3q,fzd,onzt4,gzt,bzt,h0zt,euzt
  real(kind=DP) :: bet, delt

! ------
  istart = ist;  iend = ien

  select case ( pot_type )
  case (6)                           ! pbesol
     bet = beta_pbesol
  case (7)                           ! pbeint
     bet = beta_pbeint
  case default
     bet = beta_pbe
  end select

  delt = bet /gamma
! ------


  facw = ispin

!---- Spin dependency
#ifdef NEC_TUNE_MXCP
!CDIR NOCONCUR
#endif
  if(present(ecor)) ecor = 0.d0
  do is = 1, ispin
     exc1 = 0.d0
     if ( ispin == 1 ) then
        g   = 1.d0
        g3  = g**3
!!$        g4  = g3*g
        facpon = -delt/(g3*bet)

#ifdef NEC_TUNE_MXCP
!CDIR INNER
#endif
        do i = istart, iend ! MPI
           d  = facw*chgrhr_l(i, 1)
           if(d < density_minimum) cycle
! gradient of charge density 95/12/3 H.S.
           dd = facw* grad_trho(i)

           rs = (0.75d0/(PAI*d))**thrd
           fk = (3.d0*PAI*PAI*d)**thrd
           sk = dsqrt(4.d0*fk/PAI)
           if(d > density_minimum2) then
              t = dd/(d*sk*2.d0)
           else
              t = 0.d0
           endif

           q0 = -2.d0 *a*(1.d0 + a1*rs)
           rs12 = dsqrt(rs)
           rs32 = rs12**3
           rsp = rs**p
           q1 = 2.d0 *a*(b1*rs12+b2n*rs+b3*rs32+b4*rs*rsp)
           q2 = log( 1.d0 + 1.d0/q1)
           eu = q0*q2
           q3 = a*(b1/rs12 +2.d0*b2n +3.d0*b3*rs12 +2.d0*b4*p1*rsp)
           eurs = -2.d0 *a*a1*q2 - q0*q3/(q1**2 + q1)
           ecd = eu-thrd*rs*eurs
           pon = facpon*eu
           b   = delt/(dexp(pon)-1.d0)
           b2  = b*b
           t2  = t*t

           t4  = t2*t2
           t6  = t4*t2
           q4  = 1.d0 + b*t2
           q5  = 1.d0 + b*t2+b2*t4
           h0  = g3*(bet/delt)*dlog(1d0+delt*q4*t2/q5)
           h   = h0
           q8  = q5*q5+delt*q4*q5*t2
           h0t = 2.d0 *bet*t*( 1.d0 +2d0 *b*t2)/q8 * g3
           h0b = -bet*t6*( 2.d0 *b+b2*t2)/q8 * g3
           h0rs = h0b*b*eurs*(b+delt)/ bet / g3
           ht   = h0t
           hrs  = h0rs
           ec1  = d*h/facw
           ec1d = h-thrd*rs*hrs-sixth7*t*ht

           ec1dd = 0.5d0*ht/sk

           exc0 = eu*d/facw + ec1
           excd = ecd + ec1d
           dF_drho(i, is) = dF_drho(i, is) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(dabs(grad_trho(i)) > 1.0d-9) then
              grad_trho(i) = excdd / grad_trho(i)
           else
              grad_trho(i) = 0.d0
           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d/facw
           exc1 = exc1 + exc0*wos(i)
#ifdef __EDA__
           if(sw_eda==ON) then
! -----  ascat starts modifying  -----
           exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + exc0!*f2or1(i)
! -----  ascat ceases modifying  -----
           endif
#endif
        end do

     else if ( ispin ==  2 ) then
        thrd2 = thrd * 2.d0
        thrd4 = thrd * 4.d0
        fzdd = 8.d0 /(9.d0*(2.d0**thrd4 - 2.d0 ))

        do i = ista_r,iend_r
           d = chgrhr_l(i, 1) + chgrhr_l(i, ispin)
           if(d < density_minimum ) cycle
           zeta = ( chgrhr_l(i, 1) - chgrhr_l(i, 2) ) / d
           onpzeta = 2.d0 *chgrhr_l(i, 1) / d
           onmzeta = 2.d0 *chgrhr_l(i, 2) / d
           g = ( onpzeta**thrd2 + onmzeta**thrd2 )/2.d0
           if(dabs(zeta) < zeta_minimum ) then
              fz= 4.d0/9.d0 *zeta*zeta*(1.d0 + 5.d0/54.d0 *zeta*zeta &
                   &               *(1.d0 +44.d0/135.d0*zeta*zeta )) &
                   &       /(2.d0**thrd4 -2.d0)
           else
              fz = (onpzeta**thrd4 + onmzeta**thrd4 -2.d0) / (2.d0**thrd4 - 2.d0)
           end if
           g3  = g**3
!!$           g4  = g3*g
           facpon = -delt/(g3*bet)
           if(d < density_minimum) cycle
           dd = grad_trho(i)
           rs = (0.75d0/(PAI*d))**thrd
           fk = (3.d0*PAI*PAI*d)**thrd
           sk = dsqrt(4.d0*fk/PAI)

           if(d > density_minimum2) then
              t = dd/(d*sk*2.d0*g)
           else
              t = 0.d0
           endif

           q0 = -2.d0 *a*(1.d0+a1*rs)
           rs12 = dsqrt(rs)
           rs32 = rs12**3
           rsp = rs**p
           q1 = 2.d0 *a*(b1*rs12+b2n*rs+b3*rs32+b4*rs*rsp)
           q2 = log( 1.d0 +1.d0/q1 )
           q3 = a*(b1/rs12 +2.d0*b2n +3.d0*b3*rs12 +2.d0*b4*p1*rsp)

           q0p = -2.d0 *ap*(1.d0+a1p*rs)
           q1p = 2.d0*ap*(b1p*rs12+b2np*rs+b3p*rs32+b4p*rs*rsp)
           q2p = log(1.d0 +1.d0/q1p)
           q3p = ap*(b1p/rs12 +2.d0*b2np +3.d0*b3p*rs12 +2.d0*b4p*p1*rsp)
! vwn 95/12/2 Y.M
           q0q = -2.d0*aq*(1.d0+a1q*rs)
           q1q =  2.d0*aq*(b1q*rs12+b2nq*rs+b3q*rs32+b4q*rs*rsp)
           q2q = log(1.d0 +1.d0/q1q)
           q3q = aq*(b1q/rs12 +2.d0*b2nq +3.d0*b3q*rs12 +2.d0*b4q*p1*rsp)
           if(dabs(zeta) < zeta_minimum) then
              fzd = 8.d0/9.d0 /( 2.d0**thrd4 -2.d0 ) &
                   &     *zeta *( 1.d0 + 5.d0/27.d0 *zeta*zeta &
                   &           *( 1.d0 +22.d0/45.d0 *zeta*zeta ))
           else
              fzd = 4.d0/3.d0 /( 2.d0**thrd4 -2.d0) &
                   &     *( onpzeta**thrd -onmzeta**thrd )
           end if
           onzt4= onmzeta*onpzeta *( 1.d0 + zeta*zeta )
           eu = q0*q2 + ( q0p * q2p - q0 * q2 ) * fz * zeta**4 &
                &             -   q0q * q2q * fz / fzdd * onzt4
           eurs = ( -2.d0*a*a1*q2 -q0*q3/(q1**2 +q1 )) * ( 1.d0 -fz*zeta**4 )&
                &+( -2.d0*ap*a1p*q2p -q0p*q3p/( q1p**2 +q1p )) *fz*zeta**4 &
                &-( -2.d0*aq*a1q*q2q -q0q*q3q/( q1q**2 +q1q )) *fz/fzdd*onzt4
           euzt =  -q0q*q2q /fzdd *( fzd*onzt4 -4.d0*fz*zeta**3 )&
                &        +( q0p*q2p -q0*q2 )*(fzd*zeta +4.d0*fz )*zeta**3 
! vwn 95/12/2 Y.M

           zetadxd = - 2.d0 * chgrhr_l(i, 3-is) * (-1.d0)**is / d
           ecd = eu-thrd*rs*eurs + euzt * zetadxd 

           pon = facpon*eu
           b   = delt/(dexp(pon)-1.d0)
           b2  = b*b
           t2  = t*t

           t4  = t2*t2
           t6  = t4*t2
           q4  = 1.d0 +b*t2
           q5  = 1.d0 +b*t2+b2*t4
           h0  = g3*(bet/delt)*dlog(1.d0+delt*q4*t2/q5)
           h   = h0

           q8  = q5*q5+delt*q4*q5*t2
           h0t = 2.d0 *bet*t*(1.d0 +2.d0*b*t2)/q8 * g3
           h0b = -bet*t6*( 2.d0*b +b2*t2)/q8 * g3
           h0rs = h0b * b * eurs * (b+delt) / bet / g3
           ht   = h0t
           hrs  = h0rs

           if(dabs(zeta) < zeta_minimum) then
              gzt = -2.d0/9.d0*zeta *( 1.d0 +14.d0/27.d0 *zeta*zeta &
                   &             *( 1.d0 +13.d0/18.d0 *zeta*zeta ))
           else
              gzt = ( (onpzeta**2+eta)**sixthm - (onmzeta**2+eta)**sixthm)/3.d0
           end if
           bzt = b * ( b+delt) * ( euzt - 3.d0*eu*gzt/g ) / ( bet * g3 ) 
           h0zt = 3.d0 * gzt * h0 / g + h0b * bzt


           ec1  = d*h * 0.5d0

           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxd
           ec1dd = 0.5d0*ht/sk / g

10         continue

           exc0 = eu*d * 0.5d0 + ec1
           excd = ecd + ec1d
           dF_drho(i, is) = dF_drho(i, is) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(is ==  2) then
              if(dabs(grad_trho(i)) > 1.0d-9) then
                 grad_trho(i) = excdd / grad_trho(i)
              else
                 grad_trho(i) = 0.d0
              endif
           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d * 0.5
           exc1 = exc1 + exc0*wos(i)
#ifdef __EDA__
           if(sw_eda==ON) then
! -----  ascat starts modifying  -----
           exc_on_a_grid_wk(i) = exc_on_a_grid_wk(i) + exc0!*f2or1(i)
! -----  ascat ceases modifying  -----
           endif
#endif
        end do
     end if
     exc = exc + exc1
     if(present(ecor)) ecor = ecor + exc1
  end do
end subroutine cr_gga_library

subroutine init_spline_htbs( mu, kappa_rpbe, kappa_wc, coeff_wc, &
     &                       s_low_htbs, s_high_htbs, coeff_spline )
  use m_Const_Parameters,  only : DP
  implicit none

  real(kind=DP), intent(in) :: mu, kappa_rpbe, kappa_wc, coeff_wc
  real(kind=DP), intent(in) :: s_low_htbs, s_high_htbs
  real(kind=DP), intent(out) :: coeff_spline(0:5)
!
  real(kind=DP) :: f_0, df_0, d2f_0, f_1, df_1, d2f_1
  real(kind=DP) :: delta_s, ctmp1, ctmp2, ctmp3

  call set_values_endpoint0( s_low_htbs,  f_0, df_0, d2f_0 )   ! WC06 at lower limit
  call set_values_endpoint1( s_high_htbs, f_1, df_1, d2f_1 )   ! RPBE at high limit

! scaling from s=(s1,s2) to y=(0, 1)   ! y = (s-s1)/(s2-s1)
!
  delta_s = s_high_htbs -s_low_htbs
  df_0 = df_0 *delta_s;     d2f_0 = d2f_0 *delta_s**2
  df_1 = df_1 *delta_s;     d2f_1 = d2f_1 *delta_s**2
!
  coeff_spline(0) = f_0
  coeff_spline(1) = df_0
  coeff_spline(2) = d2f_0 /2.0d0
!
  ctmp1 =   f_1 -coeff_spline(0) -coeff_spline(1) -coeff_spline(2)
  ctmp2 =  df_1 -coeff_spline(1) -2.0d0 *coeff_spline(2)
  ctmp3 = d2f_1 -2.0d0 *coeff_spline(2)
!
  coeff_spline(3) = ( 20.0d0 *ctmp1 -8.0d0 *ctmp2 + ctmp3 ) /2.0d0
  coeff_spline(4) = ( 4.0d0 *ctmp2 -ctmp3 -6.0d0 *coeff_spline(3) ) /4.0d0
  coeff_spline(5) = ctmp1 -coeff_spline(3) -coeff_spline(4)

contains

  subroutine set_values_endpoint0( s, f, df_ds, d2f_ds2 )
    real(kind=DP), intent(in) :: s
    real(kind=DP), intent(out) :: f, df_ds, d2f_ds2

    real(kind=DP), parameter :: mu_gel    = 0.1234567901234568d0   !  = 10.0 /81.0

    real(kind=DP) :: s2, s4, x, dx_ds, d2x_ds2, df_dx, d2f_dx2
    real(kind=DP) :: ctmp1, ctmp2, ctmp3

    s2 = s**2;  s4 = s2**2

    ctmp1 = exp( -s2 )
    ctmp2 = coeff_WC *s4

    x = mu_gel *s2 + ( mu -mu_gel )*s2 *ctmp1 + log( 1.0d0 +ctmp2 )
    dx_ds = 2.0d0 *mu_gel *s + ( mu -mu_gel ) *2.0d0 *s *ctmp1 *( 1.0d0 -s2 ) &
         &  +4.0d0 *coeff_WC *s2 *s /( 1.0d0 +ctmp2 )
    d2x_ds2 = 2.0d0 *mu_gel &
         &   + 2.0d0 *( mu -mu_gel ) *ctmp1 *( 1.0d0 -5.0d0 *s2 +2.0d0 *s4 ) &
         &   + 4.0d0 *coeff_WC *s2 *( 3.0d0 -ctmp2 ) /( 1.0d0 +ctmp2 )**2

    ctmp3 = kappa_wc +x

    f = 1.0d0 +kappa_wc *x /ctmp3
    df_dx = kappa_wc**2 /ctmp3**2
    d2f_dx2 = -2.0d0 *df_dx /ctmp3

    df_ds = dx_ds * df_dx
    d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2

  end subroutine set_values_endpoint0

  subroutine set_values_endpoint1( s, f, df_ds, d2f_ds2 )
    real(kind=DP), intent(in) :: s
    real(kind=DP), intent(out) :: f, df_ds, d2f_ds2

    real(kind=DP) :: s2, x, dx_ds, d2x_ds2
    real(kind=DP) :: ctmp1

    s2 = s**2

    x = mu *s2;   dx_ds = 2.0d0 *mu *s;   d2x_ds2 = 2.0d0 *mu
    ctmp1 = exp( -x/ kappa_rpbe )

    f = 1.0d0 +kappa_rpbe -kappa_rpbe *ctmp1
    df_ds = dx_ds *ctmp1
    d2f_ds2 = ctmp1 * ( d2x_ds2 -dx_ds**2 /kappa_rpbe )

  end subroutine set_values_endpoint1

end subroutine init_spline_htbs
! ===================== 13.0A
