
module m_Fx_omega_PBE
 use m_Const_Parameters,  only : DP,PAI,PAIsqrt
 implicit none

 real(kind=DP), parameter :: &
                 ca =  1.0161144d0 &
               , cb = -0.37170836d0 &
               , cc = -0.077215461d0 &
               , cd =  0.57786348d0 &
               , ce = -0.051955731d0 &
               , ck = -4.d0/27.d0

 real(kind=DP), parameter :: &
                 a1 = -0.000205484d0 &
               , a2 = -0.109465240d0 &
               , a3 = -0.064078780d0 &
               , a4 = -0.008181735d0 &
               , a5 = -0.000110666d0

 real(kind=DP), parameter :: &
                 b1 =  0.006601306d0 &
               , b2 =  0.259931140d0 &
               , b3 =  0.520352224d0 &
               , b4 =  0.118551043d0 &
               , b5 =  0.046003777d0

 real(kind=DP), parameter :: pi = PAI
 real(kind=DP), parameter :: rpi = PAIsqrt

contains

 function Fx_omega(fk,s,omega)
  implicit none
  real(kind=DP) :: Fx_omega
  real(kind=DP), intent(in) :: fk,s,omega

  real(kind=DP) :: s2, alp1, alp2, bet

  s2 = s*s
  alp1 = s2 * h(s)
  alp2 = alp1 + cd
  bet = omega/fk

  Fx_omega = a1 * yjw1(alp1+b1,bet) &
           + a2 * yjw1(alp1+b2,bet) &
           + a3 * yjw2(alp1+b3,bet) &
           + a4 * yjw2(alp1+b4,bet) &
           + a5 * yjw3(alp1+b5,bet) &
           + cb * yjw1(alp2,bet) &
           + cc * (1.d0+s2*f(s)) * yjw3(alp2,bet) &
           +      (ce+s2eg(s)) * yjw5(alp2,bet)
  Fx_omega = -8.d0/9.d0 * Fx_omega
 end function Fx_omega
 
 function dFxds_omega(fk,s,omega)
  implicit none
  real(kind=DP) :: dFxds_omega
  real(kind=DP), intent(in) :: fk,s,omega

  real(kind=DP) :: s2, alp1, alp2, dalp, bet
  real(kind=DP) :: ff, gg
  real(kind=DP) :: dff, dgg

  s2 = s*s
  alp1 = s2 * h(s)
  alp2 = alp1 + cd 
  dalp = 2.d0*s*h(s) + s2*dh(s)
  bet = omega/fk
  ff = cc * (1.d0+s2*f(s))
  gg = ce+s2eg(s)
  dff = cc * ( 2.d0*s*f(s) + s2*df(s) )
  dgg = 2.d0*s*eg(s) + s2*deg(s)

  dFxds_omega = ( a1 * dyjw1a(alp1+b1,bet) &
                + a2 * dyjw1a(alp1+b2,bet) &
                + a3 * dyjw2a(alp1+b3,bet) &
                + a4 * dyjw2a(alp1+b4,bet) &
                + a5 * dyjw3a(alp1+b5,bet) &
                + cb * dyjw1a(alp2,bet) &
                + ff * dyjw3a(alp2,bet) &
                + gg * dyjw5a(alp2,bet) ) * dalp &
                + dff * yjw3(alp2,bet) &
                + dgg * yjw5(alp2,bet) 
  dFxds_omega = -8.d0/9.d0 * dFxds_omega
 end function dFxds_omega

 function kdFxdk_omega(fk,s,omega)
  implicit none
  real(kind=DP) :: kdFxdk_omega
  real(kind=DP), intent(in) :: fk,s,omega

  real(kind=DP) :: s2, alp1, alp2, bet

  s2 = s*s
  alp1 = s2 * h(s)
  alp2 = alp1 + cd
  bet = omega/fk

  kdFxdk_omega = a1 * dyjw1b(alp1+b1,bet) &
              + a2 * dyjw1b(alp1+b2,bet) &
              + a3 * dyjw2b(alp1+b3,bet) &
              + a4 * dyjw2b(alp1+b4,bet) &
              + a5 * dyjw3b(alp1+b5,bet) &
              + cb * dyjw1b(alp2,bet) &
              + cc * (1.d0+s2*f(s)) * dyjw3b(alp2,bet) &
              +      (ce+s2eg(s)) * dyjw5b(alp2,bet)
  kdFxdk_omega = 8.d0/9.d0 * bet * kdFxdk_omega
 end function kdFxdk_omega

 function k2dFxdk_omega(fk,s,omega)
  implicit none
  real(kind=DP) :: k2dFxdk_omega
  real(kind=DP), intent(in) :: fk,s,omega

  real(kind=DP) :: s2, alp1, alp2, bet

  s2 = s*s
  alp1 = s2 * h(s)
  alp2 = alp1 + cd
  bet = omega/fk

  k2dFxdk_omega = a1 * dyjw1b(alp1+b1,bet) &
              + a2 * dyjw1b(alp1+b2,bet) &
              + a3 * dyjw2b(alp1+b3,bet) &
              + a4 * dyjw2b(alp1+b4,bet) &
              + a5 * dyjw3b(alp1+b5,bet) &
              + cb * dyjw1b(alp2,bet) &
              + cc * (1.d0+s2*f(s)) * dyjw3b(alp2,bet) &
              +      (ce+s2eg(s)) * dyjw5b(alp2,bet)
  k2dFxdk_omega = 8.d0/9.d0 * omega * k2dFxdk_omega
 end function k2dFxdk_omega

 function f(s)
  implicit none
  real(kind=DP) :: f
  real(kind=DP), intent(in) :: s
 
  f = ( h(s) * (16.d0*ca**2 + 36.d0*(cb-ca*cd)) + 9.d0*ck ) / (36.d0*cc)
 end function f

 function df(s)
  implicit none
  real(kind=DP) :: df
  real(kind=DP), intent(in) :: s
 
  df = ( dh(s) * (16.d0*ca**2 + 36.d0*(cb-ca*cd)) ) / (36.d0*cc)
 end function df

 function eg(s)
  implicit none
  real(kind=DP) :: eg
  real(kind=DP), intent(in) :: s

  real(kind=DP) :: a,b,s2,t,d

  s2 = s * s

  t = cd + h(s) * s2
  d = 16.d0 * t**3.5d0

  a = sqrt(pi) * ( 15.d0*ce + t * ( 6.d0*cc*(1.d0+f(s)*s2) &
                  + 4.d0*cb*t + 8.d0*ca*t*t ) ) / d  &
      - 3.d0*pi*sqrt(ca)/4.d0 &
               * exp(9.d0*h(s)*s2/(4.d0*ca)) * erfc(1.5d0*s*sqrt(h(s)/ca))

  b = 15.d0*sqrt(pi)*s2 / d

  eg = - (3.d0*pi/4.d0+a) / b
 end function eg

 function s2eg(s)
  implicit none
  real(kind=DP) :: s2eg
  real(kind=DP), intent(in) :: s

  real(kind=DP) :: a,b,s2,t,d

  s2 = s * s

  t = cd + h(s) * s2
  d = 16.d0 * t**3.5d0

  a = sqrt(pi) * ( 15.d0*ce + t * ( 6.d0*cc*(1.d0+f(s)*s2) &
                  + 4.d0*cb*t + 8.d0*ca*t*t ) ) / d  &
      - 3.d0*pi*sqrt(ca)/4.d0 &
               * exp(9.d0*h(s)*s2/(4.d0*ca)) * erfc(1.5d0*s*sqrt(h(s)/ca))

  b = 15.d0*sqrt(pi) / d

  s2eg = - (3.d0*pi/4.d0+a) / b
 end function s2eg


 function deg(s)
  implicit none
  real(kind=DP) :: deg
  real(kind=DP), intent(in) :: s

  real(kind=DP) :: a,b,b2,s2,t,u,d
  real(kind=DP) :: da,db,dt
  real(kind=DP) :: rh

  s2 = s * s
  t = cd + h(s) * s2
  u = 1.d0 + f(s) * s2
  d = 16.d0 * t**3.5d0

  a = sqrt(pi) * ( 15.d0*ce + t * ( 6.d0*cc*(1.d0+f(s)*s2) &
                  + 4.d0*cb*t + 8.d0*ca*t*t ) ) / d  &
      - 3.d0*pi*sqrt(ca)/4.d0 &
               * exp(9.d0*h(s)*s2/(4.d0*ca)) * erfc(1.5d0*s*sqrt(h(s)/ca))
  b = 15.d0*sqrt(pi)*s2 / d
  b2 = b * b

  dt = dh(s) * s2 + 2.d0 * s * h(s)
  rh = sqrt(h(s))
  da = sqrt(pi) * ( 6.d0*cc*( df(s)*s2 + 2.d0*s*f(s) ) * t &
     &            + ( 6.d0*cc*u + 8.d0*cb*t+24.d0*ca*t*t ) * dt ) &
     &            / d &
     & - 7.d0*sqrt(pi)/32.d0 * ( 15.d0*ce + 6.d0*cc*u*t + 4.d0*cb*t*t + 8.d0*ca*t**3 ) &
     &                      * dt / t**4.5d0 &
     & - 27.d0*pi/(16.d0*sqrt(ca)) * dt * exp(9.d0*h(s)*s2/(4.d0*ca)) * erfc(1.5d0*s*sqrt(h(s)/ca)) &           
     & + 9.d0*sqrt(pi)/4.d0*( rh + 0.5d0*s*dh(s)/rh )
  db = 15.d0*sqrt(pi) * ( 2.d0*s - 3.5d0*s2*dt/t ) / d

  deg = -da/b + db*(3.d0*pi/4.d0+a)/b2 
 end function deg

 function h(s)
  implicit none
  real(kind=DP) :: h
  real(kind=DP), intent(in) :: s

  real(kind=DP) :: s2, s4
  real(kind=DP), parameter :: a1 = 0.00979681d0 &
                           , a2 = 0.0410834d0 &
                           , a3 = 0.187440d0 &
                           , a4 = 0.00120824d0 &
                           , a5 = 0.0347188d0

  s2 = s * s
  s4 = s2 * s2

  h = ( a1 * s2 + a2 * s4 ) / ( 1.d0 + s4 * (a3 + a4 * s + a5 * s2 ) )
 end function h

 function dh(s)
  implicit none
  real(kind=DP) :: dh
  real(kind=DP), intent(in) :: s

  real(kind=DP) :: s2, s4, s5, s6
  real(kind=DP), parameter :: a1 = 0.00979681d0 &
                           , a2 = 0.0410834d0 &
                           , a3 = 0.187440d0 &
                           , a4 = 0.00120824d0 &
                           , a5 = 0.0347188d0

  s2 = s * s
  s4 = s2 * s2
  s5 = s * s4
  s6 = s2 * s4

  dh = -s*( a1*( 2.d0*a3*s4 + 3.d0*a4*s5 + 4.d0*a5*s6 - 2.d0) + a2*s2*(a4*s5+2.d0*a5*s6-4.d0)) &
     &    / ( 1.d0 + s4 * (a3 + a4 * s + a5 * s2 ) )**2
 end function dh

 function yjw1(a,b)
  implicit none
  real(kind=DP) :: yjw1
  real(kind=DP), intent(in) :: a, b
  
  yjw1 = (1.d0-b/sqrt(a+b*b))/(2.d0*a)
 end function yjw1

 function dyjw1a(a,b)
  implicit none
  real(kind=DP) :: dyjw1a
  real(kind=DP), intent(in) :: a, b

  real(kind=DP) :: apb2
  real(kind=DP) :: rapb2

  apb2 = a + b*b
  rapb2 = sqrt(apb2)

  dyjw1a = (b*(2.d0/rapb2+a/(apb2*rapb2))-2.d0)/(4.d0*a*a)
 end function dyjw1a

 function dyjw1b(a,b)
  implicit none
  real(kind=DP) :: dyjw1b
  real(kind=DP), intent(in) :: a, b

  real(kind=DP) :: apb2

  apb2 = a + b*b

  dyjw1b = -0.5d0/apb2**1.5d0
 end function dyjw1b

 function yjw2(a,b)
  implicit none
  real(kind=DP) :: yjw2
  real(kind=DP), intent(in) :: a, b

  real(kind=DP) :: ra, b2

  ra = sqrt(a)
  b2 = b * b

  yjw2 = (-ra*b+(a+b2)*acot(b/ra))/(2.d0*a**(3.d0/2.d0)*(a+b2)*rpi)
 end function yjw2

 function dyjw2a(a,b)
  implicit none
  real(kind=DP) :: dyjw2a
  real(kind=DP), intent(in) :: a, b
 
  real(kind=DP) :: ra, b2

  ra = sqrt(a)
  b2 = b * b

  dyjw2a = (ra*b*(5.d0*a+3.d0*b2)-3.d0*(a+b2)**2*acot(b/ra))/(4.d0*a**(5.d0/2.d0)*(a+b2)**2*rpi)
 end function dyjw2a

 function dyjw2b(a,b)
  implicit none
  real(kind=DP) :: dyjw2b
  real(kind=DP), intent(in) :: a, b

  real(kind=DP) :: ab2

  ab2 = a + b*b

  dyjw2b = -1.d0/(sqrt(pi)*ab2*ab2)
 end function dyjw2b

 function yjw3(a,b)
  implicit none
  real(kind=DP) :: yjw3
  real(kind=DP), intent(in) :: a, b

  real(kind=DP) :: b2, ab2, rab2

  b2 = b * b
  ab2 = a+b2
  rab2 = sqrt(ab2)

  yjw3 = (2.d0*b2*(-b+rab2)+a*(-3.d0*b+2.d0*rab2))/(4.d0*a*a*ab2*rab2)
 end function yjw3

 function dyjw3a(a,b)
  implicit none
  real(kind=DP) :: dyjw3a
  real(kind=DP), intent(in) :: a, b

  real(kind=DP) :: a2, b2, b4, ab2, rab2

  a2 = a * a
  b2 = b * b
  b4 = b2 * b2
  ab2 = a+b2
  rab2 = sqrt(ab2)

  dyjw3a = (a2*(15.d0*b-8.d0*rab2) + 8.d0*b4*(b-rab2) - 4.d0*a*b2*(4.d0*rab2-5.d0*b))/(8.d0*a2*a*ab2*ab2*rab2)
 end function dyjw3a

 function dyjw3b(a,b)
  implicit none
  real(kind=DP) :: dyjw3b
  real(kind=DP), intent(in) :: a, b

  real(kind=DP) :: ab2

  ab2 = a + b*b

  dyjw3b = -0.75d0/ab2**2.5d0
 end function dyjw3b

 function yjw5(a,b)
  implicit none
  real(kind=DP) :: yjw5
  real(kind=DP), intent(in) :: a, b

  real(kind=DP) :: a2, b2, b4, ab2, rab2

  a2 = a * a
  b2 = b * b
  b4 = b2 * b2
  ab2 = a+b2
  rab2 = sqrt(ab2)

  yjw5 = (8.d0*b4*(-b+rab2)+4.d0*a*b2*(-5.d0*b+4.d0*rab2)+a2*(-15.d0*b+8.d0*rab2))/(8.d0*a2*a*ab2*ab2*rab2)
 end function yjw5

 function dyjw5a(a,b)
  implicit none
  real(kind=DP) :: dyjw5a
  real(kind=DP), intent(in) :: a, b

  real(kind=DP) :: a2, a3, a4, b2, b4, b6, ab2, rab2

  a2 = a * a
  a3 = a * a2
  a4 = a2 * a2
  b2 = b * b
  b4 = b2 * b2
  b6 = b2 * b4
  ab2 = a+b2
  rab2 = sqrt(ab2)

  dyjw5a = -3.d0*( 16.d0*b6*(rab2-b) + 8.d0*a*b4*(6.d0*rab2-7.d0*b) + a3*(16.d0*rab2-35.d0*b) &
                & + 2.d0*a2*b2*(24.d0*rab2 - 35.d0*b) ) / (16.d0*a4*ab2*ab2*ab2*rab2)
 end function dyjw5a

 function dyjw5b(a,b)
  implicit none
  real(kind=DP) :: dyjw5b
  real(kind=DP), intent(in) :: a, b

  real(kind=DP) :: ab2

  ab2 = a + b*b

  dyjw5b = -15.d0/(8.d0*ab2**3.5d0)
 end function dyjw5b

 function acot(a)
  implicit none
  real(kind=DP) :: acot
  real(kind=DP), intent(in) :: a
 
  acot = atan(1.d0/a)
 end function acot
 
end module m_Fx_omega_PBE

subroutine ex_omegapbe(amix,omega,nspin,ispin,ista_fftph,iend_fftph,chgrhr_l,grad_rho,f2or1,exc,dFx_drho,dFx_dgradrho,nfft_y)
  use m_Const_Parameters,  only : DP,PAI,PAI4
  use m_Fx_omega_PBE,      only : Fx_omega, dFxds_omega, k2dFxdk_omega
  implicit none

  real(kind=DP),intent(in)  :: amix,omega
  integer,intent(in)        :: nspin,ispin,nfft_y
  integer,intent(in)        :: ista_fftph, iend_fftph
  real(kind=DP),intent(in)  :: chgrhr_l(1:nfft_y,ispin)
  real(kind=DP),intent(in)  :: grad_rho(1:nfft_y,nspin)
  real(kind=DP),intent(in)  :: f2or1(1:nfft_y)  
  real(kind=DP),intent(inout) :: exc
  real(kind=DP),intent(inout) :: dFx_drho(1:nfft_y,nspin)
  real(kind=DP),intent(inout) :: dFx_dgradrho(1:nfft_y,nspin)

  real(kind=DP), parameter :: ax = -0.7385587663820224d0 !! ax = -(3/4)*(3/PI)^(1/3)
  real(kind=DP), parameter :: thrd = 0.33333333333d0, thrd4 = 1.333333333333333d0
  real(kind=DP), parameter :: th8pi = -0.1193662073189215d0

  real(kind=DP) :: facw,d,dd,fk,s,fac,f,ex,fs,ff,exd,exdd,exc0,exc1
  integer       :: is,i

!---- Spin dependency

  facw = ispin
#ifdef NEC_TUNE_MXCP
!CDIR NOCONCUR
#endif
  do is = 1, ispin
     exc1 = 0.d0
#ifdef NEC_TUNE_MXCP
!CDIR INNER
#endif
     do i = 1, nfft_y
        d  = facw * chgrhr_l(i, is)
        if(abs(d) < 1.d-30) d = 1.d-30
        fk = (3*PAI*PAI*d)**thrd
        dd = facw * grad_rho(i, is)
        if(abs(d) > 1.d-05) then
           s = dd/(d*fk*2)
        else
           s = 0.d0
        end if
!-------------------------------------
        f = Fx_omega(fk,s,omega)
        ff = k2dFxdk_omega(fk,s,omega)
        if(abs(s)>1.d-10) then
           fs = dFxds_omega(fk,s,omega)
        else
           fs = 0.d0
           !!write(2000,'("s=",f20.5," f=",f20.5," ff=",f20.5)') s, f, ff
        end if
        fac = ax*d**thrd
        ex  = fac*f*d
        exd = thrd4*fac*(f-s*fs)-ff/PAI4
        exdd = th8pi*fs
!------------------------------------------     
        exc0 = ex / facw
        dFx_drho(i, is) = dFx_drho(i, is) + amix * exd 
        if(dabs(grad_rho(i, is)) > 1.d-9) then
           dFx_dgradrho(i, is) = dFx_dgradrho(i, is) + amix * exdd / grad_rho(i, is)
        endif
        exc1 = exc1 + exc0*f2or1(i)
        !!write(1000,'("i,f,fs,ff,d,dd,fk,s=",i8,1x,7(1x,f20.5))') i,f,fs,ff,d,dd,fk,s
     end do
     exc = exc + amix * exc1
  end do
  !!stop 'Debug: OmegaPBE'
end subroutine ex_omegapbe

! ----------------------------------
! HJS
!
! J. Chem. Phys. 128 (2008) 194105.
! ----------------------------------
subroutine ex_omegapbe_library_HJS( amix, omega, &
     &                              nspin, ispin, ista_r, iend_r, chgrhr_l, grad_rho, &
     &                              wos, exc, dFx_drho, dFx_dgradrho, pot_type )
  use m_Const_Parameters,  only : DP,PAI
  implicit none

  real(kind=DP),intent(in)  :: amix, omega
  integer,intent(in)        :: nspin,ispin,ista_r,iend_r, pot_type
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)  :: grad_rho(ista_r:iend_r,nspin)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)
  real(kind=DP),intent(out) :: exc
  real(kind=DP),intent(out) :: dFx_drho(ista_r:iend_r,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(ista_r:iend_r,nspin)

  real(kind=DP), parameter :: ax = -0.7385587663820224d0 !! ax = -(3/4)*(3/PI)^(1/3)
  real(kind=DP), parameter :: thrd = 0.33333333333d0, thrd4 = 1.333333333333333d0
  real(kind=DP), parameter :: thpith = 3.0936677262801d0

  integer       :: is, i

  real(kind=DP) :: facw, d, dd, fac, ex, &
       &           exd, exdd, exc0, excd, excdd, exc1
  real(kind=DP) :: fx, dfx, ddfx

!---- Spin dependency

  facw = ispin
#ifdef NEC_TUNE_MXCP
!CDIR NOCONCUR
#endif
  do is = 1, ispin
     exc1 = 0.d0
#ifdef NEC_TUNE_MXCP
!CDIR INNER
#endif
     do i = ista_r, iend_r
        d  = facw * chgrhr_l(i, is)
!        if(abs(d) < 1.d-30) d = 1.d-30
!        if(abs(d) < 1.d-8) d = 1.d-8
        if(abs(d) < 1.d-12) d = 1.d-12

        dd = facw * grad_rho(i, is)
        fac = ax*d**thrd

        select case (pot_type)
        case (1)         ! PBE
           call omegapbe_analytic( 1, d, dd, omega, fx, dfx, ddfx )
        case (6)         ! PBE-sol
           call omegapbe_analytic( 2, d, dd, omega, fx, dfx, ddfx )
        end select

        ex  = fac *fx *d
        exd = thrd4 *fac* fx +fac *d *dfx
        exdd = fac * d *ddfx

        exc0 = ex / facw
        dFx_drho(i, is) = dFx_drho(i, is) + amix * exd 
        if(dabs(grad_rho(i, is)) > 1.d-9) then
           dFx_dgradrho(i, is) = dFx_dgradrho(i, is) + amix * exdd / grad_rho(i, is)
        endif
        exc1 = exc1 + exc0*wos(i)
     end do
     exc = exc + amix * exc1
  end do

end subroutine ex_omegapbe_library_HJS

subroutine omegapbe_analytic( pot_type, rho, gradrho, omega, fx, dfx_drho, dfx_dgradrho )
  use m_Const_Parameters,  only : DP,PAI
  implicit none
  
  integer, intent(in) :: pot_type
  real(kind=DP), intent(in) :: rho, gradrho, omega
  real(kind=DP), intent(out) :: fx, dfx_drho, dfx_dgradrho
!  
#if 1
  real(kind=DP), parameter :: A_HUS =  0.757211D0
  real(kind=DP), parameter :: B_HUS = -0.106304D0
  real(kind=DP), parameter :: C_HUS = -0.118649D0
  real(kind=DP), parameter :: D_HUS =  0.609650D0
#else
  real(kind=DP), parameter :: A_HUS =  0.451606D0
  real(kind=DP), parameter :: B_HUS = -0.371708D0
  real(kind=DP), parameter :: C_HUS = -0.0772155D0
  real(kind=DP), parameter :: D_HUS =  0.577863D0
#endif
!
  real(kind=DP) :: s, kf, ds_drho, ds_dgradrho
  real(kind=DP) :: zeta, dzeta_ds, dzeta_drho, dzeta_dgradrho
  real(kind=DP) ::  eta,  deta_ds,  deta_drho,  deta_dgradrho
  real(kind=DP) :: lambda, dlambda_ds, dlambda_drho, dlambda_dgradrho
  real(kind=DP) :: hs, dhs_ds
  real(kind=DP) :: fs, dfs_ds, dfs_drho, dfs_dgradrho
  real(kind=DP) :: egs, degs_ds, degs_drho, degs_dgradrho

  real(kind=DP) :: nu, nu2, dnu_drho
  real(kind=DP) :: lam2, lam3, lam4
  real(kind=DP) :: chi, chi2, chi3, chi4, chi5
  real(kind=DP) :: dchi_drho, dchi_dgradrho

  real(kind=DP) :: c0a, c0b, c0c, c1, c2, c3, c4a, c4b, c4c, c5a, c5b, c5c
  real(kind=DP) :: dc1, dc2, dc3, dc4a, dc4b, dc4c, dc5a, dc5b, dc5c
  real(kind=DP) :: ddc1, ddc2, ddc3, ddc4a, ddc4b, ddc4c, ddc5a, ddc5b, ddc5c

  real(kind=DP) :: f1, f2, f3, f4, f5, f6
  real(kind=DP) :: df1, df2, df3, df4, df5, df6
  real(kind=DP) :: ddf1, ddf2, ddf3, ddf4, ddf5, ddf6

  integer :: i

#if 0
  if ( abs(rho) > 1.d-05) then
     kf = ( 3.0d0 *PAI *PAI *rho )**(1.0d0/3.0d0)
     s = gradrho /2.0d0 /kf /rho
     ds_drho = -4.0d0 /3.0d0 *s /rho
     ds_dgradrho = 0.5d0 /kf /rho
  else
     s = 0.d0
     ds_drho = 0.0d0
     ds_dgradrho = 0.0d0
  endif
#else
  kf = ( 3.0d0 *PAI *PAI *rho )**(1.0d0/3.0d0)
  s = gradrho /2.0d0 /kf /rho
  ds_drho = -4.0d0 /3.0d0 *s /rho
  ds_dgradrho = 0.5d0 /kf /rho
#endif
  select case (pot_type)
  case (1)
     call diff_h_case_pbe( s, hs, dhs_ds )
  case (2)
     call diff_h_case_pbesol( s, hs, dhs_ds )
  end select

  call diff_zeta_eta_lamb( s, hs, dhs_ds, zeta, dzeta_ds, eta, deta_ds, &
       &                   lambda, dlambda_ds )
  call diff_fs( s, hs, dhs_ds, fs, dfs_ds )
  call diff_egs( zeta, dzeta_ds, eta, deta_ds, lambda, dlambda_ds, fs, dfs_ds, &
       &         egs, degs_ds )

  dzeta_drho = dzeta_ds *ds_drho
  deta_drho  = deta_ds  *ds_drho
  dlambda_drho = dlambda_ds *ds_drho
  dfs_drho = dfs_ds *ds_drho
  degs_drho = degs_ds *ds_drho

  dzeta_dgradrho = dzeta_ds *ds_dgradrho
  deta_dgradrho = deta_ds *ds_dgradrho
  dlambda_dgradrho = dlambda_ds *ds_dgradrho
  dfs_dgradrho = dfs_ds *ds_dgradrho
  degs_dgradrho = degs_ds *ds_dgradrho
  
  nu = omega /kf;
  nu2 = nu *nu
  c0a = zeta +nu2;     c0b = eta  +nu2;      c0c = lambda +nu2

  chi = nu /sqrt(c0c)
  chi2 = chi *chi;   chi3 = chi2 *chi;    chi4 = chi3 *chi;   chi5 = chi4 *chi

  lam2 = lambda **2;   lam3 = lambda *lam2;   lam4 = lambda *lam3

  dnu_drho = -1.0d0 /3.0d0 *nu /rho

  dchi_drho = dnu_drho /sqrt(c0c) &
       &      -nu /2.0d0 /c0c**1.5d0 *( dlambda_drho +2.0d0 *nu *dnu_drho )
  dchi_dgradrho = -nu /2.0d0 /c0c**1.5d0 *dlambda_dgradrho

  c1 = 1.0d0 -chi
  c2 = 1.0d0 -1.5d0 *chi +0.5d0 *chi3
  c3 = 1.0d0 -15.0d0/8.0d0 *chi +5.0d0/4.0d0 *chi3 -3.0d0/8.0d0 *chi5
  c4a = sqrt(c0a);     c5a = log(nu +c4a)
  c4b = sqrt(c0b);     c5b = log(nu +c4b)
  c4c = sqrt(c0c);     c5c = log(nu +c4c)

! -----------------------
  dc1 = -dchi_drho
  dc2 = ( -1.5d0 +1.5d0 *chi2 ) *dchi_drho
  dc3 = ( -15.0d0 /8.0d0 +15.0d0 /4.0d0 *chi2 -15.0d0 /8.0d0 *chi4 ) *dchi_drho
  dc4a = 0.5d0 /c4a *( dzeta_drho    +2.0d0 *nu *dnu_drho )
  dc4b = 0.5d0 /c4b *( deta_drho     +2.0d0 *nu *dnu_drho )
  dc4c = 0.5d0 /c4c *( dlambda_drho  +2.0d0 *nu *dnu_drho )
  dc5a = ( dnu_drho +dc4a ) /( nu +c4a )
  dc5b = ( dnu_drho +dc4b ) /( nu +c4b )
  dc5c = ( dnu_drho +dc4c ) /( nu +c4c )

! ----------------------
  ddc1 = -dchi_dgradrho
  ddc2 = ( -1.5d0 +1.5d0 *chi2 ) *dchi_dgradrho
  ddc3 = ( -15.0d0 /8.0d0 +15.0d0 /4.0d0 *chi2 -15.0d0 /8.0d0 *chi4 ) *dchi_dgradrho
  ddc4a = 0.5d0 /c4a *( dzeta_dgradrho )
  ddc4b = 0.5d0 /c4b *( deta_dgradrho  ) 
  ddc4c = 0.5d0 /c4c *( dlambda_dgradrho )
  ddc5a = ( ddc4a ) /( nu +c4a )
  ddc5b = ( ddc4b ) /( nu +c4b )
  ddc5c = ( ddc4c ) /( nu +c4c )

! -------------------------
  f1 = -4.0d0 /9.0d0 *B_HUS /lambda *c1
  f2 = -4.0d0 /9.0d0 *C_HUS *fs /lam2 *c2
  f3 = -8.0d0 /9.0d0 *egs /lam3 *c3
  f4 = 2.0d0 *nu *( c4a -c4b )
  f5 = 2.0d0 *zeta *( c5a -c5c )
  f6 = -2.0d0 *eta *( c5b -c5c )

  fx = A_HUS +f1 +f2 +f3 +f4 +f5 +f6

! ------------------------------------
  df1 = dc1 /lambda - c1 /lam2 *dlambda_drho
  df1 = -4.0d0 /9.0d0 *B_HUS *df1

  df2 = dfs_drho /lam2 *c2 +fs *(-2.0d0) /lam3 *dlambda_drho *c2 +fs /lam2 *dc2
  df2 = -4.0d0 /9.0d0 *C_HUS *df2

  df3 = degs_drho /lam3 *c3 -3.0d0 *egs /lam4 *c3 *dlambda_drho +egs /lam3 *dc3
  df3 = -8.0d0 /9.0d0 *df3

  df4 = dnu_drho *( c4a -c4b ) +nu *( dc4a -dc4b )
  df4 = df4 *2.0d0

  df5 = dzeta_drho *( c5a -c5c ) +zeta *( dc5a -dc5c )
  df5 = df5 *2.0d0

  df6 = deta_drho *( c5b -c5c ) +eta *( dc5b -dc5c )
  df6 = -df6 *2.0d0

  dfx_drho = df1 +df2 +df3 +df4 +df5 +df6

! ------------------------------------
  ddf1 = ddc1 /lambda - c1 /lam2 *dlambda_dgradrho
  ddf1 = -4.0d0 /9.0d0 *B_HUS *ddf1

  ddf2 = dfs_dgradrho /lam2 *c2 +fs *(-2.0d0) /lam3 *dlambda_dgradrho *c2 +fs /lam2 *ddc2
  ddf2 = -4.0d0 /9.0d0 *C_HUS *ddf2

  ddf3 = degs_dgradrho /lam3 *c3 -3.0d0 *egs /lam4 *c3 *dlambda_dgradrho +egs /lam3 *ddc3
  ddf3 = -8.0d0 /9.0d0 *ddf3

  ddf4 = nu *( ddc4a -ddc4b )
  ddf4 = ddf4 *2.0d0

  ddf5 = dzeta_dgradrho *( c5a -c5c ) +zeta *( ddc5a -ddc5c )
  ddf5 = ddf5 *2.0d0

  ddf6 = deta_dgradrho *( c5b -c5c ) +eta *( ddc5b -ddc5c )
  ddf6 = -ddf6 *2.0d0

  dfx_dgradrho = ddf1 +ddf2 +ddf3 +ddf4 +ddf5 +ddf6

contains

  subroutine diff_egs( zeta, dzeta_ds, eta, deta_ds, lambda, dlambda_ds, fs, dfs_ds, &
       &               egs, degs_ds )
    real(kind=DP), intent(in) :: zeta, dzeta_ds
    real(kind=DP), intent(in) :: eta,  deta_ds
    real(kind=DP), intent(in) :: lambda, dlambda_ds
    real(kind=DP), intent(in) :: fs, dfs_ds
    real(kind=DP), intent(out) :: egs, degs_ds
    
    real(kind=DP) :: c1, c2, c3, c4, c5, c6, f1, df1
    
    c1 = lambda **2
    c2 = c1 *lambda
    c3 = lambda **2.5d0
    c4 = lambda **3.5d0
    c5 = sqrt(zeta)
    c6 = sqrt(eta)

    f1 = -0.4d0 *C_HUS *fs *lambda -4.0d0/15.0d0 *B_HUS *c1 -1.2d0 *A_HUS *c2 &
         & -0.8d0 *sqrt(PAI) *c4 -2.4d0 *c4 *(c5 -c6)
    
    df1 = -0.40d0 *C_HUS *( dfs_ds *lambda +fs *dlambda_ds ) &
         &     -8.0d0 /15.0d0 *B_HUS *lambda *dlambda_ds &
         &     -18.0d0 /5.0d0 *A_HUS *c1 *dlambda_ds &
         &     -14.0d0 /5.0d0 *sqrt(PAI) *c3 *dlambda_ds &
         &     -2.4d0 *( (3.5d0 *c3 *dlambda_ds)*( c5 -c6 ) &
         &              +c4*( dzeta_ds /c5 -deta_ds /c6 )/2.0d0 )
    egs = f1
    degs_ds = df1
  end subroutine diff_egs

  subroutine diff_fs( s, hs, dhs_ds, fs, dfs_ds )
    real(kind=DP), intent(in) :: s, hs, dhs_ds
    real(kind=DP), intent(out) :: fs, dfs_ds
!
    real(kind=DP) :: s0, s0_sq, s2, c1, c2, c3
    
    s0 = 2.0d0;   s0_sq = s0 *s0
    s2 =  s *s

    c1 = 1.0d0 +s2 /s0_sq 
    c2 = 1.0d0 /27.0d0 /C_HUS *s2 /c1
    c3 = 2.0d0 /27.0d0 /C_HUS *s /c1**2 

    fs = 1.0d0 -c2 -0.5d0 /C_HUS *s2 *hs
    dfs_ds = -c3 -0.5d0 /C_HUS *( 2.0d0 *s *hs +s2 *dhs_ds )
!
  end subroutine diff_fs

  subroutine diff_zeta_eta_lamb( s, hs, dhs_ds, zeta, dzeta_ds, eta, deta_ds, &
       &                         lambda, dlambda_ds )
    real(kind=DP), intent(in) :: s, hs, dhs_ds
    real(kind=DP), intent(out) :: zeta, dzeta_ds
    real(kind=DP), intent(out) :: eta,  deta_ds
    real(kind=DP), intent(out) :: lambda, dlambda_ds
    
    real(kind=DP) :: s2
    
    s2 = s *s
    zeta = s2 *hs
    eta = A_HUS +zeta
    lambda = D_HUS +zeta
    
    dzeta_ds = 2.0d0 *s *hs +s2 *dhs_ds
    deta_ds =  dzeta_ds
    dlambda_ds =  dzeta_ds
  end subroutine diff_zeta_eta_lamb

  subroutine diff_h_case_pbe( s, fout, dfout )
    real(kind=DP), intent(in) :: s
    real(kind=DP), intent(out) :: fout, dfout
!
    real(kind=DP), parameter :: a2 =  0.0159941D0
    real(kind=DP), parameter :: a3 =  0.0852995D0
    real(kind=DP), parameter :: a4 = -0.1603680D0
    real(kind=DP), parameter :: a5 =  0.1526450D0
    real(kind=DP), parameter :: a6 = -0.0971263D0
    real(kind=DP), parameter :: a7 =  0.0422061D0
    
    real(kind=DP), parameter :: b1 =  5.33319D0
    real(kind=DP), parameter :: b2 =-12.47800D0
    real(kind=DP), parameter :: b3 = 11.09880D0
    real(kind=DP), parameter :: b4 = -5.11013D0
    real(kind=DP), parameter :: b5 =  1.714680D0
    real(kind=DP), parameter :: b6 = -0.610380D0
    real(kind=DP), parameter :: b7 =  0.307555D0
    real(kind=DP), parameter :: b8 = -0.0770547D0
    real(kind=DP), parameter :: b9 =  0.0334840D0
    
    real(kind=DP) :: s2, s3, s4, s5, s6, s7, s8, s9
    real(kind=DP) :: c1, c2, c3, c4

    s2 =  s *s;     s3 = s2 *s;     s4 = s3 *s;    s5 = s4 *s;
    s6 = s5 *s;     s7 = s6 *s;     s8 = s7 *s;    s9 = s8 *s
    
    c1 = a2*s2 +a3*s3 +a4*s4 +a5*s5 +a6*s6 +a7*s7
    c2 = 1.0d0 +b1 *s &
         & +b2*s2 +b3*s3 +b4*s4 +b5*s5 +b6*s6 +b7*s7 +b8*s8 +b9 *s9
    fout = c1 /c2
!
    c3 = 2.d0 *a2 *s +3.0d0 *a3 *s2 +4.0d0 *a4 *s3 +5.0d0 *a5 *s4 +6.0d0 *a6 *s5 &
         & +7.0d0 *a7 *s6
    c4 = b1 +2.d0 *b2 *s +3.0d0 *b3 *s2 +4.0d0 *b4 *s3 +5.0d0 *b5 *s4 +6.0d0 *b6 *s5 &
         & +7.0d0 *b7 *s6 +8.0d0 *b8 *s7 +9.0d0 *b9 *s8
    
    dfout = c3 /c2 -c1 /c2**2 *c4
  end subroutine diff_h_case_pbe

  subroutine diff_h_case_pbesol( s, fout, dfout )
    real(kind=DP), intent(in) :: s
    real(kind=DP), intent(out) :: fout, dfout
!
    real(kind=DP), parameter :: a2 =  0.0047333D0
    real(kind=DP), parameter :: a3 =  0.0403304D0
    real(kind=DP), parameter :: a4 = -0.0574615D0
    real(kind=DP), parameter :: a5 =  0.0435395D0
    real(kind=DP), parameter :: a6 = -0.0216251D0
    real(kind=DP), parameter :: a7 =  0.0063721D0
    
    real(kind=DP), parameter :: b1 =  8.52056D0
    real(kind=DP), parameter :: b2 =-13.9885D0
    real(kind=DP), parameter :: b3 =  9.28583D0
    real(kind=DP), parameter :: b4 = -3.27287D0
    real(kind=DP), parameter :: b5 =  0.843499D0
    real(kind=DP), parameter :: b6 = -0.235543D0
    real(kind=DP), parameter :: b7 =  0.0847074D0
    real(kind=DP), parameter :: b8 = -0.0171561D0
    real(kind=DP), parameter :: b9 =  0.0050552D0
    
    real(kind=DP) :: s2, s3, s4, s5, s6, s7, s8, s9
    real(kind=DP) :: c1, c2, c3, c4

    s2 =  s *s;     s3 = s2 *s;     s4 = s3 *s;    s5 = s4 *s;
    s6 = s5 *s;     s7 = s6 *s;     s8 = s7 *s;    s9 = s8 *s
    
    c1 = a2*s2 +a3*s3 +a4*s4 +a5*s5 +a6*s6 +a7*s7
    c2 = 1.0d0 +b1 *s &
         & +b2*s2 +b3*s3 +b4*s4 +b5*s5 +b6*s6 +b7*s7 +b8*s8 +b9 *s9
    fout = c1 /c2
!
    c3 = 2.d0 *a2 *s +3.0d0 *a3 *s2 +4.0d0 *a4 *s3 +5.0d0 *a5 *s4 +6.0d0 *a6 *s5 &
         & +7.0d0 *a7 *s6
    c4 = b1 +2.d0 *b2 *s +3.0d0 *b3 *s2 +4.0d0 *b4 *s3 +5.0d0 *b5 *s4 +6.0d0 *b6 *s5 &
         & +7.0d0 *b7 *s6 +8.0d0 *b8 *s7 +9.0d0 *b9 *s8
    
    dfout = c3 /c2 -c1 /c2**2 *c4
  end subroutine diff_h_case_pbesol
  
end subroutine omegapbe_analytic

! ---------------------------------------
! Gau-PBE
!
! J. Chem. Phys. 135 (2011) 071103.
! ---------------------------------------
subroutine ex_gaupbe( nspin, ispin, ista_r, iend_r, chgrhr_l, grad_rho, &
     &                wos, exc, dFx_drho, dFx_dgradrho )
  use m_Const_Parameters,  only : DP,PAI
  implicit none

  integer,intent(in)        :: nspin,ispin,ista_r,iend_r
  real(kind=DP),intent(in)  :: chgrhr_l(ista_r:iend_r,ispin)
  real(kind=DP),intent(in)  :: grad_rho(ista_r:iend_r,nspin)
  real(kind=DP),intent(in)  :: wos(ista_r:iend_r)
  real(kind=DP),intent(out) :: exc
  real(kind=DP),intent(out) :: dFx_drho(ista_r:iend_r,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(ista_r:iend_r,nspin)

  real(kind=DP), parameter :: ax = -0.7385587663820224d0 !! ax = -(3/4)*(3/PI)^(1/3)
  real(kind=DP), parameter :: thrd = 0.33333333333d0, thrd4 = 1.333333333333333d0
  real(kind=DP), parameter :: thpith = 3.0936677262801d0

! ----- parameters of enhanced factor --
  real(kind=DP), parameter :: mu_pbe    = 0.2195149727645171d0
  real(kind=DP), parameter :: kappa_pbe    = 0.804d0

  real(kind=DP), parameter :: alpha_gaupbe = 0.15d0
  real(kind=DP), parameter :: beta_gaupbe = 0.24d0

  integer       :: is, i

  real(kind=DP) :: facw, d, dd, fk, s, fac, s2, mu_s2, f, ex, df_ds, &
       &           exd, exdd, exc0, excd, excdd, exc1
  real(kind=DP) :: x, dx_ds, ctmp1, df_dx
  real(kind=DP) :: df_drho, df_ddrho
  real(kind=DP) :: vald, dvald_drho, dvald_ddrho
  real(kind=DP) :: valg, dvalg_drho, dvalg_ddrho

!---- Spin dependency

  facw = ispin
#if 0
  facw = 1.0d0
  if ( nspin == 1 ) facw = 0.5d0
#endif

  exc  = 0.d0
#ifdef NEC_TUNE_MXCP
!CDIR NOCONCUR
#endif
  do is = 1, ispin
     exc1 = 0.d0
#ifdef NEC_TUNE_MXCP
!CDIR INNER
#endif
     do i = ista_r, iend_r
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

        x = mu_pbe *s2;   dx_ds = 2.0d0 *mu_pbe *s
        ctmp1 = kappa_pbe + x
        
        f = 1.0d0 +kappa_pbe *x /ctmp1
        df_dx = kappa_pbe**2 /ctmp1**2
        df_ds = dx_ds *df_dx
!
        if(d > 1.d-05) then
           df_drho = -4.0d0 /3.0d0 *s /d *df_ds
           df_ddrho = 1.0d0 /(d*fk*2.d0) *df_ds
        else
           df_drho = 0.0d0
           df_ddrho = 0.0d0
        endif
! 
        call diff_d( d, dd, s, f, df_drho, df_ddrho, vald, dvald_drho, dvald_ddrho )
        call diff_g( vald, dvald_drho, dvald_ddrho, valg, dvalg_drho, dvalg_ddrho )
        
        ex  = fac *f *d *valg
        exd = thrd4 *fac* ( f -s *df_ds ) *valg &
             &   +fac *d *f *dvalg_drho
        exdd = ax * df_ds *0.5d0 /thpith *valg &
             &   +fac *d *f *dvalg_ddrho
!------------------------------------------

#if 1
        exc0 = ex / facw
#endif
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
     end do
     exc = exc + exc1
  end do

#if 0
  if ( nspin == 2 ) then
     exc = exc /2.0d0
  else
     dFx_drho = dFx_drho *2.0d0
     dFx_dgradrho = dFx_dgradrho *2.0d0
  endif
#endif

contains

  subroutine diff_g( d, dd_dn, dd_ddn, fout, dfout, ddfout )
    real(kind=DP), intent(in) :: d, dd_dn, dd_ddn
    real(kind=DP), intent(out) :: fout, dfout, ddfout
    !
    real(kind=DP) :: c1, c2, c3, c4, c5, c6, c7, c8
    real(kind=DP) :: e, de_dn, de_ddn
    !
    call diff_e( d, dd_dn, dd_ddn, e, de_dn, de_ddn )
    
    c1 = 4.0d0 /3.0d0 *beta_gaupbe *sqrt(PAI) /sqrt(alpha_gaupbe)
    
    c2 = sqrt(PAI) *erf( 0.5d0 /d )
    c3 = 2.0d0 *d -16.0d0 *d**3
    c4 = c2 +c3 *e -4.0d0 *d
    
    fout = 1.0d0 -c1 *d *c4
    !
    c5 = -exp( -(0.5d0 /d)**2 ) /d**2
    c6 = 2.0d0 -48.0d0 *d**2
    
    c7 = ( c5 +c6 *e -4.0d0 ) *dd_dn  +c3 *de_dn
    c8 = ( c5 +c6 *e -4.0d0 ) *dd_ddn +c3 *de_ddn
    
    dfout = -c1 *( c4 *dd_dn + d *c7 )
    ddfout = -c1 *( c4 *dd_ddn +d *c8 )
    
  end subroutine diff_g
  
  subroutine diff_e( d, dd_dn, dd_ddn, fout, dfout, ddfout )
    implicit none
    real(kind=DP), intent(in) :: d, dd_dn, dd_ddn
    real(kind=DP), intent(out) :: fout, dfout, ddfout
    !
    real(kind=DP) :: c1, c2, c3
    
    c1 = -0.25d0 /d**2
    c2 = exp( c1 )
    
    fout = c2 -1.0d0
    
    c3 = 0.50d0 /d**3
    dfout = c3 *c2 *dd_dn
    
    ddfout = c3 *c2 *dd_ddn
  end subroutine diff_e
  
  subroutine diff_d( n, dn, s, fs, dfs_dn, dfs_ddn, fout, dfout, ddfout )
    implicit none
    real(kind=DP), intent(in) :: n, dn, s, fs, dfs_dn, dfs_ddn
    real(kind=DP), intent(out) :: fout, dfout, ddfout
    !
    real(kind=DP) :: c1, c2, c3
    
    c1 = sqrt( alpha_gaupbe ) /6.0d0 /sqrt(PAI) *sqrt( -ax )
    
    fout = c1 *n**(-1.0d0/3.0d0) *sqrt(fs)
    
    c2 = -1.0d0 /3.d0 /n +dfs_dn /fs /2.0d0
    dfout = fout *c2
    
    c3 = dfs_ddn /2.0d0 /fs
    ddfout = fout *c3
    
  end subroutine diff_d

end subroutine ex_gaupbe

