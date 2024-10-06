!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUTINE: sphrp2_for_Berry
!
!  AUTHOR(S): T. Yamamoto   May/01/2004
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
  subroutine sphrp2_for_Berry(is,gv,ylm)
! $ID: $
    use m_Const_Parameters,    only : DP, PAI, PAI4
    integer, intent(in) :: is
    real(kind=DP), intent(in), dimension(3)  :: gv
    real(kind=DP), intent(out) :: ylm
!!$C spherical harmonics.
!!$c     Some lines are rewitted by T.Yamasaki according to Y.Morikawa's
!!$c   e-mail. at 4th May 1994

    real(kind=DP) :: a,b,c,d,e,f,g
    real(kind=DP) :: gx,gy,gz,gr
    real(kind=DP), parameter :: g_minimum = 1.d-10

    gx = gv(1); gy = gv(2); gz = gv(3); gr = dsqrt(gx**2+gy**2+gz**2)
    
    if(gr < g_minimum) then
      ylm = dsqrt(1.d0/PAI4)
      return
    end if

    if(is == 1) then
       a = dsqrt(1.d0/PAI4)
       ylm = a
    else if(is == 2) then
       a = dsqrt(3.d0/PAI4)
       ylm = a*gx/gr
    else if(is == 3) then
       a = dsqrt(3.d0/PAI4)
       ylm = a*gy/gr
    else if(is == 4) then
       a = dsqrt(3.d0/PAI4)
       ylm = a*gz/gr
    else if(is == 5) then
       a = dsqrt(5.d0/(16*PAI))
       d = gz**2
       e = gx**2+gy**2+d
       ylm = a*(3*d-e)/e
    else if(is == 6) then
       a = dsqrt(15.d0/(16*PAI))
       b = gx**2
       c = gy**2
       e = b+c+gz**2
       ylm = a*(b-c)/e
    else if(is == 7) then
       a = dsqrt(15.d0/PAI4)
       e = gr**2
       ylm = a*gx*gy/e
    else if(is == 8) then
       a = dsqrt(15.d0/PAI4)
       e = gr**2
       ylm = a*gy*gz/e
    else if(is == 9) then
       a = dsqrt(15.d0/PAI4)
       e = gr**2
       ylm = a*gz*gx/e
    else if(is == 10) then
       a = dsqrt(7.d0/(16*PAI))
       d = gz**2
       e = gr**2
       f = e * gr
       ylm = a*gz*(5*d-3*e)/f
    else if(is == 11) then
       a = dsqrt(21.d0/(32*PAI))
       d = gz**2
       e = gr**2
       f = e * gr
       ylm = a*gx*(5*d-e)/f
    else if(is == 12) then
       a = dsqrt(21.d0/(32*PAI))
       d = gz**2
       e = gr**2
       f = e * gr
       ylm = a*gy*(5*d-e)/f
    else if(is == 13) then
       a = dsqrt(105.d0/(16*PAI))
       b = gx**2
       c = gy**2
       e = gr**2
       f = e * gr
       ylm = a*gz*(b-c)/f
    else if(is == 14) then
       a = dsqrt(105.d0/PAI4)
       f = gr**3
       ylm = a*gx*gy*gz/f
    else if(is == 15) then
       a = dsqrt(35.d0/(32*PAI))
       b = gx**2
       c = gy**2
       f = gr**3
       ylm = a*gx*(b-3*c)/f
    else if(is == 16) then
       a = dsqrt(35.d0/(32*pai))
       b = gx**2
       c = gy**2
       f = gr**3
       ylm = a*gy*(3*b-c)/f
    else if(is == 17) then
       a = 3.d0/8.d0/dsqrt(PAI4)
       d = gz**2
       e = gr**2
       f = e**2
       ylm = a*(5*d*(7*d-6*e)/f+3.d0)
    else if(is == 18) then
       a = 15.d0/4.d0/dsqrt(10*PAI)
       d = gz**2
       e = gr**2
       f = e**2
       ylm = a*gz*gx*(7*d-3*e)/f
    else if(is == 19) then
       a = 15.d0/4.d0/dsqrt(10.d0*PAI)
       d = gz**2
       e = gr**2
       f = e**2
       ylm = a*gy*gz*(7*d-3*e)/f
    else if(is == 20) then
       a = 15.d0/8.d0/dsqrt(5*PAI)
       b = gx**2
       c = gy**2
       d = gz**2
       e = b+c+d
       f = e**2
       ylm = a*(7*d-e)*(b-c)/f
    else if(is == 21) then
       a = 15.d0/4.d0/dsqrt(5*PAI)
       d = gz**2
       e = gr**2
       f = e**2
       ylm = a*(7*d-e)*gx*gy/f
    else if(is == 22) then
       a = 105.d0/4.d0/dsqrt(70*pai)
       b = gx**2
       c = gy**2
       f = (b+c+gz**2)**2
       ylm = a*(b-3*c)*gz*gx/f
    else if(is == 23) then
       a = 105.d0/4.d0/dsqrt(70*PAI)
       b = gx**2
       c = gy**2
       f = (b+c+gz**2)**2
       ylm = a*(3.d0*b-c)*gy*gz/f
    else if(is == 24) then
       a = 105.d0/16.d0/dsqrt(35*PAI)
       b = gx**2
       c = gy**2
       f = (b+c+gz**2)**2
       ylm = a*((b-c)**2-4*b*c)/f
    else if(is == 25) then
       a = 105.d0/4.d0/dsqrt(35*PAI)
       b = gx**2
       c = gy**2
       f = (b+c+gz**2)**2
       ylm = a*(b-c)*gx*gy/f
    else if(is == 26) then
       a = dsqrt(11.d0/PAI)/16.d0
       b = gz*gz
       c = gx*gx+gy*gy+b
       f = gr**5
       ylm = a*gz*(63.d0*b*b-70.d0*b*c+15.d0*c*c)/f
    else if(is == 27) then
       a = dsqrt(165.d0/PAI)/16.d0
       b = gz*gz
       c = gx*gx+gy*gy+b
       f = gr**5
       ylm = a*gx*(21.d0*b*b-14.d0*b*c+c*c)/f
    else if(is == 28) then
       a = dsqrt(165.d0/PAI)/16.d0
       b = gz*gz
       c = gx*gx+gy*gy+b
       f = gr**5
       ylm = a*gy*(21.d0*b*b-14.d0*b*c+c*c)/f
    else if(is == 29) then
       a = dsqrt(1155.d0/(64*PAI))
       b = gz*gz
       c = gx*gx
       d = gy*gy
       e = c+d
       f = gr**5
       ylm = a*gz*(2.d0*b-e)*(c-d)/f
    else if(is == 30) then
       a = dsqrt(1155.d0/(16.d0*PAI))
       b = gz*gz
       c = gx*gx
       d = gy*gy
       e = c+d
       f = gr**5
       ylm = a*gx*gy*gz*(2.d0*b-e)/f
    else if(is == 31) then
       a = dsqrt(385.d0/(512.d0*PAI))
       b = gz*gz
       c = gx*gx
       d = gy*gy
       e = c+d
       f = gr**5
       ylm = a*gx*(8.d0*b-e)*(c-3.d0*d)/f
    else if(is == 32) then
       a = dsqrt(385.d0/(512.d0*PAI))
       b = gz*gz
       c = gx*gx
       d = gy*gy
       e = c+d
       f = gr**5
       ylm = a*gy*(8.d0*b-e)*(3.d0*c-d)/f
    else if(is == 33) then
       a = dsqrt(385.d0/(256.d0*PAI))*3.d0
       b = gx*gx
       c = gy*gy
       f = gr**5
       ylm = a*gz*(b*b-6.d0*b*c+c*c)/f
    else if(is == 34) then
       a = dsqrt(385.d0/(16.d0*PAI))*3.d0
       b = gx*gx-gy*gy
       f = gr**5
       ylm = a*gx*gy*gz*b/f
    else if(is == 35) then
       a = dsqrt(693.d0/(512.d0*PAI))
       b = gx*gx
       c = gy*gy
       f = gr**5
       ylm = a*gx*(b*b-10.d0*b*c+5.d0*c*c)/f
    else if(is == 36) then
       a = dsqrt(693.d0/(512.d0*PAI))
       b = gy*gy
       c = gx*gx
       f = gr**5
       ylm = a*gy*(b*b-10.d0*b*c+5.d0*c*c)/f
    else if(is == 37) then
       a = dsqrt(13.d0/1024.d0*PAI)
       b = gz*gz
       c = b*b
       d = c*b
       e = gx*gx+gy*gy+b
       c = c*e
       f = e*e
       b = b*f
       f = f*e
       ylm = a*((231.d0*d-315.d0*c+105.d0*b)/f-5.d0)
    else if(is == 38) then
       a = dsqrt(273.d0/(256.d0*PAI))
       b = gz*gz
       c = b*b
       e = gx*gx+gy*gy+b
       b = b*e
       d = e*e
       f = d*e
       ylm = a*gz*gx*(33.d0*c-30.d0*b+5.d0*d)/f
    else if(is == 39) then
       a = dsqrt(273.d0/(256.d0*PAI))
       b = gz*gz
       c = b*b
       e = gx*gx+gy*gy+b
       b = b*e
       d = e*e
       f = d*e
       ylm = a*gy*gz*(33.d0*c-30.d0*b+5.d0*d)/f
    else if(is == 40) then
       a = dsqrt(1365.d0/(2048.d0*PAI))
       b = gz*gz
       c = b*b
       d = gx*gx
       f = gy*gy
       e = d+f+b 
       d = d-f
       b = b*e
       g = e*e
       f = g*e
       ylm = a*d*(33.d0*c-18.d0*b+g)/f
    else if(is == 41) then
       a = dsqrt(1365.d0/(512.d0*PAI))
       b = gz*gz
       c = b*b
       e = gx*gx+gy*gy+b 
       b = b*e
       g = e*e
       f = g*e
       ylm = a*gx*gy*(33.d0*c-18.d0*b+g)/f
    else if(is == 42) then
       a = dsqrt(1365.d0/(512.d0*PAI))
       b = gz*gz
       c = gx*gx
       d = gy*gy
       e = c+d+b 
       f = e**3
       ylm = a*gz*gx*(11.d0*b-3.d0*e)*(c-3.d0*d)/f
    else if(is == 43) then
       a = dsqrt(1365.d0/(512.d0*PAI))
       b = gz*gz
       c = gy*gy
       d = gx*gx
       e = d+c+b 
       f = e**3
       ylm = a*gy*gz*(11.d0*b-3.d0*e)*(3.d0*d-c)/f
    else if(is == 44) then
       a = dsqrt(819.d0/(1024d0*PAI))
       b = gx*gx
       c = gy*gy
       d = gz*gz
       e = b+c+d
       f = e**3
       ylm = a*(11.d0*d-e)*(b*b-6.d0*b*c+c*c)/f
    else if(is == 45) then
       a = dsqrt(819.d0/(64.d0*PAI))
       b = gx*gx
       c = gy*gy
       d = gz*gz
       e = b+c+d
       f = e**3
       ylm = a*gx*gy*(11.d0*d-e)*(b-c)/f
    else if(is == 46) then
       a = dsqrt(9009.d0/(512.d0*PAI))
       b = gx*gx
       c = gy*gy
       d = gz*gz
       e = b+c+d
       f = e**3
       ylm = a*gz*gx*(b*b-10.d0*b*c+5.d0*c*c)/f
    else if(is == 47) then
       a = dsqrt(9009.d0/(512.d0*PAI))
       b = gy*gy
       c = gx*gx
       e = c+b+gz*gz
       f = e**3
       ylm = a*gy*gz*(b*b-10.d0*b*c+5.d0*c*c)/f
    else if(is == 48) then
       a = dsqrt(3003.d0/(2048.d0*PAI))
       b = gx*gx
       c = gy*gy
       e = b+c+gz*gz
       f = e**3
       ylm = a*(b-c)*(b*b-14.d0*b*c+c*c)/f
    else if(is == 49) then
       a = dsqrt(3003.d0/(2048.d0*PAI))
       b = gx*gx
       c = gy*gy
       d = b+c
       e = b+c+gz*gz
       f = e**3
       ylm = a*gx*gy*(6.d0*d*d-32.d0*b*c)/f
    end if
  end subroutine sphrp2_for_Berry
