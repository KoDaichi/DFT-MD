#define GGA_PAW_REVISED
#define USE_3RD_ORDER_DERIVATIVE

! from b_XC_Potential.F90

subroutine ex_ggapbe_paw_drv2(nrc,dnr,nspin,chgrhr_l,grad_rho,exc, &
                                dFx_drho,dFx_dgradrho, &
                                dFx_drr,dFx_drg,dFx_dgg, &
                                dFx_drrr,dFx_drrg,dFx_drgg,dFx_dggg)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,PAI
  implicit none

  integer,intent(in)        :: nrc,dnr,nspin
  real(kind=DP),intent(in)  :: chgrhr_l(nrc,nspin)
  real(kind=DP),intent(in)  :: grad_rho(nrc,nspin)
!  real(kind=DP),intent(in)  :: wos(nrc)  
  real(kind=DP),intent(out) :: exc(nrc)
  real(kind=DP),intent(out) :: dFx_drho(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drr(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dgg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drrr(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drrg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drgg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dggg(nrc,nspin)

#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter ::  ax = -0.7385587663820224d0
  real(kind=DP), parameter :: thrd = 0.33333333333333333333333d0, &
       &                      thrd4 = 1.3333333333333333333333d0
  real(kind=DP), parameter :: thpith = 3.0936677262801d0      ! ???
  real(kind=DP), parameter :: yum = 0.2195149727645171d0
  real(kind=DP) :: cupa

  real(kind=DP), parameter :: ninth2=0.2222222222222222222d0
  real(kind=DP), parameter :: ninth4=0.4444444444444444444d0
  real(kind=DP), parameter :: ninth8=0.8888888888888888888d0
  real(kind=DP), parameter :: thrd2 =0.6666666666666666666d0
  real(kind=DP), parameter :: eghth =0.125d0
  real(kind=DP), parameter :: twtysvnth8 = 0.29629629629629627985d0
  real(kind=DP), parameter :: twtysvnth32= 1.18518518518518511939d0

#else

  real(kind=DP), parameter :: ax = -0.7385588d0
  real(kind=DP), parameter :: thrd = 0.33333333333d0, thrd4 = 1.333333333333333d0
  real(kind=DP), parameter :: thpith = 3.0936677262801d0
  real(kind=DP), parameter :: cupa = 0.804d0, yum = 0.2195149727645171d0
  
  real(kind=DP), parameter :: ninth2=0.22222222222d0
  real(kind=DP), parameter :: ninth4=0.44444444444d0
  real(kind=DP), parameter :: ninth8=0.88888888888d0
  real(kind=DP), parameter :: thrd2 =0.66666666666d0
  real(kind=DP), parameter :: eghth =0.125d0
  real(kind=DP), parameter :: twtysvnth8 = 0.296296296296d0
  real(kind=DP), parameter :: twtysvnth32= 1.18518518518519d0
#endif

  real(kind=DP) :: facw,d,dd,fk,s,fac,s2,yums2,f,ex,fs,exd,exdd,exc0,excd,excdd,exc1
  integer       :: is,i
  real(kind=DP) :: fss,fsss,ivd,ivdthrd2

#ifdef GGA_PAW_REVISED
  cupa = 0.8040d0
#endif

!---- Spin dependency

  facw = nspin
  exc  = 0.d0
  dFx_drho = 0.d0
  dFx_dgradrho = 0.d0
  dFx_drr = 0.0d0
  dFx_drg = 0.0d0
  dFx_dgg = 0.0d0
  dFx_drrr = 0.0d0
  dFx_drrg = 0.0d0
  dFx_drgg = 0.0d0
  dFx_dggg = 0.0d0

  do is = 1, nspin
     exc1 = 0.d0
     do i = 1,nrc,dnr
        d  = facw * chgrhr_l(i, is)
        if (d.le.0.0d0) cycle
        dd = facw * grad_rho(i, is)
        if(d > 1.d-05) then
           fk = (3.d0*PAI*PAI*d)**thrd
!!$        sk = dsqrt(4*fk/PAI)
           s = dd/(d*fk*2.d0)
        else
           fk = 0.d0
           s  = 0.d0
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
!           dFx_dgradrho(i, is) = excdd / grad_rho(i, is)           ! 2010 5
           dFx_dgradrho(i, is) = excdd
        else
           dFx_dgradrho(i, is) = 0.d0
        endif
! gradient of charge density 95/12/3 H.S.
!        exc1 = exc1 + exc0*wos(i)
        exc(i)=exc(i)+exc0
        
! *****  Higher derivatives *****
       
        fss = 2.d0 * cupa ** 2 * yum * ( cupa - 3.d0*yums2 ) / ( cupa + yums2 ) ** 3
        fsss= 24.d0 * cupa **2 * yum **2 * s * ( yums2 - cupa ) /  ( cupa + yums2 ) ** 4
        
        ivd=1.d0/d
        ivdthrd2=ivd**thrd2
        
        if(nspin == 1) then
            dFx_drr(i,is)   =ninth4*ax*ivdthrd2*( f-s*fs+4.d0*fss*s2 )
            dFx_drg(i,is)   =-thrd2*ax/thpith*fss*s*ivd
            dFx_dgg(i,is)   =0.25d0*ax/thpith/thpith*fss*ivdthrd2**2
#ifdef USE_3RD_ORDER_DERIVATIVE
            dFx_drrr(i,is)  =-twtysvnth8*ax*ivdthrd2*ivd*( f-fs*s+18.d0*fss*s2+8.d0*fsss*s2*s )
            dFx_drrg(i,is)  =ninth2*ax*ivd*ivd/thpith*( 7.d0*fss*s+4.d0*fsss*s2 )
            dFx_drgg(i,is)  =-thrd*ax/thpith/thpith*ivd*ivdthrd2**2*( fss+fsss*s )
            dFx_dggg(i,is)  =eghth*ax/(thpith**3)*fsss*ivd**2*ivdthrd2
#endif
        else if(nspin == 2) then
            dFx_drr(i,is)   =ninth8*ax*ivdthrd2*( f-s*fs+4.d0*fss*s2 )
            dFx_drg(i,is)   =-thrd4*ax/thpith*fss*s*ivd
            dFx_dgg(i,is)   =0.5d0*ax/thpith/thpith*fss*ivdthrd2**2
#ifdef USE_3RD_ORDER_DERIVATIVE
            dFx_drrr(i,is)  =-twtysvnth32*ax*ivdthrd2*ivd*( f-fs*s+18.d0*fss*s2+8.d0*fsss*s2*s )
            dFx_drrg(i,is)  =ninth8*ax*ivd*ivd/thpith*( 7.d0*fss*s+4.d0*fsss*s2 )
            dFx_drgg(i,is)  =-thrd4*ax/(thpith**2)*ivd*ivdthrd2**2*( fss+fsss*s )
            dFx_dggg(i,is)  =0.5d0*ax/(thpith**3)*fsss*ivd**2*ivdthrd2
#endif
        end if
        
     end do
!     exc(i) = exc(i) + exc1
  end do
end subroutine ex_ggapbe_paw_drv2
! === For nrc decomposion. by takto 2012/12/05 =================================
subroutine ex_ggapbe_paw_drv2_3D(nrc,dnr,nspin,chgrhr_l,grad_rho,exc, &
                                dFx_drho,dFx_dgradrho, &
                                dFx_drr,dFx_drg,dFx_dgg, &
                                dFx_drrr,dFx_drrg,dFx_drgg,dFx_dggg,ista_nrc,iend_nrc,ist,ien)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,PAI
  implicit none

  integer,intent(in)        :: nrc,dnr,nspin
  real(kind=DP),intent(in)  :: chgrhr_l(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(in)  :: grad_rho(ista_nrc:iend_nrc,nspin)
!  real(kind=DP),intent(in)  :: wos(nrc)  
  real(kind=DP),intent(out) :: exc(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFx_drho(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drr(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drg(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dgg(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drrr(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drrg(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drgg(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dggg(ista_nrc:iend_nrc,nspin)
  integer :: ista_nrc, iend_nrc, ist, ien

#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter ::  ax = -0.7385587663820224d0
  real(kind=DP), parameter :: thrd = 0.33333333333333333333333d0, &
       &                      thrd4 = 1.3333333333333333333333d0
  real(kind=DP), parameter :: thpith = 3.0936677262801d0      ! ???
  real(kind=DP), parameter :: yum = 0.2195149727645171d0
  real(kind=DP) :: cupa

  real(kind=DP), parameter :: ninth2=0.2222222222222222222d0
  real(kind=DP), parameter :: ninth4=0.4444444444444444444d0
  real(kind=DP), parameter :: ninth8=0.8888888888888888888d0
  real(kind=DP), parameter :: thrd2 =0.6666666666666666666d0
  real(kind=DP), parameter :: eghth =0.125d0
  real(kind=DP), parameter :: twtysvnth8 = 0.29629629629629627985d0
  real(kind=DP), parameter :: twtysvnth32= 1.18518518518518511939d0

#else

  real(kind=DP), parameter :: ax = -0.7385588d0
  real(kind=DP), parameter :: thrd = 0.33333333333d0, thrd4 = 1.333333333333333d0
  real(kind=DP), parameter :: thpith = 3.0936677262801d0
  real(kind=DP), parameter :: cupa = 0.804d0, yum = 0.2195149727645171d0
  
  real(kind=DP), parameter :: ninth2=0.22222222222d0
  real(kind=DP), parameter :: ninth4=0.44444444444d0
  real(kind=DP), parameter :: ninth8=0.88888888888d0
  real(kind=DP), parameter :: thrd2 =0.66666666666d0
  real(kind=DP), parameter :: eghth =0.125d0
  real(kind=DP), parameter :: twtysvnth8 = 0.296296296296d0
  real(kind=DP), parameter :: twtysvnth32= 1.18518518518519d0
#endif

  real(kind=DP) :: facw,d,dd,fk,s,fac,s2,yums2,f,ex,fs,exd,exdd,exc0,excd,excdd,exc1
  integer       :: is,i
  real(kind=DP) :: fss,fsss,ivd,ivdthrd2

#ifdef GGA_PAW_REVISED
  cupa = 0.8040d0
#endif

!---- Spin dependency

  facw = nspin
  exc  = 0.d0
  dFx_drho = 0.d0
  dFx_dgradrho = 0.d0
  dFx_drr = 0.0d0
  dFx_drg = 0.0d0
  dFx_dgg = 0.0d0
  dFx_drrr = 0.0d0
  dFx_drrg = 0.0d0
  dFx_drgg = 0.0d0
  dFx_dggg = 0.0d0

  do is = 1, nspin
     exc1 = 0.d0
     do i = ist,ien,dnr
        d  = facw * chgrhr_l(i, is)
        if (d.le.0.0d0) cycle
        dd = facw * grad_rho(i, is)
        if(d > 1.d-05) then
           fk = (3*PAI*PAI*d)**thrd
!!$        sk = dsqrt(4*fk/PAI)
           s = dd/(d*fk*2)
        else
           fk = 0.d0
           s  = 0.d0
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
!           dFx_dgradrho(i, is) = excdd / grad_rho(i, is)           ! 2010 5
           dFx_dgradrho(i, is) = excdd
        else
           dFx_dgradrho(i, is) = 0.d0
        endif
! gradient of charge density 95/12/3 H.S.
!        exc1 = exc1 + exc0*wos(i)
        exc(i)=exc(i)+exc0
        
! *****  Higher derivatives *****
       
        fss = 2 * cupa ** 2 * yum * ( cupa - 3*yums2 ) / ( cupa + yums2 ) ** 3
        fsss= 24 * cupa **2 * yum **2 * s * ( yums2 - cupa ) /  ( cupa + yums2 ) ** 4
        
        ivd=1.d0/d
        ivdthrd2=ivd**thrd2
        
        if(nspin == 1) then
            dFx_drr(i,is)   =ninth4*ax*ivdthrd2*( f-s*fs+4.d0*fss*s2 )
            dFx_drg(i,is)   =-thrd2*ax/thpith*fss*s*ivd
            dFx_dgg(i,is)   =0.25d0*ax/thpith/thpith*fss*ivdthrd2**2
#ifdef USE_3RD_ORDER_DERIVATIVE
            dFx_drrr(i,is)  =-twtysvnth8*ax*ivdthrd2*ivd*( f-fs*s+18.d0*fss*s2+8.d0*fsss*s2*s )
            dFx_drrg(i,is)  =ninth2*ax*ivd*ivd/thpith*( 7.d0*fss*s+4.d0*fsss*s2 )
            dFx_drgg(i,is)  =-thrd*ax/thpith/thpith*ivd*ivdthrd2**2*( fss+fsss*s )
            dFx_dggg(i,is)  =eghth*ax/(thpith**3)*fsss*ivd**2*ivdthrd2
#endif
        else if(nspin == 2) then
            dFx_drr(i,is)   =ninth8*ax*ivdthrd2*( f-s*fs+4.d0*fss*s2 )
            dFx_drg(i,is)   =-thrd4*ax/thpith*fss*s*ivd
            dFx_dgg(i,is)   =0.5d0*ax/thpith/thpith*fss*ivdthrd2**2
#ifdef USE_3RD_ORDER_DERIVATIVE
            dFx_drrr(i,is)  =-twtysvnth32*ax*ivdthrd2*ivd*( f-fs*s+18.d0*fss*s2+8.d0*fsss*s2*s )
            dFx_drrg(i,is)  =ninth8*ax*ivd*ivd/thpith*( 7.d0*fss*s+4.d0*fsss*s2 )
            dFx_drgg(i,is)  =-thrd4*ax/(thpith**2)*ivd*ivdthrd2**2*( fss+fsss*s )
            dFx_dggg(i,is)  =0.5d0*ax/(thpith**3)*fsss*ivd**2*ivdthrd2
#endif
        end if
        
     end do
!     exc(i) = exc(i) + exc1
  end do
end subroutine ex_ggapbe_paw_drv2_3D
! ==============================================================================

subroutine ex_ggapbe_paw_drv3(nrc,num_ir,irs,nspin,chgrhr_l,grad_rho,exc, &
                                dFx_drho,dFx_dgradrho, &
                                dFx_drr,dFx_drg,dFx_dgg, &
                                dFx_drrr,dFx_drrg,dFx_drgg,dFx_dggg)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,PAI
  implicit none

  integer,intent(in)        :: nrc,num_ir,irs(num_ir),nspin
  real(kind=DP),intent(in)  :: chgrhr_l(nrc,nspin)
  real(kind=DP),intent(in)  :: grad_rho(nrc,nspin)
!  real(kind=DP),intent(in)  :: wos(nrc)  
  real(kind=DP),intent(out) :: exc(nrc)
  real(kind=DP),intent(out) :: dFx_drho(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drr(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dgg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drrr(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drrg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drgg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dggg(nrc,nspin)

  real(kind=DP), parameter :: ax = -0.7385588d0
  real(kind=DP), parameter :: thrd = 0.33333333333d0, thrd4 = 1.333333333333333d0
  real(kind=DP), parameter :: thpith = 3.0936677262801d0
  real(kind=DP), parameter :: cupa = 0.804d0, yum = 0.2195149727645171d0
  
  real(kind=DP), parameter :: ninth2=0.22222222222d0
  real(kind=DP), parameter :: ninth4=0.44444444444d0
  real(kind=DP), parameter :: ninth8=0.88888888888d0
  real(kind=DP), parameter :: thrd2 =0.66666666666d0
  real(kind=DP), parameter :: eghth =0.125d0
  real(kind=DP), parameter :: twtysvnth8 = 0.296296296296d0
  real(kind=DP), parameter :: twtysvnth32= 1.18518518518519d0

  real(kind=DP) :: facw,d,dd,fk,s,fac,s2,yums2,f,ex,fs,exd,exdd,exc0,excd,excdd,exc1
  integer       :: is,i,nr
  real(kind=DP) :: fss,fsss,ivd,ivdthrd2

!---- Spin dependency

  facw = nspin
  exc  = 0.d0
  do is = 1, nspin
     exc1 = 0.d0
!     do i = 1,nrc,dnr
     do nr=1,num_ir
        i=irs(nr)
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
!           dFx_dgradrho(i, is) = excdd / grad_rho(i, is)           ! 2010 5
           dFx_dgradrho(i, is) = excdd
        else
           dFx_dgradrho(i, is) = 0.d0
        endif
! gradient of charge density 95/12/3 H.S.
!        exc1 = exc1 + exc0*wos(i)
        exc(i)=exc(i)+exc0
        
! *****  Higher derivatives *****
       
        fss = 2 * cupa ** 2 * yum * ( cupa - 3*yums2 ) / ( cupa + yums2 ) ** 3
        fsss= 24 * cupa **2 * yum **2 * s * ( yums2 - cupa ) /  ( cupa + yums2 ) ** 4
        
        ivd=1.d0/d
        ivdthrd2=ivd**thrd2
        
        if(nspin == 1) then
            dFx_drr(i,is)   =ninth4*ax*ivdthrd2*( f-s*fs+4.d0*fss*s2 )
            dFx_drg(i,is)   =-thrd2*ax/thpith*fss*s*ivd
            dFx_dgg(i,is)   =0.25d0*ax/thpith/thpith*fss*ivdthrd2**2
            dFx_drrr(i,is)  =-twtysvnth8*ax*ivdthrd2*ivd*( f-fs*s+18.d0*fss*s2+8.d0*fsss*s2*s )
            dFx_drrg(i,is)  =ninth2*ax*ivd*ivd/thpith*( 7.d0*fss*s+4.d0*fsss*s2 )
            dFx_drgg(i,is)  =-thrd*ax/thpith/thpith*ivd*ivdthrd2**2*( fss+fsss*s )
            dFx_dggg(i,is)  =eghth*ax/(thpith**3)*fsss*ivd**2*ivdthrd2
        else if(nspin == 2) then
            dFx_drr(i,is)   =ninth8*ax*ivdthrd2*( f-s*fs+4.d0*fss*s2 )
            dFx_drg(i,is)   =-thrd4*ax/thpith*fss*s*ivd
            dFx_dgg(i,is)   =0.5d0*ax/thpith/thpith*fss*ivdthrd2**2
            dFx_drrr(i,is)  =-twtysvnth32*ax*ivdthrd2*ivd*( f-fs*s+18.d0*fss*s2+8.d0*fsss*s2*s )
            dFx_drrg(i,is)  =ninth8*ax*ivd*ivd/thpith*( 7.d0*fss*s+4.d0*fsss*s2 )
            dFx_drgg(i,is)  =-thrd4*ax/(thpith**2)*ivd*ivdthrd2**2*( fss+fsss*s )
            dFx_dggg(i,is)  =0.5d0*ax/(thpith**3)*fsss*ivd**2*ivdthrd2
        end if
        
     end do
!     exc(i) = exc(i) + exc1
  end do
end subroutine ex_ggapbe_paw_drv3

subroutine cr_ggapbe_paw_drv2(nrc,dnr,nspin,chgrhr_l,grad_trho,exc,dF_drho,dF_dgradrho &
                                ,dFc_daa,dFc_dbb,dFc_dgg,dFc_dab,dFc_dag,dFc_dbg &
                                ,dFc_daaa,dFc_dbbb,dFc_dggg &
                                ,dFc_daab,dFc_daag,dFc_dabb,dFc_dbbg,dFc_dagg,dFc_dbgg,dFc_dabg)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : PAI,DP
  implicit none

  integer,intent(in)        :: nrc,dnr,nspin
  real(kind=DP),intent(in)  :: chgrhr_l(nrc,nspin)
  real(kind=DP),intent(in) :: grad_trho(nrc)
  real(kind=DP),intent(inout) :: exc(nrc)
  real(kind=DP),intent(out) :: dF_drho(nrc,nspin)
  real(kind=DP),intent(out) :: dF_dgradrho(nrc)
  real(kind=DP),intent(out) :: dFc_daa(nrc)
  real(kind=DP),intent(out) :: dFc_dbb(nrc)
  real(kind=DP),intent(out) :: dFc_dgg(nrc)
  real(kind=DP),intent(out) :: dFc_dab(nrc)
  real(kind=DP),intent(out) :: dFc_dag(nrc)
  real(kind=DP),intent(out) :: dFc_dbg(nrc)
  real(kind=DP),intent(out) :: dFc_daaa(nrc)
  real(kind=DP),intent(out) :: dFc_dbbb(nrc)
  real(kind=DP),intent(out) :: dFc_dggg(nrc)
  real(kind=DP),intent(out) :: dFc_daab(nrc)
  real(kind=DP),intent(out) :: dFc_daag(nrc)
  real(kind=DP),intent(out) :: dFc_dabb(nrc)
  real(kind=DP),intent(out) :: dFc_dbbg(nrc)
  real(kind=DP),intent(out) :: dFc_dagg(nrc)
  real(kind=DP),intent(out) :: dFc_dbgg(nrc)
  real(kind=DP),intent(out) :: dFc_dabg(nrc)

  real(kind=DP), parameter :: zeta_minimum    = 1.d-4
  real(kind=DP), parameter :: zeta_minimum2   = 1.d-6
  real(kind=DP), parameter :: density_minimum = 1.d-20
  real(kind=DP), parameter :: density_minimum2= 1.d-5

  real(kind=DP), parameter  :: a  = 0.0310907d0,    a1   = 0.21370d0
  real(kind=DP), parameter  :: b1 = 7.5957d0,       b2n  = 3.5876d0
  real(kind=DP), parameter  :: b3 = 1.6382d0,       b4   = 0.49294d0
  real(kind=DP), parameter  :: p  = 1.00d0,         p1   = p+1.d0

#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: ap  = 0.01554535d0,    a1p  = 0.20548d0
#else
  real(kind=DP), parameter  :: ap  = 0.015545d0,    a1p  = 0.20548d0
#endif

  real(kind=DP), parameter  :: b1p = 14.1189d0,     b2np = 6.1977d0
  real(kind=DP), parameter  :: b3p = 3.3662d0,      b4p  = 0.62517d0

! vwn 95/12/2 Y.M
#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: aq  = 0.0168869d0,    a1q  = 0.11125d0
#else
  real(kind=DP), parameter  :: aq  = 0.016887d0,    a1q  = 0.11125d0
#endif

  real(kind=DP), parameter  :: b1q = 10.3570d0,     b2nq = 3.6231d0
  real(kind=DP), parameter  :: b3q = 0.88026d0,     b4q  = 0.49671d0

! vwn 95/12/2 Y.M
#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: gamma = 0.03109069086965489503494086371273d0
#else
  real(kind=DP), parameter  :: gamma = 0.031091
#endif

!--------------------------

#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: bet=0.06672455060314922d0
  real(kind=DP), parameter  :: delt=bet/gamma

  real(kind=DP), parameter  :: thrd = 0.33333333333333330,&
       &                       sixth7 = 1.1666666666666666d0
  real(kind=DP), parameter  :: thrd7= 2.333333333333333333d0
  real(kind=DP), parameter  :: sixth13= 2.166666666666666666d0
  real(kind=DP), parameter  :: sixth19= 3.166666666666666666d0
  real(kind=DP), parameter  :: ninth= 0.11111111111111111111d0

#else
  real(kind=DP), parameter  :: xnu = 15.75592d0, cc0 = 0.004235d0
  real(kind=DP) :: bet, delt

  real(kind=DP), parameter  :: thrd = 0.333333333333d0, sixth7 = 1.1666666666666666d0
  real(kind=DP), parameter  :: thrd7= 2.333333333333d0
  real(kind=DP), parameter  :: sixth13= 2.1666666666666666d0
  real(kind=DP), parameter  :: sixth19= 3.1666666666666666d0
  real(kind=DP), parameter  :: ninth= 0.111111111111111d0
#endif

  real(kind=DP), parameter  :: twntysvnth=  3.703703703703703d-2
  real(kind=DP), parameter  :: fftyfrth= 1.851851851851852d-2
  real(kind=DP), parameter  :: eityonth= 1.234567901234568d-2

  integer       :: is,i
  real(kind=DP) :: facw, exc1,g,g3,facpon,d,dd,rs,fk,sk,t,s,q0 &
       &         , rs12,rs32,rsp,q1,q2,eu,q3,eurs,ecd,pon &
       &         , b,b2,t2,t4,t6,q4,q5,h0,h,q8,h0t,h0b,h0rs,ht,hrs &
       &         , ec1,ec1d,ec1dd,exc0,excd,excdd &
       &         , thrd2,thrd4,fzdd,zeta,onpzeta,onmzeta,fz,zetadxd &
       &         , q0p,q1p,q2p,q3p,q0q,q1q,q2q,q3q,fzd,onzt4,gzt,bzt,h0zt,euzt &
       &         , zetadxda,zetadxdb,ecda,ecdb,ec1da,ec1db,excda,excdb

  real(kind=DP) :: rsdn,rsdnn,rsdnnn
  real(kind=DP) :: q1drr,q1drrr
  real(kind=DP) :: eudrr,eudrrr,invq12pq1,q12p1
  real(kind=DP) :: eudn,eudnn,eudnnn,ecdd,ecddd
  real(kind=DP) :: dtdn,dtdg,dtdnn,dtdgg,dtdng,dtdnnn,dtdnng,dtdngg,dtdggg
  real(kind=DP) :: dAde,dAdee,dAdeee
  real(kind=DP) :: dAdn,dAdnn,dAdnnn
  real(kind=DP) :: fE,fEdt,fEdtt,fEda,fEdaa
  real(kind=DP) :: fF,fFdt,fFdtt,fFda,fFdaa
  real(kind=DP) :: fG,fGdt,fGdtt,fGda,fGdaa
  real(kind=DP) :: at2,a2t4,a3t6,t3,t5,t7,t8
  real(kind=DP) :: dHda,dHdt,dHdaa,dHdat,dHdtt,dHdaaa,dHdaat,dHdatt,dHdttt,betivE,betivE2
  real(kind=DP) :: dHdn,dHdg,dHdnn,dHdng,dHdgg,dHdnnn,dHdnng,dHdngg,dHdggg
  real(kind=DP) :: zetada,zetadb,zetadaa,zetadab,zetadbb,zetadaaa,zetadaab,zetadabb,zetadbbb
  real(kind=DP) :: id2,id3,id4
  real(kind=DP) :: e1dr,e1drr,e1drrr
  real(kind=DP) :: e1pdr,e1pdrr,e1pdrrr
  real(kind=DP) :: e1qdr,e1qdrr,e1qdrrr
  real(kind=DP) :: q1pdrr,q1pdrrr,invq1p2pq1p,q1p2p1
  real(kind=DP) :: q1qdrr,q1qdrrr,invq1q2pq1q,q1q2p1
  real(kind=DP) :: fzdd0,fzddd,opzthrdm2,omzthrdm2
  real(kind=DP) :: zeta2,zeta3,zeta4,eudrz,eudzz,eudrrz,eudrzz,eudzzz
  real(kind=DP) :: e1,e2,e3
  real(kind=DP) :: euda,eudb,eudaa,eudbb,eudab,eudaaa,eudaab,eudabb,eudbbb
  real(kind=DP) :: rsdn2,rsdn3,zetada2,zetada3,zetadb2,zetadb3
  real(kind=DP) :: gzz,gzzz,duda,dudb,dudaa,dudab,dudbb,dudaaa,dudaab,dudabb,dudbbb
  real(kind=DP) :: dtdu,dtduu,dtduuu,dtdnu,dtdnnu,dtdnuu,dtdug,dtduug,dtdugg,dtdnug,invg
  real(kind=DP) :: dtda,dtdb,dtdaa,dtdab,dtdag,dtdbg,dtdbb
  real(kind=DP) :: dtdaaa,dtdaab,dtdabb,dtdbbb,dtdaag,dtdagg,dtdbbg,dtdbgg,dtdabg
  real(kind=DP) :: invg3,m3u
  real(kind=DP) :: dAdu,dAduu,dAduuu,dAdeu,dAdeeu,dAdeuu
  real(kind=DP) :: dAda,dAdb,dAdaa,dAdab,dAdbb,dAdaaa,dAdaab,dAdabb,dAdbbb
  real(kind=DP) :: dHdf,dHdff,dHdft,dHdfff,dHdfft,dHdftt
  real(kind=DP) :: dHdu,dHduu,dHduuu,dHdfu,dHdtu,dHdffu,dHdfuu,dHdttu,dHdtuu,dHdftu
  real(kind=DP) :: tiu,siu2
  real(kind=DP) :: dHdb,dHdbb,dHdag,dHdbg,dHdab,dHdbbb,dHdaab,dHdabb,dHdbbg,dHdbgg,dHdaag,dHdagg,dHdabg


#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: eta = 1.d-12
  real(kind=DP), parameter  :: sixthm = -0.1666666666666666666666d0
#endif

  facw = nspin

#ifndef GGA_PAW_REVISED
  bet = xnu*cc0
  delt = bet / gamma
#endif

     if ( nspin == 1 ) then
        g   = 1.d0
        g3  = g**3
!!$        g4  = g3*g
        facpon = -delt/(g3*bet)

        do i = 1,nrc,dnr ! MPI
           d  = facw*chgrhr_l(i, 1)
           if(d < density_minimum) cycle
! gradient of charge density 95/12/3 H.S.
           dd = facw* grad_trho(i)
           rs = (0.75d0/(PAI*d))**thrd
           fk = (3.d0*PAI*PAI*d)**thrd
           sk = dsqrt(4.d0*fk/PAI)

           if(d > density_minimum2) then
              t = dd/(d*sk*2.d0)
              s = dd/(d*fk*2.d0)
           else
              t = 0.d0
              s = 0.d0
           endif

           q0 = -2.d0 *a *( 1.d0 +a1*rs )
           rs12 = dsqrt(rs)
           rs32 = rs12**3
           rsp = rs**p
           q1 = 2.d0 *a*( b1*rs12 +b2n*rs +b3*rs32 +b4*rs*rsp )
           q2 = log( 1.d0 +1.d0/q1 )
           eu = q0*q2
           q3 = a*( b1/rs12 +2.d0*b2n +3.d0 *b3*rs12 +2.d0*b4*p1*rsp )
           eurs = -2.d0*a*a1*q2 - q0*q3/(q1**2 + q1)
           ecd = eu-thrd*rs*eurs
           pon = facpon*eu
           b   = delt/(dexp(pon)-1.d0)
           b2  = b*b
           t2  = t*t
!!$           s2  = s*s
           t4  = t2*t2
           t6  = t4*t2
           q4  = 1.d0 + b*t2
           q5  = 1.d0 + b*t2 +b2*t4
           h0  = g3*(bet/delt)*dlog( 1.d0 +delt*q4*t2/q5 )
           h   = h0
           q8  = q5*q5 +delt*q4*q5*t2
           h0t = 2.d0*bet*t*( 1.d0 +2.d0*b*t2 )/q8 * g3
           h0b = -bet*t6*( 2.d0*b +b2*t2 )/q8 * g3
           h0rs = h0b*b*eurs*(b+delt)/ bet / g3
           ht   = h0t
           hrs  = h0rs
           ec1  = d*h/facw
           ec1d = h -thrd*rs*hrs -sixth7*t*ht
           ec1dd = 0.5d0*ht/sk

           exc0 = eu*d/facw + ec1
           excd = ecd + ec1d
           dF_drho(i, 1) = dF_drho(i, 1) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(dabs(grad_trho(i)) > 1.0d-9) then
!              dF_dgradrho(i) = excdd / grad_trho(i)
              dF_dgradrho(i) = excdd
           else
              dF_dgradrho(i) = 0.d0
           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d/facw
!           exc1 = exc1 + exc0*wos(i)
           exc(i)=exc(i)+exc0
        
! ***** Calc heigher derivatives. *****
           
           rsdn=-thrd*rs/d
           rsdnn=-4.d0*thrd*rsdn/d
           rsdnnn=-thrd7*rsdnn/d
!print '(i2,5e19.10)',i,d,rs,rsdn,rsdnn,rsdnnn
!q3 = a*(b1/rs12+2*b2n+3*b3*rs12+2*b4*p1*rsp)   
           q1drr =0.5d0*a*( -b1/rs32+3.d0*b3/rs12+8.d0*b4 )
           q1drrr=0.75d0*a*( b1/rs - b3 )/rs32
!eurs = -2*a*a1*q2 - q0*q3/(q1**2 + q1)
           invq12pq1=1.d0/(q1**2 + q1)
           q12p1=2.d0*q1+1.d0
           eudrr=4.d0*a*a1*q3*invq12pq1-q0*( q1drr-q12p1*q3**2*invq12pq1 )*invq12pq1
           eudrrr=6.d0*a*a1*( q1drr - q12p1*q3**2*invq12pq1 ) * invq12pq1 &
                    -q0*( q1drrr -(2.d0*q3**3+3.d0*q12p1*q3*q1drr -  &
                    2.d0*q12p1**2*q3**3* invq12pq1)*invq12pq1 ) * invq12pq1
           eudn=eurs*rsdn
           eudnn=eudrr*rsdn**2+eurs*rsdnn
           eudnnn=eudrrr*rsdn**3+3.d0*eudrr*rsdn*rsdnn+eurs*rsdnnn
!print '(i2,5e19.12)',i,d,eu,eudn,eudnn,eudnnn
           ecdd=2.d0*eudn+d*eudnn
           ecddd=3.d0*eudnn+d*eudnnn
!print '(i2,5e19.12)',i,d,eu*d,ecdd,ecddd
           
           dtdn  = -sixth7*t/d
           dtdnn = -sixth13*dtdn/d
           dtdnnn= -sixth19*dtdnn/d
!print '(i2,5e19.12)',i,d,t,dtdn,dtdnn,dtdnnn
           if(dabs(grad_trho(i)) > 1.0d-9) then
               dtdg  = t/dd
               dtdng = -sixth7*dtdg/d
               dtdnng= -sixth13*dtdng/d
           else
               dtdg  = 0.d0
               dtdng = 0.d0
               dtdnng= 0.d0
           endif
           dtdngg = 0.d0
           dtdggg = 0.d0
           
           dAde   = b * ( b + delt ) / bet
           dAdee  = 2.d0 * dAde * ( b + 0.5d0 * delt ) / bet
           dAdeee = 2.d0 * ( dAde**2 + ( b + 0.5d0 * delt )*dAdee ) /bet
           
           dAdn   = dAde  *eudn
           dAdnn  = dAdee *eudn**2  + dAde*eudnn
           dAdnnn = dAdeee*eudn**3  + 3.d0*dAdee*eudn*eudnn + dAde*eudnnn   
!print '(i2,5e19.10)' ,i,d,b,dAdnn
           
           at2  = b * t2
           a2t4 = at2**2
           a3t6 = at2*a2t4
           t3   = t2*t
           t5   = t4*t
           t7   = t6*t
           t8   = t7*t
           fE    = q5**2 + delt*t2*q4*q5
           fEdt  = 4.d0*b*t*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt *2.d0*t*(1.d0 + 4.d0*at2 + 6.d0*a2t4 + 4.d0*a3t6 )
           fEdtt = 4.d0*b*( 1.d0 + 9.d0*at2 + 15.d0*a2t4 + 14.d0*a3t6 ) + &
                    delt*( 2.d0 + 24.d0*at2 + 60.d0*a2t4 + 56.d0*a3t6 )
           fEda  = 2.d0*t2*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt*t4*( 2.d0 + 4.d0*at2 + 3.d0*a2t4 )
           fEdaa = 6.d0*t4*( 2.d0*q5 - 1.d0) + &
                    delt*2.d0*t6*( 2.d0 + 3.d0*at2 )
           fF    = t*( 1.d0 + 2.d0*at2 )         
           fFdt  = 1.d0 + 6.d0*at2
           fFdtt = 12.d0*b*t
           fFda  = 2.d0*t3
           fFdaa = 0.d0
           fG    = 2.d0*b*t6 + b2*t8
           fGdt  = 12.d0*b*t5 + 8.d0*b2*t7
           fGdtt = 60.d0*b*t4 + 56.d0*b2*t6
           fGda  = 2.d0*t6 + 2.d0*b*t8
           fGdaa = 2.d0*t8
           
           betivE=bet/fE
           betivE2=betivE/fE
           dHda   =     -fG*betivE
           dHdt   = 2.d0*fF*betivE
           dHdaa  =     -( fGda*fE-fG*fEda )*betivE2
           dHdat  =     -( fGdt*fE-fG*fEdt )*betivE2
           dHdtt  = 2.d0*( fFdt*fE-fF*fEdt )*betivE2
           dHdaaa =     -( fGdaa*fE-fG*fEdaa-2.d0*( fGda*fE-fG*fEda )*fEda/fE )*betivE2
           dHdaat = 2.d0*( fFdaa*fE-fF*fEdaa-2.d0*( fFda*fE-fF*fEda )*fEda/fE )*betivE2
           dHdatt =     -( fGdtt*fE-fG*fEdtt-2.d0*( fGdt*fE-fG*fEdt )*fEdt/fE )*betivE2
           dHdttt = 2.d0*( fFdtt*fE-fF*fEdtt-2.d0*( fFdt*fE-fF*fEdt )*fEdt/fE )*betivE2
           
           dHdn     = dHda  *dAdn       + dHdt*dtdn
           dHdg     = dHdt  *dtdg
           dHdnn    = dHdaa *dAdn**2    + 2.d0*dHdat*dAdn*dtdn  + dHdtt*dtdn**2   + dHda*dAdnn  + dHdt*dtdnn
           dHdng    = dHdat *dAdn*dtdg  + dHdtt*dtdn*dtdg       + dHdt*dtdng
           dHdgg    = dHdtt*dtdg**2
           dHdnnn   = dHdaaa*dAdn**3    + dHdttt*dtdn**3 &
                        + 3.d0*dHdaat*dAdn**2*dtdn  + 3.d0*dHdatt*dAdn*dtdn**2 &
                        + 3.d0*dHdaa*dAdn*dAdnn + 3.d0*dHdtt*dtdn*dtdnn + 3.d0*dHdat*(dAdnn*dtdn + dAdn*dtdnn) &
                        + dHda*dAdnnn   + dHdt*dtdnnn
           dHdnng   = dHdaat*dAdn**2*dtdg   + 2.d0*dHdatt*dAdn*dtdn*dtdg    + dHdttt*dtdn**2*dtdg &
                        + dHdat*(dAdnn*dtdg + 2.d0*dAdn*dtdng)  + dHdtt*(dtdnn*dtdg + 2.d0*dtdn*dtdng) &
                        + dHdt*dtdnng
           dHdngg   = dHdatt*dAdn*dtdg**2   + dHdttt*dtdn*dtdg**2   + 2.d0*dHdtt*dtdg*dtdng
           dHdggg   = dHdttt*dtdg**3
!dF_drho(i, 1)=h+d*dHdn+ecd
           dFc_daa(i)   = 2.d0*dHdn + d*dHdnn       + ecdd
           dFc_dag(i)   = dHdg  + d*dHdng
           dFc_dgg(i)   = d*dHdgg
           dFc_daaa(i)  = 3.d0*dHdnn    + d*dHdnnn  + ecddd
           dFc_daag(i)  = 2.d0*dHdng    + d*dHdnng
           dFc_dagg(i)  = dHdgg         + d*dHdngg
           dFc_dggg(i)  = d*dHdggg  
        end do

     else if ( nspin ==  2 ) then
        thrd2 = thrd * 2.d0
        thrd4 = thrd * 4.d0
        fzdd = 8.d0 /(9.d0*(2.d0**thrd4 - 2.d0 ))

        do i = 1,nrc,dnr
           d = chgrhr_l(i, 1) + chgrhr_l(i, nspin)
           if(d < density_minimum .or. chgrhr_l(i,1) < density_minimum .or. chgrhr_l(i,nspin)<density_minimum) cycle

           zeta = ( chgrhr_l(i, 1) - chgrhr_l(i, 2) ) / d
           onpzeta = 2.d0 *chgrhr_l(i, 1) / d
           onmzeta = 2.d0 *chgrhr_l(i, 2) / d
           g = ( onpzeta**thrd2 + onmzeta**thrd2 )/2
           if(dabs(zeta) < zeta_minimum ) then
              fz= 4.d0/9 *zeta*zeta*(1.d0 + 5.d0/54.d0 *zeta*zeta &
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
           ! ****** aas mkatsu *****
!              t = dd/(d*sk*2)
              t = dd/(d*sk*2.d0)/g         !    mkatsu aas
              s = dd/(d*fk*2.d0)
           else
              t = 0.d0
              s = 0.d0
           endif

           q0 = -2.d0*a*( 1.d0 +a1*rs )
           rs12 = dsqrt(rs)
           rs32 = rs12**3
           rsp = rs**p
           q1 = 2.d0 *a*( b1*rs12 +b2n*rs +b3*rs32 +b4*rs*rsp )
           q2 = log( 1.d0 +1.d0/q1 )
           q3 = a*( b1/rs12 +2.d0*b2n +3.d0*b3*rs12 +2.d0*b4*p1*rsp )

           q0p = -2.d0 *ap *( 1.d0 +a1p*rs )
           q1p = 2.d0 *ap*( b1p*rs12 +b2np*rs+b3p*rs32+b4p*rs*rsp )
           q2p = log( 1.d0 +1.d0/q1p )
           q3p = ap*( b1p/rs12 +2.d0*b2np +3.d0*b3p*rs12 +2.d0*b4p*p1*rsp )
! vwn 95/12/2 Y.M
           q0q = -2.d0 *aq*( 1.d0 +a1q*rs )
           q1q =  2.d0 *aq*( b1q*rs12 +b2nq*rs +b3q*rs32 +b4q*rs*rsp )
           q2q = log( 1.d0 +1.d0/q1q )
           q3q = aq*( b1q/rs12 +2.d0*b2nq +3.d0*b3q*rs12 +2.d0*b4q*p1*rsp )

           if(dabs(zeta) < zeta_minimum) then
              fzd = 8.d0/9.d0 /( 2.d0**thrd4 -2.d0 ) &
                   &     *zeta *( 1.d0 + 5.d0/27.d0 *zeta*zeta &
                   &           *( 1.d0 +22.d0/45.d0 *zeta*zeta ))
           else
              fzd = 4.d0/3.d0 /( 2.d0**thrd4 -2.d0 ) &
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

!           zetadxd = - 2 * chgrhr_l(i, 3-is) * (-1.d0)**is / d
           zetadxda = - 2.d0 * chgrhr_l(i, 2) * (-1.d0) / d
           zetadxdb = - 2.d0 * chgrhr_l(i, 1) / d
!           ecd = eu-thrd*rs*eurs + euzt * zetadxd
           ecda = eu-thrd*rs*eurs + euzt * zetadxda 
           ecdb = eu-thrd*rs*eurs + euzt * zetadxdb

#ifndef GGA_PAW_REVISED
           if(onpzeta < zeta_minimum2 .or. onmzeta < zeta_minimum2) then
              ec1  = 0.d0;    ec1d = 0.d0;     ec1dd= 0.d0
              ec1da = 0.d0;     ec1db = 0.d0
              go to 10
           end if
#endif

           pon = facpon*eu
           b   = delt/(dexp(pon)-1.d0)
           b2  = b*b
           t2  = t*t
!!$           s2  = s*s
           t4  = t2*t2
           t6  = t4*t2
           q4  = 1.d0 +b*t2
           q5  = 1.d0 +b*t2+b2*t4
           h0  = g3*(bet/delt)*dlog( 1.d0 +delt*q4*t2/q5 )
           h   = h0

           q8  = q5*q5+delt*q4*q5*t2
           h0t = 2.d0 *bet*t*( 1.d0 +2.d0 *b*t2 )/q8 * g3
           h0b = -bet*t6*( 2.d0 *b+b2*t2 )/q8 * g3
           h0rs = h0b * b * eurs * (b+delt) / bet / g3
           ht   = h0t
           hrs  = h0rs

           if(dabs(zeta) < zeta_minimum) then
              gzt = -2.d0/9*zeta *( 1.d0 +14.d0/27 *zeta*zeta &
                   &             *( 1.d0 +13.d0/18 *zeta*zeta ))
           else
#ifdef GGA_PAW_REVISED
              gzt = ( (onpzeta**2+eta)**sixthm - (onmzeta**2+eta)**sixthm)/3.d0
#else
              gzt = ( 1.d0/onpzeta**thrd -1.d0/onmzeta**thrd )/3.d0
#endif
           end if

           bzt = b * ( b+delt) * ( euzt - 3*eu*gzt/g ) / ( bet * g3 ) 
           h0zt = 3.d0 * gzt * h0 / g + h0b * bzt

           ec1  = d*h * 0.5d0

! --> revised by Kato-san,  1999Jul16
!!$           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + h0zt * zetadxd
!           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxd
           ec1da = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxda
           ec1db = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxdb
           ec1dd = 0.5d0*ht/sk / g

10         continue


           exc0 = eu*d * 0.5d0 + ec1
!           excd = ecd + ec1d
           excda= ecda + ec1da
           excdb= ecdb + ec1db
!           dF_drho(i, is) = dF_drho(i, is) + excd
           dF_drho(i, 1) = dF_drho(i, 1) + excda
           dF_drho(i, 2) = dF_drho(i, 2) + excdb
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
!           if(is ==  2) then
              if(dabs(grad_trho(i)) > 1.0d-9) then
!                 dF_dgradrho(i) = excdd / grad_trho(i)
                 dF_dgradrho(i) = excdd
              else
                 dF_dgradrho(i) = 0.d0
              endif
!           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d * 0.5
!           exc1 = exc1 + exc0*wos(i)
           exc(i)=exc(i)+exc0*2.d0
           
! ***** Calc heigher order derivatives. *****

           rsdn=-thrd*rs/d
           rsdnn=-4.d0*thrd*rsdn/d
           rsdnnn=-thrd7*rsdnn/d
           
           id2=1.d0/d/d
           id3=id2/d
           id4=id3/d
           zetada   =     zetadxda/d
           zetadb   =     zetadxdb/d
           zetadaa  =    -4.d0*chgrhr_l(i, 2)*id3
           zetadab  =     2.d0*zeta*id2
           zetadbb  =     4.d0*chgrhr_l(i, 1)*id3
           zetadaaa =    12.d0*chgrhr_l(i, 2)*id4
           zetadaab =   ( -4.d0*chgrhr_l(i, 1) + 8.d0*chgrhr_l(i, 2))*id4
           zetadabb =   ( -8.d0*chgrhr_l(i, 1) + 4.d0*chgrhr_l(i, 2))*id4
           zetadbbb =   -12.d0*chgrhr_l(i, 1)*id4
           
!print '(i2,7d19.10)', i,chgrhr_l(i,1:2),zetadab,zetadaaa,zetadaab,zetadabb,zetadbbb
             
           q1drr =0.5d0*a*( -b1/rs32+3.d0*b3/rs12+8.d0*b4 )
           q1drrr=0.75d0*a*( b1/rs - b3 )/rs32
           invq12pq1=1.d0/(q1**2 + q1)
           q12p1=2.d0*q1+1.d0
           e1dr=-2.d0*a*a1*q2-q0*q3*invq12pq1
           e1drr=4.d0*a*a1*q3*invq12pq1-q0*( q1drr-q12p1*q3**2*invq12pq1 )*invq12pq1
           e1drrr=6.d0*a*a1*( q1drr - q12p1*q3**2*invq12pq1 ) * invq12pq1 &
                    -q0*( q1drrr -(2.d0*q3**3+3.d0*q12p1*q3*q1drr -  &
                    2.d0*q12p1**2*q3**3* invq12pq1)*invq12pq1 ) * invq12pq1
                    
           q1pdrr =0.5d0*ap*( -b1p/rs32+3.d0*b3p/rs12+8.d0*b4p )
           q1pdrrr=0.75d0*ap*( b1p/rs - b3p )/rs32
           invq1p2pq1p=1.d0/(q1p**2 + q1p)
           q1p2p1=2.d0*q1p+1.d0
           e1pdr=-2.d0*ap*a1p*q2p-q0p*q3p*invq1p2pq1p
           e1pdrr=4.d0*ap*a1p*q3p*invq1p2pq1p-q0p*( q1pdrr-q1p2p1*q3p**2*invq1p2pq1p )*invq1p2pq1p
           e1pdrrr=6.d0*ap*a1p*( q1pdrr - q1p2p1*q3p**2*invq1p2pq1p ) * invq1p2pq1p &
                    -q0p*( q1pdrrr -(2.d0*q3p**3+3.d0*q1p2p1*q3p*q1pdrr -  &
                    2.d0*q1p2p1**2*q3p**3* invq1p2pq1p)*invq1p2pq1p ) * invq1p2pq1p
                    
           q1qdrr =0.5d0*aq*( -b1q/rs32+3.d0*b3q/rs12+8.d0*b4q )
           q1qdrrr=0.75d0*aq*( b1q/rs - b3q )/rs32
           invq1q2pq1q=1.d0/(q1q**2 + q1q)
           q1q2p1=2.d0*q1q+1.d0
           e1qdr=-2.d0*aq*a1q*q2q-q0q*q3q*invq1q2pq1q
           e1qdrr=4.d0*aq*a1q*q3q*invq1q2pq1q-q0q*( q1qdrr-q1q2p1*q3q**2*invq1q2pq1q )*invq1q2pq1q
           e1qdrrr=6.d0*aq*a1q*( q1qdrr - q1q2p1*q3q**2*invq1q2pq1q ) * invq1q2pq1q &
                    -q0q*( q1qdrrr -(2.d0*q3q**3+3.d0*q1q2p1*q3q*q1qdrr -  &
                    2.d0*q1q2p1**2*q3q**3* invq1q2pq1q)*invq1q2pq1q ) * invq1q2pq1q
           
           zeta2=zeta**2
           zeta3=zeta2*zeta
           zeta4=zeta2*zeta2
           
           if(dabs(zeta) < zeta_minimum) then
              fzdd0 = fzdd &
                   &            *( 1 + 5.d0*ninth *zeta2 &
                   &            *( 1 +22.d0*twntysvnth *zeta2 ))
              fzddd = fzdd*10.d0*ninth*zeta*( 1.d0 + 44.d0*twntysvnth*zeta2 )
              gzz   = -2.d0*ninth*( 1.d0 + 14.d0*ninth*zeta2*( 1.d0 + 65.d0*fftyfrth*zeta2 ) )
              gzzz  = -56.d0*eityonth*zeta*( 1.d0 + 65.d0*twntysvnth*zeta2 )
           else
              opzthrdm2=1.d0/onpzeta**(thrd2)
              omzthrdm2=1.d0/onmzeta**(thrd2)
              fzdd0 =  fzdd*0.5d0*( opzthrdm2 + omzthrdm2 )
              fzddd = -fzdd*thrd*( opzthrdm2/onpzeta - omzthrdm2/onmzeta )
              gzz   = -ninth*( opzthrdm2**2 + omzthrdm2**2 )
!              gzz   = -ninth*( 1.d0/onpzeta**(thrd4) + 1.d0/onmzeta**(thrd4) )
              gzzz  = 4.d0*twntysvnth*( opzthrdm2**2/onpzeta - omzthrdm2**2/onmzeta )
           end if
           
           e1=q0*q2
           e2=q0p*q2p
           e3=q0q*q2q
           
           eudrr    = e1drr  + ( e1pdrr  -  e1drr  ) * fz * zeta4 - e1qdrr  * fz * onzt4 / fzdd
           eudrrr   = e1drrr + ( e1pdrrr -  e1drrr ) * fz * zeta4 - e1qdrrr * fz * onzt4 / fzdd
           eudzz    = ( e2 - e1 )*( fzdd0*zeta4 + 8.d0*fzd*zeta3 + 12.d0*fz*zeta2 ) &
                        - e3*( fzdd0*onzt4 - 8.d0*fzd*zeta3 - 12.d0*fz*zeta2 )/fzdd
           eudzzz   = ( e2 - e1 )*( fzddd*zeta4 + 12.d0*fzdd0*zeta3 + 36.d0*fzd*zeta2 + 24.d0*fz*zeta ) &
                        - e3*( fzddd*onzt4 - 12.d0*fzdd0*zeta3 - 36.d0*fzd*zeta2 - 24.d0*fz*zeta )/fzdd
           eudrz    = ( e1pdr - e1dr )*( fzd*zeta4 + 4.d0*fz*zeta3 ) &
                        - e1qdr *( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd 
           eudrrz   = ( e1pdrr  - e1drr  )*( fzd*zeta4 + 4.d0*fz*zeta3 ) &
                        - e1qdrr *( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd
           eudrzz   = ( e1pdr - e1dr )*( fzdd0*zeta4 + 8.d0*fzd*zeta3 + 12.d0*fz*zeta2 ) &
                        - e1qdr *( fzdd0*onzt4 - 8.d0*fzd*zeta3 - 12.d0*fz*zeta2 )/fzdd
                        
           rsdn2    = rsdn*rsdn
           rsdn3    = rsdn2*rsdn
           zetada2  = zetada**2
           zetada3  = zetada2*zetada
           zetadb2  = zetadb**2
           zetadb3  = zetadb2*zetadb
                        
           euda     = eurs*rsdn + euzt*zetada
           eudb     = eurs*rsdn + euzt*zetadb
           eudaa    = eudrr*rsdn2 + 2.d0*eudrz*rsdn*zetada + eudzz*zetada2 + eurs*rsdnn + euzt*zetadaa
           eudab    = eudrr*rsdn2 + eudrz*rsdn*( zetada + zetadb ) + eudzz*zetada*zetadb &
                        + eurs*rsdnn + euzt*zetadab
           eudbb    = eudrr*rsdn2 + 2.d0*eudrz*rsdn*zetadb + eudzz*zetadb2 + eurs*rsdnn + euzt*zetadbb
           eudaaa   = eudrrr*rsdn3 &
                        + 3.d0*eudrrz*rsdn2*zetada &
                        + 3.d0*eudrzz*rsdn*zetada2 + eudzzz*zetada3 &
                        + 3.d0*eudrr*rsdn*rsdnn + 3.d0*eudzz*zetada*zetadaa &
                        + 3.d0*eudrz*( rsdnn*zetada + rsdn*zetadaa ) &
                        + eurs*rsdnnn + euzt*zetadaaa
           eudaab   = eudrrr*rsdn3 &
                        + eudrrz*rsdn2*( zetadb + 2.d0*zetada ) &
                        + eudrzz*rsdn*zetada*( zetada + 2.d0*zetadb ) &
                        + eudzzz*zetada2*zetadb &
                        + 3.d0*eudrr*rsdn*rsdnn &
                        + eudrz*( rsdnn*( 2.d0*zetada + zetadb ) + rsdn*( 2.d0*zetadab + zetadaa ) ) &
                        + eudzz*( 2.d0*zetada*zetadab + zetadaa*zetadb ) &
                        + eurs*rsdnnn + euzt*zetadaab
           eudabb   = eudrrr*rsdn3 &
                        + eudrrz*rsdn2*( zetada + 2.d0*zetadb ) &
                        + eudrzz*rsdn*zetadb*( zetadb + 2.d0*zetada ) &
                        + eudzzz*zetada*zetadb2 &
                        + 3.d0*eudrr*rsdn*rsdnn &
                        + eudrz*( rsdnn*( zetada + 2.d0*zetadb ) + rsdn*( 2.d0*zetadab + zetadbb ) ) &
                        + eudzz*( zetada*zetadbb + 2.d0*zetadb*zetadab ) &
                        + eurs*rsdnnn + euzt*zetadabb
           eudbbb   = eudrrr*rsdn3 &
                        + 3.d0*eudrrz*rsdn2*zetadb &
                        + 3.d0*eudrzz*rsdn*zetadb2 + eudzzz*zetadb3 &
                        + 3.d0*eudrr*rsdn*rsdnn + 3.d0*eudzz*zetadb*zetadbb &
                        + 3.d0*eudrz*( rsdnn*zetadb + rsdn*zetadbb ) &
                        + eurs*rsdnnn + euzt*zetadbbb
                        
!print '(i2,12d19.10)', i,chgrhr_l(i,1:2),eu,euda,eudb,eudaa,eudbb,eudab,eudaaa,eudaab,eudabb,eudbbb
!print *,eu,     e1+(e2-e1)* fz * zeta4 - e3  * fz * onzt4 / fzdd
!print *,eurs,   e1dr  + ( e1pdr  -  e1dr  ) * fz * zeta4 - e1qdr  * fz * onzt4 / fzdd
!print *,euzt,   (e2-e1)*( fzd*zeta4 + 4.d0*fz*zeta3 )-e3*( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd 

           duda     = gzt*zetada
           dudb     = gzt*zetadb 
           dudaa    = gzz*zetada2     + gzt*zetadaa
           dudab    = gzz*zetada*zetadb + gzt*zetadab
           dudbb    = gzz*zetadb2     + gzt*zetadbb
           dudaaa   = gzzz*zetada3    + 3.d0*gzz*zetada*zetadaa   + gzt*zetadaaa
           dudaab   = gzzz*zetada2*zetadb + gzz*( zetadaa*zetadb + 2.d0*zetada*zetadab ) + gzt*zetadaab
           dudabb   = gzzz*zetada*zetadb2 + gzz*( zetada*zetadbb + 2.d0*zetadab*zetadb ) + gzt*zetadabb
           dudbbb   = gzzz*zetadb3    + 3.d0*gzz*zetadb*zetadbb   + gzt*zetadbbb 
           
!print '(i2,12d19.10)', i,chgrhr_l(i,1:2),g,duda,dudb,dudaa,dudbb,dudab,dudaaa,dudaab,dudabb,dudbbb

           
           dtdn  = -sixth7*t/d
           dtdnn = -sixth13*dtdn/d
           dtdnnn= -sixth19*dtdnn/d
!print '(i2,5e19.12)',i,d,t,dtdn,dtdnn,dtdnnn
           if(dabs(grad_trho(i)) > 1.0d-9) then
               dtdg  = t/dd
               dtdng = -sixth7*dtdg/d
               dtdnng= -sixth13*dtdng/d
!               dtdnug = -dtdn*invg/dd
           else
               dtdg  = 0.d0
               dtdng = 0.d0
               dtdnng= 0.d0
               dtdnug= 0.d0 
           endif
           dtdgg  = 0.d0
           dtdngg = 0.d0
           dtdggg = 0.d0
           invg=1.d0/g
           dtdu   = -t*invg
           dtduu  = -2.d0*dtdu*invg 
           dtduuu = -3.d0*dtduu*invg
           dtdnu  = -dtdn*invg
           dtdnnu = -dtdnn*invg
           dtdnuu = -sixth7*dtduu/d
           dtdug  = -dtdg*invg
           dtduug = 2.d0*dtdg*invg**2
           dtdugg = 0.d0
           if(dabs(grad_trho(i)) > 1.0d-9) then
                dtdnug = dtdnu/dd
           else
                dtdnug = 0.d0
           end if
           
           dtda   = dtdn + dtdu*duda
           dtdb   = dtdn + dtdu*dudb
           dtdaa  = dtdnn + 2.d0*dtdnu*duda + dtduu*duda**2 + dtdu*dudaa
           dtdbb  = dtdnn + 2.d0*dtdnu*dudb + dtduu*dudb**2 + dtdu*dudbb
           dtdab  = dtdnn + dtdnu*( duda + dudb ) + dtduu*duda*dudb + dtdu*dudab
           dtdag  = dtdng + dtdug*duda
           dtdbg  = dtdng + dtdug*dudb
           dtdaaa = dtdnnn+ 3.d0*dtdnnu*duda + 3.d0*dtdnuu*duda**2 + dtduuu*duda**3 &
                    + 3.d0*dtduu*duda*dudaa + 3.d0*dtdnu*dudaa + dtdu*dudaaa
           dtdbbb = dtdnnn+ 3.d0*dtdnnu*dudb + 3.d0*dtdnuu*dudb**2 + dtduuu*dudb**3 &
                    + 3.d0*dtduu*dudb*dudbb + 3.d0*dtdnu*dudbb + dtdu*dudbbb
           dtdaab = dtdnnn+ dtdnnu*( 2.d0*duda + dudb ) + dtdnuu*duda*( duda + 2.d0*dudb ) + dtduuu*duda**2*dudb &
                    + dtduu*( dudaa*dudb + 2.d0*duda*dudab ) + dtdnu*( dudaa + 2.d0*dudab ) + dtdu*dudaab
           dtdabb = dtdnnn+ dtdnnu*( duda + 2.d0*dudb ) + dtdnuu*dudb*( dudb + 2.d0*duda ) + dtduuu*dudb**2*duda &
                    + dtduu*( dudbb*duda + 2.d0*dudb*dudab ) + dtdnu*( dudbb + 2.d0*dudab ) + dtdu*dudabb 
           dtdaag = dtdnng+ 2.d0*dtdnug*duda + dtduug*duda**2 + dtdug*dudaa
           dtdbbg = dtdnng+ 2.d0*dtdnug*dudb + dtduug*dudb**2 + dtdug*dudbb
           dtdabg = dtdnng+ dtdnug*( duda + dudb ) + dtduug*duda*dudb + dtdug*dudab 
           
!print '(i2,20d19.10)', i,t,dtda,dtdb,dtdg,dtdaa,dtdbb,dtdgg,dtdab,dtdag,dtdbg,dtdaaa,dtdbbb,dtdggg,dtdaab,dtdabb,dtdaag,dtdagg,dtdbbg,dtdbgg,dtdabg     
!print '(i2,12d19.10)', i,t,dtdag,dtdbg,dtdaag,dtdbbg,dtdabg     

           invg3  = 1.d0  /g3
           dAde   = b * ( b + delt ) / bet * invg3
           dAdee  = 2.d0 * dAde * ( b + 0.5d0 * delt ) / bet *invg3
           dAdeee = 2.d0 * ( dAde**2 + ( b + 0.5d0 * delt )*dAdee ) /bet * invg3
           m3u    = -3.d0*invg
           dAdu     = m3u*eu*dAde
           dAdeu    = m3u*( dAde + eu*dAdee )
           dAduu    = m3u*eu*invg*( -dAde + g*dAdeu )
           dAdeeu   = m3u*( 2.d0*dAdee + eu*dAdeee )
           dAdeuu   = m3u*invg*( -dAde + g*dAdeu - eu*dAdee + g*eu*dAdeeu )
           dAduuu   = m3u*eu*invg**2*( 2.d0*dAde -2.d0*g*dAdeu + g**2*dAdeuu )
           
           dAda     = dAde*euda + dAdu*duda
           dAdb     = dAde*eudb + dAdu*dudb
           dAdaa    = dAdee*euda**2 + 2.d0*dAdeu*euda*duda + dAduu*duda**2 + dAde*eudaa + dAdu*dudaa 
           dAdab    = dAdee*euda*eudb + dAdeu*( euda*dudb + eudb*duda ) + dAduu*duda*dudb + dAde*eudab + dAdu*dudab 
           dAdbb    = dAdee*eudb**2 + 2.d0*dAdeu*eudb*dudb + dAduu*dudb**2 + dAde*eudbb + dAdu*dudbb
           dAdaaa   = dAdeee*euda**3 + 3.d0*dAdeeu*euda**2*duda + 3.d0*dAdeuu*euda*duda**2 + dAduuu*duda**3 &
                        + 3.d0*dAdee*euda*eudaa + 3.d0*dAduu*duda*dudaa &
                        + 3.d0*dAdeu*( euda*dudaa + eudaa*duda ) + dAde*eudaaa + dAdu*dudaaa
           dAdaab   = dAdeee*euda**2*eudb + dAduuu*duda**2*dudb &
                        + dAdeeu*euda*( euda*dudb + 2.d0*eudb*duda ) &
                        + dAdeuu*( 2.d0*euda*dudb + eudb*duda )*duda &
                        + dAdee*( 2.d0*euda*eudab + eudb*eudaa ) + dAduu*( 2.d0*duda*dudab + dudb*dudaa ) &
                        + dAdeu*( 2.d0*eudab*duda + 2.d0*euda*dudab + dudb*eudaa + eudb*dudaa ) &
                        + dAde*eudaab + dAdu*dudaab
           dAdabb   = dAdeee*euda*eudb**2 + dAduuu*dudb**2*duda &
                        + dAdeeu*eudb*( eudb*duda + 2.d0*euda*dudb ) &
                        + dAdeuu*( 2.d0*eudb*duda + euda*dudb )*dudb &
                        + dAdee*( eudbb*euda + 2.d0*eudb*eudab ) + dAduu*( dudbb*duda + 2.d0*dudb*dudab ) &
                        + dAdeu*( eudbb*duda + 2.d0*eudb*dudab + 2.d0*dudb*eudab + dudbb*euda ) &
                        + dAde*eudabb + dAdu*dudabb
           dAdbbb   = dAdeee*eudb**3 + 3.d0*dAdeeu*eudb**2*dudb + 3.d0*dAdeuu*eudb*dudb**2 + dAduuu*dudb**3 &
                        + 3.d0*dAdee*eudb*eudbb + 3.d0*dAduu*dudb*dudbb &
                        + 3.d0*dAdeu*( eudb*dudbb + eudbb*dudb ) + dAde*eudbbb + dAdu*dudbbb 
                        
!print '(i2,12d19.10)', i,b,dAda,dAdb,dAdaa,dAdab,dAdbb,dAdaaa,dAdaab,dAdabb,dAdbbb

           at2  = b * t2
           a2t4 = at2**2
           a3t6 = at2*a2t4
           t3   = t2*t
           t5   = t4*t
           t7   = t6*t
           t8   = t7*t
           fE    = q5**2 + delt*t2*q4*q5
           fEdt  = 4.d0*b*t*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt *2.d0*t*(1.d0 + 4.d0*at2 + 6.d0*a2t4 + 4.d0*a3t6 )
           fEdtt = 4.d0*b*( 1.d0 + 9.d0*at2 + 15.d0*a2t4 + 14.d0*a3t6 ) + &
                    delt*( 2.d0 + 24.d0*at2 + 60.d0*a2t4 + 56.d0*a3t6 )
           fEda  = 2.d0*t2*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt*t4*( 2.d0 + 4.d0*at2 + 3.d0*a2t4 )
           fEdaa = 6.d0*t4*( 2.d0*q5 - 1.d0) + &
                    delt*2.d0*t6*( 2.d0 + 3.d0*at2 )
           fF    = t*( 1.d0 + 2.d0*at2 )         
           fFdt  = 1.d0 + 6.d0*at2
           fFdtt = 12.d0*b*t
           fFda  = 2.d0*t3
           fFdaa = 0.d0
           fG    = 2.d0*b*t6 + b2*t8
           fGdt  = 12.d0*b*t5 + 8.d0*b2*t7
           fGdtt = 60.d0*b*t4 + 56.d0*b2*t6
           fGda  = 2.d0*t6 + 2.d0*b*t8
           fGdaa = 2.d0*t8
           
           betivE=bet/fE*g3
           betivE2=betivE/fE
           dHdf   =     -fG*betivE
           dHdt   = 2.d0*fF*betivE
           dHdff  =     -( fGda*fE-fG*fEda )*betivE2
           dHdft  =     -( fGdt*fE-fG*fEdt )*betivE2
           dHdtt  = 2.d0*( fFdt*fE-fF*fEdt )*betivE2
           dHdfff =     -( fGdaa*fE-fG*fEdaa-2.d0*( fGda*fE-fG*fEda )*fEda/fE )*betivE2
           dHdfft = 2.d0*( fFdaa*fE-fF*fEdaa-2.d0*( fFda*fE-fF*fEda )*fEda/fE )*betivE2
           dHdftt =     -( fGdtt*fE-fG*fEdtt-2.d0*( fGdt*fE-fG*fEdt )*fEdt/fE )*betivE2
           dHdttt = 2.d0*( fFdtt*fE-fF*fEdtt-2.d0*( fFdt*fE-fF*fEdt )*fEdt/fE )*betivE2
           
           tiu    = 3.d0*invg
           siu2   = tiu*2.d0*invg
           dHdu   = tiu*h
           dHduu  = siu2*h
           dHduuu = siu2*invg*h
           dHdfu  = tiu*dHdf
           dHdtu  = tiu*dHdt
           dHdffu = tiu*dHdff
           dHdfuu = siu2*dHdf
           dHdttu = tiu*dHdtt
           dHdtuu = siu2*dHdt
           dHdftu = tiu*dHdft
           
           dHda   = dHdf*dAda + dHdt*dtda + dHdu*duda
           dHdb   = dHdf*dAdb + dHdt*dtdb + dHdu*dudb
           dHdg   =             dHdt*dtdg
           dHdaa  = dHdff*dAda**2 + dHdtt*dtda**2 + dHduu*duda**2 &
                    + 2.d0*dHdft*dAda*dtda + 2.d0*dHdfu*dAda*duda + 2.d0*dHdtu*dtda*duda &
                    + dHdf*dAdaa + dHdt*dtdaa + dHdu*dudaa
           dHdab  = dHdff*dAda*dAdb + dHdtt*dtda*dtdb + dHduu*duda*dudb &
                    + dHdft*( dAda*dtdb + dAdb*dtda ) + dHdfu*( dAda*dudb + dAdb*duda ) + dHdtu*( dtda*dudb + dtdb*duda ) &
                    + dHdf*dAdab + dHdt*dtdab + dHdu*dudab
           dHdbb  = dHdff*dAdb**2 + dHdtt*dtdb**2 + dHduu*dudb**2 &
                    + 2.d0*dHdft*dAdb*dtdb + 2.d0*dHdfu*dAdb*dudb + 2.d0*dHdtu*dtdb*dudb &
                    + dHdf*dAdbb + dHdt*dtdbb + dHdu*dudbb
           dHdgg  = dHdtt*dtdg**2
           dHdag  = dHdft*dAda*dtdg + dHdtt*dtda*dtdg + dHdtu*duda*dtdg + dHdt*dtdag
           dHdbg  = dHdft*dAdb*dtdg + dHdtt*dtdb*dtdg + dHdtu*dudb*dtdg + dHdt*dtdbg
           dHdaaa = dHdfff*dAda**3 + 3.d0*dHdfft*dAda**2*dtda + 3.d0*dHdffu*dAda**2*duda &
                    + 3.d0*dHdftt*dAda*dtda**2 + dHdttt*dtda**3 + 3.d0*dHdttu*dtda**2*duda &
                    + 3.d0*dHdfuu*dAda*duda**2 + 3.d0*dHdtuu*dtda*duda**2 + dHduuu*duda**3 &
                    + 6.d0*dHdftu*dAda*dtda*duda &
                    + 3.d0*dHdff*dAda*dAdaa + 3.d0*dHdtt*dtda*dtdaa + 3.d0*dHduu*duda*dudaa &
                    + 3.d0*dHdft*( dAdaa*dtda + dAda*dtdaa ) &
                    + 3.d0*dHdtu*( dtdaa*duda + dtda*dudaa ) &
                    + 3.d0*dHdfu*( dAdaa*duda + dAda*dudaa ) &
                    + dHdf*dAdaaa + dHdt*dtdaaa + dHdu*dudaaa
           dHdaab = dHdfff*dAda**2*dAdb + dHdfft*dAda*( dAda*dtdb + 2.d0*dAdb*dtda ) &
                    + dHdffu*dAda*( dAda*dudb + 2.d0*dAdb*duda ) &
                    + dHdftt*dtda*( dAdb*dtda + 2.d0*dAda*dtdb ) &
                    + dHdttt*dtdb*dtda**2 + dHdttu*dtda*( dtda*dudb + 2.d0*dtdb*duda ) &
                    + dHdfuu*duda*( dAdb*duda + 2.d0*dAda*dudb ) &
                    + dHdtuu*duda*( dtdb*duda + 2.d0*dtda*dudb ) + dHduuu*duda**2*dudb &
                    + 2.d0*dHdftu*( dAdb*duda*dtda + dAda*dudb*dtda + dAda*duda*dtdb ) &
                    + dHdff*( 2.d0*dAda*dAdab + dAdb*dAdaa ) &
                    + dHdtt*( 2.d0*dtda*dtdab + dtdb*dtdaa ) &
                    + dHduu*( 2.d0*duda*dudab + dudb*dudaa ) &
                    + dHdft*( 2.d0*dAdab*dtda + 2.d0*dAda*dtdab &
                    + dAdaa*dtdb + dAdb*dtdaa ) &
                    + dHdtu*( 2.d0*dtdab*duda + 2.d0*dtda*dudab &
                    + dtdaa*dudb + dtdb*dudaa ) &
                    + dHdfu*( 2.d0*dAdab*duda + 2.d0*dAda*dudab &
                    + dAdaa*dudb + dAdb*dudaa ) &
                    + dHdf*dAdaab + dHdt*dtdaab + dHdu*dudaab
           dHdabb = dHdfff*dAda*dAdb**2 + dHdfft*dAdb*( dAdb*dtda + 2.d0*dAda*dtdb ) &
                    + dHdffu*dAdb*( dAdb*duda + 2.d0*dAda*dudb ) &
                    + dHdftt*dtdb*( dAda*dtdb + 2.d0*dAdb*dtda ) + dHdttt*dtda*dtdb**2 &
                    + dHdttu*dtdb*( dtdb*duda + 2.d0*dtda*dudb ) &
                    + dHdfuu*dudb*( dAda*dudb + 2.d0*dAdb*duda ) &
                    + dHdtuu*dudb*( dtda*dudb + 2.d0*dtdb*duda ) + dHduuu*duda*dudb**2 &
                    + 2.d0*dHdftu*( dAda*dtdb*dudb + dAdb*dtda*dudb + dAdb*dtdb*duda ) &
                    + dHdff*( 2.d0*dAdb*dAdab + dAda*dAdbb ) &
                    + dHdtt*( 2.d0*dtdb*dtdab + dtda*dtdbb ) &
                    + dHduu*( 2.d0*dudb*dudab + duda*dudbb ) &
                    + dHdft*( 2.d0*dAdab*dtdb + 2.d0*dAdb*dtdab + dAdbb*dtda + dAda*dtdbb ) &
                    + dHdtu*( 2.d0*dtdab*dudb + 2.d0*dtdb*dudab + dtdbb*duda + dtda*dudbb ) &
                    + dHdfu*( 2.d0*dAdab*dudb + 2.d0*dAdb*dudab + dAdbb*duda + dAda*dudbb ) &
                    + dHdf*dAdabb + dHdt*dtdabb + dHdu*dudabb
           dHdbbb = dHdfff*dAdb**3 + 3.d0*dHdfft*dAdb**2*dtdb + 3.d0*dHdffu*dAdb**2*dudb &
                    + 3.d0*dHdftt*dAdb*dtdb**2 + dHdttt*dtdb**3 + 3.d0*dHdttu*dtdb**2*dudb &
                    + 3.d0*dHdfuu*dAdb*dudb**2 + 3.d0*dHdtuu*dtdb*dudb**2 + dHduuu*dudb**3 &
                    + 6.d0*dHdftu*dAdb*dtdb*dudb &
                    + 3.d0*dHdff*dAdb*dAdbb + 3.d0*dHdtt*dtdb*dtdbb + 3.d0*dHduu*dudb*dudbb &
                    + 3.d0*dHdft*( dAdbb*dtdb + dAdb*dtdbb ) &
                    + 3.d0*dHdtu*( dtdbb*dudb + dtdb*dudbb ) &
                    + 3.d0*dHdfu*( dAdbb*dudb + dAdb*dudbb ) &
                    + dHdf*dAdbbb + dHdt*dtdbbb + dHdu*dudbbb
           dHdabg = dHdfft*dAda*dAdb*dtdg   + dHdftt*( dAdb*dtda + dAda*dtdb )*dtdg &
                    + dHdttt*dtda*dtdb*dtdg + dHdttu*( dtdb*duda + dtda*dudb )*dtdg &
                    + dHdtuu*duda*dudb*dtdg + dHdftu*( dAdb*duda + dAda*dudb )*dtdg &
                    + dHdft*( dAdab*dtdg + dAdb*dtdag + dAda*dtdbg ) &
                    + dHdtt*( dtdab*dtdg + dtdb*dtdag + dtda*dtdbg ) &
                    + dHdtu*( dudab*dtdg + dudb*dtdag + duda*dtdbg ) &
                    + dHdt*dtdabg
           dHdaag = dHdfft*dAda**2*dtdg + dHdttt*dtda**2*dtdg + dHdtuu*duda**2*dtdg &
                    + 2.d0*dHdttu*dtda*duda*dtdg + 2.d0*dHdftt*dAda*dtda*dtdg + 2.d0*dHdftu*dAda*duda*dtdg &
                    + dHdtt*( dtdaa*dtdg + 2.d0*dtda*dtdag ) + dHdft*( dAdaa*dtdg + 2.d0*dAda*dtdag ) &
                    + dHdtu*( dudaa*dtdg + 2.d0*duda*dtdag ) + dHdt*dtdaag         
           dHdbbg = dHdfft*dAdb**2*dtdg + dHdttt*dtdb**2*dtdg + dHdtuu*dudb**2*dtdg &
                    + 2.d0*dHdttu*dtdb*dudb*dtdg + 2.d0*dHdftt*dAdb*dtdb*dtdg + 2.d0*dHdftu*dAdb*dudb*dtdg &
                    + dHdtt*( dtdbb*dtdg + 2.d0*dtdb*dtdbg ) + dHdft*( dAdbb*dtdg + 2.d0*dAdb*dtdbg ) &
                    + dHdtu*( dudbb*dtdg + 2.d0*dudb*dtdbg ) + dHdt*dtdbbg       
           dHdagg = dHdftt*dAda*dtdg**2 + dHdttt*dtda*dtdg**2 + dHdttu*duda*dtdg**2 + 2.d0*dHdtt*dtdg*dtdag
           dHdbgg = dHdftt*dAdb*dtdg**2 + dHdttt*dtdb*dtdg**2 + dHdttu*dudb*dtdg**2 + 2.d0*dHdtt*dtdg*dtdbg 
           dHdggg = dHdttt*dtdg**3        
                    
!print '(i2,20d19.10)', i,h,dHda,dHdb,dHdg,dHdaa,dHdbb,dHdgg,dHdag,dHdbg,dHdab,dHdaaa,dHdbbb,dHdggg,dHdaab,dHdabb,dHdaag,dHdagg,dHdbbg,dHdbgg,dHdabg
    
           dFc_daa(i)   = 2.d0*euda + d*eudaa + 2.d0*dHda + d*dHdaa 
           dFc_dbb(i)   = 2.d0*eudb + d*eudbb + 2.d0*dHdb + d*dHdbb
           dFc_dgg(i)   = d*dHdgg
           dFc_dab(i)   = euda + eudb + d*eudab + dHda + dHdb + d*dHdab 
           dFc_dag(i)   = dHdg + d*dHdag
           dFc_dbg(i)   = dHdg + d*dHdbg
           dFc_daaa(i)  = 3.d0*eudaa + d*eudaaa + 3.d0*dHdaa + d*dHdaaa
           dFc_dbbb(i)  = 3.d0*eudbb + d*eudbbb + 3.d0*dHdbb + d*dHdbbb
           dFc_dggg(i)  = d*dHdggg
           dFc_daab(i)  = eudaa + 2.d0*eudab +d*eudaab + dHdaa + 2.d0*dHdab + d*dHdaab
           dFc_daag(i)  = 2.d0*dHdag + d*dHdaag
           dFc_dabb(i)  = eudbb + 2.d0*eudab + d*eudabb + dHdbb + 2.d0*dHdab + d*dHdabb
           dFc_dbbg(i)  = 2.d0*dHdbg + d*dHdbbg
           dFc_dagg(i)  = dHdgg + d*dHdagg
           dFc_dbgg(i)  = dHdgg + d*dHdbgg
           dFc_dabg(i)  = dHdag + dHdbg + d*dHdabg
                    
!print '(i2,21d19.10)', i,exc(i),d*(eu+h),dF_drho(i,1),dF_drho(i,2),dF_dgradrho(i)*grad_trho(i) &
!                      ,dFc_daa(i),dFc_dbb(i),dFc_dgg(i),dFc_dab(i),dFc_dag(i),dFc_dbg(i),dFc_daaa(i),dFc_dbbb(i) &
!                      ,dFc_dggg(i),dFc_daab(i),dFc_daag(i),dFc_dabb(i),dFc_dbbg(i),dFc_dagg(i),dFc_dbgg(i),dFc_dabg(i)
!print *,dHdag,dHdaag,d,2.d0*dHdag + d*dHdaag
!print *,2.d0*dHdbg + d*dHdbbg
!print *,dHdag + dHdbg + d*dHdabg
               
           
        end do
     end if
     
end subroutine cr_ggapbe_paw_drv2
! === For nrc decomposion. by takto 2012/12/05 =================================
subroutine cr_ggapbe_paw_drv2_3D(nrc,dnr,nspin,chgrhr_l,grad_trho,exc,dF_drho,dF_dgradrho &
                                ,dFc_daa,dFc_dbb,dFc_dgg,dFc_dab,dFc_dag,dFc_dbg &
                                ,dFc_daaa,dFc_dbbb,dFc_dggg &
                                ,dFc_daab,dFc_daag,dFc_dabb,dFc_dbbg,dFc_dagg,dFc_dbgg,dFc_dabg &
                                ,ista_nrc,iend_nrc,ist,ien)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : PAI,DP
  implicit none

  integer,intent(in)        :: nrc,dnr,nspin
  real(kind=DP),intent(in)  :: chgrhr_l(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(in) :: grad_trho(ista_nrc:iend_nrc)
  real(kind=DP),intent(inout) :: exc(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dF_drho(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dF_dgradrho(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_daa(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbb(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dgg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dab(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dag(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_daaa(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbbb(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dggg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_daab(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_daag(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dabb(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbbg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dagg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbgg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dabg(ista_nrc:iend_nrc)
  integer :: ista_nrc, iend_nrc, ist, ien

  real(kind=DP), parameter :: zeta_minimum    = 1.d-4
  real(kind=DP), parameter :: zeta_minimum2   = 1.d-6
  real(kind=DP), parameter :: density_minimum = 1.d-20
  real(kind=DP), parameter :: density_minimum2= 1.d-5

  real(kind=DP), parameter  :: a  = 0.0310907d0,    a1   = 0.21370d0
  real(kind=DP), parameter  :: b1 = 7.5957d0,       b2n  = 3.5876d0
  real(kind=DP), parameter  :: b3 = 1.6382d0,       b4   = 0.49294d0
  real(kind=DP), parameter  :: p  = 1.00d0,         p1   = p+1.d0

#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: ap  = 0.01554535d0,    a1p  = 0.20548d0
#else
  real(kind=DP), parameter  :: ap  = 0.015545d0,    a1p  = 0.20548d0
#endif

  real(kind=DP), parameter  :: b1p = 14.1189d0,     b2np = 6.1977d0
  real(kind=DP), parameter  :: b3p = 3.3662d0,      b4p  = 0.62517d0
! vwn 95/12/2 Y.M
! vwn 95/12/2 Y.M
#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: aq  = 0.0168869d0,    a1q  = 0.11125d0
#else
  real(kind=DP), parameter  :: aq  = 0.016887d0,    a1q  = 0.11125d0
#endif

  real(kind=DP), parameter  :: b1q = 10.3570d0,     b2nq = 3.6231d0
  real(kind=DP), parameter  :: b3q = 0.88026d0,     b4q  = 0.49671d0

! vwn 95/12/2 Y.M
#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: gamma = 0.03109069086965489503494086371273d0
#else
  real(kind=DP), parameter  :: gamma = 0.031091
#endif

!--------------------------

#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: bet=0.06672455060314922d0
  real(kind=DP), parameter  :: delt=bet/gamma

  real(kind=DP), parameter  :: thrd = 0.33333333333333330,&
       &                       sixth7 = 1.1666666666666666d0
  real(kind=DP), parameter  :: thrd7= 2.333333333333333333d0
  real(kind=DP), parameter  :: sixth13= 2.166666666666666666d0
  real(kind=DP), parameter  :: sixth19= 3.166666666666666666d0
  real(kind=DP), parameter  :: ninth= 0.11111111111111111111d0

#else

  real(kind=DP), parameter  :: xnu = 15.75592d0, cc0 = 0.004235d0
  real(kind=DP) :: bet, delt

  real(kind=DP), parameter  :: thrd = 0.333333333333d0, sixth7 = 1.1666666666666666d0
  real(kind=DP), parameter  :: thrd7= 2.333333333333d0
  real(kind=DP), parameter  :: sixth13= 2.1666666666666666d0
  real(kind=DP), parameter  :: sixth19= 3.1666666666666666d0
  real(kind=DP), parameter  :: ninth= 0.111111111111111d0
#endif

  real(kind=DP), parameter  :: twntysvnth=  3.703703703703703d-2
  real(kind=DP), parameter  :: fftyfrth= 1.851851851851852d-2
  real(kind=DP), parameter  :: eityonth= 1.234567901234568d-2

  integer       :: is,i
  real(kind=DP) :: facw,exc1,g,g3,facpon,d,dd,rs,fk,sk,t,s,q0 &
       &         , rs12,rs32,rsp,q1,q2,eu,q3,eurs,ecd,pon &
       &         , b,b2,t2,t4,t6,q4,q5,h0,h,q8,h0t,h0b,h0rs,ht,hrs &
       &         , ec1,ec1d,ec1dd,exc0,excd,excdd &
       &         , thrd2,thrd4,fzdd,zeta,onpzeta,onmzeta,fz,zetadxd &
       &         , q0p,q1p,q2p,q3p,q0q,q1q,q2q,q3q,fzd,onzt4,gzt,bzt,h0zt,euzt &
       &         , zetadxda,zetadxdb,ecda,ecdb,ec1da,ec1db,excda,excdb
  real(kind=DP) :: rsdn,rsdnn,rsdnnn
  real(kind=DP) :: q1drr,q1drrr
  real(kind=DP) :: eudrr,eudrrr,invq12pq1,q12p1
  real(kind=DP) :: eudn,eudnn,eudnnn,ecdd,ecddd
  real(kind=DP) :: dtdn,dtdg,dtdnn,dtdgg,dtdng,dtdnnn,dtdnng,dtdngg,dtdggg
  real(kind=DP) :: dAde,dAdee,dAdeee
  real(kind=DP) :: dAdn,dAdnn,dAdnnn
  real(kind=DP) :: fE,fEdt,fEdtt,fEda,fEdaa
  real(kind=DP) :: fF,fFdt,fFdtt,fFda,fFdaa
  real(kind=DP) :: fG,fGdt,fGdtt,fGda,fGdaa
  real(kind=DP) :: at2,a2t4,a3t6,t3,t5,t7,t8
  real(kind=DP) :: dHda,dHdt,dHdaa,dHdat,dHdtt,dHdaaa,dHdaat,dHdatt,dHdttt,betivE,betivE2
  real(kind=DP) :: dHdn,dHdg,dHdnn,dHdng,dHdgg,dHdnnn,dHdnng,dHdngg,dHdggg
  real(kind=DP) :: zetada,zetadb,zetadaa,zetadab,zetadbb,zetadaaa,zetadaab,zetadabb,zetadbbb
  real(kind=DP) :: id2,id3,id4
  real(kind=DP) :: e1dr,e1drr,e1drrr
  real(kind=DP) :: e1pdr,e1pdrr,e1pdrrr
  real(kind=DP) :: e1qdr,e1qdrr,e1qdrrr
  real(kind=DP) :: q1pdrr,q1pdrrr,invq1p2pq1p,q1p2p1
  real(kind=DP) :: q1qdrr,q1qdrrr,invq1q2pq1q,q1q2p1
  real(kind=DP) :: fzdd0,fzddd,opzthrdm2,omzthrdm2
  real(kind=DP) :: zeta2,zeta3,zeta4,eudrz,eudzz,eudrrz,eudrzz,eudzzz
  real(kind=DP) :: e1,e2,e3
  real(kind=DP) :: euda,eudb,eudaa,eudbb,eudab,eudaaa,eudaab,eudabb,eudbbb
  real(kind=DP) :: rsdn2,rsdn3,zetada2,zetada3,zetadb2,zetadb3
  real(kind=DP) :: gzz,gzzz,duda,dudb,dudaa,dudab,dudbb,dudaaa,dudaab,dudabb,dudbbb
  real(kind=DP) :: dtdu,dtduu,dtduuu,dtdnu,dtdnnu,dtdnuu,dtdug,dtduug,dtdugg,dtdnug,invg
  real(kind=DP) :: dtda,dtdb,dtdaa,dtdab,dtdag,dtdbg,dtdbb
  real(kind=DP) :: dtdaaa,dtdaab,dtdabb,dtdbbb,dtdaag,dtdagg,dtdbbg,dtdbgg,dtdabg
  real(kind=DP) :: invg3,m3u
  real(kind=DP) :: dAdu,dAduu,dAduuu,dAdeu,dAdeeu,dAdeuu
  real(kind=DP) :: dAda,dAdb,dAdaa,dAdab,dAdbb,dAdaaa,dAdaab,dAdabb,dAdbbb
  real(kind=DP) :: dHdf,dHdff,dHdft,dHdfff,dHdfft,dHdftt
  real(kind=DP) :: dHdu,dHduu,dHduuu,dHdfu,dHdtu,dHdffu,dHdfuu,dHdttu,dHdtuu,dHdftu
  real(kind=DP) :: tiu,siu2
  real(kind=DP) :: dHdb,dHdbb,dHdag,dHdbg,dHdab,dHdbbb,dHdaab,dHdabb,dHdbbg,dHdbgg,dHdaag,dHdagg,dHdabg


#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: eta = 1.d-12
  real(kind=DP), parameter  :: sixthm = -0.1666666666666666666666d0
#endif

  facw = nspin

#ifndef GGA_PAW_REVISED
  bet = xnu*cc0
  delt = bet / gamma
#endif

     if ( nspin == 1 ) then
        g   = 1.d0
        g3  = g**3
!!$        g4  = g3*g
        facpon = -delt/(g3*bet)

        do i = ist,ien,dnr ! MPI
           d  = facw*chgrhr_l(i, 1)
           if(d < density_minimum) cycle
! gradient of charge density 95/12/3 H.S.
           dd = facw* grad_trho(i)
           rs = (0.75/(PAI*d))**thrd
           fk = (3*PAI*PAI*d)**thrd
           sk = dsqrt(4*fk/PAI)

           if(d > density_minimum2) then
              t = dd/(d*sk*2)
              s = dd/(d*fk*2)
           else
              t = 0.d0
              s = 0.d0
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
           dF_drho(i, 1) = dF_drho(i, 1) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(dabs(grad_trho(i)) > 1.0d-9) then
!              dF_dgradrho(i) = excdd / grad_trho(i)
              dF_dgradrho(i) = excdd
           else
              dF_dgradrho(i) = 0.d0
           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d/facw
!           exc1 = exc1 + exc0*wos(i)
           exc(i)=exc(i)+exc0
        
! ***** Calc heigher derivatives. *****
           
           rsdn=-thrd*rs/d
           rsdnn=-4.d0*thrd*rsdn/d
           rsdnnn=-thrd7*rsdnn/d
!print '(i2,5e19.10)',i,d,rs,rsdn,rsdnn,rsdnnn
!q3 = a*(b1/rs12+2*b2n+3*b3*rs12+2*b4*p1*rsp)   
           q1drr =0.5d0*a*( -b1/rs32+3.d0*b3/rs12+8.d0*b4 )
           q1drrr=0.75d0*a*( b1/rs - b3 )/rs32
!eurs = -2*a*a1*q2 - q0*q3/(q1**2 + q1)
           invq12pq1=1.d0/(q1**2 + q1)
           q12p1=2.d0*q1+1.d0
           eudrr=4.d0*a*a1*q3*invq12pq1-q0*( q1drr-q12p1*q3**2*invq12pq1 )*invq12pq1
           eudrrr=6.d0*a*a1*( q1drr - q12p1*q3**2*invq12pq1 ) * invq12pq1 &
                    -q0*( q1drrr -(2.d0*q3**3+3.d0*q12p1*q3*q1drr -  &
                    2.d0*q12p1**2*q3**3* invq12pq1)*invq12pq1 ) * invq12pq1
           eudn=eurs*rsdn
           eudnn=eudrr*rsdn**2+eurs*rsdnn
           eudnnn=eudrrr*rsdn**3+3.d0*eudrr*rsdn*rsdnn+eurs*rsdnnn
!print '(i2,5e19.12)',i,d,eu,eudn,eudnn,eudnnn
           ecdd=2.d0*eudn+d*eudnn
           ecddd=3.d0*eudnn+d*eudnnn
!print '(i2,5e19.12)',i,d,eu*d,ecdd,ecddd
           
           dtdn  = -sixth7*t/d
           dtdnn = -sixth13*dtdn/d
           dtdnnn= -sixth19*dtdnn/d
!print '(i2,5e19.12)',i,d,t,dtdn,dtdnn,dtdnnn
           if(dabs(grad_trho(i)) > 1.0d-9) then
               dtdg  = t/dd
               dtdng = -sixth7*dtdg/d
               dtdnng= -sixth13*dtdng/d
           else
               dtdg  = 0.d0
               dtdng = 0.d0
               dtdnng= 0.d0
           endif
           dtdngg = 0.d0
           dtdggg = 0.d0
           
           dAde   = b * ( b + delt ) / bet
           dAdee  = 2.d0 * dAde * ( b + 0.5d0 * delt ) / bet
           dAdeee = 2.d0 * ( dAde**2 + ( b + 0.5d0 * delt )*dAdee ) /bet
           
           dAdn   = dAde  *eudn
           dAdnn  = dAdee *eudn**2  + dAde*eudnn
           dAdnnn = dAdeee*eudn**3  + 3.d0*dAdee*eudn*eudnn + dAde*eudnnn   
!print '(i2,5e19.10)' ,i,d,b,dAdnn
           
           at2  = b * t2
           a2t4 = at2**2
           a3t6 = at2*a2t4
           t3   = t2*t
           t5   = t4*t
           t7   = t6*t
           t8   = t7*t
           fE    = q5**2 + delt*t2*q4*q5
           fEdt  = 4.d0*b*t*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt *2.d0*t*(1.d0 + 4.d0*at2 + 6.d0*a2t4 + 4.d0*a3t6 )
           fEdtt = 4.d0*b*( 1.d0 + 9.d0*at2 + 15.d0*a2t4 + 14.d0*a3t6 ) + &
                    delt*( 2.d0 + 24.d0*at2 + 60.d0*a2t4 + 56.d0*a3t6 )
           fEda  = 2.d0*t2*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt*t4*( 2.d0 + 4.d0*at2 + 3.d0*a2t4 )
           fEdaa = 6.d0*t4*( 2.d0*q5 - 1.d0) + &
                    delt*2.d0*t6*( 2.d0 + 3.d0*at2 )
           fF    = t*( 1.d0 + 2.d0*at2 )         
           fFdt  = 1.d0 + 6.d0*at2
           fFdtt = 12.d0*b*t
           fFda  = 2.d0*t3
           fFdaa = 0.d0
           fG    = 2.d0*b*t6 + b2*t8
           fGdt  = 12.d0*b*t5 + 8.d0*b2*t7
           fGdtt = 60.d0*b*t4 + 56.d0*b2*t6
           fGda  = 2.d0*t6 + 2.d0*b*t8
           fGdaa = 2.d0*t8
           
           betivE=bet/fE
           betivE2=betivE/fE
           dHda   =     -fG*betivE
           dHdt   = 2.d0*fF*betivE
           dHdaa  =     -( fGda*fE-fG*fEda )*betivE2
           dHdat  =     -( fGdt*fE-fG*fEdt )*betivE2
           dHdtt  = 2.d0*( fFdt*fE-fF*fEdt )*betivE2
           dHdaaa =     -( fGdaa*fE-fG*fEdaa-2.d0*( fGda*fE-fG*fEda )*fEda/fE )*betivE2
           dHdaat = 2.d0*( fFdaa*fE-fF*fEdaa-2.d0*( fFda*fE-fF*fEda )*fEda/fE )*betivE2
           dHdatt =     -( fGdtt*fE-fG*fEdtt-2.d0*( fGdt*fE-fG*fEdt )*fEdt/fE )*betivE2
           dHdttt = 2.d0*( fFdtt*fE-fF*fEdtt-2.d0*( fFdt*fE-fF*fEdt )*fEdt/fE )*betivE2
           
           dHdn     = dHda  *dAdn       + dHdt*dtdn
           dHdg     = dHdt  *dtdg
           dHdnn    = dHdaa *dAdn**2    + 2.d0*dHdat*dAdn*dtdn  + dHdtt*dtdn**2   + dHda*dAdnn  + dHdt*dtdnn
           dHdng    = dHdat *dAdn*dtdg  + dHdtt*dtdn*dtdg       + dHdt*dtdng
           dHdgg    = dHdtt*dtdg**2
           dHdnnn   = dHdaaa*dAdn**3    + dHdttt*dtdn**3 &
                        + 3.d0*dHdaat*dAdn**2*dtdn  + 3.d0*dHdatt*dAdn*dtdn**2 &
                        + 3.d0*dHdaa*dAdn*dAdnn + 3.d0*dHdtt*dtdn*dtdnn + 3.d0*dHdat*(dAdnn*dtdn + dAdn*dtdnn) &
                        + dHda*dAdnnn   + dHdt*dtdnnn
           dHdnng   = dHdaat*dAdn**2*dtdg   + 2.d0*dHdatt*dAdn*dtdn*dtdg    + dHdttt*dtdn**2*dtdg &
                        + dHdat*(dAdnn*dtdg + 2.d0*dAdn*dtdng)  + dHdtt*(dtdnn*dtdg + 2.d0*dtdn*dtdng) &
                        + dHdt*dtdnng
           dHdngg   = dHdatt*dAdn*dtdg**2   + dHdttt*dtdn*dtdg**2   + 2.d0*dHdtt*dtdg*dtdng
           dHdggg   = dHdttt*dtdg**3
!dF_drho(i, 1)=h+d*dHdn+ecd
           dFc_daa(i)   = 2.d0*dHdn + d*dHdnn       + ecdd
           dFc_dag(i)   = dHdg  + d*dHdng
           dFc_dgg(i)   = d*dHdgg
#ifdef USE_3RD_ORDER_DERIVATIVE
           dFc_daaa(i)  = 3.d0*dHdnn    + d*dHdnnn  + ecddd
           dFc_daag(i)  = 2.d0*dHdng    + d*dHdnng
           dFc_dagg(i)  = dHdgg         + d*dHdngg
           dFc_dggg(i)  = d*dHdggg  
#endif
        end do
     else if ( nspin ==  2 ) then
        thrd2 = thrd * 2
        thrd4 = thrd * 4
        fzdd = 8/(9*(2**thrd4 - 2 ))

        do i = ist,ien,dnr
           d = chgrhr_l(i, 1) + chgrhr_l(i, nspin)
           if(d < density_minimum .or. chgrhr_l(i,1) < density_minimum .or. chgrhr_l(i,nspin)<density_minimum) cycle
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
           ! ****** aas mkatsu *****
!              t = dd/(d*sk*2)
              t = dd/(d*sk*2)/g         !    mkatsu aas
              s = dd/(d*fk*2)
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

!           zetadxd = - 2 * chgrhr_l(i, 3-is) * (-1.d0)**is / d
           zetadxda = - 2 * chgrhr_l(i, 2) * (-1.d0) / d
           zetadxdb = - 2 * chgrhr_l(i, 1) / d
!           ecd = eu-thrd*rs*eurs + euzt * zetadxd
           ecda = eu-thrd*rs*eurs + euzt * zetadxda 
           ecdb = eu-thrd*rs*eurs + euzt * zetadxdb

#ifndef GGA_PAW_REVISED
           if(onpzeta < zeta_minimum2 .or. onmzeta < zeta_minimum2) then
              ec1  = 0.d0;    ec1d = 0.d0;     ec1dd= 0.d0
              ec1da = 0.d0;     ec1db = 0.d0
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
#ifdef GGA_PAW_REVISED
              gzt = ( (onpzeta**2+eta)**sixthm - (onmzeta**2+eta)**sixthm)/3.d0
#else
              gzt = ( 1.d0/onpzeta**thrd -1.d0/onmzeta**thrd )/3.d0
#endif
           end if

           bzt = b * ( b+delt) * ( euzt - 3*eu*gzt/g ) / ( bet * g3 ) 
           h0zt = 3 * gzt * h0 / g + h0b * bzt

           ec1  = d*h * 0.5
! --> revised by Kato-san,  1999Jul16
!!$           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + h0zt * zetadxd
!           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxd
           ec1da = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxda
           ec1db = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxdb
           ec1dd = 0.5*ht/sk / g

10         continue
           exc0 = eu*d * 0.5 + ec1
!           excd = ecd + ec1d
           excda= ecda + ec1da
           excdb= ecdb + ec1db
!           dF_drho(i, is) = dF_drho(i, is) + excd
           dF_drho(i, 1) = dF_drho(i, 1) + excda
           dF_drho(i, 2) = dF_drho(i, 2) + excdb
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
!           if(is ==  2) then
              if(dabs(grad_trho(i)) > 1.0d-9) then
!                 dF_dgradrho(i) = excdd / grad_trho(i)
                 dF_dgradrho(i) = excdd
              else
                 dF_dgradrho(i) = 0.d0
              endif
!           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d * 0.5
!           exc1 = exc1 + exc0*wos(i)
           exc(i)=exc(i)+exc0*2.d0
           
! ***** Calc heigher order derivatives. *****

           rsdn=-thrd*rs/d
           rsdnn=-4.d0*thrd*rsdn/d
           rsdnnn=-thrd7*rsdnn/d
           
           id2=1.d0/d/d
           id3=id2/d
           id4=id3/d
           zetada   =     zetadxda/d
           zetadb   =     zetadxdb/d
           zetadaa  =    -4.d0*chgrhr_l(i, 2)*id3
           zetadab  =     2.d0*zeta*id2
           zetadbb  =     4.d0*chgrhr_l(i, 1)*id3
           zetadaaa =    12.d0*chgrhr_l(i, 2)*id4
           zetadaab =   ( -4.d0*chgrhr_l(i, 1) + 8.d0*chgrhr_l(i, 2))*id4
           zetadabb =   ( -8.d0*chgrhr_l(i, 1) + 4.d0*chgrhr_l(i, 2))*id4
           zetadbbb =   -12.d0*chgrhr_l(i, 1)*id4
           
!print '(i2,7d19.10)', i,chgrhr_l(i,1:2),zetadab,zetadaaa,zetadaab,zetadabb,zetadbbb
             
           q1drr =0.5d0*a*( -b1/rs32+3.d0*b3/rs12+8.d0*b4 )
           q1drrr=0.75d0*a*( b1/rs - b3 )/rs32
           invq12pq1=1.d0/(q1**2 + q1)
           q12p1=2.d0*q1+1.d0
           e1dr=-2.d0*a*a1*q2-q0*q3*invq12pq1
           e1drr=4.d0*a*a1*q3*invq12pq1-q0*( q1drr-q12p1*q3**2*invq12pq1 )*invq12pq1
           e1drrr=6.d0*a*a1*( q1drr - q12p1*q3**2*invq12pq1 ) * invq12pq1 &
                    -q0*( q1drrr -(2.d0*q3**3+3.d0*q12p1*q3*q1drr -  &
                    2.d0*q12p1**2*q3**3* invq12pq1)*invq12pq1 ) * invq12pq1
                    
           q1pdrr =0.5d0*ap*( -b1p/rs32+3.d0*b3p/rs12+8.d0*b4p )
           q1pdrrr=0.75d0*ap*( b1p/rs - b3p )/rs32
           invq1p2pq1p=1.d0/(q1p**2 + q1p)
           q1p2p1=2.d0*q1p+1.d0
           e1pdr=-2.d0*ap*a1p*q2p-q0p*q3p*invq1p2pq1p
           e1pdrr=4.d0*ap*a1p*q3p*invq1p2pq1p-q0p*( q1pdrr-q1p2p1*q3p**2*invq1p2pq1p )*invq1p2pq1p
           e1pdrrr=6.d0*ap*a1p*( q1pdrr - q1p2p1*q3p**2*invq1p2pq1p ) * invq1p2pq1p &
                    -q0p*( q1pdrrr -(2.d0*q3p**3+3.d0*q1p2p1*q3p*q1pdrr -  &
                    2.d0*q1p2p1**2*q3p**3* invq1p2pq1p)*invq1p2pq1p ) * invq1p2pq1p
                    
           q1qdrr =0.5d0*aq*( -b1q/rs32+3.d0*b3q/rs12+8.d0*b4q )
           q1qdrrr=0.75d0*aq*( b1q/rs - b3q )/rs32
           invq1q2pq1q=1.d0/(q1q**2 + q1q)
           q1q2p1=2.d0*q1q+1.d0
           e1qdr=-2.d0*aq*a1q*q2q-q0q*q3q*invq1q2pq1q
           e1qdrr=4.d0*aq*a1q*q3q*invq1q2pq1q-q0q*( q1qdrr-q1q2p1*q3q**2*invq1q2pq1q )*invq1q2pq1q
           e1qdrrr=6.d0*aq*a1q*( q1qdrr - q1q2p1*q3q**2*invq1q2pq1q ) * invq1q2pq1q &
                    -q0q*( q1qdrrr -(2.d0*q3q**3+3.d0*q1q2p1*q3q*q1qdrr -  &
                    2.d0*q1q2p1**2*q3q**3* invq1q2pq1q)*invq1q2pq1q ) * invq1q2pq1q
           
           zeta2=zeta**2
           zeta3=zeta2*zeta
           zeta4=zeta2*zeta2
           
           if(dabs(zeta) < zeta_minimum) then
              fzdd0 = fzdd &
                   &            *( 1 + 5.d0*ninth *zeta2 &
                   &            *( 1 +22.d0*twntysvnth *zeta2 ))
              fzddd = fzdd*10.d0*ninth*zeta*( 1.d0 + 44.d0*twntysvnth*zeta2 )
              gzz   = -2.d0*ninth*( 1.d0 + 14.d0*ninth*zeta2*( 1.d0 + 65.d0*fftyfrth*zeta2 ) )
              gzzz  = -56.d0*eityonth*zeta*( 1.d0 + 65.d0*twntysvnth*zeta2 )
           else
              opzthrdm2=1.d0/onpzeta**(thrd2)
              omzthrdm2=1.d0/onmzeta**(thrd2)
              fzdd0 =  fzdd*0.5d0*( opzthrdm2 + omzthrdm2 )
              fzddd = -fzdd*thrd*( opzthrdm2/onpzeta - omzthrdm2/onmzeta )
              gzz   = -ninth*( opzthrdm2**2 + omzthrdm2**2 )
!              gzz   = -ninth*( 1.d0/onpzeta**(thrd4) + 1.d0/onmzeta**(thrd4) )
              gzzz  = 4.d0*twntysvnth*( opzthrdm2**2/onpzeta - omzthrdm2**2/onmzeta )
           end if
           
           e1=q0*q2
           e2=q0p*q2p
           e3=q0q*q2q
           
           eudrr    = e1drr  + ( e1pdrr  -  e1drr  ) * fz * zeta4 - e1qdrr  * fz * onzt4 / fzdd
           eudrrr   = e1drrr + ( e1pdrrr -  e1drrr ) * fz * zeta4 - e1qdrrr * fz * onzt4 / fzdd
           eudzz    = ( e2 - e1 )*( fzdd0*zeta4 + 8.d0*fzd*zeta3 + 12.d0*fz*zeta2 ) &
                        - e3*( fzdd0*onzt4 - 8.d0*fzd*zeta3 - 12.d0*fz*zeta2 )/fzdd
           eudzzz   = ( e2 - e1 )*( fzddd*zeta4 + 12.d0*fzdd0*zeta3 + 36.d0*fzd*zeta2 + 24.d0*fz*zeta ) &
                        - e3*( fzddd*onzt4 - 12.d0*fzdd0*zeta3 - 36.d0*fzd*zeta2 - 24.d0*fz*zeta )/fzdd
           eudrz    = ( e1pdr - e1dr )*( fzd*zeta4 + 4.d0*fz*zeta3 ) &
                        - e1qdr *( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd 
           eudrrz   = ( e1pdrr  - e1drr  )*( fzd*zeta4 + 4.d0*fz*zeta3 ) &
                        - e1qdrr *( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd
           eudrzz   = ( e1pdr - e1dr )*( fzdd0*zeta4 + 8.d0*fzd*zeta3 + 12.d0*fz*zeta2 ) &
                        - e1qdr *( fzdd0*onzt4 - 8.d0*fzd*zeta3 - 12.d0*fz*zeta2 )/fzdd
                        
           rsdn2    = rsdn*rsdn
           rsdn3    = rsdn2*rsdn
           zetada2  = zetada**2
           zetada3  = zetada2*zetada
           zetadb2  = zetadb**2
           zetadb3  = zetadb2*zetadb
                        
           euda     = eurs*rsdn + euzt*zetada
           eudb     = eurs*rsdn + euzt*zetadb
           eudaa    = eudrr*rsdn2 + 2.d0*eudrz*rsdn*zetada + eudzz*zetada2 + eurs*rsdnn + euzt*zetadaa
           eudab    = eudrr*rsdn2 + eudrz*rsdn*( zetada + zetadb ) + eudzz*zetada*zetadb &
                        + eurs*rsdnn + euzt*zetadab
           eudbb    = eudrr*rsdn2 + 2.d0*eudrz*rsdn*zetadb + eudzz*zetadb2 + eurs*rsdnn + euzt*zetadbb
           eudaaa   = eudrrr*rsdn3 &
                        + 3.d0*eudrrz*rsdn2*zetada &
                        + 3.d0*eudrzz*rsdn*zetada2 + eudzzz*zetada3 &
                        + 3.d0*eudrr*rsdn*rsdnn + 3.d0*eudzz*zetada*zetadaa &
                        + 3.d0*eudrz*( rsdnn*zetada + rsdn*zetadaa ) &
                        + eurs*rsdnnn + euzt*zetadaaa
           eudaab   = eudrrr*rsdn3 &
                        + eudrrz*rsdn2*( zetadb + 2.d0*zetada ) &
                        + eudrzz*rsdn*zetada*( zetada + 2.d0*zetadb ) &
                        + eudzzz*zetada2*zetadb &
                        + 3.d0*eudrr*rsdn*rsdnn &
                        + eudrz*( rsdnn*( 2.d0*zetada + zetadb ) + rsdn*( 2.d0*zetadab + zetadaa ) ) &
                        + eudzz*( 2.d0*zetada*zetadab + zetadaa*zetadb ) &
                        + eurs*rsdnnn + euzt*zetadaab
           eudabb   = eudrrr*rsdn3 &
                        + eudrrz*rsdn2*( zetada + 2.d0*zetadb ) &
                        + eudrzz*rsdn*zetadb*( zetadb + 2.d0*zetada ) &
                        + eudzzz*zetada*zetadb2 &
                        + 3.d0*eudrr*rsdn*rsdnn &
                        + eudrz*( rsdnn*( zetada + 2.d0*zetadb ) + rsdn*( 2.d0*zetadab + zetadbb ) ) &
                        + eudzz*( zetada*zetadbb + 2.d0*zetadb*zetadab ) &
                        + eurs*rsdnnn + euzt*zetadabb
           eudbbb   = eudrrr*rsdn3 &
                        + 3.d0*eudrrz*rsdn2*zetadb &
                        + 3.d0*eudrzz*rsdn*zetadb2 + eudzzz*zetadb3 &
                        + 3.d0*eudrr*rsdn*rsdnn + 3.d0*eudzz*zetadb*zetadbb &
                        + 3.d0*eudrz*( rsdnn*zetadb + rsdn*zetadbb ) &
                        + eurs*rsdnnn + euzt*zetadbbb
                        
!print '(i2,12d19.10)', i,chgrhr_l(i,1:2),eu,euda,eudb,eudaa,eudbb,eudab,eudaaa,eudaab,eudabb,eudbbb
!print *,eu,     e1+(e2-e1)* fz * zeta4 - e3  * fz * onzt4 / fzdd
!print *,eurs,   e1dr  + ( e1pdr  -  e1dr  ) * fz * zeta4 - e1qdr  * fz * onzt4 / fzdd
!print *,euzt,   (e2-e1)*( fzd*zeta4 + 4.d0*fz*zeta3 )-e3*( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd 

           duda     = gzt*zetada
           dudb     = gzt*zetadb 
           dudaa    = gzz*zetada2     + gzt*zetadaa
           dudab    = gzz*zetada*zetadb + gzt*zetadab
           dudbb    = gzz*zetadb2     + gzt*zetadbb
           dudaaa   = gzzz*zetada3    + 3.d0*gzz*zetada*zetadaa   + gzt*zetadaaa
           dudaab   = gzzz*zetada2*zetadb + gzz*( zetadaa*zetadb + 2.d0*zetada*zetadab ) + gzt*zetadaab
           dudabb   = gzzz*zetada*zetadb2 + gzz*( zetada*zetadbb + 2.d0*zetadab*zetadb ) + gzt*zetadabb
           dudbbb   = gzzz*zetadb3    + 3.d0*gzz*zetadb*zetadbb   + gzt*zetadbbb 
           
!print '(i2,12d19.10)', i,chgrhr_l(i,1:2),g,duda,dudb,dudaa,dudbb,dudab,dudaaa,dudaab,dudabb,dudbbb


           
           dtdn  = -sixth7*t/d
           dtdnn = -sixth13*dtdn/d
           dtdnnn= -sixth19*dtdnn/d
!print '(i2,5e19.12)',i,d,t,dtdn,dtdnn,dtdnnn
           if(dabs(grad_trho(i)) > 1.0d-9) then
               dtdg  = t/dd
               dtdng = -sixth7*dtdg/d
               dtdnng= -sixth13*dtdng/d
!               dtdnug = -dtdn*invg/dd
           else
               dtdg  = 0.d0
               dtdng = 0.d0
               dtdnng= 0.d0
               dtdnug= 0.d0 
           endif
           dtdgg  = 0.d0
           dtdngg = 0.d0
           dtdggg = 0.d0
           invg=1.d0/g
           dtdu   = -t*invg
           dtduu  = -2.d0*dtdu*invg 
           dtduuu = -3.d0*dtduu*invg
           dtdnu  = -dtdn*invg
           dtdnnu = -dtdnn*invg
           dtdnuu = -sixth7*dtduu/d
           dtdug  = -dtdg*invg
           dtduug = 2.d0*dtdg*invg**2
           dtdugg = 0.d0
           if(dabs(grad_trho(i)) > 1.0d-9) then
                dtdnug = dtdnu/dd
           else
                dtdnug = 0.d0
           end if
           
           dtda   = dtdn + dtdu*duda
           dtdb   = dtdn + dtdu*dudb
           dtdaa  = dtdnn + 2.d0*dtdnu*duda + dtduu*duda**2 + dtdu*dudaa
           dtdbb  = dtdnn + 2.d0*dtdnu*dudb + dtduu*dudb**2 + dtdu*dudbb
           dtdab  = dtdnn + dtdnu*( duda + dudb ) + dtduu*duda*dudb + dtdu*dudab
           dtdag  = dtdng + dtdug*duda
           dtdbg  = dtdng + dtdug*dudb
           dtdaaa = dtdnnn+ 3.d0*dtdnnu*duda + 3.d0*dtdnuu*duda**2 + dtduuu*duda**3 &
                    + 3.d0*dtduu*duda*dudaa + 3.d0*dtdnu*dudaa + dtdu*dudaaa
           dtdbbb = dtdnnn+ 3.d0*dtdnnu*dudb + 3.d0*dtdnuu*dudb**2 + dtduuu*dudb**3 &
                    + 3.d0*dtduu*dudb*dudbb + 3.d0*dtdnu*dudbb + dtdu*dudbbb
           dtdaab = dtdnnn+ dtdnnu*( 2.d0*duda + dudb ) + dtdnuu*duda*( duda + 2.d0*dudb ) + dtduuu*duda**2*dudb &
                    + dtduu*( dudaa*dudb + 2.d0*duda*dudab ) + dtdnu*( dudaa + 2.d0*dudab ) + dtdu*dudaab
           dtdabb = dtdnnn+ dtdnnu*( duda + 2.d0*dudb ) + dtdnuu*dudb*( dudb + 2.d0*duda ) + dtduuu*dudb**2*duda &
                    + dtduu*( dudbb*duda + 2.d0*dudb*dudab ) + dtdnu*( dudbb + 2.d0*dudab ) + dtdu*dudabb 
           dtdaag = dtdnng+ 2.d0*dtdnug*duda + dtduug*duda**2 + dtdug*dudaa
           dtdbbg = dtdnng+ 2.d0*dtdnug*dudb + dtduug*dudb**2 + dtdug*dudbb
           dtdabg = dtdnng+ dtdnug*( duda + dudb ) + dtduug*duda*dudb + dtdug*dudab 
           
!print '(i2,20d19.10)', i,t,dtda,dtdb,dtdg,dtdaa,dtdbb,dtdgg,dtdab,dtdag,dtdbg,dtdaaa,dtdbbb,dtdggg,dtdaab,dtdabb,dtdaag,dtdagg,dtdbbg,dtdbgg,dtdabg     
!print '(i2,12d19.10)', i,t,dtdag,dtdbg,dtdaag,dtdbbg,dtdabg     

           invg3  = 1.d0  /g3
           dAde   = b * ( b + delt ) / bet * invg3
           dAdee  = 2.d0 * dAde * ( b + 0.5d0 * delt ) / bet *invg3
           dAdeee = 2.d0 * ( dAde**2 + ( b + 0.5d0 * delt )*dAdee ) /bet * invg3
           m3u    = -3.d0*invg
           dAdu     = m3u*eu*dAde
           dAdeu    = m3u*( dAde + eu*dAdee )
           dAduu    = m3u*eu*invg*( -dAde + g*dAdeu )
           dAdeeu   = m3u*( 2.d0*dAdee + eu*dAdeee )
           dAdeuu   = m3u*invg*( -dAde + g*dAdeu - eu*dAdee + g*eu*dAdeeu )
           dAduuu   = m3u*eu*invg**2*( 2.d0*dAde -2.d0*g*dAdeu + g**2*dAdeuu )
           
           dAda     = dAde*euda + dAdu*duda
           dAdb     = dAde*eudb + dAdu*dudb
           dAdaa    = dAdee*euda**2 + 2.d0*dAdeu*euda*duda + dAduu*duda**2 + dAde*eudaa + dAdu*dudaa 
           dAdab    = dAdee*euda*eudb + dAdeu*( euda*dudb + eudb*duda ) + dAduu*duda*dudb + dAde*eudab + dAdu*dudab 
           dAdbb    = dAdee*eudb**2 + 2.d0*dAdeu*eudb*dudb + dAduu*dudb**2 + dAde*eudbb + dAdu*dudbb
           dAdaaa   = dAdeee*euda**3 + 3.d0*dAdeeu*euda**2*duda + 3.d0*dAdeuu*euda*duda**2 + dAduuu*duda**3 &
                        + 3.d0*dAdee*euda*eudaa + 3.d0*dAduu*duda*dudaa &
                        + 3.d0*dAdeu*( euda*dudaa + eudaa*duda ) + dAde*eudaaa + dAdu*dudaaa
           dAdaab   = dAdeee*euda**2*eudb + dAduuu*duda**2*dudb &
                        + dAdeeu*euda*( euda*dudb + 2.d0*eudb*duda ) &
                        + dAdeuu*( 2.d0*euda*dudb + eudb*duda )*duda &
                        + dAdee*( 2.d0*euda*eudab + eudb*eudaa ) + dAduu*( 2.d0*duda*dudab + dudb*dudaa ) &
                        + dAdeu*( 2.d0*eudab*duda + 2.d0*euda*dudab + dudb*eudaa + eudb*dudaa ) &
                        + dAde*eudaab + dAdu*dudaab
           dAdabb   = dAdeee*euda*eudb**2 + dAduuu*dudb**2*duda &
                        + dAdeeu*eudb*( eudb*duda + 2.d0*euda*dudb ) &
                        + dAdeuu*( 2.d0*eudb*duda + euda*dudb )*dudb &
                        + dAdee*( eudbb*euda + 2.d0*eudb*eudab ) + dAduu*( dudbb*duda + 2.d0*dudb*dudab ) &
                        + dAdeu*( eudbb*duda + 2.d0*eudb*dudab + 2.d0*dudb*eudab + dudbb*euda ) &
                        + dAde*eudabb + dAdu*dudabb
           dAdbbb   = dAdeee*eudb**3 + 3.d0*dAdeeu*eudb**2*dudb + 3.d0*dAdeuu*eudb*dudb**2 + dAduuu*dudb**3 &
                        + 3.d0*dAdee*eudb*eudbb + 3.d0*dAduu*dudb*dudbb &
                        + 3.d0*dAdeu*( eudb*dudbb + eudbb*dudb ) + dAde*eudbbb + dAdu*dudbbb 
                        
!print '(i2,12d19.10)', i,b,dAda,dAdb,dAdaa,dAdab,dAdbb,dAdaaa,dAdaab,dAdabb,dAdbbb

           at2  = b * t2
           a2t4 = at2**2
           a3t6 = at2*a2t4
           t3   = t2*t
           t5   = t4*t
           t7   = t6*t
           t8   = t7*t
           fE    = q5**2 + delt*t2*q4*q5
           fEdt  = 4.d0*b*t*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt *2.d0*t*(1.d0 + 4.d0*at2 + 6.d0*a2t4 + 4.d0*a3t6 )
           fEdtt = 4.d0*b*( 1.d0 + 9.d0*at2 + 15.d0*a2t4 + 14.d0*a3t6 ) + &
                    delt*( 2.d0 + 24.d0*at2 + 60.d0*a2t4 + 56.d0*a3t6 )
           fEda  = 2.d0*t2*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt*t4*( 2.d0 + 4.d0*at2 + 3.d0*a2t4 )
           fEdaa = 6.d0*t4*( 2.d0*q5 - 1.d0) + &
                    delt*2.d0*t6*( 2.d0 + 3.d0*at2 )
           fF    = t*( 1.d0 + 2.d0*at2 )         
           fFdt  = 1.d0 + 6.d0*at2
           fFdtt = 12.d0*b*t
           fFda  = 2.d0*t3
           fFdaa = 0.d0
           fG    = 2.d0*b*t6 + b2*t8
           fGdt  = 12.d0*b*t5 + 8.d0*b2*t7
           fGdtt = 60.d0*b*t4 + 56.d0*b2*t6
           fGda  = 2.d0*t6 + 2.d0*b*t8
           fGdaa = 2.d0*t8
           
           betivE=bet/fE*g3
           betivE2=betivE/fE
           dHdf   =     -fG*betivE
           dHdt   = 2.d0*fF*betivE
           dHdff  =     -( fGda*fE-fG*fEda )*betivE2
           dHdft  =     -( fGdt*fE-fG*fEdt )*betivE2
           dHdtt  = 2.d0*( fFdt*fE-fF*fEdt )*betivE2
           dHdfff =     -( fGdaa*fE-fG*fEdaa-2.d0*( fGda*fE-fG*fEda )*fEda/fE )*betivE2
           dHdfft = 2.d0*( fFdaa*fE-fF*fEdaa-2.d0*( fFda*fE-fF*fEda )*fEda/fE )*betivE2
           dHdftt =     -( fGdtt*fE-fG*fEdtt-2.d0*( fGdt*fE-fG*fEdt )*fEdt/fE )*betivE2
           dHdttt = 2.d0*( fFdtt*fE-fF*fEdtt-2.d0*( fFdt*fE-fF*fEdt )*fEdt/fE )*betivE2
           
           tiu    = 3.d0*invg
           siu2   = tiu*2.d0*invg
           dHdu   = tiu*h
           dHduu  = siu2*h
           dHduuu = siu2*invg*h
           dHdfu  = tiu*dHdf
           dHdtu  = tiu*dHdt
           dHdffu = tiu*dHdff
           dHdfuu = siu2*dHdf
           dHdttu = tiu*dHdtt
           dHdtuu = siu2*dHdt
           dHdftu = tiu*dHdft
           
           dHda   = dHdf*dAda + dHdt*dtda + dHdu*duda
           dHdb   = dHdf*dAdb + dHdt*dtdb + dHdu*dudb
           dHdg   =             dHdt*dtdg
           dHdaa  = dHdff*dAda**2 + dHdtt*dtda**2 + dHduu*duda**2 &
                    + 2.d0*dHdft*dAda*dtda + 2.d0*dHdfu*dAda*duda + 2.d0*dHdtu*dtda*duda &
                    + dHdf*dAdaa + dHdt*dtdaa + dHdu*dudaa
           dHdab  = dHdff*dAda*dAdb + dHdtt*dtda*dtdb + dHduu*duda*dudb &
                    + dHdft*( dAda*dtdb + dAdb*dtda ) + dHdfu*( dAda*dudb + dAdb*duda ) + dHdtu*( dtda*dudb + dtdb*duda ) &
                    + dHdf*dAdab + dHdt*dtdab + dHdu*dudab
           dHdbb  = dHdff*dAdb**2 + dHdtt*dtdb**2 + dHduu*dudb**2 &
                    + 2.d0*dHdft*dAdb*dtdb + 2.d0*dHdfu*dAdb*dudb + 2.d0*dHdtu*dtdb*dudb &
                    + dHdf*dAdbb + dHdt*dtdbb + dHdu*dudbb
           dHdgg  = dHdtt*dtdg**2
           dHdag  = dHdft*dAda*dtdg + dHdtt*dtda*dtdg + dHdtu*duda*dtdg + dHdt*dtdag
           dHdbg  = dHdft*dAdb*dtdg + dHdtt*dtdb*dtdg + dHdtu*dudb*dtdg + dHdt*dtdbg
           dHdaaa = dHdfff*dAda**3 + 3.d0*dHdfft*dAda**2*dtda + 3.d0*dHdffu*dAda**2*duda &
                    + 3.d0*dHdftt*dAda*dtda**2 + dHdttt*dtda**3 + 3.d0*dHdttu*dtda**2*duda &
                    + 3.d0*dHdfuu*dAda*duda**2 + 3.d0*dHdtuu*dtda*duda**2 + dHduuu*duda**3 &
                    + 6.d0*dHdftu*dAda*dtda*duda &
                    + 3.d0*dHdff*dAda*dAdaa + 3.d0*dHdtt*dtda*dtdaa + 3.d0*dHduu*duda*dudaa &
                    + 3.d0*dHdft*( dAdaa*dtda + dAda*dtdaa ) &
                    + 3.d0*dHdtu*( dtdaa*duda + dtda*dudaa ) &
                    + 3.d0*dHdfu*( dAdaa*duda + dAda*dudaa ) &
                    + dHdf*dAdaaa + dHdt*dtdaaa + dHdu*dudaaa
           dHdaab = dHdfff*dAda**2*dAdb + dHdfft*dAda*( dAda*dtdb + 2.d0*dAdb*dtda ) &
                    + dHdffu*dAda*( dAda*dudb + 2.d0*dAdb*duda ) &
                    + dHdftt*dtda*( dAdb*dtda + 2.d0*dAda*dtdb ) &
                    + dHdttt*dtdb*dtda**2 + dHdttu*dtda*( dtda*dudb + 2.d0*dtdb*duda ) &
                    + dHdfuu*duda*( dAdb*duda + 2.d0*dAda*dudb ) &
                    + dHdtuu*duda*( dtdb*duda + 2.d0*dtda*dudb ) + dHduuu*duda**2*dudb &
                    + 2.d0*dHdftu*( dAdb*duda*dtda + dAda*dudb*dtda + dAda*duda*dtdb ) &
                    + dHdff*( 2.d0*dAda*dAdab + dAdb*dAdaa ) &
                    + dHdtt*( 2.d0*dtda*dtdab + dtdb*dtdaa ) &
                    + dHduu*( 2.d0*duda*dudab + dudb*dudaa ) &
                    + dHdft*( 2.d0*dAdab*dtda + 2.d0*dAda*dtdab &
                    + dAdaa*dtdb + dAdb*dtdaa ) &
                    + dHdtu*( 2.d0*dtdab*duda + 2.d0*dtda*dudab &
                    + dtdaa*dudb + dtdb*dudaa ) &
                    + dHdfu*( 2.d0*dAdab*duda + 2.d0*dAda*dudab &
                    + dAdaa*dudb + dAdb*dudaa ) &
                    + dHdf*dAdaab + dHdt*dtdaab + dHdu*dudaab
           dHdabb = dHdfff*dAda*dAdb**2 + dHdfft*dAdb*( dAdb*dtda + 2.d0*dAda*dtdb ) &
                    + dHdffu*dAdb*( dAdb*duda + 2.d0*dAda*dudb ) &
                    + dHdftt*dtdb*( dAda*dtdb + 2.d0*dAdb*dtda ) + dHdttt*dtda*dtdb**2 &
                    + dHdttu*dtdb*( dtdb*duda + 2.d0*dtda*dudb ) &
                    + dHdfuu*dudb*( dAda*dudb + 2.d0*dAdb*duda ) &
                    + dHdtuu*dudb*( dtda*dudb + 2.d0*dtdb*duda ) + dHduuu*duda*dudb**2 &
                    + 2.d0*dHdftu*( dAda*dtdb*dudb + dAdb*dtda*dudb + dAdb*dtdb*duda ) &
                    + dHdff*( 2.d0*dAdb*dAdab + dAda*dAdbb ) &
                    + dHdtt*( 2.d0*dtdb*dtdab + dtda*dtdbb ) &
                    + dHduu*( 2.d0*dudb*dudab + duda*dudbb ) &
                    + dHdft*( 2.d0*dAdab*dtdb + 2.d0*dAdb*dtdab + dAdbb*dtda + dAda*dtdbb ) &
                    + dHdtu*( 2.d0*dtdab*dudb + 2.d0*dtdb*dudab + dtdbb*duda + dtda*dudbb ) &
                    + dHdfu*( 2.d0*dAdab*dudb + 2.d0*dAdb*dudab + dAdbb*duda + dAda*dudbb ) &
                    + dHdf*dAdabb + dHdt*dtdabb + dHdu*dudabb
           dHdbbb = dHdfff*dAdb**3 + 3.d0*dHdfft*dAdb**2*dtdb + 3.d0*dHdffu*dAdb**2*dudb &
                    + 3.d0*dHdftt*dAdb*dtdb**2 + dHdttt*dtdb**3 + 3.d0*dHdttu*dtdb**2*dudb &
                    + 3.d0*dHdfuu*dAdb*dudb**2 + 3.d0*dHdtuu*dtdb*dudb**2 + dHduuu*dudb**3 &
                    + 6.d0*dHdftu*dAdb*dtdb*dudb &
                    + 3.d0*dHdff*dAdb*dAdbb + 3.d0*dHdtt*dtdb*dtdbb + 3.d0*dHduu*dudb*dudbb &
                    + 3.d0*dHdft*( dAdbb*dtdb + dAdb*dtdbb ) &
                    + 3.d0*dHdtu*( dtdbb*dudb + dtdb*dudbb ) &
                    + 3.d0*dHdfu*( dAdbb*dudb + dAdb*dudbb ) &
                    + dHdf*dAdbbb + dHdt*dtdbbb + dHdu*dudbbb
           dHdabg = dHdfft*dAda*dAdb*dtdg   + dHdftt*( dAdb*dtda + dAda*dtdb )*dtdg &
                    + dHdttt*dtda*dtdb*dtdg + dHdttu*( dtdb*duda + dtda*dudb )*dtdg &
                    + dHdtuu*duda*dudb*dtdg + dHdftu*( dAdb*duda + dAda*dudb )*dtdg &
                    + dHdft*( dAdab*dtdg + dAdb*dtdag + dAda*dtdbg ) &
                    + dHdtt*( dtdab*dtdg + dtdb*dtdag + dtda*dtdbg ) &
                    + dHdtu*( dudab*dtdg + dudb*dtdag + duda*dtdbg ) &
                    + dHdt*dtdabg
           dHdaag = dHdfft*dAda**2*dtdg + dHdttt*dtda**2*dtdg + dHdtuu*duda**2*dtdg &
                    + 2.d0*dHdttu*dtda*duda*dtdg + 2.d0*dHdftt*dAda*dtda*dtdg + 2.d0*dHdftu*dAda*duda*dtdg &
                    + dHdtt*( dtdaa*dtdg + 2.d0*dtda*dtdag ) + dHdft*( dAdaa*dtdg + 2.d0*dAda*dtdag ) &
                    + dHdtu*( dudaa*dtdg + 2.d0*duda*dtdag ) + dHdt*dtdaag         
           dHdbbg = dHdfft*dAdb**2*dtdg + dHdttt*dtdb**2*dtdg + dHdtuu*dudb**2*dtdg &
                    + 2.d0*dHdttu*dtdb*dudb*dtdg + 2.d0*dHdftt*dAdb*dtdb*dtdg + 2.d0*dHdftu*dAdb*dudb*dtdg &
                    + dHdtt*( dtdbb*dtdg + 2.d0*dtdb*dtdbg ) + dHdft*( dAdbb*dtdg + 2.d0*dAdb*dtdbg ) &
                    + dHdtu*( dudbb*dtdg + 2.d0*dudb*dtdbg ) + dHdt*dtdbbg       
           dHdagg = dHdftt*dAda*dtdg**2 + dHdttt*dtda*dtdg**2 + dHdttu*duda*dtdg**2 + 2.d0*dHdtt*dtdg*dtdag
           dHdbgg = dHdftt*dAdb*dtdg**2 + dHdttt*dtdb*dtdg**2 + dHdttu*dudb*dtdg**2 + 2.d0*dHdtt*dtdg*dtdbg 
           dHdggg = dHdttt*dtdg**3        
                    
!print '(i2,20d19.10)', i,h,dHda,dHdb,dHdg,dHdaa,dHdbb,dHdgg,dHdag,dHdbg,dHdab,dHdaaa,dHdbbb,dHdggg,dHdaab,dHdabb,dHdaag,dHdagg,dHdbbg,dHdbgg,dHdabg
    
           dFc_daa(i)   = 2.d0*euda + d*eudaa + 2.d0*dHda + d*dHdaa 
           dFc_dbb(i)   = 2.d0*eudb + d*eudbb + 2.d0*dHdb + d*dHdbb
           dFc_dgg(i)   = d*dHdgg
           dFc_dab(i)   = euda + eudb + d*eudab + dHda + dHdb + d*dHdab 
           dFc_dag(i)   = dHdg + d*dHdag
           dFc_dbg(i)   = dHdg + d*dHdbg
#ifdef USE_3RD_ORDER_DERIVATIVE
           dFc_daaa(i)  = 3.d0*eudaa + d*eudaaa + 3.d0*dHdaa + d*dHdaaa
           dFc_dbbb(i)  = 3.d0*eudbb + d*eudbbb + 3.d0*dHdbb + d*dHdbbb
           dFc_dggg(i)  = d*dHdggg
           dFc_daab(i)  = eudaa + 2.d0*eudab +d*eudaab + dHdaa + 2.d0*dHdab + d*dHdaab
           dFc_daag(i)  = 2.d0*dHdag + d*dHdaag
           dFc_dabb(i)  = eudbb + 2.d0*eudab + d*eudabb + dHdbb + 2.d0*dHdab + d*dHdabb
           dFc_dbbg(i)  = 2.d0*dHdbg + d*dHdbbg
           dFc_dagg(i)  = dHdgg + d*dHdagg
           dFc_dbgg(i)  = dHdgg + d*dHdbgg
           dFc_dabg(i)  = dHdag + dHdbg + d*dHdabg
#endif                    
!print '(i2,21d19.10)', i,exc(i),d*(eu+h),dF_drho(i,1),dF_drho(i,2),dF_dgradrho(i)*grad_trho(i) &
!                      ,dFc_daa(i),dFc_dbb(i),dFc_dgg(i),dFc_dab(i),dFc_dag(i),dFc_dbg(i),dFc_daaa(i),dFc_dbbb(i) &
!                      ,dFc_dggg(i),dFc_daab(i),dFc_daag(i),dFc_dabb(i),dFc_dbbg(i),dFc_dagg(i),dFc_dbgg(i),dFc_dabg(i)
!print *,dHdag,dHdaag,d,2.d0*dHdag + d*dHdaag
!print *,2.d0*dHdbg + d*dHdbbg
!print *,dHdag + dHdbg + d*dHdabg
               
           
        end do
     end if
     
end subroutine cr_ggapbe_paw_drv2_3D
! ==============================================================================

subroutine cr_ggapbe_paw_drv3(nrc,num_ir,irs,nspin,chgrhr_l,grad_trho,exc,dF_drho,dF_dgradrho &
                                ,dFc_daa,dFc_dbb,dFc_dgg,dFc_dab,dFc_dag,dFc_dbg &
                                ,dFc_daaa,dFc_dbbb,dFc_dggg &
                                ,dFc_daab,dFc_daag,dFc_dabb,dFc_dbbg,dFc_dagg,dFc_dbgg,dFc_dabg)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : PAI,DP
  implicit none

  integer,intent(in)        :: nrc,num_ir,irs(num_ir),nspin
  real(kind=DP),intent(in)  :: chgrhr_l(nrc,nspin)
  real(kind=DP),intent(in) :: grad_trho(nrc)
  real(kind=DP),intent(inout) :: exc(nrc)
  real(kind=DP),intent(out) :: dF_drho(nrc,nspin)
  real(kind=DP),intent(out) :: dF_dgradrho(nrc)
  real(kind=DP),intent(out) :: dFc_daa(nrc)
  real(kind=DP),intent(out) :: dFc_dbb(nrc)
  real(kind=DP),intent(out) :: dFc_dgg(nrc)
  real(kind=DP),intent(out) :: dFc_dab(nrc)
  real(kind=DP),intent(out) :: dFc_dag(nrc)
  real(kind=DP),intent(out) :: dFc_dbg(nrc)
  real(kind=DP),intent(out) :: dFc_daaa(nrc)
  real(kind=DP),intent(out) :: dFc_dbbb(nrc)
  real(kind=DP),intent(out) :: dFc_dggg(nrc)
  real(kind=DP),intent(out) :: dFc_daab(nrc)
  real(kind=DP),intent(out) :: dFc_daag(nrc)
  real(kind=DP),intent(out) :: dFc_dabb(nrc)
  real(kind=DP),intent(out) :: dFc_dbbg(nrc)
  real(kind=DP),intent(out) :: dFc_dagg(nrc)
  real(kind=DP),intent(out) :: dFc_dbgg(nrc)
  real(kind=DP),intent(out) :: dFc_dabg(nrc)

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
  real(kind=DP), parameter  :: thrd7= 2.333333333333d0
  real(kind=DP), parameter  :: sixth13= 2.1666666666666666d0
  real(kind=DP), parameter  :: sixth19= 3.1666666666666666d0
  real(kind=DP), parameter  :: ninth= 0.111111111111111d0
  real(kind=DP), parameter  :: twntysvnth=  3.703703703703703d-2
  real(kind=DP), parameter  :: fftyfrth= 1.851851851851852d-2
  real(kind=DP), parameter  :: eityonth= 1.234567901234568d-2

  integer       :: is,i,nr
  real(kind=DP) :: facw,bet,delt,exc1,g,g3,facpon,d,dd,rs,fk,sk,t,s,q0 &
       &         , rs12,rs32,rsp,q1,q2,eu,q3,eurs,ecd,pon &
       &         , b,b2,t2,t4,t6,q4,q5,h0,h,q8,h0t,h0b,h0rs,ht,hrs &
       &         , ec1,ec1d,ec1dd,exc0,excd,excdd &
       &         , thrd2,thrd4,fzdd,zeta,onpzeta,onmzeta,fz,zetadxd &
       &         , q0p,q1p,q2p,q3p,q0q,q1q,q2q,q3q,fzd,onzt4,gzt,bzt,h0zt,euzt &
       &         , zetadxda,zetadxdb,ecda,ecdb,ec1da,ec1db,excda,excdb
  real(kind=DP) :: rsdn,rsdnn,rsdnnn
  real(kind=DP) :: q1drr,q1drrr
  real(kind=DP) :: eudrr,eudrrr,invq12pq1,q12p1
  real(kind=DP) :: eudn,eudnn,eudnnn,ecdd,ecddd
  real(kind=DP) :: dtdn,dtdg,dtdnn,dtdgg,dtdng,dtdnnn,dtdnng,dtdngg,dtdggg
  real(kind=DP) :: dAde,dAdee,dAdeee
  real(kind=DP) :: dAdn,dAdnn,dAdnnn
  real(kind=DP) :: fE,fEdt,fEdtt,fEda,fEdaa
  real(kind=DP) :: fF,fFdt,fFdtt,fFda,fFdaa
  real(kind=DP) :: fG,fGdt,fGdtt,fGda,fGdaa
  real(kind=DP) :: at2,a2t4,a3t6,t3,t5,t7,t8
  real(kind=DP) :: dHda,dHdt,dHdaa,dHdat,dHdtt,dHdaaa,dHdaat,dHdatt,dHdttt,betivE,betivE2
  real(kind=DP) :: dHdn,dHdg,dHdnn,dHdng,dHdgg,dHdnnn,dHdnng,dHdngg,dHdggg
  real(kind=DP) :: zetada,zetadb,zetadaa,zetadab,zetadbb,zetadaaa,zetadaab,zetadabb,zetadbbb
  real(kind=DP) :: id2,id3,id4
  real(kind=DP) :: e1dr,e1drr,e1drrr
  real(kind=DP) :: e1pdr,e1pdrr,e1pdrrr
  real(kind=DP) :: e1qdr,e1qdrr,e1qdrrr
  real(kind=DP) :: q1pdrr,q1pdrrr,invq1p2pq1p,q1p2p1
  real(kind=DP) :: q1qdrr,q1qdrrr,invq1q2pq1q,q1q2p1
  real(kind=DP) :: fzdd0,fzddd,opzthrdm2,omzthrdm2
  real(kind=DP) :: zeta2,zeta3,zeta4,eudrz,eudzz,eudrrz,eudrzz,eudzzz
  real(kind=DP) :: e1,e2,e3
  real(kind=DP) :: euda,eudb,eudaa,eudbb,eudab,eudaaa,eudaab,eudabb,eudbbb
  real(kind=DP) :: rsdn2,rsdn3,zetada2,zetada3,zetadb2,zetadb3
  real(kind=DP) :: gzz,gzzz,duda,dudb,dudaa,dudab,dudbb,dudaaa,dudaab,dudabb,dudbbb
  real(kind=DP) :: dtdu,dtduu,dtduuu,dtdnu,dtdnnu,dtdnuu,dtdug,dtduug,dtdugg,dtdnug,invg
  real(kind=DP) :: dtda,dtdb,dtdaa,dtdab,dtdag,dtdbg,dtdbb
  real(kind=DP) :: dtdaaa,dtdaab,dtdabb,dtdbbb,dtdaag,dtdagg,dtdbbg,dtdbgg,dtdabg
  real(kind=DP) :: invg3,m3u
  real(kind=DP) :: dAdu,dAduu,dAduuu,dAdeu,dAdeeu,dAdeuu
  real(kind=DP) :: dAda,dAdb,dAdaa,dAdab,dAdbb,dAdaaa,dAdaab,dAdabb,dAdbbb
  real(kind=DP) :: dHdf,dHdff,dHdft,dHdfff,dHdfft,dHdftt
  real(kind=DP) :: dHdu,dHduu,dHduuu,dHdfu,dHdtu,dHdffu,dHdfuu,dHdttu,dHdtuu,dHdftu
  real(kind=DP) :: tiu,siu2
  real(kind=DP) :: dHdb,dHdbb,dHdag,dHdbg,dHdab,dHdbbb,dHdaab,dHdabb,dHdbbg,dHdbgg,dHdaag,dHdagg,dHdabg

  facw = nspin
  bet = xnu*cc0
  delt = bet / gamma
     if ( nspin == 1 ) then
        g   = 1.d0
        g3  = g**3
!!$        g4  = g3*g
        facpon = -delt/(g3*bet)

!        do i = 1,nrc,dnr ! MPI
        do nr=1,num_ir
           i=irs(nr)
           d  = facw*chgrhr_l(i, 1)
           if(d < density_minimum) cycle
! gradient of charge density 95/12/3 H.S.
           dd = facw* grad_trho(i)
           rs = (0.75/(PAI*d))**thrd
           fk = (3*PAI*PAI*d)**thrd
           sk = dsqrt(4*fk/PAI)
           if(d > density_minimum2) then
              t = dd/(d*sk*2)
              s = dd/(d*fk*2)
           else
              t = 0.d0
              s = 0.d0
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
           dF_drho(i, 1) = dF_drho(i, 1) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(dabs(grad_trho(i)) > 1.0d-9) then
!              dF_dgradrho(i) = excdd / grad_trho(i)
              dF_dgradrho(i) = excdd
           else
              dF_dgradrho(i) = 0.d0
           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d/facw
!           exc1 = exc1 + exc0*wos(i)
           exc(i)=exc(i)+exc0
        
! ***** Calc heigher derivatives. *****
           
           rsdn=-thrd*rs/d
           rsdnn=-4.d0*thrd*rsdn/d
           rsdnnn=-thrd7*rsdnn/d
!print '(i2,5e19.10)',i,d,rs,rsdn,rsdnn,rsdnnn
!q3 = a*(b1/rs12+2*b2n+3*b3*rs12+2*b4*p1*rsp)   
           q1drr =0.5d0*a*( -b1/rs32+3.d0*b3/rs12+8.d0*b4 )
           q1drrr=0.75d0*a*( b1/rs - b3 )/rs32
!eurs = -2*a*a1*q2 - q0*q3/(q1**2 + q1)
           invq12pq1=1.d0/(q1**2 + q1)
           q12p1=2.d0*q1+1.d0
           eudrr=4.d0*a*a1*q3*invq12pq1-q0*( q1drr-q12p1*q3**2*invq12pq1 )*invq12pq1
           eudrrr=6.d0*a*a1*( q1drr - q12p1*q3**2*invq12pq1 ) * invq12pq1 &
                    -q0*( q1drrr -(2.d0*q3**3+3.d0*q12p1*q3*q1drr -  &
                    2.d0*q12p1**2*q3**3* invq12pq1)*invq12pq1 ) * invq12pq1
           eudn=eurs*rsdn
           eudnn=eudrr*rsdn**2+eurs*rsdnn
           eudnnn=eudrrr*rsdn**3+3.d0*eudrr*rsdn*rsdnn+eurs*rsdnnn
!print '(i2,5e19.12)',i,d,eu,eudn,eudnn,eudnnn
           ecdd=2.d0*eudn+d*eudnn
           ecddd=3.d0*eudnn+d*eudnnn
!print '(i2,5e19.12)',i,d,eu*d,ecdd,ecddd
           
           dtdn  = -sixth7*t/d
           dtdnn = -sixth13*dtdn/d
           dtdnnn= -sixth19*dtdnn/d
!print '(i2,5e19.12)',i,d,t,dtdn,dtdnn,dtdnnn
           if(dabs(grad_trho(i)) > 1.0d-9) then
               dtdg  = t/dd
               dtdng = -sixth7*dtdg/d
               dtdnng= -sixth13*dtdng/d
           else
               dtdg  = 0.d0
               dtdng = 0.d0
               dtdnng= 0.d0
           endif
           dtdngg = 0.d0
           dtdggg = 0.d0
           
           dAde   = b * ( b + delt ) / bet
           dAdee  = 2.d0 * dAde * ( b + 0.5d0 * delt ) / bet
           dAdeee = 2.d0 * ( dAde**2 + ( b + 0.5d0 * delt )*dAdee ) /bet
           
           dAdn   = dAde  *eudn
           dAdnn  = dAdee *eudn**2  + dAde*eudnn
           dAdnnn = dAdeee*eudn**3  + 3.d0*dAdee*eudn*eudnn + dAde*eudnnn   
!print '(i2,5e19.10)' ,i,d,b,dAdnn
           
           at2  = b * t2
           a2t4 = at2**2
           a3t6 = at2*a2t4
           t3   = t2*t
           t5   = t4*t
           t7   = t6*t
           t8   = t7*t
           fE    = q5**2 + delt*t2*q4*q5
           fEdt  = 4.d0*b*t*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt *2.d0*t*(1.d0 + 4.d0*at2 + 6.d0*a2t4 + 4.d0*a3t6 )
           fEdtt = 4.d0*b*( 1.d0 + 9.d0*at2 + 15.d0*a2t4 + 14.d0*a3t6 ) + &
                    delt*( 2.d0 + 24.d0*at2 + 60.d0*a2t4 + 56.d0*a3t6 )
           fEda  = 2.d0*t2*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt*t4*( 2.d0 + 4.d0*at2 + 3.d0*a2t4 )
           fEdaa = 6.d0*t4*( 2.d0*q5 - 1.d0) + &
                    delt*2.d0*t6*( 2.d0 + 3.d0*at2 )
           fF    = t*( 1.d0 + 2.d0*at2 )         
           fFdt  = 1.d0 + 6.d0*at2
           fFdtt = 12.d0*b*t
           fFda  = 2.d0*t3
           fFdaa = 0.d0
           fG    = 2.d0*b*t6 + b2*t8
           fGdt  = 12.d0*b*t5 + 8.d0*b2*t7
           fGdtt = 60.d0*b*t4 + 56.d0*b2*t6
           fGda  = 2.d0*t6 + 2.d0*b*t8
           fGdaa = 2.d0*t8
           
           betivE=bet/fE
           betivE2=betivE/fE
           dHda   =     -fG*betivE
           dHdt   = 2.d0*fF*betivE
           dHdaa  =     -( fGda*fE-fG*fEda )*betivE2
           dHdat  =     -( fGdt*fE-fG*fEdt )*betivE2
           dHdtt  = 2.d0*( fFdt*fE-fF*fEdt )*betivE2
           dHdaaa =     -( fGdaa*fE-fG*fEdaa-2.d0*( fGda*fE-fG*fEda )*fEda/fE )*betivE2
           dHdaat = 2.d0*( fFdaa*fE-fF*fEdaa-2.d0*( fFda*fE-fF*fEda )*fEda/fE )*betivE2
           dHdatt =     -( fGdtt*fE-fG*fEdtt-2.d0*( fGdt*fE-fG*fEdt )*fEdt/fE )*betivE2
           dHdttt = 2.d0*( fFdtt*fE-fF*fEdtt-2.d0*( fFdt*fE-fF*fEdt )*fEdt/fE )*betivE2
           
           dHdn     = dHda  *dAdn       + dHdt*dtdn
           dHdg     = dHdt  *dtdg
           dHdnn    = dHdaa *dAdn**2    + 2.d0*dHdat*dAdn*dtdn  + dHdtt*dtdn**2   + dHda*dAdnn  + dHdt*dtdnn
           dHdng    = dHdat *dAdn*dtdg  + dHdtt*dtdn*dtdg       + dHdt*dtdng
           dHdgg    = dHdtt*dtdg**2
           dHdnnn   = dHdaaa*dAdn**3    + dHdttt*dtdn**3 &
                        + 3.d0*dHdaat*dAdn**2*dtdn  + 3.d0*dHdatt*dAdn*dtdn**2 &
                        + 3.d0*dHdaa*dAdn*dAdnn + 3.d0*dHdtt*dtdn*dtdnn + 3.d0*dHdat*(dAdnn*dtdn + dAdn*dtdnn) &
                        + dHda*dAdnnn   + dHdt*dtdnnn
           dHdnng   = dHdaat*dAdn**2*dtdg   + 2.d0*dHdatt*dAdn*dtdn*dtdg    + dHdttt*dtdn**2*dtdg &
                        + dHdat*(dAdnn*dtdg + 2.d0*dAdn*dtdng)  + dHdtt*(dtdnn*dtdg + 2.d0*dtdn*dtdng) &
                        + dHdt*dtdnng
           dHdngg   = dHdatt*dAdn*dtdg**2   + dHdttt*dtdn*dtdg**2   + 2.d0*dHdtt*dtdg*dtdng
           dHdggg   = dHdttt*dtdg**3
!dF_drho(i, 1)=h+d*dHdn+ecd
           dFc_daa(i)   = 2.d0*dHdn + d*dHdnn       + ecdd
           dFc_dag(i)   = dHdg  + d*dHdng
           dFc_dgg(i)   = d*dHdgg
#ifdef USE_3RD_ORDER_DERIVATIVE
           dFc_daaa(i)  = 3.d0*dHdnn    + d*dHdnnn  + ecddd
           dFc_daag(i)  = 2.d0*dHdng    + d*dHdnng
           dFc_dagg(i)  = dHdgg         + d*dHdngg
           dFc_dggg(i)  = d*dHdggg  
#endif
        end do
     else if ( nspin ==  2 ) then
        thrd2 = thrd * 2
        thrd4 = thrd * 4
        fzdd = 8/(9*(2**thrd4 - 2 ))

!        do i = 1,nrc,dnr
        do nr=1,num_ir
           i=irs(nr)
           d = chgrhr_l(i, 1) + chgrhr_l(i, nspin)
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
           ! ****** aas mkatsu *****
!              t = dd/(d*sk*2)
              t = dd/(d*sk*2)/g         !    mkatsu aas
              s = dd/(d*fk*2)
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

!           zetadxd = - 2 * chgrhr_l(i, 3-is) * (-1.d0)**is / d
           zetadxda = - 2 * chgrhr_l(i, 2) * (-1.d0) / d
           zetadxdb = - 2 * chgrhr_l(i, 1) / d
!           ecd = eu-thrd*rs*eurs + euzt * zetadxd
           ecda = eu-thrd*rs*eurs + euzt * zetadxda 
           ecdb = eu-thrd*rs*eurs + euzt * zetadxdb
           if(onpzeta < zeta_minimum2 .or. onmzeta < zeta_minimum2) then
              ec1  = 0.d0;    ec1d = 0.d0;     ec1dd= 0.d0
              ec1da = 0.d0;     ec1db = 0.d0
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
!           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxd
           ec1da = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxda
           ec1db = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxdb
           ec1dd = 0.5*ht/sk / g

10         continue
           exc0 = eu*d * 0.5 + ec1
!           excd = ecd + ec1d
           excda= ecda + ec1da
           excdb= ecdb + ec1db
!           dF_drho(i, is) = dF_drho(i, is) + excd
           dF_drho(i, 1) = dF_drho(i, 1) + excda
           dF_drho(i, 2) = dF_drho(i, 2) + excdb
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
!           if(is ==  2) then
              if(dabs(grad_trho(i)) > 1.0d-9) then
!                 dF_dgradrho(i) = excdd / grad_trho(i)
                 dF_dgradrho(i) = excdd
              else
                 dF_dgradrho(i) = 0.d0
              endif
!           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d * 0.5
!           exc1 = exc1 + exc0*wos(i)
           exc(i)=exc(i)+exc0*2.d0
           
! ***** Calc heigher order derivatives. *****

           rsdn=-thrd*rs/d
           rsdnn=-4.d0*thrd*rsdn/d
           rsdnnn=-thrd7*rsdnn/d
           
           id2=1.d0/d/d
           id3=id2/d
           id4=id3/d
           zetada   =     zetadxda/d
           zetadb   =     zetadxdb/d
           zetadaa  =    -4.d0*chgrhr_l(i, 2)*id3
           zetadab  =     2.d0*zeta*id2
           zetadbb  =     4.d0*chgrhr_l(i, 1)*id3
           zetadaaa =    12.d0*chgrhr_l(i, 2)*id4
           zetadaab =   ( -4.d0*chgrhr_l(i, 1) + 8.d0*chgrhr_l(i, 2))*id4
           zetadabb =   ( -8.d0*chgrhr_l(i, 1) + 4.d0*chgrhr_l(i, 2))*id4
           zetadbbb =   -12.d0*chgrhr_l(i, 1)*id4
           
!print '(i2,7d19.10)', i,chgrhr_l(i,1:2),zetadab,zetadaaa,zetadaab,zetadabb,zetadbbb
             
           q1drr =0.5d0*a*( -b1/rs32+3.d0*b3/rs12+8.d0*b4 )
           q1drrr=0.75d0*a*( b1/rs - b3 )/rs32
           invq12pq1=1.d0/(q1**2 + q1)
           q12p1=2.d0*q1+1.d0
           e1dr=-2.d0*a*a1*q2-q0*q3*invq12pq1
           e1drr=4.d0*a*a1*q3*invq12pq1-q0*( q1drr-q12p1*q3**2*invq12pq1 )*invq12pq1
           e1drrr=6.d0*a*a1*( q1drr - q12p1*q3**2*invq12pq1 ) * invq12pq1 &
                    -q0*( q1drrr -(2.d0*q3**3+3.d0*q12p1*q3*q1drr -  &
                    2.d0*q12p1**2*q3**3* invq12pq1)*invq12pq1 ) * invq12pq1
                    
           q1pdrr =0.5d0*ap*( -b1p/rs32+3.d0*b3p/rs12+8.d0*b4p )
           q1pdrrr=0.75d0*ap*( b1p/rs - b3p )/rs32
           invq1p2pq1p=1.d0/(q1p**2 + q1p)
           q1p2p1=2.d0*q1p+1.d0
           e1pdr=-2.d0*ap*a1p*q2p-q0p*q3p*invq1p2pq1p
           e1pdrr=4.d0*ap*a1p*q3p*invq1p2pq1p-q0p*( q1pdrr-q1p2p1*q3p**2*invq1p2pq1p )*invq1p2pq1p
           e1pdrrr=6.d0*ap*a1p*( q1pdrr - q1p2p1*q3p**2*invq1p2pq1p ) * invq1p2pq1p &
                    -q0p*( q1pdrrr -(2.d0*q3p**3+3.d0*q1p2p1*q3p*q1pdrr -  &
                    2.d0*q1p2p1**2*q3p**3* invq1p2pq1p)*invq1p2pq1p ) * invq1p2pq1p
                    
           q1qdrr =0.5d0*aq*( -b1q/rs32+3.d0*b3q/rs12+8.d0*b4q )
           q1qdrrr=0.75d0*aq*( b1q/rs - b3q )/rs32
           invq1q2pq1q=1.d0/(q1q**2 + q1q)
           q1q2p1=2.d0*q1q+1.d0
           e1qdr=-2.d0*aq*a1q*q2q-q0q*q3q*invq1q2pq1q
           e1qdrr=4.d0*aq*a1q*q3q*invq1q2pq1q-q0q*( q1qdrr-q1q2p1*q3q**2*invq1q2pq1q )*invq1q2pq1q
           e1qdrrr=6.d0*aq*a1q*( q1qdrr - q1q2p1*q3q**2*invq1q2pq1q ) * invq1q2pq1q &
                    -q0q*( q1qdrrr -(2.d0*q3q**3+3.d0*q1q2p1*q3q*q1qdrr -  &
                    2.d0*q1q2p1**2*q3q**3* invq1q2pq1q)*invq1q2pq1q ) * invq1q2pq1q
           
           zeta2=zeta**2
           zeta3=zeta2*zeta
           zeta4=zeta2*zeta2
           
           if(dabs(zeta) < zeta_minimum) then
              fzdd0 = fzdd &
                   &            *( 1 + 5.d0*ninth *zeta2 &
                   &            *( 1 +22.d0*twntysvnth *zeta2 ))
              fzddd = fzdd*10.d0*ninth*zeta*( 1.d0 + 44.d0*twntysvnth*zeta2 )
              gzz   = -2.d0*ninth*( 1.d0 + 14.d0*ninth*zeta2*( 1.d0 + 65.d0*fftyfrth*zeta2 ) )
              gzzz  = -56.d0*eityonth*zeta*( 1.d0 + 65.d0*twntysvnth*zeta2 )
           else
              opzthrdm2=1.d0/onpzeta**(thrd2)
              omzthrdm2=1.d0/onmzeta**(thrd2)
              fzdd0 =  fzdd*0.5d0*( opzthrdm2 + omzthrdm2 )
              fzddd = -fzdd*thrd*( opzthrdm2/onpzeta - omzthrdm2/onmzeta )
              gzz   = -ninth*( opzthrdm2**2 + omzthrdm2**2 )
!              gzz   = -ninth*( 1.d0/onpzeta**(thrd4) + 1.d0/onmzeta**(thrd4) )
              gzzz  = 4.d0*twntysvnth*( opzthrdm2**2/onpzeta - omzthrdm2**2/onmzeta )
           end if
           
           e1=q0*q2
           e2=q0p*q2p
           e3=q0q*q2q
           
           eudrr    = e1drr  + ( e1pdrr  -  e1drr  ) * fz * zeta4 - e1qdrr  * fz * onzt4 / fzdd
           eudrrr   = e1drrr + ( e1pdrrr -  e1drrr ) * fz * zeta4 - e1qdrrr * fz * onzt4 / fzdd
           eudzz    = ( e2 - e1 )*( fzdd0*zeta4 + 8.d0*fzd*zeta3 + 12.d0*fz*zeta2 ) &
                        - e3*( fzdd0*onzt4 - 8.d0*fzd*zeta3 - 12.d0*fz*zeta2 )/fzdd
           eudzzz   = ( e2 - e1 )*( fzddd*zeta4 + 12.d0*fzdd0*zeta3 + 36.d0*fzd*zeta2 + 24.d0*fz*zeta ) &
                        - e3*( fzddd*onzt4 - 12.d0*fzdd0*zeta3 - 36.d0*fzd*zeta2 - 24.d0*fz*zeta )/fzdd
           eudrz    = ( e1pdr - e1dr )*( fzd*zeta4 + 4.d0*fz*zeta3 ) &
                        - e1qdr *( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd 
           eudrrz   = ( e1pdrr  - e1drr  )*( fzd*zeta4 + 4.d0*fz*zeta3 ) &
                        - e1qdrr *( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd
           eudrzz   = ( e1pdr - e1dr )*( fzdd0*zeta4 + 8.d0*fzd*zeta3 + 12.d0*fz*zeta2 ) &
                        - e1qdr *( fzdd0*onzt4 - 8.d0*fzd*zeta3 - 12.d0*fz*zeta2 )/fzdd
                        
           rsdn2    = rsdn*rsdn
           rsdn3    = rsdn2*rsdn
           zetada2  = zetada**2
           zetada3  = zetada2*zetada
           zetadb2  = zetadb**2
           zetadb3  = zetadb2*zetadb
                        
           euda     = eurs*rsdn + euzt*zetada
           eudb     = eurs*rsdn + euzt*zetadb
           eudaa    = eudrr*rsdn2 + 2.d0*eudrz*rsdn*zetada + eudzz*zetada2 + eurs*rsdnn + euzt*zetadaa
           eudab    = eudrr*rsdn2 + eudrz*rsdn*( zetada + zetadb ) + eudzz*zetada*zetadb &
                        + eurs*rsdnn + euzt*zetadab
           eudbb    = eudrr*rsdn2 + 2.d0*eudrz*rsdn*zetadb + eudzz*zetadb2 + eurs*rsdnn + euzt*zetadbb
           eudaaa   = eudrrr*rsdn3 &
                        + 3.d0*eudrrz*rsdn2*zetada &
                        + 3.d0*eudrzz*rsdn*zetada2 + eudzzz*zetada3 &
                        + 3.d0*eudrr*rsdn*rsdnn + 3.d0*eudzz*zetada*zetadaa &
                        + 3.d0*eudrz*( rsdnn*zetada + rsdn*zetadaa ) &
                        + eurs*rsdnnn + euzt*zetadaaa
           eudaab   = eudrrr*rsdn3 &
                        + eudrrz*rsdn2*( zetadb + 2.d0*zetada ) &
                        + eudrzz*rsdn*zetada*( zetada + 2.d0*zetadb ) &
                        + eudzzz*zetada2*zetadb &
                        + 3.d0*eudrr*rsdn*rsdnn &
                        + eudrz*( rsdnn*( 2.d0*zetada + zetadb ) + rsdn*( 2.d0*zetadab + zetadaa ) ) &
                        + eudzz*( 2.d0*zetada*zetadab + zetadaa*zetadb ) &
                        + eurs*rsdnnn + euzt*zetadaab
           eudabb   = eudrrr*rsdn3 &
                        + eudrrz*rsdn2*( zetada + 2.d0*zetadb ) &
                        + eudrzz*rsdn*zetadb*( zetadb + 2.d0*zetada ) &
                        + eudzzz*zetada*zetadb2 &
                        + 3.d0*eudrr*rsdn*rsdnn &
                        + eudrz*( rsdnn*( zetada + 2.d0*zetadb ) + rsdn*( 2.d0*zetadab + zetadbb ) ) &
                        + eudzz*( zetada*zetadbb + 2.d0*zetadb*zetadab ) &
                        + eurs*rsdnnn + euzt*zetadabb
           eudbbb   = eudrrr*rsdn3 &
                        + 3.d0*eudrrz*rsdn2*zetadb &
                        + 3.d0*eudrzz*rsdn*zetadb2 + eudzzz*zetadb3 &
                        + 3.d0*eudrr*rsdn*rsdnn + 3.d0*eudzz*zetadb*zetadbb &
                        + 3.d0*eudrz*( rsdnn*zetadb + rsdn*zetadbb ) &
                        + eurs*rsdnnn + euzt*zetadbbb
                        
!print '(i2,12d19.10)', i,chgrhr_l(i,1:2),eu,euda,eudb,eudaa,eudbb,eudab,eudaaa,eudaab,eudabb,eudbbb
!print *,eu,     e1+(e2-e1)* fz * zeta4 - e3  * fz * onzt4 / fzdd
!print *,eurs,   e1dr  + ( e1pdr  -  e1dr  ) * fz * zeta4 - e1qdr  * fz * onzt4 / fzdd
!print *,euzt,   (e2-e1)*( fzd*zeta4 + 4.d0*fz*zeta3 )-e3*( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd 

           duda     = gzt*zetada
           dudb     = gzt*zetadb 
           dudaa    = gzz*zetada2     + gzt*zetadaa
           dudab    = gzz*zetada*zetadb + gzt*zetadab
           dudbb    = gzz*zetadb2     + gzt*zetadbb
           dudaaa   = gzzz*zetada3    + 3.d0*gzz*zetada*zetadaa   + gzt*zetadaaa
           dudaab   = gzzz*zetada2*zetadb + gzz*( zetadaa*zetadb + 2.d0*zetada*zetadab ) + gzt*zetadaab
           dudabb   = gzzz*zetada*zetadb2 + gzz*( zetada*zetadbb + 2.d0*zetadab*zetadb ) + gzt*zetadabb
           dudbbb   = gzzz*zetadb3    + 3.d0*gzz*zetadb*zetadbb   + gzt*zetadbbb 
           
!print '(i2,12d19.10)', i,chgrhr_l(i,1:2),g,duda,dudb,dudaa,dudbb,dudab,dudaaa,dudaab,dudabb,dudbbb

           
           dtdn  = -sixth7*t/d
           dtdnn = -sixth13*dtdn/d
           dtdnnn= -sixth19*dtdnn/d
!print '(i2,5e19.12)',i,d,t,dtdn,dtdnn,dtdnnn
           if(dabs(grad_trho(i)) > 1.0d-9) then
               dtdg  = t/dd
               dtdng = -sixth7*dtdg/d
               dtdnng= -sixth13*dtdng/d
!               dtdnug = -dtdn*invg/dd
           else
               dtdg  = 0.d0
               dtdng = 0.d0
               dtdnng= 0.d0
               dtdnug= 0.d0 
           endif
           dtdgg  = 0.d0
           dtdngg = 0.d0
           dtdggg = 0.d0
           invg=1.d0/g
           dtdu   = -t*invg
           dtduu  = -2.d0*dtdu*invg 
           dtduuu = -3.d0*dtduu*invg
           dtdnu  = -dtdn*invg
           dtdnnu = -dtdnn*invg
           dtdnuu = -sixth7*dtduu/d
           dtdug  = -dtdg*invg
           dtduug = 2.d0*dtdg*invg**2
           dtdugg = 0.d0
           if(dabs(grad_trho(i)) > 1.0d-9) then
                dtdnug = dtdnu/dd
           else
                dtdnug = 0.d0
           end if
           
           dtda   = dtdn + dtdu*duda
           dtdb   = dtdn + dtdu*dudb
           dtdaa  = dtdnn + 2.d0*dtdnu*duda + dtduu*duda**2 + dtdu*dudaa
           dtdbb  = dtdnn + 2.d0*dtdnu*dudb + dtduu*dudb**2 + dtdu*dudbb
           dtdab  = dtdnn + dtdnu*( duda + dudb ) + dtduu*duda*dudb + dtdu*dudab
           dtdag  = dtdng + dtdug*duda
           dtdbg  = dtdng + dtdug*dudb
           dtdaaa = dtdnnn+ 3.d0*dtdnnu*duda + 3.d0*dtdnuu*duda**2 + dtduuu*duda**3 &
                    + 3.d0*dtduu*duda*dudaa + 3.d0*dtdnu*dudaa + dtdu*dudaaa
           dtdbbb = dtdnnn+ 3.d0*dtdnnu*dudb + 3.d0*dtdnuu*dudb**2 + dtduuu*dudb**3 &
                    + 3.d0*dtduu*dudb*dudbb + 3.d0*dtdnu*dudbb + dtdu*dudbbb
           dtdaab = dtdnnn+ dtdnnu*( 2.d0*duda + dudb ) + dtdnuu*duda*( duda + 2.d0*dudb ) + dtduuu*duda**2*dudb &
                    + dtduu*( dudaa*dudb + 2.d0*duda*dudab ) + dtdnu*( dudaa + 2.d0*dudab ) + dtdu*dudaab
           dtdabb = dtdnnn+ dtdnnu*( duda + 2.d0*dudb ) + dtdnuu*dudb*( dudb + 2.d0*duda ) + dtduuu*dudb**2*duda &
                    + dtduu*( dudbb*duda + 2.d0*dudb*dudab ) + dtdnu*( dudbb + 2.d0*dudab ) + dtdu*dudabb 
           dtdaag = dtdnng+ 2.d0*dtdnug*duda + dtduug*duda**2 + dtdug*dudaa
           dtdbbg = dtdnng+ 2.d0*dtdnug*dudb + dtduug*dudb**2 + dtdug*dudbb
           dtdabg = dtdnng+ dtdnug*( duda + dudb ) + dtduug*duda*dudb + dtdug*dudab 
           
!print '(i2,20d19.10)', i,t,dtda,dtdb,dtdg,dtdaa,dtdbb,dtdgg,dtdab,dtdag,dtdbg,dtdaaa,dtdbbb,dtdggg,dtdaab,dtdabb,dtdaag,dtdagg,dtdbbg,dtdbgg,dtdabg     
!print '(i2,12d19.10)', i,t,dtdag,dtdbg,dtdaag,dtdbbg,dtdabg     

           invg3  = 1.d0  /g3
           dAde   = b * ( b + delt ) / bet * invg3
           dAdee  = 2.d0 * dAde * ( b + 0.5d0 * delt ) / bet *invg3
           dAdeee = 2.d0 * ( dAde**2 + ( b + 0.5d0 * delt )*dAdee ) /bet * invg3
           m3u    = -3.d0*invg
           dAdu     = m3u*eu*dAde
           dAdeu    = m3u*( dAde + eu*dAdee )
           dAduu    = m3u*eu*invg*( -dAde + g*dAdeu )
           dAdeeu   = m3u*( 2.d0*dAdee + eu*dAdeee )
           dAdeuu   = m3u*invg*( -dAde + g*dAdeu - eu*dAdee + g*eu*dAdeeu )
           dAduuu   = m3u*eu*invg**2*( 2.d0*dAde -2.d0*g*dAdeu + g**2*dAdeuu )
           
           dAda     = dAde*euda + dAdu*duda
           dAdb     = dAde*eudb + dAdu*dudb
           dAdaa    = dAdee*euda**2 + 2.d0*dAdeu*euda*duda + dAduu*duda**2 + dAde*eudaa + dAdu*dudaa 
           dAdab    = dAdee*euda*eudb + dAdeu*( euda*dudb + eudb*duda ) + dAduu*duda*dudb + dAde*eudab + dAdu*dudab 
           dAdbb    = dAdee*eudb**2 + 2.d0*dAdeu*eudb*dudb + dAduu*dudb**2 + dAde*eudbb + dAdu*dudbb
           dAdaaa   = dAdeee*euda**3 + 3.d0*dAdeeu*euda**2*duda + 3.d0*dAdeuu*euda*duda**2 + dAduuu*duda**3 &
                        + 3.d0*dAdee*euda*eudaa + 3.d0*dAduu*duda*dudaa &
                        + 3.d0*dAdeu*( euda*dudaa + eudaa*duda ) + dAde*eudaaa + dAdu*dudaaa
           dAdaab   = dAdeee*euda**2*eudb + dAduuu*duda**2*dudb &
                        + dAdeeu*euda*( euda*dudb + 2.d0*eudb*duda ) &
                        + dAdeuu*( 2.d0*euda*dudb + eudb*duda )*duda &
                        + dAdee*( 2.d0*euda*eudab + eudb*eudaa ) + dAduu*( 2.d0*duda*dudab + dudb*dudaa ) &
                        + dAdeu*( 2.d0*eudab*duda + 2.d0*euda*dudab + dudb*eudaa + eudb*dudaa ) &
                        + dAde*eudaab + dAdu*dudaab
           dAdabb   = dAdeee*euda*eudb**2 + dAduuu*dudb**2*duda &
                        + dAdeeu*eudb*( eudb*duda + 2.d0*euda*dudb ) &
                        + dAdeuu*( 2.d0*eudb*duda + euda*dudb )*dudb &
                        + dAdee*( eudbb*euda + 2.d0*eudb*eudab ) + dAduu*( dudbb*duda + 2.d0*dudb*dudab ) &
                        + dAdeu*( eudbb*duda + 2.d0*eudb*dudab + 2.d0*dudb*eudab + dudbb*euda ) &
                        + dAde*eudabb + dAdu*dudabb
           dAdbbb   = dAdeee*eudb**3 + 3.d0*dAdeeu*eudb**2*dudb + 3.d0*dAdeuu*eudb*dudb**2 + dAduuu*dudb**3 &
                        + 3.d0*dAdee*eudb*eudbb + 3.d0*dAduu*dudb*dudbb &
                        + 3.d0*dAdeu*( eudb*dudbb + eudbb*dudb ) + dAde*eudbbb + dAdu*dudbbb 
                        
!print '(i2,12d19.10)', i,b,dAda,dAdb,dAdaa,dAdab,dAdbb,dAdaaa,dAdaab,dAdabb,dAdbbb

           at2  = b * t2
           a2t4 = at2**2
           a3t6 = at2*a2t4
           t3   = t2*t
           t5   = t4*t
           t7   = t6*t
           t8   = t7*t
           fE    = q5**2 + delt*t2*q4*q5
           fEdt  = 4.d0*b*t*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt *2.d0*t*(1.d0 + 4.d0*at2 + 6.d0*a2t4 + 4.d0*a3t6 )
           fEdtt = 4.d0*b*( 1.d0 + 9.d0*at2 + 15.d0*a2t4 + 14.d0*a3t6 ) + &
                    delt*( 2.d0 + 24.d0*at2 + 60.d0*a2t4 + 56.d0*a3t6 )
           fEda  = 2.d0*t2*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt*t4*( 2.d0 + 4.d0*at2 + 3.d0*a2t4 )
           fEdaa = 6.d0*t4*( 2.d0*q5 - 1.d0) + &
                    delt*2.d0*t6*( 2.d0 + 3.d0*at2 )
           fF    = t*( 1.d0 + 2.d0*at2 )         
           fFdt  = 1.d0 + 6.d0*at2
           fFdtt = 12.d0*b*t
           fFda  = 2.d0*t3
           fFdaa = 0.d0
           fG    = 2.d0*b*t6 + b2*t8
           fGdt  = 12.d0*b*t5 + 8.d0*b2*t7
           fGdtt = 60.d0*b*t4 + 56.d0*b2*t6
           fGda  = 2.d0*t6 + 2.d0*b*t8
           fGdaa = 2.d0*t8
           
           betivE=bet/fE*g3
           betivE2=betivE/fE
           dHdf   =     -fG*betivE
           dHdt   = 2.d0*fF*betivE
           dHdff  =     -( fGda*fE-fG*fEda )*betivE2
           dHdft  =     -( fGdt*fE-fG*fEdt )*betivE2
           dHdtt  = 2.d0*( fFdt*fE-fF*fEdt )*betivE2
           dHdfff =     -( fGdaa*fE-fG*fEdaa-2.d0*( fGda*fE-fG*fEda )*fEda/fE )*betivE2
           dHdfft = 2.d0*( fFdaa*fE-fF*fEdaa-2.d0*( fFda*fE-fF*fEda )*fEda/fE )*betivE2
           dHdftt =     -( fGdtt*fE-fG*fEdtt-2.d0*( fGdt*fE-fG*fEdt )*fEdt/fE )*betivE2
           dHdttt = 2.d0*( fFdtt*fE-fF*fEdtt-2.d0*( fFdt*fE-fF*fEdt )*fEdt/fE )*betivE2
           
           tiu    = 3.d0*invg
           siu2   = tiu*2.d0*invg
           dHdu   = tiu*h
           dHduu  = siu2*h
           dHduuu = siu2*invg*h
           dHdfu  = tiu*dHdf
           dHdtu  = tiu*dHdt
           dHdffu = tiu*dHdff
           dHdfuu = siu2*dHdf
           dHdttu = tiu*dHdtt
           dHdtuu = siu2*dHdt
           dHdftu = tiu*dHdft
           
           dHda   = dHdf*dAda + dHdt*dtda + dHdu*duda
           dHdb   = dHdf*dAdb + dHdt*dtdb + dHdu*dudb
           dHdg   =             dHdt*dtdg
           dHdaa  = dHdff*dAda**2 + dHdtt*dtda**2 + dHduu*duda**2 &
                    + 2.d0*dHdft*dAda*dtda + 2.d0*dHdfu*dAda*duda + 2.d0*dHdtu*dtda*duda &
                    + dHdf*dAdaa + dHdt*dtdaa + dHdu*dudaa
           dHdab  = dHdff*dAda*dAdb + dHdtt*dtda*dtdb + dHduu*duda*dudb &
                    + dHdft*( dAda*dtdb + dAdb*dtda ) + dHdfu*( dAda*dudb + dAdb*duda ) + dHdtu*( dtda*dudb + dtdb*duda ) &
                    + dHdf*dAdab + dHdt*dtdab + dHdu*dudab
           dHdbb  = dHdff*dAdb**2 + dHdtt*dtdb**2 + dHduu*dudb**2 &
                    + 2.d0*dHdft*dAdb*dtdb + 2.d0*dHdfu*dAdb*dudb + 2.d0*dHdtu*dtdb*dudb &
                    + dHdf*dAdbb + dHdt*dtdbb + dHdu*dudbb
           dHdgg  = dHdtt*dtdg**2
           dHdag  = dHdft*dAda*dtdg + dHdtt*dtda*dtdg + dHdtu*duda*dtdg + dHdt*dtdag
           dHdbg  = dHdft*dAdb*dtdg + dHdtt*dtdb*dtdg + dHdtu*dudb*dtdg + dHdt*dtdbg
           dHdaaa = dHdfff*dAda**3 + 3.d0*dHdfft*dAda**2*dtda + 3.d0*dHdffu*dAda**2*duda &
                    + 3.d0*dHdftt*dAda*dtda**2 + dHdttt*dtda**3 + 3.d0*dHdttu*dtda**2*duda &
                    + 3.d0*dHdfuu*dAda*duda**2 + 3.d0*dHdtuu*dtda*duda**2 + dHduuu*duda**3 &
                    + 6.d0*dHdftu*dAda*dtda*duda &
                    + 3.d0*dHdff*dAda*dAdaa + 3.d0*dHdtt*dtda*dtdaa + 3.d0*dHduu*duda*dudaa &
                    + 3.d0*dHdft*( dAdaa*dtda + dAda*dtdaa ) &
                    + 3.d0*dHdtu*( dtdaa*duda + dtda*dudaa ) &
                    + 3.d0*dHdfu*( dAdaa*duda + dAda*dudaa ) &
                    + dHdf*dAdaaa + dHdt*dtdaaa + dHdu*dudaaa
           dHdaab = dHdfff*dAda**2*dAdb + dHdfft*dAda*( dAda*dtdb + 2.d0*dAdb*dtda ) &
                    + dHdffu*dAda*( dAda*dudb + 2.d0*dAdb*duda ) &
                    + dHdftt*dtda*( dAdb*dtda + 2.d0*dAda*dtdb ) &
                    + dHdttt*dtdb*dtda**2 + dHdttu*dtda*( dtda*dudb + 2.d0*dtdb*duda ) &
                    + dHdfuu*duda*( dAdb*duda + 2.d0*dAda*dudb ) &
                    + dHdtuu*duda*( dtdb*duda + 2.d0*dtda*dudb ) + dHduuu*duda**2*dudb &
                    + 2.d0*dHdftu*( dAdb*duda*dtda + dAda*dudb*dtda + dAda*duda*dtdb ) &
                    + dHdff*( 2.d0*dAda*dAdab + dAdb*dAdaa ) &
                    + dHdtt*( 2.d0*dtda*dtdab + dtdb*dtdaa ) &
                    + dHduu*( 2.d0*duda*dudab + dudb*dudaa ) &
                    + dHdft*( 2.d0*dAdab*dtda + 2.d0*dAda*dtdab + dAdaa*dtdb + dAdb*dtdaa ) &
                    + dHdtu*( 2.d0*dtdab*duda + 2.d0*dtda*dudab + dtdaa*dudb + dtdb*dudaa ) &
                    + dHdfu*( 2.d0*dAdab*duda + 2.d0*dAda*dudab + dAdaa*dudb + dAdb*dudaa ) &
                    + dHdf*dAdaab + dHdt*dtdaab + dHdu*dudaab
           dHdabb = dHdfff*dAda*dAdb**2 + dHdfft*dAdb*( dAdb*dtda + 2.d0*dAda*dtdb ) &
                    + dHdffu*dAdb*( dAdb*duda + 2.d0*dAda*dudb ) &
                    + dHdftt*dtdb*( dAda*dtdb + 2.d0*dAdb*dtda ) &
                    + dHdttt*dtda*dtdb**2 + dHdttu*dtdb*( dtdb*duda + 2.d0*dtda*dudb ) &
                    + dHdfuu*dudb*( dAda*dudb + 2.d0*dAdb*duda ) &
                    + dHdtuu*dudb*( dtda*dudb + 2.d0*dtdb*duda ) + dHduuu*duda*dudb**2 &
                    + 2.d0*dHdftu*( dAda*dtdb*dudb + dAdb*dtda*dudb + dAdb*dtdb*duda ) &
                    + dHdff*( 2.d0*dAdb*dAdab + dAda*dAdbb ) &
                    + dHdtt*( 2.d0*dtdb*dtdab + dtda*dtdbb ) &
                    + dHduu*( 2.d0*dudb*dudab + duda*dudbb ) &
                    + dHdft*( 2.d0*dAdab*dtdb + 2.d0*dAdb*dtdab + dAdbb*dtda + dAda*dtdbb ) &
                    + dHdtu*( 2.d0*dtdab*dudb + 2.d0*dtdb*dudab + dtdbb*duda + dtda*dudbb ) &
                    + dHdfu*( 2.d0*dAdab*dudb + 2.d0*dAdb*dudab + dAdbb*duda + dAda*dudbb ) &
                    + dHdf*dAdabb + dHdt*dtdabb + dHdu*dudabb
           dHdbbb = dHdfff*dAdb**3 + 3.d0*dHdfft*dAdb**2*dtdb + 3.d0*dHdffu*dAdb**2*dudb &
                    + 3.d0*dHdftt*dAdb*dtdb**2 + dHdttt*dtdb**3 + 3.d0*dHdttu*dtdb**2*dudb &
                    + 3.d0*dHdfuu*dAdb*dudb**2 + 3.d0*dHdtuu*dtdb*dudb**2 + dHduuu*dudb**3 &
                    + 6.d0*dHdftu*dAdb*dtdb*dudb &
                    + 3.d0*dHdff*dAdb*dAdbb + 3.d0*dHdtt*dtdb*dtdbb + 3.d0*dHduu*dudb*dudbb &
                    + 3.d0*dHdft*( dAdbb*dtdb + dAdb*dtdbb ) &
                    + 3.d0*dHdtu*( dtdbb*dudb + dtdb*dudbb ) &
                    + 3.d0*dHdfu*( dAdbb*dudb + dAdb*dudbb ) &
                    + dHdf*dAdbbb + dHdt*dtdbbb + dHdu*dudbbb
           dHdabg = dHdfft*dAda*dAdb*dtdg   + dHdftt*( dAdb*dtda + dAda*dtdb )*dtdg &
                    + dHdttt*dtda*dtdb*dtdg + dHdttu*( dtdb*duda + dtda*dudb )*dtdg &
                    + dHdtuu*duda*dudb*dtdg + dHdftu*( dAdb*duda + dAda*dudb )*dtdg &
                    + dHdft*( dAdab*dtdg + dAdb*dtdag + dAda*dtdbg ) &
                    + dHdtt*( dtdab*dtdg + dtdb*dtdag + dtda*dtdbg ) &
                    + dHdtu*( dudab*dtdg + dudb*dtdag + duda*dtdbg ) &
                    + dHdt*dtdabg
           dHdaag = dHdfft*dAda**2*dtdg + dHdttt*dtda**2*dtdg + dHdtuu*duda**2*dtdg &
                    + 2.d0*dHdttu*dtda*duda*dtdg + 2.d0*dHdftt*dAda*dtda*dtdg + 2.d0*dHdftu*dAda*duda*dtdg &
                    + dHdtt*( dtdaa*dtdg + 2.d0*dtda*dtdag ) + dHdft*( dAdaa*dtdg + 2.d0*dAda*dtdag ) &
                    + dHdtu*( dudaa*dtdg + 2.d0*duda*dtdag ) + dHdt*dtdaag         
           dHdbbg = dHdfft*dAdb**2*dtdg + dHdttt*dtdb**2*dtdg + dHdtuu*dudb**2*dtdg &
                    + 2.d0*dHdttu*dtdb*dudb*dtdg + 2.d0*dHdftt*dAdb*dtdb*dtdg + 2.d0*dHdftu*dAdb*dudb*dtdg &
                    + dHdtt*( dtdbb*dtdg + 2.d0*dtdb*dtdbg ) + dHdft*( dAdbb*dtdg + 2.d0*dAdb*dtdbg ) &
                    + dHdtu*( dudbb*dtdg + 2.d0*dudb*dtdbg ) + dHdt*dtdbbg       
           dHdagg = dHdftt*dAda*dtdg**2 + dHdttt*dtda*dtdg**2 + dHdttu*duda*dtdg**2 + 2.d0*dHdtt*dtdg*dtdag
           dHdbgg = dHdftt*dAdb*dtdg**2 + dHdttt*dtdb*dtdg**2 + dHdttu*dudb*dtdg**2 + 2.d0*dHdtt*dtdg*dtdbg 
           dHdggg = dHdttt*dtdg**3        
                    
!print '(i2,20d19.10)', i,h,dHda,dHdb,dHdg,dHdaa,dHdbb,dHdgg,dHdag,dHdbg,dHdab,dHdaaa,dHdbbb,dHdggg,dHdaab,dHdabb,dHdaag,dHdagg,dHdbbg,dHdbgg,dHdabg
    
           dFc_daa(i)   = 2.d0*euda + d*eudaa + 2.d0*dHda + d*dHdaa 
           dFc_dbb(i)   = 2.d0*eudb + d*eudbb + 2.d0*dHdb + d*dHdbb
           dFc_dgg(i)   = d*dHdgg
           dFc_dab(i)   = euda + eudb + d*eudab + dHda + dHdb + d*dHdab 
           dFc_dag(i)   = dHdg + d*dHdag
           dFc_dbg(i)   = dHdg + d*dHdbg
#ifdef USE_3RD_ORDER_DERIVATIVE
           dFc_daaa(i)  = 3.d0*eudaa + d*eudaaa + 3.d0*dHdaa + d*dHdaaa
           dFc_dbbb(i)  = 3.d0*eudbb + d*eudbbb + 3.d0*dHdbb + d*dHdbbb
           dFc_dggg(i)  = d*dHdggg
           dFc_daab(i)  = eudaa + 2.d0*eudab +d*eudaab + dHdaa + 2.d0*dHdab + d*dHdaab
           dFc_daag(i)  = 2.d0*dHdag + d*dHdaag
           dFc_dabb(i)  = eudbb + 2.d0*eudab + d*eudabb + dHdbb + 2.d0*dHdab + d*dHdabb
           dFc_dbbg(i)  = 2.d0*dHdbg + d*dHdbbg
           dFc_dagg(i)  = dHdgg + d*dHdagg
           dFc_dbgg(i)  = dHdgg + d*dHdbgg
           dFc_dabg(i)  = dHdag + dHdbg + d*dHdabg
#endif                    
!print '(i2,21d19.10)', i,exc(i),d*(eu+h),dF_drho(i,1),dF_drho(i,2),dF_dgradrho(i)*grad_trho(i) &
!                      ,dFc_daa(i),dFc_dbb(i),dFc_dgg(i),dFc_dab(i),dFc_dag(i),dFc_dbg(i),dFc_daaa(i),dFc_dbbb(i) &
!                      ,dFc_dggg(i),dFc_daab(i),dFc_daag(i),dFc_dabb(i),dFc_dbbg(i),dFc_dagg(i),dFc_dbgg(i),dFc_dabg(i)
!print *,dHdag,dHdaag,d,2.d0*dHdag + d*dHdaag
!print *,2.d0*dHdbg + d*dHdbbg
!print *,dHdag + dHdbg + d*dHdabg
               
           
        end do
     end if
     
end subroutine cr_ggapbe_paw_drv3

subroutine ex_ggapw91_paw_drv2(nrc,dnr,nspin,chgrhr_l,grad_rho,exc, &
                                dFx_drho,dFx_dgradrho, &
                                dFx_drr,dFx_drg,dFx_dgg, &
                                dFx_drrr,dFx_drrg,dFx_drgg,dFx_dggg)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,PAI
  implicit none

  integer,intent(in)        :: nrc,dnr,nspin
  real(kind=DP),intent(in)  :: chgrhr_l(nrc,nspin)
  real(kind=DP),intent(in)  :: grad_rho(nrc,nspin)
!  real(kind=DP),intent(in)  :: wos(nrc)  
  real(kind=DP),intent(out) :: exc(nrc)
  real(kind=DP),intent(out) :: dFx_drho(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drr(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dgg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drrr(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drrg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drgg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dggg(nrc,nspin)

  real(kind=DP), parameter :: a1 = 0.19645d0
  real(kind=DP), parameter :: a2 = 0.27430d0
  real(kind=DP), parameter :: a3 = 0.15084d0
  real(kind=DP), parameter :: a4 = 100.0d0

#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter ::  ax = -0.7385587663820224d0
  real(kind=DP), parameter :: a  =  7.7956d0
  real(kind=DP), parameter :: b1 =  0.004d0
  real(kind=DP), parameter :: thrd = 0.33333333333333333333333d0, &
       &                      thrd4 = 1.3333333333333333333333d0
  real(kind=DP), parameter :: thpith = 3.0936677262801d0      ! ???

  real(kind=DP), parameter :: ninth2=0.2222222222222222222d0
  real(kind=DP), parameter :: ninth4=0.4444444444444444444d0
  real(kind=DP), parameter :: ninth8=0.8888888888888888888d0
  real(kind=DP), parameter :: thrd2 =0.6666666666666666666d0
  real(kind=DP), parameter :: eghth =0.125d0
  real(kind=DP), parameter :: twtysvnth8 = 0.29629629629629627985d0
  real(kind=DP), parameter :: twtysvnth32= 1.18518518518518511939d0

#else
  real(kind=DP), parameter :: ax = -0.7385588d0
  real(kind=DP), parameter :: a  =  7.7956d0
  real(kind=DP), parameter :: b1 =  0.004d0
  real(kind=DP), parameter :: thrd   = 0.33333333333d0
  real(kind=DP), parameter :: thrd4  = 1.333333333333333d0
  real(kind=DP), parameter :: thpith = 3.0936677262801d0
  
  real(kind=DP), parameter :: ninth2=0.22222222222d0
  real(kind=DP), parameter :: ninth4=0.44444444444d0
  real(kind=DP), parameter :: ninth8=0.88888888888d0
  real(kind=DP), parameter :: thrd2 =0.66666666667d0
  real(kind=DP), parameter :: eghth =0.125d0
  real(kind=DP), parameter :: twtysvnth8 = 0.296296296296d0
  real(kind=DP), parameter :: twtysvnth32= 1.18518518518519d0
#endif

  real(kind=DP) :: exc1,facw,d,dd,fk,s,fac,s2,s3,s4&
       & ,p0,p1,p2,p3,p4,p5,p6,f,ex,fs,exd,exdd,exc0,excd,excdd

  integer       :: is,i
  
  
  real(kind=DP) :: fss,fsss,ivd,ivdthrd2
  real(kind=DP) :: dfa,ddfa,dddfa
  real(kind=DP) :: fb,ddfb,dddfb
  real(kind=DP) :: p03,p05,a2s2,p32,p33,p34

!---- Spin dependency

  facw = nspin
  exc  = 0.d0
  do is = 1, nspin
     exc1 = 0.d0
     do i = 1,nrc,dnr
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
        p0  = 1.d0 /(dsqrt(1.d0+a*a*s2))
        p1  = dlog(a*s +1.d0/p0)
        p2  = dexp(-a4*s2)
        p3  = 1.d0/(1.d0 +a1*s*p1+b1*s4)
        p4  = 1.d0 + a1*s*p1 + (a2-a3*p2)*s2
        f   = p3*p4
        ex  = fac*f*d
        p5  = 2.d0 *(s*(a2-a3*p2)+a3*a4*s3*p2 - 2.d0*b1*s3)
        p6  = (a1*(p1+a*s*p0)+4*b1*s3)*((a2-a3*p2)*s2-b1*s4)
        fs  = (p5*p3-p6*p3*p3)
        exd = thrd4*fac*(f-s*fs)
        exdd = ax*fs*0.5d0/thpith

        exc0 = ex / facw
        excd = exd
        dFx_drho(i, is) = excd 
        excdd = exdd

! gradient of charge density 95/12/3 H.S.
        if(dabs(grad_rho(i, is)) > 1.d-9) then
!           dFx_dgradrho(i, is) = excdd / grad_rho(i, is)           ! 2010 5
           dFx_dgradrho(i, is) = excdd
        else
           dFx_dgradrho(i, is) = 0.d0
        endif
! gradient of charge density 95/12/3 H.S.
!        exc1 = exc1 + exc0*wos(i)
        exc(i)=exc(i)+exc0
        
! *****  Higher derivatives *****

        p03=p0**3
        p05=p03*p0**2
        a2s2=a**2*s2
        dfa   = a1*(p1+a*s*p0)+4.d0*b1*s3
        ddfa  = a1*a*(2.d0+a2s2)*p03+12.d0*b1*s2
        dddfa = -a1*a**3*s*(4.d0+a2s2)*p05+24.d0*b1*s 
        fb    = (a2-a3*p2)*s2-b1*s4
        ddfb  = 2.d0*(a2-a3*(1.d0-(5.d0-2.d0*a4*s2)*a4*s2)*p2)-12.d0*b1*s2
        dddfb = 4.d0*a3*a4*s*(6.d0+(2.d0*a4*s2-9.d0)*a4*s2)*p2-24.d0*b1*s
       
        p32=p3**2
        p33=p3*p32
        p34=p32*p32       
       
        fss = ddfb*p3-(2.d0*dfa*p5+ddfa*fb)*p32+2.d0*dfa**2*fb*p33
        fsss= dddfb*p3-(3.d0*(dfa*ddfb+ddfa*p5)+dddfa*fb)*p32+ &
                        6.d0*dfa*(dfa*p5+ddfa*fb)*p33-6.d0*dfa**3*fb*p34
                                
        ivd=1.d0/d
        ivdthrd2=ivd**thrd2
        
        if(nspin == 1) then
            dFx_drr(i,is)   =ninth4*ax*ivdthrd2*( f-s*fs+4.d0*fss*s2 )
            dFx_drrr(i,is)  =-twtysvnth8*ax*ivdthrd2*ivd*( f-fs*s+18.d0*fss*s2+8.d0*fsss*s2*s )
            if(dabs(grad_rho(i, is)) > 1.d-9) then
                dFx_drg(i,is)   =-thrd2*ax/thpith*fss*s*ivd
                dFx_dgg(i,is)   =0.25d0*ax/thpith/thpith*fss*ivdthrd2**2
                dFx_drrg(i,is)  =ninth2*ax*ivd*ivd/thpith*( 7.d0*fss*s+4.d0*fsss*s2 )
                dFx_drgg(i,is)  =-thrd*ax/thpith/thpith*ivd*ivdthrd2**2*( fss+fsss*s )
                dFx_dggg(i,is)  =eghth*ax/(thpith**3)*fsss*ivd**2*ivdthrd2
            else
                dFx_drg(i,is)   = 0.d0
                dFx_dgg(i,is)   = 0.d0
                dFx_drrg(i,is)  = 0.d0
                dFx_drgg(i,is)  = 0.d0
                dFx_dggg(i,is)  = 0.d0
            endif
        else if(nspin == 2) then
            dFx_drr(i,is)   =ninth8*ax*ivdthrd2*( f-s*fs+4.d0*fss*s2 )
            dFx_drrr(i,is)  =-twtysvnth32*ax*ivdthrd2*ivd*( f-fs*s+18.d0*fss*s2+8.d0*fsss*s2*s )
            if(dabs(grad_rho(i, is)) > 1.d-9) then
                dFx_drg(i,is)   =-thrd4*ax/thpith*fss*s*ivd
                dFx_dgg(i,is)   =0.5d0*ax/thpith/thpith*fss*ivdthrd2**2
                dFx_drrg(i,is)  =ninth8*ax*ivd*ivd/thpith*( 7.d0*fss*s+4.d0*fsss*s2 )
                dFx_drgg(i,is)  =-thrd4*ax/(thpith**2)*ivd*ivdthrd2**2*( fss+fsss*s )
                dFx_dggg(i,is)  =0.5d0*ax/(thpith**3)*fsss*ivd**2*ivdthrd2
            else
                dFx_drg(i,is)   = 0.d0
                dFx_dgg(i,is)   = 0.d0
                dFx_drrg(i,is)  = 0.d0
                dFx_drgg(i,is)  = 0.d0
                dFx_dggg(i,is)  = 0.d0
            end if
        end if
        
     end do
!     exc(i) = exc(i) + exc1
  end do
end subroutine ex_ggapw91_paw_drv2
! === For nrc decomposion. by takto 2012/12/05 =================================
subroutine ex_ggapw91_paw_drv2_3D(nrc,dnr,nspin,chgrhr_l,grad_rho,exc, &
                                dFx_drho,dFx_dgradrho, &
                                dFx_drr,dFx_drg,dFx_dgg, &
                                dFx_drrr,dFx_drrg,dFx_drgg,dFx_dggg,ista_nrc,iend_nrc,ist,ien)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,PAI
  implicit none

  integer,intent(in)        :: nrc,dnr,nspin
  real(kind=DP),intent(in)  :: chgrhr_l(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(in)  :: grad_rho(ista_nrc:iend_nrc,nspin)
!  real(kind=DP),intent(in)  :: wos(nrc)  
  real(kind=DP),intent(out) :: exc(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFx_drho(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drr(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drg(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dgg(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drrr(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drrg(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drgg(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dggg(ista_nrc:iend_nrc,nspin)

  integer, intent(in)       :: ista_nrc,iend_nrc,ist,ien

  real(kind=DP), parameter :: a1 = 0.19645d0
  real(kind=DP), parameter :: a2 = 0.27430d0
  real(kind=DP), parameter :: a3 = 0.15084d0
  real(kind=DP), parameter :: a4 = 100.0d0

#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter ::  ax = -0.7385587663820224d0
  real(kind=DP), parameter :: a  =  7.7956d0
  real(kind=DP), parameter :: b1 =  0.004d0
  real(kind=DP), parameter :: thrd = 0.33333333333333333333333d0, &
       &                      thrd4 = 1.3333333333333333333333d0
  real(kind=DP), parameter :: thpith = 3.0936677262801d0      ! ???

  real(kind=DP), parameter :: ninth2=0.2222222222222222222d0
  real(kind=DP), parameter :: ninth4=0.4444444444444444444d0
  real(kind=DP), parameter :: ninth8=0.8888888888888888888d0
  real(kind=DP), parameter :: thrd2 =0.6666666666666666666d0
  real(kind=DP), parameter :: eghth =0.125d0
  real(kind=DP), parameter :: twtysvnth8 = 0.29629629629629627985d0
  real(kind=DP), parameter :: twtysvnth32= 1.18518518518518511939d0

#else
  real(kind=DP), parameter :: ax = -0.7385588d0
  real(kind=DP), parameter :: a  =  7.7956d0
  real(kind=DP), parameter :: b1 =  0.004d0
  real(kind=DP), parameter :: thrd   = 0.33333333333d0
  real(kind=DP), parameter :: thrd4  = 1.333333333333333d0
  real(kind=DP), parameter :: thpith = 3.0936677262801d0
  
  real(kind=DP), parameter :: ninth2=0.22222222222d0
  real(kind=DP), parameter :: ninth4=0.44444444444d0
  real(kind=DP), parameter :: ninth8=0.88888888888d0
  real(kind=DP), parameter :: thrd2 =0.66666666667d0
  real(kind=DP), parameter :: eghth =0.125d0
  real(kind=DP), parameter :: twtysvnth8 = 0.296296296296d0
  real(kind=DP), parameter :: twtysvnth32= 1.18518518518519d0
#endif

  real(kind=DP) :: exc1,facw,d,dd,fk,s,fac,s2,s3,s4&
       & ,p0,p1,p2,p3,p4,p5,p6,f,ex,fs,exd,exdd,exc0,excd,excdd

  integer       :: is,i
  
  
  real(kind=DP) :: fss,fsss,ivd,ivdthrd2
  real(kind=DP) :: dfa,ddfa,dddfa
  real(kind=DP) :: fb,ddfb,dddfb
  real(kind=DP) :: p03,p05,a2s2,p32,p33,p34

!---- Spin dependency

  facw = nspin
  exc  = 0.d0
  do is = 1, nspin
     exc1 = 0.d0
     do i = ist,ien,dnr
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

! gradient of charge density 95/12/3 H.S.
        if(dabs(grad_rho(i, is)) > 1.d-9) then
!           dFx_dgradrho(i, is) = excdd / grad_rho(i, is)           ! 2010 5
           dFx_dgradrho(i, is) = excdd
        else
           dFx_dgradrho(i, is) = 0.d0
        endif
! gradient of charge density 95/12/3 H.S.
!        exc1 = exc1 + exc0*wos(i)
        exc(i)=exc(i)+exc0
        
! *****  Higher derivatives *****

        p03=p0**3
        p05=p03*p0**2
        a2s2=a**2*s2
        dfa   = a1*(p1+a*s*p0)+4.d0*b1*s3
        ddfa  = a1*a*(2.d0+a2s2)*p03+12.d0*b1*s2
        dddfa = -a1*a**3*s*(4.d0+a2s2)*p05+24.d0*b1*s 
        fb    = (a2-a3*p2)*s2-b1*s4
        ddfb  = 2.d0*(a2-a3*(1.d0-(5.d0-2.d0*a4*s2)*a4*s2)*p2)-12.d0*b1*s2
        dddfb = 4.d0*a3*a4*s*(6.d0+(2.d0*a4*s2-9.d0)*a4*s2)*p2-24.d0*b1*s
       
        p32=p3**2
        p33=p3*p32
        p34=p32*p32       
       
        fss = ddfb*p3-(2.d0*dfa*p5+ddfa*fb)*p32+2.d0*dfa**2*fb*p33
        fsss= dddfb*p3-(3.d0*(dfa*ddfb+ddfa*p5)+dddfa*fb)*p32+ &
                        6.d0*dfa*(dfa*p5+ddfa*fb)*p33-6.d0*dfa**3*fb*p34
                                
        ivd=1.d0/d
        ivdthrd2=ivd**thrd2
        
        if(nspin == 1) then
            dFx_drr(i,is)   =ninth4*ax*ivdthrd2*( f-s*fs+4.d0*fss*s2 )
            dFx_drrr(i,is)  =-twtysvnth8*ax*ivdthrd2*ivd*( f-fs*s+18.d0*fss*s2+8.d0*fsss*s2*s )
            if(dabs(grad_rho(i, is)) > 1.d-9) then
                dFx_drg(i,is)   =-thrd2*ax/thpith*fss*s*ivd
                dFx_dgg(i,is)   =0.25d0*ax/thpith/thpith*fss*ivdthrd2**2
                dFx_drrg(i,is)  =ninth2*ax*ivd*ivd/thpith*( 7.d0*fss*s+4.d0*fsss*s2 )
                dFx_drgg(i,is)  =-thrd*ax/thpith/thpith*ivd*ivdthrd2**2*( fss+fsss*s )
                dFx_dggg(i,is)  =eghth*ax/(thpith**3)*fsss*ivd**2*ivdthrd2
            else
                dFx_drg(i,is)   = 0.d0
                dFx_dgg(i,is)   = 0.d0
                dFx_drrg(i,is)  = 0.d0
                dFx_drgg(i,is)  = 0.d0
                dFx_dggg(i,is)  = 0.d0
            endif
        else if(nspin == 2) then
            dFx_drr(i,is)   =ninth8*ax*ivdthrd2*( f-s*fs+4.d0*fss*s2 )
            dFx_drrr(i,is)  =-twtysvnth32*ax*ivdthrd2*ivd*( f-fs*s+18.d0*fss*s2+8.d0*fsss*s2*s )
            if(dabs(grad_rho(i, is)) > 1.d-9) then
                dFx_drg(i,is)   =-thrd4*ax/thpith*fss*s*ivd
                dFx_dgg(i,is)   =0.5d0*ax/thpith/thpith*fss*ivdthrd2**2
                dFx_drrg(i,is)  =ninth8*ax*ivd*ivd/thpith*( 7.d0*fss*s+4.d0*fsss*s2 )
                dFx_drgg(i,is)  =-thrd4*ax/(thpith**2)*ivd*ivdthrd2**2*( fss+fsss*s )
                dFx_dggg(i,is)  =0.5d0*ax/(thpith**3)*fsss*ivd**2*ivdthrd2
            else
                dFx_drg(i,is)   = 0.d0
                dFx_dgg(i,is)   = 0.d0
                dFx_drrg(i,is)  = 0.d0
                dFx_drgg(i,is)  = 0.d0
                dFx_dggg(i,is)  = 0.d0
            end if
        end if
        
     end do
!     exc(i) = exc(i) + exc1
  end do
end subroutine ex_ggapw91_paw_drv2_3D
! ==============================================================================

subroutine ex_ggapw91_paw_drv3(nrc,num_ir,irs,nspin,chgrhr_l,grad_rho,exc, &
                                dFx_drho,dFx_dgradrho, &
                                dFx_drr,dFx_drg,dFx_dgg, &
                                dFx_drrr,dFx_drrg,dFx_drgg,dFx_dggg)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : DP,PAI
  implicit none

  integer,intent(in)        :: nrc,num_ir,irs(num_ir),nspin
  real(kind=DP),intent(in)  :: chgrhr_l(nrc,nspin)
  real(kind=DP),intent(in)  :: grad_rho(nrc,nspin)
!  real(kind=DP),intent(in)  :: wos(nrc)  
  real(kind=DP),intent(out) :: exc(nrc)
  real(kind=DP),intent(out) :: dFx_drho(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drr(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dgg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drrr(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drrg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drgg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dggg(nrc,nspin)

  real(kind=DP), parameter :: a1 = 0.19645d0
  real(kind=DP), parameter :: a2 = 0.27430d0
  real(kind=DP), parameter :: a3 = 0.15084d0
  real(kind=DP), parameter :: a4 = 100.0d0
  real(kind=DP), parameter :: ax = -0.7385588d0
  real(kind=DP), parameter :: a  =  7.7956d0
  real(kind=DP), parameter :: b1 =  0.004d0
  real(kind=DP), parameter :: thrd   = 0.33333333333d0
  real(kind=DP), parameter :: thrd4  = 1.333333333333333d0
  real(kind=DP), parameter :: thpith = 3.0936677262801d0
  
  real(kind=DP), parameter :: ninth2=0.22222222222d0
  real(kind=DP), parameter :: ninth4=0.44444444444d0
  real(kind=DP), parameter :: ninth8=0.88888888888d0
  real(kind=DP), parameter :: thrd2 =0.66666666667d0
  real(kind=DP), parameter :: eghth =0.125d0
  real(kind=DP), parameter :: twtysvnth8 = 0.296296296296d0
  real(kind=DP), parameter :: twtysvnth32= 1.18518518518519d0

  real(kind=DP) :: exc1,facw,d,dd,fk,s,fac,s2,s3,s4&
       & ,p0,p1,p2,p3,p4,p5,p6,f,ex,fs,exd,exdd,exc0,excd,excdd

  integer       :: is,i,nr
  
  
  real(kind=DP) :: fss,fsss,ivd,ivdthrd2
  real(kind=DP) :: dfa,ddfa,dddfa
  real(kind=DP) :: fb,ddfb,dddfb
  real(kind=DP) :: p03,p05,a2s2,p32,p33,p34

!---- Spin dependency

  facw = nspin
  exc  = 0.d0
  do is = 1, nspin
     exc1 = 0.d0
!     do i = 1,nrc,dnr
     do nr=1,num_ir
        i=irs(nr)
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

! gradient of charge density 95/12/3 H.S.
        if(dabs(grad_rho(i, is)) > 1.d-9) then
!           dFx_dgradrho(i, is) = excdd / grad_rho(i, is)           ! 2010 5
           dFx_dgradrho(i, is) = excdd
        else
           dFx_dgradrho(i, is) = 0.d0
        endif
! gradient of charge density 95/12/3 H.S.
!        exc1 = exc1 + exc0*wos(i)
        exc(i)=exc(i)+exc0
        
! *****  Higher derivatives *****

        p03=p0**3
        p05=p03*p0**2
        a2s2=a**2*s2
        dfa   = a1*(p1+a*s*p0)+4.d0*b1*s3
        ddfa  = a1*a*(2.d0+a2s2)*p03+12.d0*b1*s2
        dddfa = -a1*a**3*s*(4.d0+a2s2)*p05+24.d0*b1*s 
        fb    = (a2-a3*p2)*s2-b1*s4
        ddfb  = 2.d0*(a2-a3*(1.d0-(5.d0-2.d0*a4*s2)*a4*s2)*p2)-12.d0*b1*s2
        dddfb = 4.d0*a3*a4*s*(6.d0+(2.d0*a4*s2-9.d0)*a4*s2)*p2-24.d0*b1*s
       
        p32=p3**2
        p33=p3*p32
        p34=p32*p32       
       
        fss = ddfb*p3-(2.d0*dfa*p5+ddfa*fb)*p32+2.d0*dfa**2*fb*p33
        fsss= dddfb*p3-(3.d0*(dfa*ddfb+ddfa*p5)+dddfa*fb)*p32+ &
                        6.d0*dfa*(dfa*p5+ddfa*fb)*p33-6.d0*dfa**3*fb*p34
                                
        ivd=1.d0/d
        ivdthrd2=ivd**thrd2
        
        if(nspin == 1) then
            dFx_drr(i,is)   =ninth4*ax*ivdthrd2*( f-s*fs+4.d0*fss*s2 )
            dFx_drrr(i,is)  =-twtysvnth8*ax*ivdthrd2*ivd*( f-fs*s+18.d0*fss*s2+8.d0*fsss*s2*s )
            if(dabs(grad_rho(i, is)) > 1.d-9) then
                dFx_drg(i,is)   =-thrd2*ax/thpith*fss*s*ivd
                dFx_dgg(i,is)   =0.25d0*ax/thpith/thpith*fss*ivdthrd2**2
                dFx_drrg(i,is)  =ninth2*ax*ivd*ivd/thpith*( 7.d0*fss*s+4.d0*fsss*s2 )
                dFx_drgg(i,is)  =-thrd*ax/thpith/thpith*ivd*ivdthrd2**2*( fss+fsss*s )
                dFx_dggg(i,is)  =eghth*ax/(thpith**3)*fsss*ivd**2*ivdthrd2
            else
                dFx_drg(i,is)   = 0.d0
                dFx_dgg(i,is)   = 0.d0
                dFx_drrg(i,is)  = 0.d0
                dFx_drgg(i,is)  = 0.d0
                dFx_dggg(i,is)  = 0.d0
            endif
        else if(nspin == 2) then
            dFx_drr(i,is)   =ninth8*ax*ivdthrd2*( f-s*fs+4.d0*fss*s2 )
            dFx_drrr(i,is)  =-twtysvnth32*ax*ivdthrd2*ivd*( f-fs*s+18.d0*fss*s2+8.d0*fsss*s2*s )
            if(dabs(grad_rho(i, is)) > 1.d-9) then
                dFx_drg(i,is)   =-thrd4*ax/thpith*fss*s*ivd
                dFx_dgg(i,is)   =0.5d0*ax/thpith/thpith*fss*ivdthrd2**2
                dFx_drrg(i,is)  =ninth8*ax*ivd*ivd/thpith*( 7.d0*fss*s+4.d0*fsss*s2 )
                dFx_drgg(i,is)  =-thrd4*ax/(thpith**2)*ivd*ivdthrd2**2*( fss+fsss*s )
                dFx_dggg(i,is)  =0.5d0*ax/(thpith**3)*fsss*ivd**2*ivdthrd2
            else
                dFx_drg(i,is)   = 0.d0
                dFx_dgg(i,is)   = 0.d0
                dFx_drrg(i,is)  = 0.d0
                dFx_drgg(i,is)  = 0.d0
                dFx_dggg(i,is)  = 0.d0
            end if
        end if
        
     end do
!     exc(i) = exc(i) + exc1
  end do
end subroutine ex_ggapw91_paw_drv3

subroutine cr_ggapw91_paw_drv2(nrc,dnr,nspin,chgrhr_l,grad_trho,exc,dF_drho,dF_dgradrho &
                                ,dFc_daa,dFc_dbb,dFc_dgg,dFc_dab,dFc_dag,dFc_dbg &
                                ,dFc_daaa,dFc_dbbb,dFc_dggg &
                                ,dFc_daab,dFc_daag,dFc_dabb,dFc_dbbg,dFc_dagg,dFc_dbgg,dFc_dabg)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : PAI,DP
  implicit none

  integer,intent(in)        :: nrc,dnr,nspin
  real(kind=DP),intent(in)  :: chgrhr_l(nrc,nspin)
  real(kind=DP),intent(in) :: grad_trho(nrc)
  real(kind=DP),intent(inout) :: exc(nrc)
  real(kind=DP),intent(out) :: dF_drho(nrc,nspin)
  real(kind=DP),intent(out) :: dF_dgradrho(nrc)
  real(kind=DP),intent(out) :: dFc_daa(nrc)
  real(kind=DP),intent(out) :: dFc_dbb(nrc)
  real(kind=DP),intent(out) :: dFc_dgg(nrc)
  real(kind=DP),intent(out) :: dFc_dab(nrc)
  real(kind=DP),intent(out) :: dFc_dag(nrc)
  real(kind=DP),intent(out) :: dFc_dbg(nrc)
  real(kind=DP),intent(out) :: dFc_daaa(nrc)
  real(kind=DP),intent(out) :: dFc_dbbb(nrc)
  real(kind=DP),intent(out) :: dFc_dggg(nrc)
  real(kind=DP),intent(out) :: dFc_daab(nrc)
  real(kind=DP),intent(out) :: dFc_daag(nrc)
  real(kind=DP),intent(out) :: dFc_dabb(nrc)
  real(kind=DP),intent(out) :: dFc_dbbg(nrc)
  real(kind=DP),intent(out) :: dFc_dagg(nrc)
  real(kind=DP),intent(out) :: dFc_dbgg(nrc)
  real(kind=DP),intent(out) :: dFc_dabg(nrc)

  real(kind=DP), parameter :: zeta_minimum    = 1.d-4
  real(kind=DP), parameter :: zeta_minimum2   = 1.d-6
  real(kind=DP), parameter :: density_minimum = 1.d-20
  real(kind=DP), parameter :: density_minimum2= 1.d-5

  real(kind=DP), parameter  :: a  = 0.0310907d0,    a1   = 0.21370d0
  real(kind=DP), parameter  :: b1 = 7.5957d0,       b2n  = 3.5876d0
  real(kind=DP), parameter  :: b3 = 1.6382d0,       b4   = 0.49294d0
  real(kind=DP), parameter  :: p  = 1.00d0,         p1   = p+1.d0

#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: ap  = 0.01554535d0,    a1p  = 0.20548d0
#else
  real(kind=DP), parameter  :: ap  = 0.015545d0,    a1p  = 0.20548d0
#endif

  real(kind=DP), parameter  :: b1p = 14.1189d0,     b2np = 6.1977d0
  real(kind=DP), parameter  :: b3p = 3.3662d0,      b4p  = 0.62517d0

! vwn 95/12/2 Y.M
#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: aq  = 0.0168869d0,    a1q  = 0.11125d0
#else
  real(kind=DP), parameter  :: aq  = 0.016887d0,    a1q  = 0.11125d0
#endif

  real(kind=DP), parameter  :: b1q = 10.3570d0,     b2nq = 3.6231d0
  real(kind=DP), parameter  :: b3q = 0.88026d0,     b4q  = 0.49671d0

! vwn 95/12/2 Y.M
#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: gamma = 0.03109069086965489503494086371273d0
#else
  real(kind=DP), parameter  :: gamma = 0.031091
#endif

!--------------------------
  real(kind=DP), parameter :: cx  = -0.001667212d0, alf =  0.09d0
  real(kind=DP), parameter :: c1  =  0.002568d0,    c2  =  0.023266d0
  real(kind=DP), parameter :: c3  =  7.389d-6,      c4  =  8.723d0
  real(kind=DP), parameter :: c5  =  0.472d0,       c6  =  7.389d-2
  real(kind=DP), parameter :: a4  = 100.d0

  real(kind=DP), parameter  :: xnu = 15.75592d0, cc0 = 0.004235d0

#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: thrd = 0.33333333333333330,&
       &                       sixth7 = 1.1666666666666666d0
  real(kind=DP), parameter  :: thrd7= 2.333333333333333333d0
  real(kind=DP), parameter  :: sixth13= 2.166666666666666666d0
  real(kind=DP), parameter  :: sixth19= 3.166666666666666666d0
  real(kind=DP), parameter  :: ninth= 0.11111111111111111111d0

#else
  real(kind=DP), parameter  :: thrd = 0.333333333333d0, sixth7 = 1.1666666666666666d0
  real(kind=DP), parameter  :: thrd7= 2.333333333333d0
  real(kind=DP), parameter  :: sixth13= 2.1666666666666666d0
  real(kind=DP), parameter  :: sixth19= 3.1666666666666666d0
  real(kind=DP), parameter  :: ninth= 0.111111111111111d0
#endif

  real(kind=DP), parameter  :: twntysvnth=  3.703703703703703d-2
  real(kind=DP), parameter  :: fftyfrth= 1.851851851851852d-2
  real(kind=DP), parameter  :: eityonth= 1.234567901234568d-2

  integer       :: is,i
  real(kind=DP) :: facw,bet,delt,exc1, g,g3,g4,facpon,d,dd,rs,fk,sk,t,s,q0 &
       &         , rs12,rs32,rsp,q1,q2,eu,q3,eurs,ecd,pon &
       &         , b,b2,t2,s2,t4,t6,rs2,rs3,q4,q5,q6,q7,cc,r0,r1,coeff,r2,r3,a4ms,h0,h1,h,q8 &
       &         , h0t,h0b,h0rs,h1t,ccrs,r1rs,h1rs,ht,hrs &
       &         , ec1,ec1d,ec1dd,exc0,excd,excdd &
       &         , thrd2,thrd4,fzdd,zeta,onpzeta,onmzeta,fz,zetadxd &
       &         , q0p,q1p,q2p,q3p,q0q,q1q,q2q,q3q,fzd,onzt4,gzt,bzt,h0zt,h1zt,euzt &
       &         , zetadxda,zetadxdb,ecda,ecdb,ec1da,ec1db,excda,excdb
  real(kind=DP) :: ccrs2,ccrs3,dq6,ddq6,dq7,ddq7,dddq7
  real(kind=DP) :: iq7,iq72,iq73,iq74
       
  real(kind=DP) :: rsdn,rsdnn,rsdnnn
  real(kind=DP) :: q1drr,q1drrr
  real(kind=DP) :: eudrr,eudrrr,invq12pq1,q12p1
  real(kind=DP) :: eudn,eudnn,eudnnn,ecdd,ecddd
  real(kind=DP) :: dtdn,dtdg,dtdnn,dtdgg,dtdng,dtdnnn,dtdnng,dtdngg,dtdggg
  real(kind=DP) :: dAde,dAdee,dAdeee
  real(kind=DP) :: dAdn,dAdnn,dAdnnn
  real(kind=DP) :: fE,fEdt,fEdtt,fEda,fEdaa
  real(kind=DP) :: fF,fFdt,fFdtt,fFda,fFdaa
  real(kind=DP) :: fG,fGdt,fGdtt,fGda,fGdaa
  real(kind=DP) :: at2,a2t4,a3t6,t3,t5,t7,t8
  real(kind=DP) :: dHda,dHdt,dHdaa,dHdat,dHdtt,dHdaaa,dHdaat,dHdatt,dHdttt,betivE,betivE2
  real(kind=DP) :: dHdn,dHdg,dHdnn,dHdng,dHdgg,dHdnnn,dHdnng,dHdngg,dHdggg
  real(kind=DP) :: zetada,zetadb,zetadaa,zetadab,zetadbb,zetadaaa,zetadaab,zetadabb,zetadbbb
  real(kind=DP) :: id2,id3,id4
  real(kind=DP) :: e1dr,e1drr,e1drrr
  real(kind=DP) :: e1pdr,e1pdrr,e1pdrrr
  real(kind=DP) :: e1qdr,e1qdrr,e1qdrrr
  real(kind=DP) :: q1pdrr,q1pdrrr,invq1p2pq1p,q1p2p1
  real(kind=DP) :: q1qdrr,q1qdrrr,invq1q2pq1q,q1q2p1
  real(kind=DP) :: fzdd0,fzddd,opzthrdm2,omzthrdm2
  real(kind=DP) :: zeta2,zeta3,zeta4,eudrz,eudzz,eudrrz,eudrzz,eudzzz
  real(kind=DP) :: e1,e2,e3
  real(kind=DP) :: euda,eudb,eudaa,eudbb,eudab,eudaaa,eudaab,eudabb,eudbbb
  real(kind=DP) :: rsdn2,rsdn3,zetada2,zetada3,zetadb2,zetadb3
  real(kind=DP) :: gzz,gzzz,duda,dudb,dudaa,dudab,dudbb,dudaaa,dudaab,dudabb,dudbbb
  real(kind=DP) :: dtdu,dtduu,dtduuu,dtdnu,dtdnnu,dtdnuu,dtdug,dtduug,dtdugg,dtdnug,invg
  real(kind=DP) :: dtda,dtdb,dtdaa,dtdab,dtdag,dtdbg,dtdbb
  real(kind=DP) :: dtdaaa,dtdaab,dtdabb,dtdbbb,dtdaag,dtdagg,dtdbbg,dtdbgg,dtdabg
  real(kind=DP) :: invg3,m3u
  real(kind=DP) :: dAdu,dAduu,dAduuu,dAdeu,dAdeeu,dAdeuu
  real(kind=DP) :: dAda,dAdb,dAdaa,dAdab,dAdbb,dAdaaa,dAdaab,dAdabb,dAdbbb
  real(kind=DP) :: dHdf,dHdff,dHdft,dHdfff,dHdfft,dHdftt
  real(kind=DP) :: dHdu,dHduu,dHduuu,dHdfu,dHdtu,dHdffu,dHdfuu,dHdttu,dHdtuu,dHdftu
  real(kind=DP) :: tiu,siu2
  real(kind=DP) :: dHdb,dHdbb,dHdag,dHdbg,dHdab,dHdbbb,dHdaab,dHdabb,dHdbbg,dHdbgg,dHdaag,dHdagg,dHdabg
  
  real(kind=DP) :: dH1dr,dH1dt,dH1drr,dH1drt,dH1dtt,dH1drrr,dH1drrt,dH1drtt,dH1dttt,dH1drtu
  real(kind=DP) :: dH1du,dH1duu,dH1dru,dH1dtu,dH1duuu,dH1drru,dH1druu,dH1dtuu,dH1dttu
  real(kind=DP) :: r1rs2,r1rs3,r1rs4,ivt,ivt2,ivtmr1t,ivu,ivu2,r1t2,r1t22,r1t23
  real(kind=DP) :: dH1dn,dH1dg,dH1dnn,dH1dng,dH1dgg,dH1dnnn,dH1dnng,dH1dngg,dH1dggg
  real(kind=DP) :: dH1da,dH1db,dH1daa,dH1dab,dH1dbb,dH1dbg,dH1dag
  real(kind=DP) :: dH1daaa,dH1daab,dH1dabb,dH1dbbb,dH1dbbg,dH1dbgg,dH1dagg,dH1daag,dH1dabg

#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: eta = 1.d-12
  real(kind=DP), parameter  :: sixthm = -0.1666666666666666666666d0
#endif

  facw = nspin

#ifdef GGA_PAW_REVISED
  bet=0.06672455060314922d0
#else
  bet = xnu*cc0
#endif

  delt = 2.d0*alf/bet
  
     if ( nspin == 1 ) then
        g   = 1.d0
        g3  = g**3
        g4  = g3*g
        facpon = -delt/(g3*bet)

        do i = 1,nrc,dnr ! MPI
           d  = facw*chgrhr_l(i, 1)
           if(d < density_minimum) cycle
! gradient of charge density 95/12/3 H.S.
           dd = facw* grad_trho(i)
           rs = (0.75d0/(PAI*d))**thrd
           fk = (3.d0*PAI*PAI*d)**thrd
           sk = dsqrt(4.d0*fk/PAI)
           if(d > density_minimum2) then
              t = dd/(d*sk*2.d0)
              s = dd/(d*fk*2.d0)
           else
              t = 0.d0
              s = 0.d0
           endif

           q0 = -2.d0 *a*( 1.d0 + a1*rs)
           rs12 = dsqrt(rs)
           rs32 = rs12**3.d0
           rsp = rs**p
           q1 = 2.d0 *a*(b1*rs12+b2n*rs+b3*rs32+b4*rs*rsp)
           q2 = log( 1.d0 + 1.d0/q1)
           eu = q0*q2
           q3 = a*( b1/rs12 +2.d0*b2n +3.d0*b3*rs12 +2.d0*b4*p1*rsp )
           eurs = -2.d0 *a*a1*q2 - q0*q3/( q1**2 + q1)
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
           q4  = 1.d0 + b*t2
           q5  = 1.d0 + b*t2 +b2*t4
           q6  = c1+c2*rs+c3*rs2
           q7  = 1.d0+c4*rs+c5*rs2+c6*rs3
           cc  = -cx + q6/q7
           r0  = (sk/fk)**2.d0
           r1  = a4*r0*g4
           coeff = cc -cc0 -3.d0*cx /7.d0
           r2  = xnu*coeff*g3
           a4ms = -a4*g4*s2
           r3  = dexp(a4ms)
           h0  = g3*(bet/delt) *dlog(1.d0 +delt*q4*t2/q5 )
           h1  = r3*r2*t2
           h   = h0+h1
!======================================================
           q8  = q5*q5+delt*q4*q5*t2
           h0t = 2.d0 *bet*t*(1.d0+2.d0*b*t2)/q8 * g3
           h0b = -bet*t6*(2.d0*b+b2*t2)/q8 * g3
           h0rs = h0b*b*eurs*(b+delt)/ bet / g3
           h1t  = 2.d0*r3*r2*t*(1.d0-r1*t2)
           ccrs = (c2+2*c3*rs)/q7-q6*(c4+2*c5*rs+3*c6*rs2)/q7**2
           r1rs = 100.d0 *r0/rs
           h1rs = xnu*t2*r3*(ccrs - coeff*t2*r1rs) * g3
           ht   = h0t+h1t
           hrs  = h0rs + h1rs

           ec1  = d*h/facw
           ec1d = h-thrd*rs*hrs-sixth7*t*ht
           ec1dd = 0.5d0*ht/sk

           exc0 = eu*d/facw + ec1
           excd = ecd + ec1d
           dF_drho(i, 1) = dF_drho(i, 1) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(dabs(grad_trho(i)) > 1.0d-9) then
!              dF_dgradrho(i) = excdd / grad_trho(i)
              dF_dgradrho(i) = excdd
           else
              dF_dgradrho(i) = 0.d0
           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d/facw
!           exc1 = exc1 + exc0*wos(i)
           exc(i)=exc(i)+exc0
        
! ***** Calc heigher derivatives. *****
           
           rsdn=-thrd*rs/d
           rsdnn=-4.d0*thrd*rsdn/d
           rsdnnn=-thrd7*rsdnn/d
!print '(i2,5e19.10)',i,d,rs,rsdn,rsdnn,rsdnnn
!q3 = a*(b1/rs12+2*b2n+3*b3*rs12+2*b4*p1*rsp)   
           q1drr =0.5d0*a*( -b1/rs32+3.d0*b3/rs12+8.d0*b4 )
           q1drrr=0.75d0*a*( b1/rs - b3 )/rs32
!eurs = -2*a*a1*q2 - q0*q3/(q1**2 + q1)
           invq12pq1=1.d0/(q1**2 + q1)
           q12p1=2.d0*q1+1.d0
           eudrr=4.d0*a*a1*q3*invq12pq1-q0*( q1drr-q12p1*q3**2*invq12pq1 )*invq12pq1
           eudrrr=6.d0*a*a1*( q1drr - q12p1*q3**2*invq12pq1 ) * invq12pq1 &
                    -q0*( q1drrr -(2.d0*q3**3+3.d0*q12p1*q3*q1drr -  &
                    2.d0*q12p1**2*q3**3* invq12pq1)*invq12pq1 ) * invq12pq1
           eudn=eurs*rsdn
           eudnn=eudrr*rsdn**2+eurs*rsdnn
           eudnnn=eudrrr*rsdn**3+3.d0*eudrr*rsdn*rsdnn+eurs*rsdnnn
!print '(i2,5e19.12)',i,d,eu,eudn,eudnn,eudnnn
           ecdd=2.d0*eudn+d*eudnn
           ecddd=3.d0*eudnn+d*eudnnn
!print '(i2,5e19.12)',i,d,eu*d,ecdd,ecddd
           
           dtdn  = -sixth7*t/d
           dtdnn = -sixth13*dtdn/d
           dtdnnn= -sixth19*dtdnn/d
!print '(i2,5e19.12)',i,d,t,dtdn,dtdnn,dtdnnn
           if(dabs(grad_trho(i)) > 1.0d-9) then
               dtdg  = t/dd
               dtdng = -sixth7*dtdg/d
               dtdnng= -sixth13*dtdng/d
           else
               dtdg  = 0.d0
               dtdng = 0.d0
               dtdnng= 0.d0
           endif
           dtdngg = 0.d0
           dtdggg = 0.d0
           
           dAde   = b * ( b + delt ) / bet
           dAdee  = 2.d0 * dAde * ( b + 0.5d0 * delt ) / bet
           dAdeee = 2.d0 * ( dAde**2 + ( b + 0.5d0 * delt )*dAdee ) /bet
           
           dAdn   = dAde  *eudn
           dAdnn  = dAdee *eudn**2  + dAde*eudnn
           dAdnnn = dAdeee*eudn**3  + 3.d0*dAdee*eudn*eudnn + dAde*eudnnn   
!print '(i2,5e19.10)' ,i,d,b,dAdnn
           
           at2  = b * t2
           a2t4 = at2**2
           a3t6 = at2*a2t4
           t3   = t2*t
           t5   = t4*t
           t7   = t6*t
           t8   = t7*t
           fE    = q5**2 + delt*t2*q4*q5
           fEdt  = 4.d0*b*t*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt *2.d0*t*(1.d0 + 4.d0*at2 + 6.d0*a2t4 + 4.d0*a3t6 )
           fEdtt = 4.d0*b*( 1.d0 + 9.d0*at2 + 15.d0*a2t4 + 14.d0*a3t6 ) + &
                    delt*( 2.d0 + 24.d0*at2 + 60.d0*a2t4 + 56.d0*a3t6 )
           fEda  = 2.d0*t2*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt*t4*( 2.d0 + 4.d0*at2 + 3.d0*a2t4 )
           fEdaa = 6.d0*t4*( 2.d0*q5 - 1.d0) + &
                    delt*2.d0*t6*( 2.d0 + 3.d0*at2 )
           fF    = t*( 1.d0 + 2.d0*at2 )         
           fFdt  = 1.d0 + 6.d0*at2
           fFdtt = 12.d0*b*t
           fFda  = 2.d0*t3
           fFdaa = 0.d0
           fG    = 2.d0*b*t6 + b2*t8
           fGdt  = 12.d0*b*t5 + 8.d0*b2*t7
           fGdtt = 60.d0*b*t4 + 56.d0*b2*t6
           fGda  = 2.d0*t6 + 2.d0*b*t8
           fGdaa = 2.d0*t8
           
           betivE=bet/fE
           betivE2=betivE/fE
           dHda   =     -fG*betivE
           dHdt   = 2.d0*fF*betivE
           dHdaa  =     -( fGda*fE-fG*fEda )*betivE2
           dHdat  =     -( fGdt*fE-fG*fEdt )*betivE2
           dHdtt  = 2.d0*( fFdt*fE-fF*fEdt )*betivE2
           dHdaaa =     -( fGdaa*fE-fG*fEdaa-2.d0*( fGda*fE-fG*fEda )*fEda/fE )*betivE2
           dHdaat = 2.d0*( fFdaa*fE-fF*fEdaa-2.d0*( fFda*fE-fF*fEda )*fEda/fE )*betivE2
           dHdatt =     -( fGdtt*fE-fG*fEdtt-2.d0*( fGdt*fE-fG*fEdt )*fEdt/fE )*betivE2
           dHdttt = 2.d0*( fFdtt*fE-fF*fEdtt-2.d0*( fFdt*fE-fF*fEdt )*fEdt/fE )*betivE2
           
           dHdn     = dHda  *dAdn       + dHdt*dtdn
           dHdg     = dHdt  *dtdg
           dHdnn    = dHdaa *dAdn**2    + 2.d0*dHdat*dAdn*dtdn  + dHdtt*dtdn**2   + dHda*dAdnn  + dHdt*dtdnn
           dHdng    = dHdat *dAdn*dtdg  + dHdtt*dtdn*dtdg       + dHdt*dtdng
           dHdgg    = dHdtt*dtdg**2
           dHdnnn   = dHdaaa*dAdn**3    + dHdttt*dtdn**3 &
                        + 3.d0*dHdaat*dAdn**2*dtdn  + 3.d0*dHdatt*dAdn*dtdn**2 &
                        + 3.d0*dHdaa*dAdn*dAdnn + 3.d0*dHdtt*dtdn*dtdnn + 3.d0*dHdat*(dAdnn*dtdn + dAdn*dtdnn) &
                        + dHda*dAdnnn   + dHdt*dtdnnn
           dHdnng   = dHdaat*dAdn**2*dtdg   + 2.d0*dHdatt*dAdn*dtdn*dtdg    + dHdttt*dtdn**2*dtdg &
                        + dHdat*(dAdnn*dtdg + 2.d0*dAdn*dtdng)  + dHdtt*(dtdnn*dtdg + 2.d0*dtdn*dtdng) &
                        + dHdt*dtdnng
           dHdngg   = dHdatt*dAdn*dtdg**2   + dHdttt*dtdn*dtdg**2   + 2.d0*dHdtt*dtdg*dtdng
           dHdggg   = dHdttt*dtdg**3
           
           dq6   = c2+2.d0*c3*rs
           ddq6  = 2.d0*c3
           dq7   = c4+2.d0*c5*rs+3.d0*c6*rs2
           ddq7  = 2.d0*c5+6.d0*c6*rs
           dddq7 = 6.d0*c6
           
           iq7=1.d0/q7
           iq72=iq7*iq7
           iq73=iq72*iq7
           iq74=iq73*iq7
           
           ccrs2 =  ddq6*iq7-(2.d0*dq6*dq7+q6*ddq7)*iq72+2.d0*q6*dq7**2*iq73
           ccrs3 =  -iq72*(3.d0*(ddq6*dq7+dq6*ddq7)+q6*dddq7) + &
                                    6.d0*dq7*(dq6*dq7+q6*ddq7)*iq73 - &
                                    6.d0*q6*dq7**3*iq74
           r1rs2=r1rs*r1rs
           r1rs3=r1rs2*r1rs
           r1rs4=r1rs3*r1rs
           
           if(dabs(grad_trho(i)) > 1.0d-9) then
               ivt=1.d0/t
               ivt2=ivt**2
               ivtmr1t=ivt-r1*t
           else
               ivt=0.d0
               ivt2=0.d0
               ivtmr1t=0.d0
           end if
           
           dH1dr=h1rs
           dH1dt=h1t                     
           dH1drr  = xnu*ccrs2*t2*r3-2.d0*h1rs*r1rs*t2-h1*r1rs2*t4
           dH1drt  = 2.d0*h1rs*ivtmr1t-2.d0*h1*t*r1rs
           dH1dtt  = 2.d0*h1t*ivtmr1t-2.d0*h1*(ivt2+r1)
           dH1drrr = xnu*ccrs3*t2*r3-3.d0*t2*r1rs*(dH1drr+h1rs*r1rs*t2)-h1*r1rs3*t6
           dH1drrt = 2.d0*dH1drr*ivtmr1t-4.d0*h1rs*r1rs*t
           dH1drtt = 2.d0*dH1drt*ivtmr1t-2.d0*h1t*r1rs*t-2.d0*ivt*h1rs*(ivt+r1*t)-2.d0*h1*r1rs
           dH1dttt = 2.d0*dH1dtt*ivtmr1t-4.d0*h1t*(ivt2+r1)+4.d0*h1*ivt2*ivt
           
           dH1dn     = dH1dr  *rsdn       + dH1dt*dtdn
           dH1dg     = dH1dt  *dtdg
           dH1dnn    = dH1drr *rsdn**2    + 2.d0*dH1drt*rsdn*dtdn  + dH1dtt*dtdn**2   + dH1dr*rsdnn  + dH1dt*dtdnn
           dH1dng    = dH1drt *rsdn*dtdg  + dH1dtt*dtdn*dtdg       + dH1dt*dtdng
           dH1dgg    = dH1dtt*dtdg**2
           dH1dnnn   = dH1drrr*rsdn**3    + dH1dttt*dtdn**3 &
                        + 3.d0*dH1drrt*rsdn**2*dtdn  + 3.d0*dH1drtt*rsdn*dtdn**2 &
                        + 3.d0*dH1drr*rsdn*rsdnn + 3.d0*dH1dtt*dtdn*dtdnn + 3.d0*dH1drt*(rsdnn*dtdn + rsdn*dtdnn) &
                        + dH1dr*rsdnnn   + dH1dt*dtdnnn
           dH1dnng   = dH1drrt*rsdn**2*dtdg   + 2.d0*dH1drtt*rsdn*dtdn*dtdg    + dH1dttt*dtdn**2*dtdg &
                        + dH1drt*(rsdnn*dtdg + 2.d0*rsdn*dtdng)  + dH1dtt*(dtdnn*dtdg + 2.d0*dtdn*dtdng) &
                        + dH1dt*dtdnng
           dH1dngg   = dH1drtt*rsdn*dtdg**2   + dH1dttt*dtdn*dtdg**2   + 2.d0*dH1dtt*dtdg*dtdng
           dH1dggg   = dH1dttt*dtdg**3
           
!dF_drho(i, 1)=h+d*dHdn+ecd
           dFc_daa(i)   = 2.d0*(dHdn + dH1dn) + d*(dHdnn + dH1dnn)       + ecdd
           dFc_dag(i)   = dHdg + dH1dg + d*(dHdng + dH1dng)
           dFc_dgg(i)   = d*(dHdgg + dH1dgg)
           dFc_daaa(i)  = 3.d0*(dHdnn + dH1dnn)    + d*(dHdnnn + dH1dnnn)  + ecddd
           dFc_daag(i)  = 2.d0*(dHdng + dH1dng)    + d*(dHdnng + dH1dnng)
           dFc_dagg(i)  = dHdgg + dH1dgg        + d*(dHdngg + dH1dngg)
           dFc_dggg(i)  = d*(dHdggg + dH1dggg)  
        end do

     else if ( nspin ==  2 ) then
        thrd2 = thrd * 2.d0
        thrd4 = thrd * 4.d0
        fzdd = 8.d0 /(9.d0*(2.d0**thrd4 - 2.d0 ))

        do i = 1,nrc,dnr
           d = chgrhr_l(i, 1) + chgrhr_l(i, nspin)
           if(d < density_minimum ) cycle
           zeta = ( chgrhr_l(i, 1) - chgrhr_l(i, 2) ) / d
           onpzeta = 2.d0 *chgrhr_l(i, 1) / d
           onmzeta = 2.d0 *chgrhr_l(i, 2) / d
           g = ( onpzeta**thrd2 + onmzeta**thrd2 )/2.d0
           if(dabs(zeta) < zeta_minimum ) then
              fz= 4.d0/9.d0 *zeta*zeta*(1.d0 + 5.d0/54.d0 *zeta*zeta &
                   &               *(1.d0 +44.d0/135.d0 *zeta*zeta )) &
                   &       /(2.d0**thrd4 -2.d0)
           else
              fz = (onpzeta**thrd4 + onmzeta**thrd4 -2.d0) / (2.d0**thrd4 - 2.d0)
           end if
           g3  = g**3
           g4  = g3*g
           facpon = -delt/(g3*bet)
           if(d < density_minimum) cycle
           dd = grad_trho(i)
           rs = (0.75d0/(PAI*d))**thrd
           fk = (3.d0*PAI*PAI*d)**thrd
           sk = dsqrt(4.d0*fk/PAI)

           if(d > density_minimum2) then
           ! ****** aas mkatsu *****
!              t = dd/(d*sk*2)
              t = dd/(d*sk*2.d0)/g         !    mkatsu aas
              s = dd/(d*fk*2.d0)
           else
              t = 0.d0
              s = 0.d0
           endif

           q0 = -2.d0*a*(1.d0+a1*rs)
           rs12 = dsqrt(rs)
           rs32 = rs12**3
           rsp = rs**p
           q1 = 2.d0 *a*( b1*rs12+b2n*rs+b3*rs32+b4*rs*rsp )
           q2 = log( 1.d0 +1.d0/q1 )
           q3 = a*( b1/rs12 +2.d0*b2n +3.d0*b3*rs12 +2.d0*b4*p1*rsp )

           q0p = -2.d0 *ap*( 1.d0 +a1p*rs )
           q1p = 2.d0 *ap*(b1p*rs12+b2np*rs+b3p*rs32+b4p*rs*rsp)
           q2p = log(1.d0 +1.d0/q1p)
           q3p = ap*(b1p/rs12 +2.d0*b2np +3.d0*b3p*rs12 +2.d0*b4p*p1*rsp)
! vwn 95/12/2 Y.M
           q0q = -2.d0*aq*(1.d0+a1q*rs)
           q1q =  2.d0*aq*(b1q*rs12+b2nq*rs+b3q*rs32+b4q*rs*rsp)
           q2q = log( 1.d0 +1.d0/q1q )
           q3q = aq*(b1q/rs12+2*b2nq+3*b3q*rs12+2*b4q*p1*rsp)

           if(dabs(zeta) < zeta_minimum) then
              fzd = 8.d0/9.d0 /( 2.d0**thrd4 -2.d0 ) &
                   &     *zeta *( 1 + 5.d0/27.d0 *zeta*zeta &
                   &           *( 1 +22.d0/45.d0 *zeta*zeta ))
           else
              fzd = 4.d0/3.d0 /( 2.d0**thrd4 -2.d0 ) &
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

!           zetadxd = - 2 * chgrhr_l(i, 3-is) * (-1.d0)**is / d
           zetadxda = - 2.d0 * chgrhr_l(i, 2) * (-1.d0) / d
           zetadxdb = - 2.d0 * chgrhr_l(i, 1) / d
!           ecd = eu-thrd*rs*eurs + euzt * zetadxd
           ecda = eu-thrd*rs*eurs + euzt * zetadxda 
           ecdb = eu-thrd*rs*eurs + euzt * zetadxdb

#ifndef GGA_PAW_REVISED
           if(onpzeta < zeta_minimum2 .or. onmzeta < zeta_minimum2) then
              ec1  = 0.d0;    ec1d = 0.d0;     ec1dd= 0.d0
              ec1da = 0.d0;     ec1db = 0.d0
              go to 10
           end if
#endif

           pon = facpon*eu
           b   = delt/(dexp(pon)-1.d0)
           b2  = b*b
           t2  = t*t
           s2  = s*s
           t4  = t2*t2
           t6  = t4*t2
           rs2 = rs*rs
           rs3 = rs2*rs
           q4  = 1.d0 +b*t2
           q5  = 1.d0 +b*t2+b2*t4
           q6  = c1+c2*rs+c3*rs2
           q7  = 1.d0+c4*rs+c5*rs2+c6*rs3
           cc  = -cx + q6/q7
           r0  = (sk/fk)**2
           r1  = a4*r0*g4
           coeff = cc-cc0-3.d0*cx/7.d0
           r2  = xnu*coeff*g3
!           a4ms = -a4*g4*s2                                !aas
           a4ms = -r1*t2                                    !aas
           r3  = dexp(a4ms)                       
           h0  = g3*(bet/delt)*dlog(1.d0+delt*q4*t2/q5)
           h1  = r3*r2*t2
           h   = h0+h1
           
           q8  = q5*q5+delt*q4*q5*t2
           h0t = 2.d0 *bet*t*(1+2*b*t2)/q8 * g3
           h0b = -bet*t6*(2.d0*b+b2*t2)/q8 * g3
!           h0rs = h0b * b2 * eurs / bet / g3               !aas
           h0rs = h0b * b * eurs * (b+delt) / bet / g3      !aas
           h1t  = 2.d0*r3*r2*t*(1.d0-r1*t2)
           ccrs = (c2+2*c3*rs)/q7-q6*(c4+2*c5*rs+3*c6*rs2)/q7**2
!           r1rs = 100*r0/rs                                ! aas
           r1rs = r1/rs                                     ! aas
           h1rs = xnu*t2*r3*(ccrs - coeff*t2*r1rs) * g3
           
           ht   = h0t+h1t
           hrs  = h0rs + h1rs

           if(dabs(zeta) < zeta_minimum) then
              gzt = -2.d0/9.d0*zeta *( 1.d0 +14.d0/27.d0 *zeta*zeta &
                   &             *( 1.d0 +13.d0/18.d0 *zeta*zeta ))
           else
#ifdef GGA_PAW_REVISED
              gzt = ( (onpzeta**2+eta)**sixthm - (onmzeta**2+eta)**sixthm)/3.d0
#else
              gzt = ( 1.d0/onpzeta**thrd -1.d0/onmzeta**thrd )/3.d0
#endif
           end if

!           bzt = pon*b**2 * (euzt-3*eu*gzt/g) / ( bet * g3 )           !aas
           bzt = b * ( b+delt) * ( euzt - 3.d0*eu*gzt/g ) / ( bet * g3 )   !aas
           h0zt = 3.d0 * gzt * h0 / g + h0b * bzt
!           h1zt =(3*gzt + 4*gzt*a4ms ) * h1 / g                        !aas
           h1zt = (3.d0 + 4.d0*a4ms) * h1 * gzt / g                    !aas

           ec1  = d*h * 0.5d0
! --> revised by Kato-san,  1999Jul16
!!$           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + h0zt * zetadxd
!           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxd
           ec1da = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*ht/g + h1zt ) * zetadxda
           ec1db = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*ht/g + h1zt ) * zetadxdb
           ec1dd = 0.5d0*ht/sk / g

10         continue
           exc0 = eu*d * 0.5d0 + ec1
!           excd = ecd + ec1d
           excda= ecda + ec1da
           excdb= ecdb + ec1db
!           dF_drho(i, is) = dF_drho(i, is) + excd
           dF_drho(i, 1) = dF_drho(i, 1) + excda
           dF_drho(i, 2) = dF_drho(i, 2) + excdb
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
!           if(is ==  2) then
              if(dabs(grad_trho(i)) > 1.0d-9) then
!                 dF_dgradrho(i) = excdd / grad_trho(i)
                 dF_dgradrho(i) = excdd
              else
                 dF_dgradrho(i) = 0.d0
              endif
!           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d * 0.5
!           exc1 = exc1 + exc0*wos(i)
           exc(i)=exc(i)+exc0*2.d0
!print '(i3,5e19.12)', i,exc0,dF_drho(i, 1:2),dF_dgradrho(i),excdd/ grad_trho(i)
      
! ***** Calc heigher order derivatives. *****

           rsdn=-thrd*rs/d
           rsdnn=-4.d0*thrd*rsdn/d
           rsdnnn=-thrd7*rsdnn/d
           
           id2=1.d0/d/d
           id3=id2/d
           id4=id3/d
           zetada   =     zetadxda/d
           zetadb   =     zetadxdb/d
           zetadaa  =    -4.d0*chgrhr_l(i, 2)*id3
           zetadab  =     2.d0*zeta*id2
           zetadbb  =     4.d0*chgrhr_l(i, 1)*id3
           zetadaaa =    12.d0*chgrhr_l(i, 2)*id4
           zetadaab =   ( -4.d0*chgrhr_l(i, 1) + 8.d0*chgrhr_l(i, 2))*id4
           zetadabb =   ( -8.d0*chgrhr_l(i, 1) + 4.d0*chgrhr_l(i, 2))*id4
           zetadbbb =   -12.d0*chgrhr_l(i, 1)*id4
           
!print '(i2,7d19.10)', i,chgrhr_l(i,1:2),zetadab,zetadaaa,zetadaab,zetadabb,zetadbbb
             
           q1drr =0.5d0*a*( -b1/rs32+3.d0*b3/rs12+8.d0*b4 )
           q1drrr=0.75d0*a*( b1/rs - b3 )/rs32
           invq12pq1=1.d0/(q1**2 + q1)
           q12p1=2.d0*q1+1.d0
           e1dr=-2.d0*a*a1*q2-q0*q3*invq12pq1
           e1drr=4.d0*a*a1*q3*invq12pq1-q0*( q1drr-q12p1*q3**2*invq12pq1 )*invq12pq1
           e1drrr=6.d0*a*a1*( q1drr - q12p1*q3**2*invq12pq1 ) * invq12pq1 &
                    -q0*( q1drrr -(2.d0*q3**3+3.d0*q12p1*q3*q1drr -  &
                    2.d0*q12p1**2*q3**3* invq12pq1)*invq12pq1 ) * invq12pq1
                    
           q1pdrr =0.5d0*ap*( -b1p/rs32+3.d0*b3p/rs12+8.d0*b4p )
           q1pdrrr=0.75d0*ap*( b1p/rs - b3p )/rs32
           invq1p2pq1p=1.d0/(q1p**2 + q1p)
           q1p2p1=2.d0*q1p+1.d0
           e1pdr=-2.d0*ap*a1p*q2p-q0p*q3p*invq1p2pq1p
           e1pdrr=4.d0*ap*a1p*q3p*invq1p2pq1p-q0p*( q1pdrr-q1p2p1*q3p**2*invq1p2pq1p )*invq1p2pq1p
           e1pdrrr=6.d0*ap*a1p*( q1pdrr - q1p2p1*q3p**2*invq1p2pq1p ) * invq1p2pq1p &
                    -q0p*( q1pdrrr -(2.d0*q3p**3+3.d0*q1p2p1*q3p*q1pdrr -  &
                    2.d0*q1p2p1**2*q3p**3* invq1p2pq1p)*invq1p2pq1p ) * invq1p2pq1p
                    
           q1qdrr =0.5d0*aq*( -b1q/rs32+3.d0*b3q/rs12+8.d0*b4q )
           q1qdrrr=0.75d0*aq*( b1q/rs - b3q )/rs32
           invq1q2pq1q=1.d0/(q1q**2 + q1q)
           q1q2p1=2.d0*q1q+1.d0
           e1qdr=-2.d0*aq*a1q*q2q-q0q*q3q*invq1q2pq1q
           e1qdrr=4.d0*aq*a1q*q3q*invq1q2pq1q-q0q*( q1qdrr-q1q2p1*q3q**2*invq1q2pq1q )*invq1q2pq1q
           e1qdrrr=6.d0*aq*a1q*( q1qdrr - q1q2p1*q3q**2*invq1q2pq1q ) * invq1q2pq1q &
                    -q0q*( q1qdrrr -(2.d0*q3q**3+3.d0*q1q2p1*q3q*q1qdrr -  &
                    2.d0*q1q2p1**2*q3q**3* invq1q2pq1q)*invq1q2pq1q ) * invq1q2pq1q
           
           zeta2=zeta**2
           zeta3=zeta2*zeta
           zeta4=zeta2*zeta2
           
           if(dabs(zeta) < zeta_minimum) then
              fzdd0 = fzdd &
                   &            *( 1 + 5.d0*ninth *zeta2 &
                   &            *( 1 +22.d0*twntysvnth *zeta2 ))
              fzddd = fzdd*10.d0*ninth*zeta*( 1.d0 + 44.d0*twntysvnth*zeta2 )
              gzz   = -2.d0*ninth*( 1.d0 + 14.d0*ninth*zeta2*( 1.d0 + 65.d0*fftyfrth*zeta2 ) )
              gzzz  = -56.d0*eityonth*zeta*( 1.d0 + 65.d0*twntysvnth*zeta2 )
           else
              opzthrdm2=1.d0/onpzeta**(thrd2)
              omzthrdm2=1.d0/onmzeta**(thrd2)
              fzdd0 =  fzdd*0.5d0*( opzthrdm2 + omzthrdm2 )
              fzddd = -fzdd*thrd*( opzthrdm2/onpzeta - omzthrdm2/onmzeta )
              gzz   = -ninth*( opzthrdm2**2 + omzthrdm2**2 )
!              gzz   = -ninth*( 1.d0/onpzeta**(thrd4) + 1.d0/onmzeta**(thrd4) )
              gzzz  = 4.d0*twntysvnth*( opzthrdm2**2/onpzeta - omzthrdm2**2/onmzeta )
           end if
           
           e1=q0*q2
           e2=q0p*q2p
           e3=q0q*q2q
           
           eudrr    = e1drr  + ( e1pdrr  -  e1drr  ) * fz * zeta4 - e1qdrr  * fz * onzt4 / fzdd
           eudrrr   = e1drrr + ( e1pdrrr -  e1drrr ) * fz * zeta4 - e1qdrrr * fz * onzt4 / fzdd
           eudzz    = ( e2 - e1 )*( fzdd0*zeta4 + 8.d0*fzd*zeta3 + 12.d0*fz*zeta2 ) &
                        - e3*( fzdd0*onzt4 - 8.d0*fzd*zeta3 - 12.d0*fz*zeta2 )/fzdd
           eudzzz   = ( e2 - e1 )*( fzddd*zeta4 + 12.d0*fzdd0*zeta3 + 36.d0*fzd*zeta2 + 24.d0*fz*zeta ) &
                        - e3*( fzddd*onzt4 - 12.d0*fzdd0*zeta3 - 36.d0*fzd*zeta2 - 24.d0*fz*zeta )/fzdd
           eudrz    = ( e1pdr - e1dr )*( fzd*zeta4 + 4.d0*fz*zeta3 ) &
                        - e1qdr *( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd 
           eudrrz   = ( e1pdrr  - e1drr  )*( fzd*zeta4 + 4.d0*fz*zeta3 ) &
                        - e1qdrr *( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd
           eudrzz   = ( e1pdr - e1dr )*( fzdd0*zeta4 + 8.d0*fzd*zeta3 + 12.d0*fz*zeta2 ) &
                        - e1qdr *( fzdd0*onzt4 - 8.d0*fzd*zeta3 - 12.d0*fz*zeta2 )/fzdd
                        
           rsdn2    = rsdn*rsdn
           rsdn3    = rsdn2*rsdn
           zetada2  = zetada**2
           zetada3  = zetada2*zetada
           zetadb2  = zetadb**2
           zetadb3  = zetadb2*zetadb
                        
           euda     = eurs*rsdn + euzt*zetada
           eudb     = eurs*rsdn + euzt*zetadb
           eudaa    = eudrr*rsdn2 + 2.d0*eudrz*rsdn*zetada + eudzz*zetada2 + eurs*rsdnn + euzt*zetadaa
           eudab    = eudrr*rsdn2 + eudrz*rsdn*( zetada + zetadb ) + eudzz*zetada*zetadb &
                        + eurs*rsdnn + euzt*zetadab
           eudbb    = eudrr*rsdn2 + 2.d0*eudrz*rsdn*zetadb + eudzz*zetadb2 + eurs*rsdnn + euzt*zetadbb
           eudaaa   = eudrrr*rsdn3 &
                        + 3.d0*eudrrz*rsdn2*zetada &
                        + 3.d0*eudrzz*rsdn*zetada2 + eudzzz*zetada3 &
                        + 3.d0*eudrr*rsdn*rsdnn + 3.d0*eudzz*zetada*zetadaa &
                        + 3.d0*eudrz*( rsdnn*zetada + rsdn*zetadaa ) &
                        + eurs*rsdnnn + euzt*zetadaaa
           eudaab   = eudrrr*rsdn3 &
                        + eudrrz*rsdn2*( zetadb + 2.d0*zetada ) &
                        + eudrzz*rsdn*zetada*( zetada + 2.d0*zetadb ) &
                        + eudzzz*zetada2*zetadb &
                        + 3.d0*eudrr*rsdn*rsdnn &
                        + eudrz*( rsdnn*( 2.d0*zetada + zetadb ) + rsdn*( 2.d0*zetadab + zetadaa ) ) &
                        + eudzz*( 2.d0*zetada*zetadab + zetadaa*zetadb ) &
                        + eurs*rsdnnn + euzt*zetadaab
           eudabb   = eudrrr*rsdn3 &
                        + eudrrz*rsdn2*( zetada + 2.d0*zetadb ) &
                        + eudrzz*rsdn*zetadb*( zetadb + 2.d0*zetada ) &
                        + eudzzz*zetada*zetadb2 &
                        + 3.d0*eudrr*rsdn*rsdnn &
                        + eudrz*( rsdnn*( zetada + 2.d0*zetadb ) + rsdn*( 2.d0*zetadab + zetadbb ) ) &
                        + eudzz*( zetada*zetadbb + 2.d0*zetadb*zetadab ) &
                        + eurs*rsdnnn + euzt*zetadabb
           eudbbb   = eudrrr*rsdn3 &
                        + 3.d0*eudrrz*rsdn2*zetadb &
                        + 3.d0*eudrzz*rsdn*zetadb2 + eudzzz*zetadb3 &
                        + 3.d0*eudrr*rsdn*rsdnn + 3.d0*eudzz*zetadb*zetadbb &
                        + 3.d0*eudrz*( rsdnn*zetadb + rsdn*zetadbb ) &
                        + eurs*rsdnnn + euzt*zetadbbb
                        
!print '(i2,12d19.10)', i,chgrhr_l(i,1:2),eu,euda,eudb,eudaa,eudbb,eudab,eudaaa,eudaab,eudabb,eudbbb
!print *,eu,     e1+(e2-e1)* fz * zeta4 - e3  * fz * onzt4 / fzdd
!print *,eurs,   e1dr  + ( e1pdr  -  e1dr  ) * fz * zeta4 - e1qdr  * fz * onzt4 / fzdd
!print *,euzt,   (e2-e1)*( fzd*zeta4 + 4.d0*fz*zeta3 )-e3*( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd 

           duda     = gzt*zetada
           dudb     = gzt*zetadb 
           dudaa    = gzz*zetada2     + gzt*zetadaa
           dudab    = gzz*zetada*zetadb + gzt*zetadab
           dudbb    = gzz*zetadb2     + gzt*zetadbb
           dudaaa   = gzzz*zetada3    + 3.d0*gzz*zetada*zetadaa   + gzt*zetadaaa
           dudaab   = gzzz*zetada2*zetadb + gzz*( zetadaa*zetadb + 2.d0*zetada*zetadab ) + gzt*zetadaab
           dudabb   = gzzz*zetada*zetadb2 + gzz*( zetada*zetadbb + 2.d0*zetadab*zetadb ) + gzt*zetadabb
           dudbbb   = gzzz*zetadb3    + 3.d0*gzz*zetadb*zetadbb   + gzt*zetadbbb 
           
!print '(i2,12d19.10)', i,chgrhr_l(i,1:2),g,duda,dudb,dudaa,dudbb,dudab,dudaaa,dudaab,dudabb,dudbbb

           
           dtdn  = -sixth7*t/d
           dtdnn = -sixth13*dtdn/d
           dtdnnn= -sixth19*dtdnn/d
!print '(i2,5e19.12)',i,d,t,dtdn,dtdnn,dtdnnn
           if(dabs(grad_trho(i)) > 1.0d-9) then
               dtdg  = t/dd
               dtdng = -sixth7*dtdg/d
               dtdnng= -sixth13*dtdng/d
!               dtdnug = -dtdn*invg/dd
           else
               dtdg  = 0.d0
               dtdng = 0.d0
               dtdnng= 0.d0
               dtdnug= 0.d0 
           endif
           dtdgg  = 0.d0
           dtdngg = 0.d0
           dtdggg = 0.d0
           invg=1.d0/g
           dtdu   = -t*invg
           dtduu  = -2.d0*dtdu*invg 
           dtduuu = -3.d0*dtduu*invg
           dtdnu  = -dtdn*invg
           dtdnnu = -dtdnn*invg
           dtdnuu = -sixth7*dtduu/d
           dtdug  = -dtdg*invg
           dtduug = 2.d0*dtdg*invg**2
           dtdugg = 0.d0
           if(dabs(grad_trho(i)) > 1.0d-9) then
                dtdnug = dtdnu/dd
           else
                dtdnug = 0.d0
           end if
           
           dtda   = dtdn + dtdu*duda
           dtdb   = dtdn + dtdu*dudb
           dtdaa  = dtdnn + 2.d0*dtdnu*duda + dtduu*duda**2 + dtdu*dudaa
           dtdbb  = dtdnn + 2.d0*dtdnu*dudb + dtduu*dudb**2 + dtdu*dudbb
           dtdab  = dtdnn + dtdnu*( duda + dudb ) + dtduu*duda*dudb + dtdu*dudab
           dtdag  = dtdng + dtdug*duda
           dtdbg  = dtdng + dtdug*dudb
           dtdaaa = dtdnnn+ 3.d0*dtdnnu*duda + 3.d0*dtdnuu*duda**2 + dtduuu*duda**3 &
                    + 3.d0*dtduu*duda*dudaa + 3.d0*dtdnu*dudaa + dtdu*dudaaa
           dtdbbb = dtdnnn+ 3.d0*dtdnnu*dudb + 3.d0*dtdnuu*dudb**2 + dtduuu*dudb**3 &
                    + 3.d0*dtduu*dudb*dudbb + 3.d0*dtdnu*dudbb + dtdu*dudbbb
           dtdaab = dtdnnn+ dtdnnu*( 2.d0*duda + dudb ) + dtdnuu*duda*( duda + 2.d0*dudb ) + dtduuu*duda**2*dudb &
                    + dtduu*( dudaa*dudb + 2.d0*duda*dudab ) + dtdnu*( dudaa + 2.d0*dudab ) + dtdu*dudaab
           dtdabb = dtdnnn+ dtdnnu*( duda + 2.d0*dudb ) + dtdnuu*dudb*( dudb + 2.d0*duda ) + dtduuu*dudb**2*duda &
                    + dtduu*( dudbb*duda + 2.d0*dudb*dudab ) + dtdnu*( dudbb + 2.d0*dudab ) + dtdu*dudabb 
           dtdaag = dtdnng+ 2.d0*dtdnug*duda + dtduug*duda**2 + dtdug*dudaa
           dtdbbg = dtdnng+ 2.d0*dtdnug*dudb + dtduug*dudb**2 + dtdug*dudbb
           dtdabg = dtdnng+ dtdnug*( duda + dudb ) + dtduug*duda*dudb + dtdug*dudab 
           
!print '(i2,20d19.10)', i,t,dtda,dtdb,dtdg,dtdaa,dtdbb,dtdgg,dtdab,dtdag,dtdbg,dtdaaa,dtdbbb,dtdggg,dtdaab,dtdabb,dtdaag,dtdagg,dtdbbg,dtdbgg,dtdabg     
!print '(i2,12d19.10)', i,t,dtdag,dtdbg,dtdaag,dtdbbg,dtdabg     

           invg3  = 1.d0  /g3
           dAde   = b * ( b + delt ) / bet * invg3
           dAdee  = 2.d0 * dAde * ( b + 0.5d0 * delt ) / bet *invg3
           dAdeee = 2.d0 * ( dAde**2 + ( b + 0.5d0 * delt )*dAdee ) /bet * invg3
           m3u    = -3.d0*invg
           dAdu     = m3u*eu*dAde
           dAdeu    = m3u*( dAde + eu*dAdee )
           dAduu    = m3u*eu*invg*( -dAde + g*dAdeu )
           dAdeeu   = m3u*( 2.d0*dAdee + eu*dAdeee )
           dAdeuu   = m3u*invg*( -dAde + g*dAdeu - eu*dAdee + g*eu*dAdeeu )
           dAduuu   = m3u*eu*invg**2*( 2.d0*dAde -2.d0*g*dAdeu + g**2*dAdeuu )
           
           dAda     = dAde*euda + dAdu*duda
           dAdb     = dAde*eudb + dAdu*dudb
           dAdaa    = dAdee*euda**2 + 2.d0*dAdeu*euda*duda + dAduu*duda**2 + dAde*eudaa + dAdu*dudaa 
           dAdab    = dAdee*euda*eudb + dAdeu*( euda*dudb + eudb*duda ) + dAduu*duda*dudb + dAde*eudab + dAdu*dudab 
           dAdbb    = dAdee*eudb**2 + 2.d0*dAdeu*eudb*dudb + dAduu*dudb**2 + dAde*eudbb + dAdu*dudbb
           dAdaaa   = dAdeee*euda**3 + 3.d0*dAdeeu*euda**2*duda + 3.d0*dAdeuu*euda*duda**2 + dAduuu*duda**3 &
                        + 3.d0*dAdee*euda*eudaa + 3.d0*dAduu*duda*dudaa &
                        + 3.d0*dAdeu*( euda*dudaa + eudaa*duda ) + dAde*eudaaa + dAdu*dudaaa
           dAdaab   = dAdeee*euda**2*eudb + dAduuu*duda**2*dudb &
                        + dAdeeu*euda*( euda*dudb + 2.d0*eudb*duda ) &
                        + dAdeuu*( 2.d0*euda*dudb + eudb*duda )*duda &
                        + dAdee*( 2.d0*euda*eudab + eudb*eudaa ) + dAduu*( 2.d0*duda*dudab + dudb*dudaa ) &
                        + dAdeu*( 2.d0*eudab*duda + 2.d0*euda*dudab + dudb*eudaa + eudb*dudaa ) &
                        + dAde*eudaab + dAdu*dudaab
           dAdabb   = dAdeee*euda*eudb**2 + dAduuu*dudb**2*duda &
                        + dAdeeu*eudb*( eudb*duda + 2.d0*euda*dudb ) &
                        + dAdeuu*( 2.d0*eudb*duda + euda*dudb )*dudb &
                        + dAdee*( eudbb*euda + 2.d0*eudb*eudab ) + dAduu*( dudbb*duda + 2.d0*dudb*dudab ) &
                        + dAdeu*( eudbb*duda + 2.d0*eudb*dudab + 2.d0*dudb*eudab + dudbb*euda ) &
                        + dAde*eudabb + dAdu*dudabb
           dAdbbb   = dAdeee*eudb**3 + 3.d0*dAdeeu*eudb**2*dudb + 3.d0*dAdeuu*eudb*dudb**2 + dAduuu*dudb**3 &
                        + 3.d0*dAdee*eudb*eudbb + 3.d0*dAduu*dudb*dudbb &
                        + 3.d0*dAdeu*( eudb*dudbb + eudbb*dudb ) + dAde*eudbbb + dAdu*dudbbb 
                        
!print '(i2,12d19.10)', i,b,dAda,dAdb,dAdaa,dAdab,dAdbb,dAdaaa,dAdaab,dAdabb,dAdbbb

           at2  = b * t2
           a2t4 = at2**2
           a3t6 = at2*a2t4
           t3   = t2*t
           t5   = t4*t
           t7   = t6*t
           t8   = t7*t
           fE    = q5**2 + delt*t2*q4*q5
           fEdt  = 4.d0*b*t*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt *2.d0*t*(1.d0 + 4.d0*at2 + 6.d0*a2t4 + 4.d0*a3t6 )
           fEdtt = 4.d0*b*( 1.d0 + 9.d0*at2 + 15.d0*a2t4 + 14.d0*a3t6 ) + &
                    delt*( 2.d0 + 24.d0*at2 + 60.d0*a2t4 + 56.d0*a3t6 )
           fEda  = 2.d0*t2*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt*t4*( 2.d0 + 4.d0*at2 + 3.d0*a2t4 )
           fEdaa = 6.d0*t4*( 2.d0*q5 - 1.d0) + &
                    delt*2.d0*t6*( 2.d0 + 3.d0*at2 )
           fF    = t*( 1.d0 + 2.d0*at2 )         
           fFdt  = 1.d0 + 6.d0*at2
           fFdtt = 12.d0*b*t
           fFda  = 2.d0*t3
           fFdaa = 0.d0
           fG    = 2.d0*b*t6 + b2*t8
           fGdt  = 12.d0*b*t5 + 8.d0*b2*t7
           fGdtt = 60.d0*b*t4 + 56.d0*b2*t6
           fGda  = 2.d0*t6 + 2.d0*b*t8
           fGdaa = 2.d0*t8
           
           betivE=bet/fE*g3
           betivE2=betivE/fE
           dHdf   =     -fG*betivE
           dHdt   = 2.d0*fF*betivE
           dHdff  =     -( fGda*fE-fG*fEda )*betivE2
           dHdft  =     -( fGdt*fE-fG*fEdt )*betivE2
           dHdtt  = 2.d0*( fFdt*fE-fF*fEdt )*betivE2
           dHdfff =     -( fGdaa*fE-fG*fEdaa-2.d0*( fGda*fE-fG*fEda )*fEda/fE )*betivE2
           dHdfft = 2.d0*( fFdaa*fE-fF*fEdaa-2.d0*( fFda*fE-fF*fEda )*fEda/fE )*betivE2
           dHdftt =     -( fGdtt*fE-fG*fEdtt-2.d0*( fGdt*fE-fG*fEdt )*fEdt/fE )*betivE2
           dHdttt = 2.d0*( fFdtt*fE-fF*fEdtt-2.d0*( fFdt*fE-fF*fEdt )*fEdt/fE )*betivE2
           
           tiu    = 3.d0*invg
           siu2   = tiu*2.d0*invg
           dHdu   = tiu*h0
           dHduu  = siu2*h0
           dHduuu = siu2*invg*h0
           dHdfu  = tiu*dHdf
           dHdtu  = tiu*dHdt
           dHdffu = tiu*dHdff
           dHdfuu = siu2*dHdf
           dHdttu = tiu*dHdtt
           dHdtuu = siu2*dHdt
           dHdftu = tiu*dHdft
           
           dHda   = dHdf*dAda + dHdt*dtda + dHdu*duda
           dHdb   = dHdf*dAdb + dHdt*dtdb + dHdu*dudb
           dHdg   =             dHdt*dtdg
           dHdaa  = dHdff*dAda**2 + dHdtt*dtda**2 + dHduu*duda**2 &
                    + 2.d0*dHdft*dAda*dtda + 2.d0*dHdfu*dAda*duda + 2.d0*dHdtu*dtda*duda &
                    + dHdf*dAdaa + dHdt*dtdaa + dHdu*dudaa
           dHdab  = dHdff*dAda*dAdb + dHdtt*dtda*dtdb + dHduu*duda*dudb &
                    + dHdft*( dAda*dtdb + dAdb*dtda ) + dHdfu*( dAda*dudb + dAdb*duda ) + dHdtu*( dtda*dudb + dtdb*duda ) &
                    + dHdf*dAdab + dHdt*dtdab + dHdu*dudab
           dHdbb  = dHdff*dAdb**2 + dHdtt*dtdb**2 + dHduu*dudb**2 &
                    + 2.d0*dHdft*dAdb*dtdb + 2.d0*dHdfu*dAdb*dudb + 2.d0*dHdtu*dtdb*dudb &
                    + dHdf*dAdbb + dHdt*dtdbb + dHdu*dudbb
           dHdgg  = dHdtt*dtdg**2
           dHdag  = dHdft*dAda*dtdg + dHdtt*dtda*dtdg + dHdtu*duda*dtdg + dHdt*dtdag
           dHdbg  = dHdft*dAdb*dtdg + dHdtt*dtdb*dtdg + dHdtu*dudb*dtdg + dHdt*dtdbg
           dHdaaa = dHdfff*dAda**3 + 3.d0*dHdfft*dAda**2*dtda + 3.d0*dHdffu*dAda**2*duda &
                    + 3.d0*dHdftt*dAda*dtda**2 + dHdttt*dtda**3 + 3.d0*dHdttu*dtda**2*duda &
                    + 3.d0*dHdfuu*dAda*duda**2 + 3.d0*dHdtuu*dtda*duda**2 + dHduuu*duda**3 &
                    + 6.d0*dHdftu*dAda*dtda*duda &
                    + 3.d0*dHdff*dAda*dAdaa + 3.d0*dHdtt*dtda*dtdaa + 3.d0*dHduu*duda*dudaa &
                    + 3.d0*dHdft*( dAdaa*dtda + dAda*dtdaa ) &
                    + 3.d0*dHdtu*( dtdaa*duda + dtda*dudaa ) &
                    + 3.d0*dHdfu*( dAdaa*duda + dAda*dudaa ) &
                    + dHdf*dAdaaa + dHdt*dtdaaa + dHdu*dudaaa
           dHdaab = dHdfff*dAda**2*dAdb + dHdfft*dAda*( dAda*dtdb + 2.d0*dAdb*dtda ) &
                    + dHdffu*dAda*( dAda*dudb + 2.d0*dAdb*duda ) &
                    + dHdftt*dtda*( dAdb*dtda + 2.d0*dAda*dtdb ) &
                    + dHdttt*dtdb*dtda**2 + dHdttu*dtda*( dtda*dudb + 2.d0*dtdb*duda ) &
                    + dHdfuu*duda*( dAdb*duda + 2.d0*dAda*dudb ) &
                    + dHdtuu*duda*( dtdb*duda + 2.d0*dtda*dudb ) + dHduuu*duda**2*dudb &
                    + 2.d0*dHdftu*( dAdb*duda*dtda + dAda*dudb*dtda + dAda*duda*dtdb ) &
                    + dHdff*( 2.d0*dAda*dAdab + dAdb*dAdaa ) &
                    + dHdtt*( 2.d0*dtda*dtdab + dtdb*dtdaa ) &
                    + dHduu*( 2.d0*duda*dudab + dudb*dudaa ) &
                    + dHdft*( 2.d0*dAdab*dtda + 2.d0*dAda*dtdab + dAdaa*dtdb + dAdb*dtdaa ) &
                    + dHdtu*( 2.d0*dtdab*duda + 2.d0*dtda*dudab + dtdaa*dudb + dtdb*dudaa ) &
                    + dHdfu*( 2.d0*dAdab*duda + 2.d0*dAda*dudab + dAdaa*dudb + dAdb*dudaa ) &
                    + dHdf*dAdaab + dHdt*dtdaab + dHdu*dudaab
           dHdabb = dHdfff*dAda*dAdb**2 + dHdfft*dAdb*( dAdb*dtda + 2.d0*dAda*dtdb ) &
                    + dHdffu*dAdb*( dAdb*duda + 2.d0*dAda*dudb ) &
                    + dHdftt*dtdb*( dAda*dtdb + 2.d0*dAdb*dtda ) + dHdttt*dtda*dtdb**2 &
                    + dHdttu*dtdb*( dtdb*duda + 2.d0*dtda*dudb ) &
                    + dHdfuu*dudb*( dAda*dudb + 2.d0*dAdb*duda ) + dHdtuu*dudb*( dtda*dudb &
                    + 2.d0*dtdb*duda ) + dHduuu*duda*dudb**2 &
                    + 2.d0*dHdftu*( dAda*dtdb*dudb + dAdb*dtda*dudb + dAdb*dtdb*duda ) &
                    + dHdff*( 2.d0*dAdb*dAdab + dAda*dAdbb ) &
                    + dHdtt*( 2.d0*dtdb*dtdab + dtda*dtdbb ) &
                    + dHduu*( 2.d0*dudb*dudab + duda*dudbb ) &
                    + dHdft*( 2.d0*dAdab*dtdb + 2.d0*dAdb*dtdab + dAdbb*dtda + dAda*dtdbb ) &
                    + dHdtu*( 2.d0*dtdab*dudb + 2.d0*dtdb*dudab + dtdbb*duda + dtda*dudbb ) &
                    + dHdfu*( 2.d0*dAdab*dudb + 2.d0*dAdb*dudab + dAdbb*duda + dAda*dudbb ) &
                    + dHdf*dAdabb + dHdt*dtdabb + dHdu*dudabb
           dHdbbb = dHdfff*dAdb**3 + 3.d0*dHdfft*dAdb**2*dtdb + 3.d0*dHdffu*dAdb**2*dudb &
                    + 3.d0*dHdftt*dAdb*dtdb**2 + dHdttt*dtdb**3 + 3.d0*dHdttu*dtdb**2*dudb &
                    + 3.d0*dHdfuu*dAdb*dudb**2 + 3.d0*dHdtuu*dtdb*dudb**2 + dHduuu*dudb**3 &
                    + 6.d0*dHdftu*dAdb*dtdb*dudb &
                    + 3.d0*dHdff*dAdb*dAdbb + 3.d0*dHdtt*dtdb*dtdbb + 3.d0*dHduu*dudb*dudbb &
                    + 3.d0*dHdft*( dAdbb*dtdb + dAdb*dtdbb ) &
                    + 3.d0*dHdtu*( dtdbb*dudb + dtdb*dudbb ) &
                    + 3.d0*dHdfu*( dAdbb*dudb + dAdb*dudbb ) &
                    + dHdf*dAdbbb + dHdt*dtdbbb + dHdu*dudbbb
           dHdabg = dHdfft*dAda*dAdb*dtdg   + dHdftt*( dAdb*dtda + dAda*dtdb )*dtdg &
                    + dHdttt*dtda*dtdb*dtdg + dHdttu*( dtdb*duda + dtda*dudb )*dtdg &
                    + dHdtuu*duda*dudb*dtdg + dHdftu*( dAdb*duda + dAda*dudb )*dtdg &
                    + dHdft*( dAdab*dtdg + dAdb*dtdag + dAda*dtdbg ) &
                    + dHdtt*( dtdab*dtdg + dtdb*dtdag + dtda*dtdbg ) &
                    + dHdtu*( dudab*dtdg + dudb*dtdag + duda*dtdbg ) &
                    + dHdt*dtdabg
           dHdaag = dHdfft*dAda**2*dtdg + dHdttt*dtda**2*dtdg + dHdtuu*duda**2*dtdg &
                    + 2.d0*dHdttu*dtda*duda*dtdg + 2.d0*dHdftt*dAda*dtda*dtdg + 2.d0*dHdftu*dAda*duda*dtdg &
                    + dHdtt*( dtdaa*dtdg + 2.d0*dtda*dtdag ) + dHdft*( dAdaa*dtdg + 2.d0*dAda*dtdag ) &
                    + dHdtu*( dudaa*dtdg + 2.d0*duda*dtdag ) + dHdt*dtdaag         
           dHdbbg = dHdfft*dAdb**2*dtdg + dHdttt*dtdb**2*dtdg + dHdtuu*dudb**2*dtdg &
                    + 2.d0*dHdttu*dtdb*dudb*dtdg + 2.d0*dHdftt*dAdb*dtdb*dtdg + 2.d0*dHdftu*dAdb*dudb*dtdg &
                    + dHdtt*( dtdbb*dtdg + 2.d0*dtdb*dtdbg ) + dHdft*( dAdbb*dtdg + 2.d0*dAdb*dtdbg ) &
                    + dHdtu*( dudbb*dtdg + 2.d0*dudb*dtdbg ) + dHdt*dtdbbg       
           dHdagg = dHdftt*dAda*dtdg**2 + dHdttt*dtda*dtdg**2 + dHdttu*duda*dtdg**2 + 2.d0*dHdtt*dtdg*dtdag
           dHdbgg = dHdftt*dAdb*dtdg**2 + dHdttt*dtdb*dtdg**2 + dHdttu*dudb*dtdg**2 + 2.d0*dHdtt*dtdg*dtdbg 
           dHdggg = dHdttt*dtdg**3   
           
! ****** Derivative of h1 *****
           
           dq6   = c2+2.d0*c3*rs
           ddq6  = 2.d0*c3
           dq7   = c4+2.d0*c5*rs+3.d0*c6*rs2
           ddq7  = 2.d0*c5+6.d0*c6*rs
           dddq7 = 6.d0*c6
           
           iq7=1.d0/q7
           iq72=iq7*iq7
           iq73=iq72*iq7
           iq74=iq73*iq7
           
           ccrs2 =  ddq6*iq7-(2.d0*dq6*dq7+q6*ddq7)*iq72+2.d0*q6*dq7**2*iq73
           ccrs3 =  -iq72*(3.d0*(ddq6*dq7+dq6*ddq7)+q6*dddq7) + &
                                    6.d0*dq7*(dq6*dq7+q6*ddq7)*iq73 - &
                                    6.d0*q6*dq7**3*iq74
           r1rs2=r1rs*r1rs
           r1rs3=r1rs2*r1rs
           r1rs4=r1rs3*r1rs
           
           if(dabs(grad_trho(i)) > 1.0d-9) then
               ivt=1.d0/t
               ivt2=ivt**2
           else
               ivt=0.d0
               ivt2=0.d0
           end if
           
!           ivtmr1t=ivt-r1*t
           ivu=1.d0/g
           ivu2=1.d0/g/g
           r1t2=r1*t2
           r1t22= r1t2**2
           r1t23= r1t2*r1t22
           
           dH1dr   =  h1rs
           dH1du   =  h1*ivu*(3.d0-4.d0*r1t2)
           dH1dt   =  h1t
           dH1drr  =  xnu*ccrs2*t2*g3*r3-2.d0*h1rs*r1rs*t2-h1*r1rs2*t4
           dH1duu  =  2.d0*h1*ivu2*(3.d0 - 18.d0*r1t2 + 8.d0*r1t22) 
           dH1dtt  =  2.d0*h1*ivt2*(1.d0 -  5.d0*r1t2 + 2.d0*r1t22)
           dH1dru  =  ivu*(dH1dr*(3.d0 - 4.d0*r1t2)-4.d0*h1*r1rs*t2)
           dH1dtu  =  2.d0*ivu*ivt*h1*(3.d0-11.d0*r1t2+4.d0*r1t22)
           dH1drt  =  2.d0*(ivt*dH1dr*(1.d0-r1t2)-h1*r1rs*t)
           dH1drrr =  xnu*ccrs3*t2*g3*r3-3.d0*t2*r1rs*(dH1drr+h1rs*r1rs*t2)-h1*r1rs3*t6
           dH1drru =  ivu*(dH1drr*(3.d0-4.d0*r1t2)-8.d0*t2*dH1dr*r1rs)
           dH1druu =  2.d0*ivu2*(dH1dr*(3.d0-18.d0*r1t2+8.d0*r1t22)+2.d0*h1*t2*r1rs*(8.d0*r1t2-9.d0))
           dH1duuu =  2.d0*h1*ivu*ivu2*(3.d0-102.d0*r1t2+144.d0*r1t22-32.d0*r1t23)
           dH1dtuu =  4.d0*h1*ivu2*ivt*(3.d0-39.d0*r1t2+42.d0*r1t22-8.d0*r1t23)
           dH1dttu =  0.5d0*g*ivt*dH1dtuu
           dH1dttt =  4.d0*h1*r1*ivt*(-6.d0+9.d0*r1t2-2.d0*r1t22)
           dH1drtt =  2.d0*dH1dr*ivt2*(1.d0-5.d0*r1t2+2.d0*r1t22)+2.d0*h1*r1rs*(4.d0*r1t2-5.d0)
           dH1drrt =  2.d0*ivt*dH1drr*(1.d0-r1t2)-4.d0*dH1dr*r1rs*t
           dH1drtu =  2.d0*ivu*( ivt*dH1dr*(3.d0-11.d0*r1t2+4.d0*r1t22)+h1*t*r1rs*(-11.d0+8.d0*r1t2))
           
           dH1da   = dH1dr*rsdn + dH1dt*dtda + dH1du*duda
           dH1db   = dH1dr*rsdn + dH1dt*dtdb + dH1du*dudb
           dH1dg   =              dH1dt*dtdg
           dH1daa  = dH1drr*rsdn**2 + dH1dtt*dtda**2 + dH1duu*duda**2 &
                    + 2.d0*dH1drt*rsdn*dtda + 2.d0*dH1dru*rsdn*duda + 2.d0*dH1dtu*dtda*duda &
                    + dH1dr*rsdnn + dH1dt*dtdaa + dH1du*dudaa
           dH1dab  = dH1drr*rsdn*rsdn + dH1dtt*dtda*dtdb + dH1duu*duda*dudb &
                    + dH1drt*rsdn*(dtdb + dtda) + dH1dru*rsdn*( dudb + duda ) + dH1dtu*( dtda*dudb + dtdb*duda ) &
                    + dH1dr*rsdnn + dH1dt*dtdab + dH1du*dudab
           dH1dbb  = dH1drr*rsdn**2 + dH1dtt*dtdb**2 + dH1duu*dudb**2 &
                    + 2.d0*dH1drt*rsdn*dtdb + 2.d0*dH1dru*rsdn*dudb + 2.d0*dH1dtu*dtdb*dudb &
                    + dH1dr*rsdnn + dH1dt*dtdbb + dH1du*dudbb
           dH1dgg  = dH1dtt*dtdg**2
           dH1dag  = dH1drt*rsdn*dtdg + dH1dtt*dtda*dtdg + dH1dtu*duda*dtdg + dH1dt*dtdag
           dH1dbg  = dH1drt*rsdn*dtdg + dH1dtt*dtdb*dtdg + dH1dtu*dudb*dtdg + dH1dt*dtdbg
           dH1daaa = dH1drrr*rsdn3 + 3.d0*dH1drrt*rsdn2*dtda + 3.d0*dH1drru*rsdn2*duda &
                    + 3.d0*dH1drtt*rsdn*dtda**2 + dH1dttt*dtda**3 + 3.d0*dH1dttu*dtda**2*duda &
                    + 3.d0*dH1druu*rsdn*duda**2 + 3.d0*dH1dtuu*dtda*duda**2 + dH1duuu*duda**3 &
                    + 6.d0*dH1drtu*rsdn*dtda*duda &
                    + 3.d0*dH1drr*rsdn*rsdnn + 3.d0*dH1dtt*dtda*dtdaa + 3.d0*dH1duu*duda*dudaa &
                    + 3.d0*dH1drt*( rsdnn*dtda + rsdn*dtdaa ) &
                    + 3.d0*dH1dtu*( dtdaa*duda + dtda*dudaa ) &
                    + 3.d0*dH1dru*( rsdnn*duda + rsdn*dudaa ) &
                    + dH1dr*rsdnnn + dH1dt*dtdaaa + dH1du*dudaaa
           dH1daab = dH1drrr*rsdn3 + dH1drrt*rsdn2*( dtdb + 2.d0*dtda ) + dH1drru*rsdn2*( dudb + 2.d0*duda ) &
                    + dH1drtt*dtda*rsdn*( dtda + 2.d0*dtdb ) + dH1dttt*dtdb*dtda**2 &
                    + dH1dttu*dtda*( dtda*dudb + 2.d0*dtdb*duda ) &
                    + dH1druu*duda*rsdn*( duda + 2.d0*dudb ) &
                    + dH1dtuu*duda*( dtdb*duda + 2.d0*dtda*dudb ) + dH1duuu*duda**2*dudb &
                    + 2.d0*dH1drtu*rsdn*( duda*dtda + dudb*dtda + duda*dtdb ) &
                    + dH1drr*3.d0*rsdn*rsdnn &
                    + dH1dtt*( 2.d0*dtda*dtdab + dtdb*dtdaa ) &
                    + dH1duu*( 2.d0*duda*dudab + dudb*dudaa ) &
                    + dH1drt*( 2.d0*rsdnn*dtda + 2.d0*rsdn*dtdab + rsdnn*dtdb + rsdn*dtdaa ) &
                    + dH1dtu*( 2.d0*dtdab*duda + 2.d0*dtda*dudab + dtdaa*dudb + dtdb*dudaa ) &
                    + dH1dru*( 2.d0*rsdnn*duda + 2.d0*rsdn*dudab + rsdnn*dudb + rsdn*dudaa ) &
                    + dH1dr*rsdnnn + dH1dt*dtdaab + dH1du*dudaab
           dH1dabb = dH1drrr*rsdn3 + dH1drrt*rsdn2*( dtda + 2.d0*dtdb ) + dH1drru*rsdn2*( duda + 2.d0*dudb ) &
                    + dH1drtt*dtdb*rsdn*( dtdb + 2.d0*dtda ) &
                    + dH1dttt*dtda*dtdb**2 + dH1dttu*dtdb*( dtdb*duda + 2.d0*dtda*dudb ) &
                    + dH1druu*dudb*rsdn*( dudb + 2.d0*duda ) &
                    + dH1dtuu*dudb*( dtda*dudb + 2.d0*dtdb*duda ) + dH1duuu*duda*dudb**2 &
                    + 2.d0*dH1drtu*rsdn*( dtdb*dudb + dtda*dudb + dtdb*duda ) &
                    + dH1drr*3.d0*rsdn*rsdnn &
                    + dH1dtt*( 2.d0*dtdb*dtdab + dtda*dtdbb ) &
                    + dH1duu*( 2.d0*dudb*dudab + duda*dudbb ) &
                    + dH1drt*( 2.d0*rsdnn*dtdb + 2.d0*rsdn*dtdab + rsdnn*dtda + rsdn*dtdbb ) &
                    + dH1dtu*( 2.d0*dtdab*dudb + 2.d0*dtdb*dudab + dtdbb*duda + dtda*dudbb ) &
                    + dH1dru*( 2.d0*rsdnn*dudb + 2.d0*rsdn*dudab + rsdnn*duda + rsdn*dudbb ) &
                    + dH1dr*rsdnnn + dH1dt*dtdabb + dH1du*dudabb
           dH1dbbb = dH1drrr*rsdn3 + 3.d0*dH1drrt*rsdn2*dtdb + 3.d0*dH1drru*rsdn2*dudb &
                    + 3.d0*dH1drtt*rsdn*dtdb**2 + dH1dttt*dtdb**3 + 3.d0*dH1dttu*dtdb**2*dudb &
                    + 3.d0*dH1druu*rsdn*dudb**2 + 3.d0*dH1dtuu*dtdb*dudb**2 + dH1duuu*dudb**3 &
                    + 6.d0*dH1drtu*rsdn*dtdb*dudb &
                    + 3.d0*dH1drr*rsdn*rsdnn + 3.d0*dH1dtt*dtdb*dtdbb + 3.d0*dH1duu*dudb*dudbb &
                    + 3.d0*dH1drt*( rsdnn*dtdb + rsdn*dtdbb ) &
                    + 3.d0*dH1dtu*( dtdbb*dudb + dtdb*dudbb ) &
                    + 3.d0*dH1dru*( rsdnn*dudb + rsdn*dudbb ) &
                    + dH1dr*rsdnnn + dH1dt*dtdbbb + dH1du*dudbbb
           dH1dabg = dH1drrt*rsdn2*dtdg   + dH1drtt*rsdn*( dtda + dtdb )*dtdg &
                    + dH1dttt*dtda*dtdb*dtdg + dH1dttu*( dtdb*duda + dtda*dudb )*dtdg &
                    + dH1dtuu*duda*dudb*dtdg + dH1drtu*rsdn*( duda + dudb )*dtdg &
                    + dH1drt*( rsdnn*dtdg + rsdn*dtdag + rsdn*dtdbg ) &
                    + dH1dtt*( dtdab*dtdg + dtdb*dtdag + dtda*dtdbg ) &
                    + dH1dtu*( dudab*dtdg + dudb*dtdag + duda*dtdbg ) &
                    + dH1dt*dtdabg
           dH1daag = dH1drrt*rsdn2*dtdg + dH1dttt*dtda**2*dtdg + dH1dtuu*duda**2*dtdg &
                    + 2.d0*dH1dttu*dtda*duda*dtdg + 2.d0*dH1drtt*rsdn*dtda*dtdg + 2.d0*dH1drtu*rsdn*duda*dtdg &
                    + dH1dtt*( dtdaa*dtdg + 2.d0*dtda*dtdag ) + dH1drt*( rsdnn*dtdg + 2.d0*rsdn*dtdag ) &
                    + dH1dtu*( dudaa*dtdg + 2.d0*duda*dtdag ) + dH1dt*dtdaag  
           dH1dbbg = dH1drrt*rsdn2*dtdg + dH1dttt*dtdb**2*dtdg + dH1dtuu*dudb**2*dtdg &
                    + 2.d0*dH1dttu*dtdb*dudb*dtdg + 2.d0*dH1drtt*rsdn*dtdb*dtdg + 2.d0*dH1drtu*rsdn*dudb*dtdg &
                    + dH1dtt*( dtdbb*dtdg + 2.d0*dtdb*dtdbg ) + dH1drt*( rsdnn*dtdg + 2.d0*rsdn*dtdbg ) &
                    + dH1dtu*( dudbb*dtdg + 2.d0*dudb*dtdbg ) + dH1dt*dtdbbg       
           dH1dagg = dH1drtt*rsdn*dtdg**2 + dH1dttt*dtda*dtdg**2 + dH1dttu*duda*dtdg**2 + 2.d0*dH1dtt*dtdg*dtdag
           dH1dbgg = dH1drtt*rsdn*dtdg**2 + dH1dttt*dtdb*dtdg**2 + dH1dttu*dudb*dtdg**2 + 2.d0*dH1dtt*dtdg*dtdbg 
           dH1dggg = dH1dttt*dtdg**3   
                    
!print '(i2,20d19.10)', i,h0,dHda,dHdb,dHdg,dHdaa,dHdbb,dHdgg,dHdag,dHdbg,dHdab,dHdaaa,dHdbbb,dHdggg,dHdaab,dHdabb,dHdaag,dHdagg,dHdbbg,dHdbgg,dHdabg
!print '(i2,20d19.10)', i,h1,dH1da,dH1db,dH1dg,dH1daa,dH1dbb,dH1dgg,dH1dag,dH1dbg,dH1dab,dH1daaa,dH1dbbb,dH1dggg,dH1daab,dH1dabb,dH1daag,dH1dagg,dH1dbbg,dH1dbgg,dH1dabg
    
           dFc_daa(i)   = 2.d0*euda + d*eudaa + 2.d0*(dHda + dH1da) + d*(dHdaa + dH1daa) 
           dFc_dbb(i)   = 2.d0*eudb + d*eudbb + 2.d0*(dHdb + dH1db) + d*(dHdbb + dH1dbb)
           dFc_dgg(i)   = d*(dHdgg + dH1dgg)
           dFc_dab(i)   = euda + eudb + d*eudab + dHda + dH1da + dHdb + dH1db + d*(dHdab + dH1dab) 
           dFc_dag(i)   = dHdg + dH1dg + d*(dHdag + dH1dag)
           dFc_dbg(i)   = dHdg + dH1dg + d*(dHdbg + dH1dbg)
           dFc_daaa(i)  = 3.d0*eudaa + d*eudaaa + 3.d0*(dHdaa + dH1daa) + d*(dHdaaa + dH1daaa)
           dFc_dbbb(i)  = 3.d0*eudbb + d*eudbbb + 3.d0*(dHdbb + dH1dbb) + d*(dHdbbb + dH1dbbb)
           dFc_dggg(i)  = d*(dHdggg + dH1dggg)
           dFc_daab(i)  = eudaa + 2.d0*eudab + d*eudaab + dHdaa + dH1daa + 2.d0*(dHdab + dH1dab) + d*(dHdaab + dH1daab)
           dFc_daag(i)  = 2.d0*(dHdag + dH1dag) + d*(dHdaag + dH1daag)
           dFc_dabb(i)  = eudbb + 2.d0*eudab + d*eudabb + dHdbb + dH1dbb + 2.d0*(dHdab + dH1dab) + d*(dHdabb + dH1dabb)
           dFc_dbbg(i)  = 2.d0*(dHdbg + dH1dbg) + d*(dHdbbg + dH1dbbg)
           dFc_dagg(i)  = dHdgg + dH1dgg + d*(dHdagg + dH1dagg)
           dFc_dbgg(i)  = dHdgg + dH1dgg + d*(dHdbgg + dH1dbgg)
           dFc_dabg(i)  = dHdag + dH1dag + dHdbg + dH1dbg + d*(dHdabg + dH1dabg)
                    
!print '(i2,21d19.10)', i,exc(i),d*(eu+h),dF_drho(i,1),dF_drho(i,2),dF_dgradrho(i)*grad_trho(i) &
!                      ,dFc_daa(i),dFc_dbb(i),dFc_dgg(i),dFc_dab(i),dFc_dag(i),dFc_dbg(i),dFc_daaa(i),dFc_dbbb(i) &
!                      ,dFc_dggg(i),dFc_daab(i),dFc_daag(i),dFc_dabb(i),dFc_dbbg(i),dFc_dagg(i),dFc_dbgg(i),dFc_dabg(i)
!print *,dHdag,dHdaag,d,2.d0*dHdag + d*dHdaag
!print *,2.d0*dHdbg + d*dHdbbg
!print *,dHdag + dHdbg + d*dHdabg
               
           
        end do
     end if
     
end subroutine cr_ggapw91_paw_drv2
! === For nrc decomposion. by takto 2012/12/05 =================================
subroutine cr_ggapw91_paw_drv2_3D(nrc,dnr,nspin,chgrhr_l,grad_trho,exc,dF_drho,dF_dgradrho &
                                ,dFc_daa,dFc_dbb,dFc_dgg,dFc_dab,dFc_dag,dFc_dbg &
                                ,dFc_daaa,dFc_dbbb,dFc_dggg &
                                ,dFc_daab,dFc_daag,dFc_dabb,dFc_dbbg,dFc_dagg,dFc_dbgg,dFc_dabg &
                                ,ista_nrc,iend_nrc,ist,ien)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : PAI,DP
  implicit none

  integer,intent(in)        :: nrc,dnr,nspin
  real(kind=DP),intent(in)  :: chgrhr_l(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(in) :: grad_trho(ista_nrc:iend_nrc)
  real(kind=DP),intent(inout) :: exc(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dF_drho(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dF_dgradrho(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_daa(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbb(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dgg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dab(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dag(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_daaa(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbbb(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dggg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_daab(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_daag(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dabb(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbbg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dagg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbgg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dabg(ista_nrc:iend_nrc)
  integer :: ista_nrc, iend_nrc, ist, ien

  real(kind=DP), parameter :: zeta_minimum    = 1.d-4
  real(kind=DP), parameter :: zeta_minimum2   = 1.d-6
  real(kind=DP), parameter :: density_minimum = 1.d-20
  real(kind=DP), parameter :: density_minimum2= 1.d-5

  real(kind=DP), parameter  :: a  = 0.0310907d0,    a1   = 0.21370d0
  real(kind=DP), parameter  :: b1 = 7.5957d0,       b2n  = 3.5876d0
  real(kind=DP), parameter  :: b3 = 1.6382d0,       b4   = 0.49294d0
  real(kind=DP), parameter  :: p  = 1.00d0,         p1   = p+1.d0

#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: ap  = 0.01554535d0,    a1p  = 0.20548d0
#else
  real(kind=DP), parameter  :: ap  = 0.015545d0,    a1p  = 0.20548d0
#endif

  real(kind=DP), parameter  :: b1p = 14.1189d0,     b2np = 6.1977d0
  real(kind=DP), parameter  :: b3p = 3.3662d0,      b4p  = 0.62517d0
! vwn 95/12/2 Y.M
#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: aq  = 0.0168869d0,    a1q  = 0.11125d0
#else
  real(kind=DP), parameter  :: aq  = 0.016887d0,    a1q  = 0.11125d0
#endif

  real(kind=DP), parameter  :: b1q = 10.3570d0,     b2nq = 3.6231d0
  real(kind=DP), parameter  :: b3q = 0.88026d0,     b4q  = 0.49671d0
! vwn 95/12/2 Y.M
#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: gamma = 0.03109069086965489503494086371273d0
#else
  real(kind=DP), parameter  :: gamma = 0.031091
#endif

!--------------------------
  real(kind=DP), parameter :: cx  = -0.001667212d0, alf =  0.09d0
  real(kind=DP), parameter :: c1  =  0.002568d0,    c2  =  0.023266d0
  real(kind=DP), parameter :: c3  =  7.389d-6,      c4  =  8.723d0
  real(kind=DP), parameter :: c5  =  0.472d0,       c6  =  7.389d-2
  real(kind=DP), parameter :: a4  = 100.d0
  
  real(kind=DP), parameter  :: xnu = 15.75592d0, cc0 = 0.004235d0

#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: thrd = 0.33333333333333330,&
       &                       sixth7 = 1.1666666666666666d0
  real(kind=DP), parameter  :: thrd7= 2.333333333333333333d0
  real(kind=DP), parameter  :: sixth13= 2.166666666666666666d0
  real(kind=DP), parameter  :: sixth19= 3.166666666666666666d0
  real(kind=DP), parameter  :: ninth= 0.11111111111111111111d0

#else
  real(kind=DP), parameter  :: thrd = 0.333333333333d0, sixth7 = 1.1666666666666666d0
  real(kind=DP), parameter  :: thrd7= 2.333333333333d0
  real(kind=DP), parameter  :: sixth13= 2.1666666666666666d0
  real(kind=DP), parameter  :: sixth19= 3.1666666666666666d0
  real(kind=DP), parameter  :: ninth= 0.111111111111111d0
#endif

  real(kind=DP), parameter  :: twntysvnth=  3.703703703703703d-2
  real(kind=DP), parameter  :: fftyfrth= 1.851851851851852d-2
  real(kind=DP), parameter  :: eityonth= 1.234567901234568d-2

  integer       :: is,i
  real(kind=DP) :: facw,bet,delt,exc1,g,g3,g4,facpon,d,dd,rs,fk,sk,t,s,q0 &
       &         , rs12,rs32,rsp,q1,q2,eu,q3,eurs,ecd,pon &
       &         , b,b2,t2,s2,t4,t6,rs2,rs3,q4,q5,q6,q7,cc,r0,r1,coeff,r2,r3,a4ms,h0,h1,h,q8 &
       &         , h0t,h0b,h0rs,h1t,ccrs,r1rs,h1rs,ht,hrs &
       &         , ec1,ec1d,ec1dd,exc0,excd,excdd &
       &         , thrd2,thrd4,fzdd,zeta,onpzeta,onmzeta,fz,zetadxd &
       &         , q0p,q1p,q2p,q3p,q0q,q1q,q2q,q3q,fzd,onzt4,gzt,bzt,h0zt,h1zt,euzt &
       &         , zetadxda,zetadxdb,ecda,ecdb,ec1da,ec1db,excda,excdb
  real(kind=DP) :: ccrs2,ccrs3,dq6,ddq6,dq7,ddq7,dddq7
  real(kind=DP) :: iq7,iq72,iq73,iq74
       
  real(kind=DP) :: rsdn,rsdnn,rsdnnn
  real(kind=DP) :: q1drr,q1drrr
  real(kind=DP) :: eudrr,eudrrr,invq12pq1,q12p1
  real(kind=DP) :: eudn,eudnn,eudnnn,ecdd,ecddd
  real(kind=DP) :: dtdn,dtdg,dtdnn,dtdgg,dtdng,dtdnnn,dtdnng,dtdngg,dtdggg
  real(kind=DP) :: dAde,dAdee,dAdeee
  real(kind=DP) :: dAdn,dAdnn,dAdnnn
  real(kind=DP) :: fE,fEdt,fEdtt,fEda,fEdaa
  real(kind=DP) :: fF,fFdt,fFdtt,fFda,fFdaa
  real(kind=DP) :: fG,fGdt,fGdtt,fGda,fGdaa
  real(kind=DP) :: at2,a2t4,a3t6,t3,t5,t7,t8
  real(kind=DP) :: dHda,dHdt,dHdaa,dHdat,dHdtt,dHdaaa,dHdaat,dHdatt,dHdttt,betivE,betivE2
  real(kind=DP) :: dHdn,dHdg,dHdnn,dHdng,dHdgg,dHdnnn,dHdnng,dHdngg,dHdggg
  real(kind=DP) :: zetada,zetadb,zetadaa,zetadab,zetadbb,zetadaaa,zetadaab,zetadabb,zetadbbb
  real(kind=DP) :: id2,id3,id4
  real(kind=DP) :: e1dr,e1drr,e1drrr
  real(kind=DP) :: e1pdr,e1pdrr,e1pdrrr
  real(kind=DP) :: e1qdr,e1qdrr,e1qdrrr
  real(kind=DP) :: q1pdrr,q1pdrrr,invq1p2pq1p,q1p2p1
  real(kind=DP) :: q1qdrr,q1qdrrr,invq1q2pq1q,q1q2p1
  real(kind=DP) :: fzdd0,fzddd,opzthrdm2,omzthrdm2
  real(kind=DP) :: zeta2,zeta3,zeta4,eudrz,eudzz,eudrrz,eudrzz,eudzzz
  real(kind=DP) :: e1,e2,e3
  real(kind=DP) :: euda,eudb,eudaa,eudbb,eudab,eudaaa,eudaab,eudabb,eudbbb
  real(kind=DP) :: rsdn2,rsdn3,zetada2,zetada3,zetadb2,zetadb3
  real(kind=DP) :: gzz,gzzz,duda,dudb,dudaa,dudab,dudbb,dudaaa,dudaab,dudabb,dudbbb
  real(kind=DP) :: dtdu,dtduu,dtduuu,dtdnu,dtdnnu,dtdnuu,dtdug,dtduug,dtdugg,dtdnug,invg
  real(kind=DP) :: dtda,dtdb,dtdaa,dtdab,dtdag,dtdbg,dtdbb
  real(kind=DP) :: dtdaaa,dtdaab,dtdabb,dtdbbb,dtdaag,dtdagg,dtdbbg,dtdbgg,dtdabg
  real(kind=DP) :: invg3,m3u
  real(kind=DP) :: dAdu,dAduu,dAduuu,dAdeu,dAdeeu,dAdeuu
  real(kind=DP) :: dAda,dAdb,dAdaa,dAdab,dAdbb,dAdaaa,dAdaab,dAdabb,dAdbbb
  real(kind=DP) :: dHdf,dHdff,dHdft,dHdfff,dHdfft,dHdftt
  real(kind=DP) :: dHdu,dHduu,dHduuu,dHdfu,dHdtu,dHdffu,dHdfuu,dHdttu,dHdtuu,dHdftu
  real(kind=DP) :: tiu,siu2
  real(kind=DP) :: dHdb,dHdbb,dHdag,dHdbg,dHdab,dHdbbb,dHdaab,dHdabb,dHdbbg,dHdbgg,dHdaag,dHdagg,dHdabg
  
  real(kind=DP) :: dH1dr,dH1dt,dH1drr,dH1drt,dH1dtt,dH1drrr,dH1drrt,dH1drtt,dH1dttt,dH1drtu
  real(kind=DP) :: dH1du,dH1duu,dH1dru,dH1dtu,dH1duuu,dH1drru,dH1druu,dH1dtuu,dH1dttu
  real(kind=DP) :: r1rs2,r1rs3,r1rs4,ivt,ivt2,ivtmr1t,ivu,ivu2,r1t2,r1t22,r1t23
  real(kind=DP) :: dH1dn,dH1dg,dH1dnn,dH1dng,dH1dgg,dH1dnnn,dH1dnng,dH1dngg,dH1dggg
  real(kind=DP) :: dH1da,dH1db,dH1daa,dH1dab,dH1dbb,dH1dbg,dH1dag
  real(kind=DP) :: dH1daaa,dH1daab,dH1dabb,dH1dbbb,dH1dbbg,dH1dbgg,dH1dagg,dH1daag,dH1dabg

#ifdef GGA_PAW_REVISED
  real(kind=DP), parameter  :: eta = 1.d-12
  real(kind=DP), parameter  :: sixthm = -0.1666666666666666666666d0
#endif

  facw = nspin

#ifdef GGA_PAW_REVISED
  bet=0.06672455060314922d0
#else
  bet = xnu*cc0
#endif

  delt = 2.d0*alf/bet
  
     if ( nspin == 1 ) then
        g   = 1.d0
        g3  = g**3
        g4  = g3*g
        facpon = -delt/(g3*bet)

        do i = ist,ien,dnr ! MPI
           d  = facw*chgrhr_l(i, 1)
           if(d < density_minimum) cycle
! gradient of charge density 95/12/3 H.S.
           dd = facw* grad_trho(i)
           rs = (0.75/(PAI*d))**thrd
           fk = (3*PAI*PAI*d)**thrd
           sk = dsqrt(4*fk/PAI)
           if(d > density_minimum2) then
              t = dd/(d*sk*2)
              s = dd/(d*fk*2)
           else
              t = 0.d0
              s = 0.d0
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
           s2  = s*s
           t4  = t2*t2
           t6  = t4*t2
           rs2 = rs*rs
           rs3 = rs2*rs
           q4  = 1 + b*t2
           q5  = 1 + b*t2+b2*t4
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
           h0t = 2*bet*t*(1+2*b*t2)/q8 * g3
           h0b = -bet*t6*(2*b+b2*t2)/q8 * g3
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

           exc0 = eu*d/facw + ec1
           excd = ecd + ec1d
           dF_drho(i, 1) = dF_drho(i, 1) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(dabs(grad_trho(i)) > 1.0d-9) then
!              dF_dgradrho(i) = excdd / grad_trho(i)
              dF_dgradrho(i) = excdd
           else
              dF_dgradrho(i) = 0.d0
           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d/facw
!           exc1 = exc1 + exc0*wos(i)
           exc(i)=exc(i)+exc0
        
! ***** Calc heigher derivatives. *****
           
           rsdn=-thrd*rs/d
           rsdnn=-4.d0*thrd*rsdn/d
           rsdnnn=-thrd7*rsdnn/d
!print '(i2,5e19.10)',i,d,rs,rsdn,rsdnn,rsdnnn
!q3 = a*(b1/rs12+2*b2n+3*b3*rs12+2*b4*p1*rsp)   
           q1drr =0.5d0*a*( -b1/rs32+3.d0*b3/rs12+8.d0*b4 )
           q1drrr=0.75d0*a*( b1/rs - b3 )/rs32
!eurs = -2*a*a1*q2 - q0*q3/(q1**2 + q1)
           invq12pq1=1.d0/(q1**2 + q1)
           q12p1=2.d0*q1+1.d0
           eudrr=4.d0*a*a1*q3*invq12pq1-q0*( q1drr-q12p1*q3**2*invq12pq1 )*invq12pq1
           eudrrr=6.d0*a*a1*( q1drr - q12p1*q3**2*invq12pq1 ) * invq12pq1 &
                    -q0*( q1drrr -(2.d0*q3**3+3.d0*q12p1*q3*q1drr -  &
                    2.d0*q12p1**2*q3**3* invq12pq1)*invq12pq1 ) * invq12pq1
           eudn=eurs*rsdn
           eudnn=eudrr*rsdn**2+eurs*rsdnn
           eudnnn=eudrrr*rsdn**3+3.d0*eudrr*rsdn*rsdnn+eurs*rsdnnn
!print '(i2,5e19.12)',i,d,eu,eudn,eudnn,eudnnn
           ecdd=2.d0*eudn+d*eudnn
           ecddd=3.d0*eudnn+d*eudnnn
!print '(i2,5e19.12)',i,d,eu*d,ecdd,ecddd
           
           dtdn  = -sixth7*t/d
           dtdnn = -sixth13*dtdn/d
           dtdnnn= -sixth19*dtdnn/d
!print '(i2,5e19.12)',i,d,t,dtdn,dtdnn,dtdnnn
           if(dabs(grad_trho(i)) > 1.0d-9) then
               dtdg  = t/dd
               dtdng = -sixth7*dtdg/d
               dtdnng= -sixth13*dtdng/d
           else
               dtdg  = 0.d0
               dtdng = 0.d0
               dtdnng= 0.d0
           endif
           dtdngg = 0.d0
           dtdggg = 0.d0
           
           dAde   = b * ( b + delt ) / bet
           dAdee  = 2.d0 * dAde * ( b + 0.5d0 * delt ) / bet
           dAdeee = 2.d0 * ( dAde**2 + ( b + 0.5d0 * delt )*dAdee ) /bet
           
           dAdn   = dAde  *eudn
           dAdnn  = dAdee *eudn**2  + dAde*eudnn
           dAdnnn = dAdeee*eudn**3  + 3.d0*dAdee*eudn*eudnn + dAde*eudnnn   
!print '(i2,5e19.10)' ,i,d,b,dAdnn
           
           at2  = b * t2
           a2t4 = at2**2
           a3t6 = at2*a2t4
           t3   = t2*t
           t5   = t4*t
           t7   = t6*t
           t8   = t7*t
           fE    = q5**2 + delt*t2*q4*q5
           fEdt  = 4.d0*b*t*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt *2.d0*t*(1.d0 + 4.d0*at2 + 6.d0*a2t4 + 4.d0*a3t6 )
           fEdtt = 4.d0*b*( 1.d0 + 9.d0*at2 + 15.d0*a2t4 + 14.d0*a3t6 ) + &
                    delt*( 2.d0 + 24.d0*at2 + 60.d0*a2t4 + 56.d0*a3t6 )
           fEda  = 2.d0*t2*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt*t4*( 2.d0 + 4.d0*at2 + 3.d0*a2t4 )
           fEdaa = 6.d0*t4*( 2.d0*q5 - 1.d0) + &
                    delt*2.d0*t6*( 2.d0 + 3.d0*at2 )
           fF    = t*( 1.d0 + 2.d0*at2 )         
           fFdt  = 1.d0 + 6.d0*at2
           fFdtt = 12.d0*b*t
           fFda  = 2.d0*t3
           fFdaa = 0.d0
           fG    = 2.d0*b*t6 + b2*t8
           fGdt  = 12.d0*b*t5 + 8.d0*b2*t7
           fGdtt = 60.d0*b*t4 + 56.d0*b2*t6
           fGda  = 2.d0*t6 + 2.d0*b*t8
           fGdaa = 2.d0*t8
           
           betivE=bet/fE
           betivE2=betivE/fE
           dHda   =     -fG*betivE
           dHdt   = 2.d0*fF*betivE
           dHdaa  =     -( fGda*fE-fG*fEda )*betivE2
           dHdat  =     -( fGdt*fE-fG*fEdt )*betivE2
           dHdtt  = 2.d0*( fFdt*fE-fF*fEdt )*betivE2
           dHdaaa =     -( fGdaa*fE-fG*fEdaa-2.d0*( fGda*fE-fG*fEda )*fEda/fE )*betivE2
           dHdaat = 2.d0*( fFdaa*fE-fF*fEdaa-2.d0*( fFda*fE-fF*fEda )*fEda/fE )*betivE2
           dHdatt =     -( fGdtt*fE-fG*fEdtt-2.d0*( fGdt*fE-fG*fEdt )*fEdt/fE )*betivE2
           dHdttt = 2.d0*( fFdtt*fE-fF*fEdtt-2.d0*( fFdt*fE-fF*fEdt )*fEdt/fE )*betivE2
           
           dHdn     = dHda  *dAdn       + dHdt*dtdn
           dHdg     = dHdt  *dtdg
           dHdnn    = dHdaa *dAdn**2    + 2.d0*dHdat*dAdn*dtdn  + dHdtt*dtdn**2   + dHda*dAdnn  + dHdt*dtdnn
           dHdng    = dHdat *dAdn*dtdg  + dHdtt*dtdn*dtdg       + dHdt*dtdng
           dHdgg    = dHdtt*dtdg**2
           dHdnnn   = dHdaaa*dAdn**3    + dHdttt*dtdn**3 &
                        + 3.d0*dHdaat*dAdn**2*dtdn  + 3.d0*dHdatt*dAdn*dtdn**2 &
                        + 3.d0*dHdaa*dAdn*dAdnn + 3.d0*dHdtt*dtdn*dtdnn + 3.d0*dHdat*(dAdnn*dtdn + dAdn*dtdnn) &
                        + dHda*dAdnnn   + dHdt*dtdnnn
           dHdnng   = dHdaat*dAdn**2*dtdg   + 2.d0*dHdatt*dAdn*dtdn*dtdg    + dHdttt*dtdn**2*dtdg &
                        + dHdat*(dAdnn*dtdg + 2.d0*dAdn*dtdng)  + dHdtt*(dtdnn*dtdg + 2.d0*dtdn*dtdng) &
                        + dHdt*dtdnng
           dHdngg   = dHdatt*dAdn*dtdg**2   + dHdttt*dtdn*dtdg**2   + 2.d0*dHdtt*dtdg*dtdng
           dHdggg   = dHdttt*dtdg**3
           
           dq6   = c2+2.d0*c3*rs
           ddq6  = 2.d0*c3
           dq7   = c4+2.d0*c5*rs+3.d0*c6*rs2
           ddq7  = 2.d0*c5+6.d0*c6*rs
           dddq7 = 6.d0*c6
           
           iq7=1.d0/q7
           iq72=iq7*iq7
           iq73=iq72*iq7
           iq74=iq73*iq7
           
           ccrs2 =  ddq6*iq7-(2.d0*dq6*dq7+q6*ddq7)*iq72+2.d0*q6*dq7**2*iq73
           ccrs3 =  -iq72*(3.d0*(ddq6*dq7+dq6*ddq7)+q6*dddq7) + &
                                    6.d0*dq7*(dq6*dq7+q6*ddq7)*iq73 - &
                                    6.d0*q6*dq7**3*iq74
           r1rs2=r1rs*r1rs
           r1rs3=r1rs2*r1rs
           r1rs4=r1rs3*r1rs
           
           if(dabs(grad_trho(i)) > 1.0d-9) then
               ivt=1.d0/t
               ivt2=ivt**2
               ivtmr1t=ivt-r1*t
           else
               ivt=0.d0
               ivt2=0.d0
               ivtmr1t=0.d0
           end if
           
           dH1dr=h1rs
           dH1dt=h1t                     
           dH1drr  = xnu*ccrs2*t2*r3-2.d0*h1rs*r1rs*t2-h1*r1rs2*t4
           dH1drt  = 2.d0*h1rs*ivtmr1t-2.d0*h1*t*r1rs
           dH1dtt  = 2.d0*h1t*ivtmr1t-2.d0*h1*(ivt2+r1)
           dH1drrr = xnu*ccrs3*t2*r3-3.d0*t2*r1rs*(dH1drr+h1rs*r1rs*t2)-h1*r1rs3*t6
           dH1drrt = 2.d0*dH1drr*ivtmr1t-4.d0*h1rs*r1rs*t
           dH1drtt = 2.d0*dH1drt*ivtmr1t-2.d0*h1t*r1rs*t-2.d0*ivt*h1rs*(ivt+r1*t)-2.d0*h1*r1rs
           dH1dttt = 2.d0*dH1dtt*ivtmr1t-4.d0*h1t*(ivt2+r1)+4.d0*h1*ivt2*ivt
           
           dH1dn     = dH1dr  *rsdn       + dH1dt*dtdn
           dH1dg     = dH1dt  *dtdg
           dH1dnn    = dH1drr *rsdn**2    + 2.d0*dH1drt*rsdn*dtdn  + dH1dtt*dtdn**2   + dH1dr*rsdnn  + dH1dt*dtdnn
           dH1dng    = dH1drt *rsdn*dtdg  + dH1dtt*dtdn*dtdg       + dH1dt*dtdng
           dH1dgg    = dH1dtt*dtdg**2
           dH1dnnn   = dH1drrr*rsdn**3    + dH1dttt*dtdn**3 &
                        + 3.d0*dH1drrt*rsdn**2*dtdn  + 3.d0*dH1drtt*rsdn*dtdn**2 &
                        + 3.d0*dH1drr*rsdn*rsdnn + 3.d0*dH1dtt*dtdn*dtdnn + 3.d0*dH1drt*(rsdnn*dtdn + rsdn*dtdnn) &
                        + dH1dr*rsdnnn   + dH1dt*dtdnnn
           dH1dnng   = dH1drrt*rsdn**2*dtdg   + 2.d0*dH1drtt*rsdn*dtdn*dtdg    + dH1dttt*dtdn**2*dtdg &
                        + dH1drt*(rsdnn*dtdg + 2.d0*rsdn*dtdng)  + dH1dtt*(dtdnn*dtdg + 2.d0*dtdn*dtdng) &
                        + dH1dt*dtdnng
           dH1dngg   = dH1drtt*rsdn*dtdg**2   + dH1dttt*dtdn*dtdg**2   + 2.d0*dH1dtt*dtdg*dtdng
           dH1dggg   = dH1dttt*dtdg**3
           
!dF_drho(i, 1)=h+d*dHdn+ecd
           dFc_daa(i)   = 2.d0*(dHdn + dH1dn) + d*(dHdnn + dH1dnn)       + ecdd
           dFc_dag(i)   = dHdg + dH1dg + d*(dHdng + dH1dng)
           dFc_dgg(i)   = d*(dHdgg + dH1dgg)
           dFc_daaa(i)  = 3.d0*(dHdnn + dH1dnn)    + d*(dHdnnn + dH1dnnn)  + ecddd
           dFc_daag(i)  = 2.d0*(dHdng + dH1dng)    + d*(dHdnng + dH1dnng)
           dFc_dagg(i)  = dHdgg + dH1dgg        + d*(dHdngg + dH1dngg)
           dFc_dggg(i)  = d*(dHdggg + dH1dggg)  
        end do
     else if ( nspin ==  2 ) then
        thrd2 = thrd * 2
        thrd4 = thrd * 4
        fzdd = 8/(9*(2**thrd4 - 2 ))

        do i = ist,ien,dnr
           d = chgrhr_l(i, 1) + chgrhr_l(i, nspin)
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
           g4  = g3*g
           facpon = -delt/(g3*bet)
           if(d < density_minimum) cycle
           dd = grad_trho(i)
           rs = (0.75d0/(PAI*d))**thrd
           fk = (3*PAI*PAI*d)**thrd
           sk = dsqrt(4*fk/PAI)

           if(d > density_minimum2) then
           ! ****** aas mkatsu *****
!              t = dd/(d*sk*2)
              t = dd/(d*sk*2)/g         !    mkatsu aas
              s = dd/(d*fk*2)
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

!           zetadxd = - 2 * chgrhr_l(i, 3-is) * (-1.d0)**is / d
           zetadxda = - 2 * chgrhr_l(i, 2) * (-1.d0) / d
           zetadxdb = - 2 * chgrhr_l(i, 1) / d
!           ecd = eu-thrd*rs*eurs + euzt * zetadxd
           ecda = eu-thrd*rs*eurs + euzt * zetadxda 
           ecdb = eu-thrd*rs*eurs + euzt * zetadxdb

#ifndef GGA_PAW_REVISED
           if(onpzeta < zeta_minimum2 .or. onmzeta < zeta_minimum2) then
              ec1  = 0.d0;    ec1d = 0.d0;     ec1dd= 0.d0
              ec1da = 0.d0;     ec1db = 0.d0
              go to 10
           end if
#endif

           pon = facpon*eu
           b   = delt/(dexp(pon)-1)
           b2  = b*b
           t2  = t*t
           s2  = s*s
           t4  = t2*t2
           t6  = t4*t2
           rs2 = rs*rs
           rs3 = rs2*rs
           q4  = 1+b*t2
           q5  = 1+b*t2+b2*t4
           q6  = c1+c2*rs+c3*rs2
           q7  = 1.d0+c4*rs+c5*rs2+c6*rs3
           cc  = -cx + q6/q7
           r0  = (sk/fk)**2
           r1  = a4*r0*g4
           coeff = cc-cc0-3.d0*cx/7.d0
           r2  = xnu*coeff*g3
!           a4ms = -a4*g4*s2                                !aas
           a4ms = -r1*t2                                    !aas
           r3  = dexp(a4ms)                       
           h0  = g3*(bet/delt)*dlog(1+delt*q4*t2/q5)
           h1  = r3*r2*t2
           h   = h0+h1
           
           q8  = q5*q5+delt*q4*q5*t2
           h0t = 2*bet*t*(1+2*b*t2)/q8 * g3
           h0b = -bet*t6*(2*b+b2*t2)/q8 * g3
!           h0rs = h0b * b2 * eurs / bet / g3               !aas
           h0rs = h0b * b * eurs * (b+delt) / bet / g3      !aas
           h1t  = 2*r3*r2*t*(1.d0-r1*t2)
           ccrs = (c2+2*c3*rs)/q7-q6*(c4+2*c5*rs+3*c6*rs2)/q7**2
!           r1rs = 100*r0/rs                                ! aas
           r1rs = r1/rs                                     ! aas
           h1rs = xnu*t2*r3*(ccrs - coeff*t2*r1rs) * g3
           
           ht   = h0t+h1t
           hrs  = h0rs + h1rs

           if(dabs(zeta) < zeta_minimum) then
              gzt = -2.d0/9*zeta *( 1 +14.d0/27 *zeta*zeta &
                   &             *( 1 +13.d0/18 *zeta*zeta ))
           else
#ifdef GGA_PAW_REVISED
              gzt = ( (onpzeta**2+eta)**sixthm - (onmzeta**2+eta)**sixthm)/3.d0
#else
              gzt = ( 1.d0/onpzeta**thrd -1.d0/onmzeta**thrd )/3.d0
#endif

           end if
!           bzt = pon*b**2 * (euzt-3*eu*gzt/g) / ( bet * g3 )           !aas
           bzt = b * ( b+delt) * ( euzt - 3*eu*gzt/g ) / ( bet * g3 )   !aas
           h0zt = 3 * gzt * h0 / g + h0b * bzt
!           h1zt =(3*gzt + 4*gzt*a4ms ) * h1 / g                        !aas
           h1zt = (3.d0 + 4.d0*a4ms) * h1 * gzt / g                    !aas

           ec1  = d*h * 0.5
! --> revised by Kato-san,  1999Jul16
!!$           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + h0zt * zetadxd
!           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxd
           ec1da = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*ht/g + h1zt ) * zetadxda
           ec1db = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*ht/g + h1zt ) * zetadxdb
           ec1dd = 0.5*ht/sk / g

10         continue
           exc0 = eu*d * 0.5 + ec1
!           excd = ecd + ec1d
           excda= ecda + ec1da
           excdb= ecdb + ec1db
!           dF_drho(i, is) = dF_drho(i, is) + excd
           dF_drho(i, 1) = dF_drho(i, 1) + excda
           dF_drho(i, 2) = dF_drho(i, 2) + excdb
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
!           if(is ==  2) then
              if(dabs(grad_trho(i)) > 1.0d-9) then
!                 dF_dgradrho(i) = excdd / grad_trho(i)
                 dF_dgradrho(i) = excdd
              else
                 dF_dgradrho(i) = 0.d0
              endif
!           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d * 0.5
!           exc1 = exc1 + exc0*wos(i)
           exc(i)=exc(i)+exc0*2.d0
!print '(i3,5e19.12)', i,exc0,dF_drho(i, 1:2),dF_dgradrho(i),excdd/ grad_trho(i)
      
! ***** Calc heigher order derivatives. *****

           rsdn=-thrd*rs/d
           rsdnn=-4.d0*thrd*rsdn/d
           rsdnnn=-thrd7*rsdnn/d
           
           id2=1.d0/d/d
           id3=id2/d
           id4=id3/d
           zetada   =     zetadxda/d
           zetadb   =     zetadxdb/d
           zetadaa  =    -4.d0*chgrhr_l(i, 2)*id3
           zetadab  =     2.d0*zeta*id2
           zetadbb  =     4.d0*chgrhr_l(i, 1)*id3
           zetadaaa =    12.d0*chgrhr_l(i, 2)*id4
           zetadaab =   ( -4.d0*chgrhr_l(i, 1) + 8.d0*chgrhr_l(i, 2))*id4
           zetadabb =   ( -8.d0*chgrhr_l(i, 1) + 4.d0*chgrhr_l(i, 2))*id4
           zetadbbb =   -12.d0*chgrhr_l(i, 1)*id4
           
!print '(i2,7d19.10)', i,chgrhr_l(i,1:2),zetadab,zetadaaa,zetadaab,zetadabb,zetadbbb
             
           q1drr =0.5d0*a*( -b1/rs32+3.d0*b3/rs12+8.d0*b4 )
           q1drrr=0.75d0*a*( b1/rs - b3 )/rs32
           invq12pq1=1.d0/(q1**2 + q1)
           q12p1=2.d0*q1+1.d0
           e1dr=-2.d0*a*a1*q2-q0*q3*invq12pq1
           e1drr=4.d0*a*a1*q3*invq12pq1-q0*( q1drr-q12p1*q3**2*invq12pq1 )*invq12pq1
           e1drrr=6.d0*a*a1*( q1drr - q12p1*q3**2*invq12pq1 ) * invq12pq1 &
                    -q0*( q1drrr -(2.d0*q3**3+3.d0*q12p1*q3*q1drr -  &
                    2.d0*q12p1**2*q3**3* invq12pq1)*invq12pq1 ) * invq12pq1
                    
           q1pdrr =0.5d0*ap*( -b1p/rs32+3.d0*b3p/rs12+8.d0*b4p )
           q1pdrrr=0.75d0*ap*( b1p/rs - b3p )/rs32
           invq1p2pq1p=1.d0/(q1p**2 + q1p)
           q1p2p1=2.d0*q1p+1.d0
           e1pdr=-2.d0*ap*a1p*q2p-q0p*q3p*invq1p2pq1p
           e1pdrr=4.d0*ap*a1p*q3p*invq1p2pq1p-q0p*( q1pdrr-q1p2p1*q3p**2*invq1p2pq1p )*invq1p2pq1p
           e1pdrrr=6.d0*ap*a1p*( q1pdrr - q1p2p1*q3p**2*invq1p2pq1p ) * invq1p2pq1p &
                    -q0p*( q1pdrrr -(2.d0*q3p**3+3.d0*q1p2p1*q3p*q1pdrr -  &
                    2.d0*q1p2p1**2*q3p**3* invq1p2pq1p)*invq1p2pq1p ) * invq1p2pq1p
                    
           q1qdrr =0.5d0*aq*( -b1q/rs32+3.d0*b3q/rs12+8.d0*b4q )
           q1qdrrr=0.75d0*aq*( b1q/rs - b3q )/rs32
           invq1q2pq1q=1.d0/(q1q**2 + q1q)
           q1q2p1=2.d0*q1q+1.d0
           e1qdr=-2.d0*aq*a1q*q2q-q0q*q3q*invq1q2pq1q
           e1qdrr=4.d0*aq*a1q*q3q*invq1q2pq1q-q0q*( q1qdrr-q1q2p1*q3q**2*invq1q2pq1q )*invq1q2pq1q
           e1qdrrr=6.d0*aq*a1q*( q1qdrr - q1q2p1*q3q**2*invq1q2pq1q ) * invq1q2pq1q &
                    -q0q*( q1qdrrr -(2.d0*q3q**3+3.d0*q1q2p1*q3q*q1qdrr -  &
                    2.d0*q1q2p1**2*q3q**3* invq1q2pq1q)*invq1q2pq1q ) * invq1q2pq1q
           
           zeta2=zeta**2
           zeta3=zeta2*zeta
           zeta4=zeta2*zeta2
           
           if(dabs(zeta) < zeta_minimum) then
              fzdd0 = fzdd &
                   &            *( 1 + 5.d0*ninth *zeta2 &
                   &            *( 1 +22.d0*twntysvnth *zeta2 ))
              fzddd = fzdd*10.d0*ninth*zeta*( 1.d0 + 44.d0*twntysvnth*zeta2 )
              gzz   = -2.d0*ninth*( 1.d0 + 14.d0*ninth*zeta2*( 1.d0 + 65.d0*fftyfrth*zeta2 ) )
              gzzz  = -56.d0*eityonth*zeta*( 1.d0 + 65.d0*twntysvnth*zeta2 )
           else
              opzthrdm2=1.d0/onpzeta**(thrd2)
              omzthrdm2=1.d0/onmzeta**(thrd2)
              fzdd0 =  fzdd*0.5d0*( opzthrdm2 + omzthrdm2 )
              fzddd = -fzdd*thrd*( opzthrdm2/onpzeta - omzthrdm2/onmzeta )
              gzz   = -ninth*( opzthrdm2**2 + omzthrdm2**2 )
!              gzz   = -ninth*( 1.d0/onpzeta**(thrd4) + 1.d0/onmzeta**(thrd4) )
              gzzz  = 4.d0*twntysvnth*( opzthrdm2**2/onpzeta - omzthrdm2**2/onmzeta )
           end if
           
           e1=q0*q2
           e2=q0p*q2p
           e3=q0q*q2q
           
           eudrr    = e1drr  + ( e1pdrr  -  e1drr  ) * fz * zeta4 - e1qdrr  * fz * onzt4 / fzdd
           eudrrr   = e1drrr + ( e1pdrrr -  e1drrr ) * fz * zeta4 - e1qdrrr * fz * onzt4 / fzdd
           eudzz    = ( e2 - e1 )*( fzdd0*zeta4 + 8.d0*fzd*zeta3 + 12.d0*fz*zeta2 ) &
                        - e3*( fzdd0*onzt4 - 8.d0*fzd*zeta3 - 12.d0*fz*zeta2 )/fzdd
           eudzzz   = ( e2 - e1 )*( fzddd*zeta4 + 12.d0*fzdd0*zeta3 + 36.d0*fzd*zeta2 + 24.d0*fz*zeta ) &
                        - e3*( fzddd*onzt4 - 12.d0*fzdd0*zeta3 - 36.d0*fzd*zeta2 - 24.d0*fz*zeta )/fzdd
           eudrz    = ( e1pdr - e1dr )*( fzd*zeta4 + 4.d0*fz*zeta3 ) &
                        - e1qdr *( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd 
           eudrrz   = ( e1pdrr  - e1drr  )*( fzd*zeta4 + 4.d0*fz*zeta3 ) &
                        - e1qdrr *( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd
           eudrzz   = ( e1pdr - e1dr )*( fzdd0*zeta4 + 8.d0*fzd*zeta3 + 12.d0*fz*zeta2 ) &
                        - e1qdr *( fzdd0*onzt4 - 8.d0*fzd*zeta3 - 12.d0*fz*zeta2 )/fzdd
                        
           rsdn2    = rsdn*rsdn
           rsdn3    = rsdn2*rsdn
           zetada2  = zetada**2
           zetada3  = zetada2*zetada
           zetadb2  = zetadb**2
           zetadb3  = zetadb2*zetadb
                        
           euda     = eurs*rsdn + euzt*zetada
           eudb     = eurs*rsdn + euzt*zetadb
           eudaa    = eudrr*rsdn2 + 2.d0*eudrz*rsdn*zetada + eudzz*zetada2 + eurs*rsdnn + euzt*zetadaa
           eudab    = eudrr*rsdn2 + eudrz*rsdn*( zetada + zetadb ) + eudzz*zetada*zetadb &
                        + eurs*rsdnn + euzt*zetadab
           eudbb    = eudrr*rsdn2 + 2.d0*eudrz*rsdn*zetadb + eudzz*zetadb2 + eurs*rsdnn + euzt*zetadbb
           eudaaa   = eudrrr*rsdn3 &
                        + 3.d0*eudrrz*rsdn2*zetada &
                        + 3.d0*eudrzz*rsdn*zetada2 + eudzzz*zetada3 &
                        + 3.d0*eudrr*rsdn*rsdnn + 3.d0*eudzz*zetada*zetadaa &
                        + 3.d0*eudrz*( rsdnn*zetada + rsdn*zetadaa ) &
                        + eurs*rsdnnn + euzt*zetadaaa
           eudaab   = eudrrr*rsdn3 &
                        + eudrrz*rsdn2*( zetadb + 2.d0*zetada ) &
                        + eudrzz*rsdn*zetada*( zetada + 2.d0*zetadb ) &
                        + eudzzz*zetada2*zetadb &
                        + 3.d0*eudrr*rsdn*rsdnn &
                        + eudrz*( rsdnn*( 2.d0*zetada + zetadb ) + rsdn*( 2.d0*zetadab + zetadaa ) ) &
                        + eudzz*( 2.d0*zetada*zetadab + zetadaa*zetadb ) &
                        + eurs*rsdnnn + euzt*zetadaab
           eudabb   = eudrrr*rsdn3 &
                        + eudrrz*rsdn2*( zetada + 2.d0*zetadb ) &
                        + eudrzz*rsdn*zetadb*( zetadb + 2.d0*zetada ) &
                        + eudzzz*zetada*zetadb2 &
                        + 3.d0*eudrr*rsdn*rsdnn &
                        + eudrz*( rsdnn*( zetada + 2.d0*zetadb ) + rsdn*( 2.d0*zetadab + zetadbb ) ) &
                        + eudzz*( zetada*zetadbb + 2.d0*zetadb*zetadab ) &
                        + eurs*rsdnnn + euzt*zetadabb
           eudbbb   = eudrrr*rsdn3 &
                        + 3.d0*eudrrz*rsdn2*zetadb &
                        + 3.d0*eudrzz*rsdn*zetadb2 + eudzzz*zetadb3 &
                        + 3.d0*eudrr*rsdn*rsdnn + 3.d0*eudzz*zetadb*zetadbb &
                        + 3.d0*eudrz*( rsdnn*zetadb + rsdn*zetadbb ) &
                        + eurs*rsdnnn + euzt*zetadbbb
                        
!print '(i2,12d19.10)', i,chgrhr_l(i,1:2),eu,euda,eudb,eudaa,eudbb,eudab,eudaaa,eudaab,eudabb,eudbbb
!print *,eu,     e1+(e2-e1)* fz * zeta4 - e3  * fz * onzt4 / fzdd
!print *,eurs,   e1dr  + ( e1pdr  -  e1dr  ) * fz * zeta4 - e1qdr  * fz * onzt4 / fzdd
!print *,euzt,   (e2-e1)*( fzd*zeta4 + 4.d0*fz*zeta3 )-e3*( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd 

           duda     = gzt*zetada
           dudb     = gzt*zetadb 
           dudaa    = gzz*zetada2     + gzt*zetadaa
           dudab    = gzz*zetada*zetadb + gzt*zetadab
           dudbb    = gzz*zetadb2     + gzt*zetadbb
           dudaaa   = gzzz*zetada3    + 3.d0*gzz*zetada*zetadaa   + gzt*zetadaaa
           dudaab   = gzzz*zetada2*zetadb + gzz*( zetadaa*zetadb + 2.d0*zetada*zetadab ) + gzt*zetadaab
           dudabb   = gzzz*zetada*zetadb2 + gzz*( zetada*zetadbb + 2.d0*zetadab*zetadb ) + gzt*zetadabb
           dudbbb   = gzzz*zetadb3    + 3.d0*gzz*zetadb*zetadbb   + gzt*zetadbbb 
           
!print '(i2,12d19.10)', i,chgrhr_l(i,1:2),g,duda,dudb,dudaa,dudbb,dudab,dudaaa,dudaab,dudabb,dudbbb

           
           dtdn  = -sixth7*t/d
           dtdnn = -sixth13*dtdn/d
           dtdnnn= -sixth19*dtdnn/d
!print '(i2,5e19.12)',i,d,t,dtdn,dtdnn,dtdnnn
           if(dabs(grad_trho(i)) > 1.0d-9) then
               dtdg  = t/dd
               dtdng = -sixth7*dtdg/d
               dtdnng= -sixth13*dtdng/d
!               dtdnug = -dtdn*invg/dd
           else
               dtdg  = 0.d0
               dtdng = 0.d0
               dtdnng= 0.d0
               dtdnug= 0.d0 
           endif
           dtdgg  = 0.d0
           dtdngg = 0.d0
           dtdggg = 0.d0
           invg=1.d0/g
           dtdu   = -t*invg
           dtduu  = -2.d0*dtdu*invg 
           dtduuu = -3.d0*dtduu*invg
           dtdnu  = -dtdn*invg
           dtdnnu = -dtdnn*invg
           dtdnuu = -sixth7*dtduu/d
           dtdug  = -dtdg*invg
           dtduug = 2.d0*dtdg*invg**2
           dtdugg = 0.d0
           if(dabs(grad_trho(i)) > 1.0d-9) then
                dtdnug = dtdnu/dd
           else
                dtdnug = 0.d0
           end if
           
           dtda   = dtdn + dtdu*duda
           dtdb   = dtdn + dtdu*dudb
           dtdaa  = dtdnn + 2.d0*dtdnu*duda + dtduu*duda**2 + dtdu*dudaa
           dtdbb  = dtdnn + 2.d0*dtdnu*dudb + dtduu*dudb**2 + dtdu*dudbb
           dtdab  = dtdnn + dtdnu*( duda + dudb ) + dtduu*duda*dudb + dtdu*dudab
           dtdag  = dtdng + dtdug*duda
           dtdbg  = dtdng + dtdug*dudb
           dtdaaa = dtdnnn+ 3.d0*dtdnnu*duda + 3.d0*dtdnuu*duda**2 + dtduuu*duda**3 &
                    + 3.d0*dtduu*duda*dudaa + 3.d0*dtdnu*dudaa + dtdu*dudaaa
           dtdbbb = dtdnnn+ 3.d0*dtdnnu*dudb + 3.d0*dtdnuu*dudb**2 + dtduuu*dudb**3 &
                    + 3.d0*dtduu*dudb*dudbb + 3.d0*dtdnu*dudbb + dtdu*dudbbb
           dtdaab = dtdnnn+ dtdnnu*( 2.d0*duda + dudb ) + dtdnuu*duda*( duda + 2.d0*dudb ) + dtduuu*duda**2*dudb &
                    + dtduu*( dudaa*dudb + 2.d0*duda*dudab ) + dtdnu*( dudaa + 2.d0*dudab ) + dtdu*dudaab
           dtdabb = dtdnnn+ dtdnnu*( duda + 2.d0*dudb ) + dtdnuu*dudb*( dudb + 2.d0*duda ) + dtduuu*dudb**2*duda &
                    + dtduu*( dudbb*duda + 2.d0*dudb*dudab ) + dtdnu*( dudbb + 2.d0*dudab ) + dtdu*dudabb 
           dtdaag = dtdnng+ 2.d0*dtdnug*duda + dtduug*duda**2 + dtdug*dudaa
           dtdbbg = dtdnng+ 2.d0*dtdnug*dudb + dtduug*dudb**2 + dtdug*dudbb
           dtdabg = dtdnng+ dtdnug*( duda + dudb ) + dtduug*duda*dudb + dtdug*dudab 
           
!print '(i2,20d19.10)', i,t,dtda,dtdb,dtdg,dtdaa,dtdbb,dtdgg,dtdab,dtdag,dtdbg,dtdaaa,dtdbbb,dtdggg,dtdaab,dtdabb,dtdaag,dtdagg,dtdbbg,dtdbgg,dtdabg     
!print '(i2,12d19.10)', i,t,dtdag,dtdbg,dtdaag,dtdbbg,dtdabg     

           invg3  = 1.d0  /g3
           dAde   = b * ( b + delt ) / bet * invg3
           dAdee  = 2.d0 * dAde * ( b + 0.5d0 * delt ) / bet *invg3
           dAdeee = 2.d0 * ( dAde**2 + ( b + 0.5d0 * delt )*dAdee ) /bet * invg3
           m3u    = -3.d0*invg
           dAdu     = m3u*eu*dAde
           dAdeu    = m3u*( dAde + eu*dAdee )
           dAduu    = m3u*eu*invg*( -dAde + g*dAdeu )
           dAdeeu   = m3u*( 2.d0*dAdee + eu*dAdeee )
           dAdeuu   = m3u*invg*( -dAde + g*dAdeu - eu*dAdee + g*eu*dAdeeu )
           dAduuu   = m3u*eu*invg**2*( 2.d0*dAde -2.d0*g*dAdeu + g**2*dAdeuu )
           
           dAda     = dAde*euda + dAdu*duda
           dAdb     = dAde*eudb + dAdu*dudb
           dAdaa    = dAdee*euda**2 + 2.d0*dAdeu*euda*duda + dAduu*duda**2 + dAde*eudaa + dAdu*dudaa 
           dAdab    = dAdee*euda*eudb + dAdeu*( euda*dudb + eudb*duda ) + dAduu*duda*dudb + dAde*eudab + dAdu*dudab 
           dAdbb    = dAdee*eudb**2 + 2.d0*dAdeu*eudb*dudb + dAduu*dudb**2 + dAde*eudbb + dAdu*dudbb
           dAdaaa   = dAdeee*euda**3 + 3.d0*dAdeeu*euda**2*duda + 3.d0*dAdeuu*euda*duda**2 + dAduuu*duda**3 &
                        + 3.d0*dAdee*euda*eudaa + 3.d0*dAduu*duda*dudaa &
                        + 3.d0*dAdeu*( euda*dudaa + eudaa*duda ) + dAde*eudaaa + dAdu*dudaaa
           dAdaab   = dAdeee*euda**2*eudb + dAduuu*duda**2*dudb &
                        + dAdeeu*euda*( euda*dudb + 2.d0*eudb*duda ) &
                        + dAdeuu*( 2.d0*euda*dudb + eudb*duda )*duda &
                        + dAdee*( 2.d0*euda*eudab + eudb*eudaa ) + dAduu*( 2.d0*duda*dudab + dudb*dudaa ) &
                        + dAdeu*( 2.d0*eudab*duda + 2.d0*euda*dudab + dudb*eudaa + eudb*dudaa ) &
                        + dAde*eudaab + dAdu*dudaab
           dAdabb   = dAdeee*euda*eudb**2 + dAduuu*dudb**2*duda &
                        + dAdeeu*eudb*( eudb*duda + 2.d0*euda*dudb ) &
                        + dAdeuu*( 2.d0*eudb*duda + euda*dudb )*dudb &
                        + dAdee*( eudbb*euda + 2.d0*eudb*eudab ) + dAduu*( dudbb*duda + 2.d0*dudb*dudab ) &
                        + dAdeu*( eudbb*duda + 2.d0*eudb*dudab + 2.d0*dudb*eudab + dudbb*euda ) &
                        + dAde*eudabb + dAdu*dudabb
           dAdbbb   = dAdeee*eudb**3 + 3.d0*dAdeeu*eudb**2*dudb + 3.d0*dAdeuu*eudb*dudb**2 + dAduuu*dudb**3 &
                        + 3.d0*dAdee*eudb*eudbb + 3.d0*dAduu*dudb*dudbb &
                        + 3.d0*dAdeu*( eudb*dudbb + eudbb*dudb ) + dAde*eudbbb + dAdu*dudbbb 
                        
!print '(i2,12d19.10)', i,b,dAda,dAdb,dAdaa,dAdab,dAdbb,dAdaaa,dAdaab,dAdabb,dAdbbb

           at2  = b * t2
           a2t4 = at2**2
           a3t6 = at2*a2t4
           t3   = t2*t
           t5   = t4*t
           t7   = t6*t
           t8   = t7*t
           fE    = q5**2 + delt*t2*q4*q5
           fEdt  = 4.d0*b*t*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt *2.d0*t*(1.d0 + 4.d0*at2 + 6.d0*a2t4 + 4.d0*a3t6 )
           fEdtt = 4.d0*b*( 1.d0 + 9.d0*at2 + 15.d0*a2t4 + 14.d0*a3t6 ) + &
                    delt*( 2.d0 + 24.d0*at2 + 60.d0*a2t4 + 56.d0*a3t6 )
           fEda  = 2.d0*t2*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt*t4*( 2.d0 + 4.d0*at2 + 3.d0*a2t4 )
           fEdaa = 6.d0*t4*( 2.d0*q5 - 1.d0) + &
                    delt*2.d0*t6*( 2.d0 + 3.d0*at2 )
           fF    = t*( 1.d0 + 2.d0*at2 )         
           fFdt  = 1.d0 + 6.d0*at2
           fFdtt = 12.d0*b*t
           fFda  = 2.d0*t3
           fFdaa = 0.d0
           fG    = 2.d0*b*t6 + b2*t8
           fGdt  = 12.d0*b*t5 + 8.d0*b2*t7
           fGdtt = 60.d0*b*t4 + 56.d0*b2*t6
           fGda  = 2.d0*t6 + 2.d0*b*t8
           fGdaa = 2.d0*t8
           
           betivE=bet/fE*g3
           betivE2=betivE/fE
           dHdf   =     -fG*betivE
           dHdt   = 2.d0*fF*betivE
           dHdff  =     -( fGda*fE-fG*fEda )*betivE2
           dHdft  =     -( fGdt*fE-fG*fEdt )*betivE2
           dHdtt  = 2.d0*( fFdt*fE-fF*fEdt )*betivE2
           dHdfff =     -( fGdaa*fE-fG*fEdaa-2.d0*( fGda*fE-fG*fEda )*fEda/fE )*betivE2
           dHdfft = 2.d0*( fFdaa*fE-fF*fEdaa-2.d0*( fFda*fE-fF*fEda )*fEda/fE )*betivE2
           dHdftt =     -( fGdtt*fE-fG*fEdtt-2.d0*( fGdt*fE-fG*fEdt )*fEdt/fE )*betivE2
           dHdttt = 2.d0*( fFdtt*fE-fF*fEdtt-2.d0*( fFdt*fE-fF*fEdt )*fEdt/fE )*betivE2
           
           tiu    = 3.d0*invg
           siu2   = tiu*2.d0*invg
           dHdu   = tiu*h0
           dHduu  = siu2*h0
           dHduuu = siu2*invg*h0
           dHdfu  = tiu*dHdf
           dHdtu  = tiu*dHdt
           dHdffu = tiu*dHdff
           dHdfuu = siu2*dHdf
           dHdttu = tiu*dHdtt
           dHdtuu = siu2*dHdt
           dHdftu = tiu*dHdft
           
           dHda   = dHdf*dAda + dHdt*dtda + dHdu*duda
           dHdb   = dHdf*dAdb + dHdt*dtdb + dHdu*dudb
           dHdg   =             dHdt*dtdg
           dHdaa  = dHdff*dAda**2 + dHdtt*dtda**2 + dHduu*duda**2 &
                    + 2.d0*dHdft*dAda*dtda + 2.d0*dHdfu*dAda*duda + 2.d0*dHdtu*dtda*duda &
                    + dHdf*dAdaa + dHdt*dtdaa + dHdu*dudaa
           dHdab  = dHdff*dAda*dAdb + dHdtt*dtda*dtdb + dHduu*duda*dudb &
                    + dHdft*( dAda*dtdb + dAdb*dtda ) + dHdfu*( dAda*dudb + dAdb*duda ) + dHdtu*( dtda*dudb + dtdb*duda ) &
                    + dHdf*dAdab + dHdt*dtdab + dHdu*dudab
           dHdbb  = dHdff*dAdb**2 + dHdtt*dtdb**2 + dHduu*dudb**2 &
                    + 2.d0*dHdft*dAdb*dtdb + 2.d0*dHdfu*dAdb*dudb + 2.d0*dHdtu*dtdb*dudb &
                    + dHdf*dAdbb + dHdt*dtdbb + dHdu*dudbb
           dHdgg  = dHdtt*dtdg**2
           dHdag  = dHdft*dAda*dtdg + dHdtt*dtda*dtdg + dHdtu*duda*dtdg + dHdt*dtdag
           dHdbg  = dHdft*dAdb*dtdg + dHdtt*dtdb*dtdg + dHdtu*dudb*dtdg + dHdt*dtdbg
           dHdaaa = dHdfff*dAda**3 + 3.d0*dHdfft*dAda**2*dtda + 3.d0*dHdffu*dAda**2*duda &
                    + 3.d0*dHdftt*dAda*dtda**2 + dHdttt*dtda**3 + 3.d0*dHdttu*dtda**2*duda &
                    + 3.d0*dHdfuu*dAda*duda**2 + 3.d0*dHdtuu*dtda*duda**2 + dHduuu*duda**3 &
                    + 6.d0*dHdftu*dAda*dtda*duda &
                    + 3.d0*dHdff*dAda*dAdaa + 3.d0*dHdtt*dtda*dtdaa + 3.d0*dHduu*duda*dudaa &
                    + 3.d0*dHdft*( dAdaa*dtda + dAda*dtdaa ) &
                    + 3.d0*dHdtu*( dtdaa*duda + dtda*dudaa ) &
                    + 3.d0*dHdfu*( dAdaa*duda + dAda*dudaa ) &
                    + dHdf*dAdaaa + dHdt*dtdaaa + dHdu*dudaaa
           dHdaab = dHdfff*dAda**2*dAdb + dHdfft*dAda*( dAda*dtdb + 2.d0*dAdb*dtda ) &
                    + dHdffu*dAda*( dAda*dudb + 2.d0*dAdb*duda ) &
                    + dHdftt*dtda*( dAdb*dtda + 2.d0*dAda*dtdb ) &
                    + dHdttt*dtdb*dtda**2 + dHdttu*dtda*( dtda*dudb + 2.d0*dtdb*duda ) &
                    + dHdfuu*duda*( dAdb*duda + 2.d0*dAda*dudb ) &
                    + dHdtuu*duda*( dtdb*duda + 2.d0*dtda*dudb ) + dHduuu*duda**2*dudb &
                    + 2.d0*dHdftu*( dAdb*duda*dtda + dAda*dudb*dtda + dAda*duda*dtdb ) &
                    + dHdff*( 2.d0*dAda*dAdab + dAdb*dAdaa ) &
                    + dHdtt*( 2.d0*dtda*dtdab + dtdb*dtdaa ) &
                    + dHduu*( 2.d0*duda*dudab + dudb*dudaa ) &
                    + dHdft*( 2.d0*dAdab*dtda + 2.d0*dAda*dtdab + dAdaa*dtdb + dAdb*dtdaa ) &
                    + dHdtu*( 2.d0*dtdab*duda + 2.d0*dtda*dudab + dtdaa*dudb + dtdb*dudaa ) &
                    + dHdfu*( 2.d0*dAdab*duda + 2.d0*dAda*dudab + dAdaa*dudb + dAdb*dudaa ) &
                    + dHdf*dAdaab + dHdt*dtdaab + dHdu*dudaab
           dHdabb = dHdfff*dAda*dAdb**2 + dHdfft*dAdb*( dAdb*dtda + 2.d0*dAda*dtdb ) &
                    + dHdffu*dAdb*( dAdb*duda + 2.d0*dAda*dudb ) &
                    + dHdftt*dtdb*( dAda*dtdb + 2.d0*dAdb*dtda ) + dHdttt*dtda*dtdb**2 &
                    + dHdttu*dtdb*( dtdb*duda + 2.d0*dtda*dudb ) &
                    + dHdfuu*dudb*( dAda*dudb + 2.d0*dAdb*duda ) + dHdtuu*dudb*( dtda*dudb &
                    + 2.d0*dtdb*duda ) + dHduuu*duda*dudb**2 &
                    + 2.d0*dHdftu*( dAda*dtdb*dudb + dAdb*dtda*dudb + dAdb*dtdb*duda ) &
                    + dHdff*( 2.d0*dAdb*dAdab + dAda*dAdbb ) &
                    + dHdtt*( 2.d0*dtdb*dtdab + dtda*dtdbb ) &
                    + dHduu*( 2.d0*dudb*dudab + duda*dudbb ) &
                    + dHdft*( 2.d0*dAdab*dtdb + 2.d0*dAdb*dtdab + dAdbb*dtda + dAda*dtdbb ) &
                    + dHdtu*( 2.d0*dtdab*dudb + 2.d0*dtdb*dudab + dtdbb*duda + dtda*dudbb ) &
                    + dHdfu*( 2.d0*dAdab*dudb + 2.d0*dAdb*dudab + dAdbb*duda + dAda*dudbb ) &
                    + dHdf*dAdabb + dHdt*dtdabb + dHdu*dudabb
           dHdbbb = dHdfff*dAdb**3 + 3.d0*dHdfft*dAdb**2*dtdb + 3.d0*dHdffu*dAdb**2*dudb &
                    + 3.d0*dHdftt*dAdb*dtdb**2 + dHdttt*dtdb**3 + 3.d0*dHdttu*dtdb**2*dudb &
                    + 3.d0*dHdfuu*dAdb*dudb**2 + 3.d0*dHdtuu*dtdb*dudb**2 + dHduuu*dudb**3 &
                    + 6.d0*dHdftu*dAdb*dtdb*dudb &
                    + 3.d0*dHdff*dAdb*dAdbb + 3.d0*dHdtt*dtdb*dtdbb + 3.d0*dHduu*dudb*dudbb &
                    + 3.d0*dHdft*( dAdbb*dtdb + dAdb*dtdbb ) &
                    + 3.d0*dHdtu*( dtdbb*dudb + dtdb*dudbb ) &
                    + 3.d0*dHdfu*( dAdbb*dudb + dAdb*dudbb ) &
                    + dHdf*dAdbbb + dHdt*dtdbbb + dHdu*dudbbb
           dHdabg = dHdfft*dAda*dAdb*dtdg   + dHdftt*( dAdb*dtda + dAda*dtdb )*dtdg &
                    + dHdttt*dtda*dtdb*dtdg + dHdttu*( dtdb*duda + dtda*dudb )*dtdg &
                    + dHdtuu*duda*dudb*dtdg + dHdftu*( dAdb*duda + dAda*dudb )*dtdg &
                    + dHdft*( dAdab*dtdg + dAdb*dtdag + dAda*dtdbg ) &
                    + dHdtt*( dtdab*dtdg + dtdb*dtdag + dtda*dtdbg ) &
                    + dHdtu*( dudab*dtdg + dudb*dtdag + duda*dtdbg ) &
                    + dHdt*dtdabg
           dHdaag = dHdfft*dAda**2*dtdg + dHdttt*dtda**2*dtdg + dHdtuu*duda**2*dtdg &
                    + 2.d0*dHdttu*dtda*duda*dtdg + 2.d0*dHdftt*dAda*dtda*dtdg + 2.d0*dHdftu*dAda*duda*dtdg &
                    + dHdtt*( dtdaa*dtdg + 2.d0*dtda*dtdag ) + dHdft*( dAdaa*dtdg + 2.d0*dAda*dtdag ) &
                    + dHdtu*( dudaa*dtdg + 2.d0*duda*dtdag ) + dHdt*dtdaag         
           dHdbbg = dHdfft*dAdb**2*dtdg + dHdttt*dtdb**2*dtdg + dHdtuu*dudb**2*dtdg &
                    + 2.d0*dHdttu*dtdb*dudb*dtdg + 2.d0*dHdftt*dAdb*dtdb*dtdg + 2.d0*dHdftu*dAdb*dudb*dtdg &
                    + dHdtt*( dtdbb*dtdg + 2.d0*dtdb*dtdbg ) + dHdft*( dAdbb*dtdg + 2.d0*dAdb*dtdbg ) &
                    + dHdtu*( dudbb*dtdg + 2.d0*dudb*dtdbg ) + dHdt*dtdbbg       
           dHdagg = dHdftt*dAda*dtdg**2 + dHdttt*dtda*dtdg**2 + dHdttu*duda*dtdg**2 + 2.d0*dHdtt*dtdg*dtdag
           dHdbgg = dHdftt*dAdb*dtdg**2 + dHdttt*dtdb*dtdg**2 + dHdttu*dudb*dtdg**2 + 2.d0*dHdtt*dtdg*dtdbg 
           dHdggg = dHdttt*dtdg**3   
           
! ****** Derivative of h1 *****
           
           dq6   = c2+2.d0*c3*rs
           ddq6  = 2.d0*c3
           dq7   = c4+2.d0*c5*rs+3.d0*c6*rs2
           ddq7  = 2.d0*c5+6.d0*c6*rs
           dddq7 = 6.d0*c6
           
           iq7=1.d0/q7
           iq72=iq7*iq7
           iq73=iq72*iq7
           iq74=iq73*iq7
           
           ccrs2 =  ddq6*iq7-(2.d0*dq6*dq7+q6*ddq7)*iq72+2.d0*q6*dq7**2*iq73
           ccrs3 =  -iq72*(3.d0*(ddq6*dq7+dq6*ddq7)+q6*dddq7) + &
                                    6.d0*dq7*(dq6*dq7+q6*ddq7)*iq73 - &
                                    6.d0*q6*dq7**3*iq74
           r1rs2=r1rs*r1rs
           r1rs3=r1rs2*r1rs
           r1rs4=r1rs3*r1rs
           
           if(dabs(grad_trho(i)) > 1.0d-9) then
               ivt=1.d0/t
               ivt2=ivt**2
           else
               ivt=0.d0
               ivt2=0.d0
           end if
           
!           ivtmr1t=ivt-r1*t
           ivu=1.d0/g
           ivu2=1.d0/g/g
           r1t2=r1*t2
           r1t22= r1t2**2
           r1t23= r1t2*r1t22
           
           dH1dr   =  h1rs
           dH1du   =  h1*ivu*(3.d0-4.d0*r1t2)
           dH1dt   =  h1t
           dH1drr  =  xnu*ccrs2*t2*g3*r3-2.d0*h1rs*r1rs*t2-h1*r1rs2*t4
           dH1duu  =  2.d0*h1*ivu2*(3.d0 - 18.d0*r1t2 + 8.d0*r1t22) 
           dH1dtt  =  2.d0*h1*ivt2*(1.d0 -  5.d0*r1t2 + 2.d0*r1t22)
           dH1dru  =  ivu*(dH1dr*(3.d0 - 4.d0*r1t2)-4.d0*h1*r1rs*t2)
           dH1dtu  =  2.d0*ivu*ivt*h1*(3.d0-11.d0*r1t2+4.d0*r1t22)
           dH1drt  =  2.d0*(ivt*dH1dr*(1.d0-r1t2)-h1*r1rs*t)
           dH1drrr =  xnu*ccrs3*t2*g3*r3-3.d0*t2*r1rs*(dH1drr+h1rs*r1rs*t2)-h1*r1rs3*t6
           dH1drru =  ivu*(dH1drr*(3.d0-4.d0*r1t2)-8.d0*t2*dH1dr*r1rs)
           dH1druu =  2.d0*ivu2*(dH1dr*(3.d0-18.d0*r1t2+8.d0*r1t22)+2.d0*h1*t2*r1rs*(8.d0*r1t2-9.d0))
           dH1duuu =  2.d0*h1*ivu*ivu2*(3.d0-102.d0*r1t2+144.d0*r1t22-32.d0*r1t23)
           dH1dtuu =  4.d0*h1*ivu2*ivt*(3.d0-39.d0*r1t2+42.d0*r1t22-8.d0*r1t23)
           dH1dttu =  0.5d0*g*ivt*dH1dtuu
           dH1dttt =  4.d0*h1*r1*ivt*(-6.d0+9.d0*r1t2-2.d0*r1t22)
           dH1drtt =  2.d0*dH1dr*ivt2*(1.d0-5.d0*r1t2+2.d0*r1t22)+2.d0*h1*r1rs*(4.d0*r1t2-5.d0)
           dH1drrt =  2.d0*ivt*dH1drr*(1.d0-r1t2)-4.d0*dH1dr*r1rs*t
           dH1drtu =  2.d0*ivu*( ivt*dH1dr*(3.d0-11.d0*r1t2+4.d0*r1t22)+h1*t*r1rs*(-11.d0+8.d0*r1t2))
           
           dH1da   = dH1dr*rsdn + dH1dt*dtda + dH1du*duda
           dH1db   = dH1dr*rsdn + dH1dt*dtdb + dH1du*dudb
           dH1dg   =              dH1dt*dtdg
           dH1daa  = dH1drr*rsdn**2 + dH1dtt*dtda**2 + dH1duu*duda**2 &
                    + 2.d0*dH1drt*rsdn*dtda + 2.d0*dH1dru*rsdn*duda + 2.d0*dH1dtu*dtda*duda &
                    + dH1dr*rsdnn + dH1dt*dtdaa + dH1du*dudaa
           dH1dab  = dH1drr*rsdn*rsdn + dH1dtt*dtda*dtdb + dH1duu*duda*dudb &
                    + dH1drt*rsdn*(dtdb + dtda) + dH1dru*rsdn*( dudb + duda ) + dH1dtu*( dtda*dudb + dtdb*duda ) &
                    + dH1dr*rsdnn + dH1dt*dtdab + dH1du*dudab
           dH1dbb  = dH1drr*rsdn**2 + dH1dtt*dtdb**2 + dH1duu*dudb**2 &
                    + 2.d0*dH1drt*rsdn*dtdb + 2.d0*dH1dru*rsdn*dudb + 2.d0*dH1dtu*dtdb*dudb &
                    + dH1dr*rsdnn + dH1dt*dtdbb + dH1du*dudbb
           dH1dgg  = dH1dtt*dtdg**2
           dH1dag  = dH1drt*rsdn*dtdg + dH1dtt*dtda*dtdg + dH1dtu*duda*dtdg + dH1dt*dtdag
           dH1dbg  = dH1drt*rsdn*dtdg + dH1dtt*dtdb*dtdg + dH1dtu*dudb*dtdg + dH1dt*dtdbg
           dH1daaa = dH1drrr*rsdn3 + 3.d0*dH1drrt*rsdn2*dtda + 3.d0*dH1drru*rsdn2*duda &
                    + 3.d0*dH1drtt*rsdn*dtda**2 + dH1dttt*dtda**3 + 3.d0*dH1dttu*dtda**2*duda &
                    + 3.d0*dH1druu*rsdn*duda**2 + 3.d0*dH1dtuu*dtda*duda**2 + dH1duuu*duda**3 &
                    + 6.d0*dH1drtu*rsdn*dtda*duda &
                    + 3.d0*dH1drr*rsdn*rsdnn + 3.d0*dH1dtt*dtda*dtdaa + 3.d0*dH1duu*duda*dudaa &
                    + 3.d0*dH1drt*( rsdnn*dtda + rsdn*dtdaa ) &
                    + 3.d0*dH1dtu*( dtdaa*duda + dtda*dudaa ) &
                    + 3.d0*dH1dru*( rsdnn*duda + rsdn*dudaa ) &
                    + dH1dr*rsdnnn + dH1dt*dtdaaa + dH1du*dudaaa
           dH1daab = dH1drrr*rsdn3 + dH1drrt*rsdn2*( dtdb + 2.d0*dtda ) + dH1drru*rsdn2*( dudb + 2.d0*duda ) &
                    + dH1drtt*dtda*rsdn*( dtda + 2.d0*dtdb ) + dH1dttt*dtdb*dtda**2 &
                    + dH1dttu*dtda*( dtda*dudb + 2.d0*dtdb*duda ) &
                    + dH1druu*duda*rsdn*( duda + 2.d0*dudb ) &
                    + dH1dtuu*duda*( dtdb*duda + 2.d0*dtda*dudb ) + dH1duuu*duda**2*dudb &
                    + 2.d0*dH1drtu*rsdn*( duda*dtda + dudb*dtda + duda*dtdb ) &
                    + dH1drr*3.d0*rsdn*rsdnn &
                    + dH1dtt*( 2.d0*dtda*dtdab + dtdb*dtdaa ) &
                    + dH1duu*( 2.d0*duda*dudab + dudb*dudaa ) &
                    + dH1drt*( 2.d0*rsdnn*dtda + 2.d0*rsdn*dtdab + rsdnn*dtdb + rsdn*dtdaa ) &
                    + dH1dtu*( 2.d0*dtdab*duda + 2.d0*dtda*dudab + dtdaa*dudb + dtdb*dudaa ) &
                    + dH1dru*( 2.d0*rsdnn*duda + 2.d0*rsdn*dudab + rsdnn*dudb + rsdn*dudaa ) &
                    + dH1dr*rsdnnn + dH1dt*dtdaab + dH1du*dudaab
           dH1dabb = dH1drrr*rsdn3 + dH1drrt*rsdn2*( dtda + 2.d0*dtdb ) + dH1drru*rsdn2*( duda + 2.d0*dudb ) &
                    + dH1drtt*dtdb*rsdn*( dtdb + 2.d0*dtda ) &
                    + dH1dttt*dtda*dtdb**2 + dH1dttu*dtdb*( dtdb*duda + 2.d0*dtda*dudb ) &
                    + dH1druu*dudb*rsdn*( dudb + 2.d0*duda ) &
                    + dH1dtuu*dudb*( dtda*dudb + 2.d0*dtdb*duda ) + dH1duuu*duda*dudb**2 &
                    + 2.d0*dH1drtu*rsdn*( dtdb*dudb + dtda*dudb + dtdb*duda ) &
                    + dH1drr*3.d0*rsdn*rsdnn &
                    + dH1dtt*( 2.d0*dtdb*dtdab + dtda*dtdbb ) &
                    + dH1duu*( 2.d0*dudb*dudab + duda*dudbb ) &
                    + dH1drt*( 2.d0*rsdnn*dtdb + 2.d0*rsdn*dtdab + rsdnn*dtda + rsdn*dtdbb ) &
                    + dH1dtu*( 2.d0*dtdab*dudb + 2.d0*dtdb*dudab + dtdbb*duda + dtda*dudbb ) &
                    + dH1dru*( 2.d0*rsdnn*dudb + 2.d0*rsdn*dudab + rsdnn*duda + rsdn*dudbb ) &
                    + dH1dr*rsdnnn + dH1dt*dtdabb + dH1du*dudabb
           dH1dbbb = dH1drrr*rsdn3 + 3.d0*dH1drrt*rsdn2*dtdb + 3.d0*dH1drru*rsdn2*dudb &
                    + 3.d0*dH1drtt*rsdn*dtdb**2 + dH1dttt*dtdb**3 + 3.d0*dH1dttu*dtdb**2*dudb &
                    + 3.d0*dH1druu*rsdn*dudb**2 + 3.d0*dH1dtuu*dtdb*dudb**2 + dH1duuu*dudb**3 &
                    + 6.d0*dH1drtu*rsdn*dtdb*dudb &
                    + 3.d0*dH1drr*rsdn*rsdnn + 3.d0*dH1dtt*dtdb*dtdbb + 3.d0*dH1duu*dudb*dudbb &
                    + 3.d0*dH1drt*( rsdnn*dtdb + rsdn*dtdbb ) &
                    + 3.d0*dH1dtu*( dtdbb*dudb + dtdb*dudbb ) &
                    + 3.d0*dH1dru*( rsdnn*dudb + rsdn*dudbb ) &
                    + dH1dr*rsdnnn + dH1dt*dtdbbb + dH1du*dudbbb
           dH1dabg = dH1drrt*rsdn2*dtdg   + dH1drtt*rsdn*( dtda + dtdb )*dtdg &
                    + dH1dttt*dtda*dtdb*dtdg + dH1dttu*( dtdb*duda + dtda*dudb )*dtdg &
                    + dH1dtuu*duda*dudb*dtdg + dH1drtu*rsdn*( duda + dudb )*dtdg &
                    + dH1drt*( rsdnn*dtdg + rsdn*dtdag + rsdn*dtdbg ) &
                    + dH1dtt*( dtdab*dtdg + dtdb*dtdag + dtda*dtdbg ) &
                    + dH1dtu*( dudab*dtdg + dudb*dtdag + duda*dtdbg ) &
                    + dH1dt*dtdabg
           dH1daag = dH1drrt*rsdn2*dtdg + dH1dttt*dtda**2*dtdg + dH1dtuu*duda**2*dtdg &
                    + 2.d0*dH1dttu*dtda*duda*dtdg + 2.d0*dH1drtt*rsdn*dtda*dtdg + 2.d0*dH1drtu*rsdn*duda*dtdg &
                    + dH1dtt*( dtdaa*dtdg + 2.d0*dtda*dtdag ) + dH1drt*( rsdnn*dtdg + 2.d0*rsdn*dtdag ) &
                    + dH1dtu*( dudaa*dtdg + 2.d0*duda*dtdag ) + dH1dt*dtdaag  
           dH1dbbg = dH1drrt*rsdn2*dtdg + dH1dttt*dtdb**2*dtdg + dH1dtuu*dudb**2*dtdg &
                    + 2.d0*dH1dttu*dtdb*dudb*dtdg + 2.d0*dH1drtt*rsdn*dtdb*dtdg + 2.d0*dH1drtu*rsdn*dudb*dtdg &
                    + dH1dtt*( dtdbb*dtdg + 2.d0*dtdb*dtdbg ) + dH1drt*( rsdnn*dtdg + 2.d0*rsdn*dtdbg ) &
                    + dH1dtu*( dudbb*dtdg + 2.d0*dudb*dtdbg ) + dH1dt*dtdbbg       
           dH1dagg = dH1drtt*rsdn*dtdg**2 + dH1dttt*dtda*dtdg**2 + dH1dttu*duda*dtdg**2 + 2.d0*dH1dtt*dtdg*dtdag
           dH1dbgg = dH1drtt*rsdn*dtdg**2 + dH1dttt*dtdb*dtdg**2 + dH1dttu*dudb*dtdg**2 + 2.d0*dH1dtt*dtdg*dtdbg 
           dH1dggg = dH1dttt*dtdg**3   
                    
!print '(i2,20d19.10)', i,h0,dHda,dHdb,dHdg,dHdaa,dHdbb,dHdgg,dHdag,dHdbg,dHdab,dHdaaa,dHdbbb,dHdggg,dHdaab,dHdabb,dHdaag,dHdagg,dHdbbg,dHdbgg,dHdabg
!print '(i2,20d19.10)', i,h1,dH1da,dH1db,dH1dg,dH1daa,dH1dbb,dH1dgg,dH1dag,dH1dbg,dH1dab,dH1daaa,dH1dbbb,dH1dggg,dH1daab,dH1dabb,dH1daag,dH1dagg,dH1dbbg,dH1dbgg,dH1dabg
    
           dFc_daa(i)   = 2.d0*euda + d*eudaa + 2.d0*(dHda + dH1da) + d*(dHdaa + dH1daa) 
           dFc_dbb(i)   = 2.d0*eudb + d*eudbb + 2.d0*(dHdb + dH1db) + d*(dHdbb + dH1dbb)
           dFc_dgg(i)   = d*(dHdgg + dH1dgg)
           dFc_dab(i)   = euda + eudb + d*eudab + dHda + dH1da + dHdb + dH1db + d*(dHdab + dH1dab) 
           dFc_dag(i)   = dHdg + dH1dg + d*(dHdag + dH1dag)
           dFc_dbg(i)   = dHdg + dH1dg + d*(dHdbg + dH1dbg)
           dFc_daaa(i)  = 3.d0*eudaa + d*eudaaa + 3.d0*(dHdaa + dH1daa) + d*(dHdaaa + dH1daaa)
           dFc_dbbb(i)  = 3.d0*eudbb + d*eudbbb + 3.d0*(dHdbb + dH1dbb) + d*(dHdbbb + dH1dbbb)
           dFc_dggg(i)  = d*(dHdggg + dH1dggg)
           dFc_daab(i)  = eudaa + 2.d0*eudab + d*eudaab + dHdaa + dH1daa + 2.d0*(dHdab + dH1dab) + d*(dHdaab + dH1daab)
           dFc_daag(i)  = 2.d0*(dHdag + dH1dag) + d*(dHdaag + dH1daag)
           dFc_dabb(i)  = eudbb + 2.d0*eudab + d*eudabb + dHdbb + dH1dbb + 2.d0*(dHdab + dH1dab) + d*(dHdabb + dH1dabb)
           dFc_dbbg(i)  = 2.d0*(dHdbg + dH1dbg) + d*(dHdbbg + dH1dbbg)
           dFc_dagg(i)  = dHdgg + dH1dgg + d*(dHdagg + dH1dagg)
           dFc_dbgg(i)  = dHdgg + dH1dgg + d*(dHdbgg + dH1dbgg)
           dFc_dabg(i)  = dHdag + dH1dag + dHdbg + dH1dbg + d*(dHdabg + dH1dabg)
                    
!print '(i2,21d19.10)', i,exc(i),d*(eu+h),dF_drho(i,1),dF_drho(i,2),dF_dgradrho(i)*grad_trho(i) &
!                      ,dFc_daa(i),dFc_dbb(i),dFc_dgg(i),dFc_dab(i),dFc_dag(i),dFc_dbg(i),dFc_daaa(i),dFc_dbbb(i) &
!                      ,dFc_dggg(i),dFc_daab(i),dFc_daag(i),dFc_dabb(i),dFc_dbbg(i),dFc_dagg(i),dFc_dbgg(i),dFc_dabg(i)
!print *,dHdag,dHdaag,d,2.d0*dHdag + d*dHdaag
!print *,2.d0*dHdbg + d*dHdbbg
!print *,dHdag + dHdbg + d*dHdabg
               
           
        end do
     end if
     
end subroutine cr_ggapw91_paw_drv2_3D
! ==============================================================================

subroutine cr_ggapw91_paw_drv3(nrc,num_ir,irs,nspin,chgrhr_l,grad_trho,exc,dF_drho,dF_dgradrho &
                                ,dFc_daa,dFc_dbb,dFc_dgg,dFc_dab,dFc_dag,dFc_dbg &
                                ,dFc_daaa,dFc_dbbb,dFc_dggg &
                                ,dFc_daab,dFc_daag,dFc_dabb,dFc_dbbg,dFc_dagg,dFc_dbgg,dFc_dabg)
!                           @(#)b_XC_Potential.F90 1.4 02/09/27 20:39:38
  use m_Const_Parameters,  only : PAI,DP
  implicit none

  integer,intent(in)        :: nrc,num_ir,irs(num_ir),nspin
  real(kind=DP),intent(in)  :: chgrhr_l(nrc,nspin)
  real(kind=DP),intent(in) :: grad_trho(nrc)
  real(kind=DP),intent(inout) :: exc(nrc)
  real(kind=DP),intent(out) :: dF_drho(nrc,nspin)
  real(kind=DP),intent(out) :: dF_dgradrho(nrc)
  real(kind=DP),intent(out) :: dFc_daa(nrc)
  real(kind=DP),intent(out) :: dFc_dbb(nrc)
  real(kind=DP),intent(out) :: dFc_dgg(nrc)
  real(kind=DP),intent(out) :: dFc_dab(nrc)
  real(kind=DP),intent(out) :: dFc_dag(nrc)
  real(kind=DP),intent(out) :: dFc_dbg(nrc)
  real(kind=DP),intent(out) :: dFc_daaa(nrc)
  real(kind=DP),intent(out) :: dFc_dbbb(nrc)
  real(kind=DP),intent(out) :: dFc_dggg(nrc)
  real(kind=DP),intent(out) :: dFc_daab(nrc)
  real(kind=DP),intent(out) :: dFc_daag(nrc)
  real(kind=DP),intent(out) :: dFc_dabb(nrc)
  real(kind=DP),intent(out) :: dFc_dbbg(nrc)
  real(kind=DP),intent(out) :: dFc_dagg(nrc)
  real(kind=DP),intent(out) :: dFc_dbgg(nrc)
  real(kind=DP),intent(out) :: dFc_dabg(nrc)

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
  real(kind=DP), parameter :: cx  = -0.001667212d0, alf =  0.09d0
  real(kind=DP), parameter :: c1  =  0.002568d0,    c2  =  0.023266d0
  real(kind=DP), parameter :: c3  =  7.389d-6,      c4  =  8.723d0
  real(kind=DP), parameter :: c5  =  0.472d0,       c6  =  7.389d-2
  real(kind=DP), parameter :: a4  = 100.d0
  
  real(kind=DP), parameter  :: thrd = 0.333333333333d0, sixth7 = 1.1666666666666666d0
  real(kind=DP), parameter  :: thrd7= 2.333333333333d0
  real(kind=DP), parameter  :: sixth13= 2.1666666666666666d0
  real(kind=DP), parameter  :: sixth19= 3.1666666666666666d0
  real(kind=DP), parameter  :: ninth= 0.111111111111111d0
  real(kind=DP), parameter  :: twntysvnth=  3.703703703703703d-2
  real(kind=DP), parameter  :: fftyfrth= 1.851851851851852d-2
  real(kind=DP), parameter  :: eityonth= 1.234567901234568d-2

  integer       :: is,i,nr
  real(kind=DP) :: facw,bet,delt,exc1,g,g3,g4,facpon,d,dd,rs,fk,sk,t,s,q0 &
       &         , rs12,rs32,rsp,q1,q2,eu,q3,eurs,ecd,pon &
       &         , b,b2,t2,s2,t4,t6,rs2,rs3,q4,q5,q6,q7,cc,r0,r1,coeff,r2,r3,a4ms,h0,h1,h,q8 &
       &         , h0t,h0b,h0rs,h1t,ccrs,r1rs,h1rs,ht,hrs &
       &         , ec1,ec1d,ec1dd,exc0,excd,excdd &
       &         , thrd2,thrd4,fzdd,zeta,onpzeta,onmzeta,fz,zetadxd &
       &         , q0p,q1p,q2p,q3p,q0q,q1q,q2q,q3q,fzd,onzt4,gzt,bzt,h0zt,h1zt,euzt &
       &         , zetadxda,zetadxdb,ecda,ecdb,ec1da,ec1db,excda,excdb
  real(kind=DP) :: ccrs2,ccrs3,dq6,ddq6,dq7,ddq7,dddq7
  real(kind=DP) :: iq7,iq72,iq73,iq74
       
  real(kind=DP) :: rsdn,rsdnn,rsdnnn
  real(kind=DP) :: q1drr,q1drrr
  real(kind=DP) :: eudrr,eudrrr,invq12pq1,q12p1
  real(kind=DP) :: eudn,eudnn,eudnnn,ecdd,ecddd
  real(kind=DP) :: dtdn,dtdg,dtdnn,dtdgg,dtdng,dtdnnn,dtdnng,dtdngg,dtdggg
  real(kind=DP) :: dAde,dAdee,dAdeee
  real(kind=DP) :: dAdn,dAdnn,dAdnnn
  real(kind=DP) :: fE,fEdt,fEdtt,fEda,fEdaa
  real(kind=DP) :: fF,fFdt,fFdtt,fFda,fFdaa
  real(kind=DP) :: fG,fGdt,fGdtt,fGda,fGdaa
  real(kind=DP) :: at2,a2t4,a3t6,t3,t5,t7,t8
  real(kind=DP) :: dHda,dHdt,dHdaa,dHdat,dHdtt,dHdaaa,dHdaat,dHdatt,dHdttt,betivE,betivE2
  real(kind=DP) :: dHdn,dHdg,dHdnn,dHdng,dHdgg,dHdnnn,dHdnng,dHdngg,dHdggg
  real(kind=DP) :: zetada,zetadb,zetadaa,zetadab,zetadbb,zetadaaa,zetadaab,zetadabb,zetadbbb
  real(kind=DP) :: id2,id3,id4
  real(kind=DP) :: e1dr,e1drr,e1drrr
  real(kind=DP) :: e1pdr,e1pdrr,e1pdrrr
  real(kind=DP) :: e1qdr,e1qdrr,e1qdrrr
  real(kind=DP) :: q1pdrr,q1pdrrr,invq1p2pq1p,q1p2p1
  real(kind=DP) :: q1qdrr,q1qdrrr,invq1q2pq1q,q1q2p1
  real(kind=DP) :: fzdd0,fzddd,opzthrdm2,omzthrdm2
  real(kind=DP) :: zeta2,zeta3,zeta4,eudrz,eudzz,eudrrz,eudrzz,eudzzz
  real(kind=DP) :: e1,e2,e3
  real(kind=DP) :: euda,eudb,eudaa,eudbb,eudab,eudaaa,eudaab,eudabb,eudbbb
  real(kind=DP) :: rsdn2,rsdn3,zetada2,zetada3,zetadb2,zetadb3
  real(kind=DP) :: gzz,gzzz,duda,dudb,dudaa,dudab,dudbb,dudaaa,dudaab,dudabb,dudbbb
  real(kind=DP) :: dtdu,dtduu,dtduuu,dtdnu,dtdnnu,dtdnuu,dtdug,dtduug,dtdugg,dtdnug,invg
  real(kind=DP) :: dtda,dtdb,dtdaa,dtdab,dtdag,dtdbg,dtdbb
  real(kind=DP) :: dtdaaa,dtdaab,dtdabb,dtdbbb,dtdaag,dtdagg,dtdbbg,dtdbgg,dtdabg
  real(kind=DP) :: invg3,m3u
  real(kind=DP) :: dAdu,dAduu,dAduuu,dAdeu,dAdeeu,dAdeuu
  real(kind=DP) :: dAda,dAdb,dAdaa,dAdab,dAdbb,dAdaaa,dAdaab,dAdabb,dAdbbb
  real(kind=DP) :: dHdf,dHdff,dHdft,dHdfff,dHdfft,dHdftt
  real(kind=DP) :: dHdu,dHduu,dHduuu,dHdfu,dHdtu,dHdffu,dHdfuu,dHdttu,dHdtuu,dHdftu
  real(kind=DP) :: tiu,siu2
  real(kind=DP) :: dHdb,dHdbb,dHdag,dHdbg,dHdab,dHdbbb,dHdaab,dHdabb,dHdbbg,dHdbgg,dHdaag,dHdagg,dHdabg
  
  real(kind=DP) :: dH1dr,dH1dt,dH1drr,dH1drt,dH1dtt,dH1drrr,dH1drrt,dH1drtt,dH1dttt,dH1drtu
  real(kind=DP) :: dH1du,dH1duu,dH1dru,dH1dtu,dH1duuu,dH1drru,dH1druu,dH1dtuu,dH1dttu
  real(kind=DP) :: r1rs2,r1rs3,r1rs4,ivt,ivt2,ivtmr1t,ivu,ivu2,r1t2,r1t22,r1t23
  real(kind=DP) :: dH1dn,dH1dg,dH1dnn,dH1dng,dH1dgg,dH1dnnn,dH1dnng,dH1dngg,dH1dggg
  real(kind=DP) :: dH1da,dH1db,dH1daa,dH1dab,dH1dbb,dH1dbg,dH1dag
  real(kind=DP) :: dH1daaa,dH1daab,dH1dabb,dH1dbbb,dH1dbbg,dH1dbgg,dH1dagg,dH1daag,dH1dabg

  facw = nspin
  bet = xnu*cc0
  delt = 2.d0*alf/bet
  
     if ( nspin == 1 ) then
        g   = 1.d0
        g3  = g**3
        g4  = g3*g
        facpon = -delt/(g3*bet)

!        do i = 1,nrc,dnr ! MPI
        do nr=1,num_ir
           i=irs(nr)
           d  = facw*chgrhr_l(i, 1)
           if(d < density_minimum) cycle
! gradient of charge density 95/12/3 H.S.
           dd = facw* grad_trho(i)
           rs = (0.75/(PAI*d))**thrd
           fk = (3*PAI*PAI*d)**thrd
           sk = dsqrt(4*fk/PAI)
           if(d > density_minimum2) then
              t = dd/(d*sk*2)
              s = dd/(d*fk*2)
           else
              t = 0.d0
              s = 0.d0
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
           s2  = s*s
           t4  = t2*t2
           t6  = t4*t2
           rs2 = rs*rs
           rs3 = rs2*rs
           q4  = 1 + b*t2
           q5  = 1 + b*t2+b2*t4
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
           h0t = 2*bet*t*(1+2*b*t2)/q8 * g3
           h0b = -bet*t6*(2*b+b2*t2)/q8 * g3
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

           exc0 = eu*d/facw + ec1
           excd = ecd + ec1d
           dF_drho(i, 1) = dF_drho(i, 1) + excd
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
           if(dabs(grad_trho(i)) > 1.0d-9) then
!              dF_dgradrho(i) = excdd / grad_trho(i)
              dF_dgradrho(i) = excdd
           else
              dF_dgradrho(i) = 0.d0
           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d/facw
!           exc1 = exc1 + exc0*wos(i)
           exc(i)=exc(i)+exc0
        
! ***** Calc heigher derivatives. *****
           
           rsdn=-thrd*rs/d
           rsdnn=-4.d0*thrd*rsdn/d
           rsdnnn=-thrd7*rsdnn/d
!print '(i2,5e19.10)',i,d,rs,rsdn,rsdnn,rsdnnn
!q3 = a*(b1/rs12+2*b2n+3*b3*rs12+2*b4*p1*rsp)   
           q1drr =0.5d0*a*( -b1/rs32+3.d0*b3/rs12+8.d0*b4 )
           q1drrr=0.75d0*a*( b1/rs - b3 )/rs32
!eurs = -2*a*a1*q2 - q0*q3/(q1**2 + q1)
           invq12pq1=1.d0/(q1**2 + q1)
           q12p1=2.d0*q1+1.d0
           eudrr=4.d0*a*a1*q3*invq12pq1-q0*( q1drr-q12p1*q3**2*invq12pq1 )*invq12pq1
           eudrrr=6.d0*a*a1*( q1drr - q12p1*q3**2*invq12pq1 ) * invq12pq1 &
                    -q0*( q1drrr -(2.d0*q3**3+3.d0*q12p1*q3*q1drr -  &
                    2.d0*q12p1**2*q3**3* invq12pq1)*invq12pq1 ) * invq12pq1
           eudn=eurs*rsdn
           eudnn=eudrr*rsdn**2+eurs*rsdnn
           eudnnn=eudrrr*rsdn**3+3.d0*eudrr*rsdn*rsdnn+eurs*rsdnnn
!print '(i2,5e19.12)',i,d,eu,eudn,eudnn,eudnnn
           ecdd=2.d0*eudn+d*eudnn
           ecddd=3.d0*eudnn+d*eudnnn
!print '(i2,5e19.12)',i,d,eu*d,ecdd,ecddd
           
           dtdn  = -sixth7*t/d
           dtdnn = -sixth13*dtdn/d
           dtdnnn= -sixth19*dtdnn/d
!print '(i2,5e19.12)',i,d,t,dtdn,dtdnn,dtdnnn
           if(dabs(grad_trho(i)) > 1.0d-9) then
               dtdg  = t/dd
               dtdng = -sixth7*dtdg/d
               dtdnng= -sixth13*dtdng/d
           else
               dtdg  = 0.d0
               dtdng = 0.d0
               dtdnng= 0.d0
           endif
           dtdngg = 0.d0
           dtdggg = 0.d0
           
           dAde   = b * ( b + delt ) / bet
           dAdee  = 2.d0 * dAde * ( b + 0.5d0 * delt ) / bet
           dAdeee = 2.d0 * ( dAde**2 + ( b + 0.5d0 * delt )*dAdee ) /bet
           
           dAdn   = dAde  *eudn
           dAdnn  = dAdee *eudn**2  + dAde*eudnn
           dAdnnn = dAdeee*eudn**3  + 3.d0*dAdee*eudn*eudnn + dAde*eudnnn   
!print '(i2,5e19.10)' ,i,d,b,dAdnn
           
           at2  = b * t2
           a2t4 = at2**2
           a3t6 = at2*a2t4
           t3   = t2*t
           t5   = t4*t
           t7   = t6*t
           t8   = t7*t
           fE    = q5**2 + delt*t2*q4*q5
           fEdt  = 4.d0*b*t*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt *2.d0*t*(1.d0 + 4.d0*at2 + 6.d0*a2t4 + 4.d0*a3t6 )
           fEdtt = 4.d0*b*( 1.d0 + 9.d0*at2 + 15.d0*a2t4 + 14.d0*a3t6 ) + &
                    delt*( 2.d0 + 24.d0*at2 + 60.d0*a2t4 + 56.d0*a3t6 )
           fEda  = 2.d0*t2*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt*t4*( 2.d0 + 4.d0*at2 + 3.d0*a2t4 )
           fEdaa = 6.d0*t4*( 2.d0*q5 - 1.d0) + &
                    delt*2.d0*t6*( 2.d0 + 3.d0*at2 )
           fF    = t*( 1.d0 + 2.d0*at2 )         
           fFdt  = 1.d0 + 6.d0*at2
           fFdtt = 12.d0*b*t
           fFda  = 2.d0*t3
           fFdaa = 0.d0
           fG    = 2.d0*b*t6 + b2*t8
           fGdt  = 12.d0*b*t5 + 8.d0*b2*t7
           fGdtt = 60.d0*b*t4 + 56.d0*b2*t6
           fGda  = 2.d0*t6 + 2.d0*b*t8
           fGdaa = 2.d0*t8
           
           betivE=bet/fE
           betivE2=betivE/fE
           dHda   =     -fG*betivE
           dHdt   = 2.d0*fF*betivE
           dHdaa  =     -( fGda*fE-fG*fEda )*betivE2
           dHdat  =     -( fGdt*fE-fG*fEdt )*betivE2
           dHdtt  = 2.d0*( fFdt*fE-fF*fEdt )*betivE2
           dHdaaa =     -( fGdaa*fE-fG*fEdaa-2.d0*( fGda*fE-fG*fEda )*fEda/fE )*betivE2
           dHdaat = 2.d0*( fFdaa*fE-fF*fEdaa-2.d0*( fFda*fE-fF*fEda )*fEda/fE )*betivE2
           dHdatt =     -( fGdtt*fE-fG*fEdtt-2.d0*( fGdt*fE-fG*fEdt )*fEdt/fE )*betivE2
           dHdttt = 2.d0*( fFdtt*fE-fF*fEdtt-2.d0*( fFdt*fE-fF*fEdt )*fEdt/fE )*betivE2
           
           dHdn     = dHda  *dAdn       + dHdt*dtdn
           dHdg     = dHdt  *dtdg
           dHdnn    = dHdaa *dAdn**2    + 2.d0*dHdat*dAdn*dtdn  + dHdtt*dtdn**2   + dHda*dAdnn  + dHdt*dtdnn
           dHdng    = dHdat *dAdn*dtdg  + dHdtt*dtdn*dtdg       + dHdt*dtdng
           dHdgg    = dHdtt*dtdg**2
           dHdnnn   = dHdaaa*dAdn**3    + dHdttt*dtdn**3 &
                        + 3.d0*dHdaat*dAdn**2*dtdn  + 3.d0*dHdatt*dAdn*dtdn**2 &
                        + 3.d0*dHdaa*dAdn*dAdnn + 3.d0*dHdtt*dtdn*dtdnn + 3.d0*dHdat*(dAdnn*dtdn + dAdn*dtdnn) &
                        + dHda*dAdnnn   + dHdt*dtdnnn
           dHdnng   = dHdaat*dAdn**2*dtdg   + 2.d0*dHdatt*dAdn*dtdn*dtdg    + dHdttt*dtdn**2*dtdg &
                        + dHdat*(dAdnn*dtdg + 2.d0*dAdn*dtdng)  + dHdtt*(dtdnn*dtdg + 2.d0*dtdn*dtdng) &
                        + dHdt*dtdnng
           dHdngg   = dHdatt*dAdn*dtdg**2   + dHdttt*dtdn*dtdg**2   + 2.d0*dHdtt*dtdg*dtdng
           dHdggg   = dHdttt*dtdg**3
           
           dq6   = c2+2.d0*c3*rs
           ddq6  = 2.d0*c3
           dq7   = c4+2.d0*c5*rs+3.d0*c6*rs2
           ddq7  = 2.d0*c5+6.d0*c6*rs
           dddq7 = 6.d0*c6
           
           iq7=1.d0/q7
           iq72=iq7*iq7
           iq73=iq72*iq7
           iq74=iq73*iq7
           
           ccrs2 =  ddq6*iq7-(2.d0*dq6*dq7+q6*ddq7)*iq72+2.d0*q6*dq7**2*iq73
           ccrs3 =  -iq72*(3.d0*(ddq6*dq7+dq6*ddq7)+q6*dddq7) + &
                                    6.d0*dq7*(dq6*dq7+q6*ddq7)*iq73 - &
                                    6.d0*q6*dq7**3*iq74
           r1rs2=r1rs*r1rs
           r1rs3=r1rs2*r1rs
           r1rs4=r1rs3*r1rs
           
           if(dabs(grad_trho(i)) > 1.0d-9) then
               ivt=1.d0/t
               ivt2=ivt**2
               ivtmr1t=ivt-r1*t
           else
               ivt=0.d0
               ivt2=0.d0
               ivtmr1t=0.d0
           end if
           
           dH1dr=h1rs
           dH1dt=h1t                     
           dH1drr  = xnu*ccrs2*t2*r3-2.d0*h1rs*r1rs*t2-h1*r1rs2*t4
           dH1drt  = 2.d0*h1rs*ivtmr1t-2.d0*h1*t*r1rs
           dH1dtt  = 2.d0*h1t*ivtmr1t-2.d0*h1*(ivt2+r1)
           dH1drrr = xnu*ccrs3*t2*r3-3.d0*t2*r1rs*(dH1drr+h1rs*r1rs*t2)-h1*r1rs3*t6
           dH1drrt = 2.d0*dH1drr*ivtmr1t-4.d0*h1rs*r1rs*t
           dH1drtt = 2.d0*dH1drt*ivtmr1t-2.d0*h1t*r1rs*t-2.d0*ivt*h1rs*(ivt+r1*t)-2.d0*h1*r1rs
           dH1dttt = 2.d0*dH1dtt*ivtmr1t-4.d0*h1t*(ivt2+r1)+4.d0*h1*ivt2*ivt
           
           dH1dn     = dH1dr  *rsdn       + dH1dt*dtdn
           dH1dg     = dH1dt  *dtdg
           dH1dnn    = dH1drr *rsdn**2    + 2.d0*dH1drt*rsdn*dtdn  + dH1dtt*dtdn**2   + dH1dr*rsdnn  + dH1dt*dtdnn
           dH1dng    = dH1drt *rsdn*dtdg  + dH1dtt*dtdn*dtdg       + dH1dt*dtdng
           dH1dgg    = dH1dtt*dtdg**2
           dH1dnnn   = dH1drrr*rsdn**3    + dH1dttt*dtdn**3 &
                        + 3.d0*dH1drrt*rsdn**2*dtdn  + 3.d0*dH1drtt*rsdn*dtdn**2 &
                        + 3.d0*dH1drr*rsdn*rsdnn + 3.d0*dH1dtt*dtdn*dtdnn + 3.d0*dH1drt*(rsdnn*dtdn + rsdn*dtdnn) &
                        + dH1dr*rsdnnn   + dH1dt*dtdnnn
           dH1dnng   = dH1drrt*rsdn**2*dtdg   + 2.d0*dH1drtt*rsdn*dtdn*dtdg    + dH1dttt*dtdn**2*dtdg &
                        + dH1drt*(rsdnn*dtdg + 2.d0*rsdn*dtdng)  + dH1dtt*(dtdnn*dtdg + 2.d0*dtdn*dtdng) &
                        + dH1dt*dtdnng
           dH1dngg   = dH1drtt*rsdn*dtdg**2   + dH1dttt*dtdn*dtdg**2   + 2.d0*dH1dtt*dtdg*dtdng
           dH1dggg   = dH1dttt*dtdg**3
           
!dF_drho(i, 1)=h+d*dHdn+ecd
           dFc_daa(i)   = 2.d0*(dHdn + dH1dn) + d*(dHdnn + dH1dnn)       + ecdd
           dFc_dag(i)   = dHdg + dH1dg + d*(dHdng + dH1dng)
           dFc_dgg(i)   = d*(dHdgg + dH1dgg)
           dFc_daaa(i)  = 3.d0*(dHdnn + dH1dnn)    + d*(dHdnnn + dH1dnnn)  + ecddd
           dFc_daag(i)  = 2.d0*(dHdng + dH1dng)    + d*(dHdnng + dH1dnng)
           dFc_dagg(i)  = dHdgg + dH1dgg        + d*(dHdngg + dH1dngg)
           dFc_dggg(i)  = d*(dHdggg + dH1dggg)  
        end do
     else if ( nspin ==  2 ) then
        thrd2 = thrd * 2
        thrd4 = thrd * 4
        fzdd = 8/(9*(2**thrd4 - 2 ))

!        do i = 1,nrc,dnr
        do nr=1,num_ir
           i=irs(nr)
           d = chgrhr_l(i, 1) + chgrhr_l(i, nspin)
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
           g4  = g3*g
           facpon = -delt/(g3*bet)
           if(d < density_minimum) cycle
           dd = grad_trho(i)
           rs = (0.75d0/(PAI*d))**thrd
           fk = (3*PAI*PAI*d)**thrd
           sk = dsqrt(4*fk/PAI)

           if(d > density_minimum2) then
           ! ****** aas mkatsu *****
!              t = dd/(d*sk*2)
              t = dd/(d*sk*2)/g         !    mkatsu aas
              s = dd/(d*fk*2)
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

!           zetadxd = - 2 * chgrhr_l(i, 3-is) * (-1.d0)**is / d
           zetadxda = - 2 * chgrhr_l(i, 2) * (-1.d0) / d
           zetadxdb = - 2 * chgrhr_l(i, 1) / d
!           ecd = eu-thrd*rs*eurs + euzt * zetadxd
           ecda = eu-thrd*rs*eurs + euzt * zetadxda 
           ecdb = eu-thrd*rs*eurs + euzt * zetadxdb
           if(onpzeta < zeta_minimum2 .or. onmzeta < zeta_minimum2) then
              ec1  = 0.d0;    ec1d = 0.d0;     ec1dd= 0.d0
              ec1da = 0.d0;     ec1db = 0.d0
              go to 10
           end if
           pon = facpon*eu
           b   = delt/(dexp(pon)-1)
           b2  = b*b
           t2  = t*t
           s2  = s*s
           t4  = t2*t2
           t6  = t4*t2
           rs2 = rs*rs
           rs3 = rs2*rs
           q4  = 1+b*t2
           q5  = 1+b*t2+b2*t4
           q6  = c1+c2*rs+c3*rs2
           q7  = 1.d0+c4*rs+c5*rs2+c6*rs3
           cc  = -cx + q6/q7
           r0  = (sk/fk)**2
           r1  = a4*r0*g4
           coeff = cc-cc0-3.d0*cx/7.d0
           r2  = xnu*coeff*g3
!           a4ms = -a4*g4*s2                                !aas
           a4ms = -r1*t2                                    !aas
           r3  = dexp(a4ms)                       
           h0  = g3*(bet/delt)*dlog(1+delt*q4*t2/q5)
           h1  = r3*r2*t2
           h   = h0+h1
           
           q8  = q5*q5+delt*q4*q5*t2
           h0t = 2*bet*t*(1+2*b*t2)/q8 * g3
           h0b = -bet*t6*(2*b+b2*t2)/q8 * g3
!           h0rs = h0b * b2 * eurs / bet / g3               !aas
           h0rs = h0b * b * eurs * (b+delt) / bet / g3      !aas
           h1t  = 2*r3*r2*t*(1.d0-r1*t2)
           ccrs = (c2+2*c3*rs)/q7-q6*(c4+2*c5*rs+3*c6*rs2)/q7**2
!           r1rs = 100*r0/rs                                ! aas
           r1rs = r1/rs                                     ! aas
           h1rs = xnu*t2*r3*(ccrs - coeff*t2*r1rs) * g3
           
           ht   = h0t+h1t
           hrs  = h0rs + h1rs

           if(dabs(zeta) < zeta_minimum) then
              gzt = -2.d0/9*zeta *( 1 +14.d0/27 *zeta*zeta &
                   &             *( 1 +13.d0/18 *zeta*zeta ))
           else
              gzt = ( 1.d0/onpzeta**thrd -1.d0/onmzeta**thrd )/3.d0
           end if
!           bzt = pon*b**2 * (euzt-3*eu*gzt/g) / ( bet * g3 )           !aas
           bzt = b * ( b+delt) * ( euzt - 3*eu*gzt/g ) / ( bet * g3 )   !aas
           h0zt = 3 * gzt * h0 / g + h0b * bzt
!           h1zt =(3*gzt + 4*gzt*a4ms ) * h1 / g                        !aas
           h1zt = (3.d0 + 4.d0*a4ms) * h1 * gzt / g                    !aas

           ec1  = d*h * 0.5
! --> revised by Kato-san,  1999Jul16
!!$           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + h0zt * zetadxd
!           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxd
           ec1da = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*ht/g + h1zt ) * zetadxda
           ec1db = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*ht/g + h1zt ) * zetadxdb
           ec1dd = 0.5*ht/sk / g

10         continue
           exc0 = eu*d * 0.5 + ec1
!           excd = ecd + ec1d
           excda= ecda + ec1da
           excdb= ecdb + ec1db
!           dF_drho(i, is) = dF_drho(i, is) + excd
           dF_drho(i, 1) = dF_drho(i, 1) + excda
           dF_drho(i, 2) = dF_drho(i, 2) + excdb
           excdd = ec1dd
! gradient of charge density 95/12/3 H.S.
!           if(is ==  2) then
              if(dabs(grad_trho(i)) > 1.0d-9) then
!                 dF_dgradrho(i) = excdd / grad_trho(i)
                 dF_dgradrho(i) = excdd
              else
                 dF_dgradrho(i) = 0.d0
              endif
!           endif
! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d * 0.5
!           exc1 = exc1 + exc0*wos(i)
           exc(i)=exc(i)+exc0*2.d0
!print '(i3,5e19.12)', i,exc0,dF_drho(i, 1:2),dF_dgradrho(i),excdd/ grad_trho(i)
      
! ***** Calc heigher order derivatives. *****

           rsdn=-thrd*rs/d
           rsdnn=-4.d0*thrd*rsdn/d
           rsdnnn=-thrd7*rsdnn/d
           
           id2=1.d0/d/d
           id3=id2/d
           id4=id3/d
           zetada   =     zetadxda/d
           zetadb   =     zetadxdb/d
           zetadaa  =    -4.d0*chgrhr_l(i, 2)*id3
           zetadab  =     2.d0*zeta*id2
           zetadbb  =     4.d0*chgrhr_l(i, 1)*id3
           zetadaaa =    12.d0*chgrhr_l(i, 2)*id4
           zetadaab =   ( -4.d0*chgrhr_l(i, 1) + 8.d0*chgrhr_l(i, 2))*id4
           zetadabb =   ( -8.d0*chgrhr_l(i, 1) + 4.d0*chgrhr_l(i, 2))*id4
           zetadbbb =   -12.d0*chgrhr_l(i, 1)*id4
           
!print '(i2,7d19.10)', i,chgrhr_l(i,1:2),zetadab,zetadaaa,zetadaab,zetadabb,zetadbbb
             
           q1drr =0.5d0*a*( -b1/rs32+3.d0*b3/rs12+8.d0*b4 )
           q1drrr=0.75d0*a*( b1/rs - b3 )/rs32
           invq12pq1=1.d0/(q1**2 + q1)
           q12p1=2.d0*q1+1.d0
           e1dr=-2.d0*a*a1*q2-q0*q3*invq12pq1
           e1drr=4.d0*a*a1*q3*invq12pq1-q0*( q1drr-q12p1*q3**2*invq12pq1 )*invq12pq1
           e1drrr=6.d0*a*a1*( q1drr - q12p1*q3**2*invq12pq1 ) * invq12pq1 &
                    -q0*( q1drrr -(2.d0*q3**3+3.d0*q12p1*q3*q1drr -  &
                    2.d0*q12p1**2*q3**3* invq12pq1)*invq12pq1 ) * invq12pq1
                    
           q1pdrr =0.5d0*ap*( -b1p/rs32+3.d0*b3p/rs12+8.d0*b4p )
           q1pdrrr=0.75d0*ap*( b1p/rs - b3p )/rs32
           invq1p2pq1p=1.d0/(q1p**2 + q1p)
           q1p2p1=2.d0*q1p+1.d0
           e1pdr=-2.d0*ap*a1p*q2p-q0p*q3p*invq1p2pq1p
           e1pdrr=4.d0*ap*a1p*q3p*invq1p2pq1p-q0p*( q1pdrr-q1p2p1*q3p**2*invq1p2pq1p )*invq1p2pq1p
           e1pdrrr=6.d0*ap*a1p*( q1pdrr - q1p2p1*q3p**2*invq1p2pq1p ) * invq1p2pq1p &
                    -q0p*( q1pdrrr -(2.d0*q3p**3+3.d0*q1p2p1*q3p*q1pdrr -  &
                    2.d0*q1p2p1**2*q3p**3* invq1p2pq1p)*invq1p2pq1p ) * invq1p2pq1p
                    
           q1qdrr =0.5d0*aq*( -b1q/rs32+3.d0*b3q/rs12+8.d0*b4q )
           q1qdrrr=0.75d0*aq*( b1q/rs - b3q )/rs32
           invq1q2pq1q=1.d0/(q1q**2 + q1q)
           q1q2p1=2.d0*q1q+1.d0
           e1qdr=-2.d0*aq*a1q*q2q-q0q*q3q*invq1q2pq1q
           e1qdrr=4.d0*aq*a1q*q3q*invq1q2pq1q-q0q*( q1qdrr-q1q2p1*q3q**2*invq1q2pq1q )*invq1q2pq1q
           e1qdrrr=6.d0*aq*a1q*( q1qdrr - q1q2p1*q3q**2*invq1q2pq1q ) * invq1q2pq1q &
                    -q0q*( q1qdrrr -(2.d0*q3q**3+3.d0*q1q2p1*q3q*q1qdrr -  &
                    2.d0*q1q2p1**2*q3q**3* invq1q2pq1q)*invq1q2pq1q ) * invq1q2pq1q
           
           zeta2=zeta**2
           zeta3=zeta2*zeta
           zeta4=zeta2*zeta2
           
           if(dabs(zeta) < zeta_minimum) then
              fzdd0 = fzdd &
                   &            *( 1 + 5.d0*ninth *zeta2 &
                   &            *( 1 +22.d0*twntysvnth *zeta2 ))
              fzddd = fzdd*10.d0*ninth*zeta*( 1.d0 + 44.d0*twntysvnth*zeta2 )
              gzz   = -2.d0*ninth*( 1.d0 + 14.d0*ninth*zeta2*( 1.d0 + 65.d0*fftyfrth*zeta2 ) )
              gzzz  = -56.d0*eityonth*zeta*( 1.d0 + 65.d0*twntysvnth*zeta2 )
           else
              opzthrdm2=1.d0/onpzeta**(thrd2)
              omzthrdm2=1.d0/onmzeta**(thrd2)
              fzdd0 =  fzdd*0.5d0*( opzthrdm2 + omzthrdm2 )
              fzddd = -fzdd*thrd*( opzthrdm2/onpzeta - omzthrdm2/onmzeta )
              gzz   = -ninth*( opzthrdm2**2 + omzthrdm2**2 )
!              gzz   = -ninth*( 1.d0/onpzeta**(thrd4) + 1.d0/onmzeta**(thrd4) )
              gzzz  = 4.d0*twntysvnth*( opzthrdm2**2/onpzeta - omzthrdm2**2/onmzeta )
           end if
           
           e1=q0*q2
           e2=q0p*q2p
           e3=q0q*q2q
           
           eudrr    = e1drr  + ( e1pdrr  -  e1drr  ) * fz * zeta4 - e1qdrr  * fz * onzt4 / fzdd
           eudrrr   = e1drrr + ( e1pdrrr -  e1drrr ) * fz * zeta4 - e1qdrrr * fz * onzt4 / fzdd
           eudzz    = ( e2 - e1 )*( fzdd0*zeta4 + 8.d0*fzd*zeta3 + 12.d0*fz*zeta2 ) &
                        - e3*( fzdd0*onzt4 - 8.d0*fzd*zeta3 - 12.d0*fz*zeta2 )/fzdd
           eudzzz   = ( e2 - e1 )*( fzddd*zeta4 + 12.d0*fzdd0*zeta3 + 36.d0*fzd*zeta2 + 24.d0*fz*zeta ) &
                        - e3*( fzddd*onzt4 - 12.d0*fzdd0*zeta3 - 36.d0*fzd*zeta2 - 24.d0*fz*zeta )/fzdd
           eudrz    = ( e1pdr - e1dr )*( fzd*zeta4 + 4.d0*fz*zeta3 ) &
                        - e1qdr *( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd 
           eudrrz   = ( e1pdrr  - e1drr  )*( fzd*zeta4 + 4.d0*fz*zeta3 ) &
                        - e1qdrr *( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd
           eudrzz   = ( e1pdr - e1dr )*( fzdd0*zeta4 + 8.d0*fzd*zeta3 + 12.d0*fz*zeta2 ) &
                        - e1qdr *( fzdd0*onzt4 - 8.d0*fzd*zeta3 - 12.d0*fz*zeta2 )/fzdd
                        
           rsdn2    = rsdn*rsdn
           rsdn3    = rsdn2*rsdn
           zetada2  = zetada**2
           zetada3  = zetada2*zetada
           zetadb2  = zetadb**2
           zetadb3  = zetadb2*zetadb
                        
           euda     = eurs*rsdn + euzt*zetada
           eudb     = eurs*rsdn + euzt*zetadb
           eudaa    = eudrr*rsdn2 + 2.d0*eudrz*rsdn*zetada + eudzz*zetada2 + eurs*rsdnn + euzt*zetadaa
           eudab    = eudrr*rsdn2 + eudrz*rsdn*( zetada + zetadb ) + eudzz*zetada*zetadb &
                        + eurs*rsdnn + euzt*zetadab
           eudbb    = eudrr*rsdn2 + 2.d0*eudrz*rsdn*zetadb + eudzz*zetadb2 + eurs*rsdnn + euzt*zetadbb
           eudaaa   = eudrrr*rsdn3 &
                        + 3.d0*eudrrz*rsdn2*zetada &
                        + 3.d0*eudrzz*rsdn*zetada2 + eudzzz*zetada3 &
                        + 3.d0*eudrr*rsdn*rsdnn + 3.d0*eudzz*zetada*zetadaa &
                        + 3.d0*eudrz*( rsdnn*zetada + rsdn*zetadaa ) &
                        + eurs*rsdnnn + euzt*zetadaaa
           eudaab   = eudrrr*rsdn3 &
                        + eudrrz*rsdn2*( zetadb + 2.d0*zetada ) &
                        + eudrzz*rsdn*zetada*( zetada + 2.d0*zetadb ) &
                        + eudzzz*zetada2*zetadb &
                        + 3.d0*eudrr*rsdn*rsdnn &
                        + eudrz*( rsdnn*( 2.d0*zetada + zetadb ) + rsdn*( 2.d0*zetadab + zetadaa ) ) &
                        + eudzz*( 2.d0*zetada*zetadab + zetadaa*zetadb ) &
                        + eurs*rsdnnn + euzt*zetadaab
           eudabb   = eudrrr*rsdn3 &
                        + eudrrz*rsdn2*( zetada + 2.d0*zetadb ) &
                        + eudrzz*rsdn*zetadb*( zetadb + 2.d0*zetada ) &
                        + eudzzz*zetada*zetadb2 &
                        + 3.d0*eudrr*rsdn*rsdnn &
                        + eudrz*( rsdnn*( zetada + 2.d0*zetadb ) + rsdn*( 2.d0*zetadab + zetadbb ) ) &
                        + eudzz*( zetada*zetadbb + 2.d0*zetadb*zetadab ) &
                        + eurs*rsdnnn + euzt*zetadabb
           eudbbb   = eudrrr*rsdn3 &
                        + 3.d0*eudrrz*rsdn2*zetadb &
                        + 3.d0*eudrzz*rsdn*zetadb2 + eudzzz*zetadb3 &
                        + 3.d0*eudrr*rsdn*rsdnn + 3.d0*eudzz*zetadb*zetadbb &
                        + 3.d0*eudrz*( rsdnn*zetadb + rsdn*zetadbb ) &
                        + eurs*rsdnnn + euzt*zetadbbb
                        
!print '(i2,12d19.10)', i,chgrhr_l(i,1:2),eu,euda,eudb,eudaa,eudbb,eudab,eudaaa,eudaab,eudabb,eudbbb
!print *,eu,     e1+(e2-e1)* fz * zeta4 - e3  * fz * onzt4 / fzdd
!print *,eurs,   e1dr  + ( e1pdr  -  e1dr  ) * fz * zeta4 - e1qdr  * fz * onzt4 / fzdd
!print *,euzt,   (e2-e1)*( fzd*zeta4 + 4.d0*fz*zeta3 )-e3*( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd 

           duda     = gzt*zetada
           dudb     = gzt*zetadb 
           dudaa    = gzz*zetada2     + gzt*zetadaa
           dudab    = gzz*zetada*zetadb + gzt*zetadab
           dudbb    = gzz*zetadb2     + gzt*zetadbb
           dudaaa   = gzzz*zetada3    + 3.d0*gzz*zetada*zetadaa   + gzt*zetadaaa
           dudaab   = gzzz*zetada2*zetadb + gzz*( zetadaa*zetadb + 2.d0*zetada*zetadab ) + gzt*zetadaab
           dudabb   = gzzz*zetada*zetadb2 + gzz*( zetada*zetadbb + 2.d0*zetadab*zetadb ) + gzt*zetadabb
           dudbbb   = gzzz*zetadb3    + 3.d0*gzz*zetadb*zetadbb   + gzt*zetadbbb 
           
!print '(i2,12d19.10)', i,chgrhr_l(i,1:2),g,duda,dudb,dudaa,dudbb,dudab,dudaaa,dudaab,dudabb,dudbbb

           
           dtdn  = -sixth7*t/d
           dtdnn = -sixth13*dtdn/d
           dtdnnn= -sixth19*dtdnn/d
!print '(i2,5e19.12)',i,d,t,dtdn,dtdnn,dtdnnn
           if(dabs(grad_trho(i)) > 1.0d-9) then
               dtdg  = t/dd
               dtdng = -sixth7*dtdg/d
               dtdnng= -sixth13*dtdng/d
!               dtdnug = -dtdn*invg/dd
           else
               dtdg  = 0.d0
               dtdng = 0.d0
               dtdnng= 0.d0
               dtdnug= 0.d0 
           endif
           dtdgg  = 0.d0
           dtdngg = 0.d0
           dtdggg = 0.d0
           invg=1.d0/g
           dtdu   = -t*invg
           dtduu  = -2.d0*dtdu*invg 
           dtduuu = -3.d0*dtduu*invg
           dtdnu  = -dtdn*invg
           dtdnnu = -dtdnn*invg
           dtdnuu = -sixth7*dtduu/d
           dtdug  = -dtdg*invg
           dtduug = 2.d0*dtdg*invg**2
           dtdugg = 0.d0
           if(dabs(grad_trho(i)) > 1.0d-9) then
                dtdnug = dtdnu/dd
           else
                dtdnug = 0.d0
           end if
           
           dtda   = dtdn + dtdu*duda
           dtdb   = dtdn + dtdu*dudb
           dtdaa  = dtdnn + 2.d0*dtdnu*duda + dtduu*duda**2 + dtdu*dudaa
           dtdbb  = dtdnn + 2.d0*dtdnu*dudb + dtduu*dudb**2 + dtdu*dudbb
           dtdab  = dtdnn + dtdnu*( duda + dudb ) + dtduu*duda*dudb + dtdu*dudab
           dtdag  = dtdng + dtdug*duda
           dtdbg  = dtdng + dtdug*dudb
           dtdaaa = dtdnnn+ 3.d0*dtdnnu*duda + 3.d0*dtdnuu*duda**2 + dtduuu*duda**3 &
                    + 3.d0*dtduu*duda*dudaa + 3.d0*dtdnu*dudaa + dtdu*dudaaa
           dtdbbb = dtdnnn+ 3.d0*dtdnnu*dudb + 3.d0*dtdnuu*dudb**2 + dtduuu*dudb**3 &
                    + 3.d0*dtduu*dudb*dudbb + 3.d0*dtdnu*dudbb + dtdu*dudbbb
           dtdaab = dtdnnn+ dtdnnu*( 2.d0*duda + dudb ) + dtdnuu*duda*( duda + 2.d0*dudb ) + dtduuu*duda**2*dudb &
                    + dtduu*( dudaa*dudb + 2.d0*duda*dudab ) + dtdnu*( dudaa + 2.d0*dudab ) + dtdu*dudaab
           dtdabb = dtdnnn+ dtdnnu*( duda + 2.d0*dudb ) + dtdnuu*dudb*( dudb + 2.d0*duda ) + dtduuu*dudb**2*duda &
                    + dtduu*( dudbb*duda + 2.d0*dudb*dudab ) + dtdnu*( dudbb + 2.d0*dudab ) + dtdu*dudabb 
           dtdaag = dtdnng+ 2.d0*dtdnug*duda + dtduug*duda**2 + dtdug*dudaa
           dtdbbg = dtdnng+ 2.d0*dtdnug*dudb + dtduug*dudb**2 + dtdug*dudbb
           dtdabg = dtdnng+ dtdnug*( duda + dudb ) + dtduug*duda*dudb + dtdug*dudab 
           
!print '(i2,20d19.10)', i,t,dtda,dtdb,dtdg,dtdaa,dtdbb,dtdgg,dtdab,dtdag,dtdbg,dtdaaa,dtdbbb,dtdggg,dtdaab,dtdabb,dtdaag,dtdagg,dtdbbg,dtdbgg,dtdabg     
!print '(i2,12d19.10)', i,t,dtdag,dtdbg,dtdaag,dtdbbg,dtdabg     

           invg3  = 1.d0  /g3
           dAde   = b * ( b + delt ) / bet * invg3
           dAdee  = 2.d0 * dAde * ( b + 0.5d0 * delt ) / bet *invg3
           dAdeee = 2.d0 * ( dAde**2 + ( b + 0.5d0 * delt )*dAdee ) /bet * invg3
           m3u    = -3.d0*invg
           dAdu     = m3u*eu*dAde
           dAdeu    = m3u*( dAde + eu*dAdee )
           dAduu    = m3u*eu*invg*( -dAde + g*dAdeu )
           dAdeeu   = m3u*( 2.d0*dAdee + eu*dAdeee )
           dAdeuu   = m3u*invg*( -dAde + g*dAdeu - eu*dAdee + g*eu*dAdeeu )
           dAduuu   = m3u*eu*invg**2*( 2.d0*dAde -2.d0*g*dAdeu + g**2*dAdeuu )
           
           dAda     = dAde*euda + dAdu*duda
           dAdb     = dAde*eudb + dAdu*dudb
           dAdaa    = dAdee*euda**2 + 2.d0*dAdeu*euda*duda + dAduu*duda**2 + dAde*eudaa + dAdu*dudaa 
           dAdab    = dAdee*euda*eudb + dAdeu*( euda*dudb + eudb*duda ) + dAduu*duda*dudb + dAde*eudab + dAdu*dudab 
           dAdbb    = dAdee*eudb**2 + 2.d0*dAdeu*eudb*dudb + dAduu*dudb**2 + dAde*eudbb + dAdu*dudbb
           dAdaaa   = dAdeee*euda**3 + 3.d0*dAdeeu*euda**2*duda + 3.d0*dAdeuu*euda*duda**2 + dAduuu*duda**3 &
                        + 3.d0*dAdee*euda*eudaa + 3.d0*dAduu*duda*dudaa &
                        + 3.d0*dAdeu*( euda*dudaa + eudaa*duda ) + dAde*eudaaa + dAdu*dudaaa
           dAdaab   = dAdeee*euda**2*eudb + dAduuu*duda**2*dudb &
                        + dAdeeu*euda*( euda*dudb + 2.d0*eudb*duda ) &
                        + dAdeuu*( 2.d0*euda*dudb + eudb*duda )*duda &
                        + dAdee*( 2.d0*euda*eudab + eudb*eudaa ) + dAduu*( 2.d0*duda*dudab + dudb*dudaa ) &
                        + dAdeu*( 2.d0*eudab*duda + 2.d0*euda*dudab + dudb*eudaa + eudb*dudaa ) &
                        + dAde*eudaab + dAdu*dudaab
           dAdabb   = dAdeee*euda*eudb**2 + dAduuu*dudb**2*duda &
                        + dAdeeu*eudb*( eudb*duda + 2.d0*euda*dudb ) &
                        + dAdeuu*( 2.d0*eudb*duda + euda*dudb )*dudb &
                        + dAdee*( eudbb*euda + 2.d0*eudb*eudab ) + dAduu*( dudbb*duda + 2.d0*dudb*dudab ) &
                        + dAdeu*( eudbb*duda + 2.d0*eudb*dudab + 2.d0*dudb*eudab + dudbb*euda ) &
                        + dAde*eudabb + dAdu*dudabb
           dAdbbb   = dAdeee*eudb**3 + 3.d0*dAdeeu*eudb**2*dudb + 3.d0*dAdeuu*eudb*dudb**2 + dAduuu*dudb**3 &
                        + 3.d0*dAdee*eudb*eudbb + 3.d0*dAduu*dudb*dudbb &
                        + 3.d0*dAdeu*( eudb*dudbb + eudbb*dudb ) + dAde*eudbbb + dAdu*dudbbb 
                        
!print '(i2,12d19.10)', i,b,dAda,dAdb,dAdaa,dAdab,dAdbb,dAdaaa,dAdaab,dAdabb,dAdbbb

           at2  = b * t2
           a2t4 = at2**2
           a3t6 = at2*a2t4
           t3   = t2*t
           t5   = t4*t
           t7   = t6*t
           t8   = t7*t
           fE    = q5**2 + delt*t2*q4*q5
           fEdt  = 4.d0*b*t*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt *2.d0*t*(1.d0 + 4.d0*at2 + 6.d0*a2t4 + 4.d0*a3t6 )
           fEdtt = 4.d0*b*( 1.d0 + 9.d0*at2 + 15.d0*a2t4 + 14.d0*a3t6 ) + &
                    delt*( 2.d0 + 24.d0*at2 + 60.d0*a2t4 + 56.d0*a3t6 )
           fEda  = 2.d0*t2*q5*( 1.d0 + 2.d0*at2 ) + &
                    delt*t4*( 2.d0 + 4.d0*at2 + 3.d0*a2t4 )
           fEdaa = 6.d0*t4*( 2.d0*q5 - 1.d0) + &
                    delt*2.d0*t6*( 2.d0 + 3.d0*at2 )
           fF    = t*( 1.d0 + 2.d0*at2 )         
           fFdt  = 1.d0 + 6.d0*at2
           fFdtt = 12.d0*b*t
           fFda  = 2.d0*t3
           fFdaa = 0.d0
           fG    = 2.d0*b*t6 + b2*t8
           fGdt  = 12.d0*b*t5 + 8.d0*b2*t7
           fGdtt = 60.d0*b*t4 + 56.d0*b2*t6
           fGda  = 2.d0*t6 + 2.d0*b*t8
           fGdaa = 2.d0*t8
           
           betivE=bet/fE*g3
           betivE2=betivE/fE
           dHdf   =     -fG*betivE
           dHdt   = 2.d0*fF*betivE
           dHdff  =     -( fGda*fE-fG*fEda )*betivE2
           dHdft  =     -( fGdt*fE-fG*fEdt )*betivE2
           dHdtt  = 2.d0*( fFdt*fE-fF*fEdt )*betivE2
           dHdfff =     -( fGdaa*fE-fG*fEdaa-2.d0*( fGda*fE-fG*fEda )*fEda/fE )*betivE2
           dHdfft = 2.d0*( fFdaa*fE-fF*fEdaa-2.d0*( fFda*fE-fF*fEda )*fEda/fE )*betivE2
           dHdftt =     -( fGdtt*fE-fG*fEdtt-2.d0*( fGdt*fE-fG*fEdt )*fEdt/fE )*betivE2
           dHdttt = 2.d0*( fFdtt*fE-fF*fEdtt-2.d0*( fFdt*fE-fF*fEdt )*fEdt/fE )*betivE2
           
           tiu    = 3.d0*invg
           siu2   = tiu*2.d0*invg
           dHdu   = tiu*h0
           dHduu  = siu2*h0
           dHduuu = siu2*invg*h0
           dHdfu  = tiu*dHdf
           dHdtu  = tiu*dHdt
           dHdffu = tiu*dHdff
           dHdfuu = siu2*dHdf
           dHdttu = tiu*dHdtt
           dHdtuu = siu2*dHdt
           dHdftu = tiu*dHdft
           
           dHda   = dHdf*dAda + dHdt*dtda + dHdu*duda
           dHdb   = dHdf*dAdb + dHdt*dtdb + dHdu*dudb
           dHdg   =             dHdt*dtdg
           dHdaa  = dHdff*dAda**2 + dHdtt*dtda**2 + dHduu*duda**2 &
                    + 2.d0*dHdft*dAda*dtda + 2.d0*dHdfu*dAda*duda + 2.d0*dHdtu*dtda*duda &
                    + dHdf*dAdaa + dHdt*dtdaa + dHdu*dudaa
           dHdab  = dHdff*dAda*dAdb + dHdtt*dtda*dtdb + dHduu*duda*dudb &
                    + dHdft*( dAda*dtdb + dAdb*dtda ) + dHdfu*( dAda*dudb + dAdb*duda ) + dHdtu*( dtda*dudb + dtdb*duda ) &
                    + dHdf*dAdab + dHdt*dtdab + dHdu*dudab
           dHdbb  = dHdff*dAdb**2 + dHdtt*dtdb**2 + dHduu*dudb**2 &
                    + 2.d0*dHdft*dAdb*dtdb + 2.d0*dHdfu*dAdb*dudb + 2.d0*dHdtu*dtdb*dudb &
                    + dHdf*dAdbb + dHdt*dtdbb + dHdu*dudbb
           dHdgg  = dHdtt*dtdg**2
           dHdag  = dHdft*dAda*dtdg + dHdtt*dtda*dtdg + dHdtu*duda*dtdg + dHdt*dtdag
           dHdbg  = dHdft*dAdb*dtdg + dHdtt*dtdb*dtdg + dHdtu*dudb*dtdg + dHdt*dtdbg
           dHdaaa = dHdfff*dAda**3 + 3.d0*dHdfft*dAda**2*dtda + 3.d0*dHdffu*dAda**2*duda &
                    + 3.d0*dHdftt*dAda*dtda**2 + dHdttt*dtda**3 + 3.d0*dHdttu*dtda**2*duda &
                    + 3.d0*dHdfuu*dAda*duda**2 + 3.d0*dHdtuu*dtda*duda**2 + dHduuu*duda**3 &
                    + 6.d0*dHdftu*dAda*dtda*duda &
                    + 3.d0*dHdff*dAda*dAdaa + 3.d0*dHdtt*dtda*dtdaa + 3.d0*dHduu*duda*dudaa &
                    + 3.d0*dHdft*( dAdaa*dtda + dAda*dtdaa ) &
                    + 3.d0*dHdtu*( dtdaa*duda + dtda*dudaa ) &
                    + 3.d0*dHdfu*( dAdaa*duda + dAda*dudaa ) &
                    + dHdf*dAdaaa + dHdt*dtdaaa + dHdu*dudaaa
           dHdaab = dHdfff*dAda**2*dAdb + dHdfft*dAda*( dAda*dtdb + 2.d0*dAdb*dtda ) + &
                      dHdffu*dAda*( dAda*dudb + 2.d0*dAdb*duda ) &
                    + dHdftt*dtda*( dAdb*dtda + 2.d0*dAda*dtdb ) &
                    + dHdttt*dtdb*dtda**2 + dHdttu*dtda*( dtda*dudb + 2.d0*dtdb*duda ) &
                    + dHdfuu*duda*( dAdb*duda + 2.d0*dAda*dudb ) &
                    + dHdtuu*duda*( dtdb*duda + 2.d0*dtda*dudb ) + dHduuu*duda**2*dudb &
                    + 2.d0*dHdftu*( dAdb*duda*dtda + dAda*dudb*dtda + dAda*duda*dtdb ) &
                    + dHdff*( 2.d0*dAda*dAdab + dAdb*dAdaa ) &
                    + dHdtt*( 2.d0*dtda*dtdab + dtdb*dtdaa ) &
                    + dHduu*( 2.d0*duda*dudab + dudb*dudaa ) &
                    + dHdft*( 2.d0*dAdab*dtda + 2.d0*dAda*dtdab + dAdaa*dtdb + dAdb*dtdaa ) &
                    + dHdtu*( 2.d0*dtdab*duda + 2.d0*dtda*dudab + dtdaa*dudb + dtdb*dudaa ) &
                    + dHdfu*( 2.d0*dAdab*duda + 2.d0*dAda*dudab + dAdaa*dudb + dAdb*dudaa ) &
                    + dHdf*dAdaab + dHdt*dtdaab + dHdu*dudaab
           dHdabb = dHdfff*dAda*dAdb**2 + dHdfft*dAdb*( dAdb*dtda + 2.d0*dAda*dtdb ) &
                    + dHdffu*dAdb*( dAdb*duda + 2.d0*dAda*dudb ) &
                    + dHdftt*dtdb*( dAda*dtdb + 2.d0*dAdb*dtda ) &
                    + dHdttt*dtda*dtdb**2 + dHdttu*dtdb*( dtdb*duda + 2.d0*dtda*dudb ) &
                    + dHdfuu*dudb*( dAda*dudb + 2.d0*dAdb*duda ) &
                    + dHdtuu*dudb*( dtda*dudb + 2.d0*dtdb*duda ) + dHduuu*duda*dudb**2 &
                    + 2.d0*dHdftu*( dAda*dtdb*dudb + dAdb*dtda*dudb + dAdb*dtdb*duda ) &
                    + dHdff*( 2.d0*dAdb*dAdab + dAda*dAdbb ) &
                    + dHdtt*( 2.d0*dtdb*dtdab + dtda*dtdbb ) &
                    + dHduu*( 2.d0*dudb*dudab + duda*dudbb ) &
                    + dHdft*( 2.d0*dAdab*dtdb + 2.d0*dAdb*dtdab + dAdbb*dtda + dAda*dtdbb ) &
                    + dHdtu*( 2.d0*dtdab*dudb + 2.d0*dtdb*dudab + dtdbb*duda + dtda*dudbb ) &
                    + dHdfu*( 2.d0*dAdab*dudb + 2.d0*dAdb*dudab + dAdbb*duda + dAda*dudbb ) &
                    + dHdf*dAdabb + dHdt*dtdabb + dHdu*dudabb
           dHdbbb = dHdfff*dAdb**3 + 3.d0*dHdfft*dAdb**2*dtdb + 3.d0*dHdffu*dAdb**2*dudb &
                    + 3.d0*dHdftt*dAdb*dtdb**2 + dHdttt*dtdb**3 + 3.d0*dHdttu*dtdb**2*dudb &
                    + 3.d0*dHdfuu*dAdb*dudb**2 + 3.d0*dHdtuu*dtdb*dudb**2 + dHduuu*dudb**3 &
                    + 6.d0*dHdftu*dAdb*dtdb*dudb &
                    + 3.d0*dHdff*dAdb*dAdbb + 3.d0*dHdtt*dtdb*dtdbb + 3.d0*dHduu*dudb*dudbb &
                    + 3.d0*dHdft*( dAdbb*dtdb + dAdb*dtdbb ) &
                    + 3.d0*dHdtu*( dtdbb*dudb + dtdb*dudbb ) &
                    + 3.d0*dHdfu*( dAdbb*dudb + dAdb*dudbb ) &
                    + dHdf*dAdbbb + dHdt*dtdbbb + dHdu*dudbbb
           dHdabg = dHdfft*dAda*dAdb*dtdg   + dHdftt*( dAdb*dtda + dAda*dtdb )*dtdg &
                    + dHdttt*dtda*dtdb*dtdg + dHdttu*( dtdb*duda + dtda*dudb )*dtdg &
                    + dHdtuu*duda*dudb*dtdg + dHdftu*( dAdb*duda + dAda*dudb )*dtdg &
                    + dHdft*( dAdab*dtdg + dAdb*dtdag + dAda*dtdbg ) &
                    + dHdtt*( dtdab*dtdg + dtdb*dtdag + dtda*dtdbg ) &
                    + dHdtu*( dudab*dtdg + dudb*dtdag + duda*dtdbg ) &
                    + dHdt*dtdabg
           dHdaag = dHdfft*dAda**2*dtdg + dHdttt*dtda**2*dtdg + dHdtuu*duda**2*dtdg &
                    + 2.d0*dHdttu*dtda*duda*dtdg + 2.d0*dHdftt*dAda*dtda*dtdg + 2.d0*dHdftu*dAda*duda*dtdg &
                    + dHdtt*( dtdaa*dtdg + 2.d0*dtda*dtdag ) + dHdft*( dAdaa*dtdg + 2.d0*dAda*dtdag ) &
                    + dHdtu*( dudaa*dtdg + 2.d0*duda*dtdag ) + dHdt*dtdaag         
           dHdbbg = dHdfft*dAdb**2*dtdg + dHdttt*dtdb**2*dtdg + dHdtuu*dudb**2*dtdg &
                    + 2.d0*dHdttu*dtdb*dudb*dtdg + 2.d0*dHdftt*dAdb*dtdb*dtdg + 2.d0*dHdftu*dAdb*dudb*dtdg &
                    + dHdtt*( dtdbb*dtdg + 2.d0*dtdb*dtdbg ) + dHdft*( dAdbb*dtdg + 2.d0*dAdb*dtdbg ) &
                    + dHdtu*( dudbb*dtdg + 2.d0*dudb*dtdbg ) + dHdt*dtdbbg       
           dHdagg = dHdftt*dAda*dtdg**2 + dHdttt*dtda*dtdg**2 + dHdttu*duda*dtdg**2 + 2.d0*dHdtt*dtdg*dtdag
           dHdbgg = dHdftt*dAdb*dtdg**2 + dHdttt*dtdb*dtdg**2 + dHdttu*dudb*dtdg**2 + 2.d0*dHdtt*dtdg*dtdbg 
           dHdggg = dHdttt*dtdg**3   
           
! ****** Derivative of h1 *****
           
           dq6   = c2+2.d0*c3*rs
           ddq6  = 2.d0*c3
           dq7   = c4+2.d0*c5*rs+3.d0*c6*rs2
           ddq7  = 2.d0*c5+6.d0*c6*rs
           dddq7 = 6.d0*c6
           
           iq7=1.d0/q7
           iq72=iq7*iq7
           iq73=iq72*iq7
           iq74=iq73*iq7
           
           ccrs2 =  ddq6*iq7-(2.d0*dq6*dq7+q6*ddq7)*iq72+2.d0*q6*dq7**2*iq73
           ccrs3 =  -iq72*(3.d0*(ddq6*dq7+dq6*ddq7)+q6*dddq7) + &
                                    6.d0*dq7*(dq6*dq7+q6*ddq7)*iq73 - &
                                    6.d0*q6*dq7**3*iq74
           r1rs2=r1rs*r1rs
           r1rs3=r1rs2*r1rs
           r1rs4=r1rs3*r1rs
           
           if(dabs(grad_trho(i)) > 1.0d-9) then
               ivt=1.d0/t
               ivt2=ivt**2
           else
               ivt=0.d0
               ivt2=0.d0
           end if
           
!           ivtmr1t=ivt-r1*t
           ivu=1.d0/g
           ivu2=1.d0/g/g
           r1t2=r1*t2
           r1t22= r1t2**2
           r1t23= r1t2*r1t22
           
           dH1dr   =  h1rs
           dH1du   =  h1*ivu*(3.d0-4.d0*r1t2)
           dH1dt   =  h1t
           dH1drr  =  xnu*ccrs2*t2*g3*r3-2.d0*h1rs*r1rs*t2-h1*r1rs2*t4
           dH1duu  =  2.d0*h1*ivu2*(3.d0 - 18.d0*r1t2 + 8.d0*r1t22) 
           dH1dtt  =  2.d0*h1*ivt2*(1.d0 -  5.d0*r1t2 + 2.d0*r1t22)
           dH1dru  =  ivu*(dH1dr*(3.d0 - 4.d0*r1t2)-4.d0*h1*r1rs*t2)
           dH1dtu  =  2.d0*ivu*ivt*h1*(3.d0-11.d0*r1t2+4.d0*r1t22)
           dH1drt  =  2.d0*(ivt*dH1dr*(1.d0-r1t2)-h1*r1rs*t)
           dH1drrr =  xnu*ccrs3*t2*g3*r3-3.d0*t2*r1rs*(dH1drr+h1rs*r1rs*t2)-h1*r1rs3*t6
           dH1drru =  ivu*(dH1drr*(3.d0-4.d0*r1t2)-8.d0*t2*dH1dr*r1rs)
           dH1druu =  2.d0*ivu2*(dH1dr*(3.d0-18.d0*r1t2+8.d0*r1t22)+2.d0*h1*t2*r1rs*(8.d0*r1t2-9.d0))
           dH1duuu =  2.d0*h1*ivu*ivu2*(3.d0-102.d0*r1t2+144.d0*r1t22-32.d0*r1t23)
           dH1dtuu =  4.d0*h1*ivu2*ivt*(3.d0-39.d0*r1t2+42.d0*r1t22-8.d0*r1t23)
           dH1dttu =  0.5d0*g*ivt*dH1dtuu
           dH1dttt =  4.d0*h1*r1*ivt*(-6.d0+9.d0*r1t2-2.d0*r1t22)
           dH1drtt =  2.d0*dH1dr*ivt2*(1.d0-5.d0*r1t2+2.d0*r1t22)+2.d0*h1*r1rs*(4.d0*r1t2-5.d0)
           dH1drrt =  2.d0*ivt*dH1drr*(1.d0-r1t2)-4.d0*dH1dr*r1rs*t
           dH1drtu =  2.d0*ivu*( ivt*dH1dr*(3.d0-11.d0*r1t2+4.d0*r1t22)+h1*t*r1rs*(-11.d0+8.d0*r1t2))
           
           dH1da   = dH1dr*rsdn + dH1dt*dtda + dH1du*duda
           dH1db   = dH1dr*rsdn + dH1dt*dtdb + dH1du*dudb
           dH1dg   =              dH1dt*dtdg
           dH1daa  = dH1drr*rsdn**2 + dH1dtt*dtda**2 + dH1duu*duda**2 &
                    + 2.d0*dH1drt*rsdn*dtda + 2.d0*dH1dru*rsdn*duda + 2.d0*dH1dtu*dtda*duda &
                    + dH1dr*rsdnn + dH1dt*dtdaa + dH1du*dudaa
           dH1dab  = dH1drr*rsdn*rsdn + dH1dtt*dtda*dtdb + dH1duu*duda*dudb &
                    + dH1drt*rsdn*(dtdb + dtda) + dH1dru*rsdn*( dudb + duda ) + dH1dtu*( dtda*dudb + dtdb*duda ) &
                    + dH1dr*rsdnn + dH1dt*dtdab + dH1du*dudab
           dH1dbb  = dH1drr*rsdn**2 + dH1dtt*dtdb**2 + dH1duu*dudb**2 &
                    + 2.d0*dH1drt*rsdn*dtdb + 2.d0*dH1dru*rsdn*dudb + 2.d0*dH1dtu*dtdb*dudb &
                    + dH1dr*rsdnn + dH1dt*dtdbb + dH1du*dudbb
           dH1dgg  = dH1dtt*dtdg**2
           dH1dag  = dH1drt*rsdn*dtdg + dH1dtt*dtda*dtdg + dH1dtu*duda*dtdg + dH1dt*dtdag
           dH1dbg  = dH1drt*rsdn*dtdg + dH1dtt*dtdb*dtdg + dH1dtu*dudb*dtdg + dH1dt*dtdbg
           dH1daaa = dH1drrr*rsdn3 + 3.d0*dH1drrt*rsdn2*dtda + 3.d0*dH1drru*rsdn2*duda &
                    + 3.d0*dH1drtt*rsdn*dtda**2 + dH1dttt*dtda**3 + 3.d0*dH1dttu*dtda**2*duda &
                    + 3.d0*dH1druu*rsdn*duda**2 + 3.d0*dH1dtuu*dtda*duda**2 + dH1duuu*duda**3 &
                    + 6.d0*dH1drtu*rsdn*dtda*duda &
                    + 3.d0*dH1drr*rsdn*rsdnn + 3.d0*dH1dtt*dtda*dtdaa + 3.d0*dH1duu*duda*dudaa &
                    + 3.d0*dH1drt*( rsdnn*dtda + rsdn*dtdaa ) &
                    + 3.d0*dH1dtu*( dtdaa*duda + dtda*dudaa ) &
                    + 3.d0*dH1dru*( rsdnn*duda + rsdn*dudaa ) &
                    + dH1dr*rsdnnn + dH1dt*dtdaaa + dH1du*dudaaa
           dH1daab = dH1drrr*rsdn3 + dH1drrt*rsdn2*( dtdb + 2.d0*dtda ) + dH1drru*rsdn2*( dudb + 2.d0*duda ) &
                    + dH1drtt*dtda*rsdn*( dtda + 2.d0*dtdb ) + dH1dttt*dtdb*dtda**2 + dH1dttu*dtda*( dtda*dudb + 2.d0*dtdb*duda ) &
                    + dH1druu*duda*rsdn*( duda + 2.d0*dudb ) + dH1dtuu*duda*( dtdb*duda + 2.d0*dtda*dudb ) + dH1duuu*duda**2*dudb &
                    + 2.d0*dH1drtu*rsdn*( duda*dtda + dudb*dtda + duda*dtdb ) &
                    + dH1drr*3.d0*rsdn*rsdnn &
                    + dH1dtt*( 2.d0*dtda*dtdab + dtdb*dtdaa ) &
                    + dH1duu*( 2.d0*duda*dudab + dudb*dudaa ) &
                    + dH1drt*( 2.d0*rsdnn*dtda + 2.d0*rsdn*dtdab + rsdnn*dtdb + rsdn*dtdaa ) &
                    + dH1dtu*( 2.d0*dtdab*duda + 2.d0*dtda*dudab + dtdaa*dudb + dtdb*dudaa ) &
                    + dH1dru*( 2.d0*rsdnn*duda + 2.d0*rsdn*dudab + rsdnn*dudb + rsdn*dudaa ) &
                    + dH1dr*rsdnnn + dH1dt*dtdaab + dH1du*dudaab
           dH1dabb = dH1drrr*rsdn3 + dH1drrt*rsdn2*( dtda + 2.d0*dtdb ) + dH1drru*rsdn2*( duda + 2.d0*dudb ) &
                    + dH1drtt*dtdb*rsdn*( dtdb + 2.d0*dtda ) + dH1dttt*dtda*dtdb**2 + dH1dttu*dtdb*( dtdb*duda + 2.d0*dtda*dudb ) &
                    + dH1druu*dudb*rsdn*( dudb + 2.d0*duda ) + dH1dtuu*dudb*( dtda*dudb + 2.d0*dtdb*duda ) + dH1duuu*duda*dudb**2 &
                    + 2.d0*dH1drtu*rsdn*( dtdb*dudb + dtda*dudb + dtdb*duda ) &
                    + dH1drr*3.d0*rsdn*rsdnn &
                    + dH1dtt*( 2.d0*dtdb*dtdab + dtda*dtdbb ) &
                    + dH1duu*( 2.d0*dudb*dudab + duda*dudbb ) &
                    + dH1drt*( 2.d0*rsdnn*dtdb + 2.d0*rsdn*dtdab + rsdnn*dtda + rsdn*dtdbb ) &
                    + dH1dtu*( 2.d0*dtdab*dudb + 2.d0*dtdb*dudab + dtdbb*duda + dtda*dudbb ) &
                    + dH1dru*( 2.d0*rsdnn*dudb + 2.d0*rsdn*dudab + rsdnn*duda + rsdn*dudbb ) &
                    + dH1dr*rsdnnn + dH1dt*dtdabb + dH1du*dudabb
           dH1dbbb = dH1drrr*rsdn3 + 3.d0*dH1drrt*rsdn2*dtdb + 3.d0*dH1drru*rsdn2*dudb &
                    + 3.d0*dH1drtt*rsdn*dtdb**2 + dH1dttt*dtdb**3 + 3.d0*dH1dttu*dtdb**2*dudb &
                    + 3.d0*dH1druu*rsdn*dudb**2 + 3.d0*dH1dtuu*dtdb*dudb**2 + dH1duuu*dudb**3 &
                    + 6.d0*dH1drtu*rsdn*dtdb*dudb &
                    + 3.d0*dH1drr*rsdn*rsdnn + 3.d0*dH1dtt*dtdb*dtdbb + 3.d0*dH1duu*dudb*dudbb &
                    + 3.d0*dH1drt*( rsdnn*dtdb + rsdn*dtdbb ) &
                    + 3.d0*dH1dtu*( dtdbb*dudb + dtdb*dudbb ) &
                    + 3.d0*dH1dru*( rsdnn*dudb + rsdn*dudbb ) &
                    + dH1dr*rsdnnn + dH1dt*dtdbbb + dH1du*dudbbb
           dH1dabg = dH1drrt*rsdn2*dtdg   + dH1drtt*rsdn*( dtda + dtdb )*dtdg &
                    + dH1dttt*dtda*dtdb*dtdg + dH1dttu*( dtdb*duda + dtda*dudb )*dtdg &
                    + dH1dtuu*duda*dudb*dtdg + dH1drtu*rsdn*( duda + dudb )*dtdg &
                    + dH1drt*( rsdnn*dtdg + rsdn*dtdag + rsdn*dtdbg ) &
                    + dH1dtt*( dtdab*dtdg + dtdb*dtdag + dtda*dtdbg ) &
                    + dH1dtu*( dudab*dtdg + dudb*dtdag + duda*dtdbg ) &
                    + dH1dt*dtdabg
           dH1daag = dH1drrt*rsdn2*dtdg + dH1dttt*dtda**2*dtdg + dH1dtuu*duda**2*dtdg &
                    + 2.d0*dH1dttu*dtda*duda*dtdg + 2.d0*dH1drtt*rsdn*dtda*dtdg + 2.d0*dH1drtu*rsdn*duda*dtdg &
                    + dH1dtt*( dtdaa*dtdg + 2.d0*dtda*dtdag ) + dH1drt*( rsdnn*dtdg + 2.d0*rsdn*dtdag ) &
                    + dH1dtu*( dudaa*dtdg + 2.d0*duda*dtdag ) + dH1dt*dtdaag  
           dH1dbbg = dH1drrt*rsdn2*dtdg + dH1dttt*dtdb**2*dtdg + dH1dtuu*dudb**2*dtdg &
                    + 2.d0*dH1dttu*dtdb*dudb*dtdg + 2.d0*dH1drtt*rsdn*dtdb*dtdg + 2.d0*dH1drtu*rsdn*dudb*dtdg &
                    + dH1dtt*( dtdbb*dtdg + 2.d0*dtdb*dtdbg ) + dH1drt*( rsdnn*dtdg + 2.d0*rsdn*dtdbg ) &
                    + dH1dtu*( dudbb*dtdg + 2.d0*dudb*dtdbg ) + dH1dt*dtdbbg       
           dH1dagg = dH1drtt*rsdn*dtdg**2 + dH1dttt*dtda*dtdg**2 + dH1dttu*duda*dtdg**2 + 2.d0*dH1dtt*dtdg*dtdag
           dH1dbgg = dH1drtt*rsdn*dtdg**2 + dH1dttt*dtdb*dtdg**2 + dH1dttu*dudb*dtdg**2 + 2.d0*dH1dtt*dtdg*dtdbg 
           dH1dggg = dH1dttt*dtdg**3   
                    
!print '(i2,20d19.10)', i,h0,dHda,dHdb,dHdg,dHdaa,dHdbb,dHdgg,dHdag,dHdbg,dHdab,dHdaaa,dHdbbb,dHdggg,dHdaab,dHdabb,dHdaag,dHdagg,dHdbbg,dHdbgg,dHdabg
!print '(i2,20d19.10)', i,h1,dH1da,dH1db,dH1dg,dH1daa,dH1dbb,dH1dgg,dH1dag,dH1dbg,dH1dab,dH1daaa,dH1dbbb,dH1dggg,dH1daab,dH1dabb,dH1daag,dH1dagg,dH1dbbg,dH1dbgg,dH1dabg
    
           dFc_daa(i)   = 2.d0*euda + d*eudaa + 2.d0*(dHda + dH1da) + d*(dHdaa + dH1daa) 
           dFc_dbb(i)   = 2.d0*eudb + d*eudbb + 2.d0*(dHdb + dH1db) + d*(dHdbb + dH1dbb)
           dFc_dgg(i)   = d*(dHdgg + dH1dgg)
           dFc_dab(i)   = euda + eudb + d*eudab + dHda + dH1da + dHdb + dH1db + d*(dHdab + dH1dab) 
           dFc_dag(i)   = dHdg + dH1dg + d*(dHdag + dH1dag)
           dFc_dbg(i)   = dHdg + dH1dg + d*(dHdbg + dH1dbg)
           dFc_daaa(i)  = 3.d0*eudaa + d*eudaaa + 3.d0*(dHdaa + dH1daa) + d*(dHdaaa + dH1daaa)
           dFc_dbbb(i)  = 3.d0*eudbb + d*eudbbb + 3.d0*(dHdbb + dH1dbb) + d*(dHdbbb + dH1dbbb)
           dFc_dggg(i)  = d*(dHdggg + dH1dggg)
           dFc_daab(i)  = eudaa + 2.d0*eudab + d*eudaab + dHdaa + dH1daa + 2.d0*(dHdab + dH1dab) + d*(dHdaab + dH1daab)
           dFc_daag(i)  = 2.d0*(dHdag + dH1dag) + d*(dHdaag + dH1daag)
           dFc_dabb(i)  = eudbb + 2.d0*eudab + d*eudabb + dHdbb + dH1dbb + 2.d0*(dHdab + dH1dab) + d*(dHdabb + dH1dabb)
           dFc_dbbg(i)  = 2.d0*(dHdbg + dH1dbg) + d*(dHdbbg + dH1dbbg)
           dFc_dagg(i)  = dHdgg + dH1dgg + d*(dHdagg + dH1dagg)
           dFc_dbgg(i)  = dHdgg + dH1dgg + d*(dHdbgg + dH1dbgg)
           dFc_dabg(i)  = dHdag + dH1dag + dHdbg + dH1dbg + d*(dHdabg + dH1dabg)
                    
!print '(i2,21d19.10)', i,exc(i),d*(eu+h),dF_drho(i,1),dF_drho(i,2),dF_dgradrho(i)*grad_trho(i) &
!                      ,dFc_daa(i),dFc_dbb(i),dFc_dgg(i),dFc_dab(i),dFc_dag(i),dFc_dbg(i),dFc_daaa(i),dFc_dbbb(i) &
!                      ,dFc_dggg(i),dFc_daab(i),dFc_daag(i),dFc_dabb(i),dFc_dbbg(i),dFc_dagg(i),dFc_dbgg(i),dFc_dabg(i)
!print *,dHdag,dHdaag,d,2.d0*dHdag + d*dHdaag
!print *,2.d0*dHdbg + d*dHdbbg
!print *,dHdag + dHdbg + d*dHdabg
               
           
        end do
     end if
     
end subroutine cr_ggapw91_paw_drv3

! ==== other xc functional 
subroutine ex_gga_paw_library( nrc, dnr, nspin, chgrhr_l, grad_rho, exc, &
     &                         dFx_drho, dFx_dgradrho, &
     &                         dFx_drr,  dFx_drg,  dFx_dgg, &
     &                         dFx_drrr, dFx_drrg, dFx_drgg, dFx_dggg, pot_type )
  use m_Const_Parameters,  only : DP,PAI
  implicit none

  integer,intent(in)        :: nrc, dnr, nspin, pot_type
  real(kind=DP),intent(in)  :: chgrhr_l(nrc,nspin)
  real(kind=DP),intent(in)  :: grad_rho(nrc,nspin)

  real(kind=DP),intent(out) :: exc(nrc)
  real(kind=DP),intent(out) :: dFx_drho(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drr(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dgg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drrr(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drrg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drgg(nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dggg(nrc,nspin)

  real(kind=DP), parameter ::  ax = -0.7385587663820224d0
  real(kind=DP), parameter :: thrd = 0.33333333333333333333333d0, &
       &                      thrd4 = 1.3333333333333333333333d0
  real(kind=DP), parameter :: thpith = 3.0936677262801d0      ! ???

  real(kind=DP), parameter :: ninth2=0.2222222222222222222d0
  real(kind=DP), parameter :: ninth4=0.4444444444444444444d0
  real(kind=DP), parameter :: ninth8=0.8888888888888888888d0
  real(kind=DP), parameter :: thrd2 =0.6666666666666666666d0
  real(kind=DP), parameter :: eghth =0.125d0
  real(kind=DP), parameter :: twtysvnth8 = 0.29629629629629627985d0
  real(kind=DP), parameter :: twtysvnth32= 1.18518518518518511939d0

  ! ----- parameters of enhancement factor --
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

  real(kind=DP) :: facw, d, dd, fk, s, fac, s2, f, ex, exd, exdd, &
       &           exc0, excd, excdd, exc1, pot_add
  real(kind=DP) :: df_ds, d2f_ds2, d3f_ds3, ivd, ivdthrd2, s4, s6
  real(kind=DP) :: ctmp1, ctmp2, ctmp3, ctmp4, ctmp5, ctmp6, ctmp7, &
       &           x, dx_ds, d2x_ds2, d3x_ds3
  real(kind=DP) :: df_dx, d2f_dx2, d3f_dx3
  real(kind=DP) :: dmu_dx, d2mu_dx2, d3mu_dx3
  real(kind=DP) :: alpha_pbeint, mu_this, dmu_ds, d2mu_ds2, d3mu_ds3
  real(kind=DP) :: term(4), dterm1(4), dterm2(4), dterm3(4)
  real(kind=DP) :: y, y2, y3, y4, dy_ds, d2y_ds2, d3y_ds3

  logical :: First = .true.
  real(kind=DP), save :: coeff_spline_htbs(0:5)

  if ( pot_type == 5 ) then
     if ( First ) then
        call init_spline_htbs( mu_wc, kappa_rpbe, kappa_wc, coeff_wc, &
             &                 s_low_htbs, s_high_htbs, coeff_spline_htbs )
        First = .false.
     endif
  endif

  !---- Spin dependency

  facw = nspin

  exc  = 0.d0
  dFx_drho = 0.0d0;    dFx_dgradrho = 0.0d0;    dFx_drr  = 0.0d0
  dFx_drg  = 0.0d0;    dFx_dgg      = 0.0d0;    dFx_drrr = 0.0d0
  dFx_drrg = 0.0d0;    dFx_drgg     = 0.0d0;    dFx_dggg = 0.0d0

  do is = 1, nspin
     exc1 = 0.d0

     do i = 1,nrc,dnr
        d  = facw * chgrhr_l(i, is)
        if (d.le.0.0d0) cycle

        dd = facw * grad_rho(i, is)

        if (d > 1.d-05) then
           fk = (3.d0*PAI*PAI*d)**thrd
!!$        sk = dsqrt(4*fk/PAI)
           s = dd/(d*fk*2.d0)
        else
           fk = 0.d0
           s  = 0.d0
        endif
        !-------------------------------------
        fac = ax*d**thrd
        s2  = s*s

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
           x = mu_pbe *s2;   dx_ds = 2.0d0 *mu_pbe *s;    d2x_ds2 = 2.0d0 *mu_pbe
           d3x_ds3 = 0.0d0

           ctmp1 = kappa_pbe + x

           f = 1.0d0 +kappa_pbe *x /ctmp1
           df_dx = kappa_pbe**2 /ctmp1**2
           d2f_dx2 = -2.0d0 *df_dx /ctmp1
           d3f_dx3 = -3.0d0 *d2f_dx2 /ctmp1

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (2)                        ! revpbe
           x = mu_revpbe *s2;   dx_ds = 2.0d0 *mu_revpbe *s;  d2x_ds2 = 2.0d0 *mu_revpbe
           d3x_ds3 = 0.0d0

           ctmp1 = kappa_revpbe + x

           f = 1.0d0 +kappa_revpbe *x /ctmp1
           df_dx = kappa_revpbe**2 /ctmp1**2
           d2f_dx2 = -2.0d0 *df_dx /ctmp1
           d3f_dx3 = -3.0d0 *d2f_dx2 /ctmp1

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (3)                        ! rpbe
           x = mu_rpbe *s2;   dx_ds = 2.0d0 *mu_rpbe *s;    d2x_ds2 = 2.0d0 *mu_rpbe
           d3x_ds3 = 0.0d0

           ctmp1 = exp( -x/ kappa_rpbe )

           f = 1.0d0 +kappa_rpbe -kappa_rpbe *ctmp1
           df_dx = ctmp1
           d2f_dx2 = -df_dx /kappa_rpbe
           d3f_dx3 = -d2f_dx2 /kappa_rpbe

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (4)                        ! W. Cohen
           s4 = s2 *s2
           ctmp1 = exp( -s2 )
           ctmp2 = coeff_WC *s4

           x = mu_gel *s2 + ( mu_wc -mu_gel )*s2 *ctmp1 + log( 1.0d0 +ctmp2 )
           dx_ds = 2.0d0 *mu_gel *s + ( mu_wc -mu_gel ) *2.0d0 *s *ctmp1 *( 1.0d0 -s2 ) &
                &  +4.0d0 *coeff_WC *s2 *s /( 1.0d0 +ctmp2 )
           d2x_ds2 = 2.0d0 *mu_gel &
                &   + 2.0d0 *( mu_wc -mu_gel ) *ctmp1 *( 1.0d0 -5.0d0 *s2 +2.0d0 *s4 ) &
                &   + 4.0d0 *coeff_WC *s2 *( 3.0d0 -ctmp2 ) /( 1.0d0 +ctmp2 )**2
           d3x_ds3 = -4.0d0 *s *( mu_wc -mu_gel ) *ctmp1 &
                &                   *( 6.0d0 -9.0d0 *s2 +2.0d0 *s4 ) &
                &    -8.0d0 *coeff_WC *s *( 12.0d0 *ctmp2 -ctmp2**2 -3.0d0 ) &
                &            / ( 1.0d0 +ctmp2 )**3
           
           ctmp3 = kappa_wc +x

           f = 1.0d0 +kappa_wc *x /ctmp3
           df_dx = kappa_wc**2 /ctmp3**2
           d2f_dx2 = -2.0d0 *df_dx /ctmp3
           d3f_dx3 = -3.0d0 *d2f_dx2 /ctmp3

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (5)                       ! HTBS
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
           df_ds = coeff_spline_htbs(1) + 2.0d0 *coeff_spline_htbs(2) *ctmp1 &
                &                       + 3.0d0 *coeff_spline_htbs(3) *ctmp2 &
                &                       + 4.0d0 *coeff_spline_htbs(4) *ctmp3 &
                &                       + 5.0d0 *coeff_spline_htbs(5) *ctmp4
           d2f_ds2 = 2.0d0 *coeff_spline_htbs(2) + 6.0d0 *coeff_spline_htbs(3) *ctmp1 &
                &                                +12.0d0 *coeff_spline_htbs(4) *ctmp2 &
                &                                +20.0d0 *coeff_spline_htbs(5) *ctmp3 
           d3f_ds3 = 6.0d0 *coeff_spline_htbs(3) +24.0d0 *coeff_spline_htbs(4) *ctmp1 &
                &                                +60.0d0 *coeff_spline_htbs(5) *ctmp2 

           df_ds   = df_ds   / (s_high_htbs -s_low_htbs )
           d2f_ds2 = d2f_ds2 / (s_high_htbs -s_low_htbs )
           d3f_ds3 = d3f_ds3 / (s_high_htbs -s_low_htbs )

        case (6)                       ! pbesol
           x = mu_pbesol *s2;   dx_ds = 2.0d0 *mu_pbesol *s;  d2x_ds2 = 2.0d0 *mu_pbesol
           d3x_ds3 = 0.0d0

           ctmp1 = kappa_pbesol + x

           f = 1.0d0 +kappa_pbesol *x /ctmp1
           df_dx = kappa_pbesol**2 /ctmp1**2
           d2f_dx2 = -2.0d0 *df_dx /ctmp1
           d3f_dx3 = -3.0d0 *d2f_dx2 /ctmp1

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (7)                       ! pbeint
           ctmp1 = mu_pbe - mu_gel
           
           alpha_pbeint = mu_gel**2 /kappa_pbeint /ctmp1
           ctmp2 = alpha_pbeint *s2
           ctmp3 = 1.0d0 +ctmp2
           
           mu_this = mu_gel +ctmp1 *ctmp2 /ctmp3
           dmu_ds =  2.0d0 *ctmp1 *alpha_pbeint *s  /ctmp3**2
           d2mu_ds2 = 2.0d0 *ctmp1 *alpha_pbeint *( 1.0d0 -3.0d0 *ctmp2 ) /ctmp3**3
           d3mu_ds3 = ctmp1 *alpha_pbeint**2 *( -24.0d0 * s )* ( 1.0d0 -ctmp2 ) &
                &           / ctmp3**4
           
           x = mu_this *s2;   dx_ds = dmu_ds *s2 +2.0d0 *mu_this *s
           d2x_ds2 = d2mu_ds2 *s2 +4.0d0 *dmu_ds *s + 2.0d0 *mu_this
           d3x_ds3 = d3mu_ds3 *s2 +6.0d0 *d2mu_ds2 *s + 6.0d0 *dmu_ds
           
           ctmp4= kappa_pbeint + x
           
           f = 1.0d0 +kappa_pbeint *x /ctmp4
           df_dx = kappa_pbeint**2 /ctmp1**2
           d2f_dx2 = -2.0d0 *df_dx /ctmp1
           d3f_dx3 = -3.0d0 *d2f_dx2 /ctmp1

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (11)                        ! B86r  I. Hamada
           x = mu_b86r *s2;   dx_ds = 2.0d0 *mu_b86r *s;    d2x_ds2 = 2.0d0 *mu_b86r
           d3x_ds3 = 0.0d0

           ctmp1 = 1.0d0 + x /kappa_b86r
           ctmp2 = ctmp1**(-4.0d0/5.0d0)

           ctmp3 = -0.8d0 /kappa_b86r *ctmp2 /ctmp1
           ctmp4 = -1.8d0 /kappa_b86r *ctmp3 /ctmp1
           ctmp5 = -2.8d0 /kappa_b86r *ctmp4 /ctmp1

           f = 1.0d0 + x *ctmp2

           df_dx = ctmp2 +ctmp3 *x
           d2f_dx2 = 2.0d0 *ctmp3 +ctmp4 *x
           d3f_dx3 = 3.0d0 *ctmp4 +ctmp5 *x

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (12)                      ! optpbe
           x = mu_optpbe *s2;   dx_ds = 2.0d0 *mu_optpbe *s;   d2x_ds2 = 2.0d0 *mu_optpbe
           d3x_ds3 = 0.0d0

           ctmp1 = kappa_optpbe + x

           f = 1.0d0 +kappa_optpbe *x /ctmp1
           df_dx = kappa_optpbe**2 /ctmp1**2
           d2f_dx2 = -2.0d0 *df_dx /ctmp1
           d3f_dx3 = -3.0d0 *d2f_dx2 /ctmp1

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (13)                         ! optb86b
           x = mu_optb86b *s2;   dx_ds = 2.0d0 *mu_optb86b *s;
           d2x_ds2 = 2.0d0 *mu_optb86b
           d3x_ds3 = 0.0d0

           ctmp1 = 1.0d0 + x /kappa_optb86b
           ctmp2 = ctmp1**(-4.0d0/5.0d0)

           ctmp3 = -0.8d0 /kappa_optb86b *ctmp2 /ctmp1
           ctmp4 = -1.8d0 /kappa_optb86b *ctmp3 /ctmp1
           ctmp5 = -2.8d0 /kappa_optb86b *ctmp4 /ctmp1

           f = 1.0d0 + x *ctmp2

           df_dx = ctmp2 +ctmp3 *x
           d2f_dx2 = 2.0d0 *ctmp3 +ctmp4 *x
           d3f_dx3 = 3.0d0 *ctmp4 +ctmp5 *x

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (14)                      ! pw86r
           s4 = s2 *s2;  s6 = s4 *s2

           x = 1.0d0 +15.d0 *a_pw86r *s2 +b_pw86r *s4 + c_pw86r *s6
           dx_ds = ( 30.0d0 *a_pw86r + 4.0d0 *b_pw86r *s2 +6.0d0 *c_pw86r *s4 ) *s
           d2x_ds2 = 30.0d0 *a_pw86r + 12.0d0 *b_pw86r *s2 +30.0d0 *c_pw86r *s4
           d3x_ds3 = ( 24.0d0 *b_pw86r +120.0d0 *c_pw86r *s2 ) *s

           f = x**(1.0d0/15.0d0)
           df_dx = f /15.0d0 /x
           d2f_dx2 = -14.0d0 /15.0d0 *df_dx /x
           d3f_dx3 = -29.0d0 /15.0d0 *d2f_dx2 /x

           df_ds = dx_ds *df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx + 3.0d0 *d2x_ds2 *dx_ds *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (15)                      ! c09x
           x = mu_c09 *s2;  dx_ds = 2.0d0 *mu_c09 *s
           d2x_ds2 = 2.0d0 *mu_c09
           d3x_ds3 = 0.0d0

           ctmp1 = alpha_c09 /mu_c09

           ctmp2 = exp( -ctmp1 *x )
           ctmp3 = exp( -ctmp1 *x /2.0d0 )

           f = 1.0d0 +x *ctmp2  + kappa_c09 *( 1.0d0 -ctmp3 )
           df_dx = ctmp2 *( 1.0d0 -ctmp1 *x ) +0.5d0 *kappa_c09 *ctmp1 *ctmp3
           d2f_dx2 = ctmp2 *ctmp1 *( ctmp1 *x -2.0d0 ) &
                &   -0.25d0 *kappa_c09 *ctmp1**2 *ctmp3
           d3f_dx3 = ctmp2 *ctmp1**2 *( -ctmp1 *x +3.0d0 ) &
                &   +0.125d0 *kappa_c09 *ctmp1**3 *ctmp3

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (16)                      ! lv-pw86r
           s4 = s2 *s2;  s6 = s4 *s2

           ctmp1 = alpha_lvpw86r *s6
           ctmp2 = 1.0d0 /(1.0d0 +ctmp1)
           ctmp3 = 1.0d0 /(beta_lvpw86r +ctmp1)

           ctmp4 = 1.0d0 +15.d0 *a_pw86r *s2 +b_pw86r *s4 + c_pw86r *s6
           ctmp5 = ( 30.0d0 *a_pw86r + 4.0d0 *b_pw86r *s2 +6.0d0 *c_pw86r *s4 ) *s
           ctmp6 = 30.0d0 *a_pw86r + 12.0d0 *b_pw86r *s2 +30.0d0 *c_pw86r *s4
           ctmp7 = ( 24.0d0 *b_pw86r +120.0d0 *c_pw86r *s2 ) *s

           term(1) = ctmp2
           term(2) = 1.d0 +mu_lvpw86r *s2
           term(3) = ctmp1 *ctmp3
           term(4) = ctmp4 **(1.0d0/15.0d0)

           dterm1(1) = -6.0d0 *alpha_lvpw86r *ctmp2**2 *s4 *s
           dterm1(2) = 2.0d0 *mu_lvpw86r *s
           dterm1(3) =  6.0d0 *alpha_lvpw86r *beta_lvpw86r *ctmp3**2 *s4 *s
           dterm1(4) = term(4) /ctmp4 *ctmp5 /15.0d0

           dterm2(1) = 6.0d0 *alpha_lvpw86r *ctmp2**3 *s4 *( 7.0d0 *ctmp1 -5.0d0 )
           dterm2(2) = 2.0d0 *mu_lvpw86r
           dterm2(3) =-6.0d0 *alpha_lvpw86r *beta_lvpw86r *ctmp3**3 *s4 &
                &          *( 7.0d0 *ctmp1 -5.0d0 *beta_lvpw86r )
           dterm2(4) = term(4) *( -14.0d0/15.0d0 /ctmp4**2 *ctmp5**2 +ctmp6 /ctmp4 ) &
                &      / 15.d0

           dterm3(1) = 6.0d0 *alpha_lvpw86r *ctmp2**4 *s2 *s &
                &            *( -56.0d0 *ctmp1**2 +140.0d0 *ctmp1 -20.0d0 )
           dterm3(2) = 0.0d0
           dterm3(3) =-6.0d0 *alpha_lvpw86r *beta_lvpw86r *ctmp3**4 *s2 *s &
                &            *( -56.0d0 *ctmp1**2 +140.0d0 *ctmp1 *beta_lvpw86r &
                &               -20.0d0 *beta_lvpw86r**2 )
           dterm3(4) = term(4) *( 14.0d0 /15.0d0 *29.0d0 /15.0d0 /ctmp4**3 *ctmp5**3 &
                &                -3.0d0 *14.0d0 /15.0d0 /ctmp4**2 *ctmp5 *ctmp6 &
                &                +ctmp7 /ctmp4 ) /15.0d0

           f = term(1) *term(2) +term(3) *term(4)
           df_ds   = dterm1(1) *term(2) +term(1) *dterm1(2) &
                &   +dterm1(3) *term(4) +term(3) *dterm1(4)
           d2f_ds2 = dterm2(1) *term(2) +2.0d0 *dterm1(1) *dterm1(2) &
                &   + term(1)  *dterm2(2) &
                &   +dterm2(3) *term(4) +2.0d0 *dterm1(3) *dterm1(4) &
                &   + term(3) *dterm2(4)
           d3f_ds3 = dterm3(1) *term(2) +3.0d0 *dterm2(1) *dterm1(2) &
                &   +3.0d0 *dterm1(1) *dterm2(2) + term(1) *dterm3(2) &
                &   +dterm3(3) *term(4) +3.0d0 *dterm2(3) *dterm1(4) &
                &   +3.0d0 *dterm1(3) *dterm2(4) + term(3) *dterm3(4)

        case (20)            !  Engel Vosko 93
           s4 = s2 *s2;  s6 = s4 *s2

           x = 1.0D0 +a1_ev93 *s2 +a2_ev93 *s4 +a3_ev93 *s6
           y = 1.0D0 +b1_ev93 *s2 +b2_ev93 *s4 +b3_ev93 *s6

           dx_ds = s *( 2.0D0 *a1_ev93 +4.0D0 *a2_ev93 *s2 +6.0D0 *a3_ev93 *s4 )
           dy_ds = s *( 2.0D0 *b1_ev93 +4.0D0 *b2_ev93 *s2 +6.0D0 *b3_ev93 *s4 )

           d2x_ds2 = 2.0D0 *a1_ev93 +12.0D0 *a2_ev93 *s2 +30.0D0 *a3_ev93 *s4
           d2y_ds2 = 2.0D0 *b1_ev93 +12.0D0 *b2_ev93 *s2 +30.0D0 *b3_ev93 *s4

           d3x_ds3 = s *( 24.0D0 *a2_ev93 +120.0D0 *a3_ev93 *s2 )
           d3y_ds3 = s *( 24.0D0 *b2_ev93 +120.0D0 *b3_ev93 *s2 )

           ctmp1 = y*dx_ds -x*dy_ds
           ctmp2 = y*d2x_ds2 -x*d2y_ds2
           ctmp3 = y*d3x_ds3 -x*d3y_ds3
           ctmp4 = dy_ds *d2x_ds2 -dx_ds *d2y_ds2

           y2 = y**2;    y3 = y2 *y;     y4 = y3 *y

           f = x /y
           df_ds = ctmp1 /y2
           d2f_ds2 = -2.0D0 *ctmp1 *dy_ds /y3 +ctmp2 /y2
           d3f_ds3 = ( 6.0D0 *dy_ds**2 /y4 -2.0D0 *d2y_ds2 /y3 ) *ctmp1 &
                &   -4.0D0 *dy_ds *ctmp2 /y3 &
                &   +( ctmp3 +ctmp4 ) /y2

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
           call phase_error_with_msg(6,"xc functional not supported (PAW)",__LINE__,__FILE__)
        end select

        ex  = fac *f *d
        exd = thrd4 *fac *( f - s *df_ds )
        exdd = ax *df_ds *0.5d0 /thpith
        !------------------------------------------     
        exc0 = ex / facw
        excd = exd
        dFx_drho(i, is) = excd 
        if ( pot_type == 25 ) dFx_drho(i,is) = dFx_drho(i,is) +pot_add

        excdd = exdd
        ! gradient of charge density 95/12/3 H.S.
        if(dabs(grad_rho(i, is)) > 1.d-9) then
           dFx_dgradrho(i, is) = excdd
           !           dFx_dgradrho(i, is) = excdd / grad_rho(i, is)           ! 2010 5
        else
           dFx_dgradrho(i, is) = 0.d0
        endif
        ! gradient of charge density 95/12/3 H.S.
        !        exc1 = exc1 + exc0*wos(i)
        exc(i)=exc(i)+exc0

        ! -----------------------------------
        ivd =1.d0 /d
        ivdthrd2 =ivd**thrd2

        if ( nspin == 1 ) then
           dFx_drr(i,is)   = ninth4 *ax *ivdthrd2*( f -s *df_ds +4.d0 *d2f_ds2 *s2  )
           dFx_drg(i,is)   =-thrd2 *ax /thpith *d2f_ds2 *s *ivd
           dFx_dgg(i,is)   = 0.25d0 *ax /thpith /thpith *d2f_ds2 *ivdthrd2**2
           dFx_drrr(i,is)  =-twtysvnth8 *ax *ivdthrd2 *ivd &
                &           *( f -df_ds *s +18.d0 *d2f_ds2 *s2 +8.d0 *d3f_ds3 *s2 *s )
           dFx_drrg(i,is)  = ninth2 *ax *ivd *ivd /thpith &
                &           *( 7.d0 *d2f_ds2 *s +4.d0 *d3f_ds3 *s2 )
           dFx_drgg(i,is)  =-thrd *ax /thpith /thpith *ivd *ivdthrd2**2 &
                &           *( d2f_ds2 +d3f_ds3 *s )
           dFx_dggg(i,is)  =eghth *ax /(thpith**3) *d3f_ds3 *ivd**2 *ivdthrd2
        else if ( nspin == 2 ) then
           dFx_drr(i,is)   = ninth8 *ax *ivdthrd2 *( f -s *df_ds +4.d0 *d2f_ds2 *s2 )
           dFx_drg(i,is)   =-thrd4 *ax /thpith *d2f_ds2 *s *ivd
           dFx_dgg(i,is)   = 0.5d0 *ax /thpith /thpith *d2f_ds2 *ivdthrd2**2
           dFx_drrr(i,is)  =-twtysvnth32 *ax *ivdthrd2 *ivd &
                &           *( f -df_ds*s +18.d0 *d2f_ds2 *s2 +8.d0 *d3f_ds3 *s2 *s )
           dFx_drrg(i,is)  = ninth8 *ax *ivd *ivd /thpith &
                &           *( 7.d0 *d2f_ds2 *s +4.d0 *d3f_ds3 *s2 )
           dFx_drgg(i,is)  =-thrd4 *ax/(thpith**2) *ivd *ivdthrd2**2 &
                &           *( d2f_ds2 +d3f_ds3 *s )
           dFx_dggg(i,is)  = 0.5d0 *ax/(thpith**3) *d3f_ds3 *ivd**2 *ivdthrd2
        end if

     end do
     !     exc(i) = exc(i) + exc1
  end do
end subroutine ex_gga_paw_library

subroutine cr_gga_paw_library( nrc, dnr, nspin, chgrhr_l, grad_trho, exc, &
     &                         dF_drho,  dF_dgradrho, dFc_daa,  dFc_dbb,  dFc_dgg ,&
     &                         dFc_dab,  dFc_dag,     dFc_dbg,  dFc_daaa, dFc_dbbb, &
     &                         dFc_dggg, dFc_daab,    dFc_daag, dFc_dabb, dFc_dbbg, &
     &                         dFc_dagg, dFc_dbgg, dFc_dabg, pot_type )
  use m_Const_Parameters,  only : PAI,DP
  implicit none

  integer,intent(in)        :: nrc, dnr, nspin, pot_type
  real(kind=DP),intent(in)  :: chgrhr_l(nrc,nspin)
  real(kind=DP),intent(in) :: grad_trho(nrc)
  real(kind=DP),intent(inout) :: exc(nrc)
  real(kind=DP),intent(out) :: dF_drho(nrc,nspin)
  real(kind=DP),intent(out) :: dF_dgradrho(nrc)
  real(kind=DP),intent(out) :: dFc_daa(nrc)
  real(kind=DP),intent(out) :: dFc_dbb(nrc)
  real(kind=DP),intent(out) :: dFc_dgg(nrc)
  real(kind=DP),intent(out) :: dFc_dab(nrc)
  real(kind=DP),intent(out) :: dFc_dag(nrc)
  real(kind=DP),intent(out) :: dFc_dbg(nrc)
  real(kind=DP),intent(out) :: dFc_daaa(nrc)
  real(kind=DP),intent(out) :: dFc_dbbb(nrc)
  real(kind=DP),intent(out) :: dFc_dggg(nrc)
  real(kind=DP),intent(out) :: dFc_daab(nrc)
  real(kind=DP),intent(out) :: dFc_daag(nrc)
  real(kind=DP),intent(out) :: dFc_dabb(nrc)
  real(kind=DP),intent(out) :: dFc_dbbg(nrc)
  real(kind=DP),intent(out) :: dFc_dagg(nrc)
  real(kind=DP),intent(out) :: dFc_dbgg(nrc)
  real(kind=DP),intent(out) :: dFc_dabg(nrc)

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

  ! vwn 95/12/2 Y.M
  real(kind=DP), parameter  :: aq  = 0.0168869d0,    a1q  = 0.11125d0

  real(kind=DP), parameter  :: b1q = 10.3570d0,     b2nq = 3.6231d0
  real(kind=DP), parameter  :: b3q = 0.88026d0,     b4q  = 0.49671d0

  ! vwn 95/12/2 Y.M
  real(kind=DP), parameter  :: gamma = 0.03109069086965489503494086371273d0

  real(kind=DP), parameter  :: thrd = 0.33333333333333330,&
       &                       sixth7 = 1.1666666666666666d0
  real(kind=DP), parameter  :: thrd7= 2.333333333333333333d0
  real(kind=DP), parameter  :: sixth13= 2.166666666666666666d0
  real(kind=DP), parameter  :: sixth19= 3.166666666666666666d0
  real(kind=DP), parameter  :: ninth= 0.11111111111111111111d0

  real(kind=DP), parameter  :: twntysvnth=  3.703703703703703d-2
  real(kind=DP), parameter  :: fftyfrth= 1.851851851851852d-2
  real(kind=DP), parameter  :: eityonth= 1.234567901234568d-2

! ----- parameters
  real(kind=DP), parameter :: beta_pbe   = 0.06672455060314922d0
  real(kind=DP), parameter :: beta_pbesol = 0.046d0
  real(kind=DP), parameter :: beta_pbeint = 0.052d0
!--------------------------

  integer       :: is,i
  real(kind=DP) :: facw, exc1,g,g3,facpon,d,dd,rs,fk,sk,t,s,q0 &
       &         , rs12,rs32,rsp,q1,q2,eu,q3,eurs,ecd,pon &
       &         , b,b2,t2,t4,t6,q4,q5,h0,h,q8,h0t,h0b,h0rs,ht,hrs &
       &         , ec1,ec1d,ec1dd,exc0,excd,excdd &
       &         , thrd2,thrd4,fzdd,zeta,onpzeta,onmzeta,fz,zetadxd &
       &         , q0p,q1p,q2p,q3p,q0q,q1q,q2q,q3q,fzd,onzt4,gzt,bzt,h0zt,euzt &
       &         , zetadxda,zetadxdb,ecda,ecdb,ec1da,ec1db,excda,excdb

  real(kind=DP) :: bet, delt

  real(kind=DP) :: rsdn,rsdnn,rsdnnn
  real(kind=DP) :: q1drr,q1drrr
  real(kind=DP) :: eudrr,eudrrr,invq12pq1,q12p1
  real(kind=DP) :: eudn,eudnn,eudnnn,ecdd,ecddd
  real(kind=DP) :: dtdn,dtdg,dtdnn,dtdgg,dtdng,dtdnnn,dtdnng,dtdngg,dtdggg
  real(kind=DP) :: dAde,dAdee,dAdeee
  real(kind=DP) :: dAdn,dAdnn,dAdnnn
  real(kind=DP) :: fE,fEdt,fEdtt,fEda,fEdaa
  real(kind=DP) :: fF,fFdt,fFdtt,fFda,fFdaa
  real(kind=DP) :: fG,fGdt,fGdtt,fGda,fGdaa
  real(kind=DP) :: at2,a2t4,a3t6,t3,t5,t7,t8
  real(kind=DP) :: dHda,dHdt,dHdaa,dHdat,dHdtt,dHdaaa,dHdaat,dHdatt,dHdttt,betivE,betivE2
  real(kind=DP) :: dHdn,dHdg,dHdnn,dHdng,dHdgg,dHdnnn,dHdnng,dHdngg,dHdggg
  real(kind=DP) :: zetada,zetadb,zetadaa,zetadab,zetadbb,zetadaaa,zetadaab,zetadabb,zetadbbb
  real(kind=DP) :: id2,id3,id4
  real(kind=DP) :: e1dr,e1drr,e1drrr
  real(kind=DP) :: e1pdr,e1pdrr,e1pdrrr
  real(kind=DP) :: e1qdr,e1qdrr,e1qdrrr
  real(kind=DP) :: q1pdrr,q1pdrrr,invq1p2pq1p,q1p2p1
  real(kind=DP) :: q1qdrr,q1qdrrr,invq1q2pq1q,q1q2p1
  real(kind=DP) :: fzdd0,fzddd,opzthrdm2,omzthrdm2
  real(kind=DP) :: zeta2,zeta3,zeta4,eudrz,eudzz,eudrrz,eudrzz,eudzzz
  real(kind=DP) :: e1,e2,e3
  real(kind=DP) :: euda,eudb,eudaa,eudbb,eudab,eudaaa,eudaab,eudabb,eudbbb
  real(kind=DP) :: rsdn2,rsdn3,zetada2,zetada3,zetadb2,zetadb3
  real(kind=DP) :: gzz,gzzz,duda,dudb,dudaa,dudab,dudbb,dudaaa,dudaab,dudabb,dudbbb
  real(kind=DP) :: dtdu,dtduu,dtduuu,dtdnu,dtdnnu,dtdnuu,dtdug,dtduug,dtdugg,dtdnug,invg
  real(kind=DP) :: dtda,dtdb,dtdaa,dtdab,dtdag,dtdbg,dtdbb
  real(kind=DP) :: dtdaaa,dtdaab,dtdabb,dtdbbb,dtdaag,dtdagg,dtdbbg,dtdbgg,dtdabg
  real(kind=DP) :: invg3,m3u
  real(kind=DP) :: dAdu,dAduu,dAduuu,dAdeu,dAdeeu,dAdeuu
  real(kind=DP) :: dAda,dAdb,dAdaa,dAdab,dAdbb,dAdaaa,dAdaab,dAdabb,dAdbbb
  real(kind=DP) :: dHdf,dHdff,dHdft,dHdfff,dHdfft,dHdftt
  real(kind=DP) :: dHdu,dHduu,dHduuu,dHdfu,dHdtu,dHdffu,dHdfuu,dHdttu,dHdtuu,dHdftu
  real(kind=DP) :: tiu,siu2
  real(kind=DP) :: dHdb,dHdbb,dHdag,dHdbg,dHdab,dHdbbb,dHdaab,dHdabb,dHdbbg,dHdbgg,dHdaag,dHdagg,dHdabg

  real(kind=DP), parameter  :: eta = 1.d-12
  real(kind=DP), parameter  :: sixthm = -0.1666666666666666666666d0

! -----
  select case ( pot_type )
  case (6)                           ! pbesol
     bet = beta_pbesol
  case (7)                           ! pbeint
     bet = beta_pbeint
  case default
     bet = beta_pbe
  end select

  delt = bet / gamma
! -----

  facw = nspin

  if ( nspin == 1 ) then
     g   = 1.d0
     g3  = g**3
!!$        g4  = g3*g
     facpon = -delt/(g3*bet)

     do i = 1,nrc,dnr ! MPI
        d  = facw*chgrhr_l(i, 1)
        if(d < density_minimum) cycle
        ! gradient of charge density 95/12/3 H.S.
        dd = facw* grad_trho(i)
        rs = (0.75d0/(PAI*d))**thrd
        fk = (3.d0*PAI*PAI*d)**thrd
        sk = dsqrt(4.d0*fk/PAI)

        if(d > density_minimum2) then
           t = dd/(d*sk*2.d0)
           s = dd/(d*fk*2.d0)
        else
           t = 0.d0
           s = 0.d0
        endif

        q0 = -2.d0 *a *( 1.d0 +a1*rs )
        rs12 = dsqrt(rs)
        rs32 = rs12**3
        rsp = rs**p
        q1 = 2.d0 *a*( b1*rs12 +b2n*rs +b3*rs32 +b4*rs*rsp )
        q2 = log( 1.d0 +1.d0/q1 )
        eu = q0*q2
        q3 = a*( b1/rs12 +2.d0*b2n +3.d0 *b3*rs12 +2.d0*b4*p1*rsp )
        eurs = -2.d0*a*a1*q2 - q0*q3/(q1**2 + q1)
        ecd = eu-thrd*rs*eurs
        pon = facpon*eu
        b   = delt/(dexp(pon)-1.d0)
        b2  = b*b
        t2  = t*t
!!$           s2  = s*s
        t4  = t2*t2
        t6  = t4*t2
        q4  = 1.d0 + b*t2
        q5  = 1.d0 + b*t2 +b2*t4
        h0  = g3*(bet/delt)*dlog( 1.d0 +delt*q4*t2/q5 )
        h   = h0
        q8  = q5*q5 +delt*q4*q5*t2
        h0t = 2.d0*bet*t*( 1.d0 +2.d0*b*t2 )/q8 * g3
        h0b = -bet*t6*( 2.d0*b +b2*t2 )/q8 * g3
        h0rs = h0b*b*eurs*(b+delt)/ bet / g3
        ht   = h0t
        hrs  = h0rs
        ec1  = d*h/facw
        ec1d = h -thrd*rs*hrs -sixth7*t*ht
        ec1dd = 0.5d0*ht/sk

        exc0 = eu*d/facw + ec1
        excd = ecd + ec1d
        dF_drho(i, 1) = dF_drho(i, 1) + excd
        excdd = ec1dd
        ! gradient of charge density 95/12/3 H.S.
        if(dabs(grad_trho(i)) > 1.0d-9) then
           !              dF_dgradrho(i) = excdd / grad_trho(i)
           dF_dgradrho(i) = excdd
        else
           dF_dgradrho(i) = 0.d0
        endif
        ! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d/facw
        !           exc1 = exc1 + exc0*wos(i)
        exc(i)=exc(i)+exc0

        ! ***** Calc heigher derivatives. *****

        rsdn=-thrd*rs/d
        rsdnn=-4.d0*thrd*rsdn/d
        rsdnnn=-thrd7*rsdnn/d
        !print '(i2,5e19.10)',i,d,rs,rsdn,rsdnn,rsdnnn
        !q3 = a*(b1/rs12+2*b2n+3*b3*rs12+2*b4*p1*rsp)   
        q1drr =0.5d0*a*( -b1/rs32+3.d0*b3/rs12+8.d0*b4 )
        q1drrr=0.75d0*a*( b1/rs - b3 )/rs32
        !eurs = -2*a*a1*q2 - q0*q3/(q1**2 + q1)
        invq12pq1=1.d0/(q1**2 + q1)
        q12p1=2.d0*q1+1.d0
        eudrr=4.d0*a*a1*q3*invq12pq1-q0*( q1drr-q12p1*q3**2*invq12pq1 )*invq12pq1
        eudrrr=6.d0*a*a1*( q1drr - q12p1*q3**2*invq12pq1 ) * invq12pq1 &
             -q0*( q1drrr -(2.d0*q3**3+3.d0*q12p1*q3*q1drr -  &
             2.d0*q12p1**2*q3**3* invq12pq1)*invq12pq1 ) * invq12pq1
        eudn=eurs*rsdn
        eudnn=eudrr*rsdn**2+eurs*rsdnn
        eudnnn=eudrrr*rsdn**3+3.d0*eudrr*rsdn*rsdnn+eurs*rsdnnn
        !print '(i2,5e19.12)',i,d,eu,eudn,eudnn,eudnnn
        ecdd=2.d0*eudn+d*eudnn
        ecddd=3.d0*eudnn+d*eudnnn
        !print '(i2,5e19.12)',i,d,eu*d,ecdd,ecddd

        dtdn  = -sixth7*t/d
        dtdnn = -sixth13*dtdn/d
        dtdnnn= -sixth19*dtdnn/d
        !print '(i2,5e19.12)',i,d,t,dtdn,dtdnn,dtdnnn
        if(dabs(grad_trho(i)) > 1.0d-9) then
           dtdg  = t/dd
           dtdng = -sixth7*dtdg/d
           dtdnng= -sixth13*dtdng/d
        else
           dtdg  = 0.d0
           dtdng = 0.d0
           dtdnng= 0.d0
        endif
        dtdngg = 0.d0
        dtdggg = 0.d0

        dAde   = b * ( b + delt ) / bet
        dAdee  = 2.d0 * dAde * ( b + 0.5d0 * delt ) / bet
        dAdeee = 2.d0 * ( dAde**2 + ( b + 0.5d0 * delt )*dAdee ) /bet

        dAdn   = dAde  *eudn
        dAdnn  = dAdee *eudn**2  + dAde*eudnn
        dAdnnn = dAdeee*eudn**3  + 3.d0*dAdee*eudn*eudnn + dAde*eudnnn   
        !print '(i2,5e19.10)' ,i,d,b,dAdnn

        at2  = b * t2
        a2t4 = at2**2
        a3t6 = at2*a2t4
        t3   = t2*t
        t5   = t4*t
        t7   = t6*t
        t8   = t7*t
        fE    = q5**2 + delt*t2*q4*q5
        fEdt  = 4.d0*b*t*q5*( 1.d0 + 2.d0*at2 ) + &
             delt *2.d0*t*(1.d0 + 4.d0*at2 + 6.d0*a2t4 + 4.d0*a3t6 )
        fEdtt = 4.d0*b*( 1.d0 + 9.d0*at2 + 15.d0*a2t4 + 14.d0*a3t6 ) + &
             delt*( 2.d0 + 24.d0*at2 + 60.d0*a2t4 + 56.d0*a3t6 )
        fEda  = 2.d0*t2*q5*( 1.d0 + 2.d0*at2 ) + &
             delt*t4*( 2.d0 + 4.d0*at2 + 3.d0*a2t4 )
        fEdaa = 6.d0*t4*( 2.d0*q5 - 1.d0) + &
             delt*2.d0*t6*( 2.d0 + 3.d0*at2 )
        fF    = t*( 1.d0 + 2.d0*at2 )         
        fFdt  = 1.d0 + 6.d0*at2
        fFdtt = 12.d0*b*t
        fFda  = 2.d0*t3
        fFdaa = 0.d0
        fG    = 2.d0*b*t6 + b2*t8
        fGdt  = 12.d0*b*t5 + 8.d0*b2*t7
        fGdtt = 60.d0*b*t4 + 56.d0*b2*t6
        fGda  = 2.d0*t6 + 2.d0*b*t8
        fGdaa = 2.d0*t8

        betivE=bet/fE
        betivE2=betivE/fE
        dHda   =     -fG*betivE
        dHdt   = 2.d0*fF*betivE
        dHdaa  =     -( fGda*fE-fG*fEda )*betivE2
        dHdat  =     -( fGdt*fE-fG*fEdt )*betivE2
        dHdtt  = 2.d0*( fFdt*fE-fF*fEdt )*betivE2
        dHdaaa =     -( fGdaa*fE-fG*fEdaa-2.d0*( fGda*fE-fG*fEda )*fEda/fE )*betivE2
        dHdaat = 2.d0*( fFdaa*fE-fF*fEdaa-2.d0*( fFda*fE-fF*fEda )*fEda/fE )*betivE2
        dHdatt =     -( fGdtt*fE-fG*fEdtt-2.d0*( fGdt*fE-fG*fEdt )*fEdt/fE )*betivE2
        dHdttt = 2.d0*( fFdtt*fE-fF*fEdtt-2.d0*( fFdt*fE-fF*fEdt )*fEdt/fE )*betivE2

        dHdn     = dHda  *dAdn       + dHdt*dtdn
        dHdg     = dHdt  *dtdg
        dHdnn    = dHdaa *dAdn**2    + 2.d0*dHdat*dAdn*dtdn  + dHdtt*dtdn**2   + dHda*dAdnn  + dHdt*dtdnn
        dHdng    = dHdat *dAdn*dtdg  + dHdtt*dtdn*dtdg       + dHdt*dtdng
        dHdgg    = dHdtt*dtdg**2
        dHdnnn   = dHdaaa*dAdn**3    + dHdttt*dtdn**3 &
             + 3.d0*dHdaat*dAdn**2*dtdn  + 3.d0*dHdatt*dAdn*dtdn**2 &
             + 3.d0*dHdaa*dAdn*dAdnn + 3.d0*dHdtt*dtdn*dtdnn + 3.d0*dHdat*(dAdnn*dtdn + dAdn*dtdnn) &
             + dHda*dAdnnn   + dHdt*dtdnnn
        dHdnng   = dHdaat*dAdn**2*dtdg   + 2.d0*dHdatt*dAdn*dtdn*dtdg    + dHdttt*dtdn**2*dtdg &
             + dHdat*(dAdnn*dtdg + 2.d0*dAdn*dtdng)  + dHdtt*(dtdnn*dtdg + 2.d0*dtdn*dtdng) &
             + dHdt*dtdnng
        dHdngg   = dHdatt*dAdn*dtdg**2   + dHdttt*dtdn*dtdg**2   + 2.d0*dHdtt*dtdg*dtdng
        dHdggg   = dHdttt*dtdg**3
        !dF_drho(i, 1)=h+d*dHdn+ecd
        dFc_daa(i)   = 2.d0*dHdn + d*dHdnn       + ecdd
        dFc_dag(i)   = dHdg  + d*dHdng
        dFc_dgg(i)   = d*dHdgg
        dFc_daaa(i)  = 3.d0*dHdnn    + d*dHdnnn  + ecddd
        dFc_daag(i)  = 2.d0*dHdng    + d*dHdnng
        dFc_dagg(i)  = dHdgg         + d*dHdngg
        dFc_dggg(i)  = d*dHdggg  
     end do

  else if ( nspin ==  2 ) then
     thrd2 = thrd * 2.d0
     thrd4 = thrd * 4.d0
     fzdd = 8.d0 /(9.d0*(2.d0**thrd4 - 2.d0 ))

     do i = 1,nrc,dnr
        d = chgrhr_l(i, 1) + chgrhr_l(i, nspin)
        if(d < density_minimum .or. chgrhr_l(i,1) < density_minimum .or. chgrhr_l(i,nspin)<density_minimum) cycle

        zeta = ( chgrhr_l(i, 1) - chgrhr_l(i, 2) ) / d
        onpzeta = 2.d0 *chgrhr_l(i, 1) / d
        onmzeta = 2.d0 *chgrhr_l(i, 2) / d
        g = ( onpzeta**thrd2 + onmzeta**thrd2 )/2
        if(dabs(zeta) < zeta_minimum ) then
           fz= 4.d0/9 *zeta*zeta*(1.d0 + 5.d0/54.d0 *zeta*zeta &
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
           ! ****** aas mkatsu *****
           !              t = dd/(d*sk*2)
           t = dd/(d*sk*2.d0)/g         !    mkatsu aas
           s = dd/(d*fk*2.d0)
        else
           t = 0.d0
           s = 0.d0
        endif

        q0 = -2.d0*a*( 1.d0 +a1*rs )
        rs12 = dsqrt(rs)
        rs32 = rs12**3
        rsp = rs**p
        q1 = 2.d0 *a*( b1*rs12 +b2n*rs +b3*rs32 +b4*rs*rsp )
        q2 = log( 1.d0 +1.d0/q1 )
        q3 = a*( b1/rs12 +2.d0*b2n +3.d0*b3*rs12 +2.d0*b4*p1*rsp )

        q0p = -2.d0 *ap *( 1.d0 +a1p*rs )
        q1p = 2.d0 *ap*( b1p*rs12 +b2np*rs+b3p*rs32+b4p*rs*rsp )
        q2p = log( 1.d0 +1.d0/q1p )
        q3p = ap*( b1p/rs12 +2.d0*b2np +3.d0*b3p*rs12 +2.d0*b4p*p1*rsp )
        ! vwn 95/12/2 Y.M
        q0q = -2.d0 *aq*( 1.d0 +a1q*rs )
        q1q =  2.d0 *aq*( b1q*rs12 +b2nq*rs +b3q*rs32 +b4q*rs*rsp )
        q2q = log( 1.d0 +1.d0/q1q )
        q3q = aq*( b1q/rs12 +2.d0*b2nq +3.d0*b3q*rs12 +2.d0*b4q*p1*rsp )

        if(dabs(zeta) < zeta_minimum) then
           fzd = 8.d0/9.d0 /( 2.d0**thrd4 -2.d0 ) &
                &     *zeta *( 1.d0 + 5.d0/27.d0 *zeta*zeta &
                &           *( 1.d0 +22.d0/45.d0 *zeta*zeta ))
        else
           fzd = 4.d0/3.d0 /( 2.d0**thrd4 -2.d0 ) &
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

        !           zetadxd = - 2 * chgrhr_l(i, 3-is) * (-1.d0)**is / d
        zetadxda = - 2.d0 * chgrhr_l(i, 2) * (-1.d0) / d
        zetadxdb = - 2.d0 * chgrhr_l(i, 1) / d
        !           ecd = eu-thrd*rs*eurs + euzt * zetadxd
        ecda = eu-thrd*rs*eurs + euzt * zetadxda 
        ecdb = eu-thrd*rs*eurs + euzt * zetadxdb

        pon = facpon*eu
        b   = delt/(dexp(pon)-1.d0)
        b2  = b*b
        t2  = t*t
!!$           s2  = s*s
        t4  = t2*t2
        t6  = t4*t2
        q4  = 1.d0 +b*t2
        q5  = 1.d0 +b*t2+b2*t4
        h0  = g3*(bet/delt)*dlog( 1.d0 +delt*q4*t2/q5 )
        h   = h0

        q8  = q5*q5+delt*q4*q5*t2
        h0t = 2.d0 *bet*t*( 1.d0 +2.d0 *b*t2 )/q8 * g3
        h0b = -bet*t6*( 2.d0 *b+b2*t2 )/q8 * g3
        h0rs = h0b * b * eurs * (b+delt) / bet / g3
        ht   = h0t
        hrs  = h0rs

        if(dabs(zeta) < zeta_minimum) then
           gzt = -2.d0/9*zeta *( 1.d0 +14.d0/27 *zeta*zeta &
                &             *( 1.d0 +13.d0/18 *zeta*zeta ))
        else
           gzt = ( (onpzeta**2+eta)**sixthm - (onmzeta**2+eta)**sixthm)/3.d0
        end if

        bzt = b * ( b+delt) * ( euzt - 3*eu*gzt/g ) / ( bet * g3 ) 
        h0zt = 3.d0 * gzt * h0 / g + h0b * bzt

        ec1  = d*h * 0.5d0

        ! --> revised by Kato-san,  1999Jul16
!!$           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + h0zt * zetadxd
        !           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxd
        ec1da = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxda
        ec1db = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxdb
        ec1dd = 0.5d0*ht/sk / g

10      continue


        exc0 = eu*d * 0.5d0 + ec1
        !           excd = ecd + ec1d
        excda= ecda + ec1da
        excdb= ecdb + ec1db
        !           dF_drho(i, is) = dF_drho(i, is) + excd
        dF_drho(i, 1) = dF_drho(i, 1) + excda
        dF_drho(i, 2) = dF_drho(i, 2) + excdb
        excdd = ec1dd
        ! gradient of charge density 95/12/3 H.S.
        !           if(is ==  2) then
        if(dabs(grad_trho(i)) > 1.0d-9) then
           !                 dF_dgradrho(i) = excdd / grad_trho(i)
           dF_dgradrho(i) = excdd
        else
           dF_dgradrho(i) = 0.d0
        endif
        !           endif
        ! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d * 0.5
        !           exc1 = exc1 + exc0*wos(i)
        exc(i)=exc(i)+exc0*2.d0

        ! ***** Calc heigher order derivatives. *****

        rsdn=-thrd*rs/d
        rsdnn=-4.d0*thrd*rsdn/d
        rsdnnn=-thrd7*rsdnn/d

        id2=1.d0/d/d
        id3=id2/d
        id4=id3/d
        zetada   =     zetadxda/d
        zetadb   =     zetadxdb/d
        zetadaa  =    -4.d0*chgrhr_l(i, 2)*id3
        zetadab  =     2.d0*zeta*id2
        zetadbb  =     4.d0*chgrhr_l(i, 1)*id3
        zetadaaa =    12.d0*chgrhr_l(i, 2)*id4
        zetadaab =   ( -4.d0*chgrhr_l(i, 1) + 8.d0*chgrhr_l(i, 2))*id4
        zetadabb =   ( -8.d0*chgrhr_l(i, 1) + 4.d0*chgrhr_l(i, 2))*id4
        zetadbbb =   -12.d0*chgrhr_l(i, 1)*id4

        !print '(i2,7d19.10)', i,chgrhr_l(i,1:2),zetadab,zetadaaa,zetadaab,zetadabb,zetadbbb

        q1drr =0.5d0*a*( -b1/rs32+3.d0*b3/rs12+8.d0*b4 )
        q1drrr=0.75d0*a*( b1/rs - b3 )/rs32
        invq12pq1=1.d0/(q1**2 + q1)
        q12p1=2.d0*q1+1.d0
        e1dr=-2.d0*a*a1*q2-q0*q3*invq12pq1
        e1drr=4.d0*a*a1*q3*invq12pq1-q0*( q1drr-q12p1*q3**2*invq12pq1 )*invq12pq1
        e1drrr=6.d0*a*a1*( q1drr - q12p1*q3**2*invq12pq1 ) * invq12pq1 &
             -q0*( q1drrr -(2.d0*q3**3+3.d0*q12p1*q3*q1drr -  &
             2.d0*q12p1**2*q3**3* invq12pq1)*invq12pq1 ) * invq12pq1

        q1pdrr =0.5d0*ap*( -b1p/rs32+3.d0*b3p/rs12+8.d0*b4p )
        q1pdrrr=0.75d0*ap*( b1p/rs - b3p )/rs32
        invq1p2pq1p=1.d0/(q1p**2 + q1p)
        q1p2p1=2.d0*q1p+1.d0
        e1pdr=-2.d0*ap*a1p*q2p-q0p*q3p*invq1p2pq1p
        e1pdrr=4.d0*ap*a1p*q3p*invq1p2pq1p-q0p*( q1pdrr-q1p2p1*q3p**2*invq1p2pq1p )*invq1p2pq1p
        e1pdrrr=6.d0*ap*a1p*( q1pdrr - q1p2p1*q3p**2*invq1p2pq1p ) * invq1p2pq1p &
             -q0p*( q1pdrrr -(2.d0*q3p**3+3.d0*q1p2p1*q3p*q1pdrr -  &
             2.d0*q1p2p1**2*q3p**3* invq1p2pq1p)*invq1p2pq1p ) * invq1p2pq1p

        q1qdrr =0.5d0*aq*( -b1q/rs32+3.d0*b3q/rs12+8.d0*b4q )
        q1qdrrr=0.75d0*aq*( b1q/rs - b3q )/rs32
        invq1q2pq1q=1.d0/(q1q**2 + q1q)
        q1q2p1=2.d0*q1q+1.d0
        e1qdr=-2.d0*aq*a1q*q2q-q0q*q3q*invq1q2pq1q
        e1qdrr=4.d0*aq*a1q*q3q*invq1q2pq1q-q0q*( q1qdrr-q1q2p1*q3q**2*invq1q2pq1q )*invq1q2pq1q
        e1qdrrr=6.d0*aq*a1q*( q1qdrr - q1q2p1*q3q**2*invq1q2pq1q ) * invq1q2pq1q &
             -q0q*( q1qdrrr -(2.d0*q3q**3+3.d0*q1q2p1*q3q*q1qdrr -  &
             2.d0*q1q2p1**2*q3q**3* invq1q2pq1q)*invq1q2pq1q ) * invq1q2pq1q

        zeta2=zeta**2
        zeta3=zeta2*zeta
        zeta4=zeta2*zeta2

        if(dabs(zeta) < zeta_minimum) then
           fzdd0 = fzdd &
                &            *( 1 + 5.d0*ninth *zeta2 &
                &            *( 1 +22.d0*twntysvnth *zeta2 ))
           fzddd = fzdd*10.d0*ninth*zeta*( 1.d0 + 44.d0*twntysvnth*zeta2 )
           gzz   = -2.d0*ninth*( 1.d0 + 14.d0*ninth*zeta2*( 1.d0 + 65.d0*fftyfrth*zeta2 ) )
           gzzz  = -56.d0*eityonth*zeta*( 1.d0 + 65.d0*twntysvnth*zeta2 )
        else
           opzthrdm2=1.d0/onpzeta**(thrd2)
           omzthrdm2=1.d0/onmzeta**(thrd2)
           fzdd0 =  fzdd*0.5d0*( opzthrdm2 + omzthrdm2 )
           fzddd = -fzdd*thrd*( opzthrdm2/onpzeta - omzthrdm2/onmzeta )
           gzz   = -ninth*( opzthrdm2**2 + omzthrdm2**2 )
           !              gzz   = -ninth*( 1.d0/onpzeta**(thrd4) + 1.d0/onmzeta**(thrd4) )
           gzzz  = 4.d0*twntysvnth*( opzthrdm2**2/onpzeta - omzthrdm2**2/onmzeta )
        end if

        e1=q0*q2
        e2=q0p*q2p
        e3=q0q*q2q

        eudrr    = e1drr  + ( e1pdrr  -  e1drr  ) * fz * zeta4 - e1qdrr  * fz * onzt4 / fzdd
        eudrrr   = e1drrr + ( e1pdrrr -  e1drrr ) * fz * zeta4 - e1qdrrr * fz * onzt4 / fzdd
        eudzz    = ( e2 - e1 )*( fzdd0*zeta4 + 8.d0*fzd*zeta3 + 12.d0*fz*zeta2 ) &
             - e3*( fzdd0*onzt4 - 8.d0*fzd*zeta3 - 12.d0*fz*zeta2 )/fzdd
        eudzzz   = ( e2 - e1 )*( fzddd*zeta4 + 12.d0*fzdd0*zeta3 + 36.d0*fzd*zeta2 + 24.d0*fz*zeta ) &
             - e3*( fzddd*onzt4 - 12.d0*fzdd0*zeta3 - 36.d0*fzd*zeta2 - 24.d0*fz*zeta )/fzdd
        eudrz    = ( e1pdr - e1dr )*( fzd*zeta4 + 4.d0*fz*zeta3 ) &
             - e1qdr *( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd 
        eudrrz   = ( e1pdrr  - e1drr  )*( fzd*zeta4 + 4.d0*fz*zeta3 ) &
             - e1qdrr *( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd
        eudrzz   = ( e1pdr - e1dr )*( fzdd0*zeta4 + 8.d0*fzd*zeta3 + 12.d0*fz*zeta2 ) &
             - e1qdr *( fzdd0*onzt4 - 8.d0*fzd*zeta3 - 12.d0*fz*zeta2 )/fzdd

        rsdn2    = rsdn*rsdn
        rsdn3    = rsdn2*rsdn
        zetada2  = zetada**2
        zetada3  = zetada2*zetada
        zetadb2  = zetadb**2
        zetadb3  = zetadb2*zetadb

        euda     = eurs*rsdn + euzt*zetada
        eudb     = eurs*rsdn + euzt*zetadb
        eudaa    = eudrr*rsdn2 + 2.d0*eudrz*rsdn*zetada + eudzz*zetada2 + eurs*rsdnn + euzt*zetadaa
        eudab    = eudrr*rsdn2 + eudrz*rsdn*( zetada + zetadb ) + eudzz*zetada*zetadb &
             + eurs*rsdnn + euzt*zetadab
        eudbb    = eudrr*rsdn2 + 2.d0*eudrz*rsdn*zetadb + eudzz*zetadb2 + eurs*rsdnn + euzt*zetadbb
        eudaaa   = eudrrr*rsdn3 &
             + 3.d0*eudrrz*rsdn2*zetada &
             + 3.d0*eudrzz*rsdn*zetada2 + eudzzz*zetada3 &
             + 3.d0*eudrr*rsdn*rsdnn + 3.d0*eudzz*zetada*zetadaa &
             + 3.d0*eudrz*( rsdnn*zetada + rsdn*zetadaa ) &
             + eurs*rsdnnn + euzt*zetadaaa
        eudaab   = eudrrr*rsdn3 &
             + eudrrz*rsdn2*( zetadb + 2.d0*zetada ) &
             + eudrzz*rsdn*zetada*( zetada + 2.d0*zetadb ) &
             + eudzzz*zetada2*zetadb &
             + 3.d0*eudrr*rsdn*rsdnn &
             + eudrz*( rsdnn*( 2.d0*zetada + zetadb ) + rsdn*( 2.d0*zetadab + zetadaa ) ) &
             + eudzz*( 2.d0*zetada*zetadab + zetadaa*zetadb ) &
             + eurs*rsdnnn + euzt*zetadaab
        eudabb   = eudrrr*rsdn3 &
             + eudrrz*rsdn2*( zetada + 2.d0*zetadb ) &
             + eudrzz*rsdn*zetadb*( zetadb + 2.d0*zetada ) &
             + eudzzz*zetada*zetadb2 &
             + 3.d0*eudrr*rsdn*rsdnn &
             + eudrz*( rsdnn*( zetada + 2.d0*zetadb ) + rsdn*( 2.d0*zetadab + zetadbb ) ) &
             + eudzz*( zetada*zetadbb + 2.d0*zetadb*zetadab ) &
             + eurs*rsdnnn + euzt*zetadabb
        eudbbb   = eudrrr*rsdn3 &
             + 3.d0*eudrrz*rsdn2*zetadb &
             + 3.d0*eudrzz*rsdn*zetadb2 + eudzzz*zetadb3 &
             + 3.d0*eudrr*rsdn*rsdnn + 3.d0*eudzz*zetadb*zetadbb &
             + 3.d0*eudrz*( rsdnn*zetadb + rsdn*zetadbb ) &
             + eurs*rsdnnn + euzt*zetadbbb

        !print '(i2,12d19.10)', i,chgrhr_l(i,1:2),eu,euda,eudb,eudaa,eudbb,eudab,eudaaa,eudaab,eudabb,eudbbb
        !print *,eu,     e1+(e2-e1)* fz * zeta4 - e3  * fz * onzt4 / fzdd
        !print *,eurs,   e1dr  + ( e1pdr  -  e1dr  ) * fz * zeta4 - e1qdr  * fz * onzt4 / fzdd
        !print *,euzt,   (e2-e1)*( fzd*zeta4 + 4.d0*fz*zeta3 )-e3*( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd 

        duda     = gzt*zetada
        dudb     = gzt*zetadb 
        dudaa    = gzz*zetada2     + gzt*zetadaa
        dudab    = gzz*zetada*zetadb + gzt*zetadab
        dudbb    = gzz*zetadb2     + gzt*zetadbb
        dudaaa   = gzzz*zetada3    + 3.d0*gzz*zetada*zetadaa   + gzt*zetadaaa
        dudaab   = gzzz*zetada2*zetadb + gzz*( zetadaa*zetadb + 2.d0*zetada*zetadab ) + gzt*zetadaab
        dudabb   = gzzz*zetada*zetadb2 + gzz*( zetada*zetadbb + 2.d0*zetadab*zetadb ) + gzt*zetadabb
        dudbbb   = gzzz*zetadb3    + 3.d0*gzz*zetadb*zetadbb   + gzt*zetadbbb 

        !print '(i2,12d19.10)', i,chgrhr_l(i,1:2),g,duda,dudb,dudaa,dudbb,dudab,dudaaa,dudaab,dudabb,dudbbb


        dtdn  = -sixth7*t/d
        dtdnn = -sixth13*dtdn/d
        dtdnnn= -sixth19*dtdnn/d
        !print '(i2,5e19.12)',i,d,t,dtdn,dtdnn,dtdnnn
        if(dabs(grad_trho(i)) > 1.0d-9) then
           dtdg  = t/dd
           dtdng = -sixth7*dtdg/d
           dtdnng= -sixth13*dtdng/d
           !               dtdnug = -dtdn*invg/dd
        else
           dtdg  = 0.d0
           dtdng = 0.d0
           dtdnng= 0.d0
           dtdnug= 0.d0 
        endif
        dtdgg  = 0.d0
        dtdngg = 0.d0
        dtdggg = 0.d0
        invg=1.d0/g
        dtdu   = -t*invg
        dtduu  = -2.d0*dtdu*invg 
        dtduuu = -3.d0*dtduu*invg
        dtdnu  = -dtdn*invg
        dtdnnu = -dtdnn*invg
        dtdnuu = -sixth7*dtduu/d
        dtdug  = -dtdg*invg
        dtduug = 2.d0*dtdg*invg**2
        dtdugg = 0.d0
        if(dabs(grad_trho(i)) > 1.0d-9) then
           dtdnug = dtdnu/dd
        else
           dtdnug = 0.d0
        end if

        dtda   = dtdn + dtdu*duda
        dtdb   = dtdn + dtdu*dudb
        dtdaa  = dtdnn + 2.d0*dtdnu*duda + dtduu*duda**2 + dtdu*dudaa
        dtdbb  = dtdnn + 2.d0*dtdnu*dudb + dtduu*dudb**2 + dtdu*dudbb
        dtdab  = dtdnn + dtdnu*( duda + dudb ) + dtduu*duda*dudb + dtdu*dudab
        dtdag  = dtdng + dtdug*duda
        dtdbg  = dtdng + dtdug*dudb
        dtdaaa = dtdnnn+ 3.d0*dtdnnu*duda + 3.d0*dtdnuu*duda**2 + dtduuu*duda**3 &
             + 3.d0*dtduu*duda*dudaa + 3.d0*dtdnu*dudaa + dtdu*dudaaa
        dtdbbb = dtdnnn+ 3.d0*dtdnnu*dudb + 3.d0*dtdnuu*dudb**2 + dtduuu*dudb**3 &
             + 3.d0*dtduu*dudb*dudbb + 3.d0*dtdnu*dudbb + dtdu*dudbbb
        dtdaab = dtdnnn+ dtdnnu*( 2.d0*duda + dudb ) + dtdnuu*duda*( duda + 2.d0*dudb ) + dtduuu*duda**2*dudb &
             + dtduu*( dudaa*dudb + 2.d0*duda*dudab ) + dtdnu*( dudaa + 2.d0*dudab ) + dtdu*dudaab
        dtdabb = dtdnnn+ dtdnnu*( duda + 2.d0*dudb ) + dtdnuu*dudb*( dudb + 2.d0*duda ) + dtduuu*dudb**2*duda &
             + dtduu*( dudbb*duda + 2.d0*dudb*dudab ) + dtdnu*( dudbb + 2.d0*dudab ) + dtdu*dudabb 
        dtdaag = dtdnng+ 2.d0*dtdnug*duda + dtduug*duda**2 + dtdug*dudaa
        dtdbbg = dtdnng+ 2.d0*dtdnug*dudb + dtduug*dudb**2 + dtdug*dudbb
        dtdabg = dtdnng+ dtdnug*( duda + dudb ) + dtduug*duda*dudb + dtdug*dudab 

        !print '(i2,20d19.10)', i,t,dtda,dtdb,dtdg,dtdaa,dtdbb,dtdgg,dtdab,dtdag,dtdbg,dtdaaa,dtdbbb,dtdggg,dtdaab,dtdabb,dtdaag,dtdagg,dtdbbg,dtdbgg,dtdabg     
        !print '(i2,12d19.10)', i,t,dtdag,dtdbg,dtdaag,dtdbbg,dtdabg     

        invg3  = 1.d0  /g3
        dAde   = b * ( b + delt ) / bet * invg3
        dAdee  = 2.d0 * dAde * ( b + 0.5d0 * delt ) / bet *invg3
        dAdeee = 2.d0 * ( dAde**2 + ( b + 0.5d0 * delt )*dAdee ) /bet * invg3
        m3u    = -3.d0*invg
        dAdu     = m3u*eu*dAde
        dAdeu    = m3u*( dAde + eu*dAdee )
        dAduu    = m3u*eu*invg*( -dAde + g*dAdeu )
        dAdeeu   = m3u*( 2.d0*dAdee + eu*dAdeee )
        dAdeuu   = m3u*invg*( -dAde + g*dAdeu - eu*dAdee + g*eu*dAdeeu )
        dAduuu   = m3u*eu*invg**2*( 2.d0*dAde -2.d0*g*dAdeu + g**2*dAdeuu )

        dAda     = dAde*euda + dAdu*duda
        dAdb     = dAde*eudb + dAdu*dudb
        dAdaa    = dAdee*euda**2 + 2.d0*dAdeu*euda*duda + dAduu*duda**2 + dAde*eudaa + dAdu*dudaa 
        dAdab    = dAdee*euda*eudb + dAdeu*( euda*dudb + eudb*duda ) + dAduu*duda*dudb + dAde*eudab + dAdu*dudab 
        dAdbb    = dAdee*eudb**2 + 2.d0*dAdeu*eudb*dudb + dAduu*dudb**2 + dAde*eudbb + dAdu*dudbb
        dAdaaa   = dAdeee*euda**3 + 3.d0*dAdeeu*euda**2*duda + 3.d0*dAdeuu*euda*duda**2 + dAduuu*duda**3 &
             + 3.d0*dAdee*euda*eudaa + 3.d0*dAduu*duda*dudaa &
             + 3.d0*dAdeu*( euda*dudaa + eudaa*duda ) + dAde*eudaaa + dAdu*dudaaa
        dAdaab   = dAdeee*euda**2*eudb + dAduuu*duda**2*dudb &
             + dAdeeu*euda*( euda*dudb + 2.d0*eudb*duda ) &
             + dAdeuu*( 2.d0*euda*dudb + eudb*duda )*duda &
             + dAdee*( 2.d0*euda*eudab + eudb*eudaa ) + dAduu*( 2.d0*duda*dudab + dudb*dudaa ) &
             + dAdeu*( 2.d0*eudab*duda + 2.d0*euda*dudab + dudb*eudaa + eudb*dudaa ) &
             + dAde*eudaab + dAdu*dudaab
        dAdabb   = dAdeee*euda*eudb**2 + dAduuu*dudb**2*duda &
             + dAdeeu*eudb*( eudb*duda + 2.d0*euda*dudb ) &
             + dAdeuu*( 2.d0*eudb*duda + euda*dudb )*dudb &
             + dAdee*( eudbb*euda + 2.d0*eudb*eudab ) + dAduu*( dudbb*duda + 2.d0*dudb*dudab ) &
             + dAdeu*( eudbb*duda + 2.d0*eudb*dudab + 2.d0*dudb*eudab + dudbb*euda ) &
             + dAde*eudabb + dAdu*dudabb
        dAdbbb   = dAdeee*eudb**3 + 3.d0*dAdeeu*eudb**2*dudb + 3.d0*dAdeuu*eudb*dudb**2 + dAduuu*dudb**3 &
             + 3.d0*dAdee*eudb*eudbb + 3.d0*dAduu*dudb*dudbb &
             + 3.d0*dAdeu*( eudb*dudbb + eudbb*dudb ) + dAde*eudbbb + dAdu*dudbbb 

        !print '(i2,12d19.10)', i,b,dAda,dAdb,dAdaa,dAdab,dAdbb,dAdaaa,dAdaab,dAdabb,dAdbbb

        at2  = b * t2
        a2t4 = at2**2
        a3t6 = at2*a2t4
        t3   = t2*t
        t5   = t4*t
        t7   = t6*t
        t8   = t7*t
        fE    = q5**2 + delt*t2*q4*q5
        fEdt  = 4.d0*b*t*q5*( 1.d0 + 2.d0*at2 ) + &
             delt *2.d0*t*(1.d0 + 4.d0*at2 + 6.d0*a2t4 + 4.d0*a3t6 )
        fEdtt = 4.d0*b*( 1.d0 + 9.d0*at2 + 15.d0*a2t4 + 14.d0*a3t6 ) + &
             delt*( 2.d0 + 24.d0*at2 + 60.d0*a2t4 + 56.d0*a3t6 )
        fEda  = 2.d0*t2*q5*( 1.d0 + 2.d0*at2 ) + &
             delt*t4*( 2.d0 + 4.d0*at2 + 3.d0*a2t4 )
        fEdaa = 6.d0*t4*( 2.d0*q5 - 1.d0) + &
             delt*2.d0*t6*( 2.d0 + 3.d0*at2 )
        fF    = t*( 1.d0 + 2.d0*at2 )         
        fFdt  = 1.d0 + 6.d0*at2
        fFdtt = 12.d0*b*t
        fFda  = 2.d0*t3
        fFdaa = 0.d0
        fG    = 2.d0*b*t6 + b2*t8
        fGdt  = 12.d0*b*t5 + 8.d0*b2*t7
        fGdtt = 60.d0*b*t4 + 56.d0*b2*t6
        fGda  = 2.d0*t6 + 2.d0*b*t8
        fGdaa = 2.d0*t8

        betivE=bet/fE*g3
        betivE2=betivE/fE
        dHdf   =     -fG*betivE
        dHdt   = 2.d0*fF*betivE
        dHdff  =     -( fGda*fE-fG*fEda )*betivE2
        dHdft  =     -( fGdt*fE-fG*fEdt )*betivE2
        dHdtt  = 2.d0*( fFdt*fE-fF*fEdt )*betivE2
        dHdfff =     -( fGdaa*fE-fG*fEdaa-2.d0*( fGda*fE-fG*fEda )*fEda/fE )*betivE2
        dHdfft = 2.d0*( fFdaa*fE-fF*fEdaa-2.d0*( fFda*fE-fF*fEda )*fEda/fE )*betivE2
        dHdftt =     -( fGdtt*fE-fG*fEdtt-2.d0*( fGdt*fE-fG*fEdt )*fEdt/fE )*betivE2
        dHdttt = 2.d0*( fFdtt*fE-fF*fEdtt-2.d0*( fFdt*fE-fF*fEdt )*fEdt/fE )*betivE2

        tiu    = 3.d0*invg
        siu2   = tiu*2.d0*invg
        dHdu   = tiu*h
        dHduu  = siu2*h
        dHduuu = siu2*invg*h
        dHdfu  = tiu*dHdf
        dHdtu  = tiu*dHdt
        dHdffu = tiu*dHdff
        dHdfuu = siu2*dHdf
        dHdttu = tiu*dHdtt
        dHdtuu = siu2*dHdt
        dHdftu = tiu*dHdft

        dHda   = dHdf*dAda + dHdt*dtda + dHdu*duda
        dHdb   = dHdf*dAdb + dHdt*dtdb + dHdu*dudb
        dHdg   =             dHdt*dtdg
        dHdaa  = dHdff*dAda**2 + dHdtt*dtda**2 + dHduu*duda**2 &
             + 2.d0*dHdft*dAda*dtda + 2.d0*dHdfu*dAda*duda + 2.d0*dHdtu*dtda*duda &
             + dHdf*dAdaa + dHdt*dtdaa + dHdu*dudaa
        dHdab  = dHdff*dAda*dAdb + dHdtt*dtda*dtdb + dHduu*duda*dudb &
             + dHdft*( dAda*dtdb + dAdb*dtda ) + dHdfu*( dAda*dudb + dAdb*duda ) + dHdtu*( dtda*dudb + dtdb*duda ) &
             + dHdf*dAdab + dHdt*dtdab + dHdu*dudab
        dHdbb  = dHdff*dAdb**2 + dHdtt*dtdb**2 + dHduu*dudb**2 &
             + 2.d0*dHdft*dAdb*dtdb + 2.d0*dHdfu*dAdb*dudb + 2.d0*dHdtu*dtdb*dudb &
             + dHdf*dAdbb + dHdt*dtdbb + dHdu*dudbb
        dHdgg  = dHdtt*dtdg**2
        dHdag  = dHdft*dAda*dtdg + dHdtt*dtda*dtdg + dHdtu*duda*dtdg + dHdt*dtdag
        dHdbg  = dHdft*dAdb*dtdg + dHdtt*dtdb*dtdg + dHdtu*dudb*dtdg + dHdt*dtdbg
        dHdaaa = dHdfff*dAda**3 + 3.d0*dHdfft*dAda**2*dtda + 3.d0*dHdffu*dAda**2*duda &
             + 3.d0*dHdftt*dAda*dtda**2 + dHdttt*dtda**3 + 3.d0*dHdttu*dtda**2*duda &
             + 3.d0*dHdfuu*dAda*duda**2 + 3.d0*dHdtuu*dtda*duda**2 + dHduuu*duda**3 &
             + 6.d0*dHdftu*dAda*dtda*duda &
             + 3.d0*dHdff*dAda*dAdaa + 3.d0*dHdtt*dtda*dtdaa + 3.d0*dHduu*duda*dudaa &
             + 3.d0*dHdft*( dAdaa*dtda + dAda*dtdaa ) &
             + 3.d0*dHdtu*( dtdaa*duda + dtda*dudaa ) &
             + 3.d0*dHdfu*( dAdaa*duda + dAda*dudaa ) &
             + dHdf*dAdaaa + dHdt*dtdaaa + dHdu*dudaaa
        dHdaab = dHdfff*dAda**2*dAdb + dHdfft*dAda*( dAda*dtdb + 2.d0*dAdb*dtda ) &
             + dHdffu*dAda*( dAda*dudb + 2.d0*dAdb*duda ) &
             + dHdftt*dtda*( dAdb*dtda + 2.d0*dAda*dtdb ) &
             + dHdttt*dtdb*dtda**2 + dHdttu*dtda*( dtda*dudb + 2.d0*dtdb*duda ) &
             + dHdfuu*duda*( dAdb*duda + 2.d0*dAda*dudb ) &
             + dHdtuu*duda*( dtdb*duda + 2.d0*dtda*dudb ) + dHduuu*duda**2*dudb &
             + 2.d0*dHdftu*( dAdb*duda*dtda + dAda*dudb*dtda + dAda*duda*dtdb ) &
             + dHdff*( 2.d0*dAda*dAdab + dAdb*dAdaa ) &
             + dHdtt*( 2.d0*dtda*dtdab + dtdb*dtdaa ) &
             + dHduu*( 2.d0*duda*dudab + dudb*dudaa ) &
             + dHdft*( 2.d0*dAdab*dtda + 2.d0*dAda*dtdab &
             + dAdaa*dtdb + dAdb*dtdaa ) &
             + dHdtu*( 2.d0*dtdab*duda + 2.d0*dtda*dudab &
             + dtdaa*dudb + dtdb*dudaa ) &
             + dHdfu*( 2.d0*dAdab*duda + 2.d0*dAda*dudab &
             + dAdaa*dudb + dAdb*dudaa ) &
             + dHdf*dAdaab + dHdt*dtdaab + dHdu*dudaab
        dHdabb = dHdfff*dAda*dAdb**2 + dHdfft*dAdb*( dAdb*dtda + 2.d0*dAda*dtdb ) &
             + dHdffu*dAdb*( dAdb*duda + 2.d0*dAda*dudb ) &
             + dHdftt*dtdb*( dAda*dtdb + 2.d0*dAdb*dtda ) + dHdttt*dtda*dtdb**2 &
             + dHdttu*dtdb*( dtdb*duda + 2.d0*dtda*dudb ) &
             + dHdfuu*dudb*( dAda*dudb + 2.d0*dAdb*duda ) &
             + dHdtuu*dudb*( dtda*dudb + 2.d0*dtdb*duda ) + dHduuu*duda*dudb**2 &
             + 2.d0*dHdftu*( dAda*dtdb*dudb + dAdb*dtda*dudb + dAdb*dtdb*duda ) &
             + dHdff*( 2.d0*dAdb*dAdab + dAda*dAdbb ) &
             + dHdtt*( 2.d0*dtdb*dtdab + dtda*dtdbb ) &
             + dHduu*( 2.d0*dudb*dudab + duda*dudbb ) &
             + dHdft*( 2.d0*dAdab*dtdb + 2.d0*dAdb*dtdab + dAdbb*dtda + dAda*dtdbb ) &
             + dHdtu*( 2.d0*dtdab*dudb + 2.d0*dtdb*dudab + dtdbb*duda + dtda*dudbb ) &
             + dHdfu*( 2.d0*dAdab*dudb + 2.d0*dAdb*dudab + dAdbb*duda + dAda*dudbb ) &
             + dHdf*dAdabb + dHdt*dtdabb + dHdu*dudabb
        dHdbbb = dHdfff*dAdb**3 + 3.d0*dHdfft*dAdb**2*dtdb + 3.d0*dHdffu*dAdb**2*dudb &
             + 3.d0*dHdftt*dAdb*dtdb**2 + dHdttt*dtdb**3 + 3.d0*dHdttu*dtdb**2*dudb &
             + 3.d0*dHdfuu*dAdb*dudb**2 + 3.d0*dHdtuu*dtdb*dudb**2 + dHduuu*dudb**3 &
             + 6.d0*dHdftu*dAdb*dtdb*dudb &
             + 3.d0*dHdff*dAdb*dAdbb + 3.d0*dHdtt*dtdb*dtdbb + 3.d0*dHduu*dudb*dudbb &
             + 3.d0*dHdft*( dAdbb*dtdb + dAdb*dtdbb ) &
             + 3.d0*dHdtu*( dtdbb*dudb + dtdb*dudbb ) &
             + 3.d0*dHdfu*( dAdbb*dudb + dAdb*dudbb ) &
             + dHdf*dAdbbb + dHdt*dtdbbb + dHdu*dudbbb
        dHdabg = dHdfft*dAda*dAdb*dtdg   + dHdftt*( dAdb*dtda + dAda*dtdb )*dtdg &
             + dHdttt*dtda*dtdb*dtdg + dHdttu*( dtdb*duda + dtda*dudb )*dtdg &
             + dHdtuu*duda*dudb*dtdg + dHdftu*( dAdb*duda + dAda*dudb )*dtdg &
             + dHdft*( dAdab*dtdg + dAdb*dtdag + dAda*dtdbg ) &
             + dHdtt*( dtdab*dtdg + dtdb*dtdag + dtda*dtdbg ) &
             + dHdtu*( dudab*dtdg + dudb*dtdag + duda*dtdbg ) &
             + dHdt*dtdabg
        dHdaag = dHdfft*dAda**2*dtdg + dHdttt*dtda**2*dtdg + dHdtuu*duda**2*dtdg &
             + 2.d0*dHdttu*dtda*duda*dtdg + 2.d0*dHdftt*dAda*dtda*dtdg + 2.d0*dHdftu*dAda*duda*dtdg &
             + dHdtt*( dtdaa*dtdg + 2.d0*dtda*dtdag ) + dHdft*( dAdaa*dtdg + 2.d0*dAda*dtdag ) &
             + dHdtu*( dudaa*dtdg + 2.d0*duda*dtdag ) + dHdt*dtdaag         
        dHdbbg = dHdfft*dAdb**2*dtdg + dHdttt*dtdb**2*dtdg + dHdtuu*dudb**2*dtdg &
             + 2.d0*dHdttu*dtdb*dudb*dtdg + 2.d0*dHdftt*dAdb*dtdb*dtdg + 2.d0*dHdftu*dAdb*dudb*dtdg &
             + dHdtt*( dtdbb*dtdg + 2.d0*dtdb*dtdbg ) + dHdft*( dAdbb*dtdg + 2.d0*dAdb*dtdbg ) &
             + dHdtu*( dudbb*dtdg + 2.d0*dudb*dtdbg ) + dHdt*dtdbbg       
        dHdagg = dHdftt*dAda*dtdg**2 + dHdttt*dtda*dtdg**2 + dHdttu*duda*dtdg**2 + 2.d0*dHdtt*dtdg*dtdag
        dHdbgg = dHdftt*dAdb*dtdg**2 + dHdttt*dtdb*dtdg**2 + dHdttu*dudb*dtdg**2 + 2.d0*dHdtt*dtdg*dtdbg 
        dHdggg = dHdttt*dtdg**3        

        !print '(i2,20d19.10)', i,h,dHda,dHdb,dHdg,dHdaa,dHdbb,dHdgg,dHdag,dHdbg,dHdab,dHdaaa,dHdbbb,dHdggg,dHdaab,dHdabb,dHdaag,dHdagg,dHdbbg,dHdbgg,dHdabg

        dFc_daa(i)   = 2.d0*euda + d*eudaa + 2.d0*dHda + d*dHdaa 
        dFc_dbb(i)   = 2.d0*eudb + d*eudbb + 2.d0*dHdb + d*dHdbb
        dFc_dgg(i)   = d*dHdgg
        dFc_dab(i)   = euda + eudb + d*eudab + dHda + dHdb + d*dHdab 
        dFc_dag(i)   = dHdg + d*dHdag
        dFc_dbg(i)   = dHdg + d*dHdbg
        dFc_daaa(i)  = 3.d0*eudaa + d*eudaaa + 3.d0*dHdaa + d*dHdaaa
        dFc_dbbb(i)  = 3.d0*eudbb + d*eudbbb + 3.d0*dHdbb + d*dHdbbb
        dFc_dggg(i)  = d*dHdggg
        dFc_daab(i)  = eudaa + 2.d0*eudab +d*eudaab + dHdaa + 2.d0*dHdab + d*dHdaab
        dFc_daag(i)  = 2.d0*dHdag + d*dHdaag
        dFc_dabb(i)  = eudbb + 2.d0*eudab + d*eudabb + dHdbb + 2.d0*dHdab + d*dHdabb
        dFc_dbbg(i)  = 2.d0*dHdbg + d*dHdbbg
        dFc_dagg(i)  = dHdgg + d*dHdagg
        dFc_dbgg(i)  = dHdgg + d*dHdbgg
        dFc_dabg(i)  = dHdag + dHdbg + d*dHdabg

        !print '(i2,21d19.10)', i,exc(i),d*(eu+h),dF_drho(i,1),dF_drho(i,2),dF_dgradrho(i)*grad_trho(i) &
        !                      ,dFc_daa(i),dFc_dbb(i),dFc_dgg(i),dFc_dab(i),dFc_dag(i),dFc_dbg(i),dFc_daaa(i),dFc_dbbb(i) &
        !                      ,dFc_dggg(i),dFc_daab(i),dFc_daag(i),dFc_dabb(i),dFc_dbbg(i),dFc_dagg(i),dFc_dbgg(i),dFc_dabg(i)
        !print *,dHdag,dHdaag,d,2.d0*dHdag + d*dHdaag
        !print *,2.d0*dHdbg + d*dHdbbg
        !print *,dHdag + dHdbg + d*dHdabg


     end do
  end if

end subroutine cr_gga_paw_library
! ===
! ==== other xc functional 
subroutine ex_gga_paw_library_3D( nrc, dnr, nspin, chgrhr_l, grad_rho, exc, &
     &                         dFx_drho, dFx_dgradrho, &
     &                         dFx_drr,  dFx_drg,  dFx_dgg, &
     &                         dFx_drrr, dFx_drrg, dFx_drgg, dFx_dggg, &
     &                         ista_nrc, iend_nrc, ist, ien, pot_type )
  use m_Const_Parameters,  only : DP,PAI
  implicit none

  integer,intent(in)        :: nrc, dnr, nspin, pot_type
  real(kind=DP),intent(in)  :: chgrhr_l(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(in)  :: grad_rho(ista_nrc:iend_nrc,nspin)

  real(kind=DP),intent(out) :: exc(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFx_drho(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dgradrho(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drr(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drg(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dgg(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drrr(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drrg(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_drgg(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dFx_dggg(ista_nrc:iend_nrc,nspin)
  integer :: ista_nrc, iend_nrc, ist, ien

  real(kind=DP), parameter ::  ax = -0.7385587663820224d0
  real(kind=DP), parameter :: thrd = 0.33333333333333333333333d0, &
       &                      thrd4 = 1.3333333333333333333333d0
  real(kind=DP), parameter :: thpith = 3.0936677262801d0      ! ???

  real(kind=DP), parameter :: ninth2=0.2222222222222222222d0
  real(kind=DP), parameter :: ninth4=0.4444444444444444444d0
  real(kind=DP), parameter :: ninth8=0.8888888888888888888d0
  real(kind=DP), parameter :: thrd2 =0.6666666666666666666d0
  real(kind=DP), parameter :: eghth =0.125d0
  real(kind=DP), parameter :: twtysvnth8 = 0.29629629629629627985d0
  real(kind=DP), parameter :: twtysvnth32= 1.18518518518518511939d0

  ! ----- parameters of enhancement factor --
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

  real(kind=DP) :: facw, d, dd, fk, s, fac, s2, f, ex, exd, exdd, &
       &           exc0, excd, excdd, exc1, pot_add
  real(kind=DP) :: df_ds, d2f_ds2, d3f_ds3, ivd, ivdthrd2, s4, s6
  real(kind=DP) :: ctmp1, ctmp2, ctmp3, ctmp4, ctmp5, ctmp6, ctmp7, &
       &           x, dx_ds, d2x_ds2, d3x_ds3
  real(kind=DP) :: df_dx, d2f_dx2, d3f_dx3
  real(kind=DP) :: dmu_dx, d2mu_dx2, d3mu_dx3
  real(kind=DP) :: alpha_pbeint, mu_this, dmu_ds, d2mu_ds2, d3mu_ds3
  real(kind=DP) :: term(4), dterm1(4), dterm2(4), dterm3(4)
  real(kind=DP) :: y, y2, y3, y4, dy_ds, d2y_ds2, d3y_ds3

  logical :: First = .true.
  real(kind=DP), save :: coeff_spline_htbs(0:5)

  if ( pot_type == 5 ) then
     if ( First ) then
        call init_spline_htbs( mu_wc, kappa_rpbe, kappa_wc, coeff_wc, &
             &                 s_low_htbs, s_high_htbs, coeff_spline_htbs )
        First = .false.
     endif
  endif

  !---- Spin dependency

  facw = nspin

  exc  = 0.d0
  dFx_drho = 0.0d0;    dFx_dgradrho = 0.0d0;    dFx_drr  = 0.0d0
  dFx_drg  = 0.0d0;    dFx_dgg      = 0.0d0;    dFx_drrr = 0.0d0
  dFx_drrg = 0.0d0;    dFx_drgg     = 0.0d0;    dFx_dggg = 0.0d0

  do is = 1, nspin
     exc1 = 0.d0

     do i = ist,ien,dnr
        d  = facw * chgrhr_l(i, is)
        if (d.le.0.0d0) cycle

        dd = facw * grad_rho(i, is)

        if (d > 1.d-05) then
           fk = (3.d0*PAI*PAI*d)**thrd
!!$        sk = dsqrt(4*fk/PAI)
           s = dd/(d*fk*2.d0)
        else
           fk = 0.d0
           s  = 0.d0
        endif
        !-------------------------------------
        fac = ax*d**thrd
        s2  = s*s

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
           x = mu_pbe *s2;   dx_ds = 2.0d0 *mu_pbe *s;    d2x_ds2 = 2.0d0 *mu_pbe
           d3x_ds3 = 0.0d0

           ctmp1 = kappa_pbe + x

           f = 1.0d0 +kappa_pbe *x /ctmp1
           df_dx = kappa_pbe**2 /ctmp1**2
           d2f_dx2 = -2.0d0 *df_dx /ctmp1
           d3f_dx3 = -3.0d0 *d2f_dx2 /ctmp1

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (2)                        ! revpbe
           x = mu_revpbe *s2;   dx_ds = 2.0d0 *mu_revpbe *s;  d2x_ds2 = 2.0d0 *mu_revpbe
           d3x_ds3 = 0.0d0

           ctmp1 = kappa_revpbe + x

           f = 1.0d0 +kappa_revpbe *x /ctmp1
           df_dx = kappa_revpbe**2 /ctmp1**2
           d2f_dx2 = -2.0d0 *df_dx /ctmp1
           d3f_dx3 = -3.0d0 *d2f_dx2 /ctmp1

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (3)                        ! rpbe
           x = mu_rpbe *s2;   dx_ds = 2.0d0 *mu_rpbe *s;    d2x_ds2 = 2.0d0 *mu_rpbe
           d3x_ds3 = 0.0d0

           ctmp1 = exp( -x/ kappa_rpbe )

           f = 1.0d0 +kappa_rpbe -kappa_rpbe *ctmp1
           df_dx = ctmp1
           d2f_dx2 = -df_dx /kappa_rpbe
           d3f_dx3 = -d2f_dx2 /kappa_rpbe

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (4)                        ! W. Cohen
           s4 = s2 *s2
           ctmp1 = exp( -s2 )
           ctmp2 = coeff_WC *s4

           x = mu_gel *s2 + ( mu_wc -mu_gel )*s2 *ctmp1 + log( 1.0d0 +ctmp2 )
           dx_ds = 2.0d0 *mu_gel *s + ( mu_wc -mu_gel ) *2.0d0 *s *ctmp1 *( 1.0d0 -s2 ) &
                &  +4.0d0 *coeff_WC *s2 *s /( 1.0d0 +ctmp2 )
           d2x_ds2 = 2.0d0 *mu_gel &
                &   + 2.0d0 *( mu_wc -mu_gel ) *ctmp1 *( 1.0d0 -5.0d0 *s2 +2.0d0 *s4 ) &
                &   + 4.0d0 *coeff_WC *s2 *( 3.0d0 -ctmp2 ) /( 1.0d0 +ctmp2 )**2
           d3x_ds3 = -4.0d0 *s *( mu_wc -mu_gel ) *ctmp1 &
                &                   *( 6.0d0 -9.0d0 *s2 +2.0d0 *s4 ) &
                &    -8.0d0 *coeff_WC *s *( 12.0d0 *ctmp2 -ctmp2**2 -3.0d0 ) &
                &            / ( 1.0d0 +ctmp2 )**3
           
           ctmp3 = kappa_wc +x

           f = 1.0d0 +kappa_wc *x /ctmp3
           df_dx = kappa_wc**2 /ctmp3**2
           d2f_dx2 = -2.0d0 *df_dx /ctmp3
           d3f_dx3 = -3.0d0 *d2f_dx2 /ctmp3

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (5)                       ! HTBS
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
           df_ds = coeff_spline_htbs(1) + 2.0d0 *coeff_spline_htbs(2) *ctmp1 &
                &                       + 3.0d0 *coeff_spline_htbs(3) *ctmp2 &
                &                       + 4.0d0 *coeff_spline_htbs(4) *ctmp3 &
                &                       + 5.0d0 *coeff_spline_htbs(5) *ctmp4
           d2f_ds2 = 2.0d0 *coeff_spline_htbs(2) + 6.0d0 *coeff_spline_htbs(3) *ctmp1 &
                &                                +12.0d0 *coeff_spline_htbs(4) *ctmp2 &
                &                                +20.0d0 *coeff_spline_htbs(5) *ctmp3 
           d3f_ds3 = 6.0d0 *coeff_spline_htbs(3) +24.0d0 *coeff_spline_htbs(4) *ctmp1 &
                &                                +60.0d0 *coeff_spline_htbs(5) *ctmp2 

           df_ds   = df_ds   / (s_high_htbs -s_low_htbs )
           d2f_ds2 = d2f_ds2 / (s_high_htbs -s_low_htbs )
           d3f_ds3 = d3f_ds3 / (s_high_htbs -s_low_htbs )

        case (6)                       ! pbesol
           x = mu_pbesol *s2;   dx_ds = 2.0d0 *mu_pbesol *s;  d2x_ds2 = 2.0d0 *mu_pbesol
           d3x_ds3 = 0.0d0

           ctmp1 = kappa_pbesol + x

           f = 1.0d0 +kappa_pbesol *x /ctmp1
           df_dx = kappa_pbesol**2 /ctmp1**2
           d2f_dx2 = -2.0d0 *df_dx /ctmp1
           d3f_dx3 = -3.0d0 *d2f_dx2 /ctmp1

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (7)                       ! pbeint
           ctmp1 = mu_pbe - mu_gel
           
           alpha_pbeint = mu_gel**2 /kappa_pbeint /ctmp1
           ctmp2 = alpha_pbeint *s2
           ctmp3 = 1.0d0 +ctmp2
           
           mu_this = mu_gel +ctmp1 *ctmp2 /ctmp3
           dmu_ds =  2.0d0 *ctmp1 *alpha_pbeint *s  /ctmp3**2
           d2mu_ds2 = 2.0d0 *ctmp1 *alpha_pbeint *( 1.0d0 -3.0d0 *ctmp2 ) /ctmp3**3
           d3mu_ds3 = ctmp1 *alpha_pbeint**2 *( -24.0d0 * s )* ( 1.0d0 -ctmp2 ) &
                &           / ctmp3**4
           
           x = mu_this *s2;   dx_ds = dmu_ds *s2 +2.0d0 *mu_this *s
           d2x_ds2 = d2mu_ds2 *s2 +4.0d0 *dmu_ds *s + 2.0d0 *mu_this
           d3x_ds3 = d3mu_ds3 *s2 +6.0d0 *d2mu_ds2 *s + 6.0d0 *dmu_ds
           
           ctmp4= kappa_pbeint + x
           
           f = 1.0d0 +kappa_pbeint *x /ctmp4
           df_dx = kappa_pbeint**2 /ctmp1**2
           d2f_dx2 = -2.0d0 *df_dx /ctmp1
           d3f_dx3 = -3.0d0 *d2f_dx2 /ctmp1

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (11)                        ! B86r  I. Hamada
           x = mu_b86r *s2;   dx_ds = 2.0d0 *mu_b86r *s;    d2x_ds2 = 2.0d0 *mu_b86r
           d3x_ds3 = 0.0d0

           ctmp1 = 1.0d0 + x /kappa_b86r
           ctmp2 = ctmp1**(-4.0d0/5.0d0)

           ctmp3 = -0.8d0 /kappa_b86r *ctmp2 /ctmp1
           ctmp4 = -1.8d0 /kappa_b86r *ctmp3 /ctmp1
           ctmp5 = -2.8d0 /kappa_b86r *ctmp4 /ctmp1

           f = 1.0d0 + x *ctmp2

           df_dx = ctmp2 +ctmp3 *x
           d2f_dx2 = 2.0d0 *ctmp3 +ctmp4 *x
           d3f_dx3 = 3.0d0 *ctmp4 +ctmp5 *x

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (12)                      ! optpbe
           x = mu_optpbe *s2;   dx_ds = 2.0d0 *mu_optpbe *s;   d2x_ds2 = 2.0d0 *mu_optpbe
           d3x_ds3 = 0.0d0

           ctmp1 = kappa_optpbe + x

           f = 1.0d0 +kappa_optpbe *x /ctmp1
           df_dx = kappa_optpbe**2 /ctmp1**2
           d2f_dx2 = -2.0d0 *df_dx /ctmp1
           d3f_dx3 = -3.0d0 *d2f_dx2 /ctmp1

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (13)                         ! optb86b
           x = mu_optb86b *s2;   dx_ds = 2.0d0 *mu_optb86b *s;
           d2x_ds2 = 2.0d0 *mu_optb86b
           d3x_ds3 = 0.0d0

           ctmp1 = 1.0d0 + x /kappa_optb86b
           ctmp2 = ctmp1**(-4.0d0/5.0d0)

           ctmp3 = -0.8d0 /kappa_optb86b *ctmp2 /ctmp1
           ctmp4 = -1.8d0 /kappa_optb86b *ctmp3 /ctmp1
           ctmp5 = -2.8d0 /kappa_optb86b *ctmp4 /ctmp1

           f = 1.0d0 + x *ctmp2

           df_dx = ctmp2 +ctmp3 *x
           d2f_dx2 = 2.0d0 *ctmp3 +ctmp4 *x
           d3f_dx3 = 3.0d0 *ctmp4 +ctmp5 *x

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (14)                      ! pw86r
           s4 = s2 *s2;  s6 = s4 *s2

           x = 1.0d0 +15.d0 *a_pw86r *s2 +b_pw86r *s4 + c_pw86r *s6
           dx_ds = ( 30.0d0 *a_pw86r + 4.0d0 *b_pw86r *s2 +6.0d0 *c_pw86r *s4 ) *s
           d2x_ds2 = 30.0d0 *a_pw86r + 12.0d0 *b_pw86r *s2 +30.0d0 *c_pw86r *s4
           d3x_ds3 = ( 24.0d0 *b_pw86r +120.0d0 *c_pw86r *s2 ) *s

           f = x**(1.0d0/15.0d0)
           df_dx = f /15.0d0 /x
           d2f_dx2 = -14.0d0 /15.0d0 *df_dx /x
           d3f_dx3 = -29.0d0 /15.0d0 *d2f_dx2 /x

           df_ds = dx_ds *df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx + 3.0d0 *d2x_ds2 *dx_ds *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (15)                      ! c09x
           x = mu_c09 *s2;  dx_ds = 2.0d0 *mu_c09 *s
           d2x_ds2 = 2.0d0 *mu_c09
           d3x_ds3 = 0.0d0

           ctmp1 = alpha_c09 /mu_c09

           ctmp2 = exp( -ctmp1 *x )
           ctmp3 = exp( -ctmp1 *x /2.0d0 )

           f = 1.0d0 +x *ctmp2  + kappa_c09 *( 1.0d0 -ctmp3 )
           df_dx = ctmp2 *( 1.0d0 -ctmp1 *x ) +0.5d0 *kappa_c09 *ctmp1 *ctmp3
           d2f_dx2 = ctmp2 *ctmp1 *( ctmp1 *x -2.0d0 ) &
                &   -0.25d0 *kappa_c09 *ctmp1**2 *ctmp3
           d3f_dx3 = ctmp2 *ctmp1**2 *( -ctmp1 *x +3.0d0 ) &
                &   +0.125d0 *kappa_c09 *ctmp1**3 *ctmp3

           df_ds = dx_ds * df_dx
           d2f_ds2 = d2x_ds2 *df_dx +dx_ds**2 *d2f_dx2
           d3f_ds3 = d3x_ds3 *df_dx +3.0d0 *dx_ds *d2x_ds2 *d2f_dx2 +dx_ds**3 *d3f_dx3

        case (16)                      ! lv-pw86r
           s4 = s2 *s2;  s6 = s4 *s2

           ctmp1 = alpha_lvpw86r *s6
           ctmp2 = 1.0d0 /(1.0d0 +ctmp1)
           ctmp3 = 1.0d0 /(beta_lvpw86r +ctmp1)

           ctmp4 = 1.0d0 +15.d0 *a_pw86r *s2 +b_pw86r *s4 + c_pw86r *s6
           ctmp5 = ( 30.0d0 *a_pw86r + 4.0d0 *b_pw86r *s2 +6.0d0 *c_pw86r *s4 ) *s
           ctmp6 = 30.0d0 *a_pw86r + 12.0d0 *b_pw86r *s2 +30.0d0 *c_pw86r *s4
           ctmp7 = ( 24.0d0 *b_pw86r +120.0d0 *c_pw86r *s2 ) *s

           term(1) = ctmp2
           term(2) = 1.d0 +mu_lvpw86r *s2
           term(3) = ctmp1 *ctmp3
           term(4) = ctmp4 **(1.0d0/15.0d0)

           dterm1(1) = -6.0d0 *alpha_lvpw86r *ctmp2**2 *s4 *s
           dterm1(2) = 2.0d0 *mu_lvpw86r *s
           dterm1(3) =  6.0d0 *alpha_lvpw86r *beta_lvpw86r *ctmp3**2 *s4 *s
           dterm1(4) = term(4) /ctmp4 *ctmp5 /15.0d0

           dterm2(1) = 6.0d0 *alpha_lvpw86r *ctmp2**3 *s4 *( 7.0d0 *ctmp1 -5.0d0 )
           dterm2(2) = 2.0d0 *mu_lvpw86r
           dterm2(3) =-6.0d0 *alpha_lvpw86r *beta_lvpw86r *ctmp3**3 *s4 &
                &          *( 7.0d0 *ctmp1 -5.0d0 *beta_lvpw86r )
           dterm2(4) = term(4) *( -14.0d0/15.0d0 /ctmp4**2 *ctmp5**2 +ctmp6 /ctmp4 ) &
                &      / 15.d0

           dterm3(1) = 6.0d0 *alpha_lvpw86r *ctmp2**4 *s2 *s &
                &            *( -56.0d0 *ctmp1**2 +140.0d0 *ctmp1 -20.0d0 )
           dterm3(2) = 0.0d0
           dterm3(3) =-6.0d0 *alpha_lvpw86r *beta_lvpw86r *ctmp3**4 *s2 *s &
                &            *( -56.0d0 *ctmp1**2 +140.0d0 *ctmp1 *beta_lvpw86r &
                &               -20.0d0 *beta_lvpw86r**2 )
           dterm3(4) = term(4) *( 14.0d0 /15.0d0 *29.0d0 /15.0d0 /ctmp4**3 *ctmp5**3 &
                &                -3.0d0 *14.0d0 /15.0d0 /ctmp4**2 *ctmp5 *ctmp6 &
                &                +ctmp7 /ctmp4 ) /15.0d0

           f = term(1) *term(2) +term(3) *term(4)
           df_ds   = dterm1(1) *term(2) +term(1) *dterm1(2) &
                &   +dterm1(3) *term(4) +term(3) *dterm1(4)
           d2f_ds2 = dterm2(1) *term(2) +2.0d0 *dterm1(1) *dterm1(2) &
                &   + term(1)  *dterm2(2) &
                &   +dterm2(3) *term(4) +2.0d0 *dterm1(3) *dterm1(4) &
                &   + term(3) *dterm2(4)
           d3f_ds3 = dterm3(1) *term(2) +3.0d0 *dterm2(1) *dterm1(2) &
                &   +3.0d0 *dterm1(1) *dterm2(2) + term(1) *dterm3(2) &
                &   +dterm3(3) *term(4) +3.0d0 *dterm2(3) *dterm1(4) &
                &   +3.0d0 *dterm1(3) *dterm2(4) + term(3) *dterm3(4)

        case (20)            !  Engel Vosko 93
           s4 = s2 *s2;  s6 = s4 *s2

           x = 1.0D0 +a1_ev93 *s2 +a2_ev93 *s4 +a3_ev93 *s6
           y = 1.0D0 +b1_ev93 *s2 +b2_ev93 *s4 +b3_ev93 *s6

           dx_ds = s *( 2.0D0 *a1_ev93 +4.0D0 *a2_ev93 *s2 +6.0D0 *a3_ev93 *s4 )
           dy_ds = s *( 2.0D0 *b1_ev93 +4.0D0 *b2_ev93 *s2 +6.0D0 *b3_ev93 *s4 )

           d2x_ds2 = 2.0D0 *a1_ev93 +12.0D0 *a2_ev93 *s2 +30.0D0 *a3_ev93 *s4
           d2y_ds2 = 2.0D0 *b1_ev93 +12.0D0 *b2_ev93 *s2 +30.0D0 *b3_ev93 *s4

           d3x_ds3 = s *( 24.0D0 *a2_ev93 +120.0D0 *a3_ev93 *s2 )
           d3y_ds3 = s *( 24.0D0 *b2_ev93 +120.0D0 *b3_ev93 *s2 )

           ctmp1 = y*dx_ds -x*dy_ds
           ctmp2 = y*d2x_ds2 -x*d2y_ds2
           ctmp3 = y*d3x_ds3 -x*d3y_ds3
           ctmp4 = dy_ds *d2x_ds2 -dx_ds *d2y_ds2

           y2 = y**2;    y3 = y2 *y;     y4 = y3 *y

           f = x /y
           df_ds = ctmp1 /y2
           d2f_ds2 = -2.0D0 *ctmp1 *dy_ds /y3 +ctmp2 /y2
           d3f_ds3 = ( 6.0D0 *dy_ds**2 /y4 -2.0D0 *d2y_ds2 /y3 ) *ctmp1 &
                &   -4.0D0 *dy_ds *ctmp2 /y3 &
                &   +( ctmp3 +ctmp4 ) /y2

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
           call phase_error_with_msg(6,"xc functional not supported (PAW)",__LINE__,__FILE__)
        end select

        ex  = fac *f *d
        exd = thrd4 *fac *( f - s *df_ds )
        exdd = ax *df_ds *0.5d0 /thpith
        !------------------------------------------     
        exc0 = ex / facw
        excd = exd
        dFx_drho(i, is) = excd 
        if ( pot_type == 25 ) dFx_drho(i,is) = dFx_drho(i,is) +pot_add

        excdd = exdd
        ! gradient of charge density 95/12/3 H.S.
        if(dabs(grad_rho(i, is)) > 1.d-9) then
           dFx_dgradrho(i, is) = excdd
           !           dFx_dgradrho(i, is) = excdd / grad_rho(i, is)           ! 2010 5
        else
           dFx_dgradrho(i, is) = 0.d0
        endif
        ! gradient of charge density 95/12/3 H.S.
        !        exc1 = exc1 + exc0*wos(i)
        exc(i)=exc(i)+exc0

        ! -----------------------------------
        ivd =1.d0 /d
        ivdthrd2 =ivd**thrd2

        if ( nspin == 1 ) then
           dFx_drr(i,is)   = ninth4 *ax *ivdthrd2*( f -s *df_ds +4.d0 *d2f_ds2 *s2  )
           dFx_drg(i,is)   =-thrd2 *ax /thpith *d2f_ds2 *s *ivd
           dFx_dgg(i,is)   = 0.25d0 *ax /thpith /thpith *d2f_ds2 *ivdthrd2**2
           dFx_drrr(i,is)  =-twtysvnth8 *ax *ivdthrd2 *ivd &
                &           *( f -df_ds *s +18.d0 *d2f_ds2 *s2 +8.d0 *d3f_ds3 *s2 *s )
           dFx_drrg(i,is)  = ninth2 *ax *ivd *ivd /thpith &
                &           *( 7.d0 *d2f_ds2 *s +4.d0 *d3f_ds3 *s2 )
           dFx_drgg(i,is)  =-thrd *ax /thpith /thpith *ivd *ivdthrd2**2 &
                &           *( d2f_ds2 +d3f_ds3 *s )
           dFx_dggg(i,is)  =eghth *ax /(thpith**3) *d3f_ds3 *ivd**2 *ivdthrd2
        else if ( nspin == 2 ) then
           dFx_drr(i,is)   = ninth8 *ax *ivdthrd2 *( f -s *df_ds +4.d0 *d2f_ds2 *s2 )
           dFx_drg(i,is)   =-thrd4 *ax /thpith *d2f_ds2 *s *ivd
           dFx_dgg(i,is)   = 0.5d0 *ax /thpith /thpith *d2f_ds2 *ivdthrd2**2
           dFx_drrr(i,is)  =-twtysvnth32 *ax *ivdthrd2 *ivd &
                &           *( f -df_ds*s +18.d0 *d2f_ds2 *s2 +8.d0 *d3f_ds3 *s2 *s )
           dFx_drrg(i,is)  = ninth8 *ax *ivd *ivd /thpith &
                &           *( 7.d0 *d2f_ds2 *s +4.d0 *d3f_ds3 *s2 )
           dFx_drgg(i,is)  =-thrd4 *ax/(thpith**2) *ivd *ivdthrd2**2 &
                &           *( d2f_ds2 +d3f_ds3 *s )
           dFx_dggg(i,is)  = 0.5d0 *ax/(thpith**3) *d3f_ds3 *ivd**2 *ivdthrd2
        end if

     end do
     !     exc(i) = exc(i) + exc1
  end do
end subroutine ex_gga_paw_library_3D

subroutine cr_gga_paw_library_3D( nrc, dnr, nspin, chgrhr_l, grad_trho, exc, &
     &                         dF_drho,  dF_dgradrho, dFc_daa,  dFc_dbb,  dFc_dgg ,&
     &                         dFc_dab,  dFc_dag,     dFc_dbg,  dFc_daaa, dFc_dbbb, &
     &                         dFc_dggg, dFc_daab,    dFc_daag, dFc_dabb, dFc_dbbg, &
     &                         dFc_dagg, dFc_dbgg, dFc_dabg, &
     &                         ista_nrc,iend_nrc,ist,ien,pot_type )
  use m_Const_Parameters,  only : PAI,DP
  implicit none

  integer,intent(in)        :: nrc, dnr, nspin, pot_type
  real(kind=DP),intent(in)  :: chgrhr_l(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(in) :: grad_trho(ista_nrc:iend_nrc)
  real(kind=DP),intent(inout) :: exc(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dF_drho(ista_nrc:iend_nrc,nspin)
  real(kind=DP),intent(out) :: dF_dgradrho(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_daa(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbb(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dgg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dab(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dag(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_daaa(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbbb(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dggg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_daab(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_daag(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dabb(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbbg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dagg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dbgg(ista_nrc:iend_nrc)
  real(kind=DP),intent(out) :: dFc_dabg(ista_nrc:iend_nrc)
  integer :: ista_nrc, iend_nrc, ist, ien

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

  ! vwn 95/12/2 Y.M
  real(kind=DP), parameter  :: aq  = 0.0168869d0,    a1q  = 0.11125d0

  real(kind=DP), parameter  :: b1q = 10.3570d0,     b2nq = 3.6231d0
  real(kind=DP), parameter  :: b3q = 0.88026d0,     b4q  = 0.49671d0

  ! vwn 95/12/2 Y.M
  real(kind=DP), parameter  :: gamma = 0.03109069086965489503494086371273d0

  real(kind=DP), parameter  :: thrd = 0.33333333333333330,&
       &                       sixth7 = 1.1666666666666666d0
  real(kind=DP), parameter  :: thrd7= 2.333333333333333333d0
  real(kind=DP), parameter  :: sixth13= 2.166666666666666666d0
  real(kind=DP), parameter  :: sixth19= 3.166666666666666666d0
  real(kind=DP), parameter  :: ninth= 0.11111111111111111111d0

  real(kind=DP), parameter  :: twntysvnth=  3.703703703703703d-2
  real(kind=DP), parameter  :: fftyfrth= 1.851851851851852d-2
  real(kind=DP), parameter  :: eityonth= 1.234567901234568d-2

! ----- parameters
  real(kind=DP), parameter :: beta_pbe   = 0.06672455060314922d0
  real(kind=DP), parameter :: beta_pbesol = 0.046d0
  real(kind=DP), parameter :: beta_pbeint = 0.052d0
!--------------------------

  integer       :: is,i
  real(kind=DP) :: facw, exc1,g,g3,facpon,d,dd,rs,fk,sk,t,s,q0 &
       &         , rs12,rs32,rsp,q1,q2,eu,q3,eurs,ecd,pon &
       &         , b,b2,t2,t4,t6,q4,q5,h0,h,q8,h0t,h0b,h0rs,ht,hrs &
       &         , ec1,ec1d,ec1dd,exc0,excd,excdd &
       &         , thrd2,thrd4,fzdd,zeta,onpzeta,onmzeta,fz,zetadxd &
       &         , q0p,q1p,q2p,q3p,q0q,q1q,q2q,q3q,fzd,onzt4,gzt,bzt,h0zt,euzt &
       &         , zetadxda,zetadxdb,ecda,ecdb,ec1da,ec1db,excda,excdb

  real(kind=DP) :: bet, delt

  real(kind=DP) :: rsdn,rsdnn,rsdnnn
  real(kind=DP) :: q1drr,q1drrr
  real(kind=DP) :: eudrr,eudrrr,invq12pq1,q12p1
  real(kind=DP) :: eudn,eudnn,eudnnn,ecdd,ecddd
  real(kind=DP) :: dtdn,dtdg,dtdnn,dtdgg,dtdng,dtdnnn,dtdnng,dtdngg,dtdggg
  real(kind=DP) :: dAde,dAdee,dAdeee
  real(kind=DP) :: dAdn,dAdnn,dAdnnn
  real(kind=DP) :: fE,fEdt,fEdtt,fEda,fEdaa
  real(kind=DP) :: fF,fFdt,fFdtt,fFda,fFdaa
  real(kind=DP) :: fG,fGdt,fGdtt,fGda,fGdaa
  real(kind=DP) :: at2,a2t4,a3t6,t3,t5,t7,t8
  real(kind=DP) :: dHda,dHdt,dHdaa,dHdat,dHdtt,dHdaaa,dHdaat,dHdatt,dHdttt,betivE,betivE2
  real(kind=DP) :: dHdn,dHdg,dHdnn,dHdng,dHdgg,dHdnnn,dHdnng,dHdngg,dHdggg
  real(kind=DP) :: zetada,zetadb,zetadaa,zetadab,zetadbb,zetadaaa,zetadaab,zetadabb,zetadbbb
  real(kind=DP) :: id2,id3,id4
  real(kind=DP) :: e1dr,e1drr,e1drrr
  real(kind=DP) :: e1pdr,e1pdrr,e1pdrrr
  real(kind=DP) :: e1qdr,e1qdrr,e1qdrrr
  real(kind=DP) :: q1pdrr,q1pdrrr,invq1p2pq1p,q1p2p1
  real(kind=DP) :: q1qdrr,q1qdrrr,invq1q2pq1q,q1q2p1
  real(kind=DP) :: fzdd0,fzddd,opzthrdm2,omzthrdm2
  real(kind=DP) :: zeta2,zeta3,zeta4,eudrz,eudzz,eudrrz,eudrzz,eudzzz
  real(kind=DP) :: e1,e2,e3
  real(kind=DP) :: euda,eudb,eudaa,eudbb,eudab,eudaaa,eudaab,eudabb,eudbbb
  real(kind=DP) :: rsdn2,rsdn3,zetada2,zetada3,zetadb2,zetadb3
  real(kind=DP) :: gzz,gzzz,duda,dudb,dudaa,dudab,dudbb,dudaaa,dudaab,dudabb,dudbbb
  real(kind=DP) :: dtdu,dtduu,dtduuu,dtdnu,dtdnnu,dtdnuu,dtdug,dtduug,dtdugg,dtdnug,invg
  real(kind=DP) :: dtda,dtdb,dtdaa,dtdab,dtdag,dtdbg,dtdbb
  real(kind=DP) :: dtdaaa,dtdaab,dtdabb,dtdbbb,dtdaag,dtdagg,dtdbbg,dtdbgg,dtdabg
  real(kind=DP) :: invg3,m3u
  real(kind=DP) :: dAdu,dAduu,dAduuu,dAdeu,dAdeeu,dAdeuu
  real(kind=DP) :: dAda,dAdb,dAdaa,dAdab,dAdbb,dAdaaa,dAdaab,dAdabb,dAdbbb
  real(kind=DP) :: dHdf,dHdff,dHdft,dHdfff,dHdfft,dHdftt
  real(kind=DP) :: dHdu,dHduu,dHduuu,dHdfu,dHdtu,dHdffu,dHdfuu,dHdttu,dHdtuu,dHdftu
  real(kind=DP) :: tiu,siu2
  real(kind=DP) :: dHdb,dHdbb,dHdag,dHdbg,dHdab,dHdbbb,dHdaab,dHdabb,dHdbbg,dHdbgg,dHdaag,dHdagg,dHdabg

  real(kind=DP), parameter  :: eta = 1.d-12
  real(kind=DP), parameter  :: sixthm = -0.1666666666666666666666d0

! -----
  select case ( pot_type )
  case (6)                           ! pbesol
     bet = beta_pbesol
  case (7)                           ! pbeint
     bet = beta_pbeint
  case default
     bet = beta_pbe
  end select

  delt = bet / gamma
! -----

  facw = nspin

  if ( nspin == 1 ) then
     g   = 1.d0
     g3  = g**3
!!$        g4  = g3*g
     facpon = -delt/(g3*bet)

     do i = ist,ien,dnr ! MPI
        d  = facw*chgrhr_l(i, 1)
        if(d < density_minimum) cycle
        ! gradient of charge density 95/12/3 H.S.
        dd = facw* grad_trho(i)
        rs = (0.75d0/(PAI*d))**thrd
        fk = (3.d0*PAI*PAI*d)**thrd
        sk = dsqrt(4.d0*fk/PAI)

        if(d > density_minimum2) then
           t = dd/(d*sk*2.d0)
           s = dd/(d*fk*2.d0)
        else
           t = 0.d0
           s = 0.d0
        endif

        q0 = -2.d0 *a *( 1.d0 +a1*rs )
        rs12 = dsqrt(rs)
        rs32 = rs12**3
        rsp = rs**p
        q1 = 2.d0 *a*( b1*rs12 +b2n*rs +b3*rs32 +b4*rs*rsp )
        q2 = log( 1.d0 +1.d0/q1 )
        eu = q0*q2
        q3 = a*( b1/rs12 +2.d0*b2n +3.d0 *b3*rs12 +2.d0*b4*p1*rsp )
        eurs = -2.d0*a*a1*q2 - q0*q3/(q1**2 + q1)
        ecd = eu-thrd*rs*eurs
        pon = facpon*eu
        b   = delt/(dexp(pon)-1.d0)
        b2  = b*b
        t2  = t*t
!!$           s2  = s*s
        t4  = t2*t2
        t6  = t4*t2
        q4  = 1.d0 + b*t2
        q5  = 1.d0 + b*t2 +b2*t4
        h0  = g3*(bet/delt)*dlog( 1.d0 +delt*q4*t2/q5 )
        h   = h0
        q8  = q5*q5 +delt*q4*q5*t2
        h0t = 2.d0*bet*t*( 1.d0 +2.d0*b*t2 )/q8 * g3
        h0b = -bet*t6*( 2.d0*b +b2*t2 )/q8 * g3
        h0rs = h0b*b*eurs*(b+delt)/ bet / g3
        ht   = h0t
        hrs  = h0rs
        ec1  = d*h/facw
        ec1d = h -thrd*rs*hrs -sixth7*t*ht
        ec1dd = 0.5d0*ht/sk

        exc0 = eu*d/facw + ec1
        excd = ecd + ec1d
        dF_drho(i, 1) = dF_drho(i, 1) + excd
        excdd = ec1dd
        ! gradient of charge density 95/12/3 H.S.
        if(dabs(grad_trho(i)) > 1.0d-9) then
           !              dF_dgradrho(i) = excdd / grad_trho(i)
           dF_dgradrho(i) = excdd
        else
           dF_dgradrho(i) = 0.d0
        endif
        ! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d/facw
        !           exc1 = exc1 + exc0*wos(i)
        exc(i)=exc(i)+exc0

        ! ***** Calc heigher derivatives. *****

        rsdn=-thrd*rs/d
        rsdnn=-4.d0*thrd*rsdn/d
        rsdnnn=-thrd7*rsdnn/d
        !print '(i2,5e19.10)',i,d,rs,rsdn,rsdnn,rsdnnn
        !q3 = a*(b1/rs12+2*b2n+3*b3*rs12+2*b4*p1*rsp)   
        q1drr =0.5d0*a*( -b1/rs32+3.d0*b3/rs12+8.d0*b4 )
        q1drrr=0.75d0*a*( b1/rs - b3 )/rs32
        !eurs = -2*a*a1*q2 - q0*q3/(q1**2 + q1)
        invq12pq1=1.d0/(q1**2 + q1)
        q12p1=2.d0*q1+1.d0
        eudrr=4.d0*a*a1*q3*invq12pq1-q0*( q1drr-q12p1*q3**2*invq12pq1 )*invq12pq1
        eudrrr=6.d0*a*a1*( q1drr - q12p1*q3**2*invq12pq1 ) * invq12pq1 &
             -q0*( q1drrr -(2.d0*q3**3+3.d0*q12p1*q3*q1drr -  &
             2.d0*q12p1**2*q3**3* invq12pq1)*invq12pq1 ) * invq12pq1
        eudn=eurs*rsdn
        eudnn=eudrr*rsdn**2+eurs*rsdnn
        eudnnn=eudrrr*rsdn**3+3.d0*eudrr*rsdn*rsdnn+eurs*rsdnnn
        !print '(i2,5e19.12)',i,d,eu,eudn,eudnn,eudnnn
        ecdd=2.d0*eudn+d*eudnn
        ecddd=3.d0*eudnn+d*eudnnn
        !print '(i2,5e19.12)',i,d,eu*d,ecdd,ecddd

        dtdn  = -sixth7*t/d
        dtdnn = -sixth13*dtdn/d
        dtdnnn= -sixth19*dtdnn/d
        !print '(i2,5e19.12)',i,d,t,dtdn,dtdnn,dtdnnn
        if(dabs(grad_trho(i)) > 1.0d-9) then
           dtdg  = t/dd
           dtdng = -sixth7*dtdg/d
           dtdnng= -sixth13*dtdng/d
        else
           dtdg  = 0.d0
           dtdng = 0.d0
           dtdnng= 0.d0
        endif
        dtdngg = 0.d0
        dtdggg = 0.d0

        dAde   = b * ( b + delt ) / bet
        dAdee  = 2.d0 * dAde * ( b + 0.5d0 * delt ) / bet
        dAdeee = 2.d0 * ( dAde**2 + ( b + 0.5d0 * delt )*dAdee ) /bet

        dAdn   = dAde  *eudn
        dAdnn  = dAdee *eudn**2  + dAde*eudnn
        dAdnnn = dAdeee*eudn**3  + 3.d0*dAdee*eudn*eudnn + dAde*eudnnn   
        !print '(i2,5e19.10)' ,i,d,b,dAdnn

        at2  = b * t2
        a2t4 = at2**2
        a3t6 = at2*a2t4
        t3   = t2*t
        t5   = t4*t
        t7   = t6*t
        t8   = t7*t
        fE    = q5**2 + delt*t2*q4*q5
        fEdt  = 4.d0*b*t*q5*( 1.d0 + 2.d0*at2 ) + &
             delt *2.d0*t*(1.d0 + 4.d0*at2 + 6.d0*a2t4 + 4.d0*a3t6 )
        fEdtt = 4.d0*b*( 1.d0 + 9.d0*at2 + 15.d0*a2t4 + 14.d0*a3t6 ) + &
             delt*( 2.d0 + 24.d0*at2 + 60.d0*a2t4 + 56.d0*a3t6 )
        fEda  = 2.d0*t2*q5*( 1.d0 + 2.d0*at2 ) + &
             delt*t4*( 2.d0 + 4.d0*at2 + 3.d0*a2t4 )
        fEdaa = 6.d0*t4*( 2.d0*q5 - 1.d0) + &
             delt*2.d0*t6*( 2.d0 + 3.d0*at2 )
        fF    = t*( 1.d0 + 2.d0*at2 )         
        fFdt  = 1.d0 + 6.d0*at2
        fFdtt = 12.d0*b*t
        fFda  = 2.d0*t3
        fFdaa = 0.d0
        fG    = 2.d0*b*t6 + b2*t8
        fGdt  = 12.d0*b*t5 + 8.d0*b2*t7
        fGdtt = 60.d0*b*t4 + 56.d0*b2*t6
        fGda  = 2.d0*t6 + 2.d0*b*t8
        fGdaa = 2.d0*t8

        betivE=bet/fE
        betivE2=betivE/fE
        dHda   =     -fG*betivE
        dHdt   = 2.d0*fF*betivE
        dHdaa  =     -( fGda*fE-fG*fEda )*betivE2
        dHdat  =     -( fGdt*fE-fG*fEdt )*betivE2
        dHdtt  = 2.d0*( fFdt*fE-fF*fEdt )*betivE2
        dHdaaa =     -( fGdaa*fE-fG*fEdaa-2.d0*( fGda*fE-fG*fEda )*fEda/fE )*betivE2
        dHdaat = 2.d0*( fFdaa*fE-fF*fEdaa-2.d0*( fFda*fE-fF*fEda )*fEda/fE )*betivE2
        dHdatt =     -( fGdtt*fE-fG*fEdtt-2.d0*( fGdt*fE-fG*fEdt )*fEdt/fE )*betivE2
        dHdttt = 2.d0*( fFdtt*fE-fF*fEdtt-2.d0*( fFdt*fE-fF*fEdt )*fEdt/fE )*betivE2

        dHdn     = dHda  *dAdn       + dHdt*dtdn
        dHdg     = dHdt  *dtdg
        dHdnn    = dHdaa *dAdn**2    + 2.d0*dHdat*dAdn*dtdn  + dHdtt*dtdn**2   + dHda*dAdnn  + dHdt*dtdnn
        dHdng    = dHdat *dAdn*dtdg  + dHdtt*dtdn*dtdg       + dHdt*dtdng
        dHdgg    = dHdtt*dtdg**2
        dHdnnn   = dHdaaa*dAdn**3    + dHdttt*dtdn**3 &
             + 3.d0*dHdaat*dAdn**2*dtdn  + 3.d0*dHdatt*dAdn*dtdn**2 &
             + 3.d0*dHdaa*dAdn*dAdnn + 3.d0*dHdtt*dtdn*dtdnn + 3.d0*dHdat*(dAdnn*dtdn + dAdn*dtdnn) &
             + dHda*dAdnnn   + dHdt*dtdnnn
        dHdnng   = dHdaat*dAdn**2*dtdg   + 2.d0*dHdatt*dAdn*dtdn*dtdg    + dHdttt*dtdn**2*dtdg &
             + dHdat*(dAdnn*dtdg + 2.d0*dAdn*dtdng)  + dHdtt*(dtdnn*dtdg + 2.d0*dtdn*dtdng) &
             + dHdt*dtdnng
        dHdngg   = dHdatt*dAdn*dtdg**2   + dHdttt*dtdn*dtdg**2   + 2.d0*dHdtt*dtdg*dtdng
        dHdggg   = dHdttt*dtdg**3
        !dF_drho(i, 1)=h+d*dHdn+ecd
        dFc_daa(i)   = 2.d0*dHdn + d*dHdnn       + ecdd
        dFc_dag(i)   = dHdg  + d*dHdng
        dFc_dgg(i)   = d*dHdgg
        dFc_daaa(i)  = 3.d0*dHdnn    + d*dHdnnn  + ecddd
        dFc_daag(i)  = 2.d0*dHdng    + d*dHdnng
        dFc_dagg(i)  = dHdgg         + d*dHdngg
        dFc_dggg(i)  = d*dHdggg  
     end do

  else if ( nspin ==  2 ) then
     thrd2 = thrd * 2.d0
     thrd4 = thrd * 4.d0
     fzdd = 8.d0 /(9.d0*(2.d0**thrd4 - 2.d0 ))

     do i = ist,ien,dnr ! MPI
        d = chgrhr_l(i, 1) + chgrhr_l(i, nspin)
        if(d < density_minimum .or. chgrhr_l(i,1) < density_minimum .or. chgrhr_l(i,nspin)<density_minimum) cycle

        zeta = ( chgrhr_l(i, 1) - chgrhr_l(i, 2) ) / d
        onpzeta = 2.d0 *chgrhr_l(i, 1) / d
        onmzeta = 2.d0 *chgrhr_l(i, 2) / d
        g = ( onpzeta**thrd2 + onmzeta**thrd2 )/2
        if(dabs(zeta) < zeta_minimum ) then
           fz= 4.d0/9 *zeta*zeta*(1.d0 + 5.d0/54.d0 *zeta*zeta &
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
           ! ****** aas mkatsu *****
           !              t = dd/(d*sk*2)
           t = dd/(d*sk*2.d0)/g         !    mkatsu aas
           s = dd/(d*fk*2.d0)
        else
           t = 0.d0
           s = 0.d0
        endif

        q0 = -2.d0*a*( 1.d0 +a1*rs )
        rs12 = dsqrt(rs)
        rs32 = rs12**3
        rsp = rs**p
        q1 = 2.d0 *a*( b1*rs12 +b2n*rs +b3*rs32 +b4*rs*rsp )
        q2 = log( 1.d0 +1.d0/q1 )
        q3 = a*( b1/rs12 +2.d0*b2n +3.d0*b3*rs12 +2.d0*b4*p1*rsp )

        q0p = -2.d0 *ap *( 1.d0 +a1p*rs )
        q1p = 2.d0 *ap*( b1p*rs12 +b2np*rs+b3p*rs32+b4p*rs*rsp )
        q2p = log( 1.d0 +1.d0/q1p )
        q3p = ap*( b1p/rs12 +2.d0*b2np +3.d0*b3p*rs12 +2.d0*b4p*p1*rsp )
        ! vwn 95/12/2 Y.M
        q0q = -2.d0 *aq*( 1.d0 +a1q*rs )
        q1q =  2.d0 *aq*( b1q*rs12 +b2nq*rs +b3q*rs32 +b4q*rs*rsp )
        q2q = log( 1.d0 +1.d0/q1q )
        q3q = aq*( b1q/rs12 +2.d0*b2nq +3.d0*b3q*rs12 +2.d0*b4q*p1*rsp )

        if(dabs(zeta) < zeta_minimum) then
           fzd = 8.d0/9.d0 /( 2.d0**thrd4 -2.d0 ) &
                &     *zeta *( 1.d0 + 5.d0/27.d0 *zeta*zeta &
                &           *( 1.d0 +22.d0/45.d0 *zeta*zeta ))
        else
           fzd = 4.d0/3.d0 /( 2.d0**thrd4 -2.d0 ) &
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

        !           zetadxd = - 2 * chgrhr_l(i, 3-is) * (-1.d0)**is / d
        zetadxda = - 2.d0 * chgrhr_l(i, 2) * (-1.d0) / d
        zetadxdb = - 2.d0 * chgrhr_l(i, 1) / d
        !           ecd = eu-thrd*rs*eurs + euzt * zetadxd
        ecda = eu-thrd*rs*eurs + euzt * zetadxda 
        ecdb = eu-thrd*rs*eurs + euzt * zetadxdb

        pon = facpon*eu
        b   = delt/(dexp(pon)-1.d0)
        b2  = b*b
        t2  = t*t
!!$           s2  = s*s
        t4  = t2*t2
        t6  = t4*t2
        q4  = 1.d0 +b*t2
        q5  = 1.d0 +b*t2+b2*t4
        h0  = g3*(bet/delt)*dlog( 1.d0 +delt*q4*t2/q5 )
        h   = h0

        q8  = q5*q5+delt*q4*q5*t2
        h0t = 2.d0 *bet*t*( 1.d0 +2.d0 *b*t2 )/q8 * g3
        h0b = -bet*t6*( 2.d0 *b+b2*t2 )/q8 * g3
        h0rs = h0b * b * eurs * (b+delt) / bet / g3
        ht   = h0t
        hrs  = h0rs

        if(dabs(zeta) < zeta_minimum) then
           gzt = -2.d0/9*zeta *( 1.d0 +14.d0/27 *zeta*zeta &
                &             *( 1.d0 +13.d0/18 *zeta*zeta ))
        else
           gzt = ( (onpzeta**2+eta)**sixthm - (onmzeta**2+eta)**sixthm)/3.d0
        end if

        bzt = b * ( b+delt) * ( euzt - 3*eu*gzt/g ) / ( bet * g3 ) 
        h0zt = 3.d0 * gzt * h0 / g + h0b * bzt

        ec1  = d*h * 0.5d0

        ! --> revised by Kato-san,  1999Jul16
!!$           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + h0zt * zetadxd
        !           ec1d = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxd
        ec1da = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxda
        ec1db = h-thrd*rs*hrs-sixth7*t*ht  + (h0zt-gzt*t*h0t/g) * zetadxdb
        ec1dd = 0.5d0*ht/sk / g

10      continue


        exc0 = eu*d * 0.5d0 + ec1
        !           excd = ecd + ec1d
        excda= ecda + ec1da
        excdb= ecdb + ec1db
        !           dF_drho(i, is) = dF_drho(i, is) + excd
        dF_drho(i, 1) = dF_drho(i, 1) + excda
        dF_drho(i, 2) = dF_drho(i, 2) + excdb
        excdd = ec1dd
        ! gradient of charge density 95/12/3 H.S.
        !           if(is ==  2) then
        if(dabs(grad_trho(i)) > 1.0d-9) then
           !                 dF_dgradrho(i) = excdd / grad_trho(i)
           dF_dgradrho(i) = excdd
        else
           dF_dgradrho(i) = 0.d0
        endif
        !           endif
        ! gradient of charge density 95/12/3 H.S.
!!$           exec1 = eu*d * 0.5
        !           exc1 = exc1 + exc0*wos(i)
        exc(i)=exc(i)+exc0*2.d0

        ! ***** Calc heigher order derivatives. *****

        rsdn=-thrd*rs/d
        rsdnn=-4.d0*thrd*rsdn/d
        rsdnnn=-thrd7*rsdnn/d

        id2=1.d0/d/d
        id3=id2/d
        id4=id3/d
        zetada   =     zetadxda/d
        zetadb   =     zetadxdb/d
        zetadaa  =    -4.d0*chgrhr_l(i, 2)*id3
        zetadab  =     2.d0*zeta*id2
        zetadbb  =     4.d0*chgrhr_l(i, 1)*id3
        zetadaaa =    12.d0*chgrhr_l(i, 2)*id4
        zetadaab =   ( -4.d0*chgrhr_l(i, 1) + 8.d0*chgrhr_l(i, 2))*id4
        zetadabb =   ( -8.d0*chgrhr_l(i, 1) + 4.d0*chgrhr_l(i, 2))*id4
        zetadbbb =   -12.d0*chgrhr_l(i, 1)*id4

        !print '(i2,7d19.10)', i,chgrhr_l(i,1:2),zetadab,zetadaaa,zetadaab,zetadabb,zetadbbb

        q1drr =0.5d0*a*( -b1/rs32+3.d0*b3/rs12+8.d0*b4 )
        q1drrr=0.75d0*a*( b1/rs - b3 )/rs32
        invq12pq1=1.d0/(q1**2 + q1)
        q12p1=2.d0*q1+1.d0
        e1dr=-2.d0*a*a1*q2-q0*q3*invq12pq1
        e1drr=4.d0*a*a1*q3*invq12pq1-q0*( q1drr-q12p1*q3**2*invq12pq1 )*invq12pq1
        e1drrr=6.d0*a*a1*( q1drr - q12p1*q3**2*invq12pq1 ) * invq12pq1 &
             -q0*( q1drrr -(2.d0*q3**3+3.d0*q12p1*q3*q1drr -  &
             2.d0*q12p1**2*q3**3* invq12pq1)*invq12pq1 ) * invq12pq1

        q1pdrr =0.5d0*ap*( -b1p/rs32+3.d0*b3p/rs12+8.d0*b4p )
        q1pdrrr=0.75d0*ap*( b1p/rs - b3p )/rs32
        invq1p2pq1p=1.d0/(q1p**2 + q1p)
        q1p2p1=2.d0*q1p+1.d0
        e1pdr=-2.d0*ap*a1p*q2p-q0p*q3p*invq1p2pq1p
        e1pdrr=4.d0*ap*a1p*q3p*invq1p2pq1p-q0p*( q1pdrr-q1p2p1*q3p**2*invq1p2pq1p )*invq1p2pq1p
        e1pdrrr=6.d0*ap*a1p*( q1pdrr - q1p2p1*q3p**2*invq1p2pq1p ) * invq1p2pq1p &
             -q0p*( q1pdrrr -(2.d0*q3p**3+3.d0*q1p2p1*q3p*q1pdrr -  &
             2.d0*q1p2p1**2*q3p**3* invq1p2pq1p)*invq1p2pq1p ) * invq1p2pq1p

        q1qdrr =0.5d0*aq*( -b1q/rs32+3.d0*b3q/rs12+8.d0*b4q )
        q1qdrrr=0.75d0*aq*( b1q/rs - b3q )/rs32
        invq1q2pq1q=1.d0/(q1q**2 + q1q)
        q1q2p1=2.d0*q1q+1.d0
        e1qdr=-2.d0*aq*a1q*q2q-q0q*q3q*invq1q2pq1q
        e1qdrr=4.d0*aq*a1q*q3q*invq1q2pq1q-q0q*( q1qdrr-q1q2p1*q3q**2*invq1q2pq1q )*invq1q2pq1q
        e1qdrrr=6.d0*aq*a1q*( q1qdrr - q1q2p1*q3q**2*invq1q2pq1q ) * invq1q2pq1q &
             -q0q*( q1qdrrr -(2.d0*q3q**3+3.d0*q1q2p1*q3q*q1qdrr -  &
             2.d0*q1q2p1**2*q3q**3* invq1q2pq1q)*invq1q2pq1q ) * invq1q2pq1q

        zeta2=zeta**2
        zeta3=zeta2*zeta
        zeta4=zeta2*zeta2

        if(dabs(zeta) < zeta_minimum) then
           fzdd0 = fzdd &
                &            *( 1 + 5.d0*ninth *zeta2 &
                &            *( 1 +22.d0*twntysvnth *zeta2 ))
           fzddd = fzdd*10.d0*ninth*zeta*( 1.d0 + 44.d0*twntysvnth*zeta2 )
           gzz   = -2.d0*ninth*( 1.d0 + 14.d0*ninth*zeta2*( 1.d0 + 65.d0*fftyfrth*zeta2 ) )
           gzzz  = -56.d0*eityonth*zeta*( 1.d0 + 65.d0*twntysvnth*zeta2 )
        else
           opzthrdm2=1.d0/onpzeta**(thrd2)
           omzthrdm2=1.d0/onmzeta**(thrd2)
           fzdd0 =  fzdd*0.5d0*( opzthrdm2 + omzthrdm2 )
           fzddd = -fzdd*thrd*( opzthrdm2/onpzeta - omzthrdm2/onmzeta )
           gzz   = -ninth*( opzthrdm2**2 + omzthrdm2**2 )
           !              gzz   = -ninth*( 1.d0/onpzeta**(thrd4) + 1.d0/onmzeta**(thrd4) )
           gzzz  = 4.d0*twntysvnth*( opzthrdm2**2/onpzeta - omzthrdm2**2/onmzeta )
        end if

        e1=q0*q2
        e2=q0p*q2p
        e3=q0q*q2q

        eudrr    = e1drr  + ( e1pdrr  -  e1drr  ) * fz * zeta4 - e1qdrr  * fz * onzt4 / fzdd
        eudrrr   = e1drrr + ( e1pdrrr -  e1drrr ) * fz * zeta4 - e1qdrrr * fz * onzt4 / fzdd
        eudzz    = ( e2 - e1 )*( fzdd0*zeta4 + 8.d0*fzd*zeta3 + 12.d0*fz*zeta2 ) &
             - e3*( fzdd0*onzt4 - 8.d0*fzd*zeta3 - 12.d0*fz*zeta2 )/fzdd
        eudzzz   = ( e2 - e1 )*( fzddd*zeta4 + 12.d0*fzdd0*zeta3 + 36.d0*fzd*zeta2 + 24.d0*fz*zeta ) &
             - e3*( fzddd*onzt4 - 12.d0*fzdd0*zeta3 - 36.d0*fzd*zeta2 - 24.d0*fz*zeta )/fzdd
        eudrz    = ( e1pdr - e1dr )*( fzd*zeta4 + 4.d0*fz*zeta3 ) &
             - e1qdr *( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd 
        eudrrz   = ( e1pdrr  - e1drr  )*( fzd*zeta4 + 4.d0*fz*zeta3 ) &
             - e1qdrr *( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd
        eudrzz   = ( e1pdr - e1dr )*( fzdd0*zeta4 + 8.d0*fzd*zeta3 + 12.d0*fz*zeta2 ) &
             - e1qdr *( fzdd0*onzt4 - 8.d0*fzd*zeta3 - 12.d0*fz*zeta2 )/fzdd

        rsdn2    = rsdn*rsdn
        rsdn3    = rsdn2*rsdn
        zetada2  = zetada**2
        zetada3  = zetada2*zetada
        zetadb2  = zetadb**2
        zetadb3  = zetadb2*zetadb

        euda     = eurs*rsdn + euzt*zetada
        eudb     = eurs*rsdn + euzt*zetadb
        eudaa    = eudrr*rsdn2 + 2.d0*eudrz*rsdn*zetada + eudzz*zetada2 + eurs*rsdnn + euzt*zetadaa
        eudab    = eudrr*rsdn2 + eudrz*rsdn*( zetada + zetadb ) + eudzz*zetada*zetadb &
             + eurs*rsdnn + euzt*zetadab
        eudbb    = eudrr*rsdn2 + 2.d0*eudrz*rsdn*zetadb + eudzz*zetadb2 + eurs*rsdnn + euzt*zetadbb
        eudaaa   = eudrrr*rsdn3 &
             + 3.d0*eudrrz*rsdn2*zetada &
             + 3.d0*eudrzz*rsdn*zetada2 + eudzzz*zetada3 &
             + 3.d0*eudrr*rsdn*rsdnn + 3.d0*eudzz*zetada*zetadaa &
             + 3.d0*eudrz*( rsdnn*zetada + rsdn*zetadaa ) &
             + eurs*rsdnnn + euzt*zetadaaa
        eudaab   = eudrrr*rsdn3 &
             + eudrrz*rsdn2*( zetadb + 2.d0*zetada ) &
             + eudrzz*rsdn*zetada*( zetada + 2.d0*zetadb ) &
             + eudzzz*zetada2*zetadb &
             + 3.d0*eudrr*rsdn*rsdnn &
             + eudrz*( rsdnn*( 2.d0*zetada + zetadb ) + rsdn*( 2.d0*zetadab + zetadaa ) ) &
             + eudzz*( 2.d0*zetada*zetadab + zetadaa*zetadb ) &
             + eurs*rsdnnn + euzt*zetadaab
        eudabb   = eudrrr*rsdn3 &
             + eudrrz*rsdn2*( zetada + 2.d0*zetadb ) &
             + eudrzz*rsdn*zetadb*( zetadb + 2.d0*zetada ) &
             + eudzzz*zetada*zetadb2 &
             + 3.d0*eudrr*rsdn*rsdnn &
             + eudrz*( rsdnn*( zetada + 2.d0*zetadb ) + rsdn*( 2.d0*zetadab + zetadbb ) ) &
             + eudzz*( zetada*zetadbb + 2.d0*zetadb*zetadab ) &
             + eurs*rsdnnn + euzt*zetadabb
        eudbbb   = eudrrr*rsdn3 &
             + 3.d0*eudrrz*rsdn2*zetadb &
             + 3.d0*eudrzz*rsdn*zetadb2 + eudzzz*zetadb3 &
             + 3.d0*eudrr*rsdn*rsdnn + 3.d0*eudzz*zetadb*zetadbb &
             + 3.d0*eudrz*( rsdnn*zetadb + rsdn*zetadbb ) &
             + eurs*rsdnnn + euzt*zetadbbb

        !print '(i2,12d19.10)', i,chgrhr_l(i,1:2),eu,euda,eudb,eudaa,eudbb,eudab,eudaaa,eudaab,eudabb,eudbbb
        !print *,eu,     e1+(e2-e1)* fz * zeta4 - e3  * fz * onzt4 / fzdd
        !print *,eurs,   e1dr  + ( e1pdr  -  e1dr  ) * fz * zeta4 - e1qdr  * fz * onzt4 / fzdd
        !print *,euzt,   (e2-e1)*( fzd*zeta4 + 4.d0*fz*zeta3 )-e3*( fzd*onzt4 - 4.d0*fz*zeta3 )/fzdd 

        duda     = gzt*zetada
        dudb     = gzt*zetadb 
        dudaa    = gzz*zetada2     + gzt*zetadaa
        dudab    = gzz*zetada*zetadb + gzt*zetadab
        dudbb    = gzz*zetadb2     + gzt*zetadbb
        dudaaa   = gzzz*zetada3    + 3.d0*gzz*zetada*zetadaa   + gzt*zetadaaa
        dudaab   = gzzz*zetada2*zetadb + gzz*( zetadaa*zetadb + 2.d0*zetada*zetadab ) + gzt*zetadaab
        dudabb   = gzzz*zetada*zetadb2 + gzz*( zetada*zetadbb + 2.d0*zetadab*zetadb ) + gzt*zetadabb
        dudbbb   = gzzz*zetadb3    + 3.d0*gzz*zetadb*zetadbb   + gzt*zetadbbb 

        !print '(i2,12d19.10)', i,chgrhr_l(i,1:2),g,duda,dudb,dudaa,dudbb,dudab,dudaaa,dudaab,dudabb,dudbbb


        dtdn  = -sixth7*t/d
        dtdnn = -sixth13*dtdn/d
        dtdnnn= -sixth19*dtdnn/d
        !print '(i2,5e19.12)',i,d,t,dtdn,dtdnn,dtdnnn
        if(dabs(grad_trho(i)) > 1.0d-9) then
           dtdg  = t/dd
           dtdng = -sixth7*dtdg/d
           dtdnng= -sixth13*dtdng/d
           !               dtdnug = -dtdn*invg/dd
        else
           dtdg  = 0.d0
           dtdng = 0.d0
           dtdnng= 0.d0
           dtdnug= 0.d0 
        endif
        dtdgg  = 0.d0
        dtdngg = 0.d0
        dtdggg = 0.d0
        invg=1.d0/g
        dtdu   = -t*invg
        dtduu  = -2.d0*dtdu*invg 
        dtduuu = -3.d0*dtduu*invg
        dtdnu  = -dtdn*invg
        dtdnnu = -dtdnn*invg
        dtdnuu = -sixth7*dtduu/d
        dtdug  = -dtdg*invg
        dtduug = 2.d0*dtdg*invg**2
        dtdugg = 0.d0
        if(dabs(grad_trho(i)) > 1.0d-9) then
           dtdnug = dtdnu/dd
        else
           dtdnug = 0.d0
        end if

        dtda   = dtdn + dtdu*duda
        dtdb   = dtdn + dtdu*dudb
        dtdaa  = dtdnn + 2.d0*dtdnu*duda + dtduu*duda**2 + dtdu*dudaa
        dtdbb  = dtdnn + 2.d0*dtdnu*dudb + dtduu*dudb**2 + dtdu*dudbb
        dtdab  = dtdnn + dtdnu*( duda + dudb ) + dtduu*duda*dudb + dtdu*dudab
        dtdag  = dtdng + dtdug*duda
        dtdbg  = dtdng + dtdug*dudb
        dtdaaa = dtdnnn+ 3.d0*dtdnnu*duda + 3.d0*dtdnuu*duda**2 + dtduuu*duda**3 &
             + 3.d0*dtduu*duda*dudaa + 3.d0*dtdnu*dudaa + dtdu*dudaaa
        dtdbbb = dtdnnn+ 3.d0*dtdnnu*dudb + 3.d0*dtdnuu*dudb**2 + dtduuu*dudb**3 &
             + 3.d0*dtduu*dudb*dudbb + 3.d0*dtdnu*dudbb + dtdu*dudbbb
        dtdaab = dtdnnn+ dtdnnu*( 2.d0*duda + dudb ) + dtdnuu*duda*( duda + 2.d0*dudb ) + dtduuu*duda**2*dudb &
             + dtduu*( dudaa*dudb + 2.d0*duda*dudab ) + dtdnu*( dudaa + 2.d0*dudab ) + dtdu*dudaab
        dtdabb = dtdnnn+ dtdnnu*( duda + 2.d0*dudb ) + dtdnuu*dudb*( dudb + 2.d0*duda ) + dtduuu*dudb**2*duda &
             + dtduu*( dudbb*duda + 2.d0*dudb*dudab ) + dtdnu*( dudbb + 2.d0*dudab ) + dtdu*dudabb 
        dtdaag = dtdnng+ 2.d0*dtdnug*duda + dtduug*duda**2 + dtdug*dudaa
        dtdbbg = dtdnng+ 2.d0*dtdnug*dudb + dtduug*dudb**2 + dtdug*dudbb
        dtdabg = dtdnng+ dtdnug*( duda + dudb ) + dtduug*duda*dudb + dtdug*dudab 

        !print '(i2,20d19.10)', i,t,dtda,dtdb,dtdg,dtdaa,dtdbb,dtdgg,dtdab,dtdag,dtdbg,dtdaaa,dtdbbb,dtdggg,dtdaab,dtdabb,dtdaag,dtdagg,dtdbbg,dtdbgg,dtdabg     
        !print '(i2,12d19.10)', i,t,dtdag,dtdbg,dtdaag,dtdbbg,dtdabg     

        invg3  = 1.d0  /g3
        dAde   = b * ( b + delt ) / bet * invg3
        dAdee  = 2.d0 * dAde * ( b + 0.5d0 * delt ) / bet *invg3
        dAdeee = 2.d0 * ( dAde**2 + ( b + 0.5d0 * delt )*dAdee ) /bet * invg3
        m3u    = -3.d0*invg
        dAdu     = m3u*eu*dAde
        dAdeu    = m3u*( dAde + eu*dAdee )
        dAduu    = m3u*eu*invg*( -dAde + g*dAdeu )
        dAdeeu   = m3u*( 2.d0*dAdee + eu*dAdeee )
        dAdeuu   = m3u*invg*( -dAde + g*dAdeu - eu*dAdee + g*eu*dAdeeu )
        dAduuu   = m3u*eu*invg**2*( 2.d0*dAde -2.d0*g*dAdeu + g**2*dAdeuu )

        dAda     = dAde*euda + dAdu*duda
        dAdb     = dAde*eudb + dAdu*dudb
        dAdaa    = dAdee*euda**2 + 2.d0*dAdeu*euda*duda + dAduu*duda**2 + dAde*eudaa + dAdu*dudaa 
        dAdab    = dAdee*euda*eudb + dAdeu*( euda*dudb + eudb*duda ) + dAduu*duda*dudb + dAde*eudab + dAdu*dudab 
        dAdbb    = dAdee*eudb**2 + 2.d0*dAdeu*eudb*dudb + dAduu*dudb**2 + dAde*eudbb + dAdu*dudbb
        dAdaaa   = dAdeee*euda**3 + 3.d0*dAdeeu*euda**2*duda + 3.d0*dAdeuu*euda*duda**2 + dAduuu*duda**3 &
             + 3.d0*dAdee*euda*eudaa + 3.d0*dAduu*duda*dudaa &
             + 3.d0*dAdeu*( euda*dudaa + eudaa*duda ) + dAde*eudaaa + dAdu*dudaaa
        dAdaab   = dAdeee*euda**2*eudb + dAduuu*duda**2*dudb &
             + dAdeeu*euda*( euda*dudb + 2.d0*eudb*duda ) &
             + dAdeuu*( 2.d0*euda*dudb + eudb*duda )*duda &
             + dAdee*( 2.d0*euda*eudab + eudb*eudaa ) + dAduu*( 2.d0*duda*dudab + dudb*dudaa ) &
             + dAdeu*( 2.d0*eudab*duda + 2.d0*euda*dudab + dudb*eudaa + eudb*dudaa ) &
             + dAde*eudaab + dAdu*dudaab
        dAdabb   = dAdeee*euda*eudb**2 + dAduuu*dudb**2*duda &
             + dAdeeu*eudb*( eudb*duda + 2.d0*euda*dudb ) &
             + dAdeuu*( 2.d0*eudb*duda + euda*dudb )*dudb &
             + dAdee*( eudbb*euda + 2.d0*eudb*eudab ) + dAduu*( dudbb*duda + 2.d0*dudb*dudab ) &
             + dAdeu*( eudbb*duda + 2.d0*eudb*dudab + 2.d0*dudb*eudab + dudbb*euda ) &
             + dAde*eudabb + dAdu*dudabb
        dAdbbb   = dAdeee*eudb**3 + 3.d0*dAdeeu*eudb**2*dudb + 3.d0*dAdeuu*eudb*dudb**2 + dAduuu*dudb**3 &
             + 3.d0*dAdee*eudb*eudbb + 3.d0*dAduu*dudb*dudbb &
             + 3.d0*dAdeu*( eudb*dudbb + eudbb*dudb ) + dAde*eudbbb + dAdu*dudbbb 

        !print '(i2,12d19.10)', i,b,dAda,dAdb,dAdaa,dAdab,dAdbb,dAdaaa,dAdaab,dAdabb,dAdbbb

        at2  = b * t2
        a2t4 = at2**2
        a3t6 = at2*a2t4
        t3   = t2*t
        t5   = t4*t
        t7   = t6*t
        t8   = t7*t
        fE    = q5**2 + delt*t2*q4*q5
        fEdt  = 4.d0*b*t*q5*( 1.d0 + 2.d0*at2 ) + &
             delt *2.d0*t*(1.d0 + 4.d0*at2 + 6.d0*a2t4 + 4.d0*a3t6 )
        fEdtt = 4.d0*b*( 1.d0 + 9.d0*at2 + 15.d0*a2t4 + 14.d0*a3t6 ) + &
             delt*( 2.d0 + 24.d0*at2 + 60.d0*a2t4 + 56.d0*a3t6 )
        fEda  = 2.d0*t2*q5*( 1.d0 + 2.d0*at2 ) + &
             delt*t4*( 2.d0 + 4.d0*at2 + 3.d0*a2t4 )
        fEdaa = 6.d0*t4*( 2.d0*q5 - 1.d0) + &
             delt*2.d0*t6*( 2.d0 + 3.d0*at2 )
        fF    = t*( 1.d0 + 2.d0*at2 )         
        fFdt  = 1.d0 + 6.d0*at2
        fFdtt = 12.d0*b*t
        fFda  = 2.d0*t3
        fFdaa = 0.d0
        fG    = 2.d0*b*t6 + b2*t8
        fGdt  = 12.d0*b*t5 + 8.d0*b2*t7
        fGdtt = 60.d0*b*t4 + 56.d0*b2*t6
        fGda  = 2.d0*t6 + 2.d0*b*t8
        fGdaa = 2.d0*t8

        betivE=bet/fE*g3
        betivE2=betivE/fE
        dHdf   =     -fG*betivE
        dHdt   = 2.d0*fF*betivE
        dHdff  =     -( fGda*fE-fG*fEda )*betivE2
        dHdft  =     -( fGdt*fE-fG*fEdt )*betivE2
        dHdtt  = 2.d0*( fFdt*fE-fF*fEdt )*betivE2
        dHdfff =     -( fGdaa*fE-fG*fEdaa-2.d0*( fGda*fE-fG*fEda )*fEda/fE )*betivE2
        dHdfft = 2.d0*( fFdaa*fE-fF*fEdaa-2.d0*( fFda*fE-fF*fEda )*fEda/fE )*betivE2
        dHdftt =     -( fGdtt*fE-fG*fEdtt-2.d0*( fGdt*fE-fG*fEdt )*fEdt/fE )*betivE2
        dHdttt = 2.d0*( fFdtt*fE-fF*fEdtt-2.d0*( fFdt*fE-fF*fEdt )*fEdt/fE )*betivE2

        tiu    = 3.d0*invg
        siu2   = tiu*2.d0*invg
        dHdu   = tiu*h
        dHduu  = siu2*h
        dHduuu = siu2*invg*h
        dHdfu  = tiu*dHdf
        dHdtu  = tiu*dHdt
        dHdffu = tiu*dHdff
        dHdfuu = siu2*dHdf
        dHdttu = tiu*dHdtt
        dHdtuu = siu2*dHdt
        dHdftu = tiu*dHdft

        dHda   = dHdf*dAda + dHdt*dtda + dHdu*duda
        dHdb   = dHdf*dAdb + dHdt*dtdb + dHdu*dudb
        dHdg   =             dHdt*dtdg
        dHdaa  = dHdff*dAda**2 + dHdtt*dtda**2 + dHduu*duda**2 &
             + 2.d0*dHdft*dAda*dtda + 2.d0*dHdfu*dAda*duda + 2.d0*dHdtu*dtda*duda &
             + dHdf*dAdaa + dHdt*dtdaa + dHdu*dudaa
        dHdab  = dHdff*dAda*dAdb + dHdtt*dtda*dtdb + dHduu*duda*dudb &
             + dHdft*( dAda*dtdb + dAdb*dtda ) + dHdfu*( dAda*dudb + dAdb*duda ) + dHdtu*( dtda*dudb + dtdb*duda ) &
             + dHdf*dAdab + dHdt*dtdab + dHdu*dudab
        dHdbb  = dHdff*dAdb**2 + dHdtt*dtdb**2 + dHduu*dudb**2 &
             + 2.d0*dHdft*dAdb*dtdb + 2.d0*dHdfu*dAdb*dudb + 2.d0*dHdtu*dtdb*dudb &
             + dHdf*dAdbb + dHdt*dtdbb + dHdu*dudbb
        dHdgg  = dHdtt*dtdg**2
        dHdag  = dHdft*dAda*dtdg + dHdtt*dtda*dtdg + dHdtu*duda*dtdg + dHdt*dtdag
        dHdbg  = dHdft*dAdb*dtdg + dHdtt*dtdb*dtdg + dHdtu*dudb*dtdg + dHdt*dtdbg
        dHdaaa = dHdfff*dAda**3 + 3.d0*dHdfft*dAda**2*dtda + 3.d0*dHdffu*dAda**2*duda &
             + 3.d0*dHdftt*dAda*dtda**2 + dHdttt*dtda**3 + 3.d0*dHdttu*dtda**2*duda &
             + 3.d0*dHdfuu*dAda*duda**2 + 3.d0*dHdtuu*dtda*duda**2 + dHduuu*duda**3 &
             + 6.d0*dHdftu*dAda*dtda*duda &
             + 3.d0*dHdff*dAda*dAdaa + 3.d0*dHdtt*dtda*dtdaa + 3.d0*dHduu*duda*dudaa &
             + 3.d0*dHdft*( dAdaa*dtda + dAda*dtdaa ) &
             + 3.d0*dHdtu*( dtdaa*duda + dtda*dudaa ) &
             + 3.d0*dHdfu*( dAdaa*duda + dAda*dudaa ) &
             + dHdf*dAdaaa + dHdt*dtdaaa + dHdu*dudaaa
        dHdaab = dHdfff*dAda**2*dAdb + dHdfft*dAda*( dAda*dtdb + 2.d0*dAdb*dtda ) &
             + dHdffu*dAda*( dAda*dudb + 2.d0*dAdb*duda ) &
             + dHdftt*dtda*( dAdb*dtda + 2.d0*dAda*dtdb ) &
             + dHdttt*dtdb*dtda**2 + dHdttu*dtda*( dtda*dudb + 2.d0*dtdb*duda ) &
             + dHdfuu*duda*( dAdb*duda + 2.d0*dAda*dudb ) &
             + dHdtuu*duda*( dtdb*duda + 2.d0*dtda*dudb ) + dHduuu*duda**2*dudb &
             + 2.d0*dHdftu*( dAdb*duda*dtda + dAda*dudb*dtda + dAda*duda*dtdb ) &
             + dHdff*( 2.d0*dAda*dAdab + dAdb*dAdaa ) &
             + dHdtt*( 2.d0*dtda*dtdab + dtdb*dtdaa ) &
             + dHduu*( 2.d0*duda*dudab + dudb*dudaa ) &
             + dHdft*( 2.d0*dAdab*dtda + 2.d0*dAda*dtdab &
             + dAdaa*dtdb + dAdb*dtdaa ) &
             + dHdtu*( 2.d0*dtdab*duda + 2.d0*dtda*dudab &
             + dtdaa*dudb + dtdb*dudaa ) &
             + dHdfu*( 2.d0*dAdab*duda + 2.d0*dAda*dudab &
             + dAdaa*dudb + dAdb*dudaa ) &
             + dHdf*dAdaab + dHdt*dtdaab + dHdu*dudaab
        dHdabb = dHdfff*dAda*dAdb**2 + dHdfft*dAdb*( dAdb*dtda + 2.d0*dAda*dtdb ) &
             + dHdffu*dAdb*( dAdb*duda + 2.d0*dAda*dudb ) &
             + dHdftt*dtdb*( dAda*dtdb + 2.d0*dAdb*dtda ) + dHdttt*dtda*dtdb**2 &
             + dHdttu*dtdb*( dtdb*duda + 2.d0*dtda*dudb ) &
             + dHdfuu*dudb*( dAda*dudb + 2.d0*dAdb*duda ) &
             + dHdtuu*dudb*( dtda*dudb + 2.d0*dtdb*duda ) + dHduuu*duda*dudb**2 &
             + 2.d0*dHdftu*( dAda*dtdb*dudb + dAdb*dtda*dudb + dAdb*dtdb*duda ) &
             + dHdff*( 2.d0*dAdb*dAdab + dAda*dAdbb ) &
             + dHdtt*( 2.d0*dtdb*dtdab + dtda*dtdbb ) &
             + dHduu*( 2.d0*dudb*dudab + duda*dudbb ) &
             + dHdft*( 2.d0*dAdab*dtdb + 2.d0*dAdb*dtdab + dAdbb*dtda + dAda*dtdbb ) &
             + dHdtu*( 2.d0*dtdab*dudb + 2.d0*dtdb*dudab + dtdbb*duda + dtda*dudbb ) &
             + dHdfu*( 2.d0*dAdab*dudb + 2.d0*dAdb*dudab + dAdbb*duda + dAda*dudbb ) &
             + dHdf*dAdabb + dHdt*dtdabb + dHdu*dudabb
        dHdbbb = dHdfff*dAdb**3 + 3.d0*dHdfft*dAdb**2*dtdb + 3.d0*dHdffu*dAdb**2*dudb &
             + 3.d0*dHdftt*dAdb*dtdb**2 + dHdttt*dtdb**3 + 3.d0*dHdttu*dtdb**2*dudb &
             + 3.d0*dHdfuu*dAdb*dudb**2 + 3.d0*dHdtuu*dtdb*dudb**2 + dHduuu*dudb**3 &
             + 6.d0*dHdftu*dAdb*dtdb*dudb &
             + 3.d0*dHdff*dAdb*dAdbb + 3.d0*dHdtt*dtdb*dtdbb + 3.d0*dHduu*dudb*dudbb &
             + 3.d0*dHdft*( dAdbb*dtdb + dAdb*dtdbb ) &
             + 3.d0*dHdtu*( dtdbb*dudb + dtdb*dudbb ) &
             + 3.d0*dHdfu*( dAdbb*dudb + dAdb*dudbb ) &
             + dHdf*dAdbbb + dHdt*dtdbbb + dHdu*dudbbb
        dHdabg = dHdfft*dAda*dAdb*dtdg   + dHdftt*( dAdb*dtda + dAda*dtdb )*dtdg &
             + dHdttt*dtda*dtdb*dtdg + dHdttu*( dtdb*duda + dtda*dudb )*dtdg &
             + dHdtuu*duda*dudb*dtdg + dHdftu*( dAdb*duda + dAda*dudb )*dtdg &
             + dHdft*( dAdab*dtdg + dAdb*dtdag + dAda*dtdbg ) &
             + dHdtt*( dtdab*dtdg + dtdb*dtdag + dtda*dtdbg ) &
             + dHdtu*( dudab*dtdg + dudb*dtdag + duda*dtdbg ) &
             + dHdt*dtdabg
        dHdaag = dHdfft*dAda**2*dtdg + dHdttt*dtda**2*dtdg + dHdtuu*duda**2*dtdg &
             + 2.d0*dHdttu*dtda*duda*dtdg + 2.d0*dHdftt*dAda*dtda*dtdg + 2.d0*dHdftu*dAda*duda*dtdg &
             + dHdtt*( dtdaa*dtdg + 2.d0*dtda*dtdag ) + dHdft*( dAdaa*dtdg + 2.d0*dAda*dtdag ) &
             + dHdtu*( dudaa*dtdg + 2.d0*duda*dtdag ) + dHdt*dtdaag         
        dHdbbg = dHdfft*dAdb**2*dtdg + dHdttt*dtdb**2*dtdg + dHdtuu*dudb**2*dtdg &
             + 2.d0*dHdttu*dtdb*dudb*dtdg + 2.d0*dHdftt*dAdb*dtdb*dtdg + 2.d0*dHdftu*dAdb*dudb*dtdg &
             + dHdtt*( dtdbb*dtdg + 2.d0*dtdb*dtdbg ) + dHdft*( dAdbb*dtdg + 2.d0*dAdb*dtdbg ) &
             + dHdtu*( dudbb*dtdg + 2.d0*dudb*dtdbg ) + dHdt*dtdbbg       
        dHdagg = dHdftt*dAda*dtdg**2 + dHdttt*dtda*dtdg**2 + dHdttu*duda*dtdg**2 + 2.d0*dHdtt*dtdg*dtdag
        dHdbgg = dHdftt*dAdb*dtdg**2 + dHdttt*dtdb*dtdg**2 + dHdttu*dudb*dtdg**2 + 2.d0*dHdtt*dtdg*dtdbg 
        dHdggg = dHdttt*dtdg**3        

        !print '(i2,20d19.10)', i,h,dHda,dHdb,dHdg,dHdaa,dHdbb,dHdgg,dHdag,dHdbg,dHdab,dHdaaa,dHdbbb,dHdggg,dHdaab,dHdabb,dHdaag,dHdagg,dHdbbg,dHdbgg,dHdabg

        dFc_daa(i)   = 2.d0*euda + d*eudaa + 2.d0*dHda + d*dHdaa 
        dFc_dbb(i)   = 2.d0*eudb + d*eudbb + 2.d0*dHdb + d*dHdbb
        dFc_dgg(i)   = d*dHdgg
        dFc_dab(i)   = euda + eudb + d*eudab + dHda + dHdb + d*dHdab 
        dFc_dag(i)   = dHdg + d*dHdag
        dFc_dbg(i)   = dHdg + d*dHdbg
        dFc_daaa(i)  = 3.d0*eudaa + d*eudaaa + 3.d0*dHdaa + d*dHdaaa
        dFc_dbbb(i)  = 3.d0*eudbb + d*eudbbb + 3.d0*dHdbb + d*dHdbbb
        dFc_dggg(i)  = d*dHdggg
        dFc_daab(i)  = eudaa + 2.d0*eudab +d*eudaab + dHdaa + 2.d0*dHdab + d*dHdaab
        dFc_daag(i)  = 2.d0*dHdag + d*dHdaag
        dFc_dabb(i)  = eudbb + 2.d0*eudab + d*eudabb + dHdbb + 2.d0*dHdab + d*dHdabb
        dFc_dbbg(i)  = 2.d0*dHdbg + d*dHdbbg
        dFc_dagg(i)  = dHdgg + d*dHdagg
        dFc_dbgg(i)  = dHdgg + d*dHdbgg
        dFc_dabg(i)  = dHdag + dHdbg + d*dHdabg

        !print '(i2,21d19.10)', i,exc(i),d*(eu+h),dF_drho(i,1),dF_drho(i,2),dF_dgradrho(i)*grad_trho(i) &
        !                      ,dFc_daa(i),dFc_dbb(i),dFc_dgg(i),dFc_dab(i),dFc_dag(i),dFc_dbg(i),dFc_daaa(i),dFc_dbbb(i) &
        !                      ,dFc_dggg(i),dFc_daab(i),dFc_daag(i),dFc_dabb(i),dFc_dbbg(i),dFc_dagg(i),dFc_dbgg(i),dFc_dabg(i)
        !print *,dHdag,dHdaag,d,2.d0*dHdag + d*dHdaag
        !print *,2.d0*dHdbg + d*dHdbbg
        !print *,dHdag + dHdbg + d*dHdabg


     end do
  end if

end subroutine cr_gga_paw_library_3D
! ===
