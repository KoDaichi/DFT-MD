!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  PROGRAM: TDLRMAIN
!
!  AUTHOR(S): K. Tagami et al   Aug. 1 2011
!
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)
!
!=======================================================================

! ======================================================================
!  This is a collection of subroutines for calculating exchange correlation kernels.
! ======================================================================

! ======================= history ======================================
!  ver 1.0 :  2010/3/31
!               applicable to the crystal such as Si.  
!  ver 2.0 :  2011/3/31
!               applicable to the isolated molecule such as C6H6.
!
! ======================================================================

!------------------------------------------------------------------
!!
!!!        Calc Ex on the mesh   (GGAPBE96)
!!
!-------------------------------------------------------------------
subroutine calc_kt_ex_ggapbe( nspin, ispin, chgrhr_l, grad_rho, f2or1, &
     &                        ktmp_ex )

  use m_Const_Parameters,  only : PAI, PAI4, DP
  use m_Parallelization,   only : ista_fftph, iend_fftph
  use m_Crystal_Structure, only : univol

  implicit none

  integer,intent(in)        :: nspin,ispin
  real(kind=DP),intent(in)  :: chgrhr_l(ista_fftph:iend_fftph,ispin)
  real(kind=DP),intent(in)  :: grad_rho(ista_fftph:iend_fftph,nspin)
  real(kind=DP),intent(in)  :: f2or1(ista_fftph:iend_fftph)
  real(kind=DP),intent(out) :: ktmp_ex(ista_fftph:iend_fftph,nspin)
!-------------------------
  if ( ispin==1 ) then
     call ktmp_exchange_spinless
  else
     call ktmp_exchange_polarized
  endif

contains
  
  subroutine ktmp_exchange_spinless
    integer        :: is, i
    real(kind=DP)  :: ctmp1, ctmp2

!    real(kind=DP)  :: density_minimum = 1.0d-20
!    real(kind=DP)  :: density_minimum = 1.0d-6
!    real(kind=DP)  :: density_minimum = 1.0d-8
!    real(kind=DP)  :: density_minimum = 1.0d-12
    real(kind=DP)  :: density_minimum = 1.0d-5
!!
    real(kind=DP) :: fac1, fac2, rho, grho, kf, s, exunif
    real(kind=DP) :: Fx, dFx_ds, d2Fx_ds2
    real(kind=DP) :: ds_drho, d2s_drho2
! ------------------------------ start -------------
    fac1 = ( 3.0d0 *PAI**2 )**(1./3.)
    fac2 = -3.0d0 / PAI4
! ---------------------------- main ---------------
    do is = 1, ispin
       do i = ista_fftph, iend_fftph
          if ( chgrhr_l(i,is) < density_minimum ) cycle
          rho = chgrhr_l(i, is)
          grho = grad_rho(i,is)

          kf = fac1 *rho**(1./3.)
          s = grho / 2.0d0 / kf /rho

          exunif = fac2 *kf
          call derivative_Fx_s( s, Fx, dFx_ds, d2Fx_ds2 )
          call derivaitive_s_rho( rho, s, ds_drho, d2s_drho2 )

          ktmp_ex(i,is ) = exunif *Fx
       End do
    End do
  end subroutine ktmp_exchange_spinless

  subroutine ktmp_exchange_polarized
    integer        :: is, i
    real(kind=DP)  :: ctmp1, ctmp2

!    real(kind=DP)  :: density_minimum = 1.0d-20
!    real(kind=DP)  :: density_minimum = 1.0d-6
!    real(kind=DP)  :: density_minimum = 1.0d-8
    real(kind=DP)  :: density_minimum = 1.0d-12
!!
    real(kind=DP) :: fac1, fac2, rho, grho, kf, s, exunif
    real(kind=DP) :: Fx, dFx_ds, d2Fx_ds2
    real(kind=DP) :: ds_drho, d2s_drho2
! ---------------------------------start ----------------    
    fac1 = ( 3.0d0 *PAI**2 )**(1./3.)
    fac2 = -3.0d0 / PAI4
! ---------------------------------- main ---------
    do is = 1, ispin
       do i = ista_fftph, iend_fftph
          if ( chgrhr_l(i,is) < density_minimum ) cycle
          rho = chgrhr_l(i, is)*2.0d0
          grho = grad_rho(i,is)*2.0d0

          kf = fac1 *rho**(1./3.)
          s = grho / 2.0d0 / kf /rho

          exunif = fac2 *kf
          call derivative_Fx_s( s, Fx, dFx_ds, d2Fx_ds2 )
          call derivaitive_s_rho( rho, s, ds_drho, d2s_drho2 )

          ktmp_ex(i,is ) = exunif* Fx
       End do
    End do
  end subroutine ktmp_exchange_polarized

  subroutine derivaitive_s_rho( rho, s, fn1, fn2 )
    Real(kind=DP), intent(in) :: s, rho
    Real(kind=DP), intent(out) ::fn1, fn2

    fn1 = -4.0D0 /3.0D0  * s /rho               ! ds / drho
    fn2 = -7.0D0 /3.0D0  *fn1 /rho          ! d2s / drho2
  end subroutine derivaitive_s_rho

  subroutine derivative_Fx_s( s, fn0, fn1, fn2 )
    Real(kind=DP), intent(in) :: s
    Real(kind=DP), intent(out) ::fn0, fn1, fn2
    
    Real(kind=DP) :: stmp1, c1, c2
    Real(kind=DP) :: kappa, mu
    Parameter( kappa = 0.804D0,  mu = 0.21951D0 )
! --------------------------------------------    
    stmp1 = s*s
    c1 = 1.0D0 + mu *stmp1 / kappa
    c2 = 1.0D0 - 3.0d0 *mu *stmp1 / kappa
    
    fn0 = 1.0D0 + kappa - kappa / c1
    fn1 = 2.0D0 *mu *s /c1**2                         ! dFx / ds
    fn2 = 2.0D0 *mu *c2 /c1**3                        ! d2Fx / ds2
  end subroutine derivative_Fx_s

end subroutine calc_kt_ex_ggapbe

!------------------------------------------------------------------
!!
!!!        Calc Ec on the mesh   (GGAPBE96)
!!
!-------------------------------------------------------------------
subroutine calc_kt_cr_ggapbe( nspin, ispin, chgrhr_l, grad_trho, f2or1, &
     &                          ktmp_cr )
  use m_Const_Parameters,  only : PAI, PAI4, DP
  use m_Parallelization,   only : ista_fftph, iend_fftph
  use m_Crystal_Structure, only : univol

  implicit none

  integer,intent(in)        :: nspin,ispin
  real(kind=DP),intent(in)  :: chgrhr_l(ista_fftph:iend_fftph,ispin)
  real(kind=DP),intent(in)  :: grad_trho(ista_fftph:iend_fftph)
  real(kind=DP),intent(in)  :: f2or1(ista_fftph:iend_fftph)
  real(kind=DP),intent(out) :: ktmp_cr(ista_fftph:iend_fftph,nspin,nspin)

! --------------------------------
  Real(kind=DP) :: beta, gamma
  Parameter( beta = 0.066725D0, gamma = 0.031091D0 )
!-------------------------
  if ( ispin==1 ) then
     call ktmp_correlation_spinless
  else
     call ktmp_correlation_polarized
  endif

contains

  subroutine ktmp_correlation_spinless
!    real(kind=DP)  :: density_minimum = 1.0d-20
!    real(kind=DP)  :: density_minimum = 1.0d-6
!    real(kind=DP)  :: density_minimum = 1.0d-8
!    real(kind=DP)  :: density_minimum = 1.0d-12
    real(kind=DP)  :: density_minimum = 1.0d-5
! ---------------
    Real(kind=DP) :: rho, grho, ks, kf, rs, t
    real(kind=DP) :: fac1, fac2
    real(kind=DP) :: ec, dec_drs, d2ec_drs2
    real(kind=DP) :: H, dH_drs_tot, d2H_drs2_tot
    Real(kind=DP) :: c1, ctmp1
! 
    integer       :: i
! ------------------------------ start ---------
    fac1 = ( 3.0d0 *PAI**2 )**(1./3.)
    fac2 = -3.0d0 / PAI4
! ----------------------------- main ------------
    Do i = ista_fftph, iend_fftph
       rho  = chgrhr_l(i,1)
       grho = grad_trho(i)

       if ( rho < density_minimum ) cycle

       kf = fac1 *rho**(1.00/3.0d0)
       ks = ( 4.0d0*kf /PAI )**(1.0d0/2.0d0)
       t = abs(grho) / 2.0d0 / ks /rho
       rs   = ( 3.0d0 / PAI4 / rho )**(1.0d0/3.0d0)

       call differential_ec_rs_only( rs, ec, dec_drs, d2ec_drs2 )
       call differeintial_HPBE_rs_only( rs, t, ec, dec_drs, d2ec_drs2, &
            &                           H, dH_drs_tot, d2H_drs2_tot )

       c1 = rs /3.0d0
       ctmp1 = -c1*( dec_drs + dH_drs_tot - c1*d2ec_drs2 - c1*d2H_drs2_tot )
       ktmp_cr( i,1,1 ) = ec +H
    End do
  end subroutine ktmp_correlation_spinless
  
  subroutine ktmp_correlation_polarized
    write(*,*) '!!!!!! Correlation Kernel for GGAPBE96 is not supported '
    stop
  end subroutine ktmp_correlation_polarized

  subroutine differeintial_HPBE_rs_only( rs, t, ec, dec_drs, d2ec_drs2, &
       &                                 f0, f1, f2 )
     Real(kind=DP), intent(in) :: rs, t, ec, dec_drs, d2ec_drs2
     Real(kind=DP), intent(out) :: f0, f1, f2

     Real(kind=DP) :: t2, t4, At2, A2t4
     Real(kind=DP) :: dt_drs, d2t_drs2
     Real(kind=DP) :: A, dA_drs, d2A_drs2
!
     Real(kind=DP) :: dH_drs, d2H_drs2, dH_dt, d2H_dt2
     Real(kind=DP) :: d2H_drsdt
     Real(kind=DP) :: kfn, dkfn_drs, d2kfn_drs2
     Real(kind=DP) :: dkfn_dt, d2kfn_dt2
     Real(kind=DP) :: d2kfn_drsdt
!
     Real(kind=DP) :: ctmp1
     Real(kind=DP) :: phi, gphi3
! ----------------------------- experimental---------------------
     Real(kind=DP) :: prefactor1
     prefactor1 = 0.5d0
! ------------------------- start -----------------------------
     phi = 1.0d0;    gphi3 = gamma *phi**3
! -------------------------------------------------------------
     call derivative_t_rs( rs, t, dt_drs, d2t_drs2 )
     call derivative_A_rs( rs, ec, dec_drs, d2ec_drs2, A, dA_drs, d2A_drs2 )
     call derivative_kfunc1_rs_t( rs, t, A, dA_drs, d2A_drs2, &
          &                       kfn, dkfn_drs, d2kfn_drs2, &
          &                       dkfn_dt,  d2kfn_dt2, &
          &                       d2kfn_drsdt ) 
     ctmp1 = 1.0d0 +kfn
! 
     dH_drs = gphi3 * dkfn_drs / ctmp1
     d2H_drs2 = gphi3 *( -dkfn_drs**2 /ctmp1**2 + d2kfn_drs2 /ctmp1 )

     dH_dt = gphi3 *dkfn_dt / ctmp1
     d2H_dt2 = gphi3 *( -dkfn_dt**2 /ctmp1**2 + d2kfn_dt2 /ctmp1 )

     d2H_drsdt = gphi3*( -dkfn_drs *dkfn_dt /ctmp1**2 + d2kfn_drsdt /ctmp1 )
! ---------
     f0  = gphi3 *log( ctmp1 )
     f1 = dH_drs + dt_drs * dH_dt              ! total derivative
     f2 = d2H_drs2 + 2.0d0 *dt_drs *d2H_drsdt &
          & + d2t_drs2 *dH_dt + dt_drs**2 *d2H_dt2 
! ---------------------------- experimental --------
     f2 = f2 + d2t_drs2 *dH_dt *prefactor1            ! naze ? 0.5
! --------------------------------------------------
   end subroutine differeintial_HPBE_rs_only

   subroutine derivative_kfunc1_rs( rs, t, A, dA_drs, d2A_drs2, f0, f1, f2 )
     Real(kind=DP), intent(in) :: rs, t, A, dA_drs, d2A_drs2
     Real(kind=DP), intent(out) :: f0, f1, f2
!
     Real(kind=DP) :: t2, t4, At2, A2t4, A3t6
     Real(kind=DP) :: q1, q2
     Real(kind=DP) :: ctmp1, ctmp2, ctmp3
!----------------------------- start ------------
     t2 = t*t;  t4 = t2**2
     At2 = A *t2;    A2t4 = At2 **2;  A3t6 = A2t4 *At2

     q1 = 1.0D0 + At2;  q2 = q1 + A2t4
! ----
     ctmp1 = beta /gamma *t2
     ctmp2 = t2*( -2.0d0 *At2 - A2t4 )
     ctmp3 = t4*( -2.0d0 + 6.0d0 *A2t4 + 2.0d0*A3t6 )
! ---------------------------------------------
     f0 = ctmp1 * q1/q2
     f1 = ctmp1 *ctmp2 / q2**2 *dA_drs
     f2 = ctmp1 *( ctmp3 /q2**3 *dA_drs**2 + ctmp2 *d2A_drs2 )
   end subroutine derivative_kfunc1_rs

   subroutine derivative_kfunc1_rs_t( rs, t, A, dA_drs, d2A_drs2, &
        &                               kfn, dkfn_drs, d2kfn_drs2, &
        &                               dkfn_dt, d2kfn_dt2, &
        &                               d2kfn_drsdt )
     Real(kind=DP), intent(in) :: rs, t, A, dA_drs, d2A_drs2
     Real(kind=DP), intent(out) :: kfn, dkfn_drs, d2kfn_drs2
     Real(kind=DP), intent(out) :: dkfn_dt, d2kfn_dt2
     Real(kind=DP), intent(out) :: d2kfn_drsdt
!
     Real(kind=DP) :: t2, t4, At2, A2t4, A3t6
     Real(kind=DP) :: q1, q2
     Real(kind=DP) :: ctmp0, ctmp1, ctmp2, ctmp3, ctmp4, ctmp5, ctmp6
! ------------------------- start ----------------------------
     t2 = t*t;  t4 = t2**2
     At2 = A *t2;    A2t4 = At2 **2;  A3t6 = A2t4 *At2
     q1 = 1.0D0 + At2;  q2 = q1 + A2t4
! -------------
     ctmp0 = beta /gamma
     ctmp1 = ctmp0 *t2
     ctmp2 = t2*( -2.0d0 *At2 - A2t4 )
     ctmp3 = t4*( -2.0d0 + 6.0d0 *A2t4 + 2.0d0*A3t6 )
     ctmp4 = t *( 2.0d0 + 4.0d0 *At2 )
     ctmp5 = 2.0d0 + 6.0d0 *At2 - 18.0d0 *A2t4 - 20.0d0*A3t6
     ctmp6 = -12.0d0 *t *( At2 + A2t4 )
! -----------------------
     kfn = ctmp1 * q1/q2
     dkfn_drs = ctmp1 *ctmp2 / q2**2 *dA_drs
     d2kfn_drs2 = ctmp1 *( ctmp3 /q2**3 *dA_drs**2 + ctmp2 /q2**2 *d2A_drs2 )
!
     dkfn_dt = ctmp0 *ctmp4 / q2**2
     d2kfn_dt2 = ctmp0 *ctmp5 / q2**3
!
     d2kfn_drsdt = ctmp1 *ctmp6 / q2**3 *dA_drs
   end subroutine derivative_kfunc1_rs_t

   subroutine derivative_A_rs( rs, ec, dec_drs, d2ec_drs2, f0, f1, f2 )
     Real(kind=DP), intent(in) :: rs, ec, dec_drs, d2ec_drs2
     Real(kind=DP), intent(out) :: f0, f1, f2

     Real(kind=DP) :: c1, c2, c3
     Real(kind=DP) :: ctmp1, ctmp2, ctmp3
     Real(kind=DP) :: phi, gphi3
! ------------------------------ start ------------
     phi = 1.0d0;       gphi3 = gamma * phi
! ---------------------------------------------
     c1 = ec / gphi3
     c2 = exp( -c1 );    c3 = c2 -1.0D0

     ctmp1 = -beta / gamma /c3**2
     ctmp2 = -1.0d0 /gphi3 *c2
     ctmp3 = -2.0d0 *ctmp1 / c3
! --------------------------------------------------
     f0 = beta / gamma /c3            ! value of A
     f1 = ctmp1 *ctmp2 *dec_drs       ! dA/drs
     f2 = ctmp3 *ctmp2**2 *dec_drs**2  &
          & -ctmp1 *ctmp2 /gphi3 *dec_drs**2 &
          & +ctmp1 *ctmp2 *d2ec_drs2 
   end subroutine derivative_A_rs

   subroutine derivative_t_rs( rs, t, f1, f2 )
     Real(kind=DP), intent(in) :: t, rs
     Real(kind=DP), intent(out) :: f1, f2

     f1 = 3.5D0 *t / rs
     f2 = 2.5D0 *f1 /rs
   end subroutine derivative_t_rs

  subroutine differential_ec_rs_only( rs, ec, dec_drs, d2ec_drs2 )
    real(kind=DP), intent(in) :: rs
    real(kind=DP), intent(out):: ec, dec_drs
    real(kind=DP), intent(out):: d2ec_drs2
    
    real(kind=DP) :: eu, deu_1, deu_2

    call derivative_ec_rs_zeta0( rs, eu, deu_1, deu_2 )
    ec = eu;  dec_drs = deu_1;   d2ec_drs2 = deu_2
  end subroutine differential_ec_rs_only

  subroutine derivative_ec_rs_zeta0( rs, ec0, ec1, ec2 )
    real(kind=DP), intent(in) :: rs
    real(kind=DP), intent(out):: ec0, ec1, ec2
! ---------------------------- const -------------------
    real(kind=DP), parameter :: p = 1.0d0
    real(kind=DP), parameter :: A  = 0.0310907d0,     alpha1 = 0.21370d0
    real(kind=DP), parameter :: beta1 = 7.5957d0,     beta2  = 3.5876d0
    real(kind=DP), parameter :: beta3 = 1.6382d0,     beta4  = 0.49294d0
! --
    call derivative_G_rs( rs, ec0, ec1, ec2, &
       &                  A, alpha1, beta1, beta2, beta3, beta4, p )
  end subroutine derivative_ec_rs_zeta0

  subroutine derivative_G_rs( rs, g0, g1, g2, &
       &           A, alpha1, beta1, beta2, beta3, beta4, p )
    implicit none

    real(kind=DP), intent(in) :: rs, A, alpha1, p
    real(kind=DP), intent(in) :: beta1, beta2, beta3, beta4
    real(kind=DP), intent(out):: g0, g1, g2

    real(kind=DP) :: ctmp1, ctmp2, ctmp3, ctmp4, ctmp5, ctmp6, ctmp7
    real(kind=DP) :: c1, c2
    real(kind=DP) :: Q0, Q1, Q0p, Q1p, Q1p2

    ctmp1 = rs **(1.0d0/2.0d0 )  ! rs **(1/2)
    ctmp2 = ctmp1 *rs            ! rs **(3/2)
    ctmp3 = ctmp1 /rs            ! rs **(-1/2)
    ctmp4 = ctmp3 /rs            ! rs **(-3/2)
    ctmp5 = rs **(p+1.0d0)       ! rs **(p+1)
    ctmp6 = ctmp5 / rs           ! rs **(p)
    ctmp7 = ctmp6 / rs           ! rs **(p-1)

    Q0   = -2.0d0 *A *( 1.0d0 + alpha1 *rs )
    Q1   =  2.0d0 *A *( beta1*ctmp1 + beta2*rs + beta3*ctmp2 +beta4*ctmp5 )
    Q1p  =  A *( beta1*ctmp3 + 2.0d0*beta2 + 3.0d0*beta3*ctmp1 &
         &       + 2.0d0*(p+1.0d0)*beta4*ctmp6 )
    Q0p  = -2.0d0 *A *alpha1

    Q1p2 =  A *( -beta1*ctmp4 + 3.0d0*beta3*ctmp3 &
         &        + 4.0d0 *p*(p+1.0d0)*beta4*ctmp7 ) / 2.0d0

    c1 = log( 1.0d0 + 1.0d0 /Q1 );   c2 = Q1**2 + Q1
    g0 = Q0 *c1
    g1 = -2.0d0 *A *alpha1 *c1 - Q0 *Q1p / c2

    g2 = 2.0d0 *A *alpha1 *Q1p / c2 -( Q0p*Q1p +Q0*Q1p2 )/c2 &
         & + Q0*( 2.0d0*Q1 + 1.0d0 ) *Q1p**2 / c2**2
  end subroutine derivative_G_rs
  
end subroutine calc_kt_cr_ggapbe

