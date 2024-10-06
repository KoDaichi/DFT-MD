! =========== Contributions ===================================
!
! Through the courtesy of contributors, the following functions are added.
!
! Company:  ASMS Co.,Ltd.
!
! Functions:  BeckeRoussel ex functional
!             modified BJ  ex functional
!
! =============================================================

! ********************************************
!
!  Becke Roussel meta-gga functinal
!                       [ A. D. Becke et al, Phys. Rev. A 39 (1989) 3761 ]
!
! ********************************************
subroutine ex_mgga_br89( nspin, ispin, ista_r, iend_r, iskip_r, &
     &                   rho, norm_grad_rho, grad2_rho, ekin_dens, vx )
  use m_Const_Parameters,  only : DP,PAI
  implicit none

  integer, intent(in) :: ispin, nspin, ista_r, iend_r, iskip_r
  real(kind=DP), intent(in) :: rho( ista_r:iend_r, ispin )
  real(kind=DP), intent(in) :: norm_grad_rho( ista_r:iend_r, nspin )
  real(kind=DP), intent(in) :: grad2_rho( ista_r:iend_r, nspin )
  real(kind=DP), intent(in) :: ekin_dens( ista_r:iend_r, nspin )
  real(kind=DP), intent(out) :: vx( ista_r:iend_r, nspin )

  integer :: is, i
  real(kind=DP) :: c1, c2, c3, c4, d1, d2, factor, coeff1
  real(kind=DP) :: val_d, val_q, val_x, val_y, val_b
  real(kind=DP) :: ctmp

  real(kind=DP), parameter :: delta1 = 1.0E-14
!  real(kind=DP), parameter :: delta1 = 1.0E-18
!  real(kind=DP), parameter :: delta1 = 1.0E-20

  real(kind=DP), parameter :: delta2 = 1.0E-20

  factor = ispin /2.0d0

  coeff1 = 2.0d0 /3.0d0 *PAI**(2.0d0/3.0d0)

  vx = 0.0d0

  Do is=1, ispin
     Do i=ista_r, iend_r, iskip_r
        c1 = ekin_dens(i,is) *factor
        c2 = norm_grad_rho(i,is) *factor
        c3 = rho(i,is) *factor
        c4 = grad2_rho(i,is) *factor

        if ( c3 > delta1 ) then
           val_d = 2.0d0*c1 -c2**2 /c3 /4.0d0
        else
           val_d = 2.0d0*c1
        endif

        val_q = ( c4 -2.0d0 *val_d ) /6.0d0
        if ( c3 > delta1 .and. abs(val_q) > delta1 ) then
           val_y = coeff1 *c3**(5.0d0/3.0d0) /val_q
        else
           val_y = 0.0d0
        endif

        call determine_x_from_y_analytic( val_y, val_x )

#if 0
        if ( c3 > delta1 ) then
           d1 = val_x**3 *exp( -val_x ) /( 8.0d0*PAI*c3 )

           if ( d1 > delta2 ) then
              val_b = d1 **(1.0d0/3.0d0)
           else
              val_b = 0.0d0
           endif
           if ( val_b > delta2 ) then
              d2 = exp( -val_x )
              vx(i,is) = -( 1.0d0 -d2 -d2*val_x/2.0d0  )/val_b
           endif
        endif
#else
        if ( abs(val_x) < 1.0D-8 ) then
           ctmp = 0.50d0
        else
           ctmp = exp( val_x/3.0d0 ) /val_x &
                &      *( 1.0d0 -exp(-val_x) -0.50d0*val_x*exp(-val_x) )
        endif
        vx(i,is) = -2.0*PAI**(1./3.)*c3**(1./3.) *ctmp
#endif
     End Do
  End Do

contains

! ********************************************
!
!  Analytical solution of nonlinear eq.
!           [ E. Proynov et al, Chem. Phy. Lett. 455 (2008) 103. ]
!
! ********************************************
  subroutine determine_x_from_y_analytic( val_y, val_x )
    real(kind=DP), intent(in) :: val_y
    real(kind=DP), intent(out) :: val_x

    real(kind=DP), parameter :: a1 =  1.5255251812009530D+00
    real(kind=DP), parameter :: a2 =  0.4576575543602858D+00
    real(kind=DP), parameter :: a3 =  0.4292036732051034D+00

    real(kind=DP), parameter :: b0 =  0.4771976183772063D+00
    real(kind=DP), parameter :: b1 = -1.7799813494556270D+00
    real(kind=DP), parameter :: b2 =  3.8433841862302150D+00
    real(kind=DP), parameter :: b3 = -9.5912050880518490D+00
    real(kind=DP), parameter :: b4 =  2.1730180285916720D+00
    real(kind=DP), parameter :: b5 = -30.425133851603660D+00
    real(kind=DP), parameter :: c0 =  0.7566445420735584D+00
    real(kind=DP), parameter :: c1 = -2.6363977871370960D+00
    real(kind=DP), parameter :: c2 =  5.4745159964232880D+00
    real(kind=DP), parameter :: c3 = -12.657308127108290D+00
    real(kind=DP), parameter :: c4 =  4.1250584725121360D+00
    real(kind=DP), parameter :: c5 = -30.425133957163840D+00

    real(kind=DP), parameter :: d0 = 0.00004435009886795587D+00
    real(kind=DP), parameter :: d1 = 0.58128653604457910D+00
    real(kind=DP), parameter :: d2 = 66.742764515940610D+00
    real(kind=DP), parameter :: d3 = 434.26780897229770D+00
    real(kind=DP), parameter :: d4 = 824.7765766052239000D+00
    real(kind=DP), parameter :: d5 = 1657.9652731582120D+00
    real(kind=DP), parameter :: e0 = 0.00003347285060926091D+00
    real(kind=DP), parameter :: e1 = 0.47917931023971350D+00
    real(kind=DP), parameter :: e2 = 62.392268338574240D+00
    real(kind=DP), parameter :: e3 = 463.14816427938120D+00
    real(kind=DP), parameter :: e4 = 785.2360350104029000D+00
    real(kind=DP), parameter :: e5 = 1657.962968223273000000D+00

    real(kind=DP), parameter :: BB = 2.085749716493756D+00

    real(kind=DP) :: y1, y2, y3, y4, y5
    real(kind=DP) :: ctmp1, gy, p1y, p2y

    y1 = val_y
    y2 = val_y *y1
    y3 = val_y *y2
    y4 = val_y *y3
    y5 = val_y *y4

    if ( val_y <= 0.0d0 ) then
       ctmp1 = a1 *y1 +a2
       gy = -atan( ctmp1 ) +a3
       p1y = c0 +c1*y1 +c2*y2 +c3*y3 +c4*y4 +c5*y5
       p2y = b0 +b1*y1 +b2*y2 +b3*y3 +b4*y4 +b5*y5
    else
       ctmp1 = BB *y1
       gy = arccsch( ctmp1 ) +2.0d0
       p1y = d0 +d1*y1 +d2*y2 +d3*y3 +d4*y4 +d5*y5
       p2y = e0 +e1*y1 +e2*y2 +e3*y3 +e4*y4 +e5*y5
    endif

    val_x = gy *p1y /p2y
    if ( val_x < 0.0 ) val_x = 0.0d0

  end subroutine determine_x_from_y_analytic

  real*8 function arccsch(x)
    real(kind=DP), intent(in) :: x

    real(kind=DP) :: c1, c2

    c1 = 1.0d0 +sqrt( 1.0d0 +x**2 )
    c2 = c1 /abs(x)
    arccsch = log(c2)

    return
  end function arccsch

end subroutine ex_mgga_br89

! ********************************************
!
!  modified Becke Johnson meta-gga functinal
!                       [ F. Tran et al, Phys. Rev. Lett. 102 (2009) 226401 ]
!
! ********************************************
subroutine set_cval_tb09( val_g, val_c )
  use m_Const_Parameters,     only : DP
  implicit none

  real(kind=DP), intent(in) :: val_g
  real(kind=DP), intent(out) :: val_c
!
  real(kind=DP), parameter :: alpha = -0.012D+00,  beta = 1.023D+00
!
  val_c = alpha + beta *sqrt( val_g )

end subroutine set_cval_tb09

subroutine ex_mgga_tb09( nspin, ispin, ista_r, iend_r, iskip_r, &
     &                    rho, norm_grad_rho, grad2_rho, ekin_dens, vx, val_c )
  use m_Const_Parameters,     only : DP,PAI
  use m_Crystal_Structure,    only : univol
  use m_Parallelization,      only : MPI_CommGroup, npes, mype
  use mpi

  implicit none
!  include 'mpif.h'

  integer, intent(in) :: ispin, nspin, ista_r, iend_r, iskip_r
  real(kind=DP), intent(in) :: rho( ista_r:iend_r, ispin )
  real(kind=DP), intent(in) :: norm_grad_rho( ista_r:iend_r, nspin)
  real(kind=DP), intent(in) :: grad2_rho( ista_r:iend_r, nspin )
  real(kind=DP), intent(in) :: ekin_dens( ista_r:iend_r, nspin )
  real(kind=DP), intent(in) :: val_c
  real(kind=DP), intent(out) :: vx( ista_r:iend_r, nspin )

!  real(kind=DP), parameter :: delta1 = 1.0E-14
!  real(kind=DP), parameter :: delta1 = 1.0E-18

  real(kind=DP), parameter :: delta1 = 1.0E-20
  real(kind=DP), parameter :: delta2 = 1.0E-20
!  real(kind=DP), parameter :: delta1 = 1.0E-10
!  real(kind=DP), parameter :: delta2 = 1.0E-10

  integer :: i, is
  real(kind=DP) :: ctmp1, fac

  vx = 0.0d0
  call ex_mgga_br89( nspin, ispin, ista_r, iend_r, iskip_r, &
     &               rho, norm_grad_rho, grad2_rho, ekin_dens, vx )

  fac = sqrt(5.0d0/12.0d0) /PAI

  Do is=1, nspin
     Do i=ista_r, iend_r, iskip_r
        if ( rho(i,is) > delta1 ) then
           ctmp1 = 2.0d0 *ekin_dens(i,is) /rho(i,is)
           if ( ctmp1 > delta2 ) then
              ctmp1 = fac *sqrt(ctmp1)
           else
              ctmp1 = 0.0d0
           endif
        else
           ctmp1 = 0.0d0
        endif

        vx(i,is) = val_c *vx(i,is) + ( 3.0d0*val_c -2.0d0 ) *ctmp1
     End do
  End Do

end subroutine ex_mgga_tb09

subroutine set_gval_tb09( nspin, ista_r, iend_r, rho, norm_grad_rho, weight, gval )
  use m_Const_Parameters,     only : DP,PAI


  use m_FFT,                  only : fft_box_size_CD_3D
  use m_Parallelization,      only : MPI_CommGroup, npes, nrank_g, mpi_ke_world
  use mpi

  implicit none
!  include 'mpif.h'

  integer, intent(in) :: nspin, ista_r, iend_r
  real(kind=DP), intent(in) :: rho( ista_r:iend_r, nspin )
  real(kind=DP), intent(in) :: norm_grad_rho( ista_r:iend_r, nspin)
  real(kind=DP), intent(in) :: weight( ista_r:iend_r )
  real(kind=DP), intent(out) :: gval

  real(kind=DP), parameter :: alpha = -0.012D+00,  beta = 1.023D+00
  real(kind=DP), parameter :: delta1 = 1.0E-20

  integer :: i, is, ierr
  real(kind=DP) :: csum_mpi, csum

  csum = 0.0d0
  Do is=1, nspin
     Do i=ista_r, iend_r
        if ( abs( rho(i,is) ) > delta1 ) then
           csum = csum + norm_grad_rho(i,is) /rho(i,is) *weight(i)
        end if
     End Do
  End Do

  if ( nspin == 2 ) csum = csum /2.0d0


#ifdef CD_FFT_ALL
    if ( npes > 1 ) then
     call mpi_allreduce( csum, csum_mpi, 1, mpi_double_precision, mpi_sum, &
          &              MPI_CommGroup, ierr )
       csum = csum_mpi
    end if
#else
    if ( nrank_g > 1 ) then
       call mpi_allreduce( csum, csum_mpi, 1, mpi_double_precision, mpi_sum, &
            &              mpi_ke_world, ierr )
       csum = csum_mpi
    end if
#endif
  csum = csum /product(fft_box_size_CD_3D(1:3,1))
  gval = csum

end subroutine set_gval_tb09
