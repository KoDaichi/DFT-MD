module m_SpinOrbit_Potential
! $Id: m_SpinOrbit_Potential.F90 606 2020-04-15 06:45:49Z ktagami $
  use m_Control_Parameters,  only : noncol, ndim_spinor, ndim_chgpot, proj_attribute, &
       &                            SpinOrbit_mode
  use m_Const_Parameters,    only : DP, CMPLDP, zi, ON, Hartree, Neglected
  use m_PseudoPotential,     only : nlmt, ltp, mtp, taup, prodphi, nlmta, ilmt, &
       &                            nloc, ntau
  use m_Ionic_System,        only : natm, itab_spinorbit_addition, ityp, ntyp, &
       &                            speciesname

! == KT_add === 13.2S
  use m_Crystal_Structure,  only : sw_fix_global_quantz_axis,  &
       &                           Global_Quantz_Axis_Fixed, &
       &                           sw_spinorbit_second_variation
  use m_Files,  only : nfout
! ============= 13.2S

  use m_Parallelization,   only : mype
  use mpi

  implicit none
!  include 'mpif.h'

! -----------------------------------
  complex(kind=CMPLDP), target :: MatU_ylm_RC_L0( 2*0+1, -0:0 )
  complex(kind=CMPLDP), target :: MatU_ylm_RC_L1( 2*1+1, -1:1 )
  complex(kind=CMPLDP), target :: MatU_ylm_RC_L2( 2*2+1, -2:2 )
  complex(kind=CMPLDP), target :: MatU_ylm_RC_L3( 2*3+1, -3:3 )
  complex(kind=CMPLDP), target :: MatU_ylm_RC_L4( 2*4+1, -4:4 )
  complex(kind=CMPLDP), target :: MatU_ylm_RC_L5( 2*5+1, -5:5 )
  complex(kind=CMPLDP), target :: MatU_ylm_RC_L6( 2*6+1, -6:6 )
!
  complex(kind=CMPLDP) :: Mat_LS_with_cmplx_ylm_L0(  0:0, 0:0, 4 )
  complex(kind=CMPLDP) :: Mat_LS_with_cmplx_ylm_L1( -1:1,-1:1, 4 )
  complex(kind=CMPLDP) :: Mat_LS_with_cmplx_ylm_L2( -2:2,-2:2, 4 )
  complex(kind=CMPLDP) :: Mat_LS_with_cmplx_ylm_L3( -3:3,-3:3, 4 )
!
  complex(kind=CMPLDP) :: Mat_LS_with_real_ylm_L0( 1, 1, 4 )
  complex(kind=CMPLDP) :: Mat_LS_with_real_ylm_L1( 3, 3, 4 )
  complex(kind=CMPLDP) :: Mat_LS_with_real_ylm_L2( 5, 5, 4 )
  complex(kind=CMPLDP) :: Mat_LS_with_real_ylm_L3( 7, 7, 4 )
! -----------------------------------
  complex(kind=CMPLDP) :: EigenWfns_MatLS_L0( 1*2,1*2 )
  complex(kind=CMPLDP) :: EigenWfns_MatLS_L1( 3*2,3*2 )
  complex(kind=CMPLDP) :: EigenWfns_MatLS_L2( 5*2,5*2 )
  complex(kind=CMPLDP) :: EigenWfns_MatLS_L3( 7*2,7*2 )
! -----------------------------------
  real(kind=DP) :: EigenVals_MatLS_L0( 1*2 )
  real(kind=DP) :: EigenVals_MatLS_L1( 3*2 )
  real(kind=DP) :: EigenVals_MatLS_L2( 5*2 )
  real(kind=DP) :: EigenVals_MatLS_L3( 7*2 )
! -----
!
  Complex(kind=CMPLDP), allocatable  :: dsoc(:,:,:,:)

  real(kind=DP), allocatable :: Mat_SOC_Strength(:,:,:,:)

contains

! ============================================================
!!
!!!  Allocate and deallocate arrays
!!
!===============================================================-
  subroutine m_SO_alloc_Dsoc
    if (nlmt <= 0) call phase_error_with_msg(nfout,' nlmt <=  0 (in m_ES_alloc_Dsoc)',__LINE__,__FILE__)
    allocate( dsoc(nlmt,nlmt,natm,ndim_chgpot )); dsoc = 0.d0
  end subroutine m_SO_alloc_Dsoc

  subroutine m_SO_dealloc_Dsoc
    if (allocated(dsoc)) deallocate(dsoc)
  end subroutine m_SO_dealloc_Dsoc
  
  subroutine m_SO_alloc_Mat_SOC_Strenth
    allocate( Mat_SOC_Strength(nloc,ntau,ntau,natm) )
    Mat_SOC_Strength = 0.0d0
  end subroutine m_SO_alloc_Mat_SOC_Strenth

  subroutine m_SO_dealloc_Mat_SOC_Strenth
    if ( allocated( Mat_SOC_Strength) ) deallocate( Mat_SOC_Strength )
  end subroutine m_SO_dealloc_Mat_SOC_Strenth

! ========================================================================
!!
!!!   calc < Ylm Chi_s | L.S | Chi_s' Ylm' >, where Ylm is real spherical harmonics
!!
!=========================================================================
  subroutine m_SO_calc_MatLS_orb_s_to_f
    real(kind=DP) :: theta, phi
    
    theta = 0.0; phi = 0.0d0

    if ( noncol ) then
    else 
       if ( sw_spinorbit_second_variation == ON &
            &     .and. sw_fix_global_quantz_axis == ON ) then
          theta = acos( Global_Quantz_Axis_Fixed(3))
          phi = atan2( Global_Quantz_Axis_Fixed(2), Global_Quantz_Axis_Fixed(1) )
!          write(nfout,*) "theta = ", theta, phi
       endif
    endif

    call m_SO_calc_MatLS_with_Cmplx_ylm( 0, theta, phi, &
         &                               Mat_LS_with_cmplx_ylm_L0 )
    call m_SO_calc_MatLS_with_Real_ylm( 0, MatU_ylm_RC_L0, &
         &                              Mat_LS_with_cmplx_ylm_L0, &
         &                              Mat_LS_with_real_ylm_L0 )
!
    call m_SO_calc_MatLS_with_Cmplx_ylm( 1, theta, phi, &
         &                               Mat_LS_with_cmplx_ylm_L1 )
    call m_SO_calc_MatLS_with_Real_ylm( 1, MatU_ylm_RC_L1, &
         &                              Mat_LS_with_cmplx_ylm_L1, &
         &                              Mat_LS_with_real_ylm_L1 )
!
    call m_SO_calc_MatLS_with_Cmplx_ylm( 2, theta, phi, &
         &                               Mat_LS_with_cmplx_ylm_L2 )
    call m_SO_calc_MatLS_with_Real_ylm( 2, MatU_ylm_RC_L2, &
         &                              Mat_LS_with_cmplx_ylm_L2, &
         &                              Mat_LS_with_real_ylm_L2 )
!
    call m_SO_calc_MatLS_with_Cmplx_ylm( 3, theta, phi, &
         &                            Mat_LS_with_cmplx_ylm_L3 )
    call m_SO_calc_MatLS_with_Real_ylm( 3, MatU_ylm_RC_L3, &
         &                              Mat_LS_with_cmplx_ylm_L3, &
         &                              Mat_LS_with_real_ylm_L3 )
  end subroutine m_SO_calc_MatLS_orb_s_to_f

  subroutine m_SO_diagonalize_MatLS
    integer :: my_l, lwork, info
    integer :: nsize, nsize0
    real(kind=DP), allocatable :: rwork(:), MatW(:)
    complex(kind=CMPLDP), allocatable :: work(:), MatA(:,:)
    integer :: m1, m2
    
    Do my_l = 1, 3                 ! p-orbital to f-orbital
       nsize0 = 2 *my_l + 1
       nsize = nsize0 *ndim_spinor

       allocate( MatA( nsize, nsize ) ); allocate( MatW( nsize ) )
       MatA = 0.0d0

       call set_Mat_A( my_l )

       lwork = max( 1, 2*nsize -1 )
       allocate( rwork( max( 1, 3*nsize -2 ) )); rwork = 0.0d0
       allocate( work(  max( 1, lwork ) )); work = 0.0d0
!
       call zheev( 'V', 'U', nsize, MatA, nsize, MatW, &
            &       work, lwork, rwork, info )
!
       if ( my_l == 1 ) EigenVals_MatLS_L1 = MatW
       if ( my_l == 2 ) EigenVals_MatLS_L2 = MatW
       if ( my_l == 3 ) EigenVals_MatLS_L3 = MatW

!       Do m1=1, nsize
!          write(*,*) 'Eigen Val m1 ', m1, MatW(m1)
!          Do m2=1, nsize
!             write(*,*) m2, m1, MatA( m2,m1 )
!          End do
!       End do

       call get_eigenwfns_with_real_ylm( my_l )
       deallocate( rwork, work, MatA, MatW )

    End Do
! --
    EigenWfns_MatLS_L0(1,1) = 1.0d0;      EigenWfns_MatLS_L0(2,2) = 1.0d0
!    stop "AAA"
    !    
  contains

    subroutine get_eigenwfns_with_real_ylm( my_l )
      integer, intent(in) :: my_l
      integer :: m1, m2, ma
      integer :: is1
      complex(kind=CMPLDP) :: z1

      Do m2=1, nsize
         Do is1=1, ndim_spinor
            Do m1=1, 2*my_l +1
               z1 = 0.0d0
               if ( my_l == 1 ) then
                  Do ma=-my_l, my_l
                     z1 = z1 +  MatU_ylm_RC_L1( m1, ma )  &
                          &    *MatA( nsize0*(is1-1)+ma+my_l +1, m2 )
                  End do
                  EigenWfns_MatLS_L1( nsize0*(is1-1) +m1, m2 ) = z1
               else if ( my_l == 2 ) then
                  Do ma=-my_l, my_l
                     z1 = z1 + MatU_ylm_RC_L2( m1, ma ) &
                          &    *MatA( nsize0*(is1-1)+ma+my_l +1, m2 )
                  End do
                  EigenWfns_MatLS_L2( nsize0*(is1-1) +m1, m2 ) = z1
               else if ( my_l == 3 ) then
                  Do ma=-my_l, my_l
                     z1 = z1 + MatU_ylm_RC_L3( m1,ma ) &
                          &    *MatA( nsize0*(is1-1)+ma+my_l +1, m2 )
                  End do
                  EigenWfns_MatLS_L3( nsize0*(is1-1) +m1, m2 ) = z1
               endif
            End do
         End do
      End Do
    end subroutine get_eigenwfns_with_real_ylm

    subroutine set_mat_A( my_l )
      integer, intent(in) :: my_l

      integer :: is1, is2, istmp
      integer :: m1, m2

      Do is1=1, ndim_spinor
         Do is2=1, ndim_spinor
            istmp = ( is1 -1 )*ndim_spinor + is2
            Do m1=1, nsize0
               Do m2=1, nsize0
                  
                  if ( my_l == 1 ) then
                     MatA( nsize0*(is1-1)+m1, nsize0*(is2-1)+m2 ) &
                          &    = Mat_LS_with_cmplx_ylm_L1( m1-my_l-1,m2-my_l-1,istmp )
                  else if ( my_l == 2 ) then
                     MatA( nsize0*(is1-1)+m1, nsize0*(is2-1)+m2 ) &
                          &    = Mat_LS_with_cmplx_ylm_L2( m1-my_l-1,m2-my_l-1,istmp )
                  else if ( my_l == 3 ) then
                     MatA( nsize0*(is1-1)+m1, nsize0*(is2-1)+m2 ) &
                          &    = Mat_LS_with_cmplx_ylm_L3( m1-my_l-1,m2-my_l-1,istmp )
                  endif
               End Do
            End Do
         End Do
      End do
    end subroutine set_mat_A

  end subroutine m_SO_diagonalize_MatLS

  subroutine m_SO_print_SOC_on_atoms
    integer :: l_min, l_max, t1 ,t2, it, ia, l1, num, lun
    real(kind=DP), allocatable :: guzai(:)

    if ( SpinOrbit_mode == Neglected ) return

    l_min = 2;   l_max = 4
    t1 = 1; t2 = 1

    allocate( guzai(l_min:l_max) );  guzai = 0.0d0

    lun = 1000
    if ( mype == 0 ) then
       open( unit=lun, file="SOI_phase", status="unknown", form = "formatted" )
       write(lun,'(A,8X,A,9X,A,9X,A)') "# no.  name", "p-orbital", &
            &           "d-orbital", "f-orbital  (unit eV)"
       Do ia=1, natm
          it = ityp(ia)
          Do l1=l_min, nloc
             guzai(l1) = Mat_Soc_Strength(l1,t1,t2,ia)
          End Do
          write(lun,'(I4,2X,A5,3F18.8)') ia, speciesname(it), guzai(l_min:l_max) *Hartree
       End Do
       close(lun)
! --
    end if
    deallocate( guzai )
  end subroutine m_SO_print_SOC_on_atoms

! ===================================================================
!!
!!!  set array "dsoc", zeta *< Ylm Chi_s | L.S | Chi_s' Ylm' > 
!!
!===================================================================
  subroutine m_SO_set_Dsoc_potential1     ! Dsoc  is hermitian.

    if ( noncol ) then
       call goto_noncollinear
    else
       call phase_error_with_msg(nfout,'Not supported: Spin-Orbit potential in the collinear case',__LINE__,__FILE__)
    end if

  contains 
    
    subroutine goto_noncollinear
      integer :: ia, ih, it, il, is
      integer :: iyy, ilmt1, ilmt2
      integer :: l1, l2, m1, m2, t1, t2

      real(kind=DP) :: LScoupling
      Complex(kind=CMPLDP) :: z1

      dsoc = 0.d0

      Do ia=1, natm
         ih = itab_spinorbit_addition(ia)
         if (ih == 0) cycle

         LScoupling = proj_attribute(ih)%LScoupling0 &
              &      *proj_attribute(ih)%LScoupling_scaling_factor

         it = proj_attribute(ih)%ityp
         il = proj_attribute(ih)%l
!
         Do ilmt1=1, ilmt(it)
            Do ilmt2=1, ilmt(it)
               l1 = ltp(ilmt1,it);  m1 = mtp(ilmt1,it);  t1 = taup(ilmt1,it)
               l2 = ltp(ilmt2,it);  m2 = mtp(ilmt2,it);  t2 = taup(ilmt2,it)

!!               if ( t1 /= 1 .or. t2 /= 1 ) cycle
!               if ( t1 /= t2 ) cycle

               if ( l1 == il+1 .and. l2 == l1 ) then

                  Do is=1, ndim_chgpot
                     if ( il == 0 ) then
                        z1 = Mat_LS_with_real_ylm_L0( m1, m2, is )
                     else if ( il == 1 ) then
                        z1 = Mat_LS_with_real_ylm_L1( m1, m2, is )
                     else if ( il == 2 ) then
                        z1 = Mat_LS_with_real_ylm_L2( m1, m2, is )
                     else if ( il == 3 ) then
                        z1 = Mat_LS_with_real_ylm_L3( m1, m2, is )
                     endif
                     dsoc( ilmt1,ilmt2, ia, is ) &
                          &   = LScoupling *prodphi( ih, t1, t2 ) *z1

!                     write(854,*) ilmt1, ilmt2, ia, is, dsoc( ilmt1, ilmt2, ia, is)
!                     write(855,*) ilmt1, ilmt2, ia, is, z1
                  End do
               endif
            End do
         End do
      End do

!      stop 'IIIIIIIIIIIII '

    end subroutine goto_noncollinear

  end subroutine m_SO_set_Dsoc_potential1

  subroutine m_SO_set_Dsoc_potential2
    integer :: ia, it
    integer :: ilmt1, ilmt2, l1, l2, t1, t2, m1, m2
    integer :: is, il
    complex(kind=CMPLDP) :: z1
    
    dsoc = 0.0d0
    
    Do ia=1, natm
       it = ityp(ia)
       
       Do ilmt1=1, ilmt(it)
          Do ilmt2=1, ilmt(it)
             l1 = ltp(ilmt1,it);  m1 = mtp(ilmt1,it);  t1 = taup(ilmt1,it)
             l2 = ltp(ilmt2,it);  m2 = mtp(ilmt2,it);  t2 = taup(ilmt2,it)
             
             if ( l1 /=l2  ) cycle
             il = l1 - 1
             
             Do is=1, ndim_chgpot
                if ( il == 0 ) then
                   z1 = Mat_LS_with_real_ylm_L0( m1, m2, is )
                else if ( il == 1 ) then
                   z1 = Mat_LS_with_real_ylm_L1( m1, m2, is )
                else if ( il == 2 ) then
                   z1 = Mat_LS_with_real_ylm_L2( m1, m2, is )
                else if ( il == 3 ) then
                   z1 = Mat_LS_with_real_ylm_L3( m1, m2, is )
                endif
                
                dsoc( ilmt1, ilmt2, ia, is ) &
                     &     = Mat_Soc_Strength( l1, t1, t2, ia ) *z1
             End do
          End do
       End Do
    End Do

  end subroutine m_SO_set_Dsoc_potential2

! ====================
!
!  calc < Ylm a | LS | Yl'm' b >
!
! =======================
  subroutine m_SO_calc_MatLS_with_Real_ylm( l_in, MatU_ylm_RC, &
       &                                    Mat_LS_with_cmplx_ylm, &
       &                                    Mat_LS_with_real_ylm )
    integer, intent(in) :: l_in
    Complex(kind=CMPLDP), intent(in) :: MatU_ylm_RC( 2*l_in+1, -l_in:l_in )
    complex(kind=CMPLDP), intent(in) :: &
         &    Mat_LS_with_cmplx_ylm( -l_in:l_in, -l_in:l_in, ndim_chgpot )
    complex(kind=CMPLDP), intent(out) :: &
         &    Mat_LS_with_real_ylm( 2*l_in+1, 2*l_in+1, ndim_chgpot )

    integer :: m_min, m_max
    !
    integer :: is
    integer :: m1, m2, ma, mb
!
    Complex(kind=CMPLDP) :: ztmp( ndim_chgpot )
    !
!
!    write(860,*) 'lin = ', l_in
!    write(860,*) Mat_LS_with_cmplx_ylm

    m_min = -l_in;  m_max = l_in
!
!    write(870,*) 'lin = ', l_in
!    Do m1=1, 2*l_in + 1
!       Do ma=m_min, m_max
!          write(870,*) m1, ma, MatU_ylm_RC( m1, ma )
!       End do
!    End do
!
    Do m1=1, 2 *l_in +1
       Do m2=1, 2 *l_in +1
          
          ztmp = 0.0d0
          
          Do ma = m_min, m_max
             Do mb = m_min, m_max
                Do is=1, ndim_chgpot
#if 0
                   ztmp(is) = ztmp(is) + conjg(MatU_ylm_RC( m1,ma ))&
                        &                *Mat_LS_with_cmplx_ylm( ma,mb,is ) &
                        &                *( MatU_ylm_RC( m2,mb ) )
#else
                   ztmp(is) = ztmp(is) + MatU_ylm_RC( m1,ma )&
                        &                *Mat_LS_with_cmplx_ylm( ma,mb,is ) &
                        &                *conjg( MatU_ylm_RC( m2,mb ) )
#endif

                End do
             End do
          End do
          Mat_LS_with_real_ylm( m1,m2,: ) = ztmp(:)
       End do
    End do
  end subroutine m_SO_calc_MatLS_with_Real_ylm
  
  subroutine m_SO_calc_MatLS_with_Cmplx_ylm( l_in, theta, phi, Mat_LS_term )
    integer, intent(in) :: l_in
    real(kind=DP), intent(in) :: theta, phi
    complex(kind=CMPLDP), intent(out) :: &
         &          Mat_LS_term( -l_in:l_in, -l_in:l_in, ndim_chgpot )
    
    real(kind=DP) :: cos_th, sin_th, cos2_th_h, sin2_th_h
    real(kind=DP) :: ctmp_m, ctmp_p, ctmp_0
    
    complex(kind=CMPLDP) :: zi_phi_p, zi_phi_m
    complex(kind=CMPLDP) :: ztmp1, ztmp2, ztmp3
    
    integer :: m1, m2
    
    cos_th = cos( theta );  sin_th = sin( theta )
    cos2_th_h = cos( theta/2.0d0 )**2
    sin2_th_h = 1.0d0 - cos2_th_h
    
    zi_phi_p = exp( zi*phi );    zi_phi_m = exp( -zi*phi )

    Do m1=-l_in, l_in
       Do m2=-l_in, l_in

          if ( m2 == m1 +1 ) then
             ctmp_m = sqrt( (l_in +m1 +1.0d0) *( l_in -m1 ) )
          else
             ctmp_m = 0.0d0
          endif

          if ( m2 == m1 -1 ) then
             ctmp_p = sqrt( (l_in -m1 +1.0d0) *( l_in +m1 ) )
          else
             ctmp_p = 0.0d0
          endif
          !
          if ( m1 == m2 ) then
             ctmp_0 = dble(m1)
          else
             ctmp_0 = 0.0d0
          endif

          !
          ztmp1 =  ctmp_0 *cos_th + ctmp_p *sin_th *zi_phi_m / 2.0d0 &
               &                  + ctmp_m *sin_th *zi_phi_p / 2.0d0 

          ztmp2 = -ctmp_0 *sin_th - ctmp_p *sin2_th_h *zi_phi_m &
               &                  + ctmp_m *cos2_th_h *zi_phi_p

          ztmp3 = -ctmp_0 *sin_th + ctmp_p *cos2_th_h *zi_phi_m &
               &                  - ctmp_m *sin2_th_h *zi_phi_p
          
          Mat_LS_term( m1,m2,1 ) =  ztmp1 /2.0d0         ! (11) component
          Mat_LS_term( m1,m2,4 ) = -ztmp1 /2.0d0         ! (22) 
          Mat_LS_term( m1,m2,2 ) =  ztmp2 /2.0d0         ! (12)
          Mat_LS_term( m1,m2,3 ) =  ztmp3 /2.0d0         ! (21)
       End do
    End Do
  end subroutine m_SO_calc_MatLS_with_Cmplx_ylm

! -------------------------------------
!
!            matrix U ( index1, index2 )  
!                           index1 : px, py, pz ....
!!                          index2 : -1, 0, 1 ....
!
!             px =  1/sqrt(2) [ |1,-1> - |1, 1> ]
!             py = zi/sqrt(2) [ |1,-1> + |1, 1> ]
!             pz =  1,0>
! ---------------------------------------
  subroutine m_SO_set_MatU_ylm_RC
    Real(kind=DP) :: ctmp
    Complex(kind=CMPLDP) :: ztmp
    !
    MatU_ylm_RC_L0 = 0.0d0
    MatU_ylm_RC_L1 = 0.0d0
    MatU_ylm_RC_L2 = 0.0d0
    MatU_ylm_RC_L3 = 0.0d0
    MatU_ylm_RC_L4 = 0.0d0
    ! -------
    ctmp = 1.0d0 / sqrt(2.0d0)

! ------------------------- ktDEBUG ------------ 2012/11/12
    ztmp = -ctmp * zi
!!    ztmp = ctmp * zi
! ------------------------- ktDEBUG ------------ 2012/11/12
    !
! s-orbital
    MatU_ylm_RC_L0(  1,0 ) = 1.0d0
! p-orbital
    MatU_ylm_RC_L1(  1,-1 ) =  ctmp;    MatU_ylm_RC_L1( 1,1 ) = -ctmp
    MatU_ylm_RC_L1(  2,-1 ) =  ztmp;    MatU_ylm_RC_L1( 2,1 ) =  ztmp
    MatU_ylm_RC_L1(  3, 0 ) = 1.0d0
! d-orbital
    MatU_ylm_RC_L2(  1, 0 ) = 1.0d0
    MatU_ylm_RC_L2(  2,-2 ) =  ctmp;    MatU_ylm_RC_L2( 2,2 ) =  ctmp
    MatU_ylm_RC_L2(  3,-2 ) =  ztmp;    MatU_ylm_RC_L2( 3,2 ) = -ztmp
    MatU_ylm_RC_L2(  4,-1 ) =  ztmp;    MatU_ylm_RC_L2( 4,1 ) =  ztmp
    MatU_ylm_RC_L2(  5,-1 ) =  ctmp;    MatU_ylm_RC_L2( 5,1 ) = -ctmp
! f-orbital
    MatU_ylm_RC_L3(  1, 0 ) = 1.0d0
    MatU_ylm_RC_L3(  2,-1 ) =  ctmp;   MatU_ylm_RC_L3( 2,1 ) = -ctmp
    MatU_ylm_RC_L3(  3,-1 ) =  ztmp;   MatU_ylm_RC_L3( 3,1 ) =  ztmp
    MatU_ylm_RC_L3(  4,-2 ) =  ctmp;   MatU_ylm_RC_L3( 4,2 ) =  ctmp
    MatU_ylm_RC_L3(  5,-2 ) =  ztmp;   MatU_ylm_RC_L3( 5,2 ) = -ztmp
    MatU_ylm_RC_L3(  6,-3 ) =  ctmp;   MatU_ylm_RC_L3( 6,3 ) = -ctmp
    MatU_ylm_RC_L3(  7,-3 ) =  ztmp;   MatU_ylm_RC_L3( 7,3 ) =  ztmp
! g-orbital
    MatU_ylm_RC_L4(  1, 0 ) = 1.0d0
    MatU_ylm_RC_L4(  2,-1 ) =  ctmp;   MatU_ylm_RC_L4( 2,1 ) = -ctmp
    MatU_ylm_RC_L4(  3,-1 ) =  ztmp;   MatU_ylm_RC_L4( 3,1 ) =  ztmp
    MatU_ylm_RC_L4(  4,-2 ) =  ctmp;   MatU_ylm_RC_L4( 4,2 ) =  ctmp
    MatU_ylm_RC_L4(  5,-2 ) =  ztmp;   MatU_ylm_RC_L4( 5,2 ) = -ztmp
    MatU_ylm_RC_L4(  6,-3 ) =  ctmp;   MatU_ylm_RC_L4( 6,3 ) = -ctmp
    MatU_ylm_RC_L4(  7,-3 ) =  ztmp;   MatU_ylm_RC_L4( 7,3 ) =  ztmp
    MatU_ylm_RC_L4(  8,-4 ) =  ctmp;   MatU_ylm_RC_L4( 8,4 ) =  ctmp
    MatU_ylm_RC_L4(  9,-4 ) =  ztmp;   MatU_ylm_RC_L4( 9,4 ) = -ztmp
! L=5
    MatU_ylm_RC_L5(  1, 0 ) = 1.0d0
    MatU_ylm_RC_L5(  2,-1 ) =  ctmp;   MatU_ylm_RC_L5(  2,1 ) = -ctmp
    MatU_ylm_RC_L5(  3,-1 ) =  ztmp;   MatU_ylm_RC_L5(  3,1 ) =  ztmp
    MatU_ylm_RC_L5(  4,-2 ) =  ctmp;   MatU_ylm_RC_L5(  4,2 ) =  ctmp
    MatU_ylm_RC_L5(  5,-2 ) =  ztmp;   MatU_ylm_RC_L5(  5,2 ) = -ztmp
    MatU_ylm_RC_L5(  6,-3 ) =  ctmp;   MatU_ylm_RC_L5(  6,3 ) = -ctmp
    MatU_ylm_RC_L5(  7,-3 ) =  ztmp;   MatU_ylm_RC_L5(  7,3 ) =  ztmp
    MatU_ylm_RC_L5(  8,-4 ) =  ctmp;   MatU_ylm_RC_L5(  8,4 ) =  ctmp
    MatU_ylm_RC_L5(  9,-4 ) =  ztmp;   MatU_ylm_RC_L5(  9,4 ) = -ztmp
    MatU_ylm_RC_L5( 10,-5 ) =  ctmp;   MatU_ylm_RC_L5( 10,5 ) = -ctmp
    MatU_ylm_RC_L5( 11,-5 ) =  ztmp;   MatU_ylm_RC_L5( 11,5 ) =  ztmp
! L=6
    MatU_ylm_RC_L6(  1, 0 ) = 1.0d0
    MatU_ylm_RC_L6(  2,-1 ) =  ctmp;   MatU_ylm_RC_L6(  2,1 ) = -ctmp
    MatU_ylm_RC_L6(  3,-1 ) =  ztmp;   MatU_ylm_RC_L6(  3,1 ) =  ztmp
    MatU_ylm_RC_L6(  4,-2 ) =  ctmp;   MatU_ylm_RC_L6(  4,2 ) =  ctmp
    MatU_ylm_RC_L6(  5,-2 ) =  ztmp;   MatU_ylm_RC_L6(  5,2 ) = -ztmp
    MatU_ylm_RC_L6(  6,-3 ) =  ctmp;   MatU_ylm_RC_L6(  6,3 ) = -ctmp
    MatU_ylm_RC_L6(  7,-3 ) =  ztmp;   MatU_ylm_RC_L6(  7,3 ) =  ztmp
    MatU_ylm_RC_L6(  8,-4 ) =  ctmp;   MatU_ylm_RC_L6(  8,4 ) =  ctmp
    MatU_ylm_RC_L6(  9,-4 ) =  ztmp;   MatU_ylm_RC_L6(  9,4 ) = -ztmp
    MatU_ylm_RC_L6( 10,-5 ) =  ctmp;   MatU_ylm_RC_L6( 10,5 ) = -ctmp
    MatU_ylm_RC_L6( 11,-5 ) =  ztmp;   MatU_ylm_RC_L6( 11,5 ) =  ztmp
    MatU_ylm_RC_L6( 12,-6 ) =  ctmp;   MatU_ylm_RC_L6( 12,6 ) =  ctmp
    MatU_ylm_RC_L6( 13,-6 ) =  ztmp;   MatU_ylm_RC_L6( 13,6 ) = -ztmp
  end subroutine m_SO_set_MatU_ylm_RC

end module m_SpinOrbit_Potential

