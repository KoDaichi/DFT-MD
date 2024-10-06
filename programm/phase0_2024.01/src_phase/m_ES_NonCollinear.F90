module m_ES_NonCollinear
! $Id: m_ES_NonCollinear.f90 587 2018-09-27 02:25:59Z ktagami $
  use m_Const_Parameters,   only : DP, CMPLDP, zi, yes, PAI, &
       &                            BuiltIn, ByProjector, ByPawPot, ZeffApprox, &
       &                           BUCS, CARTS, CRDTYP, DELTA10, ReadFromPP

  use m_Control_Parameters,  only : noncol, ndim_spinor, ndim_magmom, ndim_chgpot, &
       &                            kimg, nspin, ipri, iprixc, SpinOrbit_Mode

  use m_Ionic_System,       only : natm, ntyp, ityp, mag_direction0_atomtyp, &
       &                           natm2, pos, iwei, magmom_local_now

  use m_PseudoPotential,     only : dion_so, dion_scr_noncl, q_noncl, flg_soc, &
       &                            pot_has_soc, dion0_noncl, &
       &                            nlmt, ipaw, dion, dion_paw, q, ltp, mtp, jtp, &
       &                            nac, nac_p, nlmta_phi, &
       &                            ilmt

  use m_Parallelization,      only : ista_kngp, iend_kngp, ista_fftp, iend_fftp, &
       &                             ista_fftph, iend_fftph, &
       &                             myrank_ggacmp, nrank_ggacmp, npes, ierr, &
       &                             MPI_CommGroup

  use m_Files,                 only : nfout

!!!

  use m_SpinOrbit_Potential,    only : MatU_ylm_RC_L0, MatU_ylm_RC_L1, &
       &                              MatU_ylm_RC_L2, MatU_ylm_RC_L3, &
       &                              dsoc

!  use m_Crystal_Structure,      only : mag_direction0_global

!  use m_Crystal_Structure,      only : nopr, op, ig01, tau, m_CS_op_in_PUCD

! === KT_add === 2014/08/04
  use m_Crystal_Structure,   only : threshold_vorticity, threshold_helicity
! ============== 2014/08/04
  use mpi

  implicit none

!  include 'mpif.h'
  integer istatus(mpi_status_size)

  complex(kind=CMPLDP), allocatable :: Spinor_EigenWfn0_atomtyp(:,:,:)
  complex(kind=CMPLDP), allocatable :: Spinor_EigenWfn0_global(:,:)
  complex(kind=CMPLDP), allocatable :: factor_f( :,:,:,:,: )

  real(kind=DP) :: Global_Quantz_Axis_now(3)

contains

! ------------------------------------------------------------------------
!!
!!!   Spinor eigenwfns calculated from the initial mag. orientaion.
!!
! -----------------------------------------------------------------------
  subroutine m_ES_alloc_spinor_eigenwfn_0
    allocate( Spinor_EigenWfn0_atomtyp( ntyp, ndim_spinor, ndim_spinor ) )
    Spinor_EigenWfn0_atomtyp = 0.0d0
  end subroutine m_ES_alloc_spinor_eigenwfn_0

  subroutine m_ES_dealloc_spinor_eigenwfn_0
    deallocate( Spinor_EigenWfn0_atomtyp )
  end subroutine m_ES_dealloc_spinor_eigenwfn_0

  subroutine m_ES_set_spinor_eigenwfn_0
    integer :: it
    real(kind=DP) :: c_mx, c_my, c_mz
    real(kind=DP) :: theta, phi
    real(kind=DP) :: c1, s1
    complex(kind=CMPLDP) :: z1

    Do it=1, ntyp
       c_mx = mag_direction0_atomtyp( it,1 )
       c_my = mag_direction0_atomtyp( it,2 )
       c_mz = mag_direction0_atomtyp( it,3 )
!
       theta = atan2( sqrt( c_mx**2 +c_my**2), c_mz )
       phi = atan2( c_my, c_mx )
!
!       write(*,*) 'cmx = ', c_mx, c_my, c_mz
!       write(*,*) 'Th phi = ', theta, phi
!
       c1 = cos( theta/2.0d0 );  s1 = sin( theta/2.0d0 )
       z1 = exp( zi *phi /2.0d0 )
!!! ------------------------------------------------------ Kai_+ ---
!!!       Spinor_EigenWfn0_atomtyp( it,1,1 ) = conjg(z1)*c1
!!!       Spinor_EigenWfn0_atomtyp( it,2,1 ) = z1 *s1
!!! ------------------------------------------------------ Kai_- ---
!!!       Spinor_EigenWfn0_atomtyp( it,1,2 ) = -conjg(z1)*s1
!!!       Spinor_EigenWfn0_atomtyp( it,2,2 ) = z1 *c1
! ------------------------------------------------------ Kai_+ ---
       Spinor_EigenWfn0_atomtyp( it,1,1 ) =  z1 *c1
       Spinor_EigenWfn0_atomtyp( it,2,1 ) = -z1 *s1
! ------------------------------------------------------ Kai_- ---
       Spinor_EigenWfn0_atomtyp( it,1,2 ) = conjg(z1) *s1
       Spinor_EigenWfn0_atomtyp( it,2,2 ) = conjg(z1) *c1
    End do

  end subroutine m_ES_set_spinor_eigenwfn_0

! ------------------------------------------------------------------------
!!
!!!        Factor F defined defined by A. Coroso et al
!!!                                    see. PRB 71 (2005) 115106.
!!
! -------------------------------------------------------------------------
  subroutine m_ES_alloc_factor_fss
    allocate( factor_f( nlmt, nlmt, ndim_spinor, ndim_spinor, ntyp ) )
    factor_f = 0.0d0
  end subroutine m_ES_alloc_factor_fss

  subroutine m_ES_dealloc_factor_fss
    deallocate( factor_f )
  end subroutine m_ES_dealloc_factor_fss

  subroutine m_ES_set_factor_fss
    integer :: it
    integer :: ilmt1, ilmt2, l1, l2, m1, m2
    integer :: is1, is2

    integer :: ni, number_of_mj

    real(kind=DP) :: jval1, jval2, mj
    real(kind=DP) :: c1, c2
    complex(kind=CMPLDP) :: z1, z2, zsum

    Do it=1, ntyp

       if ( .not. pot_has_soc(it)  ) cycle

       Do ilmt1=1, nlmt
          l1 = ltp( ilmt1, it )
          m1 = mtp( ilmt1, it )
          jval1 = jtp( ilmt1, it ) -0.50d0

          Do ilmt2=1, nlmt
             l2 = ltp( ilmt2, it )
             m2 = mtp( ilmt2, it )
             jval2 = jtp( ilmt2, it ) -0.5d0

             if ( l1 /= l2 ) cycle
             if ( abs(jval1-jval2) < 1.0D-4 )  cycle

             number_of_mj = nint( 2*jval1 +1 )

             write(*,*) 'mj = ', number_of_mj


             Do is1=1, ndim_spinor
                Do is2=1, ndim_spinor

                   zsum = 0.0d0
                   Do ni=1, number_of_mj
                      mj = -jval1 + ni -1
                      c1 = coeff_ClebshGordan( l1, jval1, mj, is1 )
                      c2 = coeff_ClebshGordan( l1, jval1, mj, is2 )
                      z1 = coeff_MatUs( l1, jval1, mj, m1, is1 )
                      z2 = coeff_MatUs( l1, jval1, mj, m2, is2 )
                      zsum = zsum + c1 *z1 *c2 *conjg( z2 )
                   End do
!
                   factor_f( ilmt1, ilmt2, is1, is2, it ) = zsum
                End do
             End do
          End do
       End do
    End Do

  end subroutine m_ES_set_factor_fss

!------------------------------------------------------------------
!!
!!!     Generalizatin of Dion and q
!!
! ----------------------------------------------------------------
  subroutine m_ES_set_Mat_q_noncl
    integer :: it

    q_noncl = 0.0d0

    if ( SpinOrbit_mode == BuiltIn ) then
       Do it=1, ntyp
          call set_q_so( it )
       End do
    else
       Do it=1, ntyp
          q_noncl( 1:nlmt, 1:nlmt, 1, it ) = q( 1:nlmt, 1:nlmt, it )
          q_noncl( 1:nlmt, 1:nlmt, ndim_chgpot, it ) = q( 1:nlmt, 1:nlmt, it )
       End do
    endif

  contains

    subroutine set_q_so( it )
      integer, intent(in) :: it

      integer :: lmt1, lmt2, lmt3, lmt4
      integer :: is1, is2, is3, is4
      integer :: is_tmp
      complex(kind=CMPLDP) :: z1, zsum

      Do lmt1 = 1, nlmt
         Do lmt2 = 1, nlmt

            Do is_tmp = 1, ndim_chgpot
               is1 = ( is_tmp -1 ) /ndim_spinor + 1
               is2 = mod( is_tmp -1,2 ) +1

               zsum = 0.0d0

               Do lmt3 = 1, nlmt
                  Do lmt4 = 1, nlmt
                     Do is3=1, ndim_spinor
                        is4 = is3
                        z1 = q( lmt1, lmt2, it ) &
                             &   *factor_f( lmt1, lmt3, is1, is3, it ) &
                             &   *factor_f( lmt4, lmt2, is4, is2, it )
                        zsum = zsum + z1
                     End do
                  End do
               End do
                  !
               q_noncl( lmt1, lmt2, is_tmp, it ) = zsum
               !
            End do
         End do
      End do
    end subroutine set_q_so

  end subroutine m_ES_set_Mat_q_noncl

  subroutine m_ES_init_Mat_dion0_noncl
    integer :: it, ia

    if ( SpinOrbit_Mode == BuiltIn ) then      ! full-relativisitic
       Do ia=1, natm
          it = ityp(ia)
          call set_dion_so( ia, it )

       End do

    else if ( SpinOrbit_Mode == ByProjector .or. &
         &    SpinOrbit_Mode == ZeffApprox ) then      ! soc by parameter
       Do ia=1, natm
          it = ityp(ia)
          dion0_noncl( 1:nlmt, 1:nlmt, 1, ia ) = dion( 1:nlmt, 1:nlmt, it )
          dion0_noncl( 1:nlmt, 1:nlmt, ndim_chgpot, ia ) &
               &                 = dion( 1:nlmt, 1:nlmt, it )

          call add_Dsoc_term_noncl( ia )
       End do

    else if ( SpinOrbit_Mode == ReadFromPP .or. &
         &    SpinOrbit_Mode == ByPawPot ) then

       Do ia=1, natm
          it = ityp(ia)
          dion0_noncl( 1:nlmt, 1:nlmt, 1, ia ) = dion( 1:nlmt, 1:nlmt, it )
          dion0_noncl( 1:nlmt, 1:nlmt, ndim_chgpot, ia ) &
               &                 = dion( 1:nlmt, 1:nlmt, it )

          call add_Dsoc_term_noncl( ia )
       End do

    else                                 ! neglect soc
       Do ia=1, natm
          it = ityp(ia)
          dion0_noncl( 1:nlmt, 1:nlmt, 1, ia ) = dion( 1:nlmt, 1:nlmt, it )
          dion0_noncl( 1:nlmt, 1:nlmt, ndim_chgpot, ia ) &
               &                 = dion( 1:nlmt, 1:nlmt, it )
       End do
    end if

    dion_scr_noncl = dion0_noncl

  contains

    subroutine add_Dsoc_term_noncl( ia )
      integer, intent (in ) :: ia
      integer is

      Do is=1, ndim_chgpot
         dion0_noncl(:,:,is,ia) = dion0_noncl(:,:,is,ia) + dsoc(:,:,ia,is)
      End do
    end subroutine add_Dsoc_term_noncl

! --
    subroutine set_dion_so( ia, it )
      integer, intent (in ) :: ia, it
!
      integer :: is1, is2, is_tmp
      integer :: ilmt1, ilmt2, l1, m1, l2, m2
      real(kind=DP) :: jval1, jval2
      complex(kind=CMPLDP) :: z1, z2
      !
      Do ilmt1=1, nlmt
         l1 = ltp( ilmt1, it )
         m1 = mtp( ilmt1, it )
         jval1 = jtp( ilmt1, it )

         Do ilmt2=1, nlmt
            l2 = ltp( ilmt2, it )
            m2 = mtp( ilmt2, it )
            jval2 = jtp( ilmt2, it )

            if ( l1 /= l2 ) cycle
            if ( abs(jval1-jval2) < 1.0D-4 )  cycle

            z1 = dion( ilmt1, ilmt2, it )
            Do is1=1, ndim_spinor
               Do is2=1, ndim_spinor
                  is_tmp = ( is1 -1 )*ndim_spinor + is2

                  z2 = z1 *factor_f( ilmt1, ilmt2, is1, is2, it )
                  dion0_noncl( ilmt1, ilmt2, is_tmp, ia ) = z2
               End do
            End do
         End do
      End do

    end subroutine set_dion_so

  end subroutine m_ES_init_Mat_dion0_noncl

  subroutine m_ES_update_Mat_dion0_noncl         ! for PAW
    integer :: it, ia
    complex(kind=CMPLDP), allocatable :: dion_paw_ssrep(:,:,:,:)

    allocate( dion_paw_ssrep(nlmt,nlmt,ndim_magmom,natm) )
    call m_ES_MagMom_To_DensMat_DionPaw( dion_paw, dion_paw_ssrep )

    if ( SpinOrbit_Mode == BuiltIn ) then      ! full-relativisitic
       call phase_error_with_msg(nfout, 'kt : Not suppported spin-orbit + paw',__LINE__,__FILE__)

       Do ia=1, natm
          it = ityp(ia)
          call set_dion_so( ia, it )      !  ?????????????
       End do

    else if ( SpinOrbit_Mode == ByProjector .or. &
         &    SpinOrbit_Mode == ZeffApprox ) then      ! soc by parameter
       Do ia=1, natm
          dion0_noncl( 1:nlmt, 1:nlmt, :, ia ) = dion_paw_ssrep( 1:nlmt, 1:nlmt, :, ia )
          call add_Dsoc_term_noncl( ia )
       End do

    else if ( SpinOrbit_Mode == ReadFromPP .or. &
         &    SpinOrbit_Mode == ByPawPot ) then
       Do ia=1, natm
          dion0_noncl( 1:nlmt, 1:nlmt, :, ia ) = dion_paw_ssrep( 1:nlmt, 1:nlmt, :, ia )
          call add_Dsoc_term_noncl( ia )
       End do

    else                                 ! neglect soc
       Do ia=1, natm
          dion0_noncl( 1:nlmt, 1:nlmt, :, ia ) = dion_paw_ssrep( 1:nlmt, 1:nlmt, :, ia )
       End do
    end if

    dion_scr_noncl = dion0_noncl
    deallocate( dion_paw_ssrep )

  contains

    subroutine add_Dsoc_term_noncl( ia )
      integer, intent (in ) :: ia
      integer is

      Do is=1, ndim_chgpot
         dion0_noncl(:,:,is,ia) = dion0_noncl(:,:,is,ia) + dsoc(:,:,ia,is)
      End do
    end subroutine add_Dsoc_term_noncl

! --
    subroutine set_dion_so( ia, it )
      integer, intent (in ) :: ia, it
!
      integer :: is1, is2, is_tmp
      integer :: ilmt1, ilmt2, l1, m1, l2, m2
      real(kind=DP) :: jval1, jval2
      complex(kind=CMPLDP) :: z1, z2
      !
      Do ilmt1=1, nlmt
         l1 = ltp( ilmt1, it )
         m1 = mtp( ilmt1, it )
         jval1 = jtp( ilmt1, it )

         Do ilmt2=1, nlmt
            l2 = ltp( ilmt2, it )
            m2 = mtp( ilmt2, it )
            jval2 = jtp( ilmt2, it )

            if ( l1 /= l2 ) cycle
            if ( abs(jval1-jval2) < 1.0D-4 )  cycle

            z1 = dion( ilmt1, ilmt2, it )
            Do is1=1, ndim_spinor
               Do is2=1, ndim_spinor
                  is_tmp = ( is1 -1 )*ndim_spinor + is2

                  z2 = z1 *factor_f( ilmt1, ilmt2, is1, is2, it )
                  dion0_noncl( ilmt1, ilmt2, is_tmp, ia ) = z2
               End do
            End do
         End do
      End do

    end subroutine set_dion_so

  end subroutine m_ES_update_Mat_dion0_noncl

  subroutine m_ES_set_Mat_dion_scr_noncl( VlhxcQ_l, dhub_aimag )
    real(kind=DP), intent(in) :: VlhxcQ_l( nlmt, nlmt, natm, ndim_magmom )
    real(kind=DP), optional, intent(in) :: dhub_aimag( nlmt, nlmt, natm, ndim_magmom )

    complex(kind=CMPLDP), allocatable :: vlhxcQ_ssrep( :,:,:,: )
    complex(kind=CMPLDP) :: PauliMatrix( ndim_magmom, ndim_spinor, ndim_spinor )
!
    integer :: ia, it

    allocate( vlhxcQ_ssrep( nlmt,nlmt,natm,ndim_magmom ) )
    vlhxcQ_ssrep = 0.0d0
!
    if ( present(dhub_aimag) ) then
       call m_ES_MagMom_to_DensMat_vlhxcQ( vlhxcQ_l, vlhxcQ_ssrep, dhub_aimag )
    else
       call m_ES_MagMom_to_DensMat_vlhxcQ( vlhxcQ_l, vlhxcQ_ssrep )
    endif
!
    call m_ES_set_Pauli_matrix( PauliMatrix )

!------------------------------------------
    dion_scr_noncl = dion0_noncl
! -------------------------------------

    if ( SpinOrbit_Mode == BuiltIn ) then      ! full-relativisitic
       Do ia=1, natm
          it = ityp(ia)
          call add_VlhxcQ_term_noncl_B( ia )
       End do

    else if ( SpinOrbit_Mode == ByProjector .or. &
         &    SpinOrbit_Mode == ZeffApprox ) then      ! soc by parameter
       Do ia=1, natm
          it = ityp(ia)
          call add_VlhxcQ_term_noncl_A( ia )
       End do

    else if ( SpinOrbit_Mode == ReadFromPP .or. &
         &    SpinOrbit_Mode == ByPawPot ) then
       Do ia=1, natm
          it = ityp(ia)
          call add_VlhxcQ_term_noncl_A( ia )
       End do

    else                                 ! neglect soc
       Do ia=1, natm
          it = ityp(ia)
          call add_VlhxcQ_term_noncl_A( ia )
       End do
    end if

    deallocate( vlhxcQ_ssrep )

  contains

    subroutine add_VlhxcQ_term_noncl_A( ia )
      integer, intent (in ) :: ia
      integer is

      Do is=1, ndim_chgpot
         dion_scr_noncl(:,:,is,ia) = dion0_noncl(:,:,is,ia) +vlhxcQ_ssrep( :,:,ia,is)
      End do
    end subroutine add_VlhxcQ_term_noncl_A

    subroutine add_VlhxcQ_term_noncl_B( ia )
      integer, intent (in ) :: ia

      integer :: is1, is2, is3, is4
      integer :: lmt1, lmt2, lmt3, lmt4
      integer :: ialph, is_tmp

      complex(kind=CMPLDP) :: zsum, z1

      Do lmt1 = 1, nlmt
         Do lmt2 = 1, nlmt

            Do is_tmp = 1, ndim_chgpot
               is1 = ( is_tmp -1 ) /ndim_spinor + 1
               is2 = mod( is_tmp -1,2 ) +1

               zsum = 0.0d0

               Do lmt3 = 1, nlmt
                  Do lmt4 = 1, nlmt
                     Do is3=1, ndim_spinor
                        Do is4=1, ndim_spinor

                           Do ialph=1, ndim_magmom
                              z1 = PauliMatrix( ialph, is3, is4 ) &
                                   &     *VlhxcQ_ssrep( lmt3, lmt4, ia, ialph ) &
                                   &     *factor_f( lmt1, lmt3, is1, is3, it ) &
                                   &     *factor_f( lmt4, lmt2, is4, is2, it )
                              zsum = zsum + z1
                           End do
                        End do
                     End do
                  End do
               End do
!
               dion_scr_noncl(:,:,is_tmp,ia) = dion0_noncl(:,:,is_tmp,ia) +zsum
!
            End do
         End do
      End do

    end subroutine add_VlhxcQ_term_noncl_B

  end subroutine m_ES_set_Mat_dion_scr_noncl

!--------------------------------------------------------------
!!
!!!   Clebsh Gordan coeff etc defined by A. Corso et al
!!!                                    see. PRB 71 (2005) 115106.
!!
! -------------------------------------------------------------
  complex(kind=CMPLDP) function coeff_MatUs( l, j, mj, mprime, is)
    integer, intent(in) :: l, is, mprime
    real(kind=DP), intent(in) :: j, mj

    integer :: m
    complex(kind=CMPLDP) :: z1

    call phase_error_with_msg(nfout, "Not checked",__LINE__,__FILE__)

    z1 = 0.0d0
    if ( abs( j -l -0.5d0 ) < 1.0D-4 ) then
       m = nint( mj - 0.5d0 )
       if ( is == 1 .and. m >=-l ) then
          if ( l == 0 ) z1 = MatU_ylm_RC_L0( mprime, m )
          if ( l == 1 ) z1 = MatU_ylm_RC_L1( mprime, m )
          if ( l == 2 ) z1 = MatU_ylm_RC_L2( mprime, m )
          if ( l == 3 ) z1 = MatU_ylm_RC_L3( mprime, m )
       else if ( is == 2 .and. m <= l-1 ) then
          if ( l == 0 ) z1 = MatU_ylm_RC_L0( mprime, m+1 )
          if ( l == 1 ) z1 = MatU_ylm_RC_L1( mprime, m+1 )
          if ( l == 2 ) z1 = MatU_ylm_RC_L2( mprime, m+1 )
          if ( l == 3 ) z1 = MatU_ylm_RC_L3( mprime, m+1 )
       endif
    else if ( abs ( j -l +0.5d0 ) < 1.0D-4 ) then
       m = nint( mj + 0.5d0 )
       if ( is == 1 .and. m >=-l+1 ) then
          if ( l == 0 ) z1 = MatU_ylm_RC_L0( mprime, m-1 )
          if ( l == 1 ) z1 = MatU_ylm_RC_L1( mprime, m-1 )
          if ( l == 2 ) z1 = MatU_ylm_RC_L2( mprime, m-1 )
          if ( l == 3 ) z1 = MatU_ylm_RC_L3( mprime, m-1 )
       else if ( is == 2 .and. m <= l ) then
          if ( l == 0 ) z1 = MatU_ylm_RC_L0( mprime, m )
          if ( l == 1 ) z1 = MatU_ylm_RC_L1( mprime, m )
          if ( l == 2 ) z1 = MatU_ylm_RC_L2( mprime, m )
          if ( l == 3 ) z1 = MatU_ylm_RC_L3( mprime, m )
       endif
    end if

    coeff_MatUs = conjg(z1)

  end function coeff_MatUs

  real(kind=DP) function coeff_ClebshGordan( l, j, mj, is )
    integer, intent(in) :: l, is
    real(kind=DP), intent(in) :: j, mj

    real(kind=DP) :: c1
    integer :: m

    c1 = 0.0d0
    if ( abs( j -l -0.5d0 ) < 1.0D-4 ) then
       m = nint( mj - 0.5d0 )
       if ( is == 1 ) c1 = sqrt( ( l+m+1)/(2.d0*l+1.d0) )
       if ( is == 2 ) c1 = sqrt( ( l-m)/(2.d0*l+1.d0) )
    else if ( abs ( j -l +0.5d0 ) < 1.0D-4 ) then
       m = nint( mj + 0.5d0 )
       if ( is == 1 ) c1 = sqrt( ( l-m+1)/(2.d0*l+1.d0) )
       if ( is == 2 ) c1 =-sqrt( ( l+m)/(2.d0*l+1.d0) )
    endif

    coeff_ClebshGordan = c1
    return

  end function coeff_ClebshGordan

  subroutine m_ES_set_Pauli_matrix( PauliMatrix )
    complex(kind=CMPLDP), intent(out) :: &
         &                PauliMatrix( ndim_magmom, ndim_spinor, ndim_spinor )

    PauliMatrix = 0.0d0

    PauliMatrix( 1,1,1 ) = 1.0d0;     PauliMatrix( 1,2,2 ) = 1.0d0   ! Sigma 0
    PauliMatrix( 2,1,2 ) = 1.0d0;     PauliMatrix( 2,2,1 ) = 1.0d0   ! Sigma X
    PauliMatrix( 3,1,2 ) =   -zi;     PauliMatrix( 3,2,1 ) =    zi   ! Sigma Y
    PauliMatrix( 4,1,1 ) = 1.0d0;     PauliMatrix( 4,2,2 ) =-1.0d0   ! Sigma Z

  end subroutine m_ES_set_Pauli_matrix

! -------------------------------
  subroutine m_ES_set_Mat_hsr_with_soc( hsr_ssrep, hsi_ssrep, &
       &                                hsr_ssrep_with_soc, hsi_ssrep_with_soc )
    real(kind=DP), intent(in) :: hsr_ssrep( natm, nlmt, nlmt, ndim_chgpot )
    real(kind=DP), intent(in) :: hsi_ssrep( natm, nlmt, nlmt, ndim_chgpot )
    real(kind=DP), intent(out):: hsr_ssrep_with_soc( natm, nlmt, nlmt, ndim_chgpot )
    real(kind=DP), intent(out):: hsi_ssrep_with_soc( natm, nlmt, nlmt, ndim_chgpot )

    integer :: ia, it
    integer :: is1, is2, is3, is4, is_tmp, is_tmp2
    integer :: ilmt1, ilmt2, ilmt3, ilmt4
    integer :: l1, l2, l3, l4, m1, m2, m3, m4
    real(kind=DP) :: jval1, jval2, jval3, jval4
    complex(kind=CMPLDP) :: z1, z2, z3, zsum

    Do ia=1, natm
       it = ityp(ia)

       Do ilmt1=1, nlmt
          l1 = ltp( ilmt1, it )
          m1 = mtp( ilmt1, it )
          jval1 = jtp( ilmt1, it )

          Do ilmt2=1, nlmt
             l2 = ltp( ilmt2, it )
             m2 = mtp( ilmt2, it )
             jval2 = jtp( ilmt2, it )

             Do is_tmp = 1, ndim_chgpot
                is1 = ( is_tmp -1 ) /ndim_spinor + 1
                is2 = mod( is_tmp -1, ndim_spinor ) +1

                zsum = 0.0d0
                Do ilmt3=1, nlmt
                   l3 = ltp( ilmt3, it )
                   m3 = mtp( ilmt3, it )
                   jval3 = jtp( ilmt3, it )

                   if ( l1 /= l3 ) cycle
                   if ( abs(jval1-jval3) > 1.0E-4 ) cycle

                   Do ilmt4=1, nlmt
                      l4 = ltp( ilmt4, it )
                      m4 = mtp( ilmt4, it )
                      jval4 = jtp( ilmt4, it )

                      if ( l2 /= l4 ) cycle
                      if ( abs(jval2-jval4) > 1.0E-4 ) cycle

                      Do is3=1, ndim_spinor
                         Do is4=1, ndim_spinor
                            is_tmp2 = ( is3 -1 )*ndim_spinor + is4

                            z1 = factor_f( ilmt3, ilmt1, is3, is1, it )
                            z2 = factor_f( ilmt2, ilmt4, is2, is4, it )

                            z3 = dcmplx( hsr_ssrep( ia, ilmt3, ilmt4, is_tmp2 ), &
                                 &       hsi_ssrep( ia, ilmt3, ilmt4, is_tmp2 ) )
                            zsum = zsum + z1 *z3 *z2
                         End do
                      End do
                   End do
                End do

                hsr_ssrep_with_soc( ia, ilmt1, ilmt2, is_tmp ) = real(zsum)
                hsi_ssrep_with_soc( ia, ilmt1, ilmt2, is_tmp ) = aimag(zsum)
             End do
          End do
       End do
    End Do
 end subroutine m_ES_set_Mat_hsr_with_soc
!

! ----------------------------------------------------------------------
!!
!!!   Converters between the ss-repsentaion and magmom-representaion
!!
! ----------------------------------------------------------------------
  subroutine m_ES_MagMom_to_DensMat_Gspace( Rho_MagMom, DensMat )
    real(kind=DP), intent(in)  :: Rho_MagMom( ista_kngp:iend_kngp,kimg,ndim_magmom )
    real(kind=DP), intent(out) :: DensMat(    ista_kngp:iend_kngp,kimg,ndim_chgpot )
    ! --
    !                                       note : ndim_magmom = 4
! --
    integer i, k
    Real(kind=DP) :: c_rho(kimg), c_mx (kimg), c_my (kimg), c_mz (kimg)
    Real(kind=DP) :: c_n11(kimg), c_n12(kimg), c_n21(kimg), c_n22(kimg)
!
    DensMat = 0.0d0
    Do i=ista_kngp, iend_kngp

       Do k=1, kimg
          c_rho(k) = Rho_MagMom( i,k,1 )
          c_mx(k)  = Rho_MagMom( i,k,2 )
          c_my(k)  = Rho_MagMom( i,k,3 )
          c_mz(k)  = Rho_MagMom( i,k,4 )
       End do
       ! -------------------- real part ------
       c_n11(1) = ( c_rho(1) + c_mz(1) ) /2.0d0
       c_n22(1) = ( c_rho(1) - c_mz(1) ) /2.0d0
       c_n12(1) = ( c_mx(1)  + c_my(2) ) /2.0d0
       c_n21(1) = ( c_mx(1)  - c_my(2) ) /2.0d0
       ! ---------------------- imaginary part --
       c_n11(2) = ( c_rho(2) + c_mz(2) ) /2.0d0
       c_n22(2) = ( c_rho(2) - c_mz(2) ) /2.0d0
       c_n12(2) = ( c_mx(2)  - c_my(1) ) /2.0d0
       c_n21(2) = ( c_mx(2)  + c_my(1) ) /2.0d0
       ! ----
       Do k=1,  kimg
          DensMat( i,k,1 ) = c_n11(k)
          DensMat( i,k,2 ) = c_n12(k)
          DensMat( i,k,3 ) = c_n21(k)
          DensMat( i,k,4 ) = c_n22(k)
       End do
    End Do
  end subroutine m_ES_MagMom_to_DensMat_Gspace

! -----------------------------------------------------------------
!!
!!!  n11(G),n12(G), n21(G), n22(G) ---> n(G), mx(G), my(G), mz(G)
!!
! -----------------------------------------------------------------
  subroutine m_ES_DensMat_To_MagMom_Gspace( DensMat, Rho_MagMom )

    real(kind=DP), intent(in)   :: DensMat(    ista_kngp:iend_kngp,kimg,ndim_chgpot )
    real(kind=DP), intent(out)  :: Rho_MagMom( ista_kngp:iend_kngp,kimg,ndim_magmom )
!
    integer i, k
    Real(kind=DP) :: c_rho(kimg), c_mx (kimg), c_my (kimg), c_mz (kimg)
    Real(kind=DP) :: c_n11(kimg), c_n12(kimg), c_n21(kimg), c_n22(kimg)
    !
    Rho_MagMom = 0.0d0
    Do i=ista_kngp, iend_kngp

       Do k=1,  kimg
          c_n11(k) = DensMat( i,k,1 )
          c_n12(k) = DensMat( i,k,2 )
          c_n21(k) = DensMat( i,k,3 )
          c_n22(k) = DensMat( i,k,4 )
       End do
! -------------------- real part ------
       c_rho(1) = c_n11(1) + c_n22(1)
       c_mx(1)  = c_n12(1) + c_n21(1)
       c_my(1)  = c_n21(2) - c_n12(2)
       c_mz(1)  = c_n11(1) - c_n22(1)
! ---------------------- imaginary part --
       c_rho(2) = c_n11(2) + c_n22(2)
       c_mx(2)  = c_n12(2) + c_n21(2)
       c_my(2)  =-c_n21(1) + c_n12(1)
       c_mz(2)  = c_n11(2) - c_n22(2)
! ----
       Do k=1, kimg
          Rho_MagMom( i,k,1 ) = c_rho(k)
          Rho_MagMom( i,k,2 ) = c_mx(k)
          Rho_MagMom( i,k,3 ) = c_my(k)
          Rho_MagMom( i,k,4 ) = c_mz(k)
       End do
! --
    End Do
  end subroutine m_ES_DensMat_To_MagMom_Gspace

  subroutine m_ES_DensMat_To_MagMom_hsr( natm, nlmt, hsr_densmat, hsi_densmat, &
       &                                 hsr_magmom, hsi_magmom )
    integer, intent(in) :: natm, nlmt
    real(kind=DP), intent(in)   :: hsr_densmat( natm, nlmt, nlmt, ndim_chgpot )
    real(kind=DP), intent(in)   :: hsi_densmat( natm, nlmt, nlmt, ndim_chgpot )
    real(kind=DP), intent(out)  :: hsr_magmom(  natm, nlmt, nlmt, ndim_magmom )
    real(kind=DP), intent(out)  :: hsi_magmom(  natm, nlmt, nlmt, ndim_magmom )
!
    integer i, k
    integer :: ia, lmt1, lmt2
    Real(kind=DP) :: c_rho(kimg), c_mx (kimg), c_my (kimg), c_mz (kimg)
    Real(kind=DP) :: c_n11(kimg), c_n12(kimg), c_n21(kimg), c_n22(kimg)
!
    hsr_magmom = 0.0d0
    Do lmt2=1, nlmt
       Do lmt1=1, nlmt
          Do ia=1, natm
             c_n11(1) = hsr_densmat( ia, lmt1, lmt2, 1 )
             c_n12(1) = hsr_densmat( ia, lmt1, lmt2, 2 )
             c_n21(1) = hsr_densmat( ia, lmt1, lmt2, 3 )
             c_n22(1) = hsr_densmat( ia, lmt1, lmt2, 4 )
             !
             c_n11(2) = hsi_densmat( ia, lmt1, lmt2, 1 )
             c_n12(2) = hsi_densmat( ia, lmt1, lmt2, 2 )
             c_n21(2) = hsi_densmat( ia, lmt1, lmt2, 3 )
             c_n22(2) = hsi_densmat( ia, lmt1, lmt2, 4 )
!
! -------------------- real part ------
             c_rho(1) = c_n11(1) + c_n22(1)
             c_mx(1)  = c_n12(1) + c_n21(1)
             c_my(1)  = c_n21(2) - c_n12(2)
             c_mz(1)  = c_n11(1) - c_n22(1)
! ---------------------- imaginary part --
             c_rho(2) = c_n11(2) + c_n22(2)
             c_mx(2)  = c_n12(2) + c_n21(2)
             c_my(2)  =-c_n21(1) + c_n12(1)
             c_mz(2)  = c_n11(2) - c_n22(2)
!
             hsr_magmom( ia, lmt1, lmt2, 1 ) = c_rho(1)
             hsr_magmom( ia, lmt1, lmt2, 2 ) = c_mx(1)
             hsr_magmom( ia, lmt1, lmt2, 3 ) = c_my(1)
             hsr_magmom( ia, lmt1, lmt2, 4 ) = c_mz(1)
!
             hsi_magmom( ia, lmt1, lmt2, 1 ) = c_rho(2)
             hsi_magmom( ia, lmt1, lmt2, 2 ) = c_mx(2)
             hsi_magmom( ia, lmt1, lmt2, 3 ) = c_my(2)
             hsi_magmom( ia, lmt1, lmt2, 4 ) = c_mz(2)
          End do
       End do
! --
    End Do
  end subroutine m_ES_DensMat_To_MagMom_hsr

  subroutine m_ES_DensMat_To_MagMom_hsr0( hsr_densmat, hsi_densmat, hsr_magmom )
    real(kind=DP), intent(in)   :: hsr_densmat( natm, nlmt, nlmt, ndim_chgpot )
    real(kind=DP), intent(in)   :: hsi_densmat( natm, nlmt, nlmt, ndim_chgpot )
    real(kind=DP), intent(out)  :: hsr_magmom(  natm, nlmt, nlmt, ndim_magmom )
!
    integer i, k
    integer :: ia, lmt1, lmt2
    Real(kind=DP) :: c_rho(kimg), c_mx (kimg), c_my (kimg), c_mz (kimg)
    Real(kind=DP) :: c_n11(kimg), c_n12(kimg), c_n21(kimg), c_n22(kimg)
!
    hsr_magmom = 0.0d0
    Do lmt2=1, nlmt
       Do lmt1=1, nlmt
          Do ia=1, natm
             c_n11(1) = hsr_densmat( ia, lmt1, lmt2, 1 )
             c_n12(1) = hsr_densmat( ia, lmt1, lmt2, 2 )
             c_n21(1) = hsr_densmat( ia, lmt1, lmt2, 3 )
             c_n22(1) = hsr_densmat( ia, lmt1, lmt2, 4 )
             !
             c_n11(2) = hsi_densmat( ia, lmt1, lmt2, 1 )
             c_n12(2) = hsi_densmat( ia, lmt1, lmt2, 2 )
             c_n21(2) = hsi_densmat( ia, lmt1, lmt2, 3 )
             c_n22(2) = hsi_densmat( ia, lmt1, lmt2, 4 )
!
! -------------------- real part ------
             c_rho(1) = c_n11(1) + c_n22(1)
             c_mx(1)  = c_n12(1) + c_n21(1)
             c_my(1)  = c_n21(2) - c_n12(2)
             c_mz(1)  = c_n11(1) - c_n22(1)
! ---------------------- imaginary part --
             c_rho(2) = c_n11(2) + c_n22(2)
             c_mx(2)  = c_n12(2) + c_n21(2)
             c_my(2)  =-c_n21(1) + c_n12(1)
             c_mz(2)  = c_n11(2) - c_n22(2)

             hsr_magmom( ia, lmt1, lmt2, 1 ) = c_rho(1)
             hsr_magmom( ia, lmt1, lmt2, 2 ) = c_mx(1)
             hsr_magmom( ia, lmt1, lmt2, 3 ) = c_my(1)
             hsr_magmom( ia, lmt1, lmt2, 4 ) = c_mz(1)
          End do
       End do
! --
    End Do
  end subroutine m_ES_DensMat_To_MagMom_hsr0

  subroutine m_ES_MagMom_To_DensMat_hsr( natm, nlmt, hsr_magmom, hsi_magmom, &
       &                                 hsr_densmat, hsi_densmat )
    integer, intent(in) :: natm, nlmt
    real(kind=DP), intent(in)  :: hsr_magmom(  natm, nlmt, nlmt, ndim_magmom )
    real(kind=DP), intent(in)  :: hsi_magmom(  natm, nlmt, nlmt, ndim_magmom )
    real(kind=DP), intent(out)   :: hsr_densmat( natm, nlmt, nlmt, ndim_chgpot )
    real(kind=DP), intent(out)   :: hsi_densmat( natm, nlmt, nlmt, ndim_chgpot )
!
    integer i, k
    integer :: ia, lmt1, lmt2
    Real(kind=DP) :: c_rho(2), c_mx(2), c_my(2), c_mz(2)
    Real(kind=DP) :: c_n11(kimg), c_n12(kimg), c_n21(kimg), c_n22(kimg)
!
    hsr_densmat = 0.0d0
    hsi_densmat = 0.0d0
    Do lmt2=1, nlmt
       Do lmt1=1, nlmt
          Do ia=1, natm

             c_rho(1) = hsr_magmom( ia, lmt1, lmt2, 1 )
             c_mx(1)  = hsr_magmom( ia, lmt1, lmt2, 2 )
             c_my(1)  = hsr_magmom( ia, lmt1, lmt2, 3 )
             c_mz(1)  = hsr_magmom( ia, lmt1, lmt2, 4 )

             c_rho(2) = hsi_magmom( ia, lmt1, lmt2, 1 )
             c_mx(2)  = hsi_magmom( ia, lmt1, lmt2, 2 )
             c_my(2)  = hsi_magmom( ia, lmt1, lmt2, 3 )
             c_mz(2)  = hsi_magmom( ia, lmt1, lmt2, 4 )

       ! -------------------- real part ------
             c_n11(1) = ( c_rho(1) + c_mz(1) ) /2.0d0
             c_n22(1) = ( c_rho(1) - c_mz(1) ) /2.0d0
             c_n12(1) = ( c_mx(1)  + c_my(2) ) /2.0d0
             c_n21(1) = ( c_mx(1)  - c_my(2) ) /2.0d0
       ! ---------------------- imaginary part --
             c_n11(2) = ( c_rho(2) + c_mz(2) ) /2.0d0
             c_n22(2) = ( c_rho(2) - c_mz(2) ) /2.0d0
             c_n12(2) = ( c_mx(2)  - c_my(1) ) /2.0d0
             c_n21(2) = ( c_mx(2)  + c_my(1) ) /2.0d0
       ! ----
             hsr_densmat( ia, lmt1, lmt2, 1 ) = c_n11(1)
             hsr_densmat( ia, lmt1, lmt2, 2 ) = c_n12(1)
             hsr_densmat( ia, lmt1, lmt2, 3 ) = c_n21(1)
             hsr_densmat( ia, lmt1, lmt2, 4 ) = c_n22(1)

             hsi_densmat( ia, lmt1, lmt2, 1 ) = c_n11(2)
             hsi_densmat( ia, lmt1, lmt2, 2 ) = c_n12(2)
             hsi_densmat( ia, lmt1, lmt2, 3 ) = c_n21(2)
             hsi_densmat( ia, lmt1, lmt2, 4 ) = c_n22(2)

          End do
       End do
! --
    End Do
  end subroutine m_ES_MagMom_To_DensMat_hsr

  subroutine m_ES_MagMom_To_DensMat_hsr0( hsr_magmom, hsr_densmat, hsi_densmat )
    real(kind=DP), intent(in)  :: hsr_magmom(  natm, nlmt, nlmt, ndim_magmom )
    real(kind=DP), intent(out)   :: hsr_densmat( natm, nlmt, nlmt, ndim_chgpot )
    real(kind=DP), intent(out)   :: hsi_densmat( natm, nlmt, nlmt, ndim_chgpot )
!
    integer i, k
    integer :: ia, lmt1, lmt2
    Real(kind=DP) :: c_rho, c_mx, c_my, c_mz
    Real(kind=DP) :: c_n11(kimg), c_n12(kimg), c_n21(kimg), c_n22(kimg)
!
    hsr_densmat = 0.0d0
    hsi_densmat = 0.0d0
    Do lmt2=1, nlmt
       Do lmt1=1, nlmt
          Do ia=1, natm

             c_rho = hsr_magmom( ia, lmt1, lmt2, 1 )
             c_mx  = hsr_magmom( ia, lmt1, lmt2, 2 )
             c_my  = hsr_magmom( ia, lmt1, lmt2, 3 )
             c_mz  = hsr_magmom( ia, lmt1, lmt2, 4 )
       ! -------------------- real part ------
             c_n11(1) = ( c_rho + c_mz ) /2.0d0
             c_n22(1) = ( c_rho - c_mz ) /2.0d0
             c_n12(1) = c_mx /2.0d0
             c_n21(1) = c_mx /2.0d0
       ! ---------------------- imaginary part --
             c_n11(2) = 0.0d0
             c_n22(2) = 0.0d0
             c_n12(2) = -c_my /2.0d0
             c_n21(2) =  c_my /2.0d0
       ! ----
             hsr_densmat( ia, lmt1, lmt2, 1 ) = c_n11(1)
             hsr_densmat( ia, lmt1, lmt2, 2 ) = c_n12(1)
             hsr_densmat( ia, lmt1, lmt2, 3 ) = c_n21(1)
             hsr_densmat( ia, lmt1, lmt2, 4 ) = c_n22(1)

             hsi_densmat( ia, lmt1, lmt2, 1 ) = c_n11(2)
             hsi_densmat( ia, lmt1, lmt2, 2 ) = c_n12(2)
             hsi_densmat( ia, lmt1, lmt2, 3 ) = c_n21(2)
             hsi_densmat( ia, lmt1, lmt2, 4 ) = c_n22(2)

          End do
       End do
! --
    End Do
  end subroutine m_ES_MagMom_To_DensMat_hsr0

  subroutine m_ES_MagMom_to_DensMat_vlhxcQ0( vlhxcQ_magmom, vlhxcQ_densmat )
    real(kind=DP),        intent(in) :: vlhxcQ_magmom( nlmt, nlmt, natm, ndim_magmom )
    complex(kind=CMPLDP), intent(out):: vlhxcQ_densmat( nlmt, nlmt, natm, ndim_chgpot )

    Real(kind=DP) :: c_rho, c_mx, c_my, c_mz
    Real(kind=DP) :: c_n11(kimg), c_n12(kimg), c_n21(kimg), c_n22(kimg)
!
    integer :: ia, lmt1, lmt2
!
!    c_rho(2) = 0.0d0;  c_mx(2) = 0.0d0;   c_my(2) = 0.0d0;   c_mz(2) = 0.0d0
!
    vlhxcQ_densmat = 0.0d0
    Do ia=1, natm
       Do lmt2=1, nlmt
          Do lmt1=1, nlmt
             c_rho = vlhxcQ_magmom( lmt1, lmt2, ia, 1 )
             c_mx  = vlhxcQ_magmom( lmt1, lmt2, ia, 2 )
             c_my  = vlhxcQ_magmom( lmt1, lmt2, ia, 3 )
             c_mz  = vlhxcQ_magmom( lmt1, lmt2, ia, 4 )
             !
             c_n11(1) = c_rho + c_mz
             c_n22(1) = c_rho - c_mz
             c_n12(1) = c_mx
             c_n21(1) = c_mx
! ---------------------- imaginary part --
             c_n11(2) = 0.0d0
             c_n22(2) = 0.0d0
             c_n12(2) = -c_my
             c_n21(2) =  c_my
! ----
             vlhxcQ_densmat( lmt1, lmt2, ia, 1 ) = dcmplx( c_n11(1), c_n11(2) )
             vlhxcQ_densmat( lmt1, lmt2, ia, 2 ) = dcmplx( c_n12(1), c_n12(2) )
             vlhxcQ_densmat( lmt1, lmt2, ia, 3 ) = dcmplx( c_n21(1), c_n21(2) )
             vlhxcQ_densmat( lmt1, lmt2, ia, 4 ) = dcmplx( c_n22(1), c_n22(2) )
          End do
       End do
    End Do
  end subroutine m_ES_MagMom_to_DensMat_vlhxcQ0

  subroutine m_ES_MagMom_to_DensMat_vlhxcQ( vlhxcQ_r_magmom, vlhxcQ_densmat, &
       &                                     vlhxcQ_i_magmom )
    real(kind=DP),         intent(in) :: vlhxcQ_r_magmom( nlmt, nlmt, natm, ndim_magmom )
    real(kind=DP),optional,intent(in) :: vlhxcQ_i_magmom( nlmt, nlmt, natm, ndim_magmom )
    complex(kind=CMPLDP), intent(out):: vlhxcQ_densmat( nlmt, nlmt, natm, ndim_chgpot )

    Real(kind=DP) :: c_rho(kimg), c_mx(kimg), c_my(kimg), c_mz(kimg)
    Real(kind=DP) :: c_n11(kimg), c_n12(kimg), c_n21(kimg), c_n22(kimg)
!
    integer :: ia, lmt1, lmt2

    vlhxcQ_densmat = 0.0d0
    Do ia=1, natm
       Do lmt2=1, nlmt
          Do lmt1=1, nlmt
             c_rho(1) = vlhxcQ_r_magmom( lmt1, lmt2, ia, 1 )
             c_mx(1)  = vlhxcQ_r_magmom( lmt1, lmt2, ia, 2 )
             c_my(1)  = vlhxcQ_r_magmom( lmt1, lmt2, ia, 3 )
             c_mz(1)  = vlhxcQ_r_magmom( lmt1, lmt2, ia, 4 )
!
             if ( present(vlhxcQ_i_magmom) ) then
                c_rho(2) = vlhxcQ_i_magmom( lmt1, lmt2, ia, 1 )
                c_mx(2)  = vlhxcQ_i_magmom( lmt1, lmt2, ia, 2 )
                c_my(2)  = vlhxcQ_i_magmom( lmt1, lmt2, ia, 3 )
                c_mz(2)  = vlhxcQ_i_magmom( lmt1, lmt2, ia, 4 )
             else
                c_rho(2) = 0.0d0
                c_mx(2) = 0.0d0
                c_my(2) = 0.0d0
                c_mz(2) = 0.0d0
             endif
       ! -------------------- real part ------
             c_n11(1) = c_rho(1) + c_mz(1)
             c_n22(1) = c_rho(1) - c_mz(1)
             c_n12(1) = c_mx(1)  + c_my(2)
             c_n21(1) = c_mx(1)  - c_my(2)
       ! ---------------------- imaginary part --
             c_n11(2) = c_rho(2) + c_mz(2)
             c_n22(2) = c_rho(2) - c_mz(2)
             c_n12(2) = c_mx(2)  - c_my(1)
             c_n21(2) = c_mx(2)  + c_my(1)
       ! ----
             vlhxcQ_densmat( lmt1, lmt2, ia, 1 ) = dcmplx( c_n11(1), c_n11(2) )
             vlhxcQ_densmat( lmt1, lmt2, ia, 2 ) = dcmplx( c_n12(1), c_n12(2) )
             vlhxcQ_densmat( lmt1, lmt2, ia, 3 ) = dcmplx( c_n21(1), c_n21(2) )
             vlhxcQ_densmat( lmt1, lmt2, ia, 4 ) = dcmplx( c_n22(1), c_n22(2) )
          End do
       End do
    End Do

  end subroutine m_ES_MagMom_to_DensMat_vlhxcQ

  subroutine m_ES_MagMom_to_DensMat_vlhxcl( vlhxc_l, vlhxc_ssrep )

    real(kind=DP), intent(in):: vlhxc_l(ista_kngp:iend_kngp,kimg,ndim_magmom)
    real(kind=DP), intent(out) :: vlhxc_ssrep(ista_kngp:iend_kngp,kimg,ndim_chgpot)
!
    integer :: i, k
    Real(kind=DP) :: c_rho(kimg), c_mx(kimg), c_my(kimg), c_mz(kimg)
    Real(kind=DP) :: c_n11(kimg), c_n12(kimg), c_n21(kimg), c_n22(kimg)

    vlhxc_ssrep = 0.0d0
    Do i=ista_kngp, iend_kngp

       Do k=1, kimg
          c_rho(k) = vlhxc_l( i,k,1 )
          c_mx(k)  = vlhxc_l( i,k,2 )
          c_my(k)  = vlhxc_l( i,k,3 )
          c_mz(k)  = vlhxc_l( i,k,4 )
       End do
       ! -------------------- real part ------
       c_n11(1) = c_rho(1) + c_mz(1)
       c_n22(1) = c_rho(1) - c_mz(1)
       c_n12(1) = c_mx(1)  + c_my(2)
       c_n21(1) = c_mx(1)  - c_my(2)
       ! ---------------------- imaginary part --
       c_n11(2) = c_rho(2) + c_mz(2)
       c_n22(2) = c_rho(2) - c_mz(2)
       c_n12(2) = c_mx(2)  - c_my(1)
       c_n21(2) = c_mx(2)  + c_my(1)
       ! ----
       Do k=1, kimg
          vlhxc_ssrep( i,k,1 ) = c_n11(k)
          vlhxc_ssrep( i,k,2 ) = c_n12(k)
          vlhxc_ssrep( i,k,3 ) = c_n21(k)
          vlhxc_ssrep( i,k,4 ) = c_n22(k)
       End do
! --
    End Do

  end subroutine m_ES_MagMom_to_DensMat_vlhxcl

  subroutine m_ES_DensMat_to_MagMom_vlhxcl( vlhxc_ssrep, vlhxc_l )

    real(kind=DP), intent(in) :: vlhxc_ssrep(ista_kngp:iend_kngp,kimg,ndim_chgpot)
    real(kind=DP), intent(out):: vlhxc_l(ista_kngp:iend_kngp,kimg,ndim_magmom)
!
    integer :: i, k
    Real(kind=DP) :: c_rho(kimg), c_mx(kimg), c_my(kimg), c_mz(kimg)
    Real(kind=DP) :: c_n11(kimg), c_n12(kimg), c_n21(kimg), c_n22(kimg)

    vlhxc_l = 0.0d0
    Do i=ista_kngp, iend_kngp

       Do k=1, kimg
          c_n11(k) = vlhxc_ssrep( i,k,1 )
          c_n12(k) = vlhxc_ssrep( i,k,2 )
          c_n21(k) = vlhxc_ssrep( i,k,3 )
          c_n22(k) = vlhxc_ssrep( i,k,4 )
       End do
! -------------------- real part ------
       c_rho(1) = ( c_n11(1) + c_n22(1) ) /2.0d0
       c_mx(1)  = ( c_n12(1) + c_n21(1) ) /2.0d0
       c_my(1)  = ( c_n21(2) - c_n12(2) ) /2.0d0
       c_mz(1)  = ( c_n11(1) - c_n22(1) ) /2.0d0
! ---------------------- imaginary part --
       c_rho(2) = ( c_n11(2) + c_n22(2) ) /2.0d0
       c_mx(2)  = ( c_n12(2) + c_n21(2) ) /2.0d0
       c_my(2)  = (-c_n21(1) + c_n12(1) ) /2.0d0
       c_mz(2)  = ( c_n11(2) - c_n22(2) ) /2.0d0
! ----
       Do k=1, kimg
          vlhxc_l( i,k,1 ) = c_rho(k)
          vlhxc_l( i,k,2 ) = c_mx(k)
          vlhxc_l( i,k,3 ) = c_my(k)
          vlhxc_l( i,k,4 ) = c_mz(k)
       End do
! --
    End Do

  end subroutine m_ES_DensMat_to_MagMom_vlhxcl


  subroutine m_ES_DensMat_To_MagMom_dm( dm_ssrep, dm_magmom, max2lp, max_projs )
    integer, intent(in) :: max2lp, max_projs
    complex(kind=CMPLDP), intent(in)  :: dm_ssrep( max2lp, max2lp, max_projs, &
         &                                         natm, ndim_chgpot )
    real(kind=DP), intent(out)  :: dm_magmom( max2lp, max2lp, max_projs, &
         &                                    natm, ndim_magmom )
!
    integer i, k
    integer :: ia, i1, i2, ip
    Real(kind=DP) :: c_rho(kimg), c_mx (kimg), c_my (kimg), c_mz (kimg)
    Real(kind=DP) :: c_n11(kimg), c_n12(kimg), c_n21(kimg), c_n22(kimg)
!
    dm_magmom = 0.0d0
    Do ia=1, natm
       Do ip=1, max_projs
          Do i2=1, max2lp
             Do i1=1, max2lp

                c_n11(1) = real( dm_ssrep( i1, i2, ip, ia, 1 ) )
                c_n12(1) = real( dm_ssrep( i1, i2, ip, ia, 2 ) )
                c_n21(1) = real( dm_ssrep( i1, i2, ip, ia, 3 ) )
                c_n22(1) = real( dm_ssrep( i1, i2, ip, ia, 4 ) )

                c_n11(2) = aimag( dm_ssrep( i1, i2, ip, ia, 1 ) )
                c_n12(2) = aimag( dm_ssrep( i1, i2, ip, ia, 2 ) )
                c_n21(2) = aimag( dm_ssrep( i1, i2, ip, ia, 3 ) )
                c_n22(2) = aimag( dm_ssrep( i1, i2, ip, ia, 4 ) )

! -------------------- real part ------
                c_rho(1) = c_n11(1) + c_n22(1)
                c_mx(1)  = c_n12(1) + c_n21(1)
                c_my(1)  = c_n21(2) - c_n12(2)
                c_mz(1)  = c_n11(1) - c_n22(1)
! ---------------------- imaginary part --
                c_rho(2) = c_n11(2) + c_n22(2)
                c_mx(2)  = c_n12(2) + c_n21(2)
                c_my(2)  =-c_n21(1) + c_n12(1)
                c_mz(2)  = c_n11(2) - c_n22(2)
!
                dm_magmom( i1, i2, ip, ia, 1 )  = c_rho(1)
                dm_magmom( i1, i2, ip, ia, 2 )  = c_mx(1)
                dm_magmom( i1, i2, ip, ia, 3 )  = c_my(1)
                dm_magmom( i1, i2, ip, ia, 4 )  = c_mz(1)
             End Do
          End do
       End do
! --
    End Do
  end subroutine m_ES_DensMat_To_MagMom_dm

  subroutine m_ES_MagMom_To_DensMat_porb0( nsize, porb_r_magmom, &
       &                                  porb_ssrep )

    integer, intent(in) :: nsize
    real(kind=DP), intent(in)  :: porb_r_magmom( nsize, ndim_magmom )
!    real(kind=DP), intent(in)  :: porb_i_magmom( nsize, ndim_magmom )
    complex(kind=CMPLDP), intent(out)  :: porb_ssrep( nsize, ndim_chgpot )
!
    integer :: i
    Real(kind=DP) :: c_rho(kimg), c_mx(kimg), c_my(kimg), c_mz(kimg)
    Real(kind=DP) :: c_n11(kimg), c_n12(kimg), c_n21(kimg), c_n22(kimg)
!
    porb_ssrep = 0.0d0

    Do i=1, nsize
       c_rho(1) = porb_r_magmom( i,1 )
       c_mx(1)  = porb_r_magmom( i,2 )
       c_my(1)  = porb_r_magmom( i,3 )
       c_mz(1)  = porb_r_magmom( i,4 )

!       c_rho(2) = porb_i_magmom( i,1 )
!       c_mx(2)  = porb_i_magmom( i,2 )
!       c_my(2)  = porb_i_magmom( i,3 )
!       c_mz(2)  = porb_i_magmom( i,4 )
       c_rho(2) = 0.0d0
       c_mx(2)  = 0.0d0
       c_my(2)  = 0.0d0
       c_mz(2)  = 0.0d0

       ! -------------------- real part ------
       c_n11(1) = ( c_rho(1) + c_mz(1) ) /2.0d0
       c_n22(1) = ( c_rho(1) - c_mz(1) ) /2.0d0
       c_n12(1) = ( c_mx(1)  + c_my(2) ) /2.0d0
       c_n21(1) = ( c_mx(1)  - c_my(2) ) /2.0d0
       ! ---------------------- imaginary part --
       c_n11(2) = ( c_rho(2) + c_mz(2) ) /2.0d0
       c_n22(2) = ( c_rho(2) - c_mz(2) ) /2.0d0
       c_n12(2) = ( c_mx(2)  - c_my(1) ) /2.0d0
       c_n21(2) = ( c_mx(2)  + c_my(1) ) /2.0d0
!
       porb_ssrep( i, 1 ) = dcmplx( c_n11(1), c_n11(2) )
       porb_ssrep( i, 2 ) = dcmplx( c_n12(1), c_n12(2) )
       porb_ssrep( i, 3 ) = dcmplx( c_n21(1), c_n21(2) )
       porb_ssrep( i, 4 ) = dcmplx( c_n22(1), c_n22(2) )
! --
    End Do
  end subroutine m_ES_MagMom_To_DensMat_porb0

  subroutine m_ES_MagMom_To_DensMat_porb( nsize, porb_r_magmom, porb_i_magmom, &
       &                                  porb_ssrep )

    integer, intent(in) :: nsize
    real(kind=DP), intent(in)  :: porb_r_magmom( nsize, ndim_magmom )
    real(kind=DP), intent(in)  :: porb_i_magmom( nsize, ndim_magmom )
    complex(kind=CMPLDP), intent(out)  :: porb_ssrep( nsize, ndim_chgpot )
!
    integer :: i
    Real(kind=DP) :: c_rho(kimg), c_mx(kimg), c_my(kimg), c_mz(kimg)
    Real(kind=DP) :: c_n11(kimg), c_n12(kimg), c_n21(kimg), c_n22(kimg)
!
    porb_ssrep = 0.0d0

    Do i=1, nsize
       c_rho(1) = porb_r_magmom( i,1 )
       c_mx(1)  = porb_r_magmom( i,2 )
       c_my(1)  = porb_r_magmom( i,3 )
       c_mz(1)  = porb_r_magmom( i,4 )

       c_rho(2) = porb_i_magmom( i,1 )
       c_mx(2)  = porb_i_magmom( i,2 )
       c_my(2)  = porb_i_magmom( i,3 )
       c_mz(2)  = porb_i_magmom( i,4 )

       ! -------------------- real part ------
       c_n11(1) = ( c_rho(1) + c_mz(1) ) /2.0d0
       c_n22(1) = ( c_rho(1) - c_mz(1) ) /2.0d0
       c_n12(1) = ( c_mx(1)  + c_my(2) ) /2.0d0
       c_n21(1) = ( c_mx(1)  - c_my(2) ) /2.0d0
       ! ---------------------- imaginary part --
       c_n11(2) = ( c_rho(2) + c_mz(2) ) /2.0d0
       c_n22(2) = ( c_rho(2) - c_mz(2) ) /2.0d0
       c_n12(2) = ( c_mx(2)  - c_my(1) ) /2.0d0
       c_n21(2) = ( c_mx(2)  + c_my(1) ) /2.0d0
!
       porb_ssrep( i, 1 ) = dcmplx( c_n11(1), c_n11(2) )
       porb_ssrep( i, 2 ) = dcmplx( c_n12(1), c_n12(2) )
       porb_ssrep( i, 3 ) = dcmplx( c_n21(1), c_n21(2) )
       porb_ssrep( i, 4 ) = dcmplx( c_n22(1), c_n22(2) )
! --
    End Do
  end subroutine m_ES_MagMom_To_DensMat_porb

  subroutine m_ES_DensMat_To_MagMom_porb( nsize, porb_ssrep, porb_magmom )

    integer, intent(in) :: nsize
    complex(kind=CMPLDP), intent(in)  :: porb_ssrep( nsize, ndim_chgpot )
    real(kind=DP), intent(out)  :: porb_magmom( nsize, ndim_magmom )
!
    integer :: i
    Real(kind=DP) :: c_rho(kimg), c_mx (kimg), c_my (kimg), c_mz (kimg)
    Real(kind=DP) :: c_n11(kimg), c_n12(kimg), c_n21(kimg), c_n22(kimg)
!
    porb_magmom = 0.0d0
    Do i=1, nsize
       c_n11(1) = real( porb_ssrep( i, 1 ) )
       c_n12(1) = real( porb_ssrep( i, 2 ) )
       c_n21(1) = real( porb_ssrep( i, 3 ) )
       c_n22(1) = real( porb_ssrep( i, 4 ) )

       c_n11(2) = aimag( porb_ssrep( i, 1 ) )
       c_n12(2) = aimag( porb_ssrep( i, 2 ) )
       c_n21(2) = aimag( porb_ssrep( i, 3 ) )
       c_n22(2) = aimag( porb_ssrep( i, 4 ) )

! -------------------- real part ------
       c_rho(1) = c_n11(1) + c_n22(1)
       c_mx(1)  = c_n12(1) + c_n21(1)
       c_my(1)  = c_n21(2) - c_n12(2)
       c_mz(1)  = c_n11(1) - c_n22(1)
! ---------------------- imaginary part --
       c_rho(2) = c_n11(2) + c_n22(2)
       c_mx(2)  = c_n12(2) + c_n21(2)
       c_my(2)  =-c_n21(1) + c_n12(1)
       c_mz(2)  = c_n11(2) - c_n22(2)
!
       porb_magmom( i, 1 ) = c_rho(1)
       porb_magmom( i, 2 ) = c_mx(1)
       porb_magmom( i, 3 ) = c_my(1)
       porb_magmom( i, 4 ) = c_mz(1)
! --
    End Do
  end subroutine m_ES_DensMat_To_MagMom_porb

  subroutine m_ES_MagMom_to_DensMat_DionPaw( dion_magmom, dion_ssrep )
    real(kind=DP), intent(in):: dion_magmom(nlmt,nlmt,ndim_magmom,natm)
    complex(kind=CMPLDP), intent(out) :: dion_ssrep(nlmt,nlmt,ndim_chgpot,natm)
!
    integer :: ia, lmt1, lmt2
    Real(kind=DP) :: c_rho, c_mx, c_my, c_mz
    Real(kind=DP) :: c_n11(kimg), c_n12(kimg), c_n21(kimg), c_n22(kimg)

    dion_ssrep = 0.0d0
    Do ia=1, natm
       Do lmt2=1, nlmt
          Do lmt1=1, nlmt
             c_rho = dion_magmom( lmt1, lmt2, 1, ia )
             c_mx  = dion_magmom( lmt1, lmt2, 2, ia )
             c_my  = dion_magmom( lmt1, lmt2, 3, ia )
             c_mz  = dion_magmom( lmt1, lmt2, 4, ia )
             ! -------------------- real part ------
             c_n11(1) = c_rho + c_mz
             c_n22(1) = c_rho - c_mz
             c_n12(1) = c_mx
             c_n21(1) = c_mx
       ! ---------------------- imaginary part --
             c_n11(2) = 0.0d0
             c_n22(2) = 0.0d0
             c_n12(2) = -c_my
             c_n21(2) =  c_my
       ! ----
             dion_ssrep( lmt1, lmt2, 1, ia ) = dcmplx( c_n11(1), c_n11(2) )
             dion_ssrep( lmt1, lmt2, 2, ia ) = dcmplx( c_n12(1), c_n12(2) )
             dion_ssrep( lmt1, lmt2, 3, ia ) = dcmplx( c_n21(1), c_n21(2) )
             dion_ssrep( lmt1, lmt2, 4, ia ) = dcmplx( c_n22(1), c_n22(2) )
          End do
       End Do
! --
    End Do

  end subroutine m_ES_MagMom_to_DensMat_DionPaw

  subroutine m_ES_DensMat_to_MagMom_Dhub0( dhub_ssrep, dhub_magmom )
    complex(kind=CMPLDP), intent(in) :: dhub_ssrep(nlmt,nlmt,natm,ndim_chgpot)
    real(kind=DP), intent(out):: dhub_magmom(nlmt,nlmt,natm,ndim_magmom)
!
    integer :: ia, lmt1, lmt2
    Real(kind=DP) :: c_rho(kimg), c_mx(kimg), c_my(kimg), c_mz(kimg)
    Real(kind=DP) :: c_n11(kimg), c_n12(kimg), c_n21(kimg), c_n22(kimg)
!
    dhub_magmom = 0.0d0

    Do ia=1, natm
       Do lmt2=1, nlmt
          Do lmt1=1, nlmt
             c_n11(1) = real( dhub_ssrep( lmt1, lmt2, ia, 1 ) )
             c_n12(1) = real( dhub_ssrep( lmt1, lmt2, ia, 2 ) )
             c_n21(1) = real( dhub_ssrep( lmt1, lmt2, ia, 3 ) )
             c_n22(1) = real( dhub_ssrep( lmt1, lmt2, ia, 4 ) )

             c_n11(2) = aimag( dhub_ssrep( lmt1, lmt2, ia, 1 ) )
             c_n12(2) = aimag( dhub_ssrep( lmt1, lmt2, ia, 2 ) )
             c_n21(2) = aimag( dhub_ssrep( lmt1, lmt2, ia, 3 ) )
             c_n22(2) = aimag( dhub_ssrep( lmt1, lmt2, ia, 4 ) )

! -------------------- real part ------
             c_rho(1) = ( c_n11(1) + c_n22(1) ) /2.0d0
             c_mx(1)  = ( c_n12(1) + c_n21(1) ) /2.0d0
             c_my(1)  = ( c_n21(2) - c_n12(2) ) /2.0d0
             c_mz(1)  = ( c_n11(1) - c_n22(1) ) /2.0d0
! ---------------------- imaginary part --
             c_rho(2) = ( c_n11(2) + c_n22(2) ) /2.0d0
             c_mx(2)  = ( c_n12(2) + c_n21(2) ) /2.0d0
             c_my(2)  = (-c_n21(1) + c_n12(1) ) /2.0d0
             c_mz(2)  = ( c_n11(2) - c_n22(2) ) /2.0d0
! ----
             dhub_magmom( lmt1, lmt2, ia, 1 ) = c_rho(1)
             dhub_magmom( lmt1, lmt2, ia, 2 ) = c_mx(1)
             dhub_magmom( lmt1, lmt2, ia, 3 ) = c_my(1)
             dhub_magmom( lmt1, lmt2, ia, 4 ) = c_mz(1)
!
!             write(500,*) 'ilmt1 ilmt2 ia ', lmt1, lmt2, ia
!             write(500,*) c_rho(2), c_mx(2), c_my(2), c_mz(2)
          End do
       End Do
! --
    End Do

  end subroutine m_ES_DensMat_to_MagMom_Dhub0

  subroutine m_ES_DensMat_to_MagMom_Dhub( dhub_ssrep, dhub_r_magmom, dhub_i_magmom )
    complex(kind=CMPLDP), intent(in) :: dhub_ssrep(nlmt,nlmt,natm,ndim_chgpot)
    real(kind=DP), intent(out):: dhub_r_magmom(nlmt,nlmt,natm,ndim_magmom)
    real(kind=DP), intent(out):: dhub_i_magmom(nlmt,nlmt,natm,ndim_magmom)
!
    integer :: ia, lmt1, lmt2
    Real(kind=DP) :: c_rho(kimg), c_mx(kimg), c_my(kimg), c_mz(kimg)
    Real(kind=DP) :: c_n11(kimg), c_n12(kimg), c_n21(kimg), c_n22(kimg)
!
    dhub_r_magmom = 0.0d0
    dhub_i_magmom = 0.0d0

    Do ia=1, natm
       Do lmt2=1, nlmt
          Do lmt1=1, nlmt
             c_n11(1) = real( dhub_ssrep( lmt1, lmt2, ia, 1 ) )
             c_n12(1) = real( dhub_ssrep( lmt1, lmt2, ia, 2 ) )
             c_n21(1) = real( dhub_ssrep( lmt1, lmt2, ia, 3 ) )
             c_n22(1) = real( dhub_ssrep( lmt1, lmt2, ia, 4 ) )

             c_n11(2) = aimag( dhub_ssrep( lmt1, lmt2, ia, 1 ) )
             c_n12(2) = aimag( dhub_ssrep( lmt1, lmt2, ia, 2 ) )
             c_n21(2) = aimag( dhub_ssrep( lmt1, lmt2, ia, 3 ) )
             c_n22(2) = aimag( dhub_ssrep( lmt1, lmt2, ia, 4 ) )

! -------------------- real part ------
             c_rho(1) = ( c_n11(1) + c_n22(1) ) /2.0d0
             c_mx(1)  = ( c_n12(1) + c_n21(1) ) /2.0d0
             c_my(1)  = ( c_n21(2) - c_n12(2) ) /2.0d0
             c_mz(1)  = ( c_n11(1) - c_n22(1) ) /2.0d0
! ---------------------- imaginary part --
             c_rho(2) = ( c_n11(2) + c_n22(2) ) /2.0d0
             c_mx(2)  = ( c_n12(2) + c_n21(2) ) /2.0d0
             c_my(2)  = (-c_n21(1) + c_n12(1) ) /2.0d0
             c_mz(2)  = ( c_n11(2) - c_n22(2) ) /2.0d0
! ----
             dhub_r_magmom( lmt1, lmt2, ia, 1 ) = c_rho(1)
             dhub_r_magmom( lmt1, lmt2, ia, 2 ) = c_mx(1)
             dhub_r_magmom( lmt1, lmt2, ia, 3 ) = c_my(1)
             dhub_r_magmom( lmt1, lmt2, ia, 4 ) = c_mz(1)
! ----
             dhub_i_magmom( lmt1, lmt2, ia, 1 ) = c_rho(2)
             dhub_i_magmom( lmt1, lmt2, ia, 2 ) = c_mx(2)
             dhub_i_magmom( lmt1, lmt2, ia, 3 ) = c_my(2)
             dhub_i_magmom( lmt1, lmt2, ia, 4 ) = c_mz(2)
!
!             write(500,*) 'ilmt1 ilmt2 ia ', lmt1, lmt2, ia
!             write(500,*) c_rho(2), c_mx(2), c_my(2), c_mz(2)
          End do
       End Do
! --
    End Do

  end subroutine m_ES_DensMat_to_MagMom_Dhub

  subroutine m_ES_MagMom_to_DensMat_Dhub( dhub_r_magmom, dhub_i_magmom, dhub_ssrep )
    real(kind=DP), intent(in):: dhub_r_magmom(nlmt,nlmt,natm,ndim_magmom)
    real(kind=DP), intent(in):: dhub_i_magmom(nlmt,nlmt,natm,ndim_magmom)
    complex(kind=CMPLDP), intent(out) :: dhub_ssrep(nlmt,nlmt,natm,ndim_chgpot)
!
    integer :: ia, lmt1, lmt2
    Real(kind=DP) :: c_rho(kimg), c_mx(kimg), c_my(kimg), c_mz(kimg)
    Real(kind=DP) :: c_n11(kimg), c_n12(kimg), c_n21(kimg), c_n22(kimg)
!
    dhub_ssrep = 0.0d0

    Do ia=1, natm
       Do lmt2=1, nlmt
          Do lmt1=1, nlmt

             c_rho(1) = dhub_r_magmom( lmt1, lmt2, ia, 1 )
             c_mx(1)  = dhub_r_magmom( lmt1, lmt2, ia, 2 )
             c_my(1)  = dhub_r_magmom( lmt1, lmt2, ia, 3 )
             c_mz(1)  = dhub_r_magmom( lmt1, lmt2, ia, 4 )
!
             c_rho(2) = dhub_i_magmom( lmt1, lmt2, ia, 1 )
             c_mx(2)  = dhub_i_magmom( lmt1, lmt2, ia, 2 )
             c_my(2)  = dhub_i_magmom( lmt1, lmt2, ia, 3 )
             c_mz(2)  = dhub_i_magmom( lmt1, lmt2, ia, 4 )

       ! -------------------- real part ------
             c_n11(1) = c_rho(1) + c_mz(1)
             c_n22(1) = c_rho(1) - c_mz(1)
             c_n12(1) = c_mx(1)  + c_my(2)
             c_n21(1) = c_mx(1)  - c_my(2)
       ! ---------------------- imaginary part --
             c_n11(2) = c_rho(2) + c_mz(2)
             c_n22(2) = c_rho(2) - c_mz(2)
             c_n12(2) = c_mx(2)  - c_my(1)
             c_n21(2) = c_mx(2)  + c_my(1)
       ! ----
             dhub_ssrep( lmt1, lmt2, ia, 1 ) = dcmplx( c_n11(1), c_n11(2) )
             dhub_ssrep( lmt1, lmt2, ia, 2 ) = dcmplx( c_n12(1), c_n12(2) )
             dhub_ssrep( lmt1, lmt2, ia, 3 ) = dcmplx( c_n21(1), c_n21(2) )
             dhub_ssrep( lmt1, lmt2, ia, 4 ) = dcmplx( c_n22(1), c_n22(2) )

          End do
       End Do
! --
    End Do

  end subroutine m_ES_MagMom_to_DensMat_Dhub

! ------------------------------------------------------------------------
!!
!!!   Diagonalize the density matrix, obtain the local quantization axis,
!!!                and up/down spin desity along this axis
!!
! -----------------------------------------------------------------------
  subroutine m_ES_get_Angles_MagMom_Rspace( RhoMag_R, Rot_Angles )
    Real(kind=DP), intent(in)  :: RhoMag_R(    ista_fftph:iend_fftph,ndim_magmom )
    Real(kind=DP), intent(out) :: Rot_Angles( ista_fftph:iend_fftph,2 )

    integer :: i
    Real(kind=DP) :: c_n, c_mx, c_my, c_mz
    Real(kind=DP) :: theta, phi, m_norm

!    real(kind=DP) :: m_norm_criteria = 5.0D-3
!    real(kind=DP) :: m_norm_criteria = 1.0D-6
!    real(kind=DP) :: m_norm_criteria = 1.0D-12
    real(kind=DP) :: m_norm_criteria = 1.0D-20

! --
    real(kind=DP) :: theta_neighbor, phi_neighbor

! --------- initialize --
    theta_neighbor = 0.0d0
    phi_neighbor = 0.0d0
! --
!    write(940,*) 'aa '

    Do i=ista_fftph, iend_fftph
       c_n  = RhoMag_R( i,1 );  c_mx = RhoMag_R( i,2 )
       c_my = RhoMag_R( i,3 );  c_mz = RhoMag_R( i,4 )
!
       m_norm = sqrt( c_mx**2 + c_my**2 + c_mz**2 )
!
       if ( m_norm > m_norm_criteria ) then

! ---------------- original ------
!!!!!!!!          theta = atan2( sqrt( c_mx**2 +c_my**2), c_mz )
!!!!!          phi = atan2(  c_my, c_mx )

! --- for test -
          theta = atan2( sqrt( c_mx**2 +c_my**2), abs(c_mz) )
          phi = atan2( abs(c_my), abs(c_mx) )

!          goto 100
          goto 120

! --- for test2 -

          if ( c_mz > 0.0d0 ) then
             theta = atan2( sqrt( c_mx**2 +c_my**2), c_mz )
             phi = atan2( c_my, c_mx )
          else if ( c_mz < 0.0d0 ) then
             theta = atan2( sqrt( c_mx**2 +c_my**2), -c_mz )
             phi = atan2( -c_my, -c_mx )
          else
             theta = atan2( sqrt( c_mx**2 +c_my**2), -c_mz )
             if ( c_mx > 0.0d0 ) then
                phi = atan2( c_my,  c_mx )
             else if ( c_mx < 0.0d0 ) then
                phi = atan2( -c_my,  -c_mx )
             else
                phi = atan2( abs(c_my), c_mx )
             endif
          endif

120       continue
! --- for test3 -

!          write(940,*) i, c_n
          if ( c_mz >= 0.0d0 ) then
             theta = atan2( sqrt( c_mx**2 +c_my**2), c_mz )
             phi = atan2( c_my, c_mx )
          else if ( c_mz < 0.0d0 ) then
             theta = atan2( sqrt( c_mx**2 +c_my**2), -c_mz )
             phi = atan2( -c_my, -c_mx )
          endif

100       continue

       else
!!!!!!          theta = 0.0d0;       phi =0.0d0
          theta = theta_neighbor;  phi  = phi_neighbor
       endif
!
       theta_neighbor = theta;  phi_neighbor = phi
!!
       Rot_Angles( i, 1 ) = theta
       Rot_Angles( i, 2 ) = phi
    End do
  end subroutine m_ES_get_Angles_MagMom_Rspace

  subroutine m_ES_get_Angles_MagMom_Rspace2( RhoMag_R, Rot_Angles )
    Real(kind=DP), intent(in)  :: RhoMag_R(    ista_fftph:iend_fftph,ndim_magmom )
    Real(kind=DP), intent(out) :: Rot_Angles( ista_fftph:iend_fftph,2 )

    integer :: i
    Real(kind=DP) :: c_n, c_mx, c_my, c_mz
    Real(kind=DP) :: mmx, mmy, mmz, ctmp
    Real(kind=DP) :: theta, phi, m_norm

!    real(kind=DP) :: m_norm_criteria = 5.0D-3
!    real(kind=DP) :: m_norm_criteria = 1.0D-6
!    real(kind=DP) :: m_norm_criteria = 1.0D-12
    real(kind=DP) :: m_norm_criteria = 1.0D-20
! --
    real(kind=DP) :: theta_neighbor, phi_neighbor
    real(kind=DP) :: theta0, phi0

    real(kind=DP) :: MagDirec_OnMesh_avg(3)

! --------- initialize --
    theta_neighbor = 0.0d0
    phi_neighbor = 0.0d0
! --
!    write(940,*) 'aa '

    call calc_averaged_mag_direc_now( MagDirec_OnMesh_avg )
!
    theta0 = atan2( sqrt( MagDirec_OnMesh_avg(1)**2 &
         &             + MagDirec_OnMesh_avg(2)**2 ), MagDirec_OnMesh_avg(3) )
    phi0 = atan2( MagDirec_OnMesh_avg(2), MagDirec_OnMesh_avg(1) )
!
!    write(990,*) 'th0 ph0 = ', theta0, phi0
!    write(990,*) MagDirec_OnMesh_Avg(1), MagDirec_OnMesh_Avg(2), MagDirec_OnMesh_Avg(3)
!
    Do i=ista_fftph, iend_fftph
       c_n  = RhoMag_R( i,1 );  c_mx = RhoMag_R( i,2 )
       c_my = RhoMag_R( i,3 );  c_mz = RhoMag_R( i,4 )
!
       m_norm = sqrt( c_mx**2 + c_my**2 + c_mz**2 )
!
       if ( m_norm > m_norm_criteria ) then
!
          mmx = c_mx / m_norm
          mmy = c_my / m_norm
          mmz = c_mz / m_norm

          ctmp = mmx *MagDirec_OnMesh_avg(1) &
               & + mmy *MagDirec_OnMesh_avg(2) &
               & + mmz *MagDirec_OnMesh_avg(3)

          theta = atan2( sqrt( c_mx**2 +c_my**2), c_mz )
          phi = atan2( c_my, c_mx )

          if ( ctmp >= 0.0d0 ) then
          else
!             theta = atan2( sqrt( c_mx**2 +c_my**2), -c_mz )
!             phi = atan2( -c_my, -c_mx )
             theta = PAI - theta
             phi =   PAI + phi
          endif

       else
!          theta = theta0
!          phi = phi0
          theta = theta_neighbor
          phi = phi_neighbor
       endif

       theta_neighbor = theta;  phi_neighbor = phi
!!
       Rot_Angles( i, 1 ) = theta
       Rot_Angles( i, 2 ) = phi
    End do


  contains

    subroutine calc_averaged_mag_direc_now( MagDirec_OnMesh_avg )
      real(kind=DP), intent(out) :: MagDirec_OnMesh_avg(3)

      real(kind=DP) :: MagVectorSum(3), MagVectorSum_mpi(3)
      real(kind=DP) :: cnorm
      integer :: i

      MagDirec_OnMesh_avg = 0.0d0

      MagVectorSum = 0.0d0
      Do i=ista_fftph, iend_fftph
         MagVectorSum(1) = MagVectorSum(1) +  abs( RhoMag_R( i,2 ) )
         MagVectorSum(2) = MagVectorSum(2) +  abs( RhoMag_R( i,3 ) )
         MagVectorSum(3) = MagVectorSum(3) +  abs( RhoMag_R( i,4 ) )
      End do

      if ( npes > 1 ) then
         call mpi_allreduce( MagVectorSum, MagVectorSum_mpi, 3, &
              &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
         MagVectorSum = MagVectorSum_mpi
      endif
!
      cnorm = MagVectorSum(1)**2 + MagVectorSum(2)**2 + MagVectorSum(3)**2
      cnorm = sqrt( cnorm )
!
      if ( cnorm > m_norm_criteria ) then
         MagDirec_OnMesh_avg = MagVectorSum / cnorm
      else
         MagDirec_OnMesh_avg(3) = 1.0d0
      endif
    end subroutine calc_averaged_mag_direc_now

  end subroutine m_ES_get_Angles_MagMom_Rspace2


  subroutine m_ES_SpinDens_Along_QuantzAxis( RhoMag_R, Rot_Angles, chgrhr_l )
    Real(kind=DP), intent(in)  :: RhoMag_R(   ista_fftph:iend_fftph,ndim_magmom )
    Real(kind=DP), intent(in)  :: Rot_Angles( ista_fftph:iend_fftph,nspin )
    Real(kind=DP), intent(out) :: chgrhr_l(   ista_fftph:iend_fftph,nspin )

    integer :: i
    Real(kind=DP) :: c_n, c_mx, c_my, c_mz, c_tmp
    Real(kind=DP) :: theta, phi
    Real(kind=DP) :: cos_th, sin_th, cos_ph, sin_ph
    !
    chgrhr_l = 0.0d0
    Do i=ista_fftph, iend_fftph
       c_n  = RhoMag_R( i,1 );      c_mx = RhoMag_R( i,2 )
       c_my = RhoMag_R( i,3 );      c_mz = RhoMag_R( i,4 )
!
       theta = Rot_Angles( i,1 );    phi   = Rot_Angles( i,2 )
!
       cos_th = cos( theta );  sin_th = sin( theta )
       cos_ph = cos( phi );    sin_ph = sin( phi )
       !
!!!!!!!!!!!!!!!!!       c_tmp = ( c_mx *cos_ph - c_my *sin_ph ) /2.0d0
       c_tmp = ( c_mx *cos_ph + c_my *sin_ph ) /2.0d0
       c_tmp = c_tmp *sin_th
!
       chgrhr_l( i,1 ) = ( c_n + c_mz *cos_th ) /2.0d0  + c_tmp       ! up
       chgrhr_l( i,2 ) = ( c_n - c_mz *cos_th ) /2.0d0  - c_tmp       ! down
! --
!       write(*,*) "  A ", i, theta, phi
!       write(*,*) "  B ", i, c_mz

    End do
!
    if ( iprixc >=2 ) then
       write(nfout,*) '!ES XC =='
       write(nfout,*) '!ES XC :         ispin = 1 '
       write(nfout,*) '!ES XC -- chgrhr_l -- '
       write(nfout,'(" !ES XC ",6f12.6)') (chgrhr_l(i,1),i=ista_fftph,ista_fftph+17)
!
       write(nfout,*) '!ES XC =='
       write(nfout,*) '!ES XC :         ispin = 2 '
       write(nfout,*) '!ES XC -- chgrhr_l -- '
       write(nfout,'(" !ES XC ",6f12.6)') (chgrhr_l(i,2),i=ista_fftph,ista_fftph+17)

    endif
!
  end subroutine m_ES_SpinDens_Along_QuantzAxis

  subroutine m_ES_SpinDens_Along_QuantzAxis2( RhoMag_R, chgrhr_l, &
       &                                      quantz_axis_inversion_flg_mesh, &
       &                                      skip_flag )
    Real(kind=DP), intent(in)  :: RhoMag_R(   ista_fftph:iend_fftph,ndim_magmom )
    Real(kind=DP), intent(out) :: chgrhr_l(   ista_fftph:iend_fftph,nspin )
    integer, intent(inout) :: quantz_axis_inversion_flg_mesh(ista_fftph:iend_fftph)
    logical, intent(in) :: skip_flag

    integer :: i
    Real(kind=DP) :: c_n, c_mx, c_my, c_mz, ctmp
    real(kind=DP) :: d1, d2
    real(kind=DP) :: MagDirec_OnMesh_avg(3)
    real(kind=DP) :: mmx, mmy, mmz, m_norm
    real(kind=DP) :: m_norm_criteria = 1.0D-20
!    real(kind=DP) :: m_norm_lowlimit = 1.0D-40
    real(kind=DP) :: m_norm_lowlimit = 1.0D-30

    integer :: exit_flag
    integer :: count1, count1_sum
    integer :: count2, count2_sum

! ============
!    call calc_averaged_mag_direc_now( MagDirec_OnMesh_avg )

    call calc_averaged_mag_direc_now_A( Global_Quantz_axis_now, exit_flag )
    if ( exit_flag == 1 ) then
       call calc_averaged_mag_direc_now_B( Global_Quantz_axis_now )
    endif
! ===========

    if ( iprixc >=2 ) then
       write(nfout,'(A,3E20.10)') "!ES XC : GlobalAxis = ", Global_Quantz_axis_now(1:3)
    endif

    chgrhr_l = 0.0d0
    Do i=ista_fftph, iend_fftph
       c_n  = RhoMag_R( i,1 );      c_mx = RhoMag_R( i,2 )
       c_my = RhoMag_R( i,3 );      c_mz = RhoMag_R( i,4 )
!
       m_norm = sqrt( c_mx**2 + c_my**2 + c_mz**2 )
       d1 = (c_n +m_norm) /2.0d0;  d2 = ( c_n - m_norm )/2.0d0

! ---- only test, but may be important ------
       if ( skip_flag ) goto 130

!       if ( c_n < 0.0 ) then
!          d1 = m_norm_lowlimit;    d2 = m_norm_lowlimit
!       else
          if ( d1 < 0.0 ) d1 = m_norm_lowlimit
          if ( d2 < 0.0 ) d2 = m_norm_lowlimit
!       endif

130    continue
! -----------------------------------------
       if ( m_norm > m_norm_criteria ) then
!
!          mmx = c_mx / m_norm
!          mmy = c_my / m_norm
!         mmz = c_mz / m_norm
          mmx = c_mx
          mmy = c_my
          mmz = c_mz

!          ctmp = mmx *MagDirec_OnMesh_avg(1) &
!               & + mmy *MagDirec_OnMesh_avg(2) &
!               & + mmz *MagDirec_OnMesh_avg(3)
          ctmp = mmx *Global_Quantz_Axis_now(1) &
               & + mmy *Global_Quantz_Axis_now(2) &
               & + mmz *Global_Quantz_Axis_now(3)

          if ( ctmp > 0.0 ) then
             chgrhr_l( i,1 ) = d1;  chgrhr_l( i,2 ) = d2
             quantz_axis_inversion_flg_mesh(i) = 1
          else
             chgrhr_l( i,1 ) = d2;  chgrhr_l( i,2 ) = d1
             quantz_axis_inversion_flg_mesh(i) = -1
          endif
!
       else
!!!!          d1 = ( d1 +d2 )/2.0d0;  d2 = d1
          chgrhr_l( i,1 ) = d1;  chgrhr_l( i,2 ) = d2
          quantz_axis_inversion_flg_mesh(i) = 1
       endif

    End do
!
    if ( iprixc >=2 ) then
       write(nfout,*) '!ES XC =='
       write(nfout,*) '!ES XC :         ispin = 1 '
       write(nfout,*) '!ES XC -- chgrhr_l -- '
       write(nfout,'(" !ES XC ",6f12.6)') (chgrhr_l(i,1),i=ista_fftph,ista_fftph+17)
!
       write(nfout,*) '!ES XC =='
       write(nfout,*) '!ES XC :         ispin = 2 '
       write(nfout,*) '!ES XC -- chgrhr_l -- '
       write(nfout,'(" !ES XC ",6f12.6)') (chgrhr_l(i,2),i=ista_fftph,ista_fftph+17)
    endif
!
    if ( iprixc >=2 ) then
       count1 = 0;  count2 = 0
       Do i=ista_fftph, iend_fftph
          if ( quantz_axis_inversion_flg_mesh(i) ==  1 ) count1 = count1 +1
          if ( quantz_axis_inversion_flg_mesh(i) == -1 ) count2 = count2 +1
       Enddo
       call mpi_allreduce( count1, count1_sum, 1, mpi_integer, &
            &              mpi_sum, MPI_CommGroup, ierr )
       call mpi_allreduce( count2, count2_sum, 1, mpi_integer, &
            &              mpi_sum, MPI_CommGroup, ierr )
       write(nfout,*) 'count_sum = ', count1_sum, count2_sum
    endif

  contains

    subroutine calc_averaged_mag_direc_now_A( MagDirec_OnMesh_avg, exit_flag )
      real(kind=DP), intent(out) :: MagDirec_OnMesh_avg(3)
      integer, intent(out) :: exit_flag

      real(kind=DP) :: MagVectorSum(3), MagVectorSum_mpi(3)
      real(kind=DP) :: cnorm
      integer :: i

      MagDirec_OnMesh_avg = 0.0d0

      MagVectorSum = 0.0d0
      Do i=ista_fftph, iend_fftph
         MagVectorSum(1) = MagVectorSum(1) +  RhoMag_R( i,2 )
         MagVectorSum(2) = MagVectorSum(2) +  RhoMag_R( i,3 )
         MagVectorSum(3) = MagVectorSum(3) +  RhoMag_R( i,4 )
      End do

      if ( npes > 1 ) then
         call mpi_allreduce( MagVectorSum, MagVectorSum_mpi, 3, &
              &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
         MagVectorSum = MagVectorSum_mpi
      endif
!
      cnorm = MagVectorSum(1)**2 + MagVectorSum(2)**2 + MagVectorSum(3)**2
      cnorm = sqrt( cnorm )
!
      if ( cnorm > m_norm_criteria ) then
         MagDirec_OnMesh_avg = MagVectorSum / cnorm
         exit_flag = 0
      else
         MagDirec_OnMesh_avg(3) = 1.0d0
         exit_flag = 1
      endif
    end subroutine calc_averaged_mag_direc_now_A

    subroutine calc_averaged_mag_direc_now_B( MagDirec_OnMesh_avg )
      real(kind=DP), intent(out) :: MagDirec_OnMesh_avg(3)

      real(kind=DP) :: MagVectorSum(3), MagVectorSum_mpi(3)
      real(kind=DP) :: cnorm
      integer :: i

      MagDirec_OnMesh_avg = 0.0d0

      MagVectorSum = 0.0d0
      Do i=ista_fftph, iend_fftph
         MagVectorSum(1) = MagVectorSum(1) +  abs( RhoMag_R( i,2 ) )
         MagVectorSum(2) = MagVectorSum(2) +  abs( RhoMag_R( i,3 ) )
         MagVectorSum(3) = MagVectorSum(3) +  abs( RhoMag_R( i,4 ) )
      End do

      if ( npes > 1 ) then
         call mpi_allreduce( MagVectorSum, MagVectorSum_mpi, 3, &
              &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
         MagVectorSum = MagVectorSum_mpi
      endif
!
      cnorm = MagVectorSum(1)**2 + MagVectorSum(2)**2 + MagVectorSum(3)**2
      cnorm = sqrt( cnorm )
!
      if ( cnorm > m_norm_criteria ) then
         MagDirec_OnMesh_avg = MagVectorSum / cnorm
      else
         MagDirec_OnMesh_avg(3) = 1.0d0
      endif
    end subroutine calc_averaged_mag_direc_now_B

  end subroutine m_ES_SpinDens_Along_QuantzAxis2

! -----------------------------------------------------------------------
!!
!!!    convert "up/down"-xc potential calculated along the quantization axis
!!!            to the original xyz-axis basis.
!!
! -----------------------------------------------------------------------
  subroutine m_ES_cp_VxcR_to_afft( iloop, VxcPot_R, Rot_Angles, afft )
    integer, intent(in) :: iloop
    Real(kind=DP), intent(in)  :: VxcPot_R(   ista_fftph:iend_fftph,nspin )
    Real(kind=DP), intent(in)  :: Rot_Angles( ista_fftph:iend_fftph,nspin )
    Real(kind=DP), intent(out) :: afft( ista_fftp:iend_fftp )

    integer :: i
    Real(kind=DP) :: theta, phi
    Real(kind=DP) :: V0, dV
    Real(kind=DP) :: cos_th, sin_th, cos_ph, sin_ph

    afft = 0.0d0
    Do i=ista_fftph, iend_fftph
       theta = Rot_Angles( i,1 ); phi = Rot_Angles( i,2 )
       cos_th = cos( theta );  sin_th = sin( theta )
       cos_ph = cos( phi );    sin_ph = sin( phi )
!
       V0 = ( VxcPot_R( i,1 ) + VxcPot_R(i,2) ) /2.0d0
       dV = ( VxcPot_R( i,1 ) - VxcPot_R(i,2) ) /2.0d0
!
       select case (iloop)
       case (1)
          afft( 2*i-1 ) = v0 + dV *cos_th
       case (2)
          afft( 2*i-1 ) =  dV *cos_ph *sin_th
          afft( 2*i   ) = -dV *sin_ph *sin_th
       case (3)
          afft( 2*i-1 ) =  dV *cos_ph *sin_th
          afft( 2*i   ) =  dV *sin_ph *sin_th
       case (4)
          afft( 2*i-1 ) = v0 - dV *cos_th
       end select
    End Do
  end subroutine m_ES_cp_VxcR_to_afft

  subroutine m_ES_cp_VxcR_to_afft2( iloop, quantz_axis_inversion_flg_mesh, &
       &                            RhoMag_R, VxcPot_R, afft )
                                         ! output is in magmom rep.
    integer, intent(in) :: iloop
    integer, intent(in) :: quantz_axis_inversion_flg_mesh(ista_fftph:iend_fftph)
    Real(kind=DP), intent(in)  :: RhoMag_R(   ista_fftph:iend_fftph,ndim_magmom )
    Real(kind=DP), intent(in)  :: VxcPot_R(   ista_fftph:iend_fftph,nspin )
    Real(kind=DP), intent(out) :: afft( ista_fftp:iend_fftp )

    integer :: i
    Real(kind=DP) :: V0, dV
    real(kind=DP) :: c_n, c_mx, c_my, c_mz
    real(kind=DP) :: m_norm
    real(kind=DP) :: mmx, mmy, mmz
    real(kind=DP) :: m_norm_criteria = 1.0D-20
!!    real(kind=DP) :: m_norm_criteria = 1.0D-16
!!    real(kind=DP) :: m_norm_criteria = 1.0D-18

    afft = 0.0d0
    Do i=ista_fftph, iend_fftph
       V0 = ( VxcPot_R( i,1 ) + VxcPot_R(i,2) ) /2.0d0
       dV = ( VxcPot_R( i,1 ) - VxcPot_R(i,2) ) /2.0d0
!
       c_n  = RhoMag_R( i,1 );      c_mx = RhoMag_R( i,2 )
       c_my = RhoMag_R( i,3 );      c_mz = RhoMag_R( i,4 )
!
       m_norm = sqrt( c_mx**2 + c_my**2 + c_mz**2 )

       if ( m_norm > m_norm_criteria ) then
          mmx = c_mx / m_norm
          mmy = c_my / m_norm
          mmz = c_mz / m_norm
       else
!          if ( quantz_axis_inversion_flg_mesh(i) == 1 ) then
!             mmx = Global_Quantz_axis_now(1)
!             mmy = Global_Quantz_axis_now(2)
!             mmz = Global_Quantz_axis_now(3)
!          else
!             mmx = -Global_Quantz_axis_now(1)
!             mmy = -Global_Quantz_axis_now(2)
!             mmz = -Global_Quantz_axis_now(3)
!          endif
          mmx = 0.0d0; mmy = 0.0d0;  mmz = 0.0d0
       endif

       select case (iloop)
       case (1)
          afft( 2*i-1 ) = v0
       case (2)
          afft( 2*i-1 ) =  dV *mmx *quantz_axis_inversion_flg_mesh(i)
       case (3)
          afft( 2*i-1 ) =  dV *mmy *quantz_axis_inversion_flg_mesh(i)
       case (4)
          afft( 2*i-1 ) =  dV *mmz *quantz_axis_inversion_flg_mesh(i)
       end select
    End Do
  end subroutine m_ES_cp_VxcR_to_afft2

  subroutine m_ES_Chgsoft_Along_QuantzAxis( RhoMag_R, chgrhr_l, &
       &                                    quantz_axis_inversion_flg_mesh )

    Real(kind=DP), intent(in)  :: RhoMag_R( ista_fftph:iend_fftph,ndim_magmom )
    Real(kind=DP), intent(out) :: chgrhr_l( ista_fftph:iend_fftph,nspin )
    integer, intent(in) :: quantz_axis_inversion_flg_mesh(ista_fftph:iend_fftph)

    integer :: i
    Real(kind=DP) :: c_n, c_mx, c_my, c_mz, m_norm
    real(kind=DP) :: d1, d2

    Do i=ista_fftph, iend_fftph
       c_n  = RhoMag_R( i,1 );    c_mx = RhoMag_R( i,2 )
       c_my = RhoMag_R( i,3 );    c_mz = RhoMag_R( i,4 )
!
       m_norm = sqrt( c_mx**2 + c_my**2 + c_mz**2 )
       d1 = (c_n +m_norm) /2.0d0;  d2 = ( c_n - m_norm )/2.0d0
!
       if ( quantz_axis_inversion_flg_mesh(i) == 1 ) then
          chgrhr_l(i,1) = d1
          chgrhr_l(i,2) = d2
       else
          chgrhr_l(i,1) = d2
          chgrhr_l(i,2) = d1
       endif
    End do
  end subroutine m_ES_Chgsoft_Along_QuantzAxis

  subroutine m_ES_QuantzAxis_inv_flg_atm( magmom_local, &
       &                                  quantz_axis_inversion_flg_atm )
    Real(kind=DP), intent(in)  :: magmom_local(natm,3)
    integer, intent(out) :: quantz_axis_inversion_flg_atm(natm)

    integer :: ia
    Real(kind=DP) :: c_mx, c_my, c_mz, ctmp
    real(kind=DP) :: mmx, mmy, mmz, m_norm
    real(kind=DP) :: m_norm_criteria = 1.0D-20
   !

    Do ia=1, natm
       c_mx = magmom_local( ia,1 )
       c_my = magmom_local( ia,2 )
       c_mz = magmom_local( ia,3 )
!
       m_norm = sqrt( c_mx**2 + c_my**2 + c_mz**2 )
!
       if ( m_norm > m_norm_criteria ) then
!
          mmx = c_mx / m_norm
          mmy = c_my / m_norm
          mmz = c_mz / m_norm

          ctmp = mmx *Global_Quantz_Axis_now(1) &
               & + mmy *Global_Quantz_Axis_now(2) &
               & + mmz *Global_Quantz_Axis_now(3)

          if ( ctmp > 0.0 ) then
             quantz_axis_inversion_flg_atm(ia) = 1
          else
             quantz_axis_inversion_flg_atm(ia) = -1
          endif
       else
          quantz_axis_inversion_flg_atm(ia) = 1
       endif
    End do
  end subroutine m_ES_QuantzAxis_inv_flg_atm

  subroutine m_ES_Chghard_Along_QuantzAxis( hsr, hsrd,  hsr_along_QuantzAxis, &
       &                                    hsrd_along_QuantzAxis )
!
    real(kind=DP), intent(in) :: hsr( natm,nlmt,nlmt,ndim_magmom )
    real(kind=DP), intent(in) :: hsrd( natm,nlmt,nlmt,ndim_magmom,3,3 )
    real(kind=DP), intent(out) :: hsr_along_QuantzAxis( natm,nlmt,nlmt,nspin )
   real(kind=DP), intent(out) :: hsrd_along_QuantzAxis( natm,nlmt,nlmt,nspin,3,3 )

! ------
!!!!    integer, intent(in) :: quantz_axis_inversion_flg_atm(natm)
! -----
    real(kind=DP) :: cos_th, sin_th, cos_ph, sin_ph
    real(kind=DP) :: c_n, c_mx, c_my, c_mz, c_tmp, c1
    real(kind=DP) :: d1, d2
!
    integer :: ia, it, lmt1, lmt2, i1, i2
!
    cos_th = Global_Quantz_Axis_now(3)
    sin_th = sqrt( 1.0D0 - cos_th**2 )
!
    if ( sin_th > 0.0d0 ) then
       c1 = Global_Quantz_Axis_now(1)**2 + Global_Quantz_Axis_now(2)**2
       cos_ph = Global_Quantz_Axis_now(1) / c1
       sin_ph = Global_Quantz_Axis_now(2) / c1
    else
       cos_ph = 0.0d0;  sin_ph = 0.0d0
    endif

    Do ia=1, natm
       it = ityp(ia)
       Do lmt1=1, ilmt(it)
          Do lmt2=1, ilmt(it)

             c_n =  hsr( ia, lmt1, lmt2, 1 )
             c_mx = hsr( ia, lmt1, lmt2, 2 )
             c_my = hsr( ia, lmt1, lmt2, 3 )
             c_mz = hsr( ia, lmt1, lmt2, 4 )
!
             c_tmp = ( c_mx *cos_ph + c_my *sin_ph ) /2.0d0
             c_tmp = c_tmp *sin_th
!
             d1 = ( c_n + c_mz *cos_th ) /2.0d0  + c_tmp
             d2 = ( c_n - c_mz *cos_th ) /2.0d0  - c_tmp
!
             hsr_along_QuantzAxis(ia,lmt1,lmt2,1) = d1
             hsr_along_QuantzAxis(ia,lmt1,lmt2,2) = d2
          End do
       End Do
    End do
!
    Do ia=1, natm
       Do lmt1=1, ilmt(it)
          Do lmt2=1, ilmt(it)

             Do i1=1, 3
                Do i2=1, 3
                   c_n =  hsrd( ia, lmt1, lmt2, 1, i1, i2 )
                   c_mx = hsrd( ia, lmt1, lmt2, 2, i1, i2 )
                   c_my = hsrd( ia, lmt1, lmt2, 3, i1, i2 )
                   c_mz = hsrd( ia, lmt1, lmt2, 4, i1, i2 )
!
                   c_tmp = ( c_mx *cos_ph + c_my *sin_ph ) /2.0d0
                   c_tmp = c_tmp *sin_th
!
                   d1 = ( c_n + c_mz *cos_th ) /2.0d0  + c_tmp
                   d2 = ( c_n - c_mz *cos_th ) /2.0d0  - c_tmp
!
                   hsrd_along_QuantzAxis(ia,lmt1,lmt2,1,i1,i2) = d1
                   hsrd_along_QuantzAxis(ia,lmt1,lmt2,2,i1,i2) = d2
                End do
             End Do
          End do
       End Do
    End do

  end subroutine m_ES_Chghard_Along_QuantzAxis

! ------
  subroutine m_ES_GradCorr_Along_QuantzAxis( RhoMag_R, bfft_kt, &
       &                                    quantz_axis_inversion_flg_mesh )
    real(kind=DP), intent(in) :: RhoMag_R( ista_fftph:iend_fftph,ndim_magmom )
    real(kind=DP), intent(inout) :: bfft_kt( ista_fftp:iend_fftp,ndim_magmom )
    integer, intent(in) :: quantz_axis_inversion_flg_mesh(ista_fftph:iend_fftph)

! --
    integer :: i
    real(kind=DP) :: c_n, c_mx, c_my, c_mz
    real(kind=DP) :: ctarget_n, ctarget_mx, ctarget_my, ctarget_mz
    real(kind=DP) :: cos_th, sin_th, cos_ph, sin_ph
    real(kind=DP) :: d1, d2, c_tmp
    real(kind=DP) :: m_norm
!    real(kind=DP), parameter :: m_norm_criteria = 1.0D-20
    real(kind=DP), parameter :: m_norm_criteria = 1.0D-16

    Do i = ista_fftph, iend_fftph
       c_n  = RhoMag_R( i,1 );    c_mx = RhoMag_R( i,2 )
       c_my = RhoMag_R( i,3 );    c_mz = RhoMag_R( i,4 )
!
       m_norm = sqrt( c_mx**2 + c_my**2 + c_mz**2 )
!
       ctarget_n = bfft_kt(2*i,1)
       ctarget_mx = bfft_kt(2*i,2)
       ctarget_my = bfft_kt(2*i,3)
       ctarget_mz = bfft_kt(2*i,4)

       if ( m_norm < m_norm_criteria ) then
          d1 = ctarget_n /2.0d0;   d2 = d1
       else
          cos_th = c_mz / m_norm
          sin_th = sqrt( 1.0D0 - cos_th**2 )

          if ( sin_th > 1.0D-16 ) then
             cos_ph = ( c_mx /m_norm ) / sin_th
             sin_ph = ( c_my /m_norm ) / sin_th
          else
             cos_ph = 0.0d0;  sin_ph = 0.0d0
          endif
!
          c_tmp = ( ctarget_mx *cos_ph + ctarget_my *sin_ph )/2.0d0
          c_tmp = c_tmp *sin_th
          d1 = ( ctarget_n + ctarget_mz *cos_th ) /2.0d0 +c_tmp
          d2 = ( ctarget_n - ctarget_mz *cos_th ) /2.0d0 -c_tmp
       endif

       if ( quantz_axis_inversion_flg_mesh(i) == 1 ) then
          bfft_kt(2*i,1) = d1;   bfft_kt(2*i,2) = d2
       else
          bfft_kt(2*i,1) = d2;   bfft_kt(2*i,2) = d1
       endif
    End do

  end subroutine m_ES_GradCorr_Along_QuantzAxis

  subroutine m_ES_GradCorr_Along_QuantzAxis2( RhoMag_R, bfft_kt, &
       &                                    quantz_axis_inversion_flg_mesh )
    real(kind=DP), intent(in) :: RhoMag_R( ista_fftph:iend_fftph,ndim_magmom )
    real(kind=DP), intent(inout) :: bfft_kt( ista_fftp:iend_fftp,ndim_magmom )
    integer, intent(in) :: quantz_axis_inversion_flg_mesh(ista_fftph:iend_fftph)

! --
    integer :: i
    real(kind=DP) :: c_n, c_mx, c_my, c_mz
    real(kind=DP) :: ctarget_n, ctarget_mx, ctarget_my, ctarget_mz
    real(kind=DP) :: d1, d2, ctmp1, da1, da2
    real(kind=DP) :: m_norm, mmx, mmy, mmz

    real(kind=DP), parameter :: m_norm_criteria = 1.0D-20
    real(kind=DP), parameter :: m_norm_criteria2 = 1.0D-12

    Do i = ista_fftph, iend_fftph
       c_n  = RhoMag_R( i,1 );    c_mx = RhoMag_R( i,2 )
       c_my = RhoMag_R( i,3 );    c_mz = RhoMag_R( i,4 )
       m_norm = sqrt( c_mx**2 + c_my**2 + c_mz**2 )
!
       da1 = (c_n +m_norm) /2.0d0;  da2 = ( c_n - m_norm )/2.0d0
!
       if ( m_norm < m_norm_criteria ) then
          mmx = 0.0d0;  mmy = 0.0d0;  mmz = 0.0d0

       else if ( m_norm < m_norm_criteria2 ) then
          if ( quantz_axis_inversion_flg_mesh(i) == 1 ) then
             mmx = Global_Quantz_Axis_now(1)
             mmy = Global_Quantz_Axis_now(2)
             mmz = Global_Quantz_Axis_now(3)
          else
             mmx = -Global_Quantz_Axis_now(1)
             mmy = -Global_Quantz_Axis_now(2)
             mmz = -Global_Quantz_Axis_now(3)
          endif

       else
          mmx = c_mx / m_norm
          mmy = c_my / m_norm
          mmz = c_mz / m_norm
       endif
!
       ctarget_n = bfft_kt(2*i,1)
       ctarget_mx = bfft_kt(2*i,2)
       ctarget_my = bfft_kt(2*i,3)
       ctarget_mz = bfft_kt(2*i,4)
!
       ctmp1 = ctarget_mx *mmx + ctarget_my *mmy + ctarget_mz *mmz
!
       d1 = ( ctarget_n +ctmp1 ) /2.0d0
       d2 = ( ctarget_n -ctmp1 )/2.0d0

       if ( quantz_axis_inversion_flg_mesh(i) == 1 ) then
          bfft_kt(2*i,1) = d1;   bfft_kt(2*i,2) = d2
       else
          bfft_kt(2*i,1) = d2;   bfft_kt(2*i,2) = d1
       endif

       bfft_kt(2*i-1,1:4) = 0.0d0
       bfft_kt(2*i,3:4) = 0.0d0;

    End do

  end subroutine m_ES_GradCorr_Along_QuantzAxis2

! ===== KT_add === 2014/08/04
  subroutine m_ES_add_contrib_to_Vorticity( ixyz, bfft_kt, MagVorticity )
    integer, intent(in) :: ixyz
    real(kind=DP), intent(in)  :: bfft_kt( ista_fftp:iend_fftp,ndim_magmom )
    real(kind=DP), intent(out) :: MagVorticity( ista_fftph:iend_fftph, 3 )

    integer :: i
    real(kind=DP) :: ctarget_n, ctarget_mx, ctarget_my, ctarget_mz

    Do i = ista_fftph, iend_fftph
!       ctarget_n = bfft_kt(2*i,1)
       ctarget_mx = bfft_kt(2*i,2)
       ctarget_my = bfft_kt(2*i,3)
       ctarget_mz = bfft_kt(2*i,4)
!
       if ( ixyz == 1 ) then
          MagVorticity(i,2) = MagVorticity(i,2) -ctarget_mz
          MagVorticity(i,3) = MagVorticity(i,3) +ctarget_my
       else if ( ixyz == 2 ) then
          MagVorticity(i,3) = MagVorticity(i,3) -ctarget_mx
          MagVorticity(i,1) = MagVorticity(i,1) +ctarget_mz
       else if ( ixyz == 3 ) then
          MagVorticity(i,1) = MagVorticity(i,1) -ctarget_my
          MagVorticity(i,2) = MagVorticity(i,2) +ctarget_mx
       endif
    End Do
  end subroutine m_ES_add_contrib_to_Vorticity

  subroutine m_ES_neglect_low_Helicity( RhoMag_R, MagVorticity, Mat )
    real(kind=DP), intent(in) :: RhoMag_R( ista_fftph:iend_fftph, ndim_magmom )
    real(kind=DP), intent(in) :: MagVorticity( ista_fftph:iend_fftph, 3 )
    real(kind=DP), intent(inout) :: Mat( ista_fftph:iend_fftph )

    integer :: i, k
    real(kind=DP) :: ctmp1

    Do i=ista_fftph, iend_fftph
       ctmp1 = 0.0d0
       Do k=1, 3
          ctmp1 = ctmp1 +RhoMag_R(i,k+1)*MagVorticity(i,k)
       End do
       if ( abs(ctmp1) < threshold_helicity ) Mat(i) = 0.0d0
    End Do

  end subroutine m_ES_neglect_low_Helicity

  subroutine m_ES_neglect_low_Vorticity( RhoMag_R, MagVorticity, Mat )
    real(kind=DP), intent(in) :: RhoMag_R( ista_fftph:iend_fftph, ndim_magmom )
    real(kind=DP), intent(in) :: MagVorticity( ista_fftph:iend_fftph, 3 )
    real(kind=DP), intent(inout) :: Mat( ista_fftph:iend_fftph )

    integer :: i, k
    real(kind=DP) :: ctmp1

    Do i=ista_fftph, iend_fftph
       ctmp1 = 0.0d0
       Do k=1, 3
          ctmp1 = ctmp1 +MagVorticity(i,k)**2
       End do
       ctmp1 = sqrt(ctmp1)

       if ( ctmp1 < threshold_vorticity ) Mat(i) = 0.0d0
    End Do

  end subroutine m_ES_neglect_low_Vorticity
! ================ 2014/08/04

end module m_ES_NonCollinear

