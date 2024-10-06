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
!  This is a module for calculating Chi0 for the non-intearcting system.
! ======================================================================

! ======================= history ======================================
!  ver 1.0 :  2010/3/31
!               applicable to the crystal such as Si.
!  ver 2.0 :  2011/3/31
!               applicable to the isolated molecule such as C6H6.
!
! ======================================================================

module m_LinearResponse_NonInt

  use m_Kpoints,                    only : kv3_ek, kv3, vkxyz_ek, qwgt_ek, &
       &                                   k_symmetry, np0, np2
  use m_Const_Parameters,           only : DP, GAMMA, ELECTRON, DIRECT, &
       &                                   INVERSE, ON, OFF, CMPLDP, DELTA, &
       &                                   BUCS, &
       &                                   Hartree, NO, YES, &
       &                                   TETRAHEDRON, PARABOLIC

  use m_Control_Parameters,         only : ipri, nspin, kimg, af, neg, printable, &
       &                                   way_of_smearing
  use m_LinearResponse_Control,     only  : nrd_efermi, nmax_G_LR, &
       &                                    vqxyz,  &
       &                                    scissor, eta, sw_LongWaveLimit

  use m_LinearResponse_Density,     only : RhoTilde, occup_lkt_ek

  use m_Control_Parameters,         only : Num_q_Points
  use m_Electronic_Structure,       only : eko_ek
  use m_Crystal_Structure,    only :  rltv, univol
  use m_Timing,                     only : tstatc0_begin, tstatc0_end

  use m_Parallelization,          only : map_e,  map_z, myrank_e, MPI_CommGroup, &
       &                                 npes, ierr, mype

  use m_LinearResponse_Control,     only  : tddft_eqn_type, DYSON, BS

  Implicit None
  include 'mpif.h'

! ----------------------------- public -------------
  integer :: nband_LR, band_start_LR, band_end_LR
  Complex(kind=CMPLDP), allocatable :: Mat_Chi0(:,:,:)

contains

!------------------------------------------------------------------
!!
!!!            Alloc and Dealloc Chi0
!!
!-------------------------------------------------------------------
  subroutine m_LR_alloc_Mat_Chi0
    if ( tddft_eqn_type == DYSON ) then
       if ( mype == 0 ) then
          allocate( Mat_Chi0( nmax_G_LR, nmax_G_LR, nspin ) ); Mat_Chi0 = 0.0D0
       endif
    endif
    if ( tddft_eqn_type == BS ) then
       if ( mype == 0 ) then
          allocate( Mat_Chi0( 1,1,nspin ) ); Mat_Chi0 = 0.0D0
       endif
    endif
  end subroutine m_LR_alloc_Mat_Chi0

  subroutine m_LR_dealloc_Mat_Chi0
    if ( mype == 0 ) deallocate( Mat_Chi0 )
  end subroutine m_LR_dealloc_Mat_Chi0

!------------------------------------------------------------------
!!
!!!            Calc Chi0 in the case of Q > 0
!!
!-------------------------------------------------------------------
  subroutine Calc_Chi0_General( e_value )
    real(kind=DP), intent(in) :: e_value

    real(kind=DP)         :: weight
    real(kind=DP)         :: occ1, occ2
    complex(kind=CMPLDP)  :: ztmp, omega1, omega2
    integer               :: i,j, k, ispin
    integer               :: ik1, ik2, ik3
    integer               :: ib1, ib2

    integer               :: k_max, nq_max
    integer :: id_sname = -1
! ------------------------------ start ------------
    call tstatc0_begin('Calc_Chi0_General ', id_sname)
!
    nq_max = Num_q_points
    call set_value_kmax
    omega1 = cmplx( e_value, eta );  omega2 = cmplx( e_value, -eta );
! ---
    if ( mype==0 ) Mat_Chi0 = 0.0d0
! ------------------------------- main -----------
    Do i=1, nmax_G_LR
       Do j=1, nmax_G_LR

          Do ispin=1, nspin
             ztmp = 0.0d0

             Do k=1, k_max
                ik1 = ( nq_max +1 )*nspin*(k-1) + ispin
                ik2 = ik1 + nspin
                ik3 =  nspin*(k-1) + ispin

                weight = qwgt_ek(ik1) *nspin

                Do ib1=band_start_LR, band_end_LR           ! unoccpuied state
                   occ1 = Occup_lkt_ek( ib1, ik1 );
                   if ( occ1 > DELTA ) cycle

                   Do ib2=band_start_LR, band_end_LR          ! occupied state
                      occ2 = Occup_lkt_ek( ib2, ik2 )
                      if ( occ2 < DELTA ) cycle
!
                      if ( map_e(ib2) /= myrank_e ) cycle

                      call summation_ztmp
                   End do
                End do
             End do
             call add_ztmp_to_Chi0
          End do
       End do
    End do
! -------------------------------- end ----------
    if ( mype == 0 ) then
       if ( nspin==1 ) Mat_Chi0 = Mat_Chi0 *2.0d0
    endif
    call tstatc0_end(id_sname)

  contains

    subroutine set_value_kmax
      if ( way_of_smearing == TETRAHEDRON ) then
         k_max = np2
      else
         k_max = kv3_ek / ( nspin*( nq_max + 1 ))
      endif
    end subroutine set_value_kmax

    subroutine summation_ztmp
      integer :: jb2
      Real(kind=DP) :: ebi, ebj, ediff, c1
      Complex(kind=CMPLDP) :: zd1, zd2, z1, z2, gfn

      jb2 = map_z(ib2)

      ebi = Eko_ek( ib1,ik1 ); ebj = Eko_ek( ib2,ik2 )
      ediff = ebi - ebj + scissor

      zd1 = omega1 - ediff;  zd2 = omega2 + ediff
      gfn = 1.0d0 / zd1 - 1.0d0 / zd2
!                      c1 = occ1 - occ2
      z1 = RhoTilde( i, ib1, jb2, ik3 ); z2 = RhoTilde( j, ib1, jb2, ik3 )
!
      c1 = ( 1.0D0 - occ1 ) *occ2
      ztmp = ztmp + c1 *conjg( z1 )* z2 *gfn *weight
    end subroutine summation_ztmp

    subroutine add_ztmp_to_Chi0
      Complex(kind=CMPLDP) :: zsum

      if ( npes > 1 ) then
         call mpi_allreduce( ztmp, zsum, 2, MPI_DOUBLE_PRECISION, &
              &              MPI_SUM, MPI_CommGroup, ierr )
         ztmp = zsum
      endif
      if ( mype == 0 ) Mat_Chi0( i,j,ispin ) = ztmp / univol
    end subroutine add_ztmp_to_Chi0

  end subroutine Calc_Chi0_General

!------------------------------------------------------------------
!!
!!!            Calc Chi0 in the case of LongWaveLimit
!!
!-------------------------------------------------------------------
  subroutine Calc_Chi0_LWLimit( e_value )
    real(kind=DP), intent(in) :: e_value

    real(kind=DP)         :: occ1, occ2
    real(kind=DP)         :: weight
    complex(kind=CMPLDP)  :: ztmp, omega1, omega2

    integer               :: i,j, k, ispin
    integer               :: ik1, ik2, ik3
    integer               :: ib1, ib2
    integer               :: k_max, nq_max
    integer :: id_sname = -1
! ------------------------------- start ------------
    call tstatc0_begin('Calc_Mat_Chi0_LWLimit ', id_sname)
!
    nq_max = 0
    call set_value_kmax
    omega1 = cmplx( e_value, eta );  omega2 = cmplx( e_value, -eta );
! --
    if ( mype==0 ) Mat_Chi0 = 0.0d0
! -------------------------------- main -------------
    Do i=1, nmax_G_LR
!       write(*,*) ' i = ', i
       Do j=1, nmax_G_LR

          Do ispin=1, nspin
             ztmp = 0.0d0

             Do k=1, k_max
                ik1 = ( nq_max +1 )*nspin*(k-1) + ispin
                ik2 = ik1 + nspin
                ik3 =  nspin*(k-1) + ispin
                !
                weight = qwgt_ek(ik1) *nspin

                Do ib1=band_start_LR, band_end_LR           ! unoccpuied state
                   occ1 = Occup_lkt_ek( ib1, ik1 );
                   if ( occ1 > DELTA ) cycle

                   Do ib2=band_start_LR, band_end_LR          ! occupied state
                      occ2 = Occup_lkt_ek( ib2, ik1 )
                      if ( occ2 < DELTA ) cycle

                      if ( map_e(ib2) /= myrank_e ) cycle

                      call summation_ztmp
                   End do
                End do
             End do
             call add_ztmp_to_Chi0
          End do
       End do
    End do
! --------------------------------- end -------------
    if ( mype == 0 ) then
       if ( nspin==1 ) Mat_Chi0 = Mat_Chi0 *2.0d0
    endif
    call tstatc0_end(id_sname)

  contains

    subroutine set_value_kmax
      if ( way_of_smearing == TETRAHEDRON ) then
         k_max = np2
      else
         k_max = kv3_ek / ( nspin*( nq_max + 1 ))
      endif
    end subroutine set_value_kmax

    subroutine summation_ztmp
      integer :: jb2
      Real(kind=DP) :: ebi, ebj, ediff, c1
      Complex(kind=CMPLDP) :: zd1, zd2, z1, z2, gfn

      jb2 = map_z(ib2)

      ebi = Eko_ek( ib1,ik1 ); ebj = Eko_ek( ib2,ik1 )
      ediff = ebi - ebj + scissor

      zd1 = omega1 - ediff;  zd2 = omega2 + ediff
      gfn = 1.0d0 / zd1 - 1.0d0 / zd2
!
      z1 = RhoTilde( i, ib1,jb2, ik3 );  z2 = RhoTilde( j, ib1,jb2, ik3 )
!                      c1 = occ1 - occ2
      c1 = ( 1.0D0 - occ1 ) *occ2
      ztmp = ztmp + c1 *conjg( z1 )* z2 *gfn *weight
    end subroutine summation_ztmp

    subroutine add_ztmp_to_Chi0
      Complex(kind=CMPLDP) :: zsum

      if ( npes > 1 ) then
         call mpi_allreduce( ztmp, zsum, 2, MPI_DOUBLE_PRECISION, &
              &              MPI_SUM, MPI_CommGroup, ierr )
         ztmp = zsum
      endif
      if ( mype == 0 ) Mat_Chi0( i,j,ispin ) = ztmp / univol
    end subroutine add_ztmp_to_Chi0

  end subroutine Calc_Chi0_LWLimit

!------------------------------------------------------------------
!!
!!!            Calc Total Charge for DEBUG
!!
!-------------------------------------------------------------------
  subroutine Check_Total_Charge
    real(kind=DP)         :: occ1
    complex(kind=CMPLDP)  :: ztmp, z1

    integer               :: i,j, k, ispin
    integer               :: ik1, ik2, ik3
    integer               :: ib1, jb1
    integer               :: k_max, nq_max

    integer :: id_sname = -1
! ----------------------- start ----------------
    call set_value_nqmax
    call set_value_kmax
! ------------------------------ main ----------
    Do i=1, nmax_G_LR
       ztmp = 0.0d0

       Do ispin=1, nspin
          Do k=1, k_max
             ik1 = ( nq_max +1 )*nspin*(k-1) + ispin
             ik2 = ik1 + nspin
             ik3 =  nspin*(k-1) + ispin
             !
             Do ib1=1, neg           ! occpuied state
                occ1 = Occup_lkt_ek( ib1, ik1 );
                if ( occ1 < 0.5 ) cycle

                if ( map_e(ib1) /= myrank_e ) cycle
                jb1 = map_z(ib1)

                z1 = RhoTilde( i, ib1, jb1, ik3 )
                ztmp = ztmp + occ1 * z1 *qwgt_ek(ik1)
             End do
          End do
       End do
       if ( i==1 ) call print_total_charge
    End do

  contains

    subroutine set_value_nqmax
      if ( sw_LongWaveLimit==ON ) then
         nq_max = 0
      else
         nq_max = Num_Q_Points
      endif
    end subroutine set_value_nqmax

    subroutine set_value_kmax
      if ( way_of_smearing == TETRAHEDRON ) then
         k_max = np2
      else
         k_max = kv3_ek / ( nspin*( nq_max + 1 ))
      endif
    end subroutine set_value_kmax

    subroutine print_total_charge
      Complex(kind=CMPLDP) :: zsum

      if ( npes > 1 ) then
         call mpi_allreduce( ztmp, zsum, 1, MPI_DOUBLE_PRECISION, &
              &              MPI_SUM, MPI_CommGroup, ierr )
         ztmp = zsum
      endif
      write(*,*) 'Charge i ', i, ztmp / univol, ztmp
    end subroutine print_total_charge

  end subroutine Check_Total_Charge

end module m_LinearResponse_NonInt
