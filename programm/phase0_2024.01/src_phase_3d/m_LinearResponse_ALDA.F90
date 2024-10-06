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
!  This is a module for calculating exchange kernels on the real grid
!  based on ALDA.
! ======================================================================

! ======================= history ======================================
!  ver 1.0 :  2010/3/31
!               applicable to the crystal such as Si.
!  ver 2.0 :  2011/3/31
!               applicable to the isolated molecule such as C6H6.
!
! ======================================================================

module m_LinearResponse_ALDA

  use m_Const_Parameters,   only : ON, OFF, DP, CMPLDP, DELTA, &
       &                           Partial_Core_Charge, &
       &                           Valence_plus_PC_Charge, &
       &                           TETRAHEDRON, PARABOLIC, &
       &                           GGA, LDA, PAI2, PAI

  use m_Control_Parameters,  only : Num_q_Points, way_of_smearing, neg, &
       &                            ipri, nspin, printable, xctype, kimg

  use m_Kpoints,                    only : kv3_ek, kv3, vkxyz_ek, qwgt_ek, &
       &                                   k_symmetry, np0, np2
  use m_LinearResponse_Control,     only  : nrd_efermi, nstep, e, nmax_G_LR, &
       &                                    vqxyz,  &
       &                                    xc_kernel_type, RPA, ALDA_G, LRC,  &
       &                                    LRC_alpha, sw_NLF, &
       &                                    spectrum_type, OPTICS, EELS, &
       &                                    scissor, eta, &
       &                                    sw_LongWaveLimit

  use m_LinearResponse_tools,  only  : occup_k, occup_kmq, &
       &                                  eko_k, &
       &                                    fsval_k, fsval_kmq, &
       &                                    wfn_k, wfn_kmq, &
       &                                    map_z_k, map_z_kmq, &
       &                                   nbase_k, nbase_kmq, &
       &                                   igf_k, igf_kmq, iba_k, iba_kmq

  use m_LinearResponse_Density,  only : RhoTilde, occup_lkt_ek

  use m_Electronic_Structure,       only : eko_ek
  use m_Crystal_Structure,    only :  rltv, univol

  use m_FFT,                  only : fft_box_size_CD, fft_box_size_CD_c, &
       &                             nfft, nfftp, &
       &                             m_FFT_CD_inverse_c, &
       &                             m_FFT_check_of_negative_CD, &
       &                             m_FFT_CD_direct_c, &
       &                             m_FFT_alloc_CD_box, &
       &                             m_FFT_dealloc_CD_box


#ifdef _MPIFFT_
  use m_FFT,                  only : m_FFT_set_cdata  &
       &                           , lx,ly,lz,ly_d,lz_d, ny_d,nz_d
  use m_PlaneWaveBasisSet,    only : igfp_l_c
#else
  use m_FFT,                  only : m_FFT_alloc_CD_box &
       &                           , m_FFT_dealloc_CD_box
#endif

  use m_Files,                 only : nfout
  use m_PlaneWaveBasisSet,    only : ngabc, igf, kg1,kg, igfp_l, kgp

  use m_Parallelization,      only : ista_kngp,iend_kngp,npes,mype, ierr, &
       & npes_cdfft, nrank_ggacmp, myrank_cdfft, &
       & myrank_ggacmp, ista_fftp, iend_fftp, ista_fftph, iend_fftp, &
       & nis_fftp, nie_fftp, nel_fftp, idisp_fftp,np_fftp,mp_fftp,  &
       & nel_fftph, idisp_fftph, nrest_cdfft, &
       & mpi_ggacmp_cross_world,mpi_cdfft_world, &
       & map_ggacmp, map_pe2ggacmp, map_pe2cdfft, MPI_CommGroup, &
       & ista_fftph, iend_fftph, &
       & np_e, ista_k, iend_k, map_e, map_z, myrank_e

  use m_Charge_Density,        only : chgq_l
  use m_XC_Potential,         only : check_of_xctype
  use m_Ionic_System,         only : ntyp,natm,iwei,ityp,pos,zfm3_l
  use m_PseudoPotential,      only : itpcc,ilmt,ltp,taup,il2p,isph,iqitg, &
       & rhpcg_l, qitg_l,rhpcg_diff_l,qitg_diff_l,dl2p, &
       & m_PP_include_vanderbilt_pot, m_PP_find_maximum_l,nlmt

  use m_Timing,               only : tstatc0_begin, tstatc0_end


  use m_PlaneWaveBasisSet,    only : kg0, nbase
  use m_FFT,                  only : fft_box_size_WF, &
       &                             m_FFT_alloc_wf_work, &
       &                             m_FFT_dealloc_wf_work, m_FFT_wf
  use m_Const_Parameters,     only : DIRECT, ELECTRON

  use m_LinearResponse_tools,  only  : Get_WF_in_Rspace, Get_WF_in_Rspace_CD, &
       &                               fsval_k, fsval_kmq

  use m_Const_Parameters,      only : PAI4, BUCS, Bohr
  use m_IterationNumbers,     only : nk_in_the_process

  use m_Electronic_Structure,  only : zaj_l
  use m_Const_Parameters ,    only : INVERSE, DIRECT, OFF, GAMMA

  use m_LinearResponse_NonInt,   only : nband_LR, band_start_LR, band_end_LR
!
  use m_LinearResponse_Control,  only : sw_Fxc_enhanced
!
  implicit none
  include 'mpif.h'

! --------------------------- public -------------
  integer, allocatable :: igfp_LR(:)                ! d(kgp)
! -------------
  real(kind=DP), allocatable :: ker_ex(:,:), ker_cr(:,:,:)
  real(kind=DP), allocatable :: ker_xc_triplet(:,:,:),  ker_xc_singlet(:,:,:)
  real(kind=DP), allocatable :: zcos_LR(:), zsin_LR(:)
! ----------------------------------------- DYSON && ALDA_R ---
  Complex(kind=CMPLDP), allocatable :: MatrixT(:,:,:,:)
  Complex(kind=CMPLDP), allocatable :: MatrixG(:,:,:,:,:)
!                                  Matrix G =  RhoTilde * Fxc
! ---------------------------- temporary ----------
  real(kind=DP), allocatable,dimension(:):: afft     ! d(ista_fftp:iend_fftp)
!
  real(kind=DP), private, allocatable, dimension(:,:)  :: chden_l
  real(kind=DP), private, allocatable, dimension(:,:)  :: grad_rho
  real(kind=DP), private, allocatable, dimension(:)  :: grad_trho
#ifndef _XC_SAVE_MEMORY_
  real(kind=DP),private,allocatable,dimension(:,:,:):: cgrad_rho ! MPI d(ista_fftph:iend_fftph,3,nspin)
#else
  real(kind=DP),private,allocatable,dimension(:):: cggawk13 ! MPI d(ista_fftph:iend_fftph)
#endif
! -----------------------
  real(kind=DP),private,pointer,dimension(:,:)  :: chgrhr_l    ! MPI d(ista_fftph:iend_fftph,ispin)
! GGA arrays
  integer, private,allocatable,dimension(:)   :: inx,jnx,knx ! MPI d(ista_fftp:iend_fftp
  real(kind=DP),private, pointer, dimension(:)   :: f2or1 ! MPI d(ista_fftph,iend_fftph)
  logical,      private, dimension(3)               :: lmn_even

  integer istatus( mpi_status_size )

! ------------
contains

!------------------------------------------------------------------
!!
!!!            Alloc and Dealloc Kernels evaluated on R-space
!!
!-------------------------------------------------------------------
  subroutine m_LR_alloc_RKernels_ALDA
    Allocate( ker_ex( ista_fftph:iend_fftph,nspin ) ); ker_ex = 0.0d0
    Allocate( ker_cr( ista_fftph:iend_fftph,nspin,nspin ) ); ker_cr = 0.0d0
  end subroutine m_LR_alloc_RKernels_ALDA

  subroutine m_LR_dealloc_RKernels_ALDA
    deallocate( ker_ex, ker_cr )
  end subroutine m_LR_dealloc_RKernels_ALDA

  subroutine m_LR_alloc_phase_factor
    Allocate( zcos_LR(ista_fftph:iend_fftph) );  zcos_LR = 0.0d0
    Allocate( zsin_LR(ista_fftph:iend_fftph) );  zsin_LR = 0.0d0
  end subroutine m_LR_alloc_phase_factor

  subroutine m_LR_dealloc_phase_factor
    Deallocate( zcos_LR, zsin_LR )
  end subroutine m_LR_dealloc_phase_factor

!------------------------------------------------------------------
!!
!!!            Alloc and Dealloc in the case of DYSON && ALDA_R
!!
!-------------------------------------------------------------------
  subroutine m_LR_alloc_Matrix_G
    if ( sw_LongWaveLimit == OFF ) then
       Allocate( MatrixG( kg1, neg, np_e, kv3_ek/(Num_q_Points+1), nspin ) )
    else
       Allocate( MatrixG( kg1, neg, np_e, kv3_ek, nspin ) )
    endif
  end subroutine m_LR_alloc_Matrix_G

  subroutine m_LR_dealloc_Matrix_G
    Deallocate( MatrixG )
  end subroutine m_LR_dealloc_Matrix_G

  subroutine m_LR_alloc_Matrix_T
    if ( mype == 0 ) then
       allocate( MatrixT( nmax_G_LR, nmax_G_LR, nspin, nspin ) )
       MatrixT = 0.0d0
    endif
  end subroutine m_LR_alloc_Matrix_T

  subroutine m_LR_dealloc_Matrix_T
    if ( mype == 0 ) deallocate( MatrixT )
  end subroutine m_LR_dealloc_Matrix_T

!-------------------------------------------------------
!!
!!!       igfp for LR
!!
!--------------------------------------------------------
  subroutine m_LR_alloc_igfp
    allocate( igfp_LR(kgp) ); igfp_LR = 0
  end subroutine m_LR_alloc_igfp

  subroutine m_LR_dealloc_igfp
    deallocate( igfp_LR )
  end subroutine m_LR_dealloc_igfp

  subroutine m_LR_set_igfp
    integer, allocatable     :: igfp_tmp(:)

    if ( npes > 1 ) then
       allocate( igfp_tmp(kgp) ); igfp_tmp = 0
       igfp_tmp( ista_kngp:iend_kngp ) = igfp_l( ista_kngp:iend_kngp )
       call mpi_allreduce( igfp_tmp, igfp_LR, kgp, mpi_integer, mpi_sum, &
            &              MPI_CommGroup, ierr )
       deallocate( igfp_tmp )
    else
       igfp_LR( ista_kngp:iend_kngp ) = igfp_l( ista_kngp:iend_kngp )
    endif
  end subroutine m_LR_set_igfp

!--------------------------------------------------------
!!
!!!        Matrix T  ( DYSON && ALDA_R ) in the case of Q > 0
!!
!---------------------------------------------------------
  subroutine Calc_MatrixT_General( e_value )
    real(kind=DP), intent(in) :: e_value

    integer               :: i,j, k, ispin, jspin
    integer               :: ik1, ik2, ik3, ib1, ib2
    integer               :: k_max, nq_max

    real(kind=DP)         :: occ1, occ2, weight
    complex(kind=CMPLDP)  :: omega1, omega2, ztmp
!
    integer :: id_sname = -1
! -------------------------------- start -----------
    call tstatc0_begin('Calc_MatrixT ', id_sname)

    nq_max = Num_q_points
    call set_value_kmax
    omega1 = cmplx( e_value, eta );  omega2 = cmplx( e_value, -eta );
! --
    if ( mype == 0) MatrixT = 0.0d0
! --------------------------------- main ----------
    Do i=1, nmax_G_LR
       Do j=1, nmax_G_LR

          Do ispin=1, nspin
             Do jspin=1, nspin
                ztmp = 0.0d0

                Do k=1, k_max
                   ik1 = ( nq_max +1 )*nspin*(k-1) + ispin
                   ik2 = ik1 + nspin
                   ik3 =  nspin*(k-1) + ispin

                   weight = qwgt_ek(ik1) *nspin

                   Do ib1=band_start_LR, band_end_LR      ! unoccpuied state
                      occ1 = Occup_lkt_ek( ib1, ik1 );
                      if ( occ1 > DELTA ) cycle

                      Do ib2=band_start_LR, band_end_LR     ! occupied state
                         occ2 = Occup_lkt_ek( ib2, ik2 )
                         if ( occ2 < DELTA ) cycle
                         if ( map_e(ib2) /= myrank_e ) cycle

                         call summation_ztmp
                      End do
                   End do
                End do
                call add_ztmp_to_MatrixT
             End do
          End do
       End do
    End do
!! ----------------------------- end ------------
    if ( mype == 0 ) then
       if ( nspin==1 ) MatrixT = MatrixT *2.0d0
!!    if ( nspin==1 ) MatrixT = MatrixT *2.0d0 * 2.0d0
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

      z1 = RhoTilde( i, ib1, jb2, ik3 )
      z2 = MatrixG( j, ib1, jb2,  ik3, jspin )
!
      c1 = ( 1.0D0 - occ1 ) *occ2
!      ztmp = ztmp + c1 *conjg( z1 )* z2 *gfn *weight
      ztmp = ztmp + c1 *z1 * z2 *gfn *weight
    end subroutine summation_ztmp

    subroutine add_ztmp_to_MatrixT
      Complex(kind=CMPLDP) :: zsum

      if ( npes > 1 ) then
         call mpi_allreduce( ztmp, zsum, 2, MPI_DOUBLE_PRECISION, &
              &              MPI_SUM, MPI_CommGroup, ierr )
         ztmp = zsum
      endif
      if ( mype == 0 ) MatrixT( i,j,ispin,jspin ) = ztmp / univol
    end subroutine add_ztmp_to_MatrixT

  end subroutine Calc_MatrixT_General

!--------------------------------------------------------
!!
!!!        Matrix T  ( DYSON && ALDA_R ) in the case of LongWaveLimit
!!
!---------------------------------------------------------
  subroutine Calc_MatrixT_LWLimit( e_value )
    real(kind=DP), intent(in) :: e_value

    integer               :: i,j, k, ispin, jspin
    integer               :: ik1, ik2, ik3, ib1, ib2
    integer               :: k_max, nq_max
    real(kind=DP)         :: occ1, occ2, weight
    complex(kind=CMPLDP)  :: omega1, omega2, ztmp
! ----
    integer :: id_sname = -1
! ------------------------------- start ---------
    call tstatc0_begin('Calc_MatrixT_LWLimit ', id_sname)
!
    nq_max = 0
    call set_value_kmax
    omega1 = cmplx( e_value, eta );  omega2 = cmplx( e_value, -eta );
! --
    if ( mype == 0 ) MatrixT = 0.0d0
! ------------------------------ main -----------------
    Do i=1, nmax_G_LR
       Do j=1, nmax_G_LR

          Do ispin=1, nspin
             Do jspin=1, nspin
                ztmp = 0.0d0

                Do k=1, k_max
                   ik1 = ( nq_max +1 )*nspin*(k-1) + ispin
                   ik2 = ik1 + nspin
                   ik3 =  nspin*(k-1) + ispin
                !
                   weight = qwgt_ek(ik1) *nspin

                   Do ib1=band_start_LR, band_end_LR         ! unoccpuied state
                      occ1 = Occup_lkt_ek( ib1, ik1 );
                      if ( occ1 > DELTA ) cycle

                      Do ib2=band_start_LR, band_end_LR         ! occupied state
                         occ2 = Occup_lkt_ek( ib2, ik1 )
                         if ( occ2 < DELTA ) cycle
                         if ( map_e(ib2) /= myrank_e ) cycle

                         call summation_ztmp
                      End do
                   End do
                End do
                call add_ztmp_to_MatrixT
             End do
          End do
       End do
    End do
!! ----------------------------- end ------------
    if ( mype == 0 ) then
       if ( nspin==1 ) MatrixT = MatrixT *2.0d0
!!    if ( nspin==1 ) MatrixT = MatrixT *2.0d0 * 2.0d0
    endif
!    MatrixT = 0.0d0
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
      integer       :: jb2
      Real(kind=DP) :: ebi, ebj, ediff, c1
      Complex(kind=CMPLDP) :: zd1, zd2, z1, z2, gfn

      jb2 = map_z(ib2)

      ebi = Eko_ek( ib1,ik1 ); ebj = Eko_ek( ib2,ik1 )
      ediff = ebi - ebj + scissor

      zd1 = omega1 - ediff;  zd2 = omega2 + ediff
      gfn = 1.0d0 / zd1 - 1.0d0 / zd2

      z1 = RhoTilde( i, ib1, jb2, ik3 )
      z2 = MatrixG(  j, ib1, jb2, ik3, jspin )

      c1 = ( 1.0D0 - occ1 ) *occ2
      ztmp = ztmp + c1 *conjg( z1 )* z2 *gfn *weight
!!!!!!      ztmp = ztmp + c1 * z1 * z2 *gfn *weight
    end subroutine summation_ztmp

    subroutine add_ztmp_to_MatrixT
      Complex(kind=CMPLDP) :: zsum

      if ( npes > 1 ) then
         call mpi_allreduce( ztmp, zsum, 2, MPI_DOUBLE_PRECISION, &
              &              MPI_SUM, MPI_CommGroup, ierr )
         ztmp = zsum
      endif
      if ( mype == 0 ) MatrixT( i,j,ispin,jspin ) = ztmp / univol
    end subroutine add_ztmp_to_MatrixT

  end subroutine Calc_MatrixT_LWLimit

!--------------------------------------------------------
!!
!!!        Matrix G   ( DYSON && ALDA_R ) in the case of Q > 0
!!
!---------------------------------------------------------
  subroutine Add_SoftPart_MatG_General
    Real(kind=DP), allocatable :: wf1(:), wf2(:), chg1r(:), chg1i(:)
    Real(kind=DP), allocatable :: wf2_mpi(:), wf_tmp(:)

    Real(kind=DP) :: occ1, occ2
    integer ngrid, nffth, nfftph
    integer ik, jk
    integer ib1, ib2
    integer :: id_sname = -1
! --------------------------- start ----------
    call tstatc0_begin('Add_SoftPart_MatG_General', id_sname)
!
    ngrid = product(fft_box_size_CD(1:3,1))
    nffth  = nfft/2;    nfftph = nfftp/2
    call m_FFT_alloc_CD_box()
    call set_work_arrays()
! ---------------------------- main ------------------
    Do ik=1, kv3

       Do ib1=1, neg                  ! unOccupied Band

          if ( nrd_efermi ==1 ) then
             occ1 = occup_k( ib1,ik )
             if ( occ1 > DELTA ) cycle
          endif
          call set_wf1()

          Do ib2=1, neg                ! occupied Band

             if ( nrd_efermi==1 ) then
                occ2 = occup_kmq( ib2,ik )
                if ( occ2 < DELTA ) cycle
             end if
             call set_wf2

             Do jk=1, kv3              ! kv3 should be equal to nspin
                call set_chg1
                call set_each_element_in_MatG
             End do
          End do
       End do
    End do
! -----------------------------------end --------------------
    call unset_work_arrays
    call m_FFT_dealloc_CD_box()
    call tstatc0_end(id_sname)

  contains

    subroutine set_work_arrays
      allocate( wf1(ista_fftp:iend_fftp) ); wf1 = 0.0d0
      allocate( wf2(ista_fftp:iend_fftp) ); wf2 = 0.0d0
      allocate( chg1r(ista_fftph:iend_fftph) ); chg1r = 0.0d0
      allocate( chg1i(ista_fftph:iend_fftph) ); chg1i = 0.0d0
      allocate( wf2_mpi(nfftp) ); wf2_mpi= 0.0d0
      if ( npes > 1 ) then
         allocate( wf_tmp(nfftp) ); wf_tmp = 0.0d0
      endif
    end subroutine set_work_arrays

    subroutine unset_work_arrays
      deallocate( wf1, wf2 ); deallocate( chg1r, chg1i )
      deallocate( wf2_mpi )
      if ( npes > 1 ) deallocate( wf_tmp )
    end subroutine unset_work_arrays

    subroutine set_wf1()
      call  Get_WF_in_Rspace_CD( ik, ib1, wfn_k, wf1, &
           &                     nbase_k, iba_k, igfp_LR, map_z_k )
      wf1 = wf1 / sqrt(univol)
    end subroutine set_wf1

    subroutine set_wf2()
      call  Get_WF_in_Rspace_CD( ik, ib2, wfn_kmq, wf2, &
           &                     nbase_kmq, iba_kmq, igfp_LR, map_z_kmq )
      wf2 = wf2 / sqrt(univol)
    end subroutine set_wf2

    subroutine set_chg1()
      integer :: i, ip

      Do i=ista_fftph, iend_fftph
         ip = 2*i-1
         chg1r(i) = wf1(ip) *wf2(ip)   + wf1(ip+1) *wf2(ip+1)
         chg1i(i) = wf1(ip) *wf2(ip+1) - wf1(ip+1) *wf2(ip)
      End do
    end subroutine set_chg1

    subroutine set_each_element_in_MatG
      integer :: i, ip, ri
      integer :: ic1, jb2
      Real(kind=DP) :: weight

      Do i=ista_fftph,iend_fftph
         ip = 2*i-1
         weight = ker_cr( i,ik,jk )
         chg1r(i) =  weight *chg1r(i);   chg1i(i) =  weight *chg1i(i)
      End do
      wf2 = 0.0d0
      do i=ista_fftph,iend_fftph
         ip = 2*i-1
         wf2(ip)   =  chg1r(i);  wf2(ip+1) =  chg1i(i)
      end do
      call m_FFT_CD_direct_c( nfout, wf2 )   !  Rspace -> G-space
!--
      wf2_mpi = 0.0d0
      wf2_mpi(ista_fftp:iend_fftp) = wf2(ista_fftp:iend_fftp) / dble(ngrid) *univol
      if ( npes > 1 ) then
         call mpi_allreduce( wf2_mpi, wf_tmp, nfftp, MPI_DOUBLE_PRECISION, &
              &              MPI_SUM, MPI_CommGroup, ierr )
         wf2_mpi = wf_tmp
      endif
! --------------------
      if ( map_e(ib2) == myrank_e ) then
         jb2 = map_z_k(ib2)

         ic1 = ( nk_in_the_process -1 )/ ( nspin*(Num_q_Points+1) )
         Do i=1, iba_k(ik)
            ip = 2 *igfp_LR( nbase_k(i,ik) ) -1
            MatrixG( i, ib1,jb2, ik +ic1*nspin, jk ) &
                 &     = dcmplx( wf2_mpi(ip), wf2_mpi(ip+1) )
!               MatrixG( i, ib1,jb2, ik +ic1*nspin, jk ) &
!                    &     = dcmplx( wf2_mpi(ip), -wf2_mpi(ip+1) )
         Enddo
      endif
    end subroutine set_each_element_in_MatG

  end subroutine Add_SoftPart_MatG_General

!--------------------------------------------------------
!!
!!!        Matrix G   ( DYSON && ALDA_R ) in the case of LongWaveLimit
!!
!---------------------------------------------------------
  subroutine Add_SoftPart_MatG_LWLimit
    Real(kind=DP), allocatable :: wf1(:), wf2(:), chg1r(:), chg1i(:)
    Real(kind=DP), allocatable :: wf2_mpi(:), wf_tmp(:)

    Real(kind=DP) ::  occ1, occ2
    integer ngrid, nffth, nfftph
    integer ik, jk
    integer ib1, ib2
    integer :: id_sname = -1
!
    integer iii
! ------------------------------ start -------------
    call tstatc0_begin('Add_SoftPart_MatG_LWLimit ', id_sname)
! -
    ngrid = product(fft_box_size_CD(1:3,1))
    nffth  = nfft/2;    nfftph = nfftp/2
    call m_FFT_alloc_CD_box()
    call set_work_arrays()
! ------------------------------- main --------------
    Do ik=1, kv3

       Do ib1=1, neg               ! unOccupied Band

          if ( nrd_efermi ==1 ) then
             occ1 = occup_k( ib1,ik )
             if ( occ1 > DELTA ) cycle
          endif
          call set_wf1()
!            Do iii=1, nfftp
!               write(*,*) ' wf1 i ', iii, wf1(iii)
!            End Do
!            stop

          Do ib2=1, neg                ! occupied Band

             if ( nrd_efermi==1 ) then
                occ2 = occup_k( ib2,ik )
                if ( occ2 < DELTA ) cycle
             end if
             call set_wf2

             Do jk=1, kv3              ! kv3 should be equal to nspin
                call set_chg1
                call set_each_element_in_MatG
             End do
          End do
       End do
    End do
! ------------------------------end ---------------
    call unset_work_arrays
    call m_FFT_dealloc_CD_box()
    call tstatc0_end(id_sname)

      Do jk=1, nmax_G_LR
         Do ib1=1, neg
            Do ib2=1, neg
               write(*,*) 'AA A ', jk, ib1, ib2, MatrixG( jk,ib1,ib2,1,1)
            End Do
         End Do
      End do
      stop
! ------------

  contains

    subroutine set_work_arrays()
      allocate( wf1(ista_fftp:iend_fftp) ); wf1 = 0.0d0
      allocate( wf2(ista_fftp:iend_fftp) ); wf2 = 0.0d0
      allocate( chg1r(ista_fftph:iend_fftph) ); chg1r = 0.0d0
      allocate( chg1i(ista_fftph:iend_fftph) ); chg1i = 0.0d0
      allocate( wf2_mpi(nfftp) ); wf2_mpi= 0.0d0
      if ( npes > 1 ) then
         allocate( wf_tmp(nfftp) ); wf_tmp = 0.0d0
      endif
    end subroutine set_work_arrays

    subroutine unset_work_arrays()
      deallocate( wf1, wf2 ); deallocate( chg1r, chg1i )
      deallocate( wf2_mpi )
      if ( npes > 1 ) deallocate( wf_tmp )
    end subroutine unset_work_arrays

    subroutine set_wf1
      call  Get_WF_in_Rspace_CD( ik, ib1, wfn_k, wf1, &
           &                     nbase_k, iba_k, igfp_LR, map_z_k )
    end subroutine set_wf1

    subroutine set_wf2
      call  Get_WF_in_Rspace_CD( ik, ib2, wfn_k, wf2, &
           &                     nbase_k, iba_k, igfp_LR,  map_z_k )
      wf2 = wf2 / sqrt(univol)
    end subroutine set_wf2

    subroutine set_chg1
      integer :: i, ip

      Do i=ista_fftph, iend_fftph
         ip = 2*i-1
         chg1r(i) = wf1(ip) *wf2(ip)   + wf1(ip+1) *wf2(ip+1)
         chg1i(i) = wf1(ip) *wf2(ip+1) - wf1(ip+1) *wf2(ip)
      End do
    end subroutine set_chg1

    subroutine set_each_element_in_MatG
      integer :: i, ip, ri
      integer :: ic1, jb2
      Real(kind=DP) :: weight

      Do i=ista_fftph,iend_fftph
         ip = 2*i-1
         weight = ker_cr( i,ik,jk )
         chg1r(i) =  weight *chg1r(i);   chg1i(i) =  weight *chg1i(i)
      End do
      wf2 = 0.0d0
      do i=ista_fftph,iend_fftph
         ip = 2*i-1
         wf2(ip)   =  chg1r(i);  wf2(ip+1) =  chg1i(i)
      end do
      call m_FFT_CD_direct_c( nfout, wf2 )   !  Rspace -> G-space
!--
      wf2_mpi = 0.0d0
      wf2_mpi(ista_fftp:iend_fftp) = wf2(ista_fftp:iend_fftp) / dble(ngrid) *univol
      if ( npes > 1 ) then
         call mpi_allreduce( wf2_mpi, wf_tmp, nfftp, MPI_DOUBLE_PRECISION, &
              &              MPI_SUM, MPI_CommGroup, ierr )
         wf2_mpi = wf_tmp
      endif
! --------------------
      if ( map_e(ib2) == myrank_e ) then
         jb2 = map_z_k(ib2)

         ic1 = ( nk_in_the_process -1 )/ nspin
         Do i=1, iba_k(ik)
            ip = 2 *igfp_LR( nbase_k(i,ik) ) -1
            MatrixG( i, ib1,jb2, ik +ic1*nspin, jk ) &
                 &     = dcmplx( wf2_mpi(ip), wf2_mpi(ip+1) )
!            MatrixG( i, ib1,jb2, ik +ic1*nspin, jk ) &
!                 &     = dcmplx( wf2_mpi(ip), -wf2_mpi(ip+1) )
         Enddo
      endif
    end subroutine set_each_element_in_MatG

  end subroutine Add_SoftPart_MatG_LWLimit

!-----------------------------------------------------
!!
!!!          phase factor
!!
!------------------------------------------------------
  subroutine calc_phase_factor
    integer i, ip, ngrid
    Real(kind=DP) :: ph, da(3), tmp_vec(3)
!
    ngrid = product( fft_box_size_CD(1:3,1) )
    da(1:3) = 1.d0 / fft_box_size_CD(1:3,1)
!
    write(*,*) 'UUUU', da(1), da(2), da(3)
!    stop
    Do i = ista_fftph, iend_fftph
       ip = 2*i - 1
        write(*,*) 'ZZZ', i
       tmp_vec(1) = dble( inx(ip) ) *da(1)
       tmp_vec(2) = dble( jnx(ip) ) *da(2)
       tmp_vec(3) = dble( knx(ip) ) *da(3)
!         write(*,*) 'YY', ip, inx(ip), jnx(ip), knx(ip)
!         write(*,*) 'YY', ip, tmp_vec(1), tmp_vec(2), tmp_vec(3)

       ph = -PAI2 * dot_product( vqxyz(1,:,BUCS), tmp_vec(:) )
       zcos_LR(i) = cos(ph);      zsin_LR(i) = sin(ph)
    End do
  end subroutine calc_phase_factor

  subroutine set_ispin_kt(ispin, input_Charge )
    integer, intent(out) :: ispin
    integer, intent(in)  :: input_Charge
    if( nspin == 2 .and. input_charge == Valence_plus_PC_Charge ) then
       ispin = 2
    else
       ispin = 1
    endif
  end subroutine set_ispin_kt

  subroutine Calc_Density_in_Rspace(ispin, input_Charge )
    integer, intent(in) :: ispin, input_Charge

    integer :: iloop
    integer jj
    Real(kind=DP) :: csum

    Do iloop = 1, ispin
       if(ipri >= 2) write(nfout,'(" ! (m_XC_cal_potential) , iloop = ",i5)') iloop
#ifdef _OLD_MAP_CHARGE_
       call map_charge_onto_a_fft_box(chgq_l)     !-(m_XC_Pot.) -> afft(*) (xcchg2)
#else
       call map_charge_onto_a_fft_box2(nfout,chgq_l)    !-(m_XC_Pot.) -> afft(*) (xcchg2)
#endif
       if(ipri >= 2) then
          write(nfout,'(" ista_fftp, iend_fftp, np_fftp, mp_fftp = ",4i8)' ) &
               & ista_fftp, iend_fftp, np_fftp, mp_fftp
       end if
       if (check_of_xctype() == GGA) then
          chden_l(ista_fftp:iend_fftp,iloop) = afft(ista_fftp:iend_fftp)
       end if

       if(myrank_ggacmp < nrank_ggacmp) then
!          write(*,*) 'DEBUG A', myrank_ggacmp, nrank_ggacmp, afft(1)*univol
          call m_FFT_CD_inverse_c(nfout,afft)  ! afft(G_sp.) -> afft(R_sp.)
          call m_FFT_check_of_negative_CD(npes_cdfft,ista_fftp,iend_fftp&
               & ,ista_fftph,iend_fftph,afft,nfout,nspin,iloop)
       end if
       if(ipri >= 2) write(nfout,'(" ! out of m_FFT_check_of_negative_CD <<m_XC_cal_potential>>")')
       call cp_afft_to_chgrhr    ! -(contained here) afft -> chgrhr_l
    End do
! ------------------------ check --
!    csum = 0.0d0
!    Do jj=1, nfftp /2
!       csum = csum + chgrhr_l(jj,1)
!    End do
!    write(*,*) 'csum = ', csum, csum * univol / dble(nfftp/2)
!    stop                       ! csum *univol /(nffph) == Nelec
! ------------------------------ gga -----
    call check_lmn_even  ! -(contained in subr. m_XC_cal_potential) ->lmn_even

    if(xctype == 'ldapw91' .or. xctype == 'ldapbe ') then
       grad_rho = 0.d0; grad_trho = 0.d0
#ifndef _XC_SAVE_MEMORY_
       cgrad_rho = 0.d0
#endif
    else
       call abs_grad_rho_up_down_total(ispin)
      !   chden_l ->
      !   |grad(rho_up(r))|,|grad(rho_down(r))|  --> grad_rho
      !   |grad(rho_up(r) + rho_down(r))|        --> grad_trho
    end if
! ----------------------------------------------
    !
   contains

    subroutine cp_afft_to_chgrhr
      integer  :: i

      if(myrank_ggacmp < nrank_ggacmp) then
         do i = ista_fftph, iend_fftph    ! MPI
            chgrhr_l(i,iloop) = afft(i*2-1)
         end do
      else
         do i = ista_fftph, iend_fftph    ! MPI
            chgrhr_l(i,iloop) = 0.d0
         end do
      end if

      if(ipri >= 2) then
         write(nfout,'(" !XC -- chgrhr_l -- <<cp_afft_to_chgrhr>>")')
         write(nfout,'(" !XC ",6f12.6)') (chgrhr_l(i,1),i=ista_fftph,ista_fftph+17)
         write(nfout,'(" !XC -- afft --     <<cp_afft_to_chgrhr>>")')
         write(nfout,'(" !XC ",6f12.6)') (afft(i*2-1),i=ista_fftph,ista_fftph+17)
      end if
    end subroutine cp_afft_to_chgrhr

    subroutine map_charge_onto_a_fft_box2(nfout,chgq_l)
      integer, intent(in) :: nfout
      real(DP),intent(in) :: chgq_l(ista_kngp:iend_kngp,kimg,nspin)
      real(kind=DP)      :: fac
      integer            :: mm, it, i, j
      real(kind=DP),allocatable,dimension(:)    :: tmp_s, tmp_r ! MPI d(nfftp/npes)
      real(kind=DP),allocatable,dimension(:,:)  :: chgqplus_l ! MPI d(ista_kngp:iend_kngp)
      integer, allocatable, dimension(:)        :: ip ! d(ista_kngp:iend_kngp,3)
      integer, allocatable, dimension(:)        :: igfp_s, igfp_r ! d(nfftp/npes)
      integer, allocatable, dimension(:,:)      :: igfp_ijk ! d(ista_kngp:iend_kngp+(0,1),3)
      integer, allocatable, dimension(:)        :: ipout_s, ipout_r ! d(nfftp/npes)
      integer, allocatable, dimension(:)        :: np ! d(0:npes_cdfft-1)
      integer, allocatable, dimension(:,:)      :: np_send, np_tmp ! d(ista_kngp:iend_kngp,npes)
      integer :: pe_s, pe_r, datasize_s, datasize_r, nfp, igp, nb, npc, maxdatasize
      integer :: nstart, nend, lx, lxy, ipout, j0, k, l, np0, idp, icdfft
      integer, allocatable, dimension(:) :: req_r,req_s

      integer :: id_sname = -1
      call tstatc0_begin('map_charge_onto_a_fft_box2 ',id_sname)

      allocate(chgqplus_l(ista_kngp:iend_kngp,kimg))

      if(input_charge == Valence_plus_PC_Charge) then
         fac = 1.d0/nspin
      else if(input_charge == Partial_Core_Charge) then
         fac = 1.d0
      else
        stop ' Error of ixc in <<<map_charge_onto_a_fft_box2>>>'
      endif

      chgqplus_l = 0.d0
      mm = 0
      do it = 1, ntyp
         if(itpcc(it) == 0) cycle
         mm = mm + 1
         do j = 1, kimg
#ifdef NEC_TUNE_MXCP
!CDIR NODEP
!CDIR INNER
#endif
            do i = ista_kngp, iend_kngp  !for mpi
               chgqplus_l(i,j) = chgqplus_l(i,j) &
                    & + fac*zfm3_l(i,it,j)*rhpcg_l(i,mm)   !mpi
            end do
         end do
      end do

      if(input_charge == Valence_plus_PC_Charge) then
         do j = 1, kimg
#ifdef NEC_TUNE_MXCP
!CDIR NODEP
!CDIR INNER
#endif
            do i = ista_kngp, iend_kngp  !for mpi
               chgqplus_l(i,j) = chgqplus_l(i,j) + chgq_l(i,j,iloop)
            end do
         end do
      end if

      afft = 0.d0
      if(npes > 1) then
         allocate(ip(ista_kngp:iend_kngp)); ip = -1
         allocate(np(0:npes-1)); np = 0
         allocate(np_send(0:npes-1,0:npes-1))
      end if

#ifdef _MPIFFT_
      if(mod(iend_kngp-ista_kngp+1,2) == 1) then
         allocate(igfp_ijk(ista_kngp:iend_kngp,3))
      else
         allocate(igfp_ijk(ista_kngp:iend_kngp+1,3))
      end if
      idp = fft_box_size_CD_c(1,0)
      do l = ista_kngp, iend_kngp
         np0 = igfp_l_c(l)
         i   = mod(np0-1,idp)+1
         j   = mod((np0-i)/idp, ly) + 1
         k   = (np0-i-(j-1)*idp)/(idp*ly)+1
         igfp_ijk(l,1)  = i
         igfp_ijk(l,2)  = j
         igfp_ijk(l,3)  = k
      end do
      lx     = fft_box_size_CD_c(1,0)
      lxy    = lx*ly_d
#endif

      if(npes > 1) then
         np_send = 0
         np = 0
         do k = 0, npes-1
            if(map_pe2ggacmp(k) >= nrank_ggacmp) cycle
#ifdef _MPIFFT_
            icdfft = map_pe2cdfft(k)
            nstart = icdfft*ny_d+1
            nend   = min(fft_box_size_CD(2,1),(icdfft+1)*ny_d)
            do l = ista_kngp, iend_kngp
               j   = igfp_ijk(l,2)
               if(j >= nstart .and. j <= nend) then
                  if(map_pe2ggacmp(k) == 0) ip(l) = k
                  np(k) = np(k) + 1
               end if
            end do
#else
            icdfft = mod(k,npes_cdfft)
            do l = ista_kngp, iend_kngp
               igp = (igfp_l(l)-1)*kimg+1
               if(igp >= nis_fftp(icdfft) .and. igp <= nie_fftp(icdfft)) then
                  ip(l) = map_pe2cdfft(k)
                  np(k) = np(k) + 1
               end if
            end do
#endif
            if(ipri>=2)  write(nfout,'(" np(",i5,") = ",i8)') k, np(k)
         end do

         np_send(:,mype) = np(:)  ! number element sent from mype to j
         if(npes>1) then
            allocate(np_tmp(0:npes-1,0:npes-1))
            call mpi_allreduce(np_send,np_tmp,npes*npes,mpi_integer, mpi_sum, MPI_CommGroup,ierr)
            np_send = np_tmp
            deallocate(np_tmp)
         end if
         maxdatasize = 0
         do j = 0, npes-1
            do i = 0, npes-1
               if(maxdatasize < np_send(i,j)) maxdatasize = np_send(i,j)
            end do
         end do
         if(ipri>=2) write(nfout,'(" maxdatasize = ",i8)') maxdatasize
      end if

#ifdef _MPIFFT_
      allocate(ipout_s(ista_kngp:iend_kngp))
      do l = ista_kngp, iend_kngp
         i = igfp_ijk(l,1)
         j = igfp_ijk(l,2)
         k = igfp_ijk(l,3)
         icdfft = (j-1)/ny_d
         j0 = j - ny_d*icdfft
         ipout_s(l) = kimg*(i + (j0-1)*lx + (k-1)*lxy - 1)  + nis_fftp(icdfft)
      end do
#endif

      ! mype -> mype
      if(myrank_ggacmp < nrank_ggacmp) then
#ifdef _MPIFFT_
         nstart = myrank_cdfft*ny_d+1
         nend   = min(fft_box_size_CD(2,1),(myrank_cdfft+1)*ny_d)
         if(kimg == 1) then
            if(npes > 1) then
               do l = ista_kngp, iend_kngp
                  i   = igfp_ijk(l,1)
                  if(i < 1 .or. i > fft_box_size_CD(1,1)) cycle
                  j   = igfp_ijk(l,2)
                  if(j >= nstart .and. j <= nend) then
                     ipout = ipout_s(l)
                     afft(ipout) = chgqplus_l(l,1)
                  end if
               end do
            else
               do l = ista_kngp, iend_kngp
                  ipout = ipout_s(l)
                  afft(ipout) = chgqplus_l(l,1)
               end do
            end if
         else
            if(npes > 1) then
               do l = ista_kngp, iend_kngp
                  i   = igfp_ijk(l,1)
                  if(i < 1 .or. i > fft_box_size_CD(1,1)) cycle
                  j   = igfp_ijk(l,2)
                  if(j >= nstart .and. j <= nend) then
                     ipout = ipout_s(l)
                     afft(ipout)   = chgqplus_l(l,1)
                     afft(ipout+1) = chgqplus_l(l,2)
                  end if
               end do
            else
               do l = ista_kngp, iend_kngp
                  ipout = ipout_s(l)
                  afft(ipout)   = chgqplus_l(l,1)
                  afft(ipout+1) = chgqplus_l(l,2)
               end do
            end if
         end if
#else
         if(kimg == 1) then
            do l = ista_kngp, iend_kngp
               igp = igfp_l(l)
               if(igp >= ista_fftp .and. igp <= iend_fftp) then
                  afft(igp) = chgqplus_l(l,1)
               end if
            end do
         else
#ifdef NEC_TUNE_MXCP
!CDIR NODEP
#endif
            do l = ista_kngp, iend_kngp
               igp = (igfp_l(l)-1)*kimg+1
               if(igp >= ista_fftp .and. igp <= iend_fftp) then
                  afft(igp)   = chgqplus_l(l,1)
                  afft(igp+1) = chgqplus_l(l,2)
               end if
            end do
         end if
#endif
      end if

      if(npes > 1) then
         allocate(tmp_s(maxdatasize*kimg))
         allocate(tmp_r(maxdatasize*kimg))
         allocate(igfp_s(maxdatasize),igfp_r(maxdatasize))
         allocate(req_r(npes-1))
         allocate(req_s(npes-1))
         do nb = 1, npes-1
            pe_s = mod(mype+nb,      npes)
            pe_r = mod(mype-nb+npes, npes)
#ifdef _MPIFFT_
            icdfft = map_pe2cdfft(pe_s)
#endif
            if(map_pe2ggacmp(pe_s) < nrank_ggacmp .and. np_send(pe_s,mype) > 0) then
               tmp_s = 0.d0
               igfp_s = 0
               nfp = 0
          if(kimg==1) then
               do l = ista_kngp, iend_kngp
                  if(ip(l) == map_pe2cdfft(pe_s)) then
                     nfp = nfp+1
#ifdef _MPIFFT_
                     igfp_s(nfp) = ipout_s(l)
#else
                     igfp_s(nfp) = kimg*(igfp_l(l)-1)+1
#endif
                        tmp_s(nfp) = chgqplus_l(l,1)
                  end if
               end do
          else
               do l = ista_kngp, iend_kngp
                  if(ip(l) == map_pe2cdfft(pe_s)) then
                     nfp = nfp+1
#ifdef _MPIFFT_
                     igfp_s(nfp) = ipout_s(l)
#else
                     igfp_s(nfp) = kimg*(igfp_l(l)-1)+1
#endif
                        tmp_s(2*nfp-1) = chgqplus_l(l,1)
                        tmp_s(2*nfp)   = chgqplus_l(l,2)
                  end if
               end do
          endif

               datasize_s = np_send(pe_s,mype)*kimg
               if(ipri>=2) then
                  write(nfout,'(" nb, pe_s,pe_r = ",3i8)') nb, pe_s,pe_r
                  write(nfout,'(" nfp*kimg   = ", i8)') nfp*kimg
                  write(nfout,'(" datasize_s = ",i8)') datasize_s
               end if
               if(datasize_s /= nfp*kimg) stop ' datasize_s /= nfp*kimg'
               call mpi_isend(tmp_s,datasize_s,mpi_double_precision, &
                    & pe_s, mype, MPI_CommGroup,req_s(nb),ierr)
            end if
            if(myrank_ggacmp < nrank_ggacmp .and. np_send(mype,pe_r) > 0) then
               datasize_r = np_send(mype,pe_r)*kimg
               if(ipri>=2) write(nfout,'(" datasize_r = ",i8)') datasize_r
               call mpi_irecv(tmp_r,datasize_r,mpi_double_precision, &
                    & pe_r, pe_r,  MPI_CommGroup,req_r(nb),ierr)
            end if
            if(map_pe2ggacmp(pe_s) < nrank_ggacmp .and. np_send(pe_s,mype) > 0) &
                 & call mpi_wait(req_s(nb),istatus,ierr)
            if(myrank_ggacmp < nrank_ggacmp .and. np_send(mype,pe_r) > 0) &
                 & call mpi_wait(req_r(nb),istatus,ierr)
            if(ipri>=2) write(nfout,'( " tmp_s , tmp_r  completed")')

            datasize_s = np_send(pe_s,mype)
            datasize_r = np_send(mype,pe_r)
            if(map_pe2ggacmp(pe_s) < nrank_ggacmp .and. np_send(pe_s,mype) > 0) then
               call mpi_isend(igfp_s,datasize_s,mpi_integer, &
                    & pe_s, mype, MPI_CommGroup, req_s(nb),ierr)
            end if
            if(myrank_ggacmp < nrank_ggacmp .and. np_send(mype,pe_r) > 0) then
               call mpi_irecv(igfp_r,datasize_r,mpi_integer, &
                    & pe_r, pe_r,  MPI_CommGroup, req_r(nb),ierr)
            end if
            if(map_pe2ggacmp(pe_s) <nrank_ggacmp .and. np_send(pe_s,mype) > 0) &
                 & call mpi_wait(req_s(nb),istatus,ierr)
            if(myrank_ggacmp < nrank_ggacmp .and. np_send(mype,pe_r) > 0) &
                 & call mpi_wait(req_r(nb),istatus,ierr)
            if(ipri>=2) write(nfout,'( " igfp_s , igfp_r  completed")')

            if(myrank_ggacmp < nrank_ggacmp) then
               if(kimg == 1) then
                  do l = 1, np_send(mype,pe_r)
                     ipout = igfp_r(l)
                     afft(ipout) = tmp_r(l)
                  end do
               else
#ifdef NEC_TUNE_MXCP
!CDIR NODEP
#endif
                  do l = 1, np_send(mype,pe_r)
                     ipout = igfp_r(l)
                     afft(ipout)   = tmp_r(2*l-1)
                     afft(ipout+1) = tmp_r(2*l)
                  end do
               end if
            end if
         end do
         deallocate(req_s, req_r)
         deallocate(igfp_r, igfp_s)
         deallocate(tmp_r,tmp_s)
      end if
#ifdef _MPIFFT_
      deallocate(ipout_s)
      deallocate(igfp_ijk)
#endif
      if(npes > 1) then
         deallocate(np_send,np, ip)
      end if
      deallocate(chgqplus_l)

      call tstatc0_end(id_sname)

    end subroutine map_charge_onto_a_fft_box2

  end subroutine Calc_Density_in_Rspace

  subroutine afft_allgatherv(afft_mpi0,afft)
    !  Upgraded on 19 Sep. 2006 by T. Yamasaki
    !    * MPIFFT
    real(kind=DP),intent(in),dimension(:) :: afft_mpi0(ista_fftp:iend_fftp)
    real(kind=DP),intent(out),dimension(:):: afft(nfftp)
#ifdef _MPIFFT_
    real(kind=DP), allocatable, dimension(:) :: afft_tmp
    integer  :: nstart, nend, lx, lxy, llx, lly, llxy, i, j, jj, k, ipin, ipout
    integer  :: nb, datasize, pe_s, pe_r
    integer, allocatable, dimension(:) :: req_s, req_r
    real(kind=DP), allocatable, dimension(:) :: tmp_s, tmp_r
#endif
    integer  :: id_sname = -1
!----------------------------- start ------------
    call tstatc0_begin('afft_allgatherv(in m_XC_Potential) ',id_sname)
    if(npes>=2)  call mpi_barrier(MPI_CommGroup,ierr)

    if(npes_cdfft >= 2) then
       if(myrank_ggacmp < nrank_ggacmp) then
#ifdef _MPIFFT_
          afft = 0.d0
#ifndef _ALLREDUCE_AFFT_ALLGAHTERV_

          datasize = nfftp/npes_cdfft
          allocate(tmp_s(datasize)); tmp_s = 0.d0
          allocate(tmp_r(datasize)); tmp_r = 0.d0
          allocate(req_r(npes-1),req_s(npes-1))

          lx     = fft_box_size_CD_c(1,0)
          lxy    = lx*ly_d
          llx    = fft_box_size_CD_c(1,0)
          lly    = fft_box_size_CD_c(2,0)
          llxy   = llx*lly

          do i = 1, iend_fftp-ista_fftp+1
             tmp_s(i) = afft_mpi0(i+ista_fftp-1)
          end do

          do nb = 1, npes_cdfft-1
             pe_s = mod(myrank_cdfft+nb,           npes_cdfft)
             pe_r = mod(myrank_cdfft-nb+npes_cdfft,npes_cdfft)
             call mpi_irecv(tmp_r,datasize,mpi_double_precision, &
                  & pe_r,pe_r,        mpi_cdfft_world(myrank_ggacmp),req_r(nb),ierr)
             call mpi_isend(tmp_s,datasize,mpi_double_precision, &
                  & pe_s,myrank_cdfft,mpi_cdfft_world(myrank_ggacmp),req_s(nb),ierr)
             call mpi_wait(req_r(nb),istatus,ierr)
             call mpi_wait(req_s(nb),istatus,ierr)

             nstart = pe_r*ny_d+1
             nend   = min(fft_box_size_CD(2,1),(pe_r+1)*ny_d)

             if(kimg == 1) then
                do k = 1, fft_box_size_CD(3,1)
                   do j = nstart, nend
                      jj = j - nstart + 1
                      do i = 1, fft_box_size_CD(1,1)
                         ipin  = i + (jj-1)* lx + (k-1)* lxy
                         ipout = i + ( j-1)*llx + (k-1)*llxy
                         afft(ipout)   = tmp_r(ipin)
                      end do
                   end do
                end do
             else if(kimg == 2) then
                do k = 1, fft_box_size_CD(3,1)
                   do j = nstart, nend
                      jj = j - nstart + 1
                      do i = 1, fft_box_size_CD(1,1)
                         ipin  = 2*(i + (jj-1)* lx + (k-1)* lxy)-1
                         ipout = 2*(i + ( j-1)*llx + (k-1)*llxy)-1
                         afft(ipout)   = tmp_r(ipin)
                         afft(ipout+1) = tmp_r(ipin+1)
                      end do
                   end do
                end do
             end if
          end do
          deallocate(req_s,req_r)
          deallocate(tmp_s)
          deallocate(tmp_r)

! diagonal part
          nstart = myrank_cdfft*ny_d+1
          nend   = min(fft_box_size_CD(2,1),(myrank_cdfft+1)*ny_d)
!!$          allocate(afft_tmp(nfftp))

          if(kimg == 1) then
             do k = 1, fft_box_size_CD(3,1)
                do j = nstart, nend
                   jj = j - nstart + 1
                   do i = 1, fft_box_size_CD(1,1)
                      ipin  = ista_fftp-1 + (i + (jj-1)* lx + (k-1)* lxy)
                      ipout =               (i + ( j-1)*llx + (k-1)*llxy)
                      afft(ipout)   = afft_mpi0(ipin)
                   end do
                end do
             end do
          else if(kimg == 2) then
             do k = 1, fft_box_size_CD(3,1)
                do j = nstart, nend
                   jj = j - nstart + 1
                   do i = 1, fft_box_size_CD(1,1)
                      ipin  = ista_fftp-1 + 2*(i + (jj-1)* lx + (k-1)* lxy)-1
                      ipout =               2*(i + ( j-1)*llx + (k-1)*llxy)-1
                      afft(ipout)   = afft_mpi0(ipin)
                      afft(ipout+1) = afft_mpi0(ipin+1)
                   end do
                end do
             end do
          end if
#else
          nstart = myrank_cdfft*ny_d+1
          nend   = min(fft_box_size_CD(2,1),(myrank_cdfft+1)*ny_d)
          lx     = fft_box_size_CD_c(1,0)
          lxy    = lx*ly_d
          llx    = fft_box_size_CD_c(1,0)
          lly    = fft_box_size_CD_c(2,0)
          llxy   = llx*lly
          allocate(afft_tmp(nfftp))

          if(kimg == 1) then
             do k = 1, fft_box_size_CD(3,1)
                do j = nstart, nend
                   jj = j - nstart + 1
                   do i = 1, fft_box_size_CD(1,1)
                      ipin  = ista_fftp-1 + (i + (jj-1)* lx + (k-1)* lxy)
                      ipout =               (i + ( j-1)*llx + (k-1)*llxy)
                      afft(ipout)   = afft_mpi0(ipin)
                   end do
                end do
             end do
          else if(kimg == 2) then
             do k = 1, fft_box_size_CD(3,1)
                do j = nstart, nend
                   jj = j - nstart + 1
                   do i = 1, fft_box_size_CD(1,1)
                      ipin  = ista_fftp-1 + 2*(i + (jj-1)* lx + (k-1)* lxy)-1
                      ipout =               2*(i + ( j-1)*llx + (k-1)*llxy)-1
                      afft(ipout)   = afft_mpi0(ipin)
                      afft(ipout+1) = afft_mpi0(ipin+1)
                   end do
                end do
             end do
          end if
          call mpi_allreduce(afft, afft_tmp, nfftp, mpi_double_precision, mpi_sum, mpi_cdfft_world(myrank_ggacmp),ierr)
          afft = afft_tmp
#endif

#else
          call mpi_allgatherv(afft_mpi0,nel_fftp(myrank_cdfft),mpi_double_precision, afft,nel_fftp,idisp_fftp,&
          & mpi_double_precision,mpi_cdfft_world(myrank_ggacmp),ierr)
#endif
       end if
    else
       afft = afft_mpi0
    end if

    if(nrest_cdfft >= 1) then
       if(mype > npes_cdfft*nrank_ggacmp-1) then
          call mpi_recv(afft,nfftp,mpi_double_precision &
               & , mype-npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,istatus,ierr)
       end if
       if(mype < nrest_cdfft) then
          call mpi_send(afft,nfftp,mpi_double_precision &
               & , mype+npes_cdfft*nrank_ggacmp,1,MPI_CommGroup,ierr)
       end if
    end if

    call tstatc0_end(id_sname)

  end subroutine afft_allgatherv

  subroutine xc_allocate(ispin,nfout)
    integer, intent(in)      :: ispin,nfout

    integer                  :: idp,nlp,nmp,nnp,nd2p,nd3p,ip, idph, nlph
    integer                  :: n, nn, i, j, k
#ifdef _MPIFFT_
    integer                  :: np0, j0, i2, j2, k2, n0
#endif

#ifndef _MPIFFT_
    call m_FFT_alloc_CD_box()
#endif
    allocate(afft(ista_fftp:iend_fftp));afft=0.d0
    allocate(chgrhr_l(ista_fftph:iend_fftph,ispin))    ! MPI
    allocate(f2or1(ista_fftph:iend_fftph))

    nlp  = fft_box_size_CD(1,1)
    nmp  = fft_box_size_CD(2,1)
    nnp  = fft_box_size_CD(3,1)
#ifdef _MPIFFT_
    idp  = fft_box_size_CD_c(1,0)
    nd2p = fft_box_size_CD_c(2,0)
    nd3p = fft_box_size_CD_c(3,0)
#else
    idp  = fft_box_size_CD(1,0)
    nd2p = fft_box_size_CD(2,0)
    nd3p = fft_box_size_CD(3,0)
#endif

    call set_f2or1(npes_cdfft,ista_fftph,iend_fftph,f2or1)

    if(check_of_xctype() == GGA) then
       allocate(inx(ista_fftp:iend_fftp)); inx = 0
       allocate(jnx(ista_fftp:iend_fftp)); jnx = 0
       allocate(knx(ista_fftp:iend_fftp)); knx = 0

!
!       write(*,*) 'Maji Maji ', inx
!       stop
#ifdef _MPIFFTTESTTEST_
       if(kimg == 1) then
          do n = ista_fftp, iend_fftp
             nn = n
             if(nn == ista_fftp-1) write(nfout,'(" n == ista_fftp-1 = ",i8)') nn
             np0 = nn - idp*ly_d*lz*myrank_cdfft
             i  = mod(np0-1,idp)+1
             j0 = mod((np0-i)/idp,ly_d) + 1
             j  = j0 + ny_d*myrank_cdfft
             k  = (np0-i-(j0-1)*idp)/(idp*ly_d) + 1
             if(i >= idp) i = i - idp
             inx(n) = i - 1;  if(2*inx(n) > nlp) inx(n) = inx(n) - nlp
             jnx(n) = j - 1;  if(2*jnx(n) > nmp) jnx(n) = jnx(n) - nmp
             knx(n) = k - 1;  if(2*knx(n) > nnp) knx(n) = knx(n) - nnp
          end do

!!$          do k = 0, nnp-1
!!$             k2 = k; if(2*k2 > nnp) k2 = k2 - nnp
!!$             do j0 = 0, min(ny_d, nmp-ny_d*myrank_cdfft)-1
!!$                j = j0 + ny_d*myrank_cdfft
!!$                j2 = j; if(2*j2 > nmp) j2 = j2 - nmp
!!$                n0 = 1+j0*idp + k*idp*ly_d + ista_fftp-1
!!$                do i = 0, nlp-1
!!$                   n = i+n0
!!$                   inx(n) = i; if(2*inx(n) > nlp) inx(n) = inx(n) - nlp
!!$                   if(n == ista_fftp-1) write(nfout,'(" n = ista_fftp-1 = ",i8)') n
!!$                   jnx(n) = j2
!!$                   knx(n) = k2
!!$                end do
!!$             end do
!!$          end do

!!$          do i = 0, nlp-1
!!$             i2 = i; if(2*i2 > nlp) i2 = i2 - nlp
!!$             do j0 = 0, min(ny_d, nmp-ny_d*myrank_cdfft)-1
!!$                j = j0 + ny_d*myrank_cdfft
!!$                j2 = j; if(2*j2 > nmp) j2 = j2 - nmp
!!$                do k = 0, nnp/2
!!$                   n = i+1 + j0*idp + k*idp*ly_d + ista_fftp-1
!!$                   inx(n) = i2
!!$                   jnx(n) = j2
!!$                   knx(n) = k
!!$                end do
!!$                do k = nnp/2+1, nnp-1
!!$                   n = i+1 + j0*idp + k*idp*ly_d + ista_fftp-1
!!$                   inx(n) = i2
!!$                   jnx(n) = j2
!!$                   knx(n) = k - nnp
!!$                end do
!!$             end do
!!$          end do
       else if(kimg == 2) then
          do k = 0, nnp-1
             k2 = k; if(2*k2 > nnp) k2 = k2 - nnp
             do j0 = 0, min(ny_d, nmp-ny_d*myrank_cdfft)-1
                j = j0 + ny_d*myrank_cdfft
                j2 = j; if(2*j2 > nmp) j2 = j2 - nmp
                do i = 0, nlp/2
                   n = i+1 + j0*idp + k*idp*ly_d
                   nn = n*2-1 + ista_fftp-1
                   inx(nn) = i ; inx(nn+1) = i
                   jnx(nn) = j2; jnx(nn+1) = j2
                   knx(nn) = k2; knx(nn+1) = k2
                end do
                do i = nlp/2+1, nlp-1
                   n = i+1 + j0*idp + k*idp*ly_d
                   nn = n*2-1 + ista_fftp-1
                   inx(nn) = i-nlp ; inx(nn+1) = i-nlp
                   jnx(nn) = j2; jnx(nn+1) = j2
                   knx(nn) = k2; knx(nn+1) = k2
                end do
             end do
          end do
       end if
#else
       do n = ista_fftp, iend_fftp
          nn = (n+kimg-1)/kimg    ! nn = n (when kimg=1), nn = (n+1)/2 (when kimg=2)
#ifdef _MPIFFT_
!!$          if(nn == ista_fftp-1) write(nfout,'(" n == ista_fftp-1 = ",i8)') nn
          np0 = nn - idp*ly_d*lz*myrank_cdfft
          i  = mod(np0-1,idp)+1
          j0 = mod((np0-i)/idp,ly_d) + 1
          j  = j0 + ny_d*myrank_cdfft
          k  = (np0-i-(j0-1)*idp)/(idp*ly_d) + 1
          if(i >= idp) i = i - idp
#else
          i  = mod(nn,idp)
          j  = mod((nn-1)/idp,nd2p) + 1
          k  = (nn - (j-1)*idp - i)/(idp*nd2p) + 1
#endif
          inx(n) = i - 1;  if(2*inx(n) > nlp) inx(n) = inx(n) - nlp
          jnx(n) = j - 1;  if(2*jnx(n) > nmp) jnx(n) = jnx(n) - nmp
          knx(n) = k - 1;  if(2*knx(n) > nnp) knx(n) = knx(n) - nnp
       end do
#endif

!       write(*,*) 'Ua ', mype, inx(1), jnx(1), knx(1)
!       write(*,*) 'Va ', mype, inx(2), jnx(2), knx(2)
!       write(*,*) 'Ub ', mype, inx(3), jnx(3), knx(3)
!       write(*,*) 'Vb ', mype, inx(4), jnx(4), knx(4)
!       stop

       allocate(chden_l(ista_fftp:iend_fftp,nspin)); chden_l = 0.d0
       allocate(grad_trho(ista_fftph:iend_fftph));   grad_trho = 0.0d0
       allocate(grad_rho(ista_fftph:iend_fftph,nspin)); grad_rho = 0.d0

#ifndef _XC_SAVE_MEMORY_
       if(nrank_ggacmp > 1) then
          allocate(cgrad_rho(ista_fftph:iend_fftph,myrank_ggacmp+1:myrank_ggacmp+1,nspin))  ! MPI
       else
          allocate(cgrad_rho(ista_fftph:iend_fftph,3,nspin))  ! MPI
       end if
#else
       allocate(cggawk13(ista_fftph:iend_fftph))           ! MPI
#endif
!!       allocate(dF_drho(ista_fftph:iend_fftph,nspin))      ! MPI
!!       allocate(dF_dgradrho(ista_fftph:iend_fftph,nspin))  ! MPI
    end if

  contains
    subroutine set_f2or1(npes,ista,iend,f2or1)
      integer, intent(in) :: npes,ista, iend
      real(kind=DP), intent(out), dimension(ista:iend) :: f2or1
      integer :: idph,nlph,ip,i,j,k

      if(kimg == 1) then
         idph = idp/2
         nlph = nlp/2
#ifdef _MPIFFT_
         f2or1 = 0.d0
!!$         do j = 1, min(nz_d,nnp-nz_d*myrank_cdfft)*nmp
!!$         do j = 1, min(nz_d,nnp-nz_d*myrank_cdfft)*nd2p
         do k = 1, min(nz_d, nnp-nz_d*myrank_cdfft)
            do j = 1, nmp
               do i = 1, nlph
!!$                  ip = i + idph*(j-1) + idph*ly*lz_d*myrank_cdfft
                  ip = i + idph*(j-1) + idph*ly*(k-1) + idph*ly*lz_d*myrank_cdfft
                  f2or1(ip) = 2.d0
               end do
               ip = 1 + idph*(j-1) + idph*ly*(k-1) + idph*ly*lz_d*myrank_cdfft
               f2or1(ip) = 1.d0
               ip = nlph+1 + idph*(j-1) + idph*ly*(k-1) + idph*ly*lz_d*myrank_cdfft
               f2or1(ip) = 1.d0
            end do
         end do
!!$            ip = idph*(j-1) + 1 + idph*ly*lz_d*myrank_cdfft
!!$            f2or1(ip) = 1.d0
!!$            ip = idph*(j-1)+ nlph + 1 + idph*ly*lz_d*myrank_cdfft
!!$            f2or1(ip) = 1.d0
!!$         end do
#else
         f2or1 = 2.d0
         if(npes >= 2) then
            do i = 1, nd2p*nnp
               ip = idph*(i-1) + 1
               if(ip>= ista .and. ip <= iend) f2or1(ip) = 1.d0
            end do
            do i = 1, nd2p*nnp
               ip = idph*(i-1) + nlph + 1
               if(ip>= ista .and. ip <= iend) f2or1(ip) = 1.d0
            end do
            do j = nlph+2, idph
               do i = 1,nd2p*nnp
                  ip = idph*(i-1)+j
                  if(ip>= ista .and. ip <= iend) f2or1(ip) = 0.d0
               end do
            end do
            do j = nmp+1, nd2p
               do k = 1, nnp
                  do i = 1, nlph
                     ip = i + idph*(j-1) + idph*nd2p*(k-1)
                     if(ip >= ista .and. ip <= iend) f2or1(ip) = 0.d0
                  end do
               end do
            end do
! for SX6 ASL 20040817
            do k = nnp+1, nd3p
               do i = 1, idph*nd2p
                  ip = i + idph*nd2p*(k-1)
                  if(ip >= ista .and. ip <= iend) f2or1(ip) = 0.d0
               end do
            end do
         else
            do i = 1, nd2p*nnp
               ip = idph*(i-1) + 1
               f2or1(ip) = 1.d0
            end do
            do i = 1, nd2p*nnp
               ip = idph*(i-1) + nlph + 1
               f2or1(ip) = 1.d0
            end do
            do j = nlph+2, idph
               do i = 1,nd2p*nnp
                  ip = idph*(i-1)+j
                  f2or1(ip) = 0.d0
               end do
            end do
            do j = nmp+1, nd2p
               do k = 1, nnp
                  do i = 1, nlph
                     ip = i + idph*(j-1) + idph*nd2p*(k-1)
                     f2or1(ip) = 0.d0
                  end do
               end do
            end do
! for SX6 ASL 20040817
            do k = nnp+1, nd3p
               do i = 1, idph*nd2p
                  ip = i + idph*nd2p*(k-1)
                  f2or1(ip) = 0.d0
               end do
            end do
         end if
!!$       do i = ista, iend
!!$          if(mod(i*2,idp) == 2 .or. mod(i*2,idp) == 0) f2or1(i) = 1.d0
!!$       end do
#endif
      else
#ifdef _MPIFFT_
!!$         f2or1 = 0.d0                               ! f2or1 works to the fft data in R space.
!!$         do k = 1, min(nz_d,nnp-nz_d*myrank_cdfft)
!!$            do j = 1, nmp     ! nmp = fft_box_size_CD(2,1)
!!$               do i = 1, nlp  ! nlp = fft_box_size_CD(1,1)
!!$                  ip = i+(j-1)*idp+(k-1)*idp*nd2p+idp*nd2p*lz_d*myrank_cdfft
!!$                  f2or1(ip) = 1.d0
!!$               end do
!!$            end do
!!$         end do
         if(ipri >= 2 ) write(nfout,'(" ix kimg = 2 <<set_f2or1>>")')
         f2or1 = 1.d0
         do j = nlp+1, idp      ! x
            do i = 1, ly*nz_d
               ip = idp*(i-1)+j+ista-1
               f2or1(ip) = 0.d0
            end do
         end do
         if(ipri >= 2 ) write(nfout,'(" iy kimg = 2 <<set_f2or1>>")')
         do j = nmp+1,ly         ! y
            do k = 1, nz_d
               do i = 1, nlp
                  ip = i + idp*(j-1) + idp*ly*(k-1) + ista-1
                  f2or1(ip) = 0.d0
               end do
            end do
         end do
         if(ipri >= 2 ) write(nfout,'(" iz kimg = 2 <<set_f2or1>>")')
         do  k = nz_d+1, lz_d   ! z
            do i = 1, idp*ly
               ip = i + idp*ly*(k-1) + ista-1
               f2or1(ip) = 0.d0
            end do
         end do
#else
         f2or1 = 1.d0
         if(npes >= 2) then
            do j = nlp+1, idp    ! x
               do i = 1, nd2p*nnp
                  ip = idp*(i-1)+j
                  if(ip>= ista .and. ip <= iend) f2or1(ip) = 0.d0
               end do
            end do
            do j = nmp+1, nd2p   ! y
               do k = 1, nnp
                  do i = 1, nlp
                     ip = i + idp*(j-1) + idp*nd2p*(k-1)
                     if(ip >= ista .and. ip <= iend) f2or1(ip) = 0.d0
                  end do
               end do
            end do
! for SX6 ASL 20040817
            do k = nnp+1, nd3p   ! z
               do i = 1, idp*nd2p
                  ip = i + idp*nd2p*(k-1)
                  if(ip >= ista .and. ip <= iend) f2or1(ip) = 0.d0
               end do
            end do
         else
            do j = nlp+1, idp    ! x
               do i = 1, nd2p*nnp
                  ip = idp*(i-1)+j
                  f2or1(ip) = 0.d0
               end do
            end do
            do j = nmp+1, nd2p   ! y
               do k = 1, nnp
                  do i = 1, nlp
                     ip = i + idp*(j-1) + idp*nd2p*(k-1)
                     f2or1(ip) = 0.d0
                  end do
               end do
            end do
! for SX6 ASL 20040817
            do k = nnp+1, nd3p   ! z
               do i = 1, idp*nd2p
                  ip = i + idp*nd2p*(k-1)
                  f2or1(ip) = 0.d0
               end do
            end do
         end if
#endif
      end if
    end subroutine set_f2or1
  end subroutine xc_allocate

  subroutine xc_deallocate
#ifndef _MPIFFT_
    call m_FFT_dealloc_CD_box()
#endif
    deallocate(afft)
    deallocate(chgrhr_l)
    deallocate(f2or1)

    if(check_of_xctype() == GGA) then
       deallocate(inx); deallocate(jnx); deallocate(knx)
       deallocate(chden_l)
       deallocate(grad_trho)
#ifndef _XC_SAVE_MEMORY_
       deallocate(cgrad_rho)
#else
       deallocate(cggawk13)
#endif
       deallocate(grad_rho);
!!        deallocate(dF_drho)
!!       deallocate(dF_dgradrho)

    end if
  end subroutine xc_deallocate

! ------------------------- GGA --------------------------
  subroutine abs_grad_rho_up_down_total(ispin)
    integer, intent(in) :: ispin
    real(kind=DP),allocatable,dimension(:)    :: grad_rho_c2  ! MPI d(ista_fftph:iend_fftph)
    integer  :: is, in, i
    real(kind=DP) :: x,y,z

#ifdef _DETAIL_GGA_TIMING_
    integer :: id_sname = -1
#endif
    if(nrank_ggacmp > 1) allocate(grad_rho_c2(ista_fftph:iend_fftph))    ! MPI

    do is = 1, ispin
       grad_rho(:,is) = 0.d0
       do in = 1, 3
          if(map_ggacmp(in) /= myrank_ggacmp) cycle
          call g_xyz_chden_l(in,is)        ! G_xyz * rho_{up|down}(G)-->afft
          if(ipri >= 2) call wd_small_part_of_afft_kt(34,'G space before<m_FFT_CD_inverse_c)',120, is )
          call m_FFT_CD_inverse_c(nfout,afft)! (-i)*d(rho_{up|down}(r))/d(x|y|z)
          if(ipri >= 2) call wd_small_part_of_afft_kt(34,'R space after <m_FFT_CD_inverse_c)',120, is)
#ifndef _XC_SAVE_MEMORY_
          call cp_afft_to_cgrad_rho(is,in) ! -> cgrad_rho(i,in,is) = -afft(i*2)
#endif
          call add_sq_afft_to_grad_rho(is) ! grad_rho <--  + afft**2
          if(ipri >= 2) then
             write(nfout,'(" !XC after add_sq_afft_to_grad_rho: is, in = ",2i8)') is, in
             write(nfout,'(" !XC -- grad_rho -- <<abs_grad_rho_up_down_total>>")')
             write(nfout,'(" !XC ",6f12.6)') (grad_rho(i,1),i=ista_fftph,ista_fftph+17)
             write(nfout,'(" !XC ",6f12.6)') (grad_rho(i,1),i=iend_fftph-17,iend_fftph)
#ifndef _XC_SAVE_MEMORY_
             write(nfout,'(" !XC -- cgrad_rho -- <<abs_grad_rho_up_down_total>>")')
             write(nfout,'(" !XC ",6f12.6)') (cgrad_rho(i,in,1),i=ista_fftph,ista_fftph+17)
             write(nfout,'(" !XC ",6f12.6)') (cgrad_rho(i,in,1),i=iend_fftph-17,iend_fftph)
#endif
          end if

       end do

       if(ipri >= 2) write(nfout,'(" !XC after add_sq_afft_to_grad_rho")')
! --> grad_rho_c
       if(nrank_ggacmp > 1) then
          grad_rho_c2(:) = grad_rho(:,is)
#ifdef _DETAIL_GGA_TIMING_
          call tstatc0_begin('mpi_allreduce(ggaxcp) ',id_sname)
#endif
          call mpi_allreduce(grad_rho_c2,grad_trho,iend_fftph-ista_fftph+1 &
               & ,mpi_double_precision, mpi_sum, mpi_ggacmp_cross_world(myrank_cdfft),ierr)
#ifdef _DETAIL_GGA_TIMING_
          call tstatc0_end(id_sname)
#endif
          grad_rho(:,is) = dsqrt(grad_trho(:))
       else
          grad_rho(:,is) = dsqrt(grad_rho(:,is)) ! grad_rho(s) <-- |grad(rho_s(r))|
       end if
    enddo
    if(ispin == 2) then
#ifndef _XC_SAVE_MEMORY_
       if(nrank_ggacmp > 1) then
          grad_rho_c2 = 0.d0
!!$            if(myrank_ggacmp >= 0 .and. myrank_ggacmp <=2) then
          do i = ista_fftph, iend_fftph
             x = cgrad_rho(i,myrank_ggacmp+1,1) + cgrad_rho(i,myrank_ggacmp+1,2)
             grad_rho_c2(i) = grad_rho_c2(i) + x*x
          end do
#ifdef _DETAIL_GGA_TIMING_
          call tstatc0_begin('mpi_allreduce(ggaxcp) ',id_sname)
#endif
          call mpi_allreduce(grad_rho_c2,grad_trho,iend_fftph-ista_fftph+1 &
               & , mpi_double_precision, mpi_sum, mpi_ggacmp_cross_world(myrank_cdfft),ierr)
#ifdef _DETAIL_GGA_TIMING_
          call tstatc0_end(id_sname)
#endif
          grad_trho = dsqrt(grad_trho)
       else
          do in = ista_fftph, iend_fftph
             x = cgrad_rho(in,1,1) + cgrad_rho(in,1,2)
             y = cgrad_rho(in,2,1) + cgrad_rho(in,2,2)
             z = cgrad_rho(in,3,1) + cgrad_rho(in,3,2)
             grad_trho(in) = dsqrt(x*x+y*y+z*z)
          end do
       end if
#else
       do in = 1, 3
          call g_xyz_total_chden_l(in)     ! G_xyz*(rho(G)up+rho(G)down) -> afft
          call m_FFT_CD_inverse_c(nfout,afft)!(-i)*d(rho_total(r))/d(x|y|z)
          call add_sq_afft_to_grad_trho
       end do
       grad_trho = dsqrt(grad_trho)     ! grad_trho <-- |grad(rho_total(r))|
#endif
    else
       grad_trho(:) = grad_rho(:,1)
    end if
    if(nrank_ggacmp > 1) deallocate(grad_rho_c2)
  end subroutine abs_grad_rho_up_down_total

  subroutine g_xyz_chden_l(in,is)
    integer, intent(in) :: in,is
    integer       :: n
    real(kind=DP) :: gxyz

!!$      deallocate(afft)
!!$      allocate(afft(ista_fftp:iend_fftp))

    afft = 0.d0
    do n = ista_fftp, iend_fftp  ! MPI
       gxyz = rltv(in,1)*inx(n)+rltv(in,2)*jnx(n)+rltv(in,3)*knx(n)
       afft(n) = gxyz*chden_l(n,is)
    enddo
    call boundary_zero_into_afft(in)  ! -(contained in subr. m_XC_cal_potential)
  end subroutine g_xyz_chden_l

  subroutine check_lmn_even
    integer             :: nlmn, i

    do i = 1, 3
       nlmn = fft_box_size_CD(i,1)/2
       if(2*nlmn == fft_box_size_CD(i,1)) then
          lmn_even(i) = .true.
       else
          lmn_even(i) = .false.
       end if
    end do
  end subroutine check_lmn_even

  subroutine boundary_zero_into_afft(in)
    integer, intent(in) :: in
    integer             :: i,j,k,nn,n,idp,nlp,nmp,nnp,nd2p, j0

    nlp = fft_box_size_CD(1,1)
    nmp = fft_box_size_CD(2,1)
    nnp = fft_box_size_CD(3,1)
#ifdef _MPIFFT_
    idp = fft_box_size_CD_c(1,0)
    nd2p = fft_box_size_CD_c(2,0)
#else
    idp = fft_box_size_CD(1,0)
    nd2p = fft_box_size_CD(2,0)
#endif

    if(kimg == 1) then
       if(lmn_even(in)) then
          if( in == 1) then
#ifdef _MPIFFT_
             do j = 1, ly_d*lz
!!$               do j0 = 1, ly_d
!!$                  do k = 1, nnp
!!$                     nn = nlp/2 + 1 + idp*(j0-1) + idp*ly_d*(k-1) + idp*ly_d*lz*myrank_cdfft
                nn = nlp/2 + 1 + idp*(j-1) + idp*ly_d*lz*myrank_cdfft
                afft(nn) = 0.d0
!!$                  end do
             end do
#else
             do j = 1, nmp      ! y
                do k = 1, nnp   ! z
!!$                  nn = nlp/2 + 1 + idp*(j-1) + idp*nmp*(k-1)
                   nn = nlp/2 + 1 + idp*(j-1) + idp*nd2p*(k-1)
                   if(nn >= ista_fftp .and. nn <= iend_fftp) afft(nn) = 0.d0
                end do
             end do
#endif
          else if( in == 2) then
#ifdef _MPIFFT_
             j0 = nmp/2 + 1 - ny_d*myrank_cdfft
             if(j0 >= 1 .and. j0 <= ny_d) then
                do i = 1, idp
                   do k = 1, nnp
                      nn = i + idp*(j0-1) + idp*ly_d*(k-1) + idp*ly_d*lz*myrank_cdfft
                      afft(nn) = 0.d0
                   end do
                end do
             end if
#else
             do i = 1, idp      ! x
                do k = 1, nnp   ! z
                   nn = i + idp*(nmp/2) + idp*nd2p*(k-1)
                   if(nn >= ista_fftp .and. nn <= iend_fftp) afft(nn) = 0.d0
                end do
             end do
#endif
          else if(in == 3) then
#ifdef _MPIFFT_
             do i = 1, idp
                do j = 1+ny_d*myrank_cdfft, min(nmp,ny_d*(myrank_cdfft+1))
                   nn = i + idp*(j-1-ny_d*myrank_cdfft) + idp*ly_d*(nnp/2) &
                        &                          + idp*ly_d*lz*myrank_cdfft
                   afft(nn) = 0.d0
                end do
             end do
#else
             do i = 1, idp      ! x
                do j = 1, nmp   ! y
                   nn = i + idp*(j-1) + idp*nd2p*(nnp/2)
                   if(nn >= ista_fftp .and. nn <= iend_fftp) afft(nn) = 0.d0
                end do
             end do
#endif
          end if
       end if
    else if(kimg == 2) then
       if(lmn_even(in)) then
          if( in == 1) then
#ifdef _MPIFFT_
             do j = 1, ly_d*lz
                nn = nlp/2 + 1 + idp*(j-1) + idp*ly_d*lz*myrank_cdfft
                n = nn*2 - 1
                afft(n) = 0.d0
                afft(n+1) = 0.d0
             end do
#else
             do j = 1, nmp
                do k = 1, nnp
!!$                  nn = nlp/2 + 1 + idp*(j-1) + idp*nmp*(k-1)
                   nn = nlp/2 + 1 + idp*(j-1) + idp*nd2p*(k-1)
!!$                     nn = nlp/2 + 1 + idp*(j-1) -idp*ly_d*myrank_cdfft &
!!$                          & +idp*ly_d*(k-1+myrank_cdfft)
                   n  = nn*2 - 1
                   if(n >= ista_fftp .and. n <= iend_fftp) afft(n) = 0.d0
                   n  = nn*2
                   if(n >= ista_fftp .and. n <= iend_fftp) afft(n) = 0.d0
                end do
             end do
#endif
          else if( in == 2) then
#ifdef _MPIFFT_
             j0 = nmp/2 + 1 - ny_d*myrank_cdfft
             if(j0 >= 1 .and. j0 <= ny_d) then
                do i = 1, idp
                   do k = 1, nnp
                      nn = i + idp*(j0-1) + idp*ly_d*(k-1) + idp*ly_d*lz*myrank_cdfft
                      n  = nn*2 - 1
                      afft(n) = 0.d0; afft(n+1) = 0.d0
                   end do
                end do
             end if
#else
             do i = 1, idp
                do k = 1, nnp
!!$                     nn = i + idp*(nmp/2) + idp*nmp*(k-1)
                   nn = i + idp*(nmp/2) + idp*nd2p*(k-1)
                   n  = nn*2 - 1
                   if(n >= ista_fftp .and. n <= iend_fftp) afft(n) = 0.d0
                   n  = nn*2
                   if(n >= ista_fftp .and. n <= iend_fftp) afft(n) = 0.d0
                end do
             end do
#endif
          else if(in == 3) then
#ifdef _MPIFFT_
             do i = 1, idp
                do j = 1+ny_d*myrank_cdfft, min(nmp,ny_d*(myrank_cdfft+1))
                   nn = i + idp*(j-1-ny_d*myrank_cdfft) + idp*ly_d*(nnp/2) &
                        &                          + idp*ly_d*lz*myrank_cdfft
                   n  = nn*2 - 1
                   afft(n) = 0.d0
                   afft(n+1) = 0.d0
                end do
             end do
#else
             do i = 1, idp
                do j = 1, nmp
!!$                  nn = i + idp*(j-1) + idp*nmp*(nnp/2)
                   nn = i + idp*(j-1) + idp*nd2p*(nnp/2)
                   n  = nn*2 - 1
                   if(n >= ista_fftp .and. n <= iend_fftp) afft(n) = 0.d0
                   n  = nn*2
                   if(n >= ista_fftp .and. n <= iend_fftp) afft(n) = 0.d0
                end do
             end do
#endif
          end if
       end if
    end if
  end subroutine boundary_zero_into_afft

#ifndef _XC_SAVE_MEMORY_
  subroutine cp_afft_to_cgrad_rho(is,in)
    integer, intent(in) :: is, in
    integer :: i
    do i = ista_fftph, iend_fftph
       cgrad_rho(i,in,is) = -afft(i*2)
    end do
  end subroutine cp_afft_to_cgrad_rho

#endif

  subroutine wd_small_part_of_afft_kt( len_str,str,ne,iloop )
    integer, intent(in)        :: len_str, iloop
    character(len=len_str),intent(in) :: str
    integer, intent(in)        :: ne

    integer                    :: i
!!$      if(mype == 0) then
    write(nfout, *) 'afft (',str,')'
    if(nspin == 2) write(nfout,*) ' #spin = ', iloop
    if(ne*2 > np_fftp) then
       write(nfout,*)  ' all elements'
       write(nfout,'(6f12.6)') (afft(i),i=ista_fftp, iend_fftp)
    else
       write(nfout,*) ' first ',ne,' elements'
!!$            write(nfout,'(6f12.6)') (afft(i),i=1,ne)
       write(nfout,'(6f12.6)') (afft(ista_fftp-1+i),i=1,ne)
       write(nfout,*) ' last ', ne,' elements'
       write(nfout,'(6f12.6)') (afft(iend_fftp-ne+i),i=1,ne)
    endif
!!$         write(nfout,'(" iend_fftp = ",i8)') iend_fftp
!!$         write(nfout,*) ' last ', lx*2,' elements'
!!$         write(nfout,'(" j, k = ", 2i8)') ly_d+ly_d*myrank_cdfft, lz
!!$         n = 2*((ly_d-1)*lx + (lz-1)*lx*ly_d) + ista_fftp - 1
!!$         write(nfout,'(8d12.4)') (afft( n+i),i=1,lx*2)
!!$         write(nfout,'(" n+lx*2 = ",i8)') n+lx*2
!!$            end do
!!$         end do
!!$         n = 0
!!$         do i = ista_fftp, iend_fftp
!!$            if(abs(afft(i)) > 1.d-8) n = n+1
!!$         end do
!!$         write(nfout,'(" number of non-zero term of afft = ",i8)') n
!!$      end if
  end subroutine wd_small_part_of_afft_kt

  subroutine add_sq_afft_to_grad_trho
    integer   :: i
    do i = ista_fftph, iend_fftph     ! MPI
       grad_trho(i) = grad_trho(i) + afft(i*2)**2
    end do
  end subroutine add_sq_afft_to_grad_trho

  subroutine add_sq_afft_to_grad_rho(is)
    integer, intent(in) :: is
    integer :: i
    do i = ista_fftph, iend_fftph     ! MPI
       grad_rho(i,is) = grad_rho(i,is) + afft(i*2)**2
    end do
  end subroutine add_sq_afft_to_grad_rho
! --------------------------------------------

  subroutine Calc_Kernels_ALDA_Rspace(ispin)
    integer, intent(in) :: ispin

    integer :: i, ni, nj
! --------------------------------- start ----------
    if ( xctype == 'ldapw91' ) then
       Call calc_kernel_ex_ldapw91( nspin, ispin, chgrhr_l, ker_ex )
       Call calc_kernel_cr_ldapw91( nspin, ispin, chgrhr_l, ker_cr )
    else if ( xctype == 'ggapbe' ) then
       Call calc_kernel_ex_ggapbe( nspin, ispin, chgrhr_l, grad_rho, f2or1, ker_ex )
       Call calc_kernel_cr_ggapbe( nspin, ispin, chgrhr_l, grad_trho, f2or1, ker_cr )
    endif
! ------------------------- Exc ni shita baao ------
!    if ( xctype == 'ggapbe' ) then
!       Call calc_kt_ex_ggapbe( nspin, ispin, chgrhr_l, grad_rho, f2or1, ker_ex )
!       Call calc_kt_cr_ggapbe( nspin, ispin, chgrhr_l, grad_trho, f2or1, ker_cr )
!    endif
! ---------------------------------------------------
!      ker_ex = 0.0d0; ker_cr = 0.0d0
    Do ni=1, ispin
       Do nj=1, ispin
          if ( ni==nj ) then
             Do i=ista_fftph, iend_fftph
                ker_cr(i,ni,nj) = ker_ex(i,ni) + ker_cr(i,ni,nj)
             End do
          endif
       End do
    End do
! ---------------------- debug ----
!    ker_cr(:,1,1) = chgrhr_l(:,1)     ! atteru.
!    return
! ----------------------- only for experimtal ---
!            see. A. Marini et al, Comp.  Phys. Comm. 180 (2009) 1392.
!
! -----------------------------------------------
    if ( sw_Fxc_enhanced == ON ) ker_cr = ker_cr * 2.0d0
!    if ( sw_Fxc_enhanced ) ker_cr = ker_cr * PAI
  end subroutine Calc_Kernels_ALDA_Rspace

end module m_LinearResponse_ALDA

