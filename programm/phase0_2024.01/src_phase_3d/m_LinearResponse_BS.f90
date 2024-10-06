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
!  This is a module for setting arrays used in the BS-like equation
!  for the LR-TDDFT .
! ======================================================================

! ======================= history ======================================
!  ver 1.0 :  2010/3/31
!               applicable to the crystal such as Si.
!  ver 2.0 :  2011/3/31
!               applicable to the isolated molecule such as C6H6.
!
! ======================================================================

module m_LinearResponse_BS

  use m_Const_Parameters,   only : ON, OFF, DP, CMPLDP, PAI, PAI4, Bohr, &
       &                           TETRAHEDRON, PARABOLIC, DELTA

  use m_LinearResponse_Control,   only :  sw_NODA
  use m_LinearResponse_NonInt,   only : nband_LR, band_start_LR, band_end_LR
  use m_Parallelization,      only : mype, ierr, npes, &
       &                  ista_fftp, iend_fftp, &
       &                  ista_fftph, iend_fftph, &
       &                  np_e, ista_k, iend_k, map_e, map_z, myrank_e, &
       &                  MPI_CommGroup

  use m_Control_Parameters,  only : Num_q_Points, way_of_smearing, neg, &
       &                            ipri, nspin, printable, xctype, kimg

 use m_FFT,                  only : fft_box_size_CD, fft_box_size_CD_c, &
       &                             nfft, nfftp, &
       &                             m_FFT_CD_inverse_c, &
       &                             m_FFT_check_of_negative_CD, &
       &                             m_FFT_CD_direct_c, &
       &                             m_FFT_alloc_CD_box, &
       &                             m_FFT_dealloc_CD_box
  use m_Files,                 only : nfout
  use m_Electronic_Structure,       only : eko_ek
  use m_Crystal_Structure,    only :  rltv, univol
  use m_PlaneWaveBasisSet,    only : ngabc, igf, kg1,kg, igfp_l, kgp

  use m_LinearResponse_tools,  only  : occup_k, occup_kmq, &
       &                                  eko_k, &
       &                                    fsval_k, fsval_kmq, &
       &                                    wfn_k, wfn_kmq, &
       &                                    map_z_k, map_z_kmq, &
       &                                   nbase_k, nbase_kmq, &
       &                                   igf_k, igf_kmq, iba_k, iba_kmq

  use m_Kpoints,                    only : kv3_ek, kv3, vkxyz_ek, qwgt_ek, &
       &                                   k_symmetry, np0, np2
  use m_LinearResponse_ALDA,   only : ker_cr, igfp_LR
  use m_LinearResponse_Density,  only : RhoTilde, occup_lkt_ek

  use m_LinearResponse_Control,     only  : nrd_efermi, nstep, e, nmax_G_LR, &
       &                                    vqxyz,  &
       &                                    xc_kernel_type, RPA, ALDA_R, &
       &                                    scissor, eta, &
       &                                    sw_LongWaveLimit

  use m_LinearResponse_Kernel,       only : Kernel_Coulomb

  use m_Timing,                     only : tstatc0_begin, tstatc0_end
  use m_LinearResponse_tools,  only  : Get_WF_in_Rspace, Get_WF_in_Rspace_CD

  implicit none
  include 'mpif.h'

! ------------------------------------------ BS && ALDA_R --
  Complex(kind=CMPLDP), allocatable :: MatV(:,:), MatK(:,:), MatXi(:,:)
  Complex(kind=CMPLDP), allocatable :: MatL0(:,:,:)

contains

!------------------------------------------------------------------
!!
!!!            Alloc and Dealloc in the case of BS && ALDA_R
!!
!-------------------------------------------------------------------
 subroutine m_LR_alloc_MatK_ALDA
    integer :: matsize
    matsize = nband_LR**2

    if ( mype == 0 ) then
       allocate( MatK(matsize, matsize) ); MatK = 0.0D0
    endif
  end subroutine m_LR_alloc_MatK_ALDA

  subroutine m_LR_dealloc_MatK_ALDA
    if ( mype == 0 ) deallocate( MatK )
  end subroutine m_LR_dealloc_MatK_ALDA

  subroutine m_LR_alloc_MatV
    integer :: matsize
    matsize = nband_LR **2

    if ( mype == 0 ) then
       allocate( MatV(matsize, matsize) ); MatV = 0.0D0
    endif
  end subroutine m_LR_alloc_MatV

  subroutine m_LR_dealloc_MatV
    if ( mype == 0 ) deallocate( MatV )
  end subroutine m_LR_dealloc_MatV

  subroutine m_LR_Alloc_MatXi
    integer :: matsize
    matsize = nband_LR **2

    if ( mype == 0 ) then
       allocate( MatXi(matsize, matsize) ); MatXi = 0.0D0
    endif
  end subroutine m_LR_Alloc_MatXi

  subroutine m_LR_dealloc_MatXi
    if ( mype == 0 ) deallocate( MatXi )
  end subroutine m_LR_dealloc_MatXi

  subroutine m_LR_alloc_MatL0
    if ( mype == 0 ) then
       allocate( MatL0( nband_LR, nband_LR, nspin ) );  MatL0 = 0.0D0
    endif
  end subroutine m_LR_alloc_MatL0

  subroutine m_LR_dealloc_MatL0
    if ( mype == 0 ) deallocate( MatL0 )
  end subroutine m_LR_dealloc_MatL0

!--------------------------------------------------------
!!
!!!           Matrix K, V, Xi      for   BS && ALDA_R
!!
!---------------------------------------------------------
  subroutine Calc_MatK_ALDA
    Real(kind=DP), allocatable :: wf1_tmp(:), wf2_tmp(:), wf3_tmp(:), wf4_tmp(:)
    Real(kind=DP), allocatable ::  chg1r(:), chg1i(:)
    Real(kind=DP), allocatable ::  chg2r(:), chg2i(:)
    Real(kind=DP), allocatable :: Wfns_on_Rspace(:,:)

    integer ib1, ib2, ib3, ib4
    integer ix1, ix2

    integer ik, jk
    integer ngrid
    integer matsize
! ---------------------------------- start ---------------------
    matsize = nband_LR**2
    ngrid = product(fft_box_size_CD(1:3,1))
    ik = 1; jk = 1

    call m_FFT_alloc_CD_box()
    call set_work_arrays
! ----------------------------- main ---------------------
    Do ix1=1, matsize
       ib1 = ( ix1-1 ) /nband_LR +1
       ib2 = mod( ix1-1,nband_LR ) + 1
!
       ib1 = ib1 + band_start_LR - 1;    ib2 = ib2 + band_start_LR - 1

       if ( mype == 0 ) then
          write( nfout,* ) 'Matrix K for ', ix1, '-th element is set'
       endif
       call set_chg1
! ---
       Do ix2=1, matsize
          if ( sw_NODA == ON .and. ix1 /= ix2 ) cycle

          ib3 = ( ix2-1 ) /nband_LR +1
          ib4 = mod( ix2-1,nband_LR ) + 1

          ib3 = ib3 + band_start_LR - 1;      ib4 = ib4 + band_start_LR - 1

          call set_chg2
          call set_each_element_in_MatK
       End do
    End do
! ----------------------------- end ------------
    if ( mype == 0 ) then
       MatK = MatK * univol / dble(ngrid)
!!!!!!!!!!!!!!!!!!!!!!!       MatK = MatK*2.0d0        ! ???? bug ??
    endif
! --
    call unset_work_arrays
    call m_FFT_dealloc_CD_box()

  contains

    subroutine set_work_arrays
      integer jj
      Real*8 c1
      allocate( wf1_tmp(ista_fftp:iend_fftp) ); wf1_tmp = 0.0d0
      allocate( wf2_tmp(ista_fftp:iend_fftp) ); wf2_tmp = 0.0d0
!
      allocate( chg1r(ista_fftph:iend_fftph) ); chg1r = 0.0d0
      allocate( chg1i(ista_fftph:iend_fftph) ); chg1i = 0.0d0
      allocate( chg2r(ista_fftph:iend_fftph) ); chg2r = 0.0d0
      allocate( chg2i(ista_fftph:iend_fftph) ); chg2i = 0.0d0
!
      allocate( Wfns_on_Rspace( ista_fftp:iend_fftp, band_start_LR: band_end_LR ) );
      Wfns_on_Rspace = 0.0d0
      Do ib1 = band_start_LR, band_end_LR
         Call  Get_WF_in_Rspace_CD( ik, ib1, wfn_k, wf1_tmp, &
              &                     nbase_k, iba_k, igfp_LR, map_z_k )
         Wfns_on_Rspace( ista_fftp:iend_fftp, ib1 ) = wf1_tmp( ista_fftp:iend_fftp )
      End do
      Wfns_on_Rspace = Wfns_on_Rspace / sqrt( univol )
! ----
!      c1 = 0.0d0
!      DO jj=ista_fftph,iend_fftph
!         c1 = c1 + Wfns_on_Rspace( 2*jj-1, band_start_LR )**2 &
!              &  + Wfns_on_Rspace( 2*jj,   band_start_LR )**2
!      End do
!      write(*,*) 'c1 = ', c1, c1 / dble(iend_fftph-ista_fftph+1), c1 * univol, c1 * univol /dble(iend_fftph-ista_fftph+1)           ! c1/(iend-isa+1) == 1 deshita.
! ---------------
!             wfns_on_rspace niha sqrt(V) ga kakatteiru
!      stop
! --
    end subroutine set_work_arrays

    subroutine unset_work_arrays
      deallocate( wf1_tmp, wf2_tmp );
      deallocate( chg1r, chg1i ); deallocate( chg2r, chg2i )
      deallocate( Wfns_On_Rspace )
    end subroutine unset_work_arrays

    subroutine set_chg1
      integer ::i, ip

      wf1_tmp( ista_fftp:iend_fftp ) = Wfns_on_Rspace( ista_fftp:iend_fftp, ib1 )
      wf2_tmp( ista_fftp:iend_fftp ) = Wfns_on_Rspace( ista_fftp:iend_fftp, ib2 )
      Do i=ista_fftph, iend_fftph
         ip = 2*i-1
         chg1r(i) = wf1_tmp(ip) *wf2_tmp(ip)   + wf1_tmp(ip+1) *wf2_tmp(ip+1)
         chg1i(i) = wf1_tmp(ip) *wf2_tmp(ip+1) - wf1_tmp(ip+1) *wf2_tmp(ip)
      End do
    end subroutine set_chg1

    subroutine set_chg2
      integer ::i, ip

      wf1_tmp( ista_fftp:iend_fftp ) = Wfns_on_Rspace( ista_fftp:iend_fftp, ib3 )
      wf2_tmp( ista_fftp:iend_fftp ) = Wfns_on_Rspace( ista_fftp:iend_fftp, ib4 )
      Do i=ista_fftph, iend_fftph
         ip = 2*i-1
         chg2r(i) = wf1_tmp(ip) *wf2_tmp(ip)   + wf1_tmp(ip+1) *wf2_tmp(ip+1)
         chg2i(i) = wf1_tmp(ip) *wf2_tmp(ip+1) - wf1_tmp(ip+1) *wf2_tmp(ip)
      End do
    end subroutine set_chg2

    subroutine set_each_element_in_MatK
      integer :: i
      Complex(kind=CMPLDP) :: z1, ztmp, zsum
      Real(kind=DP) :: c1, c2

      ztmp = 0.0D0
      Do i=ista_fftph, iend_fftph
         c1 = chg1r(i)*chg2r(i) + chg1i(i)*chg2i(i)
         c2 = chg1i(i)*chg2r(i) - chg1r(i)*chg2i(i)
         z1 = dcmplx( c1,c2 )
         ztmp = ztmp + z1*ker_cr( i,ik,jk )
      End do
      if ( npes > 1 ) then
         call mpi_allreduce( ztmp, zsum, 2, MPI_DOUBLE_PRECISION, &
              &              MPI_SUM, MPI_CommGroup, ierr )
         ztmp = zsum
      endif
      if ( mype==0 ) MatK(ix1,ix2) = ztmp
    end subroutine set_each_element_in_MatK

  end subroutine Calc_MatK_ALDA

  subroutine Calc_MatV
    Complex(kind=CMPLDP), allocatable :: rho1_tmp(:,:), rho2_tmp(:,:)
    Complex(kind=CMPLDP), allocatable :: rho2_tmp_mpi(:,:)

    integer ib1, ib2, ib3, ib4
    integer ix1, ix2
    integer matsize
!
    real(kind=DP) :: ttr(6)
! ------------------------------------ start- -------
    matsize = nband_LR**2
    call set_work_arrays
    Call getttr (rltv,ttr)
! --
!!    write(*,*) 'nband_LR = ', nband_LR
!!    write(*,*) 'nmax_G_LR = ', nmax_G_LR
! --------------------------- main ----------
    Do ix1=1, matsize
       ib1 = ( ix1-1 ) /nband_LR +1
       ib2 = mod( ix1-1,nband_LR ) + 1

       ib1 = ib1 + band_start_LR - 1;     ib2 = ib2 + band_start_LR - 1

       if ( mype == 0 ) then
          write( nfout,* ) 'Matrix V for ', ix1, '-th element is set'
       endif

       call set_rho1_tmp

       Do ix2=1, matsize
          if ( sw_NODA == ON .and. ix1 /= ix2 ) cycle

          ib3 = ( ix2-1 ) /nband_LR +1
          ib4 = mod( ix2-1,nband_LR ) + 1

          ib3 = ib3 + band_start_LR - 1;     ib4 = ib4 + band_start_LR - 1

          call set_rho2_tmp
          call set_each_element_in_MatV
       End do
    End do
! -------------------------- end -----------------
    if ( mype== 0 ) then
       MatV = MatV / univol
    endif
    call unset_work_arrays
  contains

    subroutine set_work_arrays
      allocate( rho1_tmp( kg1, kv3_ek ) ); rho1_tmp = 0.0d0
      allocate( rho2_tmp( kg1, kv3_ek ) ); rho2_tmp = 0.0d0
    end subroutine set_work_arrays

    subroutine unset_work_arrays
      deallocate( rho1_tmp, rho2_tmp )
    end subroutine unset_work_arrays

    subroutine set_rho1_tmp
      integer :: jb2

      if ( map_e(ib2) == myrank_e ) then
         jb2 = map_z(ib2)
         rho1_tmp( 1:kg1, 1:kv3_ek ) = RhoTilde( 1:kg1, ib1, jb2, 1:kv3_ek )
      endif
      if ( npes > 1 ) then
         call mpi_bcast( rho1_tmp, 2*kg1*kv3_ek, MPI_DOUBLE_PRECISION, &
              &         map_e(ib2), MPI_CommGroup, ierr )
      endif
    end subroutine set_rho1_tmp

    subroutine set_rho2_tmp
      integer :: jb4

      if ( map_e(ib4) == myrank_e ) then
         jb4 = map_z(ib4)
         rho2_tmp( 1:kg1, 1:kv3_ek ) = RhoTilde( 1:kg1, ib3, jb4, 1:kv3_ek )
      endif
      if ( npes > 1 ) then
         call mpi_bcast( rho2_tmp, 2*kg1*kv3_ek, MPI_DOUBLE_PRECISION, &
              &         map_e(ib4), MPI_CommGroup, ierr )
      endif
    end subroutine set_rho2_tmp

    subroutine set_each_element_in_MatV
      integer     :: i, ig
      Real(kind=DP) :: ga, gb, gc, g2
      Complex(kind=CMPLDP) :: z1, z2, ztmp
! -
      if ( mype /= 0 ) return
!
      ztmp = 0.0d0
      Do i=2, nmax_G_LR
         z1 = rho1_tmp( i,1 );  z2 = rho2_tmp( i,1 )
         ztmp = ztmp + z1* conjg(z2) * Kernel_Coulomb(i)
      End Do
      MatV( ix1,ix2 ) = ztmp
    end subroutine set_each_element_in_MatV

  end subroutine Calc_MatV

  subroutine Calc_MatXi
    if ( mype == 0 )  then
       select case( xc_kernel_type )
       case( RPA )
          MatXi = -MatV
       case ( ALDA_R )
          MatXi = -MatK - MatV
       end select
    endif
  end subroutine Calc_MatXi

!--------------------------------------------------------
!!
!!!        Matrix L0    for   BS && ALDA_R,   LongWaveLimit
!!
!---------------------------------------------------------
  subroutine Calc_MatL0( e_value )
    real(kind=DP), intent(in) :: e_value

    real(kind=DP)         :: occ1, occ2
    complex(kind=CMPLDP)  :: omega1, omega2

    integer               :: i,j, k, ispin
    integer               :: ik1, ik2, ik3, ib1, ib2
    integer               :: k_max, nq_max
    integer :: id_sname = -1
! ------------------------------- start ------------
    if ( mype /= 0 ) return
    nq_max = 0
    call set_value_kmax
    omega1 = cmplx( e_value, eta );  omega2 = cmplx( e_value, -eta );
!
    MatL0 = 0.0d0
! -------------------------------- main -------------
    Do ispin=1, nspin
       Do k=1, k_max
          ik1 = ( nq_max +1 )*nspin*(k-1) + ispin
          ik2 = ik1 + nspin
          ik3 =  nspin*(k-1) + ispin
          !
          Do ib1=band_start_LR, band_end_LR           ! unoccpuied state
             occ1 = Occup_lkt_ek( ib1, ik1 );

             Do ib2=band_start_LR, band_end_LR          ! occupied state
                occ2 = Occup_lkt_ek( ib2, ik1 )
                call set_each_element_in_MatL0
             End do
          End do
       End do
    End do
    if ( nspin==1 ) MatL0 = MatL0 *2.0d0
! --------------------------------- end -------------
    call tstatc0_end(id_sname)

  contains

    subroutine set_value_kmax
      if ( way_of_smearing == TETRAHEDRON ) then
         k_max = np2
      else
         k_max = kv3_ek / ( nspin*( nq_max + 1 ))
      endif
    end subroutine set_value_kmax

    subroutine set_each_element_in_MatL0
      integer ::  itmp1, itmp2
      Real(kind=DP) :: ebi, ebj, ediff, c1
      Complex(kind=CMPLDP) :: zd1, zd2, z1, z2, gfn

      ebi = Eko_ek( ib1,ik1 );  ebj = Eko_ek( ib2,ik1 )
      ediff = ebi - ebj + scissor

      zd1 = omega1 - ediff;  zd2 = omega2 + ediff
      gfn = 1.0d0 / zd1 - 1.0d0 / zd2

      itmp1 = ib1 - band_start_LR +1; itmp2 = ib2 - band_start_LR +1

      c1 = ( 1.0D0 - occ1 ) *occ2
      MatL0( itmp1,itmp2,ispin ) = c1*gfn *cmplx( 0.0D0, 1.0D0 )
    end subroutine set_each_element_in_MatL0

  end subroutine Calc_MatL0

!!!!!!!!!!!!!!! DEBUG -
  subroutine Check_Charge_on_mesh
    Real(kind=DP), allocatable :: wf1_tmp(:), wf2_tmp(:), wf3_tmp(:), wf4_tmp(:)
    Real(kind=DP), allocatable :: Wfns_on_Rspace(:,:)

    Real(kind=DP) :: c1
    integer :: ik, jk
    real(kind=DP)         :: occ1
    complex(kind=CMPLDP)  :: ztmp, z1

    integer :: ngrid, nfftph
    integer               :: i, ib1
    integer :: id_sname = -1

    ngrid = product(fft_box_size_CD(1:3,1))
    nfftph = nfftp / 2
    ik = 1; jk = 1

    call m_FFT_alloc_CD_box()
    call set_work_arrays
! ----------------------------- main ---------------------
    Do i=1, nfftph
       ztmp = 0.0d0
       Do ib1=1, neg
          occ1 = Occup_lkt_ek( ib1, ik );
          if ( occ1 < Delta ) cycle
          z1 = cmplx( Wfns_On_Rspace( 2*i-1,ib1 ), Wfns_On_Rspace( 2*i,ib1 ) )
          ztmp = ztmp + occ1 *z1 *conjg(z1)
       End do
!
       if ( ker_cr(i,1,1) > 1.0E-4 ) then
          c1 = 2.0 *real(ztmp) / ker_cr(i,1,1)
          write(*,*) i, ker_cr(i,1,1), c1
       endif
    End do
! --
    call unset_work_arrays
    call m_FFT_dealloc_CD_box()

  contains

    subroutine set_work_arrays
      allocate( wf1_tmp(ista_fftp:iend_fftp) ); wf1_tmp = 0.0d0
      allocate( Wfns_on_Rspace( ista_fftp:iend_fftp, 1:neg ) );
      Wfns_on_Rspace = 0.0d0
      Do ib1 = 1, neg
         Call  Get_WF_in_Rspace_CD( ik, ib1, wfn_k, wf1_tmp, &
              &                     nbase_k, iba_k, igfp_LR, map_z_k )
         Wfns_on_Rspace( ista_fftp:iend_fftp, ib1 ) = wf1_tmp( ista_fftp:iend_fftp )
      End do
      Wfns_on_Rspace = Wfns_on_Rspace / sqrt( univol )
    end subroutine set_work_arrays

    subroutine unset_work_arrays
      deallocate( wf1_tmp )
      deallocate( Wfns_On_Rspace )
    end subroutine unset_work_arrays

  End subroutine Check_Charge_on_mesh

  subroutine Check_Exc_energy
    Real(kind=DP), allocatable :: wf1_tmp(:), wf2_tmp(:), wf3_tmp(:), wf4_tmp(:)
    Real(kind=DP), allocatable :: Wfns_on_Rspace(:,:)

    Real(kind=DP) :: c1, weight
    integer :: ik, jk
    real(kind=DP)         :: occ1
    complex(kind=CMPLDP)  :: ztmp, z1

    integer :: ngrid, nfftph
    integer               :: i, ib1
    integer :: id_sname = -1

    ngrid = product(fft_box_size_CD(1:3,1))
    nfftph = nfftp / 2
    ik = 1; jk = 1

    call m_FFT_alloc_CD_box()
    call set_work_arrays
! ----------------------------- main ---------------------
    ztmp = 0.0d0
    Do i=1, nfftph
       Do ib1=1, neg
          occ1 = Occup_lkt_ek( ib1, ik );
          if ( occ1 < Delta ) cycle
          z1 = cmplx( Wfns_On_Rspace( 2*i-1,ib1 ), Wfns_On_Rspace( 2*i,ib1 ) )
          weight = ker_cr( i,1,1 )
          ztmp = ztmp + occ1 *z1 *conjg(z1) *weight
       End do
    End do
    ! --
    write(*,*) 'Exc = ', ztmp *univol / dble(ngrid) * 2.0d0            ! OK
! ----------------------------------------------- 2.0 : spin degeneracy
! --
    call unset_work_arrays
    call m_FFT_dealloc_CD_box()

  contains

    subroutine set_work_arrays
      allocate( wf1_tmp(ista_fftp:iend_fftp) ); wf1_tmp = 0.0d0
      allocate( Wfns_on_Rspace( ista_fftp:iend_fftp, 1:neg ) );
      Wfns_on_Rspace = 0.0d0
      Do ib1 = 1, neg
         Call  Get_WF_in_Rspace_CD( ik, ib1, wfn_k, wf1_tmp, &
              &                     nbase_k, iba_k, igfp_LR, map_z_k )
         Wfns_on_Rspace( ista_fftp:iend_fftp, ib1 ) = wf1_tmp( ista_fftp:iend_fftp )
      End do
      Wfns_on_Rspace = Wfns_on_Rspace / sqrt( univol )
    end subroutine set_work_arrays

    subroutine unset_work_arrays
      deallocate( wf1_tmp )
      deallocate( Wfns_On_Rspace )
    end subroutine unset_work_arrays

  End subroutine Check_Exc_energy

! --
  subroutine Print_MatK
    integer ib1, ib2, ib3, ib4
    integer ix1, ix2
    integer matsize
!
    matsize = nband_LR**2
! --------------------------- main ----------
    Do ix1=1, matsize
       ib1 = ( ix1-1 ) /nband_LR +1
       ib2 = mod( ix1-1,nband_LR ) + 1

       ib1 = ib1 + band_start_LR - 1;     ib2 = ib2 + band_start_LR - 1

!!!!       if ( ib1 == band_start_LR .and. ib2 == band_end_LR ) then
          if ( mype == 0 ) then
             write(*,*) 'K element = ', ib1, ib2, MatK( ix1,ix1 )
          endif
!!!       endif
    End do
   if ( npes > 1 ) call mpi_barrier( MPI_CommGroup,ierr )
  end subroutine Print_MatK

  subroutine Print_MatV
    integer ib1, ib2, ib3, ib4
    integer ix1, ix2
    integer matsize
!
    matsize = nband_LR**2
! --------------------------- main ----------
    Do ix1=1, matsize
       ib1 = ( ix1-1 ) /nband_LR +1
       ib2 = mod( ix1-1,nband_LR ) + 1

       ib1 = ib1 + band_start_LR - 1;     ib2 = ib2 + band_start_LR - 1

!!       if ( ib1 == band_start_LR .and. ib2 == band_end_LR ) then
          if ( mype == 0 ) then
             write(*,*) 'V element = ', ib1, ib2, MatV( ix1,ix1 )
          endif
!!!!       endif
    End do
    if ( npes > 1 ) call mpi_barrier( MPI_CommGroup,ierr )
  end subroutine Print_MatV
! -----------------
end module m_LinearResponse_BS
