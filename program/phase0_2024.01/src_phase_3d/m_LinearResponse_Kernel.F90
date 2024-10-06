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
!  This is a module for setting Coulomb and exchange kernels for the LR-TDDFT .
! ======================================================================

! ======================= history ======================================
!  ver 1.0 :  2010/3/31
!               applicable to the crystal such as Si.
!  ver 2.0 :  2011/3/31
!               applicable to the isolated molecule such as C6H6.
!
! ======================================================================

module m_LinearResponse_Kernel

  use m_Charge_Density,        only : chgq_l
  use m_Const_Parameters,        only : Bohr
  use m_PseudoPotential,      only : itpcc,ilmt,ltp,taup,il2p,isph,iqitg, &
      & rhpcg_l, qitg_l,rhpcg_diff_l,qitg_diff_l,dl2p, &
      & m_PP_include_vanderbilt_pot, m_PP_find_maximum_l,nlmt

  use m_Ionic_System,         only : ntyp,natm,iwei,ityp,pos,zfm3_l
  use m_Const_Parameters,     only : DP, CMPLDP, BUCS, PAI, PAI2, PAI4, &
       &                             GGA, LDA, DOWN,UP, &
       &                             Partial_Core_Charge,Valence_plus_PC_Charge

  use m_PlaneWaveBasisSet,    only : ngabc, igf, kg1,kg, igfp_l, n_rGpv, kgp
  use m_Crystal_Structure,    only :  rltv, univol

  use m_Control_Parameters,   only : ipri, nspin, kimg,af, neg, printable, &
       & xctype
  use m_Files,                only : nfout

  use m_Timing,               only : tstatc0_begin, tstatc0_end
  use m_Parallelization,      only : ista_kngp,iend_kngp,npes,mype, ierr, &
       & npes_cdfft, nrank_ggacmp, myrank_cdfft, &
       & myrank_ggacmp, ista_fftp, iend_fftp, ista_fftph, iend_fftp, &
       & nis_fftp, nie_fftp, nel_fftp, idisp_fftp,np_fftp,mp_fftp,  &
       & nel_fftph, idisp_fftph, nrest_cdfft, &
       & mpi_ggacmp_cross_world,mpi_cdfft_world, &
       & map_ggacmp, map_pe2ggacmp, map_pe2cdfft, MPI_CommGroup, &
       & ista_fftph, iend_fftph

  use m_FFT,                  only : fft_box_size_CD, fft_box_size_CD_c, &
       &                             nfftp, &
       &                             m_FFT_CD_inverse_c, &
       &                             m_FFT_check_of_negative_CD, &
       &                             m_FFT_CD_direct_c

#ifdef _MPIFFT_
  use m_FFT,                  only : m_FFT_set_cdata  &
       &                           , lx,ly,lz,ly_d,lz_d, ny_d,nz_d
  use m_PlaneWaveBasisSet,    only : igfp_l_c
#else
  use m_FFT,                  only : m_FFT_alloc_CD_box &
       &                           , m_FFT_dealloc_CD_box
#endif

  use m_FFT,                   only : nfft, fft_box_size_WF
  use m_XC_Potential,         only : check_of_xctype
! ----------------
  use m_LinearResponse_Control,  only : nmax_G_LR, vqxyz, xc_kernel_type, &
       &                                RPA, ALDA_G, LRC, LRC_alpha
  use m_LinearResponse_ALDA,    only :  ker_ex, ker_cr, set_ispin_kt, &
       &                                xc_deallocate, xc_allocate, &
       &                                m_LR_alloc_phase_factor, &
       &                                m_LR_dealloc_phase_factor, &
       &                                calc_phase_factor, afft_allgatherv, &
       &                                Calc_Density_in_Rspace, &
       &                                m_LR_alloc_RKernels_ALDA, &
       &                                m_LR_dealloc_RKernels_ALDA, afft, &
       &                                igfp_LR, &
       &                                Calc_Kernels_ALDA_Rspace

  use m_LinearResponse_Control,  only : rcut_yukawa

!
  Implicit None
  include 'mpif.h'

! ------------------ public ---------
  real(kind=DP), allocatable :: Kernel_Coulomb(:)
  complex(kind=CMPLDP), allocatable :: Kernel_XC(:,:,:,:)

contains

!------------------------------------------------------------------
!!
!!!            Alloc and Dealloc Kernels
!!
!-------------------------------------------------------------------
  subroutine m_LR_alloc_Coulomb_Kernels
    if ( mype == 0 ) then
       allocate( Kernel_Coulomb( nmax_G_LR ) ); Kernel_Coulomb = 0.0D0
    endif
  end subroutine m_LR_alloc_Coulomb_Kernels

  subroutine m_LR_dealloc_Coulomb_Kernels
    if ( mype == 0 ) deallocate( Kernel_Coulomb )
  end subroutine m_LR_dealloc_Coulomb_Kernels

  subroutine m_LR_alloc_XC_Kernels
    if ( mype == 0 ) then
       allocate( Kernel_xc( nmax_G_LR,nmax_G_LR,nspin,nspin ) ); Kernel_XC = 0.0D0
    endif
  end subroutine m_LR_alloc_XC_Kernels

  subroutine m_LR_dealloc_XC_Kernels
    if ( mype == 0 ) deallocate( Kernel_xc )
  end subroutine m_LR_dealloc_XC_Kernels

!-----------------------------------------------------------
!!
!!!              Coulomb Kernels
!!
!-----------------------------------------------------------
  subroutine set_Coulomb_Kernel
    integer i
    real(kind=DP) :: ga, gb, gc, g2
    real(kind=DP) :: ttr(6)
! ----------------------------- start -------------------
    if ( mype /= 0 ) return
    Call getttr (rltv,ttr)
! ------------------------------ main ------------------
    Do i=1, nmax_G_LR
       ga = ngabc(i,1) + vqxyz( 1,1,BUCS )
       gb = ngabc(i,2) + vqxyz( 1,2,BUCS )
       gc = ngabc(i,3) + vqxyz( 1,3,BUCS )
       g2 = ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
            &   +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga
       Kernel_Coulomb(i) = PAI4 / g2
    End Do
  end subroutine set_Coulomb_Kernel

  subroutine set_Screened_Coulomb_Kernel
    integer i
    real(kind=DP) :: ga, gb, gc, g2, gm
    real(kind=DP) :: ttr(6)
! ----------------------------- start -------------------
    if ( mype /= 0 ) return
    Call getttr (rltv,ttr)
    gm = 1.0d0 / rcut_yukawa
    gm = gm **2
! ------------------------------ main ------------------
    Do i=1, nmax_G_LR
       ga = ngabc(i,1) + vqxyz( 1,1,BUCS )
       gb = ngabc(i,2) + vqxyz( 1,2,BUCS )
       gc = ngabc(i,3) + vqxyz( 1,3,BUCS )
       g2 = ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
            &   +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga
       g2 = g2 + gm
       Kernel_Coulomb(i) = PAI4 / g2
    End Do
  end subroutine set_Screened_Coulomb_Kernel

! ---------------------------------------------------------
!!
!!!          XC Kernels ( RPA, LRC, ALDA )
!!
! ---------------------------------------------------------
  subroutine set_XC_Kernel_RPA
    if ( mype == 0 )  Kernel_XC = 0.0d0
  end subroutine set_XC_Kernel_RPA

  subroutine set_XC_Kernel_LRC( alpha_LRC )
    Real(kind=DP), intent(in) :: alpha_LRC
!
    integer i, ni, nj
    real(kind=DP) :: ga, gb, gc, g2
    real(kind=DP) :: ttr(6)
! ---------------------------------- start -------
    if ( mype /= 0 ) return
    Call getttr (rltv,ttr)
! ---------------------------------- main --------
    Do i=1, nmax_G_LR
       ga = ngabc(i,1) + vqxyz( 1,1,BUCS )
       gb = ngabc(i,2) + vqxyz( 1,2,BUCS )
       gc = ngabc(i,3) + vqxyz( 1,3,BUCS )
       g2 = ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
               &   +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga
       Do ni=1, nspin
          Do nj=1, nspin
             Kernel_XC(i,i,ni,nj) = -alpha_LRC / g2
          End Do
       End Do
    End Do
  end subroutine set_XC_Kernel_LRC

  subroutine set_XC_Kernel_LRC_Dynamic( alpha_LRC, beta_LRC, ene )
    Real(kind=DP), intent(in) :: alpha_LRC, beta_LRC, ene
!
    integer i, ni, nj
    real(kind=DP) :: ga, gb, gc, g2, c1
    real(kind=DP) :: ttr(6)
! ------------------------------- start --------------
    if ( mype /= 0 ) return
    Call getttr (rltv,ttr)
! ------------------------------ main --------------
    c1 = alpha_LRC + beta_LRC * ene**2
!
    Do i=1, nmax_G_LR
       ga = ngabc(i,1) + vqxyz( 1,1,BUCS )
       gb = ngabc(i,2) + vqxyz( 1,2,BUCS )
       gc = ngabc(i,3) + vqxyz( 1,3,BUCS )
       g2 = ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
            &   +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga
       Do ni=1, nspin
          Do nj=1, nspin
             Kernel_XC(i,i,ni,nj) = -c1 / g2
          End do
       End do
    End Do
  end subroutine set_XC_Kernel_LRC_Dynamic

  subroutine set_XC_kernel_ALDA_G( input_Charge )
    integer, intent(in) :: input_charge
!
    Real(kind=DP), allocatable :: afft_mpi1(:)
    integer, allocatable       :: GVec_Table( :,:,: )

    integer ngrid, iloop
    integer ispin
! --------------------------- start ---------------
    ngrid = product( fft_box_size_CD(1:3,1) )
! --------------------------- main ---------------
    call set_ispin_kt( ispin, input_charge )
    if(ipri >= 2) write(nfout,*) ' ! ispin = ', ispin
    call xc_allocate( ispin,nfout )
!
    call m_LR_alloc_Phase_Factor
    Call Calc_Phase_Factor

    Call Calc_Density_in_Rspace(ispin, input_Charge )
    call Calc_Kernels_ALDA_Rspace(ispin)

    call Alloc_Tables_For_Gspace
    call Set_Tables_For_Gspace
    call Map_Kernels_To_Gspace
    call Dealloc_Tables_For_Gspace

    call m_LR_dealloc_Phase_Factor
    call xc_deallocate

  contains

    subroutine Alloc_Tables_For_Gspace
      integer :: i,j,k
      i = n_rGpv(1); j = n_rGpv(2); k = n_rGpv(3)
      allocate( GVec_Table(-i:i,-j:j,-k:k) );  GVec_Table = 0
    end subroutine Alloc_Tables_For_Gspace

    subroutine Dealloc_Tables_For_Gspace
      Deallocate( GVec_Table )
    end subroutine Dealloc_Tables_For_Gspace

    subroutine Set_Tables_For_Gspace
      integer :: i, igf1, igf2, igf3
      do i = 1, kgp
         igf1 = ngabc(i,1); igf2 = ngabc(i,2); igf3 = ngabc(i,3)
         if(        igf1 >= -n_rGpv(1) .and. igf1 <= n_rGpv(1) &
            & .and. igf2 >= -n_rGpv(2) .and. igf2 <= n_rGpv(2) &
            & .and. igf3 >= -n_rGpv(3) .and. igf3 <= n_rGpv(3)) then
            GVec_Table(igf1,igf2,igf3) = i
         end if
      end do
    end subroutine Set_Tables_For_Gspace

    subroutine Search_Tables_For_Gspace( ix, iy, iz, G_found )
      integer, intent(in)  :: ix, iy, iz
      integer, intent(out) :: G_found

      G_found = -1
      if ( ix < -n_rGpv(1) .or. ix > n_rGpv(1) ) return
      if ( iy < -n_rGpv(2) .or. iy > n_rGpv(2) ) return
      if ( iz < -n_rGpv(3) .or. iz > n_rGpv(3) ) return
      G_found = GVec_Table( ix, iy, iz )
    end subroutine Search_Tables_For_Gspace

    subroutine Map_Kernels_To_Gspace
      integer ntmp_x, ntmp_y, ntmp_z, G_to_be_mapped
      integer i, j, ni, nj
      integer ip, i1
      real(kind=DP) :: normalization_factor
! ----
      if ( mype == 0 ) Kernel_XC = 0.0d0
      allocate( afft_mpi1(nfftp) )
      call m_FFT_alloc_CD_box()
      normalization_factor = 1.0d0/ product(fft_box_size_CD(1:3,1))
! ----------
      Do ni=1, ispin
         Do nj=1, ispin
            afft = 0.0d0
            Do i=ista_fftph, iend_fftph
               ip = 2 *i - 1
               afft(ip)   = ker_cr(i,ni,nj);     afft(ip+1) = 0.0d0
            End do
            ! ------------- Map to Kernel_XC --
            call m_FFT_CD_direct_c( nfout,afft )     ! afft(R sp) -> afft(G sp)
            ! -------
            call afft_allgatherv( afft, afft_mpi1 )
            afft_mpi1 = afft_mpi1 * normalization_factor
            !
            Do i=1, nmax_G_LR
               Do j=1, nmax_G_LR
!                  ntmp_x = ngabc(i,1) -ngabc(j,1);
!                  ntmp_y = ngabc(i,2) -ngabc(j,2);
!                  ntmp_z = ngabc(i,3) -ngabc(j,3);
                  ntmp_x = ngabc(j,1) -ngabc(i,1);
                  ntmp_y = ngabc(j,2) -ngabc(i,2);
                  ntmp_z = ngabc(j,3) -ngabc(i,3);

                  call Search_Tables_For_Gspace( ntmp_x, ntmp_y, ntmp_z, &
                       &                         G_to_be_mapped  )
                  if ( G_to_be_mapped <= 0 ) then
                     write(*,*) 'Error Code ', i,j, G_to_be_mapped
                     goto 1000
                  endif
                  i1 = igfp_LR( G_to_be_mapped )
                  ip = 2 * i1 -1
! --------------------
                  if ( mype == 0 ) then
                     Kernel_XC(i,j,ni,nj) = dcmplx( afft_mpi1(ip), afft_mpi1(ip+1) ) &
                          &         / univol
                  endif
1000              continue

               End do
            End do
         End do
      End do
! -------------------------------
      call m_FFT_dealloc_CD_box()
      deallocate( afft_mpi1 )
    end subroutine Map_Kernels_To_Gspace

  end subroutine set_XC_kernel_ALDA_G

  subroutine set_XC_Kernel_ALDA_R( input_charge )
    integer, intent(in) :: input_charge

    integer ngrid, iloop
    integer ispin
    Real(kind=DP), allocatable :: afft_mpi1(:)

    call set_ispin_kt(ispin, input_Charge )
    if (ipri >= 2) write(nfout,*) ' ! ispin = ', ispin
    call xc_allocate( ispin,nfout )
    Call Calc_Density_in_Rspace( ispin, input_Charge )
    call Calc_Kernels_ALDA_Rspace( ispin )
    call xc_deallocate
  end subroutine set_XC_Kernel_ALDA_R

end module m_LinearResponse_Kernel

