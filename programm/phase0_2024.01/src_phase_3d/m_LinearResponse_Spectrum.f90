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
!  This is a module for calculating and printing spectra.
! ======================================================================

! ======================= history ======================================
!  ver 1.0 :  2010/3/31
!               applicable to the crystal such as Si.
!  ver 2.0 :  2011/3/31
!               applicable to the isolated molecule such as C6H6.
!
! ======================================================================

module m_LinearResponse_Spectrum

  use m_Const_Parameters,           only : DP, ON, OFF, CMPLDP, DELTA, &
       &                                   BUCS, Hartree, PAI4, speed_of_light, AU_VELOCITY
  use m_Control_Parameters,         only : nspin, neg
  use m_Timing,                     only : tstatc0_begin, tstatc0_end
!
!  use m_Parallelization,        only : mype, myrank_e, npes, MPI_CommGroup, &
!       &                               ierr, map_e, map_z
  use m_Parallelization
  use m_Files,                      only  : nfout, nf_LR_spectra
  use m_Crystal_Structure,   only : univol, rltv
  use m_Electronic_Structure,    only : eko_ek
  use m_LinearResponse_Control,     only  : e, nmax_G_LR, &
       &                                    vqxyz,  &
       &                                    xc_kernel_type, ALDA_R,  &
       &                                    sw_NLF, &
       &                                    spectrum_type, &
       &                                    OPTICS, PACS, EELS, IXSS, &
       &                                    e_high
  use m_LinearResponse_Density,     only :  RhoTilde, occup_lkt_ek

!  use m_IterationNumbers,     only : nk_in_the_process, nk_converged
!  use m_PlaneWaveBasisSet,    only : kg1

  use m_LinearResponse_Kernel, only : Kernel_Coulomb, Kernel_XC


  use m_LinearResponse_NonInt,  only : Mat_Chi0, &
       &                               nband_LR, band_start_LR, band_end_LR


  use m_LinearResponse_ALDA, only : MatrixT
  use m_PlaneWaveBasisSet,    only : ngabc
  use m_LinearResponse_BS,   only : MatL0, MatXi

  Implicit None
!  include 'mpif.h'

  Complex(kind=CMPLDP), allocatable :: Chi_Head(:,:)
  Complex(kind=CMPLDP), allocatable :: MatL(:,:)

contains

!-------------------------------------------------------------
!!
!!!               Alloc /Dealloc Chi_Head
!!
! ------------------------------------------------------------
  subroutine m_LR_alloc_Chi_Head
    if ( mype == 0 ) then
       allocate( Chi_Head( nspin,nspin )  ); Chi_Head = 0.0d0
    endif
  end subroutine m_LR_alloc_Chi_Head

  subroutine m_LR_dealloc_Chi_Head
    if ( mype == 0 ) deallocate( Chi_Head )
  end subroutine m_LR_dealloc_Chi_Head

!-------------------------------------------------------------
!!
!!!               Alloc /Dealloc Matrix L of BS equation
!!
! ------------------------------------------------------------
  subroutine m_LR_Alloc_MatL
    integer :: matsize
    matsize = nband_LR**2
    if ( mype == 0 ) then
       allocate( MatL( matsize, matsize) ); MatL = 0.0D0
    endif
  end subroutine m_LR_Alloc_MatL

  subroutine m_LR_dealloc_MatL
    if ( mype == 0 ) deallocate( MatL )
  end subroutine m_LR_dealloc_MatL
!---------------------------------------------------------
!!
!!!          Solver for Dyson-like Equation
!!
!----------------------------------------------------------
  subroutine Calc_Chi_Head( ene )
    real(kind=DP), intent(in)  :: ene

    Complex(kind=CMPLDP), allocatable :: MatA( :,: ), TmpVector(:,:)
    integer :: matsize, i, j
    integer :: id_sname = -1
! ---------------------------- start --------
    if ( mype /= 0 ) return

    call tstatc0_begin('Calc_Chi_Head ', id_sname)
    call Alloc_Arrays_For_LinEq

    if ( xc_kernel_type == ALDA_R ) then
       call set_XC_contrib_GatOnida
    else
       call set_XC_contrib_default
    endif

    if ( sw_NLF == OFF ) call add_LF_effects
    if ( spectrum_type /= OPTICS ) call add_hartree_contrib_G0

    call Add_Diagonal_Elements
    call Initial_Guess_VecX
    Call Calc_LinearEq( nspin, matsize, MatA, TmpVector )
    Call Extract_ChiHead
    call Dealloc_Arrays_For_LinEq
!
!    call tstatc0_end(id_sname)

  contains

    subroutine Extract_ChiHead
      Chi_Head(1:nspin,1:nspin) = TmpVector(1:nspin,1:nspin)
    end subroutine Extract_ChiHead

    subroutine Alloc_Arrays_For_LinEq
      matsize = nspin *nmax_G_LR
      Allocate( MatA( matsize, matsize ) );   MatA = 0.0d0
      Allocate( TmpVector( matsize,nspin ) ); TmpVector = 0.0d0
    end subroutine Alloc_Arrays_For_LinEq

    subroutine Dealloc_Arrays_For_LinEq
      Deallocate( MatA, TmpVector )
    end subroutine Dealloc_Arrays_For_LinEq

    subroutine Set_XC_Contrib_GatOnida
      integer :: i, j, ni, nj, tmp_i, tmp_j
      Do i=1, nmax_G_LR
         Do j=1, nmax_G_LR
            Do ni=1, nspin
               Do nj=1, nspin
                  tmp_i = ( i-1 )*nspin + ni
                  tmp_j = ( j-1 )*nspin + nj
                  MatA( tmp_i, tmp_j ) = -MatrixT( i,j,ni,nj )
               End do
!               write(*,*) 'Y', i,j,ni,nj, MatA( tmp_i, tmp_j )
            End Do
         End Do
      End Do
!      stop
    end subroutine Set_XC_Contrib_GatOnida

    subroutine set_XC_contrib_default
      integer       :: i, j, k, ni, nj, nk
      integer       :: tmp_i, tmp_j, tmp_k
      complex(kind=CMPLDP) :: d1

      Do i=1, nmax_G_LR
         Do j=1, nmax_G_LR
            Do ni=1, nspin
               Do nj=1, nspin
!
                  d1 = 0.0d0
                  Do k=1, nmax_G_LR
                     Do nk=1, nspin

                        if ( ni == nk ) then
                           tmp_i = (i-1)*nspin + ni
                           tmp_j = (j-1)*nspin + nj
                           tmp_k = (k-1)*nspin + nk
                           d1 = d1 + Mat_Chi0( i,k,ni ) *Kernel_XC( k,j,nk,nj )
                        endif
                     End do
                  End do
                  MatA( tmp_i,tmp_j ) = -d1
!                  write(*,*) 'Y', i,j,ni,nj, MatA( tmp_i, tmp_j )
               End do
            End do
         End do
      End do
!      stop
    end subroutine set_XC_contrib_default

    subroutine add_LF_effects
      integer       :: i, j, ni, nj
      integer       :: tmp_i, tmp_j

      Do i=1, nmax_G_LR
         Do j=2, nmax_G_LR
!            if ( i /= j ) cycle

            Do ni=1, nspin
               Do nj=1, nspin
                  tmp_i = (i-1)*nspin + ni
                  tmp_j = (j-1)*nspin + nj

                  MatA( tmp_i,tmp_j ) = MatA( tmp_i, tmp_j ) &
                       &         - Mat_Chi0( i,j,ni ) *Kernel_Coulomb(j)

!                 write(*,*) tmp_i, tmp_j, Chi0_Full(i,j,ni)
!                 write(*,*) tmp_i, tmp_j, Kernel_Coulomb(j)
               End do
            End do
         End do
      End do
!      stop
    end subroutine add_LF_effects

    subroutine add_hartree_contrib_G0
      integer :: i, j, ni, nj, tmp_i, tmp_j
      Do i=1, nmax_G_LR
         Do j=1, 1
!            if ( i /= j ) cycle

            Do ni=1, nspin
               Do nj=1, nspin
                  tmp_i = (i-1)*nspin + ni
                  tmp_j = (j-1)*nspin + nj
                  MatA( tmp_i,tmp_j ) = MatA( tmp_i,tmp_j ) &
                       &              - Mat_Chi0( i,j,ni ) *Kernel_Coulomb(j)
               End do
            End do
         End do
      End do
    end subroutine add_hartree_contrib_G0

    subroutine Add_Diagonal_Elements
      integer  :: i, ni, tmp_i
      Do i=1, nmax_G_LR
         Do ni=1, nspin
            tmp_i = (i-1)*nspin + ni
            MatA( tmp_i,tmp_i ) = MatA( tmp_i,tmp_i )  +1.0d0
!            write(*,*) 'F ',tmp_i, MatA( tmp_i,tmp_i)
         End do
      End do
!      stop
    end subroutine Add_Diagonal_Elements

    subroutine Initial_Guess_VecX
      integer :: i, j, ni, nj, tmp_i, tmp_j

      TmpVector = 0.0d0
      Do j=1, 1
         Do nj=1, nspin
            Do i=1, nmax_G_LR
               Do ni=1, nspin
                  if ( ni==nj ) then
                     tmp_i = (i-1)*nspin + ni
                     tmp_j = (j-1)*nspin + nj
                     TmpVector(tmp_i,nj) = Mat_Chi0( i,j,ni )
                  endif
               End do
            End do
         End do
      End do
    end subroutine Initial_Guess_VecX

  end subroutine Calc_Chi_Head

  subroutine Calc_LinearEq( ncolumn, matsize, MatA, VecX )
    integer, intent(in)                 :: ncolumn, matsize
    complex(kind=CMPLDP), intent(inout) :: MatA( matsize, matsize )
    complex(kind=CMPLDP), intent(inout) :: VecX( matsize, ncolumn )

    integer                 ::  kt_info
    integer,    allocatable :: iwork(:)
    integer :: id_sname = -1
    !
    call tstatc0_begin('Calc_Linear_Equation ', id_sname)
    !
    Allocate( iwork( matsize ) )
    Call ZGESV( matsize, ncolumn, MatA, matsize, iwork, VecX, &
         &      matsize, kt_info )
    if ( kt_info /= 0 ) write(nfout,*) 'kt_info = ', kt_info

    Deallocate( iwork )
    call tstatc0_end(id_sname)

  end subroutine Calc_LinearEq

!---------------------------------------------------------
!!
!!!          Limit band indices in consideraiton
!!
!----------------------------------------------------------
 subroutine set_nband_LR
    integer ispin, k, ik1, ib1
    Real(kind=DP) :: width, emin, emax, ebi
    Real(kind=DP) :: occ1, e_homo, e_lumo
    integer :: k_max
! ------------
    k_max = 1
    e_homo = -1.0D+3;   e_lumo= 1.0D+3
! --------------------------
    Do ispin=1, nspin
       Do k=1, k_max
          ik1 = nspin*(k-1) + ispin
          Do ib1=1, neg           ! unoccpuied state
             ebi = Eko_ek( ib1,ik1 )
             occ1 = Occup_lkt_ek( ib1, ik1 );
             if ( occ1 > Delta )  then
                e_homo = max( ebi, e_homo )
             else
                e_lumo = min( ebi, e_lumo )
             endif
          End do
       End do
    End do
! ----------------
!    Do ib1=1, neg
!       write(80+mype,*) ib1, eko_ek( ib1,1 ), occup_lkt_ek( ib1,1 )
!    End do
!    stop
    if ( mype == 0 ) then
       write(nfout,*) '! Highest  occupied level ', e_homo
       write(nfout,*) '! Lowest unoccupied level ', e_lumo
    endif
! ------------------------
    width = e_high
    emin = e_lumo - width; emax = e_homo + width
! --------------------------
    Do ispin=1, nspin
       Do k=1, k_max
          ik1 = nspin*(k-1) + ispin
          Do ib1=1, neg           ! unoccpuied state
             ebi = Eko_ek( ib1,ik1 )
             if ( ebi < emin ) then
                band_start_LR = ib1
             else if ( ebi < emax ) then
                band_end_LR = ib1
             endif
          End do
       End do
    End do
! -------------------- temp --------
!    band_start_LR = 6
!    band_end_LR = 7
! ----------------------------------
    nband_LR = band_end_LR - band_start_LR +1
    if ( mype == 0 ) then
       write(nfout,*) '!* band LR : start, end = ', band_start_LR, band_end_LR
    endif
  end subroutine set_nband_LR

!---------------------------------------------------------
!!
!!!          Solver for BS-like Equation
!!
!----------------------------------------------------------
  subroutine Calc_MatL
    Complex(kind=CMPLDP), allocatable :: MatA(:,:), work(:)
    integer, allocatable :: ipiv(:)
    integer matsize, lwork, info
!
    Complex(kind=CMPLDP) :: ztmp
    integer :: ix1, ix2, ib1, ib2
! ---
    matsize = nband_LR**2
    lwork = matsize
! --
    if ( mype /= 0 ) return

    allocate( MatA(matsize,matsize) ); MatA = 0.0D0
    allocate( ipiv(matsize) ); ipiv = 0
    allocate( work(lwork) ); work = 0.0d0
!
    Do ix1=1, matsize
       ib1 = ( ix1-1 ) /nband_LR +1
       ib2 = mod( ix1-1,nband_LR ) + 1
       ztmp = 0.0D0
       Do ix2=1, matsize
          MatA(ix1,ix2) = -MatL0( ib1,ib2,1 ) *MatXi(ix1,ix2) *dcmplx(0.0D0,1.0D0)
       End do
    End do
    Do ix1=1, matsize
       MatA(ix1,ix1) = MatA(ix1,ix1) + 1.0D0
    End do
!
    Call ZGETRF( matsize, matsize, MatA, matsize, ipiv, info )
    if ( info /= 0 ) write(*,*) 'INFO A2 = ', info
    Call ZGETRI( matsize, MatA, matsize, ipiv, work, lwork, info )
    if ( info /= 0 ) write(*,*) 'INFO B2 = ', info
!
    Do ix1=1,matsize
       Do ix2=1,matsize
          ib1 = ( ix2-1 ) /nband_LR +1
          ib2 = mod( ix2-1,nband_LR ) + 1
          MatL(ix1,ix2) = MatA(ix1,ix2) *MatL0(ib1,ib2,1)
       End do
    End do
    deallocate( ipiv, work )
    deallocate( MatA )

  end subroutine Calc_MatL

  subroutine Calc_Chi_Head_Mol( ene )
    real(kind=DP), intent(in) :: ene

    Complex(kind=CMPLDP) :: ztmp1, ztmp2, z0, z1, z2
    integer ix1, ix2, ib1, ib2, ib3, ib4
    integer ispin
! ---------------
    if ( npes > 1 ) call mpi_barrier( MPI_CommGroup,ierr )
! --
    ztmp1 = 0.0D0
    ztmp2 = 0.0D0
    z0 = 0.0D0
    Do ix1=1, nband_LR**2
       ib1 = ( ix1-1 ) /nband_LR +1
       ib2 = mod( ix1-1,nband_LR ) + 1

       if ( mype == 0 ) z0 = MatL0( ib1, ib2, 1 )

       ib1 = ib1 + band_start_LR -1;     ib2 = ib2 + band_start_LR -1

       call set_value_z1
       ztmp1 = ztmp1 -conjg(z1) *z1 *z0 *dcmplx(0.0,1.0)

       Do ix2=1, nband_LR**2
!!!!          if ( ix1 /= ix2 ) cycle

          ib3 = ( ix2-1 ) /nband_LR +1
          ib4 = mod( ix2-1,nband_LR ) + 1

          ib3 = ib3 + band_start_LR -1;     ib4 = ib4 + band_start_LR -1

          call set_value_z2
          if ( mype == 0 ) then
             ztmp2 = ztmp2 -z1 *conjg(z2) * MatL( ix1,ix2 ) *dcmplx( 0.0D0, 1.0D0 )
          endif
       End do
    End do
!    if ( mype == 0 ) then
!!       write(91,*) ene*Hartree, aimag(ztmp1) /univol, aimag(ztmp2) /univol
!    endif
!--
    if ( mype == 0 ) then
       Mat_Chi0(1,1,1) = ztmp1 / univol
       Chi_Head(1,1) = ztmp2 / univol
    endif

  contains

    subroutine set_value_z1
      integer :: jb2

      if ( map_e(ib2) == myrank_e ) then
         jb2 = map_z(ib2)
         z1 = RhoTilde( 1, ib1, jb2, 1 )
      endif
      if ( npes > 1 ) then
         call mpi_bcast( z1, 2, MPI_DOUBLE_PRECISION, &
              &          map_e(ib2), MPI_CommGroup, ierr )
      endif
    end subroutine set_value_z1

    subroutine set_value_z2
      integer :: jb4

      if ( map_e(ib4) == myrank_e ) then
         jb4 = map_z(ib4)
         z2 = RhoTilde( 1, ib3, jb4, 1 )
      endif
      if ( npes > 1 ) then
         call mpi_bcast( z2, 2, MPI_DOUBLE_PRECISION, &
                  &          map_e(ib4), MPI_CommGroup, ierr )
      end if
    end subroutine Set_value_z2

  end subroutine Calc_Chi_Head_Mol

! ------------------------------------------------------------
!!
!!!         Output for Spectrum Data
!!
! ------------------------------------------------------------
  subroutine Write_Spectrum_Header
    if ( mype /= 0 ) return

    select case(spectrum_type)
    case (OPTICS)
       write(nf_LR_spectra,'("#",16x," Optical spectrum ")')
       write(nf_LR_spectra,'("#",27x," NonInteracting ",14x, "Interacting")')
       write(nf_LR_spectra,'("#",6x,"Energy[eV]", 8x, "Real",8x,"Imaginary",9x,"Real",8x,"Imaginary")')
    case (PACS)
       write(nf_LR_spectra,'("#",16x," Photo Absorption Cross Section")')
       write(nf_LR_spectra,'("#",6x,"Energy[eV]", 2x," NonInteracting ", 2x, "Interacting")')
    case (EELS)
       write(nf_LR_spectra,'("#",16x," EELS")')
       write(nf_LR_spectra,'("#",6x,"Energy[eV]", 2x," NonInteracting ", 2x, "Interacting")')
    case (IXSS)
       write(nf_LR_spectra,'("#",16x," IXSS")')
       write(nf_LR_spectra,'("#",6x,"Energy[eV]", 2x," NonInteracting ", 2x, "Interacting")')
    end select

  end subroutine Write_Spectrum_Header

  subroutine Write_Spectrum_Data( ene )
    real(kind=DP), intent(in) :: ene

    integer :: ni, nj, i
    real(kind=DP) :: c1
    real(kind=DP) :: au_of_speed_of_light
    complex(kind=CMPLDP)      :: z0, z1, eps0, eps1
!
    real(kind=DP) :: pot_v0, ttr(6), ga, gb, gc, g2
! --------------------------------- start ------------
    if ( mype /= 0 ) return
    Call getttr (rltv,ttr)
! --------------------------------- main -------------
    pot_v0 = Kernel_Coulomb(1)
! ---------------
    Do i=1, 1
       ga = ngabc(i,1) + vqxyz( 1,1,BUCS )
       gb = ngabc(i,2) + vqxyz( 1,2,BUCS )
       gc = ngabc(i,3) + vqxyz( 1,3,BUCS )
       g2 = ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
              &   +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga
       c1 = 1.0d0
       pot_v0 = PAI4 / g2 *c1
    End Do
    !
    z0 = 0.0d0;    z1 = 0.0d0
    Do ni=1, nspin
       z0 = z0 + Mat_Chi0(1,1,ni)
    End do
    Do ni=1, nspin
       Do nj=1, nspin
          z1 = z1 + Chi_Head(ni,nj)
       End do
    End do
    eps0 = 1.0d0 - pot_v0 *z0;      eps1 = 1.0d0 - pot_v0 *z1

    select case (spectrum_type)
    case (OPTICS)
       write(nf_LR_spectra,'(5F15.6)') ene *Hartree, &
            &                          real(eps0), aimag(eps0), &
            &                          real(eps1), aimag(eps1)
    case (PACS)
       au_of_speed_of_light = speed_of_light/ AU_VELOCITY
       c1 = univol / au_of_speed_of_light
       write(nf_LR_spectra,'(3F15.6)') ene *Hartree, &
            &                          c1 *ene *aimag(eps0), &
            &                          c1 *ene *aimag(eps1)
    case (EELS)
       write(nf_LR_spectra,'(3F15.6)') ene *Hartree, &
            &                          aimag(eps0), &
            &                          aimag(eps1)
    case (IXSS)
       c1 = -2.0d0 *univol
       write(nf_LR_spectra,'(3F15.6)') ene *Hartree, &
            &                          c1 *aimag(z0), &
            &                          c1 *aimag(z1)
    end select
  end subroutine Write_Spectrum_Data

end module m_LinearResponse_Spectrum
