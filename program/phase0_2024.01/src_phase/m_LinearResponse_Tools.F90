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
!  This is a module for tools used in the LR-TDDFT.
! ======================================================================

! ======================= history ======================================
!  ver 1.0 :  2010/3/31
!               applicable to the crystal such as Si.
!  ver 2.0 :  2011/3/31
!               applicable to the isolated molecule such as C6H6.
!
! ======================================================================

module m_LinearResponse_tools

  use m_Const_Parameters,               only : DP, CARTS, BUCS, PAI, PAI2, PAI4,  &
       &                                       ELECTRON, INVERSE, ON, OFF, CMPLDP, &
       &                                       GAMMA

  use m_Files,                          only : nfout
  use m_FFT,                            only : nfft, m_FFT_WF

  use m_Kpoints,                        only : kv3_ek, kv3, vkxyz_ek, qwgt_ek, &
       &                                       vkxyz, qwgt, &
       &                                       k_symmetry
  use m_Control_Parameters,             only : ipri, nspin, kimg, af, neg, printable

  use m_Parallelization,                only : ista_e,iend_e,istep_e,map_z,np_e, &
    &                                          ista_k, iend_k, map_k, map_ek,  &
    &                                          myrank_k, myrank_e, map_e, npes, &
    &                                          MPI_CommGroup, ierr, mype, istatus

  use m_Electronic_Structure, only : zaj_l, occup_l, eko_l, eko_ek, neordr, efermi
  use m_PlaneWaveBasisSet,    only : ngabc,igf,kg1,kg, iba, nbase, nbase_gamma
  use m_IterationNumbers,     only : nk_in_the_process, nk_converged

  use m_Timing,                     only : tstatc0_begin, tstatc0_end

!  use m_LinearResponse_Control,     only : Num_Q_Points
  use m_Control_Parameters,     only : Num_Q_Points, sw_LinearResponse

! ----------------------------
  use m_LinearResponse_Control,         only : vqxyz, vec_q, nrd_efermi
  use m_LinearResponse_Control,         only : sw_tddft,  sw_LongWaveLimit
! ------------------- USPP ---
  use m_PseudoPotential,        only : modnrm, nlmta
  use m_Const_Parameters,       only : EXECUT
  use m_Electronic_Structure,    only : fsr_l, fsi_l

  use m_FFT,                    only : nfftp, m_FFT_CD_INVERSE_c
  use m_Parallelization,         only : ista_kngp, iend_kngp, &
       &                                ista_fftp, iend_fftp
  use m_PlaneWaveBasisSet,    only : kgp
  use mpi
! -----------------------------
  Implicit None
!  include 'mpif.h'

  ! --------------------- temporary work arrays ------------------------
  real(kind=DP), allocatable :: wfn_k(:,:,:,:), wfn_kmq(:,:,:,:)
  integer,  allocatable      :: map_z_k(:), map_z_kmq(:)
  integer,  allocatable      :: nbase_k( :,: ), nbase_kmq(:,:)
  integer,  allocatable      :: igf_k(:), igf_kmq(:), iba_k(:), iba_kmq(:)
  !
  real(kind=DP), allocatable :: occup_k(:,:), occup_kmq(:,:)
  real(kind=DP), allocatable :: eko_k(:,:), eko_kmq(:,:)
! -------------------- For USPP -------------
  real(kind=DP), allocatable :: fsval_k(:,:,:,:), fsval_kmq(:,:,:,:)
 ! ---------------------------------------------------------

contains

!----------------------------------------------------
!!
!!!         Alloc and Dealloc temporary arrays
!!
!----------------------------------------------------
  subroutine m_LR_alloc_Array_WFns
    if ( sw_LongWaveLimit == OFF ) then
       Allocate( wfn_k  ( kg1, np_e, kv3, kimg) ); wfn_k = 0.0d0
       Allocate( wfn_kmq( kg1, np_e, kv3, kimg) ); wfn_kmq = 0.0d0
    else
       Allocate( wfn_k  ( kg1, np_e, kv3, kimg) ); wfn_k = 0.0d0
    endif
  end subroutine m_LR_alloc_Array_WFns

  subroutine m_LR_dealloc_Array_WFns
    if ( sw_LongWaveLimit == OFF ) then
       Deallocate( wfn_k, wfn_kmq )
    else
       Deallocate( wfn_k )
    endif
  end subroutine m_LR_dealloc_Array_WFns

  subroutine m_LR_alloc_Array_MapFns
    if ( sw_LongWaveLimit == OFF ) then
       Allocate( map_z_k  ( neg ) );      map_z_k = 0
       Allocate( map_z_kmq( neg ) );      map_z_kmq = 0
       Allocate( nbase_k  ( kg1, kv3 ) ); nbase_k = 0
       Allocate( nbase_kmq( kg1, kv3 ) ); nbase_kmq = 0
       Allocate( igf_k  ( kg ) );         igf_k = 0
       Allocate( igf_kmq( kg ) );         igf_kmq = 0
       Allocate( iba_k  ( kv3 ) );        iba_k = 0
       Allocate( iba_kmq( kv3 ) );        iba_kmq = 0
    else
       Allocate( map_z_k  ( neg ) );      map_z_k = 0
       Allocate( nbase_k  ( kg1, kv3 ) ); nbase_k = 0
       Allocate( igf_k  ( kg ) );         igf_k = 0
       Allocate( iba_k  ( kv3 ) );        iba_k = 0
    endif
  end subroutine m_LR_alloc_Array_MapFns

  subroutine m_LR_dealloc_Array_MapFns
    if ( sw_LongWaveLimit == OFF ) then
       Deallocate( map_z_k, map_z_kmq ); Deallocate( nbase_k, nbase_kmq )
       Deallocate( igf_k, igf_kmq );     Deallocate( iba_k, iba_kmq )
    else
       Deallocate( map_z_k );        Deallocate( nbase_k )
       Deallocate( igf_k );          Deallocate( iba_k )
    endif
  end subroutine m_LR_dealloc_Array_MapFns

  subroutine m_LR_alloc_Array_Occups
    if ( sw_LongWaveLimit == OFF ) then
       Allocate( occup_k( neg, ista_k:iend_k) );   occup_k = 0.0d0
       Allocate( occup_kmq( neg, ista_k:iend_k) ); occup_kmq = 0.0d0
    else
       Allocate( occup_k( neg, ista_k:iend_k) ); occup_k = 0.0d0
    endif
  end subroutine m_LR_alloc_Array_Occups

  subroutine m_LR_dealloc_Array_Occups
    if ( sw_LongWaveLimit == OFF ) then
       deallocate( occup_k, occup_kmq )
    else
       deallocate( occup_k )
    endif
  end subroutine m_LR_dealloc_Array_Occups

  subroutine m_LR_alloc_Array_EigenVals
    if ( sw_LongWaveLimit == OFF ) then
    else
       Allocate( eko_k( neg, ista_k:iend_k ) );     eko_k = 0.0d0
    endif
  end subroutine m_LR_alloc_Array_EigenVals

  subroutine m_LR_dealloc_Array_EigenVals
    if ( sw_LongWaveLimit == OFF ) then
    else
       deallocate( eko_k )
    endif
  end subroutine m_LR_dealloc_Array_EigenVals

  subroutine m_LR_alloc_Array_CorePart
    if ( sw_LongWaveLimit == OFF ) then
       Allocate( fsval_k  ( np_e, nlmta, ista_k:iend_k, kimg ) ); fsval_k   = 0.0d0
       Allocate( fsval_kmq( np_e, nlmta, ista_k:iend_k, kimg ) ); fsval_kmq = 0.0d0
    else
       Allocate( fsval_k  ( np_e, nlmta, ista_k:iend_k, kimg ) ); fsval_k   = 0.0d0
    endif
  end subroutine m_LR_alloc_Array_CorePart

  subroutine m_LR_dealloc_Array_CorePart
    if ( sw_LongWaveLimit == OFF ) then
       Deallocate( fsval_k, fsval_kmq )
    else
       Deallocate( fsval_k )
    endif
  end subroutine m_LR_dealloc_Array_CorePart
!-------------------------------------------------------------
!!
!!!      obtain a wavefunction in the real space
!!
!--------------------------------------------------------------
  subroutine Get_WF_in_Rspace( ik, ib, wfn_in, bfft, &
       &                       nbase_in, iba_in, igf_in, map_z_in )
    real(kind=DP), intent(out) :: bfft(nfft)
    real(kind=DP), intent(in) :: wfn_in( kg1, np_e, kv3, kimg )

    integer, intent(in) :: ik, ib
    integer, intent(in) :: nbase_in( kg1,kv3 )
    integer, intent(in) :: igf_in( kg )
    integer, intent(in) :: iba_in( kv3 )
    integer, intent(in) :: map_z_in( neg )

    integer :: i,i1, ri, j, i2, ii
! --------------------------------------- start --------------
    bfft = 0.d0
! --------------------------------------- main ---------------
    if ( k_symmetry(ik) == GAMMA ) then
       if ( kimg == 1 ) then
!         write(*,*) 'YYYYYYYYYYYYYYYYYYYYYYYYYYY 7'
          i1 = igf_in(1)
          bfft(i1) = wfn_in(1, ib, ik, 1)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba_in(ik)
             i = nbase_in(ii,ik)
             i1 = igf_in(i)
             bfft(i1) = wfn_in(ii, ib, ik, 1)
             j = nbase_gamma(ii,2)
             i2 = igf_in(j)
             bfft(i2) = wfn_in(ii, ib, ik, 2)
          end do
       else ! kimg == 2
!         write(*,*) 'YYYYYYYYYYYYYYYYYYYYYYYYYYY 8'
          i1 = 2*igf_in(1) - 1
          bfft(i1)   = wfn_in(1, ib, ik, 1)
          bfft(i1+1) = wfn_in(1, ib, ik, 2)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba_in(ik)
             i = nbase_in( ii,ik )
             i1 = 2*igf_in(i)-1
             bfft(i1  ) = wfn_in(ii, ib, ik, 1)
             bfft(i1+1) = wfn_in(ii, ib, ik, 2)
             j = nbase_gamma(ii,2)
             i2 = 2*igf_in(j)-1
             bfft(i2  ) =  wfn_in(ii, ib, ik, 1)
             bfft(i2+1) = -wfn_in(ii, ib, ik, 2)
          end do
       end if
    else
      if(kimg == 1) then
!         write(*,*) 'YYYYYYYYYYYYYYYYYYYYYYYYYYY 9'
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
        do i = 1, iba_in(ik)
           i1 = igf_in( nbase_in(i,ik) )
           bfft(i1) = wfn_in( i, ib, ik, 1 )
        end do
      else ! kimg == 2
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
!         write(*,*) 'YYYYYYYYYYYYYYYYYYYYYYYYYYY 10'
        do i = 1, iba_in(ik)
           i1 = 2*igf_in( nbase_in( i,ik ) ) - 1
           bfft(i1) = wfn_in( i, ib, ik, 1 )
           i2 = 2*igf_in( nbase_in(i,ik) )
           bfft(i2) = wfn_in( i, ib, ik, 2 )
        end do
      end if
    end if
! ----------------------
    call m_FFT_WF( ELECTRON,nfout,bfft,INVERSE,OFF )
!p    call m_FFT_WF( ELECTRON,nfout,bfft,INVERSE,ON )
  end subroutine Get_WF_in_Rspace

  subroutine Get_WF_in_Rspace_CD( ik, ib, wfn_in, bfft, &
       &                          nbase_in, iba_in, igfp_LR_in, map_z_in )
! --
    real(kind=DP), intent(in) :: wfn_in( kg1, np_e, ista_k:iend_k, kimg )
    real(kind=DP), intent(inout) :: bfft(ista_fftp:iend_fftp)

    integer, intent(in) :: ik, ib
    integer, intent(in) :: nbase_in( kg1,kv3 )
    integer, intent(in) :: iba_in( kv3 )
    integer, intent(in) :: map_z_in( neg )
    integer, intent(in) :: igfp_LR_in( kgp )

    integer :: i,i1, ri, j, i2, ii
    integer :: jb

    real(kind=DP), allocatable :: bfft_mpi(:)
! -------
    bfft = 0.0d0
    if ( npes > 1 ) then
       allocate( bfft_mpi( nfftp )); bfft_mpi = 0.0d0
       if ( map_e( ib ) == myrank_e ) then
          jb = map_z_in( ib )
          do ri = 1, kimg
             do i = 1, iba_in(ik)
                i1 = kimg *igfp_LR_in( nbase_in(i,ik) ) + (ri - kimg)
                bfft_mpi(i1) = wfn_in( i, jb, ik, ri )   ! MPI
             end do
          end do
       endif
       call mpi_bcast( bfft_mpi, nfftp, MPI_DOUBLE_PRECISION, &
            &          map_e(ib), MPI_CommGroup, ierr )
       if ( ierr /= 0 ) write(*,*) 'MPI error ', mype, ierr
       bfft(ista_fftp:iend_fftp) = bfft_mpi(ista_fftp:iend_fftp)

       deallocate( bfft_mpi )
    else
       do ri = 1, kimg
          do i = 1, iba_in(ik)
             i1 = kimg *igfp_LR_in( nbase_in(i,ik) ) + (ri - kimg)
             bfft(i1) = wfn_in( i, ib, ik, ri )   ! MPI
          end do
       end do
    endif
! --
    call m_FFT_CD_inverse_c( nfout, bfft )   ! From G- to R-space

  end subroutine Get_WF_in_Rspace_CD
!--------------------------------------------------------
!!
!!!            Copying tools
!!
!--------------------------------------------------------

  subroutine Copy_Zaj_To( Wfn_target )
    real(kind=DP), intent(inout) ::  Wfn_target ( kg1, np_e, ista_k:iend_k, kimg )

! ------------------------------ start -----------
    Wfn_target( 1:kg1, 1:np_e, ista_k:iend_k, 1:kimg ) &
         & = zaj_l( 1:kg1, 1:np_e, ista_k:iend_k, 1:kimg )
! ------------------------------ end -----------
  end subroutine Copy_Zaj_To

  subroutine Copy_HardPart_To( fs_target )
    real(kind=DP), intent(inout) :: fs_target( np_e, nlmta, ista_k:iend_k, kimg )
! ---------------------------------- main ---------
    fs_target( 1:np_e, 1:nlmta, ista_k:iend_k, 1 ) &
         & = fsr_l( 1:np_e, 1:nlmta, ista_k:iend_k )
    if ( kimg ==2 ) then
       fs_target( 1:np_e, 1:nlmta, ista_k:iend_k, 2 ) &
            & = fsi_l( 1:np_e, 1:nlmta, ista_k:iend_k )
    endif
  end subroutine Copy_HardPart_To

  subroutine Copy_EigenVals_To( Eko_target, occup_tmp )
    real(kind=DP), intent(out)           :: Eko_target( neg, kv3 )
    real(kind=DP), intent(out), optional :: occup_tmp( neg, kv3 )
!
    real(kind=DP), allocatable :: e_mpi(:,:), e2_mpi(:,:)
    integer :: ik, ie
! -------------------------- start -------
    allocate( e_mpi(neg,kv3));   e_mpi = 0.0d0
    allocate(e2_mpi(neg,kv3));  e2_mpi = 0.0d0
! --------------------------- main -------
    Do ik = 1, kv3
       if ( map_k(ik) /= myrank_k ) cycle
       Do ie = 1, neg
          if ( map_e(ie) /= myrank_e ) cycle
          e_mpi( ie,ik ) = eko_l( map_z(ie),ik )
       End do
    End do
    if( npes >= 2 ) then
       call mpi_allreduce(e_mpi,e2_mpi,neg*kv3,mpi_double_precision &
            &               ,mpi_sum,MPI_CommGroup,ierr)
    else
       e2_mpi = e_mpi
    end if
!
    Eko_target = e2_mpi
    if ( nrd_efermi==ON ) then
       Call Estimate_Ocuup_Simply( e2_mpi, occup_tmp )
    endif
! ------------------- end ------------
    deallocate( e_mpi, e2_mpi )
  end subroutine Copy_EigenVals_To

  subroutine Copy_MapZ_To( MapFunction )
    integer, intent(inout) :: MapFunction( neg )
    MapFunction  = map_z
  end subroutine Copy_MapZ_To

  subroutine Copy_MapE_To( MapFunction )
    integer, intent(inout) :: MapFunction( neg )
    MapFunction  = map_e
  end subroutine Copy_MapE_To

  subroutine Copy_MapEK_To( MapFunction )
    integer, intent(inout) :: MapFunction( neg, kv3 )
    MapFunction  = map_ek
  end subroutine Copy_MapEK_To

  subroutine Copy_NBase_To( MapFunction )
    integer, intent(inout) :: MapFunction( kg1, kv3 )
    MapFunction  = nbase
  end subroutine Copy_NBase_To

  subroutine Copy_Igf_To( MapFunction )
    integer, intent(inout) :: MapFunction( kg  )
    MapFunction  = igf
  end Subroutine Copy_Igf_To

  subroutine Copy_Iba_To( MapFunction )
    integer, intent(inout) :: MapFunction( kv3 )
    MapFunction  = iba
  end subroutine Copy_Iba_To

  subroutine Copy_Eko_Or_Occ_To( MatA, MatB )
    real(kind=DP), intent(in)    ::  MatA( neg, kv3 )
    real(kind=DP), intent(inout) ::  MatB( neg, kv3_ek )
    integer i, ik

    Do ik=1, kv3
       Do i=1, neg
          MatB( i, ik + nk_in_the_process -1 ) = MatA( i,ik )
       End do
    End do
  end subroutine Copy_Eko_Or_Occ_To

  subroutine Copy_Eko_To( MatA, MatB )
    real(kind=DP), intent(in)    ::  MatA( neg, kv3 )
    real(kind=DP), intent(inout) ::  MatB( neg, kv3_ek )
    integer i, ik

    Do ik=1, kv3
       Do i=1, neg
          MatB( i, ik + nk_in_the_process -1 ) = MatA( i,ik )
       End do
    End do
  end subroutine Copy_Eko_To

  subroutine Copy_EigenVals_To2( Eko_target, Occup_target, occup_tmp )
    real(kind=DP), intent(out)           :: Eko_target( neg, kv3_ek )
    real(kind=DP), intent(out), optional :: Occup_target( neg, kv3_ek )
    real(kind=DP), intent(out), optional :: occup_tmp( neg, kv3 )
!
    real(kind=DP), allocatable :: e_mpi(:,:), e2_mpi(:,:)
    integer ik, ie, i
! -------------------------- start -----------
    allocate( e_mpi(neg,kv3));   e_mpi = 0.0d0
    allocate(e2_mpi(neg,kv3));  e2_mpi = 0.0d0
! -------------------------- main -----------
    Do ik = 1, kv3
       if ( map_k(ik) /= myrank_k ) cycle
       Do ie = 1, neg
          if ( map_e(ie) /= myrank_e ) cycle
          e_mpi( ie,ik ) = eko_l( map_z(ie),ik )
       End do
    End do
    if( npes >= 2 ) then
       call mpi_allreduce(e_mpi,e2_mpi,neg*kv3,mpi_double_precision &
            &               ,mpi_sum,MPI_CommGroup,ierr)
    else
       e2_mpi = e_mpi
    end if
!
    Do ik=1, kv3
       Do i=1, neg
          Eko_target( i, ik + nk_in_the_process -1 ) = e2_mpi( i,ik )
       End do
    End do
! -
    if ( nrd_efermi == ON ) then
       Call Estimate_Ocuup_Simply( e2_mpi, occup_tmp )
       Do ik=1, kv3
          Do i=1, neg
             Occup_target( i, ik + nk_in_the_process -1 ) = occup_tmp( i,ik )
          End do
       End do
    endif
! ------------------------- end ---------------------
    deallocate( e_mpi, e2_mpi )
  end subroutine Copy_EigenVals_To2

  subroutine Copy_Occups_To( occup_tmp )
    real(kind=DP), intent(out), optional :: occup_tmp( neg, kv3 )
!
    real(kind=DP), allocatable :: e_mpi(:,:), e2_mpi(:,:)
    integer ik, ie, i
! -------------------------- start ----------
    allocate( e_mpi(neg,kv3));   e_mpi = 0.0d0
    allocate(e2_mpi(neg,kv3));  e2_mpi = 0.0d0
! --------------------------- main ----------
    Do ik = 1, kv3
       if ( map_k(ik) /= myrank_k ) cycle
       Do ie = 1, neg
          if ( map_e(ie) /= myrank_e ) cycle
          e_mpi( ie,ik ) = eko_l( map_z(ie),ik )
       End do
    End do
    if( npes >= 2 ) then
       call mpi_allreduce(e_mpi,e2_mpi,neg*kv3,mpi_double_precision &
            &               ,mpi_sum,MPI_CommGroup,ierr)
    else
       e2_mpi = e_mpi
    end if
!
    if ( nrd_efermi == ON ) then
       Call Estimate_Ocuup_Simply( e2_mpi, occup_tmp )
    endif
! ------------------------- end --------------
    deallocate( e_mpi, e2_mpi )
  end subroutine Copy_Occups_To

  subroutine Estimate_Ocuup_Simply( eko_in, occup_out )
    real(kind=DP), intent(in)  :: eko_in( neg, kv3 )
    real(kind=DP), intent(out) :: occup_out( neg, kv3 )

    integer ik, j

    Do ik=1, kv3
       Do j=1, neg
          if ( eko_in( j,ik) < efermi ) then
             occup_out( j,ik ) = 1.0d0
          else
            occup_out( j,ik ) = 0.0d0
          endif
       End do
    End do
  end subroutine Estimate_Ocuup_Simply

end module m_LinearResponse_tools

