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
!  This is a module for calculating < n'k | exp(i (q+G) r ) | n k-q >
!  for the LR-TDDFT .
! ======================================================================

! ======================= history ======================================
!  ver 1.0 :  2010/3/31
!               applicable to the crystal such as Si.
!  ver 2.0 :  2011/3/31
!               applicable to the isolated molecule such as C6H6.
!
! ======================================================================

module m_LinearResponse_Density

  use m_Files,                      only : nfout, nfpot, &
       &                                   m_Files_open_ps_files, &
       &                                   m_Files_close_ps_files, &
       &                                   m_Files_open_ps_file, &
       &                                   m_Files_close_ps_file

  use m_FFT,                        only : nfft, fft_box_size_WF, m_FFT_wf, &
       &                                   m_fft_alloc_wf_work, &
       &                                   m_fft_dealloc_wf_work
  use m_Const_Parameters,           only : DP, GAMMA, ELECTRON, DIRECT, &
       &                                   INVERSE, ON, OFF, CMPLDP, DELTA, &
       &                                   CARTS, BUCS, Valence_Plus_PC_Charge, &
       &                                   Hartree, NO, YES, &
       &                                   TETRAHEDRON, PARABOLIC, &
       &                                   PAI, PAI2, CMPLDP, EXECUT, SKIP, &
       &                                   DELTA_FermiSearchRange

  use m_Kpoints,                    only : kv3_ek, kv3, vkxyz_ek, qwgt_ek, &
       &                                   k_symmetry, np0, np2
  use m_Control_Parameters,         only : ipri, nspin, kimg, af, neg, printable, &
       &                                   way_of_smearing,width_tetra

  use m_LinearResponse_Control,     only  : nrd_efermi, nstep, e, nmax_G_LR, &
       &                                    vqxyz,  &
       &                                    scissor, eta

  use m_LinearResponse_tools,  only : wfn_k, wfn_kmq, map_z_k, map_z_kmq,  &
       &                                   nbase_k, nbase_kmq, &
       &                                   igf_k, igf_kmq, iba_k, iba_kmq

  use m_LinearResponse_tools,  only  : occup_k, occup_kmq, &
       &                                   Get_WF_in_Rspace, &
       &                                    eko_k, &
       &                                    fsval_k, fsval_kmq

  use m_PlaneWaveBasisSet,    only : ngabc,igf, kg0, kg1,kg, kgp, iba, &
       &                             nbase, nbase_gamma,   nbmx

  use m_IterationNumbers,     only : nk_in_the_process, nk_converged

  use m_Kpoints,                        only : kv3_ek, kv3, vkxyz_ek, qwgt_ek, &
       &                                       vkxyz, qwgt, nxyz_tetra, ip20, iwt, &
       &                                       ip2cub

  use m_Crystal_Structure,              only : altv,rltv, univol
  use m_Parallelization,                only : ista_e,iend_e,istep_e,map_z,np_e, &
    &                                          ista_k, iend_k, map_k, map_ek,  &
    &                                          myrank_k, myrank_e, map_e, npes, &
    &                                          MPI_CommGroup, ierr, mype, istatus, &
    &                                          ista_snl, iend_snl

  use m_Electronic_Structure,       only : fsr_l,fsi_l,fsr_add_l,fsi_add_l,vlhxcQ, &
       &                                    zaj_l
  use m_ES_nonlocal,                only : m_ES_add_betar_dot_WFs
  use m_Ionic_System,               only : pos, cps, ntyp,ityp,iwei,natm,natm2 ,ivan,&
       &                                   iatomn
  use m_PseudoPotential,            only : ival,ilmt,nlmt,nlmtt,nlmta,lmta,lmtt, &
       &                                   ltp,mtp, q, dion, modnrm,   &
       &                                   nac,fqwei,ilmt_add,nlmta_add, &
       &                                   ltp_add,mtp_add,lmta_add,lpsmax, &
       &                                   m_PP_include_vanderbilt_pot, &
       &                                   m_PP_tell_lmtt_l_m_tau, &
       &                                   m_PP_make_index_lmtt_add, &
       &                                   m_PP_tell_lmtt_l_m_tau_add

  use m_Timing,                     only : tstatc0_begin, tstatc0_end

  use m_Control_Parameters,         only : Num_q_Points, m_CtrlP_way_of_smearing, &
       &                                   iprioccup, width, sw_use_add_proj, &
       &                                   paramset, icond
  use m_PseudoPotential,            only : vec_q_plus_G_LR
  use m_Parallelization,            only : ista_kngp, iend_kngp

  use m_Electronic_Structure,       only : efermi, eko_ek, totch, metalic_system, vbm
  use m_LinearResponse_Control,     only : sw_LongWaveLimit
!
  use m_LinearResponse_Qpt,         only :ftqval

!! UVSOR
  use m_ES_occup,             only : check_totch, check_if_metalic
!! UVSOR

  implicit none
  include 'mpif.h'

! ------------------------------
  complex(kind=CMPLDP), allocatable :: RhoTilde( :,:,:,: )
! -------------------------------------------
  real(kind=DP),        allocatable :: occup_lkt_ek( :,: )
  real(kind=DP),        allocatable :: qwgt_All(:)

! ----------------------------------------------
  real(kind=DP), allocatable :: dipole_dxyz_us(:,:,:)
  real(kind=DP), allocatable :: ptrans(:,:,:,:,:)

  integer,  allocatable :: nppc_data(:)
  integer,  allocatable :: phase_ylm1(:,:),phase_ylm2(:,:)
  integer,  allocatable :: dipole_tau1(:,:),dipole_tau2(:,:)

  integer,  allocatable :: ilocal_l(:)
  integer,  allocatable :: PP_norm_type(:), PP_local_type(:)

  integer               :: nppcorr, mnppc
!
  integer  :: nonlocal=0, n_check_ts = 1
  integer  :: NC_PP = 1, US_PP = 2
  integer  :: BHS_POLY=1, ORBITAL=2

contains

!----------------------------------------------------
!!
!!!         Alloc and Dealloc RhoTilde
!!
!----------------------------------------------------
  subroutine m_LR_alloc_Array_Densities
    if ( sw_LongWaveLimit == OFF ) then
       Allocate( RhoTilde( kg1, neg, np_e, kv3_ek / (Num_q_Points+1) ) )
    else
       Allocate( RhoTilde( kg1, neg, np_e, kv3_ek ) )
    endif
    RhoTilde = 0.0d0
  end subroutine m_LR_alloc_Array_Densities

  subroutine m_LR_dealloc_Array_Densities
    Deallocate( RhoTilde )
  end subroutine m_LR_dealloc_Array_Densities


  subroutine m_LR_alloc_Array_Occups_ek
    Allocate( Occup_lkt_ek( neg,kv3_ek) ) ; occup_lkt_ek = 0.0d0
  end subroutine m_LR_alloc_Array_Occups_ek

  subroutine m_LR_dealloc_Array_Occups_ek
    Deallocate( Occup_lkt_ek )
  end subroutine m_LR_dealloc_Array_Occups_ek

  subroutine m_LR_alloc_Array_Eko_ek
    Allocate( eko_ek( neg,kv3_ek) ) ; eko_ek = 0.0d0
  end subroutine m_LR_alloc_Array_Eko_ek

  subroutine m_LR_dealloc_Array_Eko_ek
    Deallocate( eko_ek )
  end subroutine m_LR_dealloc_Array_Eko_ek

!------------------------------------------------------------------
!!
!!!            Calc RhoTilde in the case of Q > 0   ( soft part )
!!
!-------------------------------------------------------------------
  subroutine Add_SoftPart_General
    Real(kind=DP), allocatable :: wfn_k_Rspace(:,:), Wfn_kmq_Rspace(:,:)
    Real(kind=DP), allocatable :: chgr(:), chgi(:)
    Real(kind=DP), allocatable :: wf1_tmp(:), wf2_tmp(:)

    Real(kind=DP) occ1, occ2, c1
    integer ngrid, nffth
    integer ispin, ik
    integer ib1, ib2, jb2
    integer :: id_sname = -1
! -------------------------------- start ------------
    call tstatc0_begin( 'Add_SoftPart_General', id_sname )
!
    ngrid = product(fft_box_size_WF(1:3,1));   nffth = nfft/2
    call m_FFT_alloc_WF_work()
    call set_work_arrays
! ---------------------------------- main --------------------
! ------------------ < Phi_ib1 | exp(igr) | Phi_ib2 > --------
! ------------------------------------------------------------
    Do ik=1, kv3
!
       Do ib1=1, neg                 ! unoccupied Band

          if ( nrd_efermi ==1 ) then
             occ1 = occup_k( ib1,ik )
             if ( occ1 > DELTA ) cycle
          endif
          call set_wf1_tmp

          Do ib2=1, neg                ! occupied Band

             if ( nrd_efermi==1 ) then
                occ2 = occup_kmq( ib2,ik )
                if ( occ2 < DELTA ) cycle
             end if

             if ( map_e(ib2) /= myrank_e ) cycle
             jb2 = map_z_kmq(ib2)

             call set_wf2_tmp
             call set_RhoTilde_in_Gspace
          End do
       End do
    End do
! --------------------------- end ------------
    call unset_work_arrays
    Call m_FFT_dealloc_WF_work()
    call tstatc0_end(id_sname)

  contains
! ---------
    subroutine set_work_arrays
      allocate( wfn_k_Rspace(nfft,np_e) );   wfn_k_Rspace = 0.0d0
      allocate( wfn_kmq_Rspace(nfft,np_e) ); wfn_kmq_Rspace = 0.0d0
      allocate( chgr(nffth)); allocate(chgi(nffth)); chgr = 0.0d0; chgi = 0.0d0
      allocate( wf1_tmp(nfft) ); wf1_tmp = 0.0d0
      allocate( wf2_tmp(nfft) ); wf2_tmp = 0.0d0

      Do ik=1, kv3
         Do ib1=1, np_e
            Call  Get_WF_in_Rspace( ik, ib1, wfn_k, wf1_tmp, &
                 &                  nbase_k, iba_k, igf_k, map_z_k )
            wfn_k_Rspace(1:nfft,ib1) = wf1_tmp(1:nfft)
         End do
      End do
      Do ik=1, kv3
         Do ib1=1, np_e
            Call  Get_WF_in_Rspace( ik, ib1, wfn_kmq, wf1_tmp, &
                 &                  nbase_kmq, iba_kmq, igf_kmq, map_z_kmq )
            wfn_kmq_Rspace(1:nfft,ib1) = wf1_tmp(1:nfft)
         End do
      End do
!
      wfn_k_Rspace   = wfn_k_Rspace / sqrt(univol)
      wfn_kmq_Rspace = wfn_kmq_Rspace / sqrt(univol)
    end subroutine set_work_arrays

    subroutine unset_work_arrays
      deallocate( chgr )
      if ( allocated( chgi) ) deallocate( chgi )
      deallocate( wfn_k_Rspace, wfn_kmq_Rspace )
      deallocate( wf1_tmp, wf2_tmp )
    end subroutine unset_work_arrays

    subroutine set_wf1_tmp
      integer :: jb1

      wf1_tmp = 0.0d0
      if ( map_e(ib1) == myrank_e ) then
         jb1 = map_z_k(ib1)
         wf1_tmp(1:nfft) = wfn_k_Rspace(1:nfft,jb1)
      endif
      if ( npes > 1 ) then
         call mpi_bcast( wf1_tmp, nfft, MPI_DOUBLE_PRECISION, &
              &         map_e(ib1), MPI_CommGroup, ierr )
         if ( ierr /= 0 ) write(*,*) 'MPI error ', mype, ierr
      endif
    end subroutine set_wf1_tmp

    subroutine set_wf2_tmp
      wf2_tmp(1:nfft) = wfn_kmq_Rspace(1:nfft,jb2)
    end subroutine set_wf2_tmp

    subroutine set_RhoTilde_in_Gspace
      integer i, ip
      integer ic1, i1

      if ( k_symmetry(ik) == GAMMA ) then
         do i=1,nffth
            ip = 2*i-1
            chgr(i) = wf1_tmp(ip)  *wf2_tmp(ip)
            chgi(i) =-wf1_tmp(ip+1)*wf2_tmp(ip)
         end do
      else
         do i=1,nffth
            ip = 2*i-1
            chgr(i) = wf1_tmp(ip)*wf2_tmp(ip)  +wf1_tmp(ip+1)*wf2_tmp(ip+1)
            chgi(i) = wf1_tmp(ip)*wf2_tmp(ip+1)-wf1_tmp(ip+1)*wf2_tmp(ip)
         end do
      end if
! ---
      wf2_tmp = 0.0d0
      do i=1,nffth
         ip = 2*i-1
         wf2_tmp(ip)   =  chgr(i);  wf2_tmp(ip+1) =  chgi(i)
      end do
! ----------------------------- FFT : Rspace -> G-space ---------
      call m_FFT_WF( ELECTRON,nfout, wf2_tmp, DIRECT,OFF )
      wf2_tmp = wf2_tmp / dble(ngrid) *univol
! ----------------------------------------------------------------
      ic1 = ( nk_in_the_process -1 )/ ( nspin*(Num_q_Points+1) )
      Do i=1, kg1
         i1 = igf_k( nbase_k(i,ik) )
         ip = 2*i1 - 1
         RhoTilde( i, ib1, jb2, ik + ic1*nspin ) &
              & = dcmplx( wf2_tmp(ip), wf2_tmp(ip+1) )
      End do
    end subroutine set_RhoTilde_in_Gspace

  end subroutine Add_SoftPart_General

!------------------------------------------------------------------
!!
!!!            Calc RhoTilde in the case of Q > 0   ( hard part )
!!
!-------------------------------------------------------------------
  Subroutine Add_HardPart_General
    Real(kind=DP), allocatable :: fs1_tmp( :,: ), fs2_tmp(:,:)
    complex(kind=CMPLDP),   allocatable  :: zsum(:,:)

    real(kind=DP) :: occ1, occ2
!
    integer :: ik, ib1, ib2, jb2
    integer :: id_sname = -1
! -------------------------------- start ------------
    call tstatc0_begin( 'Add_HardPart_General', id_sname )
    call set_work_arrays
!------------------------------------------------------
!------------------------------- main -----------------
!------------------------------------------------------
    Do ik=1, kv3

       Do ib1= 1, neg

          if ( nrd_efermi==1 ) then
             occ1 = occup_k( ib1,ik )
             if ( occ1 > DELTA ) cycle
          endif
          call set_fs1_tmp

          Do ib2=1, neg                ! occupied Band

             if ( nrd_efermi==1 ) then
                occ2 = occup_kmq( ib2,ik )
                if ( occ2 < DELTA ) cycle
             end if

             if ( map_e(ib2) /= myrank_e ) cycle
             jb2 = map_z_kmq(ib2)

             call set_fs2_tmp
             zsum = 0.0d0
             Call part_of_USPP_correction
             call add_zsum_to_RhoTilde

          End do
       End do
    End do
! ------------------------------- end -----------
    call unset_work_arrays
    call tstatc0_end(id_sname)

  contains

    subroutine set_work_arrays
      allocate( zsum( Num_q_Points, kg1 ) ); zsum = 0.0d0
      allocate( fs1_tmp(nlmta,kimg)); fs1_tmp = 0.0d0
      allocate( fs2_tmp(nlmta,kimg)); fs2_tmp = 0.0d0
    end subroutine set_work_arrays

    subroutine unset_work_arrays
      deallocate( fs1_tmp, fs2_tmp )
      deallocate( zsum )
    end subroutine unset_work_arrays

    subroutine set_fs1_tmp
      integer :: jb1

      if ( map_e(ib1) == myrank_e ) then
         jb1 = map_z_k(ib1)
         fs1_tmp( 1:nlmta,1:kimg ) = fsval_k( jb1, 1:nlmta, ik, 1:kimg )
      endif
! --
      if ( npes > 1 ) then
         call mpi_bcast( fs1_tmp, nlmta*kimg, MPI_DOUBLE_PRECISION, &
              &          map_e(ib1), MPI_CommGroup, ierr )
         if ( ierr /= 0 ) write(*,*) 'MPI error ', mype, ierr
      endif
    end subroutine set_fs1_tmp

    subroutine set_fs2_tmp
      fs2_tmp( 1:nlmta,1:kimg ) = fsval_kmq( jb2, 1:nlmta, ik, 1:kimg )
    end subroutine set_fs2_tmp

    subroutine part_of_USPP_correction
      integer       :: ia, it, mdvdb
      integer       :: ig, nq
      integer       :: lmt1, lmt2, u, v
      real(kind=DP) :: ctmp_r, ctmp_i, ph
      complex(kind=CMPLDP) :: z1, z2, ztmp
! -------------
      Do ia=1, natm
         it = ityp(ia)
         mdvdb = m_PP_include_vanderbilt_pot(it)
         if ( mdvdb == SKIP ) cycle

         Do lmt1=1,ilmt(it)
            u = lmta(lmt1,ia)

            Do lmt2=1,ilmt(it)
               v = lmta(lmt2,ia)

               if ( kimg==1 ) then
                  ctmp_r = fs1_tmp( u,1 ) *fs2_tmp( v,1 )
                  ctmp_i = 0.0d0
               else
                  ctmp_r = fs1_tmp( u,1 ) *fs2_tmp( v,1 ) &
                       &  +fs1_tmp( u,2 ) *fs2_tmp( v,2 )
                  ctmp_i = fs1_tmp( u,1 ) *fs2_tmp( v,2 ) &
                       &  -fs1_tmp( u,2 ) *fs2_tmp( v,1 )
               endif
               Do nq=1, Num_q_Points
                  Do ig=1, kg1
                     ph = vec_q_plus_G_LR(ig,1,nq) *cps(ia,1) &
                          & +vec_q_plus_G_LR(ig,2,nq)*cps(ia,2) &
                          & +vec_q_plus_G_LR(ig,3,nq)*cps(ia,3)
!                     ph = dot_product( vec_q_plus_G_LR(ig,:,nq), &
!                       & cps(ia,:) )
                     z1 = dcmplx( 0.0d0, ph )
                     z2 = ftqval( lmt1,lmt2,it,nq,ig )
!                                                conjg wo toruhitsuyou ?
                     ztmp = exp(-z1) *z2 *iwei(ia)
!!!                     ztmp = z2 *iwei(ia)
!                     if ( ig==1 ) then
!                        write(*,*) 'ao ', exp(z1), z2
!                     endif
                     zsum(nq,ig) = zsum(nq,ig) + ztmp *dcmplx( ctmp_r, ctmp_i )
                  End do
               End do

            End do
         End do
      End do
    end subroutine part_of_USPP_correction

    subroutine add_zsum_to_RhoTilde           ! nq == 1 is supported
      integer :: ic1, nq, ig

      ic1 = ( nk_in_the_process -1 )/ (nspin*(Num_q_Points+1))
      Do nq=1, Num_q_Points
!!!                Do ig=1, nmax_G
         Do ig=1, kg1
            RhoTilde( ig, ib1,jb2, ik+ic1*nspin ) = &
                 & RhoTilde( ig, ib1,jb2, ik+ic1*nspin ) &
                 & + zsum(nq,ig)
         End do
      End do
    end subroutine add_zsum_to_RhoTilde

  end Subroutine Add_HardPart_General

  subroutine Delete_G0_Contrib
    RhoTilde(1,:,:,:) = 0.0d0
  end subroutine Delete_G0_Contrib

!------------------------------------------------------------------
!!
!!!            Calc RhoTilde in the case of LongWaveLimit   ( soft part )
!!
!-------------------------------------------------------------------
  subroutine Add_SoftPart_LWLimit
    Real(kind=DP), allocatable :: wfn_k_Rspace(:,:)
    Real(kind=DP), allocatable :: chgr(:), chgi(:)
    Real(kind=DP), allocatable :: wf1_tmp(:), wf2_tmp(:)

    Real(kind=DP) :: occ1, occ2, c1

    integer ngrid, nffth
    integer ispin, ik
    integer ib1, ib2, jb2

    integer :: id_sname = -1
! ---------------------------------- start ---------
    call tstatc0_begin('Add_SoftPart_LWLimit ', id_sname)
!
    ngrid = product(fft_box_size_WF(1:3,1));   nffth = nfft/2
    call m_FFT_alloc_WF_work()
    call set_work_arrays
! ---------------------------------- main --------------------
! ------------------ < Phi_ib1 | exp(igr) | Phi_ib2 > --------
! ------------------------------------------------------------
    Do ik=1, kv3
       !
       Do ib1=1, neg                  ! unoccupied Band

          if ( nrd_efermi==1 ) then
             occ1 = occup_k( ib1,ik )
             if ( occ1 > DELTA ) cycle
          endif
          call set_wf1_tmp

          Do ib2=1, neg                ! occupied Band

             if ( nrd_efermi==1 ) then
                occ2 = occup_k( ib2,ik )
                if ( occ2 < DELTA ) cycle
             end if

             if ( map_e(ib2) /= myrank_e ) cycle
             jb2 = map_z_k(ib2)

             call set_wf2_tmp
             call set_RhoTilde_in_Gspace
!
          End do
       End do
    End do
! --------------------------------- end -------------
    call unset_work_arrays
    Call m_FFT_dealloc_WF_work()
    call tstatc0_end(id_sname)

  contains
!
    subroutine set_work_arrays
      allocate( wfn_k_Rspace(nfft,np_e) ); wfn_k_Rspace = 0.0d0
      allocate( chgr(nffth)); allocate(chgi(nffth)); chgr = 0.0d0; chgi = 0.0d0
      allocate( wf1_tmp(nfft) ); wf1_tmp = 0.0d0
      allocate( wf2_tmp(nfft) ); wf2_tmp = 0.0d0

      Do ik=1, kv3
         Do ib1=1, np_e
            Call  Get_WF_in_Rspace( ik, ib1, wfn_k, wf1_tmp, &
                 &                  nbase_k, iba_k, igf_k, map_z_k )
            wfn_k_Rspace(1:nfft,ib1) = wf1_tmp(1:nfft)
         End do
      End do
      wfn_k_Rspace = wfn_k_Rspace / sqrt(univol)
    end subroutine set_work_arrays
!
    subroutine unset_work_arrays
      deallocate( chgr )
      if ( allocated( chgi) ) deallocate( chgi )
      deallocate( wfn_k_Rspace )
      deallocate( wf1_tmp, wf2_tmp )
    end subroutine unset_work_arrays

    subroutine set_wf1_tmp
      integer :: jb1

      wf1_tmp = 0.0d0
      if ( map_e(ib1) == myrank_e ) then
         jb1 = map_z_k(ib1)
         wf1_tmp(1:nfft) = wfn_k_Rspace(1:nfft,jb1)
      endif
      if ( npes > 1 ) then
         call mpi_bcast( wf1_tmp, nfft, MPI_DOUBLE_PRECISION, &
              &         map_e(ib1), MPI_CommGroup, ierr )
         if ( ierr /= 0 ) write(*,*) 'MPI error ', mype, ierr
      endif
    end subroutine set_wf1_tmp

    subroutine set_wf2_tmp
      wf2_tmp(1:nfft) = wfn_k_Rspace(1:nfft,jb2)
    end subroutine set_wf2_tmp

    subroutine set_RhoTilde_in_Gspace
      integer i, ip, i1
      integer ic1

      if ( k_symmetry(ik) == GAMMA ) then
         do i=1,nffth
            ip = 2*i-1
            chgr(i) = wf1_tmp(ip)  *wf2_tmp(ip)
            chgi(i) =-wf1_tmp(ip+1)*wf2_tmp(ip)
         end do
      else
         do i=1,nffth
            ip = 2*i-1
            chgr(i) = wf1_tmp(ip)*wf2_tmp(ip)  +wf1_tmp(ip+1)*wf2_tmp(ip+1)
            chgi(i) = wf1_tmp(ip)*wf2_tmp(ip+1)-wf1_tmp(ip+1)*wf2_tmp(ip)
         end do
      end if
      ! ---
      wf2_tmp = 0.0d0
      do i=1,nffth
         ip = 2*i-1
         wf2_tmp(ip)   =  chgr(i);  wf2_tmp(ip+1) =  chgi(i)
      end do
! ----------------------------- FFT : Rspace -> G-space ---------
      call m_FFT_WF( ELECTRON,nfout, wf2_tmp, DIRECT,OFF )
      wf2_tmp = wf2_tmp / dble(ngrid) *univol
! ----------------------------------------------------------------
      ic1 = ( nk_in_the_process -1 )/ ( nspin*(Num_q_Points+1) )
      Do i=1, kg1
         if(nbase_k(i,ik)<1) cycle
         i1 = igf_k( nbase_k(i,ik) )
         ip = 2*i1 - 1
         RhoTilde( i, ib1, jb2, ik + ic1*nspin ) &
              & = dcmplx( wf2_tmp(ip), wf2_tmp(ip+1) )
      End do
    end subroutine set_RhoTilde_in_Gspace

  end subroutine Add_SoftPart_LWLimit

!------------------------------------------------------------------
!!
!!!            Calc RhoTilde in the case of LongWaveLimit   ( hard part )
!!
!-------------------------------------------------------------------
  subroutine Add_HardPart_LWLimit
    Real(kind=DP), allocatable :: fs1_tmp( :,: ), fs2_tmp(:,:)
    Real(kind=DP), allocatable :: fs2_tmp_mpi( :,: )
    complex(kind=CMPLDP),   allocatable  :: zsum(:,:)

    real(kind=DP) :: occ1, occ2

    integer :: ik, ib1, ib2
    integer :: jb2
    integer :: id_sname = -1
! -------------------------------- start ------------
    call tstatc0_begin( 'Add_HardPart_LWLimit', id_sname )
    call set_work_arrays
!------------------------------------------------------
!------------------------------- main -----------------
!------------------------------------------------------
    Do ik=1, kv3

       Do ib1= 1, neg

          if ( nrd_efermi==1 ) then
             occ1 = occup_k( ib1,ik )
             if ( occ1 > DELTA ) cycle
          endif
          call set_fs1_tmp

          Do ib2=1, neg                ! occupied Band

             if ( nrd_efermi==1 ) then
                occ2 = occup_k( ib2,ik )
                if ( occ2 < DELTA ) cycle
             end if

             if ( map_e(ib2) /= myrank_e ) cycle
             jb2 = map_z_k(ib2)

             call set_fs2_tmp
             zsum = 0.0d0
             Call part_of_USPP_correction
             call add_zsum_to_RhoTilde

          End do
       End do
    End do
! ----------------------------------end -------------------
    call unset_work_arrays
!    stop

    call tstatc0_end(id_sname)

  contains

    subroutine set_work_arrays
      allocate( zsum( Num_q_Points, kg1 ) ); zsum = 0.0d0
      allocate( fs1_tmp(nlmta,kimg)); fs1_tmp = 0.0d0
      allocate( fs2_tmp(nlmta,kimg)); fs2_tmp = 0.0d0
    end subroutine set_work_arrays

    subroutine unset_work_arrays
      if ( allocated(fs1_tmp) ) deallocate( fs1_tmp )
      if ( allocated(fs2_tmp) ) deallocate( fs2_tmp )
      if ( allocated(zsum   ) ) deallocate( zsum )
    end subroutine unset_work_arrays

    subroutine set_fs1_tmp
      integer :: jb1

      if ( map_e(ib1) == myrank_e ) then
         jb1 = map_z_k(ib1)
         fs1_tmp( 1:nlmta,1:kimg ) = fsval_k( jb1, 1:nlmta, ik, 1:kimg )
      endif
! --
      if ( npes > 1 ) then
         call mpi_bcast( fs1_tmp, nlmta*kimg, MPI_DOUBLE_PRECISION, &
              &          map_e(ib1), MPI_CommGroup, ierr )
         if ( ierr /= 0 ) write(*,*) 'MPI error ', mype, ierr
      endif
    end subroutine set_fs1_tmp

    subroutine set_fs2_tmp
      fs2_tmp( 1:nlmta,1:kimg ) = fsval_k( jb2, 1:nlmta, ik, 1:kimg )
    end subroutine set_fs2_tmp

    subroutine part_of_USPP_correction
      integer       :: ia, it, mdvdb
      integer       :: nq, ig
      integer       :: lmt1, lmt2, u, v
      real(kind=DP) :: ctmp_r, ctmp_i, ph
      complex(kind=CMPLDP) :: z1, z2, ztmp
! -------------
      Do ia=1, natm
         it = ityp(ia)
         mdvdb = m_PP_include_vanderbilt_pot(it)
         if ( mdvdb == SKIP ) cycle

         Do lmt1=1,ilmt(it)
            u = lmta(lmt1,ia)

            Do lmt2=1,ilmt(it)
               v = lmta(lmt2,ia)

               if ( kimg==1 ) then
                  ctmp_r = fs1_tmp( u,1 ) *fs2_tmp( v,1 )
                  ctmp_i = 0.0d0
               else
                  ctmp_r = fs1_tmp( u,1 ) *fs2_tmp( v,1 ) &
                       &  +fs1_tmp( u,2 ) *fs2_tmp( v,2 )
                  ctmp_i = fs1_tmp( u,1 ) *fs2_tmp( v,2 ) &
                       &  -fs1_tmp( u,2 ) *fs2_tmp( v,1 )
               endif

!               write(*,*) 'ctmp_r, i =', ctmp_r, ctmp_i, iwei(ia)
!
               Do nq=1, Num_q_Points
                  Do ig=1, nmax_G_LR
!!!                  Do ig=1, min( nmax_G_LR, kg0 )

                     ph = vec_q_plus_G_LR(ig,1,nq) *cps(ia,1) &
                          & +vec_q_plus_G_LR(ig,2,nq)*cps(ia,2) &
                          & +vec_q_plus_G_LR(ig,3,nq)*cps(ia,3)

!                     ph = dot_product( vec_q_plus_G_LR(ig,:,nq), &
!                       & cps(ia,:) )
                     z1 = cmplx( 0.0d0, ph )
                     z2 = ftqval( lmt1,lmt2,it,nq,ig )
                     ztmp = exp(-z1) *z2 *iwei(ia)
!                     write(*,*) 'ztmp = ', ztmp

!!!                     ztmp = z2 *iwei(ia)
!                     if ( ig==1 ) then
!                        write(*,*) 'ao ', ztmp, ctmp_r
!                     endif
                     zsum(nq,ig) = zsum(nq,ig) + ztmp *dcmplx( ctmp_r, ctmp_i )
                  End do
               End do

            End do
         End do

      End do

    end subroutine part_of_USPP_correction

    subroutine add_zsum_to_RhoTilde           ! nq == 1 is supported
      integer :: ic1, nq, ig

      ic1 = ( nk_in_the_process -1 )/ nspin
      Do nq=1, Num_q_Points
!!!                Do ig=1, nmax_G
         Do ig=1, kg1
            RhoTilde( ig, ib1,jb2, ik+ic1*nspin ) = &
                 & RhoTilde( ig, ib1,jb2, ik+ic1*nspin ) &
                 & + zsum(nq,ig)
         End do
      End do
    end subroutine add_zsum_to_RhoTilde

  end Subroutine Add_HardPart_LWLimit

!------------------------------------------------------------------
!!
!!!            Calc RhoTilde in the case of LongWaveLimit   ( soft, correction )
!!
!-------------------------------------------------------------------
  subroutine Add_Correction_To_G0
    Real(kind=DP), allocatable :: wf1_tmp(:,:)
    Real(kind=DP), allocatable :: wf2_tmp(:,:)
    Real(kind=DP), allocatable :: k_plus_G(:,:)

    real(kind=DP) :: csum_r(3), csum_i(3)
    integer       :: ik
    integer       :: ib1, ib2, jb2
    real(kind=DP) :: occ1, occ2
! -------------------------------------- start ---------
    call set_work_arrays
! ----------------------------------------- main -------
    Do ik=1, kv3
       call set_vectors_kplusG

       Do ib1=1, neg                  ! unoccupied Band

          if ( nrd_efermi ==1 ) then
             occ1 = occup_k( ib1,ik )
             if ( occ1 > DELTA ) cycle
          endif

          call set_wf1_tmp

          Do ib2=1, neg                ! occupied Band

             if ( nrd_efermi ==1 ) then
                occ2 = occup_k( ib2,ik )
                if ( occ2 < DELTA ) cycle
             end if
             if ( ib1 == ib2 ) cycle

             if ( map_e(ib2) /= myrank_e ) cycle
             jb2 = map_z(ib2)

             call set_wf2_tmp
!
             csum_r = 0.0d0; csum_i = 0.0d0
             call part_of_local_contrib
             call add_csum_to_RhoTilde

          End do
       End do

    End do
! ---------------------------- end ---------------
    call unset_work_arrays
  contains

    subroutine set_work_arrays
      allocate( k_plus_G( kg1,3 ) ); k_plus_G = 0.0d0
      allocate( wf1_tmp(kg1,kimg) ); wf1_tmp = 0.0d0
      allocate( wf2_tmp(kg1,kimg) ); wf2_tmp = 0.0d0
    end subroutine set_work_arrays

    subroutine unset_work_arrays
      deallocate( k_plus_G )
      deallocate( wf1_tmp, wf2_tmp )
    end subroutine unset_work_arrays

    subroutine set_vectors_kplusG
      integer :: i, ip
      Real(kind=DP) :: ga, gb, gc

      Do i = 1, iba_k(ik)
         ip = nbase_k(i,ik)
         ga = vkxyz(ik,1,BUCS) + real(ngabc(ip,1),kind=DP)
         gb = vkxyz(ik,2,BUCS) + real(ngabc(ip,2),kind=DP)
         gc = vkxyz(ik,3,BUCS) + real(ngabc(ip,3),kind=DP)

         k_plus_G(i,1)  = rltv(1,1)*ga + rltv(1,2)*gb + rltv(1,3)*gc
         k_plus_G(i,2)  = rltv(2,1)*ga + rltv(2,2)*gb + rltv(2,3)*gc
         k_plus_G(i,3)  = rltv(3,1)*ga + rltv(3,2)*gb + rltv(3,3)*gc
      End do
    end subroutine set_vectors_kplusG

    subroutine set_wf1_tmp
      integer :: jb1
!
      if ( map_e(ib1) == myrank_e ) then
         jb1 = map_z(ib1)
         wf1_tmp( 1:kg1, 1:kimg ) = wfn_k( 1:kg1, jb1, ik, 1:kimg )
      endif
      if ( npes > 1 ) then
         call mpi_bcast( wf1_tmp, kg1*kimg, MPI_DOUBLE_PRECISION, &
              &         map_e(ib1), MPI_CommGroup, ierr )
         if ( ierr /= 0 ) write(*,*) 'MPI error ', mype, ierr
      endif
    end subroutine set_wf1_tmp

    subroutine set_wf2_tmp
      wf2_tmp(1:kg1,1:kimg) = wfn_k(1:kg1,jb2,ik,1:kimg)
    end subroutine set_wf2_tmp

    subroutine part_of_local_contrib
      integer       :: i
      real(kind=DP) :: ctmp_r, ctmp_i
! ----------------------------------
      if ( kimg <= 1 ) then
         Do i = 1, iba(ik)
            ctmp_r = wf1_tmp( i,1 ) *wf2_tmp( i,1 )
            csum_r(1) = csum_r(1) + ctmp_r *k_plus_G(i,1)
            csum_r(2) = csum_r(2) + ctmp_r *k_plus_G(i,2)
            csum_r(3) = csum_r(3) + ctmp_r *k_plus_G(i,3)
         End do
      else
         Do i=1, iba(ik)
            ctmp_r = wf1_tmp( i,1 ) *wf2_tmp( i,1 ) &
                 &  +wf1_tmp( i,2 ) *wf2_tmp( i,2 )
            ctmp_i = wf1_tmp( i,1 ) *wf2_tmp( i,2 ) &
                 &  -wf1_tmp( i,2 ) *wf2_tmp( i,1 )

            csum_r(1) = csum_r(1) + ctmp_r *k_plus_G(i,1)
            csum_r(2) = csum_r(2) + ctmp_r *k_plus_G(i,2)
            csum_r(3) = csum_r(3) + ctmp_r *k_plus_G(i,3)
            !
            csum_i(1) = csum_i(1) + ctmp_i *k_plus_G(i,1)
            csum_i(2) = csum_i(2) + ctmp_i *k_plus_G(i,2)
            csum_i(3) = csum_i(3) + ctmp_i *k_plus_G(i,3)
         end do
      end if
    end subroutine part_of_local_contrib

    subroutine add_csum_to_RhoTilde
      integer :: ic1
      Complex(kind=CMPLDP) :: ztmp(3), zz1
      Real(kind=DP) :: ebi, ebj, ediff

      ic1 = ( nk_in_the_process -1 ) / nspin
!
      ztmp(1) = dcmplx( csum_r(1), csum_i(1) ) *vqxyz( 1,1,CARTS )
      ztmp(2) = dcmplx( csum_r(2), csum_i(2) ) *vqxyz( 1,2,CARTS )
      ztmp(3) = dcmplx( csum_r(3), csum_i(3) ) *vqxyz( 1,3,CARTS )
      !
      ebi = eko_k( ib1,ik ); ebj = eko_k( ib2,ik )
      ediff = ebi - ebj
      !
      zz1 = ( ztmp(1)+ztmp(2)+ztmp(3) ) / ediff
      !
! --------------- test-
      if ( mype == 0 ) then
	write(*,'(I5,2I3,6F12.7)') ik+ic1, ib1, ib2, vkxyz(ik,1,CARTS), &
        & vkxyz(ik,2,CARTS), vkxyz(ik,3,CARTS), csum_r(1), csum_r(2), csum_r(3)
      endif
! ---------------

      RhoTilde( 1, ib1, jb2, ik + ic1*nspin ) = &
           & RhoTilde( 1, ib1, jb2, ik + ic1*nspin )  + zz1
    end subroutine add_csum_to_RhoTilde

  end subroutine Add_Correction_To_G0

!------------------------------------------------------------------
!!
!!!            Correction by Kageshima Shiraishi
!!
!-------------------------------------------------------------------
  subroutine Add_KageShira_To_G0
    real(kind=DP) :: occ1, occ2
    real(kind=DP) :: csum_r(3), csum_i(3)

    integer       :: i, ip, ik
    integer       :: ib1, ib2, jb2
! -------------------------------- start -------------
    nppcorr = 2
!    call check_PP(nfout)
    if ( nppcorr>0  ) call calc_ptrans_ek
    if ( nppcorr==2 ) call calc_ptrans_TM_PP_ek
! -------------------------------- main -----------
    Do ik=1, kv3

       Do ib1=1, neg                  ! unoccupied Band

          if ( nrd_efermi == ON ) then
             occ1 = occup_k( ib1,ik )
             if ( occ1 > DELTA ) cycle
          endif

          Do ib2=1, neg                ! occupied Band
             if ( ib1 == ib2 ) cycle

             if ( map_e(ib2) /= myrank_e ) cycle
             jb2 = map_z(ib2)
!
             csum_r = 0.0d0; csum_i = 0.0d0
             call part_of_nonlocal_contrib
             call add_csum_to_RhoTilde
          End do
       End do
    End do
!
  contains

    subroutine part_of_nonlocal_contrib
      csum_r(1) = ptrans( ik, ib1, ib2, 1, 1 )
      csum_r(2) = ptrans( ik, ib1, ib2, 2, 1 )
      csum_r(3) = ptrans( ik, ib1, ib2, 3, 1 )

      csum_i(1) = ptrans( ik, ib1, ib2, 1, 2 )
      csum_i(2) = ptrans( ik, ib1, ib2, 2, 2 )
      csum_i(3) = ptrans( ik, ib1, ib2, 3, 2 )
    end subroutine part_of_nonlocal_contrib

    subroutine add_csum_to_RhoTilde
      integer :: ic1
      Complex(kind=CMPLDP) :: ztmp(3), zz1
      Real(kind=DP) :: ebi, ebj, ediff

      ic1 = ( nk_in_the_process -1 )/  nspin
!
      ztmp(1) = dcmplx( csum_r(1), csum_i(1) ) *vqxyz( 1,1,CARTS )
      ztmp(2) = dcmplx( csum_r(2), csum_i(2) ) *vqxyz( 1,2,CARTS )
      ztmp(3) = dcmplx( csum_r(3), csum_i(3) ) *vqxyz( 1,3,CARTS )
!
      ebi = eko_k( ib1,ik ); ebj = eko_k( ib2,ik )
      ediff = ebi - ebj

      zz1 = ( ztmp(1)+ztmp(2)+ztmp(3) ) / ediff

      RhoTilde( 1, ib1, jb2, ik + ic1*nspin ) = &
           & + RhoTilde( 1, ib1, jb2, ik + ic1*nspin ) + zz1
    end subroutine add_csum_to_RhoTilde

  end subroutine Add_KageShira_To_G0

!---------------------------------------------------------------------
!!
!!!           calc Fermi level
!!                            from WriteDownData_onto_Files_ek.f90
!----------------------------------------------------------------------

  subroutine Calc_FermiLevel_ek(nrd_efermi)
!                        calculate fermi energy (efermi) if nrd_efermi=0
    integer             :: way_of_smearing
    integer, intent(in) :: nrd_efermi
    real(kind=DP)            :: efermi1
! ----------------------------------------------
    if ( nrd_efermi== 1 ) then
       efermi1 = efermi;     efermi = 0.0d0
    end if
    way_of_smearing = m_CtrlP_way_of_smearing()

    if ( way_of_smearing==PARABOLIC ) then
       call Fermi_parabolic_ek(nfout)
    else if ( way_of_smearing == TETRAHEDRON ) then
       !       write(*,*) ' ------------------ Under construction ---------'
       !       stop
       call Fermi_tetra_ek(nfout)
    end if
    !
    if ( nrd_efermi==0 ) then
       if (printable) write(nfout,10) efermi
    else
       if (printable) write(nfout,20) efermi1-efermi
       efermi=efermi1
    end if
10  format(1x,"!*--- efermi = ",f10.5)
20  format(1x,"!*--- difference between read and calculated efermi = ",f10.5)
  end subroutine Calc_FermiLevel_ek

  subroutine Fermi_parabolic_ek(nfout)
!
!   Coded from <m_ESoc_fermi_parabolic>
!                      by T. Yamasaki (FUJITSU LABORATORIES Ltd.), 28th Jun. 2003
    integer, intent(in) :: nfout

    integer             :: jcount
    real(kind=DP)       :: emin, emax, e1, e2, tot
    real(kind=DP), parameter :: DELTA_TOTCH = 1.d-10
    integer,       parameter :: MAXITR = 100000
    integer :: id_sname = -1
!! UVSOR
    integer             :: ik, ie
    integer             :: ik1, ik2, k, ispin, nq, nq_max
    real(kind=DP),allocatable,dimension(:,:) :: occup_mpi

    allocate(occup_mpi(neg,kv3_ek)); occup_mpi =0.0d0
    !! UVSOR

    call tstatc0_begin('Fermi_parabolic_ek ', id_sname)
    call check_totch(nfout)

    emin = minval(eko_ek);     emax = maxval(eko_ek)
    efermi = emax
    e1 = emin - dabs(DELTA_FermiSearchRange)
    e2 = emax + dabs(DELTA_FermiSearchRange)

    jcount = 1
    occupied_ch_equals_totch : do

       call get_tot_eps()             ! -(contained here) ->tot
!          ~~~~~~~~~~~~~~~~~~~~
       if(jcount == 1 .and. tot < totch) &
            call wd_fermi_error1(nfout,emin,emax,tot,totch) ! -(b_Fermi)
       if(dabs(tot - totch) < DELTA_TOTCH) exit occupied_ch_equals_totch
       if( tot < totch) then
          e1 = efermi; efermi = efermi + (e2-efermi)/2
       else
          e2 = efermi; efermi = efermi + (e1-efermi)/2
       end if
       jcount = jcount + 1
       if(jcount > MAXITR) then
          call wd_fermi_error2 &   !-(b_Fermi)
               & (nfout,e1,e2,efermi,emin,emax,tot,totch,neg,MAXITR)
       end if
    end do occupied_ch_equals_totch
    if(iprioccup >= 2) &
         &write(nfout,'( " efermi = ", f10.4, " tot = ", f10.4)') efermi, tot

    call check_if_metalic(nfout,metalic_system)
    if(iprioccup >= 2 .and. .not.metalic_system) &
         &write(nfout,'( " The highest occupied band enrgy = ", f10.4)') vbm

    call tstatc0_end(id_sname)
!!  -----------------------------------------
    if ( sw_LongWaveLimit ==ON ) then
       nq_max = 0
    else
       nq_max = Num_q_Points
    endif

    Do k=1, kv3_ek / ( nspin*(nq_max+1 ))
       Do ispin=1, nspin, af+1
          ik1 = ( nq_max +1 )*nspin*(k-1) + ispin
          Do nq=0, nq_max
             ik2 = ik1 + nspin *nq

!             Do ie = 1, neg                                   ! MPI
!                if(map_e(ie) == myrank_e) then                ! MPI
!                   occup_l_ek(map_z(ie),ik2) = occup_mpi(ie,ik2)! MPI
!                end if                                        ! MPI
!             end do                                           ! MPI
! ----------
             occup_lkt_ek(1:neg, ik2) = occup_mpi(1:neg,ik2)
! ------------
          End do
       end do                                              ! MPI
   End do
   deallocate(occup_mpi)
   !! UVSOR

 contains
    subroutine get_tot_eps
      integer       :: k, i, ik1, ik2, ispin, nq
      integer       :: nq_max

      real(kind=DP) :: wspin = 1.d0, e, dos, weight, totw

      tot = 0.d0
! ------------------------
      if ( sw_LongWaveLimit ==ON ) then
         nq_max = 0
      else
         nq_max = Num_q_Points
      endif
! ---
      Do k=1, kv3_ek / ( nspin*(nq_max+1 ))
         Do ispin=1, nspin, af+1
            ik1 = ( nq_max +1 )*nspin*(k-1) + ispin
            Do nq=0, nq_max
               ik2 = ik1 + nspin *nq

               Do i = 1, neg
                  e = eko_ek(i,ik2)

                  call width2(e,efermi,width,dos,weight)  ! -(b_Fermi)
                  totw = weight *wspin *kv3_ek *qwgt_ek(k)
                  !! LR-TDDDFT
                  totw = totw / dble( nq_max + 1 )
                  occup_mpi(i,ik2) = totw
                  !            write(*,*) 'occ ' , i, k, totw
                  !!
                  tot = tot + 2*totw
               End do
            End do
         End do
      End do

      tot = tot/kv3_ek * (af+1)
      if(iprioccup >= 2) write(nfout,'(" efermi, tot = ",d20.12,d20.12," jcount = ",i7)') efermi,tot, jcount
!      write(*,*) 'aaa ', efermi, tot, jcount
!      stop
    end subroutine get_tot_eps

  end subroutine Fermi_parabolic_ek

  subroutine Fermi_tetra_ek(nfout)
    integer, intent(in) :: nfout

#ifndef NO_TETRAHEDRON
    real(kind=DP), parameter :: delta = 1.d-12
    integer, parameter       :: idim = 3
    integer        :: neig,nengy,ispin,ip2,ik,instts1,ikee,nxx,nyy,nzz, ib, jb
    real(kind=DP)  :: efermi2,eval,totind, et
    real(kind=DP), pointer, dimension(:,:,:,:) :: eig2
    real(kind=DP), pointer, dimension(:) :: eawk,cdwk,cswk,cdos,cind,valud
    integer             :: id_sname = -1

    integer        :: ip_mpi, ieig
    integer, allocatable, dimension(:,:,:) :: neordr_ek
    real(kind=DP), pointer, dimension(:,:,:,:) :: occup2

    integer nq, ik_kt1, ik_kt2, nq_max
! --------------------------
    call tstatc0_begin('m_ESoc_fermi_tetra_ek ', id_sname)

    if ( sw_LongWaveLimit ==ON ) then
       nq_max = 0
    else
       nq_max = Num_q_Points
    endif

    allocate( eig2( np2, neg, nspin, 0:nq_max )); eig2 = 0.d0
    allocate( eawk(np0) ); eawk = 0.d0
    allocate( cdwk(np0) ); cdwk = 0.d0
    allocate( cswk(np0) ); cswk = 0.d0
    allocate( cdos(np2*neg) ); cdos = 0.d0
    allocate( cind(np2*neg) ); cind = 0.d0
    allocate( valud(nspin) ); valud = 0.d0
!! UVSOR
    allocate( neordr_ek( neg, kv3_ek, 0:nq_max ) ); neordr_ek = 0
    allocate( occup2( neg, np2, nspin, 0:nq_max )); occup2 = 0.0d0
!! UVSOR

    call check_totch(nfout)           ! totch is checked

    nxx = nxyz_tetra(1)
    nyy = nxyz_tetra(2)
    nzz = nxyz_tetra(3)
    neig=neg
    nengy=0
!
    Do nq=0, nq_max
       Do ispin=1,nspin
          Do ip2=1,np2
             ik = nspin*(ip2-1)+ispin

             ik_kt1 = ( nq_max +1)*nspin*(ip2-1) + ispin
             ik_kt2 = ik_kt1 + nspin *nq

             neordr_ek(1:neg,ik,nq) = (/(ib,ib = 1, neg)/)

             Do ib=1,neg
                eig2(ip2,ib,ispin,nq) = eko_ek(ib,ik_kt2)
             End Do
             Do ib = 1,neg-1
                Do jb = ib+1, neg
                   if ( eig2(ip2,jb,ispin,nq) < eig2(ip2,ib,ispin,nq)-delta ) then
                      et = eig2(ip2,ib,ispin,nq)
                      eig2(ip2,ib,ispin,nq) = eig2(ip2,jb,ispin,nq)
                      eig2(ip2,jb,ispin,nq) = et

                      neordr_ek(ib,ik,nq) =jb
                      neordr_ek(jb,ik,nq) =ib

                   end if
                End do
             End do
          Enddo
       Enddo
    End do
!
    if(iprioccup>=2) write(nfout,*) ' === tetrahedron method', &
                &  ' for k-space integration ==='
    if(iprioccup>=2) then
       write(nfout,'(" --- eig2 ---")')
       do ip2 = 1, np2
          write(nfout,'(" (ip2 = ",i8,")",99(/,8d12.4))') ip2,(eig2(ip2,ib,ispin,0),ib=1,neg)
       end do
    end if

    efermi2 = efermi
    call fermi1(nxx,nyy,nzz,np2,np2,neig,neg,nspin, &
             &  eig2(1,1,1,0), ip20, np0,totch, efermi2,eval,valud, &
             &  iwt,ip2cub,iprioccup)
    if(iprioccup>=2) write(nfout,*) 'eval=',eval

! ------------------------------ nq == 0 ------ evaluated @ q=0 --
    nq = 0
    Do ispin=1,nspin
       instts1 = 0
       call nsdos3(nfout,idim,efermi2,nxx,nyy,nzz,np2,1,neig, &
            &  eig2(1,1,ispin,nq),ip20,np0,eawk,instts1,np2, &
            &  iwt,ip2cub,iprioccup)
       call nstt3i(idim,nengy,efermi2,nxx,nyy,nzz, &
            &  np2,np2,neig,eig2(1,1,ispin,nq), &
            &  ip20,np0,eawk,cdwk,cswk,np2,neig,cdos,cind,width_tetra )
       efermi=efermi2

       totind=0.d0
       Do ip2=1,np2
          Do ib=1,neig
             ikee=ip2+np2*(ib-1)

             ik = nspin * (ip2-1) + ispin
             ip_mpi = neordr_ek(ib,ik,nq)
             occup2(ip_mpi,ip2,ispin,nq) = cind(ikee)*dble(np2)

!             write(*,*) 'FFA ', cind(ikee)*dble(np2)
             totind=totind+cind(ikee)
          End do
       End do
       if(iprioccup>=2) write(nfout,*) 'for spin=',ispin, &
                  &  ' ** TOTAL CHARGE after fermi1 = ',totind
!       write(*,*) 'aaaaa', totind
    Enddo
! ---------------------------- nq > 0 --------
    Do nq=1, nq_max
       Do ispin=1,nspin
          instts1 = 0
          call nsdos3(nfout,idim,efermi2,nxx,nyy,nzz,np2,1,neig, &
               &  eig2(1,1,ispin,nq),ip20,np0,eawk,instts1,np2, &
               &  iwt,ip2cub,iprioccup)
          call nstt3i(idim,nengy,efermi2,nxx,nyy,nzz, &
               &  np2,np2,neig,eig2(1,1,ispin,nq), &
               &  ip20,np0,eawk,cdwk,cswk,np2,neig,cdos,cind,width_tetra )

          Do ip2=1,np2
             Do ib=1,neig
                ikee=ip2+np2*(ib-1)

                ik = nspin * (ip2-1) + ispin
                ip_mpi = neordr_ek(ib,ik,nq)
                occup2(ip_mpi,ip2,ispin,nq) = cind(ikee)*dble(np2)
!!                write(*,*) 'FFF ', cind(ikee)*dble(np2)
             End do
          End do
       End do
    Enddo
! -------------------------------------
!! UVSOR
    Do ieig = 1, neg
!!!!!!!       if(map_e(ieig) /= myrank_e) cycle                          ! MPI
       Do nq=0, nq_max
          Do ispin=1,nspin
             Do ip2=1,np2
                ik = nspin*(ip2-1)+ispin
                ik_kt1 = ( nq_max+1 )*nspin*(ip2-1) + ispin
                ik_kt2 = ik_kt1 + nspin *nq

!!!!!!!!                occup_l_ek(map_z(ieig),ik_kt2) = occup2(ieig,ip2,ispin,nq) ! MPI
! --------------------------
                occup_lkt_ek(ieig,ik_kt2) = occup2(ieig,ip2,ispin,nq) ! MPI
! -------------------------
             End do
          End do
       End do
    End do
    deallocate(neordr_ek)
    deallocate(occup2)
!! UVSOR

    deallocate(valud); deallocate(cind); deallocate(cdos); deallocate(cswk)
    deallocate(cdwk);  deallocate(eawk)
    deallocate(eig2)

    call check_if_metalic(nfout,metalic_system)
    if(iprioccup >= 2 .and. .not.metalic_system) &
         &write(nfout,'( " The highest occupied band enrgy = ", f10.4)') vbm

    call tstatc0_end(id_sname)
#endif
  end subroutine Fermi_tetra_ek

! ------------------------------------------------------
!! -----------------------------
!!!              other subroutines inspired from m_Epsilon_ek.f90
!! -----------------------------
! ------------------------------------------------------

  subroutine Dealloc_Ptrans_Nppc_ilocal_etc
    call Dealloc_Ptrans
    if ( nppcorr>=1 )  then
       call Dealloc_Nppc_data
       call Dealloc_ptrans_data_array_ek
    end if
    if(nppcorr==2)  call Dealloc_ilocal_l
  end subroutine Dealloc_Ptrans_Nppc_ilocal_etc

  subroutine Alloc_ptrans
    allocate( ptrans( kv3, neg,neg,3,2 ) ); ptrans = 0.0d0
  end subroutine Alloc_ptrans

  subroutine Dealloc_ptrans
    deallocate( ptrans )
  end subroutine Dealloc_ptrans

  subroutine Alloc_Nppc_data
    allocate(nppc_data(ntyp)) ; nppc_data=0
  end subroutine Alloc_Nppc_data

  subroutine Dealloc_Nppc_data
    deallocate(nppc_data)
  end subroutine Dealloc_Nppc_data

  subroutine Alloc_ilocal_l
    allocate(ilocal_l(ntyp))
  end subroutine Alloc_ilocal_l

  subroutine Dealloc_ilocal_l
    deallocate( ilocal_l )
  end subroutine Dealloc_ilocal_l

  subroutine alloc_ptrans_data_array_ek    ! allocate data array for core-repair term
    integer :: it
    ! find maximum of nppc_data(it): it=1-ntyp
    if(mype == 0) then
       if(ntyp==1) then
          mnppc=nppc_data(1)
       else
          mnppc=nppc_data(1)
          do it=1,ntyp
             if(mnppc<nppc_data(it)) mnppc=nppc_data(it)
          end do
       end if
    end if
    call mpi_bcast(mnppc,1,mpi_integer,0,MPI_CommGroup, ierr) ! MPI
    if (printable) write(nfout,'(1x,"!* mnppc = ",i3)') mnppc
    ! allocate data array
    allocate(dipole_dxyz_us(ntyp,mnppc,3)); dipole_dxyz_us=0.0d0
    allocate(dipole_tau1(ntyp,mnppc));      dipole_tau1=0
    allocate(dipole_tau2(ntyp,mnppc));      dipole_tau2=0
    allocate(phase_ylm1(ntyp,mnppc));       phase_ylm1=0
    allocate(phase_ylm2(ntyp,mnppc));       phase_ylm2=0
  end subroutine alloc_ptrans_data_array_ek

  subroutine dealloc_ptrans_data_array_ek
    deallocate( dipole_tau1, dipole_tau2 )
    deallocate( dipole_dxyz_us )
    deallocate( phase_ylm1, phase_ylm2 )
  end subroutine dealloc_ptrans_data_array_ek

  subroutine set_ppc_data(nfout)     ! set core-repair param for Kagashima-Shiraishi
!
    integer, intent(in) :: nfout
    integer             :: it, nfpp, ierr

    call Alloc_Nppc_data
    call m_Files_open_ps_files(ivan,iatomn,ntyp,ierr)
    if(mype == 0) then
       if(nppcorr==2) then
          if(printable) &
           & write(nfout,'(1x,"!* KS Transition Moment Correction for Troullier-Martin pseudopotential")')
          if(sw_use_add_proj == ON) then
             paramset=.false.
             call m_PP_make_index_lmtt_add(nfout,paramset)
          end if
       end if
    end if
    call find_ppc_data_number(nfout)
    call Alloc_ptrans_data_array_ek
    call read_ptrans_data_ek
    if (nppcorr==2) then
       call Alloc_ilocal_l
       call prepare_for_TM_PP_ek
    endif

  end subroutine set_ppc_data

  subroutine set_ppc_data_it(nfout)
    integer, intent(in) :: nfout
    integer             :: it, nfpp, ierr
    call Alloc_Nppc_data
    if(mype==0)then
       if(nppcorr==2) then
          if(printable) &
           & write(nfout,'(1x,"!* KS Transition Moment Correction for Troullier-Martin pseudopotential")')
          if(sw_use_add_proj == ON) then
             paramset=.false.
             call m_PP_make_index_lmtt_add(nfout,paramset)
          end if
       end if
    endif
    do it=1,ntyp
       call m_Files_open_ps_file(ivan,iatomn,ntyp,it,ierr)
       call find_ppc_data_number_it(nfout,it)
       call m_Files_close_ps_file(it)
    enddo
    call Alloc_ptrans_data_array_ek()
    do it=1,ntyp
       call m_Files_open_ps_file(ivan,iatomn,ntyp,it,ierr)
       call find_ppc_data_number_it(nfout,it)
       call read_ptrans_data_ek_it(it)
       call m_Files_close_ps_file(it)
    enddo
    if (nppcorr==2) then
       call Alloc_ilocal_l
       call prepare_for_TM_PP_ek
    endif
  end subroutine set_ppc_data_it

  subroutine Alloc_PP_norm_type
    allocate(PP_norm_type(ntyp))  ; PP_norm_type = 0
  end subroutine Alloc_PP_norm_type

  subroutine Dealloc_PP_norm_type
    deallocate(PP_norm_type)
  end subroutine Dealloc_PP_norm_type

  subroutine Alloc_PP_local_type
    allocate(PP_local_type(ntyp))  ; PP_local_type = 0
  end subroutine Alloc_PP_local_type

  subroutine Dealloc_PP_local_type
    deallocate(PP_local_type)
  end subroutine Dealloc_PP_local_type

  subroutine check_PP(nfout)
!
!   find pseudopotantial norm and local potential types
!
!   norm type        NC_PP : norm conserving
!                    US_PP : ultrasoft
!   local potential  BHS_POLY : BHS or polynomial type
!                    ORBITAL  : orbital local potential for Troullier-Martin pseudopotentials

    integer, intent(in)                      :: nfout
    integer                                  :: it, ntype, ltype, ilocal
    integer, allocatable, dimension(:)       :: PP_local_orbital
    character(len=1), dimension(7)           :: local_potential
    data local_potential/'s','p','d','f',' ',' ',' '/


    allocate( PP_local_orbital(ntyp) );
    PP_local_orbital = 0

    do it = 1, ntyp
       call find_norm_type(it,ntype)
       PP_norm_type(it) = ntype
       call find_local_orbital(it,ilocal)
       if(ilocal == 0) then
          PP_local_type(it) = BHS_POLY
       else
          PP_local_type(it) = ORBITAL
       end if
       PP_local_orbital(it) = ilocal
    end do


    if(printable) write(nfout,'(1x,"!* ---------- pseudopotential type ----------")')
    if(printable) write(nfout,'(1x,4x,"it",10x,"norm",12x,"local potential")')

    do it = 1, ntyp
       if(printable) then
          if(PP_norm_type(it)==NC_PP.and.PP_local_type(it)==BHS_POLY) &
         & write(nfout,'(1x,3x,i3,5x,"norm conserving",5x,"BHS or polynomial")') it
          if(PP_norm_type(it)==NC_PP.and.PP_local_type(it)==ORBITAL) &
         & write(nfout,'(1x,3x,i3,5x,"norm conserving",9x,a1,1x,"orbital")') it, local_potential(PP_local_orbital(it))
          if(PP_norm_type(it)==US_PP.and.PP_local_type(it)==BHS_POLY) &
         & write(nfout,'(1x,3x,i3,8x,"ultrasoft",8x,"BHS or polynomial")') it
          if(PP_norm_type(it)==US_PP.and.PP_local_type(it)==ORBITAL) &
         & write(nfout,'(1x,3x,i3,8x,"ultrasoft",12x,a1,1x,"orbital")') it, local_potential(PP_local_orbital(it))
       end if
    end do
!
!   check transition moment option and stop if it is wrong
!
    if(n_check_ts/=0) then
       call check_ts_option
    else
       if(printable) write(nfout,'(/1x,"!* transition moment option is not checked",/)')
    end if

    deallocate(PP_local_orbital)

  contains
    subroutine find_norm_type(it,ntype)
      implicit none
      integer, intent(in)  :: it
      integer, intent(out) :: ntype
      integer              :: lmt1, lmt2
      real(kind=DP)        :: deficit_ch
      ntype = NC_PP
      do lmt1 = 1, ilmt(it)
         do lmt2 = lmt1, ilmt(it)
            deficit_ch=dabs(q(lmt1,lmt2,it))
            if(deficit_ch>0.0d0) then
               ntype=US_PP
               exit
            end if
         end do
         if(ntype==US_PP) exit
      end do
    end subroutine find_norm_type

    subroutine find_local_orbital(it,ilocal)
      implicit none
      integer, intent(in)  :: it
      integer, intent(out) :: ilocal
      integer :: lmt, il1, il2
      ilocal = 0
      do il1=1,lpsmax(it)
         do lmt=1,ilmt(it)
            il2=ltp(lmt,it)
            if(il1==il2) exit
         end do
         if(il1/=il2) then
            ilocal=il1
         end if
      end do
    end subroutine find_local_orbital

    subroutine check_ts_option
      implicit none
      integer :: it
      if(icond==3) return
      if(printable) write(nfout,'(/1x,"!* enter check of transition moment option")')
      if(nonlocal==0.and.nppcorr==0) then
         if(printable) write(nfout,'(1x,"!* local transition moment : check is skipped",/)')
         return
      end if

      if(nonlocal/=0) then
         do it = 1, ntyp
            if(printable.and.PP_norm_type(it)/=NC_PP) then
               write(nfout,'(1x,"!* A pseudopotential in use is the ultra-soft type")')
               write(nfout,'(1x,"!* Read and Needs method cannot be used for the ultra-soft type.   LinearResponse Calc. stop")')
               stop
            end if
         end do
      end if
      if(nppcorr/=0) then
         do it = 1, ntyp
            if(PP_local_type(it)/=BHS_POLY) then
               if(printable.and.sw_use_add_proj /= ON) then
                  write(nfout,'(1x,"!* A pseudopotential in use is the Troullier-Martin type")')
                  write(nfout,'(1x,"!* set use_additional_projector = on in Control tag.   LinearRespone Calc. stop")')
                  stop
               end if
            end if
         end do
      end if
      if(printable) write(nfout,'(1x,"!* no problem with transition moment option setting",/)')

    end subroutine check_ts_option
  end subroutine check_PP

  subroutine find_ppc_data_number(nfout)
    !
    !   find data number of core-repair term
    !   nppc_data(it) : data number for it-th type pseudopotential
    !
    integer,intent(in) :: nfout
    integer            :: it, nfpp, nptrans
    nfpp=0
    if(mype == 0) then
       do it=1, ntyp
          nfpp=nfpp+1

          call find_dipole_section(nfpot(nfpp),nfout,it,nptrans)
          if(printable) &
               & write(nfout,'(1x,"! PP transition moment correction data : it = ",i3, &
               & 2x,"number of data read from PP file = ",i3)') it, nptrans
          nppc_data(it)=nptrans
       end do
    end if
    call mpi_bcast(nppc_data,ntyp,mpi_integer,0,MPI_CommGroup,ierr) ! MPI
  end subroutine find_ppc_data_number

  subroutine find_ppc_data_number_it(nfout,it)
    integer,intent(in) :: nfout,it
    integer, save :: nfpp=0
    integer :: nptrans
    if(mype == 0) then
       nfpp=nfpp+1
       call find_dipole_section(nfpot(nfpp),nfout,it,nptrans)
       if(printable) &
            & write(nfout,'(1x,"! PP transition moment correction data : it = ",i3, &
            & 2x,"number of data read from PP file = ",i3)') it, nptrans
       nppc_data(it)=nptrans
    end if
    if(it==ntyp) then
       call mpi_bcast(nppc_data,ntyp,mpi_integer,0,MPI_CommGroup,ierr) ! MPI
       nfpp=0
    endif
  end subroutine find_ppc_data_number_it

  subroutine find_dipole_section(nfp,nfout,it,nptrans)
    !
    !   find dipole section in gncpp2 potential file
    !
    integer,intent(in)   :: nfp, nfout,it
    integer, intent(out) :: nptrans
    integer              :: idipole
    character(len=10)    :: line
    if(mype /= 0) return
10  read(nfp,'(a10)',end=20) line
    idipole=index(line,'DIPOLE')

    if ( idipole/=0 ) then
       read(nfp,*) nptrans
       if ( nptrans==0 ) then
          if (printable) &
               & write(nfout,'(1x,"!* Warning : number of dipole data is zero for atom type ",i3 )') it
!!!          stop
       end if
       return
    end if
    goto 10
! ------
20  if (printable) then
       write(*,'(1x,"!* Warning : dipole section is not found in pseudopotential file")')
       nptrans = 0
!!    stop
    endif

! --
  end subroutine find_dipole_section


  subroutine read_ptrans_data_ek
    implicit none
    !
    !   read core-repair term
    !
    integer :: it,nfpp
    nfpp=0
    if(mype == 0)  then
       do it =1, ntyp
          nfpp=nfpp+1
          if(printable) &
               & write(nfout,'(1x,"! ppc data for it = ",i4," are read: nppc = ",i4)') it,nppc_data(it)
          call read_ptrans_data_ek_core(nfpot(nfpp),nfout,it)
       end do
    end if

    if (printable) write(nfout,'(1x,"!* ntyp = ",i3)') ntyp
    call mpi_bcast(dipole_dxyz_us,ntyp*mnppc*3,mpi_double_precision,0,MPI_CommGroup,ierr) ! MPI
    call mpi_bcast(dipole_tau1,ntyp*mnppc,mpi_integer,0,MPI_CommGroup,ierr)      ! MPI
    call mpi_bcast(dipole_tau2,ntyp*mnppc,mpi_integer,0,MPI_CommGroup,ierr)      ! MPI
    call mpi_bcast(phase_ylm1,ntyp*mnppc,mpi_integer,0,MPI_CommGroup,ierr)       ! MPI
    call mpi_bcast(phase_ylm2,ntyp*mnppc,mpi_integer,0,MPI_CommGroup,ierr)       ! MPI

 end subroutine read_ptrans_data_ek

 subroutine read_ptrans_data_ek_it(it)
    integer, intent(in) :: it
    integer, save :: nfpp=0
    if(mype == 0)  then
       nfpp=nfpp+1
       if(printable) &
            & write(nfout,'(1x,"! ppc data for it = ",i4," are read: nppc = ",i4)') it,nppc_data(it)
       call read_ptrans_data_ek_core(nfpot(nfpp),nfout,it)
    end if

    if (printable) write(nfout,'(1x,"!* ntyp = ",i3)') ntyp
    if(it==ntyp) then
       call mpi_bcast(dipole_dxyz_us,ntyp*mnppc*3,mpi_double_precision,0,MPI_CommGroup,ierr) ! MPI
       call mpi_bcast(dipole_tau1,ntyp*mnppc,mpi_integer,0,MPI_CommGroup,ierr)      ! MPI
       call mpi_bcast(dipole_tau2,ntyp*mnppc,mpi_integer,0,MPI_CommGroup,ierr)      ! MPI
       call mpi_bcast(phase_ylm1,ntyp*mnppc,mpi_integer,0,MPI_CommGroup,ierr)       ! MPI
       call mpi_bcast(phase_ylm2,ntyp*mnppc,mpi_integer,0,MPI_CommGroup,ierr)       ! MPI
       nfpp=0
    endif
 end subroutine read_ptrans_data_ek_it

 subroutine read_ptrans_data_ek_core(nfp,nfout,it)
   implicit none
!
!      read core-repair term pij(I)
!      cf. H. Kageshima and K. Shiraishi,Phys. Rev. B vol.56, 14985 (1997)
!      pij(I)
!      i, j :: atomic orbital index (CIAO)
!      i ->  ltmltm(n1,l1,t1,m1)  j-> ltmltm(n2,l2,t1,m1)  ltm: compound index
!      n :: principal quantum number    l :: azimuthal quantum number
!      m :: magnetic quantum number     t :: energy reference index
!      dipole_dxyz_us(it,ltmltm,ixyz)     :: core repair term for it-th type pseudopotential in ixyz direction
!                                         :: ixyz=1 -> x; ixyz=2 -> y;  ixyz= 3 -> z
!      phase_ylm  :: Ylm index for phase
!      dipole_tau :: reference index
!
!
   integer, intent(in) :: nfp,nfout,it
   integer :: n1,l1,t1,m1,n2,l2,t2,m2
   integer :: ltmltm,lmax
   do ltmltm = 1, nppc_data(it)
!         n1 = n1_dipole_lm_us(ltmltm) ; n2 = n2_dipole_lm_us(ltmltm)
!         l1 = l1_dipole_lm_us(ltmltm) ; l2 = l2_dipole_lm_us(ltmltm)
!         t1 = t1_dipole_lm_us(ltmltm) ; t2 = t2_dipole_lm_us(ltmltm)
!         m1 = m1_dipole_lm_us(ltmltm) ; m2 = m2_dipole_lm_us(ltmltm)
      read (nfp,53) n1,l1,t1,m1,n2,l2,t2,m2, &
      dipole_dxyz_us(it,ltmltm,1),dipole_dxyz_us(it,ltmltm,2),dipole_dxyz_us(it,ltmltm,3), &
      phase_ylm1(it,ltmltm), phase_ylm2(it,ltmltm)

      dipole_tau1(it,ltmltm)=t1;  dipole_tau2(it,ltmltm)=t2

      if ( printable ) then
         write(nfout,53) n1,l1,dipole_tau1(it,ltmltm),m1,n2,l2,dipole_tau2(it,ltmltm),m2, &
       & dipole_dxyz_us(it,ltmltm,1),dipole_dxyz_us(it,ltmltm,2),dipole_dxyz_us(it,ltmltm,3), &
       & phase_ylm1(it,ltmltm),phase_ylm2(it,ltmltm)
      end if
   end do
53 format(1x,8i3,3e18.10,2i3)
 end subroutine read_ptrans_data_ek_core

 subroutine prepare_for_TM_PP_ek
   implicit none
    integer :: lmax,ilocal,il1,il2,lmt
    integer :: it,lmt1,lmt2

    if(mype == 0) then
       ilocal=0
       do it = 1, ntyp
          lmax=lpsmax(it)
          call check_local_orbital
          if(printable) then
             if(ilocal==0) then
                write(nfout,'(1x,"!* all non-local projectors are given for it =",i3)') it
                write(nfout,'(1x,"!* skip correction for Troullier-Martin pseudopotential for it =",i3)') it
             else
                write(nfout,'(1x,"!* ptrans correction for Troullier-Matrtin pseudopotential for it =",i3)') it
                write(nfout,'(1x,"!* non-local projector for iloc =",i3," is not given")') ilocal
                write(nfout,'(1x,"!* correction for iloc =",i3," is made")') ilocal
             end if
          end if
          ilocal_l(it)=ilocal
       end do
    end if
    call mpi_bcast(ilocal_l,ntyp,mpi_integer,0,MPI_CommGroup,ierr)
    contains
     subroutine check_local_orbital
       do il1=1,lmax
          do lmt=1,ilmt(it)
             il2=ltp(lmt,it)
             if(il1==il2) exit
          end do
          if(il1/=il2) then
             ilocal=il1
          end if
       end do
     end subroutine check_local_orbital
 end subroutine prepare_for_TM_PP_ek


 subroutine calc_ptrans_ek
    implicit none
!
!   calculate KS correction term
!
    integer                                :: id_sname = -1
    integer                                :: ispin, it, lmt1, lmt2, il1, im1, il2, im2, ia
    integer                                :: ik, ii, ib, ib1, ilmta, p, p1, index, ifact
    integer                                :: nspher1,nspher2
    real(kind=DP) :: fac, eib, eib1
    real(kind=DP),pointer,dimension(:,:,:) :: wkfsr, wkfsi
! --> T. Yamasaki 2008/02/21
    real(kind=DP),pointer,dimension(:)     :: wkfsr_tmp, wkfsi_tmp
! <-- T. Yamasaki 2008/02/21


    call tstatc0_begin('calc_ptans_ek ',id_sname)
    allocate(wkfsr(neg,nlmta,kv3)); allocate(wkfsi(neg,nlmta,kv3))

    if(npes >= 2) call mpi_barrier(MPI_CommGroup,ierr)

!   make copy of <beta|WF>
! --> T. Yamasaki 2008/02/21
    allocate(wkfsr_tmp(nlmta),wkfsi_tmp(nlmta))
!!$    do ik = 1, kv3, af+1
!!$       do ib = 1, neg
!!$          if(map_ek(ib,ik) == mype) then
!!$               do ilmta=1, nlmta
!!$                  wkfsr(ib,ilmta,ik) = fsr_l(map_z(ib),ilmta,ik)
!!$                  wkfsi(ib,ilmta,ik) = fsi_l(map_z(ib),ilmta,ik)
!!$               end do
!!$               if(map_ek(ib,ik) /= 0) then
!!$                   call mpi_send(wkfsr,neg*nlmta*kv3,mpi_double_precision,0,1,MPI_CommGroup,ierr) ! MPI
!!$                   call mpi_send(wkfsi,neg*nlmta*kv3,mpi_double_precision,0,1,MPI_CommGroup,ierr) ! MPI
!!$               end if
!!$          else if(mype == 0 .and. map_ek(ib,ik) /= 0) then
!!$            call mpi_recv(wkfsr,neg*nlmta*kv3,mpi_double_precision,map_ek(ib,ik),1,MPI_CommGroup,istatus,ierr)!MPI
!!$            call mpi_recv(wkfsi,neg*nlmta*kv3,mpi_double_precision,map_ek(ib,ik),1,MPI_CommGroup,istatus,ierr)!MPI
!!$         end if
!!$         if(npes >= 2)  then
!!$              call mpi_bcast(wkfsr,neg*nlmta*kv3,mpi_double_precision,0,MPI_CommGroup,ierr)!MPI
!!$              call mpi_bcast(wkfsi,neg*nlmta*kv3,mpi_double_precision,0,MPI_CommGroup,ierr)!MPI
!!$         end if
!!$       end do
!!$    end do
    do ik = 1, kv3, af+1
       do ib = 1, neg
          if(map_ek(ib,ik) == mype) then
             if(mype == 0) then
                do ilmta=1, nlmta
                   wkfsr(ib,ilmta,ik) = fsr_l(map_z(ib),ilmta,ik)
                   wkfsi(ib,ilmta,ik) = fsi_l(map_z(ib),ilmta,ik)
                end do
             else
                do ilmta=1, nlmta
                   wkfsr_tmp(ilmta) = fsr_l(map_z(ib),ilmta,ik)
                   wkfsi_tmp(ilmta) = fsi_l(map_z(ib),ilmta,ik)
                end do
                call mpi_send(wkfsr_tmp,nlmta,mpi_double_precision,0,1,MPI_CommGroup,ierr) ! MPI
                call mpi_send(wkfsi_tmp,nlmta,mpi_double_precision,0,1,MPI_CommGroup,ierr) ! MPI
             end if
          else if(mype == 0 .and. map_ek(ib,ik) /= 0) then
             call mpi_recv(wkfsr_tmp,nlmta,mpi_double_precision,map_ek(ib,ik),1,MPI_CommGroup,istatus,ierr)!MPI
             call mpi_recv(wkfsi_tmp,nlmta,mpi_double_precision,map_ek(ib,ik),1,MPI_CommGroup,istatus,ierr)!MPI
             do ilmta=1, nlmta
                wkfsr(ib,ilmta,ik) = wkfsr_tmp(ilmta)
                wkfsi(ib,ilmta,ik) = wkfsi_tmp(ilmta)
             end do
         end if
       end do
    end do
    if(npes >= 2)  then
       call mpi_bcast(wkfsr,neg*nlmta*kv3,mpi_double_precision,0,MPI_CommGroup,ierr)!MPI
       call mpi_bcast(wkfsi,neg*nlmta*kv3,mpi_double_precision,0,MPI_CommGroup,ierr)!MPI
    end if

    deallocate(wkfsr_tmp,wkfsi_tmp)
! <-- T. Yamasaki 2008/02/21

    ptrans=0.0d0

!   calculate sum[<WF1|beta(i)>pij<beta(j)|WF2>
    do ii = 1,3
       do ispin = 1, nspin, af+1
          do ik = ispin, kv3-nspin+ispin, nspin
             do ib = 1, neg
                do ib1 = 1, neg

!                   if(nrd_efermi==1) then
!                      eib=e2_mpi(ib,ik)
!                      eib1=e2_mpi(ib1,ik)
!                      if(eib.gt.efermi.and.eib1.le.efermi) call calc_ptrans_ek_core(ib,ib1)                   else
                      if(ib/=ib1) call calc_ptrans_ek_core(ib,ib1)
!                   end if
                end do
             end do
          end do
       end do
    end do

!!$    if(nk_in_the_process + kv3-1 >= kv3_ek) stop ' m_Epsilon_ek (1)'

    ptrans = ptrans*(af+1)

    deallocate(wkfsr)
    deallocate(wkfsi)
    call tstatc0_end(id_sname)

    contains
     subroutine find_ptrans_index_ek(it,lmt1,lmt2,nptrans,index,ifact,nspher1,nspher2)
       implicit none
!
!      find core repair term pij for <WF1|beta(i)>pij<beta(j)|WF2>
!      i -> (it,lmt1)   j -> (it,lmt2)
!
       integer, intent(in)  :: it, lmt1, lmt2, nptrans
       integer, intent(out) :: index, ifact
       integer              :: lmtt1, il1, im1, tau1, nspher1
       integer              :: lmtt2, il2, im2, tau2, nspher2
       integer              :: nspher10, nspher20
       integer              :: tau10, tau20
       integer              :: iptrans
       index=0
       ifact=0
       call m_PP_tell_lmtt_l_m_tau(lmt1,it,lmtt1,il1,im1,tau1,nspher1)
       call m_PP_tell_lmtt_l_m_tau(lmt2,it,lmtt2,il2,im2,tau2,nspher2)
       if(nspher1>nspher2) then
          nspher20=nspher1
          nspher10=nspher2
          tau10=tau2
          tau20=tau1
          ifact=-1
       else
          nspher10=nspher1
          nspher20=nspher2
          tau10=tau1
          tau20=tau2
          ifact=1
       end if

!      find core-repair term with dipole_tau = tau, phase_ylm = nspher
       do iptrans=1,nptrans
          if(phase_ylm1(it,iptrans)==nspher10.and.phase_ylm2(it,iptrans)==nspher20) then
             if(dipole_tau1(it,iptrans)==tau10.and.dipole_tau2(it,iptrans)==tau20) then
                index=iptrans
                exit
             end if
          end if
       end do
     end subroutine find_ptrans_index_ek

     subroutine calc_ptrans_ek_core(ib,ib1)
!
!      calculate <WF1|beta(i)>pij<beta(j)|WF2>
!
       integer,intent(in) :: ib,ib1
       integer            :: index,iv,ic,ifind

!       call find_ind_vb_and_cb2(ib1,ib,iv,ic,nk_in_the_process+ik-1,ifind)
!       if(ifind==0.and.printable) then
!          write(nfout,'(1x,"!!* valence or conduction band index is not found   UVSOR-Epsilon STOP at calc_ptrans")')
!       end if

! --
       iv = ib1
       ic = ib
! ---

       do it=1,ntyp
          do lmt1 = 1, ilmt(it)
             do lmt2 = 1, ilmt(it)
                il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
                call find_ptrans_index_ek(it,lmt1,lmt2,nppc_data(it),index,ifact,nspher1,nspher2)
                if(index==0) cycle
                do ia = 1, natm
                   if(ityp(ia) /= it) cycle
                   p = lmta(lmt1,ia)
                   p1 = lmta(lmt2,ia)
                   fac=real(iwei(ia),kind=DP)
                   fac = real(ifact,kind=DP)*fac
! --- debug
!                   write(*,*) it, lmt1, lmt2, ia, fac, dipole_dxyz_us(it,index,ii)
! --

                   ptrans(ik,ic,iv,ii,1) = ptrans(ik,ic,iv,ii,1) &
                & + fac*(dipole_dxyz_us(it,index,ii)) &
                & *(wkfsr(ib,p,ik)*wkfsi(ib1,p1,ik) - wkfsi(ib,p,ik)*wkfsr(ib1,p1,ik))
                   ptrans(ik,ic,iv,ii,2) = ptrans(ik,ic,iv,ii,2) &
                & -1.0d0* fac*(dipole_dxyz_us(it,index,ii)) &
                & *(wkfsr(ib,p,ik)*wkfsr(ib1,p1,ik) + wkfsi(ib,p,ik)*wkfsi(ib1,p1,ik))
!
!                   write(*,*) 'G ', it, lmt1, lmt2, ia, ptrans(ik,ic,iv,ii,1), ptrans(ik,ic,iv,ii,2)
!
                end do
             end do
          end do
       end do

     end subroutine calc_ptrans_ek_core
 end subroutine calc_ptrans_ek

 subroutine calc_ptrans_TM_PP_ek
!
!   calculate KS correction term for Troullier-Martin PP
!
    implicit none
    integer                                :: id_sname = -1
    integer                                :: ispin, it, lmt, lmt1, lmt2, il, il1, il2, ia
    integer                                :: ik, ii, ib, ib1, ilmta, p, p1
    integer                                :: ilmta_add, nspher1, nspher2
    integer                                :: icount, im1
    real(kind=DP)                          :: fac, eib, eib1
    real(kind=DP),pointer,dimension(:,:,:) :: wkfsr, wkfsi, wkfsr_add, wkfsi_add

    call tstatc0_begin('calc_ptrans_TM_PP_ek ',id_sname)
    allocate(wkfsr(neg,nlmta,kv3)); wkfsr=0.0d0
    allocate(wkfsi(neg,nlmta,kv3)); wkfsi=0.0d0
    allocate(wkfsr_add(neg,nlmta_add,kv3)); wkfsr_add=0.0d0
    allocate(wkfsi_add(neg,nlmta_add,kv3)); wkfsi_add=0.0d0

    if(sw_use_add_proj==ON) then
       call m_ES_add_betar_dot_WFs(nfout)
       if(printable) &
       & write(nfout,'(1x,"!* fsr_add_l and fsi_add_l have been calculated")')
    else
       if(printable) then
          write(nfout,'(1x,"!* fsr_add_l and fsi_add_l cannot be callculated")')
          write(nfout,'(1x,"!* skip correction for additional projector ")')
       end if
       return
    end if

    if(npes >= 2) call mpi_barrier(MPI_CommGroup,ierr)

! make copy of fsr_l and fsi_l
    do ik = 1, kv3, af+1
       do ib = 1, neg
          if(map_ek(ib,ik) == mype) then
               do ilmta=1, nlmta
                  wkfsr(ib,ilmta,ik) = fsr_l(map_z(ib),ilmta,ik)
                  wkfsi(ib,ilmta,ik) = fsi_l(map_z(ib),ilmta,ik)
               end do
               if(map_ek(ib,ik) /= 0) then
                   call mpi_send(wkfsr,neg*nlmta*kv3,mpi_double_precision,0,1,MPI_CommGroup,ierr) ! MPI
                   call mpi_send(wkfsi,neg*nlmta*kv3,mpi_double_precision,0,1,MPI_CommGroup,ierr) ! MPI
               end if
          else if(mype == 0 .and. map_ek(ib,ik) /= 0) then
            call mpi_recv(wkfsr,neg*nlmta*kv3,mpi_double_precision,map_ek(ib,ik),1,MPI_CommGroup,istatus,ierr)!MPI
            call mpi_recv(wkfsi,neg*nlmta*kv3,mpi_double_precision,map_ek(ib,ik),1,MPI_CommGroup,istatus,ierr)!MPI
         end if
         if(npes >= 2)  then
              call mpi_bcast(wkfsr,neg*nlmta*kv3,mpi_double_precision,0,MPI_CommGroup,ierr)!MPI
              call mpi_bcast(wkfsi,neg*nlmta*kv3,mpi_double_precision,0,MPI_CommGroup,ierr)!MPI
         end if
       end do
    end do

! make copy of fsr_add_l and fsi_add_l
    do ik = 1, kv3, af+1
       do ib = 1, neg
          if(map_ek(ib,ik) == mype) then
               do ilmta_add = 1, nlmta_add
                  wkfsr_add(ib,ilmta_add,ik) = fsr_add_l(map_z(ib),ilmta_add,ik)
                  wkfsi_add(ib,ilmta_add,ik) = fsi_add_l(map_z(ib),ilmta_add,ik)
               end do
               if(map_ek(ib,ik) /= 0) then
                  call mpi_send(wkfsr_add,neg*nlmta_add*kv3,mpi_double_precision,0,1,MPI_CommGroup,ierr) ! MPI
                  call mpi_send(wkfsi_add,neg*nlmta_add*kv3,mpi_double_precision,0,1,MPI_CommGroup,ierr) ! MPI
               end if
          else if(mype == 0 .and. map_ek(ib,ik) /= 0) then
            call mpi_recv(wkfsr_add,neg*nlmta_add*kv3,mpi_double_precision,map_ek(ib,ik),1,MPI_CommGroup,istatus,ierr)!MPI
            call mpi_recv(wkfsi_add,neg*nlmta_add*kv3,mpi_double_precision,map_ek(ib,ik),1,MPI_CommGroup,istatus,ierr)!MPI
         end if
         if(npes >= 2)  then
              call mpi_bcast(wkfsr_add,neg*nlmta_add*kv3,mpi_double_precision,0,MPI_CommGroup,ierr)!MPI
              call mpi_bcast(wkfsi_add,neg*nlmta_add*kv3,mpi_double_precision,0,MPI_CommGroup,ierr)!MPI
         end if
       end do
    end do

! skip if PP has no local orbital
    do it = 1, ntyp
       !!if(printable) then
       !!  write(nfout,*) 'local orbital: it=',it
       !!  write(nfout,*) 'local orbital: l=',ilocal_l(it)
       !!end if
       if(ilocal_l(it)==0) then
          if(printable) then
             write(nfout,'(1x,"!* ilocal =0 for it = ",i3)') it
             write(nfout,'(1x,"!* additional projector correction is skipped")')
          end if
          cycle
       end if

! check <additional projector|WF> for l=ilocal
       do lmt=1,ilmt_add(it)
          il=ltp_add(lmt,it)
          if(il==ilocal_l(it)) then
             exit
          else
             if(lmt== ilmt_add(it)) then
                if(printable) then
                   write(nfout,'(1x,"!* additional projector for for iloc is not given for it = ",i3)') it
                   write(nfout,'(1x,"!* calculate additional projector for iloc")')
                end if
                stop
             end if
          end if
       end do

! calculate correction for local orbital
! add sum[<WF1|additional projector>p<beta|WF2>] term
       !!if(printable) write(nfout,*) 'nrd_efermi=',nrd_efermi
       icount=0
       do ii = 1, 3
          do ispin = 1, nspin, af+1
             do ik = ispin, kv3-nspin+ispin, nspin
                do ib = 1, neg
                   do ib1 = 1, neg

!                      if(nrd_efermi==1) then
!                         eib=e2_mpi(ib,ik)
!                         eib1=e2_mpi(ib1,ik)
                         !!if(printable) write(nfout,*) 'calc_ptrans_TM_PP_ek_core_1, eib,eib1=',eib,eib1
!                         if(eib.gt.efermi.and.eib1.le.efermi) call calc_ptrans_TM_PP_ek_core_1
!                      else
                         !!if(printable) write(nfout,*) 'calc_ptrans_TM_PP_ek_core_1, ib,ib1=',ib,ib1
                         if(ib/=ib1) call calc_ptrans_TM_PP_ek_core_1
!                      end if
                   end do
                end do
             end do
          end do
       end do
       if(printable) write(nfout,'(1x,"!* number of additional ppc1 terms = ", i8)') icount
! add sum sum[<WF1|beta>p<additional projector|WF2>] term
       icount=0
       do ii = 1, 3
          do ispin = 1, nspin, af+1
             do ik = ispin, kv3-nspin+ispin, nspin
                do ib = 1, neg
                   do ib1 = 1, neg
!                      if(nrd_efermi==1) then
!                         eib=e2_mpi(ib,ik)
!                         eib1=e2_mpi(ib1,ik)
!                         if(eib.gt.efermi.and.eib1.le.efermi) call calc_ptrans_TM_PP_ek_core_2
!                      else
                         if(ib/=ib1) call calc_ptrans_TM_PP_ek_core_2
!                      end if
                   end do
                end do
             end do
          end do
       end do
       if(printable) write(nfout,'(1x,"!* number of additional ppc2 terms = ", i8)') icount
    end do

    deallocate(wkfsr)
    deallocate(wkfsi)
    deallocate(wkfsr_add)
    deallocate(wkfsi_add)

    call tstatc0_end(id_sname)

    contains
     subroutine calc_ptrans_TM_PP_ek_core_1
       integer :: iv, ic, ifind
       integer :: index, ifact

!       call find_ind_vb_and_cb2(ib1,ib,iv,ic,nk_in_the_process+ik-1,ifind)
!       if(ifind==0.and.printable) then
!          write(nfout,'(1x,"!!* valence or conduction band index is not found   UVSOR-Epsilon STOP at calc_ptrans")')
!       end if

! -----------------
       iv = ib1
       ic = ib
! ----------------

       do lmt1 = 1, ilmt_add(it)
          il1 = ltp_add(lmt1,it)
          if(il1/=ilocal_l(it)) cycle
          do lmt2 = 1, ilmt(it)
! lmt1 -> index for additional projector
! lmt2 -> index for beta
             !!if(printable) &
             !!& write(nfout,*) 'find_ptrans_index_add_beta: it,lmt1,il1,lmt2,nppc_data(it)=',it,lmt1,il1,lmt2,nppc_data(it)
             call find_ptrans_index_add_beta(it,lmt1,lmt2,nppc_data(it),index,ifact)
             !!if(printable) write(nfout,*) 'find_ptrans_index_add_beta: index,ifact=',index,ifact
             if(index==0) cycle
             do ia = 1, natm
                if(ityp(ia) /= it) cycle
                fac = real(ifact,kind=DP)*real(iwei(ia),kind=DP)
                p =  lmta_add(lmt1,ia)
                p1 = lmta(lmt2,ia)
                ptrans(ik,ic,iv,ii,1) = ptrans(ik,ic,iv,ii,1) &
             & + fac*dipole_dxyz_us(it,index,ii) &
             & *(wkfsr_add(ib,p,ik)*wkfsi(ib1,p1,ik) - wkfsi_add(ib,p,ik)*wkfsr(ib1,p1,ik))
                ptrans(ik,ic,iv,ii,2) = ptrans(ik,ic,iv,ii,2) &
             & -1.0d0*fac*dipole_dxyz_us(it,index,ii) &
             & *(wkfsr_add(ib,p,ik)*wkfsr(ib1,p1,ik) + wkfsi_add(ib,p,ik)*wkfsi(ib1,p1,ik))
                icount=icount+1
             end do
          end do
       end do
     end subroutine calc_ptrans_TM_PP_ek_core_1

     subroutine calc_ptrans_TM_PP_ek_core_2
       integer :: iv, ic, ifind
       integer :: index, ifact

!       call find_ind_vb_and_cb2(ib1,ib,iv,ic,nk_in_the_process+ik-1,ifind)
!       if(ifind==0.and.printable) then
!          write(nfout,'(1x,"!!* valence or conduction band index is not found   UVSOR-Epsilon STOP at calc_ptrans")')
!       end if

! ------------------
       iv = ib1
       ic = ib
! -----------------
       do lmt1 = 1, ilmt(it)
          do lmt2 = 1, ilmt_add(it)
! lmt1 --> index for beta
! lmt2 --> index for additional projector
             il2 = ltp_add(lmt2,it)
             if(il2/=ilocal_l(it)) cycle
             call find_ptrans_index_add_beta(it,lmt2,lmt1,nppc_data(it),index,ifact)
             if(index==0) cycle
             ifact=-1*ifact
             do ia = 1, natm
                if(ityp(ia) /= it) cycle
                fac = real(ifact,kind=DP)*real(iwei(ia),kind=DP)
                p = lmta(lmt1,ia)
                p1 = lmta_add(lmt2,ia)
! p  --> index for beta
! p1 --> index for addtional projector
                ptrans(ik,ic,iv,ii,1) = ptrans(ik,ic,iv,ii,1) &
             & + fac*dipole_dxyz_us(it,index,ii) &
             & *(wkfsr(ib,p,ik)*wkfsi_add(ib1,p1,ik) - wkfsi(ib,p,ik)*wkfsr_add(ib1,p1,ik))
                ptrans(ik,ic,iv,ii,2) = ptrans(ik,ic,iv,ii,2) &
             & -1.0d0*fac*dipole_dxyz_us(it,index,ii) &
             & *(wkfsr(ib,p,ik)*wkfsr_add(ib1,p1,ik) + wkfsi(ib,p,ik)*wkfsi_add(ib1,p1,ik))
                icount=icount+1
             end do
          end do
       end do
     end subroutine calc_ptrans_TM_PP_ek_core_2

     subroutine find_ptrans_index_add_beta(it,lmt1,lmt2,nptrans,index,ifact)
       implicit none
!
!      find core-repair term pij for <WF|addotional projector(i)>pij<beta(j)|WF2>
!
!      lmt1      : lmt index for additinal projector
!      lmt2      : lmt index for beta
!      nspher1   : spherical harmonics index for lmt1
!      nspher2   : spherical harmonics index for lmt2
!      ifact = 1 : nspher1 < nspher2
!      ifact = -1: nspher1 > nspher2
!      index     : dipole_dxyz_us index
!            = 0 : for nspher1=nspher2 case
!
       integer, intent(in)  :: it,lmt1,lmt2, nptrans
       integer, intent(out) :: index,ifact
       integer              :: lmtt1,il1,im1,tau1,nspher1
       integer              :: lmtt2,il2,im2,tau2,nspher2
       integer              :: nspher10, nspher20
       integer              :: tau10,tau20
       integer              :: iptrans

       index=0
       ifact=0

       !!if(printable) write(nfout,*) 'm_PP_tell_lmtt_l_m_tau_add:lmt1=',lmt1
       call m_PP_tell_lmtt_l_m_tau_add(lmt1,it,lmtt1,il1,im1,tau1,nspher1)
       !!if(printable) then
       !!   write(nfout,*) 'm_PP_tell_lmtt_l_m_tau_add:tau1,nspher1=',tau1,nspher1
       !!   write(nfout,*) 'm_PP_tell_lmtt_l_m_tau:lmt2=',lmt2
       !!end if
       call m_PP_tell_lmtt_l_m_tau(lmt2,it,lmtt2,il2,im2,tau2,nspher2)
       !!if(printable) write(nfout,*) 'm_PP_tell_lmtt_l_m_tau:tau2,nspher2=',tau2,nspher2

       if(nspher1>nspher2) then
          nspher20=nspher1
          nspher10=nspher2
          tau20=tau1
          tau10=tau2
          ifact=-1
       else
          nspher10=nspher1
          nspher20=nspher2
          tau10=tau1
          tau20=tau2
         ifact=1
       end if
       !!if(printable) write(nfout,*) 'loop iptans'
!      find core repair term with dipole_tau=tau, phase_ylm=nspher
       do iptrans=1,nptrans
          if(phase_ylm1(it,iptrans)==nspher10.and.phase_ylm2(it,iptrans)==nspher20) then
             if(dipole_tau1(it,iptrans)==tau10.and.dipole_tau2(it,iptrans)==tau20) then
                index=iptrans
                exit
              end if
          end if
       end do
       !!if(printable) write(nfout,*) 'end loop iptans'
     end subroutine find_ptrans_index_add_beta

   end subroutine calc_ptrans_TM_PP_ek

 end module m_LinearResponse_Density
