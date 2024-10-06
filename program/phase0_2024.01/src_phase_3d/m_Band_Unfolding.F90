module m_Band_Unfolding
  use m_Const_Parameters,   only : PAI2, CARTS, DP, OFF, MAPPED, NOTMAPPED, SKIP, &
       &                           GAMMA, BUCS, CMPLDP, ON, Neglected
  use m_Control_Parameters,  only : ekmode, neg, noncol, num_extra_bands, &
       &                            ndim_spinor, kimg, flag_mpi_g_dot_r, nspin, af, &
       &                            sw_orb_popu, orb_popu_method, SpinOrbit_Mode, &
       &                            wf_orb_proj_print_format, num_proj_elems, &
       &                            proj_group, proj_attribute, &
       &                            sw_calc_wf_orb_projection, band_unfolding_active, &
       &                            tolerance_Gvec_matching

  use m_PseudoPotential,    only : modnrm, nlmta, ilmt, lmta, q, ilmt_phi, ltp_phi, &
       &                           mtp_phi, taup_phi, lmta_phi, lmtt_phi, &
       &                           ltp, mtp, taup, nlmtt_phi, nlmta_phi, &
       &                           q_phirt_pw, qorb, m_PP_tell_iorb_lmtt, &
       &                           m_PP_tell_iorb_ia_l_m_tau
  use m_Crystal_Structure,  only : altv_refcell, rltv, univol
  use m_Kpoints,  only : kv3, vkxyz, k_symmetry, kv3_ek, vkxyz_refcell
  use m_Files,        only : nfout, nfband_spectr_wght, &
       &                     m_Files_open_nfband_spwt, &
       &                     m_Files_close_nfband_spwt, nfwfk_orb_proj, &
       &                     m_Files_open_nfwfk_orb_proj, &
       &                     m_Files_close_nfwfk_orb_proj

  use m_PlaneWaveBasisSet,    only : kg, ngabc, iba, nbase, nbmx, GVec_on_refcell

  use m_Parallelization,  only : mype, npes, myrank_k, map_k, ista_e, iend_e, map_z, &
       &                         MPI_CommGroup, ista_nbmx, mp_nbmx, np_nbmx, ng_nbmx, &
       &                         mpi_nbmx_world, ierr, ista_k, iend_k, np_e, istep_e, &
       &                         nrank_e, myrank_e, map_e
  use m_Parallelization,  only : ista_g1k, iend_g1k, np_g1k, mpi_ke_world, &
       &                         nrank_g, np_fs, ista_fs, iend_fs, ista_atm, iend_atm

  use m_IterationNumbers,  only : nk_in_the_process
  use m_Ionic_System,   only : pos, natm, ityp, altv, iproj_group
  use m_Electronic_Structure,  only : zaj_l, neordr, fsr_l, fsi_l
  use m_ES_nonlocal,     only : zfcos, zfsin, zfcos_mpi, zfsin_mpi, &
       &                        m_ES_alloc_zfsincos, m_ES_dealloc_zfsincos, &
       &                        alloc_zfsincos_mpi, dealloc_zfsincos_mpi, &
       &                        m_ES_phir_dot_WFs_3D, &
       &                        m_ES_betar_dot_Psi_4_each_k_3D
  use m_NonLocal_Potential,  only : snl, phig, norm_phig

  use m_Kpoints,  only :  sw_force_kpt_inside_bz,  GvecTrans_kpt

  use m_SpinOrbit_Potential,  only :  MatU_ylm_RC_L0,  MatU_ylm_RC_L1,  MatU_ylm_RC_L2, &
       &                              MatU_ylm_RC_L3
  use mpi

  implicit none
!  include 'mpif.h'

  real(kind=DP), allocatable :: fsr_filtered(:,:,:), fsi_filtered(:,:,:)
  real(kind=DP), allocatable :: compr_filtered(:,:,:,:), compi_filtered(:,:,:,:)
contains

  subroutine m_BU_set_GVec_Flag_refcell
    integer :: i, j, k, ik
    real(kind=DP) :: vec1(3), cx, cy, cz, delta

    if ( allocated( Gvec_on_refcell ) ) deallocate( Gvec_on_refcell )
    if ( sw_force_kpt_inside_bz == ON ) then
       allocate( GVec_on_refcell(kg,ista_k:iend_k) );
    else
       allocate( GVec_on_refcell(kg,1) );
    endif
    GVec_on_refcell = 0

    delta = tolerance_Gvec_matching

    if ( sw_force_kpt_inside_bz == ON ) then
       Do ik=ista_k, iend_k
          Do i=1, kg
             vec1(1) = ( ngabc(i,1) +GVecTrans_kpt(ik,1) )*rltv(1,1) &
                  &   +( ngabc(i,2) +GVecTrans_kpt(ik,2) )*rltv(1,2) &
                  &   +( ngabc(i,3) +GVecTrans_kpt(ik,3) )*rltv(1,3)
             vec1(2) = ( ngabc(i,1) +GVecTrans_kpt(ik,1) )*rltv(2,1) &
                  &   +( ngabc(i,2) +GVecTrans_kpt(ik,2) )*rltv(2,2) &
                  &   +( ngabc(i,3) +GVecTrans_kpt(ik,3) )*rltv(2,3)
             vec1(3) = ( ngabc(i,1) +GVecTrans_kpt(ik,1) )*rltv(3,1) &
                  &   +( ngabc(i,2) +GVecTrans_kpt(ik,2) )*rltv(3,2) &
                  &   +( ngabc(i,3) +GVecTrans_kpt(ik,3) )*rltv(3,3)
             cx = 0.0d0; cy = 0.0d0; cz = 0.0d0
             Do k=1, 3
                cx = cx +altv_refcell(k,1) *vec1(k)
                cy = cy +altv_refcell(k,2) *vec1(k)
                cz = cz +altv_refcell(k,3) *vec1(k)
             End Do

             cx = cx /PAI2;     cy = cy /PAI2;      cz = cz /PAI2
             cx = abs( cx -nint(cx) );
             cy = abs( cy -nint(cy) );
             cz = abs( cz -nint(cz) );

             if ( cx < delta .and. cy < delta .and. cz < delta ) then
                GVec_on_refcell(i,ik) = 1
             endif
          ENd Do
       End Do
    else
       Do i=1, kg
          vec1(1) = ngabc(i,1) *rltv(1,1) +ngabc(i,2) *rltv(1,2) +ngabc(i,3) *rltv(1,3)
          vec1(2) = ngabc(i,1) *rltv(2,1) +ngabc(i,2) *rltv(2,2) +ngabc(i,3) *rltv(2,3)
          vec1(3) = ngabc(i,1) *rltv(3,1) +ngabc(i,2) *rltv(3,2) +ngabc(i,3) *rltv(3,3)

          cx = 0.0d0; cy = 0.0d0; cz = 0.0d0
          Do k=1, 3
             cx = cx +altv_refcell(k,1) *vec1(k)
             cy = cy +altv_refcell(k,2) *vec1(k)
             cz = cz +altv_refcell(k,3) *vec1(k)
          End Do

          cx = cx /PAI2;     cy = cy /PAI2;      cz = cz /PAI2
          cx = abs( cx -nint(cx) );
          cy = abs( cy -nint(cy) );
          cz = abs( cz -nint(cz) );

          if ( cx < delta .and. cy < delta .and. cz < delta ) then
             GVec_on_refcell(i,1) = 1
          endif
       ENd Do
    endif

  end subroutine m_BU_set_GVec_Flag_refcell


  subroutine m_BU_calc_spectral_weight_3D
    real(kind=DP), allocatable :: SpectralWeight(:,:)
    integer :: ik, ie, neg_t
    integer :: lun, j, k

    real(kind=DP) :: c1
    real(kind=DP), allocatable :: work(:,:), work2(:)
    integer, allocatable :: iwork(:,:), neordr_t(:,:)

    band_unfolding_active = .true.

    if ( modnrm /= SKIP ) then
       allocate( fsr_filtered(np_e,np_fs,ista_k:iend_k) )
       if( .not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
          allocate( fsi_filtered(np_e,np_fs,ista_k:iend_k));
       else
          allocate( fsi_filtered(1,1,1) )
       endif
       fsr_filtered = 0.0d0;     fsi_filtered = 0.d0
       call m_BU_betar_dot_WFs_3D(nfout)
    end if

    allocate( SpectralWeight(neg,kv3) ); SpectralWeight = 0.0d0

    call softpart
    if ( modnrm /= SKIP ) call hardpart

    if ( npes > 1 ) then
       allocate( work(neg,kv3) ); work = 0.0d0
       call mpi_allreduce( SpectralWeight, work, kv3 *neg, &
            &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
       SpectralWeight = work
       deallocate( work )
    endif

    neg_t = neg -num_extra_bands

    allocate( work2(neg) )
#if 1
    allocate( neordr_t(neg,kv3) ); neordr_t = 0
    neordr_t(:,ista_k:iend_k) = neordr(:,ista_k:iend_k)
    if ( npes > 1 ) then
       allocate( iwork(neg,kv3) )
       call mpi_allreduce( neordr_t, iwork, neg*kv3, mpi_integer, mpi_sum, &
            &              MPI_CommGroup, ierr )
       neordr_t = iwork /nrank_e /nrank_g
       deallocate( iwork )
    endif
#endif

    if ( mype == 0 ) then
       if ( ekmode == ON ) then
          if ( nk_in_the_process == 1 ) then
             call m_Files_open_nfband_spwt(2)
             write(nfband_spectr_wght,'(A,I8)') 'num_kpoints = ', kv3_ek /ndim_spinor
             write(nfband_spectr_wght,'(A,I8)') 'num_bands   = ', neg_t
             write(nfband_spectr_wght,'(A,I8)') 'nspin       = ', nspin /ndim_spinor
             write(nfband_spectr_wght,*)
          else
             call m_Files_open_nfband_spwt(3)
          endif
          Do ik=1, kv3, ndim_spinor
             if ( ik+nk_in_the_process -1 > kv3_ek ) cycle
             write( nfband_spectr_wght,'(A,I6,A,3F16.8,A)') &
                  &             "ik = ", ik +nk_in_the_process -1, " ( ", &
                  &             vkxyz_refcell(ik+nk_in_the_process-1,1:3,BUCS), " )"
             if ( noncol ) then
                Do ie=1, neg
#if 1
                   work2(ie) = SpectralWeight(neordr_t(ie,ik),ik) &
                        &     +SpectralWeight(neordr_t(ie,ik+1),ik+1)
#else
                   work2(ie) = SpectralWeight(ie,ik) &
                        &     +SpectralWeight(ie,ik+1)
#endif
                End Do
             else
                Do ie=1, neg
#if 1
                   work2(ie) = SpectralWeight(neordr_t(ie,ik),ik)
#else
                   work2(ie) = SpectralWeight(ie,ik)
#endif
                End Do
             endif
             write( nfband_spectr_wght,'(4f18.10)') ( work2(ie),ie=1, neg_t )
             write( nfband_spectr_wght,* )
          End Do
          call m_Files_close_nfband_spwt

       else
          call m_Files_open_nfband_spwt(2)
          write(nfband_spectr_wght,'(A,I8)') 'num_kpoints = ', kv3 /ndim_spinor
          write(nfband_spectr_wght,'(A,I8)') 'num_bands   = ', neg_t
          write(nfband_spectr_wght,'(A,I8)') 'nspin       = ', nspin /ndim_spinor
          write(nfband_spectr_wght,*)
          Do ik=1, kv3, ndim_spinor
             write( nfband_spectr_wght,'(A,I6,A,3F16.8,A)') &
                  &             "ik = ", ik, " ( ", &
                  &                    vkxyz_refcell(ik,1:3,BUCS), " )"
             if ( noncol ) then
                Do ie=1, neg
#if 1
                   work2(ie) = SpectralWeight(neordr_t(ie,ik),ik) &
                        &     +SpectralWeight(neordr_t(ie,ik+1),ik+1)
#else
                   work2(ie) = SpectralWeight(ie,ik ) &
                        &     +SpectralWeight(ie,ik+1)
#endif
                End Do
             else
                Do ie=1, neg
#if 1
                   work2(ie) = SpectralWeight(neordr_t(ie,ik),ik)
#else
                   work2(ie) = SpectralWeight(ie,ik)
#endif
                End Do
             endif
             write( nfband_spectr_wght,'(4f18.10)') ( work2(ie),ie=1, neg_t )
             write( nfband_spectr_wght,* )
          End Do
          call m_Files_close_nfband_spwt
       endif
    endif
    deallocate( work2 )
#if 1
    deallocate( neordr_t )
#endif

    deallocate( SpectralWeight )
    if ( allocated( fsr_filtered ) ) deallocate( fsr_filtered )
    if ( allocated( fsi_filtered ) ) deallocate( fsi_filtered )

    band_unfolding_active = .false.

  contains

   subroutine hardpart
      integer :: ik, i, ib1, ia, it, lmt1, lmt2, p1, p2, iadd
      real(kind=DP) :: csum
      real(kind=DP), allocatable :: work1(:,:), work2(:,:), fs(:,:)

      allocate( work1(neg,kv3) ) ; work1 = 0.0d0

      Do ik=1, kv3, ndim_spinor
         if ( map_k(ik) /= myrank_k ) cycle

!         Do ie=ista_e, iend_e
         Do ie=1, neg
            if (map_e(ie) /= myrank_e) cycle
            ib1 = map_z(ie)

            if ( k_symmetry(ik) == GAMMA) then
               allocate(fs(nlmta,1));   fs = 0.d0
               do i = ista_fs, iend_fs
                  iadd = i - ista_fs + 1
                  fs(i,1) = fsr_filtered(ib1,iadd,ik)
               end do
               call mpi_allreduce( MPI_IN_PLACE, fs, nlmta, mpi_double_precision, &
                    &              mpi_sum, mpi_ke_world, ierr )
            else
               allocate(fs(nlmta,2));    fs = 0.d0
               do i = ista_fs, iend_fs
                  iadd = i - ista_fs + 1
                  fs(i,1) = fsr_filtered(ib1,iadd,ik)
                  fs(i,2) = fsi_filtered(ib1,iadd,ik)
               end do
               call mpi_allreduce( MPI_IN_PLACE, fs, nlmta*2, mpi_double_precision, &
                    &              mpi_sum, mpi_ke_world, ierr )
            endif
            csum = 0.0d0

            do ia = ista_atm, iend_atm
               it = ityp(ia)
               do lmt1 = 1, ilmt(it)
                  p1 = lmta(lmt1,ia)
                  do lmt2 = 1, ilmt(it)
                     p2 = lmta(lmt2,ia)

                     if ( k_symmetry(ik) == GAMMA) then
                        csum = csum +fs( p1,1 ) *fs( p2,1 ) &
                             &      *q( lmt1, lmt2, it )
                     else
                        csum = csum +( fs( p1,1 ) *fs( p2,1 ) &
                             &        +fs( p1,2 ) *fs( p2,2 ) ) &
                             &      *q( lmt1, lmt2, it )
                     endif
                  end do
               end do
            end do
            deallocate(fs)
            work1(ie,ik) = work1(ie,ik) +csum
         End Do
      End Do
      SpectralWeight = SpectralWeight +work1
      deallocate(work1)
    end subroutine hardpart

    subroutine softpart
      integer :: ik, ie, ib1, ig, ii, ierr
      real(kind=DP) :: csum, c1, c2

      Do ik=1, kv3, ndim_spinor
         if ( map_k(ik) /= myrank_k ) cycle

!         Do ie=ista_e, iend_e
         Do ie=1, neg
            if (map_e(ie) /= myrank_e) cycle

            ib1 = map_z(ie)

            csum = 0.0d0

            if ( kimg == 1 ) then
               Do ig=ista_g1k(ik), iend_g1k(ik)
                  ii = nbase(ig, ik )
                  if ( sw_force_kpt_inside_bz == ON ) then
                     if ( GVec_on_refcell(ii,ik) == 0 ) cycle
                  else
                     if ( GVec_on_refcell(ii,1) == 0 ) cycle
                  endif
                  c1 = zaj_l( ig-ista_g1k(ik)+1, ib1, ik, 1 )
                  csum = csum + c1 **2
               End Do
            else
               Do ig=ista_g1k(ik), iend_g1k(ik)
                  ii = nbase(ig, ik )
                  if ( sw_force_kpt_inside_bz == ON ) then
                     if ( GVec_on_refcell(ii,ik) == 0 ) cycle
                  else
                     if ( GVec_on_refcell(ii,1) == 0 ) cycle
                  endif
                  c1 = zaj_l( ig-ista_g1k(ik)+1, ib1, ik, 1 )
                  c2 = zaj_l( ig-ista_g1k(ik)+1, ib1, ik, kimg )
                  csum = csum + c1 **2 +c2**2
               End Do
            endif
            SpectralWeight(ie,ik) = csum
         End Do
      End Do

    end subroutine softpart

  end subroutine m_BU_calc_spectral_weight_3D

  subroutine m_BU_betar_dot_WFs_3D(nfout)
    integer, intent(in)    :: nfout
    integer :: ik
    integer :: mod_ball = ON

    Do ik=ista_k, iend_k
       call m_ES_betar_dot_Psi_4_each_k_3D( nfout, zaj_l, ista_k, iend_k, ik, &
            &                               fsr_filtered, fsi_filtered, mod_ball )
    End Do
  end subroutine m_BU_betar_dot_WFs_3D

  subroutine m_BU_set_unfolding_active( mode )
    logical, intent(in) :: mode

    if ( mode ) then
       band_unfolding_active = .true.
    else
       band_unfolding_active = .false.
    endif
  end subroutine m_BU_set_unfolding_active

  subroutine m_BU_betar_dot_WFs2(nfout)
    integer, intent(in)    :: nfout
    integer :: ik
    integer :: mod_ball = ON

    Do ik=ista_k, iend_k
       call m_ES_betar_dot_Psi_4_each_k_3D( nfout, zaj_l, ista_k, iend_k, ik, &
            &                               fsr_l, fsi_l, mod_ball )
    End Do
  end subroutine m_BU_betar_dot_WFs2

  subroutine m_BU_phir_dot_WFs_3D
    band_unfolding_active = .true.
    call m_ES_phir_dot_WFs_3D( nfout )
    band_unfolding_active = .false.

  end subroutine m_BU_phir_dot_WFs_3D

end module m_Band_Unfolding
