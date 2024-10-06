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
       &                         nrank_e

  use m_IterationNumbers,  only : nk_in_the_process
  use m_Ionic_System,   only : pos, natm, ityp, altv, iproj_group
  use m_Electronic_Structure,  only : zaj_l, neordr, fsr_l, fsi_l
  use m_ES_nonlocal,     only : zfcos, zfsin, zfcos_mpi, zfsin_mpi, &
       &                        m_ES_alloc_zfsincos, m_ES_dealloc_zfsincos, &
       &                        alloc_zfsincos_mpi, dealloc_zfsincos_mpi, &
       &                        m_ES_phir_dot_WFs, &
       &                        m_ES_betar_dot_Psi_4_each_k
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

  subroutine m_BU_calc_spectral_weight
    real(kind=DP), allocatable :: SpectralWeight(:,:)
    integer :: ik, ie, neg_t
    integer :: lun, j, k

    real(kind=DP) :: c1
    real(kind=DP), allocatable :: work(:,:), work2(:)
    integer, allocatable :: iwork(:,:), neordr_t(:,:)

    band_unfolding_active = .true.

    if ( modnrm /= SKIP ) then
       allocate( fsr_filtered(np_e,nlmta,ista_k:iend_k) )
       if( .not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
          allocate( fsi_filtered(np_e,nlmta,ista_k:iend_k));
       else
          allocate( fsi_filtered(1,1,1) )
       endif
       fsr_filtered = 0.0d0;     fsi_filtered = 0.d0
       call m_BU_betar_dot_WFs(nfout)
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
       neordr_t = iwork /nrank_e
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
                  &              vkxyz_refcell(ik+nk_in_the_process-1,1:3,BUCS), " )"
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
                  &                      vkxyz_refcell(ik,1:3,BUCS), " )"
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

! -------------------------------------------
#if 0
    if ( sw_orb_popu == ON .and. sw_calc_wf_orb_projection == ON &
         &                 .and. orb_popu_method == 2 ) then
       allocate( compr_filtered(np_e,nlmta_phi,1,ista_k:iend_k));
       if( .not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
          allocate( compi_filtered(np_e,nlmta_phi,1,ista_k:iend_k));
       else
          allocate( compi_filtered(1,1,1,1) )
       endif
       compr_filtered = 0.0d0;  compi_filtered = 0.0d0
       call m_BU_phir_dot_WFs(nfout)
       call m_BU_wd_Wfn_orb_proj
       deallocate( compr_filtered );    deallocate( compi_filtered )
    endif
#endif

    deallocate( SpectralWeight )
    if ( allocated( fsr_filtered ) ) deallocate( fsr_filtered )
    if ( allocated( fsi_filtered ) ) deallocate( fsi_filtered )

    band_unfolding_active = .false.

  contains

    subroutine hardpart
      integer :: ik, i, ib1, ia, it, lmt1, lmt2, p1, p2
      real(kind=DP) :: csum

      Do ik=1, kv3, ndim_spinor
         if ( map_k(ik) /= myrank_k ) cycle

         Do ie=ista_e, iend_e
            ib1 = map_z(ie)

            csum = 0.0d0
            do ia = 1, natm
               it = ityp(ia)
               do lmt1 = 1, ilmt(it)
                  p1 = lmta(lmt1,ia)
                  do lmt2 = 1, ilmt(it)
                     p2 = lmta(lmt2,ia)

                     if ( k_symmetry(ik) == GAMMA) then
                        csum = csum +fsr_filtered( ib1, p1, ik ) &
                             &      *fsr_filtered( ib1, p2, ik ) &
                             &      *q( lmt1, lmt2, it )
                     else
                        csum = csum +( fsr_filtered( ib1, p1, ik ) &
                             &        *fsr_filtered( ib1, p2, ik ) &
                             &        +fsi_filtered( ib1, p1, ik ) &
                             &        *fsi_filtered( ib1, p2, ik ) ) &
                             &      *q( lmt1, lmt2, it )
                     endif
                  end do
               end do
            end do
            SpectralWeight(ie,ik) = SpectralWeight(ie,ik) +csum
         End Do
      End Do

    end subroutine hardpart

    subroutine softpart
      integer :: ik, ie, ib1, ig, ii, ierr
      real(kind=DP) :: csum, c1, c2

      Do ik=1, kv3, ndim_spinor
         if ( map_k(ik) /= myrank_k ) cycle

         Do ie=ista_e, iend_e
            ib1 = map_z(ie)

            csum = 0.0d0

            if ( kimg == 1 ) then
               Do ig=1, iba(ik)
                  ii = nbase(ig, ik )
                  if ( sw_force_kpt_inside_bz == ON ) then
                     if ( GVec_on_refcell(ii,ik) == 0 ) cycle
                  else
                     if ( GVec_on_refcell(ii,1) == 0 ) cycle
                  endif
                  c1 = zaj_l( ig, ib1, ik, 1 )
                  csum = csum + c1 **2
               End Do
            else
               Do ig=1, iba(ik)
                  ii = nbase(ig, ik )
                  if ( sw_force_kpt_inside_bz == ON ) then
                     if ( GVec_on_refcell(ii,ik) == 0 ) cycle
                  else
                     if ( GVec_on_refcell(ii,1) == 0 ) cycle
                  endif
                  c1 = zaj_l( ig, ib1, ik, 1 )
                  c2 = zaj_l( ig, ib1, ik, kimg )
                  csum = csum + c1 **2 +c2**2
               End Do
            endif
            SpectralWeight(ie,ik) = csum
         End Do
      End Do

    end subroutine softpart

  end subroutine m_BU_calc_spectral_weight

  subroutine m_BU_betar_dot_WFs(nfout)
    integer, intent(in)    :: nfout
    integer :: ik

    Do ik=ista_k, iend_k
       call m_ES_betar_dot_Psi_4_each_k( nfout, zaj_l, ista_k, iend_k, ik, &
            &                            fsr_filtered, fsi_filtered )
    End Do
  end subroutine m_BU_betar_dot_WFs

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

!!    band_unfolding_active = .true.
#if 1
    Do ik=ista_k, iend_k
       call m_ES_betar_dot_Psi_4_each_k( nfout, zaj_l, ista_k, iend_k, ik, &
            &                            fsr_l, fsi_l )
    End Do
#endif
!!    band_unfolding_active = .false.

  end subroutine m_BU_betar_dot_WFs2

  subroutine m_BU_phir_dot_WFs
    band_unfolding_active = .true.
    call m_ES_phir_dot_WFs( nfout )
    band_unfolding_active = .false.

  end subroutine m_BU_phir_dot_WFs

#if 0
  subroutine m_BU_phir_dot_WFs(nfout)
    use m_ES_nonlocal,  only : G_dot_R_map, G_dot_R_mpi

    integer, intent(in)    :: nfout

    integer ia, ik, ikphig, mapmode
    integer     :: id_sname = -1

!    call tstatc0_begin('m_BU_phir_dot_WFs ',id_sname,level=1)

!    if(ipri >= 2) write(nfout,'(" -- m_BU_phir_dot_WFs --")')
    call m_ES_alloc_zfsincos(0)
    call m_ES_alloc_arai(0)
    call alloc_zfsincos_mpi()

    do ia = 1, natm
       if ( kv3/nspin == 1 ) then
          call G_dot_R_map(ia,1);      mapmode = MAPPED
       else
          call G_dot_R_mpi(ia) ! pos(ia,1:3), ngabc -> zfcos,zfsin
          mapmode = NOTMAPPED
       endif

       do ik = 1, kv3, af+1
          if(map_k(ik) /= myrank_k) cycle         ! MPI

          ikphig = (ik-1)/nspin + 1
          call m_ES_phir_dot_WFs_4_lmta_k( ista_k, iend_k, ik, zaj_l, ia, ikphig, &
               &                           phig, 1, compr_filtered, compi_filtered,  &
               &                           mapmode )
       end do
    end do
    call dealloc_zfsincos_mpi()
    call m_ES_dealloc_arai()
    call m_ES_dealloc_zfsincos()

!    if ( modnrm /= SKIP ) then
       if ( orb_popu_method == 2 ) call add_overlap_phirt_psirpw
!    endif

!    call tstatc0_end(id_sname)
  end subroutine m_BU_phir_dot_WFs

  subroutine add_overlap_phirt_psirpw
    integer :: ik, ia, it, lmt1, lmt2, p, q, s
    integer :: il1, il2, im1, im2, t1, t2
    integer :: ikphig
    real(kind=DP) :: ctmp, c1

    Do ik=1, kv3
       if ( map_k(ik) /= myrank_k ) cycle

       ikphig = ( ik -1 )/nspin +1

       Do ia=1, natm
          it = ityp(ia)

          Do lmt1=1, ilmt_phi(it)
             il1=ltp_phi(lmt1,it);  im1 = mtp_phi(lmt1,it);  t1 = taup_phi(lmt1,it)
             p = lmta_phi(lmt1,ia);  s = lmtt_phi(lmt1,it)
             c1 = sqrt( norm_phig( s,ikphig ) )

             Do lmt2=1, ilmt(it)
                il2=ltp(lmt2,it);  im2 = mtp(lmt2,it);  t2 = taup(lmt2,it)
                if ( il1 /= il2 ) cycle
                if ( im1 /= im2 ) cycle

                q = lmta(lmt2,ia)
                ctmp = q_phirt_pw(il1,t1,t2,it) /c1

                if ( k_symmetry(ik) == GAMMA ) then
                   compr_filtered(:,p,1,ik) &
                        &       = compr_filtered(:,p,1,ik) &
                        &        +ctmp *fsr_filtered(:,q,ik) *sqrt(2.0d0)  ! why ?
                else
                   compr_filtered(:,p,1,ik) &
                        &    = compr_filtered(:,p,1,ik) +ctmp *fsr_filtered(:,q,ik)
                   compi_filtered(:,p,1,ik) &
                        &    = compi_filtered(:,p,1,ik) +ctmp *fsi_filtered(:,q,ik)
                endif
             End Do
          End Do
       End Do
    End Do
  end subroutine add_overlap_phirt_psirpw

  subroutine m_BU_wd_Wfn_orb_proj
    integer :: neg_t
    real(kind=DP), allocatable :: compr(:,:,:,:), compi(:,:,:,:), norm_phig_mpi(:,:)

    if ( ekmode == OFF ) return

    allocate(compr(neg,nlmta_phi,1,kv3));  compr = 0.d0
    allocate(compi(neg,nlmta_phi,1,kv3));  compi = 0.d0

    if(.not.allocated(norm_phig_mpi)) allocate(norm_phig_mpi(nlmtt_phi,kv3/nspin))
    norm_phig_mpi=0.d0

    call set_array_compri_etc

    neg_t = neg -num_extra_bands

    if ( mype == 0 ) then
       if ( nk_in_the_process == 1 ) then
          call m_Files_open_nfwfk_orb_proj(2)

          write(nfwfk_orb_proj,'(A,I8)') 'num_kpoints = ', kv3_ek /ndim_spinor
          write(nfwfk_orb_proj,'(A,I8)') 'num_bands   = ', neg_t
          write(nfwfk_orb_proj,'(A,I8)') 'nspin       = ', nspin

          write(nfwfk_orb_proj,*)
          write(nfwfk_orb_proj,'(A,I8)') 'num of orbitals = ', nlmta_phi
          write(nfwfk_orb_proj,*)
       else
          call m_Files_open_nfwfk_orb_proj(3)
       endif
    endif

    if ( SpinOrbit_Mode /= Neglected .and. wf_orb_proj_print_format == 1 ) then
       call case_with_j
    else
       call case_ordinal    ! EXPERIMENTAL
    endif
    if ( mype == 0 ) call m_Files_close_nfwfk_orb_proj

    deallocate( compr ); deallocate( compi ); deallocate( norm_phig_mpi )

  contains

    subroutine set_array_compri_etc
      integer :: ik, ie, ib, iksnl
      integer :: iorb, lmt
      integer :: ia, il, im, tau, is
      real(kind=DP), allocatable :: compr_mpi(:,:,:,:), compi_mpi(:,:,:,:), &
           &                        norm_phig_mpi2(:,:)
      real(kind=DP), allocatable :: porb(:)

      do ik = 1, kv3
         if(map_k(ik) /= myrank_k) cycle
         iksnl = (ik-1)/nspin + 1

         do ie = ista_e, iend_e, istep_e
            ib = map_z(ie)
            compr(ie,1:nlmta_phi,1,ik) = compr_filtered(ib,1:nlmta_phi,1,ik)
            compi(ie,1:nlmta_phi,1,ik) = compi_filtered(ib,1:nlmta_phi,1,ik)
         end do
         norm_phig_mpi(1:nlmtt_phi,iksnl)  = norm_phig(1:nlmtt_phi,iksnl)
      end do

      if ( npes >1 ) then
         allocate( compr_mpi( neg, nlmta_phi, 1, kv3 ) ); compr_mpi = 0.0d0
         allocate( compi_mpi( neg, nlmta_phi, 1, kv3 ) ); compi_mpi = 0.0d0
         allocate( norm_phig_mpi2( nlmtt_phi, kv3/nspin ) )
         call mpi_allreduce( compr, compr_mpi, neg*nlmta_phi*1*kv3, &
              &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
         call mpi_allreduce( compi, compi_mpi, neg*nlmta_phi*1*kv3, &
              &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
         call mpi_allreduce( norm_phig_mpi, norm_phig_mpi2, nlmtt_phi*kv3/nspin, &
              &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
         compr = compr_mpi;   compi = compi_mpi

         norm_phig_mpi = norm_phig_mpi2 /dble(nrank_e)

         deallocate( compr_mpi ); deallocate( compi_mpi );
         deallocate( norm_phig_mpi2 )
      end if
    end subroutine set_array_compri_etc

    subroutine case_ordinal
      integer :: ik, iksnl, iorb, ia, is, ib
      integer :: il, im, tau, lmtt
      real(kind=DP), allocatable :: porb(:)

      if ( mype /= 0 ) return

      allocate( porb(neg) )

      do ik = 1, kv3, ndim_spinor
         iksnl = (ik-1)/nspin + 1

         write(nfwfk_orb_proj,'(A,I6,A,3F16.8,A)') 'ik = ', &
              &            ik +nk_in_the_process -1, " ( ", &
              &            vkxyz_refcell(ik,1:3,BUCS), " )"

         do iorb = 1,nlmta_phi
            call m_PP_tell_iorb_ia_l_m_tau(iorb,ia,il,im,tau)
            if ( iproj_group(ia) == 0) cycle

            write(nfwfk_orb_proj,'(I5,3I3,A)') ia, il-1, im, tau, ' : ia, l, m, tau'

            call m_PP_tell_iorb_lmtt(iorb,lmtt)

            porb = 0.0d0
            Do is=1, ndim_spinor
               do ib = 1, neg
                  porb(ib) = porb(ib) &
                       & + ( compr(ib,iorb,1,ik+is-1)**2 &
                       &    +compi(ib,iorb,1,ik+is-1)**2 ) &
                       &     *( 1.d0+qorb(iorb)/norm_phig_mpi(lmtt,iksnl) )
               end do
            End Do
!            write(nfwfk_orb_proj,'(4F18.10)') ( porb( neordr(ib,ik) ), ib=1, neg_t )
            write(nfwfk_orb_proj,'(4E18.10)') ( porb( neordr(ib,ik) ), ib=1, neg_t )
         end do
         write(nfwfk_orb_proj,*)
      end do
      deallocate( porb )

    end subroutine case_ordinal

    subroutine case_with_j
      integer :: ik, iksnl, ia, ig, ii, it, ilp, ll, tau, ip
      integer :: iorb, lmtt1, ib, m1
      real(kind=DP) :: c1, c2
      complex(kind=CMPLDP) :: z1, z2

      real(kind=DP), allocatable :: porb(:)
      complex(kind=CMPLDP), allocatable :: zcomp(:,:,:)

! ----------------
      if ( mype /= 0 ) return

      allocate( porb(neg) )

      Do ik=1, kv3, ndim_spinor
         write(nfwfk_orb_proj,'(A,I6,A,3F16.8,A)') 'ik = ', &
              &            ik +nk_in_the_process -1, " ( ", &
              &            vkxyz_refcell(ik,1:3,BUCS), " )"

         iksnl = ( ik -1 )/nspin +1

         Do ia=1, natm
            ig = iproj_group(ia)
            if ( ig == 0 ) cycle

            do ii=1,num_proj_elems(ig)
               ip = proj_group( ii, ig )
               it = proj_attribute(ip)%ityp
               ilp = proj_attribute(ip)%l +1
               ll = proj_attribute(ip)%l
               tau = proj_attribute(ip)%t
!
               allocate( zcomp( -ll:ll, neg, ndim_spinor ) ); zcomp = 0.0d0
               call tranform_compri_r2c_sph( ia, it, ll, tau, ik, &
                    &                        compr, compi, zcomp )

               if ( ll == 0 ) then
                  call find_iorb_from_lmt( ia, it, ll, 1, tau, iorb )
                  call m_PP_tell_iorb_lmtt( iorb, lmtt1 )

                  write(nfwfk_orb_proj,'(I5,F7.2,I3,F7.2,I3,A)') &
                       &         ia, ll+0.5d0, ll, 0.5d0, tau, ' : ia, j, l, mj, tau'
                  Do ib=1, neg
                     z1 = cmplx( compr(ib,iorb,1,ik),   compi(ib,iorb,1,ik) )
                     z2 = cmplx( compr(ib,iorb,1,ik+1), compi(ib,iorb,1,ik+1) )
                     porb(ib) = ( z1*conjg(z1) +z2*conjg(z2) ) &
                          &     *( 1.d0+qorb(iorb)/norm_phig_mpi( lmtt1, iksnl ) )
                  End Do
                  write(nfwfk_orb_proj,'(4F18.10)') &
                       &               ( porb( neordr(ib,ik) ), ib=1, neg_t )
               else
! j_up
                  Do m1=-ll -1, ll
                     c1 = dble( ll +m1 + 1 ) / dble( 2 *ll +1 )
                     c2 = dble( ll -m1 )     / dble( 2 *ll +1 )
                     c1 = sqrt(c1);  c2 = sqrt(c2)

                     call find_iorb_from_lmt( ia, it, ll, 1, tau, iorb )
                     call m_PP_tell_iorb_lmtt( iorb, lmtt1 )

                     write(nfwfk_orb_proj,'(I5,F7.2,I3,F7.2,I3,A)') &
                          &         ia, ll+0.5d0, ll, m1+0.5d0, tau, &
                          &             ' : ia, j, l, mj, tau'
                     Do ib=1, neg
                        z1 = 0.0d0
                        if ( m1 > -ll -1 ) z1 = z1 +c1 *zcomp( m1,    ib, 1 )
                        if ( m1 < ll     ) z1 = z1 +c2 *zcomp( m1 +1, ib, 2 )
                        porb(ib) = z1 *conjg(z1) &
                             &     *( 1.d0+qorb(iorb)/norm_phig_mpi(lmtt1,iksnl) )
                     End Do
                     write(nfwfk_orb_proj,'(4F18.10)') &
                          &         ( porb( neordr(ib,ik) ), ib=1, neg_t )
                  End Do
! j_down
                  Do m1=-ll+1, ll
                     c1 = dble( ll -m1 + 1 ) / dble( 2 *ll +1 )
                     c2 = dble( ll +m1 )     / dble( 2 *ll +1 )
                     c1 = sqrt(c1);  c2 = -sqrt(c2)

                     call find_iorb_from_lmt( ia, it, ll, 1, tau, iorb )
                     call m_PP_tell_iorb_lmtt( iorb, lmtt1 )

                     write(nfwfk_orb_proj,'(I5,F7.2,I3,F7.2,I3,A)') &
                          &         ia, ll-0.5d0, ll, m1-0.5d0, tau, &
                          &             ' : ia, j, l, mj, tau'
                     Do ib=1, neg
                        z1 = c1 *zcomp( m1-1, ib, 1 ) +c2 *zcomp( m1, ib, 2 )
                        porb(ib) = z1 *conjg(z1) &
                             &     *( 1.d0+qorb(iorb)/norm_phig_mpi(lmtt1,iksnl) )
                     End Do
                     write(nfwfk_orb_proj,'(4F18.10)') &
                          &         ( porb( neordr(ib,ik) ), ib=1, neg_t )
                  End Do
               end if
               deallocate( zcomp )

            End Do
         End Do
         write(nfwfk_orb_proj,*)

      End Do
      deallocate( porb )

    end subroutine case_with_j

    subroutine tranform_compri_r2c_sph( ia, it, ll, tau, ik, &
         &                              compr, compi, zcomp )
      integer, intent(in) :: ia, it, ll, tau, ik
      real(kind=DP), intent(in) :: compr( neg, nlmta_phi, 1, kv3 )
      real(kind=DP), intent(in) :: compi( neg, nlmta_phi, 1, kv3 )
      complex(kind=CMPLDP), intent(out) :: zcomp( -ll:ll, neg, ndim_spinor )

      integer :: m1, m2, ib, iorb
      complex(kind=CMPLDP) :: z1
      complex(kind=CMPLDP) :: ztmp( ndim_spinor )

      Do m1=-ll, ll
         Do m2=1, 2*ll +1
            call find_iorb_from_lmt( ia, it, ll, m2, tau, iorb )
            Do ib=1, neg
               ztmp(1) = cmplx( compr(ib,iorb,1,ik),   compi(ib,iorb,1,ik)   )
               ztmp(2) = cmplx( compr(ib,iorb,1,ik+1), compi(ib,iorb,1,ik+1) )

               if ( ll == 0 ) z1 = MatU_ylm_RC_L0( m2, m1 )
               if ( ll == 1 ) z1 = MatU_ylm_RC_L1( m2, m1 )
               if ( ll == 2 ) z1 = MatU_ylm_RC_L2( m2, m1 )
               if ( ll == 3 ) z1 = MatU_ylm_RC_L3( m2, m1 )

               zcomp(m1,ib,:) = zcomp(m1,ib,:) +z1 *ztmp(:)
            End Do
         End Do
      End Do
    end subroutine tranform_compri_r2c_sph

    subroutine find_iorb_from_lmt( ia, it, ll, mm, tau, iorb )
      integer, intent(in) :: ia, it, ll, mm, tau
      integer, intent(out) :: iorb

      integer :: lmt1, l1, m1, t1

      iorb = 0
      Do lmt1=1, ilmt_phi(it)
         l1 = ltp_phi(lmt1,it); m1 = mtp_phi(lmt1,it);  t1 = taup_phi(lmt1,it)
         if ( l1 == ll +1 .and. m1 == mm .and. t1 == tau ) then
            exit
         endif
      ENd Do
      iorb = lmta_phi( lmt1,ia )

    end subroutine find_iorb_from_lmt

  end subroutine m_BU_wd_Wfn_orb_proj
#endif

end module m_Band_Unfolding
