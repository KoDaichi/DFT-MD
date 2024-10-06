module m_SpinOrbit_SecondVariation
  use m_Const_Parameters,  only : DP, SP, CMPLDP, PARABOLIC, TETRAHEDRON, COLD, &
       &                          FERMI_DIRAC, MP, GAMMA, &
       &                          Neglected, ByPawPot, ByProjector, ZeffApprox, &
       &                          ReadFromPP, OFF, ON, BUCS, &
       &                          FIXED_CHARGE_CONTINUATION, FIXED_CHARGE, ONE_BY_ONE
  use m_Control_Parameters,  only : neg, nspin, m_CtrlP_way_of_smearing, &
       &                            ndim_spinor, ndim_chgpot, kimg, ekmode, &
       &                            SpinOrbit_mode, icond, fixed_charge_k_parallel, &
       &                            neg_is_enlarged, num_extra_bands, &
       &                            precision_WFfile, sw_modified_kpoint_increment, &
       &                            printable, sw_write_zaj_socsv
  use m_Files,              only : nfout
  use m_Kpoints,          only : kv3, k_symmetry, vkxyz, kv3_ek, qwgt, vkxyz_ek

  use m_Electronic_Structure,  only : eko_l, fsr_l, fsi_l, occup_l, nrvf_ordr, efermi, &
       &                              zaj_l, metalic_system, vbm, eko_ek
  use m_Ionic_System,       only : ityp, natm
  use m_PseudoPotential,   only : nlmt, ltp, taup, nlmta, ilmt, lmta

  use m_Parallelization,  only : map_k, map_e, map_z, myrank_k, myrank_e, mpi_k_world, &
       &                         ista_e, iend_e, istep_e, np_e, npes, ierr, &
       &                         MPI_CommGroup, mype, ista_k, iend_k, map_ek, &
       &                         nrank_k, nis_kv3_ek

  use m_ES_Occup,            only : m_ESoc_fermi_parabolic, m_ESoc_fermi_tetrahedron, &
       &                            m_ESoc_fermi_ColdSmearing, m_ESoc_fermi_Dirac, &
       &                            m_ESoc_methfessel_paxton
  use m_ES_occup,            only : m_ESoc_fermi_parabolic_ek, &
       &                            m_ESoc_fermi_tetra_ek, &
       &                            m_ESoc_fermi_Dirac_ek, &
       &                            m_ESoc_mp_ek

  use m_SpinOrbit_Potential,     only : m_SO_alloc_dsoc, dsoc, &
       &                                m_SO_set_Dsoc_potential1, &
       &                                m_SO_calc_MatLS_orb_s_to_f, &
       &                                m_SO_diagonalize_MatLS, &
       &                                m_SO_set_MatU_ylm_RC, &
       &                                m_SO_alloc_Mat_SOC_strenth, &
       &                                m_SO_set_Dsoc_potential2, &
       &                                m_SO_dealloc_Mat_SOC_strenth, &
       &                                m_SO_dealloc_dsoc
  use m_SpinOrbit_RadInt,       only : m_SO_calc_SOC_strength_pawpot, &
       &                               m_SO_calc_SOC_strength_zeff, &
       &                               m_SO_check_mode_Builtin, &
       &                               m_SO_check_mode_Pawpot, &
       &                               m_SO_check_mode_Zeff
  use m_SpinOrbit_FromFile,     only : m_SO_set_SOC_strength_from_PP

  use m_IterationNumbers,     only : nk_in_the_process, nkgroup, &
       &                             first_kpoint_in_this_job, nk_converged

  use m_PlaneWaveBasisSet,    only : kg1
  use mpi

  implicit none
!  include 'mpif.h'
  integer istatus(mpi_status_size)
!
  integer, parameter :: ndim_spinor_socsv = 2
  integer :: neg_doubled
  integer :: ndim_spinor_prev, ndim_chgpot_prev

  real(kind=DP), allocatable :: eko_ek_socsv(:,:)
  real(kind=DP), allocatable :: zaj_band1_kpt1(:,:,:,:)
  real(kind=DP), allocatable :: zaj_band1_kpt2(:,:,:,:)
  real(kind=DP), allocatable :: zaj_band2_kpt1(:,:,:,:)
  real(kind=DP), allocatable :: zaj_band2_kpt2(:,:,:,:)
!
contains

  subroutine m_SO_init_second_variation
    ndim_spinor_prev = ndim_spinor;   ndim_chgpot_prev = ndim_chgpot
    ndim_spinor = ndim_spinor_socsv;    ndim_chgpot = ndim_spinor **2

    call m_SO_check_mode_Pawpot
    call m_SO_check_mode_Zeff

    if ( SpinOrbit_Mode /= Neglected ) then
       call m_SO_set_MatU_ylm_RC
       call m_SO_calc_MatLS_orb_s_to_f
       call m_SO_diagonalize_MatLS
    endif
    if ( SpinOrbit_Mode == ByPawPot ) then
       call m_SO_alloc_Mat_SOC_strenth
       call m_SO_calc_SOC_strength_pawpot
       call m_SO_alloc_Dsoc
       call m_SO_set_Dsoc_potential2
    endif
    if ( SpinOrbit_Mode == ZeffApprox ) then
       call m_SO_alloc_Mat_SOC_strenth
       call m_SO_calc_SOC_strength_zeff
       call m_SO_alloc_Dsoc
       call m_SO_set_Dsoc_potential2
    endif
    if ( SpinOrbit_Mode == ReadFromPP ) then
       call m_SO_alloc_Mat_SOC_strenth
       call m_SO_set_SOC_strength_from_PP
       call m_SO_alloc_Dsoc
       call m_SO_set_Dsoc_potential2
    endif
    if ( SpinOrbit_Mode == ByProjector ) then
       call m_SO_alloc_Dsoc
       call m_SO_set_Dsoc_potential1
    endif

    if ( ekmode == ON ) then
       allocate( eko_ek_socsv(neg,kv3_ek) );  eko_ek_socsv = 0.0d0
    endif
    ndim_spinor = ndim_spinor_prev;       ndim_chgpot = ndim_chgpot_prev

  end subroutine m_SO_init_second_variation

  subroutine m_SO_finalize_second_variation
    if ( SpinOrbit_Mode == ByPawPot .or. SpinOrbit_Mode == ZeffApprox &
         &                          .or. SpinOrbit_Mode == ReadFromPP ) then
       call m_SO_dealloc_Mat_SOC_strenth
       call m_SO_dealloc_Dsoc
    endif
    if ( SpinOrbit_Mode == ByProjector ) then
       call m_SO_dealloc_Dsoc
    endif
    if ( ekmode == ON ) deallocate( eko_ek_socsv )

  end subroutine m_SO_finalize_second_variation

  subroutine m_SO_calc_band_energy_socsv
    integer :: ik
    real(kind=DP) :: efermi_save, eband
    real(kind=DP), allocatable :: eko_save(:,:), occ_save(:,:)
    complex(kind=CMPLDP), allocatable :: MatH(:,:)

    if ( SpinOrbit_mode == Neglected ) return

    neg_doubled = neg *ndim_spinor_socsv
    allocate( MatH( neg_doubled, neg_doubled ) ); MatH = 0.0d0
    allocate( eko_save( np_e, ista_k:iend_k) ); eko_save = 0.0d0
    allocate( occ_save( np_e, ista_k:iend_k) ); occ_save = 0.0d0
    eko_save = eko_l; occ_save = occup_l;  efermi_save = efermi

    if ( ekmode == OFF ) then
       if ( sw_write_zaj_socsv == ON ) then
          allocate( zaj_band1_kpt1(kg1,np_e,ista_k:iend_k,kimg) )
          allocate( zaj_band2_kpt1(kg1,np_e,ista_k:iend_k,kimg) )
          allocate( zaj_band1_kpt2(kg1,np_e,ista_k:iend_k,kimg) )
          allocate( zaj_band2_kpt2(kg1,np_e,ista_k:iend_k,kimg) )
       endif
!
       Do ik=1, kv3, nspin
          if ( map_k(ik) /= myrank_k ) cycle
          MatH = 0.0d0
          call set_contrib_from_dsoc
          call add_diagonal_elements
          call renew_level_socsv( ik, MatH, eko_l )
          if ( sw_write_zaj_socsv == ON ) call renew_wfn( ik, MatH )
       End Do

       write(nfout,*)
       write(nfout,*) "<<< Spin-orbit by Second Variation "

       call FermiEnergyLevel
       call calc_bandenergy_socsv( map_k, occup_l, eko_l, eband )
       write(nfout,*) "* band energy  = ", eband
       write(nfout,*) ">>>"
       write(nfout,*)

       call m_SO_wd_EigenValues( .true., .false. )

       if ( sw_write_zaj_socsv == ON ) then
          call m_SO_wd_WFs()
          deallocate( zaj_band1_kpt1 );   deallocate( zaj_band2_kpt1 )
          deallocate( zaj_band1_kpt2 );   deallocate( zaj_band2_kpt2 )
       endif

    else
       Do ik=1, kv3, nspin
          if ( map_k(ik) /= myrank_k ) cycle
          MatH = 0.0d0
          call set_contrib_from_dsoc
          call add_diagonal_elements
          call renew_level_socsv( ik, MatH, eko_l )
          call m_SO_cp_eko_l_to_eko_ek_socsv
       End Do
       call m_SO_wd_EigenValues( .false., .true. )
    endif
!
    eko_l = eko_save;   occup_l = occ_save;  efermi = efermi_save
    deallocate( eko_save ); deallocate( occ_save )

  contains

    subroutine renew_wfn( ik, Coeff )
      integer, intent(in) :: ik
      complex(kind=CMPLDP), intent(in) :: Coeff(neg_doubled,neg_doubled)

      integer :: ib, jb, ig, ierr
      complex(kind=CMPLDP) :: z1, z2, zsum1, zsum2, zsum3, zsum4
      real(kind=DP), allocatable :: zaj_wk(:,:,:)

      allocate( zaj_wk(kg1,neg_doubled,kimg) ); zaj_wk = 0.0d0
      if ( nspin==1 ) then
         Do ib=ista_e, iend_e, istep_e
            zaj_wk(:,ib,    1) = zaj_l(:,map_z(ib),ik,1)
            zaj_wk(:,ib+neg,1) = zaj_l(:,map_z(ib),ik,1)
            if (kimg==2) then
               zaj_wk(:,ib,    2) = zaj_l(:,map_z(ib),ik,2)
               zaj_wk(:,ib+neg,2) = zaj_l(:,map_z(ib),ik,2)
            endif
         End Do
      else
         Do ib=ista_e, iend_e, istep_e
            zaj_wk(:,ib,    1) = zaj_l(:,map_z(ib),ik,  1)
            zaj_wk(:,ib+neg,1) = zaj_l(:,map_z(ib),ik+1,1)
            if (kimg==2) then
               zaj_wk(:,ib,    2) = zaj_l(:,map_z(ib),ik,  2)
               zaj_wk(:,ib+neg,2) = zaj_l(:,map_z(ib),ik+1,2)
            endif
         End Do
      endif

      if ( npes > 1 ) then
         call mpi_allreduce( MPI_IN_PLACE, zaj_wk, kg1 *neg_doubled *kimg,&
               &             mpi_double_precision, mpi_sum, &
               &             mpi_k_world(myrank_k),ierr )
      endif

      if ( nspin==1 ) then
         Do ib=ista_e, iend_e, istep_e
            Do ig=1, kg1
               zsum1 = 0.0d0;  zsum2 = 0.0d0;  zsum3 = 0.0d0;  zsum4 = 0.0d0
               Do jb=1, neg
                  z1 = dcmplx( zaj_wk(ig,jb,    1), zaj_wk(ig,jb,    2) )
                  z2 = dcmplx( zaj_wk(ig,jb+neg,1), zaj_wk(ig,jb+neg,2) )
                  zsum1 = zsum1 +coeff(jb,    ib) *z1
                  zsum2 = zsum2 +coeff(jb+neg,ib) *z2
                  zsum3 = zsum3 +coeff(jb,    ib+neg) *z1
                  zsum4 = zsum4 +coeff(jb+neg,ib+neg) *z2
               End Do
               zaj_band1_kpt1(ig,map_z(ib),ik,1) = real(zsum1)
               zaj_band1_kpt1(ig,map_z(ib),ik,2) = aimag(zsum1)
               zaj_band1_kpt2(ig,map_z(ib),ik,1) = real(zsum2)
               zaj_band1_kpt2(ig,map_z(ib),ik,2) = aimag(zsum2)
               zaj_band2_kpt1(ig,map_z(ib),ik,1) = real(zsum3)
               zaj_band2_kpt1(ig,map_z(ib),ik,2) = aimag(zsum3)
               zaj_band2_kpt2(ig,map_z(ib),ik,1) = real(zsum4)
               zaj_band2_kpt2(ig,map_z(ib),ik,2) = aimag(zsum4)
            End Do
         End Do
      else
         Do ib=ista_e, iend_e, istep_e
            Do ig=1, kg1
               zsum1 = 0.0d0;  zsum2 = 0.0d0;  zsum3 = 0.0d0;  zsum4 = 0.0d0
               Do jb=1, neg
                  z1 = dcmplx( zaj_wk(ig,jb,    1), zaj_wk(ig,jb,    2) )
                  z2 = dcmplx( zaj_wk(ig,jb+neg,1), zaj_wk(ig,jb+neg,2) )
                  zsum1 = zsum1 +coeff(jb,    ib) *z1
                  zsum2 = zsum2 +coeff(jb+neg,ib) *z2
                  zsum3 = zsum3 +coeff(jb,    ib+neg) *z1
                  zsum4 = zsum4 +coeff(jb+neg,ib+neg) *z2
               End Do
               zaj_band1_kpt1(ig,map_z(ib),ik,  1) = real(zsum1)
               zaj_band1_kpt1(ig,map_z(ib),ik,  2) = aimag(zsum1)
               zaj_band1_kpt1(ig,map_z(ib),ik+1,1) = real(zsum2)
               zaj_band1_kpt1(ig,map_z(ib),ik+1,2) = aimag(zsum2)
               zaj_band2_kpt1(ig,map_z(ib),ik,  1) = real(zsum3)
               zaj_band2_kpt1(ig,map_z(ib),ik,  2) = aimag(zsum3)
               zaj_band2_kpt1(ig,map_z(ib),ik+1,1) = real(zsum4)
               zaj_band2_kpt1(ig,map_z(ib),ik+1,2) = aimag(zsum4)
            End Do
         End Do
      endif
      deallocate( zaj_wk )
    end subroutine renew_wfn

    subroutine renew_level_socsv( ik, MatH, eko_l )
      integer, intent(in) :: ik
      complex(kind=CMPLDP), intent(inout) :: MatH( neg_doubled, neg_doubled )
      real(kind=DP), intent(out) :: eko_l(np_e,ista_k:iend_k)

      integer :: matsize, lwork, info
      integer :: is, ib, itmp

      real(kind=DP), allocatable :: rwork(:), eigvals(:), eko_wk(:,:)
      complex(kind=CMPLDP), allocatable :: work(:)

      matsize = neg_doubled;   lwork = 2*matsize
      allocate( work( lwork ) );  allocate( rwork( 3*matsize-2 ) )
      allocate( eigvals( neg_doubled ) ); eigvals = 0.0d0

      call zheev( 'V', 'U', matsize,  MatH,  matsize, eigvals, work, lwork, rwork, info )
      if ( info /= 0 ) then
         write(*,*)'Error in zheev';  stop
      endif
      deallocate( work ); deallocate( rwork )
!
      allocate( eko_wk( neg, ndim_spinor_socsv ) ); eko_wk = 0.0d0
      Do ib=1, neg
         Do is=1, ndim_spinor_socsv
            itmp = ( ib -1 )*ndim_spinor_socsv +is
            eko_wk( ib, is ) = eigvals( itmp )
         End do
      End Do

      Do is=1, ndim_spinor_socsv
         itmp = 1
         if ( nspin == 2 ) itmp = is

         Do ib=ista_e, iend_e, istep_e
            eko_l(map_z(ib),ik+itmp-1) = eko_wk(ib,is)
         End do
      End Do
      deallocate( eigvals );   deallocate( eko_wk )
!
    end subroutine renew_level_socsv

    subroutine calc_bandenergy_socsv( map_k, occup_l, eko_l, eband )
      integer, intent(in) :: map_k( kv3 )
      real(kind=DP), intent(in) :: occup_l(np_e,ista_k:iend_k)
      real(kind=DP), intent(in) :: eko_l(np_e,ista_k:iend_k)
      real(kind=DP), intent(out) :: eband
!
      integer             :: ib, ik, ierr
      real(kind=DP)       :: eband_mpi

      eband = 0.d0
      do ik =1, kv3
         if (map_k(ik) /= myrank_k) cycle
         do ib = 1, np_e
            eband = eband + occup_l(ib,ik)*eko_l(ib,ik)
         end do
      end do

      if (npes > 1) then
         call mpi_allreduce( eband, eband_mpi, 1, mpi_double_precision, &
              &              mpi_sum, MPI_CommGroup,ierr )
         eband = eband_mpi
      end if
      eband = 2.0d0 *eband /dble(kv3)

    end subroutine calc_bandenergy_socsv

    subroutine add_diagonal_elements
      integer :: is, ib1, itmp
      real(kind=DP), allocatable :: eko_wk(:,:), eko_mpi(:,:)

      allocate( eko_wk( neg, ndim_spinor_socsv ) ); eko_wk = 0.0d0
      Do is=1, ndim_spinor_socsv
         itmp = 1
         if ( nspin == 2 ) itmp = is

         Do ib1=ista_e, iend_e, istep_e
            eko_wk( ib1, is ) = eko_l( map_z(ib1),ik +itmp-1 )
         End do
      End Do

      if ( npes > 1 ) then
         allocate( eko_mpi( neg, ndim_spinor_socsv ) ); eko_mpi = 0.0d0
         call mpi_allreduce( eko_wk, eko_mpi, neg*ndim_spinor_socsv, &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k), ierr )
         eko_wk = eko_mpi
         deallocate( eko_mpi )
      endif
!
      Do is=1, ndim_spinor_socsv
         Do ib1=1, neg
            itmp = ( is -1 )*neg +ib1
            MatH( itmp,itmp ) = MatH( itmp,itmp ) +eko_wk( ib1,is )
         End do
      End Do
      deallocate( eko_wk )

    end subroutine add_diagonal_elements

    subroutine set_contrib_from_dsoc
      integer :: ib1, ib2, is1, is2, itmp1, itmp2
      integer :: ia, it, istmp, jb1, jb2
      integer :: lmt1, lmt2, il1, il2, it1, it2, lmta1, lmta2
      complex(kind=CMPLDP) :: z1, wf1, wf2
      complex(kind=CMPLDP), allocatable :: wk_fsri(:,:), MatH_mpi(:,:)

      allocate( wk_fsri( nlmta, ndim_spinor_socsv ) ); wk_fsri = 0.0d0
      Do ib2=1, neg
         wk_fsri = 0.0d0

         if ( map_e(ib2) == myrank_e ) then
            Do is1=1, ndim_spinor_socsv
               itmp1 = 1
               if ( nspin == 2 ) itmp1 = is1

               if ( kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2 ) then
                  wk_fsri(:,is1) = fsr_l( map_z(ib2),:,ik+itmp1-1 )
               else
                  wk_fsri(:,is1) = dcmplx( fsr_l( map_z(ib2),:,ik+itmp1-1 ), &
                       &                   fsi_l( map_z(ib2),:,ik+itmp1-1 ) )
               endif
            End do
         endif
         call mpi_bcast( wk_fsri, 2*nlmta*ndim_spinor_socsv, mpi_double_precision, &
              &          map_e(ib2), mpi_k_world(myrank_k), ierr )
         !

         DO ib1=ista_e, iend_e, istep_e

            Do ia=1, natm
               it = ityp(ia)

               Do lmt1=1, ilmt(it)
                  il1 = ltp(lmt1,it); it1 = taup(lmt1,it)
                  lmta1 = lmta( lmt1,ia )

                  Do lmt2=1, ilmt(it)
                     il2 = ltp(lmt2,it); it2 = taup(lmt2,it)
                     lmta2 = lmta( lmt2,ia )

                     Do is1=1, ndim_spinor_socsv
                        itmp1 = 1
                        if ( nspin == 2 ) itmp1 = is1

                        Do is2=1, ndim_spinor_socsv
                           itmp2 = 1
                           if ( nspin == 2 ) itmp2 = is2

                           if ( kv3/nspin == 1 .and. k_symmetry(1) == GAMMA &
                                &              .and. kimg == 2 ) then
                              wf1 = fsr_l( map_z(ib1), lmta1, ik +itmp1 -1 )
                           else
                              wf1 = dcmplx( fsr_l( map_z(ib1), lmta1, ik +itmp1 -1 ), &
                                   &        fsi_l( map_z(ib1), lmta1, ik +itmp1 -1 ) )
                           endif
                           wf2 = wk_fsri( lmta2, is2 )

                           istmp = ( is1 -1 )*ndim_spinor_socsv +is2
                           z1 = dconjg(wf1) *dsoc( lmt1, lmt2, ia, istmp ) *wf2

                           jb1 = neg *( is1 -1 )+ib1
                           jb2 = neg *( is2 -1 )+ib2
                           MatH(jb1,jb2) = MatH (jb1,jb2) +z1
                        End Do
                     End do
                  End Do

               End Do
            End do
         End DO
      End Do
!
      if ( npes > 1 ) then
         allocate( MatH_mpi( neg_doubled, neg_doubled ) ); MatH_mpi = 0.0d0
         call mpi_allreduce( MatH, MatH_mpi, 2*neg_doubled**2, &
               &             mpi_double_precision, mpi_sum, &
               &             mpi_k_world(myrank_k), ierr )
         MatH = MatH_mpi
         deallocate( MatH_mpi )
      endif
      deallocate( wk_fsri )

    end subroutine set_contrib_from_dsoc

  end subroutine m_SO_calc_band_energy_socsv

  subroutine FermiEnergyLevel()
    integer :: way_of_smearing

    way_of_smearing = m_CtrlP_way_of_smearing()
    if(way_of_smearing == PARABOLIC) then
       call m_ESoc_fermi_parabolic(nfout)
    else if(way_of_smearing == TETRAHEDRON) then
       call phase_error_with_msg(nfout,"Not supported",__LINE__,__FILE__)
       call m_ESoc_fermi_tetrahedron(nfout)
    else if(way_of_smearing == COLD) then
       call m_ESoc_fermi_ColdSmearing(nfout)
    else if(way_of_smearing == FERMI_DIRAC) then
       call m_ESoc_fermi_Dirac(nfout)
    else if(way_of_smearing == MP) then
        call m_ESoc_methfessel_paxton(nfout)
    end if
  end subroutine FermiEnergyLevel

  subroutine FermiEnergyLevel_ek()
    integer :: way_of_smearing
    way_of_smearing = m_CtrlP_way_of_smearing()

    if(way_of_smearing == PARABOLIC) then
       call m_ESoc_fermi_parabolic_ek(nfout)  ! -> efermi, metalic_system
    else if(way_of_smearing == TETRAHEDRON) then
       call m_ESoc_fermi_tetra_ek(nfout)     ! -> efermi, metalic_system
    else if(way_of_smearing == Fermi_Dirac) then
       call m_ESoc_fermi_dirac_ek(nfout)     ! -> efermi, metalic_system
    else if(way_of_smearing == MP) then
       call m_ESoc_mp_ek(nfout)
    end if
  end subroutine FermiEnergyLevel_ek

  subroutine m_SO_wd_WFs()
    integer :: lun
    character*72 file1

#if 1
    if ( nspin == 2 ) call multiply_wfn_by_rotated_spinor
#endif

    lun = 2000
    file1 = "zaj_socsv.data"
    open( lun, file=file1, status="unknown", form="unformatted" )

    if (precision_WFfile==SP) then
       call case_single_precision
    else
       call case_double_precision
    endif
    close(lun)

  contains

    subroutine case_single_precision
      integer :: ik, ib, ri, is1, is2, ierr
      real(kind=SP), allocatable, dimension(:,:)  :: wf_l

      allocate(wf_l(kg1,kimg));    wf_l = 0
      call mpi_barrier(MPI_CommGroup,ierr)

      if ( nspin == 1 ) then
         do ik = 1, kv3
            Do is1=1, 2
               Do is2=1, 2
                  do ib = 1, neg
                     if (map_ek(ib,ik) == mype) then
                        do ri = 1, kimg
                           if (is1==1 .and. is2==1) then
                              wf_l(1:kg1,ri) = zaj_band1_kpt1(1:kg1,map_z(ib),ik,ri)
                           else if (is1==1 .and. is2==2) then
                              wf_l(1:kg1,ri) = zaj_band2_kpt1(1:kg1,map_z(ib),ik,ri)
                           else if (is1==2 .and. is2==1) then
                              wf_l(1:kg1,ri) = zaj_band1_kpt2(1:kg1,map_z(ib),ik,ri)
                           else if (is1==2 .and. is2==2) then
                              wf_l(1:kg1,ri) = zaj_band2_kpt2(1:kg1,map_z(ib),ik,ri)
                           endif
                        end do
                        if (map_ek(ib,ik) /= 0) then
                           call mpi_send( wf_l, kg1*kimg, mpi_real, 0, 1, &
                                &         MPI_CommGroup, ierr)
                        endif
                     else if(mype == 0 .and. map_ek(ib,ik) /= 0) then
                        call mpi_recv( wf_l, kg1*kimg, mpi_real, map_ek(ib,ik), 1, &
                             &         MPI_CommGroup, istatus, ierr )
                     end if
                     if (mype == 0) write(lun)  wf_l
                  end do
               end Do
            end Do
         End Do
      else
         do ik = 1, kv3
            Do is2=1, 2
               do ib = 1, neg
                  if (map_ek(ib,ik) == mype) then
                     do ri = 1, kimg
                        if (is2==1) then
                           wf_l(1:kg1,ri) = zaj_band1_kpt1(1:kg1,map_z(ib),ik,ri)
                        else if (is2==2) then
                           wf_l(1:kg1,ri) = zaj_band2_kpt1(1:kg1,map_z(ib),ik,ri)
                        endif
                     end do
                     if (map_ek(ib,ik) /= 0) then
                        call mpi_send( wf_l, kg1*kimg, mpi_real, 0, 1, &
                             &         MPI_CommGroup, ierr)
                     endif
                  else if(mype == 0 .and. map_ek(ib,ik) /= 0) then
                     call mpi_recv( wf_l, kg1*kimg, mpi_real, map_ek(ib,ik), 1, &
                          &         MPI_CommGroup, istatus, ierr )
                  end if
                  if (mype == 0) write(lun)  wf_l
               end do
            end Do
         End Do
      end if
      deallocate(wf_l)
    end subroutine case_single_precision

    subroutine case_double_precision
      integer :: ik, ib, ri, is1, is2, ierr
      real(kind=DP), allocatable, dimension(:,:)  :: wfdp_l

      allocate(wfdp_l(kg1,kimg));    wfdp_l = 0
      call mpi_barrier(MPI_CommGroup,ierr)

      if ( nspin == 1 ) then
         do ik = 1, kv3
            Do is1=1, 2
               Do is2=1, 2
                  do ib = 1, neg
                     if (map_ek(ib,ik) == mype) then
                        do ri = 1, kimg
                           if (is1==1 .and. is2==1) then
                              wfdp_l(1:kg1,ri) = zaj_band1_kpt1(1:kg1,map_z(ib),ik,ri)
                           else if (is1==1 .and. is2==2) then
                              wfdp_l(1:kg1,ri) = zaj_band2_kpt1(1:kg1,map_z(ib),ik,ri)
                           else if (is1==2 .and. is2==1) then
                              wfdp_l(1:kg1,ri) = zaj_band1_kpt2(1:kg1,map_z(ib),ik,ri)
                           else if (is1==2 .and. is2==2) then
                              wfdp_l(1:kg1,ri) = zaj_band2_kpt2(1:kg1,map_z(ib),ik,ri)
                           endif
                        end do
                        if (map_ek(ib,ik) /= 0) then
                           call mpi_send( wfdp_l, kg1*kimg, mpi_double_precision, &
                                &         0, 1, &
                                &         MPI_CommGroup, ierr)
                        endif
                     else if(mype == 0 .and. map_ek(ib,ik) /= 0) then
                        call mpi_recv( wfdp_l, kg1*kimg, mpi_double_precision, &
                             &         map_ek(ib,ik), 1, &
                             &         MPI_CommGroup, istatus, ierr )
                     end if
                     if (mype == 0) write(lun)  wfdp_l
                  end do
               end Do
            end Do
         End Do
      else
         do ik = 1, kv3
            Do is2=1, 2
               do ib = 1, neg
                  if (map_ek(ib,ik) == mype) then
                     do ri = 1, kimg
                        if (is2==1) then
                           wfdp_l(1:kg1,ri) = zaj_band1_kpt1(1:kg1,map_z(ib),ik,ri)
                        else if (is2==2) then
                           wfdp_l(1:kg1,ri) = zaj_band2_kpt1(1:kg1,map_z(ib),ik,ri)
                        endif
                     end do
                     if (map_ek(ib,ik) /= 0) then
                        call mpi_send( wfdp_l, kg1*kimg, mpi_real, 0, 1, &
                             &         MPI_CommGroup, ierr)
                     endif
                  else if(mype == 0 .and. map_ek(ib,ik) /= 0) then
                     call mpi_recv( wfdp_l, kg1*kimg, mpi_real, map_ek(ib,ik), 1, &
                          &         MPI_CommGroup, istatus, ierr )
                  end if
                  if (mype == 0) write(lun)  wfdp_l
               end do
            end Do
         End Do
      end if
      deallocate(wfdp_l)
    end subroutine case_double_precision

  end subroutine m_SO_wd_WFs

  subroutine multiply_wfn_by_rotated_spinor
    use m_Crystal_Structure,  only : sw_fix_global_quantz_axis, &
         &                           Global_Quantz_Axis_Fixed
    use m_Const_Parameters,   only : zi

    integer :: ik, ib, ig
    real(kind=DP) :: spn_quant_dir(3), theta, phi
    complex(kind=CMPLDP) :: rotated_spinor(2,2), z1, z2, ztmp1, ztmp2

    if ( sw_fix_global_quantz_axis == OFF ) return

    spn_quant_dir(1:3) = Global_Quantz_Axis_Fixed(1:3)
    theta = acos( spn_quant_dir(3) )
    phi = atan2( spn_quant_dir(2), spn_quant_dir(1) )

    if ( theta == 0.0D0 ) return

    rotated_spinor(1,1) =  exp( -zi *phi /2.0d0 ) *cos( theta /2.0d0 )
    rotated_spinor(2,1) =  exp(  zi *phi /2.0d0 ) *sin( theta /2.0d0 )

    rotated_spinor(1,2) = -exp( -zi *phi /2.0d0 ) *sin( theta /2.0d0 )
    rotated_spinor(2,2) =  exp(  zi *phi /2.0d0 ) *cos( theta /2.0d0 )

    do ik = 1, kv3, nspin
       if ( map_k(ik) /= myrank_k ) cycle
       Do ig=1, kg1
          Do ib=1, np_e
             z1 = dcmplx( zaj_band1_kpt1(ig,ib,ik,  1), zaj_band1_kpt1(ig,ib,ik,  kimg) )
             z2 = dcmplx( zaj_band1_kpt1(ig,ib,ik+1,1), zaj_band1_kpt1(ig,ib,ik+1,kimg) )
             ztmp1 = rotated_spinor(1,1) *z1 +rotated_spinor(2,1) *z2
             ztmp2 = rotated_spinor(1,2) *z1 +rotated_spinor(2,2) *z2
             zaj_band1_kpt1(ig,ib,ik,1)      = real(ztmp1)
             zaj_band1_kpt1(ig,ib,ik,kimg)   = aimag(ztmp1)
             zaj_band1_kpt1(ig,ib,ik+1,1)    = real(ztmp2)
             zaj_band1_kpt1(ig,ib,ik+1,kimg) = aimag(ztmp2)

             z1 = dcmplx( zaj_band2_kpt1(ig,ib,ik,  1), zaj_band2_kpt1(ig,ib,ik,  kimg) )
             z2 = dcmplx( zaj_band2_kpt1(ig,ib,ik+1,1), zaj_band2_kpt1(ig,ib,ik+1,kimg) )
             ztmp1 = rotated_spinor(1,1) *z1 +rotated_spinor(2,1) *z2
             ztmp2 = rotated_spinor(1,2) *z1 +rotated_spinor(2,2) *z2
             zaj_band2_kpt1(ig,ib,ik,1)      = real(ztmp1)
             zaj_band2_kpt1(ig,ib,ik,kimg)   = aimag(ztmp1)
             zaj_band2_kpt1(ig,ib,ik+1,1)    = real(ztmp2)
             zaj_band2_kpt1(ig,ib,ik+1,kimg) = aimag(ztmp2)
          End Do
       End Do
    end do
  end subroutine multiply_wfn_by_rotated_spinor

  subroutine m_SO_wd_EigenValues( fin_mode, no_occ )
    logical, intent(in) :: fin_mode, no_occ
    integer                          :: ie,  ipri0, kv3_i, ks
    integer                          :: hconst_min, lzero_max
    integer, parameter :: NCOLUMN = 6
    integer, parameter :: EIGEN_VALUES = 1, OCCUPATIONS = 2
    real(kind=DP),allocatable, dimension(:,:) :: e_mpi, o_mpi

    integer :: lun
    character*72 :: file1

    lun = 1000
    file1 = "nfefermi_socsv.data"

    if ( mype == 0 ) then
       if ( fin_mode ) then
          open( lun, file=file1, status="unknown", form="formatted" )
          write( lun,'(f16.8," : Efermi")') efermi
          close( lun )
       endif
    endif

    file1 = "nfenergy_socsv.data"
    if ( mype == 0 ) then
       if ( fin_mode ) then
          open( lun, file=file1, status="unknown", form="formatted" )
       else
          if ( nk_in_the_process == 1 ) then
             open( lun, file=file1, status="unknown", form="formatted" )
          else
             open( lun, file=file1, status="unknown", form="formatted", &
                  &     position='append' )
          endif
       endif
    endif
    allocate(e_mpi(neg_doubled,kv3/nspin)); e_mpi = 0.d0
    allocate(o_mpi(neg_doubled,kv3/nspin)); o_mpi = 0.d0

    call set_kv3_i_and_ks() ! -> kv3_i, ks

    call put_kpartArray_into(eko_l,e_mpi)
    if ( mype==0 ) call wd_k_and_values(EIGEN_VALUES)

    if ( .not. no_occ ) then
       call put_kpartArray_into(occup_l,o_mpi)
       if ( mype==0 ) call wd_k_and_values(OCCUPATIONS)
    endif
    deallocate(e_mpi);   deallocate(o_mpi)

    if ( mype == 0 ) close(lun)

  contains

    subroutine set_kv3_i_and_ks()
      if((icond == FIXED_CHARGE_CONTINUATION .or. icond == FIXED_CHARGE) .and. &
           & fixed_charge_k_parallel == ONE_BY_ONE) then
         kv3_i = kv3_ek - kv3*(nkgroup-1)
         if(kv3_i > kv3) kv3_i = kv3
         ks = max(1,first_kpoint_in_this_job) - 1 + kv3*(nkgroup-1)
      else
         kv3_i = kv3
         ks = 0
      end if

    end subroutine set_kv3_i_and_ks

    subroutine wd_k_and_values(mode)
      integer, intent(in) :: mode
      integer :: ik, nb, iktmp
      integer :: ie_s, ie_e, nhw, neg_t
      real(kind=DP) :: hw, hc, hv

      do ik = 1, kv3_i /nspin
         iktmp = (ik-1) *nspin +1
         if ( mode == EIGEN_VALUES) write(lun,'(" ===== energy eigen values =====")')
         if ( mode == OCCUPATIONS)  write(lun,'(" ===== occupations =====")')
         call wd_k_points_socsv(ik)
         neg_t = neg_doubled
         if(neg_is_enlarged) neg_t = neg_doubled - num_extra_bands*2

         if(mode == EIGEN_VALUES) then
            write(lun,'(5f16.8)') (e_mpi(nb,ik),nb=1,neg_t)
                                                    ! =eko(neordr(nb,ik),ik)
         else if(mode == OCCUPATIONS) then
            write(lun,'(5f16.8)') (o_mpi(nb,ik)/(qwgt(iktmp)*kv3),nb = 1,neg_t)
                                                    ! =occup(neordr(nb,ik),ik)
         end if
      end do
    end subroutine wd_k_and_values

    subroutine put_kpartArray_into(a_l,a_all)
      real(kind=DP), intent(in), dimension(np_e,ista_k:iend_k) :: a_l
      real(kind=DP), intent(out), dimension(neg_doubled,kv3/nspin) :: a_all
      integer :: ik, iktmp, ierr, ie, is, ito

      a_all = 0.d0
      if ( nspin==1 ) then
         do ik = 1, kv3
            if(map_k(ik) /= myrank_k) cycle
            do ie = 1, neg
               if(map_e(ie) /= myrank_e) cycle
               ito = nrvf_ordr(ie,ik)
               Do is=1, ndim_spinor_socsv
                  a_all( (ito-1)*ndim_spinor_socsv +is,ik) = a_l(map_z(ie),ik)
               End Do
            end do
         end do
      else
         do ik = 1, kv3, nspin
            iktmp = ( ik-1 )/nspin +1
            if(map_k(ik) /= myrank_k) cycle
            do ie = 1, neg
               if(map_e(ie) /= myrank_e) cycle
               ito = nrvf_ordr(ie,ik)
               Do is=1, ndim_spinor_socsv
                  a_all( (ito-1)*ndim_spinor_socsv +is,iktmp) = a_l(map_z(ie),ik+is-1)
               End Do
            end do
         end do
      endif

      if(npes >= 2) then
         call mpi_allreduce( MPI_IN_PLACE, a_all, neg_doubled*kv3/nspin,&
              &             mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
      end if
    end subroutine put_kpartArray_into

    subroutine wd_k_points_socsv(ik)
      integer, intent(in) :: ik
      real(kind=DP) :: vkxyz_wk(3)

      vkxyz_wk(1:3) = vkxyz(ik,1:3,BUCS)
      write(lun,'(" ik = ",i5," (",3f10.6," )")') &
           &                 ndim_spinor_socsv*(ik+ks-1)+1, vkxyz_wk(1:3)
    end subroutine wd_k_points_socsv

  end subroutine m_SO_wd_EigenValues

  subroutine m_SO_wd_EigenValues_ek_fin
    real(kind=DP), parameter :: delta = 1.d-12
    real(kind=DP), allocatable, dimension(:) :: eko_t
    integer, allocatable, dimension(:)       :: neordr_t
    integer                     :: ik, ib,jb,ibo,jbo, neg_t

    integer :: lun
    character*72 :: file1
    real(kind=DP) :: efermi_save, vbm_save
    real(kind=DP), allocatable :: eko_ek_save(:,:)

    lun = 1000
    file1 = "nfenergy_socsv.data"
    if ( mype == 0 ) then
       open( lun, file=file1, status="unknown", form="formatted" )
    endif

    vbm_save = vbm;    efermi_save = efermi
    allocate( eko_ek_save(neg,kv3_ek) );  eko_ek_save = eko_ek
    allocate(eko_t(neg_doubled));     eko_t = 0.0d0
    allocate(neordr_t(neg_doubled));  neordr_t = 0

    eko_ek = eko_ek_socsv

    call FermiEnergylevel_ek()

    if ( mype == 0 ) then
       write(lun,'(" num_kpoints = ",i6)') kv3_ek /nspin
       if(neg_is_enlarged) then
          write(lun,'(" num_bands   = ",i6)') neg_doubled -num_extra_bands*2
       else
          write(lun,'(" num_bands   = ",i6)') neg_doubled
       end if
       write(lun,'(" nspin       = ",i6)') 1

       if(metalic_system) then
          write(lun,'(" Fermi energy level = ",f10.6/)') efermi
       else
          write(lun,'(" Valence band max   = ",f10.6/)') vbm
       end if

       write(lun,'(" nk_converged = ",i8)') min(kv3_ek,nk_converged)
       do ik = 1, kv3_ek/nspin
          call wd_k_points_socsv
       end do
       write(lun,'(" -----")')
    end if

    do ik = 1, nk_converged, nspin
       if (ik > kv3_ek) cycle
       if (mype == 0) write(lun,'("=== energy_eigen_values ===")')

       if ( nspin == 1 ) then
          Do ib=1, neg
             eko_t(2*ib-1) = eko_ek_socsv(ib,ik)
             eko_t(2*ib  ) = eko_ek_socsv(ib,ik)
          End do
       else
          Do ib=1, neg
             eko_t(2*ib-1) = eko_ek_socsv(ib,ik)
             eko_t(2*ib  ) = eko_ek_socsv(ib,ik+1)
          End do
       endif

       neordr_t(1:neg_doubled) = (/(ib,ib=1,neg_doubled)/)
       do ib = 1, neg_doubled -1
          do jb = ib+1, neg_doubled
             ibo = neordr_t(ib);  jbo = neordr_t(jb)
             if(eko_t(jbo)  < eko_t(ibo)-delta) then        ! MPI
                neordr_t(jb) = ibo;   neordr_t(ib) = jbo
             end if
          end do
       end do
       if ( mype==0 ) then
          call wd_k_points_socsv
          neg_t = neg_doubled
          if(neg_is_enlarged) neg_t = neg_doubled - num_extra_bands *2
          write(lun,'(4f18.10)') (eko_t(neordr_t(ib)),ib=1,neg_t)
       end if
    end do

    deallocate(neordr_t);     deallocate(eko_t)

    eko_ek = eko_ek_save;  vbm = vbm_save;  efermi = efermi_save
    deallocate( eko_ek_save )
    if ( mype==0 ) close(lun)

  contains

    subroutine wd_k_points_socsv
      integer :: j, k
      real(kind=DP) :: c1, vkxyz_wk(3)

      vkxyz_wk(1:3) = vkxyz_ek(ik,1:3,BUCS)
      write(lun,'(" ik = ",i8," (",3f10.6," )")') &
           &                   ndim_spinor_socsv*(ik-1)+1, vkxyz_wk(1:3)

    end subroutine wd_k_points_socsv

  end subroutine m_SO_wd_EigenValues_ek_fin

  subroutine m_SO_cp_eko_l_to_eko_ek_socsv()
    real(kind=DP),allocatable,dimension(:,:) :: eko_t, eko_t2 ! d(neg,kv3)
    integer :: ik, ib, kv3_e
    integer :: ikt,is, kv3t

    allocate(eko_t(neg,kv3)); allocate(eko_t2(neg,kv3))

    eko_t = 0.d0
    do ik = 1, kv3
       if(map_k(ik) /= myrank_k) cycle
       do ib = 1, neg
          if(map_e(ib) /= myrank_e) cycle
          eko_t(ib,ik) = eko_l(map_z(ib),ik)
       end do
    end do
    if(npes >= 2) then
       call mpi_allreduce(eko_t,eko_t2,neg*kv3,mpi_double_precision &
            &              ,mpi_sum, MPI_CommGroup,ierr)
    else
       eko_t2 = eko_t
    end if

    kv3_e = min(kv3,kv3_ek-kv3*(nkgroup-1))
    if(nk_in_the_process+kv3_e-1 > kv3_ek) kv3_e = kv3_ek - nk_in_the_process + 1
    if(kv3_e <=0) kv3_e = 1

    if(sw_modified_kpoint_increment == ON)then
      do ik=0,nrank_k-1
         do is=1,nspin
            ikt = nspin*(nis_kv3_ek(ik)-1)+(nkgroup-1)*nspin+is
            if(ikt>kv3_ek) ikt = kv3_ek
            eko_ek_socsv(:,ikt) = eko_t2(:,ik*nspin+is)
         enddo
      enddo
    else
      eko_ek_socsv(:,nk_in_the_process:nk_in_the_process+kv3_e-1) = eko_t2(:,1:kv3_e)
    endif

    deallocate(eko_t); deallocate(eko_t2)
  end subroutine m_SO_cp_eko_l_to_eko_ek_socsv

  subroutine m_SO_calc_band_energy_socsv_ek
    use m_ES_occup_EPS, only : occup_l_ek

    integer :: ik, ib
    real(kind=DP) :: eband, eband_mpi

    allocate( occup_l_ek( np_e, kv3_ek ) )
    eko_ek = eko_ek_socsv

    call FermiEnergyLevel_ek2

    eband = 0.d0
    do ik =1, kv3_ek
       do ib = ista_e, iend_e, istep_e
          eband = eband + occup_l_ek( map_z(ib),ik ) *eko_ek(ib,ik)
       end do
    end do

    if (npes > 1) then
       call mpi_allreduce( eband, eband_mpi, 1, mpi_double_precision, &
            &              mpi_sum, mpi_k_world(myrank_k),ierr )
       eband = eband_mpi
    end if

    eband = 2.0d0 *eband /dble(kv3_ek)

!    write(nfout,*)
    write(nfout,*) "<<< Spin-orbit by Second Variation (ekmode)"
    write(nfout,*) "* band energy  = ", eband
   write(nfout,*) ">>>"
    write(nfout,*)

    deallocate( occup_l_ek )

  end subroutine m_SO_calc_band_energy_socsv_ek

  subroutine FermiEnergyLevel_ek2()
    use m_Control_Parameters,   only : iprioccup
    use m_ES_occup_EPS, only : m_ESoc_EPS_fermi_parabolic_ek, &
         &                     m_ESoc_EPS_fermi_tetra_ek, &
         &                     m_ESoc_EPS_fermi_dirac_ek

    integer :: way_of_smearing, iprioccup_save

    efermi = 0.0d0
    iprioccup_save = iprioccup;  iprioccup = -1

    way_of_smearing = m_CtrlP_way_of_smearing()

    if(way_of_smearing == PARABOLIC) then
       call m_ESoc_EPS_fermi_parabolic_ek(nfout)
    else if(way_of_smearing == TETRAHEDRON) then
       call m_ESoc_EPS_fermi_tetra_ek(nfout)
    else if ( way_of_smearing == Fermi_Dirac ) then
       call m_ESoc_EPS_fermi_dirac_ek(nfout)
    end if

    iprioccup = iprioccup_save

  end subroutine FermiEnergyLevel_ek2

end module m_SpinOrbit_SecondVariation
