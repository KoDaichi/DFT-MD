module m_CLS_dipquad
  use m_Const_Parameters,  only : DP, CMPLDP, zi, ON, BUCS, PAI2, PAI, PAI4, &
       &                          MacroDielectric, PACS, Hartree, InvHyperFineConst
  use m_PseudoPotential,    only : nlmt, nloc, ntau, radr_paw, iltpw, lppw, tppw, &
       &                           nmesh, mmesh, ltp, mtp, taup, phirpw, ilmt, psirpw, &
       &                           nlmt_add, ilmt_add, ltp_add, mtp_add, lpsmax, itau, &
       &                           lmta, lmta_add, nlmta
  use m_Ionic_System,   only :  ityp
  use m_CoreLevel_Spectrum,  only : qnum_l_to_probe, psir_core_ae_wfns, &
       &                            m_CLS_find_orb_index_to_probe, atom_to_probe, &
       &                            num_core_states, fsr_core_states, fsi_core_states, &
       &                            ene_core_states, mimic_soc_split_spectrum, &
       &                            ene_initial_state_splitting, qnum_n_to_probe, &
       &                            m_CLS_set_wfn_core_states, &
       &                            m_CLS_alloc_wfn_core_states

  use m_SpinOrbit_Potential,   only : MatU_ylm_RC_L0, MatU_ylm_RC_L1, &
       &                              MatU_ylm_RC_L2, MatU_ylm_RC_L3, &
       &                              MatU_ylm_RC_L4

  use m_Parallelization,    only : mype, np_e, myrank_k, map_k, map_z, &
       &                           ista_e, iend_e, istep_e, npes, MPI_CommGroup, &
       &                           map_e, myrank_e, ista_fs, iend_fs, mpi_ke_world, &
       &                           nrank_g, mpi_g_world

  use m_Control_Parameters, only : sw_use_add_proj, ndim_spinor, nspin, neg, noncol, &
       &                           sw_v2c_xes
  use m_Kpoints,            only : kv3_ek, kv3, qwgt_ek, num_star_of_k, star_of_k, &
       &                           iopr_k_fbz_to_ibz, trev_k_fbz_to_ibz, kv3_fbz
  use m_IterationNumbers,           only : nk_in_the_process
  use m_Electronic_Structure,  only : fsr_l, fsi_l, fsr_add_l, fsi_add_l, eko_ek, &
       &                              efermi
  use m_CS_Magnetic,         only : magmom_dir_inversion_opr_flag
  use m_Crystal_Structure,   only : op, tau, univol
  use m_Files,               only : nfepsout, m_Files_open_nfeps, m_Files_close_nfeps, &
       &                            nfout
  use m_ES_nonlocal,    only : m_ES_add_betar_dot_WFs
  use mpi

  implicit none
!  include 'mpif.h'

  integer, private :: lmax = 4
  integer, private, parameter :: nsph_max = 25
  integer, private, parameter :: PARABOLIC_B = 1, L_TETRAHEDRON = 2, &
       &                         GAUSSIAN_B = 3, LORENTZIAN_B = 4

  real(kind=DP), parameter :: hartree_in_eV = Hartree
! real(kind=DP), parameter :: InvHyperFineConst = 137.035999679D0

  complex(kind=CMPLDP) :: rhat_cylm(nsph_max,nsph_max,3)
  real(kind=DP) :: rhat_rylm(nsph_max,nsph_max,3)
  real(kind=DP) :: rhat2_rylm(nsph_max,nsph_max,3,3)

  real(kind=DP), allocatable :: rdipole_val_core(:,:,:)
  real(kind=DP), allocatable :: rquadrupole_val_core(:,:,:,:)

  complex(kind=CMPLDP), allocatable :: trm_dipole(:,:,:,:)
  complex(kind=CMPLDP), allocatable :: trm_quadrupole(:,:,:,:,:)

contains

! fixed_charge version
  subroutine m_CLS_alloc_transmom_ek
    allocate( trm_dipole( np_e, num_core_states, kv3_ek, 3 ) )
    allocate( trm_quadrupole( np_e, num_core_states, kv3_ek, 3, 3 ) )
    trm_dipole = 0.0d0
    trm_quadrupole = 0.0d0
  end subroutine m_CLS_alloc_transmom_ek

  subroutine m_CLS_dealloc_transmom_ek
    deallocate( trm_dipole )
    deallocate( trm_quadrupole )
  end subroutine m_CLS_dealloc_transmom_ek

  subroutine m_CLS_wd_transmom_ek
    integer :: lun

    if ( mype == 0 ) then
       lun = get_unused_unitnumber()
       open( lun, file="trans_bin.data", form="unformatted", status="unknown")
!       open( newunit=lun, file="trans_bin.data", form="unformatted", status="unknown")
    endif
    call write_trm_dipole(lun)
    call write_trm_quadrupole(lun)
    if ( mype == 0 )  close( lun )

  contains

    subroutine write_trm_dipole(lun)
      integer, intent(in) :: lun
      integer :: ik, ib1, ib2, nsize, ierr
      complex(kind=CMPLDP), allocatable :: zwk1(:,:,:,:)
      complex(kind=CMPLDP), allocatable :: zwk1_mpi(:,:,:,:)

      nsize = num_core_states *neg *kv3_ek *3
      allocate( zwk1(neg,num_core_states,kv3_ek,3) )
      zwk1 = 0.0d0

      Do ib2=1, num_core_states
         Do ik=1, kv3_ek
            Do ib1=1, neg
               if ( map_e(ib1) /= myrank_e ) cycle
               zwk1(ib1,ib2,ik,1:3) = trm_dipole(map_z(ib1),ib2,ik,1:3)
            End do
         End do
      End Do
      if ( npes > 1 ) then
         allocate( zwk1_mpi(neg,num_core_states,kv3_ek,3) )
         call mpi_allreduce( zwk1, zwk1_mpi, nsize*2, mpi_double_precision, &
              &              mpi_sum, mpi_g_world, ierr )
         zwk1 = zwk1_mpi
         deallocate( zwk1_mpi )
      endif
      if ( mype == 0 ) then
         write(lun) 2
         write(lun) num_core_states, neg, kv3_ek
         write(lun) zwk1
      endif
      deallocate( zwk1 )
    end subroutine write_trm_dipole

    subroutine write_trm_quadrupole(lun)
      integer, intent(in) :: lun

      integer :: ik, ib1, ib2, nsize, ierr
      complex(kind=CMPLDP), allocatable :: zwk1(:,:,:,:,:)
      complex(kind=CMPLDP), allocatable :: zwk1_mpi(:,:,:,:,:)

      nsize = num_core_states *neg *kv3_ek *3 *3
      allocate( zwk1(neg,num_core_states,kv3_ek,3,3) )
      zwk1 = 0.0d0

      Do ib2=1, num_core_states
         Do ik=1, kv3_ek
            Do ib1=1, neg
               if ( map_e(ib1) /= myrank_e ) cycle
               zwk1(ib1,ib2,ik,1:3,1:3) = trm_quadrupole(map_z(ib1),ib2,ik,1:3,1:3)
            End do
         End do
      End Do
      if ( npes > 1 ) then
         allocate( zwk1_mpi(neg,num_core_states,kv3_ek,3,3) )
         call mpi_allreduce( zwk1, zwk1_mpi, nsize*2, mpi_double_precision, &
              &              mpi_sum, mpi_g_world, ierr )
         zwk1 = zwk1_mpi
         deallocate( zwk1_mpi )
      endif
      if ( mype == 0 ) then
         write(lun) 4
         write(lun) num_core_states, neg, kv3_ek
         write(lun) zwk1
      endif
      deallocate( zwk1 )
    end subroutine write_trm_quadrupole

  end subroutine m_CLS_wd_transmom_ek

  subroutine m_CLS_rd_transmom_ek
    integer :: lun

    if ( mype == 0 ) then
       lun = get_unused_unitnumber()
       open( lun, file="trans_bin.data", form="unformatted", status="old" )
!       open( newunit=lun, file="trans_bin.data", form="unformatted", status="old" )
    endif
    call read_trm_dipole(lun)
    call read_trm_quadrupole(lun)
    if ( mype == 0 )  close( lun )

  contains

    subroutine read_trm_dipole(lun)
      integer, intent(in) :: lun
      integer :: ik, ib1, ib2, nsize, n1, n2, n3, n4, ierr
      complex(kind=CMPLDP), allocatable :: zwk1(:,:,:,:)
      complex(kind=CMPLDP), allocatable :: zwk1_mpi(:,:,:,:)

      nsize = num_core_states *neg *kv3_ek *3
      allocate( zwk1(neg,num_core_states,kv3_ek,3) )

      if ( mype == 0 ) then
         read(lun) n1
         read(lun) n2, n3, n4
         if ( n1 /= 2 ) then
            call phase_error_with_msg(nfout,"trans bin data (dipole) is different",__LINE__,__FILE__)
         endif
         if ( n2 /= num_core_states ) then
            call phase_error_with_msg(nfout,"trans bin data (dipole) is different",__LINE__,__FILE__)
         endif
         if ( n3 /= neg ) then
            call phase_error_with_msg(nfout,"trans bin data (dipole) is different",__LINE__,__FILE__)
         endif
         if ( n4 /= kv3_ek ) then
            call phase_error_with_msg(nfout,"trans bin data (dipole) is different",__LINE__,__FILE__)
         endif
         read(lun) zwk1
      endif
      call mpi_bcast( zwk1, nsize*2, mpi_double_precision, 0, MPI_CommGroup, ierr )

      Do ib2=1, num_core_states
         Do ik=1, kv3_ek
            Do ib1=1, neg
               if ( map_e(ib1) /= myrank_e ) cycle
               trm_dipole(map_z(ib1),ib2,ik,1:3) = zwk1(ib1,ib2,ik,1:3)
            End do
         End do
      End Do
      deallocate( zwk1 )
    end subroutine read_trm_dipole

    subroutine read_trm_quadrupole(lun)
      integer, intent(in) :: lun

      integer :: ik, ib1, ib2, nsize, n1, n2, n3, n4, ierr
      complex(kind=CMPLDP), allocatable :: zwk1(:,:,:,:,:)
      complex(kind=CMPLDP), allocatable :: zwk1_mpi(:,:,:,:,:)

      nsize = num_core_states *neg *kv3_ek *3 *3
      allocate( zwk1(neg,num_core_states,kv3_ek,3,3) )

      if ( mype == 0 ) then
         read(lun) n1
         read(lun) n2, n3, n4
         if ( n1 /= 4 ) then
            call phase_error_with_msg(nfout,"trans bin data (quadrupole) is different",__LINE__,__FILE__)
         endif
         if ( n2 /= num_core_states ) then
            call phase_error_with_msg(nfout,"trans bin data (quadrupole) is different",__LINE__,__FILE__)
         endif
         if ( n3 /= neg ) then
            call phase_error_with_msg(nfout,"trans bin data (quadrupole) is different",__LINE__,__FILE__)
         endif
         if ( n4 /= kv3_ek ) then
            call phase_error_with_msg(nfout,"trans bin data (quadrupole) is different",__LINE__,__FILE__)
         endif
         read(lun) zwk1
      endif
      call mpi_bcast( zwk1, nsize*2, mpi_double_precision, 0, MPI_CommGroup, ierr )

      Do ib2=1, num_core_states
         Do ik=1, kv3_ek
            Do ib1=1, neg
               if ( map_e(ib1) /= myrank_e ) cycle
               trm_quadrupole(map_z(ib1),ib2,ik,1:3,1:3) = zwk1(ib1,ib2,ik,1:3,1:3)
            End do
         End do
      End Do
      deallocate( zwk1 )
    end subroutine read_trm_quadrupole

  end subroutine m_CLS_rd_transmom_ek

  subroutine m_CLS_calc_transmom_ek
    integer :: nlmt_core
    integer :: ia, it, ik, ik2, ib1, ib2
    integer :: lmt1, lmt2, lmta2
    complex(kind=CMPLDP) :: z1, wf1, wf2
    complex(kind=CMPLDP), allocatable :: wk_fsri(:)

    integer :: i, iadd, ierr
    real(kind=DP), allocatable :: fs(:,:,:)
    real(kind=DP), allocatable :: work(:,:,:)

    nlmt_core = 2*qnum_l_to_probe+1
    allocate( wk_fsri(nlmt_core) ); wk_fsri = 0.0d0

    allocate(fs(np_e,nlmta,2));  fs = 0.0d0

    call m_CLS_alloc_wfn_core_states( .false. )
    call m_CLS_set_wfn_core_states( .false. )

    if ( sw_use_add_proj == ON ) call m_ES_add_betar_dot_WFs(nfout)

    Do ia = atom_to_probe, atom_to_probe
       it = ityp(ia)

       Do ik=1, kv3
          if ( map_k(ik) /= myrank_k ) cycle
          ik2 = nk_in_the_process +ik -1

          fs = 0.0d0
          DO ib2=1, np_e
             Do i=ista_fs, iend_fs
                iadd = i - ista_fs + 1
                fs(ib2,i,1) = fsr_l(ib2,iadd,ik)
                fs(ib2,i,2) = fsi_l(ib2,iadd,ik)
             end do
          End DO
          call mpi_allreduce( MPI_IN_PLACE, fs, np_e*nlmta*2, mpi_double_precision, &
               &              mpi_sum, mpi_ke_world, ierr )

          Do ib1=1, num_core_states
             wk_fsri(:) = dcmplx( fsr_core_states( ib1,:,ik ), &
                  &               fsi_core_states( ib1,:,ik ) )

             DO ib2=1, np_e

                Do lmt2=1, ilmt(it)
                   lmta2 = lmta(lmt2,ia)
                   wf2 = dcmplx( fs(ib2,lmta2,1), fs(ib2,lmta2,2) )

                   Do lmt1=1, nlmt_core
                      wf1 = wk_fsri(lmt1)
                      z1 =  conjg(wf2) *wf1
                      trm_dipole( ib2, ib1, ik2, 1:3 ) &
                           &  = trm_dipole( ib2, ib1, ik2, 1:3 ) &
                           &    + z1 *rdipole_val_core( lmt2, lmt1, 1:3 )
                      trm_quadrupole( ib2, ib1, ik2, 1:3, 1:3 ) &
                           &  = trm_quadrupole( ib2, ib1, ik2, 1:3, 1:3 ) &
                           &    +z1 *rquadrupole_val_core( lmt2, lmt1, 1:3, 1:3 )
                   End do
                End Do
!
                if ( sw_use_add_proj == ON ) then
                   Do lmt2=1, ilmt_add(it)
                      lmta2 = lmta_add(lmt2,ia)
                      wf2 = dcmplx( fsr_add_l( ib2, lmta2, ik ), &
                           &        fsi_add_l( ib2, lmta2, ik ) )

                      Do lmt1=1, nlmt_core
                         wf1 = wk_fsri(lmt1)
                         z1 = conjg(wf2) *wf1
                         trm_dipole( ib2, ib1, ik2, 1:3 ) &
                              &  = trm_dipole( ib2, ib1, ik2, 1:3 ) &
                              &    + z1 *rdipole_val_core( lmt2+nlmt, lmt1, 1:3 )
                         trm_quadrupole( ib2, ib1, ik2, 1:3, 1:3 )&
                              &  = trm_quadrupole( ib2, ib1, ik2, 1:3, 1:3 ) &
                              &    + z1 *rquadrupole_val_core( lmt2+nlmt, lmt1, 1:3, 1:3 )
                      End do
                   End Do
                endif

             End Do
          End DO
       End Do
    End Do
    deallocate( wk_fsri )
    deallocate( fs )
  end subroutine m_CLS_calc_transmom_ek

  subroutine m_CLS_calc_spectr_fn_l( nstep, evec, kvec, scissor, way_BZ_int, &
       &                             occup_l_ek, ene, e_step, width, spectrum_type )
    integer, intent(in) :: nstep, way_BZ_int, spectrum_type
    real(kind=DP), intent(in) :: evec(3), kvec(3), scissor, ene(nstep), e_step, width
    real(kind=DP), intent(in) :: occup_l_ek(neg,kv3_ek)

    integer :: i, ik, jk, iopr, trev
    integer :: ib1, ib2, is, ixyz1, ixyz2, ierr
    real(kind=DP) :: evec_rotated(3), kvec_rotated(3), weight
    real(kind=DP) :: ediff, odiff, occ1, occ2, txyz(3)
    real(kind=DP) :: Delta03 = 1.0D-3

    complex(kind=CMPLDP) :: zrho1, zrho2
    real(kind=DP) :: trm2_E1_E1, trm2_E2_E2, trm2_E1_E2
    real(kind=DP), allocatable :: work(:)
    real(kind=DP), allocatable :: spectr_E1_E1(:), spectr_E2_E2(:), spectr_E1_E2(:)

    real*8 :: fac1, fac2, wavenum_xray

    allocate( spectr_E1_E1(nstep) ); spectr_E1_E1 = 0.0d0
    allocate( spectr_E2_E2(nstep) ); spectr_E2_E2 = 0.0d0
    allocate( spectr_E1_E2(nstep) ); spectr_E1_E2 = 0.0d0

    Do ib1=1, num_core_states

       Do ik=1, kv3_ek, ndim_spinor
          weight = kv3_ek *qwgt_ek(ik) /dble(ndim_spinor)

          Do i=1, num_star_of_k(ik)
             jk = star_of_k(ik,i)
             iopr = iopr_k_fbz_to_ibz(jk)
             trev = trev_k_fbz_to_ibz(jk)

             evec_rotated = matmul( op(:,:,iopr), evec(:) )
             kvec_rotated = matmul( op(:,:,iopr), kvec(:) )
             if ( trev == 1 ) then
                evec_rotated = -evec_rotated
                kvec_rotated = -kvec_rotated
             endif

             txyz(1:3) = tau(1:3,iopr,BUCS)*PAI2

             Do ib2=1, neg
                if ( map_e(ib2) /= myrank_e ) cycle
                occ2 = occup_l_ek(map_z(ib2),ik) /weight

                if ( sw_v2c_xes == ON ) then
                   if ( occ2 < Delta03 ) cycle
                   odiff = occ2
                else
                   if ( 1.0D0 -occ2 < Delta03 ) cycle
                   odiff = 1.0d0 -occ2
                endif

                ediff = eko_ek(ib2,ik) +scissor -ene_core_states(ib1)
                wavenum_xray = ediff /InvHyperFineConst
                fac1 = wavenum_xray /2.0d0
                fac2 = fac1 *fac1

                zrho1 = 0.0d0;      zrho2 = 0.0d0
                Do is=1, ndim_spinor
                   Do ixyz1=1, 3
                      zrho1 = zrho1 &
                           & + trm_dipole( map_z(ib2), ib1, ik+is-1, ixyz1 ) &
                           &   *evec_rotated(ixyz1)
                   End Do
                   Do ixyz1=1, 3
                      Do ixyz2=1, 3
                         zrho2 = zrho2 &
                              & + trm_quadrupole( map_z(ib2), ib1, ik+is-1,ixyz1,ixyz2 ) &
                              &   *evec_rotated(ixyz1) *kvec_rotated(ixyz2)
                      End Do
                   ENd Do
                End do

                zrho2 = zrho2 *zi

                trm2_E1_E1 = conjg( zrho1 ) *zrho1
                trm2_E2_E2 = conjg( zrho2 ) *zrho2 *fac2
                trm2_E1_E2 = ( conjg( zrho1  ) *zrho2 &
                     &         +conjg( zrho2  ) *zrho1 ) *fac1


                call add_spectral_function( ediff, odiff, trm2_E1_E1, &
                     &                      way_BZ_int, ene, e_step, width, &
                     &                      nstep, Spectr_E1_E1 )
                call add_spectral_function( ediff, odiff, trm2_E2_E2, &
                     &                      way_BZ_int, ene, e_step, width, &
                     &                      nstep, Spectr_E2_E2 )
                call add_spectral_function( ediff, odiff, trm2_E1_E2, &
                     &                      way_BZ_int, ene, e_step, width, &
                     &                      nstep, Spectr_E1_E2 )
             End do
          End Do
       End do
    End Do

    if ( npes > 1 ) then
       allocate( work(nstep ) ); work = 0.0d0
       call mpi_allreduce( Spectr_E1_E1, work, nstep, &
            &              mpi_double_precision, mpi_sum, mpi_g_world, ierr )
       Spectr_E1_E1 = work

       work = 0.0d0
       call mpi_allreduce( Spectr_E2_E2, work, nstep, &
            &              mpi_double_precision, mpi_sum, mpi_g_world, ierr )
       Spectr_E2_E2 = work

       work = 0.0d0
       call mpi_allreduce( Spectr_E1_E2, work, nstep, &
            &              mpi_double_precision, mpi_sum, mpi_g_world, ierr )
       Spectr_E1_E2 = work
       deallocate( work )
    endif

    fac1 = PAI4 *PAI              ! latter : i pi delta

    Spectr_E1_E1(:) = Spectr_E1_E1(:) /dble(kv3_fbz/nspin) /univol *fac1
    Spectr_E2_E2(:) = Spectr_E2_E2(:) /dble(kv3_fbz/nspin) /univol *fac1
    Spectr_E1_E2(:) = Spectr_E1_E2(:) /dble(kv3_fbz/nspin) /univol *fac1
    if ( nspin == 1 ) then
       Spectr_E1_E1 = Spectr_E1_E1 *2.0d0
       Spectr_E2_E2 = Spectr_E2_E2 *2.0d0
       Spectr_E1_E2 = Spectr_E1_E2 *2.0d0
    endif

    if ( qnum_l_to_probe > 0 .and. mimic_soc_split_spectrum ) then
       if ( ene_initial_state_splitting /= 0.0d0 ) then
          allocate( work(nstep ) );
          work = 0.0d0
          call split_spectr_by_spinorb_real( nstep, e_step, Spectr_E1_E1, work )
          Spectr_E1_E1 = work

          work = 0.0d0
          call split_spectr_by_spinorb_real( nstep, e_step, Spectr_E2_E2, work )
          Spectr_E2_E2 = work

          work = 0.0d0
          call split_spectr_by_spinorb_real( nstep, e_step, Spectr_E1_E2, work )
          Spectr_E1_E2 = work
       endif
    end if

    if ( mype == 0 ) then
       call m_Files_open_nfeps( "E1_E1" )
       call write_data( nfepsout, nstep, ene, Spectr_E1_E1, "dipole-dipole" )
       call m_Files_close_nfeps

       call m_Files_open_nfeps( "E2_E2" )
       call write_data( nfepsout, nstep, ene, Spectr_E2_E2, &
            &           "quadrupole-quadrupole" )
       call m_Files_close_nfeps

       call m_Files_open_nfeps( "E1_E2" )
       call write_data( nfepsout, nstep, ene, Spectr_E1_E2, &
            &           "dipole-quadrupole" )
       call m_Files_close_nfeps
    endif

    deallocate( Spectr_E1_E1 ); deallocate( Spectr_E2_E2 ); deallocate( Spectr_E1_E2 )

  contains

    subroutine write_data( lun, nstep, e, spectr_in, word )
      integer, intent(in) :: lun, nstep
      real(kind=DP), intent(in) :: e(nstep)
      real(kind=DP), intent(in) :: spectr_in(nstep)
      character*(*), intent(in) :: word

      integer :: istep
      real(kind=DP) :: ene

      write(lun,'(A,A,A)') "# Core level spectrum ( ", trim(word), " )"
      write(lun,'(A,I4,A,I2,A,I2)') "#  atom:", atom_to_probe, &
           &                  ", orbital: n = ", qnum_n_to_probe, &
           &                           ", l =", qnum_l_to_probe
      write(lun,'(A)') "#"
      write(lun,'(A,3F10.5)') "# Polarization: ", evec(1:3)
      write(lun,'(A,3F10.5)') "# Wave vector : ", kvec(1:3)
      write(lun,'(A)') "#"
      write(lun,'(A)') "#    Energy [eV]      spectrum "
      write(lun,'(A)') "#"

      select case (spectrum_type)
      case (MacroDielectric)         ! epsilon
         Do istep=1, nstep
            ene = e(istep)
            write(lun,'(E18.10,E18.10E2)') ene*Hartree_in_eV, &
                 &                         spectr_in(istep)
         End Do
      case (PACS)         ! cross section
         Do istep=1, nstep
            ene = e(istep)
            write(lun,'(E18.10,E18.10E2)') ene*Hartree_in_eV, &
                 &                         spectr_in(istep) &
                 &                             *ene /InvHyperFineConst *univol
         End Do
      end select
    end subroutine write_data

  end subroutine m_CLS_calc_spectr_fn_l

  subroutine m_CLS_calc_spectr_fn_c( nstep, evec, evec2, kvec, scissor, way_BZ_int, &
       &                             occup_l_ek, ene, e_step, width, spectrum_type )
    integer, intent(in) :: nstep, way_BZ_int, spectrum_type
    real(kind=DP), intent(in) :: evec(3), kvec(3), scissor, ene(nstep), e_step, width
    real(kind=DP), intent(in) :: occup_l_ek(neg,kv3_ek)
    real(kind=DP), intent(in) :: evec2(3)

    integer :: i, ik, jk, iopr, trev
    integer :: ib1, ib2, is, ixyz1, ixyz2, ierr
    real(kind=DP) :: evec_rotated(3), kvec_rotated(3), weight
    real(kind=DP) :: evec2_rotated(3)
    real(kind=DP) :: ediff, odiff, occ1, occ2, txyz(3)
    real(kind=DP) :: Delta03 = 1.0D-3

    complex(kind=CMPLDP) :: zrho1p, zrho2p
    complex(kind=CMPLDP) :: zrho1m, zrho2m
    real(kind=DP) :: trm2p_E1_E1, trm2p_E2_E2, trm2p_E1_E2
    real(kind=DP) :: trm2m_E1_E1, trm2m_E2_E2, trm2m_E1_E2

    real(kind=DP), allocatable :: work(:,:)
    real(kind=DP), allocatable :: spectr_E1_E1(:,:), spectr_E2_E2(:,:), spectr_E1_E2(:,:)

    real*8 :: fac1, fac2, wavenum_xray

    allocate( spectr_E1_E1(nstep,2) ); spectr_E1_E1 = 0.0d0
    allocate( spectr_E2_E2(nstep,2) ); spectr_E2_E2 = 0.0d0
    allocate( spectr_E1_E2(nstep,2) ); spectr_E1_E2 = 0.0d0

    Do ib1=1, num_core_states
       Do ik=1, kv3_ek, ndim_spinor
          weight = kv3_ek *qwgt_ek(ik) /dble(ndim_spinor)

          Do i=1, num_star_of_k(ik)
             jk = star_of_k(ik,i)
             iopr = iopr_k_fbz_to_ibz(jk)
             trev = trev_k_fbz_to_ibz(jk)

             kvec_rotated = matmul( op(:,:,iopr), kvec(:) )
             evec_rotated = matmul( op(:,:,iopr), evec(:) )
             evec2_rotated = matmul( op(:,:,iopr), evec2(:) )
             if ( noncol .and. magmom_dir_inversion_opr_flag(iopr) == -1 ) then
                evec2_rotated = -evec2_rotated
             endif

             if ( trev == 1 ) then
                kvec_rotated = -kvec_rotated
                evec_rotated = -evec_rotated
                evec2_rotated = -evec2_rotated
             endif

             txyz(1:3) = tau(1:3,iopr,BUCS)*PAI2

             Do ib2=1, neg
                if ( map_e(ib2) /= myrank_e ) cycle
                occ2 = occup_l_ek(map_z(ib2),ik) /weight

                if ( sw_v2c_xes == ON ) then
                   if ( occ2 < Delta03 ) cycle
                   odiff = occ2
                else
!                 if ( eko_ek(ib2,ik) < efermi ) cycle
                   if ( 1.0D0 -occ2 < Delta03 ) cycle
                   odiff = 1.0d0 -occ2
                endif

                ediff = eko_ek(ib2,ik) +scissor -ene_core_states(ib1)
                wavenum_xray = ediff /InvHyperFineConst
                fac1 = wavenum_xray /2.0d0
                fac2 = fac1 *fac1

                zrho1p = 0.0d0;      zrho2p = 0.0d0
                zrho1m = 0.0d0;      zrho2m = 0.0d0
                Do is=1, ndim_spinor
                   Do ixyz1=1, 3
                      zrho1p  = zrho1p &
                           & + trm_dipole( map_z(ib2), ib1, ik+is-1, ixyz1 ) &
                           &   *( evec_rotated(ixyz1) +zi*evec2_rotated(ixyz1) )
                      zrho1m = zrho1m &
                           & + trm_dipole( map_z(ib2), ib1, ik+is-1, ixyz1 ) &
                           &   *( evec_rotated(ixyz1) -zi*evec2_rotated(ixyz1) )
                   End Do
                   Do ixyz1=1, 3
                      Do ixyz2=1, 3
                         zrho2p = zrho2p &
                              & + trm_quadrupole( map_z(ib2),ib1,ik+is-1,ixyz1,ixyz2 ) &
                              &   *(evec_rotated(ixyz1) +zi*evec2_rotated(ixyz1) ) &
                              &   *kvec_rotated(ixyz2)
                         zrho2m = zrho2m &
                              & + trm_quadrupole( map_z(ib2),ib1,ik+is-1,ixyz1,ixyz2 ) &
                              &   *(evec_rotated(ixyz1) -zi*evec2_rotated(ixyz1) ) &
                              &   *kvec_rotated(ixyz2)
                      End Do
                   ENd Do
                End do

                zrho2p = zrho2p *zi
                zrho2m = zrho2m *zi

                trm2p_E1_E1 = conjg( zrho1p ) *zrho1p/2.0d0
                trm2p_E2_E2 = conjg( zrho2p ) *zrho2p *fac2 /2.0d0
                trm2p_E1_E2 = ( conjg( zrho2p ) *zrho1p  &
                     &         +conjg( zrho1p ) *zrho2p ) *fac1 /2.0d0

                trm2m_E1_E1 = conjg( zrho1m ) *zrho1m /2.0d0
                trm2m_E2_E2 = conjg( zrho2m ) *zrho2m *fac2 /2.0d0
                trm2m_E1_E2 = ( conjg( zrho1m  ) *zrho2m &
                     &         +conjg( zrho2m  ) *zrho1m ) *fac1 /2.0d0

                call add_spectral_function( ediff, odiff, trm2m_E1_E1, way_BZ_int, &
                     &                      ene, e_step, width, &
                     &                      nstep, Spectr_E1_E1(:,1) )
                call add_spectral_function( ediff, odiff, trm2p_E1_E1, way_BZ_int, &
                     &                      ene, e_step, width, &
                     &                      nstep, Spectr_E1_E1(:,2) )

                call add_spectral_function( ediff, odiff, trm2m_E2_E2, way_BZ_int,&
                     &                      ene, e_step, width, &
                     &                      nstep, Spectr_E2_E2(:,1) )
                call add_spectral_function( ediff, odiff, trm2p_E2_E2, way_BZ_int,&
                     &                      ene, e_step, width, &
                     &                      nstep, Spectr_E2_E2(:,2) )

                call add_spectral_function( ediff, odiff, trm2m_E1_E2, way_BZ_int,&
                     &                      ene, e_step, width, &
                     &                      nstep, Spectr_E1_E2(:,1) )
                call add_spectral_function( ediff, odiff, trm2p_E1_E2, way_BZ_int,&
                     &                      ene, e_step, width, &
                     &                      nstep, Spectr_E1_E2(:,2) )
             End do
          End Do
       End do
    End Do

    if ( npes > 1 ) then
       allocate( work(nstep,2) ); work = 0.0d0
       call mpi_allreduce( Spectr_E1_E1, work, nstep*2, &
            &              mpi_double_precision, mpi_sum, mpi_g_world, ierr )
       Spectr_E1_E1 = work

       work = 0.0d0
       call mpi_allreduce( Spectr_E2_E2, work, nstep*2, &
            &              mpi_double_precision, mpi_sum, mpi_g_world, ierr )
       Spectr_E2_E2 = work

       work = 0.0d0
       call mpi_allreduce( Spectr_E1_E2, work, nstep*2, &
            &              mpi_double_precision, mpi_sum, mpi_g_world, ierr )
       Spectr_E1_E2 = work
       deallocate( work )
    endif

    fac1 = PAI4 *PAI              ! latter : i pi delta

    Spectr_E1_E1(:,:) = Spectr_E1_E1(:,:) /dble(kv3_fbz/nspin) /univol *fac1
    Spectr_E2_E2(:,:) = Spectr_E2_E2(:,:) /dble(kv3_fbz/nspin) /univol *fac1
    Spectr_E1_E2(:,:) = Spectr_E1_E2(:,:) /dble(kv3_fbz/nspin) /univol *fac1
    if ( nspin == 1 ) then
       Spectr_E1_E1 = Spectr_E1_E1 *2.0d0
       Spectr_E2_E2 = Spectr_E2_E2 *2.0d0
       Spectr_E1_E2 = Spectr_E1_E2 *2.0d0
    endif

    if ( qnum_l_to_probe > 0 .and. mimic_soc_split_spectrum ) then
       if ( ene_initial_state_splitting /= 0.0d0 ) then
          allocate( work(nstep,2) );
          work = 0.0d0
          call split_spectr_by_spinorb_real( nstep, e_step, Spectr_E1_E1(:,1),work(:,1))
          call split_spectr_by_spinorb_real( nstep, e_step, Spectr_E1_E1(:,2),work(:,2))
          Spectr_E1_E1 = work

          work = 0.0d0
          call split_spectr_by_spinorb_real( nstep, e_step, Spectr_E2_E2(:,1),work(:,1))
          call split_spectr_by_spinorb_real( nstep, e_step, Spectr_E2_E2(:,2),work(:,2))
          Spectr_E2_E2 = work

          work = 0.0d0
          call split_spectr_by_spinorb_real( nstep, e_step, Spectr_E1_E2(:,1),work(:,1))
          call split_spectr_by_spinorb_real( nstep, e_step, Spectr_E1_E2(:,2),work(:,2))
          Spectr_E1_E2 = work
       endif
    end if

    if ( mype == 0 ) then
       call m_Files_open_nfeps( "E1_E1" )
       call write_data( nfepsout, nstep, ene, Spectr_E1_E1, "dipole-dipole" )
       call m_Files_close_nfeps

       call m_Files_open_nfeps( "E2_E2" )
       call write_data( nfepsout, nstep, ene, Spectr_E2_E2, &
            &           "quadrupole-quadrupole" )
       call m_Files_close_nfeps

       call m_Files_open_nfeps( "E1_E2" )
       call write_data( nfepsout, nstep, ene, Spectr_E1_E2, &
            &           "dipole-quadrupole" )
       call m_Files_close_nfeps
    endif

    deallocate( Spectr_E1_E1 ); deallocate( Spectr_E2_E2 ); deallocate( Spectr_E1_E2 )

  contains

    subroutine write_data( lun, nstep, e, spectr_in, word )
      integer, intent(in) :: lun, nstep
      real(kind=DP), intent(in) :: e(nstep)
      real(kind=DP), intent(in) :: spectr_in(nstep,2)
      character*(*), intent(in) :: word

      integer :: istep
      real(kind=DP) :: ene, c1

      write(lun,'(A,A,A)') "# Core level spectrum ( ", trim(word), " )"
      write(lun,'(A,I4,A,I2,A,I2)') "#  atom:", atom_to_probe, &
           &                  ", orbital: n = ", qnum_n_to_probe, &
           &                           ", l =", qnum_l_to_probe
      write(lun,'(A)') "#"
      write(lun,'(A,3F10.5)') "# Wave vector : ", kvec(1:3)
      write(lun,'(A)') "#"
      write(lun,'(A)') "#    Energy [eV]      spectrum (left,-)     (right,+)      (CD)"
      write(lun,'(A)') "#"

      select case (spectrum_type)
      case (MacroDielectric)         ! epsilon
         Do istep=1, nstep
            ene = e(istep)
            write(lun,'(E18.10,3E18.10E2)') ene*Hartree_in_eV, &
                 &                         spectr_in(istep,1:2), &
                 &                         spectr_in(istep,2) -spectr_in(istep,1)
         End Do
      case (PACS)         ! cross section
         Do istep=1, nstep
            ene = e(istep)
            c1 = ene /InvHyperFineConst *univol
            write(lun,'(E18.10,3E18.10E2)') ene*Hartree_in_eV, &
                 &                          c1 *spectr_in(istep,1:2), &
                 &                          c1 *(spectr_in(istep,2)-spectr_in(istep,1))
         End Do
      end select
    end subroutine write_data

  end subroutine m_CLS_calc_spectr_fn_c

  subroutine split_spectr_by_spinorb_real( nstep, e_step, spectr_in, spectr_out )
    integer, intent(in) :: nstep
    real(kind=DP), intent(in) :: e_step
    real(kind=DP), intent(in) :: spectr_in(nstep)
    real(kind=DP), intent(out) :: spectr_out(nstep)

    integer :: ngeta, istep, jstep
    real(kind=DP) :: c1, c2

    ngeta = nint( ene_initial_state_splitting /e_step )
    c1 = qnum_l_to_probe /( 2.0d0*qnum_l_to_probe +1.0d0 )
    c2 = 1.0d0 -c1

    spectr_out(:) = 0.0d0
    Do istep=1, nstep
       spectr_out(istep) = spectr_out(istep) + c2 *spectr_in(istep)
       jstep = istep +ngeta
       if ( jstep <=nstep ) then
          spectr_out(jstep) = spectr_out(jstep) + c1 *spectr_in(istep)
       endif
    End Do
  end subroutine split_spectr_by_spinorb_real

  subroutine add_spectral_function( e_center, weight, trm2, way_BZ_int, &
       &                            e, e_step, width, nstep, SpectrFn )
    integer, intent(in) :: nstep, way_BZ_int
    real(kind=DP), intent(in) :: e_center, weight, e_step, width
    real(kind=DP), intent(in) :: e(nstep)
    real(kind=DP), intent(in) :: trm2
    real(kind=DP), intent(inout) :: SpectrFn( nstep )

    real(kind=DP) :: e_low
    real(kind=DP) :: gaussian_extent = 4.0d0

    e_low = e(1)
    select case (way_BZ_int)
    case (PARABOLIC_B)
       call Parabolic_broadening
!       stop "Does not work"
    case (L_TETRAHEDRON)
       call phase_error_with_msg(nfout,"Does not work",__LINE__,__FILE__)
    case (GAUSSIAN_B)
       call Gaussian_broadening
    case (LORENTZIAN_B)
       call Lorentzian_broadening
    end select

  contains

    subroutine Parabolic_broadening
      integer :: i
      real(kind=DP) :: c2, ctmp

      Do i=1, nstep
         call width2( e(i), e_center, width, c2, ctmp )
         SpectrFn(i) = SpectrFn(i) +weight *trm2 *c2
      End do

    end subroutine Parabolic_broadening

    subroutine Gaussian_broadening           ! width == sqrt(2)*sigma,
      integer :: istart, iend, i
      real(kind=DP) :: sigma, sig2, c1, c2, c3, e_min, e_max, e0

      sigma = width /sqrt(2.0d0)
      e0 = e_center -e_low
      e_min = e0 -sigma *gaussian_extent;
      e_max = e0 +sigma *gaussian_extent
      istart = nint( e_min /e_step ) +1;  iend = nint( e_max /e_step ) +1

      istart = max( 1, istart );      istart = min( nstep, istart )
      iend   = min( nstep, iend );    iend   = max( 1, iend )

      sig2 = sigma**2;   c1 = sqrt( PAI2 *sig2 )
      Do i=istart, iend
         c2 = ( e(i) -e_center )**2 /sig2 /2.0d0
         c3 = exp( -c2 ) /c1
         SpectrFn(i) = SpectrFn(i) +weight *trm2 *c3
      End do

    end subroutine Gaussian_broadening

    subroutine Lorentzian_broadening
      integer :: i
      real(kind=DP) :: e0, c1, c2, c3, eta

      eta = width /2.0d0
      !
      c3 = eta /PAI
!      e0 = e_center -e_low
      e0 = e_center

      Do i=1, nstep
         c1 = ( e(i) -e0 )**2 + eta**2
         c2 = c3 /c1
         SpectrFn(i) = SpectrFn(i) +weight *trm2 *c2
      End do
    end subroutine Lorentzian_broadening

  end subroutine add_spectral_function

  subroutine m_CLS_alloc_dipquad
    integer :: nlmt_val

    if ( sw_use_add_proj == ON ) then
       nlmt_val = nlmt +nlmt_add
    else
       nlmt_val = nlmt
    endif

    allocate( rdipole_val_core( nlmt_val, 2*qnum_l_to_probe+1, 3 ) )
    allocate( rquadrupole_val_core( nlmt_val, 2*qnum_l_to_probe+1,3,3) )
    rdipole_val_core = 0.0d0;   rquadrupole_val_core = 0.0d0;

  end subroutine m_CLS_alloc_dipquad

  subroutine m_CLS_dealloc_dipquad
    deallocate( rdipole_val_core )
    deallocate( rquadrupole_val_core )
  end subroutine m_CLS_dealloc_dipquad

  subroutine m_CLS_set_dipquad_elements
    integer :: orb_index, it, ilt1, il1, tau1, ir, ierr
    integer :: ilmt2, il2, im2, tau2, nsph_core, nsph_val, im1
    real(kind=DP) :: c1, c2, ctmp
    real(kind=DP), allocatable :: RadIntegral_order(:,:,:), wos(:)
!
    call set_rhat_cylm
    call set_rhat_rylm
    call set_rhat2_rylm
!
    orb_index = m_CLS_find_orb_index_to_probe()
    it = ityp( atom_to_probe )
    allocate( RadIntegral_order(2,nloc,ntau) )
    allocate( wos(mmesh) )
    call set_weight_exp( ierr, 1, nmesh(it), radr_paw(:,it), wos )

    Do il1=1, lpsmax(it)
       Do tau1=1, itau(il1,it)
          c1 = 0.0d0;  c2 = 0.0d0
          Do ir=1, nmesh(it)
             ctmp = psir_core_ae_wfns(ir,orb_index) *psirpw(ir,il1,tau1,it)
             c1 = c1 +wos(ir) *ctmp *radr_paw(ir,it)
             c2 = c2 +wos(ir) *ctmp *radr_paw(ir,it)**2
          End do
          RadIntegral_order(1,il1,tau1) = c1
          RadIntegral_order(2,il1,tau1) = c2
       End Do
    End Do

    Do im1=1, 2*qnum_l_to_probe +1
       nsph_core = qnum_l_to_probe**2 +im1
       Do ilmt2=1, ilmt(it)        ! valence
          il2=ltp(ilmt2,it);  im2=mtp(ilmt2,it);   tau2=taup(ilmt2,it)
          nsph_val = (il2 -1)**2 +im2
          rdipole_val_core(ilmt2,im1,1:3) = RadIntegral_order(1,il2,tau2) &
               &               *rhat_rylm(nsph_val,nsph_core,1:3)
          rquadrupole_val_core(ilmt2,im1,1:3,1:3) = RadIntegral_order(2,il2,tau2) &
               &               *rhat2_rylm(nsph_val,nsph_core,1:3,1:3)
       End do
       if ( sw_use_add_proj == ON ) then
          Do ilmt2=1, ilmt_add(it)
             il2=ltp_add(ilmt2,it);  im2=mtp_add(ilmt2,it);   tau2=1
             nsph_val = (il2 -1)**2 +im2
             rdipole_val_core(nlmt+ilmt2,im1,1:3) &
                  &              = RadIntegral_order(1,il2,tau2) &
                  &                 *rhat_rylm(nsph_val,nsph_core,1:3)
             rquadrupole_val_core(nlmt+ilmt2,im1,1:3,1:3) &
                  &              = RadIntegral_order(2,il2,tau2) &
                  &                 *rhat2_rylm(nsph_val,nsph_core,1:3,1:3)
          End do
       endif
    End Do

    deallocate( RadIntegral_order )
    deallocate( wos )
  end subroutine m_CLS_set_dipquad_elements

  subroutine set_rhat2_rylm
    integer :: ixyz1, ixyz2, i, j, k
    real(kind=DP) :: csum

    Do ixyz2=1, 3
       Do ixyz1=1, 3
          Do i=1, nsph_max
             Do j=1, nsph_max
                csum = 0.0d0
                Do k=1, nsph_max
                   csum = csum +rhat_rylm(i,k,ixyz1)*rhat_rylm(k,j,ixyz2)
                End Do
                rhat2_rylm(i,j,ixyz1,ixyz2) = csum
             End Do
          End Do
       ENd Do
    End Do
  end subroutine set_rhat2_rylm

  subroutine set_rhat_cylm         ! rhat = x^, y^, z^ ( x^=x/r )
    integer :: il, im, nsph1, nsph2, nsph3, ind1, ind2
    real(kind=DP) :: c1, c2, a1, a2, a3, a4, a5, a6
    integer :: get_nspher_cylm

    rhat_cylm = 0.0d0

    Do il=0, lmax
       c1 = (2*il+1)*(2*il+3)
       c2 = 0.0d0
       if ( il > 0 ) c2 = (2*il-1)*(2*il+1)

       Do im=-il, il
          nsph1 = get_nspher_cylm( il, im )

          if ( abs(im+1) <= il+1 ) then
             a1 =  -sqrt( (il+im+1)*(il+im+2)/c1 )
             nsph2 = get_nspher_cylm( il+1, im+1 )
             if ( nsph2 <= nsph_max ) then
                rhat_cylm(nsph2,nsph1,1) = a1 /2.0d0
                rhat_cylm(nsph2,nsph1,2) = -a1 /2.0d0 *zi
             endif
          endif

          if ( abs(im+1) <= il-1 ) then
             a2 = sqrt( (il-im-1)*(il-im)/c2 )
             nsph2 = get_nspher_cylm( il-1, im+1 )
             if ( nsph2 <= nsph_max ) then
                rhat_cylm(nsph2,nsph1,1) = a2 /2.0d0
                rhat_cylm(nsph2,nsph1,2) = -a2 /2.0d0 *zi
             endif
          endif

          if ( abs(im-1) <= il+1 ) then
             a3 =  sqrt( (il-im+1)*(il-im+2)/c1 )
             nsph2 = get_nspher_cylm( il+1, im-1 )
             if ( nsph2 <= nsph_max ) then
                rhat_cylm(nsph2,nsph1,1) = a3 /2.0d0
                rhat_cylm(nsph2,nsph1,2) = a3 /2.0d0 *zi
             endif
          endif

          if ( abs(im-1) <= il-1 ) then
             a4 =  -sqrt( (il+im-1)*(il+im)/c2 )
             nsph2 = get_nspher_cylm( il-1, im-1 )
             if ( nsph2 <= nsph_max ) then
                rhat_cylm(nsph2,nsph1,1) = a4 /2.0d0
                rhat_cylm(nsph2,nsph1,2) = a4 /2.0d0 *zi
             endif
          endif

          if ( abs(im) <= il+1 ) then
             a5 =  sqrt( (il-im+1)*(il+im+1)/c1 )
             nsph2 = get_nspher_cylm( il+1, im )
             if ( nsph2 <= nsph_max ) then
                rhat_cylm(nsph2,nsph1,3) = a5
             endif
          endif

          if ( abs(im) <= il-1 ) then
             a6 =  sqrt( (il-im)*(il+im)/c2 )
             nsph2 = get_nspher_cylm( il-1, im )
             if ( nsph2 <= nsph_max ) then
                rhat_cylm(nsph2,nsph1,3) = a6
             endif
          endif

       End Do
    End Do
  end subroutine set_rhat_cylm

  subroutine set_rhat_rylm
    integer :: il1, il2, im1, im2, im3, im4, ind1, ind2
    integer :: nsph3, nsph4
    integer :: get_nspher_cylm
!    external :: get_nspher_cylm

    complex(kind=CMPLDP) :: zsum(3), ztmp1, ztmp2

    Do il1=0, lmax
       Do im1=1, 2*il1 +1
          ind1 = il1**2 +im1

          Do il2=0, lmax
             Do im2=1, 2*il2 +1
                ind2 = il2**2 +im2
                zsum = 0.0d0

                Do im3=-il1, il1
                   if ( il1 == 0 ) then
                      ztmp1 = MatU_ylm_RC_L0(im1,im3)
                   else if ( il1 == 1  ) then
                      ztmp1 = MatU_ylm_RC_L1(im1,im3)
                   else if ( il1 == 2  ) then
                      ztmp1 = MatU_ylm_RC_L2(im1,im3)
                   else if ( il1 == 3  ) then
                      ztmp1 = MatU_ylm_RC_L3(im1,im3)
                   else if ( il1 == 4  ) then
                      ztmp1 = MatU_ylm_RC_L4(im1,im3)
                   endif

                   Do im4=-il2, il2
                      if ( il2 == 0 ) then
                         ztmp2 = MatU_ylm_RC_L0(im2,im4)
                      else if ( il2 == 1  ) then
                         ztmp2 = MatU_ylm_RC_L1(im2,im4)
                      else if ( il2 == 2  ) then
                         ztmp2 = MatU_ylm_RC_L2(im2,im4)
                      else if ( il2 == 3  ) then
                         ztmp2 = MatU_ylm_RC_L3(im2,im4)
                      else if ( il2 == 4  ) then
                         ztmp2 = MatU_ylm_RC_L4(im2,im4)
                      endif

                      nsph3 = get_nspher_cylm( il1, im3 )
                      nsph4 = get_nspher_cylm( il2, im4 )

                      zsum(1) = zsum(1) +ztmp1 *rhat_cylm(nsph3,nsph4,1) *conjg(ztmp2)
                      zsum(2) = zsum(2) +ztmp1 *rhat_cylm(nsph3,nsph4,2) *conjg(ztmp2)
                      zsum(3) = zsum(3) +ztmp1 *rhat_cylm(nsph3,nsph4,3) *conjg(ztmp2)
                   End Do
                ENd Do
                rhat_rylm(ind1,ind2,1:3) = zsum(1:3)
             End Do
          End Do
       End Do
    End Do
  end subroutine set_rhat_rylm

  integer function get_unused_unitnumber()
    integer :: i
    logical :: op
    do i=500,9999
      inquire(unit=i,opened=op)
      if( .not.op ) then
        get_unused_unitnumber = i
        return
      endif
    enddo
    call phase_error_with_msg(nfout,'  insufficient file handle',__LINE__,__FILE__)
  end function get_unused_unitnumber

end module m_CLS_dipquad
