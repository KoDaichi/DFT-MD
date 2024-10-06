!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 635 $)
!
!  MODULE: m_ES_IO
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
!=======================================================================
!
!     The original version of this set of the computer programs "PHASE"
!  was developed by the members of the Theory Group of Joint Research
!  Center for Atom Technology (JRCAT), based in Tsukuba, in the period
!  1993-2001.
!
!     Since 2002, this set has been tuned and new functions have been
!  added to it as a part of the national project "Frontier Simulation
!  Software for Industrial Science (FSIS)",  which is supported by
!  the IT program of the Ministry of Education, Culture, Sports,
!  Science and Technology (MEXT) of Japan.
!     Since 2006, this program set has been developed as a part of the
!  national project "Revolutionary Simulation Software (RSS21)", which
!  is supported by the next-generation IT program of MEXT of Japan.
!   Since 2013, this program set has been further developed centering on PHASE System
!  Consortium.
!   The activity of development of this program set has been supervised by Takahisa Ohno.
!
#ifdef __TIMER_SUB__
#   define __TIMER_SUB_START(a)  call timer_sta(a)
#   define __TIMER_SUB_STOP(a)   call timer_end(a)
#else
#   define __TIMER_SUB_START(a)
#   define __TIMER_SUB_STOP(a)
#endif
#ifdef __TIMER_IODO__
#   define __TIMER_IODO_START(a)   call timer_sta(a)
#   define __TIMER_IODO_STOP(a)    call timer_end(a)
#else
#   define __TIMER_IODO_START(a)
#   define __TIMER_IODO_STOP(a)
#endif
#ifdef __TIMER_IOCOMM__
#   define __TIMER_IOCOMM_START_w_BARRIER(str,a)   call timer_barrier(str) ;   call timer_sta(a)
#   define __TIMER_IOCOMM_START(a)       call timer_sta(a)
#   define __TIMER_IOCOMM_STOP(a)        call timer_end(a)
#else
#   define __TIMER_IOCOMM_START_w_BARRIER(str,a)
#   define __TIMER_IOCOMM_START(a)
#   define __TIMER_IOCOMM_STOP(a)
#endif

!#define _DEBUG_ESIO_

!
module m_ES_IO
! $Id: m_ES_IO.F90 635 2021-02-26 07:16:10Z jkoga $
  use m_Electronic_Structure, only : zaj_l,neordr,nrvf_ordr,eko_l,occup_l,efermi,efermi_spin,totch&
       &                            ,vnlph_l,vlhxc_l,eko_ek, zaj_l_prev, metalic_system, vbm
  use m_Electronic_Structure, only : m_ES_WF_in_Rspace
  use m_PlaneWaveBasisSet,    only : kgp,kg1,ngabc,nbase,iba
  use m_Kpoints,              only : kv3, vkxyz, vkxyz_ek, kv3_ek, k_symmetry, qwgt, &
       &                             m_Kp_get_nkmesh, m_Kp_get_kptable_bxsf, &
       &                             vkxyz_refcell
  use m_Timing,               only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters,   only : nspin,kimg,neg,num_extra_bands,af,ipri, printable, neg_previous &
       &                           , wf_filetype, wf_title, eigmin_wf, eigmax_wf, ekmode, neg_is_enlarged &
       &                           , icond, fixed_charge_k_parallel, sw_ekzaj, numk_zajsaved, Nw_Psicoef &
       &                           , precision_WFfile
  use m_Const_parameters,     only : DP, SP, CMPLDP, BUCS, OFF, YES, EK, SCF, DENSITY_ONLY &
       &                           , CUBE, VTK, BINARY, GAMMA, GAMMA_base_symmetrization, ONE_BY_ONE &
       &                           , FIXED_CHARGE, FIXED_CHARGE_CONTINUATION, DELTA &
       &                           , EFERMI_VICINITY, ALL_VALUES, GRID, ON
  use m_Parallelization,      only : MPI_CommGroup,mpi_k_world,mpi_e_world,is_kngp,ie_kngp,npes &
       &                           , mype,ierr,map_k, map_ek,ista_e,iend_e,istep_e,map_z, np_e &
       &                           , ista_k,iend_k,myrank_e,myrank_k,map_e,nrank_e &
       &                           , ista_kngp,iend_kngp, nrank_k  &
       &                           , ista_g1k,iend_g1k, np_g1k , myrank_g, nrank_g, np_g1k_prev &
       &                           , nis_g1k,nie_g1k,map_rank_gek  &
       &                           , ista_spin, iend_spin, map_s, myrank_spin, map_eks
  use m_IterationNumbers,     only : nk_in_the_process, nk_converged, nkgroup &
       &                           , first_kpoint_in_this_job, iteration_ionic, iteration_electronic
  use m_FFT,                  only : fft_box_size_WF,nfft
  use m_Crystal_Structure,    only : altv, sw_fix_total_spin, altv_refcell
  use m_Ionic_System,         only : natm,natm2,iatomn,m_IS_pack_all_ions_in_uc
  use m_PseudoPotential,      only : ival
  use m_Crystal_Structure,    only : univol, rltv

! ===================================== added by K. Tagami ============= 11.0
  use m_Control_Parameters,    only : ndim_spinor, noncol, &
       &                              previous_nspin_collinear, &
       &                              previous_nband_collinear
! ====================================================================== 11.0
  use m_ErrorMessages,        only : EOF_REACHED

  use m_Control_Parameters, only : ndim_magmom, ik_wf_squared, &
       &                           ib1_wf_squared, ib2_wf_squared, &
       &                           wf_squared_filetype,  max_projs, proj_attribute, &
       &                           ndim_chgpot, SpinOrbit_Mode, &
       &                           wf_orb_proj_print_format, proj_group, num_proj_elems, &
       &                           sw_band_unfolding
  use m_Const_parameters,     only : Neglected, PAI2, CARTS, DELTA07, HARTREE, BOHR
  use m_Files,              only :  nfwfk_sq, m_Files_open_nfwfksq_noncl, &
       &                            nfwfk_integ_mom, &
       &                            m_Files_open_nfwfk_integ_mom, &
       &                            m_Files_close_nfwfk_integ_mom, &
       &                            m_Files_open_nfwfk_orb_proj, &
       &                            m_Files_close_nfwfk_orb_proj, &
       &                            nfwfk_orb_proj
  use m_PseudoPotential,   only : nlmt, ilmt, lmta, q, &
       &                          nlmta_phi, nlmtt_phi, qorb, m_PP_tell_iorb_lmtt, &
       &                          m_PP_tell_iorb_ia_l_m_tau, ilmt_phi, &
       &                          mtp_phi, lmta_phi, ltp_phi, taup_phi
  use m_Nonlocal_Potential,   only : norm_phig
  use m_Charge_Density,    only : chgq_l, hsr, hsi, &
       &                          m_CD_softpart_ktsub_noncl, &
       &                          m_CD_hardpart_ktsub_noncl, &
       &                          m_CD_alloc_rspace_charge, &
       &                          m_CD_dealloc_rspace_charge, &
       &                          m_CD_rspace_charge_noncl
  use m_Ionic_System,      only : ityp, iproj_group
  use m_Electronic_Structure,  only : fsr_l, fsi_l, compr_l, compi_l
  use m_ES_Noncollinear,   only : m_ES_set_Pauli_Matrix
  use m_SpinOrbit_Potential,  only :  MatU_ylm_RC_L0,  MatU_ylm_RC_L1,  MatU_ylm_RC_L2, &
       &                              MatU_ylm_RC_L3

! ==== EXP_CELLOPT ==== 2015/09/24
  use m_PlaneWaveBasisSet,    only : kg1_prev
! ===================== 2015/09/24
  use mpi

  implicit none
!  include 'mpif.h'

  integer istatus(mpi_status_size)

  real(kind=SP), allocatable, dimension(:,:)  :: wf_l   ! work wave functions
  real(kind=DP), allocatable, dimension(:,:)  :: wfdp_l ! work wave functions

!  1.  m_ESIO_rd_EigenValues_etc    <-(Initial_Electronic_Structure)
!  2.  m_ESIO_wd_EigenValues_etc    <-(WriteDownData_onto_Files, Postprocessing)
!  3.  m_ESIO_wd_EigenValues        <-(WriteDownData_onto_Files, Convergence_Check, Postprocessing)
!  4.  m_ESIO_wd_EigenValues_ek     <-(WriteDownData_onto_Files)
!  5.  m_ESIO_wd_vlhxc              <-(Postprocessing)
!  6.  m_ESIO_rd_WFs                <-(Initial_Electronic_Structure, scf_rd_wf_and_chg)
!  7.  m_ESIO_rd_WFs_import_frm_collin    <-(Initial_Electronic_Structure)
!  8.  m_ESIO_wd_WFs                <-(WriteDownData_onto_Files)
!  9.  m_ESIO_wd_WFs_standardout    <-(Renewal_of_WaveFunctions)
! 10.  m_ESIO_rd_WFs_and_EVs_ek     <-(Initial_Electronic_Structure)
! 11.  m_ESIO_rd_EVs_ek             <-(Initial_Electronic_Structure)
! 12.  m_ESIO_wd_Psicoef            <-(WriteDownData_onto_Files)
! 13.  m_ESIO_wd_WFs_and_EVs_ek     <-(WriteDownData_onto_Files, WriteDownData_onto_Files_ek, Convergence_Check)
! 14.  m_ESIO_wd_WFn                <-(Postprocessing)
! 15.  m_ESIO_wd_Efermi             <-(WriteDownData_onto_Files)
! 16.  m_ESIO_rd_Eferm              <-(Initial_Electronic_Structure, scf_rd_wf_and_chg)

contains
  subroutine m_ESIO_rd_EigenValues_etc(nfout,nfcntn_bin,F_CNTN_BIN_partitioned)

    integer, intent(in) :: nfout, nfcntn_bin
    logical, intent(in) :: F_CNTN_BIN_partitioned
    integer  :: ik, ie, ispin
    integer, allocatable, dimension(:,:) :: n1_wk, n2_wk  ! MPI
    real(DP),allocatable, dimension(:,:) :: e1_wk, e2_wk  ! MPI
!!$    read(nfcntn_bin) neordr,nrvf_ordr,eko_l,occup_l,efermi,totch
                                                  __TIMER_SUB_START(1370)
    if(F_CNTN_BIN_partitioned) then
       if(neg_previous /= neg) then
          write(nfout,'(" !! neg_previous /= neg <<m_ESIO_rd_EigenValues_etc>>")')
          write(nfout,'(" !! neg_prevous = ",i8)') neg_previous
          write(nfout,'(" !! neg         = ",i8)') neg
          write(nfout,'(" neg_previous sould be neg when F_CNTN_BIN_in_partitioned is true")')
          call phase_error_with_msg(nfout,' neg_previous sould be neg when F_CNTN_BIN_in_partitioned is true',__LINE__,__FILE__)
       end if
       allocate(n1_wk(neg,iend_k-ista_k+1), n2_wk(neg,iend_k-ista_k+1))
       allocate(e1_wk(np_e,iend_k-ista_k+1),e2_wk(np_e,iend_k-ista_k+1))
       n1_wk = 0; n2_wk = 0
       e1_wk = 0; e2_wk = 0
       ! -- neordr, nrvf_ordr --
                                                  __TIMER_IODO_START(1405)
       read(nfcntn_bin) n1_wk
       read(nfcntn_bin) n2_wk
                                                  __TIMER_IODO_STOP(1405)
                                                  __TIMER_IODO_START(1406)
       do ik = ista_k, iend_k
          neordr(1:neg,ik) = n1_wk(1:neg,ik-ista_k+1)
          nrvf_ordr(1:neg,ik) = n2_wk(1:neg,ik-ista_k+1)
       end do
                                                  __TIMER_IODO_STOP(1406)
       ! -- eko_l, occup_l --
                                                  __TIMER_IODO_START(1407)
       read(nfcntn_bin) e1_wk
       read(nfcntn_bin) e2_wk
                                                  __TIMER_IODO_STOP(1407)
                                                  __TIMER_IODO_START(1408)
       do ik = ista_k, iend_k
          do ie = ista_e, iend_e
             eko_l(map_z(ie),ik) = e1_wk(ie-ista_e+1,ik-ista_k+1)
             occup_l(map_z(ie),ik) = e2_wk(ie-ista_e+1,ik-ista_k+1)
          end do
       end do
                                                  __TIMER_IODO_STOP(1408)
       deallocate(e2_wk,e1_wk,n2_wk,n1_wk)
       ! -- nfermi, totch --
       read(nfcntn_bin) efermi, totch
    else
       allocate(n1_wk(neg_previous,kv3),n2_wk(neg_previous,kv3)) ! MPI
       allocate(e1_wk(neg_previous,kv3),e2_wk(neg_previous,kv3)) ! MPI
       n1_wk = 0; n2_wk = 0
       e1_wk = 0; e2_wk = 0

       ! -- neordr, nrvf_ordr --
                                                  __TIMER_IODO_START(1409)
       if(mype == 0) read(nfcntn_bin) n1_wk                ! MPI
       if(mype == 0) read(nfcntn_bin) n2_wk                ! MPI
                                                  __TIMER_IODO_STOP(1409)
                                                  __TIMER_IOCOMM_START_w_BARRIER(MPI_CommGroup,1410)
       if(npes > 1) then
          call mpi_bcast(n1_wk,neg_previous*kv3,mpi_integer,0,MPI_CommGroup,ierr) ! MPI
          call mpi_bcast(n2_wk,neg_previous*kv3,mpi_integer,0,MPI_CommGroup,ierr) ! MPI
       endif
                                                  __TIMER_IOCOMM_STOP(1410)
                                                  __TIMER_IODO_START(1411)
       do ispin = ista_spin, iend_spin
       !do ik = ista_k, iend_k                              ! MPI
       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) /= myrank_k) cycle                 ! MPI
          neordr(1:neg_previous,ik) = n1_wk(1:neg_previous,ik)
          nrvf_ordr(1:neg_previous,ik) = n2_wk(1:neg_previous,ik)
          if(neg_previous < neg) then
             do ie = neg_previous+1, neg
                neordr(ie,ik) = ie
                nrvf_ordr(ie,ik) = ie
             end do
          end if
       end do                                              ! MPI
       end do                                              ! MPI
                                                  __TIMER_IODO_STOP(1411)
       ! -- eko_l, occup_l --
                                                  __TIMER_IODO_START(1412)
       if(mype == 0) read(nfcntn_bin) e1_wk                ! MPI
       if(mype == 0) read(nfcntn_bin) e2_wk                ! MPI
                                                  __TIMER_IODO_STOP(1412)
                                                  __TIMER_IOCOMM_START_w_BARRIER(MPI_CommGroup,1413)
       if(npes > 1) then
          call mpi_bcast(e1_wk,neg_previous*kv3,mpi_double_precision,0,MPI_CommGroup,ierr) ! MPI
          call mpi_bcast(e2_wk,neg_previous*kv3,mpi_double_precision,0,MPI_CommGroup,ierr)! MPI
       end if
                                                  __TIMER_IOCOMM_STOP(1413)
                                                  __TIMER_IODO_START(1414)
       do ispin = ista_spin, iend_spin
       !do ik = ista_k, iend_k                              ! MPI
       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) /= myrank_k) cycle                 ! MPI
          do ie = 1, neg_previous                          ! MPI
             if(map_e(ie) == myrank_e) then                ! MPI
                eko_l(map_z(ie),ik) = e1_wk(ie,ik)         ! MPI
                occup_l(map_z(ie),ik) = e2_wk(ie,ik)       ! MPI
             end if
          end do                                           ! MPI
          if(neg_previous < neg) then
             do ie = neg_previous+1, neg
                if(map_e(ie) == myrank_e) then
                   eko_l(map_z(ie),ik) = 1.d+15
                   occup_l(map_z(ie),ik) = 0.d0
                end if
            end do
          end if
       end do                                              ! MPI
       end do                                              ! MPI
                                                  __TIMER_IODO_STOP(1414)
       ! -- nfermi, totch --
       if(mype == 0) read(nfcntn_bin) efermi, totch        ! MPI
       if(npes > 1) then
          call mpi_bcast(efermi,1,mpi_double_precision,0,MPI_CommGroup,ierr) ! MPI
          call mpi_bcast(totch,1, mpi_double_precision,0,MPI_CommGroup,ierr) ! MPI
       end if

       deallocate(e2_wk,e1_wk,n2_wk,n1_wk)
    end if

    if(printable) write(nfout,'(" TOTCH (total charge) = ",f12.6 &
         & ," (= ",8d25.12,") at m_ESIO_rd_EigenValues_etc")') totch,totch

                                                  __TIMER_SUB_STOP(1370)
  end subroutine m_ESIO_rd_EigenValues_etc

  subroutine m_ESIO_wd_EigenValues_etc(nfcntn_bin,F_CNTN_BIN_partitioned,totch_flag)

    integer, intent(in) :: nfcntn_bin
    logical, intent(in) :: F_CNTN_BIN_partitioned
    integer, optional, intent(in) :: totch_flag

    integer  :: ik, ie, ispin, mpi_comm
    integer, allocatable, dimension(:,:) :: n_wk, n2_mpi  ! MPI
    real(DP),allocatable, dimension(:,:) :: e_wk, e2_mpi  ! MPI
    integer  :: id_sname = -1
                                                  __TIMER_SUB_START(1371)
    call tstatc0_begin('m_ESIO_wd_EigenValues_etc ',id_sname)

    mpi_comm = MPI_CommGroup

    if(F_CNTN_BIN_partitioned) then
       allocate(n_wk(neg,iend_k-ista_k+1))
       allocate(e_wk(np_e,iend_k-ista_k+1))
       n_wk = 0; e_wk = 0
       !  -- neordr --
                                                  __TIMER_IODO_START(1415)
       do ik = ista_k, iend_k
          n_wk(1:neg,ik-ista_k+1) = neordr(1:neg,ik)
       end do
                                                  __TIMER_IODO_STOP(1415)
                                                  __TIMER_IODO_START(1416)
       write(nfcntn_bin) n_wk
                                                  __TIMER_IODO_STOP(1416)
       !  -- nrvf_ordr --
                                                  __TIMER_IODO_START(1417)
       do ik = ista_k, iend_k
          n_wk(1:neg,ik-ista_k+1) = nrvf_ordr(1:neg,ik)
       end do
                                                  __TIMER_IODO_STOP(1417)
                                                  __TIMER_IODO_START(1418)
       write(nfcntn_bin) n_wk
                                                  __TIMER_IODO_STOP(1418)
       !  -- eko_l --
       e_wk = 0.d0
                                                  __TIMER_IODO_START(1419)
       do ik = ista_k, iend_k
          do ie = ista_e, iend_e
             e_wk(ie-ista_e+1,ik-ista_k+1) = eko_l(map_z(ie),ik)
!!$             e_wk(ie,ik-ista_k+1) = eko_l(map_z(ie),ik)
          end do
       end do
                                                  __TIMER_IODO_STOP(1419)
                                                  __TIMER_IODO_START(1420)
       write(nfcntn_bin) e_wk
                                                  __TIMER_IODO_STOP(1420)
!!$       if(npes >= 2) then
!!$          call mpi_allreduce(e_wk,e2_mpi,neg*(iend_k-ista_k+1),mpi_double_precision &
!!$               &               , mpi_sum, mpi_k_world)
!!$       else
!!$          e2_mpi = e_wk
!!$       end if
!!$       write(nfcntn_bin) e2_mpi

       !  -- occup_l --
       e_wk = 0.d0
                                                  __TIMER_IODO_START(1421)
       do ik = ista_k, iend_k
          do ie = ista_e, iend_e, istep_e
             e_wk(ie-ista_e+1,ik-ista_k+1) = occup_l(map_z(ie),ik)
!!$             e_wk(ie,ik-ista_k+1) = occup_l(map_z(ie),ik)
          end do
       end do
                                                  __TIMER_IODO_STOP(1421)
                                                  __TIMER_IODO_START(1422)
       write(nfcntn_bin) e_wk
                                                  __TIMER_IODO_STOP(1422)
!!$       if(npes >= 2) then
!!$          call mpi_allreduce(e_wk,e2_mpi,neg*(iend_k-ista_k+1),mpi_double_precision &
!!$               &               , mpi_sum, mpi_k_world)
!!$       else
!!$          e2_mpi = e_wk
!!$       end if
!!$       write(nfcntn_bin) e2_mpi
!!$       deallocate(e2_mpi, e_wk, n_wk)
       deallocate(e_wk, n_wk)
       if(totch_flag == OFF) then
          write(nfcntn_bin) efermi                     ! MPI
       else
          write(nfcntn_bin) efermi,totch               ! MPI
       end if
    else
       allocate(n_wk(neg,kv3)); allocate(n2_mpi(neg,kv3))! MPI
       allocate(e_wk(neg,kv3)); allocate(e2_mpi(neg,kv3))! MPI
       n_wk = 0; n2_mpi = 0
       e_wk = 0; e2_mpi = 0
       !  -- neordr --
       n_wk = 0                                          ! MPI
                                                  __TIMER_IODO_START(1423)
       do ispin = ista_spin, iend_spin
       !do ik = 1, kv3                                     ! MPI
       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) /= myrank_k) cycle                 ! MPI
          n_wk(1:neg,ik) = neordr(1:neg,ik)               ! MPI
       end do                                             ! MPI
       end do                                             ! MPI
                                                  __TIMER_IODO_STOP(1423)
       if(npes >= 2) then
          call mpi_allreduce(n_wk,n2_mpi,neg*kv3,mpi_integer,mpi_sum,mpi_comm,ierr)
          n2_mpi = n2_mpi/nrank_e
       else
          n2_mpi = n_wk
       end if
                                                  __TIMER_IODO_START(1425)
       if(mype == 0) write(nfcntn_bin) n2_mpi             ! MPI ; writing (neordr)
                                                  __TIMER_IODO_STOP(1425)
       !  -- nrvf_ordr --
                                                  __TIMER_IODO_START(1426)
       n_wk = 0                                           ! MPI
       do ispin = ista_spin, iend_spin
       !do ik = 1, kv3                                     ! MPI
       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) /= myrank_k) cycle                 ! MPI
          n_wk(1:neg,ik) = nrvf_ordr(1:neg,ik)            ! MPI
       end do                                             ! MPI
       end do                                             ! MPI
                                                  __TIMER_IODO_STOP(1426)
       if(npes >= 2) then
          call mpi_allreduce(n_wk,n2_mpi,neg*kv3,mpi_integer,mpi_sum,mpi_comm,ierr)
          n2_mpi = n2_mpi/nrank_e
       else
          n2_mpi = n_wk
       end if
                                                  __TIMER_IODO_START(1428)
       if(mype == 0) write(nfcntn_bin) n2_mpi             ! MPI ; writing (nrvf_ordr)
                                                  __TIMER_IODO_STOP(1428)
       !  -- eko_l --
       e_wk = 0.d0                                        ! MPI
                                                  __TIMER_IODO_START(1429)
       do ispin = ista_spin, iend_spin
       !do ik = 1, kv3                                     ! MPI
       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) /= myrank_k) cycle                 ! MPI
          do ie = 1, neg                                  ! MPI
             if(map_e(ie) /= myrank_e) cycle              ! MPI
             e_wk(ie,ik) = eko_l(map_z(ie),ik)            ! MPI
          end do
       end do
       end do
                                                  __TIMER_IODO_STOP(1429)
       if(npes >= 2) then
          call mpi_allreduce(e_wk,e2_mpi,neg*kv3,mpi_double_precision,mpi_sum,mpi_comm,ierr)
       else
          e2_mpi = e_wk
       end if
                                                  __TIMER_IODO_START(1431)
       if(mype == 0) write(nfcntn_bin) e2_mpi             ! MPI ; writing (eko_l)
                                                  __TIMER_IODO_STOP(1431)

       !  -- occup_l --
       e_wk = 0.d0
                                                  __TIMER_IODO_START(1432)
       do ispin = ista_spin, iend_spin
       !do ik = 1, kv3                                     ! MPI
       do ik = ispin, kv3-nspin+ispin, nspin
          if(map_k(ik) /= myrank_k) cycle                 ! MPI
          do ie = 1, neg                                  ! MPI
             if(map_e(ie) /= myrank_e) cycle              ! MPI
             e_wk(ie,ik) = occup_l(map_z(ie),ik)          ! MPI
          end do                                            ! MPI
       end do                                               ! MPI
       end do                                               ! MPI
                                                  __TIMER_IODO_STOP(1432)
       if(npes >= 2) then
          call mpi_allreduce(e_wk,e2_mpi,neg*kv3,mpi_double_precision,mpi_sum,mpi_comm,ierr)
       else
          e2_mpi = e_wk
       end if
                                                  __TIMER_IODO_START(1434)
       if(mype == 0) write(nfcntn_bin) e2_mpi             ! MPI ; writing (occup_l)
                                                  __TIMER_IODO_STOP(1434)
       if(mype == 0) then
          if(totch_flag == OFF) then
             if(metalic_system) then
               write(nfcntn_bin) efermi                     ! MPI
             else
               write(nfcntn_bin) vbm+1.d-10
             endif
          else
             if(metalic_system) then
               write(nfcntn_bin) efermi,totch               ! MPI
             else
               write(nfcntn_bin) vbm+1.d-10,totch
             endif
          end if
       end if
       deallocate(n_wk); deallocate(n2_mpi)              ! MPI
       deallocate(e_wk); deallocate(e2_mpi)              ! MPI
    end if
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(1371)
  end subroutine m_ESIO_wd_EigenValues_etc

  subroutine m_ESIO_wd_EigenValues(nf,iprieigen,nooccupation)

    integer, intent(in)              :: nf
    integer, intent(in)              :: iprieigen
    integer, intent(in)              :: nooccupation
    integer                          :: ie,  ipri0, kv3_i, ks
    integer                          :: hconst_min, lzero_max
    integer, parameter :: NCOLUMN = 6
    integer, parameter :: EIGEN_VALUES = 1, OCCUPATIONS = 2
    integer :: writemode
    real(kind=DP),allocatable, dimension(:,:) :: e_mpi, o_mpi
                                                  __TIMER_SUB_START(1378)
    allocate(e_mpi(neg,kv3)); e_mpi = 0.d0
    allocate(o_mpi(neg,kv3)); o_mpi = 0.d0

    call set_writemode(writemode)  ! ->(writemode) = ALL_VALUES or FERMI_VICINITY
    call get_ipri0(iprieigen,ipri0)

    if(ipri0 >= 2) then
       if(ipri0 >= 3 .and. nf == 6 .and. printable) call wd_neordr()

       call set_kv3_i_and_ks() ! -> kv3_i, ks
#ifndef _DEBUG_WRITE_
       if(writemode == EFERMI_VICINITY .and. kv3_i == kv3) &
            & call cal_vicinity_range(hconst_min,lzero_max) ! -> lzero_max, hconst_min
#endif
!     --- Energy eigen values ---
       call put_kpartArray_into(eko_l,e_mpi)
       if(printable) then
          if(ks == 0 .and. nf==6) call wd_efermi()
!!$          if(ks == 0) call wd_efermi()
          call wd_k_and_values(EIGEN_VALUES)
       end if
    end if
!     --- Occupations ---
    if(ipri0 >= 2 .and. nooccupation /= YES) then
       call put_kpartArray_into(occup_l,o_mpi)
       if(printable) call wd_k_and_values(OCCUPATIONS)
    end if
    deallocate(e_mpi)
    deallocate(o_mpi)
                                                  __TIMER_SUB_STOP(1378)
  contains
    subroutine set_writemode(writemode)
      integer, intent(out) :: writemode
      writemode = ALL_VALUES
#ifndef _DEBUG_WRITE_
      if(nf==6) then
         if((icond == FIXED_CHARGE_CONTINUATION .or. icond == FIXED_CHARGE) .and. &
              & fixed_charge_k_parallel == ONE_BY_ONE) then
            writemode = ALL_VALUES
         else
            if(neg <= NCOLUMN) then
               writemode = ALL_VALUES
            else
               writemode = EFERMI_VICINITY
            end if
         end if
      end if
#endif
    end subroutine set_writemode

    subroutine set_kv3_i_and_ks()
!!$      if(iprieigen>=2 .and. printable) then
      if((icond == FIXED_CHARGE_CONTINUATION .or. icond == FIXED_CHARGE) .and. &
           & fixed_charge_k_parallel == ONE_BY_ONE) then
         kv3_i = kv3_ek - kv3*(nkgroup-1)
         if(kv3_i > kv3) kv3_i = kv3
         ks = max(1,first_kpoint_in_this_job) - 1 + kv3*(nkgroup-1)
      else
         kv3_i = kv3
         ks = 0
      end if

      if(iprieigen>=3 .and. printable) &
           & write(nf,'(" kv3_i, kv3, ks, nkgroup = ",4i8)') kv3_i, kv3, ks, nkgroup
!!$
!!$      call mpi_bcast(kv3_i,1,mpi_integer,0,MPI_CommGroup,ierr)
!!$      call mpi_bcast(ks,1,mpi_integer,0,MPI_CommGroup,ierr)
    end subroutine set_kv3_i_and_ks

    subroutine wd_k_and_values(mode)
      integer, intent(in) :: mode
      integer :: ik, nb
      integer :: ie_s, ie_e, nhw, neg_t
      real(kind=DP) :: hw, hc, hv

#ifndef _DEBUG_WRITE_
      if(writemode == EFERMI_VICINITY) then
         hw = (lzero_max-hconst_min)*0.5d0
         hv = (lzero_max+hconst_min)*0.5d0
         hc = NCOLUMN*0.5d0
         nhw = Int(hw/hc + 1.d0)
         ie_s = max(nint(hv-nhw*hc + DELTA),1)
         ie_e = min(ie_s + nhw*NCOLUMN - 1,neg)
      end if
#endif

#if 1
      if ( nf /= 6 .and. mode == EIGEN_VALUES ) then
         write(nf,*)
         write(nf,*) "kpoint list"
         do ik = 1, kv3, ndim_spinor
            call wd_k_points(ik)
         end do
         write(nf,'(" -----")')
      end if
#endif

      if(mode == EIGEN_VALUES) then
#ifndef _DEBUG_WRITE_
         if(writemode == EFERMI_VICINITY) then
            write(nf,'(" ======  Energy Eigen Values in the vicinity of the Fermi energy level (Range=" &
                 & ,i7," :",i7,") =====")') ie_s, ie_e
         else
#endif
            if(nf==6) write(nf,'(" ======  Energy Eigen Values ======")')
#ifndef _DEBUG_WRITE_
         end if
#endif
      else
#ifndef _DEBUG_WRITE_
         if(writemode == EFERMI_VICINITY) then
            write(nf,'(" ======  Occupations in the vicinity of the Fermi energy level (Range=" &
                 & ,i7," :",i7,") =====")') ie_s, ie_e
         else
#endif
            if(nf==6) write(nf,'(" ======  Occupations ======")')
#ifndef _DEBUG_WRITE_
         end if
#endif
      end if
                                                  __TIMER_IODO_START(1463)
      do ik = 1, kv3_i, ndim_spinor
#ifndef _DEBUG_WRITE_
!!$         if(mode == OCCUPATIONS) e_mpi(:,ik) = e_mpi(:,ik)/(qwgt(ik)*kv3)
!!$         if(mode == OCCUPATIONS) o_mpi(:,ik) = o_mpi(:,ik)/(qwgt(ik)*kv3)
         if(writemode==EFERMI_VICINITY .and. kv3==kv3_i) then
            call wd_k_and_efermi_vicinities(ik,ie_s,ie_e,mode)
         else
#endif
            if(nf /= 6 .and. mode == EIGEN_VALUES) write(nf,'(" ===== energy eigenvalues =====")')
            if(nf /= 6 .and. mode == OCCUPATIONS)  write(nf,'(" ===== occupations =====")')
            if ( noncol ) then
               call wd_k_points_noncl(ik)
            else
               call wd_k_points(ik)
            endif
            neg_t = neg
            if(neg_is_enlarged) neg_t = neg - num_extra_bands
            if(mode == EIGEN_VALUES) then
               write(nf,'(5f16.8)') (e_mpi(nb,ik),nb = 1, neg_t) ! =eko(neordr(nb,ik),ik)
            else if(mode == OCCUPATIONS) then
               write(nf,'(5f16.8)') (o_mpi(nb,ik)/(qwgt(ik)*kv3/ndim_spinor),nb = 1, neg_t) ! =occup(neordr(nb,ik),ik)
            end if
#ifndef _DEBUG_WRITE_
         end if
#endif
      end do
                                                  __TIMER_IODO_STOP(1463)
    end subroutine wd_k_and_values

    subroutine put_kpartArray_into(a_l,a_all)
      real(kind=DP), intent(in), dimension(np_e,ista_k:iend_k) :: a_l
      real(kind=DP), intent(out), dimension(neg,kv3) :: a_all
      integer :: ik, ierr, ie, is
      integer :: ito

                                                  __TIMER_IODO_START(1461)
      a_all = 0.d0
      do is = ista_spin, iend_spin
      do ik = is, kv3-nspin+is, nspin
      !do ik = 1, kv3
         if(map_k(ik) /= myrank_k) cycle
         do ie = 1, neg
            if(map_e(ie) /= myrank_e) cycle
            ito = nrvf_ordr(ie,ik)
            a_all(ito,ik) = a_l(map_z(ie),ik)
         end do
      end do
      enddo
                                                  __TIMER_IODO_STOP(1461)
      if(npes >= 2) then
         call mpi_allreduce(MPI_IN_PLACE,a_all,neg*kv3,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
      end if
     end subroutine put_kpartArray_into

    subroutine wd_neordr()
      integer :: ik
      write(nf,'(" kv3 = ",i8, " neg = ",i8)') kv3,neg
      do ik= 1,kv3
         write(nf,'(" map_k(",i3,") = ",i8)') ik,map_k(ik)
         if(map_k(ik) /= myrank_k) cycle
         write(nf,'(" neordr ik=",i8)') ik
         write(nf,'(10i8)') neordr(1:neg,ik)
      end do
    end subroutine wd_neordr

    subroutine get_ipri0(ipri_in, ipri_out)
      integer, intent(in)  :: ipri_in
      integer, intent(out) :: ipri_out
      if(npes > 1) then
         if(mype == 0) ipri_out = ipri_in
         call mpi_bcast(ipri_out,1,mpi_integer,0,MPI_CommGroup,ierr)
      else
         ipri_out = ipri_in
      end if
    end subroutine get_ipri0

    subroutine wd_k_points(ik)
      integer, intent(in) :: ik
      real(kind=DP) :: vkxyz_wk(3)

      if ( sw_band_unfolding == ON ) then
         vkxyz_wk(1:3) = vkxyz_refcell(ik,1:3,BUCS)
      else
         vkxyz_wk(1:3) = vkxyz(ik,1:3,BUCS)
      endif

      if(nspin == 1) then
#ifdef _EIGENVALUES_IN_OLD_FORMAT_
          write(nf,'(i6,3f16.8)') ik+ks, vkxyz_wk(1:3)
#else
          write(nf,'(" ik = ",i8," (",3f10.6," )")') ik+ks, vkxyz_wk(1:3)
#endif
       else
#ifdef _EIGENVALUES_IN_OLD_FORMAT_
          if(mod(ik,2) == 1) then
             write(nf,'(i6,"    UP ",3f16.8)') ik+ks, vkxyz_wk(1:3)
          else
             write(nf,'(i6,"  DOWN ",3f16.8)') ik+ks, vkxyz_wk(1:3)
          end if
#else
          if(mod(ik,2) == 1) then
             write(nf,'(" ik = ",i8," (",3f10.6," )    UP ")') ik+ks, vkxyz_wk(1:3)
          else
             write(nf,'(" ik = ",i8," (",3f10.6," )  DOWN ")') ik+ks, vkxyz_wk(1:3)
          end if
#endif
       end if
     end subroutine wd_k_points

! ============================== added by K. Tagami ==================== 11.0
     subroutine wd_k_points_noncl(ik)
       integer, intent(in) :: ik
       real(kind=DP) :: vkxyz_wk(3)

       if ( sw_band_unfolding == ON ) then
          vkxyz_wk(1:3) = vkxyz_refcell(ik,1:3,BUCS)
       else
          vkxyz_wk(1:3) = vkxyz(ik,1:3,BUCS)
       endif

#ifdef _EIGENVALUES_IN_OLD_FORMAT_
       write(nf,'(i6,3f16.8)') ik+ks, vkxyz_wk(1:3)
#else
       write(nf,'(" ik = ",i5," (",3f10.6," )")') ik+ks, vkxyz_wk(1:3)
#endif
     end subroutine wd_k_points_noncl
! ====================================================================== 11.0

     subroutine cal_vicinity_range(hconst_min, lzero_max)
       integer, intent(out) :: hconst_min, lzero_max
       integer :: hconst,lzero
       integer :: ik, ie, nb
       integer :: ito

       e_mpi = 0.d0
       do ik = 1, kv3
          if(map_k(ik) /= myrank_k) cycle
          do ie = 1, neg
             if(map_e(ie) /= myrank_e) cycle
             ito = nrvf_ordr(ie,ik)
             e_mpi(ito,ik) = occup_l(map_z(ie),ik)
          end do
       end do
       if(npes >= 2) then
          call mpi_allreduce(MPI_IN_PLACE,e_mpi,neg*kv3,mpi_double_precision &
               &                  ,mpi_sum,MPI_CommGroup,ierr) ! MPI
       end if
       hconst_min = neg
       lzero_max = 0
       do ik = 1, kv3, ndim_spinor
          hconst = neg
          lzero  = 1
          do nb = neg,1,-1
             if(e_mpi(nb,ik)<DELTA) lzero = nb
!!$             if(e_mpi(nb,ik)>=qwgt(ik)*kv3-DELTA ) then
             if(e_mpi(nb,ik)>=qwgt(ik)*kv3/ndim_spinor-DELTA ) then
                hconst = nb
                exit
             end if
          end do
          if(hconst < hconst_min) hconst_min = hconst
          if(lzero_max < lzero)   lzero_max = lzero
       end do
       if(hconst_min > lzero_max) hconst_min = max(lzero_max-1,1)

       if(ipri>=2) write(nf,'(" hconst_min,lzero_max = ",2i8)') hconst_min,lzero_max
     end subroutine cal_vicinity_range

     subroutine wd_k_and_efermi_vicinities(ik,ie_s,ie_e,mode)
       integer, intent(in) :: ik, ie_s,ie_e,mode
       integer :: nb, nbloop, ie, ie1, ie2
       if(hconst_min >= 1) then
!!$          nb = lzero_max-hconst_min+1
          nb = ie_e - ie_s + 1
          nbloop = Int((nb-1)/NCOLUMN+1)
          do ie = 1, nbloop
!!$             ie1 = hconst_min+(ie-1)*NCOLUMN
!!$             ie2 = min(hconst_min+ie*NCOLUMN-1,neg)
             ie1 = max(ie_s+(ie-1)*NCOLUMN,1)
             ie2 = min(ie_s+ie*NCOLUMN-1,neg)
             if(ie == 1) then
                if(mode == EIGEN_VALUES) then
                   write(nf,'(" ik = ",i5," ",8f12.6)') ik, (e_mpi(nb,ik),nb=ie1,ie2)
                else if(mode == OCCUPATIONS) then
                   write(nf,'(" ik = ",i5," ",8f12.6)') ik, &
                        &           (o_mpi(nb,ik)/(qwgt(ik)*kv3/ndim_spinor),nb=ie1,ie2)
                end if
             else
                if(mode == EIGEN_VALUES) then
                   write(nf,'(12x,8f12.6)') (e_mpi(nb,ik),nb=ie1,ie2)
                else if(mode == OCCUPATIONS) then
                   write(nf,'(12x,8f12.6)') &
                        &      (o_mpi(nb,ik)/(qwgt(ik)*kv3/ndim_spinor),nb=ie1,ie2)
                end if
             end if
          end do
       end if
     end subroutine wd_k_and_efermi_vicinities

     subroutine wd_efermi()
       if(nf == 6) write(nf,'(" **** Eigen Values and Occupations ****")')
       write(nf,'(" ** iteration_ionic = ",i8, ", iteration_electronic = ",i8," **")') &
            & iteration_ionic, iteration_electronic
       write(nf,'(" EFermi = ",f16.8)') efermi
       if(sw_fix_total_spin == YES .and. nspin == 2) then
          write(nf,'(" Efermi_spin(1) = ",f16.8, ",  Efermi_spin(2) = ",f16.8)') &
               & efermi_spin(1), efermi_spin(2)
       end if
     end subroutine wd_efermi
  end subroutine m_ESIO_wd_EigenValues

  subroutine m_ESIO_wd_EigenValues_ek(nf,mode)
    integer, intent(in)              :: nf, mode

    real(kind=DP), parameter :: delta = 1.d-12
    real(kind=DP), allocatable, dimension(:) :: eko_t
    integer, allocatable, dimension(:)       :: neordr_t
    integer                     :: ik, ib,jb,ibo,jbo, neg_t

! =========================== added by K. Tagami ================ 11.0
    integer :: ikskip
! =============================================================== 11.0

    allocate(eko_t(neg))
    allocate(neordr_t(neg))
    eko_t = 0; neordr_t = 0

    if(mode == SCF .and. printable) write(nf,'(" ======  Energy Eigen Values ======")')
!!$    do ik = 1, kv3_ek
!!$    do ik = 1, nk_in_the_process

    if(printable) then
! ========================= modified by K. Tagami ============ 11.0
!       write(nf,'(" nk_converged = ",i8)') min(kv3_ek,nk_converged)
!       do ik = 1, kv3_ek
!          call wd_k_points
!       end do
       if ( noncol ) then
          write(nf,'(" nk_converged = ",i8)') min(kv3_ek,nk_converged) /ndim_spinor
          do ik = 1, kv3_ek, ndim_spinor
             call wd_k_points_noncl
          end do
       else
          write(nf,'(" nk_converged = ",i8)') min(kv3_ek,nk_converged)
          do ik = 1, kv3_ek
             call wd_k_points
          end do
       endif
! ============================================================= 11.0

       write(nf,'(" -----")')
    end if

! ====================== added by K. Tagami =================== 11.0
    if ( noncol ) then
       ikskip = ndim_spinor
    else
       ikskip = 1
    endif
! ============================================================ 11.0

! ============================ modified by K. Tagami ============ 11.0
!    do ik = 1, nk_converged
    do ik = 1, nk_converged, ikskip
! =============================================================== 11.0

!       if(sw_ekzaj == OFF .and. ik <= first_kpoint_in_this_job) cycle
       if(ik > kv3_ek) cycle
       if(mode == EK .and. printable) write(nf,'("=== energy eigenvalues ===")')
       eko_t = eko_ek(:,ik)
       if(nspin == 1 .or. (nspin == 2 .and. mod(ik,2) == 1)) &
            & neordr_t(1:neg) = (/(ib,ib=1,neg)/)
       do ib = 1, neg-1
          do jb = ib+1, neg
             ibo = neordr_t(ib)
             jbo = neordr_t(jb)
             if(eko_t(jbo)  < eko_t(ibo)-delta) then        ! MPI
                neordr_t(jb) = ibo
                neordr_t(ib) = jbo
             end if
          end do
       end do
       if(printable) then
! ================================ modified by K. Tagami ========== 11.0
!          call wd_k_points
          if ( noncol ) then
             call wd_k_points_noncl
          else
             call wd_k_points
          endif
! ================================================================= 11.0
          neg_t = neg

          if(neg_is_enlarged) neg_t = neg - num_extra_bands
          if(mode == SCF) then
             write(nf,'(5f16.8)') (eko_t(neordr_t(ib)),ib=1,neg_t)
          else
             write(nf,'(4f18.10)') (eko_t(neordr_t(ib)),ib=1,neg_t)
          end if
       end if
    end do

    deallocate(neordr_t);     deallocate(eko_t)

  contains

    subroutine wd_k_points
      integer :: j, k
      real(kind=DP) :: c1, vkxyz_wk(3)

      if ( sw_band_unfolding == ON ) then
         vkxyz_wk(1:3) = vkxyz_refcell(ik,1:3,BUCS)
      else
         vkxyz_wk(1:3) = vkxyz_ek(ik,1:3,BUCS)
      endif

      if(mode == SCF) then
         if(nspin == 1) then
            write(nf,'(i6,3f18.10)') ik, vkxyz_wk(1:3)
         else
            if(mod(ik,2) == 1) then
               write(nf,'(i6,"    UP ",3f18.10)') ik, vkxyz_wk(1:3)
            else
               write(nf,'(i6,"  DOWN ",3f18.10)') ik, vkxyz_wk(1:3)
            end if
         end if
      else
         if(nspin == 1) then
            write(nf,'(" ik = ",i8," (",3f10.6," )")') ik, vkxyz_wk(1:3)
         else
            if(mod(ik,2) == 1) then
               write(nf,'(" ik = ",i8," (",3f10.6," )    UP ")') ik, vkxyz_wk(1:3)
            else
               write(nf,'(" ik = ",i8," (",3f10.6," )  DOWN ")') ik, vkxyz_wk(1:3)
            end if
         end if
      end if

    end subroutine wd_k_points

! ============================== added by K. Tagami ==================== 11.0
    subroutine wd_k_points_noncl
      integer :: j, k
      real(kind=DP) :: c1, vkxyz_wk(3)

      if ( sw_band_unfolding == ON ) then
         vkxyz_wk(1:3) = vkxyz_refcell(ik,1:3,BUCS)
      else
         vkxyz_wk(1:3) = vkxyz_ek(ik,1:3,BUCS)
      endif

      if (mode == SCF) then
         write(nf,'(i6,3f18.10)') ik, vkxyz_wk(1:3)
      else
         write(nf,'(" ik = ",i8," (",3f10.6," )")') ik, vkxyz_wk(1:3)
      endif
    end subroutine wd_k_points_noncl
! ====================================================================== 11.0

  end subroutine m_ESIO_wd_EigenValues_ek

! ======================================= modified by K. Tagami ========== 11.0
!  subroutine m_ESIO_wd_vlhxc(nfvlc)
!
  subroutine m_ESIO_wd_vlhxc( nfvlc, ismax )
    integer, intent(in)              :: ismax
! ======================================================================= 11.0

    integer, intent(in)              :: nfvlc
    integer                          :: is, ik, i
    real(DP),allocatable, dimension(:,:,:):: vlhxc_mpi,vlhxc_mpi2

    if(npes >= 2) then

! ========================== modiifed by K. Tagami ============= 11.0
!       allocate(vlhxc_mpi(kgp,kimg,nspin)); vlhxc_mpi = 0.d0  ! MPI
!       allocate(vlhxc_mpi2(kgp,kimg,nspin))

       allocate(vlhxc_mpi(kgp,kimg,ismax));
       allocate(vlhxc_mpi2(kgp,kimg,ismax))
! ============================================================= 11.0
       vlhxc_mpi = 0.0d0;       vlhxc_mpi2 = 0.0d0

! ========================== modiifed by K. Tagami ============= 11.0
!       do is = 1, nspin
       do is = 1, ismax
! ============================================================== 11.0
          do ik = 1, kimg
             do i = ista_kngp, iend_kngp
                vlhxc_mpi(i,ik,is) = vlhxc_l(i,ik,is)
             end do
          end do
       end do

! ========================== modiifed by K. Tagami ====================== 11.0
!       call mpi_allreduce(vlhxc_mpi,vlhxc_mpi2,kgp*kimg*nspin &
!            &     , mpi_double_precision, mpi_sum, MPI_CommGroup,ierr)
       call mpi_allreduce( vlhxc_mpi, vlhxc_mpi2, kgp*kimg*ismax, &
            &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
! ====================================================================== 11.0

       if (mype == 0) write(nfvlc) vlhxc_mpi2

       deallocate(vlhxc_mpi); deallocate(vlhxc_mpi2)
    else
       write(nfvlc) vlhxc_l
    end if
  end subroutine m_ESIO_wd_vlhxc

  subroutine m_ESIO_rd_WFs(nfout,nfzaj, F_ZAJ_partitioned, prev_zaj)
    integer, intent(in) :: nfout, nfzaj
    logical, intent(in) :: F_ZAJ_partitioned
    logical, intent(in), optional :: prev_zaj
    logical    :: pzaj
    integer    :: ik,ib,ri, i, ispin
    integer    :: id_sname = -1
    integer    :: ierror
    integer    :: kg1t
                                                  __TIMER_SUB_START(1372)
    call tstatc0_begin('m_ESIO_rd_WFs ',id_sname)
    pzaj = .false.
    if(present(prev_zaj)) pzaj = prev_zaj

    kg1t = kg1
    if(pzaj) then
      kg1t = kg1_prev
    endif
    if(precision_WFfile==SP) then
       if(ipri >= 1) write(nfout,*) ' !D Reading zaj (single_precision)'
    else
       if(ipri >= 1) write(nfout,*) ' !D Reading zaj (double_precision)'
    end if
    if(precision_WFfile==SP) then
       allocate(wf_l(kg1t,kimg)); wf_l = 0.d0
    else
       allocate(wfdp_l(kg1t,kimg)); wfdp_l = 0.d0
    end if
    rewind nfzaj
    if(F_ZAJ_partitioned) then
       do ik = ista_k, iend_k, af+1        ! MPI
                                                  __TIMER_IODO_START(1435)
          do ib = ista_e, iend_e, istep_e  ! MPI
             if(ib > neg_previous) cycle
                                                  __TIMER_IODO_START(1436)
           if(precision_WFfile==SP) then
             read(nfzaj) wf_l
                                                  __TIMER_IODO_STOP(1436)
             if(kimg == 1) then
                if(pzaj)then
                  do i = 1, kg1t
                     zaj_l_prev(i,map_z(ib),ik,1) = wf_l(i,1)
                  end do
                else
                  do i = 1, kg1t
                     zaj_l(i,map_z(ib),ik,1) = wf_l(i,1)
                  end do
                endif
             else if(kimg==2) then
                if(pzaj)then
                  do i = 1, kg1t
                     zaj_l_prev(i,map_z(ib),ik,1) = wf_l(i,1)
                     zaj_l_prev(i,map_z(ib),ik,2) = wf_l(i,2)
                  end do
                else
                  do i = 1, kg1t
                     zaj_l(i,map_z(ib),ik,1) = wf_l(i,1)
                     zaj_l(i,map_z(ib),ik,2) = wf_l(i,2)
                  end do
                endif
             end if
           else if(precision_WFfile==DP) then
             read(nfzaj) wfdp_l
                                                  __TIMER_IODO_STOP(1436)
             if(kimg == 1) then
                if(pzaj)then
                  do i = 1, kg1t
                     zaj_l_prev(i,map_z(ib),ik,1) = wfdp_l(i,1)
                  end do
                else
                  do i = 1, kg1t
                     zaj_l(i,map_z(ib),ik,1) = wfdp_l(i,1)
                  end do
                endif
             else if(kimg==2) then
                if(pzaj)then
                  do i = 1, kg1t
                     zaj_l_prev(i,map_z(ib),ik,1) = wfdp_l(i,1)
                     zaj_l_prev(i,map_z(ib),ik,2) = wfdp_l(i,2)
                  end do
                else
                  do i = 1, kg1t
                     zaj_l(i,map_z(ib),ik,1) = wfdp_l(i,1)
                     zaj_l(i,map_z(ib),ik,2) = wfdp_l(i,2)
                  end do
                endif
             end if

           end if
          end do
                                                  __TIMER_IODO_STOP(1435)
       end do
    else
#ifdef _DEBUG_ESIO_
       if(mype == 0) write(nfout,'("### zaj reading")')
#endif
       do ispin = 1, nspin, af+1
       !do ik = 1, kv3, af+1
       do ik = ispin, kv3-nspin+ispin, nspin
                                                  __TIMER_IODO_START(1437)
          do ib = 1, neg_previous
                                                  __TIMER_IODO_START(1438)
! -----------------
           if(precision_WFfile==SP) then
             if(mype == 0) read(nfzaj, end = 9999, err = 9999) wf_l  ! MPI
#ifdef _DEBUG_ESIO_
             if(mype == 0) then
                write(nfout,'(" ik = ",i3, " ib = ",i4)')  ik, ib
                write(nfout,'(8f8.4)') (wf_l(ri,1),ri=1,8)
             end if
#endif
                                                  __TIMER_IODO_STOP(1438)
                                                  __TIMER_IOCOMM_START_w_BARRIER(mpi_comm,1439)
                                                  __TIMER_IOCOMM_STOP(1439)
             if(mype == 0 .and. map_eks(ib,ik,ispin) /= 0) then ! MPI
                call mpi_send(wf_l,kg1t*kimg,mpi_real,map_eks(ib,ik,ispin),1,MPI_CommGroup,ierr) ! MPI
             else if(map_eks(ib,ik,ispin) == mype .and. map_eks(ib,ik,ispin) /= 0) then                  ! MPI
                call mpi_recv(wf_l,kg1t*kimg,mpi_real,0,1,MPI_CommGroup,istatus,ierr)     ! MPI
             end if
             if(map_eks(ib,ik,ispin) == mype) then              ! MPI
                if(pzaj)then
                  do ri = 1, kimg
                     zaj_l_prev(1:kg1t,map_z(ib),ik,ri) = wf_l(1:kg1t,ri)  ! MPI
                  end do
                else
                  do ri = 1, kimg
                     zaj_l(1:kg1t,map_z(ib),ik,ri) = wf_l(1:kg1t,ri)  ! MPI
                  end do
                endif
             end if

! -----------------
           else if(precision_WFfile==DP) then
             if(mype == 0) read(nfzaj, end = 9999, err = 9999) wfdp_l  ! MPI
#ifdef _DEBUG_ESIO_
             if(mype == 0) then
                write(nfout,'(" ik = ",i3, " ib = ",i4)')  ik, ib
                write(nfout,'(8f8.4)') (wfdp_l(ri,1),ri=1,8)
             end if
             call flush(nfout)
#endif
                                                  __TIMER_IODO_STOP(1438)
                                                  __TIMER_IOCOMM_START_w_BARRIER(mpi_comm,1439)
                                                  __TIMER_IOCOMM_STOP(1439)
             if(mype == 0 .and. map_eks(ib,ik,ispin) /= 0) then ! MPI
                call mpi_send(wfdp_l,kg1t*kimg,mpi_double_precision,map_eks(ib,ik,ispin),1,MPI_CommGroup,ierr) ! MPI
             else if(map_eks(ib,ik,ispin) == mype .and. map_eks(ib,ik,ispin) /= 0) then                  ! MPI
                call mpi_recv(wfdp_l,kg1t*kimg,mpi_double_precision,0,1,MPI_CommGroup,istatus,ierr)     ! MPI
             end if
             if(map_eks(ib,ik,ispin) == mype) then              ! MPI
                if(pzaj)then
                  do ri = 1, kimg
                     zaj_l_prev(1:kg1t,map_z(ib),ik,ri) = wfdp_l(1:kg1t,ri)  ! MPI
                  end do
                else
                  do ri = 1, kimg
                     zaj_l(1:kg1t,map_z(ib),ik,ri) = wfdp_l(1:kg1t,ri)  ! MPI
                  end do
                endif
             end if

           endif
! -----------------
          end do
                                                  __TIMER_IODO_STOP(1437)
       end do
       end do
    end if

    if(precision_WFfile==SP) then
       deallocate(wf_l)
    else if(precision_WFfile==DP) then
       deallocate(wfdp_l)
    end if
    call tstatc0_end(id_sname)
    return
9999 continue
    ierror = EOF_REACHED
    call phase_error_wo_filename(ierror, nfout, nfzaj, __LINE__, __FILE__)
                                                  __TIMER_SUB_STOP(1372)
  end subroutine m_ESIO_rd_WFs

! ==== EXP_CELLOPT ==== 2015/09/24
  subroutine m_ESIO_import_WFs_prev_cell(nfout,nfzaj, F_ZAJ_partitioned)
    integer, intent(in) :: nfout, nfzaj
    logical, intent(in) :: F_ZAJ_partitioned
    integer    :: ik,ib,ri, i
    integer    :: id_sname = -1
    integer    :: ierror

    call tstatc0_begin('m_ESIO_import_WFs_prev_cell ',id_sname)

    if(precision_WFfile==SP) then
       if(ipri >= 1) write(nfout,*) ' !D Reading zaj (single_precision)'
    else
       if(ipri >= 1) write(nfout,*) ' !D Reading zaj (double_precision)'
    end if
    if(precision_WFfile==SP) then
       allocate(wf_l(kg1_prev,kimg)); wf_l = 0.d0
    else
       allocate(wfdp_l(kg1_prev,kimg)); wfdp_l = 0.d0
    end if
    rewind nfzaj

    zaj_l = 0.0d0

    if(F_ZAJ_partitioned) then
       do ik = ista_k, iend_k, af+1        ! MPI

          do ib = ista_e, iend_e, istep_e  ! MPI
             if(ib > neg_previous) cycle

             if(precision_WFfile==SP) then
                read(nfzaj) wf_l

                if(kimg == 1) then
                   do i = 1, min( kg1, kg1_prev )
                      zaj_l(i,map_z(ib),ik,1) = wf_l(i,1)
                   end do
                else if(kimg==2) then
                   do i = 1, min( kg1, kg1_prev )
                      zaj_l(i,map_z(ib),ik,1) = wf_l(i,1)
                      zaj_l(i,map_z(ib),ik,2) = wf_l(i,2)
                   end do
                end if
             else if(precision_WFfile==DP) then
                read(nfzaj) wfdp_l

                if(kimg == 1) then
                   do i = 1, min( kg1, kg1_prev )
                      zaj_l(i,map_z(ib),ik,1) = wfdp_l(i,1)
                   end do
                else if(kimg==2) then
                   do i = 1, min( kg1, kg1_prev )
                      zaj_l(i,map_z(ib),ik,1) = wfdp_l(i,1)
                      zaj_l(i,map_z(ib),ik,2) = wfdp_l(i,2)
                   end do
                end if

             end if
          end do
       end do
    else

       do ik = 1, kv3, af+1
          do ib = 1, neg_previous
             ! -----------------
             if(precision_WFfile==SP) then
                if(mype == 0) read(nfzaj, end = 9999, err = 9999) wf_l
                if(mype == 0 .and. map_ek(ib,ik) /= 0) then
                   call mpi_send(wf_l,kg1_prev*kimg,mpi_real,map_ek(ib,ik),1,MPI_CommGroup,ierr) ! MPI
                else if(map_ek(ib,ik) == mype .and. map_ek(ib,ik) /= 0) then
                   call mpi_recv(wf_l,kg1_prev*kimg,mpi_real,0,1,MPI_CommGroup,istatus,ierr)     ! MPI
                end if
                if(map_ek(ib,ik) == mype) then              ! MPI
                   do ri = 1, kimg
                      do i = 1, min( kg1, kg1_prev )
                         zaj_l(i,map_z(ib),ik,ri) = wf_l(i,ri)  ! MPI
                      end do
                   end do
                end if

                ! -----------------
             else if(precision_WFfile==DP) then
                if(mype == 0) read(nfzaj, end = 9999, err = 9999) wfdp_l
                if(mype == 0 .and. map_ek(ib,ik) /= 0) then
                   call mpi_send(wfdp_l,kg1_prev*kimg,mpi_double_precision,map_ek(ib,ik),1,MPI_CommGroup,ierr) ! MPI
                else if(map_ek(ib,ik) == mype .and. map_ek(ib,ik) /= 0) then
                   call mpi_recv(wfdp_l,kg1_prev*kimg,mpi_double_precision,0,1,MPI_CommGroup,istatus,ierr)     ! MPI
                end if
                if(map_ek(ib,ik) == mype) then              ! MPI
                   do i = 1, min( kg1, kg1_prev )
                      zaj_l(i,map_z(ib),ik,ri) = wfdp_l(i,ri)  ! MPI
                   end do
                end if
             endif
             ! -----------------
          end do
       end do
    end if
!
    if(precision_WFfile==SP) then
       deallocate(wf_l)
    else if(precision_WFfile==DP) then
       deallocate(wfdp_l)
    end if
    call tstatc0_end(id_sname)
    return
9999 continue
    ierror = EOF_REACHED
    call phase_error_wo_filename(ierror, nfout, nfzaj, __LINE__, __FILE__)
  end subroutine m_ESIO_import_WFs_prev_cell
! ===================== 2015/09/24

! ==================================== added by K. Tagami =============== 11.0
  subroutine m_ESIO_rd_WFs_import_frm_collin(nfout,nfzaj, F_ZAJ_partitioned)
    use m_Crystal_Structure,  only : sw_fix_global_quantz_axis, &
         &                           Global_Quantz_Axis_Fixed
    use m_Const_Parameters,   only : zi

    integer, intent(in) :: nfout, nfzaj
    logical, intent(in) :: F_ZAJ_partitioned

    integer :: neg_to_be_read
    integer    :: id_sname = -1
    call tstatc0_begin('m_ESIO_rd_WFs_import_frm_collin ',id_sname)

    if(precision_WFfile==SP) then
       allocate(wf_l(kg1,kimg))   ;       wf_l = 0
    else
       allocate(wfdp_l(kg1,kimg)) ;       wfdp_l = 0
    end if

    rewind nfzaj
    if(ipri >= 1) write(nfout,*) ' !D Reading zaj'

    zaj_l = 0.0d0

    if(F_ZAJ_partitioned) then
       write(*,*) &
            & 'Not supported : importing collinear Wfns  when F_ZAJ_partitioned = true'

    else
       neg_to_be_read = neg / 2

       write(nfout,*) '******************************** '
       write(nfout,*) '!! Collinear wavefunctions are used. '
       write(nfout,*) '!! neg_to_be_read is assumed to be ', neg_to_be_read
       write(nfout,*) '******************************** '

       if ( previous_nspin_collinear == 1 ) then
          call case_previous_nspin_eq_1
       else if ( previous_nspin_collinear == 2 ) then
          call case_previous_nspin_eq_2
       endif

    end if

    if(precision_WFfile==SP) then
       deallocate(wf_l)
    else if(precision_WFfile==DP) then
       deallocate(wfdp_l)
    end if
    call tstatc0_end(id_sname)

  contains

    subroutine case_previous_nspin_eq_2
      integer :: ik, ib, ri, i
      integer :: ib_0, is, ig
      real(kind=DP) :: spn_quant_dir(3), theta, phi
      complex(kind=CMPLDP) :: rotated_spinor(2,2), z1, z2, ztmp1, ztmp2

      do ik = 1, kv3, ndim_spinor
         Do is=1, ndim_spinor
            do ib_0 = 1, previous_nband_collinear
               ib = ( ib_0 -1 )*ndim_spinor + is
               if(precision_WFfile==SP) then
                  if(mype == 0) read(nfzaj) wf_l

                  if(mype == 0 .and. map_ek(ib,ik) /= 0) then
                     call mpi_send( wf_l, kg1*kimg, mpi_real, map_ek(ib,ik), 1, &
                          &         MPI_CommGroup, ierr )
                  else if(map_ek(ib,ik) == mype .and. map_ek(ib,ik) /= 0 ) then
                     call mpi_recv( wf_l, kg1*kimg, mpi_real, 0, 1, &
                          &         MPI_CommGroup, istatus, ierr )
                  end if

                  if ( ib_0 > neg_to_be_read ) cycle

                  if(map_ek(ib,ik) == mype) then              ! MPI
                     do ri = 1, kimg
                        zaj_l(1:kg1,map_z(ib),ik+is-1,ri) = wf_l(1:kg1,ri)  ! MPI
                     end do
                  end if
               else if(precision_WFfile==DP) then
                  if(mype == 0) read(nfzaj) wfdp_l

                  if(mype == 0 .and. map_ek(ib,ik) /= 0) then
                     call mpi_send( wfdp_l, kg1*kimg, mpi_double_precision, map_ek(ib,ik), 1, MPI_CommGroup, ierr )
                  else if(map_ek(ib,ik) == mype .and. map_ek(ib,ik) /= 0 ) then
                     call mpi_recv( wfdp_l, kg1*kimg, mpi_double_precision, 0, 1, MPI_CommGroup, istatus, ierr )
                  end if

                  if ( ib_0 > neg_to_be_read ) cycle

                  if(map_ek(ib,ik) == mype) then              ! MPI
                     do ri = 1, kimg
                        zaj_l(1:kg1,map_z(ib),ik+is-1,ri) = wfdp_l(1:kg1,ri)  ! MPI
                     end do
                  end if
               end if
            end do
         End Do
      end do
      if ( sw_fix_global_quantz_axis == ON ) then
         spn_quant_dir(1:3) = Global_Quantz_Axis_Fixed(1:3)
         theta = acos( spn_quant_dir(3) )
         phi = atan2( spn_quant_dir(2), spn_quant_dir(1) )

         rotated_spinor(1,1) =  exp( -zi *phi /2.0d0 ) *cos( theta /2.0d0 )
         rotated_spinor(2,1) =  exp(  zi *phi /2.0d0 ) *sin( theta /2.0d0 )

         rotated_spinor(1,2) = -exp( -zi *phi /2.0d0 ) *sin( theta /2.0d0 )
         rotated_spinor(2,2) =  exp(  zi *phi /2.0d0 ) *cos( theta /2.0d0 )

         do ik = 1, kv3, ndim_spinor
            if ( map_k(ik) /= myrank_k ) cycle
            Do ig=1, kg1
               Do ib=1, np_e
                  z1 = dcmplx( zaj_l(ig,ib,ik,  1), zaj_l(ig,ib,ik,  kimg) )
                  z2 = dcmplx( zaj_l(ig,ib,ik+1,1), zaj_l(ig,ib,ik+1,kimg) )
                  ztmp1 = rotated_spinor(1,1) *z1 +rotated_spinor(2,1) *z2
                  ztmp2 = rotated_spinor(1,2) *z1 +rotated_spinor(2,2) *z2
                  zaj_l(ig,ib,ik,1)      = real(ztmp1)
                  zaj_l(ig,ib,ik,kimg)   = aimag(ztmp1)
                  zaj_l(ig,ib,ik+1,1)    = real(ztmp2)
                  zaj_l(ig,ib,ik+1,kimg) = aimag(ztmp2)
               End Do
            End Do
         end do
      end if
    end subroutine case_previous_nspin_eq_2

    subroutine case_previous_nspin_eq_1
      integer :: ik, ib, ri, i
      integer :: ib_0, is

      do ik = 1, kv3, ndim_spinor
         do ib_0 = 1, previous_nband_collinear
            if(precision_WFfile==SP) then
               if(mype == 0) read(nfzaj) wf_l              ! MPI

               Do is=1, ndim_spinor
                  ib = ( ib_0 -1 )*ndim_spinor + is
                  if(mype == 0 .and. map_ek(ib,ik) /= 0) then ! MPI
                     call mpi_send( wf_l, kg1*kimg, mpi_real, map_ek(ib,ik), 1, &
                          &         MPI_CommGroup, ierr ) ! MPI
                  else if(map_ek(ib,ik) == mype .and. map_ek(ib,ik) /= 0 ) then     ! MPI
                     call mpi_recv( wf_l, kg1*kimg, mpi_real, 0, 1, &
                          &         MPI_CommGroup, istatus, ierr )     ! MPI
                  end if

                  if ( ib_0 > neg_to_be_read ) cycle

                  if(map_ek(ib,ik) == mype) then              ! MPI
                     do ri = 1, kimg
                        zaj_l(1:kg1,map_z(ib),ik+is-1,ri) = wf_l(1:kg1,ri)  ! MPI
                     end do
                  end if

               end do
            else
               Do is=1, ndim_spinor
                  ib = ( ib_0 -1 )*ndim_spinor + is
                  if(mype == 0 .and. map_ek(ib,ik) /= 0) then
                     call mpi_send( wfdp_l, kg1*kimg, mpi_double_precision, map_ek(ib,ik), 1, MPI_CommGroup, ierr )
                  else if(map_ek(ib,ik) == mype .and. map_ek(ib,ik) /= 0 ) then
                     call mpi_recv( wfdp_l, kg1*kimg, mpi_double_precision, 0, 1, MPI_CommGroup, istatus, ierr )
                  end if

                  if ( ib_0 > neg_to_be_read ) cycle

                  if(map_ek(ib,ik) == mype) then              ! MPI
                     do ri = 1, kimg
                        zaj_l(1:kg1,map_z(ib),ik+is-1,ri) = wfdp_l(1:kg1,ri)  ! MPI
                     end do
                  end if

               end do
            end if
         End Do
      end do
    end subroutine case_previous_nspin_eq_1

  end subroutine m_ESIO_rd_WFs_import_frm_collin
!===================================================================== 11.0

  subroutine m_ESIO_wd_WFs(nfout,nfzaj,F_ZAJ_partitioned,rew)
    integer, intent(in) :: nfout,nfzaj
    logical, intent(in) :: F_ZAJ_partitioned
    logical, intent(in), optional :: rew
    integer :: ik,ib,ri, ispin
    integer :: id_sname = -1
    call tstatc0_begin('m_ESIO_wd_WFs ',id_sname)

! -------------------- T. Hamada 2021.9.22 --------------------
   if(precision_WFfile == 0) return                  ! UVSOR
! -------------------------------------------------------------
   if(precision_WFfile==SP) then
    allocate(wf_l(kg1,kimg));    wf_l = 0
    call mpi_barrier(MPI_CommGroup,ierr)
    !!$ print *, ' !D Writing zaj '
    if(ipri >= 1) write(nfout,*) ' !D Writing zaj (single_precision)'
    if (present(rew)) then
      if(rew) then
        rewind nfzaj
      endif
    else
      rewind nfzaj
    endif
    if(F_ZAJ_partitioned) then
       do ik = ista_k, iend_k, af+1        ! MPI
          do ib = ista_e, iend_e, istep_e  ! MPI
             do ri = 1, kimg
                wf_l(1:kg1,ri) = zaj_l(1:kg1,map_z(ib),ik,ri)
             end do
             write(nfzaj) wf_l
          end do
       end do
    else
#ifdef _DEBUG_ESIO_
       if(mype == 0) write(nfout,'("### zaj writing")')
#endif
       do ispin = 1, nspin, af+1
       !do ik = 1, kv3, af+1
       do ik = ispin, kv3-nspin+ispin, nspin
          do ib = 1, neg
             !if(map_ek(ib,ik) == mype) then                          ! MPI
             if(map_eks(ib,ik,ispin) == mype) then                          ! MPI
                do ri = 1, kimg
                   wf_l(1:kg1,ri) = zaj_l(1:kg1,map_z(ib),ik,ri)
                end do
                if(map_eks(ib,ik,ispin) /= 0) &                             ! MPI
                     &   call mpi_send(wf_l,kg1*kimg,mpi_real,0,1,MPI_CommGroup,ierr) ! MPI
             else if(mype == 0 .and. map_eks(ib,ik,ispin) /= 0) then
                call mpi_recv(wf_l,kg1*kimg,mpi_real,map_eks(ib,ik,ispin),1,MPI_CommGroup,istatus,ierr)!MPI
             end if
             if(mype == 0) write(nfzaj)  wf_l                        ! MPI
#ifdef _DEBUG_ESIO_
             if(mype == 0) then
                write(nfout,'(" ik = ",i4, " ib = ",i5)')  ik, ib
                write(nfout,'(8f8.4)') (wf_l(ri,1),ri=1,8)
             end if
#endif
          end do
       end do
       end do
    end if
    deallocate(wf_l)
   else if(precision_WFfile==DP) then
    allocate(wfdp_l(kg1,kimg));    wfdp_l = 0
    call mpi_barrier(MPI_CommGroup,ierr)
    if(ipri >= 1) write(nfout,*) ' !D Writing zaj (double_precision) '
    if (present(rew)) then
      if(rew) then
        rewind nfzaj
      endif
    else
      rewind nfzaj
    endif
    if(F_ZAJ_partitioned) then
       do ik = ista_k, iend_k, af+1
          do ib = ista_e, iend_e, istep_e
             do ri = 1, kimg
                wfdp_l(1:kg1,ri) = zaj_l(1:kg1,map_z(ib),ik,ri)
             end do
             write(nfzaj) wfdp_l
          end do
       end do
    else
#ifdef _DEBUG_ESIO_
       if(mype == 0) write(nfout,'("### zaj writing")')
#endif
       do ispin = 1, nspin, af+1
       !do ik = 1, kv3, af+1
       do ik = ispin, kv3-nspin+ispin, nspin
          do ib = 1, neg
             if(map_eks(ib,ik,ispin) == mype) then
                do ri = 1, kimg
                   wfdp_l(1:kg1,ri) = zaj_l(1:kg1,map_z(ib),ik,ri)
                end do
                if(map_eks(ib,ik,ispin) /= 0) &
                     &   call mpi_send(wfdp_l,kg1*kimg,mpi_double_precision,0,1,MPI_CommGroup,ierr)
             else if(mype == 0 .and. map_eks(ib,ik,ispin) /= 0) then
                call mpi_recv(wfdp_l,kg1*kimg,mpi_double_precision,map_eks(ib,ik,ispin),1,MPI_CommGroup,istatus,ierr)
             end if
             if(mype == 0) write(nfzaj)  wfdp_l
#ifdef _DEBUG_ESIO_
             if(mype == 0) then
                write(nfout,'(" ik = ",i4, " ib = ",i5)')  ik, ib
                write(nfout,'(8f8.4)') (wfdp_l(ri,1),ri=1,8)
             end if
#endif
          end do
       end do
       end do
    end if
    deallocate(wfdp_l)
   end if
    call tstatc0_end(id_sname)
  end subroutine m_ESIO_wd_WFs

  subroutine m_ESIO_wd_WFs_standardout(nfout,ipriwf)
    integer, intent(in) :: nfout,ipriwf
    integer :: ik,ib,ri, i, ic, ipriwf0, icycle, icolumn, max_elements, istart, iend
    integer :: id_sname = -1
    real(kind=DP) :: phase2r, phase2i, phaser,phasei
    complex(kind=CMPLDP) :: exp2theta, exptheta
    call tstatc0_begin('m_ESIO_wd_WFs_stndout ',id_sname)

    ipriwf0 = ipriwf
    if(npes > 1) call mpi_bcast(ipriwf0,1,mpi_integer,0,MPI_CommGroup,ierr)

    if(ipriwf0 >= 2) then
       icolumn = 10
       if(precision_WFfile==SP) then
          allocate(wf_l(kg1,kimg+3)); wf_l = 0.d0
       else
          allocate(wfdp_l(kg1,kimg+3)); wfdp_l = 0.d0
       end if
       call mpi_barrier(MPI_CommGroup,ierr)
       if(mype == 0)  write(nfout,*) ' !wf Writing zaj '

       do ik = 1, kv3, af+1
          max_elements = iba(ik)
          if(mype == 0) write(nfout,'(" !wf   ik = ",i5)') ik
          do ib = 1, neg
             if(mype == 0) write(nfout,'(" !wf   ib = ",i5)') ib
           if(precision_WFfile==SP) then
             if(map_ek(ib,ik) == mype) then                          ! MPI
                do ri = 1, kimg
                   wf_l(1:kg1,ri) = zaj_l(1:kg1,map_z(ib),ik,ri)
                end do
                if(map_ek(ib,ik) /= 0) &                             ! MPI
                     &   call mpi_send(wf_l,kg1*kimg,mpi_real,0,1,MPI_CommGroup,ierr) ! MPI
             else if(mype == 0 .and. map_ek(ib,ik) /= 0) then
                call mpi_recv(wf_l,kg1*kimg,mpi_real,map_ek(ib,ik),1,MPI_CommGroup,istatus,ierr)!MPI
             end if
           else
             if(map_ek(ib,ik) == mype) then                          ! MPI
                do ri = 1, kimg
                   wfdp_l(1:kg1,ri) = zaj_l(1:kg1,map_z(ib),ik,ri)
                end do
                if(map_ek(ib,ik) /= 0) &                             ! MPI
                     &   call mpi_send(wfdp_l,kg1*kimg,mpi_double_precision,0,1,MPI_CommGroup,ierr) ! MPI
             else if(mype == 0 .and. map_ek(ib,ik) /= 0) then
                call mpi_recv(wfdp_l,kg1*kimg,mpi_double_precision,map_ek(ib,ik),1,MPI_CommGroup,istatus,ierr)!MPI
             end if
           end if
             if(mype == 0) then
                if(kimg == 2) then
                 if(precision_WFFile==SP) then
                   do i = 1, iba(ik)
                      wf_l(i,3) = wf_l(i,1)**2 + wf_l(i,2)**2
                   end do
                 else
                   do i = 1, iba(ik)
                      wfdp_l(i,3) = wfdp_l(i,1)**2 + wfdp_l(i,2)**2
                   end do
                 end if
                   if(k_symmetry(ik) == GAMMA .or. k_symmetry(ik) == GAMMA_base_symmetrization) then
                    if(precision_WFfile==SP) then
                      phase2r = (wf_l(1,1)**2 - wf_l(1,2)**2)/wf_l(1,3)
                      phase2i = -2.d0*wf_l(1,1)*wf_l(1,2)/wf_l(1,3)
                    else
                      phase2r = (wfdp_l(1,1)**2 - wfdp_l(1,2)**2)/wfdp_l(1,3)
                      phase2i = -2.d0*wfdp_l(1,1)*wfdp_l(1,2)/wfdp_l(1,3)
                    end if
                      exp2theta = cmplx(phase2r, phase2i)
                      exptheta = sqrt(exp2theta)
                      phaser = real(exptheta)
                      phasei = imag(exptheta)
                      write(nfout,'(" !wf exp2theta = ",2d20.8)') exp2theta
                      write(nfout,'(" !wf           = ",2d20.8)') phase2r, phase2i
                      write(nfout,'(" !wf |exp2theta|**2 = ",d20.8)') dsqrt(phase2r**2 + phase2i**2)
                      write(nfout,'(" !wf exptheta  = ",2d20.8)') exptheta
                      write(nfout,'(" !wf |exptheta| = ",d20.8)') abs(exptheta)
                     if(precision_WFfile==SP) then
                      do i = 1, iba(ik)
!!$                      wf_l(i,4) = real(exptheta*cmplx(wf_l(i,1),wf_l(i,2)))
!!$                      wf_l(i,5) = imag(exptheta*cmplx(wf_l(i,1),wf_l(i,2)))
                         wf_l(i,4) = phaser*wf_l(i,1) - phasei*wf_l(i,2)
                         wf_l(i,5) = phaser*wf_l(i,2) + phasei*wf_l(i,1)
                      end do
                     else
                      do i = 1, iba(ik)
                         wfdp_l(i,4) = phaser*wfdp_l(i,1) - phasei*wfdp_l(i,2)
                         wfdp_l(i,5) = phaser*wfdp_l(i,2) + phasei*wfdp_l(i,1)
                      end do
                     end if
                   end if
                else
                 if(precision_WFfile==SP) then
                   do i = 1, iba(ik)
                      wf_l(i,2) = wf_l(i,1)**2
                   end do
                 else
                   do i = 1, iba(ik)
                      wfdp_l(i,2) = wfdp_l(i,1)**2
                   end do
                 end if
                   exp2theta = 1.d0
                   exptheta = 1.d0
                end if
                icycle = ceiling(dble(min(max_elements,iba(ik)))/icolumn)
                istart = 1
                do ic = 1, icycle
                   iend = min(istart+icolumn-1,max_elements,iba(ik))
                   write(nfout,'(" !wf (nx)    ",10i10)') (ngabc(nbase(i,ik),1),i=istart,iend)
                   write(nfout,'(" !wf (ny)    ",10i10)') (ngabc(nbase(i,ik),2),i=istart,iend)
                   write(nfout,'(" !wf (nz)    ",10i10)') (ngabc(nbase(i,ik),3),i=istart,iend)
                 if(precision_WFfile==SP) then
                   write(nfout,'(" !wf (zaj-r) ",10d10.2)') (wf_l(i,1),i=istart,iend)
                   if(kimg == 2) write(nfout,'(" !wf (zaj-i) ",10d10.2)') (wf_l(i,2),i=istart,iend)
                   write(nfout,'(" !wf abs     ",10d10.2)') (wf_l(i,kimg+1),i=istart,iend)
                   if(kimg == 2) then
                      if(k_symmetry(ik) == GAMMA .or. k_symmetry(ik) == GAMMA_base_symmetrization) then
                         write(nfout,'(" !wf (zaj-r)d",10d10.2)') (wf_l(i,4),i=istart,iend)
                         write(nfout,'(" !wf (zaj-i)d",10d10.2)') (wf_l(i,5),i=istart,iend)
                      end if
                   end if
                 else
                   write(nfout,'(" !wf (zaj-r) ",10d10.2)') (wfdp_l(i,1),i=istart,iend)
                   if(kimg == 2) write(nfout,'(" !wf (zaj-i) ",10d10.2)') (wfdp_l(i,2),i=istart,iend)
                   write(nfout,'(" !wf abs     ",10d10.2)') (wfdp_l(i,kimg+1),i=istart,iend)
                   if(kimg == 2) then
                      if(k_symmetry(ik) == GAMMA .or. k_symmetry(ik) == GAMMA_base_symmetrization) then
                         write(nfout,'(" !wf (zaj-r)d",10d10.2)') (wfdp_l(i,4),i=istart,iend)
                         write(nfout,'(" !wf (zaj-i)d",10d10.2)') (wfdp_l(i,5),i=istart,iend)
                      end if
                   end if
                 end if
                   istart = iend+1
                end do
             end if
          end do
       end do
      if(precision_WFfile==SP) then
       deallocate(wf_l)
      else
       deallocate(wfdp_l)
      end if
    end if
    call tstatc0_end(id_sname)
  end subroutine m_ESIO_wd_WFs_standardout

  subroutine m_ESIO_rd_next_wfs_ek(ik,nfout,nfzaj)
    integer, intent(in) :: ik
    integer, intent(in) :: nfout,nfzaj
    integer :: ib, ri, j, ierror
    if(precision_WFfile==SP) then
       allocate(wf_l(kg1,kimg)); wf_l = 0_sp
    else
       allocate(wfdp_l(kg1,kimg)); wfdp_l = 0.d0
    end if
!    zaj_l(:,:,ik,:) = 0.d0
    do ib = 1, neg
! -----------------
       if(precision_WFfile==SP) then
          if(mype == 0) read(nfzaj, end = 9999, err = 9999) wf_l  ! MPI
          if(mype == 0 .and. map_ek(ib,ik) /= 0) then ! MPI
             call mpi_send(wf_l,kg1*kimg,mpi_real4,map_ek(ib,ik),1,MPI_CommGroup,ierr) ! MPI
          else if(map_ek(ib,ik) == mype .and. map_ek(ib,ik) /= 0) then                  ! MPI
             call mpi_recv(wf_l,kg1*kimg,mpi_real4,0,1,MPI_CommGroup,istatus,ierr)     ! MPI
          end if
          if(map_ek(ib,ik) == mype) then              ! MPI
            do ri = 1, kimg
               zaj_l(1:kg1,map_z(ib),ik,ri) = wf_l(1:kg1,ri)  ! MPI
            end do
          end if

! -----------------
       else if(precision_WFfile==DP) then
          if(mype == 0) read(nfzaj, end = 9999, err = 9999) wfdp_l  ! MPI
          if(mype == 0 .and. map_ek(ib,ik) /= 0) then ! MPI
             call mpi_send(wfdp_l,kg1*kimg,mpi_double_precision,map_ek(ib,ik),1,MPI_CommGroup,ierr) ! MPI
          else if(map_ek(ib,ik) == mype .and. map_ek(ib,ik) /= 0) then                  ! MPI
             call mpi_recv(wfdp_l,kg1*kimg,mpi_double_precision,0,1,MPI_CommGroup,istatus,ierr)     ! MPI
          end if
          if(map_ek(ib,ik) == mype) then              ! MPI
            do ri = 1, kimg
               zaj_l(1:kg1,map_z(ib),ik,ri) = wfdp_l(1:kg1,ri)  ! MPI
            end do
          end if
       endif
    end do

    if(precision_WFfile==SP) then
       deallocate(wf_l)
    else
       deallocate(wfdp_l)
    end if

    return

9999 continue
    ierror = EOF_REACHED
    call phase_error_wo_filename(ierror, nfout, nfzaj, __LINE__, __FILE__)
  end subroutine m_ESIO_rd_next_wfs_ek

  subroutine m_ESIO_rd_WFs_and_EVs_ek(nfout,nf)
    integer, intent(in) :: nfout,nf
    integer  :: ik, ie, iks, ri, ikg, ikt, ike, ike2
    integer, allocatable, dimension(:,:) :: n_mpi, n2_mpi  ! MPI
    real(DP),allocatable, dimension(:,:) :: e_mpi, e2_mpi  ! MPI
!!$    read(nf) neordr,nrvf_ordr,eko_l,occup_l,efermi,totch

    allocate(n_mpi(neg,nspin)); allocate(n2_mpi(neg,nspin)) ! MPI
    allocate(e_mpi(neg,nspin)); allocate(e2_mpi(neg,nspin)) ! MPI
    allocate(wf_l(kg1,kimg))

    n_mpi =0; n2_mpi = 0
    e_mpi =0; e2_mpi = 0; wf_l = 0
    if(ipri >= 1) write(nfout,*) ' !D Reading zaj'

    rewind nf

    eko_ek = 0.d0
!!$    do ik = 1, nk_in_the_process-kv3, kv3
    do ik = 1, nk_in_the_process - nspin, nspin
       if(ipri >= 1) write(nfout,*) ' !D     skipping ik = ', ik
!!$       do iks = 1, kv3, af+1
       do iks = 1, nspin, af+1
          do ie = 1, neg
             if(mype == 0) read(nf) wf_l
          end do
       end do
       if(mype == 0) read(nf) n_mpi
       if(mype == 0) read(nf) n2_mpi
       if(mype == 0) read(nf) e_mpi

       if(mype == 0) then
!!$          do iks = 1, kv3
          do iks = 1, nspin
             do ie = 1, neg
                eko_ek(ie,ik+iks-1) = e_mpi(n_mpi(ie,iks),iks)
             end do
          end do
          if(ipri >= 3) then
!!$             do iks=1,kv3
             do iks=1,nspin
                write(nfout,'(" ik = ",i5)') ik+iks-1
                write(nfout,'(8f8.4)') (e_mpi(n_mpi(ie,iks),iks),ie=1,neg)
             end do
          end if
       end if
    end do

    if(nk_in_the_process > kv3_ek) goto 1001

!!$    do ik = 1, kv3, af+1
    KPOINT_LOOP: do ikg = 1, nrank_k
       do ikt = 1, nspin, af+1
          ik = (ikg-1)*nspin+ikt
          if(nk_in_the_process -1 + ik > kv3_ek) exit KPOINT_LOOP
          if(nk_in_the_process -1 + ik > numk_zajsaved) exit KPOINT_LOOP
          if(ipri>=1) write(nfout,*) ' !D     reading  ik = ', ik+first_kpoint_in_this_job-1

          do ie = 1, neg
             if(ipri>=2) write(nfout,*) ' !D     ie = ', ie
             if(mype == 0) read(nf,err=2,end=2) wf_l     ! MPI
             if(mype == 0 .and. map_ek(ie,ik) /= 0) then ! MPI
                call mpi_send(wf_l,kg1*kimg,mpi_real4,map_ek(ie,ik),1,MPI_CommGroup,ierr) ! MPI
             else if(map_ek(ie,ik) == mype .and. map_ek(ie,ik) /= 0) then                  ! MPI
                call mpi_recv(wf_l,kg1*kimg,mpi_real4,0,1,MPI_CommGroup,istatus,ierr)     ! MPI
             end if
             if(map_ek(ie,ik) == mype) then              ! MPI
                do ri = 1, kimg
                   zaj_l(1:kg1,map_z(ie),ik,ri) = wf_l(1:kg1,ri)  ! MPI
                end do
             end if
          end do
       end do

       if(mype == 0) read(nf) n_mpi                ! MPI
       if(mype == 0) read(nf) n2_mpi               ! MPI
       call mpi_bcast(n_mpi,neg*nspin,mpi_integer,0,MPI_CommGroup,ierr) ! MPI
       call mpi_bcast(n2_mpi,neg*nspin,mpi_integer,0,MPI_CommGroup,ierr)! MPI

       do ikt = 1, nspin                             ! MPI
          ik = (ikg-1)*nspin + ikt
          if(map_k(ik) == myrank_k) then
             neordr(1:neg,ik) = n_mpi(1:neg,ikt)     ! MPI
             nrvf_ordr(1:neg,ik) = n2_mpi(1:neg,ikt) ! MPI
          end if
       end do                                        ! MPI

       if(mype == 0) read(nf) e_mpi                ! MPI
       call mpi_bcast(e_mpi,neg*nspin,mpi_double_precision,0,MPI_CommGroup,ierr) ! MPI

       do ikt = 1, nspin                                   ! MPI
          ik = (ikg-1)*nspin + ikt
          do ie = 1, neg                                   ! MPI
             if(map_ek(ie,ik) == mype) then                ! MPI
                eko_l(map_z(ie),ik) = e_mpi(ie,ikt)        ! MPI
             end if                                        ! MPI
          end do                                           ! MPI
          do ie = 1, neg
             eko_ek(ie,nk_in_the_process+ik-1) = e_mpi(n_mpi(ie,ikt),ikt)
          end do
       end do                                              ! MPI
       goto 3
2      continue
       call phase_error_with_msg(nfout, ' eof from nf <<m_ESIO_rd_WFs_and_EVs_ek>>',__LINE__,__FILE__)
3      continue
    end do KPOINT_LOOP

    if(ipri>=2) write(nfout,*) ' !D     ikg  = ', ikg
!!$    if(ikg <= 1) stop ' ikg <= 1 <<m_ESIO_rd_WFs_and_EVs_ek>>'
    if(ikg < nrank_k .and. ikg > 1 ) then
       do ike = ikg, nrank_k
          if(ipri>=2) write(nfout,*) ' !D     ike  = ', ike
          if(ipri>=2) write(nfout,*) ' !D     zaj_l'
          do ikt = 1, nspin, af+1
             ik   = (ikg-2)*nspin+ikt
             ike2 = (ike-1)*nspin+ikt
             if(ipri>=2) write(nfout,*) ' !D     ik, ike2 = ',ik,ike2
             do ie = 1, neg
                if(map_ek(ie,ik) == mype) then
                   do ri = 1, kimg
                      wf_l(1:kg1,ri) = zaj_l(1:kg1,map_z(ie),ik,ri)
                   end do
                end if

                if(map_ek(ie,ik) == mype) then
                   call mpi_send(wf_l,kg1*kimg,mpi_real4,map_ek(ie,ike2),1,MPI_CommGroup,ierr)
                else if(map_ek(ie,ike2) == mype) then
                   call mpi_recv(wf_l,kg1*kimg,mpi_real4,map_ek(ie,ik),1,MPI_CommGroup,istatus,ierr)
                end if
                if(map_ek(ie,ike2) == mype) then
                   do ri = 1, kimg
                      zaj_l(1:kg1,map_z(ie),ike2,ri) = wf_l(1:kg1,ri)
                   end do
                end if
             end do
          end do

          ! ---> neordr, nrvf_ordr
          do ikt = 1, nspin
             ik = (ikg-2)*nspin+ikt
             if(map_k(ik) == myrank_k) then
                n_mpi(1:neg,ikt) = neordr(1:neg,ik)
                n2_mpi(1:neg,ikt) = nrvf_ordr(1:neg,ik)
             end if
          end do
          call mpi_bcast(n_mpi,neg*nspin,mpi_integer,map_k(ik),MPI_CommGroup,ierr)
          call mpi_bcast(n2_mpi,neg*nspin,mpi_integer,map_k(ik),MPI_CommGroup,ierr)

          if(ipri>=2) write(nfout,*) ' !D     neordr and nrvf_ordr'
          do ikt = 1, nspin
             ike2 = (ike-1)*nspin+ikt
             if(map_k(ike2) == myrank_k) then
                neordr(1:neg,ike2) = n_mpi(1:neg,ikt)     ! MPI
                nrvf_ordr(1:neg,ike2) = n2_mpi(1:neg,ikt) ! MPI
             end if
          end do

          ! ---> eko_l
          if(ipri>=2) write(nfout,*) ' !D     eko_l'
          e_mpi = 0.d0
          do ikt = 1, nspin
             ik = (ikg-2)*nspin+ikt
             if(ipri>=2) write(nfout,'(" !D eko_l ik = ",i6)') ik
             do ie = 1, neg
                if(map_ek(ie,ik) == mype) then
                   e_mpi(ie,ikt) = eko_l(map_z(ie),ik)
                end if
             end do
          end do
          if(npes >= 2) then
             call mpi_allreduce(e_mpi,e2_mpi,neg*nspin,mpi_double_precision &
                  & ,mpi_sum,MPI_CommGroup,ierr)
          else
             e2_mpi  = e_mpi
          end if

          if(ipri>=2) then
             write(nfout,'(" <<m_ESIO_rd_WFs_and_EVs_ek>>")')
             do ikt = 1, nspin
                ik = (ikg-2)*nspin+ikt
                write(nfout,'(" ik = ",i5)') ik
                write(nfout,'(10f8.4)') (e2_mpi(ie,ikt),ie=1,neg)
             end do
          end if

          do ikt = 1, nspin
             ike2 = (ike-1)*nspin+ikt
             if(ipri>=2) write(nfout,'(" !D eko_l ike2 = ",i6)') ike2
             do ie = 1, neg
                if(map_ek(ie,ike2) == mype) then
                   eko_l(map_z(ie),ike2) = e2_mpi(ie,ikt)
                end if
             end do
          end do

       end do
    end if

1001 continue

    if(npes >= 2) &
         & call mpi_bcast(eko_ek,neg*kv3_ek,mpi_double_precision,0,MPI_CommGroup,ierr)

    if(nk_in_the_process > kv3_ek) goto 1002
    call m_ESIO_wd_EigenValues(nfout,2,nooccupation=YES)

1002 continue

!!$    ik = nk_in_the_process
!!$    do iks = 1, kv3
!!$          do ie = 1, neg
!!$             eko_ek(ie,ik+iks-1) = e_mpi(n_mpi(ie,iks),iks)
!!$          end do
!!$       end do
!!$    end if
!!$
    if(ipri >= 1) then
       write(nfout,'(" <<m_ESIO_rd_WFs_and_EVs_ek>>")')
       do iks = 1, kv3_ek
          if(iks > numk_zajsaved) cycle
          write(nfout,'(" ik = ",i5)') iks
          write(nfout,'(10f8.4)') (eko_ek(ie,iks),ie=1,neg)
       end do
    end if


    rewind nf

    do ik = 1, nk_in_the_process-nspin, nspin
       do iks = 1, nspin, af+1
          do ie = 1, neg
             if(mype == 0) read(nf) wf_l                 ! MPI
          end do
       end do
       if(mype == 0) read(nf) n_mpi                   ! MPI
       if(mype == 0) read(nf) n2_mpi                  ! MPI
       if(mype == 0) read(nf) e_mpi                   ! MPI
    end do

    deallocate(wf_l)                                    ! MPI
    deallocate(n_mpi); deallocate(n2_mpi)               ! MPI
    deallocate(e_mpi); deallocate(e2_mpi)               ! MPI

!!$    write(nfout,'(" ---<< m_ESIO_rd_WFs_and_EVs_ek>>---")')
  end subroutine m_ESIO_rd_WFs_and_EVs_ek

  subroutine m_ESIO_rd_EVs_ek(nfout,nf)
    integer, intent(in) :: nfout,nf
    integer  :: ik, ie, iks
    integer, allocatable, dimension(:,:) :: n_mpi, n2_mpi  ! MPI
    real(DP),allocatable, dimension(:,:) :: e_mpi, e2_mpi  ! MPI

    allocate(n_mpi(neg,nspin)); allocate(n2_mpi(neg,nspin)) ! MPI
    allocate(e_mpi(neg,nspin)); allocate(e2_mpi(neg,nspin)) ! MPI
    allocate(wf_l(kg1,kimg))
    n_mpi =0; n2_mpi = 0
    e_mpi =0; e2_mpi = 0; wf_l = 0

    if(ipri >= 1) write(nfout,*) ' !D Reading zaj <<m_ESIO_rd_EVs_ek>>'

    rewind nf

    eko_ek = 0.d0
!!$    do ik = 1, nk_in_the_process-kv3, kv3
    do ik = 1, nk_in_the_process - nspin, nspin
       if(ipri >= 1) write(nfout,*) ' !D     skipping ik = ', ik
!!$       do iks = 1, kv3, af+1
       do iks = 1, nspin, af+1
          do ie = 1, neg
             if(mype == 0) read(nf) wf_l
          end do
       end do
       if(mype == 0) read(nf) n_mpi
       if(mype == 0) read(nf) n2_mpi
       if(mype == 0) read(nf) e_mpi

       if(mype == 0) then
          do iks = 1, nspin
             do ie = 1, neg
                eko_ek(ie,ik+iks-1) = e_mpi(n_mpi(ie,iks),iks)
             end do
          end do
          if(ipri >= 3) then
             do iks=1,nspin
                write(nfout,'(" ik = ",i5)') ik+iks-1
                write(nfout,'(8f8.4)') (e_mpi(n_mpi(ie,iks),iks),ie=1,neg)
             end do
          end if
       end if
    end do

    write(nfout,'(" ---<< m_ESIO_rd_EVs_ek>>---")')
  end subroutine m_ESIO_rd_EVs_ek

  subroutine m_ESIO_wd_Psicoef(ipri,nfout,nf)
    integer, intent(in) :: ipri,nfout, nf
    integer :: ii, jj
    integer, parameter :: Ncol = 5
    integer :: ik, ie, ri,  nel, ig, ib,ib1,ib2,ibt,ibsize
    integer, allocatable, dimension(:) :: n_mpi
    real(DP),allocatable, dimension(:) :: e_mpi
    real(DP),allocatable, dimension(:,:,:) :: wf

    integer :: id_sname = -1
    call tstatc0_begin('m_ESIO_wd_Psicoef ',id_sname)

    allocate(e_mpi(neg)); e_mpi = 0.d0
    allocate(n_mpi(neg)); n_mpi = 0

    if(mype == 0) write(nf,'(" !!COEFFICIENTS of WAVE functions")')

    KPOINT: do ik = 1, kv3, af+1
       nel = min(Nw_Psicoef,iba(ik))
       allocate(wf_l(nel,kimg)); wf_l = 0.d0
!-----------------------------Modified by T.A.Ariasoca-------------------------------------------------
       if ( noncol ) then
             call wd_k_points_noncl()
          else
             call wd_k_points()
          endif
!--------------------------------------------------------------
       e_mpi = 0.d0
       n_mpi = 0
       if(map_k(ik) == myrank_k) then
!!$          if(ipri>=1) then
          if(mype == 0) then
             write(nfout,'(" ik = ", i8, " neordr ")') ik
             write(nfout,'(8i8)') (neordr(ie,ik),ie=1,neg)
          end if
          do ie = 1, neg
!!$             if(map_ek(ie,ik) == mype) n_mpi(ie) = neordr(ie,ik)
!-----------------------------Modified by T.A.Ariasoca-------------------------------------------------
             if ( noncol ) then
		        jj = ik-1
		        ii = ik-mod(jj,2)
                n_mpi(ie) = neordr(ie,ii)
                if(map_e(ie) /= myrank_e) cycle
                e_mpi(ie) = eko_l(map_z(ie),ii)
             else
                n_mpi(ie) = neordr(ie,ik)
                if(map_e(ie) /= myrank_e) cycle
                e_mpi(ie) = eko_l(map_z(ie),ik)
             endif
!--------------------------------------------------------------
          end do
       end if
       if(npes>=2) call mpi_allreduce(MPI_IN_PLACE,e_mpi,neg,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       if(nrank_k>=2) then
          if(map_ek(1,ik) == mype .and. map_ek(1,ik) /= 0) then
             call mpi_send(n_mpi,neg,mpi_integer,0,1,MPI_CommGroup,ierr)
          else if(mype == 0 .and. map_ek(1,ik) /= 0) then
             call mpi_recv(n_mpi,neg,mpi_integer,map_ek(1,ik),1,MPI_CommGroup,istatus,ierr)
          end if
!!$          call mpi_allreduce(MPI_IN_PLACE,n_mpi,neg,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
       end if
       do ie = 1, neg, Ncol
          ib1 = ie; ib2 = min(ie+Ncol-1,neg)
          ibsize = ib2-ib1+1
          if(mype == 0) then
             allocate(wf(nel,kimg,ibsize)); wf=0.d0
          end if
          do ib = ib1, ib2
!!$             ibt = n_mpi(ib)
             ibt = ib
             if(map_ek(ibt,ik) == mype) then                          ! MPI
                do ri = 1, kimg
                   wf_l(1:nel,ri) = zaj_l(1:nel,map_z(ibt),ik,ri)
                end do
                if(map_ek(ibt,ik) /= 0) &                             ! MPI
                     &   call mpi_send(wf_l,nel*kimg,mpi_real4,0,1,MPI_CommGroup,ierr) ! MPI
             else if(mype == 0 .and. map_ek(ibt,ik) /= 0) then
                call mpi_recv(wf_l,nel*kimg,mpi_real4,map_ek(ibt,ik),1,MPI_CommGroup,istatus,ierr)!MPI
             end if
             if(mype == 0) wf(:,:,ib-ib1+1) = wf_l(:,:)
          end do
          call wd_eko(ib1,ib2)
!!$          if(ik <= 2) then
             do ig = 1, nel
                call wd_coef(ibsize,ig)
             end do
!!$          end if
          if(mype == 0) deallocate(wf)
       end do
       deallocate(wf_l)
    end do KPOINT

    call tstatc0_end(id_sname)
  contains
    subroutine wd_coef(ibsize,ig)
      integer, intent(in) :: ibsize,ig
      integer :: ib
!!$      character(3) :: nan
      if(mype == 0) then
         if(kimg==2) then
            write(nf,'(i4," ( ",3i4," )", 5(" (",2f11.5," )"))') &
                 & ig, ngabc(nbase(ig,ik),1:3),(wf(ig,1,ib),wf(ig,2,ib),ib=1,ibsize)
!!$                 & ig, ngabc(nbase(ig,ik),1:3),(wf(ig,1,ie),wf(ig,2,ie),ie=1,3)
         else
!!$            nan="---"
            write(nf,'(i4," ( ",3i4," )", 5(" (",f11.5,4x,"---",4x," )"))') &
                 & ig, ngabc(nbase(ig,ik),1:3),(wf(ig,1,ib),ib=1,ibsize)
!!$            write(nf,'(i4," ( ",3i4," )", 5f11.5)') &
!!$                 & ig, ngabc(nbase(ig,ik),1:3),(wf(ig,1,ie),ie=1,ibsize)
         end if
      end if
    end subroutine wd_coef

    subroutine wd_eko(ib1,ib2)
      integer, intent(in) :: ib1,ib2
      integer :: ie
      if(mype == 0) then
         write(nf,'(a12,5x,a5,5(i7,2x,f10.5,7x))') "ig", "\ e: ", (ie,e_mpi(n_mpi(ie)),ie=ib1,ib2)
      end if
    end subroutine wd_eko

    subroutine wd_k_points()
      if(mype == 0) then
         if(nspin == 1) then
            write(nf,'(" ik = ",i6,"    ( ",3f14.6," )")') ik,(vkxyz(ik,1:3,BUCS))
         else
            if(mod(ik,2) == 1) then
               write(nf,'(" ik = ",i6,"    UP ","    ( ",3f14.6," )")') ik,(vkxyz(ik,1:3,BUCS))
            else
               write(nf,'(" ik = ",i6,"  DOWN ","    ( ",3f14.6," )")') ik,(vkxyz(ik,1:3,BUCS))
            end if
         end if
      end if
    end subroutine wd_k_points
!-----------------------------Modified by T.A.Ariasoca-------------------------------------------------
    subroutine wd_k_points_noncl()
      if(mype == 0) then
         write(nf,'(" ik = ",i6,"    ( ",3f14.6," )")') ik,(vkxyz(ik,1:3,BUCS))
      end if
    end subroutine wd_k_points_noncl
!------------------------------------------------------------------------------
  end subroutine m_ESIO_wd_Psicoef

  subroutine m_ESIO_wd_BandSymInput(ipri,nfout,nf)
    integer, intent(in) :: ipri,nfout, nf
    integer :: ii, jj
    integer :: ik, ie, ri,  nel, ig, ib,ib1,ib2,ibt,ibsize
    integer, allocatable, dimension(:) :: n_mpi
    real(DP),allocatable, dimension(:) :: e_mpi

    integer :: id_sname = -1
    call tstatc0_begin('m_ESIO_wd_BandSymInput ',id_sname)

    allocate(e_mpi(neg)); e_mpi = 0.d0
    allocate(n_mpi(neg)); n_mpi = 0

    if(mype == 0) then
       write(nf,'("##PSIINPSTART")')
       write(nf,'("#MAGNETIC_STATE")')

!-----------------------------Modified by T.A.Ariasoca-------------------------------------------------
       if( noncol ) then
          write(nf,'("2  1  2 ! nspin = 2, af = 1, ndim_spinor = 2 (NONCOLLINEAR)"/"#")')
	      else if(nspin == 2 .and. af == 0) then
          write(nf,'("2  0  1 ! nspin = 2, af = 0, ndim_spinor = 1 (FERRO)"/"#")')
       else if(nspin == 2 .and. af == 1) then
          write(nf,'("2  1  1 ! nspin = 2, af = 1, ndim_spinor = 1 (ANTIFERRO)"/"#")')
       else if(nspin == 1) then
          write(nf,'("1  0  1 ! nspin = 1, af = 0, ndim_spinor = 1 (PARAMAGNETIC)"/"#")')
       else
          write(nf,'(2i5," ! magnetic state = unknown"/"#")') nspin, af
       end if
!-----------------------------------------------------------------------------------
    end if

    KPOINT: do ik = 1, kv3, af+1

       nel = min(Nw_Psicoef,iba(ik))
       allocate(wf_l(nel,kimg)); wf_l = 0.d0
!-----------------------------Modified by T.A.Ariasoca-------------------------------------------------
       if ( noncol ) then
             call wd_k_points_noncl()
          else
             call wd_k_points()
          endif
!----------------------------------------------
       ! --- Eigen Energies --->
       e_mpi = 0.d0
       n_mpi = 0
       if(map_k(ik) == myrank_k) then
          if(ipri>=1) then
             write(nfout,'(" ik = ", i4, " neordr ")') ik
             write(nfout,'(8i8)') (neordr(ie,ik),ie=1,neg)
          end if
          do ie = 1, neg
!!$             if(map_ek(ie,ik) == mype) n_mpi(ie) = neordr(ie,ik)
!-----------------------------Modified by T.A.Ariasoca-------------------------------------------------
             if ( noncol ) then
		        jj = ik-1
		        ii = ik-mod(jj,2)
                n_mpi(ie) = neordr(ie,ii)
                if(map_e(ie) /= myrank_e) cycle
                e_mpi(ie) = eko_l(map_z(ie),ii)
             else
                n_mpi(ie) = neordr(ie,ik)
                if(map_e(ie) /= myrank_e) cycle
                e_mpi(ie) = eko_l(map_z(ie),ik)
             endif
!-----------------------------------------------
          end do
       end if
       if(npes>=2) call mpi_allreduce(MPI_IN_PLACE,e_mpi,neg,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       if(nrank_k>=2) then
          if(map_ek(1,ik) == mype .and. map_ek(1,ik) /= 0) then
             call mpi_send(n_mpi,neg,mpi_integer,0,1,MPI_CommGroup,ierr)
          else if(mype == 0 .and. map_ek(1,ik) /= 0) then
             call mpi_recv(n_mpi,neg,mpi_integer,map_ek(1,ik),1,MPI_CommGroup,istatus,ierr)
          end if
!!$          call mpi_allreduce(MPI_IN_PLACE,n_mpi,neg,mpi_integer,mpi_sum,MPI_CommGroup,ierr)
       end if
       call wd_eko(1, neg)
       ! <---
       if(mype == 0) write(nf,'("#Number_of_Gvectors"/,i0,/,"#")') nel
       call wd_gvectors(ik)

       ! --- Wave Function Coefficients --->
      if(mype == 0)  write(nf,'("#PSICOEF")')
       do ie = 1, neg
          if(mype == 0)  write(nf,'(i0,"  #EIGEN_ENERGY NUMBER")') ie
          if(map_ek(ie,ik) == mype) then
             wf_l(1:nel,1:kimg) = zaj_l(1:nel,map_z(ie),ik,1:kimg)
             if(map_ek(ie,ik) /= 0) &
                  &   call mpi_send(wf_l,nel*kimg,mpi_real,0,1,MPI_CommGroup,ierr) ! MPI
          else if(mype == 0 .and. map_ek(ie,ik) /= 0) then
             call mpi_recv(wf_l,nel*kimg,mpi_real,map_ek(ie,ik),1,MPI_CommGroup,istatus,ierr)!MPI
          end if
          call wd_coef2()
       end do
       if(mype == 0) write(nf,'("#")')
       deallocate(wf_l)
    end do KPOINT
    if(mype == 0) write(nf,'("##")')

    call tstatc0_end(id_sname)
  contains
    subroutine wd_coef2()
      integer :: ig
      if(mype == 0) then
         if(kimg==2) then
            write(nf,'(8f11.5)') (wf_l(ig,1),wf_l(ig,2),ig=1,nel)
         else
            write(nf,'(8f11.5)') (wf_l(ig,1), 0.d0, ig=1,nel)
         end if
      end if
    end subroutine wd_coef2

    subroutine wd_gvectors(ik)
      integer, intent(in) :: ik
      integer :: i
      if(mype == 0) then
         write(nf,'("#GVECTOR")')
         do i = 1, nel
            write(nf,'(i0,1x,3i6)') i,ngabc(nbase(i,ik),1:3)
         end do
         write(nf,'("#")')
      end if
    end subroutine wd_gvectors

    subroutine wd_nelements()
      if(mype == 0) write(nf,'("#Number_of_Gvectors"/,i0,/,"#")') nel
    end subroutine wd_nelements

    subroutine wd_eko(ib1,ib2)
      integer, intent(in) :: ib1,ib2
      integer :: ie
      if(mype == 0) then
         write(nf,'("#EIGEN_ENERGY")')
         do ie = ib1,ib2
            write(nf,'(i0,2x,f10.5)') ie,e_mpi(n_mpi(ie))
         end do
         write(nf,'("#")')
      end if
    end subroutine wd_eko

    subroutine wd_k_points()
      character(len=4),dimension(3) :: spinstate
      data spinstate/"    ","  UP","DOWN"/
      integer :: ip
      if(mype == 0) then
         if(nspin == 1) then
            ip = 1
         else
            if(mod(ik,2) == 1) then
               ip = 2
            else
               ip = 3
            end if
         end if
         write(nf,'("#KPOINT"/,i06,3f14.6,2x,a4,/"#")') ik,(vkxyz(ik,1:3,BUCS)),spinstate(ip)
      end if
    end subroutine wd_k_points
!-----------------------------Modified by T.A.Ariasoca-------------------------------------------------
    subroutine wd_k_points_noncl()
      write(nf,'("#KPOINT"/,i06,3f14.6,2x,a4,/"#")') ik,(vkxyz(ik,1:3,BUCS))
    end subroutine wd_k_points_noncl
!------------------------------------------------------------------------------

  end subroutine m_ESIO_wd_BandSymInput

  subroutine m_ESIO_wd_WFs_and_EVs_ek(nfout,nf)
    integer, intent(in) :: nfout, nf
    integer  :: ik, ie, ri, ikg, ikt
    integer, allocatable, dimension(:,:) :: n_mpi
    real(DP),allocatable, dimension(:,:) :: e_mpi
    integer :: id_sname = -1
    call tstatc0_begin('m_ESIO_wd_WFs_and_EVs_ek ',id_sname)

    if(precision_WFfile==SP) then
       allocate(wf_l(kg1,kimg));    wf_l = 0
    else
       allocate(wfdp_l(kg1,kimg));    wfdp_l = 0
    end if

    call mpi_barrier(MPI_CommGroup,ierr)
    if(ipri>=2) write(nfout,'(" !D Writing WaveFunctions ")')

    allocate(n_mpi(neg,nspin)); n_mpi = 0
    allocate(e_mpi(neg,nspin)); e_mpi = 0.d0

    KPOINT: do ikg = 1, nrank_k
       ! ---> zaj_l
       if((ikg-1)*nspin + 1 > kv3) then
          if(ipri >= 1) write(nfout,'(" !D ik > kv3")')
          exit KPOINT
       end if

       do ikt = 1, nspin, af+1
          ik = (ikg-1)*nspin + ikt
          if(ipri>=1) write(nfout,'(" !D Writing WaveFunctions ik = ",i5)') ik
          do ie = 1, neg
           if(precision_WFfile==SP) then
             if(map_ek(ie,ik) == mype) then                          ! MPI
                do ri = 1, kimg
                   wf_l(1:kg1,ri) = zaj_l(1:kg1,map_z(ie),ik,ri)
                end do
                if(map_ek(ie,ik) /= 0) &                             ! MPI
                     &   call mpi_send(wf_l,kg1*kimg,mpi_real4,0,1,MPI_CommGroup,ierr) ! MPI
             else if(mype == 0 .and. map_ek(ie,ik) /= 0) then
                call mpi_recv(wf_l,kg1*kimg,mpi_real4,map_ek(ie,ik),1,MPI_CommGroup,istatus,ierr)!MPI
             end if
             if(mype == 0) write(nf)  wf_l                        ! MPI
           else if(precision_WFfile==DP) then
             if(map_ek(ie,ik) == mype) then
                do ri = 1, kimg
                   wfdp_l(1:kg1,ri) = zaj_l(1:kg1,map_z(ie),ik,ri)
                end do
                if(map_ek(ie,ik) /= 0) &
                     &   call mpi_send(wfdp_l,kg1*kimg,mpi_double_precision,0,1,MPI_CommGroup,ierr)
             else if(mype == 0 .and. map_ek(ie,ik) /= 0) then
                call mpi_recv(wfdp_l,kg1*kimg,mpi_double_precision,map_ek(ie,ik),1,MPI_CommGroup,istatus,ierr)
             end if
             if(mype == 0) write(nf)  wfdp_l
           end if
          end do
       end do

       ! --->  neordr
       n_mpi = 0                                          ! MPI
       do ikt = 1, nspin                                  ! MPI
          ik = (ikg-1)*nspin + ikt
          if(map_k(ik) /= myrank_k) cycle                 ! MPI
          n_mpi(1:neg,ikt) = neordr(1:neg,ik)             ! MPI
       end do                                             ! MPI
       if(npes >= 2) then
          call mpi_allreduce(MPI_IN_PLACE,n_mpi,neg*nspin,mpi_integer,mpi_sum &
               &                      ,MPI_CommGroup,ierr)  ! MPI
          n_mpi = n_mpi/nrank_e
       end if
       if(ipri>=2) write(nfout,'(" !D Writing neordr ik = ",i5)') ik
       if(mype == 0) write(nf) n_mpi             ! MPI ; writing (neordr)

       ! --->  nrvf_ordr
       n_mpi = 0                                          ! MPI
       do ikt = 1, nspin                                  ! MPI
          ik = (ikg-1)*nspin + ikt
          if(map_k(ik) /= myrank_k) cycle                 ! MPI
          n_mpi(1:neg,ikt) = nrvf_ordr(1:neg,ik)          ! MPI
       end do                                             ! MPI
       if(npes >= 2) then
          call mpi_allreduce(MPI_IN_PLACE,n_mpi,neg*nspin,mpi_integer,mpi_sum &
               &                      ,MPI_CommGroup,ierr)  ! MPI
          n_mpi = n_mpi/nrank_e
       end if
       if(ipri>=2) write(nfout,'(" !D Writing nrvf_ordr ik = ",i5)') ik
       if(mype == 0) write(nf) n_mpi             ! MPI ; writing (nrvf_ordr)

       e_mpi = 0.d0                                       ! MPI
       do ikt = 1, nspin                                  ! MPI
          ik = (ikg-1)*nspin + ikt
          if(map_k(ik) /= myrank_k) cycle                 ! MPI
          do ie = 1, neg                                  ! MPI
             if(map_e(ie) /= myrank_e) cycle              ! MPI
             e_mpi(ie,ikt) = eko_l(map_z(ie),ik)           ! MPI
          end do
       end do
       if(npes >= 2) then
          call mpi_allreduce(MPI_IN_PLACE,e_mpi,neg*nspin,mpi_double_precision &
               &               ,mpi_sum,MPI_CommGroup,ierr) ! MPI
       end if
       if(ipri>=2) write(nfout,'(" !D Writing eko_l ik = ",i5)') ik
       if(mype == 0) write(nf) e_mpi             ! MPI ; writing (eko_l)
    end do KPOINT

!!$    do ik = 1, kv3, af+1
!!$       do ie = 1, neg
!!$          if(map_ek(ie,ik) == mype) then                          ! MPI
!!$             do ri = 1, kimg
!!$                wf_l(1:kg1,ri) = zaj_l(1:kg1,map_z(ie),ik,ri)
!!$             end do
!!$             if(map_ek(ie,ik) /= 0) &                             ! MPI
!!$            &   call mpi_send(wf_l,kg1*kimg,mpi_real,0,1,MPI_CommGroup,ierr) ! MPI
!!$          else if(mype == 0 .and. map_ek(ie,ik) /= 0) then
!!$             call mpi_recv(wf_l,kg1*kimg,mpi_real,map_ek(ie,ik),1,MPI_CommGroup,istatus,ierr)!MPI
!!$          end if
!!$          if(mype == 0) write(nf)  wf_l                        ! MPI
!!$       end do
!!$    end do
!!$    deallocate(wf_l)
!!$
!!$    allocate(n_mpi(neg,kv3)); allocate(n2_mpi(neg,kv3))! MPI
!!$    allocate(e_mpi(neg,kv3)); allocate(e2_mpi(neg,kv3))! MPI
!!$
!!$    n_mpi = 0                                          ! MPI
!!$    do ik = 1, kv3                                     ! MPI
!!$       if(map_k(ik) /= myrank_k) cycle                 ! MPI
!!$       n_mpi(1:neg,ik) = neordr(1:neg,ik)              ! MPI
!!$    end do                                             ! MPI
!!$    if(npes >= 2) then
!!$       call mpi_allreduce(n_mpi,n2_mpi,neg*kv3,mpi_integer,mpi_sum &
!!$            &                      ,MPI_CommGroup,ierr)  ! MPI
!!$       n2_mpi = n2_mpi/nrank_e
!!$    else
!!$       n2_mpi = n_mpi
!!$    end if
!!$    if(mype == 0) write(nf) n2_mpi             ! MPI ; writing (neordr)
!!$
!!$    n_mpi = 0                                          ! MPI
!!$    do ik = 1, kv3                                     ! MPI
!!$       if(map_k(ik) /= myrank_k) cycle                 ! MPI
!!$       n_mpi(1:neg,ik) = nrvf_ordr(1:neg,ik)           ! MPI
!!$    end do                                             ! MPI
!!$    if(npes >= 2) then
!!$       call mpi_allreduce(n_mpi,n2_mpi,neg*kv3,mpi_integer,mpi_sum &
!!$            &                      ,MPI_CommGroup,ierr)  ! MPI
!!$       n2_mpi = n2_mpi/nrank_e
!!$    else
!!$       n2_mpi = n_mpi
!!$    end if
!!$    if(mype == 0) write(nf) n2_mpi             ! MPI ; writing (nrvf_ordr)
!!$
!!$    e_mpi = 0.d0                                       ! MPI
!!$    do ik = 1, kv3                                     ! MPI
!!$       if(map_k(ik) /= myrank_k) cycle                 ! MPI
!!$       do ie = 1, neg                                  ! MPI
!!$          if(map_e(ie) /= myrank_e) cycle              ! MPI
!!$          e_mpi(ie,ik) = eko_l(map_z(ie),ik)           ! MPI
!!$       end do
!!$    end do
!!$    if(npes >= 2) then
!!$       call mpi_allreduce(e_mpi,e2_mpi,neg*kv3,mpi_double_precision &
!!$            &               ,mpi_sum,MPI_CommGroup,ierr) ! MPI
!!$    else
!!$       e2_mpi = e_mpi
!!$    end if
!!$    if(mype == 0) write(nf) e2_mpi             ! MPI ; writing (eko_l)

!!$    e_mpi = 0.d0                                       ! MPI
!!$    do ik = 1, kv3                                     ! MPI
!!$       if(map_k(ik) /= myrank_k) cycle                 ! MPI
!!$       do ie = 1, neg                                  ! MPI
!!$          if(map_e(ie) /= myrank_e) cycle              ! MPI
!!$          e_mpi(ie,ik) = occup_l(map_z(ie),ik)         ! MPI
!!$       end do                                          ! MPI
!!$    end do                                             ! MPI
!!$    if(npes >= 2) then
!!$       call mpi_allreduce(e_mpi,e2_mpi,neg*kv3,mpi_double_precision &
!!$            &                  ,mpi_sum,MPI_CommGroup,ierr) ! MPI
!!$    else
!!$       e2_mpi = e_mpi
!!$    end if
!!$    if(mype == 0) write(nf) e2_mpi             ! MPI ; writing (occup_l)

    if(precision_WFfile==SP) then
       deallocate(wf_l)
    else if(precision_WFfile==DP) then
       deallocate(wfdp_l)
    end if
    deallocate(n_mpi)
    deallocate(e_mpi)

    call tstatc0_end(id_sname)

  end subroutine m_ESIO_wd_WFs_and_EVs_ek

  subroutine m_ESIO_wd_WFn(nfout,nfwfk,ik,ib)
    integer, intent(in) :: nfout,nfwfk
    integer, intent(in) :: ik,ib
    integer :: jb
    real(kind=DP), allocatable :: bfft(:)
    real(kind=DP) :: eig
    integer :: id_sname = -1
    call tstatc0_begin('m_ESIO_wd_WFn ',id_sname)

    allocate(bfft(nfft));    bfft = 0
    call mpi_barrier(MPI_CommGroup,ierr)
    if(ipri >= 1) write(nfout,*) ' !D Writing Wavefunctions '
    rewind nfwfk

    if(map_k(ik) == myrank_k) jb = neordr(ib,ik)
    if(nrank_k > 1) call mpi_bcast(jb,1,mpi_integer,map_k(ik),mpi_e_world(myrank_e),ierr)

    if(map_ek(jb,ik) == mype) then
       call m_ES_WF_in_Rspace(ik,jb,bfft)
       eig = eko_l(map_z(jb),ik)
       if(map_ek(jb,ik) /= 0) then
          call mpi_send(bfft,nfft,mpi_double_precision,0,1,MPI_CommGroup,ierr)
          call mpi_send(eig,1,mpi_double_precision,0,1,MPI_CommGroup,ierr)
       end if
    else if(mype == 0 .and. map_ek(jb,ik) /= 0) then
       call mpi_recv(bfft,nfft,mpi_double_precision,map_ek(jb,ik),1,MPI_CommGroup,istatus,ierr)
       call mpi_recv(eig,1,mpi_double_precision,map_ek(jb,ik),1,MPI_CommGroup,istatus,ierr)
    end if
    if(mype == 0) then
       call wd_wf(nfout,nfwfk,ik,ib,eig)
    end if
    deallocate(bfft)
    call tstatc0_end(id_sname)
  contains
    subroutine wd_wf(nfout,nfwfk,ik,ib,eig)
      integer, intent(in) :: nfout,nfwfk,ik,ib
      real(kind=DP), intent(in) :: eig
      integer :: i,j,k, id, nl, nm, nn, nlhf,inew,jnew,knew,ip,mm
      real(kind=DP),allocatable,dimension(:,:,:,:) :: wkwf
      real(kind=DP) ::      x,y,z
      integer, parameter :: UP = 1 , DOWN = 2
      integer ::            up_down
      real(kind=DP),allocatable,dimension(:,:) :: cps_full
      integer, allocatable,dimension(:) :: ityp_full
      integer :: m, nk
      real(kind=DP) :: norm
      real(kind=DP) :: normr,normi
      integer :: n1,n2,n3
      real(kind=DP) :: dn1,dn2,dn3
      integer :: icomp

      id = fft_box_size_WF(1,0)
      mm = fft_box_size_WF(2,0)
      nl = fft_box_size_WF(1,1)
      nm = fft_box_size_WF(2,1)
      nn = fft_box_size_WF(3,1)

      if(kimg == 1) then
         nlhf = id/2
      else
         nlhf = id
      end if

      if(wf_filetype == DENSITY_ONLY &
      & .or. wf_filetype == VTK &
      & .or. wf_filetype == BINARY) then
         allocate(wkwf(nl,nm,nn,2)); wkwf = 0.d0
      else if(wf_filetype == CUBE) then
         allocate(wkwf(nn,nm,nl,2)); wkwf = 0.d0
      end if

      if(nspin==2) then
         if(mod(ik/nspin,2) == 1) then
            up_down = UP
         else
            up_down = DOWN
         end if
      end if
      nk = (ik-1)/nspin+1

      if(ipri >= 2) write(nfout,9001) nl*nm*nn, nl, nm, nn
9001  format(' Wavefunction ',i8,'(',3i5,')')

      if(ipri >= 2) write(nfout,*) ' !D FFT cube mapping start'
      do i = 1, nm
         do j = 1, nn
            do k = 1, nl
               if(kimg == 1 .and. k > nlhf) then
                  knew = id - k
                  jnew = nn+2 - j
                  inew = nm+2 - i
                  if(jnew > nn) then
                     jnew = jnew - nn
                  end if
                  if(inew > nm) then
                     inew = inew - nm
                  end if
               else
                  knew = k; jnew = j; inew = i
               end if
               ip = nlhf*mm*(jnew-1) + nlhf*(inew-1) + knew
               if(wf_filetype == DENSITY_ONLY &
               & .or. wf_filetype == VTK &
               & .or. wf_filetype == BINARY) then
                  wkwf(k,i,j,1) = bfft(ip*2-1)
                  wkwf(k,i,j,2) = bfft(ip*2)
               else if(wf_filetype == CUBE) then
                  wkwf(j,i,k,1) = bfft(ip*2-1)
                  wkwf(j,i,k,2) = bfft(ip*2)
               end if
            end do
         end do
      end do
! Normalization
      normr = 0.d0
      normi = 0.d0
      do i = 1, nm
         do j = 1, nn
            do k = 1, nl
               if(wf_filetype == DENSITY_ONLY &
               & .or. wf_filetype == VTK &
               & .or. wf_filetype == BINARY) then
                  normr = normr + wkwf(k,i,j,1)*wkwf(k,i,j,1)
                  normi = normi + wkwf(k,i,j,2)*wkwf(k,i,j,2)
               else if(wf_filetype == CUBE) then
                  normr = normr + wkwf(j,i,k,1)*wkwf(j,i,k,1)
                  normi = normi + wkwf(j,i,k,2)*wkwf(j,i,k,2)
               end if
            end do
         end do
      end do
      norm = normr + normi
      write(nfout,*) 'Real and imaginary parts of wf = ',normr/norm,normi/norm
      if(wf_filetype == VTK .or. wf_filetype == BINARY) then
         if(normr>normi) then
            norm = normr
            icomp = 1
            write(nfout,*) 'Real part of wf will be outputed.'
         else
            norm = normi
            icomp = 2
            write(nfout,*) 'Imaginary part of wf will be outputed.'
         end if
      end if
      norm = univol*norm/dble(nm*nn*nl)
      norm = 1.d0/dsqrt(norm)
      do i = 1, nm
         do j = 1, nn
            do k = 1, nl
               if(wf_filetype == DENSITY_ONLY &
               & .or. wf_filetype == VTK &
               & .or. wf_filetype == BINARY) then
                  wkwf(k,i,j,1) = norm*wkwf(k,i,j,1)
                  wkwf(k,i,j,2) = norm*wkwf(k,i,j,2)
               else if(wf_filetype == CUBE) then
                  wkwf(j,i,k,1) = norm*wkwf(j,i,k,1)
                  wkwf(j,i,k,2) = norm*wkwf(j,i,k,2)
               end if
            end do
         end do
      end do

      if(wf_filetype == DENSITY_ONLY) then
         write(nfwfk,9001) nl*nm*nn, nl, nm, nn
         write(nfwfk,'(6e13.5)') wkwf
      else if(wf_filetype == BINARY) then
         write(nfwfk) nl*nm*nn, nl, nm, nn
         write(nfwfk) altv, nspin, up_down, nk, ib, eig
         write(nfwfk) wkwf(:,:,:,icomp)
      else if(wf_filetype == VTK) then
         write(nfwfk,'("# vtk DataFile Version 2.0")')
         if(nspin == 2) then
            if(up_down == 1) then
               write(nfwfk,'(" SCF Wavefunction UP : k=",i7," n=",i7," eig=",f20.5)') nk,ib,eig
            else
               write(nfwfk,'(" SCF Wavefunction DOWN : k=",i7," n=",i7," eig=",f20.5)') nk,ib,eig
            end if
         else
            write(nfwfk,'(" SCF Wavefunction : k=",i7," n=",i7," eig=",f20.5)') nk, ib, eig
         end if
         write(nfwfk,'("ASCII")')
         write(nfwfk,'("DATASET STRUCTURED_GRID")')
         write(nfwfk,'("DIMENSIONS",3(1x,i5))') nl+1,nm+1,nn+1
         write(nfwfk,'("POINTS",1x,i10,1x,"float")') (nl+1)*(nm+1)*(nn+1)
         do n1=0,nl
            do n2=0,nm
               do n3=0,nn
                  dn1 = n1/dble(nl)
                  dn2 = n2/dble(nm)
                  dn3 = n3/dble(nn)
                  x = altv(1,1)*dn1 + altv(1,2)*dn2 + altv(1,3)*dn3
                  y = altv(2,1)*dn1 + altv(2,2)*dn2 + altv(2,3)*dn3
                  z = altv(3,1)*dn1 + altv(3,2)*dn2 + altv(3,3)*dn3
                  write(nfwfk,'(3(1x,e13.5))') x,y,z
               end do
            end do
         end do
         write(nfwfk,'("")')
         write(nfwfk,'("POINT_DATA",1x,i10)') (nl+1)*(nm+1)*(nn+1)
         write(nfwfk,'("SCALARS scalars float")')
         write(nfwfk,'("LOOKUP_TABLE default")')
         do n1=0,nl
            i=n1+1
             if(n1==nl) i=1
            do n2=0,nm
               j=n2+1
               if(n2==nm) j=1
               do n3=0,nn
                  k=n3+1
                  if(n3==nn) k=1
                  write(nfwfk,'(e13.5)') wkwf(i,j,k,icomp)
               end do
            end do
         end do
      else if(wf_filetype == CUBE) then
         if(len_trim(wf_title) >= 1) then
            write(nfwfk,*) trim(wf_title)
         else
            write(nfwfk,'(" Calculated by phase")')
         end if
         if(nspin == 2) then
            if(up_down == 1) then
               write(nfwfk,'(" SCF Wavefunction UP : k=",i7," n=",i7," eig=",f20.5)') nk,ib,eig
            else
               write(nfwfk,'(" SCF Wavefunction DOWN : k=",i7," n=",i7," eig=",f20.5)') nk,ib,eig
            end if
         else
            write(nfwfk,'(" SCF Wavefunction : k=",i7," n=",i7," eig=",f20.5)') nk, ib, eig
         end if
         x = 0.d0; y = 0.d0; z = 0.d0
         write(nfwfk,'(i6,3f10.4)') -natm2, x,y,z
         do i = 1, 3
            write(nfwfk,'(i6,3f10.6)') fft_box_size_WF(i,1), altv(1:3,i)/dble(fft_box_size_WF(i,1))
         end do

         allocate(cps_full(natm2,3))
         allocate(ityp_full(natm2))
         cps_full = 0; ityp_full = 0
         call m_IS_pack_all_ions_in_uc(ityp_full,cps_full)
         do i = 1, natm2
            m = ityp_full(i)
            write(nfwfk,'(f8.4,4f10.6)') iatomn(m), ival(m), cps_full(i,1:3)
         end do
         deallocate(ityp_full,cps_full)

         write(nfwfk,'(10i5)') 2,1,2
         write(nfwfk,'(6e13.5)') wkwf(:,:,:,1)
         write(nfwfk,'(6e13.5)') wkwf(:,:,:,2)

      end if
      if(allocated(wkwf)) deallocate(wkwf)
    end subroutine wd_wf
  end subroutine m_ESIO_wd_WFn

  logical function m_ESIO_check_energy(ik,ib)
    integer, intent(in) :: ik,ib
    integer :: jb
    real(kind=DP) :: eig
    if(map_k(ik) == myrank_k) jb = neordr(ib,ik)
    if(nrank_k > 1) call mpi_bcast(jb,1,mpi_integer,map_k(ik),mpi_e_world(myrank_e),ierr)

    if(map_ek(jb,ik) == mype) then
       eig = eko_l(map_z(jb),ik)
       if(eig >= eigmin_wf .and. eig <= eigmax_wf) then
          m_ESIO_check_energy = .true.
       else
          m_ESIO_check_energy = .false.
       end if
    end if
    call mpi_bcast(m_ESIO_check_energy,1,mpi_logical,map_ek(jb,ik),MPI_CommGroup,ierr)
  end function m_ESIO_check_energy

  subroutine m_ESIO_wd_Efermi(nfout,nfefermi)
    integer, intent(in) :: nfout, nfefermi
                                                  __TIMER_SUB_START(1377)
    if(mype == 0) then
                                                  __TIMER_IODO_START(1460)
       write(nfefermi,'(f16.8," : Efermi")') efermi
       if(sw_fix_total_spin == YES .and. nspin == 2) then
          write(nfefermi,'(2f16.8," : Ffermi_spin(1), Efermi_spin(2)")') &
               & efermi_spin(1),efermi_spin(2)
       end if
                                                  __TIMER_IODO_STOP(1460)
    end if
                                                  __TIMER_SUB_STOP(1377)
  end subroutine m_ESIO_wd_Efermi

  subroutine m_ESIO_rd_Efermi(nfout,nfefermi)
    integer, intent(in) :: nfout, nfefermi
                                                  __TIMER_SUB_START(1376)
    if(mype == 0) then
                                                  __TIMER_IODO_START(1458)
       read(nfefermi,*,err=1001,end=1001) efermi
       write(nfout,'(" ! efermi = ",f16.8," : this is read from nfefermi")') efermi
       if(sw_fix_total_spin == YES .and. nspin == 2) then
          read(nfefermi,*,err=1002,end=1002) efermi_spin(1),efermi_spin(2)
          write(nfout,'(" ! efermi_spin = ",2f16.8," : these are read from nfefermi")') efermi_spin(1:2)
       end if
                                                  __TIMER_IODO_STOP(1458)
       goto 1010
1001   continue
       efermi = 0.d0
1002   continue
       if(sw_fix_total_spin == YES .and. nspin == 2) then
          efermi_spin(1) = 0.d0; efermi_spin(2) = 0.d0
       end if
1010   continue
    end if
    if(npes > 1) then
                                                  __TIMER_IOCOMM_START_w_BARRIER(MPI_CommGroup,1459)
       call mpi_bcast(efermi,1,mpi_double_precision,0,MPI_CommGroup,ierr)
       if(sw_fix_total_spin == YES .and. nspin == 2) then
          call mpi_bcast(efermi_spin,2,mpi_double_precision,0,MPI_CommGroup,ierr)
       end if
                                                  __TIMER_IOCOMM_STOP(1459)
    end if
                                                  __TIMER_SUB_STOP(1376)
  end subroutine m_ESIO_rd_Efermi


! === KT_add === 2015/05/16
  subroutine m_ESIO_wd_Wfn_squared_noncl( nfout )
    integer, intent(in) :: nfout

    integer :: iloop, ib
    logical :: goflag
    real(kind=DP), allocatable :: chgq_l_save(:,:,:)
    real(kind=DP), allocatable :: hsr_save(:,:,:,:), hsi_save(:,:,:,:)

    if ( ekmode == ON ) then
       if ( nk_in_the_process <= ik_wf_squared &
            &  .and. nk_in_the_process -1 +kv3 > ik_wf_squared ) then
          goflag = .true.
       else
          goflag = .false.
       endif
    else
       goflag = .true.
    endif

    if ( .not. goflag ) return

    allocate( chgq_l_save(ista_kngp:iend_kngp,kimg,ndim_magmom) )
    allocate( hsr_save(natm,nlmt,nlmt,ndim_magmom) )
    allocate( hsi_save(natm,nlmt,nlmt,ndim_magmom) )

    chgq_l_save = chgq_l;     hsr_save = hsr;     hsi_save = hsi

    Do ib=ib1_wf_squared, ib2_wf_squared
       if ( ekmode == ON ) then
          call m_CD_softpart_ktsub_noncl( nfout, kv3, &
               &                          ik_wf_squared -nk_in_the_process+1, ib )
          call m_CD_hardpart_ktsub_noncl( nfout, ik_wf_squared -nk_in_the_process +1, &
               &                          ib )
       else
          call m_CD_softpart_ktsub_noncl( nfout, kv3, ik_wf_squared, ib )
          call m_CD_hardpart_ktsub_noncl( nfout, ik_wf_squared, ib )
       endif

       call m_CD_alloc_rspace_charge()
       Do iloop=1, ndim_magmom
          call m_Files_open_nfwfksq_noncl( iloop, ik_wf_squared, ib )
          call m_CD_rspace_charge_noncl( iloop, nfwfk_sq, nfout, wf_squared_filetype )
       End do
       call m_CD_dealloc_rspace_charge()

    End Do

    chgq_l = chgq_l_save;     hsr_save = hsr;     hsi_save = hsi
    deallocate( chgq_l_save ); deallocate( hsr_save );  deallocate( hsi_save )

  end subroutine m_ESIO_wd_Wfn_squared_noncl
! ============== 2015/05/16

  subroutine m_ESIO_wd_Wfn_integ_magmom       ! noncl
    integer :: ik, ie, ib1, ig, is1, is2, ik2
    integer :: ia, it, lmt1, lmt2, p1, p2
    integer :: neg_t
    real(kind=DP) :: csum(3)
    complex(kind=CMPLDP) :: z1, z2
    complex(kind=CMPLDP) :: PauliMatrix( ndim_magmom, ndim_spinor, ndim_spinor )

    real(kind=DP), allocatable :: work(:,:,:), magmom_each_wfn(:,:,:)

    call m_ES_set_Pauli_matrix( PauliMatrix )

    allocate( magmom_each_wfn( kv3/ndim_spinor, neg, 3 ) )
    magmom_each_wfn = 0.0d0

    Do ik=1, kv3, ndim_spinor
       if ( map_k(ik) /= myrank_k ) cycle! MPI

       Do ie=ista_e, iend_e
          ib1 = map_z(ie)

          csum = 0.0d0

          Do ig=1, iba(ik)
             Do is1=1, ndim_spinor
                Do is2=1, ndim_spinor
                   z1 = cmplx( zaj_l( ig, ib1, ik+is1-1, 1 ), &
                        &      zaj_l( ig, ib1, ik+is1-1, kimg ) )
                   z2 = cmplx( zaj_l( ig, ib1, ik+is2-1, 1 ), &
                        &      zaj_l( ig, ib1, ik+is2-1, kimg ) )

                   csum(1) = csum(1) +conjg(z1) *PauliMatrix(2,is1,is2) *z2
                   csum(2) = csum(2) +conjg(z1) *PauliMatrix(3,is1,is2) *z2
                   csum(3) = csum(3) +conjg(z1) *PauliMatrix(4,is1,is2) *z2
                End do
             End Do
          End Do
!
          Do ia=1, natm
             it = ityp(ia)
             do lmt1 = 1, ilmt(it)
                p1 = lmta(lmt1,ia)
                do lmt2 = 1, ilmt(it)
                   p2 = lmta(lmt2,ia)
                   Do is1=1, ndim_spinor
                      Do is2=1, ndim_spinor
                         z1 = cmplx( fsr_l(ib1,p1,ik+is1-1), fsi_l(ib1,p1,ik+is1-1) )
                         z2 = cmplx( fsr_l(ib1,p2,ik+is2-1), fsi_l(ib1,p2,ik+is2-1) )

                         csum(1) = csum(1) +conjg(z1) *q(lmt1,lmt2,it) &
                              &                       *PauliMatrix(2,is1,is2) *z2
                         csum(2) = csum(2) +conjg(z1) *q(lmt1,lmt2,it) &
                              &                       *PauliMatrix(3,is1,is2) *z2
                         csum(3) = csum(3) +conjg(z1) *q(lmt1,lmt2,it) &
                              &                       *PauliMatrix(4,is1,is2) *z2
                      End do
                   End do
                End do
             End do
          End Do
!
          ik2 = (ik-1)/ndim_spinor +1
          magmom_each_wfn( ik2, ie, 1:3 ) = csum(1:3)
       End Do
    End Do
!
    if ( npes > 1 ) then
       allocate( work (kv3/ndim_spinor, neg, 3 ) ); work = 0.0d0
       call mpi_allreduce( magmom_each_wfn, work, kv3/ndim_spinor *neg *3, &
            &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
       magmom_each_wfn = work
       deallocate( work )
    endif

!
    neg_t = neg -num_extra_bands

    if ( ekmode == on ) then
       if ( mype == 0 ) then
          if ( nk_in_the_process == 1 ) then
             call m_Files_open_nfwfk_integ_mom(2)
             write(nfwfk_integ_mom,'(A,I8)') 'num_kpoints = ', kv3_ek /ndim_spinor
             write(nfwfk_integ_mom,'(A,I8)') 'num_bands   = ', neg_t
             write(nfwfk_integ_mom,'(A,I8)') 'nspin       = ', nspin /ndim_spinor
                                                ! in order to work band**pl properly
             write(nfwfk_integ_mom,*)
          else
             call m_Files_open_nfwfk_integ_mom(3)
          endif
       endif
    else
       if ( mype == 0 ) then
          call m_Files_open_nfwfk_integ_mom(icond)
          write(nfwfk_integ_mom,'(A,I8)') 'num_kpoints = ', kv3 /ndim_spinor
          write(nfwfk_integ_mom,'(A,I8)') 'num_bands   = ', neg_t
          write(nfwfk_integ_mom,'(A,I8)') 'nspin       = ', nspin
          write(nfwfk_integ_mom,*)
       endif
    endif

    if ( ekmode == ON ) then
       if ( mype == 0 ) then
          Do ik=1, kv3, ndim_spinor
             write(nfwfk_integ_mom,'(A)') "================= "
             write(nfwfk_integ_mom,'(" ik = ",i4," (",3f10.6," )")') &
                  &                   ik +nk_in_the_process-1, (vkxyz(ik,1:3,BUCS))
             Do ie=1, neg_t
                ik2 = ( ik -1 )/ndim_spinor +1
                write(nfwfk_integ_mom,'(I5,3F16.8)') ie, magmom_each_wfn( ik2, ie, 1:3 )
             End Do
          End Do
          write(nfwfk_integ_mom,*)
       endif
    else
       if ( mype == 0 ) then
          Do ik=1, kv3, ndim_spinor
             write(nfwfk_integ_mom,'(A)') "================= "
             write(nfwfk_integ_mom,'(" ik = ",i4," (",3f10.6," )")') &
                  &                   ik, (vkxyz(ik,1:3,BUCS))
             Do ie=1, neg_t
                ik2 = ( ik -1 )/ndim_spinor +1
                write(nfwfk_integ_mom,'(I5,3F16.8)') ie, magmom_each_wfn( ik2, ie, 1:3 )
             End Do
          End Do
          write(nfwfk_integ_mom,*)
       endif
    endif

    deallocate( magmom_each_wfn )

    call m_Files_close_nfwfk_integ_mom

  end subroutine m_ESIO_wd_Wfn_integ_magmom

  subroutine m_ESIO_wd_EigenValues_bxsf
    integer :: ik, ib,jb,ibo,jbo, neg_t
    integer :: lun, nbz_mesh(3), ispin
    integer :: kv3_wk

    integer, allocatable :: neordr_t(:,:), kvtab(:)
    real(kind=DP), parameter :: delta = 1.d-12

    if ( ekmode == ON ) then
       allocate(neordr_t(neg,kv3_ek));   neordr_t = 0
       neg_t = neg
       if(neg_is_enlarged) neg_t = neg -num_extra_bands

       do ik = 1, kv3_ek, ndim_spinor
          if (nspin == 1 .or. (nspin == 2 .and. mod(ik,2) == 1)) &
               & neordr_t(1:neg,ik) = (/(ib,ib=1,neg)/)
          do ib = 1, neg-1
             do jb = ib+1, neg
                ibo = neordr_t(ib,ik);  jbo = neordr_t(jb,ik)
                if ( eko_ek(jbo,ik) < eko_ek(ibo,ik)-delta ) then        ! MPI
                   neordr_t(jb,ik) = ibo;    neordr_t(ib,ik) = jbo
                end if
             end do
          end do
       end do
    else
       !write(*,*) "ekmode == Off is not supported "
       call phase_error_with_msg(6,"ekmode == Off is not supported ",__LINE__,__FILE__)
    endif

    call m_Kp_get_nkmesh( nbz_mesh )
    kv3_wk = ( nbz_mesh(1) +1 )*( nbz_mesh(2) +1 )*( nbz_mesh(3) +1 ) *nspin

    allocate(kvtab(kv3_wk));  kvtab = 0
    call m_Kp_get_kptable_bxsf( kv3_wk, kvtab, 1 )

    lun = 360
    Do ispin=1, nspin, ndim_spinor
       if ( nspin == 1 .or. ndim_spinor == 2 ) then
          Open( unit=lun, file="./myband.bxsf", status="unknown", form="formatted" )
       else
          if ( ispin == 1 ) then
             Open( unit=lun, file="./myband_up.bxsf", status="unknown", &
                  &          form="formatted" )
          else
             Open( unit=lun, file="./myband_down.bxsf", status="unknown", &
                  &          form="formatted" )
          endif
       endif
       call print_header( lun );  call print_body( lun );   call print_tail( lun )
       close( lun )
    End Do

    deallocate(kvtab)
    if ( allocated(neordr_t) ) deallocate(neordr_t);

  contains

    subroutine print_header( lun )
      integer, intent(in) :: lun

      write( lun, '(A)') "BEGIN_INFO"
      write( lun, '(A)') "  #"
      write( lun, '(A)') "  # this is a Band-XCRYSDEN-Structure-File"
      write( lun, '(A)') "  #  aimed for Visualization of Fermi Surface"
      write( lun, '(A)') "  #"
      write( lun, '(A,F20.15)') "    Fermi Energy:",  0.0d0
      write( lun, '(A)') "END_INFO"
      write( lun, * )
      write( lun, '(A)') "BEGIN_BLOCK_BANDGRID_3D"
      write( lun, '(A)') "  from_phase/0"
      write( lun, '(A)') "  BEGIN_BANDGRID_3D_fermi"
      write( lun, '(I6)') neg_t
      write( lun, '(3I6)') nbz_mesh(1:3) +1
      write( lun, '(3F10.2)') 0.0, 0.0, 0.0
      write( lun, '(3F12.8)') rltv(1:3,1) /Bohr
      write( lun, '(3F12.8)') rltv(1:3,2) /Bohr
      write( lun, '(3F12.8)') rltv(1:3,3) /Bohr
    end subroutine print_header

    subroutine print_body( lun )
      integer, intent(in) :: lun
      integer ib, ik

      if ( ekmode == ON ) then
         if ( ndim_spinor == 1 ) then
            Do ib=1, neg_t
               write( lun, '(A,I8)') "    BAND: ", ib
               write( lun, '(5F15.8)') &
                    &   ( ( eko_ek( neordr_t(ib,kvtab(ik)),kvtab(ik) )-efermi ) &
                    &         *Hartree,  ik=ispin,kv3_wk,nspin )
            End Do
         else
            Do ib=1, neg_t
               write( lun, '(A,I8)') "    BAND: ", ib
               write( lun, '(5F15.8)') &
                    &   ( ( eko_ek( neordr_t(ib,kvtab(ik)),kvtab(ik) ) -efermi ) &
                    &          *Hartree, ik=ispin,kv3_wk,nspin )
            End Do
         endif
      endif
    end subroutine print_body

    subroutine print_tail( lun )
      integer, intent(in) :: lun
       write( lun, '(A)') "  END_BANDGRID_3D_fermi"
       write( lun, '(A)') "END_BLOCK_BANDGRID_3D"
    end subroutine print_tail

  end subroutine m_ESIO_wd_EigenValues_bxsf

  subroutine m_ESIO_wd_EigenValues_frmsf
    integer :: ik, ib,jb,ibo,jbo, neg_t
    integer :: lun, nbz_mesh(3), ispin
    integer :: kv3_wk

    integer, allocatable :: neordr_t(:,:), kvtab(:)
    real(kind=DP), parameter :: delta = 1.d-12

    if ( ekmode == ON ) then
       allocate(neordr_t(neg,kv3_ek));   neordr_t = 0
       neg_t = neg
       if(neg_is_enlarged) neg_t = neg -num_extra_bands

       do ik = 1, kv3_ek, ndim_spinor
          if (nspin == 1 .or. (nspin == 2 .and. mod(ik,2) == 1)) &
               & neordr_t(1:neg,ik) = (/(ib,ib=1,neg)/)
          do ib = 1, neg-1
             do jb = ib+1, neg
                ibo = neordr_t(ib,ik);  jbo = neordr_t(jb,ik)
                if ( eko_ek(jbo,ik) < eko_ek(ibo,ik)-delta ) then        ! MPI
                   neordr_t(jb,ik) = ibo;    neordr_t(ib,ik) = jbo
                end if
             end do
          end do
       end do
    else
       !write(*,*) "ekmode == Off is not supported "
       call phase_error_with_msg(6,"ekmode == Off is not supported ",__LINE__,__FILE__)
    endif

    call m_Kp_get_nkmesh( nbz_mesh )
    kv3_wk = ( nbz_mesh(1) )*( nbz_mesh(2) )*( nbz_mesh(3) ) *nspin

    allocate(kvtab(kv3_wk));  kvtab = 0
    call m_Kp_get_kptable_bxsf( kv3_wk, kvtab, 0 )

    lun = 360
    Do ispin=1, nspin, ndim_spinor
       if ( nspin == 1 .or. ndim_spinor == 2 ) then
          Open( unit=lun, file="./myband.frmsf", status="unknown", form="formatted" )
       else
          if ( ispin == 1 ) then
             Open( unit=lun, file="./myband_up.frmsf", status="unknown", &
                  &          form="formatted" )
          else
             Open( unit=lun, file="./myband_down.frmsf", status="unknown", &
                  &          form="formatted" )
          endif
       endif
       call print_header( lun );  call print_body( lun );   call print_tail( lun )
       close( lun )
    End Do

    deallocate(kvtab)
    if ( allocated(neordr_t) ) deallocate(neordr_t);

  contains

    subroutine print_header( lun )
      integer, intent(in) :: lun

      write( lun, '(3I6)') nbz_mesh(1:3)
      write( lun, '(I6)') 1
      write( lun, '(I6)') neg_t
      write( lun, '(3F12.8)') rltv(1:3,1) /Bohr
      write( lun, '(3F12.8)') rltv(1:3,2) /Bohr
      write( lun, '(3F12.8)') rltv(1:3,3) /Bohr
    end subroutine print_header

    subroutine print_body( lun )
      integer, intent(in) :: lun
      integer ib, ik

      if ( ekmode == ON ) then
         if ( ndim_spinor == 1 ) then
            Do ib=1, neg_t
               write( lun, '(5F15.8)') &
                    &   ( ( eko_ek( neordr_t(ib,kvtab(ik)),kvtab(ik) )-efermi ) &
                    &         *Hartree,  ik=ispin,kv3_wk,nspin )
            End Do
         else
            Do ib=1, neg_t
               write( lun, '(5F15.8)') &
                    &   ( ( eko_ek( neordr_t(ib,kvtab(ik)),kvtab(ik) ) -efermi ) &
                    &          *Hartree, ik=ispin,kv3_wk,nspin )
            End Do
         endif
      endif
    end subroutine print_body

    subroutine print_tail( lun )
      integer, intent(in) :: lun
    end subroutine print_tail

  end subroutine m_ESIO_wd_EigenValues_frmsf

  subroutine m_ESIO_wd_wfn_parity
    use m_Control_Parameters, only : eval_parity_on_fftmesh
    use m_Const_Parameters,   only : ELECTRON, INVERSE, INITIAL, FIXED_CHARGE, zi
    use m_Crystal_Structure,  only : tau, nopr, op, Global_Quantz_Axis_Fixed
    use m_Ionic_System,       only : pos
    use m_PseudoPotential,    only : ltp
    use m_PlaneWaveBasisSet,  only : kg, igf, ngpt_l, nbase_gamma
    use m_Files,              only : nfout
    use m_FFT,                only : m_FFT_alloc_WF_work, m_FFT_dealloc_WF_work, &
         &                           m_FFT_WF

    integer :: iopr_inv
    real(kind=DP) :: inv_center(3)

    integer, allocatable :: ngpt_inv(:), igf_shift(:,:,:,:)
    real(kind=DP), allocatable :: eko_wk(:,:), moment_wk(:,:,:)
    complex(kind=CMPLDP), allocatable, target :: parity_bands(:,:)

    call chkif_system_has_inv_symm
    if ( iopr_inv == 0 ) then
       !write(*,*) "The system does not have inversion symmetry."
       call phase_error_with_msg(nfout,"The system does not have inversion symmetry.",__LINE__,__FILE__)
    endif

    allocate( eko_wk(neg,kv3/ndim_spinor) );  eko_wk = 0.0d0
    allocate( moment_wk(ndim_magmom,neg,kv3/ndim_spinor) ); moment_wk = 0.0d0
    call set_eko_wk
    call calc_moment_wk

    allocate( parity_bands( ista_e:iend_e, ista_k:iend_k ) )
    parity_bands = 0.0d0

    if ( eval_parity_on_fftmesh == YES ) then
       call softpart_eval_on_Rspace
    else
       allocate( ngpt_inv(kg) );  ngpt_inv = 0
       call set_ngpt_inv( iopr_inv, ngpt_inv )
       call set_igf_with_shift
       call softpart_eval_on_Gspace
    endif

    call hardpart
    call print_parity

    deallocate( parity_bands );
    deallocate( eko_wk );    deallocate( moment_wk )
    if ( allocated( ngpt_inv ) ) deallocate( ngpt_inv )
    if ( allocated( igf_shift ) ) deallocate( igf_shift )

  contains

    subroutine set_eko_wk
      integer :: ik, ie, ierr, iksnl
      real(kind=DP), allocatable :: eko_mpi(:,:), occ_mpi(:,:)

      if ( noncol ) then
         do ik = 1, kv3, ndim_spinor
            if (map_k(ik) /= myrank_k) cycle
            iksnl = (ik-1)/nspin + 1
            do ie = 1, neg
               if (map_e(ie) /= myrank_e) cycle
               eko_wk(ie,iksnl) = eko_l(map_z(ie),ik)
            end do
         end do
      else
         do ik = 1, kv3
            if (map_k(ik) /= myrank_k) cycle
            do ie = 1, neg
               if (map_e(ie) /= myrank_e) cycle
               eko_wk(ie,ik) = eko_l(map_z(ie),ik)
            end do
         end do
      endif

      if ( npes > 1 ) then
         allocate( eko_mpi(neg,kv3/ndim_spinor) ); eko_mpi = 0.0d0
         call mpi_allreduce( eko_wk, eko_mpi, neg*kv3/ndim_spinor, mpi_double_precision,&
              &              mpi_sum, MPI_CommGroup, ierr )
         eko_wk = eko_mpi;
         deallocate( eko_mpi );
      endif
    end subroutine set_eko_wk

    subroutine calc_moment_wk
      integer :: ik, iksnl, ib, ie, is, ig
      integer :: ia, it, lmt1, lmt2, p1, p2
      complex(kind=CMPLDP) :: z1, z2, z3, z4, ztmp1, ztmp2, ztmp3, ztmp4
      real(kind=DP) :: csum(ndim_magmom)
      real(kind=DP), allocatable :: mom_mpi(:,:,:)

      if ( noncol ) then
         Do ik=ista_k, iend_k, ndim_spinor
            iksnl = (ik-1)/nspin + 1
            Do ie =ista_e, iend_e, istep_e
               ib = map_z(ie)
               csum = 0.0d0
               Do ig=1, iba(ik)
                  z1 = cmplx( zaj_l(ig,ib,ik,1),   zaj_l(ig,ib,ik,2) )
                  z2 = cmplx( zaj_l(ig,ib,ik+1,1), zaj_l(ig,ib,ik+1,2) )
                  csum(1) = csum(1) +conjg(z1)*z1 +conjg(z2)*z2
                  csum(2) = csum(2) +conjg(z1)*z2 +conjg(z2)*z1
                  csum(3) = csum(3) +( -conjg(z1)*z2 +conjg(z2)*z1 )*zi
                  csum(4) = csum(4) +conjg(z1)*z1 -conjg(z2)*z2
               End Do
               !
               do ia = 1, natm
                  it = ityp(ia)
                  do lmt1 = 1, ilmt(it)
                     p1 = lmta(lmt1,ia)
                     do lmt2 = 1, ilmt(it)
                        p2 = lmta(lmt2,ia)
                        z1 = cmplx( fsr_l(ib,p1,ik),   fsi_l(ib,p1,ik) )
                        z2 = cmplx( fsr_l(ib,p1,ik+1), fsi_l(ib,p1,ik+1) )
                        z3 = cmplx( fsr_l(ib,p2,ik),   fsi_l(ib,p2,ik) )
                        z4 = cmplx( fsr_l(ib,p2,ik+1), fsi_l(ib,p2,ik+1) )
                        ztmp1 = conjg(z1)*z3 +conjg(z2)*z4
                        ztmp2 = conjg(z1)*z4 +conjg(z2)*z3
                        ztmp3 = ( -conjg(z1)*z4 +conjg(z2)*z3 ) *zi
                        ztmp4 = conjg(z1)*z3 -conjg(z2)*z4
                        csum(1) = csum(1) +q(lmt1,lmt2,it) *ztmp1
                        csum(2) = csum(2) +q(lmt1,lmt2,it) *ztmp2
                        csum(3) = csum(3) +q(lmt1,lmt2,it) *ztmp3
                        csum(4) = csum(4) +q(lmt1,lmt2,it) *ztmp4
                     End do
                  End do
               End do
               moment_wk(1:4,ie,iksnl) = csum(1:4)
            End Do
         End Do
      else
      endif
      if ( npes > 1 ) then
         allocate( mom_mpi(ndim_magmom,neg,kv3/ndim_spinor) ); mom_mpi = 0.0d0
         call mpi_allreduce( moment_wk, mom_mpi, ndim_magmom*neg*kv3/ndim_spinor, &
              &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
         moment_wk = mom_mpi;
         deallocate( mom_mpi );
      endif
    end subroutine calc_moment_wk

    subroutine print_parity
      integer :: ierr, ik, ib, is, lun, iksnl
      real(kind=DP) :: c1, cmul(ndim_spinor)
      complex(kind=CMPLDP) :: z1
      complex(kind=CMPLDP), allocatable, target :: work(:,:)
      complex(kind=CMPLDP), allocatable, target :: work2(:,:)
      complex(kind=CMPLDP), pointer :: wk_array(:,:)
      logical, save :: First = .true.

      lun = 1000

      if ( npes > 1 ) then
         allocate( work(neg,kv3) );  work  = 0.0d0
         allocate( work2(neg,kv3) ); work2 = 0.0d0
         work(ista_e:iend_e,ista_k:iend_k) = parity_bands(ista_e:iend_e,ista_k:iend_k)
         call mpi_allreduce( work, work2, 2*neg*kv3, mpi_double_precision, &
              &              mpi_sum, MPI_CommGroup, ierr )
         deallocate(work)
         wk_array => work2
      else
         wk_array => parity_bands
      endif

      if ( mype == 0 ) then
         if ( ekmode == ON ) then
            if ( nk_in_the_process == 1 ) then
               open( lun, file="WF_parity", status="unknown", form="formatted")
               write(lun,'(A)') '####### Parity Eigenvalues #######'
            else
               open( lun, file="WF_parity", status="unknown", form="formatted", &
                    &     position = 'append' )
            endif

            Do ik=1, kv3, ndim_spinor
               iksnl = (ik-1)/ndim_spinor + 1
               write(lun,'(A,I6,A,3F16.8,A)') 'ik = ', &
                    &       ik +nk_in_the_process -1, " ( ", vkxyz(ik,1:3,BUCS), " )"
               write(lun,'(A,5X,A,12X,A)') &
                    &     "band ", "energy", "parity"

               cmul = 1.0d0
               Do ib=1, neg
                  if ( noncol ) then
                     c1 = moment_wk(2,ib,iksnl) *Global_Quantz_Axis_Fixed(1) &
                          & +moment_wk(3,ib,iksnl) *Global_Quantz_Axis_Fixed(2) &
                          & +moment_wk(4,ib,iksnl) *Global_Quantz_Axis_Fixed(3)
                  else
                     c1 = 1.0d0
                  endif

                  z1 = 0.0d0
                  Do is=1, ndim_spinor
                     z1 = z1 +wk_array(ib,ik+is-1)
                  End Do

                  if ( c1 > 0.0 ) then
                     write(lun,'(I5,2X,F16.12,2X,F10.5)') &
                          &     ib, ( eko_wk(ib,iksnl) -efermi )*Hartree, real(z1)
                  else
                     write(lun,'(I5,2X,F16.12,2X,10X,F10.5)') &
                          &     ib, ( eko_wk(ib,iksnl) -efermi )*Hartree, real(z1)
                  endif

                  if ( eko_wk(ib,iksnl) <= efermi ) then
                     if ( noncol ) then
                        if ( c1 > 0.0 ) then
                           cmul(1) = cmul(1) *real(z1)
                        else
                           cmul(2) = cmul(2) *real(z1)
                        endif
                     else
                        cmul = cmul *real(z1)
                     endif
                  endif
               End Do
               write(lun,'(A)')
               write(lun,'(A)') " product of eigenvalues (below EF) : "
               if ( noncol ) then
                  write(lun,'(25X,2F10.5)') cmul
               else
                  write(lun,'(25X,F10.5)') cmul
               endif
               write(lun,'(A)')
            End Do

            close(lun)

         else
            if ( icond == INITIAL .or. icond == FIXED_CHARGE ) then
               open( lun, file="WF_parity", status="unknown", form="formatted")
               write(lun,'(A)') '####### Parity Eigenvalues #######'
            else
               open( lun, file="WF_parity", status="unknown", form="formatted", &
                    &     position = 'append' )
            endif

            Do ik=1, kv3, ndim_spinor
               iksnl = (ik-1)/ndim_spinor + 1
               write(lun,'(A,I6,A,3F16.8,A)') 'ik = ', &
                    &       ik, " ( ", vkxyz(ik,1:3,BUCS), " )"
               write(lun,'(A,5X,A,12X,A)') &
                    &     "band ", "energy", "parity"

               cmul = 1.0d0
               Do ib=1, neg
                  c1 = moment_wk(2,ib,iksnl) *Global_Quantz_Axis_Fixed(1) &
                       & +moment_wk(3,ib,iksnl) *Global_Quantz_Axis_Fixed(2) &
                       & +moment_wk(4,ib,iksnl) *Global_Quantz_Axis_Fixed(3)

                  z1 = 0.0d0
                  Do is=1, ndim_spinor
                     z1 = z1 +wk_array(ib,ik+is-1)
                  End Do

                  if ( c1 > 0.0 ) then
                     write(lun,'(I5,2X,F16.12,2X,F10.5)') &
                          &     ib, ( eko_wk(ib,iksnl) -efermi )*Hartree, real(z1)
                  else
                     write(lun,'(I5,2X,F16.12,2X,10X,F10.5)') &
                          &     ib, ( eko_wk(ib,iksnl) -efermi )*Hartree, real(z1)
                  endif

                  if ( eko_wk(ib,iksnl) <= efermi ) then
                     if ( noncol ) then
                        if ( c1 > 0.0 ) then
                           cmul(1) = cmul(1) *real(z1)
                        else
                           cmul(2) = cmul(2) *real(z1)
                        endif
                     else
                        cmul = cmul *real(z1)
                     endif
                  endif
               End Do
               write(lun,'(A)')
               write(lun,'(A)') " product of eigenvalues (below EF) : "
               if ( noncol ) then
                  write(lun,'(25X,2F10.5)') cmul
               else
                  write(lun,'(25X,F10.5)') cmul
               endif
               write(lun,'(A)')
            End Do
            close(lun)
         end if
      endif

      if ( allocated(work2) ) deallocate( work2 )
    end subroutine print_parity

    subroutine chkif_system_has_inv_symm
      integer :: iopr
      real(kind=DP) :: c1

      iopr_inv = 0;     inv_center = 0.0d0
      Do iopr=1, nopr
         c1 = op(1,1,iopr)+op(2,2,iopr)+op(3,3,iopr)
         if ( abs( c1 +3.d0 ) < 1.0D-3 ) then
            iopr_inv = iopr
            inv_center(:) = tau(:,iopr,BUCS) /2.0d0
         endif
      End Do
    end subroutine chkif_system_has_inv_symm

    subroutine softpart_eval_on_Gspace
      integer :: ib1, ik, is, ig, ip1, i, nshift(3)
      real(kind=DP) :: vec(3), vec2(3), c1
      real(kind=DP) :: dnorm, c2r, c2i
      complex(kind=CMPLDP) :: zsum

      real(kind=DP), allocatable :: wf1(:), wf2(:), zfcos(:), zfsin(:)

      dnorm = sqrt( inv_center(1)**2 +inv_center(2)**2 +inv_center(3)**2 )
      if ( dnorm > 1.0D-3 ) then
         allocate( zfcos(kg1) );    allocate( zfsin(kg1) )
      endif

      allocate( wf1(nfft) ); wf1 = 0.0d0;
      allocate( wf2(nfft) ); wf2 = 0.0d0;

      Do ik = 1, kv3
         if ( map_k(ik) /= myrank_k ) cycle! MPI
         vec(:) = vkxyz(ik,:,BUCS)*2.d0
         vec2(:) = vec(:) -nint(vec(:))
         c1 = sqrt( vec2(1)**2 +vec2(2)**2 +vec2(3)**2 )
         if ( c1 > 1.0D-3 ) cycle

         nshift(:) = -nint(vec(:))

         if ( dnorm > 1.0D-3 ) then
            zfcos = 0.0d0;  zfsin = 0.0d0
            do ig = 1, iba(ik)
               ip1 = nbase(ig,ik)
               vec(:) = vkxyz(ik,:,BUCS) + ngabc(ip1,:)
               c1 = vec(1) *inv_center(1) +vec(2) *inv_center(2) +vec(3) *inv_center(3)
               c1 = 2.d0 *c1 *PAI2
               zfcos(ig) = cos(c1);    zfsin(ig) = sin(c1)
            End Do
         endif
         Do ib1 = ista_e, iend_e, istep_e     ! MPI
            if ( dnorm > 1.0D-3 ) then
               call map_WF_on_fftmesh( ik, ik, ik, zaj_l(:,map_z(ib1),ik,:), wf1, &
                    &                  igf_shift(:,nshift(1),nshift(2),nshift(3)), &
                    &                  zfcos, zfsin )
            else
               call map_WF_on_fftmesh( ik, ik, ik, zaj_l(:,map_z(ib1),ik,:), wf1, &
                    &                  igf_shift(:,nshift(1),nshift(2),nshift(3)) )
            endif
            call map_WF_on_fftmesh( ik, ik, ik, zaj_l(:,map_z(ib1),ik,:), wf2, &
                 &                  igf )
            zsum = 0.0d0
            Do i=1, nfft, 2
               c2r = wf2(i) *wf1(i)  +wf2(i+1) *wf1(i+1)
               c2i = wf2(i) *wf1(i+1)-wf2(i+1) *wf1(i)
               zsum = zsum +cmplx(c2r,c2i)
            End Do
            parity_bands( ib1,ik ) = zsum
         End Do
      End Do

      deallocate( wf1 );    deallocate( wf2 )
      if ( allocated(zfcos) ) deallocate( zfcos )
      if ( allocated(zfsin) ) deallocate( zfsin )

    end subroutine softpart_eval_on_Gspace

    subroutine softpart_eval_on_Rspace
      integer :: ib1, ik, is, ig
      integer :: ix1, iy1, iz1, ix2, iy2, iz2, nxy, ip1, ip2
      real(kind=DP) :: c1, c2r, c2i, c3r, c3i, dnorm
      real(kind=DP) :: da(3), damul, vec(3), cos1, sin1, vec2(3)
      complex(kind=CMPLDP) :: zsum
      real(kind=DP), allocatable :: wf1(:), zfcos(:), zfsin(:)

      nxy = fft_box_size_WF(1,0) *fft_box_size_WF(2,0)
      da(1:3) = 1.d0/fft_box_size_WF(1:3,1)
      damul = da(1) *da(2) *da(3)

      allocate( wf1(nfft) ); wf1 = 0.0d0;
      call m_FFT_alloc_WF_work()
      !
      dnorm = sqrt( inv_center(1)**2 +inv_center(2)**2 +inv_center(3)**2 )
      if ( dnorm > 1.0D-3 ) then
         allocate( zfcos(kg1) );    allocate( zfsin(kg1) )
      endif

      Do ik = 1, kv3
         if ( map_k(ik) /= myrank_k ) cycle! MPI
         !
         vec(:) = vkxyz(ik,:,BUCS)*2.d0
         vec2(:) = vec(:) -nint(vec(:))
         c1 = sqrt( vec2(1)**2 +vec2(2)**2 +vec2(3)**2 )
         if ( c1 > 1.0D-3 ) cycle

         if ( dnorm > 1.0D-3 ) then
            zfcos = 0.0d0;  zfsin = 0.0d0
            do ig = 1, iba(ik)
               ip1 = nbase(ig,ik)
               vec(:) = vkxyz(ik,:,BUCS) + ngabc(ip1,:)
               c1 = vec(1) *inv_center(1) +vec(2) *inv_center(2) +vec(3) *inv_center(3)
               c1 = c1 *PAI2
               zfcos(ig) = cos(c1);    zfsin(ig) = sin(c1)
            End Do
         endif

         Do ib1 = ista_e, iend_e, istep_e     ! MPI
            if ( dnorm > 1.0D-3 ) then
               call map_WF_on_fftmesh( ik, ik, ik, zaj_l(:,map_z(ib1),ik,:), &
                    &                  wf1, igf, zfcos, zfsin )
               call m_FFT_WF( ELECTRON, nfout, wf1, INVERSE, ON )
            else
               call m_ES_WF_in_Rspace( ik, ib1, wf1 )  ! unk(r)
            endif

            ! unk(-r)
            zsum = 0.0d0

            Do ix1=1, fft_box_size_WF(1,1)
               Do iy1=1, fft_box_size_WF(2,1)
                  Do iz1=1, fft_box_size_WF(3,1)
                     ip1 = (iz1-1) *nxy +(iy1-1) *fft_box_size_WF(1,0) +ix1
                     ix2 = mod( 1-ix1 +fft_box_size_WF(1,1),fft_box_size_WF(1,1) )+1
                     iy2 = mod( 1-iy1 +fft_box_size_WF(2,1),fft_box_size_WF(2,1) )+1
                     iz2 = mod( 1-iz1 +fft_box_size_WF(3,1),fft_box_size_WF(3,1) )+1
                     ip2 = (iz2-1) *nxy +(iy2-1) *fft_box_size_WF(1,0) +ix2
! exp( -2*i *k*r )
                     c1 = vkxyz(ik,1,BUCS) *( ix1 -1 )*da(1) &
                          & +vkxyz(ik,2,BUCS) *( iy1 -1 )*da(2) &
                          & +vkxyz(ik,3,BUCS) *( iz1 -1 )*da(3)
                     cos1 = cos( -2.0 *c1 *PAI2 )
                     sin1 = sin( -2.0 *c1 *PAI2 )

                     c2r = wf1( 2*ip1-1 )*wf1( 2*ip2-1 ) +wf1( 2*ip1 )*wf1( 2*ip2 )
                     c2i = wf1( 2*ip1-1 )*wf1( 2*ip2 )   -wf1( 2*ip1 )*wf1( 2*ip2-1 )
                     !
                     c3r = c2r *cos1 -c2i *sin1
                     c3i = c2r *sin1 +c2i *cos1
                     zsum = zsum +cmplx(c3r,c3i)
                  end Do
               end Do
            end Do
            parity_bands( ib1,ik ) = zsum *damul
         End Do
      End Do

      call m_FFT_dealloc_WF_work()

      deallocate( wf1 )
      if ( allocated(zfcos) ) deallocate( zfcos )
      if ( allocated(zfsin) ) deallocate( zfsin )

    end subroutine softpart_eval_on_Rspace

    subroutine hardpart
      integer :: ia, ja, it, ib1, ik, i
      integer :: lmt1, lmt2, p1, p2, il2
      real(kind=DP) :: dx, dy, dz, c1r, c1i, dist, fac, pos_a_inv(3)
      real(kind=DP) :: ph, c3r, c3i, vec(3), vec2(3)
      complex(kind=CMPLDP) :: zsum

      integer, allocatable :: inv_ia(:)
      real(kind=DP), allocatable :: pos_t(:,:), wkcos(:), wksin(:)

      allocate( pos_t(natm,3) )
      allocate( inv_ia(natm) ); inv_ia = 0

      pos_t(1:natm,1:3) = pos(1:natm,1:3)
      do ia = 1, natm
         do i = 1, 3
            pos_t(ia,i) = pos_t(ia,i) - floor(pos_t(ia,i))
         end do
      end do
      Do ia=1, natm
         pos_a_inv(:) = -pos_t(ia,:) +inv_center(:) *2.d0
         pos_a_inv(:) = pos_a_inv(:) - floor(pos_a_inv(:))
         loop_ja: Do ja=1, natm
            dx = pos_t(ja,1) -pos_a_inv(1)
            dy = pos_t(ja,2) -pos_a_inv(2)
            dz = pos_t(ja,3) -pos_a_inv(3)
            dist = sqrt( dx**2 +dy**2 +dz**2 )
            if ( dist < 1.0D-5 ) then
               inv_ia(ia) = ja;     exit loop_ja
            endif
         End Do loop_ja
      End Do
      !
      Do ia=1, natm
         if ( inv_ia(ia) == 0 ) then
            write(*,*) "Inverted atom is not found for atom", ia
            call phase_error_with_msg(nfout,"Inverted atom is not found ",__LINE__,__FILE__)
         endif
      End Do

      allocate( wkcos(natm) );    allocate( wksin(natm) )

      Do ik = 1, kv3
         if ( map_k(ik) /= myrank_k ) cycle! MPI

         vec(:) = vkxyz(ik,:,BUCS)*2.d0
         vec2(:) = vec(:) -nint(vec(:))
         c1r = sqrt( vec2(1)**2 +vec2(2)**2 +vec2(3)**2 )
         if ( c1r > 1.0D-3 ) cycle
         !
         Do ia=1, natm
            ph = vkxyz(ik,1,BUCS) *( inv_center(1) -pos_t(ia,1) ) &
                 & +vkxyz(ik,2,BUCS) *( inv_center(2) -pos_t(ia,2) ) &
                 & +vkxyz(ik,3,BUCS) *( inv_center(3) -pos_t(ia,3) )
            c1r = 2.0 *ph *PAI2
            wkcos(ia) = cos(c1r);  wksin(ia) = sin(c1r)
         End Do

         Do ib1 = ista_e, iend_e, istep_e     ! MPI
            zsum = 0.0d0

            Do ia=1, natm
               it = ityp(ia)
               do lmt1 = 1, ilmt(it)
                  p1 = lmta(lmt1,ia)
                  if (k_symmetry(ik) == GAMMA) then
                     do lmt2 = 1, ilmt(it)
                        il2 = ltp(lmt2,it) -1
                        fac = 1.0d0
                        if ( mod(il2, 2) == 1 ) fac =-fac

                        p2 = lmta(lmt2,inv_ia(ia))
                        c1r = fsr_l(map_z(ib1),p1,ik) *fsr_l(map_z(ib1),p2,ik)

                        zsum = zsum +fac *c1r *q(lmt1,lmt2,it)
                     end do

                  else
                     do lmt2 = 1, ilmt(it)
                        il2 = ltp(lmt2,it) -1
                        fac = 1.0d0
                        if ( mod(il2, 2) == 1 ) fac =-fac

                        p2 = lmta(lmt2,inv_ia(ia))
                        c1r = fsr_l(map_z(ib1),p1,ik) *fsr_l(map_z(ib1),p2,ik) &
                             & +fsi_l(map_z(ib1),p1,ik) *fsi_l(map_z(ib1),p2,ik)
                        c1i = fsr_l(map_z(ib1),p1,ik) *fsi_l(map_z(ib1),p2,ik) &
                             &  -fsi_l(map_z(ib1),p1,ik) *fsr_l(map_z(ib1),p2,ik)

                        c3r = c1r *wkcos(ia) -c1i *wksin(ia)
                        c3i = c1r *wksin(ia) +c1i *wkcos(ia)

                        zsum = zsum +fac *cmplx(c3r,c3i) *q(lmt1,lmt2,it)
                     end do
                  endif
               End do
            End Do
            parity_bands(ib1,ik) = parity_bands(ib1,ik) +zsum
         End Do
      End Do
      deallocate( pos_t );  deallocate( inv_ia )
      deallocate( wkcos );  deallocate( wksin )
    end subroutine hardpart

    subroutine set_igf_with_shift
      integer :: id, i, i2
      integer :: igf1, igf2, igf3, nx, ny, nz
      integer :: nxmin, nymin, nzmin, nxmax, nymax, nzmax
      !
      nxmin = -1;   nxmax = 1
      nymin = -1;   nymax = 1
      nzmin = -1;   nzmax = 1

      allocate( igf_shift( kg, nxmin:nxmax, nymin:nymax, nzmin:nzmax ) )
      igf_shift = 0
      !
      id = fft_box_size_WF(1,0)

      Do nx=nxmin, nxmax
         Do ny=nymin, nymax
            Do nz=nzmin, nzmax
               do i = 1, kg
                  i2 = ngpt_inv(i)
                  igf1 = ngabc(i2,1) + nx +1
                  igf2 = ngabc(i2,2) + ny +1
                  igf3 = ngabc(i2,3) + nz +1
                  !
                  if ( igf1 <= 0 ) igf1 = igf1 + fft_box_size_WF(1,1)
                  if ( igf2 <= 0 ) igf2 = igf2 + fft_box_size_WF(2,1)
                  if ( igf3 <= 0 ) igf3 = igf3 + fft_box_size_WF(3,1)

                  igf_shift(i,nx,ny,nz) = igf1 + (igf2-1)*id &
                       &                 + (igf3-1)*id*fft_box_size_WF(2,0)
               enddo
            End Do
         End Do
      ENd Do
    end subroutine set_igf_with_shift

    subroutine map_WF_on_fftmesh( k1, k2, ik, psi_l, bfft, igf_in, zfcos, zfsin )
      integer, intent(in) :: k1, k2, ik
      integer, intent(in) :: igf_in( kg )
      real(kind=DP), intent(in) :: psi_l( kg1, 1, k1:k2, kimg )
      real(kind=DP), intent(in), optional :: zfcos( kg1 ), zfsin( kg1 )
      real(kind=DP), intent(inout) :: bfft( nfft )

      integer :: i,i1,ri, j, i2, ii
      real(kind=DP) :: c1, c2

      bfft = 0.d0
      if(k_symmetry(ik) == GAMMA) then
         if(kimg == 1) then
            i1 = igf_in(1)
            bfft(i1) = psi_l(1,1,ik,1)
#ifdef NEC_TUNE_SMP
            !CDIR NODEP
#endif
            if ( present(zfcos) ) then
               do ii = 2, iba(ik)
!!$             i = nbase(ii,1)
                  i = nbase(ii,ik)
                  i1 = igf_in(i)
                  bfft(i1) = psi_l(ii,1,ik,1) *zfcos(ii)
                  j = nbase_gamma(ii,2)
                  i2 = igf_in(j)
                  bfft(i2) = psi_l(ii,1,ik,1) *zfcos(ii)
               end do
            else
               do ii = 2, iba(ik)
!!$             i = nbase(ii,1)
                  i = nbase(ii,ik)
                  i1 = igf_in(i)
                  bfft(i1) = psi_l(ii,1,ik,1)
                  j = nbase_gamma(ii,2)
                  i2 = igf_in(j)
                  bfft(i2) = psi_l(ii,1,ik,1)
               end do
            endif

         else if(kimg == 2) then
            i1 = 2*igf_in(1) - 1
            bfft(i1)   = psi_l(1,1,ik,1)
            bfft(i1+1) = psi_l(1,1,ik,2)
#ifdef NEC_TUNE_SMP
            !CDIR NODEP
#endif
            if ( present(zfcos) ) then
               do ii = 2, iba(ik)
!!$             i = nbase(ii,1)
                  i = nbase(ii,ik)
                  i1 = 2*igf_in(i)-1
                  c1 = psi_l(ii,1,ik,1) *zfcos(ii) -psi_l(ii,1,ik,2) *zfsin(ii)
                  c2 = psi_l(ii,1,ik,1) *zfsin(ii) +psi_l(ii,1,ik,2) *zfcos(ii)
                  bfft(i1  ) = c1;       bfft(i1+1) = c2

                  j = nbase_gamma(ii,2)
                  i2 = 2*igf_in(j)-1
                  c1 = psi_l(ii,1,ik,1) *zfcos(ii) -psi_l(ii,1,ik,2) *zfsin(ii)
                  c2 = psi_l(ii,1,ik,1) *zfsin(ii) +psi_l(ii,1,ik,2) *zfcos(ii)
                  bfft(i2  ) = c1;       bfft(i2+1) = -c2
               end do
            else
               do ii = 2, iba(ik)
!!$             i = nbase(ii,1)
                  i = nbase(ii,ik)
                  i1 = 2*igf_in(i)-1
                  c1 = psi_l(ii,1,ik,1)
                  c2 = psi_l(ii,1,ik,2)
                  bfft(i1  ) = c1;       bfft(i1+1) = c2

                  j = nbase_gamma(ii,2)
                  i2 = 2*igf_in(j)-1
                  c1 = psi_l(ii,1,ik,1)
                  c2 = psi_l(ii,1,ik,2)
                  bfft(i2  ) = c1;       bfft(i2+1) = -c2
               end do
            end if
         endif
      else
#ifdef NEC_TUNE_SMP
         !CDIR NOLOOPCHG
#endif
         if ( present(zfcos) ) then
            if ( kimg == 1 ) then
               do i = 1, iba(ik)
                  i1 = igf_in(nbase(i,ik))
                  c1 = psi_l(i,1,ik,1) *zfcos(i)
                  bfft(i1) = c1
               end do
            else
               do i = 1, iba(ik)
                  i1 = kimg*igf_in(nbase(i,ik)) -1
                  c1 = psi_l(i,1,ik,1) *zfcos(i) -psi_l(i,1,ik,2) *zfsin(i)
                  c2 = psi_l(i,1,ik,1) *zfsin(i) +psi_l(i,1,ik,2) *zfcos(i)
                  bfft(i1) = c1;        bfft(i1+1) = c2
               end do
            endif
         else
            if ( kimg == 1 ) then
               do i = 1, iba(ik)
                  i1 = igf_in(nbase(i,ik))
                  c1 = psi_l(i,1,ik,1)
                  bfft(i1) = c1;
               end do
            else
               do i = 1, iba(ik)
                  i1 = kimg*igf_in(nbase(i,ik)) -1
                  c1 = psi_l(i,1,ik,1)
                  c2 = psi_l(i,1,ik,2)
                  bfft(i1) = c1;        bfft(i1+1) = c2
               end do
            end if
         endif
      end if
    end subroutine map_WF_on_fftmesh

    subroutine set_ngpt_inv( iopr_inv, ngpt_inv )
      integer, intent(in) :: iopr_inv
      integer, intent(out) :: ngpt_inv(kg)

      integer :: iopr, i, ierr
      integer, allocatable :: ngpt_tmp(:)

      iopr = iopr_inv

      if (npes > 1) then
         allocate(ngpt_tmp(kgp));   ngpt_tmp = 0
         do i = ista_kngp, iend_kngp
            ngpt_tmp(i) = ngpt_l(i,iopr)
         end do
         call mpi_allreduce( MPI_IN_PLACE, ngpt_tmp, kgp, mpi_integer, &
              &              mpi_sum, MPI_CommGroup, ierr )
         ngpt_inv(1:kg) = ngpt_tmp(1:kg)
         deallocate( ngpt_tmp )
      else
         do i = 1, kg
            ngpt_inv(i) = ngpt_l(i,iopr)
         end do
      end if
    end subroutine set_ngpt_inv

  end subroutine m_ESIO_wd_wfn_parity


end module m_ES_IO
