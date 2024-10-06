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

#define _USE_ALLREDUCE_IN_WD_WFS_

!
module m_ES_IO
! $Id: m_ES_IO.F90 635 2021-02-26 07:16:10Z jkoga $
  use m_Electronic_Structure, only : zaj_l,neordr,nrvf_ordr,eko_l,occup_l,efermi,efermi_spin,totch&
       &                            ,vnlph_l,vlhxc_l,eko_ek, zaj_l_prev, metalic_system, vbm
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
       &                           , nis_g1k,nie_g1k,map_rank_gek, myrank_ke, mpi_ke_world  &
       &                           , ista_spin, iend_spin, map_s, myrank_spin
  use m_IterationNumbers,     only : nk_in_the_process, nk_converged, nkgroup &
       &                           , first_kpoint_in_this_job, iteration_ionic, iteration_electronic
  use m_FFT,                  only : fft_box_size_WF,nfft
  use m_Crystal_Structure,    only : altv, sw_fix_total_spin, altv_refcell
  use m_Ionic_System,         only : natm,natm2,iatomn,m_IS_pack_all_ions_in_uc
  use m_PseudoPotential,      only : ival
  use m_Crystal_Structure,    only : univol, rltv
  use m_Parallelization,      only : neg_g,ista_g1k_prev,iend_g1k_prev

! ===================================== added by K. Tagami ============= 11.0
  use m_Control_Parameters,    only : ndim_spinor, noncol, &
       &                              previous_nspin_collinear, &
       &                              previous_nband_collinear
! ====================================================================== 11.0
  use m_ErrorMessages,        only : EOF_REACHED
  use m_Parallelization,      only : mpi_kg_world, mpi_ge_world

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

  use m_ES_initialWF, only : m_ESIW_by_randomnumbers0_3D
  use m_ES_ortho, only : m_ES_modified_gram_schmidt
  use m_ES_nonlocal, only : m_ES_betar_dot_WFs_3D
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
! 17.  m_ESIO_wd_WFs_3D             <-(WriteDownData_onto_Files)
! 18.  m_ESIO_rd_WFs_dp_3D
! 19.  m_ESIO_wd_WFs_dp_3D

contains
  subroutine m_ESIO_rd_EigenValues_etc(nfout,nfcntn_bin,F_CNTN_BIN_partitioned)
   use m_Parallelization,     only : mpi_ke_world, mpi_kg_world, mpi_chg_world

    integer, intent(in) :: nfout, nfcntn_bin
    logical, intent(in) :: F_CNTN_BIN_partitioned
    integer  :: ik, ie, ispin
    integer  :: ie_from
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
          do ie = 1, np_e
             eko_l(ie,ik)   = e1_wk(ie,ik-ista_k+1)
             occup_l(ie,ik) = e2_wk(ie,ik-ista_k+1)
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
       do ik = ispin, kv3-nspin+ispin, nspin
       !do ik = ista_k, iend_k                              ! MPI
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
       do ik = ispin, kv3-nspin+ispin, nspin
       !do ik = ista_k, iend_k                              ! MPI
          if(map_k(ik) /= myrank_k) cycle                 ! MPI
          do ie = 1, np_e                          ! MPI
             ie_from = neordr(neg_g(ie),ik)
             if( ie_from > neg_previous ) cycle
              eko_l(ie,ik) = e1_wk(ie_from,ik)            ! MPI
              occup_l(ie,ik) = e2_wk(ie_from,ik)            ! MPI
          end do                                           ! MPI
          if(neg_previous < neg) then
            do ie = 1, np_e
             ie_from = neordr(neg_g(ie),ik)
             if( ie_from < neg_previous ) cycle
               eko_l(ie,ik) = 1.d+15
               occup_l(ie,ik) = 0.d0
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
   use m_Parallelization,     only : mpi_ke_world, mpi_kg_world, mpi_ge_world, mpi_chg_world

    integer, intent(in) :: nfcntn_bin
    logical, intent(in) :: F_CNTN_BIN_partitioned
    integer, optional, intent(in) :: totch_flag

    integer  :: ik, ie, mpi_comm1,mpi_comm2
    integer, allocatable, dimension(:,:) :: n_wk, n2_mpi  ! MPI
    real(DP),allocatable, dimension(:,:) :: e_wk, e2_mpi  ! MPI
    integer    :: ispin
    integer  :: id_sname = -1
                                                  __TIMER_SUB_START(1371)
    call tstatc0_begin('m_ESIO_wd_EigenValues_etc ',id_sname)

    mpi_comm1 = mpi_kg_world
    mpi_comm2 = mpi_ge_world

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
          do ie = 1, np_e
             e_wk(ie,ik-ista_k+1) = eko_l(ie,ik)
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
          do ie = 1, np_e
             e_wk(ie,ik-ista_k+1) = occup_l(ie,ik)
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
                                                  __TIMER_IOCOMM_START_w_BARRIER(mpi_comm1,1424)
          call mpi_allreduce(MPI_IN_PLACE,n_wk,neg*kv3,mpi_integer,mpi_sum,mpi_comm1,ierr)
          call mpi_allreduce(MPI_IN_PLACE,n_wk,neg*kv3,mpi_integer,mpi_sum,mpi_comm2,ierr)
          n2_mpi = n_wk/nrank_e
                                                  __TIMER_IOCOMM_STOP(1424)
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
                                                  __TIMER_IOCOMM_START_w_BARRIER(mpi_comm1,1427)
          call mpi_allreduce(MPI_IN_PLACE,n_wk,neg*kv3,mpi_integer,mpi_sum,mpi_comm1,ierr)
          call mpi_allreduce(MPI_IN_PLACE,n_wk,neg*kv3,mpi_integer,mpi_sum,mpi_comm2,ierr)
          n2_mpi = n_wk/nrank_e
                                                  __TIMER_IOCOMM_STOP(1427)
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
          do ie = 1, np_e                                   ! MPI
             e_wk(neordr(neg_g(ie),ik),ik) = eko_l(ie,ik)           ! MPI
          end do
       end do
       end do
                                                  __TIMER_IODO_STOP(1429)
       if(npes >= 2) then
                                                  __TIMER_IOCOMM_START_w_BARRIER(mpi_comm1,1430)
          call mpi_allreduce(MPI_IN_PLACE,e_wk,neg*kv3,mpi_double_precision,mpi_sum,mpi_comm1,ierr)
          call mpi_allreduce(MPI_IN_PLACE,e_wk,neg*kv3,mpi_double_precision,mpi_sum,mpi_comm2,ierr)
          e2_mpi = e_wk
                                                  __TIMER_IOCOMM_STOP(1430)
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
          do ie = 1, np_e                                   ! MPI
             e_wk(neordr(neg_g(ie),ik),ik) = occup_l(ie,ik)         ! MPI
          end do                                            ! MPI
       end do                                               ! MPI
       end do                                               ! MPI
                                                  __TIMER_IODO_STOP(1432)
       if(npes >= 2) then
                                                  __TIMER_IOCOMM_START_w_BARRIER(mpi_comm1,1433)
          call mpi_allreduce(MPI_IN_PLACE,e_wk,neg*kv3,mpi_double_precision,mpi_sum,mpi_comm1,ierr)
          call mpi_allreduce(MPI_IN_PLACE,e_wk,neg*kv3,mpi_double_precision,mpi_sum,mpi_comm2,ierr)
          e2_mpi = e_wk
                                                  __TIMER_IOCOMM_STOP(1433)
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
   use m_Parallelization,     only : mpi_ke_world, mpi_kg_world, mpi_chg_world
   use m_Parallelization,     only : mpi_ge_world

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
               call wd_k_points(ik)
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

                                                  __TIMER_IODO_START(1461)
      a_all = 0.d0
      do is = ista_spin, iend_spin
      do ik = is, kv3-nspin+is, nspin
      !do ik = 1, kv3
         if(map_k(ik) /= myrank_k) cycle
         do ie = 1, np_e
            a_all(neg_g(ie),ik) = a_l(ie,ik)
         end do
      end do
      enddo
                                                  __TIMER_IODO_STOP(1461)
      if(npes >= 2) then
                                                  __TIMER_IOCOMM_START_w_BARRIER(mpi_comm,1433)
         call mpi_allreduce(MPI_IN_PLACE,a_all,neg*kv3,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
         call mpi_allreduce(MPI_IN_PLACE,a_all,neg*kv3,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)
                                                  __TIMER_IOCOMM_STOP(1462)
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

       e_mpi = 0.d0
       do ik = 1, kv3
          if(map_k(ik) /= myrank_k) cycle
          do ie = 1, np_e
             e_mpi(neg_g(ie),ik) = occup_l(ie,ik)
          end do
       end do
       if(npes >= 2) then
          call mpi_allreduce(MPI_IN_PLACE,e_mpi,neg*kv3,mpi_double_precision &
               &                  ,mpi_sum,mpi_kg_world,ierr) ! MPI
! === DEBUG by tkato 2012/12/19 ================================================
          call mpi_allreduce(MPI_IN_PLACE,e_mpi,neg*kv3,mpi_double_precision &
                            ,mpi_sum,mpi_ge_world,ierr)
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
    use m_Parallelization, only : mpi_ke_world, mpi_chg_world
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
            &              mpi_double_precision, mpi_sum, mpi_chg_world, ierr )
! ====================================================================== 11.0

       if (mype == 0) write(nfvlc) vlhxc_mpi2

       deallocate(vlhxc_mpi); deallocate(vlhxc_mpi2)
    else
       write(nfvlc) vlhxc_l
    end if
  end subroutine m_ESIO_wd_vlhxc

  subroutine m_ESIO_rd_WFs(nfout,nfzaj, F_ZAJ_partitioned, prev_zaj)
   use m_Parallelization,     only : mpi_ke_world, mpi_kg_world
    integer, intent(in) :: nfout, nfzaj
    logical, intent(in) :: F_ZAJ_partitioned
    logical, intent(in), optional :: prev_zaj
    logical    :: pzaj
    integer    :: ik,ib,ri, i
    integer    :: j
    integer    :: id_sname = -1
    integer    :: ierror
    integer    :: kg1t
    integer    :: ispin
    integer, allocatable, dimension(:) :: npg1kt,istag1kt,iendg1kt
                                                  __TIMER_SUB_START(1372)
    call tstatc0_begin('m_ESIO_rd_WFs ',id_sname)
    pzaj = .false.
    if(present(prev_zaj)) pzaj = prev_zaj

    kg1t = kg1
    allocate(npg1kt(kv3))
    allocate(istag1kt(kv3))
    allocate(iendg1kt(kv3))
    npg1kt = np_g1k
    istag1kt = ista_g1k
    iendg1kt = iend_g1k
    if(pzaj) then
      kg1t = kg1_prev
      npg1kt = np_g1k_prev
      istag1kt = ista_g1k_prev
      iendg1kt = iend_g1k_prev
    endif
    if(precision_WFfile==SP) then
       if(ipri >= 1) write(nfout,*) ' !D Reading zaj (single_precision)'
    else
       if(ipri >= 1) write(nfout,*) ' !D Reading zaj (double_precision)'
    end if
    rewind nfzaj
    if(F_ZAJ_partitioned) then
       if(precision_WFfile==SP) then
          allocate(wf_l(maxval(npg1kt),kimg)); wf_l = 0.d0
       else
          allocate(wfdp_l(maxval(npg1kt),kimg)); wfdp_l = 0.d0
       end if
       do ik = ista_k, iend_k, af+1        ! MPI
                                                  __TIMER_IODO_START(1435)
          do ib = 1, np_e  ! MPI
             if(ista_e+ib-1 > neg_previous) cycle
                                                  __TIMER_IODO_START(1436)
           if(precision_WFfile==SP) then
             read(nfzaj) wf_l
                                                  __TIMER_IODO_STOP(1436)
             if(kimg == 1) then
                if(pzaj)then
                  do i = 1, npg1kt(ik)
                     zaj_l_prev(i,ib,ik,1) = wf_l(i,1)
                  end do
                else
                  do i = 1, npg1kt(ik)
                     zaj_l(i,ib,ik,1) = wf_l(i,1)
                  end do
                endif
             else if(kimg==2) then
                if(pzaj)then
                  do i = 1, npg1kt(ik)
                     zaj_l_prev(i,ib,ik,1) = wf_l(i,1)
                     zaj_l_prev(i,ib,ik,2) = wf_l(i,2)
                  end do
                else
                  do i = 1, npg1kt(ik)
                     zaj_l(i,ib,ik,1) = wf_l(i,1)
                     zaj_l(i,ib,ik,2) = wf_l(i,2)
                  end do
                endif
             end if
           else if(precision_WFfile==DP) then
             read(nfzaj) wfdp_l
                                                  __TIMER_IODO_STOP(1436)
             if(kimg == 1) then
                if(pzaj) then
                  do i = 1, npg1kt(ik)
                     zaj_l_prev(i,ib,ik,1) = wfdp_l(i,1)
                  end do
                else
                  do i = 1, npg1kt(ik)
                     zaj_l(i,ib,ik,1) = wfdp_l(i,1)
                  end do
                endif
             else if(kimg==2) then
                if(pzaj)then
                  do i = 1, npg1kt(ik)
                     zaj_l_prev(i,ib,ik,1) = wfdp_l(i,1)
                     zaj_l_prev(i,ib,ik,2) = wfdp_l(i,2)
                  end do
                else
                  do i = 1, npg1kt(ik)
                     zaj_l(i,ib,ik,1) = wfdp_l(i,1)
                     zaj_l(i,ib,ik,2) = wfdp_l(i,2)
                  end do
                endif
             end if

           end if
          end do
                                                  __TIMER_IODO_STOP(1435)
       end do
       if(precision_WFfile==SP) then
          deallocate(wf_l)
       else if(precision_WFfile==DP) then
          deallocate(wfdp_l)
       end if
    else
#ifdef _DEBUG_ESIO_
       if(mype == 0) write(nfout,'("### zaj reading")')
#endif
       if(precision_WFfile==SP) then
          allocate(wf_l(kg1t,kimg));  wf_l = 0
       else
          allocate(wfdp_l(kg1t,kimg));  wfdp_l = 0
       end if
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
             call mpi_bcast(wf_l,kg1t*kimg,mpi_real4,0,MPI_CommGroup,ierr)
                                                  __TIMER_IOCOMM_STOP(1439)
             if((map_k(ik) == myrank_k) .and. (map_e(ib) == myrank_e) .and. (map_s(ispin) == myrank_spin))then         ! MPI
                if(pzaj) then
                  do ri = 1, kimg
                    do j = istag1kt(ik),iendg1kt(ik)
                       zaj_l_prev(j-istag1kt(ik)+1,map_z(ib),ik,ri) = wf_l(j,ri)  ! MPI
                    end do
                  end do
                else
                  do ri = 1, kimg
                    do j = istag1kt(ik),iendg1kt(ik)
                       zaj_l(j-istag1kt(ik)+1,map_z(ib),ik,ri) = wf_l(j,ri)  ! MPI
                    end do
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
             call mpi_bcast(wfdp_l,kg1t*kimg,mpi_double_precision,0,MPI_CommGroup,ierr)
                                                  __TIMER_IOCOMM_STOP(1439)
             if((map_k(ik) == myrank_k) .and. (map_e(ib) == myrank_e) .and. (map_s(ispin) == myrank_spin))then         ! MPI
                if(pzaj)then
                  do ri = 1, kimg
                    do j = istag1kt(ik),iendg1kt(ik)
                       zaj_l_prev(j-istag1kt(ik)+1,map_z(ib),ik,ri) = wfdp_l(j,ri)  ! MPI
                    end do
                  end do
                else
                  do ri = 1, kimg
                    do j = istag1kt(ik),iendg1kt(ik)
                       zaj_l(j-istag1kt(ik)+1,map_z(ib),ik,ri) = wfdp_l(j,ri)  ! MPI
                    end do
                  end do
                endif
             end if

           endif
! -----------------
          end do
                                                  __TIMER_IODO_STOP(1437)
       end do
       end do
       if(precision_WFfile==SP) then
          deallocate(wf_l)
       else if(precision_WFfile==DP) then
          deallocate(wfdp_l)
       end if
    end if

    deallocate(npg1kt)
    deallocate(istag1kt)
    deallocate(iendg1kt)
    call tstatc0_end(id_sname)
    return
9999 continue
    ierror = EOF_REACHED
    call phase_error_wo_filename(ierror, nfout, nfzaj, __LINE__, __FILE__)
                                                  __TIMER_SUB_STOP(1372)
  end subroutine m_ESIO_rd_WFs

! ==== TY 2019/06/25 revised the same subroutine in the 2D version m_ES_IO.f90
! ==== EXP_CELLOPT ==== 2015/09/24
  subroutine m_ESIO_import_WFs_prev_cell(nfout,nfzaj, F_ZAJ_partitioned)
    integer, intent(in) :: nfout, nfzaj
    logical, intent(in) :: F_ZAJ_partitioned
    integer    :: ik,ib,ri, i, j
    integer    :: id_sname = -1
    integer    :: ierror

    call tstatc0_begin('m_ESIO_import_WFs_prev_cell ',id_sname)

    if(precision_WFfile==SP) then
       if(ipri >= 1) write(nfout,*) ' !D Reading zaj (single_precision) in m_ESIO_import_WFs_prev_cell'
    else
       if(ipri >= 1) write(nfout,*) ' !D Reading zaj (double_precision) in m_ESIO_import_WFs_prev_cell'
    end if
    rewind nfzaj
    if(F_ZAJ_partitioned) then
       if(precision_WFfile==SP) then
          allocate(wf_l(maxval(np_g1k_prev),kimg)); wf_l = 0.d0
       else
          allocate(wfdp_l(maxval(np_g1k_prev),kimg)); wfdp_l = 0.d0
       end if
!!$
!!$    zaj_l = 0.0d0
       do ik = ista_k, iend_k, af+1        ! MPI
          do ib = 1, np_e
!!$          do ib = ista_e, iend_e, istep_e  ! MPI
             if(ista_e+ib-1 >neg_previous) cycle
!!$             if(ib > neg_previous) cycle
             if(precision_WFfile==SP) then
                read(nfzaj) wf_l

                if(kimg == 1) then
                   do i = 1, min(np_g1k(ik),np_g1k_prev(ik))
                      zaj_l(i,ib,ik,1) = wf_l(i,1)
                   end do
                else if(kimg==2) then
                   do i = 1, min(np_g1k(ik),np_g1k_prev(ik))
                      zaj_l(i,ib,ik,1) = wf_l(i,1)
                      zaj_l(i,ib,ik,2) = wf_l(i,2)
                   end do
                end if
             else if(precision_WFfile==DP) then
                read(nfzaj) wfdp_l

                if(kimg == 1) then
                   do i = 1, min(np_g1k(ik),np_g1k_prev(ik))
                      zaj_l(i,ib,ik,1) = wfdp_l(i,1)
                   end do
                else if(kimg==2) then
                   do i = 1, min(np_g1k(ik),np_g1k_prev(ik))
                      zaj_l(i,ib,ik,1) = wfdp_l(i,1)
                      zaj_l(i,ib,ik,2) = wfdp_l(i,2)
                   end do
                end if

             end if
          end do

       end do
       if(precision_WFfile==SP) then
          deallocate(wf_l)
       else if(precision_WFfile==DP) then
          deallocate(wfdp_l)
       end if
    else
       zaj_l = 0.d0
       if(precision_WFfile==SP) then
          allocate(wf_l(kg1_prev,kimg));  wf_l = 0
       else
          allocate(wfdp_l(kg1_prev,kimg));  wfdp_l = 0
       end if
       do ik = 1, kv3, af+1
          do ib = 1, neg_previous
! -----------------
           if(precision_WFfile==SP) then
             if(mype == 0) read(nfzaj, end = 9999, err = 9999) wf_l  ! MPI
             call mpi_bcast(wf_l,kg1_prev*kimg,mpi_real4,0,MPI_CommGroup,ierr)

             if((map_k(ik) == myrank_k) .and. (map_e(ib) == myrank_e) )then         ! MPI
                do ri = 1, kimg
                  do j = min(ista_g1k(ik),kg1_prev),min(iend_g1k(ik),kg1_prev)
                     zaj_l(j-ista_g1k(ik)+1,map_z(ib),ik,ri) = wf_l(j,ri)  ! MPI
                  end do
                end do
             end if
! -----------------
           else if(precision_WFfile==DP) then
             if(mype == 0) read(nfzaj, end = 9999, err = 9999) wfdp_l  ! MPI
             call mpi_bcast(wfdp_l,kg1_prev*kimg,mpi_double_precision,0,MPI_CommGroup,ierr)
             if((map_k(ik) == myrank_k) .and. (map_e(ib) == myrank_e) )then         ! MPI
                do ri = 1, kimg
                  do j = min(ista_g1k(ik),kg1_prev),min(iend_g1k(ik),kg1_prev)
                     zaj_l(j-ista_g1k(ik)+1,map_z(ib),ik,ri) = wfdp_l(j,ri)  ! MPI
                  end do
                end do
             end if
           endif
! -----------------
          end do
       end do
       if(precision_WFfile==SP) then
          deallocate(wf_l)
       else if(precision_WFfile==DP) then
          deallocate(wfdp_l)
       end if
    end if
    return
9999 continue
    ierror = EOF_REACHED
    call phase_error_wo_filename(ierror, nfout, nfzaj, __LINE__, __FILE__)
  end subroutine m_ESIO_import_WFs_prev_cell
! ===================== 2015/09/24
! ==== TY 2019/06/25


  subroutine m_ESIO_wd_WFs_standardout(nfout,ipriwf)
    integer, intent(in) :: nfout,ipriwf
    integer :: ik,ib,ri, i, ic, ipriwf0, icycle, icolumn, max_elements, istart, iend
    integer :: id_sname = -1
    real(kind=DP) :: phase2r, phase2i, phaser,phasei
    real(kind=SP), allocatable, dimension(:,:)  :: wf_mpi   ! work wave functions
    real(kind=DP), allocatable, dimension(:,:)  :: wfdp_mpi
    complex(kind=CMPLDP) :: exp2theta, exptheta
                                                  __TIMER_SUB_START(1379)
    call tstatc0_begin('m_ESIO_wd_WFs_stndout ',id_sname)

    ipriwf0 = ipriwf
    if(npes > 1) call mpi_bcast(ipriwf0,1,mpi_integer,0,MPI_CommGroup,ierr)

    if(ipriwf0 >= 2) then
       icolumn = 10
       if(precision_WFfile==SP) then
          allocate(wf_l(kg1,kimg+3)); wf_l = 0.d0
          allocate(wf_mpi(kg1,kimg+3))
       else
          allocate(wfdp_l(kg1,kimg+3)); wfdp_l = 0.d0
          allocate(wfdp_mpi(kg1,kimg+3))
       end if
       call mpi_barrier(MPI_CommGroup,ierr)
       if(mype == 0)  write(nfout,*) ' !wf Writing zaj '

       do ik = 1, kv3, af+1
          max_elements = iba(ik)
          if(mype == 0) write(nfout,'(" !wf   ik = ",i5)') ik
          do ib = 1, neg
             if(mype == 0) write(nfout,'(" !wf   ib = ",i5)') ib
           if(precision_WFfile==SP) then
             wf_mpi = 0.0d0
             if((map_k(ik) == myrank_k) .and. (map_e(ib) == myrank_e) )then         ! MPI
                                                  __TIMER_IODO_START(1467)
                do ri = 1, kimg
                  do i = ista_g1k(ik),iend_g1k(ik)
                     wf_mpi(i,ri) = zaj_l(i-ista_g1k(ik)+1,map_z(ib),ik,ri)
                  end do
                end do
                                                  __TIMER_IODO_STOP(1467)
             endif
                                                  __TIMER_IOCOMM_START_w_BARRIER(MPI_CommGroup,1468)
             call mpi_allreduce(wf_mpi,wf_l,kg1*kimg,mpi_real4,mpi_sum,MPI_CommGroup,ierr)
                                                  __TIMER_IOCOMM_STOP(1468)
                                                  __TIMER_IODO_START(1469)
           else
             wfdp_mpi = 0.0d0
             if((map_k(ik) == myrank_k) .and. (map_e(ib) == myrank_e) )then         ! MPI
                                                  __TIMER_IODO_START(1467)
                do ri = 1, kimg
                  do i = ista_g1k(ik),iend_g1k(ik)
                     wfdp_mpi(i,ri) = zaj_l(i-ista_g1k(ik)+1,map_z(ib),ik,ri)
                  end do
                end do
                                                  __TIMER_IODO_STOP(1467)
             endif
                                                  __TIMER_IOCOMM_START_w_BARRIER(MPI_CommGroup,1468)
             call mpi_allreduce(wfdp_mpi,wf_l,kg1*kimg,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
                                                  __TIMER_IOCOMM_STOP(1468)
                                                  __TIMER_IODO_START(1469)
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
                                                  __TIMER_IODO_STOP(1469)
          end do
       end do
      if(precision_WFfile==SP) then
       deallocate(wf_l)
      else
       deallocate(wfdp_l)
      end if
    end if
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(1379)
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
          call mpi_bcast(wf_l,kg1*kimg,mpi_real4,0,MPI_CommGroup,ierr)
          if((map_k(ik) == myrank_k) .and. (map_e(ib) == myrank_e) )then         ! MPI
            do ri = 1, kimg
               do j = ista_g1k(ik),iend_g1k(ik)
                  zaj_l(j-ista_g1k(ik)+1,map_z(ib),ik,ri) = wf_l(j,ri)  ! MPI
               end do
            end do
          end if

! -----------------
       else if(precision_WFfile==DP) then
          if(mype == 0) read(nfzaj, end = 9999, err = 9999) wfdp_l  ! MPI
          call mpi_bcast(wfdp_l,kg1*kimg,mpi_double_precision,0,MPI_CommGroup,ierr)
          if((map_k(ik) == myrank_k) .and. (map_e(ib) == myrank_e) )then         ! MPI
             do ri = 1, kimg
               do j = ista_g1k(ik),iend_g1k(ik)
                  zaj_l(j-ista_g1k(ik)+1,map_z(ib),ik,ri) = wfdp_l(j,ri)  ! MPI
               end do
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
    integer  :: ik, ie, iks, ri, ikg, ikt, ike, ike2, ib, j, ig
    integer, allocatable, dimension(:,:) :: n_mpi, n2_mpi  ! MPI
    real(DP),allocatable, dimension(:,:) :: e_mpi, e2_mpi  ! MPI
    real(SP),allocatable, dimension(:,:) :: wf_buf
!!$    read(nf) neordr,nrvf_ordr,eko_l,occup_l,efermi,totch

    allocate(n_mpi(neg,nspin)); allocate(n2_mpi(neg,nspin)) ! MPI
    allocate(e_mpi(neg,nspin)); allocate(e2_mpi(neg,nspin)) ! MPI
    allocate(wf_l(kg1,kimg))
    allocate(wf_buf(kg1,kimg))

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

!    deallocate(wf_l)
!    allocate(wf_l(maxval(np_g1k),kimg)); wf_l = 0.d0
!!$    do ik = 1, kv3, af+1
    KPOINT_LOOP: do ikg = 1, nrank_k
       do ikt = 1, nspin, af+1
          ik = (ikg-1)*nspin+ikt
          if(nk_in_the_process -1 + ik > kv3_ek) exit KPOINT_LOOP
          if(nk_in_the_process -1 + ik > numk_zajsaved) exit KPOINT_LOOP
          if(ipri>=1) write(nfout,*) ' !D     reading  ik = ', ik+first_kpoint_in_this_job-1
          do ib = 1, neg
             if(mype == 0) read(nf,err=2,end=2) wf_l     ! MPI
#ifdef _DEBUG_ESIO_
             if(mype == 0) then
                write(nfout,'(" ik = ",i3, " ib = ",i4)')  ik, ib
                write(nfout,'(8f8.4)') (wf_l(ri,1),ri=1,8)
             end if
#endif
             call mpi_bcast(wf_l,kg1*kimg,mpi_real4,0,MPI_CommGroup,ierr)
                                                  __TIMER_IOCOMM_STOP(1439)
             if((map_k(ik) == myrank_k) .and. (map_e(ib) == myrank_e) )then         ! MPI
                do ri = 1, kimg
                  do j = ista_g1k(ik),iend_g1k(ik)
                     zaj_l(j-ista_g1k(ik)+1,map_z(ib),ik,ri) = wf_l(j,ri)  ! MPI
                  end do
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
             !if(map_ek(ie,ik) == mype) then                ! MPI
             if((map_k(ik) == myrank_k) .and. (map_e(ie) == myrank_e) )then         ! MPI
                eko_l(map_z(ie),ik) = e_mpi(ie,ikt)        ! MPI
             end if                                        ! MPI
          end do                                           ! MPI
          do ie = 1, neg
             eko_ek(ie,nk_in_the_process+ik-1) = e_mpi(n_mpi(ie,ikt),ikt)
          end do
       end do                                              ! MPI
       goto 3
2      continue
       call phase_error_with_msg(nfout,' eof from nf <<m_ESIO_rd_WFs_and_EVs_ek>>',__LINE__,__FILE__)
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
             !wf_l = 0.0
             !do ie = 1, neg
             !   if((map_k(ik) == myrank_k) .and. (map_e(ie) == myrank_e) )then         ! MPI
             !      do ri = 1, kimg
             !         do ig = ista_g1k(ik), iend_g1k(ik)
!            !           do ig = 1, np_g1k(ik)
             !            wf_l(ig,ri) = zaj_l(ig-ista_g1k(ik)+1,map_z(ie),ik,ri)
             !         enddo
             !      end do
             !   endif
             !   call mpi_allreduce(mpi_in_place, wf_l, kg1*kimg, mpi_real4, mpi_sum, mpi_ke_world,ierr)
             !   if((map_k(ike2) == myrank_k) .and. (map_e(ie) == myrank_e) )then         ! MPI
             !      do ri = 1, kimg
             !         do ig = ista_g1k(ike2), iend_g1k(ike2)
             !            zaj_l(ig-ista_g1k(ike2)+1,map_z(ie),ike2,ri) = wf_buf(ig,ri)
             !         enddo
             !      end do
             !   end if
             !end do
             if(map_k(ike2)==myrank_k) call m_ESIW_by_randomnumbers0_3D(nfout,kv3,1,neg,ike2,ike2)
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
                if((map_k(ik) == myrank_k) .and. (map_e(ie) == myrank_e) )then         ! MPI
                   e_mpi(ie,ikt) = eko_l(map_z(ie),ik)
                end if
             end do
          end do
          if(npes >= 2) then
             call mpi_allreduce(MPI_IN_PLACE,e_mpi,neg*nspin,mpi_double_precision &
                  & ,mpi_sum,mpi_kg_world,ierr)
             call mpi_allreduce(MPI_IN_PLACE,e_mpi,neg*nspin,mpi_double_precision &
                  & ,mpi_sum,mpi_ge_world,ierr)

             call mpi_allreduce(e_mpi,e2_mpi,neg*nspin,mpi_double_precision &
                  & ,mpi_sum,mpi_kg_world,ierr)
             call mpi_allreduce(e_mpi,e2_mpi,neg*nspin,mpi_double_precision &
                  & ,mpi_sum,mpi_ge_world,ierr)
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
                if((map_k(ike2) == myrank_k) .and. (map_e(ie) == myrank_e) )then         ! MPI
                   eko_l(map_z(ie),ike2) = e2_mpi(ie,ikt)
                end if
             end do
          end do

       end do
       do ik = 1, kv3, af+1
          if(map_k(ik) /= myrank_k) cycle
          call m_ES_betar_dot_WFs_3D(nfout,ik)     ! (fsrfsi,sumset)
       end do
       call m_ES_modified_gram_schmidt(nfout) ! (grmsmd)
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
    return
9999 continue
    call phase_error_wo_filename(EOF_REACHED, nfout, nf, __LINE__, __FILE__)

!!$    write(nfout,'(" ---<< m_ESIO_rd_WFs_and_EVs_ek>>---")')
  end subroutine m_ESIO_rd_WFs_and_EVs_ek


  subroutine m_ESIO_wd_WFs_and_EVs_ek(nfout,nf)
    integer, intent(in) :: nfout, nf
    integer  :: ik, ie, ri, ikg, ikt, ib, j
    integer, allocatable, dimension(:,:) :: n_mpi
    real(DP),allocatable, dimension(:,:) :: e_mpi
    integer :: id_sname = -1
    call tstatc0_begin('m_ESIO_wd_WFs_and_EVs_ek ',id_sname)

    allocate(wf_l(kg1,kimg));    wf_l = 0

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
#ifdef _USE_ALLREDUCE_IN_WD_WFS_
          do ib = 1, neg
             wf_l = 0.0d0
             if((map_k(ik) == myrank_k) .and. (map_e(ib) == myrank_e)) then         ! MPI
                do ri = 1, kimg
                  do j = ista_g1k(ik),iend_g1k(ik)
                     wf_l(j,ri) = zaj_l(j-ista_g1k(ik)+1,map_z(ib),ik,ri)
                  end do
                end do
             endif
             call mpi_allreduce(MPI_IN_PLACE,wf_l,kg1*kimg,mpi_real4,mpi_sum,MPI_CommGroup,ierr)
             if(mype == 0) write(nf)  wf_l                        ! MPI
          end do
#else
!! coded by T. Yamasaki, 2021.02.25 -->
          allocate(wf_buf(maxval(np_g1k),kimg))
          do ib = 1, neg
             wf_l = 0.0d0
             wf_buf = 0.d0
             do j = 0, nrank_g-1
                if(map_k(ik)==myrank_k .and. map_e(ib)==myrank_e .and. j==myrank_g) then
                   do ig = ista_g1k(ik), iend_g1k(ik)
                      do ri=1, kimg
                         wf_buf(ig-ista_g1k(ik)+1,ri) = zaj_l(ig-ista_g1k(ik)+1,map_z(ib),ik,ri)
                      end do
                   end do
                   call mpi_send(wf_buf,maxval(np_g1k)*kimg,mpi_real4,0,1,MPI_CommGroup,ierr)
                end if
                if(mype==0) then
                   call mpi_recv(wf_buf,maxval(np_g1k)*kimg,mpi_real4,map_rank_gek(j,map_e(ib),map_k(ik)),1,MPI_CommGroup,istatus,ierr)
                   wf_l(nis_g1k(j,ik):nie_g1k(j,ik),1:kimg) = wf_buf(1:nie_g1k(j,ik)-nis_g1k(j,ik)+1,1:kimg)
                end if
             end do
             if(mype==0) write(nf) wf_l
          end do
          deallocate(wf_buf)
!! <--
#endif
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
          n_mpi = n_mpi/nrank_e/nrank_g
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
          n_mpi = n_mpi/nrank_e/nrank_g
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
               &               ,mpi_sum,mpi_kg_world,ierr) ! MPI
          call mpi_allreduce(MPI_IN_PLACE,e_mpi,neg*nspin,mpi_double_precision &
               &               ,mpi_sum,mpi_ge_world,ierr) ! MPI
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

  subroutine m_ESIO_wd_WFs_3D(nfout,nfzaj,F_ZAJ_partitioned,rew)
   use m_Parallelization,     only : mpi_ke_world, mpi_kg_world

    integer, intent(in) :: nfout,nfzaj
    logical, intent(in) :: F_ZAJ_partitioned
    logical, intent(in), optional :: rew
    real(kind=SP), allocatable, dimension(:,:) :: wf_buf
    integer :: ik,ib,ri,j,i,ig, ispin
    integer :: id_sname = -1
                                                  __TIMER_SUB_START(1373)
    call tstatc0_begin('m_ESIO_wd_WFs ',id_sname)

    call mpi_barrier(MPI_CommGroup,ierr)
    if (present(rew)) then
      if(rew) then
        rewind nfzaj
      endif
    else
      rewind nfzaj
    endif
   if(precision_WFfile==SP) then
    if(ipri >= 1) write(nfout,*) ' !D Writing zaj (single_precision)'
    if(F_ZAJ_partitioned) then
       allocate(wf_l(maxval(np_g1k),kimg))
                                                  __TIMER_IODO_START(1440)
       do ik = ista_k, iend_k, af+1
          do ib = 1, np_e
             do ri = 1, kimg
                wf_l(1:np_g1k(ik),ri) = zaj_l(1:np_g1k(ik),ib,ik,ri)
             end do
                                                  __TIMER_IODO_START(1441)
             write(nfzaj) wf_l
                                                  __TIMER_IODO_STOP(1441)
          end do
       end do
                                                  __TIMER_IODO_STOP(1440)
       deallocate(wf_l)
    else
       allocate(wf_l(kg1,kimg))
       do ispin = 1, nspin, af+1
       !do ik = 1, kv3, af+1
       do ik = ispin, kv3-nspin+ispin, nspin
                                                  __TIMER_IODO_START(1442)
#ifdef _USE_ALLREDUCE_IN_WD_WFS_
          do ib = 1, neg
             wf_l = 0.0d0
             if((map_k(ik) == myrank_k) .and. (map_e(ib) == myrank_e) .and. (map_s(ispin) == myrank_spin) ) then         ! MPI
                do ri = 1, kimg
                  do j = ista_g1k(ik),iend_g1k(ik)
                     wf_l(j,ri) = zaj_l(j-ista_g1k(ik)+1,map_z(ib),ik,ri)
                  end do
                end do
             endif
                                                  __TIMER_IOCOMM_START_w_BARRIER(MPI_CommGroup,1443)
             call mpi_allreduce(MPI_IN_PLACE,wf_l,kg1*kimg,mpi_real4,mpi_sum,MPI_CommGroup,ierr)
                                                  __TIMER_IOCOMM_STOP(1443)
                                                  __TIMER_IODO_START(1444)
             if(mype == 0) write(nfzaj)  wf_l                        ! MPI
                                                  __TIMER_IODO_STOP(1444)
          end do
#else
!! coded by T. Yamasaki, 2021.02.25 -->
          allocate(wf_buf(maxval(np_g1k),kimg))
          do ib = 1, neg
             wf_l = 0.0d0
             wf_buf = 0.d0
             do j = 0, nrank_g-1
                if(map_k(ik)==myrank_k .and. map_e(ib)==myrank_e .and. j==myrank_g) then
                   do ig = ista_g1k(ik), iend_g1k(ik)
                      do ri=1, kimg
                         wf_buf(ig-ista_g1k(ik)+1,ri) = zaj_l(ig-ista_g1k(ik)+1,map_z(ib),ik,ri)
                      end do
                   end do
                   call mpi_send(wf_buf,maxval(np_g1k)*kimg,mpi_real4,0,1,MPI_CommGroup,ierr)
                end if
                if(mype==0) then
                   call mpi_recv(wf_buf,maxval(np_g1k)*kimg,mpi_real4,map_rank_gek(j,map_e(ib),map_k(ik)),1,MPI_CommGroup,istatus,ierr)
                   wf_l(nis_g1k(j,ik):nie_g1k(j,ik),1:kimg) = wf_buf(1:nie_g1k(j,ik)-nis_g1k(j,ik)+1,1:kimg)
                end if
             end do
             if(mype==0) write(nfzaj) wf_l
          end do
          deallocate(wf_buf)
!! <--
#endif
                                                  __TIMER_IODO_STOP(1442)
       end do
       enddo
       deallocate(wf_l)
    end if
   else if(precision_WFfile==DP) then
    if(ipri >= 1) write(nfout,*) ' !D Writing zaj (double_precision)'
    if(F_ZAJ_partitioned) then
       allocate(wfdp_l(maxval(np_g1k),kimg))
                                                  __TIMER_IODO_START(1440)
       do ik = ista_k, iend_k, af+1
          do ib = 1, np_e
             do ri = 1, kimg
                wfdp_l(1:np_g1k(ik),ri) = zaj_l(1:np_g1k(ik),ib,ik,ri)
             end do
                                                  __TIMER_IODO_START(1441)
             write(nfzaj) wfdp_l
                                                  __TIMER_IODO_STOP(1441)
          end do
       end do
                                                  __TIMER_IODO_STOP(1440)
       deallocate(wfdp_l)
    else
       allocate(wfdp_l(kg1,kimg))
       do ispin = 1, nspin, af+1
       !do ik = 1, kv3, af+1
       do ik = ispin, kv3-nspin+ispin, nspin
                                                  __TIMER_IODO_START(1442)
          do ib = 1, neg
             wfdp_l = 0.0d0
             if((map_k(ik) == myrank_k) .and. (map_e(ib) == myrank_e) .and. (map_s(ispin) == myrank_spin))then         ! MPI
                do ri = 1, kimg
                  do j = ista_g1k(ik),iend_g1k(ik)
                     wfdp_l(j,ri) = zaj_l(j-ista_g1k(ik)+1,map_z(ib),ik,ri)
                  end do
                end do
             endif
                                                  __TIMER_IOCOMM_START_w_BARRIER(MPI_CommGroup,1443)
             call mpi_allreduce(MPI_IN_PLACE,wfdp_l,kg1*kimg,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
                                                  __TIMER_IOCOMM_STOP(1443)
                                                  __TIMER_IODO_START(1444)
             if(mype == 0) write(nfzaj)  wfdp_l
                                                  __TIMER_IODO_STOP(1444)
          end do
                                                  __TIMER_IODO_STOP(1442)
       end do
       end do
       deallocate(wfdp_l)
    end if
   end if
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(1373)
  end subroutine m_ESIO_wd_WFs_3D

!-- for test
  subroutine m_ESIO_rd_WFs_dp_3D(nfout,nfzaj, F_ZAJ_partitioned)
   use m_Parallelization,     only : mpi_ke_world, mpi_kg_world

    real(kind=DP), allocatable, dimension(:,:)  :: wf_ldp   ! work wave functions
    integer, intent(in) :: nfout, nfzaj
    logical, intent(in) :: F_ZAJ_partitioned
    integer    :: ik,ib,ri, i, j
    integer    :: id_sname = -1
    call tstatc0_begin('m_ESIO_rd_WFs ',id_sname)

!f    allocate(wf_ldp(kg1,kimg))
!!$ASASASASAS
!f    wf_l = 0
!!$ASASASASAS
    rewind nfzaj
    if(ipri >= 1) write(nfout,*) ' !D Reading zaj'
    if(F_ZAJ_partitioned) then
       allocate(wf_ldp(maxval(np_g1k),kimg))
       wf_ldp = 0
       do ik = ista_k, iend_k, af+1        ! MPI
          do ib = 1, np_e  ! MPI
             if(ista_e+ib-1 > neg_previous) cycle
             read(nfzaj) wf_ldp
             if(kimg == 1) then
                do i = 1, np_g1k(ik)
                   zaj_l(i,ib,ik,1) = wf_ldp(i,1)
                end do
             else if(kimg==2) then
                do i = 1, np_g1k(ik)
                   zaj_l(i,ib,ik,1) = wf_ldp(i,1)
                   zaj_l(i,ib,ik,2) = wf_ldp(i,2)
                end do
             end if
          end do
       end do
       deallocate(wf_ldp)
    else
       allocate(wf_ldp(kg1,kimg))
       wf_ldp = 0
       do ik = 1, kv3, af+1
          do ib = 1, neg_previous
             if(mype == 0) read(nfzaj) wf_ldp              ! MPI
             call mpi_bcast(wf_ldp,kg1*kimg,mpi_real8,0,MPI_CommGroup,ierr)
             if((map_k(ik) == myrank_k) .and. (map_e(ib) == myrank_e) )then         ! MPI
                do ri = 1, kimg
                  do j = ista_g1k(ik),iend_g1k(ik)
                     zaj_l(j-ista_g1k(ik)+1,map_z(ib),ik,ri) = wf_ldp(j,ri)  ! MPI
                  end do
                end do
             end if
          end do
       end do
       deallocate(wf_ldp)
    end if

    call tstatc0_end(id_sname)
  end subroutine m_ESIO_rd_WFs_dp_3D

  subroutine m_ESIO_wd_WFs_dp_3D(nfout,nfzaj,F_ZAJ_partitioned)
   use m_Parallelization,     only : mpi_ke_world, mpi_kg_world

    real(kind=DP), allocatable, dimension(:,:)  :: wf_ldp   ! work wave functions
    integer, intent(in) :: nfout,nfzaj
    logical, intent(in) :: F_ZAJ_partitioned
    integer :: ik,ib,ri,j,i
    integer :: id_sname = -1
    real(kind=DP), allocatable, dimension(:,:)  :: wf_mpi   ! work wave functions
    call tstatc0_begin('m_ESIO_wd_WFs ',id_sname)

!f    allocate(wf_ldp(kg1,kimg))
    call mpi_barrier(MPI_CommGroup,ierr)
    !!$ print *, ' !D Writing zaj '
    if(ipri >= 1) write(nfout,*) ' !D Writing zaj '
    rewind nfzaj
    if(F_ZAJ_partitioned) then
       allocate(wf_ldp(maxval(np_g1k),kimg))
       allocate(wf_mpi(kg1,kimg))
       do ik = ista_k, iend_k, af+1        ! MPI
          do ib = 1, np_e
             do ri = 1, kimg
                wf_ldp(1:np_g1k(ik),ri) = zaj_l(1:np_g1k(ik),ib,ik,ri)
             end do
             write(nfzaj) wf_ldp
          end do
       end do
       deallocate(wf_ldp)
       deallocate(wf_mpi)
    else
       allocate(wf_ldp(kg1,kimg))
       allocate(wf_mpi(kg1,kimg))
       do ik = 1, kv3, af+1
          do ib = 1, neg
             wf_ldp = 0.0d0
             if((map_k(ik) == myrank_k) .and. (map_e(ib) == myrank_e) )then         ! MPI
                do ri = 1, kimg
                  do j = ista_g1k(ik),iend_g1k(ik)
                     wf_ldp(j,ri) = zaj_l(j-ista_g1k(ik)+1,map_z(ib),ik,ri)
                  end do
                end do
             endif
             call mpi_allreduce(wf_ldp,wf_mpi,kg1*kimg,mpi_real8,mpi_sum,MPI_CommGroup,ierr)
             if(mype == 0) write(nfzaj)  wf_mpi                        ! MPI
          end do
       end do
       deallocate(wf_ldp)
       deallocate(wf_mpi)
    end if
    call tstatc0_end(id_sname)
  end subroutine m_ESIO_wd_WFs_dp_3D


end module m_ES_IO
