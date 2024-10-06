!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 626 $)
!
!  SUBROUINE: WriteDownData_onto_Files_ek, FermiEnergyLevel_ek()
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  FURTHER MODIFICATION: T. Yamasaki, December/01/2003
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
!$$#ifndef PARA3D
subroutine WriteDownData_onto_Files_ek
! $Id: WriteDownData_onto_Files_ek.f90 626 2020-06-03 02:18:38Z jkoga $
  use m_Files,               only : m_Files_open_nfcntn,m_Files_close_all &
       &                     ,nfcntn,nfchgt,nfzaj,m_Files_open_nfcntn_bin &
       &                     ,nfcntn_bin,nfout,nfeng, nfstatus &
       &                     ,F_CNTN_BIN_partitioned &
       &                     , m_Files_open_nfstatus, m_Files_skiptoend &
       &                     , m_Files_open_nfzaj_with_check, m_Files_close_logfile &
       &                     , F_ZAJ
  use m_Const_Parameters,    only: EK,SCF,YES,EK_CONVERGED,PARABOLIC,TETRAHEDRON,ON,FINISH, OFF
  use m_Timing,              only: m_Timing_wd_status
  use m_Control_parameters,  only: ipri,ipriekzaj, neg, num_extra_bands, nspin &
       &                         , iconvergence_previous_job, iprijobstatus &
       &                         , jobstatus_series, jobstatus_format &
       &                         , continuation_using_ppdata, ekmode, neg_is_enlarged&
       &                         , m_CtrlP_wd_cpu_total, m_CtrlP_wd_iconvergence &
       &                         , m_CtrlP_wd_iconv_ek &
       &                         , m_CtrlP_way_of_smearing &
       &                         , m_CtrlP_wd_numk_zajsaved &
       &                         , sw_berry_phase, sw_write_bxsf_file, sw_write_zaj
  use m_IterationNumbers,    only: iteration, iteration_ionic, iteration_electronic &
       &                         , nk_in_the_process &
       &                         , m_Iter_wd_iters_and_nk, m_Iter_electronic_reset
  use m_Parallelization,     only: mype
  use m_Kpoints,             only: kv3, kv3_ek
  use m_PseudoPotential,     only: m_PP_wd_PP_parameters_3D, m_PP_wd_betar
  use m_PlaneWaveBasisSet,   only: kgp
  use m_NonLocal_Potential,  only: m_NLP_wd_snl_3D
  use m_ES_LHXC,             only: m_ESlhxc_potential_3D
  use m_ES_IO,               only: m_ESIO_wd_EigenValues_ek &
       &                         , m_ESIO_wd_WFs_and_EVS_ek
  use m_ES_occup,            only: m_ESoc_fermi_parabolic_ek, m_ESoc_fermi_tetra_ek &
       &                         , m_ESoc_check_num_bands
  use m_Electronic_Structure,only: efermi, vbm, metalic_system, iconv_ek &
       &                         , m_ES_wd_zaj_small_portion0
  use m_ES_nonlocal,         only: m_ES_wd_fsr_fsi

! ============================== added by K. Tagami =============== 11.0
  use m_Control_Parameters,  only : noncol, ndim_spinor
! ================================================================== 11.0

! ============================== KT_Add ===================== 13.0E
  use m_Const_Parameters,     only : Fermi_Dirac
  use m_ES_occup,            only: m_ESoc_fermi_dirac_ek
! =========================================================== 13.0E

  use m_BerryPhase,          only : m_BP_wd_cntn_data

  use m_Control_Parameters,  only : sw_corelevel_spectrum,  sw_local_approx_trans_moment
  use m_CLS_dipquad,        only : m_CLS_wd_transmom_ek

  use m_Crystal_Structure,  only : sw_spinorbit_second_variation
  use m_SpinOrbit_SecondVariation,  only : m_SO_wd_EigenValues_ek_fin, &
       &                                   m_SO_calc_band_energy_socsv_ek

  implicit none
  integer :: status_wdmode, numk_zajsaved
  logical :: Allkpoints_are_Calculated, enough_bands

  logical :: flg_wd_zaj = .true.

  flg_wd_zaj = .true.
  if ( F_ZAJ == "/dev/null" .or. sw_write_zaj == OFF ) flg_wd_zaj = .false.

  if(iconvergence_previous_job < EK_CONVERGED) then
     call m_Files_open_nfcntn                        ! -(m_Files)
     if(mype==0) rewind nfcntn
     if(ipriekzaj<=0) call m_Iter_electronic_reset()
     call m_Iter_wd_iters_and_nk(nfcntn)
     call m_CtrlP_wd_iconvergence(nfcntn)
     numk_zajsaved = min(kv3_ek,nk_in_the_process+kv3-1)
     call m_CtrlP_wd_numk_zajsaved(nfcntn,numk_zajsaved)
     call m_CtrlP_wd_iconv_ek(kv3_ek,iconv_ek,nfcntn)

     call m_ESIO_wd_EigenValues_ek(nfout,mode=SCF)

     if(continuation_using_ppdata == YES) then
        call m_Files_open_nfcntn_bin
        if(mype==0) rewind nfcntn_bin
        call m_PP_wd_PP_parameters_3D(nfout,nfcntn_bin,F_CNTN_BIN_partitioned,kgp)
        call m_PP_wd_betar(nfcntn_bin)
        call m_NLP_wd_snl_3D(nfout,nfcntn_bin,F_CNTN_BIN_partitioned,kv3)
     end if

     if(Allkpoints_are_Calculated()) then
        enough_bands = m_ESoc_check_num_bands()
        if(enough_bands) then
           call FermiEnergylevel_ek()   ! -> metalic_system?, vbm, efermi
!!$           stop ' WriteDownData_onto_Files_ek (1)'
           if(ipri >= 1) write(nfout,'(" --- efermi = ",f8.4)') efermi
        else
           if(ipri >= 1) write(nfout,'(" --- efermi can not be estimated --")')
        end if
        rewind nfeng

        if(ipri >= 1) then
! -------------------------------------
!
! --- The following modification is intended to work band.pl correctly.
! 
! =============== modified by K. Tagami ======================= 11.0
!           write(nfeng,'(" num_kpoints = ",i6)') kv3_ek
!           write(nfeng,'(" num_bands   = ",i6)') neg-num_extra_bands
!           write(nfeng,'(" nspin       = ",i6)') nspin

!!$           if ( noncol ) then
           write(nfeng,'(" num_kpoints = ",i6)') kv3_ek /ndim_spinor
           if(neg_is_enlarged) then
              write(nfeng,'(" num_bands   = ",i6)') neg-num_extra_bands
           else
              write(nfeng,'(" num_bands   = ",i6)') neg
           end if
           write(nfeng,'(" nspin       = ",i6)') nspin / ndim_spinor
!!$           else
!!$              write(nfeng,'(" num_kpoints = ",i6)') kv3_ek
!!$              write(nfeng,'(" num_bands   = ",i6)') neg-num_extra_bands
!!$              write(nfeng,'(" nspin       = ",i6)') nspin
!!$           endif
! =============================================================== 11.0
        end if

        if(enough_bands) then
           if(metalic_system) then
              if(ipri>=1) write(nfeng,'(" Fermi energy level = ",f10.6/)') efermi
           else
!!$           write(nfeng,'(" The Highest occupied band energy = ",f10.6/)') vbm
              if(ipri>=1) write(nfeng,'(" Valence band max   = ",f10.6/)') vbm
           end if
        else
           if(ipri>=1) write(nfeng,'(" Fermi energy level = unknown")')
        end if
!!$        call m_ES_wd_eko(nfeng,mode=EK)
        call m_ESIO_wd_EigenValues_ek(nfeng,mode=EK)

        if ( .not. noncol .and. sw_spinorbit_second_variation == ON ) then
           call m_SO_wd_EigenValues_ek_fin
           call m_SO_calc_band_energy_socsv_ek
        endif

     else
!!$     if(.not.Allkpoints_are_Calculated() .and. ipriekzaj >= 1) &
        if(ipriekzaj >= 1) then
           if(ipri >= 1) write(nfout,'(" ipriekzaj = ",i8," <<WriteDownData_onto_Files_ek>>")') ipriekzaj
           call m_Files_open_nfzaj_with_check()
           if ( flg_wd_zaj ) then
              call m_ESIO_wd_WFs_and_EVs_ek(nfout,nfzaj)
           endif
        end if
     end if
  else
     if(ipri>=1) write(nfout,'(" --- iconvergence_previous >= EK_CONVERGED ---")')
  end if
  if(sw_berry_phase==ON)then
     call m_BP_wd_cntn_data()
  endif

  if ( sw_corelevel_spectrum == ON .and. sw_local_approx_trans_moment == ON ) then
     call m_CLS_wd_transmom_ek
  endif

  call m_ES_wd_fsr_fsi()                    ! if(ipri >= 2) 
  call m_ES_wd_zaj_small_portion0(" -- WDonFiles_ek --",18)
!!$  stop ' WriteDownData_onto_Files_ek (2)'

  call m_CtrlP_wd_cpu_total()
  if(iprijobstatus >= 1) then
     call m_Files_open_nfstatus()
     if(jobstatus_series == ON) then
        call m_Files_skiptoend(nfstatus)
     else
     end if
     status_wdmode = FINISH
     call m_Timing_wd_status(nfstatus,jobstatus_format,jobstatus_series,status_wdmode &
          &                , iteration,iteration_ionic,iteration_electronic)
  end if 

  call m_Files_close_all()                    ! -(m_Files)
  call PrintStatus()
  call m_Files_close_logfile()
!!$  stop ' WriteDownData_onto_Files_ek (3)'
contains

  subroutine FermiEnergyLevel_ek()
    integer :: way_of_smearing
    way_of_smearing = m_CtrlP_way_of_smearing()
    if(way_of_smearing == PARABOLIC) then
       call m_ESoc_fermi_parabolic_ek(nfout)  ! -> efermi, metalic_system
!!$  else if(way_of_smearing == MP) then
!!$     call fermi_mesfessel_paxton(nfout)
     else if(way_of_smearing == TETRAHEDRON) then
        call m_ESoc_fermi_tetra_ek(nfout)     ! -> efermi, metalic_system

! ========================= KT_add ================= 13.0E
     else if(way_of_smearing == Fermi_Dirac) then
        call m_ESoc_fermi_dirac_ek(nfout)     ! -> efermi, metalic_system
! ================================================== 13.0E
    end if
  end subroutine FermiEnergyLevel_ek

end subroutine WriteDownData_onto_Files_ek
!$$#endif

subroutine WriteDown_totalcpu()
  use m_Control_parameters,  only: m_CtrlP_wd_cpu_total

  call m_CtrlP_wd_cpu_total()
end subroutine WriteDown_totalcpu

subroutine write_energies(mnk,nk,kxyz,ek,mnt,nt,tk,ie)
  use m_Const_Parameters, only : DP
  use m_Files, only :            nfegrid, m_Files_open_nfegrid
  implicit none

  integer, intent(in) :: mnk,nk,mnt,nt,ie
  real(kind=DP), intent(in), dimension(3,mnk) :: kxyz
  real(kind=DP), intent(in), dimension(mnk) ::   ek
  integer, intent(in), dimension(4,mnt) ::       tk

  integer :: it, ik

  call m_Files_open_nfegrid()  ! open nfegrid

  write (nfegrid,*) '$ DATA=CURVE3D NAME = ', nt

  do it = 1, nt
     if( kxyz(3,tk(1,it)) == 0.d0 .and. kxyz(3,tk(2,it)) == 0.d0 .and. &
          &       kxyz(3,tk(3,it)) == 0d0) then
        write (nfegrid,*) sngl(kxyz(1,tk(1,it))),sngl(kxyz(2,tk(1,it))),sngl(ek(tk(1,it))),tk(1,it)
        write (nfegrid,*) sngl(kxyz(1,tk(2,it))),sngl(kxyz(2,tk(2,it))),sngl(ek(tk(2,it))),tk(2,it)
        write (nfegrid,*) sngl(kxyz(1,tk(3,it))),sngl(kxyz(2,tk(3,it))),sngl(ek(tk(3,it))),tk(3,it)
     endif

     if (kxyz(3,tk(1,it)) == 0d0 .and. kxyz(3,tk(2,it)) == 0d0 .and. &
          &       kxyz(3,tk(4,it)) == 0d0) then
        write (nfegrid,*)
        write (nfegrid,*) sngl(kxyz(1,tk(1,it))),sngl(kxyz(2,tk(1,it))),sngl(ek(tk(1,it))),tk(1,it)
        write (nfegrid,*) sngl(kxyz(1,tk(2,it))),sngl(kxyz(2,tk(2,it))),sngl(ek(tk(2,it))),tk(2,it)
        write (nfegrid,*) sngl(kxyz(1,tk(4,it))),sngl(kxyz(2,tk(4,it))),sngl(ek(tk(4,it))),tk(4,it)
     endif

     if (kxyz(3,tk(1,it)) == 0d0 .and. kxyz(3,tk(3,it)) == 0d0 .and. &
          &      kxyz(3,tk(4,it)) == 0d0) then
        write (nfegrid,*)
        write (nfegrid,*) sngl(kxyz(1,tk(1,it))),sngl(kxyz(2,tk(1,it))),sngl(ek(tk(1,it))),tk(1,it)
        write (nfegrid,*) sngl(kxyz(1,tk(3,it))),sngl(kxyz(2,tk(3,it))),sngl(ek(tk(3,it))),tk(3,it)
        write (nfegrid,*) sngl(kxyz(1,tk(4,it))),sngl(kxyz(2,tk(4,it))),sngl(ek(tk(4,it))),tk(4,it)
     endif

     if (kxyz(3,tk(2,it)) == 0d0 .and. kxyz(3,tk(3,it)) == 0d0 .and. &
          &     kxyz(3,tk(4,it)) == 0d0) then
        write (nfegrid,*)
        write (nfegrid,*) sngl(kxyz(1,tk(2,it))),sngl(kxyz(2,tk(2,it))),sngl(ek(tk(2,it))),tk(2,it)
        write (nfegrid,*) sngl(kxyz(1,tk(3,it))),sngl(kxyz(2,tk(3,it))),sngl(ek(tk(3,it))),tk(3,it)
        write (nfegrid,*) sngl(kxyz(1,tk(4,it))),sngl(kxyz(2,tk(4,it))),sngl(ek(tk(4,it))),tk(4,it)
     endif
     write (nfegrid,*)
  enddo

  write (nfegrid,*) '$ end  ', ie, nk

!_____write output
  write (nfegrid,*) nk
  do ik = 1, nk
     write (nfegrid,*) kxyz(1,ik), kxyz(2,ik), kxyz(3,ik), ek(ik)
  enddo

  write (nfegrid,*) nt
  do it = 1, nt
     write (nfegrid,*) tk(1,it), tk(2,it), tk(3,it), tk(4,it)
  enddo
end subroutine write_energies
