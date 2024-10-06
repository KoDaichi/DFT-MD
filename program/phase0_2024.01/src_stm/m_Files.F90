!================================================
!  Software name : STM
!  Module : m_Files
!  Subroutine(s) : open_nfcntn_bin, m_Files_set_default_filenames,
!                  m_Files_read_file_names_data, m_Files_open_files_initially,
!                  open0, m_Files_close_all_files
!  Author(s)     : Takahiro Yamasaki and Koichi Kato (June 7, 2004)
!
!  Contact address :  IIS,The University of Tokyo RSS21 project
!  
!  "Multiscale Simulation System for Function Analysis of Nanomaterials"  
!
!================================================
!
!     The original version of this program "STM" was developed in the
!  period 1998-2001 at JRCAT by Koichi Kato (Toshiba Corporate
!  Research and Development Center) and Takahiro Yamasaki (Fujitsu 
!  Laboratories Ltd.) through a joint research between JRCAT and
!  Toshiba Corporate Research and Development Center and another joint
!  research between JRCAT and FUJITSU Laboratories Ltd.
!     Since 2003, this program has been tuned and revised as it
!  matches to the first-principles molecular dynamics program "PHASE"
!  as a part of the national project "Frontier Simulation Software for
!  Industrial Science (FSIS)", which is supported by the IT program of
!  the Ministry of Education, Culture, Sports, Science and Technology
!  (MEXT) of Japan.
!     Since 2006, this program set has been developed as a part of the
!  national project "Revolutionary Simulation Software (RSS21)", which
!  is supported by the next-generation IT program of MEXT of Japan.

module m_Files
! $Id: m_Files.F90,v 1.1.1.1 2004/06/26 11:21:37 yamasaki Exp $
  use m_Const_Parameters, only : OLD, NEW, UNKNOWN, FORMATTED &
       & ,UNFORMATTED, check_file_name_on, check_file_name_off, DP &
       & ,ON, GENERAL

  implicit none

  integer :: nfinp, nfout, nfcntn_bin, nfchgu, nfchgd, nfzaj, nfvlc, nfchgu_p, nfchgd_p, nfvlcr
  data nfinp, nfout, nfcntn_bin, nfchgu, nfchgd, nfzaj, nfvlc, nfchgu_p, nfchgd_p, nfvlcr  &
       & /31,     6,     55,     60,     61,    44,    45,  46,  47,48 /
!!$  integer :: nfinp, nfopgr, nfmatbp, nfkpoint, nfkindex, nfgpt &
!!$       &, nfout, nfcntn, nfcntn_bin, nfchgt, nfchgu, nfchgd, nfzaj, nfvlc
!!$  data nfinp, nfopgr, nfmatbp, nfkpoint, nfkindex, nfgpt &
!!$       &, nfout, nfcntn, nfcntn_bin, nfchgt, nfchgu, nfchgd, nfzaj, nfvlc  &
!!$       & /31,21,22,23,24,54 &
!!$       &,  6,42,55,43,60,61,44,45 /

  integer, parameter :: number_of_all_files = 10
  integer, dimension(number_of_all_files) :: n_file
  data n_file &
       & /31,6,55,60,61,44,45,46,47,48 /
!!$       & /31,21,22,23,24,54 &
!!$       &,  6,42,55,43,60,61,44,45 /

  character*80 &
       &  F_INP, F_CNTN_BIN, F_CHGU, F_CHGD, F_ZAJ, F_VLC, F_CHGU_P, F_CHGD_P, F_VLCR &
!!$       &  F_INP, F_OPGR, F_MATBP, F_KPOINT, F_KINDEX, F_GPT &
!!$       &, F_CNTN, F_CNTN_BIN, F_CHGT, F_CHGU, F_CHGD, F_ZAJ, F_VLC &
       &, F_file_names_data
  namelist/fnames/&
       &  F_INP, F_CNTN_BIN, F_CHGU, F_CHGD, F_ZAJ, F_VLC, F_CHGU_P, F_CHGD_P, F_VLCR
!!$       &  F_INP, F_OPGR, F_MATBP, F_KPOINT, F_KINDEX, F_GPT &
!!$       &, F_CNTN, F_CNTN_BIN, F_CHGT, F_CHGU, F_CHGD, F_ZAJ, F_VLC

contains
!!$  subroutine open_nfcntn
!!$    logical open
!!$    inquire(unit = nfcntn, opened = open)
!!$    if(.not.open) &
!!$         & call open0(nfcntn,F_CNTN,'F_CNTN    ',unknown,formatted&
!!$         &     ,check_file_name_on)
!!$  end subroutine open_nfcntn

  subroutine open_nfcntn_bin
    logical open
    inquire(unit = nfcntn_bin, opened = open)
    if(.not.open) &
         & call open0(nfcntn_bin,F_CNTN_BIN,'F_CNTN_BIN',unknown,unformatted&
         &     ,check_file_name_on)
  end subroutine open_nfcntn_bin

  subroutine m_Files_set_default_filenames
    F_file_names_data   = './stm_file_names.data'

    F_INP       = './nfstminput.data'
!c---------- nbztype >= 100
!!$    F_OPGR      = './opgr.data'
!!$    F_MATBP     = './matrix.BP'
!!$    F_KPOINT    = './kpoint.dat'
!!$    F_KINDEX    = './f.kp0'
!c----------- Continue file
!!$    F_GPT       = './nfgpt.data'
!!$    F_CHGT      = './nfchgt.data'
    F_CHGU      = './nfchgu.cube'
    F_CHGD      = './nfchgd.cube'
    F_CHGU_P    = './nfchgu_asis.cube'
    F_CHGD_P    = './nfchgd_asis.cube'
!!$    F_CNTN      = './continue.data'
    F_CNTN_BIN  = './continue_bin_stm.data'
    F_ZAJ       = './zaj.data'
!c---------------------
    F_VLC       = './nfvlc.data'
    F_VLCR      = './nfvlcr.data'
  end subroutine m_Files_set_default_filenames

  subroutine m_Files_read_file_names_data
    integer, parameter ::  nffile = 5
    logical :: ex
    inquire(file=F_file_names_data,exist=ex)
    if (ex) then
      call open0(nffile, F_file_names_data, 'FFILENAMES', old, formatted,check_file_name_on)
      rewind nffile
      read(nffile,NML = fnames, err = 1004, end = 1004)
1004  continue
      close(nffile,status='keep')
    endif

    print *,' F_INP      = ', F_INP
!!$    print *,' F_CHGT     = ', F_CHGT
    print *,' F_CHGU     = ', F_CHGU
    print *,' F_CHGD     = ', F_CHGD
    print *,' F_CHGU_P   = ', F_CHGU_P
    print *,' F_CHGD_P   = ', F_CHGD_P
!!$    print *,' F_CNTN     = ', F_CNTN
    print *,' F_CNTN_BIN = ', F_CNTN_BIN
    print *,' F_ZAJ      = ', F_ZAJ
    print *,' F_VLC      = ', F_VLC
    print *,' F_VLCR     = ', F_VLCR
  end subroutine m_Files_read_file_names_data

  subroutine m_Files_open_files_initially
    print *,' open_files_initially start'
    call open0(nfinp, F_INP, 'F_INP     ',    old,   formatted,check_file_name_on)
!    call open0(nfchgt,F_CHGT,'F_CHGT    ',unknown, unformatted,check_file_name_on)
    call open0(nfchgu,F_CHGU,'F_CHGU    ',unknown,   formatted,check_file_name_on)
    call open0(nfchgd,F_CHGD,'F_CHGD    ',unknown,   formatted,check_file_name_on)
    call open0(nfchgu_p,F_CHGU_P,'F_CHGU_P  ',unknown,   formatted,check_file_name_on)
    call open0(nfchgd_p,F_CHGD_P,'F_CHGD_P  ',unknown,   formatted,check_file_name_on)
    call open0(nfzaj, F_ZAJ, 'F_ZAJ     ',unknown, unformatted,check_file_name_on)
!!$    call open0(nfcntn,F_CNTN,'F_CNTN    ',    old,   formatted,check_file_name_on)
    call open0(nfcntn_bin,F_CNTN_BIN,'F_CNTN_BIN',    old, unformatted,check_file_name_on)
    call open0(nfvlc, F_VLC, 'F_VLC     ',unknown, unformatted,check_file_name_on)
    call open0(nfvlcr,F_VLCR,'F_VLCR    ',unknown, formatted,check_file_name_on)
    print *,' open_files_initially end'
  end subroutine m_Files_open_files_initially

  subroutine open0(nfile,filename,tagname,status, form, write_or_no)
    integer, intent(in) :: nfile
    character*(*), intent(in) :: filename
    character*(10),intent(in) ::  tagname
    integer, intent(in) :: status, form, write_or_no

    if(status==old .and. form==formatted) then
       open(nfile,file=filename, status='old', form='formatted')
    else if(status==old .and. form==unformatted) then
       open(nfile,file=filename, status='old', form='unformatted')
    else if(status==unknown .and. form==formatted) then
       open(nfile,file=filename, status='unknown', form='formatted')
    else if(status==unknown .and. form==unformatted) then
       open(nfile,file=filename, status='unknown', form='unformatted')
    endif

    if(write_or_no == check_file_name_on) then
       write(nfout,9001) tagname, ' = ', filename
    endif
9001 format(' ',a10,a3,a60)
  end subroutine open0

  subroutine m_Files_close_all_files
    integer i
    logical open
    do i = 1, number_of_all_files
       if(n_file(i) == 6) cycle
       inquire(unit = n_file(i), opened = open)
       if(open) then
          close(n_file(i),status='keep')
          write(nfout,*) ' closed filenumber = ', n_file(i)
       endif
    enddo
    close(6, status='keep')
  end subroutine m_Files_close_all_files
end module m_Files
