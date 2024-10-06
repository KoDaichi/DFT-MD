!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  MODULE: m_Const_Parameters
!
!  AUTHOR(S): T. Yamasaki and K. Mae   August/20/2003
!  
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
!=======================================================================
module m_ErrorMessages
! $Id: m_ErrorMessages.F90 570 2017-04-21 20:34:50Z yamasaki $
  implicit none

  integer,parameter :: NO_ERROR                   = 0
! ---
  integer,parameter :: FORMAT_ERROR               = 1101
  integer,parameter :: EOF_REACHED                = 1103
! ---
  integer,parameter :: FILE_NOT_EXIST             = 1200
  integer,parameter :: FILENAMES_NOT_EXIST        = 1201
  integer,parameter :: CONT_FILES_NOT_EXIST       = 1202
  integer,parameter :: F_POT_FILE_NOT_EXIST       = 1203
  integer,parameter :: F_ZAJ_FILE_NOT_EXIST       = 1206
  integer,parameter :: F_CHGT_FILE_NOT_EXIST      = 1207
  integer,parameter :: F_CNTN_FILE_NOT_EXIST      = 1208
  integer,parameter :: FILENAMES_FORMAT_ERROR     = 1209
  integer,parameter :: FILENAMES_FORMAT_ERROR_NEB = 1210
  integer,parameter :: ERROR_IN_INPUTFILE_OPENING = 1211
  integer,parameter :: F_CHGT_FILE_NOT_EXIST_EK   = 1212
  integer,parameter :: INVALID_CHARGE_MIXING      = 1302
  integer,parameter :: INVALID_ATOMIC_NUMBER      = 1401
  integer,parameter :: PARALLELIZATION_INVALID_2D = 2101
  integer,parameter :: PARALLELIZATION_INVALID_3D = 2102
  integer,parameter :: PARALLELIZATION_INVALID_NK = 2103
  integer,parameter :: PARALLELIZATION_INVALID_NE = 2104
! ---
  integer,parameter :: len_1101_error = 37
  character(len=len_1101_error) :: msg_format_error = "Input data format error when reading."
!!$                                                    12345678901234567890123456789012345678901234567890  
  integer,parameter :: len_1103_error = 31
  character(len=len_1103_error) :: msg_1103_error   = "Encounter the EOF when reading."
!!$                                                    12345678901234567890123456789012345678901234567890  
  integer,parameter :: len_1206_error = 48
  character(len=len_1206_error) :: msg_1204_error   = "The file for initial wavefunctions is not exist."
!!$                                                    12345678901234567890123456789012345678901234567890  
  integer,parameter :: len_1207_error = 49
  character(len=len_1207_error) :: msg_1205_error   = "The file for initial charge density is not exist."
!!$                                                    12345678901234567890123456789012345678901234567890  
  integer,parameter :: len_1211_error = 35
  character(len=len_1211_error) :: msg_1211_error   = "Error in opening of the inputfile: "
!!$                                                    12345678901234567890123456789012345678901234567890  
  integer,parameter :: len_1212_error = 46
  character(len=len_1212_error) :: msg_1212_error   = "Fixed_charge calc require charge-density file."
!!$                                                    12345678901234567890123456789012345678901234567890  
!  integer,parameter :: len_1208_error = 35
!  character(len=len_1208_error) :: msg_1208_error     = "Error in opening of the inputfile: "
! ---
  integer,parameter :: len_6100_error = 29
  character(len=len_6100_error) :: msg_6100_error     = "A CPP defintion is not proper"
  integer,parameter :: CPP_DEFINE_ERROR_1 = 6101
  integer,parameter :: len_6101_error = 31
  character(len=len_6101_error) :: msg_6101_error     = "TRANSPOSE_WITHOUT_REARRANGEMENT"
!!$                                                      1234567890123456789012345678901234567  

contains
  subroutine m_EMsg_abort01()
  end subroutine m_EMsg_abort01

  subroutine m_EMsg_Warning(iwarning, nfout)
    implicit none
    integer,intent(in) :: iwarning, nfout

    select case (iwarning)
    case (FILENAMES_NOT_EXIST)
       write(*,*) '### Warning(1201): "file_names.data" does not exist.'
       write(*,*) '                   A default set of file names is adapted.'
    case default
       write(nfout,*) '### Warning(9999)'
    end select
  end subroutine m_EMsg_Warning

end module m_ErrorMessages
  
