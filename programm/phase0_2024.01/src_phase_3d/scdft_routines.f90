function isSCDFT() result (ret)
  use m_Const_Parameters, only : DRIVER_SC_DFT
  use m_Control_Parameters, only : driver
  implicit none
  logical :: ret
  if (driver /= DRIVER_SC_DFT) then
    ret = .false.
    return
  endif
  ret = .true.
end function isSCDFT

subroutine wd_eps0(nfcntn)
   use m_Control_Parameters, only : epsilon0, epsilon0_previous
   use m_Parallelization, only : mype
   implicit none
   integer, intent(in) :: nfcntn
   character(len('epsilon0')), parameter :: tag_epsilon0 = 'epsilon0'
   if(mype==0) then
      write(nfcntn,*) tag_epsilon0
      write(nfcntn,'(2f10.5)') epsilon0,epsilon0_previous
   end if
end subroutine wd_eps0

subroutine rd_eps0(nfcntn)
   use m_Control_Parameters, only : DP, epsilon0, epsilon0_previous
   use m_Parallelization, only : mype,MPI_CommGroup,npes
   use m_Control_Parameters, only : printable, m_CtrlP_set_alpha_exx
   use m_Files, only : nfout
   use mpi
   implicit none
!   include 'mpif.h'
   integer, intent(in) :: nfcntn
   logical             :: EOF_reach, tag_is_found
   integer             :: ierr
   real(kind=DP)       :: alp
   integer, parameter  :: len_str = 132
   character(len('epsilon0')), parameter :: tag_epsilon0 = 'epsilon0'
   character*(len_str) :: str
   if(mype == 0)then
      call rewind_to_tag0(nfcntn,len(tag_epsilon0),tag_epsilon0 &
           &, EOF_reach, tag_is_found, str,len_str)
      if(.not.tag_is_found) then
         stop ' tag_epsilon0 is not found'
      else
         read(nfcntn,*) epsilon0,epsilon0_previous
      endif
   endif
   if(npes>1) then
      call mpi_bcast(epsilon0,1,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(epsilon0_previous,1,mpi_double_precision,0,MPI_CommGroup,ierr)
   endif
   if(epsilon0>1e-8) then
      if(printable) &
      & write(nfout,'(a,2f10.5)') ' !!** read epsilon0 and epsilon0_previous : ',epsilon0,epsilon0_previous
      alp = 1.d0/epsilon0
      call m_CtrlP_set_alpha_exx(nfout,alp)
   endif
end subroutine rd_eps0

