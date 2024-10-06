!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  MODULE: m_Total_Energy
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
!=======================================================================
module m_Total_Energy
! $Id: m_Total_Energy_ep.f90 570 2017-04-21 20:34:50Z yamasaki $
  use m_Force,            only : etotal
  use m_Const_Parameters, only : tag_total_energy
  implicit none

contains
  subroutine m_TE_rd_total_energy(nfcntn)
    integer, intent(in) :: nfcntn
    read(nfcntn,*)
    read(nfcntn,*) etotal
  end subroutine m_TE_rd_total_energy

  subroutine m_TE_wd_total_energy(nfcntn)
    integer, intent(in) :: nfcntn
    print *,' tag_total_energy'
    write(nfcntn,*) tag_total_energy
    write(nfcntn,'(2d24.16)') etotal
  end subroutine m_TE_wd_total_energy
end module m_Total_Energy
