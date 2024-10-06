!================================================
!  Software name : STM
!  Module : m_ArraySize_Parameters
!  Subroutine(s) : m_ArraySize_Parameters_rd
!  Author(s)     : Takahiro Yamasaki (June 7, 2004)
!
!  FURTHER MODIFICATION: Junichiro  Koga (June 24, 2004)
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

module m_ArraySize_Parameters
! $Id: m_ArraySize_Parameters.f90,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
  implicit none
  integer :: &
       & kng, kngp&
       &,kng1 &
       &,keg ,   kimg &
       &, knl , knm , knn , kid , kfft&
       &, knlp , knmp , knnp , kidp, kfftp &
       &,knv3 &
       &,kspin

contains
  subroutine m_ArraySize_Parameters_rd(nfcntn_bin)
    integer, intent(in) :: nfcntn_bin

    read(nfcntn_bin) kng
    read(nfcntn_bin) kngp
    read(nfcntn_bin) kng1
    read(nfcntn_bin) keg
    read(nfcntn_bin) kimg
    read(nfcntn_bin) knl
    read(nfcntn_bin) knm
    read(nfcntn_bin) knn
    read(nfcntn_bin) kid
    read(nfcntn_bin) knlp
    read(nfcntn_bin) knmp
    read(nfcntn_bin) knnp
    read(nfcntn_bin) kidp
    read(nfcntn_bin) knv3
    read(nfcntn_bin) kspin
#ifdef DEBUG    
    write(0,*) 'kng',kng
    write(0,*) 'kngp',kngp
    write(0,*) 'kng1',kng1
    write(0,*) 'keg',keg
    write(0,*) 'kimg',kimg
    write(0,*) 'knl',knl
    write(0,*) 'knm',knm
    write(0,*) 'knn',knn
    write(0,*) 'kid',kid
    write(0,*) 'knlp',knlp
    write(0,*) 'knmp',knmp
    write(0,*) 'knnp',knnp
    write(0,*) 'kidp',kidp
    write(0,*) 'knv3',knv3
    write(0,*) 'kspin',kspin
#endif

    kfft  = kid *knm *knn *kimg
    kfftp = kidp*knmp*knnp*kimg

  end subroutine m_ArraySize_Parameters_rd
end module m_ArraySize_Parameters
