!================================================
!  Software name : STM
!  Subroutine(s) : Initilization, aavers
!  Author(s)     : Takahiro Yamasaki (June 7, 2004)
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

subroutine Initialization
! $Id: Initialization.F90,v 1.3 2004/06/26 11:43:44 yamasaki Exp $
  use m_Files       ,only : nfout &
       &                  , m_Files_set_default_filenames &
       &                  , m_Files_read_file_names_data &
       &                  , m_Files_open_files_initially
  use m_Timing      ,only : tstatc_init
  implicit none

!    -----------
  call m_Files_set_default_filenames()
  call aavers                            ! -(here)
  call tstatc_init                       ! -(m_Timing)
  call m_Files_read_file_names_data()
  call m_Files_open_files_initially()
contains
  subroutine aavers
    character*72 vers, system
    vers = 'STM program version 2024.01'
#ifdef VPP
    system = '@(#)system=vpp'
#elif DEC
    system = '@(#)system=dec'
#elif HP
    system = '@(#)system=hp'
#elif SUN
    system = '@(#)system=sun'
#elif ONYX
    system = '@(#)system=onyx'
#elif IRIX64
    system = '@(#)system=irix64'
#elif CRAY
    system = '@(#)system=crayxmp'
#elif HIUX
    system = '@(#)system=hi-ux'
#elif SX
    system = '@(#)system=sx'
#elif SP2
    system = '@(#)system=aix'
#elif Linux
#ifdef PGI
    system = '@(#)system=linux_pgi'
#else
    system = '@(#)system=linux'
#endif
#endif
    write(nfout,*) vers
    write(nfout,*) system
  end subroutine aavers

end subroutine Initialization
