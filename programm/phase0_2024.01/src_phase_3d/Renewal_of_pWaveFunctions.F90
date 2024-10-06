!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: Renewal_of_pWaveFunctions
!
!  AUTHOR(S): M. Saito, T. Yamasaki   November/04/2003
!  
!  FURTHER MODIFICATION: M. Saito, March/26/2007
!
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
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
!

logical function pEigenValues_are_Converged()
  use m_Files, only :                   nfout
  use m_IterationNumbers, only :        iteration_positron_wf
  use m_Control_Parameters, only :      ipripositron, pev_max_iteration
  use m_Positron_Wave_Functions, only : m_pWF_pevdff, m_pWF_wd_pev
  implicit none
  pEigenValues_are_Converged = m_pWF_pevdff()

  if(pEigenValues_are_Converged) then
     if(ipripositron >= 1) &
          & write(nfout,'(" Positron Wave Functions have converged at ",i7," -th iteration")') &
          & iteration_positron_wf
     call m_pWF_wd_pev(nfout)
  else if(iteration_positron_wf >= pev_max_iteration) then
     if(ipripositron >= 1) then
        write(nfout,'(" Positron Wave Functions have not converged")')
        write(nfout,'(" iteration_positrion_wf have reached to pev_max_iteration (= ",i5," )")') &
             & pev_max_iteration
     end if
     call m_pWF_wd_pev(nfout)
     pEigenValues_are_Converged = .true.
  end if

end function pEigenValues_are_Converged


integer function what_is_the_pcore_solver(isolver)
  use m_Const_Parameters,  only : lmSD, SD, lmMSD, MSD, lmCG, CG
  integer, intent(in) :: isolver

  if(isolver == lmSD) then
     what_is_the_pcore_solver = SD
  else if(isolver == lmMSD) then
     what_is_the_pcore_solver = MSD
!!$  else if(isolver == lmCG) then
!!$     what_is_the_core_solver = CG
  else
     what_is_the_pcore_solver = MSD
  end if
end function what_is_the_pcore_solver
