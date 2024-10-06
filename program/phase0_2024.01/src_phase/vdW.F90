#define VDW_ONESHOT_NEW

!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 552 $)
!
!  "First-principles Electronic Structure Calculation Program"
!
!  PROGRAM: vdW-Soler
!
!  AUTHOR(S): Y. Ono
!  
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
!=======================================================================
!
!   The original version of this set of the computer programs "PHASE" was developed by 
!  the members of the Theory Group of Joint Research Center for Atom Technology 
!  (JRCAT), based in Tsukuba, in the period 1993-2001.  
!   Since 2002, this program set had been intensively developed as a part of the following 
!  national projects supported by the Ministry of Education, Culture, Sports, Science and 
!  Technology (MEXT) of Japan; "Frontier Simulation Software for Industrial Science 
!  (FSIS)" from 2002 to 2005, "Revolutionary Simulation Software (RSS21)" from 2006 to 
!  2008. "Research and Development of Innovative Simulation Software (RISS)" from 2008 
!  to 2013. These projects is lead by the Center for Research on Innovative Simulation 
!  Software (CISS), the Institute of Industrial Science (IIS), the University of Tokyo.
!   Since 2013, this program set has been further developed centering on PHASE System 
!  Consortium. 
!   The activity of development of this program set has been supervised by Takahisa Ohno.
!
#ifndef DISABLE_VDWDF
subroutine vdW_scf( nspin, ispin, na, nb, nc, chgr, grad_rho, ecnl, &
     &              dFdrho, dFddrho, version_no, vflag )

  use m_Control_Parameters, only : printable, sw_use_WuGygi_method
  use m_Const_Parameters, only : DP, EXC_ONLY, ON
  use m_Files, only : nfout
  use m_Timing, only : tstatc0_begin,tstatc0_end
  use m_vdWDF, only : initialize_vdwdf_scf,build_theta,vdWdf_core,corrections, &
       &              get_dFdrho_dFddrho,finalize_vdwdf,Ecnl_12,Ecnl_3,Ecnl_3s, &
       &              univol,rinplw, ecnl_vdwdf, oneshot
  implicit none

  integer, intent(in) :: nspin,ispin,na,nb,nc
  integer, intent(in) :: version_no
  integer, intent(in) :: vflag

  real(kind=DP), dimension(na*nb*nc), intent(in) :: chgr,grad_rho
  real(kind=DP), intent(out) :: ecnl
  real(kind=DP), dimension(na,nb,nc), intent(out) :: dFdrho,dFddrho
  integer :: i,i1,i2,i3
  logical, save :: initialized=.false.
  integer :: id_sname=-1

  call tstatc0_begin('vdW_scf ',id_sname,1)
  call initialize_vdwdf_scf( nspin, ispin, na, nb, nc, chgr, grad_rho, version_no )
  call build_theta()

  if ( vflag == EXC_ONLY ) oneshot = .true.
  call vdWdf_core()
  if ( vflag == EXC_ONLY ) oneshot = .false.

  if ( sw_use_WuGygi_method == ON ) then
     ecnl = Ecnl_12
  else
     call corrections()
     ecnl = Ecnl_12 + Ecnl_3 - Ecnl_3s
  endif

  ecnl_vdwdf = ecnl

  call get_dFdrho_dFddrho(na,nb,nc,dFdrho,dFddrho)
  call finalize_vdwdf()
  call tstatc0_end(id_sname)

end subroutine vdW_scf

subroutine vdW_scf_stress( nspin, ispin, na, nb, nc, chgr, grad_rho, cgrad_rho, &
     &                     version_no, input_charge )
  use m_Const_Parameters, only : DP, Partial_Core_Charge
  use m_Timing, only : tstatc0_begin, tstatc0_end
  use m_vdWDF, only : initialize_vdwdf_stress, build_theta, vdWdf_stress_tensor_core, &
       &              finalize_vdwdf, univol,rinplw, s_cnl1, s_cnl2, s_cnl1_pc, s_cnl2_pc
  use m_Files,  only : nfout

  implicit none
  integer, intent(in) :: nspin,ispin,na,nb,nc, input_charge, version_no
  real(kind=DP), intent(in) :: chgr(na*nb*nc), grad_rho(na*nb*nc)
  real(kind=DP), intent(in) :: cgrad_rho(na*nb*nc,3)

  integer :: i,i1,i2,i3
  logical, save :: initialized=.false.
  integer :: id_sname=-1

  call tstatc0_begin('vdW_scf_stress ',id_sname,1)
  call initialize_vdwdf_stress( nspin, ispin, na, nb, nc, chgr, grad_rho, cgrad_rho, &
       &                        version_no )
  call build_theta()
  call vdWdf_stress_tensor_core()
  call finalize_vdwdf()
!
  if ( input_charge == Partial_Core_Charge ) then
     s_cnl1_pc = s_cnl1
     s_cnl2_pc = s_cnl2
  endif
  call tstatc0_end(id_sname)

end subroutine vdW_scf_stress

#ifdef VDW_ONESHOT_NEW
subroutine vdW_oneshot()
  use m_Const_Parameters,  only : Valence_plus_PC_Charge, EXC_ONLY, VXC_AND_EXC
  use m_Control_Parameters, only : len_xctype, xctype, printable, oneshot
  use m_Files,   only : nfout
  use m_Charge_Density,  only : chgq_l
  use m_XC_Potential ,  only : m_XC_cal_potential
  use m_Timing, only : tstatc0_begin,tstatc0_end
  use m_Parallelization, only : mype
  use m_Total_Energy, only : etotal
  use m_vdWDF, only:  ecnl_vdwdf

  implicit none

  integer :: id_sname = -1

  call tstatc0_begin('vdW_oneshot ',id_sname,1)
  if(mype==0) then
     write(nfout,*)
     write(nfout,'(a)')    "** 'oneshot' calculation of the vdW-interaction **"
  endif
  write(nfout,*)

  oneshot = .false.
  call m_XC_cal_potential( nfout, Valence_plus_PC_Charge, chgq_l, EXC_ONLY )
  oneshot = .true.


  etotal = etotal + ecnl_vdwdf

  if(printable) then
     write(nfout,'(a,f20.10,a)') '--> total energy : ',etotal,' hartree'
  endif
  write(nfout,*)

  call tstatc0_end(id_sname)

end subroutine vdW_oneshot

#else
subroutine vdW_oneshot()
  use m_Const_Parameters, only : DP
  use m_Control_Parameters, only : printable,nspin
  use m_Files, only : nfout
  use m_Parallelization, only : mype
  use m_Timing, only : tstatc0_begin,tstatc0_end
  use m_Total_Energy, only : etotal
  use m_vdWDF, only: initialize_vdwdf_oneshot,build_theta,vdWdf_core,corrections,Ecnl_12,Ecnl_3,Ecnl_3s,finalize_vdwdf
  implicit none
  integer :: is
  integer :: id_sname = -1
  real(kind=DP) :: Ecnl
  call tstatc0_begin('vdW_oneshot ',id_sname,1)
  if(mype==0) then
     write(nfout,'(a)')    "** 'oneshot' calculation of the vdW-interaction **"
  endif
  Ecnl=0.d0
  is=1
  if(nspin>1) is=-1

  call initialize_vdwdf_oneshot(is)
  call build_theta()
  call vdWdf_core()
  Ecnl = Ecnl+Ecnl_12
  call corrections()
  Ecnl = Ecnl+Ecnl_3 - Ecnl_3s
  call finalize_vdwdf()
  etotal = etotal + Ecnl
  if(printable) then
     write(nfout,'(a,f20.10,a)') 'vdW energy       : ',Ecnl,  ' hartree'
     write(nfout,'(a,f20.10,a)') '--> total energy : ',etotal,' hartree'
  endif
  call tstatc0_end(id_sname)
end subroutine vdW_oneshot
#endif

#endif
