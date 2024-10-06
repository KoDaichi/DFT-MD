
subroutine UnitaryTransform_WF(imode,ik,zaj_l,update_fsr,wf_mode)
  use m_Const_Parameters, only : DP,ORTHONORMALIZATION
  use m_PlaneWaveBasisSet, only : kg1
  use m_Parallelization, only : np_e,ista_k,iend_k
  use m_Control_Parameters, only : kimg
  use m_ES_WF_by_submat, only : m_ESsubmat_utransform_wf
  use m_Files, only : nfout
  use m_ES_ortho, only : m_ES_MGS_4_each_k
  implicit none
  integer, intent(in) :: imode,ik
  real(kind=DP), dimension(kg1,np_e,ista_k:iend_k,kimg),intent(inout) :: zaj_l
  logical, intent(in) :: update_fsr,wf_mode
  call m_ESsubmat_utransform_wf(imode,ik,zaj_l,update_fsr,wf_mode)
  call m_ES_MGS_4_each_k(nfout,ik,ORTHONORMALIZATION)
end subroutine UnitaryTransform_WF

