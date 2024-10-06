!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 590 $)
!
!  SUBROUINE: Ewald_and_Structure_Factor
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
!=======================================================================
subroutine Ewald_and_Structure_Factor
! $Id: Ewald_and_Structure_Factor.F90 590 2018-11-06 04:20:58Z jkoga $
  use m_Const_Parameters, only : ON,OFF,Valence_plus_PC_Charge,VXC_AND_EXC,DP,VDW_DFTD3
  use m_Files,            only : nfout
  use m_PseudoPotential,  only : ival
#ifndef ENABLE_ESM_PACK
  use m_Control_Parameters,only : skip_alloc_phonon,printable,sw_wf_predictor,sw_charge_predictor &
    &                           , m_CtrlP_in_initialization,sw_extrapolate_charge, noncol
#else
  use m_Control_Parameters,only : skip_alloc_phonon &
    &                           , sw_esm, printable,sw_wf_predictor,sw_charge_predictor,sw_extrapolate_charge &
    &                           , m_CtrlP_in_initialization,kimg, noncol
#endif
  use m_Control_Parameters,only : vdw_method,sw_vdw_correction
  use m_PlaneWaveBasisSet,only : kgp,ngabc,kg,gr_l,igfp_l
  use m_Ionic_System, only : num_regions, m_IS_regions
!$$#ifndef PARA3D
  use m_Ionic_System,     only : m_IS_structure_factor, m_IS_ewald &
     &                         , ntyp_vdw, m_IS_vdw,m_IS_vdwdf3,eewald,fxyzew_l,natm &
     &                         , m_IS_get_extpl_factor, iatomn, ityp
!$$#endif
  use m_IterationNumbers, only : iteration_ionic

  use m_Charge_Density,   only : m_CD_predictor_pre,m_CD_predictor_post,chgq_l
  use m_ES_ortho,             only : m_ES_modified_gram_schmidt
  use m_XC_Potential,         only : m_XC_cal_potential,vxc_l
  use m_Electronic_Structure, only : m_ES_energy_eigen_values
  use m_ES_LHXC,              only : m_ESlhxc_potential
  use m_ES_Intgr_VlhxcQlm,    only : m_ESiVQ_integrate_VlhxcQlm
  use m_ES_wf_extrpl,         only : m_ES_wf_extrpl_doit
  use m_Electronic_Structure, only : vloc_esm

! === KT_add === 13.0PP
  use m_ES_LHXC,         only : m_ESlhxc_potential_noncl
  use m_Electronic_Structure,   only : m_ES_energy_eigen_vals_noncl
  use m_ES_Ortho,               only : m_ES_modified_gramschmidt_noncl
! ============== 13.0PP

  use m_FFT,                  only : fft_box_size_CD
  implicit none

  real(kind=DP) :: alpha,beta,rms
  integer :: nextpl,iat
  real(kind=DP) :: eewald_g,eewald_r
  real(kind=DP), allocatable, dimension(:,:) :: forc_g
  real(kind=DP) :: ew_alpha
  integer :: i,j
  real(kind=DP), allocatable, dimension(:,:) :: agauss,bgauss
  real(kind=DP), dimension(2) :: alp,cc
  real(kind=DP) :: ivaltmp,alf
  real(kind=DP), allocatable, dimension(:,:) :: few_esm
  integer :: itpcc,nfftcd

!$$#ifndef PARA3D
  if(.not.m_CtrlP_in_initialization().and.sw_charge_predictor==ON) &
       & call m_CD_predictor_pre(nfout,printable)
  call m_IS_structure_factor(nfout,kgp,ngabc)

  if(.not.m_CtrlP_in_initialization().and.&
       & (sw_charge_predictor==ON.or.sw_wf_predictor==on)) then

     if(sw_extrapolate_charge==ON.or.sw_wf_predictor==ON) &
          &  call m_IS_get_extpl_factor(alpha,beta,rms,nextpl)
     if(sw_charge_predictor==ON) &
          & call m_CD_predictor_post(alpha,beta,rms,nextpl,nfout,printable)
     if(sw_wf_predictor == ON) then
        if ( noncol ) then
        else
           call m_ES_wf_extrpl_doit(alpha,beta,rms,nextpl)
        endif
     endif

     call m_XC_cal_potential(nfout,Valence_plus_PC_Charge,chgq_l, VXC_AND_EXC)  ! -> vxc_l
     if ( noncol ) then
        call m_ESlhxc_potential_noncl(nfout,chgq_l,vxc_l)
        call m_ES_modified_gramschmidt_noncl(nfout)
        call m_ESiVQ_integrate_VlhxcQlm(nfout) ! (lclchh) -> vlhxcQ
        call m_ES_energy_eigen_vals_noncl(nfout)
     else
        call m_ESlhxc_potential(nfout,chgq_l,vxc_l) ! (stlhxc) ->vlhxc_l
        call m_ES_modified_gram_schmidt(nfout)
        call m_ESiVQ_integrate_VlhxcQlm(nfout) ! (lclchh) -> vlhxcQ
        call m_ES_energy_eigen_values(nfout)   ! (eigen0) -> eko_l,neordr
     endif
  endif

#ifdef ENABLE_ESM_PACK
   if(sw_esm==ON)then
      nfftcd = fft_box_size_CD(1,0)*fft_box_size_CD(2,0)*fft_box_size_CD(3,0)
      if(.not.allocated(vloc_esm)) allocate(vloc_esm(nfftcd))
      vloc_esm(:) = dcmplx(0.d0,0.d0)
      allocate(agauss(natm,2))
      allocate(bgauss(natm,2))
      do i=1,natm
         call psbhs0(nfout,nint(iatomn(ityp(i))),ivaltmp,itpcc,alp,cc)  ! -(b_PseudoPotential)
         agauss(i,1) = alp(1)
         agauss(i,2) = alp(2)
         bgauss(i,1) = cc(1)
         bgauss(i,2) = cc(2)
      enddo
      call esm_local_(nfftcd,vloc_esm,natm,2,agauss,bgauss)
      vloc_esm(:) = vloc_esm(:)*0.5d0 !Ry -> Ha
      deallocate(agauss)
      deallocate(bgauss)

      call m_IS_ewald(nfout,kg,gr_l,kgp,ngabc,ival,alf) ! R only
      ew_alpha = (1.0d0/alf)**2
      call esm_ewald_g(ew_alpha,eewald_g)
      eewald = eewald+0.5d0*eewald_g  !Ry -> Ha

      allocate(few_esm(3,natm));few_esm=0.d0
      call esm_force_ew(ew_alpha,few_esm)
      do i=1,natm
         do j=1,3
            fxyzew_l(i,j) = fxyzew_l(i,j) + 0.5d0*few_esm(j,i) !Ry -> Ha
         enddo
      enddo
      deallocate(few_esm)
   else
     call m_IS_ewald(nfout,kg,gr_l,kgp,ngabc,ival)
   endif
#else
  call m_IS_ewald(nfout,kg,gr_l,kgp,ngabc,ival)
#endif
  if(sw_vdw_correction == ON .and. vdw_method == VDW_DFTD3) call m_IS_vdwdf3(nfout)
  if(ntyp_vdw>0.and.vdw_method /= VDW_DFTD3) call m_IS_vdw(nfout)
!$$#endif
  if(num_regions>0) call m_IS_regions(nfout)

  if(skip_alloc_phonon.and. &
 &  (sw_charge_predictor==OFF.and.sw_extrapolate_charge==OFF.and.sw_wf_predictor==OFF)) &
 &   call Initial_Electronic_Structure

end subroutine Ewald_and_Structure_Factor
