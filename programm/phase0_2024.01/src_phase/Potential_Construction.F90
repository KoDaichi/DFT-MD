subroutine Potential_Construction

  use m_Const_Parameters,   only : ON, OFF
  use m_Control_Parameters,  only : sw_modified_TFW_functional, sw_potential_mixing, &
       &                            sw_hubbard, use_deltaV, kimg, nspin
  use m_ThomasFermiW_Potential,  only : Vin_TFW, Vlda_TFW

  use m_Parallelization, only : mype, ista_kngp, iend_kngp
  use m_Electronic_Structure,  only : vlhxc_l, vlhxc_l_old, m_ES_cp_vlhxc_to_old
  use m_Orbital_Population,   only : m_OP_cp_om_to_ommix
  use m_Files,  only : nfout

  use m_Potential_Mixing,  only : m_Pot_cp_vnonlocal_to_old

  implicit none

  integer :: i, ri, is

  if ( sw_modified_TFW_functional == OFF .and. sw_potential_mixing == OFF ) return

! ---
  call m_ES_cp_vlhxc_to_old
!  if ( sw_mix_charge_hardpart == ON ) call m_Pot_cp_vnonlocal_to_old
! ----
  call Renewal_of_Potential

  if ( sw_potential_mixing == ON ) then
     if (sw_hubbard==ON) then
        call Renewal_of_OccMat(.false., ON, .false. )  ! hsr --> om
        call m_OP_cp_om_to_ommix( nfout, 0.0d0 )       ! om --> ommix
        call Renewal_of_Hubbard_Parameters
        call Renewal_of_Hubbard_Potential
     endif
  end if

  if ( sw_modified_TFW_functional == ON ) then
     Vin_TFW(ista_kngp:iend_kngp,1:kimg,1:nspin) &
          &    = vlhxc_l_old(ista_kngp:iend_kngp,1:kimg,1:nspin)
     Vlda_TFW(ista_kngp:iend_kngp,1:kimg,1:nspin) &
          &    = vlhxc_l(   ista_kngp:iend_kngp,1:kimg,1:nspin)
  endif

end subroutine Potential_Construction
