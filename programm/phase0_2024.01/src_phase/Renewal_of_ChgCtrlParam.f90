! =================================== added by K. Tagami =========== 5.0
subroutine Renewal_of_Chg_Ctrl_Param

  use m_Const_Parameters,    only  : DP, ON, OFF, YES, WHOLE
  use m_Total_Energy,         only : m_TE_what_is_edeltb_now
  use m_Ionic_System,        only : natm
  use m_Control_Parameters,  only : sw_eval_energy_before_charge, &
       &                            sw_update_charge_total, spin_density_mixfactor
!
  use m_Const_Parameters,  only : MSD, CG, SD, LMSD, LMCG,lmeazyCG, lmmsd
  use m_Control_Parameters,  only :  m_CtrlP_solver_for_WFs_now, intzaj

  use  m_IterationNumbers,    only : iteration, iteration_ionic, iteration_electronic, &
       &                             iteration_unit_cell

  use m_Control_Parameters,       only : sw_update_charge_hsr, eval_energy_before_charge
  use m_Control_Parameters,     only : sw_recomposing, sw_force_simple_mixing
  use m_Control_Parameters,     only : sw_recomposing_hsr, sw_force_simple_mixing_hsr

  use m_Files,               only : nfout
  use m_Control_Parameters,   only : printable

! =========================== KT_Test ============================== 12.5Exp
  use m_Control_Parameters,   only : sw_hybrid_functional, &
       &                             truncate_vxw_updating, sw_update_vxw, &
       &                             edelta_for_hyb_chgfix, edelta_for_hyb_convgd
! ================================================================== 12.5Exp

! ================= KT_add =================== 13.0F
  use m_Control_Parameters,   only : use_hybrid_functional, &
       &                             hybrid_calc_is_active, edelta_change_to_hybrid
! ============================================ 13.0F

  use m_CD_Mag_Moment,   only : sw_monitor_atomcharge,  m_CD_set_rad_cov_default, &
       &                        m_CD_set_rad_cov_now

  use m_Const_Parameters,  only : LINEAR, ABRUPT, STEPWISE

  use m_Crystal_Structure,  only : sw_magnetic_constraint, mag_constraint_lambda, &
       &                           damping_method_mag_constraint, &
       &                           num_intermid_lambda, &
       &                           nmax_intermid_lambda, &
       &                           edelta_change_lambda_first, &
       &                           edelta_change_lambda_last, &
       &                           max_iter_elec_mag_constraint, &
       &                           max_iter_ion_mag_constraint, &
       &                           max_iter_cell_mag_constraint, &
       &                           sw_fix_charge_after_constraint

! ================= KT_add =================== 13.0XX
  use m_Control_Parameters,   only : sw_calc_ekin_density, &
       &                             ekin_density_is_active, ekmode
! ============================================ 13.0XX

! ===== KT_add ============== 13.0U3
  use m_Control_Parameters,   only : smearing_width => width, &
       &                             smearing_width_initial => width_initial, &
       &                             method_change_smearing_width, &
       &                             nmax_intermid_width, num_intermid_width, &
       &                             edelta_change_width_first, &
       &                             edelta_change_width_last, &
       &                             sw_wf_mixing, wf_mixing_is_active, &
       &                             edelta_start_wf_mixing
! =========================== 13.0U3

  use m_Parallelization,  only : mype

  implicit none

  real(kind=DP) :: edeltab_per_atom, edeltb_per_atom, edelta

! -----------------------
!  if ( sw_monitor_atomcharge == ON ) then
!     if ( iteration_ionic >1 .and. iteration_electronic ==1 ) then
!        call m_CD_set_rad_cov_default
!        call m_CD_set_rad_cov_now
!     endif
!  endif

! ------------------
  edelta = m_TE_what_is_edeltb_now()

  edeltab_per_atom = m_TE_what_is_edeltb_now() / natm
  edeltab_per_atom = abs(edeltab_per_atom )
!
!  call update_or_fix_charge
!  call update_or_fix_charge_hardpart
!

  if ( sw_magnetic_constraint == ON ) then
     call update_lambda_mag_constraint
  endif

! =========================== KT_Test ============================== 12.5Exp
  if ( sw_hybrid_functional == ON ) then
     if ( truncate_vxw_updating ) call update_or_fix_charge_hybrid
  endif
! ================================================================== 12.5Exp
    
! ================= KT_add =================== 13.0F
  if ( use_hybrid_functional ) then

     if ( sw_hybrid_functional == ON .and. (.not. hybrid_calc_is_active) ) then
        sw_hybrid_functional = off
     endif
     if ( sw_hybrid_functional == OFF ) then
        if ( edeltab_per_atom < edelta_change_to_hybrid ) then
           sw_hybrid_functional = on
           hybrid_calc_is_active = .true.

           write(nfout,*) '****************'
           write(nfout,*) ' hybrid_calc is turned active now'
           write(nfout,*) '****************'
        endif
     endif
  end if
! ============================================ 13.0F

  if ( eval_energy_before_charge == YES ) then
    call order_energy_chg_construction
  endif

!  call update_spin_density_mixfactor
!   call update_spin_recomposing

! === KT_add ==== 13.0XX
  if ( sw_calc_ekin_density == ON ) then
     if ( edeltab_per_atom < 1.0D-2 .or. ekmode == ON ) then
!     if ( iteration_electronic > 1 ) then

        if ( .not. ekin_density_is_active ) then
           ekin_density_is_active = .true.
           write(nfout,*) "** ekin_density_is turned active"
        endif
     endif
  endif
! =============== 13.0XX

! ====== KT_add ======= 13.0U3
  if ( sw_wf_mixing == ON ) call checkif_wf_mixig_active
!
  call update_smearing_width
! ===================== 13.0U3

contains

  subroutine update_spin_recomposing
    logical, save :: first = .true.
    integer :: nstep_threshold = 15 

    integer :: sw_recomposing0 = OFF
    integer :: sw_force_simple_mixing0 = OFF
!
    if ( first ) then
       sw_recomposing0 = sw_recomposing
       sw_force_simple_mixing0 = sw_force_simple_mixing
       first = .false.
    endif
!
    if ( sw_recomposing0 == ON ) then
       if ( iteration <= nstep_threshold ) then
          sw_recomposing = ON
          sw_force_simple_mixing = ON
          sw_recomposing_hsr = ON
          sw_force_simple_mixing_hsr = ON
       else
          sw_recomposing = OFF
          sw_force_simple_mixing = OFF
          sw_recomposing_hsr = OFF
          sw_force_simple_mixing_hsr = OFF
       endif
    end if
  end subroutine update_spin_recomposing

  subroutine update_spin_density_mixfactor
    logical, save :: first = .true.
    real(kind=DP),save ::  spin_dens_mixfactor_0

    if ( first ) then
       spin_dens_mixfactor_0 = spin_density_mixfactor
       first = .false.
    endif
!
    if ( iteration <= 5 ) then
       spin_density_mixfactor = 4.0
    else
       spin_density_mixfactor = spin_dens_mixfactor_0
    endif
!
    if ( printable ) then
       write(nfout,*) '*** spin_density_mixfacotr now = ', &
            & iteration, spin_density_mixfactor
    end if
! --------------------------------------
  end subroutine update_spin_density_mixfactor

  subroutine order_energy_chg_construction
!
    real(kind=DP) :: edeltb_now
    integer :: isolver, sw_submat
    
    integer :: nsteps_for_transition = 3
    
    edeltb_now = edeltab_per_atom
    isolver = m_CtrlP_solver_for_WFs_now(iteration_electronic,iteration_ionic,intzaj &
         &                              ,edeltb_now, sw_submat)
!
    select case(isolver)
    case (SD,MSD)       
       sw_eval_energy_before_charge = OFF
    case (lmSD, lmMSD, lmCG, lmeazyCG)
       sw_eval_energy_before_charge = OFF
    case default
       if ( iteration <= nsteps_for_transition ) then
          sw_eval_energy_before_charge = OFF
       else
          sw_eval_energy_before_charge = ON
       endif
    end select
!    
    if ( printable ) then
       write(nfout,*) '*-- sw_eval_ene_before_charge now = ',  &
            &               iteration, sw_eval_energy_before_charge
    endif
  end subroutine order_energy_chg_construction

  subroutine update_or_fix_charge
    integer, save :: counter = 0

    real(kind=DP) :: ene_threshold1
    integer :: nsteps_for_fixing_charge = -2

    ene_threshold1 = 1.0E-3
    if ( iteration_electronic == 0 &
         & .or. iteration_electronic <= nsteps_for_fixing_charge ) then
       sw_update_charge_total = ON
    else
       if ( edeltab_per_atom > ene_threshold1 ) then
          sw_update_charge_total = OFF
       else
          sw_update_charge_total = ON
       endif
    endif
!
    if ( iteration <= 1 ) then
       sw_update_charge_total = ON
    else if ( iteration_electronic <= nsteps_for_fixing_charge ) then
       sw_update_charge_total = OFF
    else
       sw_update_charge_total = ON
    endif
! 
  end subroutine update_or_fix_charge

  subroutine update_or_fix_charge_hardpart
    if ( sw_update_charge_total == OFF ) then
       sw_update_charge_hsr = OFF
    else

    endif

  end subroutine update_or_fix_charge_hardpart

! =========================  KT_Test ===================== 12.5Exp
  subroutine update_or_fix_charge_hybrid
    integer :: count1 = 0
    integer :: count2 = 0

    if ( abs(edelta) < edelta_for_hyb_chgfix ) then
       count1 = count1 + 1
       if ( count1 == 3 ) then
          sw_update_charge_total = OFF
          if ( printable ) write(nfout,*) '*----- sw_update_charge_total is set OFF --'
       endif
    endif
       !
    if ( abs(edelta) < edelta_for_hyb_convgd ) then
       count2 = count2 + 1
       if ( count2 == 3 ) then
          sw_update_vxw = OFF
          if ( printable ) write(nfout,*) '*----- sw_update_vxw is set OFF --'
       endif
    endif

!
  end subroutine update_or_fix_charge_hybrid
! ========================================================= 12.5Exp

  subroutine update_lambda_mag_constraint

    integer, parameter :: succession = 3
    real(kind=DP) :: threshold( nmax_intermid_lambda +1 )
!
    integer, save :: count = 0, istep = 1

    logical, save :: First = .true.
    logical, save :: mag_constraint_is_over = .false.
    logical, save :: lambda_is_changed 
    integer, save :: sw_fix_charge_after_constr_org

    real(kind=DP), save :: lambda_org, lambda_old, lambda_00

    integer :: i, nn
    real(kind=DP) :: c1, ratio

    if ( First ) then
       lambda_org = mag_constraint_lambda; 
       lambda_old = mag_constraint_lambda
       sw_fix_charge_after_constr_org = sw_fix_charge_after_constraint
       istep = 1
       First = .false.
    else 
!       if ( iteration_unit_cell > 1 .and. iteration_ionic == 1 &
!            &                       .and. iteration_electronic == 1 ) then
!          mag_constraint_is_over = .false.
!          mag_constraint_lambda = lambda_org
!!          istep = 1;  count = 0
!       endif
       if ( iteration_electronic == 1 ) then
          mag_constraint_is_over = .false.
          mag_constraint_lambda = lambda_org
          sw_fix_charge_after_constraint = sw_fix_charge_after_constr_org
          istep = 1;  count = 0
       endif
    endif

    if ( damping_method_mag_constraint == 0 )  return

    if ( iteration_ionic > max_iter_ion_mag_constraint ) then
       mag_constraint_is_over = .true.
       mag_constraint_lambda = 0.0d0;  
    endif
    if ( iteration_unit_cell > max_iter_cell_mag_constraint ) then
       mag_constraint_is_over = .true.
       mag_constraint_lambda = 0.0d0;  
    endif

    if ( mag_constraint_is_over ) then
       if ( sw_fix_charge_after_constraint == ON ) then
          sw_update_charge_total = OFF
       endif
       return
    endif

    lambda_is_changed = .false.

! ------------------
    select case( damping_method_mag_constraint )

    case( STEPWISE )
       nn = num_intermid_lambda +1
       threshold(1)  = edelta_change_lambda_first
       threshold(nn) = edelta_change_lambda_last

       Do i=2, nn-1
          c1 = threshold(nn) *(i-1) + threshold(1) *(nn-i)
          c1 = c1 /dble(nn-1)
          threshold(i) = c1
       End do

       c1 = abs( edeltab_per_atom )
!    c1 = abs( edelta )
!
       if ( c1 < threshold( istep ) ) then
          count = count + 1
          if ( count == succession *istep ) then
             ratio = 1.0D0 - dble(istep) /dble(nn)
             mag_constraint_lambda = lambda_org *ratio

             istep = istep +1
             lambda_is_changed = .true.
          endif
       endif

       if ( istep > nn ) then
          mag_constraint_is_over = .true.
       endif

    case (ABRUPT)
       if ( iteration_electronic > max_iter_elec_mag_constraint ) then
          mag_constraint_lambda = 0.0d0
          lambda_is_changed = .true.
          mag_constraint_is_over = .true.
       endif

       c1 = abs( edeltab_per_atom )

       if ( c1 <  edelta_change_lambda_last ) then
          count = count + 1
          if ( count == succession ) then
             mag_constraint_lambda = 0.0d0
             mag_constraint_is_over = .true.
             lambda_is_changed = .true.
          endif
       endif

    case (LINEAR)
       ratio = dble( iteration_electronic ) / dble(max_iter_elec_mag_constraint)
       ratio = 1.0D0 -ratio

       if ( iteration_electronic <= max_iter_elec_mag_constraint ) then
          mag_constraint_lambda = lambda_org *ratio
          lambda_is_changed = .true.
       else
          mag_constraint_lambda = 0.0d0
          mag_constraint_is_over = .true.
       endif

       c1 = abs( edeltab_per_atom )
       if ( c1 <  edelta_change_lambda_first ) then
          count = count + 1
          if ( count == succession ) then
             mag_constraint_is_over = .true.
             lambda_is_changed = .true.
          endif
       endif
    end select

    if ( lambda_is_changed ) then
       write(nfout,*) "------ mag_constraint_lambda is changed to ", &
            &                       mag_constraint_lambda
!       write(nfout,*) 'lambda_old = ', lambda_old
    endif

    lambda_old = mag_constraint_lambda
!
  end subroutine update_lambda_mag_constraint

! === KT_add === 13.0U3
  subroutine update_smearing_width

    integer, parameter :: ntimes = 3
    real(kind=DP) :: threshold( nmax_intermid_width +1 )
!
    integer, save :: count = 0, istep = 1

    logical, save :: First = .true.
    logical, save :: width_is_changed 
    real(kind=DP), save :: width_initial, width_final, width_old

    integer :: i, nn
    real(kind=DP) :: c1, ratio

    if ( method_change_smearing_width == 0 )  return
    if ( iteration_ionic > 1 ) return

    if ( First ) then
       width_initial = smearing_width_initial
       width_final = smearing_width

       width_old = smearing_width_initial
       smearing_width = width_initial

       istep = 1
       First = .false.
    endif

    width_is_changed = .false.

! ------------------
    select case( method_change_smearing_width )

    case( STEPWISE )
       nn = num_intermid_width +1
       threshold(1)  = edelta_change_width_first
       threshold(nn) = edelta_change_width_last

       Do i=2, nn-1
          c1 = threshold(nn) *(i-1) + threshold(1) *(nn-i)
          c1 = c1 /dble(nn-1)
          threshold(i) = c1
       End do

       c1 = abs( edeltab_per_atom )
!    c1 = abs( edelta )
!
       if ( c1 < threshold( istep ) ) then
          count = count + 1
          if ( count == ntimes *istep ) then
             ratio = 1.0D0 - dble(istep) /dble(nn)
             smearing_width = width_initial *ratio &
                  &          +width_final *( 1.0d0 -ratio )

             istep = istep +1
             width_is_changed = .true.
          endif
       endif
       if ( istep > nn )  method_change_smearing_width = 0

    end select

    if ( width_is_changed ) then
       write(nfout,*) "------ smearing_width is changed to ", &
            &                       smearing_width
    endif

    width_old = smearing_width

  end subroutine update_smearing_width

  subroutine checkif_wf_mixig_active
    integer, parameter :: ntimes = 3
    real(kind=DP), parameter :: threshold1 = 1.0D1

    integer, save :: count = 0
    logical, save :: GoFlag = .false.

    if ( edelta_start_wf_mixing > threshold1 )  then
       wf_mixing_is_active = .true.
       return
    endif

    if ( .not. GoFlag ) then
       if ( edeltab_per_atom < edelta_start_wf_mixing ) then
          count = count +1
          if ( count == ntimes ) then
             GoFlag = .true.
             write(nfout,*) "** wf_mixing is trurned active"
          endif
       endif
    endif

    if ( GoFlag ) then
       wf_mixing_is_active = .true.
    else
       wf_mixing_is_active = .false.
    endif
    
  end subroutine checkif_wf_mixig_active

end subroutine Renewal_of_Chg_Ctrl_Param
