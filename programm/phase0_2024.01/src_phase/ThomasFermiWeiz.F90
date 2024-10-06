subroutine ThomasFermiWeiz_loop

  use m_Control_Parameters, only : sw_modified_TFW_functional, &
       &                           sw_spherical_averaged_nonlocal, &
       &                           printable, ipritfwfunc, &
       &                           max_iteration_loop, &
       &                           threshold_exit_loop, &
       &                           threshold_skip_loop, &
       &                           threshold_enter_loop, &
       &                           use_deltaV, use_deltaW, use_averaged_nonlocal, &
       &                           kimg, nspin, sw_potential_mixing, &
       &                           m_CtrlP_rmx_now, noncol, ekmode, sw_hubbard, &
       &                           use_preconditioning, &
       &                           max_iter_linmin => max_iteration_linmin, &
       &                           threshold_linmin, iprinegativecharge

  use m_Const_Parameters,  only : ON, OFF, DP, PAI
  use m_ThomasFermiW_Potential,  only : m_TFW_calc_nonlocal_pot_avg_A, &
       &                                m_TFW_calc_nonlocal_pot_avg_B, &
       &                                m_TFW_CGoptimize_init, &
       &                                m_TFW_CGoptimize_finalize, &
       &                                m_TFW_calc_SD_direction, &
       &                                m_TFW_calc_CG_direction, &
       &                                m_TFW_calc_deltaW, &
       &                                m_TFW_store_Psi, &
       &                                m_TFW_try_new_Psi_and_charge, &
       &                                m_TFW_calc_energy_terms, &
       &                                m_TFW_remove_negative_charge, &
       &                                deltaV, &
       &                                m_TFW_alloc_deltaV, &
       &                                m_TFW_alloc_deltaW, &
       &                                m_TFW_dealloc_deltaV, &
       &                                m_TFW_dealloc_deltaW, &
       &                                m_TFW_init_Psi, &
       &                                Vin_tfw, Vlda_tfw, Vout_tfw, &
       &                                m_TFW_Kerker_mixing, &
       &                                m_TFW_calc_phi_MatH_phi

  use m_Total_Energy,  only : get_local_potential_energy, get_hartree_energy, &
       &                      ene_hartree => ehartr, &
       &                      ene_local  => elocal, m_TE_what_is_edeltb_now
  use m_XC_Potential,  only : ene_exc => exc

  use m_IterationNumbers,    only : iteration_electronic
  use m_Ionic_System,       only : natm
  use m_Electronic_Structure,  only : vlhxc_l, vlhxc_l_old, vlhxcQ, totch, dhub_aimag

  use m_Files,       only : nfout
  use m_Parallelization, only : mype, ista_kngp, iend_kngp
  use m_Charge_Density,  only : m_CD_cp_chgq_to_chgqo, chgq_l
  use m_CD_Mixing,       only : rmxtrc

  use m_IterationNumbers,    only : iteration_electronic, iteration_ionic
  use m_PlaneWaveBasisSet,  only : kg, kg1, kg_tfw
  use m_PseudoPotential,   only : flg_paw
  use m_Potential_Mixing,  only : vnonlocal_i,  m_Pot_add_vnonlocal_to_vlhxcQ
  use m_ES_Intgr_VlhxcQlm,   only : m_ESiVQ_integrate_VlhxcQlm
  use m_ES_NonCollinear,   only : m_ES_update_Mat_dion0_noncl, &
       &                          m_ES_set_Mat_dion_scr_noncl
  use m_Orbital_Population,    only :  m_OP_cp_om_to_ommix

  implicit none

! -----------
! variables
! -----------
  logical :: reset_cg
  integer :: niter, succession
  integer :: max_iter, i, flag, is
  integer :: iprinegativecharge_save

  real(kind=DP) :: zeta_product, zeta_product_old
  real(kind=DP) :: Ene(2), Ene_old, deltaE, E_first
  real(kind=DP) :: ekin_TF, ekin_Weiz, ene_deltaW, ene_deltaV
  real(kind=DP) :: ene_nonlocal_avg, dx
  real(kind=DP) :: psi_MatH_psi, phi_MatH_psi, phi_MatH_phi, theta
  real(kind=DP) :: threshold0, threshold1, threshold2
  real(kind=DP) :: edeltab_per_atom, edelta, zeta_norm, rmxt

  real(kind=DP), allocatable :: chgq_l_save(:,:,:)

! -----------
! parameters 
! -----------
  integer, parameter :: linmin_method = 2
  real(kind=DP), parameter :: min_zeta_norm = 1.0D-8

! --------------------------
!!
!!! initialization
!!
! --------------------------
  if ( sw_modified_TFW_functional == OFF ) return

!-----------------------------------------------------
!  if ( iteration_electronic < 5 ) return
  if ( iteration_electronic < 3 ) return
! --------------------------------------------------

  max_iter   = max_iteration_loop
  threshold0 = threshold_enter_loop
  threshold1 = threshold_skip_loop
  threshold2 = threshold_exit_loop

  edelta = m_TE_what_is_edeltb_now()
  edeltab_per_atom = m_TE_what_is_edeltb_now() / natm
  edeltab_per_atom = abs(edeltab_per_atom )

  if ( edeltab_per_atom > threshold0 ) return
  if ( edeltab_per_atom < threshold1 ) return

! --------------------------
!!
!!! set arrays for TFW 
!!
! --------------------------
  call m_TFW_CGoptimize_init
  if ( use_deltaW ) call m_TFW_alloc_deltaW
  if ( use_deltaV ) call m_TFW_alloc_deltaV

!!!!!!!!!  call m_TFW_remove_negative_charge
  call m_TFW_init_Psi

  if ( use_averaged_nonlocal ) then
     if ( sw_spherical_averaged_nonlocal == OFF ) then
        call m_TFW_calc_nonlocal_pot_avg_A
     else
        call m_TFW_calc_nonlocal_pot_avg_B
     endif
  endif

! 
! ---- Now, chgq_l (vlhxc_l) is the variable set after charge/potential mixing --
!
  if ( sw_potential_mixing == ON ) then
     Vout_tfw = vlhxc_l
  else
     call Renewal_of_Potential
     Vout_tfw = vlhxc_l
  endif

  if ( use_deltaW ) call m_TFW_calc_deltaW
  if ( use_deltaV ) call set_deltaV

  if ( sw_potential_mixing == OFF ) then
     allocate( chgq_l_save(ista_kngp:iend_kngp,kimg,nspin) )
     chgq_l_save = chgq_l
  endif

  iprinegativecharge_save = iprinegativecharge
  iprinegativecharge = 0
! --------------------------
!!
!!!  start scf loop 
!!
! --------------------------
  niter = 0;   succession = 0

  zeta_product_old = 0.0d0;  
  Ene_Old = 0.0d0;  deltaE = 0.0d0;

  vlhxc_l = Vin_tfw
! -----------
!
  do while (.true.)
     call m_TFW_calc_SD_direction( psi_MatH_psi, zeta_product )
     zeta_norm = zeta_product /totch
     if ( zeta_norm < min_zeta_norm ) then
        if ( mype == 0 ) then
           write(nfout,'(A,E12.5,A,E12.5)') "!** zeta_norm ", zeta_norm, &
                &                           " is lower than ", min_zeta_norm
        endif
        goto 1000
     endif

     call m_TFW_store_Psi
     call m_TFW_try_new_Psi_and_charge( 0.0d0 )        ! for numerical accuracy

     call m_TFW_calc_energy_terms( ekin_TF, ekin_Weiz, &
          &                        ene_nonlocal_avg, ene_deltaV, ene_deltaW )
     call Renewal_of_Potential
     call get_local_potential_energy
     call get_hartree_energy


     Ene(1) = ekin_TF +ekin_Weiz +ene_local +ene_hartree +ene_exc &
          &       +ene_nonlocal_avg  +ene_deltaV +ene_deltaW
     deltaE = ( Ene(1) -Ene_Old ) /dble(natm)

     if ( printable .and. ipritfwfunc >=2 ) then
        if ( mype == 0 )then
           write(nfout,*) "**** TFW functional info ***"
           write(nfout,*) "!- iter.   TF        Weiz     local    Hartree     Exc     AvgNonl    dV     dW"
           write(nfout,'(I7,8F10.4)') niter, ekin_TF, ekin_Weiz, &
                &                     ene_local, ene_hartree, &
                &                     ene_exc, ene_nonlocal_avg, ene_deltaV, ene_deltaW

           write(nfout,*) "!  iter.   Total      deltaE"
           write(nfout,'(I7,F15.6,E15.6)') niter, Ene(1), deltaE
           write(nfout,*) "**** *** ***"
        endif
     endif

     if ( abs(deltaE) < threshold2 )  succession = succession  +1
     if ( succession == 2 ) exit 

#if 0
     if ( niter == 0 ) E_first = Ene(1)
     if ( Ene(1) > E_first ) exit
#endif

! ------------
     reset_cg = .false.
!     if ( mod(niter,20) == 0 ) reset_cg = .true.
     if ( niter == 0 ) reset_cg = .true.
     
     call m_TFW_calc_CG_direction( reset_cg, zeta_product, zeta_product_old, &
          &                        phi_MatH_psi )
#if 0
     call check_differential
     stop
#endif

     select case (linmin_method)
     case(1)
        call line_minimization_rough( dx )
     case(2)
        call line_minimization_accurate( dx )
     end select

     if ( printable .and. ipritfwfunc >=2 ) then
        if ( mype == 0 ) write(nfout,*) '!! dx for psi update = ', dx
     endif
     call m_TFW_try_new_Psi_and_charge( dx )

     call Renewal_of_Potential

! -------------
     zeta_product_old = zeta_product
     ene_old = Ene(1)

     niter = niter + 1
     if ( niter == max_iter ) exit
  End do

1000  continue

  call m_TFW_CGoptimize_finalize
  call m_TFW_dealloc_deltaV
  call m_TFW_dealloc_deltaW

  call Renewal_of_Potential

! --------------------------
!!
!!!  set next chgq_l/vlhxc_l 
!!
! --------------------------
  rmxt = m_CtrlP_rmx_now(iteration_electronic, iteration_ionic)

  if ( sw_potential_mixing == ON ) then
     call m_TFW_Kerker_mixing( Vout_tfw, vlhxc_l, rmxt )   
  else
     call m_TFW_Kerker_mixing( chgq_l_save, chgq_l, rmxt )  
  endif

  call recover_higer_G_terms

  if ( sw_potential_mixing == OFF ) then

     call Renewal_of_Potential
  endif

! --------------------------
!!
!!!  set next vlhxcQ and related variables
!!
! --------------------------
  if ( sw_potential_mixing == ON ) then
     call update_nonlocal_for_potmix
  else
     call update_nonlocal_for_chgmix
  endif

  if ( allocated( chgq_l_save ) ) deallocate( chgq_l_save )

  iprinegativecharge = iprinegativecharge_save

contains

  subroutine line_minimization1( x )
    real(kind=DP), intent(out) :: x
    
    real(kind=DP) :: c1, c2
    
    call m_TFW_calc_phi_MatH_phi( phi_MatH_phi )
    c1 = psi_MatH_psi -phi_MatH_phi
    c2 = 2.0d0 *phi_MatH_psi
    x = atan(c2/c1) /2.0d0

  end subroutine line_minimization1
  
  subroutine line_minimization_rough( x )
    real(kind=DP), intent(out) :: x
    
    real(kind=DP) :: x0, c1, c2
    
    x0 = 1.0d0
    
    call m_TFW_try_new_Psi_and_charge( x0 )
    call m_TFW_calc_energy_terms( ekin_TF, ekin_Weiz, &
         &                        ene_nonlocal_avg, ene_deltaV, ene_deltaW )
    call Renewal_of_Potential
    call get_local_potential_energy
    call get_hartree_energy
    Ene(2) = ekin_TF +ekin_Weiz +ene_local +ene_hartree +ene_exc &
         &       +ene_nonlocal_avg +ene_deltaV +ene_deltaW

    c2 = phi_MatH_psi *2.0d0
    c1 = ( Ene(2) - c2 *x0 - Ene(1) ) / x0**2
    x = -c2 /c1/2.0d0
  end subroutine line_minimization_rough
  
  subroutine line_minimization_accurate( x )
    real(kind=DP), intent(out) :: x

    real(kind=DP) :: dx, df, d2f
    real(kind=DP), parameter :: x0 = 0.10d0

    integer :: num

    num = 0;  x = x0;   dx = 0.0d0
    Do while (.true.)
       call calc_df_and_d2f( x, df, d2f )
       dx = -df / abs(d2f);    x = x +dx

       if ( ipritfwfunc >=3 ) then
          if ( mype == 0 ) write(nfout,*) 'xnew = ', x
       endif
       num = num +1
       if ( num > max_iter_linmin ) exit
       if ( abs(dx) < threshold_linmin ) exit
    End do
    if ( dx > x0)  dx =  x0
    if ( dx <-x0 ) dx = -x0

  end subroutine line_minimization_accurate

  subroutine calc_df_and_d2f( x, df, d2f )
    real(kind=DP), intent(in) :: x
    real(kind=DP), intent(out) :: df, d2f
    real(kind=DP) :: dx, x1, x2, x3, e1, e2, e3

    dx = 1.0D-3
    x1 = x -dx;  x2 = x;  x3 = x + dx

    call m_TFW_try_new_Psi_and_charge( x1 )
    call m_TFW_calc_energy_terms( ekin_TF, ekin_Weiz, &
         &                        ene_nonlocal_avg, ene_deltaV, ene_deltaW )
    call Renewal_of_Potential
    call get_local_potential_energy
    call get_hartree_energy
    e1 = ekin_TF +ekin_Weiz +ene_local +ene_hartree +ene_exc &
         &       +ene_nonlocal_avg +ene_deltaV +ene_deltaW

    call m_TFW_try_new_Psi_and_charge( x2 )
    call m_TFW_calc_energy_terms( ekin_TF, ekin_Weiz, &
         &                        ene_nonlocal_avg, ene_deltaV, ene_deltaW )
    call Renewal_of_Potential
    call get_local_potential_energy
    call get_hartree_energy
    e2 = ekin_TF +ekin_Weiz +ene_local +ene_hartree +ene_exc &
         &       +ene_nonlocal_avg +ene_deltaV +ene_deltaW

    call m_TFW_try_new_Psi_and_charge( x3 )
    call m_TFW_calc_energy_terms( ekin_TF, ekin_Weiz, &
         &                        ene_nonlocal_avg, ene_deltaV, ene_deltaW )
    call Renewal_of_Potential
    call get_local_potential_energy
    call get_hartree_energy
    e3 = ekin_TF +ekin_Weiz +ene_local +ene_hartree +ene_exc &
         &       +ene_nonlocal_avg +ene_deltaV +ene_deltaW

    df = ( e3 -e1 ) / (2.0d0*dx)
    d2f = ( e1 +e3 -2.d0 *e2 ) /dx**2

  end subroutine calc_df_and_d2f

  subroutine check_differential
    real(kind=DP) :: x0

    x0 = 1.0D-5
    
    call m_TFW_try_new_Psi_and_charge( x0 )
    call m_TFW_calc_energy_terms( ekin_TF, ekin_Weiz, &
         &                        ene_nonlocal_avg, ene_deltaV, ene_deltaW )
    call Renewal_of_Potential
    call get_local_potential_energy
    call get_hartree_energy
    Ene(2) = ekin_TF +ekin_Weiz +ene_local +ene_hartree +ene_exc &
         &       +ene_nonlocal_avg +ene_deltaV +ene_deltaW
!
    if ( mype == 0 ) then
       write(nfout,*) 'EDiff numerical, analtical =', &
            &                  (Ene(2)-Ene(1))/x0, 2*phi_MatH_Psi
    endif

  end subroutine check_differential

  subroutine update_nonlocal_for_chgmix

    call Renewal_of_Potential

    if (sw_hubbard==ON) then
!       call Renewal_of_Hubbard_Parameters
       call Renewal_of_Hubbard_Potential
    end if

  end subroutine update_nonlocal_for_chgmix

  subroutine update_nonlocal_for_potmix
    call m_ESiVQ_integrate_VlhxcQlm(nfout) ! (lclchh) -> vlhxcQ
    call m_Pot_add_vnonlocal_to_vlhxcQ

    if ( noncol ) then
       if ( flg_paw ) then
          if ( ekmode==ON .or. iteration_electronic >= 1 ) then
             call m_ES_update_Mat_dion0_noncl
          endif
       endif
       call m_ES_set_Mat_dion_scr_noncl( VlhxcQ, vnonlocal_i )
    endif
  end subroutine update_nonlocal_for_potmix

  subroutine recover_higer_G_terms
    integer :: i

    if ( sw_potential_mixing == ON ) then
       Do i=ista_kngp, iend_kngp
          if ( i > kg_tfw ) then
             vlhxc_l(i,:,:) = vout_tfw(i,:,:)
          endif
       End do
    else
       Do i=ista_kngp, iend_kngp
          if ( i > kg_tfw ) then
             chgq_l(i,:,:) = chgq_l_save(i,:,:)
          endif
       End do
    endif
  end subroutine recover_higer_G_terms

  subroutine set_deltaV
    integer :: which_deltaV

    if ( sw_potential_mixing == ON ) then
       which_deltaV = 1
    else
       which_deltaV = 3
    endif

    select case (which_deltaV)
    case (1)
       deltaV = vout_tfw -vlda_tfw     
    case (2)
       deltaV = -( vout_tfw -vlda_tfw )
    case (3)
       deltaV = vout_tfw -vin_tfw  
    case (4) 
       deltaV = -( vout_tfw -vin_tfw )
    end select

  end subroutine set_deltaV

end subroutine ThomasFermiWeiz_loop
