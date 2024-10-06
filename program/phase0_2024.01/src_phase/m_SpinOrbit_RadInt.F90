module m_SpinOrbit_RadInt
! $Id: m_SpinOrbit_RadInt.F90 633 2020-12-01 05:11:03Z jkoga $
  use m_Const_Parameters,    only : DP, CMPLDP, PAI4, ON, &
       &                             SphericalHarmonicsExpansion, &
       &                             BUILTIN, ByPawPot, ZeffApprox, yes, SKIP, Hartree, InvHyperFineConst
  use m_Files,               only : nfout

  use m_Control_Parameters,  only : nspin ,ndim_magmom, noncol, ndim_chgpot, &
       &                            SpinOrbit_MassCorrection, &
       &                            SpinOrbit_Mode, sw_write_soi_on_atoms

  use m_PseudoPotential,   only : ipaw, nlmt, ilmt, ltp, mtp, taup, &
       &                          lppw, tppw, wf_mnrc, radr_paw, &
       &                          iltpw, psirpw, ival, mmesh, nmesh, xh, &
       &                          lpsmax, nloc, ntau, &
       &                          pot_has_soc, iloc, &
       &                          m_PP_include_vanderbilt_pot, dion_paw

  use m_Ionic_System,          only : ityp, ntyp, natm, iwei, iatomn, &
       &                              scaling_so, magmom_local_now

  use m_PAW_XC_Potential,    only : vxc_ae_k, vxc_ps_k

  use m_PAW_ChargeDensity,    only : paw_dnr, surface_integral_method, &
       &                             m_PAWCD_set_ae_cd_sphex2, &
       &                             m_PAWCD_ae_cd_sphex2_nonclA, &
       &                             m_PAWCD_set_ps_cd_sphex2, &
       &                             m_PAWCD_ps_cd_sphex2_nonclA

  use m_SpinOrbit_Potential,         only : Mat_LS_with_real_ylm_L0, &
       &                             Mat_LS_with_real_ylm_L1, &
       &                             Mat_LS_with_real_ylm_L2, &
       &                             Mat_LS_with_real_ylm_L3, dsoc

  use m_SpinOrbit_Potential,           only :  Mat_SOC_Strength

  use m_PseudoPotential,   only : flg_paw, Mat_SOC_Strength_nonpaw, rhcorpw, rhpcrpw
  use m_Parallelization,  only : mype
  use m_Charge_Density,  only : hsr

! === KT_add === 2014/08/11
  use m_Control_Parameters, only : sw_use_ival_for_paw_ps_soc, &
       &                           sw_use_rphi_Hsoc_rphi
  use m_PseudoPotential,  only : vlocr_pw, phirpw
! ============== 2014/08/11
  use mpi


  implicit none
!  include 'mpif.h'

  integer max_sph_expansion
  parameter( max_sph_expansion = 25 )

contains

  subroutine m_SO_calc_corelevel_splitting( ia, it, num_core_wfns, qnum_l, &
       &                                    num_mesh, psir_core_wfns, &
       &                                    level_splitting )
    integer, intent(in) :: ia, it, num_mesh, num_core_wfns
    integer, intent(in) :: qnum_l( num_core_wfns )
    real(kind=DP), intent(in) :: psir_core_wfns( num_mesh, num_core_wfns )
    real(kind=DP), intent(out) :: level_splitting( num_core_wfns )
!
    integer :: i, ir, ier, nrc, nrmax
    real(kind=DP) :: csum
    real(kind=DP), allocatable :: rho_00_ae(:), vtot_rad_ae(:), dV_rad_ae(:)
    real(kind=DP), allocatable :: rho_00_ps(:), vtot_rad_ps(:), dV_rad_ps(:)
    real(kind=DP), allocatable :: wos(:)

    level_splitting = 0.0d0

    if ( .not. flg_paw ) return

! ---- init ---
    allocate( rho_00_ae( mmesh ) );   rho_00_ae = 0.0d0
    allocate( vtot_rad_ae( mmesh ) ); vtot_rad_ae = 0.0d0
    allocate( dV_rad_ae( mmesh ) )  ; dV_rad_ae = 0.0d0
    allocate( wos( mmesh ) );  wos = 0.0d0

    if ( sw_use_rphi_Hsoc_rphi == ON  ) then
       allocate( rho_00_ps( mmesh ) );   rho_00_ps = 0.0d0
       allocate( vtot_rad_ps( mmesh ) ); vtot_rad_ps = 0.0d0
       allocate( dV_rad_ps( mmesh ) )  ; dV_rad_ps = 0.0d0
    endif

! -----------
    call calc_radial_charge_Y00( ia, it, nrc, rho_00_ae, 1 )      ! AE
    call calc_vh( nmesh(it), radr_paw(:,it), rho_00_ae, vtot_rad_ae )
    call calc_total_potential_Y00( ia, it, nrc, vtot_rad_ae, 1 )     ! AE
    call calc_central_potential( it, nrc, vtot_rad_ae, dV_rad_ae )

    if ( sw_use_rphi_Hsoc_rphi == ON ) then
       call calc_radial_charge_Y00( ia, it, nrc, rho_00_ps, 2 )      ! PS
       call calc_vh( nmesh(it), radr_paw(:,it), rho_00_ps, vtot_rad_ps )
       call calc_total_potential_Y00( ia, it, nrc, vtot_rad_ps, 2 )     ! PS
       call calc_central_potential( it, nrc, vtot_rad_ps, dV_rad_ps )
    endif

    call set_weight_exp( ier, 1, nrc, radr_paw(:,it), wos )

#if 0
    nrmax = nmesh(it)
#else
    nrmax = nrc
#endif

    Do i=1, num_core_wfns
       if ( qnum_l(i) == 0 ) cycle

       csum = 0.0d0
       if ( sw_use_rphi_Hsoc_rphi == ON ) then
          Do ir=1, nrmax
             csum = csum + wos(ir) *( dV_rad_ae(ir) -dV_rad_ps(ir) ) &
                  &                *psir_core_wfns(ir,i)**2
          End do
       else
          Do ir=1, nrmax
             csum = csum + wos(ir) *dV_rad_ae(ir) *psir_core_wfns(ir,i)**2
          End do
       endif
       level_splitting(i) = csum *( 2.0d0*qnum_l(i) +1.0d0 ) /2.0d0
    End do

! --- finalize ---
    deallocate( rho_00_ae ); deallocate( vtot_rad_ae );  deallocate( dV_rad_ae )
    if ( sw_use_rphi_Hsoc_rphi == ON ) then
       deallocate( rho_00_ps ); deallocate( vtot_rad_ps );  deallocate( dV_rad_ps )
    endif
    deallocate( wos )
! -----------

  end subroutine m_SO_calc_corelevel_splitting

  subroutine m_SO_calc_contrib_corelevels( ia, it, num_core_wfns, &
       &                                   orb_index, qnum_l_to_probe, &
       &                                   num_mesh, psir_core_wfns, &
       &                                   e_level, e_so_splitting )
    integer, intent(in) :: ia, it, num_mesh, num_core_wfns
    integer, intent(in) :: orb_index, qnum_l_to_probe
    real(kind=DP), intent(in) :: psir_core_wfns( num_mesh, num_core_wfns )
    real(kind=DP), intent(out) :: e_level, e_so_splitting
!
    integer :: i, ir, ier, nrc, nrmax, ll
    real(kind=DP) :: e_kin, e_pot_ae, e_pot_ps, e_soc, c1, c2
    real(kind=DP), allocatable :: rho_00_ae(:), vtot_rad_ae(:), dV_rad_ae(:)
    real(kind=DP), allocatable :: rho_00_ps(:), vtot_rad_ps(:), dV_rad_ps(:)
    real(kind=DP), allocatable :: vh_core(:)
    real(kind=DP), allocatable :: dpsir_tmp(:), wos(:)

! ---- init ---
    if ( .not. flg_paw ) return

    allocate( vh_core( mmesh ) ); vh_core = 0.0d0

    allocate( rho_00_ae( mmesh ) );   rho_00_ae = 0.0d0
    allocate( vtot_rad_ae( mmesh ) ); vtot_rad_ae = 0.0d0
    allocate( dV_rad_ae( mmesh ) )  ; dV_rad_ae = 0.0d0
    allocate( dpsir_tmp( mmesh ) );  dpsir_tmp = 0.0d0
    allocate( wos( mmesh ) );  wos = 0.0d0

    if ( sw_use_rphi_Hsoc_rphi == ON  ) then
       allocate( rho_00_ps( mmesh ) );   rho_00_ps = 0.0d0
       allocate( vtot_rad_ps( mmesh ) ); vtot_rad_ps = 0.0d0
       allocate( dV_rad_ps( mmesh ) )  ; dV_rad_ps = 0.0d0
    endif

    write(*,*) allocated( vxc_ae_k )
    write(*,*) allocated( vxc_ps_k )
!    stop

#if 0
    Do ir=1, nmesh(it)
       rho_00_ae(ir) = rhcorpw(ir,it) /radr_paw(ir,it)**2 /PAI4
    End Do
    call calc_vh( nmesh(it), radr_paw(:,it), rho_00_ae, vh_core )

    Do ir=1, nmesh(it)
       vlocr_pw(ir,it) = -dble(iatomn(it)) /radr_paw(ir,it) +vh_core(ir)
    End Do
#endif

    call calc_radial_charge_Y00( ia, it, nrc, rho_00_ae, 1 )      ! AE
    call calc_vh( nmesh(it), radr_paw(:,it), rho_00_ae, vtot_rad_ae )
    call calc_total_potential_Y00( ia, it, nrc, vtot_rad_ae, 1 )     ! AE
    call calc_central_potential( it, nrc, vtot_rad_ae, dV_rad_ae )
    
    if ( sw_use_rphi_Hsoc_rphi == ON ) then
       call calc_radial_charge_Y00( ia, it, nrc, rho_00_ps, 2 )      ! PS
       call calc_vh( nmesh(it), radr_paw(:,it), rho_00_ps, vtot_rad_ps )
       call calc_total_potential_Y00( ia, it, nrc, vtot_rad_ps, 2 )     ! PS
       call calc_central_potential( it, nrc, vtot_rad_ps, dV_rad_ps )
    endif
    
    call set_weight_exp( ier, 1, nrc, radr_paw(:,it), wos )

#if 0
    nrmax = nmesh(it)
#else
    nrmax = nrc
#endif

    i = orb_index;  ll = qnum_l_to_probe

    call calc_diff_exp( ier, 4, num_mesh, radr_paw(:,it), psir_core_wfns(:,i), &
         &              dpsir_tmp )
    
    e_kin = 0.0d0;   e_pot_ae = 0.0d0;   e_pot_ps = 0.0d0;   e_soc = 0.0d0
    Do ir=1, nrmax
       c1 = dpsir_tmp(ir)**2
       c2 = ( psir_core_wfns(ir,i)/radr_paw(ir,it) )**2
       e_kin = e_kin + wos(ir) *( c1 + c2 *ll*( ll +1 ) )
    End do
    e_kin = e_kin /2.0d0
    
    if ( sw_use_rphi_Hsoc_rphi == ON ) then
       Do ir=1, nrmax
          e_pot_ae = e_pot_ae + wos(ir) *Vtot_rad_ae(ir) *psir_core_wfns(ir,i)**2
          e_pot_ps = e_pot_ps + wos(ir) *Vtot_rad_ps(ir) *psir_core_wfns(ir,i)**2
          e_soc = e_soc + wos(ir) *( dV_rad_ae(ir) -dV_rad_ps(ir) ) &
               &                  *psir_core_wfns(ir,i)**2
       End do
    else
       Do ir=1, nrmax
          e_pot_ae = e_pot_ae + wos(ir) *Vtot_rad_ae(ir) *psir_core_wfns(ir,i)**2
          e_pot_ps = e_pot_ps + wos(ir) *Vtot_rad_ps(ir) *psir_core_wfns(ir,i)**2
          e_soc = e_soc + wos(ir) *dV_rad_ae(ir) *psir_core_wfns(ir,i)**2
       End do
    endif
    
    call add_kinetic_correction( e_kin, e_pot_ae )

    if ( sw_use_rphi_Hsoc_rphi == ON ) then
!       e_level = e_kin +e_pot_ae -e_pot_ps 
       e_level = e_pot_ae -e_pot_ps 
!       e_level = e_kin +e_pot_ae 
    else
       e_level = e_kin +e_pot_ae
    endif
    e_so_splitting = e_soc *( 2.0d0*ll +1.0d0 ) /2.0d0

! --- finalize ---
    deallocate( rho_00_ae ); deallocate( vtot_rad_ae );  deallocate( dV_rad_ae )
    if ( sw_use_rphi_Hsoc_rphi == ON ) then
       deallocate( rho_00_ps ); deallocate( vtot_rad_ps );  deallocate( dV_rad_ps )
    endif
    deallocate( dpsir_tmp );     deallocate( wos )
! -----------

  contains

    subroutine add_kinetic_correction( ene_kin, ene_pot )
      real(kind=DP), intent(inout) :: ene_kin
      real(kind=DP), intent(in) :: ene_pot

      real(kind=DP) :: HyperFineConst, fac1
      real(kind=DP) :: ene_corr, ene1, ene2

      HyperFineConst = 1.0d0 / InvHyperFineConst
      fac1 = 0.5d0 * HyperFineConst**2

      ene_corr = 0.0d0
      Do ir=1, nmesh(it)
         ene_corr = ene_corr + wos(ir) *Vtot_rad_ae(ir)**2 *psir_core_wfns(ir,i)**2
      End do

      ene1 = ene_kin +ene_pot
      ene2 = ene1**2 -2.0D0 *ene1 *ene_pot + ene_corr

      ene2 = -ene2 /2.d0 /( 511.0*1.0E3/Hartree )
      ene_kin = ene_kin +ene2

    end subroutine add_kinetic_correction

  end subroutine m_SO_calc_contrib_corelevels

  subroutine m_SO_calc_SOC_strength_pawpot

    integer :: ia, it, nrc, mesh_t, ir
    real(kind=DP), allocatable :: rho_00_ae(:), vtot_rad_ae(:), dV_rad_ae(:)
    real(kind=DP), allocatable :: rho_00_ps(:), vtot_rad_ps(:), dV_rad_ps(:)

! ---- init ---
    allocate( rho_00_ae( mmesh ) );   rho_00_ae = 0.0d0
    allocate( vtot_rad_ae( mmesh ) ); vtot_rad_ae = 0.0d0
    allocate( dV_rad_ae( mmesh ) )  ; dV_rad_ae = 0.0d0

    if ( sw_use_rphi_Hsoc_rphi == ON  ) then
       allocate( rho_00_ps( mmesh ) );   rho_00_ps = 0.0d0
       allocate( vtot_rad_ps( mmesh ) ); vtot_rad_ps = 0.0d0
       allocate( dV_rad_ps( mmesh ) )  ; dV_rad_ps = 0.0d0
    endif

! ----------
    Mat_SOC_Strength = 0.0d0

    do ia = 1, natm
       it = ityp(ia)

!       if( ipaw(it)/=1 .or. &
!              & surface_integral_method(it).ne.SphericalHarmonicsExpansion) cycle

       call calc_radial_charge_Y00( ia, it, nrc, rho_00_ae, 1 )      ! AE
       call calc_vh( nmesh(it), radr_paw(:,it), rho_00_ae, vtot_rad_ae )
       call calc_total_potential_Y00( ia, it, nrc, vtot_rad_ae, 1 )     ! AE
       call calc_central_potential( it, nrc, vtot_rad_ae, dV_rad_ae )

       if ( sw_use_rphi_Hsoc_rphi == ON ) then
          call calc_radial_charge_Y00( ia, it, nrc, rho_00_ps, 2 )      ! PS
          call calc_vh( nmesh(it), radr_paw(:,it), rho_00_ps, vtot_rad_ps )
          call calc_total_potential_Y00( ia, it, nrc, vtot_rad_ps, 2 )     ! PS
          call calc_central_potential( it, nrc, vtot_rad_ps, dV_rad_ps )
       endif

       if ( sw_use_rphi_Hsoc_rphi == ON ) then
          call integrate_central_potential_2( it, nrc, dV_rad_ae, dV_rad_ps, &
               &                              Mat_SOC_Strength(:,:,:,ia) )
       else
          call integrate_central_potential_1( it, nrc, dV_rad_ae, &
               &                              Mat_SOC_Strength(:,:,:,ia) )
       endif
    End Do

! --- finalize ---
    deallocate( rho_00_ae ); deallocate( vtot_rad_ae );  deallocate( dV_rad_ae )
    if ( sw_use_rphi_Hsoc_rphi == ON ) then
       deallocate( rho_00_ps ); deallocate( vtot_rad_ps );  deallocate( dV_rad_ps )
    endif
! -----------

  end subroutine m_SO_calc_SOC_strength_pawpot

  subroutine calc_central_potential( it, nrc, vtot_rad, dV_rad )        !  dV/dr *1/r
    integer, intent(in) :: it, nrc
    real(kind=DP), intent(in) :: vtot_rad(mmesh)
    real(kind=DP), intent(out) :: dV_rad(mmesh)
!
    real(kind=DP) :: HyperFineConst
    real(kind=DP) :: fac1, fac2
    integer :: ir, ier
!
    HyperFineConst = 1.0d0 / InvHyperFineConst
!
    dV_rad = 0.0d0
    call calc_diff_exp( ier, 4, nrc, radr_paw(:,it), vtot_rad, dV_rad )
!
    fac1 = 0.5d0 * HyperFineConst**2
!
    if ( SpinOrbit_MassCorrection == 0 ) then
       Do ir=1, nrc
          dV_rad(ir) = fac1 *dV_rad(ir) / radr_paw(ir,it)
       End do

    else if ( SpinOrbit_MassCorrection == 1 ) then
       Do ir=1, nrc
          fac2 = 1.0d0 - vtot_rad(ir) *HyperFineConst **2
          dV_rad(ir) = fac1 *dV_rad(ir) / radr_paw(ir,it) /fac2
       End do

    else if ( SpinOrbit_MassCorrection == 2 ) then
       Do ir=1, nrc
          fac2 = 1.0d0 - vtot_rad(ir) *HyperFineConst **2 /2.0d0
          dV_rad(ir) = fac1 *dV_rad(ir) / radr_paw(ir,it) /fac2
       End do
    endif
!
    if ( allocated( scaling_so ) ) then
       dV_rad(:) = dV_rad(:) *scaling_so(it)
    endif
  end subroutine calc_central_potential

  subroutine calc_total_potential_Y00( ia, it, nrc, vtot_rad, ae_or_ps )
    integer, intent(in) :: ia, it, nrc, ae_or_ps
    real(kind=DP), intent(inout) :: vtot_rad(mmesh)

    integer :: ir, ierr, itry
    real(kind=DP) :: factor, c1, c2, c3
    real(kind=DP) :: f1, f2, f3, r1, r2, r3
    real(kind=DP) :: ctmp1, ctmp2, rtmp1, rtmp2, cn
    
    factor = sqrt(PAI4)
    
    if ( ae_or_ps == 1 ) then              ! AE potential
!      Do ir=1, nrc
       Do ir=1, nmesh(it)
          if ( nspin==2 ) then
             c1 = vxc_ae_k( ir,1,1,ia ) + vxc_ae_k( ir,2,1,ia )
             c1 = c1 / 2.0d0
          else
             c1 = vxc_ae_k( ir,1,1,ia )
          endif
       
          c2 = -dble(iatomn(it)) / radr_paw(ir,it)
          c3 = ( c2 + vtot_rad(ir) ) *factor
          vtot_rad(ir) = c1 + c3
       End do

    else if ( ae_or_ps == 2 ) then          ! PS potential
       Do ir=1, nmesh(it)
          if ( nspin==2 ) then
             c1 = vxc_ps_k( ir,1,1,ia ) + vxc_ps_k( ir,2,1,ia )
             c1 = c1 / 2.0d0
          else
             c1 = vxc_ps_k( ir,1,1,ia )
          endif

          if ( sw_use_ival_for_paw_ps_soc == ON ) then
             c2 = -dble(ival(it)) / radr_paw(ir,it)
          else
             c2 = vlocr_pw(ir,it)
          endif

          c3 = ( c2 + vtot_rad(ir) ) *factor
          vtot_rad(ir) = c1 + c3
       End do
    endif

    vtot_rad = vtot_rad /factor

    return

#ifdef _POT_SMOOTHING_
    itry = 1
    Do while (.true.)
       r1 = radr_paw(nrc,  it)
       r2 = radr_paw(nrc- 5*itry,it)
       r3 = radr_paw(nrc-10*itry,it)
       f1 = vtot_rad( nrc )
       f2 = vtot_rad( nrc -5*itry )
       f3 = vtot_rad( nrc-10*itry )
!
       if ( f2/f1 > 1.0 .and. f3/f1 > 1.0 ) exit
       itry = itry + 1

       write(*,*) 'itry = ', itry
    End Do

    ctmp1 = log( -log( f1/f2)  ); ctmp2 = log( -log( f1/f3 ) )
    rtmp1 = log( r1-r2 ); rtmp2 = log( r1-r3 )

    cn = ( ctmp1 -ctmp2 ) / ( rtmp1 -rtmp2 )
    c1 = -log( f1/f2 ) /(r1-r2)**cn

    write(*,*) 'nrc = ', nrc
    write(*,*) 'ctmp1 ctmp2 = ', ctmp1, ctmp2
    write(*,*) 'rtmp1 rtmp2 = ', rtmp1, rtmp2
    
    write(*,*) 'c1, cn = ', c1, cn
    
    if ( c1 < 0.0 ) c1 = -c1
!
    Do ir=nrc+1, nmesh(it)
       c2 = radr_paw(ir,it) -radr_paw(nrc,it)
       vtot_rad(ir) = vtot_rad(nrc) *exp( -c1*c2**cn )
    End do
#endif

  end subroutine calc_total_potential_Y00
  

  subroutine calc_radial_charge_Y00( ia, it, nrc, rho_00, ae_or_ps )
    integer, intent(in) :: ia, it, ae_or_ps
    integer, intent(out) :: nrc
    real(kind=DP), intent(out) :: rho_00(mmesh)
    
    integer :: ier, ir, itry, nrmax
    integer :: nrc0, dnr, zz, nrc_tmp
    
    integer :: msphmx, msph, msphmx_chg
    integer :: num_isph_chg, isph_chg(max_sph_expansion)
    integer        :: id_sname = -1
!
    real(kind=DP) :: factor
    real(kind=DP) :: c1, c2, f1, f2, f3, r1, r2, r3
    real(kind=DP) :: ctmp1, ctmp2, rtmp1, rtmp2, cn

    real(kind=DP), allocatable :: dummy1(:,:)
    real(kind=DP), allocatable :: dummy2(:,:,:)
    real(kind=DP), allocatable :: rho_work(:,:,:)
    real(kind=DP), allocatable :: wos(:)
    
    allocate( wos(mmesh) ); wos = 0.0d0
   
    factor = sqrt(PAI4)

    dnr = paw_dnr(it)
!      dnr = 1

    rho_00 = 0.0d0

    if (dnr.gt.1) then
       nrc0 = wf_mnrc(it)
       nrc = 1 +int((nrc0-1)/dnr)*dnr
       zz = dble(nrc0-nrc)/dble(dnr)
       nrc = nrc +2*dnr                                                 ! 3rd
       call set_weight_exp3(ier,1,nrc,dnr,radr_paw(:,it),zz,wos)     !  3rd
    else
       nrc=wf_mnrc(it)
       call set_weight_exp(ier,1,nrc,radr_paw(:,it),wos)
    end if
    
    do ir=1,nrc,dnr
       wos(ir) = wos(ir)*radr_paw(ir,it)**2
    end do
    wos = wos*iwei(ia)
    
    msphmx_chg = 0
    msph = max_sph_expansion
    
    if ( noncol ) then
       allocate( dummy1(3,msph) );  allocate( dummy2(nrc,ndim_magmom,msph) )
       allocate( rho_work(mmesh,ndim_magmom,max_sph_expansion) )
       rho_work = 0.0d0

       if ( ae_or_ps == 1 ) then
          call m_PAWCD_ae_cd_sphex2_nonclA( ia, nspin, nrc, dnr, &
               &                            msph, rho_work(1:nrc,:,:), msphmx_chg, &
               &                            num_isph_chg, isph_chg, wos, &
               &                            1, dummy1, dummy2 )
       else if ( ae_or_ps == 2 ) then
          call m_PAWCD_ps_cd_sphex2_nonclA( ia, nspin, nrc, dnr, &
               &                            msph, rho_work(1:nrc,:,:), msphmx_chg, &
               &                            num_isph_chg, isph_chg, wos, &
               &                            1, dummy1, dummy2 )
       endif
       deallocate( dummy1 );  deallocate( dummy2 )

    else
       allocate( rho_work(mmesh,nspin,max_sph_expansion) ); rho_work = 0.0d0
#if 1
       nrmax = nrc
#else
       nrmax = nmesh(it)
#endif

       if ( ae_or_ps == 1 ) then
          call m_PAWCD_set_ae_cd_sphex2( ia, nspin, nrmax, dnr, &
               &                         msph, rho_work(1:nrmax,:,:), msphmx_chg, &
               &                         num_isph_chg, isph_chg )
       else if ( ae_or_ps == 2 ) then
          call m_PAWCD_set_ps_cd_sphex2( ia, nspin, nrmax, dnr, &
               &                         msph, rho_work(1:nrmax,:,:), msphmx_chg, &
               &                         num_isph_chg, isph_chg )
       endif

    endif
         
    if ( noncol ) then
       rho_00(:) = rho_work( :,1,1 ) + rho_work(:,2,1)
    else
       if ( nspin==1 ) then
          rho_00(:) = rho_work(:,1,1)
       else if ( nspin==2 ) then
          rho_00(:) = rho_work(:,1,1) + rho_work(:,2,1)
       endif
    end if
    
    deallocate( rho_work ); deallocate( wos )

    return

#ifdef _POT_SMOOTHING_
! -------------------------
!  in the vicinity of r=nrc, assuming rho = A*exp(-B*R**n )
!  f1 = A*exp(-B*r1**n); f2 = A*exp(-B*r2**n);
!  f1/f2 = exp(-B*(r1-r2)**n );  log(f1/f2) = -B*(r1-r2)**n;
!  f1/f3 = exp(-B*(r1-r3)**n );  log(f1/f3) = -B*(r1-r3)**n;
!
!  log( -log(f1/f2) ) = log(B) + n*log(r1-r2)
!  log( -log(f1/f3) ) = log(B) + n*log(r1-r3)
! -------------------------

    if ( m_PP_include_vanderbilt_pot(it) /= SKIP ) then
       nrc_tmp = nrc
       Do while ( .true. )
          if ( rho_00(nrc_tmp-1) < rho_00(nrc_tmp) ) then
             nrc_tmp = nrc_tmp -1
          else
             exit
          endif
       End Do

       itry = 1
       Do while (.true.)
          r1 = radr_paw(nrc_tmp,  it)
          r2 = radr_paw(nrc_tmp- 5*itry,it)
          r3 = radr_paw(nrc_tmp-10*itry,it)
          f1 = rho_00( nrc_tmp )
          f2 = rho_00( nrc_tmp-5 )
          f3 = rho_00( nrc_tmp-10 )
          !
          if ( f2/f1 > 1.0 .and. f3/f1 > 1.0 ) exit
          itry = itry + 1

!          write(*,*) 'itry = ', itry
       End Do

       ctmp1 = log( -log(f1/f2));  ctmp2 = log( -log(f1/f3) )
       rtmp1 = log( r1-r2 ); rtmp2 = log( r1-r3 )
       
       cn = ( ctmp1 -ctmp2 ) / ( rtmp1 -rtmp2 )

       c1 = -log(f1/f2) /(r1-r2)**cn

!         write(*,*) 'c1, cn = ', c1, cn

       if ( c1 < 0.0 ) c1 = -c1
       !
       Do ir=nrc_tmp+1, nmesh(it)
          c2 = radr_paw(ir,it) -radr_paw(nrc_tmp,it)
          rho_00(ir) = rho_00(nrc_tmp) *exp( -c1*c2**cn )
       End do
       
    endif
#endif

!  ---------
!      if ( SpinOrbit_MassCorrection == 2 ) then
!!!         rho_00 = rho_00 *factor
!      endif
! ---
  end subroutine calc_radial_charge_Y00

  subroutine calc_vh( nmesh, rpos, rho, vh )
    integer, intent(in) :: nmesh
    real(kind=DP), intent(in) :: rpos(mmesh), rho(mmesh)
    real(kind=DP), intent(out) :: vh(mmesh)
    
    integer :: ir, ii, i0, is, j, jr, ier
    real(kind=DP) :: sum1, sum2
    real(kind=DP), allocatable :: wt(:)
    
    allocate( wt (nmesh) ); wt = 0.0d0
    
    do ir = 1,nmesh
       sum1 = 0.d0;  sum2 = 0.d0
       if (ir == 1) then
          sum1 = 0.d0
       else if ((ir >= 2).and.(ir <= 5)) then
          do ii = 2,ir
             i0 = ii-1; is = 1
             call set_open_weight_exp(ier,i0,is,rpos,wt)
             
             do j = 1,4
                sum1 = sum1 +rpos(i0+j*is)**2 *rho(i0+j*is) *wt(i0+j*is)
             end do
          end do
       else
          call set_weight_exp(ier,1,ir,rpos,wt)
          do jr = 1,ir
             sum1 = sum1 + rpos(jr)**2*rho(jr)*wt(jr)
          end do
       end if
       sum1 = sum1 *PAI4 /rpos(ir)
       
       if (ir == nmesh) then
          sum2 = 0.d0
       else if ((ir <= nmesh-1).and.(ir >= nmesh-4)) then
          do ii = ir,nmesh-1
             i0 = ii+1; is = -1
             call set_open_weight_exp(ier,i0,is,rpos,wt)
             do j = 1,4
                sum2 = sum2 -rpos(i0+j*is)**2 *rho(i0+j*is) *wt(i0+j*is)
             end do
          end do
       else
          call set_weight_exp(ier,ir,nmesh,rpos,wt)
          do jr = ir,nmesh
             sum2 = sum2 + rpos(jr)*rho(jr)*wt(jr)
          end do
       end if
       sum2 = sum2 *PAI4
       vh(ir) = sum1 + sum2
    end do
    
    deallocate( wt )
  end subroutine calc_vh

  subroutine integrate_central_potential_1( it, nrc, dV_rad, Mat )
    integer, intent(in) :: nrc, it
    real(kind=DP), intent(in) :: dV_rad(mmesh)
    real(kind=DP), intent(out) :: Mat( nloc, ntau, ntau )
    
    integer :: ier, ir
    integer :: ilt1, ilt2
    integer :: il1, il2, it1, it2
    real(kind=DP) :: csum, tmp1
    real(kind=DP), allocatable :: wos(:)
    
    allocate( wos(mmesh) ); wos = 0.0d0
    call set_weight_exp( ier, 1, nrc, radr_paw(:,it), wos )

    Mat = 0.0d0
    
    Do ilt1=1, iltpw(it)
       il1 = lppw(ilt1,it);   it1 = tppw(ilt1,it)
       Do ilt2=1, iltpw(it)
          il2 = lppw(ilt2,it);  it2 = tppw(ilt2,it)
          if ( il1 /= il2 ) cycle
          
          csum = 0.0d0
          Do ir=1, nrc
             tmp1 = dV_rad(ir) *psirpw( ir, il1, it1, it )  &
                  &            *psirpw( ir, il2, it2, it )  &
                  &            *wos(ir)
             csum  = csum + tmp1
          End do
          Mat( il1, it1, it2 ) = csum
       End do
    End Do

    if ( sw_write_soi_on_atoms == ON ) then
       do il1 = 1, lpsmax(it)
          if ( il1 == iloc(it) ) then
             il2 = il1; it1 = 1; it2 = 1

             csum = 0.0d0
             Do ir=1, nrc
                tmp1 = dV_rad(ir) *psirpw( ir, il1, it1, it )  &
                     &            *psirpw( ir, il2, it2, it )  &
                     &            *wos(ir)
                csum  = csum + tmp1
             End do
             Mat( il1, it1, it2 ) = csum
          endif
       end do
    endif
    
    deallocate( wos )

  end subroutine integrate_central_potential_1

  subroutine integrate_central_potential_2( it, nrc, dV_rad_ae, dV_rad_ps, Mat )
    integer, intent(in) :: nrc, it
    real(kind=DP), intent(in) :: dV_rad_ae(mmesh), dV_rad_ps(mmesh)
    real(kind=DP), intent(out) :: Mat( nloc, ntau, ntau )
    
    integer :: ier, ir
    integer :: ilt1, ilt2
    integer :: il1, il2, it1, it2
    real(kind=DP) :: csum, tmp1, tmp2
    real(kind=DP), allocatable :: wos(:)
    
    allocate( wos(mmesh) ); wos = 0.0d0
    call set_weight_exp( ier, 1, nrc, radr_paw(:,it), wos )

    Mat = 0.0d0
    
    Do ilt1=1, iltpw(it)
       il1 = lppw(ilt1,it);   it1 = tppw(ilt1,it)
       Do ilt2=1, iltpw(it)
          il2 = lppw(ilt2,it);  it2 = tppw(ilt2,it)
          if ( il1 /= il2 ) cycle
          
          csum = 0.0d0
          Do ir=1, nrc
             tmp1 = dV_rad_ae(ir) *psirpw( ir, il1, it1, it )  &
                  &               *psirpw( ir, il2, it2, it )  &
                  &               *wos(ir)
             tmp2 = dV_rad_ps(ir) *phirpw( ir, il1, it1, it )  &
                  &               *phirpw( ir, il2, it2, it )  &
                  &               *wos(ir)
             csum  = csum + tmp1 -tmp2
          End do
          Mat( il1, it1, it2 ) = csum
       End do
    End Do

    if ( sw_write_soi_on_atoms == ON ) then
       do il1 = 1, lpsmax(it)
          if ( il1 == iloc(it) ) then
             il2 = il1; it1 = 1; it2 = 1

             csum = 0.0d0
             Do ir=1, nrc
                tmp1 = dV_rad_ae(ir) *psirpw( ir, il1, it1, it )  &
                     &               *psirpw( ir, il2, it2, it )  &
                     &               *wos(ir)
                tmp2 = dV_rad_ps(ir) *phirpw( ir, il1, it1, it )  &
                     &               *phirpw( ir, il2, it2, it )  &
                     &               *wos(ir)
                csum  = csum + tmp1 -tmp2
             End do
             Mat( il1, it1, it2 ) = csum
          endif
       end do
    endif

    deallocate( wos )

  end subroutine integrate_central_potential_2

  subroutine m_SO_check_mode_Builtin
    integer it
!
    Do it=1, 1
       if ( pot_has_soc(it) ) then
          
          if ( SpinOrbit_Mode /= BUILTIN ) then
             SpinOrbit_Mode = BUILTIN
             write(nfout,*) '** ------------------------- **'
             write(nfout,*) '** SpinOrbit mode is force to set to BUILTIN'
          endif
       else
          if ( SpinOrbit_Mode == BUILTIN ) then
             write(nfout,*) '** ------------------------- **'
             write(nfout,*) &
                  & '** Please use the spin-orbit splitted pseudopotential '
             write(nfout,*) &
                  & '** for atom type ', it
             call phase_error_with_msg(nfout,'pp is not spin-orbit splitted ',__LINE__,__FILE__)
          endif
       endif
    End do

    if ( SpinOrbit_Mode == BUILTIN ) then
       Do it=2, ntyp
          if ( .not. pot_has_soc(it) ) then
             write(nfout,*) '** ------------------------- **'
             write(nfout,*) &
                  & '** Please use the spin-orbit splitted pseudopotential '
             write(nfout,*) &
                  & '** for atom type ', it
             call phase_error_with_msg(nfout,'pp is not spin-orbit splitted ',__LINE__,__FILE__)
          endif
       End do
    end if
!
  end subroutine m_SO_check_mode_Builtin

  subroutine m_SO_check_mode_Pawpot
    integer it
!
    if ( SpinOrbit_Mode == BUILTIN ) return
    if ( SpinOrbit_Mode == ByPawPot ) then
       Do it=1, ntyp
          if ( ipaw(it) /= yes ) then
             write(nfout,*) '** ------------------------- **'
             write(nfout,*) &
                  & '** Please use the paw pseudopotential for atom type ', it
             write(nfout,*) '*** and confirm that "paw = off" is commeted out in the accuracy tag '
             call phase_error_with_msg(nfout,'non-paw potenttial unsupported ',__LINE__,__FILE__)
          endif
       End Do
    endif
!
  end subroutine m_SO_check_mode_Pawpot

  subroutine m_SO_check_mode_Zeff
    integer it
!
    if ( SpinOrbit_Mode == BUILTIN ) return
    if ( SpinOrbit_Mode == ByPawPot ) return
    if ( SpinOrbit_Mode == ZeffApprox ) then

    endif
!
  end subroutine m_SO_check_mode_Zeff

  subroutine m_SO_calc_SOC_strength_zeff
    integer :: ia, it, nrc
    real(kind=DP), allocatable :: vtot_rad(:), dV_rad(:)

    Mat_SOC_Strength = 0.0d0

    if ( flg_paw ) then
       allocate( vtot_rad( mmesh ) ); vtot_rad = 0.0d0
       allocate( dV_rad( mmesh ) )  ; dV_rad = 0.0d0

       do ia = 1, natm
          it = ityp(ia)
          nrc = wf_mnrc(it)
          call calc_central_pot_zeff( it, nrc, dV_rad )
          call integrate_central_potential_1( it, nrc, dV_rad, &
               &                              Mat_SOC_Strength(:,:,:,ia) )
       End Do
       deallocate( vtot_rad ); deallocate( dV_rad )

    else
       Do ia=1, natm
          it = ityp(ia)
          Mat_SOC_Strength(:,:,:,ia) = Mat_SOC_Strength_nonpaw(:,:,:,it)  &
               &                       * scaling_so(it)
       End do
    endif 

  contains

    subroutine calc_central_pot_zeff( it, nrc, dV_rad )        !  dV/dr *1/r
      integer, intent(in) :: it, nrc
      real(kind=DP), intent(out) :: dV_rad(mmesh)
!
      real(kind=DP) :: HyperFineConst
      real(kind=DP) :: fac1, fac2
      integer :: ir, ier
!
      HyperFineConst = 1.0d0 / InvHyperFineConst
!
      dV_rad = 0.0d0
      fac1 = 0.5d0 * HyperFineConst**2

      if ( SpinOrbit_MassCorrection == 0 ) then
         Do ir=1, nrc
            dV_rad(ir) = fac1 *dble(iatomn(it)) / radr_paw(ir,it)**3
         End do

      else if ( SpinOrbit_MassCorrection == 1 ) then
         Do ir=1, nrc
            fac2 = 1.0d0 + HyperFineConst**2 &
                 &         *dble(iatomn(it))/ radr_paw(ir,it) 
            dV_rad(ir) = fac1 *dble(iatomn(it)) / radr_paw(ir,it)**3 /fac2
         End do

      else if ( SpinOrbit_MassCorrection == 2 ) then
         Do ir=1, nrc
            fac2 = 1.0d0 + HyperFineConst**2 /2.0d0 &
                 &         *dble(iatomn(it))/ radr_paw(ir,it) 
            dV_rad(ir) = fac1 *dble(iatomn(it)) / radr_paw(ir,it)**3 /fac2
         End do

      end if
!
      dV_rad(:) = dV_rad(:) * scaling_so(it)
!
    end subroutine calc_central_pot_zeff

  end subroutine m_SO_calc_SOC_strength_zeff

#ifdef EVAL_SOC_PREVIOUS
  subroutine set_size_mesht( it, nrc, rho, mesh_t )
    integer, intent(in) :: it, nrc
    integer, intent(out) :: mesh_t
    real(kind=DP), intent(in) :: rho(mmesh)
    
    integer     :: i
    real(kind=DP), parameter :: CRDAMP = 1.d0
    real(kind=DP), parameter :: CRDIST = 10.d0
    
    mesh_t = nrc
!      return

    do i = 10, nmesh(it)-1
       if ( rho(i) - rho(i+1) > CRDAMP .and. radr_paw(i,it) < CRDIST) then
          mesh_t = i
!            if(iprippex>=1) write(nfout,'(" LMTO pot. r_ws=",i5,f12.6)') i, radr(i)

          write(*,*) 'mesht nrc = ', mesh_t, nrc, nmesh(it)
          
          return
       end if
    enddo
    mesh_t = nmesh(it)

  end subroutine set_size_mesht

  subroutine set_rho_for_poisson_eq( it, nrc, rho ) 
                                     ! 4 *pi *r**2 *rho_lm( r,l=0,ispin=1 )
    integer, intent(in) :: it, nrc
    real(kind=DP), intent(inout) :: rho(mmesh)

    integer :: ir

!      Do ir=1, nrc
    Do ir=1, nmesh(it)
       rho(ir) = PAI4 *radr_paw(ir,it)**2 * rho(ir)
    End Do
    
  end subroutine set_rho_for_poisson_eq

  subroutine calc_Hartree_pot_Y00( it, nsize, mesh_t, radr, rhvr, vvv )    
    integer, intent(in) :: nsize, mesh_t, it
    real(kind=DP), intent(in) :: radr( mmesh )
    real(kind=DP), intent(in) :: rhvr( mmesh )
    real(kind=DP), intent(out) :: vvv( mmesh )
    
    real(kind=DP)       :: s2, rhs, rhs1, bm
    real(kind=DP),allocatable,dimension(:) :: da, db ! d(nsize)
    real(kind=DP), allocatable :: wkx(:), wky(:), wkz(:)
    
    integer             :: i
    
    real(kind=DP) ::   hh
    
    hh = 1.0 / xh(it)
    allocate(wkx(mmesh)); wkx = 0.d0
    allocate(wky(mmesh)); wky = 0.d0
    allocate(wkz(mmesh)); wkz = 0.d0
    
    !+++++++++++++++++++++++++++++++
    allocate(da(nsize)); da = 0.d0
    allocate(db(nsize)); db = 0.d0
    !+++++++++++++++++++++++++++++++
    s2 = dlog(rhvr(2)/rhvr(1))/ hh
    rhs = rhvr(1)
    wkx(1)  = rhs*radr(1)/(s2+1)
    wky(1)  = rhs/s2
    db(1)   = hh*rhs*3.d0
    da(1)   = db(1)*radr(1)
    rhs1    = rhs
    do i = 2,3
       rhs     = rhvr(i)
       wkx(i)  = wkx(i-1) + hh *(rhs*radr(i)+rhs1*radr(i-1))*0.5d0
       wky(i)  = wky(i-1) + hh *(rhs        +rhs1          )*0.5d0
       db(i)   = hh *rhs*3.d0
       da(i)   = db(i)*radr(i)
       rhs1    = rhs
    enddo
    do i = 4,mesh_t
       rhs    = rhvr(i)
       db(4)  = hh *rhs*3.d0
       da(4)  = db(4)*radr(i)
       wkx(i)=(9*wkx(i-1)-wkx(i-3)+da(4)+2.d0*da(3)-da(2))/8.d0
       wky(i)=(9*wky(i-1)-wky(i-3)+db(4)+2.d0*db(3)-db(2))/8.d0
       da(1)  = da(2)
       db(1)  = db(2)
       da(2)  = da(3)
       db(2)  = db(3)
       da(3)  = da(4)
       db(3)  = db(4)
    enddo
    bm               = wky(mesh_t)
     !C--*--COULOMB POTENTIAL RVC
    vvv = 0.d0
    do i = 1,mesh_t
!!!         vvv(i) = wkx(i) + radr(i)*(bm-wky(i))        ! r *Vh(r)
       vvv(i) = wkx(i) + (bm-wky(i))                   ! Vh(r)
    enddo
     !+++++++++++++++++++++++++++++++
    deallocate(da); deallocate(db)
    deallocate(wkx)
    deallocate(wky)
    deallocate(wkz)
    !+++++++++++++++++++++++++++++++

  end subroutine Calc_Hartree_pot_Y00
#endif

end module m_SpinOrbit_RadInt
