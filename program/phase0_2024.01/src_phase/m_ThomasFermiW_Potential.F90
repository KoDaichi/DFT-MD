! --------------------------------------------------------
!
!  hybridized_kinetic_energy_fuctional
!
!   T_kin[rho] = a *T_TF[rho] + b *T_vW[rho];
!                             TF = Thomas Fermi,  vW = von Weizsacker
!
!    a: weight_TF_functional,  b: weight_Weiz_functional
!
! --------------------------------------------------------

module m_ThomasFermiW_Potential

  use m_Const_Parameters,  only : DP, CMPLDP, zi, PAI4, PAI2, PAI, ON, OFF

  use m_Control_Parameters,  only : kimg, nspin, sw_hubbard, &
       &                            use_averaged_nonlocal, use_deltaW, use_deltaV, &
       &                            weight_TF_functional, weight_Weiz_functional, &
       &                            ipritfwfunc, amin, use_preconditioning

  use m_Parallelization, only :  ista_kngp,iend_kngp, ista_fftp, iend_fftp, &
       &                         nel_fftp, nis_fftp, nie_fftp, ista_fftph, iend_fftph, &
       &                         ierr, MPI_CommGroup, mype, npes, mp_fftp, idisp_fftp, &
       &                         ista_fftp, iend_fftp, ista_k

  use m_Ionic_System,    only : ntyp, natm, ityp, pos, iatomn
  use m_Crystal_Structure,  only : rltv, univol
  use m_PseudoPotential,  only : n_non0_lmtxlmt, m_PP_find_maximum_l, dion, dion_paw, &
       &                         ipaw, il2p, dl2p, isph, mmesh, nmesh, ltp, taup, &
       &                         index_lmt1_lmt2, psir_val, radr, wos, mtp, ilmt, &
       &                         nloc, iloc, ntau, itau, lpsmax
  use m_PlaneWaveBasisSet, only : ngabc, gr_l, igfp_l, kgp_reduced, kg1, kg_tfw, &
       &                          gr_l
  use m_NonLocal_Potential, only : new_radr_and_wos
  use m_Electronic_Structure, only : dhub, vlhxc_l, vlhxc_l_old, totch, efermi

  use m_Charge_Density,  only : chgq_l, hsr, hsro
  use m_CD_Mag_Moment,   only : rad_cov_default, sw_monitor_atomcharge, &
       &                        m_CD_set_rad_cov_default

  use m_FFT,                   only : m_FFT_CD_direct_c, m_FFT_CD_inverse_c, &
       &                              m_FFT_alloc_CD_box, m_FFT_dealloc_CD_box, &
       &                              fft_box_size_CD, nfftp

  use m_Files,     only : nfout

  use m_Potential_Mixing,  only :   c_p_pot => c_p, precon_4_pot_mix
  use mpi

  implicit none
!  include 'mpif.h'

! --------------------------------------------------------
! nonlocal avgeraged potential
!
  integer, allocatable :: Table_for_ChargeDensityBasisFn(:,:,:,:)
  real(kind=DP), allocatable :: ChargeDensityBasisFn(:,:,:)
  real(kind=DP), allocatable :: vnonlocal_avg(:,:,:)
! --------------------------------------------------------
! Thomas Fermi von Weizsacker functional
!
  real(kind=DP) :: Coeff_TF
!
  real(kind=DP), allocatable :: Psi(:,:,:), Phi(:,:,:), Phi_Old(:,:,:)
  real(kind=DP), allocatable :: Zeta(:,:,:), deltaW(:,:,:), deltaV(:,:,:)
!
  real(kind=DP), allocatable :: Psi_Old(:,:,:), Vkin_TF(:,:,:)
  real(kind=DP), allocatable :: Vin_tfw(:,:,:), Vlda_tfw(:,:,:), Vout_tfw(:,:,:)
! --------------------------------------------------------
contains

! --------------------------------------------------------
!
! The "atomic wave funcions" Psir_val are normalized.
!
! --------------------------------------------------------
  subroutine m_TFW_normalize_psir_val
    integer :: it, ik, ilmt1, il1, tau1, ir
    real(kind=DP) :: csum, factor, rcut
    real(kind=DP), allocatable :: tmp_fn(:)

    allocate( tmp_fn(mmesh) );  allocate( radr(mmesh) );  allocate( wos(mmesh) )

    Do it=1, ntyp
       ik = ista_k
       call new_radr_and_wos(ik,it)                 ! --> radr, wos

       rcut = rad_cov_default( nint(iatomn(it)) ) ! Revised according to a report from ASMS Co.ltd, 10 March 2016.

       Do il1=1, lpsmax(it)
          if ( il1 == iloc(it) ) cycle

          Do tau1=1, itau(il1,it)
             tmp_fn = 0.0d0
             Do ir=1, nmesh(it)
                tmp_fn(ir) = psir_val(ir,il1,tau1,it) *exp( -radr(ir)/rcut)
             End do

             csum = 0.0d0
             Do ir=1, nmesh(it)
                csum = csum + wos(ir) *tmp_fn(ir)**2
             End do

             factor = 1.0d0 /sqrt(csum)

             psir_val(:,il1,tau1,it) = 0.0d0
             do ir = 1,nmesh(it)
                psir_val(ir,il1,tau1,it) = tmp_fn(ir) *factor
             end do
          End Do
       End Do
    End Do

    call mpi_barrier( MPI_CommGroup, ierr )
    deallocate(tmp_fn); deallocate(radr);  deallocate(wos)

  end subroutine m_TFW_normalize_psir_val

! --------------------------------------------------------
!
! ChargeDensityBasisFn is Fourier transform of
!
! (A-type) :  Psi_lt(r)*Y_lm  *Psi_lt'(r)*Y_lm'
! (B-type) :  Psi_lt(r)*Y_00  *Psi_lt'(r)*Y_00
!                           ( spherically averaged )
!
! --------------------------------------------------------
  subroutine m_TFW_alloc_ChgDensityBasisFn_A
    integer :: num, it

    num = 0
    Do it=1, ntyp
       num = num +n_non0_lmtxlmt(it)
    End do

    allocate( ChargeDensityBasisFn( ista_kngp:iend_kngp, num, 2 ) )
    ChargeDensityBasisFn = 0.0d0
  end subroutine m_TFW_alloc_ChgDensityBasisFn_A

  subroutine m_TFW_alloc_ChgDensityBasisFn_B
    integer :: num, it, il1, tau1, tau2

    allocate( Table_for_ChargeDensityBasisFn( ntyp, nloc, ntau, ntau ) )
    Table_for_ChargeDensityBasisFn = 0

    num = 0
    Do it=1, ntyp
       Do il1=1, lpsmax(it)
          if ( il1 == iloc(it) ) cycle
          do tau1=1, itau(il1,it)
             do tau2=1, itau(il1,it)
                num = num +1
                Table_for_ChargeDensityBasisFn(it,il1,tau1,tau2) = num
             end do
          End do
       End do
    End do

    allocate( ChargeDensityBasisFn( ista_kngp:iend_kngp, num, 1 ) )
    ChargeDensityBasisFn = 0.0d0
  end subroutine m_TFW_alloc_ChgDensityBasisFn_B

  subroutine m_TFW_dealloc_ChgDensityBasisFn
    deallocate( ChargeDensityBasisFn )
    if ( allocated(Table_for_ChargeDensityBasisFn) ) then
       deallocate( Table_for_ChargeDensityBasisFn )
    endif
  end subroutine m_TFW_dealloc_ChgDensityBasisFn

  subroutine m_TFW_set_ChgDensityBasisFn_A
    integer :: lmax, n, ig, ir, num
    integer :: it, ip, ik, nspher
    integer :: ilmt1, ilmt2, il1, il2, tau1, tau2, im1, im2
    real(kind=DP) :: fac, facr, csum
    complex(kind=CMPLDP) :: zfact

    integer, allocatable :: il3(:)
    real(kind=DP), allocatable :: qx(:), qy(:), qz(:), vlength(:)
    real(kind=DP), allocatable :: ylm(:), wka(:), wkb(:)
    complex(kind=CMPLDP), allocatable :: zsnl2(:)

! init ----
    call m_PP_find_maximum_l( lmax )
    n = 2*lmax -1
    allocate(il3(n**2)); call substitute_il3(n**2,il3)
!
    call alloc_arrays

    Do ig=ista_kngp, iend_kngp
       qx(ig) = rltv(1,1)*ngabc(ig,1) +rltv(1,2)*ngabc(ig,2) +rltv(1,3)*ngabc(ig,3)
       qy(ig) = rltv(2,1)*ngabc(ig,1) +rltv(2,2)*ngabc(ig,2) +rltv(2,3)*ngabc(ig,3)
       qz(ig) = rltv(3,1)*ngabc(ig,1) +rltv(3,2)*ngabc(ig,2) +rltv(3,3)*ngabc(ig,3)
       vlength(ig) = sqrt( qx(ig)**2 + qy(ig)**2 + qz(ig)**2 )
    End do

! begin --
!!    fac = PAI4 / univol     !! ??
    fac = sqrt( PAI4 )/ univol

    num = 0
    Do it=1, ntyp
       ik = ista_k
       call new_radr_and_wos(ik,it)                 ! --> radr, wos

       Do ip=1, n_non0_lmtxlmt(it)
          ilmt1 = index_lmt1_lmt2(ip,it,1);  ilmt2 = index_lmt1_lmt2(ip,it,2)

          il1 = ltp(ilmt1,it);  il2 = ltp(ilmt2,it)
          tau1 = taup(ilmt1,it);  tau2 = taup(ilmt2,it)
!
          num = num +1;    zsnl2 = 0.d0

          Do n=1, il2p( ilmt1,ilmt2, it )
             nspher = isph(ilmt1,ilmt2,n,it)
             zfact = zi**il3( nspher) *fac

             call sphr( iend_kngp -ista_kngp +1, nspher, &
                  &     qx(ista_kngp:iend_kngp), &
                  &     qy(ista_kngp:iend_kngp), &
                  &     qz(ista_kngp:iend_kngp), ylm(ista_kngp:iend_kngp) )
                                                        ! -(bottom_Subr.)
             do ir = 1,nmesh(it)
                facr = wos(ir) *psir_val(ir,il1,tau1,it) &
                     &         *psir_val(ir,il2,tau2,it)

                Do ig=ista_kngp, iend_kngp
                   wka(ig) = vlength(ig) *radr(ir)
                End do

                call dsjnv( il3(nspher), iend_kngp -ista_kngp +1, &
                     &      wka(ista_kngp:iend_kngp), &
                     &      wkb(ista_kngp:iend_kngp) )
                                                     ! -(bottom_Subr.)
                Do ig=ista_kngp, iend_kngp
                   zsnl2(ig) = zsnl2(ig) + zfact *facr *wkb(ig) *ylm(ig) &
                        &                        *dl2p(ilmt1,ilmt2,n,it)
                End do
             end do
          End Do

          Do ig=ista_kngp, iend_kngp
             ChargeDensityBasisFn( ig, num, 1 ) = real( zsnl2(ig) )
             ChargeDensityBasisFn( ig, num, 2 ) = aimag( zsnl2(ig) )
          End do
       End Do
    End Do

    call dealloc_arrays;   deallocate( il3 )

  contains

    subroutine alloc_arrays
      allocate( radr(mmesh) ); allocate( wos(mmesh) )
      allocate( qx(ista_kngp:iend_kngp) );
      allocate( qy(ista_kngp:iend_kngp) )
      allocate( qz(ista_kngp:iend_kngp) );
      allocate( ylm(ista_kngp:iend_kngp) ); allocate( vlength(ista_kngp:iend_kngp) )
      allocate( wka(ista_kngp:iend_kngp) );  allocate( wkb(ista_kngp:iend_kngp) )
      allocate( zsnl2(ista_kngp:iend_kngp) )
    end subroutine alloc_arrays

    subroutine dealloc_arrays
      deallocate( radr ); deallocate( wos )
      deallocate( qx ); deallocate( qy );  deallocate( qz ); deallocate( vlength )
      deallocate( wka ); deallocate( wkb ); deallocate( ylm )
      deallocate( zsnl2 )
    end subroutine dealloc_arrays

  end subroutine m_TFW_set_ChgDensityBasisFn_A

  subroutine m_TFW_set_ChgDensityBasisFn_B
    integer :: lmax, n, ig, ir, num
    integer :: it, ip, ik, nspher
    integer :: il1, tau1, tau2
    real(kind=DP) :: fac, facr, csum

    integer, allocatable :: il3(:)
    real(kind=DP), allocatable :: qx(:), qy(:), qz(:), vlength(:)
    real(kind=DP), allocatable :: ylm(:), wka(:), wkb(:)
    real(kind=DP), allocatable :: snl2(:)

! init ----
    call alloc_arrays

    Do ig=ista_kngp, iend_kngp
       qx(ig) = rltv(1,1)*ngabc(ig,1) +rltv(1,2)*ngabc(ig,2) +rltv(1,3)*ngabc(ig,3)
       qy(ig) = rltv(2,1)*ngabc(ig,1) +rltv(2,2)*ngabc(ig,2) +rltv(2,3)*ngabc(ig,3)
       qz(ig) = rltv(3,1)*ngabc(ig,1) +rltv(3,2)*ngabc(ig,2) +rltv(3,3)*ngabc(ig,3)
       vlength(ig) = sqrt( qx(ig)**2 + qy(ig)**2 + qz(ig)**2 )
    End do

! begin --
!!    fac = PAI4 / univol    !???
    fac = sqrt(PAI4) / univol

    Do it=1, ntyp
       ik = ista_k
       call new_radr_and_wos(ik,it)                 ! --> radr, wos

       Do il1=1, lpsmax(it)
          if ( il1 == iloc(it) ) cycle

          do tau1=1, itau(il1,it)
             do tau2=1, itau(il1,it)

                snl2 = 0.0d0;   nspher = 1

                call sphr( iend_kngp -ista_kngp +1, nspher, &
                     &     qx(ista_kngp:iend_kngp), &
                     &     qy(ista_kngp:iend_kngp), &
                     &     qz(ista_kngp:iend_kngp), ylm(ista_kngp:iend_kngp) )
                                                        ! -(bottom_Subr.)
                do ir = 1,nmesh(it)
                   facr = wos(ir) *psir_val(ir,il1,tau1,it) &
                        &         *psir_val(ir,il1,tau2,it)

                   Do ig=ista_kngp, iend_kngp
                      wka(ig) = vlength(ig) *radr(ir)
                   End do

                   call dsjnv( nspher -1, iend_kngp -ista_kngp +1, &
                        &      wka(ista_kngp:iend_kngp), &
                        &      wkb(ista_kngp:iend_kngp) )
                                                     ! -(bottom_Subr.)
                   Do ig=ista_kngp, iend_kngp
                      snl2(ig) = snl2(ig) + fac *facr *wkb(ig) *ylm(ig)
                   End do
                end do

                num = Table_for_ChargeDensityBasisFn( it, il1, tau1, tau2 )
                if ( num == 0 ) stop "PPPPPPP"

                Do ig=ista_kngp, iend_kngp
                   ChargeDensityBasisFn( ig, num, 1 ) = snl2(ig)
                End Do
             end do
          end do

       End Do
    End Do

    call dealloc_arrays
!
  contains

    subroutine alloc_arrays
      allocate( radr(mmesh) ); allocate( wos(mmesh) );
      allocate( qx(ista_kngp:iend_kngp) );
      allocate( qy(ista_kngp:iend_kngp) )
      allocate( qz(ista_kngp:iend_kngp) );
      allocate( ylm(ista_kngp:iend_kngp) ); allocate( vlength(ista_kngp:iend_kngp) )
      allocate( wka(ista_kngp:iend_kngp) );  allocate( wkb(ista_kngp:iend_kngp) )
      allocate( snl2(ista_kngp:iend_kngp) )
    end subroutine alloc_arrays

    subroutine dealloc_arrays
      deallocate( radr ); deallocate( wos );
      deallocate( qx ); deallocate( qy );  deallocate( qz ); deallocate( vlength )
      deallocate( wka ); deallocate( wkb ); deallocate( ylm )
      deallocate( snl2 )
    end subroutine dealloc_arrays

  end subroutine m_TFW_set_ChgDensityBasisFn_B

! --------------------------------------------------------
!
! Calc "averaged" nonlocal potential
!
! --------------------------------------------------------
  subroutine m_TFW_calc_nonlocal_pot_avg_A
    integer :: num, ig, it, ia, is
    integer :: ip, ilmt1, ilmt2
    real(kind=DP) :: weight, factor
    real(kind=DP) :: f1, f2, f3, ga, gb, gc, ph
    real(kind=DP), allocatable :: zfcos(:), zfsin(:)

    vnonlocal_avg = 0.0d0

    allocate( zfcos(ista_kngp:iend_kngp) ); allocate( zfsin(ista_kngp:iend_kngp) )
!
    num = 0
    Do it=1, ntyp
       Do ip=1, n_non0_lmtxlmt(it)
          ilmt1 = index_lmt1_lmt2(ip,it,1);  ilmt2 = index_lmt1_lmt2(ip,it,2)
          num = num +1

          if ( ilmt1 == ilmt2 ) then
             factor = 1.0d0
          else
             factor = 2.0d0
          endif

          Do ia=1, natm
             if ( ityp(ia) /= it ) cycle
             f1 = pos(ia,1)*PAI2; f2 = pos(ia,2)*PAI2; f3 = pos(ia,3)*PAI2
             Do ig=ista_kngp, iend_kngp
                ga = real(ngabc(ig,1),kind=DP)
                gb = real(ngabc(ig,2),kind=DP)
                gc = real(ngabc(ig,3),kind=DP)
                ph = ga *f1 +gb*f2 + gc*f3
                zfcos(ig) = cos(ph);  zfsin(ig) = sin(ph)
             End Do

             Do is=1, nspin
                if ( ipaw(it) == 0 ) then
                   weight = dion(ilmt1, ilmt2, it )
                else
                   weight = dion_paw(ilmt1, ilmt2, is, ia)
                endif
                if ( sw_hubbard == ON ) then
                   weight = weight + dhub( ilmt1,ilmt2,ia,is )
                endif

                weight = weight *hsr(ia, ilmt1,ilmt2,is )
!                weight = weight *hsro(ia, ilmt1,ilmt2,is )

                Do ig=ista_kngp, iend_kngp
                   vnonlocal_avg(ig,1,is) = vnonlocal_avg(ig,1,is) &
                        &      + factor *weight *( &
                        &                zfcos(ig) *ChargeDensityBasisFn( ig,num,1 ) &
                        &               -zfsin(ig) *ChargeDensityBasisFn( ig,num,2 ) )
                   if ( kimg == 2 ) then
                      vnonlocal_avg(ig,2,is) = vnonlocal_avg(ig,2,is) &
                           &   + factor *weight *( &
                           &             zfsin(ig) *ChargeDensityBasisFn( ig,num,1 ) &
                           &            +zfcos(ig) *ChargeDensityBasisFn( ig,num,2 ) )
                   endif
                Enddo

             End Do
          End Do
       End Do
    End Do

    deallocate( zfcos ); deallocate( zfsin )

  end subroutine m_TFW_calc_nonlocal_pot_avg_A

  subroutine m_TFW_calc_nonlocal_pot_avg_B
    integer :: num, ig, it, ia, is
    integer :: ip, ilmt1, ilmt2, il1, il2, tau1, tau2
    real(kind=DP) :: weight, factor
    real(kind=DP) :: f1, f2, f3, ga, gb, gc, ph
    real(kind=DP), allocatable :: zfcos(:), zfsin(:)

    vnonlocal_avg = 0.0d0

    allocate( zfcos(ista_kngp:iend_kngp) ); allocate( zfsin(ista_kngp:iend_kngp) )
!
    Do it=1, ntyp

       Do ia=1, natm
          if ( ityp(ia) /= it ) cycle
          f1 = pos(ia,1)*PAI2; f2 = pos(ia,2)*PAI2; f3 = pos(ia,3)*PAI2

          Do ig=ista_kngp, iend_kngp
             ga = real(ngabc(ig,1),kind=DP)
             gb = real(ngabc(ig,2),kind=DP)
             gc = real(ngabc(ig,3),kind=DP)
             ph = ga *f1 +gb*f2 + gc*f3
             zfcos(ig) = cos(ph);  zfsin(ig) = sin(ph)
          End Do

          Do ilmt1=1, ilmt(it)
             il1 = ltp(ilmt1,it);  tau1 = taup(ilmt1,it)
             Do ilmt2=1, ilmt(it)
                il2 = ltp(ilmt2,it);  tau2 = taup(ilmt2,it)

                if ( il1 /= il2 ) cycle

!                if ( ilmt1 == ilmt2 ) then
!                   factor = 1.0d0
!                else
!                   factor = 2.0d0
!                endif

                Do is=1, nspin
                   if ( ipaw(it) == 0 ) then
                      weight = dion(ilmt1, ilmt2, it )
                   else
                      weight = dion_paw(ilmt1, ilmt2, is, ia)
                   endif
                   if ( sw_hubbard == ON ) then
                      weight = weight + dhub( ilmt1,ilmt2,ia,is )
                   endif

                   weight = weight *hsr(ia, ilmt1,ilmt2,is )
!                   weight = weight *hsro(ia, ilmt1,ilmt2,is )

                   num = Table_for_ChargeDensityBasisFn( it, il1, tau1, tau2 )
                   Do ig=ista_kngp, iend_kngp
                      vnonlocal_avg(ig,1,is) = vnonlocal_avg(ig,1,is) &
                           &      + weight *( &
                           &             zfcos(ig) *ChargeDensityBasisFn( ig,num,1 ) )
                      if ( kimg == 2 ) then
                         vnonlocal_avg(ig,2,is) = vnonlocal_avg(ig,2,is) &
                              &   + weight *( &
                              &          zfsin(ig) *ChargeDensityBasisFn( ig,num,1 ) )
                      endif
                   Enddo
                End do
             End Do
          End Do
       End Do
    End Do

    deallocate( zfcos ); deallocate( zfsin )
  end subroutine m_TFW_calc_nonlocal_pot_avg_B

! ----------------------------
!
!  subroutines for Conjugate Gradient minimization
!
! ---------------------------
  subroutine m_TFW_alloc_deltaW
    allocate( deltaW(ista_kngp:iend_kngp,kimg,nspin) ); deltaW = 0.0d0
  end subroutine m_TFW_alloc_deltaW

  subroutine m_TFW_alloc_deltaV
    allocate( deltaV(ista_kngp:iend_kngp,kimg,nspin) ); deltaV = 0.0d0
  end subroutine m_TFW_alloc_deltaV

  subroutine m_TFW_alloc_Vnonlocal_avg
    allocate( vnonlocal_avg(ista_kngp:iend_kngp,kimg,nspin) )
    vnonlocal_avg = 0.0d0
  end subroutine m_TFW_alloc_Vnonlocal_avg

  subroutine m_TFW_dealloc_deltaW
    if (allocated(deltaW) ) deallocate(deltaW)
  end subroutine m_TFW_dealloc_deltaW

  subroutine m_TFW_dealloc_deltaV
    if (allocated(deltaV) ) deallocate(deltaV)
  end subroutine m_TFW_dealloc_deltaV

  subroutine m_TFW_dealloc_Vnonlocal_avg
    if (allocated(vnonlocal_avg) ) deallocate(vnonlocal_avg)
  end subroutine m_TFW_dealloc_Vnonlocal_avg

  subroutine m_TFW_alloc_Vin_and_Vlda
    allocate( Vin_tfw(ista_kngp:iend_kngp,kimg,nspin) );  Vin_tfw = 0.0d0
    allocate( Vlda_tfw(ista_kngp:iend_kngp,kimg,nspin) ); Vlda_tfw = 0.0d0
    allocate( Vout_tfw(ista_kngp:iend_kngp,kimg,nspin) ); Vout_tfw = 0.0d0
  end subroutine m_TFW_alloc_Vin_and_Vlda

  subroutine m_TFW_dealloc_Vin_and_Vlda
    if( allocated( Vin_tfw ) ) deallocate( Vin_tfw )
    if( allocated( Vlda_tfw ) ) deallocate( Vlda_tfw )
  end subroutine m_TFW_dealloc_Vin_and_Vlda

  subroutine m_TFW_CGoptimize_init
    Coeff_TF = ( 3.0*PAI*PAI )**(2.0d0/3.0d0) *3.0d0/10.0d0

    allocate( Psi(ista_kngp:iend_kngp,kimg,nspin) ); Psi = 0.0d0
    allocate( Phi(ista_kngp:iend_kngp,kimg,nspin) ); Phi = 0.0d0
    allocate( Phi_old(ista_kngp:iend_kngp,kimg,nspin) ); Phi_Old = 0.0d0
    allocate( Zeta(ista_kngp:iend_kngp,kimg,nspin) ); Zeta = 0.0d0

    allocate( Psi_Old(ista_kngp:iend_kngp,kimg,nspin) ); Psi_Old = 0.0d0
    allocate( Vkin_TF(ista_kngp:iend_kngp,kimg,nspin) ); Vkin_TF = 0.0d0

  end subroutine m_TFW_CGoptimize_init

  subroutine m_TFW_CGoptimize_finalize
    deallocate( Psi );   deallocate( Phi );   deallocate( Phi_old )
    deallocate( Zeta );

    deallocate( Psi_old ); deallocate( Vkin_TF )

  end subroutine m_TFW_CGoptimize_finalize

! ----------------------------
!
!  tools for using meshes
!
! ----------------------------
  subroutine map_matrix_on_FFTcd_mesh( matrix, afft )
    real(kind=DP), intent(in) :: matrix(ista_kngp:iend_kngp,kimg)
    real(kind=DP), intent(out) :: afft(ista_fftp:iend_fftp)

    integer :: is, j, i, ip, ilast
    real(kind=DP), allocatable :: afft_mpi1(:), afft_mpi2(:), afft_mpi3(:)

    ilast = min(iend_kngp,kgp_reduced)

    allocate( afft_mpi1(nfftp) )
    if ( npes > 1 ) then
       allocate( afft_mpi2(mp_fftp) );  allocate( afft_mpi3(mp_fftp) )
    end if

    afft = 0.0d0;   afft_mpi1=0.d0
    Do j=1, kimg
       do i=ista_kngp,ilast
          ip = ( igfp_l(i) -1 )*kimg +j
          afft_mpi1(ip) = matrix( i, j )
       end do
    End Do

    if ( npes > 1) then
       call mpi_barrier(MPI_CommGroup,ierr)
       do j=0,npes-1
          do i=nis_fftp(j),nie_fftp(j)
             afft_mpi2( i-nis_fftp(j) +1 ) = afft_mpi1(i)
          end do
          call mpi_allreduce( afft_mpi2, afft_mpi3, mp_fftp, &
               &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
          if ( j == mype ) then
             do i=ista_fftp,iend_fftp
                afft(i) = afft_mpi3(i-ista_fftp+1)
             end do
          end if
       end do
    else
       afft = afft_mpi1
    end if

    deallocate( afft_mpi1 );
    if ( npes > 1 ) then
       deallocate( afft_mpi2 ); deallocate( afft_mpi3 )
    endif

  end subroutine map_matrix_on_FFTcd_mesh

  subroutine extract_matrix_from_FFTcd_mesh( afft, matrix )
    real(kind=DP), intent(in) :: afft(ista_fftp:iend_fftp)
    real(kind=DP), intent(out) :: matrix(ista_kngp:iend_kngp,kimg)

    integer :: ik, i, ip
    real(kind=DP) :: rinplw
    real(kind=DP), allocatable :: afft_mpi1(:)

    matrix = 0.0d0

    allocate(afft_mpi1(nfftp))
    if ( npes >1 ) then
       call mpi_allgatherv( afft, nel_fftp(mype), mpi_double_precision, &  ! MPI
            &               afft_mpi1, nel_fftp, idisp_fftp, mpi_double_precision, &
            &               MPI_CommGroup, ierr )
    else
       afft_mpi1 = afft
    end if

    rinplw = 1.d0 /product(fft_box_size_CD(1:3,1))
    do ik=1, kimg
       do i=ista_kngp, iend_kngp          !for mpi
          ip = (igfp_l(i)-1)*kimg + ik
          matrix(i,ik) = afft_mpi1(ip)*rinplw
       end do
    end do

  end subroutine extract_matrix_from_FFTcd_mesh

  subroutine dot_product_of_vectors_Rspace( MatA, MatB, csum )
    real(kind=DP), intent(in) :: MatA(ista_kngp:iend_kngp,kimg,nspin)
    real(kind=DP), intent(in) :: MatB(ista_kngp:iend_kngp,kimg,nspin)
    real(kind=DP), intent(out) :: csum

    integer :: is, i
    real(kind=DP) :: csum_mpi
    real(kind=DP), allocatable :: afft(:), bfft(:)

    csum = 0.0d0

    call m_FFT_alloc_CD_box

    allocate( afft(ista_fftp:iend_fftp) )
    allocate( bfft(ista_fftp:iend_fftp) )

    Do is=1, nspin
       call map_matrix_on_FFTcd_mesh( MatA(:,:,is), afft )
       call m_FFT_CD_inverse_c( nfout, afft )       ! afft(G) -> afft(R)
       call map_matrix_on_FFTcd_mesh( MatB(:,:,is), bfft )
       call m_FFT_CD_inverse_c( nfout, bfft )       ! afft(G) -> afft(R)

       if ( kimg == 1 ) then
          Do i=ista_fftp, iend_fftp
             csum = csum +afft(i) *bfft(i)
          End do
       else
          Do i=ista_fftph, iend_fftph
             csum = csum +afft(2*i-1) *bfft(2*i-1) &
                  &      +afft(2*i) *bfft(2*i)
          End do
       endif
    End do

    if ( npes > 1 ) then
       call mpi_allreduce( csum, csum_mpi, 1, mpi_double_precision, &
            &              mpi_sum, MPI_CommGroup, ierr )
       csum = csum_mpi
    endif
    csum = csum *univol /product(fft_box_size_CD(1:3,1))

    call m_FFT_dealloc_CD_box
    deallocate( afft );  deallocate( bfft )

  end subroutine dot_product_of_vectors_Rspace

  subroutine dot_product_of_vectors_Gspace( MatA, MatB, csum )
    real(kind=DP), intent(in) :: MatA(ista_kngp:iend_kngp,kimg,nspin)
    real(kind=DP), intent(in) :: MatB(ista_kngp:iend_kngp,kimg,nspin)
    real(kind=DP), intent(out) :: csum

    integer :: is, i, ri
    real(kind=DP) :: csum_mpi

    csum = 0.0d0
    Do is=1, nspin
       Do ri=1, kimg
          Do i=ista_kngp, iend_kngp
             csum = csum + MatA(i,ri,is) *MatB(i,ri,is)
          End Do
       End do
    End do
    if ( npes > 1 ) then
       call mpi_allreduce( csum, csum_mpi, 1, mpi_double_precision, &
            &              mpi_sum, MPI_CommGroup, ierr )
       csum = csum_mpi
    endif

    csum = csum *univol

  end subroutine dot_product_of_vectors_Gspace

! ---------------------
!
!  Psi(r) is rho(r)^(1/2)
!
! --------------------
  subroutine m_TFW_init_Psi
    integer :: is, i, ri
    real(kind=DP) :: c1
    real(kind=DP), allocatable :: afft(:), bfft(:), cfft(:)

    allocate( afft(ista_fftp:iend_fftp) ); afft = 0.0d0
    allocate( cfft(ista_fftp:iend_fftp) ); cfft = 0.0d0

    call m_FFT_alloc_CD_box

    Do is=1, nspin
       afft = 0.0d0;  cfft = 0.0d0

! -- calc charge density on FFT-mesh (Rspace)
       call map_matrix_on_FFTcd_mesh( chgq_l(:,:,is), cfft )
       call m_FFT_CD_inverse_c( nfout, cfft )       ! afft(G) -> afft(R)

! -- calc "wavefunction" on FFT-mesh (Rspace)
       if ( kimg == 1 ) then
          Do i=ista_fftp, iend_fftp
             c1 = cfft(i)
             if ( c1 > 0.0 ) afft(i) = sqrt(c1)
          End do
       else
          Do i=ista_fftph, iend_fftph
             c1 = cfft(2*i-1)
             if ( c1 > 0.0 ) then
                afft(2*i-1) = sqrt(c1)
#if 0
             else
                afft(2*i) = sqrt(-c1)
#endif
             endif
          End do
       endif

! -- Back to G-space
       call m_FFT_CD_direct_c( nfout, afft )         ! afft(R) -> afft(G)
       call extract_matrix_from_FFTcd_mesh( afft, Psi(:,:,is) )
    End Do

    call m_FFT_dealloc_CD_box

    deallocate( afft );  deallocate( cfft )

#if 0
    call cut_high_G_component( Psi )
#endif

  end subroutine m_TFW_init_Psi

  subroutine m_TFW_store_Psi
    Psi_Old = Psi
  end subroutine m_TFW_store_Psi

! ----------------------------
!
!  deltaW term
!
! ----------------------------
  subroutine m_TFW_calc_deltaW
    integer :: is, i, ri
    real(kind=DP) :: c0, c1, c2, gx, gy, gz, g2
    real(kind=DP), allocatable :: afft(:), bfft(:), cfft(:)
    real(kind=DP), allocatable :: vlocal_tmp(:,:)

    allocate( afft(ista_fftp:iend_fftp) ); afft = 0.0d0
    allocate( bfft(ista_fftp:iend_fftp) ); bfft = 0.0d0
    allocate( cfft(ista_fftp:iend_fftp) ); cfft = 0.0d0
    allocate( vlocal_tmp(ista_kngp:iend_kngp,kimg) ); vlocal_tmp = 0.0d0

    Zeta = 0.0d0;  deltaW = 0.0d0

    call m_FFT_alloc_CD_box

    Do is=1, nspin
       afft = 0.0d0;  bfft = 0.0d0;  cfft = 0.0d0;  vlocal_tmp = 0.0d0

! -- calc "wavefunction" on FFT-mesh (Rspace)
       call map_matrix_on_FFTcd_mesh( Psi(:,:,is), afft )
       call m_FFT_CD_inverse_c( nfout, afft )       ! afft(G) -> afft(R)

! -- calc charge density on FFT-mesh (Rspace)
       if ( kimg == 1 ) then
          Do i=ista_fftp, iend_fftp
             cfft(i) = afft(i)**2
          End do
       else
          Do i=ista_fftph, iend_fftph
#if 0
             cfft(2*i-1) = afft(2*i-1)**2 +afft(2*i)**2
#else
             cfft(2*i-1) = afft(2*i-1)**2
#endif
          End do
       endif

! -- local +averaged nonlocal pot on FFT mesh (Rspace)
       vlocal_tmp(:,:) = vin_tfw(:,:,is)
       if ( use_averaged_nonlocal ) then
          vlocal_tmp(:,:) = vlocal_tmp(:,:) +vnonlocal_avg(:,:,is)
       endif

       call map_matrix_on_FFTcd_mesh( vlocal_tmp, bfft )
       call m_FFT_CD_inverse_c( nfout, bfft )       ! afft(G) -> afft(R)

! -- add Thomas Fermi term on FFT mesh (Rspace)
       c0 = 5.0d0 /3.0d0 *Coeff_TF *weight_TF_functional
       if ( kimg == 1 ) then
          Do i=ista_fftp, iend_fftp
             c1 = cfft(i) *nspin
             if ( c1 > 0.0 ) then
                c1 = c1**(2.0d0/3.0d0)
             else
                c1 = 0.0d0
             endif
             c2 = bfft(i)
             cfft(i) = ( c0 *c1 +c2 )*afft(i)
          End do

       else
          Do i=ista_fftph, iend_fftph
             c1 = cfft(2*i-1) *nspin
             if ( c1 > 0.0 ) then
                c1 = c1**(2.0d0/3.0d0)
             else
                c1 = 0.0d0
             endif
             c2 = bfft(2*i-1)

             cfft(2*i-1) = ( c0 *c1 +c2 )*afft(2*i-1)
#if 0
             cfft(2*i) = ( c0 *c1 +c2 )*afft(2*i)
#endif
          End do
       endif
! -- Back to G-space
       call m_FFT_CD_direct_c( nfout, cfft )         ! afft(R) -> afft(G)
       call extract_matrix_from_FFTcd_mesh( cfft, Zeta(:,:,is) )
                                                      ! zeta is used as temporay array

!       call m_FFT_CD_direct_c( nfout, afft )         ! afft(R) -> afft(G)
!       call extract_matrix_from_FFTcd_mesh( afft, Psi(:,:,is) )

! -- add von Weizsacker term
       Do ri=1, kimg
          Do i=ista_kngp, iend_kngp
             gx = ngabc(i,1)*rltv(1,1) +ngabc(i,2)*rltv(1,2) +ngabc(i,3)*rltv(1,3)
             gy = ngabc(i,1)*rltv(2,1) +ngabc(i,2)*rltv(2,2) +ngabc(i,3)*rltv(2,3)
             gz = ngabc(i,1)*rltv(3,1) +ngabc(i,2)*rltv(3,2) +ngabc(i,3)*rltv(3,3)
             g2 = ( gx**2 +gy**2 +gz**2 ) /2.0d0
             Zeta(i,ri,is) = Zeta(i,ri,is) +weight_Weiz_functional *g2 *Psi(i,ri,is)
          End do
       End Do

       Zeta(:,:,is) = Zeta(:,:,is) -efermi *Psi(:,:,is)
       deltaW(:,:,is) = -Zeta(:,:,is)
    End Do

    call m_FFT_dealloc_CD_box

    deallocate( afft );  deallocate( bfft );  deallocate( cfft )
    deallocate( vlocal_tmp )

#if 0
    call cut_high_G_component( deltaW )
#endif

  end subroutine m_TFW_calc_deltaW

! ------------------------
!
!  Zeta: steepest descent vector
!  Phi : Conjugate Gradient vector
!
! ------------------------
  subroutine m_TFW_calc_SD_direction( psi_MatH_psi, zeta_product )
    real(kind=DP), intent(out) :: psi_MatH_psi, zeta_product

    integer :: is, i, ri
    real(kind=DP) :: c0, c1, c2, gx, gy, gz, g2
    real(kind=DP), allocatable :: afft(:), bfft(:), cfft(:), dfft(:)
    real(kind=DP), allocatable :: vlocal_tmp(:,:), HPsi_tmp(:,:)

    allocate( afft(ista_fftp:iend_fftp) ); afft = 0.0d0
    allocate( bfft(ista_fftp:iend_fftp) ); bfft = 0.0d0
    allocate( cfft(ista_fftp:iend_fftp) ); cfft = 0.0d0
    allocate( dfft(ista_fftp:iend_fftp) ); dfft = 0.0d0
    allocate( vlocal_tmp(ista_kngp:iend_kngp,kimg) ); vlocal_tmp = 0.0d0

    psi_MatH_psi = 0.0d0
    Zeta = 0.0d0

    call m_FFT_alloc_CD_box

    Do is=1, nspin
       afft = 0.0d0;  bfft = 0.0d0;  cfft = 0.0d0;  dfft = 0.0d0

! -- calc "wavefunction" on FFT-mesh (Rspace)
       call map_matrix_on_FFTcd_mesh( Psi(:,:,is), afft )
       call m_FFT_CD_inverse_c( nfout, afft )       ! afft(G) -> afft(R)

! -- calc charge density on FFT-mesh (Rspace)
       if ( kimg == 1 ) then
          Do i=ista_fftp, iend_fftp
             cfft(i) = afft(i)**2
          End do
       else
          Do i=ista_fftph, iend_fftph
#if 0
             cfft(2*i-1) = afft(2*i-1)**2 +afft(2*i)**2
#else
             cfft(2*i-1) = afft(2*i-1)**2
#endif
          End do
       endif

! -- calc Thomas Fermi kinetic functional
       if ( kimg == 1 ) then
          Do i=ista_fftp, iend_fftp
             c1 = cfft(i) *nspin
             if ( c1 > 0.0 ) dfft(i) = Coeff_TF *c1**(2.0d0/3.0d0)
          End do
       else
          Do i=ista_fftph, iend_fftph
             c1 = cfft(2*i-1) *nspin
             if ( c1 > 0.0 ) dfft(2*i-1) = Coeff_TF *c1**(2.0d0/3.0d0)
          End do
       endif

! -- local +averaged nonlocal pot on FFT mesh (Rspace)
       if ( use_averaged_nonlocal ) then
          vlocal_tmp(:,:) = vlhxc_l(:,:,is) +vnonlocal_avg(:,:,is)
       else
          vlocal_tmp(:,:) = vlhxc_l(:,:,is)
       endif

       if ( use_deltaV ) vlocal_tmp(:,:) = vlocal_tmp(:,:) +deltaV(:,:,is)

       call map_matrix_on_FFTcd_mesh( vlocal_tmp, bfft )
       call m_FFT_CD_inverse_c( nfout, bfft )       ! afft(G) -> afft(R)

! -- add Thomas Fermi term on FFT mesh (Rspace)
       c1 = 5.0d0 /3.0d0 *weight_TF_functional
       if ( kimg == 1 ) then
          Do i=ista_fftp, iend_fftp
             cfft(i) = ( bfft(i) + c1*dfft(i) )*afft(i)
          End do
       else
          Do i=ista_fftph, iend_fftph
             cfft(2*i-1) = ( bfft(2*i-1) + c1*dfft(2*i-1) )*afft(2*i-1)
#if 0
             cfft(2*i) = ( bfft(2*i-1) + c1*dfft(2*i-1) )*afft(2*i)
#endif
          End do
       endif

! -- Back to G-space
       call m_FFT_CD_direct_c( nfout, cfft )         ! afft(R) -> afft(G)
       call extract_matrix_from_FFTcd_mesh( cfft, Zeta(:,:,is) )

!       call m_FFT_CD_direct_c( nfout, afft )         ! afft(R) -> afft(G)
!       call extract_matrix_from_FFTcd_mesh( afft, Psi(:,:,is) )

       call m_FFT_CD_direct_c( nfout, dfft )         ! afft(R) -> afft(G)
       call extract_matrix_from_FFTcd_mesh( dfft, Vkin_TF(:,:,is) )

! -- add von Weizsacker term
       Do ri=1, kimg
          Do i=ista_kngp, iend_kngp
             gx = ngabc(i,1)*rltv(1,1) +ngabc(i,2)*rltv(1,2) +ngabc(i,3)*rltv(1,3)
             gy = ngabc(i,1)*rltv(2,1) +ngabc(i,2)*rltv(2,2) +ngabc(i,3)*rltv(2,3)
             gz = ngabc(i,1)*rltv(3,1) +ngabc(i,2)*rltv(3,2) +ngabc(i,3)*rltv(3,3)
             g2 = ( gx**2 +gy**2 +gz**2 ) /2.0d0
             Zeta(i,ri,is) = Zeta(i,ri,is) +Weight_Weiz_functional *g2 *Psi(i,ri,is)
          End do
       End Do
       if ( use_deltaW ) Zeta(:,:,is) = Zeta(:,:,is) + deltaW(:,:,is)
    End do

    call m_FFT_dealloc_CD_box
    deallocate( afft ); deallocate( bfft );  deallocate( cfft );  deallocate( dfft )
    deallocate( vlocal_tmp )

#if 0
    call cut_high_G_component( Zeta )
#endif

! ------------------
! if normalization is required
    call dot_product_of_vectors_Gspace( Psi, Psi, c1 )
    if ( ipritfwfunc >= 3 ) then
       if ( mype == 0 ) then
          write(nfout,*) "**** TFW functional info ***"
          write(nfout,*) "dot product of Psi-Psi before normalization: ", c1
       endif
    endif

    c1 = sqrt( totch /c1 )
    Psi = Psi *c1
! ------------------

    call dot_product_of_vectors_Gspace( Psi, Zeta, psi_MatH_psi )

    Zeta = psi_MatH_psi/totch *Psi -Zeta

    call dot_product_of_vectors_Gspace( zeta, zeta, zeta_product )
    if ( ipritfwfunc >= 3 ) then
       if ( mype == 0 ) write(nfout,*) "zeta_product is ", zeta_product
    endif
  end subroutine m_TFW_calc_SD_direction

  subroutine m_TFW_calc_CG_direction( FirstFlag, zeta_product, zeta_product_old, &
       &                              phi_MatH_psi )
    logical, intent(in) ::FirstFlag
    real(kind=DP), intent(in) :: zeta_product, zeta_product_old
    real(kind=DP), intent(out) :: phi_MatH_psi

    real(kind=DP) :: gamma, c1, c2

    if ( FirstFlag ) then
       gamma = 0.0d0
       Phi = zeta
    else
       gamma = zeta_product / zeta_product_old
       Phi = zeta +gamma *Phi_old
    endif
    Phi_old = Phi

#if 0
    call dot_product_of_vectors_Gspace( Phi, Phi, c1 )
    if ( ipritfwfunc >= 3 ) then
       if ( mype == 0 ) then
          write(nfout,*) "dot product of Phi-Phi before normalization: ", c1
       endif
    endif

    call dot_product_of_vectors_Gspace( Psi, Phi, c1 )
    Phi = Phi - c1 / totch *Psi
#endif

    call dot_product_of_vectors_Gspace( Phi, Phi, c1 )
    if ( ipritfwfunc >= 3 ) then
       if ( mype == 0 ) then
          write(nfout,*) "dot product of Phi-Phi after orthogonalization: ", c1
       endif
    endif

#if 1
    if ( abs(c1) < 1.0E-10 ) then
       if ( mype == 0) write(nfout,*) 'Norm of Phi is too small, ', abs(c1)
       Phi = Phi_old
    else
       c2 = totch /c1
       Phi = Phi *sqrt(c2)
    endif
#endif

    call dot_product_of_vectors_Gspace( Phi, Zeta, c1 )
    phi_MatH_psi = -c1

  end subroutine m_TFW_calc_CG_direction

  subroutine m_TFW_calc_phi_MatH_phi( phi_MatH_phi )
    real(kind=DP), intent(out) :: phi_MatH_phi

    integer :: is, i, ri
    real(kind=DP) :: c0, c1, c2, gx, gy, gz, g2
    real(kind=DP), allocatable :: afft(:), bfft(:), cfft(:)
    real(kind=DP), allocatable :: vlocal_tmp(:,:), HPhi_tmp(:,:)

    allocate( afft(ista_fftp:iend_fftp) ); afft = 0.0d0
    allocate( bfft(ista_fftp:iend_fftp) ); bfft = 0.0d0
    allocate( cfft(ista_fftp:iend_fftp) ); cfft = 0.0d0
    allocate( vlocal_tmp(ista_kngp:iend_kngp,kimg) ); vlocal_tmp = 0.d0

    phi_MatH_phi = 0.0d0

    call m_FFT_alloc_CD_box

    Do is=1, nspin
       afft = 0.0d0;  bfft = 0.0d0;  cfft = 0.0d0

! -- calc "wavefunction" on FFT-mesh (Rspace)
       call map_matrix_on_FFTcd_mesh( Phi(:,:,is), afft )
       call m_FFT_CD_inverse_c( nfout, afft )       ! afft(G) -> afft(R)

! -- calc charge density on FFT-mesh (Rspace)
       call map_matrix_on_FFTcd_mesh( chgq_l(:,:,is), cfft )
       call m_FFT_CD_inverse_c( nfout, cfft )       ! afft(G) -> afft(R)

! -- local +averaged nonlocal pot on FFT mesh (Rspace)
       if ( use_averaged_nonlocal ) then
          vlocal_tmp(:,:) = vlhxc_l(:,:,is) +vnonlocal_avg(:,:,is)
       else
          vlocal_tmp(:,:) = vlhxc_l(:,:,is)
       endif

       if ( use_deltaV ) vlocal_tmp(:,:) = vlocal_tmp(:,:) +deltaV(:,:,is)

       call map_matrix_on_FFTcd_mesh( vlocal_tmp, bfft )
       call m_FFT_CD_inverse_c( nfout, bfft )       ! afft(G) -> afft(R)

! -- add Thomas Fermi term on FFT mesh (Rspace)
       if ( kimg == 1 ) then
          Do i=ista_fftp, iend_fftp
             c0 = cfft(i) *nspin
             c1 = 5.0d0 /3.0d0 *Coeff_TF *weight_TF_functional *c0**(2.0d0/3.0d0)
             c2 = bfft(i)
             cfft(i) = ( c1 +c2 )*afft(i)
          End do
       else
          Do i=ista_fftph, iend_fftph
             c0 = cfft(2*i-1) *nspin
             c1 = 5.0d0 /3.0d0 *Coeff_TF *weight_TF_functional *c0**(2.0d0/3.0d0)
             c2 = bfft(2*i-1)
             cfft(2*i-1) = ( c1 +c2 )*afft(2*i-1)
          End do
       endif
! -- Back to G-space
       call m_FFT_CD_direct_c( nfout, cfft )         ! afft(R) -> afft(G)
       call extract_matrix_from_FFTcd_mesh( cfft, Zeta(:,:,is) )

! -- add von Weizsacker term
       Do ri=1, kimg
          Do i=ista_kngp, iend_kngp
             gx = ngabc(i,1)*rltv(1,1) +ngabc(i,2)*rltv(1,2) +ngabc(i,3)*rltv(1,3)
             gy = ngabc(i,1)*rltv(2,1) +ngabc(i,2)*rltv(2,2) +ngabc(i,3)*rltv(2,3)
             gz = ngabc(i,1)*rltv(3,1) +ngabc(i,2)*rltv(3,2) +ngabc(i,3)*rltv(3,3)
             g2 = ( gx**2 +gy**2 +gz**2 ) /2.0d0
             Zeta(i,ri,is) = Zeta(i,ri,is) +Weight_Weiz_functional *g2 *Phi(i,ri,is)
          End do
       End Do

       if ( use_deltaW ) Zeta(:,:,is) = Zeta(:,:,is) + deltaW(:,:,is)
    End do

    call m_FFT_dealloc_CD_box

    deallocate( afft );  deallocate( bfft );  deallocate( cfft )
    deallocate( vlocal_tmp )

    call dot_product_of_vectors_Gspace( Phi, Zeta, phi_MatH_phi )

  end subroutine m_TFW_calc_phi_MatH_phi

  subroutine m_TFW_try_new_Psi_and_charge( dx )
    real(kind=DP), intent(in) :: dx

    integer :: is, i
    real(kind=DP) :: csum, c1, cos_th, sin_th
    real(kind=DP), allocatable :: afft(:), bfft(:), cfft(:), dfft(:)
    real(kind=DP), allocatable :: K_phi(:,:,:)

    if ( use_preconditioning ) then
       allocate( K_phi( ista_kngp:iend_kngp, kimg, nspin ) )

       call m_TFW_set_preconditioned_Phi( Phi, K_Phi )
       Psi = Psi_Old + dx *K_Phi

       deallocate( K_Phi )
    else
       Psi = Psi_Old + dx *Phi
    endif

    if ( ipritfwfunc >= 3 ) then
       call dot_product_of_vectors_Gspace( Psi_Old, Psi_Old, csum )
       if ( mype == 0 ) write(nfout,*) "dot product of PsiOld-PsiOld : ", csum
    endif

    if ( ipritfwfunc >= 3 ) then
       call dot_product_of_vectors_Gspace( Phi, Phi, csum )
       if ( mype == 0 ) write(nfout,*) "dot product of Phi-Phi : ", csum
    endif

! ------------------
! if normalization is required
    call dot_product_of_vectors_Gspace( Psi, Psi, c1 )
    if ( ipritfwfunc >= 3 ) then
       if ( mype == 0 ) then
          write(nfout,*) "dot product of Psi-Psi before normalization: ", c1
          write(nfout,*) "*** ***"
       endif
    endif

    c1 = sqrt( totch /c1 )
    Psi = Psi *c1
! ------------------

    allocate( afft(ista_fftp:iend_fftp) ); afft = 0.0d0
    allocate( cfft(ista_fftp:iend_fftp) ); cfft = 0.0d0
    allocate( dfft(ista_fftp:iend_fftp) ); dfft = 0.0d0

    call m_FFT_alloc_CD_box

    csum = 0.0d0
    Do is=1, nspin
       afft = 0.0d0;  cfft = 0.0d0;   dfft = 0.0d0

! -- calc "wavefunction" on FFT-mesh (Rspace)
       call map_matrix_on_FFTcd_mesh( Psi(:,:,is), cfft )
       call m_FFT_CD_inverse_c( nfout, cfft )       ! afft(G) -> afft(R)

! -- calc charge density on FFT-mesh (Rspace)
       if ( kimg == 1 ) then
          Do i=ista_fftp, iend_fftp
             afft(i) = cfft(i)**2
          End do
       else
          Do i=ista_fftph, iend_fftph
#if 0
             afft(2*i-1) = cfft(2*i-1)**2 +cfft(2*i)**2
#else
             afft(2*i-1) = cfft(2*i-1)**2
#endif
          End do
       endif

! -- calc Thomas Fermi functional
       if ( kimg == 1 ) then
          Do i=ista_fftp, iend_fftp
             c1 = afft(i) *nspin
             dfft(i) = Coeff_TF *c1**(2.0d0/3.0d0)
          End do
       else
          Do i=ista_fftph, iend_fftph
             c1 = afft(2*i-1) *nspin
             dfft(2*i-1) = Coeff_TF *c1**(2.0d0/3.0d0)
          End do
       endif

! Back to G-space
       call m_FFT_CD_direct_c( nfout, afft )         ! afft(R) -> afft(G)
       call extract_matrix_from_FFTcd_mesh( afft, chgq_l(:,:,is) )

       call m_FFT_CD_direct_c( nfout, dfft )         ! afft(R) -> afft(G)
       call extract_matrix_from_FFTcd_mesh( dfft, Vkin_TF(:,:,is) )
    End do

    deallocate(afft); deallocate(cfft);  deallocate(dfft)
    call m_FFT_dealloc_CD_box

#if 0
    call m_TFW_remove_negative_charge
#endif

  end subroutine m_TFW_try_new_Psi_and_charge

  subroutine m_TFW_calc_energy_terms( ekin_TF, ekin_Weiz, &
       &                              ene_nonlocal_avg, ene_deltaV, ene_deltaW )
    real(kind=DP), intent(out) :: ekin_TF, ekin_Weiz
    real(kind=DP), intent(out) :: ene_nonlocal_avg, ene_deltaV, ene_deltaW

    integer :: is, i, ri
    real(kind=DP) :: csum, csum_mpi, gx, gy, gz, g2

    ene_deltaW = 0.0d0;  ene_deltaV = 0.0d0
    ene_nonlocal_avg = 0.0d0;  ekin_TF = 0.0d0;  ekin_Weiz = 0.0d0

! -- Thomas Fermi kinetic energy
    call dot_product_of_vectors_Gspace( chgq_l, Vkin_TF, csum )
!!!    ekin_TF = weight_TF_functional *csum /dble(nspin)
    ekin_TF = weight_TF_functional *csum

! -- von Weizsacker kinetic energy
    csum = 0.0d0
    Do is=1, nspin
       Do ri=1, kimg
          Do i=ista_kngp, iend_kngp
             gx = ngabc(i,1)*rltv(1,1) +ngabc(i,2)*rltv(1,2) +ngabc(i,3)*rltv(1,3)
             gy = ngabc(i,1)*rltv(2,1) +ngabc(i,2)*rltv(2,2) +ngabc(i,3)*rltv(2,3)
             gz = ngabc(i,1)*rltv(3,1) +ngabc(i,2)*rltv(3,2) +ngabc(i,3)*rltv(3,3)
             g2 = ( gx**2 +gy**2 +gz**2 ) /2.0d0
             csum = csum + g2 *Psi(i,ri,is)**2
          End do
       End Do
    End do
    if ( npes > 1 ) then
       call mpi_allreduce( csum, csum_mpi, 1, mpi_double_precision, mpi_sum, &
            &              MPI_CommGroup, ierr )
       csum = csum_mpi
    endif
    ekin_Weiz = Weight_Weiz_functional *csum *univol

! nonlocal avg energy
    if ( use_averaged_nonlocal ) then
       call dot_product_of_vectors_Gspace( chgq_l, vnonlocal_avg, csum )
       ene_nonlocal_avg = csum
    endif

! deltaV eergy
    if ( use_deltaV ) then
       call dot_product_of_vectors_Gspace( chgq_l, deltaV, csum )
       ene_deltaV = -csum
    endif

! correction energy
    if ( use_deltaW ) then
       call dot_product_of_vectors_Gspace( deltaW, Psi, csum )
       ene_deltaW = csum *2.0d0
    endif

  end subroutine m_TFW_calc_energy_terms

! ---------------------
!
! others
!
! ---------------------
  subroutine cut_high_G_component( mat )
    real(kind=DP), intent(inout) :: mat( ista_kngp:iend_kngp, kimg, nspin )

    integer :: i

    Do i=ista_kngp, iend_kngp
       if ( i> kg_tfw ) mat(i,:,:) = 0.0d0
    End do

  end subroutine cut_high_G_component

  subroutine m_TFW_remove_negative_charge
    integer :: is, i
    real(kind=DP) :: c1
    real(kind=DP), allocatable :: cfft(:)

    allocate( cfft(ista_fftp:iend_fftp) )

    call m_FFT_alloc_CD_box

    Do is=1, nspin
       cfft = 0.0d0

! -- calc charge density on FFT-mesh (Rspace)
       call map_matrix_on_FFTcd_mesh( chgq_l(:,:,is), cfft )
       call m_FFT_CD_inverse_c( nfout, cfft )       ! afft(G) -> afft(R)

! -- calc "wavefunction" on FFT-mesh (Rspace)
       if ( kimg == 1 ) then
          Do i=ista_fftph, iend_fftph
             c1 = cfft(2*i-1)
             if ( c1 < 0.0 ) cfft(2*i-1) = 0.0d0
          End do
       else
          Do i=ista_fftp, iend_fftp
             c1 = cfft(i)
             if ( c1 < 0.0 ) cfft(i) = 0.0d0
          End do
       endif
! -- Back to G-space
       call m_FFT_CD_direct_c( nfout, cfft )         ! afft(R) -> afft(G)
       call extract_matrix_from_FFTcd_mesh( cfft, chgq_l(:,:,is) )
    End do

    call m_FFT_dealloc_CD_box
    deallocate( cfft )

    if ( mype == 0 ) then
       c1 = 0.0d0
       Do is=1, nspin
          c1 = c1 + chgq_l(1,1,is)
       End do
       c1 = totch /(c1*univol)
    endif
    if ( npes > 1 ) then
       call mpi_bcast( c1, 1, mpi_double_precision, 0, MPI_CommGroup, ierr )
    endif
    chgq_l = c1 *chgq_l
  end subroutine m_TFW_remove_negative_charge

! -------------------------------
!
! Sometimes, the variable chgq_l/vlhxc_l obtained by TFW equations
! are not well screened. It means, the lower G components have non-negligible
! amplitudes.
!
! Using such variables as new inputs for Kohn-Sham equation,
! the calculations become unstable due to the charge sloshing.
! In order to prevent this, Kerker mixing is adopted.
!
! ------------------------------
  subroutine m_TFW_Kerker_mixing( V_in, V_out, rmxt )
    real(kind=DP), intent(in) :: rmxt
    real(kind=DP), intent(in)    :: V_in(  ista_kngp:iend_kngp,kimg,nspin )
    real(kind=DP), intent(inout) :: V_out( ista_kngp:iend_kngp,kimg,nspin )
!
    integer :: i, is, ri
    real(kind=DP) :: gg, q0
    real(kind=DP) :: rmxtrc(nspin)
!
    allocate( c_p_pot( ista_kngp:iend_kngp,nspin) );  c_p_pot = 0.0d0

    rmxtrc = rmxt

#if 0
    call precon_4_pot_mix( rmxtrc, c_p_pot )
#else
    q0 = 10.0d0
!    q0 = 1.0d0

    Do i=ista_kngp, iend_kngp
       gg = gr_l(i)*gr_l(i)
       Do is=1, nspin
!          c_p_pot(i,is) = rmxtrc(is) *max(gg/(gg+q0),amin)
!          c_p_pot(i,is) = rmxtrc(is) *gg/(gg+q0)
          c_p_pot(i,is) = gg/(gg+q0)
!          c_p_pot(i,is) = rmxtrc(is)
       End do
    End Do
#endif

    Do is=1, nspin
       Do ri=1, kimg
          Do i=ista_kngp, iend_kngp
             V_out(i,ri,is) = V_in(i,ri,is) &
                  &         + c_p_pot(i,is)*( V_out(i,ri,is) -V_in(i,ri,is) )
          End do
       End do
    End Do

    deallocate( c_p_pot )

  end subroutine m_TFW_Kerker_mixing

  subroutine m_TFW_set_preconditioned_Phi( V_in, V_out )
    real(kind=DP), intent(in)    :: V_in(  ista_kngp:iend_kngp,kimg,nspin )
    real(kind=DP), intent(inout) :: V_out( ista_kngp:iend_kngp,kimg,nspin )

    integer :: i, is, ri
    real(kind=DP) :: gg, q0
    real(kind=DP) :: rmxtrc(nspin)
!
    allocate( c_p_pot( ista_kngp:iend_kngp,nspin) );  c_p_pot = 0.0d0

#if 0
    call precon_4_pot_mix( rmxtrc, c_p_pot )
#else
    q0 = 10.0d0
!    q0 = 1.0d0

    Do i=ista_kngp, iend_kngp
       gg = gr_l(i)*gr_l(i)
       Do is=1, nspin
          c_p_pot(i,is) = gg/(gg+q0)
       End do
    End Do
#endif

    Do is=1, nspin
       Do ri=1, kimg
          Do i=ista_kngp, iend_kngp
             V_out(i,ri,is) = c_p_pot(i,is)*V_in(i,ri,is)
          End do
       End do
    End Do
    deallocate( c_p_pot )

  end subroutine m_TFW_set_preconditioned_Phi

end module m_ThomasFermiW_Potential
