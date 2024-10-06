#define MOMENT_AS_PSEUDO_VECTOR

module m_OP_Moment
! $Id: m_OP_Moment.F90 479 2016-03-12 12:30:51Z ktagami $

! =========== Contributions ===================================
!
! Through the courtesy of contributors, the following functions are added.
!
! Company:  ASMS Co.,Ltd.
!       Calculation of Orbital Moment
!
! Note:
!       In order to work the whole functions, LIB_ASMS_INTEGRAL is required.
!
! =============================================================

  use m_Control_Parameters,  only : proj_attribute, proj_group, num_proj_elems, &
       &                            max_projs, ndim_spinor, ndim_magmom, ndim_chgpot, &
       &                            iprimagmom

  use m_Const_Parameters,     only : DP, CMPLDP, zi
  use m_Ionic_System,         only : natm, iproj_group

  use m_SpinOrbit_Potential,  only :  MatU_ylm_RC_L0,  MatU_ylm_RC_L1, &
       &                              MatU_ylm_RC_L2,  MatU_ylm_RC_L3
  use m_Orbital_Population,   only :  om, om_aimag, ommix, ommix_aimag, i2lp
  use m_Parallelization,       only : mype

  use m_ES_NonCollinear,       only :  m_ES_MagMom_To_DensMat_porb, &
       &                               m_ES_set_Pauli_Matrix
  use m_Files,                 only : nfout

! ==== KT_add === 2014/08/26
  use m_KPoints,           only : kv3, vkxyz
  use m_Const_Parameters,   only : PAI4, CARTS, ON, BOHR
  use m_Control_Parameters,  only : kimg, noncol, neg, sw_mix_charge_hardpart, &
       &                            sw_use_add_proj, sw_use_contracted_psir
  use m_PseudoPotential,    only : ltp, mtp, taup, nlmt, isph, il2p, dl2p, ilmt, &
       &                            q, q_noncl, nlmta, lmta, nmesh, mmesh, radr, wos, &
       &                            flg_paw, ia2ia_symmtry_op_inv, &
       &                            crotylm_paw, nylm_paw, iylm_paw, wf_mnrc, &
       &                            psir_val, &
       &                            nloc, ntau, lpsmax, itau, &
       &                            ilmt_add, nlmt_add, ltp_add, mtp_add, &
       &                            lmta_add, nlmta_add
  use m_PlaneWaveBasisSet,    only : ngabc, kgp, kg1, iba, nbase
  use m_Parallelization,    only : ista_kngp, iend_kngp, npes, MPI_CommGroup, ierr, &
       &                           map_k, map_e, map_z, myrank_k, myrank_e, np_e, &
       &                           ista_k, iend_k
  use m_NonLocal_Potential, only : new_radr_and_wos
  use m_Ionic_System,        only : ityp, pos, cps, iatomn, ntyp, speciesname
  use m_Crystal_Structure,   only : rltv, univol, nopr, op, altv
  use m_CS_Magnetic,     only : magmom_dir_inversion_opr_flag, determinant_op
  use m_Charge_Density,      only : hsr, hsi, hsr_add, hsi_add
  use m_CD_Mag_Moment,        only : rad_cov, rad_cov_default, product_psir_val, &
       &                             RhoMag_on_atom
  use m_ES_NonCollinear,     only : m_ES_MagMom_to_DensMat_Gspace, &
       &                            m_ES_MagMom_To_DensMat_hsr, &
       &                            m_ES_DensMat_To_MagMom_hsr
  use m_Electronic_Structure,  only : zaj_l, occup_l
! =============== 2014/08/26

  use m_Control_Parameters,   only : charge_symm_mode
  use m_Kpoints,             only : nopr_from_fbz_to_ibz, flg_opr_from_fbz_to_ibz
  use m_Const_Parameters,    only : chg_symm_level1
  use mpi

  implicit none
!  include 'mpif.h'

  complex(kind=CMPLDP) :: Mat_L_with_cmplx_ylm_L0(  0:0, 0:0, 3 )
  complex(kind=CMPLDP) :: Mat_L_with_cmplx_ylm_L1( -1:1,-1:1, 3 )
  complex(kind=CMPLDP) :: Mat_L_with_cmplx_ylm_L2( -2:2,-2:2, 3 )
  complex(kind=CMPLDP) :: Mat_L_with_cmplx_ylm_L3( -3:3,-3:3, 3 )
!
  complex(kind=CMPLDP) :: Mat_L_with_real_ylm_L0( 1, 1, 3 )
  complex(kind=CMPLDP) :: Mat_L_with_real_ylm_L1( 3, 3, 3 )
  complex(kind=CMPLDP) :: Mat_L_with_real_ylm_L2( 5, 5, 3 )
  complex(kind=CMPLDP) :: Mat_L_with_real_ylm_L3( 7, 7, 3 )

! ==== KT_add === 2014/08/26
  real(kind=DP), allocatable :: OrbMag_on_atom(:,:)
  real(kind=DP), allocatable :: rho_ylm1_ylm2_r(:,:,:,:), rho_ylm1_ylm2_i(:,:,:,:)
  complex(kind=CMPLDP), allocatable :: wfn_sph_decomposed(:,:,:)
! =============== 2014/08/26

contains

  subroutine m_OP_calc_orbmom_from_OCC           ! From occupation marix
    integer :: ia, it, i, ig, ie, ip, l, m1, m2
    integer :: is1, is2, istmp, size1
!
    real(kind=DP) :: orbmom(3), spinmom(3)
    complex(kind=CMPLDP) :: ztmp(3)
    complex(kind=CMPLDP) :: PauliMatrix( ndim_magmom, ndim_spinor, ndim_spinor )
!
    real(kind=DP), allocatable :: dmmat_r_magmom( :,:,: )
    real(kind=DP), allocatable :: dmmat_i_magmom( :,:,: )
    complex(kind=CMPLDP), allocatable :: dmmat_ssrep( :,:,: )

    logical :: FirstFlag = .true.
    logical :: print_header

! -----------
    if ( iprimagmom < 3 ) return

    if ( FirstFlag ) then
       call m_OP_calc_MatL_orb_s_to_f
       FirstFlag = .false.
    endif
    print_header = .true.

    call m_ES_set_Pauli_Matrix( PauliMatrix )

    do ia=1, natm
       ig = iproj_group(ia)
       if ( ig <1 ) cycle

       do i=1,num_proj_elems(ig)
          ip = proj_group(i,ig)
          it = proj_attribute(ip)%ityp
          ie = proj_attribute(ip)%ielem
          l  = proj_attribute(ip)%l

          size1 = i2lp(ip)
          allocate( dmmat_r_magmom(size1,size1,ndim_magmom) ); dmmat_r_magmom = 0.0d0
          allocate( dmmat_i_magmom(size1,size1,ndim_magmom) ); dmmat_i_magmom = 0.0d0
          allocate( dmmat_ssrep(size1,size1,ndim_chgpot) );  dmmat_ssrep = 0.0d0

!          dmmat_r_magmom(:,:,:) = om(:,:,ie,ia,:)
!          dmmat_i_magmom(:,:,:) = om_aimag(:,:,ie,ia,:)
          dmmat_r_magmom(:,:,:) = ommix(:,:,ie,ia,:)
          dmmat_i_magmom(:,:,:) = ommix_aimag(:,:,ie,ia,:)

          call m_ES_MagMom_To_DensMat_porb( size1**2, dmmat_r_magmom, dmmat_i_magmom, &
               &                            dmmat_ssrep )

          orbmom = 0.0d0; spinmom = 0.0d0

          Do is1=1, ndim_spinor
             Do is2=1, ndim_spinor
                if ( is1 /= is2 ) cycle

                istmp = ( is1 -1 )*ndim_spinor + is2

                Do m1=1, size1
                   Do m2=1, size1
                      ztmp = 0.0d0
                      if ( l == 1 ) then
                         ztmp(:) = Mat_L_with_real_ylm_L1( m2, m1, : )
                      else if ( l == 2 ) then
                         ztmp(:) = Mat_L_with_real_ylm_L2( m2, m1, : )
                      else if ( l == 3 ) then
                         ztmp(:) = Mat_L_with_real_ylm_L3( m2, m1, : )
                      endif

                      orbmom(1) = orbmom(1) +dmmat_ssrep(m1,m2,istmp) *ztmp(1)
                      orbmom(2) = orbmom(2) +dmmat_ssrep(m1,m2,istmp) *ztmp(2)
                      orbmom(3) = orbmom(3) +dmmat_ssrep(m1,m2,istmp) *ztmp(3)
                   End do
                End do
             End do
          End do

          Do is1=1, ndim_spinor
             Do is2=1, ndim_spinor
                istmp = ( is1 -1 )*ndim_spinor + is2
                Do m1=1, size1
                   Do m2=1, size1
                      if ( m1 /= m2 ) cycle

                      spinmom(1) = spinmom(1) +dmmat_ssrep(m1,m2,istmp) &
                           &                  *PauliMatrix(2,is2,is1)
                      spinmom(2) = spinmom(2) +dmmat_ssrep(m1,m2,istmp) &
                           &                  *PauliMatrix(3,is2,is1)
                      spinmom(3) = spinmom(3) +dmmat_ssrep(m1,m2,istmp) &
                           &                  *PauliMatrix(4,is2,is1)
                   End do
                End do
             End do
          End do
!
          deallocate( dmmat_r_magmom ); deallocate( dmmat_i_magmom )
          deallocate( dmmat_ssrep )
! --
          call  print_spin_orb_mom_local( ia, l, spinmom, orbmom, print_header )
          print_header = .false.
       end do
    End do

  contains

    subroutine print_spin_orb_mom_local( ia, l, spinmom, orbmom, print_header )
      integer, intent(in) :: ia, l
      real(kind=DP), intent(in) :: spinmom(3), orbmom(3)
      logical, intent(in) :: print_header

      if ( mype /= 0 ) return

      if (print_header) then
         write(nfout,*) '! ------------ Local Momemnt ( occmat )  --- '
      endif

      write(nfout,*) '!   ia   l   type          mx             my             mz  '
!
      write(nfout,'(I7,I4,3X,A5,3F15.8)') ia, l, &
           &                              "spin ", spinmom(1), spinmom(2), spinmom(3)
      write(nfout,'(14X,     A5,3F15.8)') &
           &                              "orb  ", orbmom(1), orbmom(2), orbmom(3)
!
      write(nfout,*) '! ------ '

    end subroutine print_spin_orb_mom_local

  end subroutine m_OP_calc_orbmom_from_OCC


! ------------------
  subroutine m_OP_calc_MatL_orb_s_to_f
    real(kind=DP) :: theta, phi

    theta = 0.0; phi = 0.0d0
    call m_OP_calc_MatL_with_Cmplx_ylm( 0, theta, phi, &
         &                              Mat_L_with_cmplx_ylm_L0 )
    call m_OP_calc_MatL_with_Real_ylm( 0, MatU_ylm_RC_L0, &
         &                             Mat_L_with_cmplx_ylm_L0, &
         &                             Mat_L_with_real_ylm_L0 )
!
    call m_OP_calc_MatL_with_Cmplx_ylm( 1, theta, phi, &
         &                              Mat_L_with_cmplx_ylm_L1 )
    call m_OP_calc_MatL_with_Real_ylm( 1, MatU_ylm_RC_L1, &
         &                             Mat_L_with_cmplx_ylm_L1, &
         &                             Mat_L_with_real_ylm_L1 )
!
    call m_OP_calc_MatL_with_Cmplx_ylm( 2, theta, phi, &
         &                              Mat_L_with_cmplx_ylm_L2 )
    call m_OP_calc_MatL_with_Real_ylm( 2, MatU_ylm_RC_L2, &
         &                              Mat_L_with_cmplx_ylm_L2, &
         &                              Mat_L_with_real_ylm_L2 )
!
    call m_OP_calc_MatL_with_Cmplx_ylm( 3, theta, phi, &
         &                              Mat_L_with_cmplx_ylm_L3 )
    call m_OP_calc_MatL_with_Real_ylm( 3, MatU_ylm_RC_L3, &
         &                             Mat_L_with_cmplx_ylm_L3, &
         &                             Mat_L_with_real_ylm_L3 )
  end subroutine m_OP_calc_MatL_orb_s_to_f


  subroutine m_OP_calc_MatL_with_Real_ylm( l_in, MatU_ylm_RC, &
       &                                   Mat_L_with_cmplx_ylm, &
       &                                   Mat_L_with_real_ylm )
    integer, intent(in) :: l_in
    Complex(kind=CMPLDP), intent(in) :: MatU_ylm_RC( 2*l_in+1, -l_in:l_in )
    complex(kind=CMPLDP), intent(in) :: &
         &    Mat_L_with_cmplx_ylm( -l_in:l_in, -l_in:l_in, 3 )
    complex(kind=CMPLDP), intent(out) :: &
         &    Mat_L_with_real_ylm( 2*l_in+1, 2*l_in+1, 3 )

    integer :: m_min, m_max
    !
    integer :: ixyz
    integer :: m1, m2, ma, mb
!
    Complex(kind=CMPLDP) :: ztmp( 3 )
!
    m_min = -l_in;  m_max = l_in
!
    Do m1=1, 2 *l_in +1
       Do m2=1, 2 *l_in +1

          ztmp = 0.0d0

          Do ma = m_min, m_max
             Do mb = m_min, m_max
                Do ixyz=1, 3
#if 0
                   ztmp(ixyz) = ztmp(ixyz) + conjg(MatU_ylm_RC( m1,ma ))&
                        &                *Mat_L_with_cmplx_ylm( ma,mb,ixyz ) &
                        &                *( MatU_ylm_RC( m2,mb ) )
#else
                   ztmp(ixyz) = ztmp(ixyz) + MatU_ylm_RC( m1,ma )&
                        &                *Mat_L_with_cmplx_ylm( ma,mb,ixyz ) &
                        &                *conjg( MatU_ylm_RC( m2,mb ) )
#endif

                End do
             End do
          End do
          Mat_L_with_real_ylm( m1,m2,: ) = ztmp(:)
       End do
    End do
  end subroutine m_OP_calc_MatL_with_Real_ylm

  subroutine m_OP_calc_MatL_with_Cmplx_ylm( l_in, theta, phi, Mat_L_term )
    integer, intent(in) :: l_in
    real(kind=DP), intent(in) :: theta, phi
    complex(kind=CMPLDP), intent(out) :: &
         &          Mat_L_term( -l_in:l_in, -l_in:l_in, 3 )

    real(kind=DP) :: cos_th, sin_th, cos2_th_h, sin2_th_h
    real(kind=DP) :: ctmp_m, ctmp_p, ctmp_0

    integer :: m1, m2

    Do m1=-l_in, l_in
       Do m2=-l_in, l_in

          if ( m2 == m1 +1 ) then
             ctmp_m = sqrt( (l_in +m1 +1.0d0) *( l_in -m1 ) )
          else
             ctmp_m = 0.0d0
          endif

          if ( m2 == m1 -1 ) then
             ctmp_p = sqrt( (l_in -m1 +1.0d0) *( l_in +m1 ) )
          else
            ctmp_p = 0.0d0
          endif
          !
          if ( m1 == m2 ) then
             ctmp_0 = dble(m1)
          else
             ctmp_0 = 0.0d0
          endif

          Mat_L_term( m1,m2,1 ) =  ( ctmp_p +ctmp_m ) /2.0d0
          Mat_L_term( m1,m2,2 ) =  ( ctmp_p -ctmp_m ) /2.0d0 /zi
          Mat_L_term( m1,m2,3 ) =  ctmp_0

       End do
    End Do

  end subroutine m_OP_calc_MatL_with_Cmplx_ylm

  subroutine m_OP_print_OrbMagMom_on_atom(nfout)
    integer, intent(in) :: nfout

    integer :: ia, it

    if ( mype /= 0 ) return

    write(nfout,*) &
         & '! ------------ Local Orbital Momemnt (in sphere) at this scf step --- '

    write(nfout,*) &
         & '    id      name         mx             my             mz'
    Do ia=1, natm
       it = ityp(ia)
       write(nfout,'(I7,6X,A5,3F15.8)') ia, speciesname(it), &
            &                          OrbMag_on_atom(ia,1), &
            &                          OrbMag_on_atom(ia,2), &
            &                          OrbMag_on_atom(ia,3)
    End Do
    write(nfout,*) '! ---------------------------------------------'

  end subroutine m_OP_print_OrbMagMom_on_atom

! ==== KT_add === 2014/08/26
  subroutine m_OP_calc_OrbMagMom_in_sphere
    logical :: FirstFlag = .true.

    if ( FirstFlag ) then
       call m_OP_calc_MatL_orb_s_to_f
       FirstFlag = .false.
    endif

    if ( kimg == 1 ) return
    if ( .not. noncol ) return

    call alloc_arrays

#ifdef USE_ASMS_INTEGRAL
    call calc_ylm_dot_wfns_fast
#else
    call calc_ylm_dot_wfns_exact
#endif

    call calc_rho_ylm1_ylm2
    if ( flg_paw .or. sw_mix_charge_hardpart == ON ) then
       call symmetrz_rho_ylm1_ylm2
    endif

    if ( .not. allocated( OrbMag_on_atom ) ) allocate( OrbMag_on_atom(natm,3) )
    OrbMag_on_atom = 0.0d0

    call calc_soft_part_contrib
    call calc_hard_part_contrib

    call dealloc_arrays
!    deallocate( OrbMag_on_atom )

!    call m_OP_print_OrbMagMom_on_atom(nfout)

  contains

    subroutine alloc_arrays
      allocate( wfn_sph_decomposed( np_e, nlmta, ista_k:iend_k ) )
      wfn_sph_decomposed = 0.0d0

      allocate( rho_ylm1_ylm2_r( natm, nlmt, nlmt, ndim_magmom ) )
      allocate( rho_ylm1_ylm2_i( natm, nlmt, nlmt, ndim_magmom ) )
      rho_ylm1_ylm2_r = 0.0d0;  rho_ylm1_ylm2_i = 0.0d0

    end subroutine alloc_arrays

    subroutine dealloc_arrays
      deallocate( wfn_sph_decomposed )
      deallocate( rho_ylm1_ylm2_r );  deallocate( rho_ylm1_ylm2_i )
    end subroutine dealloc_arrays

    subroutine calc_ylm_dot_wfns_exact
      integer :: ia, it, ik, ib, ie, ig
      integer :: lmt1, il1, im1, it1, ilmta
      integer :: nspher, ir
      real(kind=DP) :: rad1, c1, facr
      complex(kind=CMPLDP) :: zfac1, zsum, z2, csum

      real(kind=DP), allocatable :: qx(:), qy(:), qz(:), vlength(:)
      real(kind=DP), allocatable :: ylm(:), wka(:), wkb(:)
      real(kind=DP), allocatable :: rwork(:)
      complex(kind=CMPLDP), allocatable :: zph(:)

      allocate( radr(mmesh) );  allocate( wos(mmesh) )

      allocate( qx(kg1) ); allocate( qy(kg1) ); allocate( qz(kg1) );
      allocate( vlength(kg1) ); allocate( ylm(kg1) )
      allocate( wka(kg1) ); allocate( wkb(kg1) )

      allocate( zph(kg1) ); allocate( rwork(kg1) )

      Do ia=1, natm
         it = ityp(ia);  rad1 = rad_cov(ia)
         call new_radr_and_wos(ista_k,it)

         csum = 0.0d0
         Do ir=1, nmesh(it)
            if ( radr(ir) > rad_cov(ia) ) exit
!            csum = csum +wos(ir) *radr(ir)
            csum = csum +wos(ir) *radr(ir)**2
         End do
         csum = sqrt( csum )

         do lmt1=1,ilmt(it)
            il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
            if ( it1 /=1 ) cycle

            ilmta = lmta(lmt1,ia)

            nspher = (il1-1)**2 + im1
            zfac1 = zi **(il1-1)

            Do ik=1, kv3
               if ( map_k(ik) /= myrank_k ) cycle

               call k_plus_G_vectors( ik, kgp, kg1, kv3, iba, nbase, vkxyz, &
                    &                 ngabc, rltv, qx, qy, qz, vlength )
               call sphr( iba(ik), nspher, qx, qy, qz, ylm )

               Do ig=1, iba(ik)
                  c1 = qx(ig)*cps(ia,1) +qy(ig)*cps(ia,2) +qz(ig)*cps(ia,3)
                  zph(ig) = dcmplx( cos(c1), sin(c1) )
               End do

               rwork = 0.0d0

               Do ir=1, nmesh(it)
                  if ( radr(ir) > rad_cov(ia) ) exit
!                  if ( radr(ir) > 2.5 ) exit
!!                  if ( ir > wf_mnrc(it) ) exit

                  facr = wos(ir) *radr(ir)**2
!!!                  facr = wos(ir) *radr(ir)

                  Do ig=1, iba(ik)
                     wka(ig) = vlength(ig) *radr(ir)
                  End do
                  call dsjnv( il1-1, iba(ik), wka, wkb )

                  Do ig=1, iba(ik)
                     rwork(ig) = rwork(ig) + facr *wkb(ig) *ylm(ig)
                  End do
               End do
!
               Do ib=1, neg
                  if ( map_e(ib) /= myrank_e ) cycle
                  ie = map_z(ib)

                  zsum = 0.0d0
                  Do ig=1, iba(ik)
                     z2 = dcmplx( zaj_l(ig,ie,ik,1), zaj_l(ig,ie,ik,kimg) )
                     zsum = zsum + z2 *zph(ig) *rwork(ig)
                  End do
                  wfn_sph_decomposed(ie,ilmta,ik) = zsum *zfac1 /csum
               End do
            End Do
         End do
      End Do
!      stop
!
      wfn_sph_decomposed = wfn_sph_decomposed *PAI4 /sqrt(univol)

      deallocate( zph ); deallocate( radr ); deallocate( wos )
      deallocate( qx ); deallocate( qy ); deallocate( qz ); deallocate( vlength )
      deallocate( ylm )
      deallocate( wka ); deallocate( wkb )
    end subroutine calc_ylm_dot_wfns_exact

#ifdef USE_ASMS_INTEGRAL
    subroutine calc_ylm_dot_wfns_fast
      interface
         real(8) function ASMS_integ_spherical_bessel_j(a,b,c)
           integer :: a
           real(8) :: b, c
         end function ASMS_integ_spherical_bessel_j
      end interface

      integer :: ia, it, ik, ib, ie, ig
      integer :: lmt1, il1, im1, it1, ilmta
      integer :: nspher, ir
      real(kind=DP) :: rad1, c1, facr, bfact
      complex(kind=CMPLDP) :: zfac1, zsum, z2

      real(kind=DP), allocatable :: qx(:), qy(:), qz(:), vlength(:)
      real(kind=DP), allocatable :: ylm(:)
      real(kind=DP), allocatable :: sph_besselj_integ(:)
      complex(kind=CMPLDP), allocatable :: zph(:)

      allocate( qx(kg1) ); allocate( qy(kg1) ); allocate( qz(kg1) );
      allocate( vlength(kg1) );
      allocate( ylm(kg1) );     allocate( zph(kg1) );
      allocate( sph_besselj_integ(kg1) )

      Do ia=1, natm
         it = ityp(ia);  rad1 = rad_cov(ia)

         do lmt1=1,ilmt(it)
            il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
            ilmta = lmta(lmt1,ia)

            nspher = (il1-1)**2 + im1
            zfac1 = zi **(il1-1)

            Do ik=1, kv3
               if ( map_k(ik) /= myrank_k ) cycle

               call k_plus_G_vectors( ik, kgp, kg1, kv3, iba, nbase, vkxyz, &
                    &                 ngabc, rltv, qx, qy, qz, vlength )
               call sphr( iba(ik), nspher, qx, qy, qz, ylm )

               Do ig=1, iba(ik)
                  c1 = qx(ig)*cps(ia,1) +qy(ig)*cps(ia,2) +qz(ig)*cps(ia,3)
                  zph(ig) = dcmplx( cos(c1), sin(c1) )

                  sph_besselj_integ(ig) &
                       &   = ASMS_integ_spherical_bessel_j( il1 -1, vlength(ig), rad1 )
               End do

               Do ib=1, neg
                  if ( map_e(ib) /= myrank_e ) cycle
                  ie = map_z(ib)

                  zsum = 0.0d0
                  Do ig=1, iba(ik)
                     z2 = dcmplx( zaj_l(ig,ie,ik,1), zaj_l(ig,ie,ik,kimg) )
                     zsum = zsum + z2 *zph(ig) *ylm(ig) *sph_besselj_integ(ig)
                  End do
                  wfn_sph_decomposed(ie,ilmta,ik) = zsum *zfac1
               End do
            End Do
         End do
      End Do
!
      wfn_sph_decomposed = wfn_sph_decomposed *PAI4 /sqrt(univol)

      deallocate( qx ); deallocate( qy ); deallocate( qz ); deallocate( vlength )
      deallocate( ylm );deallocate( zph );
      deallocate( sph_besselj_integ )

    end subroutine calc_ylm_dot_wfns_fast
#endif

    subroutine calc_soft_part_contrib
      integer :: ia, it, lmt1, lmt2, is1, is2, istmp
      integer :: il1, il2, im1, im2, it1, it2
      real(kind=DP), allocatable :: hsr_ssrep(:,:,:,:), hsi_ssrep(:,:,:,:)
      real(kind=DP) :: orbmom(3)
      complex(kind=CMPLDP) :: ztmp(3), z1, z2

      allocate( hsr_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsr_ssrep = 0.0d0
      allocate( hsi_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsi_ssrep = 0.0d0

      call m_ES_MagMom_To_DensMat_hsr( natm, nlmt, rho_ylm1_ylm2_r, rho_ylm1_ylm2_i, &
           &                           hsr_ssrep, hsi_ssrep )

      Do ia=1, natm
         it = ityp(ia)

         orbmom = 0.0d0

         Do is1=1, ndim_spinor
            Do is2=1, ndim_spinor
               if ( is1 /= is2 ) cycle
               istmp = (is1 -1)*ndim_spinor +is2

               do lmt1=1,ilmt(it)
                  il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
                  if ( it1 > 1 ) cycle
                  if ( il1 <= 1 .or. il1 >4 ) cycle

                  do lmt2=1,ilmt(it)
                     il2=ltp(lmt2,it);  im2=mtp(lmt2,it); it2=taup(lmt2,it)
                     if ( it2 > 1 ) cycle
                     if ( il1 /= il2 ) cycle

                     ztmp = 0.0d0
                     if ( il1 == 2 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L1( im2, im1, : )
                     else if ( il1 == 3 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L2( im2, im1, : )
                     else if ( il1 == 4 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L3( im2, im1, : )
                     endif

                     z1 = dcmplx( hsr_ssrep(ia,lmt1,lmt2,istmp), &
                          &       hsi_ssrep(ia,lmt1,lmt2,istmp) )

                     orbmom(1) = orbmom(1) +z1*ztmp(1)
                     orbmom(2) = orbmom(2) +z1*ztmp(2)
                     orbmom(3) = orbmom(3) +z1*ztmp(3)
                  End do
               End Do
!
            End do
         End do

         OrbMag_on_atom(ia,1) = OrbMag_on_atom(ia,1) +orbmom(1)
         OrbMag_on_atom(ia,2) = OrbMag_on_atom(ia,2) +orbmom(2)
         OrbMag_on_atom(ia,3) = OrbMag_on_atom(ia,3) +orbmom(3)
      End Do
      deallocate( hsr_ssrep, hsi_ssrep )
    end subroutine calc_soft_part_contrib

    subroutine calc_hard_part_contrib
      integer :: ia, it, lmt1, lmt2, is1, is2, istmp
      integer :: il1, il2, im1, im2, it1, it2
      real(kind=DP), allocatable :: hsr_ssrep(:,:,:,:), hsi_ssrep(:,:,:,:)
      real(kind=DP), allocatable :: ovp(:,:,:,:)
      real(kind=DP) :: orbmom(3)
      complex(kind=CMPLDP) :: ztmp(3), z1, z2

      allocate( hsr_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsr_ssrep = 0.0d0
      allocate( hsi_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsi_ssrep = 0.0d0

      call m_ES_MagMom_To_DensMat_hsr( natm, nlmt, hsr, hsi, hsr_ssrep, hsi_ssrep )

      allocate( ovp(ntau,ntau,nloc,ntyp ) );  ovp = 0.0d0
      Do it=1, ntyp
         do lmt1=1,ilmt(it)
            il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
            do lmt2=1,ilmt(it)
               il2=ltp(lmt2,it);  im2=mtp(lmt2,it); it2=taup(lmt2,it)
               if ( il1 == il2 .and. im1 == im2 ) then
                  ovp(it1,it2,il1,it) = q_noncl(lmt1,lmt2,1,it)
               endif
            end do
         end do
      end Do

      Do ia=1, natm
         it = ityp(ia)

         orbmom = 0.0d0

         Do is1=1, ndim_spinor
            Do is2=1, ndim_spinor
               if ( is1 /= is2 ) cycle
               istmp = (is1 -1)*ndim_spinor +is2

               do lmt1=1,ilmt(it)
                  il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
                  do lmt2=1,ilmt(it)
                     il2=ltp(lmt2,it);  im2=mtp(lmt2,it); it2=taup(lmt2,it)

                     if ( il1 <= 1 .or. il1 >4 ) cycle
                     if ( il1 /= il2 ) cycle

                     ztmp = 0.0d0
                     if ( il1 == 2 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L1( im2, im1, : )
                     else if ( il1 == 3 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L2( im2, im1, : )
                     else if ( il1 == 4 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L3( im2, im1, : )
                     endif

                     z1 = dcmplx( hsr_ssrep(ia,lmt1,lmt2,istmp), &
                          &       hsi_ssrep(ia,lmt1,lmt2,istmp) )
#if 0
!                     z2 = z1 *q(lmt1,lmt2,it)
                     z2 = z1 *q_noncl(lmt1,lmt2,istmp,it)
#else
                     z2 = z1 *ovp(it1,it2,il1,it)
#endif

                     orbmom(1) = orbmom(1) +z2*ztmp(1)
                     orbmom(2) = orbmom(2) +z2*ztmp(2)
                     orbmom(3) = orbmom(3) +z2*ztmp(3)
                  End do
               End Do
!
            End do
         End do
         OrbMag_on_atom(ia,1) = OrbMag_on_atom(ia,1) +orbmom(1)
         OrbMag_on_atom(ia,2) = OrbMag_on_atom(ia,2) +orbmom(2)
         OrbMag_on_atom(ia,3) = OrbMag_on_atom(ia,3) +orbmom(3)
      End Do
      deallocate( hsr_ssrep, hsi_ssrep )
      deallocate( ovp )

    end subroutine calc_hard_part_contrib

  end subroutine m_OP_calc_OrbMagMom_in_sphere

  subroutine m_OP_calc_OrbMagMom_method1
    logical :: FirstFlag = .true.

    if ( FirstFlag ) then
       call m_OP_calc_MatL_orb_s_to_f
       FirstFlag = .false.
    endif

    if ( kimg == 1 ) return
    if ( .not. noncol ) return

    call alloc_arrays

    call calc_psir_ylm_dot_wfns

    call calc_rho_ylm1_ylm2
    if ( flg_paw .or. sw_mix_charge_hardpart == ON ) then
       call symmetrz_rho_ylm1_ylm2
    endif

    if ( .not. allocated( OrbMag_on_atom ) ) allocate( OrbMag_on_atom(natm,3) )
    OrbMag_on_atom = 0.0d0

    call calc_soft_part_contrib
    call calc_hard_part_contrib

    call print_orb_mom

    call dealloc_arrays
!    deallocate( OrbMag_on_atom )

  contains

    subroutine alloc_arrays
      allocate( wfn_sph_decomposed( np_e, nlmta, ista_k:iend_k ) )
      wfn_sph_decomposed = 0.0d0

      allocate( rho_ylm1_ylm2_r( natm, nlmt, nlmt, ndim_magmom ) )
      allocate( rho_ylm1_ylm2_i( natm, nlmt, nlmt, ndim_magmom ) )
      rho_ylm1_ylm2_r = 0.0d0;  rho_ylm1_ylm2_i = 0.0d0

    end subroutine alloc_arrays

    subroutine dealloc_arrays
      deallocate( wfn_sph_decomposed )
      deallocate( rho_ylm1_ylm2_r );  deallocate( rho_ylm1_ylm2_i )
    end subroutine dealloc_arrays

    subroutine calc_psir_ylm_dot_wfns
      integer :: ia, it, ik, ib, ie, ig
      integer :: lmt1, il1, im1, it1, ilmta
      integer :: nspher, ir
      real(kind=DP) :: rad1, c1, facr
      complex(kind=CMPLDP) :: zfac1, zsum, z2

      real(kind=DP), allocatable :: qx(:), qy(:), qz(:), vlength(:)
      real(kind=DP), allocatable :: ylm(:), wka(:), wkb(:)
      real(kind=DP), allocatable :: rwork(:)
      complex(kind=CMPLDP), allocatable :: zph(:)

      allocate( radr(mmesh) );  allocate( wos(mmesh) )

      allocate( qx(kg1) ); allocate( qy(kg1) ); allocate( qz(kg1) );
      allocate( vlength(kg1) ); allocate( ylm(kg1) )
      allocate( wka(kg1) ); allocate( wkb(kg1) )

      allocate( zph(kg1) ); allocate( rwork(kg1) )

      Do it=1, ntyp
         call new_radr_and_wos(ista_k,it)
         rad1 = rad_cov_default( nint(iatomn(it)) ) ! Revised according to a report from ASMS Co.ltd, 10 March 2016.

         do lmt1=1,ilmt(it)
            il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)

            nspher = (il1-1)**2 +im1;
            zfac1 = zi **(il1-1)

            Do ik=1, kv3
               if ( map_k(ik) /= myrank_k ) cycle

               call k_plus_G_vectors( ik, kgp, kg1, kv3, iba, nbase, vkxyz, &
                    &                 ngabc, rltv, qx, qy, qz, vlength )
               call sphr( iba(ik), nspher, qx, qy, qz, ylm )

               rwork = 0.0d0

               Do ir=1, nmesh(it)
!                  if ( ir > wf_mnrc(it) ) exit
                  if ( radr(ir) > rad1 ) exit

!                  facr = wos(ir) *psirpw(ir,il1,it1,it) *radr(ir)
                  facr = wos(ir) *psir_val(ir,il1,it1,it) *radr(ir)
!                  facr = wos(ir) *phirpw(ir,il1,it1,it) *radr(ir)

                  Do ig=1, iba(ik)
                     wka(ig) = vlength(ig) *radr(ir)
                  End do
                  call dsjnv( il1-1, iba(ik), wka, wkb )

                  Do ig=1, iba(ik)
                     rwork(ig) = rwork(ig) + facr *wkb(ig) *ylm(ig)
                  End do
               End do

               Do ia=1, natm
                  if ( ityp(ia) /= it ) cycle
                  ilmta = lmta(lmt1,ia)

                  Do ig=1, iba(ik)
                     c1 = qx(ig)*cps(ia,1) +qy(ig)*cps(ia,2) +qz(ig)*cps(ia,3)
                     zph(ig) = dcmplx( cos(c1), sin(c1) )
                  End do

!
                  Do ib=1, neg
                     if ( map_e(ib) /= myrank_e ) cycle
                     ie = map_z(ib)

                     zsum = 0.0d0
                     Do ig=1, iba(ik)
                        z2 = dcmplx( zaj_l(ig,ie,ik,1), zaj_l(ig,ie,ik,kimg) )
                        zsum = zsum + z2 *zph(ig) *rwork(ig)
                     End do
                     wfn_sph_decomposed(ie,ilmta,ik) = zsum *zfac1
                  End do
               End do
            End Do
         End do
      End Do
!
      wfn_sph_decomposed = wfn_sph_decomposed *PAI4 /sqrt(univol)

      deallocate( zph ); deallocate( radr ); deallocate( wos )
      deallocate( qx ); deallocate( qy ); deallocate( qz ); deallocate( vlength )
      deallocate( ylm )
      deallocate( wka ); deallocate( wkb )
    end subroutine calc_psir_ylm_dot_wfns

    subroutine calc_soft_part_contrib
      integer :: ia, it, lmt1, lmt2, is1, is2, istmp
      integer :: il1, il2, im1, im2, it1, it2
      real(kind=DP), allocatable :: hsr_ssrep(:,:,:,:), hsi_ssrep(:,:,:,:)
      real(kind=DP) :: orbmom(3)
      complex(kind=CMPLDP) :: ztmp(3), z1, z2

      allocate( hsr_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsr_ssrep = 0.0d0
      allocate( hsi_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsi_ssrep = 0.0d0

      call m_ES_MagMom_To_DensMat_hsr( natm, nlmt, rho_ylm1_ylm2_r, rho_ylm1_ylm2_i, &
           &                           hsr_ssrep, hsi_ssrep )

      Do ia=1, natm
         it = ityp(ia)

         orbmom = 0.0d0

         Do is1=1, ndim_spinor
            Do is2=1, ndim_spinor
               if ( is1 /= is2 ) cycle
               istmp = (is1 -1)*ndim_spinor +is2

               do lmt1=1,ilmt(it)
                  il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
                  if ( il1 <= 1 .or. il1 >4 ) cycle
!                  if (it1 >1 ) cycle

                  do lmt2=1,ilmt(it)
                     il2=ltp(lmt2,it);  im2=mtp(lmt2,it); it2=taup(lmt2,it)
                     if ( il1 /= il2 ) cycle
!                     if (it2 >1 ) cycle

                     ztmp = 0.0d0
                     if ( il1 == 2 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L1( im2, im1, : )
                     else if ( il1 == 3 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L2( im2, im1, : )
                     else if ( il1 == 4 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L3( im2, im1, : )
                     endif

                     z1 = dcmplx( hsr_ssrep(ia,lmt1,lmt2,istmp), &
                          &       hsi_ssrep(ia,lmt1,lmt2,istmp) )
                     z1 = z1 *product_psir_val(il1,it1,it2,it)

                     orbmom(1) = orbmom(1) +z1*ztmp(1)
                     orbmom(2) = orbmom(2) +z1*ztmp(2)
                     orbmom(3) = orbmom(3) +z1*ztmp(3)
                  End do
               End Do
!
            End do
         End do

         OrbMag_on_atom(ia,1) = OrbMag_on_atom(ia,1) +orbmom(1)
         OrbMag_on_atom(ia,2) = OrbMag_on_atom(ia,2) +orbmom(2)
         OrbMag_on_atom(ia,3) = OrbMag_on_atom(ia,3) +orbmom(3)
      End Do
      deallocate( hsr_ssrep, hsi_ssrep )
    end subroutine calc_soft_part_contrib

    subroutine calc_hard_part_contrib
      integer :: ia, it, lmt1, lmt2, is1, is2, istmp
      integer :: il1, il2, im1, im2, it1, it2
      real(kind=DP), allocatable :: hsr_ssrep(:,:,:,:), hsi_ssrep(:,:,:,:)
      real(kind=DP), allocatable :: ovp(:,:,:,:)
      real(kind=DP) :: orbmom(3)
      complex(kind=CMPLDP) :: ztmp(3), z1, z2

      allocate( hsr_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsr_ssrep = 0.0d0
      allocate( hsi_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsi_ssrep = 0.0d0

      call m_ES_MagMom_To_DensMat_hsr( natm, nlmt, hsr, hsi, hsr_ssrep, hsi_ssrep )

      allocate( ovp(ntau,ntau,nloc,ntyp ) );  ovp = 0.0d0
      Do it=1, ntyp
         do lmt1=1,ilmt(it)
            il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
            do lmt2=1,ilmt(it)
               il2=ltp(lmt2,it);  im2=mtp(lmt2,it); it2=taup(lmt2,it)
               if ( il1 == il2 .and. im1 == im2 ) then
                  ovp(it1,it2,il1,it) = q_noncl(lmt1,lmt2,1,it)
               endif
            end do
         end do
      end Do

      Do ia=1, natm
         it = ityp(ia)

         orbmom = 0.0d0

         Do is1=1, ndim_spinor
            Do is2=1, ndim_spinor
               if ( is1 /= is2 ) cycle
               istmp = (is1 -1)*ndim_spinor +is2

               do lmt1=1,ilmt(it)
                  il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
                  do lmt2=1,ilmt(it)
                     il2=ltp(lmt2,it);  im2=mtp(lmt2,it); it2=taup(lmt2,it)

                     if ( il1 <= 1 .or. il1 >4 ) cycle
                     if ( il1 /= il2 ) cycle

                     ztmp = 0.0d0
                     if ( il1 == 2 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L1( im2, im1, : )
                     else if ( il1 == 3 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L2( im2, im1, : )
                     else if ( il1 == 4 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L3( im2, im1, : )
                     endif

                     z1 = dcmplx( hsr_ssrep(ia,lmt1,lmt2,istmp), &
                          &       hsi_ssrep(ia,lmt1,lmt2,istmp) )
#if 0
!                     z2 = z1 *q(lmt1,lmt2,it)
                     z2 = z1 *q_noncl(lmt1,lmt2,istmp,it)
#else
                     z2 = z1 *ovp(it1,it2,il1,it)
#endif
                     orbmom(1) = orbmom(1) +z2*ztmp(1)
                     orbmom(2) = orbmom(2) +z2*ztmp(2)
                     orbmom(3) = orbmom(3) +z2*ztmp(3)
                  End do
               End Do
!
            End do
         End do
         OrbMag_on_atom(ia,1) = OrbMag_on_atom(ia,1) +orbmom(1)
         OrbMag_on_atom(ia,2) = OrbMag_on_atom(ia,2) +orbmom(2)
         OrbMag_on_atom(ia,3) = OrbMag_on_atom(ia,3) +orbmom(3)
      End Do
      deallocate( hsr_ssrep, hsi_ssrep )
      deallocate( ovp )

    end subroutine calc_hard_part_contrib

    subroutine print_orb_mom
      integer :: ia, it, k

      if ( mype /= 0 ) return

      write(nfout,*) &
           & '! ------------ Local Orbital Momemnt (method1) at this scf step --- '

       write(nfout,*) &
            & '!   id   atom no.        mx             my             mz'
       Do ia=1, natm
          it = ityp(ia)
          write(nfout,'(I7,F10.4,3F15.8)') ia, iatomn(it), &
               &                          OrbMag_on_atom(ia,1), &
               &                          OrbMag_on_atom(ia,2), &
               &                          OrbMag_on_atom(ia,3)
       End Do
       write(nfout,*) '! ---------------------------------------------'

    end subroutine print_orb_mom

  end subroutine m_OP_calc_OrbMagMom_method1

  subroutine m_OP_calc_OrbMagMom_method2
    logical :: FirstFlag = .true.

    if ( iprimagmom < 3 ) return

    if ( FirstFlag ) then
       call m_OP_calc_MatL_orb_s_to_f
       FirstFlag = .false.
    endif

    if ( kimg == 1 ) return
    if ( .not. noncol ) return

    if ( .not. allocated( OrbMag_on_atom ) ) allocate( OrbMag_on_atom(natm,3) )
    OrbMag_on_atom = 0.0d0

    call hardpart_contrib
    if ( sw_use_add_proj == ON ) call hardpart_contrib_add_proj

    call print_orb_mom

!    deallocate( OrbMag_on_atom )

  contains

    subroutine hardpart_contrib
      integer :: ia, it, lmt1, lmt2, is1, is2, istmp
      integer :: il1, il2, im1, im2, it1, it2
      real(kind=DP), allocatable :: hsr_ssrep(:,:,:,:), hsi_ssrep(:,:,:,:)
      real(kind=DP) :: orbmom(3)
      complex(kind=CMPLDP) :: ztmp(3), z1, z2

      allocate( hsr_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsr_ssrep = 0.0d0
      allocate( hsi_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsi_ssrep = 0.0d0

      call m_ES_MagMom_To_DensMat_hsr( natm, nlmt, hsr, hsi, hsr_ssrep, hsi_ssrep )

      Do ia=1, natm
         it = ityp(ia)

         orbmom = 0.0d0

         Do is1=1, ndim_spinor
            Do is2=1, ndim_spinor
               if ( is1 /= is2 ) cycle
               istmp = (is1 -1)*ndim_spinor +is2

               do lmt1=1,ilmt(it)
                  il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
                  do lmt2=1,ilmt(it)
                     il2=ltp(lmt2,it);  im2=mtp(lmt2,it); it2=taup(lmt2,it)

                     if ( il1 <= 1 .or. il1 >4 ) cycle
                     if ( il1 /= il2 ) cycle

                     ztmp = 0.0d0
                     if ( il1 == 2 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L1( im2, im1, : )
                     else if ( il1 == 3 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L2( im2, im1, : )
                     else if ( il1 == 4 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L3( im2, im1, : )
                     endif

                     z1 = dcmplx( hsr_ssrep(ia,lmt1,lmt2,istmp), &
                          &       hsi_ssrep(ia,lmt1,lmt2,istmp) )
                     if ( sw_use_contracted_psir == ON ) then
                        z1 = z1 *product_psir_val( il1, it1, it2, it )
                     endif

                     orbmom(1) = orbmom(1) +z1*ztmp(1)
                     orbmom(2) = orbmom(2) +z1*ztmp(2)
                     orbmom(3) = orbmom(3) +z1*ztmp(3)
                  End do
               End Do
!
            End do
         End do

         OrbMag_on_atom(ia,1) = OrbMag_on_atom(ia,1) +orbmom(1)
         OrbMag_on_atom(ia,2) = OrbMag_on_atom(ia,2) +orbmom(2)
         OrbMag_on_atom(ia,3) = OrbMag_on_atom(ia,3) +orbmom(3)
      End Do
      deallocate( hsr_ssrep, hsi_ssrep )

    end subroutine hardpart_contrib

    subroutine hardpart_contrib_add_proj
      integer :: ia, it, lmt1, lmt2, is1, is2, istmp
      integer :: il1, il2, im1, im2, it1, it2
      real(kind=DP), allocatable :: hsr_ssrep(:,:,:,:), hsi_ssrep(:,:,:,:)
      real(kind=DP) :: orbmom(3)
      complex(kind=CMPLDP) :: ztmp(3), z1, z2

      allocate( hsr_ssrep( natm, nlmt_add, nlmt_add, ndim_chgpot ) ); hsr_ssrep = 0.0d0
      allocate( hsi_ssrep( natm, nlmt_add, nlmt_add, ndim_chgpot ) ); hsi_ssrep = 0.0d0

      call m_ES_MagMom_To_DensMat_hsr( natm, nlmt_add, hsr_add, hsi_add, &
           &                           hsr_ssrep, hsi_ssrep )

      Do ia=1, natm
         it = ityp(ia)

         orbmom = 0.0d0

         Do is1=1, ndim_spinor
            Do is2=1, ndim_spinor
               if ( is1 /= is2 ) cycle
               istmp = (is1 -1)*ndim_spinor +is2

               do lmt1=1,ilmt_add(it)
                  il1=ltp_add(lmt1,it); im1=mtp_add(lmt1,it); it1=1
                  do lmt2=1,ilmt_add(it)
                     il2=ltp_add(lmt2,it);  im2=mtp_add(lmt2,it); it2=1

                     if ( il1 <= 1 .or. il1 >4 ) cycle
                     if ( il1 /= il2 ) cycle

                     ztmp = 0.0d0
                     if ( il1 == 2 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L1( im2, im1, : )
                     else if ( il1 == 3 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L2( im2, im1, : )
                     else if ( il1 == 4 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L3( im2, im1, : )
                     endif

                     z1 = dcmplx( hsr_ssrep(ia,lmt1,lmt2,istmp), &
                          &       hsi_ssrep(ia,lmt1,lmt2,istmp) )

                     orbmom(1) = orbmom(1) +z1*ztmp(1)
                     orbmom(2) = orbmom(2) +z1*ztmp(2)
                     orbmom(3) = orbmom(3) +z1*ztmp(3)
                  End do
               End Do
!
            End do
         End do

         OrbMag_on_atom(ia,1) = OrbMag_on_atom(ia,1) +orbmom(1)
         OrbMag_on_atom(ia,2) = OrbMag_on_atom(ia,2) +orbmom(2)
         OrbMag_on_atom(ia,3) = OrbMag_on_atom(ia,3) +orbmom(3)
      End Do
      deallocate( hsr_ssrep, hsi_ssrep )

    end subroutine hardpart_contrib_add_proj

    subroutine print_orb_mom
      integer :: ia, it, k

      if ( mype /= 0 ) return

      write(nfout,*) &
           & '! ------------ Local Orbital Momemnt (method2) at this scf step --- '
       write(nfout,*) &
            & '    id      name          mx             my             mz'
       Do ia=1, natm
          it = ityp(ia)
          write(nfout,'(I7,6X,A5,3F15.8)') ia, speciesname(it), &
               &                          OrbMag_on_atom(ia,1), &
               &                          OrbMag_on_atom(ia,2), &
               &                          OrbMag_on_atom(ia,3)
       End Do

       write(nfout,*)

    end subroutine print_orb_mom

  end subroutine m_OP_calc_OrbMagMom_method2

  subroutine calc_rho_ylm1_ylm2
    integer         :: ia, is, k, i, lmt1, lmt2, p, q, it
    integer :: is1, is2, is_tmp, k1, k2

    real(kind=DP)   :: w_n, d_factor
    complex(kind=CMPLDP) :: z1

    integer :: il1, il2, im1, im2, it1, it2, istmp, jstmp

    real(kind=DP), allocatable, dimension(:,:,:,:) :: hsr_or_hsi_mpi
    real(kind=DP), allocatable, dimension(:,:,:,:) :: hsr_ssrep
    real(kind=DP), allocatable, dimension(:,:,:,:) :: hsi_ssrep
!
! ---------------------- kt : uncertain ---------
!!!!!    d_factor = 2.d0/ ( kv3 /ndim_spinor )
!
    d_factor = 1.d0/ ( kv3 /ndim_spinor )
! -------------------------------------------

    allocate( hsr_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsr_ssrep = 0.0d0
    allocate( hsi_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsi_ssrep = 0.0d0

    do ia = 1, natm
       it = ityp(ia)

       do is1 = 1, ndim_spinor
          do is2 = 1, ndim_spinor
             is_tmp = ( is1 -1 )*ndim_spinor + is2
             do i = 1, np_e                                  ! MPI

                do k = 1, kv3, ndim_spinor
                   if ( map_k(k) /= myrank_k ) cycle            ! MPI
                   w_n = occup_l(i,k) *d_factor

                   k1 = k + is1 -1;  k2 = k + is2 -1
                   do lmt1 = 1, ilmt(it)
                      p = lmta(lmt1,ia)
                      do lmt2 = 1, ilmt(it)
                         q = lmta(lmt2,ia)
                         z1 = wfn_sph_decomposed(i,p,k1) &
                              &      *conjg( wfn_sph_decomposed(i,q,k2) )

                         hsr_ssrep(ia,lmt1,lmt2,is_tmp) &
                              &  = hsr_ssrep(ia,lmt1,lmt2,is_tmp) &
                              &    + w_n * real(z1)
                         hsi_ssrep(ia,lmt1,lmt2,is_tmp) &
                              &  = hsi_ssrep(ia,lmt1,lmt2,is_tmp) &
                              &    + w_n * aimag(z1)
                      end do! lmt2
                   end do! lmt1
                end do! k
             end do! i
          end do! is2
       end do! is1
    end do! ia
!
    if (npes >= 2) then
       allocate(hsr_or_hsi_mpi(natm,nlmt,nlmt,ndim_chgpot))
       hsr_or_hsi_mpi = 0.0d0

       call mpi_allreduce( hsr_ssrep, hsr_or_hsi_mpi, natm*nlmt*nlmt*ndim_chgpot, &
            &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
       hsr_ssrep = hsr_or_hsi_mpi

       call mpi_allreduce( hsi_ssrep, hsr_or_hsi_mpi, natm*nlmt*nlmt*ndim_chgpot, &
            &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
       hsi_ssrep = hsr_or_hsi_mpi

       deallocate(hsr_or_hsi_mpi)
    end if
!
    call m_ES_DensMat_To_MagMom_hsr( natm, nlmt, hsr_ssrep, hsi_ssrep, &
         &                           rho_ylm1_ylm2_r, rho_ylm1_ylm2_i )
!
    deallocate( hsr_ssrep, hsi_ssrep )
!
  end subroutine calc_rho_ylm1_ylm2

  subroutine symmetrz_rho_ylm1_ylm2
    integer         :: ia, is, iopr, i, lmt1, lmt2, lmt3, lmt4, it
    real(kind=DP), allocatable, dimension(:,:,:,:) :: hsr_tmp
    real(kind=DP), allocatable, dimension(:,:,:,:) :: hsi_tmp

    integer :: il1,im1,it1,il2,im2,it2,il3,im3,it3,il4,im4,it4
    integer :: ii,jj,kk,ll,n,m,iii,jjj
    integer :: ja

    integer :: ixyz1, ixyz2, is_tmp
    real(kind=DP) :: ctmp1, weight, ctmp2, weight2, fi

    fi = 1.d0/nopr
    if ( charge_symm_mode >= chg_symm_level1 ) then
       fi = 1.0d0 /dble(nopr_from_fbz_to_ibz)
    endif

    allocate(hsr_tmp(natm,nlmt,nlmt,ndim_magmom)); hsr_tmp = 0.0d0
    allocate(hsi_tmp(natm,nlmt,nlmt,ndim_magmom)); hsi_tmp = 0.0d0
!
    hsr_tmp = rho_ylm1_ylm2_r;  hsi_tmp = rho_ylm1_ylm2_i

    do ia=1,natm
       it=ityp(ia)
       do is =1, ndim_magmom
          do lmt2=1,ilmt(it)
             do lmt1=lmt2+1,ilmt(it)
                hsr_tmp(ia,lmt1,lmt2,is) =  hsr_tmp(ia,lmt2,lmt1,is)
                hsi_tmp(ia,lmt1,lmt2,is) = -hsi_tmp(ia,lmt2,lmt1,is)
             end do
          end do
       end do
    end do

    rho_ylm1_ylm2_r = 0.d0; rho_ylm1_ylm2_i = 0.0d0

    do iopr=1,nopr
       if ( charge_symm_mode >= chg_symm_level1 )then
         if(flg_opr_from_fbz_to_ibz(iopr) == 0 ) cycle
       endif

       if ( magmom_dir_inversion_opr_flag(iopr) == -1 ) then
          weight = -1.0d0
       else
          weight = 1.0d0
       endif

       if ( determinant_op(iopr) > 0 ) then
          weight2 = 1.0d0
       else
          weight2 = -1.0d0
       endif

       do ia = 1, natm
          it = ityp(ia)
          ja=abs(ia2ia_symmtry_op_inv(ia,iopr))

          do lmt1 = 1, ilmt(it)
             il1=ltp(lmt1,it)
             im1=mtp(lmt1,it)
             it1=taup(lmt1,it)
             ii=(il1-1)**2+im1

             do lmt2 = lmt1, ilmt(it)
!             do lmt2 = 1, ilmt(it)
                il2=ltp(lmt2,it)
                im2=mtp(lmt2,it)
                it2=taup(lmt2,it)
                jj=(il2-1)**2+im2

                do n=1,nylm_paw(ii,iopr,ia)
                   iii=iylm_paw(n,ii,iopr,ia)
                   do m=1,nylm_paw(jj,iopr,ia)
                      jjj=iylm_paw(m,jj,iopr,ia)

                      do lmt3=1,ilmt(it)
                         il3=ltp(lmt3,it)
                         im3=mtp(lmt3,it)
                         it3=taup(lmt3,it)
                         kk=(il3-1)**2+im3

                         if(kk.ne.iii .or. it1.ne.it3) cycle

                         do lmt4=1,ilmt(it)
                            il4=ltp(lmt4,it)
                            im4=mtp(lmt4,it)
                            it4=taup(lmt4,it)
                            ll=(il4-1)**2+im4
                            if(ll.ne.jjj .or. it2.ne.it4) cycle

                            ctmp1 = 1.0d0;   ctmp2 = weight

                            rho_ylm1_ylm2_r(ia,lmt1,lmt2,1) = &
                                 rho_ylm1_ylm2_r(ia,lmt1,lmt2,1) + &
                                 ctmp1 * &
                                 hsr_tmp(ja,lmt3,lmt4,1)* &
                                 crotylm_paw(n,ii,iopr,ia)* &
                                 crotylm_paw(m,jj,iopr,ia)
                            rho_ylm1_ylm2_i(ia,lmt1,lmt2,1) = &
                                 rho_ylm1_ylm2_i(ia,lmt1,lmt2,1) + &
                                 ctmp2 * &
                                 hsi_tmp(ja,lmt3,lmt4,1)* &
                                 crotylm_paw(n,ii,iopr,ia)* &
                                 crotylm_paw(m,jj,iopr,ia)

                            Do ixyz1=1, 3
                               Do ixyz2=1, 3
                                  ctmp1 = op(ixyz2, ixyz1, iopr) *weight2 *weight
                                  ctmp2 = op(ixyz2, ixyz1, iopr) *weight2
!
                                  rho_ylm1_ylm2_r(ia,lmt1,lmt2,ixyz2+1) &
                                       & = rho_ylm1_ylm2_r(ia,lmt1,lmt2,ixyz2+1)  &
                                       &  + ctmp1 &
                                       &    *hsr_tmp(ja,lmt3,lmt4,ixyz1+1) &
                                       &    *crotylm_paw(n,ii,iopr,ia)  &
                                       &    *crotylm_paw(m,jj,iopr,ia)
                                  rho_ylm1_ylm2_i(ia,lmt1,lmt2,ixyz2+1) &
                                       & = rho_ylm1_ylm2_i(ia,lmt1,lmt2,ixyz2+1)  &
                                       &  + ctmp2 &
                                       &    *hsi_tmp(ja,lmt3,lmt4,ixyz1+1) &
                                       &    *crotylm_paw(n,ii,iopr,ia)  &
                                       &    *crotylm_paw(m,jj,iopr,ia)
                               End do
                            End do

                         end do! lmt4
                      end do! lmt3

                   end do! jjj
                end do! iii

             end do! lmt2
          end do! lmt1
       end do! ia
    end do! iopr

!    rho_ylm1_ylm2_r = rho_ylm1_ylm2_r/nopr; rho_ylm1_ylm2_i = rho_ylm1_ylm2_i/nopr;
    rho_ylm1_ylm2_r = rho_ylm1_ylm2_r *fi; rho_ylm1_ylm2_i = rho_ylm1_ylm2_i *fi

    do ia=1,natm
       it=ityp(ia)
       do is =1, ndim_magmom
          do lmt2=1,ilmt(it)
             do lmt1=lmt2+1,ilmt(it)
                rho_ylm1_ylm2_r(ia,lmt1,lmt2,is) = rho_ylm1_ylm2_r(ia,lmt2,lmt1,is)
                rho_ylm1_ylm2_i(ia,lmt1,lmt2,is) =-rho_ylm1_ylm2_i(ia,lmt2,lmt1,is)
             end do
          end do
       end do
    end do

    deallocate(hsr_tmp);  deallocate(hsi_tmp)

  end subroutine symmetrz_rho_ylm1_ylm2
! ============= 2014/08/26

! ====== KT_add === 2014/09/01
  subroutine m_OP_orb_decomposed_orb_magmom
    logical :: FirstFlag = .true.

    real(kind=DP), allocatable :: orb_magmom_orb_decomposed(:,:)
    real(kind=DP), allocatable :: orb_magmom_orb_decomposed_add(:,:)

    if ( FirstFlag ) then
       call m_OP_calc_MatL_orb_s_to_f
       FirstFlag = .false.
    endif

    if ( kimg == 1 ) return
    if ( .not. noncol ) return

    allocate( orb_magmom_orb_decomposed( nlmta, 3 ) )
    orb_magmom_orb_decomposed = 0.0d0
    if ( sw_use_add_proj == ON ) then
       allocate( orb_magmom_orb_decomposed_add( nlmta_add, 3 ) )
       orb_magmom_orb_decomposed_add = 0.0d0
    endif

    call contrib_hardpart
    if ( sw_use_add_proj == ON ) call contrib_hardpart_add_proj

    write(nfout,*) &
         &   "! --------- Orbital decomposed orbital magnetic moment (method2) ------"
    write(nfout,*) "< Detail >"
    call print_orb_decomposed_orb_mag_1

    write(nfout,*) "< Summed over m >"
    call print_orb_decomposed_orb_mag_2

    write(nfout,*) "< Summed over l,m >"
    call print_orb_decomposed_orb_mag_3
    write(nfout,*)

    deallocate( orb_magmom_orb_decomposed )
    if ( sw_use_add_proj == ON ) then
       deallocate( orb_magmom_orb_decomposed_add )
    endif

  contains

    subroutine contrib_hardpart
      integer :: ia, it, lmt1, lmt2, is1, is2, istmp
      integer :: il1, il2, im1, im2, it1, it2, ilmta
      real(kind=DP), allocatable :: hsr_ssrep(:,:,:,:), hsi_ssrep(:,:,:,:)
      real(kind=DP) :: ctmp(3)
      complex(kind=CMPLDP) :: ztmp(3), z1, z2

      allocate( hsr_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsr_ssrep = 0.0d0
      allocate( hsi_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsi_ssrep = 0.0d0

      call m_ES_MagMom_To_DensMat_hsr( natm, nlmt, hsr, hsi, hsr_ssrep, hsi_ssrep )

      Do ia=1, natm
         it = ityp(ia)

         Do is1=1, ndim_spinor
            Do is2=1, ndim_spinor
               if ( is1 /= is2 ) cycle
               istmp = (is1 -1)*ndim_spinor +is2

               do lmt1=1,ilmt(it)
                  il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
                  ilmta = lmta(lmt1,ia)

                  do lmt2=1,ilmt(it)
                     il2=ltp(lmt2,it);  im2=mtp(lmt2,it); it2=taup(lmt2,it)

                     if ( il1 <= 1 .or. il1 >4 ) cycle
                     if ( il1 /= il2 ) cycle

                     ztmp = 0.0d0
                     if ( il1 == 2 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L1( im2, im1, : )
                     else if ( il1 == 3 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L2( im2, im1, : )
                     else if ( il1 == 4 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L3( im2, im1, : )
                     endif

                     z1 = dcmplx( hsr_ssrep(ia,lmt1,lmt2,istmp), &
                          &       hsi_ssrep(ia,lmt1,lmt2,istmp) )
                     if ( sw_use_contracted_psir == ON ) then
                        z1 = z1 *product_psir_val( il1, it1, it2, it )
                     endif

                     ctmp(:) = z1*ztmp(:)

                     orb_magmom_orb_decomposed( ilmta,: ) &
                             & = orb_magmom_orb_decomposed( ilmta,: ) +ctmp(:)
                  End do
               End Do
            End do
         End do
      End Do
      deallocate( hsr_ssrep, hsi_ssrep )

    end subroutine contrib_hardpart

    subroutine contrib_hardpart_add_proj
      integer :: ia, it, lmt1, lmt2, is1, is2, istmp
      integer :: il1, il2, im1, im2, it1, it2, ilmta
      real(kind=DP), allocatable :: hsr_ssrep(:,:,:,:), hsi_ssrep(:,:,:,:)
      real(kind=DP) :: ctmp(3)
      complex(kind=CMPLDP) :: ztmp(3), z1, z2

      allocate( hsr_ssrep( natm, nlmt_add, nlmt_add, ndim_chgpot ) ); hsr_ssrep = 0.0d0
      allocate( hsi_ssrep( natm, nlmt_add, nlmt_add, ndim_chgpot ) ); hsi_ssrep = 0.0d0

      call m_ES_MagMom_To_DensMat_hsr( natm, nlmt_add, hsr_add, hsi_add, &
           &                           hsr_ssrep, hsi_ssrep )

      Do ia=1, natm
         it = ityp(ia)

         Do is1=1, ndim_spinor
            Do is2=1, ndim_spinor
               if ( is1 /= is2 ) cycle
               istmp = (is1 -1)*ndim_spinor +is2

               do lmt1=1,ilmt_add(it)
                  il1=ltp_add(lmt1,it); im1=mtp_add(lmt1,it); it1=1
                  ilmta = lmta_add(lmt1,ia)
                  do lmt2=1,ilmt_add(it)
                     il2=ltp_add(lmt2,it);  im2=mtp_add(lmt2,it); it2=1

                     if ( il1 <= 1 .or. il1 >4 ) cycle
                     if ( il1 /= il2 ) cycle

                     ztmp = 0.0d0
                     if ( il1 == 2 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L1( im2, im1, : )
                     else if ( il1 == 3 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L2( im2, im1, : )
                     else if ( il1 == 4 ) then
                        ztmp(:) = Mat_L_with_real_ylm_L3( im2, im1, : )
                     endif

                     z1 = dcmplx( hsr_ssrep(ia,lmt1,lmt2,istmp), &
                          &       hsi_ssrep(ia,lmt1,lmt2,istmp) )
                     ctmp(:) = z1*ztmp(:)

                     orb_magmom_orb_decomposed_add( ilmta,: ) &
                          & = orb_magmom_orb_decomposed_add( ilmta,: ) +ctmp(:)

                  End do
               End Do
            End do
         End do
      End Do
      deallocate( hsr_ssrep, hsi_ssrep )

    end subroutine contrib_hardpart_add_proj

    subroutine print_orb_decomposed_orb_mag_1
      integer :: ia, it, il1, im1, it1, lmt1, ilmta
      integer :: nspher
      real(kind=DP), allocatable :: orb_magmom_ylm(:,:)

      if ( mype /= 0 ) return

      allocate( orb_magmom_ylm(25,3) )

      write(nfout,'(A,5X,A)') &
           &         '    id      name     l    m', &
           &         '    mx             my             mz'

      Do ia=1, natm
         it = ityp(ia)
         orb_magmom_ylm = 0.0d0

         Do lmt1=1, ilmt(it)
            il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
            ilmta = lmta(lmt1,ia)

            nspher = ( il1 -1 )**2 +im1
            orb_magmom_ylm( nspher,: ) = orb_magmom_ylm( nspher,: ) &
                 &                     + orb_magmom_orb_decomposed( ilmta,: )
         End do
         if ( sw_use_add_proj == ON ) then
            Do lmt1=1, ilmt_add(it)
               il1=ltp_add(lmt1,it); im1=mtp_add(lmt1,it); it1=1
               ilmta = lmta_add(lmt1,ia)

               nspher = ( il1 -1 )**2 +im1
               orb_magmom_ylm( nspher,: ) = orb_magmom_ylm( nspher,: ) &
                    &                     + orb_magmom_orb_decomposed_add( ilmta,: )
            End do
         endif
         Do il1=1, lpsmax(it)
            Do im1=1, 2*il1 -1
               nspher = ( il1 -1 )**2 +im1
               write(nfout,'(I6,6X,A5,2I5,3F15.8)') &
                    &              ia, speciesname(it), il1 -1, im1, &
                    &              orb_magmom_ylm(nspher,1:3)
            End do
         End do
      End do

      deallocate( orb_magmom_ylm )

    end subroutine print_orb_decomposed_orb_mag_1

    subroutine print_orb_decomposed_orb_mag_2
      integer :: ia, it, il1, im1, it1, lmt1, ilmta
      integer :: nspher
      real(kind=DP) :: csum(3)
      real(kind=DP), allocatable :: orb_magmom_ylm(:,:)

      if ( mype /= 0 ) return

      allocate( orb_magmom_ylm(25,3) )

      write(nfout,'(A,5X,A)') &
           & '    id      name     l    ', &
           & '     mx             my             mz'

      Do ia=1, natm
         it = ityp(ia)
         orb_magmom_ylm = 0.0d0

         Do lmt1=1, ilmt(it)
            il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
            ilmta = lmta(lmt1,ia)

            nspher = ( il1 -1 )**2 +im1
            orb_magmom_ylm( nspher,: ) = orb_magmom_ylm( nspher,: ) &
                 &                     + orb_magmom_orb_decomposed( ilmta,: )
         End do
         if ( sw_use_add_proj == ON ) then
            Do lmt1=1, ilmt_add(it)
               il1=ltp_add(lmt1,it); im1=mtp_add(lmt1,it); it1=1
               ilmta = lmta_add(lmt1,ia)

               nspher = ( il1 -1 )**2 +im1
               orb_magmom_ylm( nspher,: ) = orb_magmom_ylm( nspher,: ) &
                    &                     + orb_magmom_orb_decomposed_add( ilmta,: )
            End do
         endif
         Do il1=1, lpsmax(it)
            csum = 0.0d0
            Do im1=1, 2*il1 -1
               nspher = ( il1 -1 )**2 +im1
               csum(:) = csum(:) +orb_magmom_ylm(nspher,:)
            End do
!
            write(nfout,'(I6,6X,A5,I5,5X,3F15.8)') &
                 &              ia, speciesname(it), il1 -1, &
                 &              csum(1:3)
         End do
      End do

      deallocate( orb_magmom_ylm )

    end subroutine print_orb_decomposed_orb_mag_2

    subroutine print_orb_decomposed_orb_mag_3
      integer :: ia, it, il1, im1, it1, lmt1, ilmta
      integer :: nspher
      real(kind=DP) :: csum(3)
      real(kind=DP), allocatable :: orb_magmom_ylm(:,:)

      if ( mype /= 0 ) return

      allocate( orb_magmom_ylm(25,3) )

      write(nfout,'(A,5X,A)') &
           & '    id      name          ', &
           & '     mx             my             mz'

      Do ia=1, natm
         it = ityp(ia)
         orb_magmom_ylm = 0.0d0

         Do lmt1=1, ilmt(it)
            il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
            ilmta = lmta(lmt1,ia)

            nspher = ( il1 -1 )**2 +im1
            orb_magmom_ylm( nspher,: ) = orb_magmom_ylm( nspher,: ) &
                 &                     + orb_magmom_orb_decomposed( ilmta,: )
         End do
         if ( sw_use_add_proj == ON ) then
            Do lmt1=1, ilmt_add(it)
               il1=ltp_add(lmt1,it); im1=mtp_add(lmt1,it); it1=1
               ilmta = lmta_add(lmt1,ia)

               nspher = ( il1 -1 )**2 +im1
               orb_magmom_ylm( nspher,: ) = orb_magmom_ylm( nspher,: ) &
                    &                     + orb_magmom_orb_decomposed_add( ilmta,: )
            End do
         endif

         csum = 0.0d0
         Do il1=1, lpsmax(it)
            Do im1=1, 2*il1 -1
               nspher = ( il1 -1 )**2 +im1
               csum(:) = csum(:) +orb_magmom_ylm(nspher,:)
            End do
         End do
!
         write(nfout,'(I6,6X,A5,10X,3F15.8)') &
              &              ia, speciesname(it), &
              &              csum(1:3)
      End do

      deallocate( orb_magmom_ylm )

    end subroutine print_orb_decomposed_orb_mag_3

  end subroutine m_OP_orb_decomposed_orb_magmom
! ==================== 2014/09/01

  subroutine m_OP_wd_orbmom_xsf
    character*64 file1
    integer :: lun = 4000

    integer :: i, iloop
    real(kind=DP) :: ctmp(3)

    Do iloop=1, 2
       if ( iloop == 1 ) then
          file1 = "orb_mag_moment.xsf"
       else
          file1 = "total_mag_moment.xsf"
       endif
       open( unit=lun, file=file1, status="unknown", form="formatted" )

       write(lun,'(A)') 'CRYSTAL'
       write(lun,'(A)') 'PRIMVEC'
       write(lun,'(3F20.10)') altv(1:3,1) *Bohr
       write(lun,'(3F20.10)') altv(1:3,2) *Bohr
       write(lun,'(3F20.10)') altv(1:3,3) *Bohr
       write(lun,'(A)') 'CONVVEC'
       write(lun,'(3F20.10)') altv(1:3,1) *Bohr
       write(lun,'(3F20.10)') altv(1:3,2) *Bohr
       write(lun,'(3F20.10)') altv(1:3,3) *Bohr
       write(lun,'(A)') 'PRIMCOORD'
       write(lun,*) natm, 1
       Do i=1, natm
          if ( iloop == 1 ) then
             ctmp(1:3) = OrbMag_on_atom(i,1:3)
          else
             ctmp(1:3) = RhoMag_on_atom(i,2:4) +OrbMag_on_atom(i,1:3)
          endif
          write(lun,'(I5,6F20.10)') nint(iatomn(ityp(i))), cps(i,1:3)*Bohr, ctmp(1:3)
       End Do
       close(lun)
    End Do

  end subroutine m_OP_wd_orbmom_xsf

end module m_OP_Moment
