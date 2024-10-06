module m_CD_Mag_Moment
! $Id: m_CD_Mag_Moment.f90 634 2020-12-22 12:17:00Z ktagami $
  use m_Control_Parameters,    only : noncol, ndim_magmom, kimg, iprimagmom, ON, OFF
  use m_Const_Parameters,     only : DP, PAI4, Bohr

  use m_Ionic_System,         only : ntyp, natm, ityp, pos, cps, &
       &                             magmom_local_now, iatomn, speciesname
  use m_Crystal_Structure,    only : rltv, altv

  use m_PlaneWaveBasisSet,    only : ngabc, kgp, ngabc_kngp_l
  use m_Parallelization,      only : ista_kngp, iend_kngp, ierr, npes, mype, &
       &                             MPI_CommGroup, mpi_chg_world
  use m_PseudoPotential,      only : dl2p, ilmt

  use m_Charge_Density,       only : chgq_l, hsr

  use m_Files,                only : nfout

  use m_Crystal_Structure,   only : sw_magnetic_constraint
  use m_Control_Parameters,  only : sw_modified_TFW_functional, sw_opencore, &
       &                            sw_calc_extfnv_correction

! ==== KT_add ==== 2014/08/29, 2016/09/14
  use m_Const_Parameters,   only : CMPLDP
  use m_Control_Parameters,  only : ndim_spinor, ndim_chgpot, sw_use_add_proj, &
       &                            sw_use_contracted_psir, nspin, proj_attribute, &
       &                            num_projectors
  use m_PseudoPotential,   only : nlmt, ltp, mtp, taup, nlmt_add, ltp_add, mtp_add, &
       &                          ilmt_add, mmesh, radr, wos, nmesh, &
       &                          psir_val, lpsmax, nloc, itau, ntau, &
       &                          nlmta, nlmta_add, lmta, lmta_add
  use m_ES_Noncollinear,  only :  m_ES_set_Pauli_Matrix, m_ES_MagMom_To_DensMat_hsr
  use m_Charge_Density,       only : hsi, hsr_add, hsi_add
  use m_Parallelization,    only : ista_k
  use m_NonLocal_Potential, only : new_radr_and_wos
  use m_Orbital_QuantumNum,  only :  max_num_orb_index, num_orb_index_data, &
       &                             qnum_n_orb_index, qnum_l_orb_index, &
       &                             qnum_tau_orb_index
! ================ 2014/08/29, 2016/09/14
  use mpi

  implicit none

!  include 'mpif.h'
!  integer istatus(mpi_status_size)

  integer :: sw_monitor_atomcharge = OFF

  real(kind=DP), allocatable :: rad_cov(:)
  real(kind=DP), allocatable :: RhoMag_on_atom(:,:)

  real(kind=DP) :: rad_cov_default(120)

! ==== KT_add ==== 2014/08/29
  real(kind=DP), allocatable :: product_psir_val(:,:,:,:)
! ================ 2014/08/29

contains

  subroutine m_CD_set_sw_monitor_atomcharge
    if ( noncol ) then
       sw_monitor_atomcharge = ON
    else
       if ( sw_magnetic_constraint == ON ) then
          sw_monitor_atomcharge = ON
       endif
    endif

    if ( sw_calc_extfnv_correction == ON ) sw_monitor_atomcharge = ON
    if ( sw_opencore == ON ) sw_monitor_atomcharge = ON
    if ( sw_modified_TFW_functional == ON ) sw_monitor_atomcharge = ON

    if ( npes > 1 ) then
       call mpi_bcast( iprimagmom, 1, mpi_integer, 0, MPI_CommGroup, ierr )
    endif
    if ( iprimagmom >= 2 ) sw_monitor_atomcharge = ON

  end subroutine m_CD_set_sw_monitor_atomcharge

  subroutine m_CD_set_proj_radius
    integer :: i, it

    do i = 1, num_projectors
       it = proj_attribute(i)%ityp
       if ( it <= 0 ) cycle

       proj_attribute(i)%radius = rad_cov_default( nint(iatomn(it)) )
       proj_attribute(i)%radius_was_defined = .true.

       if ( mype == 0 ) then
          write(nfout,'(A,I5,A,F15.8)') "** The radius of ", i, &
               &         "-th projector is set to ", proj_attribute(i)%radius
       endif
    end do
  end subroutine m_CD_set_proj_radius

  subroutine m_CD_alloc_rad_cov
    if ( .not. allocated( rad_cov) ) then
       allocate( rad_cov(natm) ); rad_cov = 0.0d0
    endif
  end subroutine m_CD_alloc_rad_cov

  subroutine m_CD_set_rad_cov_default
    integer :: ia, it, inum

! -------------------
!  Ref. "Molecular Single-Bond Covalent Radii for Elements 1-118".
!        Chemistry: A European Journal 15 (2009) 186.
! -------------------

    call set_table_rad_cov(   1,  32,  46, 133, 102,  85,  75,  71,  63,  64,  67 )
    call set_table_rad_cov(  11, 155, 139, 126, 116, 111, 103,  99,  96, 196, 171 )
    call set_table_rad_cov(  21, 148, 136, 134, 122, 119, 116, 111, 110, 112, 118 )
    call set_table_rad_cov(  31, 124, 124, 121, 116, 114, 117, 210, 185, 163, 154 )
    call set_table_rad_cov(  41, 147, 138, 128, 125, 125, 120, 128, 136, 142, 140 )
    call set_table_rad_cov(  51, 140, 136, 133, 131, 232, 196, 180, 163, 176, 174 )
    call set_table_rad_cov(  61, 173, 172, 168, 169, 168, 167, 166, 165, 164, 170 )
    call set_table_rad_cov(  71, 162, 152, 146, 137, 131, 129, 122, 123, 124, 133 )
    call set_table_rad_cov(  81, 144, 144, 151, 145, 147, 142, 223, 201, 186, 175 )
    call set_table_rad_cov(  91, 169, 170, 171, 172, 166, 166, 166, 168, 165, 167 )
    call set_table_rad_cov( 101, 173, 176, 161, 157, 149, 143, 141, 134, 129, 128 )
    call set_table_rad_cov( 111, 121, 122, 136, 143, 162, 175, 165, 157, 500, 500 )

    call convert_unit_to_bohr

  contains

    subroutine convert_unit_to_bohr
      rad_cov_default = rad_cov_default / 1.0D2 / Bohr
                                 ! Note : original data is in the unit of pm.
    end subroutine convert_unit_to_bohr

    subroutine set_table_rad_cov( istart, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10 )
      integer, intent(in) :: istart
      integer, intent(in) :: d1, d2, d3, d4, d5, d6, d7, d8, d9, d10

      rad_cov_default( istart    ) = dble(d1)
      rad_cov_default( istart +1 ) = dble(d2)
      rad_cov_default( istart +2 ) = dble(d3)
      rad_cov_default( istart +3 ) = dble(d4)
      rad_cov_default( istart +4 ) = dble(d5)
      rad_cov_default( istart +5 ) = dble(d6)
      rad_cov_default( istart +6 ) = dble(d7)
      rad_cov_default( istart +7 ) = dble(d8)
      rad_cov_default( istart +8 ) = dble(d9)
      rad_cov_default( istart +9 ) = dble(d10)
    end subroutine set_table_rad_cov

  end subroutine m_CD_set_rad_cov_default

  subroutine m_CD_alloc_RhoMag_on_atom
    if ( .not. allocated( RhoMag_on_atom) ) allocate( RhoMag_on_atom(natm,ndim_magmom) )
    RhoMag_on_atom = 0.0d0
  end subroutine m_CD_alloc_RhoMag_on_atom

  subroutine m_CD_set_rad_cov_now
    integer :: ia, ja, it1, it2, inum
    integer :: nx, ny, nz

    real(kind=DP) :: x1, y1, z1, dist
    real(kind=DP) :: c1, ctmp
    real(kind=DP), allocatable :: dist_min_among_atomtypes(:,:)

    logical, save :: First = .true.

    if ( First ) then
       Do ia=1, natm
          it1 = ityp(ia)
          inum = iatomn(it1)
          rad_cov(ia) = rad_cov_default(inum)
       End do
       First = .false.
    endif

    if ( .not. allocated(dist_min_among_atomtypes) ) then
       allocate( dist_min_among_atomtypes(ntyp,ntyp) )
    endif
    dist_min_among_atomtypes = 10.0d0

    Do ia=1, natm-1
       it1 = ityp(ia)
       Do ja =ia+1, natm
          it2 = ityp(ja)
!
          c1 = rad_cov(ia) + rad_cov(ja)

          Do nx=-1, 1
             Do ny=-1, 1
                Do nz=-1, 1
                   x1 = cps(ja,1) -cps(ia,1) +altv(1,1)*nx +altv(1,2)*ny +altv(1,3)*nz
                   y1 = cps(ja,2) -cps(ia,2) +altv(2,1)*nx +altv(2,2)*ny +altv(2,3)*nz
                   z1 = cps(ja,3) -cps(ia,3) +altv(3,1)*nx +altv(3,2)*ny +altv(3,3)*nz
!
                   dist = sqrt( x1**2 +y1**2 +z1**2 )
                   if ( dist < c1 ) then
                      ctmp = dist /c1
                      rad_cov(ia) = rad_cov(ia) *ctmp
                      rad_cov(ja) = rad_cov(ja) *ctmp
                   endif
                End do
             End do
          End do
       End do
    End Do

    if ( iprimagmom >= 3 ) then
       write(nfout,*) '! ------ info : Rad_cov ------ '
       write(nfout,*) '  atom no.   rad_cov '
       Do ia=1, natm
          write(nfout,'(I7,F15.8)')   ia, rad_cov(ia)
       End do
       write(nfout,*) '! ---------------------------- '
    endif

  end subroutine m_CD_set_rad_cov_now

  subroutine m_CD_calc_ChgMagMom_in_sphere
    integer :: i, ia, is, it, ist
    real(kind=DP) :: rad1, fac1r, fac1i, fac2
    real(kind=DP) :: VecG(3), normG, normG3, gr, d1

    real(kind=DP), allocatable :: zfcos(:), zfsin(:)
    real(kind=DP), allocatable :: RhoMag_on_atom_mpi(:,:)

    RhoMag_on_atom = 0.0d0

    allocate(zfcos(ista_kngp:iend_kngp)); zfcos = 0.d0
    allocate(zfsin(ista_kngp:iend_kngp)); zfsin = 0.d0

    Do ia=1, natm
       it = ityp(ia)
       rad1 = rad_cov(ia)

       call calc_phase2( natm, pos, ia, kgp, ngabc_kngp_l, ista_kngp, iend_kngp, &
            &            zfcos, zfsin )

       if ( ista_kngp == 1 ) then
          ist = 2
       else
          ist = ista_kngp
       endif

       Do i=ist, iend_kngp

          VecG(1) = rltv(1,1)*ngabc_kngp_l(i,1) +rltv(1,2)*ngabc_kngp_l(i,2) +rltv(1,3)*ngabc_kngp_l(i,3)
          VecG(2) = rltv(2,1)*ngabc_kngp_l(i,1) +rltv(2,2)*ngabc_kngp_l(i,2) +rltv(2,3)*ngabc_kngp_l(i,3)
          VecG(3) = rltv(3,1)*ngabc_kngp_l(i,1) +rltv(3,2)*ngabc_kngp_l(i,2) +rltv(3,3)*ngabc_kngp_l(i,3)

          normG = sqrt( VecG(1)**2 +VecG(2)**2 +VecG(3)**2 )
          normG3 = normG**3

          fac1r = zfcos(i); fac1i = zfsin(i)
!
          gr = normG *rad1
          fac2 = -gr *cos(gr) + sin(gr)
          fac2 = fac2 *PAI4 /normG3
!
          Do is=1, ndim_magmom
             if ( kimg == 1 ) then
                d1 = chgq_l(i,1,is)*fac1r
             else
                d1 = chgq_l(i,1,is)*fac1r -chgq_l(i,kimg,is)*fac1i
             endif
             RhoMag_on_atom(ia,is) = RhoMag_on_atom(ia,is) + d1 *fac2
          End do

       End do

       if ( mype == 0 ) then
          fac2 = PAI4 /3.0d0 *rad1**3

          Do is=1, ndim_magmom
             d1 = chgq_l(1,1,is)
             RhoMag_on_atom(ia,is) = RhoMag_on_atom(ia,is) + d1 *fac2
          End do
       endif

    End Do
!
    deallocate( zfcos, zfsin )

    if ( npes > 1 ) then
       allocate( RhoMag_on_atom_mpi(natm,ndim_magmom) ); RhoMag_on_atom_mpi = 0.0d0
       call mpi_allreduce( RhoMag_on_atom, RhoMag_on_atom_mpi, natm*ndim_magmom, &
            &              mpi_double_precision, mpi_sum, mpi_chg_world, ierr )
       RhoMag_on_atom = RhoMag_on_atom_mpi
       deallocate( RhoMag_on_atom_mpi )
    endif

  end subroutine m_CD_calc_ChgMagMom_in_sphere

  subroutine m_CD_print_ChgMagmom_on_atom( nfout )
    integer, intent(in) :: nfout

    integer :: ia, it

    if ( iprimagmom < 2 ) return
    if ( mype /= 0 ) return

    write(nfout,*) &
         & '! ------------ Local Charge/Spin Moment (in sphere) at this scf step --- '

    if ( noncol ) then
       write(nfout,*) &
            & '    id      name          tot            mx             my             mz'
       Do ia=1, natm
          it = ityp(ia)
          write(nfout,'(I7,6X,A5,4F15.8)') ia, speciesname(it), &
               &                          RhoMag_on_atom(ia,1), &
               &                          RhoMag_on_atom(ia,2), &
               &                          RhoMag_on_atom(ia,3), &
               &                          RhoMag_on_atom(ia,4)
       End Do

    else
       if ( ndim_magmom == 1 ) then
          write(nfout,*) &
               & '    id      name          tot'
          Do ia=1, natm
             it = ityp(ia)
             write(nfout,'(I7,6X,A5,F15.8)') ia, speciesname(it), &
                  &                          RhoMag_on_atom(ia,1)
          End Do

       else if ( ndim_magmom == 2 ) then
          write(nfout,*) &
               & '    id      name          tot            mz'
          Do ia=1, natm
             it = ityp(ia)
             write(nfout,'(I7,6X,A5,2F15.8)') ia, speciesname(it), &
                  &                          RhoMag_on_atom(ia,1)+RhoMag_on_atom(ia,2),&
                  &                          RhoMag_on_atom(ia,1)-RhoMag_on_atom(ia,2)
          End Do
       endif

    endif

!    write(nfout,*) '! ---------------------------------------------'
    write(nfout,*)

  end subroutine m_CD_print_ChgMagmom_on_atom

! ------------------------------------------
  subroutine m_CD_estim_magmom_local( nfout )
    integer, intent(in) :: nfout

    integer :: ia, it, lmt1, lmt2, is
    real(kind=DP) :: sum, fac

!--
    Do ia=1, natm
       it = ityp(ia)

       Do is=2, ndim_magmom
          sum = 0.0d0
          Do lmt1=1, ilmt(it)
             Do lmt2=lmt1, ilmt(it)
                fac=2.d0;if(lmt1.eq.lmt2) fac=1.d0
                sum = sum + fac *hsr( ia, lmt1, lmt2, is )*dl2p(lmt1,lmt2,1,it)
             End do
          End do
!
          magmom_local_now(ia,is-1) = sum
       End do
    End Do

  end subroutine m_CD_estim_magmom_local

! === KT_add === 2014/08/25
  subroutine m_CD_calc_SpinMagMom_method2
    complex(kind=CMPLDP) :: PauliMatrix( ndim_magmom, ndim_spinor, ndim_spinor )

    if ( iprimagmom < 3 ) return

    call m_ES_set_Pauli_Matrix( PauliMatrix )

    call contrib_hardpart
    if ( sw_use_add_proj == ON ) call contrib_hardpart_add_proj

    call print_magmom

  contains

    subroutine contrib_hardpart
      integer :: ia, it, lmt1, lmt2, is1, is2, istmp
      integer :: il1, il2, im1, im2, it1, it2
      real(kind=DP) :: sum, fac, spinmom(3)
!
      complex(kind=CMPLDP) :: z1
!
      real(kind=DP), allocatable, dimension(:,:,:,:) :: hsr_ssrep
      real(kind=DP), allocatable, dimension(:,:,:,:) :: hsi_ssrep
!
      allocate( hsr_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsr_ssrep = 0.0d0
      allocate( hsi_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsi_ssrep = 0.0d0

      call m_ES_MagMom_To_DensMat_hsr( natm, nlmt, hsr, hsi, hsr_ssrep, hsi_ssrep )
!--
      Do ia=1, natm
         it = ityp(ia)

         spinmom = 0.0d0

         Do is1=1, ndim_spinor
            Do is2=1, ndim_spinor
               istmp = ( is1 -1 )*ndim_spinor +is2

               do lmt1=1,ilmt(it)
                  il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
                  do lmt2=1,ilmt(it)
                     il2=ltp(lmt2,it);  im2=mtp(lmt2,it); it2=taup(lmt2,it)

                     if ( il1 /= il2 ) cycle
                     if ( im1 /= im2 ) cycle

                     z1 = dcmplx( hsr_ssrep(ia,lmt1,lmt2,istmp), &
                          &       hsi_ssrep(ia,lmt1,lmt2,istmp) )
                     if ( sw_use_contracted_psir == ON ) then
                        z1 = z1 *product_psir_val( il1, it1, it2, it )
                     endif

                     spinmom(1) = spinmom(1) + z1 *PauliMatrix(2,is2,is1)
                     spinmom(2) = spinmom(2) + z1 *PauliMatrix(3,is2,is1)
                     spinmom(3) = spinmom(3) + z1 *PauliMatrix(4,is2,is1)
                  End do
               End Do
            End do
         End do
         magmom_local_now(ia,1) = spinmom(1)
         magmom_local_now(ia,2) = spinmom(2)
         magmom_local_now(ia,3) = spinmom(3)
      End Do

      deallocate( hsr_ssrep ); deallocate( hsi_ssrep )
    end subroutine contrib_hardpart

    subroutine contrib_hardpart_add_proj
      integer :: ia, it, lmt1, lmt2, is1, is2, istmp
      integer :: il1, il2, im1, im2, it1, it2
      real(kind=DP) :: sum, fac, spinmom(3)
!
      complex(kind=CMPLDP) :: z1
!
      real(kind=DP), allocatable, dimension(:,:,:,:) :: hsr_ssrep
      real(kind=DP), allocatable, dimension(:,:,:,:) :: hsi_ssrep
!
      allocate( hsr_ssrep( natm, nlmt_add, nlmt_add, ndim_chgpot ) ); hsr_ssrep = 0.0d0
      allocate( hsi_ssrep( natm, nlmt_add, nlmt_add, ndim_chgpot ) ); hsi_ssrep = 0.0d0

      call m_ES_MagMom_To_DensMat_hsr( natm, nlmt_add, hsr_add, hsi_add, &
           &                           hsr_ssrep, hsi_ssrep )
!--
      Do ia=1, natm
         it = ityp(ia)

         spinmom = 0.0d0

         Do is1=1, ndim_spinor
            Do is2=1, ndim_spinor
               istmp = ( is1 -1 )*ndim_spinor +is2

               do lmt1=1,ilmt_add(it)
                  il1=ltp_add(lmt1,it); im1=mtp_add(lmt1,it); it1=1
                  do lmt2=1,ilmt_add(it)
                     il2=ltp_add(lmt2,it);  im2=mtp_add(lmt2,it); it2=1

                     if ( il1 /= il2 ) cycle
                     if ( im1 /= im2 ) cycle

                     z1 = dcmplx( hsr_ssrep(ia,lmt1,lmt2,istmp), &
                          &       hsi_ssrep(ia,lmt1,lmt2,istmp) )

                     spinmom(1) = spinmom(1) + z1 *PauliMatrix(2,is2,is1)
                     spinmom(2) = spinmom(2) + z1 *PauliMatrix(3,is2,is1)
                     spinmom(3) = spinmom(3) + z1 *PauliMatrix(4,is2,is1)
                  End do
               End Do
            End do
         End do
         magmom_local_now(ia,1) = magmom_local_now(ia,1) +spinmom(1)
         magmom_local_now(ia,2) = magmom_local_now(ia,2) +spinmom(2)
         magmom_local_now(ia,3) = magmom_local_now(ia,3) +spinmom(3)
      End Do

      deallocate( hsr_ssrep ); deallocate( hsi_ssrep )
    end subroutine contrib_hardpart_add_proj

    subroutine print_magmom
      integer :: ia, it

      write(nfout,*) &
           & '! ------------ Local Spin Moment (method2) at this scf step --- '
      write(nfout,*) &
           & '    id      name          mx             my             mz'
      Do ia=1, natm
         it = ityp(ia)
         write(nfout,'(I7,6X,A5,3F15.8)') ia, speciesname(it), &
              &                         magmom_local_now(ia,1), &
              &                         magmom_local_now(ia,2), &
              &                         magmom_local_now(ia,3)
      End Do
      write(nfout,*)

    end subroutine print_magmom

  end subroutine m_CD_calc_SpinMagMom_method2
! =============== 2014/08/25

  subroutine m_CD_print_magmom_local( nfout )
    integer, intent(in) :: nfout

    integer :: ia, it

    if ( iprimagmom < 3 ) return
    if ( mype /= 0 ) return

    write(nfout,*) &
         & '! ------------ Local Spin Moment (hard part) at this scf step --- '
    write(nfout,*) &
         & '    id      name          mx             my             mz'
    Do ia=1, natm
       it = ityp(ia)
       write(nfout,'(I7,6X,A5,3F15.8)') ia, speciesname(it), &
            &                         magmom_local_now(ia,1), &
            &                         magmom_local_now(ia,2), &
            &                         magmom_local_now(ia,3)
    End Do
    write(nfout,*)

  end subroutine m_CD_print_magmom_local

! ===== KT_add ===== 2014/08/29
  subroutine m_CD_calc_product_psir_val
    integer :: it, ik, ilmt1, il1, tau1, tau2, ir
    real(kind=DP) :: csum, factor, rcut
    real(kind=DP), allocatable :: tmp_fn(:)

    allocate( tmp_fn(mmesh) );  allocate( radr(mmesh) );  allocate( wos(mmesh) )

    Do it=1, ntyp
       ik = ista_k
       call new_radr_and_wos(ik,it)                 ! --> radr, wos

       rcut = rad_cov_default( nint(iatomn(it) ))  ! Revised according to a report from ASMS Co.ltd, 10 March 2016.

       Do il1=1, lpsmax(it)
!          if ( il1 == iloc(it) ) cycle

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
! ---
       Do il1=1, lpsmax(it)
          Do tau1=1, itau(il1,it)
             Do tau2=1, itau(il1,it)
                csum = 0.0d0
                Do ir=1, nmesh(it)
                   if ( radr(ir) > rcut ) exit
                   csum = csum +wos(ir) *psir_val(ir,il1,tau1,it) &
                        &               *psir_val(ir,il1,tau2,it)
                End Do
                product_psir_val( il1, tau1, tau2, it ) = csum
             End do
          End do
       End Do
    End Do
!
    if ( mype == 0 .and. iprimagmom>2 ) then
       write(nfout,*) '!** product_psir_val'
       write(nfout,*) "   it     l    tau1  tau2     product"
       DO it=1, ntyp
          Do il1=1, lpsmax(it)
             Do tau1=1, itau(il1,it)
                Do tau2=1, itau(il1,it)
                   write(nfout,'(4I6,F15.8)') it, il1, tau1, tau2, &
                        &         product_psir_val( il1, tau1, tau2, it )
                End do
             End do
          End do
       End do
    end if

    call mpi_barrier( MPI_CommGroup, ierr )
    deallocate(tmp_fn); deallocate(radr);  deallocate(wos)

  end subroutine m_CD_calc_product_psir_val

  subroutine m_CD_alloc_product_psir_val
    allocate( product_psir_val(nloc,ntau,ntau,ntyp))
    product_psir_val = 0.0d0
  end subroutine m_CD_alloc_product_psir_val

  subroutine m_CD_dealloc_product_psir_val
    if ( allocated( product_psir_val ) )  deallocate( product_psir_val )
  end subroutine m_CD_dealloc_product_psir_val
! ============ 2014/08/29

! ====== KT_add === 2014/09/01
  subroutine m_CD_orb_decomposed_chg_magmom
    complex(kind=CMPLDP) :: PauliMatrix( ndim_magmom, ndim_spinor, ndim_spinor )

    real(kind=DP), allocatable :: chg_magmom_orb_decomposed(:,:)
    real(kind=DP), allocatable :: chg_magmom_orb_decomposed_add(:,:)

    call m_ES_set_Pauli_Matrix( PauliMatrix )

    allocate( chg_magmom_orb_decomposed( nlmta, ndim_magmom ) )
    chg_magmom_orb_decomposed = 0.0d0
    if ( sw_use_add_proj == ON ) then
       allocate( chg_magmom_orb_decomposed_add( nlmta_add, ndim_magmom ) )
       chg_magmom_orb_decomposed_add = 0.0d0
    endif

    call contrib_hardpart
    if ( sw_use_add_proj == ON ) call contrib_hardpart_add_proj

    write(nfout,*)
    write(nfout,*) &
         &  "! --------- Orbital decomposed charge/spin magnetic oment (method2) ------"
    write(nfout,*) "< Detail >"
    call print_orb_decomposed_chg_mag_1

    write(nfout,*) "< Summed over m >"
    call print_orb_decomposed_chg_mag_2

    write(nfout,*) "< Summed over l,m >"
    call print_orb_decomposed_chg_mag_3

!    if ( max_num_orb_index > 0 ) then
!       write(nfout,*) "< Classified by n,l >"
!       call print_orb_decomposed_chg_mag_4
!    endif

    write(nfout,*)

    deallocate( chg_magmom_orb_decomposed )
    if ( sw_use_add_proj == ON ) then
       deallocate( chg_magmom_orb_decomposed_add )
    endif

  contains

    subroutine contrib_hardpart
      integer :: ia, it, lmt1, lmt2, is1, is2, istmp
      integer :: il1, il2, im1, im2, it1, it2, ilmta
      real(kind=DP) :: sum, fac, ctmp(ndim_magmom)
!
      real(kind=DP) :: c1
      complex(kind=CMPLDP) :: z1
!
      real(kind=DP), allocatable, dimension(:,:,:,:) :: hsr_ssrep
      real(kind=DP), allocatable, dimension(:,:,:,:) :: hsi_ssrep
!
      if ( noncol ) then
         allocate( hsr_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsr_ssrep = 0.0d0
         allocate( hsi_ssrep( natm, nlmt, nlmt, ndim_chgpot ) ); hsi_ssrep = 0.0d0

         call m_ES_MagMom_To_DensMat_hsr( natm, nlmt, hsr, hsi, hsr_ssrep, hsi_ssrep )

         Do ia=1, natm
            it = ityp(ia)

            Do is1=1, ndim_spinor
               Do is2=1, ndim_spinor
                  istmp = ( is1 -1 )*ndim_spinor +is2

                  do lmt1=1,ilmt(it)
                     il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
                     ilmta = lmta(lmt1,ia)

                     do lmt2=1,ilmt(it)
                        il2=ltp(lmt2,it);  im2=mtp(lmt2,it); it2=taup(lmt2,it)

                        if ( il1 /= il2 ) cycle
                        if ( im1 /= im2 ) cycle

                        z1 = dcmplx( hsr_ssrep(ia,lmt1,lmt2,istmp), &
                             &       hsi_ssrep(ia,lmt1,lmt2,istmp) )
                        if ( sw_use_contracted_psir == ON ) then
                           z1 = z1 *product_psir_val( il1, it1, it2, it )
                        endif

                        ctmp(:) = z1 *PauliMatrix(:,is2,is1)

                        chg_magmom_orb_decomposed( ilmta,: ) &
                             & = chg_magmom_orb_decomposed( ilmta,: ) +ctmp(:)
                     End do
                  End Do
               End do
            End do
         End Do
         deallocate( hsr_ssrep ); deallocate( hsi_ssrep )

      else

         Do ia=1, natm
            it = ityp(ia)
            Do is1=1, nspin
               do lmt1=1,ilmt(it)
                  il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
                  ilmta = lmta(lmt1,ia)

                  do lmt2=1,ilmt(it)
                     il2=ltp(lmt2,it);  im2=mtp(lmt2,it); it2=taup(lmt2,it)

                     if ( il1 /= il2 ) cycle
                     if ( im1 /= im2 ) cycle

                     c1 = hsr(ia,lmt1,lmt2,is1)
                     if ( sw_use_contracted_psir == ON ) then
                        c1 = c1 *product_psir_val( il1, it1, it2, it )
                     endif

                     chg_magmom_orb_decomposed( ilmta,is1 ) &
                          & = chg_magmom_orb_decomposed( ilmta,is1 ) +c1

                  End Do
               End do
            End do
         End Do
      endif
    end subroutine contrib_hardpart

    subroutine contrib_hardpart_add_proj
      integer :: ia, it, lmt1, lmt2, is1, is2, istmp
      integer :: il1, il2, im1, im2, it1, it2, ilmta
      real(kind=DP) :: sum, fac, ctmp(ndim_magmom)
!
      real(kind=DP) :: c1
      complex(kind=CMPLDP) :: z1
!
      real(kind=DP), allocatable, dimension(:,:,:,:) :: hsr_ssrep
      real(kind=DP), allocatable, dimension(:,:,:,:) :: hsi_ssrep
!
      if ( noncol ) then
         allocate( hsr_ssrep( natm, nlmt_add, nlmt_add, ndim_chgpot ) );
         allocate( hsi_ssrep( natm, nlmt_add, nlmt_add, ndim_chgpot ) );
         hsr_ssrep = 0.0d0; hsi_ssrep = 0.0d0

         call m_ES_MagMom_To_DensMat_hsr( natm, nlmt_add, hsr_add, hsi_add, &
              &                           hsr_ssrep, hsi_ssrep )
!--
         Do ia=1, natm
            it = ityp(ia)

            Do is1=1, ndim_spinor
               Do is2=1, ndim_spinor
                  istmp = ( is1 -1 )*ndim_spinor +is2

                  do lmt1=1,ilmt_add(it)
                     il1=ltp_add(lmt1,it); im1=mtp_add(lmt1,it); it1=1
                     ilmta = lmta_add(lmt1,ia)
                     do lmt2=1,ilmt_add(it)
                        il2=ltp_add(lmt2,it);  im2=mtp_add(lmt2,it); it2=1

                        if ( il1 /= il2 ) cycle
                        if ( im1 /= im2 ) cycle

                        z1 = dcmplx( hsr_ssrep(ia,lmt1,lmt2,istmp), &
                             &       hsi_ssrep(ia,lmt1,lmt2,istmp) )

                        ctmp(:) = z1 *PauliMatrix(:,is2,is1)

                        chg_magmom_orb_decomposed_add( ilmta,: ) &
                             & = chg_magmom_orb_decomposed_add( ilmta,: ) +ctmp(:)
                     End do
                  End Do
               End do
            End do
         End Do
         deallocate( hsr_ssrep ); deallocate( hsi_ssrep )

      else

         Do ia=1, natm
            it = ityp(ia)
            Do is1=1, nspin
               do lmt1=1,ilmt_add(it)
                  il1=ltp_add(lmt1,it); im1=mtp_add(lmt1,it); it1=1
                  ilmta = lmta_add(lmt1,ia)
                  do lmt2=1,ilmt_add(it)
                     il2=ltp_add(lmt2,it);  im2=mtp_add(lmt2,it); it2=1

                     if ( il1 /= il2 ) cycle
                     if ( im1 /= im2 ) cycle

                     c1 = hsr_add(ia,lmt1,lmt2,is1)
                     chg_magmom_orb_decomposed_add( ilmta,is1 ) &
                          & = chg_magmom_orb_decomposed_add( ilmta,is1 ) +c1
                  End Do
               End do
            End do
         End Do
      endif
    end subroutine contrib_hardpart_add_proj

    subroutine print_orb_decomposed_chg_mag_1
      integer :: ia, it, il1, im1, it1, lmt1, ilmta
      integer :: nspher
      real(kind=DP), allocatable :: chg_magmom_ylm(:,:)

      if ( mype /= 0 ) return

      allocate( chg_magmom_ylm(25,ndim_magmom) )

      if ( noncol ) then
         write(nfout,'(A,5X,A)') &
              &         '    id      name     l    m', &
              &         '    tot            mx             my             mz'
      else
         if ( nspin == 1 ) then
            write(nfout,'(A,5X,A)') &
                 &         '    id      name     l    m', '    tot'
         else
            write(nfout,'(A,5X,A)') &
                 &         '    id      name     l    m', '    up            down'
         endif
      endif

      Do ia=1, natm
         it = ityp(ia)
         chg_magmom_ylm = 0.0d0

         Do lmt1=1, ilmt(it)
            il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
            ilmta = lmta(lmt1,ia)

            nspher = ( il1 -1 )**2 +im1
            chg_magmom_ylm( nspher,: ) = chg_magmom_ylm( nspher,: ) &
                 &                     + chg_magmom_orb_decomposed( ilmta,: )
         End do
         if ( sw_use_add_proj == ON ) then
            Do lmt1=1, ilmt_add(it)
               il1=ltp_add(lmt1,it); im1=mtp_add(lmt1,it); it1=1
               ilmta = lmta_add(lmt1,ia)

               nspher = ( il1 -1 )**2 +im1
               chg_magmom_ylm( nspher,: ) = chg_magmom_ylm( nspher,: ) &
                    &                     + chg_magmom_orb_decomposed_add( ilmta,: )
            End do
         endif
         Do il1=1, lpsmax(it)
            Do im1=1, 2*il1 -1
               nspher = ( il1 -1 )**2 +im1
               if ( noncol ) then
                  write(nfout,'(I6,6X,A5,2I5,4F15.8)') &
                       &              ia, speciesname(it), il1 -1, im1, &
                       &              chg_magmom_ylm(nspher,1:ndim_magmom)
               else
                  if ( nspin == 1 ) then
                     write(nfout,'(I6,6X,A5,2I5,F15.8)') &
                          &              ia, speciesname(it), il1 -1, im1, &
                          &              chg_magmom_ylm(nspher,1:ndim_magmom)
                  else
                     write(nfout,'(I6,6X,A5,2I5,2F15.8)') &
                          &              ia, speciesname(it), il1 -1, im1, &
                          &              chg_magmom_ylm(nspher,1:ndim_magmom)
                  endif
               endif
            End do
         End do
      End do

      deallocate( chg_magmom_ylm )

    end subroutine print_orb_decomposed_chg_mag_1

    subroutine print_orb_decomposed_chg_mag_2
      integer :: ia, it, il1, im1, it1, lmt1, ilmta
      integer :: nspher
      real(kind=DP) :: csum( ndim_magmom )
      real(kind=DP), allocatable :: chg_magmom_ylm(:,:)

      if ( mype /= 0 ) return

      allocate( chg_magmom_ylm(25,ndim_magmom) )

      if ( noncol ) then
         write(nfout,'(A,5X,A)') &
              & '    id      name     l    ', &
              & '     tot            mx             my             mz'
      else
         if ( nspin == 1 ) then
            write(nfout,'(A,5X,A)') &
                 & '    id      name     l    ', '     tot'
         else
            write(nfout,'(A,5X,A)') &
                 & '    id      name     l    ', '     up            down'
         endif
      endif

      Do ia=1, natm
         it = ityp(ia)
         chg_magmom_ylm = 0.0d0

         Do lmt1=1, ilmt(it)
            il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
            ilmta = lmta(lmt1,ia)

            nspher = ( il1 -1 )**2 +im1
            chg_magmom_ylm( nspher,: ) = chg_magmom_ylm( nspher,: ) &
                 &                     + chg_magmom_orb_decomposed( ilmta,: )
         End do
         if ( sw_use_add_proj == ON ) then
            Do lmt1=1, ilmt_add(it)
               il1=ltp_add(lmt1,it); im1=mtp_add(lmt1,it); it1=1
               ilmta = lmta_add(lmt1,ia)

               nspher = ( il1 -1 )**2 +im1
               chg_magmom_ylm( nspher,: ) = chg_magmom_ylm( nspher,: ) &
                    &                     + chg_magmom_orb_decomposed_add( ilmta,: )
            End do
         endif
         Do il1=1, lpsmax(it)
            csum = 0.0d0
            Do im1=1, 2*il1 -1
               nspher = ( il1 -1 )**2 +im1
               csum(:) = csum(:) +chg_magmom_ylm(nspher,:)
            End do
!
            if ( noncol ) then
               write(nfout,'(I6,6X,A5,I5,5X,4F15.8)') &
                    &              ia, speciesname(it), il1 -1, &
                    &              csum(1:ndim_magmom)
            else
               if ( nspin == 1 ) then
                  write(nfout,'(I6,6X,A5,I5,5X,F15.8)') &
                       &              ia, speciesname(it), il1 -1, &
                       &              csum(1:ndim_magmom)
               else
                  write(nfout,'(I6,6X,A5,I5,5X,2F15.8)') &
                       &              ia, speciesname(it), il1 -1, &
                       &              csum(1:ndim_magmom)
               endif
            endif
         End do
      End do

      deallocate( chg_magmom_ylm )

    end subroutine print_orb_decomposed_chg_mag_2

    subroutine print_orb_decomposed_chg_mag_3
      integer :: ia, it, il1, im1, it1, lmt1, ilmta
      integer :: nspher
      real(kind=DP) :: csum( ndim_magmom )
      real(kind=DP), allocatable :: chg_magmom_ylm(:,:)

      if ( mype /= 0 ) return

      allocate( chg_magmom_ylm(25,ndim_magmom) )

      if ( noncol ) then
         write(nfout,'(A,5X,A)') &
              & '    id      name          ', &
              & '     tot            mx             my             mz'
      else
         if ( nspin == 1 ) then
            write(nfout,'(A,5X,A)') &
                 & '    id      name          ', '     tot'
         else
            write(nfout,'(A,5X,A)') &
                 & '    id      name          ', '     up            down'
         endif
      endif

      Do ia=1, natm
         it = ityp(ia)
         chg_magmom_ylm = 0.0d0

         Do lmt1=1, ilmt(it)
            il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
            ilmta = lmta(lmt1,ia)

            nspher = ( il1 -1 )**2 +im1
            chg_magmom_ylm( nspher,: ) = chg_magmom_ylm( nspher,: ) &
                 &                     + chg_magmom_orb_decomposed( ilmta,: )
         End do
         if ( sw_use_add_proj == ON ) then
            Do lmt1=1, ilmt_add(it)
               il1=ltp_add(lmt1,it); im1=mtp_add(lmt1,it); it1=1
               ilmta = lmta_add(lmt1,ia)

               nspher = ( il1 -1 )**2 +im1
               chg_magmom_ylm( nspher,: ) = chg_magmom_ylm( nspher,: ) &
                    &                     + chg_magmom_orb_decomposed_add( ilmta,: )
            End do
         endif

         csum = 0.0d0
         Do il1=1, lpsmax(it)
            Do im1=1, 2*il1 -1
               nspher = ( il1 -1 )**2 +im1
               csum(:) = csum(:) +chg_magmom_ylm(nspher,:)
            End do
         End do
!
         if ( noncol ) then
            write(nfout,'(I6,6X,A5,10X,4F15.8)') &
                 &              ia, speciesname(it),  &
                 &              csum(1:ndim_magmom)
         else
            if ( nspin == 1 ) then
               write(nfout,'(I6,6X,A5,10X,F15.8)') &
                    &              ia, speciesname(it), &
                    &              csum(1:ndim_magmom)
            else
               write(nfout,'(I6,6X,A5,10X,2F15.8)') &
                    &              ia, speciesname(it), &
                    &              csum(1:ndim_magmom)
            endif
         endif
      End do

      deallocate( chg_magmom_ylm )

    end subroutine print_orb_decomposed_chg_mag_3

    subroutine print_orb_decomposed_chg_mag_4
      integer :: ia, it, il1, im1, it1, lmt1, ilmta
      integer :: j1, n1
      real(kind=DP) :: csum( ndim_magmom )
      real(kind=DP), allocatable :: chg_magmom_nl(:,:,:)

      if ( mype /= 0 ) return

      allocate( chg_magmom_nl(7,4,ndim_magmom) )

      if ( noncol ) then
         write(nfout,'(A,5X,A)') &
              & '    id      name     n    l', &
              & '     tot            mx             my             mz'
      else
         if ( nspin == 1 ) then
            write(nfout,'(A,5X,A)') &
                 & '    id      name     n    l', '     tot'
         else
            write(nfout,'(A,5X,A)') &
                 & '    id      name     n    l', '     up            down'
         endif
      endif

      Do ia=1, natm
         it = ityp(ia)
         if ( num_orb_index_data(it) == 0 ) cycle

         chg_magmom_nl = 0.0d0

         Do lmt1=1, ilmt(it)
            il1=ltp(lmt1,it); im1=mtp(lmt1,it); it1=taup(lmt1,it)
            ilmta = lmta(lmt1,ia)

            Do j1=1, num_orb_index_data(it)
               if ( il1 == qnum_l_orb_index(it,j1) +1 &
                    &  .and. it1 == qnum_tau_orb_index(it,j1) ) then
                  n1 = qnum_n_orb_index(it,j1)
                  exit
               endif
            End do
            chg_magmom_nl( n1,il1,: ) = chg_magmom_nl( n1,il1,: ) &
                 &                     + chg_magmom_orb_decomposed( ilmta,: )
         End do
         if ( sw_use_add_proj == ON ) then
            Do lmt1=1, ilmt_add(it)
               il1=ltp_add(lmt1,it); im1=mtp_add(lmt1,it); it1=1
               ilmta = lmta_add(lmt1,ia)

               Do j1=1, num_orb_index_data(it)
                  if ( il1 == qnum_l_orb_index(it,j1) +1 &
                       &  .and. it1 == qnum_tau_orb_index(it,j1) ) then
                     n1 = qnum_n_orb_index(it,j1)
                     exit
                  endif
               End do
               chg_magmom_nl( n1,il1,: ) = chg_magmom_nl( n1,il1,: ) &
                 &                     + chg_magmom_orb_decomposed( ilmta,: )
            End do
         endif

         Do j1=1, num_orb_index_data(it)
            n1  = qnum_n_orb_index(it,j1)
            il1 = qnum_l_orb_index(it,j1) +1

            if ( noncol ) then
               write(nfout,'(I6,6X,A5,2I5,4F15.8)') &
                    &              ia, speciesname(it), n1,il1 -1, &
                    &              chg_magmom_nl(n1,il1,1:ndim_magmom)
            else
               if ( nspin == 1 ) then
                  write(nfout,'(I6,6X,A5,2I5,F15.8)') &
                       &              ia, speciesname(it), n1, il1 -1, &
                       &              chg_magmom_nl(n1,il1,1:ndim_magmom)
               else
                  write(nfout,'(I6,6X,A5,2I5,2F15.8)') &
                       &              ia, speciesname(it), n1, il1 -1, &
                       &              chg_magmom_nl(n1,il1,1:ndim_magmom)
               endif
            endif
         End do
      End do

      deallocate( chg_magmom_nl )

    end subroutine print_orb_decomposed_chg_mag_4

  end subroutine m_CD_orb_decomposed_chg_magmom
! ======== 2014/09/01

  subroutine m_CD_wd_Magmom_xsf
    character*64 file1
    integer :: lun = 4000

    integer :: i

    file1 = "spin_mag_moment.xsf"
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
       write(lun,'(I5,6F20.10)') nint(iatomn(ityp(i))), cps(i,1:3)*Bohr, &
            &                    RhoMag_on_atom(i,2:4)
    End Do
    close(lun)

  end subroutine m_CD_wd_Magmom_xsf

end module m_CD_Mag_Moment
