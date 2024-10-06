#define MOMENT_AS_PSEUDO_VECTOR

module m_CS_Magnetic
! $Id: m_CS_Magnetic.F90 250 2012-11-16 05:29:11Z ktagami $
  use m_Const_Parameters,   only : DP, BUCS, CARTS, CRDTYP, DELTA10, &
       &                           oh=>oh_symbol, d6h=>d6h_symbol, ON, DELTA07, &
       &                           CMPLDP, zi, PAI, OFF, PAI2

  use m_Crystal_Structure,      only : nopr, op, ig01, tau, m_CS_op_in_PUCD, &
       &                               pg_symbol_system, ngen, igen, jgen, &
       &                               alloc_igen_jgen, dealloc_igen_jgen, &
       &                               sw_reduce_sym_by_magmom, &
       &                               sw_reduce_sym_by_orbital, sw_neglect_magmom, &
       &                               sw_allow_mag_sym_inversion, &
       &                               sw_allow_improper_rotation

  use m_Ionic_System,       only : natm, ntyp, ityp, mag_direction0_atomtyp, &
       &                           natm2, pos, iwei, magmom_local_now, &
       &                           has_partially_filled_lcore

  use m_Files,                 only : nfout

! == KT_add === 2014/08/14
  use m_Ionic_System,  only :   lattice_system_from_m_CS_SG
! ============= 2014/08/14

! == KT_add === 2014/08/26
  use m_Control_Parameters,   only : SpinOrbit_mode, noncol, iprisym
  use m_Const_Parameters,    only : Neglected
  use m_Ionic_System,  only :   mag_moment0_atoms, ionic_charge_atoms, &
       &                        mag_moment0_atoms_is_defined
! ============= 2014/08/26

! === KT_add === 13.2S
  use m_Crystal_Structure,   only : sw_use_magnetic_symmetry
! ============== 13.2S
  use mpi

  implicit none

!  include 'mpif.h'
  integer istatus(mpi_status_size)

! ---------------------------
  integer :: nopr_before_sym_reduction
  real(kind=DP), allocatable :: op_before_sym_reduction(:,:,:)
  integer :: invop(48)
  integer :: determinant_op(48)
! -------------------------

! == KT_add === 2014/08/14
  integer, allocatable :: magmom_dir_inversion_opr_flag(:)
! ============= 2014/08/14

  complex(kind=CMPLDP), allocatable :: op_spinor(:,:,:)

contains

! ------------------------------------------------------------------------
!!
!!!   Removing sym. operations by considering magnetic moment
!!
! ------------------------------------------------------------------------
  subroutine m_CS_set_Magnetic_Sym
    integer :: nsym_with_magmom
    real(kind=DP), allocatable :: mag_loc(:,:), pos2(:,:), ion_chg(:)

    if (  ( .not. noncol ) .and. ( sw_use_magnetic_symmetry == OFF ) ) then
       if ( .not. mag_moment0_atoms_is_defined ) return
    endif

    nopr_before_sym_reduction = nopr

    allocate( pos2(natm2,3) );  pos2 = 0.0d0
    call set_val_pos2( pos2 )

    if ( mag_moment0_atoms_is_defined ) then
       allocate( ion_chg(natm2) ); ion_chg = 0.0d0
       call set_val_ion_charge( ion_chg )
    endif

    if ( noncol .or. sw_use_magnetic_symmetry == ON ) then
       allocate( mag_loc(natm2,3) ); mag_loc = 0.0d0
       call set_val_mag_loc_noncol( mag_loc )

       if ( noncol ) then
          if ( .not. allocated( op_before_sym_reduction ) ) then
             allocate( op_before_sym_reduction(3,3,nopr) )
          endif
          op_before_sym_reduction = op
       endif

! == KT_add === 2014/08/14
       if ( allocated( magmom_dir_inversion_opr_flag ) ) then
          deallocate( magmom_dir_inversion_opr_flag )
       endif
       allocate( magmom_dir_inversion_opr_flag(nopr) )
! ============= 2014/08/14

! ==== KT_add === 2014/09/26
       if ( sw_neglect_magmom == ON ) then
          magmom_dir_inversion_opr_flag = 1;  return
       endif
! =============== 2014/09/26

       call set_magnetic_symm_noncl( natm2, pos2, nopr, ig01, op, tau, mag_loc, &
            &                  nsym_with_magmom, magmom_dir_inversion_opr_flag, &
            &                  ion_chg )

    else
       allocate( mag_loc(natm2,1) ); mag_loc = 0.0d0
       call set_val_mag_loc_col( mag_loc )
       call set_magnetic_symm_col( natm2, pos2, nopr, ig01, op, tau, mag_loc, &
            &                      nsym_with_magmom, ion_chg )
    endif

    call resize_matrix_size( nopr, nsym_with_magmom )
#if 0
    call set_magnetic_tspace_gnerators( nopr, ig01, op, tau )
#endif

    deallocate( pos2 ); deallocate( mag_loc )
    if ( allocated( ion_chg ) ) deallocate( ion_chg )

  contains

    subroutine set_val_pos2( pos2 )
      real(kind=DP), intent(out) :: pos2(natm2,3)
      integer :: i, ia

      pos2(1:natm,1:3) = pos(1:natm,1:3)

      i = natm
      do ia=1,natm
         if (iwei(ia)/=1) then
            i = i+1
            pos2(i,1:3) = -pos(ia,1:3)
         end if
      end do
    end subroutine set_val_pos2

    subroutine set_val_ion_charge( ion_chg )
      real(kind=DP), optional, intent(out) :: ion_chg(natm2)

      integer :: i, ia

      i = natm
      ion_chg(1:natm) = ionic_charge_atoms(1:natm)
      do ia=1,natm
         if (iwei(ia)/=1) then
            i = i+1
            ion_chg(i) = ionic_charge_atoms(ia)
         endif
      end do
    end subroutine set_val_ion_charge

    subroutine set_val_mag_loc_noncol( mag_loc )
      real(kind=DP), intent(out) :: mag_loc(natm2,3)

      integer :: i, ia

      i = natm
      mag_loc(1:natm,:) = magmom_local_now(1:natm,:)
      do ia=1,natm
         if (iwei(ia)/=1) then
            i = i+1
            mag_loc(i,1:3) = magmom_local_now(ia,1:3)
         end if
      end do
!
!      write(nfout,*) 'mag_loc now '
!      Do i=1, natm
!         write(nfout,'(I4,3F20.15)') i, mag_loc(i,1:3)
!      End Do

    end subroutine set_val_mag_loc_noncol

    subroutine set_val_mag_loc_col( mag_loc )
      real(kind=DP), intent(out) :: mag_loc(natm2,1)

      integer :: i, ia

      i = natm
      Do ia=1, natm
         mag_loc(ia,1) = mag_moment0_atoms(ia,1)
      End Do
      do ia=1,natm
         if (iwei(ia)/=1) then
            i = i+1
            mag_loc(i,1) = mag_moment0_atoms(ia,1)
         end if
      end do
    end subroutine set_val_mag_loc_col

    subroutine set_magnetic_symm_noncl( natom, apos, nsym, nopr, rot, tau, mag_loc, &
         &                              nsym_with_mag, magmom_dir_inversion_opr_flag, &
         &                              ion_chg )
      integer, intent(in) :: natom
      integer, intent(inout) :: nsym
      integer, intent(out) :: nopr(48), nsym_with_mag
      integer, intent(inout) :: magmom_dir_inversion_opr_flag(nsym)

! --------------------- kDEBUG --------------------------------- 20121020
!      real(kind=DP), intent(inout) :: rot(3,3,48), tau(3,48,CRDTYP)
      real(kind=DP), intent(inout) :: rot(3,3,nsym), tau(3,nsym,CRDTYP)
! --------------------- kDEBUG --------------------------------- 20121020

      real(kind=DP), intent(in) :: apos(natom,3)
      real(kind=DP), intent(in) :: mag_loc(natom,3)
      real(kind=DP), optional, intent(in) :: ion_chg(natom)

      logical :: mag_sym_flag(nsym)
      integer :: isym, count
      real(kind=DP) :: determinant
! ----
!      write(*,*) 'sw_reduce_ = ',  sw_reduce_sym_by_magmom
      if ( sw_reduce_sym_by_magmom==ON ) then
#if 0
         call check_if_op_has_magnetic_sym1( natom, apos, nsym, rot, tau, &
              &                              mag_loc, mag_sym_flag, &
              &                              magmom_dir_inversion_opr_flag, ion_chg )
#else
         call check_if_op_has_magnetic_sym2( natom, apos, nsym, rot, tau, &
              &                              mag_loc, mag_sym_flag, &
              &                              magmom_dir_inversion_opr_flag, ion_chg )
#endif
      else
         mag_sym_flag = .true.
      endif

      if ( sw_reduce_sym_by_orbital==ON ) then
         call check_if_op_has_sym_for_orb( natom, nsym, rot, mag_loc, mag_sym_flag )
      endif
!
      if ( sw_allow_improper_rotation == OFF ) then
         Do isym=1, nsym
            if ( mag_sym_flag(isym) ) then
               call calc_determinant( rot(:,:,isym), determinant )
               if ( determinant < 0 ) mag_sym_flag(isym) = .false.
            endif
         end Do
      endif
      if ( sw_allow_mag_sym_inversion == OFF ) then
         Do isym=1, nsym
            if ( mag_sym_flag(isym) ) then
               if ( magmom_dir_inversion_opr_flag(isym) == -1 ) then
                  mag_sym_flag(isym) = .false.
               endif
            endif
         End Do
      endif
! -----
      write(nfout,*) '!*** magnetic symmetry info'
      Do isym=1, nsym
         write(nfout,*) isym, mag_sym_flag(isym), magmom_dir_inversion_opr_flag(isym)
      End do

      count = 0
      Do isym=1, nsym
         if ( mag_sym_flag(isym) ) then
            count = count + 1
            nopr(count) = nopr(isym)
            rot(:,:,count) = rot(:,:,isym)
            tau(:,count,:) = tau(:,isym,:)
            magmom_dir_inversion_opr_flag(count) = magmom_dir_inversion_opr_flag(isym)
         endif
      End do

      nsym_with_mag = count

    end subroutine set_magnetic_symm_noncl

    subroutine set_magnetic_symm_col( natom, apos, nsym, nopr, rot, tau, mag_loc, &
         &                            nsym_with_mag, ion_chg )
      integer, intent(in) :: natom
      integer, intent(inout) :: nsym
      integer, intent(out) :: nopr(48), nsym_with_mag

! --------------------- kDEBUG --------------------------------- 20121020
!      real(kind=DP), intent(inout) :: rot(3,3,48), tau(3,48,CRDTYP)
      real(kind=DP), intent(inout) :: rot(3,3,nsym), tau(3,nsym,CRDTYP)
! --------------------- kDEBUG --------------------------------- 20121020

      real(kind=DP), intent(in) :: apos(natom,3)
      real(kind=DP), intent(in) :: mag_loc(natom,3)
      real(kind=DP), optional, intent(in) :: ion_chg(natm)

      logical :: mag_sym_flag(nsym)
      integer :: isym, count
      
      call check_if_op_has_magnetic_sym0( natom, apos, nsym, rot, tau, &
              &                           mag_loc, mag_sym_flag, ion_chg )

      write(nfout,*) '!*** magnetic symmetry info'
      Do isym=1, nsym
         write(nfout,*) isym, mag_sym_flag(isym)
      End do

      count = 0
      Do isym=1, nsym
         if ( mag_sym_flag(isym) ) then
            count = count + 1
            nopr(count) = nopr(isym)
            rot(:,:,count) = rot(:,:,isym)
            tau(:,count,:) = tau(:,isym,:)
         endif
      End do

      nsym_with_mag = count

    end subroutine set_magnetic_symm_col

    subroutine check_if_op_has_sym_for_orb( natom, nsym, rot, &
         &                                  mag_loc, mag_sym_flag )
      integer, intent(in) :: natom, nsym
      logical, intent(inout) :: mag_sym_flag(48)
      real(kind=DP), intent(inout) :: rot(3,3,nsym)
      real(kind=DP), intent(in) :: mag_loc(natom,3)

      integer :: isym, it, ia
      real(kind=DP) :: c1, c2
      real(kind=DP) :: vec1(3), vec_tmp(3), dvec_m(3), dvec_p(3)
      real(kind=DP), parameter :: criterion = 1.0D-5

      Do isym=1, nsym
         if ( .not. mag_sym_flag(isym) ) cycle

         Do it=1, ntyp
            if ( has_partially_filled_lcore(it) == 0 ) cycle

            Do ia=1, natom
               if ( ityp(ia) /=it ) cycle
               
               call calc_normal_vector( mag_loc(ia,1:3), vec1(1:3) )
               vec_tmp(1:3) = matmul( rot(:,:,isym), vec1(1:3) )
!
               dvec_m(1:3) = vec_tmp(1:3) - vec1(1:3)
               dvec_p(1:3) = vec_tmp(1:3) + vec1(1:3)
!
               c1 = abs( dvec_m(1) ) +abs( dvec_m(2) ) +abs( dvec_m(3) )
               c2 = abs( dvec_p(1) ) +abs( dvec_p(2) ) +abs( dvec_p(3) )
!
               if ( c1 < criterion .or. c2 < criterion ) then
                  mag_sym_flag(isym) = .true.
               else
                  mag_sym_flag(isym) = .false.
               endif

               exit
            End do
         End Do

      End Do
!
    end subroutine check_if_op_has_sym_for_orb

    subroutine calc_normal_vector( vec_in, vec_out )
      real(kind=DP), intent(in)  :: vec_in(3)
      real(kind=DP), intent(out) :: vec_out(3)
!
      integer :: k
      real(kind=DP) :: c1, cnorm, vec_tmp(3)
      real(kind=DP), parameter :: criterion = 1.0D-5

! --- trial 1 --
      vec_tmp(1) = 1.0;  vec_tmp(2) = 2.0;    vec_tmp(3) = 3.0
!
      c1 = 0.0d0; cnorm = 0.0d0
      Do k=1, 3
         c1 = c1 + vec_tmp(k)*vec_in(k)
      End do
      vec_tmp(:) = vec_tmp(:) - c1*vec_in(:)
      Do k=1, 3
         cnorm = cnorm + vec_tmp(k)**2
      End do
      cnorm = sqrt(cnorm)
!
      if ( cnorm > criterion ) then
         vec_out(:) = vec_tmp(:) / cnorm
!         write(*,*) 'vec_out = ', vec_out
         return
      endif

! --- trial 2 --
      vec_tmp(1) = 3.0;  vec_tmp(2) = 2.0;    vec_tmp(3) = 1.0
!
      c1 = 0.0d0; cnorm = 0.0d0
      Do k=1, 3
         c1 = c1 + vec_tmp(k)*vec_in(k)
      End do
      vec_tmp(:) = vec_tmp(:) - c1*vec_in(:)
      Do k=1, 3
         cnorm = cnorm + vec_tmp(k)**2
      End do
      cnorm = sqrt(cnorm)
!
      if ( cnorm > criterion ) then
         vec_out(:) = vec_tmp(:) / cnorm
      endif

    end subroutine calc_normal_vector

    subroutine check_if_op_has_magnetic_sym1( natom, apos, nsym, rot, &
         &                                    tau, mag_loc, mag_sym_flag, &
         &                                    magmom_dir_inversion_opr_flag, &
         &                                    ion_chg )
      integer, intent(in) :: natom, nsym
      logical, intent(inout) :: mag_sym_flag(nsym)
      integer, intent(inout) :: magmom_dir_inversion_opr_flag(nsym)
      real(kind=DP), intent(in) :: apos(natom,3), mag_loc(natom,3)
      real(kind=DP), optional, intent(in) :: ion_chg(natm)

! --------------------- kDEBUG --------------------------------- 20121020
!      real(kind=DP), intent(inout) :: rot(3,3,48), tau(3,48,CRDTYP)
      real(kind=DP), intent(inout) :: rot(3,3,nsym), tau(3,nsym,CRDTYP)
! --------------------- kDEBUG --------------------------------- 20121020

      real(kind=DP), parameter :: criterion = 1.0D-5
      real(kind=DP), allocatable :: apos_t(:,:)

      integer :: i, ia, ja, ja_found, it, isym, count
      real(kind=DP) :: coord_tmp(3), dcoord(3), dist
      real(kind=DP) :: mag_tmp(3), dmag_p(3), dmag_m(3)
      real(kind=DP) :: c1, c2, determinant
!
      real(kind=DP), allocatable :: rot_pr(:,:,:)
!
! -- init --
      allocate( rot_pr(3,3,nsym) ); rot_pr = 0.0d0
      call m_CS_op_in_PUCD( nfout, rot_pr, nsym, .true., &
           &                '***** Before considering magnetic symmetry ****** ')

      allocate(apos_t(natom,3))
      apos_t(1:natm,1:3) = apos(1:natom,1:3)
      do ia = 1, natom
         do i = 1, 3
            apos_t(ia,i) = apos_t(ia,i) - floor(apos_t(ia,i))
         end do
      end do

      mag_sym_flag = .false.

      magmom_dir_inversion_opr_flag = 1

! -- begin --
      Do isym=1, nsym

! ================== KT_add ======= 2014/08/14
         if ( lattice_system_from_m_CS_SG == "cubic" ) then
            if ( ig01(isym) > 24 ) cycle
         else if ( lattice_system_from_m_CS_SG == "hexagonal" ) then
            if ( ig01(isym) > 12 ) cycle
         endif
! ================================= 2014/08/14
         call calc_determinant( rot_pr(:,:,isym), determinant )

         Loop_it : Do it=1, ntyp
            Do ia=1, natom
               if ( ityp(ia) /=it ) cycle

               coord_tmp(1:3) = matmul( rot_pr(:,:,isym), apos_t(ia,1:3) ) &
                    &         + tau(1:3,isym,BUCS)
               coord_tmp(1:3) = coord_tmp(1:3) - floor(coord_tmp(1:3))

               ja_found = 0
               Do ja=1, natom
                  dcoord(:) = apos_t(ja,1:3) - coord_tmp(:)

! === KT_mod === 2015/01/05
!                  dist = dcoord(1)**2 + dcoord(2)**2 + dcoord(3)**2
!                  dist = sqrt(dist)
!                  dist = dist - floor( dist +DELTA10 )

                  dcoord = abs( dcoord )
                  dcoord = dcoord - floor(dcoord + DELTA10)
                  dist = sum(dcoord)
! ============== 2015/01/05

                  if ( dist < criterion ) then
                     ja_found = ja
                     !  goto 100
                     exit
                  endif
               End do

100            continue
               !
               if ( ja_found ==0 ) then
                  write(*,*) 'kt: Not found symmetry', isym, ia
                  !stop
                  call phase_error_with_msg(nfout,'kt: Not found symmetry',__LINE__,__FILE__)
               endif
! ----
               if ( present(ion_chg) ) then
                  if ( abs( ion_chg(ia)-ion_chg(ja_found) ) > 1.0d-7 ) then
                     mag_sym_flag(isym) = .false.;  exit Loop_it
                  endif
               endif
! ----
               mag_tmp(1:3) = matmul( rot(:,:,isym), mag_loc(ia,1:3) )
#ifdef MOMENT_AS_PSEUDO_VECTOR
               if ( determinant < 0 ) mag_tmp = -mag_tmp
#endif
               dmag_m(1:3) = mag_tmp(1:3) - mag_loc(ja_found,1:3)
               dmag_p(1:3) = mag_tmp(1:3) + mag_loc(ja_found,1:3)
!
               c1 = abs( dmag_m(1) ) +abs( dmag_m(2) ) +abs( dmag_m(3) )
               c2 = abs( dmag_p(1) ) +abs( dmag_p(2) ) +abs( dmag_p(3) )
!
! ====== KT_mod === 2014/08/14
!               if ( c1 < criterion ) then
!                  mag_sym_flag(isym) = .true.
!                  exit
!               endif
!               if ( mag_sym_flag(isym) ) then
!                  if ( c1 > criterion ) then
!                     mag_sym_flag(isym) = .false.
!                  endif
!               endif
! --
               if ( c1 < criterion ) then
                  mag_sym_flag(isym) = .true.
                  magmom_dir_inversion_opr_flag(isym) = 1
                  exit
               endif
               if ( c2 < criterion ) then
                  mag_sym_flag(isym) = .true.
                  magmom_dir_inversion_opr_flag(isym) = -1
                  exit
               endif
               if ( mag_sym_flag(isym) ) then
                  if ( c1 > criterion .and. c2 > criterion ) then
                     mag_sym_flag(isym) = .false.
                  endif
               endif
! =============== 2014/08/14
            End do
         End do Loop_it
      End do
      deallocate( rot_pr )
    end subroutine check_if_op_has_magnetic_sym1

    subroutine check_if_op_has_magnetic_sym2( natom, apos, nsym, rot, &
         &                                    tau, mag_loc, mag_sym_flag, &
         &                                    magmom_dir_inversion_opr_flag, &
         &                                    ion_chg )
      integer, intent(in) :: natom, nsym
      logical, intent(inout) :: mag_sym_flag(nsym)
      integer, intent(inout) :: magmom_dir_inversion_opr_flag(nsym)
      real(kind=DP), intent(in) :: apos(natom,3), mag_loc(natom,3)
      real(kind=DP), optional, intent(in) :: ion_chg(natm)

! --------------------- kDEBUG --------------------------------- 20121020
!      real(kind=DP), intent(inout) :: rot(3,3,48), tau(3,48,CRDTYP)
      real(kind=DP), intent(inout) :: rot(3,3,nsym), tau(3,nsym,CRDTYP)
! --------------------- kDEBUG --------------------------------- 20121020

      real(kind=DP), parameter :: criterion = 1.0D-5
      real(kind=DP), allocatable :: apos_t(:,:)

      integer :: i, ia, ja, ja_found, it, isym, count
      real(kind=DP) :: coord_tmp(3), dcoord(3), dist
      real(kind=DP) :: mag_tmp(3), dmag_p(3), dmag_m(3)
      real(kind=DP) :: c1, c2, m_norm, determinant
      real(kind=DP), allocatable :: rot_pr(:,:,:)
!
      logical :: first
!
! -- init --
      allocate( rot_pr(3,3,nsym) ); rot_pr = 0.0d0
      call m_CS_op_in_PUCD( nfout, rot_pr, nsym, .true., &
           &                '***** Before considering magnetic symmetry (2) ****** ')

      allocate(apos_t(natom,3))
      apos_t(1:natm,1:3) = apos(1:natom,1:3)
      do ia = 1, natom
         do i = 1, 3
            apos_t(ia,i) = apos_t(ia,i) - floor(apos_t(ia,i))
         end do
      end do

      mag_sym_flag = .true.
      magmom_dir_inversion_opr_flag = 1

! -- begin --
      Do isym=1, nsym
#if 0
         if ( SpinOrbit_mode /= Neglected ) then
            if ( lattice_system_from_m_CS_SG == "cubic" ) then
               if ( ig01(isym) > 24 ) then
                  mag_sym_flag(isym) = .false. ;  cycle
               endif
            else if ( lattice_system_from_m_CS_SG == "hexagonal" ) then
               if ( ig01(isym) > 12 ) then
                  mag_sym_flag(isym) = .false. ;  cycle
               endif
            endif
         endif
#endif
         call calc_determinant( rot_pr(:,:,isym), determinant )

         first = .true.

         Loop_it : Do it=1, ntyp
            Do ia=1, natom
               if ( ityp(ia) /=it ) cycle

               m_norm = sqrt( mag_loc(ia,1)**2 +mag_loc(ia,2)**2 +mag_loc(ia,3)**2 )
               if ( m_norm < 1.0D-8 ) cycle

               coord_tmp(1:3) = matmul( rot_pr(:,:,isym), apos_t(ia,1:3) ) &
                    &         + tau(1:3,isym,BUCS)
               coord_tmp(1:3) = coord_tmp(1:3) - floor(coord_tmp(1:3))

               ja_found = 0
               Loop_ja : Do ja=1, natom
                  dcoord(:) = apos_t(ja,1:3) - coord_tmp(:)

! === KT_mod === 2015/01/05
!                  dist = dcoord(1)**2 + dcoord(2)**2 + dcoord(3)**2
!                  dist = sqrt(dist)
!                  dist = dist - floor( dist +DELTA10 )

                  dcoord = abs( dcoord )
                  dcoord = dcoord - floor(dcoord + DELTA10)
                  dist = sum(dcoord)
! ============== 2015/01/05

                  if ( dist < criterion ) then
                     ja_found = ja;  exit Loop_ja
                  endif
               End do Loop_ja

100            continue
               !
               if ( ja_found ==0 ) then
                  write(*,*) 'kt: Not found symmetry', isym, ia
                  !stop
                  call phase_error_with_msg(nfout,'kt: Not found symmetry',__LINE__,__FILE__)
               endif
! ----
               if ( present(ion_chg) ) then
                  if ( abs( ion_chg(ia)-ion_chg(ja_found) ) > 1.0d-5 ) then
                     mag_sym_flag(isym) = .false.;  exit Loop_it
                  endif
               endif
#if 0
               if ( abs(mag_loc(ia,1)-mag_loc(ja_found,1) )  > 1.0d-5 ) then
                  mag_sym_flag(isym) = .false.;  exit Loop_it
               endif
               if ( abs(mag_loc(ia,2)-mag_loc(ja_found,2) )  > 1.0d-5 ) then
                  mag_sym_flag(isym) = .false.;  exit Loop_it
               endif
               if ( abs(mag_loc(ia,3)-mag_loc(ja_found,3) )  > 1.0d-5 ) then
                  mag_sym_flag(isym) = .false.;  exit Loop_it
               endif
#endif

! --
               mag_tmp(1:3) = matmul( rot(:,:,isym), mag_loc(ia,1:3) )
#ifdef MOMENT_AS_PSEUDO_VECTOR
               if ( determinant < 0 ) mag_tmp = -mag_tmp
#endif
               dmag_m(1:3) = mag_tmp(1:3) - mag_loc(ja_found,1:3)
               dmag_p(1:3) = mag_tmp(1:3) + mag_loc(ja_found,1:3)
!
               c1 = abs( dmag_m(1) ) +abs( dmag_m(2) ) +abs( dmag_m(3) )
               c2 = abs( dmag_p(1) ) +abs( dmag_p(2) ) +abs( dmag_p(3) )

               if ( c1 > criterion .and. c2 > criterion ) then
                  mag_sym_flag(isym) = .false.;  exit Loop_it
               endif
               if ( first ) then
                  if ( c1 < criterion ) then
                     magmom_dir_inversion_opr_flag(isym) = 1
                  else if ( c2 < criterion ) then
                     magmom_dir_inversion_opr_flag(isym) = -1
                  endif
                  first = .false.
               endif

               if ( magmom_dir_inversion_opr_flag(isym) == 1 &
                    &      .and. c2 < criterion ) then
                  mag_sym_flag(isym) = .false.;  exit Loop_it
               endif
               if ( magmom_dir_inversion_opr_flag(isym) == -1 &
                    &      .and. c1 < criterion ) then
                  mag_sym_flag(isym) = .false.;   exit Loop_it
               endif

            End do
         End do Loop_it
      End do
!
      deallocate( rot_pr )

!      write(nfout,*) '!*** magnetic symmetry info'
!      Do isym=1, nsym
!         write(nfout,*) isym, mag_sym_flag(isym), magmom_dir_inversion_opr_flag(isym)
!      End do

    end subroutine check_if_op_has_magnetic_sym2

    subroutine check_if_op_has_magnetic_sym0( natom, apos, nsym, rot, &
         &                                    tau, mag_loc, mag_sym_flag, &
         &                                    ion_chg )
      integer, intent(in) :: natom, nsym
      logical, intent(inout) :: mag_sym_flag(nsym)
      real(kind=DP), intent(in) :: apos(natom,3), mag_loc(natom,1)
      real(kind=DP), optional, intent(in) :: ion_chg(natm)

! --------------------- kDEBUG --------------------------------- 20121020
!      real(kind=DP), intent(inout) :: rot(3,3,48), tau(3,48,CRDTYP)
      real(kind=DP), intent(inout) :: rot(3,3,nsym), tau(3,nsym,CRDTYP)
! --------------------- kDEBUG --------------------------------- 20121020

      real(kind=DP), parameter :: criterion = 1.0D-5
      real(kind=DP), allocatable :: apos_t(:,:)

      integer :: i, ia, ja, ja_found, it, isym, count
      real(kind=DP) :: coord_tmp(3), dcoord(3), dist
      real(kind=DP) :: c1, c2, m_norm
      real(kind=DP), allocatable :: rot_pr(:,:,:)
!
      logical :: first
!
! -- init --
      allocate( rot_pr(3,3,nsym) ); rot_pr = 0.0d0
      call m_CS_op_in_PUCD( nfout, rot_pr, nsym, .true., &
           &                '***** Before considering magnetic symmetry ****** ')

      allocate(apos_t(natom,3))
      apos_t(1:natm,1:3) = apos(1:natom,1:3)
      do ia = 1, natom
         do i = 1, 3
            apos_t(ia,i) = apos_t(ia,i) - floor(apos_t(ia,i))
         end do
      end do

      mag_sym_flag = .true.

! -- begin --
      Do isym=1, nsym
         first = .true.

         Loop_it : Do it=1, ntyp
            Do ia=1, natom
               if ( ityp(ia) /=it ) cycle

!               m_norm = sqrt( mag_loc(ia,1)**2 +mag_loc(ia,2)**2 +mag_loc(ia,3)**2 )
               m_norm = abs( mag_loc(ia,1) )
               if ( m_norm < 1.0D-8 ) cycle

               coord_tmp(1:3) = matmul( rot_pr(:,:,isym), apos_t(ia,1:3) ) &
                    &         + tau(1:3,isym,BUCS)
               coord_tmp(1:3) = coord_tmp(1:3) - floor(coord_tmp(1:3))

               ja_found = 0
               Loop_ja : Do ja=1, natom
                  dcoord(:) = apos_t(ja,1:3) - coord_tmp(:)

! === KT_mod === 2015/01/05
!                  dist = dcoord(1)**2 + dcoord(2)**2 + dcoord(3)**2
!                  dist = sqrt(dist)
!                  dist = dist - floor( dist +DELTA10 )

                  dcoord = abs( dcoord )
                  dcoord = dcoord - floor(dcoord + DELTA10)
                  dist = sum(dcoord)
! ============== 2015/01/05

                  if ( dist < criterion ) then
                     ja_found = ja;  exit Loop_ja
                  endif
               End do Loop_ja

100            continue
               !
               if ( ja_found ==0 ) then
                  write(*,*) 'kt: Not found symmetry', isym, ia
                  !stop
                  call phase_error_with_msg(nfout,'kt: Not found symmetry',__LINE__,__FILE__)
               endif
! ----
               if ( present(ion_chg) ) then
                  if ( abs( ion_chg(ia)-ion_chg(ja_found) ) > 1.0d-7 ) then
                     mag_sym_flag(isym) = .false.;  exit Loop_it
                  endif
               endif
! --
               if ( abs( mag_loc(ia,1)-mag_loc(ja_found,1) ) > 1.0d-7 ) then
                  mag_sym_flag(isym) = .false.;  exit Loop_it
               endif

            End do
         End do Loop_it
      End do
!
      deallocate( rot_pr )

!      write(nfout,*) '!*** magnetic symmetry info'
!      Do isym=1, nsym
!         write(nfout,*) isym, mag_sym_flag(isym)
!      End do

    end subroutine check_if_op_has_magnetic_sym0

    subroutine resize_matrix_size( nsym, nsym_with_mag )
      integer, intent(inout) :: nsym, nsym_with_mag
      real(kind=DP), allocatable :: rot_tmp(:,:,:), tau_tmp(:,:,:)
      integer, allocatable :: magmom_dir_inversion_opr_flag_tmp(:)
!
      allocate( rot_tmp(3,3,nsym) );      rot_tmp = 0.0d0
      allocate( tau_tmp(3,nsym,CRDTYP) ); tau_tmp = 0.0d0
!
      rot_tmp = op;  tau_tmp = tau
      deallocate( op );  deallocate( tau )
!
      allocate( op(3,3,nsym_with_mag ) );       op = 0.0d0
      allocate( tau(3,nsym_with_mag,CRDTYP ) ); tau = 0.0d0
!
      op(:,:,1:nsym_with_mag) = rot_tmp(:,:,1:nsym_with_mag)
      tau(:,1:nsym_with_mag,:) = tau_tmp(:,1:nsym_with_mag,:)
!
      deallocate( rot_tmp ); deallocate( tau_tmp )
!
      if ( noncol ) then
! === KT_add === 2014/08/14
         allocate( magmom_dir_inversion_opr_flag_tmp( nsym ) )
         magmom_dir_inversion_opr_flag_tmp = 0
         
         magmom_dir_inversion_opr_flag_tmp = magmom_dir_inversion_opr_flag
         deallocate( magmom_dir_inversion_opr_flag )

         allocate( magmom_dir_inversion_opr_flag(nsym_with_mag) )
         magmom_dir_inversion_opr_flag = 0
         
         magmom_dir_inversion_opr_flag(1:nsym_with_mag) &
              &         = magmom_dir_inversion_opr_flag_tmp(1:nsym_with_mag)
         deallocate( magmom_dir_inversion_opr_flag_tmp )
      endif
! ============== 2014/08/14

      nsym = nsym_with_mag

    end subroutine resize_matrix_size

    subroutine set_magnetic_tspace_gnerators( nsym, nopr, rot, tau )
      integer, intent(in) :: nsym
      integer, intent(in) :: nopr(48)
      real(kind=DP), intent(in) :: rot(3,3,nsym), tau(3,nsym,CRDTYP)
      
      integer :: iptab(nsym,nsym)
      integer :: ng, ig(3)
      character(len=9) :: system

      system = pg_symbol_system

      call set_mini_product_table( nsym, rot, tau, iptab )
      call search_generators( nsym, iptab, ng, ig )
      call set_generators( ng, ig, nopr )
      call write_tspace_generators(system)

    end subroutine set_magnetic_tspace_gnerators

! --
    subroutine set_mini_product_table( nsym, rot, tau, iptab )
      integer, intent(in) :: nsym
      integer, intent(out) :: iptab(nsym,nsym)
      real(kind=DP), intent(in) :: rot(3,3,nsym), tau(3,nsym,CRDTYP)
      
      real(kind=DP) :: rot_pr(3,3,48)
      real(kind=DP) :: op_tmp(3,3), tau_tmp(3)
      real(kind=DP) :: c1, c2, cx, cy, cz
      real(kind=DP), parameter :: eps = 1.d-4

      integer :: isym, jsym, ksym
      integer :: ksym_found
      integer :: m1, m2, m3
      integer :: nx, ny, nz

      iptab = 0

      call m_CS_op_in_PUCD( nfout, rot_pr, nsym, .false. )

      Do jsym=1, nsym
         Do isym=1, nsym
            
            op_tmp = 0.0d0; tau_tmp = 0.0d0
            Do m1=1, 3
               Do m2=1, 3
                  Do m3=1, 3
                     op_tmp(m1,m2) = op_tmp(m1,m2) &
                          & + rot_pr(m1,m3,isym) *rot_pr(m3,m2,jsym) 
                  End do
               End Do
            End do
!            Do m1=1, 3
!               Do m2=1, 3
!                  tau_tmp(m1) = tau_tmp(m1) &
!                       &       + rot_pr(m1,m2,isym)*tau(m2,jsym,BUCS)
!               End do
!               tau_tmp(m1) = tau_tmp(m1) + tau(m1,isym,BUCS)
!            End do
! ---------------------------------
            ksym_found = 0
! ---------------------------------

!            write(*,*) 'isym jsym ', isym, jsym
!            write(*,*) 'op_tmp ', op_tmp(1,1), op_tmp(1,2), op_tmp(1,3)
!            write(*,*) 'op_tmp ', op_tmp(2,1), op_tmp(2,2), op_tmp(2,3)
!            write(*,*) 'op_tmp ', op_tmp(3,1), op_tmp(3,2), op_tmp(3,3)
!            write(*,*) 'tau_tmp ', tau_tmp(1), tau_tmp(2), tau_tmp(3)


            Do ksym=1, nsym
               c1 = 0.0d0
               Do m1=1, 3
                  Do m2=1, 3
                     c1 = c1 + abs( op_tmp(m1,m2) - rot_pr(m1,m2,ksym) )
                  End do
               End do

               if ( c1 < eps ) then
                  Do nx=-1, 1
                     Do ny=-1, 1
                        Do nz=-1, 1
!                           cx = tau_tmp(1) - tau(1,ksym,BUCS) +nx
!                           cy = tau_tmp(2) - tau(2,ksym,BUCS) +ny
!                           cz = tau_tmp(3) - tau(3,ksym,BUCS) +nz

                           cx = tau(1,jsym,BUCS) +nx
                           cy = tau(2,jsym,BUCS) +ny
                           cz = tau(3,jsym,BUCS) +nz
                           Do m1=1, 3
                              tau_tmp(m1) = rot_pr(m1,1,isym)*cx &
                                   &       +rot_pr(m1,2,isym)*cy &
                                   &       +rot_pr(m1,3,isym)*cz &
                                   &       +tau(m1,isym,BUCS)
                           End do
                           cx = tau_tmp(1) - tau(1,ksym,BUCS)
                           cy = tau_tmp(2) - tau(2,ksym,BUCS)
                           cz = tau_tmp(3) - tau(3,ksym,BUCS)
                           
                           c2 = abs(cx) + abs(cy) + abs(cz)

                           if ( c2 < eps ) then
                              ksym_found = ksym
                              goto 1000
                           endif
                        End Do
                     End Do
                  End Do
               endif
            End Do

1000        continue
            iptab(isym,jsym) = ksym_found
!            write(*,*) 'isym jsym ', isym, jsym, ksym_found
         End do
      End do

    end subroutine set_mini_product_table

    subroutine search_generators( nsym, iptab, ng, ig )
      integer, intent(in) :: nsym, iptab(nsym,nsym)
      integer, intent(out) :: ng, ig(3)
      
      integer :: i,j,k
      logical :: fexist(1:nsym)
      integer :: ii

! --------------
      ng = 1
      do i=1,nsym
         fexist(1:nsym) = .false.
         fexist(i) = .true.
         if (group_is_the_same(nsym,fexist,iptab)) then
            ig(ng) = i
            return
         end if
      end do
! --------------
      ng = 2
      do i=2,nsym
         do j=i+1,nsym
            fexist(1:nsym) = .false.
            fexist(i) = .true.
            fexist(j) = .true.
            if (group_is_the_same(nsym,fexist,iptab)) then
               ig(1)  = i
               ig(ng) = j
               return
            end if
         end do
      end do
! ------------
      ng = 3
      do i=2,nsym
         do j=i+1,nsym
            do k=j+1,nsym
               fexist(1:nsym) = .false.
               fexist(i) = .true.
               fexist(j) = .true.
               fexist(k) = .true.
               if (group_is_the_same(nsym,fexist,iptab)) then
                  ig(1)  = i
                  ig(2)  = j
                  ig(ng) = k
                  return
               end if
            end do
         end do
      end do

      !stop 'Mag: Set of generators of the space group was not found.'
      call phase_error_with_msg(nfout,'Mag: Set of generators of the space group was not found.',__LINE__,__FILE__)

    end subroutine search_generators

    logical function group_is_the_same( nsym, fexist, iptab)
      integer, intent(in) :: nsym
      integer, intent(in) :: iptab(nsym,nsym)
      logical, intent(inout) :: fexist(nsym)

      integer :: i,j
      integer :: n,no
      logical :: flag(nsym)

      n=0
      no=-1 
      do while(n > no)
         flag(1:nsym) = .false.
         do j=1,nsym
            do i=1,nsym
               if(fexist(i).and.fexist(j)) flag(iptab(i,j)) = .true.
            end do
         end do
         do i=1,nsym
            if(flag(i)) fexist(i) = .true.
         end do
         no = n
         n = 0
         do i=1,nsym
            if(fexist(i)) n=n+1
         end do
      end do
      if(n==nsym) then
         group_is_the_same = .true.
      else
         group_is_the_same = .false.
      end if

    end function group_is_the_same

    subroutine set_generators( ng, ig, nopr )
      integer, intent(in) :: ng, ig(3)
      integer, intent(in) :: nopr(48)

      integer :: i,j,k,n
      real(kind=DP) :: t,r,tb(3)
      real(kind=DP), parameter :: eps = 1.d-4

      ngen = ng
      call dealloc_igen_jgen
      call alloc_igen_jgen

      do i=1,ngen
         igen(i) = nopr(ig(i))

         tb(:) = tau(:,ig(i),BUCS)

         do j=1,3
            t = tb(j)
            do k=1,20
               r=t*k
               if(abs(nint(r)-r)<eps) then
                  jgen(1,j,i) = nint(r)
                  jgen(2,j,i) = k
                  exit
               end if
            end do
         end do
      end do
      
    end subroutine set_generators

    subroutine set_translations_zero
      integer :: i,j
      do i=1,ngen
         do j=1,3
            jgen(1,j,i) = 0
            jgen(2,j,i) = 1
         end do
      end do
    end subroutine set_translations_zero

    subroutine write_tspace_generators(system)
      character(len=9), intent(in) :: system
      
      integer :: i
      
      write(nfout,*) '----------------------------------- '
      write(nfout,'("TSPACE Generators with Magnetic symmetry :")')
      
      if (system == 'cubic') then
         do i=1,ngen
            write(nfout,'("igen,jgen(2,3)=",i2,"(",a5,")",3(3x,i2,"/",i2))') igen(i),oh(igen(i)),jgen(:,:,i)
         end do
      else
         do i=1,ngen
            write(nfout,'("igen,jgen(2,3)=",i2,"(",a5,")",3(3x,i2,"/",i2))') igen(i),d6h(igen(i)),jgen(:,:,i)
         end do
      end if

    end subroutine write_tspace_generators

  end subroutine m_CS_set_Magnetic_Sym

! ---------------------------------------------------------------------------
  subroutine m_CS_set_inverse_operation
    integer :: iopr1, iopr2
    real(kind=DP) :: ss(3,3)

    do iopr1=1,nopr
       do iopr2=1,nopr
          ss = matmul(op(:,:,iopr1),op(:,:,iopr2))
          if(abs(ss(1,1)-1.d0)<DELTA07 .and. &
          & abs(ss(2,2)-1.d0)<DELTA07 .and. &
          & abs(ss(3,3)-1.d0)<DELTA07 .and. &
          & abs(ss(1,2))<DELTA07 .and. &
          & abs(ss(1,3))<DELTA07 .and. &
          & abs(ss(2,3))<DELTA07) then
             invop(iopr1) = iopr2
             exit
          end if
       end do
    end do

  end subroutine m_CS_set_inverse_operation

  subroutine m_CS_chk_determinant_op
    integer :: iopr
    real(kind=DP) :: determinant

    do iopr=1,nopr
       call calc_determinant( op(:,:,iopr), determinant )
       determinant_op(iopr) = nint(determinant)
    end do
  end subroutine m_CS_chk_determinant_op

#if 0
  subroutine m_CS_set_op_spinor
    integer :: i, j
    real(kind=DP) :: k1, k2, k3, ux, uy, uz, c1, c2, c3
    real(kind=DP) :: sinth, costh, theta, s1, determinant
!
    real(kind=DP), parameter :: delta = 1.0D-4
    complex(kind=CMPLDP), parameter :: zi = ( 0.0d0, 1.0d0 )
!
    real(kind=DP), allocatable :: op_work(:,:,:)
!
    allocate( op_work(3,3,nopr) );     op_work = op
    if ( .not. allocated( op_spinor ) ) allocate( op_spinor(2,2,nopr) )

    Do i=1, nopr
       call calc_determinant( op(:,:,i), determinant )
       if ( determinant < 0 ) then
          op_work(1,1,i) = -op_work(1,1,i)
          op_work(2,2,i) = -op_work(2,2,i)
          op_work(3,3,i) = -op_work(3,3,i)
       endif

       costh = ( op_work(1,1,i) +op_work(2,2,i) +op_work(3,3,i) -1.0d0 ) /2.0d0
       sinth = sqrt( 1.0d0 -costh**2 )

       theta = acos( costh )
!
       if ( sinth > delta ) then
          k1 = ( op_work(3,2,i) -op_work(2,3,i) ) /2.0d0
          k2 = ( op_work(1,3,i) -op_work(3,1,i) ) /2.0d0
          k3 = ( op_work(2,1,i) -op_work(1,2,i) ) /2.0d0

          ux = k1 /sinth;       uy = k2 /sinth;         uz = k3 /sinth
       else
          if ( costh > 1.0 -delta ) then
             ux = 0.0d0;  uy = 0.0d0;   uz = 1.0d0
             theta = 0.0d0
          else if ( costh < -1.0 +delta ) then
             k1 = ( op_work(1,1,i) + 1.d0 )/2.0d0
             k2 = ( op_work(2,2,i) + 1.d0 )/2.0d0
             k3 = ( op_work(3,3,i) + 1.d0 )/2.0d0

             ux = sqrt(k1);   uy = sqrt(k2);     uz = sqrt(k3)

             if ( uz > delta ) then
                if ( op_work(1,3,i) < 0 )  ux = -ux
                if ( op_work(2,3,i) < 0 )  uy = -uy
             else
                if ( op_work(1,2,i) < 0 )  uy = -uy
             endif

          endif
       endif

#if 0
       write(*,*) "i = ", i, "costh = ", costh, "sinth = ", sinth
       write(*,'(4(A,F10.2))') "Axis: ux = ", ux,  ", uy = ", uy,  ", uz = ", uz, &
            &                  ", Angle = ", theta /PAI *180.0d0
       write(*,*)
#endif

       c1 = cos( theta /2.0d0 );   s1 = sin( theta /2.0d0 );

       op_spinor(1,1,i) = c1 -zi *uz *s1
       op_spinor(1,2,i) = ( -zi *ux -uy ) *s1
       op_spinor(2,1,i) = ( -zi *ux +uy ) *s1
       op_spinor(2,2,i) = c1 +zi *uz *s1
!
       if ( determinant < 0 ) op_spinor(:,:,i) = zi *op_spinor(:,:,i)
    End Do
!
    if ( iprisym > 1 ) then
       write(nfout,*) '*** symmetry operation for spinor ***'
       Do i=1, nopr
          write(nfout,*) ' #symmetry op. = ', i
          write(nfout,'(A,F8.4,A,F8.4,2A,F8.4,A,F8.4,A)') &
               &                '( ', real(op_spinor(1,1,i)), ', ', &
               &                      aimag(op_spinor(1,1,i)), ' I ) ',  &
               &                '( ',real(op_spinor(1,2,i)), ', ', &
               &                      aimag(op_spinor(1,2,i)), ' I )'
          write(nfout,'(A,F8.4,A,F8.4,2A,F8.4,A,F8.4,A)') &
               &                '( ', real(op_spinor(2,1,i)), ', ', &
               &                      aimag(op_spinor(2,1,i)), ' I ) ',  &
               &                '( ',real(op_spinor(2,2,i)), ', ', &
               &                      aimag(op_spinor(2,2,i)), ' I )'
       End Do
    endif
    deallocate( op_work )

  contains

    subroutine calc_determinant( a, determinant )
      real(kind=DP), intent(in) :: a(3,3)
      real(kind=DP), intent(out) :: determinant
      determinant = a(1,1)*( a(2,2)*a(3,3) -a(2,3)*a(3,2) ) &
          &       -a(1,2)*( a(2,1)*a(3,3) -a(2,3)*a(3,1) ) &
           &       +a(1,3)*( a(2,1)*a(3,2) -a(2,2)*a(3,1) )
    end subroutine calc_determinant

  end subroutine m_CS_set_op_spinor
#endif
  
  subroutine calc_determinant( a, determinant )
    real(kind=DP), intent(in) :: a(3,3)
    real(kind=DP), intent(out) :: determinant
    determinant = a(1,1)*( a(2,2)*a(3,3) -a(2,3)*a(3,2) ) &
         &       -a(1,2)*( a(2,1)*a(3,3) -a(2,3)*a(3,1) ) &
         &       +a(1,3)*( a(2,1)*a(3,2) -a(2,2)*a(3,1) )
  end subroutine calc_determinant

  subroutine m_CS_gen_opr_spinor_full( system, nsym, op_spinor )
    character(len=9), intent(in) :: system
    integer, intent(in) :: nsym
    complex(kind=CMPLDP), intent(out) :: op_spinor(2,2,nsym)

    if ( system .eq. 'cubic' ) then
       call cubic_point_group
    else if ( system .eq. 'hexagonal' ) then
       call hexagonal_point_group
    end if

  contains    

    subroutine cubic_point_group
      integer :: n
      real(kind=DP) :: direc(3), angle

      op_spinor = 0.0d0
! E
      n = 1;
      direc(1) = 0.0d0;  direc(2) = 0.0d0;   direc(3) = 1.0d0;  angle = 0.0d0
      call set_rot_mat( direc, angle, op_spinor(:,:,n) )
! C2xyz
      Do n=2, 4
         direc = 0.0d0; angle = PAI2 /2.0d0
         if ( n==2 ) direc(1) = 1.0d0
         if ( n==3 ) direc(2) = 1.0d0
         if ( n==4 ) direc(3) = 1.0d0
         call set_rot_mat( direc, angle, op_spinor(:,:,n) )
      End Do
! C3+
      n = 5
      direc = 1.0d0;  angle = PAI2 /3.0d0
      call set_rot_mat( direc, angle, op_spinor(:,:,n) )
      Do n=6, 8
         direc = -1.0d0;  angle = PAI2 /3.0d0
         if ( n==6 ) direc(3) = 1.0d0
         if ( n==7 ) direc(1) = 1.0d0
         if ( n==8 ) direc(2) = 1.0d0
         call set_rot_mat( direc, angle, op_spinor(:,:,n) )
      End Do
! C3-
      n = 9
      direc = 1.0d0;  angle = -PAI2 /3.0d0
      call set_rot_mat( direc, angle, op_spinor(:,:,n) )
      Do n=10, 12
         direc = -1.0d0;  angle = -PAI2 /3.0d0
         if ( n==10 ) direc(3) = 1.0d0
         if ( n==11 ) direc(1) = 1.0d0
         if ( n==12 ) direc(2) = 1.0d0
         call set_rot_mat( direc, angle, op_spinor(:,:,n) )
      End Do
! C2_*
      Do n=13, 14
         direc(1) = 1.0d0;  direc(3) = 0.0d0;  angle = PAI2 /2.0
         if ( n==13 ) direc(2) =  1.0d0
         if ( n==14 ) direc(2) = -1.0d0
         call set_rot_mat( direc, angle, op_spinor(:,:,n) )
      End Do
      Do n=15, 17, 2
         direc(3) = 1.0d0;  direc(2) = 0.0d0;  angle = PAI2 /2.0
         if ( n==15 ) direc(1) =  1.0d0
         if ( n==17 ) direc(1) = -1.0d0
         call set_rot_mat( direc, angle, op_spinor(:,:,n) )
      End Do
      Do n=16, 18, 2
         direc(1) = 0.0d0;  direc(2) = 1.0d0;  angle = PAI2 /2.0
         if ( n==16 ) direc(3) =  1.0d0
         if ( n==18 ) direc(3) = -1.0d0
         call set_rot_mat( direc, angle, op_spinor(:,:,n) )
      End Do
! C4+       
      Do n=19, 21
         direc = 0.0d0;  angle = PAI2 /4.0d0
         if ( n==19 ) direc(1) = 1.0d0
         if ( n==20 ) direc(2) = 1.0d0
         if ( n==21 ) direc(3) = 1.0d0
         call set_rot_mat( direc, angle, op_spinor(:,:,n) )
      End Do
! C4-
      Do n=22, 24
         direc = 0.0d0;  angle = -PAI2 /4.0d0
         if ( n==22 ) direc(1) = 1.0d0
         if ( n==23 ) direc(2) = 1.0d0
         if ( n==24 ) direc(3) = 1.0d0
         call set_rot_mat( direc, angle, op_spinor(:,:,n) )
      End Do
    !
      Do n=1, 24
         op_spinor(:,:,n+24) = op_spinor(:,:,n)
      End Do
    end subroutine cubic_point_group

    subroutine hexagonal_point_group
      integer :: n
      real(kind=DP) :: direc(3), angle

      op_spinor = 0.0d0
! E
      n = 1;
      direc(1) = 0.0d0;  direc(2) = 0.0d0;   direc(3) = 1.0d0;  angle = 0.0d0
      call set_rot_mat( direc, angle, op_spinor(:,:,n) )
! C6+
      n = 2;
      direc(1) = 0.0d0;  direc(2) = 0.0d0;   direc(3) = 1.0d0;  angle = PAI2 /6.0d0
      call set_rot_mat( direc, angle, op_spinor(:,:,n) )
! C3+
      n = 3;
      direc(1) = 0.0d0;  direc(2) = 0.0d0;   direc(3) = 1.0d0;  angle = PAI2 /3.0d0
      call set_rot_mat( direc, angle, op_spinor(:,:,n) )
! C2
      n = 4;
      direc(1) = 0.0d0;  direc(2) = 0.0d0;   direc(3) = 1.0d0;  angle = PAI2 /2.0d0
      call set_rot_mat( direc, angle, op_spinor(:,:,n) )
! C3-
      n = 5;
      direc(1) = 0.0d0;  direc(2) = 0.0d0;   direc(3) = 1.0d0;  angle = -PAI2 /3.0d0
      call set_rot_mat( direc, angle, op_spinor(:,:,n) )
! C6-
      n = 6;
      direc(1) = 0.0d0;  direc(2) = 0.0d0;   direc(3) = 1.0d0;  angle = -PAI2 /6.0d0
      call set_rot_mat( direc, angle, op_spinor(:,:,n) )
! C2 (ab plane, 90deg. direc)
      n = 7;
      direc(1) = 0.0d0;  direc(2) = 1.0d0;   direc(3) = 0.0d0;  angle = PAI2 /2.0d0
      call set_rot_mat( direc, angle, op_spinor(:,:,n) )
! C2 (ab plane, 30deg. direc)
      n = 8;
      direc(1) = sqrt(3.d0);  direc(2) = 1.d0;   direc(3) = 0.0d0;  angle = PAI2 /2.0d0
      call set_rot_mat( direc, angle, op_spinor(:,:,n) )
! C2 (ab plane, 150deg. direc)
      n = 9;
      direc(1) =-sqrt(3.d0);  direc(2) = 1.d0;   direc(3) = 0.0d0;  angle = PAI2 /2.0d0
      call set_rot_mat( direc, angle, op_spinor(:,:,n) )
! C2 (ab plane, 0deg. direc)
      n = 10;
      direc(1) = 1.0d0;  direc(2) = 0.d0;   direc(3) = 0.0d0;  angle = PAI2 /2.0d0
      call set_rot_mat( direc, angle, op_spinor(:,:,n) )
! C2 (ab plane, 120deg. direc)
      n = 11;
      direc(1) = -1.d0;  direc(2) = sqrt(3.d0);  direc(3) = 0.0d0;  angle = PAI2 /2.0d0
      call set_rot_mat( direc, angle, op_spinor(:,:,n) )
! C2 (ab plane, 60deg. direc)
      n = 12;
      direc(1) =  1.d0;  direc(2) = sqrt(3.d0);  direc(3) = 0.0d0;  angle = PAI2 /2.0d0
      call set_rot_mat( direc, angle, op_spinor(:,:,n) )

      Do n=1, 12
         op_spinor(:,:,n+12) = op_spinor(:,:,n)
      End Do
    end subroutine hexagonal_point_group

    subroutine set_rot_mat( direc, angle, mat )
      real(kind=DP), intent(in) :: direc(3), angle
      complex(kind=CMPLDP), intent(out) :: mat(2,2)

      real(kind=DP) :: s1, c1, norm, vec(3)
      
      norm = sqrt( direc(1)**2 +direc(2)**2 +direc(3)**2 )
      vec(:) = direc(:) /norm

      c1 = cos(angle/2.d0);  s1 = sin(angle/2.d0)
      mat(1,1) = c1 +zi *vec(3) *s1
      mat(1,2) = ( zi *vec(1) +vec(2) ) *s1
      mat(2,1) = ( zi *vec(1) -vec(2) ) *s1
      mat(2,2) = c1 -zi *vec(3) *s1
    end subroutine set_rot_mat

  end subroutine m_CS_gen_opr_spinor_full
  
end module m_CS_Magnetic
