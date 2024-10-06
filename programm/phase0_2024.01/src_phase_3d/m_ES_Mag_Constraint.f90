module m_ES_Mag_Constraint
! $Id: m_ES_Mag_Constraint.f90 574 2017-05-31 03:00:48Z jkoga $

  use m_Parallelization,      only : ista_kngp, iend_kngp, ierr, npes, mype, &
       &                             MPI_CommGroup

  use m_PlaneWaveBasisSet,    only : ngabc, kgp, ngabc_kngp_l

  use m_PseudoPotential,    only : ilmt, nlmt, dl2p, ltp, mtp, taup, prodphi

  use m_Control_Parameters,    only : neg, ndim_spinor, ndim_magmom, noncol, kimg, &
       &                              proj_attribute, num_proj_elems, &
       &                              proj_group

  use m_Const_Parameters,  only : DP, CMPLDP, yes, PAI4, &
       &                          MAG_MOMENT_VALS_GLOBAL, MAG_MOMENT_DIREC_GLOBAL, &
       &                          MAG_MOMENT_VALS_LOCAL,  MAG_MOMENT_DIREC_LOCAL, &
       &                          MAG_MOMENT_VALS_OCCMAT, &
       &                          MAG_MOMENT_DIREC_HARDPART

  use m_Charge_Density,        only : chgq_l, chgqo_l, hsr, hsi
  use m_Ionic_System,         only : natm, ityp,  mag_direction0_atomtyp, &
       &                             magmom_local_now, mag_moment0_atomtyp, &
       &                             pos, iproj_group, mag_moment0_atoms, &
       &                             mag_moment0_atoms_is_defined, &
       &                             mag_direction0_atoms

  use m_CD_Mag_Moment,       only : rad_cov, RhoMag_on_atom

  use m_Crystal_Structure,  only : univol, rltv, altv, &
       &                           mag_constraint_type, &
       &                           mag_moment0_global, &
       &                           mag_direc0_global, &
       &                           mag_constraint_lambda

  use m_IterationNumbers,  only : iteration_electronic
  use m_Electronic_Structure,  only : vlhxc_l, vlhxcQ

  use m_Orbital_Population,    only : om, ommix

! === KT_add === 2014/08/09
  use m_Crystal_Structure,   only : Global_Quantz_Axis_Fixed
! ============== 2014/08/09
  use mpi

  implicit none

!  include 'mpif.h'
  integer istatus(mpi_status_size)

  real(kind=DP) :: MagField_constrain_global(3)
  real(kind=DP), allocatable :: MagField_constrain_local(:,:,:)
  real(kind=DP), allocatable :: MagField_constrain_hardpart(:,:,:,:)

contains

  subroutine m_ES_add_MagConstraintPot_chgql
    if ( ndim_magmom == 1 ) return

    select case ( mag_constraint_type )
    case (MAG_MOMENT_VALS_GLOBAL)
       call case_constraint_moment_global
    case (MAG_MOMENT_DIREC_GLOBAL)
       call case_constraint_direc_global
    case (MAG_MOMENT_VALS_LOCAL)
       call case_constraint_moment_local
    case (MAG_MOMENT_DIREC_LOCAL)
       call case_constraint_direc_local
    end select

  contains

    subroutine case_constraint_moment_global
      real(kind=DP) :: MagMom(3), c1
      integer :: is, ia, it, lmt1

      if ( noncol ) then
         if ( mype == 0 ) then
            Do is=2, ndim_magmom
               MagMom(is-1) = chgq_l( 1,1,is ) *univol
            End do

            MagField_constrain_global(:) = mag_constraint_lambda &
                 &                         *( MagMom(:) -mag_moment0_global(:) )

            Do is=2, ndim_magmom
               vlhxc_l(1,1,is) = vlhxc_l(1,1,is) + MagField_constrain_global(is-1)
            End do
         endif

      else
         if ( mype == 0 ) then
            MagMom(1) = ( chgq_l( 1,1,1 ) -chgq_l( 1,1,2 ) )*univol

            MagField_constrain_global(1) = mag_constraint_lambda &
                 &                        *( MagMom(1) -mag_moment0_global(1) )

            vlhxc_l(1,1,1) = vlhxc_l(1,1,1) +MagField_constrain_global(1)
            vlhxc_l(1,1,2) = vlhxc_l(1,1,2) -MagField_constrain_global(1)
         endif

      endif

    end subroutine case_constraint_moment_global

    subroutine case_constraint_direc_global            ! only for noncol
      real(kind=DP) :: MagMom(3), MagDirec(3)
      real(kind=DP) :: c1, c2, cnorm, cnorm2
      integer :: is, ia, it, ixyz, lmt1

      real(kind=DP), parameter :: cnorm_lower_limit = 1.0D-10

      if ( .not. noncol ) return

      if ( mype == 0 ) then
         Do is=2, ndim_magmom
            MagMom(is-1) = chgq_l( 1,1,is ) *univol
         End do

         cnorm2 = 0.0d0
         Do ixyz=1, 3
            cnorm2 = cnorm2 + MagMom(ixyz)**2
         End do
         cnorm = sqrt( cnorm2 )

         if ( cnorm < cnorm_lower_limit ) then
            MagField_constrain_global = 0.0d0
            return
         endif

         MagDirec = MagMom / cnorm

         Do ixyz=1, 3
            c1 = cnorm2 - MagMom(ixyz)**2
            c2 = cnorm**3
            MagField_constrain_global(ixyz) = &
                 &  mag_constraint_lambda *c1 /c2 &
                 & *( MagDirec(ixyz) - mag_direc0_global(ixyz) )
         End do

         Do is=2, ndim_magmom
            vlhxc_l(1,1,is) = vlhxc_l(1,1,is) + MagField_constrain_global(is-1)
         End do
      endif

    end subroutine case_constraint_direc_global

    subroutine case_constraint_moment_local
      real(kind=DP) :: MagMom(3), cfactor(3)
      real(kind=DP) :: rad1, fac1r, fac1i, fac2
      real(kind=DP) :: VecG(3), normG, normG3, gr, d1

      real(kind=DP), allocatable :: zfcos(:), zfsin(:)
      real(kind=DP), allocatable :: RhoMag_on_atom_mpi(:,:)

      integer :: i, j, ia, it, is, ixyz, ixyz_max, ist

      allocate(zfcos(ista_kngp:iend_kngp)); zfcos = 0.d0
      allocate(zfsin(ista_kngp:iend_kngp)); zfsin = 0.d0

      if ( noncol ) then
         ixyz_max = 3
      else
         ixyz_max = 1
      endif

      if ( allocated( MagField_constrain_local ) ) then
         deallocate( MagField_constrain_local )
      endif

      allocate( MagField_constrain_local(ista_kngp:iend_kngp, kimg, ixyz_max ))
      MagField_constrain_local = 0.0d0
!
      Do ia=1, natm
         it = ityp(ia)
         rad1 = rad_cov(ia)

         Do ixyz=1, ixyz_max
            if ( noncol ) then
               MagMom( ixyz ) = RhoMag_on_atom( ia, ixyz +1 )
            else
               MagMom( ixyz ) = RhoMag_on_atom( ia, 1 ) -RhoMag_on_atom( ia, 2 )
            endif

            if ( mag_moment0_atoms_is_defined ) then
               cfactor( ixyz ) = mag_constraint_lambda /univol &
                    &               *( MagMom(ixyz) -mag_moment0_atoms(ia,ixyz) )
            else
               cfactor( ixyz ) = mag_constraint_lambda /univol &
                    &               *( MagMom(ixyz) -mag_moment0_atomtyp(it,ixyz) )
            endif

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
               MagField_constrain_local( i, 1, ixyz ) &
                    & = MagField_constrain_local( i, 1, ixyz ) &
                    &   +fac1r *fac2 *cfactor( ixyz )
               if ( kimg == 2 ) then
                  MagField_constrain_local( i, 2, ixyz ) &
                       & = MagField_constrain_local( i, 2, ixyz ) &
                       &  -fac1i *fac2 *cfactor( ixyz )
               endif
            End Do

            if ( mype == 0 ) then
               fac2 = PAI4 /3.0d0 *rad1**3

               MagField_constrain_local( 1, 1, ixyz ) &
                    & = MagField_constrain_local( 1, 1, ixyz ) &
                    &   +fac2 *cfactor( ixyz )
            end if
         End Do
      End Do
!
      if ( noncol ) then
         Do i=ista_kngp, iend_kngp
            Do is=2, ndim_magmom
               vlhxc_l(i,1,is) = vlhxc_l(i,1,is) +MagField_constrain_local(i,1,is-1)
               vlhxc_l(i,2,is) = vlhxc_l(i,2,is) +MagField_constrain_local(i,2,is-1)
            End do
         End Do
      else
         Do i=ista_kngp, iend_kngp
            Do j=1, kimg
               vlhxc_l(i,j,1) = vlhxc_l(i,j,1) +MagField_constrain_local(i,j,1)
               vlhxc_l(i,j,2) = vlhxc_l(i,j,2) -MagField_constrain_local(i,j,1)
            End Do
         End Do
      endif

    end subroutine case_constraint_moment_local

    subroutine case_constraint_direc_local
      real(kind=DP) :: MagMom(3), cfactor(3), MagDirec(3)
      real(kind=DP) :: rad1, fac1r, fac1i, fac2
      real(kind=DP) :: VecG(3), normG, normG3, gr, d1
      real(kind=DP) :: c1, c2, c3, cnorm, cnorm2,  cnorm4

      real(kind=DP), parameter :: cnorm_lower_limit = 1.0D-10
!      real(kind=DP), parameter :: cnorm_lower_limit = 1.0D-4

      real(kind=DP), allocatable :: zfcos(:), zfsin(:)
      real(kind=DP), allocatable :: RhoMag_on_atom_mpi(:,:)

      integer :: i, j, ia, it, is, ixyz, ixyz_max, ist

      allocate(zfcos(ista_kngp:iend_kngp)); zfcos = 0.d0
      allocate(zfsin(ista_kngp:iend_kngp)); zfsin = 0.d0

      if ( noncol ) then
         ixyz_max = 3
      else
         ixyz_max = 1
      endif

      if ( allocated( MagField_constrain_local ) ) then
         deallocate( MagField_constrain_local )
      endif

      allocate( MagField_constrain_local(ista_kngp:iend_kngp, kimg, ixyz_max ))
      MagField_constrain_local = 0.0d0
!
      Do ia=1, natm
         it = ityp(ia)
         rad1 = rad_cov(ia)

         MagMom = 0.0d0
         Do ixyz=1, ixyz_max
            if ( noncol ) then
               MagMom( ixyz ) = RhoMag_on_atom( ia, ixyz +1 )
            else
               MagMom( ixyz ) = RhoMag_on_atom( ia, 1 ) -RhoMag_on_atom( ia, 2 )
            endif
         End do

         cnorm2 = 0.0d0
         Do ixyz=1, ixyz_max
            cnorm2 = cnorm2 + MagMom(ixyz)**2
         End do
         cnorm = sqrt( cnorm2 )

         if ( cnorm < cnorm_lower_limit ) cycle

         MagDirec = MagMom / cnorm
!

         cnorm4 = 0.0d0
         if ( mag_moment0_atoms_is_defined ) then
            Do ixyz=1, ixyz_max
               cnorm4 = cnorm4 +mag_moment0_atoms(ia,ixyz)**2
            End do
         else
            Do ixyz=1, ixyz_max
               cnorm4 = cnorm4 +mag_moment0_atomtyp(it,ixyz)**2
            End do
         endif
         cnorm4 = sqrt( cnorm4 )
!
         Do ixyz=1, ixyz_max
            if ( noncol ) then
               c1 = cnorm2 - MagMom(ixyz)**2;    c2 = cnorm**3
               if ( cnorm4 > 0.0 ) then
                  if ( mag_moment0_atoms_is_defined ) then
                     c3 = mag_constraint_lambda /univol *c1 /c2 &
                          &         *( MagDirec(ixyz) -mag_direction0_atoms(ia,ixyz) )
                  else
                     c3 = mag_constraint_lambda /univol *c1 /c2 &
                          &         *( MagDirec(ixyz) -mag_direction0_atomtyp(it,ixyz) )
                  endif
               else
                  c3 = 0.0d0
               endif
               cfactor( ixyz ) = c3
            else
               if ( cnorm4 > 0.0 ) then
                  if ( mag_moment0_atoms_is_defined ) then
                     if ( mag_moment0_atoms(ia,ixyz) > 0.0 ) then
                        c3 = 1.0d0
                     else
                        c3 = -1.0d0
                     endif
                  else
                     if ( mag_moment0_atomtyp(it,ixyz) > 0.0 ) then
                        c3 = 1.0d0
                     else
                        c3 = -1.0d0
                     endif
                  endif
                  cfactor( ixyz ) = mag_constraint_lambda /univol &
                       &               *( MagDirec(ixyz) -c3 )
               else
                  cfactor( ixyz ) = 0.0d0
               endif
            endif

!            write(2000+mype,*) "ia cfac = ", ia, cfactor(ixyz), MagDirec(1), c3
!            write(2100+mype,*) "ia  = ", ia, MagMom(ixyz), cnorm, cnorm2

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
               MagField_constrain_local( i, 1, ixyz ) &
                    & = MagField_constrain_local( i, 1, ixyz ) &
                    &   +fac1r *fac2 *cfactor( ixyz )
               if ( kimg == 2 ) then
                  MagField_constrain_local( i, 2, ixyz ) &
                       & = MagField_constrain_local( i, 2, ixyz ) &
                       &  -fac1i *fac2 *cfactor( ixyz )
               endif
            End Do

            if ( mype == 0 ) then
               fac2 = PAI4 /3.0d0 *rad1**3

               MagField_constrain_local( 1, 1, ixyz ) &
                    & = MagField_constrain_local( 1, 1, ixyz ) &
                    &   +fac2 *cfactor( ixyz )
            end if
         End Do
      End Do
!
      if ( noncol ) then
         Do i=ista_kngp, iend_kngp
            Do is=2, ndim_magmom
               vlhxc_l(i,1,is) = vlhxc_l(i,1,is) +MagField_constrain_local(i,1,is-1)
               vlhxc_l(i,2,is) = vlhxc_l(i,2,is) +MagField_constrain_local(i,2,is-1)
            End do
         End Do
      else
         Do i=ista_kngp, iend_kngp
            Do j=1, kimg
               vlhxc_l(i,j,1) = vlhxc_l(i,j,1) +MagField_constrain_local(i,j,1)
               vlhxc_l(i,j,2) = vlhxc_l(i,j,2) -MagField_constrain_local(i,j,1)
            End Do
         End Do
      endif

    end subroutine case_constraint_direc_local

  end subroutine m_ES_add_MagConstraintPot_chgql

  subroutine m_ES_add_MagConstraintPot_hsr
    if ( ndim_magmom == 1 ) return

    select case ( mag_constraint_type )
    case (MAG_MOMENT_VALS_OCCMAT)
       call case_constraint_moment_occmat
    case (MAG_MOMENT_DIREC_HARDPART)
       call case_constraint_direc_hardpart
    end select

  contains

    subroutine case_constraint_moment_occmat
      integer :: ia, ig, it, i, ip, ilp
      integer :: ilmt1, ilmt2, l1, l2, m1, m2, t1, t2
      integer :: ixyz, ixyz_max, is

      real(kind=DP) :: MagMom(3), cfactor(3), c1

      if ( noncol ) then
         ixyz_max = 3
      else
         ixyz_max = 1
      endif

      if ( allocated( MagField_constrain_hardpart ) ) then
         deallocate( MagField_constrain_hardpart )
      endif

      allocate( MagField_constrain_hardpart(nlmt,nlmt,natm,ixyz_max) )
      MagField_constrain_hardpart = 0.0d0
!
      do ia=1,natm
         it = ityp(ia)

         ig = iproj_group(ia)
         if(ig<1) cycle

         Do ixyz=1, ixyz_max
            if ( noncol ) then
               c1 = 0.0d0
               do i=1,num_proj_elems(ig)
                  ip = proj_group(i,ig)
                  Do m1=1, 2*proj_attribute(ip)%l +1
                     c1 = c1 +ommix(m1,m1,i,ia,ixyz+1)
                  End Do
               End do
               MagMom(ixyz) = c1
            else
               c1 = 0.0d0
               do i=1,num_proj_elems(ig)
                  ip = proj_group(i,ig)
                  Do m1=1, 2*proj_attribute(ip)%l +1
                     c1 = c1 +ommix(m1,m1,i,ia,1) -ommix(m1,m1,i,ia,2)
                  End Do
               End do
               MagMom(ixyz) = c1
            endif
            if ( mag_moment0_atoms_is_defined ) then
               cfactor(ixyz) = mag_constraint_lambda &
                    &               *( MagMom(ixyz) -mag_moment0_atoms(ia,ixyz) )
            else
               cfactor(ixyz) = mag_constraint_lambda &
                    &               *( MagMom(ixyz) -mag_moment0_atomtyp(it,ixyz) )
            endif
         End do

         Do ixyz=1, ixyz_max

            do i=1,num_proj_elems(ig)
               ip = proj_group(i,ig)
               it = proj_attribute(ip)%ityp
               ilp = proj_attribute(ip)%l+1

               do ilmt1 = 1, ilmt(it)
                  l1 = ltp(ilmt1,it); m1 = mtp(ilmt1,it); t1 = taup(ilmt1,it)
                  if ( l1 /= ilp ) cycle

                  do ilmt2 = 1, ilmt(it)
                     l2 = ltp(ilmt2,it); m2 = mtp(ilmt2,it); t2 = taup(ilmt2,it)
                     if( l2 /= ilp ) cycle

                     if ( m1 == m2 ) then
                        MagField_constrain_hardpart(ilmt1,ilmt2,ia,ixyz) &
                             & = MagField_constrain_hardpart(ilmt1,ilmt2,ia,ixyz) &
                             & + cfactor(ixyz) *prodphi(ip,t1,t2)
                     endif
                  end do
               end do
            end do
         End Do
      End do

      if ( noncol ) then
         Do is=2, ndim_magmom
            vlhxcQ(:,:,:,is) = vlhxcQ(:,:,:,is) &
                 &               +MagField_constrain_hardpart(:,:,:,is-1)
         End Do
      else
         vlhxcQ(:,:,:,1) = vlhxcQ(:,:,:,1) +MagField_constrain_hardpart(:,:,:,1)
         vlhxcQ(:,:,:,2) = vlhxcQ(:,:,:,2) -MagField_constrain_hardpart(:,:,:,1)
      end if

    end subroutine case_constraint_moment_occmat

    subroutine case_constraint_direc_hardpart        ! noncol only
      real(kind=DP) :: MagMom(3), MagDirec(3)
      real(kind=DP) :: c1, c2, cnorm, cnorm2
      integer :: ia, it, lmt1, lmt2, is, ixyz

      real(kind=DP), parameter :: cnorm_lower_limit = 1.0D-10

      if ( .not. noncol ) return

      if ( .not. allocated( MagField_constrain_hardpart ) ) then
         allocate( MagField_constrain_hardpart(nlmt,nlmt,natm,3 ) );
      endif
!
      MagField_constrain_local = 0.0d0
!
      Do ia=1, natm
         it = ityp(ia)
         Do ixyz=1, 3
            MagMom( ixyz ) = magmom_local_now( ia,ixyz )
         End do

         cnorm2 = 0.0d0
         Do ixyz=1, 3
            cnorm2 = cnorm2 + MagMom(ixyz)**2
         End do
         cnorm = sqrt( cnorm2 )

         if ( cnorm < cnorm_lower_limit ) then
            MagField_constrain_hardpart(:,:,ia,:) = 0.0d0
            cycle
         endif

         MagDirec = MagMom / cnorm

         Do lmt1=1, ilmt(it)
            Do lmt2=1, ilmt(it)

               Do ixyz=1, 3
                  c1 = cnorm2 - MagMom(ixyz)**2
                  c2 = cnorm**3
                  if ( mag_moment0_atoms_is_defined ) then
                     MagField_constrain_hardpart(lmt1,lmt2,ia,ixyz) = &
                          &  mag_constraint_lambda *c1 /c2 &
                          & *( MagDirec(ixyz) -mag_direction0_atoms(ia,ixyz) ) &
                          & *dl2p(lmt1,lmt2,1,it)             ! approx
                  else
                     MagField_constrain_hardpart(lmt1,lmt2,ia,ixyz) = &
                          &  mag_constraint_lambda *c1 /c2 &
                          & *( MagDirec(ixyz) -mag_direction0_atomtyp(it,ixyz) ) &
                          & *dl2p(lmt1,lmt2,1,it)             ! approx
                  endif
               End do
            End do
         End do

         Do lmt1=1, ilmt(it)
            Do lmt2=1, ilmt(it)
               Do is=2, ndim_magmom
                  vlhxcQ(lmt1,lmt2,ia,is) = vlhxcQ(lmt1,lmt2,ia,is) &
                       &       + MagField_constrain_hardpart(lmt1,lmt2,ia,is-1)
               End do
            End do
         End Do

      end Do

    end subroutine case_constraint_direc_hardpart

  end subroutine m_ES_add_MagConstraintPot_hsr

  subroutine m_ES_calc_MagConstraint_Energy( ene_double_counting, ene_mag_constraint )
    real(kind=DP), intent(out) :: ene_double_counting
    real(kind=DP), intent(out) :: ene_mag_constraint

    ene_double_counting = 0.0d0;      ene_mag_constraint = 0.0d0

    select case ( mag_constraint_type )

    case (MAG_MOMENT_VALS_GLOBAL)
       call case_constraint_moment_global
    case (MAG_MOMENT_DIREC_GLOBAL)
       call case_constraint_direc_global
    case (MAG_MOMENT_VALS_LOCAL)
       call case_constraint_moment_local
    case (MAG_MOMENT_DIREC_LOCAL)
       call case_constraint_direc_local

    case (MAG_MOMENT_VALS_OCCMAT)
       call case_constraint_moment_occmat
    case (MAG_MOMENT_DIREC_HARDPART)
       call case_constraint_direc_hardpart

    end select

  contains

    subroutine case_constraint_moment_global
      real(kind=DP) :: MagMom(3), c1, c2
      real(kind=DP) :: MagMomOld(3)
      integer :: is, ixyz

      if ( mype == 0 ) then

         if ( noncol ) then
            Do is=2, ndim_magmom
               MagMom(is-1) = chgq_l(  1,1,is ) *univol
!               MagMomOld(is-1) = chgqo_l(  1,1,is ) *univol
            End do

            c1 = 0.0d0;  c2 = 0.0d0
            Do ixyz=1, 3
               c1 = c1 + MagField_constrain_global(ixyz) *MagMom(ixyz)
               c2 = c2 + ( MagMom(ixyz) -mag_moment0_global(ixyz) )**2
            End do
         else
            MagMom(1) = ( chgq_l( 1,1,1 ) -chgq_l( 1,1,2 ) )*univol

            c1 = MagField_constrain_global(1) *MagMom(1)
            c2 = ( MagMom(1) -mag_moment0_global(1) )**2
         endif

         ene_double_counting = c1
         ene_mag_constraint = c2 *mag_constraint_lambda / 2.0d0
      endif

      if (npes > 1) then
         call mpi_bcast( ene_double_counting, 1, mpi_double_precision, &
              &          0, MPI_CommGroup, ierr )
         call mpi_bcast( ene_mag_constraint,  1, mpi_double_precision, &
              &          0, MPI_CommGroup, ierr )
      end if

    end subroutine case_constraint_moment_global

    subroutine case_constraint_direc_global
      real(kind=DP) :: MagMom(3), MagDirec(3), c1, c2
      integer :: is, ixyz

      real(kind=DP) :: cnorm, cnorm2
      real(kind=DP), parameter :: cnorm_lower_limit = 1.0D-10

      if ( mype == 0 ) then
         Do is=2, ndim_magmom
            MagMom(is-1) = chgq_l(  1,1,is ) *univol
         End do

         cnorm2 = 0.0d0
         Do ixyz=1, 3
            cnorm2 = cnorm2 + MagMom(ixyz)**2
         End do
         cnorm = sqrt( cnorm2 )

         if ( cnorm > cnorm_lower_limit ) then
            MagDirec = MagMom / cnorm

            c1 = 0.0d0;  c2 = 0.0d0
            Do ixyz=1, 3
               c1 = c1 + MagField_constrain_global(ixyz) *MagMom(ixyz)
               c2 = c2 + ( MagDirec(ixyz) -mag_direc0_global(ixyz) )**2
            End do

            ene_double_counting = c1
            ene_mag_constraint = c2 *mag_constraint_lambda /2.0d0
         endif
      endif

      if (npes > 1) then
         call mpi_bcast( ene_double_counting, 1, mpi_double_precision, &
              &          0, MPI_CommGroup, ierr )
         call mpi_bcast( ene_mag_constraint,  1, mpi_double_precision, &
              &          0, MPI_CommGroup, ierr )
      end if

    end subroutine case_constraint_direc_global

    subroutine case_constraint_direc_local
      real(kind=DP) :: MagMom(3), MagDirec(3), c1, c2, c3
      integer :: i, j, ia, it, is, ixyz

      real(kind=DP) :: cnorm, cnorm2, cnorm4
      real(kind=DP) :: csum1, csum2

      real(kind=DP), parameter :: cnorm_lower_limit = 1.0D-10

      csum1 = 0.0d0; csum2 = 0.0d0

      if ( noncol ) then
         Do i=ista_kngp, iend_kngp
            Do j=1, kimg
               Do ixyz=1, 3
                  csum1 = csum1 + MagField_constrain_local(i,j,ixyz)*chgq_l(i,j,ixyz+1)
               End Do
            End Do
         End Do

         Do ia=1, natm
            it = ityp(ia)
            Do ixyz=1, 3
               MagMom(ixyz) = RhoMag_on_Atom( ia,ixyz+1 )
            End do

            cnorm2 = 0.0d0
            Do ixyz=1, 3
               cnorm2 = cnorm2 + MagMom(ixyz)**2
            End do
            cnorm = sqrt( cnorm2 )

            cnorm4 = 0.0d0
            if ( mag_moment0_atoms_is_defined ) then
               Do ixyz=1, 3
                  cnorm4 = cnorm4 +mag_moment0_atoms(ia,ixyz)**2
               End do
            else
               Do ixyz=1, 3
                  cnorm4 = cnorm4 +mag_moment0_atomtyp(it,ixyz)**2
               End do
            endif
            cnorm4 = sqrt( cnorm4 )

            if ( cnorm4 > 0.0 .and. cnorm > cnorm_lower_limit ) then
               MagDirec = MagMom / cnorm
               if ( mag_moment0_atoms_is_defined ) then
                  Do ixyz=1, 3
                     csum2 = csum2 +( MagDirec(ixyz)-mag_direction0_atoms(ia,ixyz) )**2
                  End Do
               else
                  Do ixyz=1, 3
                     csum2 = csum2 +( MagDirec(ixyz)-mag_direction0_atomtyp(it,ixyz) )**2
                  End Do
               endif
            endif
         End do

      else
         Do i=ista_kngp, iend_kngp
            Do j=1, kimg
               csum1 = csum1 + MagField_constrain_local(i,j,1) &
                    &         *( chgq_l(i,j,1) -chgq_l(i,j,2) )
            End Do
         End Do

         Do ia=1, natm
            it = ityp(ia)
            MagMom(1) = RhoMag_on_Atom( ia,1 ) -RhoMag_on_Atom( ia,2 )

            cnorm2 = 0.0d0
            Do ixyz=1, 1
               cnorm2 = cnorm2 + MagMom(ixyz)**2
            End do
            cnorm = sqrt( cnorm2 )

            cnorm4 = 0.0d0
            if ( mag_moment0_atoms_is_defined ) then
               cnorm4 = cnorm4 +mag_moment0_atoms(ia,1)**2
            else
               cnorm4 = cnorm4 +mag_moment0_atomtyp(it,1)**2
            endif

            cnorm4 = sqrt( cnorm4 )

            if ( cnorm4 > 0.0 .and. cnorm > cnorm_lower_limit ) then
               MagDirec = MagMom / cnorm
               if ( mag_moment0_atoms_is_defined ) then
                  if ( mag_moment0_atoms(ia,1) > 0.0 ) then
                     c3 = 1.0d0
                  else
                     c3 = -1.0d0
                  endif
               else
                  if ( mag_moment0_atomtyp(it,1) > 0.0 ) then
                     c3 = 1.0d0
                  else
                     c3 = -1.0d0
                  endif
               endif
               csum2 = csum2 +( MagDirec(1) -c3 )**2
            endif
         End do
      endif

      if ( npes > 1 ) then
         call mpi_allreduce( csum1, c1, 1, mpi_double_precision, &
              &              mpi_sum, MPI_CommGroup, ierr )
         csum1 = c1
      endif

      ene_double_counting = csum1 *univol
      ene_mag_constraint = csum2 *mag_constraint_lambda/ 2.0d0

    end subroutine case_constraint_direc_local

    subroutine case_constraint_moment_local
      real(kind=DP) :: MagMom(3), c1, c2
      integer :: i, j, ia, it, is, ixyz

      real(kind=DP) :: csum1, csum2

      csum1 = 0.0d0; csum2 = 0.0d0

      if ( noncol ) then
         Do i=ista_kngp, iend_kngp
            Do j=1, kimg
               Do ixyz=1, 3
                  csum1 = csum1 + MagField_constrain_local(i,j,ixyz)*chgq_l(i,j,ixyz+1)
               End Do
            End Do
         End Do

         if ( mag_moment0_atoms_is_defined ) then
            Do ia=1, natm
               Do ixyz=1, 3
                  MagMom(ixyz) = RhoMag_on_Atom( ia,ixyz+1 )
                  csum2 = csum2 + ( MagMom(ixyz) -mag_moment0_atoms(ia,ixyz) )**2
               End do
            End do
         else
            Do ia=1, natm
               it = ityp(ia)
               Do ixyz=1, 3
                  MagMom(ixyz) = RhoMag_on_Atom( ia,ixyz+1 )
                  csum2 = csum2 + ( MagMom(ixyz) -mag_moment0_atomtyp(it,ixyz) )**2
               End do
            End do
         endif
      else
         Do i=ista_kngp, iend_kngp
            Do j=1, kimg
               csum1 = csum1 + MagField_constrain_local(i,j,1) &
                    &         *( chgq_l(i,j,1) -chgq_l(i,j,2) )
            End Do
         End Do

         if ( mag_moment0_atoms_is_defined ) then
            Do ia=1, natm
               MagMom(1) = RhoMag_on_Atom( ia,1 ) -RhoMag_on_Atom( ia,2 )
               csum2 = csum2 +( MagMom(1) -mag_moment0_atoms(ia,1) )**2
            End do
         else
            Do ia=1, natm
               it = ityp(ia)
               MagMom(1) = RhoMag_on_Atom( ia,1 ) -RhoMag_on_Atom( ia,2 )
               csum2 = csum2 +( MagMom(1) -mag_moment0_atomtyp(it,1) )**2
            End do
         endif
      endif

      if ( npes > 1 ) then
         call mpi_allreduce( csum1, c1, 1, mpi_double_precision, &
              &              mpi_sum, MPI_CommGroup, ierr )
         csum1 = c1
      endif

      ene_double_counting = csum1 *univol
      ene_mag_constraint = csum2 *mag_constraint_lambda/ 2.0d0

    end subroutine case_constraint_moment_local

    subroutine case_constraint_moment_occmat
      integer :: ia, it, ilmt1, ilmt2, m1
      integer :: i, ig, ip, ixyz, ixyz_max
      real(kind=DP) :: MagMom(3), c1, csum1, csum2

      if ( noncol ) then
         ixyz_max = 3
      else
         ixyz_max = 1
      endif

      csum1 = 0.0d0; csum2 = 0.0d0

      if ( noncol ) then
         Do ia=1, natm
            it = ityp(ia)
            Do ixyz=1, ixyz_max
               Do ilmt1=1, ilmt(it)
                  Do ilmt2=1, ilmt(it)
                     csum1 = csum1 +MagField_constrain_hardpart(ilmt1,ilmt2,ia,ixyz) &
                          &          *hsr(ia,ilmt1,ilmt2,ixyz+1)
                  End do
               End do
            End do
         End Do
      else
         Do ia=1, natm
            it = ityp(ia)
            Do ilmt1=1, ilmt(it)
               Do ilmt2=1, ilmt(it)
                  csum1 = csum1 +MagField_constrain_hardpart(ilmt1,ilmt2,ia,1) &
                       &          *( hsr(ia,ilmt1,ilmt2,1) -hsr(ia,ilmt1,ilmt2,2) )
               End Do
            End do
         End do
      endif

      Do ia=1, natm
         it = ityp(ia)

         ig = iproj_group(ia)
         if(ig<1) cycle

         if ( noncol ) then
            Do ixyz=1, ixyz_max
               c1 = 0.0d0
               do i=1,num_proj_elems(ig)
                  ip = proj_group(i,ig)
                  Do m1=1, 2*proj_attribute(ip)%l +1
                     c1 = c1 +om(m1,m1,i,ia,ixyz+1)
                  End Do
               End do
               MagMom(ixyz) = c1
            End do
         else
            c1 = 0.0d0
            do i=1,num_proj_elems(ig)
               ip = proj_group(i,ig)
               Do m1=1, 2*proj_attribute(ip)%l +1
                  c1 = c1 +om(m1,m1,i,ia,1) -om(m1,m1,i,ia,2)
               End Do
            End do
            MagMom(1) = c1
         endif

         if ( mag_moment0_atoms_is_defined ) then
            Do ixyz=1, ixyz_max
               csum2 = csum2 +( MagMom(ixyz) -mag_moment0_atoms(ia,ixyz) )**2
            End do
         else
            Do ixyz=1, ixyz_max
               csum2 = csum2 +( MagMom(ixyz) -mag_moment0_atomtyp(it,ixyz) )**2
            End do
         endif
      End Do
      ene_double_counting = csum1
      ene_mag_constraint = csum2 *mag_constraint_lambda/ 2.0d0

    end subroutine case_constraint_moment_occmat

    subroutine case_constraint_direc_hardpart
      real(kind=DP) :: MagMom(3), MagDirec(3), c1, c2
      integer :: ia, it, lmt1, lmt2, is, ixyz

      real(kind=DP) :: csum1, csum2
      real(kind=DP) :: cnorm, cnorm2
      real(kind=DP), parameter :: cnorm_lower_limit = 1.0D-10

      csum1 = 0.0d0; csum2 = 0.0d0

      Do ia=1, natm
         it = ityp(ia)
         Do ixyz=1, 3
            MagMom(ixyz) = magmom_local_now( ia,ixyz )
         End do

         cnorm2 = 0.0d0
         Do ixyz=1, 3
            cnorm2 = cnorm2 + MagMom(ixyz)**2
         End do
         cnorm = sqrt( cnorm2 )

         if ( cnorm > cnorm_lower_limit ) then
            MagDirec = MagMom / cnorm

            Do lmt1=1, ilmt(it)
               Do lmt2=1, ilmt(it)

                  Do ixyz=1, 3
                     csum1 = csum1 + MagField_constrain_hardpart(lmt1,lmt2,ia,ixyz) &
                          &          *hsr(ia,lmt1,lmt2,ixyz +1)
                  End do
               End do
            End do

            if ( mag_moment0_atoms_is_defined ) then
               Do ixyz=1, 3
                  csum2 = csum2 + ( MagDirec(ixyz)-mag_direction0_atoms(ia,ixyz) )**2
               End do
            else
               Do ixyz=1, 3
                  csum2 = csum2 + ( MagDirec(ixyz)-mag_direction0_atomtyp(it,ixyz) )**2
               End do
            endif
         endif

      End do

      ene_double_counting = csum1
      ene_mag_constraint = csum2 *mag_constraint_lambda/ 2.0d0

    end subroutine case_constraint_direc_hardpart

  end subroutine m_ES_calc_MagConstraint_Energy

! === KT_add === 2014/08/09
  subroutine m_ES_proj_magmom_G_quantz_axis
    integer :: i, j, k, is, ia
    real(kind=DP) :: ctmp1, ctmp2

! soft part
    Do k=1, kimg
       Do i=ista_kngp, iend_kngp
          ctmp1 = 0.0d0
          Do is=1, 3
             ctmp1 = ctmp1 +chgq_l(i,k,is+1) *Global_Quantz_Axis_Fixed(is)
          End do
          Do is=1, 3
             chgq_l(i,k,is+1) = ctmp1 *Global_Quantz_Axis_Fixed(is)
          End do
       End do
    End do
!
! hard part
    Do ia=1, natm
       Do i=1, nlmt
          Do j=1, nlmt
             ctmp1 = 0.0d0
             ctmp2 = 0.0d0
             Do is=1, 3
                ctmp1 = ctmp1 +hsr(ia,i,j,is+1) *Global_Quantz_Axis_Fixed(is)
                ctmp2 = ctmp2 +hsi(ia,i,j,is+1) *Global_Quantz_Axis_Fixed(is)
             End do
             Do is=1, 3
                hsr(ia,i,j,is+1) = ctmp1 *Global_Quantz_Axis_Fixed(is)
                hsi(ia,i,j,is+1) = ctmp2 *Global_Quantz_Axis_Fixed(is)
             End do
          End do
       End do
    End Do

  end subroutine m_ES_proj_magmom_G_quantz_axis
! ============== 2014/08/09

! === KT_add === 2014/09/26
  subroutine m_ES_set_magmom_zero
    integer :: is

    Do is=2, ndim_magmom
       chgq_l(:,:,is) = 0.0d0
       hsr(:,:,:,is) = 0.0d0;
!!!!       hsi(:,:,:,is) = 0.0d0
    End do
  end subroutine m_ES_set_magmom_zero
! ============== 2014/09/26

end module m_ES_Mag_Constraint
