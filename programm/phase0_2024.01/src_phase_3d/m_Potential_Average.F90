module m_Potential_Average
  use m_Control_Parameters,  only : nspin, noncol, kimg
  use m_Const_Parameters,  only : DP, PAI4, OFF, ON, FMAXVALLEN, Hartree
  use m_Parallelization,  only : ista_kngp, iend_kngp, npes, MPI_CommGroup, mype, &
       &                         myrank_g, mpi_chg_world
  use m_PlaneWaveBasisSet,  only : ngabc, gr_l, kgp, ngabc_kngp_l
  use m_Ionic_System,      only : natm, pos, ntyp, ityp, zfm3_l
  use m_Crystal_Structure,  only : rltv
  use m_PseudoPotential,      only : psc_l
  use m_Charge_Density,       only : chgq_l
  use m_XC_Potential,         only : vxc_l
  use m_Files,             only : nfout
  use mpi

  implicit none

!  include 'mpif.h'
  integer istatus(mpi_status_size)

! ---
  character(len("Postprocessing")),private,parameter :: tag_postprocessing    = "postprocessing"
  character(len("potential_average")),private,parameter :: &
       &         tag_potential_average = "potential_average"
  character(len("sw_calc_pot_avg")),private,parameter :: &
       &         tag_sw_calc_pot_avg = "sw_calc_pot_avg"
  character(len("sw_add_xc_pot")),private,parameter :: &
       &         tag_sw_add_xc_pot = "sw_add_xc_pot"
  character(len("cutoff_radius")),private,parameter :: &
       &         tag_cutoff_radius = "cutoff_radius"
  character(len("weight_function")),private,parameter :: &
       &         tag_weight_function = "weight_function"

  integer :: sw_calc_pot_avg = off
  integer :: sw_add_xc_pot = off
  real(kind=DP) :: cutoff_radius = 2.0d0
  real(kind=DP) :: rad_for_pot_avg

  integer :: weight_function = 1

contains

  subroutine m_PotAvg_read_param
    integer, parameter :: LEN_TITLE = 80
    integer :: f_selectBlock, f_getStringValue, f_getRealValue, f_getIntValue
    integer :: f_selectParentBlock, f_selectTop

    integer :: iret
    real(kind=DP) :: dret
    character(len=FMAXVALLEN) :: rstr

    iret = f_selectTop()

    if ( f_selectBlock( tag_postprocessing ) == 0 ) then

       if ( f_selectBlock( tag_potential_average ) == 0 ) then
          if ( f_getIntValue( tag_sw_calc_pot_avg, iret ) == 0 ) then
             sw_calc_pot_avg = iret
          endif
          if ( f_getIntValue( tag_sw_add_xc_pot, iret ) == 0 ) then
             sw_add_xc_pot = iret
          endif
          if ( f_getIntValue( tag_weight_function, iret ) == 0 ) then
             weight_function = iret
          endif
          if ( f_getRealValue( tag_cutoff_radius, dret, 'bohr' ) == 0 ) then
             cutoff_radius = dret
          endif
          if ( mype == 0 .and. sw_calc_pot_avg == ON ) then
             write(nfout,*) "=== Info: Potential Average ==="
             write(nfout,*) "** sw_add_xc_pot is ", sw_add_xc_pot
             write(nfout,*) "** weight_function is ", weight_function
             write(nfout,*) "** cutoff_radius is ", cutoff_radius
          endif
          iret = f_selectParentBlock()
       endif
       iret = f_selectParentBlock()
    endif

  end subroutine m_PotAvg_read_param

  subroutine m_PotAvg_calc_pot_on_atoms
    integer :: ismax
    real(kind=DP), allocatable :: pot_wk(:,:,:)
    real(kind=DP), allocatable :: pot_on_atoms(:,:)

    ismax = nspin
    if ( noncol ) ismax = 1

    allocate( pot_wk(ista_kngp:iend_kngp,kimg,ismax) );  pot_wk = 0.0d0
    allocate( pot_on_atoms(natm,ismax) );                pot_on_atoms = 0.0d0

    call set_pot_Gspace( ismax, pot_wk )
    call average_pot_on_atoms( ismax, pot_wk, pot_on_atoms )
    call print_pot_on_atoms( ismax, pot_on_atoms )

    deallocate( pot_wk )
    deallocate( pot_on_atoms )
  end subroutine m_PotAvg_calc_pot_on_atoms

  subroutine average_pot_on_atoms( ismax, pot_wk, pot_on_atoms )
    integer, intent(in) :: ismax
    real(kind=DP), intent(in) ::  pot_wk(ista_kngp:iend_kngp,kimg,ismax)
    real(kind=DP), intent(out) :: pot_on_atoms(natm,ismax)

    integer :: ia, ist, i, it, ierr, is
    real(kind=DP) :: csum, d1, fac1r, fac1i
    real(kind=DP) :: rad1, fac2, gr, normG, normG3, vecG(3), vol

    real(kind=DP), allocatable :: zfcos(:), zfsin(:)
    real(kind=DP), allocatable :: pot_on_atom(:,:), pot_mpi(:,:)

    allocate(zfcos(ista_kngp:iend_kngp)); zfcos = 0.d0
    allocate(zfsin(ista_kngp:iend_kngp)); zfsin = 0.d0

    Do ia=1, natm
       it = ityp(ia)
       call calc_phase2( natm, pos, ia, kgp, ngabc_kngp_l, ista_kngp, iend_kngp, &
            &            zfcos, zfsin )
       if ( ista_kngp == 1 ) then
          ist = 2
       else
          ist = ista_kngp
       endif
!
       select case (weight_function)
       case (0)       ! just on atoms
          Do is=1, ismax
             csum = 0.0d0
             Do i=ist, iend_kngp
                fac1r = zfcos(i);   fac1i = zfsin(i)
                if ( kimg == 1 ) then
                   d1 = pot_wk(i,1,is) *fac1r
                else
                   d1 = pot_wk(i,1,is) *fac1r -pot_wk(i,kimg,is)*fac1i
                endif
                csum = csum +d1
             End Do
             pot_on_atoms(ia,is) = csum
          End Do

       case (1)       ! rad-cut
!          rad1 = rad_for_pot_avg(it)
          rad1 = cutoff_radius

          Do is=1, ismax
             csum = 0.0d0
             Do i=ist, iend_kngp
                VecG(1) = rltv(1,1)*ngabc_kngp_l(i,1) +rltv(1,2)*ngabc_kngp_l(i,2) &
                     &   +rltv(1,3)*ngabc_kngp_l(i,3)
                VecG(2) = rltv(2,1)*ngabc_kngp_l(i,1) +rltv(2,2)*ngabc_kngp_l(i,2) &
                     &   +rltv(2,3)*ngabc_kngp_l(i,3)
                VecG(3) = rltv(3,1)*ngabc_kngp_l(i,1) +rltv(3,2)*ngabc_kngp_l(i,2) &
                     &   +rltv(3,3)*ngabc_kngp_l(i,3)
                normG = sqrt( VecG(1)**2 +VecG(2)**2 +VecG(3)**2 )
                normG3 = normG**3
                fac1r = zfcos(i);   fac1i = zfsin(i)
                if ( kimg == 1 ) then
                   d1 = pot_wk(i,1,is) *fac1r
                else
                   d1 = pot_wk(i,1,is) *fac1r -pot_wk(i,kimg,is)*fac1i
                endif
                gr = normG *rad1
                fac2 = -gr *cos(gr) + sin(gr)
                fac2 = fac2 *PAI4 /normG3
                csum = csum +d1 *fac2
             End Do
             if ( myrank_g == 0 ) then
                fac2 = PAI4 /3.0d0 *rad1**3
!!!                csum = csum +pot_wk(1,kimg,is) *fac2   ! ASMS 2021/03/26
                csum = csum +pot_wk(1,1,is) *fac2
             endif

             vol = PAI4 /3.0 *rad1**3
             pot_on_atoms(ia,is) = csum /vol
          End do

       case (2)       ! gaussian
!          rad1 = rad_for_pot_avg(it)
          rad1 = cutoff_radius
          Do is=1, ismax
             csum = 0.0d0
             Do i=ista_kngp, iend_kngp
                VecG(1) = rltv(1,1)*ngabc_kngp_l(i,1) +rltv(1,2)*ngabc_kngp_l(i,2) &
                     &   +rltv(1,3)*ngabc_kngp_l(i,3)
                VecG(2) = rltv(2,1)*ngabc_kngp_l(i,1) +rltv(2,2)*ngabc_kngp_l(i,2) &
                     &   +rltv(2,3)*ngabc_kngp_l(i,3)
                VecG(3) = rltv(3,1)*ngabc_kngp_l(i,1) +rltv(3,2)*ngabc_kngp_l(i,2) &
                     &   +rltv(3,3)*ngabc_kngp_l(i,3)
                normG = sqrt( VecG(1)**2 +VecG(2)**2 +VecG(3)**2 )
                fac1r = zfcos(i);   fac1i = zfsin(i)
                if ( kimg == 1 ) then
                   d1 = pot_wk(i,1,is) *fac1r
                else
                   d1 = pot_wk(i,1,is) *fac1r -pot_wk(i,kimg,is)*fac1i
                endif
                gr = normG *rad1
                fac2 = exp( -gr**2 /2.0d0 )
                csum = csum +d1 *fac2
             End Do
             pot_on_atoms(ia,is) = csum
          End Do
       end select
    End Do

    deallocate( zfcos, zfsin )

    if ( npes > 1 ) then
       allocate( pot_mpi(natm,ismax) ); pot_mpi = 0.0d0
       call mpi_allreduce( pot_on_atoms, pot_mpi, natm*ismax, &
            &              mpi_double_precision, mpi_sum, mpi_chg_world, ierr )
       pot_on_atoms = pot_mpi
       deallocate( pot_mpi )
    endif
  end subroutine average_pot_on_atoms

  subroutine print_pot_on_atoms( ismax, pot_on_atoms )
    integer, intent(in) :: ismax
    real(kind=DP), intent(in) :: pot_on_atoms(natm,ismax)

    integer :: ia, lun
    character*72 file1

    if ( mype /= 0 ) return

    lun = 9000
    file1 = "potential_on_atoms.data"
    open( unit=lun, file=file1, status="unknown", form="formatted" )

    if ( ismax == 1 ) then
       write(lun,'(A)') "# Potential on atoms "
       write(lun,'(A)') "#   id     pot (eV)"
       Do ia=1, natm
          write(lun,'(I5,1F15.6)') ia, pot_on_atoms(ia,1)*Hartree
       End Do
    else
       write(lun,*) "# Potential on atoms "
       write(lun,*) "#   id     pot_up (eV)    pot_dn(eV)"
       Do ia=1, natm
          write(lun,'(I5,2F15.6)') ia, pot_on_atoms(ia,1:2)*Hartree
       End Do
    endif
    close(lun)

  end subroutine print_pot_on_atoms

  subroutine set_pot_Gspace( ismax, pot_wk )
    integer, intent(in) :: ismax
    real(kind=DP), intent(out) :: pot_wk(ista_kngp:iend_kngp,kimg,ismax)
    integer :: ist, is, ik, i, it
    real(kind=DP) :: c1

    ist = ista_kngp
    if(ist == 1) ist = 2

    do ik = 1, kimg
       if ( nspin == 1 .or. noncol ) then
          do i = ist, iend_kngp  !for mpi
             c1 = PAI4*chgq_l(i,ik,1)/gr_l(i)**2
             pot_wk(i,ik,1) = c1
          end do
       else if (nspin == 2) then
          do i = ist, iend_kngp  !for mpi
             c1 = PAI4*(chgq_l(i,ik,1)+chgq_l(i,ik,2))/gr_l(i)**2
             pot_wk(i,ik,1) = c1
          end do
          pot_wk(:,ik,2) = pot_wk(:,ik,1)
       endif
       do is=1, ismax
          do it = 1, ntyp
             do i = ista_kngp, iend_kngp  !for mpi
                pot_wk(i,ik,is) = pot_wk(i,ik,is) +psc_l(i,it)*zfm3_l(i,it,ik)
             end do
          end do
       end do
!
       if ( sw_add_xc_pot == ON ) then
          do is=1, ismax
             do i = ista_kngp, iend_kngp  !for mpi
                pot_wk(i,ik,is) = pot_wk(i,ik,is) +vxc_l(i,ik,is)
             end do
          end do
       endif
    end do
  end subroutine set_pot_Gspace

end module m_Potential_Average
