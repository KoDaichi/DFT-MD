module m_ES_ChargeState
  use m_Const_Parameters,   only : DP, PAI, PAI2, PAI4, FMAXVALLEN, NOCONV, ON, OFF, &
       &                           LOWER, Hartree, Bohr

  use m_Control_Parameters, only : kimg, noncol, nspin, iprichgdefect, &
       &                           sw_calc_extfnv_correction
  use m_Ionic_System,       only : natm, natm2, pos, ntyp, zfm3_l, ityp, iatomn, &
       &                           m_IS_pack_all_ions_in_uc, iatom
  use m_Crystal_Structure,  only : altv, rltv, univol
  use m_Parallelization,    only : ista_kngp, iend_kngp, ista_atm, iend_atm, &
       &                            mype, npes, MPI_CommGroup, ierr
  use m_PseudoPotential,    only : psc_l, ival
  use m_Charge_Density,     only : chgq_l
  use m_PlaneWaveBasisSet,  only : kgp, kg, gr_l, ngabc, kgp_reduced, igfp_l
  use m_Electronic_Structure, only : totch

  use m_FFT,    only : m_FFT_CD_inverse0, nfftp_nonpara, fft_box_size_CD_nonpara, &
       &               fft_box_size_CD

  use m_Files,  only : nfout, m_Files_open_electrosta_pot, &
       &                      m_Files_close_electrosta_pot, nf_elecpot_bin, &
       &                      m_Files_open_electrosta_pot_ref, &
       &                      m_Files_close_electrosta_pot_ref, nf_elecpot_bin_ref
  use m_CD_Mag_Moment,  only : rad_cov

  use m_PlaneWaveBasisSet, only : igfp_nonpara
  use mpi


  implicit none

!  include 'mpif.h'
  integer istatus(mpi_status_size)

! -------------------------------------------------------

  character(len("Postprocessing")),private,parameter :: tag_postprocessing    = "postprocessing"
  character(len("charged_defect")),private,parameter :: &
       &         tag_charged_defect    = "charged_defect"

  character(len("dielectric_constant")),private,parameter :: &
       & tag_dielectric_const = "dielectric_constant"
  character(len("exx")),private,parameter :: tag_exx = "exx"
  character(len("eyy")),private,parameter :: tag_eyy = "eyy"
  character(len("ezz")),private,parameter :: tag_ezz = "ezz"
  character(len("exy")),private,parameter :: tag_exy = "exy"
  character(len("eyz")),private,parameter :: tag_eyz = "eyz"
  character(len("ezx")),private,parameter :: tag_ezx = "ezx"

  character(len("position")),private,parameter :: tag_position  = "position"
  character(len("atom_id")),private,parameter :: tag_atom_id = "atom_id"
  character(len("x")),private,parameter :: tag_x = "x"
  character(len("y")),private,parameter :: tag_y = "y"
  character(len("z")),private,parameter :: tag_z = "z"

  integer :: atomid_defect
  logical :: atomid_is_set = .false.
  real(kind=DP) :: q, eps(3,3), pos_defect(3)

  character(len("electrostatic_potential")),private,parameter ::  &
       &       tag_electrostatic_potential  = "electrostatic_potential"
  character(len("electrostatic_pot")),private,parameter ::  &
       &       tag_electrostatic_pot  = "electrostatic_pot"
  character(len("sw_write_electrostatic_pot")),private,parameter :: &
       &        tag_sw_write_electrostatic_pot = "sw_write_electrostatic_pot"
  integer :: sw_write_electrostatic_pot = OFF

  character(len("format")),private,parameter :: tag_pot_format  = "format"
  character(len("cube")),private,parameter :: tag_cube  = "cube"
  character(len("atoms")),private,parameter :: tag_atoms = "atoms"

  character(len("unit")),private,parameter :: tag_pot_unit  = "unit"
  character(len("Hartree")),private,parameter :: tag_Hartree  = "hartree"
  character(len("Rydberg")),private,parameter :: tag_Rydberg  = "rydberg"
  character(len("eV")),private,parameter :: tag_eV  = "ev"

  integer, parameter :: CUBE = 1, ATOMS = 2
  integer, parameter :: in_Hartree = 0, in_Rydberg = 1, in_eV = 2
  integer :: pot_format = CUBE, pot_unit = in_eV

  character(len("correction")),private,parameter ::  tag_correction = "correction"

  character(len("sw_calc_extfnv_correction")),private,parameter :: &
       &        tag_sw_calc_extfnv_correction = "sw_calc_extfnv_correction"
!  integer :: sw_calc_extfnv_correction = OFF

  real(kind=DP), allocatable :: potdiff_on_grid(:,:,:), potdiff_on_atoms(:)
  real(kind=DP), allocatable :: Vpc_on_grid(:,:,:), Vpc_on_atoms(:)

contains

  subroutine m_ESCS_read_param
    integer, parameter :: LEN_TITLE = 80
    integer :: f_selectBlock, f_getStringValue, f_getRealValue, f_getIntValue
    integer :: f_selectParentBlock, f_selectTop

    integer :: iret
    real(kind=DP) :: dret
    character(len=FMAXVALLEN) :: rstr

    iret = f_selectTop()

    if ( f_selectBlock( tag_postprocessing ) == 0 ) then

       if ( f_selectBlock( tag_electrostatic_pot ) == 0 &
            & .or. f_selectBlock( tag_electrostatic_potential ) == 0 ) then
          if ( f_getIntValue( tag_sw_write_electrostatic_pot, iret ) == 0 ) then
             sw_write_electrostatic_pot = iret
             if ( mype == 0 ) &
                  &  write(nfout,*) "** sw_write_electrostatic_pot is ", iret
          endif

          iret = f_getStringValue( tag_pot_format,rstr,LOWER )
          if ( rstr == tag_cube ) then
             pot_format = CUBE
          else if ( rstr == tag_atoms ) then
             pot_format = ATOMS
          endif

          iret = f_getStringValue( tag_pot_unit,rstr,LOWER )
          if ( rstr == tag_Hartree ) then
             pot_unit = in_Hartree
          else if ( rstr == tag_Rydberg ) then
             pot_unit = in_Rydberg
          else if ( rstr == tag_eV ) then
             pot_unit = in_eV
          endif

          if ( mype == 0 ) then
             write(nfout,*) "** pot_format is ", pot_format
             write(nfout,*) "** pot_unit is ", pot_unit
          endif
          iret = f_selectParentBlock()
       endif

       if ( f_selectBlock( tag_charged_defect ) == 0 ) then

          if ( mype == 0 ) write(nfout,*) "=== Info: Charged Defect Analysis ==="

          if ( f_selectBlock( tag_correction ) == 0 ) then
             if ( f_getIntValue( tag_sw_calc_extfnv_correction, iret ) == 0 ) then
                sw_calc_extfnv_correction = iret
                if ( mype == 0 ) &
                     &  write(nfout,*) "** sw_calc_extfnv_correction is ", iret
             endif

             if ( f_selectBlock( tag_dielectric_const ) == 0 ) then
                if ( f_getRealValue( tag_exx, dret, '') == 0 ) eps(1,1) = dret
                if ( f_getRealValue( tag_eyy, dret, '') == 0 ) eps(2,2) = dret
                if ( f_getRealValue( tag_ezz, dret, '') == 0 ) eps(3,3) = dret
                if ( f_getRealValue( tag_exy, dret, '') == 0 ) eps(1,2) = dret
                if ( f_getRealValue( tag_eyz, dret, '') == 0 ) eps(2,3) = dret
                if ( f_getRealValue( tag_ezx, dret, '') == 0 ) eps(3,1) = dret
                eps(2,1) = eps(1,2);  eps(3,2) = eps(2,3);    eps(1,3) = eps(3,1)
                iret = f_selectParentBlock()
             else
                eps = 0.0d0;
                eps(1,1) = 1.0d0;  eps(2,2) = 1.0d0;  eps(3,3) = 1.0d0
             endif

             if ( f_selectBlock( tag_position ) == 0 ) then
                if ( f_getIntValue( tag_atom_id, iret ) == 0 ) then
                   atomid_defect = iret;      atomid_is_set = .true.
                endif
                if ( f_getRealValue( tag_x, dret, '') == 0 ) pos_defect(1) = dret
                if ( f_getRealValue( tag_y, dret, '') == 0 ) pos_defect(2) = dret
                if ( f_getRealValue( tag_z, dret, '') == 0 ) pos_defect(3) = dret
                iret = f_selectParentBlock()
             else
                pos_defect = 0.0d0       ! origin
             endif
! bin nomi saiyou ?

             if ( mype == 0 ) then
                if ( sw_calc_extfnv_correction == ON ) then
                   if ( atomid_is_set ) then
                      write(nfout,*) "** defect atom_id      :", atomid_defect
                   else
                      write(nfout,*) "** defect position (x) :", pos_defect(1)
                      write(nfout,*) "**                 (y) :", pos_defect(2)
                      write(nfout,*) "**                 (z) :", pos_defect(3)
                   endif
                   write(nfout,*) "** dielectric tensor (xx) :", eps(1,1)
                   write(nfout,*) "**                   (yy) :", eps(2,2)
                   write(nfout,*) "**                   (zz) :", eps(3,3)
                endif
             endif

             iret = f_selectParentBlock()
          endif

!          q = additional_charge

          iret = f_selectParentBlock()
       end if
       iret = f_selectParentBlock()

    end if

  end subroutine m_ESCS_read_param

! =================
  subroutine m_ESCS_wd_electrostatic_pot
    integer :: nlp, nmp, nnp, lun
    character*72 :: file1
    real(kind=DP), allocatable :: pot_now(:,:), pot_on_grid(:,:,:), pot_on_atoms(:)

    allocate( pot_now(ista_kngp:iend_kngp,kimg) );   pot_now = 0.0d0
    call set_electrocstatic_pot_Gspace( pot_now )

!    lun = 1000;  file1 = "./Potential.bin"
!    call wd_pot_bin( lun, file1, pot_now )

    call m_Files_open_electrosta_pot
    lun = nf_elecpot_bin
    call wd_pot_bin( lun, pot_now )

    call m_Files_close_electrosta_pot

    if ( pot_format == CUBE ) then
       nlp = fft_box_size_CD(1,1)
       nmp = fft_box_size_CD(2,1)
       nnp = fft_box_size_CD(3,1)
       allocate( pot_on_grid(nnp,nmp,nlp) ); pot_on_grid = 0.0d0
       call set_electrostatic_pot_on_grid( pot_now, nlp, nmp, nnp, pot_on_grid )

       lun = 1000;  file1 = "./electrostatic_pot.cube"
       call wd_pot_cube( lun, file1, nlp, nmp, nnp, pot_on_grid, &
            &            "Local+Hartree potential" )
       deallocate( pot_on_grid )

    else if ( pot_format == ATOMS ) then
       allocate( pot_on_atoms(natm) ); pot_on_atoms = 0.0d0
       call set_electrostatic_pot_on_atoms( pot_now, pot_on_atoms )

       lun = 1000;  file1 = "./electrostatic_pot.atoms"
       call wd_pot_on_atoms( lun, file1, pot_on_atoms, &
            &                "Local+Hartree potential" )
       deallocate( pot_on_atoms )

    endif

    deallocate( pot_now )

  end subroutine m_ESCS_wd_electrostatic_pot

  subroutine m_ESCS_wd_Madelung_pot      ! assuming Madelung pot is already calculated.
    integer :: nlp, nmp, nnp, lun
    character*72 :: file1

    if ( pot_format == CUBE ) then
       nlp = fft_box_size_CD(1,1)
       nmp = fft_box_size_CD(2,1)
       nnp = fft_box_size_CD(3,1)

       lun = 1000;  file1 = "./Madelung_pot.cube"
       call wd_pot_cube( lun, file1, nlp, nmp, nnp, Vpc_on_grid, &
            &            "Madelung pot." )

    else if ( pot_format == ATOMS ) then
       lun = 1000;  file1 = "./Madelung_pot.atoms"
       call wd_pot_on_atoms( lun, file1, Vpc_on_atoms, "Madelung pot." )
    endif
  end subroutine m_ESCS_wd_Madelung_pot

  subroutine set_electrocstatic_pot_Gspace( pot_wk )
    real(kind=DP), intent(out) :: pot_wk(ista_kngp:iend_kngp,kimg)
    integer :: ist, is, ik, i, it

    ist = ista_kngp
    if(ist == 1) ist = 2

    do ik = 1, kimg
       if ( noncol .or. nspin == 1 ) then
          do i = ist, iend_kngp  !for mpi
             pot_wk(i,ik) = PAI4*chgq_l(i,ik,1)/gr_l(i)**2
          end do
       else if (nspin == 2) then
          do i = ist, iend_kngp  !for mpi
             pot_wk(i,ik) = PAI4*(chgq_l(i,ik,1)+chgq_l(i,ik,2))/gr_l(i)**2
          end do
       endif
       do it = 1, ntyp
          do i = ista_kngp, iend_kngp  !for mpi
             pot_wk(i,ik) = pot_wk(i,ik)+psc_l(i,it)*zfm3_l(i,it,ik)
          end do
       end do
    end do
  end subroutine set_electrocstatic_pot_Gspace

  subroutine set_electrostatic_pot_on_grid( pot_wk, nlp, nmp, nnp, pot_on_grid )
    integer, intent(in) :: nlp, nmp, nnp
    real(kind=DP), intent(in) ::  pot_wk(ista_kngp:iend_kngp,kimg)
    real(kind=DP), intent(out) :: pot_on_grid(nnp,nmp,nlp)

    real(kind=DP), allocatable, dimension(:)       :: afft, bfft
    real(kind=DP), allocatable, dimension(:)       :: afft_mpi1

    allocate( afft(nfftp_nonpara) );      afft = 0.0d0
    allocate( afft_mpi1(nfftp_nonpara) ); afft_mpi1 = 0.0d0

    call map_density_to_fft_box
    call m_FFT_CD_inverse0(nfout,afft)
    call set_data_on_grid( nlp, nmp, nnp, pot_on_grid )

    deallocate( afft );   deallocate( afft_mpi1 )

  contains

    subroutine map_density_to_fft_box
      integer :: j, i, ip, it

      afft_mpi1 = 0.d0
      do j = 1, kimg
         do i = ista_kngp, iend_kngp
            if(kgp_reduced < i) cycle
#ifdef _MPIFFT_
            ip = (igfp_nonpara(i)-1)*kimg + j
#else
            ip = (igfp_l(i)-1)*kimg + j
#endif
            afft_mpi1(ip) = afft_mpi1(ip) +pot_wk(i,j)
         end do
      end do

      if(npes >= 2) then
         call mpi_allreduce(afft_mpi1,afft,nfftp_nonpara,mpi_double_precision &
              &  ,mpi_sum,MPI_CommGroup,ierr)
      else
         afft = afft_mpi1
      end if
    end subroutine map_density_to_fft_box

    subroutine set_data_on_grid( nlp, nmp, nnp, wkchr )
      integer, intent(in) :: nlp, nmp, nnp
      real(kind=DP), intent(out) :: wkchr(nnp,nmp,nlp)

      integer :: idp, mmp, up_down
      integer :: inew, jnew, knew, i, j, k, ip, nlphf
      real(kind=DP) :: s1, s2, sratio
      integer, parameter :: UP = 1 , DOWN = 2

      idp = fft_box_size_CD_nonpara(1,0)
      mmp = fft_box_size_CD_nonpara(2,0)

      if(kimg == 1) then
         nlphf = idp/2
      else
         nlphf = idp
      end if

! -- checking of symmetric or anti-symmetric about the charge densities
      s1 = 0.d0; s2 = 0.d0
      do ip = 1, nlphf*mmp*nnp,2
         s1 = s1 +dabs(afft(ip));     s2 = s2 +dabs(afft(ip+1))
      end do
!      if(ipri >= 2) write(nfout,'(" s1 = ",d14.6," s2 = ",d14.6)') s1, s2
      if(s1 > s2) then
!         if(ipri >= 2) &
!              & write(nfout,*) ' The function of pot. density is analyzed as symmetric.'
         up_down = UP;      sratio = s2/s1
      else
!         if(ipri >= 2) &
!              & write(nfout,*) ' The function of pot. density is analyzed as anti-symmetric.'
         up_down = DOWN;     sratio = s1/s2
      endif
!      if(ipri >= 2) write(nfout,*) ' !ratio = ', sratio
!      if(ipri >= 2) write(nfout,9001) nlp*nmp*nnp, nlp, nmp, nnp
!9001  format(' Pot. DENSITY NE = ',i8,'(',3i5,')')

!      if(ipri >= 2) write(nfout,*) ' !D FFT cube mapping start'
      do i = 1, nmp
         do j = 1, nnp
            do k = 1, nlp
               if(kimg == 1 .and. k > nlphf) then
                  knew = idp - k;   jnew = nnp+2 -j;   inew = nmp+2 -i
                  if(jnew > nnp) jnew = jnew -nnp
                  if(inew > nmp) inew = inew -nmp
               else
                  knew = k; jnew = j; inew = i
               end if
               ip = nlphf*mmp*(jnew-1) +nlphf*(inew-1) +knew
               wkchr(j,i,k) = afft(ip*2-2+up_down)
            end do
         end do
      end do
    end subroutine set_data_on_grid

  end subroutine set_electrostatic_pot_on_grid

  subroutine wd_pot_cube( lun, file1, nlp, nmp, nnp, wkchr, comment )
    integer, intent(in) :: lun, nlp, nmp, nnp
    character*72, intent(in) :: file1
    character*(*), intent(in) :: comment
    real(kind=DP), intent(in) :: wkchr(nnp,nmp,nlp)

    integer :: natm3, i, j, k, m
    integer :: boxsize(3), factor
    real(kind=DP) :: x, y, z, trans(3)

    integer, allocatable :: ityp_full(:)
    real(kind=DP), allocatable ::  cps_full(:,:) ,pos_full(:,:)

    if (mype /= 0) return

    open( lun, file=file1, status="unknown", form="formatted" )

    write(lun,'(A)') "Calculated by phase"
    if ( pot_unit == in_Rydberg ) then
       write(lun,'(A,A)') trim(adjustl(comment)), " in Rydberg"
    else if ( pot_unit == in_eV ) then
       write(lun,'(A,A)') trim(adjustl(comment)), " in eV"
    else
       write(lun,'(A,A)') trim(adjustl(comment)), " in Hartree"
    endif

    allocate(cps_full(natm2,3))
    allocate(ityp_full(natm2))
    call m_IS_pack_all_ions_in_uc(ityp_full,cps_full)

    natm3 = natm2

    x = 0.d0; y = 0.d0; z = 0.d0
    write(lun,'(i6,3f15.4)') natm3, x,y,z
    do i=1,3
       boxsize(i) = fft_box_size_CD(i,1)
    enddo
    do i = 1, 3
       write(lun,'(i6,3f12.6)') boxsize(i), altv(1:3,i)/dble(fft_box_size_CD(i,1))
    end do

    trans = 0.d0
    do i = 1, natm2
       m = ityp_full(i)
       cps_full(i,:) = cps_full(i,:)+trans(:)
       write(lun,'(i6,1x,f10.6,3(1x,f12.6))') &
            &        nint(iatomn(m)), iatomn(m), cps_full(i,1:3)
    end do
    deallocate(ityp_full,cps_full)

    factor = 1.0d0
    if ( pot_unit == in_Rydberg ) then
       factor = 2.0d0
    else if ( pot_unit == in_eV ) then
       factor = Hartree
    endif

    do i = 1, nlp
       do j = 1, nmp
          write(lun,'(6e13.5)') (wkchr(k,j,i)*factor,k=1,nnp)
       end do
    end do

  end subroutine wd_pot_cube

  subroutine set_electrostatic_pot_on_atoms( pot_wk, pot_on_atoms )
    real(kind=DP), intent(in) ::  pot_wk(ista_kngp:iend_kngp,kimg)
    real(kind=DP), intent(out) :: pot_on_atoms(natm)

    integer :: ia, ist, i, it
    real(kind=DP) :: csum, d1, fac1r, fac1i
    real(kind=DP) :: rad1, fac2, gr, normG, normG3, vecG(3), vol

    real(kind=DP), allocatable :: zfcos(:), zfsin(:)
    real(kind=DP), allocatable :: pot_on_atom(:), pot_mpi(:)

    allocate(zfcos(ista_kngp:iend_kngp)); zfcos = 0.d0
    allocate(zfsin(ista_kngp:iend_kngp)); zfsin = 0.d0

    Do ia=1, natm
       it = ityp(ia)
       rad1 = rad_cov(ia)
       call calc_phase2( natm, pos, ia, kgp, ngabc, ista_kngp, iend_kngp, &
            &            zfcos, zfsin )
       if ( ista_kngp == 1 ) then
          ist = 2
       else
          ist = ista_kngp
       endif

       csum = 0.0d0
       Do i=ist, iend_kngp
          VecG(1) = rltv(1,1)*ngabc(i,1) +rltv(1,2)*ngabc(i,2) +rltv(1,3)*ngabc(i,3)
          VecG(2) = rltv(2,1)*ngabc(i,1) +rltv(2,2)*ngabc(i,2) +rltv(2,3)*ngabc(i,3)
          VecG(3) = rltv(3,1)*ngabc(i,1) +rltv(3,2)*ngabc(i,2) +rltv(3,3)*ngabc(i,3)

          normG = sqrt( VecG(1)**2 +VecG(2)**2 +VecG(3)**2 )
          normG3 = normG**3

          fac1r = zfcos(i);   fac1i = zfsin(i)

          gr = normG *rad1
          fac2 = -gr *cos(gr) + sin(gr)
          fac2 = fac2 *PAI4 /normG3

          if ( kimg == 1 ) then
             d1 = pot_wk(i,1) *fac1r
          else
             d1 = pot_wk(i,1) *fac1r -pot_wk(i,kimg)*fac1i
          endif
          csum = csum +d1 *fac2
       End Do
       vol = PAI4 /3.0 *rad1**3
       Pot_on_atoms(ia) = csum /vol
    End Do

    deallocate( zfcos, zfsin )

    if ( npes > 1 ) then
       allocate( pot_mpi(natm) ); pot_mpi = 0.0d0
       call mpi_allreduce( pot_on_atoms, pot_mpi, natm, &
            &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
       pot_on_atoms = pot_mpi
       deallocate( pot_mpi )
    endif
  end subroutine set_electrostatic_pot_on_atoms

  subroutine wd_pot_on_atoms( lun, file1, pot_on_atoms, comment )
    integer, intent(in) :: lun
    character*72, intent(in) :: file1
    character*(*), intent(in) :: comment
    real(kind=DP), intent(in) :: pot_on_atoms(natm)

    integer :: i
    real(kind=DP) :: factor

    factor = 1.0d0
    if ( pot_unit == in_Rydberg ) then
       factor = 2.0d0
    else if ( pot_unit == in_eV ) then
       factor = Hartree
    endif

    if ( mype == 0 ) then
       open( unit=lun, file=file1, status="unknown", form="formatted" )
       if ( pot_unit == in_Rydberg ) then
          write(lun,'(A,A)') trim(adjustl(comment)), " in Rydberg"
       else if ( pot_unit == in_eV ) then
          write(lun,'(A,A)') trim(adjustl(comment)), " in eV"
       else
          write(lun,'(A,A)') trim(adjustl(comment)), " in Hartree"
       endif
       Do i=1, natm
          write(lun,'(I5,4F20.10)') i, pos(i,1:3),pot_on_atoms(i) *factor
       End Do
       close(lun)
    endif

  end subroutine wd_pot_on_atoms
! =================

  subroutine wd_pot_bin( lun, pot_wk )
    integer, intent(in) :: lun
    real(kind=DP), intent(in) :: pot_wk(ista_kngp:iend_kngp,kimg)

    integer :: i, ik, ierr
    real(kind=DP), allocatable :: pot_mpi1(:,:), pot_mpi2(:,:)

    if ( npes > 1 ) then
       allocate( pot_mpi1(kgp,kimg) );  pot_mpi1 = 0.0d0
       allocate( pot_mpi2(kgp,kimg) );  pot_mpi2 = 0.0d0
       do ik = 1, kimg
          do i = ista_kngp, iend_kngp
             pot_mpi1(i,ik) = pot_wk(i,ik)
          end do
       end do
       call mpi_allreduce( pot_mpi1, pot_mpi2, kgp*kimg, &
            &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
       if (mype == 0) then
          write(lun) pot_mpi2
       endif
       deallocate(pot_mpi1);  deallocate(pot_mpi2)
    else
       write(lun) pot_wk
    endif
  end subroutine wd_pot_bin

  subroutine rd_pot_bin( lun, pot_wk )
    integer, intent(in) :: lun
    real(kind=DP), intent(out) :: pot_wk(ista_kngp:iend_kngp,kimg)

    integer :: i, ik, ierr
    real(kind=DP), allocatable :: pot_mpi1(:,:), pot_mpi2(:,:)

    if ( npes > 1 ) then
       allocate( pot_mpi1(kgp,kimg) );  pot_mpi1 = 0.0d0
       if (mype == 0) then
          read(lun) pot_mpi1
       endif
       call mpi_bcast( pot_mpi1, kgp*kimg, mpi_double_precision, 0, MPI_CommGroup, ierr )
       do ik = 1, kimg
          do i = ista_kngp, iend_kngp
             pot_wk(i,ik) = pot_mpi1(i,ik)
          end do
       end do
       deallocate(pot_mpi1)
    else
       if (mype == 0) then
          read(lun) pot_wk
       endif
    endif
  end subroutine rd_pot_bin

! ===========================
  subroutine m_ESCS_calc_pot_diff
    integer :: nlp, nmp, nnp, lun
    character*72 :: file1

    real(kind=DP), allocatable :: pot_ref(:,:), pot_wk(:,:)

    allocate( pot_wk(ista_kngp:iend_kngp,kimg) );    pot_wk = 0.0d0
    call set_electrocstatic_pot_Gspace( pot_wk )

    allocate( pot_ref(ista_kngp:iend_kngp,kimg) );    pot_ref = 0.0d0

    call m_Files_open_electrosta_pot_ref()
    call rd_pot_bin( nf_elecpot_bin_ref, pot_ref )
    call m_Files_close_electrosta_pot_ref()

    pot_wk = pot_wk -pot_ref
    deallocate(pot_ref)

    nlp = fft_box_size_CD(1,1);
    nmp = fft_box_size_CD(2,1);
    nnp = fft_box_size_CD(3,1)
    allocate( potdiff_on_grid(nnp,nmp,nlp) );  potdiff_on_grid = 0.0d0
    call set_electrostatic_pot_on_grid( pot_wk, nlp, nmp, nnp, potdiff_on_grid )

    allocate( potdiff_on_atoms(natm) );  potdiff_on_atoms = 0.0d0
    call set_electrostatic_pot_on_atoms( pot_wk, potdiff_on_atoms )

    deallocate(pot_wk)

  end subroutine m_ESCS_calc_pot_diff

! ================================
  subroutine m_ESCS_finalize
    if ( allocated( potdiff_on_grid )  ) deallocate( potdiff_on_grid )
    if ( allocated( potdiff_on_atoms ) ) deallocate( potdiff_on_atoms )
    if ( allocated( Vpc_on_grid )  ) deallocate( Vpc_on_grid )
    if ( allocated( Vpc_on_atoms ) ) deallocate( Vpc_on_atoms )
  end subroutine m_ESCS_finalize

  subroutine m_ESCS_wd_extfnv_correction  ! auto
    integer :: nlp, nmp, nnp, lun, it
    character*72 :: file1
    real(kind=DP) :: Epc, zion_sum

    if ( atomid_is_set ) then
       pos_defect(1:3) = pos( atomid_defect,1:3 )
    endif

    zion_sum = 0.0d0
    Do it=1, ntyp
       zion_sum = zion_sum +ival(it) *iatom(it)
    End Do
    q = zion_sum -totch

    call m_ESCS_calc_pot_diff         ! pot_diff in Gspace

    nlp = fft_box_size_CD(1,1)
    nmp = fft_box_size_CD(2,1)
    nnp = fft_box_size_CD(3,1)
    allocate( Vpc_on_grid(nnp,nmp,nlp) ); Vpc_on_grid = 0.0d0
    allocate( Vpc_on_atoms(natm) ); Vpc_on_atoms = 0.0d0

    call calc_extended_fnv_correction( Epc )

    if ( mype == 0 ) then
       write(nfout,*)
       write(nfout,*) "=== Result: Charged Defect Analysis ==="
       write(nfout,'(A,F10.5)')    "    Ecorr     = Epc -q *dV"
       write(nfout,'(A,F10.5)')    "       q      = ", q
       write(nfout,'(A,F20.10,A)') "    Epc       = ", Epc *Hartree, " eV"
    endif

    lun = 1000;  file1 = "./defect_pot_correction.direction"
    if ( mype == 0 ) call calc_correction_cube( lun, file1 )

    lun = 1000;  file1 = "./defect_pot_correction.atoms"
    if ( mype == 0 ) call calc_correction_atoms( lun, file1 )

    write(nfout,*)

  contains

    subroutine calc_correction_cube( lun, file1 )
      integer, intent(in) :: lun
      character*72, intent(in) :: file1

      integer :: i, j, k
      integer :: nlp, nmp, nnp, ic(3)
      character*72 :: file_out
      real(kind=DP) :: vec_A(3), vec_B(3), vec_C(3), c1, c2, c3, r1, emax, emin
      real(kind=DP), allocatable :: pot_dir1(:,:), pot_dir2(:,:), pot_dir3(:,:)

      vec_A(1:3) = altv(1:3,1)
      vec_B(1:3) = altv(1:3,2)
      vec_C(1:3) = altv(1:3,3)
!
      nlp = fft_box_size_CD(1,1)
      nmp = fft_box_size_CD(2,1)
      nnp = fft_box_size_CD(3,1)
!
      allocate( pot_dir1(nlp,2) );  pot_dir1 = 0.0d0
      allocate( pot_dir2(nmp,2) );  pot_dir2 = 0.0d0
      allocate( pot_dir3(nnp,2) );  pot_dir3 = 0.0d0
!
      Do i=1, nlp
         Do j=1, nmp
            Do k=1, nnp
               c1 = -potdiff_on_grid(k,j,i) *Hartree    ! -1 is from sign of electron
               c2 = Vpc_on_grid(k,j,i) *Hartree

               pot_dir1(i,1) = pot_dir1(i,1) +c1
               pot_dir2(j,1) = pot_dir2(j,1) +c1
               pot_dir3(k,1) = pot_dir3(k,1) +c1
               pot_dir1(i,2) = pot_dir1(i,2) +c2
               pot_dir2(j,2) = pot_dir2(j,2) +c2
               pot_dir3(k,2) = pot_dir3(k,2) +c2
            End Do
         End Do
      End Do
!
      pot_dir1 = pot_dir1 /dble(nmp) /dble(nnp)
      pot_dir2 = pot_dir2 /dble(nnp) /dble(nlp)
      pot_dir3 = pot_dir3 /dble(nlp) /dble(nmp)

! direc vecA
      r1 = sqrt( Vec_A(1)**2 +Vec_A(2)**2 +Vec_A(3)**2 ) *Bohr
      file_out = trim(adjustl(file1)) // "_1"
      open( unit=lun, file=file_out, status="unknown", form="formatted" )
      write(lun,'(A,4X,A,5X,A,10X,A)') "#  dist. (Ang)", "pot_diff (eV)", &
           &         "Vpc (eV)", "pot_diff -Vpc (eV)"
      Do i=1, nlp
         c1 = r1 *(i-1) /dble(nlp)
         c2 = pot_dir1(i,1);    c3 = pot_dir1(i,2)
         write(lun,'(F10.5,3F18.8)') c1, c2, c3, c2-c3
      End Do
      close(lun)

! direc vecB
      r1 = sqrt( Vec_B(1)**2 +Vec_B(2)**2 +Vec_B(3)**2 ) *Bohr
      file_out = trim(adjustl(file1)) // "_2"
      open( unit=lun, file=file_out, status="unknown", form="formatted" )
      write(lun,'(A,4X,A,5X,A,10X,A)') "#  dist. (Ang)", "pot_diff (eV)", &
           &         "Vpc (eV)", "pot_diff -Vpc (eV)"
      Do i=1, nmp
         c1 = r1 *(i-1) /dble(nmp)
         c2 = pot_dir2(i,1);   c3 = pot_dir2(i,2)
         write(lun,'(F10.5,3F18.8)') c1, c2, c3, c2-c3
      End Do
      close(lun)

! direc vecC
      r1 = sqrt( Vec_C(1)**2 +Vec_C(2)**2 +Vec_C(3)**2 ) *Bohr
      file_out = trim(adjustl(file1)) // "_3"
      open( unit=lun, file=file_out, status="unknown", form="formatted" )
      write(lun,'(A,4X,A,5X,A,10X,A)') "#  dist. (Ang)", "pot_diff (eV)", &
           &         "Vpc (eV)", "pot_diff -Vpc (eV)"
      Do i=1, nnp
         c1 = r1 *(i-1) /dble(nnp)
         c2 = pot_dir3(i,1);     c3 = pot_dir3(i,2)
         write(lun,'(4F10.5,3F18.8)') c1, c2, c3, c2-c3
      End Do
      close(lun)
!
      if ( q > 0 ) then
         emin = 1.0D10;   ic(1) = 1
         Do i=1, nlp
            if ( pot_dir1(i,2) < emin ) then
               emin = pot_dir1(i,2);      ic(1) = i
            endif
         End Do
         emin = 1.0D10;   ic(2) = 1
         Do i=1, nmp
            if ( pot_dir2(i,2) < emin ) then
               emin = pot_dir2(i,2);      ic(2) = i
            endif
         End Do
         emin = 1.0D10;   ic(3) = 1
         Do i=1, nnp
            if ( pot_dir3(i,2) < emin ) then
               emin = pot_dir3(i,2);      ic(3) = i
            endif
         End Do
      else
         emax = -1.0D10;   ic(1) = 1
         Do i=1, nlp
            if ( pot_dir1(i,2) > emax ) then
               emax = pot_dir1(i,2);      ic(1) = i
            endif
         End Do
         emax = -1.0D10;   ic(2) = 1
         Do i=1, nmp
            if ( pot_dir2(i,2) > emax ) then
               emax = pot_dir2(i,2);      ic(2) = i
            endif
         End Do
         emax = -1.0D10;   ic(3) = 1
         Do i=1, nnp
            if ( pot_dir3(i,2) > emax ) then
               emax = pot_dir3(i,2);      ic(3) = i
            endif
         End Do
      endif

      c1 = ( pot_dir1(ic(1),1) -pot_dir1(ic(1),2) ) &
           & + ( pot_dir2(ic(2),1) -pot_dir2(ic(2),2) ) &
           & + ( pot_dir3(ic(3),1) -pot_dir3(ic(3),2) )
      c2 = c1 /3.0d0
      c3 = Epc *Hartree -q *c2

!      if ( iprichgdefect > 1 ) then
         write(nfout,'(A)') "  << if planar avarage is used >>"
         write(nfout,'(A,F20.10,A)') "    dV        = ", c2, " eV"
         write(nfout,'(A,F20.10,A)') "    Ecorr     = ", c3, " eV"
!      endif

      deallocate( pot_dir1 );   deallocate( pot_dir2 );  deallocate( pot_dir3 )

    end subroutine calc_correction_cube

    subroutine calc_correction_atoms( lun, file1 )
      integer, intent(in) :: lun
      character*72, intent(in) :: file1

      integer :: i, nx, ny, nz, num
      real(kind=DP) :: dist_min, dist, vec(3), rws
      real(kind=DP) :: x1, y1, z1, c1, c2, c3, csum
      real(kind=DP) :: vec_A(3), vec_B(3), vec_C(3)

      dist_min = 1.0D10
      Do i=1, 3
         vec(1:3) = altv(1:3,i) /2.0d0
         dist = sqrt( vec(1)**2 +vec(2)**2 +vec(3)**2 )
         dist_min = min( dist, dist_min )
      End Do
      rws = dist_min
!
      vec_A(1:3) = altv(1:3,1)
      vec_B(1:3) = altv(1:3,2)
      vec_C(1:3) = altv(1:3,3)
!
      num = 0 ; csum = 0.0d0

      open( unit=lun, file=file1, status="unknown", form="formatted" )
      write(lun,'(A,3X,A,4X,A,5X,A,10X,A)') "# no.", "dist. (Ang)", "pot_diff (eV)", &
           &         "Vpc (eV)", "pot_diff -Vpc (eV)"

      Do i=1, natm
         dist_min = 1.0D10
         Do nx=-1, 1
            Do ny=-1, 1
               Do nz=-1, 1
                  x1 = ( pos(i,1) -pos_defect(1) +nx )*vec_A(1) &
                       & +( pos(i,2) -pos_defect(2) +ny )*vec_B(1) &
                       & +( pos(i,3) -pos_defect(3) +nz )*vec_C(1)
                  y1 = ( pos(i,1) -pos_defect(1) +nx )*vec_A(2) &
                       & +( pos(i,2) -pos_defect(2) +ny )*vec_B(2) &
                       & +( pos(i,3) -pos_defect(3) +nz )*vec_C(2)
                  z1 = ( pos(i,1) -pos_defect(1) +nx )*vec_A(3) &
                       & +( pos(i,2) -pos_defect(2) +ny )*vec_B(3) &
                       & +( pos(i,3) -pos_defect(3) +nz )*vec_C(3)
                  dist = sqrt( x1**2 +y1**2 +z1**2 )
                  dist_min = min( dist, dist_min )
               End do
            End Do
         End Do

#if 0
         c1 = ( Epc +q*Vpc_on_atoms(i) ) *Hartree
         c2 =  q *potdiff_on_atoms(i) *Hartree      ! (-q) *(-1);
                                                    !    -1 is from sign of electron
         c3 = c1 +c2
#else
         c3 = ( -potdiff_on_atoms(i) -Vpc_on_atoms(i) ) *Hartree
#endif
         if ( dist_min > rws +0.1 ) then
            num = num +1;   csum = csum +c3
         endif
         write(lun,'(I5,F10.5,3F18.8)') i, dist_min*Bohr, &
              &                          -potdiff_on_atoms(i) *Hartree, &
              &                           Vpc_on_atoms(i) *Hartree, &
              &                         (-potdiff_on_atoms(i)-Vpc_on_atoms(i) )*Hartree
      End Do

!      write(nfout,*)
      write(nfout,'(A)') "  << if site avarage is used >>"
      write(nfout,'(A,F20.10,A)') "      ( Rad_ws = ", rws *Bohr, " Ang. )"

      if ( num > 0 ) then
         c2 = csum /dble(num)
         c3 = Epc *Hartree -q *c2
!         if ( iprichgdefect > 1 ) then
            write(nfout,'(A,F20.10,A)') "    dV        = ", c2, " eV"
!         endif
         write(nfout,'(A,F20.10,A)')    "    Ecorr     = ", c3, " eV"

         write(lun,*)
         write(lun,'(A,F20.10)') "# Correction energy (eV):", c3
      endif
      close(lun)

    end subroutine calc_correction_atoms

  end subroutine m_ESCS_wd_extfnv_correction

  subroutine calc_extended_fnv_correction( Epc_aniso )
    real(kind=DP), intent(out) :: Epc_aniso

    integer :: neibrd, newldg, alen(3)
    real(kind=DP) :: gamma, alf
    real(kind=DP) :: epsinv(3,3), det_eps
    real(kind=DP) :: Vpc_r, Vpc_g0
!
    real(kind=DP), parameter :: rsphere_radius = 12.5d0
    real(kind=DP), parameter :: phi            =  6.0d0
!
    real(kind=DP), allocatable :: rr(:), rxyz(:,:)
!
    integer :: i, nlp, nmp, nnp
    real(kind=DP) :: alpha_r, c1

    call calc_epsinv( eps, det_eps, epsinv )
    call decide_rxyz_size( rsphere_radius, alen, neibrd )
    allocate( rxyz(neibrd,3) );    allocate( rr(neibrd) )

    call substitute_rxyz( alen, neibrd, rxyz, rr )
    deallocate(rr)

    call decide_alf( alen, rsphere_radius, phi, alf, det_eps )   ! -> alf
    call decide_newldg( alf, phi, newldg, det_eps ) ! -> newldg = #Gvectors for summation
!
    gamma = 1.0d0 /alf
!
    call calc_Vpc_aniso_Rspace( q, gamma, neibrd, rxyz, det_eps, epsinv, Vpc_r )
    call calc_Vpc_aniso_Gspace0( q, gamma, newldg, eps, det_eps, Vpc_g0 )
!
    Epc_aniso = -q /2.0d0 *( Vpc_r +Vpc_g0 )

    nlp = fft_box_size_CD(1,1)
    nmp = fft_box_size_CD(2,1)
    nnp = fft_box_size_CD(3,1)
    call calc_Vpc_aniso_Gspace1_grid( q, gamma, newldg, eps, pos_defect, &
         &                            nlp, nmp, nnp, Vpc_on_grid )
    Vpc_on_grid = Vpc_r +Vpc_on_grid

    call calc_Vpc_aniso_Gspace1_atoms( q, gamma, newldg, eps, det_eps, pos_defect, &
         &                             Vpc_on_atoms )
    Vpc_on_atoms = Vpc_r +Vpc_on_atoms

    deallocate(rxyz);

  end subroutine calc_extended_fnv_correction

  subroutine decide_alf( alen, rsphere_radius, phi, alf, det_eps )
    integer, intent(in) :: alen(3)
    real(kind=DP), intent(in) :: rsphere_radius, phi, det_eps
    real(kind=DP), intent(out) ::alf

    real(kind=DP) :: xalen(3), aamin, c1

    c1 = det_eps **(1./3.)
    xalen = alen*(abs(int(rsphere_radius/alen)) + 1)
#if 0
    aamin = minval(xalen) /sqrt(c1)
#else
    aamin = minval(xalen)
#endif
    alf = aamin/phi
!    if(printable) write(nfout,'("  alf = ",f12.6," aamin = ",f12.6)') alf, aamin
  end subroutine decide_alf

  subroutine decide_newldg( alf, phi, newldg, det_eps )
    real(kind=DP), intent(in) :: alf, phi, det_eps
    integer, intent(out) :: newldg

    integer :: i
    integer :: iend, newldg_mpi  !mpi
    real(kind=DP) :: c1

    newldg_mpi= 1;    newldg= 1
!    if(printable) write(nfout,'(" ! kg = ",i9)') kg

    c1 = det_eps **(1./3.)

    iend = iend_kngp
    if( iend > kg ) iend = kg
    if( ista_kngp <= iend ) then
       do i = ista_kngp, iend  !for mpi
#if 0
         if ( alf*gr_l(i)*c1 < phi*2.d0 ) newldg_mpi = i
#else
         if ( alf*gr_l(i) < phi*2.d0 ) newldg_mpi = i
#endif
       end do
    endif
    if (npes > 1) then
       call mpi_allreduce(newldg_mpi,newldg,1,mpi_integer,mpi_max,MPI_CommGroup,ierr)
    else
       newldg = newldg_mpi
    end if
    if (newldg.eq.kg+1) then
!       if(printable) write(nfout,'(" **warn alf is too small: alf=",d20.10)') alf
       stop
    endif
!    if (printable) write(nfout,440) newldg
440 format(' ',' newldg = ',i8)
  end subroutine decide_newldg

  subroutine decide_rxyz_size( rsphere_radius, alen, neibrd )
    real(kind=DP), intent(in)  :: rsphere_radius
    integer,       intent(out) :: alen(3), neibrd

    integer :: i
    real(kind=DP) :: a

    do i = 1, 3
       a = dabs((rltv(1,i)*altv(1,i)+rltv(2,i)*altv(2,i)+rltv(3,i)*altv(3,i))&
            & /dsqrt(rltv(1,i)*rltv(1,i)+rltv(2,i)*rltv(2,i)+rltv(3,i)*rltv(3,i)))
       alen(i) = abs(int(rsphere_radius/a)) +1
    enddo
    neibrd = (alen(1)*2+1)*(alen(2)*2+1)*(alen(3)*2+1)

  end subroutine decide_rxyz_size

  subroutine substitute_rxyz( alen, neibrd, rxyz, rr )
    integer, intent(in) :: neibrd, alen(3)
    real(kind=DP), intent(out) :: rxyz(neibrd,3), rr(neibrd)

    integer :: i, j, k, mm
    real(kind=DP) :: f(3)

    mm = 0
    do i = -alen(1), alen(1)
       do j = -alen(2), alen(2)
          do k = -alen(3), alen(3)
             f(1) = i; f(2) = j; f(3) = k
             mm = mm + 1
             rxyz(mm,1:3) = matmul(altv,f)
             rr(mm) = dsqrt(dot_product(rxyz(mm,1:3),rxyz(mm,1:3)))
          enddo
       enddo
    enddo
#ifdef LIBRARY_BUILD
    call hpsort_phase0(neibrd,neibrd,rxyz,rr)
#else
    call hpsort(neibrd,neibrd,rxyz,rr)
#endif

  end subroutine substitute_rxyz

  subroutine calc_epsinv( eps, det, epsinv )
    real(kind=DP), intent(in) :: eps(3,3)
    real(kind=DP), intent(out) :: epsinv(3,3), det
!
    real(kind=DP) :: c11, c12, c13, c21, c22, c23, c31, c32, c33
    real(kind=DP) :: m11, m12, m13, m21, m22, m23, m31, m32, m33

    c11 = eps(1,1);   c12 = eps(1,2);   c13 = eps(1,3)
    c21 = eps(2,1);   c22 = eps(2,2);   c23 = eps(2,3)
    c31 = eps(3,1);   c32 = eps(3,2);   c33 = eps(3,3)

    det = c11*c22*c33 +c12*c23*c31 +c13*c21*c32 &
         &   -c13*c22*c31 -c12*c21*c33 -c11*c23*c32

    m11 = ( c22*c33 -c23*c32 ) /det
    m12 =-( c12*c33 -c13*c32 ) /det
    m13 = ( c12*c23 -c13*c22 ) /det
    m21 =-( c21*c33 -c23*c31 ) /det
    m22 = ( c11*c33 -c13*c31 ) /det
    m23 =-( c11*c23 -c13*c21 ) /det
    m31 = ( c21*c32 -c22*c31 ) /det
    m32 =-( c11*c32 -c12*c31 ) /det
    m33 = ( c11*c22 -c12*c21 ) /det
!
    epsinv(1,1) = m11;    epsinv(1,2) = m12;     epsinv(1,3) = m13
    epsinv(2,1) = m21;    epsinv(2,2) = m22;     epsinv(2,3) = m23
    epsinv(3,1) = m31;    epsinv(3,2) = m32;     epsinv(3,3) = m33

  end subroutine calc_epsinv

  subroutine calc_Vpc_aniso_Rspace( q, gamma, neibrd, rxyz, det_eps, epsinv, pot )
    integer, intent(in) :: neibrd
    real(kind=DP), intent(in) :: rxyz(neibrd,3)
    real(kind=DP), intent(in) :: q, gamma, det_eps, epsinv(3,3)
    real(kind=DP), intent(out) :: pot

    integer :: in, n1, n2
    real(kind=DP) :: csum1, csum2, r1, c1

    csum1 = 0.0d0
    do in = 2, neibrd
       r1 = 0.0d0
       Do n1=1, 3
          Do n2=1, 3
             r1 = r1 +rxyz(in,n1) *epsinv(n1,n2) *rxyz(in,n2)
          End do
       End Do
       c1 = erfc( gamma*r1 ) /r1
       csum1 = csum1 +c1
    End do
    csum1 = csum1 /sqrt( det_eps )
!
    csum2 = -PAI /univol /gamma**2
!
    pot = q *( csum1 +csum2 )

  end subroutine calc_Vpc_aniso_Rspace

  subroutine calc_Vpc_aniso_Gspace0( q, gamma, newldg, eps, det_eps, pot )
    real(kind=DP), intent(in) :: q, gamma, eps(3,3), det_eps
    real(kind=DP), intent(out) :: pot
    integer, intent(in) :: newldg
!
    integer :: ist, iend, in, k, n1, n2
    real(kind=DP) :: csum, csum_mpi, c1, gvec(3)
!
    ist = ista_kngp
    if (ist == 1) ist = 2

    iend = iend_kngp
    if ( iend > newldg ) iend = newldg

    csum = 0.0d0

    if( ist <= iend ) then
       do in = ist, iend  !for mpi
          Do k=1, 3
             gvec(k) = rltv(k,1) *ngabc(in,1) +rltv(k,2)*ngabc(in,2) &
                  &  + rltv(k,3) *ngabc(in,3)
          End Do
          c1 = 0.0d0
          Do n1=1, 3
             Do n2=1, 3
                c1 = c1 +gvec(n1) *eps(n1,n2) *gvec(n2)
             End Do
          End do
          csum = csum +exp( -c1 /gamma**2 /4.0d0 ) /c1
       End do
    endif

    if (npes > 1) then
       call mpi_allreduce( csum, csum_mpi, 1, mpi_double_precision, mpi_sum, &
            &              MPI_CommGroup, ierr )
       csum = csum_mpi
    end if
    pot = PAI4 *q /univol *csum -2.0d0 *gamma *q /sqrt( PAI *det_eps )

  end subroutine calc_Vpc_aniso_Gspace0

  subroutine calc_Vpc_aniso_Gspace1_atoms( q, gamma, newldg, eps, det_eps, pos_defect, &
       &                                   pot_on_atom )
    real(kind=DP), intent(in) :: q, gamma, eps(3,3), pos_defect(3), det_eps
    real(kind=DP), intent(inout) :: pot_on_atom(natm)
    integer, intent(in) :: newldg
!
    integer :: ist, iend, in, k, n1, n2, ia
    real(kind=DP) :: csum, c1, phs, gvec(3), vec_r(3)
    real(kind=DP), allocatable :: zfcos(:), zfsin(:)
    real(kind=DP), allocatable :: weight(:), work(:)
!
    ist = ista_kngp
    if (ist == 1) ist = 2

    iend = iend_kngp
    if ( iend > newldg ) iend = newldg

    allocate(zfcos(ista_kngp:iend_kngp)); zfcos = 0.d0
    allocate(zfsin(ista_kngp:iend_kngp)); zfsin = 0.d0

    allocate( weight(ista_kngp:iend_kngp) );  weight = 0.0d0

    if ( ist <= iend ) then
       do in = ist, iend  !for mpi
          Do k=1, 3
             gvec(k) = rltv(k,1) *ngabc(in,1) +rltv(k,2)*ngabc(in,2) &
                  &  + rltv(k,3) *ngabc(in,3)
          End Do
          c1 = 0.0d0
          Do n1=1, 3
             Do n2=1, 3
                c1 = c1 +gvec(n1) *eps(n1,n2) *gvec(n2)
             End Do
          End do
          weight(in) = exp( -c1 /gamma**2 /4.0d0 ) /c1
       end do
    end if

    Do ia=1, natm
       vec_r(:) = pos(ia,:) -pos_defect(:)
       call calc_phase2( 1, vec_r, 1, kgp, ngabc, ista_kngp, iend_kngp, &
            &            zfcos, zfsin )
       csum = 0.0d0
       do in=ist, iend
          csum = csum +weight(in) *zfcos(in)
       end do
#if 0
       pot_on_atom(ia) = q *PAI4 /univol *csum &
            &         -2.0d0 *gamma *q /sqrt( PAI *det_eps )
#else
       pot_on_atom(ia) = q *PAI4 /univol *csum
#endif
    End Do
    deallocate( zfcos, zfsin )

    if ( npes > 1 ) then
       allocate( work(natm) ); work = 0.0d0
       call mpi_allreduce( pot_on_atom, work, natm, mpi_double_precision, mpi_sum, &
            &              MPI_CommGroup, ierr )
       pot_on_atom = work
       deallocate( work )
    endif
    deallocate( weight )

  end subroutine calc_Vpc_aniso_Gspace1_atoms

  subroutine calc_Vpc_aniso_Gspace1_grid( q, gamma, newldg, eps, pos_defect, &
       &                                  nlp, nmp, nnp, pot_on_grid )
    real(kind=DP), intent(in) :: q, gamma, eps(3,3), pos_defect(3)
    real(kind=DP), intent(inout) :: pot_on_grid(nnp,nmp,nlp)
    integer, intent(in) :: newldg, nlp, nmp, nnp
!
    integer :: i, j, k
    integer :: ist, iend, in, n1, n2
    real(kind=DP) :: csum, c1, c2, phs, gvec(3)
    real(kind=DP), allocatable :: zfcos(:), zfsin(:)
    real(kind=DP), allocatable :: pot_wk(:,:), work(:,:,:)
!
    ist = ista_kngp
    if (ist == 1) ist = 2

    iend = iend_kngp
    if ( iend > newldg ) iend = newldg

    allocate(zfcos(ista_kngp:iend_kngp)); zfcos = 0.d0
    allocate(zfsin(ista_kngp:iend_kngp)); zfsin = 0.d0
    allocate( pot_wk(ista_kngp:iend_kngp,kimg) ); pot_wk = 0.0d0

    call calc_phase2( 1, pos_defect, 1, kgp, ngabc, ista_kngp, iend_kngp, &
         &            zfcos, zfsin )

    if ( ist <= iend ) then
       do in = ist, iend  !for mpi
          Do k=1, 3
             gvec(k) = rltv(k,1) *ngabc(in,1) +rltv(k,2)*ngabc(in,2) &
                  &  + rltv(k,3) *ngabc(in,3)
          End Do
          c1 = 0.0d0
          Do n1=1, 3
             Do n2=1, 3
                c1 = c1 +gvec(n1) *eps(n1,n2) *gvec(n2)
             End Do
          End do
          c2 = exp( -c1 /gamma**2 /4.0d0 ) /c1
          pot_wk(in,1) =  c2 *zfcos(in)
          pot_wk(in,2) = -c2 *zfsin(in)
       end do
    end if
!
    call set_electrostatic_pot_on_grid( pot_wk, nlp, nmp, nnp, pot_on_grid )
    pot_on_grid = pot_on_grid *q *PAI4 /univol

    deallocate( zfcos, zfsin )
    deallocate( pot_wk )

  end subroutine calc_Vpc_aniso_Gspace1_grid

end module m_ES_ChargeState
