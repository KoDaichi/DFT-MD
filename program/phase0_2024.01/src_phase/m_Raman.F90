!
! =========== Contributions ===================================
!
! Through the courtesy of contributors, this module is added.
!
! Company:  ASMS Co.,Ltd.
! Functions:  [Identifier: 13.1R]
!       Raman tensor/spectrum calculation
!
! Note:
!       In order to work the whole functions, LIB_ASMS_RAMAN is required.
!
! =============================================================

module m_Raman
  use m_Control_Parameters, only : sw_lo_to_splitting, sw_phonon_oneshot, &
       &                           printable, ipriphonon, &
       &                           sw_phonon_with_epsilon, sw_calc_dielectric_tensor, &
       &                           num_phonon_calc_mode, sw_excitation, &
       &                           sw_phonon, sw_calc_force, sw_use_add_proj, sw_raman

  use m_Const_Parameters, only : DP, ON, OFF, FMAXVALLEN, LOWER, PAI4, CMPLDP

  use m_Crystal_Structure, only : univol, op
  use m_Files,  only : nfout, nfeps_ph, nfoptical_coeff, nframan_spectra, &
       &               m_Files_open_nfeps_phonon, &
       &               m_Files_open_nfoptical_coeff, &
       &               m_Files_open_nframan_spectra

  use m_Ionic_System,       only : num_force_data, &
       &                           displaced_atom, displacement, u, &
       &                           natm,ionic_mass, speciesname, &
       &                           istart_phonon, iend_phonon, &
       &                           natm_prim, norder, sw_polynomial_fit, &
       &                           num_force_calc, iconf, &
       &                           napt_phonon, iequconf, &
       &                           iopr_equconf, force_was_read, phonon_iteration, &
       &                           natm_mobile, phonon_atom, imdtyp

  use m_Parallelization,    only : mype,npes,ierr, MPI_CommGroup
  use m_BP_Properties,      only : zeff

  use m_IterationNumbers, only : iteration_ionic, iteration_electronic

  implicit none
  include 'mpif.h'

! ----
! method
!
  character(len("raman") ), private, parameter :: &
       &                    tag_raman = "raman"
  character(len("scheme") ), private, parameter :: &
       &                     tag_raman_calc_scheme = "scheme"
  character(len("classical") ), private, parameter :: &
       &                     tag_classical = "classical"
!
  character(len("sw_raman") ), private, parameter :: &
       &                    tag_sw_raman = "sw_raman"
  character(len("sw_phonon_with_epsilon") ), private, parameter :: &
       &                    tag_sw_phonon_with_epsilon = "sw_phonon_with_epsilon"
  character(len("sw_calc_dielectric_tensor") ), private, parameter :: &
       &                    tag_sw_calc_dielectric_tensor = "sw_calc_dielectric_tensor"
  character(len("sw_froehlich_correction") ), private, parameter :: &
       &                    tag_sw_froehlich_correction = "sw_froehlich_correction"
  character(len("optical_coeff_from_file") ), private, parameter :: &
       &                    tag_optical_coeff_from_file = "optical_coeff_from_file"
!
!!!  integer :: sw_phonon_with_epsilon = OFF
!!  integer :: sw_calc_dielectric_tensor = OFF
  integer :: sw_froehlich_correction = OFF
  integer :: optical_coeff_from_file = OFF
!
  real(kind=DP), allocatable :: dielectric_tensors(:,:,:)
  real(kind=DP), allocatable :: dielectric_tensor0(:)
  real(kind=DP), allocatable :: dielectric_tensor_data(:,:)
  real(kind=DP), allocatable :: Mat_dEps_dR(:,:,:)
  real(kind=DP), allocatable :: Raman_susceptibility(:,:,:)
  real(kind=DP) :: nonl_op_coeff_rchi2(3,3,3)
!
! ----
! crystal type
!
  integer, parameter :: SingleCrystal = 1, PolyCrystal = 2
  integer :: Crystal_Type = SingleCrystal
!
  character(len("crystal_type") ), private, parameter :: &
       &                    tag_crystal_type = "crystal_type"
  character(len("single") ), private, parameter :: tag_single_crystal = "single"
  character(len("poly") ), private, parameter ::   tag_poly_crystal = "poly"
  character(len("powder") ), private, parameter ::   tag_powder = "powder"
!
! ----
! propagation of photon
!
  character(len("tag_x") ), private, parameter ::   tag_x = "x"
  character(len("tag_y") ), private, parameter ::   tag_y = "y"
  character(len("tag_z") ), private, parameter ::   tag_z = "z"

  character(len("propagation") ), private, parameter ::  &
       &                                 tag_propagation = "propagation"
  character(len("incoming") ), private, parameter ::  &
       &                                 tag_incoming = "incoming"
  character(len("incident") ), private, parameter ::  &
       &                                 tag_incident = "incident"
  character(len("outgoing") ), private, parameter ::  &
       &                                 tag_outgoing = "outgoing"
  character(len("scattered") ), private, parameter ::  &
       &                                 tag_scattered = "scattered"

  character(len("direction_in") ), private, parameter ::  &
       &                                 tag_direction_in = "direction_in"
  character(len("direction_out") ), private, parameter :: &
       &                                 tag_direction_out="direction_out"
  character(len("kx") ), private, parameter :: tag_kx = "kx"
  character(len("ky") ), private, parameter :: tag_ky = "ky"
  character(len("kz") ), private, parameter :: tag_kz = "kz"
!
  real(kind=DP) :: photon_k_in(3)  = (/ 0.0d0, 0.0d0, 1.0d0 /)
  real(kind=DP) :: photon_k_out(3) = (/ 0.0d0, 0.0d0,-1.0d0 /)

! ----
! laser
!
  character(len=FMAXVALLEN) :: rstr
  character(len("tag_laser") ), private, parameter ::   tag_laser = "laser"
  character(len("tag_photon") ), private, parameter ::  tag_photon = "photon"

  character(len("wavelength") ), private, parameter :: tag_wavelength = "wavelength"
!
  real(kind=DP), parameter :: wavelength_green1 = 514.5d0         ! nm
  real(kind=DP) :: wavelength = wavelength_green1
  real(kind=DP) :: omega_L

! ----
! polarization vectors
!
  character(len("polarization") ), private, parameter :: &
       &                  tag_polarization = "polarization"
!
  character(len("nonorthogonal_component") ), private, parameter :: &
       &                  tag_nonorthogonal_component = "nonorthogonal_component"
!
  character(len("rotation") ), private, parameter ::&
       &                   tag_rotation = "rotation"
  character(len("angle_dependence") ), private, parameter ::&
       &                   tag_angle_dependence = "angle_dependence"
  character(len("num_rot_angles") ), private, parameter ::&
       &                   tag_num_rot_angles = "num_rot_angles"
  character(len("angle_0") ), private, parameter ::&
       &                     tag_angle0 = "angle_0"
  character(len("angle_1") ), private, parameter ::&
       &                     tag_angle1 = "angle_1"
  character(len("parallel_or_normal") ), private, parameter ::&
       &                   tag_parallel_or_normal = "parallel_or_normal"
!
  integer :: polar_parallel_or_normal = 0
  integer :: num_rot_angles = 0
  integer, parameter :: max_num_rot_angles = 10
!
  real(kind=DP) :: nonorthogonal_component = 0.0d0
  real(kind=DP) :: rot_angle0 = 0.0d0, rot_angle1 = 0.0d0
  real(kind=DP) :: rot_angles( max_num_rot_angles )
! --
!  unpolarized config (single crystal)
!
  character(len("calc_unpolarized_config") ), private, parameter :: &
       &                   tag_calc_unpolarized_config = "calc_unpolarized_config"

  logical :: calc_unpolarized_config = .false.

! ----
! spectrum
!
  character(len("spectrum") ), private, parameter :: tag_spectrum = "spectrum"
  character(len("linewidth") ), private, parameter :: tag_linewidth = "linewidth"
  character(len("hwhm") ), private, parameter :: tag_hwhm = "hwhm"
!
  character(len("freq_pitch") ), private, parameter :: tag_freq_pitch = "freq_pitch"
  character(len("resolution") ), private, parameter :: tag_freq_resol = "resolution"
!
  real(kind=DP) :: raman_spectra_hwhm = 5.0d0                ! in cm-1
  real(kind=DP) :: raman_spectra_freq_pitch = 1.0d0          ! in cm-1

! ---
! temperature
!
  character(len("temperature") ), private, parameter :: tag_temperature = "temperature"
  character(len("temp") ), private, parameter :: tag_temp = "temp"
!
  real(kind=DP) :: raman_spectra_temperature = 0.0d0         ! in Kelvin


! ----
! general
!
  integer, parameter :: nqvec = 1

  integer :: nmodes
  real(kind=DP) :: qdir(nqvec,3)
  real(kind=DP) :: dielectric(3,3)

contains

  subroutine m_Raman_read_nfinp
    integer :: iret, f_selectBlock, f_getIntValue, f_getRealValue, f_getStringValue
    integer :: f_selectParentBlock, f_selectTop
    integer :: i
    real(kind=DP) :: dret, c1
    character(len=FMAXVALLEN) :: cret

    logical :: tf

    if( f_selectBlock( tag_raman ) == 0) then
! === 2015/10/19
       if ( f_getIntValue( tag_sw_raman, iret) == 0 ) sw_raman = iret
       if ( sw_raman == OFF ) return

       if ( sw_phonon == ON ) sw_phonon_with_epsilon = on
! === 2015/10/19

       if ( f_getIntValue( tag_sw_phonon_with_epsilon, iret) == 0) &
            &  sw_phonon_with_epsilon = iret

       if ( sw_phonon_with_epsilon == ON ) then
          sw_excitation = ON
          write(nfout,*) "!** sw_excitation is turned on"
! === 2015/10/19
          sw_use_add_proj = ON
          write(nfout,*) "!** sw_use_add_proj is turned on"
! === 2015/10/19
       endif

       if ( f_getStringValue( tag_raman_calc_scheme, rstr, LOWER) == 0 ) then
          call strncmp0( tag_classical, trim(rstr), tf )
          if(tf) sw_phonon_with_epsilon = on
       endif

! === 2015/10/19
       if ( sw_phonon_with_epsilon == ON ) then
          if ( sw_calc_force == ON ) sw_calc_dielectric_tensor = ON
       endif
! === 2015/10/19

       if ( f_getIntValue( tag_sw_calc_dielectric_tensor, iret) == 0) &
            &  sw_calc_dielectric_tensor = iret

       if ( f_getIntValue( tag_optical_coeff_from_file, iret) == 0) &
            &  optical_coeff_from_file = iret

       if ( f_getIntValue( tag_sw_froehlich_correction, iret) == 0) &
            &  sw_froehlich_correction = iret
!
       if( f_getStringValue( tag_crystal_type, rstr, LOWER) == 0 ) then
          call set_crystal_type( rstr, crystal_type )
       endif

       if( f_getRealValue( tag_temperature, dret,'') == 0 &
            &  .or. f_getRealValue( tag_temp, dret,'') == 0 ) then
          raman_spectra_temperature = dret            ! in Kelvin
          if ( raman_spectra_temperature < 0.0 ) then
             raman_spectra_temperature = 0.0d0
          endif
       endif

       if ( f_selectBlock( tag_propagation ) == 0 ) then
          if( f_getIntValue( tag_direction_in, iret ) == 0 ) then
             call set_direction_of_light( iret, photon_k_in )
          endif
          if( f_getIntValue( tag_direction_out, iret ) == 0 ) then
             call set_direction_of_light( iret, photon_k_out )
          endif
          if ( f_selectBlock( tag_incoming ) == 0 &
               & .or. f_selectBlock( tag_incident ) == 0 ) then
             if( f_getRealValue( tag_kx, dret, '' ) == 0 ) photon_k_in(1) = dret
             if( f_getRealValue( tag_ky, dret, '' ) == 0 ) photon_k_in(2) = dret
             if( f_getRealValue( tag_kz, dret, '' ) == 0 ) photon_k_in(3) = dret
             call set_direction_of_light2( photon_k_in )
             iret = f_selectParentBlock()
          end if
          if ( f_selectBlock( tag_outgoing ) == 0 &
               & .or. f_selectBlock( tag_scattered ) == 0 ) then
             if( f_getRealValue( tag_kx, dret, '' ) == 0 ) photon_k_out(1) = dret
             if( f_getRealValue( tag_ky, dret, '' ) == 0 ) photon_k_out(2) = dret
             if( f_getRealValue( tag_kz, dret, '' ) == 0 ) photon_k_out(3) = dret
             call set_direction_of_light2( photon_k_out )
             iret = f_selectParentBlock()
          endif
          iret = f_selectParentBlock()
       endif

       if ( f_selectBlock( tag_polarization ) == 0 ) then

          if ( f_getRealValue( tag_nonorthogonal_component, dret, '' ) == 0 ) then
             if ( dret > -1 .and. dret < 1 ) then
                nonorthogonal_component = dret
             endif
          endif

          if ( f_selectBlock( tag_angle_dependence ) == 0 &
               &       .or. f_selectBlock( tag_rotation ) == 0 ) then

             if ( f_getIntValue( tag_num_rot_angles, iret ) == 0 ) then
                if ( iret > 0 .and. iret <= max_num_rot_angles ) num_rot_angles = iret
             endif

             if ( num_rot_angles > 1 ) then
                if( f_getRealValue( tag_angle0, dret,'') == 0 ) rot_angle0 = dret
                if( f_getRealValue( tag_angle1, dret,'') == 0 ) rot_angle1 = dret

                c1 = ( rot_angle1 -rot_angle0 ) /dble( num_rot_angles -1 )
                Do i=1, num_rot_angles
                   rot_angles(i) = rot_angle0 +c1*(i-1)
                End do

             else if ( num_rot_angles == 1 ) then
                if ( f_getRealValue( tag_angle0, dret, '' ) == 0 )  rot_angle0 = dret
                rot_angles(1) = rot_angle0
             endif

             if ( f_getIntValue( tag_parallel_or_normal, iret ) == 0 ) then
                if ( iret == 0 .or. iret == 1 ) then
                   polar_parallel_or_normal = iret
                endif
             endif
             iret = f_selectParentBlock()

          endif

          if ( f_getIntValue( tag_calc_unpolarized_config, iret ) == 0 ) then
             if ( iret == 1 ) calc_unpolarized_config = .true.
          endif
          iret = f_selectParentBlock()
       end if
! ---
       if ( f_selectBlock( tag_laser ) == 0 ) then
          if( f_getRealValue( tag_wavelength, dret,'') == 0 ) then
             wavelength = dret         ! in nm
          endif

          if ( wavelength < 0.0 )  wavelength = wavelength_green1
! --
          iret = f_selectParentBlock()
       endif

       if( f_selectBlock( tag_spectrum ) == 0) then
          if( f_getRealValue( tag_hwhm, dret,'') == 0 &
               &  .or. f_getRealValue( tag_linewidth, dret,'') == 0 ) then
             raman_spectra_hwhm = dret                                   ! in cm-1
             if ( raman_spectra_hwhm < 0.0 ) raman_spectra_hwhm = 5.0d0
          endif

          if( f_getRealValue( tag_freq_pitch, dret,'') == 0 &
               &  .or. f_getRealValue( tag_freq_resol, dret,'') == 0 ) then

             raman_spectra_freq_pitch = dret            ! in cm-1
             if ( raman_spectra_freq_pitch < 0.0 ) then
                raman_spectra_freq_pitch = raman_spectra_hwhm *0.1d0
             endif
          else
             raman_spectra_freq_pitch = raman_spectra_hwhm *0.1d0
          endif
          iret = f_selectParentBlock()
       endif

       omega_L = 45.56335D0 /wavelength       ! in hartree

       iret = f_selectParentBlock()
    end if

  end subroutine m_Raman_read_nfinp

  subroutine m_Raman_initialize( kvec_dir )
    real(kind=DP), intent(out) :: kvec_dir(3)

    real(kind=DP) :: dnorm

    if ( sw_phonon_with_epsilon == OFF ) return

    if ( crystal_type == SingleCrystal ) then
       kvec_dir = photon_k_out -photon_k_in
       dnorm = sqrt(sum(kvec_dir(1:3)**2))

       if(abs(dnorm) < 1.d-10) then
          kvec_dir = -photon_k_in;      dnorm = 1.d0
       end if
       kvec_dir(1:3) = kvec_dir(1:3)/dnorm
    endif

    qdir(1,:) = kvec_dir(:)
! ----
    if ( sw_froehlich_correction == ON ) optical_coeff_from_file = ON
    if ( optical_coeff_from_file == ON ) call m_Raman_read_optical_coeff

    if ( sw_froehlich_correction == ON )  sw_lo_to_splitting = ON     ! ???????
!
  end subroutine m_Raman_initialize

  subroutine m_Raman_print_param
    if ( sw_raman == OFF ) return
    if ( sw_phonon_with_epsilon == OFF ) return

    write(nfout,'(A)') '!** ----- RAMAN setup ----- '
    write(nfout,*) '!** sw_phonon_with_epsilon = ', sw_phonon_with_epsilon
    write(nfout,*) '!** sw_calc_dielectric_tensor = ', sw_calc_dielectric_tensor
    write(nfout,*) '!** sw_froehlich_correction = ', sw_froehlich_correction

    write(nfout,*) '!** crystal_type = ', crystal_type
    write(nfout,'(A,F15.8)') '!** temperature [K] = ', &
         &                     raman_spectra_temperature

    write(nfout,'(A)') '!** propagation :'
    write(nfout,'(A,3F15.8)') '!** photon_k_in  = ', photon_k_in
    write(nfout,'(A,3F15.8)') '!** photon_k_out = ', photon_k_out
    write(nfout,'(A,3F15.8)') '!** qdir = ', qdir(1,:)

    if ( crystal_type == SingleCrystal ) then
       write(nfout,'(A)') '!** polarization :'
       write(nfout,'(A,F15.8)') '!** non-orthogonal component = ', &
            &                   nonorthogonal_component
       write(nfout,'(A,I0)') '!** polar_parallel_or_normal = ', polar_parallel_or_normal
       write(nfout,'(A,I0)') '!** num_rot_angles = ', num_rot_angles
       if ( num_rot_angles > 0 ) then
          write(nfout,'(A,2F15.8)') '!** First, Last = ', rot_angle0, rot_angle1
       endif
       if ( num_rot_angles == 0 ) then
          write(nfout,*) '!** calc_unpolarized_config = ', calc_unpolarized_config
       endif
    endif

    write(nfout,'(A)') '!** laser :'
    write(nfout,'(A,F15.8)') '!** wavelength [nm] = ', wavelength

    write(nfout,'(A)') '!** spectrum :'
    write(nfout,'(A,F15.8)') '!** Half width half maximum of spectr [cm-1] = ', &
         &                     raman_spectra_hwhm
    write(nfout,'(A,F15.8)') '!** Freq. pitch of spectra [cm-1] = ', &
         &                     raman_spectra_freq_pitch
    write(nfout,'(A)') '!** -----------------------'

  end subroutine m_Raman_print_param

  subroutine set_crystal_type( rstr, crystal_type )
    character(len=FMAXVALLEN),intent(in) :: rstr
    integer, intent(out) :: crystal_type
    logical :: tf

    call strncmp0( tag_single_crystal, trim(rstr), tf )
    if(tf) crystal_type = SingleCrystal

    call strncmp0( tag_poly_crystal, trim(rstr), tf )
    if(tf) crystal_type = PolyCrystal

    call strncmp0( tag_powder, trim(rstr), tf )
    if(tf) crystal_type = PolyCrystal

  end subroutine set_crystal_type

  subroutine set_direction_of_light( nn, direction )
    integer :: nn
    real(kind=DP), intent(out) :: direction(3)

    direction = 0.0d0

    if ( nn ==  1 ) direction(1) =  1.0d0
    if ( nn ==  2 ) direction(2) =  1.0d0
    if ( nn ==  3 ) direction(3) =  1.0d0
    if ( nn == -1 ) direction(1) = -1.0d0
    if ( nn == -2 ) direction(2) = -1.0d0
    if ( nn == -3 ) direction(3) = -1.0d0

    if ( direction(1)**2 +direction(2)**2 + direction(3)**2 < 1.0E-4 ) then
       direction(3) = 1.0d0
    endif
  end subroutine set_direction_of_light

  subroutine set_direction_of_light2( direction )
    real(kind=DP), intent(inout) :: direction(3)
    real(kind=DP) :: c1

    c1 = direction(1)**2 +direction(2)**2 + direction(3)**2
    if ( c1 < 1.0E-4 ) then
       direction = 0.0d0;  direction(3) = 1.0d0
    else
       direction = direction /sqrt(c1)
    endif
  end subroutine set_direction_of_light2

  subroutine m_Raman_write_dielec_tensors
    integer :: i,ic
    logical :: ini

    if (.not.allocated(dielectric_tensor_data)) then
       allocate(dielectric_tensor_data(6,num_force_data))
       ic = iconf(iteration_ionic)
       dielectric_tensor_data(1,ic) = dielectric(1,1)
       dielectric_tensor_data(2,ic) = dielectric(2,2)
       dielectric_tensor_data(3,ic) = dielectric(3,3)
       dielectric_tensor_data(4,ic) = dielectric(1,2)
       dielectric_tensor_data(5,ic) = dielectric(1,3)
       dielectric_tensor_data(6,ic) = dielectric(2,3)
    endif

    if ( mype /= 0 ) return

    ini = .false.
    do i=1,3*natm_prim*norder*2
       if ( force_was_read(i) ) then
         ini = .true.
         exit
       endif
    enddo

!    if(iteration_ionic == 1) then
    if (phonon_iteration == 1 .and. .not. ini) then
       call m_Files_open_nfeps_phonon(.false.)     ! new file (asms)
       if(istart_phonon == 1) then
          if(sw_phonon_oneshot == ON) then
             write(nfeps_ph,'(4(1x,i7),3A)') num_force_calc, norder, sw_polynomial_fit, &
                  &                          num_phonon_calc_mode, &
                  &                          "  : num_force_calc, norder, ", &
                  &                          " sw_polynomial_fit, ", &
                  &                          " num_phonon_calc_mode"
          else
             write(nfeps_ph,'(3(1x,i7),A)') num_force_calc, norder, sw_polynomial_fit, &
                  &  "  : num_force_calc, norder, sw_polynomial_fit"
          end if
       end if
    else
       call m_Files_open_nfeps_phonon(.true.)     ! append (asms)
    endif

    write(nfeps_ph,'(i4,3(1x,e25.12),2i6,2A)') displaced_atom, displacement(1:3),&
         &                                    iteration_ionic, iteration_electronic, &
         &                                    " : displaced_atom, displacement(1:3), ", &
         &                                    "iteration_ionic,iteration_electronic"

#if 0
    write(nfeps_ph,'(3F20.15)') dielectric(1,1), dielectric(2,2), &
         &                      dielectric(3,3), dielectric(1,2), &
         &                      dielectric(1,3), dielectric(2,3)
#else
    write(nfeps_ph,'(5X,3F22.15)') dielectric(1,1), dielectric(1,2), dielectric(1,3)
    write(nfeps_ph,'(5X,3F22.15)') dielectric(2,1), dielectric(2,2), dielectric(2,3)
    write(nfeps_ph,'(5X,3F22.15)') dielectric(3,1), dielectric(3,2), dielectric(3,3)
#endif
    call flush(nfeps_ph)

  end subroutine m_Raman_write_dielec_tensors

  subroutine m_Raman_read_dielec_tensors()
#ifdef WINDOWS
    use ifport
#endif
    implicit none

    ! local variables
    integer :: n,i,ia,ja,ind,j,ic,iopr
    integer :: istart,idummy, iter, jcount
    character(len=5) :: num

    integer :: n1, n2, m1, m2
    real(kind=DP) :: cmat1(3,3), cmat2(3,3), ctmp

    logical :: exi

    if(ipriphonon>=1) write(nfout,*) '<< m_Raman_read_dielectric_tensors >> Reading '
    call flush(nfout)

!!!    nmodes = natm_prim*3
    nmodes = natm_mobile*3

    allocate(dielectric_tensors(6,nmodes,norder*2))
    allocate(dielectric_tensor0(6))

!!!    num_force_data = nmodes*norder*2
    num_force_data = natm_prim*3 *norder*2

    if(sw_polynomial_fit == ON) num_force_data = num_force_data + 1
    if(.not.allocated(dielectric_tensor_data)) then
       allocate(dielectric_tensor_data(6,num_force_data))
    endif

    if(mype == 0) then
       call m_Files_open_nfeps_phonon(.false.)
       inquire( unit=nfeps_ph, exist=exi )
       if ( exi ) then
         rewind nfeps_ph
         read(nfeps_ph,*,end=100,err=100) num_force_calc, norder, sw_polynomial_fit
         do while(.true.)
            read(nfeps_ph,*,end=100,err=100) displaced_atom,displacement(1:3),iter
            ic = iconf(iter)
            read(nfeps_ph,*,end=100,err=100) cmat1(1,1), cmat1(1,2), cmat1(1,3)
            read(nfeps_ph,*,end=100,err=100) cmat1(2,1), cmat1(2,2), cmat1(2,3)
            read(nfeps_ph,*,end=100,err=100) cmat1(3,1), cmat1(3,2), cmat1(3,3)
            dielectric_tensor_data(1,ic) = cmat1(1,1)
            dielectric_tensor_data(2,ic) = cmat1(2,2)
            dielectric_tensor_data(3,ic) = cmat1(3,3)
            dielectric_tensor_data(4,ic) = cmat1(1,2)
            dielectric_tensor_data(5,ic) = cmat1(1,3)
            dielectric_tensor_data(6,ic) = cmat1(2,3)
!            dielectic_tensor_was_read(ic) = .true.
         end do
100      continue
       endif

    end if

    if(ipriphonon>=2) write(nfout,'(" num_force_data = ",i8)') num_force_data
    do i=1,num_force_data
       if(imdtyp(phonon_atom(i))==OFF) cycle

       ic = iequconf(i)
       if(ipriphonon>=2) write(nfout,'(" i, ic = ",2i8)') i, ic

       if (ic == 0) cycle

       ic = abs(ic)
       iopr = iopr_equconf(i)
       if(ipriphonon>=2) write(nfout,'(" i, ic, iopr = ",3i8)') i, ic,iopr

       cmat1(1,1) = dielectric_tensor_data(1,ic)
       cmat1(2,2) = dielectric_tensor_data(2,ic)
       cmat1(3,3) = dielectric_tensor_data(3,ic)
       cmat1(1,2) = dielectric_tensor_data(4,ic)
       cmat1(1,3) = dielectric_tensor_data(5,ic)
       cmat1(2,3) = dielectric_tensor_data(6,ic)

       cmat1(2,1) = cmat1(1,2)
       cmat1(3,1) = cmat1(1,3)
       cmat1(3,2) = cmat1(2,3)

       Do n1=1, 3
          Do n2=1, 3
             ctmp = 0.0d0
             Do m1=1, 3
                Do m2=1, 3
                   ctmp = ctmp + op(n1,m1,iopr) *cmat1(m1,m2) *op(n2,m2,iopr)
                End do
             End do
             cmat2(n1,n2) = ctmp
          End do
       End Do

       dielectric_tensor_data(1,i) = cmat2(1,1)
       dielectric_tensor_data(2,i) = cmat2(2,2)
       dielectric_tensor_data(3,i) = cmat2(3,3)
       dielectric_tensor_data(4,i) = cmat2(1,2)
       dielectric_tensor_data(5,i) = cmat2(1,3)
       dielectric_tensor_data(6,i) = cmat2(2,3)
    end do

    if ( ipriphonon >=2 ) then
       if ( mype == 0) then
          write(nfout,*) '** dielectroc_tensor_data after rotation'
          Do ic=1, num_force_data
             write(nfout,*) '** ic = ',ic
             write(nfout,*) dielectric_tensor_data(1,ic)
             write(nfout,*) dielectric_tensor_data(2,ic)
             write(nfout,*) dielectric_tensor_data(3,ic)
             write(nfout,*) dielectric_tensor_data(4,ic)
             write(nfout,*) dielectric_tensor_data(5,ic)
             write(nfout,*) dielectric_tensor_data(6,ic)
          End Do
       endif
    endif

    if(sw_polynomial_fit == ON) then
       ind=1
       dielectric_tensor0(1:6) = dielectric_tensor_data(1:6,ind)
    else
       ind=0
    end if

!!!    do i=1,nmodes
    jcount = 0
    do i=1, 3*natm_prim
       if (imdtyp((i-1)/3+1) /= OFF) jcount = jcount+1
       do n=1,norder*2
          ind = ind+1
          if (imdtyp((i-1)/3+1) == OFF) cycle

!          dielectric_tensors(1:6,i,n) = dielectric_tensor_data(1:6,ind)
          dielectric_tensors(1:6,jcount,n) = dielectric_tensor_data(1:6,ind)
       end do
    end do

    if(npes>1) then
       call mpi_bcast(dielectric_tensors,6*nmodes*norder*2,mpi_double_precision, &
            &         0,MPI_CommGroup,ierr)
       call mpi_bcast(dielectric_tensor0,6,mpi_double_precision,0,MPI_CommGroup,ierr)
    end if

    if(ipriphonon>=1) write(nfout,*) '<< m_Raman_read_dielectric_tensors >>: Read'

    call flush(nfout)

    return
  end subroutine m_Raman_read_dielec_tensors

  subroutine m_Raman_calc_Mat_dEps_dR          ! @PHONON-GAMMA
    integer :: is,nc,ind,i,j, k, icount
    real(kind=DP) :: ui,x, csum(3,3)
    real(kind=DP), allocatable :: cmat(:,:) !dim(nc,nc)
    real(kind=DP), allocatable :: fvec(:) !dim(nc)

    real(kind=DP), allocatable :: TmpMat(:,:)

!    write(*,*) "nmodes = ", nmodes

    if ( .not. allocated( Mat_dEps_dR ) ) then
       allocate( Mat_dEps_dR( 3, 3, nmodes ) )
    endif

    allocate( TmpMat( 6, nmodes) ) ; TmpMat = 0.0d0

    if(sw_polynomial_fit == ON) then
       nc = norder*2+1
       allocate(cmat(nc,nc)); cmat = 0.d0
       allocate(fvec(nc))
       cmat(1:nc,1) = 1.d0
       ind = 1
       do is = -norder,norder
          if(is==0) cycle
          ind=ind+1
          x = u/dble(is)
          do i=2,nc
             cmat(ind,i) = x**(i-1)
          end do
       end do

       do j=1, 6
          do i=1,nmodes
             fvec(1) = dielectric_tensor0(j)
             fvec(2:nc) = dielectric_tensors(j,i,1:nc-1)
             TmpMat(j,i) = polyfit_diff(nc,cmat,fvec)
          end do
       end do
    else
       ui=dble(norder)/u
       if(printable .and. ipriphonon>=2) then
          write(nfout,*) 'u=',u
          write(nfout,*) 'ui= norder/u =',ui
       endif
       if(norder == 1) then
          ui = 0.5d0*ui
          do j=1, 6
             do i=1,nmodes
                TmpMat(j,i) = ( dielectric_tensors(j,i,2) &
                     &              -dielectric_tensors(j,i,1) )*ui
             end do
          end do
       else if(norder == 2) then
          ui = ui/12.d0
          do j=1, 6
             do i=1,nmodes
                TmpMat(j,i) = ( -8.d0 *dielectric_tensors(j,i,1) &
                     &               +      dielectric_tensors(j,i,2) &
                     &               -      dielectric_tensors(j,i,3) &
                     &               +8.d0 *dielectric_tensors(j,i,4) )*ui
             end do
          end do
       else
          call phase_error_with_msg(nfout,'Your inputed norder is not supported for the differencial approximation.'&
                                   ,__LINE__,__FILE__)
       end if
    end if

    Do i=1, nmodes
       Mat_dEps_dR(1,1,i) = TmpMat(1,i)
       Mat_dEps_dR(2,2,i) = TmpMat(2,i)
       Mat_dEps_dR(3,3,i) = TmpMat(3,i)
       Mat_dEps_dR(1,2,i) = TmpMat(4,i)
       Mat_dEps_dR(1,3,i) = TmpMat(5,i)
       Mat_dEps_dR(2,3,i) = TmpMat(6,i)

       Mat_dEps_dR(2,1,i) = TmpMat(4,i)
       Mat_dEps_dR(3,1,i) = TmpMat(5,i)
       Mat_dEps_dR(3,2,i) = TmpMat(6,i)
    End do

#ifdef USE_ASMS_RAMAN
    call ASMS_Raman_setup( nfout, nframan_spectra, natm, univol, &
     &                     crystal_type, omega_L, &
     &                     photon_k_in, photon_k_out, &
     &                     raman_spectra_hwhm, raman_spectra_freq_pitch, &
     &                     raman_spectra_temperature, &
     &                     nonorthogonal_component, &
     &                     calc_unpolarized_config, &
     &                     num_rot_angles, rot_angles, &
     &                     polar_parallel_or_normal, &
     &                     qdir, dielectric, nonl_op_coeff_rchi2 )
#endif

#ifdef USE_ASMS_RAMAN
    if ( sw_froehlich_correction == ON ) then
       call ASMS_Raman_add_Froelich_term( ipriphonon, Zeff, nmodes, Mat_dEps_dR )
    endif
#endif

! --- accoustic sum rule --
    Do k=1, 3
       csum = 0.0d0
       icount = 0
       Do i=1, natm_prim
          if(imdtyp(i)==OFF) cycle
          icount = icount +1
          j = (icount -1)*3 +k
          csum(:,:) = csum(:,:) + Mat_dEps_dR(:,:,j)
       End do
       csum = csum /dble(natm)

       icount = 0
       Do i=1, natm_prim
          if(imdtyp(i)==OFF) cycle
          icount = icount +1
          j = (icount -1)*3 +k
          Mat_dEps_dR(:,:,j) = Mat_dEps_dR(:,:,j) -csum(:,:)
       End do
    End Do

#if 0
    if ( mype == 0 ) then
       Do i=1, nmodes
          write(*,*) 'BMode : ', i
          write(*,*) Mat_dEps_dR(1,i)
          write(*,*) Mat_dEps_dR(2,i)
          write(*,*) Mat_dEps_dR(3,i)
          write(*,*) Mat_dEps_dR(4,i)
          write(*,*) Mat_dEps_dR(5,i)
          write(*,*) Mat_dEps_dR(6,i)
       End do
    end if
#endif

    deallocate( TmpMat )

  end subroutine m_Raman_calc_Mat_dEps_dR

  subroutine m_Raman_read_optical_coeff
    integer :: i, j, k, n1, n2, n3

    if ( mype /= 0 ) goto 100

    call m_Files_open_nfoptical_coeff( .false. )
    read(nfoptical_coeff,*)
    read(nfoptical_coeff,*)
    Do i=1, 3
       Do j=1, 3
          read(nfoptical_coeff,*) n1, n2, dielectric(i,j)
       End do
    End do
!
    read(nfoptical_coeff,*,end=2)
    read(nfoptical_coeff,*)
    Do i=1, 3
       Do j=1, 3
          Do k=1, 3
             read(nfoptical_coeff,*) n1, n2, n3, nonl_op_coeff_rchi2(i,j,k)
          End do
       End do
    End do

2   close(nfoptical_coeff)

100 continue

    if(npes>1) then
       call mpi_bcast( dielectric, 3*3, mpi_double_precision, &
            &          0, MPI_CommGroup, ierr )
       call mpi_bcast( nonl_op_coeff_rchi2, 3*3*3, mpi_double_precision, &
            &          0, MPI_CommGroup, ierr )
    end if

  end subroutine m_Raman_read_optical_coeff

  subroutine m_Raman_calc_susceptibility( nmodes, omega, modes, &
       &                                  mode_active, mode_irr_rep )
    integer, intent(in) :: nmodes
    real(kind=DP), intent(in) :: omega( nmodes, nqvec )
    complex(kind=CMPLDP), intent(in) :: modes( nmodes, nmodes, nqvec )
    character(len=3), intent(in) :: mode_irr_rep( nmodes, nqvec )
    character(len=4), intent(in) :: mode_active( nmodes, nqvec )

!    if(sw_lo_to_splitting == ON) return

#ifdef USE_ASMS_RAMAN
    if ( mype == 0 ) then
       call m_Files_open_nframan_spectra()
       call ASMS_Raman_calc_susceptibility( nmodes, omega, modes, &
            &                                 mode_active, mode_irr_rep, &
            &                                 Mat_dEps_dR, ionic_mass )
       close( nframan_spectra )
    endif
#else
    call phase_error_with_msg(nfout,"Raman: Not supported kt",__LINE__,__FILE__)
#endif

  end subroutine m_Raman_calc_susceptibility

! ===- utility ===
  function polyfit_diff(nd,mat,vec)
    real(kind=DP) :: polyfit_diff

    integer, intent(in) :: nd
    real(kind=DP), intent(in) :: mat(nd,nd)
    real(kind=DP), intent(in) :: vec(nd)

! local variables
    integer :: ipiv(nd),info,i
    real(kind=DP) :: a(nd,nd),b(nd)
    real(kind=DP) :: x0,df
    real(kind=DP) :: xold
    real(kind=DP) :: eps = 1.d-10

    real(kind=DP) :: ux

    a = mat;   b = vec

    call dgesv(nd,1,a,nd,ipiv,b,nd,info)
    if(info .ne. 0 ) then
       if(printable) write(nfout,*) 'LAPACK routine DGETRS failure: info=',info
       call phase_error_with_msg(nfout,'LAPACK routine DGETRS failure',__LINE__,__FILE__)
    end if

    x0 = 0.d0
    df = diff_polyfunc(nd,b,x0)
    polyfit_diff = df

  end function polyfit_diff

  function polyfunc(nd,b,x)
    real(kind=DP) :: polyfunc

    integer, intent(in) :: nd
    real(kind=DP), intent(in) :: b(nd),x

    integer :: i
    real(kind=DP) :: f

    f = b(1)
    do i=2,nd
       f = f + b(i)*x**(i-1)
    end do
    polyfunc = f

  end function polyfunc

  function diff_polyfunc(nd,b,x)
    real(kind=DP) :: diff_polyfunc

    integer, intent(in) :: nd
    real(kind=DP), intent(in) :: b(nd),x

    integer :: i
    real(kind=DP) :: df

    df = b(2)
    do i=3,nd
       df = df + b(i)*dble(i-1)*x**(i-2)
    end do
    diff_polyfunc = df

  end function diff_polyfunc
! =====

end module m_Raman
