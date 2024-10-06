module m_Crystal_Field
  use m_Control_Parameters,   only : nspin, kimg, noncol, ndim_spinor, sw_hubbard, &
       &                             proj_attribute, SpinOrbit_mode
  use m_Const_Parameters, only : DP, PAI4, PAI, CMPLDP, zi, ON, OFF, Neglected, &
       &                         Hartree, CONST_EV, CONST_kB, PLANCK, speed_of_light
  use m_Crystal_Structure, only : univol, rltv
  use m_Ionic_System,     only : natm, ityp, pos, iatomn, ntyp, ivan, ival, ihubbard
  use m_PseudoPotential,  only : lpsmax, flg_paw, nmesh, mmesh, radr_paw, &
       &                         vlocr_pw, rhcorpw, nloc, psirpw, ltp, nlmt, mtp, taup
  use m_Parallelization,   only : ista_kngp, iend_kngp, MPI_CommGroup, ierr, &
       &                          ista_k, mype, npes, mype

  use m_Charge_Density,   only : chgq_l, hsr
  use m_PlaneWaveBasisSet,  only : kgp, ngabc, gr_l
  use m_Electronic_Structure,  only : vlhxc_l
  use m_XC_Potential,   only : vxc_l
!  use m_Nonlocal_Potential,  only : new_radr_and_wos
  use m_Files,  only : nfpot, nfout, m_Files_open_ps_files, m_Files_close_ps_file

  use m_CD_Mag_Moment,  only : rad_cov_default, rad_cov
  use m_SpinOrbit_Potential,  only : MatU_ylm_RC_L0, MatU_ylm_RC_L1, MatU_ylm_RC_L2, &
       &                             MatU_ylm_RC_L3, MatU_ylm_RC_L4, MatU_ylm_RC_L5, &
       &                             MatU_ylm_RC_L6, Mat_LS_with_real_ylm_L3, &
       &                             m_SO_set_MatU_ylm_RC, &
       &                             m_SO_calc_MatLS_orb_s_to_f, Mat_SOC_Strength
  use m_Orbital_Population,  only : i2lp, om
  use m_FFT,                   only : m_FFT_CD_direct_c, m_FFT_CD_inverse_c, &
       &                              m_FFT_alloc_CD_box, m_FFT_dealloc_CD_box, &
       &                              fft_box_size_CD, nfftp
  use m_Parallelization,  only :  nel_fftp, idisp_fftp, ista_fftp, iend_fftp, mp_fftp, &
       &                          nis_fftp, nie_fftp
  use m_PlaneWaveBasisSet, only : igfp_l
  use m_Ionic_System, only : cps

  use m_FFT, only : m_FFT_coef_CD_integration_kt
  use m_Parallelization,  only :  ista_fftph, iend_fftph
  use mpi

  implicit none
!  include 'mpif.h'

  integer :: max_sph = 49

  integer :: sw_contract_forb = OFF
!  integer :: sw_exclude_vh_selfcharge = ON
  integer :: sw_exclude_vh_selfcharge = OFF
  integer :: sw_exclude_vh_occmat = OFF
  integer :: sw_exclude_vxc = ON

  real(kind=DP), parameter :: Hartree_2_Kelvin = 1.0d0 / CONST_kB  ! 315777D0
  real(kind=DP), parameter :: Hartree_2_cminv = (Hartree * CONST_EV) / (1.0d2 * speed_of_light * PLANCK)  ! 2.194746D5

  real(kind=DP), allocatable :: gaunt_coeff(:,:,:)

contains

  subroutine set_psir_f( has_open_forb, psir_f, guzai_f )
    integer, intent(out) :: has_open_forb(ntyp)
    real(kind=DP), intent(out) :: psir_f( mmesh, ntyp )
    real(kind=DP), intent(out) :: guzai_f( ntyp )

    integer :: it, ierr, nfp, i, ir
    integer :: num_core_ae_wfns
    integer, allocatable :: qnum_n_core_ae_wfns(:)
    integer, allocatable :: qnum_l_core_ae_wfns(:)
    real(kind=DP), allocatable :: psir_core_ae_wfns(:,:)
    real(kind=DP), allocatable :: enelevel_core_ae_wfns(:)
    real(kind=DP), allocatable :: focc_core_ae_wfns(:)
    real(kind=DP), allocatable :: guzai_core_ae_wfns(:)
    real(kind=DP) :: weight

! == KT_DEBUG === 2015/06/15 ==
!    call m_Files_open_ps_file(ivan,iatomn,ntyp,it,ierr)
!                                      ! this does not work properly by unknown reason.
    call m_Files_open_ps_files(ivan,iatomn,ntyp,ierr)
! =============== 2015/06/15

    if (ierr/=0) call mpi_stop(nfout)

    has_open_forb = 0

    Do it=1, ntyp
       nfp = nfpot(it)

       if ( mype == 0 ) call read_num_core_ae_wfns( nfp, num_core_ae_wfns, it )
       if ( npes > 1 ) then
          call mpi_bcast( num_core_ae_wfns, 1, mpi_integer, 0, mpi_comm_world, ierr )
       endif

       if ( num_core_ae_wfns == 0 ) cycle

       allocate( qnum_n_core_ae_wfns( num_core_ae_wfns ) )
       allocate( qnum_l_core_ae_wfns( num_core_ae_wfns ) )
       allocate( psir_core_ae_wfns( nmesh(it),  num_core_ae_wfns ) )
       allocate( enelevel_core_ae_wfns( num_core_ae_wfns ) )
       allocate( focc_core_ae_wfns( num_core_ae_wfns ) )
       allocate( guzai_core_ae_wfns( num_core_ae_wfns ) )

       if ( mype == 0 ) then
          call read_data_core_ae_wfns( nfp, num_core_ae_wfns, nmesh(it), &
               &                       qnum_n_core_ae_wfns, qnum_l_core_ae_wfns, &
               &                       psir_core_ae_wfns, &
               &                       enelevel_core_ae_wfns, focc_core_ae_wfns )
          call read_data_core_soc( nfp, it, num_core_ae_wfns, &
               &                   qnum_n_core_ae_wfns, qnum_l_core_ae_wfns, &
               &                   guzai_core_ae_wfns )
       endif

       if ( npes > 1 ) then
          call mpi_bcast( qnum_n_core_ae_wfns, num_core_ae_wfns, &
               &          mpi_integer, 0, mpi_comm_world, ierr )
          call mpi_bcast( qnum_l_core_ae_wfns, num_core_ae_wfns, &
               &          mpi_integer, 0, mpi_comm_world, ierr )
          call mpi_bcast( enelevel_core_ae_wfns, num_core_ae_wfns, &
               &          mpi_double_precision, 0, mpi_comm_world, ierr )
          call mpi_bcast( focc_core_ae_wfns, num_core_ae_wfns, &
               &          mpi_double_precision, 0, mpi_comm_world, ierr )
          call mpi_bcast( psir_core_ae_wfns, nmesh(it)*num_core_ae_wfns, &
               &          mpi_double_precision, 0, mpi_comm_world, ierr )
          call mpi_bcast( guzai_core_ae_wfns, num_core_ae_wfns, &
               &          mpi_double_precision, 0, mpi_comm_world, ierr )
       endif

!
       Do i=1, num_core_ae_wfns
          if ( qnum_l_core_ae_wfns(i) == 3 ) then
             if ( focc_core_ae_wfns(i) > 0.1 .and. focc_core_ae_wfns(i) < 13.9 ) then
                has_open_forb(it) = 1
                Do ir=1, nmesh(it)
                   psir_f(ir,it) = psir_core_ae_wfns(ir,i)
                End Do
             endif
             guzai_f(it) = guzai_core_ae_wfns(i)
          endif
       End Do

       deallocate( qnum_n_core_ae_wfns )
       deallocate( qnum_l_core_ae_wfns )
       deallocate( psir_core_ae_wfns )
       deallocate( enelevel_core_ae_wfns )
       deallocate( focc_core_ae_wfns )
       deallocate( guzai_core_ae_wfns )

       call m_Files_close_ps_file(it)
    End Do
  end subroutine set_psir_f

  subroutine read_num_core_ae_wfns( nfp, nums, it )
    integer, intent(in) :: nfp, it
    integer, intent(out) :: nums

    integer :: length, ierr
    character(30) :: search_key

    nums = 0

    search_key = "CORE STATES";  length = len(search_key)
    call read_size_of_array_from_pp(nfp, nums, length, search_key, ierr)

    if ( ierr /= 0 ) then
       write(nfout,*) '----------------------'
       write(nfout,'(A,I2)') '!!! Keyword CORE STATES is not found in the PP', it
       write(nfout,*) '----------------------'
    endif

  end subroutine read_num_core_ae_wfns

  subroutine read_data_core_ae_wfns( nfp, num_core_ae_wfns, nmesh, &
       &                             qnum_n, qnum_l, psir_core, ene_level, focc )
    implicit none

    integer, intent(in) :: num_core_ae_wfns
    integer, intent(in) :: nmesh, nfp
    integer, intent(out) :: qnum_l(num_core_ae_wfns)
    integer, intent(out) :: qnum_n(num_core_ae_wfns)
    real(kind=8), intent(out) :: psir_core(nmesh, num_core_ae_wfns), &
         &                       ene_level(num_core_ae_wfns), focc(num_core_ae_wfns)

    integer :: i, n1, l1, k
    real(kind=DP) :: ene1, f1

    Do i=1, num_core_ae_wfns
       read(nfp,*) n1, l1, ene1, f1
       read(nfp,*) (psir_core(k,i),k=1,nmesh)
       !
       qnum_n(i) = n1;  qnum_l(i) = l1;  ene_level(i) = ene1;   focc(i) = f1
    End Do
  end subroutine read_data_core_ae_wfns

  subroutine read_size_of_array_from_pp( nfp, size_of_array, length, search_tag, ierr )
    implicit none

    integer, intent(in) :: nfp, length
    integer, intent(out) :: size_of_array, ierr
    character(length), intent(in) :: search_tag

    integer :: ifound
    character(30) :: line1

    size_of_array = 0;  ierr = 0

    Do while (.true.)
       read(nfp,'(a30)',end=10) line1
       ifound = index( line1, search_tag )
       if ( ifound /= 0 ) goto 20
    End do

10  ierr = 1; return

20  read(nfp,*) size_of_array

  end subroutine read_size_of_array_from_pp

  subroutine read_data_core_soc( nfp, it, num_core_ae_wfns, &
       &                         qnum_n, qnum_l, guzai_core )
    implicit none
    integer, intent(in) :: nfp, it, num_core_ae_wfns
    integer, intent(in) :: qnum_n(num_core_ae_wfns)
    integer, intent(in) :: qnum_l(num_core_ae_wfns)
    real(kind=DP), intent(out) :: guzai_core(num_core_ae_wfns)

    integer :: ifound, ier
    integer :: count, n1, l1, i, j
    real(kind=DP) :: c1
    character(30) :: search_tag
    character(30) :: line1

    ierr = 0;    guzai_core = 0.0d0

    search_tag = "SOC-CORE"

    Do while (.true.)
       read(nfp,'(a30)',end=10) line1
       ifound = index( line1, trim(search_tag) )
       if ( ifound /= 0 ) goto 20
    End do

10  ierr = 1;
    write(nfout,*) '----------------------'
    write(nfout,'(A,I2)') '!!! Keyword SOC-CORE is not found in the PP', it
    write(nfout,*) '----------------------'
    return

20  continue

    read(nfp,*) count
    Do i=1, count
       read(nfp,*) n1, l1, c1
       Do j=1, num_core_ae_wfns
          if ( qnum_n(j) == n1 .and. qnum_l(j) == l1 ) then
             guzai_core(j) = c1
             exit
          endif
       End Do
    End Do
  end subroutine read_data_core_soc

  subroutine contract_f_orbital( psir )
    real(kind=DP), intent(inout) :: psir(mmesh,ntyp)

    integer :: ik, ilmt1, il1, tau1, ir, ier, it
    real(kind=DP) :: csum, factor, rcut, beta
    real(kind=DP), allocatable :: tmp_fn(:), wos(:)

    allocate( tmp_fn(mmesh) );  allocate( wos(mmesh) );
    beta = 0.5d0

    Do it=1, ntyp
       rcut = rad_cov_default( nint(iatomn(it)) )
                    ! Revised according to a report from ASMS Co.ltd, 10 March 2016.

       call set_weight_exp( ier, 1, nmesh(it), radr_paw(:,it), wos )

       tmp_fn = 0.0d0
       Do ir=1, nmesh(it)
          tmp_fn(ir) = psir(ir,it) *exp( -radr_paw(ir,it)/rcut)
!          tmp_fn(ir) = psir(ir,it) /( 1.0d0 +exp( beta *(radr_paw(ir,it)-rcut) ) )
       End do

       csum = 0.0d0
       Do ir=1, nmesh(it)
          csum = csum + wos(ir) *tmp_fn(ir)**2
       End do

       factor = 1.0d0 /sqrt(csum)
       psir(:,it) = 0.0d0

       do ir = 1,nmesh(it)
          psir(ir,it) = tmp_fn(ir) *factor
       end do
    End Do
    call mpi_barrier( MPI_CommGroup, ierr )
    deallocate(tmp_fn); deallocate(wos)

  end subroutine contract_f_orbital

  subroutine m_CF_calc_CF_param
    integer :: val_l, lcmax, nsph1, nsph2, nsph3, nr
    integer :: ia, ig, it, excl_flg
    integer :: lun1, lun2
    integer, allocatable :: has_open_forb(:), ilist(:)
    real(kind=DP), allocatable :: qx(:), qy(:), qz(:), psir_f(:,:), guzai_f(:)
    real(kind=DP), allocatable :: zfcos(:), zfsin(:), wos(:)
    real(kind=DP), allocatable :: v_lm(:), v_lm_felec(:)
    real(kind=DP), allocatable :: chgq_self(:,:)
    real(kind=DP) :: rcut, norm_wf2
    real(kind=DP), allocatable :: r_avg(:)

    if ( .not. flg_paw ) return

    lun1 = 1000;   lun2 = 1001

    val_l = 3;     lcmax = 6

    if ( mype == 0 ) then
       open(lun1,file="Crystal_Field.input",status="old",form="formatted")
       read(lun1,*) sw_exclude_vxc
       read(lun1,*) sw_exclude_vh_selfcharge
       read(lun1,*) sw_exclude_vh_occmat
       close(lun1)
    endif
    call mpi_bcast( sw_exclude_vxc, 1, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( sw_exclude_vh_selfcharge, 1, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( sw_exclude_vh_occmat, 1, mpi_integer, 0, mpi_comm_world, ierr )

    call m_CF_set_gaunt_coeff
    if ( noncol ) then
       call m_SO_set_MatU_ylm_RC
       call m_SO_calc_MatLS_orb_s_to_f
    endif
!
    allocate( ilist(max_sph) ); ilist = 0

    Do nsph1=(val_l)**2+1, (val_l+1)**2
       Do nsph3=(val_l)**2+1, (val_l+1)**2
          Do nsph2=1, max_sph
             if ( abs(gaunt_coeff(nsph1,nsph2,nsph3)) > 1.0D-6 ) then
                ilist(nsph2) = 1
             endif
          End Do
       ENd Do
    End Do

    allocate( qx(ista_kngp:iend_kngp) )
    allocate( qy(ista_kngp:iend_kngp) )
    allocate( qz(ista_kngp:iend_kngp) )

    allocate( zfcos(ista_kngp:iend_kngp) )
    allocate( zfsin(ista_kngp:iend_kngp) )

    Do ig=ista_kngp, iend_kngp
       qx(ig) = rltv(1,1)*ngabc(ig,1) +rltv(1,2)*ngabc(ig,2) +rltv(1,3)*ngabc(ig,3)
       qy(ig) = rltv(2,1)*ngabc(ig,1) +rltv(2,2)*ngabc(ig,2) +rltv(2,3)*ngabc(ig,3)
       qz(ig) = rltv(3,1)*ngabc(ig,1) +rltv(3,2)*ngabc(ig,2) +rltv(3,3)*ngabc(ig,3)
    End do

    allocate( guzai_f(ntyp) ); guzai_f = 0.0d0
    allocate( psir_f(mmesh,ntyp) ); psir_f = 0.0d0
    allocate( has_open_forb(ntyp) );  has_open_forb = 0
    call set_radial_wavefunction( has_open_forb, psir_f, guzai_f )

    allocate( v_lm(max_sph) );  v_lm = 0.0d0

    if ( sw_exclude_vh_selfcharge == ON ) then
       allocate( chgq_self(ista_kngp:iend_kngp,kimg) );
    endif
    if ( sw_exclude_vh_occmat == ON ) then
       allocate( v_lm_felec(max_sph) )
    endif

    if ( mype == 0 ) call file_open( lun1, lun2 )

    Do ia=1, natm
       it = ityp(ia)
       if ( has_open_forb(it) == 0 .and. lpsmax(it) < 4 ) cycle

       rcut = rad_cov_default( nint(iatomn(it)) )
       call set_range_for_radial_integration( it, rcut, nr )

       allocate( wos(nr) )
       call set_weight_exp( ierr, 1, nr, radr_paw(:,it), wos )
       call calc_norm_wf2( nr, wos, norm_wf2 )

       call calc_phase2( natm, pos, ia, kgp, ngabc, ista_kngp, iend_kngp, &
            &            zfcos, zfsin )

       excl_flg = 0
!       if ( sw_exclude_vh_selfcharge == ON .and. has_open_forb(it) ==0 ) then
       if ( sw_exclude_vh_selfcharge == ON ) then
          excl_flg = 1
          call calc_chgq_inside_rcut( ia, rcut, zfcos, zfsin, chgq_self )
       endif

       call calc_v_lm( it, nr, wos, norm_wf2, zfcos, zfsin, v_lm, excl_flg )

       if ( sw_hubbard == ON .and. sw_exclude_vh_occmat == ON ) then
          call calc_v_lm_felec( ia, nr, v_lm_felec )
          v_lm = v_lm -v_lm_felec
       endif

       if ( mype == 0 ) then
          call print_field_param_A_real( ia, lun1, v_lm )
          call print_field_param_B_cmplx( ia, lun2, v_lm )
       endif

       call calc_hamil_eigenval( ia, it, has_open_forb, v_lm, guzai_f )

       deallocate( wos )
    ENd Do

    if ( mype == 0 ) call file_close( lun1, lun2 )

    deallocate(qx);    deallocate(qy);    deallocate(qz);
    deallocate(zfcos); deallocate(zfsin)
    deallocate(psir_f);deallocate(has_open_forb)
    deallocate(ilist); deallocate(v_lm);
    deallocate(guzai_f)

    if (allocated(chgq_self))   deallocate( chgq_self )
    if (allocated(v_lm_felec)) deallocate( v_lm_felec )

  contains

    subroutine file_open( lun1, lun2 )
      integer, intent(in) :: lun1, lun2

      character*32 :: file1, file2

      file1="Crystal_Field_Param_A.dat"
      open( lun1, file=file1, status="unknown", form="formatted")
      write(lun1,'(A)') "# Crystal Field Parameters (real)"

      file2="Crystal_Field_Param_B.dat"
      open( lun2, file=file2, status="unknown", form="formatted")
      write(lun2,'(A)') "# Crystal Field Parameters (complx)"
    end subroutine file_open

    subroutine file_close( lun1, lun2 )
      integer, intent(in) :: lun1, lun2

      close(lun1);  close(lun2)
    end subroutine file_close

    subroutine set_range_for_radial_integration( it, rcut, nr )
      integer, intent(in) :: it
      real(kind=DP), intent(in) :: rcut
      integer, intent(out) :: nr

      integer :: ir

      nr = 0
      Do ir=1, nmesh(it)
         if ( radr_paw(ir,it) < rcut ) nr = nr +1
      End Do
!       nr = nmesh(it)
    end subroutine set_range_for_radial_integration

    subroutine set_radial_wavefunction( has_open_forb, psir_f, guzai_f )
      integer, intent(out) :: has_open_forb(ntyp)
      real(kind=DP), intent(inout) :: psir_f(mmesh,ntyp), guzai_f(ntyp)

      integer :: it

      call set_psir_f( has_open_forb, psir_f, guzai_f )
      Do it=1, ntyp
         if ( lpsmax(it) == 4 ) then
            psir_f(1:nmesh(it),it) = psirpw(1:nmesh(it),4,1,it)
         endif
      End Do
      if ( sw_contract_forb == ON ) call contract_f_orbital( psir_f )
    end subroutine set_radial_wavefunction

    subroutine calc_norm_wf2( nr, wos, norm_wf2 )
      integer, intent(in) :: nr
      real(kind=DP), intent(in) :: wos(nr)
      real(kind=DP), intent(out) :: norm_wf2

      integer :: ir
      real(kind=DP) :: csum, c1

      csum = 0.0d0
      do ir = 1, nr
         c1 = wos(ir) *psir_f(ir,it)**2
         csum = csum +c1
      End do
      norm_wf2 = csum

    end subroutine calc_norm_wf2

    subroutine calc_r_avg( nr, wos, norm_wf2, r_avg )
      integer, intent(in) :: nr
      real(kind=DP), intent(in) :: norm_wf2
      real(kind=DP), intent(in) :: wos(nr)
      real(kind=DP), intent(out) :: r_avg(0:lcmax)

      integer :: il, ir
      real(kind=DP) :: csum, c1

      Do il=0, lcmax
         csum = 0.0d0
         do ir = 1, nr
            c1 = wos(ir) *psir_f(ir,it)**2 *radr_paw(ir,it)**il
            csum = csum +c1
         End do
         r_avg(il) = csum /norm_wf2
      End Do
    end subroutine calc_r_avg

    subroutine calc_chgq_inside_rcut( ia, rcut, zfcos, zfsin, chgq_filtered )
      integer, intent(in) :: ia
      real(kind=DP), intent(in) :: rcut
      real(kind=DP), intent(in) :: zfcos(ista_kngp:iend_kngp)
      real(kind=DP), intent(in) :: zfsin(ista_kngp:iend_kngp)
      real(kind=DP), intent(out) :: chgq_filtered(ista_kngp:iend_kngp,kimg)

      real(kind=DP), allocatable :: afft(:), bfft(:)
      real(kind=DP), allocatable :: wkq_l(:,:)

      allocate(afft(ista_fftp:iend_fftp)); afft =0.0d0
      allocate(bfft(ista_fftp:iend_fftp)); bfft =0.0d0

      call m_FFT_alloc_CD_box()

      allocate( wkq_l(ista_kngp:iend_kngp,kimg) )
      if ( noncol ) then
         wkq_l(ista_kngp:iend_kngp,:) = chgq_l(ista_kngp:iend_kngp,:,1)
      else
         if ( nspin == 1 ) then
            wkq_l(ista_kngp:iend_kngp,:) = chgq_l(ista_kngp:iend_kngp,:,1)
         else
            wkq_l(ista_kngp:iend_kngp,:) = chgq_l(ista_kngp:iend_kngp,:,1) &
                 &                           +chgq_l(ista_kngp:iend_kngp,:,2)
         endif
      endif
      call set_mat_on_fftmesh_rspace( wkq_l, afft )
      call set_filter( ia, rcut, zfcos, zfsin, wkq_l )
      call set_mat_on_fftmesh_rspace( wkq_l, bfft )
      call product_on_fft_mesh( afft, bfft )
      call set_chg( bfft, wkq_l )

      chgq_filtered = wkq_l
      deallocate( afft );    deallocate( bfft );    deallocate( wkq_l )
      call m_FFT_dealloc_CD_box()
    end subroutine calc_chgq_inside_rcut

    subroutine set_chg( afft, wkq_l )
      real(kind=DP), intent(inout) :: afft(ista_fftp:iend_fftp)
      real(kind=DP), intent(out) :: wkq_l(ista_kngp:iend_kngp,kimg)

      integer :: ik, i, ip
      real(kind=DP) :: rinplw
      real(kind=DP), allocatable :: afft_mpi1(:)

      call m_FFT_CD_direct_c( nfout, afft )        ! R -> G

      wkq_l = 0.0d0

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
            wkq_l(i,ik) = afft_mpi1(ip)*rinplw
         end do
      end do
      deallocate(afft_mpi1)

    end subroutine set_chg

    subroutine product_on_fft_mesh( afft, bfft )
      real(kind=DP), intent(in) :: afft(ista_fftp:iend_fftp)
      real(kind=DP), intent(inout) :: bfft(ista_fftp:iend_fftp)
      integer :: ig
      real(kind=DP) :: c1, c2

      Do ig=ista_fftp, iend_fftp, 2
#if 1
         c1 = afft(ig) *bfft(ig)
         bfft(ig) = c1
         bfft(ig+1) = 0.0d0
#else
         c1 = afft(ig) *bfft(ig) -afft(ig+1) *bfft(ig+1)
         c2 = afft(ig) *bfft(ig+1) +afft(ig+1) *bfft(ig)
         bfft(ig) = c1
         bfft(ig+1) = c2
#endif
      End Do
    end subroutine product_on_fft_mesh

    subroutine set_filter( ia, rcut, zfcos, zfsin, filter_l )
      integer, intent(in) :: ia
      real(kind=DP), intent(in) :: rcut
      real(kind=DP), intent(in) :: zfcos(ista_kngp:iend_kngp)
      real(kind=DP), intent(in) :: zfsin(ista_kngp:iend_kngp)
      real(kind=DP), intent(out) :: filter_l(ista_kngp:iend_kngp,kimg)

      integer :: ig
      real(kind=DP) :: c1, gx, gy, gz, fac, f1
      real(kind=DP), allocatable :: wka(:), wkb(:)

      fac = PAI4 *rcut**3 /univol

      allocate( wka(ista_kngp:iend_kngp) )
      allocate( wkb(ista_kngp:iend_kngp) )
      Do ig=ista_kngp, iend_kngp
         wka(ig) = gr_l(ig) *rcut
      End do

      call dsjnv( 1, iend_kngp -ista_kngp +1, wka, wkb )

      Do ig=ista_kngp, iend_kngp
         if ( ig==1 ) then
            f1 = fac /3.0d0
            filter_l(ig,1) =  f1
            filter_l(ig,2) =  0.0d0
         else
            f1 =  fac *wkb(ig) /wka(ig)      ! j1(GR)/GR
            filter_l(ig,1) =  f1 *zfcos(ig)
            filter_l(ig,2) = -f1 *zfsin(ig)
         endif
      End Do
      deallocate( wka );      deallocate( wkb )

    end subroutine set_filter

    subroutine set_mat_on_fftmesh_rspace( wkq_l, afft )
      real(kind=DP), intent(in) :: wkq_l(ista_kngp:iend_kngp,kimg)
      real(kind=DP), intent(out) :: afft(ista_fftp:iend_fftp)

      integer :: j, i, ip
      real(kind=DP) :: rinplw
      real(kind=DP), allocatable :: afft_mpi1(:), afft_mpi2(:), afft_mpi3(:)

      rinplw = 1.d0 /product(fft_box_size_CD(1:3,1))

      afft = 0.0d0;

      allocate(afft_mpi1(nfftp)); afft_mpi1 = 0.0d0
      if(npes >= 2) then
         allocate(afft_mpi2(mp_fftp)); afft_mpi2 = 0.0d0
         allocate(afft_mpi3(mp_fftp)); afft_mpi3 = 0.0d0
      end if
      do j = 1, kimg
         do i = ista_kngp, iend_kngp  !for mpi
            ip = (igfp_l(i)-1)*kimg + j
            afft_mpi1(ip) = wkq_l(i,j)
         end do
      end do

      if (npes >= 2) then
         call mpi_barrier(MPI_CommGroup,ierr)
         do j = 0, npes-1
            do i = nis_fftp(j),nie_fftp(j)
               afft_mpi2(i-nis_fftp(j)+1) = afft_mpi1(i)
            end do
            call mpi_allreduce(afft_mpi2,afft_mpi3,mp_fftp &
                 &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
            if(j == mype) then
               do i = ista_fftp, iend_fftp
                  afft(i) = afft_mpi3(i - ista_fftp + 1)
               end do
            end if
         end do
      else
         afft = afft_mpi1
      end if

      call m_FFT_CD_inverse_c(nfout,afft)        ! G-->R space
      deallocate(afft_mpi1)
      if (allocated(afft_mpi2) ) deallocate(afft_mpi2)
      if (allocated(afft_mpi3) ) deallocate(afft_mpi3)

    end subroutine set_mat_on_fftmesh_rspace

    subroutine calc_v_lm( it, nr, wos, norm_wf2, zfcos, zfsin, v_lm, excl_flg )
      integer, intent(in) :: it, nr, excl_flg
      real(kind=DP), intent(in) :: wos(nr)
      real(kind=DP), intent(in) :: norm_wf2
      real(kind=DP), intent(in) :: zfcos(ista_kngp:iend_kngp)
      real(kind=DP), intent(in) :: zfsin(ista_kngp:iend_kngp)
      real(kind=DP), intent(out) :: v_lm(max_sph)

      integer :: il1, im1, ir, ig, is, nsph1
      integer :: ierr, ismax
      real(kind=DP) :: fac, csum2, c1, csum, csum_mpi, facr, ctmp

      real(kind=DP), allocatable :: wka(:), wkb(:), ylm(:), snl2(:)

      fac = PAI4
      v_lm = 0.0d0
      ! ---
      allocate( wka(ista_kngp:iend_kngp) )
      allocate( wkb(ista_kngp:iend_kngp) )
      allocate( ylm(ista_kngp:iend_kngp) )
      allocate( snl2(ista_kngp:iend_kngp) )

      Do il1=0, lcmax
         Do im1=1, 2*il1 +1
            nsph1 = il1 **2 +im1
            if ( ilist(nsph1) == 0 ) cycle
            call sphr_general( iend_kngp -ista_kngp +1, nsph1, &
                 &     qx(ista_kngp:iend_kngp), &
                 &     qy(ista_kngp:iend_kngp), &
                 &     qz(ista_kngp:iend_kngp), ylm(ista_kngp:iend_kngp) )
            !
            snl2 = 0.0d0
            do ir = 1, nr
               facr = wos(ir) *psir_f(ir,it) *psir_f(ir,it)
               Do ig=ista_kngp, iend_kngp
                  wka(ig) = gr_l(ig) *radr_paw(ir,it)
               End do
               call dsjnv( il1, iend_kngp -ista_kngp +1, &
                    &      wka(ista_kngp:iend_kngp), &
                    &      wkb(ista_kngp:iend_kngp) )
               Do ig=ista_kngp, iend_kngp
                  snl2(ig) = snl2(ig) + fac *facr *wkb(ig) *ylm(ig)
               End do
            end do

            c1 = zi **il1;         snl2 = snl2 *c1

            csum = 0.0d0

            if ( noncol ) then
               ismax = 1
               Do is=1, ismax
                  Do ig=ista_kngp, iend_kngp
                     if ( sw_exclude_vxc == ON ) then
                        c1 = ( vlhxc_l(ig,1,1) -vxc_l(ig,1,1) )*zfcos(ig) &
                             & -(vlhxc_l(ig,2,is) -vxc_l(ig,2,is) )*zfsin(ig)
                     else
                        c1 = ( vlhxc_l(ig,1,is) ) *zfcos(ig) &
                             & -( vlhxc_l(ig,2,is) ) *zfsin(ig)
                     endif
                     csum = csum +c1 *snl2(ig)
                  End do
               End do
            else
               ismax = nspin
               Do is=1, ismax
                  if ( kimg ==2 ) then
                     Do ig=ista_kngp, iend_kngp
                        if ( sw_exclude_vxc == ON ) then
                           c1 = ( vlhxc_l(ig,1,is) -vxc_l(ig,1,is) ) *zfcos(ig) &
                                & -( vlhxc_l(ig,2,is) -vxc_l(ig,2,is) ) *zfsin(ig)
                        else
                           c1 = ( vlhxc_l(ig,1,is) ) *zfcos(ig) &
                                & -( vlhxc_l(ig,2,is) ) *zfsin(ig)
                        endif
                        csum = csum +c1 *snl2(ig)
                     End Do
                  else
                     Do ig=ista_kngp, iend_kngp
                        if ( sw_exclude_vxc == ON ) then
                           c1 = ( vlhxc_l(ig,1,is) -vxc_l(ig,1,is) )*zfcos(ig)
                        else
                           c1 = vlhxc_l(ig,1,is)*zfcos(ig)
                        endif
                        csum = csum +c1 *snl2(ig)
                     End Do
                  endif
               End do
               csum = csum /dble(nspin)
            endif

            if ( excl_flg == 1 ) then     ! kimg == 2 only
               Do ig=ista_kngp, iend_kngp
                  if ( ig /=1 ) then
                     c1 = ( chgq_self(ig,1) *zfcos(ig) -chgq_self(ig,2) *zfsin(ig) )&
                          &  /gr_l(ig)**2 *PAI4
                     csum = csum -c1 *snl2(ig)
                  endif
               End Do
            endif
            call mpi_allreduce( csum, csum_mpi, 1, mpi_double_precision, mpi_sum, &
                 &              MPI_CommGroup, ierr )
            !
            v_lm(nsph1) = csum_mpi
         end Do
      end Do
      v_lm = v_lm /norm_wf2
      deallocate(wka);      deallocate(wkb);      deallocate(ylm)
      deallocate(snl2)

    end subroutine calc_v_lm

    subroutine calc_v_lm_felec( ia, nr, v_lm_felec )
      integer, intent(in) :: ia, nr
      real(kind=DP), intent(out) :: v_lm_felec(max_sph)

      integer :: ih, ie, il1, il2, il3, immax, is, ir
      integer :: im1, im2, im3, nsph1, nsph2, nsph3
      real(kind=DP) :: c1, c2, csum, coeff, c3, c4
      real(kind=DP), allocatable :: rho(:), vh(:)
!
      integer :: ilmt1, ilmt2
      real(kind=DP) :: wk1(1), wk2(1), wk3(1), ylm(1)
!
      csum = 0.0
      it = ityp(ia)
      Do ir=1, nmesh(it)
         csum = csum +rhcorpw(ir,it) *wos(ir)
      End Do
!      write(*,*) "it csum = ", csum, csum/PAI4

      csum = 0.0d0
!
      v_lm_felec = 0.0d0

      ih = ihubbard(ia)
      if (ih == 0 ) return
!
      ie = proj_attribute(ih)%ielem
!      it=ityp(ia)
      il1 = proj_attribute(ih)%l

      if ( il1 /= 3 ) return

      immax = i2lp(ih)
      immax = 3*2 +1

      allocate( rho(nmesh(it)) );      allocate( vh(nmesh(it)) )

      Do il2=0, lcmax
         Do im2=1, 2*il2 +1
            nsph2 = il2**2 +im2
            if ( ilist(nsph2) == 0 ) cycle

            call real_spherical_harmonics( nsph2, 1, wk1, wk2, wk3, ylm, c3 )

            rho = 0.d0;     coeff = 0.0d0

            Do im1=1, immax
               Do im3=1, immax
                  c1 = 0.0d0
                  if ( noncol ) then
                     c1 = c1 +om(im1,im3,ie,ia,1)
                  else
                     Do is=1, nspin
                        c1 = c1 +om(im1,im3,ie,ia,is)
                     End Do
                     if ( nspin == 1 ) c1 = c1 *2.0d0
                  endif
!
!                  c1 = 0.d00
!                  if ( im1 == im3 ) then
!                     c1 = 3.0d0 /dble(immax)
!                  endif
!
                  nsph1 = il1**2 +im1;    nsph3 = il1**2 +im3
                  c2 = gaunt_coeff( nsph1, nsph2, nsph3 )
                  if ( abs(c2) < 1.0D-6 ) cycle
                  coeff = coeff +c1 *c2
               End do
            End Do

            Do ir=1, nmesh(it)
               rho(ir) = coeff *psir_f(ir,it)**2 /radr_paw(ir,it)**2
            End Do
!
            if ( nsph2 == 1 ) rho = rho /PAI4
!            rho = rho *c3**2

            vh = 0.0d0
            call calc_vh_asym( nmesh(it), radr_paw(:,it), rho, vh, il2 )

            csum = 0.0d0
            Do ir=1, nr
               csum = csum +wos(ir) *psir_f(ir,it)**2 *vh(ir)
            End Do
            v_lm_felec(nsph2) = csum
         End Do
      End Do
      deallocate(rho); deallocate(vh)
    end subroutine calc_v_lm_felec

    subroutine print_field_param_A_real( ia, lun, v_lm )
      integer, intent(in) :: lun, ia
      real(kind=DP), intent(in) :: v_lm(max_sph)

      integer :: il1, im1, im2, nsph1
      real(kind=DP) :: c1, c2
      real(kind=DP) :: wk1(1), wk2(1), wk3(1), ylm(1), coeff

      wk1(1) = 1.0;    wk2(1) = 1.0;      wk3(1) = 1.0

      Do il1=0, lcmax
         Do im1=1, 2*il1 +1
            nsph1 = il1 **2 +im1
            if ( ilist(nsph1) == 0 ) cycle
            call real_spherical_harmonics( nsph1, 1, wk1, wk2, wk3, ylm, coeff )

            im2 = get_m_tesseral_harmonics( il1, im1 )

            c1 = coeff *Hartree_2_Kelvin *v_lm(nsph1)
            write(lun,'(A,I4,A,I3,A,I3,A,F20.5,A)') &
                 &      "ia=", ia, ", L= ", il1, ", M= ", im2, &
                 &      ",   A_LM <r^L>=", c1, " [K]"
         End Do
      End Do
    end subroutine print_field_param_A_real

    subroutine print_field_param_B_cmplx( ia, lun, v_lm )
      integer, intent(in) :: lun, ia
      real(kind=DP), intent(in) :: v_lm(max_sph)

      integer :: il1, im1, nsph1, im2
      complex(kind=CMPLDP) :: z1, zsum

      Do il1=0, lcmax
         Do im1=-il1, il1
            zsum = 0.0d0
            Do im2=1, 2*il1 +1
               nsph1 = il1 **2 +im2
               if ( ilist(nsph1) == 0 ) cycle
               if ( il1 == 0 ) then
                  zsum = zsum +v_lm(nsph1) * MatU_ylm_RC_L0(im2,im1)
               else if ( il1 == 1  ) then
                  zsum = zsum +v_lm(nsph1) * MatU_ylm_RC_L1(im2,im1)
               else if ( il1 == 2 ) then
                  zsum = zsum +v_lm(nsph1) * MatU_ylm_RC_L2(im2,im1)
               else if ( il1 == 3 ) then
                  zsum = zsum +v_lm(nsph1) * MatU_ylm_RC_L3(im2,im1)
               else if ( il1 == 4 ) then
                  zsum = zsum +v_lm(nsph1) * MatU_ylm_RC_L4(im2,im1)
               else if ( il1 == 5 ) then
                  zsum = zsum +v_lm(nsph1) * MatU_ylm_RC_L5(im2,im1)
               else if ( il1 == 6 ) then
                  zsum = zsum +v_lm(nsph1) * MatU_ylm_RC_L6(im2,im1)
               endif
            End Do
            z1 = sqrt( (2*il1 +1)/PAI4 ) *zsum *Hartree_2_cminv

            write(lun,'(A,I4,A,I3,A,I3,A,2F18.5,A)') &
                 &      "ia=", ia, ", L= ", il1, ", M= ", im1, &
                 &      ",  B_LM =", real(z1), aimag(z1), " [cm-1]"
         End Do
      End Do
    end subroutine print_field_param_B_cmplx

    subroutine calc_hamil_eigenval( ia, it, has_open_forb, v_lm, guzai_f )
      integer, intent(in) :: ia, it
      integer, intent(in) :: has_open_forb(ntyp)
      real(kind=DP), intent(in) :: v_lm(max_sph), guzai_f(ntyp)

      integer :: im1, im3, nsph1, nsph2, nsph3, immax, is1, is3, istmp, size2
      integer :: lwork, info
      real(kind=DP) :: csum, gzi
      real(kind=DP), allocatable :: eigenval(:), rwork(:)
      complex(kind=CMPLDP), allocatable :: matH(:,:), work(:)

      immax = 2*val_l +1
      size2 = immax *ndim_spinor
      allocate( matH(size2,size2) ); matH = 0.0d0

      Do im1=1, immax
         Do im3=1, immax
            nsph1 = val_l**2 +im1;    nsph3 = val_l**2 +im3

            csum = 0.0d0
            Do nsph2=1, max_sph
               if ( ilist(nsph2) == 0 ) cycle
!               if ( nsph2 /= 5 ) cycle
               csum = csum +v_lm(nsph2) *gaunt_coeff(nsph1,nsph2,nsph3)
            End Do
            Do is1=1, ndim_spinor
               matH(immax*(is1-1)+im1,immax*(is1-1)+im3) = csum
            End Do
         End  Do
      End Do
!
      if ( noncol ) then
         gzi = 0.0d0
         if ( has_open_forb(it) == 1 ) then
            gzi = guzai_f(it)
         else if ( SpinOrbit_mode /= Neglected ) then
            gzi = Mat_Soc_Strength( 4, 1, 1, ia )
         endif

         Do im1=1, immax
            Do im3=1, immax
               Do is1=1, ndim_spinor
                  Do is3=1, ndim_spinor
                     istmp = ( is1 -1 )*ndim_spinor +is3
                     matH(immax*(is1-1)+im1,immax*(is3-1)+im3) &
                          &  = matH(immax*(is1-1)+im1,immax*(is3-1)+im3) &
                          &   +gzi *Mat_LS_with_real_ylm_L3(im1,im3,istmp)
                  End Do
               End Do
            End Do
         End Do
      endif

100   continue
!
      lwork = 2 *size2
      allocate( work(lwork) )
      allocate( eigenval(size2) ); eigenval = 0.0d0
      allocate( rwork( 3*size2 ) );
      call zheev( 'V', 'U', size2, MatH, size2, &
           &       eigenval, work, lwork, rwork, info )

      if ( mype == 0 ) then
         write(nfout,*) "** EigenValues of Hamiltonian of atom ", ia
         write(nfout,*) "  no.  value (cm-1)"
         Do im1=1, size2
!            write(nfout,*) im1, ( eigenval(im1) ) *Hartree_2_cminv
            write(nfout,*) im1, ( eigenval(im1) -eigenval(1) ) *Hartree_2_cminv
         End Do
      endif

      deallocate( work ); deallocate( eigenval );  deallocate( rwork )
      deallocate( matH )
    end subroutine calc_hamil_eigenval

  end subroutine m_CF_calc_CF_param

  subroutine m_CF_set_gaunt_coeff
    integer :: l1, l2, l3, m1, m2, m3, im1, im2, im3
    integer :: nsph1, nsph2, nsph3
    real(kind=DP) :: c0, c1, c2, c3
    complex(kind=CMPLDP) :: zsum
    complex(kind=CMPLDP), allocatable :: matwk(:,:,:)
    complex(kind=CMPLDP), pointer :: matU1(:,:), matU2(:,:), matU3(:,:)

    integer :: lcmax, n
    real(kind=DP), allocatable, dimension(:,:,:) :: cr2
    integer, allocatable, dimension(:,:,:)       :: isph2
    integer, allocatable, dimension(:,:)         :: mmt2
    real(kind=DP), allocatable :: matwk2(:,:,:)

    if ( allocated(gaunt_coeff) ) deallocate( gaunt_coeff )
    allocate( gaunt_coeff(16,49,16) ); gaunt_coeff = 0.0d0

    allocate( matwk(-3:3,-6:6,-3:3) ); matwk = 0.0d0

    Do l1=0, 3
       select case(l1)
       case (0)
          matU1 =>  MatU_ylm_RC_L0
       case (1)
          matU1 =>  MatU_ylm_RC_L1
       case (2)
          matU1 =>  MatU_ylm_RC_L2
       case (3)
          matU1 =>  MatU_ylm_RC_L3
       end select

       Do l2=0, 6
          select case(l2)
          case (0)
             matU2 =>  MatU_ylm_RC_L0
          case (1)
             matU2 =>  MatU_ylm_RC_L1
          case (2)
             matU2 =>  MatU_ylm_RC_L2
          case (3)
             matU2 =>  MatU_ylm_RC_L3
          case (4)
             matU2 =>  MatU_ylm_RC_L4
          case (5)
             matU2 =>  MatU_ylm_RC_L5
          case (6)
             matU2 =>  MatU_ylm_RC_L6
          end select

          Do l3=0, 3
             select case(l3)
             case (0)
                matU3 =>  MatU_ylm_RC_L0
             case (1)
                matU3 =>  MatU_ylm_RC_L1
             case (2)
                matU3 =>  MatU_ylm_RC_L2
             case (3)
                matU3 =>  MatU_ylm_RC_L3
             end select

             c0 = ( 2*l1+1 )*( 2*l2 +1 )*( 2*l3 +1 ) /PAI4
             c0 = sqrt( c0 )

             matwk = 0.0d0
             Do m1=-l1, l1
                Do m2=-l2, l2
                   Do m3=-l3, l3
                      c1 = wigner_3j( dble(l1), 0.0d0, dble(l2), 0.0d0, dble(l3), 0.0d0 )
                      c2 = wigner_3j( dble(l1), -dble(m1), dble(l2), dble(m2), &
                           &          dble(l3), dble(m3) )
                      c3 = (-1)**m1
                      matwk(m1,m2,m3) = c0 *c1 *c2 *c3
                   End Do
                End Do
             End Do
             !
             Do im1=1, 2*l1 +1
                Do im2=1, 2*l2 +1
                   Do im3=1, 2*l3 +1
                      nsph1 = l1**2 +im1
                      nsph2 = l2**2 +im2
                      nsph3 = l3**2 +im3
                      zsum = 0.0d0
                      Do m1=-l1, l1
                         Do m2=-l2, l2
                            Do m3=-l3, l3
                               zsum = zsum +conjg(matU1(im1,m1)) *matU2(im2,m2) &
                                    &      *matU3(im3,m3) *matwk(m1,m2,m3)
                            End Do
                         End Do
                      ENd Do
                      gaunt_coeff(nsph1,nsph2,nsph3) = zsum
                   End Do
                ENd Do
             ENd Do
          End Do
       ENd Do
    End Do
    deallocate( matwk )
  end subroutine m_CF_set_gaunt_coeff

  real(kind=8) function wigner_3j( val_j1, val_m1, val_j2, val_m2, val_j3, val_m3 )
    real(kind=8), intent(in) :: val_j1, val_j2, val_j3, val_m1, val_m2, val_m3

    integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, k, kmin, kmax
    real(kind=8) :: c1, c2, c3, c4, c5, c6, c7, c8, c9, c10
    real(kind=8) :: d1, d2, term_n, term_s, coeff

    if ( nint(val_m1 +val_m2 +val_m3) .ne. 0 ) then
       wigner_3j = 0.0d0;  return
    endif

    coeff = (-1)**( val_j1 -val_j2 -val_m3 )

    i1 = nint( val_j3 +val_j1 -val_j2 );    i2 = nint( val_j3 -val_j1 +val_j2 )
    i3 = nint( val_j1 +val_j2 -val_j3 );    i4 = nint( val_j3 -val_m3 )
    i5 = nint( val_j3 +val_m3 )
    i6 = nint( val_j1 +val_j2 +val_j3 +1 )
    i7 = nint( val_j1 -val_m1 );            i8 = nint( val_j1 +val_m1 )
    i9 = nint( val_j2 -val_m2 );            i10 = nint( val_j2 +val_m2 )

    c1 = factorial(i1);      c2 = factorial(i2);       c3 = factorial(i3)
    c4 = factorial(i4);      c5 = factorial(i5);       c6 = factorial(i6)
    c7 = factorial(i7);      c8 = factorial(i8);       c9 = factorial(i9)
    c10= factorial(i10)

    d1 = c1 *c2 *c3 *c4 *c5
    d2 = c6 *c7 *c8 *c9 *c10
    write(*,*) "d1 = ", d1, d2

    term_n = sqrt( d1 /d2 )

    kmin = max( 0, nint(val_m1 -val_j1) )
    kmin = max( kmin, nint(val_j2 -val_j1 -val_m3) )
    kmax = min( nint(val_j2 +val_j3 +val_m1), nint(val_j3 -val_j1 +val_j2) )
    kmax = min( kmax, nint(val_j3 -val_m3) )

    term_s = 0.0d0
    Do k=kmin, kmax
!    Do k=0, 100
       i2 = nint( val_j2 +val_j3 +val_m1 -k )
       i3 = nint( val_j1 -val_m1 +k )
       i5 = nint( val_j3 -val_j1 +val_j2 -k )
       i6 = nint( val_j3 -val_m3 -k )
       i7 = nint( val_j1 -val_j2 +val_m3 +k )
!       if ( i2 < 0 ) cycle
!       if ( i3 < 0 ) cycle
!       if ( i5 < 0 ) cycle
!       if ( i6 < 0 ) cycle
!       if ( i7 < 0 ) cycle

       c1 = (-1)**( val_j2 +val_m2 +k )
       c2 = factorial(i2);       c3 = factorial(i3);
       c4 = factorial(k)
       c5 = factorial(i5);       c6 = factorial(i6);       c7 = factorial(i7)
       d1 = c1 *c2 *c3
       d2 = c4 *c5 *c6 *c7
       term_s = term_s +d1 /d2
    End Do
    wigner_3j = coeff *term_n *term_s
    return
  end function wigner_3j

  real(kind=8) recursive function factorial(n)
    integer, intent(in) :: n

    integer :: i

    factorial = 1
    Do i=1, n
       factorial = factorial *i
    End do
    return
  end function factorial

  integer function get_m_tesseral_harmonics( il_in, im_in )
    integer, intent(in) :: il_in, im_in
    integer :: im_out

    select case(il_in)
    case (0)
       im_out = 0
    case (1)
       if ( im_in == 1 ) im_out =  1
       if ( im_in == 2 ) im_out = -1
       if ( im_in == 3 ) im_out =  0
    case (2)
       if ( im_in == 1 ) im_out =  0
       if ( im_in == 2 ) im_out =  2
       if ( im_in == 3 ) im_out = -2
       if ( im_in == 4 ) im_out = -1
       if ( im_in == 5 ) im_out =  1
    case (3,4,5,6)
       im_out = int( im_in /2 )
       if ( mod(im_in,2) == 1 ) im_out = -im_out
    end select
    get_m_tesseral_harmonics = im_out
  end function get_m_tesseral_harmonics


  subroutine calc_vh_asym( nmesh, rpos, rho, vh, qnum_l )
    integer, intent(in) :: nmesh, qnum_l
    real(kind=DP), intent(in) :: rpos(nmesh), rho(nmesh)
    real(kind=DP), intent(out) :: vh(nmesh)

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
                sum1 = sum1 +rpos(i0+j*is) **(qnum_l +2) *rho(i0+j*is) *wt(i0+j*is)
             end do
          end do
       else
          call set_weight_exp(ier,1,ir,rpos,wt)
          do jr = 1,ir
             sum1 = sum1 + rpos(jr) **( qnum_l +2 ) *rho(jr)*wt(jr)
          end do
       end if
       sum1 = sum1 *PAI4 /rpos(ir) **( qnum_l +1 )

       if (ir == nmesh) then
          sum2 = 0.d0
       else if ((ir <= nmesh-1).and.(ir >= nmesh-4)) then
          do ii = ir,nmesh-1
             i0 = ii+1; is = -1
             call set_open_weight_exp(ier,i0,is,rpos,wt)
             do j = 1,4
                sum2 = sum2 -rpos(i0+j*is)**(qnum_l +2) *rho(i0+j*is) *wt(i0+j*is)
             end do
          end do
       else
          call set_weight_exp(ier,ir,nmesh,rpos,wt)
          do jr = ir,nmesh
             sum2 = sum2 + rpos(jr)**( 1-qnum_l ) *rho(jr)*wt(jr)
          end do
       end if
       sum2 = sum2 *PAI4 *rpos(ir)**( qnum_l )
       vh(ir) = ( sum1 + sum2 ) / ( 2.*qnum_l +1. )
    end do

    deallocate( wt )
  end subroutine calc_vh_asym

  subroutine m_CF_calc_mutiple_moments
    integer :: lcmax
    integer :: ia, it, ig, il1, im1, nsph1, ist, is
    real(kind=DP) :: c1, fac, rcut
    complex(kind=CMPLDP) :: z1, z2, zsum, zsum_mpi
    real(kind=DP), allocatable :: wka(:), wkb(:), qx(:), qy(:), qz(:), ylm(:)
    real(kind=DP), allocatable :: zfcos(:), zfsin(:)
    real(kind=DP), allocatable :: multipole(:)
!
    lcmax = 6
    allocate( wka(ista_kngp:iend_kngp) )
    allocate( wkb(ista_kngp:iend_kngp) )
    allocate( ylm(ista_kngp:iend_kngp) )
    allocate( qx(ista_kngp:iend_kngp) )
    allocate( qy(ista_kngp:iend_kngp) )
    allocate( qz(ista_kngp:iend_kngp) )
    allocate( zfcos(ista_kngp:iend_kngp) )
    allocate( zfsin(ista_kngp:iend_kngp) )

    Do ig=ista_kngp, iend_kngp
       qx(ig) = rltv(1,1)*ngabc(ig,1) +rltv(1,2)*ngabc(ig,2) +rltv(1,3)*ngabc(ig,3)
       qy(ig) = rltv(2,1)*ngabc(ig,1) +rltv(2,2)*ngabc(ig,2) +rltv(2,3)*ngabc(ig,3)
       qz(ig) = rltv(3,1)*ngabc(ig,1) +rltv(3,2)*ngabc(ig,2) +rltv(3,3)*ngabc(ig,3)
    End do

    allocate( multipole(49))

    Do ia=1, natm
       it = ityp(ia)
       if ( lpsmax(it) /=4 ) cycle

!       rcut = rad_cov_default( iatomn(it) )
       rcut = rad_cov( ia )
       fac = rcut**3 *PAI4

       call calc_phase2( natm, pos, ia, kgp, ngabc, ista_kngp, iend_kngp, &
            &            zfcos, zfsin )
       Do ig=ista_kngp, iend_kngp
          wka(ig) = gr_l(ig) *rcut
       End do

       multipole = 0.0d0

       Do il1=0, lcmax
          call dsjnv( il1 +1, iend_kngp -ista_kngp +1, &
               &      wka(ista_kngp:iend_kngp), &
               &      wkb(ista_kngp:iend_kngp) )

          Do im1=1, 2*il1 +1
             nsph1 = il1**2 +im1
             call sphr_general( iend_kngp -ista_kngp +1, nsph1, &
                  &     qx(ista_kngp:iend_kngp), &
                  &     qy(ista_kngp:iend_kngp), &
                  &     qz(ista_kngp:iend_kngp), ylm(ista_kngp:iend_kngp) )

             ist = ista_kngp
             if ( ist == 1 ) ist = 2

             Do is=1, nspin
                zsum = 0.0d0
                Do ig=ist, iend_kngp
                   c1 = ylm(ig) *wkb(ig) /wka(ig)
                   z1 = cmplx( chgq_l(ig,1,is) , chgq_l(ig,2,is) )
                   z2 = cmplx( zfcos(ig), zfsin(ig) )
                   zsum = zsum +z1 *c1 *z2
                End Do
                zsum = zsum *zi **il1 *fac
                if ( il1 == 0 ) then
                   if ( mype == 0 ) zsum = zsum +chgq_l(1,1,is) *fac /3.d0 *ylm(1)
                endif
                multipole(nsph1) = multipole(nsph1) +zsum
             End Do
          end Do
       end Do
       call mpi_allreduce( MPI_IN_PLACE, multipole, 49, mpi_double_precision, &
            &              mpi_sum, MPI_CommGroup,ierr )
!
       Do nsph1=1, 49
          write(nfout,'(A,2I5,F20.12)') &
               &    "ia nsph1 contrib.", ia, nsph1, multipole(nsph1)
       end Do
    End Do
    deallocate( qx );    deallocate( qy );    deallocate( qz )
    deallocate( ylm );   deallocate( wka );   deallocate( wkb )
    deallocate( zfcos );    deallocate( zfsin )
    deallocate( multipole )

  end subroutine m_CF_calc_mutiple_moments

#if 0
R: 4f atom no aru site
vh_lm(|r-R|) = int vh(r) Ylm(|r-R|)
             = int sum_G vh(G)exp(iGr) Ylm(|r-R|)      ( r-R=r'; r = r'+R )
             =                exp(iGr') exp(iGR) Ylm(|r'|)
             = sum_G vh(G) exp(iGR) Ylm(G) j_l(Gr') ??
vh(r) R4f(r-R)^2 (r-R)^2 dr
#endif


end module m_Crystal_Field
