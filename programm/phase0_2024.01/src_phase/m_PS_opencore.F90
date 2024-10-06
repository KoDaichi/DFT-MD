module m_PS_opencore
  use m_Files,  only : nfpot, nfout, m_Files_open_ps_files, m_Files_close_ps_file, &
       &               m_Files_open_mag_opencore, m_Files_close_mag_opencore, &
       &               nf_mag_core
  use m_PseudoPotential,  only : nmesh, mmesh, radr_paw, radr, wos
  use m_NonLocal_Potential, only : new_radr_and_wos

  use m_PlaneWaveBasisSet,  only : gr_l
  use m_Parallelization,  only : npes, mype, ista_kngp, iend_kngp, ista_k, MPI_CommGroup
  use m_Const_Parameters,  only : DP, ON
  use m_Crystal_Structure,  only : univol
  use m_Ionic_System,  only : ivan, ntyp, iatomn, natm, ityp
  use m_Control_Parameters,  only : noncol,  core_spin_pol_factor, sw_fix_core_spin_pol
  use m_CD_Mag_Moment,   only :  RhoMag_on_atom

  use m_Control_Parameters,  only : nspin, sw_opencore, iprimagmom
  use m_Const_Parameters,  only : OFF
  use m_Ionic_System,  only : mag_moment0_atoms_is_defined, &
       &                      mag_moment0_atomtyp, mag_moment0_atoms, speciesname
  use mpi

  implicit none
!  include 'mpif.h'
!
  integer, allocatable :: has_opencore(:)
  real(kind=DP), allocatable :: mag_opencore_pol(:,:)
  real(kind=DP), allocatable :: rmag_opencore(:,:)
  real(kind=DP), allocatable :: rmag_opencore_l(:,:)
  real(kind=DP), allocatable :: m_opencore(:)

contains

  subroutine m_PS_set_rmag_opencore
    integer :: it, ierr, nfp, i, ir
    integer :: num_core_ae_wfns, my_l
    integer, allocatable :: qnum_n_core_ae_wfns(:)
    integer, allocatable :: qnum_l_core_ae_wfns(:)
    real(kind=DP), allocatable :: psir_core_ae_wfns(:,:)
    real(kind=DP), allocatable :: enelevel_core_ae_wfns(:)
    real(kind=DP), allocatable :: focc_core_ae_wfns(:)
    real(kind=DP) :: weight

! == KT_DEBUG === 2015/06/15 ==
!    call m_Files_open_ps_file(ivan,iatomn,ntyp,it,ierr)
!                                      ! this does not work properly by unknown reason.
    call m_Files_open_ps_files(ivan,iatomn,ntyp,ierr)
! =============== 2015/06/15

    if (ierr/=0) call mpi_stop(nfout)

    if ( allocated(rmag_opencore) ) deallocate(rmag_opencore)
    if ( allocated(has_opencore) )  deallocate(has_opencore)
    if ( allocated(m_opencore) )    deallocate(m_opencore)
    allocate( rmag_opencore(mmesh,ntyp) );   rmag_opencore = 0.0d0
    allocate( has_opencore(ntyp) );   has_opencore = 0
    allocate( m_opencore(ntyp) );    m_opencore = 0.0d0

    Do it=1, ntyp
       nfp = nfpot(it)

       if ( mype == 0 ) call read_num_core_ae_wfns( nfp, num_core_ae_wfns, it )
       if ( npes > 1 ) then
          call mpi_bcast(  num_core_ae_wfns, 1, mpi_integer, 0, mpi_comm_world, ierr )
       endif

       if ( num_core_ae_wfns == 0 ) cycle

       allocate( qnum_n_core_ae_wfns( num_core_ae_wfns ) )
       allocate( qnum_l_core_ae_wfns( num_core_ae_wfns ) )
       allocate( psir_core_ae_wfns( nmesh(it),  num_core_ae_wfns ) )
       allocate( enelevel_core_ae_wfns( num_core_ae_wfns ) )
       allocate( focc_core_ae_wfns( num_core_ae_wfns ) )

       if ( mype == 0 ) then
          call read_data_core_ae_wfns( nfp, num_core_ae_wfns, nmesh(it), &
               &                       qnum_n_core_ae_wfns, qnum_l_core_ae_wfns, &
               &                       psir_core_ae_wfns, &
               &                       enelevel_core_ae_wfns, focc_core_ae_wfns )
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
       endif
!
       Do i=1, num_core_ae_wfns
          my_l = qnum_l_core_ae_wfns(i)
          if ( my_l <= 2 ) then            ! for core-hole
             if ( abs( focc_core_ae_wfns(i) -2.0*(2*my_l +1) ) < 0.1 ) cycle

             has_opencore(it) = 1
             if ( focc_core_ae_wfns(i) <= 2*my_l +1  ) then
                weight = focc_core_ae_wfns(i)
             else
                weight =  2*(2*my_l +1) -focc_core_ae_wfns(i)
             endif
             Do ir=1, nmesh(it)
                rmag_opencore(ir,it) = weight *psir_core_ae_wfns(ir,i) **2
             End Do
             m_opencore(it) = weight
          endif
          if ( qnum_l_core_ae_wfns(i) == 3 ) then
             if ( abs( focc_core_ae_wfns(i) -14.0 ) < 0.1 ) cycle

             has_opencore(it) = 1
             if ( focc_core_ae_wfns(i) <=7  ) then
                weight = focc_core_ae_wfns(i)
             else
                weight =  14.d0 -focc_core_ae_wfns(i)
             endif
             Do ir=1, nmesh(it)
                rmag_opencore(ir,it) = weight *psir_core_ae_wfns(ir,i) **2
             End Do
             m_opencore(it) = weight
          endif
       End Do

       deallocate( qnum_n_core_ae_wfns )
       deallocate( qnum_l_core_ae_wfns )
       deallocate( psir_core_ae_wfns )
       deallocate( enelevel_core_ae_wfns )
       deallocate( focc_core_ae_wfns )

       call m_Files_close_ps_file(it)
    End Do
  end subroutine m_PS_set_rmag_opencore

  subroutine m_PS_set_pcc_opencore
    integer :: ierr, i, ir, n, it, ik
    real(kind=DP) :: gabs, fac
    real(kind=DP), allocatable :: wkx(:), wky(:)

    allocate( rmag_opencore_l( ista_kngp:iend_kngp, ntyp ) )
    rmag_opencore_l = 0.0d0

    allocate(wkx(mmesh)); wkx = 0.d0
    allocate(wky(mmesh)); wky = 0.d0
    allocate(wos(mmesh)); wos = 0

#if 1
    allocate( radr(mmesh)); radr = 0.0d0
#endif


    Do it=1, ntyp
#if 0
       call set_weight_exp( ierr, 1, nmesh(it), radr_paw(:,it), wos )
#else
       ik = ista_k
       call new_radr_and_wos(ik,it)
!       write(*,*) "wos", wos(1:5)
#endif

       if ( has_opencore(it) == 0 ) cycle
!       write(*,*) 'FFFF ', it
!       write(*,*) "FF2 ",  rmag_opencore(10,it)

       do i = ista_kngp, iend_kngp  !for mpi
          gabs = gr_l(i)
#if 0
          do n = 1, nmesh(it)
             wkx(n) = gabs * radr_paw(n,it)
          enddo
#else
          do n = 1, nmesh(it)
             wkx(n) = gabs * radr(n)
          enddo
#endif
          call dsjnv(0,nmesh(it),wkx, wky)

          do  n = 1, nmesh(it)
             fac = wos(n) *rmag_opencore(n,it) /univol
             rmag_opencore_l(i,it) = rmag_opencore_l(i,it) + fac *wky(n)
          enddo
       end do
    End Do
    deallocate( wos )
    deallocate( wkx )
    deallocate( wky )
#if 1
    deallocate( radr )
#endif

  end subroutine m_PS_set_pcc_opencore

  subroutine m_PS_init_mag_opencore_pol
    integer :: ia, it
    real(kind=DP) :: cnorm, fac(3)

    if ( sw_opencore == OFF ) return
    if ( nspin == 1 ) return

    call m_PS_alloc_mag_opencore_pol

    Do ia=1, natm
       it = ityp(ia)
       if ( has_opencore(it) == 0 ) cycle

       if ( noncol ) then
          fac(1) = 0.0d0;      fac(2) = 0.0d0;     fac(3) = 1.0d0

          if ( mag_moment0_atoms_is_defined ) then
             cnorm = sqrt( mag_moment0_atoms(ia,1)**2 &
                  &      +mag_moment0_atoms(ia,2)**2 &
                  &      +mag_moment0_atoms(ia,3)**2 )
             if ( cnorm > 0.0d0 ) then
                fac(1) = mag_moment0_atoms(ia,1) /cnorm
                fac(2) = mag_moment0_atoms(ia,2) /cnorm
                fac(3) = mag_moment0_atoms(ia,3) /cnorm
             endif
          else
             cnorm = sqrt( mag_moment0_atomtyp(it,1)**2 &
                  &      +mag_moment0_atomtyp(it,2)**2 &
                  &      +mag_moment0_atomtyp(it,3)**2 )
             if ( cnorm > 0.0d0 ) then
                fac(1) = mag_moment0_atomtyp(it,1) /cnorm
                fac(2) = mag_moment0_atomtyp(it,2) /cnorm
                fac(3) = mag_moment0_atomtyp(it,3) /cnorm
             endif
          endif
          mag_opencore_pol(ia,1:3) = core_spin_pol_factor *fac(1:3)
       else
          fac(1) = 1.0d0
          if ( mag_moment0_atoms_is_defined ) then
             if ( mag_moment0_atoms(ia,1) < 0.0 ) fac(1) = -1.0d0
          else
             if ( mag_moment0_atomtyp(it,1) < 0.0 ) fac(1) = -1.0d0
          endif
          mag_opencore_pol(ia,1) = core_spin_pol_factor *fac(1)
       endif
    End do
    return
#if 0
    Do ia=1, natm
       it = ityp(ia)
       if ( has_opencore(it) == 0 ) cycle
       if ( noncol ) then
          mag_opencore_pol(ia,3) = core_spin_pol_factor
       else
          mag_opencore_pol(ia,1) = core_spin_pol_factor
       endif
    End do
#endif
  end subroutine m_PS_init_mag_opencore_pol

  subroutine m_PS_set_mag_opencore_pol
    integer :: ia, it
    real(kind=DP) :: c1, c2

    if ( sw_opencore == OFF ) return
    if ( nspin == 1 ) return
    if ( sw_fix_core_spin_pol == ON ) return

    Do ia=1, natm
       it = ityp(ia)
       if ( has_opencore(it) == 0 ) cycle

       if ( noncol ) then
          c1 = RhoMag_on_atom(ia,2)**2 +RhoMag_on_atom(ia,3)**2 +RhoMag_on_atom(ia,4)**2
          c2 = sqrt(c1)
          if ( c2 >1.0D-11 ) then
             mag_opencore_pol(ia,1) = RhoMag_on_atom(ia,2) /c2 *core_spin_pol_factor
             mag_opencore_pol(ia,2) = RhoMag_on_atom(ia,3) /c2 *core_spin_pol_factor
             mag_opencore_pol(ia,3) = RhoMag_on_atom(ia,4) /c2 *core_spin_pol_factor
          else
             mag_opencore_pol(ia,3) = core_spin_pol_factor
          endif
       else
          c1 = RhoMag_on_atom(ia,1) - RhoMag_on_atom(ia,2)
          c2 = abs(c1)
          if ( c2 >1.0D-11 ) then
             mag_opencore_pol(ia,1) = core_spin_pol_factor *c1 /c2
          else
             mag_opencore_pol(ia,1) = core_spin_pol_factor
          endif
       endif
    End Do
  end subroutine m_PS_set_mag_opencore_pol

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

  subroutine m_PS_print_core_Magmom_on_atom( nfout )
    integer, intent(in) :: nfout

    integer :: ia, it
    real(kind=DP) :: zz

    if ( iprimagmom < 2 ) return
    if ( mype /= 0 ) return
    if ( sw_opencore == OFF ) return

    write(nfout,*) &
         & '! ------------ Core Charge Spin Moment at this scf step --- '

    if ( noncol ) then
       write(nfout,*) &
            & '    id      name          mx             my             mz'
       Do ia=1, natm
          it = ityp(ia)
          if ( has_opencore(it) == 0 ) cycle
          write(nfout,'(I7,6X,A5,3F15.8)') ia, speciesname(it), &
               &                           m_opencore(it) *mag_opencore_pol(ia,1:3)
       End Do

    else if ( nspin == 2 ) then
       write(nfout,*) &
            & '    id      name          mz'
       Do ia=1, natm
          it = ityp(ia)
          if ( has_opencore(it) == 0 ) cycle
          write(nfout,'(I7,6X,A5,F15.8)') ia, speciesname(it), &
               &                          m_opencore(it) *mag_opencore_pol(ia,1)
       End Do
    endif
    write(nfout,*)

  end subroutine m_PS_print_core_Magmom_on_atom

  subroutine m_PS_alloc_mag_opencore_pol
    if ( noncol ) then
       allocate( mag_opencore_pol(natm,3) )
    else
       allocate( mag_opencore_pol(natm,1) )
    endif
    mag_opencore_pol = 0.0d0
  end subroutine m_PS_alloc_mag_opencore_pol

  subroutine m_PS_dealloc_mag_opencore_pol
    deallocate( mag_opencore_pol )
  end subroutine m_PS_dealloc_mag_opencore_pol

  subroutine m_PS_wd_mag_opencore_pol
    integer :: ia, it

    if ( mype /= 0 ) return
    if ( sw_opencore == OFF ) return
    if ( nspin == 1 ) return

    call m_Files_open_mag_opencore
    if ( noncol ) then
       Do ia=1, natm
          it = ityp(ia)
          if ( has_opencore(it) == 0 ) cycle
          write(nf_mag_core,'(I7,6X,3F20.12)') ia, mag_opencore_pol(ia,1:3)
       End Do
    else
       Do ia=1, natm
          it = ityp(ia)
          if ( has_opencore(it) == 0 ) cycle
          write(nf_mag_core,'(I7,6X,F20.12)') ia, mag_opencore_pol(ia,1)
       End Do
    endif
    call m_Files_close_mag_opencore

  end subroutine m_PS_wd_mag_opencore_pol

  subroutine m_PS_rd_mag_opencore_pol
    integer :: ia, it, num, ierr, ipol_max
    if ( sw_opencore == OFF ) return
    if ( nspin == 1 ) return

    call m_PS_alloc_mag_opencore_pol
    if ( noncol ) then
       ipol_max = 3
    else
       ipol_max = 1
    endif

    if ( mype == 0 ) then
       call m_Files_open_mag_opencore
       Do ia=1, natm
          it = ityp(ia)
          if ( has_opencore(it) == 0 ) cycle
          read(nf_mag_core,*) num, mag_opencore_pol(ia,1:ipol_max)
       End Do
       call m_Files_close_mag_opencore
    endif

    call mpi_bcast( mag_opencore_pol, natm*ipol_max, mpi_double_precision, &
         &          0, MPI_CommGroup, ierr )

  end subroutine m_PS_rd_mag_opencore_pol

end module m_PS_opencore
