!
! =========== Contributions ===================================
!
! Through the courtesy of contributors, this module is added.
!
! Company:  ASMS Co.,Ltd.
! Functions:
!
!     SpinOrbit strength (splitting) is read from GNCPP2 file
!
! =============================================================

module m_SpinOrbit_FromFile

  use m_Const_Parameters,  only : DP
  use m_Control_Parameters,  only: printable
  use m_Files, only : nfout, nfpot, m_Files_open_ps_files, m_Files_close_ps_files, &
       &              m_Files_open_ps_file, m_Files_close_ps_file
  use m_Ionic_System,  only : natm, scaling_so, ivan, iatomn
  use m_Parallelization,  only : mype, npes, MPI_CommGroup, ierr
  use m_PseudoPotential,  only : ntyp, ityp, Mat_SOC_Strength_nonpaw, nloc, ntau
  use m_Orbital_QuantumNum, only : num_orb_index_data, &
       &                           qnum_n_orb_index, qnum_l_orb_index, &
       &                           qnum_t_orb_index, qnum_tau_orb_index

  use m_SpinOrbit_Potential,  only : Mat_SOC_Strength
  use mpi

  implicit none
!  include 'mpif.h'

contains

  subroutine m_SO_set_SOC_strength_from_PP
    integer :: it, nfp, ia, ierror

    if (mype == 0) then
!       call m_Files_open_ps_files(ivan,iatomn,ntyp,ierror)
       nfp = 0
       Do it=1, ntyp
          call m_Files_open_ps_file(ivan,iatomn,ntyp,it,ierror)
          nfp = nfp +1
          call set_soc_from_pp_file( it, nfpot(nfp) )
          call m_Files_close_ps_file(it)
       End Do
    endif

    if ( npes > 2 ) then
       call mpi_bcast( Mat_SOC_Strength_nonpaw, nloc*ntau*ntau*ntyp, &
            &          mpi_double_precision, 0, MPI_CommGroup, ierr )
    endif

    Do ia=1, natm
       it = ityp(ia)
       Mat_SOC_Strength(:,:,:,ia) = Mat_SOC_Strength_nonpaw(:,:,:,it)  &
            &                       * scaling_so(it)
    End do

!    call m_Files_close_ps_files

  contains

    subroutine set_soc_from_pp_file( it, nfp )
      integer, intent(in) :: it, nfp

      integer :: ikey
      character(len=30)    :: line

      integer :: num_soc_read
      integer :: i, j1, n1, l1, t1, n2, l2, t2, tau1, tau2
      real(kind=DP) :: c1

      num_soc_read = 0

      rewind(nfp)
      Do while (.true.)
         read(nfp,'(a30)',end=20) line
         ikey = index(line,'SOC-VALENCE')
         if ( ikey /=0 ) then
            read(nfp,*) num_soc_read
            goto 10
         endif
      End Do

! ------------------
10    Do i=1, num_soc_read
         read(nfp,*) n1, l1, t1, n2, l2, t2, c1

         if ( l1 /= l2 ) cycle

         tau1 = 0;  tau2 = 0
         Do j1=1, num_orb_index_data(it)
            if ( n1 == qnum_n_orb_index(it,j1) &
                 &    .and. l1 == qnum_l_orb_index(it,j1) &
                 &    .and. t1 == qnum_t_orb_index(it,j1) ) then
               tau1 = qnum_tau_orb_index(it,j1)
            endif
         End do
         Do j1=1, num_orb_index_data(it)
            if ( n2 == qnum_n_orb_index(it,j1) &
                 &    .and. l2 == qnum_l_orb_index(it,j1) &
                 &    .and. t2 == qnum_t_orb_index(it,j1) ) then
               tau2 = qnum_tau_orb_index(it,j1)
            endif
         End do

         if ( tau1 ==0 .or. tau2 == 0 ) cycle

         Mat_SOC_Strength_nonpaw( l1 +1, tau1, tau2, it ) = c1
      End Do

      return

! ------------------
20    if (printable) then
         write(nfout,*) "!* soc strength for valence are not found in PP for it=", it
      endif

    end subroutine set_soc_from_pp_file

  end subroutine m_SO_set_SOC_strength_from_PP

end module m_SpinOrbit_FromFile
