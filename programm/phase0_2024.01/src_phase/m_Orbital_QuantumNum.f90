!
! =========== Contributions ===================================
!
! Through the courtesy of contributors, this module is added.
!
! Company:  ASMS Co.,Ltd.
! Functions:
!
!     Conversion table between (n,l,t) in CIAO and (l,tau) in PHASE
!
! =============================================================

module m_Orbital_QuantumNum

  use m_Control_Parameters,  only : printable
  use m_Files,  only : nfout, nfpot
  use m_Ionic_System,   only : ntyp
  use m_Parallelization,  only : mype, MPI_CommGroup, ierr
  use mpi

  implicit none
!  include 'mpif.h'

  integer :: max_num_orb_index
  integer, allocatable :: num_orb_index_data(:)
  integer, allocatable :: qnum_n_orb_index(:,:), qnum_l_orb_index(:,:), &
       &                  qnum_t_orb_index(:,:), qnum_tau_orb_index(:,:)

contains

  subroutine m_OP_Qnum_init
    if ( allocated( num_orb_index_data ) )  deallocate( num_orb_index_data )
    allocate( num_orb_index_data(ntyp) ) ; num_orb_index_data = 0
  end subroutine m_OP_Qnum_init

  subroutine m_OP_Qnum_dealloc_array
    if ( allocated( num_orb_index_data ) ) deallocate( num_orb_index_data )
    if ( allocated( qnum_n_orb_index ) )   deallocate( qnum_n_orb_index )
    if ( allocated( qnum_l_orb_index ) )   deallocate( qnum_l_orb_index )
    if ( allocated( qnum_t_orb_index ) )   deallocate( qnum_t_orb_index )
    if ( allocated( qnum_tau_orb_index ) ) deallocate( qnum_tau_orb_index )
  end subroutine m_OP_Qnum_dealloc_array

  subroutine m_OP_Qnum_read_orbital_index
    call find_orbtital_index_data_number
    call alloc_orbital_index_data_array
    call read_orbital_index_data

  contains

    subroutine find_orbtital_index_data_number
      integer :: nfpp, it

      nfpp=0
      if(mype == 0) then
         write(nfout,*) "!** [asms extension] --- orbital index info ---"
         do it=1, ntyp
            nfpp=nfpp+1
            call find_orbital_index_section( nfpot(nfpp), nfout, it, &
                 &                           num_orb_index_data(it) )
            if (printable) then
               write(nfout,'(A,i3,A,i3)') ' !** num of orbital index. for it = ', &
                    &                     it, '  is ', num_orb_index_data(it)
            endif
         end do
      end if
      call mpi_bcast( num_orb_index_data, ntyp, mpi_integer, 0, MPI_CommGroup, ierr )

    end subroutine find_orbtital_index_data_number

    subroutine find_orbital_index_section( nfp, nfout, it, num_orb_index )
      integer, intent(in) :: nfp, nfout, it
      integer, intent(out) :: num_orb_index

      integer :: ikey
      character(len=30)    :: line

      if(mype /= 0) return
      num_orb_index = 0

      rewind(nfp)
      Do while (.true.)
         read(nfp,'(a30)',end=20) line
         ikey = index(line,'ORBITAL-INDEX')
         if ( ikey /=0 ) then
            read(nfp,*) num_orb_index
            return
         endif
      End Do
20    if (printable) then
         write(nfout,*) "!* orbital index is not found in PP for it=", it
      endif

    end subroutine find_orbital_index_section

    subroutine alloc_orbital_index_data_array
      if (mype == 0) then
         max_num_orb_index = maxval( num_orb_index_data )
      endif
      call mpi_bcast( max_num_orb_index, 1,  mpi_integer, 0, MPI_CommGroup, ierr )

      if (printable) then
         write(nfout,'(A,i3)') ' !** max_num_orb_index = ', max_num_orb_index
      endif

      if ( max_num_orb_index == 0 ) return

      allocate( qnum_n_orb_index(ntyp,max_num_orb_index) );   qnum_n_orb_index = 0
      allocate( qnum_l_orb_index(ntyp,max_num_orb_index) );   qnum_l_orb_index = 0
      allocate( qnum_t_orb_index(ntyp,max_num_orb_index) );   qnum_t_orb_index = 0
      allocate( qnum_tau_orb_index(ntyp,max_num_orb_index) ); qnum_tau_orb_index = 0

    end subroutine alloc_orbital_index_data_array

    subroutine read_orbital_index_data
      integer :: nfpp, it, j

      nfpp=0

      if ( max_num_orb_index == 0 ) return

      if (mype == 0)  then
         do it =1, ntyp
            nfpp=nfpp+1
            Do j=1, num_orb_index_data(it)
               read( nfpot(nfpp),*) qnum_n_orb_index(it,j), qnum_l_orb_index(it,j), &
                    &               qnum_t_orb_index(it,j), qnum_tau_orb_index(it,j)
            End do
         end do
      end if

      call mpi_bcast( qnum_n_orb_index, ntyp*max_num_orb_index, mpi_integer, 0, &
           &          MPI_CommGroup, ierr )
      call mpi_bcast( qnum_l_orb_index, ntyp*max_num_orb_index, mpi_integer, 0, &
           &          MPI_CommGroup, ierr )
      call mpi_bcast( qnum_t_orb_index, ntyp*max_num_orb_index, mpi_integer, 0, &
           &          MPI_CommGroup, ierr )
      call mpi_bcast( qnum_tau_orb_index, ntyp*max_num_orb_index, mpi_integer, 0, &
           &          MPI_CommGroup, ierr )

    end subroutine read_orbital_index_data

  end subroutine m_OP_Qnum_read_orbital_index

  subroutine m_OP_Qnum_orb_ind_data_num(it)
    integer, intent(in) :: it
    integer, save :: nfpp=0
    if(mype == 0) then
       write(nfout,*) "!** [asms extension] --- orbital index info ---"
       nfpp=nfpp+1
       call find_orbital_index_section( nfpot(nfpp), nfout, it, &
            &                           num_orb_index_data(it) )
       if (printable) then
          write(nfout,'(A,i3,A,i3)') ' !** num of orbital index. for it = ', &
               &                     it, '  is ', num_orb_index_data(it)
       endif
    end if

    if(it==ntyp) then
       call mpi_bcast( num_orb_index_data, ntyp, mpi_integer, 0, MPI_CommGroup, ierr )
       nfpp=0
    endif

    contains

    subroutine find_orbital_index_section( nfp, nfout, it, num_orb_index )
      integer, intent(in) :: nfp, nfout, it
      integer, intent(out) :: num_orb_index

      integer :: ikey
      character(len=30)    :: line

      if(mype /= 0) return
      num_orb_index = 0

      rewind(nfp)
      Do while (.true.)
         read(nfp,'(a30)',end=20) line
         ikey = index(line,'ORBITAL-INDEX')
         if ( ikey /=0 ) then
            read(nfp,*) num_orb_index
            return
         endif
      End Do
20    if (printable) then
         write(nfout,*) "!* orbital index is not found in PP for it=", it
      endif

    end subroutine find_orbital_index_section
  end subroutine m_OP_Qnum_orb_ind_data_num

  subroutine m_OP_Qnum_alloc_orb_ind_data()
    if (mype == 0) then
       max_num_orb_index = maxval( num_orb_index_data )
    endif
    call mpi_bcast( max_num_orb_index, 1,  mpi_integer, 0, MPI_CommGroup, ierr )

    if (printable) then
       write(nfout,'(A,i3)') ' !** max_num_orb_index = ', max_num_orb_index
    endif

    if ( max_num_orb_index == 0 ) return
    if(allocated(qnum_n_orb_index)) deallocate(qnum_n_orb_index)
    if(allocated(qnum_l_orb_index)) deallocate(qnum_l_orb_index)
    if(allocated(qnum_t_orb_index)) deallocate(qnum_t_orb_index)
    if(allocated(qnum_tau_orb_index)) deallocate(qnum_tau_orb_index)
    allocate( qnum_n_orb_index(ntyp,max_num_orb_index) );   qnum_n_orb_index = 0
    allocate( qnum_l_orb_index(ntyp,max_num_orb_index) );   qnum_l_orb_index = 0
    allocate( qnum_t_orb_index(ntyp,max_num_orb_index) );   qnum_t_orb_index = 0
    allocate( qnum_tau_orb_index(ntyp,max_num_orb_index) ); qnum_tau_orb_index = 0
  end subroutine m_OP_Qnum_alloc_orb_ind_data

  subroutine m_OP_Qnum_read_orbital_index_it(it)
    integer, intent(in) :: it
    integer,save :: nfpp=0
    integer ::  j

    if ( max_num_orb_index == 0 ) return

    if (mype == 0)  then
       nfpp=nfpp+1
       Do j=1, num_orb_index_data(it)
          read( nfpot(nfpp),*) qnum_n_orb_index(it,j), qnum_l_orb_index(it,j), &
               &               qnum_t_orb_index(it,j), qnum_tau_orb_index(it,j)
       End do
    end if
    if(it == ntyp) then
       call mpi_bcast( qnum_n_orb_index, ntyp*max_num_orb_index, mpi_integer, 0, &
            &          MPI_CommGroup, ierr )
       call mpi_bcast( qnum_l_orb_index, ntyp*max_num_orb_index, mpi_integer, 0, &
            &          MPI_CommGroup, ierr )
       call mpi_bcast( qnum_t_orb_index, ntyp*max_num_orb_index, mpi_integer, 0, &
            &          MPI_CommGroup, ierr )
       call mpi_bcast( qnum_tau_orb_index, ntyp*max_num_orb_index, mpi_integer, 0, &
            &          MPI_CommGroup, ierr )
       nfpp = 0
    endif
  end subroutine m_OP_Qnum_read_orbital_index_it

end module m_Orbital_QuantumNum

