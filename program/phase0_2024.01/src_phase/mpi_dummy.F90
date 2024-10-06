!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINES: mpi_init, mpi_finalize, mpi_comm_size, mpi_comm_rank,
!              mpi_barrier, mpi_send, mpi_recv, mpi_comm_dup, mpi_comm_split,
!              mpi_bcast, mpi_allgather, mpi_allgatherv, mpi_allreduce
!              mpi_copy, getarg, mpi_dummy
!
!  AUTHOR(S): H. Tsukioka   November/02/2003
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
!=======================================================================
!
!     The original version of this set of the computer programs "PHASE"
!  was developed by the members of the Theory Group of Joint Research
!  Center for Atom Technology (JRCAT), based in Tsukuba, in the period
!  1993-2001.
!
!     Since 2002, this set has been tuned and new functions have been
!  added to it as a part of the national project "Frontier Simulation
!  Software for Industrial Science (FSIS)",  which is supported by
!  the IT program of the Ministry of Education, Culture, Sports,
!  Science and Technology (MEXT) of Japan.
!     Since 2006, this program set has been developed as a part of the
!  national project "Revolutionary Simulation Software (RSS21)", which
!  is supported by the next-generation IT program of MEXT of Japan.
!   Since 2013, this program set has been further developed centering on PHASE System
!  Consortium.
!   The activity of development of this program set has been supervised by Takahisa Ohno.
!
!
! $Id: mpi_dummy.F90 570 2017-04-21 20:34:50Z yamasaki $
#ifdef _NO_MPI_

module mpi
        include 'mpif.h'
end module

subroutine mpi_init( ierr )
        integer :: ierr
        ierr = 0
end

subroutine mpi_finalize( ierr )
        integer :: ierr
        ierr = 0
end

subroutine mpi_abort(icomm, errorcode, ierr )
        integer :: icomm, errorcode, ierr
        errorcode=0
        ierr=0
end

subroutine mpi_comm_size( icomm, isize, ierr )
        integer :: icomm, isize, ierr
        isize = 1
        ierr  = 0
end

subroutine mpi_comm_rank( icomm, irank, ierr )
        integer :: icomm, irank, ierr
        irank = 0
        ierr  = 0
end

subroutine mpi_barrier( icomm, ierr )
        integer :: icomm, ierr
        ierr = 0
end

subroutine mpi_send( buf, icount, itype, idest, itag, icomm, ierr )
	integer :: icount, itype, idest, itag, icomm, ierr
	ierr = 0
end

subroutine mpi_isend( buf, icount, itype, idest, itag, icomm, irequest, ierr )
	integer :: icount, itype, idest, itag, icomm, irequest,ierr
	ierr = 0
end

subroutine mpi_recv( buf, icount, itype, isource, itag, icomm, istatus, ierr )
	integer :: icount, itype, isource, itag, icomm, istatus, ierr
	ierr = 0
end

subroutine mpi_irecv( buf, icount, itype, isource, itag, icomm, irequest, ierr )
	integer :: icount, itype, isource, itag, icomm, irequest, ierr
	ierr = 0
end

subroutine mpi_wait( irequest, istatus, ierr )
	integer :: irequest, istatus, ierr
	ierr = 0
end

subroutine mpi_comm_dup( icomm1, icomm2, ierr )
        integer :: icomm1, icomm2, ierr
        icomm2 = icomm1
        ierr  = 0
end

subroutine mpi_comm_split( icomm_in, icolor, ikey, icomm_out, ierr )
        integer :: icomm_in, icolor, ikey, icomm_out, ierr
        icomm_out = icomm_in
        ierr = 0
end

subroutine mpi_comm_free(  icomm1, ierr )
        integer :: icomm1, ierr
        ierr  = 0
end

subroutine mpi_bcast( buf, icount, itype, iroot, icomm, ierr )
	integer :: icount, itype, iroot, icomm, ierr
	ierr = 0
end

subroutine mpi_allgather( send_buf, isend_count, isend_type, &
           &  recv_buf, irecv_counts, irecv_type, icomm, ierr )
	integer    :: isend_count, isend_type
	integer    :: irecv_counts, irecv_type
	integer    :: icomm, ierr
	logical(4) :: send_buf(isend_count*2), recv_buf(isend_count*2)
	call mpi_copy( send_buf, recv_buf, isend_count, itype, ierr )
end

subroutine mpi_allgatherv( send_buf, isend_count, isend_type, &
           &  recv_buf, irecv_counts, idispls, irecv_type, icomm, ierr )
	integer    :: isend_count, isend_type
	integer    :: irecv_counts, idispls, irecv_type
	integer    :: icomm, ierr
	logical(4) :: send_buf(isend_count*2), recv_buf(isend_count*2)
	call mpi_copy( send_buf, recv_buf, isend_count, itype, ierr )
end

subroutine mpi_allreduce( send_buf, recv_buf, icount, itype, iop, icomm, ierr )
	integer    :: icount, itype, iop, icomm, ierr
	logical(4) :: send_buf(icount*2), recv_buf(icount*2)
	call mpi_copy( send_buf, recv_buf, icount, itype, ierr )
end

subroutine mpi_alltoall(send_buf,icount,itype,iroot,icomm,ierr)
    include 'mpif.h'
    integer    :: icount, itype, iroot, icomm, ierr
	logical(4) :: send_buf(icount*2)
end

subroutine mpi_copy( send_buf, recv_buf, icount, itype, ierr )
	include 'mpif.h'
	integer    :: icount, itype, ierr
	logical(4) :: send_buf(icount*2), recv_buf(icount*2)
	integer    :: i,ii
        ii = send_buf(1)
        if(ii .eq. MPI_IN_PLACE)then
            return
        endif

	if ( itype == MPI_DOUBLE_PRECISION ) then
		do i = 1, icount*2
			recv_buf(i) = send_buf(i)
		end do
		ierr = 0
	else if ( itype == MPI_INTEGER ) then
		do i = 1, icount
			recv_buf(i) = send_buf(i)
		end do
		ierr = 0
        else if ( itype == MPI_LOGICAL ) then
		do i = 1, icount
			recv_buf(i) = send_buf(i)
		end do
	else
		print *, 'Error in mpi_copy'
		print *, itype
		ierr = 1
	end if
end

subroutine mpi_sendrecv(sendbuf,sendcount,sendtype,dest,sendtag,recvbuf,recvcount,recvtype,source,recvtag,comm,status,ierr)
        integer :: sendcount,sendtype,dest,sendtag,recvcount,recvtype,source,recvtag,comm,status,ierr
        logical(4) :: sendbuf(sendcount*2),recvbuf(recvcount*2)
        call mpi_copy( sendbuf, recvbuf, sendcount, sendtype, ierr )
end

#ifndef _NO_ARG_DUMMY_
function iargc()
	iargc = 0
end

subroutine getarg(i,buf)
	integer      :: i
	character(*) :: buf
end
#endif

#else

subroutine mpi_dummy()
end

#endif
