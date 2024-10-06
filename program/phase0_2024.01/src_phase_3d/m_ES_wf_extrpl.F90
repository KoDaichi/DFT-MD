#define EXTRPL_DGEMM
  ! ========================================================================================
  ! extrapolate the wfs according to the algorithm described in Arias et al
  ! PRB45 1538 (1992)
  ! ========================================================================================

module m_ES_wf_extrpl

use m_Const_Parameters, only : DP,CMPLDP,ON,GAMMA,EXECUT
use m_Files, only : nfout
use m_Control_Parameters, only : ipripredictor, printable, sw_wf_predictor &
 & ,kimg, nspin, neg, rms_threshold
use m_ES_ortho, only : np_g1k_x,np_fs_x &
     & ,m_ESortho_set_np_g1k_x, m_ESortho_set_np_fs_x
use m_Kpoints, only : kv3,k_symmetry
use m_Timing, only : tstatc0_begin,tstatc0_end
use m_Electronic_Structure, only : zaj_l,fsr_l,fsi_l
use m_PlaneWaveBasisSet, only : kg1
use m_Parallelization, only : np_e, ista_k,iend_k, map_k,myrank_k,mpi_k_world, myrank_e,np_fs,np_g1k
use m_Ionic_System, only : m_IS_reset_extrpl_status
use m_PseudoPotential, only : modnrm,nlmta1_p,nlmta2_p,fqwei_p,nac_p,nlmta

! ============================ KT_Test ============================ 12.5Exp
#ifdef USE_ZAJ_HISTORY
use m_Control_Parameters,  only : af
use m_Parallelization,  only : MPI_CommGroup, ierr, mype, map_ek, map_z
#endif
! ================================================================= 12.5Exp
use m_Parallelization, only : myrank_g, neg_g, mpi_kg_world, mpi_ke_world
#ifdef _USE_SCALAPACK_
use m_Parallelization, only : nrank_g, nrank_e
#endif
use mpi
implicit none

!include 'mpif.h'
integer istatus(mpi_status_size)

real(kind=DP),allocatable, private, dimension(:,:,:,:) :: zaj_l_diff
real(kind=DP),allocatable, private, dimension(:,:,:,:) :: zaj_l_o
real(kind=DP),allocatable, private, dimension(:,:,:) :: fsr_l_o,fsi_l_o
complex(kind=CMPLDP),allocatable, private, dimension(:,:) :: almat
complex(kind=CMPLDP),private,parameter :: zero=(0.d0,0.d0)
complex(kind=CMPLDP),private,parameter :: one=(1.d0,0.d0)
#ifdef _USE_SCALAPACK_
logical, save :: first_call = .true.
integer, save :: nprow, npcol, block_size = 32, lda, occ
integer, save :: ictxt, myrow, mycol
integer, allocatable, dimension(:,:), save :: usermap
integer, dimension(9), save :: desca, descz
integer, save :: lwork_ = -1, lrwork_ = -1, liwork_ = -1
#endif

contains

  subroutine m_ES_wf_extrpl_alloc()
    if(sw_wf_predictor==ON)then
       allocate(almat(neg,neg));almat=zero
       allocate(zaj_l_o(maxval(np_g1k),np_e,ista_k:iend_k,kimg));zaj_l_o=0.d0
       allocate(zaj_l_diff(maxval(np_g1k),np_e,ista_k:iend_k,kimg));zaj_l_diff=0.d0
       allocate(fsr_l_o(np_e,np_fs,ista_k:iend_k));fsr_l_o=0.d0
       if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
           allocate(fsi_l_o(np_e,np_fs,ista_k:iend_k));fsi_l_o=0.d0
       else
           allocate(fsi_l_o(1,1,1));fsi_l_o=0.d0
       endif
    endif
  end subroutine m_ES_wf_extrpl_alloc

  subroutine m_ES_wf_extrpl_dealloc()
    if(allocated(zaj_l_o)) deallocate(zaj_l_o)
    if(allocated(zaj_l_diff)) deallocate(zaj_l_diff)
    if(allocated(almat)) deallocate(almat)
    if(allocated(fsr_l_o)) deallocate(fsr_l_o)
    if(allocated(fsi_l_o)) deallocate(fsi_l_o)
  end subroutine m_ES_wf_extrpl_dealloc

! ===================== KT_TEST ================================= 12.5Exp
#ifdef USE_ZAJ_HISTORY

  subroutine m_ES_wf_rd_zaj_history
    integer, parameter:: file_id  = 350
!
    integer :: ib, ik, ri
!
    real(kind=DP), allocatable :: tmp1_zaj(:,:), tmp2_zaj(:,:)
    real(kind=DP), allocatable :: tmp1_fs(:),    tmp2_fs(:)
!
    if ( sw_wf_predictor /=ON ) return
!
    allocate( tmp1_zaj(kg1,kimg) ); tmp1_zaj = 0.0d0
    allocate( tmp2_zaj(kg1,kimg) ); tmp1_zaj = 0.0d0
!
    allocate( tmp1_fs(1:nlmta) ); tmp1_fs = 0.0d0
    allocate( tmp2_fs(1:nlmta) ); tmp2_fs = 0.0d0

    call mpi_barrier(MPI_CommGroup,ierr)

    do ik = 1, kv3, af+1
       do ib = 1, neg
          if (mype == 0) read(file_id) tmp1_zaj, tmp2_zaj, tmp1_fs, tmp2_fs

          if (mype == 0 .and. map_ek(ib,ik) /= 0) then ! MPI
             call mpi_send( tmp1_zaj, kg1*kimg, mpi_real, map_ek(ib,ik), 1, &
                  &         MPI_CommGroup, ierr )
             call mpi_send( tmp2_zaj, kg1*kimg, mpi_real, map_ek(ib,ik), 1, &
                  &          MPI_CommGroup, ierr )
             call mpi_send( tmp1_fs,  nlmta, mpi_real, map_ek(ib,ik), 1, &
                  &          MPI_CommGroup, ierr )
             call mpi_send( tmp2_fs,  nlmta, mpi_real, map_ek(ib,ik), 1, &
                  &          MPI_CommGroup, ierr )

          else if(map_ek(ib,ik) == mype .and. map_ek(ib,ik) /= 0) then
             call mpi_recv( tmp1_zaj, kg1*kimg, mpi_real, 0, &
                  &         1, MPI_CommGroup, istatus,ierr )!MPI
             call mpi_recv( tmp2_zaj, kg1*kimg, mpi_real, 0, &
                  &         1, MPI_CommGroup, istatus,ierr )!MPI
             call mpi_recv( tmp1_fs, nlmta, mpi_real, 0, &
                  &         1, MPI_CommGroup, istatus,ierr )!MPI
             call mpi_recv( tmp2_fs, nlmta, mpi_real, 0, &
                  &         1, MPI_CommGroup, istatus,ierr )!MPI
          endif

          if(map_ek(ib,ik) == mype) then              ! MPI
             do ri = 1, kimg
                zaj_l_o(1:kg1,map_z(ib),ik,ri)    = tmp1_zaj(1:kg1,ri)  ! MPI
                zaj_l_diff(1:kg1,map_z(ib),ik,ri) = tmp2_zaj(1:kg1,ri)  ! MPI
             end do
              fsr_l_o(map_z(ib),1:nlmta,ik) = tmp1_fs(1:nlmta)
              fsi_l_o(map_z(ib),1:nlmta,ik) = tmp2_fs(1:nlmta)
          end if

       end do
    end do

    deallocate(tmp1_zaj); deallocate(tmp2_zaj)
    deallocate(tmp1_fs); deallocate(tmp2_fs)

  end subroutine m_ES_wf_rd_zaj_history

  subroutine m_ES_wf_wd_zaj_history
    integer, parameter :: file_id = 350
!
    integer :: ri, ib, ik
!
    real(kind=DP), allocatable :: tmp1_zaj(:,:), tmp2_zaj(:,:)
    real(kind=DP), allocatable :: tmp1_fs(:),    tmp2_fs(:)

    if (sw_wf_predictor/=ON) return

    allocate( tmp1_zaj(kg1,kimg) ); tmp1_zaj = 0.0d0
    allocate( tmp2_zaj(kg1,kimg) ); tmp1_zaj = 0.0d0
!
    allocate( tmp1_fs(1:nlmta) ); tmp1_fs = 0.0d0
    allocate( tmp2_fs(1:nlmta) ); tmp2_fs = 0.0d0

    call mpi_barrier(MPI_CommGroup,ierr)

    do ik = 1, kv3, af+1
       do ib = 1, neg
          if(map_ek(ib,ik) == mype) then                          ! MPI
             do ri = 1, kimg
                tmp1_zaj(1:kg1,ri) = zaj_l_o(1:kg1,map_z(ib),ik,ri)
                tmp2_zaj(1:kg1,ri) = zaj_l_diff(1:kg1,map_z(ib),ik,ri)
             end do
             tmp1_fs(1:nlmta) = fsr_l_o(map_z(ib),1:nlmta,ik)
             tmp2_fs(1:nlmta) = fsi_l_o(map_z(ib),1:nlmta,ik)

             if (map_ek(ib,ik) /= 0) then
                call mpi_send( tmp1_zaj, kg1*kimg, mpi_real, 0, 1, &
                     &          MPI_CommGroup, ierr )
                call mpi_send( tmp2_zaj, kg1*kimg, mpi_real, 0, 1, &
                     &          MPI_CommGroup, ierr )
                call mpi_send( tmp1_fs,  nlmta, mpi_real, 0, 1, &
                     &          MPI_CommGroup, ierr )
                call mpi_send( tmp2_fs,  nlmta, mpi_real, 0, 1, &
                     &          MPI_CommGroup, ierr )
             endif

          else if(mype == 0 .and. map_ek(ib,ik) /= 0) then
             call mpi_recv( tmp1_zaj, kg1*kimg, mpi_real, map_ek(ib,ik), &
                  &         1, MPI_CommGroup, istatus,ierr )!MPI
             call mpi_recv( tmp2_zaj, kg1*kimg, mpi_real, map_ek(ib,ik), &
                  &         1, MPI_CommGroup, istatus,ierr )!MPI
             call mpi_recv( tmp1_fs, nlmta, mpi_real, map_ek(ib,ik), &
                  &         1, MPI_CommGroup, istatus,ierr )!MPI
             call mpi_recv( tmp2_fs, nlmta, mpi_real, map_ek(ib,ik), &
                  &         1, MPI_CommGroup, istatus,ierr )!MPI
          end if
          if (mype == 0) then
             write(file_id) tmp1_zaj, tmp2_zaj, tmp1_fs, tmp2_fs
          endif
       end do
    end do

    deallocate(tmp1_zaj); deallocate(tmp2_zaj)
    deallocate(tmp1_fs); deallocate(tmp2_fs)

  end subroutine m_ES_wf_wd_zaj_history
#endif
! =============================================================== 12.5Exp

  subroutine m_ES_wf_extrpl_doit(alpha,beta,rms,nextpl)
     real(kind=DP),intent(in) :: alpha,beta,rms
     integer, intent(in) :: nextpl
     real(kind=DP) :: r
     integer :: is,ik
     real(kind=DP),allocatable,dimension(:,:,:) :: wf_t
     complex(kind=CMPLDP),allocatable,dimension(:,:) :: wf1,wf2_in,wf3_in,wf2,wf3
     integer :: id_sname = -1
     real(kind=DP),allocatable,dimension(:,:) :: bpr_t,bpi_t,bpr_t_o,bpi_t_o
     integer :: stat, ierr
     integer :: i, ib
#ifdef _USE_SCALAPACK_
     integer :: max_block_size, j, nb, nprocs, myrank, info
#endif
     call tstatc0_begin('m_ES_wf_extrpl ',id_sname,1)
     if(nextpl<2)then ! not ready
        zaj_l_diff(:,:,:,:) = zaj_l(:,:,:,:)-zaj_l_o(:,:,:,:)
        zaj_l_o(:,:,:,:) = zaj_l(:,:,:,:)
        if(modnrm==EXECUT)then
           fsr_l_o(:,:,:) = fsr_l(:,:,:)
           fsi_l_o(:,:,:) = fsi_l(:,:,:)
        endif
        return
     endif
#ifdef _USE_SCALAPACK_
     ! ScaLAPACK setup.
     if(first_call) then
        if(ipripredictor>1) write(nfout,'(a)') 'ScaLAPACK setup in m_ES_wf_extrpl_doit!'

        ! Define nprow and npcol.
        nprow = int(sqrt(real(nrank_e*nrank_g)))
        do
           if(nprow <= 1) exit
           if(mod(nrank_g*nrank_e,nprow) == 0) exit
           nprow = nprow - 1
        end do
        npcol = nrank_g*nrank_e/nprow

        ! Define block size.
        max_block_size = neg/max(nprow, npcol)
        if(block_size > max_block_size) block_size = max_block_size

        if(ipripredictor>1) write(nfout,'(a,3i8)') 'nprow, npcol, block_size: ', nprow, npcol, block_size

        ! Make usermap for BLACS.
        allocate(usermap(nprow, npcol))
        do j = 1, npcol
           do i = 1, nprow
              usermap(i, j) = myrank_k*nrank_e*nrank_g + (i-1)*npcol + (j-1)
           end do
        end do

        nb = block_size

        ! Define ScaLAPACK array size (lda, occ).
        ! - Number of blocks.
        if(mod(neg, nb) > 0) then
           lda = neg/nb + 1
        else
           lda = neg/nb
        end if
        occ = lda
        ! - Number of blocks/proc.
        if(mod(lda, nprow) > 0) then
           lda = lda/nprow + 1
        else
           lda = lda/nprow
        end if
        if(mod(occ, npcol) > 0) then
           occ = occ/npcol + 1
        else
           occ = occ/npcol
        end if
        lda = lda*nb
        occ = occ*nb

        if(ipripredictor>1) write(nfout,'(a,2i8)') 'ScaLAPACK array size (lda, occ): ', lda, occ

        ! Make BLACS context.
        call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
        call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)
        call blacs_setup(myrank, nprocs)

        call blacs_get(-1, 0, ictxt)
        call blacs_gridmap(ictxt, usermap, nprow, nprow, npcol)
        call blacs_gridinfo(ictxt, nprow, npcol, myrow, mycol)

        ! Make descripters.
        call descinit(desca, neg, neg, nb, nb, 0, 0, ictxt, lda, info)
        if(ipripredictor>1) write(nfout,'(a,i3)') 'descinit for amat: info = ', info
        call descinit(descz, neg, neg, nb, nb, 0, 0, ictxt, lda, info)
        if(ipripredictor>1) write(nfout,'(a,i3)') 'descinit for zmat: info = ', info
     end if
#endif

     if(modnrm==EXECUT)then
        call m_ESortho_set_np_fs_x() !-> np_fs_x
        allocate(bpr_t  (np_fs,neg))
        allocate(bpi_t  (np_fs,neg))
        allocate(bpr_t_o(np_fs,neg))
        allocate(bpi_t_o(np_fs,neg))
     endif
     do is=1,nspin
       do ik = is, kv3-nspin+is, nspin
          if(map_k(ik) /= myrank_k) cycle                   ! MPI
          call m_ESortho_set_np_g1k_x(ik) !-> np_g1k_x
          allocate(wf1   (maxval(np_g1k),neg)); wf1    = zero
          allocate(wf2   (maxval(np_g1k),neg)); wf2    = zero
          allocate(wf3   (maxval(np_g1k),neg)); wf3    = zero
          allocate(wf2_in(maxval(np_g1k),neg)); wf2_in = zero
          allocate(wf3_in(maxval(np_g1k),neg)); wf3_in = zero

          if(kimg == 1)then
             do ib = 1, np_e
                do i = 1, np_g1k(ik)
                   wf1   (i,neg_g(ib)) = dcmplx(zaj_l     (i,ib,ik,1),0.d0)
                   wf2_in(i,neg_g(ib)) = dcmplx(zaj_l_o   (i,ib,ik,1),0.d0)
                   wf3_in(i,neg_g(ib)) = dcmplx(zaj_l_diff(i,ib,ik,1),0.d0)
                end do
             end do
          else
             do ib = 1, np_e
                do i = 1, np_g1k(ik)
                   wf1   (i,neg_g(ib)) = dcmplx(zaj_l     (i,ib,ik,1),zaj_l     (i,ib,ik,2))
                   wf2_in(i,neg_g(ib)) = dcmplx(zaj_l_o   (i,ib,ik,1),zaj_l_o   (i,ib,ik,2))
                   wf3_in(i,neg_g(ib)) = dcmplx(zaj_l_diff(i,ib,ik,1),zaj_l_diff(i,ib,ik,2))
                end do
             end do
          end if
          call mpi_allreduce(MPI_IN_PLACE,wf1,   maxval(np_g1k)*neg,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_kg_world,ierr)
          call mpi_allreduce(MPI_IN_PLACE,wf2_in,maxval(np_g1k)*neg,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_kg_world,ierr)
          call mpi_allreduce(MPI_IN_PLACE,wf3_in,maxval(np_g1k)*neg,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_kg_world,ierr)

          if(modnrm == EXECUT) then
             bpr_t   = 0.d0
             bpi_t   = 0.d0
             bpr_t_o = 0.d0
             bpi_t_o = 0.d0

             do ib = 1, np_e
                do i = 1, np_fs
                   bpr_t  (i,neg_g(ib)) = fsr_l  (ib,i,ik)
                   bpr_t_o(i,neg_g(ib)) = fsr_l_o(ib,i,ik)
                end do
             end do
             call mpi_allreduce(MPI_IN_PLACE,bpr_t,  np_fs*neg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
             call mpi_allreduce(MPI_IN_PLACE,bpr_t_o,np_fs*neg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)

             if(.not. (kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
                do ib = 1, np_e
                   do i = 1, np_fs
                      bpi_t  (i,neg_g(ib)) = fsi_l  (ib,i,ik)
                      bpi_t_o(i,neg_g(ib)) = fsi_l_o(ib,i,ik)
                   end do
                end do
                call mpi_allreduce(MPI_IN_PLACE,bpi_t,  np_fs*neg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
                call mpi_allreduce(MPI_IN_PLACE,bpi_t_o,np_fs*neg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)
             end if
          end if

          ! align the wfs
          call check_wf_diff_norm(wf1,wf2_in,r)
          if(ipripredictor>=2.and.printable)then
             write(nfout,fmt='(a,i8,a,f20.10)') &
          & 'norm of the difference between the new and old wfs before alignment for kpoint no. ',ik,' : ',r
          endif
          call get_subspace_alignment_matrix(stat)
          if(stat/=0)then
             if(printable) write(nfout,'(a)') '!** WARN failed to obtain the subspace alignment matrix; &
         &   resetting extrapolation status'
!!x!!$             write(nfout,'(a)') '!** WARN failed to obtain the subspace alignment matrix; &
!!x!!$                  &   resetting extrapolation status'
!!x!!$             call flush(nfout)
!!x!!$             stop
             call m_IS_reset_extrpl_status()
             return
          endif

          ! apply the subspace alignment matrix to the wfs
          call zgemm('N','N',np_g1k(ik),neg,neg,one,wf2_in,maxval(np_g1k),almat,neg,zero,wf2,maxval(np_g1k))
          call zgemm('N','N',np_g1k(ik),neg,neg,one,wf3_in,maxval(np_g1k),almat,neg,zero,wf3,maxval(np_g1k))

          call check_wf_diff_norm(wf1,wf2,r)
          if(ipripredictor>=2.and.printable)then
             write(nfout,fmt='(a,i8,a,f20.10)') &
          & 'norm of the difference between the new and old wfs after  alignment for kpoint no. ',ik,' : ',r
          endif

          do ib = 1, np_e
             do i = 1, np_g1k(ik)
                zaj_l_o   (i,ib,ik,1) = dble(wf1(i,neg_g(ib)))
                zaj_l_diff(i,ib,ik,1) = dble(wf1(i,neg_g(ib))) - dble(wf2(i,neg_g(ib)))
             end do
          end do
          if(kimg > 1)then
             do ib = 1, np_e
                do i = 1, np_g1k(ik)
                   zaj_l_o   (i,ib,ik,2) = dimag(wf1(i,neg_g(ib)))
                   zaj_l_diff(i,ib,ik,2) = dimag(wf1(i,neg_g(ib))) - dimag(wf2(i,neg_g(ib)))
                end do
             end do
          end if

          if(modnrm == EXECUT) then
             do ib = 1, np_e
                do i = 1, np_fs
                   fsr_l_o(ib,i,ik) = bpr_t(i,neg_g(ib))
                   fsr_l  (ib,i,ik) = bpr_t(i,neg_g(ib))
                end do
             end do

             if(.not. (kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
                do ib = 1, np_e
                   do i = 1, np_fs
                      fsi_l_o(ib,i,ik) = bpi_t(i,neg_g(ib))
                      fsi_l  (ib,i,ik) = bpi_t(i,neg_g(ib))
                   end do
                end do
             end if
          end if

          if(rms < rms_threshold)then
             do ib = 1, np_e
                do i = 1, np_g1k(ik)
                   zaj_l(i,ib,ik,1) = dble(wf1(i,neg_g(ib))) &
                      + alpha*(dble(wf1(i,neg_g(ib)))-dble(wf2(i,neg_g(ib)))) &
                      + beta*dble(wf3(i,neg_g(ib)))
                end do
             end do
             if(kimg > 1)then
                do ib = 1, np_e
                   do i = 1, np_g1k(ik)
                      zaj_l(i,ib,ik,2) = dimag(wf1(i,neg_g(ib))) &
                         + alpha*(dimag(wf1(i,neg_g(ib)))-dimag(wf2(i,neg_g(ib)))) &
                         + beta*dimag(wf3(i,neg_g(ib)))
                   end do
                end do
             end if
          end if

          deallocate(wf1)
          deallocate(wf2)
          deallocate(wf3)
          deallocate(wf2_in)
          deallocate(wf3_in)
       enddo
     enddo
     if(modnrm==EXECUT)then
       deallocate(bpr_t)
       deallocate(bpi_t)
       deallocate(bpr_t_o)
       deallocate(bpi_t_o)
     endif
     call tstatc0_end(id_sname)

  contains

  ! ========================================================================================
  ! calculate the subspace alignment matrix, as described in Arias et.al. PRB 45 1538 (1992).
  ! i.   Uno =  <new|S|old>
  ! ii.  UU = Uno^* x Uno
  ! iii. diagonalize UU -> eigvec => A^*, eigval => E
  ! iv.  A' = (1/sqrt(E)) x A U^*
  ! v.   the subspace alignment matrix is A'x A^*
  !
  ! note : ^* denotes complex conjugate, x dentoes matrix-matrix multiplication
  ! ========================================================================================
  subroutine get_subspace_alignment_matrix(stat)
    integer, intent(out) :: stat
    complex(kind=CMPLDP),allocatable,dimension(:,:) :: Amat,APmat,Umat,UUmat,eigvec
    real(kind=DP), allocatable, dimension(:) :: eigval
    integer,save :: lwork=-1
    integer :: ig,ib1,ib2,il,iu,info,nn,ierr,ig1,ia,p,q
    real(kind=DP) :: vl,vu,abstol,dlamch,ar,ai
    complex(kind=CMPLDP), allocatable, dimension(:) :: work
    complex(kind=CMPLDP) :: cbp,cbpo
    real(kind=DP),allocatable,dimension(:) :: rwork
    integer,allocatable,dimension(:) :: iwork,ifail
    real(kind=DP) :: eiginv,eigmin
#ifdef EXTRPL_DGEMM
    complex(kind=CMPLDP),allocatable,dimension(:,:) :: tmp1, tmp2
#endif
#ifdef _USE_SCALAPACK_
    complex(kind=CMPLDP),allocatable,dimension(:,:) :: UUmat_, eigvec_
    complex(kind=CMPLDP), allocatable, dimension(:) :: work_
    real(kind=DP),allocatable,dimension(:) :: eigval_, rwork_
    integer,allocatable,dimension(:) :: iwork_
! === This should be removed later!!! ====================================================
    integer :: i_, j_, ist, jst, iind, jind
! ========================================================================================
#endif
    allocate(Umat(neg,neg));Umat=zero
    allocate(UUmat(neg,neg));UUmat=zero
    if(lwork<0) lwork = 2*neg
    allocate(eigvec(neg,neg));eigvec=zero
    allocate(eigval(neg));eigval=0.d0
#ifdef _USE_SCALAPACK_
    allocate(UUmat_(lda, occ))
    allocate(eigvec_(lda, occ))
    allocate(eigval_(neg))
#endif
    allocate(rwork(7*neg));rwork=0.d0
    allocate(iwork(5*neg));iwork=0
    allocate(ifail(neg));ifail=0
    allocate(work(1:lwork));work=zero
    allocate(Amat(neg,neg));Amat=zero
    allocate(APmat(neg,neg));APmat=zero
    almat = zero
    do ib1=1,neg
       almat(ib1,ib1) = one
    enddo
#ifdef EXTRPL_DGEMM
    allocate(tmp1(nac_p, neg))
    allocate(tmp2(nac_p, neg))
#endif
#ifdef _USE_SCALAPACK_
    if(first_call) then
       allocate(work_ (1))
       allocate(rwork_(1))
       allocate(iwork_(1))
       call pzheevd('V', 'U', neg, UUmat_, 1, 1, desca, eigval_, eigvec_, 1, 1, descz, &
                    work_, lwork_, rwork_, lrwork_, iwork_, liwork_, info)
       lwork_  = int(dble(work_(1)))
       lrwork_ = int(rwork_(1))
       liwork_ = int(iwork_(1))
       deallocate(work_)
       deallocate(rwork_)
       deallocate(iwork_)

       first_call = .false.
    end if
    allocate(work_ (lwork_))
    allocate(rwork_(lrwork_))
    allocate(iwork_(liwork_))
#endif
    ! construct U=<new|S|old>
#ifndef EXTRPL_DGEMM
    do ib1=1,neg
      do ib2=1,neg
         ig1=1
         if(k_symmetry(ik)==GAMMA.and.myrank_g==0) ig1 = 2
         do ig=ig1,np_g1k(ik)
            Umat(ib1,ib2) = Umat(ib1,ib2) + dconjg(wf1(ig,ib1))*wf2_in(ig,ib2)
         enddo
         if(k_symmetry(ik)==GAMMA) Umat(ib1,ib2) = 2.d0*Umat(ib1,ib2)
         if(ig1==2) Umat(ib1,ib2) = Umat(ib1,ib2)+dconjg(wf1(1,ib1))*wf2_in(1,ib2)
         if(modnrm==EXECUT)then ! contribution from the overlap matrix
            do ia=1,nac_p
               p = nlmta1_p(ia);q = nlmta2_p(ia)
               if(k_symmetry(ik)==GAMMA)then
                  cbp=dcmplx(bpr_t(p,ib1),0.d0)
                  cbpo=dcmplx(bpr_t_o(q,ib2),0.d0)
               else
                  cbp=dcmplx(bpr_t(p,ib1),bpi_t(p,ib1))
                  cbpo=dcmplx(bpr_t_o(q,ib2),bpi_t_o(q,ib2))
               endif
               Umat(ib1,ib2) = Umat(ib1,ib2) + fqwei_p(ia)*dconjg(cbp)*cbpo
            enddo
         endif
         if(k_symmetry(ik)==GAMMA) Umat(ib1,ib2) = dcmplx(dble(Umat(ib1,ib2)),0.d0)
      enddo
    enddo
#else
    if(k_symmetry(ik) == GAMMA) then
       call zgemm('C','N',neg,neg,np_g1k(ik),(2.0d0,0.0d0), &
                  wf1,maxval(np_g1k),wf2_in,maxval(np_g1k),(0.0d0,0.0d0),Umat,neg)
       if(myrank_g == 0) then
          do ib1 = 1, neg
             do ib2 = 1, neg
                Umat(ib1, ib2) = Umat(ib1, ib2) - dconjg(wf1(1, ib1))*wf2_in(1, ib2)
             end do
          end do
       end if
    else
       call zgemm('C','N',neg,neg,np_g1k(ik),(1.0d0,0.0d0), &
                  wf1,maxval(np_g1k),wf2_in,maxval(np_g1k),(0.0d0,0.0d0),Umat,neg)
    end if

    if(modnrm == EXECUT) then ! contribution from the overlap matrix
! === DEBUG by tkato 2014/07/29 ================================================
       if(nac_p > 0) then
! ==============================================================================
       if(k_symmetry(ik) == GAMMA) then
          do ib1 = 1, neg
             do ia = 1, nac_p
                p = nlmta1_p(ia); q = nlmta2_p(ia)
                tmp1(ia, ib1) = fqwei_p(ia)*dcmplx(bpr_t  (p, ib1), 0.0d0)
                tmp2(ia, ib1) =             dcmplx(bpr_t_o(q, ib1), 0.0d0)
             end do
          end do
       else
          do ib1 = 1, neg
             do ia = 1, nac_p
                p = nlmta1_p(ia); q = nlmta2_p(ia)
                tmp1(ia, ib1) = fqwei_p(ia)*dcmplx(bpr_t  (p,ib1), bpi_t  (p,ib1))
                tmp2(ia, ib1) =             dcmplx(bpr_t_o(q,ib1), bpi_t_o(q,ib1))
             end do
          end do
       end if
       call zgemm('C','N',neg,neg,nac_p,(1.0d0,0.0d0), &
                  tmp1,nac_p,tmp2,nac_p,(1.0d0,0.0d0),Umat,neg)
       if(k_symmetry(ik) == GAMMA) then
          do ib1 = 1, neg
             do ib2 = 1, neg
                Umat(ib1,ib2) = dcmplx(dble(Umat(ib1, ib2)), 0.d0)
             end do
          end do
       end if
! === DEBUG by tkato 2014/07/29 ================================================
       end if
! ==============================================================================
    end if
#endif
    call mpi_allreduce(MPI_IN_PLACE,Umat,neg*neg,mpi_double_complex,mpi_sum,mpi_ke_world,ierr)

    ! calculate UU=U^* U
    call zgemm('C','N',neg,neg,neg,one,Umat,neg,Umat,neg,zero,UUmat,neg)

#ifdef _USE_SCALAPACK_
! === This should be removed later!!! ====================================================
    ist = myrow*block_size + 1
    jst = mycol*block_size + 1
    jind = 0
    do j_ = jst, neg, block_size*npcol
       do ib2 = j_, j_+block_size-1
          if(ib2 <= neg) then
             jind = jind + 1
             iind = 0
             do i_ = ist, neg, block_size*nprow
                do ib1 = i_, i_+block_size-1
                   if(ib1 <= neg) then
                      iind = iind + 1
                      UUmat_(iind,jind) = UUmat(ib1,ib2)
                   end if
                enddo
             enddo
          end if
       enddo
    enddo
! ========================================================================================
#endif
    ! diagonalize UU
    vl=0.d0;vu=0.d0;il=0;iu=0
    abstol = 2*dlamch('S')
#ifndef _USE_SCALAPACK_
    call zheevx('V','A','U',neg,UUmat,neg,vl,vu,il,iu,abstol,nn,eigval,eigvec,neg, &
         &  work,lwork,rwork,iwork,ifail,info)
    if(info==0) then
       lwork = int(work(1))
    else
       stat = info
       return
    endif
#else
    call pzheevd('V', 'U', neg, UUmat_, 1, 1, desca, eigval_, eigvec_, 1, 1, descz, &
                 work_, lwork_, rwork_, lrwork_, iwork_, liwork_, info)
    if(info /= 0) then
       stat = info
       return
    endif
#endif
#ifdef _USE_SCALAPACK_
! === This should be removed later!!! ====================================================
    eigval = eigval_

    eigvec = (0.0d0, 0.0d0)
    ist = myrow*block_size + 1
    jst = mycol*block_size + 1
    jind = 0
    do j_ = jst, neg, block_size*npcol
       do ib2 = j_, j_+block_size-1
          if(ib2 <= neg) then
             jind = jind + 1
             iind = 0
             do i_ = ist, neg, block_size*nprow
                do ib1 = i_, i_+block_size-1
                   if(ib1 <= neg) then
                      iind = iind + 1
                      eigvec(ib1,ib2) = eigvec_(iind,jind)
                   end if
                enddo
             enddo
          end if
       enddo
    enddo
!    call mpi_allreduce(MPI_IN_PLACE,eigvec,neg*neg,mpi_double_complex,mpi_sum,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,eigvec,neg*neg,mpi_double_complex,mpi_sum,mpi_ke_world,ierr)
! ========================================================================================
#endif
    eigmin = minval(dabs(eigval))
    if(eigmin<1.d-15)then
       if(printable) write(nfout,'(a,f15.10)') ' !** minimum eigen value is too small : ',eigmin
       stat = 1
       return
    endif

    ! eigvec^* -> Amat
    do ib1=1,neg
      do ib2=1,neg
        Amat(ib2,ib1) = dconjg(eigvec(ib1,ib2))
      enddo
    enddo

    ! calculate APmat = (1/sqrt(E)) x A x U^*
    call zgemm('N','C',neg,neg,neg,one,Amat,neg,Umat,neg,zero,APmat,neg)
    do ib1=1,neg
      eiginv = 1.d0/dsqrt(eigval(ib1))
      do ib2=1,neg
        APmat(ib1,ib2) = eiginv*APmat(ib1,ib2)
      enddo
    enddo

    ! calculate Ap^* x A and store it to Umat
    Umat = zero
    call zgemm('C','N',neg,neg,neg,one,APmat,neg,Amat,neg,zero,Umat,neg)

    ! the c.c. of Umat is the subspace alignment matrix we want
    do ib1=1,neg
      do ib2=1,neg
        almat(ib1,ib2) = dconjg(Umat(ib2,ib1))
      enddo
    enddo
    stat = 0
    deallocate(Umat)
    deallocate(UUmat)
    deallocate(eigvec)
    deallocate(eigval)
    deallocate(rwork)
    deallocate(iwork)
    deallocate(ifail)
    deallocate(work)
    deallocate(Amat)
    deallocate(APmat)
#ifdef EXTRPL_DGEMM
    deallocate(tmp1)
    deallocate(tmp2)
#endif
#ifdef _USE_SCALAPACK_
    deallocate(UUmat_)
    deallocate(eigvec_)
    deallocate(eigval_)
    deallocate(rwork_)
    deallocate(iwork_)
    deallocate(work_)
#endif
  end subroutine get_subspace_alignment_matrix

  subroutine check_wf_diff_norm(wf1,wf2,rsum)
    complex(kind=CMPLDP),dimension(:,:),intent(in) :: wf1,wf2
    real(kind=DP),intent(out) :: rsum
    real(kind=DP) :: r
    complex(kind=CMPLDP) :: cbp,cbpo
    integer :: ib,ig,ierr,ig1,p,q,ia
    rsum=0.d0
    do ib=1,neg
       ig1=1
       r=0.d0
       if(k_symmetry(ik)==GAMMA.and.myrank_g==0) ig1 = 2
       do ig=ig1,np_g1k(ik)
          r = r+dconjg(wf1(ig,ib)-wf2(ig,ib))*(wf1(ig,ib)-wf2(ig,ib))
       enddo
       if(k_symmetry(ik)==GAMMA) r=2*r
       if(ig1==2) r=r+dconjg(wf1(1,ib)-wf2(1,ib))*(wf1(1,ib)-wf2(1,ib))
       rsum = rsum+r
    enddo
    call mpi_allreduce(MPI_IN_PLACE,rsum,1,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
  end subroutine check_wf_diff_norm

  ! ========================================================================================
  ! 'align' the wfs by the subspace alignment matrix
  ! this is done to the old wf and the diff of the wfs,
  ! so that they are 'aligned' to the current wf
  ! ========================================================================================
  subroutine align_wfs(wf,wf_aligned)
    complex(kind=CMPLDP),dimension(:,:),intent(in) :: wf
    complex(kind=CMPLDP),dimension(:,:),intent(out) :: wf_aligned
    complex(kind=CMPLDP) :: matelem
    integer :: ib1,ib2,ig,ierr
    wf_aligned = zero
    do ib1=1,neg
      do ib2=1,neg
         matelem = almat(ib2,ib1)
         do ig=1,np_g1k_x
           if(k_symmetry(ik)==GAMMA)then
              wf_aligned(ig,ib1) = wf_aligned(ig,ib1)+dble(matelem)*wf(ig,ib2)
           else
              wf_aligned(ig,ib1) = wf_aligned(ig,ib1)+matelem*wf(ig,ib2)
           endif
         enddo
      enddo
    enddo
  end subroutine align_wfs

  end subroutine m_ES_wf_extrpl_doit

end module m_ES_wf_extrpl
