module z_interface_3D

#ifdef PARA3D
  use m_PlaneWaveBasisSet,  only : kg1,iba,kgp,kg
  use m_PseudoPotential,    only : nlmta,nlmtt

  use m_Const_Parameters,   only : DP, CMPLDP, SKIP, EXECUT, ON, OFF, INVERSE, DELTA, PAI2 &
       &                         , ORTHOGONALIZATION, ORTHONORMALIZATION, NORMALIZATION &
       &                         , NORMCONSERVATION, VDB, VANDERBILT_TYPE &
       &                         , OTHER_BANDS, SAME_BAND, DELTAevdff, BUCS, SCF, EK &
       &                         , ELECTRON, GAMMA, MAPPED, NOTMAPPED, EK_CONVERGED

  use m_Parallelization,    only   : mpi_kg_world   &
       &                           , mpi_ke_world   &
       &                           , MPI_CommGroup &
       &                           , ista_k_3D, iend_k_3D       &
       &                           , np_g1k_3D, nis_g1k_3D, ista_g1k_3D, iend_g1k_3D  &
       &                           , myrank_g_3D      &
       &                           , myrank_e, myrank_e_3D, nrank_e_3D      &
       &                           , nis_e_3D, np_e_3D, nis_fs_3D, np_fs_3D   &
       &                           , np_e, nis_e, ista_k, iend_k, ista_e, iend_e   &
       &                           , map_e, map_z, mype, npes &
       &                           , ista_snl, iend_snl      &
       &                           , mpi_k_world,myrank_k   &
       &                           ,  nbs_sta, nbs_end, neg_g_all, neg_gg_all &
       &                           , ista_kngp, iend_kngp, np_kngp, mp_kngp &
       &                           , ista_kngp_3D, iend_kngp_3D, np_kngp_3D, mp_kngp_3D

  use m_Control_Parameters, only :   nblocksize_mgs          &
       &                           , nblocksize_mgs_is_given &
       &                           , neg, kimg, nspin
  use m_Ionic_System,       only : ntyp

  use m_Timing,              only : tstatc0_begin, tstatc0_end

!index band
  use m_Parallelization,    only   : lrank, nbsn, nbsn_sta, nbsn_end, neg_g, nbs_num, nbsn_num

  implicit none

   integer :: id_sname=-1

contains
  subroutine decomp_zaj_l_3D(zaj_l,zaj_l_3D,ik,nrvf_ordr,sort)

  integer :: ierr,err
  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

  real(kind=DP) :: zaj_l(kg1,np_e,ista_k:iend_k,kimg)
  real(kind=DP) :: zaj_l_3D(maxval(np_g1k_3D),np_e_3D,ista_k_3D:iend_k_3D,kimg)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, ik, is, js, nn, ito
  real(kind=DP), allocatable, dimension(:,:,:,:) :: t_wk1, t_wk2
  integer :: nrvf_ordr(neg,ista_k:iend_k)
  character(len=4) :: sort

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate( t_wk1(kg1,neg,ista_k:iend_k,kimg) )
   allocate( t_wk2(kg1,neg,ista_k:iend_k,kimg) )

   t_wk1   = 0.0d0
   t_wk2   = 0.0d0
   zaj_l_3D = 0.0d0

   do l = 1, kimg
   do k = ista_k, iend_k
   do j = 1,np_e
   do i = 1,iba(k)
     if ((sort == "sort") .and. (k == ik)) then
     ito = nrvf_ordr(nis_e(myrank_e) + j-1,ik)
     t_wk1(i,ito,k,l) = zaj_l(i,j,k,l)
   else
     js = nis_e(myrank_e) + j-1
     t_wk1(i,js,k,l) = zaj_l(i,j,k,l)
   endif
   enddo; enddo; enddo; enddo

   nn = kg1*neg*(iend_k-ista_k+1)*kimg
   call mpi_allreduce(t_wk1,t_wk2,nn,mpi_real8,mpi_sum, mpi_k_world(myrank_k),err)

   do l = 1, kimg
   do k = ista_k_3D, iend_k_3D
   do j = 1, np_e_3D
        js = neg_g(j)
   do i = 1,np_g1k_3D(k)
        is = nis_g1k_3D(myrank_g_3D,k) + i-1
        zaj_l_3D(i,j,k,l) = t_wk2(is,js,k,l)
   enddo; enddo; enddo; enddo


   deallocate( t_wk1, t_wk2)

   call tstatc0_end(id_sname)

  end subroutine decomp_zaj_l_3D

!===============================================================================
  subroutine decomp_zaj_l_r_3D(zaj_l,zaj_l_3D,ik,neordr,sort)

  integer :: ierr,err
  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

  real(kind=DP) :: zaj_l(kg1,np_e,ista_k:iend_k,kimg)
  real(kind=DP) :: zaj_l_3D(maxval(np_g1k_3D),np_e_3D,ista_k_3D:iend_k_3D,kimg)
  integer :: NB, nbs, lb, local_block, jj
  integer :: i, j, k, l, m, n, ik, is, js, nn, ifrom, jto
  real(kind=DP), allocatable, dimension(:,:,:,:) :: t_wk1, t_wk2,t_wk3, t_wk4
  integer :: neordr(neg,ista_k:iend_k)
  character(len=4) :: sort
  logical :: existance
  character(len=3) :: c_num
  character(len=2) :: c_img
  character(len=80) :: fname

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate( t_wk1(kg1,neg,ista_k_3D:iend_k_3D,kimg) )
   allocate( t_wk2(kg1,neg,ista_k_3D:iend_k_3D,kimg) )

   t_wk1   = 0.0d0
   t_wk2   = 0.0d0

   do l = 1, kimg
   do k = ista_k_3D, iend_k_3D
   do j = 1, np_e_3D
        if ((sort == "sort") .and. (k == ik)) then
             jto = neordr(neg_g(j),ik)
        else
             jto = neg_g(j)
        endif
   do i = 1,np_g1k_3D(k)
        is = nis_g1k_3D(myrank_g_3D,k) + i-1
        t_wk2(is,jto,k,l) = zaj_l_3D(i,j,k,l)
   enddo; enddo; enddo; enddo

   nn = kg1*neg*(iend_k_3D-ista_k_3D+1)*kimg
   call mpi_allreduce(t_wk2,t_wk1,nn,mpi_real8,mpi_sum, mpi_k_world(myrank_k),err)

   do l = 1, kimg
   do k = ista_k, iend_k
   do j = 1,np_e
   do i = 1,iba(k)
     js = nis_e(myrank_e) + j-1
     zaj_l(i,j,k,l) = t_wk1(i,js,k,l)
   enddo; enddo; enddo; enddo

   deallocate( t_wk1, t_wk2)

   call tstatc0_end(id_sname)

  end subroutine decomp_zaj_l_r_3D

!===============================================================================
  subroutine decomp_zaj_l_3D_ik(zaj_l,zaj_l_3D,ik,nrvf_ordr,sort)

  integer :: ierr,err
  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

  real(kind=DP) :: zaj_l(kg1,np_e,ista_k:iend_k,kimg)
  real(kind=DP) :: zaj_l_3D(maxval(np_g1k_3D),np_e_3D,ista_k_3D:iend_k_3D,kimg)
  integer :: NB, nbs, lb
  integer :: i, j, l, m, n, ik, is, js, nn, ito
  real(kind=DP), allocatable, dimension(:,:,:) :: t_wk1, t_wk2
  integer :: nrvf_ordr(neg,ista_k:iend_k)
  character(len=4) :: sort

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate( t_wk1(kg1,neg,kimg) )
   allocate( t_wk2(kg1,neg,kimg) )

   t_wk1   = 0.0d0
   t_wk2   = 0.0d0
!  zaj_l_3D = 0.0d0

   do l = 1, kimg
   do j = 1,np_e
   do i = 1,iba(ik)
   if ((sort == "sort")) then
     ito = nrvf_ordr(nis_e(myrank_e) + j-1,ik)
     t_wk1(i,ito,l) = zaj_l(i,j,ik,l)
   else
     js = nis_e(myrank_e) + j-1
     t_wk1(i,js,l) = zaj_l(i,j,ik,l)
   endif
   enddo; enddo; enddo

   nn = kg1*neg*kimg
   call mpi_allreduce(t_wk1,t_wk2,nn,mpi_real8,mpi_sum, mpi_k_world(myrank_k),err)

   do l = 1, kimg
   do j = 1, np_e_3D
        js = neg_g(j)
   do i = 1,np_g1k_3D(ik)
        is = nis_g1k_3D(myrank_g_3D,ik) + i-1
        zaj_l_3D(i,j,ik,l) = t_wk2(is,js,l)
   enddo; enddo; enddo


   deallocate( t_wk1, t_wk2)

   call tstatc0_end(id_sname)

  end subroutine decomp_zaj_l_3D_ik

!===============================================================================
  subroutine decomp_zaj_l_r_3D_ik(zaj_l,zaj_l_3D,ik,neordr,sort)

  integer :: ierr,err
  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

  real(kind=DP) :: zaj_l(kg1,np_e,ista_k:iend_k,kimg)
  real(kind=DP) :: zaj_l_3D(maxval(np_g1k_3D),np_e_3D,ista_k_3D:iend_k_3D,kimg)
  integer :: NB, nbs, lb, local_block, jj
  integer :: i, j, l, m, n, ik, is, js, nn, ifrom, jto
  real(kind=DP), allocatable, dimension(:,:,:) :: t_wk1, t_wk2
  integer :: neordr(neg,ista_k:iend_k)
  character(len=4) :: sort

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate( t_wk1(kg1,neg,kimg) )
   allocate( t_wk2(kg1,neg,kimg) )

   t_wk1   = 0.0d0
   t_wk2   = 0.0d0

   do l = 1, kimg
   do j = 1, np_e_3D
        if ((sort == "sort")) then
             jto = neordr(neg_g(j),ik)
        else
             jto = neg_g(j)
        endif
   do i = 1,np_g1k_3D(ik)
        is = nis_g1k_3D(myrank_g_3D,ik) + i-1
        t_wk2(is,jto,l) = zaj_l_3D(i,j,ik,l)
   enddo; enddo; enddo

   nn = kg1*neg*kimg
   call mpi_allreduce(t_wk2,t_wk1,nn,mpi_real8,mpi_sum, mpi_k_world(myrank_k),err)

   do l = 1, kimg
   do j = 1,np_e
   do i = 1,iba(ik)
     js = nis_e(myrank_e) + j-1
     zaj_l(i,j,ik,l) = t_wk1(i,js,l)
   enddo; enddo; enddo

   deallocate( t_wk1, t_wk2)

   call tstatc0_end(id_sname)

  end subroutine decomp_zaj_l_r_3D_ik

!===============================================================================
  subroutine decomp_fsr_l_3D(fsr_l,fsr_l_3D,ik,nrvf_ordr,sort)

  integer :: ierr,err
  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

  real(kind=DP) :: fsr_l(np_e,nlmta,ista_k:iend_k)
  real(kind=DP) :: fsr_l_3D(np_e_3D,np_fs_3D,ista_k_3D:iend_k_3D)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, ik, is, js, nn, ito
  real(kind=DP), allocatable, dimension(:,:,:) :: t_wk1, t_wk2
  integer :: nrvf_ordr(neg,ista_k:iend_k)
  character(len=4) :: sort

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate( t_wk1(neg,nlmta,ista_k:iend_k) )
   allocate( t_wk2(neg,nlmta,ista_k:iend_k) )

   t_wk1   = 0.0d0
   t_wk2   = 0.0d0
   fsr_l_3D = 0.0d0

   do k = ista_k, iend_k
   do j = 1,nlmta
   do i = 1,np_e
     if ((sort == "sort") .and. (k == ik) )then
     ito = nrvf_ordr(nis_e(myrank_e) + i-1,ik)
     t_wk1(ito,j,k) = fsr_l(i,j,k)
   else
     is = nis_e(myrank_e) + i-1
     t_wk1(is,j,k) = fsr_l(i,j,k)
   endif
   enddo; enddo; enddo

   nn = neg*nlmta*(iend_k-ista_k+1)
   call mpi_allreduce(t_wk1,t_wk2,nn,mpi_real8,mpi_sum, mpi_k_world(myrank_k),err)

   do k = ista_k_3D, iend_k_3D
   do i = 1, np_e_3D
       is = neg_g(i)
   do j = 1,np_fs_3D
       js = nis_fs_3D(myrank_g_3D) + j-1
       fsr_l_3D(i,j,k) = t_wk2(is,js,k)
   enddo; enddo; enddo


   deallocate( t_wk1, t_wk2)

   call tstatc0_end(id_sname)

  end subroutine decomp_fsr_l_3D

!===============================================================================
  subroutine decomp_fsr_l_r_3D(fsr_l,fsr_l_3D,ik,neordr,sort,isub)

  integer :: ierr,err
  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

  real(kind=DP) :: fsr_l(np_e,nlmta,ista_k:iend_k)
  real(kind=DP) :: fsr_l_3D(np_e_3D,np_fs_3D,ista_k_3D:iend_k_3D)
  integer :: NB, nbs, lb, local_block, jj
  integer :: i, j, k, l, m, n, ik, is, js, nn, ifrom, isub, ito
  real(kind=DP), allocatable, dimension(:,:,:) :: t_wk1, t_wk2, t_wk3, t_wk4
  integer :: neordr(neg,ista_k:iend_k)
  character(len=4) :: sort
  logical :: existance
  character(len=3) :: c_num
  character(len=2) :: c_img
  character(len=80) :: fname

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)


   allocate( t_wk1(neg,nlmta,ista_k:iend_k) )
   allocate( t_wk2(neg,nlmta,ista_k:iend_k) )

   t_wk1   = 0.0d0
   t_wk2   = 0.0d0

   do k = ista_k_3D, iend_k_3D
   do j = 1,np_fs_3D
   do i=1,np_e_3D
        if ((sort == "sort") .and. (k == ik))then
             ito = neordr(neg_g(i),ik)
        else
             ito = neg_g(i)
        endif
           js = nis_fs_3D(myrank_g_3D) + j-1
        t_wk2(ito,js,k) = fsr_l_3D(i,j,k)
   enddo; enddo; enddo

   nn = neg*nlmta*(iend_k-ista_k+1)
   call mpi_allreduce(t_wk2,t_wk1,nn,mpi_real8,mpi_sum, mpi_k_world(myrank_k),err)

   do k = ista_k, iend_k
   do i = 1,np_e
   do j = 1,nlmta
     is = nis_e(myrank_e) + i-1
     fsr_l(i,j,k) = t_wk1(is,j,k)
   enddo; enddo; enddo

   deallocate( t_wk1, t_wk2)

   call tstatc0_end(id_sname)

  end subroutine decomp_fsr_l_r_3D

!===============================================================================
 subroutine decomp_fsr_l_2D(inn,out,ik)
  integer :: ierr,err
  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI
  real(kind=DP), dimension(:,:,:) ,intent(in) :: inn
  real(kind=DP), dimension(:,:)   ,intent(out):: out
  integer(kind=4)                 ,intent(in) :: ik
  real(kind=DP), allocatable,dimension(:,:) :: send,recv
  integer(kind=4) :: i,j

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(neg,nlmta),stat=ierr)
   allocate(recv(neg,nlmta),stat=ierr)
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send = 0.0d0
   recv = 0.0d0

   do j = ista_e, iend_e
      do i = 1,nlmta
         send(j,i) = inn(j-ista_e+1,i,ik)
      enddo
!!    send(j,1:nlmta) = inn(j-ista_e+1,1:nlmta,ik)
   enddo

   call mpi_allreduce(send,recv,neg*nlmta,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)

!  icount = 0
!  do j = 1, neg
!     if(map_e_3D(j) .eq. myrank_e_3D) then
!        icount = icount + 1
!        out(icount,1:nlmta) = recv(j,1:nlmta)
!     endif
!  enddo
   do j = 1,np_e_3d
      do i = 1,nlmta
         out(j,i) = recv(neg_g(j),i)
      enddo
!!    out(j,1:nlmta) = recv(neg_g(j),1:nlmta)
   enddo

   deallocate(send)
   deallocate(recv)

   call tstatc0_end(id_sname)

 end subroutine decomp_fsr_l_2D

!===============================================================================
 subroutine decomp_fsr_l_r_2D(inn,out,ik)
  integer :: ierr,err
  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI
  real(kind=DP), dimension(:,:)   ,intent(in) :: inn
  real(kind=DP), dimension(:,:,:) ,intent(out):: out
  integer(kind=4)                 ,intent(in) :: ik
  real(kind=DP), allocatable,dimension(:,:) :: send,recv
  integer(kind=4) :: i,j

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(neg,nlmta),stat=ierr)
   allocate(recv(neg,nlmta),stat=ierr)
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send = 0.0d0
   recv = 0.0d0

   do i = 1,nlmta
      do j = 1, np_e_3D
         send(neg_g(j),i) = inn(j,i)
      enddo
   enddo

   call mpi_allreduce(send,recv,neg*nlmta,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)

   do i = 1,nlmta
      do j = ista_e,iend_e
         out(j-ista_e+1,i,ik) = recv(j,i)
      enddo
   enddo

   deallocate(send)
   deallocate(recv)

   call tstatc0_end(id_sname)

 end subroutine decomp_fsr_l_r_2D

!===============================================================================
  subroutine decomp_eko_l_3D(inn,out,ik)
  integer :: ierr,err
  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI
  real(kind=DP), dimension(:,:) ,intent(in) :: inn
  real(kind=DP), dimension(:)   ,intent(out):: out
  integer(kind=4)               ,intent(in) :: ik
  real(kind=DP), allocatable,dimension(:)   :: send,recv
  integer(kind=4) :: i

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(neg),stat=ierr)
   allocate(recv(neg),stat=ierr)
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send = 0.0d0
   recv = 0.0d0

   do i = ista_e, iend_e
      send(i) = inn(i-ista_e+1,ik)
   enddo

   call mpi_allreduce(send,recv,neg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)

   do i = 1,np_e_3D
      out(i) = recv(neg_g(i))
   enddo

   deallocate(send)
   deallocate(recv)

   call tstatc0_end(id_sname)
  endsubroutine decomp_eko_l_3D

!===============================================================================
  subroutine decomp_eko_l_3D_2(inn,out,ik,nrvf_ordr,sort)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(:,:) ,intent(in) :: inn
  real(kind=DP), dimension(:,:) ,intent(out):: out
  integer(kind=4)               ,intent(in) :: ik
  integer(kind=4),dimension(neg,ista_k_3D:iend_k_3D),intent(in) :: nrvf_ordr
  character(len=4) ,intent(in) :: sort
  real(kind=DP), allocatable,dimension(:) :: send, recv
  integer(kind=4) :: i, iadd, ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(neg),stat=ierr)
   allocate(recv(neg),stat=ierr)
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send = 0.0d0
   recv = 0.0d0

   do i = ista_e, iend_e
      if ((sort == "sort"))then
         iadd = nrvf_ordr(i,ik)
      else
         iadd = i
      endif
      send(iadd) = inn(i-ista_e+1,ik)
   enddo

   call mpi_allreduce(send,recv,neg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)

   do i = 1,np_e_3D
      out(i,ik) = recv(neg_g(i))
   enddo

   deallocate(send)
   deallocate(recv)

   call tstatc0_end(id_sname)

  endsubroutine decomp_eko_l_3D_2

!===============================================================================
  subroutine decomp_eko_l_r_3D(eko_l,eko_3D,ik)
  implicit none
  include 'mpif.h'
  real(kind=DP), intent(out), dimension(np_e,   ista_k:iend_k) :: eko_l
  real(kind=DP), intent(in),  dimension(np_e_3D)               :: eko_3D
  real(kind=DP),              dimension(neg)                   :: wk1
  integer, intent(in) :: ik
  integer :: i, ncnt, ierr

  call tstatc0_begin('z_interface_3D__________________',id_sname,1)

  ncnt = neg
  wk1 = 0.0d0

  do i = 1, np_e_3D
     wk1(neg_g(i)) = eko_3D(i)
  enddo

  call mpi_allreduce(MPI_IN_PLACE,wk1,ncnt,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_kg_world,ierr)

  do i = 1, np_e
     eko_l(i,ik) = wk1(i+ista_e-1)
  enddo

   call tstatc0_end(id_sname)
  endsubroutine decomp_eko_l_r_3D

!===============================================================================
  subroutine decomp_eko_l_r_3D_2(out,inn,ik,ordr,sort)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(:,:) ,intent(in) :: inn
  real(kind=DP), dimension(:,:) ,intent(out):: out
  integer(kind=4)               ,intent(in) :: ik
  integer(kind=4),dimension(neg,ista_k_3D:iend_k_3D),intent(in) :: ordr
  character(len=4), intent(in) :: sort
  real(kind=DP), allocatable,dimension(:) :: send, recv
  integer(kind=4) :: i, iadd, ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(neg),stat=ierr)
   allocate(recv(neg),stat=ierr)
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send = 0.0d0
   recv = 0.0d0
   if (myrank_g_3D == 0) then
      do i = 1, np_e_3D
         if (sort == "sort") then
            iadd = ordr(neg_g(i),ik)
         else
            iadd = neg_g(i)
         endif
         send(iadd) = inn(i,ik)
      enddo
   endif

   call mpi_allreduce(send,recv,neg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)

   iadd = 0
   do i = ista_e, iend_e
      iadd = iadd + 1
      out(iadd,ik) = recv(i)
   enddo

   deallocate(send)
   deallocate(recv)

   call tstatc0_end(id_sname)

  endsubroutine decomp_eko_l_r_3D_2

!===============================================================================
  subroutine decomp_vnlph_l_3D(inn,out,ik,kimg,nrvf_ordr)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(:,:,:) ,intent(in)  :: inn
  real(kind=DP), dimension(:,:,:) ,intent(out) :: out
  integer(kind=4)                 ,intent(in) :: ik,kimg
  integer(kind=4),dimension(neg,ista_k_3D:iend_k_3D),intent(in) :: nrvf_ordr
  real(kind=DP), allocatable,dimension(:,:) :: send,recv
  integer(kind=4) :: i,j,index1,index2,ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(iba(ik),neg),stat=ierr)
   allocate(recv(iba(ik),neg),stat=ierr)
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send = 0.0d0
   recv = 0.0d0

   do j = ista_e, iend_e
      do i = 1, iba(ik)
         send(i,nrvf_ordr(j,ik)) = inn(i,j-ista_e+1,kimg)
      enddo
   enddo

   call mpi_allreduce(send,recv,iba(ik)*neg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)

   do j = 1, np_e_3D
      index1= neg_g(j)      !   local -> global
      index2 = 0
      do i = ista_g1k_3D(ik),iend_g1k_3D(ik)
         index2 = index2 + 1
         out(index2,j,kimg) = recv(i,index1)
      enddo
   enddo

   deallocate(send)
   deallocate(recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_vnlph_l_3D

!===============================================================================
  subroutine decomp_vnlph_l_r_3D(inn,out,ik,kimg)
  integer :: ierr,err
  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI
  real(kind=DP), dimension(:,:,:) ,intent(in) :: inn
  real(kind=DP), dimension(:,:,:) ,intent(out):: out
  integer(kind=4)                 ,intent(in) :: ik,kimg
  real(kind=DP), allocatable,dimension(:,:) :: send,recv
  integer(kind=4) :: i,j,index1,index2
  logical :: existance
  character(len=3) :: c_num
  character(len=2) :: c_img
  character(len=80) :: fname

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(iba(ik),neg),stat=ierr)
   allocate(recv(iba(ik),neg),stat=ierr)
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send = 0.0d0
   recv = 0.0d0

   do j = 1, np_e_3D
      index1= neg_g(j)      !   local -> global
      index2 = 0
      do i = ista_g1k_3D(ik),iend_g1k_3D(ik)
         index2 = index2 + 1
         send(i,index1) = inn(index2,j,kimg)
      enddo
   enddo

   call mpi_allreduce(send,recv,iba(ik)*neg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)

   index1 = 0
   do j = ista_e, iend_e
      index1 = index1 + 1
      do i = 1, iba(ik)
         out(i,index1,kimg) = recv(i,j)
      enddo
   enddo

   deallocate(send)
   deallocate(recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_vnlph_l_r_3D

!===============================================================================
  subroutine decomp_vnlph_l_r_3D_2(inn,out,ik,kimg,nrvf_ordr,sort)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(:,:,:) ,intent(in) :: inn
  real(kind=DP), dimension(:,:,:) ,intent(out):: out
  integer(kind=4)                 ,intent(in) :: ik,kimg
  character(len=4)                ,intent(in) :: sort
  integer(kind=4),dimension(neg,ista_k_3D:iend_k_3D),intent(in) :: nrvf_ordr
  real(kind=DP), allocatable,dimension(:,:) :: send,recv
  integer(kind=4) :: i,j,index1,index2,ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(iba(ik),neg),stat=ierr)
   allocate(recv(iba(ik),neg),stat=ierr)
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send = 0.0d0
   recv = 0.0d0

   do j = 1, np_e_3D
!     index1= neg_g(j)      !   local -> global
      if (sort == "sort") then
         index1 = nrvf_ordr(neg_g(j),ik)
      else
         index1 = neg_g(j)
      endif
      index2 = 0
      do i = ista_g1k_3D(ik),iend_g1k_3D(ik)
         index2 = index2 + 1
         send(i,index1) = inn(index2,j,kimg)
      enddo
   enddo

   call mpi_allreduce(send,recv,iba(ik)*neg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)

   index1 = 0
   do j = ista_e, iend_e
      index1 = index1 + 1
      do i = 1, iba(ik)
         out(i,index1,kimg) = recv(i,j)
      enddo
   enddo

   deallocate(send)
   deallocate(recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_vnlph_l_r_3D_2

!===============================================================================
  subroutine decomp_snl_l_3D(inn,out,ik,nlmtt,iksnl)
  real(kind=DP), dimension(:,:,:) ,intent(in) :: inn
  real(kind=DP), dimension(:,:,:) ,intent(out):: out
  integer(kind=4)               ,intent(in) :: ik, nlmtt, iksnl
  integer(kind=4) :: i,j,index1

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   B_NLMTT : do i = 1,nlmtt
      index1 = 0
      B_BAND : do j = ista_g1k_3D(ik), iend_g1k_3D(ik)
        index1 = index1 + 1
!xx     out(index1,i,ik) = inn(j,i,ik)
        out(index1,i,iksnl) = inn(j,i,iksnl)
      enddo B_BAND
   enddo B_NLMTT

   call tstatc0_end(id_sname)

  end subroutine decomp_snl_l_3D

!===============================================================================
  subroutine decomp_snl_l_3D_2(inn,out,ik)
  real(kind=DP), dimension(kg1,nlmtt,ista_snl:iend_snl) ,intent(in) :: inn
  real(kind=DP), dimension(maxval(np_g1k_3D),nlmtt,ista_snl:iend_snl) ,intent(out):: out

  integer(kind=4) :: ik, iksnl
  integer(kind=4) :: i,j,index1

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

  out = 0.0d0
   K       : do iksnl = ista_snl,iend_snl
   B_NLMTT : do i = 1,nlmtt
      index1 = 0
      B_BAND : do j = ista_g1k_3D(iksnl), iend_g1k_3D(iksnl)
        index1 = index1 + 1
        out(index1,i,iksnl) = inn(j,i,iksnl)
      enddo B_BAND
   enddo B_NLMTT
   enddo K

   call tstatc0_end(id_sname)

  end subroutine decomp_snl_l_3D_2

!===============================================================================
  subroutine decomp_snl_l_r_3D(inn,out)
  real(kind=DP), dimension(:,:,:) ,intent(in) :: inn
  real(kind=DP), dimension(:,:,:) ,intent(out):: out
  real(kind=DP), allocatable,dimension(:,:,:) :: send,recv
  integer(kind=4) :: i,j,k,nn,ierr,err
  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(kg1,nlmtt,ista_snl:iend_snl),stat=ierr)
   allocate(recv(kg1,nlmtt,ista_snl:iend_snl),stat=ierr)
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send = 0.0d0
   recv = 0.0d0

   do k = ista_snl,iend_snl
    do j = 1,nlmtt
      do i = 1, np_g1k_3D(k)
        send(ista_g1k_3D(k)+i-1,j,k) = inn(i,j,k)
      enddo
    enddo
   enddo


   nn = kg1*nlmtt*(iend_snl-ista_snl+1)
   call mpi_allreduce(send,recv,nn,mpi_real8,mpi_sum, mpi_ke_world,err)

   do k = ista_snl,iend_snl
    do j = 1,nlmtt
      do i = 1, kg1
        out(i,j,k) = recv(i,j,k)
      enddo
    enddo
   enddo

   deallocate(send)
   deallocate(recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_snl_l_r_3D

!===============================================================================
  subroutine decomp_fsr_l_3D_ik(fsr_l,fsr_l_3D,ik,nrvf_ordr,sort)

  integer :: ierr,err
  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

  real(kind=DP) :: fsr_l(np_e,nlmta,ista_k:iend_k)
  real(kind=DP) :: fsr_l_3D(np_e_3D,np_fs_3D,ista_k_3D:iend_k_3D)
  integer :: NB, nbs, lb
  integer :: i, j, k, l, m, n, ik, is, js, nn, ito
  real(kind=DP), allocatable, dimension(:,:) :: t_wk1, t_wk2
  integer :: nrvf_ordr(neg,ista_k:iend_k)
  character(len=4) :: sort

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate( t_wk1(neg,nlmta) )
   allocate( t_wk2(neg,nlmta) )

   t_wk1   = 0.0d0
   t_wk2   = 0.0d0
!  fsr_l_3D = 0.0d0

   do j = 1,nlmta
   do i = 1,np_e
   if ((sort == "sort"))then
     ito = nrvf_ordr(nis_e(myrank_e) + i-1,ik)
     t_wk1(ito,j) = fsr_l(i,j,ik)
   else
     is = nis_e(myrank_e) + i-1
     t_wk1(is,j) = fsr_l(i,j,ik)
   endif
   enddo; enddo

   nn = neg*nlmta
   call mpi_allreduce(t_wk1,t_wk2,nn,mpi_real8,mpi_sum, mpi_k_world(myrank_k),err)

   do i = 1, np_e_3D
       is = neg_g(i)
   do j = 1,np_fs_3D
       js = nis_fs_3D(myrank_g_3D) + j-1
       fsr_l_3D(i,j,ik) = t_wk2(is,js)
   enddo; enddo


   deallocate( t_wk1, t_wk2)

   call tstatc0_end(id_sname)

  end subroutine decomp_fsr_l_3D_ik

!===============================================================================
  subroutine decomp_fsr_l_r_3D_ik(fsr_l,fsr_l_3D,ik,neordr,sort,isub)

  integer :: ierr,err
  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

  real(kind=DP) :: fsr_l(np_e,nlmta,ista_k:iend_k)
  real(kind=DP) :: fsr_l_3D(np_e_3D,np_fs_3D,ista_k_3D:iend_k_3D)
  integer :: NB, nbs, lb, local_block, jj
  integer :: i, j, ik, is, js, nn, isub, ito
  real(kind=DP), allocatable, dimension(:,:) :: t_wk1, t_wk2
  integer :: neordr(neg,ista_k:iend_k)
  character(len=4) :: sort

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate( t_wk1(neg,nlmta) )
   allocate( t_wk2(neg,nlmta) )

   t_wk1   = 0.0d0
   t_wk2   = 0.0d0

   do j = 1,np_fs_3D
   do i=1,np_e_3D
        if ((sort == "sort"))then
             ito = neordr(neg_g(i),ik)
        else
             ito = neg_g(i)
        endif
        js = nis_fs_3D(myrank_g_3D) + j-1
        t_wk2(ito,js) = fsr_l_3D(i,j,ik)
   enddo; enddo

   nn = neg*nlmta
   call mpi_allreduce(t_wk2,t_wk1,nn,mpi_real8,mpi_sum, mpi_k_world(myrank_k),err)

   do i = 1,np_e
   do j = 1,nlmta
     is = nis_e(myrank_e) + i-1
     fsr_l(i,j,ik) = t_wk1(is,j)
   enddo; enddo

   deallocate( t_wk1, t_wk2)

   call tstatc0_end(id_sname)

  end subroutine decomp_fsr_l_r_3D_ik

!===============================================================================
  subroutine replacement_zaj_ball_sequence(zaj_ball,ik,neordr,nblocksize_mgs_default)

   use m_Control_Parameters,  only : nblocksize_mgs, nblocksize_mgs_is_given

   real(kind=8),intent(inout),dimension(maxval(np_g1k_3D),neg,ista_k_3D:iend_k_3D,kimg)::zaj_ball
   integer(kind=4),intent(in)::ik, nblocksize_mgs_default
   integer(kind=4),intent(in),dimension(neg,ista_k_3D:iend_k_3D)::neordr

   real(kind=8),dimension(np_g1k_3D(ik),neg,kimg) :: work
   integer(kind=4) :: NB, nbs_num, nbsn_num, ival, icnt, i, j, k
!! integer(kind=4),dimension(neg) :: neg_g_all

    if(nblocksize_mgs_is_given) then
      NB = nblocksize_mgs
    else
      NB = nblocksize_mgs_default
    end if

    nbs_num = (neg-1) / NB + 1
    nbsn_num = (nbs_num - 1) / nrank_e_3D + 1

!!  icnt = 0
!!  do i = 1, nrank_e_3D
!!    do j = 1, nbsn_num
!!      do k = 1, NB
!!        ival = nbsn_num*NB*(j-1) + NB*(i-1) + k
!!        ival = nrank_e_3D*NB*(j-1) + NB*(i-1) + k
!!        if (ival > neg)  cycle
!!        icnt = icnt + 1
!!        if (icnt > neg)  exit
!!        neg_g_all(icnt) = ival
!!      enddo
!!    enddo
!!  enddo
!!
!!  order of eigenvalue block cyclic  ->  order of eigenvalue
!!
    do k = 1, kimg
      do j = 1, neg
        do i = 1, np_g1k_3D(ik)
          work(i,neg_g_all(j),k) = zaj_ball(i,j,ik,k)
        enddo
      enddo
    enddo
!!
!!  order of eigenvalue  ->  order of sequence
!!
    do k = 1, kimg
      do j = 1, neg
        do i = 1, np_g1k_3D(ik)
          zaj_ball(i,neordr(j,ik),ik,k) = work(i,j,k)
        enddo
      enddo
    enddo

  end subroutine replacement_zaj_ball_sequence


!===============================================================================
  subroutine replacement_zaj_ball_eigenvalue(zaj_ball,ik,nrvf_ordr)

   real(kind=8),intent(inout),dimension(maxval(np_g1k_3D),neg,ista_k_3D:iend_k_3D,kimg)::zaj_ball
   integer(kind=4),intent(in)::ik
   integer(kind=4),intent(in),dimension(neg,ista_k_3D:iend_k_3D)::nrvf_ordr

   real(kind=8),dimension(np_g1k_3D(ik),neg,kimg) :: work
   integer(kind=4) :: i, j, k

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

!!
!!  order of sequence  ->  order of eigenvalue
!!
    do k = 1, kimg
      do j = 1, neg
        do i = 1, np_g1k_3D(ik)
          work(i,nrvf_ordr(j,ik),k) = zaj_ball(i,j,ik,k)
        enddo
      enddo
    enddo
!!
!!  order of eigenvalue  ->  order of eigenvalue block cyclic
!!
    do k = 1, kimg
      do j = 1, neg
        do i = 1, np_g1k_3D(ik)
!         zaj_ball(i,j,ik,k) = work(i,neg_gg_all(j),k)
          zaj_ball(i,j,ik,k) = work(i,neg_g_all(j),k)
        enddo
      enddo
    enddo

   call tstatc0_end(id_sname)

  end subroutine replacement_zaj_ball_eigenvalue

!===============================================================================
  subroutine decomp_occup_l_3D(inn,out,ik,nrvf_ordr,sort)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(:,:) ,intent(in) :: inn
  real(kind=DP), dimension(:,:) ,intent(out):: out
  integer(kind=4)               ,intent(in) :: ik
  integer(kind=4),dimension(neg,ista_k_3D:iend_k_3D),intent(in) :: nrvf_ordr
  character(len=4) ,intent(in) :: sort
  real(kind=DP), allocatable,dimension(:) :: send, recv
  integer(kind=4) :: i, iadd, ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(neg),stat=ierr)
   allocate(recv(neg),stat=ierr)
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send = 0.0d0
   recv = 0.0d0

   do i = ista_e, iend_e
      if ((sort == "sort"))then
         iadd = nrvf_ordr(i,ik)
      else
         iadd = i
      endif
      send(iadd) = inn(i-ista_e+1,ik)
   enddo

   call mpi_allreduce(send,recv,neg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)

   do i = 1,np_e_3D
      out(i,ik) = recv(neg_g(i))
   enddo

   deallocate(send)
   deallocate(recv)

   call tstatc0_end(id_sname)

  endsubroutine decomp_occup_l_3D

!===============================================================================
  subroutine decomp_occup_l_r_3D(out,inn,ik,ordr,sort)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(:,:) ,intent(in) :: inn
  real(kind=DP), dimension(:,:) ,intent(out):: out
  integer(kind=4)               ,intent(in) :: ik
  integer(kind=4),dimension(neg,ista_k_3D:iend_k_3D),intent(in) :: ordr
  character(len=4), intent(in) :: sort
  real(kind=DP), allocatable,dimension(:) :: send, recv
  integer(kind=4) :: i, iadd, ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(neg),stat=ierr)
   allocate(recv(neg),stat=ierr)
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send = 0.0d0
   recv = 0.0d0
   if (myrank_g_3D == 0) then
      do i = 1, np_e_3D
         if (sort == "sort") then
            iadd = ordr(neg_g(i),ik)
         else
            iadd = neg_g(i)
         endif
         send(iadd) = inn(i,ik)
      enddo
   endif

   call mpi_allreduce(send,recv,neg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)

   iadd = 0
   do i = ista_e, iend_e
      iadd = iadd + 1
      out(iadd,ik) = recv(i)
   enddo

   deallocate(send)
   deallocate(recv)

   call tstatc0_end(id_sname)

  endsubroutine decomp_occup_l_r_3D

!===============================================================================
  subroutine decomp_vlhxc_l_3D(vlhxc_l,vlhxc_l_3D,ispin)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(ista_kngp:iend_kngp,kimg,nspin) ,intent(in) :: vlhxc_l
  real(kind=DP), dimension(ista_kngp_3D:iend_kngp_3D,kimg,nspin) ,intent(out):: vlhxc_l_3D
  integer(kind=4)               ,intent(in) :: ispin
  real(kind=DP), allocatable,dimension(:,:) :: send, recv
  integer(kind=4) :: i, ri, ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(1:kgp,kimg), stat=ierr)
   allocate(recv(1:kgp,kimg), stat=ierr)
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send(:,:) = 0.0d0
   do ri = 1, kimg
      do i = ista_kngp, iend_kngp
         if (i > kgp) exit
         send(i,ri) = vlhxc_l(i,ri,ispin)
      end do
   end do

   call mpi_allreduce(send,recv,kgp*kimg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)

   do ri = 1, kimg
      do i = ista_kngp_3D, iend_kngp_3D
         if (i > kgp) exit
         vlhxc_l_3D(i,ri,ispin) = recv(i,ri)
      end do
   end do

   deallocate(send,recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_vlhxc_l_3D

!===============================================================================
  subroutine decomp_vlhxc_l_r_3D(vlhxc_l,vlhxc_l_3D,ispin)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(ista_kngp:iend_kngp,kimg,nspin) ,intent(out) :: vlhxc_l
  real(kind=DP), dimension(ista_kngp_3D:iend_kngp_3D,kimg,nspin) ,intent(in):: vlhxc_l_3D
  integer(kind=4)               ,intent(in) :: ispin
  real(kind=DP), allocatable,dimension(:,:) :: send, recv
  integer(kind=4) :: i, j, ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(1:kgp,kimg))
   allocate(recv(1:kgp,kimg))
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send(:,:) = 0.0d0
!  if (myrank_e_3D == 0) then
      do j = 1, kimg
         do i = ista_kngp_3D, iend_kngp_3D
            if (i > kgp) exit
            send(i,j) = vlhxc_l_3D(i,j,ispin)
         end do
      end do
!  endif

!  call mpi_allreduce(send,recv,kgp*kimg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)
   call mpi_allreduce(send,recv,kgp*kimg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)

   vlhxc_l(:,:,ispin) = 0.0d0

   do j = 1, kimg
      do i = ista_kngp, iend_kngp
         if (i > kgp) exit
         vlhxc_l(i,j,ispin) = recv(i,j)
      end do
   end do

   deallocate(send,recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_vlhxc_l_r_3D

!===============================================================================
  subroutine decomp_ylm_l_3D(ylm_l,ylm_l_3D,nel_Ylm)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(ista_kngp:iend_kngp,nel_Ylm) ,intent(in) :: ylm_l
  real(kind=DP), dimension(ista_kngp_3D:iend_kngp_3D,nel_Ylm) ,intent(out):: ylm_l_3D
  integer(kind=4)               ,intent(in) :: nel_Ylm
  real(kind=DP), allocatable,dimension(:,:) :: send, recv
  integer(kind=4) :: i, j, ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(1:kgp,nel_Ylm), stat=ierr)
   allocate(recv(1:kgp,nel_Ylm), stat=ierr)
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send(:,:) = 0.0d0
   do j = 1, nel_Ylm
      do i = ista_kngp, iend_kngp
         if (i > kgp) exit
         send(i,j) = ylm_l(i,j)
      end do
   end do

   call mpi_allreduce(send,recv,kgp*nel_Ylm,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)

   do j = 1, nel_Ylm
      do i = ista_kngp_3D, iend_kngp_3D
         if (i > kgp) exit
         ylm_l_3D(i,j) = recv(i,j)
      end do
   end do

   deallocate(send,recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_ylm_l_3D

!===============================================================================
  subroutine decomp_qitg_l_3D(qitg_l,qitg_l_3D,nqitg)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(ista_kngp:iend_kngp,nqitg) ,intent(in) :: qitg_l
  real(kind=DP), dimension(ista_kngp_3D:iend_kngp_3D,nqitg) ,intent(out):: qitg_l_3D
  integer(kind=4)               ,intent(in) :: nqitg
  real(kind=DP), allocatable,dimension(:,:) :: send, recv
  integer(kind=4) :: i, j, ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   if (nqitg == 0) return

   allocate(send(1:kgp,nqitg))
   allocate(recv(1:kgp,nqitg))
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send(:,:) = 0.0d0
   do j = 1, nqitg
      do i = ista_kngp, iend_kngp
         if (i > kgp) exit
         send(i,j) = qitg_l(i,j)
      end do
   end do

   call mpi_allreduce(send,recv,kgp*nqitg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)

   do j = 1, nqitg
      do i = ista_kngp_3D, iend_kngp_3D
         if (i > kgp) exit
         qitg_l_3D(i,j) = recv(i,j)
      end do
   end do

   deallocate(send,recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_qitg_l_3D

!===============================================================================
  subroutine decomp_qitg_l_r_3D(qitg_l,qitg_l_3D,nqitg)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(ista_kngp:iend_kngp,nqitg)       ,intent(out) :: qitg_l
  real(kind=DP), dimension(ista_kngp_3D:iend_kngp_3D,nqitg) ,intent(in)  :: qitg_l_3D
  integer(kind=4)               ,intent(in) :: nqitg
  real(kind=DP), allocatable,dimension(:,:) :: send, recv
  integer(kind=4) :: i, j, ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(1:kgp,nqitg))
   allocate(recv(1:kgp,nqitg))
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send(:,:) = 0.0d0
   do j = 1, nqitg
      do i = ista_kngp_3D, iend_kngp_3D
         if (i > kgp) exit
         send(i,j) = qitg_l_3D(i,j)
      end do
   end do

   call mpi_allreduce(send,recv,kgp*nqitg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)

   qitg_l(:,:) = 0.0d0

   do j = 1, nqitg
      do i = ista_kngp, iend_kngp
         if (i > kgp) exit
         qitg_l(i,j) = recv(i,j)
      end do
   end do

   deallocate(send,recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_qitg_l_r_3D

!===============================================================================
  subroutine decomp_zfm3_l_3D(zfm3_l,zfm3_l_3D)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(ista_kngp:iend_kngp,ntyp,kimg) ,intent(in) :: zfm3_l
  real(kind=DP), dimension(ista_kngp_3D:iend_kngp_3D,ntyp,kimg) ,intent(out):: zfm3_l_3D
  real(kind=DP), allocatable,dimension(:,:,:) :: send, recv
  integer(kind=4) :: i, j, ri, ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(1:kgp,ntyp,kimg))
   allocate(recv(1:kgp,ntyp,kimg))
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send(:,:,:) = 0.0d0
   do ri = 1, kimg
      do j = 1, ntyp
         do i = ista_kngp, iend_kngp
            if (i > kgp) exit
            send(i,j,ri) = zfm3_l(i,j,ri)
         end do
      end do
   end do

   call mpi_allreduce(send,recv,kgp*ntyp*kimg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)

   do ri = 1, kimg
      do j = 1, ntyp
         do i = ista_kngp_3D, iend_kngp_3D
            if (i > kgp) exit
            zfm3_l_3D(i,j,ri) = recv(i,j,ri)
         end do
      end do
   end do

   deallocate(send,recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_zfm3_l_3D

!===============================================================================
  subroutine decomp_zfm3_l_r_3D(zfm3_l,zfm3_l_3D)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(ista_kngp:iend_kngp,ntyp,kimg) ,intent(out) :: zfm3_l
  real(kind=DP), dimension(ista_kngp_3D:iend_kngp_3D,ntyp,kimg) ,intent(in):: zfm3_l_3D
  real(kind=DP), allocatable,dimension(:,:,:) :: send, recv
  integer(kind=4) :: i, j, ri, ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(1:kgp,ntyp,kimg))
   allocate(recv(1:kgp,ntyp,kimg))
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send(:,:,:) = 0.0d0
   do ri = 1, kimg
      do j = 1, ntyp
         do i = ista_kngp_3D, iend_kngp_3D
            if (i > kgp) exit
            send(i,j,ri) = zfm3_l_3D(i,j,ri)
         end do
      end do
   end do

   call mpi_allreduce(send,recv,kgp*ntyp*kimg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)

   zfm3_l(:,:,:) = 0.0d0

   do ri = 1, kimg
      do j = 1, ntyp
         do i = ista_kngp, iend_kngp
            if (i > kgp) exit
            zfm3_l(i,j,ri) = recv(i,j,ri)
         end do
      end do
   end do

   deallocate(send,recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_zfm3_l_r_3D

!===============================================================================
  subroutine decomp_psc_l_3D(psc_l,psc_l_3D)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(ista_kngp:iend_kngp,ntyp) ,intent(in) :: psc_l
  real(kind=DP), dimension(ista_kngp_3D:iend_kngp_3D,ntyp) ,intent(out):: psc_l_3D
  real(kind=DP), allocatable,dimension(:,:) :: send, recv
  integer(kind=4) :: i, j, ri, ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(1:kgp,ntyp))
   allocate(recv(1:kgp,ntyp))
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send(:,:) = 0.0d0
   do j = 1, ntyp
      do i = ista_kngp, iend_kngp
         if (i > kgp) exit
         send(i,j) = psc_l(i,j)
      end do
   end do

   call mpi_allreduce(send,recv,kgp*ntyp,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)

   do j = 1, ntyp
      do i = ista_kngp_3D, iend_kngp_3D
         if (i > kgp) exit
         psc_l_3D(i,j) = recv(i,j)
      end do
   end do

   deallocate(send,recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_psc_l_3D

!===============================================================================
  subroutine decomp_psc_l_r_3D(psc_l,psc_l_3D)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(ista_kngp:iend_kngp,ntyp) ,intent(out) :: psc_l
  real(kind=DP), dimension(ista_kngp_3D:iend_kngp_3D,ntyp) ,intent(in):: psc_l_3D
  real(kind=DP), allocatable,dimension(:,:) :: send, recv
  integer(kind=4) :: i, j, ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(1:kgp,ntyp))
   allocate(recv(1:kgp,ntyp))
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send(:,:) = 0.0d0
   do j = 1, ntyp
      do i = ista_kngp_3D, iend_kngp_3D
         if (i > kgp) exit
         send(i,j) = psc_l_3D(i,j)
      end do
   end do

   call mpi_allreduce(send,recv,kgp*ntyp,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)

   psc_l(:,:) = 0.0d0

   do j = 1, ntyp
      do i = ista_kngp, iend_kngp
         if (i > kgp) exit
         psc_l(i,j) = recv(i,j)
      end do
   end do

   deallocate(send,recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_psc_l_r_3D
!===============================================================================

  subroutine decomp_vdip_l_3D(vdip_l,vdip_l_3D)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(ista_kngp:iend_kngp,kimg) ,intent(in) :: vdip_l
  real(kind=DP), dimension(ista_kngp_3D:iend_kngp_3D,kimg) ,intent(out):: vdip_l_3D
  real(kind=DP), allocatable,dimension(:,:) :: send, recv
  integer(kind=4) :: i, ri, ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(1:kgp,kimg))
   allocate(recv(1:kgp,kimg))
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send(:,:) = 0.0d0
   do ri = 1, kimg
      do i = ista_kngp, iend_kngp
         if (i > kgp) exit
         send(i,ri) = vdip_l(i,ri)
      end do
   end do

   call mpi_allreduce(send,recv,kgp*kimg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)

   do ri = 1, kimg
      do i = ista_kngp_3D, iend_kngp_3D
         if (i > kgp) exit
         vdip_l_3D(i,ri) = recv(i,ri)
      end do
   end do

   deallocate(send,recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_vdip_l_3D

!===============================================================================
  subroutine decomp_vdip_l_r_3D(vdip_l,vdip_l_3D)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(ista_kngp:iend_kngp,kimg) ,intent(out) :: vdip_l
  real(kind=DP), dimension(ista_kngp_3D:iend_kngp_3D,kimg) ,intent(in):: vdip_l_3D
  real(kind=DP), allocatable,dimension(:,:) :: send, recv
  integer(kind=4) :: i, ri, ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(1:kgp,kimg))
   allocate(recv(1:kgp,kimg))
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send(:,:) = 0.0d0
   do ri = 1, kimg
      do i = ista_kngp_3D, iend_kngp_3D
         if (i > kgp) exit
         send(i,ri) = vdip_l_3D(i,ri)
      end do
   end do

   call mpi_allreduce(send,recv,kgp*kimg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)

   vdip_l(:,:) = 0.0d0

   do ri = 1, kimg
      do i = ista_kngp, iend_kngp
         if (i > kgp) exit
         vdip_l(i,ri) = recv(i,ri)
      end do
   end do

   deallocate(send,recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_vdip_l_r_3D

!===============================================================================
  subroutine decomp_gr_l_3D(gr_l,gr_l_3D)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(ista_kngp:iend_kngp) ,intent(in) :: gr_l
  real(kind=DP), dimension(ista_kngp_3D:iend_kngp_3D) ,intent(out):: gr_l_3D
  real(kind=DP), allocatable,dimension(:) :: send, recv
  integer(kind=4) :: i, ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(1:kgp))
   allocate(recv(1:kgp))
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send(:) = 0.0d0
   do i = ista_kngp, iend_kngp
      if (i > kgp) exit
      send(i) = gr_l(i)
   end do

   call mpi_allreduce(send,recv,kgp,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)

   do i = ista_kngp_3D, iend_kngp_3D
      if (i > kgp) exit
      gr_l_3D(i) = recv(i)
   end do

   deallocate(send,recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_gr_l_3D

!===============================================================================
  subroutine decomp_gr_l_r_3D(gr_l,gr_l_3D)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(ista_kngp:iend_kngp) ,intent(out) :: gr_l
  real(kind=DP), dimension(ista_kngp_3D:iend_kngp_3D) ,intent(in):: gr_l_3D
  real(kind=DP), allocatable,dimension(:) :: send, recv
  integer(kind=4) :: i, ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(1:kgp))
   allocate(recv(1:kgp))
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send(:) = 0.0d0
   do i = ista_kngp_3D, iend_kngp_3D
      if (i > kgp) exit
      send(i) = gr_l_3D(i)
   end do

   call mpi_allreduce(send,recv,kgp,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)

   gr_l(:) = 0.0d0

   do i = ista_kngp, iend_kngp
      if (i > kgp) exit
      gr_l(i) = recv(i)
   end do

   deallocate(send,recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_gr_l_r_3D

!===============================================================================
  subroutine decomp_rr_is_over_or_under_3D(rr_is_over_or_under, rr_is_over_or_under_3D, ik, nrvf_ordr, sort)
  implicit none
  include 'mpif.h'
  integer, intent(in)  :: rr_is_over_or_under   (np_e,   ista_k:iend_k)
  integer, intent(out) :: rr_is_over_or_under_3D(np_e_3D,ista_k_3D:iend_k)
  integer              :: wk1                   (neg)
  integer, intent(in) :: ik
  integer, intent(in) :: nrvf_ordr(neg, ista_k:iend_k)
  character(len=4), intent(in) :: sort
  integer :: j, l, js, jto, ncnt
  integer :: ierr

  ncnt = neg
  wk1 = 0

  if (sort == "sort") then
     do j = 1, np_e
        js = nrvf_ordr(nis_e(myrank_e) + j - 1, ik)
        wk1(js) = rr_is_over_or_under(j,ik)
     enddo
  else
     do j = 1, np_e
        js = nis_e(myrank_e) + j - 1
        wk1(js) = rr_is_over_or_under(j,ik)
     enddo
  endif

  call mpi_allreduce(MPI_IN_PLACE,wk1,ncnt,MPI_INTEGER,MPI_SUM,MPI_CommGroup,ierr)

  do j = 1, np_e_3D
     jto = neg_g(j)
     rr_is_over_or_under_3D(j,ik) = wk1(jto)
  enddo

  end subroutine decomp_rr_is_over_or_under_3D

!===============================================================================
  subroutine decomp_rr_is_over_or_under_r_3D(rr_is_over_or_under, rr_is_over_or_under_3D, ik, neordr, sort)
  implicit none
  include 'mpif.h'
  integer, intent(out) :: rr_is_over_or_under   (np_e,   ista_k:iend_k)
  integer, intent(in)  :: rr_is_over_or_under_3D(np_e_3D,ista_k_3D:iend_k)
  integer              :: wk1                   (neg)
  integer, intent(in) :: ik
  integer, intent(in) :: neordr(neg, ista_k:iend_k)
  character(len=4), intent(in) :: sort
  integer :: j, l, js, jto, ncnt
  integer :: ierr

  ncnt = neg
  wk1 = 0

  if (sort == "sort") then
     do j = 1, np_e_3D
        jto =  neordr(neg_g(j), ik)
        wk1(jto) = rr_is_over_or_under_3D(j,ik)
     enddo
  else
     do j = 1, np_e_3D
        jto = neg_g(j)
        wk1(jto) = rr_is_over_or_under_3D(j,ik)
     enddo
  endif

  call mpi_allreduce(MPI_IN_PLACE,wk1,ncnt,MPI_INTEGER,MPI_SUM,mpi_kg_world,ierr)

  do j = 1, np_e
     js = nis_e(myrank_e) + j - 1
     rr_is_over_or_under(j,ik) = wk1(js)
  enddo

  end subroutine decomp_rr_is_over_or_under_r_3D

!===============================================================================
  subroutine decomp_wfsd_l_3D(wfsd_l, wfsd_l_3D, ik)
  implicit none
  include 'mpif.h'
  real(kind=DP), intent(in)  :: wfsd_l   (kg1,              np_e,   ik:ik,kimg)
  real(kind=DP), intent(out) :: wfsd_l_3D(maxval(np_g1k_3D),np_e_3D,ik:ik,kimg)
  real(kind=DP)              :: wk1      (kg1,              neg,    ik:ik,kimg)
  integer, intent(in) :: ik

  integer :: i, j, l, is, js, jto, ncnt
  integer :: ierr

  call tstatc0_begin('z_interface_3D__________________',id_sname,1)

  ncnt = kg1*neg*kimg
  wk1 = 0.0d0

  do l = 1, kimg
     do j = 1, np_e
        js = nis_e(myrank_e) + j - 1
        do i = 1, iba(ik)
            wk1(i,js,ik,l) = wfsd_l(i,j,ik,l)
  enddo; enddo; enddo

  call mpi_allreduce(MPI_IN_PLACE,wk1,ncnt,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_CommGroup,ierr)

  do l = 1, kimg
     do j = 1, np_e_3D
        jto = neg_g(j)
        do i = 1, np_g1k_3D(ik)
           is = nis_g1k_3D(myrank_g_3D,ik) + i - 1
           wfsd_l_3D(i,j,ik,l) = wk1(is,jto,ik,l)
  enddo; enddo; enddo

   call tstatc0_end(id_sname)

  end subroutine decomp_wfsd_l_3D

!===============================================================================
  subroutine decomp_wfsd_l_r_3D(wfsd_l, wfsd_l_3D, ik)
  implicit none
  include 'mpif.h'
  real(kind=DP), intent(out) :: wfsd_l   (kg1,              np_e,   ik:ik,kimg)
  real(kind=DP), intent(in)  :: wfsd_l_3D(maxval(np_g1k_3D),np_e_3D,ik:ik,kimg)
  real(kind=DP)              :: wk1      (kg1,              neg,    ik:ik,kimg)
  integer, intent(in) :: ik

  integer :: i, j, l, is, js, jto, ncnt
  integer :: ierr

  call tstatc0_begin('z_interface_3D__________________',id_sname,1)

  ncnt = kg1*neg*kimg
  wk1 = 0.0d0

  do l = 1, kimg
     do j = 1, np_e_3D
        jto = neg_g(j)
        do i = 1, np_g1k_3D(ik)
           is = nis_g1k_3D(myrank_g_3D,ik) + i - 1
           wk1(is,jto,ik,l) = wfsd_l_3D(i,j,ik,l)
  enddo; enddo; enddo

  call mpi_allreduce(MPI_IN_PLACE,wk1,ncnt,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_CommGroup,ierr)

  do l = 1, kimg
     do j = 1, np_e
        js = nis_e(myrank_e) + j - 1
        do i = 1, iba(ik)
            wfsd_l(i,j,ik,l) = wk1(i,js,ik,l)
  enddo; enddo; enddo

   call tstatc0_end(id_sname)

  end subroutine decomp_wfsd_l_r_3D

!===============================================================================
  subroutine decomp_bsdr_l_3D(bsdr_l, bsdr_l_3D, ik)
  implicit none
  include 'mpif.h'
  real(kind=DP), intent(in)  :: bsdr_l   (np_e,   nlmta,ik:ik)
  real(kind=DP), intent(out) :: bsdr_l_3D(np_e_3D,nlmta,ik:ik)
  real(kind=DP)              :: wk1      (neg,    nlmta)
  integer, intent(in) :: ik
  integer :: j, l, js, jto, ncnt
  integer :: ierr

  call tstatc0_begin('z_interface_3D__________________',id_sname,1)

  ncnt = neg*nlmta
  wk1 = 0.0d0

  do j = 1, np_e
     js = nis_e(myrank_e) + j - 1
     wk1(js,:) = bsdr_l(j,:,ik)
  enddo

  call mpi_allreduce(MPI_IN_PLACE,wk1,ncnt,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_CommGroup,ierr)

  do j = 1, np_e_3D
     jto = neg_g(j)
     bsdr_l_3D(j,:,ik) = wk1(jto,:)
  enddo

   call tstatc0_end(id_sname)

  end subroutine decomp_bsdr_l_3D

!===============================================================================
  subroutine decomp_zajold_l_r_3D(zajold_l, zajold_l_3D, ik, ndavid)
  implicit none
  include 'mpif.h'
  real(kind=DP), intent(out) :: zajold_l   (kg1,              np_e,   kimg,ndavid)
  real(kind=DP), intent(in)  :: zajold_l_3D(maxval(np_g1k_3D),np_e_3D,kimg,ndavid)
  real(kind=DP)              :: wk1        (kg1,              neg,    kimg,ndavid)
  integer, intent(in) :: ik, ndavid

  integer :: i, j, l, is, js, jto, ncnt
  integer :: ierr

  call tstatc0_begin('z_interface_3D__________________',id_sname,1)

  ncnt = kg1*neg*kimg*ndavid
  wk1 = 0.0d0

  do l = 1, kimg
     do j = 1, np_e_3D
        jto = neg_g(j)
        do i = 1, np_g1k_3D(ik)
           is = nis_g1k_3D(myrank_g_3D,ik) + i - 1
           wk1(is,jto,l,:) = zajold_l_3D(i,j,l,:)
  enddo; enddo; enddo

  call mpi_allreduce(MPI_IN_PLACE,wk1,ncnt,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_CommGroup,ierr)

  do l = 1, kimg
     do j = 1, np_e
        js = nis_e(myrank_e) + j - 1
        do i = 1, iba(ik)
            zajold_l(i,j,l,:) = wk1(i,js,l,:)
  enddo; enddo; enddo

   call tstatc0_end(id_sname)

  end subroutine decomp_zajold_l_r_3D

!===============================================================================
 subroutine decomp_ngpt_l_r_3D(ngpt_l,ngpt_l_3D,n)
  include 'mpif.h'                                      ! MPI
  integer(kind=4), dimension(ista_kngp:iend_kngp,n) ,intent(out) :: ngpt_l
  integer(kind=4), dimension(ista_kngp_3D:iend_kngp_3D,n) ,intent(in):: ngpt_l_3D
  integer(kind=4), allocatable,dimension(:,:) :: send, recv
  integer(kind=4) :: i, j, ierr,n

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(1:kgp,n))
   allocate(recv(1:kgp,n))
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send(:,:) = 0.0d0
   do j = 1, n
      do i = ista_kngp_3D, iend_kngp_3D
         if (i > kgp) exit
         send(i,j) = ngpt_l_3D(i,j)
      end do
   end do

   call mpi_allreduce(send,recv,kgp*n,MPI_integer,MPI_SUM,mpi_ke_world,ierr)

   ngpt_l(:,:) = 0.0d0

   do j = 1, n
      do i = ista_kngp, iend_kngp
         if (i > kgp) exit
         ngpt_l(i,j) = recv(i,j)
      end do
   end do

   deallocate(send,recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_ngpt_l_r_3D

!===============================================================================
 subroutine decomp_ngpt_l_3D(ngpt_l,ngpt_l_3D,n)
  include 'mpif.h'                                      ! MPI
  integer(kind=4), dimension(ista_kngp:iend_kngp,n) ,intent(in) :: ngpt_l
  integer(kind=4), dimension(ista_kngp_3D:iend_kngp_3D,n) ,intent(out):: ngpt_l_3D
  integer(kind=4), allocatable,dimension(:,:) :: send, recv
  integer(kind=4) :: i, j, ierr,n

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(1:kgp,n))
   allocate(recv(1:kgp,n))
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send(:,:) = 0.0d0
   do j = 1, n
      do i = ista_kngp, iend_kngp
         if (i > kgp) exit
         send(i,j) = ngpt_l(i,j)
      end do
   end do

!   call mpi_allreduce(send,recv,kgp*n,MPI_integer,MPI_SUM,mpi_ke_world,ierr)
   call mpi_allreduce(send,recv,kgp*n,MPI_integer,MPI_SUM,mpi_k_world(myrank_k),ierr)

   ngpt_l_3D(:,:) = 0.0d0

   do j = 1, n
      do i = ista_kngp_3D, iend_kngp_3D
         if (i > kgp) exit
         ngpt_l_3D(i,j) = recv(i,j)
      end do
   end do

   deallocate(send,recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_ngpt_l_3D

!===============================================================================
  subroutine decomp_igfp_l_r_3D(igfp_l,igfp_l_3D)
  include 'mpif.h'                                      ! MPI
  integer(kind=4), dimension(ista_kngp:iend_kngp) ,intent(out) :: igfp_l
  integer(kind=4), dimension(ista_kngp_3D:iend_kngp_3D) ,intent(in):: igfp_l_3D
  integer(kind=4), allocatable,dimension(:) :: send, recv
  integer(kind=4) :: i, ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(1:kgp))
   allocate(recv(1:kgp))
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send(:) = 0.0d0
   do i = ista_kngp_3D, iend_kngp_3D
      if (i > kgp) exit
      send(i) = igfp_l_3D(i)
   end do

   call mpi_allreduce(send,recv,kgp,MPI_integer,MPI_SUM,mpi_ke_world,ierr)

   igfp_l(:) = 0.0d0

   do i = ista_kngp, iend_kngp
      if (i > kgp) exit
      igfp_l(i) = recv(i)
   end do

   deallocate(send,recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_igfp_l_r_3D

!===============================================================================
  subroutine decomp_igfp_l_3D(igfp_l,igfp_l_3D)
  include 'mpif.h'                                      ! MPI
  integer(kind=4), dimension(ista_kngp:iend_kngp) ,intent(in) :: igfp_l
  integer(kind=4), dimension(ista_kngp_3D:iend_kngp_3D) ,intent(out):: igfp_l_3D
  integer(kind=4), allocatable,dimension(:) :: send, recv
  integer(kind=4) :: i, ierr

   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(1:kgp))
   allocate(recv(1:kgp))
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send(:) = 0.0d0
   do i = ista_kngp, iend_kngp
      if (i > kgp) exit
      send(i) = igfp_l(i)
   end do

!   call mpi_allreduce(send,recv,kgp,MPI_integer,MPI_SUM,mpi_ke_world,ierr)
   call mpi_allreduce(send,recv,kgp,MPI_integer,MPI_SUM,mpi_k_world(myrank_k),ierr)

   igfp_l_3D(:) = 0.0d0

   do i = ista_kngp_3D, iend_kngp_3D
      if (i > kgp) exit
      igfp_l_3D(i) = recv(i)
   end do

   deallocate(send,recv)

   call tstatc0_end(id_sname)

  end subroutine decomp_igfp_l_3D

!===============================================================================
  subroutine decomp_eko_l_3D_new(inn,out,ik,nrvf_ordr,sort)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(:,:) ,intent(in) :: inn
  real(kind=DP), dimension(:,:) ,intent(out):: out
  integer(kind=4)               ,intent(in) :: ik
  integer(kind=4),dimension(neg,ista_k_3D:iend_k_3D),intent(in) :: nrvf_ordr
  character(len=4) ,intent(in) :: sort
  real(kind=DP), allocatable,dimension(:) :: send, recv
  integer(kind=4) :: i, iadd, ierr

   integer :: id_sname=-1
   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(neg),stat=ierr)
   allocate(recv(neg),stat=ierr)
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send = 0.0d0
   recv = 0.0d0

   do i = ista_e, iend_e
      if ((sort == "sort"))then
         iadd = nrvf_ordr(i,ik)
      else
         iadd = i
      endif
      send(iadd) = inn(i-ista_e+1,ik)
   enddo

   call mpi_allreduce(send,recv,neg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)

   do i = 1,np_e_3D
      out(i,ik) = recv(neg_g(i))
   enddo

   deallocate(send)
   deallocate(recv)

   call tstatc0_end(id_sname)

  endsubroutine decomp_eko_l_3D_new
!===============================================================================
  subroutine decomp_eko_l_r_3D_new(out,inn,ik,ordr,sort)
  include 'mpif.h'                                      ! MPI
  real(kind=DP), dimension(:,:) ,intent(in) :: inn
  real(kind=DP), dimension(:,:) ,intent(out):: out
  integer(kind=4)               ,intent(in) :: ik
  integer(kind=4),dimension(neg,ista_k_3D:iend_k_3D),intent(in) :: ordr
  character(len=4), intent(in) :: sort
  real(kind=DP), allocatable,dimension(:) :: send, recv
  integer(kind=4) :: i, iadd, ierr

   integer :: id_sname=-1
   call tstatc0_begin('z_interface_3D__________________',id_sname,1)

   allocate(send(neg),stat=ierr)
   allocate(recv(neg),stat=ierr)
   if (ierr /= 0) then
      print *,' ERROR Not allocated : ',ierr
      stop
   endif

   send = 0.0d0
   recv = 0.0d0
   if (myrank_g_3D == 0) then
      do i = 1, np_e_3D
         if (sort == "sort") then
            iadd = ordr(neg_g(i),ik)
         else
            iadd = neg_g(i)
         endif
         send(iadd) = inn(i,ik)
      enddo
   endif

   call mpi_allreduce(send,recv,neg,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_k_world(myrank_k),ierr)

   iadd = 0
   do i = ista_e, iend_e
      iadd = iadd + 1
      out(iadd,ik) = recv(i)
   enddo

   deallocate(send)
   deallocate(recv)

   call tstatc0_end(id_sname)

  endsubroutine decomp_eko_l_r_3D_new
!===============================================================================

#else
  contains
  subroutine dummy()
  end subroutine dummy
#endif

end module z_interface_3D
