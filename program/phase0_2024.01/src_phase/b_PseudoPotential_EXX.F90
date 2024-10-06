subroutine alloc_qitg_exx()
  use m_PseudoPotential,    only : nqitg
  use m_PlaneWaveBasisSet,  only : kgp
  use m_ES_ExactExchange,   only : nqmk, qitg_exx
  use m_Control_Parameters, only : nmax_G_hyb
  use m_Const_Parameters,   only : ON
  implicit none
#if defined(MEMORY_SAVE_EXX) && defined(MEMORY_SAVE_MORE_EXX)
  if(.not.allocated(qitg_exx)) allocate(qitg_exx(nmax_G_hyb,nqitg))
#else
  if(.not.allocated(qitg_exx)) allocate(qitg_exx(nmax_G_hyb,nqitg,nqmk))
#endif
end subroutine alloc_qitg_exx

#if defined(MEMORY_SAVE_EXX) && defined(MEMORY_SAVE_MORE_EXX)
subroutine qitgft_qmk(it,nmm_il3,mm_il3,qrsps_mm,lcmax,h)
   use m_Const_Parameters, only : DP
   use m_PseudoPotential, only : mmesh, nqitg_sp
   use m_ES_ExactExchange, only : qitgft_qmk_each_k_setup
   implicit none
   integer, intent(in) :: it
   integer, intent(in) :: nmm_il3(lcmax+1)
   integer, intent(in) :: mm_il3(nqitg_sp(it),lcmax+1)
   real(kind=DP), intent(in) :: qrsps_mm(mmesh,nqitg_sp(it))
   integer, intent(in) :: lcmax
   real(kind=DP), intent(in) :: h(it)
   call qitgft_qmk_each_k_setup(it,nmm_il3,mm_il3,qrsps_mm,lcmax,h)
   return
end subroutine qitgft_qmk
#else
subroutine qitgft_qmk(it,nmm_il3,mm_il3,qrsps_mm,lcmax,h)
  use m_Const_Parameters,   only : DP, PAI4, DELTA, ON
  use m_Control_Parameters, only : nmax_G_hyb
  use m_Crystal_Structure,  only : rltv
  use m_PlaneWaveBasisSet,  only : ngabc,kgp
  use m_PseudoPotential,    only : mmesh,nmesh,rmax,radr,wos,nqitg_sp
  use m_Parallelization,    only : MPI_CommGroup,ista_kngp,iend_kngp,npes,ierr,ista_kngp_exx,iend_kngp_exx
  use m_ES_ExactExchange,   only : nqmk, qmk, qitg_exx
  use m_Timing,             only : tstatc0_begin, tstatc0_end
  use mpi
  implicit none
  integer, intent(in) :: it
  integer, intent(in) :: nmm_il3(lcmax+1)
  integer, intent(in) :: mm_il3(nqitg_sp(it),lcmax+1)
  real(kind=DP), intent(in) :: qrsps_mm(mmesh,nqitg_sp(it))
  integer, intent(in) :: lcmax
  real(kind=DP), intent(in) :: h(it)

  integer :: mm0, i, il3, mm, mmp, n, idp, iq, ik
  real(kind=DP) :: ttr(6), kg(3)
  real(kind=DP) :: qitg_sh, g2, gabs
  real(kind=DP) :: wkx(mmesh), wky(mmesh)
  real(kind=DP), allocatable :: gqmk_l(:) ! d(ista_kngp:iend_kngp)
  real(kind=DP), allocatable :: qitg_exx_l(:,:) ! d(ista_kngp:iend_kngp,nqitg_sp(it))
  real(kind=DP), allocatable :: qitg_t(:,:) ! d(kgp,nqitg_sp(it))

  integer :: iend_kngp0

!  include 'mpif.h'                                      ! MPI

  integer,save  :: id_sname = -1
  call tstatc0_begin('qitgft_qmk ',id_sname)

  call getttr(rltv,ttr)

  allocate(qitg_exx_l(ista_kngp:iend_kngp,nqitg_sp(it))); qitg_exx_l = 0.0d0
  allocate(gqmk_l(ista_kngp:iend_kngp))

  mm0 = 0
  do i = 1, it-1
     mm0 = mm0 + nqitg_sp(i)
  end do
  iend_kngp0 = iend_kngp
  if(iend_kngp0.gt.nmax_G_hyb) iend_kngp0 = nmax_G_hyb

  do ik=1,nqmk
     do i = ista_kngp, iend_kngp0
        kg(1:3) = qmk(ik,1:3) + ngabc(i,1:3)
        g2          = ttr(1)*kg(1)*kg(1) &
        &           + ttr(2)*kg(2)*kg(2) &
        &           + ttr(3)*kg(3)*kg(3) &
        &           + ttr(4)*kg(1)*kg(2) &
        &           + ttr(5)*kg(2)*kg(3) &
        &           + ttr(6)*kg(3)*kg(1)
        gqmk_l(i) = sqrt(g2)
     end do
     do il3 = 1, lcmax+1
        if(nmm_il3(il3) <= 0) cycle
        do i = ista_kngp, iend_kngp0
           gabs = gqmk_l(i)
           if(gabs < DELTA) then
              idp = nmesh(it)+1
           else
              idp = ceiling(nmesh(it) - dlog(rmax(it)*gabs)/h(it))
           end if

! ==== ASMS ====
           if ( idp > nmesh(it) +1 ) idp = nmesh(it) +1
! ==== ASMS ====

           wkx(1:nmesh(it)) = gabs*radr(1:nmesh(it))
           call dsjnvn(il3-1,nmesh(it),idp,wkx,wky)

           do mm = 1, nmm_il3(il3)
              mmp = mm_il3(mm,il3)
              if ( mmp <= 0 ) cycle

              qitg_sh = 0.d0
              do n = 1, nmesh(it)
                 qitg_sh = qitg_sh + wos(n)*qrsps_mm(n,mmp)*wky(n)
              end do
              qitg_exx_l(i,mmp) = qitg_sh*PAI4
           end do
        end do
     end do
     if(npes>1) then
        allocate(qitg_t(nmax_G_hyb,nqitg_sp(it))); qitg_t = 0.0d0
        do iq=1,nqitg_sp(it)
           do i=ista_kngp,iend_kngp0
              qitg_t(i,iq) = qitg_exx_l(i,iq)
           end do
        end do
        call mpi_allreduce(qitg_t,qitg_exx(1,mm0+1,ik),nmax_G_hyb*nqitg_sp(it) &
         &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
        do iq=1,nqitg_sp(it)
           do i=ista_kngp,iend_kngp0
              qitg_exx(i,mm0+iq,ik) = qitg_t(i,iq)
           enddo
        enddo
        deallocate(qitg_t)
     else
        do iq=1,nqitg_sp(it)
           do i=1,nmax_G_hyb
              qitg_exx(i,mm0+iq,ik) = qitg_exx_l(i,iq)
           end do
        end do
     end if

  end do ! ik

  deallocate(qitg_exx_l)
  deallocate(gqmk_l)

  call tstatc0_end(id_sname)
end subroutine qitgft_qmk
#endif
