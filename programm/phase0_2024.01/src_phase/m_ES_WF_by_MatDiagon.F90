!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 581 $)
!
!  MODULE: m_ES_WF_by_MatDiagon
!
!  AUTHOR(S): T. Yamasaki, M. Okamoto,   August/22/2003
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
#ifdef __TIMER_SUB__
#   define __TIMER_SUB_START(a)  call timer_sta(a)
#   define __TIMER_SUB_STOP(a)   call timer_end(a)
#else
#   define __TIMER_SUB_START(a)
#   define __TIMER_SUB_STOP(a)
#endif
#ifdef __TIMER_DO__
#   define __TIMER_DO_START(a)   call timer_sta(a)
#   define __TIMER_DO_STOP(a)    call timer_end(a)
#else
#   define __TIMER_DO_START(a)
#   define __TIMER_DO_STOP(a)
#endif
#ifdef __TIMER_COMM__
#   define __TIMER_COMM_START_w_BARRIER(str,a)   call timer_barrier(str) ;   call timer_sta(a)
#   define __TIMER_COMM_START(a)       call timer_sta(a)
#   define __TIMER_COMM_STOP(a)        call timer_end(a)
#else
#   define __TIMER_COMM_START_w_BARRIER(str,a)
#   define __TIMER_COMM_START(a)
#   define __TIMER_COMM_STOP(a)
#endif

#ifndef SX
#define DGEMM__       DGEMM
#endif

#ifndef NO_MATDIAGON_DGEMM
#define MATDIAGON_DGEMM
#endif

module m_ES_WF_by_MatDiagon
! $Id: m_ES_WF_by_MatDiagon.F90 581 2018-08-01 08:38:42Z jkoga $
  use m_Electronic_Structure,only : zaj_l, eko_l, vlhxcQ,vlhxc_l &
#ifdef SAVE_FFT_TIMES
       &                          , status_saved_phifftr &
#endif
       &                          , m_ES_sort_eigen_vals_noncl &  ! <<== added by K. Tagami ============= 11.0
       &                          , m_ES_wd_zaj_small_portion
  use m_ES_nonlocal,         only : m_ES_betar_dot_WFs
  use m_NonLocal_Potential,  only : snl
  use m_PlaneWaveBasisSet,   only : kg1_ext,kg1,iba2,iba,nmatsz,nmatsz2,nbmat,nbmat2,nbase,ngabc,ttr,igpo &
       &                          , n_rGsv, nbase_gamma, kg2_gamma &
       &                          , m_pwBS_alloc_igpo &
       &                          , m_pwBS_dealloc_igpo &
       &                          , m_pwBS_GminusGmapfunction &
       &                          , m_pwBS_set_gmaxs &
       &                          , m_pwBS_alloc_nbmat_and_iba2 &
       &                          , m_pwBS_dealloc_nbmat_and_iba2 &
       &                          , m_pwBS_mat_for_each_WF
  use m_PseudoPotential,     only : modnrm,n_non0_lmtxlmt,index_lmt1_lmt2,lmtt &
       &                          , ltp,mtp,ilmt,dion,q,ivanl, nlmtt &
       &                          , ipaw,dion_paw
  use m_Kpoints,             only : kv3, vkxyz, k_symmetry
  use m_Ionic_System,        only : ntyp, natm, ityp, iwei, pos, ivan
  use m_Timing,              only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters,  only : nspin,kimg,af,neg,iprimatdiagon,intzaj &
       &                          , n_matrix_size,gmax,eps_solve_Hx_eq_ex,gmaxs_given,ekmode &
#ifdef SAVE_FFT_TIMES
       &                          , skip_alloc_phonon, sw_save_fft
#else
       &                          , skip_alloc_phonon
#endif
  use m_Const_Parameters,    only : DP,PAI2,EXECUT,CMPLDP,BUCS,SmallestPositiveNumber,by_matrix_diagon,ON, GAMMA &
       &                          , GAMMA_base_symmetrization, OLD
  use m_Parallelization,     only : MPI_CommGroup,myrank_k,map_k,ierr,npes, np_e &
       &                          , ista_kngp,iend_kngp,ista_e,iend_e &
       &                          , istep_e,map_z

! =========================== added by K. Tagami =================== 11.0
  use m_Const_Parameters,     only : Neglected, BuiltIn, OFF, SKIP
  use m_Control_Parameters,  only : noncol, ndim_chgpot, ndim_spinor, SpinOrbit_mode, sw_hubbard
  use m_PseudoPotential,     only : dion_scr_noncl, q_noncl, nlmt
  use m_ES_NonCollinear,      only : m_ES_MagMom_to_DensMat_Gspace, &
       &                             m_ES_MagMom_to_DensMat_vlhxcl
! ================================================================== 11.0
  use mpi

  implicit none

  real(kind=DP),private,allocatable,target,dimension(:,:,:):: hsmat !d(nmatsz,nmatsz,kimg)
  real(kind=DP),private,allocatable,target,dimension(:,:,:):: hsmat_s !d(nmatsz,nmatsz,kimg)
  real(kind=DP),private,allocatable,target,dimension(:,:):: hsmat1
  real(kind=DP),private,allocatable,target,dimension(:,:):: hsmat2
  real(kind=DP),private,allocatable,dimension(:,:,:):: Linv  !d(nmatsz,nmatsz,kimg)
  real(kind=DP),private,allocatable,dimension(:,:)  :: ldag  !d(nmatsz,kimg)
  real(kind=DP),private,allocatable,dimension(:,:)  :: vlhxc   !d(nmatsz2,kimg)
  real(kind=DP),private,allocatable,dimension(:)    :: zfcos,zfsin !d(nmatsz2)
  real(kind=DP),private,allocatable,dimension(:,:,:)  :: zaj_mat !d(nmatsz,np_e,kimg)


  real(kind=DP),private,allocatable :: vlhxc_ssrep(:,:,:)   ! <<== added by K. Tagami =============== 11.0
  integer :: nmatsz_noncl                                   ! <<== added by K. Tagami =============== 11.0
  integer :: sw_use_zhegvx = on                             ! <<== added by K. Tagami =============== 11.0

  logical, private :: reduced_basis_mode = .false.

!  include 'mpif.h'

contains

  subroutine m_ESmat_set_reduced_basis_mode(reduced)
      logical, intent(in) :: reduced
      reduced_basis_mode = reduced
  end subroutine m_ESmat_set_reduced_basis_mode

  subroutine alloc_matrices()
!f    allocate(hsmat(nmatsz,nmatsz,kimg)); hsmat = 0.d0
    allocate(vlhxc(nmatsz2,kimg))
    allocate(zfcos(0:nmatsz2))
    allocate(zfsin(0:nmatsz2))
!f    if(modnrm == EXECUT) then
!f       allocate(linv(nmatsz,nmatsz,kimg))
!f       allocate(ldag(nmatsz,kimg))
!f    end if
!f    allocate(zaj_mat(nmatsz,np_e,kimg))
! ============================ modified by K. Tagami ================= 11.0
!    allocate(hsmat(nmatsz,nmatsz,kimg)); hsmat = 0.d0
!    if(modnrm == EXECUT) then
!       allocate(linv(nmatsz,nmatsz,kimg))
!       allocate(ldag(nmatsz,kimg))
!    end if
!    allocate(zaj_mat(nmatsz,np_e,kimg))
!
    if ( noncol ) then
      nmatsz_noncl = nmatsz *ndim_spinor
      allocate( hsmat( nmatsz_noncl, nmatsz_noncl, kimg )); hsmat =0.d0

      if ( modnrm == EXECUT ) then
         if ( sw_use_zhegvx == ON ) then
            allocate( hsmat_s( nmatsz_noncl, nmatsz_noncl, kimg) ); hsmat_s = 0.0d0
         else
            allocate( linv( nmatsz_noncl, nmatsz_noncl, kimg) )
            allocate( ldag( nmatsz_noncl, kimg ) )
         end if
      endif

      allocate( zaj_mat( nmatsz_noncl, np_e, kimg ) )

      allocate( vlhxc_ssrep( ista_kngp:iend_kngp,kimg,ndim_chgpot) )
      vlhxc_ssrep = 0.0d0

    else
      allocate(hsmat(nmatsz,nmatsz,kimg)); hsmat = 0.d0
      if ( modnrm == EXECUT ) then
         allocate(linv(nmatsz,nmatsz,kimg))
         allocate(ldag(nmatsz,kimg))
      end if
      allocate(zaj_mat(nmatsz,np_e,kimg))
    endif
! ===================================================================== 11.0
  end subroutine alloc_matrices

  subroutine dealloc_matrices()
    if(allocated(zaj_mat)) deallocate(zaj_mat)
    if(allocated(ldag))  deallocate(ldag)
    if(allocated(linv))  deallocate(linv)
    if(allocated(zfsin)) deallocate(zfsin)
    if(allocated(zfcos)) deallocate(zfcos)
    if(allocated(vlhxc)) deallocate(vlhxc)
    if(allocated(hsmat)) deallocate(hsmat)
    if ( allocated(vlhxc_ssrep) ) deallocate( vlhxc_ssrep )  ! <<== added by K. Tagami ========== 11.0
    if(allocated(hsmat_s)) deallocate(hsmat_s)
  end subroutine dealloc_matrices

  subroutine m_ESmat_solve_Hx_eq_eSx(nfout,iteration,iteration_electronic)
    integer, intent(in) :: nfout, iteration, iteration_electronic

    integer,parameter :: PRINT_LEVEL = 2
    integer   :: ispin, ik, iksnl, n=0
    integer   :: id_sname = -1
    logical   :: unitcell_can_change
    call tstatc0_begin('m_ESmat_solve_Hx_eq_eSx ',id_sname,1)

    if(iprimatdiagon >= 2) write(nfout,'(" !! iteration = ",i6, " m_ESmat_solve_Hx_eq_eSx")') iteration

    if(ekmode == ON ) then
       if(iteration_electronic == 1) then
          call m_pwBS_set_gmaxs(n_matrix_size,gmaxs_given)
       else
          call m_pwBS_set_gmaxs(n,gmax)
       end if
    else
       if(iteration == 1 .and. intzaj == by_matrix_diagon) then
          call m_pwBS_set_gmaxs(n_matrix_size,gmaxs_given)
       else if(skip_alloc_phonon .and. iteration_electronic == 1 .and. intzaj == by_matrix_diagon) then
          call m_pwBS_set_gmaxs(n_matrix_size,gmaxs_given)
       else if(unitcell_can_change() .and. iteration_electronic == 1 &
         &    .and. intzaj == by_matrix_diagon) then
          call m_pwBS_set_gmaxs(n_matrix_size,gmaxs_given)
       else if (reduced_basis_mode) then
          call m_pwBS_set_gmaxs(n_matrix_size,gmaxs_given)
          reduced_basis_mode = .false.
       else
          call m_pwBS_set_gmaxs(n,gmax)
       end if
    end if

    call m_pwBS_alloc_nbmat_and_iba2()
    call m_pwBS_mat_for_each_WF() ! -> nbmat,nbmat2,iba2

    call m_pwBS_alloc_igpo()
    call m_pwBS_GminusGmapfunction()  ! -> igpo, nmatsz2
    if(iprimatdiagon >= 3) write(nfout,'(" nmatsz2 = ",i9)') nmatsz2
    if(iprimatdiagon >= PRINT_LEVEL) write(nfout,'(" nmatsz, kg1_ext = ", 2i9)') nmatsz, kg1_ext

    call alloc_matrices()

    do ispin = 1, nspin, af+1
       call prjvlhxc_l2vlhxc_t()   ! vlhxc_l -> vlhxc
       do ik = ispin, kv3-nspin+ispin,nspin
          if(map_k(ik) /= myrank_k) cycle
          if(iprimatdiagon >= PRINT_LEVEL) write(nfout,'(" iba2(",i3,") = ", i9)') ik, iba2(ik)

          iksnl = (ik-1)/nspin + 1
          hsmat = 0.d0
          if(modnrm == EXECUT) then
             if(iprimatdiagon >=3 ) write(nfout,'(" !! modnrm == EXECUT (m_ESmat_solve_Hx_eq_eSx)")')
             call ovmtrx() ; if(iprimatdiagon >= 3) call wd_hsmat(" -- ovmtrx --")
             call cholde() ; if(iprimatdiagon >= 3) call wd_hsmat(" -- cholde --")   ! -> linv
          end if
          call nlmtrx()           ! (nlmtrx) -> hsmat <- snl
                          if(iprimatdiagon >= 3) call wd_hsmat(" -- nlmtrx --")
          call amtrxkinetic()     !  -> hsmat
                          if(iprimatdiagon >= 3) call wd_hsmat(" -- amtrxkinetic --")
          call amtrxlocalpt()     !  -> hsmat
                          if(iprimatdiagon >= 3) call wd_hsmat(" -- amtrxlocalpt --")
          if(modnrm == EXECUT) call lhl() ! -> hsmat
                          if(modnrm == EXECUT .and. iprimatdiagon >= 3) call wd_hsmat(" -- lhl --")
      ! [LAPACK] with LAPACK library
      !***** LAPACK [begin] ***** M.Okamoto
          call solve_Hx_eq_ex_LAPACK()
      !***** LAPACK [end] ***** M.Okamoto
          if(iprimatdiagon >= 3) call wd_hsmat(" -- after diagonalization --")
          if(modnrm == EXECUT) call lx() ! -> zaj_mat
          ! Debug
!!$          call zaj_mapping()
          if(k_symmetry(ik) == GAMMA .or. k_symmetry(ik) == GAMMA_base_symmetrization) call phase_mult()
          call cp_zaj_mat_to_zaj_l() ! zaj_mat -> zaj_l
       end do
    end do

    if(iprimatdiagon >= 3) then
       do ik = 1, kv3
          if(map_k(ik) /= myrank_k) cycle
          call m_ES_wd_zaj_small_portion(nfout,ik," (m_ESmat_solve_Hx_eq_eSx)",26)
!!$          write(nfout,'(" ik = ", i5)') ik
       end do
    end if

    call dealloc_matrices()
    call m_pwBS_dealloc_igpo()
    call m_pwBS_dealloc_nbmat_and_iba2()

    call m_ES_betar_dot_WFs(nfout) ! (fsrfsi) -> fsr_l, fsi_l

    call tstatc0_end(id_sname)
  contains
    ! Contained subroutines
    !  1)    subroutine wd_hsmat
    !  1)    subroutine wd_hsmat
    !  2)    subroutine solve_Hx_eq_ex_LAPACK
    !  3)    subroutine solve_Hx_eq_ex
    !  4)    subroutine calphase
    !  5)    subroutine nlmtrx
    !  6)    subroutine ovmtrx
    !  7)    subroutine cholde
    !  8)    subroutine amtrxkinetic
    !  9)    subroutine prjvlhxc_l2vlhxc_t
    ! 10)    subroutine amtrxlocalpt
    ! 11)    subroutine lhl
    ! 12)    subroutine lx
    ! 13)    subroutine phase_mult
    ! 14)    subroutine cp_zaj_mat_to_zaj_l

    subroutine wd_hsmat(str)
      character(len=*),intent(in) :: str
      integer :: i, imax, j, jmax
      integer, parameter :: NMATSIZE = 10
      write(nfout,'(a30," ik= ",i5)') str,ik
      imax = min(NMATSIZE, nmatsz)
      jmax = min(NMATSIZE, nmatsz)
      do j = 1, jmax
         write(nfout,'(" hsmat(1:*,",i3,",1): ",10f8.4,99(/19x,10f8.4))') j &
              & , (hsmat(i,j,1),i=1,imax)
      end do
      if(kimg == 2) then
         do j = 1, jmax
            write(nfout,'(" hsmat(1:*,",i3,",2): ",10f8.4,99(/19x,10f8.4))') j &
                 & , (hsmat(i,j,2),i=1,imax)
         end do
      end if
    end subroutine wd_hsmat

!!$    subroutine zaj_mapping()
!!$      integer :: ib, i, j, ib2
!!$      if(iba2(ik) < iba(ik)) then
!!$         write(nfout,'(" << zaj_mapping >> (nbmat -> nbase)")')
!!$         do ib = ista_e, iend_e, istep_e
!!$            ib2 = map_z(ib)
!!$            zfcos = 0.d0
!!$            do j = 1, iba2(ik)
!!$               do i = 1, iba(ik)
!!$                  if(nbase(i,ik) == nbmat(j,ik)) zfcos(i) = zaj_l(j,ib2,ik,1)
!!$               end do
!!$            end do
!!$            do i = 1, iba(ik)
!!$               zaj_l(i,ib2,ik,1) = zfcos(i)
!!$            end do
!!$         end do
!!$         if(kimg == 2) then
!!$            do ib = ista_e, iend_e, istep_e
!!$               ib2 = map_z(ib)
!!$               zfcos = 0.d0
!!$               do j = 1, iba2(ik)
!!$                  do i = 1, iba(ik)
!!$                     if(nbase(i,ik) == nbmat(j,ik)) zfcos(i) = zaj_l(j,ib2,ik,2)
!!$                  end do
!!$               end do
!!$               do i = 1, iba(ik)
!!$                  zaj_l(i,ib2,ik,2) = zfcos(i)
!!$               end do
!!$            end do
!!$         end if
!!$      end if
!!$    end subroutine zaj_mapping

   ! [LAPACK] with LAPACK library
   !***** LAPACK [begin] ***** M.Okamoto
    subroutine solve_Hx_eq_ex_LAPACK()
      integer :: n, lda, ne, ne_found, il, iu, info, ipos, i, j, ib, lwork, k
      real(kind=DP) :: abstol, vl, vu, sum, dlamch
      integer,allocatable ::              ifail(:), iwork(:)
      real(kind=DP),allocatable ::        amat(:,:), evec(:), rvmat(:,:), rwork(:)
      complex(kind=CMPLDP),allocatable :: cmat(:,:), cvmat(:,:), cwork(:)
      complex(kind=CMPLDP),parameter ::   IMG = (0.d0,1.d0)
      complex(kind=CMPLDP) :: ctemp(1)
      zaj_mat(:,:,:) = 0.d0
      lda = nmatsz
      vl = -1.d10 ; vu = 1.d10 ; il = 1 ; iu = neg   !! vl = -1.d10 ; vu = 1.d10 ; il = ista_e ; iu = iend_e
      ne = iu - il + 1
      abstol = 2*dlamch('S')
      n = iba2(ik)
      if (kimg == 1) then
         lwork = max(1,12*n)
        !+++++++++++++++++++++++++++++++++++++++++++++++++
         allocate(evec(n),iwork(5*n),ifail(n))
         allocate(amat(lda,n),rvmat(lda,ne),rwork(lwork)); amat = 0.d0
        !+++++++++++++++++++++++++++++++++++++++++++++++++
         amat(1:n,1:n) = hsmat(1:n,1:n,1)  ! initialization of amat

         if(iprimatdiagon>= 2)  call check_LAPACKinputs("solve_Hx_eq_ex_LAPACK "  &
              &                                   ,nfout,n,lda,vl,vu,il,iu,lwork,abstol,RMat=amat)
         call dsyevx('V','I','U',n,amat,lda,vl,vu,il,iu,abstol &
              &    , ne_found,evec,rvmat,lda,rwork,lwork,iwork,ifail,info)

         if(iprimatdiagon>= 2) call check_LAPACKoutputs("solve_Hx_eq_ex_LAPACK ",nfout,ik,n,ne,lda,evec,RMat=rvmat)
         call check_LAPACKerror(ne_found,ne,info)

         do ib = ista_e, iend_e, istep_e
            j = map_z(ib)
            do i = 1,n
               zaj_mat(i,j,1) = rvmat(i,ib) !!$ zaj_l(i,j,ik,1) = rvmat(i,ib-ista_e+1) ! M.Okamoto (August 22,2003)
            end do
            eko_l(j,ik) = evec(ib)          !!$ eko_l(j,ik) = evec(ib-ista_e+1) ! M.Okamoto (August 22,2003)
         end do
         if(af/=0) eko_l(:,ik+1) = eko_l(:,ik)
        !+++++++++++++++++++++++++++++
         deallocate(amat,rvmat,rwork)
        !+++++++++++++++++++++++++++++
      else if (kimg == 2) then
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         allocate(evec(n),iwork(5*n),ifail(n))
         allocate(cmat(lda,n),cvmat(lda,ne),rwork(7*n))
         evec = 0.d0; iwork = 0; ifail = 0
         cmat = dcmplx(0.d0,0.d0); cvmat = dcmplx(0.d0,0.d0); rwork = 0.d0;
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         do i = 1,n
         do j = 1,n
            cmat(i,j) = dcmplx(hsmat(i,j,1),hsmat(i,j,2))
         end do
         end do
! === DEBUG by tkato 2011/09/06 ================================================
         call zheevx('V','I','U',n,cmat,lda,vl,vu,il,iu,abstol, &
            ne_found,evec,cvmat,lda,ctemp,-1,rwork,iwork,ifail,info)
         lwork = int(real(ctemp(1)))
         allocate(cwork(lwork)); cwork = dcmplx(0.d0,0.d0)
! ==============================================================================
         if(iprimatdiagon>= 2)  call check_LAPACKinputs("solve_Hx_eq_ex_LAPACK "  &
              &                                   ,nfout,n,lda,vl,vu,il,iu,lwork,abstol,CMat=cmat)
         call zheevx('V','I','U',n,cmat,lda,vl,vu,il,iu,abstol, &
            ne_found,evec,cvmat,lda,cwork,lwork,rwork,iwork,ifail,info)

         if(iprimatdiagon>= 2) call check_LAPACKoutputs("solve_Hx_eq_ex_LAPACK ",nfout,ik,n,ne,lda,evec,CMat=cvmat)
         call check_LAPACKerror(ne_found,ne,info)

         do ib = ista_e, iend_e, istep_e
            j = map_z(ib)
            do i = 1,n
               zaj_mat(i,j,1) = real(cvmat(i,ib))
               zaj_mat(i,j,2) = aimag(cvmat(i,ib))
               !!$ zaj_l(i,j,ik,1)    = real(cvmat(i,ib-ista_e+1))  ! M.Okamoto (August 22,2003)
               !!$ zaj_l(i,j,ik,kimg) = aimag(cvmat(i,ib-ista_e+1)) ! M.Okamoto (August 22,2003)
            end do
            eko_l(j,ik) = evec(ib)    !! eko_l(j,ik) = evec(ib-ista_e+1) ! M.Okamoto (August 22,2003)
         end do
         if(af/=0) eko_l(:,ik+1) = eko_l(:,ik)
        !+++++++++++++++++++++++++++++++++++
         deallocate(cmat,cvmat,cwork,rwork)
        !+++++++++++++++++++++++++++++++++++
      end if
      if(iprimatdiagon >= 2) call wd_CalculatedEigenvalues(nfout,ik,n,ne,evec,"solve_Hx_eq_ex")

      sum = 0.d0
      do i = 1,ne
         sum = sum + abs(evec(i))
      end do
      if (sum < SmallestPositiveNumber*1.d5) call phase_error_with_msg(nfout, ' ! illegal evec (solve_Hx_eq_ex)'&
                                                                      , __LINE__, __FILE__)
      if (iprimatdiagon >= 2) write(nfout,'(" -- sum = ",f20.8)') sum
     !+++++++++++++++++++++++++++++
      deallocate(evec,iwork,ifail)
     !+++++++++++++++++++++++++++++
    end subroutine solve_Hx_eq_ex_LAPACK
   !***** LAPACK [end] ***** M.Okamoto

    subroutine calphase(ia)
      integer, intent(in) :: ia
      integer       :: i
      real(kind=DP) :: ph
      do i = 1, nmatsz2
         ph = (pos(ia,1)*ngabc(i,1)+pos(ia,2)*ngabc(i,2)+pos(ia,3)*ngabc(i,3))*PAI2
         zfcos(i) = cos(ph)
         zfsin(i) = sin(ph)
      end do
      zfcos(0) = 0.d0; zfsin(0) = 0.d0
    end subroutine calphase

    subroutine nlmtrx
      real(kind=DP),pointer,dimension(:,:,:) :: a
      real(kind=DP),allocatable,dimension(:,:)  :: s
      integer, pointer, dimension(:,:)       :: ng

      complex(kind=CMPLDP),parameter :: zi = (0.d0,1.d0)
      integer :: i,j,it,ia,u,v,n,lmt1,il1,lmt2,il2,j1,l1,l2,l3,ip,y,i1
      integer :: i2,j2,m
      real(kind=DP) :: fr,f21,f12,t

      a => hsmat
      ng => ngabc
!!$      s => snl(:,:,iksnl)
      allocate(s(kg1_ext,nlmtt)); s = 0.d0

      if(k_symmetry(ik) == GAMMA) then
         if(iprimatdiagon >= 2) write(nfout,'(" k_symmetry(ik) == GAMMA <<nlmtrx>>")')
         if(iprimatdiagon >= 2) write(nfout,'(" kg2_gamma, iba(ik), iba2(ik) = ",3i8)') kg2_gamma,iba(ik),iba2(ik)
         do j = 1, nlmtt
            ! --- finding it and lmt
            Finding: do u = 1, ntyp
               do ip = 1, ilmt(u)
                  if(lmtt(ip,u) == j) then
                     lmt1 = ip
                     it = u
                     exit Finding
                  end if
               end do
            end do Finding
            il1 = ltp(lmt1,it)
            if(iprimatdiagon >= 2) write(nfout,'(" nlmtt, j, lmt1, it, il1 = ",5i8)') nlmtt,j,lmt1,it,il1
            s(1,j) = snl(1,j,iksnl)
            if(mod(il1,2) == 0) then
               do i = 2, kg2_gamma
                  i1 = nbase(i,ik)
                  i2 = nbase_gamma(i,2)
                  s(i1,j) = snl(i,j,iksnl)
                  s(i2,j) = -snl(i,j,iksnl)
               end do
            else
               do i = 2, kg2_gamma
                  i1 = nbase(i,ik)
                  i2 = nbase_gamma(i,2)
                  s(i1,j) = snl(i,j,iksnl)
                  s(i2,j) = snl(i,j,iksnl)
               end do
            end if
         end do
      else
         if(iprimatdiagon >= 2) write(nfout,'(" k_symmetry(ik) /= GAMMA <<nlmtrx>>")')
         do j = 1, nlmtt
            do i = 1, iba(ik)
               s(i,j) = snl(i,j,iksnl)
            end do
         end do
      end if

      a = 0.d0
      do ia = 1, natm
         it = ityp(ia)
         call calphase(ia)  ! -> zfcos, zfsin
         do ip = 1, n_non0_lmtxlmt(it)
            lmt1  = index_lmt1_lmt2(ip,it,1)
            lmt2  = index_lmt1_lmt2(ip,it,2)
            u     = lmtt(lmt1,it)
            v     = lmtt(lmt2,it)
            il1   = ltp(lmt1,it)
            il2   = ltp(lmt2,it)
            if(lmt2 <= lmt1) cycle
!!$            if(ivanl(il1,it) /= 1 .and. &
            if(ivan(it) /= 1 .and. &
                 & (il1 /= il2 .or. mtp(lmt1,it) /= mtp(lmt2,it))) cycle
            if(mod(il1+il2,2) == 0) then
               if(ipaw(it)==0) then
                   fr  = dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)
               else
                   fr  = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
               endif
               if(il1==il2) then
                  f21 = fr*iwei(ia)
                  do j = 1, iba2(ik)
                     j1 = nbmat(j,ik)
                     j2 = nbmat2(j,ik)
                     l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)
                     if(kimg == 1) then
                        do i = j, iba2(ik)
                           n = nbmat(i,ik)
                           m = nbmat2(i,ik)
                           y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                           a(i,j,1) = a(i,j,1)+f21*(s(m,u)*s(j2,v)+s(m,v)*s(j2,u))*zfcos(y)
                        end do
                     else
                        do i = j, iba2(ik)
                           n = nbmat(i,ik)
                           m = nbmat2(i,ik)
                           y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                           t  = f21*(s(m,u)*s(j2,v)+s(m,v)*s(j2,u))
                           a(i,j,1) = a(i,j,1) + t*zfcos(y)
                           a(i,j,2) = a(i,j,2) - t*zfsin(y)
                        end do
                     end if
                  end do
               else
                  f21 = real(zi**(il2-il1))*fr*iwei(ia)
                  f12 = real(zi**(il1-il2))*fr*iwei(ia)
!!$                  write(nfout,'(" zi**(il1-il2) = ",2f16.8)') zi**(il1-il2)
                  do j = 1, iba2(ik)
                     j1 = nbmat(j,ik)
                     j2 = nbmat2(j,ik)
                     l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)
                     if(kimg == 1) then
                        do i = j, iba2(ik)
                           n = nbmat(i,ik)
                           m = nbmat2(i,ik)
                           y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                           a(i,j,1) = a(i,j,1)+(f21*s(m,u)*s(j2,v)+f12*s(m,v)*s(j2,u))*zfcos(y)
                        end do
                     else
                        do i = j, iba2(ik)
                           n = nbmat(i,ik)
                           m = nbmat2(i,ik)
                           y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                           t  = f21*s(m,u)*s(j2,v)+f12*s(m,v)*s(j2,u)
                           a(i,j,1) = a(i,j,1) + t*zfcos(y)
                           a(i,j,2) = a(i,j,2) - t*zfsin(y)
                        end do
                     end if
                  end do
               end if
            else if(mod(il1+il2,2) == 1) then
!!$               fr = vlhxcQ(lmt1,lmt2,ia,ispin)
               if(ipaw(it)==0) then
                   fr = vlhxcQ(lmt1,lmt2,ia,ispin)
               else
                   fr = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
               end if
               f21 = -aimag(zi**(il2-il1))*fr*iwei(ia)
               f12 = -aimag(zi**(il1-il2))*fr*iwei(ia)
               do j = 1, iba2(ik)
                  j1 = nbmat(j,ik)
                  j2 = nbmat2(j,ik)
                  l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)
                  if(kimg == 1) then
                     do i = j, iba2(ik)
                        n = nbmat(i,ik)
                        m = nbmat2(i,ik)
                        y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                        a(i,j,1) = a(i,j,1)+ &
                            & (f21*s(m,u)*s(j2,v)+f12*s(m,v)*s(j2,u))*zfsin(y)
                     end do
                  else
                     do i = j, iba2(ik)
                        n = nbmat(i,ik)
                        m = nbmat2(i,ik)
                        y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                        t = f21*s(m,u)*s(j2,v)+f12*s(m,v)*s(j2,u)
                        a(i,j,1) = a(i,j,1) + t*zfsin(y)
                        a(i,j,2) = a(i,j,2) + t*zfcos(y)
                     end do
                  end if
               end do
            end if
         end do
         do lmt1 = 1, ilmt(it)
            if(ipaw(it)==0) then
                fr = (dion(lmt1,lmt1,it) + vlhxcQ(lmt1,lmt1,ia,ispin))*iwei(ia)
            else
                fr = (dion_paw(lmt1,lmt1,ispin,ia) + vlhxcQ(lmt1,lmt1,ia,ispin))*iwei(ia)
            endif
            u  = lmtt(lmt1,it)
            do j = 1, iba2(ik)
               j1 = nbmat(j,ik)
               j2 = nbmat2(j,ik)
               l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)
               if(kimg == 1) then
                  do i = j, iba2(ik)
                     i1 = nbmat(i,ik)
                     i2 = nbmat2(i,ik)
                     y = igpo(ng(i1,1)-l1,ng(i1,2)-l2,ng(i1,3)-l3)
                     a(i,j,1) = a(i,j,1) + fr*s(i2,u)*s(j2,u)*zfcos(y)
                  end do
               else if(kimg == 2) then
                  do i = j, iba2(ik)
                     i1 = nbmat(i,ik)
                     i2 = nbmat2(i,ik)
                     y = igpo(ng(i1,1)-l1,ng(i1,2)-l2,ng(i1,3)-l3)
                     t = fr*s(i2,u)*s(j2,u)
                     a(i,j,1) = a(i,j,1) + t*zfcos(y)
                     a(i,j,2) = a(i,j,2) - t*zfsin(y)
                  end do
               end if
            end do
         end do
      end do
      deallocate(s)

    end subroutine nlmtrx

    subroutine ovmtrx
      integer :: i,j,it,ia,u,v,n,ip,y,jn,j1,l1,l2,l3,i1,lmt1,lmt2
      integer :: i2,j2,m,jm
      real(kind=DP) :: fr,t
      real(kind=DP),pointer,dimension(:,:,:) :: a
      real(kind=DP),allocatable,dimension(:,:)   :: s
      integer,pointer,dimension(:,:)         :: ng

      a => hsmat
      ng => ngabc

      allocate(s(kg1_ext,nlmtt)); s = 0.0d0

      if(k_symmetry(ik) == GAMMA) then
         do j = 1, nlmtt
            ! --- finding it and lmt
            Finding: do u = 1, ntyp
               do ip = 1, ilmt(u)
                  if(lmtt(ip,u) == j) then
                     lmt1 = ip
                     it = u
                     exit Finding
                  end if
               end do
            end do Finding
            l1 = ltp(lmt1,it)
            s(1,j) = snl(1,j,iksnl)
            if(mod(l1,2) == 0) then
               do i = 2, iba(ik)
                  i1 = nbase(i,ik)
                  i2 = nbase_gamma(i,2)
                  s(i1,j) =  snl(i,j,iksnl)
                  s(i2,j) = -snl(i,j,iksnl)
               end do
            else
               do i = 2, iba(ik)
                  i1 = nbase(i,ik)
                  i2 = nbase_gamma(i,2)
                  s(i1,j) = snl(i,j,iksnl)
                  s(i2,j) = snl(i,j,iksnl)
               end do
            end if
         end do
         if(iprimatdiagon >= 3) then
            write(nfout,'(" --- ovmtrx ---")')
            write(nfout,'(" ik = ",i8)') ik
            do j = 1, nlmtt
               write(nfout,'(" s(1:*,",i3,")      : ",10f8.4)')j, (s(i,j),i=1,10)
            end do
         end if
      else
         do j = 1, nlmtt
            do i = 1, kg1
               s(i,j) = snl(i,j,iksnl)
            end do
         end do
      end if

      a = 0.d0
      do i = 1, iba2(ik)
         a(i,i,1) = 1.d0
      end do
      do ia = 1, natm
         it = ityp(ia)
         call calphase(ia)  ! -> zfcos, zfsin
         do ip = 1, n_non0_lmtxlmt(it)
            lmt1 = index_lmt1_lmt2(ip,it,1)
            u = lmtt(lmt1,it)
            lmt2 = index_lmt1_lmt2(ip,it,2)
            if(lmt2 <= lmt1) cycle
            v = lmtt(lmt2,it)
            if(ltp(lmt1,it) /= ltp(lmt2,it) .or. mtp(lmt1,it) /= mtp(lmt2,it) ) cycle
            fr = q(lmt1,lmt2,it)*iwei(ia)
            do j = 1, iba2(ik)
               jn = nbmat(j,ik)
               jm = nbmat2(j,ik)
               l1 = ng(jn,1); l2 = ng(jn,2); l3 = ng(jn,3)
               if(kimg == 1) then
                  do i = j, iba2(ik)
                     n = nbmat(i,ik)
                     m = nbmat2(i,ik)
                     y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                     a(i,j,1) = a(i,j,1)+fr*(s(m,u)*s(jm,v)+s(m,v)*s(jm,u))*zfcos(y)
                  end do
               else
                  do i = j, iba2(ik)
                     n = nbmat(i,ik)
                     m = nbmat2(i,ik)
                     y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                     t  = fr*(s(m,u)*s(jm,v)+s(m,v)*s(jm,u))
                     a(i,j,1) = a(i,j,1) + t*zfcos(y)
                     a(i,j,2) = a(i,j,2) - t*zfsin(y)
                  end do
               end if
            end do
         end do
         do lmt1 = 1, ilmt(it)
            u  = lmtt(lmt1,it)
            fr = q(lmt1,lmt1,it)*iwei(ia)
            do j = 1, iba2(ik)
               j1 = nbmat(j,ik)
               j2 = nbmat2(j,ik)
               l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)
               if(kimg == 1) then
                  do i = j, iba2(ik)
                     i1 = nbmat(i,ik)
                     i2 = nbmat2(i,ik)
                     ip = igpo(ng(i1,1)-l1,ng(i1,2)-l2,ng(i1,3)-l3)
                     a(i,j,1) = a(i,j,1)+fr*s(i2,u)*s(j2,u)*zfcos(ip)
                  end do
               else if(kimg == 2) then
                  do i = j, iba2(ik)
                     i1 = nbmat(i,ik)
                     i2 = nbmat2(i,ik)
                     ip = igpo(ng(i1,1)-l1,ng(i1,2)-l2,ng(i1,3)-l3)
                     t  = fr*s(i2,u)*s(j2,u)
                     a(i,j,1) = a(i,j,1) + t*zfcos(ip)
                     a(i,j,2) = a(i,j,2) - t*zfsin(ip)
                  end do
               end if
            end do
         end do
      end do

      if(kimg == 1) then
         do j = 2, iba2(ik)
            do i = 1, j-1
               a(i,j,1) = a(j,i,1)
            end do
         end do
      else
         do j = 2, iba2(ik)
            do i = 1, j-1
               a(i,j,1) =  a(j,i,1)
               a(i,j,2) = -a(j,i,2)
            end do
         end do
      end if
      deallocate(s)

    end subroutine ovmtrx

    subroutine cholde
! Choleski decomposition.
! ref: S.G.Louie,K.M.Ho,and M.L.Cohen, PRB19,(1979),p1774.
      integer       :: n,i,j
      real(kind=DP) :: t
      integer  :: imax, jmax
      integer, parameter :: NMATSIZE = 10

      linv = 0.d0
      linv(1,1,1) = 1.d0/sqrt(hsmat(1,1,1))
      if(kimg == 1) then
         do n = 1, iba2(ik)-1
            ldag = 0.d0
            do  i = 1, n
               do j = 1, i
                  ldag(i,1) = ldag(i,1) + linv(i,j,1)*hsmat(j,n+1,1)
               end do
            end do
            t = 0.d0
            do i = 1, n
               t = t + ldag(i,1)**2
            end do
            ldag(n+1,1) = sqrt(hsmat(n+1,n+1,1)-t)
            linv(n+1,n+1,1) = 1.d0/sqrt(hsmat(n+1,n+1,1)-t)
            do i = 1, n
               do j = 1, n
                  linv(n+1,i,1)=linv(n+1,i,1)-linv(n+1,n+1,1)*ldag(j,1)*linv(j,i,1)
               end do
            end do
         end do
      else if(kimg == 2) then
         do n = 1, iba2(ik)-1
            ldag = 0.d0
            do i = 1, n
               do j = 1, i
                  ldag(i,1) = ldag(i,1) &
                       & + linv(i,j,1)*hsmat(j,n+1,1)-linv(i,j,2)*hsmat(j,n+1,2)
                  ldag(i,2) = ldag(i,2) &
                       & + linv(i,j,1)*hsmat(j,n+1,2)+linv(i,j,2)*hsmat(j,n+1,1)
               end do
            end do
            t = 0.d0
            do i = 1, n
               t = t + ldag(i,1)**2 + ldag(i,2)**2
            end do
            ldag(n+1,1) = sqrt(hsmat(n+1,n+1,1)-t)
            linv(n+1,n+1,1) = 1.d0/sqrt(hsmat(n+1,n+1,1)-t)
            do i = 1, n
               do j = 1, n
                  linv(n+1,i,1) = linv(n+1,i,1) &
                       &       - linv(n+1,n+1,1)*(ldag(j,1)*linv(j,i,1) &
                       &                         +ldag(j,2)*linv(j,i,2))&
                       &       + linv(n+1,n+1,2)*(ldag(j,1)*linv(j,i,2)&
                       &                         -ldag(j,2)*linv(j,i,1) )
                  linv(n+1,i,2) = linv(n+1,i,2) &
                       &       - linv(n+1,n+1,1)*(ldag(j,1)*linv(j,i,2) &
                       &                         -ldag(j,2)*linv(j,i,1))&
                       &       - linv(n+1,n+1,2)*(ldag(j,1)*linv(j,i,1)&
                       &                         -ldag(j,2)*linv(j,i,2) )
               end do
            end do
         end do
      end if
      if(iprimatdiagon >= 2) then
         imax = NMATSIZE; if(imax > nmatsz) imax = nmatsz
         jmax = NMATSIZE; if(jmax > nmatsz) jmax = nmatsz
         do j = 1, jmax
            write(nfout,'(" linv(1:*,",i3,",1): ",10f8.4,99(/18x,10f8.4))') j &
                 & , (linv(i,j,1),i=1,imax)
         end do
         if(kimg == 2) then
            do j = 1, jmax
               write(nfout,'(" linv(1:*,",i3,",2): ",10f8.4,99(/18x,10f8.4))') j &
                    & , (linv(i,j,2),i=1,imax)
            end do
         end if
      end if
    end subroutine cholde

    subroutine amtrxkinetic
      integer :: i,j
      real(kind=DP) :: ga,gb,gc
      do i = 1, iba2(ik)
         j = nbmat(i,ik)
         ga = vkxyz(ik,1,BUCS) + ngabc(j,1)
         gb = vkxyz(ik,2,BUCS) + ngabc(j,2)
         gc = vkxyz(ik,3,BUCS) + ngabc(j,3)
         hsmat(i,i,1) = hsmat(i,i,1) + 0.5 &
              & * (ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
              & +  ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga)
      end do
    end subroutine amtrxkinetic

    subroutine prjvlhxc_l2vlhxc_t()
      integer :: i,irc
      vlhxc = 0.d0
      do irc = 1, kimg
         do i = ista_kngp, iend_kngp
            if(i > nmatsz2) cycle
            vlhxc(i,irc) = vlhxc_l(i,irc,ispin)
         end do
      end do
      if(npes >= 2) then
         call mpi_allreduce(MPI_IN_PLACE,vlhxc,nmatsz2*kimg,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
      end if
    end subroutine prjvlhxc_l2vlhxc_t

    subroutine amtrxlocalpt()
      integer :: j,j1,l1,l2,l3,i,n,ip
      real(kind=DP),pointer,dimension(:,:,:) :: a
      integer,      pointer,dimension(:,:)   :: ng
      a  => hsmat
      ng => ngabc

      do j = 1, iba2(ik)
         j1 = nbmat(j,ik)
         l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)
         if(kimg == 1) then
            do i = j, iba2(ik)
               n = nbmat(i,ik)
               ip = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
               a(i,j,1) = a(i,j,1) + vlhxc(ip,1)
            end do
         else
            do i = j, iba2(ik)
               n = nbmat(i,ik)
               ip = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
               a(i,j,1) = a(i,j,1) + vlhxc(ip,1)
               a(i,j,2) = a(i,j,2) + vlhxc(ip,2)
            end do
         end if
      end do
      if(kimg == 1) then
         do j = 2, iba2(ik)
            do i = 1, j-1
               a(i,j,1) = a(j,i,1)
            end do
         end do
      else
         do j = 2, iba2(ik)
            do i = 1, j-1
               a(i,j,1) =  a(j,i,1)
               a(i,j,2) = -a(j,i,2)
            end do
         end do
      end if
    end subroutine amtrxlocalpt

    subroutine lhl
      real(kind=DP),pointer,dimension(:,:,:)     :: h
      real(kind=DP),allocatable,dimension(:,:,:) :: ww
      integer :: i,j,k, M, N

      h => hsmat
      allocate(ww(nmatsz,nmatsz,kimg)); ww = 0.d0

      if(kimg == 1) then
#ifdef MATDIAGON_DGEMM
         M = iba2(ik); N = iba2(ik); K = iba2(ik)
         call DGEMM__('N','N',M,N,K, 1.d0, linv(1,1,1),nmatsz,h(1,1,1),nmatsz,0.d0, ww(1,1,1),nmatsz)
#else
         do i = 1, iba2(ik)
            do j = 1, iba2(ik)
               do k = 1, i
                  ww(i,j,1) = ww(i,j,1) + linv(i,k,1)*h(k,j,1)
               end do
            end do
         end do
#endif
#ifdef MATDIAGON_DGEMM
         call DGEMM__('N','T',M,N,K, 1.d0, ww(1,1,1),nmatsz,linv(1,1,1),nmatsz,0.d0, h(1,1,1),nmatsz)
#else
         h = 0.d0
         do i = 1, iba2(ik)
            do j = 1, iba2(ik)
               do k = 1, j
                  h(i,j,1) = h(i,j,1) + ww(i,k,1)*linv(j,k,1)
              end do
            end do
         end do
#endif
      else if(kimg == 2) then
#ifdef MATDIAGON_DGEMM
         M = iba2(ik); N = iba2(ik); K = iba2(ik)
         call DGEMM__('N','N',M,N,K, 1.d0, linv(1,1,1),nmatsz,h(1,1,1),nmatsz,0.d0, ww(1,1,1),nmatsz)
         call DGEMM__('N','N',M,N,K,-1.d0, linv(1,1,2),nmatsz,h(1,1,2),nmatsz,1.d0, ww(1,1,1),nmatsz)
         call DGEMM__('N','N',M,N,K, 1.d0, linv(1,1,1),nmatsz,h(1,1,2),nmatsz,0.d0, ww(1,1,2),nmatsz)
         call DGEMM__('N','N',M,N,K, 1.d0, linv(1,1,2),nmatsz,h(1,1,1),nmatsz,1.d0, ww(1,1,2),nmatsz)
#else
         do i = 1, iba2(ik)
            do j = 1, iba2(ik)
               do k = 1, i
                  ww(i,j,1) = ww(i,j,1)+linv(i,k,1)*h(k,j,1)-linv(i,k,2)*h(k,j,2)
                  ww(i,j,2) = ww(i,j,2)+linv(i,k,1)*h(k,j,2)+linv(i,k,2)*h(k,j,1)
               end do
            end do
         end do
#endif
#ifdef MATDIAGON_DGEMM
         M = iba2(ik); N = iba2(ik); K = iba2(ik)
         call DGEMM__('N','T',M,N,K, 1.d0, ww(1,1,1),nmatsz,linv(1,1,1),nmatsz,0.d0, h(1,1,1),nmatsz)
         call DGEMM__('N','T',M,N,K, 1.d0, ww(1,1,2),nmatsz,linv(1,1,2),nmatsz,1.d0, h(1,1,1),nmatsz)
         call DGEMM__('N','T',M,N,K, 1.d0, ww(1,1,2),nmatsz,linv(1,1,1),nmatsz,0.d0, h(1,1,2),nmatsz)
         call DGEMM__('N','T',M,N,K,-1.d0, ww(1,1,1),nmatsz,linv(1,1,2),nmatsz,1.d0, h(1,1,2),nmatsz)
#else
         h = 0.d0
         do i = 1, iba2(ik)
            do j = 1, iba2(ik)
               do k = 1, j
                  h(i,j,1) = h(i,j,1)+ww(i,k,1)*linv(j,k,1)+ww(i,k,2)*linv(j,k,2)
                  h(i,j,2) = h(i,j,2)+ww(i,k,2)*linv(j,k,1)-ww(i,k,1)*linv(j,k,2)
               end do
            end do
         end do
#endif
      end if
      deallocate(ww)
    end subroutine lhl

    subroutine lx()
#ifdef HIUX
! *** 'poption''s have been inserted by Dr. Kino (National Institute  ***
! *** for Materials Science, Japan).  12 July 2005                    ***
#endif
      integer :: ib, ib2, i, j
      integer :: i1
      if(kimg == 1) then
         do ib = ista_e, iend_e, istep_e
            ib2 = map_z(ib)
            do i = 1, iba2(ik)
               ldag(i,1) = zaj_mat(i,ib2,1)
            end do
            zaj_mat(:,ib2,1) = 0.d0
            do i = 1, iba2(ik)
               do j = 1, iba2(ik)
                  zaj_mat(i,ib2,1) = zaj_mat(i,ib2,1) + linv(j,i,1)*ldag(j,1)
               end do
            end do
         end do
      else if(kimg == 2) then
         do ib = ista_e, iend_e, istep_e
            ib2 = map_z(ib)
            if(iba2(ik) > nmatsz)  call phase_error_with_msg(nfout, ' iba2(ik) > nmatsz' &
                                                            , __LINE__, __FILE__)
            do i = 1, iba2(ik)
               ldag(i,1) = zaj_mat(i,ib2,1)
               ldag(i,2) = zaj_mat(i,ib2,2)
            end do
            zaj_mat(:,ib2,1:2) = 0.d0
            do i = 1, iba2(ik)
               do j = i, iba2(ik)
                  zaj_mat(i,ib2,1) = zaj_mat(i,ib2,1) &
                       &        + linv(j,i,1)*ldag(j,1) + linv(j,i,2)*ldag(j,2)
                  zaj_mat(i,ib2,2) = zaj_mat(i,ib2,2) &
                       &        + linv(j,i,1)*ldag(j,2) - linv(j,i,2)*ldag(j,1)
               end do
            end do
         end do
      end if
    end subroutine lx

    subroutine phase_mult()
      integer :: ib, ib2, ig
      real(kind=DP) :: wfr,wfi,sqrwf, pcos,psin
#ifndef _PHASE_MULT_0_
      integer :: ig1, ig2, ig_max
      complex(kind=CMPLDP) :: exp2theta, exptheta
      real(kind=DP) :: sqrwf0, wfr2, wfi2, sqrwfi, sqrwf1, sqrwf2,f &
           &          ,phase2r, phase2i, cos2theta, pcos0,psin0, norm
      real(kind=DP), parameter :: abswf_min = 1.d-16
#endif

!!$      do ig = 1, nmatsz
      if(kimg == 2) then
         do ib = ista_e, iend_e, istep_e
            ib2 = map_z(ib)
#ifdef _PHASE_MULT_0_
            wfr = zaj_mat(1,ib2,1); wfi = zaj_mat(1,ib2,2)
            if(abs(wfi*wfi) > 1.d-30) then
               sqrwf = wfr**2 + wfi**2
               f = dsqrt(sqrwf)
               pcos = -wfr/f; psin = wfi/f
               if(iprimatdiagon >= 2) then
                  write(nfout,'(4x," pcos, psin = ",2e14.6)') pcos, psin
               end if
               do ig = 1, iba2(ik)
                  wfr = zaj_mat(ig,ib2,1); wfi = zaj_mat(ig,ib2,2)
                  zaj_mat(ig,ib2,1) = pcos*wfr - psin*wfi
                  zaj_mat(ig,ib2,2) = pcos*wfi + psin*wfr
               end do
            end if
#else
            wfi = zaj_mat(1,ib2,2)
            sqrwfi = wfi*wfi
            if(sqrwfi > 1.d-35) then
               wfr = zaj_mat(1,ib2,1)
               sqrwf0 = wfr**2 + sqrwfi
!!$               if(sqrwf0 > abswf_min) then
               if(sqrwfi > abswf_min) then
                  ig_max = 1
               else
                  ig_max = 1
                  g_search: do ig = 2, kg2_gamma
                     ig1 = nbase_gamma(ig,1)
                     wfr = zaj_mat(ig1,ib2,1); wfi = zaj_mat(ig1,ib2,2)
                     sqrwfi = wfi*wfi
                     sqrwf = wfr**2 + sqrwfi
                     if(sqrwfi > abswf_min) then
                        ig_max = ig
                        exit g_search
                     else
                        if(sqrwf > sqrwf0) then
                           ig_max = ig
                           sqrwf0 = sqrwf
                        end if
                     end if
                  end do g_search
               end if
               if(ig_max == 1) then
                  wfr = zaj_mat(1,ib2,1); wfi = zaj_mat(1,ib2,2)
                  sqrwf0 = wfr**2 + wfi**2
                  f = dsqrt(sqrwf0)
                  pcos = -wfr/f; psin = wfi/f
                  if(iprimatdiagon >= 2) then
                     write(nfout,'(4x," pcos, psin = ",2e14.6, " ig_max = ",i5 )') pcos, psin,ig_max
                     write(nfout,'(4x," sqrwf = ",e14.6)') sqrwf0
                  end if
               else if(ig_max >= 2) then
                  wfr = zaj_mat(1,ib2,1); wfi = zaj_mat(1,ib2,2)
                  sqrwf0 = wfr**2 + wfi**2
                  f = dsqrt(sqrwf0)
                  pcos0 = -wfr/f; psin0 = wfi/f
                  if(iprimatdiagon >= 2) then
                     write(nfout,'(4x," pcos, psin = ",2e14.6, " ig =  0")') pcos0, psin0
                  end if
                  ig1 = nbase_gamma(ig_max,1)
                  ig2 = nbase_gamma(ig_max,2)
                  wfr  = zaj_mat(ig1,ib2,1);  wfi  = zaj_mat(ig1,ib2,2)
                  wfr2 = zaj_mat(ig2,ib2,1);  wfi2 = zaj_mat(ig2,ib2,2)
                  sqrwf1 = wfr**2 + wfi**2
                  sqrwf2 = wfr2**2 + wfi2**2
                  sqrwf = sqrt(sqrwf1*sqrwf2)
                  phase2r =  (wfr*wfr2 - wfi*wfi2)/sqrwf
                  phase2i = -(wfr*wfi2 + wfr2*wfi)/sqrwf
                  norm = sqrt(phase2r*phase2r + phase2i*phase2i)
                  phase2r = phase2r/norm
                  phase2i = phase2i/norm
                  exp2theta = dcmplx(phase2r,phase2i)
                  exptheta = sqrt(exp2theta)
                  pcos = real(exptheta)
                  psin = imag(exptheta)
                  if(iprimatdiagon >= 2) then
                     write(nfout,'(4x," pcos, psin = ",2e14.6, " ig_max = ",i5 )') pcos, psin,ig_max
                     write(nfout,'(4x," sqrwf1, sqrwf2, sqrwf = ",3e14.6)') sqrwf1, sqrwf2, sqrwf
                  end if
                  cos2theta = (wfr*wfr2 - wfi*wfi2)/sqrwf
                  cos2theta = cos2theta/norm
                  if(cos2theta > 1.0) cos2theta = 1.d0
                  if(cos2theta < -1.0) then
                     pcos = 0.d0
                     psin = 1.d0
                  else
                     pcos =  sqrt((1+cos2theta)/2.0)
                     psin =  sqrt((1-cos2theta)/2.0)
                  end if
                  if(pcos0 < 0.0) pcos = - pcos
                  if(psin0 < 0.0) psin = - psin
                  if(iprimatdiagon >= 2) then
                     write(nfout,'(4x," pcos, psin = ",2e14.6, " ig_max = ",i5 )') pcos, psin,ig_max
                  end if
               end if
               do ig = 1, iba2(ik)
                  wfr = zaj_mat(ig,ib2,1); wfi = zaj_mat(ig,ib2,2)
                  zaj_mat(ig,ib2,1) = pcos*wfr - psin*wfi
                  zaj_mat(ig,ib2,2) = pcos*wfi + psin*wfr
               end do
               zaj_mat(1,ib2,2) = 0.d0
            !$else
            !$   do ig = 2, kg2_gamma
            !$      ig1 = nbase_gamma(ig,1)
            !$      wfr = zaj_mat(ig1,ib2,1)
            !$      sqrwf = sqrwf + wfr**2
            !$   end do
            !$   sqrwf = 2.d0*sqrwf + zaj_mat(1,ib2,1)**2
            !$   sqrwf = 1.d0/sqrt(sqrwf)
            !$   do ig = 1, iba2(ik)
            !$      wfr = zaj_mat(ig,ib2,1)
            !$      zaj_mat(ig,ib2,1) = wfr*sqrwf
            !$      zaj_mat(ig,ib2,2) = 0.d0
            !$   end do
            end if
#endif
            if(iprimatdiagon >=2 ) then
               write(nfout,'(i4," Re ",5e14.6)') ib, (zaj_mat(ig,map_z(ib),1),ig=1,5)
               write(nfout,'(8x,5e14.6)') (zaj_mat(ig,map_z(ib),1),ig=6,40)
               write(nfout,'(i4," Im ",5e14.6)') ib, (zaj_mat(ig,map_z(ib),2),ig=1,5)
               write(nfout,'(8x,5e14.6)') (zaj_mat(ig,map_z(ib),2),ig=6,40)
            end if
         end do

         if(k_symmetry(ik) == GAMMA_base_symmetrization) then
            do ib = ista_e, iend_e, istep_e
               ib2 = map_z(ib)
               do ig = 2, kg2_gamma
                  ig1 = nbase_gamma(ig,1)
                  ig2 = nbase_gamma(ig,2)
                  wfr = zaj_mat(ig1,ib2,1); wfi = zaj_mat(ig1,ib2,2)
                  zaj_mat(ig2,ib2,1) = wfr
                  zaj_mat(ig2,ib2,2) = -wfi
               end do
            end do
         end if
      end if
      if(iprimatdiagon >= 2) write(nfout,'(" out of <<phase_mult>>")')
    end subroutine phase_mult

    subroutine cp_zaj_mat_to_zaj_l()
      integer :: ib,ig, i, l
      zaj_l(:,:,ik,:) = 0.d0
      if(k_symmetry(ik) == GAMMA) then
         do l = 1, kimg
            do ib = 1, np_e
               do i = 1, kg2_gamma
                  ig = nbase(i,ik)
                  zaj_l(i,ib,ik,l) = zaj_mat(ig,ib,l)
               end do
            end do
         end do
      else
         do l = 1, kimg
            do ib = 1, np_e
               do i = 1, iba2(ik)
                  ig = nbmat2(i,ik)
                  zaj_l(ig,ib,ik,l) = zaj_mat(i,ib,l)
               end do
            end do
         end do
      end if
#ifdef SAVE_FFT_TIMES
      if(sw_save_fft == ON) then
         do ib = 1, np_e
            status_saved_phifftr(ib,ik) = OLD
         end do
      end if
#endif
    end subroutine cp_zaj_mat_to_zaj_l

  end subroutine m_ESmat_solve_Hx_eq_eSx

! ============================== added by Tagami ======================== 11.0
  subroutine m_ESmat_solve_Hx_eq_eSx_noncl( nfout,iteration,iteration_electronic )
    integer, intent(in) :: nfout, iteration, iteration_electronic

    integer,parameter :: PRINT_LEVEL = 2
    integer   :: ik, iksnl, n=0
    integer   :: is1, is2
    integer :: is_tmp


    integer   :: id_sname = -1
    call tstatc0_begin('m_ESmat_solve_Hx_eq_eSx_noncl ',id_sname,1)

    if(iprimatdiagon >= 2) then
       write(nfout,'(" !! iteration = ",i6, " m_ESmat_solve_Hx_eq_eSx_noncl")') &
	&        iteration
    end if

    if(ekmode == ON ) then
       if(iteration_electronic == 1) then
          call m_pwBS_set_gmaxs(n_matrix_size,gmaxs_given)
       else
          call m_pwBS_set_gmaxs(n,gmax)
       end if
    else
       if(iteration == 1 .and. intzaj == by_matrix_diagon) then
          call m_pwBS_set_gmaxs(n_matrix_size,gmaxs_given)
       else if ( skip_alloc_phonon .and. iteration_electronic == 1 &
                    	         & .and. intzaj == by_matrix_diagon) then
          call m_pwBS_set_gmaxs(n_matrix_size,gmaxs_given)
       else
          call m_pwBS_set_gmaxs(n,gmax)
       end if
    end if

    call m_pwBS_alloc_nbmat_and_iba2()
    call m_pwBS_mat_for_each_WF() ! -> nbmat,nbmat2,iba2

    call m_pwBS_alloc_igpo()
    call m_pwBS_GminusGmapfunction()  ! -> igpo, nmatsz2

    if(iprimatdiagon >= 3) write(nfout,'(" nmatsz2 = ",i9)') nmatsz2
    if(iprimatdiagon >= PRINT_LEVEL) &
         & write(nfout,'(" nmatsz, kg1_ext = ", 2i9)') nmatsz, kg1_ext

    call alloc_matrices()
!
    call m_ES_MagMom_to_DensMat_vlhxcl( vlhxc_l, vlhxc_ssrep )
!
    do ik = 1, kv3, ndim_spinor

       if(map_k(ik) /= myrank_k) cycle
       if(iprimatdiagon >= PRINT_LEVEL) &
               & write(nfout,'(" iba2(",i3,") = ", i9)') ik, iba2(ik)

       iksnl = (ik-1)/ndim_spinor + 1
       hsmat = 0.d0

       if ( modnrm == EXECUT ) then
          if ( iprimatdiagon >=3 ) then
            write(nfout,'(" !! modnrm == EXECUT (m_ESmat_solve_Hx_eq_eSx_noncl)")')
          endif

          Do is1=1, ndim_spinor
             call ovmtrx_noncl( is1, is1 )
             if ( iprimatdiagon >= 3 ) then
                call wd_hsmat_noncl(" -- ovmtrx --", is1, is1 )
             end if
          End do

          call complete_hermite_matrix_noncl

          if ( sw_use_zhegvx == ON ) then
             hsmat_s = hsmat
          else
             call cholde_noncl()
             if ( iprimatdiagon >= 3 ) then
                Do is1=1, ndim_spinor
                   call wd_hsmat_noncl(" -- cholde --",is1,is1 )
                End do
             endif
          endif

       endif

       hsmat = 0.d0
       Do is1=1, ndim_spinor
          Do is2=1, ndim_spinor
             call nlmtrx_noncl( is1, is2 )           ! (nlmtrx) -> hsmat <- snl
             if ( iprimatdiagon >= 3 ) then
                call wd_hsmat_noncl(" -- nlmtrx --",is1,is2 )
             end if

             if ( is1==is2 ) then
                call amtrxkinetic_noncl( is1, is2 )     !  -> hsmat
                if ( iprimatdiagon >= 3 ) then
                   call wd_hsmat_noncl(" -- amtrxkinetic --", is1, is2 )
                endif
             endif

             is_tmp = ( is1 -1 )*ndim_spinor + is2
             call prjvlhxc_l2vlhxc_t_noncl( is_tmp )    ! vlhxc_l -> vlhxc
             call amtrxlocalpt_noncl( is1, is2 )                !  -> hsmat
             if ( iprimatdiagon >= 3 ) then
	        call wd_hsmat_noncl(" -- amtrxlocalpt --", is1, is2 )
             endif
          End do
       End do

       call complete_hermite_matrix_noncl

       if ( modnrm == EXECUT ) then
          if ( sw_use_zhegvx == OFF ) then
             call lhl_noncl() ! -> hsmat
             if ( iprimatdiagon >= 3 ) then
                Do is1=1, ndim_spinor
                   Do is2=1, ndim_spinor
                      call wd_hsmat_noncl(" -- lhl --", is1, is2 )
                   End do
                End do
             endif
          endif
       endif

       if ( modnrm == EXECUT .and. sw_use_zhegvx == ON ) then
          call solve_Hx_eq_ex_LAPACK_noncl2()
       else
          call solve_Hx_eq_ex_LAPACK_noncl()
       endif

       if ( iprimatdiagon >= 3 ) then
          Do is1=1, ndim_spinor
             call wd_hsmat_noncl(" -- after diagonalization --",is1,1 )
          End do
       endif

       if ( modnrm == EXECUT ) then
          if ( sw_use_zhegvx == OFF ) then
             call lx_noncl() ! -> zaj_mat
          endif
       endif

!       if ( k_symmetry(ik) == GAMMA .or. &
!	&     k_symmetry(ik) == GAMMA_base_symmetrization ) then
!	  call phase_mult_noncl()
!       endif

       call cp_zaj_mat_to_zaj_l_noncl()     ! zaj_mat -> zaj_l

    end do


    if ( iprimatdiagon >= 3 ) then
       do ik = 1, kv3, ndim_spinor
          if ( map_k(ik) /= myrank_k ) cycle
          call m_ES_wd_zaj_small_portion( nfout, ik, &
	&                         " (m_ESmat_solve_Hx_eq_eSx_noncl)",31 )
       end do
    end if

    call dealloc_matrices()
    call m_pwBS_dealloc_igpo()
    call m_pwBS_dealloc_nbmat_and_iba2()

    call m_ES_betar_dot_WFs(nfout) ! (fsrfsi) -> fsr_l, fsi_l

! ****************************************** This is important ******
    call m_ES_sort_eigen_vals_noncl()
! *******************************************************************

    call tstatc0_end(id_sname)

  contains

    subroutine wd_hsmat_noncl( str, is1, is2 )
      character(len=*),intent(in) :: str
      integer, intent(in) :: is1, is2

      integer :: i, imax, j, jmax
      integer :: ishift, jshift

      integer, parameter :: NMATSIZE = 10
      write(nfout,'(a30," ik= ",i5)') str,ik

      imax = min(NMATSIZE, nmatsz)
      jmax = min(NMATSIZE, nmatsz)
!
      ishift = ( is1 -1 )*iba2(ik)
      jshift = ( is2 -1 )*iba2(ik)
!
      write(nfout,*) 'istart = ', ishift +1
      do j = 1, jmax
         write(nfout,'(" hsmat(istart:*,",i3,",1): ",10f8.4,99(/19x,10f8.4))') j+jshift &
              & , (hsmat( i+ishift, j+jshift, 1 ),i=1,imax)
      end do
      if(kimg == 2) then
         do j = 1, jmax
            write(nfout,'(" hsmat(istart:*,",i3,",2): ",10f8.4,99(/19x,10f8.4))') j+jshift &
                 & , (hsmat( i+ishift,j+jshift,2 ),i=1,imax)
         end do
      end if
    end subroutine wd_hsmat_noncl

   ! [LAPACK] with LAPACK library
   !***** LAPACK [begin] ***** M.Okamoto
    subroutine solve_Hx_eq_ex_LAPACK_noncl()
      integer :: n, lda, ne, ne_found, il, iu, info, ipos, &
                 i, j, ib, lwork, k
      real(kind=DP) ::                     abstol, vl, vu, sum, dlamch
      integer,allocatable ::               ifail(:), iwork(:)
      real(kind=DP),allocatable ::         amat(:,:), evec(:), rvmat(:,:), rwork(:)
      complex(kind=CMPLDP),allocatable ::  cmat(:,:), cvmat(:,:), cwork(:)
      complex(kind=CMPLDP),parameter ::    IMG = (0.d0,1.d0)
      complex(kind=CMPLDP) :: ctemp(1)

      zaj_mat(:,:,:) = 0.d0
      lda = nmatsz *ndim_spinor

      vl = -1.d10 ; vu = 1.d10 ; il = 1 ; iu = neg

      ne = iu - il + 1
      !abstol = eps_solve_Hx_eq_ex
!!$      abstol = 1.d-15
      !abstol = -1.d0

      abstol = 2*dlamch('S')
      if ( kimg == 1 ) then
         n = iba2(ik) *ndim_spinor

         lwork = max(1,12*n)
        !+++++++++++++++++++++++++++++++++++++++++++++++++
         allocate(evec(n),iwork(5*n),ifail(n))
         allocate(amat(lda,n),rvmat(lda,ne),rwork(lwork)); amat = 0.d0
        !+++++++++++++++++++++++++++++++++++++++++++++++++
         do i = 1,n
            do j = 1,n
               amat(i,j) = hsmat(i,j,1)
            end do
         end do
         if ( iprimatdiagon>= 2 )  call check_LAPACKinputs("solve_Hx_eq_ex_LAPACK_noncl "  &
              &                                   ,nfout,n,lda,vl,vu,il,iu,lwork,abstol,RMat=amat)
         call dsyevx( 'V', 'I', 'U', n, amat, lda, vl, vu, il, iu, abstol, &
                       ne_found, evec, rvmat, lda, rwork, lwork, iwork, ifail, info )

         if(iprimatdiagon>=2) call  check_LAPACKoutputs("solve_Hx_eq_ex_LAPACK_noncol " &
              &                                   ,nfout,ik,n,ne,lda,evec,Rmat=rvmat)
         call check_LAPACKerror(ne_found,ne,info)
         do ib = ista_e, iend_e, istep_e
            j = map_z(ib)
            do i = 1,n
               zaj_mat(i,j,1) = rvmat(i,ib)
            end do
            eko_l(j,ik) = evec(ib)
         end do
        !+++++++++++++++++++++++++++++
         deallocate(amat,rvmat,rwork)
        !+++++++++++++++++++++++++++++
      else if (kimg == 2) then
         n = iba2(ik) *ndim_spinor

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         allocate(evec(n),iwork(5*n),ifail(n))
         allocate(cmat(lda,n),cvmat(lda,ne),rwork(7*n))
! -------------------------------- ktDEBUG --------
         evec = 0.d0; iwork = 0; ifail = 0; cmat = cmplx(0.d0,0.d0)
         cvmat = dcmplx(0.d0,0.d0); rwork = 0.d0;
! -------------------------------- ktDEBUG --------
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         do i = 1,n
           do j = 1,n
              cmat(i,j) = cmplx(hsmat(i,j,1),hsmat(i,j,2))
           end do
         end do

! -------------------------------- ktDEBUG --------
         call zheevx( 'V', 'I', 'U', n, cmat, lda, vl, vu, il, iu, abstol, &
                       ne_found, evec, cvmat, lda, ctemp, -1, rwork, iwork, &
                       ifail, info )
         lwork = int(real(ctemp(1)))
         allocate(cwork(lwork)); cwork = dcmplx(0.d0,0.d0)
! -------------------------------- ktDEBUG --------

         if(iprimatdiagon>= 2)  call check_LAPACKinputs("solve_Hx_eq_ex_LAPACK_noncl "  &
              &                                   ,nfout,n,lda,vl,vu,il,iu,lwork,abstol,CMat=cmat)

         call zheevx( 'V', 'I', 'U', n, cmat, lda, vl, vu, il, iu, abstol, &
                       ne_found, evec, cvmat, lda, cwork, lwork, rwork, iwork, &
	               ifail, info )

         if(iprimatdiagon>=2) call  check_LAPACKoutputs("solve_Hx_eq_ex_LAPACK_noncol " &
              &                                   ,nfout,ik,n,ne,lda,evec,Cmat=cvmat)

         call check_LAPACKerror(ne_found,ne,info)

         do ib = ista_e, iend_e, istep_e
            j = map_z(ib)
            do i = 1,n
               zaj_mat(i,j,1) = real(cvmat(i,ib))
               zaj_mat(i,j,2) = aimag(cvmat(i,ib))
            end do
            eko_l(j,ik) = evec(ib)
         end do
        !+++++++++++++++++++++++++++++++++++
         deallocate(cmat,cvmat,cwork,rwork)
        !+++++++++++++++++++++++++++++++++++
      end if

      if(iprimatdiagon >= 2) call wd_CalculatedEigenvalues(nfout,ik,n,ne,evec,"solve_Hx_eq_ex_noncl")

      sum = 0.d0
      do i = 1,ne
         sum = sum + abs(evec(i))
      end do
      if (sum < SmallestPositiveNumber*1.d5) then
         call phase_error_with_msg(nfout, ' ! illegal evec (solve_Hx_eq_ex_noncl)',__LINE__,__FILE__)
      endif
      if (iprimatdiagon >= 2) write(nfout,'(" -- sum = ",f20.8)') sum

     !+++++++++++++++++++++++++++++
      deallocate(evec,iwork,ifail)
     !+++++++++++++++++++++++++++++

    end subroutine solve_Hx_eq_ex_LAPACK_noncl

    subroutine solve_Hx_eq_ex_LAPACK_noncl2()
      integer :: n, lda, ne, ne_found, il, iu, info, ipos, &
                 i, j, ib, lwork, k
      real(kind=DP) :: abstol, vl, vu, sum, dlamch
      integer,allocatable :: ifail(:), iwork(:)
      real(kind=DP),allocatable :: amat(:,:), evec(:), rvmat(:,:), rwork(:)
      complex(kind=CMPLDP),allocatable :: cmat(:,:), cmat2(:,:), cvmat(:,:), cwork(:)

      complex(kind=CMPLDP),parameter :: IMG = (0.d0,1.d0)

      complex(kind=CMPLDP) :: ctemp(1)

      zaj_mat(:,:,:) = 0.d0
      lda = nmatsz *ndim_spinor

      vl = -1.d10 ; vu = 1.d10 ; il = 1 ; iu = neg
      ne = iu - il + 1

      abstol = 2*dlamch('S')

      if ( kimg == 1 ) then

      else if (kimg == 2) then
         n = iba2(ik) *ndim_spinor

        !+++++++++++++++++++++++++++++++++++
         allocate(evec(n),iwork(5*n),ifail(n))
         allocate(cmat(lda,n),cmat2(lda,n),cvmat(lda,ne),rwork(7*n))
         evec = 0.d0; iwork = 0; ifail = 0;
         cmat = dcmplx(0.d0,0.d0); cmat2 = dcmplx(0.d0,0.d0); cvmat = dcmplx(0.d0,0.d0); rwork = 0.d0;
        !+++++++++++++++++++++++++++++++++++

         do i = 1,n
           do j = 1,n
              cmat(i,j) = cmplx(hsmat(i,j,1),hsmat(i,j,2))
              cmat2(i,j) = cmplx(hsmat_s(i,j,1),hsmat_s(i,j,2))
           end do
         end do

         call zhegvx( 1, 'V', 'I', 'U', n, cmat, lda, cmat2, lda, &
              &       vl, vu, il, iu, abstol, &
                      ne_found, evec, cvmat, lda, ctemp, -1, rwork, iwork, &
                      ifail, info )
         lwork = int(real(ctemp(1)))
         allocate(cwork(lwork)); cwork = dcmplx(0.d0,0.d0)

         if(iprimatdiagon>=2) call check_LAPACKinputs("solve_Hx_eq_ex_LAPACK_noncl2 "  &
              &                                   ,nfout,n,lda,vl,vu,il,iu,lwork,abstol,CMat=cmat)

         call zhegvx( 1, 'V', 'I', 'U', n, cmat, lda, cmat2, lda, &
              &       vl, vu, il, iu, abstol, &
                      ne_found, evec, cvmat, lda, cwork, lwork, rwork, iwork, &
                      ifail, info )

         if(iprimatdiagon>=2) call check_LAPACKoutputs("solve_Hx_eq_ex_LAPACK_noncl2 " &
              &                                   ,nfout,ik,n,ne,lda,evec,Cmat=cvmat)

         call check_LAPACKerror(ne_found,ne,info)

         do ib = ista_e, iend_e, istep_e
            j = map_z(ib)
            do i = 1,n
               zaj_mat(i,j,1) = real(cvmat(i,ib))
               zaj_mat(i,j,2) = aimag(cvmat(i,ib))
            end do
            eko_l(j,ik) = evec(ib)
         end do
        !+++++++++++++++++++++++++++++++++++
         deallocate(cmat,cmat2,cvmat,cwork,rwork)
        !+++++++++++++++++++++++++++++++++++
      end if

      if(iprimatdiagon >= 2) call wd_CalculatedEigenvalues(nfout,ik,n,ne,evec,"solve_Hx_eq_ex_noncl2")

      sum = 0.d0
      do i = 1,ne
         sum = sum + abs(evec(i))
      end do
      if (sum < SmallestPositiveNumber*1.d5) then
         call phase_error_with_msg(nfout, ' ! illegal evec (solve_Hx_eq_ex_noncl)',__LINE__,__FILE__)
      endif
      if (iprimatdiagon >= 2) write(nfout,'(" -- sum = ",f20.8)') sum

     !+++++++++++++++++++++++++++++
      deallocate(evec,iwork,ifail)
     !+++++++++++++++++++++++++++++

    end subroutine solve_Hx_eq_ex_LAPACK_noncl2

    subroutine calphase(ia)
      integer, intent(in) :: ia
      integer       :: i
      real(kind=DP) :: ph
      do i = 1, nmatsz2
         ph = (pos(ia,1)*ngabc(i,1)+pos(ia,2)*ngabc(i,2)+pos(ia,3)*ngabc(i,3))*PAI2
         zfcos(i) = cos(ph)
         zfsin(i) = sin(ph)
      end do
      zfcos(0) = 0.d0; zfsin(0) = 0.d0
    end subroutine calphase

    subroutine nlmtrx_noncl( is1, is2 )
      integer, intent(in) :: is1, is2

      real(kind=DP),pointer,dimension(:,:,:) :: a
      real(kind=DP),allocatable,dimension(:,:)  :: s
      integer, pointer, dimension(:,:)       :: ng

      complex(kind=CMPLDP),parameter :: zi = (0.d0,1.d0)
      integer :: i,j,it,ia,u,v,n,lmt1,il1,lmt2,il2,j1,l1,l2,l3,ip,y,i1
      integer :: i2,j2,m
      integer :: ishift, jshift

      integer :: is_tmp

      complex(kind=CMPLDP) :: fr,f21,f12,t, fr_tmp1, fr_tmp2
      real(kind=DP) :: t_r, t_i

!
      ishift = ( is1 -1 )*iba2(ik)
      jshift = ( is2 -1 )*iba2(ik)
!
      is_tmp = ( is1 -1 )*ndim_spinor + is2
!
      a => hsmat
      ng => ngabc

      allocate(s(kg1_ext,nlmtt)); s = 0.d0

      if(k_symmetry(ik) == GAMMA) then
         if(iprimatdiagon >= 2) write(nfout,'(" k_symmetry(ik) == GAMMA <<nlmtrx>>")')
         if(iprimatdiagon >= 2) write(nfout,'(" kg2_gamma, iba(ik), iba2(ik) = ",3i8)') kg2_gamma,iba(ik),iba2(ik)
         do j = 1, nlmtt
            ! --- finding it and lmt
            Finding: do u = 1, ntyp
               do ip = 1, ilmt(u)
                  if(lmtt(ip,u) == j) then
                     lmt1 = ip
                     it = u
                     exit Finding
                  end if
               end do
            end do Finding
            il1 = ltp(lmt1,it)
            if(iprimatdiagon >= 2) write(nfout,'(" nlmtt, j, lmt1, it, il1 = ",5i8)') nlmtt,j,lmt1,it,il1
            s(1,j) = snl(1,j,iksnl)
            if(mod(il1,2) == 0) then
               do i = 2, kg2_gamma
                  i1 = nbase(i,ik)
                  i2 = nbase_gamma(i,2)
                  s(i1,j) = snl(i,j,iksnl)
                  s(i2,j) = -snl(i,j,iksnl)
               end do
            else
               do i = 2, kg2_gamma
                  i1 = nbase(i,ik)
                  i2 = nbase_gamma(i,2)
                  s(i1,j) = snl(i,j,iksnl)
                  s(i2,j) = snl(i,j,iksnl)
               end do
            end if
         end do
      else
         if(iprimatdiagon >= 2) write(nfout,'(" k_symmetry(ik) /= GAMMA <<nlmtrx>>")')
         do j = 1, nlmtt
            do i = 1, iba(ik)
               s(i,j) = snl(i,j,iksnl)
            end do
         end do
      end if

      do ia = 1, natm
         it = ityp(ia)
         call calphase(ia)  ! -> zfcos, zfsin

! ===============
         if ( modnrm == SKIP .and. SpinOrbit_mode == Neglected &
              &              .and. sw_hubbard == OFF ) then
            do ip = 1, n_non0_lmtxlmt(it)
                                             ! This counts only m1==m2, l1==l2
               lmt1  = index_lmt1_lmt2(ip,it,1)
               lmt2  = index_lmt1_lmt2(ip,it,2)
               u     = lmtt(lmt1,it)
               v     = lmtt(lmt2,it)
               il1   = ltp(lmt1,it)
               il2   = ltp(lmt2,it)

               if(lmt2 <= lmt1) cycle

               if(ivan(it) /= 1 .and. &
                    & (il1 /= il2 .or. mtp(lmt1,it) /= mtp(lmt2,it))) cycle

               if(mod(il1+il2,2) == 0) then

                  fr = dion_scr_noncl(lmt1,lmt2,is_tmp,ia)

                  if(il1==il2) then
                     f21 = fr *iwei(ia)

                     do j = 1, iba2(ik)
                        j1 = nbmat(j,ik)
                        j2 = nbmat2(j,ik)
                        l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)
                        if(kimg == 1) then

                           do i = j, iba2(ik)
                              n = nbmat(i,ik)
                              m = nbmat2(i,ik)
                              y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                              a( i+ishift,j+jshift,1 ) = a( i+ishift,j+jshift,1 ) &
                                   &                   + f21*( s(m,u)*s(j2,v) &
                                   &                          +s(m,v)*s(j2,u) ) &
                                   &                    *zfcos(y)
                           end do
                        else

                           do i = j, iba2(ik)
                              n = nbmat(i,ik)
                              m = nbmat2(i,ik)
                              y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                              t  = f21*( s(m,u)*s(j2,v)+s(m,v)*s(j2,u) )

                              t_r = real(t);  t_i = aimag(t)

                              a( i+ishift,j+jshift,1 ) = a( i+ishift,j+jshift,1 ) &
                                   &                   + t_r *zfcos(y) + t_i *zfsin(y)
                              a( i+ishift,j+jshift,2 ) = a( i+ishift,j+jshift,2 ) &
                                   &                   - t_r *zfsin(y) + t_i *zfcos(y)
                           end do
                        end if
                     end do
                  else
                     f21 = real(zi**(il2-il1)) *fr *iwei(ia)
                     f12 = real(zi**(il1-il2)) *fr *iwei(ia)

                     write(nfout,'(" zi**(il1-il2) = ",2f16.8)') zi**(il1-il2)

                     do j = 1, iba2(ik)
                        j1 = nbmat(j,ik)
                        j2 = nbmat2(j,ik)
                        l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)

                        if(kimg == 1) then

                           do i = j, iba2(ik)
                              n = nbmat(i,ik)
                              m = nbmat2(i,ik)
                              y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                              a( i+ishift, j+jshift, 1 ) = a( i+ishift, j+jshift,1 )&
                                   &                      +( f21*s(m,u)*s(j2,v) &
                                   &                        +f12*s(m,v)*s(j2,u) ) &
                                   &                       *zfcos(y)
                           end do
                        else

                           do i = j, iba2(ik)
                              n = nbmat(i,ik)
                              m = nbmat2(i,ik)
                              y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)

                              t  = f21*s(m,u)*s(j2,v)+f12*s(m,v)*s(j2,u)

                              t_r = real(t);  t_i = aimag(t)

                              a( i+ishift,j+jshift,1 ) = a( i+ishift,j+jshift,1) &
                                   &                   + t_r *zfcos(y) + t_i *zfsin(y)
                              a( i+ishift,j+jshift,2 ) = a( i+ishift,j+jshift,2 ) &
                                   &                   - t_r *zfsin(y) + t_i *zfcos(y)

                           end do
                        end if
                     end do
                  end if
               else if(mod(il1+il2,2) == 1) then

                  fr = dion_scr_noncl(lmt1,lmt2,is_tmp,ia)

                  f21 = aimag(zi**(il2-il1))*fr*iwei(ia)    !    f21 = -aimag(zi**(il2-il1))*fr*iwei(ia)
                  f12 = aimag(zi**(il1-il2))*fr*iwei(ia)    !    f12 = -aimag(zi**(il1-il2))*fr*iwei(ia)

                  do j = 1, iba2(ik)
                     j1 = nbmat(j,ik)
                     j2 = nbmat2(j,ik)
                     l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)
                     if ( kimg == 1 ) then

                        do i = j, iba2(ik)
                           n = nbmat(i,ik)
                           m = nbmat2(i,ik)
                           y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                           a( i+ishift,j+jshift,1 ) = a( i+ishift,j+jshift,1 )&
                                &                    + ( f21*s(m,u)*s(j2,v) &
                                &                       +f12*s(m,v)*s(j2,u) ) &
                                &                     *zfsin(y)
                        end do
                     else

                        do i = j, iba2(ik)
                           n = nbmat(i,ik)
                           m = nbmat2(i,ik)
                           y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                           t = f21*s(m,u)*s(j2,v)+f12*s(m,v)*s(j2,u)

                           t_r = real(t);  t_i = aimag(t)

                           a( i+ishift,j+jshift,1 ) = a( i+ishift,j+jshift,1 ) &
                                &                   + t_r *zfsin(y) - t_i *zfcos(y)
                           a( i+ishift,j+jshift,2 ) = a( i+ishift,j+jshift,2 ) &
                                &                   + t_r *zfcos(y) + t_i *zfsin(y)
                        end do
                     end if
                  end do
               end if
            end do
! ----------------------------------------------- when Spin-Orbit is included --
         else

            do lmt1=1, ilmt(it)
               Do lmt2=1, ilmt(it)

                  u     = lmtt(lmt1,it)
                  v     = lmtt(lmt2,it)
                  il1   = ltp(lmt1,it)
                  il2   = ltp(lmt2,it)

                  if (lmt2 <= lmt1 ) cycle

!                  if(ivan(it) /= 1 .and. (il1 /= il2) ) cycle
!                  if ( il1 /= il2 ) cycle

                  if(mod(il1+il2,2) == 0) then

                     fr_tmp1 = dion_scr_noncl(lmt1,lmt2,is_tmp,ia)
                     fr_tmp2 = dion_scr_noncl(lmt2,lmt1,is_tmp,ia)

                     if(il1==il2) then
                        f21 = fr_tmp1 *iwei(ia)
                        f12 = fr_tmp2 *iwei(ia)

                        do j = 1, iba2(ik)
                           j1 = nbmat(j,ik)
                           j2 = nbmat2(j,ik)
                           l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)
                           if(kimg == 1) then

                              do i = j, iba2(ik)
                                 n = nbmat(i,ik)
                                 m = nbmat2(i,ik)
                                 y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                                 a( i+ishift,j+jshift,1 ) = a( i+ishift,j+jshift,1 ) &
                                      &                   + ( f21 *s(m,u)*s(j2,v) &
                                      &                      +f12 *s(m,v)*s(j2,u) ) &
                                      &                    *zfcos(y)
                              end do
                           else

                              do i = j, iba2(ik)
                                 n = nbmat(i,ik)
                                 m = nbmat2(i,ik)
                                 y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)

                                 t  = f21*s(m,u)*s(j2,v) + f12*s(m,v)*s(j2,u)

                                 t_r = real(t);  t_i = aimag(t)

                                 a( i+ishift,j+jshift,1 ) = a( i+ishift,j+jshift,1 ) &
                                      &                   + t_r *zfcos(y) + t_i *zfsin(y)
                                 a( i+ishift,j+jshift,2 ) = a( i+ishift,j+jshift,2 ) &
                                      &                   - t_r *zfsin(y) + t_i *zfcos(y)
                              end do
                           end if
                        end do
                     else

                        f21 = real(zi**(il2-il1)) *fr_tmp1 *iwei(ia)
                        f12 = real(zi**(il1-il2)) *fr_tmp2 *iwei(ia)

!!!!                        write(nfout,'(" zi**(il1-il2) = ",2f16.8)') zi**(il1-il2)

                        do j = 1, iba2(ik)
                           j1 = nbmat(j,ik)
                           j2 = nbmat2(j,ik)
                           l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)

                           if(kimg == 1) then

                              do i = j, iba2(ik)
                                 n = nbmat(i,ik)
                                 m = nbmat2(i,ik)
                                 y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                                 a( i+ishift, j+jshift, 1 ) = a( i+ishift, j+jshift,1 )&
                                      &                      +( f21*s(m,u)*s(j2,v) &
                                      &                        +f12*s(m,v)*s(j2,u) ) &
                                      &                       *zfcos(y)
                              end do
                           else

                              do i = j, iba2(ik)
                                 n = nbmat(i,ik)
                                 m = nbmat2(i,ik)
                                 y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)

                                 t  = f21*s(m,u)*s(j2,v)+f12*s(m,v)*s(j2,u)

                                 t_r = real(t);  t_i = aimag(t)

                                 a( i+ishift,j+jshift,1 ) = a( i+ishift,j+jshift,1) &
                                      &                   + t_r *zfcos(y) + t_i *zfsin(y)
                                 a( i+ishift,j+jshift,2 ) = a( i+ishift,j+jshift,2 ) &
                                      &                   - t_r *zfsin(y) + t_i *zfcos(y)

                              end do
                           end if
                        end do
                     end if

                  else if(mod(il1+il2,2) == 1) then

                     fr_tmp1 = dion_scr_noncl(lmt1,lmt2,is_tmp,ia)
                     fr_tmp2 = dion_scr_noncl(lmt2,lmt1,is_tmp,ia)

!                     f21 = -aimag(zi**(il2-il1))*fr_tmp1*iwei(ia)
!                     f12 = -aimag(zi**(il1-il2))*fr_tmp2*iwei(ia)
                     f21 = aimag(zi**(il2-il1))*fr_tmp1*iwei(ia)
                     f12 = aimag(zi**(il1-il2))*fr_tmp2*iwei(ia)

                     do j = 1, iba2(ik)
                        j1 = nbmat(j,ik)
                        j2 = nbmat2(j,ik)
                        l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)
                        if ( kimg == 1 ) then

                           do i = j, iba2(ik)
                              n = nbmat(i,ik)
                              m = nbmat2(i,ik)
                              y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                              a( i+ishift,j+jshift,1 ) = a( i+ishift,j+jshift,1 )&
                                   &                    + ( f21*s(m,u)*s(j2,v) &
                                   &                       +f12*s(m,v)*s(j2,u) ) &
                                   &                     *zfsin(y)
                           end do
                        else

                           do i = j, iba2(ik)
                              n = nbmat(i,ik)
                              m = nbmat2(i,ik)
                              y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                              t = f21*s(m,u)*s(j2,v)+f12*s(m,v)*s(j2,u)

                              t_r = real(t);  t_i = aimag(t)

                              a( i+ishift,j+jshift,1 ) = a( i+ishift,j+jshift,1 ) &
                                   &                   + t_r *zfsin(y) - t_i *zfcos(y)
                              a( i+ishift,j+jshift,2 ) = a( i+ishift,j+jshift,2 ) &
                                   &                   + t_r *zfcos(y) + t_i *zfsin(y)
                           end do
                        end if
                     end do
                  end if
               end do
            end do
         endif
! ------------------------------------------ diagonal term ----
         do lmt1 = 1, ilmt(it)

            fr = ( dion_scr_noncl(lmt1,lmt1,is_tmp,ia) )*iwei(ia)

            u  = lmtt(lmt1,it)
            do j = 1, iba2(ik)
               j1 = nbmat(j,ik)
               j2 = nbmat2(j,ik)
               l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)
               if(kimg == 1) then

                  do i = j, iba2(ik)
                     i1 = nbmat(i,ik)
                     i2 = nbmat2(i,ik)
                     y = igpo(ng(i1,1)-l1,ng(i1,2)-l2,ng(i1,3)-l3)
                     a( i+ishift,j+jshift,1 ) = a( i+ishift,j+jshift,1 ) &
	&                         + fr*s(i2,u)*s(j2,u)*zfcos(y)
                  end do
               else if(kimg == 2) then

                  do i = j, iba2(ik)
                     i1 = nbmat(i,ik)
                     i2 = nbmat2(i,ik)
                     y = igpo(ng(i1,1)-l1,ng(i1,2)-l2,ng(i1,3)-l3)
                     t = fr*s(i2,u)*s(j2,u)

                     t_r = real(t);  t_i = aimag(t)

                     a( i+ishift,j+jshift,1 ) = a( i+ishift,j+jshift,1 ) &
	&                                       + t_r *zfcos(y) + t_i *zfsin(y)
                     a( i+ishift,j+jshift,2 ) = a( i+ishift,j+jshift,2 ) &
	&                                       - t_r *zfsin(y) + t_i *zfcos(y)
                  end do
               end if
            end do
         end do
      end do
      deallocate(s)

    end subroutine nlmtrx_noncl

    subroutine ovmtrx_noncl( is1, is2 )
      integer, intent(in) :: is1, is2

      integer :: i,j,it,ia,u,v,n,ip,y,jn,j1,l1,l2,l3,i1,lmt1,lmt2
      integer :: i2,j2,m,jm
      integer :: il1, il2

      integer :: ishift, jshift

      integer :: is_tmp

      complex(kind=CMPLDP) :: fr,t, fr_tmp1, fr_tmp2
      real(kind=DP) :: t_r, t_i

      real(kind=DP),pointer,dimension(:,:,:) :: a
      real(kind=DP),allocatable,dimension(:,:)   :: s
      integer,pointer,dimension(:,:)         :: ng
!
      ishift = ( is1 -1 )*iba2(ik)
      jshift = ( is2 -1 )*iba2(ik)

      is_tmp = ( is1 -1 )*ndim_spinor + is2

      a => hsmat;        ng => ngabc
      allocate(s(kg1_ext,nlmtt))

      if(k_symmetry(ik) == GAMMA) then
         do j = 1, nlmtt
            ! --- finding it and lmt
            Finding: do u = 1, ntyp
               do ip = 1, ilmt(u)
                  if(lmtt(ip,u) == j) then
                     lmt1 = ip
                     it = u
                     exit Finding
                  end if
               end do
            end do Finding
            l1 = ltp(lmt1,it)
            s(1,j) = snl(1,j,iksnl)
            if(mod(l1,2) == 0) then
               do i = 2, iba(ik)
                  i1 = nbase(i,ik)
                  i2 = nbase_gamma(i,2)
                  s(i1,j) =  snl(i,j,iksnl)
                  s(i2,j) = -snl(i,j,iksnl)
               end do
            else
               do i = 2, iba(ik)
                  i1 = nbase(i,ik)
                  i2 = nbase_gamma(i,2)
                  s(i1,j) = snl(i,j,iksnl)
                  s(i2,j) = snl(i,j,iksnl)
               end do
            end if
         end do
         if(iprimatdiagon >= 3) then
            write(nfout,'(" --- ovmtrx ---")')
            write(nfout,'(" ik = ",i8)') ik
            do j = 1, nlmtt
               write(nfout,'(" s(1:*,",i3,")      : ",10f8.4)')j, (s(i,j),i=1,10)
            end do
         end if
      else
         do j = 1, nlmtt
            do i = 1, kg1
               s(i,j) = snl(i,j,iksnl)
            end do
         end do
      end if

      do i = 1, iba2(ik)
         a( i+ishift,i+jshift,1 ) = 1.d0
      end do

      do ia = 1, natm
         it = ityp(ia)
         call calphase(ia)  ! -> zfcos, zfsin

! ===============
         if ( SpinOrbit_mode /= BuiltIn  ) then
            do ip = 1, n_non0_lmtxlmt(it)
               lmt1 = index_lmt1_lmt2(ip,it,1)
               u = lmtt(lmt1,it)
               lmt2 = index_lmt1_lmt2(ip,it,2)

               if(lmt2 <= lmt1) cycle
               v = lmtt(lmt2,it)

               if(ltp(lmt1,it) /= ltp(lmt2,it) .or. mtp(lmt1,it) /= mtp(lmt2,it) ) cycle

               fr = q_noncl( lmt1, lmt2, is_tmp, it ) *iwei(ia)

               do j = 1, iba2(ik)
                  jn = nbmat(j,ik)
                  jm = nbmat2(j,ik)
                  l1 = ng(jn,1); l2 = ng(jn,2); l3 = ng(jn,3)

                  if ( kimg == 1 ) then

                     do i = j, iba2(ik)
                        n = nbmat(i,ik)
                        m = nbmat2(i,ik)
                        y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                        a( i+ishift, j+jshift, 1 ) = a( i+ishift, j+jshift, 1 ) &
                             &                      + fr *( s(m,u)*s(jm,v) &
                             &                             +s(m,v)*s(jm,u) ) &
                             &                        *zfcos(y)
                     end do
                  else

                     do i = j, iba2(ik)
                        n = nbmat(i,ik)
                        m = nbmat2(i,ik)
                        y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)

                        t  = fr *( s(m,u)*s(jm,v)+s(m,v)*s(jm,u) )
                        t_r = real(t);  t_i = aimag(t)

                        a( i+ishift, j+jshift, 1 ) = a( i+ishift,j+jshift, 1 ) &
                             &                      + t_r *zfcos(y) + t_i *zfsin(y)
                        a( i+ishift, j+jshift, 2 ) = a( i+ishift,j+jshift, 2 ) &
                             &                      - t_r *zfsin(y) + t_i *zfcos(y)
                     end do
                  end if
               end do
            end do

         else
! -------------------------------------- when Spin-Orbit is included -- ( builtin )
!            do lmt1=1, nlmt
!               Do lmt2=1, nlmt
            do lmt1=1, ilmt(it)
               Do lmt2=1, ilmt(it)

                  u     = lmtt(lmt1,it)
                  v     = lmtt(lmt2,it)
                  il1   = ltp(lmt1,it)
                  il2   = ltp(lmt2,it)

                  if(lmt2 <= lmt1) cycle

                  v = lmtt(lmt2,it)

                  if ( ltp(lmt1,it) /= ltp(lmt2,it) ) cycle

                  fr_tmp1 = q_noncl( lmt1, lmt2, is_tmp, it ) *iwei(ia)
                  fr_tmp2 = q_noncl( lmt2, lmt1, is_tmp, it ) *iwei(ia)

                  do j = 1, iba2(ik)
                     jn = nbmat(j,ik)
                     jm = nbmat2(j,ik)
                     l1 = ng(jn,1); l2 = ng(jn,2); l3 = ng(jn,3)

                     if ( kimg == 1 ) then

                        do i = j, iba2(ik)
                           n = nbmat(i,ik)
                           m = nbmat2(i,ik)
                           y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                           a( i+ishift, j+jshift, 1 ) = a( i+ishift, j+jshift, 1 ) &
                                &                      + ( fr_tmp1*s(m,u)*s(jm,v) &
                                &                       +fr_tmp2*s(m,v)*s(jm,u) ) &
                                &                       *zfcos(y)
                        end do
                     else

                        do i = j, iba2(ik)
                           n = nbmat(i,ik)
                           m = nbmat2(i,ik)
                           y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)

                           t  = fr_tmp1*s(m,u)*s(jm,v) +fr_tmp2*s(m,v)*s(jm,u)

                           t_r = real(t);  t_i = aimag(t)

                           a( i+ishift, j+jshift, 1 ) = a( i+ishift,j+jshift, 1 ) &
                                &                      + t_r *zfcos(y) + t_i *zfsin(y)
                           a( i+ishift, j+jshift, 2 ) = a( i+ishift,j+jshift, 2 ) &
                                &                      - t_r *zfsin(y) + t_i *zfcos(y)
                        end do
                     end if
                  end do
               End do
            End do

         endif
! ------------------------------------------ diagonal term ----
         do lmt1 = 1, ilmt(it)
            u  = lmtt(lmt1,it)

            fr = q_noncl(lmt1,lmt1,is_tmp,it) *iwei(ia)

            do j = 1, iba2(ik)
               j1 = nbmat(j,ik)
               j2 = nbmat2(j,ik)
               l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)

               if ( kimg == 1 ) then

                  do i = j, iba2(ik)
                     i1 = nbmat(i,ik)
                     i2 = nbmat2(i,ik)
                     ip = igpo(ng(i1,1)-l1,ng(i1,2)-l2,ng(i1,3)-l3)
                     a( i+ishift,j+jshift,1 ) = a( i+ishift,j+jshift,1 ) &
	&                         + fr *s(i2,u)*s(j2,u)*zfcos(ip)
                  end do
               else if (kimg == 2) then

                  do i = j, iba2(ik)
                     i1 = nbmat(i,ik)
                     i2 = nbmat2(i,ik)
                     ip = igpo(ng(i1,1)-l1,ng(i1,2)-l2,ng(i1,3)-l3)

                     t  = fr *s(i2,u)*s(j2,u)
                     t_r = real(t);  t_i = aimag(t)

                     a( i+ishift,j+jshift,1 ) = a( i+ishift,j+jshift,1 ) &
	&                         + t_r *zfcos(ip) + t_i *zfsin(ip)
                     a( i+ishift,j+jshift,2 ) = a( i+ishift,j+jshift,2 ) &
	&                         - t_r *zfsin(ip) + t_i *zfcos(ip)
                  end do
               end if
            end do
         end do
      end do

!      if(kimg == 1) then
!         do j = 2, iba2(ik)
!            do i = 1, j-1
!               a(i+isfhit,j+jshift,1) = a(j+jshift,i+isfhit,1)
!            end do
!         end do
!      else
!         do j = 2, iba2(ik)
!            do i = 1, j-1
!               a(i+isfhit,j+jsfhit,1) =  a(j+jshift,i+ishift,1)
!               a(i+isfhit,j+jshift,2) = -a(j+jshift,i+isfhit,2)
!            end do
!         end do
!      end if

      deallocate(s)

    end subroutine ovmtrx_noncl

    subroutine cholde_noncl
! Choleski decomposition.
! ref: S.G.Louie,K.M.Ho,and M.L.Cohen, PRB19,(1979),p1774.
      integer       :: n,i,j
      real(kind=DP) :: t
      integer  :: imax, jmax
      integer, parameter :: NMATSIZE = 10

      linv = 0.d0
      linv(1,1,1) = 1.d0/sqrt(hsmat(1,1,1))
      if ( kimg == 1 ) then
         do n = 1, iba2(ik)*ndim_spinor -1
            ldag = 0.d0
            do  i = 1, n
               do j = 1, i
                  ldag(i,1) = ldag(i,1) + linv(i,j,1)*hsmat(j,n+1,1)
               end do
            end do
            t = 0.d0
            do i = 1, n
               t = t + ldag(i,1)**2
            end do
            ldag(n+1,1) = sqrt(hsmat(n+1,n+1,1)-t)
            linv(n+1,n+1,1) = 1.d0/sqrt(hsmat(n+1,n+1,1)-t)
            do i = 1, n
               do j = 1, n
                  linv(n+1,i,1)=linv(n+1,i,1)-linv(n+1,n+1,1)*ldag(j,1)*linv(j,i,1)
               end do
            end do
         end do
      else if(kimg == 2) then

         do n = 1, iba2(ik)*ndim_spinor -1
            ldag = 0.d0
            do i = 1, n
               do j = 1, i
                  ldag(i,1) = ldag(i,1) &
                       & + linv(i,j,1)*hsmat(j,n+1,1)-linv(i,j,2)*hsmat(j,n+1,2)
                  ldag(i,2) = ldag(i,2) &
                       & + linv(i,j,1)*hsmat(j,n+1,2)+linv(i,j,2)*hsmat(j,n+1,1)
               end do
            end do
            t = 0.d0
            do i = 1, n
               t = t + ldag(i,1)**2 + ldag(i,2)**2
            end do

!
            ldag(n+1,1) = sqrt(hsmat(n+1,n+1,1)-t)
            linv(n+1,n+1,1) = 1.d0/sqrt(hsmat(n+1,n+1,1)-t)
            do i = 1, n
               do j = 1, n
                  linv(n+1,i,1) = linv(n+1,i,1) &
                       &       - linv(n+1,n+1,1)*(ldag(j,1)*linv(j,i,1) &
                       &                         +ldag(j,2)*linv(j,i,2))&
                       &       + linv(n+1,n+1,2)*(ldag(j,1)*linv(j,i,2)&
                       &                         -ldag(j,2)*linv(j,i,1) )
                  linv(n+1,i,2) = linv(n+1,i,2) &
                       &       - linv(n+1,n+1,1)*(ldag(j,1)*linv(j,i,2) &
                       &                         -ldag(j,2)*linv(j,i,1))&
                       &       - linv(n+1,n+1,2)*(ldag(j,1)*linv(j,i,1)&
                       &                         -ldag(j,2)*linv(j,i,2) )
               end do
            end do
         end do
      end if

      if(iprimatdiagon >= 2) then
         imax = NMATSIZE; if(imax > nmatsz) imax = nmatsz
         jmax = NMATSIZE; if(jmax > nmatsz) jmax = nmatsz

         do j = 1, jmax
            write(nfout,'(" linv(1:*,",i3,",1): ",10f8.4,99(/18x,10f8.4))') j &
                 & , (linv(i,j,1),i=1,imax)
         end do

         if(kimg == 2) then
            do j = 1, jmax
               write(nfout,'(" linv(1:*,",i3,",2): ",10f8.4,99(/18x,10f8.4))') j &
                    & , (linv(i,j,2),i=1,imax)
            end do
         end if

      end if
    end subroutine cholde_noncl

    subroutine amtrxkinetic_noncl( is1, is2 )
      integer, intent(in) :: is1, is2

      integer :: i,j
      integer :: ishift, jshift
      real(kind=DP) :: ga,gb,gc
!
      ishift = ( is1 -1 )*iba2(ik)
      jshift = ( is2 -1 )*iba2(ik)
!
      do i = 1, iba2(ik)
         j = nbmat(i,ik)
         ga = vkxyz(ik,1,BUCS) + ngabc(j,1)
         gb = vkxyz(ik,2,BUCS) + ngabc(j,2)
         gc = vkxyz(ik,3,BUCS) + ngabc(j,3)
         hsmat( i+ishift,i+jshift,1 ) = hsmat( i+ishift,i+jshift, 1 ) &
	      & + 0.5d0 * (ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
              &           +ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga )
      end do
    end subroutine amtrxkinetic_noncl

    subroutine prjvlhxc_l2vlhxc_t_noncl( is )
      integer, intent(in) :: is

      integer :: i,irc
      vlhxc = 0.d0
      do irc = 1, kimg
         do i = ista_kngp, iend_kngp
            if(i > nmatsz2) cycle
            vlhxc(i,irc) = vlhxc_ssrep(i,irc,is)
         end do
      end do
      if(npes >= 2) then
         call mpi_allreduce(MPI_IN_PLACE,vlhxc,nmatsz2*kimg,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
      end if
    end subroutine prjvlhxc_l2vlhxc_t_noncl

    subroutine amtrxlocalpt_noncl( is1, is2 )
      integer, intent(in) :: is1, is2

      integer :: j,j1,l1,l2,l3,i,n,ip
      integer :: ishift, jshift
      real(kind=DP),pointer,dimension(:,:,:) :: a
      integer,      pointer,dimension(:,:)   :: ng
!
      ishift = ( is1 -1 )*iba2(ik)
      jshift = ( is2 -1 )*iba2(ik)
!
      a  => hsmat
      ng => ngabc

      do j = 1, iba2(ik)
         j1 = nbmat(j,ik)
         l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)
         if(kimg == 1) then
!!!!!!!            do i = j, iba2(ik)
            do i = 1, iba2(ik)
               n = nbmat(i,ik)
               ip = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
               a( i+ishift,j+jshift,1 ) = a( i+ishift,j+jshift,1 ) + vlhxc(ip,1)
            end do
         else
!!!!!!!            do i = j, iba2(ik)
            do i = 1, iba2(ik)
               n = nbmat(i,ik)
               ip = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
               a( i+ishift,j+jshift,1 ) = a( i+ishift,j+jshift,1 ) + vlhxc(ip,1)
               a( i+ishift,j+jshift,2 ) = a( i+ishift,j+jshift,2 ) + vlhxc(ip,2)
            end do
         end if
      end do
! ============= dangerous -
!      if(kimg == 1) then
!         do j = 2, iba2(ik)
!            do i = 1, j-1
!               a( i+isfhit,j+jshift,1)  = a( j+jshift,i+isfhit,1 )
!!            end do
!         end do
!      else
!         do j = 2, iba2(ik)
!            do i = 1, j-1
!               a( i+ishift,j+jshift,1 ) =  a( j+jshift,i+ishift,1 )
!               a( i+ishift,j+jshift,2 ) = -a( j+jshift,i+ishift,2 )
!            end do
!         end do
!      end if
! ======================
    end subroutine amtrxlocalpt_noncl

    subroutine complete_hermite_matrix_noncl
      integer :: is1, is2, i, j
      integer :: ishift, jshift
      real(kind=DP),pointer,dimension(:,:,:) :: a

      a  => hsmat

      Do is1=1, ndim_spinor
         Do is2=1, ndim_spinor
            ishift = ( is1 -1 )*iba2(ik)
            jshift = ( is2 -1 )*iba2(ik)

!            istmp1 = ( is1 -1 )*ndim_spinor + is2
!            istmp2 = ( is2 -1 )*ndim_spinor + is1
!
            Do j=2, iba2(ik)
               Do i=1, j-1
                  a( i+ishift,j+jshift,1 ) =  a( j+jshift,i+ishift,1 )
                  a( i+ishift,j+jshift,2 ) = -a( j+jshift,i+ishift,2 )
               End do
            End do

         End do
      End do

    end subroutine complete_hermite_matrix_noncl

    subroutine lhl_noncl
      real(kind=DP),pointer,dimension(:,:,:)     :: h
      real(kind=DP),allocatable,dimension(:,:,:) :: ww
      integer :: i,j,k, M, N

      h => hsmat
      allocate( ww( nmatsz*ndim_spinor, nmatsz*ndim_spinor, kimg ) ); ww = 0.d0

      if (kimg == 1) then
#ifdef MATDIAGON_DGEMM
         M = iba2(ik) *ndim_spinor;
	 N = iba2(ik) *ndim_spinor;
	 K = iba2(ik) *ndim_spinor
         call DGEMM__(' N', 'N', M, N, K, 1.d0, linv(1,1,1), nmatsz*ndim_spinor, &
	&              h(1,1,1), nmatsz*ndim_spinor, 0.d0, &
	&              ww(1,1,1), nmatsz*ndim_spinor )
#else
         do i = 1, iba2(ik)*ndim_spinor
            do j = 1, iba2(ik)*ndim_spinor
               do k = 1, i
                  ww(i,j,1) = ww(i,j,1) + linv(i,k,1)*h(k,j,1)
               end do
            end do
         end do
#endif
#ifdef MATDIAGON_DGEMM
         call DGEMM__( 'N','T', M, N, K, 1.d0, ww(1,1,1), nmatsz*ndim_spinor, &
	&              linv(1,1,1), nmatsz*ndim_spinor, 0.d0, &
	&              h(1,1,1), nmatsz*ndim_spinor )
#else
         h = 0.d0
         do i = 1, iba2(ik) *ndim_spinor
            do j = 1, iba2(ik) *ndim_spinor
               do k = 1, j
                  h(i,j,1) = h(i,j,1) + ww(i,k,1)*linv(j,k,1)
              end do
            end do
         end do
#endif
      else if (kimg == 2) then
#ifdef MATDIAGON_DGEMM
         M = iba2(ik) *ndim_spinor
	 N = iba2(ik) *ndim_spinor
	 K = iba2(ik) *ndim_spinor
         call DGEMM__( 'N','N', M, N, K, 1.d0, linv(1,1,1), nmatsz*ndim_spinor,&
	&               h(1,1,1), nmatsz*ndim_spinor, 0.d0, &
	&               ww(1,1,1), nmatsz*ndim_spinor )
         call DGEMM__( 'N','N', M, N, K,-1.d0, linv(1,1,2), nmatsz*ndim_spinor,&
	&               h(1,1,2), nmatsz*ndim_spinor, 1.d0, &
	&               ww(1,1,1), nmatsz*ndim_spinor )
         call DGEMM__( 'N','N', M, N, K, 1.d0, linv(1,1,1), nmatsz*ndim_spinor,&
	&               h(1,1,2), nmatsz*ndim_spinor, 0.d0, &
	&               ww(1,1,2), nmatsz*ndim_spinor )
         call DGEMM__( 'N','N', M, N, K, 1.d0, linv(1,1,2), nmatsz*ndim_spinor,&
	&               h(1,1,1), nmatsz*ndim_spinor, 1.d0, &
	&               ww(1,1,2), nmatsz*ndim_spinor )
#else
         do i = 1, iba2(ik) *ndim_spinor
            do j = 1, iba2(ik) *ndim_spinor
               do k = 1, i
                  ww(i,j,1) = ww(i,j,1)+linv(i,k,1)*h(k,j,1)-linv(i,k,2)*h(k,j,2)
                  ww(i,j,2) = ww(i,j,2)+linv(i,k,1)*h(k,j,2)+linv(i,k,2)*h(k,j,1)
               end do
            end do
         end do
#endif
#ifdef MATDIAGON_DGEMM
         M = iba2(ik) *ndim_spinor;
	 N = iba2(ik) *ndim_spinor;
	 K = iba2(ik) *ndim_spinor
         call DGEMM__( 'N','T', M, N, K, 1.d0, ww(1,1,1), nmatsz*ndim_spinor, &
	&              linv(1,1,1), nmatsz*ndim_spinor, 0.d0, &
	&              h(1,1,1), nmatsz*ndim_spinor )
         call DGEMM__( 'N','T', M, N, K, 1.d0, ww(1,1,2), nmatsz*ndim_spinor, &
	&              linv(1,1,2), nmatsz*ndim_spinor, 1.d0, &
	&              h(1,1,1), nmatsz*ndim_spinor )
         call DGEMM__( 'N','T', M, N, K, 1.d0, ww(1,1,2), nmatsz*ndim_spinor, &
	&              linv(1,1,1), nmatsz*ndim_spinor, 0.d0, &
	&              h(1,1,2), nmatsz*ndim_spinor )
         call DGEMM__( 'N','T', M, N, K,-1.d0, ww(1,1,1), nmatsz*ndim_spinor, &
	&              linv(1,1,2), nmatsz*ndim_spinor, 1.d0, &
	&              h(1,1,2), nmatsz*ndim_spinor )
#else
         h = 0.d0
         do i = 1, iba2(ik) *ndim_spinor
            do j = 1, iba2(ik) *ndim_spinor
               do k = 1, j
                  h(i,j,1) = h(i,j,1)+ww(i,k,1)*linv(j,k,1)+ww(i,k,2)*linv(j,k,2)
                  h(i,j,2) = h(i,j,2)+ww(i,k,2)*linv(j,k,1)-ww(i,k,1)*linv(j,k,2)
               end do
            end do
         end do
#endif
      end if
      deallocate(ww)
    end subroutine lhl_noncl

    subroutine lx_noncl()
#ifdef HIUX
! *** 'poption''s have been inserted by Dr. Kino (National Institute  ***
! *** for Materials Science, Japan).  12 July 2005                    ***
#endif
      integer :: ib, ib2, i, j
      if(kimg == 1) then
         do ib = ista_e, iend_e, istep_e
            ib2 = map_z(ib)
            do i = 1, iba2(ik) *ndim_spinor
               ldag(i,1) = zaj_mat(i,ib2,1)
            end do
            zaj_mat(:,ib2,1) = 0.d0
            do i = 1, iba2(ik) *ndim_spinor
               do j = 1, iba2(ik) *ndim_spinor
                  zaj_mat(i,ib2,1) = zaj_mat(i,ib2,1) + linv(j,i,1)*ldag(j,1)
               end do
            end do
         end do
      else if(kimg == 2) then
         do ib = ista_e, iend_e, istep_e
            ib2 = map_z(ib)
            if(iba2(ik) > nmatsz)  call  phase_error_with_msg(nfout, ' iba2(ik) > nmatsz',__LINE__,__FILE__)
            do i = 1, iba2(ik) *ndim_spinor
               ldag(i,1) = zaj_mat(i,ib2,1)
               ldag(i,2) = zaj_mat(i,ib2,2)
            end do
            zaj_mat(:,ib2,1:2) = 0.d0
            do i = 1, iba2(ik) *ndim_spinor
               do j = i, iba2(ik)*ndim_spinor
                  zaj_mat(i,ib2,1) = zaj_mat(i,ib2,1) &
                       &        + linv(j,i,1)*ldag(j,1) + linv(j,i,2)*ldag(j,2)
                  zaj_mat(i,ib2,2) = zaj_mat(i,ib2,2) &
                       &        + linv(j,i,1)*ldag(j,2) - linv(j,i,2)*ldag(j,1)
               end do
            end do
         end do
      end if
    end subroutine lx_noncl

    subroutine phase_mult_noncl()
      integer :: ib, ib2, ig
      real(kind=DP) :: wfr, wfi, sqrwf, pcos,psin
#ifndef _PHASE_MULT_0_
      integer :: ig1, ig2, ig_max
      complex(kind=CMPLDP) :: exp2theta, exptheta
      real(kind=DP) :: sqrwf0, wfr2, wfi2, sqrwfi, sqrwf1, sqrwf2,f &
           &          ,phase2r, phase2i, cos2theta, pcos0,psin0, norm
      real(kind=DP), parameter :: abswf_min = 1.d-16
#endif

      integer :: is, igshift

!!$      do ig = 1, nmatsz
      if(kimg == 2) then
         do ib = ista_e, iend_e, istep_e
            ib2 = map_z(ib)
#ifdef _PHASE_MULT_0_

            Do is=1, ndim_spinor
              igshift = ( is-1 )*iba2(ik)
              wfr = zaj_mat( igshift+1,ib2,1 ); wfi = zaj_mat( igshift+1,ib2,2 )

              if (abs( wfi*wfi ) > 1.d-30) then
                 sqrwf = wfr1**2 + wfi**2
                 f = dsqrt(sqrwf)
                 pcos = -wfr/f; psin = wfi/f
                 if(iprimatdiagon >= 2) then
                    write(nfout,'(4x," pcos, psin = ",2e14.6)') pcos, psin
                 end if
                 do ig = 1, iba2(ik)
                    wfr = zaj_mat( ig+igshift,ib2,1 )
	            wfi = zaj_mat( ig+igshift,ib2,2 )
                    zaj_mat( ig+igshift,ib2,1 ) = pcos*wfr - psin*wfi
                    zaj_mat( ig+igshift,ib2,2 ) = pcos*wfi + psin*wfr
                 end do

              end if
            end do
#else
            Do is=1, ndim_spinor
              igshift = ( is-1 )*iba2(ik)
	      wfi = zaj_mat( igshift+1,ib2,2 )
              sqrwfi = wfi*wfi

              if (sqrwfi > 1.d-35) then
	         wfr = zaj_mat( igshift+1,ib2,1 );
                 sqrwf0 = wfr**2 + sqrwfi

                 if (sqrwfi > abswf_min) then
                    ig_max = 1
                 else
                    ig_max = 1
                    g_search: do ig = 2, kg2_gamma
                       ig1 = nbase_gamma(ig,1)
                       wfr = zaj_mat( ig1+igshift,ib2,1 );
                       wfi = zaj_mat( ig1+igshift,ib2,2 )

                       sqrwfi = wfi*wfi
                       sqrwf = wfr**2 + sqrwfi

                       if(sqrwfi > abswf_min) then
                          ig_max = ig
                          exit g_search
                       else
                          if(sqrwf > sqrwf0) then
                            ig_max = ig
                            sqrwf0 = sqrwf
                          end if
                       end if
                    end do g_search
                 end if

                 if (ig_max == 1) then
                    wfr = zaj_mat( igshift+1,ib2,1 ); wfi = zaj_mat( igshift+1,ib2,2 )
                    sqrwf0 = wfr**2 + wfi**2
                    f = dsqrt(sqrwf0)
                    pcos = -wfr/f; psin = wfi/f
                    if (iprimatdiagon >= 2 ) then
                        write(nfout,'(4x," pcos, psin = ",2e14.6, " ig_max = ",i5 )') &
	&                        pcos, psin,ig_max
                        write(nfout,'(4x," sqrwf = ",e14.6)') sqrwf0
                    end if
                 else if(ig_max >= 2) then
                    wfr = zaj_mat( igshift+1,ib2,1 ); wfi = zaj_mat( igshift+1,ib2,2 )
                    sqrwf0 = wfr**2 + wfi**2
                    f = dsqrt(sqrwf0)
                    pcos0 = -wfr/f; psin0 = wfi/f
                    if (iprimatdiagon >= 2) then
                       write(nfout,'(4x," pcos, psin = ",2e14.6, " ig =  0")') &
	&                       pcos0, psin0
                    end if

                    ig1 = nbase_gamma(ig_max,1);  ig2 = nbase_gamma(ig_max,2)

                    wfr  = zaj_mat( ig1+igshift,ib2,1 )
	            wfi  = zaj_mat( ig1+igshift,ib2,2 )
                    wfr2 = zaj_mat( ig2+igshift,ib2,1 )
	            wfi2 = zaj_mat( ig2+igshift,ib2,2 )

                    sqrwf1 = wfr**2 + wfi**2
                    sqrwf2 = wfr2**2 + wfi2**2
                    sqrwf = sqrt(sqrwf1*sqrwf2)

                    phase2r =  (wfr*wfr2 - wfi*wfi2)/sqrwf
                    phase2i = -(wfr*wfi2 + wfr2*wfi)/sqrwf
                    norm = sqrt(phase2r*phase2r + phase2i*phase2i)
                    phase2r = phase2r/norm;  phase2i = phase2i/norm
                    exp2theta = cmplx(phase2r,phase2i)
                    exptheta = sqrt(exp2theta)

                    pcos = real(exptheta);      psin = imag(exptheta)

                    if (iprimatdiagon >= 2) then
                       write(nfout,'(4x," pcos, psin = ",2e14.6, " ig_max = ",i5 )') &
	&                            pcos, psin,ig_max
                       write(nfout,'(4x," sqrwf1, sqrwf2, sqrwf = ",3e14.6)') &
	&                            sqrwf1, sqrwf2, sqrwf
                    end if

                    cos2theta = (wfr*wfr2 - wfi*wfi2)/sqrwf
                    cos2theta = cos2theta/norm
                    if (cos2theta > 1.0 ) cos2theta = 1.d0
                    if (cos2theta < -1.0 ) then
                       pcos = 0.d0;  psin = 1.d0
                    else
                       pcos = sqrt((1+cos2theta)/2.0); psin = sqrt((1-cos2theta)/2.0)
                    end if
                    if (pcos0 < 0.0 ) pcos = - pcos
                    if (psin0 < 0.0 ) psin = - psin
                    if (iprimatdiagon >= 2 ) then
                       write(nfout,'(4x," pcos, psin = ",2e14.6, " ig_max = ",i5 )') &
	&                                 pcos, psin,ig_max
                    end if
                 end if

                 do ig = 1, iba2(ik)
                    wfr = zaj_mat( ig+igshift,ib2,1 )
                    wfi = zaj_mat( ig+igshift,ib2,2 )
                    zaj_mat( ig+igshift,ib2,1 ) = pcos*wfr - psin*wfi
                    zaj_mat( ig+igshift,ib2,2 ) = pcos*wfi + psin*wfr
                 end do
                 zaj_mat( igshift+1,ib2,2 ) = 0.d0

              end if
            end do
#endif
            if ( iprimatdiagon >=2 ) then
               write(nfout,'(i4," Re ",5e14.6)') ib, (zaj_mat(ig,map_z(ib),1),ig=1,5)
               write(nfout,'(8x,5e14.6)') (zaj_mat(ig,map_z(ib),1),ig=6,40)
               write(nfout,'(i4," Im ",5e14.6)') ib, (zaj_mat(ig,map_z(ib),2),ig=1,5)
               write(nfout,'(8x,5e14.6)') (zaj_mat(ig,map_z(ib),2),ig=6,40)
            end if
         end do

         if (k_symmetry(ik) == GAMMA_base_symmetrization) then
            do ib = ista_e, iend_e, istep_e
               ib2 = map_z(ib)
               Do is=1, ndim_spinor
!!!!!!!                  igshift = ( is-1 ) *kg2_gamma
                  igshift = ( is-1 ) *iba2(ik)           ! ?????? kt : uncertain !!!

                  do ig = 2, kg2_gamma
                    ig1 = nbase_gamma(ig,1)
                    ig2 = nbase_gamma(ig,2)
                    wfr = zaj_mat( ig1+igshift,ib2,1 )
	            wfi = zaj_mat( ig1+igshift,ib2,2 )
                    zaj_mat( ig2+igshift,ib2,1 ) = wfr
                    zaj_mat( ig2+igshift,ib2,2 ) = -wfi
                  end do
               end do
            end do
         end if
      end if
      if(iprimatdiagon >= 2) write(nfout,'(" out of <<phase_mult_noncl>>")')
    end subroutine phase_mult_noncl

    subroutine cp_zaj_mat_to_zaj_l_noncl()
      integer :: ib,ig, i, is
      integer :: igshift

      if(k_symmetry(ik) == GAMMA) then
         if(kimg == 1) then
	    Do is=1, ndim_spinor
!!!!!!!	      igshift = kg2_gamma *(is-1)
	      igshift = iba2(ik) *(is-1)

              zaj_l( :, :, ik+is-1, 1 ) = 0.d0
              do ib = 1, np_e
                 do i = 1, kg2_gamma
                    ig = nbase(i,ik)
                    zaj_l( i, ib, ik+is-1, 1 ) = zaj_mat( ig+igshift, ib, 1 )
                 end do
               end do
            end do
         else if(kimg == 2) then
	    Do is=1, ndim_spinor
!!!!!!!	      igshift = kg2_gamma *(is-1)

	      igshift = iba2(ik) *(is-1)

              zaj_l( : ,:, ik+is-1, 1 ) = 0.d0
              zaj_l( : ,:, ik+is-1, 2 ) = 0.d0
              do ib = 1, np_e
                 do i = 1, kg2_gamma
                    ig = nbase(i,ik)
                    zaj_l( i, ib, ik+is-1, 1 ) = zaj_mat( ig+igshift, ib, 1 )
                    zaj_l( i, ib, ik+is-1, 2 ) = zaj_mat( ig+igshift, ib, 2 )
                 end do
              end do
           end do
         end if
      else
         if(kimg == 1) then
	    Do is=1, ndim_spinor
	      igshift = iba2(ik) *(is-1)
              zaj_l( :, :, ik+is-1, 1 ) = 0.d0

              do ib = 1, np_e
                 do i = 1, iba2(ik)
                    ig = nbmat2(i,ik)
                    zaj_l( ig, ib, ik+is-1, 1 ) = zaj_mat( i+igshift, ib, 1 )
                 end do
              end do
            end do
         else if(kimg == 2) then
	    Do is=1, ndim_spinor
	      igshift = iba2(ik) *(is-1)
              zaj_l( :, :, ik+is-1, 1 ) = 0.d0
              zaj_l( :, :, ik+is-1, 2 ) = 0.d0
              do ib = 1, np_e
                 do i = 1, iba2(ik)
                   ig = nbmat2(i,ik)
                   zaj_l( ig, ib, ik+is-1, 1 ) = zaj_mat( i+igshift, ib, 1 )
                   zaj_l( ig, ib, ik+is-1, 2 ) = zaj_mat( i+igshift, ib, 2 )
                 end do
               end do
            end do
         end if
      end if
#ifdef SAVE_FFT_TIMES
      if(sw_save_fft == ON) then
         do is=1,ndim_spinor
            do ib=1,np_e
               status_saved_phifftr(ib,ik+is-1) = OLD
            end do
         end do
      end if
#endif
    end subroutine cp_zaj_mat_to_zaj_l_noncl

  end subroutine m_ESmat_solve_Hx_eq_eSx_noncl
!=============================================================================== 11.0

!===============================================================================
  subroutine check_LAPACKinputs(str,nfout,n,lda,vl,vu,il,iu,lwork,abstol,RMat,CMat)
    character(len=*),intent(in) :: str
    integer, intent(in) ::         nfout,n,lda,il,iu,lwork
    real(DP), intent(in) ::        vl,vu,abstol
    real(DP),        optional, intent(in), dimension(lda,n) :: RMat
    complex(CMPLDP), optional, intent(in), dimension(lda,n) :: CMat
    integer :: i
    character(len=32) :: strp
    strp = "<<"//trim(str)//">>"
    write(nfout,'(" --- n   = ", i8,1x,a32)') n,strp
    write(nfout,'(" --- lda = ",i8)') lda
    write(nfout,'(" --- vl  = ",d20.8, " vu = ",d20.8)') vl, vu
    write(nfout,'(" --- il  = ",i8, " iu = ", i8)') il, iu
    write(nfout,'(" --- lwork = ",i8)') lwork
    write(nfout,'(" --- abstol = ", d20.8)') abstol
    if(present(RMat))  write(nfout,'(" --- amat(diagonal) = ",8f8.4)') (RMat(i,i),i=1,n)
    if(present(CMat)) then
       write(nfout,'(" --- Real(cmat(diagonal)) = ",8f8.4)') (real(CMat(i,i)),i=1,n)
       write(nfout,'(" --- Imag(cmat(diagonal)) = ",8f8.4)') (aimag(CMat(i,i)),i=1,n)
    end if
  end subroutine check_LAPACKinputs

  subroutine check_LAPACKoutputs(str,nfout,ik,n,ne,lda,evec,Rmat,Cmat)
    character(len=*),intent(in) ::       str
    integer, intent(in) ::               nfout,ik,n,ne,lda
    real(DP),intent(in), dimension(n) :: evec
    real(DP),        optional, intent(in), dimension(lda,n) :: RMat
    complex(CMPLDP), optional, intent(in), dimension(lda,n) :: CMat
    character(len=32) :: strp

    integer :: ib, ri, i
    character(len=4) :: a

    strp = "<<"//trim(str)//">>"
    write(nfout,'(" --- eigen values ik = ",i8,1x,a32," ---")')ik,strp
    write(nfout,'(8f10.6)') (evec(i),i=1,ne)

    if(Present(Rmat)) then
       do ib = ista_e, iend_e, istep_e               ! MPI
          write(nfout,'(" evec(",i4,",",i3,")= ",e14.6," ",4x,5e14.6)') &
               & ib,ik,evec(ib),(Rmat(i,ib),i=1,5)
       end do
    end if
    if(Present(Cmat)) then
       a = "     "
       do ib = ista_e, iend_e, istep_e               ! MPI
          do ri = 1, kimg
             if(ri == 1 .and. kimg == 2) a = "(Re) "
             if(ri == 2) a = "(Im) "
             if(ri == 1) write(nfout,'(" evec(",i4,",",i3,")= ",e14.6," ",a4,5e14.6)') &
                  & ib,ik,evec(ib),a,(real(Cmat(i,ib)),i=1,5)
             if(ri == 2) write(nfout,'(32x,a4,5e14.6)') a,(aimag(Cmat(i,ib)),i=1,5)
          end do
       end do
    end if

  end subroutine check_LAPACKoutputs

  subroutine check_LAPACKerror(ne_found,ne,info)
    integer, intent(in) :: ne_found,ne,info
    if (ne_found /= ne) then
       write(*,*) '### ERROR ### ne_found != ne'
       write(*,*) '   ne_found ...',ne_found
       write(*,*) '   ne       ...',ne
       call phase_error_with_msg(6, ' ne_found is not ne, stopped at <<check_LAPACKerror>>',__LINE__,__FILE__)
    end if
    if (info /= 0) then
       write(*,*) '### ERROR ### info /= 0'
       write(*,*) '   info ...',info
       call phase_error_with_msg(6, ' info is not zero, stopped at <<check_LAPACKerror>>',__LINE__,__FILE__)
    end if
  end subroutine check_LAPACKerror

  subroutine wd_CalculatedEigenvalues(nfout,ik,nmatsz,neg,ekt,str)
    integer, intent(in) :: nfout,ik,nmatsz,neg
    real(kind=DP),intent(in),dimension(nmatsz) :: ekt
    character(len=*),intent(in) :: str

    integer :: i, ib, ri
    character(len=5) :: ac

    write(nfout,'(" -- eigen values : ik = ",i8, 1x,a32," --")') ik,trim(str)
    write(nfout,'(8f10.6)') (ekt(i),i=1,neg)
    ac = "     "
    do ib = ista_e, iend_e, istep_e
       do ri = 1, kimg
          if(ri == 1 .and. kimg == 2) ac = "(Re) "
          if(ri == 2) ac = "(Im) "
          if(ri == 1) write(nfout,'(" eko(",i4,",",i3,")= ",e14.6," ",a4,5e14.6)') &
               & ib,ik,eko_l(map_z(ib),ik),ac,(zaj_mat(i,map_z(ib),ri),i=1,5)
          if(ri == 2) write(nfout,'(31x,a4,5e14.6)') ac,(zaj_mat(i,map_z(ib),ri),i=1,5)
          write(nfout,'(35x,5e14.6)') (zaj_mat(i,map_z(ib),ri),i=6,20)
       end do
    end do
  end subroutine wd_CalculatedEigenvalues

end module m_ES_WF_by_MatDiagon
