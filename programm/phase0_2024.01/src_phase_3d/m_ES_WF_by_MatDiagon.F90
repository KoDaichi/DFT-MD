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
       &                          , nrvf_ordr, nblocksize_mgs_default, fsr_l,fsi_l,neordr &
       &                          , m_ES_wd_zaj_small_portion_3D
  use m_ES_nonlocal,         only : m_ES_betar_dot_WFs_4_each_k_3D
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
       &                          , istep_e,map_z, nrank_k
  use m_Parallelization,   only : nis_G_g1k,nie_G_g1k,nel_G_g1k,ista_G_g1k,iend_G_g1k,np_G_g1k,mp_G_g1k,map_G_g1k &
       &                        , nis_B_g1k,nie_B_g1k,nel_B_g1k,ista_B_g1k,iend_B_g1k,np_B_g1k,mp_B_g1k,map_B_g1k &
       &                          , mpi_kg_world, mpi_ke_world          &
       &                          , myrank_g, myrank_e, nrank_g &
       &                          , mpi_chg_world, myrank_chg, nrank_chg &
       &                          , nis_g1k,nie_g1k,nel_g1k,ista_g1k,iend_g1k,np_g1k,mp_g1k  &
       &                          , neg_g,nrank_e,mype,mpi_k_world,ista_k,iend_k,neg_g_all,map_e &
       &                          , m_Parallel_init_mpi_G_iba2_3D, m_Parallel_init_mpi_B_iba2_3D
  use m_Control_Parameters,  only : ipriparallel,printable &
       &                          , sw_scalapack_md, method_scalapack, block_size, nprow, npcol, meg &
       &                          , nblocksize_mgs_is_given                    &
       &                          , nblocksize_mgs, divide_square
  use m_PlaneWaveBasisSet,   only : kgp, ngabc_kngp_l, ngabc_kngp_B_l, kg
  use m_Files,               only : nfout
  use m_Const_Parameters,    only : HOUSEHOLDER
  use m_Parallelization,     only : make_index_band_for_scalapack_md                &
       &                          , dealloate_index_band_for_scalapack_md           &
       &                          , nmatsz_col, nel_col_md, nis_col_md, nie_col_md, scl_col_md  &
       &                          , nmatsz_row, nel_row_md, nis_row_md, nie_row_md, scl_row_md  &
       &                          , scl_md_comm_max, scl_md_comm_max_r                 &
       &                          , scl_md_comm_rank, scl_md_comm_rank_r               &
       &                          , scl_md_comm_rno,  scl_md_comm_rno_r                &
       &                          , scl_md2_comm_max, scl_md2_comm_max_r                 &
       &                          , scl_md2_comm_rank, scl_md2_comm_rank_r               &
       &                          , scl_md2_comm_rno,  scl_md2_comm_rno_r                &
       &                          , nis_e, nie_e

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

  integer :: npcol_md, nprow_md
  integer :: nsclcol_md, nsclrow_md
!ScaLAPACK
    logical, save :: nsame


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
    allocate(hsmat(maxval(np_B_g1k),maxval(np_G_g1k),kimg)); hsmat = 0.d0
    if(modnrm == EXECUT) then
       allocate(hsmat_s(maxval(np_B_g1k),maxval(np_G_g1k),kimg)); hsmat_s = 0.d0
       allocate(linv(maxval(np_G_g1k),maxval(np_B_g1k),kimg)); linv = 0.d0
       allocate(ldag(maxval(np_G_g1k),kimg))
    end if
    allocate(zaj_mat(nmatsz,np_e,kimg))
  end subroutine alloc_matrices

  subroutine dealloc_matrices()
    if(allocated(zaj_mat)) deallocate(zaj_mat)
    if(allocated(ldag))  deallocate(ldag)
    if(allocated(linv))  deallocate(linv)
    if(allocated(zfsin)) deallocate(zfsin)
    if(allocated(zfcos)) deallocate(zfcos)
    if(allocated(vlhxc)) deallocate(vlhxc)
    if(allocated(hsmat)) deallocate(hsmat)
    if(allocated(hsmat_s)) deallocate(hsmat_s)
  end subroutine dealloc_matrices


!===============================================================================
  subroutine m_ESmat_solve_Hx_eq_eSx_3D(nfout,iteration,iteration_electronic)
    integer, intent(in) :: nfout, iteration, iteration_electronic

    integer   :: ispin, ik, iksnl, n=0
    integer   :: id_sname = -1
!---- ScaLAPACK
    integer :: is,ie, nb_mgs
    integer :: lda,lwork1, lwork2
    integer, save :: lrwork1,lrwork2,liwork1,liwork2,occ
    integer, save :: ictxt, myrow, mycol
    integer, dimension(9), save :: desca,descz
    integer, allocatable, save :: usermap(:,:)
    integer :: nmatsz_opt
!---- ScaLAPACK
    logical :: unitcell_can_change
                                                  __TIMER_SUB_START(1201)
    call tstatc0_begin('m_ESmat_solve_Hx_eq_eSx ',id_sname,1)
                                                  __TIMER_DO_START(1476)
    if(iprimatdiagon >= 2) write(nfout,'(" !! iteration = ",i6, " m_ESmat_solve_Hx_eq_eSx")') iteration

    if(ekmode == ON ) then
       if(iteration_electronic == 1) then
          call m_pwBS_set_gmaxs(n_matrix_size,gmaxs_given)  ! -> gmaxs,nmatsz,n_rGv
       else
          call m_pwBS_set_gmaxs(n,gmax)                     ! -> gmaxs,nmatsz,n_rGv
       end if
    else
       if(iteration == 1 .and. intzaj == by_matrix_diagon) then
          call m_pwBS_set_gmaxs(n_matrix_size,gmaxs_given)
       else if(skip_alloc_phonon .and. iteration_electronic == 1 .and. intzaj == by_matrix_diagon) then
          call m_pwBS_set_gmaxs(n_matrix_size,gmaxs_given)
       else if(unitcell_can_change() .and. iteration_electronic == 1 &
         &    .and. intzaj == by_matrix_diagon) then
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

!
    call m_Parallel_init_mpi_G_iba2_3D(nfout,ipriparallel,printable,kv3,iba2)
    call m_Parallel_init_mpi_B_iba2_3D(nfout,ipriparallel,printable,kv3,iba2)

    call alloc_matrices()

!fj --------------------
          if(nblocksize_mgs_is_given) then
             nb_mgs = nblocksize_mgs
          else
             nb_mgs = nblocksize_mgs_default
          end if
!fj    -----------------
                                                  __TIMER_DO_STOP(1476)
                                                  __TIMER_DO_START(1477)
    do ispin = 1, nspin, af+1
       call prjvlhxc_l2vlhxc_t()   ! vlhxc_l -> vlhxc
       do ik = ispin, kv3-nspin+ispin,nspin
          if(map_k(ik) /= myrank_k) cycle
          iksnl = (ik-1)/nspin + 1
          hsmat = 0.d0
          if(modnrm == EXECUT) then
             hsmat_s = 0.d0
             if(iprimatdiagon >=3 ) write(nfout,'(" !! modnrm == EXECUT (m_ESmat_solve_Hx_eq_eSx)")')
             call ovmtrx() ; if(iprimatdiagon >= 3) call wd_hsmat(" -- ovmtrx --")
!x          if(sw_scalapack /= ON) then
!x                call cholde() ; if(iprimatdiagon >= 3) call wd_hsmat(" -- cholde --")   ! -> linv
!x          endif
          end if
          call nlmtrx()           ! (nlmtrx) -> hsmat <- snl
                          if(iprimatdiagon >= 3) call wd_hsmat(" -- nlmtrx --")
          call amtrxkinetic()     !  -> hsmat
                          if(iprimatdiagon >= 3) call wd_hsmat(" -- amtrxkinetic --")
          call amtrxlocalpt()     !  -> hsmat
                          if(iprimatdiagon >= 3) call wd_hsmat(" -- amtrxlocalpt --")
!x          if(sw_scalapack /= ON) then
!x             if(modnrm == EXECUT) call lhl() ! -> hsmat
!x                             if(modnrm == EXECUT .and. iprimatdiagon >= 3) call wd_hsmat(" -- lhl --")
!x          endif

!---- ScaLAPACK
         if(sw_scalapack_md == ON) then
            nprow_md = nprow
            npcol_md = npcol
            if (divide_square > 0) then
               if(nprow==0 .or. npcol==0) then
                  nprow_md = int(sqrt(real(nrank_e*nrank_g)))
                  npcol_md = nprow_md
               else
                  if(nprow*npcol > nrank_g*nrank_e) then
                     nprow_md = nrank_e
                     npcol_md = nrank_g
                  end if
               end if
            else
               if(nprow==0 .or. npcol==0) then
                  nprow_md = nrank_e
                  npcol_md = nrank_g
               end if
               if(nprow*npcol > nrank_g*nrank_e) then
                  nprow_md = nrank_e
                  npcol_md = nrank_g
               end if
            end if
            if(block_size == 0) then
               block_size = nb_mgs
            end if
            if(nprow==nrank_g .and. npcol==nrank_e .and. block_size==nb_mgs) then
               nsame = .true.
            else
               nsame = .false.
            endif
            nsclrow_md = nprow_md
            nsclcol_md = npcol_md

            call set_nprow_npcol(nsclrow_md,nsclcol_md)
            allocate(usermap(nsclrow_md,nsclcol_md))

            nmatsz_opt = min(iba2(ik),nmatsz)
            call scalapack_setup_3D(nmatsz_opt,lwork1,lrwork1,liwork1,lwork2,&
               & lrwork2,liwork2,lda,occ,ictxt,myrow,mycol,desca,descz,usermap)
            call make_index_band_for_scalapack_md(neg, nmatsz_opt, ik, nb_mgs, block_size &
                 &                              , nsclrow_md, nsclcol_md)
               ! -> scl_md_comm_max, scl_md_comm_max_r, scl_md_comm_rank, scl_md_comm_rank_r
               ! -> scl_md_comm_rno, scl_md_comm_rno_r, scl_md2_comm_max, scl_md2_comm_max_r
               ! -> scl_md2_comm_rank, scl_md2_comm_rank_r, scl_md2_comm_rno,  scl_md2_comm_rno_r
            deallocate(usermap)
         end if

!---- ScaLAPACK
         if(sw_scalapack_md == ON) then
              call solve_Hx_eq_ex_ScaLAPACK()
         else
#ifdef _USE_LAPACK_
      ! [LAPACK] with LAPACK library
      !***** LAPACK [begin] ***** M.Okamoto
          call solve_Hx_eq_ex_LAPACK()
      !***** LAPACK [end] ***** M.Okamoto
#else
      ! [LAPACK] without LAPACK library
          call solve_Hx_eq_ex()
#endif
         endif
          if(iprimatdiagon >= 3) call wd_hsmat(" -- after diagonalization --")
!x          if(sw_scalapack /= ON) then
!x             if(modnrm == EXECUT) call lx() ! -> zaj_mat
!x          endif
          ! Debug
          if(k_symmetry(ik) == GAMMA .or. k_symmetry(ik) == GAMMA_base_symmetrization) call phase_mult()
          call cp_zaj_mat_to_zaj_l() ! zaj_mat -> zaj_l
          if(sw_scalapack_md == ON) then
            call dealloate_index_band_for_scalapack_md
          endif
       end do
    end do
                                                  __TIMER_DO_STOP(1477)
                                                  __TIMER_DO_START(1478)
    if(iprimatdiagon >= 3) then
       do ik = 1, kv3
          if(map_k(ik) /= myrank_k) cycle
          call m_ES_wd_zaj_small_portion_3D(nfout,ik," (m_ESmat_solve_Hx_eq_eSx)",26)
       end do
    end if

    call dealloc_matrices()
    call m_pwBS_dealloc_igpo()
    call m_pwBS_dealloc_nbmat_and_iba2()

!f    call m_ES_betar_dot_WFs(nfout) ! (fsrfsi) -> fsr_l, fsi_l
       do ik = 1, kv3, af+1
          if(map_k(ik) /= myrank_k) cycle         ! MPI
!          call m_ES_betar_dot_WFs_3D(nfout,ik) ! (fsrfsi) -> fsr_l, fsi_l
          call m_ES_betar_dot_WFs_4_each_k_3D(nfout,ik)
       end do
                                                  __TIMER_DO_STOP(1478)
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(1201)
  contains
    ! Contained subroutines
    !  1)    subroutine wd_hsmat
    !  2)    subroutine solve_Hx_eq_ex_LAPACK
    !  3)*   subroutine solve_Hx_eq_ex         * is same with the corresponding subroutine contained in
    !  4)*   subroutine calphase                        a subroutine m_ESmat_solve_Hx_eq_eSx
    !  5)    subroutine nlmtrx
    !  6)    subroutine ovmtrx
    !  7)    subroutine cholde
    !  8)    subroutine amtrxkinetic
    !  9)$   subroutine prjvlhxc_l2vlhxc_t     $ is almost same
    ! 10)    subroutine amtrxlocalpt
    ! 11)    subroutine lhl
    ! 12)    subroutine lx
    ! 13)    subroutine phase_mult
    ! 14)    subroutine cp_zaj_mat_to_zaj_l
    ! 15)    subroutine solve_Hx_eq_ex_ScaLAPACK
    ! 16)    subroutine set_nprow_npcol
    ! 17)    subroutine scalapack_setup_3D
    ! 18)    subroutine trans_scalapack
    ! 19)    subroutine trans_scalapack_r
    ! 20)    subroutine eigsend

    subroutine wd_hsmat(str)
      character(len=*),intent(in) :: str
      integer :: i, imax, j, jmax
      integer, parameter :: NMATSIZE = 10
                                                  __TIMER_SUB_START(1202)
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
                                                  __TIMER_SUB_STOP(1202)
    end subroutine wd_hsmat

#ifdef _USE_LAPACK_
   ! [LAPACK] with LAPACK library
   !***** LAPACK [begin] ***** M.Okamoto
    subroutine solve_Hx_eq_ex_LAPACK()
      integer :: n, lda, ne, ne_found, il, iu, info, i, j, ib, lwork, ib_wk, index_G, index_B, isize, ierr
      real(kind=DP) :: abstol, vl, vu, sum, dlamch
      integer,allocatable ::              ifail(:), iwork(:)
      real(kind=DP),allocatable ::        amat(:,:), evec(:), rvmat(:,:), rwork(:)
      complex(kind=CMPLDP),allocatable :: cmat(:,:), cmat2(:,:), cvmat(:,:), cwork(:)
      complex(kind=CMPLDP),parameter :: IMG = (0.d0,1.d0)
      integer ::  mpi_comm
#ifdef _FJ_DBG0_
      integer ::                          iiii,jjjj
#endif
!! === DEBUG by tkato 2011/09/06 ================================================
!      complex(kind=CMPLDP) :: ctemp
!! ==============================================================================
                                                  __TIMER_SUB_START(1203)
      mpi_comm = mpi_k_world(myrank_k)
      zaj_mat(:,:,:) = 0.d0
      lda = nmatsz
      !!$ vl = -1.d10 ; vu = 1.d10 ; il = ista_e ; iu = iend_e
      vl = -1.d10 ; vu = 1.d10 ; il = 1 ; iu = neg
      ne = iu - il + 1
      abstol = 2*dlamch('S')
      n = iba2(ik)
      isize = lda*n
      if (kimg == 1) then
         lwork = max(1,12*n)
        !+++++++++++++++++++++++++++++++++++++++++++++++++
         allocate(evec(n),iwork(5*n),ifail(n))
! === DEBUG by tkato 2012/11/06 ================================================
!        allocate(amat(lda,n),rvmat(lda,ne),rwork(lwork))
         allocate(amat(lda,n),rvmat(lda,n),rwork(lwork))
! ==============================================================================
        !+++++++++++++++++++++++++++++++++++++++++++++++++
         amat(:,:) = 0.d0
                                                  __TIMER_DO_START(1286)
         do i = ista_B_g1k(ik), iend_B_g1k(ik)
         do j = ista_G_g1k(ik), iend_G_g1k(ik)
            index_B = i-ista_B_g1k(ik)+1
            index_G = j-ista_G_g1k(ik)+1
            amat(i,j) = hsmat(index_B,index_G,1)
         end do
         end do
                                                  __TIMER_DO_STOP(1286)
                                                  __TIMER_COMM_START_w_BARRIER(mpi_comm,1263)
         call mpi_allreduce(MPI_IN_PLACE,amat,isize,mpi_double_precision,mpi_sum, mpi_comm,ierr)
                                                  __TIMER_COMM_STOP(1263)

         if(iprimatdiagon>= 2)  call check_LAPACKinputs("solve_Hx_eq_ex_LAPACK "  &
              &                                   ,nfout,n,lda,vl,vu,il,iu,lwork,abstol,RMat=amat)
         if(modnrm == EXECUT) then
            rvmat(:,:) = 0.d0
            do i = ista_B_g1k(ik), iend_B_g1k(ik)
               do j = ista_G_g1k(ik), iend_G_g1k(ik)
                  index_B = i-ista_B_g1k(ik)+1
                  index_G = j-ista_G_g1k(ik)+1
                  rvmat(i,j) = hsmat_s(index_B,index_G,1)
               end do
            end do
            call mpi_allreduce(MPI_IN_PLACE,rvmat,isize,mpi_double_precision,mpi_sum,mpi_comm,ierr)
                                                  __TIMER_DO_START(1287)
            Call dsygv(1,'V','U',n,amat,lda,rvmat,lda,evec,rwork,lwork,info)
            ne_found = ne   !!<<== DEBUG by tkato 2012/11/06 =====
                                                  __TIMER_DO_STOP(1287)
         else
                                                  __TIMER_DO_START(1287)
            call dsyevx('V','I','U',n,amat,lda,vl,vu,il,iu,abstol, &
                 & ne_found,evec,rvmat,lda,rwork,lwork,iwork,ifail,info)
                                                  __TIMER_DO_STOP(1287)
         endif
         if(iprimatdiagon>=2) call check_LAPACKoutputs("solve_Hx_eq_ex_LAPACK_noncl2 " &
              &                                   ,nfout,ik,n,ne,lda,evec,Rmat=rvmat)
         call check_LAPACKerror(ne_found,ne,info)
                                                  __TIMER_DO_START(1288)
         do ib = ista_e, iend_e, istep_e
            j = ib - ista_e +1                 !  j = map_z(ib)
            ib_wk = neg_g(j)
            do i = 1,n
               zaj_mat(i,j,1) = rvmat(i,ib_wk) ! zaj_l(i,j,ik,1)=rvmat(i,ib-ista_e+1)! M.Okamoto (August 22,2003)
            end do
            eko_l(j,ik) = evec(ib_wk)          !  eko_l(j,ik) = evec(ib-ista_e+1) ! M.Okamoto (August 22,2003)
         end do
                                                  __TIMER_DO_STOP(1288)
         if(af/=0) eko_l(:,ik+1) = eko_l(:,ik)
        !+++++++++++++++++++++++++++++
         deallocate(amat,rvmat,rwork)
        !+++++++++++++++++++++++++++++
      else if (kimg == 2) then
         lwork = max(1,6*n)
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         allocate(evec(n),iwork(5*n),ifail(n))
! === DEBUG by tkato 2011/09/06 ================================================
!        allocate(cmat(lda,n),cvmat(lda,ne),cwork(lwork),rwork(7*n))
         allocate(cmat(lda,n),cmat2(lda,n),cvmat(lda,ne),cwork(lwork),rwork(7*n))
         evec = 0.d0; iwork = 0; ifail = 0; cmat = dcmplx(0.d0,0.d0)
         cvmat = dcmplx(0.d0,0.d0); rwork = 0.d0;
! ==============================================================================
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                  __TIMER_DO_START(1289)
         do i = ista_B_g1k(ik), iend_B_g1k(ik)          !f       do i = 1,n
         do j = ista_G_g1k(ik), iend_G_g1k(ik)          !f       do j = 1,n
            index_B = i-ista_B_g1k(ik)+1
            index_G = j-ista_G_g1k(ik)+1
            cmat(i,j) = dcmplx(hsmat(index_B,index_G,1),hsmat(index_B,index_G,2))
!f            cmat(i,j) = cmplx(hsmat(i,j,1),hsmat(i,j,2))
         end do
         end do
                                                  __TIMER_DO_STOP(1289)
                                                  __TIMER_COMM_START_w_BARRIER(mpi_comm,1264)
         call mpi_allreduce(MPI_IN_PLACE,cmat,isize,mpi_double_complex,mpi_sum, mpi_comm,ierr)
                                                  __TIMER_COMM_STOP(1264)

! === DEBUG by tkato 2011/09/06 ================================================
!        call zheevx('V','I','U',n,cmat,lda,vl,vu,il,iu,abstol, &
!           ne_found,evec,cvmat,lda,ctemp,-1,rwork,iwork,ifail,info)
!        lwork = int(real(ctemp))
!        allocate(cwork(lwork)); cwork = dcmplx(0.d0,0.d0)
! ==============================================================================
        if(iprimatdiagon>= 2)  call check_LAPACKinputs("solve_Hx_eq_ex_LAPACK "  &
              &                                   ,nfout,n,lda,vl,vu,il,iu,lwork,abstol,CMat=cmat)
        if(modnrm == EXECUT) then
         cmat2 = dcmplx(0.d0,0.d0)
         do i = ista_B_g1k(ik), iend_B_g1k(ik)
         do j = ista_G_g1k(ik), iend_G_g1k(ik)
            index_B = i-ista_B_g1k(ik)+1
            index_G = j-ista_G_g1k(ik)+1
            cmat2(i,j) = dcmplx(hsmat_s(index_B,index_G,1),hsmat_s(index_B,index_G,2))
         end do
         end do
         call mpi_allreduce(MPI_IN_PLACE,cmat2,isize,mpi_double_complex,mpi_sum,mpi_comm,ierr)
                                                  __TIMER_DO_START(1287)
         call zhegv(1,'V','U',n,cmat,lda,cmat2,lda,evec,cwork,lwork,rwork,info)
         cvmat(:,1:ne) = cmat(:,1:ne)
                                                  __TIMER_DO_STOP(1287)
        else
                                                  __TIMER_DO_START(1287)
           call zheevx('V','I','U',n,cmat,lda,vl,vu,il,iu,abstol, &
                & ne_found,evec,cvmat,lda,cwork,lwork,rwork,iwork,ifail,info)
                                                  __TIMER_DO_STOP(1287)
        endif
         if(iprimatdiagon>=2) call check_LAPACKoutputs("solve_Hx_eq_ex_LAPACK_noncl2 " &
              &                                   ,nfout,ik,n,ne,lda,evec,Cmat=cvmat)
         if (info /= 0) then
            write(*,*) '### ERROR ### info /= 0'
            write(*,*) '   info ...',info
            call phase_error_with_msg(nfout,'### ERROR ### info /= 0',__LINE__,__FILE__)
         end if
                                                  __TIMER_DO_START(1290)
         do ib = ista_e, iend_e, istep_e
            j = ib - ista_e +1                        !f         j = map_z(ib)
            ib_wk = neg_g(j)
            do i = 1,n
               zaj_mat(i,j,1) = real(cvmat(i,ib_wk))
               zaj_mat(i,j,2) = aimag(cvmat(i,ib_wk))
               !!$ zaj_l(i,j,ik,1)    = real(cvmat(i,ib-ista_e+1))  ! M.Okamoto (August 22,2003)
               !!$ zaj_l(i,j,ik,kimg) = aimag(cvmat(i,ib-ista_e+1)) ! M.Okamoto (August 22,2003)
            end do
            eko_l(j,ik) = evec(ib_wk)  !!$ eko_l(j,ik) = evec(ib-ista_e+1) ! M.Okamoto (August 22,2003)
         end do
                                                  __TIMER_DO_STOP(1290)
           if(af/=0) eko_l(:,ik+1) = eko_l(:,ik)
        !+++++++++++++++++++++++++++++++++++
         deallocate(cmat,cvmat,cwork,rwork)
        !+++++++++++++++++++++++++++++++++++
      end if
#ifdef _FJ_DBG0_
       if(mype .eq. 0) then
            write(1002,'("evec")')
            write(1002,'(5f10.7)') (evec(iiii),iiii=1,neg)

            write(1002,'("zaj_mat:kimg=1")')
            do iiii = 1, nb_mgs
               write(1002,'("iiii=",i4)') iiii
               write(1002,'(8(f10.7))') (zaj_mat(jjjj,iiii,1),jjjj=1,200)
            enddo

            write(1002,'("zaj_mat:kimg=2")')
            do iiii = 1, nb_mgs
               write(1002,'("iiii=",i4)') iiii
               write(1002,'(8(f10.7))') (zaj_mat(jjjj,iiii,2),jjjj=1,200)
            enddo
        endif
#endif

      if (iprimatdiagon >= 2) call wd_CalculatedEigenvalues(nfout,ik,n,ne,evec,"solve_Hx_eq_ex")

      sum = 0.d0
      do i = 1,ne
         sum = sum + abs(evec(i))
      end do
      if (sum < SmallestPositiveNumber*1.d5)  call phase_error_with_msg(nfout,' ! illegal evec (solve_Hx_eq_ex)'&
                                                                       ,__LINE__,__FILE__)
      if (iprimatdiagon >= 2) write(nfout,'(" -- sum = ",f20.8)') sum
     !+++++++++++++++++++++++++++++
      deallocate(evec,iwork,ifail)
     !+++++++++++++++++++++++++++++
                                                  __TIMER_SUB_STOP(1203)
    end subroutine solve_Hx_eq_ex_LAPACK
   !***** LAPACK [end] ***** M.Okamoto
#else
   ! [LAPACK] without LAPACK library
    subroutine solve_Hx_eq_ex()
      real(kind=DP),allocatable,dimension(:,:) :: a
      real(kind=DP),allocatable,dimension(:)   :: ekt
      real(kind=DP),allocatable,dimension(:,:) :: v
      real(kind=DP),allocatable,dimension(:,:) :: w
      complex(kind=CMPLDP),allocatable,dimension(:,:) :: zw

      real(kind=DP) :: eps
      integer :: i,j,ierr,ib
      allocate(ekt(nmatsz))
      allocate(v(nmatsz*kimg,neg))

      zaj_mat(:,:,:) = 0.d0
      eps = eps_solve_Hx_eq_ex    !      default value of (eps_solve_Hx_eq_ex) = 1.d-15
      if(kimg == 1) then
         allocate(w(iba2(ik),6))
         call hobsvw(hsmat,nmatsz,iba2(ik),ekt,-neg,v,neg,eps,w,ierr)
         do ib = ista_e, iend_e, istep_e
            j = map_z(ib)
            do i = 1, iba2(ik)
               zaj_mat(i,j,1) = v(i,ib)
            end do
            eko_l(j,ik) = ekt(ib)
         end do
         deallocate(w)
! --------------
         do i = 1, nmatsz
            do j = 1, nmatsz
               hsmat(i,j,1) = v(i,j)
            end do
         end do
! -------------
      else if(kimg == 2) then
         allocate(w(3,iba2(ik)))
         allocate(zw(5,iba2(ik)))
         allocate(a(nmatsz*kimg,nmatsz)); a = 0.d0
         do j=1,nmatsz
            do i=1,nmatsz
               a(i*2-1,j) = hsmat(i,j,1)
               a(i*2,  j) = hsmat(i,j,2)
            end do
         end do
         call chobsd(a,nmatsz,iba2(ik),ekt,-neg,v,neg,eps,w,zw,ierr)
         do ib = ista_e, iend_e, istep_e
            j = map_z(ib)
            do i = 1, iba2(ik)
               zaj_mat(i,j,1) = v(i*2-1,ib)
               zaj_mat(i,j,2) = v(i*2, ib)
            end do
            eko_l(j,ik) = ekt(ib)
         end do
         deallocate(a,zw,w)
! --------------
         do i = 1, nmatsz
            do j = 1, nmatsz
               hsmat(i,j,1) = v(i*2-1,j)
               hsmat(i,j,2) = v(i*2,j)
            end do
         end do
! -------------
      end if

      if(iprimatdiagon >= 2) call wd_CalculatedEigenvalues(nfout,ik,nmatsz,neg,ekt,"solve_Hx_eq_ex")
!!$!!$         call m_ES_wd_zaj_small_portion(nfout,ik," (solve_Hx_eq_ex)",17)

      eps = 0.d0
      do ib = 1, neg
         eps = eps + abs(ekt(ib))
      end do
      if(eps < SmallestPositiveNumber*1.d5) call phase_error_with_msg(nfout,' ! illegal ekt (solve_Hx_eq_ex)',__LINE__,__FILE__)
      if (iprimatdiagon >= 2) write(nfout,'(" -- sum = ",f20.8)') eps

      deallocate(v); deallocate(ekt)
    end subroutine solve_Hx_eq_ex
#endif
   ! [LAPACK]

    subroutine calphase(ia)
      integer, intent(in) :: ia
      integer       :: i
      real(kind=DP) :: ph
                                                  __TIMER_SUB_START(1204)
                                                  __TIMER_DO_START(1291)
      do i = 1, nmatsz2
         ph = (pos(ia,1)*ngabc(i,1)+pos(ia,2)*ngabc(i,2)+pos(ia,3)*ngabc(i,3))*PAI2
         zfcos(i) = cos(ph)
         zfsin(i) = sin(ph)
      end do
                                                  __TIMER_DO_STOP(1291)
      zfcos(0) = 0.d0; zfsin(0) = 0.d0
                                                  __TIMER_SUB_STOP(1204)
    end subroutine calphase

    subroutine nlmtrx

      real(kind=DP),pointer,dimension(:,:,:) :: a
      real(kind=DP),allocatable,dimension(:,:)  :: s
      integer, pointer, dimension(:,:)       :: ng

      complex(kind=CMPLDP),parameter :: zi = (0.d0,1.d0)
      integer :: i,j,it,ia,u,v,n,lmt1,il1,lmt2,il2,j1,l1,l2,l3,ip,y,i1, i2,j2,m
      real(kind=DP) :: fr,f21,f12,t
      integer :: ista,  index_B, index_G
                                                  __TIMER_SUB_START(1205)
      a => hsmat
      ng => ngabc
!!$      s => snl(:,:,iksnl)
      allocate(s(kg1_ext,nlmtt)); s = 0.d0

      if(k_symmetry(ik) == GAMMA) then
         if(iprimatdiagon >= 2) write(nfout,'(" k_symmetry(ik) == GAMMA <<nlmtrx>>")')
         if(iprimatdiagon >= 2) write(nfout,'(" kg2_gamma, iba(ik), iba2(ik) = ",3i8)') kg2_gamma,iba(ik),iba2(ik)
                                                  __TIMER_DO_START(1292)
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
!f            s(1,j) = snl(1,j,iksnl)
            if( myrank_g == 0) s(1,j) = snl(1,j,iksnl)
            if(mod(il1,2) == 0) then
               do i = ista_g1k(ik), iend_g1k(ik)              !f    do i = 2, kg2_gamma
                  if ( i == 1) cycle
                  if ( i > kg2_gamma) cycle
                  i1 = nbase(i,ik)
                  i2 = nbase_gamma(i,2)
                  s(i1,j) =  snl(i-ista_g1k(ik)+1,j,iksnl)    !f      s(i1,j) = snl(i,j,iksnl)
                  s(i2,j) = -snl(i-ista_g1k(ik)+1,j,iksnl)    !f      s(i2,j) = -snl(i,j,iksnl)
               end do
            else
               do i = ista_g1k(ik), iend_g1k(ik)              !f    do i = 2, kg2_gamma
                  if ( i == 1) cycle
                  if ( i > kg2_gamma) cycle
                  i1 = nbase(i,ik)
                  i2 = nbase_gamma(i,2)
                  s(i1,j) = snl(i-ista_g1k(ik)+1,j,iksnl)     !f      s(i1,j) = snl(i,j,iksnl)
                  s(i2,j) = snl(i-ista_g1k(ik)+1,j,iksnl)     !f      s(i2,j) = snl(i,j,iksnl)
               end do
            end if
         end do
                                                  __TIMER_DO_STOP(1292)
      else
         if(iprimatdiagon >= 2) write(nfout,'(" k_symmetry(ik) /= GAMMA <<nlmtrx>>")')
                                                  __TIMER_DO_START(1293)
         do j = 1, nlmtt
            do i = ista_g1k(ik), iend_g1k(ik)                 !f      do i = 1, iba(ik)
               s(i,j) = snl(i-ista_g1k(ik)+1,j,iksnl)         !f        s(i,j) = snl(i,j,iksnl)
            end do
         end do
                                                  __TIMER_DO_STOP(1293)
      end if
                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,1265)
      call mpi_allreduce(MPI_IN_PLACE,s,kg1_ext*nlmtt, mpi_double_precision, mpi_sum, mpi_ke_world, ierr)
                                                  __TIMER_COMM_STOP(1265)
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
                                                  __TIMER_DO_START(1294)
                  do j = ista_G_g1k(ik), iend_G_g1k(ik)           !f        do j = 1, iba2(ik)
                     j1 = nbmat(j,ik)
                     j2 = nbmat2(j,ik)
                     l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)

                     ista = max(j,ista_B_g1k(ik))

                     if(kimg == 1) then
                        do i = ista, iend_B_g1k(ik)               !f          do i = j, iba2(ik)
                           n = nbmat(i,ik)
                           m = nbmat2(i,ik)
                           y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                           index_B = i-ista_B_g1k(ik)+1
                           index_G = j-ista_G_g1k(ik)+1
                           a(index_B,index_G,1) = a(index_B,index_G,1)+f21*(s(m,u)*s(j2,v)+s(m,v)*s(j2,u))*zfcos(y)
!f                           a(i,j,1) = a(i,j,1)+f21*(s(m,u)*s(j2,v)+s(m,v)*s(j2,u))*zfcos(y)
                        end do
                     else
                        do i = ista, iend_B_g1k(ik)               !f          do i = j, iba2(ik)
                           n = nbmat(i,ik)
                           m = nbmat2(i,ik)
                           y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                           t  = f21*(s(m,u)*s(j2,v)+s(m,v)*s(j2,u))
                           index_B = i-ista_B_g1k(ik)+1
                           index_G = j-ista_G_g1k(ik)+1
                           a(index_B,index_G,1) = a(index_B,index_G,1) + t*zfcos(y) !f a(i,j,1)=a(i,j,1)+t*zfcos(y)
                           a(index_B,index_G,2) = a(index_B,index_G,2) - t*zfsin(y) !f a(i,j,2)=a(i,j,2)-t*zfsin(y)
                        end do
                     end if
                  end do
                                                  __TIMER_DO_STOP(1294)
               else
                  f21 = real(zi**(il2-il1))*fr*iwei(ia)
                  f12 = real(zi**(il1-il2))*fr*iwei(ia)
!!$                  write(nfout,'(" zi**(il1-il2) = ",2f16.8)') zi**(il1-il2)
                                                  __TIMER_DO_START(1295)
                  do j = ista_G_g1k(ik), iend_G_g1k(ik)    !f        do j = 1, iba2(ik)
                     j1 = nbmat(j,ik)
                     j2 = nbmat2(j,ik)
                     l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)

                     ista = max(j,ista_B_g1k(ik))
                     if(kimg == 1) then
                        do i = ista, iend_B_g1k(ik)        !f          do i = j, iba2(ik)
                           n = nbmat(i,ik)
                           m = nbmat2(i,ik)
                           y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                           index_B = i-ista_B_g1k(ik)+1
                           index_G = j-ista_G_g1k(ik)+1
                         a(index_B,index_G,1)=a(index_B,index_G,1)+(f21*s(m,u)*s(j2,v)+f12*s(m,v)*s(j2,u))*zfcos(y)
!f                           a(i,j,1) = a(i,j,1)+(f21*s(m,u)*s(j2,v)+f12*s(m,v)*s(j2,u))*zfcos(y)
                        end do
                     else
                        do i = ista, iend_B_g1k(ik)        !f          do i = j, iba2(ik)
                           n = nbmat(i,ik)
                           m = nbmat2(i,ik)
                           y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                           t  = f21*s(m,u)*s(j2,v)+f12*s(m,v)*s(j2,u)
                           index_B = i-ista_B_g1k(ik)+1
                           index_G = j-ista_G_g1k(ik)+1
                           a(index_B,index_G,1) = a(index_B,index_G,1) + t*zfcos(y) !f a(i,j,1)=a(i,j,1)+t*zfcos(y)
                           a(index_B,index_G,2) = a(index_B,index_G,2) - t*zfsin(y) !f a(i,j,2)=a(i,j,2)-t*zfsin(y)
                        end do
                     end if
                  end do
                                                  __TIMER_DO_STOP(1295)
               end if
            else if(mod(il1+il2,2) == 1) then
               if(ipaw(it)==0) then
                   fr = vlhxcQ(lmt1,lmt2,ia,ispin)
               else
                   fr = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
               end if
               f21 = -aimag(zi**(il2-il1))*fr*iwei(ia)
               f12 = -aimag(zi**(il1-il2))*fr*iwei(ia)
                                                  __TIMER_DO_START(1296)
               do j = ista_G_g1k(ik), iend_G_g1k(ik)       !f          do j = 1, iba2(ik)
                  j1 = nbmat(j,ik)
                  j2 = nbmat2(j,ik)
                  l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)

                  ista = max(j,ista_B_g1k(ik))
                  if(kimg == 1) then
                     do i = ista, iend_B_g1k(ik)           !f            do i = j, iba2(ik)
                        n = nbmat(i,ik)
                        m = nbmat2(i,ik)
                        y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                        index_B = i-ista_B_g1k(ik)+1
                        index_G = j-ista_G_g1k(ik)+1
!f                        a(i,j,1) = a(i,j,1)+ &
!f                            & (f21*s(m,u)*s(j2,v)+f12*s(m,v)*s(j2,u))*zfsin(y)
                        a(index_B,index_G,1) = a(index_B,index_G,1)+ &
                            & (f21*s(m,u)*s(j2,v)+f12*s(m,v)*s(j2,u))*zfsin(y)
                     end do
                  else
                     do i = ista, iend_B_g1k(ik)           !f            do i = j, iba2(ik)
                        n = nbmat(i,ik)
                        m = nbmat2(i,ik)
                        y = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
                        t = f21*s(m,u)*s(j2,v)+f12*s(m,v)*s(j2,u)
                        index_B = i-ista_B_g1k(ik)+1
                        index_G = j-ista_G_g1k(ik)+1
                        a(index_B,index_G,1) = a(index_B,index_G,1) + t*zfsin(y) !f  a(i,j,1)=a(i,j,1)+t*zfsin(y)
                        a(index_B,index_G,2) = a(index_B,index_G,2) + t*zfcos(y) !f  a(i,j,2)=a(i,j,2)+t*zfcos(y)
                     end do
                  end if
               end do
                                                  __TIMER_DO_STOP(1296)
            end if
         end do
         do lmt1 = 1, ilmt(it)
            if(ipaw(it)==0) then
                fr = (dion(lmt1,lmt1,it) + vlhxcQ(lmt1,lmt1,ia,ispin))*iwei(ia)
            else
                fr = (dion_paw(lmt1,lmt1,ispin,ia) + vlhxcQ(lmt1,lmt1,ia,ispin))*iwei(ia)
            endif
            u  = lmtt(lmt1,it)
                                                  __TIMER_DO_START(1297)
            do j = ista_G_g1k(ik), iend_G_g1k(ik)       !f            do j = 1, iba2(ik)
               j1 = nbmat(j,ik)
               j2 = nbmat2(j,ik)
               l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)

               ista = max(j,ista_B_g1k(ik))
               if(kimg == 1) then
                  do i = ista, iend_B_g1k(ik)          !f          do i = j, iba2(ik)
                     i1 = nbmat(i,ik)
                     i2 = nbmat2(i,ik)
                     y = igpo(ng(i1,1)-l1,ng(i1,2)-l2,ng(i1,3)-l3)
                     index_B = i-ista_B_g1k(ik)+1
                     index_G = j-ista_G_g1k(ik)+1
                     a(index_B,index_G,1) = a(index_B,index_G,1) + fr*s(i2,u)*s(j2,u)*zfcos(y)
!f                     a(i,j,1) = a(i,j,1) + fr*s(i2,u)*s(j2,u)*zfcos(y)
                  end do
               else if(kimg == 2) then
                  do i = ista, iend_B_g1k(ik)          !f          do i = j, iba2(ik)
                     i1 = nbmat(i,ik)
                     i2 = nbmat2(i,ik)
                     y = igpo(ng(i1,1)-l1,ng(i1,2)-l2,ng(i1,3)-l3)
                     t = fr*s(i2,u)*s(j2,u)
                     index_B = i-ista_B_g1k(ik)+1
                     index_G = j-ista_G_g1k(ik)+1
                     a(index_B,index_G,1) = a(index_B,index_G,1) + t*zfcos(y)  !f a(i,j,1) = a(i,j,1)+t*zfcos(y)
                     a(index_B,index_G,2) = a(index_B,index_G,2) - t*zfsin(y)  !f a(i,j,2) = a(i,j,2)-t*zfsin(y)
                  end do
               end if
            end do
                                                  __TIMER_DO_STOP(1297)
         end do
      end do
      deallocate(s)
                                                  __TIMER_SUB_STOP(1205)
    end subroutine nlmtrx

    subroutine ovmtrx

      integer :: i,j,it,ia,u,v,n,ip,y,jn,j1,l1,l2,l3,i1,lmt1,lmt2
      integer :: i2,j2,m,jm
      real(kind=DP) :: fr,t
      real(kind=DP),allocatable,dimension(:,:)   :: s

      real(kind=DP),allocatable,dimension(:,:,:)   :: a_wk
      integer :: ista, index_B, index_G, k, mpi_comm
                                                  __TIMER_SUB_START(1206)
      mpi_comm = mpi_k_world(myrank_k)
      allocate(s(kg1_ext,nlmtt)); s = 0.0d0
      allocate(a_wk(nmatsz,nmatsz,kimg)); a_wk = 0.d0

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

            if( myrank_g == 0) s(1,j) = snl(1,j,iksnl)
            ip = 1
            if(mod(l1,2) == 0)  ip = -1
                                                  __TIMER_DO_START(1298)
            do i = ista_g1k(ik), iend_g1k(ik)            !f        do i = 2, iba(ik)
               if ( i == 1) cycle
               i1 = nbase(i,ik)
               i2 = nbase_gamma(i,2)
               s(i1,j) =    snl(i-ista_g1k(ik)+1,j,iksnl)  !f           s(i1,j) =             snl(i,j,iksnl)
               s(i2,j) = ip*snl(i-ista_g1k(ik)+1,j,iksnl)  !f           s(i2,j) = (-1)^{l1+1}*snl(i,j,iksnl)
            end do
                                                  __TIMER_DO_STOP(1298)
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
                                                  __TIMER_DO_START(1300)
            do i = ista_g1k(ik), iend_g1k(ik)               !f         do i = 1, kg1
               s(i,j) = snl(i-ista_g1k(ik)+1,j,iksnl)
            end do
                                                  __TIMER_DO_STOP(1300)
         end do
      end if
                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,1266)
      call mpi_allreduce(MPI_IN_PLACE, s, kg1_ext*nlmtt, mpi_double_precision, mpi_sum, mpi_ke_world, ierr)
                                                  __TIMER_COMM_STOP(1266)
      hsmat = 0.d0
!!$      do j = ista_G_g1k(ik), iend_G_g1k(ik)
!!$      do i = ista_B_g1k(ik), iend_B_g1k(ik)
!!$         index_B = i-ista_B_g1k(ik)+1
!!$         index_G = j-ista_G_g1k(ik)+1
!!$         if( i == j ) hsmat(index_B,index_G,1) = 1.d0       !f     a(i,i,1) = 1.d0
!!$      end do; end do
      do i = ista_G_g1k(ik), iend_G_g1k(ik)       !f      do i = 1, iba2(ik)
         if(ista_B_g1k(ik) <= i .and. i <= iend_B_g1k(ik)) then
            index_B = i-ista_B_g1k(ik)+1
            index_G = i-ista_G_g1k(ik)+1
            hsmat(index_B,index_G,1) = 1.d0       !f     a(i,i,1) = 1.d0
         end if
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
                                                  __TIMER_DO_START(1301)
            do j = ista_G_g1k(ik), iend_G_g1k(ik)      !f            do j = 1, iba2(ik)
               jn = nbmat(j,ik)
               jm = nbmat2(j,ik)
               l1 = ngabc(jn,1); l2 = ngabc(jn,2); l3 = ngabc(jn,3)
               ista = max(j,ista_B_g1k(ik))
               if(kimg == 1) then
                  do i = ista, iend_B_g1k(ik)          !f              do i = j, iba2(ik)
                     n = nbmat(i,ik)
                     m = nbmat2(i,ik)
                     y = igpo(ngabc(n,1)-l1,ngabc(n,2)-l2,ngabc(n,3)-l3)
                     index_B = i-ista_B_g1k(ik)+1
                     index_G = j-ista_G_g1k(ik)+1
                     hsmat(index_B,index_G,1) = hsmat(index_B,index_G,1)+fr*(s(m,u)*s(jm,v)+s(m,v)*s(jm,u))*zfcos(y)
!f                     a(i,j,1) = a(i,j,1)+fr*(s(m,u)*s(jm,v)+s(m,v)*s(jm,u))*zfcos(y)
                  end do
               else
                  do i = ista, iend_B_g1k(ik)            !f                  do i = j, iba2(ik)
                     n = nbmat(i,ik)
                     m = nbmat2(i,ik)
                     y = igpo(ngabc(n,1)-l1,ngabc(n,2)-l2,ngabc(n,3)-l3)
                     t  = fr*(s(m,u)*s(jm,v)+s(m,v)*s(jm,u))
                     index_B = i-ista_B_g1k(ik)+1
                     index_G = j-ista_G_g1k(ik)+1
                     hsmat(index_B,index_G,1) = hsmat(index_B,index_G,1)+t*zfcos(y) !f a(i,j,1)=a(i,j,1)+t*zfcos(y)
                     hsmat(index_B,index_G,2) = hsmat(index_B,index_G,2)-t*zfsin(y) !f a(i,j,2)=a(i,j,2)-t*zfsin(y)
                  end do
               end if
            end do
                                                  __TIMER_DO_STOP(1301)
         end do
         do lmt1 = 1, ilmt(it)
            u  = lmtt(lmt1,it)
            fr = q(lmt1,lmt1,it)*iwei(ia)
                                                  __TIMER_DO_START(1302)
            do j = ista_G_g1k(ik), iend_G_g1k(ik)       !f        do j = 1, iba2(ik)
               j1 = nbmat(j,ik)
               j2 = nbmat2(j,ik)
               l1 = ngabc(j1,1); l2 = ngabc(j1,2); l3 = ngabc(j1,3)
               ista = max(j,ista_B_g1k(ik))
               if(kimg == 1) then
                  do i = ista, iend_B_g1k(ik)           !f          do i = j, iba2(ik)
                     i1 = nbmat(i,ik)
                     i2 = nbmat2(i,ik)
                     ip = igpo(ngabc(i1,1)-l1,ngabc(i1,2)-l2,ngabc(i1,3)-l3)
                     index_B = i-ista_B_g1k(ik)+1
                     index_G = j-ista_G_g1k(ik)+1
                     hsmat(index_B,index_G,1) = hsmat(index_B,index_G,1)+fr*s(i2,u)*s(j2,u)*zfcos(ip)
!f                     a(i,j,1) = a(i,j,1)+fr*s(i2,u)*s(j2,u)*zfcos(ip)
                  end do
               else if(kimg == 2) then
                  do i = ista, iend_B_g1k(ik)           !f          do i = j, iba2(ik)
                     i1 = nbmat(i,ik)
                     i2 = nbmat2(i,ik)
                     ip = igpo(ngabc(i1,1)-l1,ngabc(i1,2)-l2,ngabc(i1,3)-l3)
                     t  = fr*s(i2,u)*s(j2,u)
                     index_B = i-ista_B_g1k(ik)+1
                     index_G = j-ista_G_g1k(ik)+1
                     hsmat(index_B,index_G,1)=hsmat(index_B,index_G,1)+t*zfcos(ip) !f a(i,j,1)=a(i,j,1)+t*zfcos(ip)
                     hsmat(index_B,index_G,2)=hsmat(index_B,index_G,2)-t*zfsin(ip) !f a(i,j,2)=a(i,j,2)-t*zfsin(ip)
                  end do
               end if
            end do
                                                  __TIMER_DO_STOP(1302)
         end do
      end do

      do k = 1, kimg
      do j = ista_G_g1k(ik), iend_G_g1k(ik)
      do i = ista_B_g1k(ik), iend_B_g1k(ik)
         index_G = j-ista_G_g1k(ik)+1
         index_B = i-ista_B_g1k(ik)+1
         a_wk(i,j,k) = hsmat(index_B,index_G,k)
      end do; end do; end do
                                                  __TIMER_COMM_START_w_BARRIER(mpi_comm,1267)
      call mpi_allreduce(MPI_IN_PLACE,a_wk,nmatsz*nmatsz*kimg,mpi_double_precision,mpi_sum,mpi_comm,ierr)
                                                  __TIMER_COMM_STOP(1267)

                                                  __TIMER_DO_START(1303)
      do k = 1, kimg
      do j = 2, iba2(ik)
      do i = 1, j-1
         a_wk(i,j,k) = a_wk(j,i,k)
      end do;  end do;  end do
                                                  __TIMER_DO_STOP(1303)
                                                  __TIMER_DO_START(1304)
      do k = 1, kimg
      do j = ista_G_g1k(ik), iend_G_g1k(ik)
      do i = ista_B_g1k(ik), iend_B_g1k(ik)
         index_G = j-ista_G_g1k(ik)+1
         index_B = i-ista_B_g1k(ik)+1
         hsmat(index_B,index_G,k) = a_wk(i,j,k)
         if(modnrm == EXECUT) hsmat_s(index_B,index_G,k) = a_wk(i,j,k)
      end do; end do; end do
                                                  __TIMER_DO_STOP(1304)
      deallocate(s)
      deallocate(a_wk)
                                                  __TIMER_SUB_STOP(1206)
    end subroutine ovmtrx

    subroutine cholde

! Choleski decomposition.
! ref: S.G.Louie,K.M.Ho,and M.L.Cohen, PRB19,(1979),p1774.
      integer       :: n,i,j
      real(kind=DP) :: t
      integer  :: imax, jmax
      integer, parameter :: NMATSIZE = 10

      integer  :: iend, jend, index_G, index_B, isize, nend, index_n, mpi_comm
      real(kind=DP),allocatable,dimension(:,:) :: ldag_wk
      real(kind=DP),allocatable,dimension(:) :: linv_wk1
      real(kind=DP),allocatable,dimension(:,:) :: linv_wk2, hsmat_wk
                                                  __TIMER_SUB_START(1207)
      mpi_comm = mpi_k_world(myrank_k)
      allocate(ldag_wk(np_G_g1k(ik),kimg))
      allocate(linv_wk1(kimg))
      allocate(linv_wk2(nmatsz,kimg))
      allocate(hsmat_wk(nmatsz,kimg))

      linv = 0.d0
!f      linv(1,1,1) = 1.d0/sqrt(hsmat(1,1,1))
      if((myrank_g == 0) .and. (myrank_e ==0)) linv(1,1,1) = 1.d0/sqrt(hsmat(1,1,1))
      if(kimg == 1) then
         do n = 1, iba2(ik)-1
!---
          index_n = n-ista_G_g1k(ik)+1
          isize = np_B_g1k(ik)*kimg
          hsmat_wk = 0.0d0
                                                  __TIMER_DO_START(1305)
          if((nis_G_g1k(myrank_g,ik) <= n) .and.      &
   &             (nie_G_g1k(myrank_g,ik) > n)) then
             do i = ista_B_g1k(ik), iend_B_g1k(ik)
                index_B = i-ista_B_g1k(ik)+1
                hsmat_wk(i,1)=hsmat(index_B,index_n+1,1)
             end do
          endif
                                                  __TIMER_DO_STOP(1305)
                                                  __TIMER_DO_START(1306)
          if(nis_G_g1k(myrank_g,ik) == n+1) then
             do i = ista_B_g1k(ik), iend_B_g1k(ik)
                index_B = i-ista_B_g1k(ik)+1
                hsmat_wk(i,1)=hsmat(index_B,1,1)
             end do
          endif
                                                  __TIMER_DO_STOP(1306)
          isize = nmatsz*kimg
                                                  __TIMER_COMM_START_w_BARRIER(mpi_comm,1268)
          call mpi_allreduce(MPI_IN_PLACE,hsmat_wk,isize,mpi_double_precision,mpi_sum, mpi_comm,ierr)
                                                  __TIMER_COMM_STOP(1268)
!---
            ldag = 0.d0
            ldag_wk = 0.d0
                                                  __TIMER_DO_START(1307)
            do i = ista_G_g1k(ik), min(n,iend_G_g1k(ik))      !f            do  i = 1, n
               jend = min(i, iend_B_g1k(ik))
               do j = ista_B_g1k(ik), jend                    !f               do j = 1, i
                  index_G = i-ista_G_g1k(ik)+1
                  index_B = j-ista_B_g1k(ik)+1
                  ldag_wk(index_G,1) = ldag_wk(index_G,1) + linv(index_G,index_B,1)*hsmat_wk(j,1)
               end do
            end do
                                                  __TIMER_DO_STOP(1307)
        isize = np_G_g1k(ik)*kimg
                                                  __TIMER_COMM_START_w_BARRIER(mpi_kg_world,1269)
         call mpi_allreduce(ldag_wk, ldag, isize, mpi_double_precision, mpi_sum, mpi_kg_world, ierr)
                                                  __TIMER_COMM_STOP(1269)

                                                  __TIMER_DO_START(1308)
            t = 0.d0
            do i = ista_G_g1k(ik),  min(n,iend_G_g1k(ik))     !f            do i = 1, n
               index_G = i-ista_G_g1k(ik)+1
               t = t + ldag(index_G,1)**2                     !f               t = t + ldag(i,1)**2
            end do
                                                  __TIMER_DO_STOP(1308)

                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,1270)
            call mpi_allreduce(MPI_IN_PLACE, t, 1, mpi_double_precision, mpi_sum, mpi_ke_world, ierr)
                                                  __TIMER_COMM_STOP(1270)

          if(nie_G_g1k(myrank_g,ik) == n+1) then
             ldag(index_n+1,1) = sqrt(hsmat_wk(n+1,1)-t)
          endif

                                                  __TIMER_DO_START(1309)
            linv_wk2 = 0.0d0
            linv_wk1(1) = 1.d0/sqrt(hsmat_wk(n+1,1)-t)
            do i = ista_B_g1k(ik), min(n,iend_B_g1k(ik))        !f            do i = 1, n
               do j = ista_G_g1k(ik), min(n,iend_G_g1k(ik))     !f               do j = 1, n
                  index_G = j-ista_G_g1k(ik)+1
                  index_B = i-ista_B_g1k(ik)+1
                  linv_wk2(i,1)=linv_wk2(i,1)-linv_wk1(1)*ldag(index_G,1)*linv(index_G,index_B,1)
               end do
            end do
                                                  __TIMER_DO_STOP(1309)

          isize = nmatsz*kimg
                                                  __TIMER_COMM_START_w_BARRIER(mpi_comm,1271)
         call mpi_allreduce(MPI_IN_PLACE, linv_wk2, isize, mpi_double_precision, mpi_sum, mpi_comm,ierr)
                                                  __TIMER_COMM_STOP(1271)

          linv_wk2(n+1,1) = linv_wk1(1)
    nend = min(n+1,iend_B_g1k(ik))
    if(nrank_g > 1) then
                                                  __TIMER_DO_START(1310)
          if((nis_G_g1k(myrank_g,ik) < n+1) .and.      &
   &             (nie_G_g1k(myrank_g,ik) >= n+1)) then
            do i = ista_B_g1k(ik), nend
              index_B = i-ista_B_g1k(ik)+1
              index_G = n-ista_G_g1k(ik)+1
              linv(index_G+1,index_B,1)=linv_wk2(i,1)
            end do
          else if (nis_G_g1k(myrank_g,ik) == n+1) then
            do i = ista_B_g1k(ik), nend
              index_B = i-ista_B_g1k(ik)+1
              index_G = n-ista_G_g1k(ik)+1
              linv(1,index_B,1)=linv_wk2(i,1)
            end do
          endif
                                                  __TIMER_DO_STOP(1310)
    else
                                                  __TIMER_DO_START(1311)
            do i = ista_B_g1k(ik), nend
              index_B = i-ista_B_g1k(ik)+1
              index_G = n-ista_G_g1k(ik)+1
                linv(index_G+1,index_B,1)=linv_wk2(i,1)
            end do
                                                  __TIMER_DO_STOP(1311)
    endif

         end do
      else if(kimg == 2) then
         do n = 1, iba2(ik)-1
!---
          index_n = n-ista_G_g1k(ik)+1
          isize = np_B_g1k(ik)*kimg
          hsmat_wk = 0.0d0
                                                  __TIMER_DO_START(1312)
          if((nis_G_g1k(myrank_g,ik) <= n) .and.      &
   &             (nie_G_g1k(myrank_g,ik) > n)) then
             do i = ista_B_g1k(ik), iend_B_g1k(ik)
                index_B = i-ista_B_g1k(ik)+1
                hsmat_wk(i,1)=hsmat(index_B,index_n+1,1)
                hsmat_wk(i,2)=hsmat(index_B,index_n+1,2)
             end do
          endif
                                                  __TIMER_DO_STOP(1312)
                                                  __TIMER_DO_START(1313)
          if(nis_G_g1k(myrank_g,ik) == n+1) then
             do i = ista_B_g1k(ik), iend_B_g1k(ik)
                index_B = i-ista_B_g1k(ik)+1
                hsmat_wk(i,1)=hsmat(index_B,1,1)
                hsmat_wk(i,2)=hsmat(index_B,1,2)
             end do
          endif
                                                  __TIMER_DO_STOP(1313)

          isize = nmatsz*kimg
                                                  __TIMER_COMM_START_w_BARRIER(mpi_comm,1272)
          call mpi_allreduce(MPI_IN_PLACE,hsmat_wk,isize,mpi_double_precision,mpi_sum, mpi_comm,ierr)
                                                  __TIMER_COMM_STOP(1272)
!---
            ldag = 0.d0
            ldag_wk = 0.d0
            iend = min(n,iend_G_g1k(ik))
                                                  __TIMER_DO_START(1314)
            do i = ista_G_g1k(ik), iend                      !f            do i = 1, n
               jend = min(i,iend_B_g1k(ik))
               do j = ista_B_g1k(ik), jend                   !f               do j = 1, i
                  index_G = i-ista_G_g1k(ik)+1
                  index_B = j-ista_B_g1k(ik)+1
                  ldag_wk(index_G,1) = ldag_wk(index_G,1) &
                       & + linv(index_G,index_B,1)*hsmat_wk(j,1)-linv(index_G,index_B,2)*hsmat_wk(j,2)
                  ldag_wk(index_G,2) = ldag_wk(index_G,2) &
                       & + linv(index_G,index_B,1)*hsmat_wk(j,2)+linv(index_G,index_B,2)*hsmat_wk(j,1)
               end do
            end do
                                                  __TIMER_DO_STOP(1314)
         isize = np_G_g1k(ik)*kimg
                                                  __TIMER_COMM_START_w_BARRIER(mpi_kg_world,1273)
         call mpi_allreduce(ldag_wk, ldag, isize, mpi_double_precision, mpi_sum, mpi_kg_world, ierr)
                                                  __TIMER_COMM_STOP(1273)

            t = 0.d0
                                                  __TIMER_DO_START(1315)
            do i = ista_G_g1k(ik), min(n,iend_G_g1k(ik))      !f            do i = 1, n
               index_G = i-ista_G_g1k(ik)+1
               t = t + ldag(index_G,1)**2 + ldag(index_G,2)**2
!f               t = t + ldag(i,1)**2 + ldag(i,2)**2
            end do
                                                  __TIMER_DO_STOP(1315)
                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,1274)
            call mpi_allreduce(MPI_IN_PLACE, t, 1, mpi_double_precision, mpi_sum, mpi_ke_world, ierr)
                                                  __TIMER_COMM_STOP(1274)

!f            ldag(n+1,1) = sqrt(hsmat(n+1,n+1,1)-t)
!f            linv(n+1,n+1,1) = 1.d0/sqrt(hsmat(n+1,n+1,1)-t)
          if(nie_G_g1k(myrank_g,ik) == n+1) then
             ldag(index_n+1,1) = sqrt(hsmat_wk(n+1,1)-t)
          endif
!--          if((nis_G_g1k(myrank_g,ik) < n+1) .and.      &
!--   &             (nie_G_g1k(myrank_g,ik) >= n+1)) then
!--             ldag(index_n+1,1) = sqrt(hsmat_wk(n+1,1)-t)
!--          else if (nis_G_g1k(myrank_g,ik) == n+1) then
!--             ldag(1,1) = sqrt(hsmat_wk(n+1,1)-t)
!--          endif

                                                  __TIMER_DO_START(1316)
            linv_wk2 = 0.0d0
            linv_wk1(1) = 1.d0/sqrt(hsmat_wk(n+1,1)-t)
            linv_wk1(2) = 0.0d0
            iend = min(n,iend_B_g1k(ik))
            do i = ista_B_g1k(ik), iend                   !f            do i = 1, n
               jend = iend_G_g1k(ik) ; if( n < iend_G_g1k(ik)) jend = n
               do j = ista_G_g1k(ik), jend                !f               do j = 1, n
                  index_G = j-ista_G_g1k(ik)+1
                  index_B = i-ista_B_g1k(ik)+1

                  linv_wk2(i,1) = linv_wk2(i,1) &
                       &       - linv_wk1(1)*(ldag(index_G,1)*linv(index_G,index_B,1) &
                       &                     +ldag(index_G,2)*linv(index_G,index_B,2))&
                       &       + linv_wk1(2)*(ldag(index_G,1)*linv(index_G,index_B,2)&
                       &                     -ldag(index_G,2)*linv(index_G,index_B,1) )
                  linv_wk2(i,2) = linv_wk2(i,2) &
                       &       - linv_wk1(1)*(ldag(index_G,1)*linv(index_G,index_B,2) &
                       &                     -ldag(index_G,2)*linv(index_G,index_B,1))&
                       &       - linv_wk1(2)*(ldag(index_G,1)*linv(index_G,index_B,1)&
                       &                     -ldag(index_G,2)*linv(index_G,index_B,2) )
               end do
            end do
                                                  __TIMER_DO_STOP(1316)
          isize = nmatsz*kimg
                                                  __TIMER_COMM_START_w_BARRIER(mpi_comm,1275)
         call mpi_allreduce(MPI_IN_PLACE, linv_wk2, isize, mpi_double_precision, mpi_sum, mpi_comm,ierr)
                                                  __TIMER_COMM_STOP(1275)

          linv_wk2(n+1,1) = linv_wk1(1)
    nend = min(n+1,iend_B_g1k(ik))
    if(nrank_g > 1) then

                                                  __TIMER_DO_START(1317)
          if((nis_G_g1k(myrank_g,ik) < n+1) .and.      &
   &             (nie_G_g1k(myrank_g,ik) >= n+1)) then
            do i = ista_B_g1k(ik), nend
              index_B = i-ista_B_g1k(ik)+1
              index_G = n-ista_G_g1k(ik)+1
              linv(index_G+1,index_B,1)=linv_wk2(i,1)
              linv(index_G+1,index_B,2)=linv_wk2(i,2)
            end do
          else if (nis_G_g1k(myrank_g,ik) == n+1) then
            do i = ista_B_g1k(ik), nend
              index_B = i-ista_B_g1k(ik)+1
              index_G = n-ista_G_g1k(ik)+1
              linv(1,index_B,1)=linv_wk2(i,1)
              linv(1,index_B,2)=linv_wk2(i,2)
            end do
          endif
                                                  __TIMER_DO_STOP(1317)
    else
                                                  __TIMER_DO_START(1318)
            do i = ista_B_g1k(ik), nend
              index_B = i-ista_B_g1k(ik)+1
              index_G = n-ista_G_g1k(ik)+1
                linv(index_G+1,index_B,1)=linv_wk2(i,1)
                linv(index_G+1,index_B,2)=linv_wk2(i,2)
            end do
                                                  __TIMER_DO_STOP(1318)
    endif

         end do
      end if
      if(iprimatdiagon >= 2) then
         imax = NMATSIZE; if(imax > nmatsz) imax = nmatsz
         jmax = NMATSIZE; if(jmax > nmatsz) jmax = nmatsz
         do j = 1, jmax
            write(nfout,'(" linv(1:*,",i3,",1): ",10f8.4,99(/18x,10f8.4))') j, (linv(i,j,1),i=1,imax)
         end do
         if(kimg == 2) then
            do j = 1, jmax
               write(nfout,'(" linv(1:*,",i3,",2): ",10f8.4,99(/18x,10f8.4))') j, (linv(i,j,2),i=1,imax)
            end do
         end if
      end if
      deallocate(ldag_wk)
      deallocate(linv_wk1,linv_wk2)
      deallocate(hsmat_wk)
                                                  __TIMER_SUB_STOP(1207)
    end subroutine cholde

    subroutine amtrxkinetic
      integer :: i,j
      real(kind=DP) :: ga,gb,gc
      integer :: index_B,index_G
                                                  __TIMER_SUB_START(1208)
                                                  __TIMER_DO_START(1319)
      do i = ista_G_g1k(ik), iend_G_g1k(ik)     !f      do i = 1, iba2(ik)
         if(ista_B_g1k(ik) <= i .and. i <= iend_B_g1k(ik)) then
            index_B = i-ista_B_g1k(ik)+1
            index_G = i-ista_G_g1k(ik)+1
            j = nbmat(i,ik)
            ga = vkxyz(ik,1,BUCS) + ngabc(j,1)
            gb = vkxyz(ik,2,BUCS) + ngabc(j,2)
            gc = vkxyz(ik,3,BUCS) + ngabc(j,3)
            hsmat(index_B,index_G,1) = hsmat(index_B,index_G,1) + 0.5 & !f  hsmat(i,i,1) = hsmat(i,i,1) + 0.5 &
                 & * (ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
                 & +  ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga)
         end if
      end do
!!$      do jj = ista_G_g1k(ik), iend_G_g1k(ik)
!!$      do i  = ista_B_g1k(ik), iend_B_g1k(ik)
!!$       if( i == jj) then
!!$         index_B =  i-ista_B_g1k(ik)+1
!!$         index_G = jj-ista_G_g1k(ik)+1
!!$         j = nbmat(i,ik)
!!$         ga = vkxyz(ik,1,BUCS) + ngabc(j,1)
!!$         gb = vkxyz(ik,2,BUCS) + ngabc(j,2)
!!$         gc = vkxyz(ik,3,BUCS) + ngabc(j,3)
!!$         hsmat(index_B,index_G,1) = hsmat(index_B,index_G,1) + 0.5 &
!!$              & * (ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
!!$              & +  ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga)
!!$       endif
!!$      end do; end do
                                                  __TIMER_DO_STOP(1319)
                                                  __TIMER_SUB_STOP(1208)
    end subroutine amtrxkinetic

    subroutine prjvlhxc_l2vlhxc_t()
      integer :: i,irc
                                                  __TIMER_SUB_START(1209)
      vlhxc = 0.d0
                                                  __TIMER_DO_START(1320)
      do irc = 1, kimg
         do i = ista_kngp, min(iend_kngp,nmatsz2)
            vlhxc(i,irc) = vlhxc_l(i,irc,ispin)
         end do
      end do
                                                  __TIMER_DO_STOP(1320)
      if(npes >= 2) then
                                                  __TIMER_COMM_START_w_BARRIER(mpi_chg_world,1276)
         call mpi_allreduce(MPI_IN_PLACE,vlhxc,nmatsz2*kimg,mpi_double_precision &
!f              &            ,mpi_sum,MPI_CommGroup,ierr)
              &            ,mpi_sum,mpi_chg_world,ierr)
                                                  __TIMER_COMM_STOP(1276)
      end if
                                                  __TIMER_SUB_STOP(1209)
    end subroutine prjvlhxc_l2vlhxc_t

    subroutine amtrxlocalpt()
      integer :: j,j1,l1,l2,l3,i,n,ip
      real(kind=DP),pointer,dimension(:,:,:) :: a
      integer,      pointer,dimension(:,:)   :: ng

      real(kind=DP),allocatable,dimension(:,:,:)   :: a_wk
      integer :: ista, index_B, index_G, k, mpi_comm

      a  => hsmat
      ng => ngabc
                                                  __TIMER_SUB_START(1210)

      mpi_comm = mpi_k_world(myrank_k)
      allocate(a_wk(nmatsz,nmatsz,kimg)); a_wk = 0.d0
                                                  __TIMER_DO_START(1321)
      do j = ista_G_g1k(ik), iend_G_g1k(ik)                              !f do j = 1, iba2(ik)
         j1 = nbmat(j,ik)
         l1 = ng(j1,1); l2 = ng(j1,2); l3 = ng(j1,3)
         ista = max(j,ista_B_g1k(ik))
         if(kimg == 1) then
            do i = ista, iend_B_g1k(ik)                                  !f  do i = j, iba2(ik)
               n = nbmat(i,ik)
               ip = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
               index_B = i-ista_B_g1k(ik)+1
               index_G = j-ista_G_g1k(ik)+1
               a(index_B,index_G,1) = a(index_B,index_G,1) + vlhxc(ip,1) !f    a(i,j,1) = a(i,j,1) + vlhxc(ip,1)
            end do
         else
            do i = ista, iend_B_g1k(ik)                                  !f  do i = j, iba2(ik)
               n = nbmat(i,ik)
               ip = igpo(ng(n,1)-l1,ng(n,2)-l2,ng(n,3)-l3)
               index_B = i-ista_B_g1k(ik)+1
               index_G = j-ista_G_g1k(ik)+1
               a(index_B,index_G,1) = a(index_B,index_G,1) + vlhxc(ip,1) !f    a(i,j,1) = a(i,j,1) + vlhxc(ip,1)
               a(index_B,index_G,2) = a(index_B,index_G,2) + vlhxc(ip,2) !f    a(i,j,2) = a(i,j,2) + vlhxc(ip,2)
            end do
         end if
      end do
                                                  __TIMER_DO_STOP(1321)

                                                  __TIMER_DO_START(1322)
      do k = 1, kimg
      do j = ista_G_g1k(ik), iend_G_g1k(ik)
      do i = ista_B_g1k(ik), iend_B_g1k(ik)
         index_G = j-ista_G_g1k(ik)+1
         index_B = i-ista_B_g1k(ik)+1
         a_wk(i,j,k) = a(index_B,index_G,k)
      end do; end do; end do
                                                  __TIMER_DO_STOP(1322)
                                                  __TIMER_COMM_START_w_BARRIER(mpi_comm,1277)
      call mpi_allreduce(MPI_IN_PLACE,a_wk,nmatsz*nmatsz*kimg, mpi_double_precision, mpi_sum, mpi_comm,ierr)
                                                  __TIMER_COMM_STOP(1277)

                                                  __TIMER_DO_START(1323)
      do k = 1, kimg
         l1 = 3 -2*k
         do j = 2, iba2(ik)
         do i = 1, j-1
            a_wk(i,j,k) = l1*a_wk(j,i,k)
      end do; end do; end do
                                                  __TIMER_DO_STOP(1323)

                                                  __TIMER_DO_START(1324)
      do k = 1, kimg
      do j = ista_G_g1k(ik), iend_G_g1k(ik)
      do i = ista_B_g1k(ik), iend_B_g1k(ik)
         index_G = j-ista_G_g1k(ik)+1
         index_B = i-ista_B_g1k(ik)+1
         a(index_B,index_G,k) = a_wk(i,j,k)
      end do; end do; end do
                                                  __TIMER_DO_STOP(1324)
      deallocate(a_wk)
                                                  __TIMER_SUB_STOP(1210)
    end subroutine amtrxlocalpt

    subroutine lhl

      integer :: i,j,k, index_B, index_G, n, mpi_comm
      real(kind=DP),allocatable,dimension(:,:,:) :: ww, h_wk1, linv_wk

                                                  __TIMER_SUB_START(1211)
      mpi_comm = mpi_k_world(myrank_k)
      allocate(ww(nmatsz,nmatsz,kimg)); ww = 0.d0
      allocate(h_wk1(nmatsz,np_G_g1k(ik),kimg)); h_wk1 = 0.0d0
      allocate(linv_wk(nmatsz,np_B_g1k(ik),kimg)); linv_wk = 0.0d0

                                                  __TIMER_DO_START(1325)
      do k = 1, kimg
      do j = ista_B_g1k(ik), iend_B_g1k(ik)
      do i = ista_G_g1k(ik), iend_G_g1k(ik)
         index_B = j-ista_B_g1k(ik)+1
         index_G = i-ista_G_g1k(ik)+1
         linv_wk(i,index_B,k) = linv(index_G,index_B,k)
      end do; end do; end do
                                                  __TIMER_DO_STOP(1325)
      n = nmatsz*np_B_g1k(ik)*kimg
                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_world,1278)
      call mpi_allreduce(MPI_IN_PLACE, linv_wk, n, mpi_double_precision, mpi_sum, mpi_ke_world, ierr)
                                                  __TIMER_COMM_STOP(1278)

      if(kimg == 1) then
                                                  __TIMER_DO_START(1326)
         do i = 1, iba2(ik)
            do j = ista_G_g1k(ik), iend_G_g1k(ik)                !f            do j = 1, iba2(ik)
               do k = ista_B_g1k(ik), min(i,iend_B_g1k(ik))      !f               do k = 1, i
                  index_B = k-ista_B_g1k(ik)+1
                  index_G = j-ista_G_g1k(ik)+1
                  ww(i,j,1) = ww(i,j,1) + linv_wk(i,index_B,1)*hsmat(index_B,index_G,1)
!f                                                           ww(i,j,1) = ww(i,j,1) + linv(i,k,1)*h(k,j,1)
               end do
            end do
         end do
                                                  __TIMER_DO_STOP(1326)

         n = nmatsz*nmatsz*kimg
                                                  __TIMER_COMM_START_w_BARRIER(mpi_comm,1279)
         call mpi_allreduce(MPI_IN_PLACE, ww, n, mpi_double_precision, mpi_sum, mpi_comm, ierr)
                                                  __TIMER_COMM_STOP(1279)

         h_wk1 = 0.0d0                                            !f        h = 0.d0
                                                  __TIMER_DO_START(1327)
         do i = 1, iba2(ik)
           do j = ista_G_g1k(ik), iend_G_g1k(ik)                  !f        do j = 1, iba2(ik)
               do k = ista_B_g1k(ik), min(j,iend_B_g1k(ik) )      !f          do k = 1, j
                  index_B = k-ista_B_g1k(ik)+1
                  index_G = j-ista_G_g1k(ik)+1
                  h_wk1(i,index_G,1) = h_wk1(i,index_G,1) + ww(i,k,1)*linv(index_G,index_B,1)
!f                  h(i,j,1) = h(i,j,1) + ww(i,k,1)*linv(j,k,1)
               end do
            end do
         end do
                                                  __TIMER_DO_STOP(1327)
         n = nmatsz*np_G_g1k(ik)*kimg
                                                  __TIMER_COMM_START_w_BARRIER(mpi_kg_world,1280)
         call mpi_allreduce(MPI_IN_PLACE, h_wk1, n, mpi_double_precision, mpi_sum, mpi_kg_world, ierr)
                                                  __TIMER_COMM_STOP(1280)

      else if(kimg == 2) then
                                                  __TIMER_DO_START(1328)
      do j = ista_B_g1k(ik), iend_B_g1k(ik)
      do i = ista_G_g1k(ik), iend_G_g1k(ik)
         index_B = j-ista_B_g1k(ik)+1
         index_G = i-ista_G_g1k(ik)+1
         h_wk1(j,index_G,1) = hsmat(index_B,index_G,1)
      end do; end do
                                                  __TIMER_DO_STOP(1328)
      n = nmatsz*np_G_g1k(ik)
                                                  __TIMER_COMM_START_w_BARRIER(mpi_kg_world,1281)
      call mpi_allreduce(MPI_IN_PLACE, h_wk1, n, mpi_double_precision, mpi_sum, mpi_kg_world, ierr)
                                                  __TIMER_COMM_STOP(1281)

                                                  __TIMER_DO_START(1329)
         do i = 1, iba2(ik)
            do j = ista_G_g1k(ik), iend_G_g1k(ik)                 !f        do j = 1, iba2(ik)
               do k = ista_B_g1k(ik), min(i,iend_B_g1k(ik))       !f           do k = 1, i
                  index_B = k-ista_B_g1k(ik)+1
                  index_G = j-ista_G_g1k(ik)+1
                  ww(i,j,1) = ww(i,j,1) + linv_wk(i,index_B,1)*hsmat(index_B,index_G,1)  &
     &                                  - linv_wk(i,index_B,2)*hsmat(index_B,index_G,2)
                  ww(i,j,2) = ww(i,j,2) + linv_wk(i,index_B,1)*hsmat(index_B,index_G,2)  &
!     &                                 + linv_wk(i,index_B,2)*hsmat(index_B,index_G,1)
     &                                  + linv_wk(i,index_B,2)*h_wk1(i,index_G,1)
!f                  ww(i,j,1) = ww(i,j,1)+linv(i,k,1)*h(k,j,1)-linv(i,k,2)*h(k,j,2)
!f                  ww(i,j,2) = ww(i,j,2)+linv(i,k,1)*h(k,j,2)+linv(i,k,2)*h(i,j,1)
               end do
            end do
         end do
                                                  __TIMER_DO_STOP(1329)
         n = nmatsz*nmatsz*kimg
                                                  __TIMER_COMM_START_w_BARRIER(mpi_comm,1282)
         call mpi_allreduce(MPI_IN_PLACE, ww, n, mpi_double_precision, mpi_sum, mpi_comm, ierr)
                                                  __TIMER_COMM_STOP(1282)
         h_wk1 = 0.0d0
                                                  __TIMER_DO_START(1330)
         do i = 1, iba2(ik)
            do j = ista_G_g1k(ik), iend_G_g1k(ik)                 !f        do j = 1, iba2(ik)
               do k = ista_B_g1k(ik), min(j,iend_B_g1k(ik))       !f           do k = 1, j
                  index_B = k-ista_B_g1k(ik)+1
                  index_G = j-ista_G_g1k(ik)+1
                  h_wk1(i,index_G,1) = h_wk1(i,index_G,1) + ww(i,k,1)*linv(index_G,index_B,1) &
       &                                                  + ww(i,k,2)*linv(index_G,index_B,2)
                  h_wk1(i,index_G,2) = h_wk1(i,index_G,2) + ww(i,k,2)*linv(index_G,index_B,1) &
       &                                                  - ww(i,k,1)*linv(index_G,index_B,2)
!f                  h(i,j,1) = h(i,j,1)+ww(i,k,1)*linv(j,k,1)+ww(i,k,2)*linv(j,k,2)
!f                  h(i,j,2) = h(i,j,2)+ww(i,k,2)*linv(j,k,1)-ww(i,k,1)*linv(j,k,2)
               end do
            end do
         end do
                                                  __TIMER_DO_STOP(1330)

         n = nmatsz*np_G_g1k(ik)*kimg
                                                  __TIMER_COMM_START_w_BARRIER(mpi_kg_world,1283)
         call mpi_allreduce(MPI_IN_PLACE, h_wk1, n, mpi_double_precision, mpi_sum, mpi_kg_world, ierr)
                                                  __TIMER_COMM_STOP(1283)
      end if

      hsmat = 0.0d0
                                                  __TIMER_DO_START(1331)
      do k = 1, kimg
      do j = ista_G_g1k(ik), iend_G_g1k(ik)
      do i = ista_B_g1k(ik), iend_B_g1k(ik)
         index_G = j-ista_G_g1k(ik)+1
         index_B = i-ista_B_g1k(ik)+1
         hsmat(index_B,index_G,k) = h_wk1(i,index_G,k)
      end do; end do; end do
                                                  __TIMER_DO_STOP(1331)

      deallocate(ww)
      deallocate(h_wk1)
      deallocate(linv_wk)
                                                  __TIMER_SUB_STOP(1211)
    end subroutine lhl

    subroutine lx()

#ifdef HIUX
! *** 'poption''s have been inserted by Dr. Kino (National Institute  ***
! *** for Materials Science, Japan).  12 July 2005                    ***
#endif
      integer :: ib, ib2, i, j, k, isize
      real(kind=DP),allocatable,dimension(:,:,:)  :: linv_wk
      integer :: index_G, index_B
                                                  __TIMER_SUB_START(1212)
      allocate(linv_wk(np_G_g1k(ik),nmatsz,kimg));linv_wk = 0.0d0
!
                                                  __TIMER_DO_START(1332)
         do k = 1, kimg
            do i = ista_G_g1k(ik), iend_G_g1k(ik)
               do j = ista_B_g1k(ik), iend_B_g1k(ik)
                  index_G = i-ista_G_g1k(ik)+1
                  index_B = j-ista_B_g1k(ik)+1
                  linv_wk(index_G,j,k) = linv(index_G,index_B,k)
               end do
            end do
         end do
                                                  __TIMER_DO_STOP(1332)
         isize = np_G_g1k(ik)*nmatsz*kimg
                                                  __TIMER_COMM_START_w_BARRIER(mpi_kg_world,1264)
         call mpi_allreduce(MPI_IN_PLACE,linv_wk,isize,mpi_double_precision,mpi_sum, mpi_kg_world,ierr)
                                                  __TIMER_COMM_STOP(1264)

      if(kimg == 1) then
         do ib = ista_e, iend_e, istep_e
            ib2 = ib - ista_e +1                       !f        ib2 = map_z(ib)
                                                  __TIMER_DO_START(1333)
            do i = ista_G_g1k(ik), iend_G_g1k(ik)      !f        do i = 1, iba2(ik)
               index_G = i-ista_G_g1k(ik)+1
               ldag(index_G,1) = zaj_mat(i,ib2,1)      !f           ldag(i,1) = zaj_mat(i,ib2,1)
            end do
                                                  __TIMER_DO_STOP(1333)
            zaj_mat(:,ib2,1) = 0.d0
                                                  __TIMER_DO_START(1334)
            do i = 1, iba2(ik)
!            do i = ista_B_g1k(ik), iend_B_g1k(ik)
               do j = ista_G_g1k(ik), iend_G_g1k(ik)   !f        do j = 1, iba2(ik)
                  index_B = i-ista_B_g1k(ik)+1
                  index_G = j-ista_G_g1k(ik)+1
                  zaj_mat(i,ib2,1) = zaj_mat(i,ib2,1) + linv_wk(index_G,i,1)*ldag(index_G,1)
               end do
            end do
                                                  __TIMER_DO_STOP(1334)
         end do

      else if(kimg == 2) then
         do ib = ista_e, iend_e, istep_e
            ib2 = ib - ista_e +1                       !f           ib2 = map_z(ib)
            if(iba2(ik) > nmatsz)  call phase_error_with_msg(nfout,' iba2(ik) > nmatsz',__LINE__,__FILE__)
                                                  __TIMER_DO_START(1335)

            do i = ista_G_g1k(ik), iend_G_g1k(ik)      !f        do i = 1, iba2(ik)
               index_G = i-ista_G_g1k(ik)+1
               ldag(index_G,1) = zaj_mat(i,ib2,1)      !f           ldag(i,1) = zaj_mat(i,ib2,1)
               ldag(index_G,2) = zaj_mat(i,ib2,2)      !f           ldag(i,2) = zaj_mat(i,ib2,2)
            end do
                                                  __TIMER_DO_STOP(1335)
            zaj_mat(:,ib2,1:2) = 0.d0
                                                  __TIMER_DO_START(1336)
            do i = 1, iba2(ik)
               do j = ista_G_g1k(ik), iend_G_g1k(ik)   !f        do j = i, iba2(ik)
                  index_B = i-ista_B_g1k(ik)+1
                  index_G = j-ista_G_g1k(ik)+1
                  zaj_mat(i,ib2,1) = zaj_mat(i,ib2,1) &
                       &        + linv_wk(index_G,i,1)*ldag(index_G,1) + linv_wk(index_G,i,2)*ldag(index_G,2)
                  zaj_mat(i,ib2,2) = zaj_mat(i,ib2,2) &
                       &        + linv_wk(index_G,i,1)*ldag(index_G,2) - linv_wk(index_G,i,2)*ldag(index_G,1)
!f                  zaj_mat(i,ib2,1) = zaj_mat(i,ib2,1) + linv(j,i,1)*ldag(j,1) + linv(j,i,2)*ldag(j,2)
!f                  zaj_mat(i,ib2,2) = zaj_mat(i,ib2,2) + linv(j,i,1)*ldag(j,2) - linv(j,i,2)*ldag(j,1)
               end do
            end do
                                                  __TIMER_DO_STOP(1336)
         end do
      end if

         isize = nmatsz*np_e*kimg
                                                  __TIMER_COMM_START_w_BARRIER(mpi_ke_group,1285)
         call mpi_allreduce(MPI_IN_PLACE, zaj_mat, isize, mpi_double_precision, mpi_sum, mpi_ke_world, ierr)
                                                  __TIMER_COMM_STOP(1285)

      deallocate(linv_wk)
                                                  __TIMER_SUB_STOP(1212)
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
                                                  __TIMER_SUB_START(1213)

!!$      do ig = 1, nmatsz
      if(kimg == 2) then
         do ib = ista_e, iend_e, istep_e
!f            ib2 = map_z(ib)
            ib2 = ib - ista_e +1
#ifdef _PHASE_MULT_0_
            wfr = zaj_mat(1,ib2,1); wfi = zaj_mat(1,ib2,2)
            if(abs(wfi*wfi) > 1.d-30) then
               sqrwf = wfr**2 + wfi**2
               f = dsqrt(sqrwf)
               pcos = -wfr/f; psin = wfi/f
               if(iprimatdiagon >= 2) write(nfout,'(4x," pcos, psin = ",2e14.6)') pcos, psin
                                                  __TIMER_DO_START(1337)
               do ig = 1, iba2(ik)
                  wfr = zaj_mat(ig,ib2,1); wfi = zaj_mat(ig,ib2,2)
                  zaj_mat(ig,ib2,1) = pcos*wfr - psin*wfi
                  zaj_mat(ig,ib2,2) = pcos*wfi + psin*wfr
               end do
                                                  __TIMER_DO_STOP(1337)
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
                                                  __TIMER_DO_START(1338)
               do ig = 1, iba2(ik)
                  wfr = zaj_mat(ig,ib2,1); wfi = zaj_mat(ig,ib2,2)
                  zaj_mat(ig,ib2,1) = pcos*wfr - psin*wfi
                  zaj_mat(ig,ib2,2) = pcos*wfi + psin*wfr
               end do
                                                  __TIMER_DO_STOP(1338)
               zaj_mat(1,ib2,2) = 0.d0
!!            !$else
!!            !$   do ig = 2, kg2_gamma
!!            !$      ig1 = nbase_gamma(ig,1)
!!            !$      wfr = zaj_mat(ig1,ib2,1)
!!            !$      sqrwf = sqrwf + wfr**2
!!            !$   end do
!!            !$   sqrwf = 2.d0*sqrwf + zaj_mat(1,ib2,1)**2
!!            !$   sqrwf = 1.d0/sqrt(sqrwf)
!!            !$   do ig = 1, iba2(ik)
!!            !$      wfr = zaj_mat(ig,ib2,1)
!!            !$      zaj_mat(ig,ib2,1) = wfr*sqrwf
!!            !$      zaj_mat(ig,ib2,2) = 0.d0
!!            !$   end do
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
                                                  __TIMER_DO_START(1339)
            do ib = ista_e, iend_e, istep_e
!f               ib2 = map_z(ib)
               ib2 = ib - ista_e + 1
               do ig = 2, kg2_gamma
                  ig1 = nbase_gamma(ig,1)
                  ig2 = nbase_gamma(ig,2)
                  wfr = zaj_mat(ig1,ib2,1); wfi = zaj_mat(ig1,ib2,2)
                  zaj_mat(ig2,ib2,1) = wfr
                  zaj_mat(ig2,ib2,2) = -wfi
               end do
            end do
                                                  __TIMER_DO_STOP(1339)
         end if
      end if
      if(iprimatdiagon >= 2) write(nfout,'(" out of <<phase_mult>>")')
                                                  __TIMER_SUB_STOP(1213)
    end subroutine phase_mult

    subroutine cp_zaj_mat_to_zaj_l()
      integer :: ib,ig, i, l, j, is
      real(kind=DP),allocatable, dimension(:,:,:):: zaj_l_wk
                                                  __TIMER_SUB_START(1214)
      allocate(zaj_l_wk(kg1,np_e,kimg)) ; zaj_l_wk = 0.d0
      if(k_symmetry(ik) == GAMMA) then
                                                  __TIMER_DO_START(1340)
         do l = 1, kimg
            do ib = 1, np_e
               do i = 1, kg2_gamma
                  ig = nbase(i,ik)
                  zaj_l_wk(i,ib,l) = zaj_mat(ig,ib,l)     !f       zaj_l(i,ib,ik,l) = zaj_mat(ig,ib,l)
               end do
            end do
         end do
                                                  __TIMER_DO_STOP(1340)
      else
                                                  __TIMER_DO_START(1342)
         do l = 1, kimg
            do ib = 1, np_e
               do i = 1, iba2(ik)
                  ig = nbmat2(i,ik)
                  zaj_l_wk(ig,ib,l) = zaj_mat(i,ib,l)     !f       zaj_l(ig,ib,ik,l) = zaj_mat(i,ib,l)
               end do
            end do
         end do
                                                  __TIMER_DO_STOP(1342)
      end if
                                                  __TIMER_DO_START(1344)
   do l = 1, kimg
      zaj_l(:,:,ik,l) = 0.d0
   do j = 1, np_e
   do i = 1, np_g1k(ik)
        is = nis_g1k(myrank_g,ik) + i-1
        zaj_l(i,j,ik,l) = zaj_l_wk(is,j,l)
   enddo; enddo; enddo
                                                  __TIMER_DO_STOP(1344)

#ifdef SAVE_FFT_TIMES
  if(sw_save_fft == ON) then
     do j = 1, np_e
        status_saved_phifftr(j,ik) = OLD
     end do
  end if
#endif
    deallocate(zaj_l_wk)
                                                  __TIMER_SUB_STOP(1214)
    end subroutine cp_zaj_mat_to_zaj_l

    subroutine solve_Hx_eq_ex_ScaLAPACK()

      integer ::                          n, ne, ne_found, il, iu, info, i, j, ib, k, ib_wk, NZ
      real(kind=DP) ::                    abstol, vl, vu, sum, dlamch, ORFAC
      integer,allocatable ::              ifail(:), iwork(:),    ICLUSTR(:)
      real(kind=DP),allocatable ::        evec(:),  rwork(:),    rvmat_l(:,:)
      complex(kind=CMPLDP),allocatable :: cmat(:,:),  cwork(:)
      complex(kind=CMPLDP),parameter ::   IMG = (0.d0,1.d0)

      real(kind=DP),allocatable ::        WORK(:),GAP(:)
#ifdef _FJ_DBG0_
      integer ::                          iiii,jjjj
#endif
                                                  __TIMER_SUB_START(1471)
      zaj_mat(:,:,:) = 0.d0
!      lda = nmatsz
      !!$ vl = -1.d10 ; vu = 1.d10 ; il = ista_e ; iu = iend_e
      vl = -1.d10 ; vu = 1.d10 ; il = 1 ; iu = neg
      ne = iu - il + 1

      abstol = 2*dlamch('S')

      if (kimg == 1) then
         n = iba2(ik)
        !+++++++++++++++++++++++++++++++++++++++++++++++++
         allocate(evec(n),ifail(n))
        !+++++++++++++++++++++++++++++++++++++++++++++++++
        is = 1; ie = nmatsz

        allocate(RWORK(lwork1))
        allocate(IWORK(liwork1))
        allocate(hsmat1(lda,occ));hsmat1 = 0.0d0
        allocate(hsmat2(lda,occ));hsmat1 = 0.0d0
        allocate(rvmat_l(lda,occ)); rvmat_l = 0.0d0
        allocate(GAP(nsclrow_md*nsclcol_md))
        allocate(ICLUSTR(2*nsclrow_md*nsclcol_md))

        if(modnrm == EXECUT) then
           call trans_scalapack( block_size, hsmat, hsmat1, &
       &                        maxval(np_B_g1k), maxval(np_G_g1k), lda, occ)
           call trans_scalapack( block_size, hsmat_s, hsmat2, &
       &                        maxval(np_B_g1k), maxval(np_G_g1k), lda, occ)
           ORFAC = 0.001
           ABSTOL = 2*DLAMCH('S')
                                                  __TIMER_DO_START(1479)
        if (myrow /= -1) then
           CALL PDSYGVX( 1,'V', 'I', 'U', nmatsz_opt, hsmat1, 1, 1, DESCA, hsmat2, 1, 1,                   &
     &                 DESCA, VL, VU, 1, neg, ABSTOL, ne_found, NZ, evec, ORFAC, rvmat_l, 1, 1, &
     &                 DESCZ, RWORK, lwork1, IWORK, liwork1, IFAIL,                   &
     &                 ICLUSTR, GAP, INFO)
        endif
                                                  __TIMER_DO_STOP(1479)
        else
           call trans_scalapack( block_size, hsmat, hsmat1, &
       &                        maxval(np_B_g1k), maxval(np_G_g1k), lda, occ)
           ORFAC = 0.001
           ABSTOL = 2*DLAMCH('S')
                                                  __TIMER_DO_START(1480)
        if (myrow /= -1) then
           CALl PDSYEVX( 'V', 'I', 'U', nmatsz_opt, hsmat1, 1, 1, DESCA, VL,                &
     &                    VU, 1, neg, ABSTOL, ne_found, NZ, evec, ORFAC, rvmat_l, 1,        &
     &                    1, DESCZ, RWORK, lwork1, IWORK, liwork1, IFAIL,               &
     &                    ICLUSTR, GAP, INFO )
        endif
                                                  __TIMER_DO_STOP(1480)
        endif

        call trans_scalapack_r( block_size, rvmat_l, zaj_mat, &
     &                     maxval(np_B_g1k), maxval(np_G_g1k), lda, occ)

            if (nrank_e*nrank_g /= nsclcol_md*nsclrow_md) then
               call eigsend(evec,n)
            end if

         if(iprimatdiagon >= 2) then
            write(nfout,'(" --- eigen values (solve_Hx_eq_ex_ScaLAPACK): ik = ",i8)') ik
            write(nfout,'(8f10.6)') (evec(i),i=1,ne)
!!$            a = "     "
!!$            do ib = ista_e, iend_e, istep_e
!f               write(nfout,'(" evec(",i4,",",i3,")= ",e14.6," ",a4,5e14.6)') &
!f                    & ib,ik,evec(ib),a,(rvmat(i,ib),i=1,5)
!x               write(nfout,'(" evec(",i4,",",i3,")= ",e14.6," ",a4,5e14.6)') &
!x                    & ib,ik,evec(ib),a,(amat(i,ib),i=1,5)
!!$            end do
         end if
!x         if (ne_found /= ne) then
!x            write(*,*) '### ERROR ### ne_found != ne'
!x            write(*,*) '   ne_found ...',ne_found
!x            write(*,*) '   ne       ...',ne
!x            stop
!x         end if
         if (info /= 0) then
            write(*,*) '### ERROR ### info /= 0'
            write(*,*) '   info ...',info
            call phase_error_with_msg(nfout,'### ERROR ### info /= 0',__LINE__,__FILE__)
         end if
                                                  __TIMER_DO_START(1481)
         do ib = ista_e, iend_e, istep_e
!f            j = map_z(ib)
            j = ib - ista_e +1
            ib_wk = neg_g(j)
            do i = 1,n
!f               zaj_mat(i,j,1) = rvmat(i,ib)
!f               zaj_mat(i,j,1) = rvmat(i,ib_wk)
!aa               zaj_mat(i,j,1) = amat(i,ib_wk)
               !!$ zaj_l(i,j,ik,1) = rvmat(i,ib-ista_e+1) ! M.Okamoto (August 22,2003)
            end do
!f            eko_l(j,ik) = evec(ib)
            eko_l(j,ik) = evec(ib_wk)
            !!$ eko_l(j,ik) = evec(ib-ista_e+1) ! M.Okamoto (August 22,2003)
         end do
                                                  __TIMER_DO_STOP(1481)
         if(af/=0) eko_l(:,ik+1) = eko_l(:,ik)
        !+++++++++++++++++++++++++++++
         deallocate(RWORK)
         deallocate(hsmat1,hsmat2)
         deallocate(rvmat_l)
         deallocate(GAP,ICLUSTR)
        !+++++++++++++++++++++++++++++
      else if (kimg == 2) then
         n = iba2(ik)
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         allocate(evec(n),ifail(n))
         evec = 0.d0;  ifail = 0
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         is = 1; ie = nmatsz

        allocate(CWORK(lwork2))
        allocate(RWORK(lrwork2))
        allocate(IWORK(liwork2))
        allocate(hsmat1(lda*kimg,occ)); hsmat1 = 0.0d0
        allocate(hsmat2(lda*kimg,occ)); hsmat2 = 0.0d0
        allocate(rvmat_l(lda*kimg,occ)); rvmat_l = 0.0d0
        allocate(GAP(nsclrow_md*nsclcol_md))
        allocate(ICLUSTR(2*nsclrow_md*nsclcol_md))

        if(modnrm == EXECUT) then
           call trans_scalapack( block_size, hsmat, hsmat1, &
       &                        maxval(np_B_g1k), maxval(np_G_g1k), lda, occ)
           call trans_scalapack( block_size, hsmat_s, hsmat2, &
       &                     maxval(np_B_g1k), maxval(np_G_g1k), lda, occ)
           ORFAC = 0.001
           ABSTOL = 2*DLAMCH('S')
                                                  __TIMER_DO_START(1482)
        if (myrow /= -1) then
           CALL PZHEGVX( 1,'V', 'I', 'U', nmatsz_opt, hsmat1, 1, 1, DESCA, hsmat2, 1, 1,                   &
     &              DESCA, VL, VU, 1, neg, ABSTOL, ne_found, NZ, evec, ORFAC, rvmat_l, 1, 1, &
     &              DESCZ, CWORK, lwork2, RWORK, lrwork2, IWORK, liwork2, IFAIL,                   &
     &              ICLUSTR, GAP, INFO)
        endif
                                                  __TIMER_DO_STOP(1482)
        else
           call trans_scalapack( block_size, hsmat, hsmat1, &
       &                        maxval(np_B_g1k), maxval(np_G_g1k), lda, occ)
                                                  __TIMER_DO_START(1483)
           ORFAC = 0.001
           ABSTOL = 2*DLAMCH('S')
        if (myrow /= -1) then
           CALL PZHEEVX( 'V', 'I', 'U', nmatsz_opt, hsmat1, 1, 1, DESCA, VL,                         &
     &                    VU, 1, neg, ABSTOL, ne_found, NZ, evec, ORFAC, rvmat_l, 1,            &
     &                    1, DESCZ, CWORK, lwork2, RWORK, lrwork2, IWORK,            &
     &                    liwork2, IFAIL, ICLUSTR, GAP, INFO )
        endif
                                                  __TIMER_DO_STOP(1483)
        endif

        call trans_scalapack_r( block_size, rvmat_l, zaj_mat, &
     &                     maxval(np_B_g1k), maxval(np_G_g1k), lda, occ)

            if (nrank_e*nrank_g /= nsclcol_md*nsclrow_md) then
               call eigsend(evec,n)
            end if

#ifdef _FJ_DBG0_
       if(mype .eq. 0) then
            write(1002,'("evec")')
            write(1002,'(5f10.7)') (evec(iiii),iiii=1,neg)

            write(1002,'("zaj_mat:kimg=1")')
            do iiii = 1, nb_mgs
               write(1002,'("iiii=",i4)') iiii
               write(1002,'(8(f10.7))') (zaj_mat(jjjj,iiii,1),jjjj=1,200)
            enddo

            write(1002,'("zaj_mat:kimg=2")')
            do iiii = 1, nb_mgs
               write(1002,'("iiii=",i4)') iiii
               write(1002,'(8(f10.7))') (zaj_mat(jjjj,iiii,2),jjjj=1,200)
            enddo
!            write(1002,'("rvmat_l:kimg=1")')
!            do iiii = 1, occ
!               write(1002,'("iiii=",i4)') iiii
!               write(1002,'(8(f10.8))') (rvmat_l(jjjj,iiii),jjjj=1,lda*kimg,1)
!            enddo
!            write(1002,'("rvmat_l:kimg=2")')
!            do iiii = 1, occ
!               write(1002,'("iiii=",i4)') iiii
!               write(1002,'(8(f10.8))') (rvmat_l(jjjj,iiii),jjjj=1,lda*kimg,2)
!            enddo
        endif
#endif

         if(iprimatdiagon >= 2) then
            write(nfout,'(" --- eigen values (solve_Hx_eq_ex_ScaLAPACK): ik = ",i8)') ik
            write(nfout,'(8f10.6)') (evec(i),i=1,ne)
         end if
!x         if (ne_found /= ne) then
!x            write(*,*) '### ERROR ### ne_found != ne'
!x            write(*,*) '   ne_found ...',ne_found
!x            write(*,*) '   ne       ...',ne
!x            stop
!x         end if
      if (myrow /= -1) then
         if (info /= 0) then
            write(*,*) '### ERROR ### info /= 0'
            write(*,*) '   info ...',info
            call phase_error_with_msg(nfout,'### ERROR ### info /= 0',__LINE__,__FILE__)
         end if
      endif
                                                  __TIMER_DO_START(1484)
         do ib = ista_e, iend_e, istep_e
            j = ib - ista_e +1              !f      j = map_z(ib)
            ib_wk = neg_g(j)
            eko_l(j,ik) = evec(ib_wk)       !f      eko_l(j,ik) = evec(ib)
            !!$ eko_l(j,ik) = evec(ib-ista_e+1) ! M.Okamoto (August 22,2003)
         end do
                                                  __TIMER_DO_STOP(1484)
           if(af/=0) eko_l(:,ik+1) = eko_l(:,ik)
        !+++++++++++++++++++++++++++++++++++
          deallocate(RWORK,CWORK)
          deallocate(hsmat1,hsmat2)
          deallocate(rvmat_l)
          deallocate(GAP,ICLUSTR)
        !+++++++++++++++++++++++++++++++++++
      end if
      if(iprimatdiagon >= 2) call wd_CalculatedEigenvalues(nfout,ik,n,ne,evec,"solve_Hx_eq_ex")

      sum = 0.d0
      do i = 1,ne
         sum = sum + abs(evec(i))
      end do
      if (sum < SmallestPositiveNumber*1.d5) call phase_error_with_msg(nfout,' ! illegal evec (solve_Hx_eq_ex)'&
                                                                      ,__LINE__,__FILE__)
      if (iprimatdiagon >= 2) write(nfout,'(" -- sum = ",f20.8)') sum
     !+++++++++++++++++++++++++++++
      deallocate(evec,iwork,ifail)
     !+++++++++++++++++++++++++++++
                                                  __TIMER_SUB_STOP(1471)
    end subroutine solve_Hx_eq_ex_ScaLAPACK

  subroutine set_nprow_npcol(nsclrow_md,nsclcol_md)
    integer, intent(inout) :: nsclrow_md,nsclcol_md

!      if(printable) write(nfout,'("set nprow,npcol,block_size=",3i5)') nsclrow_md,nsclcol_md,block_size
      if (iprimatdiagon>=2) write(nfout,'("set nprow,npcol,block_size=",3i5)') nsclrow_md,nsclcol_md,block_size;call flush(nfout)
!      nprow = nrank_g
!      npcol = nrank_e
!      nprow = 1
!      npcol = nrank_e*nrank_g

!   if(nprow /= nrank_e) then
!      nprow = 1
!      npcol = nrank_e
!      if(printable) write(nfout,'("set nprow,npcol=",2i5)') nprow,npcol
!   else
!       if(iprisubmat >= 2) write(nfout,'("nprow,npcol=",2i5)') nprow,npcol
!   end if
  end subroutine set_nprow_npcol

  subroutine scalapack_setup_3D(ndim,lwork1,lrwork1,liwork1,lwork2,lrwork2,liwork2,lda,occ,ictxt,myrow,mycol,desca,descz,usermap)
    integer, intent(in) :: ndim
    integer, intent(out) :: lwork1,lrwork1,liwork1,lwork2,lrwork2,liwork2,lda,occ,ictxt,myrow,mycol
    integer, dimension(9), intent(out) :: desca,descz
    integer, intent(inout) :: usermap(nprow_md,npcol_md)

    integer :: i,j, nb, np0,nq0, iam,nprocs,info
    integer, external :: indxg2l, indxg2p, numroc, iceil, pjlaenv
    integer, parameter :: BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1, &
     &                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
     &                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9
    integer :: id_sname = -1
!fj+
    integer     :: lwork0,lrwork0,liwork0
    integer     :: itmp(1)
    real(DP)    :: rtmp(1)
    complex(DP) :: ctmp(1)
    complex(DP) :: htmp1(1),htmp2(1),rvtmp1(1)
    real(DP)    :: htmp1r(1),htmp2r(1),rvtmp1r(1)
    real(DP)    :: GAP(1),evec(1)
    integer     :: ICLUSTR(1),IFAIL(1),NZ,ne_found
    real(DP)    :: ORFAC,ABSTOL,VL,VU
!fj-
    integer      :: nn
    integer :: sys2blacs_handle
                                                  __TIMER_SUB_START(1472)
    call tstatc0_begin('scalapack_setup ', id_sname)

    lwork1  = 0
    lrwork1 = 0
    liwork1 = 0
    lwork2  = 0
    lrwork2 = 0
    liwork2 = 0
    desca = 0
    descz = 0

    nb = block_size

    if(mod(ndim,nb)>0) then
       lda = ndim/nb +1
    else
       lda = ndim/nb
    end if
    occ = lda
    if(mod(lda,nprow_md)>0) then
       lda = lda/nprow_md+1
    else
       lda = lda/nprow_md
    end if
    if(mod(occ,npcol_md)>0) then
       occ = occ/npcol_md+1
    else
       occ = occ/npcol_md
    end if
    lda = lda*nb
    occ = occ*nb
!    if(iprisubmat>=2) write(nfout,'("lda,occ=",2i5)') lda,occ
    do j=1,npcol_md
       do i=1,nprow_md
!fj --------------------
!fj       usermap(i,j) = myrank_k*nrank_e + (j-1)*nprow + i - 1
    !      usermap(i,j) = myrank_k*nrank_e*nrank_g + (j-1)*nprow_md + i - 1
          usermap(i,j) = myrank_k*nrank_e*nrank_g + (j-1)*nprow_md + i - 1
!fj --------------------
       end do
    end do
!    if(iprimatdiagon>=1) then
!       write(nfout,'("USERMPAP")')
!       do i=1,nprow
!          write(nfout,'(1x,i3)') (usermap(i,j),j=1,npcol)
!       end do
!    end if
    call blacs_get(-1,0,ictxt)
!    ictxt = sys2blacs_handle(mpi_k_world(myrank_k))
!    call blacs_setup(iam,nprocs)
    call blacs_gridmap(ictxt,usermap,nprow_md,nprow_md,npcol_md)
    call blacs_gridinfo(ictxt,nprow_md,npcol_md,myrow,mycol)
    if(myrow == -1) then
!       if(printable) write(nfout,*) 'BLACS init failed.'
!       stop 'BLACS init failed.'
      return
    else
    end if
    call descinit(desca,ndim,ndim,nb,nb,0,0,ictxt,lda,info)
    call descinit(descz,ndim,ndim,nb,nb,0,0,ictxt,lda,info)

 if(kimg==1) then
       if(method_scalapack==HOUSEHOLDER) then
          if(iprimatdiagon>=2) write(nfout,*) 'HOUSEHOLDER:yes'
          nn  = max(ndim, nb, 2)
          np0 = NUMROC( max(ndim,nb,2), nb, 0, 0, nprow_md )
          nq0 = NUMROC( max(ndim,nb,2), nb, 0, 0, npcol_md )
          if(iprimatdiagon>=2) write(nfout,*) 'np0,nq0=',np0,nq0
          lwork1 = 5*ndim + max(5*nn, np0*nq0+2*nb*nb) + iceil(ndim,nprow_md*npcol_md)*nn
          liwork1= 6*max(ndim,nprow_md*npcol_md+1,4)
          if(iprimatdiagon>=2) write(nfout,*) 'lwork1 =',lwork1
          if(iprimatdiagon>=2) write(nfout,*) 'liwork1=',liwork1
       else
          if(iprimatdiagon>=2) write(nfout,*) 'HOUSEHOLDER:no'
          nn  = max(ndim, nb, 2)
          np0 = NUMROC( max(ndim,nb,2), nb, 0, 0, nprow_md )
          nq0 = NUMROC( max(ndim,nb,2), nb, 0, 0, npcol_md )
          if(iprimatdiagon>=2) write(nfout,*) 'np0,nq0=',np0,nq0
          lwork0 = 5*ndim + max(5*nn, np0*nq0+2*nb*nb) + iceil(ndim,nprow_md*npcol_md)*nn
          liwork0= 6*max(ndim,nprow_md*npcol_md+1,4)
          if(iprimatdiagon>=2) write(nfout,*) 'lwork0 =',lwork0
          if(iprimatdiagon>=2) write(nfout,*) 'liwork0=',liwork0
          if(modnrm == EXECUT) then
              CALL PDSYGVX( 1,'V', 'I', 'U', ndim, htmp1r, 1, 1, DESCA, htmp2r, 1, 1,                   &
     &                 DESCA, VL, VU, 1, neg, ABSTOL, ne_found, NZ, evec, ORFAC, rvtmp1r, 1, 1, &
     &                 DESCZ, rtmp, -1, itmp, -1, IFAIL,                   &
     &                 ICLUSTR, GAP, INFO)
          else
              CALl PDSYEVX( 'V', 'I', 'U', ndim, htmp1r, 1, 1, DESCA, VL,                &
     &                    VU, 1, neg, ABSTOL, ne_found, NZ, evec, ORFAC, rvtmp1r, 1,        &
     &                    1, DESCZ, rtmp, -1, itmp, -1, IFAIL,               &
     &                    ICLUSTR, GAP, INFO )
          endif
          lwork1  = max( lwork0, nint(rtmp(1)) )
          liwork1 = max( liwork0, itmp(1) )

          if(iprimatdiagon>=2) write(nfout,*) 'lrwork1=',lrwork1
          if(iprimatdiagon>=2) write(nfout,*) 'liwork1=',liwork1

       end if
 else
       if(method_scalapack==HOUSEHOLDER) then
          if(iprimatdiagon>=2) write(nfout,*) 'HOUSEHOLDER:yes'
          nn  = max(ndim, nb, 2)
          np0 = NUMROC( max(ndim,nb,2), nb, 0, 0, nprow_md )
          nq0 = NUMROC( max(ndim,nb,2), nb, 0, 0, npcol_md )
          if(iprimatdiagon>=2) write(nfout,*) 'np0,nq0=',np0,nq0
          lwork2 = ndim + ( np0 + nq0 + nb ) * nb
          lrwork2 = 4*ndim + max(5*nn, np0*nq0) + iceil(ndim,nprow_md*npcol_md)*nn
          liwork2 = 6*MAX( ndim, NPROW_md*NPCOL_md + 1, 4 )
          if(iprimatdiagon>=2) write(nfout,*) 'lwork2 =',lwork2
          if(iprimatdiagon>=2) write(nfout,*) 'lrwork2=',lrwork2
          if(iprimatdiagon>=2) write(nfout,*) 'liwork2=',liwork2
       else
          if(iprimatdiagon>=2) write(nfout,*) 'HOUSEHOLDER:no'
          nn  = max(ndim, nb, 2)
          np0 = NUMROC( max(ndim,nb,2), nb, 0, 0, nprow_md )
          nq0 = NUMROC( max(ndim,nb,2), nb, 0, 0, npcol_md )
          if(iprimatdiagon>=2) write(nfout,*) 'np0,nq0=',np0,nq0
          lwork0 = ndim + ( np0 + nq0 + nb ) * nb
          lrwork0 = 4*ndim + max(5*nn, np0*nq0) + iceil(ndim,nprow_md*npcol_md)*nn
          liwork0 = 6*MAX( ndim, NPROW_md*NPCOL_md + 1, 4 )
          if(iprimatdiagon>=2) write(nfout,*) 'lwork0 =',lwork0
          if(iprimatdiagon>=2) write(nfout,*) 'lrwork0=',lrwork0
          if(iprimatdiagon>=2) write(nfout,*) 'liwork0=',liwork0

          if(modnrm == EXECUT) then
              CALL PZHEGVX( 1,'V', 'I', 'U', ndim, htmp1, 1, 1, DESCA, htmp2, 1, 1,                   &
     &              DESCA, VL, VU, 1, neg, ABSTOL, ne_found, NZ, evec, ORFAC, rvtmp1, 1, 1, &
     &              DESCZ, ctmp, -1, rtmp, -1, itmp, -1, IFAIL,                   &
     &              ICLUSTR, GAP, INFO)
          else
              CALL PZHEEVX( 'V', 'I', 'U', ndim, htmp1, 1, 1, DESCA, VL,                         &
     &                    VU, 1, neg, ABSTOL, ne_found, NZ, evec, ORFAC, rvtmp1, 1,            &
     &                    1, DESCZ, ctmp, -1, rtmp, -1, itmp,            &
     &                    -1, IFAIL, ICLUSTR, GAP, INFO )
          endif

          lwork2  = max( lwork0 , nint(real(ctmp(1))) )
          lrwork2 = max( lrwork0, nint(rtmp(1)) )
          liwork2 = max( liwork0, itmp(1) )
       end if
!-- _1x10
!          lrwork2 = lrwork2*100
!-- _1x10_1
!          lwork2  = lwork2*10
!          lrwork2 = lrwork2*10
!-- _1x10_2
!          lwork2  = lwork2*10
!          lrwork2 = lrwork2*100
!-- _1x10_3
!          lwork2  = lwork2*10
!          lrwork2 = lrwork2*1000
!-- _1x10_4
!          lwork2  = lwork2*100
!          lrwork2 = lrwork2*10
!-- _1x10_5
!          lwork2  = lwork2*100
!          lrwork2 = lrwork2*100
!-- _1x10_6
!          lwork2  = lwork2*100
!          lrwork2 = lrwork2*1000
!-- _1x10_7
!          lwork2  = lwork2*1000
!          lrwork2 = lrwork2*10
!-- _1x10_8
!          lwork2  = lwork2*1000
!          lrwork2 = lrwork2*100
!-- _1x10_9
!          lwork2  = lwork2*1000
!          lrwork2 = lrwork2*1000
!--
!           lwork2  = lwork2*10
!           lrwork2 = lrwork2*10
!!$           lrwork2 = lrwork2*1000
!          lrwork2 = lrwork2*10
!          lrwork2 = lrwork2*100
!          liwork2 = liwork2*10
!--
          if(iprimatdiagon>=2) write(nfout,*) 'lwork2 =',lwork2
          if(iprimatdiagon>=2) write(nfout,*) 'lrwork2=',lrwork2
          if(iprimatdiagon>=2) write(nfout,*) 'liwork2=',liwork2

 end if
 call flush(nfout)
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(1472)
  end subroutine scalapack_setup_3D

    subroutine trans_scalapack(scl_nb,zmat_l,amat_l,maxeg,maxe,lda,occ)
     integer, intent(in) :: scl_nb, maxeg, maxe, lda, occ
     real(kind=DP), intent(in) :: zmat_l(maxeg,maxe,kimg)
     real(kind=DP), intent(inout) :: amat_l(lda*kimg,occ)

     integer , allocatable, dimension(:) :: wk_rk
     integer :: i, k, j,  nn, lrk
     integer :: lno_r, ng_r, lrk_r, rank_r
     integer :: lno_c, ng_c, lrk_c, rank_c

     real(kind=DP), allocatable, dimension(:,:) :: send_buf, recv_buf

     integer, allocatable, dimension(:)   :: req_r, req_s
     integer, allocatable, dimension(:,:) :: sta_r, sta_s
     integer :: lrank, itag, icnt_send, icnt_recv, mpi_comm, ierr
                                                  __TIMER_SUB_START(1474)

      mpi_comm = mpi_k_world(myrank_k)
      itag = 10001
      allocate(req_r(scl_md_comm_rank_r), stat=ierr)
      allocate(req_s(scl_md_comm_rank  ), stat=ierr)
      allocate(sta_r(MPI_STATUS_SIZE,scl_md_comm_rank_r), stat=ierr)
      allocate(sta_s(MPI_STATUS_SIZE,scl_md_comm_rank  ), stat=ierr)
      allocate(recv_buf(scl_md_comm_max*kimg,scl_md_comm_rank_r), stat=ierr)
      allocate(send_buf(scl_md_comm_max*kimg,scl_md_comm_rank  ), stat=ierr)
      icnt_recv = 0
      icnt_send = 0

                                                  __TIMER_DO_START(1485)
      do i = 1, scl_md_comm_rank_r
         lrank = scl_md_comm_rno_r(i)
         if (lrank /= mype) then
            icnt_recv = icnt_recv + 1
            call mpi_irecv(recv_buf(1,i), scl_md_comm_max*kimg, &
           &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
             if (ierr /= 0) then
                write(nfout,*)' trans_scalapack :  mpi_irecv error'
                call flush(nfout)
                call mpi_abort(mpi_comm_world, 10001, ierr)
             endif
         endif
      enddo
                                                  __TIMER_DO_STOP(1485)

      allocate(wk_rk(max(scl_md_comm_rank,scl_md_comm_rank_r)))
      wk_rk = 0
      rank_c = myrank_g
      rank_r = myrank_e
                                                  __TIMER_DO_START(1486)
      do k = nis_G_g1k(rank_c,ik), nie_G_g1k(rank_c,ik)
         if (mod(k,scl_nb) > 0) then
            nn = k / scl_nb + 1
         else
            nn = k / scl_nb
         end if
         lrk_c = mod((nn-1),scl_col_md)
         do j = nis_B_g1k(rank_r,ik), nie_B_g1k(rank_r,ik)
            if (mod(j,scl_nb) > 0) then
               nn = j / scl_nb + 1
            else
               nn = j / scl_nb
            end if
            lrk_r = mod((nn-1),scl_row_md)
            lrk = scl_row_md*lrk_c+lrk_r
!!          if (lrk /= mype) then
               do lrank = 1, scl_md_comm_rank
                  if (scl_md_comm_rno(lrank) == lrk) then
                     wk_rk(lrank) = wk_rk(lrank) + 1
                     if (kimg == 1) then
                        send_buf(wk_rk(lrank),lrank) = &
 &                                                 zmat_l(j-nis_B_g1k(rank_r,ik)+1,k-nis_G_g1k(rank_c,ik)+1,1)
                     else
                        send_buf(wk_rk(lrank)*2-1,lrank) = &
 &                                                 zmat_l(j-nis_B_g1k(rank_r,ik)+1,k-nis_G_g1k(rank_c,ik)+1,1)
                        send_buf(wk_rk(lrank)*2  ,lrank) = &
 &                                                 zmat_l(j-nis_B_g1k(rank_r,ik)+1,k-nis_G_g1k(rank_c,ik)+1,2)
                     end if
                     exit
                  end if
               end do
!!          end if
         end do
      end do
                                                  __TIMER_DO_STOP(1486)

                                                  __TIMER_COMM_START_w_BARRIER(mpi_comm,1489)
                                                  __TIMER_DO_START(1485)
      do i = 1, scl_md_comm_rank
         lrank = scl_md_comm_rno(i)
         if (lrank /= mype) then
            icnt_send = icnt_send + 1
            call mpi_isend(send_buf(1,i), scl_md_comm_max*kimg, &
           &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
             if (ierr /= 0) then
                write(nfout,*)' trans_scalapack :  mpi_isend error'
                call flush(nfout)
                call mpi_abort(mpi_comm_world, 10002, ierr)
             endif
         endif
      enddo
                                                  __TIMER_DO_STOP(1485)

   if ((mype-nrank_e*nrank_g*myrank_k) < (scl_col_md*scl_row_md)) then
      rank_r = mod((mype-nrank_e*nrank_g*myrank_k),scl_row_md)
      rank_c = (mype-nrank_e*nrank_g*myrank_k)/scl_row_md
                                                  __TIMER_DO_START(1487)
      do k = nis_col_md(rank_c), nie_col_md(rank_c)
         ng_c = nmatsz_col(k)
         lrk_c = map_G_g1k(ng_c,ik)
         lno_c = ng_c - nis_G_g1k(lrk_c,ik) + 1
         do j = nis_row_md(rank_r), nie_row_md(rank_r)
            ng_r = nmatsz_row(j)
            lrk_r = map_B_g1k(ng_r,ik)
            lno_r = ng_r - nis_B_g1k(lrk_r,ik) + 1
            lrk = nrank_e*lrk_c+lrk_r
            if (lrk == mype) then

               if(ng_r> ng_c) cycle

               if (kimg == 1) then
                  amat_l(j-nis_row_md(rank_r)+1,k-nis_col_md(rank_c)+1) = zmat_l(lno_r,lno_c,1)
               else
                  amat_l((j-nis_row_md(rank_r)+1)*2-1,k-nis_col_md(rank_c)+1) = zmat_l(lno_r,lno_c,1)
                  amat_l((j-nis_row_md(rank_r)+1)*2  ,k-nis_col_md(rank_c)+1) = zmat_l(lno_r,lno_c,2)
               end if

            end if
         end do
      end do
                                                  __TIMER_DO_STOP(1487)
   end if

      call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
       if (ierr /= 0) then
          write(nfout,*)' trans_scalapack :  mpi_waitall error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 10003, ierr)
       endif

      call mpi_waitall(icnt_send, req_s, sta_s, ierr)
       if (ierr /= 0) then
          write(nfout,*)' trans_scalapack :  mpi_waitall error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 10004, ierr)
       endif
                                                  __TIMER_COMM_STOP(1489)

      if ((mype-nrank_e*nrank_g*myrank_k) < (scl_col_md*scl_row_md)) then

         wk_rk = 0
         rank_r = mod((mype-nrank_e*nrank_g*myrank_k),scl_row_md)
         rank_c = (mype-nrank_e*nrank_g*myrank_k)/scl_row_md
                                                  __TIMER_DO_START(1488)
         do k = nis_col_md(rank_c), nie_col_md(rank_c)
            ng_c = nmatsz_col(k)
            lrk_c = map_G_g1k(ng_c,ik)
            do j = nis_row_md(rank_r), nie_row_md(rank_r)
               ng_r = nmatsz_row(j)
               lrk_r = map_B_g1k(ng_r,ik)
               lrk = nrank_e*lrk_c+lrk_r
!!             if (lrk /= mype) then
                  do lrank = 1, scl_md_comm_rank_r
                     if (scl_md_comm_rno_r(lrank) == lrk) then
                        wk_rk(lrank) = wk_rk(lrank) + 1

                        if(ng_r> ng_c) cycle

                        if (kimg == 1) then
                   amat_l(j-nis_row_md(rank_r)+1,k-nis_col_md(rank_c)+1) = recv_buf(wk_rk(lrank),lrank)
                        else
                   amat_l((j-nis_row_md(rank_r)+1)*2-1,k-nis_col_md(rank_c)+1) = recv_buf(wk_rk(lrank)*2-1,lrank)
                   amat_l((j-nis_row_md(rank_r)+1)*2  ,k-nis_col_md(rank_c)+1) = recv_buf(wk_rk(lrank)*2  ,lrank)
                        end if

                        exit
                     end if
                  end do
!!             end if
            end do
         end do
                                                  __TIMER_DO_STOP(1488)
      end if

      deallocate(req_r, stat=ierr)
      deallocate(req_s, stat=ierr)
      deallocate(sta_r, stat=ierr)
      deallocate(sta_s, stat=ierr)
      deallocate(send_buf, stat=ierr)
      deallocate(recv_buf, stat=ierr)
      deallocate(wk_rk, stat=ierr)
                                                  __TIMER_SUB_STOP(1474)

    end subroutine trans_scalapack

    subroutine trans_scalapack_r(scl_nb,zmat_l,z_mat,maxeg,maxe,lda,occ)
     integer, intent(in) :: scl_nb, maxeg, maxe, lda, occ
     real(kind=DP), intent(in) :: zmat_l(lda*kimg,occ)
     real(kind=DP), intent(inout) :: z_mat(nmatsz,np_e,kimg)

     integer , allocatable, dimension(:) :: wk_rk
     integer :: i, k, j, nn, mm, j_ng
     integer :: lrk, lad
     integer :: lno_r, lrk_r, rank_r
     integer :: lno_c, lrk_c, rank_c
     integer :: lrk2, lrk_r2, lrk_c2
     integer :: lrk1, lrk_r1, lrk_c1

     real(kind=DP), allocatable, dimension(:,:) :: send_buf, recv_buf
     real(kind=DP), allocatable, dimension(:,:,:) :: amat_l

     integer, allocatable, dimension(:)   :: req_r, req_s
     integer, allocatable, dimension(:,:) :: sta_r, sta_s
     integer :: ib, lrank, itag, icnt_send, icnt_recv, mpi_comm, ierr

      mpi_comm = mpi_k_world(myrank_k)
      itag = 10002
                                                  __TIMER_SUB_START(1475)
      allocate(amat_l(nmatsz,np_e,kimg), stat=ierr); amat_l = 0.0d0
      allocate(req_r(scl_md2_comm_rank_r), stat=ierr)
      allocate(req_s(scl_md2_comm_rank  ), stat=ierr)
      allocate(sta_r(MPI_STATUS_SIZE,scl_md2_comm_rank_r), stat=ierr)
      allocate(sta_s(MPI_STATUS_SIZE,scl_md2_comm_rank  ), stat=ierr)
      allocate(recv_buf(scl_md2_comm_max_r*kimg,scl_md2_comm_rank_r), stat=ierr)
      allocate(send_buf(scl_md2_comm_max_r*kimg,scl_md2_comm_rank  ), stat=ierr)
      icnt_recv = 0
      icnt_send = 0
!     send_buf=0.0d0
!     recv_buf=0.0d0

                                                  __TIMER_DO_START(1490)
      do i = 1, scl_md2_comm_rank_r
         lrank = scl_md2_comm_rno_r(i)
         if (lrank /= mype) then
            icnt_recv = icnt_recv + 1
            call mpi_irecv(recv_buf(1,i), scl_md2_comm_max_r*kimg, &
           &               mpi_double_precision, lrank, itag, mpi_comm, req_r(icnt_recv), ierr)
             if (ierr /= 0) then
                write(nfout,*)' trans_scalapack_r :  mpi_irecv error'
                call flush(nfout)
                call mpi_abort(mpi_comm_world, 10001, ierr)
             endif
         endif
      enddo
                                                  __TIMER_DO_STOP(1490)

      allocate(wk_rk(max(scl_md2_comm_rank,scl_md2_comm_rank_r)))

      if ((mype-nrank_e*nrank_g*myrank_k) < (scl_col_md*scl_row_md)) then

         wk_rk = 0

         rank_r = mod((mype-nrank_e*nrank_g*myrank_k),scl_row_md)
         rank_c = (mype-nrank_e*nrank_g*myrank_k)/scl_row_md
                                                  __TIMER_DO_START(1491)
         do k = 1, neg
            ib = neg_g_all(k)
            if (mod(ib,nb_mgs) > 0) then
               nn = ib / nb_mgs + 1
            else
               nn = ib / nb_mgs
            end if
            lrk_c1 = mod((nn-1),nrank_e)
            lrk_c2 = mod((ib-1)/scl_nb,scl_col_md)
            mm = (ib-1)/(scl_col_md*scl_nb)
            lad = mod(ib-1,scl_nb)+1
            lno_c = mm*scl_nb+lad
            do j = nis_row_md(rank_r), nie_row_md(rank_r)
               j_ng = nmatsz_row(j)
               lrk_r1 = map_G_g1k(j_ng,ik)
               lrk_r2 = mod((j_ng-1)/scl_nb,scl_row_md)
               nn = (j_ng-1)/(scl_row_md*scl_nb)
               lad = mod(j_ng-1,scl_nb)+1
               lno_r = nn*scl_nb+lad
               lrk1 = nrank_e*lrk_r1+lrk_c1
               lrk2 = scl_row_md*lrk_c2+lrk_r2
               if (lrk2 == mype) then
                  do lrank = 1, scl_md2_comm_rank
                     if(scl_md2_comm_rno(lrank) == lrk1) then
                        wk_rk(lrank) = wk_rk(lrank) + 1
                        if (kimg == 1) then
                           send_buf(wk_rk(lrank)    ,lrank) =   zmat_l(lno_r,lno_c)
                        else
                           send_buf(wk_rk(lrank)*2-1,lrank) =   zmat_l(lno_r*2-1,lno_c)
                           send_buf(wk_rk(lrank)*2  ,lrank) =   zmat_l(lno_r*2  ,lno_c)
                        endif
                        exit
                     end if
                  end do
               end if
            end do
         end do
                                                  __TIMER_DO_STOP(1491)
      end if

                                                  __TIMER_COMM_START_w_BARRIER(mpi_comm,1494)
                                                  __TIMER_DO_START(1490)
      do i = 1, scl_md2_comm_rank
         lrank = scl_md2_comm_rno(i)
         if (lrank /= mype) then
            icnt_send = icnt_send + 1
            call mpi_isend(send_buf(1,i), scl_md2_comm_max_r*kimg, &
           &               mpi_double_precision, lrank, itag, mpi_comm, req_s(icnt_send), ierr)
             if (ierr /= 0) then
                write(nfout,*)' trans_scalapack :  mpi_isend error'
                call flush(nfout)
                call mpi_abort(mpi_comm_world, 10002, ierr)
             endif
         endif
      enddo
                                                  __TIMER_DO_STOP(1490)

      rank_c = myrank_e
      rank_r = myrank_g
                                                  __TIMER_DO_START(1492)
      do k = nis_e(rank_c), nie_e(rank_c)
         ib = neg_g_all(k)
         lrk_c = mod((ib-1)/scl_nb,scl_col_md)
         mm = (ib-1)/(scl_col_md*scl_nb)
         lad = mod(ib-1,scl_nb)+1
         lno_c = mm*scl_nb+lad
         do j = nis_G_g1k(rank_r,ik), nie_G_g1k(rank_r,ik)
            lrk_r = mod((j-1)/scl_nb,scl_row_md)
            nn = (j-1)/(scl_row_md*scl_nb)
            lad = mod(j-1,scl_nb)+1
            lno_r = nn*scl_nb+lad
            lrk = scl_row_md*lrk_c+lrk_r
            if (lrk == mype) then
               if (kimg == 1) then
                  amat_l(j,k-nis_e(rank_c)+1,1) = zmat_l(lno_r,lno_c)
               else
                  amat_l(j,k-nis_e(rank_c)+1,1) = zmat_l(lno_r*2-1,lno_c)
                  amat_l(j,k-nis_e(rank_c)+1,2) = zmat_l(lno_r*2  ,lno_c)
               endif
            end if
         end do
      end do
                                                  __TIMER_DO_STOP(1492)

      call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
       if (ierr /= 0) then
          write(nfout,*)' trans_scalapack_r :  mpi_waitall error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 10003, ierr)
       endif

      call mpi_waitall(icnt_send, req_s, sta_s, ierr)
       if (ierr /= 0) then
          write(nfout,*)' trans_scalapack_r :  mpi_waitall error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 10004, ierr)
       endif
                                                  __TIMER_COMM_STOP(1494)

      wk_rk = 0
      rank_c = myrank_e
      rank_r = myrank_g
                                                  __TIMER_DO_START(1493)
      do k = nis_e(rank_c), nie_e(rank_c)
         ib = neg_g_all(k)
         lrk_c = mod((ib-1)/scl_nb,scl_col_md)
         do j = nis_G_g1k(rank_r,ik), nie_G_g1k(rank_r,ik)
            lrk_r = mod((j-1)/scl_nb,scl_row_md)
            lrk = scl_row_md*lrk_c+lrk_r
               do lrank = 1, scl_md2_comm_rank_r
                  if (scl_md2_comm_rno_r(lrank) == lrk) then
                     wk_rk(lrank) = wk_rk(lrank) + 1
                     if (kimg == 1) then
                        amat_l(j,k-nis_e(rank_c)+1,1) = recv_buf(wk_rk(lrank),lrank)
                     else
                        amat_l(j,k-nis_e(rank_c)+1,1) = recv_buf(wk_rk(lrank)*2-1,lrank)
                        amat_l(j,k-nis_e(rank_c)+1,2) = recv_buf(wk_rk(lrank)*2  ,lrank)
                     endif
                     exit
                  end if
               end do
         end do
      end do

      call mpi_allreduce(amat_l,z_mat,nmatsz*np_e*kimg,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
                                                  __TIMER_DO_STOP(1493)

      deallocate(amat_l, stat=ierr)
      deallocate(req_r, stat=ierr)
      deallocate(req_s, stat=ierr)
      deallocate(sta_r, stat=ierr)
      deallocate(sta_s, stat=ierr)
      deallocate(send_buf, stat=ierr)
      deallocate(recv_buf, stat=ierr)
      deallocate(wk_rk, stat=ierr)
                                                  __TIMER_SUB_STOP(1475)
    end subroutine trans_scalapack_r

    subroutine eigsend(evec_w,n)
    integer :: i, n, ierr, itag=100
    integer, dimension(MPI_STATUS_SIZE) :: stat
    real(kind=DP) :: evec_w(n)

    if ((myrank_e==0) .and. (myrank_g==0)) then
       do i = nrank_e*nrank_g-1, nsclcol_md*nsclrow_md, -1
          call mpi_send(evec_w, n, mpi_double_precision, i, &
         &              itag, mpi_k_world(myrank_k), ierr)
           if (ierr /= 0) then
              write(nfout,*)' eigsend :  mpi_send error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world, 10010, ierr)
           end if
       end do
    end if

    if ((nsclcol_md*nsclrow_md-1) < (nrank_e*myrank_g+myrank_e)) then
       call mpi_recv(evec_w, n, mpi_double_precision, 0, &
      &              itag, mpi_k_world(myrank_k), stat, ierr)
        if (ierr /= 0) then
           write(nfout,*)' eigsend :  mpi_irecv error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world, 10012, ierr)
        end if
    end if
    end subroutine eigsend

  end subroutine m_ESmat_solve_Hx_eq_eSx_3D
  subroutine check_LAPACKinputs(str,nfout,n,lda,vl,vu,il,iu,lwork,abstol,RMat,CMat)
    character(len=*),intent(in) :: str
    integer, intent(in) ::         nfout,n,lda,il,iu,lwork
    real(DP), intent(in) ::        vl,vu,abstol
    real(DP),        optional, intent(in), dimension(lda,n) :: RMat
    complex(CMPLDP), optional, intent(in), dimension(lda,n) :: CMat
    integer :: i
    character(len=32) :: strp
    strp = "<<"//trim(str)//">>"
    if(iprimatdiagon>=2) then
    write(nfout,'(" --- n   = ", i8,1x,a32)') n,strp
    write(nfout,'(" --- lda = ",i8)') lda
    write(nfout,'(" --- vl  = ",d20.8, " vu = ",d20.8)') vl, vu
    write(nfout,'(" --- il  = ",i8, " iu = ", i8)') il, iu
    write(nfout,'(" --- lwork = ",i8)') lwork
    write(nfout,'(" --- abstol = ", d20.8)') abstol
    endif
    if(present(RMat))  then
      if(iprimatdiagon>=2) write(nfout,'(" --- amat(diagonal) = ",8f8.4)') (RMat(i,i),i=1,n)
    endif
    if(present(CMat)) then
       if(iprimatdiagon>=2) then
         write(nfout,'(" --- Real(cmat(diagonal)) = ",8f8.4)') (real(CMat(i,i)),i=1,n)
         write(nfout,'(" --- Imag(cmat(diagonal)) = ",8f8.4)') (aimag(CMat(i,i)),i=1,n)
       endif
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
       call phase_error_with_msg(nfout,' ne_found is not ne, stopped at <<check_LAPACKerror>>',__LINE__,__FILE__)
    end if
    if (info /= 0) then
       write(*,*) '### ERROR ### info /= 0'
       write(*,*) '   info ...',info
       call phase_error_with_msg(nfout,' info is not zero, stopped at <<check_LAPACKerror>>',__LINE__,__FILE__)
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
