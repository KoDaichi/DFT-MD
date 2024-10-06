!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 606 $)
!
!  MODULE: m_Force
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!    Further modification by T. Yamasaki   Apr/10/2007
!    Further modification by T. Yamasaki   Aug/31/2007
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
#ifdef VPP
#define _VECTOR_TUNING_
#endif
#ifdef SX
#define _VECTOR_TUNING_
#endif
#ifdef HIUX
#define _VECTOR_TUNING_
#endif

#ifndef NO_FORCE_DGEMM
#define FORCE_DGEMM
#endif

module m_Force
! $Id: m_Force.F90 606 2020-04-15 06:45:49Z ktagami $
  use m_IterationNumbers,     only : iteration_ionic
  use m_Charge_Density,       only : chgq_l, hsr
  use m_Electronic_Structure, only : zaj_l,vlhxcQ,occup_l,eko_l &
       &                           , fsr_l,fsi_l,vlhxc_l
  use m_Realspace,            only : nmesh_rs_max
  use m_NonLocal_Potential,   only : snl,snld_rs
  use m_XC_Potential,         only : vxc_l, vxcpc_l
  use m_XC_Potential,         only : m_XC_cal_potential
  use m_Kpoints,              only : k_symmetry
  use m_PlaneWaveBasisSet,    only : ngabc, nbmx, nbase, iba, kgp, kg1, ylm_l &
       &                           , m_pwBS_sphrp2,igfp_l
  use m_PseudoPotential,      only : nlmta,q,ilmt,lmtt,lmta,ltp,mtp,isph,dl2p,modnrm&
       &                           , dion,taup,il2p,iqitg,qitg_l, psc_l, rhpcg_l, itpcc &
       &                           , ival, nqitg, nlmtt &
       &                           , m_PP_find_maximum_l,m_PP_include_vanderbilt_pot &
       &                           , m_PP_set_index_arrays1,m_PP_set_index_arrays2 &
       &                           , ipaw, dion_paw
  use m_Crystal_Structure,    only : rltv, nopr, op, univol, altv, sw_supercell_symmetry
  use m_Ionic_System,         only : natm,natm2,iwei,imdtyp,ityp,ntyp,nfcatm,cps &
       &                           , ntyp_vdw,fxyzvdw_l,imdtypxyz,iatomn &
       &                           , nopr_supercell, napt_supercell, iop_supercell
  use m_Ionic_System,         only : num_regions, regions
  use m_Ionic_System,         only : sw_fix_bond, bond_forc
  use m_Timing,               only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters,   only : nspin, ipriforce, forccr, kimg, neg, af, nel_Ylm &
       &                           , f_tolerable_norm_error, f_tolerable_angle_error, f_tolerable_hyper_angle_error &
       &                           , sw_dipole_correction, sw_screening_correction, elec_field, sw_fef, sw_rspace &
       &                           , sw_vdw_correction, vdw_method &
#ifdef ENABLE_ESM_PACK
       &                           , m_CtrlP_cachesize, sw_hybrid_functional,sw_esm
#else
       &                           , m_CtrlP_cachesize, sw_hybrid_functional
#endif
  use m_Const_Parameters,     only : DP, zi, EXECUT, SKIP, PAI, PAI2, PAI4, DELTA, ON&
       &                           , Partial_Core_Charge, SIMPLE_CUBIC, GDIIS&
       &                           , RELAX, BONDLENGTH_FIX, BONDLENGTH_FIX_1 &
       &                           , BONDLENGTH_FIX_2, COG_FIX, HEAT_BATH, VXC_AND_EXC &
       &                           , WITHOUTTAG, WITHTAG, GAMMA &
       &                           , VDW_DFTD3
  use m_Parallelization,      only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype,ierr &
       &                           , map_k,myrank_k,np_e, ista_fs, iend_fs, ista_k, iend_k &
       &                           , ista_spin, iend_spin
  use m_Dipole,               only : dipole,field
  use m_Screening,            only : screening
#ifdef FORCE_DGEMM
! <-- T.Kokubo & D.Fukata, Feb. 2010
  use m_Control_Parameters,   only : nblocksize_force_is_given, nblocksize_force, nb_force_default
! -->
#endif

! =============================== added by K. Tagami ======================== 11.0
  use m_Control_Parameters,     only : noncol, ndim_spinor, ndim_magmom
  use m_PseudoPotential,        only : dion_scr_noncl, q_noncl
  use m_Const_Parameters,       only : CMPLDP
! ============================================================================ 11.0

  use m_FFT,                  only : fft_box_size_CD,nfftps
  use m_ES_nonlocal,          only : betar_dot_Psi_4_each_k_in_rs0

! === Positron ==== 2015/11/28
  use m_Control_Parameters,  only : sw_positron, positron_method
  use m_Const_Parameters,    only : Positron_GGGC, OFF
  use  m_Positron_Wave_Functions,  only : pchg_l
! ================= 2015/11/28

  use m_Control_Parameters,  only : sw_opencore
  use m_PS_opencore,  only : has_opencore, rmag_opencore_l, mag_opencore_pol

  use m_Files, only : nfout
  use mpi

  implicit none
!  include 'mpif.h'

  real(kind=DP), allocatable, dimension(:,:)       :: forc_l
  real(kind=DP)                                    :: forcmx
  real(kind=DP), private,allocatable, dimension(:,:)   :: fnlxyz_l, flhxcq_l
  real(kind=DP), private,allocatable, dimension(:,:)   :: force_previous
  integer, private ::                                 force_error_check_counter = 0
!  real(kind=DP), private,pointer, dimension(:,:)   :: fdip_l,fext_l
!  real(kind=DP), private,pointer, dimension(:,:)   :: ffef_l
  real(kind=DP), private,allocatable, dimension(:,:)   :: fdip_l,fext_l
  real(kind=DP), private,allocatable, dimension(:,:)   :: ffef_l
  real(kind=DP), allocatable, dimension(:,:)   :: fexx_l

!!$  real(kind=DP), private,allocatable,target,dimension(:):: zfcos, zfsin
  integer,       private,pointer, dimension(:)     :: il3
  real(kind=DP), private,pointer, dimension(:,:)   :: work
  real(kind=DP), private,pointer, dimension(:)     :: ylm
  real(kind=DP), private,pointer, dimension(:)     :: ar, ai
  real(kind=DP), private,pointer, dimension(:,:,:) :: gr,gi ! derivative of (fsr,fsi)
  real(kind=DP), private,pointer, dimension(:,:)   :: forc_up, forc_down, forcw
  integer,       private,pointer, dimension(:)     :: ipfrc
  integer,       private  :: iteration_ionic_wd_forcmx = 0
contains
  subroutine m_Force_alloc
    if(allocated(forc_l)) deallocate(forc_l)
    allocate(forc_l(natm,3));forc_l=0.d0
    allocate(fnlxyz_l(natm,3))
    allocate(flhxcq_l(natm,3))
    if(sw_dipole_correction==ON) then
      allocate(fdip_l(natm,3))
      allocate(fext_l(natm,3))
    end if
    if(sw_fef == ON) allocate(ffef_l(natm,3))
    if(sw_hybrid_functional == ON) allocate(fexx_l(natm,3))
  end subroutine m_Force_alloc

  subroutine m_Force_sumup_and_symmetrize(nfout,fxyzew_l,napt)
    integer, intent(in) :: nfout
    real(kind=DP), intent(in), dimension(natm,3) :: fxyzew_l
    integer,intent(in), dimension(natm,nopr+af)     :: napt

    integer :: ia, j, k, iat
    integer,save :: id_sname = -1

    call tstatc0_begin('m_Force_sumup_and_symmetrize ',id_sname)
#ifdef __TIMER_SUB__
    call timer_sta(1506)
#endif
    if(ipriforce >= 2) then
       write(nfout,'(" force elements ")')
       write(nfout,'(" -- forc_l (local+pcc) --")')
       do ia = 1, natm
          write(nfout,'(i5,3f16.8)') ia,(forc_l(ia,j),j=1,3)
       end do
       write(nfout,'(" -- fnlxyz_l --")')
       do ia = 1, natm
          write(nfout,'(i5,3f16.8)') ia,(fnlxyz_l(ia,j),j=1,3)
       end do
       write(nfout,'(" -- fxyzew_l --")')
       do ia = 1, natm
          write(nfout,'(i5,3f16.8)') ia,(fxyzew_l(ia,j),j=1,3)
       end do
       write(nfout,'(" -- flhxcq_l --")')
       do ia = 1, natm
          write(nfout,'(i5,3f16.8)') ia,(flhxcq_l(ia,j),j=1,3)
       end do
       if(sw_dipole_correction==ON) then
          write(nfout,'(" -- fdip_l --")')
          do ia = 1, natm
             write(nfout,'(i5,3f16.8)') ia,(fdip_l(ia,j),j=1,3)
          end do
          write(nfout,'(" -- fext_l --")')
          do ia = 1, natm
             write(nfout,'(i5,3f16.8)') ia,(fext_l(ia,j),j=1,3)
          end do
       end if
       if(sw_fef == ON) then
          write(nfout,'(" -- ffef_l --")')
          do ia = 1, natm
             write(nfout,'(i5,3f16.8)') ia,(ffef_l(ia,j),j=1,3)
          end do
       end if

       if ( sw_vdw_correction == ON ) then
          if ( ntyp_vdw > 0 .or. vdw_method == VDW_DFTD3 ) then
             write(nfout,'(" -- fxyzvdw_l --")')
             do ia = 1, natm
                write(nfout,'(i5,3f16.8)') ia,(fxyzvdw_l(ia,j),j=1,3)
             end do
          endif
       endif

       if(num_regions>0) then
          write(nfout,'(" -- force from regions --")')
          do k=1,num_regions
          write(nfout,'(a,i5)') ' from region ',k
          do iat = 1, regions(k)%ntarget_atoms
             ia = regions(k)%target_atoms(iat)
             write(nfout,'(i5,3f16.8)') ia,(regions(k)%forc(ia,j),j=1,3)
          end do
          enddo
       endif
       if(sw_hybrid_functional == ON) then
          write(nfout,'(" -- fexx_l --")')
          do ia = 1, natm
             write(nfout,'(i5,3f16.8)') ia,(fexx_l(ia,j),j=1,3)
          end do
       end if
    end if

    forc_l = forc_l + fnlxyz_l + fxyzew_l + flhxcq_l
    if(sw_dipole_correction==ON) then
       forc_l = forc_l + fdip_l + fext_l
    end if
    if(sw_fef==ON) forc_l = forc_l + ffef_l
    if(sw_hybrid_functional==ON) forc_l = forc_l + fexx_l

    if ( sw_vdw_correction == ON ) then
       if ( ntyp_vdw > 0 .or. vdw_method == VDW_DFTD3 ) then
          forc_l = forc_l + fxyzvdw_l
       endif
    endif

    if(num_regions>0) then
        do k=1,num_regions
           forc_l = forc_l + regions(k)%forc
        enddo
    endif

    if(ipriforce >= 2) then
       write(nfout,'(" -- total forc_l --")')
       do ia = 1, natm
          write(nfout,'(i5,3f16.8)') ia,(forc_l(ia,j),j=1,3)
          call flush(nfout)
       end do
    end if
!!$    if(nbztyp >= SIMPLE_CUBIC) then
    if(nopr >= 2) then
       allocate(forcw(natm,3)); allocate(ipfrc(natm2))
!!$       write(nfout,'(" sw_supercell_symmetry = ",i8)') sw_supercell_symmetry
       call flush(nfout)
       if(nopr_supercell <= nopr .or. sw_supercell_symmetry /= ON) then
          call fd_symmetrize(natm2,natm,natm,napt,nopr+af,nopr,op,iwei &
               &, forc_l, forcw, ipfrc) ! -(bottom_Subroutines_para)
!!$          call fd_symmetrize(natm2,natm,natm,napt,nopr+af,op,nopr+af,nopr &
!!$               & ,iwei, forc_l, forcw, ipfrc) ! -(bottom_Subroutines_para)
       else
!!$          call fd_symmetrize(natm2,natm,natm,napt_supercell,nopr_supercell,op,nopr+af,nopr_supercell &
!!$               & ,iwei, forc_l, forcw, ipfrc, iop_supercell)
          call fd_supercell_symmetrize(natm2,natm,natm,napt_supercell,nopr_supercell,op,nopr+af,iwei &
               & ,forc_l, forcw, ipfrc, iop_supercell)
       end if
       deallocate(forcw); deallocate(ipfrc)
       if(ipriforce >= 2) then
          write(nfout,'(" -- total forc_l symmetrized --")')
          do ia = 1, natm
             write(nfout,'(i5,3f16.8)') ia,(forc_l(ia,j),j=1,3)
          end do
       end if
    end if
#ifdef __TIMER_SUB__
    call timer_end(1506)
#endif

  end subroutine m_Force_sumup_and_symmetrize

  subroutine m_Force_term_Elocal_and_Epc(nfout,pos)
    integer, intent(in) :: nfout
    real(kind=DP), intent(in), dimension(natm,3) :: pos
    real(kind=DP), allocatable,dimension(:):: zfcos, zfsin
    integer :: ia, it, ipc
    integer :: id_sname = -1
    complex(kind=CMPLDP),allocatable,dimension(:) :: aux
    integer :: ig,ikm,kimgtmp
    real(kind=DP), allocatable, dimension(:,:) :: forc_l_esm
    real(kind=DP), allocatable, dimension(:,:) :: agauss,bgauss
    real(kind=DP), dimension(2) :: alp,cc
    real(kind=DP) :: ivaltmp
    integer :: i,j,itpcctmp,nfftcd
    call tstatc0_begin('m_Force_term_Elocal_and_Epc ', id_sname,1)
    allocate(forcw(natm,3))
    allocate(zfsin(ista_kngp:iend_kngp))
    allocate(zfcos(ista_kngp:iend_kngp))
#ifndef REMOVE_PC_FROM_FORCE
    call m_XC_cal_potential(nfout,Partial_Core_Charge,chgq_l,VXC_AND_EXC)
                                               !  -> vxcpc_l
#endif
    forcw = 0.d0
    do ia = 1, natm
       it = ityp(ia); ipc = itpcc(it)
       call calc_phase2(natm,pos,ia,kgp,ngabc,ista_kngp,iend_kngp,zfcos,zfsin)
       call force_from_Elocal_4_each_atom ! -(contained here) ->forcw
       if(ipc /= 0 ) call force_from_Epc_for_each_atom    ! -> forcw
    end do

    do ia = 1, natm
       call pucv2cart(rltv,forcw(ia,1),forcw(ia,2),forcw(ia,3))
       forc_l(ia,1:3) = forc_l(ia,1:3) + forcw(ia,1:3)
    end do

    deallocate(zfcos,zfsin)
    deallocate(forcw)

#ifdef ENABLE_ESM_PACK
    if(sw_esm==ON)then
       nfftcd = fft_box_size_CD(1,0)*fft_box_size_CD(2,0)*fft_box_size_CD(3,0)
       allocate(aux(nfftcd));aux(:) = (0.d0,0.d0)
       allocate(agauss(natm,2))
       allocate(bgauss(natm,2))
       do i=1,natm
          call psbhs0(nfout,nint(iatomn(ityp(i))),ivaltmp,itpcctmp,alp,cc)  ! -(b_PseudoPotential)
          agauss(i,1) = alp(1)
          agauss(i,2) = alp(2)
          bgauss(i,1) = cc(1)
          bgauss(i,2) = cc(2)
       enddo
       do ig=ista_kngp,iend_kngp
         if(kimg==1)then
           if(nspin==1)then
             aux(igfp_l(ig)) = dcmplx(chgq_l(ig,1,1),0.0d0)
           else
             aux(igfp_l(ig)) = dcmplx(chgq_l(ig,1,1)+chgq_l(ig,1,2),0.0d0)
           endif
         else
           if(nspin==1 .or. noncol )then
             aux(igfp_l(ig)) = dcmplx(chgq_l(ig,1,1),chgq_l(ig,2,1))
           else
             aux(igfp_l(ig)) = dcmplx(chgq_l(ig,1,1)+chgq_l(ig,1,2),chgq_l(ig,2,1)+chgq_l(ig,2,2))
           endif
         endif
       enddo
       allocate(forc_l_esm(3,natm));forc_l_esm=0.d0
       call mpi_allreduce(MPI_IN_PLACE,aux,nfftps/kimg,mpi_double_complex,mpi_sum,MPI_CommGroup,ierr)
       call esm_force_lc_(nfftcd,aux,natm,forc_l_esm,2,agauss,bgauss)
       do i=1,natm
          do j=1,3
             forc_l(i,j) = forc_l(i,j) + forc_l_esm(j,i)*0.5d0
          enddo
       enddo
       deallocate(aux)
       deallocate(forc_l_esm)
       deallocate(agauss,bgauss)
    endif
#endif

    call tstatc0_end(id_sname)
  contains
    subroutine force_from_Epc_for_each_atom
      integer :: is, i
      real(kind=DP) :: tp, weight
      real(kind=DP), dimension(3) :: f1, f2
      real(kind=DP), dimension(3) :: f1_mpi, f2_mpi ! MPI

! ============================ added by K. Tagami ===================== 11.0
      integer :: is_max
! ===================================================================== 11.0

      f1 = 0.d0
      f2 = 0.d0
      f1_mpi = 0.d0
      f2_mpi = 0.d0

#ifdef NEC_TUNE_SMP
!CDIR NOCONCUR
#endif

! =============================== added by K. Tagami ==================== 11.0
      if ( noncol ) then
         is_max = 1
      else
         is_max = nspin
      endif
! ======================================================================== 11.0

! ================================== modified by K. Tagami ============= 11.0
!      do is = 1, nspin
      do is = 1, is_max
! ======================================================================== 11.0

         if(kimg == 1) then
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
            do i = ista_kngp, iend_kngp  !for mpi
               tp  = zfsin(i)*vxc_l(i,1,is)*rhpcg_l(i,ipc)
               f1_mpi(1) = f1_mpi(1) + tp*ngabc(i,1)
               f1_mpi(2) = f1_mpi(2) + tp*ngabc(i,2)
               f1_mpi(3) = f1_mpi(3) + tp*ngabc(i,3)
            end do
         else if(kimg == 2) then
!!$#ifdef NEC_TUNE_SMP
!!$!CDIR INNER
!!$#endif
            do i = ista_kngp, iend_kngp  !for mpi
               tp  = (zfsin(i)*vxc_l(i,1,is)+zfcos(i)*vxc_l(i,2,is))&
                    &*rhpcg_l(i,ipc)
               f1_mpi(1) = f1_mpi(1) + tp*ngabc(i,1)
               f1_mpi(2) = f1_mpi(2) + tp*ngabc(i,2)
               f1_mpi(3) = f1_mpi(3) + tp*ngabc(i,3)
            end do
         end if
      end do

      if ( sw_opencore == ON .and. nspin == 2 ) then
         if ( has_opencore(it) == 1 ) then
            if ( noncol ) then
            else
               do is = 1, is_max
                  if ( is == 1 ) then
                     weight = 0.5d0 *mag_opencore_pol(ia,1)
                  else
                     weight = -0.5d0 *mag_opencore_pol(ia,1)
                  endif

                  if ( kimg == 1 ) then
                     do i = ista_kngp, iend_kngp  !for mpi
                        tp  = zfsin(i)*vxc_l(i,1,is) *rmag_opencore_l(i,it) *weight
                        f1_mpi(1) = f1_mpi(1) + tp*ngabc(i,1)
                        f1_mpi(2) = f1_mpi(2) + tp*ngabc(i,2)
                        f1_mpi(3) = f1_mpi(3) + tp*ngabc(i,3)
                     end do
                  else if ( kimg == 2 ) then
                     do i = ista_kngp, iend_kngp  !for mpi
                        tp  = (zfsin(i)*vxc_l(i,1,is)+zfcos(i)*vxc_l(i,2,is))&
                             &   *rmag_opencore_l(i,it) *weight
                        f1_mpi(1) = f1_mpi(1) + tp*ngabc(i,1)
                        f1_mpi(2) = f1_mpi(2) + tp*ngabc(i,2)
                        f1_mpi(3) = f1_mpi(3) + tp*ngabc(i,3)
                     end do
                  endif
               end do
            endif
         end if
      endif

      call mpi_allreduce(f1_mpi,f1,3 &
           &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)

! ================================ modified by K. Tagami =================== 11.0
!!      if(nspin == 2) f1 = f1/2.d0
!
      if ( noncol ) then
      else
         if(nspin == 2) f1 = f1/2.d0
      endif
! =========================================================================== 11.0

#ifndef REMOVE_PC_FROM_FORCE
      if(kimg == 1) then
         do i = ista_kngp, iend_kngp  !for mpi
            tp  = zfsin(i)*vxcpc_l(i,1)*rhpcg_l(i,ipc)
            f2_mpi(1) = f2_mpi(1) + tp*ngabc(i,1)
            f2_mpi(2) = f2_mpi(2) + tp*ngabc(i,2)
            f2_mpi(3) = f2_mpi(3) + tp*ngabc(i,3)
         end do
      else if(kimg == 2) then
         do i = ista_kngp, iend_kngp  !for mpi
            tp  = (zfsin(i)*vxcpc_l(i,1)+zfcos(i)*vxcpc_l(i,2))*rhpcg_l(i,ipc)
            f2_mpi(1) = f2_mpi(1) + tp*ngabc(i,1)
            f2_mpi(2) = f2_mpi(2) + tp*ngabc(i,2)
            f2_mpi(3) = f2_mpi(3) + tp*ngabc(i,3)
         end do
!xocl end spread
      end if
      call mpi_allreduce(f2_mpi,f2,3 &
           &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
!!$      forcw(ia,1) = forcw(ia,1) + (f1(1)/nspin - f2(1))*univol
!!$      forcw(ia,2) = forcw(ia,2) + (f1(2)/nspin - f2(2))*univol
!!$      forcw(ia,3) = forcw(ia,3) + (f1(3)/nspin - f2(3))*univol
      forcw(ia,1) = forcw(ia,1) + (f1(1) - f2(1))*univol
      forcw(ia,2) = forcw(ia,2) + (f1(2) - f2(2))*univol
      forcw(ia,3) = forcw(ia,3) + (f1(3) - f2(3))*univol
#else
!!$      forcw(ia,1) = forcw(ia,1) + (f1(1)/nspin)*univol
!!$      forcw(ia,2) = forcw(ia,2) + (f1(2)/nspin)*univol
!!$      forcw(ia,3) = forcw(ia,3) + (f1(3)/nspin)*univol
      forcw(ia,1) = forcw(ia,1) + f1(1)*univol
      forcw(ia,2) = forcw(ia,2) + f1(2)*univol
      forcw(ia,3) = forcw(ia,3) + f1(3)*univol
#endif

    end subroutine force_from_Epc_for_each_atom

    subroutine force_from_Elocal_4_each_atom
      integer :: is, i
      real(kind=DP) :: tmp
      real(kind=DP) :: fx,fy,fz
      real(kind=DP), pointer, dimension(:) :: f1_mpi,f2_mpi ! MPI

! ============================ added by K. Tagami ===================== 11.0
      integer :: is_max
! ===================================================================== 11.0

      fx = 0.d0; fy = 0.d0; fz = 0.d0
      allocate(f1_mpi(3)); allocate(f2_mpi(3))              ! MPI

#ifdef NEC_TUNE_SMP
!CDIR NOCONCUR
#endif

! =============================== added by K. Tagami ==================== 11.0
      if ( noncol ) then
         is_max = 1
      else
         is_max = nspin
      endif
! ======================================================================== 11.0

! ================================== modified by K. Tagami ============= 11.0
!      do is = 1, nspin
      do is = 1, is_max
! ======================================================================== 11.0

         if(kimg == 1) then
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
            do i = ista_kngp, iend_kngp  !for mpi
               tmp = zfsin(i)*chgq_l(i,1,is)*psc_l(i,it)
               if(sw_screening_correction==ON) then
                 tmp = zfsin(i)*chgq_l(i,1,is)*(psc_l(i,it)-ival(it)*screening%phik(i)/univol)
               end if
               fx = fx + tmp*ngabc(i,1)
               fy = fy + tmp*ngabc(i,2)
               fz = fz + tmp*ngabc(i,3)
            end do
         else if(kimg == 2) then
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
            do i = ista_kngp, iend_kngp  !for mpi
               tmp = (zfsin(i)*chgq_l(i,1,is)+zfcos(i)*chgq_l(i,2,is))&
                    &*psc_l(i,it)
               if(sw_screening_correction==ON) then
                 tmp = (zfsin(i)*chgq_l(i,1,is)+zfcos(i)*chgq_l(i,2,is))&
                      &*(psc_l(i,it)-ival(it)*Screening%phik(i)/univol)
               end if
               fx = fx + tmp*ngabc(i,1)
               fy = fy + tmp*ngabc(i,2)
               fz = fz + tmp*ngabc(i,3)
            end do
         end if
      end do

! == POSITRON SCF === 2015/11/28
      if ( sw_positron /= OFF ) then
         if ( positron_method == Positron_GGGC ) then
            if(kimg == 1) then
               do i = ista_kngp, iend_kngp  !for mpi
                  tmp = -zfsin(i) *pchg_l(i,1) *psc_l(i,it)
                  fx = fx + tmp*ngabc(i,1)
                  fy = fy + tmp*ngabc(i,2)
                  fz = fz + tmp*ngabc(i,3)
               end do
            else
               do i = ista_kngp, iend_kngp  !for mpi
                  tmp = -( zfsin(i) *pchg_l(i,1) +zfcos(i)*pchg_l(i,2) )&
                       &  *psc_l(i,it)
                  fx = fx + tmp*ngabc(i,1)
                  fy = fy + tmp*ngabc(i,2)
                  fz = fz + tmp*ngabc(i,3)
               end do
            endif
         endif
      endif
! ============= 2015/11/28

      f1_mpi(1) = fx; f1_mpi(2) = fy; f1_mpi(3) = fz
      call mpi_allreduce(f1_mpi,f2_mpi,3 &
           &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
      forcw(ia,1:3) = f2_mpi(1:3)*univol
      deallocate(f2_mpi); deallocate(f1_mpi)
    end subroutine force_from_Elocal_4_each_atom
  end subroutine m_Force_term_Elocal_and_Epc

  subroutine m_Force_initialize
#ifdef __TIMER_SUB__
    call timer_sta(1501)
#endif
    forc_l = 0.d0
#ifdef __TIMER_SUB__
    call timer_end(1501)
#endif
  end subroutine m_Force_initialize

  subroutine m_Force_term_drv_of_vlhxcQ(nfout,pos,napt)
    !  ----
    ! (Rev) T. Yamaskai, 31, Aug, 2007
    !     1. 'call set_index_arrays1' that included a bug is replaced
    !       by 'call m_PP_set_index_arrays1', whose bug is fixed.
    !     2. 'call set_index_arrays2' is also replaced by 'call
    !       m_PP_set_index_arrays2' that can be referred from other modules.
    !     3. contained subroutines, set_index_arrays1 and set_index_arrays2 were
    !       deleted.
    integer,intent(in)                           :: nfout
    real(kind=DP), intent(in), dimension(natm,3) :: pos
    integer,intent(in), dimension(natm,nopr+af)  :: napt

    real(kind=DP), allocatable, dimension(:)                  :: ylm_t
    real(kind=DP), allocatable, target, dimension(:,:)        :: ylm_ext
    real(kind=DP), allocatable, dimension(:):: zfcos, zfsin

#ifndef _VECTOR_TUNING_
    real(kind=DP), allocatable, dimension(:,:) :: flhxcq_mpi
    real(kind=DP), allocatable, dimension(:) :: ylm_sum
    integer :: ncache, ibsize, iwidth, ibl1, ibl2, iwidthmax
#endif

#ifndef _m_Force_no_loop_exchange_
    integer :: m, maxm, ip, np, iq
    integer, parameter :: mcritical = 4*2+1
    integer, allocatable, dimension(:) :: nqitg_sp, nqitg_sp0 !d(ntyp)
    integer, allocatable, dimension(:) :: iq2l3 ! d(nqitg)
    integer, allocatable, dimension(:,:) :: nc  ! d(maxm,nqitg)
    integer :: mc ! maxval(nc)
    integer, allocatable, dimension(:,:,:) :: nc2lmt1, nc2lmt2, nc2n ! d(mc,maxm,nqitg)
#endif

! ============================ added by K. Tagami ===================== 11.0
    integer :: is_max
! ===================================================================== 11.0

    integer :: is,it,lmt1,lmt2,n,ia,mdvdb,il1,tau1,il2,tau2,ilm3,l3,iiqitg
    real(kind=DP) :: fx,fy,fz,shdg_x
    integer       :: id_sname = -1
    call tstatc0_begin('m_Force_term_drv_of_vlhxcQ ',id_sname,1)

    call m_PP_find_maximum_l(n)
    n = (n-1) + (n-1) + 1

    if(modnrm == EXECUT) then
       allocate(il3(n**2)); call substitute_il3(n**2,il3) ! -(b_Elec..)
       if(n**2 > nel_Ylm) then
          allocate(ylm_ext(ista_kngp:iend_kngp,nel_Ylm+1:n**2)); ylm_ext = 0.d0
       end if
       allocate(ylm_t(ista_kngp:iend_kngp)); ylm_t = 0.d0
       do ilm3 = nel_Ylm+1, n**2
          call m_pwBS_sphrp2(ilm3,rltv,ista_kngp,iend_kngp,ylm_t)
          ylm_ext(:,ilm3) = ylm_t(:)
       end do
       deallocate(ylm_t)

#ifndef _m_Force_no_loop_exchange_
       allocate(nqitg_sp(ntyp)); allocate(nqitg_sp0(ntyp))
       allocate(iq2l3(nqitg))
       allocate(nc(mcritical,nqitg));nc=0
       call m_PP_set_index_arrays1(nfout,ntyp,nqitg,mcritical,n**2,il3 &
            & ,maxm,mc,nqitg_sp,nqitg_sp0,iq2l3,nc)  ! -> nqitg_sp, nqitg_sp0, iq2l3,nc, maxm, mc
       allocate(nc2lmt1(mc,maxm,nqitg))
       allocate(nc2lmt2(mc,maxm,nqitg))
       allocate(nc2n(mc,maxm,nqitg))
       call m_PP_set_index_arrays2(nfout,mc,maxm,nqitg,mcritical,n**2,il3,iq2l3 &
            & ,nc2lmt1,nc2lmt2,nc2n,nc) ! -> nc2lmt1, nc2lmt2, nc2n, nc
#endif

    end if

! =============================== added by K. Tagami ==================== 11.0
    if ( noncol ) then
       is_max = ndim_magmom
    else
       is_max = nspin
    endif
! ======================================================================== 11.0

    flhxcq_l = 0.d0
#ifdef _VECTOR_TUNING_
    allocate(zfcos(ista_kngp:iend_kngp))
    allocate(zfsin(ista_kngp:iend_kngp))
    do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(mdvdb == SKIP) cycle

       do ia= 1,natm
          if(ityp(ia) /= it) cycle
          fx = 0.d0; fy = 0.d0; fz = 0.d0
          call calc_phase2(natm,pos,ia,kgp,ngabc,ista_kngp,iend_kngp,zfcos,zfsin)

! ========================================= modified by K. Tagami ========= 11.0
!          do is = 1, nspin, af+1
          do is = 1, is_max, af+1
! ========================================================================= 11.0

#ifndef _m_Force_no_loop_exchange_
             do iq = nqitg_sp0(it), nqitg_sp(it)
                l3 = iq2l3(iq)
                do m = 1, 2*l3+1
                   ilm3 = l3*l3+m
                   call sum_hsr_dot_gauntc(is,it,ia,iq,m)
                   if(ilm3 <= nel_Ylm) then
                      ylm => ylm_l(ista_kngp:iend_kngp,ilm3)
                   else
                      ylm => ylm_ext(ista_kngp:iend_kngp,ilm3)
                   end if
                   call add_hardpart_to_flhxcq_l_core(iq,ista_kngp,iend_kngp)
                   !    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   !   --> fx,fy,fz
                end do
             end do
#else
             do lmt1 = 1, ilmt(it)
                il1 = ltp(lmt1,it); tau1 = taup(lmt1,it)
                do lmt2 = lmt1, ilmt(it)
                   il2 = ltp(lmt2,it); tau2 = taup(lmt2,it)
                   do n = 1, il2p(lmt1,lmt2,it)
                      ilm3 = isph(lmt1,lmt2,n,it);    l3   =  il3(ilm3)
                      iiqitg = iqitg(il1,tau1,il2,tau2,l3+1,it)
                      if(iiqitg == 0) cycle
                      call hsr_dot_gaunt(is,it,ia,lmt1,lmt2,n)
                      if(ilm3 <= nel_Ylm) then
                         ylm => ylm_l(ista_kngp:iend_kngp,ilm3)
                      else
                         ylm => ylm_ext(ista_kngp:iend_kngp,ilm3)
                      end if
                      call add_hardpart_to_flhxcq_l_core(iiqitg,ista_kngp,iend_kngp)
                      !    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      !   --> fx,fy,fz
                   end do ! n
                end do ! lmt2
             end do
#endif
          end do
          call pucv2cart(rltv,fx,fy,fz)
!!$       call fxyz_gsum(fx,fy,fz)
          flhxcq_l(ia,1) = fx; flhxcq_l(ia,2) = fy; flhxcq_l(ia,3) = fz
       end do
    end do
    deallocate(zfsin, zfcos)
#else

#ifndef _m_Force_no_loop_exchange_
    ncache = (m_CtrlP_cachesize()*1024)*3/4
    if(ncache == 0) ibsize = iend_kngp - ista_kngp + 1

    do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(mdvdb == SKIP) cycle
       if(ipriforce >= 2) write(nfout,'(" !mForce: it = ",i8)') it

       if(ncache/=0) then
          iwidth = nqitg_sp(it) - nqitg_sp0(it) + 1
          if(kimg == 1) then ! qitg_l(i,iq)*ylm(iy)*zfsc1_1(iy)
             ibsize=ncache/(8*(iwidth + 1))
          else ! qitg_l(i,iq),ylm(iy),zfsc1_1(iy),zfsc2_1(iy)
             ibsize=ncache/(8*(iwidth + 2))
          endif
          if(ipriforce >= 2) write(nfout,'(" !mForce: ibsize, iwidth, ncache = ",3i8)') ibsize, iwidth, ncache
       end if

       IBLOCK: do ibl1=ista_kngp,iend_kngp,ibsize
          ibl2=ibl1+ibsize-1
          if(ibl2.gt.iend_kngp) ibl2=iend_kngp
          if(ibl2.gt.kgp) ibl2 = kgp
          allocate(ylm_sum(ibl1:ibl2))
          allocate(zfcos(ibl1:ibl2)); allocate(zfsin(ibl1:ibl2))

       do ia= 1,natm
          if(ityp(ia) /= it) cycle
          fx = 0.d0; fy = 0.d0; fz = 0.d0
          call calc_phase_div(ia)  ! -> zfsin, zfcos

! ========================================= modified by K. Tagami ========= 11.0
!          do is = 1, nspin, af+1
             do is = 1, is_max, af+1
! ========================================================================= 11.0

             do iq = nqitg_sp0(it), nqitg_sp(it)
                l3 = iq2l3(iq)
                ylm_sum = 0.d0
                do m = 1, 2*l3+1
                   ilm3 = l3*l3+m
                   call sum_hsr_dot_gauntc(is,it,ia,iq,m)
                   if(ilm3 <= nel_Ylm) then
                      ylm_sum(ibl1:ibl2) = ylm_sum(ibl1:ibl2) + shdg_x*ylm_l(ibl1:ibl2,ilm3)
!!!$                      ylm => ylm_l(:,ilm3)
                   else
                      ylm_sum(ibl1:ibl2) = ylm_sum(ibl1:ibl2) + shdg_x*ylm_ext(ibl1:ibl2,ilm3)
!!!$                      ylm => ylm_ext(:,ilm3)
                   end if
!!!$                   call add_hardpart_to_flhxcq_l_core(iq,ibl1,ibl2)
!!!$                   !    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!!$                   !   --> fx,fy,fz
                end do
                call add_hardpart_to_flhxcq_l_div(iq,ibl1,ibl2)
             end do
          end do
          call pucv2cart(rltv,fx,fy,fz)
          flhxcq_l(ia,1)=flhxcq_l(ia,1)+fx; flhxcq_l(ia,2)=flhxcq_l(ia,2)+fy; flhxcq_l(ia,3)=flhxcq_l(ia,3)+fz
       end do
       deallocate(zfsin,zfcos)
       deallocate(ylm_sum)
       end do IBLOCK
    end do
    if(npes>1) then
       allocate(flhxcq_mpi(natm,3)); flhxcq_mpi = 0.d0
       call mpi_allreduce(flhxcq_l,flhxcq_mpi,3*natm,mpi_double_precision &
            & , mpi_sum, MPI_CommGroup,ierr)
       flhxcq_l = flhxcq_mpi
       deallocate(flhxcq_mpi)
    end if
#else
    call phase_error_with_msg(nfout,' you cannot define _m_Force_no_loop_exchange_',__LINE__,__FILE__)
#endif
#endif

!!$!xocl spread do/ind_natm
!!$    do ia = 1, natm
!!$       fx = flhxcq_l(ia,1); fy = flhxcq_l(ia,2); fz = flhxcq_l(ia,3)
!!$       call pucv2cart(rltv,fx,fy,fz)
!!$       flhxcq_l(ia,1) = fx; flhxcq_l(ia,2) = fy; flhxcq_l(ia,3) = fz
!!$    end do
!!$!xocl end spread

    if(af /= 0) then
       allocate(forc_up(natm2,3));allocate(forc_down(natm2,3))
       call forcaf2(natm,natm2,nopr+af,natm,nopr,napt&
            &, iwei,op,ipriforce,af,flhxcq_l,forc_up,forc_down)
       deallocate(forc_up); deallocate(forc_down)
    endif
    if(modnrm == EXECUT) then
       if(allocated(ylm_ext)) deallocate(ylm_ext)
       deallocate(il3)
#ifndef _m_Force_no_loop_exchange_
       deallocate(nc2n,nc2lmt2,nc2lmt1,nc,iq2l3,nqitg_sp,nqitg_sp0)
#endif
    end if

    call tstatc0_end(id_sname)
  contains
#ifndef SX
    subroutine calc_phase_div(ia)
      integer, intent(in) :: ia
      integer :: ig
      real(kind=DP) :: fx,fy,fz, ph
      fx = pos(ia,1)*PAI2
      fy = pos(ia,2)*PAI2
      fz = pos(ia,3)*PAI2
      do ig = ibl1, ibl2
         ph = ngabc(ig,1)*fx + ngabc(ig,2)*fy + ngabc(ig,3)*fz
         zfcos(ig) = dcos(ph)
         zfsin(ig) = dsin(ph)
      end do
    end subroutine calc_phase_div
#endif

#ifndef _m_Force_no_loop_exchange_
    subroutine sum_hsr_dot_gauntc(is,it,ia,iq,m)
      integer, intent(in) :: is,it,ia,iq,m
      integer :: ip, lmt1, lmt2, np
      real(kind=DP) :: fac
      shdg_x = 0.d0
      do ip = 1, nc(m,iq)
         lmt1 = nc2lmt1(ip,m,iq)
         lmt2 = nc2lmt2(ip,m,iq)
         np = nc2n(ip,m,iq)
         fac = 2.d0; if(lmt1 == lmt2) fac = 1.d0
         shdg_x = shdg_x + &
              & fac*hsr(ia,lmt1,lmt2,is)*dl2p(lmt1,lmt2,np,it)
      end do
    end subroutine sum_hsr_dot_gauntc
#else
    subroutine hsr_dot_gaunt(is,it,ia,lmt1,lmt2,n)
      integer, intent(in) :: is, it, ia, lmt1,lmt2,n
      real(kind=DP) :: fac

      fac = 2.d0; if(lmt1 == lmt2) fac = 1.d0
      shdg_x = fac*hsr(ia,lmt1,lmt2,is)*dl2p(lmt1,lmt2,n,it)
    end subroutine hsr_dot_gaunt

#endif
    subroutine add_hardpart_to_flhxcq_l_core(iq,ibl1,ibl2)
      integer, intent(in) :: iq,ibl1,ibl2
      integer       :: i, iy
      real(kind=DP) :: flchgq, f, tx,ty,tz, f2
      real(kind=DP),pointer,dimension(:) :: f1_mpi, f2_mpi
!!$      real(kind=DP),pointer,dimension(:) :: zfsc1, zfsc2

!!$      if(mod(l3,2) == 0) then
!!$         flchgq = fac*real(zi**(-l3))*dga*hsr(ia,lmt1,lmt2,is)
!!$         flchgq = real(zi**(-l3))*shdg_x
!!$         if(kimg == 1) then
!!$            zfsc1 => zfsin
!!$         else
!!$            f2 = 1;  zfsc1 => zfsin; zfsc2 => zfcos
!!$         end if
!!$      else
!!$         flchgq = - aimag(zi**(-l3))*shdg_x
!!$         if(kimg == 1) then
!!$            zfsc1 => zfcos
!!$         else
!!$            f2 = -1; zfsc1 => zfcos; zfsc2 => zfsin
!!$         end if
!!$      end if

      tx = 0.d0; ty = 0.d0; tz = 0.d0
      if(mod(l3,2) == 0) then
         flchgq = real(zi**(-l3))*shdg_x
         if(kimg == 1) then
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
            do i = ibl1, ibl2  !for mpi
               iy = i-ista_kngp+1
!!$               f  = qitg_l(i,iq)*ylm(iy)*vlhxc_l(i,1,is)*zfsc1(iy)
               f  = qitg_l(i,iq)*ylm(iy)*vlhxc_l(i,1,is)*zfsin(i)
               tx = tx + f*ngabc(i,1);  ty = ty + f*ngabc(i,2)
               tz = tz + f*ngabc(i,3)
            end do
         else if(kimg == 2) then
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
            do i = ibl1, ibl2  !for mpi
               iy = i-ista_kngp+1
               f = qitg_l(i,iq)*ylm(iy)*&
                    &(vlhxc_l(i,1,is)*zfsin(i)+vlhxc_l(i,2,is)*zfcos(i))
!!$                    &(vlhxc_l(i,1,is)*zfsc1(iy)+vlhxc_l(i,2,is)*zfsc2(iy))
               tx = tx + f*ngabc(i,1);  ty = ty + f*ngabc(i,2)
               tz = tz + f*ngabc(i,3)
            end do
         end if
      else
         flchgq = - aimag(zi**(-l3))*shdg_x
         if(kimg == 1) then
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
            do i = ibl1, ibl2  !for mpi
               iy = i-ista_kngp+1
!!$               f  = qitg_l(i,iq)*ylm(iy)*vlhxc_l(i,1,is)*zfsc1(iy)
               f  = qitg_l(i,iq)*ylm(iy)*vlhxc_l(i,1,is)*zfcos(i)
               tx = tx + f*ngabc(i,1);  ty = ty + f*ngabc(i,2)
               tz = tz + f*ngabc(i,3)
            end do
         else if(kimg == 2) then
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
            do i = ibl1, ibl2  !for mpi
               iy = i-ista_kngp+1
               f = qitg_l(i,iq)*ylm(iy)*&
                    &(vlhxc_l(i,1,is)*zfcos(i)-vlhxc_l(i,2,is)*zfsin(i))
!!$                    &(vlhxc_l(i,1,is)*zfsc1(iy)-vlhxc_l(i,2,is)*zfsc2(iy))
               tx = tx + f*ngabc(i,1);  ty = ty + f*ngabc(i,2)
               tz = tz + f*ngabc(i,3)
            end do
         end if
      end if

      if(npes > 1) then
         allocate(f1_mpi(3)); allocate(f2_mpi(3))
         f1_mpi(1) = tx; f1_mpi(2) = ty; f1_mpi(3) = tz  ! MPI
         call mpi_allreduce(f1_mpi,f2_mpi,3 &
              &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr) ! MPI
         tx = f2_mpi(1); ty = f2_mpi(2); tz = f2_mpi(3)
         deallocate(f2_mpi); deallocate(f1_mpi)
      end if
      fx = fx + flchgq*tx*univol
      fy = fy + flchgq*ty*univol
      fz = fz + flchgq*tz*univol

    end subroutine add_hardpart_to_flhxcq_l_core

#ifndef SX
    subroutine add_hardpart_to_flhxcq_l_div(iq,ibl1,ibl2)
      integer, intent(in) :: iq,ibl1,ibl2
      integer       :: i
      real(kind=DP) :: flchgq, f, tx,ty,tz, f2
      real(kind=DP),pointer,dimension(:) :: f1_mpi, f2_mpi

      tx = 0.d0; ty = 0.d0; tz = 0.d0
      if(mod(l3,2) == 0) then
         flchgq = real(zi**(-l3))
         if(kimg == 1) then
            do i = ibl1, ibl2  !for mpi
               f  = qitg_l(i,iq)*ylm_sum(i)*vlhxc_l(i,1,is)*zfsin(i)
               tx = tx + f*ngabc(i,1);  ty = ty + f*ngabc(i,2)
               tz = tz + f*ngabc(i,3)
            end do
         else if(kimg == 2) then
            do i = ibl1, ibl2  !for mpi
               f = qitg_l(i,iq)*ylm_sum(i)*&
                    &(vlhxc_l(i,1,is)*zfsin(i)+vlhxc_l(i,2,is)*zfcos(i))
               tx = tx + f*ngabc(i,1);  ty = ty + f*ngabc(i,2)
               tz = tz + f*ngabc(i,3)
            end do
         end if
      else
         flchgq = - aimag(zi**(-l3))
         if(kimg == 1) then
            do i = ibl1, ibl2  !for mpi
               f  = qitg_l(i,iq)*ylm_sum(i)*vlhxc_l(i,1,is)*zfcos(i)
               tx = tx + f*ngabc(i,1);  ty = ty + f*ngabc(i,2)
               tz = tz + f*ngabc(i,3)
            end do
         else if(kimg == 2) then
            do i = ibl1, ibl2  !for mpi
               f = qitg_l(i,iq)*ylm_sum(i)*&
                    &(vlhxc_l(i,1,is)*zfcos(i)-vlhxc_l(i,2,is)*zfsin(i))
               tx = tx + f*ngabc(i,1);  ty = ty + f*ngabc(i,2)
               tz = tz + f*ngabc(i,3)
            end do
         end if
      end if

!!$      if(npes > 1) then
!!$         allocate(f1_mpi(3)); allocate(f2_mpi(3))
!!$         f1_mpi(1) = tx; f1_mpi(2) = ty; f1_mpi(3) = tz  ! MPI
!!$         call mpi_allreduce(f1_mpi,f2_mpi,3 &
!!$              &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr) ! MPI
!!$         tx = f2_mpi(1); ty = f2_mpi(2); tz = f2_mpi(3)
!!$         deallocate(f2_mpi); deallocate(f1_mpi)
!!$      end if
      fx = fx + flchgq*tx*univol
      fy = fy + flchgq*ty*univol
      fz = fz + flchgq*tz*univol

    end subroutine add_hardpart_to_flhxcq_l_div
#endif
  end subroutine m_Force_term_drv_of_vlhxcQ

#ifdef FORCE_DGEMM
! <-- T.Kokubo & D.Fukata, Feb. 2010
  subroutine m_Force_term_drv_of_flmt(kv3,pos,napt)
    integer, intent(in)                          :: kv3
    real(kind=DP), intent(in), dimension(natm,3) :: pos
    integer,intent(in), dimension(natm,nopr+af)     :: napt

    integer       :: id_sname = -1
    integer :: ibsize,ibl1,ibl2, is,ik,iksnl,ia,it,mdvdb,aoffset
    real(kind=DP) :: fac, fx, fy, fz

    integer, allocatable, dimension(:)       :: mil_tmp
    integer, allocatable, dimension(:,:)     :: ngabc_red
    real(kind=DP),allocatable,dimension(:,:) :: zfsin,  zfcos   ! MPI
    real(kind=DP),allocatable,dimension(:,:,:) :: ar_tmp, ai_tmp  ! MPI
    real(kind=DP),allocatable,dimension(:,:) :: fnlxyz_mpi
    real(kind=DP),allocatable,dimension(:,:,:) :: fr,fi
    real(kind=DP),dimension(3) :: tmpsum
    integer :: ii,ncount

      call tstatc0_begin('m_Force_term_drv_of_flmt ',id_sname,1)
      allocate( gr(np_e,nlmta,3) )    ! MPI
      allocate( gi(np_e,nlmta,3) )    ! MPI
      fac = 2.d0/kv3
      fnlxyz_l = 0.d0

      if(nblocksize_force_is_given) then
         ibsize = nblocksize_force
      else
         ibsize= nb_force_default
      end if
      if(ipriforce >=2) write(nfout,'(" ! ibsize = ",i8," <<m_Force_term_drv_of_flmt>>")') ibsize
      allocate(   mil_tmp(nlmta)    );   mil_tmp = 0
      allocate( ngabc_red(ibsize,3) ); ngabc_red = 0
      allocate(     zfcos(ibsize,natm),      zfsin(ibsize,natm)    )
      allocate(    ar_tmp(ibsize,nlmta,3),  ai_tmp(ibsize,nlmta,3) )
      if(sw_rspace==ON)then
         allocate(fr(np_e,nlmta,ista_k:iend_k));fr=0.d0
         allocate(fi(np_e,nlmta,ista_k:iend_k));fi=0.d0
      endif
      !do is = 1, nspin, af+1
      do is = ista_spin, iend_spin, af+1
         do ik = is, kv3+is-nspin,nspin
            if(map_k(ik) /= myrank_k) cycle
            iksnl = (ik-1)/nspin+1

            if(sw_rspace==ON)then
               do ii=1,3
                  fr=0.d0;fi=0.d0
                  call betar_dot_Psi_4_each_k_in_rs0(6,ista_k,iend_k,ik,zaj_l, &
               &  fr,fi,nmesh_rs_max,snld_rs(1:nmesh_rs_max,1:nlmta,ii))
                  gr(1:np_e,1:nlmta,ii) = fr(1:np_e,1:nlmta,ik)
                  gi(1:np_e,1:nlmta,ii) = fi(1:np_e,1:nlmta,ik)
               enddo
               aoffset = 0
               do ia = 1, natm
                  it = ityp(ia)
                  fx = 0.d0; fy = 0.d0; fz = 0.d0
                  mdvdb = m_PP_include_vanderbilt_pot(it)
                  if(mdvdb == EXECUT) then
                     call term_related_to_drv_of_flmt_VDB ! -> fx,fy,fz
                  else if(mdvdb == SKIP) then
                     call term_related_to_drv_of_flmt_NC  ! -> fx,fy,fz
                  end if
                  aoffset = aoffset + ilmt(it)
                  fnlxyz_l(ia,1) = fnlxyz_l(ia,1) + fx
                  fnlxyz_l(ia,2) = fnlxyz_l(ia,2) + fy
                  fnlxyz_l(ia,3) = fnlxyz_l(ia,3) + fz
               end do
            else
            do ibl1 = 1, iba(ik), ibsize
               ibl2 = ibl1+ibsize-1
               if(ibl2 .gt. iba(ik)) ibl2 = iba(ik)
               gr=0.d0
               gi=0.d0

               call substitute_ngabcred(ik)         ! ngabc,nbase                   --> ngabc_red
               call calc_phase_div()                ! pos, ngabc_red                --> zfcos, zfsin
               call pre_drv_of_betar_dot_WFs_div()  ! zfsin,zfcos,snl               --> ar_tmp,ai_tmp
               call drv_of_betar_dot_WFs_div(gr,gi) ! ar_tmp,ai_tmp,zaj_l,ngabc_red --> gr,gi

               aoffset = 0
               do ia = 1, natm
                  it = ityp(ia)
                  fx = 0.d0; fy = 0.d0; fz = 0.d0
                  mdvdb = m_PP_include_vanderbilt_pot(it)
                  if(mdvdb == EXECUT) then
                     call term_related_to_drv_of_flmt_VDB ! -> fx,fy,fz
                  else if(mdvdb == SKIP) then
                     call term_related_to_drv_of_flmt_NC  ! -> fx,fy,fz
                  end if
                  aoffset = aoffset + ilmt(it)
                  call pucv2cart(rltv,fx,fy,fz)
                  fnlxyz_l(ia,1) = fnlxyz_l(ia,1) + fx
                  fnlxyz_l(ia,2) = fnlxyz_l(ia,2) + fy
                  fnlxyz_l(ia,3) = fnlxyz_l(ia,3) + fz
                end do
            end do
            endif
         end do
      end do

      deallocate( mil_tmp, ngabc_red, zfcos, zfsin, ar_tmp, ai_tmp )
      if(sw_rspace==ON)then
         deallocate(fr)
         deallocate(fi)
      endif

      if(npes>1) then
         allocate(fnlxyz_mpi(natm,3)); fnlxyz_mpi = 0.d0
         call mpi_allreduce(fnlxyz_l,fnlxyz_mpi,3*natm,mpi_double_precision &
             & , mpi_sum, MPI_CommGroup,ierr)
         fnlxyz_l = fnlxyz_mpi
      end if
      if(npes>1)deallocate(fnlxyz_mpi)

      if(af /= 0) then
         allocate(forc_up(natm2,3));allocate(forc_down(natm2,3))
         call forcaf2(natm,natm2,nopr+af,natm,nopr,napt &
           & ,iwei,op,ipriforce,af,fnlxyz_l ,forc_up,forc_down)
         deallocate(forc_down);deallocate(forc_up)
      end if

      deallocate(gi); deallocate(gr)

      if(ipriforce >= 2)  then
         write(nfout,'(" --- fnlxyz_l ->>")')
         do ia = 1, natm
            write(nfout,'(i8,3f20.8)') ia, fnlxyz_l(ia,1:3)
         end do
         write(nfout,'(" <<- fnlxyz_l ---")')
      end if

      call tstatc0_end(id_sname)

    contains
      subroutine term_related_to_drv_of_flmt_NC
        integer :: lmt1, lmta1, lmt2, lmta2, l1,m1,l2,m2, i, lmt11, lmt22
        real(kind=DP) :: fac_n, a, a_n, fr1,fi1,fr2,fi2

        do lmt1 = 1, ilmt(it)
           lmt11 = lmt1 + aoffset
           lmta1 = lmta(lmt1,ia); l1 = ltp(lmt1,it); m1 = mtp(lmt1,it)
           do lmt2 = lmt1, ilmt(it)
              lmt22 = lmt2 + aoffset
              l2 = ltp(lmt2,it); m2 = mtp(lmt2,it); lmta2 = lmta(lmt2,ia)

              if ( ipaw(it) == 0 ) then
                 if(l1 /= l2 .or. m1 /= m2) cycle
                 a = dion(lmt1,lmt2,it)
              else
                 a = dion_paw(lmt1,lmt2,is,ia)
              endif
              fac_n = 2*fac; if(lmt1 == lmt2) fac_n = fac

              if(k_symmetry(ik) == GAMMA) then
                 do i = 1, np_e                                 ! MPI
                    a_n = fac_n * a *occup_l(i,ik)
                    fr1 = fsr_l(i,lmta1,ik); fr2 = fsr_l(i,lmta2,ik)
                    fx  = fx - a_n * (gr(i,lmt11,1)*fr2 + fr1*gr(i,lmt22,1))
                    fy  = fy - a_n * (gr(i,lmt11,2)*fr2 + fr1*gr(i,lmt22,2))
                    fz  = fz - a_n * (gr(i,lmt11,3)*fr2 + fr1*gr(i,lmt22,3))
                 end do
              else
                 do i = 1, np_e                                 ! MPI
                    a_n = fac_n * a *occup_l(i,ik)
                    fr1 = fsr_l(i,lmta1,ik); fr2 = fsr_l(i,lmta2,ik)
                    fi1 = fsi_l(i,lmta1,ik); fi2 = fsi_l(i,lmta2,ik)
                    fx  = fx - a_n * (gr(i,lmt11,1)*fr2+gi(i,lmt11,1)*fi2&
             +fr1*gr(i,lmt22,1)+fi1*gi(i,lmt22,1))
                    fy  = fy - a_n * (gr(i,lmt11,2)*fr2+gi(i,lmt11,2)*fi2&
             +fr1*gr(i,lmt22,2)+fi1*gi(i,lmt22,2))
                    fz  = fz - a_n * (gr(i,lmt11,3)*fr2+gi(i,lmt11,3)*fi2&
             +fr1*gr(i,lmt22,3)+fi1*gi(i,lmt22,3))
                    end do
                 end if
           end do
        end do
      end subroutine term_related_to_drv_of_flmt_NC

      subroutine term_related_to_drv_of_flmt_VDB
        integer :: lmt1, lmta1, lmt2, lmta2, i, lmt11, lmt22
        real(kind=DP) :: fac_n, a, a_n, fr1,fi1,fr2,fi2
        integer       :: id_sname = -1
        call tstatc0_begin('term_related_to_drv_of_flmt_VDB ',id_sname)

        do lmt1 = 1, ilmt(it)
           lmt11 = lmt1 + aoffset
           lmta1 = lmta(lmt1,ia)
           do lmt2 = lmt1, ilmt(it)
              lmt22 = lmt2 + aoffset
              lmta2 = lmta(lmt2,ia)
              fac_n = 2*fac; if(lmt1 == lmt2) fac_n = fac
              if(ipaw(it)==0) then
                  a = vlhxcq(lmt1,lmt2,ia,is)+dion(lmt1,lmt2,it)
              else
                  a = vlhxcq(lmt1,lmt2,ia,is)+dion_paw(lmt1,lmt2,is,ia)
              end if
              if(k_symmetry(ik) == GAMMA) then
                 do i = 1, np_e                                 ! MPI
                    a_n = fac_n*(a - eko_l(i,ik)*q(lmt1,lmt2,it))*occup_l(i,ik)
                    fr1 = fsr_l(i,lmta1,ik); fr2 = fsr_l(i,lmta2,ik)
                    fx  = fx - a_n * (gr(i,lmt11,1)*fr2 + fr1*gr(i,lmt22,1))
                    fy  = fy - a_n * (gr(i,lmt11,2)*fr2 + fr1*gr(i,lmt22,2))
                    fz  = fz - a_n * (gr(i,lmt11,3)*fr2 + fr1*gr(i,lmt22,3))
                 end do
              else
                 do i = 1, np_e                                 ! MPI
                    a_n = fac_n*(a - eko_l(i,ik)*q(lmt1,lmt2,it))*occup_l(i,ik)
                    fr1 = fsr_l(i,lmta1,ik); fr2 = fsr_l(i,lmta2,ik)
                    fi1 = fsi_l(i,lmta1,ik); fi2 = fsi_l(i,lmta2,ik)
                    fx  = fx - a_n * (gr(i,lmt11,1)*fr2+gi(i,lmt11,1)*fi2&
             +fr1*gr(i,lmt22,1)+fi1*gi(i,lmt22,1))
                    fy  = fy - a_n * (gr(i,lmt11,2)*fr2+gi(i,lmt11,2)*fi2&
             +fr1*gr(i,lmt22,2)+fi1*gi(i,lmt22,2))
                    fz  = fz - a_n * (gr(i,lmt11,3)*fr2+gi(i,lmt11,3)*fi2&
             +fr1*gr(i,lmt22,3)+fi1*gi(i,lmt22,3))
                 end do
              end if
           end do
        end do
        call tstatc0_end(id_sname)
      end subroutine term_related_to_drv_of_flmt_VDB

      subroutine substitute_ngabcred(ik)
        integer, intent(in) :: ik
        integer :: ig, i1
        if(ibl2-ibl1+1 > ibsize) call phase_error_with_msg(nfout,' ibl2-ibl1+1 > ibsize',__LINE__,__FILE__)
           do ig = 1, ibl2-ibl1+1
              i1 = nbase(ig+ibl1-1,ik)
              ngabc_red(ig,1) = ngabc(i1,1)
              ngabc_red(ig,2) = ngabc(i1,2)
              ngabc_red(ig,3) = ngabc(i1,3)
           end do
      end subroutine substitute_ngabcred

      subroutine calc_phase_div()
          integer :: ia, ig
          real(kind=DP) :: fx,fy,fz, ph

          do ia = 1, natm
             fx = pos(ia,1)*PAI2
             fy = pos(ia,2)*PAI2
             fz = pos(ia,3)*PAI2
             do ig = 1, ibl2-ibl1+1
                ph = ngabc_red(ig,1)*fx + ngabc_red(ig,2)*fy + ngabc_red(ig,3)*fz
                zfcos(ig,ia) = dcos(ph)
                zfsin(ig,ia) = dsin(ph)
             end do
          end do
      end subroutine calc_phase_div

      subroutine pre_drv_of_betar_dot_WFs_div()
          integer :: ia, i, t, lmt1, lmtt1, icnt, mil, fac_tmp

          icnt = 0
          if(k_symmetry(ik) == GAMMA) then
             if(kimg == 1) then   ! only ar
               do ia = 1, natm
                  it = ityp(ia)
                  do lmt1 = 1, ilmt(it)
                     icnt = icnt + 1
                     lmtt1 = lmtt(lmt1,it)
                     do i = ibl1, ibl2
                        ar_tmp(i-ibl1+1,icnt,1) = zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,1)
                        ar_tmp(i-ibl1+1,icnt,2) = zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,2)
                        ar_tmp(i-ibl1+1,icnt,3) = zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,3)
                     end do
                  enddo
               enddo
             else if(kimg == 2) then
               do ia = 1, natm
                  it = ityp(ia)
                  do lmt1 = 1, ilmt(it)
                     icnt = icnt + 1
                     lmtt1 = lmtt(lmt1,it)
                     mil = mod( ltp(lmt1,it),4 )
                     if( mil.eq.1 .or. mil.eq.2 ) then
                       fac_tmp= -2.d0
                     else if( mil.eq.0 .or. mil.eq.3 ) then
                       fac_tmp=  2.d0
                     endif

                     if( mil.eq.1 .or. mil.eq.3 ) then
                         do i = ibl1, ibl2
                            ar_tmp(i-ibl1+1,icnt,1) =  fac_tmp*zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,1)
                            ar_tmp(i-ibl1+1,icnt,2) =  fac_tmp*zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,2)
                            ar_tmp(i-ibl1+1,icnt,3) =  fac_tmp*zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,3)
                            ai_tmp(i-ibl1+1,icnt,1) =  fac_tmp*zfsin(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,1)
                            ai_tmp(i-ibl1+1,icnt,2) =  fac_tmp*zfsin(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,2)
                            ai_tmp(i-ibl1+1,icnt,3) =  fac_tmp*zfsin(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,3)
                         end do
                     else if( mil.eq.0 .or. mil.eq.2 ) then
                         do i = ibl1, ibl2
                            ar_tmp(i-ibl1+1,icnt,1) = -fac_tmp*zfsin(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,1)
                            ar_tmp(i-ibl1+1,icnt,2) = -fac_tmp*zfsin(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,2)
                            ar_tmp(i-ibl1+1,icnt,3) = -fac_tmp*zfsin(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,3)
                            ai_tmp(i-ibl1+1,icnt,1) =  fac_tmp*zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,1)
                            ai_tmp(i-ibl1+1,icnt,2) =  fac_tmp*zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,2)
                            ai_tmp(i-ibl1+1,icnt,3) =  fac_tmp*zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,3)
                         end do
                     endif
                  enddo
                enddo
             endif
          else
               do ia = 1, natm
                  it = ityp(ia)
                  do lmt1 = 1, ilmt(it)
                     icnt = icnt + 1
                     lmtt1 = lmtt(lmt1,it)
                     mil_tmp(icnt) = mod( ltp(lmt1,it),4 )
                      do i = ibl1, ibl2
                         ar_tmp(i-ibl1+1,icnt,1) = zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,1)
                         ar_tmp(i-ibl1+1,icnt,2) = zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,2)
                         ar_tmp(i-ibl1+1,icnt,3) = zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,3)
                         ai_tmp(i-ibl1+1,icnt,1) = zfsin(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,1)
                         ai_tmp(i-ibl1+1,icnt,2) = zfsin(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,2)
                         ai_tmp(i-ibl1+1,icnt,3) = zfsin(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,3)
                      end do
                  enddo
                enddo
          endif
      end subroutine pre_drv_of_betar_dot_WFs_div

      subroutine drv_of_betar_dot_WFs_div(gr,gi)
        real(kind=DP),intent(inout), dimension(np_e,nlmta,3) :: gr,gi
        integer :: lmt1, ib, i, mil, ibl1_2
        real(kind=DP) :: br, bi, tmp1, tmp2, tmp3
        integer :: id_sname = -1
        call tstatc0_begin('drv_of_betar_dot_WFs_div ',id_sname)

         if(k_symmetry(ik) == GAMMA) then
            if(kimg == 1) then
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ar_tmp(1,1,1),ibsize, 0.d0,gr(1,1,1),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ar_tmp(1,1,2),ibsize, 0.d0,gr(1,1,2),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ar_tmp(1,1,3),ibsize, 0.d0,gr(1,1,3),np_e)
            else if(kimg == 2) then
               ibl1_2 = max(2,ibl1)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1_2+1, 1.d0,zaj_l(ibl1_2,1,ik,1),kg1,ai_tmp(ibl1_2-ibl1+1,1,1),ibsize, &
              &           0.d0,gr(1,1,1),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1_2+1, 1.d0,zaj_l(ibl1_2,1,ik,2),kg1,ar_tmp(ibl1_2-ibl1+1,1,1),ibsize, &
              &           1.d0,gr(1,1,1),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1_2+1, 1.d0,zaj_l(ibl1_2,1,ik,1),kg1,ai_tmp(ibl1_2-ibl1+1,1,2),ibsize, &
              &           0.d0,gr(1,1,2),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1_2+1, 1.d0,zaj_l(ibl1_2,1,ik,2),kg1,ar_tmp(ibl1_2-ibl1+1,1,2),ibsize, &
              &           1.d0,gr(1,1,2),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1_2+1, 1.d0,zaj_l(ibl1_2,1,ik,1),kg1,ai_tmp(ibl1_2-ibl1+1,1,3),ibsize, &
              &           0.d0,gr(1,1,3),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1_2+1, 1.d0,zaj_l(ibl1_2,1,ik,2),kg1,ar_tmp(ibl1_2-ibl1+1,1,3),ibsize, &
              &           1.d0,gr(1,1,3),np_e)
            endif
         else
            if(kimg == 1) then
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ar_tmp(1,1,1),ibsize, 0.d0,gr(1,1,1),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ar_tmp(1,1,2),ibsize, 0.d0,gr(1,1,2),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ar_tmp(1,1,3),ibsize, 0.d0,gr(1,1,3),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ai_tmp(1,1,1),ibsize, 0.d0,gi(1,1,1),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ai_tmp(1,1,2),ibsize, 0.d0,gi(1,1,2),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ai_tmp(1,1,3),ibsize, 0.d0,gi(1,1,3),np_e)
            else if(kimg == 2) then
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ar_tmp(1,1,1),ibsize, 0.d0,gr(1,1,1),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1,-1.d0,zaj_l(ibl1,1,ik,2),kg1,ai_tmp(1,1,1),ibsize, 1.d0,gr(1,1,1),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ar_tmp(1,1,2),ibsize, 0.d0,gr(1,1,2),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1,-1.d0,zaj_l(ibl1,1,ik,2),kg1,ai_tmp(1,1,2),ibsize, 1.d0,gr(1,1,2),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ar_tmp(1,1,3),ibsize, 0.d0,gr(1,1,3),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1,-1.d0,zaj_l(ibl1,1,ik,2),kg1,ai_tmp(1,1,3),ibsize, 1.d0,gr(1,1,3),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ai_tmp(1,1,1),ibsize, 0.d0,gi(1,1,1),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,2),kg1,ar_tmp(1,1,1),ibsize, 1.d0,gi(1,1,1),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ai_tmp(1,1,2),ibsize, 0.d0,gi(1,1,2),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,2),kg1,ar_tmp(1,1,2),ibsize, 1.d0,gi(1,1,2),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ai_tmp(1,1,3),ibsize, 0.d0,gi(1,1,3),np_e)
               call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,2),kg1,ar_tmp(1,1,3),ibsize, 1.d0,gi(1,1,3),np_e)
            end if
            do lmt1 = 1, nlmta
               mil = mil_tmp(lmt1)
               if( mil == 1 ) then
                  do ib = 1, np_e                              ! MPI
                     tmp1 = gr(ib,lmt1,1)
                     tmp2 = gr(ib,lmt1,2)
                     tmp3 = gr(ib,lmt1,3)
                     gr(ib,lmt1,1) = -gi(ib,lmt1,1)
                     gr(ib,lmt1,2) = -gi(ib,lmt1,2)
                     gr(ib,lmt1,3) = -gi(ib,lmt1,3)
                     gi(ib,lmt1,1) = tmp1
                     gi(ib,lmt1,2) = tmp2
                     gi(ib,lmt1,3) = tmp3
                  end do
               elseif (mil == 2 ) then
                  do ib = 1, np_e                              ! MPI
                     gr(ib,lmt1,1) = - gr(ib,lmt1,1)
                     gr(ib,lmt1,2) = - gr(ib,lmt1,2)
                     gr(ib,lmt1,3) = - gr(ib,lmt1,3)
                     gi(ib,lmt1,1) = - gi(ib,lmt1,1)
                     gi(ib,lmt1,2) = - gi(ib,lmt1,2)
                     gi(ib,lmt1,3) = - gi(ib,lmt1,3)
                  end do
               elseif (mil == 3 ) then
                  do ib = 1, np_e                              ! MPI
                     tmp1 = gr(ib,lmt1,1)
                     tmp2 = gr(ib,lmt1,2)
                     tmp3 = gr(ib,lmt1,3)
                     gr(ib,lmt1,1) = gi(ib,lmt1,1)
                     gr(ib,lmt1,2) = gi(ib,lmt1,2)
                     gr(ib,lmt1,3) = gi(ib,lmt1,3)
                     gi(ib,lmt1,1) = -tmp1
                     gi(ib,lmt1,2) = -tmp2
                     gi(ib,lmt1,3) = -tmp3
                  end do
               endif
            enddo
         end if
         call tstatc0_end(id_sname)
      end subroutine drv_of_betar_dot_WFs_div
  end subroutine m_Force_term_drv_of_flmt
! -->

! ============================ added by K. Tagami ==================== 11.0
  subroutine m_Force_term_drv_of_flmt_noncl( kv3,pos,napt )
    integer, intent(in)                          :: kv3
    real(kind=DP), intent(in), dimension(natm,3) :: pos
    integer,intent(in), dimension(natm,nopr+af)     :: napt

    integer       :: id_sname = -1
    integer :: ibsize,ibl1,ibl2, is,ik,iksnl,ia,it,mdvdb,aoffset
    real(kind=DP) :: fac, fx, fy, fz

    integer, allocatable, dimension(:)       :: mil_tmp
    integer, allocatable, dimension(:,:)     :: ngabc_red
    real(kind=DP),allocatable,dimension(:,:) :: zfsin,  zfcos   ! MPI
    real(kind=DP),allocatable,dimension(:,:,:) :: ar_tmp, ai_tmp  ! MPI
    real(kind=DP),allocatable,dimension(:,:) :: fnlxyz_mpi

    real(kind=DP), allocatable :: gr_noncl( :,:,:,: )
    real(kind=DP), allocatable :: gi_noncl( :,:,:,: )

    call tstatc0_begin('m_Force_term_drv_of_flmt_noncl ',id_sname,1)

    allocate( gr_noncl(np_e,nlmta,3, ndim_spinor ) )    ! MPI
    allocate( gi_noncl(np_e,nlmta,3, ndim_spinor ) )    ! MPI

!!!!!!!!!!!    fac = 2.d0 / (kv3 /ndim_spinor)
    fac = 1.d0 / (kv3 /ndim_spinor)

    fnlxyz_l = 0.d0

    if(nblocksize_force_is_given) then
       ibsize = nblocksize_force
    else
       ibsize= nb_force_default
    end if
    write(nfout,'(" ! ibsize = ",i8," <<m_Force_term_drv_of_flmt_noncl>>")') ibsize
    allocate(   mil_tmp(nlmta)    );   mil_tmp = 0
    allocate( ngabc_red(ibsize,3) ); ngabc_red = 0
    allocate(     zfcos(ibsize,natm),      zfsin(ibsize,natm)    )
    allocate(    ar_tmp(ibsize,nlmta,3),  ai_tmp(ibsize,nlmta,3) )

    do ik = 1, kv3, ndim_spinor
       if(map_k(ik) /= myrank_k) cycle
       iksnl = (ik-1)/nspin+1

       do ibl1 = 1, iba(ik), ibsize
          ibl2 = ibl1+ibsize-1
          if(ibl2 .gt. iba(ik)) ibl2 = iba(ik)

          gr_noncl = 0.d0;     gi_noncl = 0.d0

          call substitute_ngabcred(ik)         ! ngabc,nbase      --> ngabc_red
          call calc_phase_div()                ! pos, ngabc_red  --> zfcos, zfsin
          call pre_drv_of_betar_dot_WFs_div()  ! zfsin,zfcos,snl               --> ar_tmp,ai_tmp
          Do is=1, ndim_spinor
             call drv_of_betar_dot_WFs_div( ik+is-1, gr_noncl(:,:,:,is), &
                  &                                  gi_noncl(:,:,:,is) )
                            ! ar_tmp,ai_tmp,zaj_l,ngabc_red --> gr,gi
          End do

          aoffset = 0
          do ia = 1, natm
             it = ityp(ia)
             fx = 0.d0; fy = 0.d0; fz = 0.d0
             mdvdb = m_PP_include_vanderbilt_pot(it)

             if(mdvdb == EXECUT) then
                call term_related_to_drv_of_flmt_VDB ! -> fx,fy,fz
             else if(mdvdb == SKIP) then
                call term_related_to_drv_of_flmt_NC  ! -> fx,fy,fz
             end if

             aoffset = aoffset + ilmt(it)
             call pucv2cart(rltv,fx,fy,fz)
             fnlxyz_l(ia,1) = fnlxyz_l(ia,1) + fx
             fnlxyz_l(ia,2) = fnlxyz_l(ia,2) + fy
             fnlxyz_l(ia,3) = fnlxyz_l(ia,3) + fz
          end do

       end do
    end do

    deallocate( mil_tmp, ngabc_red, zfcos, zfsin, ar_tmp, ai_tmp )

    if(npes>1) then
       allocate(fnlxyz_mpi(natm,3)); fnlxyz_mpi = 0.d0
       call mpi_allreduce(fnlxyz_l,fnlxyz_mpi,3*natm,mpi_double_precision &
            & , mpi_sum, MPI_CommGroup,ierr)
       fnlxyz_l = fnlxyz_mpi
    end if
    if(npes>1)deallocate(fnlxyz_mpi)

    if(af /= 0) then
       allocate(forc_up(natm2,3));allocate(forc_down(natm2,3))
       call forcaf2(natm,natm2,nopr+af,natm,nopr,napt &
            & ,iwei,op,ipriforce,af,fnlxyz_l ,forc_up,forc_down)
       deallocate(forc_down);deallocate(forc_up)
    end if

    deallocate(gi_noncl); deallocate(gr_noncl)

    if(ipriforce >= 2)  then
       write(nfout,'(" --- fnlxyz_l ->>")')
       do ia = 1, natm
          write(nfout,'(i8,3f20.8)') ia, fnlxyz_l(ia,1:3)
       end do
       write(nfout,'(" <<- fnlxyz_l ---")')
    end if

    call tstatc0_end(id_sname)

  contains

    subroutine term_related_to_drv_of_flmt_NC
      integer :: lmt1, lmta1, lmt2, lmta2, l1,m1,l2,m2, i, lmt11, lmt22
      real(kind=DP) :: fac_n, a, a_n, fr1,fi1,fr2,fi2

      integer :: is1, is2, istmp
      complex(kind=CMPLDP) :: z1, z2, z3
      complex(kind=CMPLDP) :: zg1(3), zg2(3), zsum(3)

      do lmt1 = 1, ilmt(it)
         lmt11 = lmt1 + aoffset
         lmta1 = lmta(lmt1,ia);

         l1 = ltp(lmt1,it); m1 = mtp(lmt1,it)

! ==== KT_mod ======================== 2013/10/31
!!!!!!!!         do lmt2 = lmt1, ilmt(it)
         do lmt2 = 1, ilmt(it)
! ==================================== 2013/10/31

            lmt22 = lmt2 + aoffset
            lmta2 = lmta(lmt2,ia)

            l2 = ltp(lmt2,it); m2 = mtp(lmt2,it);

!!!!!!            if(l1 /= l2 .or. m1 /= m2) cycle

! ==== KT_mod ======================== 2013/10/31
!!!            fac_n = 2*fac; if(lmt1 == lmt2) fac_n = fac
            fac_n = fac
! ==================================== 2013/10/31

            if (k_symmetry(ik) == GAMMA) then
               call phase_error_with_msg(nfout,'Not supported : Gamma symmetry in non-collinear system',__LINE__,__FILE__)

            else

               do i = 1, np_e                                 ! MPI
                  zsum = 0.0d0
                  Do is1=1, ndim_spinor
                     Do is2=1, ndim_spinor
                        istmp = ( is1 -1 )*ndim_spinor + is2

                        z1 = dcmplx( fsr_l(i,lmta1,ik+is1-1), &
                             &       fsi_l(i,lmta1,ik+is1-1 ) )
                        z2 = dcmplx( fsr_l(i,lmta2,ik+is2-1), &
                             &       fsi_l(i,lmta2,ik+is2-1 ) )
                        zg1(:) = dcmplx( gr_noncl(i,lmt11,:,is1), &
                             &           gi_noncl(i,lmt11,:,is1) )
                        zg2(:) = dcmplx( gr_noncl(i,lmt22,:,is2), &
                             &           gi_noncl(i,lmt22,:,is2) )

                        zsum(:) = zsum(:) + (conjg(z1)*zg2(:) + conjg(zg1(:))*z2 )&
                             &              *dion_scr_noncl( lmt1,lmt2,istmp,ia )
                     End do
                  End Do

                  zsum = zsum *fac_n *occup_l(i,ik)

                  fx  = fx - real( zsum(1) )
                  fy  = fy - real( zsum(2) )
                  fz  = fz - real( zsum(3) )
               end do
            end if
         end do
      end do
    end subroutine term_related_to_drv_of_flmt_NC

    subroutine term_related_to_drv_of_flmt_VDB
      integer :: lmt1, lmta1, lmt2, lmta2, i, lmt11, lmt22
      real(kind=DP) :: fac_n, a, a_n, fr1,fi1,fr2,fi2

      integer :: is1, is2, istmp
      complex(kind=CMPLDP) :: z1, z2, z3
      complex(kind=CMPLDP) :: zg1(3), zg2(3), zsum(3)
      integer       :: id_sname = -1

      call tstatc0_begin('term_related_to_drv_of_flmt_VDB ',id_sname)

      do lmt1 = 1, ilmt(it)
         lmt11 = lmt1 + aoffset
         lmta1 = lmta(lmt1,ia)

! ==== KT_mod ======================== 2013/10/31
!!!         do lmt2 = lmt1, ilmt(it)
         do lmt2 = 1, ilmt(it)
! ==================================== 2013/10/31

            lmt22 = lmt2 + aoffset
            lmta2 = lmta(lmt2,ia)

! ==== KT_mod ======================== 2013/10/31
!!!            fac_n = 2*fac; if(lmt1 == lmt2) fac_n = fac
            fac_n = fac
! ==================================== 2013/10/31

            if(k_symmetry(ik) == GAMMA) then
               call phase_error_with_msg(nfout,'Not supported : Gamma symmetry in non-collinear system',__LINE__,__FILE__)
            else

               do i = 1, np_e                                 ! MPI
                  zsum = 0.0d0
                  Do is1=1, ndim_spinor
                     Do is2=1, ndim_spinor
                        istmp = ( is1 -1 )*ndim_spinor + is2

                        z1 = dcmplx( fsr_l(i,lmta1,ik+is1-1), &
                             &       fsi_l(i,lmta1,ik+is1-1 ) )
                        z2 = dcmplx( fsr_l(i,lmta2,ik+is2-1), &
                             &       fsi_l(i,lmta2,ik+is2-1) )
                        zg1(:) = dcmplx( gr_noncl(i,lmt11,:,is1), &
                             &           gi_noncl(i,lmt11,:,is1) )
                        zg2(:) = dcmplx( gr_noncl(i,lmt22,:,is2), &
                             &           gi_noncl(i,lmt22,:,is2) )

                        z3 = dion_scr_noncl( lmt1,lmt2,istmp,ia ) &
                             &  -eko_l(i,ik) *q_noncl(lmt1,lmt2,istmp,it )

                        zsum(:) = zsum(:) + (conjg(z1)*zg2(:) + conjg(zg1(:))*z2 ) *z3
                     End do
                  End Do

                  zsum = zsum *fac_n *occup_l(i,ik)

                  fx  = fx - real( zsum(1) )
                  fy  = fy - real( zsum(2) )
                  fz  = fz - real( zsum(3) )
               end do

            end if
         end do
      end do
      call tstatc0_end(id_sname)
    end subroutine term_related_to_drv_of_flmt_VDB

    subroutine substitute_ngabcred(ik)
      integer, intent(in) :: ik
      integer :: ig, i1
      if(ibl2-ibl1+1 > ibsize) call phase_error_with_msg(nfout,' ibl2-ibl1+1 > ibsize',__LINE__,__FILE__)
      do ig = 1, ibl2-ibl1+1
         i1 = nbase(ig+ibl1-1,ik)
         ngabc_red(ig,1) = ngabc(i1,1)
         ngabc_red(ig,2) = ngabc(i1,2)
         ngabc_red(ig,3) = ngabc(i1,3)
      end do
    end subroutine substitute_ngabcred

    subroutine calc_phase_div()
      integer :: ia, ig
      real(kind=DP) :: fx,fy,fz, ph

      do ia = 1, natm
         fx = pos(ia,1)*PAI2
         fy = pos(ia,2)*PAI2
         fz = pos(ia,3)*PAI2
         do ig = 1, ibl2-ibl1+1
            ph = ngabc_red(ig,1)*fx + ngabc_red(ig,2)*fy + ngabc_red(ig,3)*fz
            zfcos(ig,ia) = dcos(ph)
            zfsin(ig,ia) = dsin(ph)
         end do
      end do
    end subroutine calc_phase_div

    subroutine pre_drv_of_betar_dot_WFs_div()
      integer :: ia, i, t, lmt1, lmtt1, icnt, mil, fac_tmp

      icnt = 0
      if(k_symmetry(ik) == GAMMA) then
         if(kimg == 1) then   ! only ar
            do ia = 1, natm
               it = ityp(ia)
               do lmt1 = 1, ilmt(it)
                  icnt = icnt + 1
                  lmtt1 = lmtt(lmt1,it)
                  do i = ibl1, ibl2
                     ar_tmp(i-ibl1+1,icnt,1) = zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,1)
                     ar_tmp(i-ibl1+1,icnt,2) = zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,2)
                     ar_tmp(i-ibl1+1,icnt,3) = zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,3)
                  end do
               enddo
            enddo
         else if(kimg == 2) then
            do ia = 1, natm
               it = ityp(ia)
               do lmt1 = 1, ilmt(it)
                  icnt = icnt + 1
                  lmtt1 = lmtt(lmt1,it)
                  mil = mod( ltp(lmt1,it),4 )
                  if( mil.eq.1 .or. mil.eq.2 ) then
                     fac_tmp= -2.d0
                  else if( mil.eq.0 .or. mil.eq.3 ) then
                     fac_tmp=  2.d0
                  endif

                  if( mil.eq.1 .or. mil.eq.3 ) then
                     do i = ibl1, ibl2
                        ar_tmp(i-ibl1+1,icnt,1) =  fac_tmp*zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,1)
                        ar_tmp(i-ibl1+1,icnt,2) =  fac_tmp*zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,2)
                        ar_tmp(i-ibl1+1,icnt,3) =  fac_tmp*zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,3)
                        ai_tmp(i-ibl1+1,icnt,1) =  fac_tmp*zfsin(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,1)
                        ai_tmp(i-ibl1+1,icnt,2) =  fac_tmp*zfsin(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,2)
                        ai_tmp(i-ibl1+1,icnt,3) =  fac_tmp*zfsin(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,3)
                     end do
                  else if( mil.eq.0 .or. mil.eq.2 ) then
                     do i = ibl1, ibl2
                        ar_tmp(i-ibl1+1,icnt,1) = -fac_tmp*zfsin(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,1)
                        ar_tmp(i-ibl1+1,icnt,2) = -fac_tmp*zfsin(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,2)
                        ar_tmp(i-ibl1+1,icnt,3) = -fac_tmp*zfsin(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,3)
                        ai_tmp(i-ibl1+1,icnt,1) =  fac_tmp*zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,1)
                        ai_tmp(i-ibl1+1,icnt,2) =  fac_tmp*zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,2)
                        ai_tmp(i-ibl1+1,icnt,3) =  fac_tmp*zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,3)
                     end do
                  endif
               enddo
            enddo
         endif
      else
         do ia = 1, natm
            it = ityp(ia)
            do lmt1 = 1, ilmt(it)
               icnt = icnt + 1
               lmtt1 = lmtt(lmt1,it)
               mil_tmp(icnt) = mod( ltp(lmt1,it),4 )
               do i = ibl1, ibl2
                  ar_tmp(i-ibl1+1,icnt,1) = zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,1)
                  ar_tmp(i-ibl1+1,icnt,2) = zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,2)
                  ar_tmp(i-ibl1+1,icnt,3) = zfcos(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,3)
                  ai_tmp(i-ibl1+1,icnt,1) = zfsin(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,1)
                  ai_tmp(i-ibl1+1,icnt,2) = zfsin(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,2)
                  ai_tmp(i-ibl1+1,icnt,3) = zfsin(i-ibl1+1,ia)*snl(i,lmtt1,iksnl)*ngabc_red(i-ibl1+1,3)
               end do
            enddo
         enddo
      endif
    end subroutine pre_drv_of_betar_dot_WFs_div

    subroutine drv_of_betar_dot_WFs_div( ik, gr, gi )
      integer, intent(in) :: ik

      real(kind=DP),intent(inout), dimension(np_e,nlmta,3) :: gr,gi
      integer :: lmt1, ib, i, mil, ibl1_2
      real(kind=DP) :: br, bi, tmp1, tmp2, tmp3
      integer :: id_sname = -1

      call tstatc0_begin('drv_of_betar_dot_WFs_div ',id_sname)

      if(k_symmetry(ik) == GAMMA) then
         if(kimg == 1) then
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ar_tmp(1,1,1),ibsize, 0.d0,gr(1,1,1),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ar_tmp(1,1,2),ibsize, 0.d0,gr(1,1,2),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ar_tmp(1,1,3),ibsize, 0.d0,gr(1,1,3),np_e)
         else if(kimg == 2) then
            ibl1_2 = max(2,ibl1)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1_2+1, 1.d0,zaj_l(ibl1_2,1,ik,1),kg1,ai_tmp(ibl1_2-ibl1+1,1,1),ibsize, &
                 &           0.d0,gr(1,1,1),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1_2+1, 1.d0,zaj_l(ibl1_2,1,ik,2),kg1,ar_tmp(ibl1_2-ibl1+1,1,1),ibsize, &
                 &           1.d0,gr(1,1,1),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1_2+1, 1.d0,zaj_l(ibl1_2,1,ik,1),kg1,ai_tmp(ibl1_2-ibl1+1,1,2),ibsize, &
                 &           0.d0,gr(1,1,2),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1_2+1, 1.d0,zaj_l(ibl1_2,1,ik,2),kg1,ar_tmp(ibl1_2-ibl1+1,1,2),ibsize, &
                 &           1.d0,gr(1,1,2),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1_2+1, 1.d0,zaj_l(ibl1_2,1,ik,1),kg1,ai_tmp(ibl1_2-ibl1+1,1,3),ibsize, &
                 &           0.d0,gr(1,1,3),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1_2+1, 1.d0,zaj_l(ibl1_2,1,ik,2),kg1,ar_tmp(ibl1_2-ibl1+1,1,3),ibsize, &
                 &           1.d0,gr(1,1,3),np_e)
         endif
      else
         if(kimg == 1) then
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ar_tmp(1,1,1),ibsize, 0.d0,gr(1,1,1),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ar_tmp(1,1,2),ibsize, 0.d0,gr(1,1,2),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ar_tmp(1,1,3),ibsize, 0.d0,gr(1,1,3),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ai_tmp(1,1,1),ibsize, 0.d0,gi(1,1,1),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ai_tmp(1,1,2),ibsize, 0.d0,gi(1,1,2),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ai_tmp(1,1,3),ibsize, 0.d0,gi(1,1,3),np_e)
         else if(kimg == 2) then
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ar_tmp(1,1,1),ibsize, 0.d0,gr(1,1,1),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1,-1.d0,zaj_l(ibl1,1,ik,2),kg1,ai_tmp(1,1,1),ibsize, 1.d0,gr(1,1,1),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ar_tmp(1,1,2),ibsize, 0.d0,gr(1,1,2),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1,-1.d0,zaj_l(ibl1,1,ik,2),kg1,ai_tmp(1,1,2),ibsize, 1.d0,gr(1,1,2),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ar_tmp(1,1,3),ibsize, 0.d0,gr(1,1,3),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1,-1.d0,zaj_l(ibl1,1,ik,2),kg1,ai_tmp(1,1,3),ibsize, 1.d0,gr(1,1,3),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ai_tmp(1,1,1),ibsize, 0.d0,gi(1,1,1),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,2),kg1,ar_tmp(1,1,1),ibsize, 1.d0,gi(1,1,1),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ai_tmp(1,1,2),ibsize, 0.d0,gi(1,1,2),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,2),kg1,ar_tmp(1,1,2),ibsize, 1.d0,gi(1,1,2),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,1),kg1,ai_tmp(1,1,3),ibsize, 0.d0,gi(1,1,3),np_e)
            call DGEMM('T','N',np_e,nlmta,ibl2-ibl1+1, 1.d0,zaj_l(ibl1,1,ik,2),kg1,ar_tmp(1,1,3),ibsize, 1.d0,gi(1,1,3),np_e)
         end if
         do lmt1 = 1, nlmta
            mil = mil_tmp(lmt1)
            if( mil == 1 ) then
               do ib = 1, np_e                              ! MPI
                  tmp1 = gr(ib,lmt1,1)
                  tmp2 = gr(ib,lmt1,2)
                  tmp3 = gr(ib,lmt1,3)
                  gr(ib,lmt1,1) = -gi(ib,lmt1,1)
                  gr(ib,lmt1,2) = -gi(ib,lmt1,2)
                  gr(ib,lmt1,3) = -gi(ib,lmt1,3)
                  gi(ib,lmt1,1) = tmp1
                  gi(ib,lmt1,2) = tmp2
                  gi(ib,lmt1,3) = tmp3
               end do
            elseif (mil == 2 ) then
               do ib = 1, np_e                              ! MPI
                  gr(ib,lmt1,1) = - gr(ib,lmt1,1)
                  gr(ib,lmt1,2) = - gr(ib,lmt1,2)
                  gr(ib,lmt1,3) = - gr(ib,lmt1,3)
                  gi(ib,lmt1,1) = - gi(ib,lmt1,1)
                  gi(ib,lmt1,2) = - gi(ib,lmt1,2)
                  gi(ib,lmt1,3) = - gi(ib,lmt1,3)
               end do
            elseif (mil == 3 ) then
               do ib = 1, np_e                              ! MPI
                  tmp1 = gr(ib,lmt1,1)
                  tmp2 = gr(ib,lmt1,2)
                  tmp3 = gr(ib,lmt1,3)
                  gr(ib,lmt1,1) = gi(ib,lmt1,1)
                  gr(ib,lmt1,2) = gi(ib,lmt1,2)
                  gr(ib,lmt1,3) = gi(ib,lmt1,3)
                  gi(ib,lmt1,1) = -tmp1
                  gi(ib,lmt1,2) = -tmp2
                  gi(ib,lmt1,3) = -tmp3
               end do
            endif
         enddo
      end if
      call tstatc0_end(id_sname)
    end subroutine drv_of_betar_dot_WFs_div

  end subroutine m_Force_term_drv_of_flmt_noncl
! ========================================================================= 11.0

#else
  subroutine m_Force_term_drv_of_flmt(kv3,pos,napt)
    integer, intent(in)                          :: kv3
    real(kind=DP), intent(in), dimension(natm,3) :: pos
    integer,intent(in), dimension(natm,nopr+af)     :: napt

    integer  :: ia,is,ik,iksnl,it,lmt1,lmtt1,il,mdvdb
    real(kind=DP) :: fac, fx, fy, fz
#ifdef SX
    real(kind=DP),allocatable,dimension(:) :: fxyz_mpi1,fxyz_mpi2  ! MPI
    real(kind=DP),allocatable,dimension(:) :: zfsin, zfcos  ! MPI
#else
    real(kind=DP),allocatable,dimension(:) :: zfsin, zfcos  ! MPI
    integer :: ibl1,ibl2,ibsize,ncache
!!$    real(kind=DP),allocatable,dimension(:,:,:) :: zaj_red ! d(ibsize,np_e,kimg)
!!$    real(kind=DP),allocatable,dimension(:,:) :: fnlxyz_mpi, snl_red
    real(kind=DP),allocatable,dimension(:,:) :: fnlxyz_mpi
    real(kind=DP),allocatable,dimension(:) :: zfsin_red,zfcos_red  ! MPI
    integer, allocatable, dimension(:,:) :: ngabc_red
#endif
    integer       :: id_sname = -1
    call tstatc0_begin('m_Force_term_drv_of_flmt ',id_sname,1)

    allocate(gr(np_e,nlmta,3));allocate(gi(np_e,nlmta,3))  ! MPI

    fac = 2.d0/kv3
#ifdef SX
    allocate(ar(kg1))        ;allocate(ai(kg1))
    allocate(zfsin(nbmx)); zfsin = 0.d0
    allocate(zfcos(nbmx)); zfcos = 0.d0
    if(npes>1) then
       allocate(fxyz_mpi1(3));    allocate(fxyz_mpi2(3))      ! MPI
    end if
    do ia = 1, natm
       fx = 0.d0; fy = 0.d0; fz = 0.d0
!!$       call calc_phase2(natm,pos,ia,kgp,nbmx,ngabc,nbmx,zfcos,zfsin)
       it = ityp(ia)
       mdvdb = m_PP_include_vanderbilt_pot(it)
       do is = 1, nspin, af+1
          do ik = is, kv3+is-nspin,nspin
             if(map_k(ik) /= myrank_k) cycle        ! MPI
             call calc_phase(ia,natm,pos,kgp,iba(ik),ngabc,kg1 &
                  & ,nbase(1,ik), nbmx, zfcos,zfsin) !-(b_Elec.) ->(zfcos,zfsin)
             iksnl = (ik-1)/nspin+1
             do lmt1 = 1, ilmt(it)
                lmtt1=lmtt(lmt1,it); il=ltp(lmt1,it)
                call drv_of_betar_dot_WFs(gr,gi) !(sumset)
                ! --> gr,gi : derivatives of fsr, fsi
             end do
             if(mdvdb == EXECUT) then
                call term_related_to_drv_of_flmt_VDB
                ! gr,gi,fsr,fsi,vlhxcq,dion,eko_l,q,occup_l --> fx,fy,fz
             else if(mdvdb == SKIP) then
                call term_related_to_drv_of_flmt_NC
                ! gr,gi,fsr,fsi,dion,occup_l --> fx,fy,fz
             end if
          end do
       end do
       call pucv2cart(rltv,fx,fy,fz)
       if(npes >1) then
          fxyz_mpi1(1) = fx; fxyz_mpi1(2) = fy; fxyz_mpi1(3) = fz  ! MPI
          call mpi_allreduce(fxyz_mpi1,fxyz_mpi2,3,mpi_double_precision &
               & ,mpi_sum,MPI_CommGroup,ierr)          ! MPI
          fx = fxyz_mpi2(1); fy = fxyz_mpi2(2); fz = fxyz_mpi2(3)  ! MPI
       end if
       fnlxyz_l(ia,1) = fx; fnlxyz_l(ia,2) = fy; fnlxyz_l(ia,3) = fz
    end do
    if(npes>1) deallocate(fxyz_mpi2,fxyz_mpi1)
    deallocate(zfcos,zfsin)
    deallocate(ai,ar)
#else
    ncache = (m_CtrlP_cachesize()*1024)*3/4

    fnlxyz_l = 0.d0
    do is = 1, nspin, af+1
       do ik = is, kv3+is-nspin,nspin
          if(ncache == 0) then
             ibsize = iba(ik)
          else
             if(k_symmetry(ik) == GAMMA) then ! core2
                if(kimg == 1) then ! ar(i):1,zaj_l(i,ib,ik,1):np_e, ngabc(i,1:3):3
                   ibsize=ncache/(8*(np_e+1) + 4*3)
                else !  ar(i),ai(i):2,zaj_l(i,ib,ik,1),zaj_l(i,ib,ik,2):np_e*2,ngabc(i,1:3):3
                   ibsize=ncache/(8*(2+np_e*2) + 4*3)
                endif
             else
                if(kimg==1) then ! ar(i),ai(i):2,zaj_l(i,ib,ik,1):np_e,ngabc(i,1:3):3
                   ibsize = ncache/(8*(2+np_e) + 4*3)
                else if(kimg == 2) then ! ar(i),ai(i):2,zaj_l(i,ib,ik,1:2):np_e*2,ngabc(i,1:3):3
                   ibsize = ncache/(8*(2+np_e*2)+4*3)
                end if
             end if
          end if
          if(ipriforce>= 2) write(nfout,'(" ! ibsize = ",i8," <<m_Force_term_drv_of_flmt>>")') ibsize
!!$          allocate(zaj_red(ibsize,np_e,kimg))
!!$          allocate(ngabc_red(ibsize,3))
!!$          allocate(snl_red(ibsize,nlmtt))
!!$          allocate(zfcos_red(ibsize)); zfcos_red = 0.d0
!!$          allocate(zfsin_red(ibsize)); zfsin_red = 0.d0

          if(map_k(ik) /= myrank_k) cycle
          iksnl = (ik-1)/nspin+1
          IBLOCK: do ibl1 = 1, iba(ik), ibsize
             ibl2 = ibl1+ibsize-1
             if(ibl2 .gt. iba(ik)) ibl2 = iba(ik)
             if(ipriforce>= 2) write(nfout,'(" ! ibl1, ibl2 = ",2i8," <<m_Force_term_drv_of_flmt>>")') ibl1,ibl2

             allocate(ngabc_red(ibl1:ibl2,3)); ngabc_red = 0
             allocate(zfcos(ibl1:ibl2))
             allocate(zfsin(ibl1:ibl2))
             allocate(ar(ibl1:ibl2))        ;allocate(ai(ibl1:ibl2))
!!$             call substitute_zajred(ik)    ! zaj_l -> zaj_red
             call substitute_ngabcred(ik)  ! ngabc,nbase -> ngabc_red
!!$             call substitute_snlred(iksnl) ! snl -> snl_red

             do ia = 1, natm
                it = ityp(ia)
                call calc_phase_div(ia) ! pos, ngabc_red -> zfcos, zfsin

                do lmt1 = 1, ilmt(it)
                   lmtt1 = lmtt(lmt1,it); il = ltp(lmt1,it)
                   call drv_of_betar_dot_WFs_div(gr,gi)
                end do

                fx = 0.d0; fy = 0.d0; fz = 0.d0
                mdvdb = m_PP_include_vanderbilt_pot(it)
                if(mdvdb == EXECUT) then
                   call term_related_to_drv_of_flmt_VDB ! -> fx,fy,fz
                else if(mdvdb == SKIP) then
                   call term_related_to_drv_of_flmt_NC  ! -> fx,fy,fz
                end if
                call pucv2cart(rltv,fx,fy,fz)
                fnlxyz_l(ia,1) = fnlxyz_l(ia,1) + fx
                fnlxyz_l(ia,2) = fnlxyz_l(ia,2) + fy
                fnlxyz_l(ia,3) = fnlxyz_l(ia,3) + fz
             end do
          deallocate(ai,ar)
          deallocate(zfsin,zfcos)
          deallocate(ngabc_red)
          end do IBLOCK
!!$          deallocate(snl_red)
!!$          deallocate(zaj_red)
       end do
    end do
    if(npes>1) then
       allocate(fnlxyz_mpi(natm,3)); fnlxyz_mpi = 0.d0
!       call mpi_allreduce(fnlxyz_l,fnlxyz_mpi,3*natm,mpi_double_precision &
!            & , mpi_sum, MPI_CommGroup,ierr)
       call mpi_allreduce(fnlxyz_l,fnlxyz_mpi,3*natm,mpi_double_precision &
            & , mpi_sum, mpi_spin_group,ierr)
       fnlxyz_l = fnlxyz_mpi
    end if
    if(npes>1)deallocate(fnlxyz_mpi)
#endif


    if(af /= 0) then
       allocate(forc_up(natm2,3));allocate(forc_down(natm2,3))
       call forcaf2(natm,natm2,nopr+af,natm,nopr,napt &
         & ,iwei,op,ipriforce,af,fnlxyz_l ,forc_up,forc_down)
       deallocate(forc_down);deallocate(forc_up)
    end if

    deallocate(gi); deallocate(gr)

    if(ipriforce >= 2)  then
       write(nfout,'(" --- fnlxyz_l ->>")')
       do ia = 1, natm
          write(nfout,'(i8,3f20.8)') ia, fnlxyz_l(ia,1:3)
       end do
       write(nfout,'(" <<- fnlxyz_l ---")')
    end if

    call tstatc0_end(id_sname)
  contains

    subroutine term_related_to_drv_of_flmt_NC
      integer :: lmt1, lmta1, lmt2, lmta2, l1,m1,l2,m2, i
      real(kind=DP) :: fac_n, a, a_n, fr1,fi1,fr2,fi2
      do lmt1 = 1, ilmt(it)
         lmta1 = lmta(lmt1,ia); l1 = ltp(lmt1,it); m1 = mtp(lmt1,it)
         do lmt2 = lmt1, ilmt(it)
            l2 = ltp(lmt2,it); m2 = mtp(lmt2,it); lmta2 = lmta(lmt2,ia)

            if ( ipaw(it) == 0 ) then
               if(l1 /= l2 .or. m1 /= m2) cycle
               a = dion(lmt1,lmt2,it)
            else
               a = dion_paw(lmt1,lmt2,is,ia)
            endif
            fac_n = 2*fac; if(lmt1 == lmt2) fac_n = fac

            if(k_symmetry(ik) == GAMMA) then
               do i = 1, np_e                                 ! MPI
                  a_n = fac_n * a *occup_l(i,ik)
                  fr1 = fsr_l(i,lmta1,ik); fr2 = fsr_l(i,lmta2,ik)
                  fx  = fx - a_n * (gr(i,lmt1,1)*fr2 + fr1*gr(i,lmt2,1))
                  fy  = fy - a_n * (gr(i,lmt1,2)*fr2 + fr1*gr(i,lmt2,2))
                  fz  = fz - a_n * (gr(i,lmt1,3)*fr2 + fr1*gr(i,lmt2,3))
               end do
            else
               do i = 1, np_e                                 ! MPI
                  a_n = fac_n * a *occup_l(i,ik)
                  fr1 = fsr_l(i,lmta1,ik); fr2 = fsr_l(i,lmta2,ik)
                  fi1 = fsi_l(i,lmta1,ik); fi2 = fsi_l(i,lmta2,ik)
                  fx  = fx - a_n * (gr(i,lmt1,1)*fr2+gi(i,lmt1,1)*fi2&
                       &           +fr1*gr(i,lmt2,1)+fi1*gi(i,lmt2,1))
                  fy  = fy - a_n * (gr(i,lmt1,2)*fr2+gi(i,lmt1,2)*fi2&
                       &           +fr1*gr(i,lmt2,2)+fi1*gi(i,lmt2,2))
                  fz  = fz - a_n * (gr(i,lmt1,3)*fr2+gi(i,lmt1,3)*fi2&
                       &           +fr1*gr(i,lmt2,3)+fi1*gi(i,lmt2,3))
               end do
            end if
         end do
      end do
    end subroutine term_related_to_drv_of_flmt_NC

    subroutine term_related_to_drv_of_flmt_VDB
      integer :: lmt1, lmta1, lmt2, lmta2, i
      real(kind=DP) :: fac_n, a, a_n, fr1,fi1,fr2,fi2
      integer       :: id_sname = -1
      call tstatc0_begin('term_related_to_drv_of_flmt_VDB ',id_sname)

      do lmt1 = 1, ilmt(it)
         lmta1 = lmta(lmt1,ia)
         do lmt2 = lmt1, ilmt(it)
            lmta2 = lmta(lmt2,ia)
            fac_n = 2*fac; if(lmt1 == lmt2) fac_n = fac
!!$            a = vlhxcq(lmt1,lmt2,ia,is)+dion(lmt1,lmt2,it)
            if(ipaw(it)==0) then
                a = vlhxcq(lmt1,lmt2,ia,is)+dion(lmt1,lmt2,it)
            else
                a = vlhxcq(lmt1,lmt2,ia,is)+dion_paw(lmt1,lmt2,is,ia)
            end if
            if(k_symmetry(ik) == GAMMA) then
               do i = 1, np_e                                 ! MPI
                  a_n = fac_n*(a - eko_l(i,ik)*q(lmt1,lmt2,it))*occup_l(i,ik)
                  fr1 = fsr_l(i,lmta1,ik); fr2 = fsr_l(i,lmta2,ik)
                  fx  = fx - a_n * (gr(i,lmt1,1)*fr2 + fr1*gr(i,lmt2,1))
                  fy  = fy - a_n * (gr(i,lmt1,2)*fr2 + fr1*gr(i,lmt2,2))
                  fz  = fz - a_n * (gr(i,lmt1,3)*fr2 + fr1*gr(i,lmt2,3))
               end do
            else
               do i = 1, np_e                                 ! MPI
                  a_n = fac_n*(a - eko_l(i,ik)*q(lmt1,lmt2,it))*occup_l(i,ik)
                  fr1 = fsr_l(i,lmta1,ik); fr2 = fsr_l(i,lmta2,ik)
                  fi1 = fsi_l(i,lmta1,ik); fi2 = fsi_l(i,lmta2,ik)
                  fx  = fx - a_n * (gr(i,lmt1,1)*fr2+gi(i,lmt1,1)*fi2&
                       &           +fr1*gr(i,lmt2,1)+fi1*gi(i,lmt2,1))
                  fy  = fy - a_n * (gr(i,lmt1,2)*fr2+gi(i,lmt1,2)*fi2&
                       &           +fr1*gr(i,lmt2,2)+fi1*gi(i,lmt2,2))
                  fz  = fz - a_n * (gr(i,lmt1,3)*fr2+gi(i,lmt1,3)*fi2&
                       &           +fr1*gr(i,lmt2,3)+fi1*gi(i,lmt2,3))
               end do
            end if
         end do
      end do
      call tstatc0_end(id_sname)
    end subroutine term_related_to_drv_of_flmt_VDB

#ifdef SX
    subroutine drv_of_betar_dot_WFs(gr,gi)
    real(kind=DP),intent(inout), dimension(np_e,nlmta,3) :: gr,gi
      integer :: ib, i, i1, id3, mil
      real(kind=DP) :: br, bi, tmp, fac
      integer :: id_sname = -1
      call tstatc0_begin('drv_of_betar_dot_WFs ',id_sname)
#ifdef HIUX
*poption parallel
#endif
      do i = 1, iba(ik)
!!$         i1 = nbase(i,ik)
!!$         ar(i) = zfcos(i1)*snl(i,lmtt1,iksnl)
!!$         ai(i) = zfsin(i1)*snl(i,lmtt1,iksnl)
         ar(i) = zfcos(i)*snl(i,lmtt1,iksnl)
         ai(i) = zfsin(i)*snl(i,lmtt1,iksnl)
      end do

      if(k_symmetry(ik) == GAMMA) then
         if(kimg == 1) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
#endif
            do ib = 1, np_e                                    ! MPI
               gr(ib,lmt1,1:3) = 0.d0
               gi(ib,lmt1,1:3) = 0.d0
               do i = 1, iba(ik)
                  i1 = nbase(i,ik)
                  br = ar(i)*zaj_l(i,ib,ik,1)
                  gr(ib,lmt1,1) = gr(ib,lmt1,1) + br*ngabc(i1,1)
                  gr(ib,lmt1,2) = gr(ib,lmt1,2) + br*ngabc(i1,2)
                  gr(ib,lmt1,3) = gr(ib,lmt1,3) + br*ngabc(i1,3)
!!$                  gi(ib,lmt1,1) = gi(ib,lmt1,1) + bi*ngabc(i1,1)
!!$                  gi(ib,lmt1,2) = gi(ib,lmt1,2) + bi*ngabc(i1,2)
!!$                  gi(ib,lmt1,3) = gi(ib,lmt1,3) + bi*ngabc(i1,3)
               end do
            end do
         else if(kimg == 2) then
            mil = mod(il,4)
            if(mil == 1 .or. mil == 3) then ! l: even
               if(mil == 1) fac = -2.d0; if(mil == 3) fac = 2.d0
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
#endif
#ifdef NEC_TUNE3
!CDIR OUTERUNROLL
!CDIR NOLOOPCHG
#endif
!CDIR OUTERUNROLL=4
!CDIR NOLOOPCHG
               do ib = 1, np_e                                    ! MPI
                  gr(ib,lmt1,1:3) = 0.d0
                  do i = 2, iba(ik)
                     i1 = nbase(i,ik)
                     bi = (ai(i)*zaj_l(i,ib,ik,1) + ar(i)*zaj_l(i,ib,ik,2))*fac
                     gr(ib,lmt1,1) = gr(ib,lmt1,1) + bi*ngabc(i1,1)
                     gr(ib,lmt1,2) = gr(ib,lmt1,2) + bi*ngabc(i1,2)
                     gr(ib,lmt1,3) = gr(ib,lmt1,3) + bi*ngabc(i1,3)
                  end do
               end do
            else if(mil == 0 .or. mil == 2) then ! l: odd
               if(mil == 0) fac = 2.d0; if(mil == 2) fac = -2.d0;
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
#endif
#ifdef NEC_TUNE3
!CDIR OUTERUNROLL
!CDIR NOLOOPCHG
#endif
               do ib = 1, np_e
                  gr(ib,lmt1,1:3) = 0.d0
                  do i = 2, iba(ik)
                     i1 = nbase(i,ik)
                     br = (ar(i)*zaj_l(i,ib,ik,1) - ai(i)*zaj_l(i,ib,ik,2))*fac
                     gr(ib,lmt1,1) = gr(ib,lmt1,1) + br*ngabc(i1,1)
                     gr(ib,lmt1,2) = gr(ib,lmt1,2) + br*ngabc(i1,2)
                     gr(ib,lmt1,3) = gr(ib,lmt1,3) + br*ngabc(i1,3)
                  end do
               end do
            end if
         end if
      else
         if(kimg == 1) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
#endif
#ifdef NEC_TUNE3
!CDIR OUTERUNROLL
!CDIR NOLOOPCHG
#endif
            do ib = 1, np_e                                    ! MPI
               gr(ib,lmt1,1:3) = 0.d0
               gi(ib,lmt1,1:3) = 0.d0
               do i = 1, iba(ik)
                  i1 = nbase(i,ik)
                  br = ar(i)*zaj_l(i,ib,ik,1); bi = ai(i)*zaj_l(i,ib,ik,1)
                  gr(ib,lmt1,1) = gr(ib,lmt1,1) + br*ngabc(i1,1)
                  gr(ib,lmt1,2) = gr(ib,lmt1,2) + br*ngabc(i1,2)
                  gr(ib,lmt1,3) = gr(ib,lmt1,3) + br*ngabc(i1,3)
                  gi(ib,lmt1,1) = gi(ib,lmt1,1) + bi*ngabc(i1,1)
                  gi(ib,lmt1,2) = gi(ib,lmt1,2) + bi*ngabc(i1,2)
                  gi(ib,lmt1,3) = gi(ib,lmt1,3) + bi*ngabc(i1,3)
               end do
            end do
         else if(kimg == 2) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
#endif
#ifdef NEC_TUNE3
!CDIR OUTERUNROLL
!CDIR NOLOOPCHG
#endif
            do ib = 1, np_e                                    ! MPI
               gr(ib,lmt1,1:3) = 0.d0
               gi(ib,lmt1,1:3) = 0.d0
               do i = 1, iba(ik)
                  i1 = nbase(i,ik)
                  br = ar(i)*zaj_l(i,ib,ik,1) - ai(i)*zaj_l(i,ib,ik,2)
                  bi = ai(i)*zaj_l(i,ib,ik,1) + ar(i)*zaj_l(i,ib,ik,2)
                  gr(ib,lmt1,1) = gr(ib,lmt1,1) + br*ngabc(i1,1)
                  gr(ib,lmt1,2) = gr(ib,lmt1,2) + br*ngabc(i1,2)
                  gr(ib,lmt1,3) = gr(ib,lmt1,3) + br*ngabc(i1,3)
                  gi(ib,lmt1,1) = gi(ib,lmt1,1) + bi*ngabc(i1,1)
                  gi(ib,lmt1,2) = gi(ib,lmt1,2) + bi*ngabc(i1,2)
                  gi(ib,lmt1,3) = gi(ib,lmt1,3) + bi*ngabc(i1,3)
               end do
            end do
         end if
         mil = mod(il,4)
         do id3 = 1, 3
            if(mil == 1 .or. mil == 3) then
               if(mil == 1) fac = 1.d0;  if(mil == 3) fac = -1.d0
#ifdef HIUX
*poption parallel
#endif
               do ib = 1, np_e                              ! MPI
                  tmp = gr(ib,lmt1,id3)
                  gr(ib,lmt1,id3) = -fac*gi(ib,lmt1,id3)
                  gi(ib,lmt1,id3) = fac*tmp
               end do
            else if(mil == 2) then
#ifdef HIUX
*poption parallel
#endif
               do ib = 1, np_e                       ! MPI
                  gr(ib,lmt1,id3) = - gr(ib,lmt1,id3)
                  gi(ib,lmt1,id3) = - gi(ib,lmt1,id3)
               end do
            end if
         end do
      end if

      call tstatc0_end(id_sname)
    end subroutine drv_of_betar_dot_WFs
#else
!!$    subroutine substitute_zajred(ik)
!!$      integer, intent(in) :: ik
!!$      integer :: ie, ig
!!$      if(ibl2-ibl1+1 > ibsize) stop ' ibl2-ibl1+1 > ibsize'
!!$      if(kimg==1) then
!!$         do ie = 1, np_e
!!$            do ig = 1, ibl2-ibl1+1
!!$               zaj_red(ig,ie,1) = zaj_l(ig+ibl1-1,ie,ik,1)
!!$            end do
!!$         end do
!!$      else if(kimg==2) then
!!$         do ie = 1, np_e
!!$            do ig = 1, ibl2-ibl1+1
!!$               zaj_red(ig,ie,1)    = zaj_l(ig+ibl1-1,ie,ik,1)
!!$               zaj_red(ig,ie,kimg) = zaj_l(ig+ibl1-1,ie,ik,kimg)
!!$            end do
!!$         end do
!!$      end if
!!$    end subroutine substitute_zajred

    subroutine substitute_ngabcred(ik)
      integer, intent(in) :: ik
      integer :: ig, i1
      if(ibl2-ibl1+1 > ibsize) call phase_error_with_msg(nfout,' ibl2-ibl1+1 > ibsize',__LINE__,__FILE__)
!!$      if(ipriforce >= 1) write(nfout,'(" ibl1, ibl2 = ",2i8)') ibl1, ibl2
      do ig = ibl1, ibl2
         i1 = nbase(ig,ik)
         ngabc_red(ig,1) = ngabc(i1,1)
         ngabc_red(ig,2) = ngabc(i1,2)
         ngabc_red(ig,3) = ngabc(i1,3)
      end do
!!$      do ig = 1, ibl2-ibl1+1
!!$         i1 = nbase(ig+ibl1-1,ik)
!!$         ngabc_red(ig,1) = ngabc(i1,1)
!!$         ngabc_red(ig,2) = ngabc(i1,2)
!!$         ngabc_red(ig,3) = ngabc(i1,3)
!!$      end do
    end subroutine substitute_ngabcred

!!$    subroutine substitute_snlred(iksnl)
!!$      integer, intent(in) :: iksnl
!!$      integer :: ig, lmt1
!!$      if(ibl2-ibl1+1 > ibsize) stop ' ibl2-ibl1+1 > ibsize'
!!$      do lmt1 = 1, nlmtt
!!$         do ig = 1, ibl2-ibl1+1
!!$            snl_red(ig,lmt1) = snl(ig+ibl1-1,lmt1,iksnl)
!!$         end do
!!$      end do
!!$    end subroutine substitute_snlred

    subroutine calc_phase_div(ia)
      integer, intent(in) :: ia
      integer :: ig
      real(kind=DP) :: fx,fy,fz, ph
      fx = pos(ia,1)*PAI2
      fy = pos(ia,2)*PAI2
      fz = pos(ia,3)*PAI2
      do ig = ibl1, ibl2
         ph = ngabc_red(ig,1)*fx + ngabc_red(ig,2)*fy + ngabc_red(ig,3)*fz
         zfcos(ig) = dcos(ph)
         zfsin(ig) = dsin(ph)
      end do
!!$      do ig = 1, ibl2-ibl1+1
!!$         ph = ngabc_red(ig,1)*fx + ngabc_red(ig,2)*fy + ngabc_red(ig,3)*fz
!!$         zfcos_red(ig) = dcos(ph)
!!$         zfsin_red(ig) = dsin(ph)
!!$      end do
    end subroutine calc_phase_div

    subroutine drv_of_betar_dot_WFs_div(gr,gi)
      real(kind=DP),intent(inout), dimension(np_e,nlmta,3) :: gr,gi
      integer :: ib, i, i1, id3, mil, ibl1_2
      real(kind=DP) :: br, bi, tmp, fac
      integer :: id_sname = -1
      call tstatc0_begin('drv_of_betar_dot_WFs_div ',id_sname)

      do i = ibl1, ibl2
         ar(i) = zfcos(i)*snl(i,lmtt1,iksnl)
         ai(i) = zfsin(i)*snl(i,lmtt1,iksnl)
      end do

      if(k_symmetry(ik) == GAMMA) then
         if(kimg == 1) then
            do ib = 1, np_e                                    ! MPI
               gr(ib,lmt1,1:3) = 0.d0
               gi(ib,lmt1,1:3) = 0.d0
               do i = ibl1, ibl2
                  br = ar(i)*zaj_l(i,ib,ik,1)
                  gr(ib,lmt1,1) = gr(ib,lmt1,1) + br*ngabc_red(i,1)
                  gr(ib,lmt1,2) = gr(ib,lmt1,2) + br*ngabc_red(i,2)
                  gr(ib,lmt1,3) = gr(ib,lmt1,3) + br*ngabc_red(i,3)
               end do
            end do
         else if(kimg == 2) then
            ibl1_2 = max(2,ibl1)
            mil = mod(il,4)
            if(mil == 1 .or. mil == 3) then ! l: even
               if(mil == 1) fac = -2.d0; if(mil == 3) fac = 2.d0
               do ib = 1, np_e                                    ! MPI
                  gr(ib,lmt1,1:3) = 0.d0
                  do i = ibl1_2, ibl2
                     bi = (ai(i)*zaj_l(i,ib,ik,1) + ar(i)*zaj_l(i,ib,ik,2))*fac
                     gr(ib,lmt1,1) = gr(ib,lmt1,1) + bi*ngabc_red(i,1)
                     gr(ib,lmt1,2) = gr(ib,lmt1,2) + bi*ngabc_red(i,2)
                     gr(ib,lmt1,3) = gr(ib,lmt1,3) + bi*ngabc_red(i,3)
                  end do
               end do
            else if(mil == 0 .or. mil == 2) then ! l: odd
               if(mil == 0) fac = 2.d0; if(mil == 2) fac = -2.d0;
               do ib = 1, np_e
                  gr(ib,lmt1,1:3) = 0.d0
                  do i = ibl1_2, ibl2
                     br = (ar(i)*zaj_l(i,ib,ik,1) - ai(i)*zaj_l(i,ib,ik,2))*fac
                     gr(ib,lmt1,1) = gr(ib,lmt1,1) + br*ngabc_red(i,1)
                     gr(ib,lmt1,2) = gr(ib,lmt1,2) + br*ngabc_red(i,2)
                     gr(ib,lmt1,3) = gr(ib,lmt1,3) + br*ngabc_red(i,3)
                  end do
               end do
            end if
         end if
      else
         if(kimg == 1) then
            do ib = 1, np_e                                    ! MPI
               gr(ib,lmt1,1:3) = 0.d0
               gi(ib,lmt1,1:3) = 0.d0
               do i = ibl1, ibl2
                  br = ar(i)*zaj_l(i,ib,ik,1); bi = ai(i)*zaj_l(i,ib,ik,1)
                  gr(ib,lmt1,1) = gr(ib,lmt1,1) + br*ngabc_red(i,1)
                  gr(ib,lmt1,2) = gr(ib,lmt1,2) + br*ngabc_red(i,2)
                  gr(ib,lmt1,3) = gr(ib,lmt1,3) + br*ngabc_red(i,3)
                  gi(ib,lmt1,1) = gi(ib,lmt1,1) + bi*ngabc_red(i,1)
                  gi(ib,lmt1,2) = gi(ib,lmt1,2) + bi*ngabc_red(i,2)
                  gi(ib,lmt1,3) = gi(ib,lmt1,3) + bi*ngabc_red(i,3)
               end do
            end do
         else if(kimg == 2) then
            do ib = 1, np_e                                    ! MPI
               gr(ib,lmt1,1:3) = 0.d0
               gi(ib,lmt1,1:3) = 0.d0
               do i = ibl1, ibl2
                  br = ar(i)*zaj_l(i,ib,ik,1) - ai(i)*zaj_l(i,ib,ik,2)
                  bi = ai(i)*zaj_l(i,ib,ik,1) + ar(i)*zaj_l(i,ib,ik,2)
                  gr(ib,lmt1,1) = gr(ib,lmt1,1) + br*ngabc_red(i,1)
                  gr(ib,lmt1,2) = gr(ib,lmt1,2) + br*ngabc_red(i,2)
                  gr(ib,lmt1,3) = gr(ib,lmt1,3) + br*ngabc_red(i,3)
                  gi(ib,lmt1,1) = gi(ib,lmt1,1) + bi*ngabc_red(i,1)
                  gi(ib,lmt1,2) = gi(ib,lmt1,2) + bi*ngabc_red(i,2)
                  gi(ib,lmt1,3) = gi(ib,lmt1,3) + bi*ngabc_red(i,3)
               end do
            end do
         end if
         mil = mod(il,4)
         do id3 = 1, 3
            if(mil == 1 .or. mil == 3) then
               if(mil == 1) fac = 1.d0;  if(mil == 3) fac = -1.d0
               do ib = 1, np_e                              ! MPI
                  tmp = gr(ib,lmt1,id3)
                  gr(ib,lmt1,id3) = -fac*gi(ib,lmt1,id3)
                  gi(ib,lmt1,id3) = fac*tmp
               end do
            else if(mil == 2) then
               do ib = 1, np_e                       ! MPI
                  gr(ib,lmt1,id3) = - gr(ib,lmt1,id3)
                  gi(ib,lmt1,id3) = - gi(ib,lmt1,id3)
               end do
            end if
         end do
      end if

      call tstatc0_end(id_sname)
    end subroutine drv_of_betar_dot_WFs_div
#endif
  end subroutine m_Force_term_drv_of_flmt

! ============================== added by K. Tagami ===================== 11.0
  subroutine m_Force_term_drv_of_flmt_noncl(kv3,pos,napt)
    integer, intent(in)                          :: kv3
    real(kind=DP), intent(in), dimension(natm,3) :: pos
    integer,intent(in), dimension(natm,nopr+af)     :: napt

    integer  :: ia,is,ik,iksnl,it,lmt1,lmtt1,il,mdvdb
    integer  :: is1
    real(kind=DP) :: fac, fx, fy, fz
#ifdef SX
    real(kind=DP),allocatable,dimension(:) :: fxyz_mpi1,fxyz_mpi2  ! MPI
    real(kind=DP),allocatable,dimension(:) :: zfsin, zfcos  ! MPI
#else
    real(kind=DP),allocatable,dimension(:) :: zfsin, zfcos  ! MPI
    integer :: ibl1,ibl2,ibsize,ncache
!!$    real(kind=DP),allocatable,dimension(:,:,:) :: zaj_red ! d(ibsize,np_e,kimg)
!!$    real(kind=DP),allocatable,dimension(:,:) :: fnlxyz_mpi, snl_red
    real(kind=DP),allocatable,dimension(:,:) :: fnlxyz_mpi
    real(kind=DP),allocatable,dimension(:) :: zfsin_red,zfcos_red  ! MPI
    integer, allocatable, dimension(:,:) :: ngabc_red
#endif

    real(kind=DP), allocatable :: gr_noncl( :,:,:,: )
    real(kind=DP), allocatable :: gi_noncl( :,:,:,: )

    integer       :: id_sname = -1
    call tstatc0_begin('m_Force_term_drv_of_flmt ',id_sname,1)

    allocate(gr_noncl(np_e,nlmta,3,ndim_spinor)); gr_noncl = 0.0d0
    allocate(gi_noncl(np_e,nlmta,3,ndim_spinor)); gi_noncl = 0.0d0

!!!!!!!!    fac = 2.d0 / (kv3/ndim_spinor )
    fac = 1.d0 / (kv3/ndim_spinor )

#ifdef SX
    allocate(ar(kg1))        ;allocate(ai(kg1))
    allocate(zfsin(nbmx)); zfsin = 0.d0
    allocate(zfcos(nbmx)); zfcos = 0.d0
    if(npes>1) then
       allocate(fxyz_mpi1(3));    allocate(fxyz_mpi2(3))      ! MPI
    end if

    do ia = 1, natm
       fx = 0.d0; fy = 0.d0; fz = 0.d0
!!$       call calc_phase2(natm,pos,ia,kgp,nbmx,ngabc,nbmx,zfcos,zfsin)
       it = ityp(ia)
       mdvdb = m_PP_include_vanderbilt_pot(it)

       do ik = 1, kv3, ndim_spinor
          if(map_k(ik) /= myrank_k) cycle        ! MPI
          call calc_phase(ia,natm,pos,kgp,iba(ik),ngabc,kg1 &
               & ,nbase(1,ik), nbmx, zfcos,zfsin) !-(b_Elec.) ->(zfcos,zfsin)
          iksnl = (ik-1)/nspin+1

          do lmt1 = 1, ilmt(it)
             lmtt1=lmtt(lmt1,it); il=ltp(lmt1,it)

             Do is1=1, ndim_spinor
                call drv_of_betar_dot_WFs( ik+is1-1, gr_noncl(:,:,:,is1), &
                     &                               gi_noncl(:,:,:,is1) )
                                       ! --> gr,gi : derivatives of fsr, fsi
             End do
          end do
          if(mdvdb == EXECUT) then
             call term_related_to_drv_of_flmt_VDB
                  ! gr,gi,fsr,fsi,vlhxcq,dion,eko_l,q,occup_l --> fx,fy,fz
          else if(mdvdb == SKIP) then
             call term_related_to_drv_of_flmt_NC
                   ! gr,gi,fsr,fsi,dion,occup_l --> fx,fy,fz
          end if
       end do

       call pucv2cart(rltv,fx,fy,fz)
       if(npes >1) then
          fxyz_mpi1(1) = fx; fxyz_mpi1(2) = fy; fxyz_mpi1(3) = fz  ! MPI
!          call mpi_allreduce(fxyz_mpi1,fxyz_mpi2,3,mpi_double_precision &
!               & ,mpi_sum,MPI_CommGroup,ierr)          ! MPI
          call mpi_allreduce(fxyz_mpi1,fxyz_mpi2,3,mpi_double_precision &
               & ,mpi_sum,mpi_spin_group,ierr)          ! MPI
          fx = fxyz_mpi2(1); fy = fxyz_mpi2(2); fz = fxyz_mpi2(3)  ! MPI
       end if
       fnlxyz_l(ia,1) = fx; fnlxyz_l(ia,2) = fy; fnlxyz_l(ia,3) = fz
    end do
    if(npes>1) deallocate(fxyz_mpi2,fxyz_mpi1)
    deallocate(zfcos,zfsin)
    deallocate(ai,ar)
#else
    ncache = (m_CtrlP_cachesize()*1024)*3/4

    fnlxyz_l = 0.d0

    do ik = 1, kv3, ndim_spinor
       if(ncache == 0) then
          ibsize = iba(ik)
       else
          if(k_symmetry(ik) == GAMMA) then ! core2
             if(kimg == 1) then ! ar(i):1,zaj_l(i,ib,ik,1):np_e, ngabc(i,1:3):3
                ibsize=ncache/(8*(np_e+1) + 4*3)
             else !  ar(i),ai(i):2,zaj_l(i,ib,ik,1),zaj_l(i,ib,ik,2):np_e*2,ngabc(i,1:3):3
                ibsize=ncache/(8*(2+np_e*2) + 4*3)
             endif
          else
             if(kimg==1) then ! ar(i),ai(i):2,zaj_l(i,ib,ik,1):np_e,ngabc(i,1:3):3
                ibsize = ncache/(8*(2+np_e) + 4*3)
             else if(kimg == 2) then ! ar(i),ai(i):2,zaj_l(i,ib,ik,1:2):np_e*2,ngabc(i,1:3):3
                ibsize = ncache/(8*(2+np_e*2)+4*3)
             end if
          end if
       end if

       if(ipriforce>= 2) write(nfout,'(" ! ibsize = ",i8," <<m_Force_term_drv_of_flmt>>")') ibsize
!!$          allocate(zaj_red(ibsize,np_e,kimg))
!!$          allocate(ngabc_red(ibsize,3))
!!$          allocate(snl_red(ibsize,nlmtt))
!!$          allocate(zfcos_red(ibsize)); zfcos_red = 0.d0
!!$          allocate(zfsin_red(ibsize)); zfsin_red = 0.d0

       if(map_k(ik) /= myrank_k) cycle
       iksnl = (ik-1)/nspin+1

       IBLOCK: do ibl1 = 1, iba(ik), ibsize
          ibl2 = ibl1+ibsize-1
          if(ibl2 .gt. iba(ik)) ibl2 = iba(ik)
          if(ipriforce>= 2) write(nfout,'(" ! ibl1, ibl2 = ",2i8," <<m_Force_term_drv_of_flmt>>")') ibl1,ibl2

          allocate(ngabc_red(ibl1:ibl2,3)); ngabc_red = 0
          allocate(zfcos(ibl1:ibl2))
          allocate(zfsin(ibl1:ibl2))
          allocate(ar(ibl1:ibl2))        ;allocate(ai(ibl1:ibl2))
!!$             call substitute_zajred(ik)    ! zaj_l -> zaj_red
          call substitute_ngabcred(ik)  ! ngabc,nbase -> ngabc_red
!!$             call substitute_snlred(iksnl) ! snl -> snl_red

          do ia = 1, natm
             it = ityp(ia)
             call calc_phase_div(ia) ! pos, ngabc_red -> zfcos, zfsin

             do lmt1 = 1, ilmt(it)
                lmtt1 = lmtt(lmt1,it); il = ltp(lmt1,it)

                Do is1=1, ndim_spinor
                   call drv_of_betar_dot_WFs_div( ik+is1-1, gr_noncl(:,:,:,is1), &
                        &                                   gi_noncl(:,:,:,is1) )
                End do

             end do

             fx = 0.d0; fy = 0.d0; fz = 0.d0
             mdvdb = m_PP_include_vanderbilt_pot(it)
             if(mdvdb == EXECUT) then
                call term_related_to_drv_of_flmt_VDB ! -> fx,fy,fz
             else if(mdvdb == SKIP) then
                call term_related_to_drv_of_flmt_NC  ! -> fx,fy,fz
             end if

             call pucv2cart(rltv,fx,fy,fz)
             fnlxyz_l(ia,1) = fnlxyz_l(ia,1) + fx
             fnlxyz_l(ia,2) = fnlxyz_l(ia,2) + fy
             fnlxyz_l(ia,3) = fnlxyz_l(ia,3) + fz
          end do
          deallocate(ai,ar)
          deallocate(zfsin,zfcos)
          deallocate(ngabc_red)
       end do IBLOCK
!!$          deallocate(snl_red)
!!$          deallocate(zaj_red)
    end do
    if(npes>1) then
       allocate(fnlxyz_mpi(natm,3)); fnlxyz_mpi = 0.d0
!       call mpi_allreduce(fnlxyz_l,fnlxyz_mpi,3*natm,mpi_double_precision &
!            & , mpi_sum, MPI_CommGroup,ierr)
       call mpi_allreduce(fnlxyz_l,fnlxyz_mpi,3*natm,mpi_double_precision &
            & , mpi_sum, mpi_spin_group,ierr)
       fnlxyz_l = fnlxyz_mpi
    end if
    if(npes>1)deallocate(fnlxyz_mpi)
#endif


    if(af /= 0) then
       allocate(forc_up(natm2,3));allocate(forc_down(natm2,3))
       call forcaf2(natm,natm2,nopr+af,natm,nopr,napt &
         & ,iwei,op,ipriforce,af,fnlxyz_l ,forc_up,forc_down)
       deallocate(forc_down);deallocate(forc_up)
    end if

    deallocate(gi); deallocate(gr)

    if(ipriforce >= 2)  then
       write(nfout,'(" --- fnlxyz_l ->>")')
       do ia = 1, natm
          write(nfout,'(i8,3f20.8)') ia, fnlxyz_l(ia,1:3)
       end do
       write(nfout,'(" <<- fnlxyz_l ---")')
    end if

    call tstatc0_end(id_sname)
  contains

    subroutine term_related_to_drv_of_flmt_NC
      integer :: lmt1, lmta1, lmt2, lmta2, l1,m1,l2,m2, i
      real(kind=DP) :: fac_n, a, a_n, fr1,fi1,fr2,fi2

      integer :: is1, is2, istmp
      complex(kind=CMPLDP) :: z1, z2, z3, zg1(3), zg2(3), zsum(3)

      do lmt1 = 1, ilmt(it)
         lmta1 = lmta(lmt1,ia); l1 = ltp(lmt1,it); m1 = mtp(lmt1,it)

! === KT_mod ================= 2013/10/31
!!!         do lmt2 = lmt1, ilmt(it)
         do lmt2 = 1, ilmt(it)
! ============================ 2013/10/31

            l2 = ltp(lmt2,it); m2 = mtp(lmt2,it); lmta2 = lmta(lmt2,ia)

!!!!            if(l1 /= l2 .or. m1 /= m2) cycle

! === KT_mod ================= 2013/10/31
!!!            fac_n = 2*fac; if(lmt1 == lmt2) fac_n = fac
            fac_n = fac
! ============================ 2013/10/31

            if(k_symmetry(ik) == GAMMA) then
               call phase_error_with_msg(nfout,'Not supported : Gamma symmetry in non-collinear system',__LINE__,__FILE__)
            else

               do i = 1, np_e                                 ! MPI
                  zsum = 0.0d0
                  Do is1=1, ndim_spinor
                     Do is2=1, ndim_spinor
                        istmp = ( is1 -1 )*ndim_spinor + is2

                        z1 = cmplx( fsr_l(i,lmta1,ik+is1-1), fsi_l(i,lmta1,ik+is1-1 ) )
                        z2 = cmplx( fsr_l(i,lmta2,ik+is2-1), fsi_l(i,lmta2,ik+is2-1 ) )
                        zg1(:) = cmplx( gr_noncl(i,lmt1,:,is1), gi_noncl(i,lmt1,:,is1) )
                        zg2(:) = cmplx( gr_noncl(i,lmt2,:,is2), gi_noncl(i,lmt2,:,is2) )

                        zsum(:) = zsum(:) + (conjg(z1)*zg2(:) + conjg(zg1(:))*z2 )&
                             &              *dion_scr_noncl( lmt1,lmt2,istmp,ia )
                     End do
                  End Do

                  zsum = zsum *fac_n *occup_l(i,ik)

                  fx  = fx - real( zsum(1) )
                  fy  = fy - real( zsum(2) )
                  fz  = fz - real( zsum(3) )

               end do
            end if
         end do
      end do
    end subroutine term_related_to_drv_of_flmt_NC

    subroutine term_related_to_drv_of_flmt_VDB
      integer :: lmt1, lmta1, lmt2, lmta2, i
      real(kind=DP) :: fac_n, a, a_n, fr1,fi1,fr2,fi2
      integer :: is1, is2, istmp
      complex(kind=CMPLDP) :: z1, z2, z3, zg1(3), zg2(3), zsum(3)

      integer       :: id_sname = -1
      call tstatc0_begin('term_related_to_drv_of_flmt_VDB ',id_sname)

      do lmt1 = 1, ilmt(it)
         lmta1 = lmta(lmt1,ia)

! ============ KT_mod ================= 2013/10/31
!!!         do lmt2 = lmt1, ilmt(it)
         do lmt2 = 1, ilmt(it)
! ===================================== 2013/10/31

            lmta2 = lmta(lmt2,ia)

! ============ KT_mod ================= 2013/10/31
!!!            fac_n = 2*fac; if(lmt1 == lmt2) fac_n = fac
            fac_n = fac
! ===================================== 2013/10/31

            if(k_symmetry(ik) == GAMMA) then
               call phase_error_with_msg(nfout,'Not supported : Gamma symmetry in non-collinear system',__LINE__,__FILE__)
            else

              do i = 1, np_e                                 ! MPI
                  zsum = 0.0d0
                  Do is1=1, ndim_spinor
                     Do is2=1, ndim_spinor
                        istmp = ( is1 -1 )*ndim_spinor + is2

                        z1 = cmplx( fsr_l(i,lmta1,ik+is1-1), fsi_l(i,lmta1,ik+is1-1 ) )
                        z2 = cmplx( fsr_l(i,lmta2,ik+is2-1), fsi_l(i,lmta2,ik+is2-1 ) )
                        zg1(:) = cmplx( gr_noncl(i,lmt1,:,is1), gi_noncl(i,lmt1,:,is1) )
                        zg2(:) = cmplx( gr_noncl(i,lmt2,:,is2), gi_noncl(i,lmt2,:,is2) )

                        z3 = dion_scr_noncl( lmt1,lmt2,istmp,ia ) &
                             &  -eko_l(i,ik) *q_noncl(lmt1,lmt2,istmp,it )

                        zsum(:) = zsum(:) + (conjg(z1)*zg2(:) + conjg(zg1(:))*z2 ) *z3
                     End do
                  End Do

                  zsum = zsum *fac_n *occup_l(i,ik)

                  fx  = fx - real( zsum(1) )
                  fy  = fy - real( zsum(2) )
                  fz  = fz - real( zsum(3) )
               end do

            end if
         end do
      end do
      call tstatc0_end(id_sname)
    end subroutine term_related_to_drv_of_flmt_VDB

#ifdef SX
    subroutine drv_of_betar_dot_WFs(ik,gr,gi)
      integer, intent(in) :: ik

      real(kind=DP),intent(inout), dimension(np_e,nlmta,3) :: gr,gi
      integer :: ib, i, i1, id3, mil
      real(kind=DP) :: br, bi, tmp, fac
      integer :: id_sname = -1
      call tstatc0_begin('drv_of_betar_dot_WFs ',id_sname)
#ifdef HIUX
*poption parallel
#endif
      do i = 1, iba(ik)
!!$         i1 = nbase(i,ik)
!!$         ar(i) = zfcos(i1)*snl(i,lmtt1,iksnl)
!!$         ai(i) = zfsin(i1)*snl(i,lmtt1,iksnl)
         ar(i) = zfcos(i)*snl(i,lmtt1,iksnl)
         ai(i) = zfsin(i)*snl(i,lmtt1,iksnl)
      end do

      if(k_symmetry(ik) == GAMMA) then
         if(kimg == 1) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
#endif
            do ib = 1, np_e                                    ! MPI
               gr(ib,lmt1,1:3) = 0.d0
               gi(ib,lmt1,1:3) = 0.d0
               do i = 1, iba(ik)
                  i1 = nbase(i,ik)
                  br = ar(i)*zaj_l(i,ib,ik,1)
                  gr(ib,lmt1,1) = gr(ib,lmt1,1) + br*ngabc(i1,1)
                  gr(ib,lmt1,2) = gr(ib,lmt1,2) + br*ngabc(i1,2)
                  gr(ib,lmt1,3) = gr(ib,lmt1,3) + br*ngabc(i1,3)
!!$                  gi(ib,lmt1,1) = gi(ib,lmt1,1) + bi*ngabc(i1,1)
!!$                  gi(ib,lmt1,2) = gi(ib,lmt1,2) + bi*ngabc(i1,2)
!!$                  gi(ib,lmt1,3) = gi(ib,lmt1,3) + bi*ngabc(i1,3)
               end do
            end do
         else if(kimg == 2) then
            mil = mod(il,4)
            if(mil == 1 .or. mil == 3) then ! l: even
               if(mil == 1) fac = -2.d0; if(mil == 3) fac = 2.d0
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
#endif
#ifdef NEC_TUNE3
!CDIR OUTERUNROLL
!CDIR NOLOOPCHG
#endif
!CDIR OUTERUNROLL=4
!CDIR NOLOOPCHG
               do ib = 1, np_e                                    ! MPI
                  gr(ib,lmt1,1:3) = 0.d0
                  do i = 2, iba(ik)
                     i1 = nbase(i,ik)
                     bi = (ai(i)*zaj_l(i,ib,ik,1) + ar(i)*zaj_l(i,ib,ik,2))*fac
                     gr(ib,lmt1,1) = gr(ib,lmt1,1) + bi*ngabc(i1,1)
                     gr(ib,lmt1,2) = gr(ib,lmt1,2) + bi*ngabc(i1,2)
                     gr(ib,lmt1,3) = gr(ib,lmt1,3) + bi*ngabc(i1,3)
                  end do
               end do
            else if(mil == 0 .or. mil == 2) then ! l: odd
               if(mil == 0) fac = 2.d0; if(mil == 2) fac = -2.d0;
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
#endif
#ifdef NEC_TUNE3
!CDIR OUTERUNROLL
!CDIR NOLOOPCHG
#endif
               do ib = 1, np_e
                  gr(ib,lmt1,1:3) = 0.d0
                  do i = 2, iba(ik)
                     i1 = nbase(i,ik)
                     br = (ar(i)*zaj_l(i,ib,ik,1) - ai(i)*zaj_l(i,ib,ik,2))*fac
                     gr(ib,lmt1,1) = gr(ib,lmt1,1) + br*ngabc(i1,1)
                     gr(ib,lmt1,2) = gr(ib,lmt1,2) + br*ngabc(i1,2)
                     gr(ib,lmt1,3) = gr(ib,lmt1,3) + br*ngabc(i1,3)
                  end do
               end do
            end if
         end if
      else
         if(kimg == 1) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
#endif
#ifdef NEC_TUNE3
!CDIR OUTERUNROLL
!CDIR NOLOOPCHG
#endif
            do ib = 1, np_e                                    ! MPI
               gr(ib,lmt1,1:3) = 0.d0
               gi(ib,lmt1,1:3) = 0.d0
               do i = 1, iba(ik)
                  i1 = nbase(i,ik)
                  br = ar(i)*zaj_l(i,ib,ik,1); bi = ai(i)*zaj_l(i,ib,ik,1)
                  gr(ib,lmt1,1) = gr(ib,lmt1,1) + br*ngabc(i1,1)
                  gr(ib,lmt1,2) = gr(ib,lmt1,2) + br*ngabc(i1,2)
                  gr(ib,lmt1,3) = gr(ib,lmt1,3) + br*ngabc(i1,3)
                  gi(ib,lmt1,1) = gi(ib,lmt1,1) + bi*ngabc(i1,1)
                  gi(ib,lmt1,2) = gi(ib,lmt1,2) + bi*ngabc(i1,2)
                  gi(ib,lmt1,3) = gi(ib,lmt1,3) + bi*ngabc(i1,3)
               end do
            end do
         else if(kimg == 2) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
#endif
#ifdef NEC_TUNE3
!CDIR OUTERUNROLL
!CDIR NOLOOPCHG
#endif
            do ib = 1, np_e                                    ! MPI
               gr(ib,lmt1,1:3) = 0.d0
               gi(ib,lmt1,1:3) = 0.d0
               do i = 1, iba(ik)
                  i1 = nbase(i,ik)
                  br = ar(i)*zaj_l(i,ib,ik,1) - ai(i)*zaj_l(i,ib,ik,2)
                  bi = ai(i)*zaj_l(i,ib,ik,1) + ar(i)*zaj_l(i,ib,ik,2)
                  gr(ib,lmt1,1) = gr(ib,lmt1,1) + br*ngabc(i1,1)
                  gr(ib,lmt1,2) = gr(ib,lmt1,2) + br*ngabc(i1,2)
                  gr(ib,lmt1,3) = gr(ib,lmt1,3) + br*ngabc(i1,3)
                  gi(ib,lmt1,1) = gi(ib,lmt1,1) + bi*ngabc(i1,1)
                  gi(ib,lmt1,2) = gi(ib,lmt1,2) + bi*ngabc(i1,2)
                  gi(ib,lmt1,3) = gi(ib,lmt1,3) + bi*ngabc(i1,3)
               end do
            end do
         end if
         mil = mod(il,4)
         do id3 = 1, 3
            if(mil == 1 .or. mil == 3) then
               if(mil == 1) fac = 1.d0;  if(mil == 3) fac = -1.d0
#ifdef HIUX
*poption parallel
#endif
               do ib = 1, np_e                              ! MPI
                  tmp = gr(ib,lmt1,id3)
                  gr(ib,lmt1,id3) = -fac*gi(ib,lmt1,id3)
                  gi(ib,lmt1,id3) = fac*tmp
               end do
            else if(mil == 2) then
#ifdef HIUX
*poption parallel
#endif
               do ib = 1, np_e                       ! MPI
                  gr(ib,lmt1,id3) = - gr(ib,lmt1,id3)
                  gi(ib,lmt1,id3) = - gi(ib,lmt1,id3)
               end do
            end if
         end do
      end if

      call tstatc0_end(id_sname)
    end subroutine drv_of_betar_dot_WFs
#else
!!$    subroutine substitute_zajred(ik)
!!$      integer, intent(in) :: ik
!!$      integer :: ie, ig
!!$      if(ibl2-ibl1+1 > ibsize) stop ' ibl2-ibl1+1 > ibsize'
!!$      if(kimg==1) then
!!$         do ie = 1, np_e
!!$            do ig = 1, ibl2-ibl1+1
!!$               zaj_red(ig,ie,1) = zaj_l(ig+ibl1-1,ie,ik,1)
!!$            end do
!!$         end do
!!$      else if(kimg==2) then
!!$         do ie = 1, np_e
!!$            do ig = 1, ibl2-ibl1+1
!!$               zaj_red(ig,ie,1)    = zaj_l(ig+ibl1-1,ie,ik,1)
!!$               zaj_red(ig,ie,kimg) = zaj_l(ig+ibl1-1,ie,ik,kimg)
!!$            end do
!!$         end do
!!$      end if
!!$    end subroutine substitute_zajred

    subroutine substitute_ngabcred(ik)
      integer, intent(in) :: ik
      integer :: ig, i1
      if(ibl2-ibl1+1 > ibsize) call phase_error_with_msg(nfout,' ibl2-ibl1+1 > ibsize',__LINE__,__FILE__)
!!$      if(ipriforce >= 1) write(nfout,'(" ibl1, ibl2 = ",2i8)') ibl1, ibl2
      do ig = ibl1, ibl2
         i1 = nbase(ig,ik)
         ngabc_red(ig,1) = ngabc(i1,1)
         ngabc_red(ig,2) = ngabc(i1,2)
         ngabc_red(ig,3) = ngabc(i1,3)
      end do
!!$      do ig = 1, ibl2-ibl1+1
!!$         i1 = nbase(ig+ibl1-1,ik)
!!$         ngabc_red(ig,1) = ngabc(i1,1)
!!$         ngabc_red(ig,2) = ngabc(i1,2)
!!$         ngabc_red(ig,3) = ngabc(i1,3)
!!$      end do
    end subroutine substitute_ngabcred

!!$    subroutine substitute_snlred(iksnl)
!!$      integer, intent(in) :: iksnl
!!$      integer :: ig, lmt1
!!$      if(ibl2-ibl1+1 > ibsize) stop ' ibl2-ibl1+1 > ibsize'
!!$      do lmt1 = 1, nlmtt
!!$         do ig = 1, ibl2-ibl1+1
!!$            snl_red(ig,lmt1) = snl(ig+ibl1-1,lmt1,iksnl)
!!$         end do
!!$      end do
!!$    end subroutine substitute_snlred

    subroutine calc_phase_div(ia)
      integer, intent(in) :: ia
      integer :: ig
      real(kind=DP) :: fx,fy,fz, ph
      fx = pos(ia,1)*PAI2
      fy = pos(ia,2)*PAI2
      fz = pos(ia,3)*PAI2
      do ig = ibl1, ibl2
         ph = ngabc_red(ig,1)*fx + ngabc_red(ig,2)*fy + ngabc_red(ig,3)*fz
         zfcos(ig) = dcos(ph)
         zfsin(ig) = dsin(ph)
      end do
!!$      do ig = 1, ibl2-ibl1+1
!!$         ph = ngabc_red(ig,1)*fx + ngabc_red(ig,2)*fy + ngabc_red(ig,3)*fz
!!$         zfcos_red(ig) = dcos(ph)
!!$         zfsin_red(ig) = dsin(ph)
!!$      end do
    end subroutine calc_phase_div

    subroutine drv_of_betar_dot_WFs_div(ik,gr,gi)
      integer, intent(in) :: ik

      real(kind=DP),intent(inout), dimension(np_e,nlmta,3) :: gr,gi
      integer :: ib, i, i1, id3, mil, ibl1_2
      real(kind=DP) :: br, bi, tmp, fac
      integer :: id_sname = -1
      call tstatc0_begin('drv_of_betar_dot_WFs_div ',id_sname)

      do i = ibl1, ibl2
         ar(i) = zfcos(i)*snl(i,lmtt1,iksnl)
         ai(i) = zfsin(i)*snl(i,lmtt1,iksnl)
      end do

      if(k_symmetry(ik) == GAMMA) then
         if(kimg == 1) then
            do ib = 1, np_e                                    ! MPI
               gr(ib,lmt1,1:3) = 0.d0
               gi(ib,lmt1,1:3) = 0.d0
               do i = ibl1, ibl2
                  br = ar(i)*zaj_l(i,ib,ik,1)
                  gr(ib,lmt1,1) = gr(ib,lmt1,1) + br*ngabc_red(i,1)
                  gr(ib,lmt1,2) = gr(ib,lmt1,2) + br*ngabc_red(i,2)
                  gr(ib,lmt1,3) = gr(ib,lmt1,3) + br*ngabc_red(i,3)
               end do
            end do
         else if(kimg == 2) then
            ibl1_2 = max(2,ibl1)
            mil = mod(il,4)
            if(mil == 1 .or. mil == 3) then ! l: even
               if(mil == 1) fac = -2.d0; if(mil == 3) fac = 2.d0
               do ib = 1, np_e                                    ! MPI
                  gr(ib,lmt1,1:3) = 0.d0
                  do i = ibl1_2, ibl2
                     bi = (ai(i)*zaj_l(i,ib,ik,1) + ar(i)*zaj_l(i,ib,ik,2))*fac
                     gr(ib,lmt1,1) = gr(ib,lmt1,1) + bi*ngabc_red(i,1)
                     gr(ib,lmt1,2) = gr(ib,lmt1,2) + bi*ngabc_red(i,2)
                     gr(ib,lmt1,3) = gr(ib,lmt1,3) + bi*ngabc_red(i,3)
                  end do
               end do
            else if(mil == 0 .or. mil == 2) then ! l: odd
               if(mil == 0) fac = 2.d0; if(mil == 2) fac = -2.d0;
               do ib = 1, np_e
                  gr(ib,lmt1,1:3) = 0.d0
                  do i = ibl1_2, ibl2
                     br = (ar(i)*zaj_l(i,ib,ik,1) - ai(i)*zaj_l(i,ib,ik,2))*fac
                     gr(ib,lmt1,1) = gr(ib,lmt1,1) + br*ngabc_red(i,1)
                     gr(ib,lmt1,2) = gr(ib,lmt1,2) + br*ngabc_red(i,2)
                     gr(ib,lmt1,3) = gr(ib,lmt1,3) + br*ngabc_red(i,3)
                  end do
               end do
            end if
         end if
      else
         if(kimg == 1) then
            do ib = 1, np_e                                    ! MPI
               gr(ib,lmt1,1:3) = 0.d0
               gi(ib,lmt1,1:3) = 0.d0
               do i = ibl1, ibl2
                  br = ar(i)*zaj_l(i,ib,ik,1); bi = ai(i)*zaj_l(i,ib,ik,1)
                  gr(ib,lmt1,1) = gr(ib,lmt1,1) + br*ngabc_red(i,1)
                  gr(ib,lmt1,2) = gr(ib,lmt1,2) + br*ngabc_red(i,2)
                  gr(ib,lmt1,3) = gr(ib,lmt1,3) + br*ngabc_red(i,3)
                  gi(ib,lmt1,1) = gi(ib,lmt1,1) + bi*ngabc_red(i,1)
                  gi(ib,lmt1,2) = gi(ib,lmt1,2) + bi*ngabc_red(i,2)
                  gi(ib,lmt1,3) = gi(ib,lmt1,3) + bi*ngabc_red(i,3)
               end do
            end do
         else if(kimg == 2) then
            do ib = 1, np_e                                    ! MPI
               gr(ib,lmt1,1:3) = 0.d0
               gi(ib,lmt1,1:3) = 0.d0
               do i = ibl1, ibl2
                  br = ar(i)*zaj_l(i,ib,ik,1) - ai(i)*zaj_l(i,ib,ik,2)
                  bi = ai(i)*zaj_l(i,ib,ik,1) + ar(i)*zaj_l(i,ib,ik,2)
                  gr(ib,lmt1,1) = gr(ib,lmt1,1) + br*ngabc_red(i,1)
                  gr(ib,lmt1,2) = gr(ib,lmt1,2) + br*ngabc_red(i,2)
                  gr(ib,lmt1,3) = gr(ib,lmt1,3) + br*ngabc_red(i,3)
                  gi(ib,lmt1,1) = gi(ib,lmt1,1) + bi*ngabc_red(i,1)
                  gi(ib,lmt1,2) = gi(ib,lmt1,2) + bi*ngabc_red(i,2)
                  gi(ib,lmt1,3) = gi(ib,lmt1,3) + bi*ngabc_red(i,3)
               end do
            end do
         end if
         mil = mod(il,4)
         do id3 = 1, 3
            if(mil == 1 .or. mil == 3) then
               if(mil == 1) fac = 1.d0;  if(mil == 3) fac = -1.d0
               do ib = 1, np_e                              ! MPI
                  tmp = gr(ib,lmt1,id3)
                  gr(ib,lmt1,id3) = -fac*gi(ib,lmt1,id3)
                  gi(ib,lmt1,id3) = fac*tmp
               end do
            else if(mil == 2) then
               do ib = 1, np_e                       ! MPI
                  gr(ib,lmt1,id3) = - gr(ib,lmt1,id3)
                  gi(ib,lmt1,id3) = - gi(ib,lmt1,id3)
               end do
            end if
         end do
      end if

      call tstatc0_end(id_sname)
    end subroutine drv_of_betar_dot_WFs_div
#endif
  end subroutine m_Force_term_drv_of_flmt_noncl
! ========================================================================= 11.0

#endif

  subroutine m_Force_term_dipole()
    integer :: i,ia
    real(kind=DP) :: ed(3),ef(3),a2
#ifdef __TIMER_SUB__
    call timer_sta(1505)
#endif
    ed = 0.d0
    do i=1,3
       a2 = dot_product(altv(1:3,i),altv(1:3,i))
       ed(1:3) = ed(1:3)-PAI4*dipole(i)*altv(1:3,i)/a2
    end do
    ef = 0.d0
    do i=1,3
       ef(1:3) = ef(1:3)-2.d0*field(i)*rltv(1:3,i)
    end do
    do ia=1,natm
       fdip_l(ia,1:3) = ival(ityp(ia))*ed(1:3)
       fext_l(ia,1:3) = ival(ityp(ia))*ef(1:3)
    end do
#ifdef __TIMER_SUB__
    call timer_end(1505)
#endif
  end subroutine m_Force_term_dipole

  subroutine m_Force_term_fef()
    integer :: ia
    do ia=1,natm
       ffef_l(ia,1:3) = ival(ityp(ia))*elec_field(1:3)
    end do
  end subroutine m_Force_term_fef


!!$  subroutine m_Force_alloc_zfsin_zfcos
!!$    allocate(zfcos(ista_kngp:iend_kngp))
!!$    allocate(zfsin(ista_kngp:iend_kngp))
!!$  end subroutine m_Force_alloc_zfsin_zfcos

!!$  subroutine m_Force_dealloc_zfsin_zfcos
!!$    deallocate(zfsin); deallocate(zfcos)
!!$  end subroutine m_Force_dealloc_zfsin_zfcos

  subroutine m_Force_cal_forcmx()
    integer       :: ia,j
    real(kind=DP) :: fca
    real(kind=DP) :: fcatmp
    real(kind=DP), dimension(3) :: frc
#ifdef __TIMER_SUB__
    call timer_sta(1507)
#endif
    forcmx = 0.d0
    if(nfcatm <= 0) then
!xocl spread do/ind_natm
!!       do ia = 1, natm
!!          if(imdtyp(ia) == RELAX .or. imdtyp(ia) > HEAT_BATH) then
!!             fca = dsqrt(forc_l(ia,1)**2+forc_l(ia,2)**2+forc_l(ia,3)**2)
!!             if(forcmx < fca) forcmx = fca
!!          end if
!!       end do
       do ia = 1, natm
          fcatmp = 0.d0
          if(imdtypxyz(ia,1) ==0 .and. imdtypxyz(ia,2)==0 .and. imdtypxyz(ia,3)==0) cycle
          frc(:) = forc_l(ia,:)
          if(sw_fix_bond==ON) frc(:) = frc(:)+bond_forc(ia,:)
          do j=1,3
             if(imdtypxyz(ia,j) == RELAX .or. imdtypxyz(ia,j) > HEAT_BATH) fcatmp = fcatmp + frc(j)**2
          enddo
          fca = dsqrt(fcatmp)
          if(forcmx < fca) forcmx = fca
       end do
!xocl end spread sum(forcmx)
    else if(nfcatm > 0) then
!xocl spread do/ind_natm
       do ia = 1, natm
          if(imdtyp(ia) == RELAX .or. imdtyp(ia) > HEAT_BATH) then
             fca = dsqrt(forc_l(ia,1)**2+forc_l(ia,2)**2+forc_l(ia,3)**2)
             if(forcmx < fca) forcmx = fca
          else if(imdtyp(ia) == BONDLENGTH_FIX .or. imdtyp(ia) == BONDLENGTH_FIX_1 &
               & .or. imdtyp(ia) == BONDLENGTH_FIX_2 .or. imdtyp(ia) == COG_FIX) then
             fca = dsqrt(forc_l(ia,1)**2+forc_l(ia,2)**2+forc_l(ia,3)**2)
             if(forcmx < fca) forcmx = fca
          end if
       end do
!xocl end spread sum(forcmx)
    end if
    if(ipriforce >= 1) then
       if(iteration_ionic_wd_forcmx /= iteration_ionic) &
            & write(nfout,'(" !D forcmx = ",d20.12, " at the iteration_ionic of ",i8)') &
            & forcmx, iteration_ionic
       call flush(nfout)
       iteration_ionic_wd_forcmx = iteration_ionic
    end if
#ifdef __TIMER_SUB__
    call timer_end(1507)
#endif
  end subroutine m_Force_cal_forcmx

  logical function m_Forces_are_Converged_core()
    call m_Force_cal_forcmx()
    if(forcmx < forccr) then
       m_Forces_are_Converged_core = .true.
    else
       m_Forces_are_Converged_core = .false.
    end if
  end function m_Forces_are_Converged_core

  subroutine extract_largerforce_atoms(nforc, ndim, iatom_larger)
    integer, intent(in) :: nforc, ndim
    integer, intent(out), dimension(ndim) :: iatom_larger
    integer :: ia, ia2, ic,ic2, it
    real(kind=DP) :: sqr_forc, sqrforc_lower, tmp
    real(kind=DP), allocatable, dimension(:) :: dforce

    allocate(dforce(natm)); dforce = 0.d0
    do ia = 1, natm
       dforce(ia) = forc_l(ia,1)*forc_l(ia,1)+forc_l(ia,2)*forc_l(ia,2)+forc_l(ia,3)*forc_l(ia,3)
    end do
    do ia = 1, natm-1
       do ia2 = ia+1, natm
          if(dforce(ia) < dforce(ia2)) then
             sqr_forc = dforce(ia)
             dforce(ia) = dforce(ia2)
             dforce(ia2) = sqr_forc
          end if
       end do
    end do
    sqrforc_lower = dforce(min(nforc,natm))
    deallocate(dforce)
    allocate(dforce(ndim)); dforce = 0.d0
    ic = 0
    do ia = 1, natm
       sqr_forc = forc_l(ia,1)*forc_l(ia,1)+forc_l(ia,2)*forc_l(ia,2)+forc_l(ia,3)*forc_l(ia,3)
       if(sqr_forc >= sqrforc_lower) then
          ic = ic + 1
          if(ic > ndim) exit
          iatom_larger(ic) = ia
          dforce(ic) = sqr_forc
       end if
    end do
    do ia = 1, ndim-1
       do ia2 = ia+1, ndim
          if(dforce(ia) < dforce(ia2)) then
             tmp = dforce(ia)
             dforce(ia) = dforce(ia2)
             dforce(ia2) = tmp
             it = iatom_larger(ia)
             iatom_larger(ia) = iatom_larger(ia2)
             iatom_larger(ia2)  = it
          end if
       end do
    end do
    deallocate(dforce)
  end subroutine extract_largerforce_atoms

  subroutine m_Force_cal_absforc(nforc, forc_lower)
    integer, intent(in) :: nforc
    real(kind=DP), intent(out) :: forc_lower
    real(kind=DP), allocatable, dimension(:) :: dforce
!!$    integer*2 compar
    real(kind=DP) :: tmp
    integer :: ia, ia2
    allocate(dforce(natm)); dforce = 0.d0
    do ia = 1, natm
       dforce(ia) = forc_l(ia,1)*forc_l(ia,1)+forc_l(ia,2)*forc_l(ia,2)+forc_l(ia,3)*forc_l(ia,3)
    end do
    do ia = 1, natm-1
       do ia2 = ia+1, natm
          if(dforce(ia) < dforce(ia2)) then
             tmp = dforce(ia)
             dforce(ia) = dforce(ia2)
             dforce(ia2) = tmp
          end if
       end do
    end do
!!$    call qsort( dforce, natm, DP, compar)
    if(nforc < natm) then
       forc_lower = sqrt(dforce(nforc))
    else
       forc_lower = sqrt(dforce(natm))
    end if
    deallocate(dforce)
  end subroutine m_Force_cal_absforc

!!$  integer*2 function compar(a,b)
!!$    real(kind=DP) :: a,b
!!$    if(a  > b) compar = -1
!!$    if(a == b) compar = 0
!!$    if(a  < b) compar = 1
!!$    return
!!$  end function compar

  subroutine m_Force_wd_force_and_cps(nf, withorwithout, cps, ndim)
    integer, intent(in) :: nf, withorwithout,ndim
    real(kind=DP), intent(in) :: cps(ndim,3)
    integer :: ia

    if(mype==0) then
       if(withorwithout == WITHOUTTAG) then
          do ia = 1, natm
             write(nf,'(" ",i4,3f15.9,3f12.6)') ia &
                  &, cps(ia,1),cps(ia,2),cps(ia,3) &
                  &, forc_l(ia,1), forc_l(ia,2), forc_l(ia,3)
          end do
       else if(withorwithout == WITHTAG) then
          do ia = 1, natm
             write(nf,'(" !forc ",i4,3f12.6,3f12.6)') ia &
                  &, cps(ia,1),cps(ia,2),cps(ia,3) &
                  &, forc_l(ia,1), forc_l(ia,2), forc_l(ia,3)
          end do
       else
          write(nf,'(" keyword in the argument list of m_Force_wd_force_and_cps is invalid")')
       end if
    end if
  end subroutine m_Force_wd_force_and_cps

  subroutine m_Force_wd_force_and_cps_lim(nf, withorwithout, cps, ndim, N_FORC)
    integer, intent(in) :: nf, withorwithout,ndim, N_FORC
    real(kind=DP), intent(in) :: cps(ndim,3)
!!$    real(kind=DP), intent(in) :: forc_lower
!!$    real(kind=DP) :: sqr_forc, dforce, forc_lower
    real(kind=DP) :: dforce, forc_lower
    integer,       allocatable, dimension(:) :: iatm_larger ! d(min(ndim,N_FORC))
!!$    real(kind=DP), allocatable, dimension(:) :: dforce
    integer :: ia,ic

    if(withorwithout /= WITHOUTTAG .and. withorwithout /= WITHTAG) then
       if(mype == 0) write(nf,'(" keyword in the argument list of m_Force_wd_force_and_cps is invalid")')
       return
    end if
    if(mype==0) then
       allocate(iatm_larger(min(natm,n_forc))); iatm_larger = 0
       call extract_largerforce_atoms(n_forc,min(natm,n_forc),iatm_larger)
!!$       sqr_forc = forc_lower*forc_lower
       write(nf,'(" !forc *** ",i4," largest-force atoms ***")') min(n_forc,natm)
!!$       write(nf,'(" !forc   no    iatm     cps(1)      cps(2)      cps(3)    forc_l(1)   forc_l(2)   forc_l(3)  abs(forc_l)")')
       write(nf,'(" !forc   no    iatm",3(6x,"cps(",i1,")"),3(3x,"forc_l(",i1,")"),"  abs(forc_l)")') &
            & 1,2,3,1,2,3
       do ic = 1, min(natm,n_forc)
          ia = iatm_larger(ic)
          dforce = sqrt(forc_l(ia,1)*forc_l(ia,1)+forc_l(ia,2)*forc_l(ia,2)+forc_l(ia,3)*forc_l(ia,3))
          if(withorwithout == WITHOUTTAG) then
             write(nf,'(" ",i4,i8,3f15.9,4f12.6)') ic,ia,cps(ia,1:3),forc_l(ia,1:3),dforce
          else if(withorwithout == WITHTAG) then
             write(nf,'(" !forc ",i4,i8,3f12.6,4f12.6)') ic,ia,cps(ia,1:3),forc_l(ia,1:3),dforce
          end if
       end do
!!$       do ia = 1, natm
!!$          if(dforce > sqr_forc) then
!!$             ic = ic + 1
!!$             if(withorwithout == WITHOUTTAG) then
!!$                write(nf,'(" ",i4,i8,3f15.9,4f12.6)') ic,ia,cps(ia,1:3),forc_l(ia,1:3),sqrt(dforce)
!!$             else if(withorwithout == WITHTAG) then
!!$                write(nf,'(" !forc ",i4,i8,3f12.6,4f12.6)') ic,ia,cps(ia,1:3),forc_l(ia,1:3),sqrt(dforce)
!!$             end if
!!$          end if
!!$       end do
    end if
  end subroutine m_Force_wd_force_and_cps_lim

  subroutine m_Force_wd_force_cps_cpd(nf, withorwithout, cps, cpd,ndim)
    integer, intent(in) :: nf, withorwithout,ndim
    real(kind=DP), intent(in) :: cps(ndim,3), cpd(ndim,3)
    integer :: ia

    if(mype==0) then
       if(withorwithout == WITHOUTTAG) then
          do ia = 1, natm
             write(nf,'(" ",i4,3f15.9,3f12.6, 3f12.6)') ia &
                  &, cps(ia,1:3), cpd(ia,1:3), forc_l(ia,1:3)
          end do
       else if(withorwithout == WITHTAG) then
          do ia = 1, natm
             write(nf,'(" !forc ",i4,3f12.6,3f12.6,3f12.6)') ia &
                  &, cps(ia,1:3),cpd(ia,1:3),forc_l(ia,1:3)
          end do
       else
          write(nf,'(" keyword in the argument list of m_Force_wd_force_cps_cpd is invalid")')
       end if
    end if
  end subroutine m_Force_wd_force_cps_cpd

!!$  subroutine m_Force_wd_force_cps_cpd_lim(nf, withorwithout, cps, cpd,ndim, forc_lower)
  subroutine m_Force_wd_force_cps_cpd_lim(nf, withorwithout, cps, cpd,ndim, N_FORC)
    integer, intent(in) :: nf, withorwithout,ndim, N_FORC
    real(kind=DP), intent(in) :: cps(ndim,3), cpd(ndim,3)
!!$    real(kind=DP), intent(in) :: forc_lower
!!$    real(kind=DP) :: sqr_forc, dforce, forc_lower
    real(kind=DP) :: dforce, forc_lower
    integer,       allocatable, dimension(:) :: iatm_larger ! d(min(ndim,N_FORC))
!!$    real(kind=DP), allocatable, dimension(:) :: dforce
    integer :: ia, ic

    if(withorwithout /= WITHOUTTAG .and. withorwithout /= WITHTAG) then
       if(mype == 0) write(nf,'(" keyword in the argument list of m_Force_wd_force_and_cps is invalid")')
       return
    end if
    if(mype==0) then
       allocate(iatm_larger(min(natm,n_forc))); iatm_larger = 0
       call extract_largerforce_atoms(n_forc,min(natm,n_forc),iatm_larger)
!!$       sqr_forc = forc_lower*forc_lower
       write(nf,'(" !forc ",i4," largest force atoms")') min(n_forc,natm)
!!$       write(nf,'(" !forc no iatm  cps(1)  cps(2)  cps(3)   forc_l(1)  forc_l(2)  forc_l(3) abs(forc_l)")')
       write(nf,'(" !forc   no    iatm",3(6x,"cps(",i1,")"),3(6x,"cpd(",i1,")") &
            & , 3(3x,"forc_l(",i1,")"),"  abs(forc_l)")') 1,2,3,1,2,3,1,2,3
!!$       ic = 0
       do ic = 1, min(natm,n_forc)
          ia = iatm_larger(ic)
          dforce = sqrt(forc_l(ia,1)*forc_l(ia,1)+forc_l(ia,2)*forc_l(ia,2)+forc_l(ia,3)*forc_l(ia,3))
          if(withorwithout == WITHOUTTAG) then
             write(nf,'(" ",i4,i8,3f15.9,3f12.6, 4f12.6)') ic,ia,cps(ia,1:3) &
                  & ,cpd(ia,1:3), forc_l(ia,1:3),dforce
          else if(withorwithout == WITHTAG) then
             write(nf,'(" !forc ",i4,i8,3f12.6,3f12.6,4f12.6)') ic,ia,cps(ia,1:3) &
                  & ,cpd(ia,1:3), forc_l(ia,1:3),dforce
          end if
       end do
!!$       do ia = 1, natm
!!$          dforce = forc_l(ia,1)*forc_l(ia,1)+forc_l(ia,2)*forc_l(ia,2)+forc_l(ia,3)*forc_l(ia,3)
!!$          if(dforce > sqr_forc) then
!!$             ic = ic + 1
!!$             if(withorwithout == WITHOUTTAG) then
!!$                write(nf,'(" ",2i4,3f15.9,3f12.6, 4f12.6)') ic,ia,cps(ia,1:3),cpd(ia,1:3),forc_l(ia,1:3), sqrt(dforce)
!!$             else if(withorwithout == WITHTAG) then
!!$                write(nf,'(" !forc ",2i4,3f12.6,3f12.6,4f12.6)') ic,ia &
!!$                     &, cps(ia,1:3),cpd(ia,1:3),forc_l(ia,1:3), sqrt(dforce)
!!$             end if
!!$          end if
!!$       end do
    end if
  end subroutine m_Force_wd_force_cps_cpd_lim

  subroutine m_Force_wd_force_sum(nf)
    integer, intent(in) :: nf
    real(kind=DP) :: force_sum(3)
    integer :: ia,j
    do j=1,3
       force_sum(j) = 0.d0
       do ia=1,natm
          force_sum(j) = force_sum(j) + forc_l(ia,j)
       end do
    end do
    write(nf,'(" !forc Sum of forces =",3(1x,f12.6))') force_sum(1:3)
  end subroutine m_Force_wd_force_sum

  subroutine alloc_force_previous()
!!$    if(.not.allocated(force_error)) then
!!$       allocate(force_error(natm,3)); force_error = 0.d0
!!$    end if
    if(.not.allocated(force_previous)) then
       allocate(force_previous(natm,3)); force_previous = 0.d0
    end if
  end subroutine alloc_force_previous

  subroutine dealloc_force_previous()
    if(allocated(force_previous)) deallocate(force_previous)
  end subroutine dealloc_force_previous

  subroutine m_Force_get_fec_count(fec_counter)
    integer, intent(out) :: fec_counter
    fec_counter = force_error_check_counter
  end subroutine m_Force_get_fec_count

  subroutine m_Force_clear_fec_counter()
    force_error_check_counter = 0
  end subroutine m_Force_clear_fec_counter

  logical function m_Force_error_tolerable(nfout)
    integer, intent(in) :: nfout
    real(kind=DP) ::       max_norm_error, max_angle_error, force_hyper_angle_error
    if(ipriforce >= 2)  write(nfout,'("!fec  force_error_check_counter = ",i5," forcmx = ",f16.8)') &
         & force_error_check_counter, forcmx

    if(force_error_check_counter == 0) then
       call alloc_force_previous()
       m_Force_error_tolerable = .false.
    else
       call diff_force(max_norm_error,max_angle_error,force_hyper_angle_error,forccr)
       if(          max_norm_error          <= f_tolerable_norm_error &
            & .and. max_angle_error         <= f_tolerable_angle_error &
            & .and. force_hyper_angle_error <= f_tolerable_hyper_angle_error) then
          m_Force_error_tolerable = .true.
       else
          m_Force_error_tolerable = .false.
       end if
    end if

    if(m_Force_error_tolerable) then
       call dealloc_force_previous()
       force_error_check_counter = 0
    else
       call cp_to_force_previous()
       force_error_check_counter = force_error_check_counter+1
    end if
    if(ipriforce == 1)  write(nfout,'("!fec  force_error_check_counter = ",i5," forcmx = ",f16.8)') &
         & force_error_check_counter, forcmx
  contains
    subroutine cp_to_force_previous()
      force_previous = forc_l
    end subroutine cp_to_force_previous

    subroutine diff_force(max_norm_error,max_angle_error,hyper_angle_error,f_critical)
      real(kind=DP), intent(out) :: max_norm_error,max_angle_error,hyper_angle_error
      real(kind=DP), intent(in)  :: f_critical
      real(kind=DP) :: fca, f_diff, f_error,f_angle_error,f_angle_error_deg
      real(kind=DP) :: f1,f2,f3, fnrom, fnorm_p, fprod0, fprod, fnorm_p0, fnorm, fnorm0
      real(kind=DP) :: farg
      integer :: ia, j
      max_norm_error = 0.d0
      max_angle_error = 0.d0
      hyper_angle_error = 0.d0
      if(ipriforce >= 2) write(nfout,'("!fec  natom, ",9x,"  fca, ",3x," f_diff,",3x &
           & ," f_error",3x,"angle_error((rad), (deg))")')
      fprod = 0.d0; fnrom = 0.d0; fnorm_p = 0.d0
      do ia = 1, natm
         if(imdtyp(ia) == RELAX .or. imdtyp(ia) == BONDLENGTH_FIX .or. imdtyp(ia) > HEAT_BATH) then
            fnorm0 = forc_l(ia,1)**2+forc_l(ia,2)**2+forc_l(ia,3)**2
            fca = dsqrt(fnorm0)
            f_error = 0.d0; f_diff = 0.d0; f_angle_error = 0.d0; f_angle_error_deg = 0.d0
            if(f_critical > DELTA .and. fca > f_critical ) then
               f1 = force_previous(ia,1) - forc_l(ia,1)
               f2 = force_previous(ia,2) - forc_l(ia,2)
               f3 = force_previous(ia,3) - forc_l(ia,3)
               f_diff  = dsqrt(f1*f1 + f2*f2 + f3*f3)
               f_error = f_diff/fca
               if(max_norm_error < f_error) max_norm_error = f_error

               fprod0 = force_previous(ia,1)*forc_l(ia,1) &
                    &  + force_previous(ia,2)*forc_l(ia,2) + force_previous(ia,3)*forc_l(ia,3)
               fnorm_p0 = force_previous(ia,1)**2 + force_previous(ia,2)**2 + force_previous(ia,3)**2
               fnorm_p  = fnorm_p + fnorm_p0
               fnorm    = fnorm   + fnorm0
               if(fnorm_p0 > DELTA) then
                  fnorm_p0 = dsqrt(fnorm_p0)
                  farg = min(fprod0/(fca*fnorm_p0),1.d0)
!!$                  f_angle_error = dacos(fprod0/(fca*fnorm_p0))
                  f_angle_error = dacos(farg)
                  f_angle_error_deg = f_angle_error*180.d0/PAI
                  if(max_angle_error < abs(f_angle_error)) max_angle_error = abs(f_angle_error)
               end if
               fprod = fprod + fprod0
            end if
            if(ipriforce >= 2) then
               if(f_critical > DELTA .and. fca > f_critical ) then
                  if(f_error > f_tolerable_norm_error .or. abs(f_angle_error) > f_tolerable_angle_error) then
                     write(nfout,'("!fec  ",i5,f16.8,4f12.4)') &
                          & ia, fca, f_diff, f_error, f_angle_error, f_angle_error_deg
                  else
                     write(nfout,'("!fec  ",i5,f16.8,4f12.4,"   clear")') &
                          & ia, fca, f_diff, f_error, f_angle_error, f_angle_error_deg
                  end if
               else
                  write(nfout,'("!fec  ",i5,f16.8," ****    **** : fca < f_critical")') ia,fca
               end if
            end if
         end if
      end do
      if(fnorm > DELTA .and. fnorm_p > DELTA) then
         hyper_angle_error = abs(dacos(fprod/(dsqrt(fnorm)*dsqrt(fnorm_p))))
      else
         hyper_angle_error = 0.d0
      end if
      if(ipriforce >= 2) write(nfout,'("!fec  hyper_angle_error = ",f12.4," (radian) = " &
           &                  ,f12.4," (degree)")') hyper_angle_error, hyper_angle_error*180.00/PAI
    end subroutine diff_force

  end function m_Force_error_tolerable

  subroutine m_Force_dealloc
!    if(allocated(forc_l)) deallocate(forc_l)
    if(allocated(fnlxyz_l)) deallocate(fnlxyz_l)
    if(allocated(flhxcq_l)) deallocate(flhxcq_l)
    if(allocated(fdip_l)) deallocate(fdip_l)
    if(allocated(fext_l)) deallocate(fext_l)
    if(allocated(ffef_l)) deallocate(ffef_l)
    if(allocated(fexx_l)) deallocate(fexx_l)
  end subroutine m_Force_dealloc

end module m_Force
