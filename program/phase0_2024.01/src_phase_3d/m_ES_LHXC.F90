!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 581 $)
!
!  MODULE: m_ES_LHXC
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
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

module m_ES_LHXC
! $Id: m_ES_LHXC.F90 581 2018-08-01 08:38:42Z jkoga $
  use m_Electronic_Structure, only : vlhxc_l, vloc_esm
  use m_PlaneWaveBasisSet,    only : kg,kgp,gr_l
  use m_PseudoPotential,      only : psc_l, ival
  use m_Ionic_System,         only : ntyp, zfm3_l
  use m_Timing,               only : tstatc0_begin, tstatc0_end
#ifdef ENABLE_ESM_PACK
  use m_Control_Parameters,   only : nspin, ipri,iprivlhxcq, kimg, sw_dipole_correction &
                                   , sw_screening_correction, sw_external_potential &
  &                                , sw_esm, convergence_criteria
#else
  use m_Control_Parameters,   only : nspin, ipri,iprivlhxcq, kimg, sw_dipole_correction &
                                   , sw_screening_correction, sw_external_potential, convergence_criteria
#endif
  use m_Const_Parameters,     only : DP, PAI4, UP, DOWN, ON, OFF, CMPLDP, DELTA_V
  use m_Parallelization,      only : ista_kngp, iend_kngp, mype, myrank_g, mpi_ke_world, mpi_kg_world &
  &                                , mpi_chg_world, myrank_chg
  use m_Crystal_Structure,    only : univol
  use m_Dipole,               only : vdip_l, vext_l
  use m_Screening,            only : screening
  use m_External_Potential,   only : espot_g

! ========================== added by K. Tagami ===================== 11.0
  use m_Control_Parameters,   only : ndim_chgpot, ndim_magmom
! =================================================================== 11.0

  use m_Charge_Density,       only : chgq_l

  use m_FFT,                  only : fft_box_size_CD
  use m_PlaneWaveBasisSet,    only : igfp_l

! === POSITRON SCF === 2015/11/28
  use m_Control_Parameters,   only : sw_positron, positron_method
  use m_Const_Parameters,     only : positron_GGGC
  use m_Positron_Wave_Functions,        only :  pchg_l
! ==================== 2015/11/28
  use mpi

  implicit none
  real(kind=DP) :: vmax
!  include 'mpif.h'
!  61. m_ESlhxc_potential
contains

!=======================================================================
!=======================================================================
!=======================================================================
  subroutine m_ESlhxc_potential_3D(nfout,chg,vxc)
   use m_Dipole,               only : m_Dipole_potential_3D, m_Dipole_calc_3D 

    integer, intent(in)       :: nfout
    real(kind=DP), intent(in) :: chg(ista_kngp:iend_kngp,kimg,nspin)
    real(kind=DP), intent(in) :: vxc(ista_kngp:iend_kngp,kimg,nspin)
    real(kind=DP), save, allocatable :: vlhxc_old(:,:,:)

    complex(kind=CMPLDP),allocatable :: vhar(:)
    complex(kind=CMPLDP),allocatable, dimension(:,:) :: chgc
    real(kind=DP) :: ehar
    real(kind=DP) :: vv
    integer :: is,ik,i,it
    integer :: ist
    integer :: nfftcd,ierr,ig
    integer :: id_sname = -1

                                                  __TIMER_SUB_START(1050)
    call tstatc0_begin('m_ESlhxc_potential_3D ',id_sname)
    if (.not.allocated(vlhxc_old)) then
       allocate(vlhxc_old(ista_kngp:iend_kngp,kimg,nspin));vlhxc_old = 0.d0
    endif
    vlhxc_l = 0.d0
    ist = ista_kngp
    if(ist == 1) ist = 2
#ifdef ENABLE_ESM_PACK
    if(sw_esm==ON) then
       nfftcd = fft_box_size_CD(1,0)*fft_box_size_CD(2,0)*fft_box_size_CD(3,0)
       allocate(vhar(nfftcd));vhar=(0.d0,0.d0)
       allocate(chgc(kgp,nspin));chgc=(0.d0,0.d0)
       if(kimg==1)then
          do ig=ista_kngp,iend_kngp
!             chgc(ig-ista_kngp+1,1:nspin) = dcmplx(chg(ig,1,1:nspin),0.d0)
             chgc(ig,1:nspin) = dcmplx(chg(ig,1,1:nspin),0.d0)
          enddo
       else
          do ig=ista_kngp,iend_kngp
!             chgc(ig-ista_kngp+1,1:nspin) = dcmplx(chg(ig,1,1:nspin),chg(ig,2,1:nspin))
             chgc(ig,1:nspin) = dcmplx(chg(ig,1,1:nspin),chg(ig,2,1:nspin))
          enddo
       endif
       call mpi_allreduce(mpi_in_place,chgc,kgp*nspin,mpi_double_complex,mpi_sum,mpi_chg_world,ierr)
       call esm_hartree(chgc,ehar,vhar)
       vhar(:) = 0.5d0*vhar(:) !Ry -> Ha
       deallocate(chgc)
    endif
#endif
    do is = 1, nspin
       do ik = 1, kimg
                                                  __TIMER_DO_START(1064)
          if(myrank_chg==0) vlhxc_l(1,ik,is)   = vxc(1,ik,is)
#ifdef ENABLE_ESM_PACK
          if(nspin == 1 .and. sw_esm==OFF) then
#else
          if(nspin == 1) then
#endif
             do i = ist, iend_kngp  !for mpi
                vlhxc_l(i,ik,is) = vxc(i,ik,is) &
                     & +PAI4*chg(i,ik,is)/gr_l(i)**2
             end do
#ifdef ENABLE_ESM_PACK
          else if(nspin == 2.and.sw_esm==OFF) then
#else
          else if(nspin == 2) then
#endif
             do i = ist, iend_kngp  !for mpi
                vlhxc_l(i,ik,is) = vxc(i,ik,is) &
                     & +PAI4*(chg(i,ik,UP)+chg(i,ik,DOWN))/gr_l(i)**2
             end do
#ifdef ENABLE_ESM_PACK
          else if (sw_esm == ON)then
             do i = ist, iend_kngp  !for mpi
                vlhxc_l(i,ik,is) = vxc(i,ik,is) 
             end do
#endif
          endif
                                                  __TIMER_DO_STOP(1064)
                                                  __TIMER_DO_START(1065)
          do it    = 1,ntyp
             do i = ista_kngp, iend_kngp  !for mpi
                vlhxc_l(i,ik,is) &
                     & = vlhxc_l(i,ik,is)+psc_l(i,it)*zfm3_l(i,it,ik)
             end do
          end do
#ifdef ENABLE_ESM_PACK
          if(sw_esm==ON)then
             if(ik==1)then
                do i=ista_kngp,iend_kngp
                   vlhxc_l(i,1,is) = vlhxc_l(i,1,is)+dble(vhar(igfp_l(i)))+dble(vloc_esm(igfp_l(i)))
                enddo
             else if(ik==2) then
                do i=ista_kngp,iend_kngp
                   vlhxc_l(i,2,is) = vlhxc_l(i,2,is)+aimag(vhar(igfp_l(i)))+aimag(vloc_esm(igfp_l(i)))
                enddo
             endif
          endif
#endif
                                                  __TIMER_DO_STOP(1065)
       end do
    end do
#ifdef ENABLE_ESM_PACK
    if(sw_esm==ON) then
       deallocate(vhar)
    endif
#endif
    if(sw_dipole_correction==ON) then
       call m_Dipole_calc_3D(nfout)
       call m_Dipole_potential_3D(nfout,chg)
                                                  __TIMER_DO_START(1066)
       do is = 1, nspin
          do ik = 1, kimg
             do i = ista_kngp, iend_kngp  !for mpi
                vlhxc_l(i,ik,is) = vlhxc_l(i,ik,is) + vdip_l(i,ik) + vext_l(i,ik)
             end do
          end do
       end do
                                                  __TIMER_DO_STOP(1066)
    end if
    if(convergence_criteria == DELTA_V)then
      vmax = -1e+30
      do is = 1,nspin
         do ik = 1,kimg
            do i = ista_kngp, iend_kngp
               vv = dabs(vlhxc_old(i,ik,is)-vlhxc_l(i,ik,is))
               if(vv>vmax) vmax = vv
            enddo
         enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,vmax,1,mpi_double_precision,mpi_max,mpi_chg_world,ierr)
      if(iprivlhxcq >= 2) write(nfout,'(a,f20.10)') '!** vmax : ',vmax
      vlhxc_old = vlhxc_l
    endif
    if(iprivlhxcq >= 2) call m_ESlhxc_wd_vlhxc_3D(nfout)
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(1050)
  end subroutine m_ESlhxc_potential_3D

  function m_ESlhxc_delta_vmax() result (res)
    real(kind=DP) :: res 
    res = vmax
  end function m_ESlhxc_delta_vmax

  subroutine m_ESlhxc_wd_vlhxc_3D(nfout)

    integer, intent(in) :: nfout
    integer :: ispin, i, is, ie, ri
                                                  __TIMER_SUB_START(1056)
    is = ista_kngp
    ie = min(is+20, iend_kngp)
    do ispin = 1, nspin
       if(nspin == 2) write(nfout,'(" !lhxc ispin = ",i5)') ispin
       write(nfout,'(" vlhxc_l")')
       do ri = 1, kimg
          if(kimg == 1 .and. ri == 1) write(nfout,*) '       real part'
          if(ri == 2)                 write(nfout,*) '       imaginary part'
          write(nfout,'(" ",4d20.12)') (vlhxc_l(i,ri,ispin),i=is,ie) ! MPI
       end do
    end do
                                                  __TIMER_SUB_STOP(1056)
  end subroutine m_ESlhxc_wd_vlhxc_3D
!=======================================================================
end module m_ES_LHXC
