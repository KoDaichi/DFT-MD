!#define DBG_SPIN_FIX
!=======================================================================
!
!  PROGRAM  PHASE/0 2023.01
!
!  MODULE: m_ES_occup
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!  FURTHER MODIFICATION: T. Yamasaki,  April/27/2007
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
! *************************************************************
!
! =========== Contributions ===================================
!
! Through the courtesy of contributors, the following functions are added.
!
! Company:  ASMS Co.,Ltd.
! Functions:  [Identifier: 13.0E]
!                    Fermi-Dirac smearing method is added
!
! =============================================================
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

#define USE_MIDDLE_EIGENVAL

module m_ES_occup
!     (m_ESoc)
! $Id: m_ES_occup.F90 633 2020-12-01 05:11:03Z jkoga $
  use m_Electronic_Structure, only : efermi, efermi_spin, occup_l, eko_l, totch, neordr, eko_ek &
       &                            , vbm, metalic_system, check_if_metalic_flag &
       &                            , band_entropy
  use m_Kpoints,              only : kv3, kv3_ek, qwgt, qwgt_ek, vkxyz &
       &                            ,np0,np2,ip20,iwt,ip2cub,nxyz_tetra
  use m_Timing,               only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters,   only : ipri, iprioccup, width, nspin, neg, af, ekmode, idimtetra &
       &                            ,partial_charge_Emin, partial_charge_Emax, partial_charge_deltaE &
       &                            ,fixed_charge_k_parallel, icond &
       &                            ,m_CtrlP_way_of_smearing, printable, width_tetra, order_mp, esearch_factor_mp,esearch
!!$       &                            ,deltaE_dos_Gaussdistrib, variance_dos_Gaussdistrib
  use m_Crystal_Structure,    only : sw_fix_total_spin, total_spin, imag
  use m_Const_Parameters,     only : DP, DELTA_FermiSearchRange, ON, OFF, NO, YES, HARTREE, DELTA10 &
       &                            ,PARABOLIC, TETRAHEDRON, COLD, FERRO, FIXED_CHARGE, ONE_BY_ONE, BUCS
  use m_Parallelization,      only : MPI_CommGroup,npes,map_ek,mype,map_e,map_k,myrank_e,myrank_k &
       &                            ,ierr,np_e,map_z,ista_e,ista_k,iend_k, nrank_k, nrank_e &
       &                            ,mp_e, nel_e, nie_e, nis_e &
       &                            ,ista_spin, iend_spin, nrank_s, MPI_CommGroup
! ==== DEBUG by tkato 2012/04/03 ===============================================
  use m_Parallelization,      only : mpi_ge_world, nrank_g
! ==============================================================================

! ======================================= added by K. Tagami ============= 11.0
   use m_Control_Parameters,   only : noncol, ndim_spinor, ndim_magmom
! ========================================================================= 11.0

  use m_Const_Parameters,     only : Fermi_Dirac, MP, LOWEST_AT_EACH_KPT
  use mpi

  implicit none
!   1. m_ESoc_fermi_parabolic       <-(ChargeDensity_Construction) ->(3)
!      - get_occup_l_and_tot, - get_entropic_term, - get_entropy, - find_eko_minimum_and_maximum
!   2. m_ESoc_fermi_tetrahedron     <-(ChargeDensity_Construction) ->(3)
!   3. check_totch           <-(1),(2)
!   4. m_ESoc_check_num_bands
!   5. m_ESoc_fermi_parabolic_ek
!      - get_tot
!   6. m_ESoc_fermi_tetra_ek
!   7. m_ESoc_check_if_metalic
!   8. check_if_metalic
!   9. m_ESoc_fermi_ColdSmearing
!      - get_occup_l_and_tot, - coldsmearing, - get_entropic_term, - get_entropy
!  10. m_ESoc_set_nEwindows_pc
!  11. find_emin_emax
!  12. m_ESoc_keep_occup
!  13. m_ESoc_retrieve_occup
!  14. m_ESoc_substitute_occup
!  15. m_ESoc_free_nEwindows
!  16. m_ESoc_if_elec_state
!
!    The subroutine "m_ESoc_fermi_parabolic" was originally coded by
!   Noriaki Hamada, 1984.05.17
!
  integer :: nEwindows_pc
  real(kind=DP), private, allocatable, dimension(:,:) :: Ewindows_partial_charge  ! d(nEwindows_pc,2)
  real(kind=DP), private, allocatable, dimension(:,:) :: occup_tmp  ! d(np_e, ista_k:iend_k)
  integer, private, allocatable, dimension(:) ::         if_elec_state ! d(nEwindows_pc)
  real(kind=DP), private ::                              e_origin, total_spin0
  logical, private, save :: efermi_is_defined = .false.

! parameters for manual occupation
  character(len("sw_manual_occupation")),private,parameter :: tag_sw_manual_occupation="sw_manual_occupation"
  integer :: sw_manual_occupation = OFF

  real(kind=DP),allocatable,dimension(:,:) :: occ_ext
  integer,allocatable,dimension(:) :: band_index
  character(len("occupation")),private,parameter :: tag_occupation = "occupation"
  character(len("band_index")),private,parameter :: tag_band_index = "band_index"
  character(len("occ")),private,parameter :: tag_occ = "occ"
  character(len("occ_up")),private,parameter :: tag_occ_up = "occ_up"
  character(len("occ_down")),private,parameter :: tag_occ_down = "occ_down"

#ifdef DBG_SPIN_FIX
  integer, private,parameter :: PRINTOUTLEVEL_TOTCH_SPIN = 1
#else
  integer, private,parameter :: PRINTOUTLEVEL_TOTCH_SPIN = 2
#endif
! -------------------------------- added by T. Hamada --------------------------
  logical                    :: system_is_already_determined = .false.
  logical                    :: system_is_not_metallic = .false.
! ------------------------------------------------------------------------------

!  include 'mpif.h'
contains
! --> T. Yamasaki, 06 Aug. 2009
  real(kind=DP) function m_ESoc_efermi_diff()
! ============================================= modified by K. Tagami ========= 11.0
!!    if(nspin == 1) then
    if ( nspin == 1 .or. noncol ) then
! ============================================================================= 11.0
       m_ESoc_efermi_diff = 0.d0
       return
    else
       m_ESoc_efermi_diff = efermi_spin(1) - efermi_spin(2)
    end if
  end function m_ESoc_efermi_diff

  real(kind=DP) function m_ESoc_new_total_spin()
    m_ESoc_new_total_spin = (total_spin + total_spin0)*0.50
  end function m_ESoc_new_total_spin

  real(kind=DP) function m_ESoc_total_spin0()
    m_ESoc_total_spin0 =  total_spin0
  end function m_ESoc_total_spin0


  subroutine check_totch(nfout)
    integer, intent(in) :: nfout
    real(kind=DP), parameter :: minimum_critical = 1.d-14
    if(totch <= minimum_critical ) then
       write(nfout,*) '!D totch is too small: totch = ', totch
       call phase_error_with_msg(nfout,'totch is too small',__LINE__,__FILE__)
    endif
  end subroutine check_totch

!!$  subroutine m_ESoc_dos_with_GaussDistrib()
!!$
!!$  end subroutine m_ESoc_dos_with_GaussDistrib

  logical function m_ESoc_check_num_bands()
    real(kind=DP) :: c

! =========================== modified by K. Tagami ================ 11.0
!    c = neg*2.d0
!
    if ( noncol ) then
       c = neg
    else
       c = neg*2.d0
    endif
! ================================================================== 11.0

    if( c > totch) then
       m_ESoc_check_num_bands = .true.
    else
       m_ESoc_check_num_bands = .false.
    end if
  end function m_ESoc_check_num_bands

  subroutine m_ESoc_fermi_parabolic_ek(nfout)
!
!   Coded from <m_ESoc_fermi_parabolic>
!                      by T. Yamasaki (FUJITSU LABORATORIES Ltd.), 28th Jun. 2003
    integer, intent(in) :: nfout

    integer             :: jcount
    real(kind=DP)       :: emin, emax, e1, e2, tot
    real(kind=DP), parameter :: DELTA_TOTCH = 1.d-10
    integer,       parameter :: MAXITR = 100000

! ================================ added by K. Tagami ================ 11.0
    integer :: ik, is, iksnl
    real(kind=DP), allocatable :: eko_tmp(:,:)
! ==================================================================== 11.0

    integer :: id_sname = -1

    call tstatc0_begin('m_ESoc_fermi_parabolic_ek ', id_sname)

    call check_totch(nfout)

! ============================ modified by K. Tagami =================  11.0
!    emin = minval(eko_ek)
!    emax = maxval(eko_ek)
!
     if ( noncol ) then
       allocate( eko_tmp( neg, kv3_ek/ndim_spinor ) ); eko_tmp = 0.0d0
       Do ik=1, kv3_ek, ndim_spinor
         iksnl = ( ik -1 )/ndim_spinor + 1
         eko_tmp( :,iksnl ) = eko_ek( :,ik )
       End do
       emin = minval(eko_tmp);   emax = maxval(eko_tmp)
       deallocate( eko_tmp )
     else
       emin = minval(eko_ek)
       emax = maxval(eko_ek)
     endif
! ===================================================================== 11.0

    efermi = emax
    e1 = emin - dabs(DELTA_FermiSearchRange)
    e2 = emax + dabs(DELTA_FermiSearchRange)

    jcount = 1
    occupied_ch_equals_totch : do

! ======================================== modified by K. Tagami ========== 11.0
!!       call get_tot()             ! -(contained here) ->tot
!
       if ( noncol ) then
         call get_tot_noncl()
       else
         call get_tot()
       endif
! ========================================================================== 11.0

!          ~~~~~~~~~~~~~~~~~~~~
       if(jcount == 1 .and. tot < totch) &
            call wd_fermi_error1(nfout,emin,emax,tot,totch) ! -(b_Fermi)
       if(dabs(tot - totch) < DELTA_TOTCH) exit occupied_ch_equals_totch
       if( tot < totch) then
          e1 = efermi; efermi = efermi + (e2-efermi)/2
       else
          e2 = efermi; efermi = efermi + (e1-efermi)/2
       end if
       jcount = jcount + 1
       if(jcount > MAXITR) then
          call wd_fermi_error2 &   !-(b_Fermi)
               & (nfout,e1,e2,efermi,emin,emax,tot,totch,neg,MAXITR)
       end if
    end do occupied_ch_equals_totch
    if(iprioccup >= 2) &
         &write(nfout,'( " efermi = ", f10.4, " tot = ", f10.4)') efermi, tot


! ===================================== modified by K. Tagami ========== 11.0
!    call check_if_metalic(nfout,metalic_system)

    if ( noncol ) then
      call check_if_metalic_noncl(nfout,metalic_system)
    else
      call check_if_metalic(nfout,metalic_system)
    endif
! ======================================================================= 11.0

    if(iprioccup >= 2 .and. .not.metalic_system) then
       if(dabs(vbm) < 1.d6) then
          write(nfout,'( " The highest occupied band energy = ", f10.4)') vbm
       else
          write(nfout,'( " The highest occupied band energy is not calculated")')
       end if
    end if

    call tstatc0_end(id_sname)
  contains
    subroutine get_tot
      integer       :: k, i
      real(kind=DP) :: wspin = 1.d0, e, dos, weight, totw
      tot = 0.d0
      do k = 1, kv3_ek, af+1
         do i = 1, neg
            e = eko_ek(i,k)
            call width2(e,efermi,width,dos,weight)  ! -(b_Fermi)
            totw = weight*wspin*kv3_ek*qwgt_ek(k)
            tot = tot + 2*totw
         end do
      end do
      tot = tot/kv3_ek * (af+1)
      if(iprioccup >= 2) write(nfout,'(" efermi, tot = ",d20.12,d20.12," jcount = ",i7)') efermi,tot, jcount
    end subroutine get_tot

! ======================================== added by K. Tagami ============ 11.0
    subroutine get_tot_noncl
      integer       :: k, i
      real(kind=DP) :: wspin = 1.d0, e, dos, weight, totw, wktmp

      tot = 0.d0
      wktmp = kv3_ek / ndim_spinor

      do k = 1, kv3_ek, ndim_spinor
         do i = 1, neg
            e = eko_ek(i,k)
            call width2(e,efermi,width,dos,weight)  ! -(b_Fermi)
            totw = weight *wspin *wktmp *qwgt_ek(k)
            tot = tot + totw
         end do
      end do
      tot = tot /wktmp
      if(iprioccup >= 2) write(nfout,'(" efermi, tot = ",d20.12,d20.12," jcount = ",i7)') efermi,tot, jcount

    end subroutine get_tot_noncl
! ======================================================================= 11.0

  end subroutine m_ESoc_fermi_parabolic_ek

  subroutine m_ESoc_fermi_tetra_ek(nfout)
!  Coded from <m_ESoc_fermi_tetrahedron>
!          by Takahiro Yamasaki(FUJITSU LABORATORIES Ltd.), 28th Jun. 2003
!
    integer, intent(in) :: nfout
#ifndef NO_TETRAHEDRON
    real(kind=DP), parameter :: delta = 1.d-12
    integer, parameter       :: idim = 3
    integer        :: neig,nengy,ispin,ip2,ik,instts1,ikee,nxx,nyy,nzz, ib, jb
    real(kind=DP)  :: efermi2,eval,totind, et
    real(kind=DP), pointer, dimension(:,:,:) :: eig2
    real(kind=DP), pointer, dimension(:) :: eawk,cdwk,cswk,cdos,cind,valud
    integer             :: id_sname = -1

! =================================== added by K. Tagami ==============- 11.0
    integer :: is
! ====================================================================== 11.0

    call tstatc0_begin('m_ESoc_fermi_tetra_ek ', id_sname)

! ================================= modified by K. Tagami ============= 11.0
!    allocate(eig2(np2,neg,nspin)); eig2 = 0.d0
    if ( noncol ) then
      allocate(eig2(np2,neg,1)); eig2 = 0.d0
    else
      allocate(eig2(np2,neg,nspin)); eig2 = 0.d0
    endif
! ===================================================================== 11.0

    allocate(eawk(np0)); eawk = 0.d0
    allocate(cdwk(np0)); cdwk = 0.d0
    allocate(cswk(np0)); cswk = 0.d0
    allocate(cdos(np2*neg)); cdos = 0.d0
    allocate(cind(np2*neg)); cind = 0.d0
! ================================= modified by K. Tagami ============= 11.0
!    allocate(valud(nspin)); valud = 0.d0
    if ( noncol ) then
      allocate(valud(1)); valud = 0.d0
    else
      allocate(valud(nspin)); valud = 0.d0
    endif
! ====================================================================== 11.0

    call check_totch(nfout)           ! totch is checked

    nxx = nxyz_tetra(1)
    nyy = nxyz_tetra(2)
    nzz = nxyz_tetra(3)
    neig=neg
    nengy=0

! ======================================= modified by K. Tagami ======== 11.0
!    do ispin=1,nspin
!       do ip2=1,np2
!          ik=nspin*(ip2-1)+ispin
!          do ib=1,neg
!             eig2(ip2,ib,ispin)=eko_ek(ib,ik)
!          enddo
!          do ib = 1,neg-1
!             do jb = ib+1, neg
!                if(eig2(ip2,jb,ispin) < eig2(ip2,ib,ispin)-delta) then
!                   et = eig2(ip2,ib,ispin)
!                   eig2(ip2,ib,ispin) = eig2(ip2,jb,ispin)
!                   eig2(ip2,jb,ispin) = et
!                end if
!             end do
!          end do
!       enddo
!    enddo
    if ( noncol ) then
      do is=1, 1
         do ip2=1,np2
            ik = ndim_spinor *(ip2-1)+is
            do ib=1,neg
               eig2(ip2,ib,is)=eko_ek(ib,ik)
            enddo
            do ib = 1,neg-1
               do jb = ib+1, neg
                  if(eig2(ip2,jb,is) < eig2(ip2,ib,is)-delta) then
                     et = eig2(ip2,ib,is)
                     eig2(ip2,ib,is) = eig2(ip2,jb,is)
                     eig2(ip2,jb,is) = et
                  end if
               end do
            end do
         enddo
      enddo
    else
      do ispin=1,nspin
         do ip2=1,np2
            ik=nspin*(ip2-1)+ispin
            do ib=1,neg
               eig2(ip2,ib,ispin)=eko_ek(ib,ik)
            enddo
            do ib = 1,neg-1
               do jb = ib+1, neg
                  if(eig2(ip2,jb,ispin) < eig2(ip2,ib,ispin)-delta) then
                     et = eig2(ip2,ib,ispin)
                     eig2(ip2,ib,ispin) = eig2(ip2,jb,ispin)
                     eig2(ip2,jb,ispin) = et
                  end if
               end do
            end do
         enddo
      enddo
    endif
! ==================================================================== 11.0

    if(iprioccup>=2) then
       write(nfout,'(" --- eko_ek --- <<m_ESoc_fermi_tetra_ek>>")')
! ================================== modified by K. Tagami =========== 11.0
!       do ik = 1, kv3_ek
       do ik = 1, kv3_ek, ndim_spinor
! ===================================================================== 11.0
          write(nfout,'(" ik = ",i8)') ik
          write(nfout,'(10f8.4)') (eko_ek(ib,ik),ib=1,neg)
       end do
    end if

    if(iprioccup>=2) write(nfout,*) ' === tetrahedron method', &
                &  ' for k-space integration ==='
    efermi2 = efermi

! ============================= modified by K. Tagami ============== 11.0
!    call fermi1(nxx,nyy,nzz,np2,np2,neig,neg,nspin, &
!             &  eig2,ip20,np0,totch, efermi2,eval,valud, &
!             &  iwt,ip2cub,iprioccup)
    if ( noncol ) then
      call fermi1( nxx, nyy, nzz, np2, np2, neig, neg, 1, &
               &   eig2, ip20, np0, totch, efermi2, eval, valud, &
               &   iwt, ip2cub, iprioccup, .true. )
    else
      call fermi1( nxx, nyy, nzz, np2, np2, neig, neg, nspin, &
               &   eig2, ip20, np0, totch, efermi2, eval, valud, &
               &   iwt, ip2cub, iprioccup, .false. )
    endif
! ================================================================== 11.0

    if(iprioccup>=2) write(nfout,*) 'eval=',eval

! ============================= modified by K. Tagami ============== 11.0
!    do ispin=1,nspin
!!!$      write(nfout,*) ' ispin=',ispin
!      instts1 = 0
!      call nsdos3(nfout,idim,efermi2,nxx,nyy,nzz,np2,1,neig, &
!               &  eig2(1,1,ispin),ip20,np0,eawk,instts1,np2,iwt,ip2cub,iprioccup)
!      call nstt3i(idim,nengy,efermi2,nxx,nyy,nzz, &
!               &  np2,np2,neig,eig2(1,1,ispin), &
!               &  ip20,np0,eawk,cdwk,cswk,np2,neig,cdos,cind )
!      efermi=efermi2
!
!      totind=0.d0
!      do ip2=1,np2
!         do ib=1,neig
!            ikee=ip2+np2*(ib-1)
!            totind=totind+cind(ikee)
!         enddo
!      enddo
!      if(iprioccup>=2) write(nfout,*) 'for spin=',ispin, &
!                  &  ' ** TOTAL CHARGE after fermi1 = ',totind
!   enddo
!
   if ( noncol ) then
      do is=1,1
        instts1 = 0
        call nsdos3( nfout, idim, efermi2, nxx, nyy, nzz, np2, 1, neig, &
                 &  eig2(1,1,is), ip20, np0, eawk, instts1, np2, &
 	         &  iwt, ip2cub, iprioccup )
        call nstt3i( idim, nengy, efermi2, nxx, nyy, nzz, &
                 &   np2, np2, neig, eig2(1,1,is), &
                 &   ip20, np0, eawk, cdwk, cswk, np2, neig, cdos, cind, width_tetra )
        efermi=efermi2

        totind=0.d0
        do ip2=1,np2
           do ib=1,neig
              ikee=ip2+np2*(ib-1)
              totind=totind+cind(ikee)
           enddo
        enddo
        if(iprioccup>=2) write(nfout,*) 'for spin=',is, &
                    &  ' ** TOTAL CHARGE after fermi1 = ',totind
     enddo
   else
      do ispin=1,nspin
        instts1 = 0
        call nsdos3( nfout, idim, efermi2, nxx, nyy, nzz, np2, 1, neig, &
                 &  eig2(1,1,ispin), ip20, np0, eawk, instts1, np2, &
 	         &  iwt, ip2cub, iprioccup )
        call nstt3i( idim, nengy, efermi2, nxx, nyy, nzz, &
                 &   np2, np2, neig, eig2(1,1,ispin), &
                 &   ip20, np0, eawk, cdwk, cswk, np2, neig, cdos, cind, width_tetra )
        efermi=efermi2

        totind=0.d0
        do ip2=1,np2
           do ib=1,neig
              ikee=ip2+np2*(ib-1)
              totind=totind+cind(ikee)
           enddo
        enddo
        if(iprioccup>=2) write(nfout,*) 'for spin=',ispin, &
                    &  ' ** TOTAL CHARGE after fermi1 = ',totind
     enddo
    endif
! =================================================================== 11.0

    deallocate(valud); deallocate(cind); deallocate(cdos); deallocate(cswk)
    deallocate(cdwk);  deallocate(eawk)
    deallocate(eig2)

! ========================================= modified by K. Tagami ====== 11.0
!    call check_if_metalic(nfout,metalic_system)
    if ( noncol ) then
      call check_if_metalic_noncl(nfout,metalic_system)
    else
      call check_if_metalic(nfout,metalic_system)
    endif
! ===================================================================== 11.0

    if(iprioccup >= 2 .and. .not.metalic_system) &
         &write(nfout,'( " The highest occupied band energy = ", f10.4)') vbm

    call tstatc0_end(id_sname)
#endif
  end subroutine m_ESoc_fermi_tetra_ek

! ============================== KT_add ====================== 13.0E
  subroutine m_ESoc_fermi_dirac_ek(nfout)
    integer, intent(in) :: nfout

    integer             :: jcount
    real(kind=DP)       :: emin, emax, e1, e2, tot
    real(kind=DP), parameter :: DELTA_TOTCH = 1.d-10
    integer,       parameter :: MAXITR = 100000

    integer :: ik, is, iksnl
    real(kind=DP), allocatable :: eko_tmp(:,:)

    integer :: id_sname = -1

    call tstatc0_begin('m_ESoc_fermi_dirac_ek ', id_sname)

    call check_totch(nfout)

    if ( noncol ) then
       allocate( eko_tmp( neg, kv3_ek/ndim_spinor ) ); eko_tmp = 0.0d0
       Do ik=1, kv3_ek, ndim_spinor
          iksnl = ( ik -1 )/ndim_spinor + 1
          eko_tmp( :,iksnl ) = eko_ek( :,ik )
       End do
       emin = minval(eko_tmp);   emax = maxval(eko_tmp)
       deallocate( eko_tmp )
    else
       emin = minval(eko_ek)
       emax = maxval(eko_ek)
    endif

    efermi = emax
    e1 = emin - dabs(DELTA_FermiSearchRange)
    e2 = emax + dabs(DELTA_FermiSearchRange)

    jcount = 1

! --------------------------------------------
    occupied_ch_equals_totch : do

       if ( noncol ) then
          call get_tot_noncl()
       else
         call get_tot()
       endif

!          ~~~~~~~~~~~~~~~~~~~~
       if(jcount == 1 .and. tot < totch) &
            call wd_fermi_error1(nfout,emin,emax,tot,totch) ! -(b_Fermi)
       if(dabs(tot - totch) < DELTA_TOTCH) exit occupied_ch_equals_totch
       if( tot < totch) then
          e1 = efermi; efermi = efermi + (e2-efermi)/2
       else
          e2 = efermi; efermi = efermi + (e1-efermi)/2
       end if
       jcount = jcount + 1
       if(jcount > MAXITR) then
          call wd_fermi_error2 &   !-(b_Fermi)
               & (nfout,e1,e2,efermi,emin,emax,tot,totch,neg,MAXITR)
       end if
    end do occupied_ch_equals_totch
    if(iprioccup >= 2) &
         &write(nfout,'( " efermi = ", f10.4, " tot = ", f10.4)') efermi, tot


    if ( noncol ) then
      call check_if_metalic_noncl(nfout,metalic_system)
    else
      call check_if_metalic(nfout,metalic_system)
    endif

    if(iprioccup >= 2 .and. .not.metalic_system) then
       if(dabs(vbm) < 1.d6) then
          write(nfout,'( " The highest occupied band energy = ", f10.4)') vbm
       else
          write(nfout,'( " The highest occupied band energy is not calculated")')
       end if
    end if

    call tstatc0_end(id_sname)
  contains
    subroutine get_tot
      integer       :: k, i
      real(kind=DP) :: wspin = 1.d0, e, dos, weight, totw
      tot = 0.d0
      do k = 1, kv3_ek, af+1
         do i = 1, neg
            e = eko_ek(i,k)
            call width_fermi_dirac(e,efermi,width,dos,weight)  ! -(b_Fermi)
            totw = weight*wspin*kv3_ek*qwgt_ek(k)
            tot = tot + 2*totw
         end do
      end do
      tot = tot/kv3_ek * (af+1)
      if(iprioccup >= 2) write(nfout,'(" efermi, tot = ",d20.12,d20.12," jcount = ",i7)') efermi,tot, jcount
    end subroutine get_tot

    subroutine get_tot_noncl
      integer       :: k, i
      real(kind=DP) :: wspin = 1.d0, e, dos, weight, totw, wktmp

      tot = 0.d0
      wktmp = kv3_ek / ndim_spinor

      do k = 1, kv3_ek, ndim_spinor
         do i = 1, neg
            e = eko_ek(i,k)
            call width_fermi_dirac(e,efermi,width,dos,weight)  ! -(b_Fermi)
            totw = weight *wspin *wktmp *qwgt_ek(k)
            tot = tot + totw
         end do
      end do
      tot = tot /wktmp
      if(iprioccup >= 2) write(nfout,'(" efermi, tot = ",d20.12,d20.12," jcount = ",i7)') efermi,tot, jcount

    end subroutine get_tot_noncl

  end subroutine m_ESoc_fermi_dirac_ek
! ========================================================== 13.0E

  subroutine m_ESoc_mp_ek(nfout)
!
!   Coded from <m_ESoc_fermi_parabolic>
!                      by T. Yamasaki (FUJITSU LABORATORIES Ltd.), 28th Jun. 2003
    integer, intent(in) :: nfout

    integer             :: jcount
    real(kind=DP)       :: emin, emax, e1, e2, tot
    real(kind=DP), parameter :: DELTA_TOTCH = 1.d-10
    integer,       parameter :: MAXITR = 100000

! ================================ added by K. Tagami ================ 11.0
    integer :: ik, is, iksnl
    real(kind=DP), allocatable :: eko_tmp(:,:)
! ==================================================================== 11.0

    integer :: id_sname = -1

    call tstatc0_begin('m_ESoc_mp_ek ', id_sname)

    call check_totch(nfout)

! ============================ modified by K. Tagami =================  11.0
!    emin = minval(eko_ek)
!    emax = maxval(eko_ek)
!
     if ( noncol ) then
       allocate( eko_tmp( neg, kv3_ek/ndim_spinor ) ); eko_tmp = 0.0d0
       Do ik=1, kv3_ek, ndim_spinor
         iksnl = ( ik -1 )/ndim_spinor + 1
         eko_tmp( :,iksnl ) = eko_ek( :,ik )
       End do
       emin = minval(eko_tmp);   emax = maxval(eko_tmp)
       deallocate( eko_tmp )
     else
       emin = minval(eko_ek)
       emax = maxval(eko_ek)
     endif
! ===================================================================== 11.0

    efermi = emax
    e1 = emin - dabs(DELTA_FermiSearchRange)
    e2 = emax + dabs(DELTA_FermiSearchRange)

    jcount = 1
    occupied_ch_equals_totch : do

! ======================================== modified by K. Tagami ========== 11.0
!!       call get_tot()             ! -(contained here) ->tot
!
       if ( noncol ) then
         call get_tot_noncl()
       else
         call get_tot()
       endif
! ========================================================================== 11.0

!          ~~~~~~~~~~~~~~~~~~~~
       if(jcount == 1 .and. tot < totch) &
            call wd_fermi_error1(nfout,emin,emax,tot,totch) ! -(b_Fermi)
       if(dabs(tot - totch) < DELTA_TOTCH) exit occupied_ch_equals_totch
       if( tot < totch) then
          e1 = efermi; efermi = efermi + (e2-efermi)/2
       else
          e2 = efermi; efermi = efermi + (e1-efermi)/2
       end if
       jcount = jcount + 1
       if(jcount > MAXITR) then
          call wd_fermi_error2 &   !-(b_Fermi)
               & (nfout,e1,e2,efermi,emin,emax,tot,totch,neg,MAXITR)
       end if
    end do occupied_ch_equals_totch
    if(iprioccup >= 2) &
         &write(nfout,'( " efermi = ", f10.4, " tot = ", f10.4)') efermi, tot


! ===================================== modified by K. Tagami ========== 11.0
!    call check_if_metalic(nfout,metalic_system)

    if ( noncol ) then
      call check_if_metalic_noncl(nfout,metalic_system)
    else
      call check_if_metalic(nfout,metalic_system)
    endif
! ======================================================================= 11.0

    if(iprioccup >= 2 .and. .not.metalic_system) then
       if(dabs(vbm) < 1.d6) then
          write(nfout,'( " The highest occupied band energy = ", f10.4)') vbm
       else
          write(nfout,'( " The highest occupied band energy is not calculated")')
       end if
    end if

    call tstatc0_end(id_sname)
  contains
    subroutine get_tot
      integer       :: k, i
      real(kind=DP) :: wspin = 1.d0, e, dos, weight, totw,ent
      tot = 0.d0
      do k = 1, kv3_ek, af+1
         do i = 1, neg
            e = eko_ek(i,k)
            call width_methfessel_paxton(order_mp,e,efermi,width,dos,weight,ent)  ! -(b_Fermi)
            totw = weight*wspin*kv3_ek*qwgt_ek(k)
            tot = tot + 2*totw
         end do
      end do
      tot = tot/kv3_ek * (af+1)
      if(iprioccup >= 2) write(nfout,'(" efermi, tot = ",d20.12,d20.12," jcount = ",i7)') efermi,tot, jcount
    end subroutine get_tot

! ======================================== added by K. Tagami ============ 11.0
    subroutine get_tot_noncl
      integer       :: k, i
      real(kind=DP) :: wspin = 1.d0, e, dos, weight, totw, wktmp, ent

      tot = 0.d0
      wktmp = kv3_ek / ndim_spinor

      do k = 1, kv3_ek, ndim_spinor
         do i = 1, neg
            e = eko_ek(i,k)
            call width_methfessel_paxton(order_mp,e,efermi,width,dos,weight,ent)  ! -(b_Fermi)
            totw = weight *wspin *wktmp *qwgt_ek(k)
            tot = tot + totw
         end do
      end do
      tot = tot /wktmp
      if(iprioccup >= 2) write(nfout,'(" efermi, tot = ",d20.12,d20.12," jcount = ",i7)') efermi,tot, jcount

    end subroutine get_tot_noncl
! ======================================================================= 11.0

  end subroutine m_ESoc_mp_ek

  subroutine m_ESoc_check_if_metalic(nfout)
!
!    Modified by T. Yamasaki   April/27/2007
!
    integer, intent(in) :: nfout
    integer :: way_of_smearing

    if(ekmode == ON .or. &
         & (icond>=FIXED_CHARGE .and. ekmode==OFF .and. fixed_charge_k_parallel == ONE_BY_ONE)) then
       way_of_smearing = m_CtrlP_way_of_smearing()
       if(way_of_smearing == PARABOLIC) then
          call m_ESoc_fermi_parabolic_ek(nfout)  ! -> efermi, metalic_system
!!$  else if(way_of_smearing == MP) then
!!$     call fermi_mesfessel_paxton(nfout)
       else if(way_of_smearing == TETRAHEDRON) then
          call m_ESoc_fermi_tetra_ek(nfout)     ! -> efermi, metalic_system

! ============================ KT_add ===================== 13.0E
       else if(way_of_smearing == Fermi_Dirac) then
          call m_ESoc_fermi_dirac_ek(nfout)
! ========================================================= 13.0E
       else if(way_of_smearing == MP) then
          call m_ESoc_mp_ek(nfout)
       else if(way_of_smearing == LOWEST_AT_EACH_KPT) then
         call phase_error_with_msg( nfout, &
              &  ' m_ESoc_check_if_metalic does not work with LOWEST_AT_EACH_KPT'  &
              &   ,__LINE__,__FILE__)
       end if
    else
       way_of_smearing = m_CtrlP_way_of_smearing()
       if(way_of_smearing == PARABOLIC) then
          call m_ESoc_fermi_parabolic_3D(nfout)
!!$  else if(way_of_smearing == MP) then
!!$     call fermi_mesfessel_paxton(nfout)
       else if(way_of_smearing == TETRAHEDRON) then
          call m_ESoc_fermi_tetrahedron_3D(nfout)
       else if(way_of_smearing == COLD) then
          call m_ESoc_fermi_ColdSmearing_3D(nfout)

! =========== KT_add=============================== 13.0E
       else if( way_of_smearing == Fermi_Dirac ) then
          call m_ESoc_fermi_Dirac_3D(nfout)
! ============== ================================== 13.0E
       else if(way_of_smearing == MP) then
          call m_ESoc_methfessel_paxton(nfout)
       else if(way_of_smearing == LOWEST_AT_EACH_KPT) then
          call m_ESoc_occup_fix_3D(nfout)
       end if
    end if

! ====================================== modified by K. Tagami ========== 11.0
!    call check_if_metalic(nfout,metalic_system)
    if ( noncol ) then
      call check_if_metalic_noncl(nfout,metalic_system)
    else
      call check_if_metalic(nfout,metalic_system)
    endif
! ====================================================================== 11.0

  end subroutine m_ESoc_check_if_metalic

  subroutine check_if_metalic(nfout,metalic)
    use m_Parallelization, only : mpi_kg_world
!
!     Modified by T. Yamasaki   May/08/2007
!
    integer, intent(in) :: nfout
    logical, intent(out) :: metalic
    integer :: ieig, ie, ik, is, ihob0, kv3_t, way_of_smearing
    real(kind=DP) :: totw, emax, cbm, vbm0, weight0
    real(kind=DP), allocatable, dimension(:) :: ehob
    real(kind=DP), allocatable, dimension(:,:) :: eko_t, eko_mpi     ! d(neg,kv3|kv3_ek)
    real(kind=DP), allocatable, dimension(:,:) :: occup_t, occup_mpi ! d(neg,kv3|kv3_ek)
    integer, allocatable, dimension(:) :: ihob ! d(kv3|kv3_ek)

    if(ekmode == ON) then
       kv3_t = kv3_ek
    else
       kv3_t = kv3
    end if

    allocate(eko_t(neg,kv3_t))

    if(ekmode == ON) then
       eko_t = eko_ek
    else
       eko_t = 0.d0
       do is = ista_spin, iend_spin
       do ik = is, kv3-nspin+is, nspin
          if(map_k(ik) /= myrank_k) cycle
          do ie = 1, neg
             if(map_e(ie) /= myrank_e) cycle
             eko_t(ie,ik) = eko_l(map_z(ie),ik)
          end do
       end do
       end do
       if(npes > 1) then
          allocate(eko_mpi(neg,kv3_t))
          call mpi_allreduce(eko_t,eko_mpi,neg*kv3,mpi_double_precision &
               &               ,mpi_sum,mpi_kg_world,ierr)
          eko_t = eko_mpi
          call mpi_allreduce(eko_t,eko_mpi,neg*kv3,mpi_double_precision &
               &               ,mpi_sum,mpi_ge_world,ierr)
          eko_t = eko_mpi
          deallocate(eko_mpi)
       end if
    end if

    metalic = .false.

    allocate(ehob(kv3_t)); ehob = -999.d+20

    if(kv3_t/nspin >= 2) then
       allocate(ihob(kv3_t)); ihob = 0
       do is = 1, nspin
          do ik = is, kv3_t, nspin
             totw = 0.d0
             ihob(ik) = 0
             firstKpoint: do ieig = 1, neg
!!$                if(eko_t(ieig,is) < efermi) then
                if(eko_t(ieig,ik) <= efermi) then
                   totw = totw + 1.d0
                   ihob(ik) = ihob(ik) + 1
                   if(ehob(ik) < eko_t(ieig,ik)) ehob(ik) = eko_t(ieig,ik)
                end if
             end do firstKpoint
             totw = totw*2.d0
             if(iprioccup >= 2) write(nfout,'(" ---  ik = ",i3," ihob = ",i6," totw = ",f12.4," totch = ",f12.4)') &
                  & ik,ihob(ik),totw,totch
          end do
          do ik = nspin+is ,kv3_t, nspin
             if(ihob(ik) /= ihob(is)) then
                metalic = .true.
                exit
             end if
          end do
       end do
       deallocate(ihob)
    else
       do is = 1, nspin
          totw = 0.d0
          do ieig=1,neg
             if(eko_t(ieig,is) <= efermi) then
                if(ehob(is) < eko_t(ieig,is)) ehob(is) = eko_t(ieig,is)
             end if
          end do
       end do
    end if

    if(.not.metalic .and. ekmode == OFF) then
       allocate(occup_t(neg,kv3_t))

       occup_t = 0.d0
       do is = ista_spin, iend_spin
       do ik = is, kv3-nspin+is, nspin
          if(map_k(ik) /= myrank_k) cycle
          do ie = 1, neg
             if(map_e(ie) /= myrank_e) cycle
             occup_t(ie,ik) = occup_l(map_z(ie),ik)
          end do
       end do
       end do
       if(npes > 1) then
          allocate(occup_mpi(neg,kv3_t))
          call mpi_allreduce(occup_t,occup_mpi,neg*kv3,mpi_double_precision &
               &               ,mpi_sum,mpi_kg_world,ierr)
          occup_t = occup_mpi
          call mpi_allreduce(occup_t,occup_mpi,neg*kv3,mpi_double_precision &
               &               ,mpi_sum,mpi_ge_world,ierr)
          occup_t = occup_mpi
          deallocate(occup_mpi)
       end if

       Kpoint2: do ik = 1 ,kv3_t
          ihob0 = 0; vbm0 = -999.d+20
          findvbm: do ieig=1, neg   ! finding Valence Band Maximum
             if(occup_t(ieig,ik) > DELTA10) then
                if(vbm0 < eko_t(ieig,ik)) then
                   vbm0 = eko_t(ieig,ik)
                   ihob0 = ieig
                end if
             end if
          end do findvbm
          weight0 = kv3*qwgt(ik)
          if(iprioccup >= 2) then
             write(nfout,'(" --- ihob0 = ",i8)') ik,ihob0
             if(ihob0 >= 1) write(nfout,'(" --- ik = ",i4," ihob0 = ",i8," vbm0, occup = "&
                  & ,2f10.5," weight0=",f8.4)') ik,ihob0,vbm0,occup_t(ihob0,ik),weight0
          end if
          if(ihob0 >= 1 .and. ihob0 <= neg) then
             if(occup_t(ihob0,ik) > DELTA10 .and. occup_t(ihob0,ik) < weight0-DELTA10) then
                metalic = .true.
                exit Kpoint2
             end if
          end if
       end do Kpoint2
       deallocate(occup_t)
    end if

    way_of_smearing = m_CtrlP_way_of_smearing()

    if(.not.metalic .or. way_of_smearing == Fermi_Dirac ) then
       cbm = 1.0D10
       Kpoint: do ik=1,kv3_t
          findcbm: do ieig = 1, neg
#if 0
             if(eko_t(ieig,ik) > efermi) then
                cbm = eko_t(ieig,ik)
                exit findcbm
             end if
#else
             if(eko_t(ieig,ik) > efermi) then
                cbm = min(eko_t(ieig,ik),cbm)
                exit findcbm
             end if
#endif
          end do findcbm

          if((way_of_smearing==PARABOLIC .and. cbm < efermi + width/2) &
               & .or. (way_of_smearing==TETRAHEDRON .and. cbm < efermi+DELTA10) &
! ========== KT_add================================= 13.0E
               & .or. (way_of_smearing==Fermi_Dirac .and. cbm < efermi +width/2) &
               & .or. (way_of_smearing==LOWEST_AT_EACH_KPT .and. cbm < efermi ) &
! ================================================== 13.0E
               & .or. (way_of_smearing==COLD .and. cbm < efermi + width/2)) then

             metalic = .true.
             if(iprioccup >= 2) then
                write(nfout,'(" --- cbm = ",f8.4," ik = ",i8)') cbm,ik
             end if
             exit Kpoint
          end if
       end do Kpoint
    end if

    if(metalic) then
       if(iprioccup>=1) write(nfout,'(" --- The system is metalic ---")')
! --------------------------- added by T. Hamada -------------------------
       system_is_not_metallic = .false.
! ------------------------------------------------------------------------
       if ( way_of_smearing == Fermi_Dirac ) then
          emax = maxval(ehob)
          if(iprioccup>=1) then
             write(nfout,'(A,f12.8)') " --- The highest band energy below Ef = ",emax
             write(nfout,'(A,f12.8)') "         lowest  band energy above Ef = ",cbm
             write(nfout,'(A,f12.8,A,f12.8,A)') &
                  &                   "                           Difference = ", &
                  &           (cbm -emax), " Ha   ( ", (cbm -emax)*Hartree, " eV )"
             write(nfout,*)
          endif
       endif
    else
       if(iprioccup>=1) write(nfout,'(" --- The system is insulating or semiconducting ---")')
       emax = maxval(ehob)
#if 0
       if(iprioccup>=1) write(nfout,'(" --- The highest occupied band energy (=vbm) = ",f12.8)') emax
#else
       if(iprioccup>=1) then
          write(nfout,'(A,f12.8)') " --- The highest  occupied band energy (=vbm) = ",emax
          write(nfout,'(A,f12.8)') "         lowest unoccupied band energy (=cbm) = ",cbm
          write(nfout,'(A,f12.8,A,f12.8,A)') &
               &           "                           estimated band gap = ", &
               &           (cbm -emax), " Ha   ( ", (cbm -emax)*Hartree, " eV )"
          write(nfout,*)
       endif
#endif
       vbm = emax
! -------------------------- added by T. Hamada --------------------------
       system_is_not_metallic = .true.
! ------------------------------------------------------------------------
    end if

    check_if_metalic_flag = .true.
! ------------------------- added by T. Hamada ---------------------------
    system_is_already_determined = .true.
!-------------------------------------------------------------------------a

    deallocate(ehob)
    deallocate(eko_t)

  end subroutine check_if_metalic

! ======================================= added by K. Tagami ========= 11.0
  subroutine check_if_metalic_noncl( nfout, metalic )
!
!     Modified by T. Yamasaki   May/08/2007
!
    integer, intent(in) :: nfout
    logical, intent(out) :: metalic
    integer :: ieig, ie, ik, ihob0, kv3_t, way_of_smearing, is

    integer  :: iksnl, iktmp

    real(kind=DP) :: totw, emax, cbm, vbm0, weight0
    real(kind=DP), allocatable, dimension(:) :: ehob
    real(kind=DP), allocatable, dimension(:,:) :: eko_t, eko_mpi     ! d(neg,kv3|kv3_ek)
    real(kind=DP), allocatable, dimension(:,:) :: occup_t, occup_mpi ! d(neg,kv3|kv3_ek)
    integer, allocatable, dimension(:) :: ihob ! d(kv3|kv3_ek)

    if(ekmode == ON) then
       kv3_t = kv3_ek
    else
       kv3_t = kv3
    end if

    allocate(eko_t(neg,kv3_t/ndim_spinor));  eko_t = 0.0d0

    if(ekmode == ON) then
! -------------------------- ktDEBUG ---------------------- 20121025
!       eko_t = eko_ek
!
       Do ik=1, kv3_t, ndim_spinor
          iksnl = ( ik -1 ) /ndim_spinor + 1
          eko_t(:,iksnl) = eko_ek(:,ik)
       End do
! -------------------------- ktDEBUG ---------------------- 20121025
    else
       eko_t = 0.d0
       do ik = 1, kv3, ndim_spinor
          if(map_k(ik) /= myrank_k) cycle
 	  iksnl = ( ik -1 )/ndim_spinor + 1
          do ie = 1, neg
             if(map_e(ie) /= myrank_e) cycle
             eko_t(ie,iksnl) = eko_l(map_z(ie),ik)
          end do
       end do
       if(npes > 1) then
          allocate(eko_mpi(neg,kv3_t/ndim_spinor))
          call mpi_allreduce( eko_t, eko_mpi, neg*kv3/ndim_spinor, &
	       &              mpi_double_precision, mpi_sum,MPI_CommGroup, ierr )
          eko_t = eko_mpi
          deallocate(eko_mpi)
       end if
    end if

    metalic = .false.

    allocate(ehob(kv3_t/ndim_spinor)); ehob = -999.d+20

    if ( kv3_t/ndim_spinor >= 2 ) then
       allocate(ihob(kv3_t/ndim_spinor)); ihob = 0

       do is = 1, 1
          do ik = 1, kv3_t /ndim_spinor
             totw = 0.d0
             ihob(ik) = 0
             firstKpoint: do ieig = 1, neg
!!$                if(eko_t(ieig,is) < efermi) then
                if(eko_t(ieig,ik) < efermi) then
                   totw = totw + 1.d0
                   ihob(ik) = ihob(ik) + 1
                   if(ehob(ik) < eko_t(ieig,ik)) ehob(ik) = eko_t(ieig,ik)
                end if
             end do firstKpoint

             if(iprioccup >= 2) write(nfout,'(" ---  ik = ",i3," ihob = ",i6," totw = ",f12.4," totch = ",f12.4)') &
                  & ik,ihob(ik),totw,totch
          end do
          do ik = 1+is, kv3_t/ndim_spinor
             if(ihob(ik) /= ihob(is)) then
                metalic = .true.
                exit
             end if
          end do
       end do
       deallocate(ihob)
    else
       do is = 1, 1
          totw = 0.d0
          do ieig=1,neg
             if(eko_t(ieig,is) < efermi) then
                if(ehob(is) < eko_t(ieig,is)) ehob(is) = eko_t(ieig,is)
             end if
          end do
       end do
    end if

    if(.not.metalic .and. ekmode == OFF) then
       allocate(occup_t(neg,kv3_t/ndim_spinor))

       occup_t = 0.d0
       do ik = 1, kv3, ndim_spinor
          if(map_k(ik) /= myrank_k) cycle
	  iksnl = ( ik -1 )/ndim_spinor + 1
          do ie = 1, neg
             if(map_e(ie) /= myrank_e) cycle
             occup_t(ie,iksnl) = occup_l(map_z(ie),ik)
          end do
       end do
       if(npes > 1) then
          allocate(occup_mpi( neg, kv3_t/ndim_spinor ))
          call mpi_allreduce( occup_t, occup_mpi, neg*kv3/ndim_spinor, &
	       &              mpi_double_precision, mpi_sum,MPI_CommGroup, ierr )
          occup_t = occup_mpi
          deallocate(occup_mpi)
       end if

       Kpoint2: do ik = 1 ,kv3_t /ndim_spinor
          ihob0 = 0; vbm0 = -999.d+20
          findvbm: do ieig=1, neg   ! finding Valence Band Maximum
             if(occup_t(ieig,ik) > DELTA10) then
                if(vbm0 < eko_t(ieig,ik)) then
                   vbm0 = eko_t(ieig,ik)
                   ihob0 = ieig
                end if
             end if
          end do findvbm

! -------------------------- ktDEBUG ---------------------- 20121027
!          weight0 = kv3 /ndim_spinor *qwgt(ik)
!
          iktmp = ndim_spinor*(ik-1) +1
          weight0 = kv3 /ndim_spinor *qwgt( iktmp )
! -------------------------- ktDEBUG ---------------------- 20121027
!
          if(iprioccup >= 2) then
             write(nfout,'(" --- ihob0 = ",i8)') ik,ihob0
             if(ihob0 >= 1) write(nfout,'(" --- ik = ",i4," ihob0 = ",i8, &
	   &              " vbm0, occup = " ,2f10.5," weight0=",f8.4)') &
	   &                ndim_spinor *ik -1, ihob0, vbm0, &
	   &                occup_t(ihob0,ik), weight0
          end if
          if(ihob0 >= 1 .and. ihob0 <= neg) then
             if(occup_t(ihob0,ik) > DELTA10 .and. occup_t(ihob0,ik) < weight0-DELTA10) then
                metalic = .true.
                exit Kpoint2
             end if
          end if
       end do Kpoint2
       deallocate(occup_t)
    end if

    way_of_smearing = m_CtrlP_way_of_smearing()

    if(.not.metalic .or. way_of_smearing == Fermi_Dirac ) then
       cbm = 1.0D10
       Kpoint: do ik= 1, kv3_t /ndim_spinor
          findcbm: do ieig = 1, neg
#if 0
             if(eko_t(ieig,ik) > efermi) then
                cbm = eko_t(ieig,ik)
                exit findcbm
             end if
#else
             if(eko_t(ieig,ik) > efermi) then
                cbm = min(eko_t(ieig,ik),cbm)
                exit findcbm
             end if
#endif
          end do findcbm

          if((way_of_smearing==PARABOLIC .and. cbm < efermi + width/2) &
               & .or. (way_of_smearing==TETRAHEDRON .and. cbm < efermi+DELTA10) &
! =================== KT_add ================ 13.0E
               & .or. (way_of_smearing==Fermi_Dirac .and. cbm < efermi +width/2 ) &
               & .or. (way_of_smearing==LOWEST_AT_EACH_KPT .and. cbm < efermi ) &
! =========================================== 13.0E
               & .or. (way_of_smearing==COLD .and. cbm < efermi + width/2)) then

             metalic = .true.
             if(iprioccup >= 2) then
                write(nfout,'(" --- cbm = ",f8.4," ik = ",i8)') cbm, ndim_spinor *ik -1
             end if
             exit Kpoint
          end if
       end do Kpoint
    end if

    if(metalic) then
       if(iprioccup>=1) write(nfout,'(" --- The system is metalic ---")')
! --------------------------- added by T. Hamada -------------------------
       system_is_not_metallic = .false.
! ------------------------------------------------------------------------
       if ( way_of_smearing == Fermi_Dirac ) then
          emax = maxval(ehob)
          if(iprioccup>=1) then
             write(nfout,'(A,f12.8)') " --- The highest band energy below Ef = ",emax
             write(nfout,'(A,f12.8)') "         lowest  band energy above Ef = ",cbm
             write(nfout,'(A,f12.8,A,f12.8,A)') &
                  &                   "                           Difference = ", &
                  &           (cbm -emax), " Ha   ( ", (cbm -emax)*Hartree, " eV )"
             write(nfout,*)
          endif
       endif
    else
       if(iprioccup>=1) write(nfout,'(" --- The system is insulating or semiconducting ---")')
       emax = maxval(ehob)
#if 0
       if(iprioccup>=1) write(nfout,'(" --- The highest occupied band energy (=vbm) = ",f12.8)') emax
#else
       if(iprioccup>=1) then
          write(nfout,'(A,f12.8)') " --- The highest  occupied band energy (=vbm) = ",emax
          write(nfout,'(A,f12.8)') "         lowest unoccupied band energy (=cbm) = ",cbm
          write(nfout,'(A,f12.8,A,f12.8,A)') &
               &           "                           estimated band gap = ", &
               &           (cbm -emax), " Ha   ( ", (cbm -emax)*Hartree, " eV )"
          write(nfout,*)
       endif
#endif
       vbm = emax
! -------------------------- added by T. Hamada --------------------------
       system_is_not_metallic = .true.
! ------------------------------------------------------------------------
    end if

    check_if_metalic_flag = .true.
! ------------------------- added by T. Hamada ---------------------------
    system_is_already_determined = .true.
!-------------------------------------------------------------------------a

    deallocate(ehob)
    deallocate(eko_t)

  end subroutine check_if_metalic_noncl
!================================================================ 11.0


  subroutine m_ESoc_set_nEwindows_pc_ek(nfout,nEwindows,eorig)
    integer, intent(in) :: nfout
    integer, intent(out) :: nEwindows
    real(kind=DP), intent(in), optional :: eorig
    real(kind=DP) :: emin, emax, emin_rev, emax_rev, e0
    integer :: nvb_windows, ncb_windows, iw, ie, ik

    emin = minval(eko_ek)
    emax = maxval(eko_ek)

    if(metalic_system)  e_origin = efermi
    if(.not.metalic_system) e_origin = vbm + 1.d-10
    if(present(eorig)) e_origin = eorig
    if(ipri >= 1) then
       write(nfout,'(" !pc e_origin = ",f9.6," <<m_ESoc_set_nEwindows_pc>>")') e_origin
       write(nfout,'(" !pc emin, emax = ",2f9.6)') emin, emax
    end if

    emin_rev = emin - e_origin
    emax_rev = emax - e_origin
    if(partial_charge_Emin > emin_rev) emin_rev = partial_charge_Emin
    if(partial_charge_Emax < emax_rev) emax_rev = partial_charge_Emax

    e0 = 0.d0
    if(emax_rev < 0.0) e0 = emax_rev
    nvb_windows = (e0-emin_rev)/partial_charge_deltaE + 1
    if(nvb_windows < 0) nvb_windows = 0
    emin_rev = e0 - nvb_windows*partial_charge_deltaE

    e0 = 0.d0
    if(emin_rev > 0.0) e0 = emin_rev
    ncb_windows = (emax_rev - e0)/partial_charge_deltaE + 1
    if(ncb_windows < 0) ncb_windows = 0

    nEwindows_pc = nvb_windows + ncb_windows
    nEwindows = nEwindows_pc
    allocate(Ewindows_partial_charge(nEwindows_pc,2))
    Ewindows_partial_charge(1,1) = emin_rev + e_origin
    do iw = 2, nEwindows_pc
       Ewindows_partial_charge(iw,1) = Ewindows_partial_charge(iw-1,1) + partial_charge_deltaE
       Ewindows_partial_charge(iw-1,2) = Ewindows_partial_charge(iw,1)
    end do
    Ewindows_partial_charge(nEwindows_pc,2) = Ewindows_partial_charge(nEwindows_pc,1) + partial_charge_deltaE

    allocate(if_elec_state(nEwindows_pc)); if_elec_state = NO

    do iw = 1, nEwindows_pc
       emin = Ewindows_partial_charge(iw,1); emax = Ewindows_partial_charge(iw,2)

       if ( noncol ) then
         Kloop: do ik = 1, kv3_ek/ndim_spinor
            do ie = 1, neg
               if(emin < eko_ek(ie,ik) .and. eko_ek(ie,ik) < emax) then
                  if_elec_state(iw) = YES
                  exit Kloop
               end if
            end do
         end do Kloop
       else
         Kloop2: do ik = 1, kv3_ek
            do ie = 1, neg
               if(emin < eko_ek(ie,ik) .and. eko_ek(ie,ik) < emax) then
                  if_elec_state(iw) = YES
                  exit Kloop2
               end if
            end do
         end do Kloop2
      endif
! =========================================================================== 11.0
    end do

    if(ipri >= 1) then
    nEwindows_pc = nvb_windows + ncb_windows
       write(nfout,'(" !pc nEwindows = ",i4,", nvb_windows = ",i4,", ncb_windows = ",i4 &
            & ," <<m_ESoc_set_nEwindows_pc>>")') nEwindows_pc, nvb_windows, ncb_windows
       write(nfout,'(" !pc",4x,"iw",2x,"if_elec_state",16x,"erange(hartree)",24x,"erange(eV)")')
       write(nfout,'(" !pc",29x,"(asis)",18x,"(shifted)",13x,"(shifted)")')
       do iw = 1, nEwindows_pc
          write(nfout,'(" !pc ",i5,7x,i1,7x,"(",2f10.6," )",1x,"(",2f10.6," )",1x,"(",2f10.6," )")') &
               &   iw, if_elec_state(iw) &
               & , Ewindows_partial_charge(iw,1),Ewindows_partial_charge(iw,2) &
               & , Ewindows_partial_charge(iw,1)-e_origin &
               & , Ewindows_partial_charge(iw,2)-e_origin &
               & , (Ewindows_partial_charge(iw,1)-e_origin)*HARTREE &
               & , (Ewindows_partial_charge(iw,2)-e_origin)*HARTREE
       end do
    end if

  end subroutine m_ESoc_set_nEwindows_pc_ek

  subroutine m_ESoc_set_nEwindows_pc(nfout,nEwindows,eorig)
    integer, intent(in) :: nfout
    integer, intent(out) :: nEwindows
    real(kind=DP), intent(in), optional :: eorig
    real(kind=DP) :: emin, emax, emin_rev, emax_rev, e0
    integer :: nvb_windows, ncb_windows, iw, ie, ik
    real(kind=DP),allocatable, dimension(:,:) :: eko_tmp

! ================================ modified by K. Tagami ============ 11.0
!    allocate(eko_tmp(neg,kv3)); eko_tmp = 0.d0
!    call find_emin_emax(emin,emax,eko_tmp)

    if ( noncol ) then
      allocate(eko_tmp(neg,kv3/ndim_spinor)); eko_tmp = 0.d0
      call find_emin_emax_noncl( emin,emax,eko_tmp )
    else
      allocate(eko_tmp(neg,kv3)); eko_tmp = 0.d0
      call find_emin_emax(emin,emax,eko_tmp)
    endif
! ================================================================== 11.0

    if(metalic_system)  e_origin = efermi
    if(.not.metalic_system) e_origin = vbm + 1.d-10
    if(present(eorig)) e_origin = eorig
    if(ipri >= 1) then
       write(nfout,'(" !pc e_origin = ",f9.6," <<m_ESoc_set_nEwindows_pc>>")') e_origin
       write(nfout,'(" !pc emin, emax = ",2f9.6)') emin, emax
    end if

    emin_rev = emin - e_origin
    emax_rev = emax - e_origin
    if(partial_charge_Emin > emin_rev) emin_rev = partial_charge_Emin
    if(partial_charge_Emax < emax_rev) emax_rev = partial_charge_Emax

    e0 = 0.d0
    if(emax_rev < 0.0) e0 = emax_rev
    nvb_windows = (e0-emin_rev)/partial_charge_deltaE + 1
    if(nvb_windows < 0) nvb_windows = 0
    emin_rev = e0 - nvb_windows*partial_charge_deltaE

    e0 = 0.d0
    if(emin_rev > 0.0) e0 = emin_rev
    ncb_windows = (emax_rev - e0)/partial_charge_deltaE + 1
    if(ncb_windows < 0) ncb_windows = 0

    nEwindows_pc = nvb_windows + ncb_windows
    nEwindows = nEwindows_pc
    allocate(Ewindows_partial_charge(nEwindows_pc,2))
    Ewindows_partial_charge(1,1) = emin_rev + e_origin
    do iw = 2, nEwindows_pc
       Ewindows_partial_charge(iw,1) = Ewindows_partial_charge(iw-1,1) + partial_charge_deltaE
       Ewindows_partial_charge(iw-1,2) = Ewindows_partial_charge(iw,1)
    end do
    Ewindows_partial_charge(nEwindows_pc,2) = Ewindows_partial_charge(nEwindows_pc,1) + partial_charge_deltaE

    allocate(if_elec_state(nEwindows_pc)); if_elec_state = NO

    do iw = 1, nEwindows_pc
       emin = Ewindows_partial_charge(iw,1); emax = Ewindows_partial_charge(iw,2)

! ========================================== modified by K. Tagami ========== 11.0
!       Kloop: do ik = 1, kv3
!          do ie = 1, neg
!             if(emin < eko_tmp(ie,ik) .and. eko_tmp(ie,ik) < emax) then
!                if_elec_state(iw) = YES
!                exit Kloop
!             end if
!          end do
!       end do Kloop
!
       if ( noncol ) then
         Kloop: do ik = 1, kv3/ndim_spinor
            do ie = 1, neg
               if(emin < eko_tmp(ie,ik) .and. eko_tmp(ie,ik) < emax) then
                  if_elec_state(iw) = YES
                  exit Kloop
               end if
            end do
         end do Kloop
       else
         Kloop2: do ik = 1, kv3
            do ie = 1, neg
               if(emin < eko_tmp(ie,ik) .and. eko_tmp(ie,ik) < emax) then
                  if_elec_state(iw) = YES
                  exit Kloop2
               end if
            end do
         end do Kloop2
      endif
! =========================================================================== 11.0
    end do

    if(ipri >= 1) then
    nEwindows_pc = nvb_windows + ncb_windows
       write(nfout,'(" !pc nEwindows = ",i4,", nvb_windows = ",i4,", ncb_windows = ",i4 &
            & ," <<m_ESoc_set_nEwindows_pc>>")') nEwindows_pc, nvb_windows, ncb_windows
       write(nfout,'(" !pc",4x,"iw",2x,"if_elec_state",16x,"erange(hartree)",24x,"erange(eV)")')
       write(nfout,'(" !pc",29x,"(asis)",18x,"(shifted)",13x,"(shifted)")')
       do iw = 1, nEwindows_pc
          write(nfout,'(" !pc ",i5,7x,i1,7x,"(",2f10.6," )",1x,"(",2f10.6," )",1x,"(",2f10.6," )")') &
               &   iw, if_elec_state(iw) &
               & , Ewindows_partial_charge(iw,1),Ewindows_partial_charge(iw,2) &
               & , Ewindows_partial_charge(iw,1)-e_origin &
               & , Ewindows_partial_charge(iw,2)-e_origin &
               & , (Ewindows_partial_charge(iw,1)-e_origin)*HARTREE &
               & , (Ewindows_partial_charge(iw,2)-e_origin)*HARTREE
       end do
    end if

    deallocate(eko_tmp)

  end subroutine m_ESoc_set_nEwindows_pc

  subroutine find_emin_emax(emin,emax,eko_tmp)
     use m_Parallelization, only : mpi_kg_world
     real(kind=DP), intent(out) :: emin, emax
     real(kind=DP), intent(out), dimension(neg,kv3) :: eko_tmp

     integer :: ik, ie, is
     real(kind=DP),allocatable, dimension(:,:) :: temp_mpi ! MPI

     if(npes >=2 ) then
        allocate(temp_mpi(neg,kv3)); temp_mpi = 0.d0
        do is = ista_spin, iend_spin
        do ik = is, kv3-nspin+is, nspin

          if(map_k(ik) /= myrank_k) cycle                 ! MPI
          do ie = 1, neg                                  ! MPI
             if(map_e(ie) /= myrank_e) cycle              ! MPI
             temp_mpi(ie,ik) = eko_l(map_z(ie),ik)        ! MPI
          end do                                          ! MPI
        end do                                             ! MPI
        end do
        call mpi_allreduce(temp_mpi,eko_tmp,neg*kv3,mpi_double_precision &  ! MPI
             &            ,mpi_sum,mpi_kg_world,ierr)                     ! MPI
        temp_mpi = eko_tmp
        call mpi_allreduce(temp_mpi,eko_tmp,neg*kv3,mpi_double_precision &  ! MPI
             &            ,mpi_sum,mpi_ge_world,ierr)                     ! MPI
        deallocate(temp_mpi )
     else
        eko_tmp = eko_l
     end if
     emin = minval(eko_tmp)
     emax = maxval(eko_tmp)
  end subroutine find_emin_emax

! ======================================= added by K. Tagami ============== 11.0
  subroutine find_emin_emax_noncl( emin,emax,eko_tmp )
     real(kind=DP), intent(out) :: emin, emax
     real(kind=DP), intent(out), dimension(neg,kv3/ndim_spinor) :: eko_tmp

     integer :: ik, ie, iksnl
     real(kind=DP),allocatable, dimension(:,:) :: temp_mpi ! MPI

     if(npes >=2 ) then
        allocate(temp_mpi(neg,kv3/ndim_spinor)); temp_mpi = 0.d0

        do ik = 1, kv3, ndim_spinor
          if(map_k(ik) /= myrank_k) cycle
	  iksnl = ( ik - 1 )/ndim_spinor + 1
          do ie = 1, neg
             if(map_e(ie) /= myrank_e) cycle
             temp_mpi(ie,iksnl) = eko_l(map_z(ie),ik)
          end do
        end do
        call mpi_allreduce( temp_mpi, eko_tmp, neg*kv3/ndim_spinor, &
	  &                 mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
        deallocate(temp_mpi )
     else
        do ik = 1, kv3, ndim_spinor
	  iksnl = ( ik -1 )/ndim_spinor + 1
          eko_tmp(:,iksnl) = eko_l(:,ik)
	end do
     end if
     emin = minval(eko_tmp);   emax = maxval(eko_tmp)
  end subroutine find_emin_emax_noncl
! ========================================================================= 11.0

  subroutine m_ESoc_keep_occup()
     integer :: ie, ik

     allocate(occup_tmp(np_e, ista_k:iend_k))
     do ik = ista_k, iend_k
        do ie = 1, np_e
           occup_tmp(ie,ik) = occup_l(ie,ik)
        end do
     end do
  end subroutine m_ESoc_keep_occup

  subroutine m_ESoc_retrieve_occup()
     integer :: ie, ik

     do ik = ista_k, iend_k
        do ie = 1, np_e
           occup_l(ie,ik) = occup_tmp(ie,ik)
        end do
     end do
     deallocate(occup_tmp)
  end subroutine m_ESoc_retrieve_occup

  subroutine m_ESoc_occup_under_ef(nfout, nk)
    use m_Electronic_Structure, only : efermi
    integer, intent(in) :: nfout
    integer, intent(in), optional :: nk
    integer :: ik, ie
    real(kind=DP) :: e, weight, eo
    integer :: noff

    noff = 1
    if(present(nk)) noff = nk
    eo = efermi
    if(.not.metalic_system) eo = vbm + 1.d-10
    do ik = ista_k, iend_k    ! MPI
       do ie = 1, np_e        ! MPI
          if(ekmode==OFF) then
             e = eko_l(ie,ik)
          else
             e = eko_ek(ista_e-1+ie,noff+ik-1)
          endif
          if(e < eo) then
             weight = 1.d0
          else
             weight = 0.d0
          end if
! =========================================== modified by K. Tagami =========== 11.0
!         occup_l(ie,ik) = weight*wspin*kv3*qwgt(ik)
!
          if ( noncol ) then
            occup_l(ie,ik) = weight * kv3/ndim_spinor *qwgt(ik)
          else
            occup_l(ie,ik) = weight * kv3 *qwgt(ik)
          endif
! ============================================================================ 11.0
       end do
    end do
  end subroutine m_ESoc_occup_under_ef

  subroutine m_ESoc_substitute_occup(nfout,iw,nk)
    use m_Parallelization, only : mpi_kg_world
    integer, intent(in) :: nfout,iw
    integer, intent(in), optional :: nk

    integer :: ie, ik, is, ipripe0, i, iset, j, is0, ie0
    integer, parameter :: NCLMN = 8
    real(kind=DP) :: wspin = 1.d0, weight, emin, emax, e
    real(kind=DP), allocatable, dimension(:,:) :: occup_mpi, temp_mpi, eko_mpi
    integer :: noff

    noff = 1
    if(present(nk)) noff = nk

    emin = Ewindows_partial_charge(iw,1)
    emax = Ewindows_partial_charge(iw,2)

    do is = ista_spin, iend_spin
    do ik = is, kv3-nspin+is, nspin
       if(map_k(ik) /= myrank_k) cycle
       do ie = 1, np_e        ! MPI
          if(ekmode==OFF) then
             e = eko_l(ie,ik)
          else
             e = eko_ek(ista_e-1+ie,noff+ik-1)
          endif
          if(e > emin .and. e < emax) then
             weight = 1.d0
          else
             weight = 0.d0
          end if
! =========================================== modified by K. Tagami =========== 11.0
!         occup_l(ie,ik) = weight*wspin*kv3*qwgt(ik)
!
          if ( noncol ) then
            occup_l(ie,ik) = weight *wspin *kv3/ndim_spinor *qwgt(ik)
          else
            occup_l(ie,ik) = weight *wspin *kv3 *qwgt(ik)
          endif
! ============================================================================ 11.0

       end do
    end do
    end do

    if(npes > 1) then
       if(mype == 0) ipripe0 = iprioccup
       call mpi_bcast(ipripe0, 1,mpi_integer,0,MPI_CommGroup,ierr)
    else
       ipripe0 = iprioccup
    end if

    if(ipripe0 >=1 ) then
       allocate(occup_mpi(neg,kv3)); occup_mpi = 0.d0
       if(npes >= 2) then
          allocate(temp_mpi(neg,kv3)); temp_mpi = 0.d0
          do is = ista_spin, iend_spin
          do ik = is, kv3-nspin+is, nspin
             if(map_k(ik) /= myrank_k) cycle
             do ie = 1, neg
                if(map_e(ie) /= myrank_e) cycle
                temp_mpi(ie,ik) = occup_l(map_z(ie),ik)
             end do
          end do
          end do
          call mpi_allreduce(temp_mpi,occup_mpi,neg*kv3,mpi_double_precision &
               &            ,mpi_sum, mpi_kg_world, ierr)
          temp_mpi = occup_mpi
          call mpi_allreduce(temp_mpi,occup_mpi,neg*kv3,mpi_double_precision &
               &            ,mpi_sum, mpi_ge_world, ierr)
          deallocate(temp_mpi)
       else
          occup_mpi = occup_l
       end if

       allocate(eko_mpi(neg,kv3)); eko_mpi = 0.d0
       if(npes >= 2) then
          allocate(temp_mpi(neg,kv3)); temp_mpi = 0.d0
          do is = ista_spin, iend_spin
          do ik = is, kv3-nspin+is, nspin
             if(map_k(ik) /= myrank_k) cycle
             do ie = 1, neg
                if(map_e(ie) /= myrank_e) cycle
                temp_mpi(ie,ik) = eko_l(map_z(ie),ik)
             end do
          end do
          end do
          call mpi_allreduce(temp_mpi,eko_mpi,neg*kv3,mpi_double_precision &
               &            ,mpi_sum, mpi_kg_world, ierr)
          temp_mpi = eko_mpi
          call mpi_allreduce(temp_mpi,eko_mpi,neg*kv3,mpi_double_precision &
               &            ,mpi_sum, mpi_ge_world, ierr)
          deallocate(temp_mpi)
       else
          eko_mpi = eko_l
       end if

    end if
    if(iprioccup >= 1) then
       write(nfout,'(" !pc iw = ",i5," <<m_ESoc_substitute_occup>>")') iw

! ======================================== modified by K. Tagami ========= 11.0
!       do ik = 1, kv3
       do ik = 1, kv3, ndim_spinor
! ======================================================================= 11.0

          is = neg+1; ie = 1-1
          Isloop: do i = 1, neg
             if(occup_mpi(i,ik) > 1.d-20) then
                is = i
                exit Isloop
             end if
          end do Isloop
          Ieloop: do i = neg, 1, -1
             if(occup_mpi(i,ik) > 1.d-20) then
                ie = i
                exit Ieloop
             end if
          end do Ieloop
          if(is > 0 .and. is <= neg .and. ie > 0 .and. ie <= neg .and. is <= ie) then
             iset = (ie-is+1+(NCLMN-1))/NCLMN
             is0 = is
             write(nfout,'(" !pc ( ik = ",i5," )")') ik
             do j = 1, iset
                ie0 = min(ie,is0+NCLMN-1)
                write(nfout,'(" !pc nband ",4x,8i10)') (i,i=is0,ie0)
                write(nfout,'(" !pc eko  ",5x,8f10.6)') (eko_mpi(i,ik),i=is0, ie0)
                write(nfout,'(" !pc eko(eV)",3x,8f10.6)') ((eko_mpi(i,ik)-e_origin)*HARTREE,i=is0, ie0)
                write(nfout,'(" !pc occup",5x,8f10.6)') (occup_mpi(i,ik),i=is0, ie0)
                is0 = is0 + NCLMN
             end do
          end if
       end do
    end if
    if(ipripe0 >=1 ) then
       deallocate(occup_mpi)
       deallocate(eko_mpi)
    end if

  end subroutine m_ESoc_substitute_occup

  subroutine m_ESoc_free_nEwindows()
    if(allocated(Ewindows_partial_charge)) deallocate(Ewindows_partial_charge)
    if(allocated(if_elec_state)) deallocate(if_elec_state)
  end subroutine m_ESoc_free_nEwindows

  integer function m_ESoc_if_elec_state(nfout,iw,emin,emax)
    integer, intent(in) :: nfout,iw
    real(kind=DP), intent(out) :: emin, emax

    if(iw < 0 .or. iw > nEwindows_pc) then
       write(nfout,'(" iw is out of the proper range. iw = ",i5," <<m_ESoc_if_elec_state>>")') iw
       m_ESoc_if_elec_state = NO
       return
    else
       emin = Ewindows_partial_charge(iw,1)
       emax = Ewindows_partial_charge(iw,2)
       m_ESoc_if_elec_state = if_elec_state(iw)
    end if
  end function m_ESoc_if_elec_state

  subroutine check_occupation(nfout,ispin,totch)
    integer, intent(in) :: nfout,ispin
    real(kind=DP),intent(in),dimension(*) :: totch
    integer  :: i
#ifdef _CHECK_OCCUPATION_STOP_
    integer :: ierror, errorcode
#endif

    if(mype==0) then
       do i =1, ispin
          if(dble(neg)-totch(i)*0.5d0 < 1.0) then
             if(iprioccup>=1) &
     &       write(nfout,'(" WARNING! neg may not be sufficient : &
     &    (neg = ",i8," ) - (totch(",i2,")*0.5 = ",f12.4,") < 1.0")') neg,i,totch(i)*0.5
          end if
       end do
    end if
#ifdef _CHECK_OCCUPATION_STOP_
    do i = 1, ispin
       if(dble(neg)-totch(i)*0.5d0 < 0.0) then
          if(iprioccup>=1) write(nfout,'(" STOP! neg is not sufficient : &
        & (neg = ",i8," ) - (totch(",i2,")*0.5 = ",f12.4,") < 0.0")') neg,i,totch(i)*0.5
          call mpi_abort(MPI_CommGroup,errorcode,ierror)
          call mpi_finalize(ierror)
       end if
    end do
#endif
  end subroutine check_occupation

! =========================== added by K. Tagami ================= 11.0
  subroutine check_occupation_noncl( nfout, totch )
    integer, intent(in) :: nfout
    real(kind=DP),intent(in),dimension(*) :: totch
    integer  :: i
#ifdef _CHECK_OCCUPATION_STOP_
    integer :: ierror, errorcode
#endif

    if(mype==0) then
       do i =1, 1
          if(dble(neg)-totch(i) < 1.0) then
             if(iprioccup>=1) &
     &       write(nfout,'(" WARNING! neg may not be sufficient : &
     &    (neg = ",i8," ) - (totch(",i2,") = ",f12.4,") < 1.0")') neg,i,totch(i)
          end if
       end do
    end if
#ifdef _CHECK_OCCUPATION_STOP_
    do i = 1, 1
       if(dble(neg)-totch(i) < 0.0) then
          if(iprioccup>=1) write(nfout,'(" STOP! neg is not sufficient : &
        & (neg = ",i8," ) - (totch(",i2,") = ",f12.4,") < 0.0")') neg,i,totch(i)
          call mpi_abort(MPI_CommGroup,errorcode,ierror)
          call mpi_finalize(ierror)
       end if
    end do
#endif
  end subroutine check_occupation_noncl

  subroutine m_ESoc_fermi_parabolic_3D(nfout)
!
!  An if-block of "if(imag == FERRO .and. sw_fix_total_spin == YES .and. nspin == 2) then",
! which defines efermi_spin and occup_l with fixing total spin, is coded.
!                        by T. Yamasaki   April/10/2007
!
    use m_Parallelization,      only : mpi_kg_world, mpi_ke_world, neg_g
    use m_Electronic_Structure, only : nrvf_ordr


    integer, intent(in) :: nfout

    integer             :: jcount, ik, ie, is, lrank, iadd
    Real(kind=DP)       :: emin, emax, e1, e2, tot, totch_spin(2), totch_spin0(2)
    real(kind=DP), parameter :: DELTA_TOTCH = 1.d-10
    integer,       parameter :: MAXITR = 100000
    real(kind=DP),allocatable, dimension(:,:) :: temp_mpi, eko_mpi, occup_mpi  ! MPI

    integer :: id_sname = -1, id_sname2 = -1

!   call tstatc0_begin('barrier(m_Eoc_fermi_parabolic) ',id_sname2)
!   if(npes > 1) call mpi_barrier(MPI_CommGroup,ierr)
!   call tstatc0_end(id_sname2)
                                                  __TIMER_SUB_START(703)

    call tstatc0_begin('m_ESoc_fermi_parabolic_3D ', id_sname)

    call check_totch(nfout)           ! -(m_ES_occup) totch will be checked

! === DEBUG by tkato 2011/12/07 ================================================
! === Initialization with zero is necessary when af = 1 or nspin = 2! ==========
    allocate(eko_mpi  (neg,kv3)); eko_mpi = 0.d0
    allocate(occup_mpi(neg,kv3)); occup_mpi = 0.d0
! ==============================================================================

    !  MPI Gather eko_l
    if(npes > 1) then
!       if(nrank_k > 1) then
          allocate(temp_mpi (neg,kv3)); temp_mpi = 0.d0
                                                  __TIMER_DO_START(800)
          do is = ista_spin, iend_spin
          !do ik = 1, kv3
          do ik = is, kv3-nspin+is, nspin
             if(map_k(ik) /= myrank_k) cycle
!            do ie = 1, neg
!               if(map_e(ie) /= myrank_e) cycle
!               temp_mpi(ie,ik) = eko_l(map_z(ie),ik)
!            end do
             do ie = 1, np_e
                temp_mpi(ista_e-1+ie,ik) = eko_l(ie,ik)
             end do
          end do
          end do
                                                  __TIMER_DO_STOP(800)
                                                  __TIMER_COMM_START_w_BARRIER(mpi_kg_world,801)
          call mpi_allreduce(temp_mpi,eko_mpi,neg*kv3,mpi_double_precision &  ! MPI
         &                   ,mpi_sum,mpi_kg_world,ierr)                     ! MPI
! ==== DEBUG by tkato 2012/04/04 ===============================================
          temp_mpi = eko_mpi
          call mpi_allreduce(temp_mpi,eko_mpi,neg*kv3,mpi_double_precision & ! MPI
         &                   ,mpi_sum,mpi_ge_world,ierr)                     ! MPI
! ==============================================================================
                                                  __TIMER_COMM_STOP(801)
          deallocate(temp_mpi )
!       else
!          allocate(temp_mpi(mp_e*kv3,0:nrank_e-1))                        ! MPI
!                                                  __TIMER_COMM_START_w_BARRIER(mpi_kg_world,802)
!          call mpi_allgather(eko_l, np_e*kv3, mpi_double_precision &
!         &                 , temp_mpi,  mp_e*kv3, mpi_double_precision, mpi_kg_world, ierr)
!                                                  __TIMER_COMM_STOP(802)
!                                                  __TIMER_DO_START(803)
!          do lrank = 0, nrank_e-1
!             do ik = 1, kv3
!                iadd = nel_e(lrank)*(ik-1)-nis_e(lrank)+1
!                do is = nis_e(lrank), nie_e(lrank)
!                   eko_mpi(is,ik) = temp_mpi(iadd+is,lrank)
!                end do
!             end do
!          end do
!                                                  __TIMER_DO_STOP(803)
!          deallocate(temp_mpi)
!       end if
    else
       eko_mpi = eko_l
    end if
    emin = minval(eko_mpi)
    emax = maxval(eko_mpi)

    efermi = emax
    e1 = emin - dabs(DELTA_FermiSearchRange)
    e2 = emax + dabs(DELTA_FermiSearchRange)

    jcount = 1
    occupied_ch_equals_totch : do

       call get_occup_l_and_tot(1,1,efermi,.true.)          ! -(contained here) ->(occup_l,tot)
!          ~~~~~~~~~~~~~~~~~~~~
       if(jcount == 1 .and. tot < totch) &
            call wd_fermi_error1(nfout,emin,emax,tot,totch) ! -(b_Fermi)
       if(dabs(tot - totch) < DELTA_TOTCH) exit occupied_ch_equals_totch
       if( tot < totch) then
          e1 = efermi; efermi = efermi + (e2-efermi)/2
       else
          e2 = efermi; efermi = efermi + (e1-efermi)/2
       end if
       jcount = jcount + 1
       if(jcount > MAXITR) then
          call wd_fermi_error2 &   !-(b_Fermi)
               & (nfout,e1,e2,efermi,emin,emax,tot,totch,neg,MAXITR)
       end if
    end do occupied_ch_equals_totch
    if(iprioccup >= 2) &
         &write(nfout,'( "== Parabolic Broadening Method ==")')
#ifdef DBG_SPIN_FIX       
    if(iprioccup >= 1) &
#else
    if(iprioccup >= 2) &
#endif
         &write(nfout,'( " efermi = ", f10.4, " tot = ", f10.4)') efermi, tot
    call get_entropic_term(band_entropy)
    if(iprioccup >= 2) &
         &write(nfout,'( " band_entropy = ", f20.14)') band_entropy

    if(imag == FERRO .and. sw_fix_total_spin == YES .and. nspin == 2) then
#ifdef DBG_SPIN_FIX       
       if(iprioccup>=1) write(nfout,'(" total_spin = ",f8.4," <<m_ESoc_fermi_parabolic_3D>>")') total_spin
#endif
       if(total_spin > totch) then
          write(nfout,'( " !WARNING total_spin should be smaller than totch <<m_ESoc_fermi_parabolic>>")')

          totch_spin(1) = totch
          totch_spin(2) = totch
       else
          totch_spin(1) = (totch + total_spin)*1.0d0
          totch_spin(2) = (totch - total_spin)*1.0d0
       end if

       do is = 1, nspin
          emin = minval(eko_mpi)
          emax = maxval(eko_mpi)

          efermi_spin(is) = emax
          e1 = emin - dabs(DELTA_FermiSearchRange)
          e2 = emax + dabs(DELTA_FermiSearchRange)

          jcount = 1
          occupied_ch_equals_totch2 : do

             call get_occup_l_and_tot(nspin,is,efermi_spin(is),.true.)     ! -(contained here) ->(occup_mpi,tot)
             !    ~~~~~~~~~~~~~~~~~~~
             if(jcount == 1 .and. tot < totch_spin(is)) &
                  call wd_fermi_error1(nfout,emin,emax,tot,totch) ! -(b_Fermi)
             if(dabs(tot - totch_spin(is)) < DELTA_TOTCH) exit occupied_ch_equals_totch2
             if( tot < totch_spin(is)) then
                e1 = efermi_spin(is); efermi_spin(is) = (efermi_spin(is) + e2)/2
             else
                e2 = efermi_spin(is); efermi_spin(is) = (efermi_spin(is) + e1)/2
             end if
             jcount = jcount + 1
             if(jcount > MAXITR) then
                call wd_fermi_error2 &   !-(b_Fermi)
                     & (nfout,e1,e2,efermi_spin(is),emin,emax,tot,totch_spin(is),neg,MAXITR)
             end if
          end do occupied_ch_equals_totch2
#ifdef DBG_SPIN_FIX
          if(iprioccup >= 1 .and. nspin==2) then
             write(nfout,'(" totch_spin = ",2f16.8)') totch_spin(1:2)
          end if
#endif
#ifdef DBG_SPIN_FIX
          if(iprioccup >= 1) &
#else
          if(iprioccup >= 2) &
#endif
               &write(nfout,'( "== Parabolic Broadening Method ==")')
#ifdef DBG_SPIN_FIX
          if(iprioccup >= 1) &
#else
          if(iprioccup >= 2) &
#endif
               &write(nfout,'( " efermi_spin(",i2," ) = ", f10.4, ", tot/2 = ", f10.4, " jcount = ",i4)') &
               & is,efermi_spin(is), tot/2.0, jcount
       end do
       call get_total_spin0
       efermi = (efermi_spin(1)+efermi_spin(2))*0.5d0
       if(iprioccup>=1) write(nfout,'(a,f10.4)') ' == new efermi = ',efermi
    else if(imag == FERRO .and. sw_fix_total_spin == NO .and. nspin == 2) then
       call get_total_spin0
    end if

!   do ik = ista_k, iend_k                              ! MPI
!      do ie = 1, neg                                   ! MPI
!         if(map_e(ie) == myrank_e) then                ! MPI
!            occup_l(map_z(ie),ik) = occup_mpi(ie,ik)   ! MPI
!         end if                                        ! MPI
!      end do                                           ! MPI
!   end do                                              ! MPI
!    do ik = ista_k, iend_k                              ! MPI
    do is = ista_spin, iend_spin
    do ik = is, kv3-nspin+is, nspin
       if(map_k(ik) /= myrank_k) cycle                 ! MPI
!      do ie = 1, np_e                                  ! MPI
!         occup_l(ie,ik) = occup_mpi(neg_g(ie),ik)      ! MPI
       do ie = 1, np_e                                  ! MPI
          occup_l(ie,ik) = occup_mpi(ista_e-1+ie,ik)    ! MPI
       end do                                           ! MPI
    end do                                              ! MPI
    end do

    deallocate(eko_mpi  )
    deallocate(occup_mpi)

    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(703)
  contains
    subroutine get_total_spin0
                                                  __TIMER_SUB_START(708)
      do is = 1, nspin
         call get_occup_l_and_tot(nspin,is,efermi,.false.)
         totch_spin0(is) = tot
      end do
      total_spin0 = (totch_spin0(1) - totch_spin0(2))*0.5d0
      if(iprioccup >= PRINTOUTLEVEL_TOTCH_SPIN) then
         if(sw_fix_total_spin == YES) then
            call wd_efermi_and_total_spin0_Plus(nfout,total_spin0,totch_spin0,totch_spin)
         else
            call wd_efermi_and_total_spin0(nfout,total_spin0,totch_spin0)
         end if
      end if
                                                  __TIMER_SUB_STOP(708)
    end subroutine get_total_spin0

    subroutine get_occup_l_and_tot(nspin,is,efermi,update_occ)
      integer, intent(in) :: nspin, is
      real(kind=DP), intent(in) :: efermi
      logical, intent(in) :: update_occ
      integer       :: k, i,iupdown
      real(kind=DP) :: wspin = 1.d0, e, dos, weight, tmp
                                                  __TIMER_SUB_START(704)
      tot = 0.d0
                                                  __TIMER_SUB_START(705)
                                                  __TIMER_SUB_STOP(705)

                                                  __TIMER_DO_START(804)
      do k = is, kv3, max(nspin,af+1)
         do i = 1, neg
            e = eko_mpi(i,k)
            call width2(e,efermi,width,dos,weight)  ! -(b_Fermi)
            if(sw_manual_occupation==ON.and..not.update_occ)then
               if(band_index(i)>0)then
                  iupdown = 1
                  if(nspin>1.and.mod(k,2)==0) iupdown = 2
                  weight = occ_ext(band_index(i),iupdown)
               endif
            endif
! ====================================================================== 12.1

            tmp = weight*wspin*kv3*qwgt(k)
            if(update_occ) occup_mpi(i,k) = tmp
            !tot = tot + 2*occup_mpi(i,k)
            tot = tot + 2*tmp
         end do
      end do
                                                  __TIMER_DO_STOP(804)
      if(af == 1) then
         tot = tot/kv3 * (af+1)
      else
         tot = tot/(kv3/nspin)
      end if
      if(iprioccup >= 2) write(nfout,'(" efermi, tot = ",d20.12,d20.12)') efermi,tot
                                                  __TIMER_SUB_STOP(704)
    end subroutine get_occup_l_and_tot


    subroutine get_entropic_term(s)
      real(kind=DP), intent(out) :: s

      integer       :: k, i
      real(kind=DP) :: e, entropy
                                                  __TIMER_SUB_START(706)
      s = 0.d0
                                                  __TIMER_DO_START(805)
      do k = 1, kv3, af+1
         do i = 1, neg
            e = eko_mpi(i,k)
            call get_entropy(e,efermi,width,entropy)
            s = s + 2.d0*entropy*qwgt(k)
         end do
      end do
                                                  __TIMER_DO_STOP(805)
      s = s * (af+1)
                                                  __TIMER_SUB_STOP(706)
    end subroutine get_entropic_term

    subroutine get_entropy(e,efermi,width,entropy)
      real(kind=DP), intent(in) :: e,efermi,width
      real(kind=DP), intent(out) :: entropy

      real(kind=DP) :: xi,xj,xx
      real(kind=DP), parameter :: d16 = 1.d0/16.d0
      real(kind=DP), parameter :: d6 = 1.d0/6.d0
      real(kind=DP), parameter :: sd24 = 7.d0/24.d0
                                                  __TIMER_SUB_START(707)
      xi = (efermi-e)/width

      if(xi < -2.d0) then
         entropy = 0.d0
      else if(xi < -1.d0) then
         xj = xi+2.d0
         xx = xj*xj
         entropy = -xx**2*d16+xx*xj*d6
      else if(xi < 1.d0) then
         xx = 0.25d0*xi*xi
         entropy = -xx*(1.d0-xx) + sd24
      else if(xi < 2.d0) then
         xj = xi-2.d0
         xx = xj*xj
         entropy = -xx**2*d16-xx*xj*d6
      else
         entropy = 0.d0
      end if
                                                  __TIMER_SUB_STOP(707)
    end subroutine get_entropy

  end subroutine m_ESoc_fermi_parabolic_3D

!===============================================================================
  subroutine m_ESoc_fermi_tetrahedron_3D(nfout)
! tetra hedron method (Dr.Hamada-san's code)
!                        by Tsuyoshi Miyazaki '94.8.9
!
!  An if-block of "if(imag == FERRO .and. sw_fix_total_spin == YES .and. nspin == 2) then",
! which defines efermi_spin and occup_l with fixing total spin, is coded.
!                        by T. Yamasaki   April/10/2007
!
    use m_Parallelization,      only : mpi_kg_world, mpi_ke_world

    integer, intent(in) :: nfout
#ifndef NO_TETRAHEDRON
    !!$integer, parameter       :: idim = 3
    integer        :: neig,nengy,ispin,ip2,ik,ieig,ikee,nxx,nyy,nzz, ip, ip_mpi
    real(kind=DP)  :: efermi2,eval,totind
    real(kind=DP), pointer, dimension(:,:,:) :: eig2,occup2,eig2_mpi,occup2_mpi
    real(kind=DP), pointer, dimension(:) :: eawk,cdwk,cswk,cdos,cind,valud
    real(kind=DP) :: totch_spin(2), totch_spin0(2)
    integer             :: id_sname = -1
                                                  __TIMER_SUB_START(709)
    call tstatc0_begin('m_ESoc_fermi_tetrahedron_3D ', id_sname)

    allocate(eig2(np2,neg,nspin)); eig2 = 0.d0
    allocate(eig2_mpi(np2,neg,nspin)); eig2_mpi = 0.d0     ! MPI
    allocate(occup2(neg,np2,nspin)); occup2 = 0.d0
    allocate(occup2_mpi(neg,np2,nspin)); occup2_mpi = 0.d0 ! MPI
    allocate(eawk(np0)); eawk = 0.d0
    allocate(cdwk(np0)); cdwk = 0.d0
    allocate(cswk(np0)); cswk = 0.d0
    allocate(cdos(np2*neg)); cdos = 0.d0
    allocate(cind(np2*neg)); cind = 0.d0
    allocate(valud(nspin)); valud = 0.d0

    call check_totch(nfout)           ! totch will be checked

    nxx = nxyz_tetra(1)
    nyy = nxyz_tetra(2)
    nzz = nxyz_tetra(3)
    neig=neg
    nengy=0
                                                  __TIMER_DO_START(806)
    do ispin = ista_spin, iend_spin
       do ip2=1,np2
          ik=nspin*(ip2-1)+ispin
          if(map_k(ik) /= myrank_k) cycle               ! MPI
          do ieig=1,neig
! === DEBUG by tkato 2011/10/05 ================================================
!            ip = neordr(ieig,ik)
!            if(map_e(ip) == myrank_e) then             ! MPI
!               eig2(ip2,ieig,ispin)=eko_l(map_z(ip),ik)! MPI
             if(map_e(ieig) == myrank_e) then             ! MPI
                eig2(ip2,ieig,ispin)=eko_l(map_z(ieig),ik)! MPI
! ==============================================================================
             end if
          enddo
       enddo
    enddo
                                                  __TIMER_DO_STOP(806)

                                                  __TIMER_COMM_START_w_BARRIER(mpi_kg_world,807)
    call mpi_allreduce(eig2,eig2_mpi,np2*neg*nspin,mpi_double_precision,mpi_sum,mpi_kg_world,ierr) ! MPI
! ==== DEBUG by tkato 2012/12/19 ===============================================
    eig2 = eig2_mpi
    call mpi_allreduce(eig2,eig2_mpi,np2*neg*nspin,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)
! ==============================================================================
                                                  __TIMER_COMM_STOP(807)
    eig2 = eig2_mpi                                   ! MPI

    if(iprioccup>=2 ) write(nfout,*) ' === tetrahedron method', &
                &  ' for k-space integration ==='

    if(imag == FERRO .and. sw_fix_total_spin == YES .and. nspin == 2) then
       if(total_spin > totch) then
          write(nfout,'( " !WARNING total_spin should be smaller than totch <<m_ESoc_fermi_tetrahedron>>")')

          totch_spin(1) = totch
          totch_spin(2) = totch
       else
          totch_spin(1) = (totch + total_spin)*1.0d0
          totch_spin(2) = (totch - total_spin)*1.0d0
       end if

       do ispin = 1, nspin
          efermi2 = (minval(eig2(:,:,ispin))+maxval(eig2(:,:,ispin)))*0.5d0

          call fermi1(nxx,nyy,nzz,np2,np2,neig,neg,1, &
               &  eig2(1:np2,1:neg,ispin),ip20,np0,totch_spin(ispin), efermi2,eval,valud, &
               &  iwt,ip2cub,iprioccup,noncol )    ! -> efermi_spin
          efermi_spin(ispin) = efermi2
       end do
! --> T. Yamasaki 06th Aug. 2009
       efermi2 = efermi
       call fermi1(nxx,nyy,nzz,np2,np2,neig,neg,nspin, &
            &  eig2,ip20,np0,totch, efermi2,eval,valud, &
            &  iwt,ip2cub,iprioccup,noncol)
       if(iprioccup >=2 ) write(nfout,*) 'eval=',eval
       efermi = efermi2
       do ispin=1, nspin
          call nstt3i(idimtetra,nengy,efermi,nxx,nyy,nzz, &
               &  np2,np2,neig,eig2(1,1,ispin), &
               &  ip20,np0,eawk,cdwk,cswk,np2,neig,cdos,cind,width_tetra )
          totind=0.d0
                                                  __TIMER_DO_START(808)
          do ip2=1,np2
             do ieig=1,neig
                ikee=ip2+np2*(ieig-1)
                ik=nspin*(ip2-1)+ispin
                if(map_k(ik) == myrank_k) then            ! MPI
                   ip_mpi = neordr(ieig,ik)
                   if(map_e(ip_mpi) == myrank_e) &        ! MPI
                        & occup2(ip_mpi,ip2,ispin)=cind(ikee)*dble(np2) ! MPI
                end if
                totind=totind+cind(ikee)
             enddo
          enddo
                                                  __TIMER_DO_STOP(808)
          totch_spin0(ispin)= totind
       end do
       total_spin0 = totch_spin0(1) - totch_spin0(2)
       totch_spin0 = totch_spin0*2.0  ! 21 Jun. 2023 by T. Y.
       if(iprioccup >= PRINTOUTLEVEL_TOTCH_SPIN) &
            & call wd_efermi_and_total_spin0_Plus(nfout,total_spin0,totch_spin0,totch_spin)
       efermi = (efermi_spin(1)+efermi_spin(2))*0.5d0
       if(iprioccup>=2) &
            & write(nfout,'(" ! efermi_spin(1:2) = ",2f16.8," efermi = ",f16.8)') efermi_spin(1:2),efermi
    else
! modified by H.Sawada on May 2, 1997
       efermi2 = efermi
! modified by H.Sawada on May 2, 1997
       call fermi1(nxx,nyy,nzz,np2,np2,neig,neg,nspin, &
            &  eig2,ip20,np0,totch, efermi2,eval,valud, &
            &  iwt,ip2cub,iprioccup,noncol)
       if(iprioccup >=2 ) write(nfout,*) 'eval=',eval
       efermi = efermi2
       efermi_spin(1) = efermi2
       if(nspin == 2) efermi_spin(2) = efermi2
    end if

    do ispin=1,nspin
      call nstt3i(idimtetra,nengy,efermi_spin(ispin),nxx,nyy,nzz, &
               &  np2,np2,neig,eig2(1,1,ispin), &
               &  ip20,np0,eawk,cdwk,cswk,np2,neig,cdos,cind,width_tetra )

      totind=0.d0
                                                  __TIMER_DO_START(809)
      do ip2=1,np2
         do ieig=1,neig
            ikee=ip2+np2*(ieig-1)
            ik=nspin*(ip2-1)+ispin
            if(map_k(ik) == myrank_k) then
               ip_mpi = neordr(ieig,ik)
               if(map_e(ip_mpi) == myrank_e) &
                    & occup2(ip_mpi,ip2,ispin)=cind(ikee)*dble(np2)
            end if
            totind=totind+cind(ikee)
         enddo
      enddo
                                                  __TIMER_DO_STOP(809)
      if(imag == FERRO .and. sw_fix_total_spin == NO .and. nspin == 2)  totch_spin0(ispin) = totind
      if(iprioccup>=2) write(nfout,'(" for spin=",i1," ** TOTAL CHARGE after fermi1 = ",f12.4)') ispin,totind
    enddo
! === DEBUG by tkato 2011/10/05 ================================================
    call mpi_allreduce(occup2,occup2_mpi,neg*np2*nspin,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
    occup2 = occup2_mpi
! ==============================================================================
    if(imag == FERRO .and. sw_fix_total_spin == NO .and. nspin == 2) then
       total_spin0 = totch_spin0(1) - totch_spin0(2)
       totch_spin0 = totch_spin0*2.0  ! 21 Jun. 2023 by T. Y.
       if(iprioccup >= PRINTOUTLEVEL_TOTCH_SPIN) call wd_efermi_and_total_spin0(nfout,total_spin0,totch_spin0)
    end if

                                                  __TIMER_DO_START(810)
    do ieig = 1, neg
       if(map_e(ieig) /= myrank_e) cycle
       do ispin=1,nspin
          do ip2=1,np2
             ik=nspin*(ip2-1)+ispin
             if(map_k(ik) == myrank_k) then
! === DEBUG by tkato 2011/10/05 ================================================
!               occup_l(map_z(ieig),ik)=occup2(ieig,ip2,ispin) ! MPI
                occup_l(map_z(ieig),ik)=occup2(neordr(ieig,ik),ip2,ispin) ! MPI
! ==============================================================================
             end if
          enddo
       enddo
    enddo
                                                  __TIMER_DO_STOP(810)
    deallocate(valud); deallocate(cind); deallocate(cdos); deallocate(cswk)
    deallocate(cdwk);  deallocate(eawk); deallocate(occup2_mpi); deallocate(occup2)
    deallocate(eig2_mpi); deallocate(eig2)
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(709)
#endif
  end subroutine m_ESoc_fermi_tetrahedron_3D

!===============================================================================
  subroutine m_ESoc_fermi_ColdSmearing_3D(nfout)
!
!  An if-block of "if(imag == FERRO .and. sw_fix_total_spin == YES .and. nspin == 2) then",
! which defines efermi_spin and occup_l with fixing total spin, is coded.
!                        by T. Yamasaki   April/10/2007
!
    use m_Parallelization,      only : mpi_kg_world, mpi_ke_world

    integer, intent(in) :: nfout

    integer             :: jcount, ik, ie, is
    real(kind=DP)       :: emin, emax, e1, e2, tot, totch_spin(2), totch_spin0(2)
    real(kind=DP), parameter :: DELTA_TOTCH = 1.d-10
    integer,       parameter :: MAXITR = 100000
    real(kind=DP),allocatable, dimension(:,:) :: temp_mpi, eko_mpi, occup_mpi  ! MPI

    integer :: id_sname = -1, id_sname2 = -1
                                                  __TIMER_SUB_START(715)
    call tstatc0_begin('barrier(m_Eoc_fermi_ColdSmearing) ',id_sname2)
    if(npes > 1) call mpi_barrier(MPI_CommGroup,ierr)
    call tstatc0_end(id_sname2)

    call tstatc0_begin('m_ESoc_fermi_ColdSmearing_3D ', id_sname)

    call check_totch(nfout)           ! -(m_ES_occup) totch will be checked

    allocate(eko_mpi  (neg,kv3)); eko_mpi = 0.d0
    allocate(occup_mpi(neg,kv3)); occup_mpi = 0.d0

    !  MPI Gather eko_l
    if(nrank_e >= 2) then
                                                  __TIMER_DO_START(821)
       do is = ista_spin, iend_spin
       do ik = 1, kv3
          if(map_k(ik) /= myrank_k) cycle
          do ie = 1, neg
             if(map_e(ie) /= myrank_e) cycle
             eko_mpi(ie,ik) = eko_l(map_z(ie),ik)
          end do
       end do
       end do
                                                  __TIMER_DO_STOP(821)
                                                  __TIMER_COMM_START_w_BARRIER(mpi_kg_world,821)
       call mpi_allreduce(MPI_IN_PLACE,eko_mpi,neg*kv3,mpi_double_precision &
            &            ,mpi_sum,mpi_kg_world,ierr)
                                                  __TIMER_COMM_STOP(821)
    else
       eko_mpi = eko_l
    end if
! ==== DEBUG by tkato 2012/12/19 ===============================================
    if(nrank_k >= 2) then
       allocate(temp_mpi(neg,kv3)); temp_mpi = 0.0d0
       temp_mpi = eko_mpi
       call mpi_allreduce(temp_mpi,eko_mpi,neg*kv3,mpi_double_precision &
                         ,mpi_sum,mpi_ge_world,ierr)
       deallocate(temp_mpi)
    end if
! ==============================================================================

    emin = minval(eko_mpi)
    emax = maxval(eko_mpi)

    efermi = emax
    e1 = emin - dabs(DELTA_FermiSearchRange)
    e2 = emax + dabs(DELTA_FermiSearchRange)

    jcount = 1
    occupied_ch_equals_totch : do

       call get_occup_l_and_tot(1,1,efermi)       ! -(contained here) ->(occup_l,tot)
!          ~~~~~~~~~~~~~~~~~~~~
       if(jcount == 1 .and. tot < totch) &
            call wd_fermi_error1(nfout,emin,emax,tot,totch) ! -(b_Fermi)
       if(dabs(tot - totch) < DELTA_TOTCH) exit occupied_ch_equals_totch
       if( tot < totch) then
          e1 = efermi; efermi = efermi + (e2-efermi)/2
       else
          e2 = efermi; efermi = efermi + (e1-efermi)/2
       end if
       jcount = jcount + 1
       if(jcount > MAXITR) then
          call wd_fermi_error2 &   !-(b_Fermi)
               & (nfout,e1,e2,efermi,emin,emax,tot,totch,neg,MAXITR)
       end if
    end do occupied_ch_equals_totch
    if(iprioccup >= 2) &
         &write(nfout,'( "== Cold Smearing Method ==")')
    if(iprioccup >= 2) &
         &write(nfout,'( " efermi = ", f10.4, " tot = ", f10.4)') efermi, tot

    call get_entropic_term(band_entropy)
    if(iprioccup >= 2) &
         &write(nfout,'( " band_entropy = ", f20.14)') band_entropy


    if(imag == FERRO .and. sw_fix_total_spin == YES .and. nspin == 2) then
       if(total_spin > totch) then
          write(nfout,'( " !WARNING total_spin should be smaller than totch <<m_ESoc_fermi_parabolic>>")')

          totch_spin(1) = totch
          totch_spin(2) = totch
       else
          totch_spin(1) = (totch + total_spin)*1.0d0
          totch_spin(2) = (totch - total_spin)*1.0d0
       end if
                                                  __TIMER_DO_START(823)
       do is = 1, nspin
          emin = minval(eko_mpi)
          emax = maxval(eko_mpi)

          efermi_spin(is) = emax
          e1 = emin - dabs(DELTA_FermiSearchRange)
          e2 = emax + dabs(DELTA_FermiSearchRange)

          jcount = 1
          occupied_ch_equals_totch2 : do

             call get_occup_l_and_tot(nspin,is,efermi_spin(is))      ! -(contained here) ->(occup_l,tot)
!                ~~~~~~~~~~~~~~~~~~~~
             if(jcount == 1 .and. tot < totch_spin(is)) &
                  call wd_fermi_error1(nfout,emin,emax,tot,totch) ! -(b_Fermi)
             if(dabs(tot - totch_spin(is)) < DELTA_TOTCH) exit occupied_ch_equals_totch2
             if( tot < totch_spin(is)) then
                e1 = efermi_spin(is); efermi_spin(is) = (efermi_spin(is) + e2)*0.5d0
             else
                e2 = efermi_spin(is); efermi_spin(is) = (efermi_spin(is) + e1)*0.5d0
             end if
             jcount = jcount + 1
             if(jcount > MAXITR) then
                call wd_fermi_error2 &   !-(b_Fermi)
                     & (nfout,e1,e2,efermi_spin(is),emin,emax,tot,totch_spin(is),neg,MAXITR)
             end if
          end do occupied_ch_equals_totch2
          if(iprioccup >= 2) &
               &write(nfout,'( "== Cold smearing Method ==")')
          if(iprioccup >= 2) &
               &write(nfout,'( " efermi_spin(",i2," ) = ", f10.4, ", tot/2 = ", f10.4)') &
               & is,efermi_spin(is), tot/2
       end do
                                                  __TIMER_DO_STOP(823)
       call get_total_spin0
       if(iprioccup >= PRINTOUTLEVEL_TOTCH_SPIN) &
            & call wd_efermi_and_total_spin0_Plus(nfout,total_spin0,totch_spin0,totch_spin)
       efermi = (efermi_spin(1)+efermi_spin(2))*0.5d0
    else if(imag == FERRO .and. sw_fix_total_spin == NO .and. nspin == 2) then
       call get_total_spin0
       if(iprioccup >= PRINTOUTLEVEL_TOTCH_SPIN) call wd_efermi_and_total_spin0(nfout,total_spin0,totch_spin0)
    end if

                                                  __TIMER_DO_START(824)
    do ik = ista_k, iend_k
       do ie = 1, neg
          if(map_e(ie) == myrank_e) then
             occup_l(map_z(ie),ik) = occup_mpi(ie,ik)
          end if
       end do
    end do
                                                  __TIMER_DO_STOP(824)

    deallocate(eko_mpi  )
    deallocate(occup_mpi)

    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(715)
  contains
    subroutine get_total_spin0(sw_occupation)
      integer, optional, intent(in) :: sw_occupation
      integer :: sw_occupation_t
      sw_occupation_t = OFF
      do is = 1, nspin
         call get_occup_l_and_tot(nspin,is,efermi,sw_occupation_t)
         totch_spin0(is) = tot
      end do
      total_spin0 = (totch_spin0(1) - totch_spin0(2))*0.5d0  ! 21th Jun. 2023 by T. Y.
    end subroutine get_total_spin0

    subroutine get_occup_l_and_tot(nspin,is,efermi,sw_occupation_t)
      integer, intent(in) :: nspin, is
      real(kind=DP), intent(in) :: efermi
      integer,optional,intent(in) :: sw_occupation_t
      integer       :: k, i
      real(kind=DP) :: wspin = 1.d0, e, dos, weight
      tot = 0.d0

      do k = is, kv3, max(nspin,af+1)   !$      do k = 1, kv3, af+1
         do i = 1, neg
            e = eko_mpi(i,k)
            call coldsmearing(e,efermi,width,dos,weight)
            tot = tot + 2*weight*wspin*kv3*qwgt(k)
            if(.not.(present(sw_occupation_t).and.sw_occupation_t == OFF)) occup_mpi(i,k) = weight*wspin*kv3*qwgt(k)
         end do
      end do
      if(af == 1) then
         tot = tot/kv3 * (af+1)
      else
         tot = tot/(kv3/nspin)
      end if
      if(ipri >= 2) write(nfout,'(" efermi, tot = ",d20.12,d20.12)') efermi,tot
    end subroutine get_occup_l_and_tot

    subroutine coldsmearing(e,efermi,width,dos,weight)
      real(kind=DP), intent(in) :: e,efermi,width
      real(kind=DP), intent(out) :: dos,weight

      real(kind=DP) :: xi,xx,xf,ex
      real(kind=DP), parameter :: msq2 = -1.4142135623730950488d0 ! -sqrt(2)
      real(kind=DP), parameter :: msq2i = -0.7071067811865475244d0  ! -1/sqrt(2)
      real(kind=DP), parameter :: sqpii = 0.56418958354775628695d0 ! 1/sqrt(pi)
      real(kind=DP), parameter :: sq2pii = 0.39894228040143267794d0 ! 1/sqrt(2pi)

      real(kind=DP) :: derfc

      xi = (efermi-e)/width
      xx = xi+msq2i
      ex = exp(-xx*xx)
      xf = (2.d0+msq2*xi)

      dos = sqpii*ex*xf/width
      weight = sq2pii*ex+0.5d0*derfc(-xx)

    end subroutine coldsmearing

    subroutine get_entropic_term(s)
      real(kind=DP), intent(out) :: s

      integer       :: k, i
      real(kind=DP) :: e, entropy
      s = 0.d0
      do k = 1, kv3, af+1
         do i = 1, neg
            e = eko_mpi(i,k)
            call get_entropy(e,efermi,width,entropy)
            s = s + 2.d0*entropy*qwgt(k)
         end do
      end do
      s = s * (af+1)
    end subroutine get_entropic_term

    subroutine get_entropy(e,efermi,width,entropy)
      real(kind=DP), intent(in) :: e,efermi,width
      real(kind=DP), intent(out) :: entropy

      real(kind=DP) :: xi,xx,xf,ex
      real(kind=DP), parameter :: msq2 = -1.4142135623730950488d0 ! -sqrt(2)
      real(kind=DP), parameter :: msq2i = -0.7071067811865475244d0  ! -1/sqrt(2)
      real(kind=DP), parameter :: sqpi2i = 0.28209479177387814347d0 ! 1/(2*sqrt(pi))

      xi = (efermi-e)/width
      xx = xi+msq2i
      ex = exp(-xx*xx)
      xf = (1.d0+msq2*xi)

      entropy = sqpi2i*ex*xf

    end subroutine get_entropy

  end subroutine m_ESoc_fermi_ColdSmearing_3D

  subroutine m_ESoc_fermi_Dirac_3D(nfout)
    use m_Parallelization,      only : mpi_kg_world

    integer, intent(in) :: nfout

    integer             :: jcount, ik, ie, is
    Real(kind=DP)       :: emin, emax, e1, e2, tot, totch_spin(2), totch_spin0(2)
    real(kind=DP), parameter :: DELTA_TOTCH = 1.d-10
    integer,       parameter :: MAXITR = 100000
    real(kind=DP),allocatable, dimension(:,:) :: temp_mpi, eko_mpi, occup_mpi  ! MPI

    integer :: id_sname = -1, id_sname2 = -1
    integer :: iksnl, iupdown

    integer             :: n1, n2
    real(kind=DP)       :: c1, c2
    logical, save       :: First = .true.

    call tstatc0_begin('barrier(m_Eoc_fermi_dirac_3D) ',id_sname2)
    if(npes > 1) call mpi_barrier(MPI_CommGroup,ierr)
    call tstatc0_end(id_sname2)

    call tstatc0_begin('m_ESoc_fermi_dirac_3D ', id_sname)

    call check_totch(nfout)           ! -(m_ES_occup) totch will be checked

    if ( noncol ) then
       allocate(eko_mpi  (neg,kv3/ndim_spinor)); eko_mpi = 0.d0
       allocate(occup_mpi(neg,kv3)); occup_mpi = 0.d0
    else
       allocate(eko_mpi  (neg,kv3)); eko_mpi = 0.d0
       allocate(occup_mpi(neg,kv3)); occup_mpi = 0.d0
    endif

    if ( noncol ) then
       do ik = 1, kv3, ndim_spinor
          Do is=2, ndim_spinor
             if ( map_k(ik+is-1) /= myrank_k ) cycle
             do ie = 1, neg
                if(map_e(ie) /= myrank_e) cycle
                !
                eko_l(map_z(ie),ik +is -1 ) = 1.0D8
             End do
          End do
       end do
    endif

    if ( noncol ) then
       if ( nrank_e >= 2) then
          allocate(temp_mpi(neg,kv3/ndim_spinor)); temp_mpi = 0.d0
          do ik = 1, kv3, ndim_spinor
             if ( map_k(ik) /= myrank_k ) cycle
             iksnl = ( ik-1 )/ndim_spinor +1

             do ie = 1, neg
                if ( map_e(ie) /= myrank_e ) cycle
                temp_mpi(ie,iksnl) = eko_l(map_z(ie),ik)
             end do
          end do
          call mpi_allreduce( temp_mpi, eko_mpi, neg*kv3/ndim_spinor, &
               &    mpi_double_precision, mpi_sum,MPI_CommGroup, ierr )
          deallocate( temp_mpi )
       else
          Do ik=1, kv3, ndim_spinor
             iksnl = ( ik-1 )/ndim_spinor +1
             eko_mpi(:,iksnl) = eko_l(:,ik)
          End do
       end if
    else
       if ( nrank_e >= 2 ) then                                 ! MPI
          do is = ista_spin, iend_spin
          do ik = is, kv3-nspin+is, nspin
             if(map_k(ik) /= myrank_k) cycle                 ! MPI
             do ie = 1, neg                                  ! MPI
                if(map_e(ie) /= myrank_e) cycle              ! MPI
                eko_mpi(ie,ik) = eko_l(map_z(ie),ik)        ! MPI
             end do                                          ! MPI
          end do                                             ! MPI
          end do                                             ! MPI
          call mpi_allreduce(MPI_IN_PLACE,eko_mpi,neg*kv3,mpi_double_precision &
               &            ,mpi_sum,mpi_kg_world,ierr)
       else                                               ! MPI
          eko_mpi = eko_l                                 ! MPI
       end if                                             ! MPI

       if(nrank_k >= 2) then
          allocate(temp_mpi(neg,kv3)); temp_mpi = 0.0d0
          temp_mpi = eko_mpi
          call mpi_allreduce(temp_mpi,eko_mpi,neg*kv3,mpi_double_precision, &
               &             mpi_sum,mpi_ge_world,ierr)
          deallocate(temp_mpi)
       end if
    endif

    emin = minval(eko_mpi);     emax = maxval(eko_mpi)
#ifdef USE_MIDDLE_EIGENVAL
    if ( First ) then
       efermi = emax;           First = .false.
    else
       if ( noncol ) then
          n1 = nint(totch);  n2 = n1+1
       else
          if ( mod(nint(totch),2) == 0 ) then
             n1 = nint(totch) /2;  n2 = n1 +1
          else
             n1 = ( nint(totch) +1 )/2;  n2 = n1
          endif
       end if
       c1 = 0.0d0; c2 = 0.0d0
       Do ik=1, kv3/ndim_spinor
          c1 = c1 +eko_mpi(n1,ik);   c2 = c2 +eko_mpi(n2,ik)
       End Do
       c1 = c1 /dble(kv3/ndim_spinor);   c2 = c2 /dble(kv3/ndim_spinor)
       efermi = ( c1+c2 )/2.0d0
    endif
#else
    efermi = emax
#endif

    e1 = emin - dabs(DELTA_FermiSearchRange)
    e2 = emax + dabs(DELTA_FermiSearchRange)

    jcount = 1
! ----------------------------
    occupied_ch_equals_totch : do

       if ( noncol ) then
          call get_occup_l_and_tot_noncl_dirac( efermi, .true. )
       else
          call get_occup_l_and_tot_dirac( 1, 1, efermi, .true. )
       endif

#ifndef USE_MIDDLE_EIGENVAL
       if(jcount == 1 .and. tot < totch) &
            call wd_fermi_error1(nfout,emin,emax,tot,totch) ! -(b_Fermi)
#endif

       if(dabs(tot - totch) < DELTA_TOTCH) exit occupied_ch_equals_totch
       if( tot < totch) then
          e1 = efermi; efermi = efermi + (e2-efermi)/2
       else
          e2 = efermi; efermi = efermi + (e1-efermi)/2
       end if
       jcount = jcount + 1
       if(jcount > MAXITR) then
          call wd_fermi_error2 &   !-(b_Fermi)
               & (nfout,e1,e2,efermi,emin,emax,tot,totch,neg,MAXITR)
       end if
    end do occupied_ch_equals_totch

    if(iprioccup >= 2) &
         &write(nfout,'( "== Fermi Dirac Broadening Method ==")')
    if(iprioccup >= 2) &
         &write(nfout,'( " efermi = ", f10.4, " tot = ", f10.4)') efermi, tot

    if ( noncol ) then
       call get_entropic_term_noncl(band_entropy)
    else
       call get_entropic_term(band_entropy)
    endif

    if(iprioccup >= 2) &
         &write(nfout,'( " band_entropy = ", f20.14)') band_entropy

    if(imag == FERRO .and. sw_fix_total_spin == YES .and. nspin == 2) then
       if(total_spin > totch) then
          write(nfout,'( " !WARNING total_spin should be smaller than totch <<m_ESoc_fermi_dirac>>")')

          totch_spin(1) = totch
          totch_spin(2) = totch
       else
          totch_spin(1) = (totch + total_spin)*1.0d0
          totch_spin(2) = (totch - total_spin)*1.0d0
       end if

       do is = 1, nspin
          emin = minval(eko_mpi)
          emax = maxval(eko_mpi)

          efermi_spin(is) = emax
          e1 = emin - dabs(DELTA_FermiSearchRange)
          e2 = emax + dabs(DELTA_FermiSearchRange)

          jcount = 1
          occupied_ch_equals_totch2 : do

             call get_occup_l_and_tot_dirac(nspin,is,efermi_spin(is),.true.)
                                      ! -(contained here) ->(occup_l,tot)
             !    ~~~~~~~~~~~~~~~~~~~
             if(jcount == 1 .and. tot < totch_spin(is)) &
                  call wd_fermi_error1(nfout,emin,emax,tot,totch) ! -(b_Fermi)
             if(dabs(tot - totch_spin(is)) < DELTA_TOTCH) exit occupied_ch_equals_totch2
             if( tot < totch_spin(is)) then
                e1 = efermi_spin(is); efermi_spin(is) = (efermi_spin(is) + e2)/2
             else
                e2 = efermi_spin(is); efermi_spin(is) = (efermi_spin(is) + e1)/2
             end if
             jcount = jcount + 1
             if(jcount > MAXITR) then
                call wd_fermi_error2 &   !-(b_Fermi)
                     & (nfout,e1,e2,efermi_spin(is),emin,emax,tot,totch_spin(is),neg,MAXITR)
             end if
          end do occupied_ch_equals_totch2

          if(iprioccup >= 2) &
               &write(nfout,'( "== Fermi Dirac Broadening Method ==")')
          if(iprioccup >= 2) &
               &write(nfout,'( " efermi_spin(",i2," ) = ", f10.4, ", tot/2 = ", f10.4)') &
               & is,efermi_spin(is), tot/2.0
       end do

       call get_total_spin0 ! -> totch_spin0

       efermi = (efermi_spin(1)+efermi_spin(2))*0.5d0
       if(iprioccup>=1)then
          write(nfout,'(a,f10.4)') ' == new efermi = ',efermi
       endif

    else if(imag == FERRO .and. sw_fix_total_spin == NO .and. nspin == 2) then
       call get_total_spin0 ! ->totch_spin0
    end if

    do ik = ista_k, iend_k                              ! MPI
       do ie = 1, neg                                   ! MPI
          if(map_e(ie) == myrank_e) then                ! MPI
             occup_l(map_z(ie),ik) = occup_mpi(ie,ik)   ! MPI
             if(sw_manual_occupation==ON)then
               if(band_index(ie)>0)then
                  iupdown = 1
                  if(nspin>1.and.mod(ik,2)==0) iupdown = 2
                  occup_l(map_z(ie),ik) = occ_ext(band_index(ie),iupdown)
               endif
             endif
          end if                                        ! MPI
       end do                                           ! MPI
    end do                                              ! MPI

    if(imag == FERRO .and. nspin == 2) then
       ie = 2
       if(sw_fix_total_spin == YES) call check_occupation(nfout,ie,totch_spin)
    else if(imag /= FERRO) then
       ie = 1
       totch_spin0(1) = totch
    end if

    if ( noncol ) then
       call check_occupation_noncl(nfout,totch_spin0)
    else
       call check_occupation(nfout,ie,totch_spin0)
    endif

    if ( af == 1 ) then
       Do ik=ista_k, iend_k, nspin
          occup_l(:,ik+1) = occup_l(:,ik)
       End Do
    endif

    deallocate(eko_mpi  )
    deallocate(occup_mpi)

    call tstatc0_end(id_sname)

  contains

    subroutine get_total_spin0
      do is = 1, nspin
         call get_occup_l_and_tot_dirac(nspin,is,efermi,.false.)
         totch_spin0(is) = tot
      end do

      total_spin0 = (totch_spin0(1) - totch_spin0(2))*0.5d0
      if(iprioccup >= PRINTOUTLEVEL_TOTCH_SPIN) then
         if(sw_fix_total_spin == YES) then
            call wd_efermi_and_total_spin0_Plus(nfout,total_spin0,totch_spin0,totch_spin)
         else
            call wd_efermi_and_total_spin0(nfout,total_spin0,totch_spin0)
         end if
      end if
    end subroutine get_total_spin0

    subroutine get_occup_l_and_tot_dirac(nspin,is,efermi,update_occ)
      integer, intent(in) :: nspin, is
      real(kind=DP), intent(in) :: efermi
      logical,intent(in) :: update_occ

      integer       :: k, i,iupdown
      real(kind=DP) :: wspin = 1.d0, e, dos, weight, tmp

      tot = 0.d0

      do k = is, kv3, max(nspin,af+1)
         do i = 1, neg
            e = eko_mpi(i,k)

            call width_fermi_dirac(e,efermi,width,dos,weight)  ! -(b_Fermi)

            if(sw_manual_occupation==ON.and..not.update_occ)then
               if(band_index(i)>0)then
                  iupdown = 1
                  if(nspin>1.and.mod(k,2)==0) iupdown = 2
                  weight = occ_ext(band_index(i),iupdown)
               endif
            endif

            tmp = weight*wspin*kv3*qwgt(k)
            if(update_occ) occup_mpi(i,k) = tmp

            tot = tot + 2*tmp
         end do
      end do
      if(af == 1) then
         tot = tot/kv3 * (af+1)
      else
         tot = tot/(kv3/nspin)
      end if
      if(ipri >= 2) write(nfout,'(" efermi, tot = ",d20.12,d20.12)') efermi,tot
    end subroutine get_occup_l_and_tot_dirac

    subroutine get_occup_l_and_tot_noncl_dirac( efermi, update_occ )
      real(kind=DP), intent(in) :: efermi
      logical,intent(in) :: update_occ

      integer       :: k, i, iksnl
      real(kind=DP) :: wspin = 1.d0, e, dos, weight, wktmp, ctmp

      tot = 0.d0;  occup_mpi = 0.0d0
      wktmp = kv3 / ndim_spinor

      do k = 1, kv3, ndim_spinor
         iksnl = ( k -1 )/ndim_spinor + 1

         do i = 1, neg
            e = eko_mpi(i,iksnl)

            call width_fermi_dirac(e,efermi,width,dos,weight)  ! -(b_Fermi)

            ctmp = weight *wspin *wktmp *qwgt(k)
            if ( update_occ ) occup_mpi(i,k) = ctmp
            tot = tot + ctmp
         end do
      end do
      tot = tot / wktmp

      if(ipri >= 2) write(nfout,'(" efermi, tot = ",d20.12,d20.12)') efermi,tot
    end subroutine get_occup_l_and_tot_noncl_dirac

    subroutine get_entropic_term(s)
      real(kind=DP), intent(out) :: s

      integer       :: k, i
      real(kind=DP) :: e, entropy
      s = 0.d0
      do k = 1, kv3, af+1
         do i = 1, neg
            e = eko_mpi(i,k)
            call get_entropy(e,efermi,width,entropy)
            s = s + 2.d0*entropy*qwgt(k)
         end do
      end do
      s = s * (af+1)
    end subroutine get_entropic_term

    subroutine get_entropic_term_noncl(s)
      real(kind=DP), intent(out) :: s

      integer       :: k, i, iksnl
      real(kind=DP) :: e, entropy
      s = 0.d0
      do k = 1, kv3, ndim_spinor
         iksnl = ( k -1 )/ndim_spinor + 1
         do i = 1, neg
            e = eko_mpi(i,iksnl)
            call get_entropy(e,efermi,width,entropy)
            s = s + 1.d0*entropy*qwgt(k)
         end do
      end do
    end subroutine get_entropic_term_noncl

    subroutine get_entropy( e, efermi, width, entropy )
      real(kind=DP), intent(in) :: e,efermi,width
      real(kind=DP), intent(out) :: entropy

      real(kind=DP) :: ee, c1, occ, unocc

      ee = ( e -efermi )/width
      c1 = exp( ee )

      occ = 1.0D0 / ( 1.0D0 +c1 )
      unocc = 1.0D0 - occ

      entropy = 0.0d0

      if ( occ > 1.0D-12 .and. unocc > 1.0D-12 ) then
         entropy = -occ*log( occ ) -unocc*log( unocc )
      endif

    end subroutine get_entropy

  end subroutine m_ESoc_fermi_Dirac_3D
!===============================================================================
!$$#endif

! ==================== KT_add ============================ 13.0E

  subroutine m_ESoc_methfessel_paxton(nfout)
    use m_Parallelization, only : mpi_kg_world
    integer, intent(in) :: nfout

    integer             :: jcount, ik, ie, is
    Real(kind=DP)       :: emin, emax, e1, e2, tot, totch_spin(2), totch_spin0(2)
    real(kind=DP)       :: emin0,emax0
    real(kind=DP), parameter :: DELTA_TOTCH = 1.d-10
    integer,       parameter :: MAXITR = 100000
    real(kind=DP),allocatable, dimension(:,:) :: temp_mpi, eko_mpi, occup_mpi  ! MPI

    integer :: id_sname = -1, id_sname2 = -1

    integer :: iksnl
    integer :: iupdown
    integer :: icount
    real(kind=DP) :: delta_energy

    call tstatc0_begin('barrier(m_Eoc_methfessel_paxton) ',id_sname2)
    if(npes > 1) call mpi_barrier(MPI_CommGroup,ierr)
    call tstatc0_end(id_sname2)

    call tstatc0_begin('m_ESoc_methfessel_paxton ', id_sname)

    call check_totch(nfout)           ! -(m_ES_occup) totch will be checked


    if ( noncol ) then
      allocate(eko_mpi  (neg,kv3/ndim_spinor)); eko_mpi = 0.d0
      allocate(occup_mpi(neg,kv3)); occup_mpi = 0.d0
    else
      allocate(eko_mpi  (neg,kv3)); eko_mpi = 0.d0
      allocate(occup_mpi(neg,kv3)); occup_mpi = 0.d0
    endif

    if ( noncol ) then
       do ik = 1, kv3, ndim_spinor
          Do is=2, ndim_spinor
            if ( map_k(ik+is-1) /= myrank_k ) cycle
             do ie = 1, neg
                if(map_e(ie) /= myrank_e) cycle
!
                eko_l(map_z(ie),ik +is -1 ) = 1.0D8
             End do
          End do
       end do
    endif

    if ( noncol ) then
      if ( npes >= 2) then
         allocate(temp_mpi(neg,kv3/ndim_spinor)); temp_mpi = 0.d0
         do ik = 1, kv3, ndim_spinor
            if ( map_k(ik) /= myrank_k ) cycle
            iksnl = ( ik-1 )/ndim_spinor +1

            do ie = 1, neg
              if ( map_e(ie) /= myrank_e ) cycle
              temp_mpi(ie,iksnl) = eko_l(map_z(ie),ik)
            end do
         end do
         call mpi_allreduce( temp_mpi, eko_mpi, neg*kv3/ndim_spinor, &
                       	&    mpi_double_precision, mpi_sum,MPI_CommGroup, ierr )
         deallocate( temp_mpi )
      else
         Do ik=1, kv3, ndim_spinor
            iksnl = ( ik-1 )/ndim_spinor +1
            eko_mpi(:,iksnl) = eko_l(:,ik)
         End do
      end if
    else
      do is = ista_spin, iend_spin
      do ik = is, kv3-nspin+is, nspin
         if(map_k(ik) /= myrank_k) cycle                 ! MPI
         do ie = 1, neg                                  ! MPI
            if(map_e(ie) /= myrank_e) cycle              ! MPI
            eko_mpi(ie,ik) = eko_l(map_z(ie),ik)        ! MPI
         end do                                          ! MPI
      end do                                             ! MPI
      end do                                             ! MPI
      call mpi_allreduce(MPI_IN_PLACE,eko_mpi,neg*kv3,mpi_double_precision &
           &            ,mpi_sum,mpi_kg_world,ierr)
      allocate(temp_mpi(neg,kv3)); temp_mpi = 0.0d0
      temp_mpi = eko_mpi
      call mpi_allreduce(temp_mpi,eko_mpi,neg*kv3,mpi_double_precision, &
           &             mpi_sum,mpi_ge_world,ierr)
      deallocate(temp_mpi)
    end if
    if(esearch==ON .and. order_mp>0)then
      delta_energy = esearch_factor_mp*width/dble(order_mp)
      emin = minval(eko_mpi);     emax = maxval(eko_mpi)
      !emin = emin+(emax-emin)*0.5d0*totch*0.5d0/dble(neg)
      emax = emin+delta_energy
      if(iprioccup>=2) write(nfout,'(a,2f10.5)') 'initial emin, emax ',emin,emax
      efermi = emax
      icount = 0
      do
         if ( noncol ) then
            call get_occup_l_and_tot_noncl_mp( efermi, .true. )
         else
            call get_occup_l_and_tot_mp( 1, 1, efermi, .true. )
         endif
         icount = icount+1
         if (tot>=totch) then
            emax = efermi
            exit
         endif
         emin = efermi
         efermi = efermi+delta_energy
      enddo
      emin0 = emin
      emax0 = emax
      efermi = emin
      e1 = emin
      e2 = emax
    else
      emin = minval(eko_mpi);     emax = maxval(eko_mpi) ; efermi = emax
      e1 = emin - dabs(DELTA_FermiSearchRange)
      e2 = emax + dabs(DELTA_FermiSearchRange)
    endif

    jcount = 1
! ----------------------------
    occupied_ch_equals_totch : do

       if ( noncol ) then
          call get_occup_l_and_tot_noncl_mp( efermi, .true. )
       else
          call get_occup_l_and_tot_mp( 1, 1, efermi, .true. )
       endif

       !if(jcount == 1 .and. tot < totch) &
       !     call wd_fermi_error1(nfout,emin,emax,tot,totch) ! -(b_Fermi)

       if(dabs(tot - totch) < DELTA_TOTCH) exit occupied_ch_equals_totch
       if( tot < totch) then
          e1 = efermi; efermi = efermi + (e2-efermi)/2
       else
          e2 = efermi; efermi = efermi + (e1-efermi)/2
       end if
       jcount = jcount + 1
       if(jcount > MAXITR) then
          call wd_fermi_error2 &   !-(b_Fermi)
               & (nfout,e1,e2,efermi,emin,emax,tot,totch,neg,MAXITR)
       end if
    end do occupied_ch_equals_totch
    if(iprioccup >= 2 .and. esearch==ON) write(nfout,'(a,i8,3f15.5)') &
       & ' !** number of iteration for bracketing, emin0, emax0 and efermi: ',icount,emin0,emax0,efermi

    if(iprioccup >= 2) &
         &write(nfout,'( "== Methfessel-Paxton Broadening Method ==")')
    if(iprioccup >= 2) &
         &write(nfout,'( " efermi = ", f10.4, " tot = ", f10.4)') efermi, tot

    if ( noncol ) then
       call get_entropic_term_noncl(band_entropy)
    else
       call get_entropic_term(band_entropy)
    endif

    if(iprioccup >= 2) &
         &write(nfout,'( " band_entropy = ", f20.14)') band_entropy

    if(imag == FERRO .and. sw_fix_total_spin == YES .and. nspin == 2) then
       if(total_spin > totch) then
          write(nfout,'( " !WARNING total_spin should be smaller than totch <<m_ESoc_methfessel_paxton>>")')

          totch_spin(1) = totch
          totch_spin(2) = totch
       else
          totch_spin(1) = (totch + total_spin)*1.0d0
          totch_spin(2) = (totch - total_spin)*1.0d0
       end if

       do is = 1, nspin
          emin = minval(eko_mpi)
          emax = maxval(eko_mpi)

          efermi_spin(is) = emax
          e1 = emin - dabs(DELTA_FermiSearchRange)
          e2 = emax + dabs(DELTA_FermiSearchRange)

          jcount = 1
          occupied_ch_equals_totch2 : do

             call get_occup_l_and_tot_mp(nspin,is,efermi_spin(is),.true.)
                                      ! -(contained here) ->(occup_l,tot)
             !    ~~~~~~~~~~~~~~~~~~~
             if(jcount == 1 .and. tot < totch_spin(is)) &
                  call wd_fermi_error1(nfout,emin,emax,tot,totch) ! -(b_Fermi)
             if(dabs(tot - totch_spin(is)) < DELTA_TOTCH) exit occupied_ch_equals_totch2
             if( tot < totch_spin(is)) then
                e1 = efermi_spin(is); efermi_spin(is) = (efermi_spin(is) + e2)/2
             else
                e2 = efermi_spin(is); efermi_spin(is) = (efermi_spin(is) + e1)/2
             end if
             jcount = jcount + 1
             if(jcount > MAXITR) then
                call wd_fermi_error2 &   !-(b_Fermi)
                     & (nfout,e1,e2,efermi_spin(is),emin,emax,tot,totch_spin(is),neg,MAXITR)
             end if
          end do occupied_ch_equals_totch2

          if(iprioccup >= 2) &
               &write(nfout,'( "== Methfessel-Paxton Broadening Method ==")')
          if(iprioccup >= 2) &
               &write(nfout,'( " efermi_spin(",i2," ) = ", f10.4, ", tot/2 = ", f10.4)') &
               & is,efermi_spin(is), tot/2.0
       end do

       call get_total_spin0 ! -> totch_spin0

       efermi = (efermi_spin(1)+efermi_spin(2))*0.5d0
       if(iprioccup>=1)then
          write(nfout,'(a,f10.4)') ' == new efermi = ',efermi
       endif

    else if(imag == FERRO .and. sw_fix_total_spin == NO .and. nspin == 2) then
       call get_total_spin0 ! ->totch_spin0
    end if

    do ik = ista_k, iend_k                              ! MPI
       do ie = 1, neg                                   ! MPI
          if(map_e(ie) == myrank_e) then                ! MPI
             occup_l(map_z(ie),ik) = occup_mpi(ie,ik)   ! MPI
             if(sw_manual_occupation==ON)then
               if(band_index(ie)>0)then
                  iupdown = 1
                  if(nspin>1.and.mod(ik,2)==0) iupdown = 2
                  occup_l(map_z(ie),ik) = occ_ext(band_index(ie),iupdown)
               endif
             endif
          end if                                        ! MPI
       end do                                           ! MPI
    end do                                              ! MPI

    if(imag == FERRO .and. nspin == 2) then
       ie = 2
       if(sw_fix_total_spin == YES) call check_occupation(nfout,ie,totch_spin)
    else if(imag /= FERRO) then
       ie = 1
       totch_spin0(1) = totch
    end if

    if ( noncol ) then
       call check_occupation_noncl(nfout,totch_spin0)
    else
       call check_occupation(nfout,ie,totch_spin0)
    endif

    if ( af == 1 ) then
       Do ik=ista_k, iend_k, nspin
          occup_l(:,ik+1) = occup_l(:,ik)
       End Do
    endif

    deallocate(eko_mpi  )
    deallocate(occup_mpi)

    call tstatc0_end(id_sname)

  contains

    subroutine get_total_spin0
      do is = 1, nspin
         call get_occup_l_and_tot_mp(nspin,is,efermi,.false.)
         totch_spin0(is) = tot
      end do

      total_spin0 = (totch_spin0(1) - totch_spin0(2))*0.5d0
      if(iprioccup >= PRINTOUTLEVEL_TOTCH_SPIN) then
         if(sw_fix_total_spin == YES) then
            call wd_efermi_and_total_spin0_Plus(nfout,total_spin0,totch_spin0,totch_spin)
         else
            call wd_efermi_and_total_spin0(nfout,total_spin0,totch_spin0)
         end if
      end if
    end subroutine get_total_spin0

    subroutine get_occup_l_and_tot_mp(nspin,is,efermi,update_occ)
      integer, intent(in) :: nspin, is
      real(kind=DP), intent(in) :: efermi
      logical,intent(in) :: update_occ

      integer       :: k, i,iupdown
      real(kind=DP) :: wspin = 1.d0, e, dos, weight, tmp, ent

      tot = 0.d0

      do k = is, kv3, max(nspin,af+1)
         do i = 1, neg
            e = eko_mpi(i,k)

            call width_methfessel_paxton(order_mp,e,efermi,width,dos,weight,ent)  ! -(b_Fermi)

            if(sw_manual_occupation==ON.and..not.update_occ)then
               if(band_index(i)>0)then
                  iupdown = 1
                  if(nspin>1.and.mod(k,2)==0) iupdown = 2
                  weight = occ_ext(band_index(i),iupdown)
               endif
            endif

            tmp = weight*wspin*kv3*qwgt(k)
            if(update_occ) occup_mpi(i,k) = tmp

            tot = tot + 2*tmp
         end do
      end do
      if(af == 1) then
         tot = tot/kv3 * (af+1)
      else
         tot = tot/(kv3/nspin)
      end if
      if(iprioccup >= 2) write(nfout,'(" efermi, tot = ",d20.12,d20.12)') efermi,tot
    end subroutine get_occup_l_and_tot_mp

    subroutine get_occup_l_and_tot_noncl_mp( efermi, update_occ )
      real(kind=DP), intent(in) :: efermi
      logical,intent(in) :: update_occ

      integer       :: k, i, iksnl
      real(kind=DP) :: wspin = 1.d0, e, dos, weight, wktmp, ctmp, ent

      tot = 0.d0;  occup_mpi = 0.0d0
      wktmp = kv3 / ndim_spinor

      do k = 1, kv3, ndim_spinor
         iksnl = ( k -1 )/ndim_spinor + 1

         do i = 1, neg
            e = eko_mpi(i,iksnl)

            call width_methfessel_paxton(order_mp,e,efermi,width,dos,weight,ent)  ! -(b_Fermi)

            ctmp = weight *wspin *wktmp *qwgt(k)
            if ( update_occ ) occup_mpi(i,k) = ctmp
            tot = tot + ctmp
         end do
      end do
      tot = tot / wktmp

      if(ipri >= 2) write(nfout,'(" efermi, tot = ",d20.12,d20.12)') efermi,tot
    end subroutine get_occup_l_and_tot_noncl_mp

    subroutine get_entropic_term(s)
      real(kind=DP), intent(out) :: s

      integer       :: k, i
      real(kind=DP) :: e, entropy
      s = 0.d0
      do k = 1, kv3, af+1
         do i = 1, neg
            e = eko_mpi(i,k)
            call get_entropy(e,efermi,width,entropy)
            s = s + 2.d0*entropy*qwgt(k)
         end do
      end do
      s = s * (af+1)
    end subroutine get_entropic_term

    subroutine get_entropic_term_noncl(s)
      real(kind=DP), intent(out) :: s

      integer       :: k, i, iksnl
      real(kind=DP) :: e, entropy
      s = 0.d0
      do k = 1, kv3, ndim_spinor
         iksnl = ( k -1 )/ndim_spinor + 1
         do i = 1, neg
            e = eko_mpi(i,iksnl)
            call get_entropy(e,efermi,width,entropy)
            s = s + 1.d0*entropy*qwgt(k)
         end do
      end do
    end subroutine get_entropic_term_noncl

    subroutine get_entropy( e, efermi, width, entropy )
      real(kind=DP), intent(in) :: e,efermi,width
      real(kind=DP), intent(out) :: entropy

      real(kind=DP) :: occ,dos
      call width_methfessel_paxton(order_mp,e,efermi,width,dos,occ,entropy)
    end subroutine get_entropy

  end subroutine m_ESoc_methfessel_paxton
! ======================================================== 13.0E

  subroutine m_ES_rd_occ_ext(nfout)
    integer, intent(in) :: nfout
    integer :: f_getIntValue,f_selectBlock, f_selectParentBlock,f_getRealValue
    integer :: iret
    if(f_selectBlock('accuracy')==0)then
       if(f_getIntValue(tag_sw_manual_occupation,iret)==0)then
          sw_manual_occupation = iret
       endif
       if(sw_manual_occupation==ON)then
          call set_external_occupation()
       endif
       iret = f_selectParentBlock()
    endif
  contains
    subroutine set_external_occupation()
      integer :: f_selectFirstTableLine, f_selectNextTableLine,f_getRealValue
      integer :: iocc,io,iret
      real(kind=DP) :: dret
      if(f_selectBlock(tag_occupation)==0)then
         iocc = 1
         do
           if(iocc==1)then
              if(f_selectFirstTableLine() /= 0 ) then
                 exit
              endif
           else
              if(f_selectNextTableLine() /= 0 )then
                 exit
              endif
           endif
           iocc = iocc+1
         enddo
         iocc = iocc-1
         allocate(occ_ext(iocc,nspin));occ_ext=0.d0
         allocate(band_index(neg));band_index=0
         iret = f_selectFirstTableLine()
         do io=1,iocc
            if(f_getIntValue(tag_band_index,iret)==0) band_index(iret) = io
            if(f_getRealValue(tag_occ,dret,'')==0) occ_ext(io,1) = dret
            if(f_getRealValue(tag_occ_up,dret,'')==0) occ_ext(io,1) = dret
            if(nspin>1)then
               if(f_getRealValue(tag_occ_down,dret,'')==0) occ_ext(io,2) = dret
            endif
            iret = f_selectNextTableLine()
         enddo
         if(printable.and.sw_manual_occupation==ON)then
            write(nfout,'(a)') ' !** manual occupation given in the input'
            do io=1,neg
               if(band_index(io)<1) cycle
               if(nspin==1) then
               write(nfout,'(a,i8,a,f10.5)') ' !** band : ',io,' occ : ',occ_ext(band_index(io),1)
               else
               write(nfout,'(a,i8,a,f10.5,a,f10.5)') ' !** band : ',io,' occ (up) : ', &
            &          occ_ext(band_index(io),1),' occ (down) : ',occ_ext(band_index(io),2)
               endif
            enddo
         endif
         iret = f_selectParentBlock()
      endif
    end subroutine set_external_occupation
  end subroutine m_ES_rd_occ_ext

  subroutine wd_efermi_and_total_spin0(nfout,total_spin0,totch_spin0)
    ! Revised by T. Y. 21th Jun. 2023
    integer, intent(in) :: nfout
    real(kind=DP), intent(in) :: total_spin0, totch_spin0(2)
    write(nfout,'(" == efermi = ",f10.4, ", totch_diff0 = ",f12.6, ", totch_spin0(1:2) = ",2f12.6)') &
         & efermi, total_spin0, totch_spin0(1:2)*0.5d0
  end subroutine wd_efermi_and_total_spin0

  subroutine wd_efermi_and_total_spin0_plus(nfout,total_spin0,totch_spin0,totch_spin)
    ! Revised by T. Y. 21th Jun. 2023
    integer, intent(in) :: nfout
    real(kind=DP), intent(in) :: total_spin0, totch_spin0(2),totch_spin(2)
    write(nfout,'(" == efermi_spin(1:2) = ",2f10.4, ", totch_diff  = ",f12.6,", totch_spin(1:2)  = ",2f12.6)') &
         & efermi_spin(1:2),total_spin, totch_spin(1:2)*0.5d0
    write(nfout,'(" == efermi = ",f10.4, 20x,", totch_diff0 = ",f12.6, ", totch_spin0(1:2) = ",2f12.6)') &
         & efermi, total_spin0, totch_spin0(1:2)*0.5d0
  end subroutine wd_efermi_and_total_spin0_plus

! =============== KT_add ========================== 13.0E
  subroutine m_ESoc_count_charge_belowEF( nfout )
    integer, intent(in) :: nfout

    integer       :: ik, ie, is
    real(kind=DP) :: wspin = 1.d0, e, dos, weight
    real(kind=DP) :: wktmp, ctmp, csum, csum_mpi
!
    csum = 0.0d0

    if ( noncol ) then
       wktmp = kv3 / ndim_spinor

       do ik = 1, kv3, ndim_spinor
          if ( map_k(ik) /= myrank_k ) cycle

          do ie = 1, neg
             if ( map_e(ie) /= myrank_e ) cycle

             e = eko_l(map_z(ie),ik)
             if ( e > efermi ) cycle

             call width_fermi_dirac(e,efermi,width,dos,weight)  ! -(b_Fermi)
             ctmp = weight *wspin *wktmp *qwgt(ik)
             csum = csum +ctmp
          end do
       end do

       if ( npes > 2 ) then
          call mpi_allreduce( csum, csum_mpi, 1, mpi_double_precision, &
               &              mpi_sum, MPI_CommGroup, ierr )
          csum = csum_mpi
       endif

       csum = csum /dble(wktmp)

    else
       do is = ista_spin, iend_spin
       do ik = is, kv3-nspin+is, nspin
          if (map_k(ik) /= myrank_k) cycle
          do ie=1, neg
             if (map_e(ie) /= myrank_e) cycle

             e = eko_l(map_z(ie),ik)
             if ( e > efermi ) cycle

             call width_fermi_dirac(e,efermi,width,dos,weight)  ! -(b_Fermi)

             ctmp = weight *wspin *kv3 *qwgt(ik)
             csum = csum + 2*ctmp
          end do
       end do
       end do

       if ( npes > 2 ) then
          call mpi_allreduce( csum, csum_mpi, 1, mpi_double_precision, &
               &              mpi_sum, MPI_CommGroup, ierr )
          csum = csum_mpi
       endif
       csum = csum/dble(kv3)/dble(nspin)

    endif
    csum = csum /dble(nrank_g)

    write(nfout,'(A,F20.12)') " !! Number of Charge below Fermi level = ", csum

  end subroutine m_ESoc_count_charge_belowEF

  subroutine m_ESoc_count_charge_belowEF_ek( nfout )
    integer, intent(in) :: nfout

    integer       :: k, i
    real(kind=DP) :: wspin = 1.d0, e, dos, weight
    real(kind=DP) :: totw, ctmp, csum, wktmp

    csum = 0.0d0

    if ( noncol ) then
       wktmp = kv3_ek /ndim_spinor

       do k=1, kv3_ek, ndim_spinor
          do i=1, neg
             e = eko_ek(i,k)
             if ( e > efermi ) cycle

             call width_fermi_dirac(e,efermi,width,dos,weight)  ! -(b_Fermi)
             totw = weight *wspin *wktmp *qwgt_ek(k)
             csum = csum + totw
          end do
       end do
       csum = csum /wktmp
    else
       do k = 1, kv3_ek, af+1
          do i = 1, neg
             e = eko_ek(i,k)
             if ( e > efermi ) cycle

             call width_fermi_dirac(e,efermi,width,dos,weight)  ! -(b_Fermi)
             totw = weight*wspin*kv3_ek*qwgt_ek(k)
             csum = csum + 2*totw
          end do
       end do
       csum = csum /kv3_ek *(af+1)
    endif

    write(nfout,'(A,F20.12,A,F20.12)') &
         &     "Total Chargeum of Charge below Fermi level = ", totch, &
         &     " Num of Charge below EF = ", csum

  end subroutine m_ESoc_count_charge_belowEF_ek
! ==================================================== 13.0E

  subroutine m_ESoc_occup_fix_3D(nfout)
    use m_Parallelization, only : mpi_kg_world

    integer, intent(in) :: nfout
    real(kind=DP),allocatable, dimension(:,:) :: eko_mpi  ! MPI
    real(kind=DP),allocatable, dimension(:) :: occup_mpi
    integer :: ik, ie, jcount, nb, iksnl
    real(kind=DP) :: emin, emax, e1, e2, tot
    real(kind=DP), parameter :: DELTA_TOTCH = 1.d-10
    integer,       parameter :: MAXITR = 10000
    integer, allocatable, dimension(:,:) :: n2_mpi  ! MPI
    real(DP),allocatable, dimension(:,:) :: e2_mpi  ! MPI
    real(kind=DP),allocatable, dimension(:,:) :: temp_mpi
    integer :: itmp1, itmp2, nocc, is
    real(DP) :: maxval_vbm(2), vbm_spin(2), c1
    real(DP) :: minval_cbm(2), cbm_spin(2)

    maxval_vbm = -1.0D99
    minval_cbm =  1.0D99

    itmp1 = nint(totch);    itmp2 = nint(total_spin)
    c1 = totch -itmp1
    if ( abs(c1) > 1E-10 ) then
       call phase_error_with_msg(nfout, 'totch should be an integer',__LINE__,__FILE__)
    endif
    if ( ndim_magmom == 1 ) then
       if ( mod(itmp1,2) == 1 ) then
          call phase_error_with_msg(nfout, 'totch should be an even number',&
               &                    __LINE__,__FILE__)
       endif
    else if ( ndim_magmom == 2 ) then
       if ( mod(itmp1+itmp2,2) == 1 ) then
          call phase_error_with_msg(nfout, 'totch+total_spin should be an even number',&
               &                    __LINE__,__FILE__)
       endif
    endif
    if ( itmp2 /= 0 ) then
       if ( sw_fix_total_spin /= ON ) then
          call phase_error_with_msg(nfout, 'sw_fix_total_spin should be ON',&
               &                    __LINE__,__FILE__)
       endif
    endif

    occup_l = 0.0d0
    do ik = ista_k, iend_k, ndim_spinor
       if ( ndim_magmom == 1 ) then
          nocc = itmp1 /2
       else if ( ndim_magmom == 2 ) then
          if ( mod(ik,2)==1 ) then
             nocc = ( itmp1 +itmp2 ) /2
          else
             nocc = ( itmp1 -itmp2 ) /2
          endif
       else if ( ndim_magmom == 4 ) then
          nocc = itmp1
       endif

       is = 1
       if ( ndim_magmom == 2 .and. mod(ik,2)==0 ) is = 2

       vbm_spin = -1.0D99
       cbm_spin =  1.0D99
       do ie = 1, neg
          if(map_e(ie) == myrank_e) then
             if ( neordr(ie,ik) <= nocc ) then
                occup_l(map_z(ie),ik) = qwgt(ik) *kv3 /ndim_spinor
             else
                occup_l(map_z(ie),ik) = 0.0d0
             endif
             if ( neordr(ie,ik) == nocc )   vbm_spin(is) = eko_l(map_z(ie),ik)
             if ( neordr(ie,ik) == nocc+1 ) cbm_spin(is) = eko_l(map_z(ie),ik)
          endif
       end do
       maxval_vbm(is) = max( vbm_spin(is), maxval_vbm(is) )
       minval_cbm(is) = min( cbm_spin(is), minval_cbm(is) )
    end do
    if ( npes > 1 ) then
       call mpi_allreduce( maxval_vbm, vbm_spin, 2, mpi_double_precision, &
            &              mpi_max, MPI_CommGroup, ierr )
       call mpi_allreduce( minval_cbm, cbm_spin, 2, mpi_double_precision, &
            &              mpi_min, MPI_CommGroup, ierr )
    endif
    if ( ndim_magmom == 2 ) then
       efermi_spin(:) = ( vbm_spin(:) +cbm_spin(:) ) /2.0d0
       efermi = ( efermi_spin(1) +efermi_spin(2) ) /2.0d0
    else
       efermi = ( vbm_spin(1) +cbm_spin(1) ) /2.0d0
    endif

    if ( vbm_spin(1) > cbm_spin(1) .or. vbm_spin(2) > cbm_spin(2) ) then
       write(nfout,*) "---- WARNING ( unphysical occupations ) ----"
       if ( vbm_spin(1) > cbm_spin(1) ) then
          write(nfout,'(A,F16.12,A,F16.12)') &
               &      " vbm_spin(1):", vbm_spin(1), " > cbm_spin(1):", cbm_spin(1)
       endif
       if ( vbm_spin(2) > cbm_spin(2) ) then
          write(nfout,'(A,F16.12,A,F16.12)') &
               &      " vbm_spin(2):", vbm_spin(2), " > cbm_spin(2):", cbm_spin(2)
       endif
       write(nfout,*) "--------------------------------------------"
    endif

  end subroutine m_ESoc_occup_fix_3D

end module m_ES_occup
