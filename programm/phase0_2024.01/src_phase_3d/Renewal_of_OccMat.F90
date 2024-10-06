!=======================================================================
!
!  SOFTWARE NAME : PHASE (ver. 7.01)
!
!  SUBROUINE: Renewal_of_OccMat
!
!  AUTHOR(S): T. Yamamoto   June/08/2005
!  
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
!=======================================================================
!   patch 0.1 by K. Tagami @adv    2009/02/04
! 
!   patch 0.1 :  change from m_OP_occ_mat_ao to m_OP_occ_mat_ao_kt
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
! =========================================== modified by K. Tagami ========== 5.0
!!subroutine Renewal_of_OccMat( pmode )
subroutine Renewal_of_OccMat( flag_dm, pmode, skip )
! ============================================================================ 5.0

! $Id: Renewal_of_OccMat.F90 577 2017-12-15 05:26:08Z jkoga $
  use m_Files,               only : nfout
  use m_Const_Parameters,   only : SPHERICAL_HARMONICS, ATOMIC_ORBITAL, ON
  use m_Control_Parameters, only : num_projectors, projector_type, sw_hubbard
  use m_Kpoints,            only : kv3
  use m_Charge_Density,     only : m_CD_den_mat
! ============================================ Modified by K. Tagami ======== 0.1
!  use m_Orbital_Population,only  : m_OP_occ_mat_ylm, m_OP_occ_mat_ao &
!       &                         , m_OP_rd_occ_mat 
  use m_Orbital_Population,only  : m_OP_occ_mat_ylm, m_OP_occ_mat_ao &
       &                         , m_OP_rd_occ_mat &
       &                         , m_OP_occ_mat_ao_kt, m_OP_cp_om_to_ommix, m_OP_cp_ommix_to_omold
! =======================================

! ================================= added by K. Tagami ===================== 5.0
  use m_Control_Parameters,      only : occmat_type
  use m_Const_Parameters,       only : OccMat_Type1, OccMat_Type2
! ============================================================================ 5.0

! ================================ added by K. Tagami =============== 11.0
  use m_Control_parameters,    only : noncol
  use m_Orbital_Population,    only  : m_OP_occ_mat_ylm_noncl, &
       &                               m_OP_occ_mat_ao_kt_noncl
! ================================================================== 11.0

  implicit none
  logical, intent(in) :: flag_dm
  integer, intent(in) :: pmode
  logical, intent(in), optional :: skip
  logical :: ski
  ski = .false.
  if(present(skip)) then
     ski = skip
  endif
  if (ski) then
     call m_OP_cp_om_to_ommix(nfout,1.d0)
     call m_OP_cp_ommix_to_omold()
     return
  endif
  
#ifdef __TIMER_SUB__
  call timer_sta(735)
#endif

! ================================= modified by K. Tagami =============== 5.0
!  if(num_projectors>0 .and. sw_hubbard == ON) then
!     if(projector_type == SPHERICAL_HARMONICS) then
!        !!$write(*,*) "m_CD_occ_mat_ylm"
!        if(flag_dm) call m_CD_den_mat(nfout,kv3)
!        call m_OP_occ_mat_ylm(nfout)
!     else if(projector_type == ATOMIC_ORBITAL) then
!!!! ====================================== Modified by K. Tagami ============ 0.1
!!        call m_OP_occ_mat_ao(nfout)
!        call m_OP_occ_mat_ao_kt(nfout)
!! ============================================================
!     end if
!  end if

! ====================== modified by K. Tagami ========== 11.0
  if ( num_projectors>0 .and. sw_hubbard == ON ) then
!     select case ( occmat_type )
!     case ( OccMat_Type1 )
!        if(flag_dm) call m_CD_den_mat(nfout,kv3)
!        call m_OP_occ_mat_ylm( nfout, pmode )
!     case ( OccMat_Type2 )
!        if(flag_dm) call m_CD_den_mat(nfout,kv3)
!        call m_OP_occ_mat_ao_kt(nfout, pmode)
!     end select

     if ( noncol ) then
        select case ( occmat_type )
        case ( OccMat_Type1 )
           if(flag_dm) call m_CD_den_mat(nfout,kv3)
           call m_OP_occ_mat_ylm_noncl( nfout, pmode )
        case ( OccMat_Type2 )
           if(flag_dm) call m_CD_den_mat(nfout,kv3)
           call m_OP_occ_mat_ao_kt_noncl(nfout, pmode)
        end select
     else
        select case ( occmat_type )
        case ( OccMat_Type1 )
           if(flag_dm) call m_CD_den_mat(nfout,kv3)
           call m_OP_occ_mat_ylm( nfout, pmode )
        case ( OccMat_Type2 )
           if(flag_dm) call m_CD_den_mat(nfout,kv3)
           call m_OP_occ_mat_ao_kt(nfout, pmode)
        end select
     endif
! ========================================================= 11.0
  end if
! =========================================================================-- 5.0

#ifdef __TIMER_SUB__
  call timer_end(735)
#endif

end subroutine Renewal_of_OccMat
