!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 606 $)
!
!  MODULE: m_ES_dos
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!      Further modification by T. Yamasaki   May/09/2004
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
module m_ES_dos
!     (m_ESdos)
! $Id: m_ES_dos.F90 606 2020-04-15 06:45:49Z ktagami $
!
! This module was originally coded by T. Yamasaki (FUJITSU Laboratories) in 2001.
! And this is transferred as match to PHASE by T. Yamasaki, 18th May. 2003.
!
  use m_Kpoints, only :              kv3, kv3_ek, qwgt,vkxyz_ek &
       &                           , np0,np2,ip20,iwt,ip2cub,nxyz_tetra,trmat &
       &                           , m_Kp_sample_mesh, qwgt_ek
  use m_Files, only :                nfout
!!$  use m_Files, only :                nfdos, nfout
  use m_Timing, only :               tstatc0_begin, tstatc0_end
  use m_Control_Parameters, only   : ekmode, ipridos, nspin, neg, af &
       &                            ,nwd_dos_window_width &
       &                            ,deltaE_dos,  variance_dos_GaussD &
       &                            ,sw_pdos, pdos_method, norbital &
       &                            ,maxorb, l_orb, t_orb, rc_orb, k_orb &
       &                            ,sw_orb_popu,dos_subroutine &
       &                            ,ipriinputfile, printable, sw_dos, dos_write_format &
       &                            ,sw_modified_kpoint_increment
  use m_Const_Parameters, only :     DP,Hartree,BUCS,EK,SCF, ALDOS, LAYERDOS, ON, OFF, TOTAL, PAI2 &
       &                           , WIDE, NARROW, eVunit
  use m_Parallelization, only :      MPI_CommGroup,map_ek,mype,map_e,map_k,myrank_e,myrank_k &
       &                            , ierr,np_e,map_z,ista_e,npes,mpi_kg_world &
       &                            , mpi_ge_world, nrank_e,ista_k,iend_k &
       &                            , nis_kv3_ek, ista_spin, iend_spin, myrank_spin, nrank_s
  use m_PseudoPotential, only :      nlmta_phi,nlmtt_phi &
       &                            ,m_PP_tell_iorb_ia_l_m_tau, qorb &
       &                            ,m_PP_tell_iorb_lmtt, ltp_phi, mtp_phi, taup_phi &
       &                            ,lmta_phi, nlmt_phi, ilmt_phi, iproj_phi, lmtt_phi
  use m_Nonlocal_Potential,  only : norm_phig
  use m_Crystal_Structure, only   : nopr, il, univol, altv
  use m_Electronic_Structure, only : occup_l, eko_l,eko_ek,totch, neordr &
       &                           , vbm, check_if_metalic_flag, metalic_system &
       &                           , compr_l, compi_l

!!$ASASASASAS
  use m_Kpoints,            only : k_symmetry
  use m_Const_Parameters,   only : GAMMA
  use m_Ionic_System,       only : iproj_group, natm, ityp, cps, if_pdos
!!$ASASASASAS


! =================================== added by K. Tagami ============ 11.0
  use m_Const_Parameters,   only : CMPLDP, zi, DIAG_SPIN_DENSITY_MATRIX, &
       &                           DIAG_CHARGE_DENSITY_MATRIX, DIAG_LS_with_t2g_octa, &
       &                           DIAG_LS, LOCAL_POINT_GROUP
  use m_Control_Parameters,   only : ndim_spinor, ndim_magmom, ndim_chgpot, noncol
! =================================================================== 11.0
  use m_ES_NonCollinear,    only : m_ES_DensMat_To_MagMom_porb, m_ES_set_Pauli_Matrix

! =================================== added by K. Tagami ============ 11.0
  use m_PseudoPotential, only :    modnrm, nloc, ntau
  use m_Parallelization, only :    ista_e,iend_e,istep_e, np_e
  use m_Control_Parameters,   only : kimg, hardpart_subroutine, &
       &                             calc_dos_magmom_contrib
  use m_Const_Parameters,   only : EXECUT, ELECTRON, DIRECT, YES, NO
  use m_PlaneWaveBasisSet, only :  kgp
  use m_FFT, only :                nfft, nfftp, nfftp_nonpara, fft_box_size_WF, &
       &                           fft_box_size_CD, &
       &                           m_FFT_alloc_WF_work, m_FFT_dealloc_WF_work, &
       &                           m_FFT_CD_inverse0,   m_FFT_WF, &
       &                           m_FFT_alloc_CD_box, m_FFT_dealloc_CD_box
  use m_Charge_Density, only :     chgq_l, chgq_enl, &
!       &                           m_CD_hardpart_sub_noncl, &
!       &                           m_CD_hardpart_sub2_noncl, &
       &                           m_CD_map_chgq_to_fft_box, &
       &                           m_CD_set_ylm_enl_etc, &
       &                           m_CD_dealloc_ylm_enl_etc, &
       &                           m_CD_map_chgqenl_to_fft_box_kt
!       &                           m_CD_map_valence_charge_to_fft_box, &
!       &                           m_CD_restore_chgq
!  use m_Electronic_Structure,only :m_ES_WF_in_Rspace
! =================================================================== 11.0

! =============== KT_add =============== 13.0E
  use m_Control_Parameters,  only : smearing_width_fdirac => width, width_tetra, &
       &                            orb_popu_method, sw_diagonalize_population, &
       &                            population_diag_mode
  use m_Electronic_Structure, only : efermi
! ====================================== 13.0E

  use m_OP_rotation,   only : m_OP_alloc_compri_rot,  m_OP_dealloc_compri_rot, &
       &                      m_OP_set_compri_rot, compr_rot_l, compi_rot_l, &
       &                      porb_rot_matrix_real,  porb_rot_matrix_cmplx
  use m_IterationNumbers, only : nkgroup, nk_in_the_process
  use mpi

  implicit none

  real(kind=DP),private ::                            Eminimum, Emaximum
  integer,private ::                                  nEwindows, nETailWindows
  integer,private,parameter ::                        nETailMax = 1000
  real(kind=DP),private ::                            ValenceBandMaximum
  real(kind=DP),private,parameter ::                  DeltaD = 1.d-18
  real(kind=DP),private,parameter ::                  DeltaDVBM = 1.d-4
  real(kind=DP),private ::                            sqrtDVI

!!$  real(kind=DP),public ::  deltaE_dos_GaussD = 1.d-4
!!$  real(kind=DP),public ::  variance_dos_GaussD = 1.d-6
!!$  integer,public ::        nwd_dos_window_width = 10

  real(kind=DP),private,allocatable,dimension(:,:) :: eko
  real(kind=DP),private,allocatable,dimension(:,:) :: dos
  real(kind=DP),private,allocatable,dimension(:,:) :: sumdos
  real(kind=DP),private,allocatable,dimension(:,:) :: dos_weight

! ===================== added by K. Tagami =============================== 11.0
  real(kind=DP),private,allocatable,dimension(:,:,:) :: dos_weight_noncl
! ======================================================================== 11.0

  integer, private ::                                 ndim_dos_weight1, ndim_dos_weight2

  !-- PDOS --
  real(kind=DP),private,allocatable,dimension(:,:,:,:) :: compr
  real(kind=DP),private,allocatable,dimension(:,:,:,:) :: compi
  real(kind=DP),private,allocatable,dimension(:,:) :: norm_phig_mpi
  real(kind=DP),private,allocatable,dimension(:,:,:) :: pdos
  real(kind=DP),private,allocatable,dimension(:,:,:) :: sumpdos

  character(len("pdos")),private,parameter :: tag_pdos       = "pdos"
  character(len("sw_pdos")),private,parameter :: tag_sw_pdos = "sw_pdos"
  character(len("sw_orb_popu")),private,parameter :: tag_sw_orb_popu = "sw_orb_popu"
  character(len("method")),private,parameter :: tag_method   = "method"
  character(len("element")),private,parameter :: tag_element = "element"
  character(len("orbitals")),private,parameter :: tag_orbitals = "orbitals"
  character(len("l")),private,parameter ::        tag_l    = "l"
  character(len("t")),private,parameter ::        tag_t    = "t"
  character(len("k")),private,parameter ::        tag_k    = "k"
  character(len("rc")),private,parameter ::       tag_rc   = "rc"
  character(len("projector")),private,parameter :: tag_projector = "projector"
  character(len("wavefunction")),private,parameter :: tag_wavefunction = "wavefunction"
  character(len("mulliken")),private,parameter :: tag_mulliken = "mulliken"

  character(len("use_rotated_compri")),private,parameter :: &
       &             tag_use_rotated_compri = "use_rotated_compri"
  integer :: use_rotated_compri = OFF

! -- pCOHP --
  character(len("pcohp")),private,parameter :: &
       &         tag_pcohp = "pcohp"
  character(len("sw_calc_pcohp")),private,parameter :: &
       &         tag_sw_calc_pcohp = "sw_calc_pcohp"
  character(len("atom1")),private,parameter :: &
       &         tag_atom1 = "atom1"
  character(len("max_distance_from_atom1")),private,parameter :: &
       &         tag_max_distance_from_atom1 = "max_distance_from_atom1"
!
  integer :: sw_calc_pcohp = off
  integer :: atom1_pcohp = 0
  real(kind=DP) :: max_distance_from_atom1_pcohp = 8.0d0  ! bohr

  real(kind=DP), private, allocatable :: pcohp(:,:,:,:)
  real(kind=DP), private, allocatable :: sum_pcohp(:,:,:,:)

  ! --- POSTPROCESSING ---
  character(len("Postprocessing")),private,parameter :: tag_postprocessing    = "postprocessing"

!  include 'mpif.h'
contains
  subroutine m_ESdos_alloc_dos_weight()
    if(.not.allocated(dos_weight)) then
       if(ekmode == ON) then
          allocate(dos_weight(neg,kv3_ek))
       else if(ekmode == OFF) then
          allocate(dos_weight(neg,kv3))
       end if
    end if
  end subroutine m_ESdos_alloc_dos_weight

  subroutine m_ESdos_dealloc_dos_weight()
    if(allocated(dos_weight)) deallocate(dos_weight)
  end subroutine m_ESdos_dealloc_dos_weight

! ========================= added by K. Tagami =============== 11.0
  subroutine m_ESdos_alloc_dos_wght_noncl()
    if (.not.allocated(dos_weight_noncl)) then
       if (ekmode == ON) then
          allocate( dos_weight_noncl(neg,kv3_ek,ndim_magmom) )
       else if (ekmode == OFF) then
          allocate( dos_weight_noncl(neg,kv3,ndim_magmom) )
       end if
       dos_weight_noncl = 0.0d0
    end if
  end subroutine m_ESdos_alloc_dos_wght_noncl

  subroutine m_ESdos_dealloc_dos_wght_noncl()
    if (allocated(dos_weight_noncl)) deallocate(dos_weight_noncl)
  end subroutine m_ESdos_dealloc_dos_wght_noncl
! ================================================================ 11.0

  subroutine m_ESdos_put_dos_weight(ne,nk,dos_weight_from)
    integer, intent(in) :: ne, nk
    real(kind=DP), intent(in), dimension(ne,nk) :: dos_weight_from
    if(.not.allocated(dos_weight)) &
         & call phase_error_with_msg(nfout,' dos_weight is not allocated <<m_ESdos_put_dos_weight>>',__LINE__,__FILE__)

    if(ekmode == ON) then
       if( ne /= neg .or. nk /= kv3_ek) then
          write(nfout,'(" !ldos ne /= neg or nk /= kv3_ek ",&
               & ": ne, neg, nk, kv3_ek = ",4i6," <<m_ESdos_put_dos_weight>>")') ne,neg,nk,kv3_ek
          call phase_error_with_msg(nfout,' ne /= neg or nk /= kv3_ek <<m_ESdos_put_dos_weight>>',__LINE__,__FILE__)
       end if
    else if(ekmode == OFF) then
       if( ne /= neg .or. nk /= kv3) then
          write(nfout,'(" !ldos ne /= neg or nk /= kv3 ",&
               & ": ne, neg, nk, kv3 = ",4i6," <<m_ESdos_put_dos_weight>>")') ne,neg,nk,kv3
          call phase_error_with_msg(nfout,' ne /= neg or nk /= kv3 <<m_ESdos_put_dos_weight>>',__LINE__,__FILE__)
       end if
    end if
    dos_weight(:,:) = dos_weight_from(:,:)
  end subroutine m_ESdos_put_dos_weight

! ============================== added by K. Tagami =============== 11.0
  subroutine m_ESdos_put_dos_weight_noncl(ne,nk,dos_weight_from)
    integer, intent(in) :: ne, nk
    real(kind=DP), intent(in), dimension(ne,nk,ndim_magmom) :: dos_weight_from

    if (.not.allocated(dos_weight_noncl)) &
         & call phase_error_with_msg(nfout,' dos_weight_noncl is not allocated <<m_ESdos_put_dos_weight_noncl>>'&
         &                          ,__LINE__,__FILE__)

    if (ekmode == ON) then
       if ( ne /= neg .or. nk /= kv3_ek) then
          write(nfout,'(" !ldos ne /= neg or nk /= kv3_ek ",&
               & ": ne, neg, nk, kv3_ek = ",4i6," <<m_ESdos_put_dos_weight_nocl>>")') &
               &  ne,neg,nk,kv3_ek
          call phase_error_with_msg(nfout,' ne /= neg or nk /= kv3_ek <<m_ESdos_put_dos_weight_noncl>>',__LINE__,__FILE__)
       end if

    else if(ekmode == OFF) then
       if ( ne /= neg .or. nk /= kv3) then
          write(nfout,'(" !ldos ne /= neg or nk /= kv3 ",&
               & ": ne, neg, nk, kv3 = ",4i6," <<m_ESdos_put_dos_weight_noncl>>")') &
               &    ne,neg,nk,kv3
          call phase_error_with_msg(nfout,' ne /= neg or nk /= kv3 <<m_ESdos_put_dos_weight_noncl>>',__LINE__,__FILE__)
       end if
    end if
    dos_weight_noncl(:,:,:) = dos_weight_from(:,:,:)

  end subroutine m_ESdos_put_dos_weight_noncl
! ============================================================== 11.0

  subroutine find_Erange(ek,neg,kv)
    integer, intent(in) :: neg,kv
    real(kind=DP), intent(in) :: ek(neg,kv)
    integer :: i, is
    logical :: found
    real(kind=DP) :: DeltaE, t
    real(kind=DP) :: derfc

    DeltaE = DeltaE_dos

    Eminimum = 9999.99d0;
    Emaximum  = -9999.99d0;

    sqrtdVI = 1.d0/dsqrt(2*Variance_dos_GaussD)
    t = DeltaE * sqrtdVI  ! t = DeltaE/dsqrt(2*Variance_dos_GaussD)

    if(ipridos >= 2) then
       write(nfout,'(" Attenuation check of a Gassuian function (find_Erange)")')
       write(nfout,'(" sqrtdVI = ",f10.5," DeltaE = ",f10.5, " t = ",f10.5 )')&
         & sqrtdVI, DeltaE, t
    end if

    i = 0; found = .false.
    do while(.not.found.and.i<=nETailMax)
       i = i + 1
       found = derfc(t*i) < DeltaD
    end do
    nETailWindows = i

    if(ipridos >= 2) write(nfout,'(" nETailWindows = ",i5)') nETailWindows

    Eminimum = minval(ek)
    Emaximum = maxval(ek)
!!$    write(nfout,'(" !! minval(ek), maxval(ek) = ",2f10.6)') Eminimum, Emaximum

    i = abs(Eminimum/DeltaE) + 1
    is = sign(1.d0,Eminimum)
    Eminimum = i*DeltaE*is - nETailWindows*DeltaE

    i = abs(Emaximum/DeltaE) + 1
    is = sign(1.d0,Emaximum)
    Emaximum = i*DeltaE*is + nETailWindows*DeltaE
    nEwindows = (Emaximum - Eminimum)/DeltaE + 1

    if(ipridos >= 2) then
       write(nfout,'(" Emaximum = ",f8.4, " Eminimum = ",f8.4 &
            & , " nEwindows = ",i8)') Emaximum, Eminimum, nEwindows+2*nETailWindows
    end if

  end subroutine find_Erange

! ================================= KT_add ======================== 13.0E
  subroutine find_Erange_fermidirac(ek,neg,kv)
    integer, intent(in) :: neg,kv
    real(kind=DP), intent(in) :: ek(neg,kv)
    integer :: i, is
    logical :: found
    real(kind=DP) :: DeltaE, t
    real(kind=DP) :: xx, yy
!
    real(kind=DP) :: threshold = 1.0E-4
!
    DeltaE = DeltaE_dos

    Eminimum = 9999.99d0;    Emaximum  = -9999.99d0;
    Eminimum = minval(ek);  Emaximum = maxval(ek)

    xx = 2.0D0 *log( 2.0D0 *sqrt(1.0D0 /threshold ) )  ! width of extent of df/dx
                                                       ! f: Fermi-Dirac fn
    xx = xx *smearing_width_fdirac

    nETailWindows = xx / DeltaE +1

    i = abs(Eminimum/DeltaE) +1;  is = sign(1.d0,Eminimum)
    Eminimum = i*DeltaE*is -nETailWindows*DeltaE

    i = abs(Emaximum/DeltaE) +1;  is = sign(1.d0,Emaximum)
    Emaximum = i*DeltaE*is + nETailWindows*DeltaE

    nEwindows = (Emaximum -Eminimum)/DeltaE +1
    if(ipridos >= 2) write(nfout,'(" nETailWindows = ",i5)') nETailWindows


  end subroutine find_Erange_fermidirac
! =========================================================== 13.0E

  subroutine alloc_eko_and_substitution(kv)
    use m_Parallelization, only : mpi_kg_world, mpi_ge_world
    integer, intent(in) :: kv
    integer ::             ik,ie,ip,iksnl
    integer :: iorb,iopr, nelm

! ==================== added by K. Tagami ================== 11.0
#ifdef forsafe
    integer :: ksym_ik
#endif
! ========================================================== 11.0
    integer :: ispin

    if ( sw_diagonalize_population == ON .and. use_rotated_compri == YES ) then
       call m_OP_alloc_compri_rot
       call m_OP_set_compri_rot
    endif

    allocate(eko(neg,kv));  eko = 0.d0
    if(sw_pdos == ON) then
       allocate(compr(neg,nlmta_phi,nopr,kv));  compr = 0.d0
       allocate(compi(neg,nlmta_phi,nopr,kv));  compi = 0.d0
       if(.not.allocated(norm_phig_mpi)) allocate(norm_phig_mpi(nlmtt_phi,kv/nspin))
       norm_phig_mpi=0.d0
    end if

! Substitution eko_l for eko in the order of thier values
    do ispin = ista_spin, iend_spin
!    do ik = 1, kv
    do ik = ispin, kv-nspin+ispin, nspin
       if(map_k(ik) /= myrank_k) cycle
! ============================= added by K. Tagami ============== 11.0
#ifdef forsafe
       ksym_ik = k_symmetry(ik)
#endif
! =============================================================== 11.0
       iksnl = (ik-1)/nspin + 1
       do ie = 1, neg
          ip = neordr(ie,ik)
          if(map_e(ip) == myrank_e) then
             eko(ie,ik) = eko_l(map_z(ip),ik)
             if(sw_pdos == ON) then
                if ( sw_diagonalize_population==ON .and. use_rotated_compri==YES ) then
                   compr(ie,1:nlmta_phi,1:nopr,ik) &
                        & = compr_rot_l(map_z(ip),1:nlmta_phi,1:nopr,ik)
#ifdef forsafe
                   if ( ksym_ik  /= GAMMA ) then
#else
                   if ( k_symmetry(ik) /= GAMMA ) then
#endif
                      compi(ie,1:nlmta_phi,1:nopr,ik) &
                           & = compi_rot_l(map_z(ip),1:nlmta_phi,1:nopr,ik)
                   else
                      compi(ie,1:nlmta_phi,1:nopr,ik) = 0.0d0
                   endif
                else
                   compr(ie,1:nlmta_phi,1:nopr,ik) &
                        & = compr_l(map_z(ip),1:nlmta_phi,1:nopr,ik)
#ifdef forsafe
                   if ( ksym_ik  /= GAMMA ) then
#else
                   if ( k_symmetry(ik) /= GAMMA ) then
#endif
                      compi(ie,1:nlmta_phi,1:nopr,ik) &
                           & = compi_l(map_z(ip),1:nlmta_phi,1:nopr,ik)
                   else
                      compi(ie,1:nlmta_phi,1:nopr,ik) = 0.0d0
                   endif
                endif
                norm_phig_mpi(1:nlmtt_phi,iksnl)  = norm_phig(1:nlmtt_phi,iksnl)
             end if
          end if
       end do
    end do
    end do

    if(npes >=2 ) then
       call mpi_allreduce(MPI_IN_PLACE,eko,neg*kv,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
       call mpi_allreduce(MPI_IN_PLACE,eko,neg*kv,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)
       if(sw_pdos == ON) then
          nelm = neg*nlmta_phi*nopr*kv
          call mpi_allreduce(MPI_IN_PLACE,compr,nelm,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
          call mpi_allreduce(MPI_IN_PLACE,compr,nelm,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)
          call mpi_allreduce(MPI_IN_PLACE,compi,nelm,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
          call mpi_allreduce(MPI_IN_PLACE,compi,nelm,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)
       !call mpi_allreduce(norm_phig_mpi,norm_phig_mpi2,nlmtt_phi*kv/nspin,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
          call mpi_allreduce(MPI_IN_PLACE,norm_phig_mpi,nlmtt_phi*kv/nspin,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)
          norm_phig_mpi = norm_phig_mpi / dble(nrank_s)
       end if
    end if
    if(ipridos >= 2) then
       ! -- writing --
       write(nfout,'(" !dos: eko <<alloc_eko_and_substitution>>")')
       do ik = 1, kv
          write(nfout,'(" !dos: ik = ",i5)') ik
          write(nfout,'(" !dos: ",10f8.4)') (eko(ie,ik),ie=1,neg)
       end do
       if(sw_pdos == ON) then
          do ik = 1, kv
             iksnl = (ik-1)/nspin + 1
             do iopr=1,nopr
                do iorb=1,nlmta_phi
                   write(nfout,'(" !dos: ik=",i5," iopr=",i5," iorb=",i5)')  ik,iopr,iorb
                   write(nfout,'(" !dos compr: ",10f8.4)') compr(1:neg,iorb,iopr,ik)
!!$ASASASASAS
!!$                write(nfout,'(" !dos compi: ",10f8.4)') &
!!$                              & compi(1:neg,iorb,iopr,ik)
! =========================== modified by K. Tagami ========= 11.0
!                if ( k_symmetry(ik) /= GAMMA ) then
#ifdef forsafe
                   if ( ksym_ik /= GAMMA ) then
#else
                      if ( k_symmetry(ik) /= GAMMA ) then
#endif
! =========================================================== 11.0

                         write(nfout,'(" !dos compi: ",10f8.4)') compi(1:neg,iorb,iopr,ik)
                      endif
!!$ASASASASAS
                   end do
                end do
             write(nfout,'(" !dos norm_phig: ",10f8.4)') norm_phig_mpi(1:nlmtt_phi,iksnl)
          end do
       end if
    end if
  end subroutine alloc_eko_and_substitution

! ========================= added by K. Tagami =================== 11.0
  subroutine alloc_eko_and_substit_noncl(kv)
    integer, intent(in) :: kv
    integer ::             ik,ie,ip,iksnl
    integer :: iorb,iopr
    integer :: is, nelm
    integer :: num_orbitals, kfac

#ifdef forsafe
    integer :: ksym_ik
#endif

    kfac = 1
    if ( sw_diagonalize_population == ON ) then
       if ( population_diag_mode /= DIAG_CHARGE_DENSITY_MATRIX &
            & .and. population_diag_mode /= LOCAL_POINT_GROUP ) then
          kfac = ndim_spinor
       endif
    endif
    num_orbitals = nlmta_phi *kfac

    if ( sw_diagonalize_population == ON .and. use_rotated_compri == YES ) then
       call m_OP_alloc_compri_rot
       call m_OP_set_compri_rot
    endif

    allocate(eko(neg,kv/ndim_spinor));  eko = 0.d0

    if(sw_pdos == ON) then
       allocate(compr(neg,num_orbitals,nopr,kv));  compr = 0.d0
       allocate(compi(neg,num_orbitals,nopr,kv));  compi = 0.d0
       allocate(norm_phig_mpi(nlmtt_phi,kv/ndim_spinor));
       norm_phig_mpi=0.d0
    end if

! Substitution eko_l for eko in the order of thier values
    do ik = 1, kv, ndim_spinor

! ============================= added by K. Tagami ============== 11.0
#ifdef forsafe
       ksym_ik = k_symmetry(ik)
#endif
! =============================================================== 11.0

       if(map_k(ik) /= myrank_k) cycle
       iksnl = (ik-1)/ndim_spinor + 1

       do ie = 1, neg
          ip = neordr(ie,ik)

          if (map_e(ip) == myrank_e) then
             eko( ie,iksnl ) = eko_l( map_z(ip),ik )

             if(sw_pdos == ON) then
                if ( sw_diagonalize_population==ON .and. use_rotated_compri==YES ) then
                   Do is=1, ndim_spinor
                      compr(ie,1:num_orbitals,1:nopr,ik+is-1) &
                           & = compr_rot_l(map_z(ip),1:num_orbitals,1:nopr,ik+is-1)
                      compi(ie,1:num_orbitals,1:nopr,ik+is-1) &
                           & = compi_rot_l(map_z(ip),1:num_orbitals,1:nopr,ik+is-1)
                   End do
                else
                   Do is=1, ndim_spinor
                      compr(ie,1:num_orbitals,1:nopr,ik+is-1) &
                           & = compr_l(map_z(ip),1:num_orbitals,1:nopr,ik+is-1)
                      compi(ie,1:num_orbitals,1:nopr,ik+is-1) &
                           & = compi_l(map_z(ip),1:num_orbitals,1:nopr,ik+is-1)
                   End do
                endif
                norm_phig_mpi(1:nlmtt_phi,iksnl)  = norm_phig(1:nlmtt_phi,iksnl)
             end if

          end if
       end do
    end do

    if (npes >=2 ) then
       call mpi_allreduce( MPI_IN_PLACE,eko,neg*kv/ndim_spinor,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       if (sw_pdos == ON) then
          nelm = neg*num_orbitals*nopr*kv
          call mpi_allreduce( MPI_IN_PLACE,compr,nelm,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
          call mpi_allreduce( MPI_IN_PLACE,compi,nelm,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
          nelm = nlmtt_phi*kv/ndim_spinor
          call mpi_allreduce( MPI_IN_PLACE,norm_phig_mpi,nelm,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
          norm_phig_mpi = norm_phig_mpi / dble(nrank_s)
       end if

    end if

    if(ipridos >= 2) then
       ! -- writing --
       write(nfout,'(" !dos: eko <<alloc_eko_and_substi_noncl>>")')

       do ik = 1, kv, ndim_spinor
          iksnl = (ik-1)/ndim_spinor + 1
          write(nfout,'(" !dos: ik = ",i5)') ik
          write(nfout,'(" !dos: ",10f8.4)') (eko(ie,iksnl),ie=1,neg)
       end do

       if(sw_pdos == ON) then
          do ik = 1, kv, ndim_spinor
             iksnl = (ik-1)/ndim_spinor + 1

             do iopr=1,nopr
                do iorb=1,num_orbitals
                   write(nfout,'(" !dos: ik=",i5," iopr=",i5," iorb=",i5)')&
                        &  ik,iopr,iorb
                   write(nfout,'(" !dos compr: ",10f8.4)') compr(1:neg,iorb,iopr,ik)
                   write(nfout,'(" !dos compi: ",10f8.4)') compi(1:neg,iorb,iopr,ik)
                end do
             end do
             write(nfout,'(" !dos norm_phig: ",10f8.4)')  norm_phig_mpi(1:nlmtt_phi,iksnl)
          end do
       end if
    end if

  end subroutine alloc_eko_and_substit_noncl
! =================================================================== 11.0

  subroutine m_ES_dos_alloc_compr_compi_ek()
    allocate(compr(neg,nlmta_phi,nopr,kv3_ek));  compr = 0.d0
    allocate(compi(neg,nlmta_phi,nopr,kv3_ek));  compi = 0.d0
    if(allocated(norm_phig_mpi)) deallocate(norm_phig_mpi)
    allocate(norm_phig_mpi(nlmtt_phi,kv3_ek/nspin))
    norm_phig_mpi=0.d0
  end subroutine m_ES_dos_alloc_compr_compi_ek

  subroutine m_ES_dos_pdos_per_k()
    integer :: ik,ikn,is,ikt
    if(sw_modified_kpoint_increment == ON) then
      ik=myrank_k
      do is=1,nspin
        ikn = nspin*(nis_kv3_ek(ik)-1)+(nkgroup-1)*nspin+is
        if(ikn>kv3_ek) cycle
        call pdos_per_k(ik+1,ikn)
      enddo
    else
      do ik = 1, kv3
         if(map_k(ik) /= myrank_k) cycle
         ikn = nk_in_the_process+ik-1
         if(ikn>kv3_ek) cycle
         call pdos_per_k(ik,ikn)
      end do
    endif

    contains

    subroutine pdos_per_k(ik,ikn)
      integer, intent(in) :: ik, ikn
      integer :: iksnl, iksn, ip
      integer :: ie,ib,jb,ibo,jbo
      integer,allocatable,dimension(:) :: neordr_ek
      integer,parameter :: delta = 1.d-12
      iksn  = (ik-1)/nspin + 1
      iksnl = (ikn-1)/nspin + 1
      allocate(neordr_ek(neg))
      neordr_ek = (/(ib,ib=1,neg)/)
      do ib = 1,neg-1
         do jb = ib+1, neg
            ibo = neordr_ek(ib)
            jbo = neordr_ek(jb)
            if(eko_ek(jbo,ikn) < eko_ek(ibo,ikn)-delta) then
               neordr_ek(jb) = ibo
               neordr_ek(ib) = jbo
            end if
         end do
      end do
      do ie = 1, neg
         ip = neordr_ek(ie)
         if(map_e(ip) == myrank_e) then
            if ( sw_diagonalize_population==ON .and. use_rotated_compri==YES ) then
               compr(ie,1:nlmta_phi,1:nopr,ikn) &
                    & = compr_rot_l(map_z(ip),1:nlmta_phi,1:nopr,ik)
               if ( k_symmetry(ik) /= GAMMA ) then
                  compi(ie,1:nlmta_phi,1:nopr,ikn) &
                     & = compi_rot_l(map_z(ip),1:nlmta_phi,1:nopr,ik)
               else
                  compi(ie,1:nlmta_phi,1:nopr,ikn) = 0.0d0
               endif
            else
               compr(ie,1:nlmta_phi,1:nopr,ikn) &
                    & = compr_l(map_z(ip),1:nlmta_phi,1:nopr,ik)
               if ( k_symmetry(ik) /= GAMMA ) then
                  compi(ie,1:nlmta_phi,1:nopr,ikn) &
                       & = compi_l(map_z(ip),1:nlmta_phi,1:nopr,ik)
               else
                  compi(ie,1:nlmta_phi,1:nopr,ikn) = 0.0d0
               endif
            endif
         end if
      end do
      norm_phig_mpi(1:nlmtt_phi,iksnl)  = norm_phig(1:nlmtt_phi,iksn)
      deallocate(neordr_ek)
    end subroutine pdos_per_k
  end subroutine m_ES_dos_pdos_per_k

  subroutine m_ES_dos_gather_compr_compi_ek()
    integer :: nelm
    integer :: ib,ilm,iopr,ik


    if(npes >=2) then
       nelm = neg*nlmta_phi*nopr*kv3_ek
       call mpi_allreduce(MPI_IN_PLACE,compr,nelm,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
       call mpi_allreduce(MPI_IN_PLACE,compr,nelm,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)
       call mpi_allreduce(MPI_IN_PLACE,compi,nelm,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
       call mpi_allreduce(MPI_IN_PLACE,compi,nelm,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)
       call mpi_allreduce(MPI_IN_PLACE,norm_phig_mpi,nlmtt_phi*kv3_ek/nspin, &
       &    mpi_double_precision,mpi_sum,mpi_ge_world,ierr)
       norm_phig_mpi = norm_phig_mpi / dble(nrank_s)
    end if
  end subroutine m_ES_dos_gather_compr_compi_ek

  subroutine pdos_ek()
    use m_Control_Parameters, only : af
    use m_Files, only : m_Files_open_nfzaj_kall, nfzaj_kall
    use m_IterationNumbers, only : m_Iter_reset_nk_in_the_process
    use m_ES_IO, only : m_ESIO_rd_next_wfs_ek
    use m_ES_nonlocal, only : m_ES_betar_dot_WFs_4_each_k_3D, m_ES_phir_dot_WFs_3D
    use m_Electronic_Structure, only : m_ES_sym_comp
    integer :: nk, is, ik

    call m_Files_open_nfzaj_kall()
    rewind(nfzaj_kall)
    call m_Iter_reset_nk_in_the_process()
    compr = 0.d0
    compi = 0.d0
    norm_phig_mpi = 0.d0
    do nk=1, kv3_ek, kv3
      call KpointNumber_Setting2()
      call Preparation_ek()
      call Preparation_for_mpi_ek()
      call epsmain_reallocate()
      call PseudoPotential_ek()
      do is=1, nspin, af+1
        do ik=is, kv3+is-nspin, nspin
          call m_ESIO_rd_next_wfs_ek(ik, nfout, nfzaj_kall)
          if(map_k(ik) /= myrank_k) cycle
          call m_ES_betar_dot_WFs_4_each_k_3D(nfout, ik)
        enddo
      enddo
      call m_ES_phir_dot_WFs_3D( nfout )
      call m_ES_sym_comp(nfout)
      call m_ES_dos_pdos_per_k()
    enddo
  end subroutine pdos_ek

  subroutine alloc_eko_and_substitution_ek(kv3_ek)
    integer, intent(in) :: kv3_ek
    integer :: ik,ie,ib,jb,ibo,jbo
    integer,allocatable,dimension(:) :: neordr_ek
    real(kind=DP),parameter :: delta = 1.d-12

    allocate(eko(neg,kv3_ek))
    allocate(neordr_ek(neg))
    do ik = 1, kv3_ek
       neordr_ek = (/(ib,ib=1,neg)/)
       do ib = 1,neg-1
          do jb = ib+1, neg
             ibo = neordr_ek(ib)
             jbo = neordr_ek(jb)
             if(eko_ek(jbo,ik) < eko_ek(ibo,ik)-delta) then
                neordr_ek(jb) = ibo
                neordr_ek(ib) = jbo
             end if
          end do
       end do
! Substitution eko_ek for eko in the order of thier values
       do ie = 1, neg
          eko(ie,ik) = eko_ek(neordr_ek(ie),ik)
       end do
    end do
!!$    write(nfout,'(" !! eko (alloc_eko_and_substitution_ek)")')
!!$    do ik = 1, kv3_ek
!!$       write(nfout,'(" !! ik = ",i5, "(",3f8.4,")")') ik, (vkxyz_ek(ib,ik,BUCS),ib=1,3)
!!$       write(nfout,'("(eko_ek): ",8f10.6)') (eko_ek(ie,ik),ie=1,neg)
!!$       write(nfout,'("(eko):    ",8f10.6)') (eko(ie,ik),ie=1,neg)
!!$    end do
    deallocate(neordr_ek)
  end subroutine alloc_eko_and_substitution_ek

! ============================== added by K. Tagami ============= 11.0
  subroutine alloc_eko_and_substit_ek_noncl(kv3_ek)
    integer, intent(in) :: kv3_ek
    integer :: ik,ie,ib,jb,ibo,jbo
    integer :: iksnl
    integer :: ikt

    integer,allocatable,dimension(:) :: neordr_ek
    real(kind=DP),parameter :: delta = 1.d-12

    allocate(eko(neg,kv3_ek/ndim_spinor))
    allocate(neordr_ek(neg))

    do ik = 1, kv3_ek, ndim_spinor
       iksnl = ( ik -1 ) /ndim_spinor + 1

       neordr_ek = (/(ib,ib=1,neg)/)
       do ib = 1,neg-1
          do jb = ib+1, neg
             ibo = neordr_ek(ib)
             jbo = neordr_ek(jb)
             if(eko_ek(jbo,ik) < eko_ek(ibo,ik)-delta) then
                neordr_ek(jb) = ibo
                neordr_ek(ib) = jbo
             end if
          end do
       end do
! Substitution eko_ek for eko in the order of thier values
       if(ik+kv3-1 > kv3_ek) then
       do ie = 1, neg
          eko(ie,iksnl) = eko_ek(neordr_ek(ie),ik)
       end do
       else
          ikt = myrank_k*kv3+mod(ik,kv3)
          eko(ie,iksnl) = eko_ek(neordr_ek(ie),ik)
       endif
    end do

    deallocate(neordr_ek)
  end subroutine alloc_eko_and_substit_ek_noncl
! ===================================================================== 11.0

  subroutine dealloc_eko()
    deallocate(eko)
    if(sw_pdos == ON) then
       deallocate(compr,compi)
       if (allocated(compr_rot_l) ) deallocate(compr_rot_l)
       if (allocated(compi_rot_l) ) deallocate(compi_rot_l)
    end if
  end subroutine dealloc_eko

  subroutine alloc_dos(i0,icomponent)
    integer,intent(in) :: i0,icomponent
    integer :: num_orbitals, kfac

    if ( noncol ) then
       allocate(dos(i0:nEwindows,ndim_magmom)); dos = 0.d0
       allocate(sumdos(i0:nEwindows,ndim_magmom)); sumdos = 0.d0
       if(icomponent == TOTAL .and. sw_pdos == ON) then
          kfac = 1
          if ( sw_diagonalize_population == ON ) then
             if ( population_diag_mode /= DIAG_CHARGE_DENSITY_MATRIX &
                  & .and. population_diag_mode /= LOCAL_POINT_GROUP ) then
                kfac = ndim_spinor
             endif
          endif
          num_orbitals = nlmta_phi *kfac
          allocate(pdos(i0:nEwindows,num_orbitals,ndim_magmom)); pdos = 0.d0
          allocate(sumpdos(i0:nEwindows,num_orbitals,ndim_magmom)); sumpdos = 0.d0
       end if

    else
       allocate(dos(i0:nEwindows,nspin)); dos = 0.d0
       allocate(sumdos(i0:nEwindows,nspin)); sumdos = 0.d0
       if(icomponent == TOTAL .and. sw_pdos == ON) then
          allocate(pdos(i0:nEwindows,nlmta_phi,nspin)); pdos = 0.d0
          allocate(sumpdos(i0:nEwindows,nlmta_phi,nspin)); sumpdos = 0.d0
       end if
    endif

  end subroutine alloc_dos

  subroutine dealloc_dos()
    deallocate(dos)
    deallocate(sumdos)
    if(sw_pdos == ON) then
       if(allocated(pdos)) deallocate(pdos)
       if(allocated(sumpdos)) deallocate(sumpdos)
    end if
  end subroutine dealloc_dos

  subroutine make_dos_with_GaussianDistrib( kv, iwsc, kpt_weight )
    integer, intent(in) :: kv, iwsc
    real(kind=DP), intent(in) :: kpt_weight( kv )

    integer ::             i, ik, is, ie, id, ispin, iorb, iopr, lmtt, iksnl
    real(kind=DP) ::       Es, e, El, Eu, tl, tu, w, DeltaE
    real(kind=DP) ::       derf
    real(kind=DP), allocatable :: porb(:)

    allocate( porb(nlmta_phi) );  porb = 0.0d0

    DeltaE = DeltaE_dos
    Es = Eminimum - DeltaE*0.5d0

    if(ipridos >= 2) write(nfout,'(" !! Es, DeltaE = ",2d20.10)') Es, DeltaE
    dos = 0.d0; sumdos = 0.d0
    if(iwsc == TOTAL .and. sw_pdos == ON) then
       pdos=0.d0; sumpdos=0.d0
    end if
    do ispin = 1, nspin
       do ik = ispin, kv, nspin
          iksnl = (ik-1)/nspin + 1
          do i = 1, neg
             w = 1.d0
             if(iwsc >= 1 ) w = dos_weight(i,ik)
             e = eko(i,ik) - nETailWindows*DeltaE - Es
             is = e/DeltaE
             ie = is + 2*nETailWindows
             if(is < 0) then
                if(ipridos >= 2) write(nfout,'(" is = ",i5," < 0")') is
                is = 0
             end if
             if(ie >= nEWindows ) then
                if(ipridos >= 2) write(nfout,'(" ie = ",i5, " > nEWindows")')ie
                ie = nEWindows-1
             end if
             if(ipridos >= 2) then
                write(nfout,'(" !! eko(",i3,",",i3,") = ",f10.6," is, ie = ",2i7, 2f10.6)') &
                     & i,ik, eko(i,ik),is,ie, Es+is*DeltaE, Es+ie*DeltaE
             end if

             do id = is, ie
                El = Es + id*DeltaE
                Eu = El + DeltaE
                tl = (El - eko(i,ik))*sqrtdVI
                tu = (Eu - eko(i,ik))*sqrtdVI
                !  d = (derf(tl) - derf(tu))*0.5d0/(DeltaE*kv)
                dos(id+1,ispin) = dos(id+1,ispin) &
                     & + w *2 *(derf(tu) - derf(tl)) *0.5d0 /DeltaE *kpt_weight(ik)
             end do
             if(iwsc == TOTAL .and. sw_pdos == ON) then
                w= 1.0d0;   porb = 0.0d0

                if ( sw_diagonalize_population == ON ) then
                   if ( use_rotated_compri == YES ) then
                      call set_porb_norot( porb )
                   else
                      call set_porb_rot_mode1( porb )
                   endif
                else
                   call set_porb_norot( porb )
                endif

                do iorb = 1, nlmta_phi
                   do id = is, ie
                      El = Es +id*DeltaE;    Eu = El +DeltaE
                      tl = (El -eko(i,ik))*sqrtdVI
                      tu = (Eu -eko(i,ik))*sqrtdVI
                      !  d = (derf(tl) - derf(tu))*0.5d0/(DeltaE*kv)
                      pdos(id+1,iorb,ispin) = pdos(id+1,iorb,ispin) &
                           & + porb(iorb) *w *2 *(derf(tu) - derf(tl)) &
                           &                *0.5d0 /DeltaE *kpt_weight(ik)
                   end do
                end do
             end if
          end do
       end do
       sumdos(1,ispin) = dos(1,ispin)*DeltaE
       do id = 1, nEWindows-1
          sumdos(id+1,ispin) = sumdos(id,ispin) + dos(id+1,ispin)*DeltaE
       end do
       if(iwsc == TOTAL .and. sw_pdos == ON) then
          do iorb = 1,nlmta_phi
             sumpdos(1,iorb,ispin) = pdos(1,iorb,ispin)*DeltaE
             do id = 1, nEWindows-1
                sumpdos(id+1,iorb,ispin) = sumpdos(id,iorb,ispin) + pdos(id+1,iorb,ispin)*DeltaE
             end do
          end do
       end if
    end do
    deallocate( porb )

  contains
    subroutine set_porb_norot( porb_out )
      real(kind=DP), intent(out) :: porb_out(nlmta_phi)

      integer :: iorb, iopr, lmtt
      real(kind=DP) :: ctemp

      Do iorb = 1, nlmta_phi
         call m_PP_tell_iorb_lmtt(iorb,lmtt)

         ctemp = 0.0d0
         if ( k_symmetry(ik) == GAMMA ) then
            do iopr=1,nopr
               ctemp = ctemp +compr(i,iorb,iopr,ik)**2 /2.0 &
                    &    *(1.d0+qorb(iorb)/(norm_phig_mpi(lmtt,iksnl)*2.) )
            end do
         else
            do iopr=1,nopr
               ctemp = ctemp + ( compr(i,iorb,iopr,ik)**2 &
                    &           +compi(i,iorb,iopr,ik)**2 ) &
                    &     *( 1.d0+qorb(iorb)/norm_phig_mpi(lmtt,iksnl) )
            end do
         endif
         porb_out(iorb) = porb_out(iorb) +ctemp /dble(nopr)
      End do
    end subroutine set_porb_norot

    subroutine set_porb_rot_mode1( porb_out )
      real(kind=DP), intent(out) :: porb_out(nlmta_phi)

      integer :: lmax, mmax, iopr
      integer :: ia, it, lmt, il, im, ip, iorb, lmtt
      integer :: lmt1, lmt2, im1, im2, iorb1, iorb2, lmtt1, lmtt2, im3
      integer :: tau, tau1, tau2
      real(kind=DP) :: c1, c2, ctmp
      real(kind=DP), allocatable :: dm(:,:,:,:)
      real(kind=DP), allocatable :: porb0(:,:,:)

      lmax = nloc;   mmax = 2*lmax -1

      allocate( dm( ntau, lmax, mmax, mmax ) )
      allocate( porb0( ntau, lmax, mmax ) )

      do ia=1,natm
         it = ityp(ia)
         if(iproj_group(ia) == 0) cycle

         dm = 0.0d0;  porb0 = 0.d00

         ! diagonal part
         do lmt=1,ilmt_phi(it)
            il = ltp_phi(lmt,it); im  = mtp_phi(lmt,it);
            tau = taup_phi(lmt,it)
            ip  = iproj_phi(lmt,it)
            iorb = lmta_phi(lmt,ia); lmtt = lmtt_phi(lmt,it)

            ctmp = 0.0d0
            if ( k_symmetry(ik) == GAMMA ) then
               do iopr=1, nopr
                  ctmp = ctmp + compr(i,iorb,iopr,ik)**2 /2.0 &
                       &      *(1.d0+qorb(iorb)/(norm_phig_mpi(lmtt,iksnl)*2.) )
               end do
            else
               do iopr=1, nopr
                  ctmp = ctmp + ( compr(i,iorb,iopr,ik)**2 &
                       &         +compi(i,iorb,iopr,ik)**2 ) &
                       &      *(1.d0+qorb(iorb)/norm_phig_mpi(lmtt,iksnl))
               end do
            endif
            dm(tau,il,im,im) = dm(tau,il,im,im) +ctmp/dble(nopr)
         end do
         ! non-diagonal part
         do lmt2=1,ilmt_phi(it)
            do lmt1=1,ilmt_phi(it)
               if ( lmt1 == lmt2 ) cycle
               if (ltp_phi(lmt1,it) /= ltp_phi(lmt2,it)) cycle
               if (taup_phi(lmt1,it) /= taup_phi(lmt2,it)) cycle

               il = ltp_phi(lmt1,it);     tau = taup_phi(lmt1,it)
               ip  = iproj_phi(lmt1,it)
               im1   = mtp_phi(lmt1,it);     im2 = mtp_phi(lmt2,it);
               iorb1 = lmta_phi(lmt1,ia);  iorb2 = lmta_phi(lmt2,ia)
               lmtt1 = lmtt_phi(lmt1,it);  lmtt2 = lmtt_phi(lmt2,it)

               c1 = ( qorb(iorb1) +qorb(iorb2) ) /2.0d0      ! approx
               c2 = sqrt(norm_phig_mpi(lmtt1,iksnl)*norm_phig_mpi(lmtt2,iksnl))

               ctmp = 0.0d0
               if(k_symmetry(ik) == GAMMA) then ! dame
                  do iopr = 1,nopr
                     ctmp = ctmp &
                          & + compr(i,iorb1,iopr,ik) &
                          &  *compr(i,iorb2,iopr,ik) /2.d0 &
                          &  *( 1.0d0 +c1 /( c2 *2.d0) )
                  end do
               else
                  do iopr = 1,nopr
                     ctmp = ctmp &
                          & + ( compr(i,iorb1,iopr,ik) &
                          &    *compr(i,iorb2,iopr,ik) &
                          &    +compi(i,iorb1,iopr,ik) &
                          &    *compi(i,iorb2,iopr,ik) )&
                          &  *( 1.0d0 +c1 /c2 )
                  end do
               endif
               dm(tau,il,im1,im2) = dm(tau,il,im1,im2) +ctmp/dble(nopr)
            end do
         end do

         Do tau=1, ntau
            Do il=1, lmax
               Do im1=1, mmax
                  ctmp = 0.0d0
                  Do im2=1, mmax
                     Do im3=1, mmax
                        ctmp = ctmp &
                             &    +porb_rot_matrix_real( ia,il,im2,im1 )  &
                             &       *dm(tau,il,im2,im3) &
                             &       *porb_rot_matrix_real( ia,il,im3,im1 )
                     End Do
                  End Do
                  porb0( tau, il, im1 ) = ctmp
               End DO
            End DO
         End Do

         Do lmt=1,ilmt_phi(it)
            il = ltp_phi(lmt,it);  im  = mtp_phi(lmt,it);  iorb = lmta_phi(lmt,ia)
            tau = taup_phi(lmt,it)
            porb_out( iorb ) = porb0( tau, il,im )
         ENd do
      End do
      deallocate( dm );  deallocate( porb0 )
    end subroutine set_porb_rot_mode1

  end subroutine make_dos_with_GaussianDistrib

! ==================================== added by K. Tagami ============== 11.0
  subroutine mkdos_with_GaussDistrib_noncl( kv, iwsc, kpt_weight )
    integer, intent(in) :: kv, iwsc
    real(kind=DP), intent(in) :: kpt_weight( kv )

    integer :: ik, iksnl, i, is, ie, istmp, id, iorb
    real(kind=DP) ::  Es, e, El, Eu, tl, tu, w, DeltaE
    real(kind=DP) ::  derf

    integer :: kfac, num_orbitals, ismax
    real(kind=DP), allocatable :: porb(:,:)
    complex(kind=CMPLDP), allocatable :: porb_ssrep(:,:)

    kfac = 1
    if ( sw_diagonalize_population == ON ) then
       if ( population_diag_mode /= DIAG_CHARGE_DENSITY_MATRIX &
            & .and. population_diag_mode /= LOCAL_POINT_GROUP ) then
          kfac = ndim_spinor
       endif
    endif
    num_orbitals = nlmta_phi *kfac

    allocate( porb(num_orbitals, ndim_magmom ) );  porb = 0.0d0
    allocate( porb_ssrep(num_orbitals, ndim_chgpot ) );  porb_ssrep = 0.0d0

    DeltaE = DeltaE_dos
    Es = Eminimum - DeltaE*0.5d0

    if (ipridos >= 2) write(nfout,'(" !! Es, DeltaE = ",2d20.10)') Es, DeltaE
    dos = 0.d0; sumdos = 0.d0
    if (iwsc == TOTAL .and. sw_pdos == ON) then
       pdos=0.d0; sumpdos=0.d0
    end if

    do ik = 1, kv, ndim_spinor
       iksnl = (ik-1)/ndim_spinor + 1
       do i = 1, neg

          e = eko(i,iksnl) - nETailWindows*DeltaE - Es

          is = e/DeltaE
          ie = is + 2*nETailWindows
          if(is < 0) then
             if(ipridos >= 2) write(nfout,'(" is = ",i5," < 0")') is
             is = 0
          end if
          if(ie >= nEWindows ) then
             if(ipridos >= 2) write(nfout,'(" ie = ",i5, " > nEWindows")')ie
             ie = nEWindows-1
          end if
          if(ipridos >= 2) then
             write(nfout,'(" !! eko(",i3,",",i3,") = ",f10.6," is, ie = ",2i7, 2f10.6)') &
                  & i,ik, eko(i,iksnl),is,ie, Es+is*DeltaE, Es+ie*DeltaE
          end if

          if ( calc_dos_magmom_contrib == YES ) then
             ismax = ndim_magmom
          else
             ismax = 1
          endif

          Do istmp=1, ismax
             w = dos_weight_noncl(i,ik,istmp)

             do id = is, ie
                El = Es + id*DeltaE
                Eu = El + DeltaE
                tl = (El - eko(i,iksnl))*sqrtdVI
                tu = (Eu - eko(i,iksnl))*sqrtdVI
                !  d = (derf(tl) - derf(tu))*0.5d0/(DeltaE*kv)
                dos(id+1,istmp) = dos(id+1,istmp) &
                     & + w *(derf(tu) - derf(tl)) *0.5d0 /DeltaE *kpt_weight(ik)
             end do
          End do

          if (iwsc == TOTAL .and. sw_pdos == ON) then
             w = 1.0d0;   porb_ssrep = 0.0d0

             if ( sw_diagonalize_population == ON ) then
                if ( use_rotated_compri == YES ) then
                   if ( population_diag_mode == DIAG_CHARGE_DENSITY_MATRIX &
                        &  .or. population_diag_mode == LOCAL_POINT_GROUP ) then
                      call set_porb_default( porb )
                   else
                      call set_porb_by_spinmixed_basis( porb )
                   endif
                else
                   if ( population_diag_mode == DIAG_CHARGE_DENSITY_MATRIX &
                        &  .or. population_diag_mode == LOCAL_POINT_GROUP ) then
                      call set_porb_rot_mode1( porb )
                   else
                      call set_porb_rot_mode2( porb )
                   endif
                endif
             else
                call set_porb_default( porb )
             endif

             Do iorb = 1, num_orbitals
                do id = is, ie
                   El = Es + id*DeltaE
                   Eu = El + DeltaE
                   tl = (El - eko(i,iksnl))*sqrtdVI
                   tu = (Eu - eko(i,iksnl))*sqrtdVI
                   !  d = (derf(tl) - derf(tu))*0.5d0/(DeltaE*kv)
                   pdos(id+1,iorb,:) = pdos(id+1,iorb,:) &
                        & + porb(iorb,:) *w *(derf(tu) - derf(tl)) &
                        &                *0.5d0 /DeltaE *kpt_weight(ik)
                end do
             End do

          end if
       end do
    end do

    sumdos(1,:) = dos(1,:)*DeltaE

    do id = 1, nEWindows-1
       sumdos(id+1,:) = sumdos(id,:) + dos(id+1,:)*DeltaE
    end do

    if(iwsc == TOTAL .and. sw_pdos == ON) then
       do iorb = 1,num_orbitals
          sumpdos(1,iorb,:) = pdos(1,iorb,:)*DeltaE
          do id = 1, nEWindows-1
             sumpdos(id+1,iorb,:) = sumpdos(id,iorb,:) + pdos(id+1,iorb,:)*DeltaE
          end do
       end do
    end if
    deallocate( porb ); deallocate( porb_ssrep )

  contains

    subroutine set_porb_default( porb_out )
      real(kind=DP), intent(out) :: porb_out(num_orbitals, ndim_magmom )

      integer :: iorb, is1, is2, istmp, iopr, lmtt
      complex(kind=CMPLDP) :: ztemp, z1, z2

      Do iorb = 1, nlmta_phi
         call m_PP_tell_iorb_lmtt(iorb,lmtt)

         Do is1=1, ndim_spinor
            Do is2=1, ndim_spinor
               istmp = ( is1 -1 )*ndim_spinor + is2

               ztemp = 0.0d0
               do iopr=1,nopr
                  z1 = dcmplx( compr(i,iorb,iopr,ik+is1-1 ), &
                       &       compi(i,iorb,iopr,ik+is1-1 ) )
                  z2 = dcmplx( compr(i,iorb,iopr,ik+is2-1 ), &
                       &       compi(i,iorb,iopr,ik+is2-1 ) )
                  ztemp = ztemp + z1 *conjg(z2) &
                       &      *( 1.d0+qorb(iorb)/norm_phig_mpi(lmtt,iksnl) )
               end do

               porb_ssrep(iorb,istmp) = porb_ssrep(iorb,istmp) &
                    &                   + ztemp /dble(nopr)
            End do
         End do
      End do
      call m_ES_DensMat_To_MagMom_porb( nlmta_phi, porb_ssrep, porb_out )

    end subroutine set_porb_default

    subroutine set_porb_by_spinmixed_basis( porb_out )
      real(kind=DP), intent(out) :: porb_out(num_orbitals, ndim_magmom )

      integer :: iorb, is1, is2, istmp, iopr, lmtt, iorb0
      complex(kind=CMPLDP) :: ztemp, z1, z2

      Do iorb = 1, num_orbitals
         iorb0 = int( (iorb-1)/ndim_spinor )+1
         call m_PP_tell_iorb_lmtt(iorb0,lmtt)

         ztemp = 0.0d0
         do iopr=1,nopr
            z1 = dcmplx( compr(i,iorb,iopr,ik ), &
                 &       compi(i,iorb,iopr,ik ) )
            z2 = dcmplx( compr(i,iorb,iopr,ik+1 ), &
                 &       compi(i,iorb,iopr,ik+1 ) )
            ztemp = ztemp + (z1+z2) *conjg(z1+z2) &
                 &  *( 1.d0+qorb(iorb0)/norm_phig_mpi(lmtt,iksnl) )
         end do
         porb_out(iorb,1) = ztemp /dble(nopr)
      End Do
    end subroutine set_porb_by_spinmixed_basis

    subroutine set_porb_rot_mode1( porb_out )
      real(kind=DP), intent(out) :: porb_out(num_orbitals, ndim_magmom )

      integer :: lmax, mmax, iopr
      integer :: ia, it, lmt, il, im, ip, iorb, lmtt, is1, is2, istmp
      integer :: lmt1, lmt2, im1, im2, iorb1, iorb2, lmtt1, lmtt2, im3
      integer :: tau, tau1, tau2
      real(kind=DP) :: c1, c2, ctmp
      complex(kind=CMPLDP) :: ztmp, z1, z2
      complex(kind=CMPLDP), allocatable :: dm_ssrep(:,:,:,:,:), dm(:,:,:,:,:)
      real(kind=DP), allocatable :: porb0(:,:,:,:)

      lmax = nloc;   mmax = 2*lmax -1

      allocate( dm_ssrep( ntau, lmax, mmax, mmax, ndim_chgpot ) );
      allocate( dm( ntau, lmax, mmax, mmax, ndim_magmom ) );
      allocate( porb0( ntau, lmax, mmax, ndim_magmom ) );

      do ia=1,natm
         it = ityp(ia)
         if(iproj_group(ia) == 0) cycle

         dm_ssrep = 0.0d0;  dm = 0.0d0;  porb0 = 0.d00

         ! diagonal part
         do lmt=1,ilmt_phi(it)
            il = ltp_phi(lmt,it); im  = mtp_phi(lmt,it);
            tau = taup_phi(lmt,it)
            ip  = iproj_phi(lmt,it)
            iorb = lmta_phi(lmt,ia); lmtt = lmtt_phi(lmt,it)

            Do is1=1, ndim_spinor
               Do is2=1, ndim_spinor
                  istmp = ( is1 -1 )*ndim_spinor + is2

                  ztmp = 0.0d0
                  do iopr=1,nopr
                     z1 = dcmplx( compr(i,iorb,iopr,ik+is1-1 ), &
                          &       compi(i,iorb,iopr,ik+is1-1 ) )
                     z2 = dcmplx( compr(i,iorb,iopr,ik+is2-1 ), &
                          &       compi(i,iorb,iopr,ik+is2-1 ) )
                     ztmp = ztmp + z1 *conjg(z2) &
                          &      *( 1.d0+qorb(iorb)/norm_phig_mpi(lmtt,iksnl) )
                  end do
                  dm_ssrep(tau,il,im,im,istmp) = dm_ssrep(tau,il,im,im,istmp) &
                       &                        +ztmp/dble(nopr)
               end do
            end Do
         end do
               ! non-diagonal part
         do lmt2=1,ilmt_phi(it)
            do lmt1=1,ilmt_phi(it)
               if ( lmt1 == lmt2 ) cycle
               if (ltp_phi(lmt1,it) /= ltp_phi(lmt2,it)) cycle
               if (taup_phi(lmt1,it) /= taup_phi(lmt2,it)) cycle

               il = ltp_phi(lmt1,it);     tau = taup_phi(lmt1,it)
               ip  = iproj_phi(lmt1,it)
               im1   = mtp_phi(lmt1,it);     im2 = mtp_phi(lmt2,it);
               iorb1 = lmta_phi(lmt1,ia);  iorb2 = lmta_phi(lmt2,ia)
               lmtt1 = lmtt_phi(lmt1,it);  lmtt2 = lmtt_phi(lmt2,it)

               c1 = ( qorb(iorb1) +qorb(iorb2) ) /2.0d0      ! approx
               c2 = sqrt(norm_phig_mpi(lmtt1,iksnl)*norm_phig_mpi(lmtt2,iksnl))

               Do is1=1, ndim_spinor
                  Do is2=1, ndim_spinor
                     istmp = ( is1 -1 )*ndim_spinor + is2
                     ztmp = 0.0d0
                     do iopr = 1,nopr
                        z1 = dcmplx( compr(i,iorb1,iopr,ik+is1-1 ), &
                             &       compi(i,iorb1,iopr,ik+is1-1 ) )
                        z2 = dcmplx( compr(i,iorb2,iopr,ik+is2-1 ), &
                             &       compi(i,iorb2,iopr,ik+is2-1 ) )
                        ztmp = ztmp + z1 *conjg(z2) &
                             &  *( 1.0d0 +c1 /c2 )
                     end do

                     dm_ssrep(tau,il,im1,im2,istmp) = dm_ssrep(tau,il,im1,im2,istmp) &
                          &                         + ztmp/dble(nopr)
                  end do
               end do
            end do
         end do

         dm(:,:,:,:,1) =  dm_ssrep(:,:,:,:,1) +dm_ssrep(:,:,:,:,4)
         dm(:,:,:,:,2) =  dm_ssrep(:,:,:,:,2) +dm_ssrep(:,:,:,:,3)
         dm(:,:,:,:,3) = (dm_ssrep(:,:,:,:,2) -dm_ssrep(:,:,:,:,3)) *zi
         dm(:,:,:,:,4) =  dm_ssrep(:,:,:,:,1) -dm_ssrep(:,:,:,:,4)

         ! ----
         Do tau=1, ntau
            Do il=1, lmax
               Do im1=1, mmax
                  DO is=1, ndim_magmom
                     ctmp = 0.0d0
                     Do im2=1, mmax
                        Do im3=1, mmax
                           ctmp = ctmp &
                                &    +porb_rot_matrix_real( ia,il,im2,im1 )  &
                                &       *dm(tau,il,im2,im3,is) &
                                &       *porb_rot_matrix_real( ia,il,im3,im1 )
                        End Do
                     End Do
                     porb0( tau, il, im1, is ) = ctmp
                  End DO
               End DO
            End Do
         End Do

         Do lmt=1,ilmt_phi(it)
            il = ltp_phi(lmt,it);  im  = mtp_phi(lmt,it);  iorb = lmta_phi(lmt,ia)
            tau = taup_phi(lmt,it)
            porb_out( iorb,1:ndim_magmom ) = porb0( tau, il,im,1:ndim_magmom )
         ENd do
      End do
      deallocate( dm ); deallocate( dm_ssrep ); deallocate( porb0 )

    end subroutine set_porb_rot_mode1

    subroutine set_porb_rot_mode2( porb_out )
      real(kind=DP), intent(out) :: porb_out(num_orbitals, ndim_magmom )

      integer :: lmax, mmax, iopr
      integer :: ia, it, lmt, il, im, ip, iorb, lmtt, is1, is2, istmp
      integer :: lmt1, lmt2, im1, im2, iorb1, iorb2, lmtt1, lmtt2, im3
      integer :: tau, tau1, tau2, is3, immax
      real(kind=DP) :: c1, c2
      complex(kind=CMPLDP) :: ztmp, z1, z2, z3
      complex(kind=CMPLDP), allocatable :: dm_ssrep(:,:,:,:,:)
      real(kind=DP), allocatable :: porb0(:,:,:,:)
      complex(kind=CMPLDP), allocatable :: zwk(:,:)
      complex(kind=CMPLDP) :: PauliMatrix( ndim_magmom, ndim_spinor, ndim_spinor )

      lmax = nloc;   mmax = 2*lmax -1

      allocate( dm_ssrep( ntau, lmax, mmax, mmax, ndim_chgpot ) );
      allocate( porb0( ntau, lmax, mmax*ndim_spinor, ndim_magmom ) );
      allocate( zwk( mmax*ndim_spinor, mmax*ndim_spinor ) );

      call m_ES_set_Pauli_Matrix( PauliMatrix )

      do ia=1,natm
         it = ityp(ia)
         if(iproj_group(ia) == 0) cycle

         dm_ssrep = 0.0d0;  porb0 = 0.d00

         ! diagonal part
         do lmt=1,ilmt_phi(it)
            il = ltp_phi(lmt,it); im  = mtp_phi(lmt,it);
            tau = taup_phi(lmt,it)
            ip  = iproj_phi(lmt,it); lmtt = lmtt_phi(lmt,it)
            iorb = lmta_phi(lmt,ia)

            Do is1=1, ndim_spinor
               Do is2=1, ndim_spinor
                  istmp = ( is1 -1 )*ndim_spinor + is2

                  ztmp = 0.0d0
                  do iopr=1,nopr
                     z1 = dcmplx( compr(i,iorb,iopr,ik+is1-1 ), &
                          &       compi(i,iorb,iopr,ik+is1-1 ) )
                     z2 = dcmplx( compr(i,iorb,iopr,ik+is2-1 ), &
                          &       compi(i,iorb,iopr,ik+is2-1 ) )
                     ztmp = ztmp + z1 *conjg(z2) &
                          &      *( 1.d0+qorb(iorb)/norm_phig_mpi(lmtt,iksnl) )
                  end do
                  dm_ssrep(tau,il,im,im,istmp) = dm_ssrep(tau,il,im,im,istmp) &
                       &                       +ztmp/dble(nopr)
               end do
            end Do
         end do
               ! non-diagonal part
         do lmt2=1,ilmt_phi(it)
            do lmt1=1,ilmt_phi(it)
               if ( lmt1 == lmt2 ) cycle
               if (ltp_phi(lmt1,it) /= ltp_phi(lmt2,it)) cycle
               if (taup_phi(lmt1,it) /= taup_phi(lmt2,it)) cycle

               il = ltp_phi(lmt1,it);     tau = taup_phi(lmt1,it)
               ip  = iproj_phi(lmt1,it)
               im1   = mtp_phi(lmt1,it);     im2 = mtp_phi(lmt2,it);
               iorb1 = lmta_phi(lmt1,ia);  iorb2 = lmta_phi(lmt2,ia)
               lmtt1 = lmtt_phi(lmt1,it);  lmtt2 = lmtt_phi(lmt2,it)

               c1 = ( qorb(iorb1) +qorb(iorb2) ) /2.0d0      ! approx
               c2 = sqrt(norm_phig_mpi(lmtt1,iksnl)*norm_phig_mpi(lmtt2,iksnl))

               Do is1=1, ndim_spinor
                  Do is2=1, ndim_spinor
                     istmp = ( is1 -1 )*ndim_spinor + is2
                     ztmp = 0.0d0
                     do iopr = 1,nopr
                        z1 = dcmplx( compr(i,iorb1,iopr,ik+is1-1 ), &
                             &       compi(i,iorb1,iopr,ik+is1-1 ) )
                        z2 = dcmplx( compr(i,iorb2,iopr,ik+is2-1 ), &
                             &       compi(i,iorb2,iopr,ik+is2-1 ) )
                        ztmp = ztmp + z1 *conjg(z2) &
                             &  *( 1.0d0 +c1 /c2 )
                     end do
                     dm_ssrep(tau,il,im1,im2,istmp) = dm_ssrep(tau,il,im1,im2,istmp) &
                          &                         + ztmp/dble(nopr)
                  end do
               end do
            end do
         end do

         Do tau=1, ntau
            Do il=1, lmax
               immax = 2 *il -1
               Do is=1, ndim_magmom
                  zwk = 0.0d0

                  Do im1=1, immax
                     Do im2=1, immax
                        Do is1=1, ndim_spinor
                           Do is2=1, ndim_spinor
                              ztmp = 0.0d0
                              Do is3=1, ndim_spinor
                                 istmp = ( is1 -1 )*ndim_spinor + is3
                                 ztmp = ztmp +dm_ssrep(tau,il,im1,im2,istmp) &
                                      &      *PauliMatrix(is,is3,is2)
                              End Do
                              zwk(im1+immax*(is1-1),im2+immax*(is2-1)) = ztmp
                           End Do
                        End Do
                     End Do
                  ENd Do
                  Do im1=1, immax *ndim_spinor
                     ztmp = 0.0d0
                     DO is2=1, ndim_spinor
                        DO is3=1, ndim_spinor
                           Do im2=1, immax
                              Do im3=1, immax
                                 z1 = porb_rot_matrix_cmplx(ia,il,immax*(is2-1)+im2,im1)
                                 z2 = zwk(immax*(is2-1)+im2,immax*(is3-1)+im3)
                                 z3 = porb_rot_matrix_cmplx(ia,il,immax*(is3-1)+im3,im1)
                                 ztmp = ztmp +conjg(z1) *z2 *z3
                              End Do
                           End Do
                        End Do
                     End DO
                     porb0( tau, il, im1, is ) = ztmp
                  End Do
               End Do
            End Do
         End Do
         Do lmt=1,ilmt_phi(it)
            il = ltp_phi(lmt,it);  im  = mtp_phi(lmt,it);  iorb = lmta_phi(lmt,ia)
            tau = taup_phi(lmt,it)
            Do is=1, ndim_spinor
               porb_out( (iorb-1)*ndim_spinor+is,1:ndim_magmom ) &
                    &    = porb0( tau,il,(im-1)*ndim_spinor+is,1:ndim_magmom )
            ENd do
         End Do
      End do
      deallocate(dm_ssrep); deallocate(porb0);    deallocate(zwk)

    end subroutine set_porb_rot_mode2

  end subroutine mkdos_with_GaussDistrib_noncl
! =================================================================== 11.0

! ====================== KT_add ======================= 13.0E
  subroutine make_dos_with_FDiracDistrib( kv, iwsc, kpt_weight )
    integer, intent(in) :: kv, iwsc
    real(kind=DP), intent(in) :: kpt_weight( kv )

    integer ::             i, ik, is, ie, id, ispin, iorb, iopr, lmtt, iksnl
    real(kind=DP) ::       Es, e, Ene1, c1, c2, w, DeltaE
    real(kind=DP) ::       porb

    DeltaE = DeltaE_dos
    Es = Eminimum - DeltaE*0.5d0

!    width_fermi_dirac = width

    if(ipridos >= 2) write(nfout,'(" !! Es, DeltaE = ",2d20.10)') Es, DeltaE
    dos = 0.d0; sumdos = 0.d0
    if(iwsc == TOTAL .and. sw_pdos == ON) then
       pdos=0.d0; sumpdos=0.d0
    end if
    do ispin = 1, nspin
       do ik = ispin, kv, nspin
          iksnl = (ik-1)/nspin + 1
          do i = 1, neg
             w = 1.d0
             if(iwsc >= 1 ) w = dos_weight(i,ik)
             e = eko(i,ik) - nETailWindows*DeltaE - Es
             is = e/DeltaE
             ie = is + 2*nETailWindows
             if(is < 0) then
                if(ipridos >= 2) write(nfout,'(" is = ",i5," < 0")') is
                is = 0
             end if
             if(ie >= nEWindows ) then
                if(ipridos >= 2) write(nfout,'(" ie = ",i5, " > nEWindows")')ie
                ie = nEWindows-1
             end if
             if(ipridos >= 2) then
                write(nfout,'(" !! eko(",i3,",",i3,") = ",f10.6," is, ie = ",2i7, 2f10.6)') &
                     & i,ik, eko(i,ik),is,ie, Es+is*DeltaE, Es+ie*DeltaE
             end if

             do id = is, ie
                ene1 = Es +id*DeltaE
                call width_fermi_dirac( ene1, eko(i,ik), smearing_width_fdirac, c1, c2 )
!!!                dos(id+1,ispin) = dos(id+1,ispin) + w *c1 *2.0d0 /dble(kv)
                dos(id+1,ispin) = dos(id+1,ispin) + w *c1 *2.0d0 *kpt_weight(ik)
             end do
             if(iwsc == TOTAL .and. sw_pdos == ON) then
                do iorb = 1,nlmta_phi
                   call m_PP_tell_iorb_lmtt(iorb,lmtt)
                   porb = 0.d0
!!$ASASASASAS
!!$                do iopr=1,nopr
!!$                   porb = porb + (compr(i,iorb,iopr,ik)**2 &
!!$                        &       + compi(i,iorb,iopr,ik)**2) &
!!$                        &     *(1.d0+qorb(iorb)/norm_phig_mpi(lmt,iksnl))
!!$                end do
                   if ( k_symmetry(ik) == GAMMA ) then
                      do iopr=1,nopr
                         porb = porb + compr(i,iorb,iopr,ik)**2  /2.0 &
                              &     *(1.d0+qorb(iorb)/(norm_phig_mpi(lmtt,iksnl)*2.) )
                      end do
                   else
                      do iopr=1,nopr
                         porb = porb + (compr(i,iorb,iopr,ik)**2 &
                              &       + compi(i,iorb,iopr,ik)**2) &
                              &     *(1.d0+qorb(iorb)/norm_phig_mpi(lmtt,iksnl))
                      end do
                   endif
!!$ASASASASAS
                   porb = porb/dble(nopr)
                   do id = is, ie
                      ene1 = Es +id*DeltaE
                      call width_fermi_dirac( ene1, eko(i,ik), &
                           &                  smearing_width_fdirac, c1, c2  )
                      pdos(id+1,iorb,ispin) = pdos(id+1,iorb,ispin) &
!!!                           & + porb *w *c1 *2.0d0 /dble(kv)
                           & + porb *w *c1 *2.0d0 *kpt_weight(ik)
                   end do
                end do
             end if
          end do
       end do
       sumdos(1,ispin) = dos(1,ispin)*DeltaE
       do id = 1, nEWindows-1
          sumdos(id+1,ispin) = sumdos(id,ispin) + dos(id+1,ispin)*DeltaE
       end do
       if(iwsc == TOTAL .and. sw_pdos == ON) then
          do iorb = 1,nlmta_phi
             sumpdos(1,iorb,ispin) = pdos(1,iorb,ispin)*DeltaE
             do id = 1, nEWindows-1
                sumpdos(id+1,iorb,ispin) = sumpdos(id,iorb,ispin) + pdos(id+1,iorb,ispin)*DeltaE
             end do
          end do
       end if
    end do
  end subroutine make_dos_with_FDiracDistrib
! =============================================================== 13.0E

  subroutine get_VBM(total_charge,deltad)
    real(kind=DP), intent(in) :: total_charge, deltad
    integer :: i
    logical :: found
    found = .false.
    i = nEWindows

    if(check_if_metalic_flag .and. .not.metalic_system) then
          ValenceBandMaximum = vbm
          if(ipridos >= 1) then
             if(dabs(ValenceBandMaximum) < 1.d6) then
                write(nfout,'(" ValenceBandMaximum = vbm = ",f12.4," total_charge = ",f12.4 &
                     & ," (m_ES_dos.get_VBM)")') ValenceBandMaximum,total_charge
             else
                write(nfout,'(" ValenceBandMaximum is not calculated")')
             end if
          end if
    else
! ============================== modified by K. Tagami ================ 11.0
!       if(nspin == 1) then
!          do while(.not.found.and.i >= 1)
!             i = i - 1
!             if(sumdos(i,nspin) < total_charge - deltad) found = .true.
!          end do
!       else if(nspin == 2) then
!          do while(.not.found.and.i<= nEWindows)
!             i = i - 1
!             if(sumdos(i,1)+sumdos(i,2) < total_charge - deltad) found = .true.
!          end do
!       end if

       if ( noncol ) then
          do while(.not.found.and.i >= 1)
             i = i - 1
             if(sumdos(i,1) < total_charge - deltad) found = .true.
          end do
       else
          if(nspin == 1) then
             do while(.not.found.and.i >= 1)
                i = i - 1
                if(sumdos(i,nspin) < total_charge - deltad) found = .true.
             end do
          else if(nspin == 2) then
             do while(.not.found.and.i<= nEWindows)
                i = i - 1
                if(sumdos(i,1)+sumdos(i,2) < total_charge - deltad) found = .true.
             end do
          end if
       endif
! ========================================================================= 11.0

       ValenceBandMaximum = Eminimum + i*DeltaE_dos
       if(ipridos >= 1) then
          write(nfout,'(" ValenceBandMaximum = ",f12.4, " i = ",i5," DeltaE_dos = ",d13.5 &
               & ," total_charge = ",f12.4," (m_ES_dos.get_VBM)")') &
               & ValenceBandMaximum,i,DeltaE_dos,total_charge
       end if
    end if

  end subroutine get_VBM

  subroutine write_dos(nf)
    integer, intent(in) :: nf
    integer       :: i, ii, nwdwidth, id, nrEWindows
    real(kind=DP) :: e,e_eV, dos_hr, dos_hr2, dos_eV, dos_eV2 &
         & , sumtotal, sumdos_avr, sumdos_avr2

    if(nwd_dos_window_width >= 1 .and. nwd_dos_window_width <= nEWindows) then
       nwdwidth = nwd_dos_window_width
    else
       nwdwidth = 1
    end if
    id = nwdwidth/2
    nrEWindows = nEWindows/nwdwidth * nwdwidth


    if(nspin == 1) then
       if(dos_write_format == WIDE) then
          write(nf,'("  No.   E(hr.)        dos(hr.)         E(eV)          dos(eV)              sum")')
       else if(dos_write_format == NARROW) then
          write(nf,'("  No.   E(hr.)    dos(hr.)       E(eV)          dos(eV)    sum")')
       else if(dos_write_format == eVunit) then
          write(nf,'("  No.   E(eV)          dos(eV)    sum")')
       end if
       do i = 1, nrEWindows, nwdwidth
          e = Eminimum + (i-1+id)*DeltaE_dos
          e_eV = (e - ValenceBandMaximum)*Hartree
          dos_hr = 0.d0; sumdos_avr = 0.d0
          do ii = i, i+nwdwidth-1
             dos_hr = dos_hr + dos(ii,nspin)
             sumdos_avr = sumdos_avr + sumdos(ii,nspin)
          end do
          dos_hr = dos_hr/nwdwidth; sumdos_avr = sumdos_avr/nwdwidth
!!$          dos_eV = dos(i,nspin)/Hartree
          dos_eV = dos_hr/Hartree
          if(dos_write_format == WIDE) then
             write(nf,'(i7,f10.5,f18.10,f14.6,f18.10,f20.10)') i+id, e, dos_hr &
                  & , e_eV, dos_eV, sumdos_avr
          else if(dos_write_format == NARROW) then
             write(nf,'(i7,f9.4,f16.8,f12.4,f16.10,f10.4)') i+id, e, dos_hr &
                  & , e_eV, dos_eV, sumdos_avr
          else if(dos_write_format == eVunit) then
             write(nf,'(i7,f12.4,f18.10,f10.4)') i+id, e_eV, dos_eV, sumdos_avr
          end if
       end do
    else if(nspin == 2) then
       if(dos_write_format == WIDE) then
          write(nf,'(1x,a47,5x,a45,6x,a27)') "No.  E(hr.)    dos_up(hr.)       dos_down(hr.)" &
               & ,"E(eV)         dos_up(eV)        dos_down(eV)" &
               & ,"sum_up   sum_down sum_total"
       else if(dos_write_format == NARROW) then
          write(nf,'(1x,a39,2x,a40,2x,a27)') "No.  E(hr.) dos_up(hr.)  dos_down(hr.)" &
               & ,"E(eV)       dos_up(eV)      dos_down(eV)" &
               & ,"sum_up   sum_down sum_total"
       else if(dos_write_format == eVunit) then
          write(nf,'(1x,a54,3x,a28)') "No.       E(eV)         dos_up(eV)        dos_down(eV)" &
               & ,"sum_up   sum_down  sum_total"
       end if
       do i = 1, nrEWindows, nwdwidth
          e = Eminimum + (i-1+id)*DeltaE_dos
          e_eV = (e - ValenceBandMaximum)*Hartree
          dos_hr = 0.d0; dos_hr2 = 0.d0;  sumdos_avr = 0.d0; sumdos_avr2 = 0.d0
          do ii = i, i+nwdwidth-1
             dos_hr  = dos_hr + dos(ii,1)
             dos_hr2 = dos_hr2 + dos(ii,2)
             sumdos_avr  = sumdos_avr  + sumdos(ii,1)
             sumdos_avr2 = sumdos_avr2 + sumdos(ii,2)
          end do
          dos_hr = dos_hr/nwdwidth; dos_hr2 = dos_hr2/nwdwidth
          sumdos_avr = sumdos_avr/nwdwidth; sumdos_avr2 = sumdos_avr2/nwdwidth
!!$          dos_eV  = dos(i,1)/Hartree
!!$          dos_eV2 = dos(i,2)/Hartree
!!$          sumtotal = sumdos(i,1) + sumdos(i,2)
          dos_eV  = dos_hr/Hartree
          dos_eV2 = dos_hr2/Hartree
          sumtotal = sumdos_avr + sumdos_avr2
          if(dos_write_format== WIDE) then
             write(nf,'(i7,f9.4,2f18.10,f15.4,2f18.10,3f10.4)') &
                  & i, e, dos_hr,dos_hr2,e_eV, dos_eV, dos_eV2, sumdos_avr,sumdos_avr2,sumtotal
          else if(dos_write_format == NARROW) then
             write(nf,'(i7,f9.4,2f14.6, f12.4,2f14.8,3f10.4)') &
                  & i, e, dos_hr,dos_hr2,e_eV, dos_eV, dos_eV2, sumdos_avr,sumdos_avr2,sumtotal
          else if(dos_write_format == eVunit) then
             write(nf,'(i7,f12.4,2f18.10,3f10.4)') i, e_eV, dos_eV, dos_eV2, sumdos_avr,sumdos_avr2,sumtotal
          end if
       end do
    end if
    write(nf,'("END")')
  end subroutine write_dos

! ================================== added by K. Tagami ================= 11.0
  subroutine write_dos_noncl(nf)
    integer, intent(in) :: nf
    integer       :: i, ii, nwdwidth, id, nrEWindows

    real(kind=DP) :: dos_hr( ndim_magmom ), sumdos_avr( ndim_magmom )
    real(kind=DP) :: dos_eV( ndim_magmom )
    real(kind=DP) :: e, e_eV
    real(kind=DP) :: sumtotal

    if(nwd_dos_window_width >= 1 .and. nwd_dos_window_width <= nEWindows) then
       nwdwidth = nwd_dos_window_width
    else
       nwdwidth = 1
    end if
    id = nwdwidth/2

    nrEWindows = nEWindows/nwdwidth * nwdwidth

    write(nf,'(2x,A,5x,A,4x,A)') "No.      E(eV)",&
         &     "  dos_chg(eV)    dos_mx(eV)    dos_my(eV)    dos_mz(eV)", &
         &      "sum_chg     sum_mx      sum_my      sum_mz"

    do i = 1, nrEWindows, nwdwidth
       e = Eminimum + (i-1+id)*DeltaE_dos
       e_eV = (e - ValenceBandMaximum)*Hartree

       dos_hr = 0.d0; sumdos_avr = 0.d0
       do ii = i, i+nwdwidth-1
          dos_hr(:) = dos_hr(:) + dos(ii,:)
          sumdos_avr(:) = sumdos_avr(:) + sumdos(ii,:)
       end do

       dos_hr = dos_hr/nwdwidth; sumdos_avr = sumdos_avr/nwdwidth
       dos_eV = dos_hr/Hartree

       write(nf,'(i5,f14.8,4f14.8,4f12.6)') &
               & i, e_eV, dos_eV(1:ndim_magmom), sumdos_avr(1:ndim_magmom)
!       write(nf,*) &
!               & i, e_eV, dos_eV(1:ndim_magmom), sumdos_avr(1:ndim_magmom)
    end do

    write(nf,'("END")')
  end subroutine write_dos_noncl

  subroutine write_dos_noncl_for_totchg(nf)

    integer, intent(in) :: nf
    integer       :: i, ii, nwdwidth, id, nrEWindows
    real(kind=DP) :: e,e_eV, dos_hr, dos_hr2, dos_eV, dos_eV2 &
         & , sumtotal, sumdos_avr, sumdos_avr2

    if(nwd_dos_window_width >= 1 .and. nwd_dos_window_width <= nEWindows) then
       nwdwidth = nwd_dos_window_width
    else
       nwdwidth = 1
    end if
    id = nwdwidth/2

    nrEWindows = nEWindows/nwdwidth * nwdwidth

    write(nf,'("     No.   E(hr.)        dos(hr.)         E(eV)          dos(eV)              sum")')

    do i = 1, nrEWindows, nwdwidth
       e = Eminimum + (i-1+id)*DeltaE_dos
       e_eV = (e - ValenceBandMaximum)*Hartree
       dos_hr = 0.d0; sumdos_avr = 0.d0

       do ii = i, i+nwdwidth-1
          dos_hr = dos_hr + dos(ii,1)
          sumdos_avr = sumdos_avr + sumdos(ii,1)
       end do

       dos_hr = dos_hr/nwdwidth; sumdos_avr = sumdos_avr/nwdwidth
       dos_eV = dos_hr/Hartree

       write(nf,'(i7,f10.5,f18.10,f14.6,f18.10,f20.10)') i+id, e, dos_hr &
            & , e_eV, dos_eV, sumdos_avr
    end do

    write(nf,'("END")')
  end subroutine write_dos_noncl_for_totchg
! ======================================================================== 11.0

  subroutine write_pdos(nf)
    integer, intent(in) :: nf

    integer :: iorb

    do iorb=1,nlmta_phi
       call write_pdos_orbital(nf,iorb)
!!$ASASASASAS
!!$       write(nf,'("END")')
!!$ASASASASAS
    end do

  end subroutine write_pdos

! =================================== added by K. Tagami =============== 11.0
  subroutine write_pdos_noncl(nf)
    integer, intent(in) :: nf
    integer :: iorb, kfac, num_orbitals

    kfac = 1
    if ( sw_diagonalize_population == ON ) then
       if ( population_diag_mode /= DIAG_CHARGE_DENSITY_MATRIX &
            & .and. population_diag_mode /= LOCAL_POINT_GROUP ) then
          kfac = ndim_spinor
       endif
    endif
    num_orbitals = nlmta_phi *kfac

    do iorb=1,num_orbitals
       call write_pdos_orbital_noncl(nf,iorb)
    end do

  end subroutine write_pdos_noncl
! ========================================================================= 11.0

  subroutine write_pdos_orbital(nf,iorb)
    integer, intent(in) :: nf
    integer, intent(in) :: iorb
    integer :: ia,il,im,tau,nspher
    integer       :: i, ii, nwdwidth, id, nrEWindows
    real(kind=DP) :: e,e_eV, dos_hr, dos_hr2, dos_eV, dos_eV2 &
         & , sumtotal, sumdos_avr, sumdos_avr2

    if(nwd_dos_window_width >= 1 .and. nwd_dos_window_width <= nEWindows) then
       nwdwidth = nwd_dos_window_width
    else
       nwdwidth = 1
    end if
    id = nwdwidth/2
    nrEWindows = nEWindows/nwdwidth * nwdwidth

    call m_PP_tell_iorb_ia_l_m_tau(iorb,ia,il,im,tau)
!!$ASASASASAS
    if ( iproj_group(ia)== 0 ) return
!!$ASASASASAS
    if ( if_pdos(ia) == 0 ) return
    if ( sw_diagonalize_population == ON ) then
       write(nf,'("PDOS: ia= ",i0," l=",i3," m''=",i3," t=",i3)') ia,il-1,im,tau
    else
       write(nf,'("PDOS: ia= ",i0," l=",i3," m=",i3," t=",i3)') ia,il-1,im,tau
    endif

    if(nspin == 1) then
       if(dos_write_format == WIDE .or. dos_write_format == NARROW) then
          write(nf,'("  No.   E(hr.)        dos(hr.)         E(eV)          dos(eV)              sum")')
       else if(dos_write_format == eVunit) then
          write(nf,'("  No.   E(eV)          dos(eV)              sum")')
       end if
       do i = 1, nrEWindows, nwdwidth
          e = Eminimum + (i-1+id)*DeltaE_dos
          e_eV = (e - ValenceBandMaximum)*Hartree
          dos_hr = 0.d0; sumdos_avr = 0.d0
          do ii = i, i+nwdwidth-1
             dos_hr = dos_hr + pdos(ii,iorb,nspin)
             sumdos_avr = sumdos_avr + sumpdos(ii,iorb,nspin)
          end do
          dos_hr = dos_hr/nwdwidth; sumdos_avr = sumdos_avr/nwdwidth
          dos_eV = dos_hr/Hartree
          if(dos_write_format == WIDE) then
             write(nf,'(i7,f10.5,f18.10,f14.6,f18.10,f20.10)') i+id, e, dos_hr &
                  & , e_eV, dos_eV, sumdos_avr
          else if(dos_write_format == NARROW) then
             write(nf,'(i7,f9.4,f16.8,f12.4,f16.10,f10.4)') i+id, e, dos_hr, e_eV, dos_eV, sumdos_avr
          else if(dos_write_format == eVunit) then
             write(nf,'(i7,f14.6,f18.10,f10.4)') i+id,e_eV,dos_eV,sumdos_avr
          end if
       end do
    else if(nspin == 2) then
       if(dos_write_format == WIDE .or. dos_write_format == NARROW) then
          write(nf,'(1x,a47,5x,a45,6x,a27)') "No.  E(hr.)    dos_up(hr.)       dos_down(hr.)" &
               & ,"E(eV)         dos_up(eV)        dos_down(eV)" &
               & ,"sum_up   sum_down sum_total"
       else if(dos_write_format == eVunit) then
          write(nf,'(1x,a54,3x,       a28)') "No.       E(eV)         dos_up(eV)        dos_down(eV)" &
               & ,"sum_up   sum_down  sum_total"
      end if
       do i = 1, nrEWindows, nwdwidth
          e = Eminimum + (i-1+id)*DeltaE_dos
          e_eV = (e - ValenceBandMaximum)*Hartree
          dos_hr = 0.d0; dos_hr2 = 0.d0;  sumdos_avr = 0.d0; sumdos_avr2 = 0.d0
          do ii = i, i+nwdwidth-1
             dos_hr  = dos_hr + pdos(ii,iorb,1)
             dos_hr2 = dos_hr2 + pdos(ii,iorb,2)
             sumdos_avr  = sumdos_avr  + sumpdos(ii,iorb,1)
             sumdos_avr2 = sumdos_avr2 + sumpdos(ii,iorb,2)
          end do
          dos_hr = dos_hr/nwdwidth; dos_hr2 = dos_hr2/nwdwidth
          sumdos_avr = sumdos_avr/nwdwidth; sumdos_avr2 = sumdos_avr2/nwdwidth
          dos_eV  = dos_hr/Hartree
          dos_eV2 = dos_hr2/Hartree
          sumtotal = sumdos_avr + sumdos_avr2
          if(dos_write_format == WIDE) then
             write(nf,'(i7,f9.4,2f18.10,f15.4,2f18.10,3f10.4)') &
                  & i, e, dos_hr,dos_hr2,e_eV, dos_eV, dos_eV2, sumdos_avr,sumdos_avr2,sumtotal
          else if(dos_write_format == NARROW) then
             write(nf,'(i7,f9.4,2f14.6,f12.4,2f14.8,3f10.4)') &
                  & i, e, dos_hr,dos_hr2,e_eV, dos_eV, dos_eV2, sumdos_avr,sumdos_avr2,sumtotal
          else if(dos_write_format == eVunit) then
             write(nf,'(i7,f12.4,2f18.10,3f10.4)') i, e_eV, dos_eV, dos_eV2, sumdos_avr,sumdos_avr2,sumtotal
          end if
       end do
    end if
!!$ASASASASAS
    write(nf,'("END")')
!!$ASASASASAS
  end subroutine write_pdos_orbital

! ================================= added by K. Tagami ==================== 11.0
  subroutine write_pdos_orbital_noncl(nf,iorb)
    integer, intent(in) :: nf
    integer, intent(in) :: iorb
    integer :: ia,il,im,tau,nspher
    integer       :: i, ii, nwdwidth, id, nrEWindows

    real(kind=DP) :: dos_hr( ndim_magmom ), sumdos_avr( ndim_magmom )
    real(kind=DP) :: dos_eV( ndim_magmom )
    real(kind=DP) :: e, e_eV

    integer :: kfac, iorb0, im0, my_l
    real(kind=DP) :: val_j, val_mj

    if(nwd_dos_window_width >= 1 .and. nwd_dos_window_width <= nEWindows) then
       nwdwidth = nwd_dos_window_width
    else
       nwdwidth = 1
    end if
    id = nwdwidth/2

    nrEWindows = nEWindows/nwdwidth * nwdwidth

    kfac = 1
    if ( sw_diagonalize_population == ON ) then
       if ( population_diag_mode /= DIAG_CHARGE_DENSITY_MATRIX &
            & .and. population_diag_mode /= LOCAL_POINT_GROUP ) then
          kfac = ndim_spinor
       endif
    endif

    iorb0 = int( (iorb-1)/kfac )+1
    call m_PP_tell_iorb_ia_l_m_tau(iorb0,ia,il,im0,tau)
    im = ( im0 -1 )*kfac +mod( iorb-1, kfac ) +1

    if ( iproj_group(ia)== 0 ) return

    if ( if_pdos(ia) == 0 ) return

    if ( sw_diagonalize_population == ON ) then
       if ( population_diag_mode == DIAG_LS ) then
          my_l = il -1
          if ( my_l > 0 ) then
             if ( im <= 2*my_l ) then
                val_j = my_l -0.5d0;    val_mj = -val_j +(im -1)
             else
                val_j = my_l +0.5d0;    val_mj = -val_j +(im -2*my_l -1)
             endif
          else
             val_j = my_l +0.5d0;       val_mj = -val_j +(im -2*my_l -1)
          endif
       endif
    endif
! ---
    if ( sw_diagonalize_population == ON ) then
       if ( population_diag_mode == DIAG_CHARGE_DENSITY_MATRIX &
            &  .or. population_diag_mode == LOCAL_POINT_GROUP ) then
          write(nf,'("PDOS: ia= ",i0," l=",i3," m''=",i3," t=",i3)') ia,il-1,im,tau
       else if ( population_diag_mode == DIAG_LS ) then
          write(nf,'(A,I0,A,I3,A,I3,A,I3,A,F4.1,A,F4.1,A)') &
               &   "PDOS: ia= ", ia, " l=", il -1," ms''=",im," t=",tau, &
               &   "                 ( j= ", val_j, "  mj= ", val_mj, " )"
       else
          write(nf,'("PDOS: ia= ",i0," l=",i3," ms''=",i3," t=",i3)') ia,il-1,im,tau
       endif
    else
       write(nf,'("PDOS: ia= ",i0," l=",i3," m=",i3," t=",i3)') ia,il-1,im,tau
    endif

    if ( sw_diagonalize_population == ON .and. use_rotated_compri == YES &
         &   .and. population_diag_mode /= DIAG_CHARGE_DENSITY_MATRIX ) then
       write(nf,'(2x,A,5x,A,4x,A)') "No.      E(eV)",&
            &     "  dos(eV)", &
            &      "sum"
    else
       write(nf,'(2x,A,5x,A,4x,A)') "No.      E(eV)",&
            &     "  dos_chg(eV)    dos_mx(eV)    dos_my(eV)    dos_mz(eV)", &
            &      "sum_chg     sum_mx      sum_my      sum_mz"
    endif

    do i = 1, nrEWindows, nwdwidth
       e = Eminimum + (i-1+id)*DeltaE_dos
       e_eV = (e - ValenceBandMaximum)*Hartree

       dos_hr = 0.0d0; sumdos_avr = 0.0d0
       do ii = i, i+nwdwidth-1
          dos_hr(:)  = dos_hr(:) + pdos(ii,iorb,:)
          sumdos_avr(:)  = sumdos_avr(:)  + sumpdos(ii,iorb,:)
       end do

       dos_hr = dos_hr/nwdwidth;         sumdos_avr = sumdos_avr/nwdwidth;
       dos_eV  = dos_hr/Hartree;

       if ( sw_diagonalize_population == ON .and. use_rotated_compri == YES &
            &   .and. population_diag_mode /= DIAG_CHARGE_DENSITY_MATRIX ) then
          write(nf,'(i5,f14.8,1f14.8,1f12.6)') &
               & i, e_eV, dos_eV(1), sumdos_avr(1)
       else
          write(nf,'(i5,f14.8,4f14.8,4f12.6)') &
               & i, e_eV, dos_eV(1:ndim_magmom), sumdos_avr(1:ndim_magmom)
       endif
    end do

    write(nf,'("END")')
  end subroutine write_pdos_orbital_noncl
! ========================================================================= 11.0

  subroutine m_ESdos_gaussdistrib(nfdos,icomponent)
    integer, intent(in) ::      nfdos, icomponent

    call alloc_eko_and_substitution(kv3) ! eko_l -> eko
    call find_Erange(eko,neg,kv3)
    call alloc_dos(1,icomponent)
    call make_dos_with_GaussianDistrib( kv3, icomponent, qwgt )
    if(icomponent == TOTAL) call get_VBM(totch,DeltaDVBM)
    if(mype == 0) call write_dos(nfdos)
!!$    call write_dos(nfout)
    if(mype == 0 .and. icomponent == TOTAL .and. sw_pdos == ON) then
       call write_pdos(nfdos)
    end if
! ---
    if ( sw_calc_pcohp == ON ) call m_ES_calc_pcohp
! ---
    call dealloc_eko()
    call dealloc_dos()
  end subroutine m_ESdos_gaussdistrib

! ================================ added by K. Tagami =================== 11.0
  subroutine m_ESdos_gaussdistrib_noncl(nfdos,icomponent)
    integer, intent(in) ::      nfdos, icomponent

    call alloc_eko_and_substit_noncl(kv3) ! eko_l -> eko
    call find_Erange( eko, neg, kv3/ndim_spinor )
    call alloc_dos(1,icomponent)

    call mkdos_with_GaussDistrib_noncl( kv3, icomponent, qwgt )
    if(icomponent == TOTAL) call get_VBM(totch,DeltaDVBM)

    if(mype == 0) then
       if ( calc_dos_magmom_contrib == YES ) then
          call write_dos_noncl(nfdos)
       else
          call write_dos_noncl_for_totchg(nfdos)
       endif
    endif

    if(mype == 0 .and. icomponent == TOTAL .and. sw_pdos == ON) then
       call write_pdos_noncl(nfdos)
    end if

    call dealloc_eko()
    call dealloc_dos()
  end subroutine m_ESdos_gaussdistrib_noncl
! ======================================================================= 11.0

  subroutine m_ESdos_gaussdistrib_ek(nfdos,icomponent)
    integer, intent(in) ::      nfdos, icomponent

    call alloc_eko_and_substitution_ek(kv3_ek) ! eko_ek -> eko
    call find_Erange(eko,neg,kv3_ek)
    call alloc_dos(1,icomponent)
    if(icomponent == TOTAL .and. sw_pdos == ON) then
       call pdos_ek()
       call m_ES_dos_gather_compr_compi_ek()
    endif
    call make_dos_with_GaussianDistrib( kv3_ek, icomponent, qwgt_ek )
    if(icomponent == TOTAL) call get_VBM(totch,DeltaDVBM)
    if(mype == 0) call write_dos(nfdos)
    if(mype == 0 .and. icomponent == TOTAL .and. sw_pdos == ON) then
       call write_pdos(nfdos)
    end if
    call dealloc_eko()
    call dealloc_dos()
  end subroutine m_ESdos_gaussdistrib_ek

! ==================================== added by K. Tagami ================ 11.0
  subroutine m_ESdos_gaussdistrib_ek_noncl(nfdos,icomponent)
    integer, intent(in) ::      nfdos, icomponent

    call alloc_eko_and_substit_ek_noncl(kv3_ek) ! eko_ek -> eko
    call find_Erange(eko,neg,kv3_ek/ndim_spinor)
    call alloc_dos(1,icomponent)
    if(icomponent == TOTAL .and. sw_pdos == ON) then
       call pdos_ek()
       call m_ES_dos_gather_compr_compi_ek()
    endif

    call mkdos_with_GaussDistrib_noncl( kv3_ek, icomponent, qwgt_ek )
    if(icomponent == TOTAL) call get_VBM(totch,DeltaDVBM)

    if(mype == 0) then
       if ( calc_dos_magmom_contrib == YES ) then
          call write_dos_noncl(nfdos)
       else
          call write_dos_noncl_for_totchg(nfdos)
       endif
    endif

    call dealloc_eko()
    call dealloc_dos()
  end subroutine m_ESdos_gaussdistrib_ek_noncl
! ========================================================================= 11.0

! ==================- KT_add ================= 13.0E
  subroutine m_ESdos_FdiracDistrib(nfdos,icomponent)
    integer, intent(in) ::      nfdos, icomponent

    call alloc_eko_and_substitution(kv3) ! eko_l -> eko
    call find_Erange_fermidirac(eko,neg,kv3)

    call alloc_dos(1,icomponent)
    call make_dos_with_FDiracDistrib( kv3, icomponent, qwgt )

!!!    if(icomponent == TOTAL) call get_VBM(totch,DeltaDVBM)
    if(icomponent == TOTAL) then
       ValenceBandMaximum = efermi
       write(nfout,*) '!!!! Efermi is used as ValenceBandMaximum in the Fermi Dirac case'
    endif

    if(mype == 0) call write_dos(nfdos)
!!$    call write_dos(nfout)
    if(mype == 0 .and. icomponent == TOTAL .and. sw_pdos == ON) then
       call write_pdos(nfdos)
    end if
    if ( sw_calc_pcohp == ON ) call m_ES_calc_pcohp

    call dealloc_eko(); call dealloc_dos()

  end subroutine m_ESdos_FdiracDistrib

  subroutine m_ESdos_FDiracDistrib_ek(nfdos,icomponent)
    integer, intent(in) ::      nfdos, icomponent

    call alloc_eko_and_substitution_ek(kv3_ek) ! eko_ek -> eko
    call find_Erange_fermidirac(eko,neg,kv3_ek)

    call alloc_dos(1,icomponent)
    if(icomponent == TOTAL .and. sw_pdos == ON) then
       call pdos_ek()
       call m_ES_dos_gather_compr_compi_ek()
    endif
    call make_dos_with_FDiracDistrib( kv3_ek, icomponent, qwgt_ek )

    if(icomponent == TOTAL) then
       ValenceBandMaximum = efermi
       write(nfout,*) '!!!! Efermi is used as ValenceBandMaximum in the Fermi Dirac case'
    endif

    if(mype == 0) call write_dos(nfdos)
    if(mype == 0 .and. icomponent == TOTAL .and. sw_pdos == ON) then
       call write_pdos(nfdos)
    end if
    call dealloc_eko();     call dealloc_dos()
  end subroutine m_ESdos_FDiracDistrib_ek
! ============================================ 13.0E

  subroutine m_ESdos_write_dos_header(nfdos,aldos_or_layerdos,icomponent,tagwords)
    integer, intent(in):: nfdos, aldos_or_layerdos, icomponent
    character(len=40), intent(in) :: tagwords
    if(mype == 0) then
       if(aldos_or_layerdos == ALDOS) then
          write(nfdos,'("ALDOS     num_atom = ",i7)') icomponent
       else if(aldos_or_layerdos == LAYERDOS) then
          write(nfdos,'("LAYERDOS   num_layer = ",i7,2x,a40)') icomponent,tagwords
       end if
    end if
  end subroutine m_ESdos_write_dos_header

  ! ----------------------------------------------

  subroutine m_ESdos_tetrahedral(nfdos,icomponent,mode)
    integer, intent(in) ::  nfdos, mode,icomponent
!!$#ifndef NO_TETRAHEDRON
    real(kind=DP), parameter :: delta = 1.d-12
    integer, parameter       :: idim = 3
    integer        :: neig,ispin,ip2,ik, ip, ib, jb, ie, nelm
    real(kind=DP)  :: et,wei,clpm
    real(kind=DP), allocatable, dimension(:,:,:)   :: eeig2
    real(kind=DP), allocatable,dimension(:,:)  :: dos_weight2 ! d(neg,np2)
!    real(kind=DP), allocatable, dimension(:,:,:,:) :: compr, compi
!    real(kind=DP), allocatable, dimension(:,:) :: norm_phig_mpi
    real(kind=DP), pointer, dimension(:) :: eawk,cdwk,cswk,e
    real(kind=DP), allocatable, dimension(:,:) ::  e_mpi
    real(kind=DP), pointer, dimension(:,:,:) :: cdos,cind
    real(kind=DP), allocatable, dimension(:,:) :: doswk,dosinwk
    integer, allocatable, dimension(:,:) :: nttra ! d(mtetra,4)
    integer ::                                  id_sname = -1
    integer :: i,iorb,iopr,ikorb,lmtt, ipridos_t, nEwindows_plus,mtetra,ikt,iloop, j, iss, ies
    integer, parameter :: ncl = 8
    real(kind=DP) :: porb

    call tstatc0_begin('m_ESdos_tetrahedral ', id_sname)

    if(mode /= SCF .and. mode /= EK) then
       write(nfout,'(" !dos:  mode = ",i6," <<m_ESdos_tetrahedral>>")') mode
       return
    end if

    if(npes > 1) then
       if(mype == 0) ipridos_t = ipridos
       call mpi_bcast(ipridos_t, 1, mpi_integer, 0, MPI_CommGroup, ierr)
    else
       ipridos_t = ipridos
    end if

    if(ipridos_t >= 3) then
       if(printable) write(nfout,'(" !dos: np2, np0 = ",2i7," <<m_ESdos_tetrahedral>>")') np2, np0
       allocate(e_mpi(neg,kv3))
       if(printable) write(nfout,'(" !dos: icomponent = ",i9," <<m_ESdos_tetrahedral>>")') icomponent
       if(printable) write(nfout,'(" !dos: --energy_eigenvalue--")')
       e_mpi = 0.d0
       do ik = 1, kv3, af+1
          if(map_k(ik) /= myrank_k) cycle
          do ie = 1, neg
             ip = neordr(ie,ik)
             if(map_e(ip) == myrank_e) e_mpi(ie,ik) = eko_l(map_z(ip),ik)
          end do
       end do
       if(npes >=2) then
          call mpi_allreduce(MPI_IN_PLACE,e_mpi,neg*kv3,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
       end if
       if(printable) then
          do ik = 1, kv3
             write(nfout,'(" !dos:  ik = ",i7)') ik
             write(nfout,'(" !dos: ",8f10.6)') (e_mpi(ie,ik),ie=1,neg)
          end do
       end if
       deallocate(e_mpi)
    end if

    allocate(eeig2(np2,neg,nspin)); eeig2 = 0.d0
    if(icomponent == TOTAL .and. sw_pdos == ON .and. ekmode==OFF) then
       if(allocated(compr))     deallocate(compr)
       if(allocated(compi))     deallocate(compi)
       if(allocated(norm_phig_mpi)) deallocate(norm_phig_mpi)
       allocate(compr(neg,nlmta_phi,nopr,np2*nspin));  compr = 0.d0
       allocate(compi(neg,nlmta_phi,nopr,np2*nspin));  compi = 0.d0
       allocate(norm_phig_mpi(nlmtt_phi,np2));  norm_phig_mpi = 0.d0
    end if

    if(icomponent == TOTAL .and. sw_pdos == ON.and. ekmode==ON) then
       call pdos_ek()
       call m_ES_dos_gather_compr_compi_ek()
    endif

    if(mode==EK .and. ipridos_t >= 1) then
       write(nfout,'(" !dos eko_ek ")')
       do ispin=1,nspin
          do ip2=1,np2
             ik = nspin*(ip2-1)+ispin
             write(nfout,'(" !dos -- ik = ",i5)') ik
             write(nfout,'(" !dos ",8f8.4)')(eko_ek(ib,ik),ib=1,neg)
          end do
       end do
    end if
    neig=neg
    iss = ista_spin
    ies = iend_spin
    if(mode == EK) then
      iss = 1
      ies = nspin
    endif
    do ispin=iss, ies
       do ip2=1,np2
          ik=nspin*(ip2-1)+ispin
          if(mode == EK) then
             do ib=1,neg
                eeig2(ip2,ib,ispin)=eko_ek(ib,ik)
             enddo
             do ib = 1,neg-1
                do jb = ib+1, neg
                   if(eeig2(ip2,jb,ispin) < eeig2(ip2,ib,ispin)-delta) then
                      et = eeig2(ip2,ib,ispin)
                      eeig2(ip2,ib,ispin) = eeig2(ip2,jb,ispin)
                      eeig2(ip2,jb,ispin) = et
                   end if
                end do
             end do
          else if(mode == SCF) then
             if(map_k(ik) /= myrank_k) cycle
             do ib = 1, neg
                ip = neordr(ib,ik)
                if(map_e(ip) == myrank_e) then
                   eeig2(ip2,ib,ispin)=eko_l(map_z(ip),ik)
                   if(icomponent == TOTAL .and. sw_pdos == ON) then
                      compr(ib,1:nlmta_phi,1:nopr,ik)=compr_l(map_z(ip),1:nlmta_phi,1:nopr,ik)
!!$ASASASASAS
!!$                      compi(ib,1:nlmta_phi,1:nopr,ik)=compi_l(map_z(ip),1:nlmta_phi,1:nopr,ik)
                      if ( k_symmetry(ik) /= GAMMA ) then
                         compi(ib,1:nlmta_phi,1:nopr,ik)=compi_l(map_z(ip),1:nlmta_phi,1:nopr,ik)
                      else
!                         compi = 0.0d0
                         compi(ib,1:nlmta_phi,1:nopr,ik) = 0.0d0
                      endif
!!$ASASASASAS
                      norm_phig_mpi(1:nlmtt_phi,ip2)=norm_phig(1:nlmtt_phi,ip2)
                   end if
                end if
             enddo
          end if
       enddo
    end do

    if(mode == SCF) then
       if(npes >= 2) then
          nelm = np2*neg*nspin
          call mpi_allreduce(MPI_IN_PLACE,eeig2,nelm,mpi_double_precision,mpi_sum, mpi_kg_world,ierr)
          call mpi_allreduce(MPI_IN_PLACE,eeig2,nelm,mpi_double_precision,mpi_sum, mpi_ge_world,ierr)
          !deallocate(eeig2_mpi)
          if(icomponent == TOTAL .and. sw_pdos == ON) then
             nelm = neg*nlmta_phi*nopr*np2*nspin
             call mpi_allreduce(MPI_IN_PLACE,compr,nelm,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
             call mpi_allreduce(MPI_IN_PLACE,compr,nelm,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)
             call mpi_allreduce(MPI_IN_PLACE,compi,nelm,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
             call mpi_allreduce(MPI_IN_PLACE,compi,nelm,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)
             nelm = nlmtt_phi*np2
             call mpi_allreduce(MPI_IN_PLACE,norm_phig_mpi,nelm,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)
             norm_phig_mpi = norm_phig_mpi / dble(nrank_s)
          end if
       end if
    end if

    if(ipridos >= 2) then
!!$       if(mype == 0) then
          write(nfout,'(" !dos: --eeig2(energy_eigenvalue)--")')
          iloop = (neg-1)/ncl+1
          do ispin = 1, nspin
             write(nfout,'(" !dos: ispin = ",i5)') ispin
             do ip2 = 1, np2
                do i = 1, iloop
                   write(nfout,'(" !dos: (",i5,") ",8f9.5)') &
                        & ip2,(eeig2(ip2,ib,ispin),ib=ncl*(i-1)+1,min(neg,ncl*i))
                end do
!!$                write(nfout,'(" !dos:  ip2 = ",i7)') ip2
!!$                write(nfout,'(" !dos: ",8f9.5)') (eeig2(ip2,ib,ispin),ib=1,neg)
             end do
          end do
          if(icomponent == TOTAL .and. sw_pdos == ON) then
          do ispin=1,nspin
             write(nfout,'(" !dos: ispin = ",i5)') ispin
             do ip2 = 1, np2
                ikorb=nspin*(ik-1)+ispin
                write(nfout,'(" !dos:  ip2 = ",i7)') ip2
                write(nfout,'(" !dos norm_phig: ",10f8.4)') &
                          & norm_phig_mpi(1:nlmtt_phi,ip2)
                do iopr=1,nopr
                   do iorb=1,nlmta_phi
                      write(nfout,'(" !dos: ik=",i5," iopr=",i5," iorb=",i5)') ikorb,iopr,iorb
                      write(nfout,'(" !dos compr: ",10f8.4)') compr(1:neg,iorb,iopr,ikorb)
                      write(nfout,'(" !dos compi: ",10f8.4)') compi(1:neg,iorb,iopr,ikorb)
                   end do
                end do
             end do
          end do
          end if
!!$       end if
    end if

    if(mode == EK) then
       call find_Erange(eko_ek,neg,kv3_ek)
    else
       Eminimum = minval(eeig2) - 0.005 ! (hartree)
       Emaximum = maxval(eeig2) + 0.005 ! (hartree)
       nEWindows = (Emaximum - Eminimum)/DeltaE_dos + 1
    end if

#ifdef DEBUG_LDOS
    if(icomponent/=TOTAL) then
       write(6,'(" <<m_ESdos_tetrahedral>")')
       if (ekmode == ON) then
          ikt = kv3_ek
       else
          ikt = kv3
       end if
       write(6,'(" kv3 = ",i8, " neg = ",i8)') ikt, neg
       do j = 1, ikt
          write(6,'(" ik = ", i8)') j
          write(6,'(10f8.4)') (dos_weight(i,j),i=1,neg)
       end do
    end if
#endif
    if(ipridos >= 2) &
         & write(nfout,'(" !dos: Emaximum, Eminimum, nEwindows = ",2f10.6,i7)') &
         &       Eminimum,Emaximum, nEWindows

    call alloc_dos(0,icomponent)

    if(nspin == 1) then
       wei = 2.d0
    else
       wei = 1.d0
    end if

    if(ipridos >= 2) write(nfout,*) ' === tetrahedron method', &
                &  ' for k-space integration === <<m_ESdos_tetrahedral>>'
    if(ipridos >= 2) write(nfout,'(" !m_ES_dos dos_subroutine = ",i5)') dos_subroutine

    do ispin=1,nspin
!!$       write(nfout,*) ' ispin=',ispin
!!$       write(nfout,'(" !! nxyz_tetra = ",3i7)') nxyz_tetra(1:3)
!!$       write(nfout,'(" !! np0 = ",i7)') np0
!!$       do ik = 1, np0
!!$          write(nfout,'(" !! ik, ip20 = ",2i6," <<m_ESdos_tetrahedral>>")') ik,ip20(ik)
!!$       end do
       if(dos_subroutine == 3) then
          allocate(e(0:nEwindows))
          e(0:nEwindows) = (/(Eminimum + DeltaE_dos*ie,ie=0,nEWindows)/)
          ! e(0) = Eminimum, e(1) = Eminimum+DeltaE_dos, e(2) = Eminimum+DeltaE_dos*2,...

          allocate(cdos(np2,neg,0:nEWindows)); cdos = 0.d0
          allocate(cind(np2,neg,0:nEWindows)); cind = 0.d0
          allocate(eawk(np0)); eawk = 0.d0
          allocate(cdwk(np0)); cdwk = 0.d0
          allocate(cswk(np0)); cswk = 0.d0
          call nstt3i(idim,nEwindows,e,nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3) &
               &  ,np2,np2,neig,eeig2(1,1,ispin) &
               &  ,ip20,np0,eawk,cdwk,cswk,np2,neig,cdos,cind,width_tetra ) ! -> cdos, cind
          deallocate(cswk,cdwk,eawk)
          clpm = 1.d0
          do ik=1,np2
             ikt = nspin*(ik-1)+ispin
             if(icomponent == TOTAL .and. sw_pdos == ON) ikorb=nspin*(ik-1)+ispin
             do ib = 1, neg
!!$             if(icomponent >= 1) clpm = dos_weight(ib,ik)
                if(icomponent /= TOTAL) clpm = dos_weight(ib,ikt)

                do ie = 0, nEwindows
                   dos(ie,ispin) = dos(ie,ispin) + cdos(ik,ib,ie)*clpm*wei
                   sumdos(ie,ispin) = sumdos(ie,ispin) + cind(ik,ib,ie)*clpm*wei
                end do
                if(icomponent == TOTAL .and. sw_pdos == ON) then
                   do iorb = 1,nlmta_phi
                      call m_PP_tell_iorb_lmtt(iorb,lmtt)
                      porb = 0.d0
!!$ASASASASAS
!!$                      do iopr=1,nopr
!!$                         porb = porb + (   compr(ib,iorb,iopr,ikorb)**2 &
!!$                              &          + compi(ib,iorb,iopr,ikorb)**2) &
!!$                              &        *(1.d0+qorb(iorb)/norm_phig_mpi(lmt,ik))
!!$                      end do
                      if ( k_symmetry(ik) == GAMMA) then
                         do iopr=1,nopr
                            porb = porb + compr(ib,iorb,iopr,ikorb)**2 /2.0 &
                            &  *(1.d0+qorb(iorb)/(norm_phig_mpi(lmtt,ik)*2.) )
                         end do
                      else
                         do iopr=1,nopr
                            porb = porb + ( compr(ib,iorb,iopr,ikorb)**2 &
                                 &        + compi(ib,iorb,iopr,ikorb)**2) &
                            &    *(1.d0+qorb(iorb)/norm_phig_mpi(lmtt,ik))
                         enddo
                      endif
!!$ASASASASAS
                      porb = porb/dble(nopr)
                      ! debug
                      ! print *,'debug iorb=',iorb,' porb=',porb
                      ! end debug
                      do ie = 0, nEwindows
                         pdos(ie,iorb,ispin) = pdos(ie,iorb,ispin) + cdos(ik,ib,ie)*wei*porb
                         sumpdos(ie,iorb,ispin) = sumpdos(ie,iorb,ispin) + cind(ik,ib,ie)*wei*porb
                      end do
                   end do
                end if
             end do
          end do
          deallocate(cind,cdos)
          deallocate(e)
       else if(dos_subroutine == 4) then
          allocate(cdos(0:nEwindows,np2,neg)); cdos = 0.d0
          allocate(cind(0:nEwindows,np2,neg)); cind = 0.d0
          nEwindows_plus = nEwindows
          if(mod(nEwindows_plus,2)==1) nEwindows_plus = nEwindows_plus+1
          allocate(doswk(0:nEwindows_plus,4))
          allocate(dosinwk(0:nEwindows_plus,4))
          mtetra = product(nxyz_tetra(1:3))*6
          allocate(nttra(mtetra,4))
          call prepare_nttra(nxyz_tetra,mtetra,nttra)
          write(nfout,'(" !dos after prepare_nttra")')
          write(nfout,'(" !dos neig = ",i5)') neig
          cdos = 0.d0; cind =0.d0
          call nstt4i(ipridos,idim,nEwindows,nxyz_tetra,np2,np2,neig,eeig2(1,1,ispin) &
               &   ,  ip20,np0,np2,neig,mtetra,nttra,deltae_dos &
               &   ,  nEwindows_plus,doswk,dosinwk,cdos,cind)
          deallocate(nttra)
          deallocate(dosinwk,doswk)
          clpm = 1.d0
          do ik=1,np2
             ikt=nspin*(ik-1)+ispin
             do ib = 1, neg
                if(icomponent /= TOTAL) clpm = dos_weight(ib,ikt)
                ! clpm -> clpm(ib,ik)
                do ie = 0, nEwindows
                   dos(ie,ispin) = dos(ie,ispin) + cdos(ie,ik,ib)*clpm*wei
                   sumdos(ie,ispin) = sumdos(ie,ispin) + cind(ie,ik,ib)*clpm*wei
                end do
             end do
          end do
          if(icomponent == TOTAL .and. sw_pdos == ON) then
             do ik = 1,np2
                ikorb=nspin*(ik-1)+ispin
                do ib = 1, neg
                   do iorb = 1,nlmta_phi
                      call m_PP_tell_iorb_lmtt(iorb,lmtt)
                      porb = 0.d0
!!$ASASASASAS
!!$                      do iopr=1,nopr
!!$                         porb = porb + (   compr(ib,iorb,iopr,ikorb)**2 &
!!$                              &          + compi(ib,iorb,iopr,ikorb)**2) &
!!$                              &        *(1.d0+qorb(iorb)/norm_phig_mpi(lmt,ik))
!!$                      end do
                      if ( k_symmetry(ik) == GAMMA ) then
                         do iopr=1,nopr
                            porb = porb + compr(ib,iorb,iopr,ikorb)**2 /2.0 &
                              &   *(1.d0+qorb(iorb)/(norm_phig_mpi(lmtt,ik)*2.))
                           end do
                      else
                         do iopr=1,nopr
                            porb = porb + ( compr(ib,iorb,iopr,ikorb)**2 &
                                 &    + compi(ib,iorb,iopr,ikorb)**2) &
                                 &   *(1.d0+qorb(iorb)/norm_phig_mpi(lmtt,ik))
                         end do
                      endif
!!$ASASASASAS
                      porb = porb/dble(nopr)
                      ! porb -> porb(ib,iorb,ik)
                      ! debug
                      ! print *,'debug iorb=',iorb,' porb=',porb
                      ! end debug
                      do ie = 0, nEwindows
                         pdos(ie,iorb,ispin) = pdos(ie,iorb,ispin) + cdos(ie,ik,ib)*wei*porb
                         sumpdos(ie,iorb,ispin) = sumpdos(ie,iorb,ispin) + cind(ie,ik,ib)*wei*porb
                      end do
                   end do
                end do
             end do
          end if
          deallocate(cind,cdos)
       else if(dos_subroutine == 5) then
          nEwindows_plus = nEwindows
          if(mod(nEwindows_plus,2)==1) nEwindows_plus = nEwindows_plus+1
          allocate(doswk(0:nEwindows_plus,4))
          allocate(dosinwk(0:nEwindows_plus,4))
          mtetra = product(nxyz_tetra(1:3))*6
          allocate(nttra(mtetra,4))
          call prepare_nttra(nxyz_tetra,mtetra,nttra)
          allocate(dos_weight2(neg,np2))
          if(icomponent == TOTAL) then
             dos_weight2 = 1.d0
          else
             do ik = 1, np2
                ikt = nspin*(ik-1)+ispin
                do ib = 1, neg
                   dos_weight2(ib,ik) = dos_weight(ib,ikt)
                end do
             end do
          end if
          call care_of_degenerate_state(neig,np2,nspin,ispin,eeig2,dos_weight2) !dos_weight
          dos_weight2 = dos_weight2*wei
          dos(:,ispin) = 0.d0
          sumdos(:,ispin) = 0.d0
          call nstt5i(ipridos,idim,Eminimum,Emaximum,nEwindows,nxyz_tetra,np2,np2 &
               &   ,  neig,eeig2(1,1,ispin) &
               &   ,  ip20,np0,np2,neig,mtetra,nttra,deltae_dos,dos_weight2 &
               &   ,  nEwindows_plus,doswk,dosinwk,dos(0,ispin),sumdos(0,ispin))
!!$          call nstt5i_3D(ipridos,idim,Eminimum,Emaximum,nEwindows,nxyz_tetra,np2,np2 &
!!$               &   ,  neig,eeig2(1,1,ispin) &
!!$               &   ,  ip20,np0,np2,neig,mtetra,nttra,deltae_dos,dos_weight2 &
!!$               &   ,  nEwindows_plus,doswk,dosinwk,dos(0,ispin),sumdos(0,ispin))

          if(icomponent == TOTAL .and. sw_pdos == ON) then
             do iorb = 1, nlmta_phi
                call m_PP_tell_iorb_lmtt(iorb,lmtt)
                do ik = 1, np2
                   ikorb=nspin*(ik-1)+ispin
                   do ib = 1, neg
                      porb = 0.d0
!!$ASASASASAS
!!$                      do iopr=1,nopr
!!$                         porb = porb + (   compr(ib,iorb,iopr,ikorb)**2 &
!!$                              &          + compi(ib,iorb,iopr,ikorb)**2) &
!!$                              &        *(1.d0+qorb(iorb)/norm_phig_mpi(lmt,ik))
!!$                      end do
                      if ( k_symmetry(ik) == GAMMA ) then
                         do iopr=1,nopr
                            porb = porb + compr(ib,iorb,iopr,ikorb)**2 /2.0 &
                             &   *(1.d0+qorb(iorb)/(norm_phig_mpi(lmtt,ik)*2.) )
                         end do
                      else
                         do iopr=1,nopr
                            porb = porb + ( compr(ib,iorb,iopr,ikorb)**2 &
                                 &      + compi(ib,iorb,iopr,ikorb)**2) &
                                 &    *(1.d0+qorb(iorb)/norm_phig_mpi(lmtt,ik))
                         end do
                      endif
!!$ASASASASAS
                      dos_weight2(ib,ik) = porb/dble(nopr)*wei
                   end do
                end do
                call care_of_degenerate_state(neig,np2,nspin,ispin,eeig2,dos_weight2) !dos_weight
                call nstt5i(ipridos,idim,Eminimum,Emaximum,nEwindows,nxyz_tetra,np2 &
                     &   ,  np2,neig,eeig2(1,1,ispin) &
                     &   ,  ip20,np0,np2,neig,mtetra,nttra,deltae_dos,dos_weight2 &
                     &   ,  nEwindows_plus,doswk,dosinwk,pdos(0,iorb,ispin),sumpdos(0,iorb,ispin))
             end do
          end if
          deallocate(dos_weight2)
          deallocate(nttra)
          deallocate(dosinwk,doswk)
       !np2: #kpoints independent
       !np0: #all kpoints
!!$       call doscal_nstt(nEwindows,nxyz_tetra,np2,neig,eeig2(1,1,ispin),ip20,np0,wei,cdos,cind)

       end if
    end do

    deallocate(eeig2)

    if(icomponent == TOTAL) call get_VBM(totch,1.d-12)  ! -> ValenceBandMaximum
    if(mype == 0) call write_dos(nfdos)
    if(mype == 0 .and. icomponent == TOTAL .and. sw_pdos == ON) call write_pdos(nfdos)

    call dealloc_dos()
    call tstatc0_end(id_sname)
!!$#endif
  end subroutine m_ESdos_tetrahedral


! ============================= added by K. Tagami ======================= 11.0
  subroutine m_ESdos_tetrahedral_noncl(nfdos,icomponent,mode)

    integer, intent(in) ::  nfdos, mode,icomponent
!!$#ifndef NO_TETRAHEDRON
    real(kind=DP), parameter :: delta = 1.d-12
    integer, parameter       :: idim = 3
    integer        :: neig, ip2, ik, ip, ib, jb, ie
    real(kind=DP)  :: et,wei,clpm
    real(kind=DP), pointer, dimension(:,:,:)   :: eeig2, eeig2_mpi
    real(kind=DP), allocatable,dimension(:,:)  :: dos_weight2 ! d(neg,np2)

    real(kind=DP), pointer, dimension(:,:,:,:) :: compr, compr_mpi
    real(kind=DP), pointer, dimension(:,:,:,:) :: compi, compi_mpi
    real(kind=DP), pointer, dimension(:,:) :: norm_phig_mpi, norm_phig_mpi2
    real(kind=DP), pointer, dimension(:) :: eawk,cdwk,cswk,e
    real(kind=DP), allocatable, dimension(:,:) ::  e_mpi

    real(kind=DP), pointer, dimension(:,:,:) :: cdos,cind
    real(kind=DP), allocatable, dimension(:,:) :: doswk,dosinwk
    integer, allocatable, dimension(:,:) :: nttra ! d(mtetra,4)
    integer ::                                  id_sname = -1
    integer :: i, iorb,iopr, lmtt, ipridos_t, nEwindows_plus, mtetra, iloop
    integer, parameter :: ncl = 8

! -------
    complex(kind=CMPLDP), allocatable :: porb_ssrep(:,:)
    complex(kind=CMPLDP), allocatable :: porb_ssrep5(:,:,:)
    real(kind=DP), allocatable :: porb(:,:)
    real(kind=DP), allocatable :: dos_weight3(:,:,:)
!
    integer :: is1, is2, istmp, ismax
    complex(kind=CMPLDP) :: z1, z2, ztemp
! ---

    call tstatc0_begin('m_ESdos_tetrahedral_noncl ', id_sname)

    if(mode /= SCF .and. mode /= EK) then
       write(nfout,'(" !dos:  mode = ",i6," <<m_ESdos_tetrahedral_noncl>>")') mode
       return
    end if

    if(npes > 1) then
       if(mype == 0) ipridos_t = ipridos
       call mpi_bcast(ipridos_t, 1, mpi_integer, 0, MPI_CommGroup, ierr)
    else
       ipridos_t = ipridos
    end if

! --------------------------------------------------------
    if (ipridos_t >= 3) then
       if (printable) then
          write(nfout,'(" !dos: np2, np0 = ",2i7," <<m_ESdos_tetrahedral_noncl>>")') &
               &      np2, np0
       endif

       allocate(e_mpi(neg,kv3))

       if (printable) then
          write(nfout,'(" !dos: icomponent = ",i9," <<m_ESdos_tetrahedral_noncl>>")') &
               &      icomponent
       endif

       if(printable) write(nfout,'(" !dos: --energy_eigenvalue--")')
       e_mpi = 0.d0

       do ik = 1, kv3, ndim_spinor
          if (map_k(ik) /= myrank_k) cycle
          do ie = 1, neg
             ip = neordr(ie,ik)
             if (map_e(ip) == myrank_e) e_mpi(ie,ik) = eko_l(map_z(ip),ik)
          end do
       end do

       if (npes >=2) then
          call mpi_allreduce(MPI_IN_PLACE,e_mpi,neg*kv3,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
       end if
       if (printable) then
          do ik = 1, kv3, ndim_spinor
             write(nfout,'(" !dos:  ik = ",i7)') ik
             write(nfout,'(" !dos: ",8f10.6)') (e_mpi(ie,ik),ie=1,neg)
          end do
       end if
       deallocate(e_mpi)
    end if
! --------------------------------------------------

    allocate( eeig2(np2,neg,1) ); eeig2 = 0.d0
    if (icomponent == TOTAL .and. sw_pdos == ON) then
       allocate( compr(neg,nlmta_phi,nopr,np2*ndim_spinor));  compr = 0.d0
       allocate( compi(neg,nlmta_phi,nopr,np2*ndim_spinor));  compi = 0.d0
       allocate( norm_phig_mpi(nlmtt_phi,np2));  norm_phig_mpi = 0.d0
    end if

    if (mode==EK .and. ipridos_t >= 1) then
       write(nfout,'(" !dos eko_ek ")')
       do ip2=1,np2
          ik = ndim_spinor*( ip2 -1 ) +1
          write(nfout,'(" !dos -- ik = ",i5)') ik
          write(nfout,'(" !dos ",8f8.4)')( eko_ek(ib,ik),ib=1,neg)
       end do
    end if

    neig=neg

! -------------------------------------------------------
    Do ip2=1, np2
       ik = ndim_spinor *(ip2-1) +1

       if(mode == EK) then
          do ib=1,neg
             eeig2(ip2,ib,1)=eko_ek(ib,ik)
          enddo
          do ib = 1,neg-1
             do jb = ib+1, neg
                if ( eeig2(ip2,jb,1) < eeig2(ip2,ib,1)-delta) then
                   et = eeig2(ip2,ib,1)
                   eeig2(ip2,ib,1) = eeig2(ip2,jb,1)
                   eeig2(ip2,jb,1) = et
                end if
             end do
          end do
       else if(mode == SCF) then
          if(map_k(ik) /= myrank_k) cycle

          do ib = 1, neg
             ip = neordr(ib,ik)
             if (map_e(ip) == myrank_e) then
                eeig2(ip2,ib,1) = eko_l(map_z(ip),ik)
                if (icomponent == TOTAL .and. sw_pdos == ON) then

                   Do is1=1, ndim_spinor
                      compr(ib,1:nlmta_phi,1:nopr,ik+is1-1) &
                           &  = compr_l(map_z(ip),1:nlmta_phi,1:nopr,ik+is1-1)

                      if ( k_symmetry(ik) /= GAMMA ) then
                         compi(ib,1:nlmta_phi,1:nopr,ik+is1-1) &
                              &  = compi_l(map_z(ip),1:nlmta_phi,1:nopr,ik+is1-1)
                      else
                         call phase_error_with_msg(nfout,"Not supported : Gamma symmetry in noncollinear",__LINE__,__FILE__)
                      endif
                   End do
                   norm_phig_mpi(1:nlmtt_phi,ip2) = norm_phig(1:nlmtt_phi,ip2)
                end if
             end if
          enddo
       end if

    End do
! -------------------------------------------------------

    if (mode == SCF) then
       if (npes >= 2) then
          allocate(eeig2_mpi(np2,neg,1)); eeig2_mpi = 0.d0
          call mpi_allreduce( eeig2, eeig2_mpi, np2*neg, mpi_double_precision, &
               &              mpi_sum, mpi_kg_world, ierr )
          eeig2 = eeig2_mpi
          deallocate(eeig2_mpi)

          if (icomponent == TOTAL .and. sw_pdos == ON) then
            allocate( compr_mpi(neg,nlmta_phi,nopr,np2*ndim_spinor) ); compr_mpi=0.d0
            allocate( compi_mpi(neg,nlmta_phi,nopr,np2*ndim_spinor) ); compi_mpi=0.d0
            allocate( norm_phig_mpi2(nlmtt_phi,np2) ); norm_phig_mpi2 = 0.d0
!!$            call mpi_allreduce( compr, compr_mpi, neg*nlmta_phi*nopr*np2*ndim_spinor,&
!!$                 &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
            call mpi_allreduce( compr, compr_mpi, neg*nlmta_phi*nopr*np2*ndim_spinor,&
                 &              mpi_double_precision, mpi_sum, mpi_kg_world, ierr )
!!$            call mpi_allreduce( compi,compi_mpi, neg*nlmta_phi*nopr*np2*ndim_spinor, &
!!$                 &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
            call mpi_allreduce( compi,compi_mpi, neg*nlmta_phi*nopr*np2*ndim_spinor, &
                 &              mpi_double_precision, mpi_sum, mpi_kg_world, ierr )
!!$            call mpi_allreduce( norm_phig_mpi, norm_phig_mpi2, nlmtt_phi*np2, &
!!$                 &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
            call mpi_allreduce( norm_phig_mpi, norm_phig_mpi2, nlmtt_phi*np2, &
                 &              mpi_double_precision, mpi_sum, mpi_kg_world, ierr )

            compr = compr_mpi;    compi = compi_mpi
            norm_phig_mpi = norm_phig_mpi2 / dble(nrank_s)
            deallocate(compr_mpi,compi_mpi,norm_phig_mpi2)
          end if

       end if
    end if

! ------------------------------------------------------------
    if (ipridos >= 2) then

       write(nfout,'(" !dos: --eeig2(energy_eigenvalue)--")')
       iloop = (neg-1)/ncl+1

       do ip2 = 1, np2
          do i = 1, iloop
             write(nfout,'(" !dos: (",i5,") ",8f9.5)') &
                  & ip2,(eeig2(ip2,ib,1),ib=ncl*(i-1)+1,min(neg,ncl*i))
          end do
       end do

       if (icomponent == TOTAL .and. sw_pdos == ON) then

          do ip2 = 1, np2
             ik = ndim_spinor *(ip2-1) + 1

             write(nfout,'(" !dos:  ip2 = ",i7)') ip2
             write(nfout,'(" !dos norm_phig: ",10f8.4)') &
                  & norm_phig_mpi(1:nlmtt_phi,ip2)

             do iopr=1,nopr
                do iorb=1,nlmta_phi
                   write(nfout,'(" !dos: ik=",i5," iopr=",i5," iorb=",i5)')&
                        &  ik, iopr, iorb
                   write(nfout,'(" !dos compr: ",10f8.4)') &
                        & compr(1:neg,iorb,iopr,ik)
                   write(nfout,'(" !dos compi: ",10f8.4)') &
                        & compi(1:neg,iorb,iopr,ik)
                end do
             end do
          end do
       end if

    end if
! -----------------------------------------------------

    if (mode == EK) then

       allocate( e_mpi( neg,kv3_ek/ndim_spinor) );   e_mpi = 0.d0;

       do ip2=1, np2
          ik = ( ip2 -1 )*ndim_spinor + 1
          e_mpi( :, ip2 ) = eko_ek( :, ik )
       end do
       call find_Erange( e_mpi, neg, kv3_ek /ndim_spinor )
       deallocate( e_mpi )

    else
       Eminimum = minval(eeig2) - 0.005 ! (hartree)
       Emaximum = maxval(eeig2) + 0.005 ! (hartree)
       nEWindows = (Emaximum - Eminimum)/DeltaE_dos + 1
    end if

    if (ipridos >= 2) &
         & write(nfout,'(" !dos: Emaximum, Eminimum, nEwindows = ",2f10.6,i7)') &
         &       Eminimum,Emaximum, nEWindows

    call alloc_dos(0,icomponent)

! ----------------------------------------------------
    wei = 1.d0                ! ??????????????? 2.0 ? docchi _

    if(ipridos >= 2) write(nfout,*) ' === tetrahedron method', &
                &  ' for k-space integration === <<m_ESdos_tetrahedral_noncl>>'
    if(ipridos >= 2) write(nfout,'(" !m_ES_dos dos_subroutine = ",i5)') dos_subroutine

! ---------------------------------------------------

    if ( calc_dos_magmom_contrib == YES ) then
       ismax = ndim_magmom
    else
       ismax = 1
    endif
! -------------------------------------------------

    if (dos_subroutine == 3) then
       allocate(e(0:nEwindows))
       e(0:nEwindows) = (/(Eminimum + DeltaE_dos*ie,ie=0,nEWindows)/)
                      ! e(0) = Eminimum, e(1) = Eminimum+DeltaE_dos,
                      ! e(2) = Eminimum+DeltaE_dos*2,...

       allocate(cdos(np2,neg,0:nEWindows)); cdos = 0.d0
       allocate(cind(np2,neg,0:nEWindows)); cind = 0.d0
       allocate(eawk(np0)); eawk = 0.d0
       allocate(cdwk(np0)); cdwk = 0.d0
       allocate(cswk(np0)); cswk = 0.d0

       call nstt3i(idim,nEwindows,e,nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3) &
            &  ,np2,np2,neig,eeig2(1,1,1) &
            &  ,ip20,np0,eawk,cdwk,cswk,np2,neig,cdos,cind,width_tetra ) ! -> cdos, cind
       deallocate(cswk,cdwk,eawk)

       clpm = 1.d0

       if (icomponent == TOTAL .and. sw_pdos == ON) then
          allocate( porb_ssrep( nlmta_phi, ndim_chgpot ) ); porb_ssrep = 0.0d0
          allocate( porb( nlmta_phi, ndim_magmom ) ); porb = 0.0d0
       endif

       do ip2 = 1, np2
          ik = ndim_spinor *(ip2-1) +1

          do ib = 1, neg

             Do istmp=1, ismax
                clpm = dos_weight_noncl( ib,ik,istmp )

                do ie = 0, nEwindows
                   dos(ie,istmp) = dos(ie,istmp) + cdos(ip2,ib,ie)*clpm*wei
                   sumdos(ie,istmp) = sumdos(ie,istmp) + cind(ip2,ib,ie)*clpm*wei
                end do
             End do

             if (icomponent == TOTAL .and. sw_pdos == ON) then
                clpm = 1.0d0
                porb_ssrep = 0.d0

! -----------------
                Do iorb = 1,nlmta_phi
                   call m_PP_tell_iorb_lmtt(iorb,lmtt)

                   if ( k_symmetry(ik) == GAMMA ) then
                      call phase_error_with_msg(nfout,'Not supported ',__LINE__,__FILE__)
                   else

                      Do is1=1, ndim_spinor
                         Do is2=1, ndim_spinor
                            istmp = ( is1 -1 )*ndim_spinor + is2

                            ztemp = 0.0d0
                            do iopr=1,nopr
                               z1 = dcmplx( compr(ib,iorb,iopr,ik+is1-1 ), &
                                    &       compi(ib,iorb,iopr,ik+is1-1 ) )
                               z2 = dcmplx( compr(ib,iorb,iopr,ik+is2-1 ), &
                                    &       compi(ib,iorb,iopr,ik+is2-1 ) )
                               ztemp = ztemp + z1 *conjg(z2) &
                                    &      *( 1.d0+qorb(iorb)/norm_phig_mpi(lmtt,ip2) )
                            end do

                            porb_ssrep(iorb,istmp) = porb_ssrep(iorb,istmp) &
                                 &                   + ztemp /dble(nopr)
                         End do
                      End do
                   endif
                End Do
! ---------------------
                call m_ES_DensMat_To_MagMom_porb( nlmta_phi, porb_ssrep, porb )

                do ie = 0, nEwindows
                   pdos(ie,iorb,:) = pdos(ie,iorb,:) &
                        &           + cdos(ip2,ib,ie) *wei *porb(iorb,:)
                   sumpdos(ie,iorb,:) = sumpdos(ie,iorb,:) &
                        &           + cind(ip2,ib,ie) *wei *porb(iorb,:)
                end do

             end if
          end do
       end do

       deallocate(cind,cdos);     deallocate(e)
       if ( allocated( porb ) ) deallocate( porb )
       if ( allocated( porb_ssrep ) ) deallocate( porb_ssrep )

! ---------------------------------------------------
    else if(dos_subroutine == 4) then
       allocate(cdos(0:nEwindows,np2,neg)); cdos = 0.d0
       allocate(cind(0:nEwindows,np2,neg)); cind = 0.d0
       nEwindows_plus = nEwindows

       if(mod(nEwindows_plus,2)==1) nEwindows_plus = nEwindows_plus+1
       allocate(doswk(0:nEwindows_plus,4))
       allocate(dosinwk(0:nEwindows_plus,4))
       mtetra = product(nxyz_tetra(1:3))*6

       allocate(nttra(mtetra,4))
       call prepare_nttra(nxyz_tetra,mtetra,nttra)
       write(nfout,'(" !dos after prepare_nttra")')
       write(nfout,'(" !dos neig = ",i5)') neig

       cdos = 0.d0; cind =0.d0
       call nstt4i(ipridos,idim,nEwindows,nxyz_tetra,np2,np2,neig,eeig2(1,1,1) &
            &   ,  ip20,np0,np2,neig,mtetra,nttra,deltae_dos &
            &   ,  nEwindows_plus,doswk,dosinwk,cdos,cind)
       deallocate(nttra);     deallocate(dosinwk,doswk)

       clpm = 1.d0

       do ip2 = 1,np2
          ik = ndim_spinor*(ip2-1) +1

          do ib = 1, neg
             Do istmp=1, ismax
                clpm = dos_weight_noncl(ib,ik,istmp)
                                       ! clpm -> clpm(ib,ik)
                do ie = 0, nEwindows
                   dos(ie,istmp) = dos(ie,istmp) + cdos(ie,ip2,ib)*clpm*wei
                   sumdos(ie,istmp) = sumdos(ie,istmp) + cind(ie,ip2,ib)*clpm*wei
                end do
             End do
          end do
       end do

       if (icomponent == TOTAL .and. sw_pdos == ON) then
          allocate( porb_ssrep( nlmta_phi, ndim_chgpot ) ); porb_ssrep = 0.0d0
          allocate( porb( nlmta_phi, ndim_magmom ) ); porb = 0.0d0

          do ip2 = 1,np2
             ik = ndim_spinor *(ip2-1) +1

             do ib = 1, neg
                porb_ssrep = 0.d0

                do iorb = 1,nlmta_phi
                   call m_PP_tell_iorb_lmtt(iorb,lmtt)

                   if ( k_symmetry(ik) == GAMMA ) then
                      call phase_error_with_msg(nfout,'Not supported ',__LINE__,__FILE__)
                   else
                      Do is1=1, ndim_spinor
                         Do is2=1, ndim_spinor
                            istmp = ( is1 -1 )*ndim_spinor + is2

                            ztemp = 0.0d0
                            do iopr=1,nopr
                               z1 = dcmplx( compr(ib,iorb,iopr,ik+is1-1 ), &
                                    &       compi(ib,iorb,iopr,ik+is1-1 ) )
                               z2 = dcmplx( compr(ib,iorb,iopr,ik+is2-1 ), &
                                    &       compi(ib,iorb,iopr,ik+is2-1 ) )
                               ztemp = ztemp + z1 *conjg(z2) &
                                    &      *( 1.d0+qorb(iorb)/norm_phig_mpi(lmtt,ip2) )
                            end do

                            porb_ssrep(iorb,istmp) = porb_ssrep(iorb,istmp) &
                                 &                   + ztemp /dble(nopr)
                         End do
                      End do
                   end if
                End do
! -----
                call m_ES_DensMat_To_MagMom_porb( nlmta_phi, porb_ssrep, porb )
!
                do ie = 0, nEwindows
                   pdos(ie,iorb,:) = pdos(ie,iorb,:) &
                        &           + cdos(ie,ip2,ib) *wei *porb(iorb,:)
                   sumpdos(ie,iorb,:) = sumpdos(ie,iorb,:) &
                        &           + cind(ie,ip2,ib) *wei *porb(iorb,:)
                end do
             end do
          end do

          deallocate( porb ); deallocate( porb_ssrep )
       end if

       deallocate(cind,cdos)

! ---------------------------------------------------------------
    else if(dos_subroutine == 5) then
       nEwindows_plus = nEwindows
       if(mod(nEwindows_plus,2)==1) nEwindows_plus = nEwindows_plus+1
       allocate(doswk(0:nEwindows_plus,4))
       allocate(dosinwk(0:nEwindows_plus,4))
       mtetra = product(nxyz_tetra(1:3))*6

       allocate(nttra(mtetra,4))
       call prepare_nttra(nxyz_tetra,mtetra,nttra)
       allocate(dos_weight2(neg,np2))

       Do istmp=1, ismax

          do ip2 = 1, np2
             ik = ndim_spinor *(ip2-1) +1
             do ib = 1, neg
                dos_weight2(ib,ip2) = dos_weight_noncl(ib,ik,istmp)
             end do
          end do

          call care_of_degenerate_state( neig, np2, 1, 1, eeig2, dos_weight2 )
                                                       !dos_weight

          dos_weight2 = dos_weight2*wei

          dos(:,istmp) = 0.d0;           sumdos(:,istmp) = 0.d0

          call nstt5i(ipridos,idim,Eminimum,Emaximum,nEwindows,nxyz_tetra,np2,np2 &
               &   ,  neig,eeig2(1,1,1) &
               &   ,  ip20,np0,np2,neig,mtetra,nttra,deltae_dos,dos_weight2 &
               &   ,  nEwindows_plus,doswk,dosinwk,dos(0,istmp),sumdos(0,istmp) )
       End do

       if (icomponent == TOTAL .and. sw_pdos == ON) then

          allocate( porb_ssrep5( neg, np2, ndim_chgpot ) ); porb_ssrep5 = 0.0d0
          allocate( dos_weight3( neg, np2, ndim_magmom ) ); dos_weight3 = 0.0d0

          do iorb = 1, nlmta_phi
             call m_PP_tell_iorb_lmtt(iorb,lmtt)

             porb_ssrep5 = 0.0d0
! ---------------
             do ip2 = 1, np2
                ik = ndim_spinor*( ip2 -1 ) +1

                do ib = 1, neg
                   if ( k_symmetry(ik) == GAMMA ) then
                      call phase_error_with_msg(nfout,'Not supported ',__LINE__,__FILE__)
                   else
                      Do is1=1, ndim_spinor
                         Do is2=1, ndim_spinor
                            istmp = ( is1 -1 )*ndim_spinor + is2

                            ztemp = 0.0d0
                            do iopr=1,nopr
                               z1 = dcmplx( compr(ib,iorb,iopr,ik+is1-1 ), &
                                    &       compi(ib,iorb,iopr,ik+is1-1 ) )
                               z2 = dcmplx( compr(ib,iorb,iopr,ik+is2-1 ), &
                                    &       compi(ib,iorb,iopr,ik+is2-1 ) )
                               ztemp = ztemp + z1 *conjg(z2) &
                                    &      *( 1.d0+qorb(iorb)/norm_phig_mpi(lmtt,ip2) )
                            end do

                            porb_ssrep5(ib,ip2,istmp) = porb_ssrep5(ib,ip2,istmp) &
                                 &                   + ztemp /dble(nopr)
                         End do
                      End Do

                   endif

                end do
             end do
! ---------------
             call m_ES_DensMat_To_MagMom_porb( np2*neg, porb_ssrep5, dos_weight3 )
             dos_weight3 = dos_weight3 *wei
! ------------------------------------ - - -- -
             Do istmp=1, ndim_magmom
                call care_of_degenerate_state( neig, np2, 1, 1, eeig2, &
                     &                         dos_weight3(:,:,istmp) )
                                                               !dos_weight
                call nstt5i( ipridos, idim, Eminimum, Emaximum, nEwindows, &
                     &       nxyz_tetra, np2, np2, neig, eeig2(1,1,1), &
                     &       ip20, np0, np2, neig, mtetra, nttra, deltae_dos, &
                     &       dos_weight3(:,:,istmp), nEwindows_plus, &
                     &       doswk, dosinwk, pdos(0,iorb,istmp), &
                     &       sumpdos(0,iorb,istmp) )
             End Do

          end do

          deallocate( dos_weight3 );  deallocate( porb_ssrep5 )

       end if

       deallocate(dos_weight2);    deallocate(nttra);  deallocate(dosinwk,doswk)
    end if
! ----------------------------------------------------------

    deallocate(eeig2)

    if(icomponent == TOTAL) call get_VBM(totch,1.d-12)  ! -> ValenceBandMaximum

    if(mype == 0) then
       if ( calc_dos_magmom_contrib == YES ) then
          call write_dos_noncl(nfdos)
       else
          call write_dos_noncl_for_totchg(nfdos)
       endif
    endif

    if(mype == 0 .and. icomponent == TOTAL .and. sw_pdos == ON) then
       call write_pdos_noncl(nfdos)
    endif

    call dealloc_dos()
    call tstatc0_end(id_sname)
!!$#endif
  end subroutine m_ESdos_tetrahedral_noncl
! ==================================================================== 11.0

  subroutine doscal_nstt(nEwindows,nxyz_tetra,np2,neg,eeig2,ip20,np0,wei,cdos,cind)
    integer, intent(in) :: nEwindows, np2, neg, np0
    integer, intent(in), dimension(3) :: nxyz_tetra
    real(kind=DP), intent(in), dimension(np2,neg) ::eeig2
    integer, intent(in), dimension(np0) :: ip20
    real(kind=DP), intent(in)           :: wei
    real(kind=DP), dimension(np2,neg,0:nEwindows) :: cdos, cind

    integer :: mtetra
    integer, allocatable, dimension(:,:)     :: nttra  ! d(mtetra,4)
    integer, allocatable, dimension(:)       :: ip8    ! d(8)
    real(kind=DP), allocatable, dimension(:) :: voltt  ! d(mtetra)
    real(kind=DP), allocatable, dimension(:,:) :: vctk, vwork ! d(4,3)
    integer, allocatable, dimension(:)       :: k_sample_mesh !d(3)
    real(kind=DP), dimension(4)              :: e
    integer, parameter :: D12 = 1, D13 = 2, D14 = 3, D23 = 4, D24 = 5 &
         &         , D34 = 6, DM  = 7, DA  = 8, mdoswk = da
    real(kind=DP), allocatable, dimension(:,:) :: dos ! d(nEwindows,neg)
    real(kind=DP), allocatable, dimension(:)   :: d

    integer :: ix,iy,iz,icub,npx,npy,ni,ip0,kx,ky,kz,nx1,ny1,nz1,nd,nn,nt,i,j,nxx,nyy,nzz,nv(3),iv,ib,ie
    integer :: iswap1, iswap2, ns, ne, idos
    real(kind=DP) :: sumv,es,ee,v,x,f0,f1,f2,f3,ed

    mtetra = product(nxyz_tetra(1:3))*6
    allocate(nttra(mtetra,4))
    allocate(voltt(mtetra))

    allocate(ip8(8))
    npx = nxyz_tetra(1)+1
    npy = nxyz_tetra(2)+1
    icub = 0
    do iz = 0, nxyz_tetra(3)-1
       do iy = 0, nxyz_tetra(2)-1
          do ix = 0, nxyz_tetra(1)-1
             icub = icub+1
             ni=npx*(npy*iz+iy)+ix
             do kz=1,2
                do ky=1,2
                   do kx=1,2
                      ip0 = ni+npx*(npy*(kz-1)+ky-1)+kx
                      ip8(kx+2*(ky-1)+4*(kz-1)) = ip0
                   end do
                end do
             end do
             nttra((icub-1)*6+1, 1:4) = (/ip8(1),ip8(2),ip8(4),ip8(8)/)
             nttra((icub-1)*6+2, 1:4) = (/ip8(1),ip8(2),ip8(6),ip8(8)/)
             nttra((icub-1)*6+3, 1:4) = (/ip8(1),ip8(5),ip8(6),ip8(8)/)
             nttra((icub-1)*6+4, 1:4) = (/ip8(1),ip8(3),ip8(4),ip8(8)/)
             nttra((icub-1)*6+5, 1:4) = (/ip8(1),ip8(3),ip8(7),ip8(8)/)
             nttra((icub-1)*6+6, 1:4) = (/ip8(1),ip8(5),ip8(7),ip8(8)/)
          end do
       end do
    end do
    deallocate(ip8)

    nxx = nxyz_tetra(1)
    nyy = nxyz_tetra(2)
    nzz = nxyz_tetra(3)
    allocate(k_sample_mesh(3))
    call m_Kp_sample_mesh(k_sample_mesh)
    nx1 = max(k_sample_mesh(1),1)
    ny1 = max(k_sample_mesh(2),1)
    nz1 = max(k_sample_mesh(3),1)
    deallocate(k_sample_mesh)
    nd = nx1*ny1*nz1
    if(ipridos >=2 ) write(nfout,'(" nd = ",i8)') nd
    allocate(vctk(4,3),vwork(4,3))
    if(ipridos>=2) write(nfout,'(" trmat = ",9f8.4)') trmat
    sumv = 0.d0
    do nt = 1, mtetra
       if(ipridos>=2) write(nfout,'(" nt = ",i8," : ",4i8)') nt, nttra(nt,1:4)
       do i = 1,4
          nn = nttra(nt,i)
          ix = mod(nn-1,nxx+1)
          iy = mod((nn-ix-1)/(nxx+1),nyy+1)
          iz = (nn-ix-1-iy*(nxx+1))/((nxx+1)*(nyy+1))
          nv(1) = ix*ny1*nz1; nv(2) = nx1*iy*nz1; nv(3) = nx1*ny1*iz
          if(ipridos >=2) write(nfout,'("  nvtk(",i4,",1:3) = ",3i8)') i,nv(1:3)
          vctk(i,1:3) = matmul(trmat,nv)/dble(nd)
       end do
       if(ipridos >=2) then
          write(nfout,'("  v(1) = ",3f8.4)') vctk(1,1:3)
          write(nfout,'("  v(2) = ",3f8.4)') vctk(2,1:3)
          write(nfout,'("  v(3) = ",3f8.4)') vctk(3,1:3)
          write(nfout,'("  v(4) = ",3f8.4)') vctk(4,1:3)
       end if

       do j = 1, 3
          vwork(1,j) = vctk(4,j) - vctk(1,j)
          vwork(2,j) = vctk(4,j) - vctk(2,j)
          vwork(3,j) = vctk(4,j) - vctk(3,j)
       enddo
       !     outer production
       vwork(4,1) = vwork(1,2)*vwork(2,3) - vwork(1,3)*vwork(2,2)
       vwork(4,2) = vwork(1,3)*vwork(2,1) - vwork(1,1)*vwork(2,3)
       vwork(4,3) = vwork(1,1)*vwork(2,2) - vwork(1,2)*vwork(2,1)
       !     inner procution
       voltt(nt) = dabs(dot_product(vwork(3,1:3),vwork(4,1:3)))/6.d0
       sumv = sumv + voltt(nt)
       if(ipridos >= 2) write(nfout,'(3x,i3,2x,f18.12)') nt,voltt(nt)
    end do
    if(ipridos >= 2) write(nfout,'(" sumv = ",f16.8, " univol*sumv/(2pi)**3 = ",f16.8)') sumv,univol*sumv/(PAI2*PAI2*PAI2)
    if(ipridos >= 2) write(nfout,'(" np2 = ",i8)') np2

    es = minval(eeig2)-0.005
    ee = maxval(eeig2)+0.005

    allocate(d(mdoswk))
    allocate(dos(nEwindows,neg))
    do iv = 1, mtetra
       v = 6.0d0*voltt(iv)*univol/(PAI2*PAI2*PAI2*nspin)
       write(nfout,'(" iv = ",i8," v = ",d12.4)') iv,v
       dos = 0.d0
       do ib = 1 ,neg
          do ie = 1, 4
             if(ip20(nttra(iv,ie)) > np2) then
                write(nfout,'(" nttra(iv,ie) = ",i8)') nttra(iv,ie)
             end if
             e(ie) = eeig2(ip20(nttra(iv,ie)),ib)
          end do
          do iswap1 = 1, 3
             do iswap2 = iswap1+1,4
                if(e(iswap1) < e(iswap2)) then
                   x = e(iswap1)
                   e(iswap1) = e(iswap2)
                   e(iswap2) = x
                end if
             end do
          end do
          if(e(4) > ee) cycle
          NS = (E(4) - ES)/DELTAE_dos + 1
          IF(NS.GT.nEwindows) cycle
          NE = (E(1) - ES)/DELTAE_dos + 1
          IF(NE.GT.nEwindows) NE = nEwindows

          write(nfout,'(" (iv,ib)=(",i4,",",i4,"), e = ",4f8.4," ns, ne = ",2i8)') iv,ib,e(1:4),ns,ne
          IF((NE-NS) == 0) THEN
             DOS(NS,ib) = DOS(NS,ib) + V/DELTAE_dos
          ELSE IF((NE-NS).EQ.1 ) THEN
             DOS(NS,ib) = DOS(NS,ib) + V/(2.D0*DELTAE_dos)
             DOS(NE,ib) = DOS(NE,ib) + V/(2.D0*DELTAE_dos)
          ELSE
             D(D34) = E(3) - E(4)
             D(D24) = E(2) - E(4)
             D(D14) = E(1) - E(4)
             D(D23) = E(2) - E(3)
             D(D13) = E(1) - E(3)
             D(D12) = E(1) - E(2)
             D(DM) = D(D13) + D(D24)
             D(DA) = (E(1)*E(2)-E(3)*E(4))/D(DM)

             if(d(d34).lt.deltae_dos) then
                f0 = 0.d0
             else
                f0 = v/(d(d34)*d(d24)*d(d14))
             endif

             F1 = V/D(DM)
             if(d(d23).lt.deltae_dos) then
                f2 = 0.d0
             else
                f2 = v*d(dm)/(d(d24)*d(d14)*d(d23)*d(d13))
             endif
             if(d(d12).lt.deltae_dos) then
                f3 = 0.d0
             else
                f3 = v/(d(d14)*d(d13)*d(d12))
             endif

             DO IDOS = NS, NE
                ED = ES + DELTAE_dos*(IDOS-0.5)
                IF(ED.LE.E(3)) THEN
                   DOS(IDOS,ib) = DOS(IDOS,ib) + (ED - E(4))*(ED - E(4)) * F0
                ELSE IF(ED.LE.E(2)) THEN
                   DOS(IDOS,ib) = DOS(IDOS,ib) + F1 - F2*(ED - D(DA))*(ED - D(DA))
                ELSE IF(ED.LE.E(1)) THEN
                   DOS(IDOS,ib) = DOS(IDOS,ib) + (ED - E(1))*(ED - E(1)) * F3
                ENDIF
             end DO
          ENDIF
          write(nfout,'(" tetra = ",i8," ik = ",4i8," ns,ne=",2i8)') &
               & iv, (ip20(nttra(iv,ie)),ie=1,4), ns, ne
!!$),ip20(nttra(iv,2)),ip20(nttra(iv,3)),ip20(nttra(iv,4)),ns,ne
          write(nfout,'(" dos = ",6d12.4)') (dos(idos,ib),idos=min(ns+10,ne),min(ns+15,ne))
          do ie = 1, 4
             i = ip20(nttra(iv,ie))
             do idos = ns, ne
                cdos(i,ib,idos) = cdos(i,ib,idos) + dos(idos,ib)*0.25/wei
             end do
          end do
       end do
    end do
    do i = 1, np2
       do ib = 1, neg
          cind(i,ib,1) = cdos(i,ib,1)*deltaE_dos
          do idos = 2, nEwindows
             cind(i,ib,idos) = cind(i,ib,idos-1)+cdos(i,ib,idos)*deltaE_dos
          end do
       end do
    end do
    deallocate(d)
    deallocate(dos)
    deallocate(vctk,vwork)
    deallocate(voltt,nttra)
  end subroutine doscal_nstt

  subroutine m_ESdos_rd_pdos_param(nfout)
    use m_Const_Parameters, only : FMAXVALLEN, ON, NOCONV, LOWER &
          &                      , Projector, Wavefunction, Mulliken
    use m_Ionic_System, only : speciesname,ntyp
    integer, intent(in) :: nfout

    character(len=FMAXVALLEN) :: rstr
    integer :: iret, f_selectBlock, f_getStringValue, f_getIntValue
    integer :: f_selectTop, f_selectParentBlock
    real(kind=DP) :: dret
    logical :: prealloc, tag_is_found
    integer :: i,it,n

    iret = f_selectTop()

    if( f_selectBlock(tag_postprocessing) == 0) then

! --- pdos
       if( f_selectBlock(tag_pdos) == 0) then
          if(f_getIntValue( tag_sw_orb_popu, iret) == 0) sw_orb_popu = iret
          if(f_getIntValue( tag_sw_pdos, iret) == 0) then
             sw_pdos = iret
             if(sw_pdos == ON) sw_orb_popu = ON
!             if ( sw_pdos == ON ) sw_dos = ON
          end if
          if(f_getIntValue( tag_use_rotated_compri, iret) == 0) then
             use_rotated_compri = iret
          endif
!!$if(sw_orb_popu == ON) then
!!$   iret = f_getStringValue(tag_method,rstr,LOWER)
!!$   if( rstr == tag_projector) then
!!$      pdos_method = Projector
!!$   else if( rstr == tag_wavefunction) then
!!$      pdos_method = Wavefunction
!!$   else if( rstr == tag_mulliken) then
!!$      pdos_method = Mulliken
!!$   end if
!!$   if( f_selectBlock(tag_orbitals) == 0) then
!!$      prealloc = .true.
!!$      call set_orbitals(prealloc) ! --> norbital
!!$      prealloc = .false.
!!$      call set_orbitals(prealloc) ! --> l_orb,t_orb
!!$      iret = f_selectParentBlock()
!!$   else
!!$      stop ' orbitals are not given properly in the inputfile'
!!$   end if
!!$end if
          if(ipriinputfile >= 1) then
             write(nfout,'(" !** sw_pdos             = ",i3)') sw_pdos
             write(nfout,'(" !** sw_orb_popu         = ",i3)') sw_orb_popu
!!$write(nfout,'(" !** pdos_method         = ",i3)') pdos_method
!!$if(sw_orb_popu == ON) then
!!$   write(nfout,'(" !** === orbitals === ")')
!!$   write(nfout,'(" !**   no   l   t   rc       k         type")')
!!$   n = 0
!!$   do it=1,ntyp
!!$      do i=1,norbital(it)
!!$         n = n + 1
!!$         write(nfout,'(" !** ",3i4,2f10.5,a6)') &
!!$           & n,l_orb(i,it),t_orb(i,it),rc_orb(i,it),k_orb(i,it),speciesname(it)
!!$      end do
!!$   end do
!!$end if
          end if
          iret = f_selectParentBlock()
       end if             ! pdos
       iret = f_selectParentBlock()
    end if

  contains

    subroutine set_orbitals(prealloc)
      logical, intent(in) :: prealloc

      character(len=FMAXVALLEN) :: rstr
      integer :: i,iret,ip,rint,it
      integer :: f_selectFirstTableLine, f_selectNextTableLine &
           &, f_getIntValue, f_getRealValue, f_getStringValue
      real(kind=DP) :: dret
      logical :: first

      if(.not.prealloc) then
         maxorb = 0
         do it=1,ntyp
            maxorb = max(maxorb,norbital(it))
         end do
         allocate(l_orb(maxorb,ntyp))
         allocate(t_orb(maxorb,ntyp))
         allocate(rc_orb(maxorb,ntyp)); rc_orb=0.d0
         allocate(k_orb(maxorb,ntyp)); k_orb=0.d0
      else
         allocate(norbital(ntyp)); norbital = 0
      end if

      do it=1,ntyp
         i = 1
         first = .true.
         do while(.true.)
            if (first) then
               if(f_selectFirstTableLine() /= 0) then
                  exit
               end if
               first = .false.
            else
               if(f_selectNextTableLine() /= 0) then
                  exit
               end if
            end if
            iret = f_getStringValue(tag_element,rstr,NOCONV)
            if(rstr == speciesname(it)) then
               if(.not.prealloc) then
                  ip = i
                  if( f_getIntValue(tag_l, rint) == 0) l_orb(ip,it) = rint
                  if( f_getIntValue(tag_t, rint) == 0) t_orb(ip,it) = rint
                  if( f_getRealValue(tag_rc, dret, 'bohr') == 0) rc_orb(ip,it) = dret
                  if( f_getRealValue(tag_k, dret, '') == 0) k_orb(ip,it) = dret
               end if
               i = i + 1
            end if

         end do
         norbital(it) = i - 1
      end do

    end subroutine set_orbitals

  end subroutine m_ESdos_rd_pdos_param

  subroutine prepare_nttra(nxyz_tetra,mtetra,nttra)
    integer, intent(in),  dimension(3)        :: nxyz_tetra
    integer, intent(in)                       :: mtetra
    integer, intent(out), dimension(mtetra,4) :: nttra

    integer :: npx,npy,icub,ix,iy,iz,ni,kx,ky,kz,ip0
    integer, allocatable, dimension(:) :: ip8

    allocate(ip8(8))
    npx = nxyz_tetra(1)+1
    npy = nxyz_tetra(2)+1
    icub = 0
    do iz = 0, nxyz_tetra(3)-1
       do iy = 0, nxyz_tetra(2)-1
          do ix = 0, nxyz_tetra(1)-1
             icub = icub+1
             ni=npx*(npy*iz+iy)+ix
             do kz=1,2
                do ky=1,2
                   do kx=1,2
                      ip0 = ni+npx*(npy*(kz-1)+ky-1)+kx
                      ip8(kx+2*(ky-1)+4*(kz-1)) = ip0
                   end do
                end do
             end do
             nttra((icub-1)*6+1, 1:4) = (/ip8(1),ip8(2),ip8(4),ip8(8)/)
             nttra((icub-1)*6+2, 1:4) = (/ip8(1),ip8(2),ip8(6),ip8(8)/)
             nttra((icub-1)*6+3, 1:4) = (/ip8(1),ip8(5),ip8(6),ip8(8)/)
             nttra((icub-1)*6+4, 1:4) = (/ip8(1),ip8(3),ip8(4),ip8(8)/)
             nttra((icub-1)*6+5, 1:4) = (/ip8(1),ip8(3),ip8(7),ip8(8)/)
             nttra((icub-1)*6+6, 1:4) = (/ip8(1),ip8(5),ip8(7),ip8(8)/)
          end do
       end do
    end do
    deallocate(ip8)
  end subroutine prepare_nttra

  subroutine care_of_degenerate_state(neig,np2,nspin,ispin,eeig2,dos_weight2)
    integer, intent(in) :: neig,np2,nspin,ispin
!!$    real(kind=DP),intent(in),dimension(neig,np2,nspin) :: eeig2
    real(kind=DP),intent(in),dimension(np2,neig,nspin) :: eeig2
    real(kind=DP),intent(inout),dimension(neig,np2) :: dos_weight2(neig,np2)
    real(kind=DP) :: eps,c1
    integer :: k2, ieig, n, i, ie
! ---- following lines are the later part of subroutine nstt3i --
!     take care of a weight on a degenerate state
!
    eps = dfloat(10)**(-5)

    do k2=1,np2
       if(ipridos >= 3) write(nfout,'(" k2 = ",i5," <<care_of_degenerate_state>>")') k2
       ieig=1
40     continue
       n=1
!!$do 42 i=1,20
       do i=1,neig
          if(ieig+i.gt.neig) go to 44
          if(dabs(eeig2(k2,ieig+i,ispin)-eeig2(k2,ieig,ispin)).lt.eps) then
             n=n+1
             cycle
          end if
          go to 44
       end do

44     continue

       if(ipridos >= 3) write(nfout,'("       n = ",i5)')  n
       c1=0
       do i=0,n-1
          if(ieig+i>neig) then
             if(ipridos>=1) write(nfout,'(" ieig+i = ",i8," > neig = ",i8)') ieig+i,neig
             cycle
          end if
          c1=c1+dos_weight2(ieig+i,k2)
       end do
       c1=c1/n
       do i=0,n-1
          if(ieig+i>neig) then
             if(ipridos>=1) write(nfout,'(" ieig+i = ",i8," > neig = ",i8)') ieig+i,neig
             cycle
          end if
          dos_weight2(ieig+i,k2)=c1
       end do
       ieig=ieig+n
       if(ieig.lt.neig) go to 40
    end do
  end subroutine care_of_degenerate_state


  subroutine m_ESdos_rd_pcohp_param(nfout)
    use m_Const_Parameters, only : FMAXVALLEN, ON, LOWER
    integer, intent(in) :: nfout

    character(len=FMAXVALLEN) :: rstr
    integer :: iret, f_selectBlock, f_getStringValue, f_getIntValue, f_getRealValue
    integer :: f_selectTop, f_selectParentBlock
    real(kind=DP) :: dret
    logical :: prealloc, tag_is_found
    integer :: i,it,n

    iret = f_selectTop()

    if( f_selectBlock(tag_postprocessing) == 0) then
! ---  pcohp --
       if( f_selectBlock( tag_pcohp ) == 0) then
          if ( f_getIntValue( tag_sw_calc_pcohp, iret) == 0) sw_calc_pcohp = iret

          if ( sw_calc_pcohp == ON ) then
             if ( nopr > 1 ) call phase_error_with_msg(nfout,"Nopr > 1 is no supported",__LINE__,__FILE__)

             sw_dos = ON ;  sw_pdos = ON;  orb_popu_method = 2
             if(printable) then
                write(nfout,'(A)') "!** orb_popu_method is forced to be 2"
             endif

             if ( f_getIntValue( tag_atom1, iret) == 0 ) then
                if ( iret > 0 .and. iret <= natm ) then
                   atom1_pcohp = iret
                else
                   write(*,*) "set atom1 in the range [1: ", natm ,"]"
                   call phase_error_with_msg(nfout,'set atom1 in the range [1:natm]',__LINE__,__FILE__)
                endif
             endif
             if ( f_getRealValue( tag_max_distance_from_atom1, dret, "") == 0) then
                max_distance_from_atom1_pcohp = dret
             endif
             if(printable) then
                write(nfout,*) ' ** atom1 of pcohp is ', atom1_pcohp
                write(nfout,*) ' ** max_distance_from_atom1 is ', &
                     &              max_distance_from_atom1_pcohp
             endif
          endif
          iret = f_selectParentBlock()
       endif
! ----------
    end if
  end subroutine m_ESdos_rd_pcohp_param

  subroutine m_ES_calc_pcohp
    integer :: ia, ja, it1, it2
    integer :: ik, ispin, lmt1, lmt2
    integer :: i, is, ie, id
    real(kind=DP) :: c1, w, e, El, Eu, tl, tu, Es, DeltaE, distance
    complex(kind=CMPLDP) :: z1

    real(kind=DP) :: derf

    complex(kind=CMPLDP), allocatable :: mat_p(:,:,:)
    complex(kind=CMPLDP), allocatable :: mat_h(:,:)

    integer, parameter :: lun = 355

! --- init --
    allocate( mat_p( neg,nlmt_phi,nlmt_phi ) );
    allocate( mat_h( nlmt_phi,nlmt_phi ) );

    allocate( pcohp( 1:nEwindows,nlmt_phi,nlmt_phi,nspin ) );
    allocate( sum_pcohp( 1:nEWindows,nlmt_phi,nlmt_phi,nspin ) );

    DeltaE = DeltaE_dos
    Es = Eminimum - DeltaE*0.5d0

    if ( mype == 0 ) then
       open( unit=lun, file='pcohp.data', status='unknown', form ='formatted' )
    endif

! ---- begin -----
    ia = atom1_pcohp
    Do ja=1, natm
       if ( ia == ja ) cycle

       call calc_distance_from_atom1( ia, ja, distance )
       if ( distance > max_distance_from_atom1_pcohp ) cycle

       it1 = ityp(ia);  it2 = ityp(ja)

       pcohp = 0.0d0;     sum_pcohp = 0.0d0

       Do ik=1, kv3
          ispin = mod( ik -1, nspin ) +1
          call calc_matrix_p_nosym( ik, ia, ja, mat_p )
          call calc_matrix_h( ik, ia, ja, mat_p, mat_h )

          Do lmt1=1, ilmt_phi(it1)
             Do lmt2=1, ilmt_phi(it2)

                do i = 1, neg
                   w = 1.d0
                   z1 = mat_p(i,lmt1,lmt2) *conjg( mat_h(lmt1,lmt2) )
                   c1 = real(z1)

                   e = eko(i,ik) -nETailWindows*DeltaE -Es
                   is = e /DeltaE;      ie = is +2*nETailWindows

                   if (is < 0) is = 0
                   if (ie >= nEWindows )ie = nEWindows-1

                   do id = is, ie
                      El = Es + id*DeltaE;    Eu = El + DeltaE
                      tl = (El -eko(i,ik))*sqrtdVI;   tu = (Eu -eko(i,ik))*sqrtdVI

                      pcohp(id+1,lmt1,lmt2,ispin) = pcohp(id+1,lmt1,lmt2,ispin) &
                           & + c1 *w *2 *(derf(tu) -derf(tl)) *0.5d0 /DeltaE &
                           &            *qwgt(ik)
                   end do
                end do
             end Do
          end Do
          ! ---- integ --
          Do lmt1=1, ilmt_phi(it1)
             Do lmt2=1, ilmt_phi(it2)
                sum_pcohp(1,lmt1,lmt2,ispin) = pcohp(1,lmt1,lmt2,ispin)*DeltaE
                do id = 1, nEWindows-1
                   sum_pcohp(id+1,lmt1,lmt2,ispin) &
                        & = sum_pcohp(id,lmt1,lmt2,ispin) &
                        &   +pcohp(id+1,lmt1,lmt2,ispin) *DeltaE
                end do
             End Do
          ENd Do

       end Do   ! ik
       ! --
       if ( mype == 0 )  then
!          call write_pcohp_full( lun, ia, ja )
          call write_pcohp_sum_over_atom( lun, ia, ja )
       endif
    End Do

    ! --- finalize--
    if ( mype == 0 ) close( lun )

    deallocate( mat_p ); deallocate( mat_h )
    deallocate( pcohp ); deallocate( sum_pcohp )

  contains

    subroutine calc_distance_from_atom1( ia, ja, dist_min )
      integer, intent(in) :: ia, ja
      real(kind=DP), intent(out) :: dist_min

      integer :: nx, ny, nz
      real(kind=DP) :: x1, y1, z1, dist

      dist_min = 1.0D3
      Do nx=-1, 1
         Do ny=-1, 1
            Do nz=-1, 1
               x1 = cps(ja,1) -cps(ia,1) +nx*altv(1,1) +ny*altv(1,2) +nz*altv(1,3)
               y1 = cps(ja,2) -cps(ia,2) +nx*altv(2,1) +ny*altv(2,2) +nz*altv(2,3)
               z1 = cps(ja,3) -cps(ia,3) +nx*altv(3,1) +ny*altv(3,2) +nz*altv(3,3)
               dist = sqrt( x1**2 +y1**2 +z1**2 )
               dist_min = min( dist_min, dist )
            End Do
         ENd Do
      End Do
    end subroutine calc_distance_from_atom1

  end subroutine m_ES_calc_pcohp

  subroutine write_pcohp_full( nf, ia, ja )
    integer, intent(in) :: nf, ia, ja

    integer :: it1, it2, lmt1, lmt2, il1, il2, im1, im2, tau1, tau2
    integer :: id, nrEWindows, ii, i, nwdwidth
    real(kind=DP) :: e, e_eV, pcohp_hr, pcohp_hr2, sum_pcohp_avr, sum_pcohp_avr2
    real(kind=DP) :: pcohp_eV, pcohp_eV2, sumtotal

    if(nwd_dos_window_width >= 1 .and. nwd_dos_window_width <= nEWindows) then
       nwdwidth = nwd_dos_window_width
    else
       nwdwidth = 1
    end if
    id = nwdwidth/2
    nrEWindows = nEWindows/nwdwidth * nwdwidth

    if ( iproj_group(ia)== 0 ) return
    if ( iproj_group(ja)== 0 ) return

    if ( if_pdos(ia)== 0 ) return
    if ( if_pdos(ja)== 0 ) return

    it1 = ityp(ia);     it2 = ityp(ja)

    Do lmt1=1, ilmt_phi(it1)
       il1 = ltp_phi(lmt1,it1); im1 = mtp_phi(lmt1,it1);  tau1 = taup_phi(lmt1,it1)

       Do lmt2=1, ilmt_phi(it2)
          il2 = ltp_phi(lmt2,it2); im2 = mtp_phi(lmt2,it2);  tau2 = taup_phi(lmt2,it2)

          write(nf,'(A,I5,3(A,I3),A,I5,3(A,I3))') &
               &       "-PCOHP: ia= ", ia, " l=", il1-1," m=", im1," t=", tau1, &
               &       "        ja= ", ja, " l=", il2-1," m=", im2," t=", tau2

          if (nspin == 1) then
             write(nf,*) "     E(eV)          pcohp(eV)           sum(eV)"
             do i = 1, nrEWindows, nwdwidth
                e = Eminimum + (i-1+id)*DeltaE_dos
                e_eV = (e - ValenceBandMaximum)*Hartree
                pcohp_hr = 0.d0;  sum_pcohp_avr = 0.d0

                do ii = i, i+nwdwidth-1
                   pcohp_hr = pcohp_hr + pcohp(ii,lmt1,lmt2,nspin)
                   sum_pcohp_avr = sum_pcohp_avr + sum_pcohp(ii,lmt1,lmt2,nspin)
                end do
                pcohp_hr = pcohp_hr /nwdwidth;
                sum_pcohp_avr = sum_pcohp_avr /nwdwidth
!                pcohp_eV = pcohp_hr /Hartree    ! ???
                pcohp_eV = pcohp_hr

                write(nf,'(f14.6,f18.10,f20.10)') &
                     &     e_eV, -pcohp_eV, -sum_pcohp_avr
             end do

          else if(nspin == 2) then
             write(nf,*) "   E(eV)         pcohp_up(eV)        pcohp_down(eV)" &
                  & ,"   sum_up    sum_down  sum_total"

             do i = 1, nrEWindows, nwdwidth
                e = Eminimum + (i-1+id)*DeltaE_dos
                e_eV = (e - ValenceBandMaximum)*Hartree
                pcohp_hr = 0.d0; pcohp_hr2 = 0.d0;
                sum_pcohp_avr = 0.d0; sum_pcohp_avr2 = 0.d0

                do ii = i, i+nwdwidth-1
                   pcohp_hr  = pcohp_hr  + pcohp(ii,lmt1,lmt2,1)
                   pcohp_hr2 = pcohp_hr2 + pcohp(ii,lmt1,lmt2,2)
                   sum_pcohp_avr  = sum_pcohp_avr  + sum_pcohp(ii,lmt1,lmt2,1)
                   sum_pcohp_avr2 = sum_pcohp_avr2 + sum_pcohp(ii,lmt1,lmt2,2)
                end do
                pcohp_hr = pcohp_hr/nwdwidth; pcohp_hr2 = pcohp_hr2/nwdwidth
                sum_pcohp_avr = sum_pcohp_avr/nwdwidth;
                sum_pcohp_avr2 = sum_pcohp_avr2/nwdwidth

!                pcohp_eV  = pcohp_hr/Hartree
!                pcohp_eV2 = pcohp_hr2/Hartree
                pcohp_eV  = pcohp_hr;
                pcohp_eV2 = pcohp_hr2
                sumtotal = sum_pcohp_avr + sum_pcohp_avr2
                write(nf,'(f15.4,2f18.10,3f10.4)') &
                     &    e_eV, -pcohp_eV, -pcohp_eV2, &
                     &    -sum_pcohp_avr, -sum_pcohp_avr2, -sumtotal
             end do
          end if
          write(nf,*) 'END'
       end Do
    end Do

  end subroutine write_pcohp_full

  subroutine write_pcohp_sum_over_atom( nf, ia, ja )
    integer, intent(in) :: nf, ia, ja

    integer :: it1, it2, lmt1, lmt2, il1, il2, im1, im2, tau1, tau2
    integer :: id, nrEWindows, ii, i, nwdwidth
    real(kind=DP) :: e, e_eV, pcohp_hr, pcohp_hr2, sum_pcohp_avr, sum_pcohp_avr2
    real(kind=DP) :: pcohp_eV, pcohp_eV2, sumtotal

    if(nwd_dos_window_width >= 1 .and. nwd_dos_window_width <= nEWindows) then
       nwdwidth = nwd_dos_window_width
    else
       nwdwidth = 1
    end if
    id = nwdwidth/2
    nrEWindows = nEWindows/nwdwidth * nwdwidth

    if ( iproj_group(ia)== 0 ) return
    if ( iproj_group(ja)== 0 ) return

    if ( if_pdos(ia)== 0 ) return
    if ( if_pdos(ja)== 0 ) return

    it1 = ityp(ia);     it2 = ityp(ja)

    write(nf,'(A,I5,A,I5)') &
         &       "-PCOHP: ia= ", ia, "        ja= ", ja
    if (nspin == 1) then
       write(nf,*) "     E(eV)          pcohp(eV)           sum(eV)"
       do i = 1, nrEWindows, nwdwidth
          e = Eminimum + (i-1+id)*DeltaE_dos
          e_eV = (e - ValenceBandMaximum)*Hartree
          pcohp_hr = 0.d0;  sum_pcohp_avr = 0.d0

          Do lmt1=1, ilmt_phi(it1)
             Do lmt2=1, ilmt_phi(it2)
                do ii = i, i+nwdwidth-1
                   pcohp_hr = pcohp_hr + pcohp(ii,lmt1,lmt2,nspin)
                   sum_pcohp_avr = sum_pcohp_avr + sum_pcohp(ii,lmt1,lmt2,nspin)
                end do
             End Do
          End Do
          pcohp_hr = pcohp_hr /nwdwidth;
          sum_pcohp_avr = sum_pcohp_avr /nwdwidth
!                pcohp_eV = pcohp_hr /Hartree    ! ???
          pcohp_eV = pcohp_hr
          sum_pcohp_avr = sum_pcohp_avr *Hartree

          write(nf,'(f14.6,f18.10,f20.10)') &
               &     e_eV, -pcohp_eV, -sum_pcohp_avr
       end do
    else if(nspin == 2) then
       write(nf,*) "   E(eV)         pcohp_up(eV)        pcohp_down(eV)" &
            & ,"sum_up   sum_down sum_total"
       do i = 1, nrEWindows, nwdwidth
          e = Eminimum + (i-1+id)*DeltaE_dos
          e_eV = (e - ValenceBandMaximum)*Hartree
          pcohp_hr = 0.d0; pcohp_hr2 = 0.d0;
          sum_pcohp_avr = 0.d0; sum_pcohp_avr2 = 0.d0

          Do lmt1=1, ilmt_phi(it1)
             Do lmt2=1, ilmt_phi(it2)
                do ii = i, i+nwdwidth-1
                   pcohp_hr  = pcohp_hr  + pcohp(ii,lmt1,lmt2,1)
                   pcohp_hr2 = pcohp_hr2 + pcohp(ii,lmt1,lmt2,2)
                   sum_pcohp_avr  = sum_pcohp_avr  + sum_pcohp(ii,lmt1,lmt2,1)
                   sum_pcohp_avr2 = sum_pcohp_avr2 + sum_pcohp(ii,lmt1,lmt2,2)
                end do
             end Do
          end Do
          pcohp_hr = pcohp_hr/nwdwidth; pcohp_hr2 = pcohp_hr2/nwdwidth
          sum_pcohp_avr = sum_pcohp_avr/nwdwidth;
          sum_pcohp_avr2 = sum_pcohp_avr2/nwdwidth

!                pcohp_eV  = pcohp_hr/Hartree
!                pcohp_eV2 = pcohp_hr2/Hartree
          pcohp_eV  = pcohp_hr
          pcohp_eV2 = pcohp_hr2

          sum_pcohp_avr  = sum_pcohp_avr  *Hartree
          sum_pcohp_avr2 = sum_pcohp_avr2 *Hartree

          sumtotal = sum_pcohp_avr + sum_pcohp_avr2
          write(nf,'(f15.4,2f18.10,3f10.4)') &
               &    e_eV, -pcohp_eV, -pcohp_eV2, &
               &    -sum_pcohp_avr, -sum_pcohp_avr2, -sumtotal
       end do
    end if
    write(nf,*) 'END'
  end subroutine write_pcohp_sum_over_atom

  subroutine calc_matrix_H( ik, ia, ja, mat_p, mat_H )
    integer, intent(in) :: ia, ja, ik
    complex(kind=CMPLDP), intent(in)  :: mat_p(neg,nlmt_phi,nlmt_phi)
    complex(kind=CMPLDP), intent(out) :: mat_h(nlmt_phi,nlmt_phi)

    integer :: it1, it2, lmt1, lmt2, ie
    complex(kind=CMPLDP) :: z1

    it1 = ityp(ia);     it2 = ityp(ja)

    mat_h = 0.0d0

    Do lmt1=1, ilmt_phi(it1)
       Do lmt2=1, ilmt_phi(it2)
          z1 = 0.0d0
          Do ie=1, neg
!             z1 = z1 + ( eko(ie,ik) -ValenceBandMaximum ) *mat_p(ie,lmt1,lmt2)
             z1 = z1 + eko(ie,ik) *mat_p(ie,lmt1,lmt2)
          End Do
          mat_h(lmt1,lmt2) = z1
       End Do
    End Do
  end subroutine calc_matrix_H

  subroutine calc_matrix_P_nosym( ik, ia, ja, mat_p )
    integer, intent(in) :: ia, ja, ik
    complex(kind=CMPLDP), intent(out) :: mat_p(neg,nlmt_phi,nlmt_phi)

    integer :: it1, it2, iorb1, iorb2, lmt1, lmt2, iopr1, iopr2, ie, i
    complex(kind=CMPLDP) :: porb, z1, z2
    logical :: First = .true.

    it1 = ityp(ia);     it2 = ityp(ja)

    mat_p = 0.0d0

    Do ie=1, neg
       Do lmt1=1, ilmt_phi(it1)
          iorb1 = lmta_phi(lmt1,ia)

          Do lmt2=1, ilmt_phi(it2)
             iorb2 = lmta_phi(lmt2,ja)

             porb = 0.0d0
             if ( k_symmetry(ik) == GAMMA ) then
                do iopr1=1, 1
                   porb = porb +compr(ie,iorb1,iopr1,ik) &
                        &      *compr(ie,iorb2,iopr1,ik) /2.0d0
                end do
             else
                Do iopr1=1, 1
                   z1 = cmplx( compr(ie,iorb1,iopr1,ik), compi(ie,iorb1,iopr1,ik) )
                   z2 = cmplx( compr(ie,iorb2,iopr1,ik), compi(ie,iorb2,iopr1,ik) )
                   porb = porb +conjg(z1)* z2
                end do
             endif
             mat_p( ie, lmt1, lmt2 ) = porb

          ENd Do
       End Do
    End Do
  end subroutine calc_matrix_P_nosym

! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
!
  subroutine nstt5i_3D(ipri,idim,es,ee,newindows,nxyz_tetra, &
       &     np2,lmnp2e,neig, &
       &     eeig,ip20,np0,lmnp2c,lmneig,mtetra,nttra,deltae,dos_weight, &
       &     nep,dos,dosin,cdos,csumdos)
!
!     nstt0i, nstt1i, nstt2i and nstt3i are merged into
!  this subroutine nstt5i
!               by T. Yamasaki, Aug 2007
!
    integer, intent(in) :: ipri,idim, newindows
    real(kind=DP), intent(in) :: es,ee, deltae
    integer, intent(in) ::  nxyz_tetra(3)
    integer, intent(in) :: np2,lmnp2e,neig,lmneig,np0,mtetra,lmnp2c,nep
    real(kind=DP), intent(in), dimension(lmnp2e,lmneig) :: eeig
    integer, intent(in), dimension(np0)  :: ip20
    integer, intent(in), dimension(mtetra,4) :: nttra
    real(kind=DP), intent(in), dimension(lmneig,lmnp2c) :: dos_weight
    real(kind=DP), intent(out), dimension(0:nep,4) :: dos,dosin
    real(kind=DP), intent(out), dimension(0:newindows) :: cdos,csumdos

    integer :: ieb(4)

    real(kind=DP), dimension(4) :: eb
    integer :: ncounter,ncounter1,ncounter2,ncounter3
    real(kind=DP) :: etime1,etime2,etime3, wct_now, wct_start,etime4,etime5,etime0
    real(kind=DP) :: d21,d31,d41,d32,d42,d43,e,d1,d2,d3,d4,yy,x,y,xx,td21,yy3,d2d2,d2d2d2,d3d3 &
         &          ,x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,esum,e1,e2,e3,e4,w,tdos
    integer :: ierr,itmp,nxx,nyy,nzz,npx,npy,npz,np,ncub,ntet,iloop,iv,ne1,ns,ns2,ns3,ne2,idos,ie,ib,ip,ne, ib_t

    nxx = nxyz_tetra(1); nyy = nxyz_tetra(2); nzz = nxyz_tetra(3)
    npx = nxx+1
    npy = nyy+1
    npz = nzz+1
    np = npx*npy*npz
    ncub = nxx*nyy*nzz
    ntet = 6*ncub

    dos = 0.0d0
    dosin = 0.0d0
    cdos = 0.0d0
    csumdos = 0.0d0
    if(ipri.ge.2) write(nfout,'(" np = ",i5," np2 = ",i5," idim = ",i5," ntet = ",i5)') np,np2,idim,ntet

    ncounter = 0
    ncounter1 = 0
    ncounter2 = 0
    ncounter3 = 0

    etime0 = 0.d0
    etime1 = 0.d0
    etime2 = 0.d0
    etime3 = 0.d0
    etime4 = 0.d0
    etime5 = 0.d0

    iloop = 0
    do iv = 1, ntet
       if(mod(iv,npes)/=(mype+1)) cycle
       if(ipri.ge.2) then
          write(nfout,'(" iv = ",i8)') iv
          write(nfout,'(" -- nttra  = ",4i8)') (nttra(iv,ie),ie=1,4)
       endif
!!$       do ib_t = 1, np_e
!!$          ib = ib_t + ista_e - 1
       do ib = 1, neig
          call gettod(wct_start)
          do ie = 1,4
             ieb(ie) = ip20(nttra(iv,ie))
             eb(ie) = eeig(ieb(ie),ib)
          end do

          call nsttod(eb,ieb)

          e1 = eb(1)
          e2 = eb(2)
          e3 = eb(3)
          e4 = eb(4)
          call nstts1(e1,e2,e3,e4)
          if(ipri.ge.3) then
             write(nfout,'(" (iv,ib) = (",2i5,") e(1:4) = ",4f9.5, " ieb(1:4) = ",4i3)') &
                  & iv,ib,e1,e2,e3,e4 ,ieb(1),ieb(2),ieb(3),ieb(4)
          end if
          if(e1 > ee) cycle
          ns = (e1 - es)/deltae
          if(es+deltae*ns .lt. e1) ns = ns + 1
          ne = (e4 - es)/deltae+1
          if(es+deltae*ne .ge. e4) ne = ne - 1
          if(es+deltae*ne .ge. e4) ne = ne - 1

          iloop = iloop + (ne-ns+1)
          ne1 = (e2-es)/deltae
          if(es+deltae*ne1>e2) ne1 = ne1-1
          ns2 = ne1+1
          ns3 = (e3-es)/deltae+1
          if(es+deltae*ns3 < e3) ns3 = ns3+1
          ne2 = ns3-1

          if(ipri .ge. 2) then
             write(nfout,'(" ns, ne1 = ",2i8," e = ",f9.5," - ",f9.5)')  ns,ne1,es+deltae*ns,es+deltae*ne1
             write(nfout,'(" ns2,ne2 = ",2i8," e = ",f9.5," - ",f9.5)')  ns2,ne2,es+deltae*ns2,es+deltae*ne2
             write(nfout,'(" ns3,ne  = ",2i8," e = ",f9.5," - ",f9.5)')  ns3,ne,es+deltae*ns3,es+deltae*ne
             if(es+deltae*ns <e1)write(nfout,'(" !! es+deltae+ns  < e1")')
             if(es+deltae*ne1>e2)write(nfout,'(" !! es+deltae*ne1 > e2")')
             if(es+deltae*ns2<e2)write(nfout,'(" !! es+deltae*ns2 < e2")')
             if(es+deltae*ne2>e3)write(nfout,'(" !! es+deltae*ne2 > e3")')
             if(es+deltae*ns3<e3)write(nfout,'(" !! es+deltae*ns3 < e3")')
             if(es+deltae*ne >e4)write(nfout,'(" !! es+deltae*ne  > e4")')
          end if

          call gettod(wct_now)
          etime0 = etime0 + (wct_now-wct_start)*1.d-6

          call gettod(wct_start)

          d21=e2-e1
          d31=e3-e1
          d41=e4-e1
          d32=e3-e2
          d42=e4-e2
          d43=e4-e3
          ncounter1 = ncounter1+(ne1-ns+1)
          do idos = ns, ne1
             e = es + deltae*idos
             d1=e-e1
             d4=e4-e
             d2=e2-e
             d3=e3-e
             yy=d41*d31*d21
             x=d2/d21+d3/d31+d4/d41
             y=(d1*d1)/yy
             dos(idos,1)=x*y
             dosin(idos,1)=0.25*d1*y*(x+1.0)
             xx=d1*d1*d1
             x=xx/(d21*yy)
             dos(idos,2)=x
             dosin(idos,2)=0.25*d1*x
             x=xx/(d31*yy)
             dos(idos,3)=x
             dosin(idos,3)=0.25*d1*x
             x=xx/(d41*yy)
             dos(idos,4)=x
             dosin(idos,4)=0.25*d1*x
          end do
          call gettod(wct_now)
          etime1 = etime1 + (wct_now-wct_start)*1.d-6

          call gettod(wct_start)
          ncounter3 = ncounter3+(ne-ns3+1)
          do idos = ns3,ne
             e = es + deltae*idos
             d1=e-e1
             d4=e4-e
             d2=e-e2
             d3=e-e3
             xx=d4*d4*d4
             yy=d41*d42*d43
             x=xx/(d41*yy)
             dos(idos,1)=x
             dosin(idos,1)=0.25*(1.0-d4*x)
             x=xx/(d42*yy)
             dos(idos,2)=x
             dosin(idos,2)=0.25*(1.0-d4*x)
             x=xx/(d43*yy)
             dos(idos,3)=x
             dosin(idos,3)=0.25*(1.0-d4*x)
             x=d3/d43+d2/d42+d1/d41
             y=(d4*d4)/yy
             dos(idos,4)=x*y
             dosin(idos,4)=0.25*(1.0-d4*y*(x+1.0))
          end do
          call gettod(wct_now)
          etime3 = etime3 + (wct_now-wct_start)*1.d-6

          call gettod(wct_start)
          ncounter2 = ncounter2+(ne2-ns2+1)
          td21 = 3.0*d21
          yy3=1.0/(d31*d42)+1.0/(d41*d32)
          do idos=ns2,ne2
             e = es + deltae*idos
!!$               e = es+deltae*(idos-0.5)
             d1=e-e1
             d4=e4-e
             d2=e-e2
             d3=e3-e
             d2d2 = d2*d2
             d2d2d2 = d2d2*d2
             d3d3 = d3*d3
             x1=(d3d3)/(d31*d31*d32)*(d2/d42+d1/d41)
             x2=(d4*d4)/(d41*d41*d42)*(d2/d32+d1/d31)
             x3=(d3*d4*d1)/(d31*d41)*yy3
             dos(idos,1)=0.5*(x1+x2+x3)
             x1=d2d2*(d32*d2+3.0*d3*(d32+d3))/12.0
             y1=x1
             x1=x1/(d31*d31*d32*d42)
             x2=d2*(d2d2*(d31+td21)+3.0*d3*(d2*d3+d32*(td21+d1)))
             x2=x2/12.0
             y2=x2
             x2=x2/(d31*d31*d32*d41)
             x3=d2d2*(d42*d2+3.0*d4*(d42+d4))/12.0
             y3=x3
             x3=x3/(d41*d41*d42*d32)
             x4=d2*(d2d2*(d41+td21)+3.0*d4*(d2*d4+d42*(td21+d1)))
             x4=x4/12.0
             y4=x4
             x4=x4/(d41*d41*d42*d31)
             x5=0.5*d2*d3*d4*(d1+d21)
             x5=x5+d2d2*(2.0*d21*(d3+d42) +(d1+d21)*(2.0*d3+d4+d42))/12.0
             x5=x5*yy3/(d31*d41)
             x6=0.25*d21*d21*(d42/d41+d32/d31+1.0)/(d41*d31)
             dosin(idos,1)=0.5*(x1+x2+x3+x4+x5)+x6

             x1=(d3d3)/(d32*d32*d31)*(d2/d42+d1/d41)
             x2=(d4*d4)/(d42*d42*d41)*(d2/d32+d1/d31)
             x3=(d3*d4*d2)/(d42*d32)*yy3
             dos(idos,2)=0.5*(x1+x2+x3)
             x1=y1/(d32*d32*d31*d42)
             x2=y2/(d32*d32*d31*d41)
             x3=y3/(d42*d42*d41*d32)
             x4=y4/(d42*d42*d41*d31)
             x5=d2d2*(d3*(d42+3.0*d4)+d32*(d42+d4))/12
             x5=x5*yy3/(d42*d32)
             x6=0.25*d21/d31*d21/d41
             dosin(idos,2)=0.5*(x1+x2+x3+x4+x5)+x6

             x1=(d2d2)/(d32*d32*d42)*(d3/d31+d4/d41)
             x2=(d1*d1)/(d31*d31*d41)*(d3/d32+d4/d42)
             x3=(d1*d2*d3)/(d32*d31)*yy3
             dos(idos,3)=0.5*(x1+x2+x3)
             x1=d2d2d2*(3*d3+d32)/12
             y1=x1
             x1=x1/(d32*d32*d42*d31)
             x2=d2d2d2*(3*d4+d42)/12
             y2=x2
             x2=x2/(d32*d32*d42*d41)
             x3=d2*(d2*d31*(d2+3*d21)+3*d3*(d2d2+3*d21*d1)+ 3*d21*d21*d32)/12
             y3=x3
             x3=x3/(d31*d31*d41*d32)
             x4=d2*(d2*d41*(d2+3*d21)+3*d4*(d2d2+3*d21*d1)+3*d21*d21*d42)/12
             y4=x4
             x4=x4/(d31*d31*d41*d42)
             x5=(d2d2)*(d3*(d21+3*d1)+d32*(d21+d1))/12
             x5=x5/(d32*d31)*yy3
             x6=0.25*(d21/d31*d21/d31*d21/d41)
             dosin(idos,3)=0.5*(x1+x2+x3+x4+x5)+x6

             x1=(d2d2)/(d42*d42*d32)*(d3/d31+d4/d41)
             x2=(d1*d1)/(d41*d41*d31)*(d3/d32+d4/d42)
             x3=(d1*d2*d4)/(d41*d42)*yy3
             dos(idos,4)=0.5*(x1+x2+x3)
             x1=y1/(d42*d42*d32*d31)
             x2=y2/(d42*d42*d32*d41)
             x3=y3/(d41*d41*d31*d32)
             x4=y4/(d41*d41*d31*d42)
             x5=(d2d2)*(d4*(d21+3*d1)+d42*(d21+d1))/12
             x5=x5/(d41*d42)*yy3
             x6=0.25*(d21/d31*d21/d41*d21/d41)
             dosin(idos,4)=0.5*(x1+x2+x3+x4+x5)+x6
          end do
          call gettod(wct_now)
          etime2 = etime2 + (wct_now-wct_start)*1.d-6

          call gettod(wct_start)
          if(idim .eq. -3) then
             esum = e1+e2+e3+e4
             do idos = ns, ne
                tdos = 0.025d0*(dos(idos,1)+dos(idos,2)+dos(idos,3)+dos(idos,4))
                dosin(idos,1)=dosin(idos,1)+tdos*(esum-4.d0*e1)
                dosin(idos,2)=dosin(idos,2)+tdos*(esum-4.d0*e2)
                dosin(idos,3)=dosin(idos,3)+tdos*(esum-4.d0*e3)
                dosin(idos,4)=dosin(idos,4)+tdos*(esum-4.d0*e4)
             end do
          end if
          do ip = 1,4
             do ie = ns, ne
                cdos(ie) = cdos(ie) + dos(ie,ip)*dos_weight(ib,ieb(ip))
                csumdos(ie) = csumdos(ie) + dosin(ie,ip)*dos_weight(ib,ieb(ip))
             end do
          end do
          w = (dos_weight(ib,ieb(1))+dos_weight(ib,ieb(2)) &
               &           +dos_weight(ib,ieb(3))+dos_weight(ib,ieb(4)))*0.25d0
          do ie = ne+1,nEwindows
             csumdos(ie) = csumdos(ie) + 1.d0*w
          end do
          call gettod(wct_now)
          etime4 = etime4 + (wct_now-wct_start)*1.d-6
       end do
    end do

    call gettod(wct_start)
    do ie = 0, nEwindows
       cdos(ie) = cdos(ie)/ntet
       csumdos(ie) = csumdos(ie)/ntet
    end do

!!$    if(nrank_e .gt.1) then
!!$       itmp = (nep+1)*4
!!$       call mpi_allreduce(MPI_IN_PLACE,dos,    itmp,       mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
!!$       call mpi_allreduce(MPI_IN_PLACE,dosin,  itmp,       mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
!!$       call mpi_allreduce(MPI_IN_PLACE,cdos,   newindows+1,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
!!$       call mpi_allreduce(MPI_IN_PLACE,csumdos,newindows+1,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
!!$    end if

    if (npes .gt. 1) then
       itmp = (nep+1)*4
       call mpi_allreduce(MPI_IN_PLACE,dos,  itmp,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       call mpi_allreduce(MPI_IN_PLACE,dosin,itmp,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       call mpi_allreduce(MPI_IN_PLACE,cdos,   newindows+1, mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       call mpi_allreduce(MPI_IN_PLACE,csumdos,newindows+1, mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
    endif

    call gettod(wct_now)
    etime5 = etime5 + (wct_now-wct_start)*1.d-6

    if(ipri .ge. 2) then
       write(nfout,'(" !dos iloop    = ",i12)') iloop
       write(nfout,'(" !dos ncounter1 = ",i12)') ncounter1
       write(nfout,'(" !dos ncounter2 = ",i12)') ncounter2
       write(nfout,'(" !dos ncounter3 = ",i12)') ncounter3
       write(nfout,'(" !dos etime0    = ",f16.8)') etime0
       write(nfout,'(" !dos etime1    = ",f16.8)') etime1
       write(nfout,'(" !dos etime2    = ",f16.8)') etime2
       write(nfout,'(" !dos etime3    = ",f16.8)') etime3
       write(nfout,'(" !dos etime4    = ",f16.8)') etime4
       write(nfout,'(" !dos etime5    = ",f16.8)') etime5
    end if
  end subroutine nstt5i_3D

end module m_ES_dos

