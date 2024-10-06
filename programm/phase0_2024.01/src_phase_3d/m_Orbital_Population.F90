!#define _DEBUG_WRITE_DFTU_MPI_PROCESSES_
!=======================================================================
!
!  SOFTWARE NAME : PHASE (ver. 7.01)
!
!  MODULE: m_Orbital_Population
!
!  AUTHOR(S): T. Yamamoto   May/04/2005
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
!=======================================================================
!   patch 0.1 by K. Tagami @adv    2009/01/19
!   patch 0.2 by K. Tagami @adv    2009/01/19
!   patch 0.2a by K. Tagami @adv    2009/02/11
!   patch 0.2b by K. Tagami @adv    2009/02/15
!
! ==========================
!   patch 0.1 :  correction in m_OP_occ_mat_ylm
!   patch 0.2 :  introduction of m_OP_occ_mat_ao_kt
!                ( because the implementation of "m_ES_orbital_den_mat"
!                  in Electronic_Structure.F90 is unclear )
!   patch 0.2a : output zero values in the case of af=1
!                  in subroutine  wd_occ_mat
!   patch 0.2b : bug fix in patch level 0.2a ( if af/=0 statement )
! =========================================================
!
!
module m_Orbital_Population
! $Id: m_Orbital_Population.F90 618 2020-05-12 08:44:53Z ktagami $
  use m_Const_Parameters,     only : DP,ON,ANEW,RENEW,SIMPLE,BROYD1,BROYD2,DFP,PULAY,OFF,UNIT_MATRIX
  use m_Files,                only : nfout, nfoccmat, m_Files_open_nfoccmat
  use m_IterationNumbers,     only : iteration_for_cmix
  use m_Control_Parameters,   only : nspin,af,ipri,printable,num_projectors &
     &                             , proj_attribute &
     &                             , proj_group, num_proj_elems, max_projs &
     &                             , iprihubbard &
     &                             , hownew,nbxmix,istrbr,ipripulay, sw_force_simple_mixing_hub &
     &                             , sw_metric_diff, alpha_pulay, sw_recomposing &
     &                             , sw_mix_bothspins_sametime, sw_force_simple_mixing
  use m_Ionic_System,         only : natm,ntyp,ityp,iproj_group,zeta1
  use m_PseudoPotential,      only : prodphi,ilmt,ltp,mtp,taup,nlmt,ntau,nlmtt
  use m_Crystal_Structure,    only : op,nopr
  use m_Ionic_System,         only : napt
  use m_Charge_Density,       only : hsr
  use m_Electronic_Structure, only : m_ES_sym_comp, m_ES_orbital_den_mat
  use m_Timing,               only : tstatc0_begin, tstatc0_end
  use m_Parallelization,      only : mype, npes, ierr, MPI_CommGroup

! ================================== added by K. Tagami ====x======= 5.0
  use m_Const_Parameters,    only :  YES, NO
  use m_Control_Parameters,  only :  sw_force_simple_mixing_occdiff, &
      &                              sw_recomposing_occmat,initial_occmat, &
      &                              initial_occmat_factor
  use m_IterationNumbers,   only : iteration_electronic, iteration
  use m_Control_Parameters,  only : spin_density_mixfactor
  use m_Control_Parameters,  only : sw_mix_bothspins_sametime_hsr
! ================================================================== 5.0

! =========================== added by K. Tagami =========== 11.0
  use m_Control_Parameters,  only :  noncol, ndim_magmom, sw_mix_charge_hardpart, &
       &                             occmat_diag_mode, ndim_spinor, ndim_chgpot, &
       &                             occmat_file_format
  use m_Electronic_Structure, only : m_ES_orbital_den_mat_noncl
  use m_PseudoPotential,      only : flg_paw
  use m_ES_NonCollinear,     only : m_ES_MagMom_To_DensMat_porb
  use m_Const_Parameters,    only :  CMPLDP
  use m_Charge_Density,       only : hsi
! ========================================================== 11.0

  use m_Crystal_Structure,  only : altv
  use m_Ionic_System,     only : ihubbard, ionic_charge_atomtyp, cps, iatomn, &
       &                         mag_moment0_atoms_is_defined, ionic_charge_atoms, &
       &                         mag_moment0_atomtyp, mag_moment0_atoms, &
       &                         sw_set_initial_magmom_by_atom
  use m_PseudoPotential,  only : radr_paw, nmesh, mmesh, nloc, psirpw
  use m_Const_Parameters,  only : PAI4, SP
  use m_CD_Mag_Moment,  only : rad_cov
  use mpi

  implicit none

  integer, public, allocatable :: i2lp(:) ! d(num_projectors)

! ================================ modified by K. Tagami ================ 11.0
!  integer, private :: max2lp ! max. of i2lp
  integer, public :: max2lp ! max. of i2lp
! ======================================================================= 11.0

  integer, private :: l1max ! max. of l1
  real(kind=DP), public, allocatable, target :: om(:,:,:,:,:) ! d(max2lp,max2lp,max_projs,natm,nspin)
  real(kind=DP), public, allocatable, target :: ommix(:,:,:,:,:) ! d(max2lp,max2lp,max_projs,natm,nspin)
  real(kind=DP), public, allocatable :: omold(:,:,:,:,:) ! d(max2lp,max2lp,max_projs,natm,nspin)
                                                   ! occupation maxtrix, n(m,m)
  integer, private :: nyymax
  logical, private :: occmat_is_read = .false.
  integer, public, allocatable :: nyy(:) ! d(num_projectors)
  integer, public, allocatable :: ilmt_yy(:,:,:) ! d(2,nyymax,num_projectors)
  integer, private, allocatable :: nryy(:,:,:,:) ! d(2,nyymax,num_projectors,nopr+af)
  integer, private, allocatable :: iryy1(:,:,:,:) ! d(max2lp,nyymax,num_projectos,nopr+af)
  integer, private, allocatable :: iryy2(:,:,:,:) ! d(max2lp,nyymax,num_projectors,nopr+af)
  real(kind=DP), private, allocatable :: cryy(:,:,:,:,:) ! d(max2lp,max2lp,nyymax,num_projectors,nopr+af)

! ================================ added by K. Tagami ================ 11.0
  real(kind=DP), public, allocatable, target :: om_aimag(:,:,:,:,:)
  real(kind=DP), public, allocatable, target :: ommix_aimag(:,:,:,:,:)
  real(kind=DP), public, allocatable :: omold_aimag(:,:,:,:,:)
! ===================================================================== 11.0

  integer :: nsize_rho
  integer, private, allocatable :: imapom(:) ! d(nsize_rho)
  real(kind=DP),private,allocatable, dimension(:,:) ::       rho,rhoo ! d(nsize_rho,nspin)
  !         rho => ommix, rhoo => omold

  ! -- For Broyden and DFP mixing method --
  integer, private,parameter                :: iU = 1, iVD = 1, iW = 2, iY = 2, iV = 2
  integer, private                          :: nspin_m
  real(DP),private,allocatable,dimension(:,:) :: d0,u,v,w,dout,dd
  real(DP),private,pointer,dimension(:,:)         :: FF
  real(DP),private,allocatable,target,dimension(:,:) :: din
  !                                             d(nsize_rho,nspin_m)
  real(DP),private,allocatable,dimension(:,:)     :: dF
  !     dF_l(deltaF):= \Delta \cal F^{m} = \cal F^{m} - \cal F^{m-1}
  !              = F[\rho^{m}] - F[\rho^{m-1}] - (\rho^{m} - \rho^{m-1})
  real(DP),private,allocatable,target,dimension(:,:,:,:) :: urec
  real(DP),private,allocatable,dimension(:,:,:)     :: f      !d(nbxmix,nbxmix,nspin)
  real(DP),private,allocatable,dimension(:)         :: g      !d(nbxmix)
  !    f and g are used only when hownew == RENEW
  integer, private,allocatable,dimension(:)         :: ncrspd !d(nbxmix)
  real(DP),private,allocatable,dimension(:,:,:)     :: uuf    !d(nbxmix,nspin,2),
                                                   !only for DFP method
  real(DP),private,allocatable,dimension(:,:)       :: uuf_p
  real(DP),private,allocatable,dimension(:,:)       :: g_p

  integer, private :: previous_waymix = 0

  logical, save :: first = .true.

! ========================================================================= 5.8
  real(kind=DP),private, allocatable, dimension(:):: rmxtrck ! d(nspin_m)
  real(kind=DP), allocatable :: om_store(:,:,:,:,:)
  real(kind=DP), allocatable :: omold_store(:,:,:,:,:)
!
  real(kind=DP), allocatable :: rho_store(:,:)
  real(kind=DP), allocatable :: rhoo_store(:,:)
! ======================================================================== 5.8

! ================================= added by K. Tagami ================ 11.0
  integer :: nsize_rho_realpart
  integer, private :: sw_mix_imaginary_component = ON
! ===================================================================== 11.0

  real(DP),private,allocatable,target,dimension(:,:,:,:) :: urec_l
  real(DP),private,allocatable,dimension(:,:) :: d0_l,u_l,v_l
  real(kind=DP), allocatable, dimension(:,:) :: ynorm
  real(DP),private,allocatable,dimension(:) :: f_p !d(ista_kgpm:iend_kgpm)

!  include 'mpif.h'

contains

  logical function m_OP_occ_mat_is_not_read()
    m_OP_occ_mat_is_not_read = .false.
    if(.not.occmat_is_read) m_OP_occ_mat_is_not_read = .true.
  end function m_OP_occ_mat_is_not_read

  subroutine m_OP_set_i2lp_max2lp
    integer :: it,ip
    integer, parameter :: ntau0=2

! =========================== added by K. Tagami ====================== 11.0
    integer :: nsize
! ===================================================================== 11.0

    do ip=1,num_projectors
       i2lp(ip) = 2*proj_attribute(ip)%l+1
    end do
    max2lp = 0
    do ip=1,num_projectors
       if(i2lp(ip) > max2lp) then
          max2lp = i2lp(ip)
          l1max  = proj_attribute(ip)%l+1
       end if
    end do

! =========================== modified by K. Tagami ====================== 11.0
!!
!!    nyymax = ntau0*l1max**2*(l1max**2+1)/2
!
    nsize = ntau0*( 2*( l1max -1 )+1 )
    nyymax = nsize *( nsize +1 ) /2

! ======================================================================== 11.0

    !!$write(nfout,*) '== m_OP_set_i2lp_max2lp =='
    !!$write(nfout,*) 'ntau=',ntau
    !!$write(nfout,*) 'ntau0=',ntau0
    !!$write(nfout,*) 'i2lp=',i2lp
    !!$write(nfout,*) 'max2lp=',max2lp
    !!$write(nfout,*) 'l1max=',l1max
    !!$write(nfout,*) 'nyymax=',nyymax
  end subroutine m_OP_set_i2lp_max2lp

  subroutine m_OP_store_om()
    omold = om
  end subroutine m_OP_store_om

  subroutine m_OP_om_diff()
    integer :: is,ia,ja,it,ilmt1,ilmt2
    integer :: l1,l2,m1,m2,t1,t2,m1r,m2r
    integer :: iyy,iopr,mm1,mm2
    integer :: ig,ip,i

    integer :: ilp
    real(kind=DP) :: diff,sumdiff
    integer :: icount
    sumdiff = 0.d0
    icount=0
    do is=1,nspin,af+1

       do ia=1,natm
          ig = iproj_group(ia)
          if(ig<1) cycle

          do i=1,num_proj_elems(ig)
             ip=proj_group(i,ig)
             it = proj_attribute(ip)%ityp
             ilp = proj_attribute(ip)%l+1

             do ilmt1 = 1, ilmt(it)
                l1 = ltp(ilmt1,it); m1 = mtp(ilmt1,it)
                t1 = taup(ilmt1,it)
                if ( l1 /= ilp ) cycle

                do ilmt2 = 1, ilmt(it)
                   l2 = ltp(ilmt2,it); m2 = mtp(ilmt2,it)
                   t2 = taup(ilmt2,it)
                   if( l2 /= ilp ) cycle

                   diff = abs(omold(m1,m2,i,ia,is)-om(m1,m2,i,ia,is))
                   sumdiff = sumdiff+diff
                   icount = icount+1
                end do
             end do
! =============================================================
          end do
          !!if(printable) write(nfout,'(a,i5,a,f20.12)') 'omdiff for atom ',ia,' : ',sumdiff/dble(icount)
          sumdiff = 0.d0
          icount = 0
       end do
    end do
  end subroutine m_OP_om_diff

  subroutine m_OP_alloc
    if(.not.allocated(i2lp)) allocate(i2lp(num_projectors))
    call m_OP_set_i2lp_max2lp

! =========================== modified by K. Tagami ====================== 11.0
!    allocate(om(max2lp,max2lp,max_projs,natm,nspin)); om=0.d0
!    allocate(ommix(max2lp,max2lp,max_projs,natm,nspin)); ommix=0.d0
!    allocate(omold(max2lp,max2lp,max_projs,natm,nspin)); omold=0.d0

    if(.not.allocated(om)) then
       allocate(om(max2lp,max2lp,max_projs,natm,ndim_magmom)); om=0.d0
    endif
    if(.not.allocated(ommix))then
       allocate(ommix(max2lp,max2lp,max_projs,natm,ndim_magmom)); ommix=0.d0
    endif
    if(.not.allocated(omold)) then
       allocate(omold(max2lp,max2lp,max_projs,natm,ndim_magmom)); omold=0.d0
    endif
! ========================================================================= 11.0

! ========================= added by K. Tagami =================== 11.0
    if ( noncol ) then
       if(.not.allocated(om_aimag))then
          allocate( om_aimag(max2lp,max2lp,max_projs,natm,ndim_magmom) )
          om_aimag = 0.d0
       endif
       if(.not.allocated(ommix_aimag))then
          allocate( ommix_aimag(max2lp,max2lp,max_projs,natm,ndim_magmom) )
          ommix_aimag = 0.d0
       endif
       if(.not.allocated(omold_aimag))then
          allocate( omold_aimag(max2lp,max2lp,max_projs,natm,ndim_magmom) )
          omold_aimag = 0.d0
       endif
       !om_aimag = 0.d0; ommix_aimag = 0.0d0; omold_aimag = 0.0d0
    endif
! ================================================================ 11.0

    if(.not.allocated(nyy)) allocate(nyy(num_projectors))
    if(.not.allocated(ilmt_yy)) then
       allocate(ilmt_yy(2,nyymax,num_projectors)); ilmt_yy=100
    endif
    !!$write(nfout,*) 'nyymax',nyymax
    !!$write(nfout,*) 'num_projectors',num_projectors
    !!$write(nfout,*) 'allocate(ilmt_yy):',ilmt_yy
  end subroutine m_OP_alloc

  subroutine m_OP_dealloc
    deallocate(i2lp)
    deallocate(om)
    deallocate(ommix)
    deallocate(omold)
    deallocate(nyy)
    deallocate(ilmt_yy)
! ============================ addd by K. Tagami =========== 11.0
    if ( noncol ) then
       deallocate( om_aimag )
       deallocate( ommix_aimag )
       deallocate( omold_aimag )
    endif
! ========================================================== 11.0

! ======= KT_add ======= 13.0AS
    deallocate( nryy )
    deallocate( iryy1, iryy2 )
    deallocate( cryy )
! ====================== 13.0AS

  end subroutine m_OP_dealloc

  subroutine m_OP_ilmt_yy
    integer :: it,iyy,ilmt1,ilmt2,l1,l2,ip,ilp

    !!$write(nfout,*) 'ntyp=',ntyp
    !!$write(nfout,*) 'allocated(nyy)',allocated(nyy)
    !!$write(nfout,*) 'allocated(ilmt_yy)',allocated(ilmt_yy)
    !!$write(nfout,*) 'allocate(ilmt_yy):',ilmt_yy
    do ip=1,num_projectors
       it=proj_attribute(ip)%ityp
       ilp=proj_attribute(ip)%l+1
       if ( it==0 ) cycle
!!$write(nfout,*) 'ip,it,ilp=',ip,it,ilp
       iyy=0
       do ilmt2=1,ilmt(it)
          do ilmt1=1,ilmt2
             !!$write(nfout,*) 'ilmt1,ilmt2=',ilmt1,ilmt2
             l1=ltp(ilmt1,it)
             l2=ltp(ilmt2,it)
             !!$write(nfout,*) 'l1,l2=',l1,l2
             if(l1==ilp .and. l2 == l1) then
                iyy=iyy+1
                !!$write(nfout,*) 'iyy=',iyy
                ilmt_yy(1,iyy,ip) = ilmt1
                !!$write(nfout,*) 'ilmt_yy(1,iyy,it)=',ilmt1
                ilmt_yy(2,iyy,ip) = ilmt2
                !!$write(nfout,*) 'ilmt_yy(2,iyy,it)=',ilmt2
             end if
          end do
       end do
       nyy(ip)=iyy
    end do
  end subroutine m_OP_ilmt_yy

  subroutine m_OP_crotylm(nfout)
    integer, intent(in) :: nfout

    integer :: mmax,nsph
    integer :: lmt,ia,ja,it,mm,isph,isph1,isph2,m1,m2,mm1,mm2,l1
    integer :: iopr,ilmt1,ilmt2
    integer :: iyy,nylm1,nylm2,m1r,m2r
    integer :: ip
    real(kind=DP), allocatable :: crotylm(:,:,:)
    integer, allocatable :: iylm(:,:,:)
    integer, allocatable :: nylm(:,:)
    real(kind=DP) :: opr(3,3)

! debug
!    write(nfout,*) '=== In m_OP_crotylm ==='
! end debug

    mmax  = 2*l1max-1
    nsph = l1max**2
    allocate(crotylm(mmax,nsph,nopr+af))
    allocate(iylm(mmax,nsph,nopr+af))
    allocate(nylm(nsph,nopr+af))

    call get_crotylm(l1max,mmax,nsph,nopr+af,crotylm,iylm,nylm,op)

!debug
!!$    if(printable) then
!!$       do iopr=1,nopr+af
!!$          do isph=1,nsph
!!$             write(nfout,'(2(1x,i3))') iopr,isph
!!$             do mm=1,nylm(isph,iopr)
!!$                write(nfout,'(i3,"=>",f20.8)') iylm(mm,isph,iopr),crotylm(mm,isph,iopr)
!!$             end do
!!$          end do
!!$       end do
!!$    end if
!end debug

    if(.not.allocated(nryy)) allocate(nryy(2,nyymax,num_projectors,nopr+af))
    if(.not.allocated(iryy1)) allocate(iryy1(max2lp,nyymax,num_projectors,nopr+af))
    iryy1=-100
    if(.not.allocated(iryy2)) allocate(iryy2(max2lp,nyymax,num_projectors,nopr+af))
    iryy2=-200
    if(.not.allocated(cryy)) allocate(cryy(max2lp,max2lp,nyymax,num_projectors,nopr+af))
    cryy=-100.d0

    do ip=1,num_projectors
       it = proj_attribute(ip)%ityp
       l1 = proj_attribute(ip)%l+1
       if ( it == 0 ) cycle

       do iyy=1,nyy(ip)
          ilmt1 = ilmt_yy(1,iyy,ip);   ilmt2 = ilmt_yy(2,iyy,ip)
          m1 = mtp(ilmt1,it);          m2 = mtp(ilmt2,it)
          isph1 = get_isph(m1,l1);     isph2 = get_isph(m2,l1)

          do iopr=1,nopr+af
             nylm1 = nylm(isph1,iopr);       nylm2 = nylm(isph2,iopr)
             nryy(1,iyy,ip,iopr) = nylm1;    nryy(2,iyy,ip,iopr) = nylm2
             do mm1=1,nylm1
                iryy1(mm1,iyy,ip,iopr) = get_m(iylm(mm1,isph1,iopr),l1)
             end do
             do mm2=1,nylm2
                iryy2(mm2,iyy,ip,iopr) = get_m(iylm(mm2,isph2,iopr),l1)
             end do
             do mm2=1,nylm2
                do mm1=1,nylm1
                   cryy(mm1,mm2,iyy,ip,iopr) = crotylm(mm1,isph1,iopr) &
                        &                     *crotylm(mm2,isph2,iopr)
               !!!write(6,'("mm1,mm2,iyy,ip,iopr,cryy=",5i3,f10.5)') mm1,mm2,iyy,ip,iopr,cryy(mm1,mm2,iyy,ip,iopr)
                end do
             end do
          end do
       end do
    end do

    deallocate(crotylm);    deallocate(iylm);    deallocate(nylm)

! debug
!!$    if(printable) then
!!$       write(nfout,'("== Summary ==")')
!!$       do iopr=1,nopr+af
!!$          write(nfout,*) 'iopr=',iopr
!!$          do ip=1,num_projectors
!!$             it=proj_attribute(ip)%ityp
!!$             write(nfout,*) 'ip,it,l=',ip,it,proj_attribute(ip)%l
!!$             do iyy=1,nyy(ip)
!!$                ilmt1 = ilmt_yy(1,iyy,ip)
!!$                ilmt2 = ilmt_yy(2,iyy,ip)
!!$                m1 = mtp(ilmt1,it)
!!$                m2 = mtp(ilmt2,it)
!!$                do mm2=1,nryy(2,iyy,ip,iopr)
!!$                   do mm1=1,nryy(1,iyy,ip,iopr)
!!$                      m1r=iryy1(mm1,iyy,ip,iopr)
!!$                      m2r=iryy2(mm2,iyy,ip,iopr)
!!$                      write(nfout,'(3i3,"=>",2i3,f20.8,":",2i3)') &
!!$                  & iyy,m1,m2,m1r,m2r,cryy(mm1,mm2,iyy,ip,iopr),mm1,mm2
!!$                   end do
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$    end if
! end debug

    !! Projector rotations
    mmax  = 2*l1max-1;    nsph = l1max**2
    allocate(crotylm(mmax,nsph,1))
    allocate(iylm(mmax,nsph,1));    allocate(nylm(nsph,1))

    do ip=1,num_projectors
       if(.not.proj_attribute(ip)%frotate) cycle
       call get_opr(proj_attribute(ip)%phi,proj_attribute(ip)%theta,proj_attribute(ip)%psi,opr)
       !debug
       !!$write(nfout,*) 'Rot mat: ip=',ip
       !!$do m1=1,3
       !!$   write(nfout,'(3f10.5)') (opr(m1,m2),m2=1,3)
       !!$end do
       !end debug
       call get_crotylm(l1max,mmax,nsph,1,crotylm,iylm,nylm,opr)
       l1 = proj_attribute(ip)%l+1
       allocate(proj_attribute(ip)%crotylm(2*l1-1,2*l1-1))
       proj_attribute(ip)%crotylm = 0.d0

       do m1=1,2*l1-1
          isph1 = get_isph(m1,l1)
          do mm2=1,nylm(isph1,1)
             isph2 = iylm(mm2,isph1,1)
             m2 = get_m(isph2,l1)
             proj_attribute(ip)%crotylm(m2,m1) = crotylm(mm2,isph1,1)
          end do
       end do
       !debug
       !!$write(nfout,*) 'CrotYlm: ip=',ip,' l=',l1-1
       !!$do m1=1,2*l1-1
       !!$   write(nfout,'(7f10.5)') (proj_attribute(ip)%crotylm(m1,m2),m2=1,2*l1-1)
       !!$end do
       !end debug
    end do
    deallocate(crotylm);    deallocate(iylm);    deallocate(nylm)

  contains

    integer function get_m(isph,l1)
      integer, intent(in) :: isph, l1
      get_m = isph - (l1-1)**2
    end function get_m

    integer function get_isph(m,l1)
      integer, intent(in) :: m, l1
      get_isph = (l1-1)**2 + m
    end function get_isph

    subroutine get_opr(phi,theta,psi,opr)
      real(kind=DP), intent(in) :: phi,theta,psi
      real(kind=DP), intent(out) :: opr(3,3)
      real(kind=DP) :: cphi,sphi,ctheta,stheta,cpsi,spsi
      cphi = cos(phi);         sphi = sin(phi)
      ctheta = cos(theta);     stheta = sin(theta)
      cpsi = cos(psi);         spsi = sin(psi)
      opr(1,1) = cphi*cpsi-sphi*ctheta*sphi
      opr(1,2) = -cphi*spsi-sphi*ctheta*cpsi
      opr(1,3) = sphi*stheta
      opr(2,1) = sphi*cpsi+cphi*ctheta*spsi
      opr(2,2) = -sphi*spsi+cphi*ctheta*cpsi
      opr(2,3) = -cphi*stheta
      opr(3,1) = stheta*spsi
      opr(3,2) = stheta*cpsi
      opr(3,3) = ctheta
    end subroutine get_opr

  end subroutine m_OP_crotylm

  subroutine m_OP_occ_mat_ylm(nfout,pmode)
    integer, intent(in) :: nfout
    integer, intent(in) :: pmode

    integer :: is,ia,ja,it,ilmt1,ilmt2
    integer :: l1,l2,m1,m2,t1,t2,m1r,m2r
    integer :: iyy,iopr,mm1,mm2
    integer :: ig,ip,i

    integer :: ilp
#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    real(kind=DP) :: om2_sum, hsr2_sum
    om2_sum = 0.d0
#endif    

#ifdef __TIMER_SUB__
    call timer_sta(737)
#endif
    !!$write(nfout,*) 'Constracting Occupation matrix...'
!!$    write(nfout,'(" --- << m_OP_occ_mat_ylm >> ---")')
    om = 0.d0
#ifdef __TIMER_DO__
  call timer_sta(861)
#endif
    do is=1,nspin,af+1

       do ia=1,natm
          ig = iproj_group(ia)
          if(ig<1) cycle

          do i=1,num_proj_elems(ig)
             ip=proj_group(i,ig)
             it = proj_attribute(ip)%ityp
             ilp = proj_attribute(ip)%l+1

             !!$write(nfout,*) 'ia,it,ip=',ia,it,ip
             !!$write(nfout,*) 'nyy=',nyy(ip)
! ============================== Modified by K. Tagami ============== 0.1
!             do iyy=1,nyy(ip)
!                !!$write(nfout,*) 'iyy=',iyy
!                ilmt1=ilmt_yy(1,iyy,ip)
!                l1=ltp(ilmt1,it)
!                m1=mtp(ilmt1,it)
!                t1=taup(ilmt1,it)
!
!                ilmt2=ilmt_yy(2,iyy,ip)
!                l2=ltp(ilmt2,it)
!                m2=mtp(ilmt2,it)
!                t2=taup(ilmt2,it)

                !!$write(nfout,*) 'om,hsr,prodphi=' &
                !!$& , om(m1,m2,i,ia,is), hsr(ia,ilmt1,ilmt2,is), prodphi(ip,t1,t2)

!                write(*,*) 'FFF I = ',i
!                if ( l1 /= l2 ) then
!                   write(*,*) 'FFFFFFFFFFFFFF AAAAAAAAA'
!                   stop
!                endif
!                om(m1,m2,i,ia,is)=om(m1,m2,i,ia,is) &
!                & + hsr(ia,ilmt1,ilmt2,is)*prodphi(ip,t1,t2)
! =
!             end do
!             do m2=1,i2lp(ip)
!                do m1=m2+1,i2lp(ip)
!                   !!$write(nfout,*) 'm1,m2=',m1,m2
!                   om(m1,m2,i,ia,is)=om(m2,m1,i,ia,is)
!                end do
!             end do
! ======= === ======== ============ =========== ============= ======== ===
             do ilmt1 = 1, ilmt(it)
                l1 = ltp(ilmt1,it); m1 = mtp(ilmt1,it)
                t1 = taup(ilmt1,it)
                if ( l1 /= ilp ) cycle

                do ilmt2 = 1, ilmt(it)
                   l2 = ltp(ilmt2,it); m2 = mtp(ilmt2,it)
                   t2 = taup(ilmt2,it)
                   if( l2 /= ilp ) cycle

                   if ( ilmt1 < ilmt2 ) then
                      om(m1,m2,i,ia,is)=om(m1,m2,i,ia,is) &
                           & + hsr(ia,ilmt1,ilmt2,is)*prodphi(ip,t1,t2)
                   else
                      om(m1,m2,i,ia,is)=om(m1,m2,i,ia,is) &
                           & + hsr(ia,ilmt2,ilmt1,is)*prodphi(ip,t1,t2)
                   endif
                end do
             end do
! =============================================================
          end do
       end do
    end do
#ifdef __TIMER_DO__
  call timer_end(861)
#endif
    if(nspin==1) om=0.5d0*om
    if(iprihubbard > 2) then
       write(nfout,*) '=== <Unsymmetrized occ. mat.> ==='
       call wd_occ_mat(om)
    end if

#ifdef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    om2_sum = 0.d0
    hsr2_sum = 0.d0
    do is = 1, nspin
       do ia = 1, natm
          ig = iproj_group(ia)
          if(ig<1) cycle
          do i = 1, num_proj_elems(ig)
             do m2 = 1, max2lp
                do m1 = 1, max2lp
                   om2_sum = om2_sum + om(m1,m2,i,ia,is)*om(m1,m2,i,ia,is)
                end do
             end do
          end do
       end do
       do ia = 1, natm
          it = ityp(ia)
          do ilmt1 = 1, ilmt(it)
             do ilmt2 = 1, ilmt(it)
                hsr2_sum = hsr2_sum + hsr(ia,ilmt2,ilmt1,is)*hsr(ia,ilmt2,ilmt1,is)
             end do
          end do
       end do
    end do
    write(nfout,'(" om2_sum = ",f20.10, " hsr2_sum = ",f20.10, " <<m_OP_occ_mat_ylm>>")') om2_sum, hsr2_sum
#endif
    if ( sw_mix_charge_hardpart == OFF ) call symmetrize_occ_mat(om)
!    call symmetrize_occ_mat(om)
    if( pmode==1 .and. iprihubbard > 1) then
       write(nfout,*) '=== <Symmetrized occ. mat.> ==='
       call wd_occ_mat(om)

       call diag_occ_mat
    end if
#ifdef __TIMER_SUB__
    call timer_end(737)
#endif

  end subroutine m_OP_occ_mat_ylm


! ========================= added by K. Tagami ============= 11.0
  subroutine m_OP_occ_mat_ylm_noncl( nfout, pmode )
    integer, intent(in) :: nfout
    integer, intent(in) :: pmode

    integer :: is,ia,ja,it,ilmt1,ilmt2
    integer :: l1,l2,m1,m2,t1,t2,m1r,m2r
    integer :: iyy,iopr,mm1,mm2
    integer :: ig,ip,i

    integer :: ilp

#ifdef __TIMER_SUB__
    call timer_sta(737)
#endif
    om = 0.d0; om_aimag = 0.0d0
#ifdef __TIMER_DO__
  call timer_sta(861)
#endif

    do is=1, ndim_magmom

       do ia=1,natm
          ig = iproj_group(ia)
          if(ig<1) cycle

          do i=1,num_proj_elems(ig)
             ip=proj_group(i,ig)
             it = proj_attribute(ip)%ityp
             ilp = proj_attribute(ip)%l+1

             do ilmt1 = 1, ilmt(it)
                l1 = ltp(ilmt1,it); m1 = mtp(ilmt1,it)
                t1 = taup(ilmt1,it)
                if ( l1 /= ilp ) cycle

                do ilmt2 = 1, ilmt(it)
                   l2 = ltp(ilmt2,it); m2 = mtp(ilmt2,it)
                   t2 = taup(ilmt2,it)
                   if( l2 /= ilp ) cycle

                   om(m1,m2,i,ia,is)=om(m1,m2,i,ia,is) &
                        & + hsr(ia,ilmt1,ilmt2,is)*prodphi(ip,t1,t2)
                   om_aimag(m1,m2,i,ia,is) = om_aimag(m1,m2,i,ia,is) &
                        & + hsi(ia,ilmt1,ilmt2,is)*prodphi(ip,t1,t2)
                end do
             end do

          end do
       end do
    end do
#ifdef __TIMER_DO__
    call timer_end(861)
#endif

!    if(iprihubbard >= 2) then
!       write(nfout,*) '=== <Unsymmetrized occ. mat.> ==='
!       call wd_occ_mat(om)
!    end if

    if ( flg_paw .or. sw_mix_charge_hardpart == ON ) then

    else
       if(iprihubbard > 2) then
          write(nfout,*) '=== <Unsymmetrized occ. mat.> ==='
          call wd_occ_mat_noncl(om,om_aimag)
       end if
       call symmetrize_occ_mat_noncl_B(om,om_aimag)
    endif

    if ( pmode==1 .and. iprihubbard > 1) then
       write(nfout,*) '=== <Symmetrized occ. mat.> ==='
       call wd_occ_mat_noncl(om,om_aimag)

!       if ( occmat_diag_mode == 1 ) then
!          call diag_occ_mat
!       else
          call diag_occ_mat_noncl
!       endif
    end if

#ifdef __TIMER_SUB__
    call timer_end(737)
#endif

  end subroutine m_OP_occ_mat_ylm_noncl
! ================================================================ 11.0

  subroutine m_OP_occ_mat_ao(nfout)
    implicit none
    integer :: nfout
  end subroutine m_OP_occ_mat_ao

! ==================================== Added by K. Tagami =============== 0.2
  subroutine m_OP_occ_mat_ao_kt(nfout, pmode )
    integer, intent(in) :: nfout
    integer, intent(in) :: pmode

    integer :: is,ia,ja,it,ilmt1,ilmt2
    integer :: l1,l2,m1,m2,t1,t2,m1r,m2r
    integer :: iyy,iopr,mm1,mm2
    integer :: ig,ip,i

    integer :: ilp, t0
#ifdef __TIMER_SUB__
    call timer_sta(738)
#endif
!!$    write(nfout,'(" *** << m_OP_occ_mat_ao_kt >> ***")')
    !!$write(nfout,*) 'Constracting Occupation matrix...'
    om = 0.d0
#ifdef __TIMER_DO__
  call timer_sta(862)
#endif
    do is=1,nspin,af+1

       do ia=1,natm
          ig = iproj_group(ia)
          if(ig<1) cycle

          do i=1,num_proj_elems(ig)
             ip=proj_group(i,ig)
             it = proj_attribute(ip)%ityp
             ilp = proj_attribute(ip)%l+1

             t0 = proj_attribute(ip)%t
!
             do ilmt1 = 1, ilmt(it)
                l1 = ltp(ilmt1,it); m1 = mtp(ilmt1,it)
                t1 = taup(ilmt1,it)
                if ( l1 /= ilp ) cycle

                do ilmt2 = 1, ilmt(it)
                   l2 = ltp(ilmt2,it); m2 = mtp(ilmt2,it)
                   t2 = taup(ilmt2,it)
                   if( l2 /= ilp ) cycle

                   if ( ilmt1 < ilmt2 ) then
                      om(m1,m2,i,ia,is)=om(m1,m2,i,ia,is) &
                           & + hsr(ia,ilmt1,ilmt2,is) &
                           &   *prodphi(ip,t1,t0)*prodphi(ip,t0,t2)
                   else
                      om(m1,m2,i,ia,is)=om(m1,m2,i,ia,is) &
                           & + hsr(ia,ilmt2,ilmt1,is) &
                           &   *prodphi(ip,t1,t0)*prodphi(ip,t0,t2)
                   endif
                end do
             end do
!
          end do
       end do
    end do
#ifdef __TIMER_DO__
  call timer_end(862)
#endif
    if(nspin==1) om=0.5d0*om
    if(iprihubbard > 2) then
       write(nfout,*) '=== <Unsymmetrized occ. mat.> ==='
       call wd_occ_mat(om)
    end if

    call symmetrize_occ_mat(om)
    if ( pmode==1 .and. iprihubbard > 1 ) then
       write(nfout,*) '=== <Symmetrized occ. mat.> ==='
       call wd_occ_mat(om)

       call diag_occ_mat
    end if
#ifdef __TIMER_SUB__
    call timer_end(738)
#endif

  end subroutine m_OP_occ_mat_ao_kt
! =============================================================================


! ======================== added by K. Tagami ================== 11.0
  subroutine m_OP_occ_mat_ao_kt_noncl(nfout, pmode )
    integer, intent(in) :: nfout
    integer, intent(in) :: pmode

    integer :: is,ia,ja,it,ilmt1,ilmt2
    integer :: l1,l2,m1,m2,t1,t2,m1r,m2r
    integer :: iyy,iopr,mm1,mm2
    integer :: ig,ip,i

    integer :: ilp, t0
#ifdef __TIMER_SUB__
    call timer_sta(738)
#endif
    !!$write(nfout,*) 'Constracting Occupation matrix...'
    om = 0.d0; om_aimag = 0.0d0
#ifdef __TIMER_DO__
  call timer_sta(862)
#endif

    do is=1, ndim_magmom, af+1

       do ia=1,natm
          ig = iproj_group(ia)
          if(ig<1) cycle

          do i=1,num_proj_elems(ig)
             ip=proj_group(i,ig)
             it = proj_attribute(ip)%ityp
             ilp = proj_attribute(ip)%l+1

             t0 = proj_attribute(ip)%t
!
             do ilmt1 = 1, ilmt(it)
                l1 = ltp(ilmt1,it); m1 = mtp(ilmt1,it)
                t1 = taup(ilmt1,it)
                if ( l1 /= ilp ) cycle

                do ilmt2 = 1, ilmt(it)
                   l2 = ltp(ilmt2,it); m2 = mtp(ilmt2,it)
                   t2 = taup(ilmt2,it)
                   if( l2 /= ilp ) cycle

                   om(m1,m2,i,ia,is) = om(m1,m2,i,ia,is) &
                        & + hsr(ia,ilmt1,ilmt2,is) &
                        &   *prodphi(ip,t1,t0)*prodphi(ip,t0,t2)
                   om_aimag(m1,m2,i,ia,is) = om_aimag(m1,m2,i,ia,is) &
                        & + hsi(ia,ilmt1,ilmt2,is) &
                        &   *prodphi(ip,t1,t0)*prodphi(ip,t0,t2)
                end do
             end do
!
          end do
       end do
    end do
#ifdef __TIMER_DO__
    call timer_end(862)
#endif

!    if(iprihubbard >= 2) then
!       write(nfout,*) '=== <Unsymmetrized occ. mat.> ==='
!       call wd_occ_mat(om)
!    end if

    if ( flg_paw .or. sw_mix_charge_hardpart == ON ) then

    else
       if(iprihubbard > 2) then
          write(nfout,*) '=== <Unsymmetrized occ. mat.> ==='
          call wd_occ_mat_noncl(om,om_aimag)
       end if
       call symmetrize_occ_mat_noncl_B(om,om_aimag)
    endif

    if ( pmode==1 .and. iprihubbard > 1 ) then
       write(nfout,*) '=== <Symmetrized occ. mat.> ==='
       call wd_occ_mat_noncl(om,om_aimag)

!       if ( occmat_diag_mode == 1 ) then
!          call diag_occ_mat
!       else
          call diag_occ_mat_noncl
!       endif
    end if

#ifdef __TIMER_SUB__
    call timer_end(738)
#endif

  end subroutine m_OP_occ_mat_ao_kt_noncl
! ============================================================= 11.0


  subroutine symmetrize_occ_mat(om)
! ===================================== modified by K. Tagami ============== 11.0
!    real(kind=DP), intent(inout) :: om(max2lp,max2lp,max_projs,natm,nspin)
    real(kind=DP), intent(inout) :: om(max2lp,max2lp,max_projs,natm,ndim_magmom)
! ========================================================================== 11.0

    integer :: is,ia,ja,it,ilmt1,ilmt2
    integer :: l1,l2,m1,m2,t1,t2,m1r,m2r
    integer :: iyy,iopr,mm1,mm2
    integer :: ig,ip,i
    real(kind=DP), allocatable :: om_nosym(:,:,:,:,:)

! ======================= modified by K. Tagami ========================= 11.0
!    allocate(om_nosym(max2lp,max2lp,max_projs,natm,nspin))

    allocate(om_nosym(max2lp,max2lp,max_projs,natm,ndim_magmom ))
! ======================================================================= 11.0

    om_nosym=om
    om = 0.d0

! ======================= modified by K. Tagami ========================= 11.0
!    do is=1,nspin,af+1
    do is=1, ndim_magmom, af+1
! ====================================================================== 11.0

       do ia=1,natm
          ig = iproj_group(ia)
          if(ig<1) cycle
          do i=1,num_proj_elems(ig)
             ip=proj_group(i,ig)
             it = proj_attribute(ip)%ityp
             !!$write(nfout,*) 'ia,it,ip=',ia,it,ip
             do iyy=1,nyy(ip)
                ilmt1 = ilmt_yy(1,iyy,ip)
                ilmt2 = ilmt_yy(2,iyy,ip)
                m1 = mtp(ilmt1,it)
                t1=taup(ilmt1,it)
                m2 = mtp(ilmt2,it)
                t2=taup(ilmt2,it)
                if(t1/=1.or.t2/=1) cycle
                do iopr=1,nopr
                   ja=napt(ia,iopr)
                   do mm2=1,nryy(2,iyy,ip,iopr)
                      do mm1=1,nryy(1,iyy,ip,iopr)
                         m1r=iryy1(mm1,iyy,ip,iopr)
                         m2r=iryy2(mm2,iyy,ip,iopr)
                  !!!print *,'mm1,mm2,iyy,ip,iopr,cryy=',mm1,mm2,iyy,ip,iopr,cryy(mm1,mm2,iyy,ip,iopr)
                         om(m1r,m2r,i,ja,is) =  om(m1r,m2r,i,ja,is) + &
                   & om_nosym(m1,m2,i,ia,is)*cryy(mm1,mm2,iyy,ip,iopr)
                         if(m1 /= m2) &
                   & om(m2r,m1r,i,ja,is) =  om(m2r,m1r,i,ja,is) + &
                   & om_nosym(m2,m1,i,ia,is)*cryy(mm1,mm2,iyy,ip,iopr)
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
    om = om/dble(nopr)
    deallocate(om_nosym)

    if(af /= 0) then
       do ia=1,natm
          ig = iproj_group(ia)
          if(ig<1) cycle
          do i=1,num_proj_elems(ig)
             ip=proj_group(i,ig)
             it = proj_attribute(ip)%ityp
             do iyy=1,nyy(ip)
                ilmt1 = ilmt_yy(1,iyy,ip)
                ilmt2 = ilmt_yy(2,iyy,ip)
                m1 = mtp(ilmt1,it)
                t1=taup(ilmt1,it)
                m2 = mtp(ilmt2,it)
                t2=taup(ilmt2,it)
                if(t1/=1.or.t2/=1) cycle
                iopr=nopr+af
                ja=napt(ia,iopr)
                do mm2=1,nryy(2,iyy,ip,iopr)
                   do mm1=1,nryy(1,iyy,ip,iopr)
                      m1r=iryy1(mm1,iyy,ip,iopr)
                      m2r=iryy2(mm2,iyy,ip,iopr)
               !!print *,'mm1,mm2,iyy,ip,iopr,cryy=',mm1,mm2,iyy,ip,iopr,cryy(mm1,mm2,iyy,ip,iopr)
                      om(m1r,m2r,i,ja,2) =  om(m1r,m2r,i,ja,2) + &
                  & om(m1,m2,i,ia,1)*cryy(mm1,mm2,iyy,ip,iopr)
                      if(m1 /= m2) &
                  & om(m2r,m1r,i,ja,2) = om(m2r,m1r,i,ja,2) + &
                  & om(m2,m1,i,ia,1)*cryy(mm1,mm2,iyy,ip,iopr)
                   end do
                end do
             end do
          end do
       end do
    end if

  end subroutine symmetrize_occ_mat


! ======================================= added by K. Tagami =============== 11.0
  subroutine symmetrize_occ_mat_noncl_A(om)
    real(kind=DP), intent(inout) :: om(max2lp,max2lp,max_projs,natm,ndim_magmom)

    integer :: is,ia,ja,it,ilmt1,ilmt2
    integer :: l1,l2,m1,m2,t1,t2,m1r,m2r
    integer :: iyy,iopr,mm1,mm2
    integer :: ig,ip,i
    real(kind=DP), allocatable :: om_nosym(:,:,:,:,:)

    allocate(om_nosym(max2lp,max2lp,max_projs,natm,ndim_magmom ))

    om_nosym=om
    om = 0.d0

    do is=1, ndim_magmom

       do ia=1,natm
          ig = iproj_group(ia)
          if(ig<1) cycle

          do i=1,num_proj_elems(ig)
             ip=proj_group(i,ig)
             it = proj_attribute(ip)%ityp
             !!$write(nfout,*) 'ia,it,ip=',ia,it,ip

             do iyy=1,nyy(ip)
                ilmt1 = ilmt_yy(1,iyy,ip)
                ilmt2 = ilmt_yy(2,iyy,ip)

                m1 = mtp(ilmt1,it);    t1=taup(ilmt1,it)
                m2 = mtp(ilmt2,it);    t2=taup(ilmt2,it)

                if (t1/=1.or.t2/=1) cycle

                do iopr=1,nopr
                   ja = napt(ia,iopr)

                   do mm2 = 1,nryy(2,iyy,ip,iopr)
                      do mm1 = 1,nryy(1,iyy,ip,iopr)
                         m1r = iryy1(mm1,iyy,ip,iopr)
                         m2r = iryy2(mm2,iyy,ip,iopr)

                         om(m1r,m2r,i,ja,is) =  om(m1r,m2r,i,ja,is) &
                              &               + om_nosym(m1,m2,i,ia,is) &
                              &                *cryy(mm1,mm2,iyy,ip,iopr)
                         if ( m1 /= m2 ) then
                            om(m2r,m1r,i,ja,is) = om(m2r,m1r,i,ja,is) &
                                 &              + om_nosym(m2,m1,i,ia,is) &
                                 &               *cryy(mm1,mm2,iyy,ip,iopr)
                         endif

                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

    om = om/dble(nopr)
    deallocate(om_nosym)

  end subroutine symmetrize_occ_mat_noncl_A

  subroutine symmetrize_occ_mat_noncl_B(om, om_aimag)
    real(kind=DP), intent(inout) :: om(max2lp,max2lp,max_projs,natm,ndim_magmom)
    real(kind=DP), intent(inout) :: om_aimag(max2lp,max2lp,max_projs,natm,ndim_magmom)

    integer :: is,ia,ja,it,ilmt1,ilmt2
    integer :: l1,l2,m1,m2,t1,t2,m1r,m2r
    integer :: iyy,iopr,mm1,mm2
    integer :: ig,ip,i

    integer :: ixyz1, ixyz2
    real(kind=DP) :: ctmp1

    real(kind=DP), allocatable :: om_nosym(:,:,:,:,:)
    real(kind=DP), allocatable :: om_nosym_aimag(:,:,:,:,:)

    allocate( om_nosym(max2lp,max2lp,max_projs,natm,ndim_magmom ))
    allocate( om_nosym_aimag(max2lp,max2lp,max_projs,natm,ndim_magmom ))

    om_nosym = om;    om_nosym_aimag = om_aimag
    om = 0.d0;   om_aimag = 0.0d0

    do ia=1,natm
       ig = iproj_group(ia)
       if(ig<1) cycle

       do i=1,num_proj_elems(ig)
          ip=proj_group(i,ig)
          it = proj_attribute(ip)%ityp
             !!$write(nfout,*) 'ia,it,ip=',ia,it,ip

          do iyy=1,nyy(ip)
             ilmt1 = ilmt_yy(1,iyy,ip)
             ilmt2 = ilmt_yy(2,iyy,ip)
             m1 = mtp(ilmt1,it);    t1=taup(ilmt1,it)
             m2 = mtp(ilmt2,it);    t2=taup(ilmt2,it)

             if (t1/=1.or.t2/=1) cycle

             do iopr=1,nopr
                ja = napt(ia,iopr)

                do mm2 = 1,nryy(2,iyy,ip,iopr)
                   do mm1 = 1,nryy(1,iyy,ip,iopr)
                      m1r = iryy1(mm1,iyy,ip,iopr)
                      m2r = iryy2(mm2,iyy,ip,iopr)

                      om(m1r,m2r,i,ja,1) = om(m1r,m2r,i,ja,1) &
                           &              +om_nosym(m1,m2,i,ia,1) &
                           &               *cryy(mm1,mm2,iyy,ip,iopr)
                      om_aimag(m1r,m2r,i,ja,1) = om_aimag(m1r,m2r,i,ja,1) &
                           &              +om_nosym_aimag(m1,m2,i,ia,1) &
                           &               *cryy(mm1,mm2,iyy,ip,iopr)

                      if ( m1 /= m2 ) then
                         om(m2r,m1r,i,ja,1) = om(m2r,m1r,i,ja,1) &
                              &              +om_nosym(m2,m1,i,ia,1) &
                              &               *cryy(mm1,mm2,iyy,ip,iopr)
                         om_aimag(m2r,m1r,i,ja,1) = om_aimag(m2r,m1r,i,ja,1) &
                              &              +om_nosym_aimag(m2,m1,i,ia,1) &
                              &               *cryy(mm1,mm2,iyy,ip,iopr)
                      endif

                      Do ixyz1=1, 3
                         Do ixyz2=1, 3
                            ctmp1 = op(ixyz2, ixyz1, iopr)

                            om(m1r,m2r,i,ja,ixyz2+1) = om(m1r,m2r,i,ja,ixyz2+1) &
                                 &                    + ctmp1 &
                                 &                     *om_nosym(m1,m2,i,ia,ixyz1+1) &
                                 &                     *cryy(mm1,mm2,iyy,ip,iopr)
                            om_aimag(m1r,m2r,i,ja,ixyz2+1) &
                                 &            = om_aimag(m1r,m2r,i,ja,ixyz2+1) &
                                 &              + ctmp1 &
                                 &               *om_nosym_aimag(m1,m2,i,ia,ixyz1+1) &
                                 &               *cryy(mm1,mm2,iyy,ip,iopr)
                            if ( m1 /= m2 ) then
                               om(m2r,m1r,i,ja,ixyz2+1) = om(m2r,m1r,i,ja,ixyz2+1) &
                                    &                    + ctmp1 &
                                    &                     *om_nosym(m2,m1,i,ia,ixyz1+1) &
                                    &                     *cryy(mm1,mm2,iyy,ip,iopr)
                               om_aimag(m2r,m1r,i,ja,ixyz2+1) &
                                    &          = om_aimag(m2r,m1r,i,ja,ixyz2+1) &
                                    &           + ctmp1 &
                                    &             *om_nosym_aimag(m2,m1,i,ia,ixyz1+1) &
                                    &             *cryy(mm1,mm2,iyy,ip,iopr)
                            endif
                         End do
                      End do
                   end do
                end do
             end do
          end do
       end do
    end do

    om = om/dble(nopr)
    om_aimag = om_aimag /dble(nopr)

    deallocate(om_nosym)
    deallocate(om_nosym_aimag)

  end subroutine symmetrize_occ_mat_noncl_B
! ====================================================================== 11.0

  subroutine wd_occ_mat(om)
! =================================== modified by K. Tagami ============= 11.0
!    real(kind=DP), intent(in) :: om(max2lp,max2lp,max_projs,natm,nspin)
    real(kind=DP), intent(in) :: om(max2lp,max2lp,max_projs,natm,ndim_magmom)
! ======================================================================= 11.0

    integer :: l,is,ia,ig,i,ip,it,m1,m2

! ========================= Modified by K. Tagami ============= 11.0
!    do is=1,nspin, af+1
    do is=1, ndim_magmom, af+1
! ============================================================= 11.0

       do ia=1,natm
          ig = iproj_group(ia)
          if(ig<1) cycle
          do i=1,num_proj_elems(ig)
             ip=proj_group(i,ig)
             it = proj_attribute(ip)%ityp
             l = proj_attribute(ip)%l
             write(nfout,'("Occupation Matrix: is,ia,l=",3i5)') is,ia,l
             do m1=1,i2lp(ip)
                write(nfout,'(14f8.3)')(om(m1,m2,i,ia,is),m2=1,i2lp(ip))
             end do
          end do
       end do
    end do
! ========================== Added by K. Tagami ============= 0.2 ==
!    if ( af /=0 ) then
! ============================== Added by K. Tagami =======  0.2b
!      is = 2
! ====================================================
!       do ia=1,natm
!          ig = iproj_group(ia)
!          if(ig<1) cycle
!          do i=1,num_proj_elems(ig)
!             ip=proj_group(i,ig)
!             it = proj_attribute(ip)%ityp
!             l = proj_attribute(ip)%l
!             if(iprihubbard>1) then
!             write(nfout,'("Occupation Matrix: is,ia,l=",3i5)') is,ia,l
!             do m1=1,i2lp(ip)
!                write(nfout,'(14f8.3)')( 0.0d0, m2=1,i2lp(ip))
!             end do
!             endif
!          end do
!       end do
!    endif
! =====================================================================

  end subroutine wd_occ_mat

! ================================== added by K. Tagami ================== 11.0
  subroutine wd_occ_mat_noncl(om, om_aimag)
    real(kind=DP), intent(in) :: om(max2lp,max2lp,max_projs,natm,ndim_magmom)
    real(kind=DP), intent(in) :: om_aimag(max2lp,max2lp,max_projs,natm,ndim_magmom)

    integer :: l,is,ia,ig,i,ip,it,m1,m2

    do is=1, ndim_magmom

       do ia=1,natm
          ig = iproj_group(ia)
          if(ig<1) cycle
          do i=1,num_proj_elems(ig)
             ip=proj_group(i,ig)
             it = proj_attribute(ip)%ityp
             l = proj_attribute(ip)%l
             write(nfout,'("Occupation Matrix: is,ia,l=",3i5)') is,ia,l

             write(nfout,*) ' (Re)'
             do m1=1,i2lp(ip)
                write(nfout,'(14f9.4)')(om(m1,m2,i,ia,is),m2=1,i2lp(ip))
             end do

             write(nfout,*) ' (Im)'
             do m1=1,i2lp(ip)
                write(nfout,'(14f9.4)')(om_aimag(m1,m2,i,ia,is),m2=1,i2lp(ip))
             end do
          end do
       end do
    end do

  end subroutine wd_occ_mat_noncl
! ========================================================================= 11.0

  subroutine m_OP_occ_mat_init(nfout)
    integer, intent(in) :: nfout
    integer :: num_om, i, iproj, ia, is, ip, m1, m2, ig, it, l
    logical :: tf
    real(kind=DP) :: nn,zeta

    call m_Files_open_nfoccmat
    om = 0.d0
! ============================= added by K. Tagami ============== 11.0
    if ( noncol ) om_aimag = 0.0d0
! =============================================================== 11.0
    do is=1, ndim_magmom
       do ia=1,natm
          ig = iproj_group(ia)
          if(ig<1) cycle
          zeta = zeta1(ityp(ia))
          nn   = (1.d0+zeta)
          if(is.eq.2) nn=(1.d0-zeta)
          if(initial_occmat==UNIT_MATRIX) nn=1.d0
          if(initial_occmat==OFF)         nn=0.d0
          nn = nn*initial_occmat_factor
          do i=1,num_proj_elems(ig)
             ip=proj_group(i,ig)
             it = proj_attribute(ip)%ityp
             l = proj_attribute(ip)%l
             do m1=1,i2lp(ip)
                om(m1,m1,i,ia,is) = nn
             end do
          end do
       end do
    end do
  end subroutine m_OP_occ_mat_init

  subroutine m_OP_rd_occ_mat(nfout)
    integer, intent(in) :: nfout
    integer :: num_om, i, iproj, ia, is, ip, m1, m2, ig
    logical :: tf

    call m_Files_open_nfoccmat

    tf = .false.
    if(mype == 0) then
       write(nfout,*) '! Reading occ. mat.'
       rewind nfoccmat
       read(nfoccmat,*,err=101,end=101) num_om
       if(num_om < 1 ) goto 101
       tf = .true.

       do i=1,num_om
          read(nfoccmat,*,err=101,end=101) is,ia,iproj
          ig = iproj_group(ia)
          ip=proj_group(iproj,ig)
          do m1=1,i2lp(ip)
             read(nfoccmat,*,err=101,end=101)(om(m1,m2,iproj,ia,is),m2=1,i2lp(ip))
          end do
       end do

       occmat_is_read = .true.
    end if
101 if(npes > 1) then
       call mpi_bcast(tf,1,mpi_logical,0,MPI_CommGroup,ierr)

! ========================== modified by K. Tagami ====================== 11.0
!       if (tf) call mpi_bcast(om,max2lp*max2lp*max_projs*natm*nspin,mpi_double_precision,0,MPI_CommGroup,ierr)
       if (tf) then
          call mpi_bcast( om, max2lp *max2lp *max_projs *natm *ndim_magmom, &
               &          mpi_double_precision, 0, MPI_CommGroup, ierr )
       endif
       call mpi_bcast(occmat_is_read,1,mpi_logical,0,MPI_CommGroup,ierr)
    end if
! ======================================================================== 11.0

    if(.not.tf) then
       if(mype==0) write(nfout,*) '! no file of occ. mat.; init. om = 0'
       om = 0.d0
    end if

    return

!!$100 continue
!!$    call mpi_bcast(tf,1,mpi_logical,0,MPI_CommGroup,ierr)
!!$    if(mype==0) write(nfout,*) '! no file of occ. mat.; init. om = 0'
!!$    om = 0.d0
!!$    return

  end subroutine m_OP_rd_occ_mat

! ======================================= added by K. Tagami ============= 11.0
  subroutine m_OP_rd_occ_mat_noncl(nfout)
    integer, intent(in) :: nfout
    integer :: num_om, i, iproj, ia, is, ip, m1, m2, ig
    logical :: tf

    om = 0.0d0; om_aimag = 0.0d0

    call m_Files_open_nfoccmat

    tf = .false.
    if (mype == 0) then
       write(nfout,*) '! Reading occ. mat.'
       rewind nfoccmat
       read(nfoccmat,*,err=101,end=101) num_om

       if(num_om < 1 ) goto 101
       tf = .true.

       do i=1,num_om
          read(nfoccmat,*,err=101,end=101) is,ia,iproj
          ig = iproj_group(ia)
          ip=proj_group(iproj,ig)
          do m1=1,i2lp(ip)
             read(nfoccmat,*,err=101,end=101)(om(m1,m2,iproj,ia,is),m2=1,i2lp(ip))
          end do

          if ( occmat_file_format == 1 ) cycle

          do m1=1,i2lp(ip)
             read(nfoccmat,*,err=101,end=101)(om_aimag(m1,m2,iproj,ia,is),m2=1,i2lp(ip))
          end do
       end do
       occmat_is_read = .true.
    end if

101 if(npes > 1) then
       call mpi_bcast(tf,1,mpi_logical,0,MPI_CommGroup,ierr)

       if (tf) then
          call mpi_bcast( om, max2lp *max2lp *max_projs *natm *ndim_magmom, &
               &          mpi_double_precision, 0, MPI_CommGroup, ierr )
          call mpi_bcast( om_aimag, max2lp *max2lp *max_projs *natm *ndim_magmom, &
               &          mpi_double_precision, 0, MPI_CommGroup, ierr )
       endif
       call mpi_bcast(occmat_is_read,1,mpi_logical,0,MPI_CommGroup,ierr)
    end if

    if(.not.tf) then
       if(mype==0) write(nfout,*) '! no file of occ. mat.; init. om = 0'
       om = 0.d0; om_aimag = 0.0d0
    end if

  end subroutine m_OP_rd_occ_mat_noncl
! ====================================================================== 11.0

  subroutine m_OP_wd_occ_mat(nfout)
    integer, intent(in) :: nfout
    integer :: num_om, i, iproj, ia, is, it, ip, ig, m1, m2, l

    call m_Files_open_nfoccmat

    if(mype /= 0) return

    write(nfout,*) '! Writing occ mat'

    num_om=0

! ================================= modified by K. Tagami =========== 11.0
!    do is=1,nspin,af+1
    do is=1, ndim_magmom, af+1
! ================================================================== 11.0

       do ia=1,natm
          ig = iproj_group(ia)
          if(ig<1) cycle
          num_om = num_om + num_proj_elems(ig)
       end do
    end do
    rewind nfoccmat
    write(nfoccmat,'(i5," : num_om")') num_om

! ================================= modified by K. Tagami =========== 11.0
!    do is=1,nspin,af+1
    do is=1, ndim_magmom, af+1
! ================================================================== 11.0

       do ia=1,natm
          ig = iproj_group(ia)
          if(ig<1) cycle
          do i=1,num_proj_elems(ig)
             ip=proj_group(i,ig)
             it = proj_attribute(ip)%ityp
             l = proj_attribute(ip)%l
             !!iproj = proj_attribute(ip)%ielem
             iproj = i
             write(nfoccmat,'(5(1x,i5)," : is, ia, iproj; it, l")') is,ia,iproj,it,l
             do m1=1,i2lp(ip)
                write(nfoccmat,'(14(1x,e20.8))')(om(m1,m2,iproj,ia,is),m2=1,i2lp(ip))
             end do
          end do
       end do
    end do

  end subroutine m_OP_wd_occ_mat

! ======================================= added by K. Tagami =============== 11.0
  subroutine m_OP_wd_occ_mat_noncl(nfout)
    integer, intent(in) :: nfout
    integer :: num_om, i, iproj, ia, is, it, ip, ig, m1, m2, l

    call m_Files_open_nfoccmat
    if(mype /= 0) return

    write(nfout,*) '! Writing occ mat'

    num_om=0

    do is=1, ndim_magmom, af+1
       do ia=1,natm
          ig = iproj_group(ia)
          if(ig<1) cycle
          num_om = num_om + num_proj_elems(ig)
       end do
    end do

    rewind nfoccmat
    write(nfoccmat,'(i5," : num_om")') num_om

    do is=1, ndim_magmom, af+1
       do ia=1,natm
          ig = iproj_group(ia)
          if(ig<1) cycle
          do i=1,num_proj_elems(ig)
             ip=proj_group(i,ig)
             it = proj_attribute(ip)%ityp
             l = proj_attribute(ip)%l
             !!iproj = proj_attribute(ip)%ielem
             iproj = i

             write(nfoccmat,'(5(1x,i5)," : is, ia, iproj; it, l")') is,ia,iproj,it,l
             do m1=1,i2lp(ip)
                write(nfoccmat,'(14(1x,e20.8))') &
                     &       (om(m1,m2,iproj,ia,is),m2=1,i2lp(ip))
             end do
             do m1=1,i2lp(ip)
                write(nfoccmat,'(14(1x,e20.8))') &
                     &       (om_aimag(m1,m2,iproj,ia,is),m2=1,i2lp(ip))
             end do

          end do
       end do
    end do

  end subroutine m_OP_wd_occ_mat_noncl
! ======================================================================= 11.0

  subroutine diag_occ_mat()
    integer :: is,ia,it,ig,i,ip,l,m1,m2,m3
    integer :: lwork,info
    integer, allocatable :: iwork(:),ifail(:)
    real(kind=DP) :: amat(max2lp,max2lp),evec(max2lp),rvmat(max2lp,max2lp)
    real(kind=DP), allocatable :: rwork(:)
    real(kind=DP) :: abstol,nfound
    real(kind=DP), external :: dlamch
    lwork = 8*max2lp
    allocate(rwork(lwork))
    allocate(iwork(5*max2lp))
    allocate(ifail(max2lp))
    abstol = 2*dlamch('S')

! ============================== modified by K. Tagami =========== 11.0
!    do is=1,nspin,af+1
    do is=1,ndim_magmom,af+1
! ================================================================ 11.0

       do ia=1,natm
          ig = iproj_group(ia)
          if(ig<1) cycle
          do i=1,num_proj_elems(ig)
             ip=proj_group(i,ig)
             it = proj_attribute(ip)%ityp
             l = proj_attribute(ip)%l

!!$! ============================== modified by K. Tagami =========== 11.0
!!$!             write(nfout,'("Diagonalizing Occupation Mattrix: is,ia,l=",3i5)') is,ia,l
!!$             write(nfout,'("Diagonalizing Occupation Matrix: is,ia,l=",3i5)') is,ia,l
!!$! ================================================================ 11.0

             amat(:,:) = om(:,:,i,ia,is)
             call dsyevx('V','A','U',i2lp(ip),amat,max2lp,0.d0,0.d0,0,0, &
             & abstol,nfound,evec,rvmat,max2lp,rwork,lwork,iwork,ifail,info)
             if(info /= 0) then
                write(nfout,*) 'dsyevx: info=',info
             end if
!!$             do m1=1,i2lp(ip)
!!$! ================================= Modified by K. Tagami =============== 0.2a
!!$!                write(nfout,'(f8.3,":",14f8.3)') evec(m1),(rvmat(m2,m1),m2=1,i2lp(ip))
!!$                write(nfout,'(f13.7,":",14f8.3)') evec(m1),(rvmat(m2,m1),m2=1,i2lp(ip))
!!$! ======================================================================
!!$             end do
             !!$if(proj_attribute(ip)%strong_correlated) then
             !!$   om(:,:,i,ia,is) = 0.d0
             !!$   do m1=1,i2lp(ip)
             !!$      do m2=1,i2lp(ip)
             !!$         do m3=i2lp(ip)-proj_attribute(ip)%norbital+1,i2lp(ip)
             !!$            om(m1,m2,i,ia,is) = om(m1,m2,i,ia,is) + rvmat(m1,m3)*evec(m3)*rvmat(m2,m3)
             !!$         end do
             !!$      end do
             !!$   end do
             !!$end if
          end do
       end do
    end do
    deallocate(rwork)
    deallocate(iwork)
    deallocate(ifail)
  end subroutine diag_occ_mat

! ===================================== added by K. Tagami ================= 11.0
  subroutine diag_occ_mat_noncl()
    integer :: ia, it, i, ig, ie, ip, l, m1, m2
    integer :: lwork, info
!
    integer :: is1, is2, istmp, size1, size2
!
    real(kind=DP), allocatable :: dmmat_r_magmom( :,:,: )
    real(kind=DP), allocatable :: dmmat_i_magmom( :,:,: )
    complex(kind=CMPLDP), allocatable :: dmmat_ssrep( :,:,: )
!
    real(kind=DP), allocatable :: rwork(:), eigenval(:)
    complex(kind=CMPLDP), allocatable :: amat(:,:), work(:)
!
    do ia=1, natm
       ig = iproj_group(ia)
       if ( ig <1 ) cycle

       do i=1,num_proj_elems(ig)
          ip = proj_group(i,ig)
          it = proj_attribute(ip)%ityp
          ie = proj_attribute(ip)%ielem
          l  = proj_attribute(ip)%l

          size1 = i2lp(ip)
          size2 = size1 *ndim_spinor
          lwork = 2 *size2

          allocate( work(lwork) )
          allocate( amat(size2,size2) ); amat = 0.0d0
          allocate( eigenval(size2) ); eigenval = 0.0d0
          allocate( rwork( 3*size2 ) );

          allocate( dmmat_r_magmom(size1,size1,ndim_magmom) ); dmmat_r_magmom = 0.0d0
          allocate( dmmat_i_magmom(size1,size1,ndim_magmom) ); dmmat_i_magmom = 0.0d0
          allocate( dmmat_ssrep(size1,size1,ndim_chgpot) );  dmmat_ssrep = 0.0d0

          dmmat_r_magmom(:,:,:) = om(:,:,ie,ia,:)
          dmmat_i_magmom(:,:,:) = om_aimag(:,:,ie,ia,:)

          call m_ES_MagMom_To_DensMat_porb( size1**2, dmmat_r_magmom, dmmat_i_magmom, &
               &                            dmmat_ssrep )

          Do is1=1, ndim_spinor
             Do is2=1, ndim_spinor
                Do m1=1, size1
                   Do m2=1, size1
                      istmp = ( is1 -1 )*ndim_spinor +is2
                      amat( size1*(is1-1)+m1, size1*(is2-1)+m2 ) &
                           &           = dmmat_ssrep( m1, m2, istmp )
                   End do
                End Do
             End Do
          End Do

!!$          write(nfout,'("Diagonalizing Occupation Matrix: ia,l=",3i5)') ia, l

          call zheev( 'V', 'U', size2, amat, size2, &
               &       eigenval, work, lwork, rwork, info )

          if (info /= 0) then
             write(nfout,*) 'zheev : info=',info
          end if

          if(iprihubbard>1)then
             write(nfout,*) '(UP spin component)'
             do m1 = 1, size2
                write(nfout,'(f13.7,": (Re)",14f8.3)') eigenval(m1), &
                     &                      ( real(amat(m2,m1)),m2=1,size1 )
                write(nfout,'(15X,"(Im)",14f8.3)') &
                     &                      ( aimag(amat(m2,m1)),m2=1,size1 )
             end do

             write(nfout,*) '(DOWN spin component)'
             do m1 = 1, size2
                write(nfout,'(f13.7,": (Re)",14f8.3)') eigenval(m1), &
                     &                      ( real(amat(m2,m1)),m2=size1+1,size2 )
                write(nfout,'(15X,"(Im)",14f8.3)') &
                     &                      ( aimag(amat(m2,m1)),m2=size1+1,size2 )
             end do
          endif
          deallocate( dmmat_ssrep, dmmat_r_magmom, dmmat_i_magmom )
          deallocate( amat )
          deallocate(rwork,work)
          deallocate(eigenval)
       end do
    end do


  end subroutine diag_occ_mat_noncl
! ========================================================================= 11.0

  subroutine m_OP_mix_om(alpha)
    real(kind=DP), intent(in) :: alpha

! ======================= modified by K. Tagami ================= 11.0
!    ommix = (1.d0-alpha)*ommix + alpha*om
!    if(printable) then
!       write(nfout,*) '=== <Mixed occ. mat. 1> ==='
!       call wd_occ_mat(ommix)
!    end if
!
    if ( noncol ) then
       ommix = (1.d0-alpha)*ommix + alpha*om
       ommix_aimag = (1.d0-alpha)*ommix_aimag + alpha*om_aimag
       if(printable) then
          if(iprihubbard>1) then
             write(nfout,*) '=== <Mixed occ. mat. 1> ==='
             call wd_occ_mat_noncl(ommix,ommix_aimag)
          endif
       end if
    else
       ommix = (1.d0-alpha)*ommix + alpha*om
       if(printable) then
          if(iprihubbard>1) then
             write(nfout,*) '=== <Mixed occ. mat. 1> ==='
             call wd_occ_mat(ommix)
          endif
       end if
    endif
! ================================================================= 11.0
  end subroutine m_OP_mix_om

  function iter_from_reset()
    integer             :: n, nbox
    integer             :: iter_from_reset
    if(hownew ==  ANEW) then
       n = (iteration_for_cmix - istrbr - 1)/(nbxmix-1)
       if(n < 0) n = 0
       nbox = iteration_for_cmix - (n*(nbxmix-1) + istrbr +1) + 2
       iter_from_reset = nbox + istrbr - 1
    else
       iter_from_reset = iteration_for_cmix
    endif
  end function iter_from_reset

  function icrspd_is(iter)
    integer, intent(in) :: iter
    integer             :: icrspd_is
    if(iter-istrbr+1 < nbxmix) then
       icrspd_is = ncrspd(iter-istrbr+1)
    else
       icrspd_is = ncrspd(nbxmix)
    endif
  end function icrspd_is

  subroutine set_ncrspd_mxiter_etc(iter,iuv,mxiter)
    integer, intent(in)  :: iuv,iter
    integer, intent(out) :: mxiter
    if(hownew == RENEW) then
       if((iter-istrbr+1) >= 3) then
          if((iter-istrbr+1) > nbxmix) then   ! When the box overflows
             call rotate_cmix_arrays          !-(contained here) ->mxiter,ncrspd,urec,f,g
          else
             mxiter = (iter-istrbr+1) - 1
             ncrspd(iter-istrbr+1) = iter-istrbr+1
          endif
       else
          mxiter = (iter-istrbr+1) - 1
          ncrspd(1) = 1
          ncrspd(2) = 2
       endif
    else ! if(hownew == ANEW)
       mxiter = (iter-istrbr+1) - 1
       ncrspd(iter-istrbr+1) = iter-istrbr+1
    endif

  contains
    subroutine rotate_cmix_arrays
      integer :: is,j,i,icr,jcr,iwork

! ======================================= modified by K. Tagami ========= 11.0
!      do is = 1, nspin, af+1
      do is = 1, ndim_magmom, af+1
! ======================================================================= 11.0

         do j = 3, nbxmix
            icr = ncrspd(2)
            jcr = ncrspd(j)
            g(j) = f(icr,jcr,is)
            do i = 3, j-1
               icr = ncrspd(i)
               g(j) = g(j) - f(icr,jcr,is)*g(i)
            enddo
            icr = ncrspd(2)
            urec(:,is,jcr,iuv) &
                 & = urec(:,is,jcr,iuv) + g(j)*urec(:,is,icr,iuv)
         enddo
      enddo
      mxiter = nbxmix-1
      iwork = ncrspd(2)
      do i = 2, mxiter
         ncrspd(i)= ncrspd(i+1)
      end do
      ncrspd(mxiter+1) = iwork
    end subroutine rotate_cmix_arrays
  end subroutine set_ncrspd_mxiter_etc

  subroutine mix_dealloc_previous()
    if(previous_waymix == BROYD1) then
       call mix_broyden_deallocate()
    else if(previous_waymix == BROYD2) then
       call mix_broyden_deallocate()
    else if(previous_waymix == DFP) then
       call mix_DFP_deallocate()
    else if(previous_waymix == PULAY) then
       call mix_PULAY_deallocate()
    end if
  end subroutine mix_dealloc_previous

  subroutine mix_broyden_allocate
! ========================= modified by K. Tagami =================== 11.0
!    nspin_m  = nspin/(af+1)
!
    if ( noncol ) then
       nspin_m  = ndim_magmom
    else
       nspin_m  = nspin/(af+1)
    endif
! ====================================================================== 11.0

    allocate(din(nsize_rho,nspin_m))
    allocate(dout(nsize_rho,nspin_m))
    allocate(dF(nsize_rho,nspin_m))
    allocate(urec(nsize_rho,nspin_m,nbxmix,2))
    if(hownew == RENEW) then
! ========================= modified by K. Tagami =================== 11.0
!       allocate(f(nbxmix,nbxmix,nspin))
       allocate(f(nbxmix,nbxmix,ndim_magmom))
! ==================================================================== 11.0

       allocate(g(nbxmix))
    end if
    allocate(ncrspd(nbxmix))
  end subroutine mix_broyden_allocate

  subroutine mix_broyden_deallocate
    if(allocated(din)) deallocate(din)
    if(allocated(dout)) deallocate(dout)
    if(allocated(dF)) deallocate(dF)
    if(allocated(urec)) deallocate(urec)
    if(allocated(f)) deallocate(f)
    if(allocated(g)) deallocate(g)
    if(allocated(ncrspd)) deallocate(ncrspd)
  end subroutine mix_broyden_deallocate

  subroutine mix_broyden_alloc2
    allocate(d0(nsize_rho,nspin_m))
    allocate(u(nsize_rho,nspin_m))
    allocate(v(nsize_rho,nspin_m))
  end subroutine mix_broyden_alloc2

  subroutine mix_broyden_dealloc2
    deallocate(d0); deallocate(u); deallocate(v)
  end subroutine mix_broyden_dealloc2

  subroutine mix_broyden_alloc3
    allocate(d0(nsize_rho,nspin_m))
    allocate(u(nsize_rho,nspin_m))
    allocate(v(nsize_rho,nspin_m))
    allocate(dd(nsize_rho,nspin_m))
  end subroutine mix_broyden_alloc3

  subroutine mix_broyden_dealloc3
    deallocate(d0); deallocate(u); deallocate(v); deallocate(dd)
  end subroutine mix_broyden_dealloc3

  subroutine mix_DFP_allocate
! ========================= modified by K. Tagami =================== 11.0
!    nspin_m  = nspin/(af+1)
!
    if ( noncol ) then
       nspin_m  = ndim_magmom
    else
       nspin_m  = nspin/(af+1)
    endif
! ====================================================================== 11.0

    allocate(din(nsize_rho,nspin_m))
    allocate(dout(nsize_rho,nspin_m))
    allocate(dF(nsize_rho,nspin_m))
    allocate(urec(nsize_rho,nspin_m,nbxmix,2))
    allocate(ncrspd(nbxmix))
    allocate(uuf(nbxmix,nspin_m,2))
  end subroutine mix_DFP_allocate

  subroutine mix_DFP_deallocate()
    if(allocated(din)) deallocate(din)
    if(allocated(dout)) deallocate(dout)
    if(allocated(dF)) deallocate(dF)
    if(allocated(urec)) deallocate(urec)
    if(allocated(ncrspd)) deallocate(ncrspd)
    if(allocated(uuf)) deallocate(uuf)
  end subroutine mix_DFP_deallocate

  subroutine mix_DFP_alloc2
    allocate(d0(nsize_rho,nspin_m))
    allocate(u(nsize_rho,nspin_m))
    allocate(w(nsize_rho,nspin_m))
  end subroutine mix_DFP_alloc2

  subroutine mix_DFP_dealloc2
    deallocate(d0); deallocate(u); deallocate(w)
  end subroutine mix_DFP_dealloc2

  subroutine renew_u_br(j,i)
    integer, intent(in) :: j,i

! ================================== modified by K. Tagami =========== 11.0
!    real(DP)      :: v_dF(nspin)
    real(DP)      :: v_dF(nspin_m)
! ===================================================================== 11.0

    integer       :: is

    v_dF = 0.d0
! ================================ modified by K. Tagami ============ 5.0
!    do is = 1, nspin, af+1
!       v_dF(is) = sum(urec(:,is,j,iV)*dF(:,is))
!       u(:,is) = u(:,is) - v_dF(is)*urec(:,is,j,iU)
!    end do

! ================================== modified by K. Tagami =========== 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ===================================================================== 11.0
       v_dF(is) = sum(urec(:,is,j,iV)*dF(:,is))
    end do

! ================================== modified by K. Tagami =========== 11.0
!    if ( nspin==2 .and. sw_mix_bothspins_sametime_hsr == YES ) then
    if ( nspin_m==2 .and. sw_mix_bothspins_sametime_hsr == YES ) then
! ===================================================================== 11.0
      v_dF(1) = v_dF(1) + v_dF(2)
      v_dF(2) = v_dF(1)
    endif

! ======================== added by K. Tagami ============== 11.0
    if ( noncol ) then
       v_dF(1) = sum( v_dF(:) )
       v_dF(:) = v_dF(1)
    endif
! =========================================================== 11.0

! ================================== modified by K. Tagami =========== 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ===================================================================== 11.0
       u(:,is) = u(:,is) - v_dF(is)*urec(:,is,j,iU)
    end do
! ================================================================== 5.0

! ================================== modified by K. Tagami =========== 11.0
!    if(hownew == RENEW) f(j,i,1:nspin) = v_dF(1:nspin)
    if(hownew == RENEW) f(j,i,1:nspin_m) = v_dF(1:nspin_m)
! ==================================================================== 11.0

  end subroutine renew_u_br

  subroutine renew_d_br(j)
    integer, intent(in) :: j

! ================================== modified by K. Tagami =========== 11.0
!    real(DP)  :: vF(nspin)
    real(DP)  :: vF(nspin_m)
! ==================================================================== 11.0

    integer :: is

    vF = 0.d0
! ================================ modified by K. Tagami ============ 5.0
!    do is = 1, nspin, af+1
!       vF(is) = sum(urec(:,is,j,iV)*FF(:,is))
!       d0(:,is) = d0(:,is) - vF(is)*urec(:,is,j,iU)
!    end do

! ================================== modified by K. Tagami =========== 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ==================================================================== 11.0
       vF(is) = sum(urec(:,is,j,iV)*FF(:,is))
    end do
!
! ================================ modified by K. Tagami ============= 11.0
!    if ( nspin==2 .and. sw_mix_bothspins_sametime_hsr == YES ) then
    if ( nspin_m==2 .and. sw_mix_bothspins_sametime_hsr == YES ) then
! ===================================================================== 11.0
      vF(1) = vF(1) + vF(2)
      vF(2) = vF(1)
    endif
! ==================================================================== 5.0

! ======================== added by K. Tagami ============== 11.0
    if ( noncol ) then
       vF(1) = sum( vF(:) )
       vF(:) = vF(1)
    endif
! =========================================================== 11.0

! ================================ modified by K. Tagami ============= 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ===================================================================== 11.0
       d0(:,is) = d0(:,is) - vF(is)*urec(:,is,j,iU)
    end do
  end subroutine renew_d_br

  subroutine renew_d_last_br
    integer   :: is, i
! ================================ modified by K. Tagami ============= 11.0
!    real(DP)  :: vF(nspin)
    real(DP)  :: vF(nspin_m)
! ===================================================================== 11.0

! ================================ added by K. Tagami ============ 5.0
    integer :: ns
! ================================================================ 5.0

    vF = 0.d0
! ================================ modified by K. Tagami ============= 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ===================================================================== 11.0
       vF(is) = sum(v(:,is)*FF(:,is))
    end do

! ================================ added by K. Tagami ============ 5.0

! ========================================== modified by K. Tagami =========== 11.0
!    if ( nspin==2 .and. sw_mix_bothspins_sametime_hsr == YES ) then
    if ( nspin_m==2 .and. sw_mix_bothspins_sametime_hsr == YES ) then
! ============================================================================ 11.0
      vF(1) = vF(1) + vF(2)
      vF(2) = vF(1)
    endif

! ================================================================ 5.0

! ======================== added by K. Tagami ==================== 11.0
    if ( noncol ) then
       vF(1) = sum( vF(:) )
       vF(:) = vF(1)
    endif
! ================================================================ 11.0

! ================================ modified by K. Tagami ============= 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ===================================================================== 11.0
       din (:,is) = rhoo(:,is) ! chgqo
       dout(:,is) = rho (:,is) ! chgq
    end do

! ============================= modified by K. Tagami =============== 5.0
!    do is = 1, nspin, af+1
!       rho(:,is) = d0(:,is) - vF(is)*u(:,is)
!    end do

! ================================ modified by K. Tagami ============= 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! ==================================================================== 11.0
       rho(:,is) = d0(:,is) - vF(is)*u(:,is)
    end do

! ============================= modified by K. Tagami =============== 11.0
!    if ( sw_force_simple_mixing_occdiff ==ON .and. sw_recomposing_occmat ==ON ) then
!       call simple_mix2_kt( rmxtrck )            ! without mapping to occmix
!    endif

    if ( .not. noncol ) then
       if ( sw_force_simple_mixing_occdiff ==ON .and. sw_recomposing_occmat ==ON ) then
          call simple_mix2_kt( rmxtrck )            ! without mapping to occmix
       endif
    endif
! ======================================================================= 11.0

!!!!!!!!!!!!    call map_rho_to_om(ommix,rho)
! ===================================================================5.0


  end subroutine renew_d_last_br

! <<< Quasi-Newton Methods >>>
!  1. Broyden's 1st method
!  2. Broyden's 2nd method
!  3. DFP method
!
  subroutine m_OP_mix_broyden1(rmx)
    real(kind=DP), intent(in) :: rmx
    integer   :: iter,j,mxiter,icr,jcr
!=========================== modified by K. Tagami ================== 11.0
!    real(DP)  :: vdF(nspin)
    real(DP)  :: vdF(nspin_m)
! ==================================================================== 11.0
    integer   :: is
    integer   :: id_sname = -1
#ifdef __TIMER_SUB__
    call timer_sta(1142)
#endif
    call tstatc0_begin('m_OP_mix_broyden1 ',id_sname)

    if(previous_waymix /= BROYD1) then
       if(first) then
          call create_map_func(.true.)
          call alloc_rho
          call create_map_func(.false.)
          first = .false.
       end if
       call mix_dealloc_previous()
       call mix_broyden_allocate()
       FF => din
    end if

    iter = iter_from_reset()                 !-(m_OP)

    if((iter-istrbr+1) <= 1) then
       call simple_mix(rmx)                  !-(m_OP)
    else
       call map_om_to_rho(om,rho)
       call map_om_to_rho(omold,rhoo)
       call mix_broyden_alloc3   !-(m_OP) d0,u,v, and dd are allocated
       call dF_F_d0_u_v_and_dd   !-(c.h.) dF, FF, dd, initial u,v,d0

       call set_ncrspd_mxiter_etc(iter,iU,mxiter) !-(m_CD) ->mxiter,ncrspd
       icr = icrspd_is(iter)                 !-(m_CD) function
       do j = 2, mxiter
          jcr = ncrspd(j)
          call renew_u_br(jcr,icr) !-(m_OP) |u(m)> = |u(m)> - <v(j)|dF(m)>|u(j)>
          call renew_v(jcr)        !-(c.h.) |v(m)> = |v(m)> - <u(j)|dd(m)>|v(j)>
          call renew_d_br(jcr)     !-(m_OP) |d(m)> = |d(m)> - <v(j)|F(m)> |u(j)>
       enddo!j-loop

       urec(:,:,icr,iU) = u(:,:)
       vdF=0.d0
#ifdef __TIMER_DO__
  call timer_sta(1191)
#endif

! ====================================== modified by K. Tagami ========== 11.0
!       do is=1,nspin,(af+1)
       do is=1,ndim_magmom,(af+1)
! ======================================================================== 11.0
          vdF(is) = sum(v(:,is)*dF(:,is))
          urec(:,is,icr,iV) = v(:,is)/vdF(is)
       end do
#ifdef __TIMER_DO__
  call timer_end(1191)
#endif

       call renew_d_last_br                 !-(m_CD)
                                            ! chgq(|d(m)>) = |d(m)>-<v(m)|F(m)>|u(m)>
       call mix_broyden_dealloc3()          !-(m_CD)
    endif

    previous_waymix = BROYD1
    call tstatc0_end(id_sname)
#ifdef __TIMER_SUB__
    call timer_end(1142)
#endif
  contains
    subroutine dF_F_d0_u_v_and_dd
      !   dF(=deltaF) = (rho - dout) - (rhoo - din)
      !   FF  = rho - rhoo (=\cal F^{m}); u  = (rhoo - din) + c_p*dF;
      !   dd = rhoo - din
      !   d0 = rhoo+c_p* F;               v  = c_p* dd_l

      integer                      :: is,k,i

! ====================================== modified by K. Tagami ========== 11.0
!      do is = 1, nspin, af+1
      do is = 1, ndim_magmom, af+1
! ======================================================================== 11.0
         dF(:,is) = (rho(:,is)-rhoo(:,is)) - (dout(:,is)-din(:,is))
         d0(:,is) = rhoo(:,is)+ rmx*(rho(:,is) - rhoo(:,is))
         dd(:,is) = rhoo(:,is) - din(:,is)
         u(:,is)  = rmx*dF(:,is) + dd(:,is)
         FF(:,is)  = rho(:,is) - rhoo(:,is)
      end do

! ====================================== modified by K. Tagami ========== 11.0
!      do is = 1, nspin, af+1
      do is = 1, ndim_magmom, af+1
! ======================================================================== 11.0
         v(:,is) = rmx*dd(:,is)
      end do

    end subroutine dF_F_d0_u_v_and_dd

    subroutine renew_v(j)
      integer, intent(in) :: j
! ====================================== modified by K. Tagami ========== 11.0
!      real(DP)  :: u_dd(nspin)
      real(DP)  :: u_dd(nspin_m)
! ======================================================================== 11.0
      integer :: is

      u_dd = 0.d0
! ====================================== modified by K. Tagami ========== 11.0
!      do is=1,nspin,(af+1)
      do is=1,ndim_magmom,(af+1)
! ======================================================================= 11.0
         u_dd(is) = sum(urec(:,is,j,iU)*dd(:,is))
         v(:,is)  = v(:,is) - u_dd(is)*urec(:,is,j,iV)
      end do
    end subroutine renew_v
  end subroutine m_OP_mix_broyden1

  subroutine m_OP_mix_broyden2(rmx)
    real(kind=DP), intent(in) :: rmx
    integer   :: iter,j,mxiter,icr,jcr
    logical   :: falloc
    integer   :: id_sname = -1
#ifdef __TIMER_SUB__
    call timer_sta(1143)
#endif
    call tstatc0_begin('m_OP_mix_broyden2 ',id_sname)

    if(previous_waymix /= BROYD2) then
       if(first) then
          call create_map_func(.true.)
          call alloc_rho
          call create_map_func(.false.)
          first = .false.
       end if
       call mix_dealloc_previous()
       call mix_broyden_allocate()
       FF => din
    end if

! ============================ modified by K. Tagami ============== 11.0
!    call map_om_to_rho(om,rho)
!    call map_om_to_rho(omold,rhoo)
!
    if ( noncol ) then
       call map_om_to_rho_noncl(om,   om_aimag,   rho)
       call map_om_to_rho_noncl(omold,omold_aimag,rhoo)
    else
       call map_om_to_rho(om,rho)
       call map_om_to_rho(omold,rhoo)
    endif
! ==================================================================== 11.0

    iter = iter_from_reset()                 !-(m_OP)

! ====================================== added by K. Tagami ========= 5.0

! ================= modified by K. Tagami ====== 11.0
!    nspin_m  = nspin/(af+1)
!    allocate( rmxtrck(nspin_m) )
!    if ( sw_recomposing_occmat == YES .and. af == 0 .and. nspin == 2 ) then
!!!       call alloc_occmatstore_recomp( rmx, rmxtrck )
!       call alloc_rhostore_recomp( rmx, rmxtrck )
!    else
!       rmxtrck = rmx
!    endif

    nspin_m  = ndim_magmom/(af+1)
    allocate( rmxtrck(nspin_m) )
    if ( noncol ) then
       rmxtrck = rmx
    else
       if ( sw_recomposing_occmat == YES .and. af == 0 .and. nspin == 2 ) then
          call alloc_rhostore_recomp( rmx, rmxtrck )
       else
          rmxtrck = rmx
       endif
    end if
! ================================================== 11.0

! =================================================================== 5.0

    if((iter-istrbr+1) <= 1) then
! ===================================== modified by K. Tagami ======== 5.0
!!       call simple_mix(rmx)                  !-(m_OP)
       call simple_mix_kt( rmxtrck )                  !-(m_OP)
! ==================================================================== 5.0
    else
       call mix_broyden_alloc2   !-(m_OP) d0,u, and v are allocated
       call dF_F_d0_u_and_v      !-(c.h.)   dF, FF, initial u,v,d0

       call set_ncrspd_mxiter_etc(iter,iU,mxiter) !-(m_OP) ->mxiter,ncrspd
       !                  when hownew == RENEW: f,g,ncrspd, and urec are reset.

       icr = icrspd_is(iter)                 !-(m_OP) function
       do j = 2, mxiter
          jcr = ncrspd(j)
          call renew_u_br(jcr,icr) !-(m_OP) |u(m)> = |u(m)> - <v(j)|dF(m)>|u(j)>
          call renew_d_br(jcr)     !-(m_OP) |d(m)> = |d(m)> - <v(j)|F(m)> |u(j)>
       enddo!j-loop

       urec(:,:,icr,iU) = u(:,:)  ! storing
       urec(:,:,icr,iV) = v(:,:)  ! storing

       call renew_d_last_br            !-(m_OP) chgq_l(|d(m)>) = |d(m)>-<v(m)|F(m)>|u(m)>
       call mix_broyden_dealloc2       !-(m_OP)
    endif

! ===================================- added by K. Tagami =============== 5.0

! ============== modiifed by K. Tagami ======= 11.0
!    if ( sw_recomposing_occmat == YES .and. af == 0 .and. nspin == 2 ) then
!       call compose_rho_dealloc_store
!    end if

    if ( .not. noncol ) then
       if ( sw_recomposing_occmat == YES .and. af == 0 .and. nspin == 2 ) then
          call compose_rho_dealloc_store
       end if
    endif
! =============================================== 11.0

    deallocate(rmxtrck)

! =========================== modified by K. Tagami =========== 11.0
!    call map_rho_to_om( ommix,rho )
!
    if ( noncol ) then
       call map_rho_to_om_noncl( ommix,ommix_aimag,rho )
    else
       call map_rho_to_om( ommix,rho )
    endif
! =============================================================== 11.0

! =======================================================================5.0

    previous_waymix = BROYD2

    call tstatc0_end(id_sname)
#ifdef __TIMER_SUB__
    call timer_end(1143)
#endif
  contains
    subroutine dF_F_d0_u_and_v
      !   dF(=deltaF) = (rho - dout) - (rhoo - din)
      !   FF = rho - rhoo (=\cal F^{m}); u = (rhoo - din) + a*dF;
      !   d0 = rhoo + a*FF;            v = dF/( |dF| )

      integer                      :: is,i
      real(DP), dimension(nspin_m) :: fff

! ==================================- modified by K. Tagami ============= 5.0
!      do is = 1, nspin, af+1
!         dF(:,is) = (rho(:,is)-rhoo(:,is)) - (dout(:,is)-din(:,is))
!         d0(:,is) = rhoo(:,is) + rmx*(rho(:,is) - rhoo(:,is))
!         u(:,is)  = rmx*dF(:,is) + (rhoo(:,is) - din(:,is))
!         FF(:,is) = rho(:,is) - rhoo(:,is)
!      end do

!      do is = 1, nspin, af+1
!         dF(:,is) = (rho(:,is)-rhoo(:,is)) - (dout(:,is)-FF(:,is))
!         d0(:,is) = rhoo(:,is) + rmx*(rho(:,is) - rhoo(:,is))
!         u(:,is)  = rmx*dF(:,is) + (rhoo(:,is) - FF(:,is))
!         FF(:,is) = rho(:,is) - rhoo(:,is)
!      end do

! ========================== modified by K. Tagami ====== 11.0
!      do is = 1, nspin, af+1
      do is = 1, ndim_magmom, af+1
! ======================================================= 11.0
         dF(:,is) = (rho(:,is)-rhoo(:,is)) - (dout(:,is)-FF(:,is))
         d0(:,is) = rhoo(:,is) + rmxtrck(is) *(rho(:,is) - rhoo(:,is))
         u(:,is)  = rmxtrck(is) *dF(:,is) + (rhoo(:,is) - FF(:,is))
         FF(:,is) = rho(:,is) - rhoo(:,is)
      end do
! ========================================================================= 5.0
! debug
!!$      write(nfout,*) 'dF=',dF
!!$      write(nfout,*) 'd0=',d0
!!$      write(nfout,*) 'u =',u
!!$      write(nfout,*) 'FF=',FF
! end debug

      fff = 0.d0

! ======================= modified by K. Tagami ========= 11.0
!      do is=1,nspin,af+1
      do is=1,ndim_magmom,af+1
! ======================================================= 11.0
         fff(is) = sum(dF(:,is)*dF(:,is))
      end do
      if(sum(fff) < 1.d-40)  call phase_error_with_msg(nfout,' fmult is too small' ,__LINE__,__FILE__)

! ================================= added by K. Tagami ================== 5.0
! === DEBUG by tkato 2011/11/19 ================================================
!     if ( nspin == 2 .and. sw_mix_bothspins_sametime_hsr == YES ) then
      if ( nspin_m == 2 .and. sw_mix_bothspins_sametime_hsr == YES ) then
! ==============================================================================
        fff(1) = fff(1) + fff(2)
        fff(2) = fff(1)
      endif
! ======================================================================= 5.0

! ========================= added by K. Tagami =========================== 11.0
      if ( noncol ) then
         fff(1) = sum( fff(:) )
         fff(:) = fff(1)
      endif
! ======================================================================== 11.0

! ========================================= modified by K. Tagami ========== 11.0
!      do is = 1, nspin, af+1
      do is = 1, ndim_magmom, af+1
! =========================================================================== 11.0
         do i=1,nsize_rho
            v(i,is) = dF(i,is)/fff(is)
         end do
      end do

    end subroutine dF_F_d0_u_and_v
  end subroutine m_OP_mix_broyden2

  subroutine m_OP_mix_DFP(rmx)
    real(kind=DP), intent(in) :: rmx
  end subroutine m_OP_mix_DFP

  subroutine mix_pulay_allocate
! ==================== modified by K. Tagami ======================== 11.0
!    nspin_m  = nspin/(af+1)
    if ( noncol ) then
       nspin_m = ndim_magmom
    else
       nspin_m  = nspin/(af+1)
    endif
! ================================================================= 11.0

! =========================================== Modified by K. Tagami =========
!    allocate(f_p(ista_kgpm:iend_kgpm)); call precon_4_mult(f_p) !-(m_CD)
    allocate(f_p(1:nsize_rho)); f_p = 0
! ============================================================================

    allocate(din(1:nsize_rho,nspin_m))
    allocate(dout(1:nsize_rho,nspin_m))
    allocate(urec_l(1:nsize_rho,nspin_m,nbxmix,2))
    allocate(uuf_p(nbxmix,nspin_m))
    allocate(f(nbxmix,nbxmix,nspin_m))
    allocate(g_p(nbxmix,nspin_m))
    allocate(ncrspd(nbxmix))

    allocate(ynorm(nbxmix,nspin_m));ynorm=1.d0
! ======================================= Added by K. Tagami ===========
    din = 0.0d0; dout = 0.0d0; urec_l = 0.0d0; uuf_p = 0.0d0; f = 0.0d0
    g_p = 0.0d0;  ncrspd = 0
! ======================================================================
  end subroutine mix_pulay_allocate

  subroutine mix_pulay_deallocate
    if(allocated(f_p)) deallocate(f_p)
    if(allocated(din)) deallocate(din)
    if(allocated(dout)) deallocate(dout)
    if(allocated(urec_l)) deallocate(urec_l)
    if(allocated(uuf_p)) deallocate(uuf_p)
    if(allocated(f)) deallocate(f)
    if(allocated(g_p)) deallocate(g_p)
    if(allocated(ncrspd)) deallocate(ncrspd)
    if (allocated(ynorm)) deallocate(ynorm)
  end subroutine mix_pulay_deallocate

  subroutine mix_pulay_alloc2
    allocate(d0_l(nsize_rho,nspin_m))
    d0_l = 0.0d0
  end subroutine mix_pulay_alloc2

  subroutine mix_pulay_dealloc2
    deallocate(d0_l)
  end subroutine mix_pulay_dealloc2

  subroutine m_OP_mix_pulay(rmx)
    integer, parameter  :: iRho = 1, iResid = 2
    real(DP),intent(in) :: rmx
    integer   :: iter, mxiter
    real(DP),pointer,dimension(:)  :: e_wk, f_wk, ww1, finv
    integer, pointer,dimension(:)  :: ip
! --> T. Yamasaki  03 Aug. 2009
    real(kind=DP), allocatable, dimension(:):: rmxtrc ! d(nspin_m)
!   real(kind=DP), allocatable, dimension(:,:,:) :: chgqstore_l, chgqostore_l
! <--
    real(kind=DP) :: rmxtt
    integer   :: id_sname = -1
    call tstatc0_begin('m_OP_mix_pulay ',id_sname,1)

    if(previous_waymix /= PULAY) then
       if(first) then
          call create_map_func(.true.)
          call alloc_rho
          call create_map_func(.false.)
          first = .false.
       end if
       call mix_dealloc_previous()
       call mix_pulay_allocate()
    end if

    if ( noncol ) then
       call map_om_to_rho_noncl(om,   om_aimag,   rho)
       call map_om_to_rho_noncl(omold,omold_aimag,rhoo)
    else
       call map_om_to_rho(om,rho)
       call map_om_to_rho(omold,rhoo)
    endif

    iter = iter_from_reset()                 !-(m_OP)

    nspin_m  = ndim_magmom/(af+1)
    allocate( rmxtrck(nspin_m) )
    if ( noncol ) then
       rmxtrck = rmx
    else
       if ( sw_recomposing_occmat == YES .and. af == 0 .and. nspin == 2 ) then
          call alloc_rhostore_recomp( rmx, rmxtrck )
       else
          rmxtrck = rmx
       endif
    end if
! ========================================================================= 11.0

    if((iter-istrbr+1) <= 1) then
! ===================================== modified by K. Tagami ======== 5.0
!!       call simple_mix(rmx)                  !-(m_OP)
       call simple_mix_kt( rmxtrck )                  !-(m_OP)
! ==================================================================== 5.0
    else
       call mix_pulay_alloc2   !-(m_CD) d0_l,u_l, and w_l are allocated
       call set_ncrspd_mxiter(nbxmix,iter-istrbr,mxiter) ! -> ncrspd, mxiter
!!$       call mix_pulay_alloc3(nbxmix,iter-istrbr)   !-(c.h.) e_wk,f_wk,ww1,finv,ip
       call mix_pulay_alloc3(mxiter)   !-(c.h.) e_wk,f_wk,ww1,finv,ip

       call Resid_and_dd_into_urec(mxiter) !-(c.h.)
       !                               dF ->urec_l; dd ->urec_l; d0_l,din,dout
       call Ri_dot_Rj(mxiter)          !-(c.h.) <R(i)|R(j)>->f
       call get_finv_lapack(nbxmix,mxiter,f)  !-(c.h.) f -> f^{-1}= <R(i)|R(j)>^{-1}

       call Rj_dot_d(mxiter)           !-(c.h.) <R(j)|d>,(j=1,iter-istrb) -> uuf_p

       call get_gmatrix(mxiter)        !-(c.h.) (f,uuf_p)->g
       call renew_d_using_g(mxiter,rmxtrck)     !-(c.h.)

       call mix_pulay_dealloc3                    !-(c.h.)
       call mix_pulay_dealloc2                    !-(m_CD)
    endif

    deallocate(rmxtrck)

    if ( noncol ) then
       call map_rho_to_om_noncl( ommix,ommix_aimag,rho )
    else
       call map_rho_to_om( ommix,rho )
    endif

    previous_waymix = PULAY
    call tstatc0_end(id_sname)
  contains
    subroutine mix_pulay_alloc3(m)
      integer, intent(in) :: m
      allocate(e_wk(m*m)); allocate(f_wk(m*m)); allocate(ww1(m)); allocate(finv(m*m))
      allocate(ip(m))
! ===================================== Added by K. Tagami ============
      e_wk = 0; f_wk = 0; ww1 = 0; finv = 0; ip = 0
! =====================================================================
    end subroutine mix_pulay_alloc3

    subroutine set_ncrspd_mxiter(n,iter,m)
      integer, intent(in)  :: n, iter
      integer, intent(out) :: m
      integer :: i, nx
      if(hownew == ANEW) then
         m = iter
!!$         ncrspd(:) = (/(i,i=1,m)/)
         do i=1,iter
            ncrspd(i) = i
         end do
      else ! hownew == RENEW
         if(iter <= n) then
            m = iter
!!$            ncrspd(:) = (/(i,i=1,m)/)
            do i=1,iter
               ncrspd(i) = i
            end do
         else
            m = n
            nx = ncrspd(1)
            do i = 1, m-1
               ncrspd(i) = ncrspd(i+1)
            end do
            ncrspd(m) = nx
         end if
      end if
    end subroutine set_ncrspd_mxiter

    subroutine mix_pulay_dealloc3
      deallocate(e_wk); deallocate(f_wk); deallocate(ww1); deallocate(finv)
      deallocate(ip)
    end subroutine mix_pulay_dealloc3

    subroutine Resid_and_dd_into_urec(iter)
      integer, intent(in) :: iter
      integer             :: itc,itc0,itc1
      integer :: i,j,k,imix
      real(kind=DP) :: sum1,sum2
      itc = ncrspd(iter)
      urec_l(:,:,itc,iResid) = rho(:,:) - rhoo(:,:) - (dout(:,:) - din(:,:)) ! =dF(=delta F^i)
      urec_l(:,:,itc,iRho  ) = rhoo(:,:) - din(:,:)                ! =dd
      d0_l(:,:) = rho(:,:) - rhoo(:,:)
      din(:,:)  = rhoo(:,:)
      dout(:,:) = rho(:,:)
      ynorm(itc,:)=0.d0
      do i=1,nspin_m
         do k=1,nsize_rho
            ynorm(itc,i) = ynorm(itc,i)+urec_l(k,i,itc,iResid)*urec_l(k,i,itc,iResid)
         enddo
      enddo
      ynorm(itc,:) = 1.d0/sqrt(ynorm(itc,:))
    end subroutine Resid_and_dd_into_urec

    subroutine Ri_dot_Rj(n)
      integer, intent(in) :: n
      integer  :: it,jt,itc,jtc
      real(DP) :: ff1(nspin_m),ff1tmp

      do it = 1, n
         itc = ncrspd(it)
         do jt = it, n
            jtc = ncrspd(jt)
            if(sw_mix_bothspins_sametime==ON .and. nspin_m>1)then
               call mult1s10_reduce_spin(urec_l,nbxmix,2,itc,iResid,urec_l,jtc,iResid,f_p,ff1tmp)   ! <delta F^i|delta F^j>
               ff1(1)=ff1tmp;ff1(2)=ff1tmp
            else
               call mult1s10(urec_l,nbxmix,2,itc,iResid,urec_l,jtc,iResid,f_p,ff1)   ! <delta F^i|delta F^j>
            endif

            if ( noncol ) then
               call mult1s10_reduce_spin( urec_l, nbxmix, 2, itc, iResid, &
                    &                     urec_l, jtc, iResid, f_p, ff1tmp )
                                                        ! <delta F^i|delta F^j>
               ff1(:) = ff1tmp
            endif
            f(it,jt,1:nspin_m) = ff1(1:nspin_m)
            if(jt /= it) f(jt,it,1:nspin_m) = f(it,jt,1:nspin_m)
         end do
      end do
    end subroutine Ri_dot_Rj

    subroutine Rj_dot_d(n)
      integer, intent(in) :: n
      integer  :: jt, jtc
      real(DP) :: ff1(nspin_m),ff1tmp

      do jt = 1, n
         jtc = ncrspd(jt)
         if(sw_mix_bothspins_sametime==ON .and. nspin_m>1)then
            call mult1s5_reduce_spin(urec_l,nbxmix,2,jtc,iResid,d0_l,f_p,ff1tmp)
            ff1(1) = ff1tmp;ff1(2)=ff1tmp
         else
            call mult1s5(urec_l,nbxmix,2,jtc,iResid,d0_l,f_p,ff1)
         endif

         if ( noncol ) then
            call mult1s5_reduce_spin(urec_l,nbxmix,2,jtc,iResid,d0_l,f_p,ff1tmp)
            ff1(:) = ff1tmp
         endif

         uuf_p(jt,1:nspin_m) = ff1(1:nspin_m)
      end do
    end subroutine Rj_dot_d

    subroutine get_finv_lapack(m,n,f)
      integer,intent(in)                             :: m,n
      real(DP),intent(inout),dimension(m,m,nspin_m) :: f
      real(DP), allocatable,dimension(:,:) :: fwork
      integer :: is,inf,it,jt,kt,nnspin
      real(DP) :: div,tmp
      allocate(fwork(n,n))
      nnspin = nspin
      if(sw_mix_bothspins_sametime==ON .or. af==1) nnspin=1

      if ( noncol ) then
         nnspin = 1
      end if

      do is=1,nnspin
         if(ipripulay >= 2) then
            write(nfout,600) n,(('(',it,jt,')',f(it,jt,is),jt=1,n),it=1,n)
600         format(//11x,"**input matrix**"/12x &
                 & ,"horder=",I5/(2x,4(1x,1a,i2,",",i2,1a,e14.6)))
         end if
         fwork=0
         do it=1,n
            do jt=1,n
               fwork(jt,it) = f(jt,it,is)*ynorm(jt,is)*ynorm(it,is)
               if(it==jt) fwork(jt,it)=fwork(jt,it)+alpha_pulay
            enddo
         enddo
         call dpotrf('U',n,fwork,n,inf)
         call dpotri('U',n,fwork,n,inf)
         do it=1,n-1
            do jt=it+1,n
               fwork(jt,it) = fwork(it,jt)
            enddo
         enddo
         do it=1,n
            do jt=1,n
               f(jt,it,is) = fwork(jt,it)*ynorm(jt,is)*ynorm(it,is)
            enddo
         enddo
         if(ipripulay >= 2) then
            write(nfout,630) (('(',it,jt,')',f(it,jt,is),it=1,n),jt=1,n)
630         format(/11x, "**inverse matrix**" &
                 & ,/(2x,4(1x,1a,i2,",",i2,1a,e14.6)))
         end if
      enddo
      if(sw_mix_bothspins_sametime==ON .and. nspin_m>1)then
         do it=1,n
            do jt=1,n
               f(jt,it,2) = f(jt,it,1)
            enddo
         enddo
      endif
! ============================== added by K. Tagami ========== 11.0
      if ( noncol ) then
         do it=1,n
            do jt=1,n
               f(jt,it,:) = f(jt,it,1)
            enddo
         end do
      endif
! ============================================================ 11.0
      deallocate(fwork)

    end subroutine get_finv_lapack

    subroutine get_finv(m,n,f)
      integer,intent(in)                             :: m,n
      real(DP),intent(inout),dimension(m,m,nspin_m) :: f

      integer                        :: icount,is,jt,it,icon
      real(DP)                       :: div

      e_wk = 0.d0
      do it = 1, n
         e_wk(it*it) = 1.d0
      end do

! ======================================= modified by K. Tagami =========== 11.0
!      do is = 1, nspin, af+1
      do is = 1, ndim_magmom, af+1
! ========================================================================== 11.0
         div = 1.d0/f(1,1,is)
         icount = 1
         do jt = 1, n
            do it = 1, n
               f_wk(icount) = f(it,jt,is)*div
               icount = icount + 1
            end do
         end do
         if(ipripulay >= 1) then
            write(nfout,600) n,(('(',it,jt,')',f(it,jt,is)*div,jt=1,n),it=1,n)
600         format(//11x,"**input matrix**"/12x &
                 & ,"horder=",I5/(2x,4(1x,1a,i2,",",i2,1a,e14.6)))
         end if
         call rdecomp(n,f_wk,ww1,ip,icon)
         if(icon /= 0) then
            call phase_error_with_msg(nfout,'LU decomposition is impossible.',__LINE__,__FILE__)
         else
            call rsolve(n,n,f_wk,e_wk,finv,ip)
         endif

         icount = 1
         do jt = 1, n
            do it = 1, n
               f(it,jt,is) = finv(icount)
               icount = icount + 1
            end do
         end do
         if(ipripulay >= 1) then
            write(nfout,630) (('(',it,jt,')',f(it,jt,is),it=1,n),jt=1,n)
630         format(/11x, "**inverse matrix**" &
                 & ,/(2x,4(1x,1a,i2,",",i2,1a,e14.6)))
         end if
      end do
    end subroutine get_finv

    subroutine get_gmatrix(n)
      integer,intent(in) :: n
      integer :: is, it, jt, nnspin
      nnspin = nspin
      if(sw_mix_bothspins_sametime==ON .or. af==1) nnspin=1

! ============================ added by K. Tagami ============= 11.0
      if ( noncol ) nnspin = 1
! ============================================================== 11.0

      g_p = 0.d0
      do is = 1, nnspin
         do it = 1, n
            do jt = 1, n
               g_p(it,is) = g_p(it,is) - f(jt,it,is)*uuf_p(jt,is)
            end do
         end do
         if(ipripulay >= 2) then
            write(nfout,'(" -- g_p(1:",i3,") --")') n
            write(nfout,'(8f20.12)') (g_p(it,is),it=1,n)
         end if
      end do
      if(sw_mix_bothspins_sametime==ON .and. nspin_m>1)then
         do it = 1,n
            g_p(it,2) = g_p(it,1)
         enddo
      endif
! ============================== added by K. Tagami ============ 11.0
      if ( noncol ) then
         do it = 1,n
            g_p(it,:) = g_p(it,1)
         enddo
      endif
! ============================================================== 11.0

    end subroutine get_gmatrix

    subroutine renew_d_using_g(n,p)
      integer, intent(in)                                :: n
      real(DP),intent(in),dimension(nspin_m) :: p
      integer    :: is, k, i, it, itc, ns

!!$      do is = 1, nspin, af+1
      ns = nspin_for_qnewton()
      do is = 1, ns,af+1
         do i = 1,nsize_rho
            rho(i,is)  = rhoo(i,is) + p(is)*d0_l(i,is)
         end do
         do it = 1, n
            itc = ncrspd(it)
            do i = 1,nsize_rho
               rho(i,is) = rho(i,is) + g_p(it,is)* &
                    &        (urec_l(i,is,itc,iRho) + p(is)*urec_l(i,is,itc,iResid))
            end do
         end do
      end do

    end subroutine renew_d_using_g

    integer function nspin_for_qnewton()
      if ( noncol ) then
         nspin_for_qnewton=ndim_magmom
      else
         nspin_for_qnewton=nspin
         if (sw_force_simple_mixing==ON .and. sw_recomposing==ON) nspin_for_qnewton=1
      endif
    end function nspin_for_qnewton

  end subroutine m_OP_mix_pulay

  subroutine create_map_func(paramset)
    logical :: paramset
    integer :: n,ia,ig,i,ip,it,m1,m2
    integer :: ii

    n=0
    do ia=1,natm
       ig = iproj_group(ia)
       if(ig<1) cycle
       do i=1,num_proj_elems(ig)
          ip=proj_group(i,ig)
          it = proj_attribute(ip)%ityp
          do m2=1,i2lp(ip)
             do m1=m2,i2lp(ip)
                n=n+1
                if(.not.paramset) &
                & imapom(n) = m1 + max2lp*(m2-1 + max2lp*( i-1 + max_projs*( ia-1 ) ) )
                !!$if(.not.paramset) write(nfout,*) 'n,m1,m2,i,ia,imapom=',n,m1,m2,i,ia,imapom(n)
             end do
          end do
       end do
    end do
    nsize_rho = n
    !!$if(.not.paramset) then
    !!$   write(nfout,*) 'm_OP: nsize_rho=',nsize_rho
    !!$   write(nfout,*) 'm_OP: imapom=',imapom(1:nsize_rho)
    !!$end if

! =============================== added by K. Tagami ================= 11.0
    if ( noncol ) then
       nsize_rho_realpart = n

      if ( sw_mix_imaginary_component == ON ) then
         do ia=1,natm
            ig = iproj_group(ia)
            if(ig<1) cycle
            do i=1,num_proj_elems(ig)
               ip=proj_group(i,ig)
               it = proj_attribute(ip)%ityp
               do m2=1,i2lp(ip)
                  do m1=m2+1,i2lp(ip)
                     n=n+1
                     if (.not.paramset) &
                          & imapom(n) = m1 + max2lp*(m2-1 &
                          &                + max2lp*( i-1 + max_projs*( ia-1 ) ) )
                  end do
               end do
            end do
         end do
         nsize_rho = n
      endif
    end if
! ===================================================================== 11.0

  end subroutine create_map_func

  subroutine alloc_rho
! ========================= modified by K. Tagami ============= 11.0
!    allocate(rho(nsize_rho,nspin))
!    allocate(rhoo(nsize_rho,nspin))
    allocate(rho(nsize_rho,ndim_magmom))
    allocate(rhoo(nsize_rho,ndim_magmom))
! ============================================================= 11.0
    allocate(imapom(nsize_rho))
  end subroutine alloc_rho

  subroutine dealloc_rho
    deallocate(rho)
    deallocate(rhoo)
    deallocate(imapom)
  end subroutine dealloc_rho

  subroutine map_om_to_rho(om,rho)
    real(kind=DP), intent(in) :: om(max2lp*max2lp*max_projs*natm,nspin)
    real(kind=DP), intent(out) :: rho(nsize_rho,nspin)

    integer :: i,is

    do is=1,nspin,(af+1)
       do i=1,nsize_rho
          rho(i,is) = om(imapom(i),is)
          ! debug
          !  write(nfout,*) 'is,i,imapom,rho=',is,i,imapom(i),rho(i,is)
          ! end debug
       end do
    end do
  end subroutine map_om_to_rho

! ==================================== added by K. Tagami ============== 11.0
  subroutine map_om_to_rho_noncl(om,om_aimag,rho)
    real(kind=DP), intent(in) :: om(max2lp*max2lp*max_projs*natm,ndim_magmom)
    real(kind=DP), intent(in) :: om_aimag(max2lp*max2lp*max_projs*natm,ndim_magmom)
    real(kind=DP), intent(out) :: rho(nsize_rho,ndim_magmom)

    integer :: i,is

    do is=1,ndim_magmom
       do i=1, nsize_rho_realpart
          rho(i,is) = om(imapom(i),is)
       end do
       do i=nsize_rho_realpart+1, nsize_rho
          rho(i,is) = om_aimag(imapom(i),is)
       end do
    end do

  end subroutine map_om_to_rho_noncl
! ============================================================================ 11.0

  subroutine map_rho_to_om(om,rho)
    real(kind=DP), intent(out) :: om(max2lp*max2lp*max_projs*natm,nspin)
    real(kind=DP), intent(in) :: rho(nsize_rho,nspin)
    integer :: i,is,ia,ig,ip,it,m1,m2

    do is=1,nspin,(af+1)
       do i=1,nsize_rho
          om(imapom(i),is) = rho(i,is)
       end do
    end do
    call symmetrize(om)

  contains

    subroutine symmetrize(om)
      real(kind=DP), intent(inout) :: om(max2lp,max2lp,max_projs,natm,nspin)

      do is=1,nspin,(af+1)
         do ia=1,natm
            ig = iproj_group(ia)
            if(ig<1) cycle
            do i=1,num_proj_elems(ig)
               ip=proj_group(i,ig)
               it = proj_attribute(ip)%ityp
               do m2=1,i2lp(ip)
                  do m1=m2,i2lp(ip)
                     if(m1/=m2) om(m2,m1,i,ia,is) = om(m1,m2,i,ia,is)
                  end do
               end do
            end do
         end do
      end do
    end subroutine symmetrize

  end subroutine map_rho_to_om

! ==================================== added by K. Tagami ================ 11.0
  subroutine map_rho_to_om_noncl(om,om_aimag,rho)
    real(kind=DP), intent(out) :: om(max2lp*max2lp*max_projs*natm,ndim_magmom)
    real(kind=DP), intent(out) :: om_aimag(max2lp*max2lp*max_projs*natm,ndim_magmom)
    real(kind=DP), intent(in) :: rho(nsize_rho,ndim_magmom)

    integer :: i,is,ia,ig,ip,it,m1,m2

    om = 0.0d0;
    if ( sw_mix_imaginary_component == ON ) om_aimag = 0.0d0

    do is=1,ndim_magmom
       do i=1,nsize_rho_realpart
          om(imapom(i),is) = rho(i,is)
       end do
       do i=nsize_rho_realpart+1, nsize_rho
          om_aimag(imapom(i),is) = rho(i,is)
       end do
    end do

    call symmetrize(om,om_aimag)

  contains

    subroutine symmetrize(om,om_aimag)

      real(kind=DP), intent(inout) :: om(max2lp,max2lp,max_projs,natm,ndim_magmom)
      real(kind=DP), intent(inout) :: om_aimag(max2lp,max2lp,max_projs,natm,ndim_magmom)

      do is=1,ndim_magmom

         do ia=1,natm
            ig = iproj_group(ia)
            if(ig<1) cycle
            do i=1,num_proj_elems(ig)
               ip=proj_group(i,ig)
               it = proj_attribute(ip)%ityp
               do m2=1,i2lp(ip)
                  do m1=m2,i2lp(ip)
                     if(m1/=m2) then
                        om(m2,m1,i,ia,is) = om(m1,m2,i,ia,is)
                        if ( sw_mix_imaginary_component == ON ) then
                           om_aimag(m2,m1,i,ia,is) = -om_aimag(m1,m2,i,ia,is)
                        endif
                     endif
                  end do
               end do
            end do
         end do
      end do
    end subroutine symmetrize

  end subroutine map_rho_to_om_noncl
! ===================================================================== 11.0

!OCL SERIAL
  subroutine simple_mix(rmx)
    real(kind=DP), intent(in) :: rmx
    din   = rhoo ! chgqo
    dout  = rho  ! chgq
! ================================== modified by K. Tagami ========== 5.0
!!    rho = (rmx-1.d0)*din + rmx*dout
    rho = (1.0d0-rmx)*din + rmx*dout
! ==================================================================== 5.0
    call map_rho_to_om(ommix,rho)
  end subroutine simple_mix

! =================================== added by K. Tagami ============= 5.0
  subroutine simple_mix_kt(rmx_this)
    real(kind=DP), intent(in) :: rmx_this(nspin_m)
! === DEBUG by tkato 2011/11/19 ================================================
    integer :: i
! ==============================================================================

    din   = rhoo ! chgqo
    dout  = rho  ! chgq
! === DEBUG by tkato 2011/11/19 ================================================
!   rho(:,1) = rmx_this(1) *dout(:,1) + ( 1.0D0 - rmx_this(1) )* din(:,1)
!   rho(:,2) = rmx_this(2) *dout(:,2) + ( 1.0D0 - rmx_this(2) )* din(:,2)
    do i = 1, nspin_m
       rho(:,i) = rmx_this(i) *dout(:,i) + ( 1.0D0 - rmx_this(i) )* din(:,i)
    enddo
! ==============================================================================
!!!!!!!!!!!    call map_rho_to_om(ommix,rho)
  end subroutine simple_mix_kt

  subroutine simple_mix2_kt(rmx_this)
    real(kind=DP), intent(in) :: rmx_this(nspin_m)

    din(:,2)   = rhoo(:,2) ! chgqo
    dout(:,2)  = rho(:,2)  ! chgq
    rho(:,2) = rmx_this(2) *dout(:,2) + ( 1.0D0 - rmx_this(2) )* din(:,2)
!!!!!!!!!!    call map_rho_to_om(ommix,rho)
  end subroutine simple_mix2_kt
! ========================================================================= 5.0

  subroutine m_OP_cp_ommix_to_omold
    omold = ommix
! ========================= added by K. Tagami ============== 11.0
    if ( noncol ) omold_aimag = ommix_aimag
! =========================================================== 11.0
  end subroutine m_OP_cp_ommix_to_omold

  subroutine m_OP_cp_ommix_to_om
    om = ommix
! ========================= added by K. Tagami ============== 11.0
    if ( noncol ) om_aimag = ommix_aimag
! =========================================================== 11.0
  end subroutine m_OP_cp_ommix_to_om

  subroutine m_OP_simple_mixing(nfout,rmxt)
    integer, intent(in) :: nfout
    real(kind=DP), intent(in) :: rmxt
#ifdef __TIMER_SUB__
    call timer_sta(1141)
#endif

! =========================== modified by K. Tagami =========== 11.0
!    ommix = (1.d0-rmxt)*ommix + rmxt*om
!    if(printable) then
!       write(nfout,*) '=== <Mixed occ. mat. 2> ==='
!       call wd_occ_mat(ommix)
!    end if
!
    if ( noncol ) then
       ommix = (1.d0-rmxt)*ommix + rmxt*om
       if ( sw_mix_imaginary_component ==  ON ) then
          ommix_aimag = (1.d0-rmxt)*ommix_aimag + rmxt*om_aimag
       endif

       if(printable) then
          if(iprihubbard>1) then
             write(nfout,*) '=== <Mixed occ. mat. 2> ==='
             call wd_occ_mat_noncl(ommix,ommix_aimag)
          endif
       end if
    else
       ommix = (1.d0-rmxt)*ommix + rmxt*om
       if(printable) then
          if(iprihubbard>1) then
             write(nfout,*) '=== <Mixed occ. mat. 2> ==='
             call wd_occ_mat(ommix)
          endif
       end if
    endif
! =============================================================== 11.0

#ifdef __TIMER_SUB__
    call timer_end(1141)
#endif
  end subroutine m_OP_simple_mixing

! ===================== added by K. Tagami ============================= 5.0
  subroutine m_OP_cp_om_to_ommix(nfout,rmxt)
    integer, intent(in) :: nfout
    real(kind=DP), intent(in) :: rmxt

! ================ modifiied by K. Tagami ======== 11.0
!    ommix = om
!    if(printable) then
!       write(nfout,*) '=== <Mixed occ. mat. 2A> ==='
!       call wd_occ_mat(ommix)
!    end if
!
    if ( noncol ) then
       ommix = om;  ommix_aimag = om_aimag
       if(printable) then
          if(iprihubbard>1) then
             write(nfout,*) '=== <Mixed occ. mat. 2A> ==='
             call wd_occ_mat_noncl(ommix,ommix_aimag)
          endif
       end if
    else
       ommix = om
       if(printable) then
          if(iprihubbard>1) then
             write(nfout,*) '=== <Mixed occ. mat. 2A> ==='
             call wd_occ_mat(ommix)
          endif
       end if
    endif
! ================================================ 11.0
  end subroutine m_OP_cp_om_to_ommix

  subroutine compose_rho_dealloc_store
    rho_store = rho

    rho(:,1) = 0.5*( rho_store(:,1) + rho_store(:,2) )
    rho(:,2) = 0.5*( rho_store(:,1) - rho_store(:,2) )

    rhoo = rhoo_store
    deallocate( rho_store, rhoo_store )
  end subroutine compose_rho_dealloc_store

  subroutine alloc_rhostore_recomp( rmxt, rmxtrc )
    real(kind=DP),intent(in) :: rmxt
    real(kind=DP),intent(out),dimension(nspin_m) :: rmxtrc

    allocate( rhoo_store( nsize_rho,nspin) )
    allocate( rho_store( nsize_rho,nspin) )

    rho_store = rho;      rhoo_store = rhoo

     rho(:,1) =  rho_store(:,1) +  rho_store(:,2)
     rho(:,2) =  rho_store(:,2) -  rho_store(:,2)
    rhoo(:,1) = rhoo_store(:,1) + rhoo_store(:,2)
    rhoo(:,2) = rhoo_store(:,2) - rhoo_store(:,2)

    rmxtrc(1) = rmxt;     rmxtrc(2) = rmxt*spin_density_mixfactor

  end subroutine alloc_rhostore_recomp

  subroutine mult1s5(u,mb,muv,j,iuv,v,f_q,fmult)
    integer,intent(in) :: mb,muv,j,iuv
    real(DP),intent(in), dimension(1:nsize_rho,nspin_m,mb,muv) :: u
    real(DP),intent(in), dimension(1:nsize_rho,nspin_m) :: v
    real(DP),intent(in), dimension(1:nsize_rho) :: f_q
    real(DP),intent(out),dimension(nspin_m)            :: fmult

    real(DP) :: p,  fac
    integer  :: is,i

    fmult = 0.d0
    do is = 1, ndim_magmom, af+1
       p = 0.d0
       fac=1.0d0
       do i = 1,nsize_rho
          if ( noncol ) then
             fac=f_q(i)
          else
             if (is==1 .or. sw_recomposing==OFF .or. sw_metric_diff==ON) then
                fac=f_q(i)
             endif
          end if
          p = p + fac*u(i,is,j,iuv)*v(i,is)
       end do
       fmult(is) = p
    enddo
  end subroutine mult1s5

  subroutine mult1s5_reduce_spin(u,mb,muv,j,iuv,v,f_q,fmult)
    integer,intent(in) :: mb,muv,j,iuv
    real(DP),intent(in), dimension(1:nsize_rho,nspin_m,mb,muv) :: u
    real(DP),intent(in), dimension(1:nsize_rho,nspin_m) :: v
    real(DP),intent(in), dimension(1:nsize_rho):: f_q
    real(DP),intent(out)            :: fmult

    real(DP) :: p,  fac
    integer  :: is,i

    fmult = 0.d0
    p = 0.d0

    do is = 1, ndim_magmom, af+1
       fac = 1.0d0
       do i = 1,nsize_rho
          if ( noncol ) then
             fac=f_q(i)
          else
             if (is==1 .or. sw_recomposing==OFF .or. sw_metric_diff==ON) then
                fac=f_q(i)
             endif
          end if
          p = p + fac*u(i,is,j,iuv)*v(i,is)
       end do
    enddo
    fmult = p
  end subroutine mult1s5_reduce_spin

  subroutine mult1s10(u,mb,muv,i,iu,v,j,iv,f_q,fmult)
    integer,intent(in) :: mb,muv,i,iu,j,iv
    real(DP),intent(in), dimension(1:nsize_rho,nspin_m,mb,muv) :: u,v
    real(DP),intent(in), dimension(1:nsize_rho):: f_q
    real(DP),intent(out),dimension(nspin_m)            :: fmult

    real(DP) :: p,  fac
    integer  :: is,ig
    fmult = 0.d0

    do is = 1, ndim_magmom, af+1
       p = 0.d0
       fac = 1.0d0
       do ig = 1,nsize_rho
          if ( noncol ) then
             fac=f_q(ig)
          else
             if (is==1 .or. sw_recomposing==OFF .or. sw_metric_diff==ON) then
                fac=f_q(ig)
             endif
          end if
          p = p + fac*u(ig,is,i,iu)*v(ig,is,j,iv)
       end do
       fmult(is) = p
    enddo
  end subroutine mult1s10

  subroutine mult1s10_reduce_spin(u,mb,muv,i,iu,v,j,iv,f_q,fmult)
    integer,intent(in) :: mb,muv,i,iu,j,iv
    real(DP),intent(in), dimension(1:nsize_rho,nspin_m,mb,muv) :: u,v
    real(DP),intent(in), dimension(1:nsize_rho):: f_q
    real(DP),intent(out)            :: fmult

    real(DP) :: p,  fac
    integer  :: is,ig

    fmult = 0.d0
    p = 0.d0

    do is = 1, ndim_magmom, af+1
       fac = 1.0d0
       do ig = 1,nsize_rho

          if ( noncol ) then
             fac=f_q(ig)
          else
             if (is==1 .or. sw_recomposing==OFF .or. sw_metric_diff==ON) then
                fac=f_q(ig)
             endif
          end if
          p = p + fac*u(ig,is,i,iu)*v(ig,is,j,iv)
       end do
    enddo
    fmult = p
  end subroutine mult1s10_reduce_spin

! ================================================================= 5.0
  subroutine m_OP_occ_mat_isotropic(nfout)
    integer, intent(in) :: nfout
!
    om = 0.0d0
    if ( noncol ) then
       call case_noncol
       call wd_occ_mat_noncl(om,om_aimag)
    else
       call case_col
       call wd_occ_mat(om)
    endif

  contains

    subroutine case_noncol
      integer :: ia, ih, it, ie, l1, immax, im1
      real(kind=DP) :: fac(4), chg1, chg2, cnorm

      Do ia=1, natm
         ih = ihubbard(ia)
         if ( ih == 0 ) cycle

         it = ityp(ia)
         ie = proj_attribute(ih)%ielem
         l1 = proj_attribute(ih)%l
         immax = 2*l1 +1

         chg1 = ( proj_attribute(ih)%initial_spin_sum &
              & +proj_attribute(ih)%initial_spin_diff ) /2.0d0
         chg2 = ( proj_attribute(ih)%initial_spin_sum &
              & -proj_attribute(ih)%initial_spin_diff ) /2.0d0
!
         fac = 0.0d0
         fac(1) = 1.0d0;  fac(ndim_magmom) = 0.0d0

         if ( mag_moment0_atoms_is_defined ) then
            cnorm = sqrt( mag_moment0_atoms(ia,1)**2 &
                  &      +mag_moment0_atoms(ia,2)**2 &
                  &      +mag_moment0_atoms(ia,3)**2 )
            if ( cnorm > 0.0d0 ) then
               fac(2) = mag_moment0_atoms(ia,1) /cnorm
               fac(3) = mag_moment0_atoms(ia,2) /cnorm
               fac(4) = mag_moment0_atoms(ia,3) /cnorm
            endif
         else
            cnorm = sqrt( mag_moment0_atomtyp(it,1)**2 &
                  &      +mag_moment0_atomtyp(it,2)**2 &
                  &      +mag_moment0_atomtyp(it,3)**2 )
            if ( cnorm > 0.0d0 ) then
               fac(2) = mag_moment0_atomtyp(it,1) /cnorm
               fac(3) = mag_moment0_atomtyp(it,2) /cnorm
               fac(4) = mag_moment0_atomtyp(it,3) /cnorm
            endif
         endif

         Do im1=1, immax
            om(im1,im1,ie,ia,1) = chg1 /dble(immax) *fac(1)
            om(im1,im1,ie,ia,2) = chg2 /dble(immax) *fac(2)
            om(im1,im1,ie,ia,3) = chg2 /dble(immax) *fac(3)
            om(im1,im1,ie,ia,4) = chg2 /dble(immax) *fac(4)
         End Do
      end Do
    end subroutine case_noncol

    subroutine case_col
      integer :: ia, ih, it, ie, l1, immax, im1
      real(kind=DP) :: fac, chg1, chg2

      Do ia=1, natm
         ih = ihubbard(ia)
         if ( ih == 0 ) cycle

         it = ityp(ia)
         ie = proj_attribute(ih)%ielem
         l1 = proj_attribute(ih)%l
         immax = 2*l1 +1

         fac = 1.0d0
         if ( mag_moment0_atoms_is_defined ) then
            if ( mag_moment0_atoms(ia,1) < 0.0 ) fac = -1.0d0
         else
            if ( mag_moment0_atomtyp(it,1) < 0.0 ) fac = -1.0d0
         endif

         chg1 = ( proj_attribute(ih)%initial_spin_sum &
              & + fac *proj_attribute(ih)%initial_spin_diff ) /2.0d0
         chg2 = ( proj_attribute(ih)%initial_spin_sum &
              & - fac *proj_attribute(ih)%initial_spin_diff ) /2.0d0

         if ( nspin == 2 ) then
            Do im1=1, immax
               om(im1,im1,ie,ia,1) = chg1 /dble(immax)
               om(im1,im1,ie,ia,2) = chg2 /dble(immax)
            End Do
         else
            call phase_error_with_msg(nfout,"UUU",__LINE__,__FILE__)
         endif
      end Do

    end subroutine case_col

  end subroutine m_OP_occ_mat_isotropic

  subroutine m_OP_occ_mat_gen_kt(nfout)
    integer, intent(in) :: nfout

    integer :: lcmax, msph_max
    integer, allocatable :: isph2(:,:,:), mmt2(:,:), il3(:)

    real(kind=DP), allocatable :: A_lm(:), rr_exp(:,:)
    real(kind=DP), allocatable :: cr2(:,:,:)
    real(kind=DP), allocatable :: elec_negativity(:)
!
    integer :: sw_use_elec_negativity = ON

    lcmax = 4
#ifdef LMAX_SPH_EQ_6
    if ( nloc == 4 ) lcmax = 6
#endif
    msph_max = ( lcmax +1 )**2

#ifdef LMAX_SPH_EQ_6
    allocate( cr2(  16,16,7) ); cr2 = 0.0d0
    allocate( isph2(16,16,7) ); isph2 = 0.0d0
    allocate( mmt2( 16,16)   ); mmt2 = 0
    call sphset_upto_L6( cr2, isph2, mmt2 )
#else
    allocate( cr2(  16,16,6) ); cr2 = 0.0d0
    allocate( isph2(16,16,6) ); isph2 = 0.0d0
    allocate( mmt2( 16,16)   ); mmt2 = 0
    call sphset2(nfout,ipri,lcmax,cr2,isph2,mmt2)
#endif

    allocate( il3(msph_max));call substitute_il3( msph_max, il3 )
    allocate( A_lm(msph_max) ); A_lm = 0.0d0
    allocate( rr_exp(0:lcmax,ntyp) ); rr_exp = 0.0d0
!
    call set_r_exp( rr_exp )

    if ( sw_use_elec_negativity == ON ) then
       allocate( elec_negativity(ntyp) )
       call set_electronegativity( ntyp, iatomn, elec_negativity )
    endif

    if ( noncol ) then
       call set_occmat_noncol
    else
       call set_occmat_col
    endif

    deallocate(il3)
    deallocate(cr2); deallocate(isph2);  deallocate(mmt2)
    deallocate(A_lm); deallocate(rr_exp)
    if ( sw_use_elec_negativity == ON ) deallocate(elec_negativity)

  contains

    subroutine set_occmat_col
      integer :: ia, ih, it, l1, immax, ie, is, n1, im1, im2
      integer :: nocc, count
      real(kind=DP) :: chg1, chg2, fac
      real(kind=DP), allocatable :: mat_wk(:,:), eigval(:), eigvec(:,:)
      real(kind=DP), allocatable :: weight(:)

      Do ia=1, natm
         ih = ihubbard(ia)
         if ( ih == 0 ) cycle

         it = ityp(ia)
         call set_mat_Alm( ia, A_lm )

         l1 = proj_attribute(ih)%l
         immax = 2*l1 +1

         allocate( mat_wk(immax,immax) )
         allocate( eigval(immax) );   allocate( eigvec(immax,immax) )
         allocate( weight(immax) );

         call set_mat_crystal_field_pot( it, l1, immax, mat_wk )
         call diagonalize( immax, mat_wk, eigval, eigvec )
!
         ie = proj_attribute(ih)%ielem
!
         fac = 1.0d0
         if ( mag_moment0_atoms_is_defined ) then
            if ( mag_moment0_atoms(ia,1) < 0.0 ) fac = -1.0d0
         else
            if ( mag_moment0_atomtyp(it,1) < 0.0 ) fac = -1.0d0
         endif
            !
         chg1 = ( proj_attribute(ih)%initial_spin_sum &
              & + fac *proj_attribute(ih)%initial_spin_diff ) /2.0d0
         chg2 = ( proj_attribute(ih)%initial_spin_sum &
              & - fac *proj_attribute(ih)%initial_spin_diff ) /2.0d0

         Do is=1, nspin
            mat_wk = 0.0d0
            if ( is==1 ) then
               nocc = nint( chg1 )
            else
               nocc = nint( chg2 )
            endif

            Do n1=1, immax
               if( n1 <= nocc ) then
                  weight(n1) = 1.0d0
               else
                  weight(n1) = 0.0d0
               endif
            End Do

            Do n1=1, immax
               Do im1=1, immax
                  Do im2=1, immax
                     mat_wk(im1,im2) = mat_wk(im1,im2) +eigvec(im1,n1)*eigvec(im2,n1) &
                          &                            *weight(n1)
                  End do
               End Do
            End do
            om(1:immax,1:immax,ie,ia,is) = mat_wk(1:immax,1:immax)
         End do
         deallocate( mat_wk )
         deallocate( eigvec );    deallocate( eigval )
         deallocate( weight )
      End Do
      call wd_occ_mat(om)

    end subroutine set_occmat_col

    subroutine set_occmat_noncol
      integer :: ia, ih, it, l1, immax, ie, is, n1, im1, im2
      integer :: nocc, count
      real(kind=DP) :: fac(ndim_magmom), cnorm, chg1, chg2
      real(kind=DP), allocatable :: mat_wk(:,:,:), eigval(:), eigvec(:,:)
      real(kind=DP), allocatable :: weight(:)

      Do ia=1, natm
         ih = ihubbard(ia)
         if ( ih == 0 ) cycle

         it = ityp(ia)
         call set_mat_Alm( ia, A_lm )

         l1 = proj_attribute(ih)%l
         immax = 2*l1 +1

         allocate( mat_wk(immax,immax,nspin) )
         allocate( eigval(immax) );   allocate( eigvec(immax,immax) )
         allocate( weight(immax) );

         call set_mat_crystal_field_pot( it, l1, immax, mat_wk )
         call diagonalize( immax, mat_wk, eigval, eigvec )
!
         ie = proj_attribute(ih)%ielem
!
         chg1 = ( proj_attribute(ih)%initial_spin_sum &
              & +proj_attribute(ih)%initial_spin_diff ) /2.0d0
         chg2 = ( proj_attribute(ih)%initial_spin_sum &
              & -proj_attribute(ih)%initial_spin_diff ) /2.0d0

         mat_wk = 0.0d0
         Do is=1, nspin
            if ( is==1 ) then
               nocc = nint( chg1 )
            else
               nocc = nint( chg2 )
            endif

            Do n1=1, immax
               if( n1 <= nocc ) then
                  weight(n1) = 1.0d0
               else
                  weight(n1) = 0.0d0
               endif
            End Do

            Do n1=1, immax
               Do im1=1, immax
                  Do im2=1, immax
                     mat_wk(im1,im2,is) = mat_wk(im1,im2,is) &
                          &              +eigvec(im1,n1)*eigvec(im2,n1) &
                          &                             *weight(n1)
                  End do
               End Do
            End do
         end Do

         fac = 0.0d0
         fac(1) = 1.0d0;  fac(ndim_magmom) = 0.0d0

         if ( mag_moment0_atoms_is_defined ) then
            cnorm = sqrt( mag_moment0_atoms(ia,1)**2 &
                  &      +mag_moment0_atoms(ia,2)**2 &
                  &      +mag_moment0_atoms(ia,3)**2 )
            if ( cnorm > 0.0d0 ) then
               fac(2) = mag_moment0_atoms(ia,1) /cnorm
               fac(3) = mag_moment0_atoms(ia,2) /cnorm
               fac(4) = mag_moment0_atoms(ia,3) /cnorm
            endif
         else
            cnorm = sqrt( mag_moment0_atomtyp(it,1)**2 &
                  &      +mag_moment0_atomtyp(it,2)**2 &
                  &      +mag_moment0_atomtyp(it,3)**2 )
            if ( cnorm > 0.0d0 ) then
               fac(2) = mag_moment0_atomtyp(it,1) /cnorm
               fac(3) = mag_moment0_atomtyp(it,2) /cnorm
               fac(4) = mag_moment0_atomtyp(it,3) /cnorm
            endif
         endif

         om(1:immax,1:immax,ie,ia,1) = ( mat_wk(1:immax,1:immax,1) &
              &                         +mat_wk(1:immax,1:immax,2) )
         Do is=2, ndim_magmom
            om(1:immax,1:immax,ie,ia,is) = ( mat_wk(1:immax,1:immax,1) &
                 &                         -mat_wk(1:immax,1:immax,2) ) *fac(is)
         End Do
         deallocate( mat_wk )
         deallocate( eigvec );    deallocate( eigval )
         deallocate( weight )
      End Do
      call wd_occ_mat_noncl(om,om_aimag)

    end subroutine set_occmat_noncol

    subroutine set_r_exp( rr_exp )        ! <r^l>
      real(kind=DP), intent(out) :: rr_exp(0:lcmax,ntyp)

      integer :: ia, ih, it, il1, l1, ir, ierr
      real(kind=DP) :: rr, c1, c2
      real(kind=DP), allocatable :: wos(:)

      allocate( wos(mmesh) ); wos = 0.0d0

      Do it=1, ntyp
         AtomLoop: Do ia=1, natm
            if ( ityp(ia) /= it ) cycle

            ih = ihubbard(ia)
            if ( ih == 0 ) cycle

            l1 = proj_attribute(ih)%l +1

            call set_weight_exp( ierr, 1, nmesh(it), radr_paw(:,it), wos )
            c1 = 0.0D0
            Do ir=1, nmesh(it)
               c1 = c1 +wos(ir) *psirpw(ir,l1,1,it) **2
            End Do

            if ( c1 > 1.1D0 ) then        ! unbound
            else
               Do il1=0, lcmax
                  c2 = 0.0d0
                  Do ir=1, nmesh(it)
                     rr = radr_paw(ir,it)
                     c2 = c2 +wos(ir) *psirpw(ir,l1,1,it) **2 *rr**(il1)
                  End do
                  rr_exp(il1,it) = c2
               End Do
            endif
            exit AtomLoop
         End do AtomLoop
      End Do
      deallocate( wos )
    end subroutine set_r_exp

    subroutine set_mat_Alm( ia, A_lm )
      integer, intent(in) :: ia
      real(kind=DP) :: A_lm(msph_max)

      integer :: nxmax, nymax, nzmax
      integer :: nx, ny, nz, ja, jt, l1, ilm1, it
      real(kind=DP) :: factor, x1, y1, z1, dist, rcut, q, c1, dist_min

      factor = 1.2
      nxmax = 1;   nymax = 1;     nzmax = 1

      A_lm = 0.0d0

      it = ityp(ia)
      dist_min = 1.0D2

      Do nz=-nzmax, nzmax
         Do ny=-nymax, nymax
            Do nx=-nxmax, nxmax
               Do ja=1, natm
                  jt = ityp(ja)
                  if ( it == jt ) cycle
                  x1 = cps(ja,1) -cps(ia,1) + nx*altv(1,1) +ny*altv(1,2) +nz*altv(1,3)
                  y1 = cps(ja,2) -cps(ia,2) + nx*altv(2,1) +ny*altv(2,2) +nz*altv(2,3)
                  z1 = cps(ja,3) -cps(ia,3) + nx*altv(3,1) +ny*altv(3,2) +nz*altv(3,3)
                  dist = sqrt( x1*x1 +y1*y1 +z1*z1 )
                  dist_min = min( dist_min, dist )
               End Do
            End Do
         End Do
      End Do
      rcut = dist_min *factor

      Do nz=-nzmax, nzmax
         Do ny=-nymax, nymax
            Do nx=-nxmax, nxmax
               Do ja=1, natm
                  jt = ityp(ja)
                  if ( it == jt ) cycle
                  x1 = cps(ja,1) -cps(ia,1) + nx*altv(1,1) +ny*altv(1,2) +nz*altv(1,3)
                  y1 = cps(ja,2) -cps(ia,2) + nx*altv(2,1) +ny*altv(2,2) +nz*altv(2,3)
                  z1 = cps(ja,3) -cps(ia,3) + nx*altv(3,1) +ny*altv(3,2) +nz*altv(3,3)
                  dist = sqrt( x1*x1 +y1*y1 +z1*z1 )

                  if ( dist > 0.1 .and. dist < rcut ) then
                     if ( sw_use_elec_negativity == OFF ) then
                        if ( sw_set_initial_magmom_by_atom == ON ) then
                           if ( mag_moment0_atoms_is_defined ) then
                              q = ionic_charge_atoms(ja)
                           else
                              q = ionic_charge_atomtyp(jt)
                           endif
                        else
                           q = ionic_charge_atomtyp(jt)
                        endif
                     else
                        q = -( elec_negativity(jt) -elec_negativity(it) )
                     endif
                     Do ilm1=1, msph_max
                        l1 = il3(ilm1)
                        call sphr_general( 1, ilm1, x1, y1, z1, c1 )
                        A_lm(ilm1) = A_lm(ilm1) -q *PAI4 /(2*l1 +1) /dist**(l1 +1)*c1
                     End do
                  endif

               End Do
            End do
         End do
      End do
    end subroutine set_mat_Alm

    subroutine set_mat_crystal_field_pot( it, l1, size1, mat )
      integer, intent(in) :: it, l1, size1
      real(kind=DP) :: mat( size1,size1 )

      integer :: nsph1, nsph2, nsph_min, nsph_max
      integer :: ilm3, n, l3
      real(kind=DP) :: c1, gaunt, rad_avg

      nsph_min = l1**2 +1;   nsph_max = ( l1+1 ) **2

      mat = 0.0d0
      Do nsph1=nsph_min, nsph_max
         Do nsph2=nsph_min, nsph_max
            c1 = 0.0d0
            Do n=1, mmt2(nsph1,nsph2)
               ilm3 = isph2(nsph1,nsph2,n); l3=il3(ilm3)
               gaunt = cr2(nsph1,nsph2,n)
               rad_avg = rr_exp( l3, it )
               c1 = c1 +rad_avg *A_lm(ilm3) *gaunt
            End Do
            mat(nsph1 -nsph_min +1,nsph2 -nsph_min +1) = c1
         End Do
      ENd Do
    end subroutine set_mat_crystal_field_pot

    subroutine diagonalize( immax, mat, eigval, eigvec )
      integer, intent(in) :: immax
      real(kind=DP), intent(in) :: mat( immax, immax )
      real(kind=DP), intent(out) :: eigval(immax), eigvec(immax,immax)

      integer :: lwork,info
      real(kind=DP) :: abstol,nfound, ctmp
      integer, allocatable :: iwork(:),ifail(:)
      real(kind=DP), allocatable :: rwork(:)
      real(kind=DP), external :: dlamch

      abstol = 2*dlamch('S')
      lwork = 8 *immax
      allocate( rwork(lwork) );
      allocate( iwork(5*immax) );   allocate( ifail(immax) )

      call dsyevx( 'V', 'A', 'U', immax, mat, immax, 0.d0, 0.d0, 0, 0, &
           &       abstol, nfound, eigval, eigvec, immax, rwork, lwork, &
           &       iwork, ifail, info )
      !
      if (info /= 0) then
         write(nfout,*) 'dsyevx: info=',info
      end if
      deallocate( rwork );   deallocate( iwork );   deallocate( ifail )
    end subroutine diagonalize

  end subroutine m_OP_occ_mat_gen_kt

! ********************
! https://en.wikipedia.org/wiki/Electronegativity
! ********************
  subroutine set_electronegativity( ntyp, iatomn, elec_negativity )
    integer, intent(in) :: ntyp
    real(kind=DP), intent(in) :: iatomn(ntyp)
    real(kind=DP), intent(inout) :: elec_negativity(ntyp)
    !
    integer :: nmax_elem = 120
    !
    integer :: it, i
    real(kind=DP), allocatable :: elec_negativity_table(:)

    allocate( elec_negativity_table( nmax_elem ) )
    elec_negativity_table = 0.0d0

    call set_table(  1, 2.20, 0.00, 0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, 0.00 )
    call set_table( 11, 0.93, 1.31, 1.61, 1.90, 2.19, 2.58, 3.16, 0.00, 0.82, 1.00 )
    call set_table( 21, 1.36, 1.54, 1.63, 1.66, 1.55, 1.83, 1.88, 1.91, 1.90, 1.65 )
    call set_table( 31, 1.81, 2.01, 2.18, 2.55, 2.96, 3.00, 0.82, 0.95, 1.22, 1.33 )
    call set_table( 41, 1.60, 2.16, 1.90, 2.20, 2.28, 2.20, 1.93, 1.69, 1.78, 1.96 )
    call set_table( 51, 2.05, 2.10, 2.66, 2.60, 0.79, 0.89, 1.10, 1.12, 1.13, 1.14 )
    call set_table( 61, 1.13, 1.17, 1.20, 1.20, 1.20, 1.22, 1.23, 1.24, 1.25, 1.10 )
    call set_table( 71, 1.27, 1.30, 1.50, 2.36, 1.90, 2.20, 2.20, 2.28, 2.54, 2.00 )
    call set_table( 81, 1.62, 2.33, 2.02, 2.00, 2.20, 0.00, 0.70, 0.89, 1.10, 1.30 )
    call set_table( 91, 1.50, 1.38, 1.36, 1.28, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30 )

    Do it=1, ntyp
       i = nint( iatomn(it) )
       elec_negativity(it) = elec_negativity_table(i)
    End do
    deallocate( elec_negativity_table )

  contains

    subroutine set_table( istart, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10 )
      integer, intent(in) :: istart
      real(kind=SP), intent(in) :: c1, c2, c3, c4, c5, c6, c7, c8, c9, c10

      elec_negativity_table( istart   ) = dble(c1)
      elec_negativity_table( istart+1 ) = dble(c2)
      elec_negativity_table( istart+2 ) = dble(c3)
      elec_negativity_table( istart+3 ) = dble(c4)
      elec_negativity_table( istart+4 ) = dble(c5)
      elec_negativity_table( istart+5 ) = dble(c6)
      elec_negativity_table( istart+6 ) = dble(c7)
      elec_negativity_table( istart+7 ) = dble(c8)
      elec_negativity_table( istart+8 ) = dble(c9)
      elec_negativity_table( istart+9 ) = dble(c10)
    end subroutine set_table

  end subroutine set_electronegativity

end module m_Orbital_Population

