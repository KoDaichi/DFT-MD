!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 581 $)
!
!  MODULE: m_BP_Properties
!
!  AUTHOR(S): T. Yamamoto   Oct/01/2004
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
! ==============
!  patch 0.1 by K. Tagami @adv    2008/10/19
!
!     patch 0.1 :  correction in the case when nspin==2
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
module m_BerryPhase
  ! $Id: m_BerryPhase.F90 581 2018-08-01 08:38:42Z jkoga $
  use m_Const_Parameters,    only : DP, PAI, PAI2, PAI4, ON, OFF
  use m_Control_Parameters,  only : ipriberry, nspin, printable, ndim_spinor, noncol
  use m_Files,               only : nfout
  use m_Parallelization,    only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype &
       &              , nrank_e,myrank_e,map_e,ista_e,iend_e,istep_e &
       &              , map_z,np_e,mpi_k_world,myrank_k,map_k,ista_k,iend_k &
       &              , map_ek, myrank_ke &
       &              , ista_g1k, iend_g1k, mpi_ke_world, mpi_chg_world &
       &              , nel_fs, nis_fs, mpi_ke_world, ista_atm, iend_atm, mpi_g_world &
#ifdef TRANSPOSE
       &              , ierr,mp_e,nis_e,nie_e,nel_e &
       &              , ista_g1,iend_g1,np_g1,mp_g1,nis_g1,nie_g1,nel_g1
#else
       &              , ierr
#endif
  use mpi

  implicit none
!  include 'mpif.h'
  integer, private :: n1,n2,jpara
  integer, private :: ig
  integer, private :: nkvec,nkprep
  integer, private :: nstate = -1
  integer, private :: ikp = 0
  integer, private :: sw_symmetry = 2 ! 1->idetical, 2->inversion, 3->full
  integer, private :: sw_deficit = ON ! 0->off, 1->on
  real(kind=DP), private :: g(3)
  real(kind=DP), allocatable, private :: kvec(:,:) ! dim(3,nkvec)
  real(kind=DP), allocatable, private :: wgh(:) ! dim(nkprep)
  real(kind=DP), allocatable, private :: wfbp0(:,:,:,:) ! dim(kg,np_e,2,nspin)
  real(kind=DP), allocatable, private :: wfbp1(:,:,:,:) ! dim(kg,np_e,2,nspin)
  real(kind=DP), allocatable, private :: wfbp2(:,:,:,:) ! dim(kg,np_e,2,nspin)
  integer, private :: kgbp0,kgbp1,kgbp2
  complex(kind=DP), allocatable, private :: det(:,:,:) ! dim(jpara,nkprep,nspin)
  real(kind=DP) :: kshift_BP(2)
  logical, private :: ikp0 = .true.

! for the US correction
  real(kind=DP), allocatable, private :: fsr_k0(:,:,:),fsi_k0(:,:,:) ! dim(np_e,nlmta,nspin)
  real(kind=DP), allocatable, private :: fsr_k1(:,:,:),fsi_k1(:,:,:) ! dim(np_e,nlmta,nspin)
  real(kind=DP), allocatable, private :: fsr_k2(:,:,:),fsi_k2(:,:,:) ! dim(np_e,nlmta,nspin
  real(kind=DP), allocatable, private :: zfcos(:),zfsin(:) ! dim(natm)
  real(kind=DP), allocatable, private :: ftqer(:,:,:),ftqei(:,:,:) ! dim(nlmt,nlmt,natm)

! Tags
  character(len("berry_phase")),private,parameter :: tag_berry_phase = "berry_phase"
  character(len("sw_berry_phase")),private,parameter :: tag_sw_berry_phase = "sw_berry_phase"
  character(len("g_index")),private,parameter :: tag_g_index = "g_index"
  character(len("mesh")),private,parameter :: tag_mesh = "mesh"
  character(len("n1")),private,parameter :: tag_n1 = "n1"
  character(len("n2")),private,parameter :: tag_n2 = "n2"
  character(len("j")),private,parameter :: tag_j = "j"
  character(len("kshift")),private,parameter :: tag_kshift = "kshift"
  character(len("k1")),private,parameter :: tag_k1 = "k1"
  character(len("k2")),private,parameter :: tag_k2 = "k2"
  character(len("sw_symmetry")),private,parameter :: tag_sw_symmetry = "sw_symmetry"
  character(len("sw_deficit")),private,parameter :: tag_sw_deficit = "sw_deficit"

! MPI
  !integer, private :: istatus
  integer, private :: istatus(mpi_status_size)

contains

  subroutine m_BP_rd_parameters
    use m_Crystal_Structure, only : rltv
    use m_Control_Parameters,  only : icond,sw_berry_phase
    use m_Const_Parameters, only : FIXED_CHARGE,FIXED_CHARGE_CONTINUATION
    implicit none

    integer :: iret, f_selectBlock, f_getIntValue, f_getRealValue &
         & ,f_getRealVectorValue
    integer :: f_selectParentBlock, f_selectTop
    real(kind=DP) :: dret

    if(icond /= FIXED_CHARGE .and. icond /= FIXED_CHARGE_CONTINUATION) return

    iret = f_selectTop()
    if( f_selectBlock( tag_berry_phase) == 0) then
       if( f_getIntValue( tag_sw_berry_phase, iret) == 0) sw_berry_phase = iret
       if( f_getIntValue( tag_sw_symmetry, iret) == 0) sw_symmetry = iret
       if( f_getIntValue( tag_sw_deficit, iret) == 0) sw_deficit = iret
       if( f_getIntValue( tag_g_index, iret) == 0) then
          ig = iret
       else
          ig = 1 ! default value
       end if
       if( f_selectBlock( tag_mesh) == 0) then
          if( f_getIntValue( tag_n1, iret) == 0) n1 = iret
          if( f_getIntValue( tag_n2, iret) == 0) n2 = iret
          if( f_getIntValue( tag_j, iret) == 0)  jpara = iret
       else
          n1 = 4 ! default value
          n2 = 4 ! default value
          jpara = 20 ! default value
       end if
       iret = f_selectParentBlock()
       if( f_selectBlock( tag_kshift) == 0) then
          if( f_getRealValue( tag_k1, dret, '' ) == 0) kshift_BP(1) = dret
          if( f_getRealValue( tag_k2, dret, '' ) == 0) kshift_BP(2) = dret
       else
          kshift_BP(1:2) = 0.d0 ! default value
       end if
       iret = f_selectParentBlock()
    end if
    iret = f_selectParentBlock()

    if(sw_berry_phase == ON .and. printable) then
       write(nfout,'(1x,"<< Berry Phase")')
       write(nfout,'(1x,"ig=",i1,1x,"g=",3(1x,f10.5))') ig,rltv(1:3,ig)
       write(nfout,'(1x,"n1=",i3,1x,"n2=",i3,1x,"J=",i3)') n1,n2,jpara
       write(nfout,'(1x,"k1=",f8.5,1x,"k2=",f8.5)') kshift_BP(1:2)
       write(nfout,'(1x,"sw_symmetry = ",i1)') sw_symmetry
       write(nfout,'(1x,"sw_deficit = ",i1)') sw_deficit
       write(nfout,'(1x,"   Berry Phase >>")')
    end if

  end subroutine m_BP_rd_parameters

  subroutine m_BP_gen_Kpoints(preallocation)
    use m_Kpoints, only: special_kpoints, kshift, mp_index
    use m_Crystal_Structure, only : rltv, nopr, op, imag
    use m_Const_Parameters,  only : CARTS,BUCS,PARA
    use m_Kpoints, only : kv3, vkxyz, qwgt, m_Kp_alloc_kpoints
    implicit none

    logical, intent(in) :: preallocation

! local variables
    integer :: i,j,k
    integer :: printlevel = 1
    integer :: nkpoint,nsym
    integer, parameter :: nkpmax0 = 10000
    real(kind=DP) :: kpoint(3,nkpmax0), weight(nkpmax0)
    real(kind=DP) :: bmat(3,3)
    real(kind=DP) :: c(3,2),b(3,2)
    real(kind=DP) :: cdot,g2,dg,w
    real(kind=DP) :: a1(3),a2(3),a3(3),volg
    real(kind=DP) :: rot(3,3,48),rot_g(3)
    real(kind=DP) :: scaled_g(3),scaled_a3(3)
    real(kind=DP) :: theta
    real(kind=DP), parameter :: eps = 1.d-12
    logical :: sw_trs = .true. ! time reversal symmetry
    logical :: magnetic, use_trs

    use_trs = .true.

! === KT_add === 2015/03/23
    if ( sw_symmetry == 1 ) use_trs = .false.

    if ( noncol ) then
       use_trs = .false.
    else
       if ( imag /= PARA ) use_trs = .false.
    endif
! ============== 2015/03/23

    if(imag/=PARA) then
       magnetic = .true.
    else
       magnetic = .false.
    end if

    bmat(1:3,1:3)=rltv(1:3,1:3)/PAI2

    g(1:3) = bmat(1:3,ig)
    k=0
    do i=1,3
       if(i /= ig) then
          k=k+1
          b(1:3,k)=bmat(1:3,i)
       end if
    end do

    g2 = sum(g(1:3)**2)
    do k=1,2
       cdot = sum(b(1:3,k)*g(1:3))/g2
       c(1:3,k) = b(1:3,k) - cdot*g(1:3)
    end do

    call cross(c(1,2),g,a1)
    volg = sum(c(1:3,1)*a1(1:3))
    a1(1:3)=a1(1:3)/volg
    call cross(g,c(1,1),a2)
    a2(1:3)=a2(1:3)/volg
    call cross(c(1,1),c(1,2),a3)
    a3(1:3)=a3(1:3)/volg

    nsym=0
    do i=1,nopr
       do j=1,3
          rot_g(j) = sum(op(j,1:3,i)*g(1:3))
       end do
       if(abs(rot_g(1)-g(1)) < eps .and. &
            & abs(rot_g(2)-g(2)) < eps .and. &
            & abs(rot_g(3)-g(3)) < eps ) then
          nsym=nsym+1
          rot(1:3,1:3,nsym) = op(1:3,1:3,i)
       end if
    end do

    if(preallocation .and. printable) then
       write(nfout,'(1x,"<< k-points for an area integral in Berry phase calc.")')
       write(nfout,*) 'c1 = ',c(1:3,1)
       write(nfout,*) 'c2 = ',c(1:3,2)
       write(nfout,*) 'g  = ',g(1:3)
       write(nfout,*) 'a1 = ',a1(1:3)
       write(nfout,*) 'a2 = ',a2(1:3)
       write(nfout,*) 'a3 = ',a3(1:3)
       theta = sum(a1(1:3)*a2(1:3))/(sqrt(sum(a1(1:3)**2))*sqrt(sum(a2(1:3)**2)))
       theta = acos(theta) *360.d0 /PAI2
       write(nfout,*) 'angle between a1 and a2 (deg.) =',theta
    end if

! for identical symmetry olny
    if(sw_symmetry == 1 .or. sw_symmetry == 2) then
       nsym=1
       rot(1:3,1:3,nsym) = 0.d0
       do i=1,3
          rot(i,i,nsym) = 1.d0
       end do
    end if
    if(sw_symmetry == 1) sw_trs = .false.

    mp_index(1) = n1;     mp_index(2) = n2;     mp_index(3) = 1
    kshift(1:2) = kshift_BP(1:2);     kshift(3)=0.d0

    scaled_g(1:3)  = g(1:3)/dble(100.d0)
    scaled_a3(1:3) = a3(1:3)*dble(100.d0)

    nkprep = 0

    call special_kpoints(printlevel,nfout,nsym,rot &
         & , c(1,1),c(1,2),scaled_g,a1,a2,scaled_a3 &
                                !& , c(1,1),c(1,2),g,a1,a2,a3 &
         & , nkprep,nkpmax0,kpoint,weight,magnetic,use_trs)
    if(preallocation.and.printable) then
       write(nfout,'(1x,"   k-points for an area integral in Berry phase calc. >>")')
    end if

    nkvec=nkprep*jpara
    allocate(kvec(3,nkvec))
    allocate(wgh(nkprep))

    dg = 1.d0 /dble(jpara)
    nkvec=0
    do j=1,nkprep
!!!       wgh(j) = weight(j)/dble(nspin)
       wgh(j) = weight(j) /dble( nspin /ndim_spinor )     !!! to be checked !!

       do i=1,jpara
          nkvec=nkvec+1
          kvec(1:2,nkvec) = kpoint(1:2,j)
          kvec(3,nkvec)   = dg*dble(i-1)
       end do
    end do

    !kv3 = nkvec*nspin
    !call m_Kp_alloc_kpoints

    if(.not.preallocation) then

       kv3 = nkvec
       w=1.d0/dble(kv3)
       do i=1,nkvec
          vkxyz(i,1:3,CARTS) = 0.d0
          do j=1,2
             vkxyz(i,1:3,CARTS) = vkxyz(i,1:3,CARTS)+ c(1:3,j)*kvec(j,i)
          end do
          vkxyz(i,1:3,CARTS) = vkxyz(i,1:3,CARTS)+ g(1:3)*kvec(j,i)
          vkxyz(i,1:3,CARTS) = vkxyz(i,1:3,CARTS)*PAI2
          qwgt(i) = w
       end do

       if(printable) then
          write(nfout,'("<< generated k-points")')
          do i=1,nkvec
             write(nfout,'(1x,i4," k=",3(1x,f10.5),1x," k_cart=", &
                  & 3(1x,f10.5),1x,"w= ",f10.5)') &
                  & i,kvec(1:3,i),vkxyz(i,1:3,CARTS),qwgt(i)
          end do
          write(nfout,'("   generated k-points >>")')
       end if

       if ( noncol ) then
          allocate(det(jpara,nkprep,1));det=(0.d0,0.d0)
       else
          allocate(det(jpara,nkprep,nspin));det=(0.d0,0.d0)
       endif

       call set_dk_BP()

    else
       kv3 = nkvec *nspin
       deallocate(kvec,wgh)
    end if

    !stop 'm_BP_gen_Kpoints'

  contains

    subroutine cross(a,b,c)
      implicit none
      real(kind=DP), intent(in) :: a(3),b(3)
      real(kind=DP), intent(out) :: c(3)

!local variables
      integer :: i,j,k

      c(1) = a(2)*b(3)-a(3)*b(2)
      c(2) = a(3)*b(1)-a(1)*b(3)
      c(3) = a(1)*b(2)-a(2)*b(1)

    end subroutine cross

    subroutine set_dk_BP()
      use m_PseudoPotential, only : dk_BP

      dk_BP(1) = sqrt(sum(g(1:3)**2))/ dble(jpara) *PAI2
      !dk_BP(1) = 0.d0 ! for debug of Fourier transforms
      if(printable) then
         write(nfout,*) '!** set_dk_BP: dk_BP=',dk_BP(1)
      end if

    end subroutine set_dk_BP

  end subroutine m_BP_gen_Kpoints

  subroutine calc_phase_for_Berry()
    use m_Ionic_System,        only : pos,natm

! local variables
    integer       :: ia
    real(kind=DP) :: ph

    if(.not.allocated(zfcos)) allocate(zfcos(natm),zfsin(natm))

    ! debug
    !do ia=1,natm
    !  if(printable) write(nfout,*) ia,': pos=',pos(ia,1:3)
    !end do
    !if(printable) write(nfout,*) '--- ia,zfcos,zfsin ---'
    ! end debug

    do ia=1,natm
       ph = pos(ia,ig)*PAI2
       zfcos(ia) = cos(ph)
       zfsin(ia) = sin(ph)
       ! debug
       ! if(printable) write(nfout,*) ia,zfcos(ia),zfsin(ia)
       ! end debug
    end do

  end subroutine calc_phase_for_Berry

  subroutine m_BP_calc_det()
    use m_PseudoPotential,     only : nlmta
! === DEBUG by T.Kato 2013/07/11 ===============================================
    !use m_Electronic_Structure, only : neg, neordr
    use m_Electronic_Structure, only : neordr
    use m_Control_Parameters, only : neg
! ==============================================================================

    implicit none

! local variables
    integer :: i,k,k1,k2,j1,jj1,j2,kg,nk1,nk2
    integer :: ik,is, ismax, js, istmp
    integer :: j1_mpi,j1_ordr
    real(kind=DP) :: mr(nstate,nstate),mi(nstate,nstate)
    real(kind=DP), allocatable :: mr_mpi(:,:),mi_mpi(:,:)
    real(kind=DP), allocatable :: wfbp(:,:,:,:)
    real(kind=DP), allocatable :: fsr_k(:,:,:),fsi_k(:,:,:)

! for debug
    integer :: n

! === KT_add === 2015/03/23
    ismax = nspin /ndim_spinor
! ============== 2015/03/23

    nk1 = mod(ikp,jpara)
    nk2 = ikp/jpara +1

    if(jpara>1) then

       kg = min(kgbp1,kgbp2)

       allocate(wfbp(kgbp2,nstate,2,nspin))
       call copy_old_wfbp_to_temp_wfbp(kgbp2,wfbp2,wfbp) ! MPI

       if(sw_deficit == ON) then
          allocate(fsr_k(nstate,nlmta,nspin),fsi_k(nstate,nlmta,nspin))
          call copy_old_fs_to_temp_fs(fsr_k2,fsi_k2,fsr_k,fsi_k) ! MPI
       end if

! === KT_mod === 2015/03/23
!       SPIN: do is =1,nspin
       SPIN: do is =1, ismax
! ============== 2015/03/23

          mr = 0.d0; mi = 0.d0
          ik = is

          Do js=1, ndim_spinor
             istmp = is +js -1

             do j2=1,nstate
                do j1=1,neg
                   if (map_ek(j1,ik) /= myrank_ke) cycle
                   j1_mpi = map_z(j1)
                   j1_ordr = neordr(j1,ik)

                   if (j1_ordr > nstate) cycle
                   !if(printable) write(nfout,*) '>>>> j1_ordr is lower than nstate+1.'

                   do i=ista_g1k(ik),iend_g1k(ik)
                      mr(j1_ordr,j2) = mr(j1_ordr,j2) &
                           &          + wfbp1(i,j1_mpi,1,istmp) *wfbp(i,j2,1,istmp) &
                           &          + wfbp1(i,j1_mpi,2,istmp) *wfbp(i,j2,2,istmp)
                      mi(j1_ordr,j2) = mi(j1_ordr,j2) &
                           &          - wfbp1(i,j1_mpi,2,istmp) *wfbp(i,j2,1,istmp) &
                           &          + wfbp1(i,j1_mpi,1,istmp) *wfbp(i,j2,2,istmp)
                   end do
                   if (sw_deficit == ON) &
                        & call add_deficit_term( mr(j1_ordr,j2), mi(j1_ordr,j2), &
                        &                        fsr_k1, fsi_k1, fsr_k, fsi_k, &
                        &                        j1_mpi, j2, istmp )
                end do
             end do
          end Do   ! js

          call mpi_barrier(MPI_CommGroup,ierr)
          if (npes>1) then
             allocate(mr_mpi(nstate,nstate),mi_mpi(nstate,nstate))
             mr_mpi=0.d0; mi_mpi=0.d0
             call mpi_allreduce(mr,mr_mpi,nstate**2,mpi_double_precision,mpi_sum, &
                  &             MPI_CommGroup,ierr)
             mr = mr_mpi
             call mpi_allreduce(mi,mi_mpi,nstate**2,mpi_double_precision,mpi_sum, &
                  &             MPI_CommGroup,ierr)
             mi = mi_mpi
             deallocate(mr_mpi,mi_mpi)
          end if

          call determinant(nstate,mr,mi,det(nk1,nk2,is))

          if(ipriberry > 1 .and. printable) then
             write(nfout,'(1x,"Overlap matrix: k1=",i4," k2=",i4," is=",i4)') nk1,nk1+1,is
             write(nfout,'(1x,"ikp=",i4)') ikp
             write(nfout,'(1x,"ikprep=",i4)') nk2
             write(nfout,'(1x,"kg=",i0)') kg
             write(nfout,'(1x,"Real part:")')
             do j1=1,nstate
                write(nfout,'(50(1x,f10.5))') mr(j1,1:nstate)
             end do
             write(nfout,'(1x,"Imag part:")')
             do j1=1,nstate
                write(nfout,'(50(1x,f10.5))') mi(j1,1:nstate)
             end do
             write(nfout,'(1x,"Determinant: k1=",i4," k2=",i4," is=",i4)') nk1,nk1+1,is
             write(nfout,'("Det = ",e25.14," +i ",e25.14)') det(nk1,nk2,is)
          end if

       end do SPIN
       deallocate(wfbp)
       if(sw_deficit==ON) deallocate(fsr_k,fsi_k)

    end if

    if(nk1==jpara-1) then
       kg = min(kgbp2,kgbp0)

       allocate(wfbp(kgbp0,nstate,2,nspin))
       call copy_old_wfbp_to_temp_wfbp(kgbp0,wfbp0,wfbp) ! MPI
       if(sw_deficit == ON) then
          allocate(fsr_k(nstate,nlmta,nspin),fsi_k(nstate,nlmta,nspin))
          call copy_old_fs_to_temp_fs(fsr_k0,fsi_k0,fsr_k,fsi_k) ! MPI
       end if

! === KT_mod === 2015/03/23
!       SPIN2: do is =1,nspin
       SPIN2: do is =1, ismax
! ============== 2015/03/23

          mr =0.d0; mi=0.d0
          ik = is

          Do js=1, ndim_spinor
             istmp = is +js -1

             do j2=1,nstate
                do j1=1,neg
                   if (map_ek(j1,ik) /= myrank_ke) cycle
                   j1_mpi = map_z(j1)
                   j1_ordr = neordr(j1,ik)

                   if (j1_ordr > nstate) cycle
                !if(printable) write(nfout,*) '>>>> j1_ordr is lower than nstate+1.'

                   do i=ista_g1k(ik),iend_g1k(ik)
                      mr(j1_ordr,j2) = mr(j1_ordr,j2) &
                           &          + wfbp2(i,j1_mpi,1,istmp)*wfbp(i,j2,1,istmp) &
                           &          + wfbp2(i,j1_mpi,2,istmp)*wfbp(i,j2,2,istmp)
                      mi(j1_ordr,j2) = mi(j1_ordr,j2) &
                           &          - wfbp2(i,j1_mpi,2,istmp)*wfbp(i,j2,1,istmp) &
                           &          + wfbp2(i,j1_mpi,1,istmp)*wfbp(i,j2,2,istmp)
                   end do
                   if (sw_deficit == ON) &
                        & call add_deficit_term( mr(j1_ordr,j2), mi(j1_ordr,j2), &
                        &                        fsr_k2, fsi_k2, fsr_k, fsi_k, &
                        &                        j1_mpi, j2, istmp )
                end do
             end do
          End Do

          call mpi_barrier(MPI_CommGroup,ierr)

! ==================== modified by K. Tagami ========= 0.1
!          deallocate(wfbp)
!          if(sw_deficit==ON) deallocate(fsr_k,fsi_k)
! ==================================================== 0.1

          if(npes>1) then
             allocate(mr_mpi(nstate,nstate),mi_mpi(nstate,nstate))
             mr_mpi=0.d0; mi_mpi=0.d0
             call mpi_allreduce(mr,mr_mpi,nstate**2,mpi_double_precision,mpi_sum &
                  &,MPI_CommGroup,ierr)
             mr = mr_mpi
             call mpi_allreduce(mi,mi_mpi,nstate**2,mpi_double_precision,mpi_sum &
                  &,MPI_CommGroup,ierr)
             mi = mi_mpi
             deallocate(mr_mpi,mi_mpi)
          end if

          call determinant(nstate,mr,mi,det(jpara,nk2,is))

          if(ipriberry > 1 .and. printable) then
             write(nfout,'(1x,"Overlap matrix: k1=",i4," k2=",i4," is=",i4)') &
                  &                           jpara,jpara+1,is
             write(nfout,'(1x,"ikprep=",i4)') nk2
             write(nfout,'(1x,"Real part:")')
             do j1=1,nstate
                write(nfout,'(50(1x,f10.5))') mr(j1,1:nstate)
             end do
             write(nfout,'(1x,"Imag part:")')
             do j1=1,nstate
                write(nfout,'(50(1x,f10.5))') mi(j1,1:nstate)
             end do

! =================== modified by K. Tagami ====== 0.1
!!            write(nfout,'(1x,"Determinant: k1=",i4," k2=",i4," is=",i4)') jpara,jpara+1
              write(nfout,'(1x,"Determinant: k1=",i4," k2=",i4," is=",i4)') &
                   &                       jpara,jpara+1, is
! ================================================ 0.1
             write(nfout,'("Det = ",e25.14," +i ",e25.14)') det(jpara,nk2,is)
          end if

       end do SPIN2

! ================= Added by K. Tagami ==== 0.1
       deallocate(wfbp)
       if(sw_deficit==ON) deallocate(fsr_k,fsi_k)
! ===========================================

       if(allocated(wfbp0)) deallocate(wfbp0)
       if(allocated(wfbp1)) deallocate(wfbp1)
       if(allocated(wfbp2)) deallocate(wfbp2)
       if(allocated(fsr_k0)) deallocate(fsr_k0,fsi_k0)
       if(allocated(fsr_k1)) deallocate(fsr_k1,fsi_k1)
       if(allocated(fsr_k2)) deallocate(fsr_k2,fsi_k2)

    end if

  contains

    subroutine add_deficit_term(mr,mi,fr1,fi1,fr2,fi2,j1,j2,ispin)
      use m_Const_Parameters,    only : SKIP
      use m_Ionic_System,        only : natm, ityp
      use m_PseudoPotential,     only : index_lmt1_lmt2,lmta &
           &                        , ltp,mtp,taup,ilmt &
           &                        , qitg_BP,iqitg,isph,il2p,dl2p &
           &                        , m_PP_include_vanderbilt_pot &
           &                        , n_non0_lmtxlmt &
           &                        , nlmta
      implicit none
      real(kind=DP), intent(inout) :: mr,mi
      real(kind=DP), intent(in), dimension(np_e,nlmta,nspin) :: fr1,fi1
      real(kind=DP), intent(in), dimension(nstate,nlmta,nspin) :: fr2,fi2
      integer, intent(in) :: j1,j2,ispin

! local varialbes
      integer :: ia,it,ik,ip,lmt1,lmt2,u,v
      real(kind=DP) :: eps = 1.d-13
      real(kind=DP) :: hsr,hsi
      integer :: mdvdb

      ATOM: do ia = ista_atm,iend_atm
         it = ityp(ia)
         mdvdb = m_PP_include_vanderbilt_pot(it)
         if(mdvdb == SKIP) cycle ATOM
         ! debug
         !  if(printable) write(nfout,*) 'ia=',ia,' it=',it
         ! end debug
         LMT_1: do lmt1=1,ilmt(it)
            u = lmta(lmt1,ia)
            !LMT_2: do lmt2=lmt1,ilmt(it)
            LMT_2: do lmt2=1,ilmt(it)
               v = lmta(lmt2,ia)
               hsr = fr1(j1,u,ispin)*fr2(j2,v,ispin) + fi1(j1,u,ispin)*fi2(j2,v,ispin)
               hsi = fr1(j1,u,ispin)*fi2(j2,v,ispin) - fi1(j1,u,ispin)*fr2(j2,v,ispin)
               mr = mr + ftqer(lmt1,lmt2,ia)*hsr-ftqei(lmt1,lmt2,ia)*hsi
               mi = mi + ftqei(lmt1,lmt2,ia)*hsr+ftqer(lmt1,lmt2,ia)*hsi
               ! debug
               !if(printable) then
               !  write(nfout,*) '---- ftqe & hs ----'
               !  write(nfout,*) 'lmt1 = ',lmt1
               !  write(nfout,*) 'lmt2 = ',lmt2
               !  write(nfout,*) 'ftqe = ',ftqer(lmt1,lmt2,ia),ftqei(lmt1,lmt2,ia)
               !  write(nfout,*) 'hs   = ',hsr,hsi
               !end if
               ! end debug
            end do LMT_2
         end do LMT_1
      end do ATOM

    end subroutine add_deficit_term

  end subroutine m_BP_calc_det

  subroutine copy_old_wfbp_to_temp_wfbp(kg,old_wfbp,temp_wfbp)
! === DEBUG by T.Kato 2013/07/11 ===============================================
    !use m_Electronic_Structure, only : neg, neordr
    use m_Electronic_Structure, only : neordr
    use m_Control_Parameters, only : neg
! ==============================================================================
    integer, intent(in) :: kg
    real(kind=DP), intent(in) :: old_wfbp(kg,np_e,2,nspin)
    real(kind=DP), intent(out) :: temp_wfbp(kg,nstate,2,nspin)

! local variables
    real(kind=DP), allocatable :: tmp(:,:,:,:)
    integer :: datasize,ie,je,j
    integer :: is,ik

    !datasize = kg*mp_e*2*nspin
    !allocate(tmp(kg,mp_e,2,nspin)); tmp = 0.d0 ! MPI

    ! debug
    !if(printable) then
    !  write(nfout,*) 'copy wfbp to tmp'
    !  write(nfout,*) 'kg=',kg
    !  write(nfout,*) 'mp_e=',mp_e
    !  write(nfout,*) 'datasize = ', datasize
    !end if
    ! end debug

    !call mpi_barrier(MPI_CommGroup,ierr)
    temp_wfbp = 0.d0
    do is=1,nspin
       ik = is
       do ie=1,neg
          if (map_ek(ie,ik) /= myrank_ke) cycle
          je = neordr(ie,ik)
          if(je > nstate) cycle
          temp_wfbp(1:kg,je,1:2,is) = old_wfbp(1:kg,map_z(ie),1:2,is)
       enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,temp_wfbp,kg*nstate*2*nspin,mpi_double_precision,mpi_sum,mpi_g_world,ierr)
    !do j=0,nrank_e-1
    !   if(j==myrank_e) tmp(1:kg,1:np_e,1:2,1:nspin) = old_wfbp(1:kg,1:np_e,1:2,1:nspin)
    !   call mpi_bcast(tmp(1,1,1,1),datasize,mpi_double_precision &
    !        & ,j,mpi_k_world(myrank_k),ierr)
    !   do is=1,nspin
    !      ik=is
    !      do ie = 1, nel_e(j)
    !         je = nis_e(j)-1+ie
    !         je = neordr(je,ik)
    !         if(je > nstate) cycle
    !         temp_wfbp(1:kg,je,1:2,is) = tmp(1:kg,ie,1:2,is)
    !      end do
    !   end do
    !   call mpi_barrier(MPI_CommGroup,ierr)
    !end do

    !deallocate(tmp)

!    call mpi_barrier(MPI_CommGroup,ierr)

  end subroutine copy_old_wfbp_to_temp_wfbp

  subroutine copy_old_fs_to_temp_fs(old_fsr,old_fsi,temp_fsr,temp_fsi)
! === DEBUG by T.Kato 2013/07/11 ===============================================
    !use m_Electronic_Structure, only : neg, neordr
    use m_Electronic_Structure, only : neordr
    use m_Control_Parameters, only : neg
! ==============================================================================
    use m_PseudoPotential,     only : nlmta
    real(kind=DP), intent(in) :: old_fsr(np_e,nlmta,nspin),old_fsi(np_e,nlmta,nspin)
    real(kind=DP), intent(out) :: temp_fsr(nstate,nlmta,nspin),temp_fsi(nstate,nlmta,nspin)

    ! local variables
    real(kind=DP), allocatable :: tmp(:,:,:,:)
    integer :: datasize,ie,je,j
    integer :: is,ik

    temp_fsr = 0.d0;temp_fsi = 0.d0
    do is=1,nspin
       ik = is
       do ie=1,neg
          if (map_ek(ie,ik) /= myrank_ke) cycle
          je = neordr(ie,ik)
          if(je > nstate) cycle
          temp_fsr(je,1:nlmta,is) = old_fsr(map_z(ie),1:nlmta,is)
          temp_fsi(je,1:nlmta,is) = old_fsi(map_z(ie),1:nlmta,is)
       enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,temp_fsr,nstate*nlmta*nspin,mpi_double_precision,mpi_sum,mpi_g_world,ierr)
    call mpi_allreduce(MPI_IN_PLACE,temp_fsi,nstate*nlmta*nspin,mpi_double_precision,mpi_sum,mpi_g_world,ierr)
!    datasize = mp_e*nlmta*2*nspin
!    allocate(tmp(mp_e,nlmta,2,nspin)); tmp = 0.d0 ! MPI

    ! debug
    !if(printable) then
    !  write(nfout,*) 'copy fsr,fsi to tmp'
    !  write(nfout,*) 'mp_e=',mp_e
    !  write(nfout,*) 'nlmta =',nlmta
    !  write(nfout,*) 'datasize = ', datasize
    !end if
    ! end debug

!    call mpi_barrier(MPI_CommGroup,ierr)
!
!    do j=0,nrank_e-1
!       if(j==myrank_e) then
!          tmp(1:np_e,1:nlmta,1,1:nspin) = old_fsr(1:np_e,1:nlmta,1:nspin)
!          tmp(1:np_e,1:nlmta,2,1:nspin) = old_fsi(1:np_e,1:nlmta,1:nspin)
!       end if
!       call mpi_bcast(tmp(1,1,1,1),datasize,mpi_double_precision &
!            & ,j,mpi_k_world(myrank_k),ierr)
!       do is=1,nspin
!          ik=is
!          do ie = 1, nel_e(j)
!             je = nis_e(j)-1+ie
!             je = neordr(je,ik)
!             if(je > nstate) cycle
!             temp_fsr(je,1:nlmta,is) = tmp(ie,1:nlmta,1,is)
!             temp_fsi(je,1:nlmta,is) = tmp(ie,1:nlmta,2,is)
!          end do
!       end do
!       call mpi_barrier(MPI_CommGroup,ierr)
!!    end do
!
!    deallocate(tmp)

  end subroutine copy_old_fs_to_temp_fs

  subroutine constract_of_ftq()
    use m_Const_Parameters,    only : SKIP,zi
    use m_Ionic_System,        only : natm, ityp
    use m_PseudoPotential,     only : index_lmt1_lmt2,lmta &
         &                        , ltp,mtp,taup,ilmt,nlmt &
         &                        , qitg_BP,iqitg,isph,il2p,dl2p &
         &                        , n_non0_lmtxlmt &
         &                        , m_PP_include_vanderbilt_pot &
         &                        , m_PP_find_maximum_l

! local varialbes
    integer :: ia,it,ik,ip,lmt1,lmt2,u,v
    real(kind=DP) :: fac
    integer :: il1,il2,tau1,tau2,l3
    integer :: n,ilm3,iiqitg
    real(kind=DP) :: dk(3),ylm,dga
    real(kind=DP) :: ftqr,ftqi,ftqb
    integer :: mdvdb
    integer, allocatable :: il3(:)

    call m_PP_find_maximum_l(n)   !  n-1: maximum l
    n = (n-1) + (n-1) + 1
    allocate(il3(n**2)); call substitute_il3(n**2,il3) ! -(b_Elec..)

    allocate(ftqer(nlmt,nlmt,natm)); ftqer=0.d0
    allocate(ftqei(nlmt,nlmt,natm)); ftqei=0.d0

    dk(1:3) = g(1:3) /dble(jpara) *PAI2

    if(ipriberry > 0.and.printable) then
       write(nfout,*) '<< Construction of Qlmt(dk): start >>'
       write(nfout,*) 'dk = ', dk(1:3)
    end if
    ! debug
    !  do ia=1,natm
    !    it = ityp(ia)
    !    if(printable) write(nfout,'("ia=",i3," it=",i3)') ia,it
    !    do ip=1,n_non0_lmtxlmt(it)
    !    lmt1 = index_lmt1_lmt2(ip,it,1)
    !      il1 = ltp(lmt1,it)
    !      tau1 = taup(lmt1,it)
    !      lmt2 = index_lmt1_lmt2(ip,it,2)
    !      il2 = ltp(lmt2,it)
    !      tau2 = taup(lmt2,it)
    !      do n=1,il2p(lmt1,lmt2,it)
    !        ilm3 = isph(lmt1,lmt2,n,it); l3=il3(ilm3)
    !        iiqitg = iqitg(il1,tau1,il2,tau2,l3+1,it)
    !        if(printable) write(nfout,'("lmt1=",i3,"lmt2=",i3,"n=",i3,"q_BP=",f10.5))') lmt1,lmt2,n,qitg_BP(1,iiqitg)
    !      end do
    !    end do
    !  end do
    ! end debug

    ATOM: do ia = 1, natm
       it = ityp(ia)
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(mdvdb == SKIP) cycle ATOM
       ! debug
       !  if(printable) write(nfout,*) 'ia=',ia,' it=',it
       ! end debug

       LMT_1: do lmt1=1,ilmt(it)
          il1 = ltp(lmt1,it)
          tau1 = taup(lmt1,it)
          LMT_2: do lmt2=lmt1,ilmt(it)
             il2 = ltp(lmt2,it)
             tau2 = taup(lmt2,it)
             ! debug
             !if(printable) then
             !  write(nfout,*) 'lmt1 = ',lmt1
             !  write(nfout,*) 'lmt2 = ',lmt2
             !end if
             ! end debug

             ftqr = 0.d0
             ftqi = 0.d0
             LM3: do n=1,il2p(lmt1,lmt2,it)
                ilm3 = isph(lmt1,lmt2,n,it); l3=il3(ilm3)
                iiqitg = iqitg(il1,tau1,il2,tau2,l3+1,it)
                if(iiqitg == 0) cycle LM3
                call sphrp2_for_Berry(ilm3,dk,ylm)
                ftqb = qitg_BP(1,iiqitg)*dl2p(lmt1,lmt2,n,it)*ylm

                ! debug
                !if(printable) then
                !  write(nfout,*) '---- n = ',n,' ----'
                !  write(nfout,*) 'l3 = ',l3
                !  write(nfout,*) 'qitg_BP = ',qitg_BP(1,iiqitg)
                !  write(nfout,*) 'dl2p    = ',dl2p(lmt1,lmt2,n,it)
                !  write(nfout,*) 'ylm     = ',ylm
                !  write(nfout,*) 'ftqb    = ',ftqb
                !end if
                ! end debug
                if(mod(l3,2)==0) then
                   ftqr=ftqr+real(zi**(-l3))*ftqb
                else
                   ftqi=ftqi+dimag(zi**(-l3))*ftqb
                end if
             end do LM3
             ! remove exponetal factor
             !ftqer(lmt1,lmt2,ia) = ftqr*zfcos(ia)+ftqi*zfsin(ia)

             ftqer(lmt1,lmt2,ia) = ftqr
             ftqer(lmt2,lmt1,ia) = ftqer(lmt1,lmt2,ia)

             !ftqei(lmt1,lmt2,ia) = ftqi*zfcos(ia)-ftqr*zfsin(ia)
             ftqei(lmt1,lmt2,ia) = ftqi
             ftqei(lmt2,lmt1,ia) = ftqei(lmt1,lmt2,ia)

             ! debug
             !if(printable) then
             !  write(nfout,*) '---- ftqe ----'
             !  write(nfout,*) 'ftq  = ',ftqr,ftqi
             !  write(nfout,*) '(zfcos,zfsin) = ',zfcos(ia),zfsin(ia)
             !  write(nfout,*) 'ftqe = ',ftqer(lmt1,lmt2,ia),ftqei(lmt1,lmt2,ia)
             !end if
             ! end debug
          end do LMT_2
       end do LMT_1
    end do ATOM

    deallocate(il3)

    if(ipriberry > 0.and.printable) then
       write(nfout,*) '<< Constuction of Qlmt(dk): end >>'
    end if

  end subroutine constract_of_ftq

  subroutine BP_line_integral(phi,cphi,ikprep,ispin)
    implicit none
    real(kind=DP), intent(out) :: phi
    complex(kind=DP), intent(out) :: cphi
    integer, intent(in) :: ikprep, ispin

! local variables
    integer :: k

    cphi = cmplx(1.d0,0.d0)
    do k=1,jpara
       if(printable) write(nfout,'(1x,"Det = ",e25.12," + i ",e25.12)') &
            &                             det(k,ikprep,ispin)
       cphi = cphi*det(k,ikprep,ispin)
    end do

    phi = dimag(log(cphi))

    if(printable) then
       write(nfout,'("cphi = ",e25.12,1x,e25.12)') cphi
       write(nfout,'("phi = ",e25.12)') phi
    end if

  end subroutine BP_line_integral

  subroutine m_BP_Polarization()
    use m_Files, only : nfout, nfzaj
    use m_Electronic_Structure, only : zaj_l, occup_l
    use m_Crystal_Structure, only : univol
    use m_Control_Parameters, only : nspin
    use m_Ionic_System,        only : pos,natm

! =============================== KT_add =================== 13.0C
    use m_Ionic_System,        only : ityp
    use m_PseudoPotential, only : ival
! ========================================================== 13.0C

    implicit none

! local variables
    integer :: i,is,ia, ismax
    real(kind=DP) :: glen, polar, dir(3)
    real(kind=DP), allocatable :: phi(:,:) ! dim(nkprep)
    complex(kind=DP), allocatable :: cphi(:,:) ! dim(nkprep)
    real(kind=DP) :: fac
    real(kind=DP) :: berry_phase, berry_ion, berry_tot
    real(kind=DP) :: polar_ele, polar_ion, polar_tot, polar_max

    if(mype /= 0) return ! MPI

! === KT_add === 2015/03/23
    ismax = nspin /ndim_spinor
! ============== 2015/03/23

    allocate(phi(nkprep,ismax),cphi(nkprep,ismax))

    do is=1,ismax
       do i=1,nkprep
          if(printable) &
               & write(nfout,'("is=",i3," ik=",i3," kprep=",2(1x,f10.5))') &
               & is,i,kvec(1:2,i*jpara)
          call BP_line_integral(phi(i,is),cphi(i,is),i,is)
       end do
    end do
    glen = sqrt(sum(g(1:3)**2))

! =============================== KT_add  =================== 13.0C
    berry_phase = 0.0d0
! =========================================================== 13.0C

    do is = 1, ismax
       berry_phase = berry_phase + sum(phi(1:nkprep,is)*wgh(1:nkprep))
    end do

    fac = 2.d0/dble(nspin)
! === KT_Add === 2015/03/23
    if ( noncol ) fac = 2.d0/dble(nspin /ndim_spinor)    !!! to be checked
! ============== 2015/03/23

    berry_phase = fac*berry_phase
    berry_phase = 4.d0*atan(tan(0.25d0*berry_phase))
    polar_ele = berry_phase/(univol*glen*PAI2)
    berry_ion = 0.d0

    do ia=1,natm
! ============================ KT_mod =================== 13.0C
!       berry_ion = berry_ion + pos(ia,ig)
       berry_ion = berry_ion + ival(ityp(ia)) *pos(ia,ig)
! ======================================================= 13.0C
    end do

    berry_ion = PAI2*berry_ion
    berry_ion = 4.d0*atan(tan(0.25d0*berry_ion))
    polar_ion = berry_ion/(univol*glen*PAI2)
    berry_tot = 4.d0*atan(tan(0.25d0*(berry_ion+berry_phase)))
    polar_tot = berry_tot/(univol*glen*PAI2)
    polar_max = 1.d0/(univol*glen)

    if(printable) then
       write(nfout,'("<< Results of Berry phase calc.")')

       do i=1,nkprep
          write(nfout,'("ik=",i3," kprep=",2(1x,f10.5))') i,kvec(1:2,i*jpara)
       end do
       write(nfout,'(1x,3x,"is",3x,"k",13x,"phi",18x,"wgh")')

       do is=1, ismax
          do i=1,nkprep
             write(nfout,'(1x,i4,1x,i4,2(1x,f20.10))') is,i,phi(i,is),wgh(i)
          end do
       end do

       dir(1:3) = g(1:3)/glen
       write(nfout,'(1x,"Unit cell volume        =",1x,e25.14)') univol
       write(nfout,'(1x,"direction of G          =",3(1x,e25.14))') dir(1:3)
       write(nfout,'(1x,"Berry phase (ele)       =",1x,e25.14)') berry_phase
       write(nfout,'(1x,"Berry phase (ion)       =",1x,e25.14)') berry_ion
       write(nfout,'(1x,"Berry phase (total)     =",1x,e25.14)') berry_tot
       write(nfout,'(1x,"Electronic polarization =",1x,e25.14)') polar_ele
       write(nfout,'(1x,"Ionic polarization      =",1x,e25.14)') polar_ion
       write(nfout,'(1x,"Total polarization      =",1x,e25.14)') polar_tot
       write(nfout,'(1x,"Maximun polarization    =",1x,e25.14)') polar_max

       write(nfout,'("   Results of Berry phase calc.>>")')
    end if

    call m_BP_write_Berry_phase(cphi,phi,nkprep)

    deallocate(phi)

  end subroutine m_BP_Polarization

  subroutine m_BP_write_Berry_phase(cphi,phi,nkprep)
    use m_Files, only : nfout, nfberry, m_Files_open_nfberry
    use m_Ionic_System, only : displaced_atom,displacement,pos,natm
    use m_Crystal_Structure, only : strain, sw_strained_cell
    implicit none
    integer, intent(in) :: nkprep
    complex(kind=DP), intent(in) :: cphi(nkprep,nspin)
    real(kind=DP), intent(in) :: phi(nkprep,nspin)

! local variables
    integer :: i,ia,is, ismax

! === KT_add === 2015/03/23
    ismax = nspin /ndim_spinor
! ============== 2015/03/23

    call m_Files_open_nfberry()
    write(nfberry,'(3(1x,i5),3(1x,e25.14))') nkprep,ig,displaced_atom,displacement(1:3)

    do is=1, ismax
       do i=1,nkprep
          write(nfberry,'(i5,4(1x,e25.14))') i,cphi(i,is),phi(i,is),wgh(i)
       end do
    end do
    if(sw_strained_cell == ON) then
       write(nfberry,'(3(1x,e25.14))') strain(1,1:3)
       write(nfberry,'(3(1x,e25.14))') strain(2,1:3)
       write(nfberry,'(3(1x,e25.14))') strain(3,1:3)
       write(nfberry,*) natm
       do ia=1,natm
          write(nfberry,'(i7,3(1x,e25.14))') ia,pos(ia,1:3)
       end do
    end if

  end subroutine m_BP_write_Berry_phase

  subroutine m_BP_write_wfbp
    use m_Files, only : nfout, nfzaj
    use m_PlaneWaveBasisSet, only : nbase, iba, kgp, kg
    use m_Electronic_Structure, only : totch, zaj_l
    use m_Control_Parameters, only : nspin


    implicit none

! local variables
    integer :: jj,n,j,jmax,ik,ierr
    real(kind=DP), allocatable, dimension(:,:,:,:) :: zajtmp

    ! debug
    !integer :: is
    !real(kind=DP) :: norm,defi
    ! end debug

    if(nstate <= 0) nstate = int(totch + 1.d-13)/2

    if(sw_deficit == ON) then
       !  if(ikp==0) then
       if(ikp0) then
          call calc_phase_for_Berry()
          call constract_of_ftq()
          ikp0 = .false.
       end if
       call m_BP_save_betar_dot_WFs()
    end if

    if(allocated(wfbp2)) then
       if(allocated(wfbp1)) deallocate(wfbp1)
       kgbp1 = kgbp2
       allocate(wfbp1(kgbp1,np_e,2,nspin))
       wfbp1 = wfbp2
       deallocate(wfbp2)
    end if

    kgbp2 = nbase(iba(1),1)
    allocate(wfbp2(kgbp2,np_e,2,nspin))
    allocate(zajtmp(kg,np_e,nspin,2));zajtmp=0.d0
    do n=1,np_e
       do ik=1,nspin
          do j=ista_g1k(ik),iend_g1k(ik)
             zajtmp(j,n,ik,1) = zaj_l(j-ista_g1k(ik)+1,n,ik,1)
             zajtmp(j,n,ik,2) = zaj_l(j-ista_g1k(ik)+1,n,ik,2)
          enddo
       enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,zajtmp,kg*np_e*2*nspin,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)

    if(printable) then
       write(nfout,'("<< Writing Wavefunctions for Berry Phase calc. >>")')
       !write(nfout,*) 'debug iba=',iba(1)
       !write(nfout,*) 'debug kgbp2=',kgbp2
    end if
    do n=1,np_e
       !if(printable) write(nfout,'(1x,"n=",i3)') n
       jj=1
       do j=1,kgbp2
          if(nbase(jj,1) == j) then
    !         wfbp2(j,n,1,1:nspin) = zaj_l(jj,n,1:nspin,1)
    !         wfbp2(j,n,2,1:nspin) = zaj_l(jj,n,1:nspin,2)
             wfbp2(j,n,1,1:nspin) = zajtmp(jj,n,1:nspin,1)
             wfbp2(j,n,2,1:nspin) = zajtmp(jj,n,1:nspin,2)
             jj = jj + 1
          else
             wfbp2(j,n,1,1:nspin) = 0.d0
             wfbp2(j,n,2,1:nspin) = 0.d0
          end if
          ! debug
          !if(j < 21.and.printable) write(nfout,'("j=",i2,"wfbp=",2(1x,f10.5))') &
          !   & j,wfbp2(j,n,1,is),wfbp2(j,n,2,is)
          ! end debug
       end do
    end do

    ! debug
    !if(printable) then
    !write(nfout,*) 'berry check norm: kgbp2=',kgbp2
    !do is=1,nspin
    !do n=1,np_e
    !  defi = norm_deficit(n,is)
    !  write(nfout,'(1x,"n=",i3," defi = ",f15.10)') n,defi
    !  norm = 0.d0
    !  do j=1,kgbp2
    !    norm = norm + wfbp2(j,n,1)**2 + wfbp2(j,n,2)**2
    !  end do
    !  write(nfout,'(1x,"n=",i3," norm = ",f10.5)') n,norm
    !  write(nfout,'(1x,"n=",i3," total norm = ",f15.10)') n,norm+defi
    !end do
    !end do
    !call mpi_barrier(MPI_CommGroup,ierr)
    !end if
    ! end debug

    if((jpara>1 .and. mod(ikp,jpara) /= 0) .or. jpara==1)  &
         & call m_BP_calc_det()

    ikp = ikp + 1
    deallocate(zajtmp)
  contains

    function norm_deficit(n,is) result(norm)
      use m_PseudoPotential,     only : nlmta
      use m_Electronic_Structure,only : fsr_l,fsi_l ! dim(np_e,nlmta,ista_k:iend_k)
      use m_Ionic_System,        only : natm, ityp
      use m_PseudoPotential,     only : n_non0_lmtxlmt,index_lmt1_lmt2,lmta &
           &                          , ltp,mtp,ilmt,q

      real(kind=DP) :: norm ! a result

      integer, intent(in) :: n ! a band index
      integer, intent(in) :: is ! a spin index

! local varialbes
      integer :: ia,it,ik,ip,lmt1,lmt2,u,v
      real(kind=DP) :: eps = 1.d-13
      real(kind=DP) :: fac
      real(kind=DP), allocatable, dimension(:,:) :: fsr2d,fsi2d

      norm = 0.d0

      ik=is ! for ekcal mode
      allocate(fsr2d(np_e,nlmta));fsr2d=0.d0
      allocate(fsi2d(np_e,nlmta));fsi2d=0.d0
      call gather_f_3d_to_2d(fsr_l,fsr2d,ik)
      call gather_f_3d_to_2d(fsi_l,fsi2d,ik)
      do ia = 1, natm
         it = ityp(ia)
         ! debug
         !  if(printable) write(nfout,*) 'n_non0_lmtxlmt = ',n_non0_lmtxlmt(it)
         ! end debug
         LMT: do ip = 1, n_non0_lmtxlmt(it)
            lmt1 = index_lmt1_lmt2(ip,it,1)
            lmt2 = index_lmt1_lmt2(ip,it,2)
            ! debug
            !if(printable) then
            !  write(nfout,*) 'lmt1 = ',lmt1
            !  write(nfout,*) 'lmt2 = ',lmt2
            !end if
            ! end debug

            if(lmt2 < lmt1) cycle LMT
            u = lmta(lmt1,ia)
            v = lmta(lmt2,ia)
            if(ltp(lmt1,it) /= ltp(lmt2,it) .or. &
                 & mtp(lmt1,it) /= mtp(lmt2,it) ) cycle LMT
            fac = q(lmt1,lmt2,it)
            if(abs(fac) < eps) cycle LMT
            if(lmt1 /= lmt2) fac=2.d0*fac

            ! debug
            !if(printable) then
            !  write(nfout,*) 'u = ',u
            !  write(nfout,*) 'v = ',v
            !  write(nfout,*) 'fsr = ',fsr_l(n,u,ik)
            !  write(nfout,*) 'fsi = ',fsi_l(n,v,ik)
            !  write(nfout,*) 'q = ',q(lmt1,lmt2,it)
            !end if
            ! end debug
            norm = norm + &
                 & fac*(fsr2d(n,u)*fsr2d(n,v)+fsi2d(n,u)*fsi2d(n,v))
         end do LMT
      end do
      deallocate(fsr2d)
      deallocate(fsi2d)
    end function norm_deficit

  end subroutine m_BP_write_wfbp

  subroutine gather_f_3d_to_2d(in, out, ik)
    use m_PseudoPotential, only : nlmta
    use m_Parallelization, only : np_fs,nrank_g
    implicit none
    real(kind=DP), dimension(np_e, np_fs,ista_k:iend_k), intent(in)  :: in  !in_3D
    real(kind=DP), dimension(np_e, nlmta),                        intent(out) :: out !out_2D
    integer, intent(in) :: ik
    integer(kind=4) :: ierr, i
    integer :: sendcount, recvcounts(0:nrank_g-1), displs(0:nrank_g-1)
    sendcount = np_e*np_fs
    do i = 0, nrank_g - 1
       recvcounts(i) = np_e*nel_fs(i)
       displs(i)     = np_e*(nis_fs(i) - 1)
    enddo

    out=0.d0
    call mpi_allgatherv(in(1,1,ik),sendcount,        MPI_DOUBLE_PRECISION, &
                        out(1,1),  recvcounts,displs,MPI_DOUBLE_PRECISION,mpi_ke_world,ierr)
    end subroutine gather_f_3d_to_2d

  subroutine m_BP_write_wfbp_gshift
    use m_Files, only : nfout, nfzaj
    use m_PlaneWaveBasisSet, only : kgp, kg, ngabc_kngp_l, nbase, iba
    use m_Electronic_Structure, only : totch, zaj_l
    use m_Control_Parameters, only : nspin
    implicit none

! local variables
    integer :: jj,n,i,j,k,jmax,ik,ierr
    integer :: igs(3)
    integer :: iga,igb,igc,jga,jgb,jgc
    real(kind=DP), allocatable, dimension(:,:,:,:) :: zajtmp
    integer,allocatable,dimension(:,:) :: ngabc_kngp

    ! debug
    !real(kind=DP) :: norm
    ! end debug

    if(mod(ikp,jpara)/=0) return

    allocate(ngabc_kngp(kgp,3));ngabc_kngp=0.0
    do i=ista_kngp,iend_kngp
       ngabc_kngp(i,:) = ngabc_kngp_l(i,:)
    enddo
    call mpi_allreduce(mpi_in_place,ngabc_kngp,kgp*3,mpi_integer,mpi_sum &
            &            ,mpi_chg_world,ierr)

    if(nstate <= 0) nstate = int(totch + 1.d-13)/2

    igs(1:3) = 0
    do i=1,3
       if(i == ig) igs(i) = 1
    end do

    kgbp0 = 1
    jj=1
    if(printable) then
       write(nfout,'("<< Writing Wavefunctions for Berry Phase calc. >>")')
       !write(nfout,*) 'debug kgp=',kgp
       !write(nfout,*) 'debug iba=',iba(1)
    end if
    loop_i: do i=1,kgp
       iga = ngabc_kngp(i,1)+igs(1)
       igb = ngabc_kngp(i,2)+igs(2)
       igc = ngabc_kngp(i,3)+igs(3)
       if(jj>iba(1)) then
          kgbp0 =i-1
          exit loop_i
       end if
       loop_j: do j=1,iba(1)
          k=nbase(j,1)
          jga = ngabc_kngp(k,1)
          jgb = ngabc_kngp(k,2)
          jgc = ngabc_kngp(k,3)
          if(iga == jga .and. igb == jgb .and. igc == jgc) then
             jj=jj+1
             exit loop_j
          end if
       end do loop_j
    end do loop_i

    !if(printable) then
    !write(nfout,*) 'debug jj=',jj
    !write(nfout,*) 'debug kgbp0=',kgbp0
    !end if

    allocate(wfbp0(kgbp0,np_e,2,nspin)); wfbp0 = 0.d0
    allocate(zajtmp(kg,np_e,nspin,2));zajtmp=0.d0
    do n=1,np_e
       do ik=1,nspin
          do j=ista_g1k(ik),iend_g1k(ik)
             zajtmp(j,n,ik,1) = zaj_l(j-ista_g1k(ik)+1,n,ik,1)
             zajtmp(j,n,ik,2) = zaj_l(j-ista_g1k(ik)+1,n,ik,2)
          enddo
       enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,zajtmp,kg*np_e*2*nspin,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
    do n=1,np_e
       loop_i2: do i=1,kgbp0
          iga = ngabc_kngp(i,1)+igs(1)
          igb = ngabc_kngp(i,2)+igs(2)
          igc = ngabc_kngp(i,3)+igs(3)
          loop_j2: do j=1,iba(1)
             k=nbase(j,1)
             jga = ngabc_kngp(k,1)
             jgb = ngabc_kngp(k,2)
             jgc = ngabc_kngp(k,3)
             if(iga == jga .and. igb == jgb .and. igc == jgc) then
                !wfbp0(i,n,1,1:nspin) = zaj_l(j,n,1:nspin,1)
                !wfbp0(i,n,2,1:nspin) = zaj_l(j,n,1:nspin,2)
                wfbp0(i,n,1,1:nspin) = zajtmp(j,n,1:nspin,1)
                wfbp0(i,n,2,1:nspin) = zajtmp(j,n,1:nspin,2)
                exit loop_j2
             end if
          end do loop_j2
       end do loop_i2
       ! debug
       !if(printable) then
       !do j=1,20
       !  write(nfout,'("j=",i2,"wfbp0=",2(1x,f10.5))') &
       !     & j,wfbp0(j,n,1),wfbp0(j,n,2)
       !end do
       !end if
       ! end debug
    end do
    call mpi_barrier(MPI_CommGroup,ierr)
    deallocate(zajtmp)
    deallocate(ngabc_kngp)
    ! debug
    !if(printable) then
    !write(nfout,*) 'berry check norm (gshift): kgbp0=',kgbp0
    !do n=1,np_e
    !  norm = 0.d0
    !  do j=1,kgbp0
    !    norm = norm + wfbp0(j,n,1)**2 + wfbp0(j,n,2)**2
    !  end do
    !  write(nfout,'(1x,"n=",i3," norm = ",f10.5)') n,norm
    !end do
    !call mpi_barrier(MPI_CommGroup,ierr)
    !end if
    ! end debug

  end subroutine m_BP_write_wfbp_gshift

  subroutine m_BP_save_betar_dot_WFs
    use m_Electronic_Structure,only : fsr_l, fsi_l ! dim(np_e,nlmta,ista_k:iend_k)
    use m_PseudoPotential,     only : nlmta, lmta, ilmt
    use m_Ionic_System,        only : natm, ityp
    implicit none

! local variables
    integer :: ik,ia,i,j,n,nkp,ikprep,is
    real(kind=DP), allocatable, dimension(:,:) :: fsr2d,fsi2d

    ik=1 ! for ekcal mode
    nkp = mod(ikp,jpara) +1
    ikprep = ikp /jpara +1

    ! debug
    !if(printable) then
    !  write(nfout,*) 'debug(m_BP_save_betar_dot_WFs)'
    !  write(nfout,*) 'nstate=',nstate
    !  write(nfout,*) 'nlmta=',nlmta
    !  write(nfout,*) 'jpara=',jpara
    !  write(nfout,*) 'nkp=',nkp
    !  write(nfout,*) 'ikprep=',ikprep
    !end if
    ! end debug
    allocate(fsr2d(np_e,nlmta));fsr2d=0.d0
    allocate(fsi2d(np_e,nlmta));fsi2d=0.d0

    if(allocated(fsr_k2)) then
       if(allocated(fsr_k1)) deallocate(fsr_k1,fsi_k1)
       allocate(fsr_k1(np_e,nlmta,nspin),fsi_k1(np_e,nlmta,nspin))
       fsr_k1 = fsr_k2
       fsi_k1 = fsi_k2
       deallocate(fsr_k2,fsi_k2)
    end if

    allocate(fsr_k2(np_e,nlmta,nspin),fsi_k2(np_e,nlmta,nspin))
    do is=1,nspin
       call gather_f_3d_to_2d(fsr_l,fsr2d,is)
       call gather_f_3d_to_2d(fsi_l,fsi2d,is)
       do i=1,nlmta
          do n=1,np_e
             fsr_k2(n,i,is) = fsr2d(n,i)
             fsi_k2(n,i,is) = fsi2d(n,i)
          end do
       end do
    end do
    if(nkp == 1) then
       allocate(fsr_k0(np_e,nlmta,nspin),fsi_k0(np_e,nlmta,nspin))
       do is=1,nspin
          call gather_f_3d_to_2d(fsr_l,fsr2d,is)
          call gather_f_3d_to_2d(fsi_l,fsi2d,is)
          do ia=1,natm
             do j=1,ilmt(ityp(ia))
                i=lmta(j,ia)
                do n=1,np_e
                   fsr_k0(n,i,is) = fsr2d(n,i)*zfcos(ia)+fsi2d(n,i)*zfsin(ia)
                   fsi_k0(n,i,is) = fsi2d(n,i)*zfcos(ia)-fsr2d(n,i)*zfsin(ia)
                end do
             end do
          end do
       end do
    end if

    deallocate(fsr2d)
    deallocate(fsi2d)

  end subroutine m_BP_save_betar_dot_WFs

  subroutine determinant(n,mr,mi,det)
    implicit none
    integer, intent(in) :: n
    real(kind=DP), intent(inout) :: mr(n,n),mi(n,n)
    complex(kind=DP), intent(out) :: det

! local variables
    integer :: i,pivot(n),parity
    real(kind=DP) :: dr,di,dt

    ! debug
    !if(printable) then
    !write(nfout,'("before LU decomposition")')
    !write(nfout,'("real part:")')
    !do i=1,n
    !   write(nfout,'(50(1x,e25.12))') mr(i,1:n)
    !end do
    !write(nfout,'("imag part:")')
    !do i=1,n
    !   write(nfout,'(50(1x,e25.12))') mi(i,1:n)
    !end do
    !end if
    ! end debug

    call LU_decomposition(n,mr,mi,pivot,parity)

    ! debug
    !if(printable) then
    !write(nfout,'("after LU decomposition")')
    !write(nfout,'("real part:")')
    !do i=1,n
    !   write(nfout,'(50(1x,e25.12))') mr(i,1:n)
    !end do
    !write(nfout,'("imag part:")')
    !do i=1,n
    !   write(nfout,'(50(1x,e25.12))') mi(i,1:n)
    !end do
    !end if
    ! end debug

    dr=1.d0
    di=0.d0
    do i=1,n
       ! debug
       !if(printable) write(nfout,'(1x,i3," mr=",e25.12," mi=",e25.12)') i,mr(i,i),mi(i,i)
       ! end debug
       dt = dr*mr(i,i) - di*mi(i,i)
       di = di*mr(i,i) + dr*mi(i,i)
       dr = dt
    end do
    det = cmplx(dr,di)
    det = parity*det

  end subroutine determinant

  subroutine LU_decomposition(n,mr,mi,pivot,parity)
    implicit none
    integer, intent(in) :: n
    real(kind=DP), intent(inout) :: mr(n,n),mi(n,n)
    integer, intent(out) :: pivot(n),parity

! local variables
    integer :: i,j,k,imax
    real(kind=DP) :: big,scaling(n)
    real(kind=DP) :: t,dr,di
    real(kind=DP) :: sumr,sumi,mt

    real(kind=DP), parameter :: eps = 3.5d-15

    parity = 1

    do i=1,n
       big = 0.d0
       do j=1,n
          t = sqrt( mr(i,j)**2 + mi(i,j)**2 )
          if(t > big)  big = t
       end do
       if(big == 0.d0) then
          if(printable) then
             write(nfout,*) 'Error(LU_decomposition): Singular matrix'
             write(nfout,*) 'i=',i
          end if
          call phase_error_with_msg(nfout,'LU_decomposition',__LINE__,__FILE__)
       else
          scaling(i) = 1.d0/big
       end if
    end do

    do j=1,n
       do i=1,j-1
          sumr=0.d0
          sumi=0.d0
          do k=1,i-1
             sumr = sumr + ( mr(i,k)*mr(k,j) - mi(i,k)*mi(k,j) )
             sumi = sumi + ( mi(i,k)*mr(k,j) + mr(i,k)*mi(k,j) )
          end do
          mr(i,j) = mr(i,j) - sumr
          mi(i,j) = mi(i,j) - sumi
       end do
       big = 0.d0
       do i=j,n
          sumr=0.d0
          sumi=0.d0
          do k=1,j-1
             sumr = sumr + ( mr(i,k)*mr(k,j) - mi(i,k)*mi(k,j) )
             sumi = sumi + ( mi(i,k)*mr(k,j) + mr(i,k)*mi(k,j) )
          end do
          sumr = mr(i,j) - sumr
          sumi = mi(i,j) - sumi
          t = scaling(i)*sqrt(sumr**2+sumi**2)
          if( t >= big) then
             big = t
             imax = i
          end if
          mr(i,j) = sumr
          mi(i,j) = sumi
       end do
       if(j /= imax) then
          do k=1,n
             dr = mr(imax,k)
             di = mi(imax,k)
             mr(imax,k) = mr(j,k)
             mi(imax,k) = mi(j,k)
             mr(j,k) = dr
             mi(j,k) = di
          end do
          parity = -parity
          scaling(imax) = scaling(j)
       end if
       pivot(j) = imax
       if( dsqrt(mr(j,j)**2+mi(j,j)**2) < eps ) then
          mr(j,j) = 0.d0
          mi(j,j) = 0.d0
       end if
       if(j /= n ) then
          t  =  mr(j,j)**2 + mi(j,j)**2
          dr =  mr(j,j)/t
          di = -mi(j,j)/t
          do i=j+1,n
             mt      = mr(i,j)*dr - mi(i,j)*di
             mi(i,j) = mi(i,j)*dr + mr(i,j)*di
             mr(i,j) = mt
          end do
       end if
    end do

  end subroutine LU_decomposition

  subroutine m_BP_wd_cntn_data()
    use m_Files, only : nfcntn_berry,m_Files_open_nfcntn_berry, &
         &              m_Files_close_nfcntn_berry,nfout
    use m_ErrorMessages, only : EOF_REACHED

    logical :: logica

    if(.not.allocated(det)) return
    call m_Files_open_nfcntn_berry()
    rewind nfcntn_berry
    if(mype==0)then
       write(nfcntn_berry) ikp
       write(nfcntn_berry) kgbp0
       write(nfcntn_berry) kgbp1
       write(nfcntn_berry) kgbp2
       write(nfcntn_berry) det
    endif

! fsr, fsi related
    logica = allocated(fsr_k2)
    if(mype==0) write(nfcntn_berry) logica
    if(logica) call write_fs(fsr_k2,nfcntn_berry)

    logica = allocated(fsi_k2)
    if(mype==0) write(nfcntn_berry) logica
    if(logica) call write_fs(fsi_k2,nfcntn_berry)

    logica = allocated(fsr_k1)
    if(mype==0) write(nfcntn_berry) logica
    if(logica) call write_fs(fsr_k1,nfcntn_berry)

    logica = allocated(fsi_k1)
    if(mype==0) write(nfcntn_berry) logica
    if(logica) call write_fs(fsi_k1,nfcntn_berry)

    logica = allocated(fsr_k0)
    if(mype==0) write(nfcntn_berry) logica
    if(logica) call write_fs(fsr_k0,nfcntn_berry)

    logica = allocated(fsi_k0)
    if(mype==0) write(nfcntn_berry) logica
    if(logica) call write_fs(fsi_k0,nfcntn_berry)

! wf related
    logica = allocated(wfbp2)
    if(mype==0) write(nfcntn_berry) logica
    if(logica) call write_wf(kgbp2,wfbp2,nfcntn_berry)

    logica = allocated(wfbp1)
    if(mype==0) write(nfcntn_berry) logica
    if(logica) call write_wf(kgbp1,wfbp1,nfcntn_berry)

    logica = allocated(wfbp0)
    if(mype==0) write(nfcntn_berry) logica
    if(logica) call write_wf(kgbp0,wfbp0,nfcntn_berry)

    call m_Files_close_nfcntn_berry()

  end subroutine m_BP_wd_cntn_data

  subroutine m_BP_rd_cntn_data()
    use m_Files,           only : nfcntn_berry,m_Files_open_nfcntn_berry, &
         &                        m_Files_close_nfcntn_berry,nfout
    use m_PseudoPotential, only : nlmta
    use m_ErrorMessages,   only : EOF_REACHED

    logical :: logica

    if(.not.allocated(det)) return
    call m_Files_open_nfcntn_berry()
    rewind nfcntn_berry
    if(mype==0)then
       read(nfcntn_berry,end=9999,err=9999) ikp
       read(nfcntn_berry,end=9999,err=9999) kgbp0
       read(nfcntn_berry,end=9999,err=9999) kgbp1
       read(nfcntn_berry,end=9999,err=9999) kgbp2
       read(nfcntn_berry,end=9999,err=9999) det
    endif
    call mpi_bcast(ikp,1,mpi_double_precision,0,MPI_CommGroup,ierr)
    call mpi_bcast(kgbp0,1,mpi_double_precision,0,MPI_CommGroup,ierr)
    call mpi_bcast(kgbp1,1,mpi_double_precision,0,MPI_CommGroup,ierr)
    call mpi_bcast(kgbp2,1,mpi_double_precision,0,MPI_CommGroup,ierr)
    call mpi_bcast(det,jpara*nkprep*nspin,mpi_double_precision,0,MPI_CommGroup,ierr)

    if(mype==0) read(nfcntn_berry,end=9999,err=9999) logica
    call mpi_bcast(logica,1,mpi_logical,0,MPI_CommGroup,ierr)
    if (logica) then
       if(.not.allocated(fsr_k2)) allocate(fsr_k2(np_e,nlmta,nspin))
       call read_fs(fsr_k2,nfcntn_berry,1)
    endif

    if(mype==0) read(nfcntn_berry,end=9999,err=9999) logica
    call mpi_bcast(logica,1,mpi_logical,0,MPI_CommGroup,ierr)
    if (logica) then
       if(.not.allocated(fsi_k2)) allocate(fsi_k2(np_e,nlmta,nspin))
       call read_fs(fsi_k2,nfcntn_berry,2)
    endif

    if(mype==0) read(nfcntn_berry,end=9999,err=9999) logica
    call mpi_bcast(logica,1,mpi_logical,0,MPI_CommGroup,ierr)
    if (logica) then
       if(.not.allocated(fsr_k1)) allocate(fsr_k1(np_e,nlmta,nspin))
       call read_fs(fsr_k1,nfcntn_berry,3)
    endif

    if(mype==0) read(nfcntn_berry,end=9999,err=9999) logica
    call mpi_bcast(logica,1,mpi_logical,0,MPI_CommGroup,ierr)
    if (logica) then
       if(.not.allocated(fsi_k1)) allocate(fsi_k1(np_e,nlmta,nspin))
       call read_fs(fsi_k1,nfcntn_berry,4)
    endif

    if(mype==0) read(nfcntn_berry,end=9999,err=9999) logica
    call mpi_bcast(logica,1,mpi_logical,0,MPI_CommGroup,ierr)
    if (logica) then
       if(.not.allocated(fsr_k0)) allocate(fsr_k0(np_e,nlmta,nspin))
       call read_fs(fsr_k0,nfcntn_berry,5)
    endif

    if(mype==0) read(nfcntn_berry,end=9999,err=9999) logica
    call mpi_bcast(logica,1,mpi_logical,0,MPI_CommGroup,ierr)
    if (logica) then
       if(.not.allocated(fsi_k0)) allocate(fsi_k0(np_e,nlmta,nspin))
       call read_fs(fsi_k0,nfcntn_berry,6)
    endif

    if(mype==0) read(nfcntn_berry,end=9999,err=9999) logica
    call mpi_bcast(logica,1,mpi_logical,0,MPI_CommGroup,ierr)
    if (logica) then
       if(allocated(wfbp2)) deallocate(wfbp2)
       allocate(wfbp2(kgbp2,np_e,2,nspin))
       call read_wf(kgbp2,wfbp2,nfcntn_berry,7)
    endif

    if(mype==0) read(nfcntn_berry,end=9999,err=9999) logica
    call mpi_bcast(logica,1,mpi_logical,0,MPI_CommGroup,ierr)
    if(logica) then
       if(allocated(wfbp1)) deallocate(wfbp1)
       allocate(wfbp1(kgbp1,np_e,2,nspin))
       call read_wf(kgbp1,wfbp1,nfcntn_berry,8)
    endif

    if(mype==0) read(nfcntn_berry,end=9999,err=9999) logica
    call mpi_bcast(logica,1,mpi_logical,0,MPI_CommGroup,ierr)
    if(logica) then
       if(allocated(wfbp0)) deallocate(wfbp0)
       allocate(wfbp0(kgbp0,np_e,2,nspin))
       call read_wf(kgbp0,wfbp0,nfcntn_berry,9)
    endif

    call m_Files_close_nfcntn_berry()
    return

9999 continue
    call phase_error_wo_filename(EOF_REACHED, nfout, nfcntn_berry, __LINE__, __FILE__)
  end subroutine m_BP_rd_cntn_data

  subroutine write_fs(fs,nf)
!!$  use m_Electronic_Structure, only : neg
    use m_Control_Parameters, only : neg
    use m_PseudoPotential, only : nlmta

    real(kind=DP),intent(in),dimension(np_e,nlmta,nspin) :: fs
    integer, intent(in) :: nf
    real(kind=DP),allocatable,dimension(:,:) :: f
    integer :: is,ib
    integer :: ierr

    allocate(f(nlmta,nspin));f=0.0d0
    do is=1,nspin
       do ib=1,neg
          if(map_ek(ib,is) == mype) then
             f(:,:) = fs(map_z(ib),:,:)
             if(map_ek(ib,is)/=0) call mpi_send(f,nlmta*nspin,mpi_double_precision,0,1,MPI_CommGroup,ierr)
          else if (mype==0 .and. map_ek(ib,is) /= 0)then
             call mpi_recv(f,nlmta*nspin,mpi_double_precision,map_ek(ib,is),1,MPI_CommGroup,istatus,ierr)
          endif
          if(mype==0) write(nf) f
       enddo
    enddo
    deallocate(f)

  end subroutine write_fs

  subroutine read_fs(fs,nf,tag)
!!$  use m_Electronic_Structure, only : neg
    use m_Control_Parameters, only : neg
    use m_PseudoPotential, only : nlmta
    use m_ErrorMessages, only : EOF_REACHED
    use m_Files, only : nfout

    real(kind=DP),intent(out),dimension(np_e,nlmta,nspin) :: fs
    integer, intent(in) :: nf
    integer, intent(in) :: tag
    real(kind=DP),allocatable,dimension(:,:) :: f
    integer :: is,ib
    integer :: ierr

    allocate(f(nlmta,nspin));f=0.0d0
    fs = 0.d0
    do is=1,nspin
       do ib=1,neg
          if(mype==0) read(nf,end=9999,err=9999) f
          if(mype==0 .and. map_ek(ib,is)/=0)then
             call mpi_send(f,nlmta*nspin,mpi_double_precision,map_ek(ib,is),tag,MPI_CommGroup,ierr)
          else if (map_ek(ib,is)==mype .and. map_ek(ib,is) /= 0)then
             call mpi_recv(f,nlmta*nspin,mpi_double_precision,0,tag,MPI_CommGroup,istatus,ierr)
          endif
          if(map_ek(ib,is) == mype) fs(map_z(ib),:,:) = f(:,:)
       enddo
    enddo
    deallocate(f)
    return

9999 continue
    call phase_error_wo_filename(EOF_REACHED, nfout, nf, __LINE__, __FILE__)

  end subroutine read_fs

  subroutine write_wf(n1,wf,nf)
    use m_Control_Parameters, only : neg

    integer, intent(in) :: n1
    real(kind=DP), dimension(n1,np_e,2,nspin),intent(in) :: wf
    integer, intent(in) :: nf
    real(kind=DP), allocatable,dimension(:,:) :: wft
    integer :: is,ib
    integer :: ierr

    allocate(wft(n1,2));wft=0.d0
    do is=1,nspin
       do ib=1,neg
          if(map_ek(ib,is) == mype)then
             wft(:,1) = wf(:,map_z(ib),1,is)
             wft(:,2) = wf(:,map_z(ib),2,is)
             if(map_ek(ib,is)/=0) call mpi_send(wft,n1*2,mpi_double_precision,0,1,MPI_CommGroup,ierr)
          else if (mype==0 .and. map_ek(ib,is) /= 0) then
             call mpi_recv(wft,n1*2,mpi_double_precision,map_ek(ib,is),1,MPI_CommGroup,istatus,ierr)
          endif
          if(mype==0) write(nf) wft
       enddo
    enddo
    deallocate(wft)

  end subroutine write_wf

  subroutine read_wf(n1,wf,nf,tag)
    use m_Control_Parameters, only : neg
    use m_ErrorMessages, only : EOF_REACHED
    use m_Files, only : nfout

    integer, intent(in) :: n1
    real(kind=DP), dimension(n1,np_e,2,nspin),intent(out) :: wf
    integer, intent(in) :: nf
    integer, intent(in) :: tag
    real(kind=DP), allocatable,dimension(:,:) :: wft
    integer :: is,ib
    integer :: ierr

    wf = 0.d0
    allocate(wft(n1,2));wft=0.d0
    do is=1,nspin
       do ib=1,neg
          if(mype==0) read(nf,end=9999,err=9999) wft
          if(mype==0 .and. map_ek(ib,is)/=0)then
             call mpi_send(wft,n1*2,mpi_double_precision,map_ek(ib,is),tag,MPI_CommGroup,ierr)
          else if (map_ek(ib,is) == mype .and. map_ek(ib,is) /= 0) then
             call mpi_recv(wft,n1*2,mpi_double_precision,0,tag,MPI_CommGroup,istatus,ierr)
          endif
          if(map_ek(ib,is) == mype) then
             wf(:,map_z(ib),1,is) = wft(:,1)
             wf(:,map_z(ib),2,is) = wft(:,2)
          endif
       enddo
    enddo
    deallocate(wft)
    return

9999 continue
    call phase_error_wo_filename(EOF_REACHED, nfout, nf, __LINE__, __FILE__)

  end subroutine read_wf

end module m_BerryPhase

