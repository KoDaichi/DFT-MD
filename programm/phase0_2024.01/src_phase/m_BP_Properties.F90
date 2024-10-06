#ifdef HIUX
*option MP(P(0))
#endif
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
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
module m_BP_Properties
  ! $ID: $
  !
  ! Properties:
  ! 1. Berry-phase polarization
  !
  !    P_a = -qe/V sum_j Z_ion(j) * r(j)_a - qe/V * f/(2*pi) sum_i A_ia * phi_i
  !
  ! 2. Effective charge tensor
  !
  !    Z_ab = Z_ion + f/(2*pi) sum_i A_ia * d(phi_i)/d(u_b)
  !
  ! 3. piezoelectirc stress tensor
  !
  !    C_ab = f/(2*pi*V) sum_i A_ia * d(phi_i)/d(e_b)
  !

  use m_Files, only : nfout
  use m_Control_Parameters, only : polar_prop,nspin,printable, ndim_spinor
  use m_Const_Parameters,   only : DP, PAI, PAI2, ON, OFF &
       &      , POLARIZATION, EFFECTIVE_CHARGE &
       &      , PIEZOELECTRIC_CONST, UNIT_PIEZO_CONST
  use m_Parallelization,    only : MPI_CommGroup,mype,npes,ierr
  use mpi

  implicit none
!  include 'mpif.h'

  integer, private :: num_berry_phase,max_num_kpoints
  integer, private, allocatable :: num_kpoints(:) ! dim(num_berry_phase)
  integer, private, allocatable :: g_index(:) ! dim(num_berry_phase)
  integer, private, allocatable :: direction(:) ! dim(num_berry_phase)
  integer, private, allocatable :: displaced_atom(:) ! dim(num_berry_phase)
  ! 1 -> x, 2 -> y, 3 -> z
  real(kind=DP), private, allocatable :: u(:) ! dim(num_berry_phase)

! atomic displacement
  real(kind=DP), private, allocatable :: displacement(:,:) ! dim(3,num_berry_phase)
  complex(kind=DP), private, allocatable :: cphi(:,:,:) ! dim(nkprep_max,nspin,num_berry_phase)
  real(kind=DP), private, allocatable :: wgh(:,:) ! dim(nkprep_max,num_berry_phase)
  real(kind=DP), private, allocatable :: epsstr(:) ! dim(num_berry_phase)
  real(kind=DP), private, allocatable :: posstr(:,:,:) ! dim(natm,3,num_berry_phase)
  integer, private, allocatable :: istrain(:) ! dim(num_berry_phase)

! Properties
  real(kind=DP), allocatable :: zeff(:,:,:) ! dim(3,3,natm)
  real(kind=DP), allocatable :: piezo(:,:) ! dim(3,6)

contains

  subroutine m_BP_read_Berry_phase(paramset)
    use m_Files, only : m_Files_open_nfberry, nfberry
    use m_Ionic_System, only : natm
    implicit none

    logical, intent(in) :: paramset

! local variables
    integer :: i,k,ia,is, ismax
    integer :: idummy
    real(kind=DP) :: dummy,cphi_re,cphi_im,phi
    real(kind=DP) :: d(3),tmp,epstmp(3,3)

    ! debug
    !write(nfout,*) 'DEBUG: reading berry.data (1)'
    ! end debug

! === KT_add === 2015/03/23
       ismax = nspin /ndim_spinor
! ============== 2015/03/23

    if(.not.paramset) then
       allocate(num_kpoints(1),g_index(1),displaced_atom(1),displacement(3,1))
       allocate(wgh(1,1))
    else
       ! debug
       !   write(nfout,*) 'DEBUG: num_berry_phase = ',num_berry_phase
       !   write(nfout,*) 'DEBUG: max_num_kpoints = ',max_num_kpoints
       ! end debug
       allocate(num_kpoints(num_berry_phase))
       allocate(g_index(num_berry_phase))
       allocate(displaced_atom(num_berry_phase))
       allocate(displacement(3,num_berry_phase))


! === KT_mod === 2015/03/23
!       allocate(cphi(max_num_kpoints,nspin,num_berry_phase))
       allocate(cphi(max_num_kpoints,ismax,num_berry_phase))
! ============== 2015/03/23

       allocate(wgh(max_num_kpoints,num_berry_phase))
       allocate(direction(num_berry_phase))
       allocate(u(num_berry_phase))

       if( polar_prop == PIEZOELECTRIC_CONST ) then
          allocate(istrain(num_berry_phase))
          allocate(epsstr(num_berry_phase))
          allocate(posstr(natm,3,num_berry_phase))
       end if
    end if

    ! debug
    !write(nfout,*) 'DEBUG: reading berry.data (2)'
    ! end debug

    if(mype == 0) then

       if(.not.paramset) max_num_kpoints = 0

       call m_Files_open_nfberry()
       rewind nfberry
       read(nfberry,*) num_berry_phase

       ! debug
!!$       write(nfout,*) 'DEBUG: num_berry_phase = ',num_berry_phase;      call flush(nfout)
       ! end debug

       if(.not.paramset) then
          do i=1,num_berry_phase
             ! debug
!!$             write(nfout,'("m_BP_read_Berry_phase: "i7,"-th berry data is reading.")') i;   call flush(nfout)
             ! end debug
             read(nfberry,*) num_kpoints(1), g_index(1), displaced_atom(1), &
                  &                                     displacement(1:3,1)
             do is=1, ismax
                do k=1,num_kpoints(1)
                   read(nfberry,*) idummy,cphi_re,cphi_im,phi,wgh(1,1)
                end do
             end do
             if(num_kpoints(1) .gt. max_num_kpoints) then
                max_num_kpoints = num_kpoints(1)
             end if
             if( polar_prop == PIEZOELECTRIC_CONST ) then
                do k=1,3
                   read(nfberry,*)
                end do
                read(nfberry,*) idummy
                do k=1,idummy
                   read(nfberry,*)
                end do
             end if
          end do
          ! debug
!!$          write(nfout,*) 'DEBUG: max_num_kpoints = ',max_num_kpoints;  call flush(nfout)
          ! end debug
       else
          do i=1,num_berry_phase
             ! debug
!!$             write(nfout,'("m_BP_read_Berry_phase: "i7,"-th berry data is reading.")') i;   call flush(nfout)
             ! end debug

             read(nfberry,*) num_kpoints(i),g_index(i), displaced_atom(i), &
                  &                                     displacement(1:3,i)

             ! debug
!!$             write(nfout,'("DEBUG: the header was read.")');       call flush(nfout)
             ! end debug

             if(polar_prop == EFFECTIVE_CHARGE) then
                direction(i) = 1
                u(i) = 0.d0
                do k=1,3
                   if(abs(displacement(k,i)) > 1.d-10) then
                      ! debug
                      !write(nfout,'("m_BP_read_Berry_phase: displacement = ",e25.14)') displacement(k,i)
                      direction(i) = k
                      u(i) = displacement(k,i)
                      !write(nfout,'("m_BP_read_Berry_phase: "i7,"-th data. u=",e25.14)') i,u(i)
                      ! end debug
                   end if
                end do

                d(1:3) = displacement(1:3,i)
                d(direction(i)) = d(direction(i)) - u(i)
                if(abs(sum(d(1:3)**2)) > 1.d-10) then
                   if(printable) then
                      write(nfout,'("Displacement vector of ",i4,"-th atom =",3(1x,f10.5))') &
                           & displaced_atom(i),displacement(1:3,i)
                   end if
                   return
                   !stop 'Error of displacement in <<<m_BP_read_Berry_phase>>>'
                end if
                if(displaced_atom(i) <= 0) direction(i) = -1
             end if

             ! debug
!!$             write(nfout,'("DEBUG: berry phase is reading.")');      call flush(nfout)
             ! end debug

             do is=1, ismax
                do k=1,num_kpoints(i)
                   read(nfberry,*) idummy,cphi_re,cphi_im,phi,wgh(k,i)
                   cphi(k,is,i) = dcmplx(cphi_re,cphi_im)
                end do
             end do

             if( polar_prop == PIEZOELECTRIC_CONST ) then
                ! debug
!!$                write(nfout,'("DEBUG: strain tensor is reading.")');    call flush(nfout)
                ! end debug
                read(nfberry,*) epstmp(1,1:3)
                read(nfberry,*) epstmp(2,1:3)
                read(nfberry,*) epstmp(3,1:3)
                call set_istrain(epstmp,istrain(i),epsstr(i))
                read(nfberry,*) idummy
                if(idummy /= natm) then
                   call phase_error_with_msg(nfout,'natm in berry.data differs from natm in PHASE.',__LINE__,__FILE__)
                end if
                do ia=1,natm
                   read(nfberry,*) idummy,posstr(ia,1:3,i)
                end do
             end if
             if(printable) then
                write(nfout,'(i7,"-th berry data was read.")') i
                call flush(nfout)
             end if
          end do
       end if

    end if

    if(.not.paramset) then
       deallocate(num_kpoints,g_index,displaced_atom,displacement)
       deallocate(wgh)
    end if

    if(npes > 1) then
       if(.not.paramset) then
          call mpi_bcast(num_berry_phase,1,mpi_integer,0,MPI_CommGroup,ierr) ! MPI
          call mpi_bcast(max_num_kpoints,1,mpi_integer,0,MPI_CommGroup,ierr) ! MPI
       else
          call mpi_bcast(num_kpoints,num_berry_phase,mpi_integer,0,MPI_CommGroup,ierr)
          call mpi_bcast(g_index,num_berry_phase,mpi_double_precision,0,MPI_CommGroup,&
               &         ierr) ! MPI
          call mpi_bcast(displaced_atom,num_berry_phase,mpi_integer,0,MPI_CommGroup,&
               &         ierr) ! MPI
          call mpi_bcast(displacement,num_berry_phase*3,mpi_double_precision,0,&
               &         MPI_CommGroup,ierr) ! MPI
!!$          call mpi_bcast(cphi,max_num_kpoints*ismax*num_berry_phase*2,&
          call mpi_bcast(cphi,max_num_kpoints*ismax*num_berry_phase,&
               &         mpi_double_complex,0,MPI_CommGroup,ierr) ! MPI
          call mpi_bcast(wgh,max_num_kpoints*num_berry_phase,mpi_double_precision, &
               &         0,MPI_CommGroup,ierr) ! MPI
          call mpi_bcast(direction,num_berry_phase,mpi_integer,0,MPI_CommGroup,ierr)
          call mpi_bcast(u,num_berry_phase,mpi_double_precision,0,MPI_CommGroup,ierr)
       end if
    end if

  contains

    subroutine set_istrain(eps,istrain,str)
      real(kind=DP), intent(in) :: eps(3,3)
      integer, intent(out) :: istrain
      real(kind=DP), intent(out) :: str

      integer :: i
      real(kind=DP) :: dst(6)

      dst(1) = eps(1,1)
      dst(2) = eps(2,2)
      dst(3) = eps(3,3)
      dst(4) = eps(2,3)+eps(3,2)
      dst(5) = eps(3,1)+eps(1,3)
      dst(6) = eps(1,2)+eps(1,2)

      ! debug
      if(printable) then
         write(nfout,*) 'debug: set_istrain'
         do i=1,6
            write(nfout,*) i,dst(i)
         end do
         write(nfout,*) 'end debug: set_istrain'
      end if
      ! end debug

      istrain = 0
      do i=1,6
         if(abs(dst(i)) > 0.d0) then
            istrain = i
            str = dst(i)
            return
         end if
      end do

    end subroutine set_istrain

  end subroutine m_BP_read_Berry_phase

  subroutine m_BP_get_Berry_phase(atom_id,g_id,dir_id,cphi_out,wgh_out,nkp,u_out,exists)
    implicit none
    integer, intent(in) :: atom_id, g_id, dir_id
    complex(kind=DP), intent(out) :: cphi_out(max_num_kpoints,nspin/ndim_spinor)
    real(kind=DP), intent(out) :: wgh_out(max_num_kpoints)
    integer, intent(out) :: nkp
    real(kind=DP), intent(out) :: u_out
    logical, intent(out) :: exists

! lcoal variables
    integer :: i,k,is, ismax

! === KT_add === 2015/03/23
    ismax = nspin /ndim_spinor
! ============== 2015/03/23

    exists = .false.
    do i=1,num_berry_phase
       ! debug
       !write(6,*) 'atom_id,g_id,dir_id=',atom_id,g_id,dir_id,': i=',i
       ! end debug
       if(displaced_atom(i).eq.atom_id .and. &
            & g_index(i).eq.g_id .and. &
            & direction(i).eq.dir_id) then
          nkp = num_kpoints(i)

          do is=1,ismax
             do k=1,nkp
                cphi_out(k,is) = cphi(k,is,i)
                wgh_out(k) = wgh(k,i)
             end do
          end do
          u_out = u(i)
          exists = .true.
          return
       end if
    end do

    return
  end subroutine m_BP_get_Berry_phase

  subroutine m_BP_calc_diff_Berry_phase(dphi,present_atom)
    use m_Ionic_System, only : natm
    implicit none
    real(kind=DP), intent(out) :: dphi(3,3,natm)
    logical, intent(out) :: present_atom(natm)

! local variables
    integer :: ia,ig,id,k,ia0,is, ismax
    integer :: nkp0,nkp1
    complex(kind=DP) :: cphi0(max_num_kpoints,nspin/ndim_spinor)
    complex(kind=DP) :: cphi1(max_num_kpoints,nspin/ndim_spinor)
    real(kind=DP) :: wgh0(max_num_kpoints)
    real(kind=DP) :: wgh1(max_num_kpoints)
    real(kind=DP) :: phi(max_num_kpoints)
    real(kind=DP) :: u0,u1
!!$logical :: exists
    logical :: exists1, exists2

    present_atom(1:natm) = .false.

! === KT_add === 2015/03/23
    ismax = nspin /ndim_spinor
! ============== 2015/03/23

    ATOM: do ia=1,natm
!!$       write(6,'(" ia = ",i8," <<m_BP_calc_diff_Berry_phase>>")') ia
       do ig=1,3
          ia0 = 0
          id = -1
!!$call m_BP_get_Berry_phase(ia0,ig,id,cphi0,wgh0,nkp0,u0,exists)
!!$if(.not.exists) cycle ATOM
          call m_BP_get_Berry_phase(ia0,ig,id,cphi0,wgh0,nkp0,u0,exists1)
          do id=1,3
!!$call m_BP_get_Berry_phase(ia,ig,id,cphi1,wgh1,nkp1,u1,exists)
!!$if(.not.exists) cycle ATOM
             call m_BP_get_Berry_phase(ia,ig,id,cphi1,wgh1,nkp1,u1,exists2)
             if(.not.(exists1.and.exists2)) cycle
             present_atom(ia) = .true.
             if(nkp0.ne.nkp1) then
                if(printable) then
                   write(nfout,*) 'nkp0=',nkp0
                   write(nfout,*) 'nkp1=',nkp1
                end if
                call phase_error_with_msg(nfout,'nkp0 and nkp1 are different in <<<m_BP_calc_diff_Berry_phase>>'&
                                         ,__LINE__,__FILE__)
             end if

             ! debug
             !   write(nfout,*) 'ia,ig,id =',ia,ig,id
             ! end debug

             dphi(id,ig,ia) = 0.d0

             do is=1,ismax
                do k=1,nkp0
                   if(wgh0(k).ne.wgh1(k)) then
                      if(printable) then
                         write(nfout,*) 'wgh0=',wgh0(k)
                         write(nfout,*) 'wgh1=',wgh1(k)
                      end if
                      call phase_error_with_msg(nfout,'wgh0 and wgh1 are different in <<<m_BP_calc_diff_Berry_phase>>'&
                                               ,__LINE__,__FILE__)
                   end if
                   phi(k) = dimag(log(cphi1(k,is)/cphi0(k,is)))
                   ! debug
                   !     write(nfout,*) 'k=',k,' phi=',phi(k)
                   ! end debug
                   dphi(id,ig,ia) = dphi(id,ig,ia) + phi(k)*wgh0(k)
                end do
             end do
             ! debug
             !   write(nfout,*) 'ia,id,ig=',id,ig,ia,' dphi=',dphi(id,ig,ia)
             ! end debug
             dphi(id,ig,ia) = 2.d0*dphi(id,ig,ia)/u1
             ! debug
             !   write(nfout,*) 'ia,id,ig=',id,ig,ia,' dphi/u=',dphi(id,ig,ia)
             ! end debug
          end do
       end do
!!$present_atom(ia) = .true.
    end do ATOM
    return

  end subroutine m_BP_calc_diff_Berry_phase

  subroutine m_BP_calc_Polarization()
    use m_Crystal_Structure, only : altv
    use m_Ionic_System, only : natm,ntyp,ityp,cps
    use m_Crystal_Structure, only : univol
    use m_PseudoPotential, only : ival
    implicit none

! local variables
    integer :: i,ia,nkp,is, ismax
    real(kind=DP) :: polar_ion(3), polar_ele(3), polar_tot(3)
    real(kind=DP) :: phi(3)
    real(kind=DP) :: fac

    do ia=1,natm
       polar_ion(1:3) = polar_ion + ival(ityp(ia))*cps(ia,1:3)
    end do
    fac = 2.d0/(PAI2*univol)
    polar_ele(1:3) = 0.d0

! === KT_add === 2015/03/23
    ismax = nspin /ndim_spinor
! ============== 2015/03/23

    do i=1,3
       nkp = num_kpoints(i)
       phi(i) = 0.d0
       do is=1, ismax
          phi(i) = phi(i) +  sum(dimag(log(cphi(1:nkp,is,i)))*wgh(1:nkp,i))
       end do
       polar_ele(1:3) = polar_ele(1:3) + fac*altv(1:3,i)*phi(i)
    end do
    polar_tot(1:3) = polar_ion(1:3) + polar_ele(1:3)

! print polarization
    if(printable) then
       write(nfout,*) ' --- Berry phase ---'
       write(nfout,'("phi =",3(1x,f20.10))') phi(1:3)
       write(nfout,*) ' --- Polarization ---'
       write(nfout,'("P_ion =",3(1x,f20.10))') polar_ion(1:3)
       write(nfout,'("P_ele =",3(1x,f20.10))') polar_ele(1:3)
       write(nfout,'("P_tot =",3(1x,f20.10))') polar_tot(1:3)
    end if

    return
  end subroutine m_BP_calc_Polarization

  subroutine m_BP_calc_Effective_charge()
    use m_Crystal_Structure, only : altv,nopr,op
    use m_Ionic_System, only : natm,ntyp,ityp,cps,napt
    use m_Crystal_Structure, only : univol
    use m_PseudoPotential, only : ival
    use m_Files, only : nfeffchg, m_Files_open_nfeffchg
    use m_Control_Parameters,  only : sw_check_sumrule_born_charge

    implicit none

! local variables
    integer :: i,j,k,ia,ja,iopr
    real(kind=DP) :: dphi(3,3,natm)
    logical :: present_atom(natm)
    real(kind=DP) :: zeff_el(3,3,natm)
    !real(kind=DP) :: zeff(3,3,natm)

    real(kind=DP) :: fac
    real(kind=DP) :: ztmp(3,3)

! == ASMS DEBUG ==
!!!real(kind=DP) :: nsym
    integer :: nsym
! ===============

    real(kind=DP) :: nonzero(3,3)

!!$    write(nfout,'(" !! (-1) Effective_charge")');    call flush(nfout)
!!$    write(nfout,*) ' allocated(zeff) = ', allocated(zeff), " natm = ",natm

    if(.not.allocated(zeff)) then
       allocate(zeff(3,3,natm))
       zeff = 0.0d0
    endif
!!$    write(nfout,'(" !! (0) Effective_charge")');    call flush(nfout)

    if( polar_prop == EFFECTIVE_CHARGE ) then

       call m_BP_calc_diff_Berry_phase(dphi,present_atom)
       fac = 1.d0/PAI2
       do ia=1,natm
          if(.not.present_atom(ia)) cycle
          do j=1,3
             do i=1,3
                zeff_el(i,j,ia) = fac*sum(altv(i,1:3)*dphi(j,1:3,ia))
                zeff(i,j,ia) = zeff_el(i,j,ia)
             end do
          end do
          do i=1,3
             zeff(i,i,ia) = zeff(i,i,ia) + ival(ityp(ia))
          end do
       end do

! print electronic effective charges
       if(printable) then
          write(nfout,*) ' --- Calculated electronic effective charges ---'
          do ia=1,natm
             if(.not.present_atom(ia)) cycle
             write(nfout,'("     ",3x,"    [",3(1x,f10.5),1x,"]")') zeff_el(1,1:3,ia)
             write(nfout,'("Zel (",i3,") = [",3(1x,f10.5),1x,"]")') ia,zeff_el(2,1:3,ia)
             write(nfout,'("     ",3x,"    [",3(1x,f10.5),1x,"]")') zeff_el(3,1:3,ia)
             write(nfout,*)
          end do

! print effective charges
          write(nfout,*) ' --- Calculated effective charges ---'
          do ia=1,natm
             if(.not.present_atom(ia)) cycle
             write(nfout,'("     ",3x,"    [",3(1x,f10.5),1x,"]")') zeff(1,1:3,ia)
             write(nfout,'("Zeff(",i3,") = [",3(1x,f10.5),1x,"]")') ia,zeff(2,1:3,ia)
             write(nfout,'("     ",3x,"    [",3(1x,f10.5),1x,"]")') zeff(3,1:3,ia)
             write(nfout,*)
          end do
       end if
    else
       call m_BP_read_Effective_charge(zeff,present_atom)
    end if

!!$    write(nfout,'(" !! (1) Effective_charge")');    call flush(nfout)

! Symmetrization
    ATOM: do ia=1,natm
       if(.not.present_atom(ia)) cycle ATOM
       ztmp(1:3,1:3) = 0.d0
       nsym=0
       OPERATION: do iopr=1,nopr
          if(napt(ia,iopr).ne.ia) cycle OPERATION
          nsym=nsym+1
          do j=1,3
             do i=1,3
                do k=1,3
                   ztmp(i,j) = ztmp(i,j)+sum(op(i,1:3,iopr)*zeff(1:3,k,ia))*op(j,k,iopr)
                end do
             end do
          end do
       end do OPERATION
       zeff(1:3,1:3,ia) = ztmp(1:3,1:3)/nsym
    end do ATOM

!!$    write(nfout,'(" !! (2) Effective_charge")')

! print effective charges
    if(printable) then
       write(nfout,*) ' --- Symmetrized effective charges ---'
       do ia=1,natm
          if(.not.present_atom(ia)) cycle
          write(nfout,'("     ",3x,"    [",3(1x,f10.5),1x,"]")') zeff(1,1:3,ia)
          write(nfout,'("Zsym(",i3,") = [",3(1x,f10.5),1x,"]")') ia,zeff(2,1:3,ia)
          write(nfout,'("     ",3x,"    [",3(1x,f10.5),1x,"]")') zeff(3,1:3,ia)
          write(nfout,*)
       end do
    end if

! calculate effective chareges of equivalent atoms
    ATOM2: do ia=1,natm
       ! debug
       !  print *,'DEGUG: ia=',ia
       if(.not.present_atom(ia)) cycle ATOM2
       OPERATION2: do iopr=1,nopr
          ja = napt(ia,iopr)
          ! debug
          !  print *,'DEGUG: operation iopr=',iopr
          !  print *,'DEGUG: equivalent atom ja=',ja

          if(ja.eq.ia) cycle OPERATION2
          if(present_atom(ja)) cycle OPERATION2
          zeff(1:3,1:3,ja) = 0.d0
          do j=1,3
             do i=1,3
                do k=1,3
                   zeff(i,j,ja) = zeff(i,j,ja) &
                        &        + sum(op(i,1:3,iopr)*zeff(1:3,k,ia))*op(j,k,iopr)
                end do
             end do
          end do
          present_atom(ja) = .true.
       end do OPERATION2
    end do ATOM2

! print effective charges of all atoms
    if(printable) then
       write(nfout,*) ' --- Effective charges of all atoms ---'
       do ia=1,natm
          write(nfout,'("     ",3x,"    [",3(1x,f10.5),1x,"]")') zeff(1,1:3,ia)
          write(nfout,'("Zeff(",i3,") = [",3(1x,f10.5),1x,"]")') ia,zeff(2,1:3,ia)
          write(nfout,'("     ",3x,"    [",3(1x,f10.5),1x,"]")') zeff(3,1:3,ia)
          write(nfout,*)
       end do
    end if

! === KT_add === 2014/06/30
    if ( sw_check_sumrule_born_charge == OFF ) goto 100
! ============== 2014/06/30

! check zero sum rule
    ztmp(1:3,1:3) = 0.d0
    nonzero(1:3,1:3) = 0
    do ia=1,natm
       do j=1,3
          do i=1,3
             if(dabs(zeff(i,j,ia)).gt.1.d-10) nonzero(i,j)=nonzero(i,j)+1
             ztmp(i,j) = ztmp(i,j) + zeff(i,j,ia)
          end do
       end do
    end do
    do j=1,3
       do i=1,3
          if(nonzero(i,j)==0) nonzero(i,j)=1
       end do
    end do
    ztmp(1:3,1:3) = ztmp(1:3,1:3)/nonzero(1:3,1:3)
! print averaged effective charge
    if(printable) then
       write(nfout,*) ' --- Averaged effective charges ---'
       write(nfout,'("       [",3(1x,f10.5),1x,"]")') ztmp(1,1:3)
       write(nfout,'("Zave = [",3(1x,f10.5),1x,"]")') ztmp(2,1:3)
       write(nfout,'("       [",3(1x,f10.5),1x,"]")') ztmp(3,1:3)
       write(nfout,*)
    end if

! correct effective charges
    do ia=1,natm
       do j=1,3
          do i=1,3
             if(dabs(zeff(i,j,ia)).lt.1.d-10) cycle
             zeff(i,j,ia) = zeff(i,j,ia) - ztmp(i,j)
          end do
       end do
    end do

! print corrected effective charges
    if(printable) then
       write(nfout,*) ' --- Corrected effective charges ---'
       do ia=1,natm
          write(nfout,'("     ",3x,"    [",3(1x,f10.5),1x,"]")') zeff(1,1:3,ia)
          write(nfout,'("Zeff(",i3,") = [",3(1x,f10.5),1x,"]")') ia,zeff(2,1:3,ia)
          write(nfout,'("     ",3x,"    [",3(1x,f10.5),1x,"]")') zeff(3,1:3,ia)
          write(nfout,*)
       end do
    end if

100 continue

! ==== KT_add === 13.1R
    if ( mype == 0 ) then
       call m_Files_open_nfeffchg()

       write(nfeffchg,*) natm
       do ia=1,natm
          write(nfeffchg,'(I6)') ia
          write(nfeffchg,'(3F18.12)') zeff(1,1:3,ia)
          write(nfeffchg,'(3F18.12)') zeff(2,1:3,ia)
          write(nfeffchg,'(3F18.12)') zeff(3,1:3,ia)
       end do
       close( nfeffchg )
    endif
! =============== 13.1R
    return

  end subroutine m_BP_calc_Effective_charge

  subroutine m_BP_read_Effective_charge(zeff,present_atom)
    use m_Ionic_System, only : natm
    use m_Files, only : nfeffchg, m_Files_open_nfeffchg
    implicit none
    real(kind=DP), intent(out) :: zeff(3,3,natm)
    logical, intent(out) :: present_atom(natm)

    integer :: num_zeff
    integer :: i,ia

    if( polar_prop == EFFECTIVE_CHARGE ) return

    present_atom(1:natm) = .false.

    call m_Files_open_nfeffchg()
    if(mype == 0) then
       rewind nfeffchg
       read(nfeffchg,*) num_zeff
       do i=1,num_zeff
          read(nfeffchg,*) ia
          read(nfeffchg,*) zeff(1,1:3,ia)
          read(nfeffchg,*) zeff(2,1:3,ia)
          read(nfeffchg,*) zeff(3,1:3,ia)
          present_atom(ia) = .true.
       end do
       close( nfeffchg )  !ASMS
    end if

    if(npes>1) then
       call mpi_bcast(zeff,3*3*natm,mpi_double_precision,0,MPI_CommGroup,ierr) ! MPI
       call mpi_bcast(present_atom,natm,mpi_logical,0,MPI_CommGroup,ierr) ! MPI
    end if

    return
  end subroutine m_BP_read_Effective_charge

  subroutine m_BP_get_Berry_phase_strain(g_id,eps_id,cphi_out,wgh_out,nkp,eps_out, &
       &                                 posstr_out,exists)
    use m_Ionic_System, only : natm
    implicit none
    integer, intent(in) :: g_id, eps_id
    complex(kind=DP), intent(out) :: cphi_out(max_num_kpoints,nspin/ndim_spinor)
    real(kind=DP), intent(out) :: wgh_out(max_num_kpoints)
    integer, intent(out) :: nkp
    real(kind=DP), intent(out) :: eps_out
    real(kind=DP), intent(out) :: posstr_out(natm,3)
    logical, intent(out) :: exists

    integer :: i,k,is, ismax

! === KT_add === 2015/03/23
    ismax = nspin /ndim_spinor
! ============== 2015/03/23

    exists = .false.
    do i=1,num_berry_phase
       if(g_index(i).eq.g_id .and. &
            & istrain(i).eq.eps_id) then
          nkp = num_kpoints(i)

          do is=1,ismax
             do k=1,nkp
                cphi_out(k,is) = cphi(k,is,i)
                wgh_out(k) = wgh(k,i)
             end do
          end do
          eps_out = epsstr(i)
          posstr_out(1:natm,1:3) = posstr(1:natm,1:3,i)
          exists = .true.
          return
       end if
    end do

  end subroutine m_BP_get_Berry_phase_strain

  subroutine m_BP_calc_diff_Berry_strain(dphi,present_strain)
    use m_Ionic_System, only : natm
    implicit none
    real(kind=DP), intent(out) :: dphi(3,6)
    logical, intent(out) :: present_strain(6)

! local variables
    integer :: ia,ig,id,k,ist,is, ismax
    integer :: nkp0,nkp1
    complex(kind=DP) :: cphi0(max_num_kpoints,nspin/ndim_spinor)
    complex(kind=DP) :: cphi1(max_num_kpoints,nspin/ndim_spinor)
    real(kind=DP) :: wgh0(max_num_kpoints)
    real(kind=DP) :: wgh1(max_num_kpoints)
    real(kind=DP) :: phi(max_num_kpoints,nspin/ndim_spinor)
    real(kind=DP) :: eps0,eps1
    real(kind=DP) :: pos0(natm,3),pos1(natm,3)
    logical :: exists
    real(kind=DP) :: dphi_ion, dphi_0, dphi_1

! === KT_add === 2015/03/23
    ismax = nspin /ndim_spinor
! ============== 2015/03/23

    present_strain(1:6) = .false.

    do ig=1,3
       ist = 0
       call m_BP_get_Berry_phase_strain(ig,ist,cphi0,wgh0,nkp0,eps0,pos0,exists)
       if(.not.exists) cycle
       STRAIN: do ist=1,6
          call m_BP_get_Berry_phase_strain(ig,ist,cphi1,wgh1,nkp1,eps1,pos1,exists)
          if(.not.exists) cycle STRAIN
          present_strain(ist) = .true.

! <<< ASMS 2023.10.31
!          if(nkp0.ne.nkp1) then
!             if(printable) then
!                write(nfout,*) 'nkp0=',nkp0
!                write(nfout,*) 'nkp1=',nkp1
!             end if
!             call phase_error_with_msg(nfout,'nkp0 and nkp1 are different in <<<m_BP_calc_diff_Berry_strain>>', &
!             __LINE__,__FILE__)
!          end if
!          ! debug
!          !   write(nfout,*) 'ig,ist =',ig,ist
!          ! end debug
! ASMS 2023.10.31 >>>

          dphi(ig,ist) = 0.d0

! <<< ASMS 2023.10.31
!          do is=1, ismax
!             do k=1,nkp0
!                if(wgh0(k).ne.wgh1(k)) then
!                   if(printable) then
!                      write(nfout,*) 'wgh0=',wgh0(k)
!                      write(nfout,*) 'wgh1=',wgh1(k)
!                   end if
!                   call phase_error_with_msg(nfout,'wgh0 and wgh1 are different in <<<m_BP_calc_diff_Berry_strain>>'&
!                                            ,__LINE__,__FILE__)
!                end if
!                phi(k,is) = dimag(log(cphi1(k,is)/cphi0(k,is)))
!
!                ! debug
!                if(printable) write(nfout,*) 'k=',k,' is=',is,' phi=',phi(k,is)
!                ! end debug
!
!                dphi(ig,ist) = dphi(ig,ist) + phi(k,is)*wgh0(k)
!             end do
!          end do

          dphi_0 = 0.0d0;    dphi_1 = 0.0d0
          do is=1, ismax
             do k=1, nkp0
                dphi_0 = dphi_0 +dimag(log(cphi0(k,is))) *wgh0(k)
             end do
             do k=1, nkp1
                dphi_1 = dphi_1 +dimag(log(cphi1(k,is))) *wgh1(k)
             end do
          end do
          dphi(ig,ist) = dphi_1 -dphi_0
! ASMS 2023.10.31 >>>

          dphi(ig,ist) = 2.d0*dphi(ig,ist)

          ! debug
          if(printable) write(nfout,*) 'ig,ist=',ig,ist,' dphi=',dphi(ig,ist)
          ! end debug

          call ionic_part(pos0,pos1,ig,dphi_ion)

          ! debug
          if(printable) then
             write(nfout,*) 'ig,ist=',ig,ist,' dphi_ion=',dphi_ion
             write(nfout,*) 'ig,ist=',ig,ist,' dphi_tot=',dphi_ion+dphi(ig,ist)
             write(nfout,*) 'ig,ist=',ig,ist,' dphi_tot(mod.)=',4.d0*atan(tan(0.25d0*(dphi_ion + dphi(ig,ist))))
          end if
          ! end debug

          dphi(ig,ist) = 4.d0*atan(tan(0.25d0*(dphi_ion + dphi(ig,ist))))/eps1

          ! debug
          if(printable) then
             write(nfout,*) 'ig,ist=',ig,ist,' eps=',eps1
             write(nfout,*) 'ig,ist=',ig,ist,' dphi/eps=',dphi(ig,ist)
          end if
          ! end debug
       end do STRAIN
    end do

  contains

    subroutine ionic_part(pos0,pos1,ig,dphi)
      use m_Crystal_Structure, only : rltv
      use m_PseudoPotential, only : ival
      use m_Ionic_System, only : ityp
      real(kind=DP), intent(in) :: pos0(natm,3) ! internal coordinates
      real(kind=DP), intent(in) :: pos1(natm,3)
      integer, intent(in) :: ig
      real(kind=DP), intent(out) :: dphi

! local variables
      integer :: ia

      dphi = 0.d0
      do ia=1,natm
         dphi = dphi + ival(ityp(ia))*(pos1(ia,ig) - pos0(ia,ig))*PAI2
      end do

    end subroutine ionic_part

  end subroutine m_BP_calc_diff_Berry_strain

  subroutine m_BP_calc_piezoelectric_const
    use m_Crystal_Structure, only : altv
    use m_Ionic_System, only : natm,ntyp,ityp,cps,napt
    use m_Crystal_Structure, only : univol
    implicit none

! local variables
    integer :: i,ia,ist
    real(kind=DP) :: dphi(3,6)
    real(kind=DP) :: fac
    logical :: present_strain(6)

    if(.not.allocated(piezo)) allocate(piezo(3,6))
    piezo = 0.d0

    if( polar_prop == PIEZOELECTRIC_CONST ) then
       dphi = 0.d0
       call m_BP_calc_diff_Berry_strain(dphi,present_strain)
       fac = 1.d0/(PAI2*univol)
       ! debug
       if(printable) then
          write(nfout,*) 'fac=',fac
          write(nfout,*) 'present_strain=',present_strain(1:6)
       end if
       ! end debug
       do ist=1,6
          if(.not.present_strain(ist)) cycle
          ! debug
          if(printable) then
             write(nfout,*) 'ist=',ist
             write(nfout,'(" dphi=",3(1x,f10.5))') dphi(1:3,ist)
          end if
          ! end debug

         do i=1,3
             ! debug
             if(printable) then
                write(nfout,'(" a1i,a2i,a3i=",3(1x,f10.5))') altv(i,1:3)
             end if
             ! end debug
             piezo(i,ist) = fac*sum(altv(i,1:3)*dphi(1:3,ist))
          end do
       end do
    end if

    if(printable) then
       write(nfout,*) '=== Piezoelectric constant (a.u.) ==='
       do ist=1,6
          write(nfout,'(i5,3(1x,f18.10))') ist,piezo(1:3,ist)
       end do
       write(nfout,*) '=== Piezoelectric constant (C/m^2) ==='
       do ist=1,6
          write(nfout,'(i5,3(1x,f18.10))') ist,piezo(1:3,ist)*UNIT_PIEZO_CONST
       end do
    end if

  end subroutine m_BP_calc_piezoelectric_const

end module m_BP_Properties
