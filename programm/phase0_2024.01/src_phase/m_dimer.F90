
module m_dimer
  use m_Parallelization, only : mype, nrank_conf, mype_conf, npes
  use m_Const_Parameters, only : DP, ON, BOHR, PAI2, SKIP_SYMMCHECK, PAI
  use m_Control_Parameters, only : printable
  use m_Crystal_Structure, only : rltv
  use m_PseudoPotential, only : flg_paw, m_PP_rd_PAW_parameters
  use m_Charge_Density,    only :  m_CD_alloc_chgsoft
  use m_Files,         only :  m_Files_open_nfcntn_bin_paw, nfout, nfcntn_bin_paw
  use m_Ionic_System,  only : neb_dt
  use m_Replica, only : neb_file_type
  use mpi

  implicit none

!  include 'mpif.h'

  type dimer_type
    integer :: max_dimer_iteration
    integer :: natm
    real(kind=DP) :: dtheta, dR
    real(kind=DP) :: dt
    real(kind=DP), allocatable :: r0(:,:), r1(:,:), r2(:,:), for(:,:,:),f1(:,:), f2(:,:), velocity(:,:)
    real(kind=DP), allocatable :: f1para(:,:), f1perp(:,:), f2para(:,:), f2perp(:,:), fpara(:,:), fperp(:,:), fperpold(:,:)
    real(kind=DP), allocatable :: fR(:,:), f(:,:), fold(:,:), ftrans(:,:)
    integer, allocatable :: mobile(:)
    real(kind=DP), allocatable :: N(:,:), theta(:,:), thetaold(:,:), currangle(:), ene(:)
    real(kind=DP) :: curvature
    real(kind=DP) :: e, e1, e2, e0
    real(kind=DP) :: threshold
    real(kind=DP) :: max_force
    integer :: nrank_r, mype_r
    type(neb_file_type) :: file_type
  end type dimer_type

  type(dimer_type), target :: dimvars

  logical :: pp_generated = .false.
  logical :: dimer_cond_set = .false.
  logical :: standard_out_opened = .true.

  contains

  subroutine m_dm_set_pp_generated(ppg)
    logical, intent(in) :: ppg
    pp_generated = ppg
  end subroutine m_dm_set_pp_generated

  subroutine cal_scf_energy_force(i,itr)
    use m_Force, only : forc_l
    use m_Total_Energy, only : etotal
    use m_IterationNumbers, only : iteration_electronic, iteration, iteration_ionic, m_Iter_dimer_incre
    integer, intent(in) :: i, itr
    logical  :: ChargeDensity_is_Converged, TotalEnergy_is_Divergent
    logical  :: Already_Converged, Positron_bulk
    logical  :: Hubbard_model
    logical  :: Forces_are_Converged, Ending_Time, Ending_Time2, Force_errors_are_tolerable
    integer  :: ia

    if(mod(i-1,nrank_conf) == mype_conf) then
    call set_dimer_replica_parameter(i, itr)
    call set_dimer_coordinates(i)
    call Preparation(0)                   ! Basis set, symmetry check etc.
    call Preparation_for_mpi(ON)           ! mpi
    if (.not.pp_generated) then
       call PseudoPotential_Construction
       pp_generated = .true.
#ifdef ENABLE_ESM_PACK
       call Preparation_for_ESM()
#endif
    else
       if ( flg_paw ) then
!          if ( itr > 1 .or. neb%cond%condition == 1 ) then
         if( itr > 1 )then
           call m_PP_rd_PAW_parameters(nfout,nfcntn_bin_paw)
         endif
         call m_CD_alloc_chgsoft
       endif
    endif
! ================================================= 13.0D
    call Ewald_and_Structure_Factor
    call Initial_Electronic_Structure
    call Initial_MD_Condition

    if(.not.Already_Converged()) then
       AtomicConfiguration: do
          ChargeDensity:    do
             call IterationNumber_Setting
             call Renewal_of_WaveFunctions
             call ChargeDensity_Construction(1)
             call ChargeDensity_Mixing
             call Renewal_of_Potential
             if(Hubbard_model()) then
                call Renewal_of_Hubbard_Potential
             end if
             if(Ending_Time()) then
                exit AtomicConfiguration
             end if
             if(TotalEnergy_is_Divergent())    exit ChargeDensity
             if(ChargeDensity_is_Converged())  exit ChargeDensity
          enddo ChargeDensity
#ifdef LIBRARY_BUILD
          call Forces_phase0
#else
          call Forces
#endif
          exit AtomicConfiguration
       enddo AtomicConfiguration
    end if
    dimvars%for(i,:,:) = forc_l(:,:)
    dimvars%ene(i) = etotal
    write(nfout,'(a,i8,f20.10)') '!** dimer method energy for configuration ',i,dimvars%ene(i)
    write(nfout,'(a,i8)')        '!** dimer method forces for configuration ',i
    do ia=1, dimvars%natm
      write(nfout,'(i8,3f20.10)') ia,dimvars%for(i,ia,1:3)
    enddo
    iteration = 0
    iteration_electronic = 0
    iteration_ionic = 0
    call WriteDownData_onto_Files(.false.)
    call deallocate_array
    endif

    call m_Iter_dimer_incre(nfout)
  end subroutine cal_scf_energy_force

  subroutine cal_dimer_scf(itr)
    integer, intent(in) :: itr
    integer :: i, ierr
    dimvars%for=0.d0;dimvars%ene=0.d0
    do i=1,2
      call cal_scf_energy_force(i,itr)
    enddo
    if(nrank_conf>1) then
      call mpi_allreduce(mpi_in_place, dimvars%for, 2*dimvars%natm*3, mpi_double_precision, &
              mpi_sum, mpi_comm_world, ierr)
      dimvars%for = dimvars%for/dble(npes)
      call mpi_allreduce(mpi_in_place, dimvars%ene, 2, mpi_double_precision, &
              mpi_sum, mpi_comm_world, ierr)
      dimvars%ene = dimvars%ene/dble(npes)
    endif
  end subroutine cal_dimer_scf

  subroutine deallocate_array

    use m_Control_Parameters, only : m_CtrlP_set_init_status
    use m_Crystal_Structure, only : m_CS_dealloc
    use m_Ionic_System, only : m_IS_dealloc
    use m_PlaneWaveBasisSet, only : m_pwBS_dealloc
    use m_Parallelization, only : m_Parallel_dealloc
    use m_Kpoints, only : m_Kp_dealloc
    use m_Force, only : m_Force_dealloc
    use m_Charge_Density, only : m_CD_dealloc
    use m_XC_Potential, only : m_XC_dealloc_vxc
!    use m_PseudoPotential
!    use m_NonLocal_Potential,
    use m_Electronic_Structure, only : m_ES_dealloc
    use m_ES_WF_by_SDorCG, only : m_ESsd_dealloc
    use m_ES_wf_extrpl, only : m_ES_wf_extrpl_dealloc

!  call m_CtrlP_dealloc

    call m_CS_dealloc(neb_mode=.true.)
    call m_IS_dealloc(neb_mode=.true.)
    call m_pwBS_dealloc
    call m_Parallel_dealloc(neb_mode=.true.)
    call m_Kp_dealloc
    call m_Force_dealloc
    call m_CD_dealloc
    call m_XC_dealloc_vxc
    !call m_PP_dealloc
    !call m_NLP_dealloc
    call m_ES_dealloc
    call m_ESsd_dealloc
    call m_ES_wf_extrpl_dealloc()
    call m_CtrlP_set_init_status(.true.)
  end subroutine deallocate_array

  subroutine set_dimer_coordinates(i)
    use m_Ionic_System, only : pos, cps
    integer, intent(in) :: i
    real(kind=DP), pointer, dimension(:,:) :: r
    integer :: ia
    if(i==1) then
      r => dimvars%r1
    else
      r => dimvars%r2
    endif
    cps(:,:) = r(:,:)
    do ia = 1, dimvars%natm
       pos(ia,1:3) = matmul(transpose(rltv),r(ia,1:3))/PAI2
    end do
    if(mype == 0) then
      write(nfout,'(a,i5)') '!** dimer method coordinates for configuration ',i
      do ia=1, dimvars%natm
        write(nfout,'(i8,3f20.10)') ia,cps(ia,1:3)
      enddo
    endif
  end subroutine set_dimer_coordinates

  subroutine normalize(nsize,vec)
    integer, intent(in) :: nsize
    real(kind=DP), dimension(nsize,3), intent(inout) :: vec
    integer :: i
    real(kind=DP) :: v(3)
    vec = vec/(sqrt(sum(vec*vec)))
    return
  end subroutine normalize

  subroutine init_dimer_method()
    use m_Ionic_System, only : natm, cps_end0, cps_end1, &
            delta_theta, delta_r, dimer_max_iteration,  &
            imdtyp
    use m_Files, only : nfenf, nfdynm, F_DYNM, F_ENF, m_Files_open_files_initially
    use m_Control_Parameters, only : icond
    use m_Const_Parameters, only : CONTINUATION, INITIAL
    use m_Ionic_System, only : m_IS_wd_speciesname_etc, neb_convergence_threshold, &
            sw_random_unit_vector
    use m_Parallelization,      only : MPI_CommGroup, nrank_e, nrank_k, nrank_conf, mype_conf
    integer :: i
    call alloc_dimer(natm)
    do i=1, dimvars%natm
      dimvars%mobile(i) = imdtyp(i)
      if(mype==0) write(nfout,'(a,2i8)') '!** dimer method moblie',i,dimvars%mobile(i)
    enddo
    dimvars%dR = delta_r
    dimvars%dtheta = delta_theta
    dimvars%max_dimer_iteration = dimer_max_iteration
    dimvars%threshold = neb_convergence_threshold
    if(sw_random_unit_vector==ON) then
      call random_unit_vector(dimvars%natm, dimvars%N)
      cps_end1 = cps_end0-2*dimvars%dR*dimvars%N
      call set_r1_r2(natm, cps_end0, cps_end1)
    else
      call set_r1_r2(natm, cps_end0, cps_end1)
    endif
    call m_Files_open_files_initially()
    if (mype == 0) then
      call m_IS_wd_speciesname_etc(nfdynm)
    endif

    dimvars%nrank_r = nrank_conf
    dimvars%mype_r  = mype_conf

    call set_dimer_condition()

  end subroutine init_dimer_method

  subroutine random_unit_vector(natm, N)
    integer, intent(in) :: natm
    real(kind=DP), dimension(natm,3), intent(out) :: N
    integer :: i, j
    real(kind=DP) :: x
    call random_number(x)
    do j=1,3
      do i=1,natm
        if(dimvars%mobile(i)==ON) then
          call random_number(x)
          N(i,j) = x-0.5d0
        endif
      enddo
    enddo
    call normalize(natm,N)
    return
  end subroutine random_unit_vector

  subroutine finalize_dimer_method()
    use m_Files, only : nfdynm, nfenf, m_Files_close_all, m_Files_close_logfile
    use m_Control_parameters, only : terminated_because
    use m_Const_Parameters,   only : FORCE_CONVERGENCE_REACHED
    call dealloc_dimer()
    close(nfdynm)
    close(nfenf)
    call m_Files_close_all()
    terminated_because = FORCE_CONVERGENCE_REACHED
    call PrintStatus()
    call m_Files_close_logfile()
  end subroutine finalize_dimer_method

  subroutine alloc_dimer(natm)
    integer, intent(in) :: natm
    dimvars%natm = natm
    allocate(dimvars%r0(natm,3)) ; dimvars%r0=0.d0
    allocate(dimvars%r1(natm,3)) ; dimvars%r1=0.d0
    allocate(dimvars%r2(natm,3)) ; dimvars%r2=0.d0
    allocate(dimvars%mobile(natm)) ; dimvars%mobile=0
    allocate(dimvars%N(natm,3)) ; dimvars%N=0.d0
    allocate(dimvars%theta(natm,3)) ; dimvars%theta=0.d0
    allocate(dimvars%thetaold(natm,3)) ; dimvars%thetaold=0.d0
    allocate(dimvars%velocity(natm,3)) ; dimvars%velocity = 0.d0
    allocate(dimvars%f1perp(natm,3)) ; dimvars%f1perp = 0.d0
    allocate(dimvars%f1para(natm,3)) ; dimvars%f1para = 0.d0
    allocate(dimvars%f2perp(natm,3)) ; dimvars%f2perp = 0.d0
    allocate(dimvars%f2para(natm,3)) ; dimvars%f2para = 0.d0
    allocate(dimvars%fperp(natm,3)) ; dimvars%fperp = 0.d0
    allocate(dimvars%fperpold(natm,3)) ; dimvars%fperpold = 0.d0
    allocate(dimvars%fpara(natm,3)) ; dimvars%fpara = 0.d0
    allocate(dimvars%fR(natm,3)) ; dimvars%fR = 0.d0
    allocate(dimvars%f(natm,3)) ; dimvars%f = 0.d0
    allocate(dimvars%for(2,natm,3)) ; dimvars%for = 0.d0
    allocate(dimvars%fold(natm,3)) ; dimvars%fold = 0.d0
    allocate(dimvars%ftrans(natm,3)) ; dimvars%ftrans = 0.d0
    allocate(dimvars%currangle(natm)) ; dimvars%currangle = 0.d0
    allocate(dimvars%ene(2)) ; dimvars%ene = 0.d0
    dimvars%dR     = 0.01d0/BOHR
    dimvars%dtheta = 0.001d0
!    dimvars%dt     = 40.d0
    dimvars%dt     = neb_dt
    dimvars%max_dimer_iteration = 2000
    !pp_generated = .false.
    write(nfout,'(a,f15.5)') '!** dimer method dt',dimvars%dt
  end subroutine alloc_dimer

  subroutine set_r1_r2(natm, r1, r2)
    integer, intent(in) :: natm
    real(kind=DP), dimension(natm,3), intent(in) :: r1, r2
    integer :: i
    do i=1, dimvars%natm
      dimvars%N(i,:)  = r1(i,:)-r2(i,:)
      dimvars%r0(i,:) = 0.5d0*(r1(i,:)+r2(i,:))
    enddo
    call normalize(dimvars%natm, dimvars%N)
    if(mype==0) then
      do i=1, dimvars%natm
!        if (dimvars%mobile(i)==ON) then
          write(nfout,'(a,i10,3f10.5)') '!** dimer method vector N ',i,dimvars%N(i,1:3)
!        endif
      enddo
    endif
    call build_dimer()
  end subroutine set_r1_r2

  subroutine set_r0()
    integer :: i
    do i=1, dimvars%natm
      dimvars%r0(i,:) = 0.5d0*(dimvars%r1(i,:)+dimvars%r2(i,:))
    enddo
  end subroutine set_r0

  subroutine write_result(istep)
    use m_Files, only : nfdynm, nfenf
    integer, intent(in) :: istep
    integer :: ia
    if(mype == 0) then
      write(nfdynm,'(a,i10)') ' cps and forc at dimer iteration ',istep
      do ia=1, dimvars%natm
        write(nfdynm, '(i5, 3f15.9, 3f12.6)') ia, dimvars%r0(ia,1:3), dimvars%fR(ia,1:3)
      enddo
      write(nfenf, '(i5, f25.15, f20.10)') istep, dimvars%e0, dimvars%max_force
    endif
  end subroutine write_result

  subroutine build_dimer()
    integer :: i
    do i=1,dimvars%natm
      if(dimvars%mobile(i) == ON) then
        dimvars%r1(i,:) = dimvars%r0(i,:) + dimvars%dR * dimvars%N(i,:)
        dimvars%r2(i,:) = dimvars%r0(i,:) - dimvars%dR * dimvars%N(i,:)
      else
        dimvars%r1(i,:) = dimvars%r0(i,:)
        dimvars%r2(i,:) = dimvars%r0(i,:)
      endif
    enddo
  end subroutine build_dimer

  subroutine cal_curvature()
    integer :: i, j
    real(kind=DP) :: sm
    dimvars%curvature = 0.d0
    do i=1, dimvars%natm
!      if(dimvars%mobile(i)==ON) then
        dimvars%curvature = dimvars%curvature + dot_product(dimvars%for(2,i,:)-dimvars%for(1,i,:),dimvars%N(i,:))
!      endif
    enddo
    dimvars%curvature = 0.5d0*dimvars%curvature/dimvars%dR
    if(printable) then
      write(nfout,'(a,f15.5)') '!** dimer method curvature',dimvars%curvature
    endif
  end subroutine cal_curvature

  subroutine cal_dimer_energies()
    integer :: i
    real(kind=DP) :: fac
    dimvars%e  = dimvars%ene(1) + dimvars%ene(2)
    write(nfout,'(a,f15.5)') '!** dimer method dimer energy ',dimvars%e
    dimvars%e0 = dimvars%e*0.5d0
    fac = 0.25d0*dimvars%dR
    do i=1, dimvars%natm
!      if(dimvars%mobile(i)==ON) then
        dimvars%e0 = dimvars%e0 + fac * dot_product(dimvars%for(1,i,:)-dimvars%for(2,i,:), dimvars%N(i,:))
!      endif
    enddo
    if(printable) then
      write(nfout,'(a,f15.5)') '!** dimer method e0',dimvars%e0
    endif
  end subroutine cal_dimer_energies

  subroutine cal_dimer_forces()
    integer :: i
    do i=1, dimvars%natm
      dimvars%f1para(i,:) = dot_product(dimvars%for(1,i,:),dimvars%N(i,:))*dimvars%N(i,:)
      dimvars%f1perp(i,:) = dimvars%for(1,i,:) - dimvars%f1para(i,:)
      dimvars%f2para(i,:) = dot_product(dimvars%for(2,i,:),dimvars%N(i,:))*dimvars%N(i,:)
      dimvars%f2perp(i,:) = dimvars%for(2,i,:) - dimvars%f2para(i,:)
    enddo
    dimvars%fpara     = (dimvars%f1para + dimvars%f2para)*0.5d0
    dimvars%fperpold  = dimvars%fperp
    dimvars%fperp     = -(dimvars%f1perp - dimvars%f2perp)

    dimvars%fR(:,:) = (dimvars%for(1,:,:)+dimvars%for(2,:,:))*0.5d0

    dimvars%thetaold(:,:) = dimvars%theta(:,:)
    dimvars%theta(:,:)    = dimvars%fperp(:,:)
    call normalize(dimvars%natm, dimvars%theta)
    if(mype==0) then
      do i=1, dimvars%natm
!        if(dimvars%mobile(i)) then
          write(nfout,'(a,i8,3f20.10)') '!** dimer method theta',i,dimvars%theta(i,1:3)
!        endif
      enddo
    endif

    dimvars%fold(:,:) = dimvars%f(:,:)
    dimvars%f(:,:)    = dimvars%for(1,:,:)-dimvars%for(2,:,:)
  end subroutine cal_dimer_forces

  subroutine rotate_dimer(natm, dtheta)
    integer, intent(in) :: natm
    real(kind=DP), intent(in) :: dtheta
    integer :: i, j
    real(kind=DP) :: rot(3)
    do i=1, dimvars%natm
      if(dimvars%mobile(i) == ON) then
        rot(:) = (dimvars%N(i,:)     * cos(dtheta) &
        &      +  dimvars%theta(i,:) * sin(dtheta)) * dimvars%dR
        dimvars%r1(i,:) = dimvars%r0(i,:) + rot(:)
        dimvars%r2(i,:) = dimvars%r0(i,:) - rot(:)
      else
        dimvars%r1(i,:) = dimvars%r0(i,:)
        dimvars%r2(i,:) = dimvars%r0(i,:)
      endif
    enddo
  end subroutine rotate_dimer

  subroutine update_N()
    integer :: i
    do i=1, dimvars%natm
      dimvars%N(i,:) = dimvars%r1(i,:)-dimvars%r2(i,:)
    enddo
    call normalize(dimvars%natm, dimvars%N)
    if(printable) then
      do i=1, dimvars%natm
        write(nfout,'(a,i10,3f10.5)') '!** dimer method vector N ',i,dimvars%N(i,1:3)
      enddo
    endif
  end subroutine update_N

  subroutine translate_dimer(itr)
    use m_Ionic_System, only : ionic_mass, dimer_time_integral, imdtyp, m_IS_fire_core
    use m_Const_Parameters, only : QUENCHED_MD, FIRE, DP
    integer, intent(in) :: itr
    integer :: i, j
    real(kind=DP), dimension(3) :: ftrans, ftrans0
    real(kind=DP) :: dt, mass
    real(8), dimension(3) :: r1, vel, unitf, v1, v
    real(8) :: fmax, fabs
    real(kind=DP), allocatable, dimension(:,:) :: cpst,forct,cpdt,transfor
    integer :: ndim
    dt = dimvars%dt
    fmax = 0.d0
    allocate(transfor(dimvars%natm,3));transfor=0.d0
    do i=1, dimvars%natm
      if(dimvars%mobile(i) /= ON) cycle
      ftrans0(:) = dimvars%ftrans(i,:)
      if(dimvars%curvature>0) then
        transfor(i,:) = -dimvars%fpara(i,:)
      else
        transfor(i,:) =  dimvars%fR(i,:)-2*dimvars%fpara(i,:)
      endif
    enddo
    if (dimer_time_integral == QUENCHED_MD) then
      do i=1, dimvars%natm
        if(dimvars%mobile(i) /= ON) cycle
        ftrans0(:) = dimvars%ftrans(i,:)
        ftrans(:)  = transfor(i,:)
        unitf=0.d0
        if(dot_product(ftrans,ftrans)/=0) unitf=ftrans/dsqrt(dot_product(ftrans,ftrans))
        v = dimvars%velocity(i,:)
        mass = ionic_mass(i)
        v1 = v + 0.5d0*dt*(ftrans0+ftrans)/mass
        v = 0.d0
        if( dot_product(v1,unitf) >= 0.0d0) then
          v = unitf * dot_product(v1,unitf)
        else
          if (mype==0) write(nfout,'(a,i6)') 'quenched velocity of atom : ',i
        end if
        dimvars%r1(i,:) = dimvars%r1(i,:) + v(:)*dt + ftrans(:)*dt**2/mass
        dimvars%r2(i,:) = dimvars%r1(i,:) - 2*dimvars%N(i,:)*dimvars%dR
        dimvars%ftrans(i,:) = ftrans(:)
        dimvars%velocity(i,:) = v
      enddo
    else if (dimer_time_integral == FIRE) then
      ndim = dimvars%natm
      allocate(cpst(ndim,3));  cpst  = dimvars%r1
      allocate(cpdt(ndim,3));  cpdt  = dimvars%velocity
      allocate(forct(ndim,3)); forct = transfor
      call m_IS_fire_core(nfout, ndim, itr, cpst, cpdt, forct, dimvars%mobile)
      dimvars%r1       = cpst
      dimvars%r2       = dimvars%r1 - 2*dimvars%N*dimvars%dR
      dimvars%ftrans   = transfor
      dimvars%velocity = cpdt
      deallocate(cpst)
      deallocate(cpdt)
      deallocate(forct)
    else
      write(nfout,*) 'unsupported dimer_time_integral ', dimer_time_integral
      !stop
      call phase_error_with_msg(nfout,'unsupported dimer_time_integral',__LINE__,__FILE__)
    endif

    fmax = 0.d0
    do i=1, dimvars%natm
      fabs = sqrt(dot_product(transfor(i,:),transfor(i,:)))
      if (fabs>fmax) fmax = fabs
    enddo
    dimvars%max_force = fmax

    deallocate(transfor)
    if(mype==0) write(nfout, '(a,f15.8)') '!** dimer method max translational force ',fmax

  end subroutine translate_dimer

  subroutine cal_dtheta(natm, deltheta)
    integer, intent(in) :: natm
    real(kind=8), intent(out) :: deltheta
    real(kind=8) :: fscalar, fscalarp, invdtheta, maxfscalar, thetat, fscalarperp, maxperp, perp
    integer :: i
    invdtheta = 1.d0/dimvars%dtheta
    deltheta = 0.d0
    maxfscalar = 0.d0
    fscalar  = 0.d0
    fscalarperp  = 0.d0
    fscalarp = 0.d0
    thetat = 0.d0
    maxperp = 0.d0
    do i=1, dimvars%natm
      fscalar = fscalar + dot_product(dimvars%fperpold(i,:), dimvars%thetaold(i,:)) &
                        + dot_product(dimvars%fperp(i,:), dimvars%theta(i,:))
      fscalarp = fscalarp + (dot_product(dimvars%f(i,:),    dimvars%theta(i,:)) -   &
                             dot_product(dimvars%fold(i,:), dimvars%thetaold(i,:))) &
               * invdtheta
!      fscalarperp = fscalarperp + dot_product(dimvars%fperpold(i,:), dimvars%thetaold(i,:)) / &
!                dimvars%dR
!      thetat = thetat + (dot_product(dimvars%f(i,:),    dimvars%theta(i,:)) +   &
!                         dot_product(dimvars%fold(i,:), dimvars%thetaold(i,:)))
    enddo
    !deltheta = -0.5d0*atan(2.d0*dimvars%dR*fscalar/abs(fscalarp))
    deltheta = -0.5d0*atan(fscalar/abs(fscalarp))
  end subroutine cal_dtheta

  subroutine dealloc_dimer()
    deallocate(dimvars%r0)
    deallocate(dimvars%r1)
    deallocate(dimvars%r2)
    deallocate(dimvars%mobile)
    deallocate(dimvars%N)
    deallocate(dimvars%theta)
    deallocate(dimvars%thetaold)
    deallocate(dimvars%velocity)
    deallocate(dimvars%f1perp)
    deallocate(dimvars%f1para)
    deallocate(dimvars%f2perp)
    deallocate(dimvars%f2para)
    deallocate(dimvars%fperp)
    deallocate(dimvars%fpara)
    deallocate(dimvars%fR)
    deallocate(dimvars%f)
    deallocate(dimvars%fold)
    deallocate(dimvars%ftrans)
  end subroutine dealloc_dimer

  subroutine m_dm_set_standard_out_opened(tf)
    logical, intent(in) :: tf
    standard_out_opened = tf
  end subroutine m_dm_set_standard_out_opened

  subroutine set_dimer_condition()
    use m_Files, only : F_OUT, F_CHGT, F_CNTN, F_CNTN_BIN, F_ZAJ, F_CHR &
            , F_CNTN_BIN_PAW, F_OCCMAT
    if(dimer_cond_set) return
    dimvars%file_type%F_OUT = trim(F_OUT)
    dimvars%file_type%F_CHGT = trim(F_CHGT)
    dimvars%file_type%F_CNTN = trim(F_CNTN)
    dimvars%file_type%F_CNTN_BIN = trim(F_CNTN_BIN)
    dimvars%file_type%F_ZAJ = trim(F_ZAJ)
    dimvars%file_type%F_CHR = trim(F_CHR)

! ======================== KT_add ================= 13.0D
    dimvars%file_type%F_CNTN_BIN_PAW = trim(F_CNTN_BIN_PAW)
! ================================================= 13.0D
    dimvars%file_type%F_OCCMAT = F_OCCMAT

    dimer_cond_set = .true.

  end subroutine set_dimer_condition

  subroutine set_dimer_replica_parameter(image_id, istep)
    use m_Const_Parameters, only : ON
    use m_Control_parameters, only : icond, ipriparadeb, sw_hubbard
    use m_IterationNumbers, only : iteration,iteration_electronic,iteration_ionic
    use m_Replica
    use m_Files
    use m_Parallelization,  only : mype

    integer, intent(in) :: image_id, istep
    character(4) image_name

    if( istep > 1 .and. icond == 0) then
       icond = 1
    end if

    write(image_name,'(a1,i3.3)') 'r', image_id

!    F_OUT = trim(dimvars%file_type%F_OUT)//'_'//image_name
    F_OUT = trim(dimvars%file_type%F_OUT)

! ===== KT_add === 2014/07/20
#ifdef NEB_NEW_FILENAMES
!  F_OUT = trim(F_OUT_BASE)//'_'//image_name &
!       &              //trim(F_CONF_EXTENSION)//trim(F_PARA_EXTENSION)
!    F_OUT = trim(adjustl(F_OUT_BASE)) // '_' //image_name &
!         &                             // trim(adjustl(F_PARA_EXTENSION))
    F_OUT = trim(adjustl(F_OUT_BASE)) // trim(adjustl(F_PARA_EXTENSION))
! ================ 2014/07/20
#endif

    if(.not.standard_out_opened) then
      call m_Files_open_standardout()
      standard_out_opened = .true.
    endif

  !F_DYNM = trim(dimvars%file_type%F_DYNM)//'_'//image_name
  !F_ENF = trim(dimvars%file_type%F_ENF)//'_'//image_name
    F_CHGT = trim(dimvars%file_type%F_CHGT)//'_'//image_name
    F_CNTN = trim(dimvars%file_type%F_CNTN)//'_'//image_name
    F_CNTN_BIN = trim(dimvars%file_type%F_CNTN_BIN)//'_'//image_name
    F_ZAJ = trim(dimvars%file_type%F_ZAJ)//'_'//image_name
    F_CHR = trim(dimvars%file_type%F_CHR)//'_'//image_name

! ======================== KT_add ================= 13.0D
    F_CNTN_BIN_PAW = trim(dimvars%file_type%F_CNTN_BIN_PAW)//'_'//image_name
! ================================================= 13.0D
    F_OCCMAT = trim(dimvars%file_type%F_OCCMAT)//'_'//image_name

    F_ZAJ_in      = F_ZAJ
    F_CHGT_in     = F_CHGT
    F_CNTN_BIN_in = F_CNTN_BIN
    F_CNTN_in     = F_CNTN

! ======================== KT_add ================= 13.0D
    F_CNTN_BIN_PAW_in = F_CNTN_BIN_PAW
! ================================================= 13.0D

    close(nfcntn)
    close(nfcntn_bin)
    open(nfcntn,file=trim(F_CNTN))
    open(nfcntn_bin,file=trim(F_CNTN_BIN),form='unformatted')

! ======================== KT_add ================= 13.0D
    if ( flg_paw ) then
       call m_Files_close_nfcntn_bin_paw()
       call m_Files_open_nfcntn_bin_paw()
    endif
! ================================================= 13.0D

    iteration = 0
    iteration_electronic = 0
    iteration_ionic = 0
  end subroutine set_dimer_replica_parameter

  subroutine m_dm_do_dimer_method()
    use m_Const_Parameters, only : DP

    implicit none

    integer :: idimer, i
    real(kind=DP) :: maxf
    real(kind=DP) :: deltheta

    call init_dimer_method()

    !calculate energy and force for the initial configuration
    call cal_dimer_scf(1)

    call cal_dimer_energies()
    call cal_dimer_forces()

    do idimer=1, dimvars%max_dimer_iteration
!  do idimer=1, 40
!    call cal_curvature()
      call cal_dimer_energies()
      call cal_dimer_forces()
      deltheta = dimvars%dtheta
      call rotate_dimer(dimvars%natm, deltheta)

      !calculate energy and force after test rotation
      call cal_dimer_scf(idimer+1)

      call cal_curvature()
      call cal_dimer_energies()
      call cal_dimer_forces()
      call cal_dtheta(dimvars%natm, deltheta)
      call rotate_dimer(dimvars%natm, deltheta)
      call update_N()

      !calculate energy and force after rotation
      call cal_dimer_scf(idimer+1)
      call translate_dimer(idimer)
      call set_r0()
      call write_result(idimer)
      if(dimvars%max_force<dimvars%threshold) then
        if(printable) write(nfout,'(a,f15.8)') &
        '!** dimer method reached convergence. max translational force : ',dimvars%max_force
        exit
      endif
    enddo
    call finalize_dimer_method()
  end subroutine m_dm_do_dimer_method

end module m_dimer

subroutine do_dimer_method()
  use m_dimer, only : m_dm_do_dimer_method
  implicit none
  call m_dm_do_dimer_method()
end subroutine do_dimer_method
