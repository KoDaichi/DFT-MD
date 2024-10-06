module m_Fcp
  use m_Const_Parameters,   only : DP,T_CONTROL,VELOCITY_SCALING,CONST_kB, &
     STEEPEST_DESCENT,VERLET,OFF,QUENCHED_MD,BFGS,GDIIS,ANEW,RENEW,ON
  use m_Control_Parameters, only : fcp_mass,fcp_temperature,fcp_mu,fcp_qmass, &
     iprifcp,printable,dtio,fcp_relax_crit,imdalg,mode_fi_coefficient,kqnmditer_p, &
     c_forc_prop_region_high,c_forc_prop_region_low,factor_prop_region,gdiis_hownew, &
     eigenvalue_threshold,sw_optimize_alpha,sw_correct_eigenvalue
  use m_Files,              only : nfout
  use m_IterationNumbers,   only : iteration_ionic
  use m_Ionic_System,       only : t_ctrl_method
  use m_Electronic_Structure, only : cps => totch, efermi
  implicit none
  private

  ! md_thermo
  real(kind=DP) :: cpd_l

  real(kind=DP) :: cpqr(2)
  real(kind=DP) :: forcp

  real(kind=DP) :: ekr
  real(kind=DP) :: tkb,cprv,frsv
  real(kind=DP) :: ekina,ekinq,ekbt,ega

  logical :: lfirst = .true.

  ! CG2
  integer, save :: iter_CG = 1, iter_linmin = 1, iter_CG_max = 1

  ! md
  real(kind=DP) :: cpd_old
  real(kind=DP) :: fi_coefficient = 1.d0

  ! GDIIS, BFGS
  real(kind=DP),allocatable,dimension(:)   :: u_l ! d(kqnmditer_p)
  real(kind=DP),allocatable,dimension(:)   :: w_l ! d(kqnmditer_p)
  real(kind=DP)                            :: fc_l
  integer      ,allocatable,dimension(:)   :: ncrspd ! d(kqnmditer_p)
  real(kind=DP),allocatable,dimension(:,:) :: f_gdiis ! d(kqnmditer_p,kqnmditer_p)
  real(kind=DP),allocatable,dimension(:,:) :: f_wk ! d(kqnmditer_p**2,2)
  real(kind=DP),allocatable,dimension(:)   :: f_rslv ! d(kqnmditer_p**2)
  real(kind=DP),allocatable,dimension(:)   :: g,e_wk,ww1,etot_trial
  real(kind=DP)                            :: forc_g
  integer, allocatable, dimension(:)       :: ip ! d(kqnmditer_p)
  integer, parameter                       :: UNIT = 1
  integer, save                            :: iter_gdiis = 0
  integer, save                            :: if_allocated = 0
  real(kind=DP) :: min_alpha = 1.d0,max_alpha = 2.0d0
  real(kind=DP) :: wolfe_c1=1.d-4,wolfe_c2=0.9

  ! NEB
  integer :: num_of_images
  type image_type
     real(kind=DP) :: cps      = 0.0d0
     real(kind=DP) :: cps0     = 0.0d0
     real(kind=DP) :: force    = 0.0d0
     real(kind=DP) :: force0   = 0.0d0
     real(kind=DP) :: velocity = 0.0d0
  end type image_type
  type(image_type), allocatable :: neb_image(:)
  logical, save :: cps_initialized = .false.

  public :: m_Fcp_md_thermo, m_Fcp_cg, m_Fcp_CG_reset, m_Fcp_cg2, &
     m_Fcp_Converged, m_Fcp_md, m_Fcp_gdiis, m_Fcp_gdiis_reset, &
     m_Fcp_create_replica, m_Fcp_set_neb_replica_totch, &
     m_Fcp_set_neb_replica_force, m_Fcp_update_neb, m_Fcp_replica_Converged, &
     m_Fcp_set_neb_parameter, m_Fcp_allreduce_neb_force, &
     m_Fcp_write_down_totch_replica, m_Fcp_read_totch_replica, &
     m_Fcp_initialize_neb_totch, m_Fcp_print_status

contains

!  Nose-Hoover, Velocity Scaling

  subroutine initialize_fcp()
    implicit none
    cpqr  = 0.0d0
    forcp = 0.0d0
    cprv  = 0.0d0
    frsv  = 0.0d0
    tkb   = fcp_temperature * CONST_kB
    cpd_l = dsqrt( tkb / fcp_mass )
  end subroutine initialize_fcp

  subroutine m_Fcp_md_thermo()
    implicit none
    real(kind=DP) :: forc_l
    integer       :: mdalg

    if( lfirst ) then
       call initialize_fcp()
       lfirst = .false.
    end if

    forc_l = fcp_mu - efermi

    mdalg = T_CONTROL

    ! --> Velocities at iteration_ionic-th step
    if(iteration_ionic == 1) then
       call forcrsv(cpd_l,tkb,frsv)
       !  -(b_Ionic_System)(cpd_l,tkb)->(frsv(=force on the thermostat coordinate))
    else
        call vlcty_accrd2_vVerlet(mdalg,forc_l)    !-(m_Ionic_System) ->(cpd_l)
    end if
!!$ 2011.06.06
    if(t_ctrl_method.ne.VELOCITY_SCALING) call heatrsv()
!!$ 2011.06.06
         ! natm-(ifix_P_ifree):= #atoms in a heat bath
         ! -(m_Ionic_System) (cpd_l,cpqr,frsv,tkb)->(cpd_l,cpqr,frsv)
    forcp = forc_l
    call ekina_ekinq_ekbt_and_ega(mdalg)          !-(m_Ionic_System) ->(ekina,..)
    !     <== Kinetic and thermostat energies at iteration_ionic-th step
    call evolve_crdn_ACCRD2_vVerlet(forc_l)       !-(m_Ionic_System) ->(cps)
    !     <== Coordinates at (iteration_ionic+1)-th step
    call evolve_cprv                                    !-(m_Ionic_System) ->(cprv)

    if(t_ctrl_method == VELOCITY_SCALING)then
       call scale_velocity()
    endif

    !if(iprifcp >= 1) then
    !   write(nfout, '(" FCP : Fermi Energy = ",F20.8)') efermi
    !   write(nfout, '(" FCP : Target Mu    = ",F20.8)') fcp_mu
    !   write(nfout, '(" FCP : Total Charge = ",F20.8)') cps
    !end if

  end subroutine m_Fcp_md_thermo

  subroutine m_Fcp_print_status()
    if(iprifcp >= 1) then
       write(nfout, '("FCP : Fermi Energy = ",F20.8)') efermi
       write(nfout, '("FCP : Target Mu    = ",F20.8)') fcp_mu
       write(nfout, '("FCP : Total Charge = ",F20.8)') cps
    end if
  end subroutine m_Fcp_print_status

  subroutine forcrsv(cpd_l,tkb,frsv)
    implicit none
    real(kind=DP),intent(in)  :: cpd_l
    real(kind=DP),intent(in)  :: tkb
    real(kind=DP),intent(out) :: frsv

    frsv = cpd_l*cpd_l*fcp_mass-3.d0*tkb
    return
  end subroutine forcrsv

  subroutine vlcty_accrd2_vVerlet(mdalg,forc_l)
    implicit none
    integer      , intent(in) :: mdalg
    real(kind=DP), intent(in) :: forc_l
    real(kind=DP)        :: cpdxyz,frcxyz
    real(kind=DP)        :: dtz

    cpdxyz = cpd_l
    frcxyz = (forc_l + forcp)/2.d0
    cpd_l  = cpdxyz + dtio/fcp_mass*frcxyz

    if(mdalg == T_CONTROL) then
       dtz = dtio*cpqr(1)/2.d0
       cpd_l = cpd_l - dtz*cpdxyz
    end if

  end subroutine vlcty_accrd2_vVerlet

  subroutine heatrsv()
    ! OUTPUT : cpd_l,cpqr,frsv
    implicit none

    real(kind=DP) :: cpdcor,cpqrpre,cpqrcor,fchk

    integer,       parameter :: mmin = 5, mmax = 100
    real(kind=DP), parameter :: eps = 1.d-12
    integer                  :: ic
    real(kind=DP)            :: dtz

    cpdcor = cpd_l

    cpqrcor = cpqr(2) + 2*dtio/fcp_qmass*frsv
    cpqr(2)  = cpqr(1)
    cpqr(1)  = cpqr(1) + 0.5*dtio/fcp_qmass*frsv

! +++++ predictor-corrector method +++++
    PREDICTOR_CORRECTOR: do ic = 1, mmax
       cpqrpre= cpqrcor
       dtz = 1 + dtio*0.5*cpqrpre
       cpdcor= cpd_l/dtz

       call forcrsv(cpdcor,tkb,frsv)         ! -(b_Ionic_System) ->(frsv)

       cpqrcor= cpqr(1)+frsv*dtio/fcp_qmass*0.5

       if(ic <= mmin) cycle

       fchk = dabs(cpqrpre/cpqrcor - 1.d0)
       if(fchk >= eps) cycle PREDICTOR_CORRECTOR

       go to 300
    end do PREDICTOR_CORRECTOR

    if(printable) then
       write(nfout,'(" ! Warning in fcp_heatrsv ***")')
       write(nfout,'(" FCP cpqrpre, cpqrcor = ",2d20.12)') &
          cpqrpre,cpqrcor
    end if
    call phase_error_with_msg(nfout,'fcp_heatrsv',__LINE__,__FILE__)

300 continue
    if(printable) then
       write(nfout,'(" FCP ic, cpqrcor, fchk= ",i1,d12.4,d12.4)') ic,cpqrcor,fchk
    end if

! ++++++++++++++++++++++++++++++++++++++++++++++++++

    cpd_l   = cpdcor
    cpqr(1) = cpqrcor

  end subroutine heatrsv

  subroutine ekina_ekinq_ekbt_and_ega(mdalg)
    implicit none
    integer,intent(in) :: mdalg
    real(kind=DP)      :: tkin
    ekr   = 0.d0   !d(nrsv)
    ekina = 0.d0
    ekbt  = 0.d0
    ekinq = 0.d0

    tkin  = cpd_l*cpd_l*fcp_mass*0.5d0
    ekina = ekina + tkin

    ekr   = ekr + tkin
    ekbt  = ekbt + 3*tkb*cprv

    ekinq = ekinq + fcp_qmass/2.d0 * cpqr(1)*cpqr(1)
    if(mdalg == T_CONTROL) then
       ekr = ekr/1.5d0
    end if

    ega = ekina + ekinq + ekbt

    if(printable) then
       write(nfout,'(" ***** ",i6,"-th time step ; ega= ",1pe16.7 /&
            &," ekina, ekinq, ekbt for FCP",3e16.7)') &
            & iteration_ionic, ega, ekina, ekinq, ekbt
       write(nfout,'(" tkb, ekr= ",2e16.7)') tkb,ekr
    end if

    if(iprifcp >= 2 .and. printable) then
       write(nfout,'(" ! *** cps for FCP ***")')
       write(nfout,'(" ",1f20.10)') cps
       write(nfout,'(" ! *** cprv, cpqr for FCP ***")')
       write(nfout,'(" ",2f20.10)') cprv,cpqr(1)
    end if
  end subroutine ekina_ekinq_ekbt_and_ega

  subroutine evolve_crdn_ACCRD2_vVerlet(forc_l)
    implicit none
    real(kind=DP),intent(in) :: forc_l
    real(kind=DP) :: dtm, dtz

!---*----*----*----*----*----> Velocity Verlet
    dtm = dtio/fcp_mass/2.d0
    cps = cps+dtio*(cpd_l+dtm*forc_l)
!---*----*----*----*----*----< Velocity Verlet
    dtz= dtio*dtio*cpqr(1)/2.d0
    cps= cps - dtz*cpd_l

  end subroutine evolve_crdn_ACCRD2_vVerlet

  subroutine evolve_cprv
    implicit none
    cprv =  cprv+dtio*(cpqr(1) + 0.5d0*dtio/fcp_qmass*frsv)
    if(printable) write(nfout,'(" FCP cprv = ",1d12.4)') cprv
  end subroutine evolve_cprv

  subroutine scale_velocity()
    implicit none
    real(kind=DP) :: tkin
    integer :: imdt

    tkin = cpd_l*cpd_l * fcp_mass*0.5d0

    if(tkin.gt.1e-12) then
       cpd_l = cpd_l * dsqrt(0.5d0*tkb/tkin)
    end if

    if(iprifcp >= 1) then
       write(nfout,'(a)') 'scaled the FCP velocities'
       write(nfout,'(" !! tkin = ",d20.8)') tkin
       write(nfout,'(" !! tkb  = ",d20.8)') tkb
       write(nfout,'(" !! dsqrt(0.5d0*tkb/tkin) = ",d20.8)') dsqrt(0.5d0*tkb/tkin)
    end if
  end subroutine scale_velocity

! CG

  subroutine m_Fcp_cg()
    implicit none
    real(kind=DP) :: forc_l
    real(kind=DP) :: etotal

    logical, save :: finit = .true., fconv = .false.
    integer, save :: itime = 1
    real(kind=DP), save :: e0,e1,e2,e3
    real(kind=DP), save :: dt1,dt2,dt3
    real(kind=DP), save :: de1,de2,de3
    real(kind=DP), save :: hvec,g0vec,g1vec,cps0

    if(printable) write(nfout,*) 'FCP gCG: iter=',iteration_ionic,' itime=',itime

    forc_l = fcp_mu - efermi
    etotal = forc_l ** 2

    if(finit) then
       hvec  = forc_l
       g0vec = hvec
       finit = .false.
    end if

    if(itime == 1) then
       call first()
       itime = 2
    else if(itime == 2) then
       e2 = etotal
       call deriv_energy(de2,hvec)
       if(de2 > 0.d0 .or. e2 > e0 ) then
          call estimate_min(dt1,dt2,e1,e2,de1,de2,fconv)
          itime = 3
       else
      !debug
          if(printable) write(nfout,*) 'FCP gCG:== Second (dt2)=='
      !end debug
          e1 = e2
          de1 = de2
          dt1 = dt2
          dt2 = 2*dt1
          cpd_l = dt2 * hvec
          cps = cps0 + cpd_l
      !debug
      ! if(printable)
      ! write(nfout,*) 'gCG:cpd_l:'
      ! do ia=1,natm
      !    write(nfout,'("gCG:",i3,3(1x,f17.9))') ia,cpd_l(ia,1:3)
      ! end do
      ! write(nfout,*) 'gCG:cps:'
      ! do ia=1,natm
      !    write(nfout,'("gCG:",i3,3(1x,f17.9))') ia,cps(ia,1:3)
      ! end do
      ! end if
      !end debug
       end if
    else
       call estimate_min(dt1,dt2,e1,e2,de1,de2,fconv)
       itime = itime + 1
       if(fconv) then
          g1vec = forc_l
          call conjugate_grad(hvec,g0vec,g1vec)
          call first()
          itime = 2
          fconv = .false.
       end if
    end if

  contains

    subroutine first()
      implicit none
      !debug
           if(printable) write(nfout,*) 'FCP gCG:== First (dt1)=='
      !end debug
      cps0 = cps
      e0 = etotal
      dt1 = 0.d0
      e1 = e0
      call deriv_energy(de1,hvec)
      dt2 = 1.d0
      cpd_l = dt2 * hvec

      cps = cps0 + cpd_l
      !debug
      ! if(printable)
      ! write(nfout,*) 'gCG:cpd_l:'
      ! do ia=1,natm
      !    write(nfout,'("gCG:",i3,3(1x,f17.9))') ia,cpd_l(ia,1:3)
      ! end do
      ! write(nfout,*) 'gCG:cps:'
      ! do ia=1,natm
      !    write(nfout,'("gCG:",i3,3(1x,f17.9))') ia,cps(ia,1:3)
      ! end do
      ! end if
      !end debug
    end subroutine first

    subroutine deriv_energy(de,hvec)
      implicit none
      real(kind=DP), intent(out) :: de
      real(kind=DP), intent(in) :: hvec

      de = - forc_l*hvec

    end subroutine deriv_energy

    subroutine estimate_min(dt1,dt2,e1,e2,de1,de2,fconv)
      implicit none
      real(kind=DP), intent(inout) :: dt1,dt2,e1,e2,de1,de2
      logical, intent(out) :: fconv

      real(kind=DP), save :: dt3,e3,de3
      real(kind=DP) :: c1,c2,c3,c4
      real(kind=DP) :: ddt,sdt,pdt2,pdt3,sq,dtm,dtp
      real(kind=DP) :: dfp,dfm,fp,fm
      logical, save :: finit = .true.

      if(.not.finit) then
         e3 = etotal
         call deriv_energy(de3,hvec)
         !debug
         if(printable) then
           write(nfout,*) 'FCP gCG:== Estimate Emin (old dt1,dt2) =='
           write(nfout,*) 'FCP gCG:dt1=',dt1,' e1=',e1
           write(nfout,*) 'FCP gCG:de1=',de1
           write(nfout,*) 'FCP gCG:dt2=',dt2,' e2=',e2
           write(nfout,*) 'FCP gCG:de2=',de2
           write(nfout,*) 'FCP gCG:dt3=',dt3,' e3=',e3
           write(nfout,*) 'FCP gCG:de3=',de3
!           write(nfout,*) 'FCP gCG:etol=',etol
         end if
         !end debug
         if(abs(forc_l) .le. fcp_relax_crit) then
            fconv = .true.
            finit = .true.
            return
         end if
         if(de3 > 0.d0 .or. e3 > e1) then
            e2 = e3
            de2 = de3
            dt2 = dt3
         else
            e1 = e3
            de1 = de3
            dt1 = dt3
         end if
      end if
      finit = .false.
      !debug
      if(printable) then
        write(nfout,*) 'FCP gCG:== Estimate Emin (new dt1,dt2) =='
        write(nfout,*) 'FCP gCG:dt1=',dt1,' dt2=',dt2
        write(nfout,*) 'FCP gCG: e1=',e1,' e2=',e2
        write(nfout,*) 'FCP gCG:de1=',de1,' de2=',de2
      end if
      !end debug

      ddt = dt1 - dt2
      sdt = dt1 + dt2
      pdt2 = dt1*dt1
      pdt3 = dt1*pdt2
      c1 = ((de1+de2)*ddt-2.d0*(e1-e2))/ddt**3
      c2 = 0.5d0*(de1-de2-3.d0*c1*ddt*sdt)/ddt
      c3 = de1 - 3.d0*c1*pdt2 - 2.d0*c2*dt1
      c4 = e1 - c1*pdt3 - c2*pdt2 - c3*dt1

!!$ASASASASASAS
!!$      sq = c2+sqrt(c2*c2-3.d0*c1*c3)
      sq = c2*c2-3.d0*c1*c3
      if (sq > 0.d0) then
         sq = sqrt(sq)
      else
         write(nfout,*) 'FCP gCG: Warning! Cannot estimate the minimum.'
         sq = 0.d0
      endif
      sq = sq + c2
!!$ASASASASASAS
      dtp = -c3/sq
      dtm = -sq/(3.d0*c1)

      dfp = 3.d0*c1*dtp**2+2.d0*c2*dtp+c3
      dfm = 3.d0*c1*dtm**2+2.d0*c2*dtm+c3

      fp = c1*dtp**3+c2*dtp**2+c3*dtp+c4
      fm = c1*dtm**3+c2*dtm**2+c3*dtm+c4

      !debug
      if(printable) then
        write(nfout,*) 'FCP gCG:== Estimate Emin (c1,c2,c3) =='
        write(nfout,*) 'FCP gCG: c1=',c1,' c2=',c2
        write(nfout,*) 'FCP gCG: c3=',c3,' c4=',c4
        write(nfout,*) 'FCP gCG:dtp=',dtp,' dtm=',dtm
        write(nfout,*) 'FCP gCG:dfp=',dfp,' dfm=',dfm
        write(nfout,*) 'FCP gCG: fp=',fp,' fm=',fm
      end if
      !end debug

      if(dtp >= dt1 .and. dtp <= dt2) then
         dt3 = dtp
      else if(dtm >= dt1 .and. dtm <= dt2) then
         dt3 = dtm
      else
        !!stop 'I cant estimated dt3 at which the total energy is minimum'
        write(nfout,*) 'FCP gCG: Warning! line-minimization failed'
        dt3 = 1.d0
      end if

      cpd_l = dt3 * hvec
      cps = cps0 + cpd_l

      !debug
      if(printable) then
       write(nfout,*) 'FCP gCG:== Estimate Emin (dt3) =='
       write(nfout,*) 'FCP gCG:dt1=',dt1,' dt2=',dt2
       write(nfout,*) 'FCP gCG:dtp=',dtp
       write(nfout,*) 'FCP gCG:dtm=',dtm
       write(nfout,*) 'FCP gCG:dt3=',dt3
      ! write(nfout,*) 'FCP gCG:cpd_l:'
      ! write(nfout,'("FCP gCG:",1x,f17.9)') cpd_l
      ! write(nfout,*) 'FCP gCG:cps:'
      ! write(nfout,'("FCP gCG:",1x,f17.9)') cps
      end if
      !end debug

    end subroutine estimate_min

    subroutine conjugate_grad(hvec,g0vec,g1vec)
      implicit none
      real(kind=DP), intent(inout) :: hvec,g0vec
      real(kind=DP), intent(in) :: g1vec

      real(kind=DP) :: gg,gam

      gg = g0vec**2
      gam = (g1vec-g0vec)*g1vec
      if(gg > 1.d-10*gam) then
         gam = gam/gg
      else
         gam = 1.d0
      end if

      hvec  = g1vec + gam * hvec
      g0vec = g1vec

      !debug
      if(printable) then
       write(nfout,*) 'FCP gCG:== CG =='
       write(nfout,*) 'FCP gCG:gam=',gam
       write(nfout,*) 'FCP gCG:gg=',gg
       write(nfout,*) 'FCP gCG: H (conjugate_grad)'
       write(nfout,'("FCP gCG:",1x,f17.9)') hvec
       write(nfout,*) 'FCP gCG: G'
       write(nfout,'("FCP gCG:",1x,f17.9)') g0vec
      end if
      !end debug

    end subroutine conjugate_grad

  end subroutine m_Fcp_cg

! CG2

  subroutine m_Fcp_CG_reset()
    implicit none
    iter_CG = 1
    iter_linmin = 1
    iter_CG_max = 1
  end subroutine m_Fcp_CG_reset

  subroutine m_Fcp_cg2()
    implicit none
    real(kind=DP) :: forc_l

    real(kind=DP), save :: alpha, alpha2, gamma
    real(kind=DP), save :: f_para1, f_para2, f_para3
    real(kind=DP), save :: vec_h_norm
    real(kind=DP), save :: vec_g, vec_h
    real(kind=DP), save ::       f_total0
    real(kind=DP), save :: cps1, f_total1
    real(kind=DP), save :: cps2, f_total2

    forc_l = fcp_mu - efermi

    f_total0 = forc_l

    if ( ( iter_CG .eq. 1 ) .and. ( iter_linmin .eq. 1 ) ) then
       call first_cg_step()
    else if ( iter_linmin .eq. 2 ) then
       f_para2 = f_total0*sign(1.0d0,vec_h)
       if(printable) write(nfout,'(a18,2f10.6)') ' FCP CG2: F_parallel =', f_para1, f_para2
       if ( abs(f_para2) .lt. f_para1*0.2d0 ) then
          if(printable) write(nfout,*) 'FCP CG2: F_parallel is small. Then linmin finished.'
          if ( iter_CG .lt. iter_CG_max ) then
             if(printable) write(nfout,*) 'FCP CG2: Go next_cg_step'
             call next_cg_step()
          else
             if(printable) write(nfout,*) 'FCP CG2: Reach to max_CG_step. Go first_cg_step'
             call first_cg_step()
          endif
       else
          call linmin_set_new_cps()
       endif
    else
       f_para3 = f_total0*sign(1.0d0,vec_h)
       if(printable) write(nfout,'(a18,3f10.6)') ' CG2: F_parallel =', f_para1, f_para2, f_para3
       f_para3 = abs(f_para3)
       f_para2 = abs(f_para2)
       if ( ( f_para3 .lt. f_para2 ) .and. ( f_para3 .lt. f_para1 ) ) then
          if ( f_para3 .lt. f_para1*0.2d0 ) then
             if(printable) write(nfout,*) 'FCP CG2: F_parallel is small. Then linmin finished.'
             if ( iter_CG .lt. iter_CG_max ) then
                if(printable) write(nfout,*) 'FCP CG2: Go next_cg_step'
                call next_cg_step()
             else
                if(printable) write(nfout,*) 'FCP CG2: Reach to max_CG_step. Go first_cg_step'
                call first_cg_step()
             endif
          else
             if(printable)then
                write(nfout,*) 'FCP CG2: CG procedure is failed !!'
                write(nfout,*) 'FCP CG2: f_para3 must be less than f_para1*0.2d0.'
                write(nfout,*) 'FCP CG2: Then, new CG procedure starts.'
             endif
             call first_cg_step()
          endif
       else
          if(printable)then
             write(nfout,*) 'FCP CG2: Linmin calc is failed !!'
             write(nfout,*) 'FCP CG2: The 3rd point is not a minimum.'
             write(nfout,*) 'FCP CG2: Then, new CG procedure starts from the 2nd point.'
          endif
          call first_cg_step_from_cps2()
       endif
    endif

  contains

    subroutine first_cg_step()
      implicit none
      f_total1 = f_total0
      vec_g = f_total1
      vec_h = vec_g
      cps1 = cps
      vec_h_norm = abs( vec_h )
      call set_alpha()
      cps = cps1 + vec_h * alpha
      f_para1 = f_total1*sign(1.0d0,vec_h)
      if ( f_para1 .gt. 0.08d0 ) then
         iter_CG_max = 1
      else
         iter_CG_max = 5
      endif
      iter_CG = 1
      iter_linmin = 2
      write(nfout,'(a28,i2,a18,i2)') ' FCP CG2: Next step is iter_CG =', iter_CG, ' and iter_linmin =', iter_linmin
    end subroutine first_cg_step

    subroutine set_alpha()
      implicit none
      if      ( vec_h_norm .lt. 0.004d0 ) then
         alpha = 3.5d0
      else if ( vec_h_norm .lt. 0.007d0 ) then
         alpha = 3.0d0
      else if ( vec_h_norm .lt. 0.010d0 ) then
         alpha = 2.5d0
      else if ( vec_h_norm .lt. 0.015d0 ) then
         alpha = 2.0d0
      else if ( vec_h_norm .lt. 0.020d0 ) then
         alpha = 1.5d0
      else
         alpha = 1.0d0
      endif
    end subroutine set_alpha

    subroutine next_cg_step()
      f_total1 = f_total0
      !     gamma = sum(f_total1(:,:)**2) / sum(vec_g(:,:)**2)
      gamma = (f_total1-vec_g)*f_total1 / vec_g**2
      vec_g = f_total1
      vec_h = vec_g + gamma * vec_h
      cps1 = cps
      vec_h_norm = abs(vec_h)
      call set_alpha()
      cps = cps1 + vec_h * alpha
      f_para1 = f_total1*sign(1.0d0,vec_h)
      iter_CG = iter_CG + 1
      iter_linmin = 2
      write(nfout,'(a28,i2,a18,i2)') ' FCP CG2: Next step is iter_CG =', iter_CG, ' and iter_linmin =', iter_linmin
    end subroutine next_cg_step

    subroutine linmin_set_new_cps()
      implicit none
      cps2 = cps
      f_total2 = f_total0
      alpha2 = f_para1 * alpha / ( f_para1 - f_para2 )
      write(nfout,*) 'FCP CG2: alpha(force ) =', alpha2
      if ( alpha2 .lt. 0.0d0 ) then
         write(nfout,*) 'FCP CG2: Linmin calc is strange !!'
         write(nfout,*) 'FCP CG2: Alpha must be positive !!'
         write(nfout,'(a13,f15.10)') ' FCP CG2: Alpha =', alpha2
         write(nfout,*) 'FCP CG2: New CG step start from the 2nd point.'
         call first_cg_step()
      else
         if ( alpha2 .gt. 10.0d0*alpha ) then
            write(nfout,*) 'FCP CG2: Linmin calculated.'
            write(nfout,'(a22,f15.10)') ' FCP CG2: Original alpha =', alpha2
            write(nfout,*) 'FCP CG2: Alpha is too large, then, adjusted.'
            alpha2 = 10.0d0*alpha
         endif
         write(nfout,'(a13,f15.10)') ' FCP CG2: Alpha =', alpha2
         cps = cps1 + vec_h * alpha2
         iter_linmin = 3
         write(nfout,'(a28,i2,a18,i2)') ' FCP CG2: Next step is iter_CG =', iter_CG, ' and iter_linmin =', iter_linmin
      endif
    end subroutine linmin_set_new_cps

    subroutine first_cg_step_from_cps2()
      implicit none
      f_total1 = f_total2
      vec_g = f_total1
      vec_h = vec_g
      cps1 = cps2
      vec_h_norm = abs(vec_h)
      call set_alpha()
      cps = cps1 + vec_h * alpha
      f_para1 = f_total1*vec_h/abs(vec_h)
      iter_CG = 1
      iter_linmin = 2
      write(nfout,'(a28,i2,a18,i2)') ' FCP CG2: Next step is iter_CG =', iter_CG, ' and iter_linmin =', iter_linmin
    end subroutine first_cg_step_from_cps2

    !! jn_130104

  end subroutine m_Fcp_cg2

  logical function m_Fcp_Converged()
    implicit none
    m_Fcp_Converged = ( abs( fcp_mu - efermi ) < fcp_relax_crit )
    return
  end function m_Fcp_Converged

! QUENCH

  subroutine m_Fcp_md(mdalg)
    integer, intent(in)     :: mdalg
    real(kind=DP) :: forc_l

    forc_l = fcp_mu - efermi

    if (mdalg .ne. VERLET) then
       cpd_old = cpd_l
       if(iprifcp >= 1) then
          write(nfout,'(" -- FCP cpd_l, cpd_old (before evolve_velocities)--")')
          write(nfout,'(f18.5)') cpd_l
       end if
       if(imdalg == STEEPEST_DESCENT .and. mode_fi_coefficient == OFF) then
          fi_coefficient = get_fi_coefficient(abs(forc_l))
          if(iprifcp >= 1) write(nfout,'(a19,f20.10)') ' fi_coefficient =  ', fi_coefficient
       end if
       call evolve_velocities(mdalg,forc_l) !-(m_Ionic_System)  cpd_l = cpd_l + dtio/mass * forc_l
    end if

    ! <<--
    if(mdalg == QUENCHED_MD) call quench_velocities(forc_l) ! -->cpd_l

    cps = cps + dtio*cpd_l

    return
  end subroutine m_Fcp_md

  subroutine evolve_velocities(mdalg,forc_l)
    implicit none
    integer, intent(in) :: mdalg
    real(kind=DP), intent(inout) :: forc_l

    real(kind=DP) :: fac, frc, rm, pdot
    real(kind=DP) :: f_dot_v, vec_norm2, force_t_norm
    real(kind=DP) :: force_t_para, force_t_vert

    if(mdalg == STEEPEST_DESCENT) then
       cpd_l=forc_l/dtio*fi_coefficient                   !! moving force cal
    else
       fac = dtio/fcp_mass
       frc = forc_l
       if (mdalg == VERLET .and. iteration_ionic == 1)  then
          !!cpd_l(ia,1:3) = fac*frc(1:3) / 2.d0
          !!cpd_old(ia,1:3) = -cpd_l(ia,1:3)
          cpd_l = cpd_l + 0.5d0*fac*frc
       else
          cpd_l = cpd_l + fac*frc
       endif
    end if

  end subroutine evolve_velocities

  real(kind=DP) function get_fi_coefficient(force_t_norm)
    implicit none
    real(kind=DP), intent(in) :: force_t_norm
    if(force_t_norm < 0.008d0) then
       get_fi_coefficient = 3.5d0
    else if (force_t_norm < 0.02d0) then
       get_fi_coefficient = 2.5d0
    else if (force_t_norm < 0.05d0) then
       get_fi_coefficient = 1.5d0
    else
       get_fi_coefficient = 0.8d0
    end if
  end function get_fi_coefficient

  subroutine quench_velocities(forc_l)
    implicit none
    real(kind=DP), intent(in) :: forc_l
    real(kind=DP) :: cpd_parallel, cpd_vertical, f, fv, fcpd, v
    v = (cpd_old + cpd_l)*0.5
    fv = forc_l*v
    if(fv < 0.d0) then
       if(printable) write(nfout,'(" FCP ---<<quench_velocity>>")')
       cpd_l = 0.d0
    end if

  end subroutine quench_velocities

! GDIIS, BFGS

  subroutine m_Fcp_gdiis_alloc(mode_init,if_allocated)
    implicit none
    integer, intent(in) ::  mode_init
    integer, intent(out) :: if_allocated
    integer ::             i,j
    real(kind=DP) ::       fac
    if(kqnmditer_p <= 0) call phase_error_with_msg(nfout,' FCP kqnmditer_p is illegal (m_Fcp_gdiis_alloc)',__LINE__,__FILE__)
    if(printable) write(nfout,'(" !! FCP kqnmditer_p = ",i6," <<m_Fcp_gdiis_alloc>>")') kqnmditer_p
    if(.not.allocated(u_l)) allocate(u_l(kqnmditer_p));u_l=0.d0
    if(.not.allocated(w_l)) allocate(w_l(kqnmditer_p));w_l=0.d0
    if(mode_init == UNIT) then
       fac = 1.d0
    else
       fac = dtio*dtio/fcp_mass
    end if
    fc_l = -fac

    if(.not.allocated(ncrspd)) allocate(ncrspd(kqnmditer_p))
    ncrspd(:) = (/(i,i=1,kqnmditer_p)/)
    if(.not.allocated(f_gdiis)) allocate(f_gdiis(kqnmditer_p,kqnmditer_p))
    if(.not.allocated(g)) allocate(g(kqnmditer_p))
    if(.not.allocated(f_wk)) allocate(f_wk(kqnmditer_p*kqnmditer_p,2))
    if(.not.allocated(f_rslv)) allocate(f_rslv(kqnmditer_p*kqnmditer_p))
    if(.not.allocated(e_wk)) allocate(e_wk(kqnmditer_p*kqnmditer_p))
    if(.not.allocated(ww1)) allocate(ww1(kqnmditer_p))
    if(.not.allocated(ip)) allocate(ip(kqnmditer_p))
    if_allocated = 1
  end subroutine m_Fcp_gdiis_alloc

  subroutine m_Fcp_gdiis_reset()
    implicit none
    iter_gdiis = 0
    if(if_allocated==0) return
    call m_Fcp_gdiis_dealloc(if_allocated)
  end subroutine m_Fcp_gdiis_reset

  subroutine m_Fcp_gdiis_dealloc(if_allocated)
    implicit none
    integer, intent(out) :: if_allocated
    deallocate(ip);
    deallocate(ww1);   deallocate(e_wk);    deallocate(f_rslv)
    deallocate(f_wk);   deallocate(g);     deallocate(f_gdiis)
    deallocate(ncrspd);   deallocate(w_l)
    deallocate(u_l)
    if_allocated = 0
  end subroutine m_Fcp_gdiis_dealloc

  subroutine m_Fcp_gdiis()
    implicit none
    real(kind=DP) :: forcmx,etotal
    real(kind=DP) :: forc_l

    integer       :: m, mode_init = 0,it,nsum,icount,icon, imd_t
    real(kind=DP) :: c_dist, forcmx_mdfy

    forc_l = fcp_mu - efermi
    etotal = forc_l ** 2

    forcmx = fcp_relax_crit
    forcmx_mdfy = forcmx

    iter_gdiis = iter_gdiis + 1
    if(if_allocated == 0) call m_Fcp_gdiis_alloc(mode_init,if_allocated)

    forc_g = forc_l
    m = mod(iter_gdiis,kqnmditer_p)
    if(m == 0) m = kqnmditer_p

    call d_fc_l(mode_init,forcmx_mdfy) ! -> fc_l
    call rot_ncrspd(nsum)  ! -> ncrspd, nsum
    it = ncrspd(nsum)
    if(it == 0) it = 1
    call stor_cps_forc(it) ! cps,forc_g -> u_l, w_l
    if (imdalg==GDIIS)then
       if(nsum > 1) then
          call gdiis_mat(nsum) ! ncrspd,w_l -> f_gdiis
          call getmatinv(nsum,icon) ! -> f_gdiis := f_gdiis^(-1)
          if(icon == 0) then
             call gdiis_g(nsum)   ! -> g
             call opt_forc(nsum)  ! -> forc_g := sum_{it}(g(it)*w_l(,ncrspd(it))
!!$                call mdfy_forc(forcmx,forcmx_mdfy2)
             if(iprifcp >= 2 .and. printable) call forc_check(" --- FCP new force ---")
             if(iprifcp >= 1) call cmp_fop_fin(it)
             call opt_geom(nsum)  ! -> cps := sum_{it}(g(it)*u_l(,ncrspd(it))
          else   ! when icon /= 0, f_gdiis matrix is singular
             if(printable) write(nfout,'(" !! FCP f_gdiis matrix is singular <<m_Fcp_gdiis>>")')
             call m_Fcp_md(QUENCHED_MD)
          end if
       end if

       cpd_l = -fc_l/dtio*forc_g
       cps   = cps + dtio*cpd_l

    else if(imdalg==BFGS) then
       call do_bfgs(nsum)
    end if
    if(imdalg==GDIIS)then
       if(iprifcp >= 2 .and. printable) call cps_check(" --- FCP new coordinates ---") ! cps

       c_dist = 0.3d0
       call cps_damp(c_dist,it,icount)
       if(icount >= 1 .and. iprifcp >= 1 .and. printable) then
          call cps_check(" --- FCP new (damped) coordinates ---") ! cps
       endif
    endif

    if(iprifcp >= 2 .and. printable) call cpd_check(" --- FCP velocity ---") ! cpd_l

  contains

    subroutine do_bfgs(nsum)
      implicit none
      integer, intent(in) :: nsum
      integer :: i,j,k,i1,j1,itr0,itr1,info
      real(DP) :: xgi,gihg
      real(kind=DP) :: gdelta,xdelta
      logical :: corrected_eig
      real(DP) :: gdotinvh
      real(DP) :: tmpforc
      real(DP) :: ihess
      real(DP) :: amat
      real(DP) :: eigv
      real(DP) :: eigvec
      real(DP) :: workar
      real(DP) :: maxoptforc,tmpmaxoptforc
      real(DP),save :: maxoptforc_old=1000.d0
      real(DP),save :: energy_old = 0.0d0
      real(DP),save :: force_old  = 0.0d0
      real(DP),save :: alpha_bfgs = 1.0d0
      logical, save :: lfirst = .true.
      integer :: icounti,icountj
      real(DP) :: c1,c2,mina,maxa,incre
      real(DP) :: tmpsum1,tmpsum2,e0,f0
      c1 = wolfe_c1
      c2 = wolfe_c2
      mina = min_alpha
      maxa = max_alpha
      incre = 1.1d0
      xdelta=0.d0
      gdelta=0.d0
      gdotinvh=0.d0
      tmpforc=0.d0
      ihess=0.d0
      if(sw_correct_eigenvalue==ON)then
         amat=0.d0
         eigv=0.0d0
         eigvec=0.0d0
         workar=0.0d0
      endif

      !      build inverse of the Hessian
      ihess = 1.0d0
      do i=2,nsum
         itr1 = ncrspd(i)
         itr0 = ncrspd(i-1)
         icountj=0

         icountj=icountj+1
         xdelta =  u_l(itr1)-u_l(itr0)
         gdelta = -w_l(itr1)+w_l(itr0)
         xgi = 1.0d0/(xdelta*gdelta)
         if(xgi<0)then
            if(printable .and. iprifcp>=2)then
               write(nfout,'(a,i3)') '!** WARNING FCP dx dot dg is negative for history ',itr1
               write(nfout,'(a)') 'skipping this update.'
            endif
            cycle
         endif
         gdotinvh = ihess*gdelta
         gihg = gdelta*gdotinvh
         ihess = ihess+xgi*xgi*(1.0d0/xgi+gihg)*xdelta*xdelta &
            - (gdotinvh*xdelta+gdotinvh*xdelta)*xgi
      enddo

      !      correct bad eigenvalues present in the Hessian
      if(sw_correct_eigenvalue==ON)then
         corrected_eig=.false.
         amat = ihess
         call dspev('V','U',1,amat,eigv,eigvec,1,workar,info)
         if(printable .and. iprifcp>=2) write(nfout,'(a)') '--- FCP eigenvalues for the approximate Hessian ---'
         if(printable.and.iprifcp>=2) write(nfout,'(f20.10)') 1.0d0/eigv
         if (1.0d0/eigv<eigenvalue_threshold)then
            eigv = 1.0d0/eigenvalue_threshold
            if(printable.and.iprifcp>=2) write(nfout,'(a,f20.10)') &
               &  'corrected the eigenvalue : ',1.0d0/eigv
            corrected_eig=.true.
         endif
         if(corrected_eig)then
            ihess = ihess+eigvec*eigvec*eigv
         endif
      endif

      !      H^-1 dot g
      tmpforc = - ihess*forc_l
      tmpmaxoptforc=abs(tmpforc)
      if(tmpmaxoptforc>maxoptforc)maxoptforc=tmpmaxoptforc

      if(printable .and. iprifcp>=1) then
         write(nfout,'(a,f20.10)') 'FCP max. optimal force obtained from the BFGS update : ',maxoptforc
      endif

      if(sw_optimize_alpha == ON) then
         if( lfirst ) then
            lfirst = .false.
         else
            tmpsum1 = force_old*tmpforc
            tmpsum2 = forc_l*tmpforc
            e0 = energy_old + c1*alpha_bfgs*tmpsum1
            f0 = -c2*tmpsum1
            if(e0>etotal .and. f0>abs(tmpsum2))then
               alpha_bfgs = alpha_bfgs*incre
            else
               alpha_bfgs = alpha_bfgs/incre
            endif
            if(alpha_bfgs>maxa) alpha_bfgs = maxa
            if(alpha_bfgs<mina) alpha_bfgs = mina
            if(printable)then
               write(nfout,'(a,2f20.10)') 'FCP Wolfe condition  i) : ',e0,etotal
               write(nfout,'(a,2f20.10)') 'FCP Wolfe condition ii) : ',f0,abs(tmpsum2)
               write(nfout,'(a,f20.10)')  'FCP new alpha           : ',alpha_bfgs
            endif
         endif

         energy_old = etotal
         force_old  = forc_l
      endif

      if(maxoptforc_old*10<maxoptforc)then
         if (printable) write(nfout,'(a)') 'FCP the estimated force seems to be very large; &
            & update will be done by the steepest-descent method'
         cps = cps+forc_l
      else
         tmpsum1=0.d0
         tmpsum2=0.d0
         cps = cps-alpha_bfgs*tmpforc
      endif
      maxoptforc_old = maxoptforc

    end subroutine do_bfgs

    subroutine rot_ncrspd(nsum)
      implicit none
      integer, intent(out) :: nsum
      integer ::              nbox, istrbr, itemp
      nbox = (iter_gdiis-1)/kqnmditer_p
      if(gdiis_hownew == ANEW) then
         istrbr = nbox*kqnmditer_p + 1
      else if (gdiis_hownew == RENEW) then
         if(nbox == 0) then
            istrbr = 1
         else
            istrbr = iter_gdiis - (kqnmditer_p-1)
            itemp = ncrspd(1)
            do it = 1, kqnmditer_p-1
               ncrspd(it) = ncrspd(it+1)
            end do
            ncrspd(kqnmditer_p) = itemp
         end if
      end if
      nsum = iter_gdiis - istrbr + 1
      if(iprifcp >= 2 .and. printable) then
         write(nfout,'(" -- FCP rot_ncrspd -- ")')
         write(nfout,'("   i  : ",8i8)') (it,it=1,nsum)
         write(nfout,'("ncrspd: ",8i8)') (ncrspd(it),it=1,nsum)
      end if
    end subroutine rot_ncrspd

    subroutine d_fc_l(mode_init,forcmx)
      implicit none
      integer      , intent(in) :: mode_init
      real(kind=DP), intent(in) :: forcmx
      if(forcmx >= c_forc_prop_region_high) then
         if(mode_init == UNIT) then
            fc_l = -1.d0
         else
            fc_l = -dtio*dtio/fcp_mass
         end if
      else if(forcmx >= c_forc_prop_region_low) then
         fc_l = - factor_prop_region/forcmx
         if(dabs(fc_l) < 1.d0) then
            fc_l = -1.d0
         end if
      end if
    end subroutine d_fc_l

    subroutine stor_cps_forc(it)
      implicit none
      integer, intent(in) :: it
      u_l(it) = cps
      w_l(it) = forc_g
    end subroutine stor_cps_forc

    subroutine gdiis_mat(nsum)
      implicit none
      integer, intent(in) :: nsum

      integer ::             it, jt, itcrspd, jtcrspd
      real(kind=DP) ::       xmul

      do it = 1, nsum
         itcrspd = ncrspd(it)
         do jt = it, nsum
            jtcrspd = ncrspd(jt)
            xmul = w_l(jtcrspd)*w_l(itcrspd)
            f_gdiis(it,jt) = xmul
            if(jt /= it) f_gdiis(jt,it) = f_gdiis(it,jt)
         end do
         if(iprifcp >= 2 .and. printable) then
            if(it == 1) write(nfout,'(" -- FCP f_gdiis -- ")')
            write(nfout,'( 6d12.4)') (f_gdiis(it,jt),jt=1,nsum)
         end if
      end do
    end subroutine gdiis_mat

    subroutine getmatinv(nsum,icon)
      implicit none
      integer ,intent(in) ::  nsum
      integer, intent(out) :: icon

      real(kind=DP) :: div
      integer ::       it, jt, ipfr, ipto
#ifdef _GDIIS_MAT_CHECK_
      integer ::       kt, ik_count, kj_count
#endif

      div = 1.0/f_gdiis(1,1)
      do it = 1, nsum
         do jt = 1, nsum
            ipto = it + (jt-1)*nsum
            f_wk(ipto, 1) = f_gdiis(it,jt)*div
#ifdef _GDIIS_MAT_CHECK_
            f_wk(ipto, 2) = f_gdiis(it,jt)*div
#endif
            if(it == jt) then
               e_wk(ipto) = 1.d0
            else
               e_wk(ipto) = 0.d0
            end if
         end do
      end do

      call rdecomp(nsum,f_wk(1,1),ww1,ip,icon)

      if(icon /= 0) then
         if(printable) then
            write(nfout,'("  FCP [f_wk] after rdecomp <<m_Fcp_gdiis_getmatinv>>")')
            do it = 1, nsum
               write(nfout,'(" f_wk(:,",i3," ) = ",8d12.4)') it,(f_wk(it+(jt-1)*nsum,1),jt=1,nsum)
            end do
            write(nfout,'(" ip = ",8i6)') (ip(ipto),ipto=1,nsum)
            write(nfout,*) ' FCP LU decomposition is impossible. <<m_Fcp_gdiis.getmatinv>>'
         end if
         return
!!$         stop ' LU decompoition is impossible <<m_Fcp_gdiis_getmatinv>>'
      else
         call rsolve(nsum,nsum,f_wk(1,1),e_wk,f_rslv,ip)
      endif
#ifdef _GDIIS_MAT_CHECK_
      ! ----------- checking of inversion matrix ------------>
      if(printable) then
         write(nfout,*) ' -- FCP below should equal to a unit matrix --'
         do it = 1, nsum
            ww1 = 0.d0
            do jt = 1, nsum
               do kt = 1, nsum
                  ik_count = (kt-1)*nsum + it
                  kj_count = (jt-1)*nsum + kt
                  ww1(jt) = ww1(jt) + f_wk(ik_count,2)*f_rslv(kj_count)
               enddo
            enddo
            write(nfout,9009) it,(ww1(jt),jt=1,nsum)
         enddo
9009     format(i3,8f8.4,/,3x,8f8.4,/,3x,8f8.4,/3x,8f8.4)
      end if
      ! <--------------------------------------------------
#endif

      do jt = 1, nsum
         do it = 1, nsum
            ipfr = it + (jt-1)*nsum
            f_gdiis(it,jt) = f_rslv(ipfr)
         end do
      end do

      if(iprifcp >= 2 .and. printable) then
         write(nfout,*) ' **inverse matrix**'
         do jt = 1, nsum
            write(nfout,9008) jt,(f_gdiis(it,jt),it=1,nsum)
         enddo
      end if
9008  format(i3,(6d12.4))

    end subroutine getmatinv

    subroutine gdiis_g(nsum)
      implicit none
      integer, intent(in) :: nsum
      real(kind=DP) ::       alpha_fac
      integer ::             it,jt
      alpha_fac = 0.d0
      do it = 1, nsum
         do jt = 1, nsum
            alpha_fac = alpha_fac + f_gdiis(it,jt)
         end do
      end do
      alpha_fac = 1.d0/alpha_fac

      g = 0.d0
      do it = 1, nsum
         do jt = 1, nsum
            g(it) = g(it) + alpha_fac*f_gdiis(it,jt)
         end do
      end do

      if(iprifcp >= 2 .and. printable) then
         write(nfout,*) ' ---- FCP alpha_fac     ----'
         write(nfout,9008) alpha_fac
         write(nfout,*) ' ---- a(1:',nsum,') ----'
         write(nfout,9008) (g(it),it=1,nsum)
      end if
9008  format(8f20.12)
    end subroutine gdiis_g

    subroutine opt_forc(nsum)
      implicit none
      integer, intent(in) :: nsum
      integer ::             it, itcrspd

      forc_g = 0.d0
      do it = 1, nsum
         itcrspd = ncrspd(it)
         forc_g = forc_g + g(it)*w_l(itcrspd)
      end do
    end subroutine opt_forc

    subroutine opt_geom(nsum)
      implicit none
      integer, intent(in) :: nsum
      integer ::             ifc, itcrspd

      cps = 0.d0
      do ifc = 1, nsum
         itcrspd = ncrspd(ifc)
         cps = cps + g(ifc)*u_l(itcrspd)
      end do

    end subroutine opt_geom

    subroutine cmp_fop_fin(it)
      implicit none
      integer, intent(in) :: it
      real(kind=DP) :: xmul, xmul_input_f
      xmul = forc_g**2
      xmul_input_f = w_l(it)**2
      if(printable) write(nfout,*) ' FCP norm of optimal force = ', xmul
      if(printable) write(nfout,*) ' FCP norm of input   force = ', xmul_input_f
    end subroutine cmp_fop_fin

    subroutine cps_damp(c_dist,it,icount)
      implicit none
      real(kind=DP), intent(in) :: c_dist
      integer, intent(in) ::      it
      integer ,intent(out) ::     icount
      real(kind=DP) :: df

      icount = 0
      df = cps - u_l(it)
      if(df > c_dist) then
         icount = icount + 1
         cps = u_l(it) + c_dist
      else if(df < -c_dist) then
         icount = icount + 1
         cps = u_l(it) - c_dist
      end if

    end subroutine cps_damp

    subroutine forc_check(name)
      character(len=*),intent(in) :: name
      write(nfout,'("FCP ",a40)') name
      write(nfout,'(f20.12)') forc_g
    end subroutine forc_check

    subroutine cps_check(name)
      character(len=*),intent(in) :: name
      write(nfout,'("FCP ",a40)') name
      write(nfout,'(f20.12)') cps
    end subroutine cps_check

    subroutine cpd_check(name)
      character(len=*),intent(in) :: name
      write(nfout,'("FCP ",a40)') name
      write(nfout,'(f20.12)') cpd_l
    end subroutine cpd_check

  end subroutine m_Fcp_gdiis

! NEB

  subroutine m_Fcp_create_replica(nimage)
    use m_Control_Parameters, only : fcp_tot_charge_first, fcp_tot_charge_last
    implicit none
    integer, intent(in) :: nimage
    logical, save :: initialized = .false.

    num_of_images = nimage
    ALLOCATE( neb_image( nimage ) )

  end subroutine m_Fcp_create_replica

  subroutine m_Fcp_initialize_neb_totch()
    use m_Control_Parameters, only : fcp_tot_charge_first, fcp_tot_charge_last
    implicit none
    integer       :: n, i
    real(kind=DP) :: first, last

    if( cps_initialized ) then
       return
    end if
    cps_initialized = .true.

    first = fcp_tot_charge_first
    last  = fcp_tot_charge_last
    n = num_of_images
    DO i = 1, n
       neb_image(i)%cps = cps + ( first * (n - i) + last * (i - 1) ) / (n - 1)
    END DO
    neb_image(:)%cps0 = neb_image(:)%cps

  end subroutine m_Fcp_initialize_neb_totch

  subroutine m_Fcp_set_neb_replica_totch(i)
    implicit none
    integer, intent(in) :: i
    cps = neb_image(i)%cps
  end subroutine m_Fcp_set_neb_replica_totch

  subroutine m_Fcp_set_neb_replica_force(i)
    implicit none
    integer, intent(in) :: i
    neb_image(i)%force = fcp_mu - efermi
  end subroutine m_Fcp_set_neb_replica_force

  subroutine m_Fcp_update_neb(dt, time_integral)
    use m_Parallelization , only : mype
    use m_Files           , only : nfneb
    implicit none
    real(8)     , intent(in) :: dt
    character(*), intent(in) :: time_integral

    integer :: i
    real(8) :: mass
    real(8) :: r, r0, r1, v, v1, f, f0, unitf

    if (mype==0) write(nfneb,*) 'FCP coordinate update...'

    do i = 2, num_of_images-1

       neb_image(i)%cps0 = neb_image(i)%cps

       r  = neb_image(i)%cps
       r0 = neb_image(i)%cps0
       v  = neb_image(i)%velocity
       f  = neb_image(i)%force
       f0 = neb_image(i)%force0
       mass = fcp_mass

       select case(time_integral)
       case('steepest_descent')
          r1 = r + dt*f
          v = 0.0d0
       case('verlet')
          r1 = 2.0d0*r - r0 + dt*dt*f/mass
          v1 = (r1-r)/(2.0d0*dt)
          v = 0.0d0
          if( v1*f >= 0.0d0) then
             if( abs(f) /= 0.0d0 ) then
                v = f * v1/f
             end if
          end if
       case('velocity_verlet')
          v1 = v + 0.5d0*dt*(f0+f)/mass
          v = 0.0d0
          unitf=0.d0
          if( abs(f) /= 0 ) unitf=f/abs(f)
          if( v1*unitf >= 0.0d0) then
             v = v1
          else
             if (mype==0) write(nfneb,'(a,i6)') 'FCP quenched velocity'
          end if
          r1 = r + v*dt + f*dt**2/mass
       end select

       neb_image(i)%cps      = r1
       neb_image(i)%velocity = v

       if ( mype==0 ) then
          write(nfneb,*) ' image: ', i
          write(nfneb,*) '  cps, cps0, cps-cps0 '
          write(nfneb,'(2f10.5,e11.3)') neb_image(i)%cps, neb_image(i)%cps0, &
             neb_image(i)%cps - neb_image(i)%cps0
       endif

    end do

  end subroutine m_Fcp_update_neb

  logical function m_Fcp_replica_Converged()
    implicit none
    m_Fcp_replica_Converged = all( abs( neb_image(2:num_of_images)%force ) < fcp_relax_crit )
    return
  end function m_Fcp_replica_Converged

  subroutine m_Fcp_set_neb_parameter(itr)
    implicit none
    integer, intent(in) :: itr
    if( itr == 1 ) return
    neb_image(:)%force0 = neb_image(:)%force
    neb_image(:)%force  = 0.0d0
  end subroutine m_Fcp_set_neb_parameter

  subroutine m_Fcp_allreduce_neb_force(nrank_r)
    use m_Parallelization,      only : MPI_CommGroup
    use mpi
    implicit none
    integer, intent(in) :: nrank_r
    integer i
    integer mype, npes, mpi_err
    real(8),allocatable :: work1(:)

!    include 'mpif.h'

    if (nrank_r.lt.1) return

    call mpi_comm_size(MPI_CommGroup,npes,mpi_err)

    allocate(work1(num_of_images));work1=0.d0
    call mpi_allreduce(neb_image(:)%force,work1,num_of_images, &
       mpi_double_precision,mpi_sum,mpi_comm_world,mpi_err)
    neb_image(:)%force = work1(:)
    deallocate(work1)

    ! temporary    npes = ne*nk
    neb_image(:)%force = neb_image(:)%force/npes

  end subroutine m_Fcp_allreduce_neb_force

  subroutine m_Fcp_write_down_totch_replica(i)
    use m_Files, only : nfnebcntn
    implicit none
    integer, intent(in) :: i

    write(nfnebcntn,*) '(totch)'
    write(nfnebcntn,'(d24.16)') neb_image(i)%cps

  end subroutine m_Fcp_write_down_totch_replica

  subroutine m_Fcp_read_totch_replica(i)
    use m_Files, only : nfnebcntn
    implicit none
    integer, intent(in) :: i

    read(nfnebcntn,*)
    read(nfnebcntn,*) neb_image(i)%cps

    neb_image(i)%cps0 = neb_image(i)%cps

    cps_initialized = .true.

  end subroutine m_Fcp_read_totch_replica

end module m_Fcp
