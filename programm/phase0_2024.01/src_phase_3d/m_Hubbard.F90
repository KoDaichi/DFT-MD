!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 606 $)
!
!  MODULE: m_Hubbard
!
!  AUTHOR(S): T. Yamamoto   October/11/2004
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
!=======================================================================
!
!
module m_Hubbard
! $Id: m_Hubbard.f90 606 2020-04-15 06:45:49Z ktagami $
  use m_Const_Parameters,     only : DP,ON,ATOMIC_ORBITAL,SPHERICAL_HARMONICS
  use m_Files,                only : nfout
  use m_Control_Parameters,   only : nspin,ipri,printable,num_projectors &
       &                           , proj_attribute, af, projector_type &
       &                           , sw_constraint, const_site, const_alpha &
       &                           , iprihubbard
  use m_Ionic_System,         only : natm,ntyp,ityp,ihubbard
  use m_PseudoPotential,      only : prodphi,ilmt,ltp,mtp,taup,nlmt,ntau,nlmtt
  use m_Crystal_Structure,    only : op,nopr
  use m_Ionic_System,         only : napt
  use m_Orbital_Population,   only : i2lp,nyy,ilmt_yy,ommix
  use m_Electronic_Structure, only : dhub

! =============================== added by K. Tagami =============== 5.0
  use m_Const_Parameters, only :  OccMat_Type1, OccMat_Type2, FLL, AMF
  use m_Orbital_Population,   only : om
  use m_Control_Parameters,   only : occmat_type, sw_eval_energy_before_charge, &
       &                             dftu_type
! =================================================================== 5.0

! =============================== added by K. Tagami =============== 11.0
  use m_Const_Parameters, only :  CMPLDP
  use m_Control_Parameters, only : ndim_chgpot, ndim_magmom, ndim_spinor
  use m_ES_NonCollinear,  only : m_ES_MagMom_To_DensMat_porb, &
       &                         m_ES_DensMat_To_MagMom_dhub
  use m_Orbital_Population,  only : max2lp
! ================================================================== 11.0

  use m_Orbital_Population,   only : om_aimag, ommix_aimag
  use m_Electronic_Structure, only : dhub_aimag

  implicit none

! ======================== added by K. Tagami======================== 5.0
  real(kind=DP) :: Ueff_prefactor = 1.0D0
! =================================================================== 5.0

contains

  subroutine m_Hubbard_energy(energy)
    real(kind=DP), intent(out) :: energy

    integer :: is,ia,it,i1,i2,ih,ie
    real(kind=DP) :: omsum, c1, nsigma
    real(kind=DP), pointer ::  om_wk(:,:,:,:,:)
    real(kind=DP), allocatable :: delta_om(:,:)

    energy=0.d0
    if(sw_constraint == ON) return

    if ( sw_eval_energy_before_charge==ON ) then
       om_wk => ommix
    else
       om_wk => om
    endif

    do is=1,nspin,af+1
       do ia=1,natm
          ih = ihubbard(ia)
! === Debug by Intel "-check all" option! by T.Kato 2011/03/28 =================
!         ie = proj_attribute(ih)%ielem
!         if(ih /= 0) then
          if(ih /= 0) then
             ie = proj_attribute(ih)%ielem
! ==============================================================================
             it=ityp(ia)
             omsum=0.d0
!
             select case( dftu_type )
             case (FLL)
                omsum=0.d0
                do i2=1,i2lp(ih)
                   omsum = omsum +om_wk(i2,i2,ie,ia,is)
                   do i1=1,i2lp(ih)
                      omsum = omsum -om_wk(i2,i1,ie,ia,is)*om_wk(i1,i2,ie,ia,is)
                   end do
                end do
                energy = energy +0.5d0*proj_attribute(ih)%Ueff*omsum

             case (AMF)
                allocate( delta_om(max2lp,max2lp) ); delta_om = 0.0d0
                delta_om(1:i2lp(ih),1:i2lp(ih)) &
                     &  = om_wk(1:i2lp(ih),1:i2lp(ih),ie,ia,is)

                c1 = 0.0d0
                do i1=1,i2lp(ih)
                   c1 = c1 +om_wk(i1,i1,ie,ia,is)
                End do
                nsigma = c1 /( 2*proj_attribute(ih)%l +1 )
                do i1=1,i2lp(ih)
                   delta_om(i1,i1) = delta_om(i1,i1) -nsigma
                end do

                omsum=0.d0
                do i2=1,i2lp(ih)
                   do i1=1,i2lp(ih)
                      omsum = omsum +delta_om(i1,i2)*delta_om(i2,i1)
                   end do
                end do
                energy = energy -0.5d0*proj_attribute(ih)%Ueff*omsum
                deallocate( delta_om )
             end select

          end if
       end do
    end do
    energy=energy*(af+1)
    if(nspin==1) energy=energy*2.d0

! =============================== added by K. Tagami ==================== 5.0
    energy = energy *Ueff_prefactor
! ======================================================================= 5.0

  end subroutine m_Hubbard_energy

! ======================== added by K. Tagami ======================== 11.0
  subroutine m_Hubbard_energy_noncl(energy)
    real(kind=DP), intent(out) :: energy

    real(kind=DP), allocatable :: dmmat_r_magmom( :,:,: )
    real(kind=DP), allocatable :: dmmat_i_magmom( :,:,: )
    complex(kind=CMPLDP), allocatable :: dmmat_ssrep( :,:,: )

    integer :: ia, it,i1,i2,ih,ie
    integer :: is, size1
    real(kind=DP) :: nsigma
    complex(kind=CMPLDP) :: dsum

    real(kind=DP), pointer ::  om_wk_r(:,:,:,:,:), om_wk_i(:,:,:,:,:)
    complex(kind=CMPLDP), allocatable :: delta_dm(:,:,:)

    if(sw_constraint == ON) return
    energy=0.d0

    if ( sw_eval_energy_before_charge==ON ) then
       om_wk_r => ommix;       om_wk_i => ommix_aimag
    else
       om_wk_r => om;          om_wk_i => om_aimag
    endif

    do ia=1,natm
       ih = ihubbard(ia)
       if(ih /= 0) then
          ie = proj_attribute(ih)%ielem
          it=ityp(ia)

          size1 = i2lp(ih)
          allocate( dmmat_r_magmom(size1,size1,ndim_magmom) ); dmmat_r_magmom = 0.0d0
          allocate( dmmat_i_magmom(size1,size1,ndim_magmom) ); dmmat_i_magmom = 0.0d0
          allocate( dmmat_ssrep(size1,size1,ndim_chgpot) );    dmmat_ssrep = 0.0d0

          if ( sw_eval_energy_before_charge == ON ) then
             dmmat_r_magmom(:,:,:) = ommix(:,:,ie,ia,:)
             dmmat_i_magmom(:,:,:) = ommix_aimag(:,:,ie,ia,:)
          else
             dmmat_r_magmom(:,:,:) = om(:,:,ie,ia,:)
             dmmat_i_magmom(:,:,:) = om_aimag(:,:,ie,ia,:)
          endif
          call m_ES_MagMom_To_DensMat_porb( size1**2, dmmat_r_magmom, dmmat_i_magmom, &
               &                            dmmat_ssrep )

          select case( dftu_type )
          case (FLL)
             dsum=0.d0
             Do is=1, ndim_chgpot
                if ( is==1 .or. is==ndim_chgpot ) then
                   do i2=1,i2lp(ih)
                      dsum = dsum + dmmat_ssrep(i2,i2,is)
                   end do
                endif
             End do
             Do is=1, ndim_chgpot
                do i2=1,i2lp(ih)
                   do i1=1,i2lp(ih)
                      dsum = dsum - dmmat_ssrep(i1,i2,is) &
                           &        *conjg( dmmat_ssrep(i1,i2,is) )
                   end do
                end do
             End do
             energy = energy + 0.5d0*proj_attribute(ih)%Ueff *real(dsum)

          case (AMF)
             allocate( delta_dm(size1,size1,ndim_chgpot) );    delta_dm = 0.0d0
             delta_dm(1:size1,1:size1,1:ndim_chgpot) &
                  &  = dmmat_ssrep(1:size1,1:size1,1:ndim_chgpot)

             dsum = 0.0d0
             Do is=1, ndim_chgpot
                if ( is==1 .or. is==ndim_chgpot ) then
                   do i2=1,i2lp(ih)
                      dsum = dsum + dmmat_ssrep(i2,i2,is)
                   end do
                endif
             End do

             nsigma = dsum /( 2*proj_attribute(ih)%l +1 ) /ndim_spinor

             Do is=1, ndim_chgpot
                if ( is==1 .or. is==ndim_chgpot ) then
                   do i2=1,i2lp(ih)
                      delta_dm(i1,i1,is) = delta_dm(i1,i1,is) -nsigma
                   end do
                endif
             End do

             dsum=0.d0
             Do is=1, ndim_chgpot
                do i2=1,i2lp(ih)
                   do i1=1,i2lp(ih)
                      dsum = dsum + delta_dm(i1,i2,is) &
                           &        *conjg( delta_dm(i1,i2,is) )
                   end do
                end do
             End do
             energy = energy -0.5d0*proj_attribute(ih)%Ueff *real(dsum)
             deallocate( delta_dm )
          end select

          deallocate( dmmat_r_magmom, dmmat_i_magmom )
          deallocate( dmmat_ssrep )
       end if
    end do

    energy = energy *Ueff_prefactor

  end subroutine m_Hubbard_energy_noncl

  subroutine m_Hubbard_energy2_noncl(energy)           ! for safety
    real(kind=DP), intent(out) :: energy

    real(kind=DP), allocatable :: dmmat_r_magmom( :,:,: )
    real(kind=DP), allocatable :: dmmat_i_magmom( :,:,: )
    complex(kind=CMPLDP), allocatable :: dmmat_ssrep( :,:,: )

    integer :: ia, it,i1,i2,ih,ie
    integer :: is, size1
    complex(kind=CMPLDP) :: dsum

    if(sw_constraint == ON) return

    energy=0.d0

    allocate( dmmat_r_magmom(max2lp,max2lp,ndim_magmom) )
    allocate( dmmat_i_magmom(max2lp,max2lp,ndim_magmom) )
    allocate( dmmat_ssrep(max2lp,max2lp,ndim_chgpot) )

    do ia=1,natm
       ih = ihubbard(ia)
       if(ih /= 0) then
          ie = proj_attribute(ih)%ielem
          it=ityp(ia)

          dmmat_r_magmom = 0.0d0
          dmmat_i_magmom = 0.0d0
          dmmat_ssrep = 0.0d0

          if ( sw_eval_energy_before_charge == ON ) then
             dmmat_r_magmom(:,:,:) = ommix(:,:,ie,ia,:)
             dmmat_i_magmom(:,:,:) = ommix_aimag(:,:,ie,ia,:)
          else
             dmmat_r_magmom(:,:,:) = om(:,:,ie,ia,:)
             dmmat_i_magmom(:,:,:) = om_aimag(:,:,ie,ia,:)
          endif
          call m_ES_MagMom_To_DensMat_porb( max2lp**2, dmmat_r_magmom, dmmat_i_magmom, &
               &                            dmmat_ssrep )

          dsum=0.d0
          Do is=1, ndim_chgpot
             if ( is==1 .or. is==ndim_chgpot ) then
                do i2=1,i2lp(ih)
                   dsum = dsum + dmmat_ssrep(i2,i2,is)
                end do
             endif
          End do
          Do is=1, ndim_chgpot
             do i2=1,i2lp(ih)
                do i1=1,i2lp(ih)
                   dsum = dsum - dmmat_ssrep(i1,i2,is) &
                        &        *conjg( dmmat_ssrep(i1,i2,is) )
                end do
             end do
          End do

          energy = energy + 0.5d0*proj_attribute(ih)%Ueff *real(dsum)

       end if
    end do

    energy = energy *Ueff_prefactor

    deallocate( dmmat_r_magmom, dmmat_i_magmom )
    deallocate( dmmat_ssrep )

  end subroutine m_Hubbard_energy2_noncl

  subroutine m_Hubbard_energy3_noncl(energy)           ! for safety
    real(kind=DP), intent(out) :: energy

    real(kind=DP), allocatable :: dmmat_r_magmom( :,:,: )
    real(kind=DP), allocatable :: dmmat_i_magmom( :,:,: )
    complex(kind=CMPLDP), allocatable :: dmmat_ssrep( :,:,: )
    complex(kind=CMPLDP), allocatable :: dm_large( :,: )

    integer :: ia, it,i1,i2,ih,ie
    integer :: is1, is2,  size1
    complex(kind=CMPLDP) :: dsum

    if(sw_constraint == ON) return

    energy=0.d0

    allocate( dmmat_r_magmom(max2lp,max2lp,ndim_magmom) )
    allocate( dmmat_i_magmom(max2lp,max2lp,ndim_magmom) )
    allocate( dmmat_ssrep(max2lp,max2lp,ndim_chgpot) )

    allocate( dm_large( max2lp*ndim_spinor, max2lp*ndim_spinor ) )
    dm_large = 0.0d0

    do ia=1,natm
       ih = ihubbard(ia)
       if(ih /= 0) then
          ie = proj_attribute(ih)%ielem
          it=ityp(ia)

          dmmat_r_magmom = 0.0d0
          dmmat_i_magmom = 0.0d0
          dmmat_ssrep = 0.0d0

          if ( sw_eval_energy_before_charge == ON ) then
             dmmat_r_magmom(:,:,:) = ommix(:,:,ie,ia,:)
             dmmat_i_magmom(:,:,:) = ommix_aimag(:,:,ie,ia,:)
          else
             dmmat_r_magmom(:,:,:) = om(:,:,ie,ia,:)
             dmmat_i_magmom(:,:,:) = om_aimag(:,:,ie,ia,:)
          endif
          call m_ES_MagMom_To_DensMat_porb( max2lp**2, dmmat_r_magmom, dmmat_i_magmom, &
               &                            dmmat_ssrep )

          dm_large = 00.0d0

          Do is1=1, ndim_spinor
             Do is2=1, ndim_spinor
                do i1=1,i2lp(ih)
                   do i2=1,i2lp(ih)
                      dm_large( max2lp*(is1-1) +i1, max2lp*(is2-1)+i2 ) &
                           &  = dmmat_ssrep( i1,i2, (is1-1)*ndim_spinor +is2 )
                   end do
                end do
             end do
          end do

          dsum=0.d0
          Do i1=1, ndim_spinor *max2lp
             dsum = dsum + dm_large( i1,i1 )
          End do
          Do i1=1, ndim_spinor *max2lp
             Do i2=1, ndim_spinor *max2lp
                dsum = dsum - dm_large( i1,i2 )*dm_large(i2,i1)
             end do
          End do

          energy = energy + 0.5d0*proj_attribute(ih)%Ueff *real(dsum)

       end if
    end do

    energy = energy *Ueff_prefactor

    deallocate( dm_large )
    deallocate( dmmat_r_magmom, dmmat_i_magmom )
    deallocate( dmmat_ssrep )

  end subroutine m_Hubbard_energy3_noncl
! ===================================================================== 11.0

  subroutine m_Hubbard_Potential(nfout)
    integer, intent(in) :: nfout
    integer :: is,ia,iyy
    integer :: ilmt1,m1,t1,ilmt2,m2,t2
    integer :: ih,it,ie,t0
    real(kind=DP) :: Ueff

    dhub = 0.d0
    do is=1,nspin,af+1
       do ia=1,natm
          if(sw_constraint == ON .and. ia /= const_site) cycle
          ih = ihubbard(ia)
          if(ih == 0) cycle
          Ueff = proj_attribute(ih)%Ueff
          it = proj_attribute(ih)%ityp
          ie = proj_attribute(ih)%ielem
          if(projector_type == SPHERICAL_HARMONICS) then
             do iyy=1,nyy(ih)
                ilmt1=ilmt_yy(1,iyy,ih)
                m1=mtp(ilmt1,it)
                t1=taup(ilmt1,it)
                ilmt2=ilmt_yy(2,iyy,ih)
                m2=mtp(ilmt2,it)
                t2=taup(ilmt2,it)
                !!$write(nfout,'("m1,t1,m2,t2=",4i3)') m1,t1,m2,t2
                if(sw_constraint == ON) then
                   if(m1==m2) then
                      dhub(ilmt1,ilmt2,ia,is)=const_alpha*prodphi(ih,t1,t2)
                   end if
                else
                   if(m1==m2) then
                      dhub(ilmt1,ilmt2,ia,is)=Ueff &
                      & *(0.5d0-ommix(m1,m2,ie,ia,is))*prodphi(ih,t1,t2)
                   else
                      dhub(ilmt1,ilmt2,ia,is)=-Ueff &
                      & *ommix(m1,m2,ie,ia,is)*prodphi(ih,t1,t2)
                   end if
                end if
                dhub(ilmt2,ilmt1,ia,is)=dhub(ilmt1,ilmt2,ia,is)
                !!$write(nfout,'("ilmt1,ilmt2,ia,is,dhub=",4i3,e12.5)') ilmt1,ilmt2,ia,is,dhub(ilmt1,ilmt2,ia,is)
                !!$write(nfout,'("ueff,ommix,prodphi=",3e12.5)') Ueff,ommix(m1,m2,ie,ia,is),prodphi(ih,t1,t2)
             end do
          else if(projector_type == ATOMIC_ORBITAL) then
             t0 = proj_attribute(ih)%t
             do iyy=1,nyy(ih)
                ilmt1=ilmt_yy(1,iyy,ih)
                m1=mtp(ilmt1,it)
                t1=taup(ilmt1,it)
                ilmt2=ilmt_yy(2,iyy,ih)
                m2=mtp(ilmt2,it)
                t2=taup(ilmt2,it)
                !!$write(nfout,'("m1,t1,m2,t2=",4i3)') m1,t1,m2,t2
                if(sw_constraint == ON) then
                   if(m1==m2) then
                      dhub(ilmt1,ilmt2,ia,is)=const_alpha*prodphi(ih,t1,t0)*prodphi(ih,t0,t2)
                   end if
                else
                   if(m1==m2) then
                      dhub(ilmt1,ilmt2,ia,is)=Ueff &
                      & *(0.5d0-ommix(m1,m2,ie,ia,is))*prodphi(ih,t1,t0)*prodphi(ih,t0,t2)
                   else
                      dhub(ilmt1,ilmt2,ia,is)=-Ueff &
                      & *ommix(m1,m2,ie,ia,is)*prodphi(ih,t1,t0)*prodphi(ih,t0,t2)
                   end if
                end if
                dhub(ilmt2,ilmt1,ia,is)=dhub(ilmt1,ilmt2,ia,is)
                !!$write(nfout,'("ilmt1,ilmt2,ia,is,dhub=",4i3,e12.5)') ilmt1,ilmt2,ia,is,dhub(ilmt1,ilmt2,ia,is)
                !!$write(nfout,'("ie,ueff,ommix,prodphi=",i3,3e12.5)') ie,Ueff,ommix(m1,m2,ie,ia,is),prodphi(ih,t1,t2)
             end do
          end if
       end do
    end do

    if(iprihubbard >2) then
       write(nfout,'("m_Hubbard_Potential")')
       do is=1,nspin,af+1
          do ia=1,natm
             ih=ihubbard(ia)
             if(ih==0) cycle
             it=ityp(ia)
             write(nfout,'("is,ia,it=",3i3)') is,ia,it
             do ilmt1=1,ilmt(it)
                write(nfout,'(14f8.3)') (dhub(ilmt1,ilmt2,ia,is),ilmt2=1,ilmt(it))
             end do
          end do
       end do
    end if
  end subroutine m_Hubbard_Potential

! =============================- added by K. Tagami ================ 5.0
  subroutine m_Hubbard_Potential2(nfout)
    integer, intent(in) :: nfout
    integer :: is,ia,iyy
    integer :: ilmt1,m1,t1,ilmt2,m2,t2
    integer :: ih,it,ie,t0
    real(kind=DP) :: Ueff, ctmp, c1, nsigma
!
    integer :: ilmt3, ilmt4, i1
    integer :: l1, l2, l3, l4, m3, m4, t3, t4, ilp
!
    if(.not.allocated(dhub)) then
       write(nfout,'(" dhub is not allocated (subroutine m_Hubbard_Potential2)")')
       call flush(nfout)
       call phase_error_with_msg(nfout,' dhub is not allocated (subroutine m_Hubbard_Potential2)',__LINE__,__FILE__)
    end if
    dhub = 0.d0
    do is=1,nspin,af+1
       do ia=1,natm
          if(sw_constraint == ON .and. ia /= const_site) cycle
          ih = ihubbard(ia)
          if(ih == 0) cycle
          Ueff = proj_attribute(ih)%Ueff

          it = proj_attribute(ih)%ityp
          ie = proj_attribute(ih)%ielem
          
          t0 = proj_attribute(ih)%t

          if ( occmat_type == OCCMAT_Type1 ) then
             select case(dftu_type)
             case (0)
                do iyy=1,nyy(ih)
                   ilmt1=ilmt_yy(1,iyy,ih); m1=mtp(ilmt1,it); t1=taup(ilmt1,it)
                   ilmt2=ilmt_yy(2,iyy,ih); m2=mtp(ilmt2,it); t2=taup(ilmt2,it)
                   if(m1==m2) then
                      dhub(ilmt1,ilmt2,ia,is)=const_alpha*prodphi(ih,t1,t2)
                   end if
                   dhub(ilmt2,ilmt1,ia,is)=dhub(ilmt1,ilmt2,ia,is)
                end do

             case (FLL)
                do iyy=1,nyy(ih)
                   ilmt1=ilmt_yy(1,iyy,ih); m1=mtp(ilmt1,it); t1=taup(ilmt1,it)
                   ilmt2=ilmt_yy(2,iyy,ih); m2=mtp(ilmt2,it); t2=taup(ilmt2,it)
                   if(m1==m2) then
                      dhub(ilmt1,ilmt2,ia,is)=Ueff &
                           & *(0.5d0-ommix(m1,m2,ie,ia,is))*prodphi(ih,t1,t2)
                   else
                      dhub(ilmt1,ilmt2,ia,is)=-Ueff &
                           & *ommix(m1,m2,ie,ia,is)*prodphi(ih,t1,t2)
                   end if
                   dhub(ilmt2,ilmt1,ia,is)=dhub(ilmt1,ilmt2,ia,is)
                end do

             case (AMF)
                c1 = 0.0d0
                do i1=1,i2lp(ih)
                   c1 = c1 +ommix(i1,i1,ie,ia,is)
               End do
                nsigma = c1 /( 2*proj_attribute(ih)%l +1 )

                do iyy=1,nyy(ih)
                   ilmt1=ilmt_yy(1,iyy,ih); m1=mtp(ilmt1,it); t1=taup(ilmt1,it)
                   ilmt2=ilmt_yy(2,iyy,ih); m2=mtp(ilmt2,it); t2=taup(ilmt2,it)
!
                   if(m1==m2) then
                      dhub(ilmt1,ilmt2,ia,is)=Ueff &
                           & *( nsigma -ommix(m1,m2,ie,ia,is))*prodphi(ih,t1,t2)
                   else
                      dhub(ilmt1,ilmt2,ia,is)=-Ueff &
                           & *ommix(m1,m2,ie,ia,is)*prodphi(ih,t1,t2)
                   end if
                   dhub(ilmt2,ilmt1,ia,is)=dhub(ilmt1,ilmt2,ia,is)
                end do

             end select

          else if ( occmat_type == OCCMAT_Type2 ) then
             select case(dftu_type)
             case (0)
                t0 = proj_attribute(ih)%t
                do iyy=1,nyy(ih)
                   ilmt1=ilmt_yy(1,iyy,ih); m1=mtp(ilmt1,it); t1=taup(ilmt1,it)
                   ilmt2=ilmt_yy(2,iyy,ih); m2=mtp(ilmt2,it); t2=taup(ilmt2,it)
                   ctmp = prodphi(ih,t1,t0)*prodphi(ih,t0,t2)
                   if(m1==m2) then
                      dhub(ilmt1,ilmt2,ia,is)=const_alpha *ctmp
                   end if
                end do

             case (FLL)
                t0 = proj_attribute(ih)%t
                do iyy=1,nyy(ih)
                   ilmt1=ilmt_yy(1,iyy,ih); m1=mtp(ilmt1,it); t1=taup(ilmt1,it)
                   ilmt2=ilmt_yy(2,iyy,ih); m2=mtp(ilmt2,it); t2=taup(ilmt2,it)
                   ctmp = prodphi(ih,t1,t0)*prodphi(ih,t0,t2)

                   if(m1==m2) then
                      dhub(ilmt1,ilmt2,ia,is)=Ueff &
                           & *(0.5d0-ommix(m1,m2,ie,ia,is)) *ctmp
                   else
                      dhub(ilmt1,ilmt2,ia,is)=-Ueff &
                           & *ommix(m1,m2,ie,ia,is)*ctmp
                   end if
                   dhub(ilmt2,ilmt1,ia,is)=dhub(ilmt1,ilmt2,ia,is)
                end do

             end select

          end if
       end do
    end do
! ======================================= addd by K. Tagami ============= 5.2
    dhub = dhub *Ueff_prefactor
! ======================================================================= 5.2

    if(iprihubbard >2) then
       write(nfout,'("m_Hubbard_Potential")')
       do is=1,nspin,af+1
          do ia=1,natm
             ih=ihubbard(ia)
             if(ih==0) cycle
             it=ityp(ia)
             write(nfout,'("is,ia,it=",3i3)') is,ia,it
             do ilmt1=1,ilmt(it)
                write(nfout,'(14f8.3)') (dhub(ilmt1,ilmt2,ia,is),ilmt2=1,ilmt(it))
              end do
          end do
       end do
    end if
  end subroutine m_Hubbard_Potential2
! ============================================================== 5.0

! ================== added by K. Tagami ================= 11.0
  subroutine m_Hubbard_Potential2_noncl(nfout)          
!                                        This does not work properly.
    integer, intent(in) :: nfout
    integer :: is,ia,iyy
    integer :: ilmt1,m1,t1,ilmt2,m2,t2
    integer :: ih,it,ie,t0
    real(kind=DP) :: Ueff, ctmp
! 
    integer :: ilmt3, ilmt4
    integer :: l1, l2, l3, l4, m3, m4, t3, t4, ilp
!
    dhub = 0.d0
    do ia=1,natm
       if(sw_constraint == ON .and. ia /= const_site) cycle

       ih = ihubbard(ia)
       if(ih == 0) cycle

       Ueff = proj_attribute(ih)%Ueff

       it = proj_attribute(ih)%ityp
       ie = proj_attribute(ih)%ielem
          
       t0 = proj_attribute(ih)%t

       Do is=1, ndim_magmom
          if ( occmat_type == OCCMAT_Type1 ) then
             do iyy=1,nyy(ih)
                ilmt1=ilmt_yy(1,iyy,ih); m1=mtp(ilmt1,it); t1=taup(ilmt1,it)
                ilmt2=ilmt_yy(2,iyy,ih); m2=mtp(ilmt2,it); t2=taup(ilmt2,it)
                
                if (sw_constraint == ON) then        
                   if ( is==1 .and. m1==m2) then     ! uncertain
                      dhub(ilmt1,ilmt2,ia,is)=const_alpha*prodphi(ih,t1,t2)
                   end if
                else
                   if ( is == 1 .and. m1==m2 ) then
                      dhub(ilmt1,ilmt2,ia,is)=Ueff &
                           & *(0.5d0-ommix(m1,m2,ie,ia,is))*prodphi(ih,t1,t2)
                   else
                      dhub(ilmt1,ilmt2,ia,is)=-Ueff &
                           & *ommix(m1,m2,ie,ia,is)*prodphi(ih,t1,t2)
                   end if
                end if
                dhub(ilmt2,ilmt1,ia,is)=dhub(ilmt1,ilmt2,ia,is)
             end do
          else if ( occmat_type == OCCMAT_Type2 ) then
             t0 = proj_attribute(ih)%t
             do iyy=1,nyy(ih)
                ilmt1=ilmt_yy(1,iyy,ih); m1=mtp(ilmt1,it); t1=taup(ilmt1,it)
                ilmt2=ilmt_yy(2,iyy,ih); m2=mtp(ilmt2,it); t2=taup(ilmt2,it)
                ctmp = prodphi(ih,t1,t0)*prodphi(ih,t0,t2)

                if (sw_constraint == ON) then
                   if ( is==1 .and. m1==m2 ) then            ! uncertain
                      dhub(ilmt1,ilmt2,ia,is)=const_alpha *ctmp
                   end if
                else
                   if ( is==1 .and. m1==m2 ) then
                      dhub(ilmt1,ilmt2,ia,is)=Ueff &
                           & *(0.5d0-ommix(m1,m2,ie,ia,is)) *ctmp
                   else
                      dhub(ilmt1,ilmt2,ia,is)=-Ueff &
                           & *ommix(m1,m2,ie,ia,is)*ctmp
                   end if
                end if
                dhub(ilmt2,ilmt1,ia,is)=dhub(ilmt1,ilmt2,ia,is)
             end do
          end if
       end do
    end do

    dhub = dhub *Ueff_prefactor

! --- for debug --
!    Do is=1, ndim_magmom
!       write(920,*) 'is = ', is
!       Do ia=1, natm
!          write(920,*) 'ia = ' , ia
!          Do ilmt1=1, ilmt(it)
!             Do ilmt2=1, ilmt(it)
!                write(920,*) 'ilmt1 ilmt2 = ' ,ilmt1, ilmt2
!                write(920,*) dhub(ilmt1,ilmt2,ia,is)
!             End do
!          End do
!       End do
!    End do
! ---------------

    if (iprihubbard >2) then
       write(nfout,'("m_Hubbard_Potential2_noncl")')
       do is=1,ndim_magmom,af+1
          do ia=1,natm
             ih=ihubbard(ia)
             if(ih==0) cycle
             it=ityp(ia)
             write(nfout,'("is,ia,it=",3i3)') is,ia,it
             do ilmt1=1,ilmt(it)
                write(nfout,'(14f8.3)') (dhub(ilmt1,ilmt2,ia,is),ilmt2=1,ilmt(it))
              end do
          end do
       end do
    end if
  end subroutine m_Hubbard_Potential2_noncl

  subroutine m_Hubbard_Potential3_noncl(nfout)
    integer, intent(in) :: nfout
    integer :: is,ia,iyy
    integer :: ilmt1,m1,t1,ilmt2,m2,t2
    integer :: ih,it,ie,t0
    real(kind=DP) :: Ueff, ctmp
! 
    integer :: ilmt3, ilmt4
    integer :: l1, l2, l3, l4, m3, m4, t3, t4, ilp

    real(kind=DP), allocatable :: dmmat_r_magmom( :,:,: )
    real(kind=DP), allocatable :: dmmat_i_magmom( :,:,: )
    complex(kind=CMPLDP), allocatable :: dmmat_ssrep( :,:,: )
!
    complex(kind=CMPLDP), allocatable :: dhub_ssrep( :,:,:,: )
!

    dhub = 0.d0

    allocate( dhub_ssrep(nlmt,nlmt,natm,ndim_magmom) )
    dhub_ssrep = 0.0d0

    allocate( dmmat_r_magmom(max2lp,max2lp,ndim_magmom) )
    allocate( dmmat_i_magmom(max2lp,max2lp,ndim_magmom) )
    allocate( dmmat_ssrep(max2lp,max2lp,ndim_chgpot) )

    do ia=1,natm
       if(sw_constraint == ON .and. ia /= const_site) cycle

       ih = ihubbard(ia)
       if(ih == 0) cycle

       Ueff = proj_attribute(ih)%Ueff

       it = proj_attribute(ih)%ityp
       ie = proj_attribute(ih)%ielem
          
       t0 = proj_attribute(ih)%t
       ilp = proj_attribute(ih)%l +1

       dmmat_r_magmom = 0.0d0;   dmmat_i_magmom = 0.0d0;  dmmat_ssrep = 0.0d0
       dmmat_r_magmom(:,:,:) = ommix(:,:,ie,ia,:)
       dmmat_i_magmom(:,:,:) = ommix_aimag(:,:,ie,ia,:)

       call m_ES_MagMom_To_DensMat_porb( max2lp**2, dmmat_r_magmom, dmmat_i_magmom, &
            &                            dmmat_ssrep )

       Do is=1, ndim_chgpot
          if ( occmat_type == OCCMAT_Type1 ) then
             select case (dftu_type)
             case (0)
                call phase_error_with_msg(nfout,'kt : not supported. ',__LINE__,__FILE__)

             case (FLL)
                do ilmt1=1, ilmt(it)
                   do ilmt2=1, ilmt(it)
                      l1 = ltp(ilmt1,it); m1 = mtp(ilmt1,it); t1 = taup(ilmt1,it)
                      l2 = ltp(ilmt2,it); m2 = mtp(ilmt2,it); t2 = taup(ilmt2,it)

                      if ( l1 /= l2 .or. l1 /= ilp ) cycle

                      if ( (is == 1 .or. is==ndim_chgpot) .and. m1==m2 ) then
                         dhub_ssrep(ilmt1,ilmt2,ia,is) = Ueff &
                              & *( 0.5d0 - dmmat_ssrep(m1,m2,is) )*prodphi(ih,t1,t2)
                      else
                         dhub_ssrep(ilmt1,ilmt2,ia,is) = -Ueff &
                              & *dmmat_ssrep(m1,m2,is) *prodphi(ih,t1,t2)
                         !                           & *conjg( dmmat_ssrep(m1,m2,is) ) *prodphi(ih,t1,t2)
                      end if
                   end do
                end do

             end select

          else if ( occmat_type == OCCMAT_Type2 ) then
             t0 = proj_attribute(ih)%t

             select case (dftu_type)
             case (0)
                call phase_error_with_msg(nfout,'kt : not supported. ',__LINE__,__FILE__)

             case (FLL)
                do ilmt1=1, ilmt(it)
                   do ilmt2=1, ilmt(it)
                      l1 = ltp(ilmt1,it); m1 = mtp(ilmt1,it); t1 = taup(ilmt1,it)
                      l2 = ltp(ilmt2,it); m2 = mtp(ilmt2,it); t2 = taup(ilmt2,it)

                      if ( l1 /= l2 .or. l1 /= ilp ) cycle
                      ctmp = prodphi(ih,t1,t0)*prodphi(ih,t0,t2)

                      if ( (is == 1 .or. is==ndim_chgpot) .and. m1==m2 ) then
                         dhub_ssrep(ilmt1,ilmt2,ia,is) = Ueff &
                              & *( 0.5d0 - dmmat_ssrep(m1,m2,is) )*ctmp
                      else
                         dhub_ssrep(ilmt1,ilmt2,ia,is) = -Ueff &
                              & *dmmat_ssrep(m1,m2,is) *ctmp
                      end if
                   end do
                end do
             end select

          end if
       end do
    end do
!
    call m_ES_DensMat_to_MagMom_Dhub( dhub_ssrep, dhub, dhub_aimag )
!
    dhub = dhub *Ueff_prefactor
    dhub_aimag = dhub_aimag *Ueff_prefactor
!
    goto 1200

! --- for debug --
    Do is=1, ndim_magmom
       write(920,*) 'is = ', is
       Do ia=1, natm
          write(920,*) 'ia = ' , ia
          Do ilmt1=1, ilmt(it)
             Do ilmt2=1, ilmt(it)
                write(920,*) 'ilmt1 ilmt2 = ' ,ilmt1, ilmt2
                write(920,*) dhub_ssrep(ilmt1,ilmt2,ia,is)
             End do
          End do
       End do
    End do
! ---------------

    do ia=1,natm
       do is=1, ndim_magmom

          Do ilmt1=1, ilmt(it)
             Do ilmt2=1, ilmt(it)
                write(960,*) ' ia lmt1 lmtm2 istmp = ', ia, ilmt1, ilmt2, is
                write(960,*) 'A ', dhub(ilmt1,ilmt2,ia,is), dhub(ilmt2,ilmt1,ia,is)
                write(960,*) 'B ', dhub_aimag(ilmt1,ilmt2,ia,is), dhub_aimag(ilmt2,ilmt1,ia,is)
             end do
          end do
       end do
    end do

1200 continue

    if (iprihubbard >2) then
       write(nfout,'("m_Hubbard_Potential3_noncl")')
       do is=1,ndim_magmom,af+1
          do ia=1,natm
             ih=ihubbard(ia)
             if(ih==0) cycle
             it=ityp(ia)
             write(nfout,'("is,ia,it=",3i3)') is,ia,it
             do ilmt1=1,ilmt(it)
                write(nfout,'(14f8.3)') (dhub(ilmt1,ilmt2,ia,is),ilmt2=1,ilmt(it))
              end do
          end do
       end do
    end if

    deallocate( dhub_ssrep )
    deallocate( dmmat_ssrep )
    deallocate( dmmat_r_magmom, dmmat_i_magmom )

  end subroutine m_Hubbard_Potential3_noncl
! ============================================================== 11.0

end module m_Hubbard
