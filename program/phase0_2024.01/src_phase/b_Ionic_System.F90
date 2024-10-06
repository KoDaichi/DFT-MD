!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  SUBROUINE: forcrsv, forc_cnst, copy_to_cpd_forc_and_forcp, 
!             cnstrnt, rplcps, rbinuc, change_of_coordinate_system
!
!  AUTHOR(S): T. Uchiyama, T. Yamasaki   August/20/2003
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
! $Id: b_Ionic_System.f90 633 2020-12-01 05:11:03Z jkoga $
subroutine forcrsv(natm,ityp,imdtyp,iwei,amion,cpd_l,tkb,nrsv,imdalg &
     & ,frsv,iw_cnst)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  SUBROUINE: forcrsv
!
!  AUTHOR(S): T. Uchiyama, T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP,BLUEMOON,QUENCHED_CONSTRAINT,ON
  use m_Ionic_System,     only : nbonds_per_thermo, sw_fix_bond
  implicit none
  integer, intent(in)       :: natm,nrsv,imdalg &
       &                      ,ityp(natm),imdtyp(natm),iwei(natm)
  real(kind=DP),intent(in)  :: amion(*),tkb(nrsv),cpd_l(natm,3)
  real(kind=DP),intent(out) :: frsv(nrsv)
  integer, intent(in) :: iw_cnst
  real(kind=DP) :: vex,vey,vez
  integer :: ia, ir, icnstrnt_typ

  frsv = 0.d0
  do ia = 1, natm
     ir = icnstrnt_typ(imdtyp(ia),imdalg)
     if(ir >= 1) then
          vex = cpd_l(ia,1)
          vey = cpd_l(ia,2)
          vez = cpd_l(ia,3)
          frsv(ir)= frsv(ir) &
     &             +(( vex*vex +vey*vey +vez*vez )*amion(ityp(ia))&
     &                -3.d0*tkb(ir))*iwei(ia)
     end if
  end do
  if(sw_fix_bond==ON) then
    do ir=1,nrsv
      frsv(ir) = frsv(ir) + tkb(ir)*nbonds_per_thermo(ir)
    enddo
  endif

  if(imdalg == BLUEMOON .or. imdalg == QUENCHED_CONSTRAINT) then
     ir = 1
     frsv(ir) = frsv(ir) + iw_cnst*tkb(ir)
  end if
end subroutine forcrsv

integer function irtyp(imd)
  implicit none
  integer, intent(in) :: imd
  if(imd <= 0) irtyp = -1
  if(imd >= 1) irtyp =  0
  if(imd > 1000) irtyp = imd - 1000
end function irtyp

integer function ibath(imd)
  integer, intent(in) :: imd
  if(imd <= 0) ibath = -1
  if(imd >= 1) ibath =  0
  if(imd > 1000) ibath = imd/1000
end function ibath

integer function icnstrnt_typ(imd,imdalg)
  use m_Const_Parameters, only : HEAT_BATH, BLUEMOON, QUENCHED_CONSTRAINT, VERLET
  implicit none
  integer, intent(in) :: imd, imdalg

  if(imdalg == VERLET) then
     if(imd > HEAT_BATH) then
        icnstrnt_typ = imd - HEAT_BATH
     else
        icnstrnt_typ = imd
     end if
  else
     if(imd <= 0) icnstrnt_typ = -1
     if(imd >= 1) icnstrnt_typ =  0
     if(imd > HEAT_BATH) then
        if(imdalg == BLUEMOON .or. imdalg == QUENCHED_CONSTRAINT) then
           icnstrnt_typ = imd/HEAT_BATH
        else  ! imdalg == TEMPERATURE_CONTROL etc.
           icnstrnt_typ = imd - HEAT_BATH
        end if
     end if
  end if
end function icnstrnt_typ

subroutine forc_cnst(natm,ityp,amion,cps,cnst_typ,nfcatm,ia_cnst,imdtyp &
     &, sgmc,gca)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  SUBROUINE: forc_cnst
!
!  AUTHOR(S): T. Uchiyama, T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP,BONDLENGTH_FIX_1, BONDLENGTH_FIX_2 &
       &, FIX, COG_FIX_L, BONDLENGTH_FIX
  implicit none
  integer, intent(in) :: natm,ityp(natm),cnst_typ,nfcatm,ia_cnst(nfcatm)&
       &,  imdtyp(natm)
  real(kind=DP),intent(in) :: amion(*),cps(natm,3),sgmc(4)
  real(kind=DP),intent(out):: gca(natm,3)

  integer       :: ia, ia1, ia2, m
  real(kind=DP) :: r12, tmass,rmass

  gca = 0.d0

  ! +++ bond length (1st order) +++
  if(cnst_typ ==  BONDLENGTH_FIX_1 .or. cnst_typ == BONDLENGTH_FIX) then
     ia1= ia_cnst(1)
     ia2= ia_cnst(2)

     r12= dsqrt(( cps(ia1,1) -cps(ia2,1) )**2  &
          &    +( cps(ia1,2) -cps(ia2,2) )**2  &
          &    +( cps(ia1,3) -cps(ia2,3) )**2 )
     gca(ia1,1:3) = ( cps(ia1,1:3) - cps(ia2,1:3))/r12
     gca(ia2,1:3) = - gca(ia1,1:3)
     return
  end if

  ! +++ bond length (2nd order) +++
  if(cnst_typ ==  BONDLENGTH_FIX_2) then
     ia1= ia_cnst(1)
     ia2= ia_cnst(2)

     gca(ia1,1:3) = 2.d0 * (cps(ia1,1:3) - cps(ia2,1:3))
     gca(ia2,1:3)= -gca(ia1,1:3)
     return
  end if

  ! +++ center-of-mass +++
  if(cnst_typ == COG_FIX_L) then

     tmass= 0.d0; rmass = 0.d0
     do m = 1, nfcatm
        tmass = tmass + amion(ityp(ia_cnst(m)))
     end do
     do ia = 1, natm
        if(imdtyp(ia) /= FIX) rmass = rmass + amion(ityp(ia))
     end do
     rmass = rmass - tmass

     do ia = 1, natm
        if(imdtyp(ia) /= FIX) then
           gca(ia,1:3) = - amion(ityp(ia))/rmass * sgmc(1:3)
        end if
     end do
     do m = 1, nfcatm
        ia = ia_cnst(m)
        gca(ia,1:3) = amion(ityp(ia))/tmass * sgmc(1:3)
     end do
     return
  end if

  print *,'!(forc_cnst) Undefined cnst_typ =', cnst_typ
  call phase_error_with_msg(6, '!(forc_cnst) Undefined cnst_typ',__LINE__,__FILE__)
  !stop
end subroutine forc_cnst

subroutine copy_to_cpd_forc_and_forcp(katm,natm,cpd_l,forc_l,cpd,forc,forcp)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  SUBROUINE: copy_to_cpd_forc_and_forcp
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP
  implicit none

  integer, intent(in) :: katm,natm
  real(kind=DP),intent(in),dimension(katm,3) :: cpd_l,forc_l
  real(kind=DP),intent(out),dimension(katm,3):: cpd,forc,forcp

  integer :: ia
  cpd = 0.d0; forc = 0.d0
  do ia = 1, natm
     cpd(ia,1:3)   = cpd_l(ia,1:3)
     forc(ia,1:3)  = forc_l(ia,1:3)
     forcp(ia,1:3) = forc_l(ia,1:3)
  end do
  do ia = 1, natm
     print *,' cpd = ', cpd_l(ia,1),cpd_l(ia,2),cpd_l(ia,3)
  end do
end subroutine copy_to_cpd_forc_and_forcp

subroutine cnstrnt(natm,ityp,amion,cps,cnst_typ,nfcatm,ia_cnst,imdtyp &
     &, sgmc,sigma)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  SUBROUINE: cnstrnt
!
!  AUTHOR(S): T. Uchiyama,  T. Yamasaki   August/20/2003
!
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================


  use m_Const_Parameters, only : DP, FIX &
       & ,   BONDLENGTH_FIX, BONDLENGTH_FIX_1, BONDLENGTH_FIX_2, COG_FIX_L
  implicit none
  integer, intent(in)      :: natm,ityp(natm),cnst_typ,nfcatm,ia_cnst(nfcatm)
  real(kind=DP),intent(in) :: amion(*),cps(natm,3),sgmc(4),imdtyp(natm)
  real(kind=DP),intent(out):: sigma

  integer       :: ia1,ia2,ia,j,m
  real(kind=DP) :: tmass, rmass, gmass

  ! +++ bond length (1st order) +++
  if(cnst_typ == BONDLENGTH_FIX_1 .or. cnst_typ == BONDLENGTH_FIX) then

     ia1= ia_cnst(1)
     ia2= ia_cnst(2)
     sigma = dsqrt(( cps(ia1,1) -cps(ia2,1) )**2  &
          &       +( cps(ia1,2) -cps(ia2,2) )**2  &
          &       +( cps(ia1,3) -cps(ia2,3) )**2 ) -sgmc(1)
     return
  end if

  ! +++ bond length (2nd order) +++
  if(cnst_typ == BONDLENGTH_FIX_2) then
     ia1= ia_cnst(1)
     ia2= ia_cnst(2)
     sigma = ( cps(ia1,1) -cps(ia2,1) )**2  &
          & +( cps(ia1,2) -cps(ia2,2) )**2  &
          & +( cps(ia1,3) -cps(ia2,3) )**2 -sgmc(1)
     return
  end if

  ! +++ center-of-mass +++
  if(cnst_typ == COG_FIX_L) then
     tmass = 0.d0
     rmass = 0.d0
     do m= 1, nfcatm
        tmass = tmass + amion (ityp(ia_cnst(m)))
     end do
     do ia = 1, natm
        if(imdtyp(ia) /= FIX) rmass = rmass + amion(ityp(ia))
     end do
     rmass = rmass - tmass
     gmass = tmass*rmass /(tmass + rmass)
     sigma= 0.d0
     do ia = 1, natm
        if(imdtyp(ia) /= FIX) then
           do j = 1, 3
              sigma = sigma - amion(ityp(ia))/rmass * cps(ia,j)*sgmc(j)
           end do
        end if
     end do
     do m= 1, nfcatm
        ia= ia_cnst(m)
        do j = 1, 3
           sigma = sigma +amion(ityp(ia))/gmass *cps(ia,j)*sgmc(j)
        end do
     end do
     sigma = sigma -sgmc(4)
     return
  end if

  print *, '!(cnstrnt) Undefined cnst_typ=', cnst_typ
  call phase_error_with_msg(6, '!(cnstrnt) Undefined cnst_typ',__LINE__,__FILE__)
  !stop
end subroutine cnstrnt

subroutine rplcps(cps,ityp,mode,natm2,natm,iwei)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  SUBROUINE: rplcps
!
!  AUTHOR(S): T. Uchiyama, T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  implicit none
  integer      , intent(in)                        :: mode
  integer      , intent(in)                        :: natm2
  integer      , intent(in)                        :: natm
  integer      , intent(in)   , dimension(natm)    :: iwei
  real(kind=kind(1.d0)), intent(inout), dimension(natm2,3) :: cps
  integer,       intent(inout), dimension(natm2)   :: ityp
  
  integer nb, ia
  
  nb = 0
  do ia = 1, natm
     if(iwei(ia) ==  2) then
        nb = nb + 1
        cps(natm+nb,1) = -cps(ia,1)
        cps(natm+nb,2) = -cps(ia,2)
        cps(natm+nb,3) = -cps(ia,3)
        if(mode.eq.1) ityp(natm+nb) = ityp(ia)
     endif
  enddo
  
end subroutine rplcps

subroutine rbinuc(pos,natm)
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  SUBROUINE: rbinuc
!
!  AUTHOR(S): T. Uchiyama, T. Yamasaki   June/16/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in) ::                             natm
  real(kind=DP), intent(inout), dimension(natm,3) :: pos

  integer :: i,j
  do j = 1, 3
     do i = 1, natm
        if(pos(i,j) < 0.d0) pos(i,j) = (ceiling(-pos(i,j))) + pos(i,j)
        if(pos(i,j) > 1.d0) pos(i,j) = pos(i,j) - floor(pos(i,j))
     end do
  end do
end subroutine rbinuc

subroutine change_of_coordinate_system(t,v,k,n,w)
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in) :: k, n
  real(kind=DP), intent(in),  dimension(3,3) :: t
  real(kind=DP), intent(in),  dimension(k,3) :: v
  real(kind=DP), intent(out), dimension(k,3) :: w
  integer i, j
  do j = 1, 3
     do i = 1, n
        w(i,j) = t(j,1)*v(i,1) + t(j,2)*v(i,2) + t(j,3)*v(i,3)
     enddo
  enddo
end subroutine change_of_coordinate_system

subroutine fd_symmetrize(natm2,katm,natm,napt,kopr,nopr,op,iwei &
     &, f_l, f_wk, npfatm)
! This subroutine is for symmetrization of atomic forces.
! Original subroutine 'forces' was parallelized using VPPfortran
! according with atoms. This will be parallelized using MPI.
! 
  use m_Const_Parameters, only : DP
  ! Revised according to an indication by Momita-san,
  !     T. Yamasaki (FUJITSU Laboratories ltd.), 2003/07/28 
  !      [ip = npfatm(ia)] -> [ ip = npfatm(iaa)]
  implicit none
  integer, intent(in)                            :: natm2,katm,natm,kopr,nopr
  integer, intent(in), dimension(katm,kopr)      :: napt
  real(kind=DP),intent(in),dimension(3,3,kopr)   :: op
  integer, intent(in), dimension(katm)           :: iwei
  real(kind=DP), intent(inout), dimension(natm,3):: f_l
  real(kind=DP), dimension(natm,3)               :: f_wk
  integer,       dimension(natm2)                :: npfatm

  integer        :: ia, iop, iaa, ip
  real(kind=DP)  :: f(3)

  f_wk = f_l
  f_l  = 0.d0
  call substitute_npfatm      ! --> npfatm
!xocl spread do/ind_katm
  do ia = 1, natm
     do iop = 1, nopr
!!$        iaa = napt(ia,iop); ip = npfatm(ia)
        iaa = napt(ia,iop); ip = npfatm(iaa)
        f = f_wk(ip,1:3)
        if(iaa > natm) f = -f
        f_l(ia,1:3) = f_l(ia,1:3) + matmul(transpose(op(1:3,1:3,iop)),f)
     end do
  end do
!xocl end spread
  f_l = f_l / nopr

contains
  subroutine substitute_npfatm
    integer :: nb, ia
    nb = 0
    do ia = 1, natm
       if(iwei(ia) == 2) then
          nb = nb + 1
          npfatm(natm+nb) = ia
       end if
       npfatm(ia) = ia
    end do
  end subroutine substitute_npfatm
end subroutine fd_symmetrize

!!$subroutine fd_symmetrize(natm2,katm,natm,napt,kopr,op,kopr2,nopr &
!!$     & ,iwei, f_l, f_wk, npfatm, iop_supercell)
!!$! This subroutine is for symmetrization of atomic forces.
!!$! Original subroutine 'forces' was parallelized using VPPfortran
!!$! according with atoms. This will be parallelized using MPI.
!!$! 
!!$  use m_Const_Parameters, only : DP
!!$  ! Revised according to an indication by Momita-san,
!!$  !     T. Yamasaki (FUJITSU Laboratories ltd.), 2003/07/28 
!!$  !      [ip = npfatm(ia)] -> [ ip = npfatm(iaa)]
!!$  implicit none
!!$  integer, intent(in)                            :: natm2,katm,natm,kopr,nopr,kopr2
!!$  integer, intent(in), dimension(katm,kopr)      :: napt
!!$  real(kind=DP),intent(in),dimension(3,3,kopr2)  :: op
!!$  integer, intent(in), dimension(katm)           :: iwei
!!$  real(kind=DP), intent(inout), dimension(natm,3):: f_l
!!$  real(kind=DP), dimension(natm,3)               :: f_wk
!!$  integer,       dimension(natm2)                :: npfatm
!!$  integer, intent(in), optional, dimension(nopr) :: iop_supercell
!!$
!!$  integer        :: ia, iop, iaa, ip, iop_temp
!!$  real(kind=DP)  :: f(3)
!!$
!!$  f_wk = f_l
!!$  f_l  = 0.d0
!!$  call substitute_npfatm      ! --> npfatm
!!$!!$  write(6,*) ' -- ia, iop, iaa, ip --'
!!$  if(present(iop_supercell)) then
!!$!xocl spread do/ind_katm
!!$     do ia = 1, natm
!!$        do iop = 1, nopr
!!$           iaa = napt(ia,iop); ip = npfatm(iaa)
!!$           f = f_wk(ip,1:3)
!!$           if(iaa > natm) f = -f
!!$           iop_temp = iop_supercell(iop)
!!$           f_l(ia,1:3) = f_l(ia,1:3) + matmul(transpose(op(1:3,1:3,iop_temp)),f)
!!$        end do
!!$     end do
!!$!xocl end spread
!!$  else
!!$!xocl spread do/ind_katm
!!$     do ia = 1, natm
!!$        do iop = 1, nopr
!!$!!$        iaa = napt(ia,iop); ip = npfatm(ia)
!!$           iaa = napt(ia,iop); ip = npfatm(iaa)
!!$           f = f_wk(ip,1:3)
!!$           if(iaa > natm) f = -f
!!$           f_l(ia,1:3) = f_l(ia,1:3) + matmul(transpose(op(1:3,1:3,iop)),f)
!!$        end do
!!$     end do
!!$!xocl end spread
!!$  end if
!!$  f_l = f_l / nopr
!!$
!!$contains
!!$  subroutine substitute_npfatm
!!$    integer :: nb, ia
!!$    nb = 0
!!$    do ia = 1, natm
!!$       if(iwei(ia) == 2) then
!!$          nb = nb + 1
!!$          npfatm(natm+nb) = ia
!!$       end if
!!$       npfatm(ia) = ia
!!$    end do
!!$  end subroutine substitute_npfatm
!!$end subroutine fd_symmetrize

subroutine gnrt_supercell_symm_operations(natm2,katm,natm,napt,kopr,nopr,op &
     &  , koprpaf,tau,ngen_tl,tau_tl,napt_tl,lattice_system,mode &
     &  , dim2,  napt_supercell, iop_supercell, tau_supercell, nopr_supercell &
     &  , dim2e, nope_supercell, pope_supercell,iwei)
! ==== coded by T. Yamasaki after a provided code by Usami-san on Nov. 2013. April 2014 ===
! 
  use m_Const_Parameters,  only : DP, PAI2, DELTA10,DELTA07
  implicit none
  integer,parameter :: modeTAU = 1, modeAPT = 2
  integer, intent(in)                            :: natm2,katm,natm,kopr,nopr
  integer, intent(in), dimension(katm,kopr)      :: napt
  real(kind=DP),intent(in),dimension(3,3,kopr)   :: op
  integer, intent(in)                            :: koprpaf
  real(kind=DP),intent(in),dimension(3,koprpaf,2):: tau
  integer, intent(in)                            :: ngen_tl
  real(kind=DP),intent(in),dimension(3,ngen_tl,2):: tau_tl
  integer, intent(in), dimension(katm,ngen_tl)   :: napt_tl
  character(len=9), intent(in)                   :: lattice_system
  integer, intent(in)                            :: mode
  integer, intent(in)                            :: dim2
  integer, intent(out),dimension(katm,dim2)      :: napt_supercell
  integer, intent(out),dimension(dim2)           :: iop_supercell
  real(kind=DP),intent(out),dimension(3,dim2)    :: tau_supercell
  integer, intent(out)                           :: nopr_supercell
  integer, intent(in)                            :: dim2e
  integer, intent(out),dimension(kopr)           :: nope_supercell
  integer, intent(out),dimension(dim2e,kopr)     :: pope_supercell
  integer, intent(in), dimension(katm)           :: iwei
  
  integer        :: ia, iop, jop, ip

  integer,             allocatable, dimension(:,:)   :: iop_opxop   ! d(nopr:nopr)
  real(kind=DP), dimension(3,3) :: op_temp
  integer                       :: iop_temp, napt1_temp
  real(kind=DP), dimension(3)   :: tau_temp
  real(kind=DP) ::  delta
  integer :: i, j, nopr_local_temp
  logical :: isComplete

  integer,       dimension(natm2)                :: npfatm
  integer :: icount

  allocate(iop_opxop(nopr,nopr));
  do iop = 1, nopr
     do jop = 1, nopr
        op_temp = matmul(op(:,:,jop), op(:,:,iop))
        iop_opxop(iop,jop) = operatorExists_in_nop()
     end do
  end do

  nope_supercell = 0
  napt_supercell = 0

  ! COPY
  do iop = 1, nopr
     iop_supercell(iop) = iop
!-------------------
     napt_supercell(:,iop) = napt(:,iop)
     nope_supercell(iop) = 1
     pope_supercell(1,iop) = iop
     if(mode==modeTAU) tau_supercell(:, iop) = tau(:, iop, 2)
  end do

  if(ngen_tl >= 1) then
     OP_LOOP: do iop = 1, nopr
        delta = dabs(op(1,1,iop)-1.0d0) + dabs(op(2,1,iop)) +       dabs(op(3,1,iop)) &
             &+ dabs(op(1,2,iop)) +       dabs(op(2,2,iop)-1.0d0) + dabs(op(3,2,iop)) &
             &+ dabs(op(1,3,iop))       + dabs(op(2,3,iop)) +       dabs(op(3,3,iop)-1.0d0)
        if(delta < DELTA10) exit OP_LOOP
     end do OP_LOOP
     if(iop <= nopr) then
!!$        write(6,'(" # of A Unit operator = ",i8)') iop
     else
        !stop ' no Unit operator <gnrt_supercell_symm_operations>'
        call phase_error_with_msg(6, 'no Unit operator <gnrt_supercell_symm_operations>',__LINE__,__FILE__)
     end if
     do i = 1, ngen_tl
        iop_supercell(nopr+i) = iop
        napt_supercell(:,nopr+i) = napt_tl(:,i)
        nope_supercell(iop) = nope_supercell(iop)+1
        pope_supercell(1+i,iop) = nopr+i
        if(mode==modeTAU) tau_supercell(:,nopr+i) = tau_tl(:,i,2)
     end do
  end if
     
  nopr_supercell = nopr + ngen_tl
  call substitute_npfatm()
  icount = 0
  EXTEND: do 
     if (nopr_supercell > natm*nopr) then
        write(6,*) 'WARNING: Force symmetrization might be a failure.'
        nopr_supercell = nopr
        exit EXTEND
     end if
     icount = icount+1
!     if(mype==0) write(0,*) 'nopr_sup, icount ',nopr_supercell,icount
     ! MAKE TABLES (REDUNDANT SPACE-GROUP)
     isComplete = .true.
     nopr_local_temp = nopr_supercell

     do iop = 1, nopr_local_temp
!        if(mype==0) write(0,*) 'iop, nopr_local_temp',iop,nopr_local_temp
        do jop = 1, nopr_local_temp
           iop_temp = iop_opxop(iop_supercell(iop),iop_supercell(jop))

           if(mode==modeTAU) then
              tau_temp = matmul(op(1:3,1:3,iop_supercell(jop)),tau_supercell(1:3,iop))+tau_supercell(1:3,jop)
              if (lattice_system .eq. 'hexagonal') then
                 ! y-axis, hexagonal only
                 tau_temp(2) = mod(tau_temp(2)+sqrt(3.0d0), sqrt(0.75d0))
              else
                 tau_temp = tau_temp - floor(tau_temp)
              end if
              if(tau_temp(1) > 1.0d0 - DELTA07) tau_temp(1) = 0.d0
              if(tau_temp(2) > 1.0d0 - DELTA07) tau_temp(2) = 0.d0
              if(tau_temp(3) > 1.0d0 - DELTA07) tau_temp(3) = 0.d0
           else if(mode== modeAPT)then
              if(napt_supercell(1,iop)>0) &
              &  napt1_temp = npfatm(napt_supercell(npfatm(napt_supercell(1,iop)),jop))
           end if
           
           if (.not. operatorExists()) then
              nopr_supercell = nopr_supercell + 1
              nope_supercell(iop_temp) = nope_supercell(iop_temp) + 1
              pope_supercell(nope_supercell(iop_temp),iop_temp) = nopr_supercell
              if (nopr_supercell > dim2) cycle EXTEND
              iop_supercell(nopr_supercell) = iop_temp
              if(mode==modeTAU) then
                 tau_supercell(1:3,     nopr_supercell) = tau_temp(1:3)
              else if(mode==modeAPT) then
                 do i = 1, natm
                    if(napt_supercell(i,iop)>0)            &
                    &  napt_supercell(i, nopr_supercell) = &
                    &  npfatm(napt_supercell(npfatm(napt_supercell(i, iop)), jop))
                 end do
              end if
              isComplete = .false.
           end if
        end do
     end do
     if (isComplete) exit EXTEND
  end do EXTEND

  deallocate(iop_opxop)
!!$  mnope_supercell = maxval(nope_supercell(1:nopr))

contains
  subroutine substitute_npfatm()
    integer :: nb, ia
    nb = 0
    do ia = 1, natm
       if(iwei(ia) == 2) then
          nb = nb + 1
          npfatm(natm+nb) = ia
       end if
       npfatm(ia) = ia
    end do
  end subroutine substitute_npfatm

  integer function operatorExists_in_nop()
    real(kind=DP) :: delta
    integer :: nope
    real(kind=DP), dimension(3,3) :: opd

    do nope = 1, nopr
       opd = op_temp - op(:,:,nope)
       delta = opd(1,1)*opd(1,1) + opd(1,2)*opd(1,2) + opd(1,3)*opd(1,3) &
            &+ opd(2,1)*opd(2,1) + opd(2,2)*opd(2,2) + opd(2,3)*opd(2,3) &
            &+ opd(3,1)*opd(3,1) + opd(3,2)*opd(3,2) + opd(3,3)*opd(3,3)
       if(delta < 1.0d-3) then
          operatorExists_in_nop = nope
          return
       end if
    end do
    operatorExists_in_nop = -1
    if(operatorExists_in_nop == -1) then
       write(6,'(" iop, jop = ",2i8)') iop,jop
       write(6,'(" !Illegal value of iop_opxop(iop,jop) = ",i8)') iop_opxop(iop,jop)
       call flush(6)
       stop '!Illegal value of iop_opxop(iop,jop)'
    end if
  end function operatorExists_in_nop
!-----------------------------------------------------------------------------------

  logical function operatorExists()
    ! return .true. if the operator (op_temp, tau_temp) is found in the operator list op_local.
    real(kind=DP) :: delta
    real(kind=DP), dimension(3) :: temp_a, temp
    integer :: nope, nopei

    do nopei = 1, nope_supercell(iop_temp)
       nope = pope_supercell(nopei,iop_temp)
       if(mode == modeTAU) then
          temp_a(:) = tau_supercell(:,nope) - tau_temp(:)
          temp(1:3) = temp_a(1:3) - floor(temp_a(1:3))
          if(temp(1) > 1.0d0 - DELTA07) temp(1) = 0.d0
          if(temp(2) > 1.0d0 - DELTA07) temp(2) = 0.d0
          if(temp(3) > 1.0d0 - DELTA07) temp(3) = 0.d0
          delta = temp(1)*temp(1) + temp(2)*temp(2) + temp(3)*temp(3)

          if (delta < 1.0d-3) then
             operatorExists = .true.  ! found
             return
          end if
       else if(mode == modeAPT) then
          if(npfatm(napt_supercell(1,nope)) == napt1_temp) then
             operatorExists = .true. ! found
             return
          end if
       end if
    end do
    operatorExists = .false.  ! not found
    return
  end function operatorExists
!-----------------------------------------------------------------------------------
end subroutine gnrt_supercell_symm_operations

subroutine gnrt_supercell_symm_operations_apt(natm2,katm,natm,napt,kopr,nopr,op &
     &  , koprpaf,tau,ngen_tl,tau_tl,napt_tl,lattice_system,mode &
     &  , dim2,  napt_supercell, iop_supercell,  nopr_supercell &
     &  , dim2e, iwei)
! ==== coded by T. Yamasaki after a provided code by Usami-san on Nov. 2013. April 2014 ===
! 
  use m_Const_Parameters,  only : DP, PAI2, DELTA10,DELTA07
  use m_Parallelization,   only : mype
  implicit none
  integer,parameter :: modeTAU = 1, modeAPT = 2
  integer, intent(in)                            :: natm2,katm,natm,kopr,nopr
  integer, intent(in), dimension(katm,kopr)      :: napt
  real(kind=DP),intent(in),dimension(3,3,kopr)   :: op
  integer, intent(in)                            :: koprpaf
  real(kind=DP),intent(in),dimension(3,koprpaf,2):: tau
  integer, intent(in)                            :: ngen_tl
  real(kind=DP),intent(in),dimension(3,ngen_tl,2):: tau_tl
  integer, intent(in), dimension(katm,ngen_tl)   :: napt_tl
  character(len=9), intent(in)                   :: lattice_system
  integer, intent(in)                            :: mode
  integer, intent(in)                            :: dim2
  integer, intent(out),dimension(katm,dim2)      :: napt_supercell
  integer, intent(out),dimension(dim2)           :: iop_supercell
  integer, intent(out)                           :: nopr_supercell
  integer, intent(in)                            :: dim2e
!  integer, intent(out),dimension(kopr)           :: nope_supercell
!  integer, intent(out),dimension(dim2e,kopr)     :: pope_supercell
  integer, intent(in), dimension(katm)           :: iwei

  integer, allocatable, dimension(:)             :: nope_supercell
  integer, allocatable, dimension(:,:)           :: pope_supercell
  
  integer        :: ia, iop, jop, ip

  integer,             allocatable, dimension(:,:)   :: iop_opxop   ! d(nopr:nopr)
  real(kind=DP), dimension(3,3) :: op_temp
  integer                       :: iop_temp, napt1_temp
  real(kind=DP), dimension(3)   :: tau_temp
  real(kind=DP) ::  delta
  integer :: i, j, nopr_local_temp
  logical :: isComplete

  integer,       dimension(natm2)                :: npfatm
  integer :: icount
  allocate(nope_supercell(kopr))
  allocate(pope_supercell(dim2e,kopr))
  allocate(iop_opxop(nopr,nopr));
  do iop = 1, nopr
     do jop = 1, nopr
        op_temp = matmul(op(:,:,jop), op(:,:,iop))
        iop_opxop(iop,jop) = operatorExists_in_nop()
     end do
  end do

  nope_supercell = 0
  napt_supercell = 0

  ! COPY
  do iop = 1, nopr
     iop_supercell(iop) = iop
!-------------------
     napt_supercell(:,iop) = napt(:,iop)
     nope_supercell(iop) = 1
     pope_supercell(1,iop) = iop
  end do

  if(ngen_tl >= 1) then
     OP_LOOP: do iop = 1, nopr
        delta = dabs(op(1,1,iop)-1.0d0) + dabs(op(2,1,iop)) +       dabs(op(3,1,iop)) &
             &+ dabs(op(1,2,iop)) +       dabs(op(2,2,iop)-1.0d0) + dabs(op(3,2,iop)) &
             &+ dabs(op(1,3,iop))       + dabs(op(2,3,iop)) +       dabs(op(3,3,iop)-1.0d0)
        if(delta < DELTA10) exit OP_LOOP
     end do OP_LOOP
     if(iop <= nopr) then
!!$        write(6,'(" # of A Unit operator = ",i8)') iop
     else
        stop ' no Unit operator <gnrt_supercell_symm_operations>'
     end if
     do i = 1, ngen_tl
        iop_supercell(nopr+i) = iop
        napt_supercell(:,nopr+i) = napt_tl(:,i)
        nope_supercell(iop) = nope_supercell(iop)+1
        pope_supercell(1+i,iop) = nopr+i
     end do
  end if
     
  nopr_supercell = nopr + ngen_tl
  call substitute_npfatm()
  icount = 0
  EXTEND: do 
     if (nopr_supercell > natm*nopr) then
        write(6,*) 'WARNING: Force symmetrization might be a failure.'
        nopr_supercell = nopr
        exit EXTEND
     end if
     icount = icount+1
     isComplete = .true.
     nopr_local_temp = nopr_supercell

     do iop = 1, nopr_local_temp
        !if(mype==0) write(0,*) 'iop, nopr_local_temp',iop,nopr_local_temp
        do jop = 1, nopr_local_temp
           iop_temp = iop_opxop(iop_supercell(iop),iop_supercell(jop))
           if(napt_supercell(1,iop)>0) &
           &  napt1_temp = npfatm(napt_supercell(npfatm(napt_supercell(1,iop)),jop))
           
           if (.not. operatorExists(kopr,iop_temp,napt1_temp,dim2e, &
                     natm2,npfatm,nope_supercell,pope_supercell)) then
              nopr_supercell = nopr_supercell + 1
              nope_supercell(iop_temp) = nope_supercell(iop_temp) + 1
              pope_supercell(nope_supercell(iop_temp),iop_temp) = nopr_supercell
              if (nopr_supercell > dim2) cycle EXTEND
              iop_supercell(nopr_supercell) = iop_temp
              do i = 1, natm
                 if(napt_supercell(i,iop)>0)            &
                 &  napt_supercell(i, nopr_supercell) = &
                 &  npfatm(napt_supercell(npfatm(napt_supercell(i, iop)), jop))
              end do
              isComplete = .false.
           end if
        end do
     end do
     if (isComplete) exit EXTEND
  end do EXTEND

  deallocate(iop_opxop)
  deallocate(nope_supercell)
  deallocate(pope_supercell)
!!$  mnope_supercell = maxval(nope_supercell(1:nopr))

contains
  subroutine substitute_npfatm()
    integer :: nb, ia
    nb = 0
    do ia = 1, natm
       if(iwei(ia) == 2) then
          nb = nb + 1
          npfatm(natm+nb) = ia
       end if
       npfatm(ia) = ia
    end do
  end subroutine substitute_npfatm

  integer function operatorExists_in_nop()
    real(kind=DP) :: delta
    integer :: nope
    real(kind=DP), dimension(3,3) :: opd

    do nope = 1, nopr
       opd = op_temp - op(:,:,nope)
       delta = opd(1,1)*opd(1,1) + opd(1,2)*opd(1,2) + opd(1,3)*opd(1,3) &
            &+ opd(2,1)*opd(2,1) + opd(2,2)*opd(2,2) + opd(2,3)*opd(2,3) &
            &+ opd(3,1)*opd(3,1) + opd(3,2)*opd(3,2) + opd(3,3)*opd(3,3)
       if(delta < 1.0d-3) then
          operatorExists_in_nop = nope
          return
       end if
    end do
    operatorExists_in_nop = -1
    if(operatorExists_in_nop == -1) then
       write(6,'(" iop, jop = ",2i8)') iop,jop
       write(6,'(" !Illegal value of iop_opxop(iop,jop) = ",i8)') iop_opxop(iop,jop)
       call flush(6)
       stop '!Illegal value of iop_opxop(iop,jop)'
    end if
  end function operatorExists_in_nop
!-----------------------------------------------------------------------------------

  logical function operatorExists(kopr,iop_temp,napt1_temp,dim2e,natm2,npfatm,nope_supercell,pope_supercell)
    ! return .true. if the operator (op_temp, tau_temp) is found in the operator list op_local.
    integer, intent(in):: kopr,iop_temp,napt1_temp,dim2e,natm2
    integer, intent(in), dimension(natm2) :: npfatm
    integer, intent(in), dimension(kopr)  :: nope_supercell
    integer, intent(in), dimension(dim2e,kopr) :: pope_supercell
    real(kind=DP) :: delta
    real(kind=DP), dimension(3) :: temp_a, temp
    integer :: nope, nopei

    do nopei = 1, nope_supercell(iop_temp)
       nope = pope_supercell(nopei,iop_temp)
       if(npfatm(napt_supercell(1,nope)) == napt1_temp) then
          operatorExists = .true. ! found
          return
       end if
    end do
    operatorExists = .false.  ! not found
    return
  end function operatorExists
!-----------------------------------------------------------------------------------
end subroutine gnrt_supercell_symm_operations_apt


subroutine fd_supercell_symmetrize(natm2,katm,natm,napt,nopr,op,kopr,iwei &
     & , f_l, f_wk, npfatm, iop_supercell)
!
! ===  Changed from the subroutine fd_symmetrize by T. Yamasaki after a provided code ===
! === by Usami-san on Nov. 2013.                                        April 2014    ===
!
! This subroutine is for symmetrization of atomic forces.
! Original subroutine 'forces' was parallelized using VPPfortran
! according with atoms. This will be parallelized using MPI.
! 
  use m_Const_Parameters,  only : DP, PAI2, DELTA10,DELTA07

  implicit none
  integer, intent(in)                            :: natm2,katm,natm,kopr,nopr
  integer, intent(in), dimension(katm,nopr) :: napt
  real(kind=DP),intent(in),dimension(3,3,kopr)   :: op
  integer, intent(in), dimension(katm)           :: iwei
  real(kind=DP), intent(inout), dimension(natm,3):: f_l
  real(kind=DP), dimension(natm,3)               :: f_wk
  integer,       dimension(natm2)                :: npfatm
  integer, intent(in), dimension(nopr) :: iop_supercell

  integer        :: ia, iop, iaa, ip, iop_temp
  real(kind=DP)  :: f(3)

  f_wk = f_l
  f_l  = 0.d0
  call substitute_npfatm      ! --> npfatm
!xocl spread do/ind_katm
  do ia = 1, natm
     do iop = 1, nopr
        iaa = napt(ia,iop); ip = npfatm(iaa)
        f = f_wk(ip,1:3)
        if(iaa > natm) f = -f
        iop_temp = iop_supercell(iop)
        f_l(ia,1:3) = f_l(ia,1:3) + matmul(transpose(op(1:3,1:3,iop_temp)),f)
     end do
  end do
!xocl end spread
  f_l = f_l / nopr
contains
  subroutine substitute_npfatm
    integer :: nb, ia
    nb = 0
    do ia = 1, natm
       if(iwei(ia) == 2) then
          nb = nb + 1
          npfatm(natm+nb) = ia
       end if
       npfatm(ia) = ia
    end do
  end subroutine substitute_npfatm
end subroutine fd_supercell_symmetrize

!!$subroutine fd_symmetrize_sc(natm2,katm,natm,napt,kopr,nopr,op,iwei &
!!$     &, f_l, f_wk, npfatm)
!!$! This subroutine is for symmetrization of atomic forces.
!!$! Original subroutine 'forces' was parallelized using VPPfortran
!!$! according with atoms. This will be parallelized using MPI.
!!$! 
!!$  use m_Const_Parameters,  only : DP, PAI2, DELTA10,DELTA07
!!$  ! Revised according to an indication by Momita-san,
!!$  !     T. Yamasaki (FUJITSU Laboratories ltd.), 2003/07/28 
!!$  !      [ip = npfatm(ia)] -> [ ip = npfatm(iaa)]
!!$  implicit none
!!$  integer, intent(in)                            :: natm2,katm,natm,kopr,nopr
!!$  integer, intent(in), dimension(katm,kopr)      :: napt
!!$  real(kind=DP),intent(in),dimension(3,3,kopr)   :: op
!!$  integer, intent(in), dimension(katm)           :: iwei
!!$  real(kind=DP), intent(inout), dimension(natm,3):: f_l
!!$  real(kind=DP), dimension(natm,3)               :: f_wk
!!$  integer,       dimension(natm2)                :: npfatm
!!$
!!$  integer        :: ia, iop, jop, iaa, ip
!!$  real(kind=DP)  :: f(3)
!!$
!!$  integer,       save, allocatable, dimension(:,:)   :: napt_local  ! d(natm,nopr_local)
!!$  integer,       save, allocatable, dimension(:)     :: iop_local   ! d(nopr_local)
!!$  integer,       save                                :: nopr_local
!!$  integer,             allocatable, dimension(:,:)   :: pope_local  ! d(:,nopr)
!!$  integer,             allocatable, dimension(:)     :: nope_local  ! d(nopr)
!!$  integer,             allocatable, dimension(:,:)   :: iop_opxop   ! d(nopr:nopr)
!!$  !
!!$  real(kind=DP), dimension(3,3) :: op_temp
!!$  integer                       :: iop_temp, napt1_temp
!!$  integer :: i, j, nopr_local_temp, dim_redundancy
!!$  logical :: isComplete, isInitialized = .false.
!!$  logical :: DEBUG_WRITE = .true.
!!$
!!$  if(.not.isInitialized) then
!!$     allocate(iop_opxop(nopr,nopr));
!!$     do iop = 1, nopr
!!$        do jop = 1, nopr
!!$           op_temp = matmul(op(:,:,jop), op(:,:,iop))
!!$           iop_opxop(iop,jop) = operatorExists_in_nop()
!!$        end do
!!$     end do
!!$  end if
!!$
!!$  if(.not.isInitialized) then
!!$     dim_redundancy = natm*nopr
!!$     allocate( napt_local(natm, dim_redundancy))
!!$     allocate( iop_local(dim_redundancy), pope_local(dim_redundancy/nopr,nopr))
!!$     allocate( nope_local(nopr)); nope_local = 0
!!$
!!$     ! COPY
!!$     do iop = 1, nopr
!!$        iop_local(iop) = iop
!!$!-------------------
!!$        napt_local(:,iop) = napt(:,iop)
!!$        nope_local(iop) = 1
!!$        pope_local(1,iop) = iop
!!$     end do
!!$     nopr_local = nopr
!!$
!!$     EXTEND: do 
!!$        if (dim_redundancy > natm*nopr) then
!!$           write(6,*) 'WARNING: Force symmetrization might be a failure.'
!!$           nopr_local = nopr
!!$           exit EXTEND
!!$        end if
!!$
!!$     ! MAKE TABLES (REDUNDANT SPACE-GROUP)
!!$        isComplete = .true.
!!$        nopr_local_temp = nopr_local
!!$        do iop = 1, nopr_local_temp
!!$           do jop = 1, nopr_local_temp
!!$              iop_temp = iop_opxop(iop_local(iop),iop_local(jop))
!!$              napt1_temp = napt_local(napt_local(1,iop),jop)
!!$              if (.not. operatorExists_nap()) then
!!$                 nopr_local = nopr_local + 1
!!$                 nope_local(iop_temp) = nope_local(iop_temp) + 1
!!$                 pope_local(nope_local(iop_temp),iop_temp) = nopr_local
!!$                 if (nopr_local > dim_redundancy) cycle EXTEND
!!$                 iop_local(nopr_local) = iop_temp
!!$                 do i = 1, natm
!!$                    napt_local(i, nopr_local) = napt_local( napt_local(i, iop), jop)
!!$                 end do
!!$                 isComplete = .false.
!!$              end if
!!$           end do
!!$        end do
!!$        if (isComplete) exit EXTEND
!!$     end do EXTEND
!!$
!!$     If(dim_redundancy > nopr_local) call reduce_napt_iop_arraysize()
!!$
!!$     if(DEBUG_WRITE) call wd_nopr_local_etc()
!!$
!!$     deallocate(iop_opxop)
!!$     deallocate(nope_local)
!!$     deallocate(pope_local)
!!$  end if
!!$
!!$  isInitialized = .true.
!!$
!!$  f_wk = f_l
!!$  f_l  = 0.d0
!!$  call substitute_npfatm      ! --> npfatm
!!$!xocl spread do/ind_katm
!!$  do ia = 1, natm
!!$     do iop = 1, nopr_local
!!$!!!!$        iaa = napt(ia,iop); ip = npfatm(ia)
!!$        iaa = napt_local(ia,iop); ip = npfatm(iaa)
!!$        f = f_wk(ip,1:3)
!!$        if(iaa > natm) f = -f
!!$        iop_temp = iop_local(iop)
!!$        f_l(ia,1:3) = f_l(ia,1:3) + matmul(transpose(op(1:3,1:3,iop_temp)),f)
!!$     end do
!!$  end do
!!$!xocl end spread
!!$  f_l = f_l / nopr_local
!!$contains
!!$  subroutine wd_nopr_local_etc()
!!$    integer :: iop, i, ic, j, k
!!$!!!!$     if(nopr_local > nopr) then
!!$    write(6,*) ' --- Symmetry Operations (CARTS, PUCV) ---'
!!$    write(6,'(" !! nopr+af = ",i8)') nopr_local
!!$    do iop = 1, nopr
!!$       write(6,'(" #symmetry op. = ",i8)')  iop
!!$       do i = 1, 3
!!$          write(6,'(3f8.4)') (op(i,k,iop),k=1,3)
!!$       end do
!!$       ic = 0
!!$       do j = 1, nopr_local
!!$          if(iop_local(j) == iop) then
!!$             ic = ic+1
!!$             write(6,'(i3,2x,i6,2x,9f5.1)') ic,j,((op(k,i,iop),i=1,3),k=1,3)
!!$          end if
!!$       end do
!!$       write(6,'(" # operation",i3," = ",i8)') iop,ic
!!$    end do
!!$  end subroutine wd_nopr_local_etc
!!$
!!$  subroutine reduce_napt_iop_arraysize()
!!$    integer, allocatable, dimension(:,:) :: napt_temp
!!$    integer, allocatable, dimension(:)   :: iop_temp
!!$    integer :: iatm, iop
!!$    allocate(napt_temp(natm,nopr_local))
!!$    allocate(iop_temp(nopr_local))
!!$
!!$    do iop = 1, nopr_local
!!$       iop_temp(iop) = iop_local(iop)
!!$       do ia = 1, natm
!!$          napt_temp(ia,iop) = napt_local(ia,iop)
!!$       end do
!!$    end do
!!$
!!$    deallocate(napt_local, iop_local)
!!$
!!$    allocate(napt_local(natm,nopr_local))
!!$    allocate(iop_local(nopr_local))
!!$
!!$    iop_local = iop_temp
!!$    napt_local = napt_temp
!!$
!!$    deallocate(napt_temp, iop_temp)
!!$
!!$    write(6,'(" ! napt_local and iop_local have been reduced")')
!!$    write(6,'(" ! dim_redundancy = ", i8)') dim_redundancy
!!$    write(6,'(" ! nopr_local     = ", i8)') nopr_local
!!$  end subroutine reduce_napt_iop_arraysize
!!$
!!$  subroutine substitute_npfatm
!!$    integer :: nb, ia
!!$    nb = 0
!!$    do ia = 1, natm
!!$       if(iwei(ia) == 2) then
!!$          nb = nb + 1
!!$          npfatm(natm+nb) = ia
!!$       end if
!!$       npfatm(ia) = ia
!!$    end do
!!$  end subroutine substitute_npfatm
!!$!-----------------------------------------------------------------------------------
!!$  integer function operatorExists_in_nop()
!!$    real(kind=DP) :: delta
!!$    integer :: nope
!!$    real(kind=DP), dimension(3,3) :: opd
!!$
!!$    do nope = 1, nopr
!!$       opd = op_temp - op(:,:,nope)
!!$       delta = opd(1,1)*opd(1,1) + opd(1,2)*opd(1,2) + opd(1,3)*opd(1,3) &
!!$            &+ opd(2,1)*opd(2,1) + opd(2,2)*opd(2,2) + opd(2,3)*opd(2,3) &
!!$            &+ opd(3,1)*opd(3,1) + opd(3,2)*opd(3,2) + opd(3,3)*opd(3,3)
!!$       if(delta < 1.0d-3) then
!!$          operatorExists_in_nop = nope
!!$          return
!!$       end if
!!$    end do
!!$    operatorExists_in_nop = -1
!!$    if(operatorExists_in_nop == -1) then
!!$       write(6,'(" iop, jop = ",2i8)') iop,jop
!!$       write(6,'(" !Illegal value of iop_opxop(iop,jop) = ",i8)') iop_opxop(iop,jop)
!!$       call flush(6)
!!$       stop '!Illegal value of iop_opxop(iop,jop)'
!!$    end if
!!$  end function operatorExists_in_nop
!!$!-----------------------------------------------------------------------------------
!!$  logical function operatorExists_nap()
!!$    ! return .true. if the operator (op_temp, napt1_temp) is found in the operator list op_local.
!!$    real(kind=DP) :: delta ,temp
!!$    integer       :: nope, nopei
!!$
!!$    do nopei = 1, nope_local(iop_temp)
!!$       nope = pope_local(nopei,iop_temp)
!!$       if(napt_local(1,nope) == napt1_temp) then
!!$          operatorExists_nap = .true. ! found
!!$          return
!!$       end if
!!$    end do
!!$    operatorExists_nap = .false.  ! not found
!!$    return
!!$
!!$  end function operatorExists_nap
!!$!-----------------------------------------------------------------------------------
!!$end subroutine fd_symmetrize_sc
!!$
!!$subroutine fd_symmetrize_sc2(natm2,katm,natm,napt,kopr,nopr,op,iwei &
!!$     &, f_l, f_wk, npfatm)

!!$! Original subroutine 'forces' was parallelized using VPPfortran
!!$! according with atoms. This will be parallelized using MPI.
!!$! 
!!$  use m_Const_Parameters,  only : DP, PAI2
!!$  use m_Crystal_Structure, only : tau, b2pmat, rltv
!!$  use m_CS_SpaceGroup,     only : lattice_system
!!$  ! Revised according to an indication by Momita-san,
!!$  !     T. Yamasaki (FUJITSU Laboratories ltd.), 2003/07/28 
!!$  !      [ip = npfatm(ia)] -> [ ip = npfatm(iaa)]
!!$  implicit none
!!$  integer, intent(in)                            :: natm2,katm,natm,kopr,nopr
!!$  integer, intent(in), dimension(katm,kopr)      :: napt
!!$  real(kind=DP),intent(in),dimension(3,3,kopr)   :: op
!!$  integer, intent(in), dimension(katm)           :: iwei
!!$  real(kind=DP), intent(inout), dimension(natm,3):: f_l
!!$  real(kind=DP), dimension(natm,3)               :: f_wk
!!$  integer,       dimension(natm2)                :: npfatm
!!$
!!$  integer        :: ia, iop, jop, iaa, ip
!!$  real(kind=DP)  :: f(3)
!!$
!!$  real(kind=DP), save, allocatable, dimension(:,:,:) :: op_local
!!$  real(kind=DP), save, allocatable, dimension(:,:)   :: tau_local  ! (=ta1) d(3,nopr), nonprimitive translation vector (A system)
!!$  integer,       save, allocatable, dimension(:,:)   :: napt_local
!!$  integer,       save, allocatable, dimension(:)     :: iop_local
!!$  integer,       save                                :: nopr_local
!!$  !
!!$  real(kind=DP), dimension(3,3) :: op_temp, matrix_temp
!!$  integer                       :: iop_temp
!!$  real(kind=DP), dimension(3)   :: tau_temp
!!$  integer :: i, nopr_local_temp, dim_redundancy, j,k, ic
!!$  logical :: isComplete, isInitialized = .false.
!!$
!!$!  call fd_symmetrize(natm2,katm,natm,napt,kopr,nopr,op,iwei &
!!$!     &, f_l, f_wk, npfatm)
!!$!  return
!!$
!!$  dim_redundancy = -52  ! -52 + 100 = 48
!!$
!!$  EXTEND: do 
!!$     if (isInitialized) exit EXTEND
!!$
!!$     if (dim_redundancy > 32000) then
!!$        write(6,*) 'WARNING: Force symmetrization might be a failure.'
!!$        nopr_local = nopr
!!$        exit EXTEND
!!$     else if (dim_redundancy > 0) then
!!$        deallocate(op_local, tau_local, napt_local,iop_local)  ! not the first-time
!!$     end if
!!$     if (dim_redundancy < 600) then
!!$        dim_redundancy = dim_redundancy + 100
!!$     else
!!$        dim_redundancy = dim_redundancy * 2
!!$     end if
!!$     write(6,'(" dim_redundancy = ",i10)') dim_redundancy
!!$     call flush(6)
!!$     allocate( op_local(3, 3, dim_redundancy), tau_local(3, dim_redundancy), napt_local(natm, dim_redundancy),iop_local(dim_redundancy) )
!!$
!!$     matrix_temp(:,:) = matmul(transpose(b2pmat(:,:)), rltv(:,:)) / PAI2
!!$     ! COPY
!!$     do iop = 1, nopr
!!$!  --- by T. Y. 2014/04/04
!!$        op_local(1:3, 1:3, iop) = transpose(op(1:3, 1:3, iop))
!!$!!        op_local(1:3, 1:3, iop) = op(1:3, 1:3, iop)
!!$!  ---
!!$        iop_local(iop) = iop
!!$!-------------------
!!$        tau_local(:, iop) = matmul(transpose(matrix_temp(:,:)), tau(:, iop,1))
!!$        if (lattice_system .eq. 'hexagonal') then
!!$           tau_local(1,iop) = tau_local(1,iop) - 0.5d0 * tau_local(2,iop)
!!$           tau_local(1,iop) = mod(tau_local(1,iop)+1.0d0, 1.0d0)
!!$           tau_local(2,iop) = sqrt(0.75d0) * tau_local(2,iop)
!!$           tau_local(2,iop) = mod(tau_local(2,iop)+sqrt(3.0d0), sqrt(0.75d0))
!!$        end if
!!$!-------------------
!!$        do i = 1, natm
!!$           napt_local(i, iop) = napt(i, iop)
!!$        end do
!!$     end do
!!$     nopr_local = nopr
!!$
!!$     ! MAKE TABLES (REDUNDANT SPACE-GROUP)
!!$     j = 0
!!$     do
!!$        j = j+1
!!$        isComplete = .true.
!!$        nopr_local_temp = nopr_local
!!$        write(6,'(" j = ",i8)') j
!!$        call flush(6)
!!$        if(j >=10) then
!!$           if(j >= 3200000) then
!!$              stop
!!$           end if
!!$        end if
!!$        write(6,'(" nopr_local_temp = ",i8)') nopr_local_temp
!!$        call flush(6)
!!$        do iop = 1, nopr_local_temp
!!$           do jop = 1, nopr_local_temp
!!$              op_temp = matmul(transpose(op_local(:,:,jop)), transpose(op_local(:,:,iop)))
!!$              op_temp = transpose(op_temp)
!!$              iop_temp = operatorExists_in_nop()
!!$              tau_temp = matmul(transpose(op_local(1:3, 1:3, jop)), tau_local(1:3, iop)) + tau_local(1:3, jop)
!!$              do i = 1, 3
!!$                 ! normalize
!!$                 if (lattice_system .eq. 'hexagonal' .and. i .eq. 2) then
!!$                    ! y-axis, hexagonal only
!!$                    tau_temp(i) = mod(tau_temp(i)+sqrt(3.0d0), sqrt(0.75d0))
!!$                 else
!!$                    tau_temp(i) = mod(tau_temp(i)+2.0d0, 1.0d0)
!!$                 end if
!!$                 if (tau_temp(i) > 0.9999999d0) tau_temp(i) = 0.0d0
!!$              end do
!!$              if (.not. operatorExists()) then
!!$                 nopr_local = nopr_local + 1
!!$                 if (nopr_local > dim_redundancy) then
!!$                    write(6,'(" nopr_local > dim_redundancy : nopr_local = ",i8, " dim_redundancy = ",i8)') nopr_local, dim_redundancy
!!$                    call flush(6)
!!$                    cycle EXTEND
!!$                 end if
!!$                 op_local(1:3, 1:3, nopr_local) = op_temp(1:3, 1:3)
!!$                 iop_local(nopr_local) = iop_temp
!!$                 tau_local(1:3,     nopr_local) = tau_temp(1:3)
!!$                 do i = 1, natm
!!$                    napt_local(i, nopr_local) = napt_local( napt_local(i, iop), jop)
!!$                 end do
!!$                 isComplete = .false.
!!$              end if
!!$           end do
!!$        end do
!!$        if (isComplete) exit EXTEND
!!$     end do
!!$  end do EXTEND
!!$  isInitialized = .true.
!!$
!!$  if(nopr_local > nopr) then
!!$     write(6,*) ' --- Symmetry Operations (CARTS, PUCV) ---'
!!$     write(6,'(" !! nopr+af = ",i8)') nopr_local
!!$     do iop = 1, nopr
!!$        write(6,'(" #symmetry op. = ",i8)')  iop
!!$        do i = 1, 3
!!$           write(6,'(3f8.4)') (op(i,k,iop),k=1,3)
!!$        end do
!!$        ic = 0
!!$        do j = 1, nopr_local
!!$           if(iop_local(j) == iop) then
!!$              ic = ic+1
!!$              write(6,'(i3,2x,i5,2x,3f16.12,2x,9f6.2)') ic,j,tau_local(1:3,j),((op_local(k,i,j),i=1,3),k=1,3)
!!$           end if
!!$        end do
!!$        write(6,'(" # operation",i3," = ",i8)') iop,ic
!!$     end do
!!$  end if
!!$
!!$  write(6,'(" nopr_local = ",i8," (fd_symmetrize_sc)")') nopr_local
!!$  call flush(6)
!!$  f_wk = f_l
!!$  f_l  = 0.d0
!!$  call substitute_npfatm      ! --> npfatm
!!$!xocl spread do/ind_katm
!!$  do ia = 1, natm
!!$     do iop = 1, nopr_local
!!$!!!$        iaa = napt(ia,iop); ip = npfatm(ia)
!!$        iaa = napt_local(ia,iop); ip = npfatm(iaa)
!!$        f = f_wk(ip,1:3)
!!$        if(iaa > natm) f = -f
!!$!!!$        f_l(ia,1:3) = f_l(ia,1:3) + matmul(transpose(op_local(1:3,1:3,iop)),f)
!!$        f_l(ia,1:3) = f_l(ia,1:3) + matmul(op_local(1:3,1:3,iop) ,f)
!!$     end do
!!$  end do
!!$!xocl end spread
!!$  f_l = f_l / nopr_local
!!$contains
!!$  subroutine substitute_npfatm
!!$    integer :: nb, ia
!!$    nb = 0
!!$    do ia = 1, natm
!!$       if(iwei(ia) == 2) then
!!$          nb = nb + 1
!!$          npfatm(natm+nb) = ia
!!$       end if
!!$       npfatm(ia) = ia
!!$    end do
!!$  end subroutine substitute_npfatm
!!$!-----------------------------------------------------------------------------------
!!$  integer function operatorExists_in_nop()
!!$    real(kind=DP) :: delta,temp
!!$    integer :: i,j,nope
!!$    do nope = 1, nopr
!!$       delta = 0.0d0
!!$       do i = 1, 3
!!$          do j = 1, 3
!!$!!             temp = op_temp(i,j) - op(i,j,nope)
!!$             temp = op_temp(i,j) - op(j,i,nope)
!!$             delta = delta + temp*temp
!!$          end do
!!$       end do
!!$       if(delta < 1.0d-3) then
!!$          operatorExists_in_nop = nope
!!$          return
!!$       end if
!!$    end do
!!$    operatorExists_in_nop = -1
!!$  end function operatorExists_in_nop
!!$
!!$  logical function operatorExists()
!!$    ! return .true. if the operator (op_temp, tau_temp) is found in the operator list op_local.
!!$    real(kind=DP) delta ,temp, temp_a, temp_b
!!$    integer i, j, nope
!!$
!!$    do nope = 1, nopr_local
!!$       ! ROTATION
!!$       if(1<=iop_temp.and. iop_temp <= nopr .and. iop_temp /= iop_local(nope)) cycle
!!$       if(iop_temp <= 0)then
!!$          delta = 0.0d0
!!$          do i = 1, 3
!!$             do j = 1, 3
!!$                temp = op_temp(i, j) - op_local(i, j, nope)
!!$                delta = delta + temp * temp
!!$             end do
!!$          end do
!!$          if (delta > 1.0d-3) cycle  ! next operator
!!$       end if
!!$
!!$       ! TRANSLATION
!!$       delta = 0.0d0
!!$       do i = 1, 3
!!$          temp_a = tau_local(i, nope)
!!$          temp_b = tau_temp(i)
!!$          if (dabs(temp_a) <=1.d-8 .and. dabs(temp_b) <=1.d-8) cycle
!!$          temp = mod(temp_a - temp_b + 2.0d0, 1.0d0)
!!$          if (temp > 0.9999999d0) cycle
!!$          temp = temp / max(abs(temp_a), abs(temp_b))
!!$          delta = delta + temp * temp
!!$       end do
!!$
!!$
!!$       if (delta < 1.0d-3) then
!!$          operatorExists = .true.  ! found
!!$          if(nopr_local+1 == 4891 .or. nopr_local+1 == 4892 .or. nopr_local+1==4893  &
!!$               & .or. nopr_local+1==4900 .or.nopr_local+1==4901 .or. nopr_local+1==4902) then
!!$             write(6,'(" ** operatorExists, nope = ",i8, " (iop, jop) = (",i8,i8," )")') nope, iop, jop
!!$          end if
!!$          if(nopr_local_temp >= 5190) then
!!$             if(jop == nopr_local_temp) then
!!$                write(6,'(" operatorExists, nope = ",i8, " (iop, jop) = (",i8,i8," )")') nope, iop, jop
!!$                call flush(6)
!!$             end if
!!$          end if
!!$          return
!!$       end if
!!$    enddo
!!$    operatorExists = .false.  ! not found
!!$  end function operatorExists
!!$!-----------------------------------------------------------------------------------
!!$end subroutine fd_symmetrize_sc2

subroutine set_iatomn_default(elen,ntyp, speciesname, iatomn)
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in) :: elen,ntyp
  character(len=elen), dimension(ntyp), intent(in) :: speciesname
  real(kind=DP), dimension(ntyp), intent(out) :: iatomn
  character(len=elen), allocatable :: ename_table(:)
  integer :: it,ie,pos,ii
  integer :: nmax_elem = 120

  allocate( ename_table( nmax_elem ) );ename_table=''
!  call set_table( 4,  1, 'H', 'He', 'Li', 'Be', 'B' )
  ename_table(1)  = 'H'
  ename_table(2)  = 'He'
  ename_table(3)  = 'Li'
  ename_table(4)  = 'Be'
  ename_table(5)  = 'B'
  ename_table(6)  = 'C'
  ename_table(7)  = 'N'
  ename_table(8)  = 'O'
  ename_table(9)  = 'F'
  ename_table(10) = 'Ne'
  ename_table(11) = 'Na'
  ename_table(12) = 'Mg'
  ename_table(13) = 'Al'
  ename_table(14) = 'Si'
  ename_table(15) = 'P'
  ename_table(16) = 'S'
  ename_table(17) = 'Cl'
  ename_table(18) = 'Ar'
  ename_table(19) = 'K'
  ename_table(20) = 'Ca'
  ename_table(21) = 'Sc'
  ename_table(22) = 'Ti'
  ename_table(23) = 'V'
  ename_table(24) = 'Cr'
  ename_table(25) = 'Mn'
  ename_table(26) = 'Fe'
  ename_table(27) = 'Co'
  ename_table(28) = 'Ni'
  ename_table(29) = 'Cu'
  ename_table(30) = 'Zn'
  ename_table(31) = 'Ga'
  ename_table(32) = 'Ge'
  ename_table(33) = 'As'
  ename_table(34) = 'Se'
  ename_table(35) = 'Br'
  ename_table(36) = 'Kr'
  ename_table(37) = 'Rb'
  ename_table(38) = 'Sr'
  ename_table(39) = 'Y'
  ename_table(40) = 'Zr'
  ename_table(41) = 'Nb'
  ename_table(42) = 'Mo'
  ename_table(43) = 'Tc'
  ename_table(44) = 'Ru'
  ename_table(45) = 'Rh'
  ename_table(46) = 'Pd'
  ename_table(47) = 'Ag'
  ename_table(48) = 'Cd'
  ename_table(49) = 'In'
  ename_table(50) = 'Sn'
  ename_table(51) = 'Sb'
  ename_table(52) = 'Te'
  ename_table(53) = 'I'
  ename_table(54) = 'Xe'
  ename_table(55) = 'Cs'
  ename_table(56) = 'Ba'
  ename_table(57) = 'La'
  ename_table(58) = 'Ce'
  ename_table(59) = 'Pr'
  ename_table(60) = 'Nd'
  ename_table(61) = 'Pm'
  ename_table(62) = 'Sm'
  ename_table(63) = 'Eu'
  ename_table(64) = 'Gd'
  ename_table(65) = 'Tb'
  ename_table(66) = 'Dy'
  ename_table(67) = 'Ho'
  ename_table(68) = 'Er'
  ename_table(69) = 'Tm'
  ename_table(70) = 'Yb'
  ename_table(71) = 'Lu'
  ename_table(72) = 'Hf'
  ename_table(73) = 'Ta'
  ename_table(74) = 'W'
  ename_table(75) = 'Re'
  ename_table(76) = 'Os'
  ename_table(77) = 'Ir'
  ename_table(78) = 'Pt'
  ename_table(79) = 'Au'
  ename_table(80) = 'Hg'
  ename_table(81) = 'Tl'
  ename_table(82) = 'Pb'
  ename_table(83) = 'Bi'
  ename_table(84) = 'Po'
  ename_table(85) = 'At'
  ename_table(86) = 'Rn'
  ename_table(87) = 'Fr'
  ename_table(88) = 'Ra'
  ename_table(89) = 'Ac'
  ename_table(90) = 'Th'
  ename_table(91) = 'Pa'
  ename_table(92) = 'U'
  ename_table(93) = 'Np'
  ename_table(94) = 'Pu'
  ename_table(95) = 'Am'
  ename_table(96) = 'Cm'
  ename_table(97) = 'Bk'
  ename_table(98) = 'Cf'
  ename_table(99) = 'Es'
  ename_table(100) = 'Fm'
  ename_table(101) = 'Md'
  ename_table(102) = 'No'
  ename_table(103) = 'Lr'
  ename_table(104) = 'Rf'
  ename_table(105) = 'Db'
  ename_table(106) = 'Sg'
  ename_table(107) = 'Bh'
  ename_table(108) = 'Hs'
  ename_table(109) = 'Mt'
  ename_table(110) = 'Ds'
  ename_table(111) = 'Rg'
  ename_table(112) = 'Cn'
  ename_table(113) = 'Nh'
  ename_table(114) = 'Fl'
  ename_table(115) = 'Mc'
  ename_table(116) = 'Lv'
  ename_table(117) = 'Ts'
  ename_table(118) = 'Og'

  iatomn = 0 
  do it=1,ntyp
     do ie=1,nmax_elem 
        if(trim(speciesname(it)) == trim(ename_table(ie))) then
           iatomn(it) = real(ie)
        else
           pos = index(trim(speciesname(it)),trim(ename_table(ie)))
           if (pos.eq.1) then
              read(speciesname(it)(len(trim(ename_table(ie)))+1:len(speciesname(it))),*,err=90,end=90) ii
              iatomn(it) = real(ie)
90            continue
           endif
        endif
     enddo
  enddo
  deallocate(ename_table) 

end subroutine set_iatomn_default

subroutine set_covrad_default(ntyp, iatomn, covrads)
  use m_Const_Parameters, only : DP,BOHR
  implicit none
!
  integer, intent(in) :: ntyp
  real(kind=DP), intent(in) :: iatomn(ntyp)
  real(kind=DP), intent(inout) :: covrads(ntyp)
  real(kind=DP) :: ang2bohr
  real(kind=DP), allocatable, dimension(:) :: cov_table
!
  integer :: it,i
  integer :: nmax_elem = 120
  ang2bohr = 1.d0/BOHR
  allocate(cov_table(nmax_elem))

  cov_table(:) = 1.d0

  cov_table(1) = 0.32d0
  cov_table(2) = 0.93d0
  cov_table(3) = 1.23d0
  cov_table(4) = 0.9d0
  cov_table(5) = 0.82d0
  cov_table(6) = 0.77d0
  cov_table(7) = 0.75d0
  cov_table(8) = 0.73d0
  cov_table(9) = 0.72d0
  cov_table(10) = 0.71d0
  cov_table(11) = 1.54d0
  cov_table(12) = 1.36d0
  cov_table(13) = 1.18d0
  cov_table(14) = 1.11d0
  cov_table(15) = 1.06d0
  cov_table(16) = 1.02d0
  cov_table(17) = 0.99d0
  cov_table(18) = 0.98d0
  cov_table(19) = 2.03d0
  cov_table(20) = 1.74d0
  cov_table(21) = 1.44d0
  cov_table(22) = 1.32d0
  cov_table(23) = 1.22d0
  cov_table(24) = 1.18d0
  cov_table(25) = 1.17d0
  cov_table(26) = 1.17d0
  cov_table(27) = 1.16d0
  cov_table(28) = 1.15d0
  cov_table(29) = 1.17d0
  cov_table(30) = 1.25d0
  cov_table(31) = 1.26d0
  cov_table(32) = 1.22d0
  cov_table(33) = 1.2d0
  cov_table(34) = 1.16d0
  cov_table(35) = 1.14d0
  cov_table(36) = 1.89d0
  cov_table(37) = 2.16d0
  cov_table(38) = 1.91d0
  cov_table(39) = 1.62d0
  cov_table(40) = 1.45d0
  cov_table(41) = 1.34d0
  cov_table(42) = 1.3d0
  cov_table(43) = 1.27d0
  cov_table(44) = 1.25d0
  cov_table(45) = 1.25d0
  cov_table(46) = 1.28d0
  cov_table(47) = 1.34d0
  cov_table(48) = 1.41d0
  cov_table(49) = 1.44d0
  cov_table(50) = 1.41d0
  cov_table(51) = 1.4d0
  cov_table(52) = 1.36d0
  cov_table(53) = 1.33d0
  cov_table(54) = 1.31d0
  cov_table(55) = 2.35d0
  cov_table(56) = 1.98d0
  cov_table(57) = 1.25d0
  cov_table(58) = 1.65d0
  cov_table(59) = 1.65d0
  cov_table(60) = 1.64d0
  cov_table(61) = 1.63d0
  cov_table(62) = 1.62d0
  cov_table(63) = 1.85d0
  cov_table(64) = 1.61d0
  cov_table(65) = 1.59d0
  cov_table(66) = 1.59d0
  cov_table(67) = 1.58d0
  cov_table(68) = 1.57d0
  cov_table(69) = 1.56d0
  cov_table(70) = 1.7d0
  cov_table(71) = 1.56d0
  cov_table(72) = 1.44d0
  cov_table(73) = 1.34d0
  cov_table(74) = 1.37d0
  cov_table(75) = 1.28d0
  cov_table(76) = 1.26d0
  cov_table(77) = 1.27d0
  cov_table(78) = 1.3d0
  cov_table(79) = 1.34d0
  cov_table(80) = 1.49d0
  cov_table(81) = 1.48d0
  cov_table(82) = 1.47d0
  cov_table(83) = 1.46d0
  cov_table(84) = 1.53d0
  cov_table(85) = 1.47d0
  cov_table(86) = 1.42d0
  cov_table(87) = 2.23d0
  cov_table(88) = 2.01d0
  cov_table(89) = 1.86d0
  cov_table(90) = 1.75d0
  cov_table(91) = 1.69d0
  cov_table(92) = 1.7d0
  cov_table(93) = 1.71d0
  cov_table(94) = 1.72d0
  cov_table(95) = 1.66d0
  cov_table(96) = 1.66d0
  cov_table(97) = 1.66d0
  cov_table(98) = 1.68d0
  cov_table(99) = 1.65d0
  cov_table(100) = 1.67d0
  cov_table(101) = 1.73d0
  cov_table(102) = 1.76d0
  cov_table(103) = 1.61d0
  cov_table(104) = 1.57d0
  cov_table(105) = 1.49d0
  cov_table(106) = 1.43d0
  cov_table(107) = 1.41d0
  cov_table(108) = 1.34d0
  cov_table(109) = 1.29d0

  do it=1,ntyp
     i = nint( iatomn(it) )
     covrads(it) = cov_table(i) * ang2bohr
  enddo

  deallocate(cov_table)

end subroutine set_covrad_default

! === KT_add ======= 13.3B
!*** ref. http://www.caslab.com/List-of-Elements/ *****
!
subroutine set_amion_default( ntyp, iatomn, amion )
  use m_Const_Parameters, only : DP, AMU
  implicit none
!
  integer, intent(in) :: ntyp
  real(kind=DP), intent(in) :: iatomn(ntyp)
  real(kind=DP), intent(inout) :: amion(ntyp)
!
  integer :: nmax_elem = 120
!
  integer :: it, i
  real(kind=DP) :: coeff
  real(kind=DP), allocatable :: mass_table(:)
!
  allocate( mass_table( nmax_elem ) )

  call set_table(   1, 1.00794d0, 4.002602d0, 6.941d0, 9.012182d0, 10.811d0 )
  call set_table(   6, 12.0107d0, 14.0067d0, 15.9994d0, 18.9984032d0, 20.1797d0 )
  call set_table(  11, 22.98976928d0, 24.3050d0, 26.9815386d0, 28.0855d0, 30.973762d0 )
  call set_table(  16, 32.065d0, 35.453d0, 39.948d0, 39.0983d0, 40.078d0 )

  call set_table(  21, 44.955912d0, 47.867d0, 50.9415d0, 51.9961d0, 54.938045d0 )
  call set_table(  26, 55.845d0, 58.933195d0, 58.6934d0, 63.546d0, 65.409d0 )
  call set_table(  31, 69.723d0, 72.64d0, 74.92160d0, 78.96d0, 79.904d0 )
  call set_table(  36, 83.798d0, 85.4678d0, 87.62d0, 88.90585d0, 91.224d0 )

  call set_table(  41, 92.90638d0, 95.94d0, 97.9072d0, 101.07d0, 102.90550d0 )
  call set_table(  46, 106.42d0, 107.8682d0, 112.411d0, 114.818d0, 118.710d0 )
  call set_table(  51, 121.760d0, 127.60d0, 126.90447d0, 131.293d0, 132.9054519d0 )
  call set_table(  56, 137.327d0, 138.90547d0, 140.116d0, 140.90765d0, 144.242d0 )

  call set_table(  61, 144.9127d0, 150.36d0, 151.964d0, 157.25d0, 158.92535d0 )
  call set_table(  66, 162.500d0, 164.93032d0, 167.259d0, 168.93421d0, 173.04d0 )
  call set_table(  71, 174.967d0, 178.49d0, 180.94788d0, 183.84d0, 186.207d0 )
  call set_table(  76, 190.23d0, 192.217d0, 195.084d0, 196.966569d0, 200.59d0 )

  call set_table(  81, 204.3833d0, 207.2d0, 208.98040d0, 208.9824d0, 209.9871d0 )
  call set_table(  86, 222.0176d0, 223.0197d0, 226.0254d0, 227.0277d0, 232.03806d0 )
  call set_table(  91, 231.03588d0, 238.02891d0, 237.0482d0, 244.0642d0, 243.0614d0 )
  call set_table(  96, 247.0704d0, 247.0703d0, 251.0796d0, 252.0830d0, 257.0951d0 )

  call set_table( 101, 258.0984d0, 259.1010d0, 262.1097d0, 261.1088d0, 262d0 )
  call set_table( 106, 266d0, 264d0, 277d0, 268d0, 271d0 )
  call set_table( 111, 272d0, 285d0, 284d0, 289d0, 288d0 )
  call set_table( 116, 292d0, 293d0, 294d0, 500d0, 500d0 )

!
  coeff = 1.0d0 / AMU  ! 1.66053d-27 /9.1093897d-31
  Do it=1, ntyp
     if ( amion(it) < -99.0 ) then
        i = nint( iatomn(it) )
        amion(it) = mass_table(i) *coeff
     endif
  End do

  deallocate( mass_table )

contains

  subroutine set_table( istart, c1, c2, c3, c4, c5 )
    integer, intent(in) :: istart
    real(kind=DP), intent(in) :: c1, c2, c3, c4, c5

    mass_table( istart ) = c1
    mass_table( istart +1 ) = c2
    mass_table( istart +2 ) = c3
    mass_table( istart +3 ) = c4
    mass_table( istart +4 ) = c5
  end subroutine set_table

end subroutine set_amion_default
! ========== 13.3B

subroutine get_ac_from_dynm(nfdynm,targetframe,natm,nchar,ntyp,avec,bvec,cvec,speciesname,ityp,cps,cpd)
  implicit none
  integer, intent(in) :: nfdynm,targetframe,natm,nchar,ntyp
  real(kind=8), dimension(3), intent(out) :: avec,bvec,cvec
  character(len=nchar), dimension(ntyp), intent(out) :: speciesname
  integer, dimension(natm), intent(out) :: ityp
  real(kind=8), dimension(natm,3), intent(out) :: cps,cpd
  integer, parameter :: ncharline=256
  integer, parameter :: ncharword=20
  integer :: nwords = 20
  character(len=ncharline) :: line
  character(len=ncharword), allocatable, dimension(:) :: output
  integer :: i
  integer :: tmpint
  logical :: in_header, in_data
  integer :: ntypc,currat,iframe
  integer :: nword_out

  allocate(output(nwords))
  in_header = .false.
  in_data = .false.
  ntypc = 0
  iframe = 0
  avec=0.d0;bvec=0.d0;cvec=0.d0
  speciesname='';ityp=0.d0;cps=0.d0
  do 
    read(unit=nfdynm,end=100,fmt='(a)') line
    call split(ncharline,ncharword,nwords,line,output,nword_out)
    if(.not. in_header .and. trim(output(1)) == '#')  then
      in_header = .true.
      cycle
    endif
    if(in_header) then
      if(trim(output(2)) == 'a_vector')then
        ntypc = 0
        currat = 0
        read(output(4),*) avec(1)
        read(output(5),*) avec(2)
        read(output(6),*) avec(3)
      else if(trim(output(2)) == 'b_vector')then
        read(output(4),*) bvec(1)
        read(output(5),*) bvec(2)
        read(output(6),*) bvec(3)
      else if(trim(output(2)) == 'c_vector')then
        read(output(4),*) cvec(1)
        read(output(5),*) cvec(2)
        read(output(6),*) cvec(3)
      else if(trim(output(2)) == '(natm->type)')then
        do i=3,nwords
          if(len(trim(output(i)))==0) exit
          ntypc = ntypc + 1
          read(output(i),*) ityp(ntypc)
        enddo
      else if(trim(output(2)) == '(speciesname)')then
        read(output(3),*) tmpint 
        speciesname(tmpint)(1:nchar) = output(5)(1:nchar)
      endif
    endif
    if(output(1) == 'cps')then
      in_header = .false.
      in_data = .true.
      ntypc = 0
      currat = 0
      cycle
    endif
    if(in_data)then
      currat = currat+1
      read(output(2),*) cps(currat,1)
      read(output(3),*) cps(currat,2)
      read(output(4),*) cps(currat,3)
      if(nword_out == 10)then
        read(output(5),*) cpd(currat,1)
        read(output(6),*) cpd(currat,2)
        read(output(7),*) cpd(currat,3)
      endif
      if(currat == natm)then
        iframe = iframe+1
        in_data = .false.
        if(targetframe .eq.iframe) then
           exit
        endif
      endif
    endif
  enddo

100 continue
  deallocate(output)
end subroutine get_ac_from_dynm

subroutine get_numat_and_ntyp(nfdynm,natm,ntyp)
  implicit none
  integer, intent(in) :: nfdynm
  integer, intent(out) :: natm,ntyp
  integer, parameter :: ncharline=256
  integer, parameter :: ncharword=20
  integer :: nwords = 20
  integer :: nword_out
  character(len=ncharline) :: line
  character(len=ncharword), allocatable, dimension(:) :: output
  logical :: in_header
  integer :: i

  natm = 0;ntyp = 0
  allocate(output(nwords))
  in_header = .false.
  do 
    read(unit=nfdynm,end=100,fmt='(a)') line
    call split(ncharline,ncharword,nwords,line,output,nword_out)
    if(.not. in_header .and. trim(output(1)) == '#')  then
      in_header = .true.
      cycle
    endif
    if(in_header) then
      if(trim(output(2)) == 'ntyp')then
        read(output(4),*) ntyp
        read(output(7),*) natm
        return
      endif
    endif
  enddo
100 continue
  write(0,*) '!** failed to resolve natm and ntyp'
  return
end subroutine get_numat_and_ntyp

subroutine split(ncharline,ncharword,nwords,input,output,nword_out) 
  implicit none
  integer, intent(in) :: ncharline,ncharword,nwords
  character(len=ncharline), intent(in) :: input
  character(len=ncharword), dimension(nwords), intent(out) :: output
  integer, intent(out) :: nword_out
  integer :: i, currword,currind
  logical :: was_ws
  was_ws = .true.
  currword = 0
  currind = 0
  do i=1,nwords
    output(i) = ''
  enddo
  do i=1,ncharline
    if (input(i:i) .eq. ' ') then
      was_ws = .true.
      cycle
    endif
    if (was_ws) then
      currword = currword + 1
      was_ws = .false.
      currind = 0
    endif
    currind = currind+1
    output(currword)(currind:currind) = input(i:i)
  enddo
  nword_out = currword
end subroutine split

