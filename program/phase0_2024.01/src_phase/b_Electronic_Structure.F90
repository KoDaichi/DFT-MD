!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  SUBROUINE:  crngabc4, setglist4, zf_list_s, calc_phase2, calc_phase
!              substitute_il3
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
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
! $Id: b_Electronic_Structure.F90 633 2020-12-01 05:11:03Z jkoga $
subroutine substitute_il3(n,il3)
  implicit none
  integer, intent(in)               :: n
  integer, intent(out), dimension(n) :: il3

  integer i, l
  l = 0
  do i = 1, n
     if(i > (l+1)**2) l = l + 1
     il3(i) = l
  end do
end subroutine substitute_il3

subroutine calc_phase(ia,natm,pos,kgp,ibaik,ngabc,kg1,nbase,kd,zfcos,zfsin)
  use m_Const_Parameters, only : PAI2, DP
  implicit none
  integer, intent(in)                          :: ia   ! #atom
  integer, intent(in)                          :: natm ! 1st dim. of pos
  real(kind=DP), intent(in), dimension(natm,3) :: pos  ! positions of atoms
  integer, intent(in)                          :: kgp  ! 1st dim. of ngabc
  integer, intent(in)                          :: ibaik! range of operation
  integer, intent(in), dimension(kgp,3)        :: ngabc! g-vectors
  integer, intent(in)                          :: kg1  ! 1st dim. of nbase
  integer, intent(in), dimension(kg1)          :: nbase! pointer to a G-vector set
  integer, intent(in)                          :: kd   ! dim. of zfcos, zfsin
  real(kind=DP), intent(out), dimension(kd)    :: zfcos, zfsin ! phase

  integer       :: i, nb
  real(kind=DP) ::  fx, fy, fz, ph
  fx = pos(ia,1)*PAI2
  fy = pos(ia,2)*PAI2
  fz = pos(ia,3)*PAI2
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
  do i = 1,ibaik
     nb = nbase(i)
     ph = ngabc(nb,1)*fx+ngabc(nb,2)*fy+ngabc(nb,3)*fz
     zfcos(i) = dcos(ph)
     zfsin(i) = dsin(ph)
  end do
end subroutine calc_phase

!$$#ifndef PARA3D
subroutine calc_phasek_b(natm,pos,n_ialist0,ia_list,kgp,iiba,ngabc,kg1,nbase,kd_adj,zfcos_x,zfsin_x)
  use m_Const_Parameters, only : PAI2, DP
  implicit none
  integer, intent(in)                          :: natm       ! 1st dim. of pos
  real(kind=DP), intent(in), dimension(natm,3) :: pos        ! positions of atoms
  integer, intent(in)                          :: n_ialist0  ! #atom
  integer, intent(in), dimension(n_ialist0)    :: ia_list
  integer, intent(in)                          :: kgp        ! 1st dim. of ngabc
  integer, intent(in)                          :: iiba       ! range of operation
  integer, intent(in), dimension(kgp,3)        :: ngabc      ! g-vectors
  integer, intent(in)                          :: kg1        ! 1st dim. of nbase
  integer, intent(in), dimension(kg1)          :: nbase      ! pointer to a G-vector set
  integer, intent(in)                          :: kd_adj     ! dim. of zfcos, zfsin
  real(kind=DP), intent(out), dimension(kd_adj,n_ialist0)   :: zfcos_x, zfsin_x ! phase

  integer       :: i, nb, ia,iap
  real(kind=DP) ::  fx, fy, fz, ph

#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption parallel
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
#ifdef NEC_TUNE_SMP
!CDIR SELECT(CONCUR)
#endif
  do iap = 1, n_ialist0
     ia = ia_list(iap)
     fx = pos(ia,1)*PAI2
     fy = pos(ia,2)*PAI2
     fz = pos(ia,3)*PAI2
     do i = 1,iiba
        nb = nbase(i)
        ph = ngabc(nb,1)*fx+ngabc(nb,2)*fy+ngabc(nb,3)*fz
        zfcos_x(i,iap) = dcos(ph)
        zfsin_x(i,iap) = dsin(ph)
     end do
   end do
end subroutine calc_phasek_b
!$$#endif


subroutine calc_phase2(natm,pos,ia,kgp,ngabc,ista_kngp,iend_kngp,zfcos,zfsin)
  use m_Const_Parameters, only : PAI2, DP
  implicit none
  integer, intent(in)                          :: natm ! 1st dim. of pos
  real(kind=DP), intent(in), dimension(natm,3) :: pos  ! positions of atoms
  integer, intent(in)                          :: ia   ! #atom
  integer, intent(in)                          :: kgp  ! range of operation
  integer, intent(in), dimension(kgp,3)        :: ngabc! g-vectors
  integer, intent(in)                          :: ista_kngp, iend_kngp
  real(kind=DP),intent(out),dimension(ista_kngp:iend_kngp) :: zfcos, zfsin ! phase

  integer       ::  i
  real(kind=DP) :: fx, fy, fz, ph
  integer       ::  iend  !mpi
  fx = pos(ia,1)*PAI2
  fy = pos(ia,2)*PAI2
  fz = pos(ia,3)*PAI2
  iend = iend_kngp
  if( iend_kngp > kgp ) iend = kgp
  if( ista_kngp <= iend ) then
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
     do i = ista_kngp, iend  !for mp
        ph = ngabc(i,1)*fx+ngabc(i,2)*fy+ngabc(i,3)*fz
        zfcos(i) = dcos(ph)
        zfsin(i) = dsin(ph)
     end do
  endif
!xocl end spread
end subroutine calc_phase2

subroutine calc_phase_b(natm,pos,ia_list,n_ialist,kgp,ngabc,ista_kngp,iend_kngp,ista_kngp_adj,iend_kngp_adj,zfcos_x,zfsin_x)
  use m_Const_Parameters, only : PAI2, DP
  implicit none
  integer, intent(in) ::                          natm ! 1st dim. of pos
  real(kind=DP), intent(in), dimension(natm,3) :: pos  ! positions of atoms
  integer, intent(in) ::                          n_ialist ! #atom
  integer, intent(in), dimension(n_ialist) ::     ia_list(n_ialist) ! list of #atom
  integer, intent(in) ::                          kgp  ! range of operation
  integer, intent(in), dimension(kgp,3)        :: ngabc! g-vectors
  integer, intent(in)                          :: ista_kngp, iend_kngp, ista_kngp_adj, iend_kngp_adj
  real(kind=DP),intent(out),dimension(ista_kngp_adj:iend_kngp_adj,n_ialist) :: zfcos_x, zfsin_x ! phase

  integer       ::  i, ia, ip, is
  real(kind=DP) :: fx_1, fy_1, fz_1, ph_1, fx_2, fy_2, fz_2, ph_2, fx_3, fy_3, fz_3, ph_3, fx_4, fy_4, fz_4, ph_4
  integer       ::  iend  !mpi
  if(n_ialist >= 1) then
     fx_1 = pos(ia_list(1),1)*PAI2
     fy_1 = pos(ia_list(1),2)*PAI2
     fz_1 = pos(ia_list(1),3)*PAI2
  end if
  if(n_ialist >= 2) then
     fx_2 = pos(ia_list(2),1)*PAI2
     fy_2 = pos(ia_list(2),2)*PAI2
     fz_2 = pos(ia_list(2),3)*PAI2
  end if
  if(n_ialist >= 3) then
     fx_3 = pos(ia_list(3),1)*PAI2
     fy_3 = pos(ia_list(3),2)*PAI2
     fz_3 = pos(ia_list(3),3)*PAI2
  end if
  if(n_ialist >= 4) then
     fx_4 = pos(ia_list(4),1)*PAI2
     fy_4 = pos(ia_list(4),2)*PAI2
     fz_4 = pos(ia_list(4),3)*PAI2
  end if

  iend = iend_kngp
  if( iend_kngp > kgp ) iend = kgp
  if( ista_kngp <= iend ) then
     is = ista_kngp_adj - ista_kngp
     if(n_ialist == 1) then
!cdir parallel do private(i,ph_1,ip)
        do i = ista_kngp, iend  !for mp
           ph_1 = ngabc(i,1)*fx_1+ngabc(i,2)*fy_1+ngabc(i,3)*fz_1
           ip = i + is
           zfcos_x(ip,1) = dcos(ph_1)
           zfsin_x(ip,1) = dsin(ph_1)
        end do
     else if(n_ialist == 2) then
!cdir parallel do private(i,ph_1,ph_2,ip)
        do i = ista_kngp, iend  !for mp
           ph_1 = ngabc(i,1)*fx_1+ngabc(i,2)*fy_1+ngabc(i,3)*fz_1
           ph_2 = ngabc(i,1)*fx_2+ngabc(i,2)*fy_2+ngabc(i,3)*fz_2
           ip = i + is
           zfcos_x(ip,1) = dcos(ph_1)
           zfsin_x(ip,1) = dsin(ph_1)
           zfcos_x(ip,2) = dcos(ph_2)
           zfsin_x(ip,2) = dsin(ph_2)
        end do
     else if(n_ialist == 3) then
!cdir parallel do private(i,ph_1,ph_2,ph_3,ip)
        do i = ista_kngp, iend  !for mp
           ph_1 = ngabc(i,1)*fx_1+ngabc(i,2)*fy_1+ngabc(i,3)*fz_1
           ph_2 = ngabc(i,1)*fx_2+ngabc(i,2)*fy_2+ngabc(i,3)*fz_2
           ph_3 = ngabc(i,1)*fx_3+ngabc(i,2)*fy_3+ngabc(i,3)*fz_3
           ip = i + is
           zfcos_x(ip,1) = dcos(ph_1)
           zfsin_x(ip,1) = dsin(ph_1)
           zfcos_x(ip,2) = dcos(ph_2)
           zfsin_x(ip,2) = dsin(ph_2)
           zfcos_x(ip,3) = dcos(ph_3)
           zfsin_x(ip,3) = dsin(ph_3)
        end do
     else if(n_ialist >= 4) then
!cdir parallel do private(i,ph_1,ph_2,ph_3,ph_4,ip)
        do i = ista_kngp, iend  !for mp
           ph_1 = ngabc(i,1)*fx_1+ngabc(i,2)*fy_1+ngabc(i,3)*fz_1
           ph_2 = ngabc(i,1)*fx_2+ngabc(i,2)*fy_2+ngabc(i,3)*fz_2
           ph_3 = ngabc(i,1)*fx_3+ngabc(i,2)*fy_3+ngabc(i,3)*fz_3
           ph_4 = ngabc(i,1)*fx_4+ngabc(i,2)*fy_4+ngabc(i,3)*fz_4
           ip = i + is
           zfcos_x(ip,1) = dcos(ph_1)
           zfsin_x(ip,1) = dsin(ph_1)
           zfcos_x(ip,2) = dcos(ph_2)
           zfsin_x(ip,2) = dsin(ph_2)
           zfcos_x(ip,3) = dcos(ph_3)
           zfsin_x(ip,3) = dsin(ph_3)
           zfcos_x(ip,4) = dcos(ph_4)
           zfsin_x(ip,4) = dsin(ph_4)
        end do
     end if
     if(n_ialist > 4 ) then
        do ia = 5, n_ialist
           fx_1 = pos(ia_list(ia),1)*PAI2
           fy_1 = pos(ia_list(ia),2)*PAI2
           fz_1 = pos(ia_list(ia),3)*PAI2
!cdir parallel do private(i,ph_1,ip)
           do i = ista_kngp, iend   !MPI
              ph_1 = ngabc(i,1)*fx_1+ngabc(i,2)*fy_1+ngabc(i,3)*fz_1
              ip = i + is
              zfcos_x(ip,ia) = dcos(ph_1)
              zfsin_x(ip,ia) = dsin(ph_1)
           end do
        end do
     end if
  endif
!xocl end spread
end subroutine calc_phase_b

subroutine zf_list_s(n_min,n_max,matm,natm,f,ia,zfcos,zfsin)
  use m_Const_Parameters, only : DP
  implicit none
  integer, intent(in)      :: n_min,n_max,matm,natm
  real(kind=DP),intent(in) :: f
  integer, intent(in)      :: ia
  real(DP),intent(out),dimension(n_min:n_max,matm+1:natm)::zfcos,zfsin

  integer :: i
  do i = n_min, n_max
     zfcos(i,ia) = dcos(i*f)
     zfsin(i,ia) = dsin(i*f)
  end do
end subroutine zf_list_s

subroutine setglist4(n_min1,n_max1,n_min2,n_max2,n_min3,n_max3 &
     &     ,nbase,ngabc,kgp,iiba, nglist)
  implicit none
  integer,intent(in) :: n_min1,n_max1,n_min2,n_max2,n_min3,n_max3
  integer,intent(in) :: kgp,iiba,nbase(*),ngabc(kgp,3)
  integer,intent(out):: nglist(n_min1:n_max1,n_min2:n_max2,n_min3:n_max3)

  integer :: i, nb, nb1,nb2,nb3

  nglist = -1

  do i = 1, iiba
     nb = nbase(i)
     nb1 = ngabc(nb,1)
     nb2 = ngabc(nb,2)
     nb3 = ngabc(nb,3)
     nglist(nb1,nb2,nb3) = i
  enddo
end subroutine setglist4


subroutine crngabc4(n_min1,n_max1,n_min2,n_max2,n_min3,n_max3 &
     &     ,nglist,kg1,nngabc,newp)
  implicit none
  integer, intent(in) :: n_min1,n_max1,n_min2,n_max2,n_min3,n_max3 &
       &, nglist(n_min1:n_max1,n_min2:n_max2,n_min3:n_max3), kg1
  integer, intent(out) :: nngabc(kg1,3),newp(kg1)
  
  integer :: ngcount, jcount,kcount,jp,kp,ip

  ngcount = 0
  do jcount = n_min2, n_max2
     jp = jcount
     do kcount = n_min3, n_max3
        kp = kcount
        do ip = n_min1, n_max1
           if(nglist(ip,jp,kp) /= -1) then
              ngcount = ngcount + 1
              nngabc(ngcount,1) = ip
              nngabc(ngcount,2) = jp
              nngabc(ngcount,3) = kp
              newp(ngcount) = nglist(ip,jp,kp)
           endif
           jp = jp + 1
           if(jp.eq.n_max2+1) jp = n_min2
           kp = kp + 1
           if(kp.eq.n_max3+1) kp = n_min3
        enddo
     enddo
  enddo
end subroutine crngabc4

! ==== experimental ===
subroutine sphset_upto_L6( cr2, isph2, mmt2 )
  use m_Const_Parameters,  only : DP, CMPLDP
  implicit none
!
  integer, intent(out) :: mmt2( 16, 16 )
  integer, intent(out) :: isph2( 16, 16, 7 )
  real(kind=DP), intent(out) :: cr2( 16, 16, 7 )
  
  real(kind=DP), allocatable :: gaunt_coeff(:,:,:)
  integer :: nsph1, nsph2, nsph3, num
!
  cr2 = 0.0d0;  isph2 = 0;  mmt2 = 0

  allocate( gaunt_coeff(16,49,16) );  gaunt_coeff = 0.0d0
  call set_gaunt_coeff( gaunt_coeff )

  Do nsph1=1, 16
     Do nsph3=1, 16
        num = 0
        Do nsph2=1, 49
!           if ( abs(gaunt_coeff(nsph1,nsph2,nsph3)) < 1.0D-8 ) cycle
           if ( abs(gaunt_coeff(nsph1,nsph2,nsph3)) < 1.0D-6 ) cycle
           num = num +1
           isph2(nsph1,nsph3,num) = nsph2
           cr2(nsph1,nsph3,num) = gaunt_coeff(nsph1,nsph2,nsph3)
        End Do
        mmt2(nsph1,nsph3) = num
!        write(*,*) "Num = ", nsph1, nsph3, num
     End Do
  End Do
  deallocate( gaunt_coeff )
end subroutine sphset_upto_L6

subroutine set_gaunt_coeff( gaunt_coeff )
  use m_Const_Parameters,  only : DP, CMPLDP, PAI4, zi
  implicit none

  real(kind=DP), intent(out) :: gaunt_coeff(16,49,16)

  integer :: l1, l2, l3, m1, m2, m3, im1, im2, im3
  integer :: nsph1, nsph2, nsph3
  real(kind=DP) :: c0, c1, c2, c3
  complex(kind=CMPLDP) :: zsum
  complex(kind=CMPLDP), allocatable :: matwk(:,:,:)
  complex(kind=CMPLDP), pointer :: matU1(:,:), matU2(:,:), matU3(:,:)
  real(kind=DP) :: wigner_3j

  real(kind=DP), allocatable :: matwk2(:,:,:)
  
  complex(kind=CMPLDP), target :: MatU_ylm_RC_L0( 2*0+1, -0:0 )
  complex(kind=CMPLDP), target :: MatU_ylm_RC_L1( 2*1+1, -1:1 )
  complex(kind=CMPLDP), target :: MatU_ylm_RC_L2( 2*2+1, -2:2 )
  complex(kind=CMPLDP), target :: MatU_ylm_RC_L3( 2*3+1, -3:3 )
  complex(kind=CMPLDP), target :: MatU_ylm_RC_L4( 2*4+1, -4:4 )
  complex(kind=CMPLDP), target :: MatU_ylm_RC_L5( 2*5+1, -5:5 )
  complex(kind=CMPLDP), target :: MatU_ylm_RC_L6( 2*6+1, -6:6 )

  gaunt_coeff = 0.0d0

  call set_matU_ylm
  
  allocate( matwk(-3:3,-6:6,-3:3) ); matwk = 0.0d0
  
  Do l1=0, 3
     select case(l1)
     case (0)
        matU1 =>  MatU_ylm_RC_L0
     case (1)
        matU1 =>  MatU_ylm_RC_L1
     case (2)
        matU1 =>  MatU_ylm_RC_L2
     case (3)
        matU1 =>  MatU_ylm_RC_L3
     end select
     
     Do l2=0, 6
        select case(l2)
        case (0)
           matU2 =>  MatU_ylm_RC_L0
        case (1)
           matU2 =>  MatU_ylm_RC_L1
        case (2)
           matU2 =>  MatU_ylm_RC_L2
        case (3)
           matU2 =>  MatU_ylm_RC_L3
        case (4)
           matU2 =>  MatU_ylm_RC_L4
        case (5)
           matU2 =>  MatU_ylm_RC_L5
        case (6)
           matU2 =>  MatU_ylm_RC_L6
        end select
        
        Do l3=0, 3
           select case(l3)
           case (0)
              matU3 =>  MatU_ylm_RC_L0
           case (1)
              matU3 =>  MatU_ylm_RC_L1
           case (2)
              matU3 =>  MatU_ylm_RC_L2
           case (3)
              matU3 =>  MatU_ylm_RC_L3
           end select
           
           c0 = ( 2*l1+1 )*( 2*l2 +1 )*( 2*l3 +1 ) /PAI4
           c0 = sqrt( c0 )
           
           matwk = 0.0d0
           Do m1=-l1, l1
              Do m2=-l2, l2
                 Do m3=-l3, l3
                    c1 = wigner_3j( dble(l1), 0.0d0, dble(l2), 0.0d0, dble(l3), 0.0d0 )
                    c2 = wigner_3j( dble(l1), -dble(m1), dble(l2), dble(m2), &
                         &          dble(l3), dble(m3) )
                    c3 = (-1)**m1
                    matwk(m1,m2,m3) = c0 *c1 *c2 *c3
                 End Do
              End Do
           End Do
           !
           Do im1=1, 2*l1 +1
              Do im2=1, 2*l2 +1
                 Do im3=1, 2*l3 +1
                    nsph1 = l1**2 +im1
                    nsph2 = l2**2 +im2
                    nsph3 = l3**2 +im3
                    zsum = 0.0d0
                    Do m1=-l1, l1
                       Do m2=-l2, l2
                          Do m3=-l3, l3
                             zsum = zsum +conjg(matU1(im1,m1)) *matU2(im2,m2) &
                                  &      *matU3(im3,m3) *matwk(m1,m2,m3)
                          End Do
                       End Do
                    ENd Do
                    gaunt_coeff(nsph1,nsph2,nsph3) = zsum
                 End Do
              ENd Do
           ENd Do
        End Do
     ENd Do
  End Do
  deallocate( matwk )

contains
  subroutine set_matU_ylm
    real(kind=DP) :: ctmp
    complex(kind=CMPLDP) :: ztmp
    !
    MatU_ylm_RC_L0 = 0.0d0
    MatU_ylm_RC_L1 = 0.0d0
    MatU_ylm_RC_L2 = 0.0d0
    MatU_ylm_RC_L3 = 0.0d0
    MatU_ylm_RC_L4 = 0.0d0
    MatU_ylm_RC_L5 = 0.0d0
    MatU_ylm_RC_L6 = 0.0d0
    ! -------
    ctmp = 1.0d0 / sqrt(2.0d0)

! ------------------------- ktDEBUG ------------ 2012/11/12
    ztmp = -ctmp * zi
!!    ztmp = ctmp * zi
! ------------------------- ktDEBUG ------------ 2012/11/12
    !
! s-orbital
    MatU_ylm_RC_L0(  1,0 ) = 1.0d0
! p-orbital
    MatU_ylm_RC_L1(  1,-1 ) =  ctmp;    MatU_ylm_RC_L1( 1,1 ) = -ctmp
    MatU_ylm_RC_L1(  2,-1 ) =  ztmp;    MatU_ylm_RC_L1( 2,1 ) =  ztmp
    MatU_ylm_RC_L1(  3, 0 ) = 1.0d0
! d-orbital
    MatU_ylm_RC_L2(  1, 0 ) = 1.0d0
    MatU_ylm_RC_L2(  2,-2 ) =  ctmp;    MatU_ylm_RC_L2( 2,2 ) =  ctmp
    MatU_ylm_RC_L2(  3,-2 ) =  ztmp;    MatU_ylm_RC_L2( 3,2 ) = -ztmp
    MatU_ylm_RC_L2(  4,-1 ) =  ztmp;    MatU_ylm_RC_L2( 4,1 ) =  ztmp
    MatU_ylm_RC_L2(  5,-1 ) =  ctmp;    MatU_ylm_RC_L2( 5,1 ) = -ctmp
! f-orbital
    MatU_ylm_RC_L3(  1, 0 ) = 1.0d0
    MatU_ylm_RC_L3(  2,-1 ) =  ctmp;   MatU_ylm_RC_L3( 2,1 ) = -ctmp
    MatU_ylm_RC_L3(  3,-1 ) =  ztmp;   MatU_ylm_RC_L3( 3,1 ) =  ztmp
    MatU_ylm_RC_L3(  4,-2 ) =  ctmp;   MatU_ylm_RC_L3( 4,2 ) =  ctmp
    MatU_ylm_RC_L3(  5,-2 ) =  ztmp;   MatU_ylm_RC_L3( 5,2 ) = -ztmp
    MatU_ylm_RC_L3(  6,-3 ) =  ctmp;   MatU_ylm_RC_L3( 6,3 ) = -ctmp
    MatU_ylm_RC_L3(  7,-3 ) =  ztmp;   MatU_ylm_RC_L3( 7,3 ) =  ztmp
! g-orbital
    MatU_ylm_RC_L4(  1, 0 ) = 1.0d0
    MatU_ylm_RC_L4(  2,-1 ) =  ctmp;   MatU_ylm_RC_L4( 2,1 ) = -ctmp
    MatU_ylm_RC_L4(  3,-1 ) =  ztmp;   MatU_ylm_RC_L4( 3,1 ) =  ztmp
    MatU_ylm_RC_L4(  4,-2 ) =  ctmp;   MatU_ylm_RC_L4( 4,2 ) =  ctmp
    MatU_ylm_RC_L4(  5,-2 ) =  ztmp;   MatU_ylm_RC_L4( 5,2 ) = -ztmp
    MatU_ylm_RC_L4(  6,-3 ) =  ctmp;   MatU_ylm_RC_L4( 6,3 ) = -ctmp
    MatU_ylm_RC_L4(  7,-3 ) =  ztmp;   MatU_ylm_RC_L4( 7,3 ) =  ztmp
    MatU_ylm_RC_L4(  8,-4 ) =  ctmp;   MatU_ylm_RC_L4( 8,4 ) =  ctmp
    MatU_ylm_RC_L4(  9,-4 ) =  ztmp;   MatU_ylm_RC_L4( 9,4 ) = -ztmp
! L=5
    MatU_ylm_RC_L5(  1, 0 ) = 1.0d0
    MatU_ylm_RC_L5(  2,-1 ) =  ctmp;   MatU_ylm_RC_L5(  2,1 ) = -ctmp
    MatU_ylm_RC_L5(  3,-1 ) =  ztmp;   MatU_ylm_RC_L5(  3,1 ) =  ztmp
    MatU_ylm_RC_L5(  4,-2 ) =  ctmp;   MatU_ylm_RC_L5(  4,2 ) =  ctmp
    MatU_ylm_RC_L5(  5,-2 ) =  ztmp;   MatU_ylm_RC_L5(  5,2 ) = -ztmp
    MatU_ylm_RC_L5(  6,-3 ) =  ctmp;   MatU_ylm_RC_L5(  6,3 ) = -ctmp
    MatU_ylm_RC_L5(  7,-3 ) =  ztmp;   MatU_ylm_RC_L5(  7,3 ) =  ztmp
    MatU_ylm_RC_L5(  8,-4 ) =  ctmp;   MatU_ylm_RC_L5(  8,4 ) =  ctmp
    MatU_ylm_RC_L5(  9,-4 ) =  ztmp;   MatU_ylm_RC_L5(  9,4 ) = -ztmp
    MatU_ylm_RC_L5( 10,-5 ) =  ctmp;   MatU_ylm_RC_L5( 10,5 ) = -ctmp
    MatU_ylm_RC_L5( 11,-5 ) =  ztmp;   MatU_ylm_RC_L5( 11,5 ) =  ztmp
! L=6
    MatU_ylm_RC_L6(  1, 0 ) = 1.0d0
    MatU_ylm_RC_L6(  2,-1 ) =  ctmp;   MatU_ylm_RC_L6(  2,1 ) = -ctmp
    MatU_ylm_RC_L6(  3,-1 ) =  ztmp;   MatU_ylm_RC_L6(  3,1 ) =  ztmp
    MatU_ylm_RC_L6(  4,-2 ) =  ctmp;   MatU_ylm_RC_L6(  4,2 ) =  ctmp
    MatU_ylm_RC_L6(  5,-2 ) =  ztmp;   MatU_ylm_RC_L6(  5,2 ) = -ztmp
    MatU_ylm_RC_L6(  6,-3 ) =  ctmp;   MatU_ylm_RC_L6(  6,3 ) = -ctmp
    MatU_ylm_RC_L6(  7,-3 ) =  ztmp;   MatU_ylm_RC_L6(  7,3 ) =  ztmp
    MatU_ylm_RC_L6(  8,-4 ) =  ctmp;   MatU_ylm_RC_L6(  8,4 ) =  ctmp
    MatU_ylm_RC_L6(  9,-4 ) =  ztmp;   MatU_ylm_RC_L6(  9,4 ) = -ztmp
    MatU_ylm_RC_L6( 10,-5 ) =  ctmp;   MatU_ylm_RC_L6( 10,5 ) = -ctmp
    MatU_ylm_RC_L6( 11,-5 ) =  ztmp;   MatU_ylm_RC_L6( 11,5 ) =  ztmp
    MatU_ylm_RC_L6( 12,-6 ) =  ctmp;   MatU_ylm_RC_L6( 12,6 ) =  ctmp
    MatU_ylm_RC_L6( 13,-6 ) =  ztmp;   MatU_ylm_RC_L6( 13,6 ) = -ztmp
  end subroutine set_matU_ylm

end subroutine set_gaunt_coeff

real(kind=DP) function wigner_3j( val_j1, val_m1, val_j2, val_m2, val_j3, val_m3 )
  use m_Const_Parameters,  only : DP, CMPLDP
  implicit none
  real(kind=DP), intent(in) :: val_j1, val_j2, val_j3, val_m1, val_m2, val_m3

  integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, k, kmin, kmax
  real(kind=DP) :: c1, c2, c3, c4, c5, c6, c7, c8, c9, c10
  real(kind=DP) :: d1, d2, term_n, term_s, coeff
  integer :: factorial

  if ( nint(val_m1 +val_m2 +val_m3) .ne. 0 ) then
     wigner_3j = 0.0d0;  return
  endif
  
  coeff = (-1)**( val_j1 -val_j2 -val_m3 )
  
  i1 = nint( val_j3 +val_j1 -val_j2 );    i2 = nint( val_j3 -val_j1 +val_j2 )
  i3 = nint( val_j1 +val_j2 -val_j3 );    i4 = nint( val_j3 -val_m3 )
  i5 = nint( val_j3 +val_m3 )
  i6 = nint( val_j1 +val_j2 +val_j3 +1 )
  i7 = nint( val_j1 -val_m1 );            i8 = nint( val_j1 +val_m1 )
  i9 = nint( val_j2 -val_m2 );            i10 = nint( val_j2 +val_m2 )
  
  c1 = factorial(i1);      c2 = factorial(i2);       c3 = factorial(i3)
  c4 = factorial(i4);      c5 = factorial(i5);       c6 = factorial(i6)
  c7 = factorial(i7);      c8 = factorial(i8);       c9 = factorial(i9)
  c10= factorial(i10)
  
  d1 = c1 *c2 *c3 *c4 *c5
  d2 = c6 *c7 *c8 *c9 *c10
!  write(*,*) "d1 = ", d1, d2
  
  term_n = sqrt( d1 /d2 )
  
  kmin = max( 0, nint(val_m1 -val_j1) )
  kmin = max( kmin, nint(val_j2 -val_j1 -val_m3) )
  kmax = min( nint(val_j2 +val_j3 +val_m1), nint(val_j3 -val_j1 +val_j2) )
  kmax = min( kmax, nint(val_j3 -val_m3) )
  
  term_s = 0.0d0
  Do k=kmin, kmax
     !    Do k=0, 100
     i2 = nint( val_j2 +val_j3 +val_m1 -k )
     i3 = nint( val_j1 -val_m1 +k )
     i5 = nint( val_j3 -val_j1 +val_j2 -k )
     i6 = nint( val_j3 -val_m3 -k )
     i7 = nint( val_j1 -val_j2 +val_m3 +k )
!       if ( i2 < 0 ) cycle
!       if ( i3 < 0 ) cycle
!       if ( i5 < 0 ) cycle
!       if ( i6 < 0 ) cycle
!       if ( i7 < 0 ) cycle

     c1 = (-1)**( val_j2 +val_m2 +k )
     c2 = factorial(i2);       c3 = factorial(i3);
     c4 = factorial(k)
     c5 = factorial(i5);       c6 = factorial(i6);       c7 = factorial(i7)
     d1 = c1 *c2 *c3
     d2 = c4 *c5 *c6 *c7
     term_s = term_s +d1 /d2
  End Do
  wigner_3j = coeff *term_n *term_s
  return
end function wigner_3j

integer recursive function factorial(n)
  integer, intent(in) :: n
  
  integer :: i
  
  factorial = 1
  Do i=1, n
     factorial = factorial *i
  End do
  return
end function factorial
  

subroutine calc_small_Wigner_function( val_j, val_m1, val_m2, beta, cret )
  use m_Const_Parameters, only : DP

  real(kind=DP), intent(in) :: val_j, val_m1, val_m2, beta
  real(kind=DP), intent(out) :: cret
  
  integer :: i1, i2, i3, i4, n, nmin, nmax, u1, u2
  real(kind=DP) :: c1, c2, c3, c4, d1, th, f1, f2, csum
  integer :: factorial
  
  i1 = nint( val_j +val_m1 );    i2 = nint( val_j -val_m1 )
  i3 = nint( val_j +val_m2 );    i4 = nint( val_j -val_m2 )
  
  c1 = factorial( i1 );    c2 = factorial( i2 )
  c1 = factorial( i1 );    c2 = factorial( i2 )
  c3 = factorial( i3 );    c4 = factorial( i4 )
  
  d1 = sqrt( c1 *c2 *c3 *c4 )
  
  nmin = max( 0, nint(val_m2 -val_m1) )
  nmax = min( nint(val_j -val_m1), nint(val_j+val_m2) )
  
  th = beta/2.0d0;    f1 = cos(th);    f2 = -sin(th)
  
  csum = 0.0D0
  Do n=nmin, nmax
     i1 = nint( val_j -val_m1 -n )
     i2 = nint( val_j +val_m2 -n )
     i3 = nint( n + val_m1 -val_m2 )
     c1 = factorial( i1 );       c2 = factorial( i2 )
     c3 = factorial( i3 );       c4 = factorial( n )
     u1 = nint( 2.0 *val_j +val_m2 -val_m1 -2 *n )
     u2 = nint( val_m1 -val_m2 +2 *n )
     csum = csum +(-1)**n *f1 **u1 *f2 **u2 /( c1*c2*c3*c4 )
  End do
  cret = csum *d1
end subroutine calc_small_Wigner_function
