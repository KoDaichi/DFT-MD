!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  MODULE: m_ES_initialWF
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
!
module m_ES_initialWF
!    ( m_ESIW)
! $Id: m_ES_initialWF.F90 570 2017-04-21 20:34:50Z yamasaki $
!
! The original subroutine was "rndzaj", which was coded with
! referring to Dr. K. Kato's program on 2nd Jul. 1994.
!
! "rndzaj" was modified to improve the quality of random
! numbers by T. Yamasaki.    17th Oct. 1994
!
! "rndzaj" was translated into "m_ESIW_by_randomnumbers" by
! T. Yamasaki, in 1999.
!

  use m_Electronic_Structure,only : zaj_l, neordr, nrvf_ordr
  use m_ES_nonlocal         ,only : m_ES_PAO_WFs
  use m_PlaneWaveBasisSet,   only : nmatsz, iba, kg1, kg_gamma, nbase_gamma, &
       &                            m_pwBS_kinetic_energies
  use m_Timing,              only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters,  only : kimg, neg, ipri
  use m_Const_Parameters,    only : DP, GAMMA, GAMMA_base_symmetrization
  use m_Kpoints,             only : k_symmetry, vkxyz
  use m_Parallelization,     only : mype,map_z,map_ek,myrank_k,map_k, ista_e, iend_e, istep_e   &
 &                                 , ista_k, iend_k , ista_g1k, iend_g1k, np_g1k, nrank_g, np_e

! ============================ added by K. Tagami ======================== 11.0
  use m_Control_Parameters,    only : noncol
  use m_ES_nonlocal,           only : m_ES_PAO_WFs_noncl
! ========================================================================= 11.0


  implicit none

!  58. m_ESIW_by_randomnumbers    <-(Initial_Electronic_Structure)
contains

  subroutine m_ESIW_by_randomnumbers(nfout,kv3,ie_start, ie_end, skip)   !(rndzaj)
    integer, intent(in) :: nfout,kv3,ie_start,ie_end
    logical, intent(in), optional :: skip

    real(kind=DP), allocatable, dimension(:) :: zaj_t

#ifdef _DEBUG_INITIALIZATION_
    real(kind=DP)   :: a,b,p,xn
    data a,b,p/32771.d0,1234567891.d0,2147483648.d0/
#endif

    integer, allocatable, dimension(:) :: vlambda, vc, nlambda, nc, mlambda, mc

    integer :: iimg, ik, ieg, i, iegt, xp, ig1, ig2
    integer(kind=8) :: jump
    logical :: first, set_ordr

    integer :: id_sname = -1
    call tstatc0_begin('m_ESIW_by_randomnumbers ',id_sname)

    if(ipri >= 2) then
       write(nfout,*)
       write(nfout,'(" <<< m_ESIW_by_randomnumbers >>>")')
       write(nfout,'(" !iba(1)        = ",i6)') iba(1)
    end if

    set_ordr = .true.
    if ( present(skip) ) then
       if ( skip ) set_ordr = .false.
    endif

    do ik = ista_k, iend_k  ! mpi
       if(iend_k > kv3) cycle
       do ieg = ista_e,iend_e,istep_e  ! mpi
!          if(ista_e < ie_start) cycle
!          if(iend_e > ie_end)   cycle
          if (ieg < ie_start) cycle
          if (ieg > ie_end)   cycle
          iegt = map_z(ieg)
          if(kimg == 1) then
             do i = 1, kg1
                zaj_l(i,iegt,ik,1) = 0.d0
             end do
          else
             do i = 1, kg1
                zaj_l(i,iegt,ik,1) = 0.d0
                zaj_l(i,iegt,ik,2) = 0.d0
             end do
          end if
       end do
    end do

    allocate(vlambda(0:1023), vc(0:1023), nlambda(0:1023), nc(0:1023), mlambda(0:1023), mc(0:1023))
    first = .true.

#ifdef _DEBUG_NO_INVERSION_
    iimg = 1
#else
    do iimg = 1, kimg
#endif
#ifdef _DEBUG_INITIALIZATION_
       do ik = 1, kv3
          do ieg = ie_start, ie_end
             do i = 1,  iba(ik)
                xn = mod(xn*a+b,p)
                if(map_ek(ieg,ik) == mype) zaj_l(i,map_z(ieg),ik,iimg) = xn/p ! MPI
             enddo
          enddo
       enddo
#else
       do ik = ista_k, iend_k              ! mpi
          if(iend_k > kv3) cycle
          if(k_symmetry(ik) == GAMMA .or. k_symmetry(ik) == GAMMA_base_symmetrization) then
             allocate(zaj_t(kg_gamma))
             do ieg = ista_e,iend_e,istep_e   ! mpi
!!$                if(ista_e < ie_start) cycle
!!$                if(iend_e > ie_end)   cycle
                if(ieg < ie_start) cycle
                if(ieg > ie_end)   cycle
                xp = 0
                jump = kg_gamma*((ie_end-ie_start+1)*(kv3*(iimg-1_8) + (ik-1_8)) + (ieg-ie_start))
                call pran0(xp,jump,kg_gamma,zaj_t)
                if(k_symmetry(ik) == GAMMA) then
                   zaj_l(1:kg_gamma,map_z(ieg),ik,iimg) = zaj_t(1:kg_gamma)
                else if(k_symmetry(ik) == GAMMA_base_symmetrization) then
                   if(iimg == 1) then
                      do i = 2, kg_gamma
                         ig1 = nbase_gamma(i,1)
                         ig2 = nbase_gamma(i,2)
                         zaj_l(ig1,map_z(ieg),ik,iimg) = zaj_t(i)
                         zaj_l(ig2,map_z(ieg),ik,iimg) = zaj_t(i)
                      end do
                      zaj_l(1,map_z(ieg),ik,iimg) = zaj_t(1)
                   else
                      do i = 2, kg_gamma
                         ig1 = nbase_gamma(i,1)
                         ig2 = nbase_gamma(i,2)
                         zaj_l(ig1,map_z(ieg),ik,iimg) =  zaj_t(i)
                         zaj_l(ig2,map_z(ieg),ik,iimg) = -zaj_t(i)
                      end do
                      zaj_l(1,map_z(ieg),ik,iimg) = 0.d0
                   end if
                end if
                first = .false.
             end do
             deallocate(zaj_t)
          else
             do ieg = ista_e,iend_e,istep_e   ! mpi
                if(ieg < ie_start) cycle
                if(ieg > ie_end)   cycle
!!$             xn = 0.d0
!!$             jump = kg1*((ie_end-ie_start)*(kv3*(iimg-1) + (ik-1)) + (ieg-ie_start))
!!$             do i = 1, jump
!!$                xn = mod(xn*a+b,p)
!!$             end do
!!$             do i = 1, iba(ik)
!!$                xn = mod(xn*a+b,p)
!!$                zaj_l(i,map_z(ieg),ik,iimg) = xn/p ! MPI
!!$             end do
                xp = 0
                jump = kg1*((ie_end-ie_start+1)*(kv3*(iimg-1_8) + (ik-1_8)) + (ieg-ie_start))
                call pran0(xp,jump,iba(ik),zaj_l(1,map_z(ieg),ik,iimg))
                first = .false.
             end do
          end if
       end do
#endif
#ifndef _DEBUG_NO_INVERSION_
    enddo
#endif
    deallocate(mc, mlambda, nc, nlambda, vc, vlambda)

    if(ipri >= 2) write(nfout,'(" initialization of neordr and nrvf_ordr")')
! Initialization

    if ( set_ordr ) then
       do ik = 1, kv3
          if(map_k(ik) /= myrank_k) cycle ! MPI
!!$       do i = 1, neg
          do i = ie_start, ie_end
             neordr(i,ik)    = i
             nrvf_ordr(i,ik) = i
          end do
       end do
    endif
    call tstatc0_end(id_sname)

  contains
    subroutine pran0(iseed,nstep,n1,r)
!
!  AUTHOR: T.Kokubo   Jun/15/2005
!
      implicit none

      integer :: iseed,n1
      integer(kind=8) :: nstep
      real(kind=DP) :: r(n1)

      integer, PARAMETER :: lambda=32771,c=1234567891,m=2147483647
      real(kind=DP), PARAMETER :: Minv=1.0D+00/2147483648.d0

      integer :: k,j
      integer(kind=8) :: istep,jstep,ijmp,jjmp,ijmpm,jjmpm

      if(first) then
         vlambda(0)=lambda
         vc(0)=c
#ifdef SX
!cdir novector
#endif
#ifdef HIUX
*poption noparallel
#endif
         do k=1,1023
            vlambda(k)=lambda*vlambda(k-1)
            vc(k)=vc(k-1)+c*vlambda(k-1)
         enddo
         nlambda(0)=vlambda(1023)
         nc(0)=vc(1023)
#ifdef SX
!cdir novector
#endif
#ifdef HIUX
*poption noparallel
#endif
         do k=1,1023
            nlambda(k)=vlambda(1023)*nlambda(k-1)
            nc(k)=nc(k-1)+vc(1023)*nlambda(k-1)
         enddo
         mlambda(0)=nlambda(1023)
         mc(0)=nc(1023)
#ifdef SX
!cdir novector
#endif
#ifdef HIUX
*poption noparallel
#endif
         do k=1,1023
            mlambda(k)=nlambda(1023)*mlambda(k-1)
            mc(k)=mc(k-1)+nc(1023)*mlambda(k-1)
         enddo
!!$         first=.false.
      endif

      istep=nstep/1024_8
      jstep=nstep-istep*1024_8

      ijmp=istep/1024_8
      jjmp=istep-ijmp*1024_8

      ijmpm=ijmp/1024_8
      jjmpm=ijmp-ijmpm*1024_8

#ifdef SX
!cdir novector
#endif
#ifdef HIUX
*poption noparallel
#endif
      do j=1,ijmpm
         iseed=iseed*mlambda(1023)+mc(1023)
      enddo
      if(jjmpm.ne.0) then
         iseed=iseed*mlambda(jjmpm-1)+mc(jjmpm-1)
      endif
      if(jjmp.ne.0) then
         iseed=iseed*nlambda(jjmp-1)+nc(jjmp-1)
      endif
      if(jstep.ne.0) then
         iseed=iseed*vlambda(jstep-1)+vc(jstep-1)
      endif

      do j=1,n1,1024
#ifdef SX
!cdir nodep
#endif
#ifdef HIUX
*poption indep
#endif
         do k=0,min(1023,n1-j)
            r(j+k)=minv*iand(iseed*vlambda(k)+vc(k),m)
         enddo
         iseed=iseed*vlambda(k-1)+vc(k-1)
      enddo
    end subroutine pran0
  end subroutine m_ESIW_by_randomnumbers

  subroutine m_ESIW_mul_by_randomnumbers( nfout, kv3, ie_start, ie_end, skip )
    use m_Control_Parameters,  only : noise_amplitude, noise_mode, gmax
                                                                     !(rndzaj)
    integer, intent(in) :: nfout,kv3,ie_start,ie_end
    logical, intent(in), optional :: skip

    real(kind=DP), allocatable, dimension(:) :: zaj_t,  zaj_wk

#ifdef _DEBUG_INITIALIZATION_
    real(kind=DP)   :: a,b,p,xn
    data a,b,p/32771.d0,1234567891.d0,2147483648.d0/
#endif

    integer, allocatable, dimension(:) :: vlambda, vc, nlambda, nc, mlambda, mc

    integer :: iimg, ik, ieg, i, iegt, xp, ig1, ig2
    integer(kind=8) :: jump
    real(kind=DP) :: fac1, c1, c2
    real(kind=DP) :: sigma2
    real(kind=DP), allocatable :: ekin(:), filter(:)

    logical :: first

    integer :: id_sname = -1
    call tstatc0_begin('m_ESIW_mul_by_randomnumbers ',id_sname)

    fac1 = noise_amplitude
    if ( noise_mode == 0 ) then
!       fac1 = fac1 /dble(kg1) *0.10d0;    sigma2 = gmax**2
!       fac1 = fac1 /dble(kg1) /2.0;    sigma2 = gmax**2 *2.0
       fac1 = fac1 /dble(kg1);    sigma2 = gmax**2 *2.0
!       fac1 = fac1 /dble(kg1);    sigma2 = gmax**2 
    endif

    if(ipri >= 2) then
       write(nfout,*)
       write(nfout,'(" <<< m_ESIW_mul_by_randomnumbers >>>")')
       write(nfout,'(" !iba(1)        = ",i6)') iba(1)
    end if

    allocate(vlambda(0:1023), vc(0:1023), nlambda(0:1023), nc(0:1023), mlambda(0:1023), mc(0:1023))
    first = .true.

    do iimg = 1, kimg

       do ik = ista_k, iend_k              ! mpi
          if(iend_k > kv3) cycle
          if(k_symmetry(ik) == GAMMA .or. k_symmetry(ik) == GAMMA_base_symmetrization) then
             allocate(zaj_t(kg_gamma))
             if ( noise_mode == 0 ) then
                allocate( ekin(kg1) );  allocate( filter(kg1) )
                call m_pwBS_kinetic_energies(ik,vkxyz,ekin)
                filter(1:iba(ik)) = exp( -ekin(1:iba(ik))/sigma2 )
!                filter = 1.0d0
             endif

             do ieg = ista_e,iend_e,istep_e   ! mpi
                if(ieg < ie_start) cycle
                if(ieg > ie_end)   cycle
                xp = 0
                jump = kg_gamma*((ie_end-ie_start+1)*(kv3*(iimg-1_8) + (ik-1_8)) + (ieg-ie_start))
                call pran0(xp,jump,kg_gamma,zaj_t)

                if (k_symmetry(ik) == GAMMA) then
                   do ig1 = 1, kg_gamma
                      if ( noise_mode == 0 ) then
                         c1 = fac1 *zaj_t(ig1)
                         zaj_l(ig1,map_z(ieg),ik,iimg) &
                              &     = zaj_l(ig1,map_z(ieg),ik,iimg) *filter(ig1) +c1
!                         c1 = fac1 *zaj_t(ig1) *filter(ig1)
!                         zaj_l(ig1,map_z(ieg),ik,iimg) &
!                              &     = zaj_l(ig1,map_z(ieg),ik,iimg) +c1
                      else if ( noise_mode == 1 ) then
                         c1 = 1.0d0 +fac1 *zaj_t(ig1)
                         zaj_l(ig1,map_z(ieg),ik,iimg) &
                              &     = zaj_l(ig1,map_z(ieg),ik,iimg) *c1
                      endif
                   end do
                else if(k_symmetry(ik) == GAMMA_base_symmetrization) then
                   if(iimg == 1) then
                      do i = 2, kg_gamma
                         ig1 = nbase_gamma(i,1);  ig2 = nbase_gamma(i,2)
                         if ( noise_mode == 0 ) then
                            c1 = fac1 *zaj_t(i)
                            zaj_l(ig1,map_z(ieg),ik,iimg) &
                                 &   = zaj_l(ig1,map_z(ieg),ik,iimg) *filter(i) +c1
                            zaj_l(ig2,map_z(ieg),ik,iimg) & 
                                 &   = zaj_l(ig2,map_z(ieg),ik,iimg) *filter(i) +c1
                         else if ( noise_mode == 1 ) then
                            c1 = 1.0d0 +fac1 *zaj_t(i)
                            zaj_l(ig1,map_z(ieg),ik,iimg) &
                                 &   = zaj_l(ig1,map_z(ieg),ik,iimg) *c1
                            zaj_l(ig2,map_z(ieg),ik,iimg) & 
                                 &   = zaj_l(ig2,map_z(ieg),ik,iimg) *c1
                         endif
                      end do
                      if ( noise_mode == 0 ) then
                         c1 = fac1 *zaj_t(1)
                         zaj_l(1,map_z(ieg),ik,iimg) &
                              &    = zaj_l(1,map_z(ieg),ik,iimg) *filter(1) +c1
                      else if ( noise_mode == 1 ) then
                         c1 = 1.0d0 +fac1 *zaj_t(1)
                         zaj_l(1,map_z(ieg),ik,iimg) &
                              &    = zaj_l(1,map_z(ieg),ik,iimg) *c1
                      endif

                   else
                      do i = 2, kg_gamma
                         ig1 = nbase_gamma(i,1);  ig2 = nbase_gamma(i,2)
                         if ( noise_mode == 0 ) then
                            c1 = fac1 *zaj_t(i)
                            c2 = -fac1 *zaj_t(i)
                            zaj_l(ig1,map_z(ieg),ik,iimg) &
                                 &   = zaj_l(ig1,map_z(ieg),ik,iimg) *filter(i) +c1
                            zaj_l(ig2,map_z(ieg),ik,iimg) &
                                 &   = zaj_l(ig2,map_z(ieg),ik,iimg) *filter(i) +c2
                         else if ( noise_mode == 1 ) then
                            c1 = 1.0d0 +fac1 *zaj_t(i)
                            c2 = 1.0d0 -fac1 *zaj_t(i)
                            zaj_l(ig1,map_z(ieg),ik,iimg) &
                                 &   = zaj_l(ig1,map_z(ieg),ik,iimg) *c1
                            zaj_l(ig2,map_z(ieg),ik,iimg) &
                                 &   = zaj_l(ig2,map_z(ieg),ik,iimg) *c2
                         endif
                      end do
                      zaj_l(1,map_z(ieg),ik,iimg) = 0.0d0
                   end if
                endif
                first = .false.
             end do
             deallocate(zaj_t)
             if ( allocated( ekin ) ) deallocate( ekin )
             if ( allocated( filter ) ) deallocate( filter )

          else
             allocate( zaj_t( iba(ik) ) )
             if ( noise_mode == 0 ) then
                allocate( ekin(kg1) );  allocate( filter(kg1) )
                call m_pwBS_kinetic_energies(ik,vkxyz,ekin)
                filter(1:iba(ik)) = exp( -ekin(1:iba(ik))/sigma2 )
             endif

             do ieg = ista_e,iend_e,istep_e   ! mpi
                if(ieg < ie_start) cycle
                if(ieg > ie_end)   cycle
                xp = 0
                jump = kg1*((ie_end-ie_start+1)*(kv3*(iimg-1_8) &
                     &                       + (ik-1_8)) + (ieg-ie_start))
                call pran0(xp,jump,iba(ik),zaj_t)

                if ( iimg == 1 ) then
                   do ig1 = 1, iba(ik)
                      if ( noise_mode == 0 ) then
                         c1 = fac1 *zaj_t(ig1)
                         zaj_l(ig1,map_z(ieg),ik,iimg) &
                              &      = zaj_l(ig1,map_z(ieg),ik,iimg) *filter(ig1) +c1
                      else if ( noise_mode == 1 ) then
                         c1 = 1.0d0 +fac1 *zaj_t(ig1)
                         zaj_l(ig1,map_z(ieg),ik,iimg) &
                              &      = zaj_l(ig1,map_z(ieg),ik,iimg) *c1
                      endif
                   end do
                else
                   do ig1 = 1, iba(ik)
                      if ( noise_mode == 0 ) then
                         c1 = fac1 *zaj_t(ig1)
                         zaj_l(ig1,map_z(ieg),ik,iimg) &
                              &       = zaj_l(ig1,map_z(ieg),ik,iimg) *filter(ig1) +c1
                      else if ( noise_mode == 1 ) then
                         c1 = 1.0d0 +fac1 *zaj_t(ig1)
                         zaj_l(ig1,map_z(ieg),ik,iimg) &
                              &       = zaj_l(ig1,map_z(ieg),ik,iimg) *c1
                      endif
                   end do
                endif
                first = .false.
             end do
             deallocate( zaj_t )
             if ( allocated( ekin ) ) deallocate( ekin )
             if ( allocated( filter ) ) deallocate( filter )
          endif
       end do
    enddo

    deallocate(mc, mlambda, nc, nlambda, vc, vlambda)

    if(ipri >= 2) write(nfout,'(" initialization of neordr and nrvf_ordr")')
! Initialization

    if ( ( .not. present(skip) ) .or. ( .not. skip ) ) then
       do ik = 1, kv3
          if(map_k(ik) /= myrank_k) cycle ! MPI
!!$       do i = 1, neg
          do i = ie_start, ie_end
             neordr(i,ik)    = i
             nrvf_ordr(i,ik) = i
          end do
       end do
    endif

    call tstatc0_end(id_sname)

  contains
    subroutine pran0(iseed,nstep,n1,r)
!
!  AUTHOR: T.Kokubo   Jun/15/2005
!
      implicit none

      integer :: iseed,n1
      integer(kind=8) :: nstep
      real(kind=DP) :: r(n1)

      integer, PARAMETER :: lambda=32771,c=1234567891,m=2147483647
      real(kind=DP), PARAMETER :: Minv=1.0D+00/2147483648.d0

      integer :: k,j
      integer(kind=8) :: istep,jstep,ijmp,jjmp,ijmpm,jjmpm

      if(first) then
         vlambda(0)=lambda
         vc(0)=c
#ifdef SX
!cdir novector
#endif
#ifdef HIUX
*poption noparallel
#endif
         do k=1,1023
            vlambda(k)=lambda*vlambda(k-1)
            vc(k)=vc(k-1)+c*vlambda(k-1)
         enddo
         nlambda(0)=vlambda(1023)
         nc(0)=vc(1023)
#ifdef SX
!cdir novector
#endif
#ifdef HIUX
*poption noparallel
#endif
         do k=1,1023
            nlambda(k)=vlambda(1023)*nlambda(k-1)
            nc(k)=nc(k-1)+vc(1023)*nlambda(k-1)
         enddo
         mlambda(0)=nlambda(1023)
         mc(0)=nc(1023)
#ifdef SX
!cdir novector
#endif
#ifdef HIUX
*poption noparallel
#endif
         do k=1,1023
            mlambda(k)=nlambda(1023)*mlambda(k-1)
            mc(k)=mc(k-1)+nc(1023)*mlambda(k-1)
         enddo
!!$         first=.false.
      endif

      istep=nstep/1024_8
      jstep=nstep-istep*1024_8

      ijmp=istep/1024_8
      jjmp=istep-ijmp*1024_8

      ijmpm=ijmp/1024_8
      jjmpm=ijmp-ijmpm*1024_8

#ifdef SX
!cdir novector
#endif
#ifdef HIUX
*poption noparallel
#endif
      do j=1,ijmpm
         iseed=iseed*mlambda(1023)+mc(1023)
      enddo
      if(jjmpm.ne.0) then
         iseed=iseed*mlambda(jjmpm-1)+mc(jjmpm-1)
      endif
      if(jjmp.ne.0) then
         iseed=iseed*nlambda(jjmp-1)+nc(jjmp-1)
      endif
      if(jstep.ne.0) then
         iseed=iseed*vlambda(jstep-1)+vc(jstep-1)
      endif

      do j=1,n1,1024
#ifdef SX
!cdir nodep
#endif
#ifdef HIUX
*poption indep
#endif
         do k=0,min(1023,n1-j)
            r(j+k)=minv*iand(iseed*vlambda(k)+vc(k),m)
         enddo
         iseed=iseed*vlambda(k-1)+vc(k-1)
      enddo
    end subroutine pran0
  end subroutine m_ESIW_mul_by_randomnumbers

  subroutine m_ESIW_by_atomic_orbitals(nfout,kv3,ie_start,ie_end) !(paozaj)
    use m_Control_Parameters,  only : noise_amplitude
    use m_Const_Parameters,    only : YES
    use m_PseudoPotential,     only : num_wfc_pao

    integer, intent(in) :: nfout,kv3,ie_start,ie_end

    integer :: ik, i

    call m_ESIW_by_randomnumbers( nfout, kv3, ie_start, ie_end, skip=.true. )

    do ik = 1, kv3
       if (map_k(ik) /= myrank_k) cycle ! MPI
       do i = ie_start, ie_end
          neordr(i,ik)    = i
          nrvf_ordr(i,ik) = i
       end do
    end do

    if ( noncol ) then
       call m_ES_PAO_WFs_noncl(nfout)
    else
       call m_ES_PAO_WFs(nfout)
    endif

    if ( noise_amplitude > 0.0 ) then
       call m_ESIW_mul_by_randomnumbers( nfout, kv3, ie_start, num_wfc_pao, skip=.true. )
    endif

  end subroutine m_ESIW_by_atomic_orbitals
!===============================================================================
end module m_ES_initialWF
