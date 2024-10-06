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
  use m_ES_nonlocal         ,only : m_ES_PAO_WFs_3D
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
!  use m_ES_nonlocal,           only : m_ES_PAO_WFs_noncl
! ========================================================================= 11.0

  implicit none

!  58. m_ESIW_by_randomnumbers    <-(Initial_Electronic_Structure)
contains

  subroutine m_ESIW_by_randomnumbers_3D(nfout,kv3,ie_start, ie_end, isk, iek, skip)   !(rndzaj)
!
!f
    use m_Parallelization,     only :  nbsn_sta, nbsn_end, neg_g, nbsn_num
    integer :: l_bl, n_sta, n_end
    integer, intent(in) :: nfout,kv3,ie_start,ie_end
    integer, intent(in), optional :: isk, iek
    logical, intent(in), optional :: skip

    real(kind=DP), allocatable, dimension(:) :: zaj_t
    integer :: is, ie

#ifdef _DEBUG_INITIALIZATION_
    real(kind=DP)   :: a,b,p,xn
    data a,b,p/32771.d0,1234567891.d0,2147483648.d0/
#endif

    integer, allocatable, dimension(:) :: vlambda, vc, nlambda, nc, mlambda, mc

    integer :: iimg, ik, ieg, i, iegt, xp, ig1, ig2
    integer(kind=8) :: jump
    logical :: first, set_ordr

    integer :: id_sname = -1
    write(nfout,'(" <<< m_ESIW_by_randomnumbers_3D >>>")')
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

!3D
    do ik = ista_k, iend_k  ! mpi
       if(iend_k > kv3) cycle
       do ieg = 1, np_e
          if ( neg_g(ieg) < ie_start ) cycle
          if ( neg_g(ieg) > ie_end )   cycle
          if (kimg == 1) then
             do i = 1, np_g1k(ik)
                zaj_l(i,ieg,ik,1) = 0.d0
             end do
          else
             do i = 1, np_g1k(ik)
                zaj_l(i,ieg,ik,1) = 0.d0
                zaj_l(i,ieg,ik,2) = 0.d0
             end do
          end if
       end do
    end do

    allocate(vlambda(0:1023), vc(0:1023), nlambda(0:1023), nc(0:1023), mlambda(0:1023), mc(0:1023))
    first = .true.

    is = ista_k
    ie = iend_k
    if(present(isk)) is = isk
    if(present(iek)) ie = iek
    do iimg = 1, kimg

       do ik = is, ie              ! mpi
          if(ie > kv3) cycle
          if(k_symmetry(ik) == GAMMA .or. k_symmetry(ik) == GAMMA_base_symmetrization) then
!f             allocate(zaj_t(kg_gamma))
             allocate(zaj_t(np_g1k(ik)))
!modfy for block
#ifdef __TIMER_INIDO__
  call timer_sta(1345)
#endif
             do ieg = 1,np_e
!f             do ieg = ista_e,iend_e,istep_e   ! mpi
!!$                if(ista_e < ie_start) cycle
!!$                if(iend_e > ie_end)   cycle
                if ( neg_g(ieg) < ie_start) cycle
                if ( neg_g(ieg) > ie_end )   cycle
                xp = 0
!f                jump = kg_gamma*((ie_end-ie_start+1)*(kv3*(iimg-1_8) + (ik-1_8)) + (ieg-ie_start))
                jump = kg_gamma*((ie_end-ie_start+1)*(kv3*(iimg-1_8) + (ik-1_8)) + (neg_g(ieg)-ie_start)) + ista_g1k(ik) -1
!f                call pran0(xp,jump,kg_gamma,zaj_t)
                call pran0(xp,jump,np_g1k(ik),zaj_t)
                if(k_symmetry(ik) == GAMMA) then
!f                   zaj_l(1:kg_gamma,map_z(ieg),ik,iimg) = zaj_t(1:kg_gamma)
                   zaj_l(1:np_g1k(ik),ieg,ik,iimg) = zaj_t(1:np_g1k(ik))
                else if(k_symmetry(ik) == GAMMA_base_symmetrization) then
                   if(iimg == 1) then
!f                      do i = 2, kg_gamma
!f                         ig1 = nbase_gamma(i,1)
!f                         ig2 = nbase_gamma(i,2)
!f                         zaj_l(ig1,ieg,ik,iimg) = zaj_t(i)
!f                         zaj_l(ig2,ieg,ik,iimg) = zaj_t(i)
!f                      end do
!f                      zaj_l(1,ieg,ik,iimg) = zaj_t(1)
                      zaj_l(1:np_g1k(ik),ieg,ik,iimg) = zaj_t(1:np_g1k(ik))
                   else
!f                      do i = 2, kg_gamma
!f                         ig1 = nbase_gamma(i,1)
!f                         ig2 = nbase_gamma(i,2)
!f                         zaj_l(ig1,ieg,ik,iimg) =  zaj_t(i)
!f                         zaj_l(ig2,ieg,ik,iimg) = -zaj_t(i)
!f                      end do
!f                      zaj_l(1,ieg,ik,iimg) = 0.d0
                      zaj_l(1:np_g1k(ik),ieg,ik,iimg) = zaj_t(1:np_g1k(ik))
                   end if
                end if
                first = .false.
             end do
#ifdef __TIMER_INIDO__
  call timer_end(1345)
#endif
             deallocate(zaj_t)
          else
#ifdef __TIMER_INIDO__
  call timer_sta(1346)
#endif
             do ieg = 1,np_e
!f             do ieg = ista_e,iend_e,istep_e   ! mpi
                if( neg_g(ieg) < ie_start) cycle
                if( neg_g(ieg) > ie_end)   cycle

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
!f                jump = kg1*((ie_end-ie_start+1)*(kv3*(iimg-1_8) + (ik-1_8)) + (ieg-ie_start))
!f                call pran0(xp,jump,iba(ik),zaj_l(1,map_z(ieg),ik,iimg))
                jump = kg1*((ie_end-ie_start+1)*(kv3*(iimg-1_8) + (ik-1_8)) + (neg_g(ieg)-ie_start)) + ista_g1k(ik) -1
                call pran0(xp,jump,np_g1k(ik),zaj_l(1,ieg,ik,iimg))
                first = .false.
             end do
#ifdef __TIMER_INIDO__
  call timer_end(1346)
#endif
          end if
       end do
!

    enddo

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
  end subroutine m_ESIW_by_randomnumbers_3D
!----
  subroutine m_ESIW_by_randomnumbers0_3D(nfout,kv3,ie_start, ie_end, isk, iek)   !(rndzaj)
!3D
!f
    use m_Parallelization,     only :  nbsn_sta, nbsn_end, neg_g, nbsn_num
    integer :: l_bl, n_sta, n_end
    integer, intent(in) :: nfout,kv3,ie_start,ie_end
    integer, intent(in), optional :: isk, iek
    real(kind=DP), allocatable, dimension(:) :: zaj_t
    real(kind=DP), allocatable, dimension(:) :: zaj_wk

#ifdef _DEBUG_INITIALIZATION_
    real(kind=DP)   :: a,b,p,xn
    data a,b,p/32771.d0,1234567891.d0,2147483648.d0/
#endif

    integer, allocatable, dimension(:) :: vlambda, vc, nlambda, nc, mlambda, mc

    integer :: iimg, ik, ieg, i, iegt, xp, ig1, ig2
    integer(kind=8) :: jump
    logical :: first
    integer :: is, ie

    integer :: id_sname = -1
    write(nfout,'(" <<< m_ESIW_by_randomnumbers0_3D >>>")')
    call tstatc0_begin('m_ESIW_by_randomnumbers ',id_sname)

    if(ipri >= 2) then
       write(nfout,*)
       write(nfout,'(" <<< m_ESIW_by_randomnumbers >>>")')
       write(nfout,'(" !iba(1)        = ",i6)') iba(1)
    end if
!3D
   allocate(zaj_wk(kg1))
   zaj_l = 0.0d0

    allocate(vlambda(0:1023), vc(0:1023), nlambda(0:1023), nc(0:1023), mlambda(0:1023), mc(0:1023))
    first = .true.
    is = ista_k
    ie = iend_k
    if(present(isk)) is = isk
    if(present(iek)) ie = iek
    do iimg = 1, kimg

       do ik = is, ie              ! mpi
          if(ie > kv3) cycle
          if(k_symmetry(ik) == GAMMA .or. k_symmetry(ik) == GAMMA_base_symmetrization) then
             allocate(zaj_t(kg_gamma))
!modfy for block
             do ieg = 1,np_e
!f             do ieg = ista_e,iend_e,istep_e   ! mpi
!!$                if(ista_e < ie_start) cycle
!!$                if(iend_e > ie_end)   cycle
                if(ieg < ie_start) cycle
                if(ieg > ie_end)   cycle
                xp = 0
                zaj_wk = 0.0d0
!f                jump = kg_gamma*((ie_end-ie_start+1)*(kv3*(iimg-1_8) + (ik-1_8)) + (ieg-ie_start))
                jump = kg_gamma*((ie_end-ie_start+1)*(kv3*(iimg-1_8) + (ik-1_8)) + (neg_g(ieg)-ie_start))
                call pran0(xp,jump,kg_gamma,zaj_t)
                if(k_symmetry(ik) == GAMMA) then
!f                   zaj_l(1:kg_gamma,map_z(ieg),ik,iimg) = zaj_t(1:kg_gamma)
                   zaj_wk(1:kg_gamma) = zaj_t(1:kg_gamma)
                else if(k_symmetry(ik) == GAMMA_base_symmetrization) then
                   if(iimg == 1) then
                      do i = 2, kg_gamma
                         ig1 = nbase_gamma(i,1)
                         ig2 = nbase_gamma(i,2)
                         zaj_wk(ig1) = zaj_t(i)
                         zaj_wk(ig2) = zaj_t(i)
                      end do
                      zaj_wk(1) = zaj_t(1)
                   else
                      do i = 2, kg_gamma
                         ig1 = nbase_gamma(i,1)
                         ig2 = nbase_gamma(i,2)
                         zaj_wk(ig1) =  zaj_t(i)
                         zaj_wk(ig2) = -zaj_t(i)
                      end do
                      zaj_wk(1) = 0.d0
                   end if
                end if
                first = .false.
                do i = ista_g1k(ik),iend_g1k(ik)
                   zaj_l(i-ista_g1k(ik)+1,ieg,ik,iimg) = zaj_wk(i)
                end do
             end do
             deallocate(zaj_t)
          else
             do ieg = 1,np_e
!f             do ieg = ista_e,iend_e,istep_e   ! mpi
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
!f                jump = kg1*((ie_end-ie_start+1)*(kv3*(iimg-1_8) + (ik-1_8)) + (ieg-ie_start))
!f                call pran0(xp,jump,iba(ik),zaj_l(1,map_z(ieg),ik,iimg))
                jump = kg1*((ie_end-ie_start+1)*(kv3*(iimg-1_8) + (ik-1_8)) + (neg_g(ieg)-ie_start))
                zaj_wk = 0.0d0
                call pran0(xp,jump,iba(ik),zaj_wk)
                do i=ista_g1k(ik),iend_g1k(ik)
                   zaj_l(i-ista_g1k(ik)+1,ieg,ik,iimg) = zaj_wk(i)
                enddo
                first = .false.
             end do
          end if
       end do
!

    enddo

    deallocate(mc, mlambda, nc, nlambda, vc, vlambda)
    deallocate(zaj_wk)

    if(ipri >= 2) write(nfout,'(" initialization of neordr and nrvf_ordr")')
! Initialization

    do ik = 1, kv3
       if(map_k(ik) /= myrank_k) cycle ! MPI
!!$       do i = 1, neg
       do i = ie_start, ie_end
          neordr(i,ik)    = i
          nrvf_ordr(i,ik) = i
       end do
    end do

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
  end subroutine m_ESIW_by_randomnumbers0_3D
!$$#endif

  subroutine m_ESIW_mul_by_randomnumbers_3D( nfout, kv3, ie_start, ie_end, &
       &                                     isk, iek, skip )  !(rndzaj)
    use m_Parallelization,     only : nbsn_sta, nbsn_end, neg_g, nbsn_num
    use m_Control_Parameters,  only : noise_amplitude, noise_mode, gmax

    integer :: l_bl, n_sta, n_end
    integer, intent(in) :: nfout,kv3,ie_start,ie_end

    integer, intent(in), optional :: isk, iek
    logical, intent(in), optional :: skip

    real(kind=DP), allocatable, dimension(:) :: zaj_t
    integer :: is, ie

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

    write(nfout,'(" <<< m_ESIW_mul_by_randomnumbers_3D >>>")')
    call tstatc0_begin('m_ESIW_mul_by_randomnumbers ',id_sname)

    fac1 = noise_amplitude
    if ( noise_mode == 0 ) then
       fac1 = fac1 /dble(kg1);    sigma2 = gmax**2 *2.0
    endif

    if(ipri >= 2) then
       write(nfout,*)
       write(nfout,'(" <<< m_ESIW_by_randomnumbers >>>")')
       write(nfout,'(" !iba(1)        = ",i6)') iba(1)
    end if
!3D

    allocate(vlambda(0:1023), vc(0:1023), nlambda(0:1023), nc(0:1023), mlambda(0:1023), mc(0:1023))
    first = .true.

    is = ista_k
    ie = iend_k
    if(present(isk)) is = isk
    if(present(iek)) ie = iek

    do iimg = 1, kimg

       do ik = is, ie              ! mpi
          if(ie > kv3) cycle
          if(k_symmetry(ik) == GAMMA .or. k_symmetry(ik) == GAMMA_base_symmetrization) then
!f             allocate(zaj_t(kg_gamma))
             allocate(zaj_t(np_g1k(ik)))
!modfy for block
#ifdef __TIMER_INIDO__
  call timer_sta(1345)
#endif
             if ( noise_mode == 0 ) then
                allocate( ekin(np_g1k(ik)) );  allocate( filter(np_g1k(ik)) )
                call m_pwBS_kinetic_energies(ik,vkxyz,ekin)
                filter(1:np_g1k(ik)) = exp( -ekin(1:np_g1k(ik))/sigma2 )
             endif

             do ieg = 1,np_e
!f             do ieg = ista_e,iend_e,istep_e   ! mpi
!!$                if(ista_e < ie_start) cycle
!!$                if(iend_e > ie_end)   cycle
                if ( neg_g(ieg) < ie_start) cycle
                if ( neg_g(ieg) > ie_end)   cycle
                xp = 0
!f                jump = kg_gamma*((ie_end-ie_start+1)*(kv3*(iimg-1_8) + (ik-1_8)) + (ieg-ie_start))
                jump = kg_gamma*((ie_end-ie_start+1)*(kv3*(iimg-1_8) + (ik-1_8)) + (neg_g(ieg)-ie_start)) + ista_g1k(ik) -1
!f                call pran0(xp,jump,kg_gamma,zaj_t)
                call pran0(xp,jump,np_g1k(ik),zaj_t)

                if(k_symmetry(ik) == GAMMA) then
                   if ( noise_mode == 0 ) then
                      zaj_l(1:np_g1k(ik),ieg,ik,iimg) &
                           &   = zaj_l(1:np_g1k(ik),ieg,ik,iimg) *filter(1:np_g1k(ik)) &
                           &    +fac1 *zaj_t(1:np_g1k(ik))
                   else if ( noise_mode == 1 ) then
                      zaj_l(1:np_g1k(ik),ieg,ik,iimg) &
                           &   = zaj_l(1:np_g1k(ik),ieg,ik,iimg) &
                           &      *( 1.0d0 +fac1 *zaj_t(1:np_g1k(ik)) )
                   endif

                else if(k_symmetry(ik) == GAMMA_base_symmetrization) then
                   if ( noise_mode == 0 ) then
                      zaj_l(1:np_g1k(ik),ieg,ik,iimg) &
                           &   = zaj_l(1:np_g1k(ik),ieg,ik,iimg) *filter(1:np_g1k(ik)) &
                           &        +fac1 *zaj_t(1:np_g1k(ik)) 
                   else if ( noise_mode == 1 ) then
                      zaj_l(1:np_g1k(ik),ieg,ik,iimg) &
                           &   = zaj_l(1:np_g1k(ik),ieg,ik,iimg) &
                           &        *( 1.0d0 +fac1 *zaj_t(1:np_g1k(ik)) )
                   end if
                end if
                first = .false.
             end do
#ifdef __TIMER_INIDO__
  call timer_end(1345)
#endif
             deallocate(zaj_t)
             if ( allocated( ekin ) ) deallocate( ekin )
             if ( allocated( filter ) ) deallocate( filter )

          else
#ifdef __TIMER_INIDO__
  call timer_sta(1346)
#endif
             allocate(zaj_t(np_g1k(ik)))

             if ( noise_mode == 0 ) then
                allocate( ekin(np_g1k(ik)) );  allocate( filter(np_g1k(ik)) )
                call m_pwBS_kinetic_energies(ik,vkxyz,ekin)
                filter(1:np_g1k(ik)) = exp( -ekin(1:np_g1k(ik))/sigma2 )
             endif

             do ieg = 1,np_e
!f             do ieg = ista_e,iend_e,istep_e   ! mpi
                if ( neg_g(ieg) < ie_start) cycle
                if ( neg_g(ieg) > ie_end)   cycle
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
!f                jump = kg1*((ie_end-ie_start+1)*(kv3*(iimg-1_8) + (ik-1_8)) + (ieg-ie_start))
!f                call pran0(xp,jump,iba(ik),zaj_l(1,map_z(ieg),ik,iimg))
                jump = kg1*((ie_end-ie_start+1)*(kv3*(iimg-1_8) + (ik-1_8)) + (neg_g(ieg)-ie_start)) + ista_g1k(ik) -1
                call pran0(xp,jump,np_g1k(ik),zaj_t)

                if ( noise_mode == 0 ) then
                   zaj_l(1:np_g1k(ik),ieg,ik,iimg) &
                        &   = zaj_l(1:np_g1k(ik),ieg,ik,iimg) *filter(1:np_g1k(ik)) &
                        &        +fac1 *zaj_t(1:np_g1k(ik)) 
                else if ( noise_mode == 1 ) then
                   zaj_l(1:np_g1k(ik),ieg,ik,iimg) &
                        &   = zaj_l(1:np_g1k(ik),ieg,ik,iimg) &
                        &        *( 1.0d0 +fac1 *zaj_t(1:np_g1k(ik)) )
                endif

                first = .false.
             end do
             deallocate( zaj_t )
             if ( allocated( ekin ) ) deallocate( ekin )
             if ( allocated( filter ) ) deallocate( filter )

#ifdef __TIMER_INIDO__
  call timer_end(1346)
#endif
          end if
       end do
!

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
  end subroutine m_ESIW_mul_by_randomnumbers_3D

  subroutine m_ESIW_by_atomic_orbitals(nfout,kv3,ie_start,ie_end) !(paozaj)
    use m_Control_Parameters,  only : noise_amplitude
    use m_Const_Parameters,    only : YES
    use m_PseudoPotential,     only : num_wfc_pao

    integer, intent(in) :: nfout,kv3,ie_start,ie_end

    integer :: ik, i

    call m_ESIW_by_randomnumbers_3D( nfout, kv3, ie_start, ie_end, &
         &                           ista_k, iend_k, skip=.true. )

    do ik = 1, kv3
       if (map_k(ik) /= myrank_k) cycle ! MPI
       do i = ie_start, ie_end
          neordr(i,ik)    = i
          nrvf_ordr(i,ik) = i
       end do
    end do

    call m_ES_PAO_WFs_3D(nfout)

    if ( noise_amplitude > 0.0 ) then
       call m_ESIW_mul_by_randomnumbers_3D( nfout, kv3, ie_start, num_wfc_pao, &
            &                               ista_k, iend_k, skip=.true. )
    endif

  end subroutine m_ESIW_by_atomic_orbitals


end module m_ES_initialWF
