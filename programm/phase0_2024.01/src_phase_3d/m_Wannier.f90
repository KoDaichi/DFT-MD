!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  MODULE: m_Wannier
!
!  AUTHOR(S): T. Yamamoto   May/05/2007
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
module m_Wannier
!! Maximally localized Wannier functions
!! Ref.
!! [1] N. Marzari and D. Vanderbilt, Phys. Rev. B 56 (1997) 12847.
!! [2] G. Berghold, et.al., Phys. Rev. B 61 (2000) 10040.
  use m_Const_Parameters,   only : DP,PAI2,SD,CG,CGPRC,BINARY,CUBE,VTK,DENSITY_ONLY &
       &                         , SKIP,EXECUT,GAMMA,ON
  use m_Control_Parameters, only : printable, kimg, neg, eps_wan, dt_wan &
       &                         , wannier_opt_method, max_iter_wan &
       &                         , wannier_filetype, ipri, sw_random_wannier &
       &                         , sw_potential_wannier, sw_continue_wannier
  use m_Crystal_Structure,  only : altv, rltv, univol
  use m_FFT,                only : nfft, fft_box_size_WF &
       &                         , m_FFT_alloc_WF_work &
       &                         , m_FFT_dealloc_WF_work
  use m_Timing,             only : tstatc0_begin, tstatc0_end
  use m_Ionic_System,       only : ntyp,natm,natm2,ityp,iatomn,m_IS_pack_all_ions_in_uc &
       &                         , iwei,pos
  use m_PseudoPotential,    only : ival,modnrm,dk_wan,qitg_wan &
       &                         , m_PP_include_vanderbilt_pot &
       &                         , n_non0_lmtxlmt,index_lmt1_lmt2,lmta,nlmt,ilmt &
       &                         , ltp,taup,il2p,isph,iqitg,dl2p &
       &                         , m_PP_find_maximum_l
  use m_Files,              only : m_Files_open_nfwannier,nfwannier &
       &                         , m_Files_open_nfcntn_wannier &
       &                         , m_Files_open_nfpot_wannier &
       &                         , nfcntn_wannier,nfpot_wannier
  use m_Kpoints,            only : k_symmetry
  implicit none
  
  integer, public                     :: nwght
  real(kind=DP), private, allocatable :: wght(:) ! d(nwght)
  integer, public, allocatable        :: ghat(:,:) ! d(3,nwght) 

  integer, private                       :: nocc
  complex(kind=DP), private, allocatable :: zmat(:,:,:) ! d(nocc,nocc,nwght)
  complex(kind=DP), private, allocatable :: zn(:,:) ! d(nocc,nwght)

  real(kind=DP), private, allocatable :: umat(:,:) ! d(nocc,nocc)
  complex(kind=DP), private, allocatable :: bmat(:,:) ! d(nocc,nocc)
  real(kind=DP), private, allocatable :: mmat(:,:) ! d(nocc,nocc)

  integer, private                    :: ntri ! (nocc-1)*nocc/2
  real(kind=DP), private, allocatable :: amat(:) ! d(ntri)
  real(kind=DP), private, allocatable :: grad(:) ! d(ntri)
  real(kind=DP), private :: omega
  real(kind=DP), private :: gradmax

  complex(kind=DP), private, allocatable :: rmat(:,:) ! d(nocc,nocc)
  real(kind=DP), private, allocatable    :: lam(:) ! d(nocc)

  integer, private                    :: iteration

  complex(kind=DP), private, allocatable :: ftqe(:,:,:,:)!d(nlmt,nlmt,natm,nwght)

contains
  subroutine m_Wan_gen_weight(nfout)
    integer, intent(in) :: nfout
    integer :: i
    real(kind=DP) :: g(3,3),w
    real(kind=DP), parameter :: eps = 1.d-10

    nwght = 6
    allocate(wght(nwght),ghat(3,nwght))

    g = matmul(transpose(altv),altv)

    nwght = 0
    ! g^=(1,0,0)
    w = g(1,1)-g(1,2)-g(1,3)
    nwght = nwght + 1
    wght(nwght) = w
    ghat(1:3,nwght) = (/1,0,0/)
    ! g^=(0,1,0)
    w = g(2,2)-g(1,2)-g(2,3)
    nwght = nwght + 1
    wght(nwght) = w
    ghat(1:3,nwght) = (/0,1,0/)
    ! g^=(0,0,1)
    w = g(3,3)-g(1,3)-g(2,3)
    nwght = nwght + 1
    wght(nwght) = w
    ghat(1:3,nwght) = (/0,0,1/)
    ! g^=(1,1,0)
    w = g(1,2)
    if(abs(w) > eps) then
       nwght = nwght + 1
       wght(nwght) = w
       ghat(1:3,nwght) = (/1,1,0/)
    end if
    ! g^=(1,0,1)
    w = g(1,3)
    if(abs(w) > eps) then
       nwght = nwght + 1
       wght(nwght) = w
       ghat(1:3,nwght) = (/1,0,1/)
    end if
    ! g^=(1,1,0)
    w = g(2,3)
    if(abs(w) > eps) then
       nwght = nwght + 1
       wght(nwght) = w
       ghat(1:3,nwght) = (/0,1,1/)
    end if

    if(printable) then
       write(nfout,*) '!Wannier: nwght=',nwght
       do i=1,nwght
          write(nfout,'(" !Wannier: i=",i1," wght=",f18.5," g^=",3i2)') &
          & i,wght(i),ghat(1:3,i)
       end do
    end if

  end subroutine m_Wan_gen_weight


end module m_Wannier
