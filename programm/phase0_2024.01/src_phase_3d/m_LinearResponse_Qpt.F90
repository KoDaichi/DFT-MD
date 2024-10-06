!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  PROGRAM: TDLRMAIN
!
!  AUTHOR(S): K. Tagami et al   Aug. 1 2011
!
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)
!
!=======================================================================

! ======================================================================
!  This is a module for setting Q-points and their related variables 
!  for the LR-TDDFT.
! ======================================================================

! ======================= history ======================================
!  ver 1.0 :  2010/3/31
!               applicable to the crystal such as Si.  
!  ver 2.0 :  2011/3/31
!               applicable to the isolated molecule such as C6H6.
!
! ======================================================================

module m_LinearResponse_Qpt

  use m_Const_Parameters,               only : DP, CARTS, BUCS, PAI, PAI2, PAI4,  ELECTRON, INVERSE, ON, OFF, CMPLDP
  use m_Crystal_Structure,              only : altv,rltv  
  use m_Files,                          only : nfout

  use m_Kpoints,                        only : kv3_ek, kv3, vkxyz_ek, qwgt_ek, & 
       &                                       vkxyz, qwgt
  use m_Control_Parameters,             only : ipri, nspin, printable

  use m_PlaneWaveBasisSet,    only : ngabc,igf,kg1,kg, iba, nbase, nbase_gamma
  use m_IterationNumbers,     only : nk_in_the_process, nk_converged

  use m_Parallelization,         only : ista_kngp, iend_kngp
  use m_PlaneWaveBasisSet,        only : kg0

! --
  use m_Control_Parameters,     only : Num_q_Points, sw_LinearResponse
  use m_PseudoPotential,        only : norm_q_plus_G_LR, &
       &                               vec_q_plus_G_LR
!
  use m_LinearResponse_Control,         only : vqxyz, vec_q, &
       &                                       nrd_efermi
  use m_LinearResponse_Control,         only : sw_tddft, sw_LongWaveLimit
  use m_LinearResponse_Control, only : nmax_G_LR, sw_LongWaveLimit, ipriqpt
! ------------------ for subroutine m_LR_kp_cr_kpoints_table ----------------------
  use m_Kpoints,                only : way_ksample, k_sample_mesh, ip20, ip2cub, iwt, &
       &                               nxyz_tetra, ip01, ip21, ip10, ip02, ip12, &
       &                               iv21, iu21
  use m_Control_Parameters,     only : ipri_kp
  use m_Const_Parameters,       only : SKPS_DIRECT_IN, WHOLE_BZ, GENERAL, &
       &                               GENERAL_LARGER
  use m_Crystal_Structure,      only : nbztyp, nbztyp_spg, il, ngen, inv, igen, jgen, &
       &                               a, b, c, ca, cb, cc
! ------------------------------- for USPP ---------
  use m_Const_Parameters,       only : SKIP, zi
  use m_PseudoPotential,        only : nmax_q_plus_G_LR, nmax_q_LR, qitg_LR, &
       &                               nqitg, nlmta, nlmt, &
       &                               m_PP_include_vanderbilt_pot, &
       &                               m_PP_find_maximum_l, &
       &                               ilmt, &
       &                               ltp, taup, il2p, isph, iqitg, dl2p
  use m_Ionic_System,           only : ntyp

! -------------------------------------------------
  Implicit None
  include 'mpif.h'

! ----------------------- for USPP ---------------
  complex(kind=CMPLDP),allocatable :: ftqval(:,:,:,:,:)
! ------------------------------------------------

contains
! ----------------------------------------------------------
  subroutine m_LR_gen_q_points
    if ( sw_tddft==ON ) then
       Num_q_Points = 1;  
       allocate( vqxyz( Num_q_Points, 3, 2 ) )
       Call q_points_in_BUCS_CARTS
    endif
  end subroutine m_LR_gen_q_points

  subroutine m_LR_del_q_points
    if ( allocated( vqxyz ) ) Deallocate( vqxyz )
  end subroutine m_LR_del_q_points

  subroutine q_points_in_BUCS_CARTS
    integer                    :: i, j
    real(kind=DP)              :: fact 
! 
    vqxyz(1,:,CARTS) = vec_q(:)
    Do i = 1, 3
       Do j = 1, Num_q_Points
          vqxyz( j,i,BUCS ) = ( altv(1,i)*vqxyz( j,1,CARTS )&
               &           +  altv(2,i) *vqxyz( j,2,CARTS)&
               &           +  altv(3,i) *vqxyz( j,3,CARTS) )/PAI2
       End Do
    End Do
! ---------------------------- check of BUCS system ---------
    Do i = 1, 3
       Do j = 1, Num_q_Points
          fact =   rltv(i,1) *vqxyz( j,1,BUCS )&
               & + rltv(i,2) *vqxyz( j,2,BUCS )&
               & + rltv(i,3) *vqxyz( j,3,BUCS )
          if ( dabs(fact - vqxyz( j,i,CARTS ) ) > 1.d-5 ) then
             if ( ipriqpt >=1 ) then
                write(nfout,*) ' !D (i = ',i,',j = ',j,')'
                write(nfout,*) vqxyz( j,i,BUCS),fact, vqxyz( j,i,CARTS)
             end if
             stop ' Coordinate transformation is invaild in LR'
          end if
       End do
    End do
  end subroutine q_points_in_BUCS_CARTS

  subroutine m_LR_set_Vecs_k_minus_q             ! ( k - q ) 
    integer i, ik, kv3_org
    integer j1, k1, ispin
    real(kind=DP), allocatable :: vkxyz_tmp( :,:,: ), qwgt_tmp(:)
! -------------------------------- start --------------
    if ( sw_LongWaveLimit==ON ) then
       write(*,*) 'stop here '
       stop
    endif
    kv3_org = kv3 / ( Num_q_Points + 1 )
    Allocate( vkxyz_tmp( kv3_org, 1:3, 1:2 ) )
    Allocate( qwgt_tmp( kv3_org ) )
! -------------------------------- main --------
    Do ik=1, kv3_org
       Do i=1, 3
          vkxyz_tmp( ik, i, CARTS )  = vkxyz( ik, i, CARTS )
          vkxyz_tmp( ik, i, BUCS )   = vkxyz( ik, i, BUCS )
       End do
       qwgt_tmp( ik ) = qwgt( ik )
    End do
! ---
    j1 = 0;     k1 = 0
!
    Do ik=1, kv3_org / nspin
       Do i=0, Num_q_Points
          Do ispin =1, nspin
             if ( i==0 ) k1 = k1 + 1
             j1 = j1 + 1

             if ( i >= 1 ) then
                vkxyz( j1, 1:3, CARTS ) = vkxyz_tmp( k1, 1:3, CARTS ) &
                     &                  + vqxyz( i, 1:3, CARTS )
                vkxyz( j1, 1:3, BUCS  ) = vkxyz_tmp( k1, 1:3, BUCS  ) &
                     &                  + vqxyz( i, 1:3, BUCS )
                qwgt( j1 ) = qwgt_tmp( k1 )
            else
                vkxyz( j1, 1:3, CARTS ) = vkxyz_tmp( k1, 1:3, CARTS )
                vkxyz( j1, 1:3, BUCS  ) = vkxyz_tmp( k1, 1:3, BUCS  )
                qwgt( j1 ) = qwgt_tmp( k1 )
            endif
          End do
       End do
    End Do
! ---------------------------------- end -------------
    Deallocate( vkxyz_tmp, qwgt_tmp )

  end subroutine m_LR_set_Vecs_k_minus_q

! --------------------------
  subroutine m_LR_alloc_Vecs_q_plus_G
    allocate( norm_q_plus_G_LR( kg0, Num_q_points ) ) ;
    allocate( vec_q_plus_G_LR( kg0, 3, Num_q_points ) ) ;
    norm_q_plus_G_LR = 0.0d0
    vec_q_plus_G_LR = 0.0d0
  end subroutine m_LR_alloc_Vecs_q_plus_G

  subroutine m_LR_dealloc_Vecs_q_plus_G
    deallocate( norm_q_plus_G_LR, vec_q_plus_G_LR )
  end subroutine m_LR_dealloc_Vecs_q_plus_G

  subroutine m_LR_set_Vecs_q_plus_G
    integer :: i, nq
!!$    real(kind=DP) :: ga, gb, gc, g2, ttr(3)
    real(kind=DP) :: ga, gb, gc, g2, ttr(6)
    real(kind=DP) :: c1, c2, c3
! ------------------------------- start --------------
    Call getttr (rltv,ttr)
! ------------------------------- main --------------
    Do nq=1, Num_q_Points
!!!!       Do i=1, nmax_G
       Do i=1, kg0
          if ( sw_LongWaveLimit == OFF ) then
             ga = real(ngabc(i,1),kind=DP) + vqxyz(nq,1,BUCS)
             gb = real(ngabc(i,2),kind=DP) + vqxyz(nq,2,BUCS)
             gc = real(ngabc(i,3),kind=DP) + vqxyz(nq,3,BUCS)
          else
             ga = real(ngabc(i,1),kind=DP) 
             gb = real(ngabc(i,2),kind=DP) 
             gc = real(ngabc(i,3),kind=DP) 
          endif
          g2 = ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
               &   +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga
          g2 = sqrt( g2 )

          norm_q_plus_G_LR(i,nq) = g2

          c1 = rltv(1,1)*ga +rltv(1,2)*gb +rltv(1,3)*gc
          c2 = rltv(2,1)*ga +rltv(2,2)*gb +rltv(2,3)*gc
          c3 = rltv(3,1)*ga +rltv(3,2)*gb +rltv(3,3)*gc

          vec_q_plus_G_LR(i,1,nq) = c1
          vec_q_plus_G_LR(i,2,nq) = c2
          vec_q_plus_G_LR(i,3,nq) = c3
       End do
    End do

  end subroutine m_LR_set_Vecs_q_plus_G

! -------------------------------------------------------------
!! 
!!!                          For USPP 
!!
! -------------------------------------------------------------

  subroutine m_LR_alloc_Qitg
!!!!!!!!    nmax_G_LinearResponse = nmax_G;
    nmax_q_plus_G_LR = kg0;     
    nmax_q_LR = Num_q_Points
    allocate( qitg_LR( kg0, nqitg, Num_q_Points ) );  qitg_LR = 0.0d0;
  end subroutine m_LR_alloc_Qitg

  subroutine m_LR_dealloc_Qitg
    Deallocate( qitg_LR )
  end subroutine m_LR_dealloc_Qitg

  subroutine m_LR_alloc_Array_FTQVals
    allocate( ftqval( nlmt,nlmt,ntyp,Num_q_Points, kg0 ) )
    ftqval=( 0.d0, 0.d0 )
  end subroutine m_LR_alloc_Array_FTQVals

  subroutine m_LR_Dealloc_Array_FTQVals
    deallocate( ftqval )
  end subroutine m_LR_Dealloc_Array_FTQVals

  subroutine m_LR_calc_FTQ_times_Ylm
    real(kind=DP), allocatable :: ftqr(:,:), ftqi(:,:)
    integer, allocatable :: il3(:)

    integer              :: n, it, lmt1, lmt2
    integer              :: il1, il2, tau1, tau2, ilm3
    integer              :: l3, iiqitg
    integer              :: nq, ig
    real(kind=DP)        :: dq(3), ylm, ftqb
    integer              :: mdvdb

! --------------------------------- start ----------
    call m_PP_find_maximum_l(n)   !  n-1: maximum l
    n = (n-1) + (n-1) + 1
    allocate( il3(n**2) ); call substitute_il3(n**2,il3) ! -(b_Elec..)
    allocate( ftqr( Num_q_Points, kg0 ) ); ftqr = 0.0d0
    allocate( ftqi( Num_q_Points, kg0 ) ); ftqi = 0.0d0
! ------------------------------------ main ------------
    Do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if ( mdvdb == SKIP ) cycle

       Do lmt1 = 1, ilmt(it)
          il1 = ltp(lmt1,it);   tau1 = taup(lmt1,it)

          Do lmt2 = lmt1, ilmt(it)
             il2 = ltp(lmt2,it);    tau2 = taup(lmt2,it)

             ftqr = 0.d0;    ftqi = 0.d0

             Do n=1, il2p(lmt1,lmt2,it)
                ilm3 = isph(lmt1,lmt2,n,it); l3=il3(ilm3)
                iiqitg = iqitg(il1,tau1,il2,tau2,l3+1,it)
                if ( iiqitg==0 ) cycle

                Do nq=1, Num_q_Points
!!!                   Do ig=1, nmax_G
                   Do ig=1, kg0
                      dq(:) = vec_q_plus_G_LR( ig, :, nq )
                      call sphrp2_for_Berry( ilm3, dq, ylm )

                      ftqb = qitg_LR( ig, iiqitg, nq ) &
                           & *dl2p( lmt1,lmt2,n,it ) *ylm

                      if ( mod(l3,2)==0 ) then
                         ftqr(nq,ig) = ftqr(nq,ig) +real(zi**(-l3))*ftqb
                      else
                         ftqi(nq,ig) = ftqi(nq,ig) +dimag(zi**(-l3))*ftqb
                      end if
                   End do
                End Do
             End do

             Do nq=1, Num_q_Points
!!!                Do ig=1, nmax_G
                Do ig=1, kg0
                   ftqval(lmt1,lmt2,it,nq,ig) = dcmplx( ftqr(nq,ig),ftqi(nq,ig))
                   ftqval(lmt2,lmt1,it,nq,ig) = ftqval(lmt1,lmt2,it,nq,ig)
                End do
             End do
!
          End do
       End do
    End do
!    write(*,*) 'ftqavl = ', ftqval
!    stop
! ------------------------------- end ---------------------
    deallocate( ftqr, ftqi )
! ----------
!    Do it = 1, ntyp
!       mdvdb = m_PP_include_vanderbilt_pot(it)
!       if ( mdvdb == SKIP ) cycle

!       Do lmt1 = 1, ilmt(it)
!          Do lmt2 = 1, ilmt(it)
!             Do ig=1, nmax_G
!                write(*,*) 'copa ',lmt1,lmt2, ftqval(lmt1,lmt2,it,1,ig), ftqval(lmt2,lmt1,it,1,ig)
!             End do
!          End do
!       End do
!    End do
!    stop
  end subroutine m_LR_calc_FTQ_times_Ylm

! -------------------------------------------------
!!
!!!                        Kpoint 
!! 
! -------------------------------------------------
  subroutine m_LR_set_kp_cr_kpt_table
!    for tetrahedron method
#ifndef NO_TETRAHEDRON
!!$    integer :: nx1,ny1,nz1,nd
    integer :: nx,ny,nz
    integer :: nxx,nyy,nzz
    integer :: ill,ii,lmnp0, lmnp1, lmnp2
    real(kind=DP), pointer, dimension(:,:) :: pa0,pb0,pb
    integer,       pointer, dimension(:,:) :: ka0,ka2
    integer,       pointer, dimension(:,:) :: ip2cub_wk
    integer,       pointer, dimension(:)   :: nstar2
    
    integer :: kv3_div_tmp
    integer :: np0, np1, np2
    integer :: nx1, ny1, nz1, nd
    integer :: ipri_qp_count = 0, ipri_qp_t = 0

! ----------------------------------- start --------------
    if ( ipriqpt >= 2) write(nfout,*) '-- create k-table for tetrahedron method in LR --'
    if ( way_ksample == SKPS_DIRECT_IN ) then
        stop 'way_ksample = SKPS_DIRECT_IN case is not&
             & supportted in tetrahedron method.'
    end if
!
    kv3_div_tmp = kv3 / nspin / ( Num_q_Points + 1 )
! --------------------------------------- main ----------
    if ( nbztyp.eq.WHOLE_BZ ) then
       np0 = (k_sample_mesh(1,2)+1)*(k_sample_mesh(2,2)+1)*(k_sample_mesh(3,2)+1)
       np2 = k_sample_mesh(1,2)*k_sample_mesh(2,2)*k_sample_mesh(3,2)
       np1 = np2
       
       if ( ipriqpt>=1 ) write(nfout,'(" np0,np1,np2 in LR = ",3i6)') np0,np1,np2
       if ( np2 /= kv3_div_tmp ) then
          if ( ipriqpt >= 1 ) then
             write(nfout,*) ' np0,np1,np2 in LR',np0,np1,np2
             write(nfout,*) ' kv3_div_tmp ', kv3_div_tmp
             write(nfout,*) ' np2  /=  (kv3/nspin/(Num_q_points+1) )'
          end if
          stop
       endif
       
       allocate(ip20(np0)) ; ip20 = 0
       if ( ipriqpt>=2 ) write(nfout,'(" !kp ip20 is allocated in LR")')
       allocate(iwt(np2)) ; iwt = 0
       allocate(ip2cub(np1)) ; ip2cub = 0
       allocate(nxyz_tetra(3)) ; nxyz_tetra = 0
        
       nx=k_sample_mesh(1,2); ny=k_sample_mesh(2,2);  nz=k_sample_mesh(3,2)
       ill = 1;         lmnp0=np0;           lmnp1=np1

       allocate(ip01(np1)); ip01 = 0; 
       allocate(pa0(3,np0)) ; pa0 = 0.d0

       call nskma0(ill,nx,ny,nz,nxyz_tetra(1),nxyz_tetra(2) &
            &,nxyz_tetra(3),nx1,ny1,nz1,nd)
       call nskp00(nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3) &
            &,nx1,ny1,nz1,nd,lmnp0,np0,pa0)

!! When nbztyp=1, ip10 in FLAPW program should be ip20 in this program !!
       call nskpbm(np0,lmnp0,lmnp1,pa0,np1,ip20,ip01)
       deallocate(pa0);     deallocate(ip01)
    else

       if ( nbztyp_spg == GENERAL .or.nbztyp_spg == GENERAL_LARGER ) then
          nx=k_sample_mesh(1,1); ny=k_sample_mesh(2,1); nz=k_sample_mesh(3,1)
          
          call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
          lmnp0=(nxx+1)*(nyy+1)*(nzz+1);    lmnp1=lmnp0;     lmnp2=lmnp0
          
          allocate(nxyz_tetra(3)) ; nxyz_tetra = 0
          allocate(ip10(lmnp0))  ; ip10 = 0
          allocate(ip20(lmnp0))  ; ip20 = 0
          allocate(ip01(lmnp1))    ; ip01 = 0
          allocate(ip02(lmnp2))    ; ip02 = 0
          allocate(ip21(lmnp1))    ; ip21 = 0
          allocate(ip12(lmnp2))    ; ip12 = 0
          allocate(iu21(lmnp1))    ; iu21 = 0
          allocate(iv21(lmnp1))    ; iv21 = 0
          allocate(nstar2(lmnp2))  ; nstar2 = 0
          allocate(pa0(3,lmnp0))  ; pa0 = 0
          allocate(pb0(3,lmnp0))  ; pb0 = 0
          allocate(pb(3,lmnp2))  ; pb = 0
          allocate(ka0(4,lmnp0))  ; ka0 = 0
          allocate(ka2(4,lmnp2))  ; ka2 = 0

          ipri_qp_t = ipri_kp *(1-ipri_qp_count)
          if ( ipri_qp_count == 0 ) ipri_qp_count = 1
          call setkp0_n(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
               & ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
               & ,nx,ny,nz,nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3) &
               & ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21 &
               & ,nstar2,pa0,pb0,pb,ka0,ka2 &
               & ,ipri_qp_t)
       else
          nx=k_sample_mesh(1,1);  ny=k_sample_mesh(2,1);   nz=k_sample_mesh(3,1)
          
          if ( ipri_kp>=2 ) write(nfout,'(" !kp nx,ny,nz = ",3i6)') nx,ny,nz
!!$           if(nbztyp_spg.eq.SIMPLE_CUBIC) il = 1
!!$           if(nbztyp_spg.eq.BCC) il = 3
!!$           if(nbztyp_spg.eq.FCC) il = 2
!!$           if(nbztyp_spg.eq.DIAMOND) il = 2
!!$           if(nbztyp_spg.eq.HEXAGONAL) il = 0

          call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
          if ( ipri_kp>=2 ) write(nfout,'(" !kp nxx,nyy,nzz = ",3i6)') nxx,nyy,nzz
          
          lmnp0=(nxx+1)*(nyy+1)*(nzz+1);         lmnp1=lmnp0;      lmnp2=lmnp0
          
          allocate(nxyz_tetra(3)) ; nxyz_tetra = 0
          allocate(ip10(lmnp0))  ; ip10 = 0
          allocate(ip20(lmnp0))  ; ip20 = 0
          allocate(ip01(lmnp1))    ; ip01 = 0
          allocate(ip02(lmnp2))    ; ip02 = 0
          allocate(ip21(lmnp1))    ; ip21 = 0
          allocate(ip12(lmnp2))    ; ip12 = 0
          allocate(iu21(lmnp1))    ; iu21 = 0
          allocate(iv21(lmnp1))    ; iv21 = 0
          allocate(nstar2(lmnp2))  ; nstar2 = 0
          allocate(pa0(3,lmnp0))  ; pa0 = 0
          allocate(pb0(3,lmnp0))  ; pb0 = 0
          allocate(pb(3,lmnp2))  ; pb = 0
          allocate(ka0(4,lmnp0))  ; ka0 = 0
!!$          allocate(ka2(4,lmnp2))  ; ka2 = 0

!!$           call setkp0_default_n(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
!!$                &               ,nbztyp_spg,altv,nx,ny,nz &
!!$                &               ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
!!$                &               ,nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3)&
!!$                &               ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21 &
!!$                &               ,nstar2,pa0,pb0,pb,ka0,ka2 &
!!$                &               ,ipri_kp)
          ipri_qp_t = ipri_kp*(1-ipri_qp_count)
          if ( ipri_qp_count==0 ) ipri_qp_count = 1
          call setkp0_default_n(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
               &               ,nx,ny,nz &
               &               ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
               &               ,nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3)&
               &               ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21 &
               &               ,nstar2,pa0,pb0,pb,ka0 &
               &               ,ipri_qp_t)
          
       end if

       if ( ipri_kp>=2 ) write(nfout,'(" np0,np1,np2 = ",3i6)') np0,np1,np2
       if ( ipri_kp>=1 ) then
          write(nfout,'("!Kp nxyz_tetra(1:3) = ",3i8)') nxyz_tetra(1:3)
       end if
       
       if ( np2 /= kv3_div_tmp ) then
          if(ipri_kp>=1) then
             write(nfout,*) ' np0,np1,np2 ',np0,np1,np2
             write(nfout,*) ' kv3_div_tmp ',kv3_div_tmp
             write(nfout,*) ' np2  /=  (kv3/nspin/(Num_Qpoints+1) )'
          end if
          stop
       endif
       
       allocate(iwt(np2)) ; iwt = 0
       allocate(ip2cub(np1)) ; ip2cub = 0
       
       deallocate(ip10); deallocate(ip01); deallocate(ip02); deallocate(ip21)
       deallocate(ip12); deallocate(iu21); deallocate(iv21)
       deallocate(nstar2); deallocate(pa0); deallocate(pb0); deallocate(pb)
       deallocate(ka0)
       if ( nbztyp_spg == GENERAL .or.nbztyp_spg == GENERAL_LARGER ) &
            & deallocate(ka2)
    end if
    
    if ( ipri_kp>=2 ) then
       write(nfout,*) 'ip20'
       write(nfout,'(" np0 = ",i8)') np0
       write(nfout,*) (ip20(ii),ii=1,np0)
       write(nfout,'(" nx1,ny1,nz1,nd = ",4i6)') nx1,ny1,nz1,nd
    end if
    
    allocate(ip2cub_wk(9,nxyz_tetra(1)*nxyz_tetra(2)*nxyz_tetra(3)))
    ip2cub_wk = 0.d0
    call wtetra &
         &  (nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3),np0,np2,ip20 &
         &  ,iwt,ip2cub &
         &  ,ip2cub_wk)
    deallocate(ip2cub_wk)
#endif
  end subroutine m_LR_set_kp_cr_kpt_table
! -------------------------------------------------------------

end module m_LinearResponse_Qpt

