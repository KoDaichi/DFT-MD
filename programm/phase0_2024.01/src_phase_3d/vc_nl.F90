   Subroutine vc_nl(nx,ny,nz,nfftx,nffty,nfftz,nspin,rhoxyz,Exc,univol,rinplw,vcnlxyz,ista,iend)

     use m_Crystal_Structure,  only : altv
     use m_Total_Energy,       only : m_TE_what_is_edeltb_now
     use m_Ionic_System,       only : natm2
     use m_IterationNumbers,   only : iteration
     use m_Const_Parameters,   only : PAI

!$   use omp_lib

      Implicit None



!********************************* Note *******************************************
! This program calculates non-local exchange-correlation potential using the 
! van der Waals density functional (vdW-DF).
! The Algorism follows Thonhauser's self-consistent method.
!            (Reference; T.Thonhauser et al. PRB 76 (2007) 125112.)
! E_total = Ex_GGA + Ec_LDA + Ec_non-local
! v_total = vx_GGA + vc_LDA + vc_non-local
!
! Periodic systems are assumed.
! Only orthogonal corrdinates are applicable.
! Atomic unit (Hartree) is used.
!
!
!
! +++++++++++ Algorithm ++++++++++
!    ggaxcp0
! (m_XC_Potential.F90 in PHASE)
!
!      |
!      |
!      V 
!
!   vc_nl    --->   d_rho
!                  gauleg
!                 calPhi1  --->  d_q0
!                 calPhi2  --->  d_q0
!                 calPhi3  --->  d_q0
!                 calPhi4  --->  d_q0
!                    Phi1  --->  d_q0
!                    Phi2  --->  d_q0
!                    Phi3  --->  d_q0
!                    Phi4  --->  d_q0
!
!
!
! ++++++++ Input & Output ++++++++
! Input
!   rhoxyz(nx*ny*nz,nspin)  : Total density
!   nx,ny,nz                : Number of the grid points
!   nspin (dummy)           : Spin is not considered
!
! Output
!   vcnlxyz(nx*ny*nz,nspin) : Total Correlation potential
!     = vc(:,:,:) 
!     = vcnl + vcLDA
!         vcnl              : Non-local Correlation potential
!         vcLDA             : LDA Correlation potential
!
!
!
! ++++++ Internal parameters +++++
! Each internal parameter is set to restrict the error less than 10 meV.
!
! na(=nb) ( 30    )        : Gauss-Legendre integration is used.
! a2(=b2) ( 60    )        : Upper limit of the integral of dadb.
! etai    ( 1.3   )        : Radius of analytical integrating sphere (a.u.)
! eta1    ( 8     )        : Radius of numerical integrating sphere  (a.u.)
! eta2    ( 40    )        : Cutoff of asymptotic function     (unit cells)
! adj     ( 1     )        : Parameter to adjust error in Estimation Time
!
!
!
! ++++ Calculation procedure +++++
! 1. Calculate the coefficients of the asymptotic function.
!
! 2. Call d_rho and calculate differential of total density.
!     --- 7 points formula is used for differentiation.
!
! 3. Call Phi1 to Phi4 (and calPhi1 to calPhi4) to calculate
!             the functional derivative of the core function.
!
! 4. Calculate the local correlation potentail using LDA.
!
! 5. Output total correlation potential.
!     --- vcnlxyz = vc = vcnl + vcLDA
!
!
!
!                                              Written by Youky Ono 
!**********************************************************************************



!++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
      Integer  nx,ny,nz,nfftx,nffty,nfftz,nspin,nrxyz,iterstart
      Integer  ista,iend
      Real(8) da,db,a1,a2,ax,ay,az,rxyz,edel,sedel
      Integer  ci,cj,ck,ca,cb,na,nb,cD
      Parameter (sedel=1.e-4,iterstart=12)
      Parameter (na=30,a1=0,a2=60)

! Gauss-Legendre integration
      Real(8) xi(na),wi(na)

      Integer ndel,cdel,nphiD,cphiD,maxnphiD
      Real(8) ddel,del,phiD,phix,phiy,mphiD,maxD,minD,maxDii,minDii
      Parameter (ndel=20,nphiD=200,maxD=10,minD=0.1)
      Real(8), Allocatable :: &
&             phi1dD(:,:),phi2dD(:,:),phi3dD(:,:),phi4dD(:,:)

      Integer  cir,cix,ciy,ciz&
&             ,cjr,cjx,cjy,cjz&
&             ,ckr,ckx,cky,ckz&
&             ,cr ,cx ,cy ,cz,lcx,hcx,lcy,hcy,lcz,hcz
      Integer  c1x,c1y,c1z,c4x,c4y,c4z

      Real(8)  di,dj,dk,dx,dy,dz,dv,n,ni,nj,nk
      Integer  i,j,k,zxp,zxm,zyp,zym,zzp,zzm,tx,ty,tz

!!$      Real(8)  eta,etai,eta1,eta2
      Real(8)  eta,etai,eta1
      Integer eta2
      Parameter (eta=0.0000000001,etai=1.3,eta1=8,eta2=20)

      Real(8), Allocatable ::  C6(:,:,:)
      Real(8), Allocatable ::  dxrho(:,:,:), d2xrho(:,:,:) &
&                            , dyrho(:,:,:), d2yrho(:,:,:) &
&                            , dzrho(:,:,:), d2zrho(:,:,:) &
&                            , dxyrho(:,:,:),dxzrho(:,:,:),dyzrho(:,:,:)

      Real(8) rs,x,nnx,nny,nnz,nn2,r,zo
      Real(8) nxp,nxm,nyp,nym,nzp,nzm

      Real(8) rhoxyz(nfftx*nffty*nfftz,nspin),vcnlxyz(nfftx*nffty*nfftz,nspin)

      Real(8), Allocatable ::  rho(:,:,:) &
&               ,vcnl(:,:,:),vc(:,:,:),vcLDA(:,:,:)
      Real(8) Exc,univol,rinplw,Ecnl,Ecnl1,Ecnl2,Ecnl3,Ecl,rhomin
      Parameter (rhomin=1.e-6)

      Real(8)  kF,ecLDA,ecLDAd,pi,Zab,wp,wq,q,m,h,e,GxcLDA
      Real(8)  gamma,a,b,C,d,q0,q0i,q0k
      Real(8)  v1,v2,v3,v4
      Parameter (pi=PAI)
      Parameter (e=1.0,m=1.0)
      Parameter (Zab=-0.8491)

      Real(8) T,W,phi,psi
      Real(8) naPhi1,naPhi2,naPhi3,naPhi4
      Real(8) n2nx,n2ny,n2nz,term1,term2,term3,term4,ns,excLDA,&
&      exLDAd,excLDAd,exc0,alfa1,alfa2,alfa3,alfa4,vcik,vcik1, &
&      vcik2,vcik3,q0dx,q0dy,q0dz,LD1,LD2,vii,vcii,v1d,v1dd,v2d,v2dd,&
&      v3d,v4d,Tv1,Tv2,Tv11,Tv22,Tv12,Tv13,Tv14,Tv23,Tv24,Td,  &
&      Tdd,Tddk,iv12,iv13,iv14,iv23,iv24,iv34

! Estimate running time
      Real(8) time
      Integer hh,mm,ss

!++++++++++++++++++++++++ end VARIABLES ++++++++++++++++++++++++



!---------------------- Calculation Start ----------------------

      Allocate(C6(0:nx-1,0:ny-1,0:nz-1))
      Allocate(dxrho(nx,ny,nz))
      Allocate(d2xrho(nx,ny,nz))
      Allocate(dyrho(nx,ny,nz))
      Allocate(d2yrho(nx,ny,nz))
      Allocate(dzrho(nx,ny,nz))
      Allocate(d2zrho(nx,ny,nz))
      Allocate(dxyrho(nx,ny,nz))
      Allocate(dxzrho(nx,ny,nz))
      Allocate(dyzrho(nx,ny,nz))
      Allocate(rho(nx,ny,nz))
      Allocate(vcnl(nx,ny,nz))
      Allocate(vc(nx,ny,nz))
      Allocate(vcLDA(nx,ny,nz))

      Do cix = 1,nx
      Do ciy = 1,ny
      Do ciz = 1,nz
         vcnl(cix,ciy,ciz)=0
         vcLDA(cix,ciy,ciz)=0
         vc(cix,ciy,ciz)=0
         C6(0:cix-1,0:ciy-1,0:ciz-1)=0
         dxrho(cix,ciy,ciz)=0
         dyrho(cix,ciy,ciz)=0
         dzrho(cix,ciy,ciz)=0
         d2xrho(cix,ciy,ciz)=0
         d2yrho(cix,ciy,ciz)=0
         d2zrho(cix,ciy,ciz)=0
         dxyrho(cix,ciy,ciz)=0
         dxzrho(cix,ciy,ciz)=0
         dyzrho(cix,ciy,ciz)=0

      End Do
      End Do
      End Do
      rho=0
      vcnlxyz=0
      da=0
      nb = na
      db = da

      dx = altv(1,1)/Real(nx)
      dy = altv(2,2)/Real(ny)
      dz = altv(3,3)/Real(nz)
      dv = dx*dy*dz

      Do cix = 1,nx
      Do ciy = 1,ny
      Do ciz = 1,nz
         nrxyz=(ciz-1)*nffty*(nfftx)+(ciy-1)*(nfftx)+cix
         !!if(nrxyz<ista .or. nrxyz>iend) cycle
         rho(cix,ciy,ciz) = rhoxyz(nrxyz,1)
      End Do
      End Do
      End Do


!++++  Calculate the coefficients of the asymptotic function. ++++
      ax = dx*nx
      ay = dy*ny
      az = dz*nz

      Do cix = 1,nx
      Do ciy = 1,ny
      Do ciz = 1,nz
         C6(cix-1,ciy-1,ciz-1) = 0
      End Do
      End Do
      End Do

      Do ckx = 0,nx-1
      Do cky = 0,ny-1
      Do ckz = 0,nz-1
         Do tx=-1*eta2,eta2
         Do ty=-1*eta2,eta2
         Do tz=-1*eta2,eta2

            rxyz= &
&  ((ckx*dx+tx*ax)**2 + (cky*dy+ty*ay)**2 + (ckz*dz+tz*az)**2)**0.5

         If(rxyz.GT.eta1) Then
            C6(ckx,cky,ckz) = C6(ckx,cky,ckz) + (rxyz**(-6))
         End If

         End Do
         End Do
         End Do
      End Do
      End Do
      End Do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!########  Call d_rho and obtain differential of the rho #########
   Call d_rho(nx,ny,nz,dx,dy,dz,rho,dxrho,dyrho,dzrho&
&              ,d2xrho,d2yrho,d2zrho,dxyrho,dxzrho,dyzrho)
!#################################################################


!############## Call gauleg for Gauss-Legendre integral ##########
      Call gauleg(a1,a2,na,xi,wi)
!#################################################################


      Allocate(phi1dD(-ndel:ndel,-1:nphiD+1))
      Allocate(phi2dD(-ndel:ndel,-1:nphiD+1))
      Allocate(phi3dD(-ndel:ndel,-1:nphiD+1))
      Allocate(phi4dD(-ndel:ndel,-1:nphiD+1))

      Do cdel = -ndel,ndel
      Do cphiD = -1,nphiD+1
         phi1dD(cdel,cphiD) = 0
         phi2dD(cdel,cphiD) = 0
         phi3dD(cdel,cphiD) = 0
         phi4dD(cdel,cphiD) = 0
      End Do
      End Do

!######################  Call Phi1 to Phi4  ######################
  If(iteration.GE.1.AND.iteration.LT.iterstart) Then
  Else
   Call calPhi1(maxD,minD,etai,eta1,na,xi,wi,ndel,nphiD,phi1dD)
   Call calPhi2(maxD,minD,etai,eta1,na,xi,wi,ndel,nphiD,phi2dD)
   Call calPhi3(maxD,minD,etai,eta1,na,xi,wi,ndel,nphiD,phi3dD)
  End If
   Call calPhi4(maxD,minD,etai,eta1,na,xi,wi,ndel,nphiD,phi4dD)
!#################################################################

      Call CPU_TIME(time)
      ss=time
      hh=time/3600
      ss=Mod(ss,3600)
      mm=ss/60
      ss=Mod(ss,60)
      time=time-3600*hh-60*mm-ss

      Ecnl  = 0
!!!!!!!!!!!!!!!!! Hot Spot !!!!!!!!!!!!!!
!$OMP  PARALLEL DO &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(cix,ciy,ciz,ni,nb,db,gamma,dv,kF,rs,x,  &
!$OMP&  nnx,nny,nnz,nn2,n2nx,n2ny,n2nz,term1,term2,term3,term4,  &
!$OMP&  ns,exLDAd,ecLDAd,excLDAd,lcx,hcx,&
!$OMP&  lcy,hcy,lcz,hcz,ckx,cky,ckz,r,nk,di,dk,vcik,vcik1,vcik2, &
!$OMP&  cD,C,psi,maxDii,minDii,q0,q0i,q0k,q0dx,q0dy,q0dz,LD1,LD2,vii,vcii,alfa1,alfa2,alfa3,alfa4,  &
!$OMP&  a,b,h,v1,v1d,v1dd,v2,v2d,v2dd,v3,v3d,v4,naPhi1,naPhi2,naPhi3,naPhi4,   &
!$OMP&  v4d,W,Tv1,Tv2,Tv11,Tv22,Tv12,Tv13,Tv14,Tv23,Tv24,Td,Tdd, &
!$OMP&  Tddk,iv12,iv13,iv14,iv23,iv24,iv34,GxcLDA,excLDA,exc0,   &
!$OMP&  ddel,del,cdel,phiD,cphiD,phix,phiy,phi) &
!$OMP&  REDUCTION(+:Ecnl)
      Do cir = 1,nx*ny*nz
         cix = 1+(cir-MOD(cir-1+ny*nz,ny*nz)+1)/(ny*nz)
         ciy = 1+(cir-ny*nz*(cix-1)-(MOD(cir-1+nz,nz)+1))/nz
         ciz = (cir-ny*nz*(cix-1)-nz*(ciy-1))
         nrxyz=(ciz-1)*nffty*(nfftx)+(ciy-1)*(nfftx)+cix
         if(nrxyz<ista .or. nrxyz>iend) cycle

         naPhi1 = 0
         naPhi2 = 0
         naPhi3 = 0
         naPhi4 = 0

         ni = rho(cix,ciy,ciz)
         If(ni.GT.rhomin) Then

!######################  Call Phi1 to Phi4  ######################

   If(iteration.GE.0.AND.iteration.LT.iterstart) Then
   Else

   Call Phi1(maxD,minD,etai,eta1,na,xi,wi,nx,ny,nz,dx,dy,dz,rho,cix,ciy      &
&            ,ciz,dxrho,dyrho,dzrho,d2xrho,d2yrho,d2zrho,C6,ndel,nphiD,phi1dD,naPhi1)
   Call Phi2(maxD,minD,etai,eta1,na,xi,wi,nx,ny,nz,dx,dy,dz,rho,cix,ciy,ciz  &
&             ,dxrho,dyrho,dzrho,d2xrho,d2yrho,d2zrho                        &
&             ,dxyrho,dxzrho,dyzrho,C6,ndel,nphiD,phi2dD,naPhi2)
   Call Phi3(maxD,minD,etai,eta1,na,xi,wi,nx,ny,nz,dx,dy,dz,rho,cix,ciy,ciz  &
&             ,dxrho,dyrho,dzrho,d2xrho,d2yrho,d2zrho,C6,ndel,nphiD,phi3dD,naPhi3)
   End If

   If(iteration.GE.0.AND.iteration.LT.8) Then
   Else
   Call Phi4(maxD,minD,etai,eta1,na,xi,wi,nx,ny,nz,dx,dy,dz,rho,cix,ciy,ciz  &
&             ,dxrho,dyrho,dzrho,d2xrho,d2yrho,d2zrho,C6,ndel,nphiD,phi4dD,naPhi4)
   End If
!#################################################################

         vcnl(cix,ciy,ciz)= &
&           naPhi1 + naPhi2 + naPhi3 + naPhi4

         Ecnl = Ecnl + 0.5*ni*naPhi4

         End If
      End Do
!$OMP END PARALLEL DO
!!!!!!!!!!!!! End Hot Spot !!!!!!!!!!!!!!
      Ecnl = dv*Ecnl

      Call CPU_TIME(time)
      ss=time
      hh=time/3600
      ss=Mod(ss,3600)
      mm=ss/60
      ss=Mod(ss,60)
      time=time-3600*hh-60*mm-ss


!++++++++++++++++++++++  Calculate vcLDA  ++++++++++++++++++++++++
      Ecl = 0
      Do cjx = 1,nx
      Do cjy = 1,ny
      Do cjz = 1,nz
         nrxyz=(cjz-1)*nffty*(nfftx)+(cjy-1)*(nfftx)+cjx
         if(nrxyz<ista .or. nrxyz>iend) cycle

         nj = rho(cjx,cjy,cjz)

       If(nj.GT.rhomin) then

         rs = (3/(4*pi*nj))**0.3333333333333333
         x = rs/11.4

         term1 = (36*pi*(nj**4))**0.33333333333333333
         ecLDAd = -0.0333*(3*x*x*Log(1+1/x)-3*x+1.5-1/x)/(11.4*term1)

         GxcLDA = 0.5*((1+x**3)*Log(1+1/x)-x**2+x/2-1/3)
         ecLDA = -0.0666*GxcLDA

         vcLDA(cjx,cjy,cjz) = ecLDAd*nj + ecLDA
         Ecl = Ecl + dv*nj*ecLDA

       End If

      End Do
      End Do
      End Do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      Exc = Exc + (Ecnl + Ecl)/(univol*rinplw)

      Do cix = 1,nx
      Do ciy = 1,ny
      Do ciz = 1,nz
         nrxyz = (ciz-1)*nffty*nfftx+(ciy-1)*nfftx+cix
         if(nrxyz<ista .or. nrxyz>iend) cycle
         vc(cix,ciy,ciz) = vcnl(cix,ciy,ciz) + vcLDA(cix,ciy,ciz)
         !vcnlxyz(nrxyz,1) = vcnlxyz(nrxyz,1) + vc(cix,ciy,ciz)
         vcnlxyz(nrxyz,1) = vc(cix,ciy,ciz)
      End Do
      End Do
      End Do


      Deallocate(C6)
      Deallocate(dxrho)
      Deallocate(d2xrho)
      Deallocate(dyrho)
      Deallocate(d2yrho)
      Deallocate(dzrho)
      Deallocate(d2zrho)
      Deallocate(dxyrho)
      Deallocate(dxzrho)
      Deallocate(dyzrho)
      Deallocate(rho)
      Deallocate(vcnl)
      Deallocate(vc)
      Deallocate(vcLDA)
      Deallocate(phi1dD)
      Deallocate(phi2dD)
      Deallocate(phi3dD)
      Deallocate(phi4dD)



   End Subroutine vc_nl



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 
 
!** SUBROUTINE CPU_TIME ************************
      Subroutine CPU_TIME(time)
      Real(8) time
      Integer count,count_rate,hh,mm,ss
 
      Call System_clock(count,count_rate)
      If(count_rate.Eq.0) Then
         time=Real(count)/count_rate
      Else
         time=-1
      End If
 
      End Subroutine CPU_TIME
!** End SUBROUTINE CPU_TIME ********************





   Subroutine d_rho(nx,ny,nz,dx,dy,dz,rho,dxrho,dyrho,dzrho,d2xrho,d2yrho,d2zrho&
&                      ,dxyrho,dxzrho,dyzrho)
!$   use omp_lib
      Implicit None


!************************ Note *********************************
! This Algorism calculates the differential of the rho.
!
! Input
!   nx,ny,nz,dx,dy,dz,rho       : Information of the unit cell 
!                                 and the electron density.
!
! Output
!   drho,d2rho                  : 1st and 2nd order differential
!                                 of the rho.     
!
!                            Written by Youky Ono
!***************************************************************



!++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
      Integer  nx,ny,nz,cjx,cjy,cjz,i,j,k
      Real(8)  dx,dy,dz,nspin,dv,nj

      Real(8)  ::  rn(3,-3:3),rx(3,-3:3),ry(3,-3:3),rz(3,-3:3)
      Integer  ::  zx(-3:3),zy(-3:3),zz(-3:3)

      Real(8)  nxp,nxm,nyp,nym,nzp,nzm &
&             ,nxpyp,nxpym,nxmyp,nxmym &
&             ,nxpzp,nxpzm,nxmzp,nxmzm &
&             ,nypzp,nypzm,nymzp,nymzm 

      Real(8) :: rho(nx,ny,nz),dxrho(nx,ny,nz),d2xrho(nx,ny,nz) &
&                             ,dyrho(nx,ny,nz),d2yrho(nx,ny,nz) &
&                             ,dzrho(nx,ny,nz),d2zrho(nx,ny,nz) &
&           ,dxyrho(nx,ny,nz),dxzrho(nx,ny,nz),dyzrho(nx,ny,nz)
!++++++++++++++++++++++++ end VARIABLES ++++++++++++++++++++++++



!---------------------- Calculation Start ----------------------
      Do cjx = 1,nx
      Do cjy = 1,ny
      Do cjz = 1,nz

         Do j = -3,3
            zx(j) = MOD(2*nx+(cjx+j)-1,nx)+1
            zy(j) = MOD(2*ny+(cjy+j)-1,ny)+1
            zz(j) = MOD(2*nz+(cjz+j)-1,nz)+1

            rn(1,j) = rho(zx(j),cjy,cjz)
            rn(2,j) = rho(cjx,zy(j),cjz)
            rn(3,j) = rho(cjx,cjy,zz(j))
         End Do

         dxrho(cjx,cjy,cjz) = &
&           (rn(1,3)-9*rn(1,2)+45*rn(1,1)-45*rn(1,-1)+9*rn(1,-2)-rn(1,-3))/(60*dx)
         dyrho(cjx,cjy,cjz) = &
&           (rn(2,3)-9*rn(2,2)+45*rn(2,1)-45*rn(2,-1)+9*rn(2,-2)-rn(2,-3))/(60*dy)
         dzrho(cjx,cjy,cjz) = &
&           (rn(3,3)-9*rn(3,2)+45*rn(3,1)-45*rn(3,-1)+9*rn(3,-2)-rn(3,-3))/(60*dz)

      End Do
      End Do
      End Do



      Do cjx = 1,nx
      Do cjy = 1,ny
      Do cjz = 1,nz

         Do j = -3,3
            zx(j) = MOD(2*nx+(cjx+j)-1,nx)+1
            zy(j) = MOD(2*ny+(cjy+j)-1,ny)+1
            zz(j) = MOD(2*nz+(cjz+j)-1,nz)+1

            rn(1,j)  =  rho(zx(j),cjy,cjz)
            rn(2,j)  =  rho(cjx,zy(j),cjz)
            rn(3,j)  =  rho(cjx,cjy,zz(j))

            rx(1,j) = dxrho(zx(j),cjy,cjz)
            rx(2,j) = dxrho(cjx,zy(j),cjz)
            rx(3,j) = dxrho(cjx,cjy,zz(j))

            ry(1,j) = dyrho(zx(j),cjy,cjz)
            ry(2,j) = dyrho(cjx,zy(j),cjz)
            ry(3,j) = dyrho(cjx,cjy,zz(j))

            rz(1,j) = dzrho(zx(j),cjy,cjz)
            rz(2,j) = dzrho(cjx,zy(j),cjz)
            rz(3,j) = dzrho(cjx,cjy,zz(j))
         End Do

         d2xrho(cjx,cjy,cjz) = &
&           (rx(1,3)-9*rx(1,2)+45*rx(1,1)-45*rx(1,-1)+9*rx(1,-2)-rx(1,-3))/(60*dx)
         d2yrho(cjx,cjy,cjz) = &
&           (ry(2,3)-9*ry(2,2)+45*ry(2,1)-45*ry(2,-1)+9*ry(2,-2)-ry(2,-3))/(60*dy)
         d2zrho(cjx,cjy,cjz) = &
&           (rz(3,3)-9*rz(3,2)+45*rz(3,1)-45*rz(3,-1)+9*rz(3,-2)-rz(3,-3))/(60*dz)

         dxyrho(cjx,cjy,cjz) = &
&           (ry(1,3)-9*ry(1,2)+45*ry(1,1)-45*ry(1,-1)+9*ry(1,-2)-ry(1,-3))/(60*dx)
         dyzrho(cjx,cjy,cjz) = &
&           (rz(2,3)-9*rz(2,2)+45*rz(2,1)-45*rz(2,-1)+9*rz(2,-2)-rz(2,-3))/(60*dy)
         dxzrho(cjx,cjy,cjz) = &
&           (rx(3,3)-9*rx(3,2)+45*rx(3,1)-45*rx(3,-1)+9*rx(3,-2)-rx(3,-3))/(60*dz)

      End Do
      End Do
      End Do



   End Subroutine d_rho



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



   Subroutine gauleg(x1,x2,n,xi,wi)
     use m_Const_Parameters, only : PAI
!$   use omp_lib
      Implicit none

!************************ Note *********************************
!
!                            Written by Youky Ono
!***************************************************************



!++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
      Integer  m,j,i,n
      Real(8)  x1,x2,z1,z,xm,xl,pp,p3,p2,p1,pi,eta
      Parameter (pi=PAI , eta=0.0000000001)

      Real(8)  ::  xi(n),wi(n)
!++++++++++++++++++++++++ end VARIABLES ++++++++++++++++++++++++



!---------------------- Calculation Start ----------------------
      m=(n+1)/2
      xm=0.5*(x2+x1)
      xl=0.5*(x2-x1)
      z1=0.d0
      Do i=1,m
         z=COS(pi*(i-0.25)/(n+0.5))

         Do While (ABS(z-z1).GT.eta)
            p1=1.0
            p2=0.0

            Do j=1,n
               p3=p2
               P2=p1
               p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j
            End Do

            pp=n*(z*p1-p2)/(z*z-1.0)
            z1=z
            z=z1-p1/pp
            
         End Do

         xi(i) = xm-xl*z
         xi(n+1-i) = xm+xl*z
         wi(i) = 2.0*xl/((1.0-z*z)*pp*pp)
         wi(n+1-i) = wi(i)

      End Do


! End of Subroutine gauleg
   End Subroutine gauleg



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



   Subroutine calPhi1(maxD,minD,etai,eta1,na,xi,wi,ndel,nphiD,phi1dD)
     use m_Const_Parameters, only : PAI
!$   use omp_lib
      Implicit None


!************************ Note *********************************
!
!                            Written by Youky Ono
!***************************************************************



!++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
      Integer  ci,cj,ck,ca,cb,na,nb

! Gauss-Legendre integration
      Real(8) xi(na),wi(na)

      Real(8) ddel,dphiD,del,phiD,tmp,phix,phiy,maxD,minD,etai,eta1
      Integer ndel,cdel,nphiD,cphiD,cdelD
      Real(8) phi1dD(-ndel:ndel,-1:nphiD+1)

      Integer  cir,cix,ciy,ciz&
&             ,cjr,cjx,cjy,cjz&
&             ,ckr,ckx,cky,ckz&
&             ,cr ,cx ,cy ,cz,lcx,hcx,lcy,hcy,lcz,hcz
      Real(8)  di,dj,dk,dx,dy,dz,dv,n,ni,nj,nk,r
      Integer  i,j,k

      Real(8)  eta
      Parameter (eta=0.0000000001)

      Real(8)  kF,excLDAd,exLDAd,ecLDAd,pi,Zab,wp,wq,q,m,h,e
      Real(8)  gamma,a,b,C,d
      Real(8)  v1,v2,v3,v4,v1d,v2d
      Parameter (pi=PAI)
      Parameter (e=1,m=1)
      Parameter (Zab=-0.8491)

      Real(8) T,W,phi,psi,Tv1,Tv2
!++++++++++++++++++++++++ end VARIABLES ++++++++++++++++++++++++



!---------------------- Calculation Start ----------------------
      gamma = 4*pi/9

! phi(del,phiD)
      ddel = 1/Real(ndel+1)
      dphiD = (maxD-minD)/(Real(nphiD))


!$OMP  PARALLEL DO &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(cdel,cphiD,del,phiD,di,dk,phi,ca,cb,a,b,  &
!$OMP&  h,v1,v1d,v2,v2d,v3,v4,W,Tv1,Tv2)

   Do cdel = -ndel,ndel
   Do cphiD = -1,nphiD+1

      del = Real(cdel)*ddel
      phiD = minD + Real(cphiD)*dphiD
      
      di = phiD*(1+del)
      dk = phiD*(1-del)

      phi = 0
      Do ca=1,na
      Do cb=1,na
         a = xi(ca)
         b = xi(cb)

         h = 1-Exp(-gamma*((a/di)**2)) 
         v1 = (a**2)/(2*h)
         v1d = -4*pi*(a**4)*(h-1)/(9*(di**3)*(h**2))

         h = 1-Exp(-gamma*((b/di)**2)) 
         v2 = (b**2)/(2*h)
         v2d = -4*pi*(b**4)*(h-1)/(9*(di**3)*(h**2))

         h = 1-Exp(-gamma*((a/dk)**2)) 
         v3 = (a**2)/(2*h)

         h = 1-Exp(-gamma*((b/dk)**2)) 
         v4 = (b**2)/(2*h)

         W = 2*((3-a*a)*b*Cos(b)*Sin(a) + (3-b*b)*a*Cos(a)*Sin(b)           &
&            + (a*a+b*b-3)*Sin(a)*Sin(b)                                    &
&            -3*a*b*Cos(a)*Cos(b))/((a*b)**3)

         Tv1 = ((1/((v1+v2)**2))*(1/((v1+v3)*(v2+v4))+1/((v1+v4)*(v2+v3)))  &
&              +(1/(v1+v2)+1/(v3+v4))                                       &
&              *(1/(((v1+v3)**2)*(v2+v4))+1/(((v1+v4)**2)*(v2+v3))))*(-0.5)

         Tv2 = ((1/((v1+v2)**2))*(1/((v1+v3)*(v2+v4))+1/((v1+v4)*(v2+v3)))  &
&              +(1/(v1+v2)+1/(v3+v4))                                       &
&              *(1/((v1+v3)*((v2+v4)**2))+1/((v1+v4)*((v2+v3)**2))))*(-0.5)

         phi = phi + &
&              wi(ca)*wi(cb)*((a*b)**2)*W*(Tv1*v1d+Tv2*v2d)

      End Do
      End Do

      phi = phi * 2*m*(e**4)/(pi**2)

      phi1dD(cdel,cphiD) = phi

   End Do
   End Do
!$OMP END PARALLEL DO

   End Subroutine calPhi1



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



   Subroutine calPhi2(maxD,minD,etai,eta1,na,xi,wi,ndel,nphiD,phi2dD)
     use m_Const_Parameters, only : PAI
!$   use omp_lib
      Implicit None


!************************ Note *********************************
!
!                            Written by Youky Ono
!***************************************************************



!++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
      Integer  ci,cj,ck,ca,cb,na,nb

! Gauss-Legendre integration
      Real(8) xi(na),wi(na)

      Real(8) ddel,dphiD,del,phiD,tmp,phix,phiy,maxD,minD,etai,eta1
      Integer ndel,cdel,nphiD,cphiD,cdelD
      Real(8) phi2dD(-ndel:ndel,-1:nphiD+1)

      Integer  cir,cix,ciy,ciz&
&             ,cjr,cjx,cjy,cjz&
&             ,ckr,ckx,cky,ckz&
&             ,cr ,cx ,cy ,cz,lcx,hcx,lcy,hcy,lcz,hcz
      Real(8)  di,dj,dk,dx,dy,dz,dv,n,ni,nj,nk,r
      Integer  i,j,k

      Real(8)  eta
      Parameter (eta=0.0000000001)

      Real(8)  kF,excLDAd,exLDAd,ecLDAd,pi,Zab,wp,wq,q,m,h,e
      Real(8)  gamma,a,b,C,d
      Parameter (pi=PAI)
      Parameter (e=1.0,m=1.0)
      Parameter (Zab=-0.8491)

      Real(8) T,W,phi,psi,vcik,vcik1,vcik2,vcii,vii
      Real(8) Tv1,Tv2,Tv11,Tv22,Tv12,Tv13,Tv14,Tv23,Tv24,Tdd
      Real(8) v1,v2,v3,v4,v1d,v2d,v3d,v4d,v1dd,v2dd

      Real(8) term1,term2,term3,term4,ns
      Real(8) iv12,iv13,iv14,iv23,iv24,iv34
!++++++++++++++++++++++++ end VARIABLES ++++++++++++++++++++++++



!---------------------- Calculation Start ----------------------
      gamma = 4*pi/9

! phi(del,phiD)
      ddel = 1/Real(ndel+1)
      dphiD = (maxD-minD)/(Real(nphiD))


!$OMP  PARALLEL DO &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(cdel,cphiD,del,phiD,di,dk,phi,ca,cb,a,b,  &
!$OMP&  h,v1,v1d,v1dd,v2,v2d,v2dd,v3,v4,W,               &
!$OMP&  iv12,iv13,iv14,iv23,iv24,iv34,Tv1,Tv2,Tv11,Tv22,Tv12,Tdd)

   Do cdel = -ndel,ndel
   Do cphiD = -1,nphiD+1

      del = Real(cdel)*ddel
      phiD = minD + Real(cphiD)*dphiD
      
      di = phiD*(1+del)
      dk = phiD*(1-del)

      phi = 0
      Do ca=1,na
      Do cb=1,na
         a = xi(ca)
         b = xi(cb)

         h = 1-Exp(-gamma*((a/di)**2)) 
         v1 = (a**2)/(2*h)
         v1d = -4*pi*(a**4)*(h-1)/(9*(di**3)*(h**2))
         v1dd = 2*((4*pi*(a**3)/(9*(di**3)))**2)*(1-h)*(2-3*h)/(h**3)

         h = 1-Exp(-gamma*((b/di)**2)) 
         v2 = (b**2)/(2*h)
         v2d = -4*pi*(b**4)*(h-1)/(9*(di**3)*(h**2))
         v2dd = 2*((4*pi*(b**3)/(9*(di**3)))**2)*(1-h)*(2-3*h)/(h**3)

         h = 1-Exp(-gamma*((a/dk)**2)) 
         v3 = (a**2)/(2*h)

         h = 1-Exp(-gamma*((b/dk)**2)) 
         v4 = (b**2)/(2*h)

         W = 2*((3-a*a)*b*Cos(b)*Sin(a) + (3-b*b)*a*Cos(a)*Sin(b)           &
&            + (a*a+b*b-3)*Sin(a)*Sin(b)                                    &
&            -3*a*b*Cos(a)*Cos(b))/((a*b)**3)


         iv12 = 1/(v1+v2)
         iv13 = 1/(v1+v3)
         iv14 = 1/(v1+v4)
         iv23 = 1/(v2+v3)
         iv24 = 1/(v2+v4)
         iv34 = 1/(v3+v4)


         Tv1 = -(1/2)*(iv12**2)*(iv13*iv24+iv14*iv23)              &
&              -(1/2)*(iv12+iv34)*(iv13**2*iv24+iv14**2*iv23)


         Tv2 = -(1/2)*(iv12**2)*(iv13*iv24+iv14*iv23)              &
&              -(1/2)*(iv12+iv34)*(iv13*iv24**2+iv14*iv23**2)


         Tv11 = iv12**3*(iv13*iv24+iv14*iv23)                      &
&              +iv12**2*(iv13**2*iv24+iv14**2*iv23)                &
&              +(iv12+iv34)*(iv13**3*iv24+iv14**3*iv23)


         Tv22 = iv12**3*(iv13*iv24+iv14*iv23)                      &
&              +iv12**2*(iv13*iv24**2+iv14*iv23**2)                &
&              +(iv12+iv34)*(iv13*iv24**3+iv14*iv23**3)


         Tv12 = iv12**3*(iv13*iv24+iv14*iv23)                      &
&              +(1/2)*iv12**2*(iv13*iv24**2+iv14*iv23**2)          &
&              +(1/2)*iv12**2*(iv13**2*iv24+iv14**2*iv23)          &
&              +(iv12+iv34)*(iv13**2*iv24**2+iv14**2*iv23**2)


      Tdd = Tv11*v1d**2 + Tv1*v1dd + Tv22*v2d**2 + Tv2*v2dd + 2*Tv12*v1d*v2d


         phi = phi + &
&              wi(ca)*wi(cb)*((a*b)**2)*W*Tdd

      End Do
      End Do

      phi = phi * 2*m*(e**4)/(pi**2)

      phi2dD(cdel,cphiD) = phi

   End Do
   End Do
!$OMP END PARALLEL DO

   End Subroutine calPhi2



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



   Subroutine calPhi3(maxD,minD,etai,eta1,na,xi,wi,ndel,nphiD,phi3dD)
     use m_Const_Parameters, only : PAI
!$   use omp_lib
      Implicit None

!************************ Note *********************************
!
!                            Written by Youky Ono
!***************************************************************



!++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
      Integer  ci,cj,ck,ca,cb,na,nb

! Gauss-Legendre integration
      Real(8) xi(na),wi(na)

      Real(8) ddel,dphiD,del,phiD,tmp,phix,phiy,maxD,minD,etai,eta1
      Integer ndel,cdel,nphiD,cphiD,cdelD
      Real(8) phi3dD(-ndel:ndel,-1:nphiD+1)

      Integer  cir,cix,ciy,ciz&
&             ,cjr,cjx,cjy,cjz&
&             ,ckr,ckx,cky,ckz&
&             ,cr ,cx ,cy ,cz,lcx,hcx,lcy,hcy,lcz,hcz
      Real(8)  di,dj,dk,dx,dy,dz,dv,n,ni,nj,nk,r
      Integer  i,j,k

      Real(8)  eta
      Parameter (eta=0.0000000001)

      Real(8)  kF,excLDAd,exLDAd,ecLDAd,pi,Zab,wp,wq,q,m,h,e
      Real(8)  gamma,a,b,C,d
      Parameter (pi=PAI)
      Parameter (e=1,m=1)
      Parameter (Zab=-0.8491)

      Real(8) T,W,phi,psi,vcik,vcik1,vcik2,vcii,vii
      Real(8) Tv1,Tv2,Tv11,Tv22,Tv12,Tv13,Tv14,Tv23,Tv24,Td,Tdd,Tddk
      Real(8) v1,v2,v3,v4,v1d,v2d,v3d,v4d,v1dd,v2dd

      Real(8) naPhi3,alfa3,term1,term2,term3,term4,ns
      Real(8) iv12,iv13,iv14,iv23,iv24,iv34
!++++++++++++++++++++++++ end VARIABLES ++++++++++++++++++++++++



!---------------------- Calculation Start ----------------------
      gamma = 4*pi/9

! phi(del,phiD)
      ddel = 1/Real(ndel+1)
      dphiD = (maxD-minD)/(Real(nphiD))


!$OMP  PARALLEL DO &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(cdel,cphiD,del,phiD,di,dk,phi,ca,cb,a,b,  &
!$OMP&  h,v1,v1d,v2,v2d,v3,v3d,v4,v4d,iv12,iv13,iv14,    &
!$OMP&  iv23,iv24,iv34,W,Tv1,Tv2,Tv11,Tv22,Tv12,Tv13,Tv14,Tv23,Tv24,Td,Tdd,Tddk)

   Do cdel = -ndel,ndel
   Do cphiD = -1,nphiD+1

      del = Real(cdel)*ddel
      phiD = minD + Real(cphiD)*dphiD

      di = phiD*(1+del)
      dk = phiD*(1-del)

      phi = 0
      Do ca=1,na
      Do cb=1,na
         a = xi(ca)
         b = xi(cb)

         h = 1-Exp(-gamma*((a/di)**2)) + eta
         v1 = (a**2)/(2*h)
         v1d = -4*pi*(a**4)*(h-1)/(9*(di**3)*(h**2))
         v1dd = 2*((4*pi*(a**3)/(9*(di**3)))**2)*(1-h)*(2-3*h)/(h**3)

         h = 1-Exp(-gamma*((b/di)**2)) + eta
         v2 = (b**2)/(2*h)
         v2d = -4*pi*(b**4)*(h-1)/(9*(di**3)*(h**2))
         v2dd = 2*((4*pi*(b**3)/(9*(di**3)))**2)*(1-h)*(2-3*h)/(h**3)

         h = 1-Exp(-gamma*((a/dk)**2)) + eta
         v3 = (a**2)/(2*h)
         v3d = -4*pi*(a**4)*(h-1)/(9*(di**3)*(h**2))

         h = 1-Exp(-gamma*((b/dk)**2)) + eta
         v4 = (b**2)/(2*h)
         v4d = -4*pi*(b**4)*(h-1)/(9*(di**3)*(h**2))

         W = 2*((3-a*a)*b*Cos(b)*Sin(a) + (3-b*b)*a*Cos(a)*Sin(b)           &
&            + (a*a+b*b-3)*Sin(a)*Sin(b)                                    &
&            -3*a*b*Cos(a)*Cos(b))/((a*b)**3)

         iv12 = 1/(v1+v2)
         iv13 = 1/(v1+v3)
         iv14 = 1/(v1+v4)
         iv23 = 1/(v2+v3)
         iv24 = 1/(v2+v4)
         iv34 = 1/(v3+v4)


         Tv1 = -(1/2)*(iv12**2)*(iv13*iv24 + iv14*iv23)              &
&              -(1/2)*(iv12 + iv34)*(iv13**2*iv24 + iv14**2*iv23)


         Tv2 = -(1/2)*(iv12**2)*(iv13*iv24 + iv14*iv23)              &
&              -(1/2)*(iv12 + iv34)*(iv13*iv24**2 + iv14*iv23**2)


         Tv11 = iv12**3*(iv13*iv24 + iv14*iv23)                      &
&              +iv12**2*(iv13**2*iv24 + iv14**2*iv23)                &
&              +(iv12 + iv34)*(iv13**3*iv24 + iv14**3*iv23)


         Tv22 = iv12**3*(iv13*iv24 + iv14*iv23)                      &
&              +iv12**2*(iv13*iv24**2 + iv14*iv23**2)                &
&              +(iv12 + iv34)*(iv13*iv24**3 + iv14*iv23**3)


         Tv12 = iv12**3*(iv13*iv24 + iv14*iv23)                      &
&              +(1/2)*iv12**2*(iv13*iv24**2 + iv14*iv23**2)          &
&              +(1/2)*iv12**2*(iv13**2*iv24 + iv14**2*iv23)          &
&              +(iv12 + iv34)*(iv13**2*iv24**2 + iv14**2*iv23**2)

         Tv13 = (1/2)*( iv12**2*(iv13**2*iv14 + iv14*iv23**2)        &
&              +iv34**2*(iv13**2*iv14 + iv14**2*iv23)                &
&              +(iv12 + iv34)*(2*iv13**3*iv14 + iv14**2*iv23**2) )


         Tv14 = (1/2)*( iv12**2*(iv13*iv14**2 + iv14**2*iv23)        &
&              +iv34**2*(iv13**2*iv14 + iv14**2*iv23)                &
&              +(iv12 + iv34)*(iv13**2*iv14**2 + 2*iv14**3*iv23) )


         Tv23 = (1/2)*( iv12**2*(iv13**2*iv14 + iv14*iv23**2)        &
&              +iv34**2*(iv13*iv14**2 + iv14*iv23**2)                &
&              +(iv12 + iv34)*(iv13**2*iv14**2 + 2*iv14*iv23**3) )


         Tv24 = (1/2)*( iv12**2*(iv13*iv14**2 + iv14**2*iv23)        &
&              +iv34**2*(iv13*iv14**2 + iv14*iv23**2)                &
&              +(iv12 + iv34)*(2*iv13*iv14**3 + iv14**2*iv23**2) )


      Td   = Tv1*v1d+Tv2*v2d
      Tdd  = Tv11*v1d**2 + Tv1*v1dd + Tv22*v2d**2 + Tv2*v2dd + 2*Tv12*v1d*v2d
      Tddk = Tv13*v3d*v1d + Tv14*v4d*v1d + Tv23*v3d*v2d + Tv14*v4d*v2d

         phi = phi + &
&              wi(ca)*wi(cb)*((a*b)**2)*W*(Td + di*Tdd + dk*Tddk)

      End Do
      End Do

      phi = phi * 2*m*(e**4)/(pi**2)

      phi3dD(cdel,cphiD) = phi

   End Do
   End Do
!$OMP END PARALLEL DO

   End Subroutine calPhi3



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



   Subroutine calPhi4(maxD,minD,etai,eta1,na,xi,wi,ndel,nphiD,phi4dD)
     use m_Const_Parameters, only : PAI
!$   use omp_lib
      Implicit None


!************************ Note *********************************
!
!                            Written by Youky Ono
!***************************************************************



!++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
      Integer  ci,cj,ck,ca,cb,na,nb

! Gauss-Legendre integration
      Real(8) xi(na),wi(na)

      Real(8) ddel,dphiD,del,phiD,tmp,phix,phiy,maxD,minD,etai,eta1
      Integer ndel,cdel,nphiD,cphiD,cdelD
      Real(8)  :: phi4dD(-ndel:ndel,-1:nphiD+1)

      Integer  cir,cix,ciy,ciz&
&             ,cjr,cjx,cjy,cjz&
&             ,ckr,ckx,cky,ckz&
&             ,cr ,cx ,cy ,cz,lcx,hcx,lcy,hcy,lcz,hcz
      Real(8)  di,dj,dk,dx,dy,dz,dv,n,ni,nj,nk,r
      Integer  i,j,k

      Real(8)  eta
      Parameter (eta=0.0000000001)

      Real(8)  kF,excLDAd,exLDAd,ecLDAd,pi,Zab,wp,wq,q,m,h,e
      Real(8)  gamma,a,b,C,d
      Parameter (pi=PAI)
      Parameter (e=1.0,m=1.0)
      Parameter (Zab=-0.8491)

      Real(8) T,W,phi,psi,vcik,vcik1,vcik2,vcii,vii
      Real(8) Tv1,Tv2,Tv11,Tv22,Tv12,Tv13,Tv14,Tv23,Tv24,Tdd
      Real(8) v1,v2,v3,v4,v1d,v2d,v3d,v4d,v1dd,v2dd

      Real(8) term1,term2,term3,term4,ns
      Real(8) iv12,iv13,iv14,iv23,iv24,iv34
!++++++++++++++++++++++++ end VARIABLES ++++++++++++++++++++++++
      gamma = 4*pi/9

! phi(del,phiD)
      ddel = 1/Real(ndel+1)
      dphiD = (maxD-minD)/(Real(nphiD))


!$OMP  PARALLEL DO &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(cdel,cphiD,del,phiD,di,dk,phi,ca,cb,a,b,  &
!$OMP&  h,v1,v2,v3,v4,W,T)

   Do cdel = -ndel,ndel
   Do cphiD = -1,nphiD+1

      del = Real(cdel)*ddel
      phiD = minD + Real(cphiD)*dphiD
      
      di = phiD*(1+del)
      dk = phiD*(1-del)

      phi = 0
      Do ca=1,na
      Do cb=1,na
         a = xi(ca)
         b = xi(cb)

         h = 1-Exp(-gamma*((a/di)**2)) + eta
         v1 = (a**2)/(2*h)

         h = 1-Exp(-gamma*((b/di)**2)) + eta
         v2 = (b**2)/(2*h)

         h = 1-Exp(-gamma*((a/dk)**2)) + eta
         v3 = (a**2)/(2*h)

         h = 1-Exp(-gamma*((b/dk)**2)) + eta
         v4 = (b**2)/(2*h)

         W = 2*((3-a*a)*b*Cos(b)*Sin(a) + (3-b*b)*a*Cos(a)*Sin(b)&
&            + (a*a+b*b-3)*Sin(a)*Sin(b)&
&            -3*a*b*Cos(a)*Cos(b))/((a*b)**3)

         T = 0.5*(1/(v1+v2) + 1/(v3+v4) )&
&            *(1/((v1+v3)*(v2+v4)) + 1/((v1+v4)*(v2+v3)))

         phi = phi + &
&              wi(ca)*wi(cb)*((a*b)**2)*W*T

      End Do
      End Do

      phi = phi * 2*m*(e**4)/(pi**2)

      phi4dD(cdel,cphiD) = phi

   End Do
   End Do
!$OMP END PARALLEL DO

   End Subroutine calPhi4



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



   Subroutine d_q0(nx,ny,nz,n,nnx,nny,nnz,q0)
     use m_Const_Parameters, only : PAI
!$   use omp_lib
      Implicit None


!************************ Note *********************************
! This program calculates q0 and d.
!
! Input
!   rho(nrxyz,nsipn) : Total density
!
! Output
!   q0               : 
!
!
!
!                            Written by Youky Ono
!***************************************************************



!++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
      Integer  nx,ny,nz,nspin
      Double Precision n,rs,x,nnx,nny,nnz,nn2,r,term1,term2,term3

      Double precision  kF,exc0,exLDA,excLDA,pi,Zab,q0,m,h,e,GxcLDA,eta
      Parameter (pi=PAI)
      Parameter (e=1.0,m=1.0)
      Parameter (Zab=-0.8491)

!++++++++++++++++++++++++ end VARIABLES ++++++++++++++++++++++++



!---------------------- Calculation Start ----------------------
      nn2 = (nnx**2)+(nny**2)+(nnz**2)

      rs = (3/(4*pi*n))**0.3333333333333333

! Eq.(58), (59) (p.93-94) of Theory of the Inhomogeneous Electron gas.
! S.Lundqvist and N.H.March 1983 Plenum Press, NY
      x = rs/11.4
      GxcLDA = 0.5*((1+x**3)*Log(1+1/x)-x**2+x/2-1/3)
      excLDA = -0.458/rs-0.0666*GxcLDA

      kF = (3*pi*pi*n)**0.3333333333333333
      exLDA = -3*e*kF/(4*pi)

! Eq.(12) of Dion PRL92,246401
      term1 = exLDA*Zab*nn2/(6*kF*n)
      term2 = term1/(6*kF*n)

      exc0 = excLDA - term2


! Eq.(11) of Dion PRL92,246401
      q0 = exc0*kF/exLDA



   End Subroutine d_q0



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



   Subroutine Phi1(maxD,minD,etai,eta1,na,xi,wi,nx,ny,nz,dx,dy,dz,rho,cix,ciy &
&           ,ciz,dxrho,dyrho,dzrho,d2xrho,d2yrho,d2zrho,C6,ndel,nphiD,phi1dD,naPhi1)
     use m_Const_Parameters, only : PAI
!$   use omp_lib
      Implicit None

!************************ Note *********************************
!
!
!                            Written by Youky Ono
!***************************************************************



!++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
      Integer  nx,ny,nz,nspin
      Real(8) da,db,a1,a2
      Integer  ci,cj,ck,ca,cb,na,nb
      Parameter (a1=0,a2=60)

! Gauss-Legendre integration
      Real(8) xi(na),wi(na)

      Real(8) ddel,dphiD,del,phiD,phix,phiy,maxD,minD
      Integer ndel,cdel,nphiD,cphiD
      Real(8) phi1dD(-ndel:ndel,-1:nphiD+1)

      Integer  cir,cix,ciy,ciz&
&             ,cjr,cjx,cjy,cjz&
&             ,ckr,ckx,cky,ckz&
&             ,cr ,cx ,cy ,cz,lcx,hcx,lcy,hcy,lcz,hcz
      Real(8)  di,dj,dk,dx,dy,dz,dv,n,ni,nj,nk
      Integer  i,j,k,zxp,zxm,zyp,zym,zzp,zzm

      Real(8)  eta,etai,eta1
      Parameter (eta=0.0000000001)

      Real(8) rs,r,x,nnx,nny,nnz,nn2,n2nx,n2ny,n2nz
      Real(8) :: rho(nx,ny,nz),Eii(nx,ny,nz)
      Real(8) :: dxrho(nx,ny,nz),d2xrho(nx,ny,nz) &
&               ,dyrho(nx,ny,nz),d2yrho(nx,ny,nz) &
&               ,dzrho(nx,ny,nz),d2zrho(nx,ny,nz)
      Real(8) rhomin
      Parameter (rhomin=1.e-6)

      Real(8)  kF,excLDAd,exLDAd,ecLDAd,pi,Zab,wp,wq,q,m,h,e
      Real(8)  gamma,a,b,C,d,q0,q0i,q0k
      Real(8)  v1,v2,v3,v4
      Parameter (pi=PAI)
      Parameter (e=1.0,m=1.0)
      Parameter (Zab=-0.8491)

      Real(8) :: C6(0:nx-1,0:ny-1,0:nz-1)

      Real(8) T,W,phi,psi,vcik,vcik1,vcik2,vcii,vii
      Real(8) Tv1,Tv2,Tv11,Tv22,Tv12,Tv13,Tv14,Tv23,Tv24
      Real(8) v1d,v2d,v3d,v4dd

      Real(8) naPhi1,alfa1,term1,term2,ns

      integer cD,nD
      parameter (nD=100)
      real(8) minDii,maxDii,LD1,LD2,PLD,LDxi(nD),LDwi(nD)
      parameter (minDii=0)
!++++++++++++++++++++++++ end VARIABLES ++++++++++++++++++++++++



!---------------------- Calculation Start ----------------------
      nb = na
      da=0
      db = da
      gamma = 4*pi/9

      dv = dx*dy*dz
      ni = rho(cix,ciy,ciz)
      kF = (3*pi*pi*ni)**0.3333333333333333
      rs = (3/(4*pi*ni))**0.3333333333333333
      x = rs/11.4

      vcik1=0
      If(ni.GT.rhomin) Then

      nnx = dxrho(cix,ciy,ciz)
      nny = dyrho(cix,ciy,ciz)      
      nnz = dzrho(cix,ciy,ciz)
      nn2 = (nnx**2)+(nny**2)+(nnz**2)

      Call d_q0(nx,ny,nz,ni,nnx,nny,nnz,q0i)

! Calculate the coefficient alfa1.
      n2nx = d2xrho(cix,ciy,ciz)
      n2ny = d2yrho(cix,ciy,ciz)      
      n2nz = d2zrho(cix,ciy,ciz)
      term1 = n2nx+n2ny+n2nz 
      term2 = ((((pi/(3*ni))**0.666666666666667)/kF)+(1/ni))*nn2
      ns = (term1-term2)/(2*kF*ni)
      term1 = (36*pi*(ni**4))**0.33333333333333333
      exLDAd = -0.458/((rs**2)*term1)
      ecLDAd = -0.0333*(3*x*x*Log(1+1/x)-3*x+1.5-1/x)/(11.4*term1)
      excLDAd = exLDAd + ecLDAd
      alfa1 = ((zab*ns/9)+(7*zab*nn2/(3*9*4*kF*ni*ni))-(4*pi*ni*excLDAd/3))/q0i

      lcx = cix - Nint(eta1/dx)
      hcx = cix + Nint(eta1/dx)
      lcy = ciy - Nint(eta1/dy)
      hcy = ciy + Nint(eta1/dy)
      lcz = ciz - Nint(eta1/dz)
      hcz = ciz + Nint(eta1/dz)

      Do ckx = lcx,hcx
      Do cky = lcy,hcy
      Do ckz = lcz,hcz
      r =                              &
&        ((REAL(cix-ckx)*dx)**2        &
&       + (REAL(ciy-cky)*dy)**2        &
&       + (REAL(ciz-ckz)*dz)**2)**0.5

   If(r.Lt.eta1.AND.r.Ge.etai) Then

         cx=MOD(20*nx+ckx-1,nx)+1
         cy=MOD(20*ny+cky-1,ny)+1
         cz=MOD(20*nz+ckz-1,nz)+1
         nk = rho(cx,cy,cz)
      If(nk.GT.rhomin) Then

      nnx = dxrho(cx,cy,cz)
      nny = dyrho(cx,cy,cz)      
      nnz = dzrho(cx,cy,cz)

      Call d_q0(nx,ny,nz,nk,nnx,nny,nnz,q0k)
      di = r*q0i + eta
      dk = r*q0k + eta

      ddel = 1/Real(ndel+1)
      del = ((di-dk)/(di+dk))
      cdel = AINT(del/ddel)
      phiD = (di+dk)/2

      If(phiD.LT.maxD.AND.phiD.GT.minD) then
         cphiD = AINT(REAL(nphiD)*(phiD-minD)/(maxD-minD))

         phix = del/ddel - REAL(cdel)
         phiy = REAL(nphiD)*(phiD-minD)/(maxD-minD) - REAL(cphiD)

         phi =  (1-phix)*(1-phiy) * phi1dD(cdel  ,cphiD  )  &
&              +   phix *(1-phiy) * phi1dD(cdel+1,cphiD  )  &
&              +(1-phix)*   phiy  * phi1dD(cdel  ,cphiD+1)  &
&              +   phix *   phiy  * phi1dD(cdel+1,cphiD+1)
      End If


! vcik1 is the non-local correlation potential 
! v_c^nl summated for (etai < rik < eta1).
      vcik1 = vcik1 + dv*alfa1*di*phi*nk

      End If

   End If

      End Do
      End Do
      End Do
         End If



   vcik2 = 0
   If(ni.GT.rhomin) Then
   C = 12*m*(e**4)*((4*pi/9)**3)

   Do ckx = 1,nx
   Do cky = 1,ny
   Do ckz = 1,nz
   nk = rho(ckx,cky,ckz)
   If(nk.GT.rhomin) Then

      nnx = dxrho(ckx,cky,ckz)
      nny = dyrho(ckx,cky,ckz)      
      nnz = dzrho(ckx,cky,ckz)
      Call d_q0(nx,ny,nz,nk,nnx,nny,nnz,q0k)

      psi = 2*(3*q0i**2+q0k**2)/(q0i**3*q0k**2*(q0i**2+q0k**2)**2)

! vcik2 is the non-local correlation potential 
! v_c^nl summated for ( eta1 < rik ).
      vcik2 = vcik2 + dv*alfa1*C*q0i*psi*nk &
&      *C6(ABS(cix-ckx),ABS(ciy-cky),ABS(ciz-ckz))

   End If

   End Do
   End Do
   End Do

   End If



   vcik = vcik1 + vcik2



   vii=0
   vcii=0
   cjx=cix
   cjy=ciy
   cjz=ciz
      n = rho(cix,ciy,ciz)
   If(n.GT.rhomin) Then

      nnx = dxrho(cjx,cjy,cjz)
      nny = dyrho(cjx,cjy,cjz)      
      nnz = dzrho(cjx,cjy,cjz)
      Call d_q0(nx,ny,nz,n,nnx,nny,nnz,q0k)
      q0 = q0k

      maxDii = q0*etai

      Call gauleg(minDii,maxDii,nD,LDxi,LDwi)
!      Call gauleg(a1,a2,na,xi,wi)
      
      Do cD=1,nD
         LD1 = LDxi(cD)
         LD2 = LD1

         di = LD1
         dk = LD2

      phi = 0

      Do ca=1,na
      Do cb=1,na

         a = xi(ca)
         b = xi(cb)

         h = 1-Exp(-gamma*((a/di)**2)) + eta
         v1 = (a**2)/(2*h)
         v1d = -4*pi*(a**4)*(h-1)/(9*(di**3)*(h**2))

         h = 1-Exp(-gamma*((b/di)**2)) + eta
         v2 = (b**2)/(2*h)
         v2d = -4*pi*(b**4)*(h-1)/(9*(di**3)*(h**2))

         h = 1-Exp(-gamma*((a/dk)**2)) + eta
         v3 = (a**2)/(2*h)

         h = 1-Exp(-gamma*((b/dk)**2)) + eta
         v4 = (b**2)/(2*h)

         W = 2*((3-a*a)*b*Cos(b)*Sin(a) + (3-b*b)*a*Cos(a)*Sin(b)           &
&            + (a*a+b*b-3)*Sin(a)*Sin(b)                                    &
&            -3*a*b*Cos(a)*Cos(b))/((a*b)**3)

         Tv1 = ((1/((v1+v2)**2))*(1/((v1+v3)*(v2+v4))+1/((v1+v4)*(v2+v3)))  &
&              +(1/(v1+v2)+1/(v3+v4))                                       &
&              *(1/(((v1+v3)**2)*(v2+v4))+1/(((v1+v4)**2)*(v2+v3))))*(-0.5)

         Tv2 = ((1/((v1+v2)**2))*(1/((v1+v3)*(v2+v4))+1/((v1+v4)*(v2+v3)))  &
&              +(1/(v1+v2)+1/(v3+v4))                                       &
&              *(1/((v1+v3)*((v2+v4)**2))+1/((v1+v4)*((v2+v3)**2))))*(-0.5)

         phi = phi + &
&              wi(ca)*wi(cb)*((a*b)**2)*W*(Tv1*v1d+Tv2*v2d)

      End Do
      End Do

      phi = phi * 2*m*(e**4)/(pi**2)

      vii = vii + LDwi(cD)*4*pi*(LD1**2)*phi

      End Do


! vcii is the non-local correlation potential 
! v_c^nl summated for ( rik < etai ).
      vcii = alfa1*n*vii/(q0**3)

   End If

   naPhi1 = vcik1 + vcik2 + vcii



! End of Phi1
   End Subroutine Phi1



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



   Subroutine Phi2(maxD,minD,etai,eta1,na,xi,wi,nx,ny,nz,dx,dy,dz,rho,cix,ciy,ciz  &
&                     ,dxrho,dyrho,dzrho,d2xrho,d2yrho,d2zrho               &
&                     ,dxyrho,dxzrho,dyzrho,C6,ndel,nphiD,phi2dD,naPhi2)
     use m_Const_Parameters, only : PAI
!$   use omp_lib
      Implicit None


!************************ Note *********************************
!
!
!                            Written by Youky Ono
!***************************************************************



!++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
      Integer  nx,ny,nz,nspin
      Real(8) da,db,a1,a2
      Integer  ci,cj,ck,ca,cb,na,nb
      Parameter (a1=0,a2=60)

! Gauss-Legendre integration
      Real(8) xi(na),wi(na)

      Real(8) ddel,dphiD,del,phiD,phix,phiy,maxD,minD
      Integer ndel,cdel,nphiD,cphiD
      Real(8) phi2dD(-ndel:ndel,-1:nphiD+1)

      Integer  cir,cix,ciy,ciz&
&             ,cjr,cjx,cjy,cjz&
&             ,ckr,ckx,cky,ckz&
&             ,cr ,cx ,cy ,cz,lcx,hcx,lcy,hcy,lcz,hcz
      Real(8)  di,dj,dk,dx,dy,dz,dv,n,ni,nj,nk
      Integer  i,j,k,zxp,zxm,zyp,zym,zzp,zzm

      Real(8)  eta,etai,eta1
      Parameter (eta=0.0000000001)

      Real(8) rs,r,x,nnx,nny,nnz,nn2,n2nx,n2ny,n2nz,nxy,nxz,nyz
      Real(8) ::      rho(nx,ny,nz),Eii(nx,ny,nz)
      Real(8) ::      dxrho(nx,ny,nz),d2xrho(nx,ny,nz) &
&                    ,dyrho(nx,ny,nz),d2yrho(nx,ny,nz) &
&                    ,dzrho(nx,ny,nz),d2zrho(nx,ny,nz) &
&  ,dxyrho(nx,ny,nz),dxzrho(nx,ny,nz),dyzrho(nx,ny,nz)

      Real(8) rhomin
      Parameter (rhomin=1.e-6)

      Real(8)  kF,GxcLDA,excLDA,excLDAd,exLDAd,ecLDAd,pi,Zab,wp,wq,q,m,h,e
      Real(8)  gamma,a,b,C,d,q0,q0i,q0k,q0dx,q0dy,q0dz
      Parameter (pi=PAI)
      Parameter (e=1.0,m=1.0)
      Parameter (Zab=-0.8491)

      Real(8) :: C6(0:nx-1,0:ny-1,0:nz-1)

      Real(8) T,W,phi,psi,vcik,vcik1,vcik2,vcii,vii
      Real(8) Tv1,Tv2,Tv11,Tv22,Tv12,Tv13,Tv14,Tv23,Tv24,Tdd
      Real(8) v1,v2,v3,v4,v1d,v2d,v3d,v4d,v1dd,v2dd

      Real(8) naPhi2,alfa2,term1,term2,term3,term4,ns
      Real(8) iv12,iv13,iv14,iv23,iv24,iv34

      integer cD,nD
      parameter (nD=100)
      real(8) minDii,maxDii,LD1,LD2,PLD,LDxi(nD),LDwi(nD)
      parameter (minDii=0)
!++++++++++++++++++++++++ end VARIABLES ++++++++++++++++++++++++



!---------------------- Calculation Start ----------------------
      nb = na
      da=0
      db = da
      gamma = 4*pi/9

      dv = dx*dy*dz
      ni = rho(cix,ciy,ciz)
      kF = (3*pi*pi*ni)**0.3333333333333333
      rs = (3/(4*pi*ni))**0.3333333333333333
      x = rs/11.4

      vcik1=0
      If(ni.GT.rhomin) Then

      nnx = dxrho(cix,ciy,ciz)
      nny = dyrho(cix,ciy,ciz)      
      nnz = dzrho(cix,ciy,ciz)
      nn2 = (nnx**2)+(nny**2)+(nnz**2)

      Call d_q0(nx,ny,nz,ni,nnx,nny,nnz,q0i)



!********* Calculate the coefficient alfa2. *********
      n2nx = d2xrho(cix,ciy,ciz)
      n2ny = d2yrho(cix,ciy,ciz)      
      n2nz = d2zrho(cix,ciy,ciz)
      nxy = dxyrho(cix,ciy,ciz)
      nxz = dxzrho(cix,ciy,ciz)
      nyz = dyzrho(cix,ciy,ciz)

! Eq.(58), (59) (p.93-94) of Theory of the Inhomogeneous Electron gas.
! S.Lundqvist and N.H.March 1983 Plenum Press, NY
      GxcLDA = 0.5*((1+x**3)*Log(1+1/x)-x**2+x/2-1/3)
      excLDA = -0.458/rs-0.0666*GxcLDA
      
      term1 = (36*pi*(ni**4))**0.33333333333333333
      exLDAd = -0.458/((rs**2)*term1)
      ecLDAd = -0.0333*(3*x*x*Log(1+1/x)-3*x+1.5-1/x)/(11.4*term1)
      excLDAd = exLDAd + ecLDAd

! q0dx
      term1 = -nnx*(4*pi/3)*(excLDAd*ni+excLDA)
      term2 = nn2*nnx*(((pi*pi*ni/9)**0.33333333333333333)+2*kF)/(-4*(kF**2)*(ni**3))
      term3 = 2*( nnx*n2nx + nny*nxy + nnz*nxz )/(4*kF*(ni**2))
      term4 = -(Zab/9)*(term2+term3)
      q0dx = term1 + term4

! q0dy
      term1 = -nny*(4*pi/3)*(excLDAd*ni+excLDA)
      term2 = nn2*nny*(((pi*pi*ni/9)**0.33333333333333333)+2*kF)/(-4*(kF**2)*(ni**3))
      term3 = 2*( nnx*nxy + nny*n2ny + nnz*nyz )/(4*kF*(ni**2))
      term4 = -(Zab/9)*(term2+term3)
      q0dy = term1 + term4

! q0dz
      term1 = -nnz*(4*pi/3)*(excLDAd*ni+excLDA)
      term2 = nn2*nnz*(((pi*pi*ni/9)**0.33333333333333333)+2*kF)/(-4*(kF**2)*(ni**3))
      term3 = 2*( nnx*nxz + nny*nyz + nnz*n2nz )/(4*kF*(ni**2))
      term4 = -(Zab/9)*(term2+term3)
      q0dz = term1 + term4

      term1 = ( nnx*q0dx + nny*q0dy + nnz*q0dz )/(2*kF*ni)

      alfa2 = (Zab/9)*term1/(q0i**2)
!******** End Calculating coefficient alfa2. ********



      lcx = cix - Nint(eta1/dx)
      hcx = cix + Nint(eta1/dx)
      lcy = ciy - Nint(eta1/dy)
      hcy = ciy + Nint(eta1/dy)
      lcz = ciz - Nint(eta1/dz)
      hcz = ciz + Nint(eta1/dz)

      Do ckx = lcx,hcx
      Do cky = lcy,hcy
      Do ckz = lcz,hcz
      r = (((cix-ckx)*dx)**2 + ((ciy-cky)*dy)**2 + ((ciz-ckz)*dz)**2)**0.5

   If(r.LT.eta1.AND.r.GE.etai) Then

         cx=MOD(20*nx+ckx-1,nx)+1
         cy=MOD(20*ny+cky-1,ny)+1
         cz=MOD(20*nz+ckz-1,nz)+1
         nk = rho(cx,cy,cz)
      If(nk.GT.rhomin) Then

      nnx = dxrho(cx,cy,cz)
      nny = dyrho(cx,cy,cz)      
      nnz = dzrho(cx,cy,cz)

      Call d_q0(nx,ny,nz,nk,nnx,nny,nnz,q0k)
      di = r*q0i + eta
      dk = r*q0k + eta

      ddel = 1/Real(ndel+1)
      del = ((di-dk)/(di+dk))
      cdel = AINT(del/ddel)
      phiD = (di+dk)/2

      If(phiD.LT.maxD.AND.phiD.GT.minD) then
         cphiD = AINT(REAL(nphiD)*(phiD-minD)/(maxD-minD))

         phix = del/ddel - REAL(cdel)
         phiy = REAL(nphiD)*(phiD-minD)/(maxD-minD) - REAL(cphiD)

         phi =  (1-phix)*(1-phiy) * phi2dD(cdel  ,cphiD  )  &
&              +   phix *(1-phiy) * phi2dD(cdel+1,cphiD  )  &
&              +(1-phix)*   phiy  * phi2dD(cdel  ,cphiD+1)  &
&              +   phix *   phiy  * phi2dD(cdel+1,cphiD+1)
      End If

! vcik1 is the non-local correlation potential 
! v_c^nl summated for (etai < rik < eta1).
      vcik1 = vcik1 + dv*alfa2*(di**2)*phi*nk

      End If

   End If

      End Do
      End Do
      End Do
         End If



   vcik2 = 0
   If(ni.GT.rhomin) Then
   C = 12*m*(e**4)*((4*pi/9)**3)

   Do ckx = 1,nx
   Do cky = 1,ny
   Do ckz = 1,nz
   nk = rho(ckx,cky,ckz)
   If(nk.GT.rhomin) Then

      nnx = dxrho(ckx,cky,ckz)
      nny = dyrho(ckx,cky,ckz)      
      nnz = dzrho(ckx,cky,ckz)
      Call d_q0(nx,ny,nz,nk,nnx,nny,nnz,q0k)

      psi = -2*(15*q0i**4+10*(q0i*q0k)**2+3*q0k**4) &
&             /(q0i**4*q0k**2*(q0i**2+q0k**2)**3)

! vcik2 is the non-local correlation potential 
! v_c^nl summated for ( eta1 < rik ).
      vcik2 = vcik2 + dv*alfa2*C*(q0i**2)*psi*nk &
&      *C6(ABS(cix-ckx),ABS(ciy-cky),ABS(ciz-ckz))

   End If

   End Do
   End Do
   End Do

   End If



   vcik = vcik1 + vcik2



   vcii=0
   vii=0
   cjx=cix
   cjy=ciy
   cjz=ciz
      n = rho(cix,ciy,ciz)
   If(n.GT.rhomin) Then

      nnx = dxrho(cjx,cjy,cjz)
      nny = dyrho(cjx,cjy,cjz)      
      nnz = dzrho(cjx,cjy,cjz)
      Call d_q0(nx,ny,nz,n,nnx,nny,nnz,q0k)
      q0 = q0k

      maxDii = q0*etai

      Call gauleg(minDii,maxDii,nD,LDxi,LDwi)
!      Call gauleg(a1,a2,na,xi,wi)
      
      Do cD=1,nD
         LD1 = LDxi(cD)
         LD2 = LD1

         di = LD1
         dk = LD2

      phi = 0

      Do ca=1,na
      Do cb=1,na

         a = xi(ca)
         b = xi(cb)

         h = 1-Exp(-gamma*((a/di)**2)) + eta
         v1 = (a**2)/(2*h)
         v1d = -4*pi*(a**4)*(h-1)/(9*(di**3)*(h**2))
         v1dd = 2*((4*pi*(a**3)/(9*(di**3)))**2)*(1-h)*(2-3*h)/(h**3)

         h = 1-Exp(-gamma*((b/di)**2)) + eta
         v2 = (b**2)/(2*h)
         v2d = -4*pi*(b**4)*(h-1)/(9*(di**3)*(h**2))
         v2dd = 2*((4*pi*(b**3)/(9*(di**3)))**2)*(1-h)*(2-3*h)/(h**3)

         h = 1-Exp(-gamma*((a/dk)**2)) + eta
         v3 = (a**2)/(2*h)

         h = 1-Exp(-gamma*((b/dk)**2)) + eta
         v4 = (b**2)/(2*h)

         W = 2*((3-a*a)*b*Cos(b)*Sin(a) + (3-b*b)*a*Cos(a)*Sin(b)           &
&            + (a*a+b*b-3)*Sin(a)*Sin(b)                                    &
&            -3*a*b*Cos(a)*Cos(b))/((a*b)**3)


         iv12 = 1/(v1+v2)
         iv13 = 1/(v1+v3)
         iv14 = 1/(v1+v4)
         iv23 = 1/(v2+v3)
         iv24 = 1/(v2+v4)
         iv34 = 1/(v3+v4)


         Tv1 = -(1/2)*(iv12**2)*(iv13*iv24+iv14*iv23)              &
&              -(1/2)*(iv12+iv34)*(iv13**2*iv24+iv14**2*iv23)


         Tv2 = -(1/2)*(iv12**2)*(iv13*iv24+iv14*iv23)              &
&              -(1/2)*(iv12+iv34)*(iv13*iv24**2+iv14*iv23**2)


         Tv11 = iv12**3*(iv13*iv24+iv14*iv23)                      &
&              +iv12**2*(iv13**2*iv24+iv14**2*iv23)                &
&              +(iv12+iv34)*(iv13**3*iv24+iv14**3*iv23)


         Tv22 = iv12**3*(iv13*iv24+iv14*iv23)                      &
&              +iv12**2*(iv13*iv24**2+iv14*iv23**2)                &
&              +(iv12+iv34)*(iv13*iv24**3+iv14*iv23**3)


         Tv12 = iv12**3*(iv13*iv24+iv14*iv23)                      &
&              +(1/2)*iv12**2*(iv13*iv24**2+iv14*iv23**2)          &
&              +(1/2)*iv12**2*(iv13**2*iv24+iv14**2*iv23)          &
&              +(iv12+iv34)*(iv13**2*iv24**2+iv14**2*iv23**2)


         Tdd = Tv11*v1d**2 + Tv1*v1dd + Tv22*v2d**2 + Tv2*v2dd + 2*Tv12*v1d*v2d
         
         phi = phi + &
&              wi(ca)*wi(cb)*((a*b)**2)*W*Tdd

      End Do
      End Do

      phi = phi * 2*m*(e**4)/(pi**2)

      vii = vii + LDwi(cD)*4*pi*(LD1**2)*phi

      End Do


! vcii is the non-local correlation potential 
! v_c^nl summated for ( rik < etai ).
      vcii = alfa2*n*vii/(q0**3)

   End If

   naPhi2 = vcik1 + vcik2 + vcii



! End of Phi2
   End Subroutine Phi2



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



   Subroutine Phi3(maxD,minD,etai,eta1,na,xi,wi,nx,ny,nz,dx,dy,dz,rho,cix,ciy,ciz  &
&                     ,dxrho,dyrho,dzrho,d2xrho,d2yrho,d2zrho,C6,ndel,nphiD,phi3dD,naPhi3)
     use m_Const_Parameters, only : PAI
!$   use omp_lib
      Implicit None


!************************ Note *********************************
!
!
!                            Written by Youky Ono
!***************************************************************



!++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
      Integer  nx,ny,nz,nspin
      Real(8) da,db
      Integer  ci,cj,ck,ca,cb,na,nb

! Gauss-Legendre integration
      Real(8) xi(na),wi(na)

      Real(8) ddel,dphiD,del,phiD,phix,phiy,maxD,minD
      Integer ndel,cdel,nphiD,cphiD
      Real(8) phi3dD(-ndel:ndel,-1:nphiD+1)

      Integer  cir,cix,ciy,ciz&
&             ,cjr,cjx,cjy,cjz&
&             ,ckr,ckx,cky,ckz&
&             ,cr ,cx ,cy ,cz,lcx,hcx,lcy,hcy,lcz,hcz
      Real(8)  di,dj,dk,dx,dy,dz,dv,n,ni,nj,nk
      Integer  i,j,k,zxp,zxm,zyp,zym,zzp,zzm

      Real(8)  eta,etai,eta1
      Parameter (eta=0.0000000001)

      Real(8) rs,r,x,nnx,nny,nnz,nn2,n2nx,n2ny,n2nz,nxy,nxz,nyz
      Real(8) ::      rho(nx,ny,nz),Eii(nx,ny,nz)
      Real(8) ::      dxrho(nx,ny,nz),d2xrho(nx,ny,nz) &
&                    ,dyrho(nx,ny,nz),d2yrho(nx,ny,nz) &
&                    ,dzrho(nx,ny,nz),d2zrho(nx,ny,nz)

      Real(8) rhomin
      Parameter (rhomin=1.e-6)

      Real(8)  kF,GxcLDA,excLDA,excLDAd,exLDAd,ecLDAd,pi,Zab,wp,wq,q,m,h,e
      Real(8)  gamma,a,b,C,d,q0,q0i,q0k,q0dx,q0dy,q0dz
      Parameter (pi=PAI)
      Parameter (e=1.0,m=1.0)
      Parameter (Zab=-0.8491)

      Real(8) :: C6(0:nx-1,0:ny-1,0:nz-1)

      Real(8) T,W,phi,psi,vcik,vcik1,vcik2,vcii,vii
      Real(8) Tv1,Tv2,Tv11,Tv22,Tv12,Tv13,Tv14,Tv23,Tv24,Td,Tdd,Tddk
      Real(8) v1,v2,v3,v4,v1d,v2d,v3d,v4d,v1dd,v2dd

      Real(8) naPhi3,alfa3,term1,term2,term3,term4,ns
      Real(8) iv12,iv13,iv14,iv23,iv24,iv34
!++++++++++++++++++++++++ end VARIABLES ++++++++++++++++++++++++



!---------------------- Calculation Start ----------------------
      nb = na
      da=0
      db = da
      gamma = 4*pi/9

      dv = dx*dy*dz
      ni = rho(cix,ciy,ciz)
      kF = (3*pi*pi*ni)**0.3333333333333333
      rs = (3/(4*pi*ni))**0.3333333333333333
      x = rs/11.4

      vcik1=0
      If(ni.GT.rhomin) Then

      nnx = dxrho(cix,ciy,ciz)
      nny = dyrho(cix,ciy,ciz)      
      nnz = dzrho(cix,ciy,ciz)
      nn2 = (nnx**2)+(nny**2)+(nnz**2)

      Call d_q0(nx,ny,nz,ni,nnx,nny,nnz,q0i)

      lcx = cix - Nint(eta1/dx)
      hcx = cix + Nint(eta1/dx)
      lcy = ciy - Nint(eta1/dy)
      hcy = ciy + Nint(eta1/dy)
      lcz = ciz - Nint(eta1/dz)
      hcz = ciz + Nint(eta1/dz)



      Do ckx = lcx,hcx
      Do cky = lcy,hcy
      Do ckz = lcz,hcz
      r =                              &
&        ((REAL(cix-ckx)*dx)**2        &
&       + (REAL(ciy-cky)*dy)**2        &
&       + (REAL(ciz-ckz)*dz)**2)**0.5

   If(r.LT.eta1.AND.r.GE.etai) Then


!********* Calculate the coefficient alfa3. *********
      term1 = (REAL(cix-ckx)*dx)*nnx
      term2 = (REAL(ciy-cky)*dy)*nny
      term3 = (REAL(ciz-ckz)*dz)*nnz
      alfa3 = (Zab/9)*(1/(2*kF*ni))*(term1+term2+term3)
!******** End Calculating coefficient alfa3. ********


         cx=MOD(20*nx+ckx-1,nx)+1
         cy=MOD(20*ny+cky-1,ny)+1
         cz=MOD(20*nz+ckz-1,nz)+1
         nk = rho(cx,cy,cz)
      If(nk.GT.rhomin) Then

      nnx = dxrho(cx,cy,cz)
      nny = dyrho(cx,cy,cz)      
      nnz = dzrho(cx,cy,cz)

      Call d_q0(nx,ny,nz,nk,nnx,nny,nnz,q0k)
      di = r*q0i + eta
      dk = r*q0k + eta

      ddel = 1/Real(ndel+1)
      del = ((di-dk)/(di+dk))
      cdel = AINT(del/ddel)
      phiD = (di+dk)/2

      If(phiD.LT.maxD.AND.phiD.GT.minD) then
         cphiD = AINT(REAL(nphiD)*(phiD-minD)/(maxD-minD))

         phix = del/ddel - REAL(cdel)
         phiy = REAL(nphiD)*(phiD-minD)/(maxD-minD) - REAL(cphiD)

         phi =  (1-phix)*(1-phiy) * phi3dD(cdel  ,cphiD  )  &
&              +   phix *(1-phiy) * phi3dD(cdel+1,cphiD  )  &
&              +(1-phix)*   phiy  * phi3dD(cdel  ,cphiD+1)  &
&              +   phix *   phiy  * phi3dD(cdel+1,cphiD+1)
      End If

! vcik1 is the non-local correlation potential 
! v_c^nl summated for (etai < rik < eta1).
      vcik1 = vcik1 + dv*alfa3*phi*nk

      End If

   End If

      End Do
      End Do
      End Do
         End If



   vcik2 = 0

   vcik = vcik1 + vcik2

   vcii=0

   naPhi3 = vcik1 + vcik2 + vcii



! End of Phi3
   End Subroutine Phi3



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



   Subroutine Phi4(maxD,minD,etai,eta1,na,xi,wi,nx,ny,nz,dx,dy,dz,rho,cix,ciy,ciz  &
&                  ,dxrho,dyrho,dzrho,d2xrho,d2yrho,d2zrho,C6,ndel,nphiD,phi4dD,naPhi4)
     use m_Const_Parameters, only : PAI
!$   use omp_lib
      Implicit None

!************************ Note *********************************
!
!
!                            Written by Youky Ono
!***************************************************************



!++++++++++++++++++++++++++ VARIABLES ++++++++++++++++++++++++++
      Integer  nx,ny,nz,nspin
      Real(8) da,db,a1,a2
      Integer  ci,cj,ck,ca,cb,na,nb
      Parameter (a1=0,a2=60)

! Gauss-Legendre integration
      Real(8) xi(na),wi(na)

      Real(8) ddel,dphiD,del,phiD,phix,phiy,maxD,minD
      Integer ndel,cdel,nphiD,cphiD
      Real(8) :: phi4dD(-ndel:ndel,-1:nphiD+1)

      Integer  cir,cix,ciy,ciz&
&             ,cjr,cjx,cjy,cjz&
&             ,ckr,ckx,cky,ckz&
&             ,cr ,cx ,cy ,cz,lcx,hcx,lcy,hcy,lcz,hcz
      Real(8)  di,dj,dk,dx,dy,dz,dv,n,ni,nj,nk
      Integer  i,j,k,zxp,zxm,zyp,zym,zzp,zzm

      Real(8)  eta,etai,eta1
      Parameter (eta=0.0000000001)

      Real(8) rs,x,nnx,nny,nnz,nn2,r,zo
      Real(8) :: rho(nx,ny,nz),Eii(nx,ny,nz)
      Real(8) :: dxrho(nx,ny,nz),d2xrho(nx,ny,nz) &
&               ,dyrho(nx,ny,nz),d2yrho(nx,ny,nz) &
&               ,dzrho(nx,ny,nz),d2zrho(nx,ny,nz)
      Real(8) rhomin
      Parameter (rhomin=1.e-6)

      Real(8)  kF,exc0,exLDA,excLDA,pi,Zab,wp,wq,q,m,h,e,GxcLDA
      Real(8)  gamma,a,b,C,d,q0,q0i,q0k
      Real(8)  v1,v2,v3,v4
      Parameter (pi=PAI)
      Parameter (e=1.0,m=1.0)
      Parameter (Zab=-0.8491)

      Real(8) :: C6(0:nx-1,0:ny-1,0:nz-1)

      Real(8) T,W,phi,psi,vcik,vcik1,vcik2,vcii,vii

      Real(8) naPhi4,alfa4

      integer cD,nD
      parameter (nD=100)
      real(8) minDii,maxDii,LD1,LD2,PLD,LDxi(nD),LDwi(nD)
      parameter (minDii=0)
!++++++++++++++++++++++++ end VARIABLES ++++++++++++++++++++++++



!---------------------- Calculation Start ----------------------

      alfa4=1
      nb = na
      da=0
      db = da
      gamma = 4*pi/9

      dv = dx*dy*dz
      ni = rho(cix,ciy,ciz)

      vcik1=0
      If(ni.GT.rhomin) Then

      nnx = dxrho(cix,ciy,ciz)
      nny = dyrho(cix,ciy,ciz)      
      nnz = dzrho(cix,ciy,ciz)
      Call d_q0(nx,ny,nz,ni,nnx,nny,nnz,q0i)

      lcx = cix - Nint(eta1/dx)
      hcx = cix + Nint(eta1/dx)
      lcy = ciy - Nint(eta1/dy)
      hcy = ciy + Nint(eta1/dy)
      lcz = ciz - Nint(eta1/dz)
      hcz = ciz + Nint(eta1/dz)

!      Do ckx = cix,cix ! lcx,hcx
!      Do cky = ciy,ciy ! lcy,hcy
!      Do ckz = ciz,ciz ! lcz,hcz
      Do ckx = lcx,hcx
      Do cky = lcy,hcy
      Do ckz = lcz,hcz
      r = (((cix-ckx)*dx)**2 + ((ciy-cky)*dy)**2 + ((ciz-ckz)*dz)**2)**0.5

   If(r.LT.eta1.AND.r.GE.etai) Then

         cx=MOD(20*nx+ckx-1,nx)+1
         cy=MOD(20*ny+cky-1,ny)+1
         cz=MOD(20*nz+ckz-1,nz)+1
         nk = rho(cx,cy,cz)
      If(nk.GT.rhomin) Then

      nnx = dxrho(cx,cy,cz)
      nny = dyrho(cx,cy,cz)      
      nnz = dzrho(cx,cy,cz)
      Call d_q0(nx,ny,nz,nk,nnx,nny,nnz,q0k)
      di = r*q0i + eta
      dk = r*q0k + eta

      ddel = 1/Real(ndel+1)
      del = ((di-dk)/(di+dk))
      cdel = AINT(del/ddel)
      phiD = (di+dk)/2

      If(phiD.LT.maxD.AND.phiD.GT.minD) then
         cphiD = AINT(REAL(nphiD)*(phiD-minD)/(maxD-minD))

         phix = del/ddel - REAL(cdel)
         phiy = REAL(nphiD)*(phiD-minD)/(maxD-minD) - REAL(cphiD)

         phi =  (1-phix)*(1-phiy) * phi4dD(cdel  ,cphiD  )  &
&              +   phix *(1-phiy) * phi4dD(cdel+1,cphiD  )  &
&              +(1-phix)*   phiy  * phi4dD(cdel  ,cphiD+1)  &
&              +   phix *   phiy  * phi4dD(cdel+1,cphiD+1)
      End If


! vcik1 is the non-local correlation potential 
! v_c^nl summated for (etai < rik < eta1).
      vcik1 = vcik1 + dv*alfa4*phi*nk

      End If

   End If

      End Do
      End Do
      End Do
         End If



   vcik2 = 0
   If(ni.GT.rhomin) Then
   C = 12*m*(e**4)*((4*pi/9)**3)

   Do ckx = 1,nx
   Do cky = 1,ny
   Do ckz = 1,nz
   nk = rho(ckx,cky,ckz)
   If(nk.GT.rhomin) Then

      nnx = dxrho(ckx,cky,ckz)
      nny = dyrho(ckx,cky,ckz)      
      nnz = dzrho(ckx,cky,ckz)
      Call d_q0(nx,ny,nz,nk,nnx,nny,nnz,q0k)

      psi = (((q0i*q0k)**2)*(q0i**2+q0k**2))**(-1)

! vcik2 is the non-local correlation potential 
! v_c^nl summated for ( eta1 < rik ).
      vcik2 = vcik2 -1*dv*C*psi*nk &
&      *C6(ABS(cix-ckx),ABS(ciy-cky),ABS(ciz-ckz))

   End If

   End Do
   End Do
   End Do

   End If



!   vcik = vcik1 + vcik2



   vcii = 0
   vii=0
   cjx=cix
   cjy=ciy
   cjz=ciz
      n = rho(cix,ciy,ciz)
   If(n.GT.rhomin) Then

      nnx = dxrho(cjx,cjy,cjz)
      nny = dyrho(cjx,cjy,cjz)      
      nnz = dzrho(cjx,cjy,cjz)
      Call d_q0(nx,ny,nz,n,nnx,nny,nnz,q0k)
      q0 = q0k

      maxDii = q0*etai

      Call gauleg(minDii,maxDii,nD,LDxi,LDwi)
!      Call gauleg(a1,a2,na,xi,wi)
      
      Do cD=1,nD
         LD1 = LDxi(cD)
         LD2 = LD1

         di = LD1
         dk = LD2

      phi = 0
!      Do ca=1,nD
!      Do cb=1,nD

!         a = LDxi(ca)
!         b = LDxi(cb)

      Do ca=1,na
      Do cb=1,na

         a = xi(ca)
         b = xi(cb)

         h = 1-Exp(-gamma*((a/di)**2)) + eta
         v1 = (a**2)/(2*h)

         h = 1-Exp(-gamma*((b/di)**2)) + eta
         v2 = (b**2)/(2*h)

         h = 1-Exp(-gamma*((a/dk)**2)) + eta
         v3 = (a**2)/(2*h)

         h = 1-Exp(-gamma*((b/dk)**2)) + eta
         v4 = (b**2)/(2*h)

         W = 2*((3-a*a)*b*Cos(b)*Sin(a) + (3-b*b)*a*Cos(a)*Sin(b)&
&            + (a*a+b*b-3)*Sin(a)*Sin(b)&
&            -3*a*b*Cos(a)*Cos(b))/((a*b)**3)

         T = 0.5*(1/(v1+v2) + 1/(v3+v4) )&
&            *(1/((v1+v3)*(v2+v4)) + 1/((v1+v4)*(v2+v3)))

         phi = phi + &
&              wi(ca)*wi(cb)*((a*b)**2)*W*T

      End Do
      End Do

      phi = phi * 2*m*(e**4)/(pi**2)

      vii = vii + LDwi(cD)*4*pi*(LD1**2)*phi

      End Do


! vcii is the non-local correlation potential 
! v_c^nl summated for ( rik < etai ).
      vcii = 2*alfa4*n*vii/(q0**3)

   End If

   naPhi4 = vcik1 + vcik2 + vcii



! End of Phi4
   End Subroutine Phi4


