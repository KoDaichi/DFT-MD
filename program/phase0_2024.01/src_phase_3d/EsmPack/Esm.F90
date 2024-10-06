! Copyright (c) 2012, Minoru Otani <minoru.otani@aist.go.jp> 
! 
! Permission is hereby granted, free of charge, to any person 
! obtaining a copy of this software and associated documentation 
! files (the "Software"), to deal in the Software without restriction, 
! including without limitation the rights to use, copy, modify, merge, 
! publish, distribute, sublicense, and/or sell copies of the Software, 
! and to permit persons to whom the Software is furnished to do so, 
! subject to the following conditions:
 
! The above copyright notice and this permission notice shall be 
! included in all copies or substantial portions of the Software.
 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
! DEALINGS IN THE SOFTWARE.

Module ESM_VARS
  !
  ! constants
  !
  Integer, Parameter :: gstart = 2
  Real(8), Parameter :: dual = 8d0

  Integer, Parameter :: esm_nfit = 4
!  Real(8), Parameter :: esm_w = 0d0
!  Real(8), Parameter :: esm_efield = 0d0
  Real(8)            :: esm_w = 0d0
  Real(8)            :: esm_efield = 0d0

  REAL(8), Parameter :: pi  = 3.141592653589793238d0
  REAL(8), Parameter :: tpi = 2d0 * pi
  REAL(8), Parameter :: fpi = 4d0 * pi
  Real(8), Parameter :: eps8 = 1d-8

  !
  ! Parameter from input file
  !
  Integer              :: nat              ! Number of Atoms
  Real(8), Allocatable :: tau(:,:)         ! Coordinates of Atoms
  Real(8)              :: at(3,3)          ! Lattice Vectors
  Real(8)              :: alat             ! Lattice Parameter
  Integer              :: nr1x, nr2x, nr3x ! FFT Grid 
  CHARACTER (LEN=3)    :: esm_bc           ! Boundary Condition for ESM
  Logical              :: gamma_only       ! gamma only

  Real(8)              :: bg(3,3)          ! Reciprocal Lattice Vectors
  Real(8)              :: omega
  Real(8), Allocatable :: upf_zp(:)
  Real(8)              :: tpiba2, ecutwfc, gcutm
  Integer              :: nrxx

  Integer                 :: ngm, nspin
  Integer,    Allocatable :: mill(:,:), nl(:), nlm(:)
  Complex(8), Allocatable :: rhog_(:,:)
  !
  INTEGER              :: ngm_2d = 0
  INTEGER, Allocatable :: mill_2d(:,:), imill_2d(:,:)
  !

  Real(8)              :: gew=-1

  COMPLEX(8), Allocatable :: vbar(:)
  Real(8)                 :: z_wall
  Real(8)                 :: bar_height
  Real(8)                 :: bar_width
  Integer                 :: izwall=0
  INTEGER              :: communicator
  INTEGER              :: npes,mype

End Module ESM_VARS

Subroutine SetupEsmVars( inputfile )
  Use ESM_VARS
  Implicit None
  !
  Character(256), Intent(In) :: inputfile
  !
  Integer :: ig, n1, n2, n3, i
  Real(8) :: at12(3)
  Character(256) :: key
  Real(8), External :: InnerProd
  Character(2) :: atom

!
! read from input file
!
  gamma_only = .False.
  Open(10, file=inputfile, status='old')
  Do While( .True. )
     Read(10,'(A256)',end=999) key
     If( key == "nat" ) Then
        Read(10,*) nat
        If( .Not. Allocated(tau)    ) Allocate( tau(3,nat)  )
        If( .Not. Allocated(upf_zp) ) Allocate( upf_zp(nat) )
     Else If( key == "positions" ) Then
        Do i = 1, nat
           Read(10,*) atom, tau(1:3,i)
           If( atom == "H" ) upf_zp(i) = 1.0
           If( atom == "O" ) upf_zp(i) = 6.0
        End Do
     Else If( key == "lattice" ) Then
        Do i = 1, 3
           Read(10,*) at(1:3,i)
        End Do
     Else If( key == "alat" ) Then
        Read(10,*) alat
     Else If( key == "grids" ) Then
        Read(10,*) nr1x, nr2x, nr3x
     Else If( key == "ecutwfc" ) Then
        Read(10,*) ecutwfc
     Else If( key == "esm_bc" ) Then
        Read(10,*) esm_bc
     Else If( key == "gamma" ) Then
        gamma_only = .True.
     End If
  End Do
999 continue
  Close(10)

! calculate parameters
  Call OuterProd( at(1,1), at(1,2), at12 )
  omega = InnerProd( at12, at(1,3) ) * ( alat ** 3 )

  tpiba2 = ( tpi / alat ) ** 2
  gcutm = dual * ecutwfc / tpiba2

  Call MakeRecLatVec( at, bg )

  nrxx = nr1x * nr2x * nr3x

  Return
End Subroutine SetupEsmVars

Subroutine DeallocEsmVars()
  Use ESM_VARS
  Implicit None
  
  If( Allocated( tau      ) ) deAllocate( tau      )
  If( Allocated( upf_zp   ) ) deAllocate( upf_zp   )
  If( Allocated( mill     ) ) deAllocate( mill     )
  If( Allocated( nl       ) ) deAllocate( nl       )
  If( Allocated( nlm      ) ) deAllocate( nlm      )
  If( Allocated( rhog_    ) ) deAllocate( rhog_    )
  If( Allocated( mill_2d  ) ) deAllocate( mill_2d  )
  If( Allocated( imill_2d ) ) deAllocate( imill_2d )

  Return
End Subroutine DeallocEsmVars


Subroutine ReadRhog( rhogfile )
  Use ESM_VARS
  Implicit None
  !
  Character(256), Intent(In) :: rhogfile
  !
  Integer :: n1, n2, n3, ig
  
  Open(11,file=rhogfile,status='old')
  Read(11,*) ngm, nspin
  Allocate( rhog_(ngm,nspin), mill(3,ngm), nl(ngm), nlm(ngm) )
  Do ig = 1, ngm
     If( nspin==1 ) Then
        Read(11,*) n1,n2,n3,nl(ig),nlm(ig),rhog_(ig,1)
     Else
        Read(11,*) n1,n2,n3,nl(ig),nlm(ig),rhog_(ig,1:2)
     End If
     mill(1,ig) = n1
     mill(2,ig) = n2
     mill(3,ig) = n3
  End Do
  Close(11)
  Return
End Subroutine ReadRhog

SUBROUTINE esm_ggen_2d()
  Use ESM_VARS
  IMPLICIT NONE
  !
  INTEGER :: n1xh, n2xh, ng, n1, n2, ng_2d
  LOGICAL, ALLOCATABLE :: do_mill_2d(:,:)
  Real(8) :: tt,t(3)
  !
  !     Make g parallel array
  !
  n1xh = nr1x / 2
  n2xh = nr2x / 2
  ALLOCATE( do_mill_2d(-n1xh:n1xh,-n2xh:n2xh) )
  do_mill_2d(:,:) = .false.

  Do ng = 1, ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     do_mill_2d(n1,n2) = .true.
  End Do
  ngm_2d = COUNT( do_mill_2d )
  
  ALLOCATE( mill_2d(2,ngm_2d), imill_2d(-n1xh:n1xh,-n2xh:n2xh) )
  mill_2d(:,:) = 0
  imill_2d(:,:) = 0
  ng_2d = 1
  DO n1 = -n1xh, n1xh
  DO n2 = -n2xh, n2xh
     IF( do_mill_2d(n1,n2) ) THEN
        mill_2d(1,ng_2d) = n1
        mill_2d(2,ng_2d) = n2
        imill_2d(n1,n2) = ng_2d
        ng_2d = ng_2d + 1
     ENDIF
  ENDDO
  ENDDO
  DEALLOCATE(do_mill_2d)  

  allocate(vbar(nr3x));vbar(:) = (0.d0,0.d0)

  RETURN
END SUBROUTINE esm_ggen_2d

!-----------------------------------------------------------------------
!--------------ESM HARTREE SUBROUTINE-----------------------------------
!-----------------------------------------------------------------------
SUBROUTINE esm_hartree (rhog, ehart, aux )
  Use ESM_VARS
  IMPLICIT NONE
!#ifdef __MPI__
!  include 'mpif.h'
!#endif
  !
  COMPLEX(8) :: rhog(ngm,nspin)   !  n(G)
  REAL(8),    Intent(Out) :: ehart             !  Hartree energy
  COMPLEX(8), Intent(Out) :: aux(nrxx)         !  v_h(G)
  !
  !    here the local variables
  !
  real(8)                 :: tt, t(2), zz, gz, z0, gp, gp2, z1, kn, cc, ss, z, L, &
       z_l, z_r, eh
  integer                  :: ipol, k, k1, k2, k3, iz, ng, n1, n2, n3, &
       nz_r, nz_l, ng_2d
  complex(8),allocatable  :: rhog3(:,:), vg2(:), vg2_in(:), vg3(:,:),vg3_mpi(:,:,:)
  complex(8)              :: xc, ci, tmp, tmp1, tmp2, tmp3, tmp4, f1, f2, f3, f4, &
       a0, a1, a2, a3, c_r, c_l, s_r, s_l, rg3
  integer :: ierr,n1h,n2h
  allocate(rhog3(nr3x,ngm_2d))
  !
  ! Map to FFT mesh (nr3x,ngm_2d)
  rhog3(:,:)=(0.d0,0.d0)
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng)+1
     IF (n3<1) n3 = n3 + nr3x
     if (nspin == 2) then
        rg3 = rhog(ng,1)+rhog(ng,2)
     else
        rg3 = rhog(ng,1)
     endif
     rhog3(n3,ng_2d)=rg3
     if ( gamma_only .and. n1==0 .and. n2==0 ) then
        n3 = -mill(3,ng)+1
        IF (n3<1) n3 = n3 + nr3x
        rhog3(n3,ng_2d)=CONJG(rg3)
     endif
  enddo

  ! End mapping
  !
  allocate(vg3(nr3x,ngm_2d))
  vg3(:,:)=(0.d0,0.d0)
  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+abs(esm_w)
  ci=(0.d0,1.d0)

!****For gp!=0 case ********************
!$omp parallel private( k1, k2, gp2, ipol, t, gp, tmp1, tmp2, vg2, iz, kn, &
!$omp                   cc, ss, tmp, vg2_in, k3, z, rg3 )
  allocate(vg2(nr3x),vg2_in(nr3x))
!$omp do
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle
     t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
     gp2 = sum( t(:) * t(:) ) * tpiba2
     gp=sqrt(gp2)
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0)
     vg2(:)=(0.d0,0.d0)
     do iz=1, nr3x
        if(iz<=nr3x/2) kn=dble(iz-1)     * tpi/L
        if(iz> nr3x/2) kn=dble(iz-1-nr3x) * tpi/L
        cc=cos(kn*z0)
        ss=sin(kn*z0)
        rg3=rhog3(iz,ng_2d)
        vg2(iz)=fpi*rg3/(gp**2+kn**2)
        if (esm_bc.eq.'bc1') then
           tmp1=tmp1+rg3*(cc+ci*ss)/(gp-ci*kn)
           tmp2=tmp2+rg3*(cc-ci*ss)/(gp+ci*kn)
        else if (esm_bc.eq.'bc2') then
           tmp=((gp+ci*kn)*exp(gp*(z1-z0))+(gp-ci*kn)*exp(-gp*(z1-z0)))/(2.d0*gp)
           tmp1=tmp1+rg3*(cc+ci*ss)/(gp**2+kn**2)*tmp
           tmp=((gp-ci*kn)*exp(gp*(z1-z0))+(gp+ci*kn)*exp(-gp*(z1-z0)))/(2.d0*gp)
           tmp2=tmp2+rg3*(cc-ci*ss)/(gp**2+kn**2)*tmp
        else if (esm_bc.eq.'bc3') then
           tmp=((gp+ci*kn)*exp(gp*(z1-z0))+(gp-ci*kn)*exp(-gp*(z1-z0)))/(2.d0*gp)
           tmp1=tmp1+rg3*(cc+ci*ss)/(gp**2+kn**2)*tmp
           tmp=(gp-ci*kn)/gp
           tmp2=tmp2+rg3*(cc-ci*ss)/(gp**2+kn**2)*tmp
        endif
     enddo

     vg2_in(1:nr3x)=vg2(1:nr3x)  ! Since cft_1z is not in-place
     call cft_1z(vg2_in,1,nr3x,nr3x,1,vg2)

     do iz=1,nr3x
        k3=iz-1
        if (k3.gt.nr3x/2) k3=iz-nr3x-1
        z=dble(k3)/dble(nr3x)*L
        if (esm_bc.eq.'bc1') then
           vg2(iz)=vg2(iz)-2.d0*pi/gp*(exp(gp*(z-z0))*tmp1+exp(-gp*(z+z0))*tmp2)
        else if (esm_bc.eq.'bc2') then
           vg2(iz)=vg2(iz)-fpi*(exp(gp*(z-z1))-exp(-gp*(z+3.d0*z1)))*tmp1 &
                         /(1.d0-exp(-4.d0*gp*z1)) &
                          +fpi*(exp(gp*(z-3.d0*z1))-exp(-gp*(z+z1)))*tmp2 &
                         /(1.d0-exp(-4.d0*gp*z1))
        else if (esm_bc.eq.'bc3') then
           vg2(iz)=vg2(iz)-fpi*exp(gp*(z-z1))*tmp1 &
                +tpi*(exp(gp*(z-z0-2.d0*z1))-exp(-gp*(z+z0)))*tmp2
        endif
     enddo
        
     vg2_in(1:nr3x)=vg2(1:nr3x)  ! Since cft_1z is not in-place
     call cft_1z(vg2_in,1,nr3x,nr3x,-1,vg2)
        
     vg3(1:nr3x,ng_2d)=vg2(1:nr3x)*2.d0
  enddo
!#ifdef __MPI__
!  n1h=nr1x/2;n2h=nr2x/2
!  allocate(vg3_mpi(1:nr3x,-n2h:n2h,-n1h:n1h));vg3_mpi(:,:,:)=(0.d0,0.d0)
!  do ng_2d=1,ngm_2d
!     k1 = mill_2d(1,ng_2d)
!     k2 = mill_2d(2,ng_2d)
!     vg3_mpi(:,k2,k1) = vg3(:,ng_2d)
!  enddo
!  call mpi_allreduce(MPI_IN_PLACE,vg3_mpi,nr3x*nr2x*nr1x,mpi_double_complex,mpi_sum,communicator,ierr)
!  do ng_2d=1,ngm_2d
!     k1 = mill_2d(1,ng_2d)
!     k2 = mill_2d(2,ng_2d)
!     vg3(:,ng_2d) = vg3_mpi(:,k2,k1)
!  enddo
!  deallocate(vg3_mpi)
!#endif
  deallocate(vg2,vg2_in)
!$omp end parallel  

!****For gp=0 case ********************
  ng_2d = imill_2d(0,0)
  if( ng_2d > 0 ) then
     allocate(vg2(nr3x),vg2_in(nr3x))
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0); tmp3=(0.d0,0.d0); tmp4=(0.d0,0.d0)
     !for smoothing
     f1=(0.d0,0.d0); f2=(0.d0,0.d0); f3=(0.d0,0.d0); f4=(0.d0,0.d0)
     nz_l=nr3x/2+1+esm_nfit
     nz_r=nr3x/2+1-esm_nfit
     z_l=dble(nz_l-1)*L/dble(nr3x)-L
     z_r=dble(nz_r-1)*L/dble(nr3x)
     !
     rg3=rhog3(1,ng_2d)
     if (esm_bc.eq.'bc1') then
        vg2(1)=-tpi*z0**2*rg3
     else if (esm_bc.eq.'bc2') then
        vg2(1)= tpi*(2.d0*z1-z0)*z0*rg3
     else if (esm_bc.eq.'bc3') then
        vg2(1)= tpi*(4.d0*z1-z0)*z0*rg3
     endif
     do iz=2,nr3x
        if(iz<=nr3x/2) kn=dble(iz-1)     *tpi/L
        if(iz> nr3x/2) kn=dble(iz-1-nr3x) *tpi/L
        cc=cos(kn*z0)
        ss=sin(kn*z0)
        rg3=rhog3(iz,ng_2d)
        if (esm_bc.eq.'bc1') then
           tmp1=tmp1+rg3*ci*(cc+ci*ss)/kn
           tmp2=tmp2+rg3*ci*(cc-ci*ss)/kn
           tmp3=tmp3+rg3*cc/kn**2
           tmp4=tmp4+(0.d0,0.d0)
        else if (esm_bc.eq.'bc2') then
           tmp1=tmp1+rg3*(cc+ci*ss)/kn**2
           tmp2=tmp2+rg3*(cc-ci*ss)/kn**2
           tmp3=tmp3+rg3*ci*cc/kn
           tmp4=tmp4+rg3*ss/kn
        else if (esm_bc.eq.'bc3') then
           tmp1=tmp1+rg3*(cc+ci*ss)/kn**2
           tmp2=tmp2+rg3*(cc-ci*ss)/kn
           tmp3=tmp3+rg3*(cc+ci*ss)/kn
           tmp4=tmp4+(0.d0,0.d0)
        endif
        vg2(iz)=fpi*rg3/(kn**2)
        !for smoothing
        c_r=cos(kn*z_r)
        s_r=sin(kn*z_r)
        c_l=cos(kn*z_l)
        s_l=sin(kn*z_l)
        f1=f1+fpi*   rg3*(c_r+ci*s_r)/kn**2
        f2=f2+fpi*   rg3*(c_l+ci*s_l)/kn**2
        f3=f3+fpi*ci*rg3*(c_r+ci*s_r)/kn
        f4=f4+fpi*ci*rg3*(c_l+ci*s_l)/kn
        !
     enddo
     
     vg2_in(1:nr3x)=vg2(1:nr3x)  ! Since cft_1z is not in-place
     call cft_1z(vg2_in,1,nr3x,nr3x,1,vg2)
     
     rg3=rhog3(1,ng_2d)
     do iz=1,nr3x
        k3=iz-1
        if (k3.gt.nr3x/2) k3=iz-nr3x-1
        z=dble(k3)/dble(nr3x)*L
        if (esm_bc.eq.'bc1') then
           vg2(iz)=vg2(iz)-tpi*z**2*rg3    &
                          -tpi*(z-z0)*tmp1 &
                          -tpi*(z+z0)*tmp2 &
                          -fpi*tmp3        &
                          +vbar(iz)
        else if (esm_bc.eq.'bc2') then
           vg2(iz)=vg2(iz)-tpi*z**2*rg3          &
                          -tpi*(z+z1)*tmp1/z1    &
                          +tpi*(z-z1)*tmp2/z1    &
                          -fpi*z*(z1-z0)/z1*tmp3 &
                          +fpi*(z1-z0)*tmp4      &
                          +vbar(iz)
        else if (esm_bc.eq.'bc3') then
           vg2(iz)=vg2(iz)-tpi*(z**2+2.d0*z*z0)*rg3 &
                          -fpi*tmp1                 &
                          -fpi*ci*(z-z0)*tmp2       &
                          -fpi*ci*(z1-z0)*tmp3      &
                          +vbar(iz)
        endif
     enddo
     !for smoothing
     if (esm_bc.eq.'bc1') then
        f1=f1-tpi*z_r**2*rg3 &
             -tpi*(z_r-z0)*tmp1 &
             -tpi*(z_r+z0)*tmp2 &
             -fpi*tmp3
        f1=f1-tpi*z0**2*rg3+vbar(nz_r)
        f2=f2-tpi*z_l**2*rg3 &
             -tpi*(z_l-z0)*tmp1 &
             -tpi*(z_l+z0)*tmp2 &
             -fpi*tmp3
        f2=f2-tpi*z0**2*rg3+vbar(nz_l)
        f3=f3-tpi*tmp1-tpi*tmp2-fpi*z_r*rg3+vbar(nz_r)
        f4=f4-tpi*tmp1-tpi*tmp2-fpi*z_l*rg3+vbar(nz_l)
     else if (esm_bc.eq.'bc2') then
        f1=f1-tpi*z_r**2*rg3 &
             -tpi*(z_r+z1)*tmp1/z1 &
             +tpi*(z_r-z1)*tmp2/z1 &
             -fpi*z*(z1-z0)/z1*tmp3 &
             +fpi  *(z1-z0)   *tmp4
        f1=f1+tpi*(2.d0*z1-z0)*z0*rg3+vbar(nz_r)
        f2=f2-tpi*z_l**2*rg3 &
             -tpi*(z_l+z1)*tmp1/z1 &
             +tpi*(z_l-z1)*tmp2/z1 &
             -fpi*z*(z1-z0)/z1*tmp3 &
             +fpi  *(z1-z0)   *tmp4
        f2=f2+tpi*(2.d0*z1-z0)*z0*rg3+vbar(nz_l)
        f3=f3-fpi*z_r*rg3-tpi*tmp1/z1+tpi*tmp2/z1-fpi*(z1-z0)/z1*tmp3+vbar(nz_r)
        f4=f4-fpi*z_l*rg3-tpi*tmp1/z1+tpi*tmp2/z1-fpi*(z1-z0)/z1*tmp3+vbar(nz_l)
     else if (esm_bc.eq.'bc3') then
        f1=f1-tpi*(z_r**2+2.d0*z_r*z0)*rg3 &
             -fpi*tmp1 &
             -fpi*ci*(z_r-z1)*tmp2 &
             -fpi*ci*(z1 -z0)*tmp3
        f1=f1+tpi*(4.d0*z1-z0)*z0*rg3+vbar(nz_r)
        f2=f2-tpi*(z_l**2+2.d0*z_l*z0)*rg3 &
             -fpi*tmp1 &
             -fpi*ci*(z_l-z1)*tmp2 &
             -fpi*ci*(z1 -z0)*tmp3
        f2=f2+tpi*(4.d0*z1-z0)*z0*rg3+vbar(nz_l)
        f3=f3-tpi*(2.d0*z_r+2.d0*z0)*rg3-fpi*ci*tmp2+vbar(nz_r)
        f4=f4-tpi*(2.d0*z_l+2.d0*z0)*rg3-fpi*ci*tmp2+vbar(nz_l)
     endif
     ! for smoothing
     !factor 2 will be multiplied later (at vg3 <= vg2)
     !f1=f1*2.d0; f2=f2*2.d0; f3=f3*2.d0; f4=f4*2.d0
     z_r=z_r
     z_l=z_l+L
     a0=(f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r) &
          +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
     a1=(f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r &
          -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
     a2=(-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
          +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
     a3=(2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
     do iz=nz_r,nz_l
        z=dble(iz-1)/dble(nr3x)*L
        vg2(iz)=(a0+a1*z+a2*z**2+a3*z**3)
     enddo
     
     vg2_in(1:nr3x)=vg2(1:nr3x)  ! Since cft_1z is not in-place
     call cft_1z(vg2_in,1,nr3x,nr3x,-1,vg2)
     
     vg3(1:nr3x,ng_2d)=vg2(1:nr3x)*2.d0
     
!#ifdef __MPI__
!  call mpi_allreduce(MPI_IN_PLACE,vg3(1,ng_2d),nr3x,mpi_double_complex,mpi_sum,communicator,ierr)
!#endif
     deallocate(vg2,vg2_in)
  endif ! if( ng_2d > 0 )
  
! Hartree Energy
  ehart=0.d0
!$omp parallel private( ng_2d, k1, k2, k, eh )
  eh = 0d0
!$omp do
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     eh = eh + sum( vg3(:,ng_2d)*conjg(rhog3(:,ng_2d)) )
  enddo
!$omp atomic
  ehart=ehart+eh
!$omp end parallel
  if( gamma_only ) then
     ehart = ehart * 2d0
     ng_2d = imill_2d(0,0)
     if( ng_2d > 0 ) then
        ehart = ehart - sum( vg3(:,ng_2d)*conjg(rhog3(:,ng_2d)) )
     endif
  endif
  ehart = ehart *omega*0.5d0
#ifdef __PARA
  call mp_sum( ehart, intra_pool_comm )
#endif
!#ifdef __MPI__
!  call mpi_allreduce(MPI_IN_PLACE,ehart,1,mpi_double_precision,mpi_sum,communicator,ierr)
!#endif

! Map to FFT mesh (nrxx)
  aux=0.0d0
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng) + 1
     if (n3<1) n3 = n3 + nr3x
     aux(nl(ng))= aux(nl(ng)) + vg3(n3,ng_2d)
  enddo
  if (gamma_only) then
     do ng=1,ngm
        aux(nlm(ng))=CONJG(aux(nl(ng)))
     enddo
  endif
  deallocate (vg3,rhog3)
  
  RETURN
END SUBROUTINE esm_hartree

SUBROUTINE esm_ewald_r(alpha,ewaldr)
  Use ESM_VARS
  implicit none
  real(8),intent(in) :: alpha
  real(8),intent(out) :: ewaldr
  integer :: i
  ewaldr = 0.d0
  Call Ewald( nat, tau, upf_zp, at, alat, alpha, bg, ewaldr )
end SUBROUTINE esm_ewald_r

SUBROUTINE esm_ewald_g(alpha,ewaldg)
  Use ESM_VARS
  implicit none
  real(8),intent(in) :: alpha
  real(8),intent(out) :: ewaldg
  integer :: i
  ewaldg = 0.d0
  Call esm_ewald( upf_zp, alpha, ewaldg )
end SUBROUTINE esm_ewald_g

!-----------------------------------------------------------------------
!--------------ESM EWALD SUBROUTINE-------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE esm_ewald ( charge, alpha, ewg )
!SUBROUTINE esm_ewald (  alpha, ewg )
  Use ESM_VARS
  implicit none
  !
#ifdef __MPI__
  include 'mpif.h'
#endif
  REAL(8), Intent(in)     :: charge(nat), alpha
  REAL(8), Intent(out)    :: ewg
  !
  !    here the local variables
  !
  real(8), external      :: qe_erfc, qe_erf
  real(8)                :: gp2, t(2), gp, sa, z1, z0, L
  integer                 :: k1, k2, k3, ipol, it1, it2, ng_2d, ierr
  real(8) :: tt, z, zp, kk1, kk2, g, cc1, cc2, arg1, arg2, t1, t2, ff, argmax, ew
#ifdef __OPENMP
  INTEGER :: nth, ith, omp_get_thread_num, omp_get_num_threads
#endif
  real(8) :: v0
  
  argmax=0.9*log(huge(1.d0))
  ewg=0.d0
  
  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+abs(esm_w)
  if(gew<0)then
  g=sqrt(alpha)
  else
  g=gew
  endif
  !g = sqrt(2.6)
  sa=omega/L
  
!$omp parallel private( nth, ith, ew, it1, it2, z, zp, tt, kk1, kk2, cc1, cc2, &
!$omp                   ng_2d, k1, k2, gp2, ipol, t, gp, ff, arg1, arg2, t1, t2 )
#ifdef __OPENMP
  nth=omp_get_num_threads()
  ith=omp_get_thread_num()
#endif
  ew=0d0
  do it1=1,nat
  do it2=1,it1
#ifdef __OPENMP
     if( mod( (it1-1)*it1/2+it2-1, nth) /= ith ) cycle
#endif
#ifdef __MPI__
     if( mod( (it1-1)*nat+it2-1, npes) /= mype ) cycle
#endif

     z=tau(3,it1)
     if (z.gt.at(3,3)*0.5) z=z-at(3,3)
     z=z*alat
     
     zp=tau(3,it2)
     if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
     zp=zp*alat
     
     tt=upf_zp(it1)*upf_zp(it2)*2.0*pi/sa
     
     kk1=0.5d0*(-(z-zp)*qe_erf(g*(z-zp))-exp(-g**2*(z-zp)**2)/g/sqrt(pi))
     
     if (esm_bc.eq.'bc1') then
        kk2=0.d0
     else if (esm_bc.eq.'bc2') then
        kk2=0.5d0*(z1-z*zp/z1)
     else if (esm_bc.eq.'bc3') then
        kk2=0.5d0*(2.d0*z1-z-zp)
     endif
     
     cc1=0.d0
     cc2=0.d0
     
     if (it1.eq.it2) then
        
        do ng_2d = 1, ngm_2d
           k1 = mill_2d(1,ng_2d)
           k2 = mill_2d(2,ng_2d)
           if( k1==0 .and. k2==0 ) cycle
           
           t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
           gp2 = sum( t(:) * t(:) ) * tpiba2
           gp=sqrt(gp2)
           
           arg1=-gp*(z-zp)
           arg2= gp*(z-zp)
           arg1=min(arg1,argmax)
           arg2=min(arg2,argmax)
           t1=exp(arg1)*qe_erfc(gp/2.d0/g-g*(z-zp))
           t2=exp(arg2)*qe_erfc(gp/2.d0/g+g*(z-zp))
           cc1=cc1+(t1+t2)/4.d0/gp
           
           if (esm_bc.eq.'bc1') then
              cc2=0.d0
           else if (esm_bc.eq.'bc2') then
              cc2=cc2+(exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                      -exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                      /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp
           else if (esm_bc.eq.'bc3') then
              cc2=cc2+(-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp
           endif
           
        enddo
        
        if( gamma_only ) then
           cc1 = cc1 * 2d0
           cc2 = cc2 * 2d0
        endif
        ew=ew+tt*(cc1+cc2)
        if(gstart==2) ew=ew+tt*(kk1+kk2)
        
     else
        
        do ng_2d = 1, ngm_2d
           k1 = mill_2d(1,ng_2d)
           k2 = mill_2d(2,ng_2d)
           if( k1==0 .and. k2==0 ) cycle
           
           t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
           gp2 = sum( t(:) * t(:) ) * tpiba2
           gp=sqrt(gp2)
           
           ff = ( ( k1*bg(1,1)+k2*bg(1,2) ) * ( tau(1,it1)-tau(1,it2) )  &
              +   ( k1*bg(2,1)+k2*bg(2,2) ) * ( tau(2,it1)-tau(2,it2) ) ) * tpi
           arg1=-gp*(z-zp)
           arg2= gp*(z-zp)
           arg1=min(arg1,argmax)
           arg2=min(arg2,argmax)
           t1=exp(arg1)*qe_erfc(gp/2.d0/g-g*(z-zp))
           t2=exp(arg2)*qe_erfc(gp/2.d0/g+g*(z-zp))
           cc1=cc1+cos(ff)*(t1+t2)/4.d0/gp
           
           if (esm_bc.eq.'bc1') then
              cc2=0.d0
           else if (esm_bc.eq.'bc2') then
              cc2=cc2+cos(ff)*(exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                              -exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                              /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp
           else if (esm_bc.eq.'bc3') then
              cc2=cc2+cos(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp
           endif

        enddo
        
        if( gamma_only ) then
           cc1 = cc1 * 2d0
           cc2 = cc2 * 2d0
        endif
        ew=ew+tt*(cc1+cc2)*2d0
        if(gstart==2) ew=ew+tt*(kk1+kk2)*2d0
        
     endif
     
  enddo
  enddo
#ifdef __MPI__
  call mpi_allreduce(MPI_IN_PLACE,ew,1,mpi_double_precision,mpi_sum,communicator,ierr)
#endif
!$omp atomic
  ewg=ewg+ew
!$omp end parallel
  
  ewg=2.0*ewg
  
  if( gstart == 2 ) then
     do it1=1,nat
        ewg=ewg- upf_zp(it1) **2 * sqrt (8.d0 / tpi * alpha)
     enddo
  endif
  
  if(esm_bc.eq.'bc2')then
     v0=esm_efield*z1*2.d0/2.d0 ! factor 1/2: unit Ry. -> hartree
     do it1=1,nat
        z = tau(3,it1) 
        if (z.gt.at(3,3)*0.5) z=z-at(3,3)
        z=z*alat
        ewg = ewg-upf_zp(it1)*v0*(1.d0-z/z1)
     enddo
     ewg=ewg+v0**2*sa/(pi*8.d0)/z1
  endif

  return
end subroutine esm_ewald

subroutine esm_local(aux)
  Use ESM_VARS
  implicit none
  COMPLEX(8)             :: aux( nrxx )     ! aux contains v_loc_short(G) (input) and v_loc(G) (output)
  real(8),allocatable    :: agauss(:,:),bgauss(:,:) 
  integer :: it
  allocate(agauss(nat,1),bgauss(nat,1))
  do it=1,nat
     agauss(it,1)=1.d0
     bgauss(it,1)=1.d0
  enddo
  call esm_local_(nrxx,aux,nat,1,agauss,bgauss)
  deallocate(agauss,bgauss)
end subroutine esm_local

!-----------------------------------------------------------------------
!--------------ESM LOCAL POTENTIAL SUBROUTINE---------------------------
!-----------------------------------------------------------------------
subroutine esm_local_(nrx,aux,natm,ngauss,agauss,bgauss)
  Use ESM_VARS
  implicit none
#ifdef __MPI__
  include 'mpif.h'
#endif
  !
  integer,intent(in)       :: nrx
  COMPLEX(8),intent(inout) :: aux( nrx )     ! aux contains v_loc_short(G) (input) and v_loc(G) (output)
  integer, intent(in)      :: natm,ngauss
  real(8), intent(in)      :: agauss(natm,ngauss),bgauss(natm,ngauss)
  !
  !    here the local variables
  !
  complex(8),allocatable :: vloc3(:,:),vg2(:),vg2_in(:)
!  real(8), external      :: qe_erf, qe_erfc
  real(8)      :: qe_erf, qe_erfc
  real(8)                :: t(3),tt,gp,gp2,sa,z1,z0,pp,cc,ss,t1,t2, &
                            z,zp,arg11,arg12,arg21,arg22,v0,tmp,L,argmax, &
                            z_l,z_r
  integer                :: iz,ig,it,ipol,k1,k2,k3,ng,n1,n2,n3, &
                            nz_l,nz_r, ng_2d
  complex(8)             :: cs,cc1,cc2,ci,a0,a1,a2,a3,f1,f2,f3,f4
#ifdef __MPI__
  integer                :: ind,ierr
#endif

  argmax=0.9*log(huge(1.d0))
  L =at(3,3)*alat
  z0=L/2.d0
  z1=z0+abs(esm_w)
  allocate(vloc3(nr3x,ngm_2d))
  vloc3(:,:)=(0.d0,0.d0)
  sa=omega/L
  v0=esm_efield*z1*2.d0/2.d0 ! factor 1/2: unit Ry. -> hartree
  ci=(0.d0,1.d0)

! for gp!=0
!$omp parallel private( k1, k2, gp2, gp, vg2, it, tt, pp, cc, ss, cs, zp, iz, &
!$omp                   k3, z, cc1, ig, tmp, arg11, arg12, arg21, arg22, t1, t2, &
!$omp                   cc2, vg2_in )
  allocate(vg2(nr3x),vg2_in(nr3x))
  vg2(:) = (0.d0,0.d0);vg2_in(:) = (0.d0,0.d0)
!$omp do
  do ng_2d = 1, ngm_2d
#ifdef __MPI__
     if (mod(ng_2d,npes) /= mype) cycle
#endif
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle
     
     t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
     gp2 = sum( t(:) * t(:) ) * tpiba2
     gp=sqrt(gp2)
     
     vg2(1:nr3x)=(0.d0,0.d0)
     do it=1,nat
        tt=-fpi*upf_zp(it)/sa
        pp=-tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2))+tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
        cc=cos(pp)
        ss=sin(pp)
        cs=CMPLX ( cc, ss, kind=8 )
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,nr3x
           k3=iz-1
           if (k3.gt.nr3x/2) k3=iz-nr3x-1
           z=dble(k3)/dble(nr3x)*L
           cc1=(0.d0,0.d0)
           do ig=1,ngauss
              tmp=sqrt(agauss(it,ig))
              arg11=-gp*(z-zp)
              arg11=min(arg11,argmax)
              arg12= gp/2.d0/tmp-tmp*(z-zp)
              arg21= gp*(z-zp)
              arg21=min(arg21,argmax)
              arg22= gp/2.d0/tmp+tmp*(z-zp)
              t1=exp(arg11)*qe_erfc(arg12)
              t2=exp(arg21)*qe_erfc(arg22)
              cc1=cc1+bgauss(it,ig)*cs*(t1+t2)/4.d0/gp
           enddo
           if (esm_bc.eq.'bc1') then
              cc2=(0.d0,0.d0)
           else if (esm_bc.eq.'bc2') then
              cc2=cs*( exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                      -exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1))) &
                      /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp 
           else if (esm_bc.eq.'bc3') then
              cc2=cs*(-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp
           endif
           vg2(iz) = vg2(iz) + tt*(cc1+cc2)*2.d0 ! factor 2: hartree -> Ry.
        enddo
     enddo
     vg2_in(1:nr3x)=vg2(1:nr3x)
     call cft_1z(vg2_in,1,nr3x,nr3x,-1,vg2)
     do iz=1,nr3x
        vloc3(iz,ng_2d)=vg2(iz)
     enddo
  enddo
#ifdef __MPI__
  call mpi_allreduce(MPI_IN_PLACE,vloc3,nr3x * ngm_2d,mpi_double_complex,mpi_sum,communicator,ierr)
#endif
  deallocate(vg2,vg2_in)
!$omp end parallel

  ng_2d=imill_2d(0,0)
  if( ng_2d > 0 ) then
     allocate(vg2(nr3x),vg2_in(nr3x))

     vg2(1:nr3x)=(0.d0,0.d0)
! for smoothing
     f1=0.d0; f2=0.d0; f3=0.d0; f4=0.d0
     nz_l=nr3x/2+1+esm_nfit
     nz_r=nr3x/2+1-esm_nfit
     z_l=dble(nz_l-1)*L/dble(nr3x)-L
     z_r=dble(nz_r-1)*L/dble(nr3x)
! add constant potential (capacitor term)
     do iz=1,nr3x
        k3=iz-1
        if (k3.gt.nr3x/2) k3=iz-nr3x-1
        z=dble(k3)/dble(nr3x)*L
        vg2(iz)=-0.5d0*v0*(z-z1)/z1*2.d0 ! factor 2: hartree -> Ry.
     enddo
     f1=-0.5d0*v0*(z_r-z1)/z1 ! unit: hartree
     f2=-0.5d0*v0*(z_l-z1)/z1 ! unit: hartree
     f3=-0.5d0*v0/z1 ! unit: hartree/a.u.
     f4=-0.5d0*v0/z1 ! unit: harteee/a.u.
! for gp=0
     do it=1,nat
        tt=-fpi*upf_zp(it)/sa
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,nr3x
           k3=iz-1
           if (k3.gt.nr3x/2) k3=iz-nr3x-1
           z=dble(k3)/dble(nr3x)*L
           cc1=(0.d0,0.d0) 
           do ig=1,ngauss
              tmp=sqrt(agauss(it,ig))
              cc1=cc1+bgauss(it,ig)*0.5d0*(-(z-zp)*qe_erf(tmp*(z-zp)) &
                   -exp(-tmp**2*(z-zp)**2)/tmp/sqrt(pi))
           enddo
           if (esm_bc.eq.'bc1') then
              cc2=(0.d0,0.d0)
           else if (esm_bc.eq.'bc2') then
              cc2=0.5d0*(z1-z*zp/z1)
           else if (esm_bc.eq.'bc3') then
              cc2=0.5d0*(2.d0*z1-z-zp)
           endif
           vg2(iz) = vg2(iz) + tt*(cc1+cc2)*2.d0 ! factor 2: hartree -> Ry.
        enddo
     ! smoothing cell edge potential (avoiding unphysical oscillation)
        do ig=1,ngauss
           tmp=sqrt(agauss(it,ig))
           f1=f1+tt*bgauss(it,ig)*0.5d0*(-(z_r-zp)*qe_erf(tmp*(z_r-zp)) &
                -exp(-tmp**2*(z_r-zp)**2)/tmp/sqrt(pi))
           f2=f2+tt*bgauss(it,ig)*0.5d0*(-(z_l-zp)*qe_erf(tmp*(z_l-zp)) &
                -exp(-tmp**2*(z_l-zp)**2)/tmp/sqrt(pi))
           f3=f3-tt*bgauss(it,ig)*0.5d0*qe_erf(tmp*(z_r-zp))
           f4=f4-tt*bgauss(it,ig)*0.5d0*qe_erf(tmp*(z_l-zp))
        enddo
        if(esm_bc.eq.'bc1')then
           f1=f1+tt*0.d0
           f2=f2+tt*0.d0
           f3=f3+tt*0.d0
           f4=f4+tt*0.d0
        elseif(esm_bc.eq.'bc2')then
           f1=f1+tt*0.5d0*(z1-z_r*zp/z1)
           f2=f2+tt*0.5d0*(z1-z_l*zp/z1)
           f3=f3+tt*(-0.5d0*(zp/z1))
           f4=f4+tt*(-0.5d0*(zp/z1))
        elseif(esm_bc.eq.'bc3')then
           f1=f1+tt*0.5d0*(2.d0*z1-z_r-zp)
           f2=f2+tt*0.5d0*(2.d0*z1-z_l-zp)
           f3=f3-tt*0.5d0
           f4=f4-tt*0.5d0
        endif
     enddo
     ! for smoothing
     f1=f1*2.d0; f2=f2*2.d0; f3=f3*2.d0; f4=f4*2.d0 ! factor 2: hartree -> Ry.
     z_r=z_r
     z_l=z_l+L
     a0=(f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r) &
          +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
     a1=(f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r &
          -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
     a2=(-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
          +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
     a3=(2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
     do iz=nz_r,nz_l
        z=dble(iz-1)/dble(nr3x)*L
        vg2(iz)=(a0+a1*z+a2*z**2+a3*z**3)
     enddo
     vg2_in(1:nr3x)=vg2(1:nr3x)
     call cft_1z(vg2_in,1,nr3x,nr3x,-1,vg2)
     do iz=1,nr3x
        vloc3(iz,ng_2d)=vg2(iz)
     enddo
     
     deallocate(vg2,vg2_in)
  endif ! if( ng_2d > 0 )
  
! Map to FFT mesh (nrxx)
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng) + 1 
     IF (n3<1) n3 = n3 + nr3x
     aux(nl(ng))= aux(nl(ng)) + vloc3(n3,ng_2d)
  enddo
  if (gamma_only) then
     do ng=1,ngm
        aux (nlm(ng))=CONJG(aux(nl(ng)))
     enddo
  endif

  deallocate(vloc3)

  return
  end subroutine esm_local_

!-----------------------------------------------------------------------
!--------------ESM EWALD-DERIVED FORCE SUBROUTINE-----------------------
!-----------------------------------------------------------------------
subroutine esm_force_ew ( alpha, forceion ) 
  Use ESM_VARS
  implicit none
#ifdef __MPI__
  include 'mpif.h'
#endif
  !
  REAL(8)                :: alpha
  REAL(8)                :: forceion(3,nat) 
  !
  !    here the local variables
  !
  real(8), external      :: qe_erfc, qe_erf
  integer  :: it1, it2, ipol, k1, k2, k3, ng_2d
  integer  :: nth, ith, omp_get_num_threads, omp_get_thread_num,ierr
  real(8) :: t1_for, t2_for, z, zp, kk1_for, kk2_for, g, gp2, gp, z1, t(2), L
  real(8) :: cx1_for, cy1_for, cz1_for, cx2_for, cy2_for, cz2_for, arg1, arg2, t1, t2, ff
  real(8) :: sa, z0, g_b,tauz1,tauz2,gt,tt,gz,argmax
  real(8), allocatable :: for(:,:), for_g(:,:)
  real(8) :: v0

  argmax=0.9*log(huge(1.d0))
  allocate(for_g(3,nat))
  for_g(:,:)=0.d0
  
  forceion(:,:)=0.d0
  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+abs(esm_w) 
  sa=omega/L
  g=sqrt(alpha)

!$omp parallel private( nth, ith, for, z, zp, t1_for, t2_for, kk1_for, kk2_for, &
!$omp                   cz1_for, cz2_for, ng_2d, k1, k2, gp2, gp, arg1, arg2, t1, t2, &
!$omp                   cx1_for, cy1_for, cx2_for, cy2_for, ff, t )
#ifdef __OPENMP
  nth=omp_get_num_threads()
  ith=omp_get_thread_num()
#endif
  allocate(for(3,nat))
  for=0d0
  do it1=1,nat
  do it2=1,nat
#ifdef __OPENMP
     if( mod( (it1-1)*nat+it2-1, nth) /= ith ) cycle
#endif
#ifdef __MPI__
     if( mod( (it1-1)*nat+it2-1, npes) /= mype ) cycle
#endif
     z=tau(3,it1)
     if (z.gt.at(3,3)*0.5) z=z-at(3,3)
     z=z*alat
     zp=tau(3,it2)
     if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
     zp=zp*alat
     if (gamma_only) then 
        t1_for=upf_zp(it1)*upf_zp(it2)*fpi/sa*2.d0
     else
        t1_for=upf_zp(it1)*upf_zp(it2)*fpi/sa
     endif
     t2_for=upf_zp(it1)*upf_zp(it2)*fpi/sa

     kk1_for=0.5d0*qe_erf(g*(z-zp))
     if (esm_bc.eq.'bc1') then
        kk2_for=0.d0
     else if (esm_bc.eq.'bc2') then
        kk2_for=-0.5d0*(z/z1)
     else if (esm_bc.eq.'bc3') then
        kk2_for=-0.5d0
     endif
  
     if (it1.eq.it2) then
        cz1_for=0.d0
        cz2_for=0.d0
        
        do ng_2d = 1, ngm_2d
           k1 = mill_2d(1,ng_2d)
           k2 = mill_2d(2,ng_2d)
           if( k1==0 .and. k2==0 ) cycle

           t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
           gp2 = sum( t(:) * t(:) ) * tpiba2
           gp=sqrt(gp2)
           
           arg1=-gp*(z-zp)
           arg2= gp*(z-zp)
           arg1=min(arg1,argmax)
           arg2=min(arg2,argmax)
           t1=exp(arg1)*qe_erfc(gp/2.d0/g-g*(z-zp))
           t2=exp(arg2)*qe_erfc(gp/2.d0/g+g*(z-zp))
           cz1_for=0.d0
           if (esm_bc.eq.'bc1') then      
              cz2_for=0.d0
           else if (esm_bc.eq.'bc2') then      
              cz2_for=cz2_for - (exp(gp*(z-zp-4.d0*z1))-exp(-gp*(z-zp+4.d0*z1)) &
                                +exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                                /(1.d0-exp(-4.d0*gp*z1))/2.d0
           else if (esm_bc.eq.'bc3') then      
              cz2_for=cz2_for - exp(gp*(z+zp-2.d0*z1))/2.d0
           endif
        enddo
        for(3,it2) = for(3,it2) + t1_for*(cz1_for+cz2_for)
        if(gstart==2) then
           for(3,it2) = for(3,it2) + t2_for*(kk1_for+kk2_for)
        endif
        
     else if (it1.gt.it2) then
        
        cx1_for=0.d0
        cy1_for=0.d0
        cz1_for=0.d0
        cx2_for=0.d0
        cy2_for=0.d0
        cz2_for=0.d0
        do ng_2d = 1, ngm_2d
           k1 = mill_2d(1,ng_2d)
           k2 = mill_2d(2,ng_2d)
           if(k1==0.and.k2==0) cycle
           
           t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
           gp2 = sum( t(:) * t(:) ) * tpiba2
           gp=sqrt(gp2)
           
           ff = ( ( k1*bg(1,1)+k2*bg(1,2) ) * ( tau(1,it1)-tau(1,it2) )  &
              +   ( k1*bg(2,1)+k2*bg(2,2) ) * ( tau(2,it1)-tau(2,it2) ) ) * tpi
           arg1=-gp*(z-zp)
           arg2= gp*(z-zp)
           arg1=min(arg1,argmax)
           arg2=min(arg2,argmax)
           t1=exp(arg1)*qe_erfc(gp/2.d0/g-g*(z-zp))
           t2=exp(arg2)*qe_erfc(gp/2.d0/g+g*(z-zp))
           
           cx1_for=cx1_for+sin(ff)*(t1+t2)/4.d0/gp*k1
           cy1_for=cy1_for+sin(ff)*(t1+t2)/4.d0/gp*k2
           cz1_for=cz1_for+cos(ff)*(t1-t2)/4.d0
           if (esm_bc.eq.'bc1') then
              cx2_for=0.d0
              cy2_for=0.d0
              cz2_for=0.d0
           else if (esm_bc.eq.'bc2') then
              cx2_for=cx2_for + sin(ff)*(exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                                       - exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                                       /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp*k1
              cy2_for=cy2_for + sin(ff)*(exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                                       - exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                                       /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp*k2
              cz2_for=cz2_for - cos(ff)*(exp(gp*(z-zp-4.d0*z1))-exp(-gp*(z-zp+4.d0*z1)) &
                                       + exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                                       /(1.d0-exp(-4.d0*gp*z1))/2.d0
           else if (esm_bc.eq.'bc3') then
              cx2_for=cx2_for+sin(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp*k1
              cy2_for=cy2_for+sin(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp*k2
              cz2_for=cz2_for+cos(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0
           endif
        enddo
        for(1,it2)=for(1,it2)+t1_for*(cx1_for+cx2_for)
        for(2,it2)=for(2,it2)+t1_for*(cy1_for+cy2_for)
        for(3,it2)=for(3,it2)+t1_for*(cz1_for+cz2_for)
        if(gstart==2) then
           for(3,it2)=for(3,it2)+t2_for*(kk1_for+kk2_for)
        endif
        
     else if (it1.lt.it2) then
        
        cx1_for=0.d0
        cy1_for=0.d0
        cz1_for=0.d0
        cx2_for=0.d0
        cy2_for=0.d0
        cz2_for=0.d0
        do ng_2d = 1, ngm_2d
           k1 = mill_2d(1,ng_2d)
           k2 = mill_2d(2,ng_2d)
           if(k1==0.and.k2==0) cycle
           
           t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
           gp2 = sum( t(:) * t(:) ) * tpiba2
           gp=sqrt(gp2)
           
           ff = ( ( k1*bg(1,1)+k2*bg(1,2) ) * ( tau(1,it1)-tau(1,it2) )  &
              +   ( k1*bg(2,1)+k2*bg(2,2) ) * ( tau(2,it1)-tau(2,it2) ) ) * tpi
           arg1=-gp*(z-zp)
           arg2= gp*(z-zp)
           arg1=min(arg1,argmax)
           arg2=min(arg2,argmax)
           t1=exp(arg1)*qe_erfc(gp/2.d0/g-g*(z-zp))
           t2=exp(arg2)*qe_erfc(gp/2.d0/g+g*(z-zp))
           
           cx1_for=cx1_for+sin(ff)*(t1+t2)/4.d0/gp*k1
           cy1_for=cy1_for+sin(ff)*(t1+t2)/4.d0/gp*k2
           cz1_for=cz1_for+cos(ff)*(t1-t2)/4.d0
           if (esm_bc.eq.'bc1') then
              cx2_for=0.d0
              cy2_for=0.d0
              cz2_for=0.d0
           else if (esm_bc.eq.'bc2') then
              cx2_for=cx2_for + sin(ff)*(exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                                       - exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                                       /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp*k1
              cy2_for=cy2_for + sin(ff)*(exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                                       - exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                                       /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp*k2
              cz2_for=cz2_for - cos(ff)*(exp(gp*(z-zp-4.d0*z1))-exp(-gp*(z-zp+4.d0*z1)) &
                                       + exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1)) ) &
                                       /(1.d0-exp(-4.d0*gp*z1))/2.d0
           else if (esm_bc.eq.'bc3') then
              cx2_for=cx2_for+sin(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp*k1
              cy2_for=cy2_for+sin(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp*k2
              cz2_for=cz2_for+cos(ff)*(-exp(gp*(z+zp-2.d0*z1)))/2.d0
           endif
        enddo
        for(1,it2)=for(1,it2)+t1_for*(cx1_for+cx2_for)
        for(2,it2)=for(2,it2)+t1_for*(cy1_for+cy2_for)
        for(3,it2)=for(3,it2)+t1_for*(cz1_for+cz2_for)
        if(gstart==2) then
           for(3,it2)=for(3,it2)+t2_for*(kk1_for+kk2_for)
        endif
     endif
     
  enddo
  enddo
#ifdef __MPI__
  call mpi_allreduce(MPI_IN_PLACE,for,nat*3,mpi_double_precision,mpi_sum,communicator,ierr)
#endif
!$omp critical
  for_g(:,:) = for_g(:,:) + for(:,:)
!$omp end critical
  deallocate(for)
!$omp end parallel

  for_g(:,:)=2.0*for_g(:,:)

  do it1=1,nat
     forceion(1,it1)=sum( for_g(1:2,it1)*bg(1,1:2) )*sqrt(tpiba2)
     forceion(2,it1)=sum( for_g(1:2,it1)*bg(2,1:2) )*sqrt(tpiba2)
     forceion(3,it1)=for_g(3,it1)
  enddo

  if(esm_bc.eq.'bc2')then
     v0=esm_efield*z1*2.d0/2.d0 ! factor 1/2: unit Ry. -> hartree
     do it1=1,nat
        z=tau(3,it1)
        if (z.gt.at(3,3)*0.5) z=z-at(3,3)
        z=z*alat
        forceion(3,it1) = forceion(3,it1)+upf_zp(it1)*v0/z1
     enddo
  endif
  
  deallocate(for_g)
  forceion(:,:)=-forceion(:,:)

  return
end subroutine esm_force_ew

subroutine esm_force_lc ( aux, forcelc )
  Use ESM_VARS
  implicit none
  !
  COMPLEX(8)             :: aux(nrxx)       ! aux contains n(G) (input)   
  REAL(8)                :: forcelc(3,nat)
  real(8),allocatable    :: agauss(:,:),bgauss(:,:) 
  integer :: it
  allocate(agauss(nat,1),bgauss(nat,1))
  do it=1,nat
     agauss(it,1)=1.d0
     bgauss(it,1)=1.d0
  enddo
  call esm_force_lc_(nrxx,aux,nat,forcelc,1,agauss,bgauss)
  deallocate(agauss,bgauss)
end subroutine esm_force_lc

!-----------------------------------------------------------------------
!--------------ESM LOCAL POTENTIAL-DERIVED FORCE SUBROUTINE-------------
!-----------------------------------------------------------------------
subroutine esm_force_lc_ (nrx, aux, natm, forcelc, ngauss, agauss, bgauss )
  Use ESM_VARS
  implicit none
#ifdef __MPI__
  include 'mpif.h'
#endif
  !
  integer,intent(in)     :: nrx
  COMPLEX(8)             :: aux(nrx)       ! aux contains n(G) (input)   
  integer,intent(in)     :: natm
  REAL(8)                :: forcelc(3,natm)
  integer,intent(in)     :: ngauss
  real(8),intent(in)     :: agauss(natm,ngauss),bgauss(natm,ngauss)
  !
  !    here are the local variables
  !
  complex(8),allocatable :: vlocx(:), vlocy(:), vlocdz(:)
  real(8),allocatable    :: for(:,:),for_g(:,:)
  real(8), external      :: qe_erf, qe_erfc
  real(8)                :: t(3),tt,gp,gp2,sa,z1,z0,pp,cc,ss,t1,t2,z,zp,L,forcelc2(3,nat)
  real(8)                :: arg11,arg12,arg21,arg22,tmp,r1,r2,fx1,fy1,fz1,fx2,fy2,fz2,argmax
  integer                 :: iz,ig,it,ipol,k1,k2,k3,ng,n1,n2,n3,ng_2d
  complex(8),allocatable :: vg2(:),vg2_fx(:),vg2_fy(:),vg2_fz(:),rhog3(:,:),vg2_fxyz(:,:)
  complex(8)             :: cx1,cy1,cz1,cx2,cy2,cz2,cc1,cc2
#ifdef __MPI__
  integer                :: ind,ierr
#endif

  argmax=0.9*log(huge(1.d0))

! Mat to FULL FFT mesh (nr1x,nr2x,nr3x)
  allocate(rhog3(nr3x,ngm_2d))
  rhog3(:,:)=(0.d0,0.d0)
  do ng=1,ngm
      n1 = mill(1,ng)
      n2 = mill(2,ng)
      ng_2d = imill_2d(n1,n2)
      n3 = mill(3,ng) + 1
      IF (n3<1) n3 = n3 + nr3x
      rhog3(n3,ng_2d)=aux(nl(ng))
      if( gamma_only .and. n1==0 .and. n2==0 ) then
         n3 = -mill(3,ng) + 1
         IF (n3<1) n3 = n3 + nr3x
         rhog3(n3,ng_2d)=conjg(aux(nl(ng)))
      endif
  enddo

  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+abs(esm_w)

  allocate(for_g(3,nat))

  sa=omega/L
  for_g(:,:)=0.d0

!**** for gp!=0 *********
!$omp parallel private( k1, k2, gp2, gp, it, tt, pp, cc, ss, zp, iz, &
!$omp                   k3, z, cx1, cy1, cz1, tmp, arg11, arg12, arg21, arg22, &
!$omp                   t1, t2, cx2, cy2, cz2, vg2_fx, vg2_fy, vg2_fz, vg2, &
!$omp                   r1, r2, fx1, fy1, fz1, fx2, fy2, fz2, for, t )
  !allocate(for(3,nat),vg2(nr3x),vg2_fx(nr3x),vg2_fy(nr3x),vg2_fz(nr3x))
  allocate(for(3,nat),vg2(nr3x),vg2_fxyz(nr3x,3))
  for(:,:)=0.d0
!  vg2_fx(:)=(0.d0,0.d0)
!  vg2_fy(:)=(0.d0,0.d0)
!  vg2_fz(:)=(0.d0,0.d0)
!$omp do
  do ng_2d = 1, ngm_2d
#ifdef __MPI__
     if (mod(ng_2d,npes) /= mype) cycle
#endif
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle

     t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
     gp2 = sum( t(:) * t(:) ) * tpiba2
     gp=sqrt(gp2)
     
     do it=1,nat
        IF (gamma_only) THEN
           tt=-fpi*upf_zp(it)/sa*2.d0
        ELSE 
           tt=-fpi*upf_zp(it)/sa
        ENDIF 
        pp=-tpi*(tau(1,it)*(k1*bg(1,1)+k2*bg(1,2))+tau(2,it)*(k1*bg(2,1)+k2*bg(2,2)))
        cc=cos(pp)
        ss=sin(pp)
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
!        vg2_fx(:)=(0.d0,0.d0)
!        vg2_fy(:)=(0.d0,0.d0)
!        vg2_fz(:)=(0.d0,0.d0)
        vg2_fxyz(:,:) = (0.d0,0.d0)
        do iz=1,nr3x
           k3=iz-1
           if (k3.gt.nr3x/2) k3=iz-nr3x-1
           z=dble(k3)/dble(nr3x)*L
           cx1=(0.d0,0.d0); cy1=(0.d0,0.d0); cz1=(0.d0,0.d0)
           do ig=1,ngauss
              tmp=sqrt(agauss(it,ig))
              arg11=-gp*(z-zp)
              arg11=min(arg11,argmax) 
              arg12= gp/2.d0/tmp-tmp*(z-zp)
              arg21= gp*(z-zp)
              arg21=min(arg21,argmax)
              arg22= gp/2.d0/tmp+tmp*(z-zp)
              t1=exp(arg11)*qe_erfc(arg12)
              t2=exp(arg21)*qe_erfc(arg22)
              cx1=cx1+bgauss(it,ig)*CMPLX(ss, -cc, kind=8) &
                   *(t1+t2)/4.d0/gp*k1
              cy1=cy1+bgauss(it,ig)*CMPLX(ss, -cc, kind=8) &
                   *(t1+t2)/4.d0/gp*k2
              cz1=cz1+bgauss(it,ig)*CMPLX(cc, ss, kind=8)  &
                   *(t1-t2)/4.d0
           enddo
           if (esm_bc.eq.'bc1') then
              cx2=(0.d0,0.d0)
              cy2=(0.d0,0.d0)
              cz2=(0.d0,0.d0)
           else if (esm_bc.eq.'bc2') then
              cx2=CMPLX(ss, -cc, kind=8)* &
                   (exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                   -exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1))) &
                   /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp*k1
              cy2=CMPLX(ss, -cc, kind=8)* &
                   (exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                   -exp(gp*(z+zp-2.d0*z1))-exp(-gp*(z+zp+2.d0*z1))) &
                   /(1.d0-exp(-4.d0*gp*z1))/2.d0/gp*k2
              cz2=CMPLX(cc, ss, kind=8)* &
                   (-exp(gp*(z-zp-4.d0*z1))+exp(-gp*(z-zp+4.d0*z1)) &
                   -exp(gp*(z+zp-2.d0*z1))+exp(-gp*(z+zp+2.d0*z1))) &
                   /(1.d0-exp(-4.d0*gp*z1))/2.d0
           else if (esm_bc.eq.'bc3') then
              cx2=CMPLX(ss, -cc, kind=8)* &
                   (-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp*k1
              cy2=CMPLX(ss, -cc, kind=8)* &
                   (-exp(gp*(z+zp-2.d0*z1)))/2.d0/gp*k2
              cz2=CMPLX(cc, ss, kind=8)* &
                   (-exp(gp*(z+zp-2.d0*z1)))/2.d0
           endif
!           vg2_fx(iz) = tt*(cx1+cx2)
!           vg2_fy(iz) = tt*(cy1+cy2)
!           vg2_fz(iz) = tt*(cz1+cz2)
           vg2_fxyz(iz,1) = tt*(cx1+cx2)
           vg2_fxyz(iz,2) = tt*(cy1+cy2)
           vg2_fxyz(iz,3) = tt*(cz1+cz2)
        enddo
        vg2(1:nr3x)=vg2_fxyz(1:nr3x,1)  ! Since cft_1z is not in-place
        call cft_1z(vg2,1,nr3x,nr3x,-1,vg2_fxyz(:,1))
        vg2(1:nr3x)=vg2_fxyz(1:nr3x,2)  ! Since cft_1z is not in-place
        call cft_1z(vg2,1,nr3x,nr3x,-1,vg2_fxyz(:,2))
        vg2(1:nr3x)=vg2_fxyz(1:nr3x,3)  ! Since cft_1z is not in-place
        call cft_1z(vg2,1,nr3x,nr3x,-1,vg2_fxyz(:,3))
        do iz=1,nr3x
           r1= dble(rhog3(iz,ng_2d))
           r2=aimag(rhog3(iz,ng_2d))
           fx1=dble(  vg2_fxyz(iz,1))
           fy1=dble(  vg2_fxyz(iz,2))
           fz1=dble(  vg2_fxyz(iz,3))
           fx2=aimag( vg2_fxyz(iz,1))
           fy2=aimag( vg2_fxyz(iz,2))
           fz2=aimag( vg2_fxyz(iz,3))
           for(1,it)=for(1,it)-r1*fx1-r2*fx2
           for(2,it)=for(2,it)-r1*fy1-r2*fy2
           for(3,it)=for(3,it)-r1*fz1-r2*fz2
        enddo
     enddo
  enddo
#ifdef __MPI__
  call mpi_allreduce(MPI_IN_PLACE,for,nat*3,mpi_double_precision,mpi_sum,communicator,ierr)
#endif
!$omp critical
  for_g(:,:) = for_g(:,:) + for(:,:)
!$omp end critical
  !deallocate(for,vg2,vg2_fx,vg2_fy,vg2_fz)
  deallocate(for,vg2,vg2_fxyz)
!$omp end parallel

!***** for gp==0********
  ng_2d = imill_2d(0,0)
  if( ng_2d > 0 ) then
     allocate(vg2(nr3x),vg2_fz(nr3x))

     vg2_fz(:)=(0.d0,0.d0)
     do it=1,nat
        tt=-fpi*upf_zp(it)/sa
        zp=tau(3,it)
        if (zp.gt.at(3,3)*0.5) zp=zp-at(3,3)
        zp=zp*alat
        do iz=1,nr3x
           k3=iz-1
           if (k3.gt.nr3x/2) k3=iz-nr3x-1
           z=dble(k3)/dble(nr3x)*L
           cc1=(0.d0,0.d0)
           do ig=1,ngauss
              tmp=sqrt(agauss(it,ig))
              cc1=cc1+bgauss(it,ig)*(0.5d0*qe_erf(tmp*(z-zp)))
           enddo
           if (esm_bc.eq.'bc1') then
              cc2=(0.d0,0.d0)
           else if (esm_bc.eq.'bc2') then
              cc2=-0.5d0*(z/z1)
           else if (esm_bc.eq.'bc3') then
              cc2=-0.5d0
           endif
           vg2_fz(iz) =  tt*(cc1+cc2)
        enddo
        vg2(1:nr3x)=vg2_fz(1:nr3x)  ! Since cft_1z is not in-place
        call cft_1z(vg2,1,nr3x,nr3x,-1,vg2_fz)
        do iz=1,nr3x
           r1=dble( rhog3(iz,ng_2d))
           r2=aimag(rhog3(iz,ng_2d))
           fz1=dble( vg2_fz(iz))
           fz2=aimag(vg2_fz(iz))
           for_g(3,it)=for_g(3,it)-r1*fz1-r2*fz2
        enddo
     enddo
     
     deallocate(vg2,vg2_fz)
  endif ! if( ng_2d > 0 )

!***** sum short_range part and long_range part in local potential force at cartecian coordinate

  do it=1,nat
     forcelc(1,it)=forcelc(1,it)+sum(for_g(1:2,it)*bg(1,1:2))*sqrt(tpiba2)*omega*2.d0
     forcelc(2,it)=forcelc(2,it)+sum(for_g(1:2,it)*bg(2,1:2))*sqrt(tpiba2)*omega*2.d0
     forcelc(3,it)=forcelc(3,it)+for_g(3,it)*omega*2.d0
  enddo

  deallocate(for_g,rhog3)
  return
end subroutine esm_force_lc_

!-----------------------------------------------------------------------
!--------------ESM FINAL PRINTOUT SUBROUTINE----------------------------
!-----------------------------------------------------------------------
!
! Prints out vlocal and vhartree to stdout once electrons are converged
! Format: z, rho(r), v_hartree, v_local, (v_hartree + v_local) 
!
SUBROUTINE esm_printpot ( rhog, v_h, v_l )
  Use ESM_VARS
  !
  IMPLICIT NONE
  !
  COMPLEX(8), Intent(In)  :: rhog(ngm,nspin)   !  n(G)
  COMPLEX(8), Intent(In)  :: v_h(nrxx)
  COMPLEX(8), Intent(In)  :: v_l(nrxx)
  REAL(8)                 :: charge,ehart,bohr,rydv,L,area
  COMPLEX(8), ALLOCATABLE :: work1(:,:), work2(:,:)
  Real(8)   , ALLOCATABLE :: work3(:,:)
  INTEGER                 :: iz,i,k3,ng,n1,n2,n3,n3_

  L = alat*at(3,3)
  area = (at(1,1)*at(2,2)-at(2,1)*at(1,2))*alat**2
  bohr = 0.52917720859d0
  rydv = 13.6058d0

  ! x,y average
  allocate( work1(nr3x,3), work2(nr3x,3), work3(nr3x,1:5) )
  work1(:,:) = 0d0
  work2(:,:) = 0d0
  work3(:,:) = 0d0

  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     if( n1 == 0 .and. n2 == 0 ) then
        n3 = mill(3,ng) + 1
        if (n3<1) n3 = n3 + nr3x
        if( nspin == 2 ) then
           work1(n3,1) = rhog(ng,1) + rhog(ng,2)
        else
           work1(n3,1) = rhog(ng,1)
        endif
        if( gamma_only ) then
           n3_ = -mill(3,ng) + 1
           if (n3_<1) n3_ = n3_ + nr3x
           work1(n3_,1) = conjg( work1(n3,1) )
        endif
     endif
  enddo

  do iz = 1, nr3x
     k3 = iz - 1
     if( k3 > nr3x/2 ) k3 = k3 - nr3x
     i=1+(iz-1)*nr1x*nr2x
     work1(iz,2) = v_h(i);
     work1(iz,3) = v_l(i);
  enddo

  call cft_1z(work1(:,:),3,nr3x,nr3x,1,work2(:,:))
  work3(:,2) = DBLE( work2(:,1) ) * area
  work3(:,3) = DBLE( work2(:,2) ) * rydv / bohr
  work3(:,4) = DBLE( work2(:,3) ) * rydv / bohr

  do iz = 1, nr3x
     k3 = iz - 1
     if( k3 > nr3x/2 ) k3 = k3 - nr3x
     work3(iz,1) = dble(k3) / nr3x * L * bohr
  enddo

  work3(:,5) = work3(:,3) + work3(:,4)

! z = position along slab (A)
! rho = planar-summed charge density of slab section (e)
! v_hartree = planar-averaged hartree potential term (eV/A)
! v_local = planar-averaged local potential term (eV/A)

  write(*,*)
  write(*,'("     ESM Charge and Potential")')
  write(*,'("     ========================")' )
  write(*,*)
  write(*,'("    z (A)      rho (e)      Avg v_hartree        Avg v_local  Avg v_hart+v_loc")' )
  write(*,'("                                   (eV/A)             (eV/A)            (eV/A)")' )
  write(*,*)
  write(*,'("    ==========================================================================")' )
  do k3 = nr3x/2-nr3x+1, nr3x/2
     iz = k3 + nr3x + 1
     if( iz > nr3x ) iz = iz - nr3x
     write(*,'(f9.3,f13.5,2f19.7,f18.7)') work3(iz,1:5)
  enddo
  write(*,*)

  deallocate(work1, work2, work3)

END SUBROUTINE esm_printpot

SUBROUTINE esm_gen_vbar()
  Use ESM_VARS
  Implicit none
  integer :: i3
  real(8) :: z
  if(izwall.eq.0) bar_height=0.d0
  if(izwall.eq.1)then
     do i3=1,nr3x
        z = dble(i3-1)*at(3,3)*alat/dble(nr3x)
        if(z.gt.at(3,3)*alat/2.0d0)z=z-at(3,3)*alat
        if(z.ge.z_wall.and.z.le.z_wall+bar_width)then
           vbar(i3) = cmplx(0.5d0*(1.d0+sin(-pi/2.d0+pi/bar_width*(z-z_wall)))*bar_height,       &
           &                0.5d0*pi/bar_width*cos(-pi/2.d0+pi/bar_width*(z-z_wall))*bar_height)
        else if (z.gt.z_wall+bar_width)then
           vbar(i3) = cmplx(bar_height,0.d0)
        endif
     enddo
  endif
END SUBROUTINE esm_gen_vbar
